#include "ThreeCenterElectronRepulsionPrimRecISL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isl,
                                 size_t idx_eri_0_gsl,
                                 size_t idx_eri_1_gsl,
                                 size_t idx_eri_1_hsk,
                                 size_t idx_eri_1_hsl,
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

    /// Set up components of auxilary buffer : GSL

    auto g_xxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl);

    auto g_xxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 1);

    auto g_xxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 2);

    auto g_xxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 3);

    auto g_xxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 4);

    auto g_xxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 5);

    auto g_xxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 6);

    auto g_xxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 7);

    auto g_xxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 8);

    auto g_xxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 9);

    auto g_xxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 10);

    auto g_xxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 11);

    auto g_xxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 12);

    auto g_xxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 13);

    auto g_xxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 14);

    auto g_xxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 15);

    auto g_xxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 16);

    auto g_xxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 17);

    auto g_xxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 18);

    auto g_xxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 19);

    auto g_xxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 20);

    auto g_xxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 21);

    auto g_xxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 22);

    auto g_xxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 23);

    auto g_xxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 24);

    auto g_xxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 25);

    auto g_xxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 26);

    auto g_xxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 27);

    auto g_xxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 28);

    auto g_xxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 29);

    auto g_xxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 30);

    auto g_xxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 31);

    auto g_xxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 32);

    auto g_xxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 33);

    auto g_xxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 34);

    auto g_xxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 35);

    auto g_xxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 36);

    auto g_xxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 37);

    auto g_xxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 38);

    auto g_xxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 39);

    auto g_xxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 40);

    auto g_xxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 41);

    auto g_xxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 42);

    auto g_xxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 43);

    auto g_xxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 44);

    auto g_xxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 45);

    auto g_xxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 47);

    auto g_xxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 50);

    auto g_xxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 54);

    auto g_xxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 59);

    auto g_xxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 65);

    auto g_xxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 72);

    auto g_xxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 80);

    auto g_xxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 90);

    auto g_xxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 91);

    auto g_xxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 93);

    auto g_xxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 96);

    auto g_xxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 100);

    auto g_xxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 105);

    auto g_xxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 111);

    auto g_xxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 118);

    auto g_xxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 135);

    auto g_xxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 136);

    auto g_xxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 137);

    auto g_xxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 138);

    auto g_xxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 139);

    auto g_xxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 140);

    auto g_xxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 141);

    auto g_xxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 142);

    auto g_xxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 143);

    auto g_xxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 144);

    auto g_xxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 145);

    auto g_xxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 146);

    auto g_xxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 147);

    auto g_xxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 148);

    auto g_xxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 149);

    auto g_xxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 150);

    auto g_xxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 151);

    auto g_xxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 152);

    auto g_xxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 153);

    auto g_xxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 154);

    auto g_xxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 155);

    auto g_xxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 156);

    auto g_xxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 157);

    auto g_xxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 158);

    auto g_xxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 159);

    auto g_xxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 160);

    auto g_xxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 161);

    auto g_xxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 162);

    auto g_xxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 163);

    auto g_xxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 164);

    auto g_xxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 165);

    auto g_xxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 166);

    auto g_xxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 167);

    auto g_xxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 168);

    auto g_xxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 169);

    auto g_xxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 170);

    auto g_xxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 171);

    auto g_xxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 172);

    auto g_xxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 173);

    auto g_xxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 174);

    auto g_xxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 175);

    auto g_xxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 176);

    auto g_xxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 177);

    auto g_xxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 178);

    auto g_xxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 179);

    auto g_xxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 225);

    auto g_xxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 226);

    auto g_xxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 227);

    auto g_xxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 228);

    auto g_xxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 229);

    auto g_xxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 230);

    auto g_xxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 231);

    auto g_xxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 232);

    auto g_xxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 233);

    auto g_xxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 234);

    auto g_xxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 235);

    auto g_xxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 236);

    auto g_xxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 237);

    auto g_xxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 238);

    auto g_xxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 239);

    auto g_xxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 240);

    auto g_xxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 241);

    auto g_xxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 242);

    auto g_xxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 243);

    auto g_xxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 244);

    auto g_xxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 245);

    auto g_xxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 246);

    auto g_xxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 247);

    auto g_xxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 248);

    auto g_xxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 249);

    auto g_xxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 250);

    auto g_xxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 251);

    auto g_xxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 252);

    auto g_xxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 253);

    auto g_xxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 254);

    auto g_xxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 255);

    auto g_xxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 256);

    auto g_xxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 257);

    auto g_xxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 258);

    auto g_xxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 259);

    auto g_xxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 260);

    auto g_xxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 261);

    auto g_xxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 262);

    auto g_xxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 263);

    auto g_xxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 264);

    auto g_xxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 265);

    auto g_xxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 266);

    auto g_xxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 267);

    auto g_xxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 268);

    auto g_xxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 269);

    auto g_xyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 271);

    auto g_xyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 273);

    auto g_xyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 274);

    auto g_xyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 276);

    auto g_xyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 277);

    auto g_xyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 278);

    auto g_xyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 280);

    auto g_xyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 281);

    auto g_xyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 282);

    auto g_xyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 283);

    auto g_xyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 285);

    auto g_xyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 286);

    auto g_xyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 287);

    auto g_xyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 288);

    auto g_xyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 289);

    auto g_xyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 291);

    auto g_xyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 292);

    auto g_xyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 293);

    auto g_xyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 294);

    auto g_xyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 295);

    auto g_xyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 296);

    auto g_xyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 298);

    auto g_xyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 299);

    auto g_xyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 300);

    auto g_xyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 301);

    auto g_xyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 302);

    auto g_xyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 303);

    auto g_xyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 304);

    auto g_xyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 306);

    auto g_xyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 307);

    auto g_xyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 308);

    auto g_xyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 309);

    auto g_xyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 310);

    auto g_xyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 311);

    auto g_xyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 312);

    auto g_xyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 313);

    auto g_xyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 314);

    auto g_xzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 407);

    auto g_xzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 409);

    auto g_xzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 410);

    auto g_xzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 412);

    auto g_xzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 413);

    auto g_xzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 414);

    auto g_xzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 416);

    auto g_xzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 417);

    auto g_xzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 418);

    auto g_xzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 419);

    auto g_xzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 421);

    auto g_xzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 422);

    auto g_xzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 423);

    auto g_xzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 424);

    auto g_xzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 425);

    auto g_xzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 427);

    auto g_xzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 428);

    auto g_xzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 429);

    auto g_xzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 430);

    auto g_xzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 431);

    auto g_xzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 432);

    auto g_xzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 434);

    auto g_xzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 435);

    auto g_xzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 436);

    auto g_xzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 437);

    auto g_xzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 438);

    auto g_xzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 439);

    auto g_xzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 440);

    auto g_xzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 441);

    auto g_xzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 442);

    auto g_xzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 443);

    auto g_xzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 444);

    auto g_xzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 445);

    auto g_xzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 446);

    auto g_xzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 447);

    auto g_xzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 448);

    auto g_xzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 449);

    auto g_yyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 450);

    auto g_yyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 451);

    auto g_yyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 452);

    auto g_yyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 453);

    auto g_yyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 454);

    auto g_yyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 455);

    auto g_yyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 456);

    auto g_yyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 457);

    auto g_yyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 458);

    auto g_yyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 459);

    auto g_yyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 460);

    auto g_yyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 461);

    auto g_yyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 462);

    auto g_yyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 463);

    auto g_yyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 464);

    auto g_yyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 465);

    auto g_yyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 466);

    auto g_yyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 467);

    auto g_yyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 468);

    auto g_yyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 469);

    auto g_yyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 470);

    auto g_yyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 471);

    auto g_yyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 472);

    auto g_yyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 473);

    auto g_yyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 474);

    auto g_yyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 475);

    auto g_yyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 476);

    auto g_yyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 477);

    auto g_yyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 478);

    auto g_yyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 479);

    auto g_yyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 480);

    auto g_yyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 481);

    auto g_yyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 482);

    auto g_yyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 483);

    auto g_yyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 484);

    auto g_yyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 485);

    auto g_yyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 486);

    auto g_yyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 487);

    auto g_yyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 488);

    auto g_yyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 489);

    auto g_yyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 490);

    auto g_yyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 491);

    auto g_yyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 492);

    auto g_yyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 493);

    auto g_yyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 494);

    auto g_yyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 496);

    auto g_yyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 498);

    auto g_yyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 501);

    auto g_yyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 505);

    auto g_yyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 510);

    auto g_yyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 516);

    auto g_yyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 523);

    auto g_yyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 531);

    auto g_yyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 540);

    auto g_yyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 541);

    auto g_yyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 542);

    auto g_yyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 543);

    auto g_yyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 544);

    auto g_yyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 545);

    auto g_yyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 546);

    auto g_yyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 547);

    auto g_yyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 548);

    auto g_yyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 549);

    auto g_yyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 550);

    auto g_yyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 551);

    auto g_yyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 552);

    auto g_yyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 553);

    auto g_yyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 554);

    auto g_yyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 555);

    auto g_yyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 556);

    auto g_yyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 557);

    auto g_yyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 558);

    auto g_yyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 559);

    auto g_yyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 560);

    auto g_yyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 561);

    auto g_yyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 562);

    auto g_yyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 563);

    auto g_yyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 564);

    auto g_yyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 565);

    auto g_yyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 566);

    auto g_yyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 567);

    auto g_yyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 568);

    auto g_yyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 569);

    auto g_yyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 570);

    auto g_yyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 571);

    auto g_yyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 572);

    auto g_yyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 573);

    auto g_yyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 574);

    auto g_yyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 575);

    auto g_yyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 576);

    auto g_yyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 577);

    auto g_yyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 578);

    auto g_yyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 579);

    auto g_yyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 580);

    auto g_yyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 581);

    auto g_yyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 582);

    auto g_yyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 583);

    auto g_yyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 584);

    auto g_yzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 585);

    auto g_yzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 587);

    auto g_yzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 589);

    auto g_yzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 590);

    auto g_yzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 592);

    auto g_yzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 593);

    auto g_yzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 594);

    auto g_yzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 596);

    auto g_yzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 597);

    auto g_yzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 598);

    auto g_yzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 599);

    auto g_yzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 601);

    auto g_yzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 602);

    auto g_yzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 603);

    auto g_yzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 604);

    auto g_yzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 605);

    auto g_yzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 607);

    auto g_yzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 608);

    auto g_yzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 609);

    auto g_yzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 610);

    auto g_yzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 611);

    auto g_yzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 612);

    auto g_yzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 614);

    auto g_yzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 615);

    auto g_yzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 616);

    auto g_yzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 617);

    auto g_yzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 618);

    auto g_yzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 619);

    auto g_yzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 620);

    auto g_yzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 622);

    auto g_yzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 623);

    auto g_yzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 624);

    auto g_yzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 625);

    auto g_yzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 626);

    auto g_yzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 627);

    auto g_yzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 628);

    auto g_yzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 629);

    auto g_zzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 630);

    auto g_zzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 631);

    auto g_zzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 632);

    auto g_zzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 633);

    auto g_zzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 634);

    auto g_zzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 635);

    auto g_zzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 636);

    auto g_zzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 637);

    auto g_zzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 638);

    auto g_zzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 639);

    auto g_zzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 640);

    auto g_zzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 641);

    auto g_zzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 642);

    auto g_zzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 643);

    auto g_zzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 644);

    auto g_zzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 645);

    auto g_zzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 646);

    auto g_zzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 647);

    auto g_zzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 648);

    auto g_zzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 649);

    auto g_zzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 650);

    auto g_zzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 651);

    auto g_zzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 652);

    auto g_zzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 653);

    auto g_zzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 654);

    auto g_zzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 655);

    auto g_zzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 656);

    auto g_zzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 657);

    auto g_zzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 658);

    auto g_zzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 659);

    auto g_zzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 660);

    auto g_zzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 661);

    auto g_zzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 662);

    auto g_zzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 663);

    auto g_zzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 664);

    auto g_zzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 665);

    auto g_zzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 666);

    auto g_zzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 667);

    auto g_zzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 668);

    auto g_zzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 669);

    auto g_zzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 670);

    auto g_zzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 671);

    auto g_zzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 672);

    auto g_zzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 673);

    auto g_zzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 674);

    /// Set up components of auxilary buffer : GSL

    auto g_xxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl);

    auto g_xxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 1);

    auto g_xxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 2);

    auto g_xxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 3);

    auto g_xxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 4);

    auto g_xxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 5);

    auto g_xxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 6);

    auto g_xxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 7);

    auto g_xxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 8);

    auto g_xxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 9);

    auto g_xxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 10);

    auto g_xxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 11);

    auto g_xxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 12);

    auto g_xxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 13);

    auto g_xxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 14);

    auto g_xxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 15);

    auto g_xxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 16);

    auto g_xxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 17);

    auto g_xxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 18);

    auto g_xxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 19);

    auto g_xxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 20);

    auto g_xxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 21);

    auto g_xxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 22);

    auto g_xxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 23);

    auto g_xxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 24);

    auto g_xxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 25);

    auto g_xxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 26);

    auto g_xxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 27);

    auto g_xxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 28);

    auto g_xxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 29);

    auto g_xxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 30);

    auto g_xxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 31);

    auto g_xxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 32);

    auto g_xxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 33);

    auto g_xxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 34);

    auto g_xxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 35);

    auto g_xxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 36);

    auto g_xxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 37);

    auto g_xxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 38);

    auto g_xxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 39);

    auto g_xxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 40);

    auto g_xxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 41);

    auto g_xxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 42);

    auto g_xxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 43);

    auto g_xxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 44);

    auto g_xxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 45);

    auto g_xxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 47);

    auto g_xxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 50);

    auto g_xxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 54);

    auto g_xxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 59);

    auto g_xxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 65);

    auto g_xxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 72);

    auto g_xxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 80);

    auto g_xxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 90);

    auto g_xxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 91);

    auto g_xxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 93);

    auto g_xxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 96);

    auto g_xxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 100);

    auto g_xxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 105);

    auto g_xxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 111);

    auto g_xxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 118);

    auto g_xxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 135);

    auto g_xxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 136);

    auto g_xxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 137);

    auto g_xxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 138);

    auto g_xxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 139);

    auto g_xxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 140);

    auto g_xxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 141);

    auto g_xxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 142);

    auto g_xxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 143);

    auto g_xxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 144);

    auto g_xxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 145);

    auto g_xxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 146);

    auto g_xxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 147);

    auto g_xxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 148);

    auto g_xxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 149);

    auto g_xxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 150);

    auto g_xxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 151);

    auto g_xxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 152);

    auto g_xxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 153);

    auto g_xxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 154);

    auto g_xxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 155);

    auto g_xxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 156);

    auto g_xxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 157);

    auto g_xxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 158);

    auto g_xxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 159);

    auto g_xxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 160);

    auto g_xxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 161);

    auto g_xxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 162);

    auto g_xxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 163);

    auto g_xxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 164);

    auto g_xxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 165);

    auto g_xxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 166);

    auto g_xxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 167);

    auto g_xxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 168);

    auto g_xxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 169);

    auto g_xxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 170);

    auto g_xxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 171);

    auto g_xxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 172);

    auto g_xxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 173);

    auto g_xxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 174);

    auto g_xxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 175);

    auto g_xxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 176);

    auto g_xxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 177);

    auto g_xxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 178);

    auto g_xxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 179);

    auto g_xxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 225);

    auto g_xxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 226);

    auto g_xxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 227);

    auto g_xxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 228);

    auto g_xxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 229);

    auto g_xxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 230);

    auto g_xxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 231);

    auto g_xxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 232);

    auto g_xxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 233);

    auto g_xxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 234);

    auto g_xxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 235);

    auto g_xxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 236);

    auto g_xxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 237);

    auto g_xxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 238);

    auto g_xxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 239);

    auto g_xxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 240);

    auto g_xxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 241);

    auto g_xxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 242);

    auto g_xxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 243);

    auto g_xxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 244);

    auto g_xxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 245);

    auto g_xxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 246);

    auto g_xxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 247);

    auto g_xxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 248);

    auto g_xxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 249);

    auto g_xxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 250);

    auto g_xxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 251);

    auto g_xxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 252);

    auto g_xxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 253);

    auto g_xxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 254);

    auto g_xxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 255);

    auto g_xxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 256);

    auto g_xxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 257);

    auto g_xxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 258);

    auto g_xxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 259);

    auto g_xxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 260);

    auto g_xxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 261);

    auto g_xxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 262);

    auto g_xxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 263);

    auto g_xxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 264);

    auto g_xxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 265);

    auto g_xxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 266);

    auto g_xxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 267);

    auto g_xxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 268);

    auto g_xxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 269);

    auto g_xyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 271);

    auto g_xyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 273);

    auto g_xyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 274);

    auto g_xyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 276);

    auto g_xyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 277);

    auto g_xyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 278);

    auto g_xyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 280);

    auto g_xyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 281);

    auto g_xyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 282);

    auto g_xyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 283);

    auto g_xyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 285);

    auto g_xyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 286);

    auto g_xyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 287);

    auto g_xyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 288);

    auto g_xyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 289);

    auto g_xyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 291);

    auto g_xyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 292);

    auto g_xyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 293);

    auto g_xyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 294);

    auto g_xyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 295);

    auto g_xyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 296);

    auto g_xyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 298);

    auto g_xyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 299);

    auto g_xyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 300);

    auto g_xyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 301);

    auto g_xyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 302);

    auto g_xyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 303);

    auto g_xyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 304);

    auto g_xyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 306);

    auto g_xyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 307);

    auto g_xyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 308);

    auto g_xyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 309);

    auto g_xyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 310);

    auto g_xyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 311);

    auto g_xyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 312);

    auto g_xyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 313);

    auto g_xyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 314);

    auto g_xzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 407);

    auto g_xzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 409);

    auto g_xzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 410);

    auto g_xzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 412);

    auto g_xzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 413);

    auto g_xzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 414);

    auto g_xzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 416);

    auto g_xzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 417);

    auto g_xzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 418);

    auto g_xzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 419);

    auto g_xzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 421);

    auto g_xzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 422);

    auto g_xzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 423);

    auto g_xzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 424);

    auto g_xzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 425);

    auto g_xzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 427);

    auto g_xzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 428);

    auto g_xzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 429);

    auto g_xzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 430);

    auto g_xzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 431);

    auto g_xzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 432);

    auto g_xzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 434);

    auto g_xzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 435);

    auto g_xzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 436);

    auto g_xzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 437);

    auto g_xzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 438);

    auto g_xzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 439);

    auto g_xzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 440);

    auto g_xzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 441);

    auto g_xzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 442);

    auto g_xzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 443);

    auto g_xzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 444);

    auto g_xzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 445);

    auto g_xzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 446);

    auto g_xzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 447);

    auto g_xzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 448);

    auto g_xzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 449);

    auto g_yyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 450);

    auto g_yyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 451);

    auto g_yyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 452);

    auto g_yyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 453);

    auto g_yyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 454);

    auto g_yyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 455);

    auto g_yyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 456);

    auto g_yyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 457);

    auto g_yyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 458);

    auto g_yyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 459);

    auto g_yyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 460);

    auto g_yyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 461);

    auto g_yyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 462);

    auto g_yyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 463);

    auto g_yyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 464);

    auto g_yyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 465);

    auto g_yyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 466);

    auto g_yyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 467);

    auto g_yyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 468);

    auto g_yyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 469);

    auto g_yyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 470);

    auto g_yyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 471);

    auto g_yyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 472);

    auto g_yyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 473);

    auto g_yyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 474);

    auto g_yyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 475);

    auto g_yyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 476);

    auto g_yyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 477);

    auto g_yyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 478);

    auto g_yyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 479);

    auto g_yyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 480);

    auto g_yyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 481);

    auto g_yyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 482);

    auto g_yyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 483);

    auto g_yyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 484);

    auto g_yyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 485);

    auto g_yyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 486);

    auto g_yyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 487);

    auto g_yyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 488);

    auto g_yyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 489);

    auto g_yyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 490);

    auto g_yyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 491);

    auto g_yyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 492);

    auto g_yyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 493);

    auto g_yyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 494);

    auto g_yyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 496);

    auto g_yyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 498);

    auto g_yyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 501);

    auto g_yyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 505);

    auto g_yyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 510);

    auto g_yyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 516);

    auto g_yyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 523);

    auto g_yyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 531);

    auto g_yyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 540);

    auto g_yyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 541);

    auto g_yyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 542);

    auto g_yyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 543);

    auto g_yyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 544);

    auto g_yyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 545);

    auto g_yyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 546);

    auto g_yyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 547);

    auto g_yyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 548);

    auto g_yyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 549);

    auto g_yyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 550);

    auto g_yyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 551);

    auto g_yyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 552);

    auto g_yyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 553);

    auto g_yyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 554);

    auto g_yyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 555);

    auto g_yyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 556);

    auto g_yyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 557);

    auto g_yyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 558);

    auto g_yyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 559);

    auto g_yyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 560);

    auto g_yyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 561);

    auto g_yyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 562);

    auto g_yyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 563);

    auto g_yyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 564);

    auto g_yyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 565);

    auto g_yyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 566);

    auto g_yyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 567);

    auto g_yyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 568);

    auto g_yyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 569);

    auto g_yyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 570);

    auto g_yyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 571);

    auto g_yyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 572);

    auto g_yyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 573);

    auto g_yyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 574);

    auto g_yyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 575);

    auto g_yyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 576);

    auto g_yyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 577);

    auto g_yyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 578);

    auto g_yyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 579);

    auto g_yyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 580);

    auto g_yyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 581);

    auto g_yyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 582);

    auto g_yyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 583);

    auto g_yyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 584);

    auto g_yzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 585);

    auto g_yzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 587);

    auto g_yzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 589);

    auto g_yzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 590);

    auto g_yzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 592);

    auto g_yzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 593);

    auto g_yzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 594);

    auto g_yzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 596);

    auto g_yzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 597);

    auto g_yzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 598);

    auto g_yzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 599);

    auto g_yzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 601);

    auto g_yzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 602);

    auto g_yzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 603);

    auto g_yzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 604);

    auto g_yzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 605);

    auto g_yzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 607);

    auto g_yzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 608);

    auto g_yzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 609);

    auto g_yzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 610);

    auto g_yzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 611);

    auto g_yzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 612);

    auto g_yzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 614);

    auto g_yzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 615);

    auto g_yzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 616);

    auto g_yzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 617);

    auto g_yzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 618);

    auto g_yzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 619);

    auto g_yzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 620);

    auto g_yzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 622);

    auto g_yzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 623);

    auto g_yzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 624);

    auto g_yzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 625);

    auto g_yzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 626);

    auto g_yzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 627);

    auto g_yzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 628);

    auto g_yzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 629);

    auto g_zzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 630);

    auto g_zzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 631);

    auto g_zzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 632);

    auto g_zzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 633);

    auto g_zzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 634);

    auto g_zzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 635);

    auto g_zzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 636);

    auto g_zzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 637);

    auto g_zzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 638);

    auto g_zzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 639);

    auto g_zzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 640);

    auto g_zzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 641);

    auto g_zzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 642);

    auto g_zzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 643);

    auto g_zzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 644);

    auto g_zzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 645);

    auto g_zzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 646);

    auto g_zzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 647);

    auto g_zzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 648);

    auto g_zzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 649);

    auto g_zzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 650);

    auto g_zzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 651);

    auto g_zzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 652);

    auto g_zzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 653);

    auto g_zzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 654);

    auto g_zzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 655);

    auto g_zzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 656);

    auto g_zzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 657);

    auto g_zzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 658);

    auto g_zzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 659);

    auto g_zzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 660);

    auto g_zzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 661);

    auto g_zzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 662);

    auto g_zzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 663);

    auto g_zzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 664);

    auto g_zzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 665);

    auto g_zzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 666);

    auto g_zzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 667);

    auto g_zzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 668);

    auto g_zzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 669);

    auto g_zzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 670);

    auto g_zzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 671);

    auto g_zzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 672);

    auto g_zzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 673);

    auto g_zzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 674);

    /// Set up components of auxilary buffer : HSK

    auto g_xxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk);

    auto g_xxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 1);

    auto g_xxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 2);

    auto g_xxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 3);

    auto g_xxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 4);

    auto g_xxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 5);

    auto g_xxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 6);

    auto g_xxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 7);

    auto g_xxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 8);

    auto g_xxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 9);

    auto g_xxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 10);

    auto g_xxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 11);

    auto g_xxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 12);

    auto g_xxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 13);

    auto g_xxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 14);

    auto g_xxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 15);

    auto g_xxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 16);

    auto g_xxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 17);

    auto g_xxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 18);

    auto g_xxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 19);

    auto g_xxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 20);

    auto g_xxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 21);

    auto g_xxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 22);

    auto g_xxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 23);

    auto g_xxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 24);

    auto g_xxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 25);

    auto g_xxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 26);

    auto g_xxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 27);

    auto g_xxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 28);

    auto g_xxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 29);

    auto g_xxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 30);

    auto g_xxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 31);

    auto g_xxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 32);

    auto g_xxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 33);

    auto g_xxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 34);

    auto g_xxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 35);

    auto g_xxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 74);

    auto g_xxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 76);

    auto g_xxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 77);

    auto g_xxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 79);

    auto g_xxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 80);

    auto g_xxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 81);

    auto g_xxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 83);

    auto g_xxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 84);

    auto g_xxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 85);

    auto g_xxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 86);

    auto g_xxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 88);

    auto g_xxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 89);

    auto g_xxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 90);

    auto g_xxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 91);

    auto g_xxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 92);

    auto g_xxxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 94);

    auto g_xxxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 95);

    auto g_xxxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 96);

    auto g_xxxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 97);

    auto g_xxxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 98);

    auto g_xxxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 99);

    auto g_xxxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 101);

    auto g_xxxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 102);

    auto g_xxxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 103);

    auto g_xxxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 104);

    auto g_xxxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 105);

    auto g_xxxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 106);

    auto g_xxxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 107);

    auto g_xxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 108);

    auto g_xxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 109);

    auto g_xxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 110);

    auto g_xxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 111);

    auto g_xxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 112);

    auto g_xxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 113);

    auto g_xxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 114);

    auto g_xxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 115);

    auto g_xxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 116);

    auto g_xxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 117);

    auto g_xxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 118);

    auto g_xxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 119);

    auto g_xxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 120);

    auto g_xxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 121);

    auto g_xxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 122);

    auto g_xxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 123);

    auto g_xxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 124);

    auto g_xxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 125);

    auto g_xxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 126);

    auto g_xxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 127);

    auto g_xxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 128);

    auto g_xxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 129);

    auto g_xxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 130);

    auto g_xxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 131);

    auto g_xxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 132);

    auto g_xxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 133);

    auto g_xxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 134);

    auto g_xxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 135);

    auto g_xxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 136);

    auto g_xxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 137);

    auto g_xxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 138);

    auto g_xxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 139);

    auto g_xxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 140);

    auto g_xxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 141);

    auto g_xxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 142);

    auto g_xxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 143);

    auto g_xxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 180);

    auto g_xxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 181);

    auto g_xxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 182);

    auto g_xxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 183);

    auto g_xxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 184);

    auto g_xxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 185);

    auto g_xxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 186);

    auto g_xxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 187);

    auto g_xxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 188);

    auto g_xxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 189);

    auto g_xxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 190);

    auto g_xxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 191);

    auto g_xxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 192);

    auto g_xxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 193);

    auto g_xxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 194);

    auto g_xxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 195);

    auto g_xxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 196);

    auto g_xxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 197);

    auto g_xxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 198);

    auto g_xxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 199);

    auto g_xxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 200);

    auto g_xxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 201);

    auto g_xxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 202);

    auto g_xxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 203);

    auto g_xxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 204);

    auto g_xxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 205);

    auto g_xxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 206);

    auto g_xxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 207);

    auto g_xxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 208);

    auto g_xxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 209);

    auto g_xxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 210);

    auto g_xxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 211);

    auto g_xxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 212);

    auto g_xxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 213);

    auto g_xxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 214);

    auto g_xxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 215);

    auto g_xxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 216);

    auto g_xxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 217);

    auto g_xxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 218);

    auto g_xxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 219);

    auto g_xxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 220);

    auto g_xxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 221);

    auto g_xxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 222);

    auto g_xxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 223);

    auto g_xxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 224);

    auto g_xxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 225);

    auto g_xxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 226);

    auto g_xxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 227);

    auto g_xxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 228);

    auto g_xxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 229);

    auto g_xxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 230);

    auto g_xxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 231);

    auto g_xxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 232);

    auto g_xxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 233);

    auto g_xxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 234);

    auto g_xxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 235);

    auto g_xxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 236);

    auto g_xxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 237);

    auto g_xxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 238);

    auto g_xxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 239);

    auto g_xxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 240);

    auto g_xxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 241);

    auto g_xxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 242);

    auto g_xxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 243);

    auto g_xxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 244);

    auto g_xxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 245);

    auto g_xxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 246);

    auto g_xxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 247);

    auto g_xxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 248);

    auto g_xxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 249);

    auto g_xxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 250);

    auto g_xxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 251);

    auto g_xxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 324);

    auto g_xxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 325);

    auto g_xxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 326);

    auto g_xxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 327);

    auto g_xxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 328);

    auto g_xxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 329);

    auto g_xxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 330);

    auto g_xxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 331);

    auto g_xxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 332);

    auto g_xxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 333);

    auto g_xxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 334);

    auto g_xxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 335);

    auto g_xxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 336);

    auto g_xxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 337);

    auto g_xxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 338);

    auto g_xxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 339);

    auto g_xxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 340);

    auto g_xxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 341);

    auto g_xxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 342);

    auto g_xxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 343);

    auto g_xxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 344);

    auto g_xxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 345);

    auto g_xxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 346);

    auto g_xxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 347);

    auto g_xxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 348);

    auto g_xxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 349);

    auto g_xxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 350);

    auto g_xxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 351);

    auto g_xxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 352);

    auto g_xxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 353);

    auto g_xxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 354);

    auto g_xxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 355);

    auto g_xxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 356);

    auto g_xxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 357);

    auto g_xxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 358);

    auto g_xxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 359);

    auto g_xyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 361);

    auto g_xyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 363);

    auto g_xyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 364);

    auto g_xyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 366);

    auto g_xyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 367);

    auto g_xyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 368);

    auto g_xyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 370);

    auto g_xyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 371);

    auto g_xyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 372);

    auto g_xyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 373);

    auto g_xyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 375);

    auto g_xyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 376);

    auto g_xyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 377);

    auto g_xyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 378);

    auto g_xyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 379);

    auto g_xyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 381);

    auto g_xyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 382);

    auto g_xyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 383);

    auto g_xyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 384);

    auto g_xyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 385);

    auto g_xyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 386);

    auto g_xyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 388);

    auto g_xyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 389);

    auto g_xyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 390);

    auto g_xyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 391);

    auto g_xyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 392);

    auto g_xyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 393);

    auto g_xyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 394);

    auto g_xyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 436);

    auto g_xyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 439);

    auto g_xyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 440);

    auto g_xyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 443);

    auto g_xyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 444);

    auto g_xyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 445);

    auto g_xyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 448);

    auto g_xyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 449);

    auto g_xyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 450);

    auto g_xyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 451);

    auto g_xyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 454);

    auto g_xyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 455);

    auto g_xyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 456);

    auto g_xyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 457);

    auto g_xyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 458);

    auto g_xyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 461);

    auto g_xyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 462);

    auto g_xyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 463);

    auto g_xyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 464);

    auto g_xyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 465);

    auto g_xyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 466);

    auto g_xzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 506);

    auto g_xzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 508);

    auto g_xzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 509);

    auto g_xzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 511);

    auto g_xzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 512);

    auto g_xzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 513);

    auto g_xzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 515);

    auto g_xzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 516);

    auto g_xzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 517);

    auto g_xzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 518);

    auto g_xzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 520);

    auto g_xzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 521);

    auto g_xzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 522);

    auto g_xzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 523);

    auto g_xzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 524);

    auto g_xzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 526);

    auto g_xzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 527);

    auto g_xzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 528);

    auto g_xzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 529);

    auto g_xzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 530);

    auto g_xzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 531);

    auto g_xzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 533);

    auto g_xzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 534);

    auto g_xzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 535);

    auto g_xzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 536);

    auto g_xzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 537);

    auto g_xzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 538);

    auto g_xzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 539);

    auto g_yyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 540);

    auto g_yyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 541);

    auto g_yyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 542);

    auto g_yyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 543);

    auto g_yyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 544);

    auto g_yyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 545);

    auto g_yyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 546);

    auto g_yyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 547);

    auto g_yyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 548);

    auto g_yyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 549);

    auto g_yyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 550);

    auto g_yyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 551);

    auto g_yyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 552);

    auto g_yyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 553);

    auto g_yyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 554);

    auto g_yyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 555);

    auto g_yyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 556);

    auto g_yyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 557);

    auto g_yyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 558);

    auto g_yyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 559);

    auto g_yyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 560);

    auto g_yyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 561);

    auto g_yyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 562);

    auto g_yyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 563);

    auto g_yyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 564);

    auto g_yyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 565);

    auto g_yyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 566);

    auto g_yyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 567);

    auto g_yyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 568);

    auto g_yyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 569);

    auto g_yyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 570);

    auto g_yyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 571);

    auto g_yyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 572);

    auto g_yyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 573);

    auto g_yyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 574);

    auto g_yyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 575);

    auto g_yyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 578);

    auto g_yyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 580);

    auto g_yyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 581);

    auto g_yyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 583);

    auto g_yyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 584);

    auto g_yyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 585);

    auto g_yyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 587);

    auto g_yyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 588);

    auto g_yyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 589);

    auto g_yyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 590);

    auto g_yyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 592);

    auto g_yyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 593);

    auto g_yyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 594);

    auto g_yyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 595);

    auto g_yyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 596);

    auto g_yyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 598);

    auto g_yyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 599);

    auto g_yyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 600);

    auto g_yyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 601);

    auto g_yyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 602);

    auto g_yyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 603);

    auto g_yyyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 605);

    auto g_yyyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 606);

    auto g_yyyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 607);

    auto g_yyyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 608);

    auto g_yyyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 609);

    auto g_yyyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 610);

    auto g_yyyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 611);

    auto g_yyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 612);

    auto g_yyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 613);

    auto g_yyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 614);

    auto g_yyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 615);

    auto g_yyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 616);

    auto g_yyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 617);

    auto g_yyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 618);

    auto g_yyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 619);

    auto g_yyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 620);

    auto g_yyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 621);

    auto g_yyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 622);

    auto g_yyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 623);

    auto g_yyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 624);

    auto g_yyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 625);

    auto g_yyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 626);

    auto g_yyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 627);

    auto g_yyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 628);

    auto g_yyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 629);

    auto g_yyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 630);

    auto g_yyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 631);

    auto g_yyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 632);

    auto g_yyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 633);

    auto g_yyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 634);

    auto g_yyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 635);

    auto g_yyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 636);

    auto g_yyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 637);

    auto g_yyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 638);

    auto g_yyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 639);

    auto g_yyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 640);

    auto g_yyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 641);

    auto g_yyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 642);

    auto g_yyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 643);

    auto g_yyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 644);

    auto g_yyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 645);

    auto g_yyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 646);

    auto g_yyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 647);

    auto g_yyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 648);

    auto g_yyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 649);

    auto g_yyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 650);

    auto g_yyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 651);

    auto g_yyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 652);

    auto g_yyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 653);

    auto g_yyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 654);

    auto g_yyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 655);

    auto g_yyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 656);

    auto g_yyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 657);

    auto g_yyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 658);

    auto g_yyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 659);

    auto g_yyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 660);

    auto g_yyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 661);

    auto g_yyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 662);

    auto g_yyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 663);

    auto g_yyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 664);

    auto g_yyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 665);

    auto g_yyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 666);

    auto g_yyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 667);

    auto g_yyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 668);

    auto g_yyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 669);

    auto g_yyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 670);

    auto g_yyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 671);

    auto g_yyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 672);

    auto g_yyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 673);

    auto g_yyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 674);

    auto g_yyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 675);

    auto g_yyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 676);

    auto g_yyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 677);

    auto g_yyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 678);

    auto g_yyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 679);

    auto g_yyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 680);

    auto g_yyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 681);

    auto g_yyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 682);

    auto g_yyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 683);

    auto g_yzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 685);

    auto g_yzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 686);

    auto g_yzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 687);

    auto g_yzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 688);

    auto g_yzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 689);

    auto g_yzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 690);

    auto g_yzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 691);

    auto g_yzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 692);

    auto g_yzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 693);

    auto g_yzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 694);

    auto g_yzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 695);

    auto g_yzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 696);

    auto g_yzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 697);

    auto g_yzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 698);

    auto g_yzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 699);

    auto g_yzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 700);

    auto g_yzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 701);

    auto g_yzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 702);

    auto g_yzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 703);

    auto g_yzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 704);

    auto g_yzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 705);

    auto g_yzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 706);

    auto g_yzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 707);

    auto g_yzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 708);

    auto g_yzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 709);

    auto g_yzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 710);

    auto g_yzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 711);

    auto g_yzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 712);

    auto g_yzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 713);

    auto g_yzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 714);

    auto g_yzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 715);

    auto g_yzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 716);

    auto g_yzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 717);

    auto g_yzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 718);

    auto g_yzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 719);

    auto g_zzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 720);

    auto g_zzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 721);

    auto g_zzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 722);

    auto g_zzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 723);

    auto g_zzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 724);

    auto g_zzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 725);

    auto g_zzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 726);

    auto g_zzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 727);

    auto g_zzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 728);

    auto g_zzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 729);

    auto g_zzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 730);

    auto g_zzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 731);

    auto g_zzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 732);

    auto g_zzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 733);

    auto g_zzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 734);

    auto g_zzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 735);

    auto g_zzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 736);

    auto g_zzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 737);

    auto g_zzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 738);

    auto g_zzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 739);

    auto g_zzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 740);

    auto g_zzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 741);

    auto g_zzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 742);

    auto g_zzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 743);

    auto g_zzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 744);

    auto g_zzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 745);

    auto g_zzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 746);

    auto g_zzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 747);

    auto g_zzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 748);

    auto g_zzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 749);

    auto g_zzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 750);

    auto g_zzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 751);

    auto g_zzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 752);

    auto g_zzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 753);

    auto g_zzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 754);

    auto g_zzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 755);

    /// Set up components of auxilary buffer : HSL

    auto g_xxxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl);

    auto g_xxxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 1);

    auto g_xxxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 2);

    auto g_xxxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 3);

    auto g_xxxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 4);

    auto g_xxxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 5);

    auto g_xxxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 6);

    auto g_xxxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 7);

    auto g_xxxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 8);

    auto g_xxxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 9);

    auto g_xxxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 10);

    auto g_xxxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 11);

    auto g_xxxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 12);

    auto g_xxxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 13);

    auto g_xxxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 14);

    auto g_xxxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 15);

    auto g_xxxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 16);

    auto g_xxxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 17);

    auto g_xxxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 18);

    auto g_xxxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 19);

    auto g_xxxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 20);

    auto g_xxxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 21);

    auto g_xxxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 22);

    auto g_xxxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 23);

    auto g_xxxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 24);

    auto g_xxxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 25);

    auto g_xxxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 26);

    auto g_xxxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 27);

    auto g_xxxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 28);

    auto g_xxxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 29);

    auto g_xxxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 30);

    auto g_xxxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 31);

    auto g_xxxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 32);

    auto g_xxxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 33);

    auto g_xxxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 34);

    auto g_xxxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 35);

    auto g_xxxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 36);

    auto g_xxxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 37);

    auto g_xxxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 38);

    auto g_xxxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 39);

    auto g_xxxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 40);

    auto g_xxxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 41);

    auto g_xxxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 42);

    auto g_xxxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 43);

    auto g_xxxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 44);

    auto g_xxxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 45);

    auto g_xxxxy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 46);

    auto g_xxxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 47);

    auto g_xxxxy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 48);

    auto g_xxxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 50);

    auto g_xxxxy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 51);

    auto g_xxxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 54);

    auto g_xxxxy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 55);

    auto g_xxxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 59);

    auto g_xxxxy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 60);

    auto g_xxxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 65);

    auto g_xxxxy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 66);

    auto g_xxxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 72);

    auto g_xxxxy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 73);

    auto g_xxxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 80);

    auto g_xxxxy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 81);

    auto g_xxxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 90);

    auto g_xxxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 91);

    auto g_xxxxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 92);

    auto g_xxxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 93);

    auto g_xxxxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 94);

    auto g_xxxxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 95);

    auto g_xxxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 96);

    auto g_xxxxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 97);

    auto g_xxxxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 98);

    auto g_xxxxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 99);

    auto g_xxxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 100);

    auto g_xxxxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 101);

    auto g_xxxxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 102);

    auto g_xxxxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 103);

    auto g_xxxxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 104);

    auto g_xxxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 105);

    auto g_xxxxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 106);

    auto g_xxxxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 107);

    auto g_xxxxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 108);

    auto g_xxxxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 109);

    auto g_xxxxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 110);

    auto g_xxxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 111);

    auto g_xxxxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 112);

    auto g_xxxxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 113);

    auto g_xxxxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 114);

    auto g_xxxxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 115);

    auto g_xxxxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 116);

    auto g_xxxxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 117);

    auto g_xxxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 118);

    auto g_xxxxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 119);

    auto g_xxxxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 120);

    auto g_xxxxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 121);

    auto g_xxxxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 122);

    auto g_xxxxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 123);

    auto g_xxxxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 124);

    auto g_xxxxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 125);

    auto g_xxxxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 127);

    auto g_xxxxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 128);

    auto g_xxxxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 129);

    auto g_xxxxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 130);

    auto g_xxxxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 131);

    auto g_xxxxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 132);

    auto g_xxxxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 133);

    auto g_xxxxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 134);

    auto g_xxxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 135);

    auto g_xxxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 136);

    auto g_xxxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 137);

    auto g_xxxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 138);

    auto g_xxxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 139);

    auto g_xxxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 140);

    auto g_xxxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 141);

    auto g_xxxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 142);

    auto g_xxxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 143);

    auto g_xxxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 144);

    auto g_xxxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 145);

    auto g_xxxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 146);

    auto g_xxxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 147);

    auto g_xxxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 148);

    auto g_xxxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 149);

    auto g_xxxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 150);

    auto g_xxxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 151);

    auto g_xxxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 152);

    auto g_xxxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 153);

    auto g_xxxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 154);

    auto g_xxxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 155);

    auto g_xxxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 156);

    auto g_xxxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 157);

    auto g_xxxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 158);

    auto g_xxxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 159);

    auto g_xxxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 160);

    auto g_xxxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 161);

    auto g_xxxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 162);

    auto g_xxxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 163);

    auto g_xxxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 164);

    auto g_xxxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 165);

    auto g_xxxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 166);

    auto g_xxxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 167);

    auto g_xxxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 168);

    auto g_xxxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 169);

    auto g_xxxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 170);

    auto g_xxxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 171);

    auto g_xxxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 172);

    auto g_xxxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 173);

    auto g_xxxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 174);

    auto g_xxxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 175);

    auto g_xxxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 176);

    auto g_xxxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 177);

    auto g_xxxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 178);

    auto g_xxxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 179);

    auto g_xxxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 225);

    auto g_xxxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 226);

    auto g_xxxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 227);

    auto g_xxxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 228);

    auto g_xxxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 229);

    auto g_xxxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 230);

    auto g_xxxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 231);

    auto g_xxxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 232);

    auto g_xxxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 233);

    auto g_xxxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 234);

    auto g_xxxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 235);

    auto g_xxxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 236);

    auto g_xxxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 237);

    auto g_xxxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 238);

    auto g_xxxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 239);

    auto g_xxxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 240);

    auto g_xxxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 241);

    auto g_xxxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 242);

    auto g_xxxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 243);

    auto g_xxxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 244);

    auto g_xxxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 245);

    auto g_xxxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 246);

    auto g_xxxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 247);

    auto g_xxxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 248);

    auto g_xxxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 249);

    auto g_xxxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 250);

    auto g_xxxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 251);

    auto g_xxxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 252);

    auto g_xxxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 253);

    auto g_xxxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 254);

    auto g_xxxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 255);

    auto g_xxxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 256);

    auto g_xxxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 257);

    auto g_xxxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 258);

    auto g_xxxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 259);

    auto g_xxxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 260);

    auto g_xxxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 261);

    auto g_xxxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 262);

    auto g_xxxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 263);

    auto g_xxxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 264);

    auto g_xxxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 265);

    auto g_xxxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 266);

    auto g_xxxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 267);

    auto g_xxxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 268);

    auto g_xxxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 269);

    auto g_xxyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 270);

    auto g_xxyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 271);

    auto g_xxyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 272);

    auto g_xxyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 273);

    auto g_xxyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 274);

    auto g_xxyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 275);

    auto g_xxyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 276);

    auto g_xxyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 277);

    auto g_xxyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 278);

    auto g_xxyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 279);

    auto g_xxyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 280);

    auto g_xxyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 281);

    auto g_xxyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 282);

    auto g_xxyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 283);

    auto g_xxyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 284);

    auto g_xxyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 285);

    auto g_xxyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 286);

    auto g_xxyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 287);

    auto g_xxyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 288);

    auto g_xxyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 289);

    auto g_xxyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 290);

    auto g_xxyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 291);

    auto g_xxyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 292);

    auto g_xxyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 293);

    auto g_xxyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 294);

    auto g_xxyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 295);

    auto g_xxyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 296);

    auto g_xxyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 297);

    auto g_xxyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 298);

    auto g_xxyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 299);

    auto g_xxyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 300);

    auto g_xxyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 301);

    auto g_xxyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 302);

    auto g_xxyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 303);

    auto g_xxyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 304);

    auto g_xxyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 305);

    auto g_xxyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 306);

    auto g_xxyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 307);

    auto g_xxyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 308);

    auto g_xxyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 309);

    auto g_xxyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 310);

    auto g_xxyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 311);

    auto g_xxyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 312);

    auto g_xxyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 313);

    auto g_xxyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 314);

    auto g_xxyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 316);

    auto g_xxyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 318);

    auto g_xxyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 321);

    auto g_xxyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 325);

    auto g_xxyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 330);

    auto g_xxyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 336);

    auto g_xxyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 343);

    auto g_xxyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 360);

    auto g_xxyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 362);

    auto g_xxyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 365);

    auto g_xxyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 369);

    auto g_xxyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 374);

    auto g_xxyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 380);

    auto g_xxyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 387);

    auto g_xxyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 395);

    auto g_xxzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 405);

    auto g_xxzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 406);

    auto g_xxzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 407);

    auto g_xxzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 408);

    auto g_xxzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 409);

    auto g_xxzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 410);

    auto g_xxzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 411);

    auto g_xxzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 412);

    auto g_xxzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 413);

    auto g_xxzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 414);

    auto g_xxzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 415);

    auto g_xxzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 416);

    auto g_xxzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 417);

    auto g_xxzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 418);

    auto g_xxzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 419);

    auto g_xxzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 420);

    auto g_xxzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 421);

    auto g_xxzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 422);

    auto g_xxzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 423);

    auto g_xxzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 424);

    auto g_xxzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 425);

    auto g_xxzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 426);

    auto g_xxzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 427);

    auto g_xxzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 428);

    auto g_xxzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 429);

    auto g_xxzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 430);

    auto g_xxzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 431);

    auto g_xxzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 432);

    auto g_xxzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 433);

    auto g_xxzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 434);

    auto g_xxzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 435);

    auto g_xxzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 436);

    auto g_xxzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 437);

    auto g_xxzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 438);

    auto g_xxzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 439);

    auto g_xxzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 440);

    auto g_xxzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 441);

    auto g_xxzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 442);

    auto g_xxzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 443);

    auto g_xxzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 444);

    auto g_xxzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 445);

    auto g_xxzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 446);

    auto g_xxzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 447);

    auto g_xxzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 448);

    auto g_xxzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 449);

    auto g_xyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 450);

    auto g_xyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 451);

    auto g_xyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 453);

    auto g_xyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 454);

    auto g_xyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 456);

    auto g_xyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 457);

    auto g_xyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 458);

    auto g_xyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 460);

    auto g_xyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 461);

    auto g_xyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 462);

    auto g_xyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 463);

    auto g_xyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 465);

    auto g_xyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 466);

    auto g_xyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 467);

    auto g_xyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 468);

    auto g_xyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 469);

    auto g_xyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 471);

    auto g_xyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 472);

    auto g_xyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 473);

    auto g_xyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 474);

    auto g_xyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 475);

    auto g_xyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 476);

    auto g_xyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 478);

    auto g_xyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 479);

    auto g_xyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 480);

    auto g_xyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 481);

    auto g_xyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 482);

    auto g_xyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 483);

    auto g_xyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 484);

    auto g_xyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 486);

    auto g_xyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 487);

    auto g_xyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 488);

    auto g_xyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 489);

    auto g_xyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 490);

    auto g_xyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 491);

    auto g_xyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 492);

    auto g_xyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 493);

    auto g_xyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 494);

    auto g_xyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 544);

    auto g_xyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 547);

    auto g_xyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 548);

    auto g_xyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 551);

    auto g_xyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 552);

    auto g_xyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 553);

    auto g_xyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 556);

    auto g_xyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 557);

    auto g_xyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 558);

    auto g_xyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 559);

    auto g_xyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 562);

    auto g_xyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 563);

    auto g_xyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 564);

    auto g_xyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 565);

    auto g_xyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 566);

    auto g_xyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 569);

    auto g_xyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 570);

    auto g_xyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 571);

    auto g_xyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 572);

    auto g_xyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 573);

    auto g_xyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 574);

    auto g_xyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 576);

    auto g_xyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 577);

    auto g_xyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 578);

    auto g_xyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 579);

    auto g_xyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 580);

    auto g_xyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 581);

    auto g_xyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 582);

    auto g_xyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 583);

    auto g_xyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 584);

    auto g_xzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 630);

    auto g_xzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 632);

    auto g_xzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 634);

    auto g_xzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 635);

    auto g_xzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 637);

    auto g_xzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 638);

    auto g_xzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 639);

    auto g_xzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 641);

    auto g_xzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 642);

    auto g_xzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 643);

    auto g_xzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 644);

    auto g_xzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 646);

    auto g_xzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 647);

    auto g_xzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 648);

    auto g_xzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 649);

    auto g_xzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 650);

    auto g_xzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 652);

    auto g_xzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 653);

    auto g_xzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 654);

    auto g_xzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 655);

    auto g_xzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 656);

    auto g_xzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 657);

    auto g_xzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 659);

    auto g_xzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 660);

    auto g_xzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 661);

    auto g_xzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 662);

    auto g_xzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 663);

    auto g_xzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 664);

    auto g_xzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 665);

    auto g_xzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 666);

    auto g_xzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 667);

    auto g_xzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 668);

    auto g_xzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 669);

    auto g_xzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 670);

    auto g_xzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 671);

    auto g_xzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 672);

    auto g_xzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 673);

    auto g_xzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 674);

    auto g_yyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 675);

    auto g_yyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 676);

    auto g_yyyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 677);

    auto g_yyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 678);

    auto g_yyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 679);

    auto g_yyyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 680);

    auto g_yyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 681);

    auto g_yyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 682);

    auto g_yyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 683);

    auto g_yyyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 684);

    auto g_yyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 685);

    auto g_yyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 686);

    auto g_yyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 687);

    auto g_yyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 688);

    auto g_yyyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 689);

    auto g_yyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 690);

    auto g_yyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 691);

    auto g_yyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 692);

    auto g_yyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 693);

    auto g_yyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 694);

    auto g_yyyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 695);

    auto g_yyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 696);

    auto g_yyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 697);

    auto g_yyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 698);

    auto g_yyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 699);

    auto g_yyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 700);

    auto g_yyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 701);

    auto g_yyyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 702);

    auto g_yyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 703);

    auto g_yyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 704);

    auto g_yyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 705);

    auto g_yyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 706);

    auto g_yyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 707);

    auto g_yyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 708);

    auto g_yyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 709);

    auto g_yyyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 710);

    auto g_yyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 711);

    auto g_yyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 712);

    auto g_yyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 713);

    auto g_yyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 714);

    auto g_yyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 715);

    auto g_yyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 716);

    auto g_yyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 717);

    auto g_yyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 718);

    auto g_yyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 719);

    auto g_yyyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 721);

    auto g_yyyyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 722);

    auto g_yyyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 723);

    auto g_yyyyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 724);

    auto g_yyyyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 725);

    auto g_yyyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 726);

    auto g_yyyyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 727);

    auto g_yyyyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 728);

    auto g_yyyyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 729);

    auto g_yyyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 730);

    auto g_yyyyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 731);

    auto g_yyyyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 732);

    auto g_yyyyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 733);

    auto g_yyyyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 734);

    auto g_yyyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 735);

    auto g_yyyyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 736);

    auto g_yyyyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 737);

    auto g_yyyyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 738);

    auto g_yyyyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 739);

    auto g_yyyyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 740);

    auto g_yyyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 741);

    auto g_yyyyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 742);

    auto g_yyyyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 743);

    auto g_yyyyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 744);

    auto g_yyyyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 745);

    auto g_yyyyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 746);

    auto g_yyyyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 747);

    auto g_yyyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 748);

    auto g_yyyyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 749);

    auto g_yyyyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 750);

    auto g_yyyyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 751);

    auto g_yyyyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 752);

    auto g_yyyyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 753);

    auto g_yyyyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 754);

    auto g_yyyyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 755);

    auto g_yyyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 756);

    auto g_yyyyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 757);

    auto g_yyyyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 758);

    auto g_yyyyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 759);

    auto g_yyyyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 760);

    auto g_yyyyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 761);

    auto g_yyyyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 762);

    auto g_yyyyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 763);

    auto g_yyyyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 764);

    auto g_yyyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 765);

    auto g_yyyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 766);

    auto g_yyyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 767);

    auto g_yyyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 768);

    auto g_yyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 769);

    auto g_yyyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 770);

    auto g_yyyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 771);

    auto g_yyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 772);

    auto g_yyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 773);

    auto g_yyyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 774);

    auto g_yyyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 775);

    auto g_yyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 776);

    auto g_yyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 777);

    auto g_yyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 778);

    auto g_yyyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 779);

    auto g_yyyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 780);

    auto g_yyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 781);

    auto g_yyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 782);

    auto g_yyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 783);

    auto g_yyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 784);

    auto g_yyyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 785);

    auto g_yyyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 786);

    auto g_yyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 787);

    auto g_yyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 788);

    auto g_yyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 789);

    auto g_yyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 790);

    auto g_yyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 791);

    auto g_yyyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 792);

    auto g_yyyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 793);

    auto g_yyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 794);

    auto g_yyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 795);

    auto g_yyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 796);

    auto g_yyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 797);

    auto g_yyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 798);

    auto g_yyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 799);

    auto g_yyyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 800);

    auto g_yyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 801);

    auto g_yyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 802);

    auto g_yyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 803);

    auto g_yyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 804);

    auto g_yyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 805);

    auto g_yyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 806);

    auto g_yyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 807);

    auto g_yyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 808);

    auto g_yyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 809);

    auto g_yyzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 810);

    auto g_yyzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 811);

    auto g_yyzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 812);

    auto g_yyzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 813);

    auto g_yyzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 814);

    auto g_yyzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 815);

    auto g_yyzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 816);

    auto g_yyzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 817);

    auto g_yyzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 818);

    auto g_yyzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 819);

    auto g_yyzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 820);

    auto g_yyzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 821);

    auto g_yyzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 822);

    auto g_yyzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 823);

    auto g_yyzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 824);

    auto g_yyzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 825);

    auto g_yyzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 826);

    auto g_yyzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 827);

    auto g_yyzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 828);

    auto g_yyzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 829);

    auto g_yyzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 830);

    auto g_yyzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 831);

    auto g_yyzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 832);

    auto g_yyzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 833);

    auto g_yyzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 834);

    auto g_yyzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 835);

    auto g_yyzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 836);

    auto g_yyzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 837);

    auto g_yyzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 838);

    auto g_yyzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 839);

    auto g_yyzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 840);

    auto g_yyzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 841);

    auto g_yyzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 842);

    auto g_yyzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 843);

    auto g_yyzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 844);

    auto g_yyzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 845);

    auto g_yyzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 846);

    auto g_yyzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 847);

    auto g_yyzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 848);

    auto g_yyzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 849);

    auto g_yyzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 850);

    auto g_yyzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 851);

    auto g_yyzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 852);

    auto g_yyzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 853);

    auto g_yyzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 854);

    auto g_yzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 855);

    auto g_yzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 856);

    auto g_yzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 857);

    auto g_yzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 858);

    auto g_yzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 859);

    auto g_yzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 860);

    auto g_yzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 861);

    auto g_yzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 862);

    auto g_yzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 863);

    auto g_yzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 864);

    auto g_yzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 865);

    auto g_yzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 866);

    auto g_yzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 867);

    auto g_yzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 868);

    auto g_yzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 869);

    auto g_yzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 870);

    auto g_yzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 871);

    auto g_yzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 872);

    auto g_yzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 873);

    auto g_yzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 874);

    auto g_yzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 875);

    auto g_yzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 876);

    auto g_yzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 877);

    auto g_yzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 878);

    auto g_yzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 879);

    auto g_yzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 880);

    auto g_yzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 881);

    auto g_yzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 882);

    auto g_yzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 883);

    auto g_yzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 884);

    auto g_yzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 885);

    auto g_yzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 886);

    auto g_yzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 887);

    auto g_yzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 888);

    auto g_yzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 889);

    auto g_yzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 890);

    auto g_yzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 891);

    auto g_yzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 892);

    auto g_yzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 893);

    auto g_yzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 894);

    auto g_yzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 895);

    auto g_yzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 896);

    auto g_yzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 897);

    auto g_yzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 898);

    auto g_yzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 899);

    auto g_zzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 900);

    auto g_zzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 901);

    auto g_zzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 902);

    auto g_zzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 903);

    auto g_zzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 904);

    auto g_zzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 905);

    auto g_zzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 906);

    auto g_zzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 907);

    auto g_zzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 908);

    auto g_zzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 909);

    auto g_zzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 910);

    auto g_zzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 911);

    auto g_zzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 912);

    auto g_zzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 913);

    auto g_zzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 914);

    auto g_zzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 915);

    auto g_zzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 916);

    auto g_zzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 917);

    auto g_zzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 918);

    auto g_zzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 919);

    auto g_zzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 920);

    auto g_zzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 921);

    auto g_zzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 922);

    auto g_zzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 923);

    auto g_zzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 924);

    auto g_zzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 925);

    auto g_zzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 926);

    auto g_zzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 927);

    auto g_zzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 928);

    auto g_zzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 929);

    auto g_zzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 930);

    auto g_zzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 931);

    auto g_zzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 932);

    auto g_zzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 933);

    auto g_zzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 934);

    auto g_zzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 935);

    auto g_zzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 936);

    auto g_zzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 937);

    auto g_zzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 938);

    auto g_zzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 939);

    auto g_zzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 940);

    auto g_zzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 941);

    auto g_zzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 942);

    auto g_zzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 943);

    auto g_zzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 944);

    /// Set up 0-45 components of targeted buffer : ISL

    auto g_xxxxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl);

    auto g_xxxxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 1);

    auto g_xxxxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 2);

    auto g_xxxxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 3);

    auto g_xxxxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 4);

    auto g_xxxxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 5);

    auto g_xxxxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 6);

    auto g_xxxxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 7);

    auto g_xxxxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 8);

    auto g_xxxxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 9);

    auto g_xxxxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 10);

    auto g_xxxxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 11);

    auto g_xxxxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 12);

    auto g_xxxxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 13);

    auto g_xxxxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 14);

    auto g_xxxxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 15);

    auto g_xxxxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 16);

    auto g_xxxxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 17);

    auto g_xxxxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 18);

    auto g_xxxxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 19);

    auto g_xxxxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 20);

    auto g_xxxxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 21);

    auto g_xxxxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 22);

    auto g_xxxxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 23);

    auto g_xxxxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 24);

    auto g_xxxxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 25);

    auto g_xxxxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 26);

    auto g_xxxxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 27);

    auto g_xxxxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 28);

    auto g_xxxxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 29);

    auto g_xxxxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 30);

    auto g_xxxxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 31);

    auto g_xxxxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 32);

    auto g_xxxxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 33);

    auto g_xxxxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 34);

    auto g_xxxxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 35);

    auto g_xxxxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 36);

    auto g_xxxxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 37);

    auto g_xxxxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 38);

    auto g_xxxxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 39);

    auto g_xxxxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 40);

    auto g_xxxxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 41);

    auto g_xxxxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 42);

    auto g_xxxxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 43);

    auto g_xxxxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 44);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxxx_0, g_xxxx_0_xxxxxxxx_1, g_xxxx_0_xxxxxxxy_0, g_xxxx_0_xxxxxxxy_1, g_xxxx_0_xxxxxxxz_0, g_xxxx_0_xxxxxxxz_1, g_xxxx_0_xxxxxxyy_0, g_xxxx_0_xxxxxxyy_1, g_xxxx_0_xxxxxxyz_0, g_xxxx_0_xxxxxxyz_1, g_xxxx_0_xxxxxxzz_0, g_xxxx_0_xxxxxxzz_1, g_xxxx_0_xxxxxyyy_0, g_xxxx_0_xxxxxyyy_1, g_xxxx_0_xxxxxyyz_0, g_xxxx_0_xxxxxyyz_1, g_xxxx_0_xxxxxyzz_0, g_xxxx_0_xxxxxyzz_1, g_xxxx_0_xxxxxzzz_0, g_xxxx_0_xxxxxzzz_1, g_xxxx_0_xxxxyyyy_0, g_xxxx_0_xxxxyyyy_1, g_xxxx_0_xxxxyyyz_0, g_xxxx_0_xxxxyyyz_1, g_xxxx_0_xxxxyyzz_0, g_xxxx_0_xxxxyyzz_1, g_xxxx_0_xxxxyzzz_0, g_xxxx_0_xxxxyzzz_1, g_xxxx_0_xxxxzzzz_0, g_xxxx_0_xxxxzzzz_1, g_xxxx_0_xxxyyyyy_0, g_xxxx_0_xxxyyyyy_1, g_xxxx_0_xxxyyyyz_0, g_xxxx_0_xxxyyyyz_1, g_xxxx_0_xxxyyyzz_0, g_xxxx_0_xxxyyyzz_1, g_xxxx_0_xxxyyzzz_0, g_xxxx_0_xxxyyzzz_1, g_xxxx_0_xxxyzzzz_0, g_xxxx_0_xxxyzzzz_1, g_xxxx_0_xxxzzzzz_0, g_xxxx_0_xxxzzzzz_1, g_xxxx_0_xxyyyyyy_0, g_xxxx_0_xxyyyyyy_1, g_xxxx_0_xxyyyyyz_0, g_xxxx_0_xxyyyyyz_1, g_xxxx_0_xxyyyyzz_0, g_xxxx_0_xxyyyyzz_1, g_xxxx_0_xxyyyzzz_0, g_xxxx_0_xxyyyzzz_1, g_xxxx_0_xxyyzzzz_0, g_xxxx_0_xxyyzzzz_1, g_xxxx_0_xxyzzzzz_0, g_xxxx_0_xxyzzzzz_1, g_xxxx_0_xxzzzzzz_0, g_xxxx_0_xxzzzzzz_1, g_xxxx_0_xyyyyyyy_0, g_xxxx_0_xyyyyyyy_1, g_xxxx_0_xyyyyyyz_0, g_xxxx_0_xyyyyyyz_1, g_xxxx_0_xyyyyyzz_0, g_xxxx_0_xyyyyyzz_1, g_xxxx_0_xyyyyzzz_0, g_xxxx_0_xyyyyzzz_1, g_xxxx_0_xyyyzzzz_0, g_xxxx_0_xyyyzzzz_1, g_xxxx_0_xyyzzzzz_0, g_xxxx_0_xyyzzzzz_1, g_xxxx_0_xyzzzzzz_0, g_xxxx_0_xyzzzzzz_1, g_xxxx_0_xzzzzzzz_0, g_xxxx_0_xzzzzzzz_1, g_xxxx_0_yyyyyyyy_0, g_xxxx_0_yyyyyyyy_1, g_xxxx_0_yyyyyyyz_0, g_xxxx_0_yyyyyyyz_1, g_xxxx_0_yyyyyyzz_0, g_xxxx_0_yyyyyyzz_1, g_xxxx_0_yyyyyzzz_0, g_xxxx_0_yyyyyzzz_1, g_xxxx_0_yyyyzzzz_0, g_xxxx_0_yyyyzzzz_1, g_xxxx_0_yyyzzzzz_0, g_xxxx_0_yyyzzzzz_1, g_xxxx_0_yyzzzzzz_0, g_xxxx_0_yyzzzzzz_1, g_xxxx_0_yzzzzzzz_0, g_xxxx_0_yzzzzzzz_1, g_xxxx_0_zzzzzzzz_0, g_xxxx_0_zzzzzzzz_1, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxxx_1, g_xxxxx_0_xxxxxxxy_1, g_xxxxx_0_xxxxxxxz_1, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxyy_1, g_xxxxx_0_xxxxxxyz_1, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxxzz_1, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyyy_1, g_xxxxx_0_xxxxxyyz_1, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxyzz_1, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxxzzz_1, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyyy_1, g_xxxxx_0_xxxxyyyz_1, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyyzz_1, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxyzzz_1, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxxzzzz_1, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyyy_1, g_xxxxx_0_xxxyyyyz_1, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyyzz_1, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyyzzz_1, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxyzzzz_1, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxxzzzzz_1, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyyy_1, g_xxxxx_0_xxyyyyyz_1, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyyzz_1, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyyzzz_1, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyyzzzz_1, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxyzzzzz_1, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xxzzzzzz_1, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyyy_1, g_xxxxx_0_xyyyyyyz_1, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyyzz_1, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyyzzz_1, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyyzzzz_1, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyyzzzzz_1, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xyzzzzzz_1, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_xzzzzzzz_1, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyyy_1, g_xxxxx_0_yyyyyyyz_1, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyyzz_1, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyyzzz_1, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyyzzzz_1, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyyzzzzz_1, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yyzzzzzz_1, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_yzzzzzzz_1, g_xxxxx_0_zzzzzzz_1, g_xxxxx_0_zzzzzzzz_1, g_xxxxxx_0_xxxxxxxx_0, g_xxxxxx_0_xxxxxxxy_0, g_xxxxxx_0_xxxxxxxz_0, g_xxxxxx_0_xxxxxxyy_0, g_xxxxxx_0_xxxxxxyz_0, g_xxxxxx_0_xxxxxxzz_0, g_xxxxxx_0_xxxxxyyy_0, g_xxxxxx_0_xxxxxyyz_0, g_xxxxxx_0_xxxxxyzz_0, g_xxxxxx_0_xxxxxzzz_0, g_xxxxxx_0_xxxxyyyy_0, g_xxxxxx_0_xxxxyyyz_0, g_xxxxxx_0_xxxxyyzz_0, g_xxxxxx_0_xxxxyzzz_0, g_xxxxxx_0_xxxxzzzz_0, g_xxxxxx_0_xxxyyyyy_0, g_xxxxxx_0_xxxyyyyz_0, g_xxxxxx_0_xxxyyyzz_0, g_xxxxxx_0_xxxyyzzz_0, g_xxxxxx_0_xxxyzzzz_0, g_xxxxxx_0_xxxzzzzz_0, g_xxxxxx_0_xxyyyyyy_0, g_xxxxxx_0_xxyyyyyz_0, g_xxxxxx_0_xxyyyyzz_0, g_xxxxxx_0_xxyyyzzz_0, g_xxxxxx_0_xxyyzzzz_0, g_xxxxxx_0_xxyzzzzz_0, g_xxxxxx_0_xxzzzzzz_0, g_xxxxxx_0_xyyyyyyy_0, g_xxxxxx_0_xyyyyyyz_0, g_xxxxxx_0_xyyyyyzz_0, g_xxxxxx_0_xyyyyzzz_0, g_xxxxxx_0_xyyyzzzz_0, g_xxxxxx_0_xyyzzzzz_0, g_xxxxxx_0_xyzzzzzz_0, g_xxxxxx_0_xzzzzzzz_0, g_xxxxxx_0_yyyyyyyy_0, g_xxxxxx_0_yyyyyyyz_0, g_xxxxxx_0_yyyyyyzz_0, g_xxxxxx_0_yyyyyzzz_0, g_xxxxxx_0_yyyyzzzz_0, g_xxxxxx_0_yyyzzzzz_0, g_xxxxxx_0_yyzzzzzz_0, g_xxxxxx_0_yzzzzzzz_0, g_xxxxxx_0_zzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xxxxxxxx_0[i] = 5.0 * g_xxxx_0_xxxxxxxx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxxx_1[i] * fz_be_0 + 8.0 * g_xxxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxxx_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxxy_0[i] = 5.0 * g_xxxx_0_xxxxxxxy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxxz_0[i] = 5.0 * g_xxxx_0_xxxxxxxz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxyy_0[i] = 5.0 * g_xxxx_0_xxxxxxyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxyz_0[i] = 5.0 * g_xxxx_0_xxxxxxyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxzz_0[i] = 5.0 * g_xxxx_0_xxxxxxzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxyyy_0[i] = 5.0 * g_xxxx_0_xxxxxyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxyyz_0[i] = 5.0 * g_xxxx_0_xxxxxyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxyzz_0[i] = 5.0 * g_xxxx_0_xxxxxyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxzzz_0[i] = 5.0 * g_xxxx_0_xxxxxzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyyyy_0[i] = 5.0 * g_xxxx_0_xxxxyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyyyz_0[i] = 5.0 * g_xxxx_0_xxxxyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyyzz_0[i] = 5.0 * g_xxxx_0_xxxxyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyzzz_0[i] = 5.0 * g_xxxx_0_xxxxyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxzzzz_0[i] = 5.0 * g_xxxx_0_xxxxzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyyyy_0[i] = 5.0 * g_xxxx_0_xxxyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyyyz_0[i] = 5.0 * g_xxxx_0_xxxyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyyzz_0[i] = 5.0 * g_xxxx_0_xxxyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyzzz_0[i] = 5.0 * g_xxxx_0_xxxyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyzzzz_0[i] = 5.0 * g_xxxx_0_xxxyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxzzzzz_0[i] = 5.0 * g_xxxx_0_xxxzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyyyy_0[i] = 5.0 * g_xxxx_0_xxyyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyyyz_0[i] = 5.0 * g_xxxx_0_xxyyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyyzz_0[i] = 5.0 * g_xxxx_0_xxyyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyzzz_0[i] = 5.0 * g_xxxx_0_xxyyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyzzzz_0[i] = 5.0 * g_xxxx_0_xxyyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyzzzzz_0[i] = 5.0 * g_xxxx_0_xxyzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxzzzzzz_0[i] = 5.0 * g_xxxx_0_xxzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyyyy_0[i] = 5.0 * g_xxxx_0_xyyyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyyyz_0[i] = 5.0 * g_xxxx_0_xyyyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyyzz_0[i] = 5.0 * g_xxxx_0_xyyyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyzzz_0[i] = 5.0 * g_xxxx_0_xyyyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyzzzz_0[i] = 5.0 * g_xxxx_0_xyyyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyzzzzz_0[i] = 5.0 * g_xxxx_0_xyyzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyzzzzzz_0[i] = 5.0 * g_xxxx_0_xyzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xzzzzzzz_0[i] = 5.0 * g_xxxx_0_xzzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyyyy_0[i] = 5.0 * g_xxxx_0_yyyyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyyyz_0[i] = 5.0 * g_xxxx_0_yyyyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyyzz_0[i] = 5.0 * g_xxxx_0_yyyyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyzzz_0[i] = 5.0 * g_xxxx_0_yyyyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyzzzz_0[i] = 5.0 * g_xxxx_0_yyyyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyzzzzz_0[i] = 5.0 * g_xxxx_0_yyyzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyzzzzzz_0[i] = 5.0 * g_xxxx_0_yyzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yzzzzzzz_0[i] = 5.0 * g_xxxx_0_yzzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_zzzzzzzz_0[i] = 5.0 * g_xxxx_0_zzzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : ISL

    auto g_xxxxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 45);

    auto g_xxxxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 46);

    auto g_xxxxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 47);

    auto g_xxxxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 48);

    auto g_xxxxxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 49);

    auto g_xxxxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 50);

    auto g_xxxxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 51);

    auto g_xxxxxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 52);

    auto g_xxxxxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 53);

    auto g_xxxxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 54);

    auto g_xxxxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 55);

    auto g_xxxxxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 56);

    auto g_xxxxxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 57);

    auto g_xxxxxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 58);

    auto g_xxxxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 59);

    auto g_xxxxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 60);

    auto g_xxxxxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 61);

    auto g_xxxxxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 62);

    auto g_xxxxxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 63);

    auto g_xxxxxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 64);

    auto g_xxxxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 65);

    auto g_xxxxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 66);

    auto g_xxxxxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 67);

    auto g_xxxxxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 68);

    auto g_xxxxxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 69);

    auto g_xxxxxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 70);

    auto g_xxxxxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 71);

    auto g_xxxxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 72);

    auto g_xxxxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 73);

    auto g_xxxxxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 74);

    auto g_xxxxxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 75);

    auto g_xxxxxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 76);

    auto g_xxxxxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 77);

    auto g_xxxxxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 78);

    auto g_xxxxxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 79);

    auto g_xxxxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 80);

    auto g_xxxxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 81);

    auto g_xxxxxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 82);

    auto g_xxxxxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 83);

    auto g_xxxxxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 84);

    auto g_xxxxxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 85);

    auto g_xxxxxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 86);

    auto g_xxxxxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 87);

    auto g_xxxxxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 88);

    auto g_xxxxxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 89);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxxx_1, g_xxxxx_0_xxxxxxxy_1, g_xxxxx_0_xxxxxxxz_1, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxyy_1, g_xxxxx_0_xxxxxxyz_1, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxxzz_1, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyyy_1, g_xxxxx_0_xxxxxyyz_1, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxyzz_1, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxxzzz_1, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyyy_1, g_xxxxx_0_xxxxyyyz_1, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyyzz_1, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxyzzz_1, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxxzzzz_1, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyyy_1, g_xxxxx_0_xxxyyyyz_1, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyyzz_1, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyyzzz_1, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxyzzzz_1, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxxzzzzz_1, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyyy_1, g_xxxxx_0_xxyyyyyz_1, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyyzz_1, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyyzzz_1, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyyzzzz_1, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxyzzzzz_1, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xxzzzzzz_1, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyyy_1, g_xxxxx_0_xyyyyyyz_1, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyyzz_1, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyyzzz_1, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyyzzzz_1, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyyzzzzz_1, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xyzzzzzz_1, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_xzzzzzzz_1, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyyy_1, g_xxxxx_0_yyyyyyyz_1, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyyzz_1, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyyzzz_1, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyyzzzz_1, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyyzzzzz_1, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yyzzzzzz_1, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_yzzzzzzz_1, g_xxxxx_0_zzzzzzz_1, g_xxxxx_0_zzzzzzzz_1, g_xxxxxy_0_xxxxxxxx_0, g_xxxxxy_0_xxxxxxxy_0, g_xxxxxy_0_xxxxxxxz_0, g_xxxxxy_0_xxxxxxyy_0, g_xxxxxy_0_xxxxxxyz_0, g_xxxxxy_0_xxxxxxzz_0, g_xxxxxy_0_xxxxxyyy_0, g_xxxxxy_0_xxxxxyyz_0, g_xxxxxy_0_xxxxxyzz_0, g_xxxxxy_0_xxxxxzzz_0, g_xxxxxy_0_xxxxyyyy_0, g_xxxxxy_0_xxxxyyyz_0, g_xxxxxy_0_xxxxyyzz_0, g_xxxxxy_0_xxxxyzzz_0, g_xxxxxy_0_xxxxzzzz_0, g_xxxxxy_0_xxxyyyyy_0, g_xxxxxy_0_xxxyyyyz_0, g_xxxxxy_0_xxxyyyzz_0, g_xxxxxy_0_xxxyyzzz_0, g_xxxxxy_0_xxxyzzzz_0, g_xxxxxy_0_xxxzzzzz_0, g_xxxxxy_0_xxyyyyyy_0, g_xxxxxy_0_xxyyyyyz_0, g_xxxxxy_0_xxyyyyzz_0, g_xxxxxy_0_xxyyyzzz_0, g_xxxxxy_0_xxyyzzzz_0, g_xxxxxy_0_xxyzzzzz_0, g_xxxxxy_0_xxzzzzzz_0, g_xxxxxy_0_xyyyyyyy_0, g_xxxxxy_0_xyyyyyyz_0, g_xxxxxy_0_xyyyyyzz_0, g_xxxxxy_0_xyyyyzzz_0, g_xxxxxy_0_xyyyzzzz_0, g_xxxxxy_0_xyyzzzzz_0, g_xxxxxy_0_xyzzzzzz_0, g_xxxxxy_0_xzzzzzzz_0, g_xxxxxy_0_yyyyyyyy_0, g_xxxxxy_0_yyyyyyyz_0, g_xxxxxy_0_yyyyyyzz_0, g_xxxxxy_0_yyyyyzzz_0, g_xxxxxy_0_yyyyzzzz_0, g_xxxxxy_0_yyyzzzzz_0, g_xxxxxy_0_yyzzzzzz_0, g_xxxxxy_0_yzzzzzzz_0, g_xxxxxy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xxxxxxxx_0[i] = g_xxxxx_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxxy_0[i] = g_xxxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxxz_0[i] = g_xxxxx_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxyy_0[i] = 2.0 * g_xxxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxyz_0[i] = g_xxxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxzz_0[i] = g_xxxxx_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxyyy_0[i] = 3.0 * g_xxxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxyyz_0[i] = 2.0 * g_xxxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxyzz_0[i] = g_xxxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxzzz_0[i] = g_xxxxx_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyyyy_0[i] = 4.0 * g_xxxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyyyz_0[i] = 3.0 * g_xxxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyyzz_0[i] = 2.0 * g_xxxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyzzz_0[i] = g_xxxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxzzzz_0[i] = g_xxxxx_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyyyy_0[i] = 5.0 * g_xxxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyyyz_0[i] = 4.0 * g_xxxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyyzz_0[i] = 3.0 * g_xxxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyzzz_0[i] = 2.0 * g_xxxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyzzzz_0[i] = g_xxxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxzzzzz_0[i] = g_xxxxx_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyyyy_0[i] = 6.0 * g_xxxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyyyz_0[i] = 5.0 * g_xxxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyyzz_0[i] = 4.0 * g_xxxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyzzz_0[i] = 3.0 * g_xxxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyzzzz_0[i] = 2.0 * g_xxxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyzzzzz_0[i] = g_xxxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxzzzzzz_0[i] = g_xxxxx_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyyyy_0[i] = 7.0 * g_xxxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyyyz_0[i] = 6.0 * g_xxxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyyzz_0[i] = 5.0 * g_xxxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyzzz_0[i] = 4.0 * g_xxxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyzzzz_0[i] = 3.0 * g_xxxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyzzzzz_0[i] = 2.0 * g_xxxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyzzzzzz_0[i] = g_xxxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xzzzzzzz_0[i] = g_xxxxx_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyyyy_0[i] = 8.0 * g_xxxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyyyz_0[i] = 7.0 * g_xxxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyyzz_0[i] = 6.0 * g_xxxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyzzz_0[i] = 5.0 * g_xxxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyzzzz_0[i] = 4.0 * g_xxxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyzzzzz_0[i] = 3.0 * g_xxxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyzzzzzz_0[i] = 2.0 * g_xxxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yzzzzzzz_0[i] = g_xxxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_zzzzzzzz_0[i] = g_xxxxx_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 90-135 components of targeted buffer : ISL

    auto g_xxxxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 90);

    auto g_xxxxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 91);

    auto g_xxxxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 92);

    auto g_xxxxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 93);

    auto g_xxxxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 94);

    auto g_xxxxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 95);

    auto g_xxxxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 96);

    auto g_xxxxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 97);

    auto g_xxxxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 98);

    auto g_xxxxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 99);

    auto g_xxxxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 100);

    auto g_xxxxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 101);

    auto g_xxxxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 102);

    auto g_xxxxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 103);

    auto g_xxxxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 104);

    auto g_xxxxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 105);

    auto g_xxxxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 106);

    auto g_xxxxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 107);

    auto g_xxxxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 108);

    auto g_xxxxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 109);

    auto g_xxxxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 110);

    auto g_xxxxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 111);

    auto g_xxxxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 112);

    auto g_xxxxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 113);

    auto g_xxxxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 114);

    auto g_xxxxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 115);

    auto g_xxxxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 116);

    auto g_xxxxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 117);

    auto g_xxxxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 118);

    auto g_xxxxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 119);

    auto g_xxxxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 120);

    auto g_xxxxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 121);

    auto g_xxxxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 122);

    auto g_xxxxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 123);

    auto g_xxxxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 124);

    auto g_xxxxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 125);

    auto g_xxxxxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 126);

    auto g_xxxxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 127);

    auto g_xxxxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 128);

    auto g_xxxxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 129);

    auto g_xxxxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 130);

    auto g_xxxxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 131);

    auto g_xxxxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 132);

    auto g_xxxxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 133);

    auto g_xxxxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 134);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxxx_1, g_xxxxx_0_xxxxxxxy_1, g_xxxxx_0_xxxxxxxz_1, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxyy_1, g_xxxxx_0_xxxxxxyz_1, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxxzz_1, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyyy_1, g_xxxxx_0_xxxxxyyz_1, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxyzz_1, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxxzzz_1, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyyy_1, g_xxxxx_0_xxxxyyyz_1, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyyzz_1, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxyzzz_1, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxxzzzz_1, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyyy_1, g_xxxxx_0_xxxyyyyz_1, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyyzz_1, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyyzzz_1, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxyzzzz_1, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxxzzzzz_1, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyyy_1, g_xxxxx_0_xxyyyyyz_1, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyyzz_1, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyyzzz_1, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyyzzzz_1, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxyzzzzz_1, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xxzzzzzz_1, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyyy_1, g_xxxxx_0_xyyyyyyz_1, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyyzz_1, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyyzzz_1, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyyzzzz_1, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyyzzzzz_1, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xyzzzzzz_1, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_xzzzzzzz_1, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyyy_1, g_xxxxx_0_yyyyyyyz_1, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyyzz_1, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyyzzz_1, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyyzzzz_1, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyyzzzzz_1, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yyzzzzzz_1, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_yzzzzzzz_1, g_xxxxx_0_zzzzzzz_1, g_xxxxx_0_zzzzzzzz_1, g_xxxxxz_0_xxxxxxxx_0, g_xxxxxz_0_xxxxxxxy_0, g_xxxxxz_0_xxxxxxxz_0, g_xxxxxz_0_xxxxxxyy_0, g_xxxxxz_0_xxxxxxyz_0, g_xxxxxz_0_xxxxxxzz_0, g_xxxxxz_0_xxxxxyyy_0, g_xxxxxz_0_xxxxxyyz_0, g_xxxxxz_0_xxxxxyzz_0, g_xxxxxz_0_xxxxxzzz_0, g_xxxxxz_0_xxxxyyyy_0, g_xxxxxz_0_xxxxyyyz_0, g_xxxxxz_0_xxxxyyzz_0, g_xxxxxz_0_xxxxyzzz_0, g_xxxxxz_0_xxxxzzzz_0, g_xxxxxz_0_xxxyyyyy_0, g_xxxxxz_0_xxxyyyyz_0, g_xxxxxz_0_xxxyyyzz_0, g_xxxxxz_0_xxxyyzzz_0, g_xxxxxz_0_xxxyzzzz_0, g_xxxxxz_0_xxxzzzzz_0, g_xxxxxz_0_xxyyyyyy_0, g_xxxxxz_0_xxyyyyyz_0, g_xxxxxz_0_xxyyyyzz_0, g_xxxxxz_0_xxyyyzzz_0, g_xxxxxz_0_xxyyzzzz_0, g_xxxxxz_0_xxyzzzzz_0, g_xxxxxz_0_xxzzzzzz_0, g_xxxxxz_0_xyyyyyyy_0, g_xxxxxz_0_xyyyyyyz_0, g_xxxxxz_0_xyyyyyzz_0, g_xxxxxz_0_xyyyyzzz_0, g_xxxxxz_0_xyyyzzzz_0, g_xxxxxz_0_xyyzzzzz_0, g_xxxxxz_0_xyzzzzzz_0, g_xxxxxz_0_xzzzzzzz_0, g_xxxxxz_0_yyyyyyyy_0, g_xxxxxz_0_yyyyyyyz_0, g_xxxxxz_0_yyyyyyzz_0, g_xxxxxz_0_yyyyyzzz_0, g_xxxxxz_0_yyyyzzzz_0, g_xxxxxz_0_yyyzzzzz_0, g_xxxxxz_0_yyzzzzzz_0, g_xxxxxz_0_yzzzzzzz_0, g_xxxxxz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xxxxxxxx_0[i] = g_xxxxx_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxxy_0[i] = g_xxxxx_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxxz_0[i] = g_xxxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxyy_0[i] = g_xxxxx_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxyz_0[i] = g_xxxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxzz_0[i] = 2.0 * g_xxxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxyyy_0[i] = g_xxxxx_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxyyz_0[i] = g_xxxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxyzz_0[i] = 2.0 * g_xxxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxzzz_0[i] = 3.0 * g_xxxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyyyy_0[i] = g_xxxxx_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyyyz_0[i] = g_xxxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyyzz_0[i] = 2.0 * g_xxxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyzzz_0[i] = 3.0 * g_xxxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxzzzz_0[i] = 4.0 * g_xxxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyyyy_0[i] = g_xxxxx_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyyyz_0[i] = g_xxxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyyzz_0[i] = 2.0 * g_xxxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyzzz_0[i] = 3.0 * g_xxxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyzzzz_0[i] = 4.0 * g_xxxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxzzzzz_0[i] = 5.0 * g_xxxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyyyy_0[i] = g_xxxxx_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyyyz_0[i] = g_xxxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyyzz_0[i] = 2.0 * g_xxxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyzzz_0[i] = 3.0 * g_xxxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyzzzz_0[i] = 4.0 * g_xxxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyzzzzz_0[i] = 5.0 * g_xxxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxzzzzzz_0[i] = 6.0 * g_xxxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyyyy_0[i] = g_xxxxx_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyyyz_0[i] = g_xxxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyyzz_0[i] = 2.0 * g_xxxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyzzz_0[i] = 3.0 * g_xxxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyzzzz_0[i] = 4.0 * g_xxxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyzzzzz_0[i] = 5.0 * g_xxxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyzzzzzz_0[i] = 6.0 * g_xxxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xzzzzzzz_0[i] = 7.0 * g_xxxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyyyy_0[i] = g_xxxxx_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyyyz_0[i] = g_xxxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyyzz_0[i] = 2.0 * g_xxxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyzzz_0[i] = 3.0 * g_xxxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyzzzz_0[i] = 4.0 * g_xxxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyzzzzz_0[i] = 5.0 * g_xxxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyzzzzzz_0[i] = 6.0 * g_xxxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yzzzzzzz_0[i] = 7.0 * g_xxxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_zzzzzzzz_0[i] = 8.0 * g_xxxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 135-180 components of targeted buffer : ISL

    auto g_xxxxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 135);

    auto g_xxxxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 136);

    auto g_xxxxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 137);

    auto g_xxxxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 138);

    auto g_xxxxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 139);

    auto g_xxxxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 140);

    auto g_xxxxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 141);

    auto g_xxxxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 142);

    auto g_xxxxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 143);

    auto g_xxxxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 144);

    auto g_xxxxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 145);

    auto g_xxxxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 146);

    auto g_xxxxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 147);

    auto g_xxxxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 148);

    auto g_xxxxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 149);

    auto g_xxxxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 150);

    auto g_xxxxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 151);

    auto g_xxxxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 152);

    auto g_xxxxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 153);

    auto g_xxxxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 154);

    auto g_xxxxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 155);

    auto g_xxxxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 156);

    auto g_xxxxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 157);

    auto g_xxxxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 158);

    auto g_xxxxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 159);

    auto g_xxxxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 160);

    auto g_xxxxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 161);

    auto g_xxxxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 162);

    auto g_xxxxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 163);

    auto g_xxxxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 164);

    auto g_xxxxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 165);

    auto g_xxxxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 166);

    auto g_xxxxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 167);

    auto g_xxxxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 168);

    auto g_xxxxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 169);

    auto g_xxxxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 170);

    auto g_xxxxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 171);

    auto g_xxxxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 172);

    auto g_xxxxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 173);

    auto g_xxxxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 174);

    auto g_xxxxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 175);

    auto g_xxxxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 176);

    auto g_xxxxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 177);

    auto g_xxxxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 178);

    auto g_xxxxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 179);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxxx_0, g_xxxx_0_xxxxxxxx_1, g_xxxx_0_xxxxxxxz_0, g_xxxx_0_xxxxxxxz_1, g_xxxx_0_xxxxxxzz_0, g_xxxx_0_xxxxxxzz_1, g_xxxx_0_xxxxxzzz_0, g_xxxx_0_xxxxxzzz_1, g_xxxx_0_xxxxzzzz_0, g_xxxx_0_xxxxzzzz_1, g_xxxx_0_xxxzzzzz_0, g_xxxx_0_xxxzzzzz_1, g_xxxx_0_xxzzzzzz_0, g_xxxx_0_xxzzzzzz_1, g_xxxx_0_xzzzzzzz_0, g_xxxx_0_xzzzzzzz_1, g_xxxxy_0_xxxxxxxx_1, g_xxxxy_0_xxxxxxxz_1, g_xxxxy_0_xxxxxxzz_1, g_xxxxy_0_xxxxxzzz_1, g_xxxxy_0_xxxxzzzz_1, g_xxxxy_0_xxxzzzzz_1, g_xxxxy_0_xxzzzzzz_1, g_xxxxy_0_xzzzzzzz_1, g_xxxxyy_0_xxxxxxxx_0, g_xxxxyy_0_xxxxxxxy_0, g_xxxxyy_0_xxxxxxxz_0, g_xxxxyy_0_xxxxxxyy_0, g_xxxxyy_0_xxxxxxyz_0, g_xxxxyy_0_xxxxxxzz_0, g_xxxxyy_0_xxxxxyyy_0, g_xxxxyy_0_xxxxxyyz_0, g_xxxxyy_0_xxxxxyzz_0, g_xxxxyy_0_xxxxxzzz_0, g_xxxxyy_0_xxxxyyyy_0, g_xxxxyy_0_xxxxyyyz_0, g_xxxxyy_0_xxxxyyzz_0, g_xxxxyy_0_xxxxyzzz_0, g_xxxxyy_0_xxxxzzzz_0, g_xxxxyy_0_xxxyyyyy_0, g_xxxxyy_0_xxxyyyyz_0, g_xxxxyy_0_xxxyyyzz_0, g_xxxxyy_0_xxxyyzzz_0, g_xxxxyy_0_xxxyzzzz_0, g_xxxxyy_0_xxxzzzzz_0, g_xxxxyy_0_xxyyyyyy_0, g_xxxxyy_0_xxyyyyyz_0, g_xxxxyy_0_xxyyyyzz_0, g_xxxxyy_0_xxyyyzzz_0, g_xxxxyy_0_xxyyzzzz_0, g_xxxxyy_0_xxyzzzzz_0, g_xxxxyy_0_xxzzzzzz_0, g_xxxxyy_0_xyyyyyyy_0, g_xxxxyy_0_xyyyyyyz_0, g_xxxxyy_0_xyyyyyzz_0, g_xxxxyy_0_xyyyyzzz_0, g_xxxxyy_0_xyyyzzzz_0, g_xxxxyy_0_xyyzzzzz_0, g_xxxxyy_0_xyzzzzzz_0, g_xxxxyy_0_xzzzzzzz_0, g_xxxxyy_0_yyyyyyyy_0, g_xxxxyy_0_yyyyyyyz_0, g_xxxxyy_0_yyyyyyzz_0, g_xxxxyy_0_yyyyyzzz_0, g_xxxxyy_0_yyyyzzzz_0, g_xxxxyy_0_yyyzzzzz_0, g_xxxxyy_0_yyzzzzzz_0, g_xxxxyy_0_yzzzzzzz_0, g_xxxxyy_0_zzzzzzzz_0, g_xxxyy_0_xxxxxxxy_1, g_xxxyy_0_xxxxxxy_1, g_xxxyy_0_xxxxxxyy_1, g_xxxyy_0_xxxxxxyz_1, g_xxxyy_0_xxxxxyy_1, g_xxxyy_0_xxxxxyyy_1, g_xxxyy_0_xxxxxyyz_1, g_xxxyy_0_xxxxxyz_1, g_xxxyy_0_xxxxxyzz_1, g_xxxyy_0_xxxxyyy_1, g_xxxyy_0_xxxxyyyy_1, g_xxxyy_0_xxxxyyyz_1, g_xxxyy_0_xxxxyyz_1, g_xxxyy_0_xxxxyyzz_1, g_xxxyy_0_xxxxyzz_1, g_xxxyy_0_xxxxyzzz_1, g_xxxyy_0_xxxyyyy_1, g_xxxyy_0_xxxyyyyy_1, g_xxxyy_0_xxxyyyyz_1, g_xxxyy_0_xxxyyyz_1, g_xxxyy_0_xxxyyyzz_1, g_xxxyy_0_xxxyyzz_1, g_xxxyy_0_xxxyyzzz_1, g_xxxyy_0_xxxyzzz_1, g_xxxyy_0_xxxyzzzz_1, g_xxxyy_0_xxyyyyy_1, g_xxxyy_0_xxyyyyyy_1, g_xxxyy_0_xxyyyyyz_1, g_xxxyy_0_xxyyyyz_1, g_xxxyy_0_xxyyyyzz_1, g_xxxyy_0_xxyyyzz_1, g_xxxyy_0_xxyyyzzz_1, g_xxxyy_0_xxyyzzz_1, g_xxxyy_0_xxyyzzzz_1, g_xxxyy_0_xxyzzzz_1, g_xxxyy_0_xxyzzzzz_1, g_xxxyy_0_xyyyyyy_1, g_xxxyy_0_xyyyyyyy_1, g_xxxyy_0_xyyyyyyz_1, g_xxxyy_0_xyyyyyz_1, g_xxxyy_0_xyyyyyzz_1, g_xxxyy_0_xyyyyzz_1, g_xxxyy_0_xyyyyzzz_1, g_xxxyy_0_xyyyzzz_1, g_xxxyy_0_xyyyzzzz_1, g_xxxyy_0_xyyzzzz_1, g_xxxyy_0_xyyzzzzz_1, g_xxxyy_0_xyzzzzz_1, g_xxxyy_0_xyzzzzzz_1, g_xxxyy_0_yyyyyyy_1, g_xxxyy_0_yyyyyyyy_1, g_xxxyy_0_yyyyyyyz_1, g_xxxyy_0_yyyyyyz_1, g_xxxyy_0_yyyyyyzz_1, g_xxxyy_0_yyyyyzz_1, g_xxxyy_0_yyyyyzzz_1, g_xxxyy_0_yyyyzzz_1, g_xxxyy_0_yyyyzzzz_1, g_xxxyy_0_yyyzzzz_1, g_xxxyy_0_yyyzzzzz_1, g_xxxyy_0_yyzzzzz_1, g_xxxyy_0_yyzzzzzz_1, g_xxxyy_0_yzzzzzz_1, g_xxxyy_0_yzzzzzzz_1, g_xxxyy_0_zzzzzzzz_1, g_xxyy_0_xxxxxxxy_0, g_xxyy_0_xxxxxxxy_1, g_xxyy_0_xxxxxxyy_0, g_xxyy_0_xxxxxxyy_1, g_xxyy_0_xxxxxxyz_0, g_xxyy_0_xxxxxxyz_1, g_xxyy_0_xxxxxyyy_0, g_xxyy_0_xxxxxyyy_1, g_xxyy_0_xxxxxyyz_0, g_xxyy_0_xxxxxyyz_1, g_xxyy_0_xxxxxyzz_0, g_xxyy_0_xxxxxyzz_1, g_xxyy_0_xxxxyyyy_0, g_xxyy_0_xxxxyyyy_1, g_xxyy_0_xxxxyyyz_0, g_xxyy_0_xxxxyyyz_1, g_xxyy_0_xxxxyyzz_0, g_xxyy_0_xxxxyyzz_1, g_xxyy_0_xxxxyzzz_0, g_xxyy_0_xxxxyzzz_1, g_xxyy_0_xxxyyyyy_0, g_xxyy_0_xxxyyyyy_1, g_xxyy_0_xxxyyyyz_0, g_xxyy_0_xxxyyyyz_1, g_xxyy_0_xxxyyyzz_0, g_xxyy_0_xxxyyyzz_1, g_xxyy_0_xxxyyzzz_0, g_xxyy_0_xxxyyzzz_1, g_xxyy_0_xxxyzzzz_0, g_xxyy_0_xxxyzzzz_1, g_xxyy_0_xxyyyyyy_0, g_xxyy_0_xxyyyyyy_1, g_xxyy_0_xxyyyyyz_0, g_xxyy_0_xxyyyyyz_1, g_xxyy_0_xxyyyyzz_0, g_xxyy_0_xxyyyyzz_1, g_xxyy_0_xxyyyzzz_0, g_xxyy_0_xxyyyzzz_1, g_xxyy_0_xxyyzzzz_0, g_xxyy_0_xxyyzzzz_1, g_xxyy_0_xxyzzzzz_0, g_xxyy_0_xxyzzzzz_1, g_xxyy_0_xyyyyyyy_0, g_xxyy_0_xyyyyyyy_1, g_xxyy_0_xyyyyyyz_0, g_xxyy_0_xyyyyyyz_1, g_xxyy_0_xyyyyyzz_0, g_xxyy_0_xyyyyyzz_1, g_xxyy_0_xyyyyzzz_0, g_xxyy_0_xyyyyzzz_1, g_xxyy_0_xyyyzzzz_0, g_xxyy_0_xyyyzzzz_1, g_xxyy_0_xyyzzzzz_0, g_xxyy_0_xyyzzzzz_1, g_xxyy_0_xyzzzzzz_0, g_xxyy_0_xyzzzzzz_1, g_xxyy_0_yyyyyyyy_0, g_xxyy_0_yyyyyyyy_1, g_xxyy_0_yyyyyyyz_0, g_xxyy_0_yyyyyyyz_1, g_xxyy_0_yyyyyyzz_0, g_xxyy_0_yyyyyyzz_1, g_xxyy_0_yyyyyzzz_0, g_xxyy_0_yyyyyzzz_1, g_xxyy_0_yyyyzzzz_0, g_xxyy_0_yyyyzzzz_1, g_xxyy_0_yyyzzzzz_0, g_xxyy_0_yyyzzzzz_1, g_xxyy_0_yyzzzzzz_0, g_xxyy_0_yyzzzzzz_1, g_xxyy_0_yzzzzzzz_0, g_xxyy_0_yzzzzzzz_1, g_xxyy_0_zzzzzzzz_0, g_xxyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xxxxxxxx_0[i] = g_xxxx_0_xxxxxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxxy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxxxxy_0[i] = 3.0 * g_xxyy_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxxyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxxxz_0[i] = g_xxxx_0_xxxxxxxz_0[i] * fbe_0 - g_xxxx_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxxy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxxxyy_0[i] = 3.0 * g_xxyy_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxxyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxxyz_0[i] = 3.0 * g_xxyy_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxxzz_0[i] = g_xxxx_0_xxxxxxzz_0[i] * fbe_0 - g_xxxx_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxxy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxxyyy_0[i] = 3.0 * g_xxyy_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxxyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxyyz_0[i] = 3.0 * g_xxyy_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxyzz_0[i] = 3.0 * g_xxyy_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxzzz_0[i] = g_xxxx_0_xxxxxzzz_0[i] * fbe_0 - g_xxxx_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxxy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxyyyy_0[i] = 3.0 * g_xxyy_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxyyyz_0[i] = 3.0 * g_xxyy_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxyyzz_0[i] = 3.0 * g_xxyy_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxyzzz_0[i] = 3.0 * g_xxyy_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxzzzz_0[i] = g_xxxx_0_xxxxzzzz_0[i] * fbe_0 - g_xxxx_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxxy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxyyyyy_0[i] = 3.0 * g_xxyy_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyyyyz_0[i] = 3.0 * g_xxyy_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyyyzz_0[i] = 3.0 * g_xxyy_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyyzzz_0[i] = 3.0 * g_xxyy_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyzzzz_0[i] = 3.0 * g_xxyy_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxzzzzz_0[i] = g_xxxx_0_xxxzzzzz_0[i] * fbe_0 - g_xxxx_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxxy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxyyyyyy_0[i] = 3.0 * g_xxyy_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyyyyz_0[i] = 3.0 * g_xxyy_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyyyzz_0[i] = 3.0 * g_xxyy_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyyzzz_0[i] = 3.0 * g_xxyy_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyzzzz_0[i] = 3.0 * g_xxyy_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyzzzzz_0[i] = 3.0 * g_xxyy_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxzzzzzz_0[i] = g_xxxx_0_xxzzzzzz_0[i] * fbe_0 - g_xxxx_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxxy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xyyyyyyy_0[i] = 3.0 * g_xxyy_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyyyyz_0[i] = 3.0 * g_xxyy_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyyyzz_0[i] = 3.0 * g_xxyy_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyyzzz_0[i] = 3.0 * g_xxyy_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyzzzz_0[i] = 3.0 * g_xxyy_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyzzzzz_0[i] = 3.0 * g_xxyy_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyzzzzzz_0[i] = 3.0 * g_xxyy_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xzzzzzzz_0[i] = g_xxxx_0_xzzzzzzz_0[i] * fbe_0 - g_xxxx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxxy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_yyyyyyyy_0[i] = 3.0 * g_xxyy_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyyyyz_0[i] = 3.0 * g_xxyy_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyyyzz_0[i] = 3.0 * g_xxyy_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyyzzz_0[i] = 3.0 * g_xxyy_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyzzzz_0[i] = 3.0 * g_xxyy_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyzzzzz_0[i] = 3.0 * g_xxyy_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyzzzzzz_0[i] = 3.0 * g_xxyy_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yzzzzzzz_0[i] = 3.0 * g_xxyy_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_zzzzzzzz_0[i] = 3.0 * g_xxyy_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 180-225 components of targeted buffer : ISL

    auto g_xxxxyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 180);

    auto g_xxxxyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 181);

    auto g_xxxxyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 182);

    auto g_xxxxyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 183);

    auto g_xxxxyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 184);

    auto g_xxxxyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 185);

    auto g_xxxxyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 186);

    auto g_xxxxyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 187);

    auto g_xxxxyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 188);

    auto g_xxxxyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 189);

    auto g_xxxxyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 190);

    auto g_xxxxyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 191);

    auto g_xxxxyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 192);

    auto g_xxxxyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 193);

    auto g_xxxxyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 194);

    auto g_xxxxyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 195);

    auto g_xxxxyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 196);

    auto g_xxxxyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 197);

    auto g_xxxxyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 198);

    auto g_xxxxyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 199);

    auto g_xxxxyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 200);

    auto g_xxxxyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 201);

    auto g_xxxxyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 202);

    auto g_xxxxyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 203);

    auto g_xxxxyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 204);

    auto g_xxxxyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 205);

    auto g_xxxxyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 206);

    auto g_xxxxyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 207);

    auto g_xxxxyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 208);

    auto g_xxxxyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 209);

    auto g_xxxxyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 210);

    auto g_xxxxyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 211);

    auto g_xxxxyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 212);

    auto g_xxxxyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 213);

    auto g_xxxxyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 214);

    auto g_xxxxyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 215);

    auto g_xxxxyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 216);

    auto g_xxxxyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 217);

    auto g_xxxxyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 218);

    auto g_xxxxyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 219);

    auto g_xxxxyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 220);

    auto g_xxxxyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 221);

    auto g_xxxxyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 222);

    auto g_xxxxyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 223);

    auto g_xxxxyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 224);

    #pragma omp simd aligned(g_xxxxy_0_xxxxxxxy_1, g_xxxxy_0_xxxxxxyy_1, g_xxxxy_0_xxxxxyyy_1, g_xxxxy_0_xxxxyyyy_1, g_xxxxy_0_xxxyyyyy_1, g_xxxxy_0_xxyyyyyy_1, g_xxxxy_0_xyyyyyyy_1, g_xxxxy_0_yyyyyyyy_1, g_xxxxyz_0_xxxxxxxx_0, g_xxxxyz_0_xxxxxxxy_0, g_xxxxyz_0_xxxxxxxz_0, g_xxxxyz_0_xxxxxxyy_0, g_xxxxyz_0_xxxxxxyz_0, g_xxxxyz_0_xxxxxxzz_0, g_xxxxyz_0_xxxxxyyy_0, g_xxxxyz_0_xxxxxyyz_0, g_xxxxyz_0_xxxxxyzz_0, g_xxxxyz_0_xxxxxzzz_0, g_xxxxyz_0_xxxxyyyy_0, g_xxxxyz_0_xxxxyyyz_0, g_xxxxyz_0_xxxxyyzz_0, g_xxxxyz_0_xxxxyzzz_0, g_xxxxyz_0_xxxxzzzz_0, g_xxxxyz_0_xxxyyyyy_0, g_xxxxyz_0_xxxyyyyz_0, g_xxxxyz_0_xxxyyyzz_0, g_xxxxyz_0_xxxyyzzz_0, g_xxxxyz_0_xxxyzzzz_0, g_xxxxyz_0_xxxzzzzz_0, g_xxxxyz_0_xxyyyyyy_0, g_xxxxyz_0_xxyyyyyz_0, g_xxxxyz_0_xxyyyyzz_0, g_xxxxyz_0_xxyyyzzz_0, g_xxxxyz_0_xxyyzzzz_0, g_xxxxyz_0_xxyzzzzz_0, g_xxxxyz_0_xxzzzzzz_0, g_xxxxyz_0_xyyyyyyy_0, g_xxxxyz_0_xyyyyyyz_0, g_xxxxyz_0_xyyyyyzz_0, g_xxxxyz_0_xyyyyzzz_0, g_xxxxyz_0_xyyyzzzz_0, g_xxxxyz_0_xyyzzzzz_0, g_xxxxyz_0_xyzzzzzz_0, g_xxxxyz_0_xzzzzzzz_0, g_xxxxyz_0_yyyyyyyy_0, g_xxxxyz_0_yyyyyyyz_0, g_xxxxyz_0_yyyyyyzz_0, g_xxxxyz_0_yyyyyzzz_0, g_xxxxyz_0_yyyyzzzz_0, g_xxxxyz_0_yyyzzzzz_0, g_xxxxyz_0_yyzzzzzz_0, g_xxxxyz_0_yzzzzzzz_0, g_xxxxyz_0_zzzzzzzz_0, g_xxxxz_0_xxxxxxxx_1, g_xxxxz_0_xxxxxxxz_1, g_xxxxz_0_xxxxxxyz_1, g_xxxxz_0_xxxxxxz_1, g_xxxxz_0_xxxxxxzz_1, g_xxxxz_0_xxxxxyyz_1, g_xxxxz_0_xxxxxyz_1, g_xxxxz_0_xxxxxyzz_1, g_xxxxz_0_xxxxxzz_1, g_xxxxz_0_xxxxxzzz_1, g_xxxxz_0_xxxxyyyz_1, g_xxxxz_0_xxxxyyz_1, g_xxxxz_0_xxxxyyzz_1, g_xxxxz_0_xxxxyzz_1, g_xxxxz_0_xxxxyzzz_1, g_xxxxz_0_xxxxzzz_1, g_xxxxz_0_xxxxzzzz_1, g_xxxxz_0_xxxyyyyz_1, g_xxxxz_0_xxxyyyz_1, g_xxxxz_0_xxxyyyzz_1, g_xxxxz_0_xxxyyzz_1, g_xxxxz_0_xxxyyzzz_1, g_xxxxz_0_xxxyzzz_1, g_xxxxz_0_xxxyzzzz_1, g_xxxxz_0_xxxzzzz_1, g_xxxxz_0_xxxzzzzz_1, g_xxxxz_0_xxyyyyyz_1, g_xxxxz_0_xxyyyyz_1, g_xxxxz_0_xxyyyyzz_1, g_xxxxz_0_xxyyyzz_1, g_xxxxz_0_xxyyyzzz_1, g_xxxxz_0_xxyyzzz_1, g_xxxxz_0_xxyyzzzz_1, g_xxxxz_0_xxyzzzz_1, g_xxxxz_0_xxyzzzzz_1, g_xxxxz_0_xxzzzzz_1, g_xxxxz_0_xxzzzzzz_1, g_xxxxz_0_xyyyyyyz_1, g_xxxxz_0_xyyyyyz_1, g_xxxxz_0_xyyyyyzz_1, g_xxxxz_0_xyyyyzz_1, g_xxxxz_0_xyyyyzzz_1, g_xxxxz_0_xyyyzzz_1, g_xxxxz_0_xyyyzzzz_1, g_xxxxz_0_xyyzzzz_1, g_xxxxz_0_xyyzzzzz_1, g_xxxxz_0_xyzzzzz_1, g_xxxxz_0_xyzzzzzz_1, g_xxxxz_0_xzzzzzz_1, g_xxxxz_0_xzzzzzzz_1, g_xxxxz_0_yyyyyyyz_1, g_xxxxz_0_yyyyyyz_1, g_xxxxz_0_yyyyyyzz_1, g_xxxxz_0_yyyyyzz_1, g_xxxxz_0_yyyyyzzz_1, g_xxxxz_0_yyyyzzz_1, g_xxxxz_0_yyyyzzzz_1, g_xxxxz_0_yyyzzzz_1, g_xxxxz_0_yyyzzzzz_1, g_xxxxz_0_yyzzzzz_1, g_xxxxz_0_yyzzzzzz_1, g_xxxxz_0_yzzzzzz_1, g_xxxxz_0_yzzzzzzz_1, g_xxxxz_0_zzzzzzz_1, g_xxxxz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xxxxxxxx_0[i] = g_xxxxz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxxxy_0[i] = g_xxxxy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxxxxz_0[i] = g_xxxxz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxxyy_0[i] = g_xxxxy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxxxyz_0[i] = g_xxxxz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxxzz_0[i] = g_xxxxz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxyyy_0[i] = g_xxxxy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxxyyz_0[i] = 2.0 * g_xxxxz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxyzz_0[i] = g_xxxxz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxzzz_0[i] = g_xxxxz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxyyyy_0[i] = g_xxxxy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxyyyz_0[i] = 3.0 * g_xxxxz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxyyzz_0[i] = 2.0 * g_xxxxz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxyzzz_0[i] = g_xxxxz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxzzzz_0[i] = g_xxxxz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyyyyy_0[i] = g_xxxxy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxyyyyz_0[i] = 4.0 * g_xxxxz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyyyzz_0[i] = 3.0 * g_xxxxz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyyzzz_0[i] = 2.0 * g_xxxxz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyzzzz_0[i] = g_xxxxz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxzzzzz_0[i] = g_xxxxz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyyyyy_0[i] = g_xxxxy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxyyyyyz_0[i] = 5.0 * g_xxxxz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyyyzz_0[i] = 4.0 * g_xxxxz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyyzzz_0[i] = 3.0 * g_xxxxz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyzzzz_0[i] = 2.0 * g_xxxxz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyzzzzz_0[i] = g_xxxxz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxzzzzzz_0[i] = g_xxxxz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyyyyy_0[i] = g_xxxxy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xyyyyyyz_0[i] = 6.0 * g_xxxxz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyyyzz_0[i] = 5.0 * g_xxxxz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyyzzz_0[i] = 4.0 * g_xxxxz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyzzzz_0[i] = 3.0 * g_xxxxz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyzzzzz_0[i] = 2.0 * g_xxxxz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyzzzzzz_0[i] = g_xxxxz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xzzzzzzz_0[i] = g_xxxxz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyyyyy_0[i] = g_xxxxy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_yyyyyyyz_0[i] = 7.0 * g_xxxxz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyyyzz_0[i] = 6.0 * g_xxxxz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyyzzz_0[i] = 5.0 * g_xxxxz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyzzzz_0[i] = 4.0 * g_xxxxz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyzzzzz_0[i] = 3.0 * g_xxxxz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyzzzzzz_0[i] = 2.0 * g_xxxxz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yzzzzzzz_0[i] = g_xxxxz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_zzzzzzzz_0[i] = g_xxxxz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 225-270 components of targeted buffer : ISL

    auto g_xxxxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 225);

    auto g_xxxxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 226);

    auto g_xxxxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 227);

    auto g_xxxxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 228);

    auto g_xxxxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 229);

    auto g_xxxxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 230);

    auto g_xxxxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 231);

    auto g_xxxxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 232);

    auto g_xxxxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 233);

    auto g_xxxxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 234);

    auto g_xxxxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 235);

    auto g_xxxxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 236);

    auto g_xxxxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 237);

    auto g_xxxxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 238);

    auto g_xxxxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 239);

    auto g_xxxxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 240);

    auto g_xxxxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 241);

    auto g_xxxxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 242);

    auto g_xxxxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 243);

    auto g_xxxxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 244);

    auto g_xxxxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 245);

    auto g_xxxxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 246);

    auto g_xxxxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 247);

    auto g_xxxxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 248);

    auto g_xxxxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 249);

    auto g_xxxxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 250);

    auto g_xxxxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 251);

    auto g_xxxxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 252);

    auto g_xxxxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 253);

    auto g_xxxxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 254);

    auto g_xxxxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 255);

    auto g_xxxxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 256);

    auto g_xxxxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 257);

    auto g_xxxxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 258);

    auto g_xxxxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 259);

    auto g_xxxxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 260);

    auto g_xxxxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 261);

    auto g_xxxxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 262);

    auto g_xxxxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 263);

    auto g_xxxxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 264);

    auto g_xxxxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 265);

    auto g_xxxxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 266);

    auto g_xxxxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 267);

    auto g_xxxxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 268);

    auto g_xxxxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 269);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxxx_0, g_xxxx_0_xxxxxxxx_1, g_xxxx_0_xxxxxxxy_0, g_xxxx_0_xxxxxxxy_1, g_xxxx_0_xxxxxxyy_0, g_xxxx_0_xxxxxxyy_1, g_xxxx_0_xxxxxyyy_0, g_xxxx_0_xxxxxyyy_1, g_xxxx_0_xxxxyyyy_0, g_xxxx_0_xxxxyyyy_1, g_xxxx_0_xxxyyyyy_0, g_xxxx_0_xxxyyyyy_1, g_xxxx_0_xxyyyyyy_0, g_xxxx_0_xxyyyyyy_1, g_xxxx_0_xyyyyyyy_0, g_xxxx_0_xyyyyyyy_1, g_xxxxz_0_xxxxxxxx_1, g_xxxxz_0_xxxxxxxy_1, g_xxxxz_0_xxxxxxyy_1, g_xxxxz_0_xxxxxyyy_1, g_xxxxz_0_xxxxyyyy_1, g_xxxxz_0_xxxyyyyy_1, g_xxxxz_0_xxyyyyyy_1, g_xxxxz_0_xyyyyyyy_1, g_xxxxzz_0_xxxxxxxx_0, g_xxxxzz_0_xxxxxxxy_0, g_xxxxzz_0_xxxxxxxz_0, g_xxxxzz_0_xxxxxxyy_0, g_xxxxzz_0_xxxxxxyz_0, g_xxxxzz_0_xxxxxxzz_0, g_xxxxzz_0_xxxxxyyy_0, g_xxxxzz_0_xxxxxyyz_0, g_xxxxzz_0_xxxxxyzz_0, g_xxxxzz_0_xxxxxzzz_0, g_xxxxzz_0_xxxxyyyy_0, g_xxxxzz_0_xxxxyyyz_0, g_xxxxzz_0_xxxxyyzz_0, g_xxxxzz_0_xxxxyzzz_0, g_xxxxzz_0_xxxxzzzz_0, g_xxxxzz_0_xxxyyyyy_0, g_xxxxzz_0_xxxyyyyz_0, g_xxxxzz_0_xxxyyyzz_0, g_xxxxzz_0_xxxyyzzz_0, g_xxxxzz_0_xxxyzzzz_0, g_xxxxzz_0_xxxzzzzz_0, g_xxxxzz_0_xxyyyyyy_0, g_xxxxzz_0_xxyyyyyz_0, g_xxxxzz_0_xxyyyyzz_0, g_xxxxzz_0_xxyyyzzz_0, g_xxxxzz_0_xxyyzzzz_0, g_xxxxzz_0_xxyzzzzz_0, g_xxxxzz_0_xxzzzzzz_0, g_xxxxzz_0_xyyyyyyy_0, g_xxxxzz_0_xyyyyyyz_0, g_xxxxzz_0_xyyyyyzz_0, g_xxxxzz_0_xyyyyzzz_0, g_xxxxzz_0_xyyyzzzz_0, g_xxxxzz_0_xyyzzzzz_0, g_xxxxzz_0_xyzzzzzz_0, g_xxxxzz_0_xzzzzzzz_0, g_xxxxzz_0_yyyyyyyy_0, g_xxxxzz_0_yyyyyyyz_0, g_xxxxzz_0_yyyyyyzz_0, g_xxxxzz_0_yyyyyzzz_0, g_xxxxzz_0_yyyyzzzz_0, g_xxxxzz_0_yyyzzzzz_0, g_xxxxzz_0_yyzzzzzz_0, g_xxxxzz_0_yzzzzzzz_0, g_xxxxzz_0_zzzzzzzz_0, g_xxxzz_0_xxxxxxxz_1, g_xxxzz_0_xxxxxxyz_1, g_xxxzz_0_xxxxxxz_1, g_xxxzz_0_xxxxxxzz_1, g_xxxzz_0_xxxxxyyz_1, g_xxxzz_0_xxxxxyz_1, g_xxxzz_0_xxxxxyzz_1, g_xxxzz_0_xxxxxzz_1, g_xxxzz_0_xxxxxzzz_1, g_xxxzz_0_xxxxyyyz_1, g_xxxzz_0_xxxxyyz_1, g_xxxzz_0_xxxxyyzz_1, g_xxxzz_0_xxxxyzz_1, g_xxxzz_0_xxxxyzzz_1, g_xxxzz_0_xxxxzzz_1, g_xxxzz_0_xxxxzzzz_1, g_xxxzz_0_xxxyyyyz_1, g_xxxzz_0_xxxyyyz_1, g_xxxzz_0_xxxyyyzz_1, g_xxxzz_0_xxxyyzz_1, g_xxxzz_0_xxxyyzzz_1, g_xxxzz_0_xxxyzzz_1, g_xxxzz_0_xxxyzzzz_1, g_xxxzz_0_xxxzzzz_1, g_xxxzz_0_xxxzzzzz_1, g_xxxzz_0_xxyyyyyz_1, g_xxxzz_0_xxyyyyz_1, g_xxxzz_0_xxyyyyzz_1, g_xxxzz_0_xxyyyzz_1, g_xxxzz_0_xxyyyzzz_1, g_xxxzz_0_xxyyzzz_1, g_xxxzz_0_xxyyzzzz_1, g_xxxzz_0_xxyzzzz_1, g_xxxzz_0_xxyzzzzz_1, g_xxxzz_0_xxzzzzz_1, g_xxxzz_0_xxzzzzzz_1, g_xxxzz_0_xyyyyyyz_1, g_xxxzz_0_xyyyyyz_1, g_xxxzz_0_xyyyyyzz_1, g_xxxzz_0_xyyyyzz_1, g_xxxzz_0_xyyyyzzz_1, g_xxxzz_0_xyyyzzz_1, g_xxxzz_0_xyyyzzzz_1, g_xxxzz_0_xyyzzzz_1, g_xxxzz_0_xyyzzzzz_1, g_xxxzz_0_xyzzzzz_1, g_xxxzz_0_xyzzzzzz_1, g_xxxzz_0_xzzzzzz_1, g_xxxzz_0_xzzzzzzz_1, g_xxxzz_0_yyyyyyyy_1, g_xxxzz_0_yyyyyyyz_1, g_xxxzz_0_yyyyyyz_1, g_xxxzz_0_yyyyyyzz_1, g_xxxzz_0_yyyyyzz_1, g_xxxzz_0_yyyyyzzz_1, g_xxxzz_0_yyyyzzz_1, g_xxxzz_0_yyyyzzzz_1, g_xxxzz_0_yyyzzzz_1, g_xxxzz_0_yyyzzzzz_1, g_xxxzz_0_yyzzzzz_1, g_xxxzz_0_yyzzzzzz_1, g_xxxzz_0_yzzzzzz_1, g_xxxzz_0_yzzzzzzz_1, g_xxxzz_0_zzzzzzz_1, g_xxxzz_0_zzzzzzzz_1, g_xxzz_0_xxxxxxxz_0, g_xxzz_0_xxxxxxxz_1, g_xxzz_0_xxxxxxyz_0, g_xxzz_0_xxxxxxyz_1, g_xxzz_0_xxxxxxzz_0, g_xxzz_0_xxxxxxzz_1, g_xxzz_0_xxxxxyyz_0, g_xxzz_0_xxxxxyyz_1, g_xxzz_0_xxxxxyzz_0, g_xxzz_0_xxxxxyzz_1, g_xxzz_0_xxxxxzzz_0, g_xxzz_0_xxxxxzzz_1, g_xxzz_0_xxxxyyyz_0, g_xxzz_0_xxxxyyyz_1, g_xxzz_0_xxxxyyzz_0, g_xxzz_0_xxxxyyzz_1, g_xxzz_0_xxxxyzzz_0, g_xxzz_0_xxxxyzzz_1, g_xxzz_0_xxxxzzzz_0, g_xxzz_0_xxxxzzzz_1, g_xxzz_0_xxxyyyyz_0, g_xxzz_0_xxxyyyyz_1, g_xxzz_0_xxxyyyzz_0, g_xxzz_0_xxxyyyzz_1, g_xxzz_0_xxxyyzzz_0, g_xxzz_0_xxxyyzzz_1, g_xxzz_0_xxxyzzzz_0, g_xxzz_0_xxxyzzzz_1, g_xxzz_0_xxxzzzzz_0, g_xxzz_0_xxxzzzzz_1, g_xxzz_0_xxyyyyyz_0, g_xxzz_0_xxyyyyyz_1, g_xxzz_0_xxyyyyzz_0, g_xxzz_0_xxyyyyzz_1, g_xxzz_0_xxyyyzzz_0, g_xxzz_0_xxyyyzzz_1, g_xxzz_0_xxyyzzzz_0, g_xxzz_0_xxyyzzzz_1, g_xxzz_0_xxyzzzzz_0, g_xxzz_0_xxyzzzzz_1, g_xxzz_0_xxzzzzzz_0, g_xxzz_0_xxzzzzzz_1, g_xxzz_0_xyyyyyyz_0, g_xxzz_0_xyyyyyyz_1, g_xxzz_0_xyyyyyzz_0, g_xxzz_0_xyyyyyzz_1, g_xxzz_0_xyyyyzzz_0, g_xxzz_0_xyyyyzzz_1, g_xxzz_0_xyyyzzzz_0, g_xxzz_0_xyyyzzzz_1, g_xxzz_0_xyyzzzzz_0, g_xxzz_0_xyyzzzzz_1, g_xxzz_0_xyzzzzzz_0, g_xxzz_0_xyzzzzzz_1, g_xxzz_0_xzzzzzzz_0, g_xxzz_0_xzzzzzzz_1, g_xxzz_0_yyyyyyyy_0, g_xxzz_0_yyyyyyyy_1, g_xxzz_0_yyyyyyyz_0, g_xxzz_0_yyyyyyyz_1, g_xxzz_0_yyyyyyzz_0, g_xxzz_0_yyyyyyzz_1, g_xxzz_0_yyyyyzzz_0, g_xxzz_0_yyyyyzzz_1, g_xxzz_0_yyyyzzzz_0, g_xxzz_0_yyyyzzzz_1, g_xxzz_0_yyyzzzzz_0, g_xxzz_0_yyyzzzzz_1, g_xxzz_0_yyzzzzzz_0, g_xxzz_0_yyzzzzzz_1, g_xxzz_0_yzzzzzzz_0, g_xxzz_0_yzzzzzzz_1, g_xxzz_0_zzzzzzzz_0, g_xxzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xxxxxxxx_0[i] = g_xxxx_0_xxxxxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxxz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxxxy_0[i] = g_xxxx_0_xxxxxxxy_0[i] * fbe_0 - g_xxxx_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxxz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxxxz_0[i] = 3.0 * g_xxzz_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxxzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxxyy_0[i] = g_xxxx_0_xxxxxxyy_0[i] * fbe_0 - g_xxxx_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxxz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxxyz_0[i] = 3.0 * g_xxzz_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxxzz_0[i] = 3.0 * g_xxzz_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxxzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxyyy_0[i] = g_xxxx_0_xxxxxyyy_0[i] * fbe_0 - g_xxxx_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxxz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxyyz_0[i] = 3.0 * g_xxzz_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxyzz_0[i] = 3.0 * g_xxzz_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxzzz_0[i] = 3.0 * g_xxzz_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxxzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxyyyy_0[i] = g_xxxx_0_xxxxyyyy_0[i] * fbe_0 - g_xxxx_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxxz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxyyyz_0[i] = 3.0 * g_xxzz_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxyyzz_0[i] = 3.0 * g_xxzz_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxyzzz_0[i] = 3.0 * g_xxzz_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxzzzz_0[i] = 3.0 * g_xxzz_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyyyyy_0[i] = g_xxxx_0_xxxyyyyy_0[i] * fbe_0 - g_xxxx_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxxz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxyyyyz_0[i] = 3.0 * g_xxzz_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyyyzz_0[i] = 3.0 * g_xxzz_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyyzzz_0[i] = 3.0 * g_xxzz_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyzzzz_0[i] = 3.0 * g_xxzz_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxzzzzz_0[i] = 3.0 * g_xxzz_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyyyyy_0[i] = g_xxxx_0_xxyyyyyy_0[i] * fbe_0 - g_xxxx_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxxz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxyyyyyz_0[i] = 3.0 * g_xxzz_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyyyzz_0[i] = 3.0 * g_xxzz_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyyzzz_0[i] = 3.0 * g_xxzz_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyzzzz_0[i] = 3.0 * g_xxzz_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyzzzzz_0[i] = 3.0 * g_xxzz_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxzzzzzz_0[i] = 3.0 * g_xxzz_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyyyyy_0[i] = g_xxxx_0_xyyyyyyy_0[i] * fbe_0 - g_xxxx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxxz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xyyyyyyz_0[i] = 3.0 * g_xxzz_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyyyzz_0[i] = 3.0 * g_xxzz_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyyzzz_0[i] = 3.0 * g_xxzz_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyzzzz_0[i] = 3.0 * g_xxzz_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyzzzzz_0[i] = 3.0 * g_xxzz_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyzzzzzz_0[i] = 3.0 * g_xxzz_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xzzzzzzz_0[i] = 3.0 * g_xxzz_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyyyy_0[i] = 3.0 * g_xxzz_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyyyz_0[i] = 3.0 * g_xxzz_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyyzz_0[i] = 3.0 * g_xxzz_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyzzz_0[i] = 3.0 * g_xxzz_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyzzzz_0[i] = 3.0 * g_xxzz_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyzzzzz_0[i] = 3.0 * g_xxzz_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyzzzzzz_0[i] = 3.0 * g_xxzz_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yzzzzzzz_0[i] = 3.0 * g_xxzz_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_zzzzzzzz_0[i] = 3.0 * g_xxzz_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 270-315 components of targeted buffer : ISL

    auto g_xxxyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 270);

    auto g_xxxyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 271);

    auto g_xxxyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 272);

    auto g_xxxyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 273);

    auto g_xxxyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 274);

    auto g_xxxyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 275);

    auto g_xxxyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 276);

    auto g_xxxyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 277);

    auto g_xxxyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 278);

    auto g_xxxyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 279);

    auto g_xxxyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 280);

    auto g_xxxyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 281);

    auto g_xxxyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 282);

    auto g_xxxyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 283);

    auto g_xxxyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 284);

    auto g_xxxyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 285);

    auto g_xxxyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 286);

    auto g_xxxyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 287);

    auto g_xxxyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 288);

    auto g_xxxyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 289);

    auto g_xxxyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 290);

    auto g_xxxyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 291);

    auto g_xxxyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 292);

    auto g_xxxyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 293);

    auto g_xxxyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 294);

    auto g_xxxyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 295);

    auto g_xxxyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 296);

    auto g_xxxyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 297);

    auto g_xxxyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 298);

    auto g_xxxyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 299);

    auto g_xxxyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 300);

    auto g_xxxyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 301);

    auto g_xxxyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 302);

    auto g_xxxyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 303);

    auto g_xxxyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 304);

    auto g_xxxyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 305);

    auto g_xxxyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 306);

    auto g_xxxyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 307);

    auto g_xxxyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 308);

    auto g_xxxyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 309);

    auto g_xxxyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 310);

    auto g_xxxyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 311);

    auto g_xxxyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 312);

    auto g_xxxyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 313);

    auto g_xxxyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 314);

    #pragma omp simd aligned(g_xxxy_0_xxxxxxxx_0, g_xxxy_0_xxxxxxxx_1, g_xxxy_0_xxxxxxxz_0, g_xxxy_0_xxxxxxxz_1, g_xxxy_0_xxxxxxzz_0, g_xxxy_0_xxxxxxzz_1, g_xxxy_0_xxxxxzzz_0, g_xxxy_0_xxxxxzzz_1, g_xxxy_0_xxxxzzzz_0, g_xxxy_0_xxxxzzzz_1, g_xxxy_0_xxxzzzzz_0, g_xxxy_0_xxxzzzzz_1, g_xxxy_0_xxzzzzzz_0, g_xxxy_0_xxzzzzzz_1, g_xxxy_0_xzzzzzzz_0, g_xxxy_0_xzzzzzzz_1, g_xxxyy_0_xxxxxxxx_1, g_xxxyy_0_xxxxxxxz_1, g_xxxyy_0_xxxxxxzz_1, g_xxxyy_0_xxxxxzzz_1, g_xxxyy_0_xxxxzzzz_1, g_xxxyy_0_xxxzzzzz_1, g_xxxyy_0_xxzzzzzz_1, g_xxxyy_0_xzzzzzzz_1, g_xxxyyy_0_xxxxxxxx_0, g_xxxyyy_0_xxxxxxxy_0, g_xxxyyy_0_xxxxxxxz_0, g_xxxyyy_0_xxxxxxyy_0, g_xxxyyy_0_xxxxxxyz_0, g_xxxyyy_0_xxxxxxzz_0, g_xxxyyy_0_xxxxxyyy_0, g_xxxyyy_0_xxxxxyyz_0, g_xxxyyy_0_xxxxxyzz_0, g_xxxyyy_0_xxxxxzzz_0, g_xxxyyy_0_xxxxyyyy_0, g_xxxyyy_0_xxxxyyyz_0, g_xxxyyy_0_xxxxyyzz_0, g_xxxyyy_0_xxxxyzzz_0, g_xxxyyy_0_xxxxzzzz_0, g_xxxyyy_0_xxxyyyyy_0, g_xxxyyy_0_xxxyyyyz_0, g_xxxyyy_0_xxxyyyzz_0, g_xxxyyy_0_xxxyyzzz_0, g_xxxyyy_0_xxxyzzzz_0, g_xxxyyy_0_xxxzzzzz_0, g_xxxyyy_0_xxyyyyyy_0, g_xxxyyy_0_xxyyyyyz_0, g_xxxyyy_0_xxyyyyzz_0, g_xxxyyy_0_xxyyyzzz_0, g_xxxyyy_0_xxyyzzzz_0, g_xxxyyy_0_xxyzzzzz_0, g_xxxyyy_0_xxzzzzzz_0, g_xxxyyy_0_xyyyyyyy_0, g_xxxyyy_0_xyyyyyyz_0, g_xxxyyy_0_xyyyyyzz_0, g_xxxyyy_0_xyyyyzzz_0, g_xxxyyy_0_xyyyzzzz_0, g_xxxyyy_0_xyyzzzzz_0, g_xxxyyy_0_xyzzzzzz_0, g_xxxyyy_0_xzzzzzzz_0, g_xxxyyy_0_yyyyyyyy_0, g_xxxyyy_0_yyyyyyyz_0, g_xxxyyy_0_yyyyyyzz_0, g_xxxyyy_0_yyyyyzzz_0, g_xxxyyy_0_yyyyzzzz_0, g_xxxyyy_0_yyyzzzzz_0, g_xxxyyy_0_yyzzzzzz_0, g_xxxyyy_0_yzzzzzzz_0, g_xxxyyy_0_zzzzzzzz_0, g_xxyyy_0_xxxxxxxy_1, g_xxyyy_0_xxxxxxy_1, g_xxyyy_0_xxxxxxyy_1, g_xxyyy_0_xxxxxxyz_1, g_xxyyy_0_xxxxxyy_1, g_xxyyy_0_xxxxxyyy_1, g_xxyyy_0_xxxxxyyz_1, g_xxyyy_0_xxxxxyz_1, g_xxyyy_0_xxxxxyzz_1, g_xxyyy_0_xxxxyyy_1, g_xxyyy_0_xxxxyyyy_1, g_xxyyy_0_xxxxyyyz_1, g_xxyyy_0_xxxxyyz_1, g_xxyyy_0_xxxxyyzz_1, g_xxyyy_0_xxxxyzz_1, g_xxyyy_0_xxxxyzzz_1, g_xxyyy_0_xxxyyyy_1, g_xxyyy_0_xxxyyyyy_1, g_xxyyy_0_xxxyyyyz_1, g_xxyyy_0_xxxyyyz_1, g_xxyyy_0_xxxyyyzz_1, g_xxyyy_0_xxxyyzz_1, g_xxyyy_0_xxxyyzzz_1, g_xxyyy_0_xxxyzzz_1, g_xxyyy_0_xxxyzzzz_1, g_xxyyy_0_xxyyyyy_1, g_xxyyy_0_xxyyyyyy_1, g_xxyyy_0_xxyyyyyz_1, g_xxyyy_0_xxyyyyz_1, g_xxyyy_0_xxyyyyzz_1, g_xxyyy_0_xxyyyzz_1, g_xxyyy_0_xxyyyzzz_1, g_xxyyy_0_xxyyzzz_1, g_xxyyy_0_xxyyzzzz_1, g_xxyyy_0_xxyzzzz_1, g_xxyyy_0_xxyzzzzz_1, g_xxyyy_0_xyyyyyy_1, g_xxyyy_0_xyyyyyyy_1, g_xxyyy_0_xyyyyyyz_1, g_xxyyy_0_xyyyyyz_1, g_xxyyy_0_xyyyyyzz_1, g_xxyyy_0_xyyyyzz_1, g_xxyyy_0_xyyyyzzz_1, g_xxyyy_0_xyyyzzz_1, g_xxyyy_0_xyyyzzzz_1, g_xxyyy_0_xyyzzzz_1, g_xxyyy_0_xyyzzzzz_1, g_xxyyy_0_xyzzzzz_1, g_xxyyy_0_xyzzzzzz_1, g_xxyyy_0_yyyyyyy_1, g_xxyyy_0_yyyyyyyy_1, g_xxyyy_0_yyyyyyyz_1, g_xxyyy_0_yyyyyyz_1, g_xxyyy_0_yyyyyyzz_1, g_xxyyy_0_yyyyyzz_1, g_xxyyy_0_yyyyyzzz_1, g_xxyyy_0_yyyyzzz_1, g_xxyyy_0_yyyyzzzz_1, g_xxyyy_0_yyyzzzz_1, g_xxyyy_0_yyyzzzzz_1, g_xxyyy_0_yyzzzzz_1, g_xxyyy_0_yyzzzzzz_1, g_xxyyy_0_yzzzzzz_1, g_xxyyy_0_yzzzzzzz_1, g_xxyyy_0_zzzzzzzz_1, g_xyyy_0_xxxxxxxy_0, g_xyyy_0_xxxxxxxy_1, g_xyyy_0_xxxxxxyy_0, g_xyyy_0_xxxxxxyy_1, g_xyyy_0_xxxxxxyz_0, g_xyyy_0_xxxxxxyz_1, g_xyyy_0_xxxxxyyy_0, g_xyyy_0_xxxxxyyy_1, g_xyyy_0_xxxxxyyz_0, g_xyyy_0_xxxxxyyz_1, g_xyyy_0_xxxxxyzz_0, g_xyyy_0_xxxxxyzz_1, g_xyyy_0_xxxxyyyy_0, g_xyyy_0_xxxxyyyy_1, g_xyyy_0_xxxxyyyz_0, g_xyyy_0_xxxxyyyz_1, g_xyyy_0_xxxxyyzz_0, g_xyyy_0_xxxxyyzz_1, g_xyyy_0_xxxxyzzz_0, g_xyyy_0_xxxxyzzz_1, g_xyyy_0_xxxyyyyy_0, g_xyyy_0_xxxyyyyy_1, g_xyyy_0_xxxyyyyz_0, g_xyyy_0_xxxyyyyz_1, g_xyyy_0_xxxyyyzz_0, g_xyyy_0_xxxyyyzz_1, g_xyyy_0_xxxyyzzz_0, g_xyyy_0_xxxyyzzz_1, g_xyyy_0_xxxyzzzz_0, g_xyyy_0_xxxyzzzz_1, g_xyyy_0_xxyyyyyy_0, g_xyyy_0_xxyyyyyy_1, g_xyyy_0_xxyyyyyz_0, g_xyyy_0_xxyyyyyz_1, g_xyyy_0_xxyyyyzz_0, g_xyyy_0_xxyyyyzz_1, g_xyyy_0_xxyyyzzz_0, g_xyyy_0_xxyyyzzz_1, g_xyyy_0_xxyyzzzz_0, g_xyyy_0_xxyyzzzz_1, g_xyyy_0_xxyzzzzz_0, g_xyyy_0_xxyzzzzz_1, g_xyyy_0_xyyyyyyy_0, g_xyyy_0_xyyyyyyy_1, g_xyyy_0_xyyyyyyz_0, g_xyyy_0_xyyyyyyz_1, g_xyyy_0_xyyyyyzz_0, g_xyyy_0_xyyyyyzz_1, g_xyyy_0_xyyyyzzz_0, g_xyyy_0_xyyyyzzz_1, g_xyyy_0_xyyyzzzz_0, g_xyyy_0_xyyyzzzz_1, g_xyyy_0_xyyzzzzz_0, g_xyyy_0_xyyzzzzz_1, g_xyyy_0_xyzzzzzz_0, g_xyyy_0_xyzzzzzz_1, g_xyyy_0_yyyyyyyy_0, g_xyyy_0_yyyyyyyy_1, g_xyyy_0_yyyyyyyz_0, g_xyyy_0_yyyyyyyz_1, g_xyyy_0_yyyyyyzz_0, g_xyyy_0_yyyyyyzz_1, g_xyyy_0_yyyyyzzz_0, g_xyyy_0_yyyyyzzz_1, g_xyyy_0_yyyyzzzz_0, g_xyyy_0_yyyyzzzz_1, g_xyyy_0_yyyzzzzz_0, g_xyyy_0_yyyzzzzz_1, g_xyyy_0_yyzzzzzz_0, g_xyyy_0_yyzzzzzz_1, g_xyyy_0_yzzzzzzz_0, g_xyyy_0_yzzzzzzz_1, g_xyyy_0_zzzzzzzz_0, g_xyyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xxxxxxxx_0[i] = 2.0 * g_xxxy_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxxxxy_0[i] = 2.0 * g_xyyy_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxxxz_0[i] = 2.0 * g_xxxy_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxxxyy_0[i] = 2.0 * g_xyyy_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxxyz_0[i] = 2.0 * g_xyyy_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxxzz_0[i] = 2.0 * g_xxxy_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxxyyy_0[i] = 2.0 * g_xyyy_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxyyz_0[i] = 2.0 * g_xyyy_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxyzz_0[i] = 2.0 * g_xyyy_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxzzz_0[i] = 2.0 * g_xxxy_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxyyyy_0[i] = 2.0 * g_xyyy_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxyyyz_0[i] = 2.0 * g_xyyy_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxyyzz_0[i] = 2.0 * g_xyyy_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxyzzz_0[i] = 2.0 * g_xyyy_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxzzzz_0[i] = 2.0 * g_xxxy_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxyyyyy_0[i] = 2.0 * g_xyyy_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyyyyz_0[i] = 2.0 * g_xyyy_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyyyzz_0[i] = 2.0 * g_xyyy_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyyzzz_0[i] = 2.0 * g_xyyy_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyzzzz_0[i] = 2.0 * g_xyyy_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxzzzzz_0[i] = 2.0 * g_xxxy_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxyyyyyy_0[i] = 2.0 * g_xyyy_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyyyyz_0[i] = 2.0 * g_xyyy_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyyyzz_0[i] = 2.0 * g_xyyy_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyyzzz_0[i] = 2.0 * g_xyyy_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyzzzz_0[i] = 2.0 * g_xyyy_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyzzzzz_0[i] = 2.0 * g_xyyy_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxzzzzzz_0[i] = 2.0 * g_xxxy_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xyyyyyyy_0[i] = 2.0 * g_xyyy_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyyyyz_0[i] = 2.0 * g_xyyy_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyyyzz_0[i] = 2.0 * g_xyyy_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyyzzz_0[i] = 2.0 * g_xyyy_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyzzzz_0[i] = 2.0 * g_xyyy_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyzzzzz_0[i] = 2.0 * g_xyyy_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyzzzzzz_0[i] = 2.0 * g_xyyy_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xzzzzzzz_0[i] = 2.0 * g_xxxy_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_yyyyyyyy_0[i] = 2.0 * g_xyyy_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyyyyz_0[i] = 2.0 * g_xyyy_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyyyzz_0[i] = 2.0 * g_xyyy_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyyzzz_0[i] = 2.0 * g_xyyy_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyzzzz_0[i] = 2.0 * g_xyyy_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyzzzzz_0[i] = 2.0 * g_xyyy_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyzzzzzz_0[i] = 2.0 * g_xyyy_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yzzzzzzz_0[i] = 2.0 * g_xyyy_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_zzzzzzzz_0[i] = 2.0 * g_xyyy_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-360 components of targeted buffer : ISL

    auto g_xxxyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 315);

    auto g_xxxyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 316);

    auto g_xxxyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 317);

    auto g_xxxyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 318);

    auto g_xxxyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 319);

    auto g_xxxyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 320);

    auto g_xxxyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 321);

    auto g_xxxyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 322);

    auto g_xxxyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 323);

    auto g_xxxyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 324);

    auto g_xxxyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 325);

    auto g_xxxyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 326);

    auto g_xxxyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 327);

    auto g_xxxyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 328);

    auto g_xxxyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 329);

    auto g_xxxyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 330);

    auto g_xxxyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 331);

    auto g_xxxyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 332);

    auto g_xxxyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 333);

    auto g_xxxyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 334);

    auto g_xxxyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 335);

    auto g_xxxyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 336);

    auto g_xxxyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 337);

    auto g_xxxyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 338);

    auto g_xxxyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 339);

    auto g_xxxyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 340);

    auto g_xxxyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 341);

    auto g_xxxyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 342);

    auto g_xxxyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 343);

    auto g_xxxyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 344);

    auto g_xxxyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 345);

    auto g_xxxyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 346);

    auto g_xxxyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 347);

    auto g_xxxyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 348);

    auto g_xxxyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 349);

    auto g_xxxyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 350);

    auto g_xxxyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 351);

    auto g_xxxyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 352);

    auto g_xxxyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 353);

    auto g_xxxyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 354);

    auto g_xxxyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 355);

    auto g_xxxyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 356);

    auto g_xxxyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 357);

    auto g_xxxyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 358);

    auto g_xxxyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 359);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxxx_1, g_xxxyy_0_xxxxxxxx_1, g_xxxyy_0_xxxxxxxy_1, g_xxxyy_0_xxxxxxxz_1, g_xxxyy_0_xxxxxxy_1, g_xxxyy_0_xxxxxxyy_1, g_xxxyy_0_xxxxxxyz_1, g_xxxyy_0_xxxxxxz_1, g_xxxyy_0_xxxxxxzz_1, g_xxxyy_0_xxxxxyy_1, g_xxxyy_0_xxxxxyyy_1, g_xxxyy_0_xxxxxyyz_1, g_xxxyy_0_xxxxxyz_1, g_xxxyy_0_xxxxxyzz_1, g_xxxyy_0_xxxxxzz_1, g_xxxyy_0_xxxxxzzz_1, g_xxxyy_0_xxxxyyy_1, g_xxxyy_0_xxxxyyyy_1, g_xxxyy_0_xxxxyyyz_1, g_xxxyy_0_xxxxyyz_1, g_xxxyy_0_xxxxyyzz_1, g_xxxyy_0_xxxxyzz_1, g_xxxyy_0_xxxxyzzz_1, g_xxxyy_0_xxxxzzz_1, g_xxxyy_0_xxxxzzzz_1, g_xxxyy_0_xxxyyyy_1, g_xxxyy_0_xxxyyyyy_1, g_xxxyy_0_xxxyyyyz_1, g_xxxyy_0_xxxyyyz_1, g_xxxyy_0_xxxyyyzz_1, g_xxxyy_0_xxxyyzz_1, g_xxxyy_0_xxxyyzzz_1, g_xxxyy_0_xxxyzzz_1, g_xxxyy_0_xxxyzzzz_1, g_xxxyy_0_xxxzzzz_1, g_xxxyy_0_xxxzzzzz_1, g_xxxyy_0_xxyyyyy_1, g_xxxyy_0_xxyyyyyy_1, g_xxxyy_0_xxyyyyyz_1, g_xxxyy_0_xxyyyyz_1, g_xxxyy_0_xxyyyyzz_1, g_xxxyy_0_xxyyyzz_1, g_xxxyy_0_xxyyyzzz_1, g_xxxyy_0_xxyyzzz_1, g_xxxyy_0_xxyyzzzz_1, g_xxxyy_0_xxyzzzz_1, g_xxxyy_0_xxyzzzzz_1, g_xxxyy_0_xxzzzzz_1, g_xxxyy_0_xxzzzzzz_1, g_xxxyy_0_xyyyyyy_1, g_xxxyy_0_xyyyyyyy_1, g_xxxyy_0_xyyyyyyz_1, g_xxxyy_0_xyyyyyz_1, g_xxxyy_0_xyyyyyzz_1, g_xxxyy_0_xyyyyzz_1, g_xxxyy_0_xyyyyzzz_1, g_xxxyy_0_xyyyzzz_1, g_xxxyy_0_xyyyzzzz_1, g_xxxyy_0_xyyzzzz_1, g_xxxyy_0_xyyzzzzz_1, g_xxxyy_0_xyzzzzz_1, g_xxxyy_0_xyzzzzzz_1, g_xxxyy_0_xzzzzzz_1, g_xxxyy_0_xzzzzzzz_1, g_xxxyy_0_yyyyyyy_1, g_xxxyy_0_yyyyyyyy_1, g_xxxyy_0_yyyyyyyz_1, g_xxxyy_0_yyyyyyz_1, g_xxxyy_0_yyyyyyzz_1, g_xxxyy_0_yyyyyzz_1, g_xxxyy_0_yyyyyzzz_1, g_xxxyy_0_yyyyzzz_1, g_xxxyy_0_yyyyzzzz_1, g_xxxyy_0_yyyzzzz_1, g_xxxyy_0_yyyzzzzz_1, g_xxxyy_0_yyzzzzz_1, g_xxxyy_0_yyzzzzzz_1, g_xxxyy_0_yzzzzzz_1, g_xxxyy_0_yzzzzzzz_1, g_xxxyy_0_zzzzzzz_1, g_xxxyy_0_zzzzzzzz_1, g_xxxyyz_0_xxxxxxxx_0, g_xxxyyz_0_xxxxxxxy_0, g_xxxyyz_0_xxxxxxxz_0, g_xxxyyz_0_xxxxxxyy_0, g_xxxyyz_0_xxxxxxyz_0, g_xxxyyz_0_xxxxxxzz_0, g_xxxyyz_0_xxxxxyyy_0, g_xxxyyz_0_xxxxxyyz_0, g_xxxyyz_0_xxxxxyzz_0, g_xxxyyz_0_xxxxxzzz_0, g_xxxyyz_0_xxxxyyyy_0, g_xxxyyz_0_xxxxyyyz_0, g_xxxyyz_0_xxxxyyzz_0, g_xxxyyz_0_xxxxyzzz_0, g_xxxyyz_0_xxxxzzzz_0, g_xxxyyz_0_xxxyyyyy_0, g_xxxyyz_0_xxxyyyyz_0, g_xxxyyz_0_xxxyyyzz_0, g_xxxyyz_0_xxxyyzzz_0, g_xxxyyz_0_xxxyzzzz_0, g_xxxyyz_0_xxxzzzzz_0, g_xxxyyz_0_xxyyyyyy_0, g_xxxyyz_0_xxyyyyyz_0, g_xxxyyz_0_xxyyyyzz_0, g_xxxyyz_0_xxyyyzzz_0, g_xxxyyz_0_xxyyzzzz_0, g_xxxyyz_0_xxyzzzzz_0, g_xxxyyz_0_xxzzzzzz_0, g_xxxyyz_0_xyyyyyyy_0, g_xxxyyz_0_xyyyyyyz_0, g_xxxyyz_0_xyyyyyzz_0, g_xxxyyz_0_xyyyyzzz_0, g_xxxyyz_0_xyyyzzzz_0, g_xxxyyz_0_xyyzzzzz_0, g_xxxyyz_0_xyzzzzzz_0, g_xxxyyz_0_xzzzzzzz_0, g_xxxyyz_0_yyyyyyyy_0, g_xxxyyz_0_yyyyyyyz_0, g_xxxyyz_0_yyyyyyzz_0, g_xxxyyz_0_yyyyyzzz_0, g_xxxyyz_0_yyyyzzzz_0, g_xxxyyz_0_yyyzzzzz_0, g_xxxyyz_0_yyzzzzzz_0, g_xxxyyz_0_yzzzzzzz_0, g_xxxyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xxxxxxxx_0[i] = g_xxxyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxxy_0[i] = g_xxxyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxxz_0[i] = g_xxxyy_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxyy_0[i] = g_xxxyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxyz_0[i] = g_xxxyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxzz_0[i] = 2.0 * g_xxxyy_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxyyy_0[i] = g_xxxyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxyyz_0[i] = g_xxxyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxyzz_0[i] = 2.0 * g_xxxyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxzzz_0[i] = 3.0 * g_xxxyy_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyyyy_0[i] = g_xxxyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyyyz_0[i] = g_xxxyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyyzz_0[i] = 2.0 * g_xxxyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyzzz_0[i] = 3.0 * g_xxxyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxzzzz_0[i] = 4.0 * g_xxxyy_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyyyy_0[i] = g_xxxyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyyyz_0[i] = g_xxxyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyyzz_0[i] = 2.0 * g_xxxyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyzzz_0[i] = 3.0 * g_xxxyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyzzzz_0[i] = 4.0 * g_xxxyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxzzzzz_0[i] = 5.0 * g_xxxyy_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyyyy_0[i] = g_xxxyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyyyz_0[i] = g_xxxyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyyzz_0[i] = 2.0 * g_xxxyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyzzz_0[i] = 3.0 * g_xxxyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyzzzz_0[i] = 4.0 * g_xxxyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyzzzzz_0[i] = 5.0 * g_xxxyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxzzzzzz_0[i] = 6.0 * g_xxxyy_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyyyy_0[i] = g_xxxyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyyyz_0[i] = g_xxxyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyyzz_0[i] = 2.0 * g_xxxyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyzzz_0[i] = 3.0 * g_xxxyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyzzzz_0[i] = 4.0 * g_xxxyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyzzzzz_0[i] = 5.0 * g_xxxyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyzzzzzz_0[i] = 6.0 * g_xxxyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xzzzzzzz_0[i] = 7.0 * g_xxxyy_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyyyy_0[i] = g_xxxyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyyyz_0[i] = g_xxxyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyyzz_0[i] = 2.0 * g_xxxyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyzzz_0[i] = 3.0 * g_xxxyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyzzzz_0[i] = 4.0 * g_xxxyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyzzzzz_0[i] = 5.0 * g_xxxyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyzzzzzz_0[i] = 6.0 * g_xxxyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yzzzzzzz_0[i] = 7.0 * g_xxxyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_zzzzzzzz_0[i] = 8.0 * g_xxxyy_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 360-405 components of targeted buffer : ISL

    auto g_xxxyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 360);

    auto g_xxxyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 361);

    auto g_xxxyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 362);

    auto g_xxxyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 363);

    auto g_xxxyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 364);

    auto g_xxxyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 365);

    auto g_xxxyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 366);

    auto g_xxxyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 367);

    auto g_xxxyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 368);

    auto g_xxxyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 369);

    auto g_xxxyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 370);

    auto g_xxxyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 371);

    auto g_xxxyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 372);

    auto g_xxxyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 373);

    auto g_xxxyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 374);

    auto g_xxxyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 375);

    auto g_xxxyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 376);

    auto g_xxxyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 377);

    auto g_xxxyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 378);

    auto g_xxxyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 379);

    auto g_xxxyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 380);

    auto g_xxxyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 381);

    auto g_xxxyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 382);

    auto g_xxxyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 383);

    auto g_xxxyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 384);

    auto g_xxxyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 385);

    auto g_xxxyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 386);

    auto g_xxxyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 387);

    auto g_xxxyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 388);

    auto g_xxxyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 389);

    auto g_xxxyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 390);

    auto g_xxxyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 391);

    auto g_xxxyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 392);

    auto g_xxxyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 393);

    auto g_xxxyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 394);

    auto g_xxxyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 395);

    auto g_xxxyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 396);

    auto g_xxxyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 397);

    auto g_xxxyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 398);

    auto g_xxxyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 399);

    auto g_xxxyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 400);

    auto g_xxxyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 401);

    auto g_xxxyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 402);

    auto g_xxxyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 403);

    auto g_xxxyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 404);

    #pragma omp simd aligned(g_xxxyzz_0_xxxxxxxx_0, g_xxxyzz_0_xxxxxxxy_0, g_xxxyzz_0_xxxxxxxz_0, g_xxxyzz_0_xxxxxxyy_0, g_xxxyzz_0_xxxxxxyz_0, g_xxxyzz_0_xxxxxxzz_0, g_xxxyzz_0_xxxxxyyy_0, g_xxxyzz_0_xxxxxyyz_0, g_xxxyzz_0_xxxxxyzz_0, g_xxxyzz_0_xxxxxzzz_0, g_xxxyzz_0_xxxxyyyy_0, g_xxxyzz_0_xxxxyyyz_0, g_xxxyzz_0_xxxxyyzz_0, g_xxxyzz_0_xxxxyzzz_0, g_xxxyzz_0_xxxxzzzz_0, g_xxxyzz_0_xxxyyyyy_0, g_xxxyzz_0_xxxyyyyz_0, g_xxxyzz_0_xxxyyyzz_0, g_xxxyzz_0_xxxyyzzz_0, g_xxxyzz_0_xxxyzzzz_0, g_xxxyzz_0_xxxzzzzz_0, g_xxxyzz_0_xxyyyyyy_0, g_xxxyzz_0_xxyyyyyz_0, g_xxxyzz_0_xxyyyyzz_0, g_xxxyzz_0_xxyyyzzz_0, g_xxxyzz_0_xxyyzzzz_0, g_xxxyzz_0_xxyzzzzz_0, g_xxxyzz_0_xxzzzzzz_0, g_xxxyzz_0_xyyyyyyy_0, g_xxxyzz_0_xyyyyyyz_0, g_xxxyzz_0_xyyyyyzz_0, g_xxxyzz_0_xyyyyzzz_0, g_xxxyzz_0_xyyyzzzz_0, g_xxxyzz_0_xyyzzzzz_0, g_xxxyzz_0_xyzzzzzz_0, g_xxxyzz_0_xzzzzzzz_0, g_xxxyzz_0_yyyyyyyy_0, g_xxxyzz_0_yyyyyyyz_0, g_xxxyzz_0_yyyyyyzz_0, g_xxxyzz_0_yyyyyzzz_0, g_xxxyzz_0_yyyyzzzz_0, g_xxxyzz_0_yyyzzzzz_0, g_xxxyzz_0_yyzzzzzz_0, g_xxxyzz_0_yzzzzzzz_0, g_xxxyzz_0_zzzzzzzz_0, g_xxxzz_0_xxxxxxx_1, g_xxxzz_0_xxxxxxxx_1, g_xxxzz_0_xxxxxxxy_1, g_xxxzz_0_xxxxxxxz_1, g_xxxzz_0_xxxxxxy_1, g_xxxzz_0_xxxxxxyy_1, g_xxxzz_0_xxxxxxyz_1, g_xxxzz_0_xxxxxxz_1, g_xxxzz_0_xxxxxxzz_1, g_xxxzz_0_xxxxxyy_1, g_xxxzz_0_xxxxxyyy_1, g_xxxzz_0_xxxxxyyz_1, g_xxxzz_0_xxxxxyz_1, g_xxxzz_0_xxxxxyzz_1, g_xxxzz_0_xxxxxzz_1, g_xxxzz_0_xxxxxzzz_1, g_xxxzz_0_xxxxyyy_1, g_xxxzz_0_xxxxyyyy_1, g_xxxzz_0_xxxxyyyz_1, g_xxxzz_0_xxxxyyz_1, g_xxxzz_0_xxxxyyzz_1, g_xxxzz_0_xxxxyzz_1, g_xxxzz_0_xxxxyzzz_1, g_xxxzz_0_xxxxzzz_1, g_xxxzz_0_xxxxzzzz_1, g_xxxzz_0_xxxyyyy_1, g_xxxzz_0_xxxyyyyy_1, g_xxxzz_0_xxxyyyyz_1, g_xxxzz_0_xxxyyyz_1, g_xxxzz_0_xxxyyyzz_1, g_xxxzz_0_xxxyyzz_1, g_xxxzz_0_xxxyyzzz_1, g_xxxzz_0_xxxyzzz_1, g_xxxzz_0_xxxyzzzz_1, g_xxxzz_0_xxxzzzz_1, g_xxxzz_0_xxxzzzzz_1, g_xxxzz_0_xxyyyyy_1, g_xxxzz_0_xxyyyyyy_1, g_xxxzz_0_xxyyyyyz_1, g_xxxzz_0_xxyyyyz_1, g_xxxzz_0_xxyyyyzz_1, g_xxxzz_0_xxyyyzz_1, g_xxxzz_0_xxyyyzzz_1, g_xxxzz_0_xxyyzzz_1, g_xxxzz_0_xxyyzzzz_1, g_xxxzz_0_xxyzzzz_1, g_xxxzz_0_xxyzzzzz_1, g_xxxzz_0_xxzzzzz_1, g_xxxzz_0_xxzzzzzz_1, g_xxxzz_0_xyyyyyy_1, g_xxxzz_0_xyyyyyyy_1, g_xxxzz_0_xyyyyyyz_1, g_xxxzz_0_xyyyyyz_1, g_xxxzz_0_xyyyyyzz_1, g_xxxzz_0_xyyyyzz_1, g_xxxzz_0_xyyyyzzz_1, g_xxxzz_0_xyyyzzz_1, g_xxxzz_0_xyyyzzzz_1, g_xxxzz_0_xyyzzzz_1, g_xxxzz_0_xyyzzzzz_1, g_xxxzz_0_xyzzzzz_1, g_xxxzz_0_xyzzzzzz_1, g_xxxzz_0_xzzzzzz_1, g_xxxzz_0_xzzzzzzz_1, g_xxxzz_0_yyyyyyy_1, g_xxxzz_0_yyyyyyyy_1, g_xxxzz_0_yyyyyyyz_1, g_xxxzz_0_yyyyyyz_1, g_xxxzz_0_yyyyyyzz_1, g_xxxzz_0_yyyyyzz_1, g_xxxzz_0_yyyyyzzz_1, g_xxxzz_0_yyyyzzz_1, g_xxxzz_0_yyyyzzzz_1, g_xxxzz_0_yyyzzzz_1, g_xxxzz_0_yyyzzzzz_1, g_xxxzz_0_yyzzzzz_1, g_xxxzz_0_yyzzzzzz_1, g_xxxzz_0_yzzzzzz_1, g_xxxzz_0_yzzzzzzz_1, g_xxxzz_0_zzzzzzz_1, g_xxxzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xxxxxxxx_0[i] = g_xxxzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxxy_0[i] = g_xxxzz_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxxz_0[i] = g_xxxzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxyy_0[i] = 2.0 * g_xxxzz_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxyz_0[i] = g_xxxzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxzz_0[i] = g_xxxzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxyyy_0[i] = 3.0 * g_xxxzz_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxyyz_0[i] = 2.0 * g_xxxzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxyzz_0[i] = g_xxxzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxzzz_0[i] = g_xxxzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyyyy_0[i] = 4.0 * g_xxxzz_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyyyz_0[i] = 3.0 * g_xxxzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyyzz_0[i] = 2.0 * g_xxxzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyzzz_0[i] = g_xxxzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxzzzz_0[i] = g_xxxzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyyyy_0[i] = 5.0 * g_xxxzz_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyyyz_0[i] = 4.0 * g_xxxzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyyzz_0[i] = 3.0 * g_xxxzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyzzz_0[i] = 2.0 * g_xxxzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyzzzz_0[i] = g_xxxzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxzzzzz_0[i] = g_xxxzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyyyy_0[i] = 6.0 * g_xxxzz_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyyyz_0[i] = 5.0 * g_xxxzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyyzz_0[i] = 4.0 * g_xxxzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyzzz_0[i] = 3.0 * g_xxxzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyzzzz_0[i] = 2.0 * g_xxxzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyzzzzz_0[i] = g_xxxzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxzzzzzz_0[i] = g_xxxzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyyyy_0[i] = 7.0 * g_xxxzz_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyyyz_0[i] = 6.0 * g_xxxzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyyzz_0[i] = 5.0 * g_xxxzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyzzz_0[i] = 4.0 * g_xxxzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyzzzz_0[i] = 3.0 * g_xxxzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyzzzzz_0[i] = 2.0 * g_xxxzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyzzzzzz_0[i] = g_xxxzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xzzzzzzz_0[i] = g_xxxzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyyyy_0[i] = 8.0 * g_xxxzz_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyyyz_0[i] = 7.0 * g_xxxzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyyzz_0[i] = 6.0 * g_xxxzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyzzz_0[i] = 5.0 * g_xxxzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyzzzz_0[i] = 4.0 * g_xxxzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyzzzzz_0[i] = 3.0 * g_xxxzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyzzzzzz_0[i] = 2.0 * g_xxxzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yzzzzzzz_0[i] = g_xxxzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_zzzzzzzz_0[i] = g_xxxzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 405-450 components of targeted buffer : ISL

    auto g_xxxzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 405);

    auto g_xxxzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 406);

    auto g_xxxzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 407);

    auto g_xxxzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 408);

    auto g_xxxzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 409);

    auto g_xxxzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 410);

    auto g_xxxzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 411);

    auto g_xxxzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 412);

    auto g_xxxzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 413);

    auto g_xxxzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 414);

    auto g_xxxzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 415);

    auto g_xxxzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 416);

    auto g_xxxzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 417);

    auto g_xxxzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 418);

    auto g_xxxzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 419);

    auto g_xxxzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 420);

    auto g_xxxzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 421);

    auto g_xxxzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 422);

    auto g_xxxzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 423);

    auto g_xxxzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 424);

    auto g_xxxzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 425);

    auto g_xxxzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 426);

    auto g_xxxzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 427);

    auto g_xxxzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 428);

    auto g_xxxzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 429);

    auto g_xxxzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 430);

    auto g_xxxzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 431);

    auto g_xxxzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 432);

    auto g_xxxzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 433);

    auto g_xxxzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 434);

    auto g_xxxzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 435);

    auto g_xxxzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 436);

    auto g_xxxzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 437);

    auto g_xxxzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 438);

    auto g_xxxzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 439);

    auto g_xxxzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 440);

    auto g_xxxzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 441);

    auto g_xxxzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 442);

    auto g_xxxzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 443);

    auto g_xxxzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 444);

    auto g_xxxzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 445);

    auto g_xxxzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 446);

    auto g_xxxzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 447);

    auto g_xxxzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 448);

    auto g_xxxzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 449);

    #pragma omp simd aligned(g_xxxz_0_xxxxxxxx_0, g_xxxz_0_xxxxxxxx_1, g_xxxz_0_xxxxxxxy_0, g_xxxz_0_xxxxxxxy_1, g_xxxz_0_xxxxxxyy_0, g_xxxz_0_xxxxxxyy_1, g_xxxz_0_xxxxxyyy_0, g_xxxz_0_xxxxxyyy_1, g_xxxz_0_xxxxyyyy_0, g_xxxz_0_xxxxyyyy_1, g_xxxz_0_xxxyyyyy_0, g_xxxz_0_xxxyyyyy_1, g_xxxz_0_xxyyyyyy_0, g_xxxz_0_xxyyyyyy_1, g_xxxz_0_xyyyyyyy_0, g_xxxz_0_xyyyyyyy_1, g_xxxzz_0_xxxxxxxx_1, g_xxxzz_0_xxxxxxxy_1, g_xxxzz_0_xxxxxxyy_1, g_xxxzz_0_xxxxxyyy_1, g_xxxzz_0_xxxxyyyy_1, g_xxxzz_0_xxxyyyyy_1, g_xxxzz_0_xxyyyyyy_1, g_xxxzz_0_xyyyyyyy_1, g_xxxzzz_0_xxxxxxxx_0, g_xxxzzz_0_xxxxxxxy_0, g_xxxzzz_0_xxxxxxxz_0, g_xxxzzz_0_xxxxxxyy_0, g_xxxzzz_0_xxxxxxyz_0, g_xxxzzz_0_xxxxxxzz_0, g_xxxzzz_0_xxxxxyyy_0, g_xxxzzz_0_xxxxxyyz_0, g_xxxzzz_0_xxxxxyzz_0, g_xxxzzz_0_xxxxxzzz_0, g_xxxzzz_0_xxxxyyyy_0, g_xxxzzz_0_xxxxyyyz_0, g_xxxzzz_0_xxxxyyzz_0, g_xxxzzz_0_xxxxyzzz_0, g_xxxzzz_0_xxxxzzzz_0, g_xxxzzz_0_xxxyyyyy_0, g_xxxzzz_0_xxxyyyyz_0, g_xxxzzz_0_xxxyyyzz_0, g_xxxzzz_0_xxxyyzzz_0, g_xxxzzz_0_xxxyzzzz_0, g_xxxzzz_0_xxxzzzzz_0, g_xxxzzz_0_xxyyyyyy_0, g_xxxzzz_0_xxyyyyyz_0, g_xxxzzz_0_xxyyyyzz_0, g_xxxzzz_0_xxyyyzzz_0, g_xxxzzz_0_xxyyzzzz_0, g_xxxzzz_0_xxyzzzzz_0, g_xxxzzz_0_xxzzzzzz_0, g_xxxzzz_0_xyyyyyyy_0, g_xxxzzz_0_xyyyyyyz_0, g_xxxzzz_0_xyyyyyzz_0, g_xxxzzz_0_xyyyyzzz_0, g_xxxzzz_0_xyyyzzzz_0, g_xxxzzz_0_xyyzzzzz_0, g_xxxzzz_0_xyzzzzzz_0, g_xxxzzz_0_xzzzzzzz_0, g_xxxzzz_0_yyyyyyyy_0, g_xxxzzz_0_yyyyyyyz_0, g_xxxzzz_0_yyyyyyzz_0, g_xxxzzz_0_yyyyyzzz_0, g_xxxzzz_0_yyyyzzzz_0, g_xxxzzz_0_yyyzzzzz_0, g_xxxzzz_0_yyzzzzzz_0, g_xxxzzz_0_yzzzzzzz_0, g_xxxzzz_0_zzzzzzzz_0, g_xxzzz_0_xxxxxxxz_1, g_xxzzz_0_xxxxxxyz_1, g_xxzzz_0_xxxxxxz_1, g_xxzzz_0_xxxxxxzz_1, g_xxzzz_0_xxxxxyyz_1, g_xxzzz_0_xxxxxyz_1, g_xxzzz_0_xxxxxyzz_1, g_xxzzz_0_xxxxxzz_1, g_xxzzz_0_xxxxxzzz_1, g_xxzzz_0_xxxxyyyz_1, g_xxzzz_0_xxxxyyz_1, g_xxzzz_0_xxxxyyzz_1, g_xxzzz_0_xxxxyzz_1, g_xxzzz_0_xxxxyzzz_1, g_xxzzz_0_xxxxzzz_1, g_xxzzz_0_xxxxzzzz_1, g_xxzzz_0_xxxyyyyz_1, g_xxzzz_0_xxxyyyz_1, g_xxzzz_0_xxxyyyzz_1, g_xxzzz_0_xxxyyzz_1, g_xxzzz_0_xxxyyzzz_1, g_xxzzz_0_xxxyzzz_1, g_xxzzz_0_xxxyzzzz_1, g_xxzzz_0_xxxzzzz_1, g_xxzzz_0_xxxzzzzz_1, g_xxzzz_0_xxyyyyyz_1, g_xxzzz_0_xxyyyyz_1, g_xxzzz_0_xxyyyyzz_1, g_xxzzz_0_xxyyyzz_1, g_xxzzz_0_xxyyyzzz_1, g_xxzzz_0_xxyyzzz_1, g_xxzzz_0_xxyyzzzz_1, g_xxzzz_0_xxyzzzz_1, g_xxzzz_0_xxyzzzzz_1, g_xxzzz_0_xxzzzzz_1, g_xxzzz_0_xxzzzzzz_1, g_xxzzz_0_xyyyyyyz_1, g_xxzzz_0_xyyyyyz_1, g_xxzzz_0_xyyyyyzz_1, g_xxzzz_0_xyyyyzz_1, g_xxzzz_0_xyyyyzzz_1, g_xxzzz_0_xyyyzzz_1, g_xxzzz_0_xyyyzzzz_1, g_xxzzz_0_xyyzzzz_1, g_xxzzz_0_xyyzzzzz_1, g_xxzzz_0_xyzzzzz_1, g_xxzzz_0_xyzzzzzz_1, g_xxzzz_0_xzzzzzz_1, g_xxzzz_0_xzzzzzzz_1, g_xxzzz_0_yyyyyyyy_1, g_xxzzz_0_yyyyyyyz_1, g_xxzzz_0_yyyyyyz_1, g_xxzzz_0_yyyyyyzz_1, g_xxzzz_0_yyyyyzz_1, g_xxzzz_0_yyyyyzzz_1, g_xxzzz_0_yyyyzzz_1, g_xxzzz_0_yyyyzzzz_1, g_xxzzz_0_yyyzzzz_1, g_xxzzz_0_yyyzzzzz_1, g_xxzzz_0_yyzzzzz_1, g_xxzzz_0_yyzzzzzz_1, g_xxzzz_0_yzzzzzz_1, g_xxzzz_0_yzzzzzzz_1, g_xxzzz_0_zzzzzzz_1, g_xxzzz_0_zzzzzzzz_1, g_xzzz_0_xxxxxxxz_0, g_xzzz_0_xxxxxxxz_1, g_xzzz_0_xxxxxxyz_0, g_xzzz_0_xxxxxxyz_1, g_xzzz_0_xxxxxxzz_0, g_xzzz_0_xxxxxxzz_1, g_xzzz_0_xxxxxyyz_0, g_xzzz_0_xxxxxyyz_1, g_xzzz_0_xxxxxyzz_0, g_xzzz_0_xxxxxyzz_1, g_xzzz_0_xxxxxzzz_0, g_xzzz_0_xxxxxzzz_1, g_xzzz_0_xxxxyyyz_0, g_xzzz_0_xxxxyyyz_1, g_xzzz_0_xxxxyyzz_0, g_xzzz_0_xxxxyyzz_1, g_xzzz_0_xxxxyzzz_0, g_xzzz_0_xxxxyzzz_1, g_xzzz_0_xxxxzzzz_0, g_xzzz_0_xxxxzzzz_1, g_xzzz_0_xxxyyyyz_0, g_xzzz_0_xxxyyyyz_1, g_xzzz_0_xxxyyyzz_0, g_xzzz_0_xxxyyyzz_1, g_xzzz_0_xxxyyzzz_0, g_xzzz_0_xxxyyzzz_1, g_xzzz_0_xxxyzzzz_0, g_xzzz_0_xxxyzzzz_1, g_xzzz_0_xxxzzzzz_0, g_xzzz_0_xxxzzzzz_1, g_xzzz_0_xxyyyyyz_0, g_xzzz_0_xxyyyyyz_1, g_xzzz_0_xxyyyyzz_0, g_xzzz_0_xxyyyyzz_1, g_xzzz_0_xxyyyzzz_0, g_xzzz_0_xxyyyzzz_1, g_xzzz_0_xxyyzzzz_0, g_xzzz_0_xxyyzzzz_1, g_xzzz_0_xxyzzzzz_0, g_xzzz_0_xxyzzzzz_1, g_xzzz_0_xxzzzzzz_0, g_xzzz_0_xxzzzzzz_1, g_xzzz_0_xyyyyyyz_0, g_xzzz_0_xyyyyyyz_1, g_xzzz_0_xyyyyyzz_0, g_xzzz_0_xyyyyyzz_1, g_xzzz_0_xyyyyzzz_0, g_xzzz_0_xyyyyzzz_1, g_xzzz_0_xyyyzzzz_0, g_xzzz_0_xyyyzzzz_1, g_xzzz_0_xyyzzzzz_0, g_xzzz_0_xyyzzzzz_1, g_xzzz_0_xyzzzzzz_0, g_xzzz_0_xyzzzzzz_1, g_xzzz_0_xzzzzzzz_0, g_xzzz_0_xzzzzzzz_1, g_xzzz_0_yyyyyyyy_0, g_xzzz_0_yyyyyyyy_1, g_xzzz_0_yyyyyyyz_0, g_xzzz_0_yyyyyyyz_1, g_xzzz_0_yyyyyyzz_0, g_xzzz_0_yyyyyyzz_1, g_xzzz_0_yyyyyzzz_0, g_xzzz_0_yyyyyzzz_1, g_xzzz_0_yyyyzzzz_0, g_xzzz_0_yyyyzzzz_1, g_xzzz_0_yyyzzzzz_0, g_xzzz_0_yyyzzzzz_1, g_xzzz_0_yyzzzzzz_0, g_xzzz_0_yyzzzzzz_1, g_xzzz_0_yzzzzzzz_0, g_xzzz_0_yzzzzzzz_1, g_xzzz_0_zzzzzzzz_0, g_xzzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xxxxxxxx_0[i] = 2.0 * g_xxxz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxxxy_0[i] = 2.0 * g_xxxz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxxxz_0[i] = 2.0 * g_xzzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxxz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxxyz_0[i] = 2.0 * g_xzzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxxzz_0[i] = 2.0 * g_xzzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxyyy_0[i] = 2.0 * g_xxxz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxyyz_0[i] = 2.0 * g_xzzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxyzz_0[i] = 2.0 * g_xzzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxzzz_0[i] = 2.0 * g_xzzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxyyyy_0[i] = 2.0 * g_xxxz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxyyyz_0[i] = 2.0 * g_xzzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxyyzz_0[i] = 2.0 * g_xzzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxyzzz_0[i] = 2.0 * g_xzzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxzzzz_0[i] = 2.0 * g_xzzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyyyyy_0[i] = 2.0 * g_xxxz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxyyyyz_0[i] = 2.0 * g_xzzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyyyzz_0[i] = 2.0 * g_xzzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyyzzz_0[i] = 2.0 * g_xzzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyzzzz_0[i] = 2.0 * g_xzzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxzzzzz_0[i] = 2.0 * g_xzzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyyyyy_0[i] = 2.0 * g_xxxz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxyyyyyz_0[i] = 2.0 * g_xzzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyyyzz_0[i] = 2.0 * g_xzzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyyzzz_0[i] = 2.0 * g_xzzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyzzzz_0[i] = 2.0 * g_xzzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyzzzzz_0[i] = 2.0 * g_xzzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxzzzzzz_0[i] = 2.0 * g_xzzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyyyyy_0[i] = 2.0 * g_xxxz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xyyyyyyz_0[i] = 2.0 * g_xzzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyyyzz_0[i] = 2.0 * g_xzzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyyzzz_0[i] = 2.0 * g_xzzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyzzzz_0[i] = 2.0 * g_xzzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyzzzzz_0[i] = 2.0 * g_xzzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyzzzzzz_0[i] = 2.0 * g_xzzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xzzzzzzz_0[i] = 2.0 * g_xzzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyyyy_0[i] = 2.0 * g_xzzz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyyyz_0[i] = 2.0 * g_xzzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyyzz_0[i] = 2.0 * g_xzzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyzzz_0[i] = 2.0 * g_xzzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyzzzz_0[i] = 2.0 * g_xzzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyzzzzz_0[i] = 2.0 * g_xzzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyzzzzzz_0[i] = 2.0 * g_xzzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yzzzzzzz_0[i] = 2.0 * g_xzzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_zzzzzzzz_0[i] = 2.0 * g_xzzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 450-495 components of targeted buffer : ISL

    auto g_xxyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 450);

    auto g_xxyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 451);

    auto g_xxyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 452);

    auto g_xxyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 453);

    auto g_xxyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 454);

    auto g_xxyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 455);

    auto g_xxyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 456);

    auto g_xxyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 457);

    auto g_xxyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 458);

    auto g_xxyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 459);

    auto g_xxyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 460);

    auto g_xxyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 461);

    auto g_xxyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 462);

    auto g_xxyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 463);

    auto g_xxyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 464);

    auto g_xxyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 465);

    auto g_xxyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 466);

    auto g_xxyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 467);

    auto g_xxyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 468);

    auto g_xxyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 469);

    auto g_xxyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 470);

    auto g_xxyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 471);

    auto g_xxyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 472);

    auto g_xxyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 473);

    auto g_xxyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 474);

    auto g_xxyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 475);

    auto g_xxyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 476);

    auto g_xxyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 477);

    auto g_xxyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 478);

    auto g_xxyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 479);

    auto g_xxyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 480);

    auto g_xxyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 481);

    auto g_xxyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 482);

    auto g_xxyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 483);

    auto g_xxyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 484);

    auto g_xxyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 485);

    auto g_xxyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 486);

    auto g_xxyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 487);

    auto g_xxyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 488);

    auto g_xxyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 489);

    auto g_xxyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 490);

    auto g_xxyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 491);

    auto g_xxyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 492);

    auto g_xxyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 493);

    auto g_xxyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 494);

    #pragma omp simd aligned(g_xxyy_0_xxxxxxxx_0, g_xxyy_0_xxxxxxxx_1, g_xxyy_0_xxxxxxxz_0, g_xxyy_0_xxxxxxxz_1, g_xxyy_0_xxxxxxzz_0, g_xxyy_0_xxxxxxzz_1, g_xxyy_0_xxxxxzzz_0, g_xxyy_0_xxxxxzzz_1, g_xxyy_0_xxxxzzzz_0, g_xxyy_0_xxxxzzzz_1, g_xxyy_0_xxxzzzzz_0, g_xxyy_0_xxxzzzzz_1, g_xxyy_0_xxzzzzzz_0, g_xxyy_0_xxzzzzzz_1, g_xxyy_0_xzzzzzzz_0, g_xxyy_0_xzzzzzzz_1, g_xxyyy_0_xxxxxxxx_1, g_xxyyy_0_xxxxxxxz_1, g_xxyyy_0_xxxxxxzz_1, g_xxyyy_0_xxxxxzzz_1, g_xxyyy_0_xxxxzzzz_1, g_xxyyy_0_xxxzzzzz_1, g_xxyyy_0_xxzzzzzz_1, g_xxyyy_0_xzzzzzzz_1, g_xxyyyy_0_xxxxxxxx_0, g_xxyyyy_0_xxxxxxxy_0, g_xxyyyy_0_xxxxxxxz_0, g_xxyyyy_0_xxxxxxyy_0, g_xxyyyy_0_xxxxxxyz_0, g_xxyyyy_0_xxxxxxzz_0, g_xxyyyy_0_xxxxxyyy_0, g_xxyyyy_0_xxxxxyyz_0, g_xxyyyy_0_xxxxxyzz_0, g_xxyyyy_0_xxxxxzzz_0, g_xxyyyy_0_xxxxyyyy_0, g_xxyyyy_0_xxxxyyyz_0, g_xxyyyy_0_xxxxyyzz_0, g_xxyyyy_0_xxxxyzzz_0, g_xxyyyy_0_xxxxzzzz_0, g_xxyyyy_0_xxxyyyyy_0, g_xxyyyy_0_xxxyyyyz_0, g_xxyyyy_0_xxxyyyzz_0, g_xxyyyy_0_xxxyyzzz_0, g_xxyyyy_0_xxxyzzzz_0, g_xxyyyy_0_xxxzzzzz_0, g_xxyyyy_0_xxyyyyyy_0, g_xxyyyy_0_xxyyyyyz_0, g_xxyyyy_0_xxyyyyzz_0, g_xxyyyy_0_xxyyyzzz_0, g_xxyyyy_0_xxyyzzzz_0, g_xxyyyy_0_xxyzzzzz_0, g_xxyyyy_0_xxzzzzzz_0, g_xxyyyy_0_xyyyyyyy_0, g_xxyyyy_0_xyyyyyyz_0, g_xxyyyy_0_xyyyyyzz_0, g_xxyyyy_0_xyyyyzzz_0, g_xxyyyy_0_xyyyzzzz_0, g_xxyyyy_0_xyyzzzzz_0, g_xxyyyy_0_xyzzzzzz_0, g_xxyyyy_0_xzzzzzzz_0, g_xxyyyy_0_yyyyyyyy_0, g_xxyyyy_0_yyyyyyyz_0, g_xxyyyy_0_yyyyyyzz_0, g_xxyyyy_0_yyyyyzzz_0, g_xxyyyy_0_yyyyzzzz_0, g_xxyyyy_0_yyyzzzzz_0, g_xxyyyy_0_yyzzzzzz_0, g_xxyyyy_0_yzzzzzzz_0, g_xxyyyy_0_zzzzzzzz_0, g_xyyyy_0_xxxxxxxy_1, g_xyyyy_0_xxxxxxy_1, g_xyyyy_0_xxxxxxyy_1, g_xyyyy_0_xxxxxxyz_1, g_xyyyy_0_xxxxxyy_1, g_xyyyy_0_xxxxxyyy_1, g_xyyyy_0_xxxxxyyz_1, g_xyyyy_0_xxxxxyz_1, g_xyyyy_0_xxxxxyzz_1, g_xyyyy_0_xxxxyyy_1, g_xyyyy_0_xxxxyyyy_1, g_xyyyy_0_xxxxyyyz_1, g_xyyyy_0_xxxxyyz_1, g_xyyyy_0_xxxxyyzz_1, g_xyyyy_0_xxxxyzz_1, g_xyyyy_0_xxxxyzzz_1, g_xyyyy_0_xxxyyyy_1, g_xyyyy_0_xxxyyyyy_1, g_xyyyy_0_xxxyyyyz_1, g_xyyyy_0_xxxyyyz_1, g_xyyyy_0_xxxyyyzz_1, g_xyyyy_0_xxxyyzz_1, g_xyyyy_0_xxxyyzzz_1, g_xyyyy_0_xxxyzzz_1, g_xyyyy_0_xxxyzzzz_1, g_xyyyy_0_xxyyyyy_1, g_xyyyy_0_xxyyyyyy_1, g_xyyyy_0_xxyyyyyz_1, g_xyyyy_0_xxyyyyz_1, g_xyyyy_0_xxyyyyzz_1, g_xyyyy_0_xxyyyzz_1, g_xyyyy_0_xxyyyzzz_1, g_xyyyy_0_xxyyzzz_1, g_xyyyy_0_xxyyzzzz_1, g_xyyyy_0_xxyzzzz_1, g_xyyyy_0_xxyzzzzz_1, g_xyyyy_0_xyyyyyy_1, g_xyyyy_0_xyyyyyyy_1, g_xyyyy_0_xyyyyyyz_1, g_xyyyy_0_xyyyyyz_1, g_xyyyy_0_xyyyyyzz_1, g_xyyyy_0_xyyyyzz_1, g_xyyyy_0_xyyyyzzz_1, g_xyyyy_0_xyyyzzz_1, g_xyyyy_0_xyyyzzzz_1, g_xyyyy_0_xyyzzzz_1, g_xyyyy_0_xyyzzzzz_1, g_xyyyy_0_xyzzzzz_1, g_xyyyy_0_xyzzzzzz_1, g_xyyyy_0_yyyyyyy_1, g_xyyyy_0_yyyyyyyy_1, g_xyyyy_0_yyyyyyyz_1, g_xyyyy_0_yyyyyyz_1, g_xyyyy_0_yyyyyyzz_1, g_xyyyy_0_yyyyyzz_1, g_xyyyy_0_yyyyyzzz_1, g_xyyyy_0_yyyyzzz_1, g_xyyyy_0_yyyyzzzz_1, g_xyyyy_0_yyyzzzz_1, g_xyyyy_0_yyyzzzzz_1, g_xyyyy_0_yyzzzzz_1, g_xyyyy_0_yyzzzzzz_1, g_xyyyy_0_yzzzzzz_1, g_xyyyy_0_yzzzzzzz_1, g_xyyyy_0_zzzzzzzz_1, g_yyyy_0_xxxxxxxy_0, g_yyyy_0_xxxxxxxy_1, g_yyyy_0_xxxxxxyy_0, g_yyyy_0_xxxxxxyy_1, g_yyyy_0_xxxxxxyz_0, g_yyyy_0_xxxxxxyz_1, g_yyyy_0_xxxxxyyy_0, g_yyyy_0_xxxxxyyy_1, g_yyyy_0_xxxxxyyz_0, g_yyyy_0_xxxxxyyz_1, g_yyyy_0_xxxxxyzz_0, g_yyyy_0_xxxxxyzz_1, g_yyyy_0_xxxxyyyy_0, g_yyyy_0_xxxxyyyy_1, g_yyyy_0_xxxxyyyz_0, g_yyyy_0_xxxxyyyz_1, g_yyyy_0_xxxxyyzz_0, g_yyyy_0_xxxxyyzz_1, g_yyyy_0_xxxxyzzz_0, g_yyyy_0_xxxxyzzz_1, g_yyyy_0_xxxyyyyy_0, g_yyyy_0_xxxyyyyy_1, g_yyyy_0_xxxyyyyz_0, g_yyyy_0_xxxyyyyz_1, g_yyyy_0_xxxyyyzz_0, g_yyyy_0_xxxyyyzz_1, g_yyyy_0_xxxyyzzz_0, g_yyyy_0_xxxyyzzz_1, g_yyyy_0_xxxyzzzz_0, g_yyyy_0_xxxyzzzz_1, g_yyyy_0_xxyyyyyy_0, g_yyyy_0_xxyyyyyy_1, g_yyyy_0_xxyyyyyz_0, g_yyyy_0_xxyyyyyz_1, g_yyyy_0_xxyyyyzz_0, g_yyyy_0_xxyyyyzz_1, g_yyyy_0_xxyyyzzz_0, g_yyyy_0_xxyyyzzz_1, g_yyyy_0_xxyyzzzz_0, g_yyyy_0_xxyyzzzz_1, g_yyyy_0_xxyzzzzz_0, g_yyyy_0_xxyzzzzz_1, g_yyyy_0_xyyyyyyy_0, g_yyyy_0_xyyyyyyy_1, g_yyyy_0_xyyyyyyz_0, g_yyyy_0_xyyyyyyz_1, g_yyyy_0_xyyyyyzz_0, g_yyyy_0_xyyyyyzz_1, g_yyyy_0_xyyyyzzz_0, g_yyyy_0_xyyyyzzz_1, g_yyyy_0_xyyyzzzz_0, g_yyyy_0_xyyyzzzz_1, g_yyyy_0_xyyzzzzz_0, g_yyyy_0_xyyzzzzz_1, g_yyyy_0_xyzzzzzz_0, g_yyyy_0_xyzzzzzz_1, g_yyyy_0_yyyyyyyy_0, g_yyyy_0_yyyyyyyy_1, g_yyyy_0_yyyyyyyz_0, g_yyyy_0_yyyyyyyz_1, g_yyyy_0_yyyyyyzz_0, g_yyyy_0_yyyyyyzz_1, g_yyyy_0_yyyyyzzz_0, g_yyyy_0_yyyyyzzz_1, g_yyyy_0_yyyyzzzz_0, g_yyyy_0_yyyyzzzz_1, g_yyyy_0_yyyzzzzz_0, g_yyyy_0_yyyzzzzz_1, g_yyyy_0_yyzzzzzz_0, g_yyyy_0_yyzzzzzz_1, g_yyyy_0_yzzzzzzz_0, g_yyyy_0_yzzzzzzz_1, g_yyyy_0_zzzzzzzz_0, g_yyyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xxxxxxxx_0[i] = 3.0 * g_xxyy_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxxx_1[i] * fz_be_0 + g_xxyyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxxxxy_0[i] = g_yyyy_0_xxxxxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxxxz_0[i] = 3.0 * g_xxyy_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxxz_1[i] * fz_be_0 + g_xxyyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxxxyy_0[i] = g_yyyy_0_xxxxxxyy_0[i] * fbe_0 - g_yyyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxxyz_0[i] = g_yyyy_0_xxxxxxyz_0[i] * fbe_0 - g_yyyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxxzz_0[i] = 3.0 * g_xxyy_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxzz_1[i] * fz_be_0 + g_xxyyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxxyyy_0[i] = g_yyyy_0_xxxxxyyy_0[i] * fbe_0 - g_yyyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxyyz_0[i] = g_yyyy_0_xxxxxyyz_0[i] * fbe_0 - g_yyyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxyzz_0[i] = g_yyyy_0_xxxxxyzz_0[i] * fbe_0 - g_yyyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxzzz_0[i] = 3.0 * g_xxyy_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxzzz_1[i] * fz_be_0 + g_xxyyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxyyyy_0[i] = g_yyyy_0_xxxxyyyy_0[i] * fbe_0 - g_yyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxyyyz_0[i] = g_yyyy_0_xxxxyyyz_0[i] * fbe_0 - g_yyyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxyyzz_0[i] = g_yyyy_0_xxxxyyzz_0[i] * fbe_0 - g_yyyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxyzzz_0[i] = g_yyyy_0_xxxxyzzz_0[i] * fbe_0 - g_yyyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxzzzz_0[i] = 3.0 * g_xxyy_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxzzzz_1[i] * fz_be_0 + g_xxyyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxyyyyy_0[i] = g_yyyy_0_xxxyyyyy_0[i] * fbe_0 - g_yyyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyyyyz_0[i] = g_yyyy_0_xxxyyyyz_0[i] * fbe_0 - g_yyyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyyyzz_0[i] = g_yyyy_0_xxxyyyzz_0[i] * fbe_0 - g_yyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyyzzz_0[i] = g_yyyy_0_xxxyyzzz_0[i] * fbe_0 - g_yyyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyzzzz_0[i] = g_yyyy_0_xxxyzzzz_0[i] * fbe_0 - g_yyyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxzzzzz_0[i] = 3.0 * g_xxyy_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxzzzzz_1[i] * fz_be_0 + g_xxyyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxyyyyyy_0[i] = g_yyyy_0_xxyyyyyy_0[i] * fbe_0 - g_yyyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyyyyz_0[i] = g_yyyy_0_xxyyyyyz_0[i] * fbe_0 - g_yyyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyyyzz_0[i] = g_yyyy_0_xxyyyyzz_0[i] * fbe_0 - g_yyyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyyzzz_0[i] = g_yyyy_0_xxyyyzzz_0[i] * fbe_0 - g_yyyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyzzzz_0[i] = g_yyyy_0_xxyyzzzz_0[i] * fbe_0 - g_yyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyzzzzz_0[i] = g_yyyy_0_xxyzzzzz_0[i] * fbe_0 - g_yyyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxzzzzzz_0[i] = 3.0 * g_xxyy_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xyyyyyyy_0[i] = g_yyyy_0_xyyyyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyyyyz_0[i] = g_yyyy_0_xyyyyyyz_0[i] * fbe_0 - g_yyyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyyyzz_0[i] = g_yyyy_0_xyyyyyzz_0[i] * fbe_0 - g_yyyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyyzzz_0[i] = g_yyyy_0_xyyyyzzz_0[i] * fbe_0 - g_yyyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyzzzz_0[i] = g_yyyy_0_xyyyzzzz_0[i] * fbe_0 - g_yyyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyzzzzz_0[i] = g_yyyy_0_xyyzzzzz_0[i] * fbe_0 - g_yyyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyzzzzzz_0[i] = g_yyyy_0_xyzzzzzz_0[i] * fbe_0 - g_yyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xzzzzzzz_0[i] = 3.0 * g_xxyy_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xzzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_yyyyyyyy_0[i] = g_yyyy_0_yyyyyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyyyyz_0[i] = g_yyyy_0_yyyyyyyz_0[i] * fbe_0 - g_yyyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyyyzz_0[i] = g_yyyy_0_yyyyyyzz_0[i] * fbe_0 - g_yyyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyyzzz_0[i] = g_yyyy_0_yyyyyzzz_0[i] * fbe_0 - g_yyyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyzzzz_0[i] = g_yyyy_0_yyyyzzzz_0[i] * fbe_0 - g_yyyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyzzzzz_0[i] = g_yyyy_0_yyyzzzzz_0[i] * fbe_0 - g_yyyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyzzzzzz_0[i] = g_yyyy_0_yyzzzzzz_0[i] * fbe_0 - g_yyyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yzzzzzzz_0[i] = g_yyyy_0_yzzzzzzz_0[i] * fbe_0 - g_yyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_zzzzzzzz_0[i] = g_yyyy_0_zzzzzzzz_0[i] * fbe_0 - g_yyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xyyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 495-540 components of targeted buffer : ISL

    auto g_xxyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 495);

    auto g_xxyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 496);

    auto g_xxyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 497);

    auto g_xxyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 498);

    auto g_xxyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 499);

    auto g_xxyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 500);

    auto g_xxyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 501);

    auto g_xxyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 502);

    auto g_xxyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 503);

    auto g_xxyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 504);

    auto g_xxyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 505);

    auto g_xxyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 506);

    auto g_xxyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 507);

    auto g_xxyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 508);

    auto g_xxyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 509);

    auto g_xxyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 510);

    auto g_xxyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 511);

    auto g_xxyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 512);

    auto g_xxyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 513);

    auto g_xxyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 514);

    auto g_xxyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 515);

    auto g_xxyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 516);

    auto g_xxyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 517);

    auto g_xxyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 518);

    auto g_xxyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 519);

    auto g_xxyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 520);

    auto g_xxyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 521);

    auto g_xxyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 522);

    auto g_xxyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 523);

    auto g_xxyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 524);

    auto g_xxyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 525);

    auto g_xxyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 526);

    auto g_xxyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 527);

    auto g_xxyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 528);

    auto g_xxyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 529);

    auto g_xxyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 530);

    auto g_xxyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 531);

    auto g_xxyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 532);

    auto g_xxyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 533);

    auto g_xxyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 534);

    auto g_xxyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 535);

    auto g_xxyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 536);

    auto g_xxyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 537);

    auto g_xxyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 538);

    auto g_xxyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 539);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxxx_1, g_xxyyy_0_xxxxxxxx_1, g_xxyyy_0_xxxxxxxy_1, g_xxyyy_0_xxxxxxxz_1, g_xxyyy_0_xxxxxxy_1, g_xxyyy_0_xxxxxxyy_1, g_xxyyy_0_xxxxxxyz_1, g_xxyyy_0_xxxxxxz_1, g_xxyyy_0_xxxxxxzz_1, g_xxyyy_0_xxxxxyy_1, g_xxyyy_0_xxxxxyyy_1, g_xxyyy_0_xxxxxyyz_1, g_xxyyy_0_xxxxxyz_1, g_xxyyy_0_xxxxxyzz_1, g_xxyyy_0_xxxxxzz_1, g_xxyyy_0_xxxxxzzz_1, g_xxyyy_0_xxxxyyy_1, g_xxyyy_0_xxxxyyyy_1, g_xxyyy_0_xxxxyyyz_1, g_xxyyy_0_xxxxyyz_1, g_xxyyy_0_xxxxyyzz_1, g_xxyyy_0_xxxxyzz_1, g_xxyyy_0_xxxxyzzz_1, g_xxyyy_0_xxxxzzz_1, g_xxyyy_0_xxxxzzzz_1, g_xxyyy_0_xxxyyyy_1, g_xxyyy_0_xxxyyyyy_1, g_xxyyy_0_xxxyyyyz_1, g_xxyyy_0_xxxyyyz_1, g_xxyyy_0_xxxyyyzz_1, g_xxyyy_0_xxxyyzz_1, g_xxyyy_0_xxxyyzzz_1, g_xxyyy_0_xxxyzzz_1, g_xxyyy_0_xxxyzzzz_1, g_xxyyy_0_xxxzzzz_1, g_xxyyy_0_xxxzzzzz_1, g_xxyyy_0_xxyyyyy_1, g_xxyyy_0_xxyyyyyy_1, g_xxyyy_0_xxyyyyyz_1, g_xxyyy_0_xxyyyyz_1, g_xxyyy_0_xxyyyyzz_1, g_xxyyy_0_xxyyyzz_1, g_xxyyy_0_xxyyyzzz_1, g_xxyyy_0_xxyyzzz_1, g_xxyyy_0_xxyyzzzz_1, g_xxyyy_0_xxyzzzz_1, g_xxyyy_0_xxyzzzzz_1, g_xxyyy_0_xxzzzzz_1, g_xxyyy_0_xxzzzzzz_1, g_xxyyy_0_xyyyyyy_1, g_xxyyy_0_xyyyyyyy_1, g_xxyyy_0_xyyyyyyz_1, g_xxyyy_0_xyyyyyz_1, g_xxyyy_0_xyyyyyzz_1, g_xxyyy_0_xyyyyzz_1, g_xxyyy_0_xyyyyzzz_1, g_xxyyy_0_xyyyzzz_1, g_xxyyy_0_xyyyzzzz_1, g_xxyyy_0_xyyzzzz_1, g_xxyyy_0_xyyzzzzz_1, g_xxyyy_0_xyzzzzz_1, g_xxyyy_0_xyzzzzzz_1, g_xxyyy_0_xzzzzzz_1, g_xxyyy_0_xzzzzzzz_1, g_xxyyy_0_yyyyyyy_1, g_xxyyy_0_yyyyyyyy_1, g_xxyyy_0_yyyyyyyz_1, g_xxyyy_0_yyyyyyz_1, g_xxyyy_0_yyyyyyzz_1, g_xxyyy_0_yyyyyzz_1, g_xxyyy_0_yyyyyzzz_1, g_xxyyy_0_yyyyzzz_1, g_xxyyy_0_yyyyzzzz_1, g_xxyyy_0_yyyzzzz_1, g_xxyyy_0_yyyzzzzz_1, g_xxyyy_0_yyzzzzz_1, g_xxyyy_0_yyzzzzzz_1, g_xxyyy_0_yzzzzzz_1, g_xxyyy_0_yzzzzzzz_1, g_xxyyy_0_zzzzzzz_1, g_xxyyy_0_zzzzzzzz_1, g_xxyyyz_0_xxxxxxxx_0, g_xxyyyz_0_xxxxxxxy_0, g_xxyyyz_0_xxxxxxxz_0, g_xxyyyz_0_xxxxxxyy_0, g_xxyyyz_0_xxxxxxyz_0, g_xxyyyz_0_xxxxxxzz_0, g_xxyyyz_0_xxxxxyyy_0, g_xxyyyz_0_xxxxxyyz_0, g_xxyyyz_0_xxxxxyzz_0, g_xxyyyz_0_xxxxxzzz_0, g_xxyyyz_0_xxxxyyyy_0, g_xxyyyz_0_xxxxyyyz_0, g_xxyyyz_0_xxxxyyzz_0, g_xxyyyz_0_xxxxyzzz_0, g_xxyyyz_0_xxxxzzzz_0, g_xxyyyz_0_xxxyyyyy_0, g_xxyyyz_0_xxxyyyyz_0, g_xxyyyz_0_xxxyyyzz_0, g_xxyyyz_0_xxxyyzzz_0, g_xxyyyz_0_xxxyzzzz_0, g_xxyyyz_0_xxxzzzzz_0, g_xxyyyz_0_xxyyyyyy_0, g_xxyyyz_0_xxyyyyyz_0, g_xxyyyz_0_xxyyyyzz_0, g_xxyyyz_0_xxyyyzzz_0, g_xxyyyz_0_xxyyzzzz_0, g_xxyyyz_0_xxyzzzzz_0, g_xxyyyz_0_xxzzzzzz_0, g_xxyyyz_0_xyyyyyyy_0, g_xxyyyz_0_xyyyyyyz_0, g_xxyyyz_0_xyyyyyzz_0, g_xxyyyz_0_xyyyyzzz_0, g_xxyyyz_0_xyyyzzzz_0, g_xxyyyz_0_xyyzzzzz_0, g_xxyyyz_0_xyzzzzzz_0, g_xxyyyz_0_xzzzzzzz_0, g_xxyyyz_0_yyyyyyyy_0, g_xxyyyz_0_yyyyyyyz_0, g_xxyyyz_0_yyyyyyzz_0, g_xxyyyz_0_yyyyyzzz_0, g_xxyyyz_0_yyyyzzzz_0, g_xxyyyz_0_yyyzzzzz_0, g_xxyyyz_0_yyzzzzzz_0, g_xxyyyz_0_yzzzzzzz_0, g_xxyyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xxxxxxxx_0[i] = g_xxyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxxy_0[i] = g_xxyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxxz_0[i] = g_xxyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxyy_0[i] = g_xxyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxyz_0[i] = g_xxyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxzz_0[i] = 2.0 * g_xxyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxyyy_0[i] = g_xxyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxyyz_0[i] = g_xxyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxyzz_0[i] = 2.0 * g_xxyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxzzz_0[i] = 3.0 * g_xxyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyyyy_0[i] = g_xxyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyyyz_0[i] = g_xxyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyyzz_0[i] = 2.0 * g_xxyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyzzz_0[i] = 3.0 * g_xxyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxzzzz_0[i] = 4.0 * g_xxyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyyyy_0[i] = g_xxyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyyyz_0[i] = g_xxyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyyzz_0[i] = 2.0 * g_xxyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyzzz_0[i] = 3.0 * g_xxyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyzzzz_0[i] = 4.0 * g_xxyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxzzzzz_0[i] = 5.0 * g_xxyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyyyy_0[i] = g_xxyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyyyz_0[i] = g_xxyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyyzz_0[i] = 2.0 * g_xxyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyzzz_0[i] = 3.0 * g_xxyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyzzzz_0[i] = 4.0 * g_xxyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyzzzzz_0[i] = 5.0 * g_xxyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxzzzzzz_0[i] = 6.0 * g_xxyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyyyy_0[i] = g_xxyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyyyz_0[i] = g_xxyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyyzz_0[i] = 2.0 * g_xxyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyzzz_0[i] = 3.0 * g_xxyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyzzzz_0[i] = 4.0 * g_xxyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyzzzzz_0[i] = 5.0 * g_xxyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyzzzzzz_0[i] = 6.0 * g_xxyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xzzzzzzz_0[i] = 7.0 * g_xxyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyyyy_0[i] = g_xxyyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyyyz_0[i] = g_xxyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyyzz_0[i] = 2.0 * g_xxyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyzzz_0[i] = 3.0 * g_xxyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyzzzz_0[i] = 4.0 * g_xxyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyzzzzz_0[i] = 5.0 * g_xxyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyzzzzzz_0[i] = 6.0 * g_xxyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yzzzzzzz_0[i] = 7.0 * g_xxyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_zzzzzzzz_0[i] = 8.0 * g_xxyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 540-585 components of targeted buffer : ISL

    auto g_xxyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 540);

    auto g_xxyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 541);

    auto g_xxyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 542);

    auto g_xxyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 543);

    auto g_xxyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 544);

    auto g_xxyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 545);

    auto g_xxyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 546);

    auto g_xxyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 547);

    auto g_xxyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 548);

    auto g_xxyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 549);

    auto g_xxyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 550);

    auto g_xxyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 551);

    auto g_xxyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 552);

    auto g_xxyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 553);

    auto g_xxyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 554);

    auto g_xxyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 555);

    auto g_xxyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 556);

    auto g_xxyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 557);

    auto g_xxyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 558);

    auto g_xxyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 559);

    auto g_xxyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 560);

    auto g_xxyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 561);

    auto g_xxyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 562);

    auto g_xxyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 563);

    auto g_xxyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 564);

    auto g_xxyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 565);

    auto g_xxyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 566);

    auto g_xxyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 567);

    auto g_xxyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 568);

    auto g_xxyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 569);

    auto g_xxyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 570);

    auto g_xxyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 571);

    auto g_xxyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 572);

    auto g_xxyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 573);

    auto g_xxyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 574);

    auto g_xxyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 575);

    auto g_xxyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 576);

    auto g_xxyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 577);

    auto g_xxyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 578);

    auto g_xxyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 579);

    auto g_xxyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 580);

    auto g_xxyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 581);

    auto g_xxyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 582);

    auto g_xxyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 583);

    auto g_xxyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 584);

    #pragma omp simd aligned(g_xxyy_0_xxxxxxxy_0, g_xxyy_0_xxxxxxxy_1, g_xxyy_0_xxxxxxyy_0, g_xxyy_0_xxxxxxyy_1, g_xxyy_0_xxxxxyyy_0, g_xxyy_0_xxxxxyyy_1, g_xxyy_0_xxxxyyyy_0, g_xxyy_0_xxxxyyyy_1, g_xxyy_0_xxxyyyyy_0, g_xxyy_0_xxxyyyyy_1, g_xxyy_0_xxyyyyyy_0, g_xxyy_0_xxyyyyyy_1, g_xxyy_0_xyyyyyyy_0, g_xxyy_0_xyyyyyyy_1, g_xxyyz_0_xxxxxxxy_1, g_xxyyz_0_xxxxxxyy_1, g_xxyyz_0_xxxxxyyy_1, g_xxyyz_0_xxxxyyyy_1, g_xxyyz_0_xxxyyyyy_1, g_xxyyz_0_xxyyyyyy_1, g_xxyyz_0_xyyyyyyy_1, g_xxyyzz_0_xxxxxxxx_0, g_xxyyzz_0_xxxxxxxy_0, g_xxyyzz_0_xxxxxxxz_0, g_xxyyzz_0_xxxxxxyy_0, g_xxyyzz_0_xxxxxxyz_0, g_xxyyzz_0_xxxxxxzz_0, g_xxyyzz_0_xxxxxyyy_0, g_xxyyzz_0_xxxxxyyz_0, g_xxyyzz_0_xxxxxyzz_0, g_xxyyzz_0_xxxxxzzz_0, g_xxyyzz_0_xxxxyyyy_0, g_xxyyzz_0_xxxxyyyz_0, g_xxyyzz_0_xxxxyyzz_0, g_xxyyzz_0_xxxxyzzz_0, g_xxyyzz_0_xxxxzzzz_0, g_xxyyzz_0_xxxyyyyy_0, g_xxyyzz_0_xxxyyyyz_0, g_xxyyzz_0_xxxyyyzz_0, g_xxyyzz_0_xxxyyzzz_0, g_xxyyzz_0_xxxyzzzz_0, g_xxyyzz_0_xxxzzzzz_0, g_xxyyzz_0_xxyyyyyy_0, g_xxyyzz_0_xxyyyyyz_0, g_xxyyzz_0_xxyyyyzz_0, g_xxyyzz_0_xxyyyzzz_0, g_xxyyzz_0_xxyyzzzz_0, g_xxyyzz_0_xxyzzzzz_0, g_xxyyzz_0_xxzzzzzz_0, g_xxyyzz_0_xyyyyyyy_0, g_xxyyzz_0_xyyyyyyz_0, g_xxyyzz_0_xyyyyyzz_0, g_xxyyzz_0_xyyyyzzz_0, g_xxyyzz_0_xyyyzzzz_0, g_xxyyzz_0_xyyzzzzz_0, g_xxyyzz_0_xyzzzzzz_0, g_xxyyzz_0_xzzzzzzz_0, g_xxyyzz_0_yyyyyyyy_0, g_xxyyzz_0_yyyyyyyz_0, g_xxyyzz_0_yyyyyyzz_0, g_xxyyzz_0_yyyyyzzz_0, g_xxyyzz_0_yyyyzzzz_0, g_xxyyzz_0_yyyzzzzz_0, g_xxyyzz_0_yyzzzzzz_0, g_xxyyzz_0_yzzzzzzz_0, g_xxyyzz_0_zzzzzzzz_0, g_xxyzz_0_xxxxxxxx_1, g_xxyzz_0_xxxxxxxz_1, g_xxyzz_0_xxxxxxzz_1, g_xxyzz_0_xxxxxzzz_1, g_xxyzz_0_xxxxzzzz_1, g_xxyzz_0_xxxzzzzz_1, g_xxyzz_0_xxzzzzzz_1, g_xxyzz_0_xzzzzzzz_1, g_xxzz_0_xxxxxxxx_0, g_xxzz_0_xxxxxxxx_1, g_xxzz_0_xxxxxxxz_0, g_xxzz_0_xxxxxxxz_1, g_xxzz_0_xxxxxxzz_0, g_xxzz_0_xxxxxxzz_1, g_xxzz_0_xxxxxzzz_0, g_xxzz_0_xxxxxzzz_1, g_xxzz_0_xxxxzzzz_0, g_xxzz_0_xxxxzzzz_1, g_xxzz_0_xxxzzzzz_0, g_xxzz_0_xxxzzzzz_1, g_xxzz_0_xxzzzzzz_0, g_xxzz_0_xxzzzzzz_1, g_xxzz_0_xzzzzzzz_0, g_xxzz_0_xzzzzzzz_1, g_xyyzz_0_xxxxxxyz_1, g_xyyzz_0_xxxxxyyz_1, g_xyyzz_0_xxxxxyz_1, g_xyyzz_0_xxxxxyzz_1, g_xyyzz_0_xxxxyyyz_1, g_xyyzz_0_xxxxyyz_1, g_xyyzz_0_xxxxyyzz_1, g_xyyzz_0_xxxxyzz_1, g_xyyzz_0_xxxxyzzz_1, g_xyyzz_0_xxxyyyyz_1, g_xyyzz_0_xxxyyyz_1, g_xyyzz_0_xxxyyyzz_1, g_xyyzz_0_xxxyyzz_1, g_xyyzz_0_xxxyyzzz_1, g_xyyzz_0_xxxyzzz_1, g_xyyzz_0_xxxyzzzz_1, g_xyyzz_0_xxyyyyyz_1, g_xyyzz_0_xxyyyyz_1, g_xyyzz_0_xxyyyyzz_1, g_xyyzz_0_xxyyyzz_1, g_xyyzz_0_xxyyyzzz_1, g_xyyzz_0_xxyyzzz_1, g_xyyzz_0_xxyyzzzz_1, g_xyyzz_0_xxyzzzz_1, g_xyyzz_0_xxyzzzzz_1, g_xyyzz_0_xyyyyyyz_1, g_xyyzz_0_xyyyyyz_1, g_xyyzz_0_xyyyyyzz_1, g_xyyzz_0_xyyyyzz_1, g_xyyzz_0_xyyyyzzz_1, g_xyyzz_0_xyyyzzz_1, g_xyyzz_0_xyyyzzzz_1, g_xyyzz_0_xyyzzzz_1, g_xyyzz_0_xyyzzzzz_1, g_xyyzz_0_xyzzzzz_1, g_xyyzz_0_xyzzzzzz_1, g_xyyzz_0_yyyyyyyy_1, g_xyyzz_0_yyyyyyyz_1, g_xyyzz_0_yyyyyyz_1, g_xyyzz_0_yyyyyyzz_1, g_xyyzz_0_yyyyyzz_1, g_xyyzz_0_yyyyyzzz_1, g_xyyzz_0_yyyyzzz_1, g_xyyzz_0_yyyyzzzz_1, g_xyyzz_0_yyyzzzz_1, g_xyyzz_0_yyyzzzzz_1, g_xyyzz_0_yyzzzzz_1, g_xyyzz_0_yyzzzzzz_1, g_xyyzz_0_yzzzzzz_1, g_xyyzz_0_yzzzzzzz_1, g_xyyzz_0_zzzzzzzz_1, g_yyzz_0_xxxxxxyz_0, g_yyzz_0_xxxxxxyz_1, g_yyzz_0_xxxxxyyz_0, g_yyzz_0_xxxxxyyz_1, g_yyzz_0_xxxxxyzz_0, g_yyzz_0_xxxxxyzz_1, g_yyzz_0_xxxxyyyz_0, g_yyzz_0_xxxxyyyz_1, g_yyzz_0_xxxxyyzz_0, g_yyzz_0_xxxxyyzz_1, g_yyzz_0_xxxxyzzz_0, g_yyzz_0_xxxxyzzz_1, g_yyzz_0_xxxyyyyz_0, g_yyzz_0_xxxyyyyz_1, g_yyzz_0_xxxyyyzz_0, g_yyzz_0_xxxyyyzz_1, g_yyzz_0_xxxyyzzz_0, g_yyzz_0_xxxyyzzz_1, g_yyzz_0_xxxyzzzz_0, g_yyzz_0_xxxyzzzz_1, g_yyzz_0_xxyyyyyz_0, g_yyzz_0_xxyyyyyz_1, g_yyzz_0_xxyyyyzz_0, g_yyzz_0_xxyyyyzz_1, g_yyzz_0_xxyyyzzz_0, g_yyzz_0_xxyyyzzz_1, g_yyzz_0_xxyyzzzz_0, g_yyzz_0_xxyyzzzz_1, g_yyzz_0_xxyzzzzz_0, g_yyzz_0_xxyzzzzz_1, g_yyzz_0_xyyyyyyz_0, g_yyzz_0_xyyyyyyz_1, g_yyzz_0_xyyyyyzz_0, g_yyzz_0_xyyyyyzz_1, g_yyzz_0_xyyyyzzz_0, g_yyzz_0_xyyyyzzz_1, g_yyzz_0_xyyyzzzz_0, g_yyzz_0_xyyyzzzz_1, g_yyzz_0_xyyzzzzz_0, g_yyzz_0_xyyzzzzz_1, g_yyzz_0_xyzzzzzz_0, g_yyzz_0_xyzzzzzz_1, g_yyzz_0_yyyyyyyy_0, g_yyzz_0_yyyyyyyy_1, g_yyzz_0_yyyyyyyz_0, g_yyzz_0_yyyyyyyz_1, g_yyzz_0_yyyyyyzz_0, g_yyzz_0_yyyyyyzz_1, g_yyzz_0_yyyyyzzz_0, g_yyzz_0_yyyyyzzz_1, g_yyzz_0_yyyyzzzz_0, g_yyzz_0_yyyyzzzz_1, g_yyzz_0_yyyzzzzz_0, g_yyzz_0_yyyzzzzz_1, g_yyzz_0_yyzzzzzz_0, g_yyzz_0_yyzzzzzz_1, g_yyzz_0_yzzzzzzz_0, g_yyzz_0_yzzzzzzz_1, g_yyzz_0_zzzzzzzz_0, g_yyzz_0_zzzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzz_0_xxxxxxxx_0[i] = g_xxzz_0_xxxxxxxx_0[i] * fbe_0 - g_xxzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxyzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxxxxy_0[i] = g_xxyy_0_xxxxxxxy_0[i] * fbe_0 - g_xxyy_0_xxxxxxxy_1[i] * fz_be_0 + g_xxyyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxxxxz_0[i] = g_xxzz_0_xxxxxxxz_0[i] * fbe_0 - g_xxzz_0_xxxxxxxz_1[i] * fz_be_0 + g_xxyzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxxxyy_0[i] = g_xxyy_0_xxxxxxyy_0[i] * fbe_0 - g_xxyy_0_xxxxxxyy_1[i] * fz_be_0 + g_xxyyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxxxyz_0[i] = g_yyzz_0_xxxxxxyz_0[i] * fbe_0 - g_yyzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxxxzz_0[i] = g_xxzz_0_xxxxxxzz_0[i] * fbe_0 - g_xxzz_0_xxxxxxzz_1[i] * fz_be_0 + g_xxyzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxxyyy_0[i] = g_xxyy_0_xxxxxyyy_0[i] * fbe_0 - g_xxyy_0_xxxxxyyy_1[i] * fz_be_0 + g_xxyyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxxyyz_0[i] = g_yyzz_0_xxxxxyyz_0[i] * fbe_0 - g_yyzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxxyzz_0[i] = g_yyzz_0_xxxxxyzz_0[i] * fbe_0 - g_yyzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxxzzz_0[i] = g_xxzz_0_xxxxxzzz_0[i] * fbe_0 - g_xxzz_0_xxxxxzzz_1[i] * fz_be_0 + g_xxyzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxyyyy_0[i] = g_xxyy_0_xxxxyyyy_0[i] * fbe_0 - g_xxyy_0_xxxxyyyy_1[i] * fz_be_0 + g_xxyyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxyyyz_0[i] = g_yyzz_0_xxxxyyyz_0[i] * fbe_0 - g_yyzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxyyzz_0[i] = g_yyzz_0_xxxxyyzz_0[i] * fbe_0 - g_yyzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxyzzz_0[i] = g_yyzz_0_xxxxyzzz_0[i] * fbe_0 - g_yyzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxzzzz_0[i] = g_xxzz_0_xxxxzzzz_0[i] * fbe_0 - g_xxzz_0_xxxxzzzz_1[i] * fz_be_0 + g_xxyzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxyyyyy_0[i] = g_xxyy_0_xxxyyyyy_0[i] * fbe_0 - g_xxyy_0_xxxyyyyy_1[i] * fz_be_0 + g_xxyyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxyyyyz_0[i] = g_yyzz_0_xxxyyyyz_0[i] * fbe_0 - g_yyzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxyyyzz_0[i] = g_yyzz_0_xxxyyyzz_0[i] * fbe_0 - g_yyzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxyyzzz_0[i] = g_yyzz_0_xxxyyzzz_0[i] * fbe_0 - g_yyzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxyzzzz_0[i] = g_yyzz_0_xxxyzzzz_0[i] * fbe_0 - g_yyzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxzzzzz_0[i] = g_xxzz_0_xxxzzzzz_0[i] * fbe_0 - g_xxzz_0_xxxzzzzz_1[i] * fz_be_0 + g_xxyzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxyyyyyy_0[i] = g_xxyy_0_xxyyyyyy_0[i] * fbe_0 - g_xxyy_0_xxyyyyyy_1[i] * fz_be_0 + g_xxyyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxyyyyyz_0[i] = g_yyzz_0_xxyyyyyz_0[i] * fbe_0 - g_yyzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyyyyzz_0[i] = g_yyzz_0_xxyyyyzz_0[i] * fbe_0 - g_yyzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyyyzzz_0[i] = g_yyzz_0_xxyyyzzz_0[i] * fbe_0 - g_yyzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyyzzzz_0[i] = g_yyzz_0_xxyyzzzz_0[i] * fbe_0 - g_yyzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyzzzzz_0[i] = g_yyzz_0_xxyzzzzz_0[i] * fbe_0 - g_yyzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxzzzzzz_0[i] = g_xxzz_0_xxzzzzzz_0[i] * fbe_0 - g_xxzz_0_xxzzzzzz_1[i] * fz_be_0 + g_xxyzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xyyyyyyy_0[i] = g_xxyy_0_xyyyyyyy_0[i] * fbe_0 - g_xxyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxyyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xyyyyyyz_0[i] = g_yyzz_0_xyyyyyyz_0[i] * fbe_0 - g_yyzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyyyyzz_0[i] = g_yyzz_0_xyyyyyzz_0[i] * fbe_0 - g_yyzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyyyzzz_0[i] = g_yyzz_0_xyyyyzzz_0[i] * fbe_0 - g_yyzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyyzzzz_0[i] = g_yyzz_0_xyyyzzzz_0[i] * fbe_0 - g_yyzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyzzzzz_0[i] = g_yyzz_0_xyyzzzzz_0[i] * fbe_0 - g_yyzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyzzzzzz_0[i] = g_yyzz_0_xyzzzzzz_0[i] * fbe_0 - g_yyzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xzzzzzzz_0[i] = g_xxzz_0_xzzzzzzz_0[i] * fbe_0 - g_xxzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxyzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_yyyyyyyy_0[i] = g_yyzz_0_yyyyyyyy_0[i] * fbe_0 - g_yyzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xyyzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyyyyz_0[i] = g_yyzz_0_yyyyyyyz_0[i] * fbe_0 - g_yyzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyyyzz_0[i] = g_yyzz_0_yyyyyyzz_0[i] * fbe_0 - g_yyzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyyzzz_0[i] = g_yyzz_0_yyyyyzzz_0[i] * fbe_0 - g_yyzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyzzzz_0[i] = g_yyzz_0_yyyyzzzz_0[i] * fbe_0 - g_yyzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyzzzzz_0[i] = g_yyzz_0_yyyzzzzz_0[i] * fbe_0 - g_yyzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyzzzzzz_0[i] = g_yyzz_0_yyzzzzzz_0[i] * fbe_0 - g_yyzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yzzzzzzz_0[i] = g_yyzz_0_yzzzzzzz_0[i] * fbe_0 - g_yyzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_zzzzzzzz_0[i] = g_yyzz_0_zzzzzzzz_0[i] * fbe_0 - g_yyzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xyyzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 585-630 components of targeted buffer : ISL

    auto g_xxyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 585);

    auto g_xxyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 586);

    auto g_xxyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 587);

    auto g_xxyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 588);

    auto g_xxyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 589);

    auto g_xxyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 590);

    auto g_xxyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 591);

    auto g_xxyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 592);

    auto g_xxyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 593);

    auto g_xxyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 594);

    auto g_xxyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 595);

    auto g_xxyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 596);

    auto g_xxyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 597);

    auto g_xxyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 598);

    auto g_xxyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 599);

    auto g_xxyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 600);

    auto g_xxyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 601);

    auto g_xxyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 602);

    auto g_xxyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 603);

    auto g_xxyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 604);

    auto g_xxyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 605);

    auto g_xxyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 606);

    auto g_xxyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 607);

    auto g_xxyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 608);

    auto g_xxyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 609);

    auto g_xxyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 610);

    auto g_xxyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 611);

    auto g_xxyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 612);

    auto g_xxyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 613);

    auto g_xxyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 614);

    auto g_xxyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 615);

    auto g_xxyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 616);

    auto g_xxyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 617);

    auto g_xxyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 618);

    auto g_xxyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 619);

    auto g_xxyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 620);

    auto g_xxyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 621);

    auto g_xxyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 622);

    auto g_xxyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 623);

    auto g_xxyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 624);

    auto g_xxyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 625);

    auto g_xxyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 626);

    auto g_xxyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 627);

    auto g_xxyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 628);

    auto g_xxyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 629);

    #pragma omp simd aligned(g_xxyzzz_0_xxxxxxxx_0, g_xxyzzz_0_xxxxxxxy_0, g_xxyzzz_0_xxxxxxxz_0, g_xxyzzz_0_xxxxxxyy_0, g_xxyzzz_0_xxxxxxyz_0, g_xxyzzz_0_xxxxxxzz_0, g_xxyzzz_0_xxxxxyyy_0, g_xxyzzz_0_xxxxxyyz_0, g_xxyzzz_0_xxxxxyzz_0, g_xxyzzz_0_xxxxxzzz_0, g_xxyzzz_0_xxxxyyyy_0, g_xxyzzz_0_xxxxyyyz_0, g_xxyzzz_0_xxxxyyzz_0, g_xxyzzz_0_xxxxyzzz_0, g_xxyzzz_0_xxxxzzzz_0, g_xxyzzz_0_xxxyyyyy_0, g_xxyzzz_0_xxxyyyyz_0, g_xxyzzz_0_xxxyyyzz_0, g_xxyzzz_0_xxxyyzzz_0, g_xxyzzz_0_xxxyzzzz_0, g_xxyzzz_0_xxxzzzzz_0, g_xxyzzz_0_xxyyyyyy_0, g_xxyzzz_0_xxyyyyyz_0, g_xxyzzz_0_xxyyyyzz_0, g_xxyzzz_0_xxyyyzzz_0, g_xxyzzz_0_xxyyzzzz_0, g_xxyzzz_0_xxyzzzzz_0, g_xxyzzz_0_xxzzzzzz_0, g_xxyzzz_0_xyyyyyyy_0, g_xxyzzz_0_xyyyyyyz_0, g_xxyzzz_0_xyyyyyzz_0, g_xxyzzz_0_xyyyyzzz_0, g_xxyzzz_0_xyyyzzzz_0, g_xxyzzz_0_xyyzzzzz_0, g_xxyzzz_0_xyzzzzzz_0, g_xxyzzz_0_xzzzzzzz_0, g_xxyzzz_0_yyyyyyyy_0, g_xxyzzz_0_yyyyyyyz_0, g_xxyzzz_0_yyyyyyzz_0, g_xxyzzz_0_yyyyyzzz_0, g_xxyzzz_0_yyyyzzzz_0, g_xxyzzz_0_yyyzzzzz_0, g_xxyzzz_0_yyzzzzzz_0, g_xxyzzz_0_yzzzzzzz_0, g_xxyzzz_0_zzzzzzzz_0, g_xxzzz_0_xxxxxxx_1, g_xxzzz_0_xxxxxxxx_1, g_xxzzz_0_xxxxxxxy_1, g_xxzzz_0_xxxxxxxz_1, g_xxzzz_0_xxxxxxy_1, g_xxzzz_0_xxxxxxyy_1, g_xxzzz_0_xxxxxxyz_1, g_xxzzz_0_xxxxxxz_1, g_xxzzz_0_xxxxxxzz_1, g_xxzzz_0_xxxxxyy_1, g_xxzzz_0_xxxxxyyy_1, g_xxzzz_0_xxxxxyyz_1, g_xxzzz_0_xxxxxyz_1, g_xxzzz_0_xxxxxyzz_1, g_xxzzz_0_xxxxxzz_1, g_xxzzz_0_xxxxxzzz_1, g_xxzzz_0_xxxxyyy_1, g_xxzzz_0_xxxxyyyy_1, g_xxzzz_0_xxxxyyyz_1, g_xxzzz_0_xxxxyyz_1, g_xxzzz_0_xxxxyyzz_1, g_xxzzz_0_xxxxyzz_1, g_xxzzz_0_xxxxyzzz_1, g_xxzzz_0_xxxxzzz_1, g_xxzzz_0_xxxxzzzz_1, g_xxzzz_0_xxxyyyy_1, g_xxzzz_0_xxxyyyyy_1, g_xxzzz_0_xxxyyyyz_1, g_xxzzz_0_xxxyyyz_1, g_xxzzz_0_xxxyyyzz_1, g_xxzzz_0_xxxyyzz_1, g_xxzzz_0_xxxyyzzz_1, g_xxzzz_0_xxxyzzz_1, g_xxzzz_0_xxxyzzzz_1, g_xxzzz_0_xxxzzzz_1, g_xxzzz_0_xxxzzzzz_1, g_xxzzz_0_xxyyyyy_1, g_xxzzz_0_xxyyyyyy_1, g_xxzzz_0_xxyyyyyz_1, g_xxzzz_0_xxyyyyz_1, g_xxzzz_0_xxyyyyzz_1, g_xxzzz_0_xxyyyzz_1, g_xxzzz_0_xxyyyzzz_1, g_xxzzz_0_xxyyzzz_1, g_xxzzz_0_xxyyzzzz_1, g_xxzzz_0_xxyzzzz_1, g_xxzzz_0_xxyzzzzz_1, g_xxzzz_0_xxzzzzz_1, g_xxzzz_0_xxzzzzzz_1, g_xxzzz_0_xyyyyyy_1, g_xxzzz_0_xyyyyyyy_1, g_xxzzz_0_xyyyyyyz_1, g_xxzzz_0_xyyyyyz_1, g_xxzzz_0_xyyyyyzz_1, g_xxzzz_0_xyyyyzz_1, g_xxzzz_0_xyyyyzzz_1, g_xxzzz_0_xyyyzzz_1, g_xxzzz_0_xyyyzzzz_1, g_xxzzz_0_xyyzzzz_1, g_xxzzz_0_xyyzzzzz_1, g_xxzzz_0_xyzzzzz_1, g_xxzzz_0_xyzzzzzz_1, g_xxzzz_0_xzzzzzz_1, g_xxzzz_0_xzzzzzzz_1, g_xxzzz_0_yyyyyyy_1, g_xxzzz_0_yyyyyyyy_1, g_xxzzz_0_yyyyyyyz_1, g_xxzzz_0_yyyyyyz_1, g_xxzzz_0_yyyyyyzz_1, g_xxzzz_0_yyyyyzz_1, g_xxzzz_0_yyyyyzzz_1, g_xxzzz_0_yyyyzzz_1, g_xxzzz_0_yyyyzzzz_1, g_xxzzz_0_yyyzzzz_1, g_xxzzz_0_yyyzzzzz_1, g_xxzzz_0_yyzzzzz_1, g_xxzzz_0_yyzzzzzz_1, g_xxzzz_0_yzzzzzz_1, g_xxzzz_0_yzzzzzzz_1, g_xxzzz_0_zzzzzzz_1, g_xxzzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xxxxxxxx_0[i] = g_xxzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxxy_0[i] = g_xxzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxxz_0[i] = g_xxzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxyz_0[i] = g_xxzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxzz_0[i] = g_xxzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxyyy_0[i] = 3.0 * g_xxzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxyyz_0[i] = 2.0 * g_xxzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxyzz_0[i] = g_xxzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxzzz_0[i] = g_xxzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyyyy_0[i] = 4.0 * g_xxzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyyyz_0[i] = 3.0 * g_xxzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyyzz_0[i] = 2.0 * g_xxzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyzzz_0[i] = g_xxzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxzzzz_0[i] = g_xxzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyyyy_0[i] = 5.0 * g_xxzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyyyz_0[i] = 4.0 * g_xxzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyyzz_0[i] = 3.0 * g_xxzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyzzz_0[i] = 2.0 * g_xxzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyzzzz_0[i] = g_xxzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxzzzzz_0[i] = g_xxzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyyyy_0[i] = 6.0 * g_xxzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyyyz_0[i] = 5.0 * g_xxzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyyzz_0[i] = 4.0 * g_xxzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyzzz_0[i] = 3.0 * g_xxzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyzzzz_0[i] = 2.0 * g_xxzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyzzzzz_0[i] = g_xxzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxzzzzzz_0[i] = g_xxzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyyyy_0[i] = 7.0 * g_xxzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyyyz_0[i] = 6.0 * g_xxzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyyzz_0[i] = 5.0 * g_xxzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyzzz_0[i] = 4.0 * g_xxzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyzzzz_0[i] = 3.0 * g_xxzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyzzzzz_0[i] = 2.0 * g_xxzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyzzzzzz_0[i] = g_xxzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xzzzzzzz_0[i] = g_xxzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyyyy_0[i] = 8.0 * g_xxzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyyyz_0[i] = 7.0 * g_xxzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyyzz_0[i] = 6.0 * g_xxzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyzzz_0[i] = 5.0 * g_xxzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyzzzz_0[i] = 4.0 * g_xxzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyzzzzz_0[i] = 3.0 * g_xxzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyzzzzzz_0[i] = 2.0 * g_xxzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yzzzzzzz_0[i] = g_xxzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_zzzzzzzz_0[i] = g_xxzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 630-675 components of targeted buffer : ISL

    auto g_xxzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 630);

    auto g_xxzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 631);

    auto g_xxzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 632);

    auto g_xxzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 633);

    auto g_xxzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 634);

    auto g_xxzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 635);

    auto g_xxzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 636);

    auto g_xxzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 637);

    auto g_xxzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 638);

    auto g_xxzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 639);

    auto g_xxzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 640);

    auto g_xxzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 641);

    auto g_xxzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 642);

    auto g_xxzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 643);

    auto g_xxzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 644);

    auto g_xxzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 645);

    auto g_xxzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 646);

    auto g_xxzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 647);

    auto g_xxzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 648);

    auto g_xxzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 649);

    auto g_xxzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 650);

    auto g_xxzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 651);

    auto g_xxzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 652);

    auto g_xxzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 653);

    auto g_xxzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 654);

    auto g_xxzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 655);

    auto g_xxzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 656);

    auto g_xxzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 657);

    auto g_xxzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 658);

    auto g_xxzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 659);

    auto g_xxzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 660);

    auto g_xxzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 661);

    auto g_xxzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 662);

    auto g_xxzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 663);

    auto g_xxzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 664);

    auto g_xxzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 665);

    auto g_xxzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 666);

    auto g_xxzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 667);

    auto g_xxzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 668);

    auto g_xxzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 669);

    auto g_xxzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 670);

    auto g_xxzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 671);

    auto g_xxzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 672);

    auto g_xxzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 673);

    auto g_xxzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 674);

    #pragma omp simd aligned(g_xxzz_0_xxxxxxxx_0, g_xxzz_0_xxxxxxxx_1, g_xxzz_0_xxxxxxxy_0, g_xxzz_0_xxxxxxxy_1, g_xxzz_0_xxxxxxyy_0, g_xxzz_0_xxxxxxyy_1, g_xxzz_0_xxxxxyyy_0, g_xxzz_0_xxxxxyyy_1, g_xxzz_0_xxxxyyyy_0, g_xxzz_0_xxxxyyyy_1, g_xxzz_0_xxxyyyyy_0, g_xxzz_0_xxxyyyyy_1, g_xxzz_0_xxyyyyyy_0, g_xxzz_0_xxyyyyyy_1, g_xxzz_0_xyyyyyyy_0, g_xxzz_0_xyyyyyyy_1, g_xxzzz_0_xxxxxxxx_1, g_xxzzz_0_xxxxxxxy_1, g_xxzzz_0_xxxxxxyy_1, g_xxzzz_0_xxxxxyyy_1, g_xxzzz_0_xxxxyyyy_1, g_xxzzz_0_xxxyyyyy_1, g_xxzzz_0_xxyyyyyy_1, g_xxzzz_0_xyyyyyyy_1, g_xxzzzz_0_xxxxxxxx_0, g_xxzzzz_0_xxxxxxxy_0, g_xxzzzz_0_xxxxxxxz_0, g_xxzzzz_0_xxxxxxyy_0, g_xxzzzz_0_xxxxxxyz_0, g_xxzzzz_0_xxxxxxzz_0, g_xxzzzz_0_xxxxxyyy_0, g_xxzzzz_0_xxxxxyyz_0, g_xxzzzz_0_xxxxxyzz_0, g_xxzzzz_0_xxxxxzzz_0, g_xxzzzz_0_xxxxyyyy_0, g_xxzzzz_0_xxxxyyyz_0, g_xxzzzz_0_xxxxyyzz_0, g_xxzzzz_0_xxxxyzzz_0, g_xxzzzz_0_xxxxzzzz_0, g_xxzzzz_0_xxxyyyyy_0, g_xxzzzz_0_xxxyyyyz_0, g_xxzzzz_0_xxxyyyzz_0, g_xxzzzz_0_xxxyyzzz_0, g_xxzzzz_0_xxxyzzzz_0, g_xxzzzz_0_xxxzzzzz_0, g_xxzzzz_0_xxyyyyyy_0, g_xxzzzz_0_xxyyyyyz_0, g_xxzzzz_0_xxyyyyzz_0, g_xxzzzz_0_xxyyyzzz_0, g_xxzzzz_0_xxyyzzzz_0, g_xxzzzz_0_xxyzzzzz_0, g_xxzzzz_0_xxzzzzzz_0, g_xxzzzz_0_xyyyyyyy_0, g_xxzzzz_0_xyyyyyyz_0, g_xxzzzz_0_xyyyyyzz_0, g_xxzzzz_0_xyyyyzzz_0, g_xxzzzz_0_xyyyzzzz_0, g_xxzzzz_0_xyyzzzzz_0, g_xxzzzz_0_xyzzzzzz_0, g_xxzzzz_0_xzzzzzzz_0, g_xxzzzz_0_yyyyyyyy_0, g_xxzzzz_0_yyyyyyyz_0, g_xxzzzz_0_yyyyyyzz_0, g_xxzzzz_0_yyyyyzzz_0, g_xxzzzz_0_yyyyzzzz_0, g_xxzzzz_0_yyyzzzzz_0, g_xxzzzz_0_yyzzzzzz_0, g_xxzzzz_0_yzzzzzzz_0, g_xxzzzz_0_zzzzzzzz_0, g_xzzzz_0_xxxxxxxz_1, g_xzzzz_0_xxxxxxyz_1, g_xzzzz_0_xxxxxxz_1, g_xzzzz_0_xxxxxxzz_1, g_xzzzz_0_xxxxxyyz_1, g_xzzzz_0_xxxxxyz_1, g_xzzzz_0_xxxxxyzz_1, g_xzzzz_0_xxxxxzz_1, g_xzzzz_0_xxxxxzzz_1, g_xzzzz_0_xxxxyyyz_1, g_xzzzz_0_xxxxyyz_1, g_xzzzz_0_xxxxyyzz_1, g_xzzzz_0_xxxxyzz_1, g_xzzzz_0_xxxxyzzz_1, g_xzzzz_0_xxxxzzz_1, g_xzzzz_0_xxxxzzzz_1, g_xzzzz_0_xxxyyyyz_1, g_xzzzz_0_xxxyyyz_1, g_xzzzz_0_xxxyyyzz_1, g_xzzzz_0_xxxyyzz_1, g_xzzzz_0_xxxyyzzz_1, g_xzzzz_0_xxxyzzz_1, g_xzzzz_0_xxxyzzzz_1, g_xzzzz_0_xxxzzzz_1, g_xzzzz_0_xxxzzzzz_1, g_xzzzz_0_xxyyyyyz_1, g_xzzzz_0_xxyyyyz_1, g_xzzzz_0_xxyyyyzz_1, g_xzzzz_0_xxyyyzz_1, g_xzzzz_0_xxyyyzzz_1, g_xzzzz_0_xxyyzzz_1, g_xzzzz_0_xxyyzzzz_1, g_xzzzz_0_xxyzzzz_1, g_xzzzz_0_xxyzzzzz_1, g_xzzzz_0_xxzzzzz_1, g_xzzzz_0_xxzzzzzz_1, g_xzzzz_0_xyyyyyyz_1, g_xzzzz_0_xyyyyyz_1, g_xzzzz_0_xyyyyyzz_1, g_xzzzz_0_xyyyyzz_1, g_xzzzz_0_xyyyyzzz_1, g_xzzzz_0_xyyyzzz_1, g_xzzzz_0_xyyyzzzz_1, g_xzzzz_0_xyyzzzz_1, g_xzzzz_0_xyyzzzzz_1, g_xzzzz_0_xyzzzzz_1, g_xzzzz_0_xyzzzzzz_1, g_xzzzz_0_xzzzzzz_1, g_xzzzz_0_xzzzzzzz_1, g_xzzzz_0_yyyyyyyy_1, g_xzzzz_0_yyyyyyyz_1, g_xzzzz_0_yyyyyyz_1, g_xzzzz_0_yyyyyyzz_1, g_xzzzz_0_yyyyyzz_1, g_xzzzz_0_yyyyyzzz_1, g_xzzzz_0_yyyyzzz_1, g_xzzzz_0_yyyyzzzz_1, g_xzzzz_0_yyyzzzz_1, g_xzzzz_0_yyyzzzzz_1, g_xzzzz_0_yyzzzzz_1, g_xzzzz_0_yyzzzzzz_1, g_xzzzz_0_yzzzzzz_1, g_xzzzz_0_yzzzzzzz_1, g_xzzzz_0_zzzzzzz_1, g_xzzzz_0_zzzzzzzz_1, g_zzzz_0_xxxxxxxz_0, g_zzzz_0_xxxxxxxz_1, g_zzzz_0_xxxxxxyz_0, g_zzzz_0_xxxxxxyz_1, g_zzzz_0_xxxxxxzz_0, g_zzzz_0_xxxxxxzz_1, g_zzzz_0_xxxxxyyz_0, g_zzzz_0_xxxxxyyz_1, g_zzzz_0_xxxxxyzz_0, g_zzzz_0_xxxxxyzz_1, g_zzzz_0_xxxxxzzz_0, g_zzzz_0_xxxxxzzz_1, g_zzzz_0_xxxxyyyz_0, g_zzzz_0_xxxxyyyz_1, g_zzzz_0_xxxxyyzz_0, g_zzzz_0_xxxxyyzz_1, g_zzzz_0_xxxxyzzz_0, g_zzzz_0_xxxxyzzz_1, g_zzzz_0_xxxxzzzz_0, g_zzzz_0_xxxxzzzz_1, g_zzzz_0_xxxyyyyz_0, g_zzzz_0_xxxyyyyz_1, g_zzzz_0_xxxyyyzz_0, g_zzzz_0_xxxyyyzz_1, g_zzzz_0_xxxyyzzz_0, g_zzzz_0_xxxyyzzz_1, g_zzzz_0_xxxyzzzz_0, g_zzzz_0_xxxyzzzz_1, g_zzzz_0_xxxzzzzz_0, g_zzzz_0_xxxzzzzz_1, g_zzzz_0_xxyyyyyz_0, g_zzzz_0_xxyyyyyz_1, g_zzzz_0_xxyyyyzz_0, g_zzzz_0_xxyyyyzz_1, g_zzzz_0_xxyyyzzz_0, g_zzzz_0_xxyyyzzz_1, g_zzzz_0_xxyyzzzz_0, g_zzzz_0_xxyyzzzz_1, g_zzzz_0_xxyzzzzz_0, g_zzzz_0_xxyzzzzz_1, g_zzzz_0_xxzzzzzz_0, g_zzzz_0_xxzzzzzz_1, g_zzzz_0_xyyyyyyz_0, g_zzzz_0_xyyyyyyz_1, g_zzzz_0_xyyyyyzz_0, g_zzzz_0_xyyyyyzz_1, g_zzzz_0_xyyyyzzz_0, g_zzzz_0_xyyyyzzz_1, g_zzzz_0_xyyyzzzz_0, g_zzzz_0_xyyyzzzz_1, g_zzzz_0_xyyzzzzz_0, g_zzzz_0_xyyzzzzz_1, g_zzzz_0_xyzzzzzz_0, g_zzzz_0_xyzzzzzz_1, g_zzzz_0_xzzzzzzz_0, g_zzzz_0_xzzzzzzz_1, g_zzzz_0_yyyyyyyy_0, g_zzzz_0_yyyyyyyy_1, g_zzzz_0_yyyyyyyz_0, g_zzzz_0_yyyyyyyz_1, g_zzzz_0_yyyyyyzz_0, g_zzzz_0_yyyyyyzz_1, g_zzzz_0_yyyyyzzz_0, g_zzzz_0_yyyyyzzz_1, g_zzzz_0_yyyyzzzz_0, g_zzzz_0_yyyyzzzz_1, g_zzzz_0_yyyzzzzz_0, g_zzzz_0_yyyzzzzz_1, g_zzzz_0_yyzzzzzz_0, g_zzzz_0_yyzzzzzz_1, g_zzzz_0_yzzzzzzz_0, g_zzzz_0_yzzzzzzz_1, g_zzzz_0_zzzzzzzz_0, g_zzzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xxxxxxxx_0[i] = 3.0 * g_xxzz_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxzzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxxxy_0[i] = 3.0 * g_xxzz_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxxxz_0[i] = g_zzzz_0_xxxxxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxxyy_0[i] = 3.0 * g_xxzz_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxxyz_0[i] = g_zzzz_0_xxxxxxyz_0[i] * fbe_0 - g_zzzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxxzz_0[i] = g_zzzz_0_xxxxxxzz_0[i] * fbe_0 - g_zzzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxyyy_0[i] = 3.0 * g_xxzz_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxyyz_0[i] = g_zzzz_0_xxxxxyyz_0[i] * fbe_0 - g_zzzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxyzz_0[i] = g_zzzz_0_xxxxxyzz_0[i] * fbe_0 - g_zzzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxzzz_0[i] = g_zzzz_0_xxxxxzzz_0[i] * fbe_0 - g_zzzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxyyyy_0[i] = 3.0 * g_xxzz_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxyyyz_0[i] = g_zzzz_0_xxxxyyyz_0[i] * fbe_0 - g_zzzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxyyzz_0[i] = g_zzzz_0_xxxxyyzz_0[i] * fbe_0 - g_zzzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxyzzz_0[i] = g_zzzz_0_xxxxyzzz_0[i] * fbe_0 - g_zzzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxzzzz_0[i] = g_zzzz_0_xxxxzzzz_0[i] * fbe_0 - g_zzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyyyyy_0[i] = 3.0 * g_xxzz_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxyyyyz_0[i] = g_zzzz_0_xxxyyyyz_0[i] * fbe_0 - g_zzzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyyyzz_0[i] = g_zzzz_0_xxxyyyzz_0[i] * fbe_0 - g_zzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyyzzz_0[i] = g_zzzz_0_xxxyyzzz_0[i] * fbe_0 - g_zzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyzzzz_0[i] = g_zzzz_0_xxxyzzzz_0[i] * fbe_0 - g_zzzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxzzzzz_0[i] = g_zzzz_0_xxxzzzzz_0[i] * fbe_0 - g_zzzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyyyyy_0[i] = 3.0 * g_xxzz_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxyyyyyz_0[i] = g_zzzz_0_xxyyyyyz_0[i] * fbe_0 - g_zzzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyyyzz_0[i] = g_zzzz_0_xxyyyyzz_0[i] * fbe_0 - g_zzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyyzzz_0[i] = g_zzzz_0_xxyyyzzz_0[i] * fbe_0 - g_zzzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyzzzz_0[i] = g_zzzz_0_xxyyzzzz_0[i] * fbe_0 - g_zzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyzzzzz_0[i] = g_zzzz_0_xxyzzzzz_0[i] * fbe_0 - g_zzzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxzzzzzz_0[i] = g_zzzz_0_xxzzzzzz_0[i] * fbe_0 - g_zzzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyyyyy_0[i] = 3.0 * g_xxzz_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xyyyyyyz_0[i] = g_zzzz_0_xyyyyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyyyzz_0[i] = g_zzzz_0_xyyyyyzz_0[i] * fbe_0 - g_zzzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyyzzz_0[i] = g_zzzz_0_xyyyyzzz_0[i] * fbe_0 - g_zzzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyzzzz_0[i] = g_zzzz_0_xyyyzzzz_0[i] * fbe_0 - g_zzzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyzzzzz_0[i] = g_zzzz_0_xyyzzzzz_0[i] * fbe_0 - g_zzzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyzzzzzz_0[i] = g_zzzz_0_xyzzzzzz_0[i] * fbe_0 - g_zzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xzzzzzzz_0[i] = g_zzzz_0_xzzzzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyyyy_0[i] = g_zzzz_0_yyyyyyyy_0[i] * fbe_0 - g_zzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyyyz_0[i] = g_zzzz_0_yyyyyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyyzz_0[i] = g_zzzz_0_yyyyyyzz_0[i] * fbe_0 - g_zzzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyzzz_0[i] = g_zzzz_0_yyyyyzzz_0[i] * fbe_0 - g_zzzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyzzzz_0[i] = g_zzzz_0_yyyyzzzz_0[i] * fbe_0 - g_zzzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyzzzzz_0[i] = g_zzzz_0_yyyzzzzz_0[i] * fbe_0 - g_zzzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyzzzzzz_0[i] = g_zzzz_0_yyzzzzzz_0[i] * fbe_0 - g_zzzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yzzzzzzz_0[i] = g_zzzz_0_yzzzzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_zzzzzzzz_0[i] = g_zzzz_0_zzzzzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 675-720 components of targeted buffer : ISL

    auto g_xyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 675);

    auto g_xyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 676);

    auto g_xyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 677);

    auto g_xyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 678);

    auto g_xyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 679);

    auto g_xyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 680);

    auto g_xyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 681);

    auto g_xyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 682);

    auto g_xyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 683);

    auto g_xyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 684);

    auto g_xyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 685);

    auto g_xyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 686);

    auto g_xyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 687);

    auto g_xyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 688);

    auto g_xyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 689);

    auto g_xyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 690);

    auto g_xyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 691);

    auto g_xyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 692);

    auto g_xyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 693);

    auto g_xyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 694);

    auto g_xyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 695);

    auto g_xyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 696);

    auto g_xyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 697);

    auto g_xyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 698);

    auto g_xyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 699);

    auto g_xyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 700);

    auto g_xyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 701);

    auto g_xyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 702);

    auto g_xyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 703);

    auto g_xyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 704);

    auto g_xyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 705);

    auto g_xyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 706);

    auto g_xyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 707);

    auto g_xyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 708);

    auto g_xyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 709);

    auto g_xyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 710);

    auto g_xyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 711);

    auto g_xyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 712);

    auto g_xyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 713);

    auto g_xyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 714);

    auto g_xyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 715);

    auto g_xyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 716);

    auto g_xyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 717);

    auto g_xyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 718);

    auto g_xyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 719);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxxxxx_0, g_xyyyyy_0_xxxxxxxy_0, g_xyyyyy_0_xxxxxxxz_0, g_xyyyyy_0_xxxxxxyy_0, g_xyyyyy_0_xxxxxxyz_0, g_xyyyyy_0_xxxxxxzz_0, g_xyyyyy_0_xxxxxyyy_0, g_xyyyyy_0_xxxxxyyz_0, g_xyyyyy_0_xxxxxyzz_0, g_xyyyyy_0_xxxxxzzz_0, g_xyyyyy_0_xxxxyyyy_0, g_xyyyyy_0_xxxxyyyz_0, g_xyyyyy_0_xxxxyyzz_0, g_xyyyyy_0_xxxxyzzz_0, g_xyyyyy_0_xxxxzzzz_0, g_xyyyyy_0_xxxyyyyy_0, g_xyyyyy_0_xxxyyyyz_0, g_xyyyyy_0_xxxyyyzz_0, g_xyyyyy_0_xxxyyzzz_0, g_xyyyyy_0_xxxyzzzz_0, g_xyyyyy_0_xxxzzzzz_0, g_xyyyyy_0_xxyyyyyy_0, g_xyyyyy_0_xxyyyyyz_0, g_xyyyyy_0_xxyyyyzz_0, g_xyyyyy_0_xxyyyzzz_0, g_xyyyyy_0_xxyyzzzz_0, g_xyyyyy_0_xxyzzzzz_0, g_xyyyyy_0_xxzzzzzz_0, g_xyyyyy_0_xyyyyyyy_0, g_xyyyyy_0_xyyyyyyz_0, g_xyyyyy_0_xyyyyyzz_0, g_xyyyyy_0_xyyyyzzz_0, g_xyyyyy_0_xyyyzzzz_0, g_xyyyyy_0_xyyzzzzz_0, g_xyyyyy_0_xyzzzzzz_0, g_xyyyyy_0_xzzzzzzz_0, g_xyyyyy_0_yyyyyyyy_0, g_xyyyyy_0_yyyyyyyz_0, g_xyyyyy_0_yyyyyyzz_0, g_xyyyyy_0_yyyyyzzz_0, g_xyyyyy_0_yyyyzzzz_0, g_xyyyyy_0_yyyzzzzz_0, g_xyyyyy_0_yyzzzzzz_0, g_xyyyyy_0_yzzzzzzz_0, g_xyyyyy_0_zzzzzzzz_0, g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxxx_1, g_yyyyy_0_xxxxxxxy_1, g_yyyyy_0_xxxxxxxz_1, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxyy_1, g_yyyyy_0_xxxxxxyz_1, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxxzz_1, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyyy_1, g_yyyyy_0_xxxxxyyz_1, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxyzz_1, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxxzzz_1, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyyy_1, g_yyyyy_0_xxxxyyyz_1, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyyzz_1, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxyzzz_1, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxxzzzz_1, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyyy_1, g_yyyyy_0_xxxyyyyz_1, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyyzz_1, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyyzzz_1, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxyzzzz_1, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxxzzzzz_1, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyyy_1, g_yyyyy_0_xxyyyyyz_1, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyyzz_1, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyyzzz_1, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyyzzzz_1, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxyzzzzz_1, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xxzzzzzz_1, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyyy_1, g_yyyyy_0_xyyyyyyz_1, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyyzz_1, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyyzzz_1, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyyzzzz_1, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyyzzzzz_1, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xyzzzzzz_1, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_xzzzzzzz_1, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyyy_1, g_yyyyy_0_yyyyyyyz_1, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyyzz_1, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyyzzz_1, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyyzzzz_1, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyyzzzzz_1, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yyzzzzzz_1, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_yzzzzzzz_1, g_yyyyy_0_zzzzzzz_1, g_yyyyy_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xxxxxxxx_0[i] = 8.0 * g_yyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxxy_0[i] = 7.0 * g_yyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxxz_0[i] = 7.0 * g_yyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxyy_0[i] = 6.0 * g_yyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxyz_0[i] = 6.0 * g_yyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxzz_0[i] = 6.0 * g_yyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxyyy_0[i] = 5.0 * g_yyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxyyz_0[i] = 5.0 * g_yyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxyzz_0[i] = 5.0 * g_yyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxzzz_0[i] = 5.0 * g_yyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyyyy_0[i] = 4.0 * g_yyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyyyz_0[i] = 4.0 * g_yyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyyzz_0[i] = 4.0 * g_yyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyzzz_0[i] = 4.0 * g_yyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxzzzz_0[i] = 4.0 * g_yyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyyyy_0[i] = 3.0 * g_yyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyyyz_0[i] = 3.0 * g_yyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyyzz_0[i] = 3.0 * g_yyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyzzz_0[i] = 3.0 * g_yyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyzzzz_0[i] = 3.0 * g_yyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxzzzzz_0[i] = 3.0 * g_yyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyyyy_0[i] = 2.0 * g_yyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyyyz_0[i] = 2.0 * g_yyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyyzz_0[i] = 2.0 * g_yyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyzzz_0[i] = 2.0 * g_yyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyzzzz_0[i] = 2.0 * g_yyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyzzzzz_0[i] = 2.0 * g_yyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxzzzzzz_0[i] = 2.0 * g_yyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyyyy_0[i] = g_yyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyyyz_0[i] = g_yyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyyzz_0[i] = g_yyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyzzz_0[i] = g_yyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyzzzz_0[i] = g_yyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyzzzzz_0[i] = g_yyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyzzzzzz_0[i] = g_yyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xzzzzzzz_0[i] = g_yyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyyyy_0[i] = g_yyyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyyyz_0[i] = g_yyyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyyzz_0[i] = g_yyyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyzzz_0[i] = g_yyyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyzzzz_0[i] = g_yyyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyzzzzz_0[i] = g_yyyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyzzzzzz_0[i] = g_yyyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yzzzzzzz_0[i] = g_yyyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_zzzzzzzz_0[i] = g_yyyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 720-765 components of targeted buffer : ISL

    auto g_xyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 720);

    auto g_xyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 721);

    auto g_xyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 722);

    auto g_xyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 723);

    auto g_xyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 724);

    auto g_xyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 725);

    auto g_xyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 726);

    auto g_xyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 727);

    auto g_xyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 728);

    auto g_xyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 729);

    auto g_xyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 730);

    auto g_xyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 731);

    auto g_xyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 732);

    auto g_xyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 733);

    auto g_xyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 734);

    auto g_xyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 735);

    auto g_xyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 736);

    auto g_xyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 737);

    auto g_xyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 738);

    auto g_xyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 739);

    auto g_xyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 740);

    auto g_xyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 741);

    auto g_xyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 742);

    auto g_xyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 743);

    auto g_xyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 744);

    auto g_xyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 745);

    auto g_xyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 746);

    auto g_xyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 747);

    auto g_xyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 748);

    auto g_xyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 749);

    auto g_xyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 750);

    auto g_xyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 751);

    auto g_xyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 752);

    auto g_xyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 753);

    auto g_xyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 754);

    auto g_xyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 755);

    auto g_xyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 756);

    auto g_xyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 757);

    auto g_xyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 758);

    auto g_xyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 759);

    auto g_xyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 760);

    auto g_xyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 761);

    auto g_xyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 762);

    auto g_xyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 763);

    auto g_xyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 764);

    #pragma omp simd aligned(g_xyyyy_0_xxxxxxxx_1, g_xyyyy_0_xxxxxxxy_1, g_xyyyy_0_xxxxxxyy_1, g_xyyyy_0_xxxxxyyy_1, g_xyyyy_0_xxxxyyyy_1, g_xyyyy_0_xxxyyyyy_1, g_xyyyy_0_xxyyyyyy_1, g_xyyyy_0_xyyyyyyy_1, g_xyyyyz_0_xxxxxxxx_0, g_xyyyyz_0_xxxxxxxy_0, g_xyyyyz_0_xxxxxxxz_0, g_xyyyyz_0_xxxxxxyy_0, g_xyyyyz_0_xxxxxxyz_0, g_xyyyyz_0_xxxxxxzz_0, g_xyyyyz_0_xxxxxyyy_0, g_xyyyyz_0_xxxxxyyz_0, g_xyyyyz_0_xxxxxyzz_0, g_xyyyyz_0_xxxxxzzz_0, g_xyyyyz_0_xxxxyyyy_0, g_xyyyyz_0_xxxxyyyz_0, g_xyyyyz_0_xxxxyyzz_0, g_xyyyyz_0_xxxxyzzz_0, g_xyyyyz_0_xxxxzzzz_0, g_xyyyyz_0_xxxyyyyy_0, g_xyyyyz_0_xxxyyyyz_0, g_xyyyyz_0_xxxyyyzz_0, g_xyyyyz_0_xxxyyzzz_0, g_xyyyyz_0_xxxyzzzz_0, g_xyyyyz_0_xxxzzzzz_0, g_xyyyyz_0_xxyyyyyy_0, g_xyyyyz_0_xxyyyyyz_0, g_xyyyyz_0_xxyyyyzz_0, g_xyyyyz_0_xxyyyzzz_0, g_xyyyyz_0_xxyyzzzz_0, g_xyyyyz_0_xxyzzzzz_0, g_xyyyyz_0_xxzzzzzz_0, g_xyyyyz_0_xyyyyyyy_0, g_xyyyyz_0_xyyyyyyz_0, g_xyyyyz_0_xyyyyyzz_0, g_xyyyyz_0_xyyyyzzz_0, g_xyyyyz_0_xyyyzzzz_0, g_xyyyyz_0_xyyzzzzz_0, g_xyyyyz_0_xyzzzzzz_0, g_xyyyyz_0_xzzzzzzz_0, g_xyyyyz_0_yyyyyyyy_0, g_xyyyyz_0_yyyyyyyz_0, g_xyyyyz_0_yyyyyyzz_0, g_xyyyyz_0_yyyyyzzz_0, g_xyyyyz_0_yyyyzzzz_0, g_xyyyyz_0_yyyzzzzz_0, g_xyyyyz_0_yyzzzzzz_0, g_xyyyyz_0_yzzzzzzz_0, g_xyyyyz_0_zzzzzzzz_0, g_yyyyz_0_xxxxxxxz_1, g_yyyyz_0_xxxxxxyz_1, g_yyyyz_0_xxxxxxz_1, g_yyyyz_0_xxxxxxzz_1, g_yyyyz_0_xxxxxyyz_1, g_yyyyz_0_xxxxxyz_1, g_yyyyz_0_xxxxxyzz_1, g_yyyyz_0_xxxxxzz_1, g_yyyyz_0_xxxxxzzz_1, g_yyyyz_0_xxxxyyyz_1, g_yyyyz_0_xxxxyyz_1, g_yyyyz_0_xxxxyyzz_1, g_yyyyz_0_xxxxyzz_1, g_yyyyz_0_xxxxyzzz_1, g_yyyyz_0_xxxxzzz_1, g_yyyyz_0_xxxxzzzz_1, g_yyyyz_0_xxxyyyyz_1, g_yyyyz_0_xxxyyyz_1, g_yyyyz_0_xxxyyyzz_1, g_yyyyz_0_xxxyyzz_1, g_yyyyz_0_xxxyyzzz_1, g_yyyyz_0_xxxyzzz_1, g_yyyyz_0_xxxyzzzz_1, g_yyyyz_0_xxxzzzz_1, g_yyyyz_0_xxxzzzzz_1, g_yyyyz_0_xxyyyyyz_1, g_yyyyz_0_xxyyyyz_1, g_yyyyz_0_xxyyyyzz_1, g_yyyyz_0_xxyyyzz_1, g_yyyyz_0_xxyyyzzz_1, g_yyyyz_0_xxyyzzz_1, g_yyyyz_0_xxyyzzzz_1, g_yyyyz_0_xxyzzzz_1, g_yyyyz_0_xxyzzzzz_1, g_yyyyz_0_xxzzzzz_1, g_yyyyz_0_xxzzzzzz_1, g_yyyyz_0_xyyyyyyz_1, g_yyyyz_0_xyyyyyz_1, g_yyyyz_0_xyyyyyzz_1, g_yyyyz_0_xyyyyzz_1, g_yyyyz_0_xyyyyzzz_1, g_yyyyz_0_xyyyzzz_1, g_yyyyz_0_xyyyzzzz_1, g_yyyyz_0_xyyzzzz_1, g_yyyyz_0_xyyzzzzz_1, g_yyyyz_0_xyzzzzz_1, g_yyyyz_0_xyzzzzzz_1, g_yyyyz_0_xzzzzzz_1, g_yyyyz_0_xzzzzzzz_1, g_yyyyz_0_yyyyyyyy_1, g_yyyyz_0_yyyyyyyz_1, g_yyyyz_0_yyyyyyz_1, g_yyyyz_0_yyyyyyzz_1, g_yyyyz_0_yyyyyzz_1, g_yyyyz_0_yyyyyzzz_1, g_yyyyz_0_yyyyzzz_1, g_yyyyz_0_yyyyzzzz_1, g_yyyyz_0_yyyzzzz_1, g_yyyyz_0_yyyzzzzz_1, g_yyyyz_0_yyzzzzz_1, g_yyyyz_0_yyzzzzzz_1, g_yyyyz_0_yzzzzzz_1, g_yyyyz_0_yzzzzzzz_1, g_yyyyz_0_zzzzzzz_1, g_yyyyz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xxxxxxxx_0[i] = g_xyyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxxxy_0[i] = g_xyyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxxxz_0[i] = 7.0 * g_yyyyz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxxyy_0[i] = g_xyyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxxyz_0[i] = 6.0 * g_yyyyz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxxzz_0[i] = 6.0 * g_yyyyz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxyyy_0[i] = g_xyyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxyyz_0[i] = 5.0 * g_yyyyz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxyzz_0[i] = 5.0 * g_yyyyz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxzzz_0[i] = 5.0 * g_yyyyz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxyyyy_0[i] = g_xyyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxyyyz_0[i] = 4.0 * g_yyyyz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxyyzz_0[i] = 4.0 * g_yyyyz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxyzzz_0[i] = 4.0 * g_yyyyz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyyyz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyyyyy_0[i] = g_xyyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxyyyyz_0[i] = 3.0 * g_yyyyz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyyyzz_0[i] = 3.0 * g_yyyyz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyyyz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyzzzz_0[i] = 3.0 * g_yyyyz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxzzzzz_0[i] = 3.0 * g_yyyyz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyyyyy_0[i] = g_xyyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxyyyyyz_0[i] = 2.0 * g_yyyyz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyyyz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyyzzz_0[i] = 2.0 * g_yyyyz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyzzzz_0[i] = 2.0 * g_yyyyz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyzzzzz_0[i] = 2.0 * g_yyyyz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxzzzzzz_0[i] = 2.0 * g_yyyyz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyyyyy_0[i] = g_xyyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xyyyyyyz_0[i] = g_yyyyz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyyyzz_0[i] = g_yyyyz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyyzzz_0[i] = g_yyyyz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyzzzz_0[i] = g_yyyyz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyzzzzz_0[i] = g_yyyyz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyzzzzzz_0[i] = g_yyyyz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xzzzzzzz_0[i] = g_yyyyz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyyyy_0[i] = g_yyyyz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyyyz_0[i] = g_yyyyz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyyzz_0[i] = g_yyyyz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyzzz_0[i] = g_yyyyz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyzzzz_0[i] = g_yyyyz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyzzzzz_0[i] = g_yyyyz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyzzzzzz_0[i] = g_yyyyz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yzzzzzzz_0[i] = g_yyyyz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_zzzzzzzz_0[i] = g_yyyyz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 765-810 components of targeted buffer : ISL

    auto g_xyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 765);

    auto g_xyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 766);

    auto g_xyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 767);

    auto g_xyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 768);

    auto g_xyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 769);

    auto g_xyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 770);

    auto g_xyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 771);

    auto g_xyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 772);

    auto g_xyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 773);

    auto g_xyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 774);

    auto g_xyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 775);

    auto g_xyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 776);

    auto g_xyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 777);

    auto g_xyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 778);

    auto g_xyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 779);

    auto g_xyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 780);

    auto g_xyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 781);

    auto g_xyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 782);

    auto g_xyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 783);

    auto g_xyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 784);

    auto g_xyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 785);

    auto g_xyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 786);

    auto g_xyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 787);

    auto g_xyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 788);

    auto g_xyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 789);

    auto g_xyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 790);

    auto g_xyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 791);

    auto g_xyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 792);

    auto g_xyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 793);

    auto g_xyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 794);

    auto g_xyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 795);

    auto g_xyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 796);

    auto g_xyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 797);

    auto g_xyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 798);

    auto g_xyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 799);

    auto g_xyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 800);

    auto g_xyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 801);

    auto g_xyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 802);

    auto g_xyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 803);

    auto g_xyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 804);

    auto g_xyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 805);

    auto g_xyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 806);

    auto g_xyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 807);

    auto g_xyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 808);

    auto g_xyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 809);

    #pragma omp simd aligned(g_xyyyzz_0_xxxxxxxx_0, g_xyyyzz_0_xxxxxxxy_0, g_xyyyzz_0_xxxxxxxz_0, g_xyyyzz_0_xxxxxxyy_0, g_xyyyzz_0_xxxxxxyz_0, g_xyyyzz_0_xxxxxxzz_0, g_xyyyzz_0_xxxxxyyy_0, g_xyyyzz_0_xxxxxyyz_0, g_xyyyzz_0_xxxxxyzz_0, g_xyyyzz_0_xxxxxzzz_0, g_xyyyzz_0_xxxxyyyy_0, g_xyyyzz_0_xxxxyyyz_0, g_xyyyzz_0_xxxxyyzz_0, g_xyyyzz_0_xxxxyzzz_0, g_xyyyzz_0_xxxxzzzz_0, g_xyyyzz_0_xxxyyyyy_0, g_xyyyzz_0_xxxyyyyz_0, g_xyyyzz_0_xxxyyyzz_0, g_xyyyzz_0_xxxyyzzz_0, g_xyyyzz_0_xxxyzzzz_0, g_xyyyzz_0_xxxzzzzz_0, g_xyyyzz_0_xxyyyyyy_0, g_xyyyzz_0_xxyyyyyz_0, g_xyyyzz_0_xxyyyyzz_0, g_xyyyzz_0_xxyyyzzz_0, g_xyyyzz_0_xxyyzzzz_0, g_xyyyzz_0_xxyzzzzz_0, g_xyyyzz_0_xxzzzzzz_0, g_xyyyzz_0_xyyyyyyy_0, g_xyyyzz_0_xyyyyyyz_0, g_xyyyzz_0_xyyyyyzz_0, g_xyyyzz_0_xyyyyzzz_0, g_xyyyzz_0_xyyyzzzz_0, g_xyyyzz_0_xyyzzzzz_0, g_xyyyzz_0_xyzzzzzz_0, g_xyyyzz_0_xzzzzzzz_0, g_xyyyzz_0_yyyyyyyy_0, g_xyyyzz_0_yyyyyyyz_0, g_xyyyzz_0_yyyyyyzz_0, g_xyyyzz_0_yyyyyzzz_0, g_xyyyzz_0_yyyyzzzz_0, g_xyyyzz_0_yyyzzzzz_0, g_xyyyzz_0_yyzzzzzz_0, g_xyyyzz_0_yzzzzzzz_0, g_xyyyzz_0_zzzzzzzz_0, g_yyyzz_0_xxxxxxx_1, g_yyyzz_0_xxxxxxxx_1, g_yyyzz_0_xxxxxxxy_1, g_yyyzz_0_xxxxxxxz_1, g_yyyzz_0_xxxxxxy_1, g_yyyzz_0_xxxxxxyy_1, g_yyyzz_0_xxxxxxyz_1, g_yyyzz_0_xxxxxxz_1, g_yyyzz_0_xxxxxxzz_1, g_yyyzz_0_xxxxxyy_1, g_yyyzz_0_xxxxxyyy_1, g_yyyzz_0_xxxxxyyz_1, g_yyyzz_0_xxxxxyz_1, g_yyyzz_0_xxxxxyzz_1, g_yyyzz_0_xxxxxzz_1, g_yyyzz_0_xxxxxzzz_1, g_yyyzz_0_xxxxyyy_1, g_yyyzz_0_xxxxyyyy_1, g_yyyzz_0_xxxxyyyz_1, g_yyyzz_0_xxxxyyz_1, g_yyyzz_0_xxxxyyzz_1, g_yyyzz_0_xxxxyzz_1, g_yyyzz_0_xxxxyzzz_1, g_yyyzz_0_xxxxzzz_1, g_yyyzz_0_xxxxzzzz_1, g_yyyzz_0_xxxyyyy_1, g_yyyzz_0_xxxyyyyy_1, g_yyyzz_0_xxxyyyyz_1, g_yyyzz_0_xxxyyyz_1, g_yyyzz_0_xxxyyyzz_1, g_yyyzz_0_xxxyyzz_1, g_yyyzz_0_xxxyyzzz_1, g_yyyzz_0_xxxyzzz_1, g_yyyzz_0_xxxyzzzz_1, g_yyyzz_0_xxxzzzz_1, g_yyyzz_0_xxxzzzzz_1, g_yyyzz_0_xxyyyyy_1, g_yyyzz_0_xxyyyyyy_1, g_yyyzz_0_xxyyyyyz_1, g_yyyzz_0_xxyyyyz_1, g_yyyzz_0_xxyyyyzz_1, g_yyyzz_0_xxyyyzz_1, g_yyyzz_0_xxyyyzzz_1, g_yyyzz_0_xxyyzzz_1, g_yyyzz_0_xxyyzzzz_1, g_yyyzz_0_xxyzzzz_1, g_yyyzz_0_xxyzzzzz_1, g_yyyzz_0_xxzzzzz_1, g_yyyzz_0_xxzzzzzz_1, g_yyyzz_0_xyyyyyy_1, g_yyyzz_0_xyyyyyyy_1, g_yyyzz_0_xyyyyyyz_1, g_yyyzz_0_xyyyyyz_1, g_yyyzz_0_xyyyyyzz_1, g_yyyzz_0_xyyyyzz_1, g_yyyzz_0_xyyyyzzz_1, g_yyyzz_0_xyyyzzz_1, g_yyyzz_0_xyyyzzzz_1, g_yyyzz_0_xyyzzzz_1, g_yyyzz_0_xyyzzzzz_1, g_yyyzz_0_xyzzzzz_1, g_yyyzz_0_xyzzzzzz_1, g_yyyzz_0_xzzzzzz_1, g_yyyzz_0_xzzzzzzz_1, g_yyyzz_0_yyyyyyy_1, g_yyyzz_0_yyyyyyyy_1, g_yyyzz_0_yyyyyyyz_1, g_yyyzz_0_yyyyyyz_1, g_yyyzz_0_yyyyyyzz_1, g_yyyzz_0_yyyyyzz_1, g_yyyzz_0_yyyyyzzz_1, g_yyyzz_0_yyyyzzz_1, g_yyyzz_0_yyyyzzzz_1, g_yyyzz_0_yyyzzzz_1, g_yyyzz_0_yyyzzzzz_1, g_yyyzz_0_yyzzzzz_1, g_yyyzz_0_yyzzzzzz_1, g_yyyzz_0_yzzzzzz_1, g_yyyzz_0_yzzzzzzz_1, g_yyyzz_0_zzzzzzz_1, g_yyyzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xxxxxxxx_0[i] = 8.0 * g_yyyzz_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxxy_0[i] = 7.0 * g_yyyzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxxz_0[i] = 7.0 * g_yyyzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxyy_0[i] = 6.0 * g_yyyzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxyz_0[i] = 6.0 * g_yyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxzz_0[i] = 6.0 * g_yyyzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxyyy_0[i] = 5.0 * g_yyyzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxyyz_0[i] = 5.0 * g_yyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxyzz_0[i] = 5.0 * g_yyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxzzz_0[i] = 5.0 * g_yyyzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyyyy_0[i] = 4.0 * g_yyyzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyyyz_0[i] = 4.0 * g_yyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyyzz_0[i] = 4.0 * g_yyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyzzz_0[i] = 4.0 * g_yyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxzzzz_0[i] = 4.0 * g_yyyzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyyyy_0[i] = 3.0 * g_yyyzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyyyz_0[i] = 3.0 * g_yyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyyzz_0[i] = 3.0 * g_yyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyzzz_0[i] = 3.0 * g_yyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyzzzz_0[i] = 3.0 * g_yyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxzzzzz_0[i] = 3.0 * g_yyyzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyyyy_0[i] = 2.0 * g_yyyzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyyyz_0[i] = 2.0 * g_yyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyyzz_0[i] = 2.0 * g_yyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyzzz_0[i] = 2.0 * g_yyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyzzzz_0[i] = 2.0 * g_yyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyzzzzz_0[i] = 2.0 * g_yyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxzzzzzz_0[i] = 2.0 * g_yyyzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyyyy_0[i] = g_yyyzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyyyz_0[i] = g_yyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyyzz_0[i] = g_yyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyzzz_0[i] = g_yyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyzzzz_0[i] = g_yyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyzzzzz_0[i] = g_yyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyzzzzzz_0[i] = g_yyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xzzzzzzz_0[i] = g_yyyzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyyyy_0[i] = g_yyyzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyyyz_0[i] = g_yyyzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyyzz_0[i] = g_yyyzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyzzz_0[i] = g_yyyzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyzzzz_0[i] = g_yyyzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyzzzzz_0[i] = g_yyyzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyzzzzzz_0[i] = g_yyyzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yzzzzzzz_0[i] = g_yyyzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_zzzzzzzz_0[i] = g_yyyzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 810-855 components of targeted buffer : ISL

    auto g_xyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 810);

    auto g_xyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 811);

    auto g_xyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 812);

    auto g_xyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 813);

    auto g_xyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 814);

    auto g_xyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 815);

    auto g_xyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 816);

    auto g_xyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 817);

    auto g_xyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 818);

    auto g_xyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 819);

    auto g_xyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 820);

    auto g_xyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 821);

    auto g_xyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 822);

    auto g_xyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 823);

    auto g_xyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 824);

    auto g_xyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 825);

    auto g_xyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 826);

    auto g_xyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 827);

    auto g_xyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 828);

    auto g_xyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 829);

    auto g_xyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 830);

    auto g_xyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 831);

    auto g_xyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 832);

    auto g_xyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 833);

    auto g_xyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 834);

    auto g_xyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 835);

    auto g_xyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 836);

    auto g_xyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 837);

    auto g_xyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 838);

    auto g_xyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 839);

    auto g_xyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 840);

    auto g_xyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 841);

    auto g_xyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 842);

    auto g_xyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 843);

    auto g_xyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 844);

    auto g_xyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 845);

    auto g_xyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 846);

    auto g_xyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 847);

    auto g_xyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 848);

    auto g_xyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 849);

    auto g_xyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 850);

    auto g_xyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 851);

    auto g_xyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 852);

    auto g_xyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 853);

    auto g_xyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 854);

    #pragma omp simd aligned(g_xyyzzz_0_xxxxxxxx_0, g_xyyzzz_0_xxxxxxxy_0, g_xyyzzz_0_xxxxxxxz_0, g_xyyzzz_0_xxxxxxyy_0, g_xyyzzz_0_xxxxxxyz_0, g_xyyzzz_0_xxxxxxzz_0, g_xyyzzz_0_xxxxxyyy_0, g_xyyzzz_0_xxxxxyyz_0, g_xyyzzz_0_xxxxxyzz_0, g_xyyzzz_0_xxxxxzzz_0, g_xyyzzz_0_xxxxyyyy_0, g_xyyzzz_0_xxxxyyyz_0, g_xyyzzz_0_xxxxyyzz_0, g_xyyzzz_0_xxxxyzzz_0, g_xyyzzz_0_xxxxzzzz_0, g_xyyzzz_0_xxxyyyyy_0, g_xyyzzz_0_xxxyyyyz_0, g_xyyzzz_0_xxxyyyzz_0, g_xyyzzz_0_xxxyyzzz_0, g_xyyzzz_0_xxxyzzzz_0, g_xyyzzz_0_xxxzzzzz_0, g_xyyzzz_0_xxyyyyyy_0, g_xyyzzz_0_xxyyyyyz_0, g_xyyzzz_0_xxyyyyzz_0, g_xyyzzz_0_xxyyyzzz_0, g_xyyzzz_0_xxyyzzzz_0, g_xyyzzz_0_xxyzzzzz_0, g_xyyzzz_0_xxzzzzzz_0, g_xyyzzz_0_xyyyyyyy_0, g_xyyzzz_0_xyyyyyyz_0, g_xyyzzz_0_xyyyyyzz_0, g_xyyzzz_0_xyyyyzzz_0, g_xyyzzz_0_xyyyzzzz_0, g_xyyzzz_0_xyyzzzzz_0, g_xyyzzz_0_xyzzzzzz_0, g_xyyzzz_0_xzzzzzzz_0, g_xyyzzz_0_yyyyyyyy_0, g_xyyzzz_0_yyyyyyyz_0, g_xyyzzz_0_yyyyyyzz_0, g_xyyzzz_0_yyyyyzzz_0, g_xyyzzz_0_yyyyzzzz_0, g_xyyzzz_0_yyyzzzzz_0, g_xyyzzz_0_yyzzzzzz_0, g_xyyzzz_0_yzzzzzzz_0, g_xyyzzz_0_zzzzzzzz_0, g_yyzzz_0_xxxxxxx_1, g_yyzzz_0_xxxxxxxx_1, g_yyzzz_0_xxxxxxxy_1, g_yyzzz_0_xxxxxxxz_1, g_yyzzz_0_xxxxxxy_1, g_yyzzz_0_xxxxxxyy_1, g_yyzzz_0_xxxxxxyz_1, g_yyzzz_0_xxxxxxz_1, g_yyzzz_0_xxxxxxzz_1, g_yyzzz_0_xxxxxyy_1, g_yyzzz_0_xxxxxyyy_1, g_yyzzz_0_xxxxxyyz_1, g_yyzzz_0_xxxxxyz_1, g_yyzzz_0_xxxxxyzz_1, g_yyzzz_0_xxxxxzz_1, g_yyzzz_0_xxxxxzzz_1, g_yyzzz_0_xxxxyyy_1, g_yyzzz_0_xxxxyyyy_1, g_yyzzz_0_xxxxyyyz_1, g_yyzzz_0_xxxxyyz_1, g_yyzzz_0_xxxxyyzz_1, g_yyzzz_0_xxxxyzz_1, g_yyzzz_0_xxxxyzzz_1, g_yyzzz_0_xxxxzzz_1, g_yyzzz_0_xxxxzzzz_1, g_yyzzz_0_xxxyyyy_1, g_yyzzz_0_xxxyyyyy_1, g_yyzzz_0_xxxyyyyz_1, g_yyzzz_0_xxxyyyz_1, g_yyzzz_0_xxxyyyzz_1, g_yyzzz_0_xxxyyzz_1, g_yyzzz_0_xxxyyzzz_1, g_yyzzz_0_xxxyzzz_1, g_yyzzz_0_xxxyzzzz_1, g_yyzzz_0_xxxzzzz_1, g_yyzzz_0_xxxzzzzz_1, g_yyzzz_0_xxyyyyy_1, g_yyzzz_0_xxyyyyyy_1, g_yyzzz_0_xxyyyyyz_1, g_yyzzz_0_xxyyyyz_1, g_yyzzz_0_xxyyyyzz_1, g_yyzzz_0_xxyyyzz_1, g_yyzzz_0_xxyyyzzz_1, g_yyzzz_0_xxyyzzz_1, g_yyzzz_0_xxyyzzzz_1, g_yyzzz_0_xxyzzzz_1, g_yyzzz_0_xxyzzzzz_1, g_yyzzz_0_xxzzzzz_1, g_yyzzz_0_xxzzzzzz_1, g_yyzzz_0_xyyyyyy_1, g_yyzzz_0_xyyyyyyy_1, g_yyzzz_0_xyyyyyyz_1, g_yyzzz_0_xyyyyyz_1, g_yyzzz_0_xyyyyyzz_1, g_yyzzz_0_xyyyyzz_1, g_yyzzz_0_xyyyyzzz_1, g_yyzzz_0_xyyyzzz_1, g_yyzzz_0_xyyyzzzz_1, g_yyzzz_0_xyyzzzz_1, g_yyzzz_0_xyyzzzzz_1, g_yyzzz_0_xyzzzzz_1, g_yyzzz_0_xyzzzzzz_1, g_yyzzz_0_xzzzzzz_1, g_yyzzz_0_xzzzzzzz_1, g_yyzzz_0_yyyyyyy_1, g_yyzzz_0_yyyyyyyy_1, g_yyzzz_0_yyyyyyyz_1, g_yyzzz_0_yyyyyyz_1, g_yyzzz_0_yyyyyyzz_1, g_yyzzz_0_yyyyyzz_1, g_yyzzz_0_yyyyyzzz_1, g_yyzzz_0_yyyyzzz_1, g_yyzzz_0_yyyyzzzz_1, g_yyzzz_0_yyyzzzz_1, g_yyzzz_0_yyyzzzzz_1, g_yyzzz_0_yyzzzzz_1, g_yyzzz_0_yyzzzzzz_1, g_yyzzz_0_yzzzzzz_1, g_yyzzz_0_yzzzzzzz_1, g_yyzzz_0_zzzzzzz_1, g_yyzzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xxxxxxxx_0[i] = 8.0 * g_yyzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxxy_0[i] = 7.0 * g_yyzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxxz_0[i] = 7.0 * g_yyzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxyy_0[i] = 6.0 * g_yyzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxyz_0[i] = 6.0 * g_yyzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxzz_0[i] = 6.0 * g_yyzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxyyy_0[i] = 5.0 * g_yyzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxyyz_0[i] = 5.0 * g_yyzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxyzz_0[i] = 5.0 * g_yyzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxzzz_0[i] = 5.0 * g_yyzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyyyy_0[i] = 4.0 * g_yyzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyyyz_0[i] = 4.0 * g_yyzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyyzz_0[i] = 4.0 * g_yyzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyzzz_0[i] = 4.0 * g_yyzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxzzzz_0[i] = 4.0 * g_yyzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyyyy_0[i] = 3.0 * g_yyzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyyyz_0[i] = 3.0 * g_yyzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyyzz_0[i] = 3.0 * g_yyzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyzzz_0[i] = 3.0 * g_yyzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyzzzz_0[i] = 3.0 * g_yyzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxzzzzz_0[i] = 3.0 * g_yyzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyyyy_0[i] = 2.0 * g_yyzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyyyz_0[i] = 2.0 * g_yyzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyyzz_0[i] = 2.0 * g_yyzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyzzz_0[i] = 2.0 * g_yyzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyzzzz_0[i] = 2.0 * g_yyzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyzzzzz_0[i] = 2.0 * g_yyzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxzzzzzz_0[i] = 2.0 * g_yyzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyyyy_0[i] = g_yyzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyyyz_0[i] = g_yyzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyyzz_0[i] = g_yyzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyzzz_0[i] = g_yyzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyzzzz_0[i] = g_yyzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyzzzzz_0[i] = g_yyzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyzzzzzz_0[i] = g_yyzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xzzzzzzz_0[i] = g_yyzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyyyy_0[i] = g_yyzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyyyz_0[i] = g_yyzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyyzz_0[i] = g_yyzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyzzz_0[i] = g_yyzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyzzzz_0[i] = g_yyzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyzzzzz_0[i] = g_yyzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyzzzzzz_0[i] = g_yyzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yzzzzzzz_0[i] = g_yyzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_zzzzzzzz_0[i] = g_yyzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 855-900 components of targeted buffer : ISL

    auto g_xyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 855);

    auto g_xyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 856);

    auto g_xyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 857);

    auto g_xyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 858);

    auto g_xyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 859);

    auto g_xyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 860);

    auto g_xyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 861);

    auto g_xyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 862);

    auto g_xyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 863);

    auto g_xyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 864);

    auto g_xyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 865);

    auto g_xyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 866);

    auto g_xyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 867);

    auto g_xyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 868);

    auto g_xyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 869);

    auto g_xyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 870);

    auto g_xyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 871);

    auto g_xyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 872);

    auto g_xyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 873);

    auto g_xyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 874);

    auto g_xyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 875);

    auto g_xyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 876);

    auto g_xyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 877);

    auto g_xyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 878);

    auto g_xyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 879);

    auto g_xyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 880);

    auto g_xyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 881);

    auto g_xyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 882);

    auto g_xyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 883);

    auto g_xyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 884);

    auto g_xyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 885);

    auto g_xyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 886);

    auto g_xyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 887);

    auto g_xyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 888);

    auto g_xyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 889);

    auto g_xyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 890);

    auto g_xyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 891);

    auto g_xyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 892);

    auto g_xyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 893);

    auto g_xyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 894);

    auto g_xyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 895);

    auto g_xyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 896);

    auto g_xyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 897);

    auto g_xyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 898);

    auto g_xyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 899);

    #pragma omp simd aligned(g_xyzzzz_0_xxxxxxxx_0, g_xyzzzz_0_xxxxxxxy_0, g_xyzzzz_0_xxxxxxxz_0, g_xyzzzz_0_xxxxxxyy_0, g_xyzzzz_0_xxxxxxyz_0, g_xyzzzz_0_xxxxxxzz_0, g_xyzzzz_0_xxxxxyyy_0, g_xyzzzz_0_xxxxxyyz_0, g_xyzzzz_0_xxxxxyzz_0, g_xyzzzz_0_xxxxxzzz_0, g_xyzzzz_0_xxxxyyyy_0, g_xyzzzz_0_xxxxyyyz_0, g_xyzzzz_0_xxxxyyzz_0, g_xyzzzz_0_xxxxyzzz_0, g_xyzzzz_0_xxxxzzzz_0, g_xyzzzz_0_xxxyyyyy_0, g_xyzzzz_0_xxxyyyyz_0, g_xyzzzz_0_xxxyyyzz_0, g_xyzzzz_0_xxxyyzzz_0, g_xyzzzz_0_xxxyzzzz_0, g_xyzzzz_0_xxxzzzzz_0, g_xyzzzz_0_xxyyyyyy_0, g_xyzzzz_0_xxyyyyyz_0, g_xyzzzz_0_xxyyyyzz_0, g_xyzzzz_0_xxyyyzzz_0, g_xyzzzz_0_xxyyzzzz_0, g_xyzzzz_0_xxyzzzzz_0, g_xyzzzz_0_xxzzzzzz_0, g_xyzzzz_0_xyyyyyyy_0, g_xyzzzz_0_xyyyyyyz_0, g_xyzzzz_0_xyyyyyzz_0, g_xyzzzz_0_xyyyyzzz_0, g_xyzzzz_0_xyyyzzzz_0, g_xyzzzz_0_xyyzzzzz_0, g_xyzzzz_0_xyzzzzzz_0, g_xyzzzz_0_xzzzzzzz_0, g_xyzzzz_0_yyyyyyyy_0, g_xyzzzz_0_yyyyyyyz_0, g_xyzzzz_0_yyyyyyzz_0, g_xyzzzz_0_yyyyyzzz_0, g_xyzzzz_0_yyyyzzzz_0, g_xyzzzz_0_yyyzzzzz_0, g_xyzzzz_0_yyzzzzzz_0, g_xyzzzz_0_yzzzzzzz_0, g_xyzzzz_0_zzzzzzzz_0, g_xzzzz_0_xxxxxxxx_1, g_xzzzz_0_xxxxxxxz_1, g_xzzzz_0_xxxxxxzz_1, g_xzzzz_0_xxxxxzzz_1, g_xzzzz_0_xxxxzzzz_1, g_xzzzz_0_xxxzzzzz_1, g_xzzzz_0_xxzzzzzz_1, g_xzzzz_0_xzzzzzzz_1, g_yzzzz_0_xxxxxxxy_1, g_yzzzz_0_xxxxxxy_1, g_yzzzz_0_xxxxxxyy_1, g_yzzzz_0_xxxxxxyz_1, g_yzzzz_0_xxxxxyy_1, g_yzzzz_0_xxxxxyyy_1, g_yzzzz_0_xxxxxyyz_1, g_yzzzz_0_xxxxxyz_1, g_yzzzz_0_xxxxxyzz_1, g_yzzzz_0_xxxxyyy_1, g_yzzzz_0_xxxxyyyy_1, g_yzzzz_0_xxxxyyyz_1, g_yzzzz_0_xxxxyyz_1, g_yzzzz_0_xxxxyyzz_1, g_yzzzz_0_xxxxyzz_1, g_yzzzz_0_xxxxyzzz_1, g_yzzzz_0_xxxyyyy_1, g_yzzzz_0_xxxyyyyy_1, g_yzzzz_0_xxxyyyyz_1, g_yzzzz_0_xxxyyyz_1, g_yzzzz_0_xxxyyyzz_1, g_yzzzz_0_xxxyyzz_1, g_yzzzz_0_xxxyyzzz_1, g_yzzzz_0_xxxyzzz_1, g_yzzzz_0_xxxyzzzz_1, g_yzzzz_0_xxyyyyy_1, g_yzzzz_0_xxyyyyyy_1, g_yzzzz_0_xxyyyyyz_1, g_yzzzz_0_xxyyyyz_1, g_yzzzz_0_xxyyyyzz_1, g_yzzzz_0_xxyyyzz_1, g_yzzzz_0_xxyyyzzz_1, g_yzzzz_0_xxyyzzz_1, g_yzzzz_0_xxyyzzzz_1, g_yzzzz_0_xxyzzzz_1, g_yzzzz_0_xxyzzzzz_1, g_yzzzz_0_xyyyyyy_1, g_yzzzz_0_xyyyyyyy_1, g_yzzzz_0_xyyyyyyz_1, g_yzzzz_0_xyyyyyz_1, g_yzzzz_0_xyyyyyzz_1, g_yzzzz_0_xyyyyzz_1, g_yzzzz_0_xyyyyzzz_1, g_yzzzz_0_xyyyzzz_1, g_yzzzz_0_xyyyzzzz_1, g_yzzzz_0_xyyzzzz_1, g_yzzzz_0_xyyzzzzz_1, g_yzzzz_0_xyzzzzz_1, g_yzzzz_0_xyzzzzzz_1, g_yzzzz_0_yyyyyyy_1, g_yzzzz_0_yyyyyyyy_1, g_yzzzz_0_yyyyyyyz_1, g_yzzzz_0_yyyyyyz_1, g_yzzzz_0_yyyyyyzz_1, g_yzzzz_0_yyyyyzz_1, g_yzzzz_0_yyyyyzzz_1, g_yzzzz_0_yyyyzzz_1, g_yzzzz_0_yyyyzzzz_1, g_yzzzz_0_yyyzzzz_1, g_yzzzz_0_yyyzzzzz_1, g_yzzzz_0_yyzzzzz_1, g_yzzzz_0_yyzzzzzz_1, g_yzzzz_0_yzzzzzz_1, g_yzzzz_0_yzzzzzzz_1, g_yzzzz_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xxxxxxxx_0[i] = g_xzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxxxxy_0[i] = 7.0 * g_yzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxxxz_0[i] = g_xzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxxxyy_0[i] = 6.0 * g_yzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxxyz_0[i] = 6.0 * g_yzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxxzz_0[i] = g_xzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxxyyy_0[i] = 5.0 * g_yzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxyyz_0[i] = 5.0 * g_yzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxyzz_0[i] = 5.0 * g_yzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxzzz_0[i] = g_xzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxyyyy_0[i] = 4.0 * g_yzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxyyyz_0[i] = 4.0 * g_yzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxyyzz_0[i] = 4.0 * g_yzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxyzzz_0[i] = 4.0 * g_yzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxzzzz_0[i] = g_xzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxyyyyy_0[i] = 3.0 * g_yzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyyyyz_0[i] = 3.0 * g_yzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyyyzz_0[i] = 3.0 * g_yzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyyzzz_0[i] = 3.0 * g_yzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyzzzz_0[i] = 3.0 * g_yzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxzzzzz_0[i] = g_xzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxyyyyyy_0[i] = 2.0 * g_yzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyyyzz_0[i] = 2.0 * g_yzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyyzzz_0[i] = 2.0 * g_yzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyzzzz_0[i] = 2.0 * g_yzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyzzzzz_0[i] = 2.0 * g_yzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxzzzzzz_0[i] = g_xzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xyyyyyyy_0[i] = g_yzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyyyyz_0[i] = g_yzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyyyzz_0[i] = g_yzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyyzzz_0[i] = g_yzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyzzzz_0[i] = g_yzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyzzzzz_0[i] = g_yzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyzzzzzz_0[i] = g_yzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xzzzzzzz_0[i] = g_xzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_yyyyyyyy_0[i] = g_yzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyyyyz_0[i] = g_yzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyyyzz_0[i] = g_yzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyyzzz_0[i] = g_yzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyzzzz_0[i] = g_yzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyzzzzz_0[i] = g_yzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyzzzzzz_0[i] = g_yzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yzzzzzzz_0[i] = g_yzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_zzzzzzzz_0[i] = g_yzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 900-945 components of targeted buffer : ISL

    auto g_xzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 900);

    auto g_xzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 901);

    auto g_xzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 902);

    auto g_xzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 903);

    auto g_xzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 904);

    auto g_xzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 905);

    auto g_xzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 906);

    auto g_xzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 907);

    auto g_xzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 908);

    auto g_xzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 909);

    auto g_xzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 910);

    auto g_xzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 911);

    auto g_xzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 912);

    auto g_xzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 913);

    auto g_xzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 914);

    auto g_xzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 915);

    auto g_xzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 916);

    auto g_xzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 917);

    auto g_xzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 918);

    auto g_xzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 919);

    auto g_xzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 920);

    auto g_xzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 921);

    auto g_xzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 922);

    auto g_xzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 923);

    auto g_xzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 924);

    auto g_xzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 925);

    auto g_xzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 926);

    auto g_xzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 927);

    auto g_xzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 928);

    auto g_xzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 929);

    auto g_xzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 930);

    auto g_xzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 931);

    auto g_xzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 932);

    auto g_xzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 933);

    auto g_xzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 934);

    auto g_xzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 935);

    auto g_xzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 936);

    auto g_xzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 937);

    auto g_xzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 938);

    auto g_xzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 939);

    auto g_xzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 940);

    auto g_xzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 941);

    auto g_xzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 942);

    auto g_xzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 943);

    auto g_xzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 944);

    #pragma omp simd aligned(g_xzzzzz_0_xxxxxxxx_0, g_xzzzzz_0_xxxxxxxy_0, g_xzzzzz_0_xxxxxxxz_0, g_xzzzzz_0_xxxxxxyy_0, g_xzzzzz_0_xxxxxxyz_0, g_xzzzzz_0_xxxxxxzz_0, g_xzzzzz_0_xxxxxyyy_0, g_xzzzzz_0_xxxxxyyz_0, g_xzzzzz_0_xxxxxyzz_0, g_xzzzzz_0_xxxxxzzz_0, g_xzzzzz_0_xxxxyyyy_0, g_xzzzzz_0_xxxxyyyz_0, g_xzzzzz_0_xxxxyyzz_0, g_xzzzzz_0_xxxxyzzz_0, g_xzzzzz_0_xxxxzzzz_0, g_xzzzzz_0_xxxyyyyy_0, g_xzzzzz_0_xxxyyyyz_0, g_xzzzzz_0_xxxyyyzz_0, g_xzzzzz_0_xxxyyzzz_0, g_xzzzzz_0_xxxyzzzz_0, g_xzzzzz_0_xxxzzzzz_0, g_xzzzzz_0_xxyyyyyy_0, g_xzzzzz_0_xxyyyyyz_0, g_xzzzzz_0_xxyyyyzz_0, g_xzzzzz_0_xxyyyzzz_0, g_xzzzzz_0_xxyyzzzz_0, g_xzzzzz_0_xxyzzzzz_0, g_xzzzzz_0_xxzzzzzz_0, g_xzzzzz_0_xyyyyyyy_0, g_xzzzzz_0_xyyyyyyz_0, g_xzzzzz_0_xyyyyyzz_0, g_xzzzzz_0_xyyyyzzz_0, g_xzzzzz_0_xyyyzzzz_0, g_xzzzzz_0_xyyzzzzz_0, g_xzzzzz_0_xyzzzzzz_0, g_xzzzzz_0_xzzzzzzz_0, g_xzzzzz_0_yyyyyyyy_0, g_xzzzzz_0_yyyyyyyz_0, g_xzzzzz_0_yyyyyyzz_0, g_xzzzzz_0_yyyyyzzz_0, g_xzzzzz_0_yyyyzzzz_0, g_xzzzzz_0_yyyzzzzz_0, g_xzzzzz_0_yyzzzzzz_0, g_xzzzzz_0_yzzzzzzz_0, g_xzzzzz_0_zzzzzzzz_0, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxxx_1, g_zzzzz_0_xxxxxxxy_1, g_zzzzz_0_xxxxxxxz_1, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxyy_1, g_zzzzz_0_xxxxxxyz_1, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxxzz_1, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyyy_1, g_zzzzz_0_xxxxxyyz_1, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxyzz_1, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxxzzz_1, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyyy_1, g_zzzzz_0_xxxxyyyz_1, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyyzz_1, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxyzzz_1, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxxzzzz_1, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyyy_1, g_zzzzz_0_xxxyyyyz_1, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyyzz_1, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyyzzz_1, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxyzzzz_1, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxxzzzzz_1, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyyy_1, g_zzzzz_0_xxyyyyyz_1, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyyzz_1, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyyzzz_1, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyyzzzz_1, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxyzzzzz_1, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xxzzzzzz_1, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyyy_1, g_zzzzz_0_xyyyyyyz_1, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyyzz_1, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyyzzz_1, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyyzzzz_1, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyyzzzzz_1, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xyzzzzzz_1, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_xzzzzzzz_1, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyyy_1, g_zzzzz_0_yyyyyyyz_1, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyyzz_1, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyyzzz_1, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyyzzzz_1, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyyzzzzz_1, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yyzzzzzz_1, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_yzzzzzzz_1, g_zzzzz_0_zzzzzzz_1, g_zzzzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xxxxxxxx_0[i] = 8.0 * g_zzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxxy_0[i] = 7.0 * g_zzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxxz_0[i] = 7.0 * g_zzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxyy_0[i] = 6.0 * g_zzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxyz_0[i] = 6.0 * g_zzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxzz_0[i] = 6.0 * g_zzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxyyy_0[i] = 5.0 * g_zzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxyyz_0[i] = 5.0 * g_zzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxyzz_0[i] = 5.0 * g_zzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxzzz_0[i] = 5.0 * g_zzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyyyz_0[i] = 4.0 * g_zzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyyzz_0[i] = 4.0 * g_zzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyzzz_0[i] = 4.0 * g_zzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxzzzz_0[i] = 4.0 * g_zzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyyyy_0[i] = 3.0 * g_zzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyyyz_0[i] = 3.0 * g_zzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyzzz_0[i] = 3.0 * g_zzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyzzzz_0[i] = 3.0 * g_zzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxzzzzz_0[i] = 3.0 * g_zzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyyyy_0[i] = 2.0 * g_zzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyyyz_0[i] = 2.0 * g_zzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyyzz_0[i] = 2.0 * g_zzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyzzz_0[i] = 2.0 * g_zzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyzzzzz_0[i] = 2.0 * g_zzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxzzzzzz_0[i] = 2.0 * g_zzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyyyy_0[i] = g_zzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyyyz_0[i] = g_zzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyyzz_0[i] = g_zzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyzzz_0[i] = g_zzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyzzzz_0[i] = g_zzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyzzzzz_0[i] = g_zzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyzzzzzz_0[i] = g_zzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xzzzzzzz_0[i] = g_zzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyyyy_0[i] = g_zzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyyyz_0[i] = g_zzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyyzz_0[i] = g_zzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyzzz_0[i] = g_zzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyzzzz_0[i] = g_zzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyzzzzz_0[i] = g_zzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyzzzzzz_0[i] = g_zzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yzzzzzzz_0[i] = g_zzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_zzzzzzzz_0[i] = g_zzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 945-990 components of targeted buffer : ISL

    auto g_yyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 945);

    auto g_yyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 946);

    auto g_yyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 947);

    auto g_yyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 948);

    auto g_yyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 949);

    auto g_yyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 950);

    auto g_yyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 951);

    auto g_yyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 952);

    auto g_yyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 953);

    auto g_yyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 954);

    auto g_yyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 955);

    auto g_yyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 956);

    auto g_yyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 957);

    auto g_yyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 958);

    auto g_yyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 959);

    auto g_yyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 960);

    auto g_yyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 961);

    auto g_yyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 962);

    auto g_yyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 963);

    auto g_yyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 964);

    auto g_yyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 965);

    auto g_yyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 966);

    auto g_yyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 967);

    auto g_yyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 968);

    auto g_yyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 969);

    auto g_yyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 970);

    auto g_yyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 971);

    auto g_yyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 972);

    auto g_yyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 973);

    auto g_yyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 974);

    auto g_yyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 975);

    auto g_yyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 976);

    auto g_yyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 977);

    auto g_yyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 978);

    auto g_yyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 979);

    auto g_yyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 980);

    auto g_yyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 981);

    auto g_yyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 982);

    auto g_yyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 983);

    auto g_yyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 984);

    auto g_yyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 985);

    auto g_yyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 986);

    auto g_yyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 987);

    auto g_yyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 988);

    auto g_yyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 989);

    #pragma omp simd aligned(g_yyyy_0_xxxxxxxx_0, g_yyyy_0_xxxxxxxx_1, g_yyyy_0_xxxxxxxy_0, g_yyyy_0_xxxxxxxy_1, g_yyyy_0_xxxxxxxz_0, g_yyyy_0_xxxxxxxz_1, g_yyyy_0_xxxxxxyy_0, g_yyyy_0_xxxxxxyy_1, g_yyyy_0_xxxxxxyz_0, g_yyyy_0_xxxxxxyz_1, g_yyyy_0_xxxxxxzz_0, g_yyyy_0_xxxxxxzz_1, g_yyyy_0_xxxxxyyy_0, g_yyyy_0_xxxxxyyy_1, g_yyyy_0_xxxxxyyz_0, g_yyyy_0_xxxxxyyz_1, g_yyyy_0_xxxxxyzz_0, g_yyyy_0_xxxxxyzz_1, g_yyyy_0_xxxxxzzz_0, g_yyyy_0_xxxxxzzz_1, g_yyyy_0_xxxxyyyy_0, g_yyyy_0_xxxxyyyy_1, g_yyyy_0_xxxxyyyz_0, g_yyyy_0_xxxxyyyz_1, g_yyyy_0_xxxxyyzz_0, g_yyyy_0_xxxxyyzz_1, g_yyyy_0_xxxxyzzz_0, g_yyyy_0_xxxxyzzz_1, g_yyyy_0_xxxxzzzz_0, g_yyyy_0_xxxxzzzz_1, g_yyyy_0_xxxyyyyy_0, g_yyyy_0_xxxyyyyy_1, g_yyyy_0_xxxyyyyz_0, g_yyyy_0_xxxyyyyz_1, g_yyyy_0_xxxyyyzz_0, g_yyyy_0_xxxyyyzz_1, g_yyyy_0_xxxyyzzz_0, g_yyyy_0_xxxyyzzz_1, g_yyyy_0_xxxyzzzz_0, g_yyyy_0_xxxyzzzz_1, g_yyyy_0_xxxzzzzz_0, g_yyyy_0_xxxzzzzz_1, g_yyyy_0_xxyyyyyy_0, g_yyyy_0_xxyyyyyy_1, g_yyyy_0_xxyyyyyz_0, g_yyyy_0_xxyyyyyz_1, g_yyyy_0_xxyyyyzz_0, g_yyyy_0_xxyyyyzz_1, g_yyyy_0_xxyyyzzz_0, g_yyyy_0_xxyyyzzz_1, g_yyyy_0_xxyyzzzz_0, g_yyyy_0_xxyyzzzz_1, g_yyyy_0_xxyzzzzz_0, g_yyyy_0_xxyzzzzz_1, g_yyyy_0_xxzzzzzz_0, g_yyyy_0_xxzzzzzz_1, g_yyyy_0_xyyyyyyy_0, g_yyyy_0_xyyyyyyy_1, g_yyyy_0_xyyyyyyz_0, g_yyyy_0_xyyyyyyz_1, g_yyyy_0_xyyyyyzz_0, g_yyyy_0_xyyyyyzz_1, g_yyyy_0_xyyyyzzz_0, g_yyyy_0_xyyyyzzz_1, g_yyyy_0_xyyyzzzz_0, g_yyyy_0_xyyyzzzz_1, g_yyyy_0_xyyzzzzz_0, g_yyyy_0_xyyzzzzz_1, g_yyyy_0_xyzzzzzz_0, g_yyyy_0_xyzzzzzz_1, g_yyyy_0_xzzzzzzz_0, g_yyyy_0_xzzzzzzz_1, g_yyyy_0_yyyyyyyy_0, g_yyyy_0_yyyyyyyy_1, g_yyyy_0_yyyyyyyz_0, g_yyyy_0_yyyyyyyz_1, g_yyyy_0_yyyyyyzz_0, g_yyyy_0_yyyyyyzz_1, g_yyyy_0_yyyyyzzz_0, g_yyyy_0_yyyyyzzz_1, g_yyyy_0_yyyyzzzz_0, g_yyyy_0_yyyyzzzz_1, g_yyyy_0_yyyzzzzz_0, g_yyyy_0_yyyzzzzz_1, g_yyyy_0_yyzzzzzz_0, g_yyyy_0_yyzzzzzz_1, g_yyyy_0_yzzzzzzz_0, g_yyyy_0_yzzzzzzz_1, g_yyyy_0_zzzzzzzz_0, g_yyyy_0_zzzzzzzz_1, g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxxx_1, g_yyyyy_0_xxxxxxxy_1, g_yyyyy_0_xxxxxxxz_1, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxyy_1, g_yyyyy_0_xxxxxxyz_1, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxxzz_1, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyyy_1, g_yyyyy_0_xxxxxyyz_1, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxyzz_1, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxxzzz_1, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyyy_1, g_yyyyy_0_xxxxyyyz_1, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyyzz_1, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxyzzz_1, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxxzzzz_1, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyyy_1, g_yyyyy_0_xxxyyyyz_1, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyyzz_1, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyyzzz_1, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxyzzzz_1, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxxzzzzz_1, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyyy_1, g_yyyyy_0_xxyyyyyz_1, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyyzz_1, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyyzzz_1, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyyzzzz_1, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxyzzzzz_1, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xxzzzzzz_1, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyyy_1, g_yyyyy_0_xyyyyyyz_1, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyyzz_1, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyyzzz_1, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyyzzzz_1, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyyzzzzz_1, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xyzzzzzz_1, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_xzzzzzzz_1, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyyy_1, g_yyyyy_0_yyyyyyyz_1, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyyzz_1, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyyzzz_1, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyyzzzz_1, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyyzzzzz_1, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yyzzzzzz_1, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_yzzzzzzz_1, g_yyyyy_0_zzzzzzz_1, g_yyyyy_0_zzzzzzzz_1, g_yyyyyy_0_xxxxxxxx_0, g_yyyyyy_0_xxxxxxxy_0, g_yyyyyy_0_xxxxxxxz_0, g_yyyyyy_0_xxxxxxyy_0, g_yyyyyy_0_xxxxxxyz_0, g_yyyyyy_0_xxxxxxzz_0, g_yyyyyy_0_xxxxxyyy_0, g_yyyyyy_0_xxxxxyyz_0, g_yyyyyy_0_xxxxxyzz_0, g_yyyyyy_0_xxxxxzzz_0, g_yyyyyy_0_xxxxyyyy_0, g_yyyyyy_0_xxxxyyyz_0, g_yyyyyy_0_xxxxyyzz_0, g_yyyyyy_0_xxxxyzzz_0, g_yyyyyy_0_xxxxzzzz_0, g_yyyyyy_0_xxxyyyyy_0, g_yyyyyy_0_xxxyyyyz_0, g_yyyyyy_0_xxxyyyzz_0, g_yyyyyy_0_xxxyyzzz_0, g_yyyyyy_0_xxxyzzzz_0, g_yyyyyy_0_xxxzzzzz_0, g_yyyyyy_0_xxyyyyyy_0, g_yyyyyy_0_xxyyyyyz_0, g_yyyyyy_0_xxyyyyzz_0, g_yyyyyy_0_xxyyyzzz_0, g_yyyyyy_0_xxyyzzzz_0, g_yyyyyy_0_xxyzzzzz_0, g_yyyyyy_0_xxzzzzzz_0, g_yyyyyy_0_xyyyyyyy_0, g_yyyyyy_0_xyyyyyyz_0, g_yyyyyy_0_xyyyyyzz_0, g_yyyyyy_0_xyyyyzzz_0, g_yyyyyy_0_xyyyzzzz_0, g_yyyyyy_0_xyyzzzzz_0, g_yyyyyy_0_xyzzzzzz_0, g_yyyyyy_0_xzzzzzzz_0, g_yyyyyy_0_yyyyyyyy_0, g_yyyyyy_0_yyyyyyyz_0, g_yyyyyy_0_yyyyyyzz_0, g_yyyyyy_0_yyyyyzzz_0, g_yyyyyy_0_yyyyzzzz_0, g_yyyyyy_0_yyyzzzzz_0, g_yyyyyy_0_yyzzzzzz_0, g_yyyyyy_0_yzzzzzzz_0, g_yyyyyy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xxxxxxxx_0[i] = 5.0 * g_yyyy_0_xxxxxxxx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxxx_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxxy_0[i] = 5.0 * g_yyyy_0_xxxxxxxy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxxy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxxz_0[i] = 5.0 * g_yyyy_0_xxxxxxxz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxxz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxyy_0[i] = 5.0 * g_yyyy_0_xxxxxxyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxyz_0[i] = 5.0 * g_yyyy_0_xxxxxxyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxyz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxzz_0[i] = 5.0 * g_yyyy_0_xxxxxxzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxyyy_0[i] = 5.0 * g_yyyy_0_xxxxxyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxyyz_0[i] = 5.0 * g_yyyy_0_xxxxxyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxyzz_0[i] = 5.0 * g_yyyy_0_xxxxxyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxyzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxzzz_0[i] = 5.0 * g_yyyy_0_xxxxxzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyyyy_0[i] = 5.0 * g_yyyy_0_xxxxyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyyyz_0[i] = 5.0 * g_yyyy_0_xxxxyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyyzz_0[i] = 5.0 * g_yyyy_0_xxxxyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyzzz_0[i] = 5.0 * g_yyyy_0_xxxxyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxzzzz_0[i] = 5.0 * g_yyyy_0_xxxxzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyyyy_0[i] = 5.0 * g_yyyy_0_xxxyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyyyz_0[i] = 5.0 * g_yyyy_0_xxxyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyyzz_0[i] = 5.0 * g_yyyy_0_xxxyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyzzz_0[i] = 5.0 * g_yyyy_0_xxxyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyzzzz_0[i] = 5.0 * g_yyyy_0_xxxyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxzzzzz_0[i] = 5.0 * g_yyyy_0_xxxzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyyyy_0[i] = 5.0 * g_yyyy_0_xxyyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyyyz_0[i] = 5.0 * g_yyyy_0_xxyyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyyzz_0[i] = 5.0 * g_yyyy_0_xxyyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyzzz_0[i] = 5.0 * g_yyyy_0_xxyyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyzzzz_0[i] = 5.0 * g_yyyy_0_xxyyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyzzzzz_0[i] = 5.0 * g_yyyy_0_xxyzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxzzzzzz_0[i] = 5.0 * g_yyyy_0_xxzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyyyy_0[i] = 5.0 * g_yyyy_0_xyyyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyyyz_0[i] = 5.0 * g_yyyy_0_xyyyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyyzz_0[i] = 5.0 * g_yyyy_0_xyyyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyzzz_0[i] = 5.0 * g_yyyy_0_xyyyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyzzzz_0[i] = 5.0 * g_yyyy_0_xyyyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyzzzzz_0[i] = 5.0 * g_yyyy_0_xyyzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyzzzzzz_0[i] = 5.0 * g_yyyy_0_xyzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xzzzzzzz_0[i] = 5.0 * g_yyyy_0_xzzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xzzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyyyy_0[i] = 5.0 * g_yyyy_0_yyyyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyyyz_0[i] = 5.0 * g_yyyy_0_yyyyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyyzz_0[i] = 5.0 * g_yyyy_0_yyyyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyzzz_0[i] = 5.0 * g_yyyy_0_yyyyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyzzzz_0[i] = 5.0 * g_yyyy_0_yyyyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyzzzzz_0[i] = 5.0 * g_yyyy_0_yyyzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyzzzzzz_0[i] = 5.0 * g_yyyy_0_yyzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yzzzzzzz_0[i] = 5.0 * g_yyyy_0_yzzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_zzzzzzzz_0[i] = 5.0 * g_yyyy_0_zzzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 990-1035 components of targeted buffer : ISL

    auto g_yyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 990);

    auto g_yyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 991);

    auto g_yyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 992);

    auto g_yyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 993);

    auto g_yyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 994);

    auto g_yyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 995);

    auto g_yyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 996);

    auto g_yyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 997);

    auto g_yyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 998);

    auto g_yyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 999);

    auto g_yyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 1000);

    auto g_yyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 1001);

    auto g_yyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 1002);

    auto g_yyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 1003);

    auto g_yyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 1004);

    auto g_yyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1005);

    auto g_yyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1006);

    auto g_yyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1007);

    auto g_yyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1008);

    auto g_yyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1009);

    auto g_yyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1010);

    auto g_yyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1011);

    auto g_yyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1012);

    auto g_yyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1013);

    auto g_yyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1014);

    auto g_yyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1015);

    auto g_yyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1016);

    auto g_yyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1017);

    auto g_yyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1018);

    auto g_yyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1019);

    auto g_yyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1020);

    auto g_yyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1021);

    auto g_yyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1022);

    auto g_yyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1023);

    auto g_yyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1024);

    auto g_yyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1025);

    auto g_yyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1026);

    auto g_yyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1027);

    auto g_yyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1028);

    auto g_yyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1029);

    auto g_yyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1030);

    auto g_yyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1031);

    auto g_yyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1032);

    auto g_yyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1033);

    auto g_yyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1034);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxxx_1, g_yyyyy_0_xxxxxxxy_1, g_yyyyy_0_xxxxxxxz_1, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxyy_1, g_yyyyy_0_xxxxxxyz_1, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxxzz_1, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyyy_1, g_yyyyy_0_xxxxxyyz_1, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxyzz_1, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxxzzz_1, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyyy_1, g_yyyyy_0_xxxxyyyz_1, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyyzz_1, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxyzzz_1, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxxzzzz_1, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyyy_1, g_yyyyy_0_xxxyyyyz_1, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyyzz_1, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyyzzz_1, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxyzzzz_1, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxxzzzzz_1, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyyy_1, g_yyyyy_0_xxyyyyyz_1, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyyzz_1, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyyzzz_1, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyyzzzz_1, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxyzzzzz_1, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xxzzzzzz_1, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyyy_1, g_yyyyy_0_xyyyyyyz_1, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyyzz_1, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyyzzz_1, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyyzzzz_1, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyyzzzzz_1, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xyzzzzzz_1, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_xzzzzzzz_1, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyyy_1, g_yyyyy_0_yyyyyyyz_1, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyyzz_1, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyyzzz_1, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyyzzzz_1, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyyzzzzz_1, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yyzzzzzz_1, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_yzzzzzzz_1, g_yyyyy_0_zzzzzzz_1, g_yyyyy_0_zzzzzzzz_1, g_yyyyyz_0_xxxxxxxx_0, g_yyyyyz_0_xxxxxxxy_0, g_yyyyyz_0_xxxxxxxz_0, g_yyyyyz_0_xxxxxxyy_0, g_yyyyyz_0_xxxxxxyz_0, g_yyyyyz_0_xxxxxxzz_0, g_yyyyyz_0_xxxxxyyy_0, g_yyyyyz_0_xxxxxyyz_0, g_yyyyyz_0_xxxxxyzz_0, g_yyyyyz_0_xxxxxzzz_0, g_yyyyyz_0_xxxxyyyy_0, g_yyyyyz_0_xxxxyyyz_0, g_yyyyyz_0_xxxxyyzz_0, g_yyyyyz_0_xxxxyzzz_0, g_yyyyyz_0_xxxxzzzz_0, g_yyyyyz_0_xxxyyyyy_0, g_yyyyyz_0_xxxyyyyz_0, g_yyyyyz_0_xxxyyyzz_0, g_yyyyyz_0_xxxyyzzz_0, g_yyyyyz_0_xxxyzzzz_0, g_yyyyyz_0_xxxzzzzz_0, g_yyyyyz_0_xxyyyyyy_0, g_yyyyyz_0_xxyyyyyz_0, g_yyyyyz_0_xxyyyyzz_0, g_yyyyyz_0_xxyyyzzz_0, g_yyyyyz_0_xxyyzzzz_0, g_yyyyyz_0_xxyzzzzz_0, g_yyyyyz_0_xxzzzzzz_0, g_yyyyyz_0_xyyyyyyy_0, g_yyyyyz_0_xyyyyyyz_0, g_yyyyyz_0_xyyyyyzz_0, g_yyyyyz_0_xyyyyzzz_0, g_yyyyyz_0_xyyyzzzz_0, g_yyyyyz_0_xyyzzzzz_0, g_yyyyyz_0_xyzzzzzz_0, g_yyyyyz_0_xzzzzzzz_0, g_yyyyyz_0_yyyyyyyy_0, g_yyyyyz_0_yyyyyyyz_0, g_yyyyyz_0_yyyyyyzz_0, g_yyyyyz_0_yyyyyzzz_0, g_yyyyyz_0_yyyyzzzz_0, g_yyyyyz_0_yyyzzzzz_0, g_yyyyyz_0_yyzzzzzz_0, g_yyyyyz_0_yzzzzzzz_0, g_yyyyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xxxxxxxx_0[i] = g_yyyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxxy_0[i] = g_yyyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxxz_0[i] = g_yyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxyy_0[i] = g_yyyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxyz_0[i] = g_yyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxzz_0[i] = 2.0 * g_yyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxyyy_0[i] = g_yyyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxyyz_0[i] = g_yyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxyzz_0[i] = 2.0 * g_yyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxzzz_0[i] = 3.0 * g_yyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyyyy_0[i] = g_yyyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyyyz_0[i] = g_yyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyyzz_0[i] = 2.0 * g_yyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyzzz_0[i] = 3.0 * g_yyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyyyy_0[i] = g_yyyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyyyz_0[i] = g_yyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyyzz_0[i] = 2.0 * g_yyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyzzzz_0[i] = 4.0 * g_yyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxzzzzz_0[i] = 5.0 * g_yyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyyyy_0[i] = g_yyyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyyyz_0[i] = g_yyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyzzz_0[i] = 3.0 * g_yyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyzzzz_0[i] = 4.0 * g_yyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyzzzzz_0[i] = 5.0 * g_yyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxzzzzzz_0[i] = 6.0 * g_yyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyyyy_0[i] = g_yyyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyyyz_0[i] = g_yyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyyzz_0[i] = 2.0 * g_yyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyzzz_0[i] = 3.0 * g_yyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyzzzz_0[i] = 4.0 * g_yyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyzzzzz_0[i] = 5.0 * g_yyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyzzzzzz_0[i] = 6.0 * g_yyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xzzzzzzz_0[i] = 7.0 * g_yyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyyyy_0[i] = g_yyyyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyyyz_0[i] = g_yyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyyzz_0[i] = 2.0 * g_yyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyzzz_0[i] = 3.0 * g_yyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyzzzz_0[i] = 4.0 * g_yyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyzzzzz_0[i] = 5.0 * g_yyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyzzzzzz_0[i] = 6.0 * g_yyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yzzzzzzz_0[i] = 7.0 * g_yyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_zzzzzzzz_0[i] = 8.0 * g_yyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 1035-1080 components of targeted buffer : ISL

    auto g_yyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 1035);

    auto g_yyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 1036);

    auto g_yyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 1037);

    auto g_yyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 1038);

    auto g_yyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 1039);

    auto g_yyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 1040);

    auto g_yyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 1041);

    auto g_yyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 1042);

    auto g_yyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 1043);

    auto g_yyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 1044);

    auto g_yyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 1045);

    auto g_yyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 1046);

    auto g_yyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 1047);

    auto g_yyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 1048);

    auto g_yyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 1049);

    auto g_yyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1050);

    auto g_yyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1051);

    auto g_yyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1052);

    auto g_yyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1053);

    auto g_yyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1054);

    auto g_yyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1055);

    auto g_yyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1056);

    auto g_yyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1057);

    auto g_yyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1058);

    auto g_yyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1059);

    auto g_yyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1060);

    auto g_yyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1061);

    auto g_yyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1062);

    auto g_yyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1063);

    auto g_yyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1064);

    auto g_yyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1065);

    auto g_yyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1066);

    auto g_yyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1067);

    auto g_yyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1068);

    auto g_yyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1069);

    auto g_yyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1070);

    auto g_yyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1071);

    auto g_yyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1072);

    auto g_yyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1073);

    auto g_yyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1074);

    auto g_yyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1075);

    auto g_yyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1076);

    auto g_yyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1077);

    auto g_yyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1078);

    auto g_yyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1079);

    #pragma omp simd aligned(g_yyyy_0_xxxxxxxy_0, g_yyyy_0_xxxxxxxy_1, g_yyyy_0_xxxxxxyy_0, g_yyyy_0_xxxxxxyy_1, g_yyyy_0_xxxxxyyy_0, g_yyyy_0_xxxxxyyy_1, g_yyyy_0_xxxxyyyy_0, g_yyyy_0_xxxxyyyy_1, g_yyyy_0_xxxyyyyy_0, g_yyyy_0_xxxyyyyy_1, g_yyyy_0_xxyyyyyy_0, g_yyyy_0_xxyyyyyy_1, g_yyyy_0_xyyyyyyy_0, g_yyyy_0_xyyyyyyy_1, g_yyyy_0_yyyyyyyy_0, g_yyyy_0_yyyyyyyy_1, g_yyyyz_0_xxxxxxxy_1, g_yyyyz_0_xxxxxxyy_1, g_yyyyz_0_xxxxxyyy_1, g_yyyyz_0_xxxxyyyy_1, g_yyyyz_0_xxxyyyyy_1, g_yyyyz_0_xxyyyyyy_1, g_yyyyz_0_xyyyyyyy_1, g_yyyyz_0_yyyyyyyy_1, g_yyyyzz_0_xxxxxxxx_0, g_yyyyzz_0_xxxxxxxy_0, g_yyyyzz_0_xxxxxxxz_0, g_yyyyzz_0_xxxxxxyy_0, g_yyyyzz_0_xxxxxxyz_0, g_yyyyzz_0_xxxxxxzz_0, g_yyyyzz_0_xxxxxyyy_0, g_yyyyzz_0_xxxxxyyz_0, g_yyyyzz_0_xxxxxyzz_0, g_yyyyzz_0_xxxxxzzz_0, g_yyyyzz_0_xxxxyyyy_0, g_yyyyzz_0_xxxxyyyz_0, g_yyyyzz_0_xxxxyyzz_0, g_yyyyzz_0_xxxxyzzz_0, g_yyyyzz_0_xxxxzzzz_0, g_yyyyzz_0_xxxyyyyy_0, g_yyyyzz_0_xxxyyyyz_0, g_yyyyzz_0_xxxyyyzz_0, g_yyyyzz_0_xxxyyzzz_0, g_yyyyzz_0_xxxyzzzz_0, g_yyyyzz_0_xxxzzzzz_0, g_yyyyzz_0_xxyyyyyy_0, g_yyyyzz_0_xxyyyyyz_0, g_yyyyzz_0_xxyyyyzz_0, g_yyyyzz_0_xxyyyzzz_0, g_yyyyzz_0_xxyyzzzz_0, g_yyyyzz_0_xxyzzzzz_0, g_yyyyzz_0_xxzzzzzz_0, g_yyyyzz_0_xyyyyyyy_0, g_yyyyzz_0_xyyyyyyz_0, g_yyyyzz_0_xyyyyyzz_0, g_yyyyzz_0_xyyyyzzz_0, g_yyyyzz_0_xyyyzzzz_0, g_yyyyzz_0_xyyzzzzz_0, g_yyyyzz_0_xyzzzzzz_0, g_yyyyzz_0_xzzzzzzz_0, g_yyyyzz_0_yyyyyyyy_0, g_yyyyzz_0_yyyyyyyz_0, g_yyyyzz_0_yyyyyyzz_0, g_yyyyzz_0_yyyyyzzz_0, g_yyyyzz_0_yyyyzzzz_0, g_yyyyzz_0_yyyzzzzz_0, g_yyyyzz_0_yyzzzzzz_0, g_yyyyzz_0_yzzzzzzz_0, g_yyyyzz_0_zzzzzzzz_0, g_yyyzz_0_xxxxxxxx_1, g_yyyzz_0_xxxxxxxz_1, g_yyyzz_0_xxxxxxyz_1, g_yyyzz_0_xxxxxxz_1, g_yyyzz_0_xxxxxxzz_1, g_yyyzz_0_xxxxxyyz_1, g_yyyzz_0_xxxxxyz_1, g_yyyzz_0_xxxxxyzz_1, g_yyyzz_0_xxxxxzz_1, g_yyyzz_0_xxxxxzzz_1, g_yyyzz_0_xxxxyyyz_1, g_yyyzz_0_xxxxyyz_1, g_yyyzz_0_xxxxyyzz_1, g_yyyzz_0_xxxxyzz_1, g_yyyzz_0_xxxxyzzz_1, g_yyyzz_0_xxxxzzz_1, g_yyyzz_0_xxxxzzzz_1, g_yyyzz_0_xxxyyyyz_1, g_yyyzz_0_xxxyyyz_1, g_yyyzz_0_xxxyyyzz_1, g_yyyzz_0_xxxyyzz_1, g_yyyzz_0_xxxyyzzz_1, g_yyyzz_0_xxxyzzz_1, g_yyyzz_0_xxxyzzzz_1, g_yyyzz_0_xxxzzzz_1, g_yyyzz_0_xxxzzzzz_1, g_yyyzz_0_xxyyyyyz_1, g_yyyzz_0_xxyyyyz_1, g_yyyzz_0_xxyyyyzz_1, g_yyyzz_0_xxyyyzz_1, g_yyyzz_0_xxyyyzzz_1, g_yyyzz_0_xxyyzzz_1, g_yyyzz_0_xxyyzzzz_1, g_yyyzz_0_xxyzzzz_1, g_yyyzz_0_xxyzzzzz_1, g_yyyzz_0_xxzzzzz_1, g_yyyzz_0_xxzzzzzz_1, g_yyyzz_0_xyyyyyyz_1, g_yyyzz_0_xyyyyyz_1, g_yyyzz_0_xyyyyyzz_1, g_yyyzz_0_xyyyyzz_1, g_yyyzz_0_xyyyyzzz_1, g_yyyzz_0_xyyyzzz_1, g_yyyzz_0_xyyyzzzz_1, g_yyyzz_0_xyyzzzz_1, g_yyyzz_0_xyyzzzzz_1, g_yyyzz_0_xyzzzzz_1, g_yyyzz_0_xyzzzzzz_1, g_yyyzz_0_xzzzzzz_1, g_yyyzz_0_xzzzzzzz_1, g_yyyzz_0_yyyyyyyz_1, g_yyyzz_0_yyyyyyz_1, g_yyyzz_0_yyyyyyzz_1, g_yyyzz_0_yyyyyzz_1, g_yyyzz_0_yyyyyzzz_1, g_yyyzz_0_yyyyzzz_1, g_yyyzz_0_yyyyzzzz_1, g_yyyzz_0_yyyzzzz_1, g_yyyzz_0_yyyzzzzz_1, g_yyyzz_0_yyzzzzz_1, g_yyyzz_0_yyzzzzzz_1, g_yyyzz_0_yzzzzzz_1, g_yyyzz_0_yzzzzzzz_1, g_yyyzz_0_zzzzzzz_1, g_yyyzz_0_zzzzzzzz_1, g_yyzz_0_xxxxxxxx_0, g_yyzz_0_xxxxxxxx_1, g_yyzz_0_xxxxxxxz_0, g_yyzz_0_xxxxxxxz_1, g_yyzz_0_xxxxxxyz_0, g_yyzz_0_xxxxxxyz_1, g_yyzz_0_xxxxxxzz_0, g_yyzz_0_xxxxxxzz_1, g_yyzz_0_xxxxxyyz_0, g_yyzz_0_xxxxxyyz_1, g_yyzz_0_xxxxxyzz_0, g_yyzz_0_xxxxxyzz_1, g_yyzz_0_xxxxxzzz_0, g_yyzz_0_xxxxxzzz_1, g_yyzz_0_xxxxyyyz_0, g_yyzz_0_xxxxyyyz_1, g_yyzz_0_xxxxyyzz_0, g_yyzz_0_xxxxyyzz_1, g_yyzz_0_xxxxyzzz_0, g_yyzz_0_xxxxyzzz_1, g_yyzz_0_xxxxzzzz_0, g_yyzz_0_xxxxzzzz_1, g_yyzz_0_xxxyyyyz_0, g_yyzz_0_xxxyyyyz_1, g_yyzz_0_xxxyyyzz_0, g_yyzz_0_xxxyyyzz_1, g_yyzz_0_xxxyyzzz_0, g_yyzz_0_xxxyyzzz_1, g_yyzz_0_xxxyzzzz_0, g_yyzz_0_xxxyzzzz_1, g_yyzz_0_xxxzzzzz_0, g_yyzz_0_xxxzzzzz_1, g_yyzz_0_xxyyyyyz_0, g_yyzz_0_xxyyyyyz_1, g_yyzz_0_xxyyyyzz_0, g_yyzz_0_xxyyyyzz_1, g_yyzz_0_xxyyyzzz_0, g_yyzz_0_xxyyyzzz_1, g_yyzz_0_xxyyzzzz_0, g_yyzz_0_xxyyzzzz_1, g_yyzz_0_xxyzzzzz_0, g_yyzz_0_xxyzzzzz_1, g_yyzz_0_xxzzzzzz_0, g_yyzz_0_xxzzzzzz_1, g_yyzz_0_xyyyyyyz_0, g_yyzz_0_xyyyyyyz_1, g_yyzz_0_xyyyyyzz_0, g_yyzz_0_xyyyyyzz_1, g_yyzz_0_xyyyyzzz_0, g_yyzz_0_xyyyyzzz_1, g_yyzz_0_xyyyzzzz_0, g_yyzz_0_xyyyzzzz_1, g_yyzz_0_xyyzzzzz_0, g_yyzz_0_xyyzzzzz_1, g_yyzz_0_xyzzzzzz_0, g_yyzz_0_xyzzzzzz_1, g_yyzz_0_xzzzzzzz_0, g_yyzz_0_xzzzzzzz_1, g_yyzz_0_yyyyyyyz_0, g_yyzz_0_yyyyyyyz_1, g_yyzz_0_yyyyyyzz_0, g_yyzz_0_yyyyyyzz_1, g_yyzz_0_yyyyyzzz_0, g_yyzz_0_yyyyyzzz_1, g_yyzz_0_yyyyzzzz_0, g_yyzz_0_yyyyzzzz_1, g_yyzz_0_yyyzzzzz_0, g_yyzz_0_yyyzzzzz_1, g_yyzz_0_yyzzzzzz_0, g_yyzz_0_yyzzzzzz_1, g_yyzz_0_yzzzzzzz_0, g_yyzz_0_yzzzzzzz_1, g_yyzz_0_zzzzzzzz_0, g_yyzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xxxxxxxx_0[i] = 3.0 * g_yyzz_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxxxy_0[i] = g_yyyy_0_xxxxxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxxxxz_0[i] = 3.0 * g_yyzz_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxxyy_0[i] = g_yyyy_0_xxxxxxyy_0[i] * fbe_0 - g_yyyy_0_xxxxxxyy_1[i] * fz_be_0 + g_yyyyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxxxyz_0[i] = 3.0 * g_yyzz_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxxzz_0[i] = 3.0 * g_yyzz_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxyyy_0[i] = g_yyyy_0_xxxxxyyy_0[i] * fbe_0 - g_yyyy_0_xxxxxyyy_1[i] * fz_be_0 + g_yyyyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxxyyz_0[i] = 3.0 * g_yyzz_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxyzz_0[i] = 3.0 * g_yyzz_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxzzz_0[i] = 3.0 * g_yyzz_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxyyyy_0[i] = g_yyyy_0_xxxxyyyy_0[i] * fbe_0 - g_yyyy_0_xxxxyyyy_1[i] * fz_be_0 + g_yyyyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxyyyz_0[i] = 3.0 * g_yyzz_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxyyzz_0[i] = 3.0 * g_yyzz_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxyzzz_0[i] = 3.0 * g_yyzz_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxzzzz_0[i] = 3.0 * g_yyzz_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyyyyy_0[i] = g_yyyy_0_xxxyyyyy_0[i] * fbe_0 - g_yyyy_0_xxxyyyyy_1[i] * fz_be_0 + g_yyyyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxyyyyz_0[i] = 3.0 * g_yyzz_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyyyzz_0[i] = 3.0 * g_yyzz_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyyzzz_0[i] = 3.0 * g_yyzz_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyzzzz_0[i] = 3.0 * g_yyzz_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxzzzzz_0[i] = 3.0 * g_yyzz_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyyyyy_0[i] = g_yyyy_0_xxyyyyyy_0[i] * fbe_0 - g_yyyy_0_xxyyyyyy_1[i] * fz_be_0 + g_yyyyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxyyyyyz_0[i] = 3.0 * g_yyzz_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyyyzz_0[i] = 3.0 * g_yyzz_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyyzzz_0[i] = 3.0 * g_yyzz_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyzzzz_0[i] = 3.0 * g_yyzz_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyzzzzz_0[i] = 3.0 * g_yyzz_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxzzzzzz_0[i] = 3.0 * g_yyzz_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyyyyy_0[i] = g_yyyy_0_xyyyyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_yyyyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xyyyyyyz_0[i] = 3.0 * g_yyzz_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyyyzz_0[i] = 3.0 * g_yyzz_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyyzzz_0[i] = 3.0 * g_yyzz_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyzzzz_0[i] = 3.0 * g_yyzz_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyzzzzz_0[i] = 3.0 * g_yyzz_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyzzzzzz_0[i] = 3.0 * g_yyzz_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xzzzzzzz_0[i] = 3.0 * g_yyzz_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyyyyy_0[i] = g_yyyy_0_yyyyyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_yyyyz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_yyyyyyyz_0[i] = 3.0 * g_yyzz_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyyyzz_0[i] = 3.0 * g_yyzz_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyyzzz_0[i] = 3.0 * g_yyzz_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyzzzz_0[i] = 3.0 * g_yyzz_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyzzzzz_0[i] = 3.0 * g_yyzz_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyzzzzzz_0[i] = 3.0 * g_yyzz_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yzzzzzzz_0[i] = 3.0 * g_yyzz_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_zzzzzzzz_0[i] = 3.0 * g_yyzz_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1080-1125 components of targeted buffer : ISL

    auto g_yyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 1080);

    auto g_yyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 1081);

    auto g_yyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 1082);

    auto g_yyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 1083);

    auto g_yyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 1084);

    auto g_yyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 1085);

    auto g_yyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 1086);

    auto g_yyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 1087);

    auto g_yyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 1088);

    auto g_yyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 1089);

    auto g_yyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 1090);

    auto g_yyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 1091);

    auto g_yyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 1092);

    auto g_yyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 1093);

    auto g_yyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 1094);

    auto g_yyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1095);

    auto g_yyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1096);

    auto g_yyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1097);

    auto g_yyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1098);

    auto g_yyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1099);

    auto g_yyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1100);

    auto g_yyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1101);

    auto g_yyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1102);

    auto g_yyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1103);

    auto g_yyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1104);

    auto g_yyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1105);

    auto g_yyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1106);

    auto g_yyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1107);

    auto g_yyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1108);

    auto g_yyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1109);

    auto g_yyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1110);

    auto g_yyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1111);

    auto g_yyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1112);

    auto g_yyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1113);

    auto g_yyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1114);

    auto g_yyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1115);

    auto g_yyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1116);

    auto g_yyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1117);

    auto g_yyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1118);

    auto g_yyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1119);

    auto g_yyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1120);

    auto g_yyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1121);

    auto g_yyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1122);

    auto g_yyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1123);

    auto g_yyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1124);

    #pragma omp simd aligned(g_yyyz_0_xxxxxxxy_0, g_yyyz_0_xxxxxxxy_1, g_yyyz_0_xxxxxxyy_0, g_yyyz_0_xxxxxxyy_1, g_yyyz_0_xxxxxyyy_0, g_yyyz_0_xxxxxyyy_1, g_yyyz_0_xxxxyyyy_0, g_yyyz_0_xxxxyyyy_1, g_yyyz_0_xxxyyyyy_0, g_yyyz_0_xxxyyyyy_1, g_yyyz_0_xxyyyyyy_0, g_yyyz_0_xxyyyyyy_1, g_yyyz_0_xyyyyyyy_0, g_yyyz_0_xyyyyyyy_1, g_yyyz_0_yyyyyyyy_0, g_yyyz_0_yyyyyyyy_1, g_yyyzz_0_xxxxxxxy_1, g_yyyzz_0_xxxxxxyy_1, g_yyyzz_0_xxxxxyyy_1, g_yyyzz_0_xxxxyyyy_1, g_yyyzz_0_xxxyyyyy_1, g_yyyzz_0_xxyyyyyy_1, g_yyyzz_0_xyyyyyyy_1, g_yyyzz_0_yyyyyyyy_1, g_yyyzzz_0_xxxxxxxx_0, g_yyyzzz_0_xxxxxxxy_0, g_yyyzzz_0_xxxxxxxz_0, g_yyyzzz_0_xxxxxxyy_0, g_yyyzzz_0_xxxxxxyz_0, g_yyyzzz_0_xxxxxxzz_0, g_yyyzzz_0_xxxxxyyy_0, g_yyyzzz_0_xxxxxyyz_0, g_yyyzzz_0_xxxxxyzz_0, g_yyyzzz_0_xxxxxzzz_0, g_yyyzzz_0_xxxxyyyy_0, g_yyyzzz_0_xxxxyyyz_0, g_yyyzzz_0_xxxxyyzz_0, g_yyyzzz_0_xxxxyzzz_0, g_yyyzzz_0_xxxxzzzz_0, g_yyyzzz_0_xxxyyyyy_0, g_yyyzzz_0_xxxyyyyz_0, g_yyyzzz_0_xxxyyyzz_0, g_yyyzzz_0_xxxyyzzz_0, g_yyyzzz_0_xxxyzzzz_0, g_yyyzzz_0_xxxzzzzz_0, g_yyyzzz_0_xxyyyyyy_0, g_yyyzzz_0_xxyyyyyz_0, g_yyyzzz_0_xxyyyyzz_0, g_yyyzzz_0_xxyyyzzz_0, g_yyyzzz_0_xxyyzzzz_0, g_yyyzzz_0_xxyzzzzz_0, g_yyyzzz_0_xxzzzzzz_0, g_yyyzzz_0_xyyyyyyy_0, g_yyyzzz_0_xyyyyyyz_0, g_yyyzzz_0_xyyyyyzz_0, g_yyyzzz_0_xyyyyzzz_0, g_yyyzzz_0_xyyyzzzz_0, g_yyyzzz_0_xyyzzzzz_0, g_yyyzzz_0_xyzzzzzz_0, g_yyyzzz_0_xzzzzzzz_0, g_yyyzzz_0_yyyyyyyy_0, g_yyyzzz_0_yyyyyyyz_0, g_yyyzzz_0_yyyyyyzz_0, g_yyyzzz_0_yyyyyzzz_0, g_yyyzzz_0_yyyyzzzz_0, g_yyyzzz_0_yyyzzzzz_0, g_yyyzzz_0_yyzzzzzz_0, g_yyyzzz_0_yzzzzzzz_0, g_yyyzzz_0_zzzzzzzz_0, g_yyzzz_0_xxxxxxxx_1, g_yyzzz_0_xxxxxxxz_1, g_yyzzz_0_xxxxxxyz_1, g_yyzzz_0_xxxxxxz_1, g_yyzzz_0_xxxxxxzz_1, g_yyzzz_0_xxxxxyyz_1, g_yyzzz_0_xxxxxyz_1, g_yyzzz_0_xxxxxyzz_1, g_yyzzz_0_xxxxxzz_1, g_yyzzz_0_xxxxxzzz_1, g_yyzzz_0_xxxxyyyz_1, g_yyzzz_0_xxxxyyz_1, g_yyzzz_0_xxxxyyzz_1, g_yyzzz_0_xxxxyzz_1, g_yyzzz_0_xxxxyzzz_1, g_yyzzz_0_xxxxzzz_1, g_yyzzz_0_xxxxzzzz_1, g_yyzzz_0_xxxyyyyz_1, g_yyzzz_0_xxxyyyz_1, g_yyzzz_0_xxxyyyzz_1, g_yyzzz_0_xxxyyzz_1, g_yyzzz_0_xxxyyzzz_1, g_yyzzz_0_xxxyzzz_1, g_yyzzz_0_xxxyzzzz_1, g_yyzzz_0_xxxzzzz_1, g_yyzzz_0_xxxzzzzz_1, g_yyzzz_0_xxyyyyyz_1, g_yyzzz_0_xxyyyyz_1, g_yyzzz_0_xxyyyyzz_1, g_yyzzz_0_xxyyyzz_1, g_yyzzz_0_xxyyyzzz_1, g_yyzzz_0_xxyyzzz_1, g_yyzzz_0_xxyyzzzz_1, g_yyzzz_0_xxyzzzz_1, g_yyzzz_0_xxyzzzzz_1, g_yyzzz_0_xxzzzzz_1, g_yyzzz_0_xxzzzzzz_1, g_yyzzz_0_xyyyyyyz_1, g_yyzzz_0_xyyyyyz_1, g_yyzzz_0_xyyyyyzz_1, g_yyzzz_0_xyyyyzz_1, g_yyzzz_0_xyyyyzzz_1, g_yyzzz_0_xyyyzzz_1, g_yyzzz_0_xyyyzzzz_1, g_yyzzz_0_xyyzzzz_1, g_yyzzz_0_xyyzzzzz_1, g_yyzzz_0_xyzzzzz_1, g_yyzzz_0_xyzzzzzz_1, g_yyzzz_0_xzzzzzz_1, g_yyzzz_0_xzzzzzzz_1, g_yyzzz_0_yyyyyyyz_1, g_yyzzz_0_yyyyyyz_1, g_yyzzz_0_yyyyyyzz_1, g_yyzzz_0_yyyyyzz_1, g_yyzzz_0_yyyyyzzz_1, g_yyzzz_0_yyyyzzz_1, g_yyzzz_0_yyyyzzzz_1, g_yyzzz_0_yyyzzzz_1, g_yyzzz_0_yyyzzzzz_1, g_yyzzz_0_yyzzzzz_1, g_yyzzz_0_yyzzzzzz_1, g_yyzzz_0_yzzzzzz_1, g_yyzzz_0_yzzzzzzz_1, g_yyzzz_0_zzzzzzz_1, g_yyzzz_0_zzzzzzzz_1, g_yzzz_0_xxxxxxxx_0, g_yzzz_0_xxxxxxxx_1, g_yzzz_0_xxxxxxxz_0, g_yzzz_0_xxxxxxxz_1, g_yzzz_0_xxxxxxyz_0, g_yzzz_0_xxxxxxyz_1, g_yzzz_0_xxxxxxzz_0, g_yzzz_0_xxxxxxzz_1, g_yzzz_0_xxxxxyyz_0, g_yzzz_0_xxxxxyyz_1, g_yzzz_0_xxxxxyzz_0, g_yzzz_0_xxxxxyzz_1, g_yzzz_0_xxxxxzzz_0, g_yzzz_0_xxxxxzzz_1, g_yzzz_0_xxxxyyyz_0, g_yzzz_0_xxxxyyyz_1, g_yzzz_0_xxxxyyzz_0, g_yzzz_0_xxxxyyzz_1, g_yzzz_0_xxxxyzzz_0, g_yzzz_0_xxxxyzzz_1, g_yzzz_0_xxxxzzzz_0, g_yzzz_0_xxxxzzzz_1, g_yzzz_0_xxxyyyyz_0, g_yzzz_0_xxxyyyyz_1, g_yzzz_0_xxxyyyzz_0, g_yzzz_0_xxxyyyzz_1, g_yzzz_0_xxxyyzzz_0, g_yzzz_0_xxxyyzzz_1, g_yzzz_0_xxxyzzzz_0, g_yzzz_0_xxxyzzzz_1, g_yzzz_0_xxxzzzzz_0, g_yzzz_0_xxxzzzzz_1, g_yzzz_0_xxyyyyyz_0, g_yzzz_0_xxyyyyyz_1, g_yzzz_0_xxyyyyzz_0, g_yzzz_0_xxyyyyzz_1, g_yzzz_0_xxyyyzzz_0, g_yzzz_0_xxyyyzzz_1, g_yzzz_0_xxyyzzzz_0, g_yzzz_0_xxyyzzzz_1, g_yzzz_0_xxyzzzzz_0, g_yzzz_0_xxyzzzzz_1, g_yzzz_0_xxzzzzzz_0, g_yzzz_0_xxzzzzzz_1, g_yzzz_0_xyyyyyyz_0, g_yzzz_0_xyyyyyyz_1, g_yzzz_0_xyyyyyzz_0, g_yzzz_0_xyyyyyzz_1, g_yzzz_0_xyyyyzzz_0, g_yzzz_0_xyyyyzzz_1, g_yzzz_0_xyyyzzzz_0, g_yzzz_0_xyyyzzzz_1, g_yzzz_0_xyyzzzzz_0, g_yzzz_0_xyyzzzzz_1, g_yzzz_0_xyzzzzzz_0, g_yzzz_0_xyzzzzzz_1, g_yzzz_0_xzzzzzzz_0, g_yzzz_0_xzzzzzzz_1, g_yzzz_0_yyyyyyyz_0, g_yzzz_0_yyyyyyyz_1, g_yzzz_0_yyyyyyzz_0, g_yzzz_0_yyyyyyzz_1, g_yzzz_0_yyyyyzzz_0, g_yzzz_0_yyyyyzzz_1, g_yzzz_0_yyyyzzzz_0, g_yzzz_0_yyyyzzzz_1, g_yzzz_0_yyyzzzzz_0, g_yzzz_0_yyyzzzzz_1, g_yzzz_0_yyzzzzzz_0, g_yzzz_0_yyzzzzzz_1, g_yzzz_0_yzzzzzzz_0, g_yzzz_0_yzzzzzzz_1, g_yzzz_0_zzzzzzzz_0, g_yzzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xxxxxxxx_0[i] = 2.0 * g_yzzz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxxxy_0[i] = 2.0 * g_yyyz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxxxxz_0[i] = 2.0 * g_yzzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxxyy_0[i] = 2.0 * g_yyyz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxxxyy_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxxxyz_0[i] = 2.0 * g_yzzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxxzz_0[i] = 2.0 * g_yzzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxyyy_0[i] = 2.0 * g_yyyz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxxyyy_1[i] * fz_be_0 + g_yyyzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxxyyz_0[i] = 2.0 * g_yzzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxyzz_0[i] = 2.0 * g_yzzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxzzz_0[i] = 2.0 * g_yzzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxyyyy_0[i] = 2.0 * g_yyyz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxyyyy_1[i] * fz_be_0 + g_yyyzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxyyyz_0[i] = 2.0 * g_yzzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxyyzz_0[i] = 2.0 * g_yzzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxyzzz_0[i] = 2.0 * g_yzzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxzzzz_0[i] = 2.0 * g_yzzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyyyyy_0[i] = 2.0 * g_yyyz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxyyyyy_1[i] * fz_be_0 + g_yyyzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxyyyyz_0[i] = 2.0 * g_yzzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyyyzz_0[i] = 2.0 * g_yzzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyyzzz_0[i] = 2.0 * g_yzzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyzzzz_0[i] = 2.0 * g_yzzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxzzzzz_0[i] = 2.0 * g_yzzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyyyyy_0[i] = 2.0 * g_yyyz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxyyyyyy_1[i] * fz_be_0 + g_yyyzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxyyyyyz_0[i] = 2.0 * g_yzzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyyyzz_0[i] = 2.0 * g_yzzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyyzzz_0[i] = 2.0 * g_yzzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyzzzz_0[i] = 2.0 * g_yzzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyzzzzz_0[i] = 2.0 * g_yzzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxzzzzzz_0[i] = 2.0 * g_yzzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyyyyy_0[i] = 2.0 * g_yyyz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xyyyyyyy_1[i] * fz_be_0 + g_yyyzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xyyyyyyz_0[i] = 2.0 * g_yzzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyyyzz_0[i] = 2.0 * g_yzzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyyzzz_0[i] = 2.0 * g_yzzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyzzzz_0[i] = 2.0 * g_yzzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyzzzzz_0[i] = 2.0 * g_yzzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyzzzzzz_0[i] = 2.0 * g_yzzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xzzzzzzz_0[i] = 2.0 * g_yzzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyyyyy_0[i] = 2.0 * g_yyyz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yyyyyyyy_1[i] * fz_be_0 + g_yyyzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_yyyyyyyz_0[i] = 2.0 * g_yzzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyyyzz_0[i] = 2.0 * g_yzzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyyzzz_0[i] = 2.0 * g_yzzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyzzzz_0[i] = 2.0 * g_yzzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyzzzzz_0[i] = 2.0 * g_yzzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyzzzzzz_0[i] = 2.0 * g_yzzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yzzzzzzz_0[i] = 2.0 * g_yzzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_zzzzzzzz_0[i] = 2.0 * g_yzzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1125-1170 components of targeted buffer : ISL

    auto g_yyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 1125);

    auto g_yyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 1126);

    auto g_yyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 1127);

    auto g_yyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 1128);

    auto g_yyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 1129);

    auto g_yyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 1130);

    auto g_yyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 1131);

    auto g_yyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 1132);

    auto g_yyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 1133);

    auto g_yyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 1134);

    auto g_yyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 1135);

    auto g_yyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 1136);

    auto g_yyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 1137);

    auto g_yyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 1138);

    auto g_yyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 1139);

    auto g_yyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1140);

    auto g_yyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1141);

    auto g_yyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1142);

    auto g_yyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1143);

    auto g_yyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1144);

    auto g_yyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1145);

    auto g_yyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1146);

    auto g_yyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1147);

    auto g_yyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1148);

    auto g_yyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1149);

    auto g_yyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1150);

    auto g_yyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1151);

    auto g_yyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1152);

    auto g_yyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1153);

    auto g_yyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1154);

    auto g_yyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1155);

    auto g_yyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1156);

    auto g_yyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1157);

    auto g_yyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1158);

    auto g_yyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1159);

    auto g_yyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1160);

    auto g_yyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1161);

    auto g_yyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1162);

    auto g_yyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1163);

    auto g_yyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1164);

    auto g_yyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1165);

    auto g_yyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1166);

    auto g_yyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1167);

    auto g_yyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1168);

    auto g_yyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1169);

    #pragma omp simd aligned(g_yyzz_0_xxxxxxxy_0, g_yyzz_0_xxxxxxxy_1, g_yyzz_0_xxxxxxyy_0, g_yyzz_0_xxxxxxyy_1, g_yyzz_0_xxxxxyyy_0, g_yyzz_0_xxxxxyyy_1, g_yyzz_0_xxxxyyyy_0, g_yyzz_0_xxxxyyyy_1, g_yyzz_0_xxxyyyyy_0, g_yyzz_0_xxxyyyyy_1, g_yyzz_0_xxyyyyyy_0, g_yyzz_0_xxyyyyyy_1, g_yyzz_0_xyyyyyyy_0, g_yyzz_0_xyyyyyyy_1, g_yyzz_0_yyyyyyyy_0, g_yyzz_0_yyyyyyyy_1, g_yyzzz_0_xxxxxxxy_1, g_yyzzz_0_xxxxxxyy_1, g_yyzzz_0_xxxxxyyy_1, g_yyzzz_0_xxxxyyyy_1, g_yyzzz_0_xxxyyyyy_1, g_yyzzz_0_xxyyyyyy_1, g_yyzzz_0_xyyyyyyy_1, g_yyzzz_0_yyyyyyyy_1, g_yyzzzz_0_xxxxxxxx_0, g_yyzzzz_0_xxxxxxxy_0, g_yyzzzz_0_xxxxxxxz_0, g_yyzzzz_0_xxxxxxyy_0, g_yyzzzz_0_xxxxxxyz_0, g_yyzzzz_0_xxxxxxzz_0, g_yyzzzz_0_xxxxxyyy_0, g_yyzzzz_0_xxxxxyyz_0, g_yyzzzz_0_xxxxxyzz_0, g_yyzzzz_0_xxxxxzzz_0, g_yyzzzz_0_xxxxyyyy_0, g_yyzzzz_0_xxxxyyyz_0, g_yyzzzz_0_xxxxyyzz_0, g_yyzzzz_0_xxxxyzzz_0, g_yyzzzz_0_xxxxzzzz_0, g_yyzzzz_0_xxxyyyyy_0, g_yyzzzz_0_xxxyyyyz_0, g_yyzzzz_0_xxxyyyzz_0, g_yyzzzz_0_xxxyyzzz_0, g_yyzzzz_0_xxxyzzzz_0, g_yyzzzz_0_xxxzzzzz_0, g_yyzzzz_0_xxyyyyyy_0, g_yyzzzz_0_xxyyyyyz_0, g_yyzzzz_0_xxyyyyzz_0, g_yyzzzz_0_xxyyyzzz_0, g_yyzzzz_0_xxyyzzzz_0, g_yyzzzz_0_xxyzzzzz_0, g_yyzzzz_0_xxzzzzzz_0, g_yyzzzz_0_xyyyyyyy_0, g_yyzzzz_0_xyyyyyyz_0, g_yyzzzz_0_xyyyyyzz_0, g_yyzzzz_0_xyyyyzzz_0, g_yyzzzz_0_xyyyzzzz_0, g_yyzzzz_0_xyyzzzzz_0, g_yyzzzz_0_xyzzzzzz_0, g_yyzzzz_0_xzzzzzzz_0, g_yyzzzz_0_yyyyyyyy_0, g_yyzzzz_0_yyyyyyyz_0, g_yyzzzz_0_yyyyyyzz_0, g_yyzzzz_0_yyyyyzzz_0, g_yyzzzz_0_yyyyzzzz_0, g_yyzzzz_0_yyyzzzzz_0, g_yyzzzz_0_yyzzzzzz_0, g_yyzzzz_0_yzzzzzzz_0, g_yyzzzz_0_zzzzzzzz_0, g_yzzzz_0_xxxxxxxx_1, g_yzzzz_0_xxxxxxxz_1, g_yzzzz_0_xxxxxxyz_1, g_yzzzz_0_xxxxxxz_1, g_yzzzz_0_xxxxxxzz_1, g_yzzzz_0_xxxxxyyz_1, g_yzzzz_0_xxxxxyz_1, g_yzzzz_0_xxxxxyzz_1, g_yzzzz_0_xxxxxzz_1, g_yzzzz_0_xxxxxzzz_1, g_yzzzz_0_xxxxyyyz_1, g_yzzzz_0_xxxxyyz_1, g_yzzzz_0_xxxxyyzz_1, g_yzzzz_0_xxxxyzz_1, g_yzzzz_0_xxxxyzzz_1, g_yzzzz_0_xxxxzzz_1, g_yzzzz_0_xxxxzzzz_1, g_yzzzz_0_xxxyyyyz_1, g_yzzzz_0_xxxyyyz_1, g_yzzzz_0_xxxyyyzz_1, g_yzzzz_0_xxxyyzz_1, g_yzzzz_0_xxxyyzzz_1, g_yzzzz_0_xxxyzzz_1, g_yzzzz_0_xxxyzzzz_1, g_yzzzz_0_xxxzzzz_1, g_yzzzz_0_xxxzzzzz_1, g_yzzzz_0_xxyyyyyz_1, g_yzzzz_0_xxyyyyz_1, g_yzzzz_0_xxyyyyzz_1, g_yzzzz_0_xxyyyzz_1, g_yzzzz_0_xxyyyzzz_1, g_yzzzz_0_xxyyzzz_1, g_yzzzz_0_xxyyzzzz_1, g_yzzzz_0_xxyzzzz_1, g_yzzzz_0_xxyzzzzz_1, g_yzzzz_0_xxzzzzz_1, g_yzzzz_0_xxzzzzzz_1, g_yzzzz_0_xyyyyyyz_1, g_yzzzz_0_xyyyyyz_1, g_yzzzz_0_xyyyyyzz_1, g_yzzzz_0_xyyyyzz_1, g_yzzzz_0_xyyyyzzz_1, g_yzzzz_0_xyyyzzz_1, g_yzzzz_0_xyyyzzzz_1, g_yzzzz_0_xyyzzzz_1, g_yzzzz_0_xyyzzzzz_1, g_yzzzz_0_xyzzzzz_1, g_yzzzz_0_xyzzzzzz_1, g_yzzzz_0_xzzzzzz_1, g_yzzzz_0_xzzzzzzz_1, g_yzzzz_0_yyyyyyyz_1, g_yzzzz_0_yyyyyyz_1, g_yzzzz_0_yyyyyyzz_1, g_yzzzz_0_yyyyyzz_1, g_yzzzz_0_yyyyyzzz_1, g_yzzzz_0_yyyyzzz_1, g_yzzzz_0_yyyyzzzz_1, g_yzzzz_0_yyyzzzz_1, g_yzzzz_0_yyyzzzzz_1, g_yzzzz_0_yyzzzzz_1, g_yzzzz_0_yyzzzzzz_1, g_yzzzz_0_yzzzzzz_1, g_yzzzz_0_yzzzzzzz_1, g_yzzzz_0_zzzzzzz_1, g_yzzzz_0_zzzzzzzz_1, g_zzzz_0_xxxxxxxx_0, g_zzzz_0_xxxxxxxx_1, g_zzzz_0_xxxxxxxz_0, g_zzzz_0_xxxxxxxz_1, g_zzzz_0_xxxxxxyz_0, g_zzzz_0_xxxxxxyz_1, g_zzzz_0_xxxxxxzz_0, g_zzzz_0_xxxxxxzz_1, g_zzzz_0_xxxxxyyz_0, g_zzzz_0_xxxxxyyz_1, g_zzzz_0_xxxxxyzz_0, g_zzzz_0_xxxxxyzz_1, g_zzzz_0_xxxxxzzz_0, g_zzzz_0_xxxxxzzz_1, g_zzzz_0_xxxxyyyz_0, g_zzzz_0_xxxxyyyz_1, g_zzzz_0_xxxxyyzz_0, g_zzzz_0_xxxxyyzz_1, g_zzzz_0_xxxxyzzz_0, g_zzzz_0_xxxxyzzz_1, g_zzzz_0_xxxxzzzz_0, g_zzzz_0_xxxxzzzz_1, g_zzzz_0_xxxyyyyz_0, g_zzzz_0_xxxyyyyz_1, g_zzzz_0_xxxyyyzz_0, g_zzzz_0_xxxyyyzz_1, g_zzzz_0_xxxyyzzz_0, g_zzzz_0_xxxyyzzz_1, g_zzzz_0_xxxyzzzz_0, g_zzzz_0_xxxyzzzz_1, g_zzzz_0_xxxzzzzz_0, g_zzzz_0_xxxzzzzz_1, g_zzzz_0_xxyyyyyz_0, g_zzzz_0_xxyyyyyz_1, g_zzzz_0_xxyyyyzz_0, g_zzzz_0_xxyyyyzz_1, g_zzzz_0_xxyyyzzz_0, g_zzzz_0_xxyyyzzz_1, g_zzzz_0_xxyyzzzz_0, g_zzzz_0_xxyyzzzz_1, g_zzzz_0_xxyzzzzz_0, g_zzzz_0_xxyzzzzz_1, g_zzzz_0_xxzzzzzz_0, g_zzzz_0_xxzzzzzz_1, g_zzzz_0_xyyyyyyz_0, g_zzzz_0_xyyyyyyz_1, g_zzzz_0_xyyyyyzz_0, g_zzzz_0_xyyyyyzz_1, g_zzzz_0_xyyyyzzz_0, g_zzzz_0_xyyyyzzz_1, g_zzzz_0_xyyyzzzz_0, g_zzzz_0_xyyyzzzz_1, g_zzzz_0_xyyzzzzz_0, g_zzzz_0_xyyzzzzz_1, g_zzzz_0_xyzzzzzz_0, g_zzzz_0_xyzzzzzz_1, g_zzzz_0_xzzzzzzz_0, g_zzzz_0_xzzzzzzz_1, g_zzzz_0_yyyyyyyz_0, g_zzzz_0_yyyyyyyz_1, g_zzzz_0_yyyyyyzz_0, g_zzzz_0_yyyyyyzz_1, g_zzzz_0_yyyyyzzz_0, g_zzzz_0_yyyyyzzz_1, g_zzzz_0_yyyyzzzz_0, g_zzzz_0_yyyyzzzz_1, g_zzzz_0_yyyzzzzz_0, g_zzzz_0_yyyzzzzz_1, g_zzzz_0_yyzzzzzz_0, g_zzzz_0_yyzzzzzz_1, g_zzzz_0_yzzzzzzz_0, g_zzzz_0_yzzzzzzz_1, g_zzzz_0_zzzzzzzz_0, g_zzzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xxxxxxxx_0[i] = g_zzzz_0_xxxxxxxx_0[i] * fbe_0 - g_zzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxxxy_0[i] = 3.0 * g_yyzz_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxxy_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxxxxz_0[i] = g_zzzz_0_xxxxxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxxyy_0[i] = 3.0 * g_yyzz_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxyy_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxxxyz_0[i] = g_zzzz_0_xxxxxxyz_0[i] * fbe_0 - g_zzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxxzz_0[i] = g_zzzz_0_xxxxxxzz_0[i] * fbe_0 - g_zzzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxyyy_0[i] = 3.0 * g_yyzz_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxyyy_1[i] * fz_be_0 + g_yyzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxxyyz_0[i] = g_zzzz_0_xxxxxyyz_0[i] * fbe_0 - g_zzzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxyzz_0[i] = g_zzzz_0_xxxxxyzz_0[i] * fbe_0 - g_zzzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxzzz_0[i] = g_zzzz_0_xxxxxzzz_0[i] * fbe_0 - g_zzzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxyyyy_0[i] = 3.0 * g_yyzz_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyyyy_1[i] * fz_be_0 + g_yyzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxyyyz_0[i] = g_zzzz_0_xxxxyyyz_0[i] * fbe_0 - g_zzzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxyyzz_0[i] = g_zzzz_0_xxxxyyzz_0[i] * fbe_0 - g_zzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxyzzz_0[i] = g_zzzz_0_xxxxyzzz_0[i] * fbe_0 - g_zzzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxzzzz_0[i] = g_zzzz_0_xxxxzzzz_0[i] * fbe_0 - g_zzzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyyyyy_0[i] = 3.0 * g_yyzz_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyyyy_1[i] * fz_be_0 + g_yyzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxyyyyz_0[i] = g_zzzz_0_xxxyyyyz_0[i] * fbe_0 - g_zzzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyyyzz_0[i] = g_zzzz_0_xxxyyyzz_0[i] * fbe_0 - g_zzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyyzzz_0[i] = g_zzzz_0_xxxyyzzz_0[i] * fbe_0 - g_zzzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyzzzz_0[i] = g_zzzz_0_xxxyzzzz_0[i] * fbe_0 - g_zzzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxzzzzz_0[i] = g_zzzz_0_xxxzzzzz_0[i] * fbe_0 - g_zzzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyyyyy_0[i] = 3.0 * g_yyzz_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyyyy_1[i] * fz_be_0 + g_yyzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxyyyyyz_0[i] = g_zzzz_0_xxyyyyyz_0[i] * fbe_0 - g_zzzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyyyzz_0[i] = g_zzzz_0_xxyyyyzz_0[i] * fbe_0 - g_zzzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyyzzz_0[i] = g_zzzz_0_xxyyyzzz_0[i] * fbe_0 - g_zzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyzzzz_0[i] = g_zzzz_0_xxyyzzzz_0[i] * fbe_0 - g_zzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyzzzzz_0[i] = g_zzzz_0_xxyzzzzz_0[i] * fbe_0 - g_zzzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxzzzzzz_0[i] = g_zzzz_0_xxzzzzzz_0[i] * fbe_0 - g_zzzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyyyyy_0[i] = 3.0 * g_yyzz_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyyyy_1[i] * fz_be_0 + g_yyzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xyyyyyyz_0[i] = g_zzzz_0_xyyyyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyyyzz_0[i] = g_zzzz_0_xyyyyyzz_0[i] * fbe_0 - g_zzzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyyzzz_0[i] = g_zzzz_0_xyyyyzzz_0[i] * fbe_0 - g_zzzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyzzzz_0[i] = g_zzzz_0_xyyyzzzz_0[i] * fbe_0 - g_zzzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyzzzzz_0[i] = g_zzzz_0_xyyzzzzz_0[i] * fbe_0 - g_zzzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyzzzzzz_0[i] = g_zzzz_0_xyzzzzzz_0[i] * fbe_0 - g_zzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xzzzzzzz_0[i] = g_zzzz_0_xzzzzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyyyyy_0[i] = 3.0 * g_yyzz_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyyyy_1[i] * fz_be_0 + g_yyzzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_yyyyyyyz_0[i] = g_zzzz_0_yyyyyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyyyzz_0[i] = g_zzzz_0_yyyyyyzz_0[i] * fbe_0 - g_zzzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyyzzz_0[i] = g_zzzz_0_yyyyyzzz_0[i] * fbe_0 - g_zzzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyzzzz_0[i] = g_zzzz_0_yyyyzzzz_0[i] * fbe_0 - g_zzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyzzzzz_0[i] = g_zzzz_0_yyyzzzzz_0[i] * fbe_0 - g_zzzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyzzzzzz_0[i] = g_zzzz_0_yyzzzzzz_0[i] * fbe_0 - g_zzzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yzzzzzzz_0[i] = g_zzzz_0_yzzzzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_zzzzzzzz_0[i] = g_zzzz_0_zzzzzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1170-1215 components of targeted buffer : ISL

    auto g_yzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 1170);

    auto g_yzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 1171);

    auto g_yzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 1172);

    auto g_yzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 1173);

    auto g_yzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 1174);

    auto g_yzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 1175);

    auto g_yzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 1176);

    auto g_yzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 1177);

    auto g_yzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 1178);

    auto g_yzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 1179);

    auto g_yzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 1180);

    auto g_yzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 1181);

    auto g_yzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 1182);

    auto g_yzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 1183);

    auto g_yzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 1184);

    auto g_yzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1185);

    auto g_yzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1186);

    auto g_yzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1187);

    auto g_yzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1188);

    auto g_yzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1189);

    auto g_yzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1190);

    auto g_yzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1191);

    auto g_yzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1192);

    auto g_yzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1193);

    auto g_yzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1194);

    auto g_yzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1195);

    auto g_yzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1196);

    auto g_yzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1197);

    auto g_yzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1198);

    auto g_yzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1199);

    auto g_yzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1200);

    auto g_yzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1201);

    auto g_yzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1202);

    auto g_yzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1203);

    auto g_yzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1204);

    auto g_yzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1205);

    auto g_yzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1206);

    auto g_yzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1207);

    auto g_yzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1208);

    auto g_yzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1209);

    auto g_yzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1210);

    auto g_yzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1211);

    auto g_yzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1212);

    auto g_yzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1213);

    auto g_yzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1214);

    #pragma omp simd aligned(g_yzzzzz_0_xxxxxxxx_0, g_yzzzzz_0_xxxxxxxy_0, g_yzzzzz_0_xxxxxxxz_0, g_yzzzzz_0_xxxxxxyy_0, g_yzzzzz_0_xxxxxxyz_0, g_yzzzzz_0_xxxxxxzz_0, g_yzzzzz_0_xxxxxyyy_0, g_yzzzzz_0_xxxxxyyz_0, g_yzzzzz_0_xxxxxyzz_0, g_yzzzzz_0_xxxxxzzz_0, g_yzzzzz_0_xxxxyyyy_0, g_yzzzzz_0_xxxxyyyz_0, g_yzzzzz_0_xxxxyyzz_0, g_yzzzzz_0_xxxxyzzz_0, g_yzzzzz_0_xxxxzzzz_0, g_yzzzzz_0_xxxyyyyy_0, g_yzzzzz_0_xxxyyyyz_0, g_yzzzzz_0_xxxyyyzz_0, g_yzzzzz_0_xxxyyzzz_0, g_yzzzzz_0_xxxyzzzz_0, g_yzzzzz_0_xxxzzzzz_0, g_yzzzzz_0_xxyyyyyy_0, g_yzzzzz_0_xxyyyyyz_0, g_yzzzzz_0_xxyyyyzz_0, g_yzzzzz_0_xxyyyzzz_0, g_yzzzzz_0_xxyyzzzz_0, g_yzzzzz_0_xxyzzzzz_0, g_yzzzzz_0_xxzzzzzz_0, g_yzzzzz_0_xyyyyyyy_0, g_yzzzzz_0_xyyyyyyz_0, g_yzzzzz_0_xyyyyyzz_0, g_yzzzzz_0_xyyyyzzz_0, g_yzzzzz_0_xyyyzzzz_0, g_yzzzzz_0_xyyzzzzz_0, g_yzzzzz_0_xyzzzzzz_0, g_yzzzzz_0_xzzzzzzz_0, g_yzzzzz_0_yyyyyyyy_0, g_yzzzzz_0_yyyyyyyz_0, g_yzzzzz_0_yyyyyyzz_0, g_yzzzzz_0_yyyyyzzz_0, g_yzzzzz_0_yyyyzzzz_0, g_yzzzzz_0_yyyzzzzz_0, g_yzzzzz_0_yyzzzzzz_0, g_yzzzzz_0_yzzzzzzz_0, g_yzzzzz_0_zzzzzzzz_0, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxxx_1, g_zzzzz_0_xxxxxxxy_1, g_zzzzz_0_xxxxxxxz_1, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxyy_1, g_zzzzz_0_xxxxxxyz_1, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxxzz_1, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyyy_1, g_zzzzz_0_xxxxxyyz_1, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxyzz_1, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxxzzz_1, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyyy_1, g_zzzzz_0_xxxxyyyz_1, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyyzz_1, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxyzzz_1, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxxzzzz_1, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyyy_1, g_zzzzz_0_xxxyyyyz_1, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyyzz_1, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyyzzz_1, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxyzzzz_1, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxxzzzzz_1, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyyy_1, g_zzzzz_0_xxyyyyyz_1, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyyzz_1, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyyzzz_1, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyyzzzz_1, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxyzzzzz_1, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xxzzzzzz_1, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyyy_1, g_zzzzz_0_xyyyyyyz_1, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyyzz_1, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyyzzz_1, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyyzzzz_1, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyyzzzzz_1, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xyzzzzzz_1, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_xzzzzzzz_1, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyyy_1, g_zzzzz_0_yyyyyyyz_1, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyyzz_1, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyyzzz_1, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyyzzzz_1, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyyzzzzz_1, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yyzzzzzz_1, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_yzzzzzzz_1, g_zzzzz_0_zzzzzzz_1, g_zzzzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xxxxxxxx_0[i] = g_zzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxxy_0[i] = g_zzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxxz_0[i] = g_zzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxyy_0[i] = 2.0 * g_zzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxyz_0[i] = g_zzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxzz_0[i] = g_zzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxyyy_0[i] = 3.0 * g_zzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxyyz_0[i] = 2.0 * g_zzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxyzz_0[i] = g_zzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxzzz_0[i] = g_zzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyyyz_0[i] = 3.0 * g_zzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyyzz_0[i] = 2.0 * g_zzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyzzz_0[i] = g_zzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxzzzz_0[i] = g_zzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyyyy_0[i] = 5.0 * g_zzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyyyz_0[i] = 4.0 * g_zzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyzzz_0[i] = 2.0 * g_zzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyzzzz_0[i] = g_zzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxzzzzz_0[i] = g_zzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyyyy_0[i] = 6.0 * g_zzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyyyz_0[i] = 5.0 * g_zzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyyzz_0[i] = 4.0 * g_zzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyzzz_0[i] = 3.0 * g_zzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyzzzzz_0[i] = g_zzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxzzzzzz_0[i] = g_zzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyyyy_0[i] = 7.0 * g_zzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyyyz_0[i] = 6.0 * g_zzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyyzz_0[i] = 5.0 * g_zzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyzzz_0[i] = 4.0 * g_zzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyzzzz_0[i] = 3.0 * g_zzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyzzzzz_0[i] = 2.0 * g_zzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyzzzzzz_0[i] = g_zzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xzzzzzzz_0[i] = g_zzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyyyy_0[i] = 8.0 * g_zzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyyyz_0[i] = 7.0 * g_zzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyyzz_0[i] = 6.0 * g_zzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyzzz_0[i] = 5.0 * g_zzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyzzzz_0[i] = 4.0 * g_zzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyzzzzz_0[i] = 3.0 * g_zzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyzzzzzz_0[i] = 2.0 * g_zzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yzzzzzzz_0[i] = g_zzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_zzzzzzzz_0[i] = g_zzzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1215-1260 components of targeted buffer : ISL

    auto g_zzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_isl + 1215);

    auto g_zzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_isl + 1216);

    auto g_zzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_isl + 1217);

    auto g_zzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_isl + 1218);

    auto g_zzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_isl + 1219);

    auto g_zzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_isl + 1220);

    auto g_zzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_isl + 1221);

    auto g_zzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_isl + 1222);

    auto g_zzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_isl + 1223);

    auto g_zzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_isl + 1224);

    auto g_zzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_isl + 1225);

    auto g_zzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_isl + 1226);

    auto g_zzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_isl + 1227);

    auto g_zzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_isl + 1228);

    auto g_zzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_isl + 1229);

    auto g_zzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1230);

    auto g_zzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1231);

    auto g_zzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1232);

    auto g_zzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1233);

    auto g_zzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1234);

    auto g_zzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1235);

    auto g_zzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1236);

    auto g_zzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1237);

    auto g_zzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1238);

    auto g_zzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1239);

    auto g_zzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1240);

    auto g_zzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1241);

    auto g_zzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1242);

    auto g_zzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1243);

    auto g_zzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1244);

    auto g_zzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1245);

    auto g_zzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1246);

    auto g_zzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1247);

    auto g_zzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1248);

    auto g_zzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1249);

    auto g_zzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1250);

    auto g_zzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_isl + 1251);

    auto g_zzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_isl + 1252);

    auto g_zzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_isl + 1253);

    auto g_zzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_isl + 1254);

    auto g_zzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_isl + 1255);

    auto g_zzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1256);

    auto g_zzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1257);

    auto g_zzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1258);

    auto g_zzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_isl + 1259);

    #pragma omp simd aligned(g_zzzz_0_xxxxxxxx_0, g_zzzz_0_xxxxxxxx_1, g_zzzz_0_xxxxxxxy_0, g_zzzz_0_xxxxxxxy_1, g_zzzz_0_xxxxxxxz_0, g_zzzz_0_xxxxxxxz_1, g_zzzz_0_xxxxxxyy_0, g_zzzz_0_xxxxxxyy_1, g_zzzz_0_xxxxxxyz_0, g_zzzz_0_xxxxxxyz_1, g_zzzz_0_xxxxxxzz_0, g_zzzz_0_xxxxxxzz_1, g_zzzz_0_xxxxxyyy_0, g_zzzz_0_xxxxxyyy_1, g_zzzz_0_xxxxxyyz_0, g_zzzz_0_xxxxxyyz_1, g_zzzz_0_xxxxxyzz_0, g_zzzz_0_xxxxxyzz_1, g_zzzz_0_xxxxxzzz_0, g_zzzz_0_xxxxxzzz_1, g_zzzz_0_xxxxyyyy_0, g_zzzz_0_xxxxyyyy_1, g_zzzz_0_xxxxyyyz_0, g_zzzz_0_xxxxyyyz_1, g_zzzz_0_xxxxyyzz_0, g_zzzz_0_xxxxyyzz_1, g_zzzz_0_xxxxyzzz_0, g_zzzz_0_xxxxyzzz_1, g_zzzz_0_xxxxzzzz_0, g_zzzz_0_xxxxzzzz_1, g_zzzz_0_xxxyyyyy_0, g_zzzz_0_xxxyyyyy_1, g_zzzz_0_xxxyyyyz_0, g_zzzz_0_xxxyyyyz_1, g_zzzz_0_xxxyyyzz_0, g_zzzz_0_xxxyyyzz_1, g_zzzz_0_xxxyyzzz_0, g_zzzz_0_xxxyyzzz_1, g_zzzz_0_xxxyzzzz_0, g_zzzz_0_xxxyzzzz_1, g_zzzz_0_xxxzzzzz_0, g_zzzz_0_xxxzzzzz_1, g_zzzz_0_xxyyyyyy_0, g_zzzz_0_xxyyyyyy_1, g_zzzz_0_xxyyyyyz_0, g_zzzz_0_xxyyyyyz_1, g_zzzz_0_xxyyyyzz_0, g_zzzz_0_xxyyyyzz_1, g_zzzz_0_xxyyyzzz_0, g_zzzz_0_xxyyyzzz_1, g_zzzz_0_xxyyzzzz_0, g_zzzz_0_xxyyzzzz_1, g_zzzz_0_xxyzzzzz_0, g_zzzz_0_xxyzzzzz_1, g_zzzz_0_xxzzzzzz_0, g_zzzz_0_xxzzzzzz_1, g_zzzz_0_xyyyyyyy_0, g_zzzz_0_xyyyyyyy_1, g_zzzz_0_xyyyyyyz_0, g_zzzz_0_xyyyyyyz_1, g_zzzz_0_xyyyyyzz_0, g_zzzz_0_xyyyyyzz_1, g_zzzz_0_xyyyyzzz_0, g_zzzz_0_xyyyyzzz_1, g_zzzz_0_xyyyzzzz_0, g_zzzz_0_xyyyzzzz_1, g_zzzz_0_xyyzzzzz_0, g_zzzz_0_xyyzzzzz_1, g_zzzz_0_xyzzzzzz_0, g_zzzz_0_xyzzzzzz_1, g_zzzz_0_xzzzzzzz_0, g_zzzz_0_xzzzzzzz_1, g_zzzz_0_yyyyyyyy_0, g_zzzz_0_yyyyyyyy_1, g_zzzz_0_yyyyyyyz_0, g_zzzz_0_yyyyyyyz_1, g_zzzz_0_yyyyyyzz_0, g_zzzz_0_yyyyyyzz_1, g_zzzz_0_yyyyyzzz_0, g_zzzz_0_yyyyyzzz_1, g_zzzz_0_yyyyzzzz_0, g_zzzz_0_yyyyzzzz_1, g_zzzz_0_yyyzzzzz_0, g_zzzz_0_yyyzzzzz_1, g_zzzz_0_yyzzzzzz_0, g_zzzz_0_yyzzzzzz_1, g_zzzz_0_yzzzzzzz_0, g_zzzz_0_yzzzzzzz_1, g_zzzz_0_zzzzzzzz_0, g_zzzz_0_zzzzzzzz_1, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxxx_1, g_zzzzz_0_xxxxxxxy_1, g_zzzzz_0_xxxxxxxz_1, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxyy_1, g_zzzzz_0_xxxxxxyz_1, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxxzz_1, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyyy_1, g_zzzzz_0_xxxxxyyz_1, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxyzz_1, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxxzzz_1, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyyy_1, g_zzzzz_0_xxxxyyyz_1, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyyzz_1, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxyzzz_1, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxxzzzz_1, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyyy_1, g_zzzzz_0_xxxyyyyz_1, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyyzz_1, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyyzzz_1, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxyzzzz_1, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxxzzzzz_1, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyyy_1, g_zzzzz_0_xxyyyyyz_1, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyyzz_1, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyyzzz_1, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyyzzzz_1, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxyzzzzz_1, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xxzzzzzz_1, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyyy_1, g_zzzzz_0_xyyyyyyz_1, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyyzz_1, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyyzzz_1, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyyzzzz_1, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyyzzzzz_1, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xyzzzzzz_1, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_xzzzzzzz_1, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyyy_1, g_zzzzz_0_yyyyyyyz_1, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyyzz_1, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyyzzz_1, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyyzzzz_1, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyyzzzzz_1, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yyzzzzzz_1, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_yzzzzzzz_1, g_zzzzz_0_zzzzzzz_1, g_zzzzz_0_zzzzzzzz_1, g_zzzzzz_0_xxxxxxxx_0, g_zzzzzz_0_xxxxxxxy_0, g_zzzzzz_0_xxxxxxxz_0, g_zzzzzz_0_xxxxxxyy_0, g_zzzzzz_0_xxxxxxyz_0, g_zzzzzz_0_xxxxxxzz_0, g_zzzzzz_0_xxxxxyyy_0, g_zzzzzz_0_xxxxxyyz_0, g_zzzzzz_0_xxxxxyzz_0, g_zzzzzz_0_xxxxxzzz_0, g_zzzzzz_0_xxxxyyyy_0, g_zzzzzz_0_xxxxyyyz_0, g_zzzzzz_0_xxxxyyzz_0, g_zzzzzz_0_xxxxyzzz_0, g_zzzzzz_0_xxxxzzzz_0, g_zzzzzz_0_xxxyyyyy_0, g_zzzzzz_0_xxxyyyyz_0, g_zzzzzz_0_xxxyyyzz_0, g_zzzzzz_0_xxxyyzzz_0, g_zzzzzz_0_xxxyzzzz_0, g_zzzzzz_0_xxxzzzzz_0, g_zzzzzz_0_xxyyyyyy_0, g_zzzzzz_0_xxyyyyyz_0, g_zzzzzz_0_xxyyyyzz_0, g_zzzzzz_0_xxyyyzzz_0, g_zzzzzz_0_xxyyzzzz_0, g_zzzzzz_0_xxyzzzzz_0, g_zzzzzz_0_xxzzzzzz_0, g_zzzzzz_0_xyyyyyyy_0, g_zzzzzz_0_xyyyyyyz_0, g_zzzzzz_0_xyyyyyzz_0, g_zzzzzz_0_xyyyyzzz_0, g_zzzzzz_0_xyyyzzzz_0, g_zzzzzz_0_xyyzzzzz_0, g_zzzzzz_0_xyzzzzzz_0, g_zzzzzz_0_xzzzzzzz_0, g_zzzzzz_0_yyyyyyyy_0, g_zzzzzz_0_yyyyyyyz_0, g_zzzzzz_0_yyyyyyzz_0, g_zzzzzz_0_yyyyyzzz_0, g_zzzzzz_0_yyyyzzzz_0, g_zzzzzz_0_yyyzzzzz_0, g_zzzzzz_0_yyzzzzzz_0, g_zzzzzz_0_yzzzzzzz_0, g_zzzzzz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xxxxxxxx_0[i] = 5.0 * g_zzzz_0_xxxxxxxx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxxy_0[i] = 5.0 * g_zzzz_0_xxxxxxxy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxxy_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxxz_0[i] = 5.0 * g_zzzz_0_xxxxxxxz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxxz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxyy_0[i] = 5.0 * g_zzzz_0_xxxxxxyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxyy_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxyz_0[i] = 5.0 * g_zzzz_0_xxxxxxyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxzz_0[i] = 5.0 * g_zzzz_0_xxxxxxzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxyyy_0[i] = 5.0 * g_zzzz_0_xxxxxyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxyyy_1[i] * fz_be_0 + g_zzzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxyyz_0[i] = 5.0 * g_zzzz_0_xxxxxyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxyyz_1[i] * fz_be_0 + g_zzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxyzz_0[i] = 5.0 * g_zzzz_0_xxxxxyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxzzz_0[i] = 5.0 * g_zzzz_0_xxxxxzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyyyy_0[i] = 5.0 * g_zzzz_0_xxxxyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyyyy_1[i] * fz_be_0 + g_zzzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyyyz_0[i] = 5.0 * g_zzzz_0_xxxxyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyyyz_1[i] * fz_be_0 + g_zzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyyzz_0[i] = 5.0 * g_zzzz_0_xxxxyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyzzz_0[i] = 5.0 * g_zzzz_0_xxxxyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxzzzz_0[i] = 5.0 * g_zzzz_0_xxxxzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyyyy_0[i] = 5.0 * g_zzzz_0_xxxyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyyyy_1[i] * fz_be_0 + g_zzzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyyyz_0[i] = 5.0 * g_zzzz_0_xxxyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyyyz_1[i] * fz_be_0 + g_zzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyyzz_0[i] = 5.0 * g_zzzz_0_xxxyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyzzz_0[i] = 5.0 * g_zzzz_0_xxxyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyzzzz_0[i] = 5.0 * g_zzzz_0_xxxyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxzzzzz_0[i] = 5.0 * g_zzzz_0_xxxzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyyyy_0[i] = 5.0 * g_zzzz_0_xxyyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyyyy_1[i] * fz_be_0 + g_zzzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyyyz_0[i] = 5.0 * g_zzzz_0_xxyyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyyyz_1[i] * fz_be_0 + g_zzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyyzz_0[i] = 5.0 * g_zzzz_0_xxyyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyzzz_0[i] = 5.0 * g_zzzz_0_xxyyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyzzzz_0[i] = 5.0 * g_zzzz_0_xxyyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyzzzzz_0[i] = 5.0 * g_zzzz_0_xxyzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxzzzzzz_0[i] = 5.0 * g_zzzz_0_xxzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyyyy_0[i] = 5.0 * g_zzzz_0_xyyyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyyyy_1[i] * fz_be_0 + g_zzzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyyyz_0[i] = 5.0 * g_zzzz_0_xyyyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_zzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyyzz_0[i] = 5.0 * g_zzzz_0_xyyyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyzzz_0[i] = 5.0 * g_zzzz_0_xyyyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyzzzz_0[i] = 5.0 * g_zzzz_0_xyyyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyzzzzz_0[i] = 5.0 * g_zzzz_0_xyyzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyzzzzzz_0[i] = 5.0 * g_zzzz_0_xyzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xzzzzzzz_0[i] = 5.0 * g_zzzz_0_xzzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyyyy_0[i] = 5.0 * g_zzzz_0_yyyyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_zzzzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyyyz_0[i] = 5.0 * g_zzzz_0_yyyyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_zzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyyzz_0[i] = 5.0 * g_zzzz_0_yyyyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyzzz_0[i] = 5.0 * g_zzzz_0_yyyyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyzzzz_0[i] = 5.0 * g_zzzz_0_yyyyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyzzzzz_0[i] = 5.0 * g_zzzz_0_yyyzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyzzzzzz_0[i] = 5.0 * g_zzzz_0_yyzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yzzzzzzz_0[i] = 5.0 * g_zzzz_0_yzzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_zzzzzzzz_0[i] = 5.0 * g_zzzz_0_zzzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

