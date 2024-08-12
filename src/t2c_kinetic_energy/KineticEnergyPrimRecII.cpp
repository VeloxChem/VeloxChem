#include "KineticEnergyPrimRecII.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ii(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ii,
                            const size_t              idx_ovl_gi,
                            const size_t              idx_kin_gi,
                            const size_t              idx_kin_hh,
                            const size_t              idx_kin_hi,
                            const size_t              idx_ovl_ii,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GI

    auto ts_xxxx_xxxxxx = pbuffer.data(idx_ovl_gi);

    auto ts_xxxx_xxxxxy = pbuffer.data(idx_ovl_gi + 1);

    auto ts_xxxx_xxxxxz = pbuffer.data(idx_ovl_gi + 2);

    auto ts_xxxx_xxxxyy = pbuffer.data(idx_ovl_gi + 3);

    auto ts_xxxx_xxxxyz = pbuffer.data(idx_ovl_gi + 4);

    auto ts_xxxx_xxxxzz = pbuffer.data(idx_ovl_gi + 5);

    auto ts_xxxx_xxxyyy = pbuffer.data(idx_ovl_gi + 6);

    auto ts_xxxx_xxxyyz = pbuffer.data(idx_ovl_gi + 7);

    auto ts_xxxx_xxxyzz = pbuffer.data(idx_ovl_gi + 8);

    auto ts_xxxx_xxxzzz = pbuffer.data(idx_ovl_gi + 9);

    auto ts_xxxx_xxyyyy = pbuffer.data(idx_ovl_gi + 10);

    auto ts_xxxx_xxyyyz = pbuffer.data(idx_ovl_gi + 11);

    auto ts_xxxx_xxyyzz = pbuffer.data(idx_ovl_gi + 12);

    auto ts_xxxx_xxyzzz = pbuffer.data(idx_ovl_gi + 13);

    auto ts_xxxx_xxzzzz = pbuffer.data(idx_ovl_gi + 14);

    auto ts_xxxx_xyyyyy = pbuffer.data(idx_ovl_gi + 15);

    auto ts_xxxx_xyyyyz = pbuffer.data(idx_ovl_gi + 16);

    auto ts_xxxx_xyyyzz = pbuffer.data(idx_ovl_gi + 17);

    auto ts_xxxx_xyyzzz = pbuffer.data(idx_ovl_gi + 18);

    auto ts_xxxx_xyzzzz = pbuffer.data(idx_ovl_gi + 19);

    auto ts_xxxx_xzzzzz = pbuffer.data(idx_ovl_gi + 20);

    auto ts_xxxx_yyyyyy = pbuffer.data(idx_ovl_gi + 21);

    auto ts_xxxx_yyyyyz = pbuffer.data(idx_ovl_gi + 22);

    auto ts_xxxx_yyyyzz = pbuffer.data(idx_ovl_gi + 23);

    auto ts_xxxx_yyyzzz = pbuffer.data(idx_ovl_gi + 24);

    auto ts_xxxx_yyzzzz = pbuffer.data(idx_ovl_gi + 25);

    auto ts_xxxx_yzzzzz = pbuffer.data(idx_ovl_gi + 26);

    auto ts_xxxx_zzzzzz = pbuffer.data(idx_ovl_gi + 27);

    auto ts_xxxy_xxxxxx = pbuffer.data(idx_ovl_gi + 28);

    auto ts_xxxy_xxxxxz = pbuffer.data(idx_ovl_gi + 30);

    auto ts_xxxy_xxxxzz = pbuffer.data(idx_ovl_gi + 33);

    auto ts_xxxy_xxxzzz = pbuffer.data(idx_ovl_gi + 37);

    auto ts_xxxy_xxzzzz = pbuffer.data(idx_ovl_gi + 42);

    auto ts_xxxy_xzzzzz = pbuffer.data(idx_ovl_gi + 48);

    auto ts_xxxz_xxxxxx = pbuffer.data(idx_ovl_gi + 56);

    auto ts_xxxz_xxxxxy = pbuffer.data(idx_ovl_gi + 57);

    auto ts_xxxz_xxxxyy = pbuffer.data(idx_ovl_gi + 59);

    auto ts_xxxz_xxxyyy = pbuffer.data(idx_ovl_gi + 62);

    auto ts_xxxz_xxyyyy = pbuffer.data(idx_ovl_gi + 66);

    auto ts_xxxz_xyyyyy = pbuffer.data(idx_ovl_gi + 71);

    auto ts_xxyy_xxxxxx = pbuffer.data(idx_ovl_gi + 84);

    auto ts_xxyy_xxxxxy = pbuffer.data(idx_ovl_gi + 85);

    auto ts_xxyy_xxxxxz = pbuffer.data(idx_ovl_gi + 86);

    auto ts_xxyy_xxxxyy = pbuffer.data(idx_ovl_gi + 87);

    auto ts_xxyy_xxxxyz = pbuffer.data(idx_ovl_gi + 88);

    auto ts_xxyy_xxxxzz = pbuffer.data(idx_ovl_gi + 89);

    auto ts_xxyy_xxxyyy = pbuffer.data(idx_ovl_gi + 90);

    auto ts_xxyy_xxxyyz = pbuffer.data(idx_ovl_gi + 91);

    auto ts_xxyy_xxxyzz = pbuffer.data(idx_ovl_gi + 92);

    auto ts_xxyy_xxxzzz = pbuffer.data(idx_ovl_gi + 93);

    auto ts_xxyy_xxyyyy = pbuffer.data(idx_ovl_gi + 94);

    auto ts_xxyy_xxyyyz = pbuffer.data(idx_ovl_gi + 95);

    auto ts_xxyy_xxyyzz = pbuffer.data(idx_ovl_gi + 96);

    auto ts_xxyy_xxyzzz = pbuffer.data(idx_ovl_gi + 97);

    auto ts_xxyy_xxzzzz = pbuffer.data(idx_ovl_gi + 98);

    auto ts_xxyy_xyyyyy = pbuffer.data(idx_ovl_gi + 99);

    auto ts_xxyy_xyyyyz = pbuffer.data(idx_ovl_gi + 100);

    auto ts_xxyy_xyyyzz = pbuffer.data(idx_ovl_gi + 101);

    auto ts_xxyy_xyyzzz = pbuffer.data(idx_ovl_gi + 102);

    auto ts_xxyy_xyzzzz = pbuffer.data(idx_ovl_gi + 103);

    auto ts_xxyy_xzzzzz = pbuffer.data(idx_ovl_gi + 104);

    auto ts_xxyy_yyyyyy = pbuffer.data(idx_ovl_gi + 105);

    auto ts_xxyy_yyyyyz = pbuffer.data(idx_ovl_gi + 106);

    auto ts_xxyy_yyyyzz = pbuffer.data(idx_ovl_gi + 107);

    auto ts_xxyy_yyyzzz = pbuffer.data(idx_ovl_gi + 108);

    auto ts_xxyy_yyzzzz = pbuffer.data(idx_ovl_gi + 109);

    auto ts_xxyy_yzzzzz = pbuffer.data(idx_ovl_gi + 110);

    auto ts_xxyy_zzzzzz = pbuffer.data(idx_ovl_gi + 111);

    auto ts_xxzz_xxxxxx = pbuffer.data(idx_ovl_gi + 140);

    auto ts_xxzz_xxxxxy = pbuffer.data(idx_ovl_gi + 141);

    auto ts_xxzz_xxxxxz = pbuffer.data(idx_ovl_gi + 142);

    auto ts_xxzz_xxxxyy = pbuffer.data(idx_ovl_gi + 143);

    auto ts_xxzz_xxxxyz = pbuffer.data(idx_ovl_gi + 144);

    auto ts_xxzz_xxxxzz = pbuffer.data(idx_ovl_gi + 145);

    auto ts_xxzz_xxxyyy = pbuffer.data(idx_ovl_gi + 146);

    auto ts_xxzz_xxxyyz = pbuffer.data(idx_ovl_gi + 147);

    auto ts_xxzz_xxxyzz = pbuffer.data(idx_ovl_gi + 148);

    auto ts_xxzz_xxxzzz = pbuffer.data(idx_ovl_gi + 149);

    auto ts_xxzz_xxyyyy = pbuffer.data(idx_ovl_gi + 150);

    auto ts_xxzz_xxyyyz = pbuffer.data(idx_ovl_gi + 151);

    auto ts_xxzz_xxyyzz = pbuffer.data(idx_ovl_gi + 152);

    auto ts_xxzz_xxyzzz = pbuffer.data(idx_ovl_gi + 153);

    auto ts_xxzz_xxzzzz = pbuffer.data(idx_ovl_gi + 154);

    auto ts_xxzz_xyyyyy = pbuffer.data(idx_ovl_gi + 155);

    auto ts_xxzz_xyyyyz = pbuffer.data(idx_ovl_gi + 156);

    auto ts_xxzz_xyyyzz = pbuffer.data(idx_ovl_gi + 157);

    auto ts_xxzz_xyyzzz = pbuffer.data(idx_ovl_gi + 158);

    auto ts_xxzz_xyzzzz = pbuffer.data(idx_ovl_gi + 159);

    auto ts_xxzz_xzzzzz = pbuffer.data(idx_ovl_gi + 160);

    auto ts_xxzz_yyyyyy = pbuffer.data(idx_ovl_gi + 161);

    auto ts_xxzz_yyyyyz = pbuffer.data(idx_ovl_gi + 162);

    auto ts_xxzz_yyyyzz = pbuffer.data(idx_ovl_gi + 163);

    auto ts_xxzz_yyyzzz = pbuffer.data(idx_ovl_gi + 164);

    auto ts_xxzz_yyzzzz = pbuffer.data(idx_ovl_gi + 165);

    auto ts_xxzz_yzzzzz = pbuffer.data(idx_ovl_gi + 166);

    auto ts_xxzz_zzzzzz = pbuffer.data(idx_ovl_gi + 167);

    auto ts_xyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 169);

    auto ts_xyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 171);

    auto ts_xyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 172);

    auto ts_xyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 174);

    auto ts_xyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 175);

    auto ts_xyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 176);

    auto ts_xyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 178);

    auto ts_xyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 179);

    auto ts_xyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 180);

    auto ts_xyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 181);

    auto ts_xyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 183);

    auto ts_xyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 184);

    auto ts_xyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 185);

    auto ts_xyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 186);

    auto ts_xyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 187);

    auto ts_xyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 189);

    auto ts_xyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 190);

    auto ts_xyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 191);

    auto ts_xyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 192);

    auto ts_xyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 193);

    auto ts_xyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 194);

    auto ts_xyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 195);

    auto ts_xzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 254);

    auto ts_xzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 256);

    auto ts_xzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 257);

    auto ts_xzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 259);

    auto ts_xzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 260);

    auto ts_xzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 261);

    auto ts_xzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 263);

    auto ts_xzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 264);

    auto ts_xzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 265);

    auto ts_xzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 266);

    auto ts_xzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 268);

    auto ts_xzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 269);

    auto ts_xzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 270);

    auto ts_xzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 271);

    auto ts_xzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 272);

    auto ts_xzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 273);

    auto ts_xzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 274);

    auto ts_xzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 275);

    auto ts_xzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 276);

    auto ts_xzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 277);

    auto ts_xzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 278);

    auto ts_xzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 279);

    auto ts_yyyy_xxxxxx = pbuffer.data(idx_ovl_gi + 280);

    auto ts_yyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 281);

    auto ts_yyyy_xxxxxz = pbuffer.data(idx_ovl_gi + 282);

    auto ts_yyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 283);

    auto ts_yyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 284);

    auto ts_yyyy_xxxxzz = pbuffer.data(idx_ovl_gi + 285);

    auto ts_yyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 286);

    auto ts_yyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 287);

    auto ts_yyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 288);

    auto ts_yyyy_xxxzzz = pbuffer.data(idx_ovl_gi + 289);

    auto ts_yyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 290);

    auto ts_yyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 291);

    auto ts_yyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 292);

    auto ts_yyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 293);

    auto ts_yyyy_xxzzzz = pbuffer.data(idx_ovl_gi + 294);

    auto ts_yyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 295);

    auto ts_yyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 296);

    auto ts_yyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 297);

    auto ts_yyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 298);

    auto ts_yyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 299);

    auto ts_yyyy_xzzzzz = pbuffer.data(idx_ovl_gi + 300);

    auto ts_yyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 301);

    auto ts_yyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 302);

    auto ts_yyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 303);

    auto ts_yyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 304);

    auto ts_yyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 305);

    auto ts_yyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 306);

    auto ts_yyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 307);

    auto ts_yyyz_xxxxxy = pbuffer.data(idx_ovl_gi + 309);

    auto ts_yyyz_xxxxyy = pbuffer.data(idx_ovl_gi + 311);

    auto ts_yyyz_xxxyyy = pbuffer.data(idx_ovl_gi + 314);

    auto ts_yyyz_xxyyyy = pbuffer.data(idx_ovl_gi + 318);

    auto ts_yyyz_xyyyyy = pbuffer.data(idx_ovl_gi + 323);

    auto ts_yyyz_yyyyyy = pbuffer.data(idx_ovl_gi + 329);

    auto ts_yyzz_xxxxxx = pbuffer.data(idx_ovl_gi + 336);

    auto ts_yyzz_xxxxxy = pbuffer.data(idx_ovl_gi + 337);

    auto ts_yyzz_xxxxxz = pbuffer.data(idx_ovl_gi + 338);

    auto ts_yyzz_xxxxyy = pbuffer.data(idx_ovl_gi + 339);

    auto ts_yyzz_xxxxyz = pbuffer.data(idx_ovl_gi + 340);

    auto ts_yyzz_xxxxzz = pbuffer.data(idx_ovl_gi + 341);

    auto ts_yyzz_xxxyyy = pbuffer.data(idx_ovl_gi + 342);

    auto ts_yyzz_xxxyyz = pbuffer.data(idx_ovl_gi + 343);

    auto ts_yyzz_xxxyzz = pbuffer.data(idx_ovl_gi + 344);

    auto ts_yyzz_xxxzzz = pbuffer.data(idx_ovl_gi + 345);

    auto ts_yyzz_xxyyyy = pbuffer.data(idx_ovl_gi + 346);

    auto ts_yyzz_xxyyyz = pbuffer.data(idx_ovl_gi + 347);

    auto ts_yyzz_xxyyzz = pbuffer.data(idx_ovl_gi + 348);

    auto ts_yyzz_xxyzzz = pbuffer.data(idx_ovl_gi + 349);

    auto ts_yyzz_xxzzzz = pbuffer.data(idx_ovl_gi + 350);

    auto ts_yyzz_xyyyyy = pbuffer.data(idx_ovl_gi + 351);

    auto ts_yyzz_xyyyyz = pbuffer.data(idx_ovl_gi + 352);

    auto ts_yyzz_xyyyzz = pbuffer.data(idx_ovl_gi + 353);

    auto ts_yyzz_xyyzzz = pbuffer.data(idx_ovl_gi + 354);

    auto ts_yyzz_xyzzzz = pbuffer.data(idx_ovl_gi + 355);

    auto ts_yyzz_xzzzzz = pbuffer.data(idx_ovl_gi + 356);

    auto ts_yyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 357);

    auto ts_yyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 358);

    auto ts_yyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 359);

    auto ts_yyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 360);

    auto ts_yyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 361);

    auto ts_yyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 362);

    auto ts_yyzz_zzzzzz = pbuffer.data(idx_ovl_gi + 363);

    auto ts_yzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 364);

    auto ts_yzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 366);

    auto ts_yzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 368);

    auto ts_yzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 369);

    auto ts_yzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 371);

    auto ts_yzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 372);

    auto ts_yzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 373);

    auto ts_yzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 375);

    auto ts_yzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 376);

    auto ts_yzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 377);

    auto ts_yzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 378);

    auto ts_yzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 380);

    auto ts_yzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 381);

    auto ts_yzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 382);

    auto ts_yzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 383);

    auto ts_yzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 384);

    auto ts_yzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 386);

    auto ts_yzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 387);

    auto ts_yzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 388);

    auto ts_yzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 389);

    auto ts_yzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 390);

    auto ts_yzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 391);

    auto ts_zzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 392);

    auto ts_zzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 393);

    auto ts_zzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 394);

    auto ts_zzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 395);

    auto ts_zzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 396);

    auto ts_zzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 397);

    auto ts_zzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 398);

    auto ts_zzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 399);

    auto ts_zzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 400);

    auto ts_zzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 401);

    auto ts_zzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 402);

    auto ts_zzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 403);

    auto ts_zzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 404);

    auto ts_zzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 405);

    auto ts_zzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 406);

    auto ts_zzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 407);

    auto ts_zzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 408);

    auto ts_zzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 409);

    auto ts_zzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 410);

    auto ts_zzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 411);

    auto ts_zzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 412);

    auto ts_zzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 413);

    auto ts_zzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 414);

    auto ts_zzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 415);

    auto ts_zzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 416);

    auto ts_zzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 417);

    auto ts_zzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 418);

    auto ts_zzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 419);

    // Set up components of auxiliary buffer : GI

    auto tk_xxxx_xxxxxx = pbuffer.data(idx_kin_gi);

    auto tk_xxxx_xxxxxy = pbuffer.data(idx_kin_gi + 1);

    auto tk_xxxx_xxxxxz = pbuffer.data(idx_kin_gi + 2);

    auto tk_xxxx_xxxxyy = pbuffer.data(idx_kin_gi + 3);

    auto tk_xxxx_xxxxyz = pbuffer.data(idx_kin_gi + 4);

    auto tk_xxxx_xxxxzz = pbuffer.data(idx_kin_gi + 5);

    auto tk_xxxx_xxxyyy = pbuffer.data(idx_kin_gi + 6);

    auto tk_xxxx_xxxyyz = pbuffer.data(idx_kin_gi + 7);

    auto tk_xxxx_xxxyzz = pbuffer.data(idx_kin_gi + 8);

    auto tk_xxxx_xxxzzz = pbuffer.data(idx_kin_gi + 9);

    auto tk_xxxx_xxyyyy = pbuffer.data(idx_kin_gi + 10);

    auto tk_xxxx_xxyyyz = pbuffer.data(idx_kin_gi + 11);

    auto tk_xxxx_xxyyzz = pbuffer.data(idx_kin_gi + 12);

    auto tk_xxxx_xxyzzz = pbuffer.data(idx_kin_gi + 13);

    auto tk_xxxx_xxzzzz = pbuffer.data(idx_kin_gi + 14);

    auto tk_xxxx_xyyyyy = pbuffer.data(idx_kin_gi + 15);

    auto tk_xxxx_xyyyyz = pbuffer.data(idx_kin_gi + 16);

    auto tk_xxxx_xyyyzz = pbuffer.data(idx_kin_gi + 17);

    auto tk_xxxx_xyyzzz = pbuffer.data(idx_kin_gi + 18);

    auto tk_xxxx_xyzzzz = pbuffer.data(idx_kin_gi + 19);

    auto tk_xxxx_xzzzzz = pbuffer.data(idx_kin_gi + 20);

    auto tk_xxxx_yyyyyy = pbuffer.data(idx_kin_gi + 21);

    auto tk_xxxx_yyyyyz = pbuffer.data(idx_kin_gi + 22);

    auto tk_xxxx_yyyyzz = pbuffer.data(idx_kin_gi + 23);

    auto tk_xxxx_yyyzzz = pbuffer.data(idx_kin_gi + 24);

    auto tk_xxxx_yyzzzz = pbuffer.data(idx_kin_gi + 25);

    auto tk_xxxx_yzzzzz = pbuffer.data(idx_kin_gi + 26);

    auto tk_xxxx_zzzzzz = pbuffer.data(idx_kin_gi + 27);

    auto tk_xxxy_xxxxxx = pbuffer.data(idx_kin_gi + 28);

    auto tk_xxxy_xxxxxz = pbuffer.data(idx_kin_gi + 30);

    auto tk_xxxy_xxxxzz = pbuffer.data(idx_kin_gi + 33);

    auto tk_xxxy_xxxzzz = pbuffer.data(idx_kin_gi + 37);

    auto tk_xxxy_xxzzzz = pbuffer.data(idx_kin_gi + 42);

    auto tk_xxxy_xzzzzz = pbuffer.data(idx_kin_gi + 48);

    auto tk_xxxz_xxxxxx = pbuffer.data(idx_kin_gi + 56);

    auto tk_xxxz_xxxxxy = pbuffer.data(idx_kin_gi + 57);

    auto tk_xxxz_xxxxyy = pbuffer.data(idx_kin_gi + 59);

    auto tk_xxxz_xxxyyy = pbuffer.data(idx_kin_gi + 62);

    auto tk_xxxz_xxyyyy = pbuffer.data(idx_kin_gi + 66);

    auto tk_xxxz_xyyyyy = pbuffer.data(idx_kin_gi + 71);

    auto tk_xxyy_xxxxxx = pbuffer.data(idx_kin_gi + 84);

    auto tk_xxyy_xxxxxy = pbuffer.data(idx_kin_gi + 85);

    auto tk_xxyy_xxxxxz = pbuffer.data(idx_kin_gi + 86);

    auto tk_xxyy_xxxxyy = pbuffer.data(idx_kin_gi + 87);

    auto tk_xxyy_xxxxyz = pbuffer.data(idx_kin_gi + 88);

    auto tk_xxyy_xxxxzz = pbuffer.data(idx_kin_gi + 89);

    auto tk_xxyy_xxxyyy = pbuffer.data(idx_kin_gi + 90);

    auto tk_xxyy_xxxyyz = pbuffer.data(idx_kin_gi + 91);

    auto tk_xxyy_xxxyzz = pbuffer.data(idx_kin_gi + 92);

    auto tk_xxyy_xxxzzz = pbuffer.data(idx_kin_gi + 93);

    auto tk_xxyy_xxyyyy = pbuffer.data(idx_kin_gi + 94);

    auto tk_xxyy_xxyyyz = pbuffer.data(idx_kin_gi + 95);

    auto tk_xxyy_xxyyzz = pbuffer.data(idx_kin_gi + 96);

    auto tk_xxyy_xxyzzz = pbuffer.data(idx_kin_gi + 97);

    auto tk_xxyy_xxzzzz = pbuffer.data(idx_kin_gi + 98);

    auto tk_xxyy_xyyyyy = pbuffer.data(idx_kin_gi + 99);

    auto tk_xxyy_xyyyyz = pbuffer.data(idx_kin_gi + 100);

    auto tk_xxyy_xyyyzz = pbuffer.data(idx_kin_gi + 101);

    auto tk_xxyy_xyyzzz = pbuffer.data(idx_kin_gi + 102);

    auto tk_xxyy_xyzzzz = pbuffer.data(idx_kin_gi + 103);

    auto tk_xxyy_xzzzzz = pbuffer.data(idx_kin_gi + 104);

    auto tk_xxyy_yyyyyy = pbuffer.data(idx_kin_gi + 105);

    auto tk_xxyy_yyyyyz = pbuffer.data(idx_kin_gi + 106);

    auto tk_xxyy_yyyyzz = pbuffer.data(idx_kin_gi + 107);

    auto tk_xxyy_yyyzzz = pbuffer.data(idx_kin_gi + 108);

    auto tk_xxyy_yyzzzz = pbuffer.data(idx_kin_gi + 109);

    auto tk_xxyy_yzzzzz = pbuffer.data(idx_kin_gi + 110);

    auto tk_xxyy_zzzzzz = pbuffer.data(idx_kin_gi + 111);

    auto tk_xxzz_xxxxxx = pbuffer.data(idx_kin_gi + 140);

    auto tk_xxzz_xxxxxy = pbuffer.data(idx_kin_gi + 141);

    auto tk_xxzz_xxxxxz = pbuffer.data(idx_kin_gi + 142);

    auto tk_xxzz_xxxxyy = pbuffer.data(idx_kin_gi + 143);

    auto tk_xxzz_xxxxyz = pbuffer.data(idx_kin_gi + 144);

    auto tk_xxzz_xxxxzz = pbuffer.data(idx_kin_gi + 145);

    auto tk_xxzz_xxxyyy = pbuffer.data(idx_kin_gi + 146);

    auto tk_xxzz_xxxyyz = pbuffer.data(idx_kin_gi + 147);

    auto tk_xxzz_xxxyzz = pbuffer.data(idx_kin_gi + 148);

    auto tk_xxzz_xxxzzz = pbuffer.data(idx_kin_gi + 149);

    auto tk_xxzz_xxyyyy = pbuffer.data(idx_kin_gi + 150);

    auto tk_xxzz_xxyyyz = pbuffer.data(idx_kin_gi + 151);

    auto tk_xxzz_xxyyzz = pbuffer.data(idx_kin_gi + 152);

    auto tk_xxzz_xxyzzz = pbuffer.data(idx_kin_gi + 153);

    auto tk_xxzz_xxzzzz = pbuffer.data(idx_kin_gi + 154);

    auto tk_xxzz_xyyyyy = pbuffer.data(idx_kin_gi + 155);

    auto tk_xxzz_xyyyyz = pbuffer.data(idx_kin_gi + 156);

    auto tk_xxzz_xyyyzz = pbuffer.data(idx_kin_gi + 157);

    auto tk_xxzz_xyyzzz = pbuffer.data(idx_kin_gi + 158);

    auto tk_xxzz_xyzzzz = pbuffer.data(idx_kin_gi + 159);

    auto tk_xxzz_xzzzzz = pbuffer.data(idx_kin_gi + 160);

    auto tk_xxzz_yyyyyy = pbuffer.data(idx_kin_gi + 161);

    auto tk_xxzz_yyyyyz = pbuffer.data(idx_kin_gi + 162);

    auto tk_xxzz_yyyyzz = pbuffer.data(idx_kin_gi + 163);

    auto tk_xxzz_yyyzzz = pbuffer.data(idx_kin_gi + 164);

    auto tk_xxzz_yyzzzz = pbuffer.data(idx_kin_gi + 165);

    auto tk_xxzz_yzzzzz = pbuffer.data(idx_kin_gi + 166);

    auto tk_xxzz_zzzzzz = pbuffer.data(idx_kin_gi + 167);

    auto tk_xyyy_xxxxxy = pbuffer.data(idx_kin_gi + 169);

    auto tk_xyyy_xxxxyy = pbuffer.data(idx_kin_gi + 171);

    auto tk_xyyy_xxxxyz = pbuffer.data(idx_kin_gi + 172);

    auto tk_xyyy_xxxyyy = pbuffer.data(idx_kin_gi + 174);

    auto tk_xyyy_xxxyyz = pbuffer.data(idx_kin_gi + 175);

    auto tk_xyyy_xxxyzz = pbuffer.data(idx_kin_gi + 176);

    auto tk_xyyy_xxyyyy = pbuffer.data(idx_kin_gi + 178);

    auto tk_xyyy_xxyyyz = pbuffer.data(idx_kin_gi + 179);

    auto tk_xyyy_xxyyzz = pbuffer.data(idx_kin_gi + 180);

    auto tk_xyyy_xxyzzz = pbuffer.data(idx_kin_gi + 181);

    auto tk_xyyy_xyyyyy = pbuffer.data(idx_kin_gi + 183);

    auto tk_xyyy_xyyyyz = pbuffer.data(idx_kin_gi + 184);

    auto tk_xyyy_xyyyzz = pbuffer.data(idx_kin_gi + 185);

    auto tk_xyyy_xyyzzz = pbuffer.data(idx_kin_gi + 186);

    auto tk_xyyy_xyzzzz = pbuffer.data(idx_kin_gi + 187);

    auto tk_xyyy_yyyyyy = pbuffer.data(idx_kin_gi + 189);

    auto tk_xyyy_yyyyyz = pbuffer.data(idx_kin_gi + 190);

    auto tk_xyyy_yyyyzz = pbuffer.data(idx_kin_gi + 191);

    auto tk_xyyy_yyyzzz = pbuffer.data(idx_kin_gi + 192);

    auto tk_xyyy_yyzzzz = pbuffer.data(idx_kin_gi + 193);

    auto tk_xyyy_yzzzzz = pbuffer.data(idx_kin_gi + 194);

    auto tk_xyyy_zzzzzz = pbuffer.data(idx_kin_gi + 195);

    auto tk_xzzz_xxxxxz = pbuffer.data(idx_kin_gi + 254);

    auto tk_xzzz_xxxxyz = pbuffer.data(idx_kin_gi + 256);

    auto tk_xzzz_xxxxzz = pbuffer.data(idx_kin_gi + 257);

    auto tk_xzzz_xxxyyz = pbuffer.data(idx_kin_gi + 259);

    auto tk_xzzz_xxxyzz = pbuffer.data(idx_kin_gi + 260);

    auto tk_xzzz_xxxzzz = pbuffer.data(idx_kin_gi + 261);

    auto tk_xzzz_xxyyyz = pbuffer.data(idx_kin_gi + 263);

    auto tk_xzzz_xxyyzz = pbuffer.data(idx_kin_gi + 264);

    auto tk_xzzz_xxyzzz = pbuffer.data(idx_kin_gi + 265);

    auto tk_xzzz_xxzzzz = pbuffer.data(idx_kin_gi + 266);

    auto tk_xzzz_xyyyyz = pbuffer.data(idx_kin_gi + 268);

    auto tk_xzzz_xyyyzz = pbuffer.data(idx_kin_gi + 269);

    auto tk_xzzz_xyyzzz = pbuffer.data(idx_kin_gi + 270);

    auto tk_xzzz_xyzzzz = pbuffer.data(idx_kin_gi + 271);

    auto tk_xzzz_xzzzzz = pbuffer.data(idx_kin_gi + 272);

    auto tk_xzzz_yyyyyy = pbuffer.data(idx_kin_gi + 273);

    auto tk_xzzz_yyyyyz = pbuffer.data(idx_kin_gi + 274);

    auto tk_xzzz_yyyyzz = pbuffer.data(idx_kin_gi + 275);

    auto tk_xzzz_yyyzzz = pbuffer.data(idx_kin_gi + 276);

    auto tk_xzzz_yyzzzz = pbuffer.data(idx_kin_gi + 277);

    auto tk_xzzz_yzzzzz = pbuffer.data(idx_kin_gi + 278);

    auto tk_xzzz_zzzzzz = pbuffer.data(idx_kin_gi + 279);

    auto tk_yyyy_xxxxxx = pbuffer.data(idx_kin_gi + 280);

    auto tk_yyyy_xxxxxy = pbuffer.data(idx_kin_gi + 281);

    auto tk_yyyy_xxxxxz = pbuffer.data(idx_kin_gi + 282);

    auto tk_yyyy_xxxxyy = pbuffer.data(idx_kin_gi + 283);

    auto tk_yyyy_xxxxyz = pbuffer.data(idx_kin_gi + 284);

    auto tk_yyyy_xxxxzz = pbuffer.data(idx_kin_gi + 285);

    auto tk_yyyy_xxxyyy = pbuffer.data(idx_kin_gi + 286);

    auto tk_yyyy_xxxyyz = pbuffer.data(idx_kin_gi + 287);

    auto tk_yyyy_xxxyzz = pbuffer.data(idx_kin_gi + 288);

    auto tk_yyyy_xxxzzz = pbuffer.data(idx_kin_gi + 289);

    auto tk_yyyy_xxyyyy = pbuffer.data(idx_kin_gi + 290);

    auto tk_yyyy_xxyyyz = pbuffer.data(idx_kin_gi + 291);

    auto tk_yyyy_xxyyzz = pbuffer.data(idx_kin_gi + 292);

    auto tk_yyyy_xxyzzz = pbuffer.data(idx_kin_gi + 293);

    auto tk_yyyy_xxzzzz = pbuffer.data(idx_kin_gi + 294);

    auto tk_yyyy_xyyyyy = pbuffer.data(idx_kin_gi + 295);

    auto tk_yyyy_xyyyyz = pbuffer.data(idx_kin_gi + 296);

    auto tk_yyyy_xyyyzz = pbuffer.data(idx_kin_gi + 297);

    auto tk_yyyy_xyyzzz = pbuffer.data(idx_kin_gi + 298);

    auto tk_yyyy_xyzzzz = pbuffer.data(idx_kin_gi + 299);

    auto tk_yyyy_xzzzzz = pbuffer.data(idx_kin_gi + 300);

    auto tk_yyyy_yyyyyy = pbuffer.data(idx_kin_gi + 301);

    auto tk_yyyy_yyyyyz = pbuffer.data(idx_kin_gi + 302);

    auto tk_yyyy_yyyyzz = pbuffer.data(idx_kin_gi + 303);

    auto tk_yyyy_yyyzzz = pbuffer.data(idx_kin_gi + 304);

    auto tk_yyyy_yyzzzz = pbuffer.data(idx_kin_gi + 305);

    auto tk_yyyy_yzzzzz = pbuffer.data(idx_kin_gi + 306);

    auto tk_yyyy_zzzzzz = pbuffer.data(idx_kin_gi + 307);

    auto tk_yyyz_xxxxxy = pbuffer.data(idx_kin_gi + 309);

    auto tk_yyyz_xxxxyy = pbuffer.data(idx_kin_gi + 311);

    auto tk_yyyz_xxxyyy = pbuffer.data(idx_kin_gi + 314);

    auto tk_yyyz_xxyyyy = pbuffer.data(idx_kin_gi + 318);

    auto tk_yyyz_xyyyyy = pbuffer.data(idx_kin_gi + 323);

    auto tk_yyyz_yyyyyy = pbuffer.data(idx_kin_gi + 329);

    auto tk_yyzz_xxxxxx = pbuffer.data(idx_kin_gi + 336);

    auto tk_yyzz_xxxxxy = pbuffer.data(idx_kin_gi + 337);

    auto tk_yyzz_xxxxxz = pbuffer.data(idx_kin_gi + 338);

    auto tk_yyzz_xxxxyy = pbuffer.data(idx_kin_gi + 339);

    auto tk_yyzz_xxxxyz = pbuffer.data(idx_kin_gi + 340);

    auto tk_yyzz_xxxxzz = pbuffer.data(idx_kin_gi + 341);

    auto tk_yyzz_xxxyyy = pbuffer.data(idx_kin_gi + 342);

    auto tk_yyzz_xxxyyz = pbuffer.data(idx_kin_gi + 343);

    auto tk_yyzz_xxxyzz = pbuffer.data(idx_kin_gi + 344);

    auto tk_yyzz_xxxzzz = pbuffer.data(idx_kin_gi + 345);

    auto tk_yyzz_xxyyyy = pbuffer.data(idx_kin_gi + 346);

    auto tk_yyzz_xxyyyz = pbuffer.data(idx_kin_gi + 347);

    auto tk_yyzz_xxyyzz = pbuffer.data(idx_kin_gi + 348);

    auto tk_yyzz_xxyzzz = pbuffer.data(idx_kin_gi + 349);

    auto tk_yyzz_xxzzzz = pbuffer.data(idx_kin_gi + 350);

    auto tk_yyzz_xyyyyy = pbuffer.data(idx_kin_gi + 351);

    auto tk_yyzz_xyyyyz = pbuffer.data(idx_kin_gi + 352);

    auto tk_yyzz_xyyyzz = pbuffer.data(idx_kin_gi + 353);

    auto tk_yyzz_xyyzzz = pbuffer.data(idx_kin_gi + 354);

    auto tk_yyzz_xyzzzz = pbuffer.data(idx_kin_gi + 355);

    auto tk_yyzz_xzzzzz = pbuffer.data(idx_kin_gi + 356);

    auto tk_yyzz_yyyyyy = pbuffer.data(idx_kin_gi + 357);

    auto tk_yyzz_yyyyyz = pbuffer.data(idx_kin_gi + 358);

    auto tk_yyzz_yyyyzz = pbuffer.data(idx_kin_gi + 359);

    auto tk_yyzz_yyyzzz = pbuffer.data(idx_kin_gi + 360);

    auto tk_yyzz_yyzzzz = pbuffer.data(idx_kin_gi + 361);

    auto tk_yyzz_yzzzzz = pbuffer.data(idx_kin_gi + 362);

    auto tk_yyzz_zzzzzz = pbuffer.data(idx_kin_gi + 363);

    auto tk_yzzz_xxxxxx = pbuffer.data(idx_kin_gi + 364);

    auto tk_yzzz_xxxxxz = pbuffer.data(idx_kin_gi + 366);

    auto tk_yzzz_xxxxyz = pbuffer.data(idx_kin_gi + 368);

    auto tk_yzzz_xxxxzz = pbuffer.data(idx_kin_gi + 369);

    auto tk_yzzz_xxxyyz = pbuffer.data(idx_kin_gi + 371);

    auto tk_yzzz_xxxyzz = pbuffer.data(idx_kin_gi + 372);

    auto tk_yzzz_xxxzzz = pbuffer.data(idx_kin_gi + 373);

    auto tk_yzzz_xxyyyz = pbuffer.data(idx_kin_gi + 375);

    auto tk_yzzz_xxyyzz = pbuffer.data(idx_kin_gi + 376);

    auto tk_yzzz_xxyzzz = pbuffer.data(idx_kin_gi + 377);

    auto tk_yzzz_xxzzzz = pbuffer.data(idx_kin_gi + 378);

    auto tk_yzzz_xyyyyz = pbuffer.data(idx_kin_gi + 380);

    auto tk_yzzz_xyyyzz = pbuffer.data(idx_kin_gi + 381);

    auto tk_yzzz_xyyzzz = pbuffer.data(idx_kin_gi + 382);

    auto tk_yzzz_xyzzzz = pbuffer.data(idx_kin_gi + 383);

    auto tk_yzzz_xzzzzz = pbuffer.data(idx_kin_gi + 384);

    auto tk_yzzz_yyyyyz = pbuffer.data(idx_kin_gi + 386);

    auto tk_yzzz_yyyyzz = pbuffer.data(idx_kin_gi + 387);

    auto tk_yzzz_yyyzzz = pbuffer.data(idx_kin_gi + 388);

    auto tk_yzzz_yyzzzz = pbuffer.data(idx_kin_gi + 389);

    auto tk_yzzz_yzzzzz = pbuffer.data(idx_kin_gi + 390);

    auto tk_yzzz_zzzzzz = pbuffer.data(idx_kin_gi + 391);

    auto tk_zzzz_xxxxxx = pbuffer.data(idx_kin_gi + 392);

    auto tk_zzzz_xxxxxy = pbuffer.data(idx_kin_gi + 393);

    auto tk_zzzz_xxxxxz = pbuffer.data(idx_kin_gi + 394);

    auto tk_zzzz_xxxxyy = pbuffer.data(idx_kin_gi + 395);

    auto tk_zzzz_xxxxyz = pbuffer.data(idx_kin_gi + 396);

    auto tk_zzzz_xxxxzz = pbuffer.data(idx_kin_gi + 397);

    auto tk_zzzz_xxxyyy = pbuffer.data(idx_kin_gi + 398);

    auto tk_zzzz_xxxyyz = pbuffer.data(idx_kin_gi + 399);

    auto tk_zzzz_xxxyzz = pbuffer.data(idx_kin_gi + 400);

    auto tk_zzzz_xxxzzz = pbuffer.data(idx_kin_gi + 401);

    auto tk_zzzz_xxyyyy = pbuffer.data(idx_kin_gi + 402);

    auto tk_zzzz_xxyyyz = pbuffer.data(idx_kin_gi + 403);

    auto tk_zzzz_xxyyzz = pbuffer.data(idx_kin_gi + 404);

    auto tk_zzzz_xxyzzz = pbuffer.data(idx_kin_gi + 405);

    auto tk_zzzz_xxzzzz = pbuffer.data(idx_kin_gi + 406);

    auto tk_zzzz_xyyyyy = pbuffer.data(idx_kin_gi + 407);

    auto tk_zzzz_xyyyyz = pbuffer.data(idx_kin_gi + 408);

    auto tk_zzzz_xyyyzz = pbuffer.data(idx_kin_gi + 409);

    auto tk_zzzz_xyyzzz = pbuffer.data(idx_kin_gi + 410);

    auto tk_zzzz_xyzzzz = pbuffer.data(idx_kin_gi + 411);

    auto tk_zzzz_xzzzzz = pbuffer.data(idx_kin_gi + 412);

    auto tk_zzzz_yyyyyy = pbuffer.data(idx_kin_gi + 413);

    auto tk_zzzz_yyyyyz = pbuffer.data(idx_kin_gi + 414);

    auto tk_zzzz_yyyyzz = pbuffer.data(idx_kin_gi + 415);

    auto tk_zzzz_yyyzzz = pbuffer.data(idx_kin_gi + 416);

    auto tk_zzzz_yyzzzz = pbuffer.data(idx_kin_gi + 417);

    auto tk_zzzz_yzzzzz = pbuffer.data(idx_kin_gi + 418);

    auto tk_zzzz_zzzzzz = pbuffer.data(idx_kin_gi + 419);

    // Set up components of auxiliary buffer : HH

    auto tk_xxxxx_xxxxx = pbuffer.data(idx_kin_hh);

    auto tk_xxxxx_xxxxy = pbuffer.data(idx_kin_hh + 1);

    auto tk_xxxxx_xxxxz = pbuffer.data(idx_kin_hh + 2);

    auto tk_xxxxx_xxxyy = pbuffer.data(idx_kin_hh + 3);

    auto tk_xxxxx_xxxyz = pbuffer.data(idx_kin_hh + 4);

    auto tk_xxxxx_xxxzz = pbuffer.data(idx_kin_hh + 5);

    auto tk_xxxxx_xxyyy = pbuffer.data(idx_kin_hh + 6);

    auto tk_xxxxx_xxyyz = pbuffer.data(idx_kin_hh + 7);

    auto tk_xxxxx_xxyzz = pbuffer.data(idx_kin_hh + 8);

    auto tk_xxxxx_xxzzz = pbuffer.data(idx_kin_hh + 9);

    auto tk_xxxxx_xyyyy = pbuffer.data(idx_kin_hh + 10);

    auto tk_xxxxx_xyyyz = pbuffer.data(idx_kin_hh + 11);

    auto tk_xxxxx_xyyzz = pbuffer.data(idx_kin_hh + 12);

    auto tk_xxxxx_xyzzz = pbuffer.data(idx_kin_hh + 13);

    auto tk_xxxxx_xzzzz = pbuffer.data(idx_kin_hh + 14);

    auto tk_xxxxx_yyyyy = pbuffer.data(idx_kin_hh + 15);

    auto tk_xxxxx_yyyyz = pbuffer.data(idx_kin_hh + 16);

    auto tk_xxxxx_yyyzz = pbuffer.data(idx_kin_hh + 17);

    auto tk_xxxxx_yyzzz = pbuffer.data(idx_kin_hh + 18);

    auto tk_xxxxx_yzzzz = pbuffer.data(idx_kin_hh + 19);

    auto tk_xxxxx_zzzzz = pbuffer.data(idx_kin_hh + 20);

    auto tk_xxxxz_xxxxz = pbuffer.data(idx_kin_hh + 44);

    auto tk_xxxxz_xxxyz = pbuffer.data(idx_kin_hh + 46);

    auto tk_xxxxz_xxxzz = pbuffer.data(idx_kin_hh + 47);

    auto tk_xxxxz_xxyyz = pbuffer.data(idx_kin_hh + 49);

    auto tk_xxxxz_xxyzz = pbuffer.data(idx_kin_hh + 50);

    auto tk_xxxxz_xxzzz = pbuffer.data(idx_kin_hh + 51);

    auto tk_xxxxz_xyyyz = pbuffer.data(idx_kin_hh + 53);

    auto tk_xxxxz_xyyzz = pbuffer.data(idx_kin_hh + 54);

    auto tk_xxxxz_xyzzz = pbuffer.data(idx_kin_hh + 55);

    auto tk_xxxxz_xzzzz = pbuffer.data(idx_kin_hh + 56);

    auto tk_xxxxz_yyyyz = pbuffer.data(idx_kin_hh + 58);

    auto tk_xxxxz_yyyzz = pbuffer.data(idx_kin_hh + 59);

    auto tk_xxxxz_yyzzz = pbuffer.data(idx_kin_hh + 60);

    auto tk_xxxxz_yzzzz = pbuffer.data(idx_kin_hh + 61);

    auto tk_xxxxz_zzzzz = pbuffer.data(idx_kin_hh + 62);

    auto tk_xxxyy_xxxxx = pbuffer.data(idx_kin_hh + 63);

    auto tk_xxxyy_xxxxy = pbuffer.data(idx_kin_hh + 64);

    auto tk_xxxyy_xxxxz = pbuffer.data(idx_kin_hh + 65);

    auto tk_xxxyy_xxxyy = pbuffer.data(idx_kin_hh + 66);

    auto tk_xxxyy_xxxyz = pbuffer.data(idx_kin_hh + 67);

    auto tk_xxxyy_xxxzz = pbuffer.data(idx_kin_hh + 68);

    auto tk_xxxyy_xxyyy = pbuffer.data(idx_kin_hh + 69);

    auto tk_xxxyy_xxyyz = pbuffer.data(idx_kin_hh + 70);

    auto tk_xxxyy_xxyzz = pbuffer.data(idx_kin_hh + 71);

    auto tk_xxxyy_xxzzz = pbuffer.data(idx_kin_hh + 72);

    auto tk_xxxyy_xyyyy = pbuffer.data(idx_kin_hh + 73);

    auto tk_xxxyy_xyyyz = pbuffer.data(idx_kin_hh + 74);

    auto tk_xxxyy_xyyzz = pbuffer.data(idx_kin_hh + 75);

    auto tk_xxxyy_xyzzz = pbuffer.data(idx_kin_hh + 76);

    auto tk_xxxyy_xzzzz = pbuffer.data(idx_kin_hh + 77);

    auto tk_xxxyy_yyyyy = pbuffer.data(idx_kin_hh + 78);

    auto tk_xxxyy_yyyyz = pbuffer.data(idx_kin_hh + 79);

    auto tk_xxxyy_yyyzz = pbuffer.data(idx_kin_hh + 80);

    auto tk_xxxyy_yyzzz = pbuffer.data(idx_kin_hh + 81);

    auto tk_xxxyy_yzzzz = pbuffer.data(idx_kin_hh + 82);

    auto tk_xxxyy_zzzzz = pbuffer.data(idx_kin_hh + 83);

    auto tk_xxxzz_xxxxx = pbuffer.data(idx_kin_hh + 105);

    auto tk_xxxzz_xxxxy = pbuffer.data(idx_kin_hh + 106);

    auto tk_xxxzz_xxxxz = pbuffer.data(idx_kin_hh + 107);

    auto tk_xxxzz_xxxyy = pbuffer.data(idx_kin_hh + 108);

    auto tk_xxxzz_xxxyz = pbuffer.data(idx_kin_hh + 109);

    auto tk_xxxzz_xxxzz = pbuffer.data(idx_kin_hh + 110);

    auto tk_xxxzz_xxyyy = pbuffer.data(idx_kin_hh + 111);

    auto tk_xxxzz_xxyyz = pbuffer.data(idx_kin_hh + 112);

    auto tk_xxxzz_xxyzz = pbuffer.data(idx_kin_hh + 113);

    auto tk_xxxzz_xxzzz = pbuffer.data(idx_kin_hh + 114);

    auto tk_xxxzz_xyyyy = pbuffer.data(idx_kin_hh + 115);

    auto tk_xxxzz_xyyyz = pbuffer.data(idx_kin_hh + 116);

    auto tk_xxxzz_xyyzz = pbuffer.data(idx_kin_hh + 117);

    auto tk_xxxzz_xyzzz = pbuffer.data(idx_kin_hh + 118);

    auto tk_xxxzz_xzzzz = pbuffer.data(idx_kin_hh + 119);

    auto tk_xxxzz_yyyyy = pbuffer.data(idx_kin_hh + 120);

    auto tk_xxxzz_yyyyz = pbuffer.data(idx_kin_hh + 121);

    auto tk_xxxzz_yyyzz = pbuffer.data(idx_kin_hh + 122);

    auto tk_xxxzz_yyzzz = pbuffer.data(idx_kin_hh + 123);

    auto tk_xxxzz_yzzzz = pbuffer.data(idx_kin_hh + 124);

    auto tk_xxxzz_zzzzz = pbuffer.data(idx_kin_hh + 125);

    auto tk_xxyyy_xxxxx = pbuffer.data(idx_kin_hh + 126);

    auto tk_xxyyy_xxxxy = pbuffer.data(idx_kin_hh + 127);

    auto tk_xxyyy_xxxxz = pbuffer.data(idx_kin_hh + 128);

    auto tk_xxyyy_xxxyy = pbuffer.data(idx_kin_hh + 129);

    auto tk_xxyyy_xxxyz = pbuffer.data(idx_kin_hh + 130);

    auto tk_xxyyy_xxxzz = pbuffer.data(idx_kin_hh + 131);

    auto tk_xxyyy_xxyyy = pbuffer.data(idx_kin_hh + 132);

    auto tk_xxyyy_xxyyz = pbuffer.data(idx_kin_hh + 133);

    auto tk_xxyyy_xxyzz = pbuffer.data(idx_kin_hh + 134);

    auto tk_xxyyy_xxzzz = pbuffer.data(idx_kin_hh + 135);

    auto tk_xxyyy_xyyyy = pbuffer.data(idx_kin_hh + 136);

    auto tk_xxyyy_xyyyz = pbuffer.data(idx_kin_hh + 137);

    auto tk_xxyyy_xyyzz = pbuffer.data(idx_kin_hh + 138);

    auto tk_xxyyy_xyzzz = pbuffer.data(idx_kin_hh + 139);

    auto tk_xxyyy_xzzzz = pbuffer.data(idx_kin_hh + 140);

    auto tk_xxyyy_yyyyy = pbuffer.data(idx_kin_hh + 141);

    auto tk_xxyyy_yyyyz = pbuffer.data(idx_kin_hh + 142);

    auto tk_xxyyy_yyyzz = pbuffer.data(idx_kin_hh + 143);

    auto tk_xxyyy_yyzzz = pbuffer.data(idx_kin_hh + 144);

    auto tk_xxyyy_yzzzz = pbuffer.data(idx_kin_hh + 145);

    auto tk_xxyyy_zzzzz = pbuffer.data(idx_kin_hh + 146);

    auto tk_xxzzz_xxxxx = pbuffer.data(idx_kin_hh + 189);

    auto tk_xxzzz_xxxxy = pbuffer.data(idx_kin_hh + 190);

    auto tk_xxzzz_xxxxz = pbuffer.data(idx_kin_hh + 191);

    auto tk_xxzzz_xxxyy = pbuffer.data(idx_kin_hh + 192);

    auto tk_xxzzz_xxxyz = pbuffer.data(idx_kin_hh + 193);

    auto tk_xxzzz_xxxzz = pbuffer.data(idx_kin_hh + 194);

    auto tk_xxzzz_xxyyy = pbuffer.data(idx_kin_hh + 195);

    auto tk_xxzzz_xxyyz = pbuffer.data(idx_kin_hh + 196);

    auto tk_xxzzz_xxyzz = pbuffer.data(idx_kin_hh + 197);

    auto tk_xxzzz_xxzzz = pbuffer.data(idx_kin_hh + 198);

    auto tk_xxzzz_xyyyy = pbuffer.data(idx_kin_hh + 199);

    auto tk_xxzzz_xyyyz = pbuffer.data(idx_kin_hh + 200);

    auto tk_xxzzz_xyyzz = pbuffer.data(idx_kin_hh + 201);

    auto tk_xxzzz_xyzzz = pbuffer.data(idx_kin_hh + 202);

    auto tk_xxzzz_xzzzz = pbuffer.data(idx_kin_hh + 203);

    auto tk_xxzzz_yyyyy = pbuffer.data(idx_kin_hh + 204);

    auto tk_xxzzz_yyyyz = pbuffer.data(idx_kin_hh + 205);

    auto tk_xxzzz_yyyzz = pbuffer.data(idx_kin_hh + 206);

    auto tk_xxzzz_yyzzz = pbuffer.data(idx_kin_hh + 207);

    auto tk_xxzzz_yzzzz = pbuffer.data(idx_kin_hh + 208);

    auto tk_xxzzz_zzzzz = pbuffer.data(idx_kin_hh + 209);

    auto tk_xyyyy_xxxxy = pbuffer.data(idx_kin_hh + 211);

    auto tk_xyyyy_xxxyy = pbuffer.data(idx_kin_hh + 213);

    auto tk_xyyyy_xxxyz = pbuffer.data(idx_kin_hh + 214);

    auto tk_xyyyy_xxyyy = pbuffer.data(idx_kin_hh + 216);

    auto tk_xyyyy_xxyyz = pbuffer.data(idx_kin_hh + 217);

    auto tk_xyyyy_xxyzz = pbuffer.data(idx_kin_hh + 218);

    auto tk_xyyyy_xyyyy = pbuffer.data(idx_kin_hh + 220);

    auto tk_xyyyy_xyyyz = pbuffer.data(idx_kin_hh + 221);

    auto tk_xyyyy_xyyzz = pbuffer.data(idx_kin_hh + 222);

    auto tk_xyyyy_xyzzz = pbuffer.data(idx_kin_hh + 223);

    auto tk_xyyyy_yyyyy = pbuffer.data(idx_kin_hh + 225);

    auto tk_xyyyy_yyyyz = pbuffer.data(idx_kin_hh + 226);

    auto tk_xyyyy_yyyzz = pbuffer.data(idx_kin_hh + 227);

    auto tk_xyyyy_yyzzz = pbuffer.data(idx_kin_hh + 228);

    auto tk_xyyyy_yzzzz = pbuffer.data(idx_kin_hh + 229);

    auto tk_xyyzz_xxxyz = pbuffer.data(idx_kin_hh + 256);

    auto tk_xyyzz_xxyyz = pbuffer.data(idx_kin_hh + 259);

    auto tk_xyyzz_xxyzz = pbuffer.data(idx_kin_hh + 260);

    auto tk_xyyzz_xyyyz = pbuffer.data(idx_kin_hh + 263);

    auto tk_xyyzz_xyyzz = pbuffer.data(idx_kin_hh + 264);

    auto tk_xyyzz_xyzzz = pbuffer.data(idx_kin_hh + 265);

    auto tk_xyyzz_yyyyz = pbuffer.data(idx_kin_hh + 268);

    auto tk_xyyzz_yyyzz = pbuffer.data(idx_kin_hh + 269);

    auto tk_xyyzz_yyzzz = pbuffer.data(idx_kin_hh + 270);

    auto tk_xyyzz_yzzzz = pbuffer.data(idx_kin_hh + 271);

    auto tk_xzzzz_xxxxz = pbuffer.data(idx_kin_hh + 296);

    auto tk_xzzzz_xxxyz = pbuffer.data(idx_kin_hh + 298);

    auto tk_xzzzz_xxxzz = pbuffer.data(idx_kin_hh + 299);

    auto tk_xzzzz_xxyyz = pbuffer.data(idx_kin_hh + 301);

    auto tk_xzzzz_xxyzz = pbuffer.data(idx_kin_hh + 302);

    auto tk_xzzzz_xxzzz = pbuffer.data(idx_kin_hh + 303);

    auto tk_xzzzz_xyyyz = pbuffer.data(idx_kin_hh + 305);

    auto tk_xzzzz_xyyzz = pbuffer.data(idx_kin_hh + 306);

    auto tk_xzzzz_xyzzz = pbuffer.data(idx_kin_hh + 307);

    auto tk_xzzzz_xzzzz = pbuffer.data(idx_kin_hh + 308);

    auto tk_xzzzz_yyyyz = pbuffer.data(idx_kin_hh + 310);

    auto tk_xzzzz_yyyzz = pbuffer.data(idx_kin_hh + 311);

    auto tk_xzzzz_yyzzz = pbuffer.data(idx_kin_hh + 312);

    auto tk_xzzzz_yzzzz = pbuffer.data(idx_kin_hh + 313);

    auto tk_xzzzz_zzzzz = pbuffer.data(idx_kin_hh + 314);

    auto tk_yyyyy_xxxxx = pbuffer.data(idx_kin_hh + 315);

    auto tk_yyyyy_xxxxy = pbuffer.data(idx_kin_hh + 316);

    auto tk_yyyyy_xxxxz = pbuffer.data(idx_kin_hh + 317);

    auto tk_yyyyy_xxxyy = pbuffer.data(idx_kin_hh + 318);

    auto tk_yyyyy_xxxyz = pbuffer.data(idx_kin_hh + 319);

    auto tk_yyyyy_xxxzz = pbuffer.data(idx_kin_hh + 320);

    auto tk_yyyyy_xxyyy = pbuffer.data(idx_kin_hh + 321);

    auto tk_yyyyy_xxyyz = pbuffer.data(idx_kin_hh + 322);

    auto tk_yyyyy_xxyzz = pbuffer.data(idx_kin_hh + 323);

    auto tk_yyyyy_xxzzz = pbuffer.data(idx_kin_hh + 324);

    auto tk_yyyyy_xyyyy = pbuffer.data(idx_kin_hh + 325);

    auto tk_yyyyy_xyyyz = pbuffer.data(idx_kin_hh + 326);

    auto tk_yyyyy_xyyzz = pbuffer.data(idx_kin_hh + 327);

    auto tk_yyyyy_xyzzz = pbuffer.data(idx_kin_hh + 328);

    auto tk_yyyyy_xzzzz = pbuffer.data(idx_kin_hh + 329);

    auto tk_yyyyy_yyyyy = pbuffer.data(idx_kin_hh + 330);

    auto tk_yyyyy_yyyyz = pbuffer.data(idx_kin_hh + 331);

    auto tk_yyyyy_yyyzz = pbuffer.data(idx_kin_hh + 332);

    auto tk_yyyyy_yyzzz = pbuffer.data(idx_kin_hh + 333);

    auto tk_yyyyy_yzzzz = pbuffer.data(idx_kin_hh + 334);

    auto tk_yyyyy_zzzzz = pbuffer.data(idx_kin_hh + 335);

    auto tk_yyyyz_xxxxz = pbuffer.data(idx_kin_hh + 338);

    auto tk_yyyyz_xxxyz = pbuffer.data(idx_kin_hh + 340);

    auto tk_yyyyz_xxxzz = pbuffer.data(idx_kin_hh + 341);

    auto tk_yyyyz_xxyyz = pbuffer.data(idx_kin_hh + 343);

    auto tk_yyyyz_xxyzz = pbuffer.data(idx_kin_hh + 344);

    auto tk_yyyyz_xxzzz = pbuffer.data(idx_kin_hh + 345);

    auto tk_yyyyz_xyyyz = pbuffer.data(idx_kin_hh + 347);

    auto tk_yyyyz_xyyzz = pbuffer.data(idx_kin_hh + 348);

    auto tk_yyyyz_xyzzz = pbuffer.data(idx_kin_hh + 349);

    auto tk_yyyyz_xzzzz = pbuffer.data(idx_kin_hh + 350);

    auto tk_yyyyz_yyyyz = pbuffer.data(idx_kin_hh + 352);

    auto tk_yyyyz_yyyzz = pbuffer.data(idx_kin_hh + 353);

    auto tk_yyyyz_yyzzz = pbuffer.data(idx_kin_hh + 354);

    auto tk_yyyyz_yzzzz = pbuffer.data(idx_kin_hh + 355);

    auto tk_yyyyz_zzzzz = pbuffer.data(idx_kin_hh + 356);

    auto tk_yyyzz_xxxxx = pbuffer.data(idx_kin_hh + 357);

    auto tk_yyyzz_xxxxy = pbuffer.data(idx_kin_hh + 358);

    auto tk_yyyzz_xxxxz = pbuffer.data(idx_kin_hh + 359);

    auto tk_yyyzz_xxxyy = pbuffer.data(idx_kin_hh + 360);

    auto tk_yyyzz_xxxyz = pbuffer.data(idx_kin_hh + 361);

    auto tk_yyyzz_xxxzz = pbuffer.data(idx_kin_hh + 362);

    auto tk_yyyzz_xxyyy = pbuffer.data(idx_kin_hh + 363);

    auto tk_yyyzz_xxyyz = pbuffer.data(idx_kin_hh + 364);

    auto tk_yyyzz_xxyzz = pbuffer.data(idx_kin_hh + 365);

    auto tk_yyyzz_xxzzz = pbuffer.data(idx_kin_hh + 366);

    auto tk_yyyzz_xyyyy = pbuffer.data(idx_kin_hh + 367);

    auto tk_yyyzz_xyyyz = pbuffer.data(idx_kin_hh + 368);

    auto tk_yyyzz_xyyzz = pbuffer.data(idx_kin_hh + 369);

    auto tk_yyyzz_xyzzz = pbuffer.data(idx_kin_hh + 370);

    auto tk_yyyzz_xzzzz = pbuffer.data(idx_kin_hh + 371);

    auto tk_yyyzz_yyyyy = pbuffer.data(idx_kin_hh + 372);

    auto tk_yyyzz_yyyyz = pbuffer.data(idx_kin_hh + 373);

    auto tk_yyyzz_yyyzz = pbuffer.data(idx_kin_hh + 374);

    auto tk_yyyzz_yyzzz = pbuffer.data(idx_kin_hh + 375);

    auto tk_yyyzz_yzzzz = pbuffer.data(idx_kin_hh + 376);

    auto tk_yyyzz_zzzzz = pbuffer.data(idx_kin_hh + 377);

    auto tk_yyzzz_xxxxx = pbuffer.data(idx_kin_hh + 378);

    auto tk_yyzzz_xxxxy = pbuffer.data(idx_kin_hh + 379);

    auto tk_yyzzz_xxxxz = pbuffer.data(idx_kin_hh + 380);

    auto tk_yyzzz_xxxyy = pbuffer.data(idx_kin_hh + 381);

    auto tk_yyzzz_xxxyz = pbuffer.data(idx_kin_hh + 382);

    auto tk_yyzzz_xxxzz = pbuffer.data(idx_kin_hh + 383);

    auto tk_yyzzz_xxyyy = pbuffer.data(idx_kin_hh + 384);

    auto tk_yyzzz_xxyyz = pbuffer.data(idx_kin_hh + 385);

    auto tk_yyzzz_xxyzz = pbuffer.data(idx_kin_hh + 386);

    auto tk_yyzzz_xxzzz = pbuffer.data(idx_kin_hh + 387);

    auto tk_yyzzz_xyyyy = pbuffer.data(idx_kin_hh + 388);

    auto tk_yyzzz_xyyyz = pbuffer.data(idx_kin_hh + 389);

    auto tk_yyzzz_xyyzz = pbuffer.data(idx_kin_hh + 390);

    auto tk_yyzzz_xyzzz = pbuffer.data(idx_kin_hh + 391);

    auto tk_yyzzz_xzzzz = pbuffer.data(idx_kin_hh + 392);

    auto tk_yyzzz_yyyyy = pbuffer.data(idx_kin_hh + 393);

    auto tk_yyzzz_yyyyz = pbuffer.data(idx_kin_hh + 394);

    auto tk_yyzzz_yyyzz = pbuffer.data(idx_kin_hh + 395);

    auto tk_yyzzz_yyzzz = pbuffer.data(idx_kin_hh + 396);

    auto tk_yyzzz_yzzzz = pbuffer.data(idx_kin_hh + 397);

    auto tk_yyzzz_zzzzz = pbuffer.data(idx_kin_hh + 398);

    auto tk_yzzzz_xxxxy = pbuffer.data(idx_kin_hh + 400);

    auto tk_yzzzz_xxxxz = pbuffer.data(idx_kin_hh + 401);

    auto tk_yzzzz_xxxyy = pbuffer.data(idx_kin_hh + 402);

    auto tk_yzzzz_xxxyz = pbuffer.data(idx_kin_hh + 403);

    auto tk_yzzzz_xxxzz = pbuffer.data(idx_kin_hh + 404);

    auto tk_yzzzz_xxyyy = pbuffer.data(idx_kin_hh + 405);

    auto tk_yzzzz_xxyyz = pbuffer.data(idx_kin_hh + 406);

    auto tk_yzzzz_xxyzz = pbuffer.data(idx_kin_hh + 407);

    auto tk_yzzzz_xxzzz = pbuffer.data(idx_kin_hh + 408);

    auto tk_yzzzz_xyyyy = pbuffer.data(idx_kin_hh + 409);

    auto tk_yzzzz_xyyyz = pbuffer.data(idx_kin_hh + 410);

    auto tk_yzzzz_xyyzz = pbuffer.data(idx_kin_hh + 411);

    auto tk_yzzzz_xyzzz = pbuffer.data(idx_kin_hh + 412);

    auto tk_yzzzz_xzzzz = pbuffer.data(idx_kin_hh + 413);

    auto tk_yzzzz_yyyyy = pbuffer.data(idx_kin_hh + 414);

    auto tk_yzzzz_yyyyz = pbuffer.data(idx_kin_hh + 415);

    auto tk_yzzzz_yyyzz = pbuffer.data(idx_kin_hh + 416);

    auto tk_yzzzz_yyzzz = pbuffer.data(idx_kin_hh + 417);

    auto tk_yzzzz_yzzzz = pbuffer.data(idx_kin_hh + 418);

    auto tk_yzzzz_zzzzz = pbuffer.data(idx_kin_hh + 419);

    auto tk_zzzzz_xxxxx = pbuffer.data(idx_kin_hh + 420);

    auto tk_zzzzz_xxxxy = pbuffer.data(idx_kin_hh + 421);

    auto tk_zzzzz_xxxxz = pbuffer.data(idx_kin_hh + 422);

    auto tk_zzzzz_xxxyy = pbuffer.data(idx_kin_hh + 423);

    auto tk_zzzzz_xxxyz = pbuffer.data(idx_kin_hh + 424);

    auto tk_zzzzz_xxxzz = pbuffer.data(idx_kin_hh + 425);

    auto tk_zzzzz_xxyyy = pbuffer.data(idx_kin_hh + 426);

    auto tk_zzzzz_xxyyz = pbuffer.data(idx_kin_hh + 427);

    auto tk_zzzzz_xxyzz = pbuffer.data(idx_kin_hh + 428);

    auto tk_zzzzz_xxzzz = pbuffer.data(idx_kin_hh + 429);

    auto tk_zzzzz_xyyyy = pbuffer.data(idx_kin_hh + 430);

    auto tk_zzzzz_xyyyz = pbuffer.data(idx_kin_hh + 431);

    auto tk_zzzzz_xyyzz = pbuffer.data(idx_kin_hh + 432);

    auto tk_zzzzz_xyzzz = pbuffer.data(idx_kin_hh + 433);

    auto tk_zzzzz_xzzzz = pbuffer.data(idx_kin_hh + 434);

    auto tk_zzzzz_yyyyy = pbuffer.data(idx_kin_hh + 435);

    auto tk_zzzzz_yyyyz = pbuffer.data(idx_kin_hh + 436);

    auto tk_zzzzz_yyyzz = pbuffer.data(idx_kin_hh + 437);

    auto tk_zzzzz_yyzzz = pbuffer.data(idx_kin_hh + 438);

    auto tk_zzzzz_yzzzz = pbuffer.data(idx_kin_hh + 439);

    auto tk_zzzzz_zzzzz = pbuffer.data(idx_kin_hh + 440);

    // Set up components of auxiliary buffer : HI

    auto tk_xxxxx_xxxxxx = pbuffer.data(idx_kin_hi);

    auto tk_xxxxx_xxxxxy = pbuffer.data(idx_kin_hi + 1);

    auto tk_xxxxx_xxxxxz = pbuffer.data(idx_kin_hi + 2);

    auto tk_xxxxx_xxxxyy = pbuffer.data(idx_kin_hi + 3);

    auto tk_xxxxx_xxxxyz = pbuffer.data(idx_kin_hi + 4);

    auto tk_xxxxx_xxxxzz = pbuffer.data(idx_kin_hi + 5);

    auto tk_xxxxx_xxxyyy = pbuffer.data(idx_kin_hi + 6);

    auto tk_xxxxx_xxxyyz = pbuffer.data(idx_kin_hi + 7);

    auto tk_xxxxx_xxxyzz = pbuffer.data(idx_kin_hi + 8);

    auto tk_xxxxx_xxxzzz = pbuffer.data(idx_kin_hi + 9);

    auto tk_xxxxx_xxyyyy = pbuffer.data(idx_kin_hi + 10);

    auto tk_xxxxx_xxyyyz = pbuffer.data(idx_kin_hi + 11);

    auto tk_xxxxx_xxyyzz = pbuffer.data(idx_kin_hi + 12);

    auto tk_xxxxx_xxyzzz = pbuffer.data(idx_kin_hi + 13);

    auto tk_xxxxx_xxzzzz = pbuffer.data(idx_kin_hi + 14);

    auto tk_xxxxx_xyyyyy = pbuffer.data(idx_kin_hi + 15);

    auto tk_xxxxx_xyyyyz = pbuffer.data(idx_kin_hi + 16);

    auto tk_xxxxx_xyyyzz = pbuffer.data(idx_kin_hi + 17);

    auto tk_xxxxx_xyyzzz = pbuffer.data(idx_kin_hi + 18);

    auto tk_xxxxx_xyzzzz = pbuffer.data(idx_kin_hi + 19);

    auto tk_xxxxx_xzzzzz = pbuffer.data(idx_kin_hi + 20);

    auto tk_xxxxx_yyyyyy = pbuffer.data(idx_kin_hi + 21);

    auto tk_xxxxx_yyyyyz = pbuffer.data(idx_kin_hi + 22);

    auto tk_xxxxx_yyyyzz = pbuffer.data(idx_kin_hi + 23);

    auto tk_xxxxx_yyyzzz = pbuffer.data(idx_kin_hi + 24);

    auto tk_xxxxx_yyzzzz = pbuffer.data(idx_kin_hi + 25);

    auto tk_xxxxx_yzzzzz = pbuffer.data(idx_kin_hi + 26);

    auto tk_xxxxx_zzzzzz = pbuffer.data(idx_kin_hi + 27);

    auto tk_xxxxy_xxxxxx = pbuffer.data(idx_kin_hi + 28);

    auto tk_xxxxy_xxxxxy = pbuffer.data(idx_kin_hi + 29);

    auto tk_xxxxy_xxxxxz = pbuffer.data(idx_kin_hi + 30);

    auto tk_xxxxy_xxxxyy = pbuffer.data(idx_kin_hi + 31);

    auto tk_xxxxy_xxxxzz = pbuffer.data(idx_kin_hi + 33);

    auto tk_xxxxy_xxxyyy = pbuffer.data(idx_kin_hi + 34);

    auto tk_xxxxy_xxxzzz = pbuffer.data(idx_kin_hi + 37);

    auto tk_xxxxy_xxyyyy = pbuffer.data(idx_kin_hi + 38);

    auto tk_xxxxy_xxzzzz = pbuffer.data(idx_kin_hi + 42);

    auto tk_xxxxy_xyyyyy = pbuffer.data(idx_kin_hi + 43);

    auto tk_xxxxy_xzzzzz = pbuffer.data(idx_kin_hi + 48);

    auto tk_xxxxy_yyyyyy = pbuffer.data(idx_kin_hi + 49);

    auto tk_xxxxz_xxxxxx = pbuffer.data(idx_kin_hi + 56);

    auto tk_xxxxz_xxxxxy = pbuffer.data(idx_kin_hi + 57);

    auto tk_xxxxz_xxxxxz = pbuffer.data(idx_kin_hi + 58);

    auto tk_xxxxz_xxxxyy = pbuffer.data(idx_kin_hi + 59);

    auto tk_xxxxz_xxxxyz = pbuffer.data(idx_kin_hi + 60);

    auto tk_xxxxz_xxxxzz = pbuffer.data(idx_kin_hi + 61);

    auto tk_xxxxz_xxxyyy = pbuffer.data(idx_kin_hi + 62);

    auto tk_xxxxz_xxxyyz = pbuffer.data(idx_kin_hi + 63);

    auto tk_xxxxz_xxxyzz = pbuffer.data(idx_kin_hi + 64);

    auto tk_xxxxz_xxxzzz = pbuffer.data(idx_kin_hi + 65);

    auto tk_xxxxz_xxyyyy = pbuffer.data(idx_kin_hi + 66);

    auto tk_xxxxz_xxyyyz = pbuffer.data(idx_kin_hi + 67);

    auto tk_xxxxz_xxyyzz = pbuffer.data(idx_kin_hi + 68);

    auto tk_xxxxz_xxyzzz = pbuffer.data(idx_kin_hi + 69);

    auto tk_xxxxz_xxzzzz = pbuffer.data(idx_kin_hi + 70);

    auto tk_xxxxz_xyyyyy = pbuffer.data(idx_kin_hi + 71);

    auto tk_xxxxz_xyyyyz = pbuffer.data(idx_kin_hi + 72);

    auto tk_xxxxz_xyyyzz = pbuffer.data(idx_kin_hi + 73);

    auto tk_xxxxz_xyyzzz = pbuffer.data(idx_kin_hi + 74);

    auto tk_xxxxz_xyzzzz = pbuffer.data(idx_kin_hi + 75);

    auto tk_xxxxz_xzzzzz = pbuffer.data(idx_kin_hi + 76);

    auto tk_xxxxz_yyyyyz = pbuffer.data(idx_kin_hi + 78);

    auto tk_xxxxz_yyyyzz = pbuffer.data(idx_kin_hi + 79);

    auto tk_xxxxz_yyyzzz = pbuffer.data(idx_kin_hi + 80);

    auto tk_xxxxz_yyzzzz = pbuffer.data(idx_kin_hi + 81);

    auto tk_xxxxz_yzzzzz = pbuffer.data(idx_kin_hi + 82);

    auto tk_xxxxz_zzzzzz = pbuffer.data(idx_kin_hi + 83);

    auto tk_xxxyy_xxxxxx = pbuffer.data(idx_kin_hi + 84);

    auto tk_xxxyy_xxxxxy = pbuffer.data(idx_kin_hi + 85);

    auto tk_xxxyy_xxxxxz = pbuffer.data(idx_kin_hi + 86);

    auto tk_xxxyy_xxxxyy = pbuffer.data(idx_kin_hi + 87);

    auto tk_xxxyy_xxxxyz = pbuffer.data(idx_kin_hi + 88);

    auto tk_xxxyy_xxxxzz = pbuffer.data(idx_kin_hi + 89);

    auto tk_xxxyy_xxxyyy = pbuffer.data(idx_kin_hi + 90);

    auto tk_xxxyy_xxxyyz = pbuffer.data(idx_kin_hi + 91);

    auto tk_xxxyy_xxxyzz = pbuffer.data(idx_kin_hi + 92);

    auto tk_xxxyy_xxxzzz = pbuffer.data(idx_kin_hi + 93);

    auto tk_xxxyy_xxyyyy = pbuffer.data(idx_kin_hi + 94);

    auto tk_xxxyy_xxyyyz = pbuffer.data(idx_kin_hi + 95);

    auto tk_xxxyy_xxyyzz = pbuffer.data(idx_kin_hi + 96);

    auto tk_xxxyy_xxyzzz = pbuffer.data(idx_kin_hi + 97);

    auto tk_xxxyy_xxzzzz = pbuffer.data(idx_kin_hi + 98);

    auto tk_xxxyy_xyyyyy = pbuffer.data(idx_kin_hi + 99);

    auto tk_xxxyy_xyyyyz = pbuffer.data(idx_kin_hi + 100);

    auto tk_xxxyy_xyyyzz = pbuffer.data(idx_kin_hi + 101);

    auto tk_xxxyy_xyyzzz = pbuffer.data(idx_kin_hi + 102);

    auto tk_xxxyy_xyzzzz = pbuffer.data(idx_kin_hi + 103);

    auto tk_xxxyy_xzzzzz = pbuffer.data(idx_kin_hi + 104);

    auto tk_xxxyy_yyyyyy = pbuffer.data(idx_kin_hi + 105);

    auto tk_xxxyy_yyyyyz = pbuffer.data(idx_kin_hi + 106);

    auto tk_xxxyy_yyyyzz = pbuffer.data(idx_kin_hi + 107);

    auto tk_xxxyy_yyyzzz = pbuffer.data(idx_kin_hi + 108);

    auto tk_xxxyy_yyzzzz = pbuffer.data(idx_kin_hi + 109);

    auto tk_xxxyy_yzzzzz = pbuffer.data(idx_kin_hi + 110);

    auto tk_xxxyy_zzzzzz = pbuffer.data(idx_kin_hi + 111);

    auto tk_xxxzz_xxxxxx = pbuffer.data(idx_kin_hi + 140);

    auto tk_xxxzz_xxxxxy = pbuffer.data(idx_kin_hi + 141);

    auto tk_xxxzz_xxxxxz = pbuffer.data(idx_kin_hi + 142);

    auto tk_xxxzz_xxxxyy = pbuffer.data(idx_kin_hi + 143);

    auto tk_xxxzz_xxxxyz = pbuffer.data(idx_kin_hi + 144);

    auto tk_xxxzz_xxxxzz = pbuffer.data(idx_kin_hi + 145);

    auto tk_xxxzz_xxxyyy = pbuffer.data(idx_kin_hi + 146);

    auto tk_xxxzz_xxxyyz = pbuffer.data(idx_kin_hi + 147);

    auto tk_xxxzz_xxxyzz = pbuffer.data(idx_kin_hi + 148);

    auto tk_xxxzz_xxxzzz = pbuffer.data(idx_kin_hi + 149);

    auto tk_xxxzz_xxyyyy = pbuffer.data(idx_kin_hi + 150);

    auto tk_xxxzz_xxyyyz = pbuffer.data(idx_kin_hi + 151);

    auto tk_xxxzz_xxyyzz = pbuffer.data(idx_kin_hi + 152);

    auto tk_xxxzz_xxyzzz = pbuffer.data(idx_kin_hi + 153);

    auto tk_xxxzz_xxzzzz = pbuffer.data(idx_kin_hi + 154);

    auto tk_xxxzz_xyyyyy = pbuffer.data(idx_kin_hi + 155);

    auto tk_xxxzz_xyyyyz = pbuffer.data(idx_kin_hi + 156);

    auto tk_xxxzz_xyyyzz = pbuffer.data(idx_kin_hi + 157);

    auto tk_xxxzz_xyyzzz = pbuffer.data(idx_kin_hi + 158);

    auto tk_xxxzz_xyzzzz = pbuffer.data(idx_kin_hi + 159);

    auto tk_xxxzz_xzzzzz = pbuffer.data(idx_kin_hi + 160);

    auto tk_xxxzz_yyyyyy = pbuffer.data(idx_kin_hi + 161);

    auto tk_xxxzz_yyyyyz = pbuffer.data(idx_kin_hi + 162);

    auto tk_xxxzz_yyyyzz = pbuffer.data(idx_kin_hi + 163);

    auto tk_xxxzz_yyyzzz = pbuffer.data(idx_kin_hi + 164);

    auto tk_xxxzz_yyzzzz = pbuffer.data(idx_kin_hi + 165);

    auto tk_xxxzz_yzzzzz = pbuffer.data(idx_kin_hi + 166);

    auto tk_xxxzz_zzzzzz = pbuffer.data(idx_kin_hi + 167);

    auto tk_xxyyy_xxxxxx = pbuffer.data(idx_kin_hi + 168);

    auto tk_xxyyy_xxxxxy = pbuffer.data(idx_kin_hi + 169);

    auto tk_xxyyy_xxxxxz = pbuffer.data(idx_kin_hi + 170);

    auto tk_xxyyy_xxxxyy = pbuffer.data(idx_kin_hi + 171);

    auto tk_xxyyy_xxxxyz = pbuffer.data(idx_kin_hi + 172);

    auto tk_xxyyy_xxxxzz = pbuffer.data(idx_kin_hi + 173);

    auto tk_xxyyy_xxxyyy = pbuffer.data(idx_kin_hi + 174);

    auto tk_xxyyy_xxxyyz = pbuffer.data(idx_kin_hi + 175);

    auto tk_xxyyy_xxxyzz = pbuffer.data(idx_kin_hi + 176);

    auto tk_xxyyy_xxxzzz = pbuffer.data(idx_kin_hi + 177);

    auto tk_xxyyy_xxyyyy = pbuffer.data(idx_kin_hi + 178);

    auto tk_xxyyy_xxyyyz = pbuffer.data(idx_kin_hi + 179);

    auto tk_xxyyy_xxyyzz = pbuffer.data(idx_kin_hi + 180);

    auto tk_xxyyy_xxyzzz = pbuffer.data(idx_kin_hi + 181);

    auto tk_xxyyy_xxzzzz = pbuffer.data(idx_kin_hi + 182);

    auto tk_xxyyy_xyyyyy = pbuffer.data(idx_kin_hi + 183);

    auto tk_xxyyy_xyyyyz = pbuffer.data(idx_kin_hi + 184);

    auto tk_xxyyy_xyyyzz = pbuffer.data(idx_kin_hi + 185);

    auto tk_xxyyy_xyyzzz = pbuffer.data(idx_kin_hi + 186);

    auto tk_xxyyy_xyzzzz = pbuffer.data(idx_kin_hi + 187);

    auto tk_xxyyy_xzzzzz = pbuffer.data(idx_kin_hi + 188);

    auto tk_xxyyy_yyyyyy = pbuffer.data(idx_kin_hi + 189);

    auto tk_xxyyy_yyyyyz = pbuffer.data(idx_kin_hi + 190);

    auto tk_xxyyy_yyyyzz = pbuffer.data(idx_kin_hi + 191);

    auto tk_xxyyy_yyyzzz = pbuffer.data(idx_kin_hi + 192);

    auto tk_xxyyy_yyzzzz = pbuffer.data(idx_kin_hi + 193);

    auto tk_xxyyy_yzzzzz = pbuffer.data(idx_kin_hi + 194);

    auto tk_xxyyy_zzzzzz = pbuffer.data(idx_kin_hi + 195);

    auto tk_xxyyz_xxxxxy = pbuffer.data(idx_kin_hi + 197);

    auto tk_xxyyz_xxxxyy = pbuffer.data(idx_kin_hi + 199);

    auto tk_xxyyz_xxxyyy = pbuffer.data(idx_kin_hi + 202);

    auto tk_xxyyz_xxyyyy = pbuffer.data(idx_kin_hi + 206);

    auto tk_xxyyz_xyyyyy = pbuffer.data(idx_kin_hi + 211);

    auto tk_xxyzz_xxxxxx = pbuffer.data(idx_kin_hi + 224);

    auto tk_xxyzz_xxxxxz = pbuffer.data(idx_kin_hi + 226);

    auto tk_xxyzz_xxxxzz = pbuffer.data(idx_kin_hi + 229);

    auto tk_xxyzz_xxxzzz = pbuffer.data(idx_kin_hi + 233);

    auto tk_xxyzz_xxzzzz = pbuffer.data(idx_kin_hi + 238);

    auto tk_xxyzz_xzzzzz = pbuffer.data(idx_kin_hi + 244);

    auto tk_xxzzz_xxxxxx = pbuffer.data(idx_kin_hi + 252);

    auto tk_xxzzz_xxxxxy = pbuffer.data(idx_kin_hi + 253);

    auto tk_xxzzz_xxxxxz = pbuffer.data(idx_kin_hi + 254);

    auto tk_xxzzz_xxxxyy = pbuffer.data(idx_kin_hi + 255);

    auto tk_xxzzz_xxxxyz = pbuffer.data(idx_kin_hi + 256);

    auto tk_xxzzz_xxxxzz = pbuffer.data(idx_kin_hi + 257);

    auto tk_xxzzz_xxxyyy = pbuffer.data(idx_kin_hi + 258);

    auto tk_xxzzz_xxxyyz = pbuffer.data(idx_kin_hi + 259);

    auto tk_xxzzz_xxxyzz = pbuffer.data(idx_kin_hi + 260);

    auto tk_xxzzz_xxxzzz = pbuffer.data(idx_kin_hi + 261);

    auto tk_xxzzz_xxyyyy = pbuffer.data(idx_kin_hi + 262);

    auto tk_xxzzz_xxyyyz = pbuffer.data(idx_kin_hi + 263);

    auto tk_xxzzz_xxyyzz = pbuffer.data(idx_kin_hi + 264);

    auto tk_xxzzz_xxyzzz = pbuffer.data(idx_kin_hi + 265);

    auto tk_xxzzz_xxzzzz = pbuffer.data(idx_kin_hi + 266);

    auto tk_xxzzz_xyyyyy = pbuffer.data(idx_kin_hi + 267);

    auto tk_xxzzz_xyyyyz = pbuffer.data(idx_kin_hi + 268);

    auto tk_xxzzz_xyyyzz = pbuffer.data(idx_kin_hi + 269);

    auto tk_xxzzz_xyyzzz = pbuffer.data(idx_kin_hi + 270);

    auto tk_xxzzz_xyzzzz = pbuffer.data(idx_kin_hi + 271);

    auto tk_xxzzz_xzzzzz = pbuffer.data(idx_kin_hi + 272);

    auto tk_xxzzz_yyyyyy = pbuffer.data(idx_kin_hi + 273);

    auto tk_xxzzz_yyyyyz = pbuffer.data(idx_kin_hi + 274);

    auto tk_xxzzz_yyyyzz = pbuffer.data(idx_kin_hi + 275);

    auto tk_xxzzz_yyyzzz = pbuffer.data(idx_kin_hi + 276);

    auto tk_xxzzz_yyzzzz = pbuffer.data(idx_kin_hi + 277);

    auto tk_xxzzz_yzzzzz = pbuffer.data(idx_kin_hi + 278);

    auto tk_xxzzz_zzzzzz = pbuffer.data(idx_kin_hi + 279);

    auto tk_xyyyy_xxxxxx = pbuffer.data(idx_kin_hi + 280);

    auto tk_xyyyy_xxxxxy = pbuffer.data(idx_kin_hi + 281);

    auto tk_xyyyy_xxxxyy = pbuffer.data(idx_kin_hi + 283);

    auto tk_xyyyy_xxxxyz = pbuffer.data(idx_kin_hi + 284);

    auto tk_xyyyy_xxxyyy = pbuffer.data(idx_kin_hi + 286);

    auto tk_xyyyy_xxxyyz = pbuffer.data(idx_kin_hi + 287);

    auto tk_xyyyy_xxxyzz = pbuffer.data(idx_kin_hi + 288);

    auto tk_xyyyy_xxyyyy = pbuffer.data(idx_kin_hi + 290);

    auto tk_xyyyy_xxyyyz = pbuffer.data(idx_kin_hi + 291);

    auto tk_xyyyy_xxyyzz = pbuffer.data(idx_kin_hi + 292);

    auto tk_xyyyy_xxyzzz = pbuffer.data(idx_kin_hi + 293);

    auto tk_xyyyy_xyyyyy = pbuffer.data(idx_kin_hi + 295);

    auto tk_xyyyy_xyyyyz = pbuffer.data(idx_kin_hi + 296);

    auto tk_xyyyy_xyyyzz = pbuffer.data(idx_kin_hi + 297);

    auto tk_xyyyy_xyyzzz = pbuffer.data(idx_kin_hi + 298);

    auto tk_xyyyy_xyzzzz = pbuffer.data(idx_kin_hi + 299);

    auto tk_xyyyy_yyyyyy = pbuffer.data(idx_kin_hi + 301);

    auto tk_xyyyy_yyyyyz = pbuffer.data(idx_kin_hi + 302);

    auto tk_xyyyy_yyyyzz = pbuffer.data(idx_kin_hi + 303);

    auto tk_xyyyy_yyyzzz = pbuffer.data(idx_kin_hi + 304);

    auto tk_xyyyy_yyzzzz = pbuffer.data(idx_kin_hi + 305);

    auto tk_xyyyy_yzzzzz = pbuffer.data(idx_kin_hi + 306);

    auto tk_xyyyy_zzzzzz = pbuffer.data(idx_kin_hi + 307);

    auto tk_xyyzz_xxxxyz = pbuffer.data(idx_kin_hi + 340);

    auto tk_xyyzz_xxxyyz = pbuffer.data(idx_kin_hi + 343);

    auto tk_xyyzz_xxxyzz = pbuffer.data(idx_kin_hi + 344);

    auto tk_xyyzz_xxyyyz = pbuffer.data(idx_kin_hi + 347);

    auto tk_xyyzz_xxyyzz = pbuffer.data(idx_kin_hi + 348);

    auto tk_xyyzz_xxyzzz = pbuffer.data(idx_kin_hi + 349);

    auto tk_xyyzz_xyyyyz = pbuffer.data(idx_kin_hi + 352);

    auto tk_xyyzz_xyyyzz = pbuffer.data(idx_kin_hi + 353);

    auto tk_xyyzz_xyyzzz = pbuffer.data(idx_kin_hi + 354);

    auto tk_xyyzz_xyzzzz = pbuffer.data(idx_kin_hi + 355);

    auto tk_xyyzz_yyyyyy = pbuffer.data(idx_kin_hi + 357);

    auto tk_xyyzz_yyyyyz = pbuffer.data(idx_kin_hi + 358);

    auto tk_xyyzz_yyyyzz = pbuffer.data(idx_kin_hi + 359);

    auto tk_xyyzz_yyyzzz = pbuffer.data(idx_kin_hi + 360);

    auto tk_xyyzz_yyzzzz = pbuffer.data(idx_kin_hi + 361);

    auto tk_xyyzz_yzzzzz = pbuffer.data(idx_kin_hi + 362);

    auto tk_xyyzz_zzzzzz = pbuffer.data(idx_kin_hi + 363);

    auto tk_xzzzz_xxxxxx = pbuffer.data(idx_kin_hi + 392);

    auto tk_xzzzz_xxxxxz = pbuffer.data(idx_kin_hi + 394);

    auto tk_xzzzz_xxxxyz = pbuffer.data(idx_kin_hi + 396);

    auto tk_xzzzz_xxxxzz = pbuffer.data(idx_kin_hi + 397);

    auto tk_xzzzz_xxxyyz = pbuffer.data(idx_kin_hi + 399);

    auto tk_xzzzz_xxxyzz = pbuffer.data(idx_kin_hi + 400);

    auto tk_xzzzz_xxxzzz = pbuffer.data(idx_kin_hi + 401);

    auto tk_xzzzz_xxyyyz = pbuffer.data(idx_kin_hi + 403);

    auto tk_xzzzz_xxyyzz = pbuffer.data(idx_kin_hi + 404);

    auto tk_xzzzz_xxyzzz = pbuffer.data(idx_kin_hi + 405);

    auto tk_xzzzz_xxzzzz = pbuffer.data(idx_kin_hi + 406);

    auto tk_xzzzz_xyyyyz = pbuffer.data(idx_kin_hi + 408);

    auto tk_xzzzz_xyyyzz = pbuffer.data(idx_kin_hi + 409);

    auto tk_xzzzz_xyyzzz = pbuffer.data(idx_kin_hi + 410);

    auto tk_xzzzz_xyzzzz = pbuffer.data(idx_kin_hi + 411);

    auto tk_xzzzz_xzzzzz = pbuffer.data(idx_kin_hi + 412);

    auto tk_xzzzz_yyyyyy = pbuffer.data(idx_kin_hi + 413);

    auto tk_xzzzz_yyyyyz = pbuffer.data(idx_kin_hi + 414);

    auto tk_xzzzz_yyyyzz = pbuffer.data(idx_kin_hi + 415);

    auto tk_xzzzz_yyyzzz = pbuffer.data(idx_kin_hi + 416);

    auto tk_xzzzz_yyzzzz = pbuffer.data(idx_kin_hi + 417);

    auto tk_xzzzz_yzzzzz = pbuffer.data(idx_kin_hi + 418);

    auto tk_xzzzz_zzzzzz = pbuffer.data(idx_kin_hi + 419);

    auto tk_yyyyy_xxxxxx = pbuffer.data(idx_kin_hi + 420);

    auto tk_yyyyy_xxxxxy = pbuffer.data(idx_kin_hi + 421);

    auto tk_yyyyy_xxxxxz = pbuffer.data(idx_kin_hi + 422);

    auto tk_yyyyy_xxxxyy = pbuffer.data(idx_kin_hi + 423);

    auto tk_yyyyy_xxxxyz = pbuffer.data(idx_kin_hi + 424);

    auto tk_yyyyy_xxxxzz = pbuffer.data(idx_kin_hi + 425);

    auto tk_yyyyy_xxxyyy = pbuffer.data(idx_kin_hi + 426);

    auto tk_yyyyy_xxxyyz = pbuffer.data(idx_kin_hi + 427);

    auto tk_yyyyy_xxxyzz = pbuffer.data(idx_kin_hi + 428);

    auto tk_yyyyy_xxxzzz = pbuffer.data(idx_kin_hi + 429);

    auto tk_yyyyy_xxyyyy = pbuffer.data(idx_kin_hi + 430);

    auto tk_yyyyy_xxyyyz = pbuffer.data(idx_kin_hi + 431);

    auto tk_yyyyy_xxyyzz = pbuffer.data(idx_kin_hi + 432);

    auto tk_yyyyy_xxyzzz = pbuffer.data(idx_kin_hi + 433);

    auto tk_yyyyy_xxzzzz = pbuffer.data(idx_kin_hi + 434);

    auto tk_yyyyy_xyyyyy = pbuffer.data(idx_kin_hi + 435);

    auto tk_yyyyy_xyyyyz = pbuffer.data(idx_kin_hi + 436);

    auto tk_yyyyy_xyyyzz = pbuffer.data(idx_kin_hi + 437);

    auto tk_yyyyy_xyyzzz = pbuffer.data(idx_kin_hi + 438);

    auto tk_yyyyy_xyzzzz = pbuffer.data(idx_kin_hi + 439);

    auto tk_yyyyy_xzzzzz = pbuffer.data(idx_kin_hi + 440);

    auto tk_yyyyy_yyyyyy = pbuffer.data(idx_kin_hi + 441);

    auto tk_yyyyy_yyyyyz = pbuffer.data(idx_kin_hi + 442);

    auto tk_yyyyy_yyyyzz = pbuffer.data(idx_kin_hi + 443);

    auto tk_yyyyy_yyyzzz = pbuffer.data(idx_kin_hi + 444);

    auto tk_yyyyy_yyzzzz = pbuffer.data(idx_kin_hi + 445);

    auto tk_yyyyy_yzzzzz = pbuffer.data(idx_kin_hi + 446);

    auto tk_yyyyy_zzzzzz = pbuffer.data(idx_kin_hi + 447);

    auto tk_yyyyz_xxxxxy = pbuffer.data(idx_kin_hi + 449);

    auto tk_yyyyz_xxxxxz = pbuffer.data(idx_kin_hi + 450);

    auto tk_yyyyz_xxxxyy = pbuffer.data(idx_kin_hi + 451);

    auto tk_yyyyz_xxxxyz = pbuffer.data(idx_kin_hi + 452);

    auto tk_yyyyz_xxxxzz = pbuffer.data(idx_kin_hi + 453);

    auto tk_yyyyz_xxxyyy = pbuffer.data(idx_kin_hi + 454);

    auto tk_yyyyz_xxxyyz = pbuffer.data(idx_kin_hi + 455);

    auto tk_yyyyz_xxxyzz = pbuffer.data(idx_kin_hi + 456);

    auto tk_yyyyz_xxxzzz = pbuffer.data(idx_kin_hi + 457);

    auto tk_yyyyz_xxyyyy = pbuffer.data(idx_kin_hi + 458);

    auto tk_yyyyz_xxyyyz = pbuffer.data(idx_kin_hi + 459);

    auto tk_yyyyz_xxyyzz = pbuffer.data(idx_kin_hi + 460);

    auto tk_yyyyz_xxyzzz = pbuffer.data(idx_kin_hi + 461);

    auto tk_yyyyz_xxzzzz = pbuffer.data(idx_kin_hi + 462);

    auto tk_yyyyz_xyyyyy = pbuffer.data(idx_kin_hi + 463);

    auto tk_yyyyz_xyyyyz = pbuffer.data(idx_kin_hi + 464);

    auto tk_yyyyz_xyyyzz = pbuffer.data(idx_kin_hi + 465);

    auto tk_yyyyz_xyyzzz = pbuffer.data(idx_kin_hi + 466);

    auto tk_yyyyz_xyzzzz = pbuffer.data(idx_kin_hi + 467);

    auto tk_yyyyz_xzzzzz = pbuffer.data(idx_kin_hi + 468);

    auto tk_yyyyz_yyyyyy = pbuffer.data(idx_kin_hi + 469);

    auto tk_yyyyz_yyyyyz = pbuffer.data(idx_kin_hi + 470);

    auto tk_yyyyz_yyyyzz = pbuffer.data(idx_kin_hi + 471);

    auto tk_yyyyz_yyyzzz = pbuffer.data(idx_kin_hi + 472);

    auto tk_yyyyz_yyzzzz = pbuffer.data(idx_kin_hi + 473);

    auto tk_yyyyz_yzzzzz = pbuffer.data(idx_kin_hi + 474);

    auto tk_yyyyz_zzzzzz = pbuffer.data(idx_kin_hi + 475);

    auto tk_yyyzz_xxxxxx = pbuffer.data(idx_kin_hi + 476);

    auto tk_yyyzz_xxxxxy = pbuffer.data(idx_kin_hi + 477);

    auto tk_yyyzz_xxxxxz = pbuffer.data(idx_kin_hi + 478);

    auto tk_yyyzz_xxxxyy = pbuffer.data(idx_kin_hi + 479);

    auto tk_yyyzz_xxxxyz = pbuffer.data(idx_kin_hi + 480);

    auto tk_yyyzz_xxxxzz = pbuffer.data(idx_kin_hi + 481);

    auto tk_yyyzz_xxxyyy = pbuffer.data(idx_kin_hi + 482);

    auto tk_yyyzz_xxxyyz = pbuffer.data(idx_kin_hi + 483);

    auto tk_yyyzz_xxxyzz = pbuffer.data(idx_kin_hi + 484);

    auto tk_yyyzz_xxxzzz = pbuffer.data(idx_kin_hi + 485);

    auto tk_yyyzz_xxyyyy = pbuffer.data(idx_kin_hi + 486);

    auto tk_yyyzz_xxyyyz = pbuffer.data(idx_kin_hi + 487);

    auto tk_yyyzz_xxyyzz = pbuffer.data(idx_kin_hi + 488);

    auto tk_yyyzz_xxyzzz = pbuffer.data(idx_kin_hi + 489);

    auto tk_yyyzz_xxzzzz = pbuffer.data(idx_kin_hi + 490);

    auto tk_yyyzz_xyyyyy = pbuffer.data(idx_kin_hi + 491);

    auto tk_yyyzz_xyyyyz = pbuffer.data(idx_kin_hi + 492);

    auto tk_yyyzz_xyyyzz = pbuffer.data(idx_kin_hi + 493);

    auto tk_yyyzz_xyyzzz = pbuffer.data(idx_kin_hi + 494);

    auto tk_yyyzz_xyzzzz = pbuffer.data(idx_kin_hi + 495);

    auto tk_yyyzz_xzzzzz = pbuffer.data(idx_kin_hi + 496);

    auto tk_yyyzz_yyyyyy = pbuffer.data(idx_kin_hi + 497);

    auto tk_yyyzz_yyyyyz = pbuffer.data(idx_kin_hi + 498);

    auto tk_yyyzz_yyyyzz = pbuffer.data(idx_kin_hi + 499);

    auto tk_yyyzz_yyyzzz = pbuffer.data(idx_kin_hi + 500);

    auto tk_yyyzz_yyzzzz = pbuffer.data(idx_kin_hi + 501);

    auto tk_yyyzz_yzzzzz = pbuffer.data(idx_kin_hi + 502);

    auto tk_yyyzz_zzzzzz = pbuffer.data(idx_kin_hi + 503);

    auto tk_yyzzz_xxxxxx = pbuffer.data(idx_kin_hi + 504);

    auto tk_yyzzz_xxxxxy = pbuffer.data(idx_kin_hi + 505);

    auto tk_yyzzz_xxxxxz = pbuffer.data(idx_kin_hi + 506);

    auto tk_yyzzz_xxxxyy = pbuffer.data(idx_kin_hi + 507);

    auto tk_yyzzz_xxxxyz = pbuffer.data(idx_kin_hi + 508);

    auto tk_yyzzz_xxxxzz = pbuffer.data(idx_kin_hi + 509);

    auto tk_yyzzz_xxxyyy = pbuffer.data(idx_kin_hi + 510);

    auto tk_yyzzz_xxxyyz = pbuffer.data(idx_kin_hi + 511);

    auto tk_yyzzz_xxxyzz = pbuffer.data(idx_kin_hi + 512);

    auto tk_yyzzz_xxxzzz = pbuffer.data(idx_kin_hi + 513);

    auto tk_yyzzz_xxyyyy = pbuffer.data(idx_kin_hi + 514);

    auto tk_yyzzz_xxyyyz = pbuffer.data(idx_kin_hi + 515);

    auto tk_yyzzz_xxyyzz = pbuffer.data(idx_kin_hi + 516);

    auto tk_yyzzz_xxyzzz = pbuffer.data(idx_kin_hi + 517);

    auto tk_yyzzz_xxzzzz = pbuffer.data(idx_kin_hi + 518);

    auto tk_yyzzz_xyyyyy = pbuffer.data(idx_kin_hi + 519);

    auto tk_yyzzz_xyyyyz = pbuffer.data(idx_kin_hi + 520);

    auto tk_yyzzz_xyyyzz = pbuffer.data(idx_kin_hi + 521);

    auto tk_yyzzz_xyyzzz = pbuffer.data(idx_kin_hi + 522);

    auto tk_yyzzz_xyzzzz = pbuffer.data(idx_kin_hi + 523);

    auto tk_yyzzz_xzzzzz = pbuffer.data(idx_kin_hi + 524);

    auto tk_yyzzz_yyyyyy = pbuffer.data(idx_kin_hi + 525);

    auto tk_yyzzz_yyyyyz = pbuffer.data(idx_kin_hi + 526);

    auto tk_yyzzz_yyyyzz = pbuffer.data(idx_kin_hi + 527);

    auto tk_yyzzz_yyyzzz = pbuffer.data(idx_kin_hi + 528);

    auto tk_yyzzz_yyzzzz = pbuffer.data(idx_kin_hi + 529);

    auto tk_yyzzz_yzzzzz = pbuffer.data(idx_kin_hi + 530);

    auto tk_yyzzz_zzzzzz = pbuffer.data(idx_kin_hi + 531);

    auto tk_yzzzz_xxxxxx = pbuffer.data(idx_kin_hi + 532);

    auto tk_yzzzz_xxxxxy = pbuffer.data(idx_kin_hi + 533);

    auto tk_yzzzz_xxxxxz = pbuffer.data(idx_kin_hi + 534);

    auto tk_yzzzz_xxxxyy = pbuffer.data(idx_kin_hi + 535);

    auto tk_yzzzz_xxxxyz = pbuffer.data(idx_kin_hi + 536);

    auto tk_yzzzz_xxxxzz = pbuffer.data(idx_kin_hi + 537);

    auto tk_yzzzz_xxxyyy = pbuffer.data(idx_kin_hi + 538);

    auto tk_yzzzz_xxxyyz = pbuffer.data(idx_kin_hi + 539);

    auto tk_yzzzz_xxxyzz = pbuffer.data(idx_kin_hi + 540);

    auto tk_yzzzz_xxxzzz = pbuffer.data(idx_kin_hi + 541);

    auto tk_yzzzz_xxyyyy = pbuffer.data(idx_kin_hi + 542);

    auto tk_yzzzz_xxyyyz = pbuffer.data(idx_kin_hi + 543);

    auto tk_yzzzz_xxyyzz = pbuffer.data(idx_kin_hi + 544);

    auto tk_yzzzz_xxyzzz = pbuffer.data(idx_kin_hi + 545);

    auto tk_yzzzz_xxzzzz = pbuffer.data(idx_kin_hi + 546);

    auto tk_yzzzz_xyyyyy = pbuffer.data(idx_kin_hi + 547);

    auto tk_yzzzz_xyyyyz = pbuffer.data(idx_kin_hi + 548);

    auto tk_yzzzz_xyyyzz = pbuffer.data(idx_kin_hi + 549);

    auto tk_yzzzz_xyyzzz = pbuffer.data(idx_kin_hi + 550);

    auto tk_yzzzz_xyzzzz = pbuffer.data(idx_kin_hi + 551);

    auto tk_yzzzz_xzzzzz = pbuffer.data(idx_kin_hi + 552);

    auto tk_yzzzz_yyyyyy = pbuffer.data(idx_kin_hi + 553);

    auto tk_yzzzz_yyyyyz = pbuffer.data(idx_kin_hi + 554);

    auto tk_yzzzz_yyyyzz = pbuffer.data(idx_kin_hi + 555);

    auto tk_yzzzz_yyyzzz = pbuffer.data(idx_kin_hi + 556);

    auto tk_yzzzz_yyzzzz = pbuffer.data(idx_kin_hi + 557);

    auto tk_yzzzz_yzzzzz = pbuffer.data(idx_kin_hi + 558);

    auto tk_yzzzz_zzzzzz = pbuffer.data(idx_kin_hi + 559);

    auto tk_zzzzz_xxxxxx = pbuffer.data(idx_kin_hi + 560);

    auto tk_zzzzz_xxxxxy = pbuffer.data(idx_kin_hi + 561);

    auto tk_zzzzz_xxxxxz = pbuffer.data(idx_kin_hi + 562);

    auto tk_zzzzz_xxxxyy = pbuffer.data(idx_kin_hi + 563);

    auto tk_zzzzz_xxxxyz = pbuffer.data(idx_kin_hi + 564);

    auto tk_zzzzz_xxxxzz = pbuffer.data(idx_kin_hi + 565);

    auto tk_zzzzz_xxxyyy = pbuffer.data(idx_kin_hi + 566);

    auto tk_zzzzz_xxxyyz = pbuffer.data(idx_kin_hi + 567);

    auto tk_zzzzz_xxxyzz = pbuffer.data(idx_kin_hi + 568);

    auto tk_zzzzz_xxxzzz = pbuffer.data(idx_kin_hi + 569);

    auto tk_zzzzz_xxyyyy = pbuffer.data(idx_kin_hi + 570);

    auto tk_zzzzz_xxyyyz = pbuffer.data(idx_kin_hi + 571);

    auto tk_zzzzz_xxyyzz = pbuffer.data(idx_kin_hi + 572);

    auto tk_zzzzz_xxyzzz = pbuffer.data(idx_kin_hi + 573);

    auto tk_zzzzz_xxzzzz = pbuffer.data(idx_kin_hi + 574);

    auto tk_zzzzz_xyyyyy = pbuffer.data(idx_kin_hi + 575);

    auto tk_zzzzz_xyyyyz = pbuffer.data(idx_kin_hi + 576);

    auto tk_zzzzz_xyyyzz = pbuffer.data(idx_kin_hi + 577);

    auto tk_zzzzz_xyyzzz = pbuffer.data(idx_kin_hi + 578);

    auto tk_zzzzz_xyzzzz = pbuffer.data(idx_kin_hi + 579);

    auto tk_zzzzz_xzzzzz = pbuffer.data(idx_kin_hi + 580);

    auto tk_zzzzz_yyyyyy = pbuffer.data(idx_kin_hi + 581);

    auto tk_zzzzz_yyyyyz = pbuffer.data(idx_kin_hi + 582);

    auto tk_zzzzz_yyyyzz = pbuffer.data(idx_kin_hi + 583);

    auto tk_zzzzz_yyyzzz = pbuffer.data(idx_kin_hi + 584);

    auto tk_zzzzz_yyzzzz = pbuffer.data(idx_kin_hi + 585);

    auto tk_zzzzz_yzzzzz = pbuffer.data(idx_kin_hi + 586);

    auto tk_zzzzz_zzzzzz = pbuffer.data(idx_kin_hi + 587);

    // Set up components of auxiliary buffer : II

    auto ts_xxxxxx_xxxxxx = pbuffer.data(idx_ovl_ii);

    auto ts_xxxxxx_xxxxxy = pbuffer.data(idx_ovl_ii + 1);

    auto ts_xxxxxx_xxxxxz = pbuffer.data(idx_ovl_ii + 2);

    auto ts_xxxxxx_xxxxyy = pbuffer.data(idx_ovl_ii + 3);

    auto ts_xxxxxx_xxxxyz = pbuffer.data(idx_ovl_ii + 4);

    auto ts_xxxxxx_xxxxzz = pbuffer.data(idx_ovl_ii + 5);

    auto ts_xxxxxx_xxxyyy = pbuffer.data(idx_ovl_ii + 6);

    auto ts_xxxxxx_xxxyyz = pbuffer.data(idx_ovl_ii + 7);

    auto ts_xxxxxx_xxxyzz = pbuffer.data(idx_ovl_ii + 8);

    auto ts_xxxxxx_xxxzzz = pbuffer.data(idx_ovl_ii + 9);

    auto ts_xxxxxx_xxyyyy = pbuffer.data(idx_ovl_ii + 10);

    auto ts_xxxxxx_xxyyyz = pbuffer.data(idx_ovl_ii + 11);

    auto ts_xxxxxx_xxyyzz = pbuffer.data(idx_ovl_ii + 12);

    auto ts_xxxxxx_xxyzzz = pbuffer.data(idx_ovl_ii + 13);

    auto ts_xxxxxx_xxzzzz = pbuffer.data(idx_ovl_ii + 14);

    auto ts_xxxxxx_xyyyyy = pbuffer.data(idx_ovl_ii + 15);

    auto ts_xxxxxx_xyyyyz = pbuffer.data(idx_ovl_ii + 16);

    auto ts_xxxxxx_xyyyzz = pbuffer.data(idx_ovl_ii + 17);

    auto ts_xxxxxx_xyyzzz = pbuffer.data(idx_ovl_ii + 18);

    auto ts_xxxxxx_xyzzzz = pbuffer.data(idx_ovl_ii + 19);

    auto ts_xxxxxx_xzzzzz = pbuffer.data(idx_ovl_ii + 20);

    auto ts_xxxxxx_yyyyyy = pbuffer.data(idx_ovl_ii + 21);

    auto ts_xxxxxx_yyyyyz = pbuffer.data(idx_ovl_ii + 22);

    auto ts_xxxxxx_yyyyzz = pbuffer.data(idx_ovl_ii + 23);

    auto ts_xxxxxx_yyyzzz = pbuffer.data(idx_ovl_ii + 24);

    auto ts_xxxxxx_yyzzzz = pbuffer.data(idx_ovl_ii + 25);

    auto ts_xxxxxx_yzzzzz = pbuffer.data(idx_ovl_ii + 26);

    auto ts_xxxxxx_zzzzzz = pbuffer.data(idx_ovl_ii + 27);

    auto ts_xxxxxy_xxxxxx = pbuffer.data(idx_ovl_ii + 28);

    auto ts_xxxxxy_xxxxxy = pbuffer.data(idx_ovl_ii + 29);

    auto ts_xxxxxy_xxxxxz = pbuffer.data(idx_ovl_ii + 30);

    auto ts_xxxxxy_xxxxyy = pbuffer.data(idx_ovl_ii + 31);

    auto ts_xxxxxy_xxxxyz = pbuffer.data(idx_ovl_ii + 32);

    auto ts_xxxxxy_xxxxzz = pbuffer.data(idx_ovl_ii + 33);

    auto ts_xxxxxy_xxxyyy = pbuffer.data(idx_ovl_ii + 34);

    auto ts_xxxxxy_xxxyyz = pbuffer.data(idx_ovl_ii + 35);

    auto ts_xxxxxy_xxxyzz = pbuffer.data(idx_ovl_ii + 36);

    auto ts_xxxxxy_xxxzzz = pbuffer.data(idx_ovl_ii + 37);

    auto ts_xxxxxy_xxyyyy = pbuffer.data(idx_ovl_ii + 38);

    auto ts_xxxxxy_xxyyyz = pbuffer.data(idx_ovl_ii + 39);

    auto ts_xxxxxy_xxyyzz = pbuffer.data(idx_ovl_ii + 40);

    auto ts_xxxxxy_xxyzzz = pbuffer.data(idx_ovl_ii + 41);

    auto ts_xxxxxy_xxzzzz = pbuffer.data(idx_ovl_ii + 42);

    auto ts_xxxxxy_xyyyyy = pbuffer.data(idx_ovl_ii + 43);

    auto ts_xxxxxy_xyyyyz = pbuffer.data(idx_ovl_ii + 44);

    auto ts_xxxxxy_xyyyzz = pbuffer.data(idx_ovl_ii + 45);

    auto ts_xxxxxy_xyyzzz = pbuffer.data(idx_ovl_ii + 46);

    auto ts_xxxxxy_xyzzzz = pbuffer.data(idx_ovl_ii + 47);

    auto ts_xxxxxy_xzzzzz = pbuffer.data(idx_ovl_ii + 48);

    auto ts_xxxxxy_yyyyyy = pbuffer.data(idx_ovl_ii + 49);

    auto ts_xxxxxy_yyyyyz = pbuffer.data(idx_ovl_ii + 50);

    auto ts_xxxxxy_yyyyzz = pbuffer.data(idx_ovl_ii + 51);

    auto ts_xxxxxy_yyyzzz = pbuffer.data(idx_ovl_ii + 52);

    auto ts_xxxxxy_yyzzzz = pbuffer.data(idx_ovl_ii + 53);

    auto ts_xxxxxy_yzzzzz = pbuffer.data(idx_ovl_ii + 54);

    auto ts_xxxxxy_zzzzzz = pbuffer.data(idx_ovl_ii + 55);

    auto ts_xxxxxz_xxxxxx = pbuffer.data(idx_ovl_ii + 56);

    auto ts_xxxxxz_xxxxxy = pbuffer.data(idx_ovl_ii + 57);

    auto ts_xxxxxz_xxxxxz = pbuffer.data(idx_ovl_ii + 58);

    auto ts_xxxxxz_xxxxyy = pbuffer.data(idx_ovl_ii + 59);

    auto ts_xxxxxz_xxxxyz = pbuffer.data(idx_ovl_ii + 60);

    auto ts_xxxxxz_xxxxzz = pbuffer.data(idx_ovl_ii + 61);

    auto ts_xxxxxz_xxxyyy = pbuffer.data(idx_ovl_ii + 62);

    auto ts_xxxxxz_xxxyyz = pbuffer.data(idx_ovl_ii + 63);

    auto ts_xxxxxz_xxxyzz = pbuffer.data(idx_ovl_ii + 64);

    auto ts_xxxxxz_xxxzzz = pbuffer.data(idx_ovl_ii + 65);

    auto ts_xxxxxz_xxyyyy = pbuffer.data(idx_ovl_ii + 66);

    auto ts_xxxxxz_xxyyyz = pbuffer.data(idx_ovl_ii + 67);

    auto ts_xxxxxz_xxyyzz = pbuffer.data(idx_ovl_ii + 68);

    auto ts_xxxxxz_xxyzzz = pbuffer.data(idx_ovl_ii + 69);

    auto ts_xxxxxz_xxzzzz = pbuffer.data(idx_ovl_ii + 70);

    auto ts_xxxxxz_xyyyyy = pbuffer.data(idx_ovl_ii + 71);

    auto ts_xxxxxz_xyyyyz = pbuffer.data(idx_ovl_ii + 72);

    auto ts_xxxxxz_xyyyzz = pbuffer.data(idx_ovl_ii + 73);

    auto ts_xxxxxz_xyyzzz = pbuffer.data(idx_ovl_ii + 74);

    auto ts_xxxxxz_xyzzzz = pbuffer.data(idx_ovl_ii + 75);

    auto ts_xxxxxz_xzzzzz = pbuffer.data(idx_ovl_ii + 76);

    auto ts_xxxxxz_yyyyyy = pbuffer.data(idx_ovl_ii + 77);

    auto ts_xxxxxz_yyyyyz = pbuffer.data(idx_ovl_ii + 78);

    auto ts_xxxxxz_yyyyzz = pbuffer.data(idx_ovl_ii + 79);

    auto ts_xxxxxz_yyyzzz = pbuffer.data(idx_ovl_ii + 80);

    auto ts_xxxxxz_yyzzzz = pbuffer.data(idx_ovl_ii + 81);

    auto ts_xxxxxz_yzzzzz = pbuffer.data(idx_ovl_ii + 82);

    auto ts_xxxxxz_zzzzzz = pbuffer.data(idx_ovl_ii + 83);

    auto ts_xxxxyy_xxxxxx = pbuffer.data(idx_ovl_ii + 84);

    auto ts_xxxxyy_xxxxxy = pbuffer.data(idx_ovl_ii + 85);

    auto ts_xxxxyy_xxxxxz = pbuffer.data(idx_ovl_ii + 86);

    auto ts_xxxxyy_xxxxyy = pbuffer.data(idx_ovl_ii + 87);

    auto ts_xxxxyy_xxxxyz = pbuffer.data(idx_ovl_ii + 88);

    auto ts_xxxxyy_xxxxzz = pbuffer.data(idx_ovl_ii + 89);

    auto ts_xxxxyy_xxxyyy = pbuffer.data(idx_ovl_ii + 90);

    auto ts_xxxxyy_xxxyyz = pbuffer.data(idx_ovl_ii + 91);

    auto ts_xxxxyy_xxxyzz = pbuffer.data(idx_ovl_ii + 92);

    auto ts_xxxxyy_xxxzzz = pbuffer.data(idx_ovl_ii + 93);

    auto ts_xxxxyy_xxyyyy = pbuffer.data(idx_ovl_ii + 94);

    auto ts_xxxxyy_xxyyyz = pbuffer.data(idx_ovl_ii + 95);

    auto ts_xxxxyy_xxyyzz = pbuffer.data(idx_ovl_ii + 96);

    auto ts_xxxxyy_xxyzzz = pbuffer.data(idx_ovl_ii + 97);

    auto ts_xxxxyy_xxzzzz = pbuffer.data(idx_ovl_ii + 98);

    auto ts_xxxxyy_xyyyyy = pbuffer.data(idx_ovl_ii + 99);

    auto ts_xxxxyy_xyyyyz = pbuffer.data(idx_ovl_ii + 100);

    auto ts_xxxxyy_xyyyzz = pbuffer.data(idx_ovl_ii + 101);

    auto ts_xxxxyy_xyyzzz = pbuffer.data(idx_ovl_ii + 102);

    auto ts_xxxxyy_xyzzzz = pbuffer.data(idx_ovl_ii + 103);

    auto ts_xxxxyy_xzzzzz = pbuffer.data(idx_ovl_ii + 104);

    auto ts_xxxxyy_yyyyyy = pbuffer.data(idx_ovl_ii + 105);

    auto ts_xxxxyy_yyyyyz = pbuffer.data(idx_ovl_ii + 106);

    auto ts_xxxxyy_yyyyzz = pbuffer.data(idx_ovl_ii + 107);

    auto ts_xxxxyy_yyyzzz = pbuffer.data(idx_ovl_ii + 108);

    auto ts_xxxxyy_yyzzzz = pbuffer.data(idx_ovl_ii + 109);

    auto ts_xxxxyy_yzzzzz = pbuffer.data(idx_ovl_ii + 110);

    auto ts_xxxxyy_zzzzzz = pbuffer.data(idx_ovl_ii + 111);

    auto ts_xxxxyz_xxxxxx = pbuffer.data(idx_ovl_ii + 112);

    auto ts_xxxxyz_xxxxxy = pbuffer.data(idx_ovl_ii + 113);

    auto ts_xxxxyz_xxxxxz = pbuffer.data(idx_ovl_ii + 114);

    auto ts_xxxxyz_xxxxyy = pbuffer.data(idx_ovl_ii + 115);

    auto ts_xxxxyz_xxxxyz = pbuffer.data(idx_ovl_ii + 116);

    auto ts_xxxxyz_xxxxzz = pbuffer.data(idx_ovl_ii + 117);

    auto ts_xxxxyz_xxxyyy = pbuffer.data(idx_ovl_ii + 118);

    auto ts_xxxxyz_xxxyyz = pbuffer.data(idx_ovl_ii + 119);

    auto ts_xxxxyz_xxxyzz = pbuffer.data(idx_ovl_ii + 120);

    auto ts_xxxxyz_xxxzzz = pbuffer.data(idx_ovl_ii + 121);

    auto ts_xxxxyz_xxyyyy = pbuffer.data(idx_ovl_ii + 122);

    auto ts_xxxxyz_xxyyyz = pbuffer.data(idx_ovl_ii + 123);

    auto ts_xxxxyz_xxyyzz = pbuffer.data(idx_ovl_ii + 124);

    auto ts_xxxxyz_xxyzzz = pbuffer.data(idx_ovl_ii + 125);

    auto ts_xxxxyz_xxzzzz = pbuffer.data(idx_ovl_ii + 126);

    auto ts_xxxxyz_xyyyyy = pbuffer.data(idx_ovl_ii + 127);

    auto ts_xxxxyz_xyyyyz = pbuffer.data(idx_ovl_ii + 128);

    auto ts_xxxxyz_xyyyzz = pbuffer.data(idx_ovl_ii + 129);

    auto ts_xxxxyz_xyyzzz = pbuffer.data(idx_ovl_ii + 130);

    auto ts_xxxxyz_xyzzzz = pbuffer.data(idx_ovl_ii + 131);

    auto ts_xxxxyz_xzzzzz = pbuffer.data(idx_ovl_ii + 132);

    auto ts_xxxxyz_yyyyyy = pbuffer.data(idx_ovl_ii + 133);

    auto ts_xxxxyz_yyyyyz = pbuffer.data(idx_ovl_ii + 134);

    auto ts_xxxxyz_yyyyzz = pbuffer.data(idx_ovl_ii + 135);

    auto ts_xxxxyz_yyyzzz = pbuffer.data(idx_ovl_ii + 136);

    auto ts_xxxxyz_yyzzzz = pbuffer.data(idx_ovl_ii + 137);

    auto ts_xxxxyz_yzzzzz = pbuffer.data(idx_ovl_ii + 138);

    auto ts_xxxxyz_zzzzzz = pbuffer.data(idx_ovl_ii + 139);

    auto ts_xxxxzz_xxxxxx = pbuffer.data(idx_ovl_ii + 140);

    auto ts_xxxxzz_xxxxxy = pbuffer.data(idx_ovl_ii + 141);

    auto ts_xxxxzz_xxxxxz = pbuffer.data(idx_ovl_ii + 142);

    auto ts_xxxxzz_xxxxyy = pbuffer.data(idx_ovl_ii + 143);

    auto ts_xxxxzz_xxxxyz = pbuffer.data(idx_ovl_ii + 144);

    auto ts_xxxxzz_xxxxzz = pbuffer.data(idx_ovl_ii + 145);

    auto ts_xxxxzz_xxxyyy = pbuffer.data(idx_ovl_ii + 146);

    auto ts_xxxxzz_xxxyyz = pbuffer.data(idx_ovl_ii + 147);

    auto ts_xxxxzz_xxxyzz = pbuffer.data(idx_ovl_ii + 148);

    auto ts_xxxxzz_xxxzzz = pbuffer.data(idx_ovl_ii + 149);

    auto ts_xxxxzz_xxyyyy = pbuffer.data(idx_ovl_ii + 150);

    auto ts_xxxxzz_xxyyyz = pbuffer.data(idx_ovl_ii + 151);

    auto ts_xxxxzz_xxyyzz = pbuffer.data(idx_ovl_ii + 152);

    auto ts_xxxxzz_xxyzzz = pbuffer.data(idx_ovl_ii + 153);

    auto ts_xxxxzz_xxzzzz = pbuffer.data(idx_ovl_ii + 154);

    auto ts_xxxxzz_xyyyyy = pbuffer.data(idx_ovl_ii + 155);

    auto ts_xxxxzz_xyyyyz = pbuffer.data(idx_ovl_ii + 156);

    auto ts_xxxxzz_xyyyzz = pbuffer.data(idx_ovl_ii + 157);

    auto ts_xxxxzz_xyyzzz = pbuffer.data(idx_ovl_ii + 158);

    auto ts_xxxxzz_xyzzzz = pbuffer.data(idx_ovl_ii + 159);

    auto ts_xxxxzz_xzzzzz = pbuffer.data(idx_ovl_ii + 160);

    auto ts_xxxxzz_yyyyyy = pbuffer.data(idx_ovl_ii + 161);

    auto ts_xxxxzz_yyyyyz = pbuffer.data(idx_ovl_ii + 162);

    auto ts_xxxxzz_yyyyzz = pbuffer.data(idx_ovl_ii + 163);

    auto ts_xxxxzz_yyyzzz = pbuffer.data(idx_ovl_ii + 164);

    auto ts_xxxxzz_yyzzzz = pbuffer.data(idx_ovl_ii + 165);

    auto ts_xxxxzz_yzzzzz = pbuffer.data(idx_ovl_ii + 166);

    auto ts_xxxxzz_zzzzzz = pbuffer.data(idx_ovl_ii + 167);

    auto ts_xxxyyy_xxxxxx = pbuffer.data(idx_ovl_ii + 168);

    auto ts_xxxyyy_xxxxxy = pbuffer.data(idx_ovl_ii + 169);

    auto ts_xxxyyy_xxxxxz = pbuffer.data(idx_ovl_ii + 170);

    auto ts_xxxyyy_xxxxyy = pbuffer.data(idx_ovl_ii + 171);

    auto ts_xxxyyy_xxxxyz = pbuffer.data(idx_ovl_ii + 172);

    auto ts_xxxyyy_xxxxzz = pbuffer.data(idx_ovl_ii + 173);

    auto ts_xxxyyy_xxxyyy = pbuffer.data(idx_ovl_ii + 174);

    auto ts_xxxyyy_xxxyyz = pbuffer.data(idx_ovl_ii + 175);

    auto ts_xxxyyy_xxxyzz = pbuffer.data(idx_ovl_ii + 176);

    auto ts_xxxyyy_xxxzzz = pbuffer.data(idx_ovl_ii + 177);

    auto ts_xxxyyy_xxyyyy = pbuffer.data(idx_ovl_ii + 178);

    auto ts_xxxyyy_xxyyyz = pbuffer.data(idx_ovl_ii + 179);

    auto ts_xxxyyy_xxyyzz = pbuffer.data(idx_ovl_ii + 180);

    auto ts_xxxyyy_xxyzzz = pbuffer.data(idx_ovl_ii + 181);

    auto ts_xxxyyy_xxzzzz = pbuffer.data(idx_ovl_ii + 182);

    auto ts_xxxyyy_xyyyyy = pbuffer.data(idx_ovl_ii + 183);

    auto ts_xxxyyy_xyyyyz = pbuffer.data(idx_ovl_ii + 184);

    auto ts_xxxyyy_xyyyzz = pbuffer.data(idx_ovl_ii + 185);

    auto ts_xxxyyy_xyyzzz = pbuffer.data(idx_ovl_ii + 186);

    auto ts_xxxyyy_xyzzzz = pbuffer.data(idx_ovl_ii + 187);

    auto ts_xxxyyy_xzzzzz = pbuffer.data(idx_ovl_ii + 188);

    auto ts_xxxyyy_yyyyyy = pbuffer.data(idx_ovl_ii + 189);

    auto ts_xxxyyy_yyyyyz = pbuffer.data(idx_ovl_ii + 190);

    auto ts_xxxyyy_yyyyzz = pbuffer.data(idx_ovl_ii + 191);

    auto ts_xxxyyy_yyyzzz = pbuffer.data(idx_ovl_ii + 192);

    auto ts_xxxyyy_yyzzzz = pbuffer.data(idx_ovl_ii + 193);

    auto ts_xxxyyy_yzzzzz = pbuffer.data(idx_ovl_ii + 194);

    auto ts_xxxyyy_zzzzzz = pbuffer.data(idx_ovl_ii + 195);

    auto ts_xxxyyz_xxxxxx = pbuffer.data(idx_ovl_ii + 196);

    auto ts_xxxyyz_xxxxxy = pbuffer.data(idx_ovl_ii + 197);

    auto ts_xxxyyz_xxxxxz = pbuffer.data(idx_ovl_ii + 198);

    auto ts_xxxyyz_xxxxyy = pbuffer.data(idx_ovl_ii + 199);

    auto ts_xxxyyz_xxxxyz = pbuffer.data(idx_ovl_ii + 200);

    auto ts_xxxyyz_xxxxzz = pbuffer.data(idx_ovl_ii + 201);

    auto ts_xxxyyz_xxxyyy = pbuffer.data(idx_ovl_ii + 202);

    auto ts_xxxyyz_xxxyyz = pbuffer.data(idx_ovl_ii + 203);

    auto ts_xxxyyz_xxxyzz = pbuffer.data(idx_ovl_ii + 204);

    auto ts_xxxyyz_xxxzzz = pbuffer.data(idx_ovl_ii + 205);

    auto ts_xxxyyz_xxyyyy = pbuffer.data(idx_ovl_ii + 206);

    auto ts_xxxyyz_xxyyyz = pbuffer.data(idx_ovl_ii + 207);

    auto ts_xxxyyz_xxyyzz = pbuffer.data(idx_ovl_ii + 208);

    auto ts_xxxyyz_xxyzzz = pbuffer.data(idx_ovl_ii + 209);

    auto ts_xxxyyz_xxzzzz = pbuffer.data(idx_ovl_ii + 210);

    auto ts_xxxyyz_xyyyyy = pbuffer.data(idx_ovl_ii + 211);

    auto ts_xxxyyz_xyyyyz = pbuffer.data(idx_ovl_ii + 212);

    auto ts_xxxyyz_xyyyzz = pbuffer.data(idx_ovl_ii + 213);

    auto ts_xxxyyz_xyyzzz = pbuffer.data(idx_ovl_ii + 214);

    auto ts_xxxyyz_xyzzzz = pbuffer.data(idx_ovl_ii + 215);

    auto ts_xxxyyz_xzzzzz = pbuffer.data(idx_ovl_ii + 216);

    auto ts_xxxyyz_yyyyyy = pbuffer.data(idx_ovl_ii + 217);

    auto ts_xxxyyz_yyyyyz = pbuffer.data(idx_ovl_ii + 218);

    auto ts_xxxyyz_yyyyzz = pbuffer.data(idx_ovl_ii + 219);

    auto ts_xxxyyz_yyyzzz = pbuffer.data(idx_ovl_ii + 220);

    auto ts_xxxyyz_yyzzzz = pbuffer.data(idx_ovl_ii + 221);

    auto ts_xxxyyz_yzzzzz = pbuffer.data(idx_ovl_ii + 222);

    auto ts_xxxyyz_zzzzzz = pbuffer.data(idx_ovl_ii + 223);

    auto ts_xxxyzz_xxxxxx = pbuffer.data(idx_ovl_ii + 224);

    auto ts_xxxyzz_xxxxxy = pbuffer.data(idx_ovl_ii + 225);

    auto ts_xxxyzz_xxxxxz = pbuffer.data(idx_ovl_ii + 226);

    auto ts_xxxyzz_xxxxyy = pbuffer.data(idx_ovl_ii + 227);

    auto ts_xxxyzz_xxxxyz = pbuffer.data(idx_ovl_ii + 228);

    auto ts_xxxyzz_xxxxzz = pbuffer.data(idx_ovl_ii + 229);

    auto ts_xxxyzz_xxxyyy = pbuffer.data(idx_ovl_ii + 230);

    auto ts_xxxyzz_xxxyyz = pbuffer.data(idx_ovl_ii + 231);

    auto ts_xxxyzz_xxxyzz = pbuffer.data(idx_ovl_ii + 232);

    auto ts_xxxyzz_xxxzzz = pbuffer.data(idx_ovl_ii + 233);

    auto ts_xxxyzz_xxyyyy = pbuffer.data(idx_ovl_ii + 234);

    auto ts_xxxyzz_xxyyyz = pbuffer.data(idx_ovl_ii + 235);

    auto ts_xxxyzz_xxyyzz = pbuffer.data(idx_ovl_ii + 236);

    auto ts_xxxyzz_xxyzzz = pbuffer.data(idx_ovl_ii + 237);

    auto ts_xxxyzz_xxzzzz = pbuffer.data(idx_ovl_ii + 238);

    auto ts_xxxyzz_xyyyyy = pbuffer.data(idx_ovl_ii + 239);

    auto ts_xxxyzz_xyyyyz = pbuffer.data(idx_ovl_ii + 240);

    auto ts_xxxyzz_xyyyzz = pbuffer.data(idx_ovl_ii + 241);

    auto ts_xxxyzz_xyyzzz = pbuffer.data(idx_ovl_ii + 242);

    auto ts_xxxyzz_xyzzzz = pbuffer.data(idx_ovl_ii + 243);

    auto ts_xxxyzz_xzzzzz = pbuffer.data(idx_ovl_ii + 244);

    auto ts_xxxyzz_yyyyyy = pbuffer.data(idx_ovl_ii + 245);

    auto ts_xxxyzz_yyyyyz = pbuffer.data(idx_ovl_ii + 246);

    auto ts_xxxyzz_yyyyzz = pbuffer.data(idx_ovl_ii + 247);

    auto ts_xxxyzz_yyyzzz = pbuffer.data(idx_ovl_ii + 248);

    auto ts_xxxyzz_yyzzzz = pbuffer.data(idx_ovl_ii + 249);

    auto ts_xxxyzz_yzzzzz = pbuffer.data(idx_ovl_ii + 250);

    auto ts_xxxyzz_zzzzzz = pbuffer.data(idx_ovl_ii + 251);

    auto ts_xxxzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 252);

    auto ts_xxxzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 253);

    auto ts_xxxzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 254);

    auto ts_xxxzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 255);

    auto ts_xxxzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 256);

    auto ts_xxxzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 257);

    auto ts_xxxzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 258);

    auto ts_xxxzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 259);

    auto ts_xxxzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 260);

    auto ts_xxxzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 261);

    auto ts_xxxzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 262);

    auto ts_xxxzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 263);

    auto ts_xxxzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 264);

    auto ts_xxxzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 265);

    auto ts_xxxzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 266);

    auto ts_xxxzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 267);

    auto ts_xxxzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 268);

    auto ts_xxxzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 269);

    auto ts_xxxzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 270);

    auto ts_xxxzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 271);

    auto ts_xxxzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 272);

    auto ts_xxxzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 273);

    auto ts_xxxzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 274);

    auto ts_xxxzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 275);

    auto ts_xxxzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 276);

    auto ts_xxxzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 277);

    auto ts_xxxzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 278);

    auto ts_xxxzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 279);

    auto ts_xxyyyy_xxxxxx = pbuffer.data(idx_ovl_ii + 280);

    auto ts_xxyyyy_xxxxxy = pbuffer.data(idx_ovl_ii + 281);

    auto ts_xxyyyy_xxxxxz = pbuffer.data(idx_ovl_ii + 282);

    auto ts_xxyyyy_xxxxyy = pbuffer.data(idx_ovl_ii + 283);

    auto ts_xxyyyy_xxxxyz = pbuffer.data(idx_ovl_ii + 284);

    auto ts_xxyyyy_xxxxzz = pbuffer.data(idx_ovl_ii + 285);

    auto ts_xxyyyy_xxxyyy = pbuffer.data(idx_ovl_ii + 286);

    auto ts_xxyyyy_xxxyyz = pbuffer.data(idx_ovl_ii + 287);

    auto ts_xxyyyy_xxxyzz = pbuffer.data(idx_ovl_ii + 288);

    auto ts_xxyyyy_xxxzzz = pbuffer.data(idx_ovl_ii + 289);

    auto ts_xxyyyy_xxyyyy = pbuffer.data(idx_ovl_ii + 290);

    auto ts_xxyyyy_xxyyyz = pbuffer.data(idx_ovl_ii + 291);

    auto ts_xxyyyy_xxyyzz = pbuffer.data(idx_ovl_ii + 292);

    auto ts_xxyyyy_xxyzzz = pbuffer.data(idx_ovl_ii + 293);

    auto ts_xxyyyy_xxzzzz = pbuffer.data(idx_ovl_ii + 294);

    auto ts_xxyyyy_xyyyyy = pbuffer.data(idx_ovl_ii + 295);

    auto ts_xxyyyy_xyyyyz = pbuffer.data(idx_ovl_ii + 296);

    auto ts_xxyyyy_xyyyzz = pbuffer.data(idx_ovl_ii + 297);

    auto ts_xxyyyy_xyyzzz = pbuffer.data(idx_ovl_ii + 298);

    auto ts_xxyyyy_xyzzzz = pbuffer.data(idx_ovl_ii + 299);

    auto ts_xxyyyy_xzzzzz = pbuffer.data(idx_ovl_ii + 300);

    auto ts_xxyyyy_yyyyyy = pbuffer.data(idx_ovl_ii + 301);

    auto ts_xxyyyy_yyyyyz = pbuffer.data(idx_ovl_ii + 302);

    auto ts_xxyyyy_yyyyzz = pbuffer.data(idx_ovl_ii + 303);

    auto ts_xxyyyy_yyyzzz = pbuffer.data(idx_ovl_ii + 304);

    auto ts_xxyyyy_yyzzzz = pbuffer.data(idx_ovl_ii + 305);

    auto ts_xxyyyy_yzzzzz = pbuffer.data(idx_ovl_ii + 306);

    auto ts_xxyyyy_zzzzzz = pbuffer.data(idx_ovl_ii + 307);

    auto ts_xxyyyz_xxxxxx = pbuffer.data(idx_ovl_ii + 308);

    auto ts_xxyyyz_xxxxxy = pbuffer.data(idx_ovl_ii + 309);

    auto ts_xxyyyz_xxxxxz = pbuffer.data(idx_ovl_ii + 310);

    auto ts_xxyyyz_xxxxyy = pbuffer.data(idx_ovl_ii + 311);

    auto ts_xxyyyz_xxxxyz = pbuffer.data(idx_ovl_ii + 312);

    auto ts_xxyyyz_xxxxzz = pbuffer.data(idx_ovl_ii + 313);

    auto ts_xxyyyz_xxxyyy = pbuffer.data(idx_ovl_ii + 314);

    auto ts_xxyyyz_xxxyyz = pbuffer.data(idx_ovl_ii + 315);

    auto ts_xxyyyz_xxxyzz = pbuffer.data(idx_ovl_ii + 316);

    auto ts_xxyyyz_xxxzzz = pbuffer.data(idx_ovl_ii + 317);

    auto ts_xxyyyz_xxyyyy = pbuffer.data(idx_ovl_ii + 318);

    auto ts_xxyyyz_xxyyyz = pbuffer.data(idx_ovl_ii + 319);

    auto ts_xxyyyz_xxyyzz = pbuffer.data(idx_ovl_ii + 320);

    auto ts_xxyyyz_xxyzzz = pbuffer.data(idx_ovl_ii + 321);

    auto ts_xxyyyz_xxzzzz = pbuffer.data(idx_ovl_ii + 322);

    auto ts_xxyyyz_xyyyyy = pbuffer.data(idx_ovl_ii + 323);

    auto ts_xxyyyz_xyyyyz = pbuffer.data(idx_ovl_ii + 324);

    auto ts_xxyyyz_xyyyzz = pbuffer.data(idx_ovl_ii + 325);

    auto ts_xxyyyz_xyyzzz = pbuffer.data(idx_ovl_ii + 326);

    auto ts_xxyyyz_xyzzzz = pbuffer.data(idx_ovl_ii + 327);

    auto ts_xxyyyz_xzzzzz = pbuffer.data(idx_ovl_ii + 328);

    auto ts_xxyyyz_yyyyyy = pbuffer.data(idx_ovl_ii + 329);

    auto ts_xxyyyz_yyyyyz = pbuffer.data(idx_ovl_ii + 330);

    auto ts_xxyyyz_yyyyzz = pbuffer.data(idx_ovl_ii + 331);

    auto ts_xxyyyz_yyyzzz = pbuffer.data(idx_ovl_ii + 332);

    auto ts_xxyyyz_yyzzzz = pbuffer.data(idx_ovl_ii + 333);

    auto ts_xxyyyz_yzzzzz = pbuffer.data(idx_ovl_ii + 334);

    auto ts_xxyyyz_zzzzzz = pbuffer.data(idx_ovl_ii + 335);

    auto ts_xxyyzz_xxxxxx = pbuffer.data(idx_ovl_ii + 336);

    auto ts_xxyyzz_xxxxxy = pbuffer.data(idx_ovl_ii + 337);

    auto ts_xxyyzz_xxxxxz = pbuffer.data(idx_ovl_ii + 338);

    auto ts_xxyyzz_xxxxyy = pbuffer.data(idx_ovl_ii + 339);

    auto ts_xxyyzz_xxxxyz = pbuffer.data(idx_ovl_ii + 340);

    auto ts_xxyyzz_xxxxzz = pbuffer.data(idx_ovl_ii + 341);

    auto ts_xxyyzz_xxxyyy = pbuffer.data(idx_ovl_ii + 342);

    auto ts_xxyyzz_xxxyyz = pbuffer.data(idx_ovl_ii + 343);

    auto ts_xxyyzz_xxxyzz = pbuffer.data(idx_ovl_ii + 344);

    auto ts_xxyyzz_xxxzzz = pbuffer.data(idx_ovl_ii + 345);

    auto ts_xxyyzz_xxyyyy = pbuffer.data(idx_ovl_ii + 346);

    auto ts_xxyyzz_xxyyyz = pbuffer.data(idx_ovl_ii + 347);

    auto ts_xxyyzz_xxyyzz = pbuffer.data(idx_ovl_ii + 348);

    auto ts_xxyyzz_xxyzzz = pbuffer.data(idx_ovl_ii + 349);

    auto ts_xxyyzz_xxzzzz = pbuffer.data(idx_ovl_ii + 350);

    auto ts_xxyyzz_xyyyyy = pbuffer.data(idx_ovl_ii + 351);

    auto ts_xxyyzz_xyyyyz = pbuffer.data(idx_ovl_ii + 352);

    auto ts_xxyyzz_xyyyzz = pbuffer.data(idx_ovl_ii + 353);

    auto ts_xxyyzz_xyyzzz = pbuffer.data(idx_ovl_ii + 354);

    auto ts_xxyyzz_xyzzzz = pbuffer.data(idx_ovl_ii + 355);

    auto ts_xxyyzz_xzzzzz = pbuffer.data(idx_ovl_ii + 356);

    auto ts_xxyyzz_yyyyyy = pbuffer.data(idx_ovl_ii + 357);

    auto ts_xxyyzz_yyyyyz = pbuffer.data(idx_ovl_ii + 358);

    auto ts_xxyyzz_yyyyzz = pbuffer.data(idx_ovl_ii + 359);

    auto ts_xxyyzz_yyyzzz = pbuffer.data(idx_ovl_ii + 360);

    auto ts_xxyyzz_yyzzzz = pbuffer.data(idx_ovl_ii + 361);

    auto ts_xxyyzz_yzzzzz = pbuffer.data(idx_ovl_ii + 362);

    auto ts_xxyyzz_zzzzzz = pbuffer.data(idx_ovl_ii + 363);

    auto ts_xxyzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 364);

    auto ts_xxyzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 365);

    auto ts_xxyzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 366);

    auto ts_xxyzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 367);

    auto ts_xxyzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 368);

    auto ts_xxyzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 369);

    auto ts_xxyzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 370);

    auto ts_xxyzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 371);

    auto ts_xxyzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 372);

    auto ts_xxyzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 373);

    auto ts_xxyzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 374);

    auto ts_xxyzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 375);

    auto ts_xxyzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 376);

    auto ts_xxyzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 377);

    auto ts_xxyzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 378);

    auto ts_xxyzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 379);

    auto ts_xxyzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 380);

    auto ts_xxyzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 381);

    auto ts_xxyzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 382);

    auto ts_xxyzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 383);

    auto ts_xxyzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 384);

    auto ts_xxyzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 385);

    auto ts_xxyzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 386);

    auto ts_xxyzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 387);

    auto ts_xxyzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 388);

    auto ts_xxyzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 389);

    auto ts_xxyzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 390);

    auto ts_xxyzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 391);

    auto ts_xxzzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 392);

    auto ts_xxzzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 393);

    auto ts_xxzzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 394);

    auto ts_xxzzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 395);

    auto ts_xxzzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 396);

    auto ts_xxzzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 397);

    auto ts_xxzzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 398);

    auto ts_xxzzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 399);

    auto ts_xxzzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 400);

    auto ts_xxzzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 401);

    auto ts_xxzzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 402);

    auto ts_xxzzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 403);

    auto ts_xxzzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 404);

    auto ts_xxzzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 405);

    auto ts_xxzzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 406);

    auto ts_xxzzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 407);

    auto ts_xxzzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 408);

    auto ts_xxzzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 409);

    auto ts_xxzzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 410);

    auto ts_xxzzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 411);

    auto ts_xxzzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 412);

    auto ts_xxzzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 413);

    auto ts_xxzzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 414);

    auto ts_xxzzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 415);

    auto ts_xxzzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 416);

    auto ts_xxzzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 417);

    auto ts_xxzzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 418);

    auto ts_xxzzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 419);

    auto ts_xyyyyy_xxxxxx = pbuffer.data(idx_ovl_ii + 420);

    auto ts_xyyyyy_xxxxxy = pbuffer.data(idx_ovl_ii + 421);

    auto ts_xyyyyy_xxxxxz = pbuffer.data(idx_ovl_ii + 422);

    auto ts_xyyyyy_xxxxyy = pbuffer.data(idx_ovl_ii + 423);

    auto ts_xyyyyy_xxxxyz = pbuffer.data(idx_ovl_ii + 424);

    auto ts_xyyyyy_xxxxzz = pbuffer.data(idx_ovl_ii + 425);

    auto ts_xyyyyy_xxxyyy = pbuffer.data(idx_ovl_ii + 426);

    auto ts_xyyyyy_xxxyyz = pbuffer.data(idx_ovl_ii + 427);

    auto ts_xyyyyy_xxxyzz = pbuffer.data(idx_ovl_ii + 428);

    auto ts_xyyyyy_xxxzzz = pbuffer.data(idx_ovl_ii + 429);

    auto ts_xyyyyy_xxyyyy = pbuffer.data(idx_ovl_ii + 430);

    auto ts_xyyyyy_xxyyyz = pbuffer.data(idx_ovl_ii + 431);

    auto ts_xyyyyy_xxyyzz = pbuffer.data(idx_ovl_ii + 432);

    auto ts_xyyyyy_xxyzzz = pbuffer.data(idx_ovl_ii + 433);

    auto ts_xyyyyy_xxzzzz = pbuffer.data(idx_ovl_ii + 434);

    auto ts_xyyyyy_xyyyyy = pbuffer.data(idx_ovl_ii + 435);

    auto ts_xyyyyy_xyyyyz = pbuffer.data(idx_ovl_ii + 436);

    auto ts_xyyyyy_xyyyzz = pbuffer.data(idx_ovl_ii + 437);

    auto ts_xyyyyy_xyyzzz = pbuffer.data(idx_ovl_ii + 438);

    auto ts_xyyyyy_xyzzzz = pbuffer.data(idx_ovl_ii + 439);

    auto ts_xyyyyy_xzzzzz = pbuffer.data(idx_ovl_ii + 440);

    auto ts_xyyyyy_yyyyyy = pbuffer.data(idx_ovl_ii + 441);

    auto ts_xyyyyy_yyyyyz = pbuffer.data(idx_ovl_ii + 442);

    auto ts_xyyyyy_yyyyzz = pbuffer.data(idx_ovl_ii + 443);

    auto ts_xyyyyy_yyyzzz = pbuffer.data(idx_ovl_ii + 444);

    auto ts_xyyyyy_yyzzzz = pbuffer.data(idx_ovl_ii + 445);

    auto ts_xyyyyy_yzzzzz = pbuffer.data(idx_ovl_ii + 446);

    auto ts_xyyyyy_zzzzzz = pbuffer.data(idx_ovl_ii + 447);

    auto ts_xyyyyz_xxxxxx = pbuffer.data(idx_ovl_ii + 448);

    auto ts_xyyyyz_xxxxxy = pbuffer.data(idx_ovl_ii + 449);

    auto ts_xyyyyz_xxxxxz = pbuffer.data(idx_ovl_ii + 450);

    auto ts_xyyyyz_xxxxyy = pbuffer.data(idx_ovl_ii + 451);

    auto ts_xyyyyz_xxxxyz = pbuffer.data(idx_ovl_ii + 452);

    auto ts_xyyyyz_xxxxzz = pbuffer.data(idx_ovl_ii + 453);

    auto ts_xyyyyz_xxxyyy = pbuffer.data(idx_ovl_ii + 454);

    auto ts_xyyyyz_xxxyyz = pbuffer.data(idx_ovl_ii + 455);

    auto ts_xyyyyz_xxxyzz = pbuffer.data(idx_ovl_ii + 456);

    auto ts_xyyyyz_xxxzzz = pbuffer.data(idx_ovl_ii + 457);

    auto ts_xyyyyz_xxyyyy = pbuffer.data(idx_ovl_ii + 458);

    auto ts_xyyyyz_xxyyyz = pbuffer.data(idx_ovl_ii + 459);

    auto ts_xyyyyz_xxyyzz = pbuffer.data(idx_ovl_ii + 460);

    auto ts_xyyyyz_xxyzzz = pbuffer.data(idx_ovl_ii + 461);

    auto ts_xyyyyz_xxzzzz = pbuffer.data(idx_ovl_ii + 462);

    auto ts_xyyyyz_xyyyyy = pbuffer.data(idx_ovl_ii + 463);

    auto ts_xyyyyz_xyyyyz = pbuffer.data(idx_ovl_ii + 464);

    auto ts_xyyyyz_xyyyzz = pbuffer.data(idx_ovl_ii + 465);

    auto ts_xyyyyz_xyyzzz = pbuffer.data(idx_ovl_ii + 466);

    auto ts_xyyyyz_xyzzzz = pbuffer.data(idx_ovl_ii + 467);

    auto ts_xyyyyz_xzzzzz = pbuffer.data(idx_ovl_ii + 468);

    auto ts_xyyyyz_yyyyyy = pbuffer.data(idx_ovl_ii + 469);

    auto ts_xyyyyz_yyyyyz = pbuffer.data(idx_ovl_ii + 470);

    auto ts_xyyyyz_yyyyzz = pbuffer.data(idx_ovl_ii + 471);

    auto ts_xyyyyz_yyyzzz = pbuffer.data(idx_ovl_ii + 472);

    auto ts_xyyyyz_yyzzzz = pbuffer.data(idx_ovl_ii + 473);

    auto ts_xyyyyz_yzzzzz = pbuffer.data(idx_ovl_ii + 474);

    auto ts_xyyyyz_zzzzzz = pbuffer.data(idx_ovl_ii + 475);

    auto ts_xyyyzz_xxxxxx = pbuffer.data(idx_ovl_ii + 476);

    auto ts_xyyyzz_xxxxxy = pbuffer.data(idx_ovl_ii + 477);

    auto ts_xyyyzz_xxxxxz = pbuffer.data(idx_ovl_ii + 478);

    auto ts_xyyyzz_xxxxyy = pbuffer.data(idx_ovl_ii + 479);

    auto ts_xyyyzz_xxxxyz = pbuffer.data(idx_ovl_ii + 480);

    auto ts_xyyyzz_xxxxzz = pbuffer.data(idx_ovl_ii + 481);

    auto ts_xyyyzz_xxxyyy = pbuffer.data(idx_ovl_ii + 482);

    auto ts_xyyyzz_xxxyyz = pbuffer.data(idx_ovl_ii + 483);

    auto ts_xyyyzz_xxxyzz = pbuffer.data(idx_ovl_ii + 484);

    auto ts_xyyyzz_xxxzzz = pbuffer.data(idx_ovl_ii + 485);

    auto ts_xyyyzz_xxyyyy = pbuffer.data(idx_ovl_ii + 486);

    auto ts_xyyyzz_xxyyyz = pbuffer.data(idx_ovl_ii + 487);

    auto ts_xyyyzz_xxyyzz = pbuffer.data(idx_ovl_ii + 488);

    auto ts_xyyyzz_xxyzzz = pbuffer.data(idx_ovl_ii + 489);

    auto ts_xyyyzz_xxzzzz = pbuffer.data(idx_ovl_ii + 490);

    auto ts_xyyyzz_xyyyyy = pbuffer.data(idx_ovl_ii + 491);

    auto ts_xyyyzz_xyyyyz = pbuffer.data(idx_ovl_ii + 492);

    auto ts_xyyyzz_xyyyzz = pbuffer.data(idx_ovl_ii + 493);

    auto ts_xyyyzz_xyyzzz = pbuffer.data(idx_ovl_ii + 494);

    auto ts_xyyyzz_xyzzzz = pbuffer.data(idx_ovl_ii + 495);

    auto ts_xyyyzz_xzzzzz = pbuffer.data(idx_ovl_ii + 496);

    auto ts_xyyyzz_yyyyyy = pbuffer.data(idx_ovl_ii + 497);

    auto ts_xyyyzz_yyyyyz = pbuffer.data(idx_ovl_ii + 498);

    auto ts_xyyyzz_yyyyzz = pbuffer.data(idx_ovl_ii + 499);

    auto ts_xyyyzz_yyyzzz = pbuffer.data(idx_ovl_ii + 500);

    auto ts_xyyyzz_yyzzzz = pbuffer.data(idx_ovl_ii + 501);

    auto ts_xyyyzz_yzzzzz = pbuffer.data(idx_ovl_ii + 502);

    auto ts_xyyyzz_zzzzzz = pbuffer.data(idx_ovl_ii + 503);

    auto ts_xyyzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 504);

    auto ts_xyyzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 505);

    auto ts_xyyzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 506);

    auto ts_xyyzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 507);

    auto ts_xyyzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 508);

    auto ts_xyyzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 509);

    auto ts_xyyzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 510);

    auto ts_xyyzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 511);

    auto ts_xyyzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 512);

    auto ts_xyyzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 513);

    auto ts_xyyzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 514);

    auto ts_xyyzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 515);

    auto ts_xyyzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 516);

    auto ts_xyyzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 517);

    auto ts_xyyzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 518);

    auto ts_xyyzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 519);

    auto ts_xyyzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 520);

    auto ts_xyyzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 521);

    auto ts_xyyzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 522);

    auto ts_xyyzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 523);

    auto ts_xyyzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 524);

    auto ts_xyyzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 525);

    auto ts_xyyzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 526);

    auto ts_xyyzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 527);

    auto ts_xyyzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 528);

    auto ts_xyyzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 529);

    auto ts_xyyzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 530);

    auto ts_xyyzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 531);

    auto ts_xyzzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 532);

    auto ts_xyzzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 533);

    auto ts_xyzzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 534);

    auto ts_xyzzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 535);

    auto ts_xyzzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 536);

    auto ts_xyzzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 537);

    auto ts_xyzzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 538);

    auto ts_xyzzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 539);

    auto ts_xyzzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 540);

    auto ts_xyzzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 541);

    auto ts_xyzzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 542);

    auto ts_xyzzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 543);

    auto ts_xyzzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 544);

    auto ts_xyzzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 545);

    auto ts_xyzzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 546);

    auto ts_xyzzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 547);

    auto ts_xyzzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 548);

    auto ts_xyzzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 549);

    auto ts_xyzzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 550);

    auto ts_xyzzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 551);

    auto ts_xyzzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 552);

    auto ts_xyzzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 553);

    auto ts_xyzzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 554);

    auto ts_xyzzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 555);

    auto ts_xyzzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 556);

    auto ts_xyzzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 557);

    auto ts_xyzzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 558);

    auto ts_xyzzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 559);

    auto ts_xzzzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 560);

    auto ts_xzzzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 561);

    auto ts_xzzzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 562);

    auto ts_xzzzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 563);

    auto ts_xzzzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 564);

    auto ts_xzzzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 565);

    auto ts_xzzzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 566);

    auto ts_xzzzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 567);

    auto ts_xzzzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 568);

    auto ts_xzzzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 569);

    auto ts_xzzzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 570);

    auto ts_xzzzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 571);

    auto ts_xzzzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 572);

    auto ts_xzzzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 573);

    auto ts_xzzzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 574);

    auto ts_xzzzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 575);

    auto ts_xzzzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 576);

    auto ts_xzzzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 577);

    auto ts_xzzzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 578);

    auto ts_xzzzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 579);

    auto ts_xzzzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 580);

    auto ts_xzzzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 581);

    auto ts_xzzzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 582);

    auto ts_xzzzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 583);

    auto ts_xzzzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 584);

    auto ts_xzzzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 585);

    auto ts_xzzzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 586);

    auto ts_xzzzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 587);

    auto ts_yyyyyy_xxxxxx = pbuffer.data(idx_ovl_ii + 588);

    auto ts_yyyyyy_xxxxxy = pbuffer.data(idx_ovl_ii + 589);

    auto ts_yyyyyy_xxxxxz = pbuffer.data(idx_ovl_ii + 590);

    auto ts_yyyyyy_xxxxyy = pbuffer.data(idx_ovl_ii + 591);

    auto ts_yyyyyy_xxxxyz = pbuffer.data(idx_ovl_ii + 592);

    auto ts_yyyyyy_xxxxzz = pbuffer.data(idx_ovl_ii + 593);

    auto ts_yyyyyy_xxxyyy = pbuffer.data(idx_ovl_ii + 594);

    auto ts_yyyyyy_xxxyyz = pbuffer.data(idx_ovl_ii + 595);

    auto ts_yyyyyy_xxxyzz = pbuffer.data(idx_ovl_ii + 596);

    auto ts_yyyyyy_xxxzzz = pbuffer.data(idx_ovl_ii + 597);

    auto ts_yyyyyy_xxyyyy = pbuffer.data(idx_ovl_ii + 598);

    auto ts_yyyyyy_xxyyyz = pbuffer.data(idx_ovl_ii + 599);

    auto ts_yyyyyy_xxyyzz = pbuffer.data(idx_ovl_ii + 600);

    auto ts_yyyyyy_xxyzzz = pbuffer.data(idx_ovl_ii + 601);

    auto ts_yyyyyy_xxzzzz = pbuffer.data(idx_ovl_ii + 602);

    auto ts_yyyyyy_xyyyyy = pbuffer.data(idx_ovl_ii + 603);

    auto ts_yyyyyy_xyyyyz = pbuffer.data(idx_ovl_ii + 604);

    auto ts_yyyyyy_xyyyzz = pbuffer.data(idx_ovl_ii + 605);

    auto ts_yyyyyy_xyyzzz = pbuffer.data(idx_ovl_ii + 606);

    auto ts_yyyyyy_xyzzzz = pbuffer.data(idx_ovl_ii + 607);

    auto ts_yyyyyy_xzzzzz = pbuffer.data(idx_ovl_ii + 608);

    auto ts_yyyyyy_yyyyyy = pbuffer.data(idx_ovl_ii + 609);

    auto ts_yyyyyy_yyyyyz = pbuffer.data(idx_ovl_ii + 610);

    auto ts_yyyyyy_yyyyzz = pbuffer.data(idx_ovl_ii + 611);

    auto ts_yyyyyy_yyyzzz = pbuffer.data(idx_ovl_ii + 612);

    auto ts_yyyyyy_yyzzzz = pbuffer.data(idx_ovl_ii + 613);

    auto ts_yyyyyy_yzzzzz = pbuffer.data(idx_ovl_ii + 614);

    auto ts_yyyyyy_zzzzzz = pbuffer.data(idx_ovl_ii + 615);

    auto ts_yyyyyz_xxxxxx = pbuffer.data(idx_ovl_ii + 616);

    auto ts_yyyyyz_xxxxxy = pbuffer.data(idx_ovl_ii + 617);

    auto ts_yyyyyz_xxxxxz = pbuffer.data(idx_ovl_ii + 618);

    auto ts_yyyyyz_xxxxyy = pbuffer.data(idx_ovl_ii + 619);

    auto ts_yyyyyz_xxxxyz = pbuffer.data(idx_ovl_ii + 620);

    auto ts_yyyyyz_xxxxzz = pbuffer.data(idx_ovl_ii + 621);

    auto ts_yyyyyz_xxxyyy = pbuffer.data(idx_ovl_ii + 622);

    auto ts_yyyyyz_xxxyyz = pbuffer.data(idx_ovl_ii + 623);

    auto ts_yyyyyz_xxxyzz = pbuffer.data(idx_ovl_ii + 624);

    auto ts_yyyyyz_xxxzzz = pbuffer.data(idx_ovl_ii + 625);

    auto ts_yyyyyz_xxyyyy = pbuffer.data(idx_ovl_ii + 626);

    auto ts_yyyyyz_xxyyyz = pbuffer.data(idx_ovl_ii + 627);

    auto ts_yyyyyz_xxyyzz = pbuffer.data(idx_ovl_ii + 628);

    auto ts_yyyyyz_xxyzzz = pbuffer.data(idx_ovl_ii + 629);

    auto ts_yyyyyz_xxzzzz = pbuffer.data(idx_ovl_ii + 630);

    auto ts_yyyyyz_xyyyyy = pbuffer.data(idx_ovl_ii + 631);

    auto ts_yyyyyz_xyyyyz = pbuffer.data(idx_ovl_ii + 632);

    auto ts_yyyyyz_xyyyzz = pbuffer.data(idx_ovl_ii + 633);

    auto ts_yyyyyz_xyyzzz = pbuffer.data(idx_ovl_ii + 634);

    auto ts_yyyyyz_xyzzzz = pbuffer.data(idx_ovl_ii + 635);

    auto ts_yyyyyz_xzzzzz = pbuffer.data(idx_ovl_ii + 636);

    auto ts_yyyyyz_yyyyyy = pbuffer.data(idx_ovl_ii + 637);

    auto ts_yyyyyz_yyyyyz = pbuffer.data(idx_ovl_ii + 638);

    auto ts_yyyyyz_yyyyzz = pbuffer.data(idx_ovl_ii + 639);

    auto ts_yyyyyz_yyyzzz = pbuffer.data(idx_ovl_ii + 640);

    auto ts_yyyyyz_yyzzzz = pbuffer.data(idx_ovl_ii + 641);

    auto ts_yyyyyz_yzzzzz = pbuffer.data(idx_ovl_ii + 642);

    auto ts_yyyyyz_zzzzzz = pbuffer.data(idx_ovl_ii + 643);

    auto ts_yyyyzz_xxxxxx = pbuffer.data(idx_ovl_ii + 644);

    auto ts_yyyyzz_xxxxxy = pbuffer.data(idx_ovl_ii + 645);

    auto ts_yyyyzz_xxxxxz = pbuffer.data(idx_ovl_ii + 646);

    auto ts_yyyyzz_xxxxyy = pbuffer.data(idx_ovl_ii + 647);

    auto ts_yyyyzz_xxxxyz = pbuffer.data(idx_ovl_ii + 648);

    auto ts_yyyyzz_xxxxzz = pbuffer.data(idx_ovl_ii + 649);

    auto ts_yyyyzz_xxxyyy = pbuffer.data(idx_ovl_ii + 650);

    auto ts_yyyyzz_xxxyyz = pbuffer.data(idx_ovl_ii + 651);

    auto ts_yyyyzz_xxxyzz = pbuffer.data(idx_ovl_ii + 652);

    auto ts_yyyyzz_xxxzzz = pbuffer.data(idx_ovl_ii + 653);

    auto ts_yyyyzz_xxyyyy = pbuffer.data(idx_ovl_ii + 654);

    auto ts_yyyyzz_xxyyyz = pbuffer.data(idx_ovl_ii + 655);

    auto ts_yyyyzz_xxyyzz = pbuffer.data(idx_ovl_ii + 656);

    auto ts_yyyyzz_xxyzzz = pbuffer.data(idx_ovl_ii + 657);

    auto ts_yyyyzz_xxzzzz = pbuffer.data(idx_ovl_ii + 658);

    auto ts_yyyyzz_xyyyyy = pbuffer.data(idx_ovl_ii + 659);

    auto ts_yyyyzz_xyyyyz = pbuffer.data(idx_ovl_ii + 660);

    auto ts_yyyyzz_xyyyzz = pbuffer.data(idx_ovl_ii + 661);

    auto ts_yyyyzz_xyyzzz = pbuffer.data(idx_ovl_ii + 662);

    auto ts_yyyyzz_xyzzzz = pbuffer.data(idx_ovl_ii + 663);

    auto ts_yyyyzz_xzzzzz = pbuffer.data(idx_ovl_ii + 664);

    auto ts_yyyyzz_yyyyyy = pbuffer.data(idx_ovl_ii + 665);

    auto ts_yyyyzz_yyyyyz = pbuffer.data(idx_ovl_ii + 666);

    auto ts_yyyyzz_yyyyzz = pbuffer.data(idx_ovl_ii + 667);

    auto ts_yyyyzz_yyyzzz = pbuffer.data(idx_ovl_ii + 668);

    auto ts_yyyyzz_yyzzzz = pbuffer.data(idx_ovl_ii + 669);

    auto ts_yyyyzz_yzzzzz = pbuffer.data(idx_ovl_ii + 670);

    auto ts_yyyyzz_zzzzzz = pbuffer.data(idx_ovl_ii + 671);

    auto ts_yyyzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 672);

    auto ts_yyyzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 673);

    auto ts_yyyzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 674);

    auto ts_yyyzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 675);

    auto ts_yyyzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 676);

    auto ts_yyyzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 677);

    auto ts_yyyzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 678);

    auto ts_yyyzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 679);

    auto ts_yyyzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 680);

    auto ts_yyyzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 681);

    auto ts_yyyzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 682);

    auto ts_yyyzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 683);

    auto ts_yyyzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 684);

    auto ts_yyyzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 685);

    auto ts_yyyzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 686);

    auto ts_yyyzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 687);

    auto ts_yyyzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 688);

    auto ts_yyyzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 689);

    auto ts_yyyzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 690);

    auto ts_yyyzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 691);

    auto ts_yyyzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 692);

    auto ts_yyyzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 693);

    auto ts_yyyzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 694);

    auto ts_yyyzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 695);

    auto ts_yyyzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 696);

    auto ts_yyyzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 697);

    auto ts_yyyzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 698);

    auto ts_yyyzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 699);

    auto ts_yyzzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 700);

    auto ts_yyzzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 701);

    auto ts_yyzzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 702);

    auto ts_yyzzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 703);

    auto ts_yyzzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 704);

    auto ts_yyzzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 705);

    auto ts_yyzzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 706);

    auto ts_yyzzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 707);

    auto ts_yyzzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 708);

    auto ts_yyzzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 709);

    auto ts_yyzzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 710);

    auto ts_yyzzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 711);

    auto ts_yyzzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 712);

    auto ts_yyzzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 713);

    auto ts_yyzzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 714);

    auto ts_yyzzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 715);

    auto ts_yyzzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 716);

    auto ts_yyzzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 717);

    auto ts_yyzzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 718);

    auto ts_yyzzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 719);

    auto ts_yyzzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 720);

    auto ts_yyzzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 721);

    auto ts_yyzzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 722);

    auto ts_yyzzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 723);

    auto ts_yyzzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 724);

    auto ts_yyzzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 725);

    auto ts_yyzzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 726);

    auto ts_yyzzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 727);

    auto ts_yzzzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 728);

    auto ts_yzzzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 729);

    auto ts_yzzzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 730);

    auto ts_yzzzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 731);

    auto ts_yzzzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 732);

    auto ts_yzzzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 733);

    auto ts_yzzzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 734);

    auto ts_yzzzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 735);

    auto ts_yzzzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 736);

    auto ts_yzzzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 737);

    auto ts_yzzzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 738);

    auto ts_yzzzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 739);

    auto ts_yzzzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 740);

    auto ts_yzzzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 741);

    auto ts_yzzzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 742);

    auto ts_yzzzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 743);

    auto ts_yzzzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 744);

    auto ts_yzzzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 745);

    auto ts_yzzzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 746);

    auto ts_yzzzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 747);

    auto ts_yzzzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 748);

    auto ts_yzzzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 749);

    auto ts_yzzzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 750);

    auto ts_yzzzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 751);

    auto ts_yzzzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 752);

    auto ts_yzzzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 753);

    auto ts_yzzzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 754);

    auto ts_yzzzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 755);

    auto ts_zzzzzz_xxxxxx = pbuffer.data(idx_ovl_ii + 756);

    auto ts_zzzzzz_xxxxxy = pbuffer.data(idx_ovl_ii + 757);

    auto ts_zzzzzz_xxxxxz = pbuffer.data(idx_ovl_ii + 758);

    auto ts_zzzzzz_xxxxyy = pbuffer.data(idx_ovl_ii + 759);

    auto ts_zzzzzz_xxxxyz = pbuffer.data(idx_ovl_ii + 760);

    auto ts_zzzzzz_xxxxzz = pbuffer.data(idx_ovl_ii + 761);

    auto ts_zzzzzz_xxxyyy = pbuffer.data(idx_ovl_ii + 762);

    auto ts_zzzzzz_xxxyyz = pbuffer.data(idx_ovl_ii + 763);

    auto ts_zzzzzz_xxxyzz = pbuffer.data(idx_ovl_ii + 764);

    auto ts_zzzzzz_xxxzzz = pbuffer.data(idx_ovl_ii + 765);

    auto ts_zzzzzz_xxyyyy = pbuffer.data(idx_ovl_ii + 766);

    auto ts_zzzzzz_xxyyyz = pbuffer.data(idx_ovl_ii + 767);

    auto ts_zzzzzz_xxyyzz = pbuffer.data(idx_ovl_ii + 768);

    auto ts_zzzzzz_xxyzzz = pbuffer.data(idx_ovl_ii + 769);

    auto ts_zzzzzz_xxzzzz = pbuffer.data(idx_ovl_ii + 770);

    auto ts_zzzzzz_xyyyyy = pbuffer.data(idx_ovl_ii + 771);

    auto ts_zzzzzz_xyyyyz = pbuffer.data(idx_ovl_ii + 772);

    auto ts_zzzzzz_xyyyzz = pbuffer.data(idx_ovl_ii + 773);

    auto ts_zzzzzz_xyyzzz = pbuffer.data(idx_ovl_ii + 774);

    auto ts_zzzzzz_xyzzzz = pbuffer.data(idx_ovl_ii + 775);

    auto ts_zzzzzz_xzzzzz = pbuffer.data(idx_ovl_ii + 776);

    auto ts_zzzzzz_yyyyyy = pbuffer.data(idx_ovl_ii + 777);

    auto ts_zzzzzz_yyyyyz = pbuffer.data(idx_ovl_ii + 778);

    auto ts_zzzzzz_yyyyzz = pbuffer.data(idx_ovl_ii + 779);

    auto ts_zzzzzz_yyyzzz = pbuffer.data(idx_ovl_ii + 780);

    auto ts_zzzzzz_yyzzzz = pbuffer.data(idx_ovl_ii + 781);

    auto ts_zzzzzz_yzzzzz = pbuffer.data(idx_ovl_ii + 782);

    auto ts_zzzzzz_zzzzzz = pbuffer.data(idx_ovl_ii + 783);

    // Set up 0-28 components of targeted buffer : II

    auto tk_xxxxxx_xxxxxx = pbuffer.data(idx_kin_ii);

    auto tk_xxxxxx_xxxxxy = pbuffer.data(idx_kin_ii + 1);

    auto tk_xxxxxx_xxxxxz = pbuffer.data(idx_kin_ii + 2);

    auto tk_xxxxxx_xxxxyy = pbuffer.data(idx_kin_ii + 3);

    auto tk_xxxxxx_xxxxyz = pbuffer.data(idx_kin_ii + 4);

    auto tk_xxxxxx_xxxxzz = pbuffer.data(idx_kin_ii + 5);

    auto tk_xxxxxx_xxxyyy = pbuffer.data(idx_kin_ii + 6);

    auto tk_xxxxxx_xxxyyz = pbuffer.data(idx_kin_ii + 7);

    auto tk_xxxxxx_xxxyzz = pbuffer.data(idx_kin_ii + 8);

    auto tk_xxxxxx_xxxzzz = pbuffer.data(idx_kin_ii + 9);

    auto tk_xxxxxx_xxyyyy = pbuffer.data(idx_kin_ii + 10);

    auto tk_xxxxxx_xxyyyz = pbuffer.data(idx_kin_ii + 11);

    auto tk_xxxxxx_xxyyzz = pbuffer.data(idx_kin_ii + 12);

    auto tk_xxxxxx_xxyzzz = pbuffer.data(idx_kin_ii + 13);

    auto tk_xxxxxx_xxzzzz = pbuffer.data(idx_kin_ii + 14);

    auto tk_xxxxxx_xyyyyy = pbuffer.data(idx_kin_ii + 15);

    auto tk_xxxxxx_xyyyyz = pbuffer.data(idx_kin_ii + 16);

    auto tk_xxxxxx_xyyyzz = pbuffer.data(idx_kin_ii + 17);

    auto tk_xxxxxx_xyyzzz = pbuffer.data(idx_kin_ii + 18);

    auto tk_xxxxxx_xyzzzz = pbuffer.data(idx_kin_ii + 19);

    auto tk_xxxxxx_xzzzzz = pbuffer.data(idx_kin_ii + 20);

    auto tk_xxxxxx_yyyyyy = pbuffer.data(idx_kin_ii + 21);

    auto tk_xxxxxx_yyyyyz = pbuffer.data(idx_kin_ii + 22);

    auto tk_xxxxxx_yyyyzz = pbuffer.data(idx_kin_ii + 23);

    auto tk_xxxxxx_yyyzzz = pbuffer.data(idx_kin_ii + 24);

    auto tk_xxxxxx_yyzzzz = pbuffer.data(idx_kin_ii + 25);

    auto tk_xxxxxx_yzzzzz = pbuffer.data(idx_kin_ii + 26);

    auto tk_xxxxxx_zzzzzz = pbuffer.data(idx_kin_ii + 27);

#pragma omp simd aligned(pa_x,                 \
                             tk_xxxx_xxxxxx,   \
                             tk_xxxx_xxxxxy,   \
                             tk_xxxx_xxxxxz,   \
                             tk_xxxx_xxxxyy,   \
                             tk_xxxx_xxxxyz,   \
                             tk_xxxx_xxxxzz,   \
                             tk_xxxx_xxxyyy,   \
                             tk_xxxx_xxxyyz,   \
                             tk_xxxx_xxxyzz,   \
                             tk_xxxx_xxxzzz,   \
                             tk_xxxx_xxyyyy,   \
                             tk_xxxx_xxyyyz,   \
                             tk_xxxx_xxyyzz,   \
                             tk_xxxx_xxyzzz,   \
                             tk_xxxx_xxzzzz,   \
                             tk_xxxx_xyyyyy,   \
                             tk_xxxx_xyyyyz,   \
                             tk_xxxx_xyyyzz,   \
                             tk_xxxx_xyyzzz,   \
                             tk_xxxx_xyzzzz,   \
                             tk_xxxx_xzzzzz,   \
                             tk_xxxx_yyyyyy,   \
                             tk_xxxx_yyyyyz,   \
                             tk_xxxx_yyyyzz,   \
                             tk_xxxx_yyyzzz,   \
                             tk_xxxx_yyzzzz,   \
                             tk_xxxx_yzzzzz,   \
                             tk_xxxx_zzzzzz,   \
                             tk_xxxxx_xxxxx,   \
                             tk_xxxxx_xxxxxx,  \
                             tk_xxxxx_xxxxxy,  \
                             tk_xxxxx_xxxxxz,  \
                             tk_xxxxx_xxxxy,   \
                             tk_xxxxx_xxxxyy,  \
                             tk_xxxxx_xxxxyz,  \
                             tk_xxxxx_xxxxz,   \
                             tk_xxxxx_xxxxzz,  \
                             tk_xxxxx_xxxyy,   \
                             tk_xxxxx_xxxyyy,  \
                             tk_xxxxx_xxxyyz,  \
                             tk_xxxxx_xxxyz,   \
                             tk_xxxxx_xxxyzz,  \
                             tk_xxxxx_xxxzz,   \
                             tk_xxxxx_xxxzzz,  \
                             tk_xxxxx_xxyyy,   \
                             tk_xxxxx_xxyyyy,  \
                             tk_xxxxx_xxyyyz,  \
                             tk_xxxxx_xxyyz,   \
                             tk_xxxxx_xxyyzz,  \
                             tk_xxxxx_xxyzz,   \
                             tk_xxxxx_xxyzzz,  \
                             tk_xxxxx_xxzzz,   \
                             tk_xxxxx_xxzzzz,  \
                             tk_xxxxx_xyyyy,   \
                             tk_xxxxx_xyyyyy,  \
                             tk_xxxxx_xyyyyz,  \
                             tk_xxxxx_xyyyz,   \
                             tk_xxxxx_xyyyzz,  \
                             tk_xxxxx_xyyzz,   \
                             tk_xxxxx_xyyzzz,  \
                             tk_xxxxx_xyzzz,   \
                             tk_xxxxx_xyzzzz,  \
                             tk_xxxxx_xzzzz,   \
                             tk_xxxxx_xzzzzz,  \
                             tk_xxxxx_yyyyy,   \
                             tk_xxxxx_yyyyyy,  \
                             tk_xxxxx_yyyyyz,  \
                             tk_xxxxx_yyyyz,   \
                             tk_xxxxx_yyyyzz,  \
                             tk_xxxxx_yyyzz,   \
                             tk_xxxxx_yyyzzz,  \
                             tk_xxxxx_yyzzz,   \
                             tk_xxxxx_yyzzzz,  \
                             tk_xxxxx_yzzzz,   \
                             tk_xxxxx_yzzzzz,  \
                             tk_xxxxx_zzzzz,   \
                             tk_xxxxx_zzzzzz,  \
                             tk_xxxxxx_xxxxxx, \
                             tk_xxxxxx_xxxxxy, \
                             tk_xxxxxx_xxxxxz, \
                             tk_xxxxxx_xxxxyy, \
                             tk_xxxxxx_xxxxyz, \
                             tk_xxxxxx_xxxxzz, \
                             tk_xxxxxx_xxxyyy, \
                             tk_xxxxxx_xxxyyz, \
                             tk_xxxxxx_xxxyzz, \
                             tk_xxxxxx_xxxzzz, \
                             tk_xxxxxx_xxyyyy, \
                             tk_xxxxxx_xxyyyz, \
                             tk_xxxxxx_xxyyzz, \
                             tk_xxxxxx_xxyzzz, \
                             tk_xxxxxx_xxzzzz, \
                             tk_xxxxxx_xyyyyy, \
                             tk_xxxxxx_xyyyyz, \
                             tk_xxxxxx_xyyyzz, \
                             tk_xxxxxx_xyyzzz, \
                             tk_xxxxxx_xyzzzz, \
                             tk_xxxxxx_xzzzzz, \
                             tk_xxxxxx_yyyyyy, \
                             tk_xxxxxx_yyyyyz, \
                             tk_xxxxxx_yyyyzz, \
                             tk_xxxxxx_yyyzzz, \
                             tk_xxxxxx_yyzzzz, \
                             tk_xxxxxx_yzzzzz, \
                             tk_xxxxxx_zzzzzz, \
                             ts_xxxx_xxxxxx,   \
                             ts_xxxx_xxxxxy,   \
                             ts_xxxx_xxxxxz,   \
                             ts_xxxx_xxxxyy,   \
                             ts_xxxx_xxxxyz,   \
                             ts_xxxx_xxxxzz,   \
                             ts_xxxx_xxxyyy,   \
                             ts_xxxx_xxxyyz,   \
                             ts_xxxx_xxxyzz,   \
                             ts_xxxx_xxxzzz,   \
                             ts_xxxx_xxyyyy,   \
                             ts_xxxx_xxyyyz,   \
                             ts_xxxx_xxyyzz,   \
                             ts_xxxx_xxyzzz,   \
                             ts_xxxx_xxzzzz,   \
                             ts_xxxx_xyyyyy,   \
                             ts_xxxx_xyyyyz,   \
                             ts_xxxx_xyyyzz,   \
                             ts_xxxx_xyyzzz,   \
                             ts_xxxx_xyzzzz,   \
                             ts_xxxx_xzzzzz,   \
                             ts_xxxx_yyyyyy,   \
                             ts_xxxx_yyyyyz,   \
                             ts_xxxx_yyyyzz,   \
                             ts_xxxx_yyyzzz,   \
                             ts_xxxx_yyzzzz,   \
                             ts_xxxx_yzzzzz,   \
                             ts_xxxx_zzzzzz,   \
                             ts_xxxxxx_xxxxxx, \
                             ts_xxxxxx_xxxxxy, \
                             ts_xxxxxx_xxxxxz, \
                             ts_xxxxxx_xxxxyy, \
                             ts_xxxxxx_xxxxyz, \
                             ts_xxxxxx_xxxxzz, \
                             ts_xxxxxx_xxxyyy, \
                             ts_xxxxxx_xxxyyz, \
                             ts_xxxxxx_xxxyzz, \
                             ts_xxxxxx_xxxzzz, \
                             ts_xxxxxx_xxyyyy, \
                             ts_xxxxxx_xxyyyz, \
                             ts_xxxxxx_xxyyzz, \
                             ts_xxxxxx_xxyzzz, \
                             ts_xxxxxx_xxzzzz, \
                             ts_xxxxxx_xyyyyy, \
                             ts_xxxxxx_xyyyyz, \
                             ts_xxxxxx_xyyyzz, \
                             ts_xxxxxx_xyyzzz, \
                             ts_xxxxxx_xyzzzz, \
                             ts_xxxxxx_xzzzzz, \
                             ts_xxxxxx_yyyyyy, \
                             ts_xxxxxx_yyyyyz, \
                             ts_xxxxxx_yyyyzz, \
                             ts_xxxxxx_yyyzzz, \
                             ts_xxxxxx_yyzzzz, \
                             ts_xxxxxx_yzzzzz, \
                             ts_xxxxxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_xxxxxx[i] = -10.0 * ts_xxxx_xxxxxx[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxxx[i] * fe_0 + 6.0 * tk_xxxxx_xxxxx[i] * fe_0 +
                              tk_xxxxx_xxxxxx[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxxx[i] * fz_0;

        tk_xxxxxx_xxxxxy[i] = -10.0 * ts_xxxx_xxxxxy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxxy[i] * fe_0 + 5.0 * tk_xxxxx_xxxxy[i] * fe_0 +
                              tk_xxxxx_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxxy[i] * fz_0;

        tk_xxxxxx_xxxxxz[i] = -10.0 * ts_xxxx_xxxxxz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxxz[i] * fe_0 + 5.0 * tk_xxxxx_xxxxz[i] * fe_0 +
                              tk_xxxxx_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxxz[i] * fz_0;

        tk_xxxxxx_xxxxyy[i] = -10.0 * ts_xxxx_xxxxyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxyy[i] * fe_0 + 4.0 * tk_xxxxx_xxxyy[i] * fe_0 +
                              tk_xxxxx_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxyy[i] * fz_0;

        tk_xxxxxx_xxxxyz[i] = -10.0 * ts_xxxx_xxxxyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxyz[i] * fe_0 + 4.0 * tk_xxxxx_xxxyz[i] * fe_0 +
                              tk_xxxxx_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxyz[i] * fz_0;

        tk_xxxxxx_xxxxzz[i] = -10.0 * ts_xxxx_xxxxzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxzz[i] * fe_0 + 4.0 * tk_xxxxx_xxxzz[i] * fe_0 +
                              tk_xxxxx_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxzz[i] * fz_0;

        tk_xxxxxx_xxxyyy[i] = -10.0 * ts_xxxx_xxxyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxyyy[i] * fe_0 + 3.0 * tk_xxxxx_xxyyy[i] * fe_0 +
                              tk_xxxxx_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxyyy[i] * fz_0;

        tk_xxxxxx_xxxyyz[i] = -10.0 * ts_xxxx_xxxyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxyyz[i] * fe_0 + 3.0 * tk_xxxxx_xxyyz[i] * fe_0 +
                              tk_xxxxx_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxyyz[i] * fz_0;

        tk_xxxxxx_xxxyzz[i] = -10.0 * ts_xxxx_xxxyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxyzz[i] * fe_0 + 3.0 * tk_xxxxx_xxyzz[i] * fe_0 +
                              tk_xxxxx_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxyzz[i] * fz_0;

        tk_xxxxxx_xxxzzz[i] = -10.0 * ts_xxxx_xxxzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxzzz[i] * fe_0 + 3.0 * tk_xxxxx_xxzzz[i] * fe_0 +
                              tk_xxxxx_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxzzz[i] * fz_0;

        tk_xxxxxx_xxyyyy[i] = -10.0 * ts_xxxx_xxyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyyyy[i] * fe_0 + 2.0 * tk_xxxxx_xyyyy[i] * fe_0 +
                              tk_xxxxx_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyyyy[i] * fz_0;

        tk_xxxxxx_xxyyyz[i] = -10.0 * ts_xxxx_xxyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyyyz[i] * fe_0 + 2.0 * tk_xxxxx_xyyyz[i] * fe_0 +
                              tk_xxxxx_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyyyz[i] * fz_0;

        tk_xxxxxx_xxyyzz[i] = -10.0 * ts_xxxx_xxyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyyzz[i] * fe_0 + 2.0 * tk_xxxxx_xyyzz[i] * fe_0 +
                              tk_xxxxx_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyyzz[i] * fz_0;

        tk_xxxxxx_xxyzzz[i] = -10.0 * ts_xxxx_xxyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyzzz[i] * fe_0 + 2.0 * tk_xxxxx_xyzzz[i] * fe_0 +
                              tk_xxxxx_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyzzz[i] * fz_0;

        tk_xxxxxx_xxzzzz[i] = -10.0 * ts_xxxx_xxzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxzzzz[i] * fe_0 + 2.0 * tk_xxxxx_xzzzz[i] * fe_0 +
                              tk_xxxxx_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxzzzz[i] * fz_0;

        tk_xxxxxx_xyyyyy[i] = -10.0 * ts_xxxx_xyyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyyyy[i] * fe_0 + tk_xxxxx_yyyyy[i] * fe_0 +
                              tk_xxxxx_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyyyy[i] * fz_0;

        tk_xxxxxx_xyyyyz[i] = -10.0 * ts_xxxx_xyyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyyyz[i] * fe_0 + tk_xxxxx_yyyyz[i] * fe_0 +
                              tk_xxxxx_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyyyz[i] * fz_0;

        tk_xxxxxx_xyyyzz[i] = -10.0 * ts_xxxx_xyyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyyzz[i] * fe_0 + tk_xxxxx_yyyzz[i] * fe_0 +
                              tk_xxxxx_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyyzz[i] * fz_0;

        tk_xxxxxx_xyyzzz[i] = -10.0 * ts_xxxx_xyyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyzzz[i] * fe_0 + tk_xxxxx_yyzzz[i] * fe_0 +
                              tk_xxxxx_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyzzz[i] * fz_0;

        tk_xxxxxx_xyzzzz[i] = -10.0 * ts_xxxx_xyzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyzzzz[i] * fe_0 + tk_xxxxx_yzzzz[i] * fe_0 +
                              tk_xxxxx_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyzzzz[i] * fz_0;

        tk_xxxxxx_xzzzzz[i] = -10.0 * ts_xxxx_xzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xzzzzz[i] * fe_0 + tk_xxxxx_zzzzz[i] * fe_0 +
                              tk_xxxxx_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xzzzzz[i] * fz_0;

        tk_xxxxxx_yyyyyy[i] = -10.0 * ts_xxxx_yyyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyyyy[i] * fe_0 + tk_xxxxx_yyyyyy[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_yyyyyy[i] * fz_0;

        tk_xxxxxx_yyyyyz[i] = -10.0 * ts_xxxx_yyyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyyyz[i] * fe_0 + tk_xxxxx_yyyyyz[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_yyyyyz[i] * fz_0;

        tk_xxxxxx_yyyyzz[i] = -10.0 * ts_xxxx_yyyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyyzz[i] * fe_0 + tk_xxxxx_yyyyzz[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_yyyyzz[i] * fz_0;

        tk_xxxxxx_yyyzzz[i] = -10.0 * ts_xxxx_yyyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyzzz[i] * fe_0 + tk_xxxxx_yyyzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_yyyzzz[i] * fz_0;

        tk_xxxxxx_yyzzzz[i] = -10.0 * ts_xxxx_yyzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyzzzz[i] * fe_0 + tk_xxxxx_yyzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_yyzzzz[i] * fz_0;

        tk_xxxxxx_yzzzzz[i] = -10.0 * ts_xxxx_yzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yzzzzz[i] * fe_0 + tk_xxxxx_yzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_yzzzzz[i] * fz_0;

        tk_xxxxxx_zzzzzz[i] = -10.0 * ts_xxxx_zzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_zzzzzz[i] * fe_0 + tk_xxxxx_zzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxxx_zzzzzz[i] * fz_0;
    }

    // Set up 28-56 components of targeted buffer : II

    auto tk_xxxxxy_xxxxxx = pbuffer.data(idx_kin_ii + 28);

    auto tk_xxxxxy_xxxxxy = pbuffer.data(idx_kin_ii + 29);

    auto tk_xxxxxy_xxxxxz = pbuffer.data(idx_kin_ii + 30);

    auto tk_xxxxxy_xxxxyy = pbuffer.data(idx_kin_ii + 31);

    auto tk_xxxxxy_xxxxyz = pbuffer.data(idx_kin_ii + 32);

    auto tk_xxxxxy_xxxxzz = pbuffer.data(idx_kin_ii + 33);

    auto tk_xxxxxy_xxxyyy = pbuffer.data(idx_kin_ii + 34);

    auto tk_xxxxxy_xxxyyz = pbuffer.data(idx_kin_ii + 35);

    auto tk_xxxxxy_xxxyzz = pbuffer.data(idx_kin_ii + 36);

    auto tk_xxxxxy_xxxzzz = pbuffer.data(idx_kin_ii + 37);

    auto tk_xxxxxy_xxyyyy = pbuffer.data(idx_kin_ii + 38);

    auto tk_xxxxxy_xxyyyz = pbuffer.data(idx_kin_ii + 39);

    auto tk_xxxxxy_xxyyzz = pbuffer.data(idx_kin_ii + 40);

    auto tk_xxxxxy_xxyzzz = pbuffer.data(idx_kin_ii + 41);

    auto tk_xxxxxy_xxzzzz = pbuffer.data(idx_kin_ii + 42);

    auto tk_xxxxxy_xyyyyy = pbuffer.data(idx_kin_ii + 43);

    auto tk_xxxxxy_xyyyyz = pbuffer.data(idx_kin_ii + 44);

    auto tk_xxxxxy_xyyyzz = pbuffer.data(idx_kin_ii + 45);

    auto tk_xxxxxy_xyyzzz = pbuffer.data(idx_kin_ii + 46);

    auto tk_xxxxxy_xyzzzz = pbuffer.data(idx_kin_ii + 47);

    auto tk_xxxxxy_xzzzzz = pbuffer.data(idx_kin_ii + 48);

    auto tk_xxxxxy_yyyyyy = pbuffer.data(idx_kin_ii + 49);

    auto tk_xxxxxy_yyyyyz = pbuffer.data(idx_kin_ii + 50);

    auto tk_xxxxxy_yyyyzz = pbuffer.data(idx_kin_ii + 51);

    auto tk_xxxxxy_yyyzzz = pbuffer.data(idx_kin_ii + 52);

    auto tk_xxxxxy_yyzzzz = pbuffer.data(idx_kin_ii + 53);

    auto tk_xxxxxy_yzzzzz = pbuffer.data(idx_kin_ii + 54);

    auto tk_xxxxxy_zzzzzz = pbuffer.data(idx_kin_ii + 55);

#pragma omp simd aligned(pa_y,                 \
                             tk_xxxxx_xxxxx,   \
                             tk_xxxxx_xxxxxx,  \
                             tk_xxxxx_xxxxxy,  \
                             tk_xxxxx_xxxxxz,  \
                             tk_xxxxx_xxxxy,   \
                             tk_xxxxx_xxxxyy,  \
                             tk_xxxxx_xxxxyz,  \
                             tk_xxxxx_xxxxz,   \
                             tk_xxxxx_xxxxzz,  \
                             tk_xxxxx_xxxyy,   \
                             tk_xxxxx_xxxyyy,  \
                             tk_xxxxx_xxxyyz,  \
                             tk_xxxxx_xxxyz,   \
                             tk_xxxxx_xxxyzz,  \
                             tk_xxxxx_xxxzz,   \
                             tk_xxxxx_xxxzzz,  \
                             tk_xxxxx_xxyyy,   \
                             tk_xxxxx_xxyyyy,  \
                             tk_xxxxx_xxyyyz,  \
                             tk_xxxxx_xxyyz,   \
                             tk_xxxxx_xxyyzz,  \
                             tk_xxxxx_xxyzz,   \
                             tk_xxxxx_xxyzzz,  \
                             tk_xxxxx_xxzzz,   \
                             tk_xxxxx_xxzzzz,  \
                             tk_xxxxx_xyyyy,   \
                             tk_xxxxx_xyyyyy,  \
                             tk_xxxxx_xyyyyz,  \
                             tk_xxxxx_xyyyz,   \
                             tk_xxxxx_xyyyzz,  \
                             tk_xxxxx_xyyzz,   \
                             tk_xxxxx_xyyzzz,  \
                             tk_xxxxx_xyzzz,   \
                             tk_xxxxx_xyzzzz,  \
                             tk_xxxxx_xzzzz,   \
                             tk_xxxxx_xzzzzz,  \
                             tk_xxxxx_yyyyy,   \
                             tk_xxxxx_yyyyyy,  \
                             tk_xxxxx_yyyyyz,  \
                             tk_xxxxx_yyyyz,   \
                             tk_xxxxx_yyyyzz,  \
                             tk_xxxxx_yyyzz,   \
                             tk_xxxxx_yyyzzz,  \
                             tk_xxxxx_yyzzz,   \
                             tk_xxxxx_yyzzzz,  \
                             tk_xxxxx_yzzzz,   \
                             tk_xxxxx_yzzzzz,  \
                             tk_xxxxx_zzzzz,   \
                             tk_xxxxx_zzzzzz,  \
                             tk_xxxxxy_xxxxxx, \
                             tk_xxxxxy_xxxxxy, \
                             tk_xxxxxy_xxxxxz, \
                             tk_xxxxxy_xxxxyy, \
                             tk_xxxxxy_xxxxyz, \
                             tk_xxxxxy_xxxxzz, \
                             tk_xxxxxy_xxxyyy, \
                             tk_xxxxxy_xxxyyz, \
                             tk_xxxxxy_xxxyzz, \
                             tk_xxxxxy_xxxzzz, \
                             tk_xxxxxy_xxyyyy, \
                             tk_xxxxxy_xxyyyz, \
                             tk_xxxxxy_xxyyzz, \
                             tk_xxxxxy_xxyzzz, \
                             tk_xxxxxy_xxzzzz, \
                             tk_xxxxxy_xyyyyy, \
                             tk_xxxxxy_xyyyyz, \
                             tk_xxxxxy_xyyyzz, \
                             tk_xxxxxy_xyyzzz, \
                             tk_xxxxxy_xyzzzz, \
                             tk_xxxxxy_xzzzzz, \
                             tk_xxxxxy_yyyyyy, \
                             tk_xxxxxy_yyyyyz, \
                             tk_xxxxxy_yyyyzz, \
                             tk_xxxxxy_yyyzzz, \
                             tk_xxxxxy_yyzzzz, \
                             tk_xxxxxy_yzzzzz, \
                             tk_xxxxxy_zzzzzz, \
                             ts_xxxxxy_xxxxxx, \
                             ts_xxxxxy_xxxxxy, \
                             ts_xxxxxy_xxxxxz, \
                             ts_xxxxxy_xxxxyy, \
                             ts_xxxxxy_xxxxyz, \
                             ts_xxxxxy_xxxxzz, \
                             ts_xxxxxy_xxxyyy, \
                             ts_xxxxxy_xxxyyz, \
                             ts_xxxxxy_xxxyzz, \
                             ts_xxxxxy_xxxzzz, \
                             ts_xxxxxy_xxyyyy, \
                             ts_xxxxxy_xxyyyz, \
                             ts_xxxxxy_xxyyzz, \
                             ts_xxxxxy_xxyzzz, \
                             ts_xxxxxy_xxzzzz, \
                             ts_xxxxxy_xyyyyy, \
                             ts_xxxxxy_xyyyyz, \
                             ts_xxxxxy_xyyyzz, \
                             ts_xxxxxy_xyyzzz, \
                             ts_xxxxxy_xyzzzz, \
                             ts_xxxxxy_xzzzzz, \
                             ts_xxxxxy_yyyyyy, \
                             ts_xxxxxy_yyyyyz, \
                             ts_xxxxxy_yyyyzz, \
                             ts_xxxxxy_yyyzzz, \
                             ts_xxxxxy_yyzzzz, \
                             ts_xxxxxy_yzzzzz, \
                             ts_xxxxxy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxy_xxxxxx[i] = tk_xxxxx_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxxx[i] * fz_0;

        tk_xxxxxy_xxxxxy[i] = tk_xxxxx_xxxxx[i] * fe_0 + tk_xxxxx_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxxy[i] * fz_0;

        tk_xxxxxy_xxxxxz[i] = tk_xxxxx_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxxz[i] * fz_0;

        tk_xxxxxy_xxxxyy[i] = 2.0 * tk_xxxxx_xxxxy[i] * fe_0 + tk_xxxxx_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxyy[i] * fz_0;

        tk_xxxxxy_xxxxyz[i] = tk_xxxxx_xxxxz[i] * fe_0 + tk_xxxxx_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxyz[i] * fz_0;

        tk_xxxxxy_xxxxzz[i] = tk_xxxxx_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxzz[i] * fz_0;

        tk_xxxxxy_xxxyyy[i] = 3.0 * tk_xxxxx_xxxyy[i] * fe_0 + tk_xxxxx_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxyyy[i] * fz_0;

        tk_xxxxxy_xxxyyz[i] = 2.0 * tk_xxxxx_xxxyz[i] * fe_0 + tk_xxxxx_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxyyz[i] * fz_0;

        tk_xxxxxy_xxxyzz[i] = tk_xxxxx_xxxzz[i] * fe_0 + tk_xxxxx_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxyzz[i] * fz_0;

        tk_xxxxxy_xxxzzz[i] = tk_xxxxx_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxzzz[i] * fz_0;

        tk_xxxxxy_xxyyyy[i] = 4.0 * tk_xxxxx_xxyyy[i] * fe_0 + tk_xxxxx_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyyyy[i] * fz_0;

        tk_xxxxxy_xxyyyz[i] = 3.0 * tk_xxxxx_xxyyz[i] * fe_0 + tk_xxxxx_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyyyz[i] * fz_0;

        tk_xxxxxy_xxyyzz[i] = 2.0 * tk_xxxxx_xxyzz[i] * fe_0 + tk_xxxxx_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyyzz[i] * fz_0;

        tk_xxxxxy_xxyzzz[i] = tk_xxxxx_xxzzz[i] * fe_0 + tk_xxxxx_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyzzz[i] * fz_0;

        tk_xxxxxy_xxzzzz[i] = tk_xxxxx_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxzzzz[i] * fz_0;

        tk_xxxxxy_xyyyyy[i] = 5.0 * tk_xxxxx_xyyyy[i] * fe_0 + tk_xxxxx_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyyyy[i] * fz_0;

        tk_xxxxxy_xyyyyz[i] = 4.0 * tk_xxxxx_xyyyz[i] * fe_0 + tk_xxxxx_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyyyz[i] * fz_0;

        tk_xxxxxy_xyyyzz[i] = 3.0 * tk_xxxxx_xyyzz[i] * fe_0 + tk_xxxxx_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyyzz[i] * fz_0;

        tk_xxxxxy_xyyzzz[i] = 2.0 * tk_xxxxx_xyzzz[i] * fe_0 + tk_xxxxx_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyzzz[i] * fz_0;

        tk_xxxxxy_xyzzzz[i] = tk_xxxxx_xzzzz[i] * fe_0 + tk_xxxxx_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyzzzz[i] * fz_0;

        tk_xxxxxy_xzzzzz[i] = tk_xxxxx_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xzzzzz[i] * fz_0;

        tk_xxxxxy_yyyyyy[i] = 6.0 * tk_xxxxx_yyyyy[i] * fe_0 + tk_xxxxx_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyyyy[i] * fz_0;

        tk_xxxxxy_yyyyyz[i] = 5.0 * tk_xxxxx_yyyyz[i] * fe_0 + tk_xxxxx_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyyyz[i] * fz_0;

        tk_xxxxxy_yyyyzz[i] = 4.0 * tk_xxxxx_yyyzz[i] * fe_0 + tk_xxxxx_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyyzz[i] * fz_0;

        tk_xxxxxy_yyyzzz[i] = 3.0 * tk_xxxxx_yyzzz[i] * fe_0 + tk_xxxxx_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyzzz[i] * fz_0;

        tk_xxxxxy_yyzzzz[i] = 2.0 * tk_xxxxx_yzzzz[i] * fe_0 + tk_xxxxx_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyzzzz[i] * fz_0;

        tk_xxxxxy_yzzzzz[i] = tk_xxxxx_zzzzz[i] * fe_0 + tk_xxxxx_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yzzzzz[i] * fz_0;

        tk_xxxxxy_zzzzzz[i] = tk_xxxxx_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_zzzzzz[i] * fz_0;
    }

    // Set up 56-84 components of targeted buffer : II

    auto tk_xxxxxz_xxxxxx = pbuffer.data(idx_kin_ii + 56);

    auto tk_xxxxxz_xxxxxy = pbuffer.data(idx_kin_ii + 57);

    auto tk_xxxxxz_xxxxxz = pbuffer.data(idx_kin_ii + 58);

    auto tk_xxxxxz_xxxxyy = pbuffer.data(idx_kin_ii + 59);

    auto tk_xxxxxz_xxxxyz = pbuffer.data(idx_kin_ii + 60);

    auto tk_xxxxxz_xxxxzz = pbuffer.data(idx_kin_ii + 61);

    auto tk_xxxxxz_xxxyyy = pbuffer.data(idx_kin_ii + 62);

    auto tk_xxxxxz_xxxyyz = pbuffer.data(idx_kin_ii + 63);

    auto tk_xxxxxz_xxxyzz = pbuffer.data(idx_kin_ii + 64);

    auto tk_xxxxxz_xxxzzz = pbuffer.data(idx_kin_ii + 65);

    auto tk_xxxxxz_xxyyyy = pbuffer.data(idx_kin_ii + 66);

    auto tk_xxxxxz_xxyyyz = pbuffer.data(idx_kin_ii + 67);

    auto tk_xxxxxz_xxyyzz = pbuffer.data(idx_kin_ii + 68);

    auto tk_xxxxxz_xxyzzz = pbuffer.data(idx_kin_ii + 69);

    auto tk_xxxxxz_xxzzzz = pbuffer.data(idx_kin_ii + 70);

    auto tk_xxxxxz_xyyyyy = pbuffer.data(idx_kin_ii + 71);

    auto tk_xxxxxz_xyyyyz = pbuffer.data(idx_kin_ii + 72);

    auto tk_xxxxxz_xyyyzz = pbuffer.data(idx_kin_ii + 73);

    auto tk_xxxxxz_xyyzzz = pbuffer.data(idx_kin_ii + 74);

    auto tk_xxxxxz_xyzzzz = pbuffer.data(idx_kin_ii + 75);

    auto tk_xxxxxz_xzzzzz = pbuffer.data(idx_kin_ii + 76);

    auto tk_xxxxxz_yyyyyy = pbuffer.data(idx_kin_ii + 77);

    auto tk_xxxxxz_yyyyyz = pbuffer.data(idx_kin_ii + 78);

    auto tk_xxxxxz_yyyyzz = pbuffer.data(idx_kin_ii + 79);

    auto tk_xxxxxz_yyyzzz = pbuffer.data(idx_kin_ii + 80);

    auto tk_xxxxxz_yyzzzz = pbuffer.data(idx_kin_ii + 81);

    auto tk_xxxxxz_yzzzzz = pbuffer.data(idx_kin_ii + 82);

    auto tk_xxxxxz_zzzzzz = pbuffer.data(idx_kin_ii + 83);

#pragma omp simd aligned(pa_z,                 \
                             tk_xxxxx_xxxxx,   \
                             tk_xxxxx_xxxxxx,  \
                             tk_xxxxx_xxxxxy,  \
                             tk_xxxxx_xxxxxz,  \
                             tk_xxxxx_xxxxy,   \
                             tk_xxxxx_xxxxyy,  \
                             tk_xxxxx_xxxxyz,  \
                             tk_xxxxx_xxxxz,   \
                             tk_xxxxx_xxxxzz,  \
                             tk_xxxxx_xxxyy,   \
                             tk_xxxxx_xxxyyy,  \
                             tk_xxxxx_xxxyyz,  \
                             tk_xxxxx_xxxyz,   \
                             tk_xxxxx_xxxyzz,  \
                             tk_xxxxx_xxxzz,   \
                             tk_xxxxx_xxxzzz,  \
                             tk_xxxxx_xxyyy,   \
                             tk_xxxxx_xxyyyy,  \
                             tk_xxxxx_xxyyyz,  \
                             tk_xxxxx_xxyyz,   \
                             tk_xxxxx_xxyyzz,  \
                             tk_xxxxx_xxyzz,   \
                             tk_xxxxx_xxyzzz,  \
                             tk_xxxxx_xxzzz,   \
                             tk_xxxxx_xxzzzz,  \
                             tk_xxxxx_xyyyy,   \
                             tk_xxxxx_xyyyyy,  \
                             tk_xxxxx_xyyyyz,  \
                             tk_xxxxx_xyyyz,   \
                             tk_xxxxx_xyyyzz,  \
                             tk_xxxxx_xyyzz,   \
                             tk_xxxxx_xyyzzz,  \
                             tk_xxxxx_xyzzz,   \
                             tk_xxxxx_xyzzzz,  \
                             tk_xxxxx_xzzzz,   \
                             tk_xxxxx_xzzzzz,  \
                             tk_xxxxx_yyyyy,   \
                             tk_xxxxx_yyyyyy,  \
                             tk_xxxxx_yyyyyz,  \
                             tk_xxxxx_yyyyz,   \
                             tk_xxxxx_yyyyzz,  \
                             tk_xxxxx_yyyzz,   \
                             tk_xxxxx_yyyzzz,  \
                             tk_xxxxx_yyzzz,   \
                             tk_xxxxx_yyzzzz,  \
                             tk_xxxxx_yzzzz,   \
                             tk_xxxxx_yzzzzz,  \
                             tk_xxxxx_zzzzz,   \
                             tk_xxxxx_zzzzzz,  \
                             tk_xxxxxz_xxxxxx, \
                             tk_xxxxxz_xxxxxy, \
                             tk_xxxxxz_xxxxxz, \
                             tk_xxxxxz_xxxxyy, \
                             tk_xxxxxz_xxxxyz, \
                             tk_xxxxxz_xxxxzz, \
                             tk_xxxxxz_xxxyyy, \
                             tk_xxxxxz_xxxyyz, \
                             tk_xxxxxz_xxxyzz, \
                             tk_xxxxxz_xxxzzz, \
                             tk_xxxxxz_xxyyyy, \
                             tk_xxxxxz_xxyyyz, \
                             tk_xxxxxz_xxyyzz, \
                             tk_xxxxxz_xxyzzz, \
                             tk_xxxxxz_xxzzzz, \
                             tk_xxxxxz_xyyyyy, \
                             tk_xxxxxz_xyyyyz, \
                             tk_xxxxxz_xyyyzz, \
                             tk_xxxxxz_xyyzzz, \
                             tk_xxxxxz_xyzzzz, \
                             tk_xxxxxz_xzzzzz, \
                             tk_xxxxxz_yyyyyy, \
                             tk_xxxxxz_yyyyyz, \
                             tk_xxxxxz_yyyyzz, \
                             tk_xxxxxz_yyyzzz, \
                             tk_xxxxxz_yyzzzz, \
                             tk_xxxxxz_yzzzzz, \
                             tk_xxxxxz_zzzzzz, \
                             ts_xxxxxz_xxxxxx, \
                             ts_xxxxxz_xxxxxy, \
                             ts_xxxxxz_xxxxxz, \
                             ts_xxxxxz_xxxxyy, \
                             ts_xxxxxz_xxxxyz, \
                             ts_xxxxxz_xxxxzz, \
                             ts_xxxxxz_xxxyyy, \
                             ts_xxxxxz_xxxyyz, \
                             ts_xxxxxz_xxxyzz, \
                             ts_xxxxxz_xxxzzz, \
                             ts_xxxxxz_xxyyyy, \
                             ts_xxxxxz_xxyyyz, \
                             ts_xxxxxz_xxyyzz, \
                             ts_xxxxxz_xxyzzz, \
                             ts_xxxxxz_xxzzzz, \
                             ts_xxxxxz_xyyyyy, \
                             ts_xxxxxz_xyyyyz, \
                             ts_xxxxxz_xyyyzz, \
                             ts_xxxxxz_xyyzzz, \
                             ts_xxxxxz_xyzzzz, \
                             ts_xxxxxz_xzzzzz, \
                             ts_xxxxxz_yyyyyy, \
                             ts_xxxxxz_yyyyyz, \
                             ts_xxxxxz_yyyyzz, \
                             ts_xxxxxz_yyyzzz, \
                             ts_xxxxxz_yyzzzz, \
                             ts_xxxxxz_yzzzzz, \
                             ts_xxxxxz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxz_xxxxxx[i] = tk_xxxxx_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxxx[i] * fz_0;

        tk_xxxxxz_xxxxxy[i] = tk_xxxxx_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxxy[i] * fz_0;

        tk_xxxxxz_xxxxxz[i] = tk_xxxxx_xxxxx[i] * fe_0 + tk_xxxxx_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxxz[i] * fz_0;

        tk_xxxxxz_xxxxyy[i] = tk_xxxxx_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxyy[i] * fz_0;

        tk_xxxxxz_xxxxyz[i] = tk_xxxxx_xxxxy[i] * fe_0 + tk_xxxxx_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxyz[i] * fz_0;

        tk_xxxxxz_xxxxzz[i] = 2.0 * tk_xxxxx_xxxxz[i] * fe_0 + tk_xxxxx_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxzz[i] * fz_0;

        tk_xxxxxz_xxxyyy[i] = tk_xxxxx_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxyyy[i] * fz_0;

        tk_xxxxxz_xxxyyz[i] = tk_xxxxx_xxxyy[i] * fe_0 + tk_xxxxx_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxyyz[i] * fz_0;

        tk_xxxxxz_xxxyzz[i] = 2.0 * tk_xxxxx_xxxyz[i] * fe_0 + tk_xxxxx_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxyzz[i] * fz_0;

        tk_xxxxxz_xxxzzz[i] = 3.0 * tk_xxxxx_xxxzz[i] * fe_0 + tk_xxxxx_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxzzz[i] * fz_0;

        tk_xxxxxz_xxyyyy[i] = tk_xxxxx_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyyyy[i] * fz_0;

        tk_xxxxxz_xxyyyz[i] = tk_xxxxx_xxyyy[i] * fe_0 + tk_xxxxx_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyyyz[i] * fz_0;

        tk_xxxxxz_xxyyzz[i] = 2.0 * tk_xxxxx_xxyyz[i] * fe_0 + tk_xxxxx_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyyzz[i] * fz_0;

        tk_xxxxxz_xxyzzz[i] = 3.0 * tk_xxxxx_xxyzz[i] * fe_0 + tk_xxxxx_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyzzz[i] * fz_0;

        tk_xxxxxz_xxzzzz[i] = 4.0 * tk_xxxxx_xxzzz[i] * fe_0 + tk_xxxxx_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxzzzz[i] * fz_0;

        tk_xxxxxz_xyyyyy[i] = tk_xxxxx_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyyyy[i] * fz_0;

        tk_xxxxxz_xyyyyz[i] = tk_xxxxx_xyyyy[i] * fe_0 + tk_xxxxx_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyyyz[i] * fz_0;

        tk_xxxxxz_xyyyzz[i] = 2.0 * tk_xxxxx_xyyyz[i] * fe_0 + tk_xxxxx_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyyzz[i] * fz_0;

        tk_xxxxxz_xyyzzz[i] = 3.0 * tk_xxxxx_xyyzz[i] * fe_0 + tk_xxxxx_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyzzz[i] * fz_0;

        tk_xxxxxz_xyzzzz[i] = 4.0 * tk_xxxxx_xyzzz[i] * fe_0 + tk_xxxxx_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyzzzz[i] * fz_0;

        tk_xxxxxz_xzzzzz[i] = 5.0 * tk_xxxxx_xzzzz[i] * fe_0 + tk_xxxxx_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xzzzzz[i] * fz_0;

        tk_xxxxxz_yyyyyy[i] = tk_xxxxx_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyyyy[i] * fz_0;

        tk_xxxxxz_yyyyyz[i] = tk_xxxxx_yyyyy[i] * fe_0 + tk_xxxxx_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyyyz[i] * fz_0;

        tk_xxxxxz_yyyyzz[i] = 2.0 * tk_xxxxx_yyyyz[i] * fe_0 + tk_xxxxx_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyyzz[i] * fz_0;

        tk_xxxxxz_yyyzzz[i] = 3.0 * tk_xxxxx_yyyzz[i] * fe_0 + tk_xxxxx_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyzzz[i] * fz_0;

        tk_xxxxxz_yyzzzz[i] = 4.0 * tk_xxxxx_yyzzz[i] * fe_0 + tk_xxxxx_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyzzzz[i] * fz_0;

        tk_xxxxxz_yzzzzz[i] = 5.0 * tk_xxxxx_yzzzz[i] * fe_0 + tk_xxxxx_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yzzzzz[i] * fz_0;

        tk_xxxxxz_zzzzzz[i] = 6.0 * tk_xxxxx_zzzzz[i] * fe_0 + tk_xxxxx_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_zzzzzz[i] * fz_0;
    }

    // Set up 84-112 components of targeted buffer : II

    auto tk_xxxxyy_xxxxxx = pbuffer.data(idx_kin_ii + 84);

    auto tk_xxxxyy_xxxxxy = pbuffer.data(idx_kin_ii + 85);

    auto tk_xxxxyy_xxxxxz = pbuffer.data(idx_kin_ii + 86);

    auto tk_xxxxyy_xxxxyy = pbuffer.data(idx_kin_ii + 87);

    auto tk_xxxxyy_xxxxyz = pbuffer.data(idx_kin_ii + 88);

    auto tk_xxxxyy_xxxxzz = pbuffer.data(idx_kin_ii + 89);

    auto tk_xxxxyy_xxxyyy = pbuffer.data(idx_kin_ii + 90);

    auto tk_xxxxyy_xxxyyz = pbuffer.data(idx_kin_ii + 91);

    auto tk_xxxxyy_xxxyzz = pbuffer.data(idx_kin_ii + 92);

    auto tk_xxxxyy_xxxzzz = pbuffer.data(idx_kin_ii + 93);

    auto tk_xxxxyy_xxyyyy = pbuffer.data(idx_kin_ii + 94);

    auto tk_xxxxyy_xxyyyz = pbuffer.data(idx_kin_ii + 95);

    auto tk_xxxxyy_xxyyzz = pbuffer.data(idx_kin_ii + 96);

    auto tk_xxxxyy_xxyzzz = pbuffer.data(idx_kin_ii + 97);

    auto tk_xxxxyy_xxzzzz = pbuffer.data(idx_kin_ii + 98);

    auto tk_xxxxyy_xyyyyy = pbuffer.data(idx_kin_ii + 99);

    auto tk_xxxxyy_xyyyyz = pbuffer.data(idx_kin_ii + 100);

    auto tk_xxxxyy_xyyyzz = pbuffer.data(idx_kin_ii + 101);

    auto tk_xxxxyy_xyyzzz = pbuffer.data(idx_kin_ii + 102);

    auto tk_xxxxyy_xyzzzz = pbuffer.data(idx_kin_ii + 103);

    auto tk_xxxxyy_xzzzzz = pbuffer.data(idx_kin_ii + 104);

    auto tk_xxxxyy_yyyyyy = pbuffer.data(idx_kin_ii + 105);

    auto tk_xxxxyy_yyyyyz = pbuffer.data(idx_kin_ii + 106);

    auto tk_xxxxyy_yyyyzz = pbuffer.data(idx_kin_ii + 107);

    auto tk_xxxxyy_yyyzzz = pbuffer.data(idx_kin_ii + 108);

    auto tk_xxxxyy_yyzzzz = pbuffer.data(idx_kin_ii + 109);

    auto tk_xxxxyy_yzzzzz = pbuffer.data(idx_kin_ii + 110);

    auto tk_xxxxyy_zzzzzz = pbuffer.data(idx_kin_ii + 111);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tk_xxxx_xxxxxx,   \
                             tk_xxxx_xxxxxz,   \
                             tk_xxxx_xxxxzz,   \
                             tk_xxxx_xxxzzz,   \
                             tk_xxxx_xxzzzz,   \
                             tk_xxxx_xzzzzz,   \
                             tk_xxxxy_xxxxxx,  \
                             tk_xxxxy_xxxxxz,  \
                             tk_xxxxy_xxxxzz,  \
                             tk_xxxxy_xxxzzz,  \
                             tk_xxxxy_xxzzzz,  \
                             tk_xxxxy_xzzzzz,  \
                             tk_xxxxyy_xxxxxx, \
                             tk_xxxxyy_xxxxxy, \
                             tk_xxxxyy_xxxxxz, \
                             tk_xxxxyy_xxxxyy, \
                             tk_xxxxyy_xxxxyz, \
                             tk_xxxxyy_xxxxzz, \
                             tk_xxxxyy_xxxyyy, \
                             tk_xxxxyy_xxxyyz, \
                             tk_xxxxyy_xxxyzz, \
                             tk_xxxxyy_xxxzzz, \
                             tk_xxxxyy_xxyyyy, \
                             tk_xxxxyy_xxyyyz, \
                             tk_xxxxyy_xxyyzz, \
                             tk_xxxxyy_xxyzzz, \
                             tk_xxxxyy_xxzzzz, \
                             tk_xxxxyy_xyyyyy, \
                             tk_xxxxyy_xyyyyz, \
                             tk_xxxxyy_xyyyzz, \
                             tk_xxxxyy_xyyzzz, \
                             tk_xxxxyy_xyzzzz, \
                             tk_xxxxyy_xzzzzz, \
                             tk_xxxxyy_yyyyyy, \
                             tk_xxxxyy_yyyyyz, \
                             tk_xxxxyy_yyyyzz, \
                             tk_xxxxyy_yyyzzz, \
                             tk_xxxxyy_yyzzzz, \
                             tk_xxxxyy_yzzzzz, \
                             tk_xxxxyy_zzzzzz, \
                             tk_xxxyy_xxxxxy,  \
                             tk_xxxyy_xxxxy,   \
                             tk_xxxyy_xxxxyy,  \
                             tk_xxxyy_xxxxyz,  \
                             tk_xxxyy_xxxyy,   \
                             tk_xxxyy_xxxyyy,  \
                             tk_xxxyy_xxxyyz,  \
                             tk_xxxyy_xxxyz,   \
                             tk_xxxyy_xxxyzz,  \
                             tk_xxxyy_xxyyy,   \
                             tk_xxxyy_xxyyyy,  \
                             tk_xxxyy_xxyyyz,  \
                             tk_xxxyy_xxyyz,   \
                             tk_xxxyy_xxyyzz,  \
                             tk_xxxyy_xxyzz,   \
                             tk_xxxyy_xxyzzz,  \
                             tk_xxxyy_xyyyy,   \
                             tk_xxxyy_xyyyyy,  \
                             tk_xxxyy_xyyyyz,  \
                             tk_xxxyy_xyyyz,   \
                             tk_xxxyy_xyyyzz,  \
                             tk_xxxyy_xyyzz,   \
                             tk_xxxyy_xyyzzz,  \
                             tk_xxxyy_xyzzz,   \
                             tk_xxxyy_xyzzzz,  \
                             tk_xxxyy_yyyyy,   \
                             tk_xxxyy_yyyyyy,  \
                             tk_xxxyy_yyyyyz,  \
                             tk_xxxyy_yyyyz,   \
                             tk_xxxyy_yyyyzz,  \
                             tk_xxxyy_yyyzz,   \
                             tk_xxxyy_yyyzzz,  \
                             tk_xxxyy_yyzzz,   \
                             tk_xxxyy_yyzzzz,  \
                             tk_xxxyy_yzzzz,   \
                             tk_xxxyy_yzzzzz,  \
                             tk_xxxyy_zzzzzz,  \
                             tk_xxyy_xxxxxy,   \
                             tk_xxyy_xxxxyy,   \
                             tk_xxyy_xxxxyz,   \
                             tk_xxyy_xxxyyy,   \
                             tk_xxyy_xxxyyz,   \
                             tk_xxyy_xxxyzz,   \
                             tk_xxyy_xxyyyy,   \
                             tk_xxyy_xxyyyz,   \
                             tk_xxyy_xxyyzz,   \
                             tk_xxyy_xxyzzz,   \
                             tk_xxyy_xyyyyy,   \
                             tk_xxyy_xyyyyz,   \
                             tk_xxyy_xyyyzz,   \
                             tk_xxyy_xyyzzz,   \
                             tk_xxyy_xyzzzz,   \
                             tk_xxyy_yyyyyy,   \
                             tk_xxyy_yyyyyz,   \
                             tk_xxyy_yyyyzz,   \
                             tk_xxyy_yyyzzz,   \
                             tk_xxyy_yyzzzz,   \
                             tk_xxyy_yzzzzz,   \
                             tk_xxyy_zzzzzz,   \
                             ts_xxxx_xxxxxx,   \
                             ts_xxxx_xxxxxz,   \
                             ts_xxxx_xxxxzz,   \
                             ts_xxxx_xxxzzz,   \
                             ts_xxxx_xxzzzz,   \
                             ts_xxxx_xzzzzz,   \
                             ts_xxxxyy_xxxxxx, \
                             ts_xxxxyy_xxxxxy, \
                             ts_xxxxyy_xxxxxz, \
                             ts_xxxxyy_xxxxyy, \
                             ts_xxxxyy_xxxxyz, \
                             ts_xxxxyy_xxxxzz, \
                             ts_xxxxyy_xxxyyy, \
                             ts_xxxxyy_xxxyyz, \
                             ts_xxxxyy_xxxyzz, \
                             ts_xxxxyy_xxxzzz, \
                             ts_xxxxyy_xxyyyy, \
                             ts_xxxxyy_xxyyyz, \
                             ts_xxxxyy_xxyyzz, \
                             ts_xxxxyy_xxyzzz, \
                             ts_xxxxyy_xxzzzz, \
                             ts_xxxxyy_xyyyyy, \
                             ts_xxxxyy_xyyyyz, \
                             ts_xxxxyy_xyyyzz, \
                             ts_xxxxyy_xyyzzz, \
                             ts_xxxxyy_xyzzzz, \
                             ts_xxxxyy_xzzzzz, \
                             ts_xxxxyy_yyyyyy, \
                             ts_xxxxyy_yyyyyz, \
                             ts_xxxxyy_yyyyzz, \
                             ts_xxxxyy_yyyzzz, \
                             ts_xxxxyy_yyzzzz, \
                             ts_xxxxyy_yzzzzz, \
                             ts_xxxxyy_zzzzzz, \
                             ts_xxyy_xxxxxy,   \
                             ts_xxyy_xxxxyy,   \
                             ts_xxyy_xxxxyz,   \
                             ts_xxyy_xxxyyy,   \
                             ts_xxyy_xxxyyz,   \
                             ts_xxyy_xxxyzz,   \
                             ts_xxyy_xxyyyy,   \
                             ts_xxyy_xxyyyz,   \
                             ts_xxyy_xxyyzz,   \
                             ts_xxyy_xxyzzz,   \
                             ts_xxyy_xyyyyy,   \
                             ts_xxyy_xyyyyz,   \
                             ts_xxyy_xyyyzz,   \
                             ts_xxyy_xyyzzz,   \
                             ts_xxyy_xyzzzz,   \
                             ts_xxyy_yyyyyy,   \
                             ts_xxyy_yyyyyz,   \
                             ts_xxyy_yyyyzz,   \
                             ts_xxyy_yyyzzz,   \
                             ts_xxyy_yyzzzz,   \
                             ts_xxyy_yzzzzz,   \
                             ts_xxyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxyy_xxxxxx[i] =
            -2.0 * ts_xxxx_xxxxxx[i] * fbe_0 * fz_0 + tk_xxxx_xxxxxx[i] * fe_0 + tk_xxxxy_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxxxx[i] * fz_0;

        tk_xxxxyy_xxxxxy[i] = -6.0 * ts_xxyy_xxxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxxy[i] * fe_0 + 5.0 * tk_xxxyy_xxxxy[i] * fe_0 +
                              tk_xxxyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxxxy[i] * fz_0;

        tk_xxxxyy_xxxxxz[i] =
            -2.0 * ts_xxxx_xxxxxz[i] * fbe_0 * fz_0 + tk_xxxx_xxxxxz[i] * fe_0 + tk_xxxxy_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxxxz[i] * fz_0;

        tk_xxxxyy_xxxxyy[i] = -6.0 * ts_xxyy_xxxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxyy[i] * fe_0 + 4.0 * tk_xxxyy_xxxyy[i] * fe_0 +
                              tk_xxxyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxxyy[i] * fz_0;

        tk_xxxxyy_xxxxyz[i] = -6.0 * ts_xxyy_xxxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxyz[i] * fe_0 + 4.0 * tk_xxxyy_xxxyz[i] * fe_0 +
                              tk_xxxyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxxyz[i] * fz_0;

        tk_xxxxyy_xxxxzz[i] =
            -2.0 * ts_xxxx_xxxxzz[i] * fbe_0 * fz_0 + tk_xxxx_xxxxzz[i] * fe_0 + tk_xxxxy_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxxzz[i] * fz_0;

        tk_xxxxyy_xxxyyy[i] = -6.0 * ts_xxyy_xxxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxyyy[i] * fe_0 + 3.0 * tk_xxxyy_xxyyy[i] * fe_0 +
                              tk_xxxyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxyyy[i] * fz_0;

        tk_xxxxyy_xxxyyz[i] = -6.0 * ts_xxyy_xxxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxyyz[i] * fe_0 + 3.0 * tk_xxxyy_xxyyz[i] * fe_0 +
                              tk_xxxyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxyyz[i] * fz_0;

        tk_xxxxyy_xxxyzz[i] = -6.0 * ts_xxyy_xxxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxyzz[i] * fe_0 + 3.0 * tk_xxxyy_xxyzz[i] * fe_0 +
                              tk_xxxyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxyzz[i] * fz_0;

        tk_xxxxyy_xxxzzz[i] =
            -2.0 * ts_xxxx_xxxzzz[i] * fbe_0 * fz_0 + tk_xxxx_xxxzzz[i] * fe_0 + tk_xxxxy_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxzzz[i] * fz_0;

        tk_xxxxyy_xxyyyy[i] = -6.0 * ts_xxyy_xxyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyyyy[i] * fe_0 + 2.0 * tk_xxxyy_xyyyy[i] * fe_0 +
                              tk_xxxyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyyyy[i] * fz_0;

        tk_xxxxyy_xxyyyz[i] = -6.0 * ts_xxyy_xxyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyyyz[i] * fe_0 + 2.0 * tk_xxxyy_xyyyz[i] * fe_0 +
                              tk_xxxyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyyyz[i] * fz_0;

        tk_xxxxyy_xxyyzz[i] = -6.0 * ts_xxyy_xxyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyyzz[i] * fe_0 + 2.0 * tk_xxxyy_xyyzz[i] * fe_0 +
                              tk_xxxyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyyzz[i] * fz_0;

        tk_xxxxyy_xxyzzz[i] = -6.0 * ts_xxyy_xxyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyzzz[i] * fe_0 + 2.0 * tk_xxxyy_xyzzz[i] * fe_0 +
                              tk_xxxyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyzzz[i] * fz_0;

        tk_xxxxyy_xxzzzz[i] =
            -2.0 * ts_xxxx_xxzzzz[i] * fbe_0 * fz_0 + tk_xxxx_xxzzzz[i] * fe_0 + tk_xxxxy_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxzzzz[i] * fz_0;

        tk_xxxxyy_xyyyyy[i] = -6.0 * ts_xxyy_xyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyyyy[i] * fe_0 + tk_xxxyy_yyyyy[i] * fe_0 +
                              tk_xxxyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyyyy[i] * fz_0;

        tk_xxxxyy_xyyyyz[i] = -6.0 * ts_xxyy_xyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyyyz[i] * fe_0 + tk_xxxyy_yyyyz[i] * fe_0 +
                              tk_xxxyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyyyz[i] * fz_0;

        tk_xxxxyy_xyyyzz[i] = -6.0 * ts_xxyy_xyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyyzz[i] * fe_0 + tk_xxxyy_yyyzz[i] * fe_0 +
                              tk_xxxyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyyzz[i] * fz_0;

        tk_xxxxyy_xyyzzz[i] = -6.0 * ts_xxyy_xyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyzzz[i] * fe_0 + tk_xxxyy_yyzzz[i] * fe_0 +
                              tk_xxxyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyzzz[i] * fz_0;

        tk_xxxxyy_xyzzzz[i] = -6.0 * ts_xxyy_xyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyzzzz[i] * fe_0 + tk_xxxyy_yzzzz[i] * fe_0 +
                              tk_xxxyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyzzzz[i] * fz_0;

        tk_xxxxyy_xzzzzz[i] =
            -2.0 * ts_xxxx_xzzzzz[i] * fbe_0 * fz_0 + tk_xxxx_xzzzzz[i] * fe_0 + tk_xxxxy_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xzzzzz[i] * fz_0;

        tk_xxxxyy_yyyyyy[i] = -6.0 * ts_xxyy_yyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyyyy[i] * fe_0 + tk_xxxyy_yyyyyy[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_yyyyyy[i] * fz_0;

        tk_xxxxyy_yyyyyz[i] = -6.0 * ts_xxyy_yyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyyyz[i] * fe_0 + tk_xxxyy_yyyyyz[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_yyyyyz[i] * fz_0;

        tk_xxxxyy_yyyyzz[i] = -6.0 * ts_xxyy_yyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyyzz[i] * fe_0 + tk_xxxyy_yyyyzz[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_yyyyzz[i] * fz_0;

        tk_xxxxyy_yyyzzz[i] = -6.0 * ts_xxyy_yyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyzzz[i] * fe_0 + tk_xxxyy_yyyzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_yyyzzz[i] * fz_0;

        tk_xxxxyy_yyzzzz[i] = -6.0 * ts_xxyy_yyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyzzzz[i] * fe_0 + tk_xxxyy_yyzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_yyzzzz[i] * fz_0;

        tk_xxxxyy_yzzzzz[i] = -6.0 * ts_xxyy_yzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yzzzzz[i] * fe_0 + tk_xxxyy_yzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_yzzzzz[i] * fz_0;

        tk_xxxxyy_zzzzzz[i] = -6.0 * ts_xxyy_zzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_zzzzzz[i] * fe_0 + tk_xxxyy_zzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxyy_zzzzzz[i] * fz_0;
    }

    // Set up 112-140 components of targeted buffer : II

    auto tk_xxxxyz_xxxxxx = pbuffer.data(idx_kin_ii + 112);

    auto tk_xxxxyz_xxxxxy = pbuffer.data(idx_kin_ii + 113);

    auto tk_xxxxyz_xxxxxz = pbuffer.data(idx_kin_ii + 114);

    auto tk_xxxxyz_xxxxyy = pbuffer.data(idx_kin_ii + 115);

    auto tk_xxxxyz_xxxxyz = pbuffer.data(idx_kin_ii + 116);

    auto tk_xxxxyz_xxxxzz = pbuffer.data(idx_kin_ii + 117);

    auto tk_xxxxyz_xxxyyy = pbuffer.data(idx_kin_ii + 118);

    auto tk_xxxxyz_xxxyyz = pbuffer.data(idx_kin_ii + 119);

    auto tk_xxxxyz_xxxyzz = pbuffer.data(idx_kin_ii + 120);

    auto tk_xxxxyz_xxxzzz = pbuffer.data(idx_kin_ii + 121);

    auto tk_xxxxyz_xxyyyy = pbuffer.data(idx_kin_ii + 122);

    auto tk_xxxxyz_xxyyyz = pbuffer.data(idx_kin_ii + 123);

    auto tk_xxxxyz_xxyyzz = pbuffer.data(idx_kin_ii + 124);

    auto tk_xxxxyz_xxyzzz = pbuffer.data(idx_kin_ii + 125);

    auto tk_xxxxyz_xxzzzz = pbuffer.data(idx_kin_ii + 126);

    auto tk_xxxxyz_xyyyyy = pbuffer.data(idx_kin_ii + 127);

    auto tk_xxxxyz_xyyyyz = pbuffer.data(idx_kin_ii + 128);

    auto tk_xxxxyz_xyyyzz = pbuffer.data(idx_kin_ii + 129);

    auto tk_xxxxyz_xyyzzz = pbuffer.data(idx_kin_ii + 130);

    auto tk_xxxxyz_xyzzzz = pbuffer.data(idx_kin_ii + 131);

    auto tk_xxxxyz_xzzzzz = pbuffer.data(idx_kin_ii + 132);

    auto tk_xxxxyz_yyyyyy = pbuffer.data(idx_kin_ii + 133);

    auto tk_xxxxyz_yyyyyz = pbuffer.data(idx_kin_ii + 134);

    auto tk_xxxxyz_yyyyzz = pbuffer.data(idx_kin_ii + 135);

    auto tk_xxxxyz_yyyzzz = pbuffer.data(idx_kin_ii + 136);

    auto tk_xxxxyz_yyzzzz = pbuffer.data(idx_kin_ii + 137);

    auto tk_xxxxyz_yzzzzz = pbuffer.data(idx_kin_ii + 138);

    auto tk_xxxxyz_zzzzzz = pbuffer.data(idx_kin_ii + 139);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tk_xxxxy_xxxxxy,  \
                             tk_xxxxy_xxxxyy,  \
                             tk_xxxxy_xxxyyy,  \
                             tk_xxxxy_xxyyyy,  \
                             tk_xxxxy_xyyyyy,  \
                             tk_xxxxy_yyyyyy,  \
                             tk_xxxxyz_xxxxxx, \
                             tk_xxxxyz_xxxxxy, \
                             tk_xxxxyz_xxxxxz, \
                             tk_xxxxyz_xxxxyy, \
                             tk_xxxxyz_xxxxyz, \
                             tk_xxxxyz_xxxxzz, \
                             tk_xxxxyz_xxxyyy, \
                             tk_xxxxyz_xxxyyz, \
                             tk_xxxxyz_xxxyzz, \
                             tk_xxxxyz_xxxzzz, \
                             tk_xxxxyz_xxyyyy, \
                             tk_xxxxyz_xxyyyz, \
                             tk_xxxxyz_xxyyzz, \
                             tk_xxxxyz_xxyzzz, \
                             tk_xxxxyz_xxzzzz, \
                             tk_xxxxyz_xyyyyy, \
                             tk_xxxxyz_xyyyyz, \
                             tk_xxxxyz_xyyyzz, \
                             tk_xxxxyz_xyyzzz, \
                             tk_xxxxyz_xyzzzz, \
                             tk_xxxxyz_xzzzzz, \
                             tk_xxxxyz_yyyyyy, \
                             tk_xxxxyz_yyyyyz, \
                             tk_xxxxyz_yyyyzz, \
                             tk_xxxxyz_yyyzzz, \
                             tk_xxxxyz_yyzzzz, \
                             tk_xxxxyz_yzzzzz, \
                             tk_xxxxyz_zzzzzz, \
                             tk_xxxxz_xxxxxx,  \
                             tk_xxxxz_xxxxxz,  \
                             tk_xxxxz_xxxxyz,  \
                             tk_xxxxz_xxxxz,   \
                             tk_xxxxz_xxxxzz,  \
                             tk_xxxxz_xxxyyz,  \
                             tk_xxxxz_xxxyz,   \
                             tk_xxxxz_xxxyzz,  \
                             tk_xxxxz_xxxzz,   \
                             tk_xxxxz_xxxzzz,  \
                             tk_xxxxz_xxyyyz,  \
                             tk_xxxxz_xxyyz,   \
                             tk_xxxxz_xxyyzz,  \
                             tk_xxxxz_xxyzz,   \
                             tk_xxxxz_xxyzzz,  \
                             tk_xxxxz_xxzzz,   \
                             tk_xxxxz_xxzzzz,  \
                             tk_xxxxz_xyyyyz,  \
                             tk_xxxxz_xyyyz,   \
                             tk_xxxxz_xyyyzz,  \
                             tk_xxxxz_xyyzz,   \
                             tk_xxxxz_xyyzzz,  \
                             tk_xxxxz_xyzzz,   \
                             tk_xxxxz_xyzzzz,  \
                             tk_xxxxz_xzzzz,   \
                             tk_xxxxz_xzzzzz,  \
                             tk_xxxxz_yyyyyz,  \
                             tk_xxxxz_yyyyz,   \
                             tk_xxxxz_yyyyzz,  \
                             tk_xxxxz_yyyzz,   \
                             tk_xxxxz_yyyzzz,  \
                             tk_xxxxz_yyzzz,   \
                             tk_xxxxz_yyzzzz,  \
                             tk_xxxxz_yzzzz,   \
                             tk_xxxxz_yzzzzz,  \
                             tk_xxxxz_zzzzz,   \
                             tk_xxxxz_zzzzzz,  \
                             ts_xxxxyz_xxxxxx, \
                             ts_xxxxyz_xxxxxy, \
                             ts_xxxxyz_xxxxxz, \
                             ts_xxxxyz_xxxxyy, \
                             ts_xxxxyz_xxxxyz, \
                             ts_xxxxyz_xxxxzz, \
                             ts_xxxxyz_xxxyyy, \
                             ts_xxxxyz_xxxyyz, \
                             ts_xxxxyz_xxxyzz, \
                             ts_xxxxyz_xxxzzz, \
                             ts_xxxxyz_xxyyyy, \
                             ts_xxxxyz_xxyyyz, \
                             ts_xxxxyz_xxyyzz, \
                             ts_xxxxyz_xxyzzz, \
                             ts_xxxxyz_xxzzzz, \
                             ts_xxxxyz_xyyyyy, \
                             ts_xxxxyz_xyyyyz, \
                             ts_xxxxyz_xyyyzz, \
                             ts_xxxxyz_xyyzzz, \
                             ts_xxxxyz_xyzzzz, \
                             ts_xxxxyz_xzzzzz, \
                             ts_xxxxyz_yyyyyy, \
                             ts_xxxxyz_yyyyyz, \
                             ts_xxxxyz_yyyyzz, \
                             ts_xxxxyz_yyyzzz, \
                             ts_xxxxyz_yyzzzz, \
                             ts_xxxxyz_yzzzzz, \
                             ts_xxxxyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxyz_xxxxxx[i] = tk_xxxxz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxxxx[i] * fz_0;

        tk_xxxxyz_xxxxxy[i] = tk_xxxxy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxxxxy[i] * fz_0;

        tk_xxxxyz_xxxxxz[i] = tk_xxxxz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxxxz[i] * fz_0;

        tk_xxxxyz_xxxxyy[i] = tk_xxxxy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxxxyy[i] * fz_0;

        tk_xxxxyz_xxxxyz[i] = tk_xxxxz_xxxxz[i] * fe_0 + tk_xxxxz_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxxyz[i] * fz_0;

        tk_xxxxyz_xxxxzz[i] = tk_xxxxz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxxzz[i] * fz_0;

        tk_xxxxyz_xxxyyy[i] = tk_xxxxy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxxyyy[i] * fz_0;

        tk_xxxxyz_xxxyyz[i] = 2.0 * tk_xxxxz_xxxyz[i] * fe_0 + tk_xxxxz_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxyyz[i] * fz_0;

        tk_xxxxyz_xxxyzz[i] = tk_xxxxz_xxxzz[i] * fe_0 + tk_xxxxz_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxyzz[i] * fz_0;

        tk_xxxxyz_xxxzzz[i] = tk_xxxxz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxzzz[i] * fz_0;

        tk_xxxxyz_xxyyyy[i] = tk_xxxxy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxyyyy[i] * fz_0;

        tk_xxxxyz_xxyyyz[i] = 3.0 * tk_xxxxz_xxyyz[i] * fe_0 + tk_xxxxz_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxyyyz[i] * fz_0;

        tk_xxxxyz_xxyyzz[i] = 2.0 * tk_xxxxz_xxyzz[i] * fe_0 + tk_xxxxz_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxyyzz[i] * fz_0;

        tk_xxxxyz_xxyzzz[i] = tk_xxxxz_xxzzz[i] * fe_0 + tk_xxxxz_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxyzzz[i] * fz_0;

        tk_xxxxyz_xxzzzz[i] = tk_xxxxz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxzzzz[i] * fz_0;

        tk_xxxxyz_xyyyyy[i] = tk_xxxxy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xyyyyy[i] * fz_0;

        tk_xxxxyz_xyyyyz[i] = 4.0 * tk_xxxxz_xyyyz[i] * fe_0 + tk_xxxxz_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyyyyz[i] * fz_0;

        tk_xxxxyz_xyyyzz[i] = 3.0 * tk_xxxxz_xyyzz[i] * fe_0 + tk_xxxxz_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyyyzz[i] * fz_0;

        tk_xxxxyz_xyyzzz[i] = 2.0 * tk_xxxxz_xyzzz[i] * fe_0 + tk_xxxxz_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyyzzz[i] * fz_0;

        tk_xxxxyz_xyzzzz[i] = tk_xxxxz_xzzzz[i] * fe_0 + tk_xxxxz_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyzzzz[i] * fz_0;

        tk_xxxxyz_xzzzzz[i] = tk_xxxxz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xzzzzz[i] * fz_0;

        tk_xxxxyz_yyyyyy[i] = tk_xxxxy_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_yyyyyy[i] * fz_0;

        tk_xxxxyz_yyyyyz[i] = 5.0 * tk_xxxxz_yyyyz[i] * fe_0 + tk_xxxxz_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyyyyz[i] * fz_0;

        tk_xxxxyz_yyyyzz[i] = 4.0 * tk_xxxxz_yyyzz[i] * fe_0 + tk_xxxxz_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyyyzz[i] * fz_0;

        tk_xxxxyz_yyyzzz[i] = 3.0 * tk_xxxxz_yyzzz[i] * fe_0 + tk_xxxxz_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyyzzz[i] * fz_0;

        tk_xxxxyz_yyzzzz[i] = 2.0 * tk_xxxxz_yzzzz[i] * fe_0 + tk_xxxxz_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyzzzz[i] * fz_0;

        tk_xxxxyz_yzzzzz[i] = tk_xxxxz_zzzzz[i] * fe_0 + tk_xxxxz_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yzzzzz[i] * fz_0;

        tk_xxxxyz_zzzzzz[i] = tk_xxxxz_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_zzzzzz[i] * fz_0;
    }

    // Set up 140-168 components of targeted buffer : II

    auto tk_xxxxzz_xxxxxx = pbuffer.data(idx_kin_ii + 140);

    auto tk_xxxxzz_xxxxxy = pbuffer.data(idx_kin_ii + 141);

    auto tk_xxxxzz_xxxxxz = pbuffer.data(idx_kin_ii + 142);

    auto tk_xxxxzz_xxxxyy = pbuffer.data(idx_kin_ii + 143);

    auto tk_xxxxzz_xxxxyz = pbuffer.data(idx_kin_ii + 144);

    auto tk_xxxxzz_xxxxzz = pbuffer.data(idx_kin_ii + 145);

    auto tk_xxxxzz_xxxyyy = pbuffer.data(idx_kin_ii + 146);

    auto tk_xxxxzz_xxxyyz = pbuffer.data(idx_kin_ii + 147);

    auto tk_xxxxzz_xxxyzz = pbuffer.data(idx_kin_ii + 148);

    auto tk_xxxxzz_xxxzzz = pbuffer.data(idx_kin_ii + 149);

    auto tk_xxxxzz_xxyyyy = pbuffer.data(idx_kin_ii + 150);

    auto tk_xxxxzz_xxyyyz = pbuffer.data(idx_kin_ii + 151);

    auto tk_xxxxzz_xxyyzz = pbuffer.data(idx_kin_ii + 152);

    auto tk_xxxxzz_xxyzzz = pbuffer.data(idx_kin_ii + 153);

    auto tk_xxxxzz_xxzzzz = pbuffer.data(idx_kin_ii + 154);

    auto tk_xxxxzz_xyyyyy = pbuffer.data(idx_kin_ii + 155);

    auto tk_xxxxzz_xyyyyz = pbuffer.data(idx_kin_ii + 156);

    auto tk_xxxxzz_xyyyzz = pbuffer.data(idx_kin_ii + 157);

    auto tk_xxxxzz_xyyzzz = pbuffer.data(idx_kin_ii + 158);

    auto tk_xxxxzz_xyzzzz = pbuffer.data(idx_kin_ii + 159);

    auto tk_xxxxzz_xzzzzz = pbuffer.data(idx_kin_ii + 160);

    auto tk_xxxxzz_yyyyyy = pbuffer.data(idx_kin_ii + 161);

    auto tk_xxxxzz_yyyyyz = pbuffer.data(idx_kin_ii + 162);

    auto tk_xxxxzz_yyyyzz = pbuffer.data(idx_kin_ii + 163);

    auto tk_xxxxzz_yyyzzz = pbuffer.data(idx_kin_ii + 164);

    auto tk_xxxxzz_yyzzzz = pbuffer.data(idx_kin_ii + 165);

    auto tk_xxxxzz_yzzzzz = pbuffer.data(idx_kin_ii + 166);

    auto tk_xxxxzz_zzzzzz = pbuffer.data(idx_kin_ii + 167);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tk_xxxx_xxxxxx,   \
                             tk_xxxx_xxxxxy,   \
                             tk_xxxx_xxxxyy,   \
                             tk_xxxx_xxxyyy,   \
                             tk_xxxx_xxyyyy,   \
                             tk_xxxx_xyyyyy,   \
                             tk_xxxxz_xxxxxx,  \
                             tk_xxxxz_xxxxxy,  \
                             tk_xxxxz_xxxxyy,  \
                             tk_xxxxz_xxxyyy,  \
                             tk_xxxxz_xxyyyy,  \
                             tk_xxxxz_xyyyyy,  \
                             tk_xxxxzz_xxxxxx, \
                             tk_xxxxzz_xxxxxy, \
                             tk_xxxxzz_xxxxxz, \
                             tk_xxxxzz_xxxxyy, \
                             tk_xxxxzz_xxxxyz, \
                             tk_xxxxzz_xxxxzz, \
                             tk_xxxxzz_xxxyyy, \
                             tk_xxxxzz_xxxyyz, \
                             tk_xxxxzz_xxxyzz, \
                             tk_xxxxzz_xxxzzz, \
                             tk_xxxxzz_xxyyyy, \
                             tk_xxxxzz_xxyyyz, \
                             tk_xxxxzz_xxyyzz, \
                             tk_xxxxzz_xxyzzz, \
                             tk_xxxxzz_xxzzzz, \
                             tk_xxxxzz_xyyyyy, \
                             tk_xxxxzz_xyyyyz, \
                             tk_xxxxzz_xyyyzz, \
                             tk_xxxxzz_xyyzzz, \
                             tk_xxxxzz_xyzzzz, \
                             tk_xxxxzz_xzzzzz, \
                             tk_xxxxzz_yyyyyy, \
                             tk_xxxxzz_yyyyyz, \
                             tk_xxxxzz_yyyyzz, \
                             tk_xxxxzz_yyyzzz, \
                             tk_xxxxzz_yyzzzz, \
                             tk_xxxxzz_yzzzzz, \
                             tk_xxxxzz_zzzzzz, \
                             tk_xxxzz_xxxxxz,  \
                             tk_xxxzz_xxxxyz,  \
                             tk_xxxzz_xxxxz,   \
                             tk_xxxzz_xxxxzz,  \
                             tk_xxxzz_xxxyyz,  \
                             tk_xxxzz_xxxyz,   \
                             tk_xxxzz_xxxyzz,  \
                             tk_xxxzz_xxxzz,   \
                             tk_xxxzz_xxxzzz,  \
                             tk_xxxzz_xxyyyz,  \
                             tk_xxxzz_xxyyz,   \
                             tk_xxxzz_xxyyzz,  \
                             tk_xxxzz_xxyzz,   \
                             tk_xxxzz_xxyzzz,  \
                             tk_xxxzz_xxzzz,   \
                             tk_xxxzz_xxzzzz,  \
                             tk_xxxzz_xyyyyz,  \
                             tk_xxxzz_xyyyz,   \
                             tk_xxxzz_xyyyzz,  \
                             tk_xxxzz_xyyzz,   \
                             tk_xxxzz_xyyzzz,  \
                             tk_xxxzz_xyzzz,   \
                             tk_xxxzz_xyzzzz,  \
                             tk_xxxzz_xzzzz,   \
                             tk_xxxzz_xzzzzz,  \
                             tk_xxxzz_yyyyyy,  \
                             tk_xxxzz_yyyyyz,  \
                             tk_xxxzz_yyyyz,   \
                             tk_xxxzz_yyyyzz,  \
                             tk_xxxzz_yyyzz,   \
                             tk_xxxzz_yyyzzz,  \
                             tk_xxxzz_yyzzz,   \
                             tk_xxxzz_yyzzzz,  \
                             tk_xxxzz_yzzzz,   \
                             tk_xxxzz_yzzzzz,  \
                             tk_xxxzz_zzzzz,   \
                             tk_xxxzz_zzzzzz,  \
                             tk_xxzz_xxxxxz,   \
                             tk_xxzz_xxxxyz,   \
                             tk_xxzz_xxxxzz,   \
                             tk_xxzz_xxxyyz,   \
                             tk_xxzz_xxxyzz,   \
                             tk_xxzz_xxxzzz,   \
                             tk_xxzz_xxyyyz,   \
                             tk_xxzz_xxyyzz,   \
                             tk_xxzz_xxyzzz,   \
                             tk_xxzz_xxzzzz,   \
                             tk_xxzz_xyyyyz,   \
                             tk_xxzz_xyyyzz,   \
                             tk_xxzz_xyyzzz,   \
                             tk_xxzz_xyzzzz,   \
                             tk_xxzz_xzzzzz,   \
                             tk_xxzz_yyyyyy,   \
                             tk_xxzz_yyyyyz,   \
                             tk_xxzz_yyyyzz,   \
                             tk_xxzz_yyyzzz,   \
                             tk_xxzz_yyzzzz,   \
                             tk_xxzz_yzzzzz,   \
                             tk_xxzz_zzzzzz,   \
                             ts_xxxx_xxxxxx,   \
                             ts_xxxx_xxxxxy,   \
                             ts_xxxx_xxxxyy,   \
                             ts_xxxx_xxxyyy,   \
                             ts_xxxx_xxyyyy,   \
                             ts_xxxx_xyyyyy,   \
                             ts_xxxxzz_xxxxxx, \
                             ts_xxxxzz_xxxxxy, \
                             ts_xxxxzz_xxxxxz, \
                             ts_xxxxzz_xxxxyy, \
                             ts_xxxxzz_xxxxyz, \
                             ts_xxxxzz_xxxxzz, \
                             ts_xxxxzz_xxxyyy, \
                             ts_xxxxzz_xxxyyz, \
                             ts_xxxxzz_xxxyzz, \
                             ts_xxxxzz_xxxzzz, \
                             ts_xxxxzz_xxyyyy, \
                             ts_xxxxzz_xxyyyz, \
                             ts_xxxxzz_xxyyzz, \
                             ts_xxxxzz_xxyzzz, \
                             ts_xxxxzz_xxzzzz, \
                             ts_xxxxzz_xyyyyy, \
                             ts_xxxxzz_xyyyyz, \
                             ts_xxxxzz_xyyyzz, \
                             ts_xxxxzz_xyyzzz, \
                             ts_xxxxzz_xyzzzz, \
                             ts_xxxxzz_xzzzzz, \
                             ts_xxxxzz_yyyyyy, \
                             ts_xxxxzz_yyyyyz, \
                             ts_xxxxzz_yyyyzz, \
                             ts_xxxxzz_yyyzzz, \
                             ts_xxxxzz_yyzzzz, \
                             ts_xxxxzz_yzzzzz, \
                             ts_xxxxzz_zzzzzz, \
                             ts_xxzz_xxxxxz,   \
                             ts_xxzz_xxxxyz,   \
                             ts_xxzz_xxxxzz,   \
                             ts_xxzz_xxxyyz,   \
                             ts_xxzz_xxxyzz,   \
                             ts_xxzz_xxxzzz,   \
                             ts_xxzz_xxyyyz,   \
                             ts_xxzz_xxyyzz,   \
                             ts_xxzz_xxyzzz,   \
                             ts_xxzz_xxzzzz,   \
                             ts_xxzz_xyyyyz,   \
                             ts_xxzz_xyyyzz,   \
                             ts_xxzz_xyyzzz,   \
                             ts_xxzz_xyzzzz,   \
                             ts_xxzz_xzzzzz,   \
                             ts_xxzz_yyyyyy,   \
                             ts_xxzz_yyyyyz,   \
                             ts_xxzz_yyyyzz,   \
                             ts_xxzz_yyyzzz,   \
                             ts_xxzz_yyzzzz,   \
                             ts_xxzz_yzzzzz,   \
                             ts_xxzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxzz_xxxxxx[i] =
            -2.0 * ts_xxxx_xxxxxx[i] * fbe_0 * fz_0 + tk_xxxx_xxxxxx[i] * fe_0 + tk_xxxxz_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxxxx[i] * fz_0;

        tk_xxxxzz_xxxxxy[i] =
            -2.0 * ts_xxxx_xxxxxy[i] * fbe_0 * fz_0 + tk_xxxx_xxxxxy[i] * fe_0 + tk_xxxxz_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxxxy[i] * fz_0;

        tk_xxxxzz_xxxxxz[i] = -6.0 * ts_xxzz_xxxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxxz[i] * fe_0 + 5.0 * tk_xxxzz_xxxxz[i] * fe_0 +
                              tk_xxxzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxxxz[i] * fz_0;

        tk_xxxxzz_xxxxyy[i] =
            -2.0 * ts_xxxx_xxxxyy[i] * fbe_0 * fz_0 + tk_xxxx_xxxxyy[i] * fe_0 + tk_xxxxz_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxxyy[i] * fz_0;

        tk_xxxxzz_xxxxyz[i] = -6.0 * ts_xxzz_xxxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxyz[i] * fe_0 + 4.0 * tk_xxxzz_xxxyz[i] * fe_0 +
                              tk_xxxzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxxyz[i] * fz_0;

        tk_xxxxzz_xxxxzz[i] = -6.0 * ts_xxzz_xxxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxzz[i] * fe_0 + 4.0 * tk_xxxzz_xxxzz[i] * fe_0 +
                              tk_xxxzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxxzz[i] * fz_0;

        tk_xxxxzz_xxxyyy[i] =
            -2.0 * ts_xxxx_xxxyyy[i] * fbe_0 * fz_0 + tk_xxxx_xxxyyy[i] * fe_0 + tk_xxxxz_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxyyy[i] * fz_0;

        tk_xxxxzz_xxxyyz[i] = -6.0 * ts_xxzz_xxxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxyyz[i] * fe_0 + 3.0 * tk_xxxzz_xxyyz[i] * fe_0 +
                              tk_xxxzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxyyz[i] * fz_0;

        tk_xxxxzz_xxxyzz[i] = -6.0 * ts_xxzz_xxxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxyzz[i] * fe_0 + 3.0 * tk_xxxzz_xxyzz[i] * fe_0 +
                              tk_xxxzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxyzz[i] * fz_0;

        tk_xxxxzz_xxxzzz[i] = -6.0 * ts_xxzz_xxxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxzzz[i] * fe_0 + 3.0 * tk_xxxzz_xxzzz[i] * fe_0 +
                              tk_xxxzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxzzz[i] * fz_0;

        tk_xxxxzz_xxyyyy[i] =
            -2.0 * ts_xxxx_xxyyyy[i] * fbe_0 * fz_0 + tk_xxxx_xxyyyy[i] * fe_0 + tk_xxxxz_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxyyyy[i] * fz_0;

        tk_xxxxzz_xxyyyz[i] = -6.0 * ts_xxzz_xxyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyyyz[i] * fe_0 + 2.0 * tk_xxxzz_xyyyz[i] * fe_0 +
                              tk_xxxzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxyyyz[i] * fz_0;

        tk_xxxxzz_xxyyzz[i] = -6.0 * ts_xxzz_xxyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyyzz[i] * fe_0 + 2.0 * tk_xxxzz_xyyzz[i] * fe_0 +
                              tk_xxxzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxyyzz[i] * fz_0;

        tk_xxxxzz_xxyzzz[i] = -6.0 * ts_xxzz_xxyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyzzz[i] * fe_0 + 2.0 * tk_xxxzz_xyzzz[i] * fe_0 +
                              tk_xxxzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxyzzz[i] * fz_0;

        tk_xxxxzz_xxzzzz[i] = -6.0 * ts_xxzz_xxzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxzzzz[i] * fe_0 + 2.0 * tk_xxxzz_xzzzz[i] * fe_0 +
                              tk_xxxzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxzzzz[i] * fz_0;

        tk_xxxxzz_xyyyyy[i] =
            -2.0 * ts_xxxx_xyyyyy[i] * fbe_0 * fz_0 + tk_xxxx_xyyyyy[i] * fe_0 + tk_xxxxz_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xyyyyy[i] * fz_0;

        tk_xxxxzz_xyyyyz[i] = -6.0 * ts_xxzz_xyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyyyz[i] * fe_0 + tk_xxxzz_yyyyz[i] * fe_0 +
                              tk_xxxzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyyyyz[i] * fz_0;

        tk_xxxxzz_xyyyzz[i] = -6.0 * ts_xxzz_xyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyyzz[i] * fe_0 + tk_xxxzz_yyyzz[i] * fe_0 +
                              tk_xxxzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyyyzz[i] * fz_0;

        tk_xxxxzz_xyyzzz[i] = -6.0 * ts_xxzz_xyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyzzz[i] * fe_0 + tk_xxxzz_yyzzz[i] * fe_0 +
                              tk_xxxzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyyzzz[i] * fz_0;

        tk_xxxxzz_xyzzzz[i] = -6.0 * ts_xxzz_xyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyzzzz[i] * fe_0 + tk_xxxzz_yzzzz[i] * fe_0 +
                              tk_xxxzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyzzzz[i] * fz_0;

        tk_xxxxzz_xzzzzz[i] = -6.0 * ts_xxzz_xzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xzzzzz[i] * fe_0 + tk_xxxzz_zzzzz[i] * fe_0 +
                              tk_xxxzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xzzzzz[i] * fz_0;

        tk_xxxxzz_yyyyyy[i] = -6.0 * ts_xxzz_yyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyyyy[i] * fe_0 + tk_xxxzz_yyyyyy[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_yyyyyy[i] * fz_0;

        tk_xxxxzz_yyyyyz[i] = -6.0 * ts_xxzz_yyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyyyz[i] * fe_0 + tk_xxxzz_yyyyyz[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_yyyyyz[i] * fz_0;

        tk_xxxxzz_yyyyzz[i] = -6.0 * ts_xxzz_yyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyyzz[i] * fe_0 + tk_xxxzz_yyyyzz[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_yyyyzz[i] * fz_0;

        tk_xxxxzz_yyyzzz[i] = -6.0 * ts_xxzz_yyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyzzz[i] * fe_0 + tk_xxxzz_yyyzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_yyyzzz[i] * fz_0;

        tk_xxxxzz_yyzzzz[i] = -6.0 * ts_xxzz_yyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyzzzz[i] * fe_0 + tk_xxxzz_yyzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_yyzzzz[i] * fz_0;

        tk_xxxxzz_yzzzzz[i] = -6.0 * ts_xxzz_yzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yzzzzz[i] * fe_0 + tk_xxxzz_yzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_yzzzzz[i] * fz_0;

        tk_xxxxzz_zzzzzz[i] = -6.0 * ts_xxzz_zzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_zzzzzz[i] * fe_0 + tk_xxxzz_zzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxxzz_zzzzzz[i] * fz_0;
    }

    // Set up 168-196 components of targeted buffer : II

    auto tk_xxxyyy_xxxxxx = pbuffer.data(idx_kin_ii + 168);

    auto tk_xxxyyy_xxxxxy = pbuffer.data(idx_kin_ii + 169);

    auto tk_xxxyyy_xxxxxz = pbuffer.data(idx_kin_ii + 170);

    auto tk_xxxyyy_xxxxyy = pbuffer.data(idx_kin_ii + 171);

    auto tk_xxxyyy_xxxxyz = pbuffer.data(idx_kin_ii + 172);

    auto tk_xxxyyy_xxxxzz = pbuffer.data(idx_kin_ii + 173);

    auto tk_xxxyyy_xxxyyy = pbuffer.data(idx_kin_ii + 174);

    auto tk_xxxyyy_xxxyyz = pbuffer.data(idx_kin_ii + 175);

    auto tk_xxxyyy_xxxyzz = pbuffer.data(idx_kin_ii + 176);

    auto tk_xxxyyy_xxxzzz = pbuffer.data(idx_kin_ii + 177);

    auto tk_xxxyyy_xxyyyy = pbuffer.data(idx_kin_ii + 178);

    auto tk_xxxyyy_xxyyyz = pbuffer.data(idx_kin_ii + 179);

    auto tk_xxxyyy_xxyyzz = pbuffer.data(idx_kin_ii + 180);

    auto tk_xxxyyy_xxyzzz = pbuffer.data(idx_kin_ii + 181);

    auto tk_xxxyyy_xxzzzz = pbuffer.data(idx_kin_ii + 182);

    auto tk_xxxyyy_xyyyyy = pbuffer.data(idx_kin_ii + 183);

    auto tk_xxxyyy_xyyyyz = pbuffer.data(idx_kin_ii + 184);

    auto tk_xxxyyy_xyyyzz = pbuffer.data(idx_kin_ii + 185);

    auto tk_xxxyyy_xyyzzz = pbuffer.data(idx_kin_ii + 186);

    auto tk_xxxyyy_xyzzzz = pbuffer.data(idx_kin_ii + 187);

    auto tk_xxxyyy_xzzzzz = pbuffer.data(idx_kin_ii + 188);

    auto tk_xxxyyy_yyyyyy = pbuffer.data(idx_kin_ii + 189);

    auto tk_xxxyyy_yyyyyz = pbuffer.data(idx_kin_ii + 190);

    auto tk_xxxyyy_yyyyzz = pbuffer.data(idx_kin_ii + 191);

    auto tk_xxxyyy_yyyzzz = pbuffer.data(idx_kin_ii + 192);

    auto tk_xxxyyy_yyzzzz = pbuffer.data(idx_kin_ii + 193);

    auto tk_xxxyyy_yzzzzz = pbuffer.data(idx_kin_ii + 194);

    auto tk_xxxyyy_zzzzzz = pbuffer.data(idx_kin_ii + 195);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tk_xxxy_xxxxxx,   \
                             tk_xxxy_xxxxxz,   \
                             tk_xxxy_xxxxzz,   \
                             tk_xxxy_xxxzzz,   \
                             tk_xxxy_xxzzzz,   \
                             tk_xxxy_xzzzzz,   \
                             tk_xxxyy_xxxxxx,  \
                             tk_xxxyy_xxxxxz,  \
                             tk_xxxyy_xxxxzz,  \
                             tk_xxxyy_xxxzzz,  \
                             tk_xxxyy_xxzzzz,  \
                             tk_xxxyy_xzzzzz,  \
                             tk_xxxyyy_xxxxxx, \
                             tk_xxxyyy_xxxxxy, \
                             tk_xxxyyy_xxxxxz, \
                             tk_xxxyyy_xxxxyy, \
                             tk_xxxyyy_xxxxyz, \
                             tk_xxxyyy_xxxxzz, \
                             tk_xxxyyy_xxxyyy, \
                             tk_xxxyyy_xxxyyz, \
                             tk_xxxyyy_xxxyzz, \
                             tk_xxxyyy_xxxzzz, \
                             tk_xxxyyy_xxyyyy, \
                             tk_xxxyyy_xxyyyz, \
                             tk_xxxyyy_xxyyzz, \
                             tk_xxxyyy_xxyzzz, \
                             tk_xxxyyy_xxzzzz, \
                             tk_xxxyyy_xyyyyy, \
                             tk_xxxyyy_xyyyyz, \
                             tk_xxxyyy_xyyyzz, \
                             tk_xxxyyy_xyyzzz, \
                             tk_xxxyyy_xyzzzz, \
                             tk_xxxyyy_xzzzzz, \
                             tk_xxxyyy_yyyyyy, \
                             tk_xxxyyy_yyyyyz, \
                             tk_xxxyyy_yyyyzz, \
                             tk_xxxyyy_yyyzzz, \
                             tk_xxxyyy_yyzzzz, \
                             tk_xxxyyy_yzzzzz, \
                             tk_xxxyyy_zzzzzz, \
                             tk_xxyyy_xxxxxy,  \
                             tk_xxyyy_xxxxy,   \
                             tk_xxyyy_xxxxyy,  \
                             tk_xxyyy_xxxxyz,  \
                             tk_xxyyy_xxxyy,   \
                             tk_xxyyy_xxxyyy,  \
                             tk_xxyyy_xxxyyz,  \
                             tk_xxyyy_xxxyz,   \
                             tk_xxyyy_xxxyzz,  \
                             tk_xxyyy_xxyyy,   \
                             tk_xxyyy_xxyyyy,  \
                             tk_xxyyy_xxyyyz,  \
                             tk_xxyyy_xxyyz,   \
                             tk_xxyyy_xxyyzz,  \
                             tk_xxyyy_xxyzz,   \
                             tk_xxyyy_xxyzzz,  \
                             tk_xxyyy_xyyyy,   \
                             tk_xxyyy_xyyyyy,  \
                             tk_xxyyy_xyyyyz,  \
                             tk_xxyyy_xyyyz,   \
                             tk_xxyyy_xyyyzz,  \
                             tk_xxyyy_xyyzz,   \
                             tk_xxyyy_xyyzzz,  \
                             tk_xxyyy_xyzzz,   \
                             tk_xxyyy_xyzzzz,  \
                             tk_xxyyy_yyyyy,   \
                             tk_xxyyy_yyyyyy,  \
                             tk_xxyyy_yyyyyz,  \
                             tk_xxyyy_yyyyz,   \
                             tk_xxyyy_yyyyzz,  \
                             tk_xxyyy_yyyzz,   \
                             tk_xxyyy_yyyzzz,  \
                             tk_xxyyy_yyzzz,   \
                             tk_xxyyy_yyzzzz,  \
                             tk_xxyyy_yzzzz,   \
                             tk_xxyyy_yzzzzz,  \
                             tk_xxyyy_zzzzzz,  \
                             tk_xyyy_xxxxxy,   \
                             tk_xyyy_xxxxyy,   \
                             tk_xyyy_xxxxyz,   \
                             tk_xyyy_xxxyyy,   \
                             tk_xyyy_xxxyyz,   \
                             tk_xyyy_xxxyzz,   \
                             tk_xyyy_xxyyyy,   \
                             tk_xyyy_xxyyyz,   \
                             tk_xyyy_xxyyzz,   \
                             tk_xyyy_xxyzzz,   \
                             tk_xyyy_xyyyyy,   \
                             tk_xyyy_xyyyyz,   \
                             tk_xyyy_xyyyzz,   \
                             tk_xyyy_xyyzzz,   \
                             tk_xyyy_xyzzzz,   \
                             tk_xyyy_yyyyyy,   \
                             tk_xyyy_yyyyyz,   \
                             tk_xyyy_yyyyzz,   \
                             tk_xyyy_yyyzzz,   \
                             tk_xyyy_yyzzzz,   \
                             tk_xyyy_yzzzzz,   \
                             tk_xyyy_zzzzzz,   \
                             ts_xxxy_xxxxxx,   \
                             ts_xxxy_xxxxxz,   \
                             ts_xxxy_xxxxzz,   \
                             ts_xxxy_xxxzzz,   \
                             ts_xxxy_xxzzzz,   \
                             ts_xxxy_xzzzzz,   \
                             ts_xxxyyy_xxxxxx, \
                             ts_xxxyyy_xxxxxy, \
                             ts_xxxyyy_xxxxxz, \
                             ts_xxxyyy_xxxxyy, \
                             ts_xxxyyy_xxxxyz, \
                             ts_xxxyyy_xxxxzz, \
                             ts_xxxyyy_xxxyyy, \
                             ts_xxxyyy_xxxyyz, \
                             ts_xxxyyy_xxxyzz, \
                             ts_xxxyyy_xxxzzz, \
                             ts_xxxyyy_xxyyyy, \
                             ts_xxxyyy_xxyyyz, \
                             ts_xxxyyy_xxyyzz, \
                             ts_xxxyyy_xxyzzz, \
                             ts_xxxyyy_xxzzzz, \
                             ts_xxxyyy_xyyyyy, \
                             ts_xxxyyy_xyyyyz, \
                             ts_xxxyyy_xyyyzz, \
                             ts_xxxyyy_xyyzzz, \
                             ts_xxxyyy_xyzzzz, \
                             ts_xxxyyy_xzzzzz, \
                             ts_xxxyyy_yyyyyy, \
                             ts_xxxyyy_yyyyyz, \
                             ts_xxxyyy_yyyyzz, \
                             ts_xxxyyy_yyyzzz, \
                             ts_xxxyyy_yyzzzz, \
                             ts_xxxyyy_yzzzzz, \
                             ts_xxxyyy_zzzzzz, \
                             ts_xyyy_xxxxxy,   \
                             ts_xyyy_xxxxyy,   \
                             ts_xyyy_xxxxyz,   \
                             ts_xyyy_xxxyyy,   \
                             ts_xyyy_xxxyyz,   \
                             ts_xyyy_xxxyzz,   \
                             ts_xyyy_xxyyyy,   \
                             ts_xyyy_xxyyyz,   \
                             ts_xyyy_xxyyzz,   \
                             ts_xyyy_xxyzzz,   \
                             ts_xyyy_xyyyyy,   \
                             ts_xyyy_xyyyyz,   \
                             ts_xyyy_xyyyzz,   \
                             ts_xyyy_xyyzzz,   \
                             ts_xyyy_xyzzzz,   \
                             ts_xyyy_yyyyyy,   \
                             ts_xyyy_yyyyyz,   \
                             ts_xyyy_yyyyzz,   \
                             ts_xyyy_yyyzzz,   \
                             ts_xyyy_yyzzzz,   \
                             ts_xyyy_yzzzzz,   \
                             ts_xyyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyyy_xxxxxx[i] = -4.0 * ts_xxxy_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxxxx[i] * fe_0 + tk_xxxyy_xxxxxx[i] * pa_y[i] +
                              2.0 * ts_xxxyyy_xxxxxx[i] * fz_0;

        tk_xxxyyy_xxxxxy[i] = -4.0 * ts_xyyy_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxxxy[i] * fe_0 + 5.0 * tk_xxyyy_xxxxy[i] * fe_0 +
                              tk_xxyyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxxxy[i] * fz_0;

        tk_xxxyyy_xxxxxz[i] = -4.0 * ts_xxxy_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxxxz[i] * fe_0 + tk_xxxyy_xxxxxz[i] * pa_y[i] +
                              2.0 * ts_xxxyyy_xxxxxz[i] * fz_0;

        tk_xxxyyy_xxxxyy[i] = -4.0 * ts_xyyy_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxxyy[i] * fe_0 + 4.0 * tk_xxyyy_xxxyy[i] * fe_0 +
                              tk_xxyyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxxyy[i] * fz_0;

        tk_xxxyyy_xxxxyz[i] = -4.0 * ts_xyyy_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxxyz[i] * fe_0 + 4.0 * tk_xxyyy_xxxyz[i] * fe_0 +
                              tk_xxyyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxxyz[i] * fz_0;

        tk_xxxyyy_xxxxzz[i] = -4.0 * ts_xxxy_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxxzz[i] * fe_0 + tk_xxxyy_xxxxzz[i] * pa_y[i] +
                              2.0 * ts_xxxyyy_xxxxzz[i] * fz_0;

        tk_xxxyyy_xxxyyy[i] = -4.0 * ts_xyyy_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxyyy[i] * fe_0 + 3.0 * tk_xxyyy_xxyyy[i] * fe_0 +
                              tk_xxyyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxyyy[i] * fz_0;

        tk_xxxyyy_xxxyyz[i] = -4.0 * ts_xyyy_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxyyz[i] * fe_0 + 3.0 * tk_xxyyy_xxyyz[i] * fe_0 +
                              tk_xxyyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxyyz[i] * fz_0;

        tk_xxxyyy_xxxyzz[i] = -4.0 * ts_xyyy_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxyzz[i] * fe_0 + 3.0 * tk_xxyyy_xxyzz[i] * fe_0 +
                              tk_xxyyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxyzz[i] * fz_0;

        tk_xxxyyy_xxxzzz[i] = -4.0 * ts_xxxy_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxzzz[i] * fe_0 + tk_xxxyy_xxxzzz[i] * pa_y[i] +
                              2.0 * ts_xxxyyy_xxxzzz[i] * fz_0;

        tk_xxxyyy_xxyyyy[i] = -4.0 * ts_xyyy_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyyyy[i] * fe_0 + 2.0 * tk_xxyyy_xyyyy[i] * fe_0 +
                              tk_xxyyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyyyy[i] * fz_0;

        tk_xxxyyy_xxyyyz[i] = -4.0 * ts_xyyy_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyyyz[i] * fe_0 + 2.0 * tk_xxyyy_xyyyz[i] * fe_0 +
                              tk_xxyyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyyyz[i] * fz_0;

        tk_xxxyyy_xxyyzz[i] = -4.0 * ts_xyyy_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyyzz[i] * fe_0 + 2.0 * tk_xxyyy_xyyzz[i] * fe_0 +
                              tk_xxyyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyyzz[i] * fz_0;

        tk_xxxyyy_xxyzzz[i] = -4.0 * ts_xyyy_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyzzz[i] * fe_0 + 2.0 * tk_xxyyy_xyzzz[i] * fe_0 +
                              tk_xxyyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyzzz[i] * fz_0;

        tk_xxxyyy_xxzzzz[i] = -4.0 * ts_xxxy_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxzzzz[i] * fe_0 + tk_xxxyy_xxzzzz[i] * pa_y[i] +
                              2.0 * ts_xxxyyy_xxzzzz[i] * fz_0;

        tk_xxxyyy_xyyyyy[i] = -4.0 * ts_xyyy_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyyyy[i] * fe_0 + tk_xxyyy_yyyyy[i] * fe_0 +
                              tk_xxyyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyyyy[i] * fz_0;

        tk_xxxyyy_xyyyyz[i] = -4.0 * ts_xyyy_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyyyz[i] * fe_0 + tk_xxyyy_yyyyz[i] * fe_0 +
                              tk_xxyyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyyyz[i] * fz_0;

        tk_xxxyyy_xyyyzz[i] = -4.0 * ts_xyyy_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyyzz[i] * fe_0 + tk_xxyyy_yyyzz[i] * fe_0 +
                              tk_xxyyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyyzz[i] * fz_0;

        tk_xxxyyy_xyyzzz[i] = -4.0 * ts_xyyy_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyzzz[i] * fe_0 + tk_xxyyy_yyzzz[i] * fe_0 +
                              tk_xxyyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyzzz[i] * fz_0;

        tk_xxxyyy_xyzzzz[i] = -4.0 * ts_xyyy_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyzzzz[i] * fe_0 + tk_xxyyy_yzzzz[i] * fe_0 +
                              tk_xxyyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyzzzz[i] * fz_0;

        tk_xxxyyy_xzzzzz[i] = -4.0 * ts_xxxy_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xzzzzz[i] * fe_0 + tk_xxxyy_xzzzzz[i] * pa_y[i] +
                              2.0 * ts_xxxyyy_xzzzzz[i] * fz_0;

        tk_xxxyyy_yyyyyy[i] = -4.0 * ts_xyyy_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyyyy[i] * fe_0 + tk_xxyyy_yyyyyy[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_yyyyyy[i] * fz_0;

        tk_xxxyyy_yyyyyz[i] = -4.0 * ts_xyyy_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyyyz[i] * fe_0 + tk_xxyyy_yyyyyz[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_yyyyyz[i] * fz_0;

        tk_xxxyyy_yyyyzz[i] = -4.0 * ts_xyyy_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyyzz[i] * fe_0 + tk_xxyyy_yyyyzz[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_yyyyzz[i] * fz_0;

        tk_xxxyyy_yyyzzz[i] = -4.0 * ts_xyyy_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyzzz[i] * fe_0 + tk_xxyyy_yyyzzz[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_yyyzzz[i] * fz_0;

        tk_xxxyyy_yyzzzz[i] = -4.0 * ts_xyyy_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyzzzz[i] * fe_0 + tk_xxyyy_yyzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_yyzzzz[i] * fz_0;

        tk_xxxyyy_yzzzzz[i] = -4.0 * ts_xyyy_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yzzzzz[i] * fe_0 + tk_xxyyy_yzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_yzzzzz[i] * fz_0;

        tk_xxxyyy_zzzzzz[i] = -4.0 * ts_xyyy_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_zzzzzz[i] * fe_0 + tk_xxyyy_zzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxyyy_zzzzzz[i] * fz_0;
    }

    // Set up 196-224 components of targeted buffer : II

    auto tk_xxxyyz_xxxxxx = pbuffer.data(idx_kin_ii + 196);

    auto tk_xxxyyz_xxxxxy = pbuffer.data(idx_kin_ii + 197);

    auto tk_xxxyyz_xxxxxz = pbuffer.data(idx_kin_ii + 198);

    auto tk_xxxyyz_xxxxyy = pbuffer.data(idx_kin_ii + 199);

    auto tk_xxxyyz_xxxxyz = pbuffer.data(idx_kin_ii + 200);

    auto tk_xxxyyz_xxxxzz = pbuffer.data(idx_kin_ii + 201);

    auto tk_xxxyyz_xxxyyy = pbuffer.data(idx_kin_ii + 202);

    auto tk_xxxyyz_xxxyyz = pbuffer.data(idx_kin_ii + 203);

    auto tk_xxxyyz_xxxyzz = pbuffer.data(idx_kin_ii + 204);

    auto tk_xxxyyz_xxxzzz = pbuffer.data(idx_kin_ii + 205);

    auto tk_xxxyyz_xxyyyy = pbuffer.data(idx_kin_ii + 206);

    auto tk_xxxyyz_xxyyyz = pbuffer.data(idx_kin_ii + 207);

    auto tk_xxxyyz_xxyyzz = pbuffer.data(idx_kin_ii + 208);

    auto tk_xxxyyz_xxyzzz = pbuffer.data(idx_kin_ii + 209);

    auto tk_xxxyyz_xxzzzz = pbuffer.data(idx_kin_ii + 210);

    auto tk_xxxyyz_xyyyyy = pbuffer.data(idx_kin_ii + 211);

    auto tk_xxxyyz_xyyyyz = pbuffer.data(idx_kin_ii + 212);

    auto tk_xxxyyz_xyyyzz = pbuffer.data(idx_kin_ii + 213);

    auto tk_xxxyyz_xyyzzz = pbuffer.data(idx_kin_ii + 214);

    auto tk_xxxyyz_xyzzzz = pbuffer.data(idx_kin_ii + 215);

    auto tk_xxxyyz_xzzzzz = pbuffer.data(idx_kin_ii + 216);

    auto tk_xxxyyz_yyyyyy = pbuffer.data(idx_kin_ii + 217);

    auto tk_xxxyyz_yyyyyz = pbuffer.data(idx_kin_ii + 218);

    auto tk_xxxyyz_yyyyzz = pbuffer.data(idx_kin_ii + 219);

    auto tk_xxxyyz_yyyzzz = pbuffer.data(idx_kin_ii + 220);

    auto tk_xxxyyz_yyzzzz = pbuffer.data(idx_kin_ii + 221);

    auto tk_xxxyyz_yzzzzz = pbuffer.data(idx_kin_ii + 222);

    auto tk_xxxyyz_zzzzzz = pbuffer.data(idx_kin_ii + 223);

#pragma omp simd aligned(pa_z,                 \
                             tk_xxxyy_xxxxx,   \
                             tk_xxxyy_xxxxxx,  \
                             tk_xxxyy_xxxxxy,  \
                             tk_xxxyy_xxxxxz,  \
                             tk_xxxyy_xxxxy,   \
                             tk_xxxyy_xxxxyy,  \
                             tk_xxxyy_xxxxyz,  \
                             tk_xxxyy_xxxxz,   \
                             tk_xxxyy_xxxxzz,  \
                             tk_xxxyy_xxxyy,   \
                             tk_xxxyy_xxxyyy,  \
                             tk_xxxyy_xxxyyz,  \
                             tk_xxxyy_xxxyz,   \
                             tk_xxxyy_xxxyzz,  \
                             tk_xxxyy_xxxzz,   \
                             tk_xxxyy_xxxzzz,  \
                             tk_xxxyy_xxyyy,   \
                             tk_xxxyy_xxyyyy,  \
                             tk_xxxyy_xxyyyz,  \
                             tk_xxxyy_xxyyz,   \
                             tk_xxxyy_xxyyzz,  \
                             tk_xxxyy_xxyzz,   \
                             tk_xxxyy_xxyzzz,  \
                             tk_xxxyy_xxzzz,   \
                             tk_xxxyy_xxzzzz,  \
                             tk_xxxyy_xyyyy,   \
                             tk_xxxyy_xyyyyy,  \
                             tk_xxxyy_xyyyyz,  \
                             tk_xxxyy_xyyyz,   \
                             tk_xxxyy_xyyyzz,  \
                             tk_xxxyy_xyyzz,   \
                             tk_xxxyy_xyyzzz,  \
                             tk_xxxyy_xyzzz,   \
                             tk_xxxyy_xyzzzz,  \
                             tk_xxxyy_xzzzz,   \
                             tk_xxxyy_xzzzzz,  \
                             tk_xxxyy_yyyyy,   \
                             tk_xxxyy_yyyyyy,  \
                             tk_xxxyy_yyyyyz,  \
                             tk_xxxyy_yyyyz,   \
                             tk_xxxyy_yyyyzz,  \
                             tk_xxxyy_yyyzz,   \
                             tk_xxxyy_yyyzzz,  \
                             tk_xxxyy_yyzzz,   \
                             tk_xxxyy_yyzzzz,  \
                             tk_xxxyy_yzzzz,   \
                             tk_xxxyy_yzzzzz,  \
                             tk_xxxyy_zzzzz,   \
                             tk_xxxyy_zzzzzz,  \
                             tk_xxxyyz_xxxxxx, \
                             tk_xxxyyz_xxxxxy, \
                             tk_xxxyyz_xxxxxz, \
                             tk_xxxyyz_xxxxyy, \
                             tk_xxxyyz_xxxxyz, \
                             tk_xxxyyz_xxxxzz, \
                             tk_xxxyyz_xxxyyy, \
                             tk_xxxyyz_xxxyyz, \
                             tk_xxxyyz_xxxyzz, \
                             tk_xxxyyz_xxxzzz, \
                             tk_xxxyyz_xxyyyy, \
                             tk_xxxyyz_xxyyyz, \
                             tk_xxxyyz_xxyyzz, \
                             tk_xxxyyz_xxyzzz, \
                             tk_xxxyyz_xxzzzz, \
                             tk_xxxyyz_xyyyyy, \
                             tk_xxxyyz_xyyyyz, \
                             tk_xxxyyz_xyyyzz, \
                             tk_xxxyyz_xyyzzz, \
                             tk_xxxyyz_xyzzzz, \
                             tk_xxxyyz_xzzzzz, \
                             tk_xxxyyz_yyyyyy, \
                             tk_xxxyyz_yyyyyz, \
                             tk_xxxyyz_yyyyzz, \
                             tk_xxxyyz_yyyzzz, \
                             tk_xxxyyz_yyzzzz, \
                             tk_xxxyyz_yzzzzz, \
                             tk_xxxyyz_zzzzzz, \
                             ts_xxxyyz_xxxxxx, \
                             ts_xxxyyz_xxxxxy, \
                             ts_xxxyyz_xxxxxz, \
                             ts_xxxyyz_xxxxyy, \
                             ts_xxxyyz_xxxxyz, \
                             ts_xxxyyz_xxxxzz, \
                             ts_xxxyyz_xxxyyy, \
                             ts_xxxyyz_xxxyyz, \
                             ts_xxxyyz_xxxyzz, \
                             ts_xxxyyz_xxxzzz, \
                             ts_xxxyyz_xxyyyy, \
                             ts_xxxyyz_xxyyyz, \
                             ts_xxxyyz_xxyyzz, \
                             ts_xxxyyz_xxyzzz, \
                             ts_xxxyyz_xxzzzz, \
                             ts_xxxyyz_xyyyyy, \
                             ts_xxxyyz_xyyyyz, \
                             ts_xxxyyz_xyyyzz, \
                             ts_xxxyyz_xyyzzz, \
                             ts_xxxyyz_xyzzzz, \
                             ts_xxxyyz_xzzzzz, \
                             ts_xxxyyz_yyyyyy, \
                             ts_xxxyyz_yyyyyz, \
                             ts_xxxyyz_yyyyzz, \
                             ts_xxxyyz_yyyzzz, \
                             ts_xxxyyz_yyzzzz, \
                             ts_xxxyyz_yzzzzz, \
                             ts_xxxyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyyz_xxxxxx[i] = tk_xxxyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxxx[i] * fz_0;

        tk_xxxyyz_xxxxxy[i] = tk_xxxyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxxy[i] * fz_0;

        tk_xxxyyz_xxxxxz[i] = tk_xxxyy_xxxxx[i] * fe_0 + tk_xxxyy_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxxz[i] * fz_0;

        tk_xxxyyz_xxxxyy[i] = tk_xxxyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxyy[i] * fz_0;

        tk_xxxyyz_xxxxyz[i] = tk_xxxyy_xxxxy[i] * fe_0 + tk_xxxyy_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxyz[i] * fz_0;

        tk_xxxyyz_xxxxzz[i] = 2.0 * tk_xxxyy_xxxxz[i] * fe_0 + tk_xxxyy_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxzz[i] * fz_0;

        tk_xxxyyz_xxxyyy[i] = tk_xxxyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxyyy[i] * fz_0;

        tk_xxxyyz_xxxyyz[i] = tk_xxxyy_xxxyy[i] * fe_0 + tk_xxxyy_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxyyz[i] * fz_0;

        tk_xxxyyz_xxxyzz[i] = 2.0 * tk_xxxyy_xxxyz[i] * fe_0 + tk_xxxyy_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxyzz[i] * fz_0;

        tk_xxxyyz_xxxzzz[i] = 3.0 * tk_xxxyy_xxxzz[i] * fe_0 + tk_xxxyy_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxzzz[i] * fz_0;

        tk_xxxyyz_xxyyyy[i] = tk_xxxyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyyyy[i] * fz_0;

        tk_xxxyyz_xxyyyz[i] = tk_xxxyy_xxyyy[i] * fe_0 + tk_xxxyy_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyyyz[i] * fz_0;

        tk_xxxyyz_xxyyzz[i] = 2.0 * tk_xxxyy_xxyyz[i] * fe_0 + tk_xxxyy_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyyzz[i] * fz_0;

        tk_xxxyyz_xxyzzz[i] = 3.0 * tk_xxxyy_xxyzz[i] * fe_0 + tk_xxxyy_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyzzz[i] * fz_0;

        tk_xxxyyz_xxzzzz[i] = 4.0 * tk_xxxyy_xxzzz[i] * fe_0 + tk_xxxyy_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxzzzz[i] * fz_0;

        tk_xxxyyz_xyyyyy[i] = tk_xxxyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyyyy[i] * fz_0;

        tk_xxxyyz_xyyyyz[i] = tk_xxxyy_xyyyy[i] * fe_0 + tk_xxxyy_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyyyz[i] * fz_0;

        tk_xxxyyz_xyyyzz[i] = 2.0 * tk_xxxyy_xyyyz[i] * fe_0 + tk_xxxyy_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyyzz[i] * fz_0;

        tk_xxxyyz_xyyzzz[i] = 3.0 * tk_xxxyy_xyyzz[i] * fe_0 + tk_xxxyy_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyzzz[i] * fz_0;

        tk_xxxyyz_xyzzzz[i] = 4.0 * tk_xxxyy_xyzzz[i] * fe_0 + tk_xxxyy_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyzzzz[i] * fz_0;

        tk_xxxyyz_xzzzzz[i] = 5.0 * tk_xxxyy_xzzzz[i] * fe_0 + tk_xxxyy_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xzzzzz[i] * fz_0;

        tk_xxxyyz_yyyyyy[i] = tk_xxxyy_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyyyy[i] * fz_0;

        tk_xxxyyz_yyyyyz[i] = tk_xxxyy_yyyyy[i] * fe_0 + tk_xxxyy_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyyyz[i] * fz_0;

        tk_xxxyyz_yyyyzz[i] = 2.0 * tk_xxxyy_yyyyz[i] * fe_0 + tk_xxxyy_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyyzz[i] * fz_0;

        tk_xxxyyz_yyyzzz[i] = 3.0 * tk_xxxyy_yyyzz[i] * fe_0 + tk_xxxyy_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyzzz[i] * fz_0;

        tk_xxxyyz_yyzzzz[i] = 4.0 * tk_xxxyy_yyzzz[i] * fe_0 + tk_xxxyy_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyzzzz[i] * fz_0;

        tk_xxxyyz_yzzzzz[i] = 5.0 * tk_xxxyy_yzzzz[i] * fe_0 + tk_xxxyy_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yzzzzz[i] * fz_0;

        tk_xxxyyz_zzzzzz[i] = 6.0 * tk_xxxyy_zzzzz[i] * fe_0 + tk_xxxyy_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_zzzzzz[i] * fz_0;
    }

    // Set up 224-252 components of targeted buffer : II

    auto tk_xxxyzz_xxxxxx = pbuffer.data(idx_kin_ii + 224);

    auto tk_xxxyzz_xxxxxy = pbuffer.data(idx_kin_ii + 225);

    auto tk_xxxyzz_xxxxxz = pbuffer.data(idx_kin_ii + 226);

    auto tk_xxxyzz_xxxxyy = pbuffer.data(idx_kin_ii + 227);

    auto tk_xxxyzz_xxxxyz = pbuffer.data(idx_kin_ii + 228);

    auto tk_xxxyzz_xxxxzz = pbuffer.data(idx_kin_ii + 229);

    auto tk_xxxyzz_xxxyyy = pbuffer.data(idx_kin_ii + 230);

    auto tk_xxxyzz_xxxyyz = pbuffer.data(idx_kin_ii + 231);

    auto tk_xxxyzz_xxxyzz = pbuffer.data(idx_kin_ii + 232);

    auto tk_xxxyzz_xxxzzz = pbuffer.data(idx_kin_ii + 233);

    auto tk_xxxyzz_xxyyyy = pbuffer.data(idx_kin_ii + 234);

    auto tk_xxxyzz_xxyyyz = pbuffer.data(idx_kin_ii + 235);

    auto tk_xxxyzz_xxyyzz = pbuffer.data(idx_kin_ii + 236);

    auto tk_xxxyzz_xxyzzz = pbuffer.data(idx_kin_ii + 237);

    auto tk_xxxyzz_xxzzzz = pbuffer.data(idx_kin_ii + 238);

    auto tk_xxxyzz_xyyyyy = pbuffer.data(idx_kin_ii + 239);

    auto tk_xxxyzz_xyyyyz = pbuffer.data(idx_kin_ii + 240);

    auto tk_xxxyzz_xyyyzz = pbuffer.data(idx_kin_ii + 241);

    auto tk_xxxyzz_xyyzzz = pbuffer.data(idx_kin_ii + 242);

    auto tk_xxxyzz_xyzzzz = pbuffer.data(idx_kin_ii + 243);

    auto tk_xxxyzz_xzzzzz = pbuffer.data(idx_kin_ii + 244);

    auto tk_xxxyzz_yyyyyy = pbuffer.data(idx_kin_ii + 245);

    auto tk_xxxyzz_yyyyyz = pbuffer.data(idx_kin_ii + 246);

    auto tk_xxxyzz_yyyyzz = pbuffer.data(idx_kin_ii + 247);

    auto tk_xxxyzz_yyyzzz = pbuffer.data(idx_kin_ii + 248);

    auto tk_xxxyzz_yyzzzz = pbuffer.data(idx_kin_ii + 249);

    auto tk_xxxyzz_yzzzzz = pbuffer.data(idx_kin_ii + 250);

    auto tk_xxxyzz_zzzzzz = pbuffer.data(idx_kin_ii + 251);

#pragma omp simd aligned(pa_y,                 \
                             tk_xxxyzz_xxxxxx, \
                             tk_xxxyzz_xxxxxy, \
                             tk_xxxyzz_xxxxxz, \
                             tk_xxxyzz_xxxxyy, \
                             tk_xxxyzz_xxxxyz, \
                             tk_xxxyzz_xxxxzz, \
                             tk_xxxyzz_xxxyyy, \
                             tk_xxxyzz_xxxyyz, \
                             tk_xxxyzz_xxxyzz, \
                             tk_xxxyzz_xxxzzz, \
                             tk_xxxyzz_xxyyyy, \
                             tk_xxxyzz_xxyyyz, \
                             tk_xxxyzz_xxyyzz, \
                             tk_xxxyzz_xxyzzz, \
                             tk_xxxyzz_xxzzzz, \
                             tk_xxxyzz_xyyyyy, \
                             tk_xxxyzz_xyyyyz, \
                             tk_xxxyzz_xyyyzz, \
                             tk_xxxyzz_xyyzzz, \
                             tk_xxxyzz_xyzzzz, \
                             tk_xxxyzz_xzzzzz, \
                             tk_xxxyzz_yyyyyy, \
                             tk_xxxyzz_yyyyyz, \
                             tk_xxxyzz_yyyyzz, \
                             tk_xxxyzz_yyyzzz, \
                             tk_xxxyzz_yyzzzz, \
                             tk_xxxyzz_yzzzzz, \
                             tk_xxxyzz_zzzzzz, \
                             tk_xxxzz_xxxxx,   \
                             tk_xxxzz_xxxxxx,  \
                             tk_xxxzz_xxxxxy,  \
                             tk_xxxzz_xxxxxz,  \
                             tk_xxxzz_xxxxy,   \
                             tk_xxxzz_xxxxyy,  \
                             tk_xxxzz_xxxxyz,  \
                             tk_xxxzz_xxxxz,   \
                             tk_xxxzz_xxxxzz,  \
                             tk_xxxzz_xxxyy,   \
                             tk_xxxzz_xxxyyy,  \
                             tk_xxxzz_xxxyyz,  \
                             tk_xxxzz_xxxyz,   \
                             tk_xxxzz_xxxyzz,  \
                             tk_xxxzz_xxxzz,   \
                             tk_xxxzz_xxxzzz,  \
                             tk_xxxzz_xxyyy,   \
                             tk_xxxzz_xxyyyy,  \
                             tk_xxxzz_xxyyyz,  \
                             tk_xxxzz_xxyyz,   \
                             tk_xxxzz_xxyyzz,  \
                             tk_xxxzz_xxyzz,   \
                             tk_xxxzz_xxyzzz,  \
                             tk_xxxzz_xxzzz,   \
                             tk_xxxzz_xxzzzz,  \
                             tk_xxxzz_xyyyy,   \
                             tk_xxxzz_xyyyyy,  \
                             tk_xxxzz_xyyyyz,  \
                             tk_xxxzz_xyyyz,   \
                             tk_xxxzz_xyyyzz,  \
                             tk_xxxzz_xyyzz,   \
                             tk_xxxzz_xyyzzz,  \
                             tk_xxxzz_xyzzz,   \
                             tk_xxxzz_xyzzzz,  \
                             tk_xxxzz_xzzzz,   \
                             tk_xxxzz_xzzzzz,  \
                             tk_xxxzz_yyyyy,   \
                             tk_xxxzz_yyyyyy,  \
                             tk_xxxzz_yyyyyz,  \
                             tk_xxxzz_yyyyz,   \
                             tk_xxxzz_yyyyzz,  \
                             tk_xxxzz_yyyzz,   \
                             tk_xxxzz_yyyzzz,  \
                             tk_xxxzz_yyzzz,   \
                             tk_xxxzz_yyzzzz,  \
                             tk_xxxzz_yzzzz,   \
                             tk_xxxzz_yzzzzz,  \
                             tk_xxxzz_zzzzz,   \
                             tk_xxxzz_zzzzzz,  \
                             ts_xxxyzz_xxxxxx, \
                             ts_xxxyzz_xxxxxy, \
                             ts_xxxyzz_xxxxxz, \
                             ts_xxxyzz_xxxxyy, \
                             ts_xxxyzz_xxxxyz, \
                             ts_xxxyzz_xxxxzz, \
                             ts_xxxyzz_xxxyyy, \
                             ts_xxxyzz_xxxyyz, \
                             ts_xxxyzz_xxxyzz, \
                             ts_xxxyzz_xxxzzz, \
                             ts_xxxyzz_xxyyyy, \
                             ts_xxxyzz_xxyyyz, \
                             ts_xxxyzz_xxyyzz, \
                             ts_xxxyzz_xxyzzz, \
                             ts_xxxyzz_xxzzzz, \
                             ts_xxxyzz_xyyyyy, \
                             ts_xxxyzz_xyyyyz, \
                             ts_xxxyzz_xyyyzz, \
                             ts_xxxyzz_xyyzzz, \
                             ts_xxxyzz_xyzzzz, \
                             ts_xxxyzz_xzzzzz, \
                             ts_xxxyzz_yyyyyy, \
                             ts_xxxyzz_yyyyyz, \
                             ts_xxxyzz_yyyyzz, \
                             ts_xxxyzz_yyyzzz, \
                             ts_xxxyzz_yyzzzz, \
                             ts_xxxyzz_yzzzzz, \
                             ts_xxxyzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyzz_xxxxxx[i] = tk_xxxzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxxx[i] * fz_0;

        tk_xxxyzz_xxxxxy[i] = tk_xxxzz_xxxxx[i] * fe_0 + tk_xxxzz_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxxy[i] * fz_0;

        tk_xxxyzz_xxxxxz[i] = tk_xxxzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxxz[i] * fz_0;

        tk_xxxyzz_xxxxyy[i] = 2.0 * tk_xxxzz_xxxxy[i] * fe_0 + tk_xxxzz_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxyy[i] * fz_0;

        tk_xxxyzz_xxxxyz[i] = tk_xxxzz_xxxxz[i] * fe_0 + tk_xxxzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxyz[i] * fz_0;

        tk_xxxyzz_xxxxzz[i] = tk_xxxzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxzz[i] * fz_0;

        tk_xxxyzz_xxxyyy[i] = 3.0 * tk_xxxzz_xxxyy[i] * fe_0 + tk_xxxzz_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxyyy[i] * fz_0;

        tk_xxxyzz_xxxyyz[i] = 2.0 * tk_xxxzz_xxxyz[i] * fe_0 + tk_xxxzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxyyz[i] * fz_0;

        tk_xxxyzz_xxxyzz[i] = tk_xxxzz_xxxzz[i] * fe_0 + tk_xxxzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxyzz[i] * fz_0;

        tk_xxxyzz_xxxzzz[i] = tk_xxxzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxzzz[i] * fz_0;

        tk_xxxyzz_xxyyyy[i] = 4.0 * tk_xxxzz_xxyyy[i] * fe_0 + tk_xxxzz_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyyyy[i] * fz_0;

        tk_xxxyzz_xxyyyz[i] = 3.0 * tk_xxxzz_xxyyz[i] * fe_0 + tk_xxxzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyyyz[i] * fz_0;

        tk_xxxyzz_xxyyzz[i] = 2.0 * tk_xxxzz_xxyzz[i] * fe_0 + tk_xxxzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyyzz[i] * fz_0;

        tk_xxxyzz_xxyzzz[i] = tk_xxxzz_xxzzz[i] * fe_0 + tk_xxxzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyzzz[i] * fz_0;

        tk_xxxyzz_xxzzzz[i] = tk_xxxzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxzzzz[i] * fz_0;

        tk_xxxyzz_xyyyyy[i] = 5.0 * tk_xxxzz_xyyyy[i] * fe_0 + tk_xxxzz_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyyyy[i] * fz_0;

        tk_xxxyzz_xyyyyz[i] = 4.0 * tk_xxxzz_xyyyz[i] * fe_0 + tk_xxxzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyyyz[i] * fz_0;

        tk_xxxyzz_xyyyzz[i] = 3.0 * tk_xxxzz_xyyzz[i] * fe_0 + tk_xxxzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyyzz[i] * fz_0;

        tk_xxxyzz_xyyzzz[i] = 2.0 * tk_xxxzz_xyzzz[i] * fe_0 + tk_xxxzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyzzz[i] * fz_0;

        tk_xxxyzz_xyzzzz[i] = tk_xxxzz_xzzzz[i] * fe_0 + tk_xxxzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyzzzz[i] * fz_0;

        tk_xxxyzz_xzzzzz[i] = tk_xxxzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xzzzzz[i] * fz_0;

        tk_xxxyzz_yyyyyy[i] = 6.0 * tk_xxxzz_yyyyy[i] * fe_0 + tk_xxxzz_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyyyy[i] * fz_0;

        tk_xxxyzz_yyyyyz[i] = 5.0 * tk_xxxzz_yyyyz[i] * fe_0 + tk_xxxzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyyyz[i] * fz_0;

        tk_xxxyzz_yyyyzz[i] = 4.0 * tk_xxxzz_yyyzz[i] * fe_0 + tk_xxxzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyyzz[i] * fz_0;

        tk_xxxyzz_yyyzzz[i] = 3.0 * tk_xxxzz_yyzzz[i] * fe_0 + tk_xxxzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyzzz[i] * fz_0;

        tk_xxxyzz_yyzzzz[i] = 2.0 * tk_xxxzz_yzzzz[i] * fe_0 + tk_xxxzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyzzzz[i] * fz_0;

        tk_xxxyzz_yzzzzz[i] = tk_xxxzz_zzzzz[i] * fe_0 + tk_xxxzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yzzzzz[i] * fz_0;

        tk_xxxyzz_zzzzzz[i] = tk_xxxzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_zzzzzz[i] * fz_0;
    }

    // Set up 252-280 components of targeted buffer : II

    auto tk_xxxzzz_xxxxxx = pbuffer.data(idx_kin_ii + 252);

    auto tk_xxxzzz_xxxxxy = pbuffer.data(idx_kin_ii + 253);

    auto tk_xxxzzz_xxxxxz = pbuffer.data(idx_kin_ii + 254);

    auto tk_xxxzzz_xxxxyy = pbuffer.data(idx_kin_ii + 255);

    auto tk_xxxzzz_xxxxyz = pbuffer.data(idx_kin_ii + 256);

    auto tk_xxxzzz_xxxxzz = pbuffer.data(idx_kin_ii + 257);

    auto tk_xxxzzz_xxxyyy = pbuffer.data(idx_kin_ii + 258);

    auto tk_xxxzzz_xxxyyz = pbuffer.data(idx_kin_ii + 259);

    auto tk_xxxzzz_xxxyzz = pbuffer.data(idx_kin_ii + 260);

    auto tk_xxxzzz_xxxzzz = pbuffer.data(idx_kin_ii + 261);

    auto tk_xxxzzz_xxyyyy = pbuffer.data(idx_kin_ii + 262);

    auto tk_xxxzzz_xxyyyz = pbuffer.data(idx_kin_ii + 263);

    auto tk_xxxzzz_xxyyzz = pbuffer.data(idx_kin_ii + 264);

    auto tk_xxxzzz_xxyzzz = pbuffer.data(idx_kin_ii + 265);

    auto tk_xxxzzz_xxzzzz = pbuffer.data(idx_kin_ii + 266);

    auto tk_xxxzzz_xyyyyy = pbuffer.data(idx_kin_ii + 267);

    auto tk_xxxzzz_xyyyyz = pbuffer.data(idx_kin_ii + 268);

    auto tk_xxxzzz_xyyyzz = pbuffer.data(idx_kin_ii + 269);

    auto tk_xxxzzz_xyyzzz = pbuffer.data(idx_kin_ii + 270);

    auto tk_xxxzzz_xyzzzz = pbuffer.data(idx_kin_ii + 271);

    auto tk_xxxzzz_xzzzzz = pbuffer.data(idx_kin_ii + 272);

    auto tk_xxxzzz_yyyyyy = pbuffer.data(idx_kin_ii + 273);

    auto tk_xxxzzz_yyyyyz = pbuffer.data(idx_kin_ii + 274);

    auto tk_xxxzzz_yyyyzz = pbuffer.data(idx_kin_ii + 275);

    auto tk_xxxzzz_yyyzzz = pbuffer.data(idx_kin_ii + 276);

    auto tk_xxxzzz_yyzzzz = pbuffer.data(idx_kin_ii + 277);

    auto tk_xxxzzz_yzzzzz = pbuffer.data(idx_kin_ii + 278);

    auto tk_xxxzzz_zzzzzz = pbuffer.data(idx_kin_ii + 279);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tk_xxxz_xxxxxx,   \
                             tk_xxxz_xxxxxy,   \
                             tk_xxxz_xxxxyy,   \
                             tk_xxxz_xxxyyy,   \
                             tk_xxxz_xxyyyy,   \
                             tk_xxxz_xyyyyy,   \
                             tk_xxxzz_xxxxxx,  \
                             tk_xxxzz_xxxxxy,  \
                             tk_xxxzz_xxxxyy,  \
                             tk_xxxzz_xxxyyy,  \
                             tk_xxxzz_xxyyyy,  \
                             tk_xxxzz_xyyyyy,  \
                             tk_xxxzzz_xxxxxx, \
                             tk_xxxzzz_xxxxxy, \
                             tk_xxxzzz_xxxxxz, \
                             tk_xxxzzz_xxxxyy, \
                             tk_xxxzzz_xxxxyz, \
                             tk_xxxzzz_xxxxzz, \
                             tk_xxxzzz_xxxyyy, \
                             tk_xxxzzz_xxxyyz, \
                             tk_xxxzzz_xxxyzz, \
                             tk_xxxzzz_xxxzzz, \
                             tk_xxxzzz_xxyyyy, \
                             tk_xxxzzz_xxyyyz, \
                             tk_xxxzzz_xxyyzz, \
                             tk_xxxzzz_xxyzzz, \
                             tk_xxxzzz_xxzzzz, \
                             tk_xxxzzz_xyyyyy, \
                             tk_xxxzzz_xyyyyz, \
                             tk_xxxzzz_xyyyzz, \
                             tk_xxxzzz_xyyzzz, \
                             tk_xxxzzz_xyzzzz, \
                             tk_xxxzzz_xzzzzz, \
                             tk_xxxzzz_yyyyyy, \
                             tk_xxxzzz_yyyyyz, \
                             tk_xxxzzz_yyyyzz, \
                             tk_xxxzzz_yyyzzz, \
                             tk_xxxzzz_yyzzzz, \
                             tk_xxxzzz_yzzzzz, \
                             tk_xxxzzz_zzzzzz, \
                             tk_xxzzz_xxxxxz,  \
                             tk_xxzzz_xxxxyz,  \
                             tk_xxzzz_xxxxz,   \
                             tk_xxzzz_xxxxzz,  \
                             tk_xxzzz_xxxyyz,  \
                             tk_xxzzz_xxxyz,   \
                             tk_xxzzz_xxxyzz,  \
                             tk_xxzzz_xxxzz,   \
                             tk_xxzzz_xxxzzz,  \
                             tk_xxzzz_xxyyyz,  \
                             tk_xxzzz_xxyyz,   \
                             tk_xxzzz_xxyyzz,  \
                             tk_xxzzz_xxyzz,   \
                             tk_xxzzz_xxyzzz,  \
                             tk_xxzzz_xxzzz,   \
                             tk_xxzzz_xxzzzz,  \
                             tk_xxzzz_xyyyyz,  \
                             tk_xxzzz_xyyyz,   \
                             tk_xxzzz_xyyyzz,  \
                             tk_xxzzz_xyyzz,   \
                             tk_xxzzz_xyyzzz,  \
                             tk_xxzzz_xyzzz,   \
                             tk_xxzzz_xyzzzz,  \
                             tk_xxzzz_xzzzz,   \
                             tk_xxzzz_xzzzzz,  \
                             tk_xxzzz_yyyyyy,  \
                             tk_xxzzz_yyyyyz,  \
                             tk_xxzzz_yyyyz,   \
                             tk_xxzzz_yyyyzz,  \
                             tk_xxzzz_yyyzz,   \
                             tk_xxzzz_yyyzzz,  \
                             tk_xxzzz_yyzzz,   \
                             tk_xxzzz_yyzzzz,  \
                             tk_xxzzz_yzzzz,   \
                             tk_xxzzz_yzzzzz,  \
                             tk_xxzzz_zzzzz,   \
                             tk_xxzzz_zzzzzz,  \
                             tk_xzzz_xxxxxz,   \
                             tk_xzzz_xxxxyz,   \
                             tk_xzzz_xxxxzz,   \
                             tk_xzzz_xxxyyz,   \
                             tk_xzzz_xxxyzz,   \
                             tk_xzzz_xxxzzz,   \
                             tk_xzzz_xxyyyz,   \
                             tk_xzzz_xxyyzz,   \
                             tk_xzzz_xxyzzz,   \
                             tk_xzzz_xxzzzz,   \
                             tk_xzzz_xyyyyz,   \
                             tk_xzzz_xyyyzz,   \
                             tk_xzzz_xyyzzz,   \
                             tk_xzzz_xyzzzz,   \
                             tk_xzzz_xzzzzz,   \
                             tk_xzzz_yyyyyy,   \
                             tk_xzzz_yyyyyz,   \
                             tk_xzzz_yyyyzz,   \
                             tk_xzzz_yyyzzz,   \
                             tk_xzzz_yyzzzz,   \
                             tk_xzzz_yzzzzz,   \
                             tk_xzzz_zzzzzz,   \
                             ts_xxxz_xxxxxx,   \
                             ts_xxxz_xxxxxy,   \
                             ts_xxxz_xxxxyy,   \
                             ts_xxxz_xxxyyy,   \
                             ts_xxxz_xxyyyy,   \
                             ts_xxxz_xyyyyy,   \
                             ts_xxxzzz_xxxxxx, \
                             ts_xxxzzz_xxxxxy, \
                             ts_xxxzzz_xxxxxz, \
                             ts_xxxzzz_xxxxyy, \
                             ts_xxxzzz_xxxxyz, \
                             ts_xxxzzz_xxxxzz, \
                             ts_xxxzzz_xxxyyy, \
                             ts_xxxzzz_xxxyyz, \
                             ts_xxxzzz_xxxyzz, \
                             ts_xxxzzz_xxxzzz, \
                             ts_xxxzzz_xxyyyy, \
                             ts_xxxzzz_xxyyyz, \
                             ts_xxxzzz_xxyyzz, \
                             ts_xxxzzz_xxyzzz, \
                             ts_xxxzzz_xxzzzz, \
                             ts_xxxzzz_xyyyyy, \
                             ts_xxxzzz_xyyyyz, \
                             ts_xxxzzz_xyyyzz, \
                             ts_xxxzzz_xyyzzz, \
                             ts_xxxzzz_xyzzzz, \
                             ts_xxxzzz_xzzzzz, \
                             ts_xxxzzz_yyyyyy, \
                             ts_xxxzzz_yyyyyz, \
                             ts_xxxzzz_yyyyzz, \
                             ts_xxxzzz_yyyzzz, \
                             ts_xxxzzz_yyzzzz, \
                             ts_xxxzzz_yzzzzz, \
                             ts_xxxzzz_zzzzzz, \
                             ts_xzzz_xxxxxz,   \
                             ts_xzzz_xxxxyz,   \
                             ts_xzzz_xxxxzz,   \
                             ts_xzzz_xxxyyz,   \
                             ts_xzzz_xxxyzz,   \
                             ts_xzzz_xxxzzz,   \
                             ts_xzzz_xxyyyz,   \
                             ts_xzzz_xxyyzz,   \
                             ts_xzzz_xxyzzz,   \
                             ts_xzzz_xxzzzz,   \
                             ts_xzzz_xyyyyz,   \
                             ts_xzzz_xyyyzz,   \
                             ts_xzzz_xyyzzz,   \
                             ts_xzzz_xyzzzz,   \
                             ts_xzzz_xzzzzz,   \
                             ts_xzzz_yyyyyy,   \
                             ts_xzzz_yyyyyz,   \
                             ts_xzzz_yyyyzz,   \
                             ts_xzzz_yyyzzz,   \
                             ts_xzzz_yyzzzz,   \
                             ts_xzzz_yzzzzz,   \
                             ts_xzzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzzz_xxxxxx[i] = -4.0 * ts_xxxz_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxxxx[i] * fe_0 + tk_xxxzz_xxxxxx[i] * pa_z[i] +
                              2.0 * ts_xxxzzz_xxxxxx[i] * fz_0;

        tk_xxxzzz_xxxxxy[i] = -4.0 * ts_xxxz_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxxxy[i] * fe_0 + tk_xxxzz_xxxxxy[i] * pa_z[i] +
                              2.0 * ts_xxxzzz_xxxxxy[i] * fz_0;

        tk_xxxzzz_xxxxxz[i] = -4.0 * ts_xzzz_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxxxz[i] * fe_0 + 5.0 * tk_xxzzz_xxxxz[i] * fe_0 +
                              tk_xxzzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxxxz[i] * fz_0;

        tk_xxxzzz_xxxxyy[i] = -4.0 * ts_xxxz_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxxyy[i] * fe_0 + tk_xxxzz_xxxxyy[i] * pa_z[i] +
                              2.0 * ts_xxxzzz_xxxxyy[i] * fz_0;

        tk_xxxzzz_xxxxyz[i] = -4.0 * ts_xzzz_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxxyz[i] * fe_0 + 4.0 * tk_xxzzz_xxxyz[i] * fe_0 +
                              tk_xxzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxxyz[i] * fz_0;

        tk_xxxzzz_xxxxzz[i] = -4.0 * ts_xzzz_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxxzz[i] * fe_0 + 4.0 * tk_xxzzz_xxxzz[i] * fe_0 +
                              tk_xxzzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxxzz[i] * fz_0;

        tk_xxxzzz_xxxyyy[i] = -4.0 * ts_xxxz_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxyyy[i] * fe_0 + tk_xxxzz_xxxyyy[i] * pa_z[i] +
                              2.0 * ts_xxxzzz_xxxyyy[i] * fz_0;

        tk_xxxzzz_xxxyyz[i] = -4.0 * ts_xzzz_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxyyz[i] * fe_0 + 3.0 * tk_xxzzz_xxyyz[i] * fe_0 +
                              tk_xxzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxyyz[i] * fz_0;

        tk_xxxzzz_xxxyzz[i] = -4.0 * ts_xzzz_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxyzz[i] * fe_0 + 3.0 * tk_xxzzz_xxyzz[i] * fe_0 +
                              tk_xxzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxyzz[i] * fz_0;

        tk_xxxzzz_xxxzzz[i] = -4.0 * ts_xzzz_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxzzz[i] * fe_0 + 3.0 * tk_xxzzz_xxzzz[i] * fe_0 +
                              tk_xxzzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxzzz[i] * fz_0;

        tk_xxxzzz_xxyyyy[i] = -4.0 * ts_xxxz_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxyyyy[i] * fe_0 + tk_xxxzz_xxyyyy[i] * pa_z[i] +
                              2.0 * ts_xxxzzz_xxyyyy[i] * fz_0;

        tk_xxxzzz_xxyyyz[i] = -4.0 * ts_xzzz_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxyyyz[i] * fe_0 + 2.0 * tk_xxzzz_xyyyz[i] * fe_0 +
                              tk_xxzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxyyyz[i] * fz_0;

        tk_xxxzzz_xxyyzz[i] = -4.0 * ts_xzzz_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxyyzz[i] * fe_0 + 2.0 * tk_xxzzz_xyyzz[i] * fe_0 +
                              tk_xxzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxyyzz[i] * fz_0;

        tk_xxxzzz_xxyzzz[i] = -4.0 * ts_xzzz_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxyzzz[i] * fe_0 + 2.0 * tk_xxzzz_xyzzz[i] * fe_0 +
                              tk_xxzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxyzzz[i] * fz_0;

        tk_xxxzzz_xxzzzz[i] = -4.0 * ts_xzzz_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxzzzz[i] * fe_0 + 2.0 * tk_xxzzz_xzzzz[i] * fe_0 +
                              tk_xxzzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxzzzz[i] * fz_0;

        tk_xxxzzz_xyyyyy[i] = -4.0 * ts_xxxz_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xyyyyy[i] * fe_0 + tk_xxxzz_xyyyyy[i] * pa_z[i] +
                              2.0 * ts_xxxzzz_xyyyyy[i] * fz_0;

        tk_xxxzzz_xyyyyz[i] = -4.0 * ts_xzzz_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyyyyz[i] * fe_0 + tk_xxzzz_yyyyz[i] * fe_0 +
                              tk_xxzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyyyyz[i] * fz_0;

        tk_xxxzzz_xyyyzz[i] = -4.0 * ts_xzzz_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyyyzz[i] * fe_0 + tk_xxzzz_yyyzz[i] * fe_0 +
                              tk_xxzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyyyzz[i] * fz_0;

        tk_xxxzzz_xyyzzz[i] = -4.0 * ts_xzzz_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyyzzz[i] * fe_0 + tk_xxzzz_yyzzz[i] * fe_0 +
                              tk_xxzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyyzzz[i] * fz_0;

        tk_xxxzzz_xyzzzz[i] = -4.0 * ts_xzzz_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyzzzz[i] * fe_0 + tk_xxzzz_yzzzz[i] * fe_0 +
                              tk_xxzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyzzzz[i] * fz_0;

        tk_xxxzzz_xzzzzz[i] = -4.0 * ts_xzzz_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xzzzzz[i] * fe_0 + tk_xxzzz_zzzzz[i] * fe_0 +
                              tk_xxzzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xzzzzz[i] * fz_0;

        tk_xxxzzz_yyyyyy[i] = -4.0 * ts_xzzz_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyyyy[i] * fe_0 + tk_xxzzz_yyyyyy[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_yyyyyy[i] * fz_0;

        tk_xxxzzz_yyyyyz[i] = -4.0 * ts_xzzz_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyyyz[i] * fe_0 + tk_xxzzz_yyyyyz[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_yyyyyz[i] * fz_0;

        tk_xxxzzz_yyyyzz[i] = -4.0 * ts_xzzz_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyyzz[i] * fe_0 + tk_xxzzz_yyyyzz[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_yyyyzz[i] * fz_0;

        tk_xxxzzz_yyyzzz[i] = -4.0 * ts_xzzz_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyzzz[i] * fe_0 + tk_xxzzz_yyyzzz[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_yyyzzz[i] * fz_0;

        tk_xxxzzz_yyzzzz[i] = -4.0 * ts_xzzz_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyzzzz[i] * fe_0 + tk_xxzzz_yyzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_yyzzzz[i] * fz_0;

        tk_xxxzzz_yzzzzz[i] = -4.0 * ts_xzzz_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yzzzzz[i] * fe_0 + tk_xxzzz_yzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_yzzzzz[i] * fz_0;

        tk_xxxzzz_zzzzzz[i] = -4.0 * ts_xzzz_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_zzzzzz[i] * fe_0 + tk_xxzzz_zzzzzz[i] * pa_x[i] +
                              2.0 * ts_xxxzzz_zzzzzz[i] * fz_0;
    }

    // Set up 280-308 components of targeted buffer : II

    auto tk_xxyyyy_xxxxxx = pbuffer.data(idx_kin_ii + 280);

    auto tk_xxyyyy_xxxxxy = pbuffer.data(idx_kin_ii + 281);

    auto tk_xxyyyy_xxxxxz = pbuffer.data(idx_kin_ii + 282);

    auto tk_xxyyyy_xxxxyy = pbuffer.data(idx_kin_ii + 283);

    auto tk_xxyyyy_xxxxyz = pbuffer.data(idx_kin_ii + 284);

    auto tk_xxyyyy_xxxxzz = pbuffer.data(idx_kin_ii + 285);

    auto tk_xxyyyy_xxxyyy = pbuffer.data(idx_kin_ii + 286);

    auto tk_xxyyyy_xxxyyz = pbuffer.data(idx_kin_ii + 287);

    auto tk_xxyyyy_xxxyzz = pbuffer.data(idx_kin_ii + 288);

    auto tk_xxyyyy_xxxzzz = pbuffer.data(idx_kin_ii + 289);

    auto tk_xxyyyy_xxyyyy = pbuffer.data(idx_kin_ii + 290);

    auto tk_xxyyyy_xxyyyz = pbuffer.data(idx_kin_ii + 291);

    auto tk_xxyyyy_xxyyzz = pbuffer.data(idx_kin_ii + 292);

    auto tk_xxyyyy_xxyzzz = pbuffer.data(idx_kin_ii + 293);

    auto tk_xxyyyy_xxzzzz = pbuffer.data(idx_kin_ii + 294);

    auto tk_xxyyyy_xyyyyy = pbuffer.data(idx_kin_ii + 295);

    auto tk_xxyyyy_xyyyyz = pbuffer.data(idx_kin_ii + 296);

    auto tk_xxyyyy_xyyyzz = pbuffer.data(idx_kin_ii + 297);

    auto tk_xxyyyy_xyyzzz = pbuffer.data(idx_kin_ii + 298);

    auto tk_xxyyyy_xyzzzz = pbuffer.data(idx_kin_ii + 299);

    auto tk_xxyyyy_xzzzzz = pbuffer.data(idx_kin_ii + 300);

    auto tk_xxyyyy_yyyyyy = pbuffer.data(idx_kin_ii + 301);

    auto tk_xxyyyy_yyyyyz = pbuffer.data(idx_kin_ii + 302);

    auto tk_xxyyyy_yyyyzz = pbuffer.data(idx_kin_ii + 303);

    auto tk_xxyyyy_yyyzzz = pbuffer.data(idx_kin_ii + 304);

    auto tk_xxyyyy_yyzzzz = pbuffer.data(idx_kin_ii + 305);

    auto tk_xxyyyy_yzzzzz = pbuffer.data(idx_kin_ii + 306);

    auto tk_xxyyyy_zzzzzz = pbuffer.data(idx_kin_ii + 307);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tk_xxyy_xxxxxx,   \
                             tk_xxyy_xxxxxz,   \
                             tk_xxyy_xxxxzz,   \
                             tk_xxyy_xxxzzz,   \
                             tk_xxyy_xxzzzz,   \
                             tk_xxyy_xzzzzz,   \
                             tk_xxyyy_xxxxxx,  \
                             tk_xxyyy_xxxxxz,  \
                             tk_xxyyy_xxxxzz,  \
                             tk_xxyyy_xxxzzz,  \
                             tk_xxyyy_xxzzzz,  \
                             tk_xxyyy_xzzzzz,  \
                             tk_xxyyyy_xxxxxx, \
                             tk_xxyyyy_xxxxxy, \
                             tk_xxyyyy_xxxxxz, \
                             tk_xxyyyy_xxxxyy, \
                             tk_xxyyyy_xxxxyz, \
                             tk_xxyyyy_xxxxzz, \
                             tk_xxyyyy_xxxyyy, \
                             tk_xxyyyy_xxxyyz, \
                             tk_xxyyyy_xxxyzz, \
                             tk_xxyyyy_xxxzzz, \
                             tk_xxyyyy_xxyyyy, \
                             tk_xxyyyy_xxyyyz, \
                             tk_xxyyyy_xxyyzz, \
                             tk_xxyyyy_xxyzzz, \
                             tk_xxyyyy_xxzzzz, \
                             tk_xxyyyy_xyyyyy, \
                             tk_xxyyyy_xyyyyz, \
                             tk_xxyyyy_xyyyzz, \
                             tk_xxyyyy_xyyzzz, \
                             tk_xxyyyy_xyzzzz, \
                             tk_xxyyyy_xzzzzz, \
                             tk_xxyyyy_yyyyyy, \
                             tk_xxyyyy_yyyyyz, \
                             tk_xxyyyy_yyyyzz, \
                             tk_xxyyyy_yyyzzz, \
                             tk_xxyyyy_yyzzzz, \
                             tk_xxyyyy_yzzzzz, \
                             tk_xxyyyy_zzzzzz, \
                             tk_xyyyy_xxxxxy,  \
                             tk_xyyyy_xxxxy,   \
                             tk_xyyyy_xxxxyy,  \
                             tk_xyyyy_xxxxyz,  \
                             tk_xyyyy_xxxyy,   \
                             tk_xyyyy_xxxyyy,  \
                             tk_xyyyy_xxxyyz,  \
                             tk_xyyyy_xxxyz,   \
                             tk_xyyyy_xxxyzz,  \
                             tk_xyyyy_xxyyy,   \
                             tk_xyyyy_xxyyyy,  \
                             tk_xyyyy_xxyyyz,  \
                             tk_xyyyy_xxyyz,   \
                             tk_xyyyy_xxyyzz,  \
                             tk_xyyyy_xxyzz,   \
                             tk_xyyyy_xxyzzz,  \
                             tk_xyyyy_xyyyy,   \
                             tk_xyyyy_xyyyyy,  \
                             tk_xyyyy_xyyyyz,  \
                             tk_xyyyy_xyyyz,   \
                             tk_xyyyy_xyyyzz,  \
                             tk_xyyyy_xyyzz,   \
                             tk_xyyyy_xyyzzz,  \
                             tk_xyyyy_xyzzz,   \
                             tk_xyyyy_xyzzzz,  \
                             tk_xyyyy_yyyyy,   \
                             tk_xyyyy_yyyyyy,  \
                             tk_xyyyy_yyyyyz,  \
                             tk_xyyyy_yyyyz,   \
                             tk_xyyyy_yyyyzz,  \
                             tk_xyyyy_yyyzz,   \
                             tk_xyyyy_yyyzzz,  \
                             tk_xyyyy_yyzzz,   \
                             tk_xyyyy_yyzzzz,  \
                             tk_xyyyy_yzzzz,   \
                             tk_xyyyy_yzzzzz,  \
                             tk_xyyyy_zzzzzz,  \
                             tk_yyyy_xxxxxy,   \
                             tk_yyyy_xxxxyy,   \
                             tk_yyyy_xxxxyz,   \
                             tk_yyyy_xxxyyy,   \
                             tk_yyyy_xxxyyz,   \
                             tk_yyyy_xxxyzz,   \
                             tk_yyyy_xxyyyy,   \
                             tk_yyyy_xxyyyz,   \
                             tk_yyyy_xxyyzz,   \
                             tk_yyyy_xxyzzz,   \
                             tk_yyyy_xyyyyy,   \
                             tk_yyyy_xyyyyz,   \
                             tk_yyyy_xyyyzz,   \
                             tk_yyyy_xyyzzz,   \
                             tk_yyyy_xyzzzz,   \
                             tk_yyyy_yyyyyy,   \
                             tk_yyyy_yyyyyz,   \
                             tk_yyyy_yyyyzz,   \
                             tk_yyyy_yyyzzz,   \
                             tk_yyyy_yyzzzz,   \
                             tk_yyyy_yzzzzz,   \
                             tk_yyyy_zzzzzz,   \
                             ts_xxyy_xxxxxx,   \
                             ts_xxyy_xxxxxz,   \
                             ts_xxyy_xxxxzz,   \
                             ts_xxyy_xxxzzz,   \
                             ts_xxyy_xxzzzz,   \
                             ts_xxyy_xzzzzz,   \
                             ts_xxyyyy_xxxxxx, \
                             ts_xxyyyy_xxxxxy, \
                             ts_xxyyyy_xxxxxz, \
                             ts_xxyyyy_xxxxyy, \
                             ts_xxyyyy_xxxxyz, \
                             ts_xxyyyy_xxxxzz, \
                             ts_xxyyyy_xxxyyy, \
                             ts_xxyyyy_xxxyyz, \
                             ts_xxyyyy_xxxyzz, \
                             ts_xxyyyy_xxxzzz, \
                             ts_xxyyyy_xxyyyy, \
                             ts_xxyyyy_xxyyyz, \
                             ts_xxyyyy_xxyyzz, \
                             ts_xxyyyy_xxyzzz, \
                             ts_xxyyyy_xxzzzz, \
                             ts_xxyyyy_xyyyyy, \
                             ts_xxyyyy_xyyyyz, \
                             ts_xxyyyy_xyyyzz, \
                             ts_xxyyyy_xyyzzz, \
                             ts_xxyyyy_xyzzzz, \
                             ts_xxyyyy_xzzzzz, \
                             ts_xxyyyy_yyyyyy, \
                             ts_xxyyyy_yyyyyz, \
                             ts_xxyyyy_yyyyzz, \
                             ts_xxyyyy_yyyzzz, \
                             ts_xxyyyy_yyzzzz, \
                             ts_xxyyyy_yzzzzz, \
                             ts_xxyyyy_zzzzzz, \
                             ts_yyyy_xxxxxy,   \
                             ts_yyyy_xxxxyy,   \
                             ts_yyyy_xxxxyz,   \
                             ts_yyyy_xxxyyy,   \
                             ts_yyyy_xxxyyz,   \
                             ts_yyyy_xxxyzz,   \
                             ts_yyyy_xxyyyy,   \
                             ts_yyyy_xxyyyz,   \
                             ts_yyyy_xxyyzz,   \
                             ts_yyyy_xxyzzz,   \
                             ts_yyyy_xyyyyy,   \
                             ts_yyyy_xyyyyz,   \
                             ts_yyyy_xyyyzz,   \
                             ts_yyyy_xyyzzz,   \
                             ts_yyyy_xyzzzz,   \
                             ts_yyyy_yyyyyy,   \
                             ts_yyyy_yyyyyz,   \
                             ts_yyyy_yyyyzz,   \
                             ts_yyyy_yyyzzz,   \
                             ts_yyyy_yyzzzz,   \
                             ts_yyyy_yzzzzz,   \
                             ts_yyyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyyy_xxxxxx[i] = -6.0 * ts_xxyy_xxxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxxx[i] * fe_0 + tk_xxyyy_xxxxxx[i] * pa_y[i] +
                              2.0 * ts_xxyyyy_xxxxxx[i] * fz_0;

        tk_xxyyyy_xxxxxy[i] = -2.0 * ts_yyyy_xxxxxy[i] * fbe_0 * fz_0 + tk_yyyy_xxxxxy[i] * fe_0 + 5.0 * tk_xyyyy_xxxxy[i] * fe_0 +
                              tk_xyyyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxxxy[i] * fz_0;

        tk_xxyyyy_xxxxxz[i] = -6.0 * ts_xxyy_xxxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxxz[i] * fe_0 + tk_xxyyy_xxxxxz[i] * pa_y[i] +
                              2.0 * ts_xxyyyy_xxxxxz[i] * fz_0;

        tk_xxyyyy_xxxxyy[i] = -2.0 * ts_yyyy_xxxxyy[i] * fbe_0 * fz_0 + tk_yyyy_xxxxyy[i] * fe_0 + 4.0 * tk_xyyyy_xxxyy[i] * fe_0 +
                              tk_xyyyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxxyy[i] * fz_0;

        tk_xxyyyy_xxxxyz[i] = -2.0 * ts_yyyy_xxxxyz[i] * fbe_0 * fz_0 + tk_yyyy_xxxxyz[i] * fe_0 + 4.0 * tk_xyyyy_xxxyz[i] * fe_0 +
                              tk_xyyyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxxyz[i] * fz_0;

        tk_xxyyyy_xxxxzz[i] = -6.0 * ts_xxyy_xxxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxzz[i] * fe_0 + tk_xxyyy_xxxxzz[i] * pa_y[i] +
                              2.0 * ts_xxyyyy_xxxxzz[i] * fz_0;

        tk_xxyyyy_xxxyyy[i] = -2.0 * ts_yyyy_xxxyyy[i] * fbe_0 * fz_0 + tk_yyyy_xxxyyy[i] * fe_0 + 3.0 * tk_xyyyy_xxyyy[i] * fe_0 +
                              tk_xyyyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxyyy[i] * fz_0;

        tk_xxyyyy_xxxyyz[i] = -2.0 * ts_yyyy_xxxyyz[i] * fbe_0 * fz_0 + tk_yyyy_xxxyyz[i] * fe_0 + 3.0 * tk_xyyyy_xxyyz[i] * fe_0 +
                              tk_xyyyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxyyz[i] * fz_0;

        tk_xxyyyy_xxxyzz[i] = -2.0 * ts_yyyy_xxxyzz[i] * fbe_0 * fz_0 + tk_yyyy_xxxyzz[i] * fe_0 + 3.0 * tk_xyyyy_xxyzz[i] * fe_0 +
                              tk_xyyyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxyzz[i] * fz_0;

        tk_xxyyyy_xxxzzz[i] = -6.0 * ts_xxyy_xxxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxzzz[i] * fe_0 + tk_xxyyy_xxxzzz[i] * pa_y[i] +
                              2.0 * ts_xxyyyy_xxxzzz[i] * fz_0;

        tk_xxyyyy_xxyyyy[i] = -2.0 * ts_yyyy_xxyyyy[i] * fbe_0 * fz_0 + tk_yyyy_xxyyyy[i] * fe_0 + 2.0 * tk_xyyyy_xyyyy[i] * fe_0 +
                              tk_xyyyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyyyy[i] * fz_0;

        tk_xxyyyy_xxyyyz[i] = -2.0 * ts_yyyy_xxyyyz[i] * fbe_0 * fz_0 + tk_yyyy_xxyyyz[i] * fe_0 + 2.0 * tk_xyyyy_xyyyz[i] * fe_0 +
                              tk_xyyyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyyyz[i] * fz_0;

        tk_xxyyyy_xxyyzz[i] = -2.0 * ts_yyyy_xxyyzz[i] * fbe_0 * fz_0 + tk_yyyy_xxyyzz[i] * fe_0 + 2.0 * tk_xyyyy_xyyzz[i] * fe_0 +
                              tk_xyyyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyyzz[i] * fz_0;

        tk_xxyyyy_xxyzzz[i] = -2.0 * ts_yyyy_xxyzzz[i] * fbe_0 * fz_0 + tk_yyyy_xxyzzz[i] * fe_0 + 2.0 * tk_xyyyy_xyzzz[i] * fe_0 +
                              tk_xyyyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyzzz[i] * fz_0;

        tk_xxyyyy_xxzzzz[i] = -6.0 * ts_xxyy_xxzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxzzzz[i] * fe_0 + tk_xxyyy_xxzzzz[i] * pa_y[i] +
                              2.0 * ts_xxyyyy_xxzzzz[i] * fz_0;

        tk_xxyyyy_xyyyyy[i] = -2.0 * ts_yyyy_xyyyyy[i] * fbe_0 * fz_0 + tk_yyyy_xyyyyy[i] * fe_0 + tk_xyyyy_yyyyy[i] * fe_0 +
                              tk_xyyyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyyyy[i] * fz_0;

        tk_xxyyyy_xyyyyz[i] = -2.0 * ts_yyyy_xyyyyz[i] * fbe_0 * fz_0 + tk_yyyy_xyyyyz[i] * fe_0 + tk_xyyyy_yyyyz[i] * fe_0 +
                              tk_xyyyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyyyz[i] * fz_0;

        tk_xxyyyy_xyyyzz[i] = -2.0 * ts_yyyy_xyyyzz[i] * fbe_0 * fz_0 + tk_yyyy_xyyyzz[i] * fe_0 + tk_xyyyy_yyyzz[i] * fe_0 +
                              tk_xyyyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyyzz[i] * fz_0;

        tk_xxyyyy_xyyzzz[i] = -2.0 * ts_yyyy_xyyzzz[i] * fbe_0 * fz_0 + tk_yyyy_xyyzzz[i] * fe_0 + tk_xyyyy_yyzzz[i] * fe_0 +
                              tk_xyyyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyzzz[i] * fz_0;

        tk_xxyyyy_xyzzzz[i] = -2.0 * ts_yyyy_xyzzzz[i] * fbe_0 * fz_0 + tk_yyyy_xyzzzz[i] * fe_0 + tk_xyyyy_yzzzz[i] * fe_0 +
                              tk_xyyyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyzzzz[i] * fz_0;

        tk_xxyyyy_xzzzzz[i] = -6.0 * ts_xxyy_xzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xzzzzz[i] * fe_0 + tk_xxyyy_xzzzzz[i] * pa_y[i] +
                              2.0 * ts_xxyyyy_xzzzzz[i] * fz_0;

        tk_xxyyyy_yyyyyy[i] =
            -2.0 * ts_yyyy_yyyyyy[i] * fbe_0 * fz_0 + tk_yyyy_yyyyyy[i] * fe_0 + tk_xyyyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyyyy[i] * fz_0;

        tk_xxyyyy_yyyyyz[i] =
            -2.0 * ts_yyyy_yyyyyz[i] * fbe_0 * fz_0 + tk_yyyy_yyyyyz[i] * fe_0 + tk_xyyyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyyyz[i] * fz_0;

        tk_xxyyyy_yyyyzz[i] =
            -2.0 * ts_yyyy_yyyyzz[i] * fbe_0 * fz_0 + tk_yyyy_yyyyzz[i] * fe_0 + tk_xyyyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyyzz[i] * fz_0;

        tk_xxyyyy_yyyzzz[i] =
            -2.0 * ts_yyyy_yyyzzz[i] * fbe_0 * fz_0 + tk_yyyy_yyyzzz[i] * fe_0 + tk_xyyyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyzzz[i] * fz_0;

        tk_xxyyyy_yyzzzz[i] =
            -2.0 * ts_yyyy_yyzzzz[i] * fbe_0 * fz_0 + tk_yyyy_yyzzzz[i] * fe_0 + tk_xyyyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyzzzz[i] * fz_0;

        tk_xxyyyy_yzzzzz[i] =
            -2.0 * ts_yyyy_yzzzzz[i] * fbe_0 * fz_0 + tk_yyyy_yzzzzz[i] * fe_0 + tk_xyyyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yzzzzz[i] * fz_0;

        tk_xxyyyy_zzzzzz[i] =
            -2.0 * ts_yyyy_zzzzzz[i] * fbe_0 * fz_0 + tk_yyyy_zzzzzz[i] * fe_0 + tk_xyyyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_zzzzzz[i] * fz_0;
    }

    // Set up 308-336 components of targeted buffer : II

    auto tk_xxyyyz_xxxxxx = pbuffer.data(idx_kin_ii + 308);

    auto tk_xxyyyz_xxxxxy = pbuffer.data(idx_kin_ii + 309);

    auto tk_xxyyyz_xxxxxz = pbuffer.data(idx_kin_ii + 310);

    auto tk_xxyyyz_xxxxyy = pbuffer.data(idx_kin_ii + 311);

    auto tk_xxyyyz_xxxxyz = pbuffer.data(idx_kin_ii + 312);

    auto tk_xxyyyz_xxxxzz = pbuffer.data(idx_kin_ii + 313);

    auto tk_xxyyyz_xxxyyy = pbuffer.data(idx_kin_ii + 314);

    auto tk_xxyyyz_xxxyyz = pbuffer.data(idx_kin_ii + 315);

    auto tk_xxyyyz_xxxyzz = pbuffer.data(idx_kin_ii + 316);

    auto tk_xxyyyz_xxxzzz = pbuffer.data(idx_kin_ii + 317);

    auto tk_xxyyyz_xxyyyy = pbuffer.data(idx_kin_ii + 318);

    auto tk_xxyyyz_xxyyyz = pbuffer.data(idx_kin_ii + 319);

    auto tk_xxyyyz_xxyyzz = pbuffer.data(idx_kin_ii + 320);

    auto tk_xxyyyz_xxyzzz = pbuffer.data(idx_kin_ii + 321);

    auto tk_xxyyyz_xxzzzz = pbuffer.data(idx_kin_ii + 322);

    auto tk_xxyyyz_xyyyyy = pbuffer.data(idx_kin_ii + 323);

    auto tk_xxyyyz_xyyyyz = pbuffer.data(idx_kin_ii + 324);

    auto tk_xxyyyz_xyyyzz = pbuffer.data(idx_kin_ii + 325);

    auto tk_xxyyyz_xyyzzz = pbuffer.data(idx_kin_ii + 326);

    auto tk_xxyyyz_xyzzzz = pbuffer.data(idx_kin_ii + 327);

    auto tk_xxyyyz_xzzzzz = pbuffer.data(idx_kin_ii + 328);

    auto tk_xxyyyz_yyyyyy = pbuffer.data(idx_kin_ii + 329);

    auto tk_xxyyyz_yyyyyz = pbuffer.data(idx_kin_ii + 330);

    auto tk_xxyyyz_yyyyzz = pbuffer.data(idx_kin_ii + 331);

    auto tk_xxyyyz_yyyzzz = pbuffer.data(idx_kin_ii + 332);

    auto tk_xxyyyz_yyzzzz = pbuffer.data(idx_kin_ii + 333);

    auto tk_xxyyyz_yzzzzz = pbuffer.data(idx_kin_ii + 334);

    auto tk_xxyyyz_zzzzzz = pbuffer.data(idx_kin_ii + 335);

#pragma omp simd aligned(pa_z,                 \
                             tk_xxyyy_xxxxx,   \
                             tk_xxyyy_xxxxxx,  \
                             tk_xxyyy_xxxxxy,  \
                             tk_xxyyy_xxxxxz,  \
                             tk_xxyyy_xxxxy,   \
                             tk_xxyyy_xxxxyy,  \
                             tk_xxyyy_xxxxyz,  \
                             tk_xxyyy_xxxxz,   \
                             tk_xxyyy_xxxxzz,  \
                             tk_xxyyy_xxxyy,   \
                             tk_xxyyy_xxxyyy,  \
                             tk_xxyyy_xxxyyz,  \
                             tk_xxyyy_xxxyz,   \
                             tk_xxyyy_xxxyzz,  \
                             tk_xxyyy_xxxzz,   \
                             tk_xxyyy_xxxzzz,  \
                             tk_xxyyy_xxyyy,   \
                             tk_xxyyy_xxyyyy,  \
                             tk_xxyyy_xxyyyz,  \
                             tk_xxyyy_xxyyz,   \
                             tk_xxyyy_xxyyzz,  \
                             tk_xxyyy_xxyzz,   \
                             tk_xxyyy_xxyzzz,  \
                             tk_xxyyy_xxzzz,   \
                             tk_xxyyy_xxzzzz,  \
                             tk_xxyyy_xyyyy,   \
                             tk_xxyyy_xyyyyy,  \
                             tk_xxyyy_xyyyyz,  \
                             tk_xxyyy_xyyyz,   \
                             tk_xxyyy_xyyyzz,  \
                             tk_xxyyy_xyyzz,   \
                             tk_xxyyy_xyyzzz,  \
                             tk_xxyyy_xyzzz,   \
                             tk_xxyyy_xyzzzz,  \
                             tk_xxyyy_xzzzz,   \
                             tk_xxyyy_xzzzzz,  \
                             tk_xxyyy_yyyyy,   \
                             tk_xxyyy_yyyyyy,  \
                             tk_xxyyy_yyyyyz,  \
                             tk_xxyyy_yyyyz,   \
                             tk_xxyyy_yyyyzz,  \
                             tk_xxyyy_yyyzz,   \
                             tk_xxyyy_yyyzzz,  \
                             tk_xxyyy_yyzzz,   \
                             tk_xxyyy_yyzzzz,  \
                             tk_xxyyy_yzzzz,   \
                             tk_xxyyy_yzzzzz,  \
                             tk_xxyyy_zzzzz,   \
                             tk_xxyyy_zzzzzz,  \
                             tk_xxyyyz_xxxxxx, \
                             tk_xxyyyz_xxxxxy, \
                             tk_xxyyyz_xxxxxz, \
                             tk_xxyyyz_xxxxyy, \
                             tk_xxyyyz_xxxxyz, \
                             tk_xxyyyz_xxxxzz, \
                             tk_xxyyyz_xxxyyy, \
                             tk_xxyyyz_xxxyyz, \
                             tk_xxyyyz_xxxyzz, \
                             tk_xxyyyz_xxxzzz, \
                             tk_xxyyyz_xxyyyy, \
                             tk_xxyyyz_xxyyyz, \
                             tk_xxyyyz_xxyyzz, \
                             tk_xxyyyz_xxyzzz, \
                             tk_xxyyyz_xxzzzz, \
                             tk_xxyyyz_xyyyyy, \
                             tk_xxyyyz_xyyyyz, \
                             tk_xxyyyz_xyyyzz, \
                             tk_xxyyyz_xyyzzz, \
                             tk_xxyyyz_xyzzzz, \
                             tk_xxyyyz_xzzzzz, \
                             tk_xxyyyz_yyyyyy, \
                             tk_xxyyyz_yyyyyz, \
                             tk_xxyyyz_yyyyzz, \
                             tk_xxyyyz_yyyzzz, \
                             tk_xxyyyz_yyzzzz, \
                             tk_xxyyyz_yzzzzz, \
                             tk_xxyyyz_zzzzzz, \
                             ts_xxyyyz_xxxxxx, \
                             ts_xxyyyz_xxxxxy, \
                             ts_xxyyyz_xxxxxz, \
                             ts_xxyyyz_xxxxyy, \
                             ts_xxyyyz_xxxxyz, \
                             ts_xxyyyz_xxxxzz, \
                             ts_xxyyyz_xxxyyy, \
                             ts_xxyyyz_xxxyyz, \
                             ts_xxyyyz_xxxyzz, \
                             ts_xxyyyz_xxxzzz, \
                             ts_xxyyyz_xxyyyy, \
                             ts_xxyyyz_xxyyyz, \
                             ts_xxyyyz_xxyyzz, \
                             ts_xxyyyz_xxyzzz, \
                             ts_xxyyyz_xxzzzz, \
                             ts_xxyyyz_xyyyyy, \
                             ts_xxyyyz_xyyyyz, \
                             ts_xxyyyz_xyyyzz, \
                             ts_xxyyyz_xyyzzz, \
                             ts_xxyyyz_xyzzzz, \
                             ts_xxyyyz_xzzzzz, \
                             ts_xxyyyz_yyyyyy, \
                             ts_xxyyyz_yyyyyz, \
                             ts_xxyyyz_yyyyzz, \
                             ts_xxyyyz_yyyzzz, \
                             ts_xxyyyz_yyzzzz, \
                             ts_xxyyyz_yzzzzz, \
                             ts_xxyyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyyz_xxxxxx[i] = tk_xxyyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxxx[i] * fz_0;

        tk_xxyyyz_xxxxxy[i] = tk_xxyyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxxy[i] * fz_0;

        tk_xxyyyz_xxxxxz[i] = tk_xxyyy_xxxxx[i] * fe_0 + tk_xxyyy_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxxz[i] * fz_0;

        tk_xxyyyz_xxxxyy[i] = tk_xxyyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxyy[i] * fz_0;

        tk_xxyyyz_xxxxyz[i] = tk_xxyyy_xxxxy[i] * fe_0 + tk_xxyyy_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxyz[i] * fz_0;

        tk_xxyyyz_xxxxzz[i] = 2.0 * tk_xxyyy_xxxxz[i] * fe_0 + tk_xxyyy_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxzz[i] * fz_0;

        tk_xxyyyz_xxxyyy[i] = tk_xxyyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxyyy[i] * fz_0;

        tk_xxyyyz_xxxyyz[i] = tk_xxyyy_xxxyy[i] * fe_0 + tk_xxyyy_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxyyz[i] * fz_0;

        tk_xxyyyz_xxxyzz[i] = 2.0 * tk_xxyyy_xxxyz[i] * fe_0 + tk_xxyyy_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxyzz[i] * fz_0;

        tk_xxyyyz_xxxzzz[i] = 3.0 * tk_xxyyy_xxxzz[i] * fe_0 + tk_xxyyy_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxzzz[i] * fz_0;

        tk_xxyyyz_xxyyyy[i] = tk_xxyyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyyyy[i] * fz_0;

        tk_xxyyyz_xxyyyz[i] = tk_xxyyy_xxyyy[i] * fe_0 + tk_xxyyy_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyyyz[i] * fz_0;

        tk_xxyyyz_xxyyzz[i] = 2.0 * tk_xxyyy_xxyyz[i] * fe_0 + tk_xxyyy_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyyzz[i] * fz_0;

        tk_xxyyyz_xxyzzz[i] = 3.0 * tk_xxyyy_xxyzz[i] * fe_0 + tk_xxyyy_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyzzz[i] * fz_0;

        tk_xxyyyz_xxzzzz[i] = 4.0 * tk_xxyyy_xxzzz[i] * fe_0 + tk_xxyyy_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxzzzz[i] * fz_0;

        tk_xxyyyz_xyyyyy[i] = tk_xxyyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyyyy[i] * fz_0;

        tk_xxyyyz_xyyyyz[i] = tk_xxyyy_xyyyy[i] * fe_0 + tk_xxyyy_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyyyz[i] * fz_0;

        tk_xxyyyz_xyyyzz[i] = 2.0 * tk_xxyyy_xyyyz[i] * fe_0 + tk_xxyyy_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyyzz[i] * fz_0;

        tk_xxyyyz_xyyzzz[i] = 3.0 * tk_xxyyy_xyyzz[i] * fe_0 + tk_xxyyy_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyzzz[i] * fz_0;

        tk_xxyyyz_xyzzzz[i] = 4.0 * tk_xxyyy_xyzzz[i] * fe_0 + tk_xxyyy_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyzzzz[i] * fz_0;

        tk_xxyyyz_xzzzzz[i] = 5.0 * tk_xxyyy_xzzzz[i] * fe_0 + tk_xxyyy_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xzzzzz[i] * fz_0;

        tk_xxyyyz_yyyyyy[i] = tk_xxyyy_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyyyy[i] * fz_0;

        tk_xxyyyz_yyyyyz[i] = tk_xxyyy_yyyyy[i] * fe_0 + tk_xxyyy_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyyyz[i] * fz_0;

        tk_xxyyyz_yyyyzz[i] = 2.0 * tk_xxyyy_yyyyz[i] * fe_0 + tk_xxyyy_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyyzz[i] * fz_0;

        tk_xxyyyz_yyyzzz[i] = 3.0 * tk_xxyyy_yyyzz[i] * fe_0 + tk_xxyyy_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyzzz[i] * fz_0;

        tk_xxyyyz_yyzzzz[i] = 4.0 * tk_xxyyy_yyzzz[i] * fe_0 + tk_xxyyy_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyzzzz[i] * fz_0;

        tk_xxyyyz_yzzzzz[i] = 5.0 * tk_xxyyy_yzzzz[i] * fe_0 + tk_xxyyy_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yzzzzz[i] * fz_0;

        tk_xxyyyz_zzzzzz[i] = 6.0 * tk_xxyyy_zzzzz[i] * fe_0 + tk_xxyyy_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_zzzzzz[i] * fz_0;
    }

    // Set up 336-364 components of targeted buffer : II

    auto tk_xxyyzz_xxxxxx = pbuffer.data(idx_kin_ii + 336);

    auto tk_xxyyzz_xxxxxy = pbuffer.data(idx_kin_ii + 337);

    auto tk_xxyyzz_xxxxxz = pbuffer.data(idx_kin_ii + 338);

    auto tk_xxyyzz_xxxxyy = pbuffer.data(idx_kin_ii + 339);

    auto tk_xxyyzz_xxxxyz = pbuffer.data(idx_kin_ii + 340);

    auto tk_xxyyzz_xxxxzz = pbuffer.data(idx_kin_ii + 341);

    auto tk_xxyyzz_xxxyyy = pbuffer.data(idx_kin_ii + 342);

    auto tk_xxyyzz_xxxyyz = pbuffer.data(idx_kin_ii + 343);

    auto tk_xxyyzz_xxxyzz = pbuffer.data(idx_kin_ii + 344);

    auto tk_xxyyzz_xxxzzz = pbuffer.data(idx_kin_ii + 345);

    auto tk_xxyyzz_xxyyyy = pbuffer.data(idx_kin_ii + 346);

    auto tk_xxyyzz_xxyyyz = pbuffer.data(idx_kin_ii + 347);

    auto tk_xxyyzz_xxyyzz = pbuffer.data(idx_kin_ii + 348);

    auto tk_xxyyzz_xxyzzz = pbuffer.data(idx_kin_ii + 349);

    auto tk_xxyyzz_xxzzzz = pbuffer.data(idx_kin_ii + 350);

    auto tk_xxyyzz_xyyyyy = pbuffer.data(idx_kin_ii + 351);

    auto tk_xxyyzz_xyyyyz = pbuffer.data(idx_kin_ii + 352);

    auto tk_xxyyzz_xyyyzz = pbuffer.data(idx_kin_ii + 353);

    auto tk_xxyyzz_xyyzzz = pbuffer.data(idx_kin_ii + 354);

    auto tk_xxyyzz_xyzzzz = pbuffer.data(idx_kin_ii + 355);

    auto tk_xxyyzz_xzzzzz = pbuffer.data(idx_kin_ii + 356);

    auto tk_xxyyzz_yyyyyy = pbuffer.data(idx_kin_ii + 357);

    auto tk_xxyyzz_yyyyyz = pbuffer.data(idx_kin_ii + 358);

    auto tk_xxyyzz_yyyyzz = pbuffer.data(idx_kin_ii + 359);

    auto tk_xxyyzz_yyyzzz = pbuffer.data(idx_kin_ii + 360);

    auto tk_xxyyzz_yyzzzz = pbuffer.data(idx_kin_ii + 361);

    auto tk_xxyyzz_yzzzzz = pbuffer.data(idx_kin_ii + 362);

    auto tk_xxyyzz_zzzzzz = pbuffer.data(idx_kin_ii + 363);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tk_xxyy_xxxxxy,   \
                             tk_xxyy_xxxxyy,   \
                             tk_xxyy_xxxyyy,   \
                             tk_xxyy_xxyyyy,   \
                             tk_xxyy_xyyyyy,   \
                             tk_xxyyz_xxxxxy,  \
                             tk_xxyyz_xxxxyy,  \
                             tk_xxyyz_xxxyyy,  \
                             tk_xxyyz_xxyyyy,  \
                             tk_xxyyz_xyyyyy,  \
                             tk_xxyyzz_xxxxxx, \
                             tk_xxyyzz_xxxxxy, \
                             tk_xxyyzz_xxxxxz, \
                             tk_xxyyzz_xxxxyy, \
                             tk_xxyyzz_xxxxyz, \
                             tk_xxyyzz_xxxxzz, \
                             tk_xxyyzz_xxxyyy, \
                             tk_xxyyzz_xxxyyz, \
                             tk_xxyyzz_xxxyzz, \
                             tk_xxyyzz_xxxzzz, \
                             tk_xxyyzz_xxyyyy, \
                             tk_xxyyzz_xxyyyz, \
                             tk_xxyyzz_xxyyzz, \
                             tk_xxyyzz_xxyzzz, \
                             tk_xxyyzz_xxzzzz, \
                             tk_xxyyzz_xyyyyy, \
                             tk_xxyyzz_xyyyyz, \
                             tk_xxyyzz_xyyyzz, \
                             tk_xxyyzz_xyyzzz, \
                             tk_xxyyzz_xyzzzz, \
                             tk_xxyyzz_xzzzzz, \
                             tk_xxyyzz_yyyyyy, \
                             tk_xxyyzz_yyyyyz, \
                             tk_xxyyzz_yyyyzz, \
                             tk_xxyyzz_yyyzzz, \
                             tk_xxyyzz_yyzzzz, \
                             tk_xxyyzz_yzzzzz, \
                             tk_xxyyzz_zzzzzz, \
                             tk_xxyzz_xxxxxx,  \
                             tk_xxyzz_xxxxxz,  \
                             tk_xxyzz_xxxxzz,  \
                             tk_xxyzz_xxxzzz,  \
                             tk_xxyzz_xxzzzz,  \
                             tk_xxyzz_xzzzzz,  \
                             tk_xxzz_xxxxxx,   \
                             tk_xxzz_xxxxxz,   \
                             tk_xxzz_xxxxzz,   \
                             tk_xxzz_xxxzzz,   \
                             tk_xxzz_xxzzzz,   \
                             tk_xxzz_xzzzzz,   \
                             tk_xyyzz_xxxxyz,  \
                             tk_xyyzz_xxxyyz,  \
                             tk_xyyzz_xxxyz,   \
                             tk_xyyzz_xxxyzz,  \
                             tk_xyyzz_xxyyyz,  \
                             tk_xyyzz_xxyyz,   \
                             tk_xyyzz_xxyyzz,  \
                             tk_xyyzz_xxyzz,   \
                             tk_xyyzz_xxyzzz,  \
                             tk_xyyzz_xyyyyz,  \
                             tk_xyyzz_xyyyz,   \
                             tk_xyyzz_xyyyzz,  \
                             tk_xyyzz_xyyzz,   \
                             tk_xyyzz_xyyzzz,  \
                             tk_xyyzz_xyzzz,   \
                             tk_xyyzz_xyzzzz,  \
                             tk_xyyzz_yyyyyy,  \
                             tk_xyyzz_yyyyyz,  \
                             tk_xyyzz_yyyyz,   \
                             tk_xyyzz_yyyyzz,  \
                             tk_xyyzz_yyyzz,   \
                             tk_xyyzz_yyyzzz,  \
                             tk_xyyzz_yyzzz,   \
                             tk_xyyzz_yyzzzz,  \
                             tk_xyyzz_yzzzz,   \
                             tk_xyyzz_yzzzzz,  \
                             tk_xyyzz_zzzzzz,  \
                             tk_yyzz_xxxxyz,   \
                             tk_yyzz_xxxyyz,   \
                             tk_yyzz_xxxyzz,   \
                             tk_yyzz_xxyyyz,   \
                             tk_yyzz_xxyyzz,   \
                             tk_yyzz_xxyzzz,   \
                             tk_yyzz_xyyyyz,   \
                             tk_yyzz_xyyyzz,   \
                             tk_yyzz_xyyzzz,   \
                             tk_yyzz_xyzzzz,   \
                             tk_yyzz_yyyyyy,   \
                             tk_yyzz_yyyyyz,   \
                             tk_yyzz_yyyyzz,   \
                             tk_yyzz_yyyzzz,   \
                             tk_yyzz_yyzzzz,   \
                             tk_yyzz_yzzzzz,   \
                             tk_yyzz_zzzzzz,   \
                             ts_xxyy_xxxxxy,   \
                             ts_xxyy_xxxxyy,   \
                             ts_xxyy_xxxyyy,   \
                             ts_xxyy_xxyyyy,   \
                             ts_xxyy_xyyyyy,   \
                             ts_xxyyzz_xxxxxx, \
                             ts_xxyyzz_xxxxxy, \
                             ts_xxyyzz_xxxxxz, \
                             ts_xxyyzz_xxxxyy, \
                             ts_xxyyzz_xxxxyz, \
                             ts_xxyyzz_xxxxzz, \
                             ts_xxyyzz_xxxyyy, \
                             ts_xxyyzz_xxxyyz, \
                             ts_xxyyzz_xxxyzz, \
                             ts_xxyyzz_xxxzzz, \
                             ts_xxyyzz_xxyyyy, \
                             ts_xxyyzz_xxyyyz, \
                             ts_xxyyzz_xxyyzz, \
                             ts_xxyyzz_xxyzzz, \
                             ts_xxyyzz_xxzzzz, \
                             ts_xxyyzz_xyyyyy, \
                             ts_xxyyzz_xyyyyz, \
                             ts_xxyyzz_xyyyzz, \
                             ts_xxyyzz_xyyzzz, \
                             ts_xxyyzz_xyzzzz, \
                             ts_xxyyzz_xzzzzz, \
                             ts_xxyyzz_yyyyyy, \
                             ts_xxyyzz_yyyyyz, \
                             ts_xxyyzz_yyyyzz, \
                             ts_xxyyzz_yyyzzz, \
                             ts_xxyyzz_yyzzzz, \
                             ts_xxyyzz_yzzzzz, \
                             ts_xxyyzz_zzzzzz, \
                             ts_xxzz_xxxxxx,   \
                             ts_xxzz_xxxxxz,   \
                             ts_xxzz_xxxxzz,   \
                             ts_xxzz_xxxzzz,   \
                             ts_xxzz_xxzzzz,   \
                             ts_xxzz_xzzzzz,   \
                             ts_yyzz_xxxxyz,   \
                             ts_yyzz_xxxyyz,   \
                             ts_yyzz_xxxyzz,   \
                             ts_yyzz_xxyyyz,   \
                             ts_yyzz_xxyyzz,   \
                             ts_yyzz_xxyzzz,   \
                             ts_yyzz_xyyyyz,   \
                             ts_yyzz_xyyyzz,   \
                             ts_yyzz_xyyzzz,   \
                             ts_yyzz_xyzzzz,   \
                             ts_yyzz_yyyyyy,   \
                             ts_yyzz_yyyyyz,   \
                             ts_yyzz_yyyyzz,   \
                             ts_yyzz_yyyzzz,   \
                             ts_yyzz_yyzzzz,   \
                             ts_yyzz_yzzzzz,   \
                             ts_yyzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyzz_xxxxxx[i] =
            -2.0 * ts_xxzz_xxxxxx[i] * fbe_0 * fz_0 + tk_xxzz_xxxxxx[i] * fe_0 + tk_xxyzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxxxx[i] * fz_0;

        tk_xxyyzz_xxxxxy[i] =
            -2.0 * ts_xxyy_xxxxxy[i] * fbe_0 * fz_0 + tk_xxyy_xxxxxy[i] * fe_0 + tk_xxyyz_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxxxxy[i] * fz_0;

        tk_xxyyzz_xxxxxz[i] =
            -2.0 * ts_xxzz_xxxxxz[i] * fbe_0 * fz_0 + tk_xxzz_xxxxxz[i] * fe_0 + tk_xxyzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxxxz[i] * fz_0;

        tk_xxyyzz_xxxxyy[i] =
            -2.0 * ts_xxyy_xxxxyy[i] * fbe_0 * fz_0 + tk_xxyy_xxxxyy[i] * fe_0 + tk_xxyyz_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxxxyy[i] * fz_0;

        tk_xxyyzz_xxxxyz[i] = -2.0 * ts_yyzz_xxxxyz[i] * fbe_0 * fz_0 + tk_yyzz_xxxxyz[i] * fe_0 + 4.0 * tk_xyyzz_xxxyz[i] * fe_0 +
                              tk_xyyzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxxxyz[i] * fz_0;

        tk_xxyyzz_xxxxzz[i] =
            -2.0 * ts_xxzz_xxxxzz[i] * fbe_0 * fz_0 + tk_xxzz_xxxxzz[i] * fe_0 + tk_xxyzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxxzz[i] * fz_0;

        tk_xxyyzz_xxxyyy[i] =
            -2.0 * ts_xxyy_xxxyyy[i] * fbe_0 * fz_0 + tk_xxyy_xxxyyy[i] * fe_0 + tk_xxyyz_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxxyyy[i] * fz_0;

        tk_xxyyzz_xxxyyz[i] = -2.0 * ts_yyzz_xxxyyz[i] * fbe_0 * fz_0 + tk_yyzz_xxxyyz[i] * fe_0 + 3.0 * tk_xyyzz_xxyyz[i] * fe_0 +
                              tk_xyyzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxxyyz[i] * fz_0;

        tk_xxyyzz_xxxyzz[i] = -2.0 * ts_yyzz_xxxyzz[i] * fbe_0 * fz_0 + tk_yyzz_xxxyzz[i] * fe_0 + 3.0 * tk_xyyzz_xxyzz[i] * fe_0 +
                              tk_xyyzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxxyzz[i] * fz_0;

        tk_xxyyzz_xxxzzz[i] =
            -2.0 * ts_xxzz_xxxzzz[i] * fbe_0 * fz_0 + tk_xxzz_xxxzzz[i] * fe_0 + tk_xxyzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxzzz[i] * fz_0;

        tk_xxyyzz_xxyyyy[i] =
            -2.0 * ts_xxyy_xxyyyy[i] * fbe_0 * fz_0 + tk_xxyy_xxyyyy[i] * fe_0 + tk_xxyyz_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxyyyy[i] * fz_0;

        tk_xxyyzz_xxyyyz[i] = -2.0 * ts_yyzz_xxyyyz[i] * fbe_0 * fz_0 + tk_yyzz_xxyyyz[i] * fe_0 + 2.0 * tk_xyyzz_xyyyz[i] * fe_0 +
                              tk_xyyzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxyyyz[i] * fz_0;

        tk_xxyyzz_xxyyzz[i] = -2.0 * ts_yyzz_xxyyzz[i] * fbe_0 * fz_0 + tk_yyzz_xxyyzz[i] * fe_0 + 2.0 * tk_xyyzz_xyyzz[i] * fe_0 +
                              tk_xyyzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxyyzz[i] * fz_0;

        tk_xxyyzz_xxyzzz[i] = -2.0 * ts_yyzz_xxyzzz[i] * fbe_0 * fz_0 + tk_yyzz_xxyzzz[i] * fe_0 + 2.0 * tk_xyyzz_xyzzz[i] * fe_0 +
                              tk_xyyzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxyzzz[i] * fz_0;

        tk_xxyyzz_xxzzzz[i] =
            -2.0 * ts_xxzz_xxzzzz[i] * fbe_0 * fz_0 + tk_xxzz_xxzzzz[i] * fe_0 + tk_xxyzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxzzzz[i] * fz_0;

        tk_xxyyzz_xyyyyy[i] =
            -2.0 * ts_xxyy_xyyyyy[i] * fbe_0 * fz_0 + tk_xxyy_xyyyyy[i] * fe_0 + tk_xxyyz_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xyyyyy[i] * fz_0;

        tk_xxyyzz_xyyyyz[i] = -2.0 * ts_yyzz_xyyyyz[i] * fbe_0 * fz_0 + tk_yyzz_xyyyyz[i] * fe_0 + tk_xyyzz_yyyyz[i] * fe_0 +
                              tk_xyyzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyyyyz[i] * fz_0;

        tk_xxyyzz_xyyyzz[i] = -2.0 * ts_yyzz_xyyyzz[i] * fbe_0 * fz_0 + tk_yyzz_xyyyzz[i] * fe_0 + tk_xyyzz_yyyzz[i] * fe_0 +
                              tk_xyyzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyyyzz[i] * fz_0;

        tk_xxyyzz_xyyzzz[i] = -2.0 * ts_yyzz_xyyzzz[i] * fbe_0 * fz_0 + tk_yyzz_xyyzzz[i] * fe_0 + tk_xyyzz_yyzzz[i] * fe_0 +
                              tk_xyyzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyyzzz[i] * fz_0;

        tk_xxyyzz_xyzzzz[i] = -2.0 * ts_yyzz_xyzzzz[i] * fbe_0 * fz_0 + tk_yyzz_xyzzzz[i] * fe_0 + tk_xyyzz_yzzzz[i] * fe_0 +
                              tk_xyyzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyzzzz[i] * fz_0;

        tk_xxyyzz_xzzzzz[i] =
            -2.0 * ts_xxzz_xzzzzz[i] * fbe_0 * fz_0 + tk_xxzz_xzzzzz[i] * fe_0 + tk_xxyzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xzzzzz[i] * fz_0;

        tk_xxyyzz_yyyyyy[i] =
            -2.0 * ts_yyzz_yyyyyy[i] * fbe_0 * fz_0 + tk_yyzz_yyyyyy[i] * fe_0 + tk_xyyzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyyyy[i] * fz_0;

        tk_xxyyzz_yyyyyz[i] =
            -2.0 * ts_yyzz_yyyyyz[i] * fbe_0 * fz_0 + tk_yyzz_yyyyyz[i] * fe_0 + tk_xyyzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyyyz[i] * fz_0;

        tk_xxyyzz_yyyyzz[i] =
            -2.0 * ts_yyzz_yyyyzz[i] * fbe_0 * fz_0 + tk_yyzz_yyyyzz[i] * fe_0 + tk_xyyzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyyzz[i] * fz_0;

        tk_xxyyzz_yyyzzz[i] =
            -2.0 * ts_yyzz_yyyzzz[i] * fbe_0 * fz_0 + tk_yyzz_yyyzzz[i] * fe_0 + tk_xyyzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyzzz[i] * fz_0;

        tk_xxyyzz_yyzzzz[i] =
            -2.0 * ts_yyzz_yyzzzz[i] * fbe_0 * fz_0 + tk_yyzz_yyzzzz[i] * fe_0 + tk_xyyzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyzzzz[i] * fz_0;

        tk_xxyyzz_yzzzzz[i] =
            -2.0 * ts_yyzz_yzzzzz[i] * fbe_0 * fz_0 + tk_yyzz_yzzzzz[i] * fe_0 + tk_xyyzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yzzzzz[i] * fz_0;

        tk_xxyyzz_zzzzzz[i] =
            -2.0 * ts_yyzz_zzzzzz[i] * fbe_0 * fz_0 + tk_yyzz_zzzzzz[i] * fe_0 + tk_xyyzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_zzzzzz[i] * fz_0;
    }

    // Set up 364-392 components of targeted buffer : II

    auto tk_xxyzzz_xxxxxx = pbuffer.data(idx_kin_ii + 364);

    auto tk_xxyzzz_xxxxxy = pbuffer.data(idx_kin_ii + 365);

    auto tk_xxyzzz_xxxxxz = pbuffer.data(idx_kin_ii + 366);

    auto tk_xxyzzz_xxxxyy = pbuffer.data(idx_kin_ii + 367);

    auto tk_xxyzzz_xxxxyz = pbuffer.data(idx_kin_ii + 368);

    auto tk_xxyzzz_xxxxzz = pbuffer.data(idx_kin_ii + 369);

    auto tk_xxyzzz_xxxyyy = pbuffer.data(idx_kin_ii + 370);

    auto tk_xxyzzz_xxxyyz = pbuffer.data(idx_kin_ii + 371);

    auto tk_xxyzzz_xxxyzz = pbuffer.data(idx_kin_ii + 372);

    auto tk_xxyzzz_xxxzzz = pbuffer.data(idx_kin_ii + 373);

    auto tk_xxyzzz_xxyyyy = pbuffer.data(idx_kin_ii + 374);

    auto tk_xxyzzz_xxyyyz = pbuffer.data(idx_kin_ii + 375);

    auto tk_xxyzzz_xxyyzz = pbuffer.data(idx_kin_ii + 376);

    auto tk_xxyzzz_xxyzzz = pbuffer.data(idx_kin_ii + 377);

    auto tk_xxyzzz_xxzzzz = pbuffer.data(idx_kin_ii + 378);

    auto tk_xxyzzz_xyyyyy = pbuffer.data(idx_kin_ii + 379);

    auto tk_xxyzzz_xyyyyz = pbuffer.data(idx_kin_ii + 380);

    auto tk_xxyzzz_xyyyzz = pbuffer.data(idx_kin_ii + 381);

    auto tk_xxyzzz_xyyzzz = pbuffer.data(idx_kin_ii + 382);

    auto tk_xxyzzz_xyzzzz = pbuffer.data(idx_kin_ii + 383);

    auto tk_xxyzzz_xzzzzz = pbuffer.data(idx_kin_ii + 384);

    auto tk_xxyzzz_yyyyyy = pbuffer.data(idx_kin_ii + 385);

    auto tk_xxyzzz_yyyyyz = pbuffer.data(idx_kin_ii + 386);

    auto tk_xxyzzz_yyyyzz = pbuffer.data(idx_kin_ii + 387);

    auto tk_xxyzzz_yyyzzz = pbuffer.data(idx_kin_ii + 388);

    auto tk_xxyzzz_yyzzzz = pbuffer.data(idx_kin_ii + 389);

    auto tk_xxyzzz_yzzzzz = pbuffer.data(idx_kin_ii + 390);

    auto tk_xxyzzz_zzzzzz = pbuffer.data(idx_kin_ii + 391);

#pragma omp simd aligned(pa_y,                 \
                             tk_xxyzzz_xxxxxx, \
                             tk_xxyzzz_xxxxxy, \
                             tk_xxyzzz_xxxxxz, \
                             tk_xxyzzz_xxxxyy, \
                             tk_xxyzzz_xxxxyz, \
                             tk_xxyzzz_xxxxzz, \
                             tk_xxyzzz_xxxyyy, \
                             tk_xxyzzz_xxxyyz, \
                             tk_xxyzzz_xxxyzz, \
                             tk_xxyzzz_xxxzzz, \
                             tk_xxyzzz_xxyyyy, \
                             tk_xxyzzz_xxyyyz, \
                             tk_xxyzzz_xxyyzz, \
                             tk_xxyzzz_xxyzzz, \
                             tk_xxyzzz_xxzzzz, \
                             tk_xxyzzz_xyyyyy, \
                             tk_xxyzzz_xyyyyz, \
                             tk_xxyzzz_xyyyzz, \
                             tk_xxyzzz_xyyzzz, \
                             tk_xxyzzz_xyzzzz, \
                             tk_xxyzzz_xzzzzz, \
                             tk_xxyzzz_yyyyyy, \
                             tk_xxyzzz_yyyyyz, \
                             tk_xxyzzz_yyyyzz, \
                             tk_xxyzzz_yyyzzz, \
                             tk_xxyzzz_yyzzzz, \
                             tk_xxyzzz_yzzzzz, \
                             tk_xxyzzz_zzzzzz, \
                             tk_xxzzz_xxxxx,   \
                             tk_xxzzz_xxxxxx,  \
                             tk_xxzzz_xxxxxy,  \
                             tk_xxzzz_xxxxxz,  \
                             tk_xxzzz_xxxxy,   \
                             tk_xxzzz_xxxxyy,  \
                             tk_xxzzz_xxxxyz,  \
                             tk_xxzzz_xxxxz,   \
                             tk_xxzzz_xxxxzz,  \
                             tk_xxzzz_xxxyy,   \
                             tk_xxzzz_xxxyyy,  \
                             tk_xxzzz_xxxyyz,  \
                             tk_xxzzz_xxxyz,   \
                             tk_xxzzz_xxxyzz,  \
                             tk_xxzzz_xxxzz,   \
                             tk_xxzzz_xxxzzz,  \
                             tk_xxzzz_xxyyy,   \
                             tk_xxzzz_xxyyyy,  \
                             tk_xxzzz_xxyyyz,  \
                             tk_xxzzz_xxyyz,   \
                             tk_xxzzz_xxyyzz,  \
                             tk_xxzzz_xxyzz,   \
                             tk_xxzzz_xxyzzz,  \
                             tk_xxzzz_xxzzz,   \
                             tk_xxzzz_xxzzzz,  \
                             tk_xxzzz_xyyyy,   \
                             tk_xxzzz_xyyyyy,  \
                             tk_xxzzz_xyyyyz,  \
                             tk_xxzzz_xyyyz,   \
                             tk_xxzzz_xyyyzz,  \
                             tk_xxzzz_xyyzz,   \
                             tk_xxzzz_xyyzzz,  \
                             tk_xxzzz_xyzzz,   \
                             tk_xxzzz_xyzzzz,  \
                             tk_xxzzz_xzzzz,   \
                             tk_xxzzz_xzzzzz,  \
                             tk_xxzzz_yyyyy,   \
                             tk_xxzzz_yyyyyy,  \
                             tk_xxzzz_yyyyyz,  \
                             tk_xxzzz_yyyyz,   \
                             tk_xxzzz_yyyyzz,  \
                             tk_xxzzz_yyyzz,   \
                             tk_xxzzz_yyyzzz,  \
                             tk_xxzzz_yyzzz,   \
                             tk_xxzzz_yyzzzz,  \
                             tk_xxzzz_yzzzz,   \
                             tk_xxzzz_yzzzzz,  \
                             tk_xxzzz_zzzzz,   \
                             tk_xxzzz_zzzzzz,  \
                             ts_xxyzzz_xxxxxx, \
                             ts_xxyzzz_xxxxxy, \
                             ts_xxyzzz_xxxxxz, \
                             ts_xxyzzz_xxxxyy, \
                             ts_xxyzzz_xxxxyz, \
                             ts_xxyzzz_xxxxzz, \
                             ts_xxyzzz_xxxyyy, \
                             ts_xxyzzz_xxxyyz, \
                             ts_xxyzzz_xxxyzz, \
                             ts_xxyzzz_xxxzzz, \
                             ts_xxyzzz_xxyyyy, \
                             ts_xxyzzz_xxyyyz, \
                             ts_xxyzzz_xxyyzz, \
                             ts_xxyzzz_xxyzzz, \
                             ts_xxyzzz_xxzzzz, \
                             ts_xxyzzz_xyyyyy, \
                             ts_xxyzzz_xyyyyz, \
                             ts_xxyzzz_xyyyzz, \
                             ts_xxyzzz_xyyzzz, \
                             ts_xxyzzz_xyzzzz, \
                             ts_xxyzzz_xzzzzz, \
                             ts_xxyzzz_yyyyyy, \
                             ts_xxyzzz_yyyyyz, \
                             ts_xxyzzz_yyyyzz, \
                             ts_xxyzzz_yyyzzz, \
                             ts_xxyzzz_yyzzzz, \
                             ts_xxyzzz_yzzzzz, \
                             ts_xxyzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzzz_xxxxxx[i] = tk_xxzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxxx[i] * fz_0;

        tk_xxyzzz_xxxxxy[i] = tk_xxzzz_xxxxx[i] * fe_0 + tk_xxzzz_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxxy[i] * fz_0;

        tk_xxyzzz_xxxxxz[i] = tk_xxzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxxz[i] * fz_0;

        tk_xxyzzz_xxxxyy[i] = 2.0 * tk_xxzzz_xxxxy[i] * fe_0 + tk_xxzzz_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxyy[i] * fz_0;

        tk_xxyzzz_xxxxyz[i] = tk_xxzzz_xxxxz[i] * fe_0 + tk_xxzzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxyz[i] * fz_0;

        tk_xxyzzz_xxxxzz[i] = tk_xxzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxzz[i] * fz_0;

        tk_xxyzzz_xxxyyy[i] = 3.0 * tk_xxzzz_xxxyy[i] * fe_0 + tk_xxzzz_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxyyy[i] * fz_0;

        tk_xxyzzz_xxxyyz[i] = 2.0 * tk_xxzzz_xxxyz[i] * fe_0 + tk_xxzzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxyyz[i] * fz_0;

        tk_xxyzzz_xxxyzz[i] = tk_xxzzz_xxxzz[i] * fe_0 + tk_xxzzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxyzz[i] * fz_0;

        tk_xxyzzz_xxxzzz[i] = tk_xxzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxzzz[i] * fz_0;

        tk_xxyzzz_xxyyyy[i] = 4.0 * tk_xxzzz_xxyyy[i] * fe_0 + tk_xxzzz_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyyyy[i] * fz_0;

        tk_xxyzzz_xxyyyz[i] = 3.0 * tk_xxzzz_xxyyz[i] * fe_0 + tk_xxzzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyyyz[i] * fz_0;

        tk_xxyzzz_xxyyzz[i] = 2.0 * tk_xxzzz_xxyzz[i] * fe_0 + tk_xxzzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyyzz[i] * fz_0;

        tk_xxyzzz_xxyzzz[i] = tk_xxzzz_xxzzz[i] * fe_0 + tk_xxzzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyzzz[i] * fz_0;

        tk_xxyzzz_xxzzzz[i] = tk_xxzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxzzzz[i] * fz_0;

        tk_xxyzzz_xyyyyy[i] = 5.0 * tk_xxzzz_xyyyy[i] * fe_0 + tk_xxzzz_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyyyy[i] * fz_0;

        tk_xxyzzz_xyyyyz[i] = 4.0 * tk_xxzzz_xyyyz[i] * fe_0 + tk_xxzzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyyyz[i] * fz_0;

        tk_xxyzzz_xyyyzz[i] = 3.0 * tk_xxzzz_xyyzz[i] * fe_0 + tk_xxzzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyyzz[i] * fz_0;

        tk_xxyzzz_xyyzzz[i] = 2.0 * tk_xxzzz_xyzzz[i] * fe_0 + tk_xxzzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyzzz[i] * fz_0;

        tk_xxyzzz_xyzzzz[i] = tk_xxzzz_xzzzz[i] * fe_0 + tk_xxzzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyzzzz[i] * fz_0;

        tk_xxyzzz_xzzzzz[i] = tk_xxzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xzzzzz[i] * fz_0;

        tk_xxyzzz_yyyyyy[i] = 6.0 * tk_xxzzz_yyyyy[i] * fe_0 + tk_xxzzz_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyyyy[i] * fz_0;

        tk_xxyzzz_yyyyyz[i] = 5.0 * tk_xxzzz_yyyyz[i] * fe_0 + tk_xxzzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyyyz[i] * fz_0;

        tk_xxyzzz_yyyyzz[i] = 4.0 * tk_xxzzz_yyyzz[i] * fe_0 + tk_xxzzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyyzz[i] * fz_0;

        tk_xxyzzz_yyyzzz[i] = 3.0 * tk_xxzzz_yyzzz[i] * fe_0 + tk_xxzzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyzzz[i] * fz_0;

        tk_xxyzzz_yyzzzz[i] = 2.0 * tk_xxzzz_yzzzz[i] * fe_0 + tk_xxzzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyzzzz[i] * fz_0;

        tk_xxyzzz_yzzzzz[i] = tk_xxzzz_zzzzz[i] * fe_0 + tk_xxzzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yzzzzz[i] * fz_0;

        tk_xxyzzz_zzzzzz[i] = tk_xxzzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_zzzzzz[i] * fz_0;
    }

    // Set up 392-420 components of targeted buffer : II

    auto tk_xxzzzz_xxxxxx = pbuffer.data(idx_kin_ii + 392);

    auto tk_xxzzzz_xxxxxy = pbuffer.data(idx_kin_ii + 393);

    auto tk_xxzzzz_xxxxxz = pbuffer.data(idx_kin_ii + 394);

    auto tk_xxzzzz_xxxxyy = pbuffer.data(idx_kin_ii + 395);

    auto tk_xxzzzz_xxxxyz = pbuffer.data(idx_kin_ii + 396);

    auto tk_xxzzzz_xxxxzz = pbuffer.data(idx_kin_ii + 397);

    auto tk_xxzzzz_xxxyyy = pbuffer.data(idx_kin_ii + 398);

    auto tk_xxzzzz_xxxyyz = pbuffer.data(idx_kin_ii + 399);

    auto tk_xxzzzz_xxxyzz = pbuffer.data(idx_kin_ii + 400);

    auto tk_xxzzzz_xxxzzz = pbuffer.data(idx_kin_ii + 401);

    auto tk_xxzzzz_xxyyyy = pbuffer.data(idx_kin_ii + 402);

    auto tk_xxzzzz_xxyyyz = pbuffer.data(idx_kin_ii + 403);

    auto tk_xxzzzz_xxyyzz = pbuffer.data(idx_kin_ii + 404);

    auto tk_xxzzzz_xxyzzz = pbuffer.data(idx_kin_ii + 405);

    auto tk_xxzzzz_xxzzzz = pbuffer.data(idx_kin_ii + 406);

    auto tk_xxzzzz_xyyyyy = pbuffer.data(idx_kin_ii + 407);

    auto tk_xxzzzz_xyyyyz = pbuffer.data(idx_kin_ii + 408);

    auto tk_xxzzzz_xyyyzz = pbuffer.data(idx_kin_ii + 409);

    auto tk_xxzzzz_xyyzzz = pbuffer.data(idx_kin_ii + 410);

    auto tk_xxzzzz_xyzzzz = pbuffer.data(idx_kin_ii + 411);

    auto tk_xxzzzz_xzzzzz = pbuffer.data(idx_kin_ii + 412);

    auto tk_xxzzzz_yyyyyy = pbuffer.data(idx_kin_ii + 413);

    auto tk_xxzzzz_yyyyyz = pbuffer.data(idx_kin_ii + 414);

    auto tk_xxzzzz_yyyyzz = pbuffer.data(idx_kin_ii + 415);

    auto tk_xxzzzz_yyyzzz = pbuffer.data(idx_kin_ii + 416);

    auto tk_xxzzzz_yyzzzz = pbuffer.data(idx_kin_ii + 417);

    auto tk_xxzzzz_yzzzzz = pbuffer.data(idx_kin_ii + 418);

    auto tk_xxzzzz_zzzzzz = pbuffer.data(idx_kin_ii + 419);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tk_xxzz_xxxxxx,   \
                             tk_xxzz_xxxxxy,   \
                             tk_xxzz_xxxxyy,   \
                             tk_xxzz_xxxyyy,   \
                             tk_xxzz_xxyyyy,   \
                             tk_xxzz_xyyyyy,   \
                             tk_xxzzz_xxxxxx,  \
                             tk_xxzzz_xxxxxy,  \
                             tk_xxzzz_xxxxyy,  \
                             tk_xxzzz_xxxyyy,  \
                             tk_xxzzz_xxyyyy,  \
                             tk_xxzzz_xyyyyy,  \
                             tk_xxzzzz_xxxxxx, \
                             tk_xxzzzz_xxxxxy, \
                             tk_xxzzzz_xxxxxz, \
                             tk_xxzzzz_xxxxyy, \
                             tk_xxzzzz_xxxxyz, \
                             tk_xxzzzz_xxxxzz, \
                             tk_xxzzzz_xxxyyy, \
                             tk_xxzzzz_xxxyyz, \
                             tk_xxzzzz_xxxyzz, \
                             tk_xxzzzz_xxxzzz, \
                             tk_xxzzzz_xxyyyy, \
                             tk_xxzzzz_xxyyyz, \
                             tk_xxzzzz_xxyyzz, \
                             tk_xxzzzz_xxyzzz, \
                             tk_xxzzzz_xxzzzz, \
                             tk_xxzzzz_xyyyyy, \
                             tk_xxzzzz_xyyyyz, \
                             tk_xxzzzz_xyyyzz, \
                             tk_xxzzzz_xyyzzz, \
                             tk_xxzzzz_xyzzzz, \
                             tk_xxzzzz_xzzzzz, \
                             tk_xxzzzz_yyyyyy, \
                             tk_xxzzzz_yyyyyz, \
                             tk_xxzzzz_yyyyzz, \
                             tk_xxzzzz_yyyzzz, \
                             tk_xxzzzz_yyzzzz, \
                             tk_xxzzzz_yzzzzz, \
                             tk_xxzzzz_zzzzzz, \
                             tk_xzzzz_xxxxxz,  \
                             tk_xzzzz_xxxxyz,  \
                             tk_xzzzz_xxxxz,   \
                             tk_xzzzz_xxxxzz,  \
                             tk_xzzzz_xxxyyz,  \
                             tk_xzzzz_xxxyz,   \
                             tk_xzzzz_xxxyzz,  \
                             tk_xzzzz_xxxzz,   \
                             tk_xzzzz_xxxzzz,  \
                             tk_xzzzz_xxyyyz,  \
                             tk_xzzzz_xxyyz,   \
                             tk_xzzzz_xxyyzz,  \
                             tk_xzzzz_xxyzz,   \
                             tk_xzzzz_xxyzzz,  \
                             tk_xzzzz_xxzzz,   \
                             tk_xzzzz_xxzzzz,  \
                             tk_xzzzz_xyyyyz,  \
                             tk_xzzzz_xyyyz,   \
                             tk_xzzzz_xyyyzz,  \
                             tk_xzzzz_xyyzz,   \
                             tk_xzzzz_xyyzzz,  \
                             tk_xzzzz_xyzzz,   \
                             tk_xzzzz_xyzzzz,  \
                             tk_xzzzz_xzzzz,   \
                             tk_xzzzz_xzzzzz,  \
                             tk_xzzzz_yyyyyy,  \
                             tk_xzzzz_yyyyyz,  \
                             tk_xzzzz_yyyyz,   \
                             tk_xzzzz_yyyyzz,  \
                             tk_xzzzz_yyyzz,   \
                             tk_xzzzz_yyyzzz,  \
                             tk_xzzzz_yyzzz,   \
                             tk_xzzzz_yyzzzz,  \
                             tk_xzzzz_yzzzz,   \
                             tk_xzzzz_yzzzzz,  \
                             tk_xzzzz_zzzzz,   \
                             tk_xzzzz_zzzzzz,  \
                             tk_zzzz_xxxxxz,   \
                             tk_zzzz_xxxxyz,   \
                             tk_zzzz_xxxxzz,   \
                             tk_zzzz_xxxyyz,   \
                             tk_zzzz_xxxyzz,   \
                             tk_zzzz_xxxzzz,   \
                             tk_zzzz_xxyyyz,   \
                             tk_zzzz_xxyyzz,   \
                             tk_zzzz_xxyzzz,   \
                             tk_zzzz_xxzzzz,   \
                             tk_zzzz_xyyyyz,   \
                             tk_zzzz_xyyyzz,   \
                             tk_zzzz_xyyzzz,   \
                             tk_zzzz_xyzzzz,   \
                             tk_zzzz_xzzzzz,   \
                             tk_zzzz_yyyyyy,   \
                             tk_zzzz_yyyyyz,   \
                             tk_zzzz_yyyyzz,   \
                             tk_zzzz_yyyzzz,   \
                             tk_zzzz_yyzzzz,   \
                             tk_zzzz_yzzzzz,   \
                             tk_zzzz_zzzzzz,   \
                             ts_xxzz_xxxxxx,   \
                             ts_xxzz_xxxxxy,   \
                             ts_xxzz_xxxxyy,   \
                             ts_xxzz_xxxyyy,   \
                             ts_xxzz_xxyyyy,   \
                             ts_xxzz_xyyyyy,   \
                             ts_xxzzzz_xxxxxx, \
                             ts_xxzzzz_xxxxxy, \
                             ts_xxzzzz_xxxxxz, \
                             ts_xxzzzz_xxxxyy, \
                             ts_xxzzzz_xxxxyz, \
                             ts_xxzzzz_xxxxzz, \
                             ts_xxzzzz_xxxyyy, \
                             ts_xxzzzz_xxxyyz, \
                             ts_xxzzzz_xxxyzz, \
                             ts_xxzzzz_xxxzzz, \
                             ts_xxzzzz_xxyyyy, \
                             ts_xxzzzz_xxyyyz, \
                             ts_xxzzzz_xxyyzz, \
                             ts_xxzzzz_xxyzzz, \
                             ts_xxzzzz_xxzzzz, \
                             ts_xxzzzz_xyyyyy, \
                             ts_xxzzzz_xyyyyz, \
                             ts_xxzzzz_xyyyzz, \
                             ts_xxzzzz_xyyzzz, \
                             ts_xxzzzz_xyzzzz, \
                             ts_xxzzzz_xzzzzz, \
                             ts_xxzzzz_yyyyyy, \
                             ts_xxzzzz_yyyyyz, \
                             ts_xxzzzz_yyyyzz, \
                             ts_xxzzzz_yyyzzz, \
                             ts_xxzzzz_yyzzzz, \
                             ts_xxzzzz_yzzzzz, \
                             ts_xxzzzz_zzzzzz, \
                             ts_zzzz_xxxxxz,   \
                             ts_zzzz_xxxxyz,   \
                             ts_zzzz_xxxxzz,   \
                             ts_zzzz_xxxyyz,   \
                             ts_zzzz_xxxyzz,   \
                             ts_zzzz_xxxzzz,   \
                             ts_zzzz_xxyyyz,   \
                             ts_zzzz_xxyyzz,   \
                             ts_zzzz_xxyzzz,   \
                             ts_zzzz_xxzzzz,   \
                             ts_zzzz_xyyyyz,   \
                             ts_zzzz_xyyyzz,   \
                             ts_zzzz_xyyzzz,   \
                             ts_zzzz_xyzzzz,   \
                             ts_zzzz_xzzzzz,   \
                             ts_zzzz_yyyyyy,   \
                             ts_zzzz_yyyyyz,   \
                             ts_zzzz_yyyyzz,   \
                             ts_zzzz_yyyzzz,   \
                             ts_zzzz_yyzzzz,   \
                             ts_zzzz_yzzzzz,   \
                             ts_zzzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzzz_xxxxxx[i] = -6.0 * ts_xxzz_xxxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxxx[i] * fe_0 + tk_xxzzz_xxxxxx[i] * pa_z[i] +
                              2.0 * ts_xxzzzz_xxxxxx[i] * fz_0;

        tk_xxzzzz_xxxxxy[i] = -6.0 * ts_xxzz_xxxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxxy[i] * fe_0 + tk_xxzzz_xxxxxy[i] * pa_z[i] +
                              2.0 * ts_xxzzzz_xxxxxy[i] * fz_0;

        tk_xxzzzz_xxxxxz[i] = -2.0 * ts_zzzz_xxxxxz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxxz[i] * fe_0 + 5.0 * tk_xzzzz_xxxxz[i] * fe_0 +
                              tk_xzzzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxxxz[i] * fz_0;

        tk_xxzzzz_xxxxyy[i] = -6.0 * ts_xxzz_xxxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxyy[i] * fe_0 + tk_xxzzz_xxxxyy[i] * pa_z[i] +
                              2.0 * ts_xxzzzz_xxxxyy[i] * fz_0;

        tk_xxzzzz_xxxxyz[i] = -2.0 * ts_zzzz_xxxxyz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxyz[i] * fe_0 + 4.0 * tk_xzzzz_xxxyz[i] * fe_0 +
                              tk_xzzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxxyz[i] * fz_0;

        tk_xxzzzz_xxxxzz[i] = -2.0 * ts_zzzz_xxxxzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxzz[i] * fe_0 + 4.0 * tk_xzzzz_xxxzz[i] * fe_0 +
                              tk_xzzzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxxzz[i] * fz_0;

        tk_xxzzzz_xxxyyy[i] = -6.0 * ts_xxzz_xxxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxyyy[i] * fe_0 + tk_xxzzz_xxxyyy[i] * pa_z[i] +
                              2.0 * ts_xxzzzz_xxxyyy[i] * fz_0;

        tk_xxzzzz_xxxyyz[i] = -2.0 * ts_zzzz_xxxyyz[i] * fbe_0 * fz_0 + tk_zzzz_xxxyyz[i] * fe_0 + 3.0 * tk_xzzzz_xxyyz[i] * fe_0 +
                              tk_xzzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxyyz[i] * fz_0;

        tk_xxzzzz_xxxyzz[i] = -2.0 * ts_zzzz_xxxyzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxyzz[i] * fe_0 + 3.0 * tk_xzzzz_xxyzz[i] * fe_0 +
                              tk_xzzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxyzz[i] * fz_0;

        tk_xxzzzz_xxxzzz[i] = -2.0 * ts_zzzz_xxxzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxzzz[i] * fe_0 + 3.0 * tk_xzzzz_xxzzz[i] * fe_0 +
                              tk_xzzzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxzzz[i] * fz_0;

        tk_xxzzzz_xxyyyy[i] = -6.0 * ts_xxzz_xxyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyyyy[i] * fe_0 + tk_xxzzz_xxyyyy[i] * pa_z[i] +
                              2.0 * ts_xxzzzz_xxyyyy[i] * fz_0;

        tk_xxzzzz_xxyyyz[i] = -2.0 * ts_zzzz_xxyyyz[i] * fbe_0 * fz_0 + tk_zzzz_xxyyyz[i] * fe_0 + 2.0 * tk_xzzzz_xyyyz[i] * fe_0 +
                              tk_xzzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxyyyz[i] * fz_0;

        tk_xxzzzz_xxyyzz[i] = -2.0 * ts_zzzz_xxyyzz[i] * fbe_0 * fz_0 + tk_zzzz_xxyyzz[i] * fe_0 + 2.0 * tk_xzzzz_xyyzz[i] * fe_0 +
                              tk_xzzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxyyzz[i] * fz_0;

        tk_xxzzzz_xxyzzz[i] = -2.0 * ts_zzzz_xxyzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxyzzz[i] * fe_0 + 2.0 * tk_xzzzz_xyzzz[i] * fe_0 +
                              tk_xzzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxyzzz[i] * fz_0;

        tk_xxzzzz_xxzzzz[i] = -2.0 * ts_zzzz_xxzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxzzzz[i] * fe_0 + 2.0 * tk_xzzzz_xzzzz[i] * fe_0 +
                              tk_xzzzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxzzzz[i] * fz_0;

        tk_xxzzzz_xyyyyy[i] = -6.0 * ts_xxzz_xyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyyyy[i] * fe_0 + tk_xxzzz_xyyyyy[i] * pa_z[i] +
                              2.0 * ts_xxzzzz_xyyyyy[i] * fz_0;

        tk_xxzzzz_xyyyyz[i] = -2.0 * ts_zzzz_xyyyyz[i] * fbe_0 * fz_0 + tk_zzzz_xyyyyz[i] * fe_0 + tk_xzzzz_yyyyz[i] * fe_0 +
                              tk_xzzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyyyyz[i] * fz_0;

        tk_xxzzzz_xyyyzz[i] = -2.0 * ts_zzzz_xyyyzz[i] * fbe_0 * fz_0 + tk_zzzz_xyyyzz[i] * fe_0 + tk_xzzzz_yyyzz[i] * fe_0 +
                              tk_xzzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyyyzz[i] * fz_0;

        tk_xxzzzz_xyyzzz[i] = -2.0 * ts_zzzz_xyyzzz[i] * fbe_0 * fz_0 + tk_zzzz_xyyzzz[i] * fe_0 + tk_xzzzz_yyzzz[i] * fe_0 +
                              tk_xzzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyyzzz[i] * fz_0;

        tk_xxzzzz_xyzzzz[i] = -2.0 * ts_zzzz_xyzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xyzzzz[i] * fe_0 + tk_xzzzz_yzzzz[i] * fe_0 +
                              tk_xzzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyzzzz[i] * fz_0;

        tk_xxzzzz_xzzzzz[i] = -2.0 * ts_zzzz_xzzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xzzzzz[i] * fe_0 + tk_xzzzz_zzzzz[i] * fe_0 +
                              tk_xzzzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xzzzzz[i] * fz_0;

        tk_xxzzzz_yyyyyy[i] =
            -2.0 * ts_zzzz_yyyyyy[i] * fbe_0 * fz_0 + tk_zzzz_yyyyyy[i] * fe_0 + tk_xzzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyyyy[i] * fz_0;

        tk_xxzzzz_yyyyyz[i] =
            -2.0 * ts_zzzz_yyyyyz[i] * fbe_0 * fz_0 + tk_zzzz_yyyyyz[i] * fe_0 + tk_xzzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyyyz[i] * fz_0;

        tk_xxzzzz_yyyyzz[i] =
            -2.0 * ts_zzzz_yyyyzz[i] * fbe_0 * fz_0 + tk_zzzz_yyyyzz[i] * fe_0 + tk_xzzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyyzz[i] * fz_0;

        tk_xxzzzz_yyyzzz[i] =
            -2.0 * ts_zzzz_yyyzzz[i] * fbe_0 * fz_0 + tk_zzzz_yyyzzz[i] * fe_0 + tk_xzzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyzzz[i] * fz_0;

        tk_xxzzzz_yyzzzz[i] =
            -2.0 * ts_zzzz_yyzzzz[i] * fbe_0 * fz_0 + tk_zzzz_yyzzzz[i] * fe_0 + tk_xzzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyzzzz[i] * fz_0;

        tk_xxzzzz_yzzzzz[i] =
            -2.0 * ts_zzzz_yzzzzz[i] * fbe_0 * fz_0 + tk_zzzz_yzzzzz[i] * fe_0 + tk_xzzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yzzzzz[i] * fz_0;

        tk_xxzzzz_zzzzzz[i] =
            -2.0 * ts_zzzz_zzzzzz[i] * fbe_0 * fz_0 + tk_zzzz_zzzzzz[i] * fe_0 + tk_xzzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 420-448 components of targeted buffer : II

    auto tk_xyyyyy_xxxxxx = pbuffer.data(idx_kin_ii + 420);

    auto tk_xyyyyy_xxxxxy = pbuffer.data(idx_kin_ii + 421);

    auto tk_xyyyyy_xxxxxz = pbuffer.data(idx_kin_ii + 422);

    auto tk_xyyyyy_xxxxyy = pbuffer.data(idx_kin_ii + 423);

    auto tk_xyyyyy_xxxxyz = pbuffer.data(idx_kin_ii + 424);

    auto tk_xyyyyy_xxxxzz = pbuffer.data(idx_kin_ii + 425);

    auto tk_xyyyyy_xxxyyy = pbuffer.data(idx_kin_ii + 426);

    auto tk_xyyyyy_xxxyyz = pbuffer.data(idx_kin_ii + 427);

    auto tk_xyyyyy_xxxyzz = pbuffer.data(idx_kin_ii + 428);

    auto tk_xyyyyy_xxxzzz = pbuffer.data(idx_kin_ii + 429);

    auto tk_xyyyyy_xxyyyy = pbuffer.data(idx_kin_ii + 430);

    auto tk_xyyyyy_xxyyyz = pbuffer.data(idx_kin_ii + 431);

    auto tk_xyyyyy_xxyyzz = pbuffer.data(idx_kin_ii + 432);

    auto tk_xyyyyy_xxyzzz = pbuffer.data(idx_kin_ii + 433);

    auto tk_xyyyyy_xxzzzz = pbuffer.data(idx_kin_ii + 434);

    auto tk_xyyyyy_xyyyyy = pbuffer.data(idx_kin_ii + 435);

    auto tk_xyyyyy_xyyyyz = pbuffer.data(idx_kin_ii + 436);

    auto tk_xyyyyy_xyyyzz = pbuffer.data(idx_kin_ii + 437);

    auto tk_xyyyyy_xyyzzz = pbuffer.data(idx_kin_ii + 438);

    auto tk_xyyyyy_xyzzzz = pbuffer.data(idx_kin_ii + 439);

    auto tk_xyyyyy_xzzzzz = pbuffer.data(idx_kin_ii + 440);

    auto tk_xyyyyy_yyyyyy = pbuffer.data(idx_kin_ii + 441);

    auto tk_xyyyyy_yyyyyz = pbuffer.data(idx_kin_ii + 442);

    auto tk_xyyyyy_yyyyzz = pbuffer.data(idx_kin_ii + 443);

    auto tk_xyyyyy_yyyzzz = pbuffer.data(idx_kin_ii + 444);

    auto tk_xyyyyy_yyzzzz = pbuffer.data(idx_kin_ii + 445);

    auto tk_xyyyyy_yzzzzz = pbuffer.data(idx_kin_ii + 446);

    auto tk_xyyyyy_zzzzzz = pbuffer.data(idx_kin_ii + 447);

#pragma omp simd aligned(pa_x,                 \
                             tk_xyyyyy_xxxxxx, \
                             tk_xyyyyy_xxxxxy, \
                             tk_xyyyyy_xxxxxz, \
                             tk_xyyyyy_xxxxyy, \
                             tk_xyyyyy_xxxxyz, \
                             tk_xyyyyy_xxxxzz, \
                             tk_xyyyyy_xxxyyy, \
                             tk_xyyyyy_xxxyyz, \
                             tk_xyyyyy_xxxyzz, \
                             tk_xyyyyy_xxxzzz, \
                             tk_xyyyyy_xxyyyy, \
                             tk_xyyyyy_xxyyyz, \
                             tk_xyyyyy_xxyyzz, \
                             tk_xyyyyy_xxyzzz, \
                             tk_xyyyyy_xxzzzz, \
                             tk_xyyyyy_xyyyyy, \
                             tk_xyyyyy_xyyyyz, \
                             tk_xyyyyy_xyyyzz, \
                             tk_xyyyyy_xyyzzz, \
                             tk_xyyyyy_xyzzzz, \
                             tk_xyyyyy_xzzzzz, \
                             tk_xyyyyy_yyyyyy, \
                             tk_xyyyyy_yyyyyz, \
                             tk_xyyyyy_yyyyzz, \
                             tk_xyyyyy_yyyzzz, \
                             tk_xyyyyy_yyzzzz, \
                             tk_xyyyyy_yzzzzz, \
                             tk_xyyyyy_zzzzzz, \
                             tk_yyyyy_xxxxx,   \
                             tk_yyyyy_xxxxxx,  \
                             tk_yyyyy_xxxxxy,  \
                             tk_yyyyy_xxxxxz,  \
                             tk_yyyyy_xxxxy,   \
                             tk_yyyyy_xxxxyy,  \
                             tk_yyyyy_xxxxyz,  \
                             tk_yyyyy_xxxxz,   \
                             tk_yyyyy_xxxxzz,  \
                             tk_yyyyy_xxxyy,   \
                             tk_yyyyy_xxxyyy,  \
                             tk_yyyyy_xxxyyz,  \
                             tk_yyyyy_xxxyz,   \
                             tk_yyyyy_xxxyzz,  \
                             tk_yyyyy_xxxzz,   \
                             tk_yyyyy_xxxzzz,  \
                             tk_yyyyy_xxyyy,   \
                             tk_yyyyy_xxyyyy,  \
                             tk_yyyyy_xxyyyz,  \
                             tk_yyyyy_xxyyz,   \
                             tk_yyyyy_xxyyzz,  \
                             tk_yyyyy_xxyzz,   \
                             tk_yyyyy_xxyzzz,  \
                             tk_yyyyy_xxzzz,   \
                             tk_yyyyy_xxzzzz,  \
                             tk_yyyyy_xyyyy,   \
                             tk_yyyyy_xyyyyy,  \
                             tk_yyyyy_xyyyyz,  \
                             tk_yyyyy_xyyyz,   \
                             tk_yyyyy_xyyyzz,  \
                             tk_yyyyy_xyyzz,   \
                             tk_yyyyy_xyyzzz,  \
                             tk_yyyyy_xyzzz,   \
                             tk_yyyyy_xyzzzz,  \
                             tk_yyyyy_xzzzz,   \
                             tk_yyyyy_xzzzzz,  \
                             tk_yyyyy_yyyyy,   \
                             tk_yyyyy_yyyyyy,  \
                             tk_yyyyy_yyyyyz,  \
                             tk_yyyyy_yyyyz,   \
                             tk_yyyyy_yyyyzz,  \
                             tk_yyyyy_yyyzz,   \
                             tk_yyyyy_yyyzzz,  \
                             tk_yyyyy_yyzzz,   \
                             tk_yyyyy_yyzzzz,  \
                             tk_yyyyy_yzzzz,   \
                             tk_yyyyy_yzzzzz,  \
                             tk_yyyyy_zzzzz,   \
                             tk_yyyyy_zzzzzz,  \
                             ts_xyyyyy_xxxxxx, \
                             ts_xyyyyy_xxxxxy, \
                             ts_xyyyyy_xxxxxz, \
                             ts_xyyyyy_xxxxyy, \
                             ts_xyyyyy_xxxxyz, \
                             ts_xyyyyy_xxxxzz, \
                             ts_xyyyyy_xxxyyy, \
                             ts_xyyyyy_xxxyyz, \
                             ts_xyyyyy_xxxyzz, \
                             ts_xyyyyy_xxxzzz, \
                             ts_xyyyyy_xxyyyy, \
                             ts_xyyyyy_xxyyyz, \
                             ts_xyyyyy_xxyyzz, \
                             ts_xyyyyy_xxyzzz, \
                             ts_xyyyyy_xxzzzz, \
                             ts_xyyyyy_xyyyyy, \
                             ts_xyyyyy_xyyyyz, \
                             ts_xyyyyy_xyyyzz, \
                             ts_xyyyyy_xyyzzz, \
                             ts_xyyyyy_xyzzzz, \
                             ts_xyyyyy_xzzzzz, \
                             ts_xyyyyy_yyyyyy, \
                             ts_xyyyyy_yyyyyz, \
                             ts_xyyyyy_yyyyzz, \
                             ts_xyyyyy_yyyzzz, \
                             ts_xyyyyy_yyzzzz, \
                             ts_xyyyyy_yzzzzz, \
                             ts_xyyyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyy_xxxxxx[i] = 6.0 * tk_yyyyy_xxxxx[i] * fe_0 + tk_yyyyy_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxxx[i] * fz_0;

        tk_xyyyyy_xxxxxy[i] = 5.0 * tk_yyyyy_xxxxy[i] * fe_0 + tk_yyyyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxxy[i] * fz_0;

        tk_xyyyyy_xxxxxz[i] = 5.0 * tk_yyyyy_xxxxz[i] * fe_0 + tk_yyyyy_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxxz[i] * fz_0;

        tk_xyyyyy_xxxxyy[i] = 4.0 * tk_yyyyy_xxxyy[i] * fe_0 + tk_yyyyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxyy[i] * fz_0;

        tk_xyyyyy_xxxxyz[i] = 4.0 * tk_yyyyy_xxxyz[i] * fe_0 + tk_yyyyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxyz[i] * fz_0;

        tk_xyyyyy_xxxxzz[i] = 4.0 * tk_yyyyy_xxxzz[i] * fe_0 + tk_yyyyy_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxzz[i] * fz_0;

        tk_xyyyyy_xxxyyy[i] = 3.0 * tk_yyyyy_xxyyy[i] * fe_0 + tk_yyyyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxyyy[i] * fz_0;

        tk_xyyyyy_xxxyyz[i] = 3.0 * tk_yyyyy_xxyyz[i] * fe_0 + tk_yyyyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxyyz[i] * fz_0;

        tk_xyyyyy_xxxyzz[i] = 3.0 * tk_yyyyy_xxyzz[i] * fe_0 + tk_yyyyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxyzz[i] * fz_0;

        tk_xyyyyy_xxxzzz[i] = 3.0 * tk_yyyyy_xxzzz[i] * fe_0 + tk_yyyyy_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxzzz[i] * fz_0;

        tk_xyyyyy_xxyyyy[i] = 2.0 * tk_yyyyy_xyyyy[i] * fe_0 + tk_yyyyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyyyy[i] * fz_0;

        tk_xyyyyy_xxyyyz[i] = 2.0 * tk_yyyyy_xyyyz[i] * fe_0 + tk_yyyyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyyyz[i] * fz_0;

        tk_xyyyyy_xxyyzz[i] = 2.0 * tk_yyyyy_xyyzz[i] * fe_0 + tk_yyyyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyyzz[i] * fz_0;

        tk_xyyyyy_xxyzzz[i] = 2.0 * tk_yyyyy_xyzzz[i] * fe_0 + tk_yyyyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyzzz[i] * fz_0;

        tk_xyyyyy_xxzzzz[i] = 2.0 * tk_yyyyy_xzzzz[i] * fe_0 + tk_yyyyy_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxzzzz[i] * fz_0;

        tk_xyyyyy_xyyyyy[i] = tk_yyyyy_yyyyy[i] * fe_0 + tk_yyyyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyyyy[i] * fz_0;

        tk_xyyyyy_xyyyyz[i] = tk_yyyyy_yyyyz[i] * fe_0 + tk_yyyyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyyyz[i] * fz_0;

        tk_xyyyyy_xyyyzz[i] = tk_yyyyy_yyyzz[i] * fe_0 + tk_yyyyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyyzz[i] * fz_0;

        tk_xyyyyy_xyyzzz[i] = tk_yyyyy_yyzzz[i] * fe_0 + tk_yyyyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyzzz[i] * fz_0;

        tk_xyyyyy_xyzzzz[i] = tk_yyyyy_yzzzz[i] * fe_0 + tk_yyyyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyzzzz[i] * fz_0;

        tk_xyyyyy_xzzzzz[i] = tk_yyyyy_zzzzz[i] * fe_0 + tk_yyyyy_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xzzzzz[i] * fz_0;

        tk_xyyyyy_yyyyyy[i] = tk_yyyyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyyyy[i] * fz_0;

        tk_xyyyyy_yyyyyz[i] = tk_yyyyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyyyz[i] * fz_0;

        tk_xyyyyy_yyyyzz[i] = tk_yyyyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyyzz[i] * fz_0;

        tk_xyyyyy_yyyzzz[i] = tk_yyyyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyzzz[i] * fz_0;

        tk_xyyyyy_yyzzzz[i] = tk_yyyyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyzzzz[i] * fz_0;

        tk_xyyyyy_yzzzzz[i] = tk_yyyyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yzzzzz[i] * fz_0;

        tk_xyyyyy_zzzzzz[i] = tk_yyyyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_zzzzzz[i] * fz_0;
    }

    // Set up 448-476 components of targeted buffer : II

    auto tk_xyyyyz_xxxxxx = pbuffer.data(idx_kin_ii + 448);

    auto tk_xyyyyz_xxxxxy = pbuffer.data(idx_kin_ii + 449);

    auto tk_xyyyyz_xxxxxz = pbuffer.data(idx_kin_ii + 450);

    auto tk_xyyyyz_xxxxyy = pbuffer.data(idx_kin_ii + 451);

    auto tk_xyyyyz_xxxxyz = pbuffer.data(idx_kin_ii + 452);

    auto tk_xyyyyz_xxxxzz = pbuffer.data(idx_kin_ii + 453);

    auto tk_xyyyyz_xxxyyy = pbuffer.data(idx_kin_ii + 454);

    auto tk_xyyyyz_xxxyyz = pbuffer.data(idx_kin_ii + 455);

    auto tk_xyyyyz_xxxyzz = pbuffer.data(idx_kin_ii + 456);

    auto tk_xyyyyz_xxxzzz = pbuffer.data(idx_kin_ii + 457);

    auto tk_xyyyyz_xxyyyy = pbuffer.data(idx_kin_ii + 458);

    auto tk_xyyyyz_xxyyyz = pbuffer.data(idx_kin_ii + 459);

    auto tk_xyyyyz_xxyyzz = pbuffer.data(idx_kin_ii + 460);

    auto tk_xyyyyz_xxyzzz = pbuffer.data(idx_kin_ii + 461);

    auto tk_xyyyyz_xxzzzz = pbuffer.data(idx_kin_ii + 462);

    auto tk_xyyyyz_xyyyyy = pbuffer.data(idx_kin_ii + 463);

    auto tk_xyyyyz_xyyyyz = pbuffer.data(idx_kin_ii + 464);

    auto tk_xyyyyz_xyyyzz = pbuffer.data(idx_kin_ii + 465);

    auto tk_xyyyyz_xyyzzz = pbuffer.data(idx_kin_ii + 466);

    auto tk_xyyyyz_xyzzzz = pbuffer.data(idx_kin_ii + 467);

    auto tk_xyyyyz_xzzzzz = pbuffer.data(idx_kin_ii + 468);

    auto tk_xyyyyz_yyyyyy = pbuffer.data(idx_kin_ii + 469);

    auto tk_xyyyyz_yyyyyz = pbuffer.data(idx_kin_ii + 470);

    auto tk_xyyyyz_yyyyzz = pbuffer.data(idx_kin_ii + 471);

    auto tk_xyyyyz_yyyzzz = pbuffer.data(idx_kin_ii + 472);

    auto tk_xyyyyz_yyzzzz = pbuffer.data(idx_kin_ii + 473);

    auto tk_xyyyyz_yzzzzz = pbuffer.data(idx_kin_ii + 474);

    auto tk_xyyyyz_zzzzzz = pbuffer.data(idx_kin_ii + 475);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tk_xyyyy_xxxxxx,  \
                             tk_xyyyy_xxxxxy,  \
                             tk_xyyyy_xxxxyy,  \
                             tk_xyyyy_xxxyyy,  \
                             tk_xyyyy_xxyyyy,  \
                             tk_xyyyy_xyyyyy,  \
                             tk_xyyyyz_xxxxxx, \
                             tk_xyyyyz_xxxxxy, \
                             tk_xyyyyz_xxxxxz, \
                             tk_xyyyyz_xxxxyy, \
                             tk_xyyyyz_xxxxyz, \
                             tk_xyyyyz_xxxxzz, \
                             tk_xyyyyz_xxxyyy, \
                             tk_xyyyyz_xxxyyz, \
                             tk_xyyyyz_xxxyzz, \
                             tk_xyyyyz_xxxzzz, \
                             tk_xyyyyz_xxyyyy, \
                             tk_xyyyyz_xxyyyz, \
                             tk_xyyyyz_xxyyzz, \
                             tk_xyyyyz_xxyzzz, \
                             tk_xyyyyz_xxzzzz, \
                             tk_xyyyyz_xyyyyy, \
                             tk_xyyyyz_xyyyyz, \
                             tk_xyyyyz_xyyyzz, \
                             tk_xyyyyz_xyyzzz, \
                             tk_xyyyyz_xyzzzz, \
                             tk_xyyyyz_xzzzzz, \
                             tk_xyyyyz_yyyyyy, \
                             tk_xyyyyz_yyyyyz, \
                             tk_xyyyyz_yyyyzz, \
                             tk_xyyyyz_yyyzzz, \
                             tk_xyyyyz_yyzzzz, \
                             tk_xyyyyz_yzzzzz, \
                             tk_xyyyyz_zzzzzz, \
                             tk_yyyyz_xxxxxz,  \
                             tk_yyyyz_xxxxyz,  \
                             tk_yyyyz_xxxxz,   \
                             tk_yyyyz_xxxxzz,  \
                             tk_yyyyz_xxxyyz,  \
                             tk_yyyyz_xxxyz,   \
                             tk_yyyyz_xxxyzz,  \
                             tk_yyyyz_xxxzz,   \
                             tk_yyyyz_xxxzzz,  \
                             tk_yyyyz_xxyyyz,  \
                             tk_yyyyz_xxyyz,   \
                             tk_yyyyz_xxyyzz,  \
                             tk_yyyyz_xxyzz,   \
                             tk_yyyyz_xxyzzz,  \
                             tk_yyyyz_xxzzz,   \
                             tk_yyyyz_xxzzzz,  \
                             tk_yyyyz_xyyyyz,  \
                             tk_yyyyz_xyyyz,   \
                             tk_yyyyz_xyyyzz,  \
                             tk_yyyyz_xyyzz,   \
                             tk_yyyyz_xyyzzz,  \
                             tk_yyyyz_xyzzz,   \
                             tk_yyyyz_xyzzzz,  \
                             tk_yyyyz_xzzzz,   \
                             tk_yyyyz_xzzzzz,  \
                             tk_yyyyz_yyyyyy,  \
                             tk_yyyyz_yyyyyz,  \
                             tk_yyyyz_yyyyz,   \
                             tk_yyyyz_yyyyzz,  \
                             tk_yyyyz_yyyzz,   \
                             tk_yyyyz_yyyzzz,  \
                             tk_yyyyz_yyzzz,   \
                             tk_yyyyz_yyzzzz,  \
                             tk_yyyyz_yzzzz,   \
                             tk_yyyyz_yzzzzz,  \
                             tk_yyyyz_zzzzz,   \
                             tk_yyyyz_zzzzzz,  \
                             ts_xyyyyz_xxxxxx, \
                             ts_xyyyyz_xxxxxy, \
                             ts_xyyyyz_xxxxxz, \
                             ts_xyyyyz_xxxxyy, \
                             ts_xyyyyz_xxxxyz, \
                             ts_xyyyyz_xxxxzz, \
                             ts_xyyyyz_xxxyyy, \
                             ts_xyyyyz_xxxyyz, \
                             ts_xyyyyz_xxxyzz, \
                             ts_xyyyyz_xxxzzz, \
                             ts_xyyyyz_xxyyyy, \
                             ts_xyyyyz_xxyyyz, \
                             ts_xyyyyz_xxyyzz, \
                             ts_xyyyyz_xxyzzz, \
                             ts_xyyyyz_xxzzzz, \
                             ts_xyyyyz_xyyyyy, \
                             ts_xyyyyz_xyyyyz, \
                             ts_xyyyyz_xyyyzz, \
                             ts_xyyyyz_xyyzzz, \
                             ts_xyyyyz_xyzzzz, \
                             ts_xyyyyz_xzzzzz, \
                             ts_xyyyyz_yyyyyy, \
                             ts_xyyyyz_yyyyyz, \
                             ts_xyyyyz_yyyyzz, \
                             ts_xyyyyz_yyyzzz, \
                             ts_xyyyyz_yyzzzz, \
                             ts_xyyyyz_yzzzzz, \
                             ts_xyyyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyz_xxxxxx[i] = tk_xyyyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxxxx[i] * fz_0;

        tk_xyyyyz_xxxxxy[i] = tk_xyyyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxxxy[i] * fz_0;

        tk_xyyyyz_xxxxxz[i] = 5.0 * tk_yyyyz_xxxxz[i] * fe_0 + tk_yyyyz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxxxz[i] * fz_0;

        tk_xyyyyz_xxxxyy[i] = tk_xyyyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxxyy[i] * fz_0;

        tk_xyyyyz_xxxxyz[i] = 4.0 * tk_yyyyz_xxxyz[i] * fe_0 + tk_yyyyz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxxyz[i] * fz_0;

        tk_xyyyyz_xxxxzz[i] = 4.0 * tk_yyyyz_xxxzz[i] * fe_0 + tk_yyyyz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxxzz[i] * fz_0;

        tk_xyyyyz_xxxyyy[i] = tk_xyyyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxyyy[i] * fz_0;

        tk_xyyyyz_xxxyyz[i] = 3.0 * tk_yyyyz_xxyyz[i] * fe_0 + tk_yyyyz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxyyz[i] * fz_0;

        tk_xyyyyz_xxxyzz[i] = 3.0 * tk_yyyyz_xxyzz[i] * fe_0 + tk_yyyyz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxyzz[i] * fz_0;

        tk_xyyyyz_xxxzzz[i] = 3.0 * tk_yyyyz_xxzzz[i] * fe_0 + tk_yyyyz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxzzz[i] * fz_0;

        tk_xyyyyz_xxyyyy[i] = tk_xyyyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxyyyy[i] * fz_0;

        tk_xyyyyz_xxyyyz[i] = 2.0 * tk_yyyyz_xyyyz[i] * fe_0 + tk_yyyyz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxyyyz[i] * fz_0;

        tk_xyyyyz_xxyyzz[i] = 2.0 * tk_yyyyz_xyyzz[i] * fe_0 + tk_yyyyz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxyyzz[i] * fz_0;

        tk_xyyyyz_xxyzzz[i] = 2.0 * tk_yyyyz_xyzzz[i] * fe_0 + tk_yyyyz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxyzzz[i] * fz_0;

        tk_xyyyyz_xxzzzz[i] = 2.0 * tk_yyyyz_xzzzz[i] * fe_0 + tk_yyyyz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxzzzz[i] * fz_0;

        tk_xyyyyz_xyyyyy[i] = tk_xyyyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xyyyyy[i] * fz_0;

        tk_xyyyyz_xyyyyz[i] = tk_yyyyz_yyyyz[i] * fe_0 + tk_yyyyz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyyyyz[i] * fz_0;

        tk_xyyyyz_xyyyzz[i] = tk_yyyyz_yyyzz[i] * fe_0 + tk_yyyyz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyyyzz[i] * fz_0;

        tk_xyyyyz_xyyzzz[i] = tk_yyyyz_yyzzz[i] * fe_0 + tk_yyyyz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyyzzz[i] * fz_0;

        tk_xyyyyz_xyzzzz[i] = tk_yyyyz_yzzzz[i] * fe_0 + tk_yyyyz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyzzzz[i] * fz_0;

        tk_xyyyyz_xzzzzz[i] = tk_yyyyz_zzzzz[i] * fe_0 + tk_yyyyz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xzzzzz[i] * fz_0;

        tk_xyyyyz_yyyyyy[i] = tk_yyyyz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyyyy[i] * fz_0;

        tk_xyyyyz_yyyyyz[i] = tk_yyyyz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyyyz[i] * fz_0;

        tk_xyyyyz_yyyyzz[i] = tk_yyyyz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyyzz[i] * fz_0;

        tk_xyyyyz_yyyzzz[i] = tk_yyyyz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyzzz[i] * fz_0;

        tk_xyyyyz_yyzzzz[i] = tk_yyyyz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyzzzz[i] * fz_0;

        tk_xyyyyz_yzzzzz[i] = tk_yyyyz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yzzzzz[i] * fz_0;

        tk_xyyyyz_zzzzzz[i] = tk_yyyyz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_zzzzzz[i] * fz_0;
    }

    // Set up 476-504 components of targeted buffer : II

    auto tk_xyyyzz_xxxxxx = pbuffer.data(idx_kin_ii + 476);

    auto tk_xyyyzz_xxxxxy = pbuffer.data(idx_kin_ii + 477);

    auto tk_xyyyzz_xxxxxz = pbuffer.data(idx_kin_ii + 478);

    auto tk_xyyyzz_xxxxyy = pbuffer.data(idx_kin_ii + 479);

    auto tk_xyyyzz_xxxxyz = pbuffer.data(idx_kin_ii + 480);

    auto tk_xyyyzz_xxxxzz = pbuffer.data(idx_kin_ii + 481);

    auto tk_xyyyzz_xxxyyy = pbuffer.data(idx_kin_ii + 482);

    auto tk_xyyyzz_xxxyyz = pbuffer.data(idx_kin_ii + 483);

    auto tk_xyyyzz_xxxyzz = pbuffer.data(idx_kin_ii + 484);

    auto tk_xyyyzz_xxxzzz = pbuffer.data(idx_kin_ii + 485);

    auto tk_xyyyzz_xxyyyy = pbuffer.data(idx_kin_ii + 486);

    auto tk_xyyyzz_xxyyyz = pbuffer.data(idx_kin_ii + 487);

    auto tk_xyyyzz_xxyyzz = pbuffer.data(idx_kin_ii + 488);

    auto tk_xyyyzz_xxyzzz = pbuffer.data(idx_kin_ii + 489);

    auto tk_xyyyzz_xxzzzz = pbuffer.data(idx_kin_ii + 490);

    auto tk_xyyyzz_xyyyyy = pbuffer.data(idx_kin_ii + 491);

    auto tk_xyyyzz_xyyyyz = pbuffer.data(idx_kin_ii + 492);

    auto tk_xyyyzz_xyyyzz = pbuffer.data(idx_kin_ii + 493);

    auto tk_xyyyzz_xyyzzz = pbuffer.data(idx_kin_ii + 494);

    auto tk_xyyyzz_xyzzzz = pbuffer.data(idx_kin_ii + 495);

    auto tk_xyyyzz_xzzzzz = pbuffer.data(idx_kin_ii + 496);

    auto tk_xyyyzz_yyyyyy = pbuffer.data(idx_kin_ii + 497);

    auto tk_xyyyzz_yyyyyz = pbuffer.data(idx_kin_ii + 498);

    auto tk_xyyyzz_yyyyzz = pbuffer.data(idx_kin_ii + 499);

    auto tk_xyyyzz_yyyzzz = pbuffer.data(idx_kin_ii + 500);

    auto tk_xyyyzz_yyzzzz = pbuffer.data(idx_kin_ii + 501);

    auto tk_xyyyzz_yzzzzz = pbuffer.data(idx_kin_ii + 502);

    auto tk_xyyyzz_zzzzzz = pbuffer.data(idx_kin_ii + 503);

#pragma omp simd aligned(pa_x,                 \
                             tk_xyyyzz_xxxxxx, \
                             tk_xyyyzz_xxxxxy, \
                             tk_xyyyzz_xxxxxz, \
                             tk_xyyyzz_xxxxyy, \
                             tk_xyyyzz_xxxxyz, \
                             tk_xyyyzz_xxxxzz, \
                             tk_xyyyzz_xxxyyy, \
                             tk_xyyyzz_xxxyyz, \
                             tk_xyyyzz_xxxyzz, \
                             tk_xyyyzz_xxxzzz, \
                             tk_xyyyzz_xxyyyy, \
                             tk_xyyyzz_xxyyyz, \
                             tk_xyyyzz_xxyyzz, \
                             tk_xyyyzz_xxyzzz, \
                             tk_xyyyzz_xxzzzz, \
                             tk_xyyyzz_xyyyyy, \
                             tk_xyyyzz_xyyyyz, \
                             tk_xyyyzz_xyyyzz, \
                             tk_xyyyzz_xyyzzz, \
                             tk_xyyyzz_xyzzzz, \
                             tk_xyyyzz_xzzzzz, \
                             tk_xyyyzz_yyyyyy, \
                             tk_xyyyzz_yyyyyz, \
                             tk_xyyyzz_yyyyzz, \
                             tk_xyyyzz_yyyzzz, \
                             tk_xyyyzz_yyzzzz, \
                             tk_xyyyzz_yzzzzz, \
                             tk_xyyyzz_zzzzzz, \
                             tk_yyyzz_xxxxx,   \
                             tk_yyyzz_xxxxxx,  \
                             tk_yyyzz_xxxxxy,  \
                             tk_yyyzz_xxxxxz,  \
                             tk_yyyzz_xxxxy,   \
                             tk_yyyzz_xxxxyy,  \
                             tk_yyyzz_xxxxyz,  \
                             tk_yyyzz_xxxxz,   \
                             tk_yyyzz_xxxxzz,  \
                             tk_yyyzz_xxxyy,   \
                             tk_yyyzz_xxxyyy,  \
                             tk_yyyzz_xxxyyz,  \
                             tk_yyyzz_xxxyz,   \
                             tk_yyyzz_xxxyzz,  \
                             tk_yyyzz_xxxzz,   \
                             tk_yyyzz_xxxzzz,  \
                             tk_yyyzz_xxyyy,   \
                             tk_yyyzz_xxyyyy,  \
                             tk_yyyzz_xxyyyz,  \
                             tk_yyyzz_xxyyz,   \
                             tk_yyyzz_xxyyzz,  \
                             tk_yyyzz_xxyzz,   \
                             tk_yyyzz_xxyzzz,  \
                             tk_yyyzz_xxzzz,   \
                             tk_yyyzz_xxzzzz,  \
                             tk_yyyzz_xyyyy,   \
                             tk_yyyzz_xyyyyy,  \
                             tk_yyyzz_xyyyyz,  \
                             tk_yyyzz_xyyyz,   \
                             tk_yyyzz_xyyyzz,  \
                             tk_yyyzz_xyyzz,   \
                             tk_yyyzz_xyyzzz,  \
                             tk_yyyzz_xyzzz,   \
                             tk_yyyzz_xyzzzz,  \
                             tk_yyyzz_xzzzz,   \
                             tk_yyyzz_xzzzzz,  \
                             tk_yyyzz_yyyyy,   \
                             tk_yyyzz_yyyyyy,  \
                             tk_yyyzz_yyyyyz,  \
                             tk_yyyzz_yyyyz,   \
                             tk_yyyzz_yyyyzz,  \
                             tk_yyyzz_yyyzz,   \
                             tk_yyyzz_yyyzzz,  \
                             tk_yyyzz_yyzzz,   \
                             tk_yyyzz_yyzzzz,  \
                             tk_yyyzz_yzzzz,   \
                             tk_yyyzz_yzzzzz,  \
                             tk_yyyzz_zzzzz,   \
                             tk_yyyzz_zzzzzz,  \
                             ts_xyyyzz_xxxxxx, \
                             ts_xyyyzz_xxxxxy, \
                             ts_xyyyzz_xxxxxz, \
                             ts_xyyyzz_xxxxyy, \
                             ts_xyyyzz_xxxxyz, \
                             ts_xyyyzz_xxxxzz, \
                             ts_xyyyzz_xxxyyy, \
                             ts_xyyyzz_xxxyyz, \
                             ts_xyyyzz_xxxyzz, \
                             ts_xyyyzz_xxxzzz, \
                             ts_xyyyzz_xxyyyy, \
                             ts_xyyyzz_xxyyyz, \
                             ts_xyyyzz_xxyyzz, \
                             ts_xyyyzz_xxyzzz, \
                             ts_xyyyzz_xxzzzz, \
                             ts_xyyyzz_xyyyyy, \
                             ts_xyyyzz_xyyyyz, \
                             ts_xyyyzz_xyyyzz, \
                             ts_xyyyzz_xyyzzz, \
                             ts_xyyyzz_xyzzzz, \
                             ts_xyyyzz_xzzzzz, \
                             ts_xyyyzz_yyyyyy, \
                             ts_xyyyzz_yyyyyz, \
                             ts_xyyyzz_yyyyzz, \
                             ts_xyyyzz_yyyzzz, \
                             ts_xyyyzz_yyzzzz, \
                             ts_xyyyzz_yzzzzz, \
                             ts_xyyyzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyzz_xxxxxx[i] = 6.0 * tk_yyyzz_xxxxx[i] * fe_0 + tk_yyyzz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxxx[i] * fz_0;

        tk_xyyyzz_xxxxxy[i] = 5.0 * tk_yyyzz_xxxxy[i] * fe_0 + tk_yyyzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxxy[i] * fz_0;

        tk_xyyyzz_xxxxxz[i] = 5.0 * tk_yyyzz_xxxxz[i] * fe_0 + tk_yyyzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxxz[i] * fz_0;

        tk_xyyyzz_xxxxyy[i] = 4.0 * tk_yyyzz_xxxyy[i] * fe_0 + tk_yyyzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxyy[i] * fz_0;

        tk_xyyyzz_xxxxyz[i] = 4.0 * tk_yyyzz_xxxyz[i] * fe_0 + tk_yyyzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxyz[i] * fz_0;

        tk_xyyyzz_xxxxzz[i] = 4.0 * tk_yyyzz_xxxzz[i] * fe_0 + tk_yyyzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxzz[i] * fz_0;

        tk_xyyyzz_xxxyyy[i] = 3.0 * tk_yyyzz_xxyyy[i] * fe_0 + tk_yyyzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxyyy[i] * fz_0;

        tk_xyyyzz_xxxyyz[i] = 3.0 * tk_yyyzz_xxyyz[i] * fe_0 + tk_yyyzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxyyz[i] * fz_0;

        tk_xyyyzz_xxxyzz[i] = 3.0 * tk_yyyzz_xxyzz[i] * fe_0 + tk_yyyzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxyzz[i] * fz_0;

        tk_xyyyzz_xxxzzz[i] = 3.0 * tk_yyyzz_xxzzz[i] * fe_0 + tk_yyyzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxzzz[i] * fz_0;

        tk_xyyyzz_xxyyyy[i] = 2.0 * tk_yyyzz_xyyyy[i] * fe_0 + tk_yyyzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyyyy[i] * fz_0;

        tk_xyyyzz_xxyyyz[i] = 2.0 * tk_yyyzz_xyyyz[i] * fe_0 + tk_yyyzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyyyz[i] * fz_0;

        tk_xyyyzz_xxyyzz[i] = 2.0 * tk_yyyzz_xyyzz[i] * fe_0 + tk_yyyzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyyzz[i] * fz_0;

        tk_xyyyzz_xxyzzz[i] = 2.0 * tk_yyyzz_xyzzz[i] * fe_0 + tk_yyyzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyzzz[i] * fz_0;

        tk_xyyyzz_xxzzzz[i] = 2.0 * tk_yyyzz_xzzzz[i] * fe_0 + tk_yyyzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxzzzz[i] * fz_0;

        tk_xyyyzz_xyyyyy[i] = tk_yyyzz_yyyyy[i] * fe_0 + tk_yyyzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyyyy[i] * fz_0;

        tk_xyyyzz_xyyyyz[i] = tk_yyyzz_yyyyz[i] * fe_0 + tk_yyyzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyyyz[i] * fz_0;

        tk_xyyyzz_xyyyzz[i] = tk_yyyzz_yyyzz[i] * fe_0 + tk_yyyzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyyzz[i] * fz_0;

        tk_xyyyzz_xyyzzz[i] = tk_yyyzz_yyzzz[i] * fe_0 + tk_yyyzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyzzz[i] * fz_0;

        tk_xyyyzz_xyzzzz[i] = tk_yyyzz_yzzzz[i] * fe_0 + tk_yyyzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyzzzz[i] * fz_0;

        tk_xyyyzz_xzzzzz[i] = tk_yyyzz_zzzzz[i] * fe_0 + tk_yyyzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xzzzzz[i] * fz_0;

        tk_xyyyzz_yyyyyy[i] = tk_yyyzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyyyy[i] * fz_0;

        tk_xyyyzz_yyyyyz[i] = tk_yyyzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyyyz[i] * fz_0;

        tk_xyyyzz_yyyyzz[i] = tk_yyyzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyyzz[i] * fz_0;

        tk_xyyyzz_yyyzzz[i] = tk_yyyzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyzzz[i] * fz_0;

        tk_xyyyzz_yyzzzz[i] = tk_yyyzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyzzzz[i] * fz_0;

        tk_xyyyzz_yzzzzz[i] = tk_yyyzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yzzzzz[i] * fz_0;

        tk_xyyyzz_zzzzzz[i] = tk_yyyzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_zzzzzz[i] * fz_0;
    }

    // Set up 504-532 components of targeted buffer : II

    auto tk_xyyzzz_xxxxxx = pbuffer.data(idx_kin_ii + 504);

    auto tk_xyyzzz_xxxxxy = pbuffer.data(idx_kin_ii + 505);

    auto tk_xyyzzz_xxxxxz = pbuffer.data(idx_kin_ii + 506);

    auto tk_xyyzzz_xxxxyy = pbuffer.data(idx_kin_ii + 507);

    auto tk_xyyzzz_xxxxyz = pbuffer.data(idx_kin_ii + 508);

    auto tk_xyyzzz_xxxxzz = pbuffer.data(idx_kin_ii + 509);

    auto tk_xyyzzz_xxxyyy = pbuffer.data(idx_kin_ii + 510);

    auto tk_xyyzzz_xxxyyz = pbuffer.data(idx_kin_ii + 511);

    auto tk_xyyzzz_xxxyzz = pbuffer.data(idx_kin_ii + 512);

    auto tk_xyyzzz_xxxzzz = pbuffer.data(idx_kin_ii + 513);

    auto tk_xyyzzz_xxyyyy = pbuffer.data(idx_kin_ii + 514);

    auto tk_xyyzzz_xxyyyz = pbuffer.data(idx_kin_ii + 515);

    auto tk_xyyzzz_xxyyzz = pbuffer.data(idx_kin_ii + 516);

    auto tk_xyyzzz_xxyzzz = pbuffer.data(idx_kin_ii + 517);

    auto tk_xyyzzz_xxzzzz = pbuffer.data(idx_kin_ii + 518);

    auto tk_xyyzzz_xyyyyy = pbuffer.data(idx_kin_ii + 519);

    auto tk_xyyzzz_xyyyyz = pbuffer.data(idx_kin_ii + 520);

    auto tk_xyyzzz_xyyyzz = pbuffer.data(idx_kin_ii + 521);

    auto tk_xyyzzz_xyyzzz = pbuffer.data(idx_kin_ii + 522);

    auto tk_xyyzzz_xyzzzz = pbuffer.data(idx_kin_ii + 523);

    auto tk_xyyzzz_xzzzzz = pbuffer.data(idx_kin_ii + 524);

    auto tk_xyyzzz_yyyyyy = pbuffer.data(idx_kin_ii + 525);

    auto tk_xyyzzz_yyyyyz = pbuffer.data(idx_kin_ii + 526);

    auto tk_xyyzzz_yyyyzz = pbuffer.data(idx_kin_ii + 527);

    auto tk_xyyzzz_yyyzzz = pbuffer.data(idx_kin_ii + 528);

    auto tk_xyyzzz_yyzzzz = pbuffer.data(idx_kin_ii + 529);

    auto tk_xyyzzz_yzzzzz = pbuffer.data(idx_kin_ii + 530);

    auto tk_xyyzzz_zzzzzz = pbuffer.data(idx_kin_ii + 531);

#pragma omp simd aligned(pa_x,                 \
                             tk_xyyzzz_xxxxxx, \
                             tk_xyyzzz_xxxxxy, \
                             tk_xyyzzz_xxxxxz, \
                             tk_xyyzzz_xxxxyy, \
                             tk_xyyzzz_xxxxyz, \
                             tk_xyyzzz_xxxxzz, \
                             tk_xyyzzz_xxxyyy, \
                             tk_xyyzzz_xxxyyz, \
                             tk_xyyzzz_xxxyzz, \
                             tk_xyyzzz_xxxzzz, \
                             tk_xyyzzz_xxyyyy, \
                             tk_xyyzzz_xxyyyz, \
                             tk_xyyzzz_xxyyzz, \
                             tk_xyyzzz_xxyzzz, \
                             tk_xyyzzz_xxzzzz, \
                             tk_xyyzzz_xyyyyy, \
                             tk_xyyzzz_xyyyyz, \
                             tk_xyyzzz_xyyyzz, \
                             tk_xyyzzz_xyyzzz, \
                             tk_xyyzzz_xyzzzz, \
                             tk_xyyzzz_xzzzzz, \
                             tk_xyyzzz_yyyyyy, \
                             tk_xyyzzz_yyyyyz, \
                             tk_xyyzzz_yyyyzz, \
                             tk_xyyzzz_yyyzzz, \
                             tk_xyyzzz_yyzzzz, \
                             tk_xyyzzz_yzzzzz, \
                             tk_xyyzzz_zzzzzz, \
                             tk_yyzzz_xxxxx,   \
                             tk_yyzzz_xxxxxx,  \
                             tk_yyzzz_xxxxxy,  \
                             tk_yyzzz_xxxxxz,  \
                             tk_yyzzz_xxxxy,   \
                             tk_yyzzz_xxxxyy,  \
                             tk_yyzzz_xxxxyz,  \
                             tk_yyzzz_xxxxz,   \
                             tk_yyzzz_xxxxzz,  \
                             tk_yyzzz_xxxyy,   \
                             tk_yyzzz_xxxyyy,  \
                             tk_yyzzz_xxxyyz,  \
                             tk_yyzzz_xxxyz,   \
                             tk_yyzzz_xxxyzz,  \
                             tk_yyzzz_xxxzz,   \
                             tk_yyzzz_xxxzzz,  \
                             tk_yyzzz_xxyyy,   \
                             tk_yyzzz_xxyyyy,  \
                             tk_yyzzz_xxyyyz,  \
                             tk_yyzzz_xxyyz,   \
                             tk_yyzzz_xxyyzz,  \
                             tk_yyzzz_xxyzz,   \
                             tk_yyzzz_xxyzzz,  \
                             tk_yyzzz_xxzzz,   \
                             tk_yyzzz_xxzzzz,  \
                             tk_yyzzz_xyyyy,   \
                             tk_yyzzz_xyyyyy,  \
                             tk_yyzzz_xyyyyz,  \
                             tk_yyzzz_xyyyz,   \
                             tk_yyzzz_xyyyzz,  \
                             tk_yyzzz_xyyzz,   \
                             tk_yyzzz_xyyzzz,  \
                             tk_yyzzz_xyzzz,   \
                             tk_yyzzz_xyzzzz,  \
                             tk_yyzzz_xzzzz,   \
                             tk_yyzzz_xzzzzz,  \
                             tk_yyzzz_yyyyy,   \
                             tk_yyzzz_yyyyyy,  \
                             tk_yyzzz_yyyyyz,  \
                             tk_yyzzz_yyyyz,   \
                             tk_yyzzz_yyyyzz,  \
                             tk_yyzzz_yyyzz,   \
                             tk_yyzzz_yyyzzz,  \
                             tk_yyzzz_yyzzz,   \
                             tk_yyzzz_yyzzzz,  \
                             tk_yyzzz_yzzzz,   \
                             tk_yyzzz_yzzzzz,  \
                             tk_yyzzz_zzzzz,   \
                             tk_yyzzz_zzzzzz,  \
                             ts_xyyzzz_xxxxxx, \
                             ts_xyyzzz_xxxxxy, \
                             ts_xyyzzz_xxxxxz, \
                             ts_xyyzzz_xxxxyy, \
                             ts_xyyzzz_xxxxyz, \
                             ts_xyyzzz_xxxxzz, \
                             ts_xyyzzz_xxxyyy, \
                             ts_xyyzzz_xxxyyz, \
                             ts_xyyzzz_xxxyzz, \
                             ts_xyyzzz_xxxzzz, \
                             ts_xyyzzz_xxyyyy, \
                             ts_xyyzzz_xxyyyz, \
                             ts_xyyzzz_xxyyzz, \
                             ts_xyyzzz_xxyzzz, \
                             ts_xyyzzz_xxzzzz, \
                             ts_xyyzzz_xyyyyy, \
                             ts_xyyzzz_xyyyyz, \
                             ts_xyyzzz_xyyyzz, \
                             ts_xyyzzz_xyyzzz, \
                             ts_xyyzzz_xyzzzz, \
                             ts_xyyzzz_xzzzzz, \
                             ts_xyyzzz_yyyyyy, \
                             ts_xyyzzz_yyyyyz, \
                             ts_xyyzzz_yyyyzz, \
                             ts_xyyzzz_yyyzzz, \
                             ts_xyyzzz_yyzzzz, \
                             ts_xyyzzz_yzzzzz, \
                             ts_xyyzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzzz_xxxxxx[i] = 6.0 * tk_yyzzz_xxxxx[i] * fe_0 + tk_yyzzz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxxx[i] * fz_0;

        tk_xyyzzz_xxxxxy[i] = 5.0 * tk_yyzzz_xxxxy[i] * fe_0 + tk_yyzzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxxy[i] * fz_0;

        tk_xyyzzz_xxxxxz[i] = 5.0 * tk_yyzzz_xxxxz[i] * fe_0 + tk_yyzzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxxz[i] * fz_0;

        tk_xyyzzz_xxxxyy[i] = 4.0 * tk_yyzzz_xxxyy[i] * fe_0 + tk_yyzzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxyy[i] * fz_0;

        tk_xyyzzz_xxxxyz[i] = 4.0 * tk_yyzzz_xxxyz[i] * fe_0 + tk_yyzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxyz[i] * fz_0;

        tk_xyyzzz_xxxxzz[i] = 4.0 * tk_yyzzz_xxxzz[i] * fe_0 + tk_yyzzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxzz[i] * fz_0;

        tk_xyyzzz_xxxyyy[i] = 3.0 * tk_yyzzz_xxyyy[i] * fe_0 + tk_yyzzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxyyy[i] * fz_0;

        tk_xyyzzz_xxxyyz[i] = 3.0 * tk_yyzzz_xxyyz[i] * fe_0 + tk_yyzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxyyz[i] * fz_0;

        tk_xyyzzz_xxxyzz[i] = 3.0 * tk_yyzzz_xxyzz[i] * fe_0 + tk_yyzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxyzz[i] * fz_0;

        tk_xyyzzz_xxxzzz[i] = 3.0 * tk_yyzzz_xxzzz[i] * fe_0 + tk_yyzzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxzzz[i] * fz_0;

        tk_xyyzzz_xxyyyy[i] = 2.0 * tk_yyzzz_xyyyy[i] * fe_0 + tk_yyzzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyyyy[i] * fz_0;

        tk_xyyzzz_xxyyyz[i] = 2.0 * tk_yyzzz_xyyyz[i] * fe_0 + tk_yyzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyyyz[i] * fz_0;

        tk_xyyzzz_xxyyzz[i] = 2.0 * tk_yyzzz_xyyzz[i] * fe_0 + tk_yyzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyyzz[i] * fz_0;

        tk_xyyzzz_xxyzzz[i] = 2.0 * tk_yyzzz_xyzzz[i] * fe_0 + tk_yyzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyzzz[i] * fz_0;

        tk_xyyzzz_xxzzzz[i] = 2.0 * tk_yyzzz_xzzzz[i] * fe_0 + tk_yyzzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxzzzz[i] * fz_0;

        tk_xyyzzz_xyyyyy[i] = tk_yyzzz_yyyyy[i] * fe_0 + tk_yyzzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyyyy[i] * fz_0;

        tk_xyyzzz_xyyyyz[i] = tk_yyzzz_yyyyz[i] * fe_0 + tk_yyzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyyyz[i] * fz_0;

        tk_xyyzzz_xyyyzz[i] = tk_yyzzz_yyyzz[i] * fe_0 + tk_yyzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyyzz[i] * fz_0;

        tk_xyyzzz_xyyzzz[i] = tk_yyzzz_yyzzz[i] * fe_0 + tk_yyzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyzzz[i] * fz_0;

        tk_xyyzzz_xyzzzz[i] = tk_yyzzz_yzzzz[i] * fe_0 + tk_yyzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyzzzz[i] * fz_0;

        tk_xyyzzz_xzzzzz[i] = tk_yyzzz_zzzzz[i] * fe_0 + tk_yyzzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xzzzzz[i] * fz_0;

        tk_xyyzzz_yyyyyy[i] = tk_yyzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyyyy[i] * fz_0;

        tk_xyyzzz_yyyyyz[i] = tk_yyzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyyyz[i] * fz_0;

        tk_xyyzzz_yyyyzz[i] = tk_yyzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyyzz[i] * fz_0;

        tk_xyyzzz_yyyzzz[i] = tk_yyzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyzzz[i] * fz_0;

        tk_xyyzzz_yyzzzz[i] = tk_yyzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyzzzz[i] * fz_0;

        tk_xyyzzz_yzzzzz[i] = tk_yyzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yzzzzz[i] * fz_0;

        tk_xyyzzz_zzzzzz[i] = tk_yyzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_zzzzzz[i] * fz_0;
    }

    // Set up 532-560 components of targeted buffer : II

    auto tk_xyzzzz_xxxxxx = pbuffer.data(idx_kin_ii + 532);

    auto tk_xyzzzz_xxxxxy = pbuffer.data(idx_kin_ii + 533);

    auto tk_xyzzzz_xxxxxz = pbuffer.data(idx_kin_ii + 534);

    auto tk_xyzzzz_xxxxyy = pbuffer.data(idx_kin_ii + 535);

    auto tk_xyzzzz_xxxxyz = pbuffer.data(idx_kin_ii + 536);

    auto tk_xyzzzz_xxxxzz = pbuffer.data(idx_kin_ii + 537);

    auto tk_xyzzzz_xxxyyy = pbuffer.data(idx_kin_ii + 538);

    auto tk_xyzzzz_xxxyyz = pbuffer.data(idx_kin_ii + 539);

    auto tk_xyzzzz_xxxyzz = pbuffer.data(idx_kin_ii + 540);

    auto tk_xyzzzz_xxxzzz = pbuffer.data(idx_kin_ii + 541);

    auto tk_xyzzzz_xxyyyy = pbuffer.data(idx_kin_ii + 542);

    auto tk_xyzzzz_xxyyyz = pbuffer.data(idx_kin_ii + 543);

    auto tk_xyzzzz_xxyyzz = pbuffer.data(idx_kin_ii + 544);

    auto tk_xyzzzz_xxyzzz = pbuffer.data(idx_kin_ii + 545);

    auto tk_xyzzzz_xxzzzz = pbuffer.data(idx_kin_ii + 546);

    auto tk_xyzzzz_xyyyyy = pbuffer.data(idx_kin_ii + 547);

    auto tk_xyzzzz_xyyyyz = pbuffer.data(idx_kin_ii + 548);

    auto tk_xyzzzz_xyyyzz = pbuffer.data(idx_kin_ii + 549);

    auto tk_xyzzzz_xyyzzz = pbuffer.data(idx_kin_ii + 550);

    auto tk_xyzzzz_xyzzzz = pbuffer.data(idx_kin_ii + 551);

    auto tk_xyzzzz_xzzzzz = pbuffer.data(idx_kin_ii + 552);

    auto tk_xyzzzz_yyyyyy = pbuffer.data(idx_kin_ii + 553);

    auto tk_xyzzzz_yyyyyz = pbuffer.data(idx_kin_ii + 554);

    auto tk_xyzzzz_yyyyzz = pbuffer.data(idx_kin_ii + 555);

    auto tk_xyzzzz_yyyzzz = pbuffer.data(idx_kin_ii + 556);

    auto tk_xyzzzz_yyzzzz = pbuffer.data(idx_kin_ii + 557);

    auto tk_xyzzzz_yzzzzz = pbuffer.data(idx_kin_ii + 558);

    auto tk_xyzzzz_zzzzzz = pbuffer.data(idx_kin_ii + 559);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tk_xyzzzz_xxxxxx, \
                             tk_xyzzzz_xxxxxy, \
                             tk_xyzzzz_xxxxxz, \
                             tk_xyzzzz_xxxxyy, \
                             tk_xyzzzz_xxxxyz, \
                             tk_xyzzzz_xxxxzz, \
                             tk_xyzzzz_xxxyyy, \
                             tk_xyzzzz_xxxyyz, \
                             tk_xyzzzz_xxxyzz, \
                             tk_xyzzzz_xxxzzz, \
                             tk_xyzzzz_xxyyyy, \
                             tk_xyzzzz_xxyyyz, \
                             tk_xyzzzz_xxyyzz, \
                             tk_xyzzzz_xxyzzz, \
                             tk_xyzzzz_xxzzzz, \
                             tk_xyzzzz_xyyyyy, \
                             tk_xyzzzz_xyyyyz, \
                             tk_xyzzzz_xyyyzz, \
                             tk_xyzzzz_xyyzzz, \
                             tk_xyzzzz_xyzzzz, \
                             tk_xyzzzz_xzzzzz, \
                             tk_xyzzzz_yyyyyy, \
                             tk_xyzzzz_yyyyyz, \
                             tk_xyzzzz_yyyyzz, \
                             tk_xyzzzz_yyyzzz, \
                             tk_xyzzzz_yyzzzz, \
                             tk_xyzzzz_yzzzzz, \
                             tk_xyzzzz_zzzzzz, \
                             tk_xzzzz_xxxxxx,  \
                             tk_xzzzz_xxxxxz,  \
                             tk_xzzzz_xxxxzz,  \
                             tk_xzzzz_xxxzzz,  \
                             tk_xzzzz_xxzzzz,  \
                             tk_xzzzz_xzzzzz,  \
                             tk_yzzzz_xxxxxy,  \
                             tk_yzzzz_xxxxy,   \
                             tk_yzzzz_xxxxyy,  \
                             tk_yzzzz_xxxxyz,  \
                             tk_yzzzz_xxxyy,   \
                             tk_yzzzz_xxxyyy,  \
                             tk_yzzzz_xxxyyz,  \
                             tk_yzzzz_xxxyz,   \
                             tk_yzzzz_xxxyzz,  \
                             tk_yzzzz_xxyyy,   \
                             tk_yzzzz_xxyyyy,  \
                             tk_yzzzz_xxyyyz,  \
                             tk_yzzzz_xxyyz,   \
                             tk_yzzzz_xxyyzz,  \
                             tk_yzzzz_xxyzz,   \
                             tk_yzzzz_xxyzzz,  \
                             tk_yzzzz_xyyyy,   \
                             tk_yzzzz_xyyyyy,  \
                             tk_yzzzz_xyyyyz,  \
                             tk_yzzzz_xyyyz,   \
                             tk_yzzzz_xyyyzz,  \
                             tk_yzzzz_xyyzz,   \
                             tk_yzzzz_xyyzzz,  \
                             tk_yzzzz_xyzzz,   \
                             tk_yzzzz_xyzzzz,  \
                             tk_yzzzz_yyyyy,   \
                             tk_yzzzz_yyyyyy,  \
                             tk_yzzzz_yyyyyz,  \
                             tk_yzzzz_yyyyz,   \
                             tk_yzzzz_yyyyzz,  \
                             tk_yzzzz_yyyzz,   \
                             tk_yzzzz_yyyzzz,  \
                             tk_yzzzz_yyzzz,   \
                             tk_yzzzz_yyzzzz,  \
                             tk_yzzzz_yzzzz,   \
                             tk_yzzzz_yzzzzz,  \
                             tk_yzzzz_zzzzzz,  \
                             ts_xyzzzz_xxxxxx, \
                             ts_xyzzzz_xxxxxy, \
                             ts_xyzzzz_xxxxxz, \
                             ts_xyzzzz_xxxxyy, \
                             ts_xyzzzz_xxxxyz, \
                             ts_xyzzzz_xxxxzz, \
                             ts_xyzzzz_xxxyyy, \
                             ts_xyzzzz_xxxyyz, \
                             ts_xyzzzz_xxxyzz, \
                             ts_xyzzzz_xxxzzz, \
                             ts_xyzzzz_xxyyyy, \
                             ts_xyzzzz_xxyyyz, \
                             ts_xyzzzz_xxyyzz, \
                             ts_xyzzzz_xxyzzz, \
                             ts_xyzzzz_xxzzzz, \
                             ts_xyzzzz_xyyyyy, \
                             ts_xyzzzz_xyyyyz, \
                             ts_xyzzzz_xyyyzz, \
                             ts_xyzzzz_xyyzzz, \
                             ts_xyzzzz_xyzzzz, \
                             ts_xyzzzz_xzzzzz, \
                             ts_xyzzzz_yyyyyy, \
                             ts_xyzzzz_yyyyyz, \
                             ts_xyzzzz_yyyyzz, \
                             ts_xyzzzz_yyyzzz, \
                             ts_xyzzzz_yyzzzz, \
                             ts_xyzzzz_yzzzzz, \
                             ts_xyzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzzz_xxxxxx[i] = tk_xzzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxxxx[i] * fz_0;

        tk_xyzzzz_xxxxxy[i] = 5.0 * tk_yzzzz_xxxxy[i] * fe_0 + tk_yzzzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxxxy[i] * fz_0;

        tk_xyzzzz_xxxxxz[i] = tk_xzzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxxxz[i] * fz_0;

        tk_xyzzzz_xxxxyy[i] = 4.0 * tk_yzzzz_xxxyy[i] * fe_0 + tk_yzzzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxxyy[i] * fz_0;

        tk_xyzzzz_xxxxyz[i] = 4.0 * tk_yzzzz_xxxyz[i] * fe_0 + tk_yzzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxxyz[i] * fz_0;

        tk_xyzzzz_xxxxzz[i] = tk_xzzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxxzz[i] * fz_0;

        tk_xyzzzz_xxxyyy[i] = 3.0 * tk_yzzzz_xxyyy[i] * fe_0 + tk_yzzzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxyyy[i] * fz_0;

        tk_xyzzzz_xxxyyz[i] = 3.0 * tk_yzzzz_xxyyz[i] * fe_0 + tk_yzzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxyyz[i] * fz_0;

        tk_xyzzzz_xxxyzz[i] = 3.0 * tk_yzzzz_xxyzz[i] * fe_0 + tk_yzzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxyzz[i] * fz_0;

        tk_xyzzzz_xxxzzz[i] = tk_xzzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxzzz[i] * fz_0;

        tk_xyzzzz_xxyyyy[i] = 2.0 * tk_yzzzz_xyyyy[i] * fe_0 + tk_yzzzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyyyy[i] * fz_0;

        tk_xyzzzz_xxyyyz[i] = 2.0 * tk_yzzzz_xyyyz[i] * fe_0 + tk_yzzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyyyz[i] * fz_0;

        tk_xyzzzz_xxyyzz[i] = 2.0 * tk_yzzzz_xyyzz[i] * fe_0 + tk_yzzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyyzz[i] * fz_0;

        tk_xyzzzz_xxyzzz[i] = 2.0 * tk_yzzzz_xyzzz[i] * fe_0 + tk_yzzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyzzz[i] * fz_0;

        tk_xyzzzz_xxzzzz[i] = tk_xzzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxzzzz[i] * fz_0;

        tk_xyzzzz_xyyyyy[i] = tk_yzzzz_yyyyy[i] * fe_0 + tk_yzzzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyyyy[i] * fz_0;

        tk_xyzzzz_xyyyyz[i] = tk_yzzzz_yyyyz[i] * fe_0 + tk_yzzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyyyz[i] * fz_0;

        tk_xyzzzz_xyyyzz[i] = tk_yzzzz_yyyzz[i] * fe_0 + tk_yzzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyyzz[i] * fz_0;

        tk_xyzzzz_xyyzzz[i] = tk_yzzzz_yyzzz[i] * fe_0 + tk_yzzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyzzz[i] * fz_0;

        tk_xyzzzz_xyzzzz[i] = tk_yzzzz_yzzzz[i] * fe_0 + tk_yzzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyzzzz[i] * fz_0;

        tk_xyzzzz_xzzzzz[i] = tk_xzzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xzzzzz[i] * fz_0;

        tk_xyzzzz_yyyyyy[i] = tk_yzzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyyyy[i] * fz_0;

        tk_xyzzzz_yyyyyz[i] = tk_yzzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyyyz[i] * fz_0;

        tk_xyzzzz_yyyyzz[i] = tk_yzzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyyzz[i] * fz_0;

        tk_xyzzzz_yyyzzz[i] = tk_yzzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyzzz[i] * fz_0;

        tk_xyzzzz_yyzzzz[i] = tk_yzzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyzzzz[i] * fz_0;

        tk_xyzzzz_yzzzzz[i] = tk_yzzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yzzzzz[i] * fz_0;

        tk_xyzzzz_zzzzzz[i] = tk_yzzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 560-588 components of targeted buffer : II

    auto tk_xzzzzz_xxxxxx = pbuffer.data(idx_kin_ii + 560);

    auto tk_xzzzzz_xxxxxy = pbuffer.data(idx_kin_ii + 561);

    auto tk_xzzzzz_xxxxxz = pbuffer.data(idx_kin_ii + 562);

    auto tk_xzzzzz_xxxxyy = pbuffer.data(idx_kin_ii + 563);

    auto tk_xzzzzz_xxxxyz = pbuffer.data(idx_kin_ii + 564);

    auto tk_xzzzzz_xxxxzz = pbuffer.data(idx_kin_ii + 565);

    auto tk_xzzzzz_xxxyyy = pbuffer.data(idx_kin_ii + 566);

    auto tk_xzzzzz_xxxyyz = pbuffer.data(idx_kin_ii + 567);

    auto tk_xzzzzz_xxxyzz = pbuffer.data(idx_kin_ii + 568);

    auto tk_xzzzzz_xxxzzz = pbuffer.data(idx_kin_ii + 569);

    auto tk_xzzzzz_xxyyyy = pbuffer.data(idx_kin_ii + 570);

    auto tk_xzzzzz_xxyyyz = pbuffer.data(idx_kin_ii + 571);

    auto tk_xzzzzz_xxyyzz = pbuffer.data(idx_kin_ii + 572);

    auto tk_xzzzzz_xxyzzz = pbuffer.data(idx_kin_ii + 573);

    auto tk_xzzzzz_xxzzzz = pbuffer.data(idx_kin_ii + 574);

    auto tk_xzzzzz_xyyyyy = pbuffer.data(idx_kin_ii + 575);

    auto tk_xzzzzz_xyyyyz = pbuffer.data(idx_kin_ii + 576);

    auto tk_xzzzzz_xyyyzz = pbuffer.data(idx_kin_ii + 577);

    auto tk_xzzzzz_xyyzzz = pbuffer.data(idx_kin_ii + 578);

    auto tk_xzzzzz_xyzzzz = pbuffer.data(idx_kin_ii + 579);

    auto tk_xzzzzz_xzzzzz = pbuffer.data(idx_kin_ii + 580);

    auto tk_xzzzzz_yyyyyy = pbuffer.data(idx_kin_ii + 581);

    auto tk_xzzzzz_yyyyyz = pbuffer.data(idx_kin_ii + 582);

    auto tk_xzzzzz_yyyyzz = pbuffer.data(idx_kin_ii + 583);

    auto tk_xzzzzz_yyyzzz = pbuffer.data(idx_kin_ii + 584);

    auto tk_xzzzzz_yyzzzz = pbuffer.data(idx_kin_ii + 585);

    auto tk_xzzzzz_yzzzzz = pbuffer.data(idx_kin_ii + 586);

    auto tk_xzzzzz_zzzzzz = pbuffer.data(idx_kin_ii + 587);

#pragma omp simd aligned(pa_x,                 \
                             tk_xzzzzz_xxxxxx, \
                             tk_xzzzzz_xxxxxy, \
                             tk_xzzzzz_xxxxxz, \
                             tk_xzzzzz_xxxxyy, \
                             tk_xzzzzz_xxxxyz, \
                             tk_xzzzzz_xxxxzz, \
                             tk_xzzzzz_xxxyyy, \
                             tk_xzzzzz_xxxyyz, \
                             tk_xzzzzz_xxxyzz, \
                             tk_xzzzzz_xxxzzz, \
                             tk_xzzzzz_xxyyyy, \
                             tk_xzzzzz_xxyyyz, \
                             tk_xzzzzz_xxyyzz, \
                             tk_xzzzzz_xxyzzz, \
                             tk_xzzzzz_xxzzzz, \
                             tk_xzzzzz_xyyyyy, \
                             tk_xzzzzz_xyyyyz, \
                             tk_xzzzzz_xyyyzz, \
                             tk_xzzzzz_xyyzzz, \
                             tk_xzzzzz_xyzzzz, \
                             tk_xzzzzz_xzzzzz, \
                             tk_xzzzzz_yyyyyy, \
                             tk_xzzzzz_yyyyyz, \
                             tk_xzzzzz_yyyyzz, \
                             tk_xzzzzz_yyyzzz, \
                             tk_xzzzzz_yyzzzz, \
                             tk_xzzzzz_yzzzzz, \
                             tk_xzzzzz_zzzzzz, \
                             tk_zzzzz_xxxxx,   \
                             tk_zzzzz_xxxxxx,  \
                             tk_zzzzz_xxxxxy,  \
                             tk_zzzzz_xxxxxz,  \
                             tk_zzzzz_xxxxy,   \
                             tk_zzzzz_xxxxyy,  \
                             tk_zzzzz_xxxxyz,  \
                             tk_zzzzz_xxxxz,   \
                             tk_zzzzz_xxxxzz,  \
                             tk_zzzzz_xxxyy,   \
                             tk_zzzzz_xxxyyy,  \
                             tk_zzzzz_xxxyyz,  \
                             tk_zzzzz_xxxyz,   \
                             tk_zzzzz_xxxyzz,  \
                             tk_zzzzz_xxxzz,   \
                             tk_zzzzz_xxxzzz,  \
                             tk_zzzzz_xxyyy,   \
                             tk_zzzzz_xxyyyy,  \
                             tk_zzzzz_xxyyyz,  \
                             tk_zzzzz_xxyyz,   \
                             tk_zzzzz_xxyyzz,  \
                             tk_zzzzz_xxyzz,   \
                             tk_zzzzz_xxyzzz,  \
                             tk_zzzzz_xxzzz,   \
                             tk_zzzzz_xxzzzz,  \
                             tk_zzzzz_xyyyy,   \
                             tk_zzzzz_xyyyyy,  \
                             tk_zzzzz_xyyyyz,  \
                             tk_zzzzz_xyyyz,   \
                             tk_zzzzz_xyyyzz,  \
                             tk_zzzzz_xyyzz,   \
                             tk_zzzzz_xyyzzz,  \
                             tk_zzzzz_xyzzz,   \
                             tk_zzzzz_xyzzzz,  \
                             tk_zzzzz_xzzzz,   \
                             tk_zzzzz_xzzzzz,  \
                             tk_zzzzz_yyyyy,   \
                             tk_zzzzz_yyyyyy,  \
                             tk_zzzzz_yyyyyz,  \
                             tk_zzzzz_yyyyz,   \
                             tk_zzzzz_yyyyzz,  \
                             tk_zzzzz_yyyzz,   \
                             tk_zzzzz_yyyzzz,  \
                             tk_zzzzz_yyzzz,   \
                             tk_zzzzz_yyzzzz,  \
                             tk_zzzzz_yzzzz,   \
                             tk_zzzzz_yzzzzz,  \
                             tk_zzzzz_zzzzz,   \
                             tk_zzzzz_zzzzzz,  \
                             ts_xzzzzz_xxxxxx, \
                             ts_xzzzzz_xxxxxy, \
                             ts_xzzzzz_xxxxxz, \
                             ts_xzzzzz_xxxxyy, \
                             ts_xzzzzz_xxxxyz, \
                             ts_xzzzzz_xxxxzz, \
                             ts_xzzzzz_xxxyyy, \
                             ts_xzzzzz_xxxyyz, \
                             ts_xzzzzz_xxxyzz, \
                             ts_xzzzzz_xxxzzz, \
                             ts_xzzzzz_xxyyyy, \
                             ts_xzzzzz_xxyyyz, \
                             ts_xzzzzz_xxyyzz, \
                             ts_xzzzzz_xxyzzz, \
                             ts_xzzzzz_xxzzzz, \
                             ts_xzzzzz_xyyyyy, \
                             ts_xzzzzz_xyyyyz, \
                             ts_xzzzzz_xyyyzz, \
                             ts_xzzzzz_xyyzzz, \
                             ts_xzzzzz_xyzzzz, \
                             ts_xzzzzz_xzzzzz, \
                             ts_xzzzzz_yyyyyy, \
                             ts_xzzzzz_yyyyyz, \
                             ts_xzzzzz_yyyyzz, \
                             ts_xzzzzz_yyyzzz, \
                             ts_xzzzzz_yyzzzz, \
                             ts_xzzzzz_yzzzzz, \
                             ts_xzzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzzz_xxxxxx[i] = 6.0 * tk_zzzzz_xxxxx[i] * fe_0 + tk_zzzzz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxxx[i] * fz_0;

        tk_xzzzzz_xxxxxy[i] = 5.0 * tk_zzzzz_xxxxy[i] * fe_0 + tk_zzzzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxxy[i] * fz_0;

        tk_xzzzzz_xxxxxz[i] = 5.0 * tk_zzzzz_xxxxz[i] * fe_0 + tk_zzzzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxxz[i] * fz_0;

        tk_xzzzzz_xxxxyy[i] = 4.0 * tk_zzzzz_xxxyy[i] * fe_0 + tk_zzzzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxyy[i] * fz_0;

        tk_xzzzzz_xxxxyz[i] = 4.0 * tk_zzzzz_xxxyz[i] * fe_0 + tk_zzzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxyz[i] * fz_0;

        tk_xzzzzz_xxxxzz[i] = 4.0 * tk_zzzzz_xxxzz[i] * fe_0 + tk_zzzzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxzz[i] * fz_0;

        tk_xzzzzz_xxxyyy[i] = 3.0 * tk_zzzzz_xxyyy[i] * fe_0 + tk_zzzzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxyyy[i] * fz_0;

        tk_xzzzzz_xxxyyz[i] = 3.0 * tk_zzzzz_xxyyz[i] * fe_0 + tk_zzzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxyyz[i] * fz_0;

        tk_xzzzzz_xxxyzz[i] = 3.0 * tk_zzzzz_xxyzz[i] * fe_0 + tk_zzzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxyzz[i] * fz_0;

        tk_xzzzzz_xxxzzz[i] = 3.0 * tk_zzzzz_xxzzz[i] * fe_0 + tk_zzzzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxzzz[i] * fz_0;

        tk_xzzzzz_xxyyyy[i] = 2.0 * tk_zzzzz_xyyyy[i] * fe_0 + tk_zzzzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyyyy[i] * fz_0;

        tk_xzzzzz_xxyyyz[i] = 2.0 * tk_zzzzz_xyyyz[i] * fe_0 + tk_zzzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyyyz[i] * fz_0;

        tk_xzzzzz_xxyyzz[i] = 2.0 * tk_zzzzz_xyyzz[i] * fe_0 + tk_zzzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyyzz[i] * fz_0;

        tk_xzzzzz_xxyzzz[i] = 2.0 * tk_zzzzz_xyzzz[i] * fe_0 + tk_zzzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyzzz[i] * fz_0;

        tk_xzzzzz_xxzzzz[i] = 2.0 * tk_zzzzz_xzzzz[i] * fe_0 + tk_zzzzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxzzzz[i] * fz_0;

        tk_xzzzzz_xyyyyy[i] = tk_zzzzz_yyyyy[i] * fe_0 + tk_zzzzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyyyy[i] * fz_0;

        tk_xzzzzz_xyyyyz[i] = tk_zzzzz_yyyyz[i] * fe_0 + tk_zzzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyyyz[i] * fz_0;

        tk_xzzzzz_xyyyzz[i] = tk_zzzzz_yyyzz[i] * fe_0 + tk_zzzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyyzz[i] * fz_0;

        tk_xzzzzz_xyyzzz[i] = tk_zzzzz_yyzzz[i] * fe_0 + tk_zzzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyzzz[i] * fz_0;

        tk_xzzzzz_xyzzzz[i] = tk_zzzzz_yzzzz[i] * fe_0 + tk_zzzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyzzzz[i] * fz_0;

        tk_xzzzzz_xzzzzz[i] = tk_zzzzz_zzzzz[i] * fe_0 + tk_zzzzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xzzzzz[i] * fz_0;

        tk_xzzzzz_yyyyyy[i] = tk_zzzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyyyy[i] * fz_0;

        tk_xzzzzz_yyyyyz[i] = tk_zzzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyyyz[i] * fz_0;

        tk_xzzzzz_yyyyzz[i] = tk_zzzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyyzz[i] * fz_0;

        tk_xzzzzz_yyyzzz[i] = tk_zzzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyzzz[i] * fz_0;

        tk_xzzzzz_yyzzzz[i] = tk_zzzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyzzzz[i] * fz_0;

        tk_xzzzzz_yzzzzz[i] = tk_zzzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yzzzzz[i] * fz_0;

        tk_xzzzzz_zzzzzz[i] = tk_zzzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 588-616 components of targeted buffer : II

    auto tk_yyyyyy_xxxxxx = pbuffer.data(idx_kin_ii + 588);

    auto tk_yyyyyy_xxxxxy = pbuffer.data(idx_kin_ii + 589);

    auto tk_yyyyyy_xxxxxz = pbuffer.data(idx_kin_ii + 590);

    auto tk_yyyyyy_xxxxyy = pbuffer.data(idx_kin_ii + 591);

    auto tk_yyyyyy_xxxxyz = pbuffer.data(idx_kin_ii + 592);

    auto tk_yyyyyy_xxxxzz = pbuffer.data(idx_kin_ii + 593);

    auto tk_yyyyyy_xxxyyy = pbuffer.data(idx_kin_ii + 594);

    auto tk_yyyyyy_xxxyyz = pbuffer.data(idx_kin_ii + 595);

    auto tk_yyyyyy_xxxyzz = pbuffer.data(idx_kin_ii + 596);

    auto tk_yyyyyy_xxxzzz = pbuffer.data(idx_kin_ii + 597);

    auto tk_yyyyyy_xxyyyy = pbuffer.data(idx_kin_ii + 598);

    auto tk_yyyyyy_xxyyyz = pbuffer.data(idx_kin_ii + 599);

    auto tk_yyyyyy_xxyyzz = pbuffer.data(idx_kin_ii + 600);

    auto tk_yyyyyy_xxyzzz = pbuffer.data(idx_kin_ii + 601);

    auto tk_yyyyyy_xxzzzz = pbuffer.data(idx_kin_ii + 602);

    auto tk_yyyyyy_xyyyyy = pbuffer.data(idx_kin_ii + 603);

    auto tk_yyyyyy_xyyyyz = pbuffer.data(idx_kin_ii + 604);

    auto tk_yyyyyy_xyyyzz = pbuffer.data(idx_kin_ii + 605);

    auto tk_yyyyyy_xyyzzz = pbuffer.data(idx_kin_ii + 606);

    auto tk_yyyyyy_xyzzzz = pbuffer.data(idx_kin_ii + 607);

    auto tk_yyyyyy_xzzzzz = pbuffer.data(idx_kin_ii + 608);

    auto tk_yyyyyy_yyyyyy = pbuffer.data(idx_kin_ii + 609);

    auto tk_yyyyyy_yyyyyz = pbuffer.data(idx_kin_ii + 610);

    auto tk_yyyyyy_yyyyzz = pbuffer.data(idx_kin_ii + 611);

    auto tk_yyyyyy_yyyzzz = pbuffer.data(idx_kin_ii + 612);

    auto tk_yyyyyy_yyzzzz = pbuffer.data(idx_kin_ii + 613);

    auto tk_yyyyyy_yzzzzz = pbuffer.data(idx_kin_ii + 614);

    auto tk_yyyyyy_zzzzzz = pbuffer.data(idx_kin_ii + 615);

#pragma omp simd aligned(pa_y,                 \
                             tk_yyyy_xxxxxx,   \
                             tk_yyyy_xxxxxy,   \
                             tk_yyyy_xxxxxz,   \
                             tk_yyyy_xxxxyy,   \
                             tk_yyyy_xxxxyz,   \
                             tk_yyyy_xxxxzz,   \
                             tk_yyyy_xxxyyy,   \
                             tk_yyyy_xxxyyz,   \
                             tk_yyyy_xxxyzz,   \
                             tk_yyyy_xxxzzz,   \
                             tk_yyyy_xxyyyy,   \
                             tk_yyyy_xxyyyz,   \
                             tk_yyyy_xxyyzz,   \
                             tk_yyyy_xxyzzz,   \
                             tk_yyyy_xxzzzz,   \
                             tk_yyyy_xyyyyy,   \
                             tk_yyyy_xyyyyz,   \
                             tk_yyyy_xyyyzz,   \
                             tk_yyyy_xyyzzz,   \
                             tk_yyyy_xyzzzz,   \
                             tk_yyyy_xzzzzz,   \
                             tk_yyyy_yyyyyy,   \
                             tk_yyyy_yyyyyz,   \
                             tk_yyyy_yyyyzz,   \
                             tk_yyyy_yyyzzz,   \
                             tk_yyyy_yyzzzz,   \
                             tk_yyyy_yzzzzz,   \
                             tk_yyyy_zzzzzz,   \
                             tk_yyyyy_xxxxx,   \
                             tk_yyyyy_xxxxxx,  \
                             tk_yyyyy_xxxxxy,  \
                             tk_yyyyy_xxxxxz,  \
                             tk_yyyyy_xxxxy,   \
                             tk_yyyyy_xxxxyy,  \
                             tk_yyyyy_xxxxyz,  \
                             tk_yyyyy_xxxxz,   \
                             tk_yyyyy_xxxxzz,  \
                             tk_yyyyy_xxxyy,   \
                             tk_yyyyy_xxxyyy,  \
                             tk_yyyyy_xxxyyz,  \
                             tk_yyyyy_xxxyz,   \
                             tk_yyyyy_xxxyzz,  \
                             tk_yyyyy_xxxzz,   \
                             tk_yyyyy_xxxzzz,  \
                             tk_yyyyy_xxyyy,   \
                             tk_yyyyy_xxyyyy,  \
                             tk_yyyyy_xxyyyz,  \
                             tk_yyyyy_xxyyz,   \
                             tk_yyyyy_xxyyzz,  \
                             tk_yyyyy_xxyzz,   \
                             tk_yyyyy_xxyzzz,  \
                             tk_yyyyy_xxzzz,   \
                             tk_yyyyy_xxzzzz,  \
                             tk_yyyyy_xyyyy,   \
                             tk_yyyyy_xyyyyy,  \
                             tk_yyyyy_xyyyyz,  \
                             tk_yyyyy_xyyyz,   \
                             tk_yyyyy_xyyyzz,  \
                             tk_yyyyy_xyyzz,   \
                             tk_yyyyy_xyyzzz,  \
                             tk_yyyyy_xyzzz,   \
                             tk_yyyyy_xyzzzz,  \
                             tk_yyyyy_xzzzz,   \
                             tk_yyyyy_xzzzzz,  \
                             tk_yyyyy_yyyyy,   \
                             tk_yyyyy_yyyyyy,  \
                             tk_yyyyy_yyyyyz,  \
                             tk_yyyyy_yyyyz,   \
                             tk_yyyyy_yyyyzz,  \
                             tk_yyyyy_yyyzz,   \
                             tk_yyyyy_yyyzzz,  \
                             tk_yyyyy_yyzzz,   \
                             tk_yyyyy_yyzzzz,  \
                             tk_yyyyy_yzzzz,   \
                             tk_yyyyy_yzzzzz,  \
                             tk_yyyyy_zzzzz,   \
                             tk_yyyyy_zzzzzz,  \
                             tk_yyyyyy_xxxxxx, \
                             tk_yyyyyy_xxxxxy, \
                             tk_yyyyyy_xxxxxz, \
                             tk_yyyyyy_xxxxyy, \
                             tk_yyyyyy_xxxxyz, \
                             tk_yyyyyy_xxxxzz, \
                             tk_yyyyyy_xxxyyy, \
                             tk_yyyyyy_xxxyyz, \
                             tk_yyyyyy_xxxyzz, \
                             tk_yyyyyy_xxxzzz, \
                             tk_yyyyyy_xxyyyy, \
                             tk_yyyyyy_xxyyyz, \
                             tk_yyyyyy_xxyyzz, \
                             tk_yyyyyy_xxyzzz, \
                             tk_yyyyyy_xxzzzz, \
                             tk_yyyyyy_xyyyyy, \
                             tk_yyyyyy_xyyyyz, \
                             tk_yyyyyy_xyyyzz, \
                             tk_yyyyyy_xyyzzz, \
                             tk_yyyyyy_xyzzzz, \
                             tk_yyyyyy_xzzzzz, \
                             tk_yyyyyy_yyyyyy, \
                             tk_yyyyyy_yyyyyz, \
                             tk_yyyyyy_yyyyzz, \
                             tk_yyyyyy_yyyzzz, \
                             tk_yyyyyy_yyzzzz, \
                             tk_yyyyyy_yzzzzz, \
                             tk_yyyyyy_zzzzzz, \
                             ts_yyyy_xxxxxx,   \
                             ts_yyyy_xxxxxy,   \
                             ts_yyyy_xxxxxz,   \
                             ts_yyyy_xxxxyy,   \
                             ts_yyyy_xxxxyz,   \
                             ts_yyyy_xxxxzz,   \
                             ts_yyyy_xxxyyy,   \
                             ts_yyyy_xxxyyz,   \
                             ts_yyyy_xxxyzz,   \
                             ts_yyyy_xxxzzz,   \
                             ts_yyyy_xxyyyy,   \
                             ts_yyyy_xxyyyz,   \
                             ts_yyyy_xxyyzz,   \
                             ts_yyyy_xxyzzz,   \
                             ts_yyyy_xxzzzz,   \
                             ts_yyyy_xyyyyy,   \
                             ts_yyyy_xyyyyz,   \
                             ts_yyyy_xyyyzz,   \
                             ts_yyyy_xyyzzz,   \
                             ts_yyyy_xyzzzz,   \
                             ts_yyyy_xzzzzz,   \
                             ts_yyyy_yyyyyy,   \
                             ts_yyyy_yyyyyz,   \
                             ts_yyyy_yyyyzz,   \
                             ts_yyyy_yyyzzz,   \
                             ts_yyyy_yyzzzz,   \
                             ts_yyyy_yzzzzz,   \
                             ts_yyyy_zzzzzz,   \
                             ts_yyyyyy_xxxxxx, \
                             ts_yyyyyy_xxxxxy, \
                             ts_yyyyyy_xxxxxz, \
                             ts_yyyyyy_xxxxyy, \
                             ts_yyyyyy_xxxxyz, \
                             ts_yyyyyy_xxxxzz, \
                             ts_yyyyyy_xxxyyy, \
                             ts_yyyyyy_xxxyyz, \
                             ts_yyyyyy_xxxyzz, \
                             ts_yyyyyy_xxxzzz, \
                             ts_yyyyyy_xxyyyy, \
                             ts_yyyyyy_xxyyyz, \
                             ts_yyyyyy_xxyyzz, \
                             ts_yyyyyy_xxyzzz, \
                             ts_yyyyyy_xxzzzz, \
                             ts_yyyyyy_xyyyyy, \
                             ts_yyyyyy_xyyyyz, \
                             ts_yyyyyy_xyyyzz, \
                             ts_yyyyyy_xyyzzz, \
                             ts_yyyyyy_xyzzzz, \
                             ts_yyyyyy_xzzzzz, \
                             ts_yyyyyy_yyyyyy, \
                             ts_yyyyyy_yyyyyz, \
                             ts_yyyyyy_yyyyzz, \
                             ts_yyyyyy_yyyzzz, \
                             ts_yyyyyy_yyzzzz, \
                             ts_yyyyyy_yzzzzz, \
                             ts_yyyyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyyy_xxxxxx[i] = -10.0 * ts_yyyy_xxxxxx[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxxx[i] * fe_0 + tk_yyyyy_xxxxxx[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_xxxxxx[i] * fz_0;

        tk_yyyyyy_xxxxxy[i] = -10.0 * ts_yyyy_xxxxxy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxxy[i] * fe_0 + tk_yyyyy_xxxxx[i] * fe_0 +
                              tk_yyyyy_xxxxxy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxxxy[i] * fz_0;

        tk_yyyyyy_xxxxxz[i] = -10.0 * ts_yyyy_xxxxxz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxxz[i] * fe_0 + tk_yyyyy_xxxxxz[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_xxxxxz[i] * fz_0;

        tk_yyyyyy_xxxxyy[i] = -10.0 * ts_yyyy_xxxxyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxyy[i] * fe_0 + 2.0 * tk_yyyyy_xxxxy[i] * fe_0 +
                              tk_yyyyy_xxxxyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxxyy[i] * fz_0;

        tk_yyyyyy_xxxxyz[i] = -10.0 * ts_yyyy_xxxxyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxyz[i] * fe_0 + tk_yyyyy_xxxxz[i] * fe_0 +
                              tk_yyyyy_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxxyz[i] * fz_0;

        tk_yyyyyy_xxxxzz[i] = -10.0 * ts_yyyy_xxxxzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxzz[i] * fe_0 + tk_yyyyy_xxxxzz[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_xxxxzz[i] * fz_0;

        tk_yyyyyy_xxxyyy[i] = -10.0 * ts_yyyy_xxxyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxyyy[i] * fe_0 + 3.0 * tk_yyyyy_xxxyy[i] * fe_0 +
                              tk_yyyyy_xxxyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxyyy[i] * fz_0;

        tk_yyyyyy_xxxyyz[i] = -10.0 * ts_yyyy_xxxyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxyyz[i] * fe_0 + 2.0 * tk_yyyyy_xxxyz[i] * fe_0 +
                              tk_yyyyy_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxyyz[i] * fz_0;

        tk_yyyyyy_xxxyzz[i] = -10.0 * ts_yyyy_xxxyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxyzz[i] * fe_0 + tk_yyyyy_xxxzz[i] * fe_0 +
                              tk_yyyyy_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxyzz[i] * fz_0;

        tk_yyyyyy_xxxzzz[i] = -10.0 * ts_yyyy_xxxzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxzzz[i] * fe_0 + tk_yyyyy_xxxzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_xxxzzz[i] * fz_0;

        tk_yyyyyy_xxyyyy[i] = -10.0 * ts_yyyy_xxyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyyyy[i] * fe_0 + 4.0 * tk_yyyyy_xxyyy[i] * fe_0 +
                              tk_yyyyy_xxyyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyyyy[i] * fz_0;

        tk_yyyyyy_xxyyyz[i] = -10.0 * ts_yyyy_xxyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyyyz[i] * fe_0 + 3.0 * tk_yyyyy_xxyyz[i] * fe_0 +
                              tk_yyyyy_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyyyz[i] * fz_0;

        tk_yyyyyy_xxyyzz[i] = -10.0 * ts_yyyy_xxyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyyzz[i] * fe_0 + 2.0 * tk_yyyyy_xxyzz[i] * fe_0 +
                              tk_yyyyy_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyyzz[i] * fz_0;

        tk_yyyyyy_xxyzzz[i] = -10.0 * ts_yyyy_xxyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyzzz[i] * fe_0 + tk_yyyyy_xxzzz[i] * fe_0 +
                              tk_yyyyy_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyzzz[i] * fz_0;

        tk_yyyyyy_xxzzzz[i] = -10.0 * ts_yyyy_xxzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxzzzz[i] * fe_0 + tk_yyyyy_xxzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_xxzzzz[i] * fz_0;

        tk_yyyyyy_xyyyyy[i] = -10.0 * ts_yyyy_xyyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyyyy[i] * fe_0 + 5.0 * tk_yyyyy_xyyyy[i] * fe_0 +
                              tk_yyyyy_xyyyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyyyy[i] * fz_0;

        tk_yyyyyy_xyyyyz[i] = -10.0 * ts_yyyy_xyyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyyyz[i] * fe_0 + 4.0 * tk_yyyyy_xyyyz[i] * fe_0 +
                              tk_yyyyy_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyyyz[i] * fz_0;

        tk_yyyyyy_xyyyzz[i] = -10.0 * ts_yyyy_xyyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyyzz[i] * fe_0 + 3.0 * tk_yyyyy_xyyzz[i] * fe_0 +
                              tk_yyyyy_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyyzz[i] * fz_0;

        tk_yyyyyy_xyyzzz[i] = -10.0 * ts_yyyy_xyyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyzzz[i] * fe_0 + 2.0 * tk_yyyyy_xyzzz[i] * fe_0 +
                              tk_yyyyy_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyzzz[i] * fz_0;

        tk_yyyyyy_xyzzzz[i] = -10.0 * ts_yyyy_xyzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyzzzz[i] * fe_0 + tk_yyyyy_xzzzz[i] * fe_0 +
                              tk_yyyyy_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyzzzz[i] * fz_0;

        tk_yyyyyy_xzzzzz[i] = -10.0 * ts_yyyy_xzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xzzzzz[i] * fe_0 + tk_yyyyy_xzzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_xzzzzz[i] * fz_0;

        tk_yyyyyy_yyyyyy[i] = -10.0 * ts_yyyy_yyyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyyyy[i] * fe_0 + 6.0 * tk_yyyyy_yyyyy[i] * fe_0 +
                              tk_yyyyy_yyyyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyyyy[i] * fz_0;

        tk_yyyyyy_yyyyyz[i] = -10.0 * ts_yyyy_yyyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyyyz[i] * fe_0 + 5.0 * tk_yyyyy_yyyyz[i] * fe_0 +
                              tk_yyyyy_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyyyz[i] * fz_0;

        tk_yyyyyy_yyyyzz[i] = -10.0 * ts_yyyy_yyyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyyzz[i] * fe_0 + 4.0 * tk_yyyyy_yyyzz[i] * fe_0 +
                              tk_yyyyy_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyyzz[i] * fz_0;

        tk_yyyyyy_yyyzzz[i] = -10.0 * ts_yyyy_yyyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyzzz[i] * fe_0 + 3.0 * tk_yyyyy_yyzzz[i] * fe_0 +
                              tk_yyyyy_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyzzz[i] * fz_0;

        tk_yyyyyy_yyzzzz[i] = -10.0 * ts_yyyy_yyzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyzzzz[i] * fe_0 + 2.0 * tk_yyyyy_yzzzz[i] * fe_0 +
                              tk_yyyyy_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyzzzz[i] * fz_0;

        tk_yyyyyy_yzzzzz[i] = -10.0 * ts_yyyy_yzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yzzzzz[i] * fe_0 + tk_yyyyy_zzzzz[i] * fe_0 +
                              tk_yyyyy_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yzzzzz[i] * fz_0;

        tk_yyyyyy_zzzzzz[i] = -10.0 * ts_yyyy_zzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_zzzzzz[i] * fe_0 + tk_yyyyy_zzzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyyy_zzzzzz[i] * fz_0;
    }

    // Set up 616-644 components of targeted buffer : II

    auto tk_yyyyyz_xxxxxx = pbuffer.data(idx_kin_ii + 616);

    auto tk_yyyyyz_xxxxxy = pbuffer.data(idx_kin_ii + 617);

    auto tk_yyyyyz_xxxxxz = pbuffer.data(idx_kin_ii + 618);

    auto tk_yyyyyz_xxxxyy = pbuffer.data(idx_kin_ii + 619);

    auto tk_yyyyyz_xxxxyz = pbuffer.data(idx_kin_ii + 620);

    auto tk_yyyyyz_xxxxzz = pbuffer.data(idx_kin_ii + 621);

    auto tk_yyyyyz_xxxyyy = pbuffer.data(idx_kin_ii + 622);

    auto tk_yyyyyz_xxxyyz = pbuffer.data(idx_kin_ii + 623);

    auto tk_yyyyyz_xxxyzz = pbuffer.data(idx_kin_ii + 624);

    auto tk_yyyyyz_xxxzzz = pbuffer.data(idx_kin_ii + 625);

    auto tk_yyyyyz_xxyyyy = pbuffer.data(idx_kin_ii + 626);

    auto tk_yyyyyz_xxyyyz = pbuffer.data(idx_kin_ii + 627);

    auto tk_yyyyyz_xxyyzz = pbuffer.data(idx_kin_ii + 628);

    auto tk_yyyyyz_xxyzzz = pbuffer.data(idx_kin_ii + 629);

    auto tk_yyyyyz_xxzzzz = pbuffer.data(idx_kin_ii + 630);

    auto tk_yyyyyz_xyyyyy = pbuffer.data(idx_kin_ii + 631);

    auto tk_yyyyyz_xyyyyz = pbuffer.data(idx_kin_ii + 632);

    auto tk_yyyyyz_xyyyzz = pbuffer.data(idx_kin_ii + 633);

    auto tk_yyyyyz_xyyzzz = pbuffer.data(idx_kin_ii + 634);

    auto tk_yyyyyz_xyzzzz = pbuffer.data(idx_kin_ii + 635);

    auto tk_yyyyyz_xzzzzz = pbuffer.data(idx_kin_ii + 636);

    auto tk_yyyyyz_yyyyyy = pbuffer.data(idx_kin_ii + 637);

    auto tk_yyyyyz_yyyyyz = pbuffer.data(idx_kin_ii + 638);

    auto tk_yyyyyz_yyyyzz = pbuffer.data(idx_kin_ii + 639);

    auto tk_yyyyyz_yyyzzz = pbuffer.data(idx_kin_ii + 640);

    auto tk_yyyyyz_yyzzzz = pbuffer.data(idx_kin_ii + 641);

    auto tk_yyyyyz_yzzzzz = pbuffer.data(idx_kin_ii + 642);

    auto tk_yyyyyz_zzzzzz = pbuffer.data(idx_kin_ii + 643);

#pragma omp simd aligned(pa_z,                 \
                             tk_yyyyy_xxxxx,   \
                             tk_yyyyy_xxxxxx,  \
                             tk_yyyyy_xxxxxy,  \
                             tk_yyyyy_xxxxxz,  \
                             tk_yyyyy_xxxxy,   \
                             tk_yyyyy_xxxxyy,  \
                             tk_yyyyy_xxxxyz,  \
                             tk_yyyyy_xxxxz,   \
                             tk_yyyyy_xxxxzz,  \
                             tk_yyyyy_xxxyy,   \
                             tk_yyyyy_xxxyyy,  \
                             tk_yyyyy_xxxyyz,  \
                             tk_yyyyy_xxxyz,   \
                             tk_yyyyy_xxxyzz,  \
                             tk_yyyyy_xxxzz,   \
                             tk_yyyyy_xxxzzz,  \
                             tk_yyyyy_xxyyy,   \
                             tk_yyyyy_xxyyyy,  \
                             tk_yyyyy_xxyyyz,  \
                             tk_yyyyy_xxyyz,   \
                             tk_yyyyy_xxyyzz,  \
                             tk_yyyyy_xxyzz,   \
                             tk_yyyyy_xxyzzz,  \
                             tk_yyyyy_xxzzz,   \
                             tk_yyyyy_xxzzzz,  \
                             tk_yyyyy_xyyyy,   \
                             tk_yyyyy_xyyyyy,  \
                             tk_yyyyy_xyyyyz,  \
                             tk_yyyyy_xyyyz,   \
                             tk_yyyyy_xyyyzz,  \
                             tk_yyyyy_xyyzz,   \
                             tk_yyyyy_xyyzzz,  \
                             tk_yyyyy_xyzzz,   \
                             tk_yyyyy_xyzzzz,  \
                             tk_yyyyy_xzzzz,   \
                             tk_yyyyy_xzzzzz,  \
                             tk_yyyyy_yyyyy,   \
                             tk_yyyyy_yyyyyy,  \
                             tk_yyyyy_yyyyyz,  \
                             tk_yyyyy_yyyyz,   \
                             tk_yyyyy_yyyyzz,  \
                             tk_yyyyy_yyyzz,   \
                             tk_yyyyy_yyyzzz,  \
                             tk_yyyyy_yyzzz,   \
                             tk_yyyyy_yyzzzz,  \
                             tk_yyyyy_yzzzz,   \
                             tk_yyyyy_yzzzzz,  \
                             tk_yyyyy_zzzzz,   \
                             tk_yyyyy_zzzzzz,  \
                             tk_yyyyyz_xxxxxx, \
                             tk_yyyyyz_xxxxxy, \
                             tk_yyyyyz_xxxxxz, \
                             tk_yyyyyz_xxxxyy, \
                             tk_yyyyyz_xxxxyz, \
                             tk_yyyyyz_xxxxzz, \
                             tk_yyyyyz_xxxyyy, \
                             tk_yyyyyz_xxxyyz, \
                             tk_yyyyyz_xxxyzz, \
                             tk_yyyyyz_xxxzzz, \
                             tk_yyyyyz_xxyyyy, \
                             tk_yyyyyz_xxyyyz, \
                             tk_yyyyyz_xxyyzz, \
                             tk_yyyyyz_xxyzzz, \
                             tk_yyyyyz_xxzzzz, \
                             tk_yyyyyz_xyyyyy, \
                             tk_yyyyyz_xyyyyz, \
                             tk_yyyyyz_xyyyzz, \
                             tk_yyyyyz_xyyzzz, \
                             tk_yyyyyz_xyzzzz, \
                             tk_yyyyyz_xzzzzz, \
                             tk_yyyyyz_yyyyyy, \
                             tk_yyyyyz_yyyyyz, \
                             tk_yyyyyz_yyyyzz, \
                             tk_yyyyyz_yyyzzz, \
                             tk_yyyyyz_yyzzzz, \
                             tk_yyyyyz_yzzzzz, \
                             tk_yyyyyz_zzzzzz, \
                             ts_yyyyyz_xxxxxx, \
                             ts_yyyyyz_xxxxxy, \
                             ts_yyyyyz_xxxxxz, \
                             ts_yyyyyz_xxxxyy, \
                             ts_yyyyyz_xxxxyz, \
                             ts_yyyyyz_xxxxzz, \
                             ts_yyyyyz_xxxyyy, \
                             ts_yyyyyz_xxxyyz, \
                             ts_yyyyyz_xxxyzz, \
                             ts_yyyyyz_xxxzzz, \
                             ts_yyyyyz_xxyyyy, \
                             ts_yyyyyz_xxyyyz, \
                             ts_yyyyyz_xxyyzz, \
                             ts_yyyyyz_xxyzzz, \
                             ts_yyyyyz_xxzzzz, \
                             ts_yyyyyz_xyyyyy, \
                             ts_yyyyyz_xyyyyz, \
                             ts_yyyyyz_xyyyzz, \
                             ts_yyyyyz_xyyzzz, \
                             ts_yyyyyz_xyzzzz, \
                             ts_yyyyyz_xzzzzz, \
                             ts_yyyyyz_yyyyyy, \
                             ts_yyyyyz_yyyyyz, \
                             ts_yyyyyz_yyyyzz, \
                             ts_yyyyyz_yyyzzz, \
                             ts_yyyyyz_yyzzzz, \
                             ts_yyyyyz_yzzzzz, \
                             ts_yyyyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyyz_xxxxxx[i] = tk_yyyyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxxx[i] * fz_0;

        tk_yyyyyz_xxxxxy[i] = tk_yyyyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxxy[i] * fz_0;

        tk_yyyyyz_xxxxxz[i] = tk_yyyyy_xxxxx[i] * fe_0 + tk_yyyyy_xxxxxz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxxz[i] * fz_0;

        tk_yyyyyz_xxxxyy[i] = tk_yyyyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxyy[i] * fz_0;

        tk_yyyyyz_xxxxyz[i] = tk_yyyyy_xxxxy[i] * fe_0 + tk_yyyyy_xxxxyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxyz[i] * fz_0;

        tk_yyyyyz_xxxxzz[i] = 2.0 * tk_yyyyy_xxxxz[i] * fe_0 + tk_yyyyy_xxxxzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxzz[i] * fz_0;

        tk_yyyyyz_xxxyyy[i] = tk_yyyyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxyyy[i] * fz_0;

        tk_yyyyyz_xxxyyz[i] = tk_yyyyy_xxxyy[i] * fe_0 + tk_yyyyy_xxxyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxyyz[i] * fz_0;

        tk_yyyyyz_xxxyzz[i] = 2.0 * tk_yyyyy_xxxyz[i] * fe_0 + tk_yyyyy_xxxyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxyzz[i] * fz_0;

        tk_yyyyyz_xxxzzz[i] = 3.0 * tk_yyyyy_xxxzz[i] * fe_0 + tk_yyyyy_xxxzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxzzz[i] * fz_0;

        tk_yyyyyz_xxyyyy[i] = tk_yyyyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyyyy[i] * fz_0;

        tk_yyyyyz_xxyyyz[i] = tk_yyyyy_xxyyy[i] * fe_0 + tk_yyyyy_xxyyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyyyz[i] * fz_0;

        tk_yyyyyz_xxyyzz[i] = 2.0 * tk_yyyyy_xxyyz[i] * fe_0 + tk_yyyyy_xxyyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyyzz[i] * fz_0;

        tk_yyyyyz_xxyzzz[i] = 3.0 * tk_yyyyy_xxyzz[i] * fe_0 + tk_yyyyy_xxyzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyzzz[i] * fz_0;

        tk_yyyyyz_xxzzzz[i] = 4.0 * tk_yyyyy_xxzzz[i] * fe_0 + tk_yyyyy_xxzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxzzzz[i] * fz_0;

        tk_yyyyyz_xyyyyy[i] = tk_yyyyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyyyy[i] * fz_0;

        tk_yyyyyz_xyyyyz[i] = tk_yyyyy_xyyyy[i] * fe_0 + tk_yyyyy_xyyyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyyyz[i] * fz_0;

        tk_yyyyyz_xyyyzz[i] = 2.0 * tk_yyyyy_xyyyz[i] * fe_0 + tk_yyyyy_xyyyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyyzz[i] * fz_0;

        tk_yyyyyz_xyyzzz[i] = 3.0 * tk_yyyyy_xyyzz[i] * fe_0 + tk_yyyyy_xyyzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyzzz[i] * fz_0;

        tk_yyyyyz_xyzzzz[i] = 4.0 * tk_yyyyy_xyzzz[i] * fe_0 + tk_yyyyy_xyzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyzzzz[i] * fz_0;

        tk_yyyyyz_xzzzzz[i] = 5.0 * tk_yyyyy_xzzzz[i] * fe_0 + tk_yyyyy_xzzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xzzzzz[i] * fz_0;

        tk_yyyyyz_yyyyyy[i] = tk_yyyyy_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyyyy[i] * fz_0;

        tk_yyyyyz_yyyyyz[i] = tk_yyyyy_yyyyy[i] * fe_0 + tk_yyyyy_yyyyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyyyz[i] * fz_0;

        tk_yyyyyz_yyyyzz[i] = 2.0 * tk_yyyyy_yyyyz[i] * fe_0 + tk_yyyyy_yyyyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyyzz[i] * fz_0;

        tk_yyyyyz_yyyzzz[i] = 3.0 * tk_yyyyy_yyyzz[i] * fe_0 + tk_yyyyy_yyyzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyzzz[i] * fz_0;

        tk_yyyyyz_yyzzzz[i] = 4.0 * tk_yyyyy_yyzzz[i] * fe_0 + tk_yyyyy_yyzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyzzzz[i] * fz_0;

        tk_yyyyyz_yzzzzz[i] = 5.0 * tk_yyyyy_yzzzz[i] * fe_0 + tk_yyyyy_yzzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yzzzzz[i] * fz_0;

        tk_yyyyyz_zzzzzz[i] = 6.0 * tk_yyyyy_zzzzz[i] * fe_0 + tk_yyyyy_zzzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_zzzzzz[i] * fz_0;
    }

    // Set up 644-672 components of targeted buffer : II

    auto tk_yyyyzz_xxxxxx = pbuffer.data(idx_kin_ii + 644);

    auto tk_yyyyzz_xxxxxy = pbuffer.data(idx_kin_ii + 645);

    auto tk_yyyyzz_xxxxxz = pbuffer.data(idx_kin_ii + 646);

    auto tk_yyyyzz_xxxxyy = pbuffer.data(idx_kin_ii + 647);

    auto tk_yyyyzz_xxxxyz = pbuffer.data(idx_kin_ii + 648);

    auto tk_yyyyzz_xxxxzz = pbuffer.data(idx_kin_ii + 649);

    auto tk_yyyyzz_xxxyyy = pbuffer.data(idx_kin_ii + 650);

    auto tk_yyyyzz_xxxyyz = pbuffer.data(idx_kin_ii + 651);

    auto tk_yyyyzz_xxxyzz = pbuffer.data(idx_kin_ii + 652);

    auto tk_yyyyzz_xxxzzz = pbuffer.data(idx_kin_ii + 653);

    auto tk_yyyyzz_xxyyyy = pbuffer.data(idx_kin_ii + 654);

    auto tk_yyyyzz_xxyyyz = pbuffer.data(idx_kin_ii + 655);

    auto tk_yyyyzz_xxyyzz = pbuffer.data(idx_kin_ii + 656);

    auto tk_yyyyzz_xxyzzz = pbuffer.data(idx_kin_ii + 657);

    auto tk_yyyyzz_xxzzzz = pbuffer.data(idx_kin_ii + 658);

    auto tk_yyyyzz_xyyyyy = pbuffer.data(idx_kin_ii + 659);

    auto tk_yyyyzz_xyyyyz = pbuffer.data(idx_kin_ii + 660);

    auto tk_yyyyzz_xyyyzz = pbuffer.data(idx_kin_ii + 661);

    auto tk_yyyyzz_xyyzzz = pbuffer.data(idx_kin_ii + 662);

    auto tk_yyyyzz_xyzzzz = pbuffer.data(idx_kin_ii + 663);

    auto tk_yyyyzz_xzzzzz = pbuffer.data(idx_kin_ii + 664);

    auto tk_yyyyzz_yyyyyy = pbuffer.data(idx_kin_ii + 665);

    auto tk_yyyyzz_yyyyyz = pbuffer.data(idx_kin_ii + 666);

    auto tk_yyyyzz_yyyyzz = pbuffer.data(idx_kin_ii + 667);

    auto tk_yyyyzz_yyyzzz = pbuffer.data(idx_kin_ii + 668);

    auto tk_yyyyzz_yyzzzz = pbuffer.data(idx_kin_ii + 669);

    auto tk_yyyyzz_yzzzzz = pbuffer.data(idx_kin_ii + 670);

    auto tk_yyyyzz_zzzzzz = pbuffer.data(idx_kin_ii + 671);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tk_yyyy_xxxxxy,   \
                             tk_yyyy_xxxxyy,   \
                             tk_yyyy_xxxyyy,   \
                             tk_yyyy_xxyyyy,   \
                             tk_yyyy_xyyyyy,   \
                             tk_yyyy_yyyyyy,   \
                             tk_yyyyz_xxxxxy,  \
                             tk_yyyyz_xxxxyy,  \
                             tk_yyyyz_xxxyyy,  \
                             tk_yyyyz_xxyyyy,  \
                             tk_yyyyz_xyyyyy,  \
                             tk_yyyyz_yyyyyy,  \
                             tk_yyyyzz_xxxxxx, \
                             tk_yyyyzz_xxxxxy, \
                             tk_yyyyzz_xxxxxz, \
                             tk_yyyyzz_xxxxyy, \
                             tk_yyyyzz_xxxxyz, \
                             tk_yyyyzz_xxxxzz, \
                             tk_yyyyzz_xxxyyy, \
                             tk_yyyyzz_xxxyyz, \
                             tk_yyyyzz_xxxyzz, \
                             tk_yyyyzz_xxxzzz, \
                             tk_yyyyzz_xxyyyy, \
                             tk_yyyyzz_xxyyyz, \
                             tk_yyyyzz_xxyyzz, \
                             tk_yyyyzz_xxyzzz, \
                             tk_yyyyzz_xxzzzz, \
                             tk_yyyyzz_xyyyyy, \
                             tk_yyyyzz_xyyyyz, \
                             tk_yyyyzz_xyyyzz, \
                             tk_yyyyzz_xyyzzz, \
                             tk_yyyyzz_xyzzzz, \
                             tk_yyyyzz_xzzzzz, \
                             tk_yyyyzz_yyyyyy, \
                             tk_yyyyzz_yyyyyz, \
                             tk_yyyyzz_yyyyzz, \
                             tk_yyyyzz_yyyzzz, \
                             tk_yyyyzz_yyzzzz, \
                             tk_yyyyzz_yzzzzz, \
                             tk_yyyyzz_zzzzzz, \
                             tk_yyyzz_xxxxxx,  \
                             tk_yyyzz_xxxxxz,  \
                             tk_yyyzz_xxxxyz,  \
                             tk_yyyzz_xxxxz,   \
                             tk_yyyzz_xxxxzz,  \
                             tk_yyyzz_xxxyyz,  \
                             tk_yyyzz_xxxyz,   \
                             tk_yyyzz_xxxyzz,  \
                             tk_yyyzz_xxxzz,   \
                             tk_yyyzz_xxxzzz,  \
                             tk_yyyzz_xxyyyz,  \
                             tk_yyyzz_xxyyz,   \
                             tk_yyyzz_xxyyzz,  \
                             tk_yyyzz_xxyzz,   \
                             tk_yyyzz_xxyzzz,  \
                             tk_yyyzz_xxzzz,   \
                             tk_yyyzz_xxzzzz,  \
                             tk_yyyzz_xyyyyz,  \
                             tk_yyyzz_xyyyz,   \
                             tk_yyyzz_xyyyzz,  \
                             tk_yyyzz_xyyzz,   \
                             tk_yyyzz_xyyzzz,  \
                             tk_yyyzz_xyzzz,   \
                             tk_yyyzz_xyzzzz,  \
                             tk_yyyzz_xzzzz,   \
                             tk_yyyzz_xzzzzz,  \
                             tk_yyyzz_yyyyyz,  \
                             tk_yyyzz_yyyyz,   \
                             tk_yyyzz_yyyyzz,  \
                             tk_yyyzz_yyyzz,   \
                             tk_yyyzz_yyyzzz,  \
                             tk_yyyzz_yyzzz,   \
                             tk_yyyzz_yyzzzz,  \
                             tk_yyyzz_yzzzz,   \
                             tk_yyyzz_yzzzzz,  \
                             tk_yyyzz_zzzzz,   \
                             tk_yyyzz_zzzzzz,  \
                             tk_yyzz_xxxxxx,   \
                             tk_yyzz_xxxxxz,   \
                             tk_yyzz_xxxxyz,   \
                             tk_yyzz_xxxxzz,   \
                             tk_yyzz_xxxyyz,   \
                             tk_yyzz_xxxyzz,   \
                             tk_yyzz_xxxzzz,   \
                             tk_yyzz_xxyyyz,   \
                             tk_yyzz_xxyyzz,   \
                             tk_yyzz_xxyzzz,   \
                             tk_yyzz_xxzzzz,   \
                             tk_yyzz_xyyyyz,   \
                             tk_yyzz_xyyyzz,   \
                             tk_yyzz_xyyzzz,   \
                             tk_yyzz_xyzzzz,   \
                             tk_yyzz_xzzzzz,   \
                             tk_yyzz_yyyyyz,   \
                             tk_yyzz_yyyyzz,   \
                             tk_yyzz_yyyzzz,   \
                             tk_yyzz_yyzzzz,   \
                             tk_yyzz_yzzzzz,   \
                             tk_yyzz_zzzzzz,   \
                             ts_yyyy_xxxxxy,   \
                             ts_yyyy_xxxxyy,   \
                             ts_yyyy_xxxyyy,   \
                             ts_yyyy_xxyyyy,   \
                             ts_yyyy_xyyyyy,   \
                             ts_yyyy_yyyyyy,   \
                             ts_yyyyzz_xxxxxx, \
                             ts_yyyyzz_xxxxxy, \
                             ts_yyyyzz_xxxxxz, \
                             ts_yyyyzz_xxxxyy, \
                             ts_yyyyzz_xxxxyz, \
                             ts_yyyyzz_xxxxzz, \
                             ts_yyyyzz_xxxyyy, \
                             ts_yyyyzz_xxxyyz, \
                             ts_yyyyzz_xxxyzz, \
                             ts_yyyyzz_xxxzzz, \
                             ts_yyyyzz_xxyyyy, \
                             ts_yyyyzz_xxyyyz, \
                             ts_yyyyzz_xxyyzz, \
                             ts_yyyyzz_xxyzzz, \
                             ts_yyyyzz_xxzzzz, \
                             ts_yyyyzz_xyyyyy, \
                             ts_yyyyzz_xyyyyz, \
                             ts_yyyyzz_xyyyzz, \
                             ts_yyyyzz_xyyzzz, \
                             ts_yyyyzz_xyzzzz, \
                             ts_yyyyzz_xzzzzz, \
                             ts_yyyyzz_yyyyyy, \
                             ts_yyyyzz_yyyyyz, \
                             ts_yyyyzz_yyyyzz, \
                             ts_yyyyzz_yyyzzz, \
                             ts_yyyyzz_yyzzzz, \
                             ts_yyyyzz_yzzzzz, \
                             ts_yyyyzz_zzzzzz, \
                             ts_yyzz_xxxxxx,   \
                             ts_yyzz_xxxxxz,   \
                             ts_yyzz_xxxxyz,   \
                             ts_yyzz_xxxxzz,   \
                             ts_yyzz_xxxyyz,   \
                             ts_yyzz_xxxyzz,   \
                             ts_yyzz_xxxzzz,   \
                             ts_yyzz_xxyyyz,   \
                             ts_yyzz_xxyyzz,   \
                             ts_yyzz_xxyzzz,   \
                             ts_yyzz_xxzzzz,   \
                             ts_yyzz_xyyyyz,   \
                             ts_yyzz_xyyyzz,   \
                             ts_yyzz_xyyzzz,   \
                             ts_yyzz_xyzzzz,   \
                             ts_yyzz_xzzzzz,   \
                             ts_yyzz_yyyyyz,   \
                             ts_yyzz_yyyyzz,   \
                             ts_yyzz_yyyzzz,   \
                             ts_yyzz_yyzzzz,   \
                             ts_yyzz_yzzzzz,   \
                             ts_yyzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyzz_xxxxxx[i] = -6.0 * ts_yyzz_xxxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxxx[i] * fe_0 + tk_yyyzz_xxxxxx[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_xxxxxx[i] * fz_0;

        tk_yyyyzz_xxxxxy[i] =
            -2.0 * ts_yyyy_xxxxxy[i] * fbe_0 * fz_0 + tk_yyyy_xxxxxy[i] * fe_0 + tk_yyyyz_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxxxxy[i] * fz_0;

        tk_yyyyzz_xxxxxz[i] = -6.0 * ts_yyzz_xxxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxxz[i] * fe_0 + tk_yyyzz_xxxxxz[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_xxxxxz[i] * fz_0;

        tk_yyyyzz_xxxxyy[i] =
            -2.0 * ts_yyyy_xxxxyy[i] * fbe_0 * fz_0 + tk_yyyy_xxxxyy[i] * fe_0 + tk_yyyyz_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxxxyy[i] * fz_0;

        tk_yyyyzz_xxxxyz[i] = -6.0 * ts_yyzz_xxxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxyz[i] * fe_0 + tk_yyyzz_xxxxz[i] * fe_0 +
                              tk_yyyzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxxyz[i] * fz_0;

        tk_yyyyzz_xxxxzz[i] = -6.0 * ts_yyzz_xxxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxzz[i] * fe_0 + tk_yyyzz_xxxxzz[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_xxxxzz[i] * fz_0;

        tk_yyyyzz_xxxyyy[i] =
            -2.0 * ts_yyyy_xxxyyy[i] * fbe_0 * fz_0 + tk_yyyy_xxxyyy[i] * fe_0 + tk_yyyyz_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxxyyy[i] * fz_0;

        tk_yyyyzz_xxxyyz[i] = -6.0 * ts_yyzz_xxxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxyyz[i] * fe_0 + 2.0 * tk_yyyzz_xxxyz[i] * fe_0 +
                              tk_yyyzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxyyz[i] * fz_0;

        tk_yyyyzz_xxxyzz[i] = -6.0 * ts_yyzz_xxxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxyzz[i] * fe_0 + tk_yyyzz_xxxzz[i] * fe_0 +
                              tk_yyyzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxyzz[i] * fz_0;

        tk_yyyyzz_xxxzzz[i] = -6.0 * ts_yyzz_xxxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxzzz[i] * fe_0 + tk_yyyzz_xxxzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_xxxzzz[i] * fz_0;

        tk_yyyyzz_xxyyyy[i] =
            -2.0 * ts_yyyy_xxyyyy[i] * fbe_0 * fz_0 + tk_yyyy_xxyyyy[i] * fe_0 + tk_yyyyz_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxyyyy[i] * fz_0;

        tk_yyyyzz_xxyyyz[i] = -6.0 * ts_yyzz_xxyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyyyz[i] * fe_0 + 3.0 * tk_yyyzz_xxyyz[i] * fe_0 +
                              tk_yyyzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxyyyz[i] * fz_0;

        tk_yyyyzz_xxyyzz[i] = -6.0 * ts_yyzz_xxyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyyzz[i] * fe_0 + 2.0 * tk_yyyzz_xxyzz[i] * fe_0 +
                              tk_yyyzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxyyzz[i] * fz_0;

        tk_yyyyzz_xxyzzz[i] = -6.0 * ts_yyzz_xxyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyzzz[i] * fe_0 + tk_yyyzz_xxzzz[i] * fe_0 +
                              tk_yyyzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxyzzz[i] * fz_0;

        tk_yyyyzz_xxzzzz[i] = -6.0 * ts_yyzz_xxzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxzzzz[i] * fe_0 + tk_yyyzz_xxzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_xxzzzz[i] * fz_0;

        tk_yyyyzz_xyyyyy[i] =
            -2.0 * ts_yyyy_xyyyyy[i] * fbe_0 * fz_0 + tk_yyyy_xyyyyy[i] * fe_0 + tk_yyyyz_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xyyyyy[i] * fz_0;

        tk_yyyyzz_xyyyyz[i] = -6.0 * ts_yyzz_xyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyyyz[i] * fe_0 + 4.0 * tk_yyyzz_xyyyz[i] * fe_0 +
                              tk_yyyzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyyyyz[i] * fz_0;

        tk_yyyyzz_xyyyzz[i] = -6.0 * ts_yyzz_xyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyyzz[i] * fe_0 + 3.0 * tk_yyyzz_xyyzz[i] * fe_0 +
                              tk_yyyzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyyyzz[i] * fz_0;

        tk_yyyyzz_xyyzzz[i] = -6.0 * ts_yyzz_xyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyzzz[i] * fe_0 + 2.0 * tk_yyyzz_xyzzz[i] * fe_0 +
                              tk_yyyzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyyzzz[i] * fz_0;

        tk_yyyyzz_xyzzzz[i] = -6.0 * ts_yyzz_xyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyzzzz[i] * fe_0 + tk_yyyzz_xzzzz[i] * fe_0 +
                              tk_yyyzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyzzzz[i] * fz_0;

        tk_yyyyzz_xzzzzz[i] = -6.0 * ts_yyzz_xzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xzzzzz[i] * fe_0 + tk_yyyzz_xzzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_xzzzzz[i] * fz_0;

        tk_yyyyzz_yyyyyy[i] =
            -2.0 * ts_yyyy_yyyyyy[i] * fbe_0 * fz_0 + tk_yyyy_yyyyyy[i] * fe_0 + tk_yyyyz_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_yyyyyy[i] * fz_0;

        tk_yyyyzz_yyyyyz[i] = -6.0 * ts_yyzz_yyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyyyz[i] * fe_0 + 5.0 * tk_yyyzz_yyyyz[i] * fe_0 +
                              tk_yyyzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyyyyz[i] * fz_0;

        tk_yyyyzz_yyyyzz[i] = -6.0 * ts_yyzz_yyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyyzz[i] * fe_0 + 4.0 * tk_yyyzz_yyyzz[i] * fe_0 +
                              tk_yyyzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyyyzz[i] * fz_0;

        tk_yyyyzz_yyyzzz[i] = -6.0 * ts_yyzz_yyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyzzz[i] * fe_0 + 3.0 * tk_yyyzz_yyzzz[i] * fe_0 +
                              tk_yyyzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyyzzz[i] * fz_0;

        tk_yyyyzz_yyzzzz[i] = -6.0 * ts_yyzz_yyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyzzzz[i] * fe_0 + 2.0 * tk_yyyzz_yzzzz[i] * fe_0 +
                              tk_yyyzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyzzzz[i] * fz_0;

        tk_yyyyzz_yzzzzz[i] = -6.0 * ts_yyzz_yzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yzzzzz[i] * fe_0 + tk_yyyzz_zzzzz[i] * fe_0 +
                              tk_yyyzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yzzzzz[i] * fz_0;

        tk_yyyyzz_zzzzzz[i] = -6.0 * ts_yyzz_zzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_zzzzzz[i] * fe_0 + tk_yyyzz_zzzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyyzz_zzzzzz[i] * fz_0;
    }

    // Set up 672-700 components of targeted buffer : II

    auto tk_yyyzzz_xxxxxx = pbuffer.data(idx_kin_ii + 672);

    auto tk_yyyzzz_xxxxxy = pbuffer.data(idx_kin_ii + 673);

    auto tk_yyyzzz_xxxxxz = pbuffer.data(idx_kin_ii + 674);

    auto tk_yyyzzz_xxxxyy = pbuffer.data(idx_kin_ii + 675);

    auto tk_yyyzzz_xxxxyz = pbuffer.data(idx_kin_ii + 676);

    auto tk_yyyzzz_xxxxzz = pbuffer.data(idx_kin_ii + 677);

    auto tk_yyyzzz_xxxyyy = pbuffer.data(idx_kin_ii + 678);

    auto tk_yyyzzz_xxxyyz = pbuffer.data(idx_kin_ii + 679);

    auto tk_yyyzzz_xxxyzz = pbuffer.data(idx_kin_ii + 680);

    auto tk_yyyzzz_xxxzzz = pbuffer.data(idx_kin_ii + 681);

    auto tk_yyyzzz_xxyyyy = pbuffer.data(idx_kin_ii + 682);

    auto tk_yyyzzz_xxyyyz = pbuffer.data(idx_kin_ii + 683);

    auto tk_yyyzzz_xxyyzz = pbuffer.data(idx_kin_ii + 684);

    auto tk_yyyzzz_xxyzzz = pbuffer.data(idx_kin_ii + 685);

    auto tk_yyyzzz_xxzzzz = pbuffer.data(idx_kin_ii + 686);

    auto tk_yyyzzz_xyyyyy = pbuffer.data(idx_kin_ii + 687);

    auto tk_yyyzzz_xyyyyz = pbuffer.data(idx_kin_ii + 688);

    auto tk_yyyzzz_xyyyzz = pbuffer.data(idx_kin_ii + 689);

    auto tk_yyyzzz_xyyzzz = pbuffer.data(idx_kin_ii + 690);

    auto tk_yyyzzz_xyzzzz = pbuffer.data(idx_kin_ii + 691);

    auto tk_yyyzzz_xzzzzz = pbuffer.data(idx_kin_ii + 692);

    auto tk_yyyzzz_yyyyyy = pbuffer.data(idx_kin_ii + 693);

    auto tk_yyyzzz_yyyyyz = pbuffer.data(idx_kin_ii + 694);

    auto tk_yyyzzz_yyyyzz = pbuffer.data(idx_kin_ii + 695);

    auto tk_yyyzzz_yyyzzz = pbuffer.data(idx_kin_ii + 696);

    auto tk_yyyzzz_yyzzzz = pbuffer.data(idx_kin_ii + 697);

    auto tk_yyyzzz_yzzzzz = pbuffer.data(idx_kin_ii + 698);

    auto tk_yyyzzz_zzzzzz = pbuffer.data(idx_kin_ii + 699);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tk_yyyz_xxxxxy,   \
                             tk_yyyz_xxxxyy,   \
                             tk_yyyz_xxxyyy,   \
                             tk_yyyz_xxyyyy,   \
                             tk_yyyz_xyyyyy,   \
                             tk_yyyz_yyyyyy,   \
                             tk_yyyzz_xxxxxy,  \
                             tk_yyyzz_xxxxyy,  \
                             tk_yyyzz_xxxyyy,  \
                             tk_yyyzz_xxyyyy,  \
                             tk_yyyzz_xyyyyy,  \
                             tk_yyyzz_yyyyyy,  \
                             tk_yyyzzz_xxxxxx, \
                             tk_yyyzzz_xxxxxy, \
                             tk_yyyzzz_xxxxxz, \
                             tk_yyyzzz_xxxxyy, \
                             tk_yyyzzz_xxxxyz, \
                             tk_yyyzzz_xxxxzz, \
                             tk_yyyzzz_xxxyyy, \
                             tk_yyyzzz_xxxyyz, \
                             tk_yyyzzz_xxxyzz, \
                             tk_yyyzzz_xxxzzz, \
                             tk_yyyzzz_xxyyyy, \
                             tk_yyyzzz_xxyyyz, \
                             tk_yyyzzz_xxyyzz, \
                             tk_yyyzzz_xxyzzz, \
                             tk_yyyzzz_xxzzzz, \
                             tk_yyyzzz_xyyyyy, \
                             tk_yyyzzz_xyyyyz, \
                             tk_yyyzzz_xyyyzz, \
                             tk_yyyzzz_xyyzzz, \
                             tk_yyyzzz_xyzzzz, \
                             tk_yyyzzz_xzzzzz, \
                             tk_yyyzzz_yyyyyy, \
                             tk_yyyzzz_yyyyyz, \
                             tk_yyyzzz_yyyyzz, \
                             tk_yyyzzz_yyyzzz, \
                             tk_yyyzzz_yyzzzz, \
                             tk_yyyzzz_yzzzzz, \
                             tk_yyyzzz_zzzzzz, \
                             tk_yyzzz_xxxxxx,  \
                             tk_yyzzz_xxxxxz,  \
                             tk_yyzzz_xxxxyz,  \
                             tk_yyzzz_xxxxz,   \
                             tk_yyzzz_xxxxzz,  \
                             tk_yyzzz_xxxyyz,  \
                             tk_yyzzz_xxxyz,   \
                             tk_yyzzz_xxxyzz,  \
                             tk_yyzzz_xxxzz,   \
                             tk_yyzzz_xxxzzz,  \
                             tk_yyzzz_xxyyyz,  \
                             tk_yyzzz_xxyyz,   \
                             tk_yyzzz_xxyyzz,  \
                             tk_yyzzz_xxyzz,   \
                             tk_yyzzz_xxyzzz,  \
                             tk_yyzzz_xxzzz,   \
                             tk_yyzzz_xxzzzz,  \
                             tk_yyzzz_xyyyyz,  \
                             tk_yyzzz_xyyyz,   \
                             tk_yyzzz_xyyyzz,  \
                             tk_yyzzz_xyyzz,   \
                             tk_yyzzz_xyyzzz,  \
                             tk_yyzzz_xyzzz,   \
                             tk_yyzzz_xyzzzz,  \
                             tk_yyzzz_xzzzz,   \
                             tk_yyzzz_xzzzzz,  \
                             tk_yyzzz_yyyyyz,  \
                             tk_yyzzz_yyyyz,   \
                             tk_yyzzz_yyyyzz,  \
                             tk_yyzzz_yyyzz,   \
                             tk_yyzzz_yyyzzz,  \
                             tk_yyzzz_yyzzz,   \
                             tk_yyzzz_yyzzzz,  \
                             tk_yyzzz_yzzzz,   \
                             tk_yyzzz_yzzzzz,  \
                             tk_yyzzz_zzzzz,   \
                             tk_yyzzz_zzzzzz,  \
                             tk_yzzz_xxxxxx,   \
                             tk_yzzz_xxxxxz,   \
                             tk_yzzz_xxxxyz,   \
                             tk_yzzz_xxxxzz,   \
                             tk_yzzz_xxxyyz,   \
                             tk_yzzz_xxxyzz,   \
                             tk_yzzz_xxxzzz,   \
                             tk_yzzz_xxyyyz,   \
                             tk_yzzz_xxyyzz,   \
                             tk_yzzz_xxyzzz,   \
                             tk_yzzz_xxzzzz,   \
                             tk_yzzz_xyyyyz,   \
                             tk_yzzz_xyyyzz,   \
                             tk_yzzz_xyyzzz,   \
                             tk_yzzz_xyzzzz,   \
                             tk_yzzz_xzzzzz,   \
                             tk_yzzz_yyyyyz,   \
                             tk_yzzz_yyyyzz,   \
                             tk_yzzz_yyyzzz,   \
                             tk_yzzz_yyzzzz,   \
                             tk_yzzz_yzzzzz,   \
                             tk_yzzz_zzzzzz,   \
                             ts_yyyz_xxxxxy,   \
                             ts_yyyz_xxxxyy,   \
                             ts_yyyz_xxxyyy,   \
                             ts_yyyz_xxyyyy,   \
                             ts_yyyz_xyyyyy,   \
                             ts_yyyz_yyyyyy,   \
                             ts_yyyzzz_xxxxxx, \
                             ts_yyyzzz_xxxxxy, \
                             ts_yyyzzz_xxxxxz, \
                             ts_yyyzzz_xxxxyy, \
                             ts_yyyzzz_xxxxyz, \
                             ts_yyyzzz_xxxxzz, \
                             ts_yyyzzz_xxxyyy, \
                             ts_yyyzzz_xxxyyz, \
                             ts_yyyzzz_xxxyzz, \
                             ts_yyyzzz_xxxzzz, \
                             ts_yyyzzz_xxyyyy, \
                             ts_yyyzzz_xxyyyz, \
                             ts_yyyzzz_xxyyzz, \
                             ts_yyyzzz_xxyzzz, \
                             ts_yyyzzz_xxzzzz, \
                             ts_yyyzzz_xyyyyy, \
                             ts_yyyzzz_xyyyyz, \
                             ts_yyyzzz_xyyyzz, \
                             ts_yyyzzz_xyyzzz, \
                             ts_yyyzzz_xyzzzz, \
                             ts_yyyzzz_xzzzzz, \
                             ts_yyyzzz_yyyyyy, \
                             ts_yyyzzz_yyyyyz, \
                             ts_yyyzzz_yyyyzz, \
                             ts_yyyzzz_yyyzzz, \
                             ts_yyyzzz_yyzzzz, \
                             ts_yyyzzz_yzzzzz, \
                             ts_yyyzzz_zzzzzz, \
                             ts_yzzz_xxxxxx,   \
                             ts_yzzz_xxxxxz,   \
                             ts_yzzz_xxxxyz,   \
                             ts_yzzz_xxxxzz,   \
                             ts_yzzz_xxxyyz,   \
                             ts_yzzz_xxxyzz,   \
                             ts_yzzz_xxxzzz,   \
                             ts_yzzz_xxyyyz,   \
                             ts_yzzz_xxyyzz,   \
                             ts_yzzz_xxyzzz,   \
                             ts_yzzz_xxzzzz,   \
                             ts_yzzz_xyyyyz,   \
                             ts_yzzz_xyyyzz,   \
                             ts_yzzz_xyyzzz,   \
                             ts_yzzz_xyzzzz,   \
                             ts_yzzz_xzzzzz,   \
                             ts_yzzz_yyyyyz,   \
                             ts_yzzz_yyyyzz,   \
                             ts_yzzz_yyyzzz,   \
                             ts_yzzz_yyzzzz,   \
                             ts_yzzz_yzzzzz,   \
                             ts_yzzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzzz_xxxxxx[i] = -4.0 * ts_yzzz_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxxxx[i] * fe_0 + tk_yyzzz_xxxxxx[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_xxxxxx[i] * fz_0;

        tk_yyyzzz_xxxxxy[i] = -4.0 * ts_yyyz_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxxxxy[i] * fe_0 + tk_yyyzz_xxxxxy[i] * pa_z[i] +
                              2.0 * ts_yyyzzz_xxxxxy[i] * fz_0;

        tk_yyyzzz_xxxxxz[i] = -4.0 * ts_yzzz_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxxxz[i] * fe_0 + tk_yyzzz_xxxxxz[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_xxxxxz[i] * fz_0;

        tk_yyyzzz_xxxxyy[i] = -4.0 * ts_yyyz_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxxxyy[i] * fe_0 + tk_yyyzz_xxxxyy[i] * pa_z[i] +
                              2.0 * ts_yyyzzz_xxxxyy[i] * fz_0;

        tk_yyyzzz_xxxxyz[i] = -4.0 * ts_yzzz_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxxyz[i] * fe_0 + tk_yyzzz_xxxxz[i] * fe_0 +
                              tk_yyzzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxxyz[i] * fz_0;

        tk_yyyzzz_xxxxzz[i] = -4.0 * ts_yzzz_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxxzz[i] * fe_0 + tk_yyzzz_xxxxzz[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_xxxxzz[i] * fz_0;

        tk_yyyzzz_xxxyyy[i] = -4.0 * ts_yyyz_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxxyyy[i] * fe_0 + tk_yyyzz_xxxyyy[i] * pa_z[i] +
                              2.0 * ts_yyyzzz_xxxyyy[i] * fz_0;

        tk_yyyzzz_xxxyyz[i] = -4.0 * ts_yzzz_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxyyz[i] * fe_0 + 2.0 * tk_yyzzz_xxxyz[i] * fe_0 +
                              tk_yyzzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxyyz[i] * fz_0;

        tk_yyyzzz_xxxyzz[i] = -4.0 * ts_yzzz_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxyzz[i] * fe_0 + tk_yyzzz_xxxzz[i] * fe_0 +
                              tk_yyzzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxyzz[i] * fz_0;

        tk_yyyzzz_xxxzzz[i] = -4.0 * ts_yzzz_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxzzz[i] * fe_0 + tk_yyzzz_xxxzzz[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_xxxzzz[i] * fz_0;

        tk_yyyzzz_xxyyyy[i] = -4.0 * ts_yyyz_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxyyyy[i] * fe_0 + tk_yyyzz_xxyyyy[i] * pa_z[i] +
                              2.0 * ts_yyyzzz_xxyyyy[i] * fz_0;

        tk_yyyzzz_xxyyyz[i] = -4.0 * ts_yzzz_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxyyyz[i] * fe_0 + 3.0 * tk_yyzzz_xxyyz[i] * fe_0 +
                              tk_yyzzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxyyyz[i] * fz_0;

        tk_yyyzzz_xxyyzz[i] = -4.0 * ts_yzzz_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxyyzz[i] * fe_0 + 2.0 * tk_yyzzz_xxyzz[i] * fe_0 +
                              tk_yyzzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxyyzz[i] * fz_0;

        tk_yyyzzz_xxyzzz[i] = -4.0 * ts_yzzz_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxyzzz[i] * fe_0 + tk_yyzzz_xxzzz[i] * fe_0 +
                              tk_yyzzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxyzzz[i] * fz_0;

        tk_yyyzzz_xxzzzz[i] = -4.0 * ts_yzzz_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxzzzz[i] * fe_0 + tk_yyzzz_xxzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_xxzzzz[i] * fz_0;

        tk_yyyzzz_xyyyyy[i] = -4.0 * ts_yyyz_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xyyyyy[i] * fe_0 + tk_yyyzz_xyyyyy[i] * pa_z[i] +
                              2.0 * ts_yyyzzz_xyyyyy[i] * fz_0;

        tk_yyyzzz_xyyyyz[i] = -4.0 * ts_yzzz_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyyyyz[i] * fe_0 + 4.0 * tk_yyzzz_xyyyz[i] * fe_0 +
                              tk_yyzzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyyyyz[i] * fz_0;

        tk_yyyzzz_xyyyzz[i] = -4.0 * ts_yzzz_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyyyzz[i] * fe_0 + 3.0 * tk_yyzzz_xyyzz[i] * fe_0 +
                              tk_yyzzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyyyzz[i] * fz_0;

        tk_yyyzzz_xyyzzz[i] = -4.0 * ts_yzzz_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyyzzz[i] * fe_0 + 2.0 * tk_yyzzz_xyzzz[i] * fe_0 +
                              tk_yyzzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyyzzz[i] * fz_0;

        tk_yyyzzz_xyzzzz[i] = -4.0 * ts_yzzz_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyzzzz[i] * fe_0 + tk_yyzzz_xzzzz[i] * fe_0 +
                              tk_yyzzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyzzzz[i] * fz_0;

        tk_yyyzzz_xzzzzz[i] = -4.0 * ts_yzzz_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xzzzzz[i] * fe_0 + tk_yyzzz_xzzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_xzzzzz[i] * fz_0;

        tk_yyyzzz_yyyyyy[i] = -4.0 * ts_yyyz_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_yyyyyy[i] * fe_0 + tk_yyyzz_yyyyyy[i] * pa_z[i] +
                              2.0 * ts_yyyzzz_yyyyyy[i] * fz_0;

        tk_yyyzzz_yyyyyz[i] = -4.0 * ts_yzzz_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyyyyz[i] * fe_0 + 5.0 * tk_yyzzz_yyyyz[i] * fe_0 +
                              tk_yyzzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyyyyz[i] * fz_0;

        tk_yyyzzz_yyyyzz[i] = -4.0 * ts_yzzz_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyyyzz[i] * fe_0 + 4.0 * tk_yyzzz_yyyzz[i] * fe_0 +
                              tk_yyzzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyyyzz[i] * fz_0;

        tk_yyyzzz_yyyzzz[i] = -4.0 * ts_yzzz_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyyzzz[i] * fe_0 + 3.0 * tk_yyzzz_yyzzz[i] * fe_0 +
                              tk_yyzzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyyzzz[i] * fz_0;

        tk_yyyzzz_yyzzzz[i] = -4.0 * ts_yzzz_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyzzzz[i] * fe_0 + 2.0 * tk_yyzzz_yzzzz[i] * fe_0 +
                              tk_yyzzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyzzzz[i] * fz_0;

        tk_yyyzzz_yzzzzz[i] = -4.0 * ts_yzzz_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yzzzzz[i] * fe_0 + tk_yyzzz_zzzzz[i] * fe_0 +
                              tk_yyzzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yzzzzz[i] * fz_0;

        tk_yyyzzz_zzzzzz[i] = -4.0 * ts_yzzz_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_zzzzzz[i] * fe_0 + tk_yyzzz_zzzzzz[i] * pa_y[i] +
                              2.0 * ts_yyyzzz_zzzzzz[i] * fz_0;
    }

    // Set up 700-728 components of targeted buffer : II

    auto tk_yyzzzz_xxxxxx = pbuffer.data(idx_kin_ii + 700);

    auto tk_yyzzzz_xxxxxy = pbuffer.data(idx_kin_ii + 701);

    auto tk_yyzzzz_xxxxxz = pbuffer.data(idx_kin_ii + 702);

    auto tk_yyzzzz_xxxxyy = pbuffer.data(idx_kin_ii + 703);

    auto tk_yyzzzz_xxxxyz = pbuffer.data(idx_kin_ii + 704);

    auto tk_yyzzzz_xxxxzz = pbuffer.data(idx_kin_ii + 705);

    auto tk_yyzzzz_xxxyyy = pbuffer.data(idx_kin_ii + 706);

    auto tk_yyzzzz_xxxyyz = pbuffer.data(idx_kin_ii + 707);

    auto tk_yyzzzz_xxxyzz = pbuffer.data(idx_kin_ii + 708);

    auto tk_yyzzzz_xxxzzz = pbuffer.data(idx_kin_ii + 709);

    auto tk_yyzzzz_xxyyyy = pbuffer.data(idx_kin_ii + 710);

    auto tk_yyzzzz_xxyyyz = pbuffer.data(idx_kin_ii + 711);

    auto tk_yyzzzz_xxyyzz = pbuffer.data(idx_kin_ii + 712);

    auto tk_yyzzzz_xxyzzz = pbuffer.data(idx_kin_ii + 713);

    auto tk_yyzzzz_xxzzzz = pbuffer.data(idx_kin_ii + 714);

    auto tk_yyzzzz_xyyyyy = pbuffer.data(idx_kin_ii + 715);

    auto tk_yyzzzz_xyyyyz = pbuffer.data(idx_kin_ii + 716);

    auto tk_yyzzzz_xyyyzz = pbuffer.data(idx_kin_ii + 717);

    auto tk_yyzzzz_xyyzzz = pbuffer.data(idx_kin_ii + 718);

    auto tk_yyzzzz_xyzzzz = pbuffer.data(idx_kin_ii + 719);

    auto tk_yyzzzz_xzzzzz = pbuffer.data(idx_kin_ii + 720);

    auto tk_yyzzzz_yyyyyy = pbuffer.data(idx_kin_ii + 721);

    auto tk_yyzzzz_yyyyyz = pbuffer.data(idx_kin_ii + 722);

    auto tk_yyzzzz_yyyyzz = pbuffer.data(idx_kin_ii + 723);

    auto tk_yyzzzz_yyyzzz = pbuffer.data(idx_kin_ii + 724);

    auto tk_yyzzzz_yyzzzz = pbuffer.data(idx_kin_ii + 725);

    auto tk_yyzzzz_yzzzzz = pbuffer.data(idx_kin_ii + 726);

    auto tk_yyzzzz_zzzzzz = pbuffer.data(idx_kin_ii + 727);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tk_yyzz_xxxxxy,   \
                             tk_yyzz_xxxxyy,   \
                             tk_yyzz_xxxyyy,   \
                             tk_yyzz_xxyyyy,   \
                             tk_yyzz_xyyyyy,   \
                             tk_yyzz_yyyyyy,   \
                             tk_yyzzz_xxxxxy,  \
                             tk_yyzzz_xxxxyy,  \
                             tk_yyzzz_xxxyyy,  \
                             tk_yyzzz_xxyyyy,  \
                             tk_yyzzz_xyyyyy,  \
                             tk_yyzzz_yyyyyy,  \
                             tk_yyzzzz_xxxxxx, \
                             tk_yyzzzz_xxxxxy, \
                             tk_yyzzzz_xxxxxz, \
                             tk_yyzzzz_xxxxyy, \
                             tk_yyzzzz_xxxxyz, \
                             tk_yyzzzz_xxxxzz, \
                             tk_yyzzzz_xxxyyy, \
                             tk_yyzzzz_xxxyyz, \
                             tk_yyzzzz_xxxyzz, \
                             tk_yyzzzz_xxxzzz, \
                             tk_yyzzzz_xxyyyy, \
                             tk_yyzzzz_xxyyyz, \
                             tk_yyzzzz_xxyyzz, \
                             tk_yyzzzz_xxyzzz, \
                             tk_yyzzzz_xxzzzz, \
                             tk_yyzzzz_xyyyyy, \
                             tk_yyzzzz_xyyyyz, \
                             tk_yyzzzz_xyyyzz, \
                             tk_yyzzzz_xyyzzz, \
                             tk_yyzzzz_xyzzzz, \
                             tk_yyzzzz_xzzzzz, \
                             tk_yyzzzz_yyyyyy, \
                             tk_yyzzzz_yyyyyz, \
                             tk_yyzzzz_yyyyzz, \
                             tk_yyzzzz_yyyzzz, \
                             tk_yyzzzz_yyzzzz, \
                             tk_yyzzzz_yzzzzz, \
                             tk_yyzzzz_zzzzzz, \
                             tk_yzzzz_xxxxxx,  \
                             tk_yzzzz_xxxxxz,  \
                             tk_yzzzz_xxxxyz,  \
                             tk_yzzzz_xxxxz,   \
                             tk_yzzzz_xxxxzz,  \
                             tk_yzzzz_xxxyyz,  \
                             tk_yzzzz_xxxyz,   \
                             tk_yzzzz_xxxyzz,  \
                             tk_yzzzz_xxxzz,   \
                             tk_yzzzz_xxxzzz,  \
                             tk_yzzzz_xxyyyz,  \
                             tk_yzzzz_xxyyz,   \
                             tk_yzzzz_xxyyzz,  \
                             tk_yzzzz_xxyzz,   \
                             tk_yzzzz_xxyzzz,  \
                             tk_yzzzz_xxzzz,   \
                             tk_yzzzz_xxzzzz,  \
                             tk_yzzzz_xyyyyz,  \
                             tk_yzzzz_xyyyz,   \
                             tk_yzzzz_xyyyzz,  \
                             tk_yzzzz_xyyzz,   \
                             tk_yzzzz_xyyzzz,  \
                             tk_yzzzz_xyzzz,   \
                             tk_yzzzz_xyzzzz,  \
                             tk_yzzzz_xzzzz,   \
                             tk_yzzzz_xzzzzz,  \
                             tk_yzzzz_yyyyyz,  \
                             tk_yzzzz_yyyyz,   \
                             tk_yzzzz_yyyyzz,  \
                             tk_yzzzz_yyyzz,   \
                             tk_yzzzz_yyyzzz,  \
                             tk_yzzzz_yyzzz,   \
                             tk_yzzzz_yyzzzz,  \
                             tk_yzzzz_yzzzz,   \
                             tk_yzzzz_yzzzzz,  \
                             tk_yzzzz_zzzzz,   \
                             tk_yzzzz_zzzzzz,  \
                             tk_zzzz_xxxxxx,   \
                             tk_zzzz_xxxxxz,   \
                             tk_zzzz_xxxxyz,   \
                             tk_zzzz_xxxxzz,   \
                             tk_zzzz_xxxyyz,   \
                             tk_zzzz_xxxyzz,   \
                             tk_zzzz_xxxzzz,   \
                             tk_zzzz_xxyyyz,   \
                             tk_zzzz_xxyyzz,   \
                             tk_zzzz_xxyzzz,   \
                             tk_zzzz_xxzzzz,   \
                             tk_zzzz_xyyyyz,   \
                             tk_zzzz_xyyyzz,   \
                             tk_zzzz_xyyzzz,   \
                             tk_zzzz_xyzzzz,   \
                             tk_zzzz_xzzzzz,   \
                             tk_zzzz_yyyyyz,   \
                             tk_zzzz_yyyyzz,   \
                             tk_zzzz_yyyzzz,   \
                             tk_zzzz_yyzzzz,   \
                             tk_zzzz_yzzzzz,   \
                             tk_zzzz_zzzzzz,   \
                             ts_yyzz_xxxxxy,   \
                             ts_yyzz_xxxxyy,   \
                             ts_yyzz_xxxyyy,   \
                             ts_yyzz_xxyyyy,   \
                             ts_yyzz_xyyyyy,   \
                             ts_yyzz_yyyyyy,   \
                             ts_yyzzzz_xxxxxx, \
                             ts_yyzzzz_xxxxxy, \
                             ts_yyzzzz_xxxxxz, \
                             ts_yyzzzz_xxxxyy, \
                             ts_yyzzzz_xxxxyz, \
                             ts_yyzzzz_xxxxzz, \
                             ts_yyzzzz_xxxyyy, \
                             ts_yyzzzz_xxxyyz, \
                             ts_yyzzzz_xxxyzz, \
                             ts_yyzzzz_xxxzzz, \
                             ts_yyzzzz_xxyyyy, \
                             ts_yyzzzz_xxyyyz, \
                             ts_yyzzzz_xxyyzz, \
                             ts_yyzzzz_xxyzzz, \
                             ts_yyzzzz_xxzzzz, \
                             ts_yyzzzz_xyyyyy, \
                             ts_yyzzzz_xyyyyz, \
                             ts_yyzzzz_xyyyzz, \
                             ts_yyzzzz_xyyzzz, \
                             ts_yyzzzz_xyzzzz, \
                             ts_yyzzzz_xzzzzz, \
                             ts_yyzzzz_yyyyyy, \
                             ts_yyzzzz_yyyyyz, \
                             ts_yyzzzz_yyyyzz, \
                             ts_yyzzzz_yyyzzz, \
                             ts_yyzzzz_yyzzzz, \
                             ts_yyzzzz_yzzzzz, \
                             ts_yyzzzz_zzzzzz, \
                             ts_zzzz_xxxxxx,   \
                             ts_zzzz_xxxxxz,   \
                             ts_zzzz_xxxxyz,   \
                             ts_zzzz_xxxxzz,   \
                             ts_zzzz_xxxyyz,   \
                             ts_zzzz_xxxyzz,   \
                             ts_zzzz_xxxzzz,   \
                             ts_zzzz_xxyyyz,   \
                             ts_zzzz_xxyyzz,   \
                             ts_zzzz_xxyzzz,   \
                             ts_zzzz_xxzzzz,   \
                             ts_zzzz_xyyyyz,   \
                             ts_zzzz_xyyyzz,   \
                             ts_zzzz_xyyzzz,   \
                             ts_zzzz_xyzzzz,   \
                             ts_zzzz_xzzzzz,   \
                             ts_zzzz_yyyyyz,   \
                             ts_zzzz_yyyyzz,   \
                             ts_zzzz_yyyzzz,   \
                             ts_zzzz_yyzzzz,   \
                             ts_zzzz_yzzzzz,   \
                             ts_zzzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzzz_xxxxxx[i] =
            -2.0 * ts_zzzz_xxxxxx[i] * fbe_0 * fz_0 + tk_zzzz_xxxxxx[i] * fe_0 + tk_yzzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxxxx[i] * fz_0;

        tk_yyzzzz_xxxxxy[i] = -6.0 * ts_yyzz_xxxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxxy[i] * fe_0 + tk_yyzzz_xxxxxy[i] * pa_z[i] +
                              2.0 * ts_yyzzzz_xxxxxy[i] * fz_0;

        tk_yyzzzz_xxxxxz[i] =
            -2.0 * ts_zzzz_xxxxxz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxxz[i] * fe_0 + tk_yzzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxxxz[i] * fz_0;

        tk_yyzzzz_xxxxyy[i] = -6.0 * ts_yyzz_xxxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxyy[i] * fe_0 + tk_yyzzz_xxxxyy[i] * pa_z[i] +
                              2.0 * ts_yyzzzz_xxxxyy[i] * fz_0;

        tk_yyzzzz_xxxxyz[i] = -2.0 * ts_zzzz_xxxxyz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxyz[i] * fe_0 + tk_yzzzz_xxxxz[i] * fe_0 +
                              tk_yzzzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxxyz[i] * fz_0;

        tk_yyzzzz_xxxxzz[i] =
            -2.0 * ts_zzzz_xxxxzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxzz[i] * fe_0 + tk_yzzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxxzz[i] * fz_0;

        tk_yyzzzz_xxxyyy[i] = -6.0 * ts_yyzz_xxxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxyyy[i] * fe_0 + tk_yyzzz_xxxyyy[i] * pa_z[i] +
                              2.0 * ts_yyzzzz_xxxyyy[i] * fz_0;

        tk_yyzzzz_xxxyyz[i] = -2.0 * ts_zzzz_xxxyyz[i] * fbe_0 * fz_0 + tk_zzzz_xxxyyz[i] * fe_0 + 2.0 * tk_yzzzz_xxxyz[i] * fe_0 +
                              tk_yzzzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxyyz[i] * fz_0;

        tk_yyzzzz_xxxyzz[i] = -2.0 * ts_zzzz_xxxyzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxyzz[i] * fe_0 + tk_yzzzz_xxxzz[i] * fe_0 +
                              tk_yzzzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxyzz[i] * fz_0;

        tk_yyzzzz_xxxzzz[i] =
            -2.0 * ts_zzzz_xxxzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxzzz[i] * fe_0 + tk_yzzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxzzz[i] * fz_0;

        tk_yyzzzz_xxyyyy[i] = -6.0 * ts_yyzz_xxyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyyyy[i] * fe_0 + tk_yyzzz_xxyyyy[i] * pa_z[i] +
                              2.0 * ts_yyzzzz_xxyyyy[i] * fz_0;

        tk_yyzzzz_xxyyyz[i] = -2.0 * ts_zzzz_xxyyyz[i] * fbe_0 * fz_0 + tk_zzzz_xxyyyz[i] * fe_0 + 3.0 * tk_yzzzz_xxyyz[i] * fe_0 +
                              tk_yzzzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxyyyz[i] * fz_0;

        tk_yyzzzz_xxyyzz[i] = -2.0 * ts_zzzz_xxyyzz[i] * fbe_0 * fz_0 + tk_zzzz_xxyyzz[i] * fe_0 + 2.0 * tk_yzzzz_xxyzz[i] * fe_0 +
                              tk_yzzzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxyyzz[i] * fz_0;

        tk_yyzzzz_xxyzzz[i] = -2.0 * ts_zzzz_xxyzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxyzzz[i] * fe_0 + tk_yzzzz_xxzzz[i] * fe_0 +
                              tk_yzzzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxyzzz[i] * fz_0;

        tk_yyzzzz_xxzzzz[i] =
            -2.0 * ts_zzzz_xxzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxzzzz[i] * fe_0 + tk_yzzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxzzzz[i] * fz_0;

        tk_yyzzzz_xyyyyy[i] = -6.0 * ts_yyzz_xyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyyyy[i] * fe_0 + tk_yyzzz_xyyyyy[i] * pa_z[i] +
                              2.0 * ts_yyzzzz_xyyyyy[i] * fz_0;

        tk_yyzzzz_xyyyyz[i] = -2.0 * ts_zzzz_xyyyyz[i] * fbe_0 * fz_0 + tk_zzzz_xyyyyz[i] * fe_0 + 4.0 * tk_yzzzz_xyyyz[i] * fe_0 +
                              tk_yzzzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyyyyz[i] * fz_0;

        tk_yyzzzz_xyyyzz[i] = -2.0 * ts_zzzz_xyyyzz[i] * fbe_0 * fz_0 + tk_zzzz_xyyyzz[i] * fe_0 + 3.0 * tk_yzzzz_xyyzz[i] * fe_0 +
                              tk_yzzzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyyyzz[i] * fz_0;

        tk_yyzzzz_xyyzzz[i] = -2.0 * ts_zzzz_xyyzzz[i] * fbe_0 * fz_0 + tk_zzzz_xyyzzz[i] * fe_0 + 2.0 * tk_yzzzz_xyzzz[i] * fe_0 +
                              tk_yzzzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyyzzz[i] * fz_0;

        tk_yyzzzz_xyzzzz[i] = -2.0 * ts_zzzz_xyzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xyzzzz[i] * fe_0 + tk_yzzzz_xzzzz[i] * fe_0 +
                              tk_yzzzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyzzzz[i] * fz_0;

        tk_yyzzzz_xzzzzz[i] =
            -2.0 * ts_zzzz_xzzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xzzzzz[i] * fe_0 + tk_yzzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xzzzzz[i] * fz_0;

        tk_yyzzzz_yyyyyy[i] = -6.0 * ts_yyzz_yyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyyyy[i] * fe_0 + tk_yyzzz_yyyyyy[i] * pa_z[i] +
                              2.0 * ts_yyzzzz_yyyyyy[i] * fz_0;

        tk_yyzzzz_yyyyyz[i] = -2.0 * ts_zzzz_yyyyyz[i] * fbe_0 * fz_0 + tk_zzzz_yyyyyz[i] * fe_0 + 5.0 * tk_yzzzz_yyyyz[i] * fe_0 +
                              tk_yzzzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyyyyz[i] * fz_0;

        tk_yyzzzz_yyyyzz[i] = -2.0 * ts_zzzz_yyyyzz[i] * fbe_0 * fz_0 + tk_zzzz_yyyyzz[i] * fe_0 + 4.0 * tk_yzzzz_yyyzz[i] * fe_0 +
                              tk_yzzzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyyyzz[i] * fz_0;

        tk_yyzzzz_yyyzzz[i] = -2.0 * ts_zzzz_yyyzzz[i] * fbe_0 * fz_0 + tk_zzzz_yyyzzz[i] * fe_0 + 3.0 * tk_yzzzz_yyzzz[i] * fe_0 +
                              tk_yzzzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyyzzz[i] * fz_0;

        tk_yyzzzz_yyzzzz[i] = -2.0 * ts_zzzz_yyzzzz[i] * fbe_0 * fz_0 + tk_zzzz_yyzzzz[i] * fe_0 + 2.0 * tk_yzzzz_yzzzz[i] * fe_0 +
                              tk_yzzzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyzzzz[i] * fz_0;

        tk_yyzzzz_yzzzzz[i] = -2.0 * ts_zzzz_yzzzzz[i] * fbe_0 * fz_0 + tk_zzzz_yzzzzz[i] * fe_0 + tk_yzzzz_zzzzz[i] * fe_0 +
                              tk_yzzzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yzzzzz[i] * fz_0;

        tk_yyzzzz_zzzzzz[i] =
            -2.0 * ts_zzzz_zzzzzz[i] * fbe_0 * fz_0 + tk_zzzz_zzzzzz[i] * fe_0 + tk_yzzzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 728-756 components of targeted buffer : II

    auto tk_yzzzzz_xxxxxx = pbuffer.data(idx_kin_ii + 728);

    auto tk_yzzzzz_xxxxxy = pbuffer.data(idx_kin_ii + 729);

    auto tk_yzzzzz_xxxxxz = pbuffer.data(idx_kin_ii + 730);

    auto tk_yzzzzz_xxxxyy = pbuffer.data(idx_kin_ii + 731);

    auto tk_yzzzzz_xxxxyz = pbuffer.data(idx_kin_ii + 732);

    auto tk_yzzzzz_xxxxzz = pbuffer.data(idx_kin_ii + 733);

    auto tk_yzzzzz_xxxyyy = pbuffer.data(idx_kin_ii + 734);

    auto tk_yzzzzz_xxxyyz = pbuffer.data(idx_kin_ii + 735);

    auto tk_yzzzzz_xxxyzz = pbuffer.data(idx_kin_ii + 736);

    auto tk_yzzzzz_xxxzzz = pbuffer.data(idx_kin_ii + 737);

    auto tk_yzzzzz_xxyyyy = pbuffer.data(idx_kin_ii + 738);

    auto tk_yzzzzz_xxyyyz = pbuffer.data(idx_kin_ii + 739);

    auto tk_yzzzzz_xxyyzz = pbuffer.data(idx_kin_ii + 740);

    auto tk_yzzzzz_xxyzzz = pbuffer.data(idx_kin_ii + 741);

    auto tk_yzzzzz_xxzzzz = pbuffer.data(idx_kin_ii + 742);

    auto tk_yzzzzz_xyyyyy = pbuffer.data(idx_kin_ii + 743);

    auto tk_yzzzzz_xyyyyz = pbuffer.data(idx_kin_ii + 744);

    auto tk_yzzzzz_xyyyzz = pbuffer.data(idx_kin_ii + 745);

    auto tk_yzzzzz_xyyzzz = pbuffer.data(idx_kin_ii + 746);

    auto tk_yzzzzz_xyzzzz = pbuffer.data(idx_kin_ii + 747);

    auto tk_yzzzzz_xzzzzz = pbuffer.data(idx_kin_ii + 748);

    auto tk_yzzzzz_yyyyyy = pbuffer.data(idx_kin_ii + 749);

    auto tk_yzzzzz_yyyyyz = pbuffer.data(idx_kin_ii + 750);

    auto tk_yzzzzz_yyyyzz = pbuffer.data(idx_kin_ii + 751);

    auto tk_yzzzzz_yyyzzz = pbuffer.data(idx_kin_ii + 752);

    auto tk_yzzzzz_yyzzzz = pbuffer.data(idx_kin_ii + 753);

    auto tk_yzzzzz_yzzzzz = pbuffer.data(idx_kin_ii + 754);

    auto tk_yzzzzz_zzzzzz = pbuffer.data(idx_kin_ii + 755);

#pragma omp simd aligned(pa_y,                 \
                             tk_yzzzzz_xxxxxx, \
                             tk_yzzzzz_xxxxxy, \
                             tk_yzzzzz_xxxxxz, \
                             tk_yzzzzz_xxxxyy, \
                             tk_yzzzzz_xxxxyz, \
                             tk_yzzzzz_xxxxzz, \
                             tk_yzzzzz_xxxyyy, \
                             tk_yzzzzz_xxxyyz, \
                             tk_yzzzzz_xxxyzz, \
                             tk_yzzzzz_xxxzzz, \
                             tk_yzzzzz_xxyyyy, \
                             tk_yzzzzz_xxyyyz, \
                             tk_yzzzzz_xxyyzz, \
                             tk_yzzzzz_xxyzzz, \
                             tk_yzzzzz_xxzzzz, \
                             tk_yzzzzz_xyyyyy, \
                             tk_yzzzzz_xyyyyz, \
                             tk_yzzzzz_xyyyzz, \
                             tk_yzzzzz_xyyzzz, \
                             tk_yzzzzz_xyzzzz, \
                             tk_yzzzzz_xzzzzz, \
                             tk_yzzzzz_yyyyyy, \
                             tk_yzzzzz_yyyyyz, \
                             tk_yzzzzz_yyyyzz, \
                             tk_yzzzzz_yyyzzz, \
                             tk_yzzzzz_yyzzzz, \
                             tk_yzzzzz_yzzzzz, \
                             tk_yzzzzz_zzzzzz, \
                             tk_zzzzz_xxxxx,   \
                             tk_zzzzz_xxxxxx,  \
                             tk_zzzzz_xxxxxy,  \
                             tk_zzzzz_xxxxxz,  \
                             tk_zzzzz_xxxxy,   \
                             tk_zzzzz_xxxxyy,  \
                             tk_zzzzz_xxxxyz,  \
                             tk_zzzzz_xxxxz,   \
                             tk_zzzzz_xxxxzz,  \
                             tk_zzzzz_xxxyy,   \
                             tk_zzzzz_xxxyyy,  \
                             tk_zzzzz_xxxyyz,  \
                             tk_zzzzz_xxxyz,   \
                             tk_zzzzz_xxxyzz,  \
                             tk_zzzzz_xxxzz,   \
                             tk_zzzzz_xxxzzz,  \
                             tk_zzzzz_xxyyy,   \
                             tk_zzzzz_xxyyyy,  \
                             tk_zzzzz_xxyyyz,  \
                             tk_zzzzz_xxyyz,   \
                             tk_zzzzz_xxyyzz,  \
                             tk_zzzzz_xxyzz,   \
                             tk_zzzzz_xxyzzz,  \
                             tk_zzzzz_xxzzz,   \
                             tk_zzzzz_xxzzzz,  \
                             tk_zzzzz_xyyyy,   \
                             tk_zzzzz_xyyyyy,  \
                             tk_zzzzz_xyyyyz,  \
                             tk_zzzzz_xyyyz,   \
                             tk_zzzzz_xyyyzz,  \
                             tk_zzzzz_xyyzz,   \
                             tk_zzzzz_xyyzzz,  \
                             tk_zzzzz_xyzzz,   \
                             tk_zzzzz_xyzzzz,  \
                             tk_zzzzz_xzzzz,   \
                             tk_zzzzz_xzzzzz,  \
                             tk_zzzzz_yyyyy,   \
                             tk_zzzzz_yyyyyy,  \
                             tk_zzzzz_yyyyyz,  \
                             tk_zzzzz_yyyyz,   \
                             tk_zzzzz_yyyyzz,  \
                             tk_zzzzz_yyyzz,   \
                             tk_zzzzz_yyyzzz,  \
                             tk_zzzzz_yyzzz,   \
                             tk_zzzzz_yyzzzz,  \
                             tk_zzzzz_yzzzz,   \
                             tk_zzzzz_yzzzzz,  \
                             tk_zzzzz_zzzzz,   \
                             tk_zzzzz_zzzzzz,  \
                             ts_yzzzzz_xxxxxx, \
                             ts_yzzzzz_xxxxxy, \
                             ts_yzzzzz_xxxxxz, \
                             ts_yzzzzz_xxxxyy, \
                             ts_yzzzzz_xxxxyz, \
                             ts_yzzzzz_xxxxzz, \
                             ts_yzzzzz_xxxyyy, \
                             ts_yzzzzz_xxxyyz, \
                             ts_yzzzzz_xxxyzz, \
                             ts_yzzzzz_xxxzzz, \
                             ts_yzzzzz_xxyyyy, \
                             ts_yzzzzz_xxyyyz, \
                             ts_yzzzzz_xxyyzz, \
                             ts_yzzzzz_xxyzzz, \
                             ts_yzzzzz_xxzzzz, \
                             ts_yzzzzz_xyyyyy, \
                             ts_yzzzzz_xyyyyz, \
                             ts_yzzzzz_xyyyzz, \
                             ts_yzzzzz_xyyzzz, \
                             ts_yzzzzz_xyzzzz, \
                             ts_yzzzzz_xzzzzz, \
                             ts_yzzzzz_yyyyyy, \
                             ts_yzzzzz_yyyyyz, \
                             ts_yzzzzz_yyyyzz, \
                             ts_yzzzzz_yyyzzz, \
                             ts_yzzzzz_yyzzzz, \
                             ts_yzzzzz_yzzzzz, \
                             ts_yzzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzzz_xxxxxx[i] = tk_zzzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxxx[i] * fz_0;

        tk_yzzzzz_xxxxxy[i] = tk_zzzzz_xxxxx[i] * fe_0 + tk_zzzzz_xxxxxy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxxy[i] * fz_0;

        tk_yzzzzz_xxxxxz[i] = tk_zzzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxxz[i] * fz_0;

        tk_yzzzzz_xxxxyy[i] = 2.0 * tk_zzzzz_xxxxy[i] * fe_0 + tk_zzzzz_xxxxyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxyy[i] * fz_0;

        tk_yzzzzz_xxxxyz[i] = tk_zzzzz_xxxxz[i] * fe_0 + tk_zzzzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxyz[i] * fz_0;

        tk_yzzzzz_xxxxzz[i] = tk_zzzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxzz[i] * fz_0;

        tk_yzzzzz_xxxyyy[i] = 3.0 * tk_zzzzz_xxxyy[i] * fe_0 + tk_zzzzz_xxxyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxyyy[i] * fz_0;

        tk_yzzzzz_xxxyyz[i] = 2.0 * tk_zzzzz_xxxyz[i] * fe_0 + tk_zzzzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxyyz[i] * fz_0;

        tk_yzzzzz_xxxyzz[i] = tk_zzzzz_xxxzz[i] * fe_0 + tk_zzzzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxyzz[i] * fz_0;

        tk_yzzzzz_xxxzzz[i] = tk_zzzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxzzz[i] * fz_0;

        tk_yzzzzz_xxyyyy[i] = 4.0 * tk_zzzzz_xxyyy[i] * fe_0 + tk_zzzzz_xxyyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyyyy[i] * fz_0;

        tk_yzzzzz_xxyyyz[i] = 3.0 * tk_zzzzz_xxyyz[i] * fe_0 + tk_zzzzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyyyz[i] * fz_0;

        tk_yzzzzz_xxyyzz[i] = 2.0 * tk_zzzzz_xxyzz[i] * fe_0 + tk_zzzzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyyzz[i] * fz_0;

        tk_yzzzzz_xxyzzz[i] = tk_zzzzz_xxzzz[i] * fe_0 + tk_zzzzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyzzz[i] * fz_0;

        tk_yzzzzz_xxzzzz[i] = tk_zzzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxzzzz[i] * fz_0;

        tk_yzzzzz_xyyyyy[i] = 5.0 * tk_zzzzz_xyyyy[i] * fe_0 + tk_zzzzz_xyyyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyyyy[i] * fz_0;

        tk_yzzzzz_xyyyyz[i] = 4.0 * tk_zzzzz_xyyyz[i] * fe_0 + tk_zzzzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyyyz[i] * fz_0;

        tk_yzzzzz_xyyyzz[i] = 3.0 * tk_zzzzz_xyyzz[i] * fe_0 + tk_zzzzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyyzz[i] * fz_0;

        tk_yzzzzz_xyyzzz[i] = 2.0 * tk_zzzzz_xyzzz[i] * fe_0 + tk_zzzzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyzzz[i] * fz_0;

        tk_yzzzzz_xyzzzz[i] = tk_zzzzz_xzzzz[i] * fe_0 + tk_zzzzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyzzzz[i] * fz_0;

        tk_yzzzzz_xzzzzz[i] = tk_zzzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xzzzzz[i] * fz_0;

        tk_yzzzzz_yyyyyy[i] = 6.0 * tk_zzzzz_yyyyy[i] * fe_0 + tk_zzzzz_yyyyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyyyy[i] * fz_0;

        tk_yzzzzz_yyyyyz[i] = 5.0 * tk_zzzzz_yyyyz[i] * fe_0 + tk_zzzzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyyyz[i] * fz_0;

        tk_yzzzzz_yyyyzz[i] = 4.0 * tk_zzzzz_yyyzz[i] * fe_0 + tk_zzzzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyyzz[i] * fz_0;

        tk_yzzzzz_yyyzzz[i] = 3.0 * tk_zzzzz_yyzzz[i] * fe_0 + tk_zzzzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyzzz[i] * fz_0;

        tk_yzzzzz_yyzzzz[i] = 2.0 * tk_zzzzz_yzzzz[i] * fe_0 + tk_zzzzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyzzzz[i] * fz_0;

        tk_yzzzzz_yzzzzz[i] = tk_zzzzz_zzzzz[i] * fe_0 + tk_zzzzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yzzzzz[i] * fz_0;

        tk_yzzzzz_zzzzzz[i] = tk_zzzzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 756-784 components of targeted buffer : II

    auto tk_zzzzzz_xxxxxx = pbuffer.data(idx_kin_ii + 756);

    auto tk_zzzzzz_xxxxxy = pbuffer.data(idx_kin_ii + 757);

    auto tk_zzzzzz_xxxxxz = pbuffer.data(idx_kin_ii + 758);

    auto tk_zzzzzz_xxxxyy = pbuffer.data(idx_kin_ii + 759);

    auto tk_zzzzzz_xxxxyz = pbuffer.data(idx_kin_ii + 760);

    auto tk_zzzzzz_xxxxzz = pbuffer.data(idx_kin_ii + 761);

    auto tk_zzzzzz_xxxyyy = pbuffer.data(idx_kin_ii + 762);

    auto tk_zzzzzz_xxxyyz = pbuffer.data(idx_kin_ii + 763);

    auto tk_zzzzzz_xxxyzz = pbuffer.data(idx_kin_ii + 764);

    auto tk_zzzzzz_xxxzzz = pbuffer.data(idx_kin_ii + 765);

    auto tk_zzzzzz_xxyyyy = pbuffer.data(idx_kin_ii + 766);

    auto tk_zzzzzz_xxyyyz = pbuffer.data(idx_kin_ii + 767);

    auto tk_zzzzzz_xxyyzz = pbuffer.data(idx_kin_ii + 768);

    auto tk_zzzzzz_xxyzzz = pbuffer.data(idx_kin_ii + 769);

    auto tk_zzzzzz_xxzzzz = pbuffer.data(idx_kin_ii + 770);

    auto tk_zzzzzz_xyyyyy = pbuffer.data(idx_kin_ii + 771);

    auto tk_zzzzzz_xyyyyz = pbuffer.data(idx_kin_ii + 772);

    auto tk_zzzzzz_xyyyzz = pbuffer.data(idx_kin_ii + 773);

    auto tk_zzzzzz_xyyzzz = pbuffer.data(idx_kin_ii + 774);

    auto tk_zzzzzz_xyzzzz = pbuffer.data(idx_kin_ii + 775);

    auto tk_zzzzzz_xzzzzz = pbuffer.data(idx_kin_ii + 776);

    auto tk_zzzzzz_yyyyyy = pbuffer.data(idx_kin_ii + 777);

    auto tk_zzzzzz_yyyyyz = pbuffer.data(idx_kin_ii + 778);

    auto tk_zzzzzz_yyyyzz = pbuffer.data(idx_kin_ii + 779);

    auto tk_zzzzzz_yyyzzz = pbuffer.data(idx_kin_ii + 780);

    auto tk_zzzzzz_yyzzzz = pbuffer.data(idx_kin_ii + 781);

    auto tk_zzzzzz_yzzzzz = pbuffer.data(idx_kin_ii + 782);

    auto tk_zzzzzz_zzzzzz = pbuffer.data(idx_kin_ii + 783);

#pragma omp simd aligned(pa_z,                 \
                             tk_zzzz_xxxxxx,   \
                             tk_zzzz_xxxxxy,   \
                             tk_zzzz_xxxxxz,   \
                             tk_zzzz_xxxxyy,   \
                             tk_zzzz_xxxxyz,   \
                             tk_zzzz_xxxxzz,   \
                             tk_zzzz_xxxyyy,   \
                             tk_zzzz_xxxyyz,   \
                             tk_zzzz_xxxyzz,   \
                             tk_zzzz_xxxzzz,   \
                             tk_zzzz_xxyyyy,   \
                             tk_zzzz_xxyyyz,   \
                             tk_zzzz_xxyyzz,   \
                             tk_zzzz_xxyzzz,   \
                             tk_zzzz_xxzzzz,   \
                             tk_zzzz_xyyyyy,   \
                             tk_zzzz_xyyyyz,   \
                             tk_zzzz_xyyyzz,   \
                             tk_zzzz_xyyzzz,   \
                             tk_zzzz_xyzzzz,   \
                             tk_zzzz_xzzzzz,   \
                             tk_zzzz_yyyyyy,   \
                             tk_zzzz_yyyyyz,   \
                             tk_zzzz_yyyyzz,   \
                             tk_zzzz_yyyzzz,   \
                             tk_zzzz_yyzzzz,   \
                             tk_zzzz_yzzzzz,   \
                             tk_zzzz_zzzzzz,   \
                             tk_zzzzz_xxxxx,   \
                             tk_zzzzz_xxxxxx,  \
                             tk_zzzzz_xxxxxy,  \
                             tk_zzzzz_xxxxxz,  \
                             tk_zzzzz_xxxxy,   \
                             tk_zzzzz_xxxxyy,  \
                             tk_zzzzz_xxxxyz,  \
                             tk_zzzzz_xxxxz,   \
                             tk_zzzzz_xxxxzz,  \
                             tk_zzzzz_xxxyy,   \
                             tk_zzzzz_xxxyyy,  \
                             tk_zzzzz_xxxyyz,  \
                             tk_zzzzz_xxxyz,   \
                             tk_zzzzz_xxxyzz,  \
                             tk_zzzzz_xxxzz,   \
                             tk_zzzzz_xxxzzz,  \
                             tk_zzzzz_xxyyy,   \
                             tk_zzzzz_xxyyyy,  \
                             tk_zzzzz_xxyyyz,  \
                             tk_zzzzz_xxyyz,   \
                             tk_zzzzz_xxyyzz,  \
                             tk_zzzzz_xxyzz,   \
                             tk_zzzzz_xxyzzz,  \
                             tk_zzzzz_xxzzz,   \
                             tk_zzzzz_xxzzzz,  \
                             tk_zzzzz_xyyyy,   \
                             tk_zzzzz_xyyyyy,  \
                             tk_zzzzz_xyyyyz,  \
                             tk_zzzzz_xyyyz,   \
                             tk_zzzzz_xyyyzz,  \
                             tk_zzzzz_xyyzz,   \
                             tk_zzzzz_xyyzzz,  \
                             tk_zzzzz_xyzzz,   \
                             tk_zzzzz_xyzzzz,  \
                             tk_zzzzz_xzzzz,   \
                             tk_zzzzz_xzzzzz,  \
                             tk_zzzzz_yyyyy,   \
                             tk_zzzzz_yyyyyy,  \
                             tk_zzzzz_yyyyyz,  \
                             tk_zzzzz_yyyyz,   \
                             tk_zzzzz_yyyyzz,  \
                             tk_zzzzz_yyyzz,   \
                             tk_zzzzz_yyyzzz,  \
                             tk_zzzzz_yyzzz,   \
                             tk_zzzzz_yyzzzz,  \
                             tk_zzzzz_yzzzz,   \
                             tk_zzzzz_yzzzzz,  \
                             tk_zzzzz_zzzzz,   \
                             tk_zzzzz_zzzzzz,  \
                             tk_zzzzzz_xxxxxx, \
                             tk_zzzzzz_xxxxxy, \
                             tk_zzzzzz_xxxxxz, \
                             tk_zzzzzz_xxxxyy, \
                             tk_zzzzzz_xxxxyz, \
                             tk_zzzzzz_xxxxzz, \
                             tk_zzzzzz_xxxyyy, \
                             tk_zzzzzz_xxxyyz, \
                             tk_zzzzzz_xxxyzz, \
                             tk_zzzzzz_xxxzzz, \
                             tk_zzzzzz_xxyyyy, \
                             tk_zzzzzz_xxyyyz, \
                             tk_zzzzzz_xxyyzz, \
                             tk_zzzzzz_xxyzzz, \
                             tk_zzzzzz_xxzzzz, \
                             tk_zzzzzz_xyyyyy, \
                             tk_zzzzzz_xyyyyz, \
                             tk_zzzzzz_xyyyzz, \
                             tk_zzzzzz_xyyzzz, \
                             tk_zzzzzz_xyzzzz, \
                             tk_zzzzzz_xzzzzz, \
                             tk_zzzzzz_yyyyyy, \
                             tk_zzzzzz_yyyyyz, \
                             tk_zzzzzz_yyyyzz, \
                             tk_zzzzzz_yyyzzz, \
                             tk_zzzzzz_yyzzzz, \
                             tk_zzzzzz_yzzzzz, \
                             tk_zzzzzz_zzzzzz, \
                             ts_zzzz_xxxxxx,   \
                             ts_zzzz_xxxxxy,   \
                             ts_zzzz_xxxxxz,   \
                             ts_zzzz_xxxxyy,   \
                             ts_zzzz_xxxxyz,   \
                             ts_zzzz_xxxxzz,   \
                             ts_zzzz_xxxyyy,   \
                             ts_zzzz_xxxyyz,   \
                             ts_zzzz_xxxyzz,   \
                             ts_zzzz_xxxzzz,   \
                             ts_zzzz_xxyyyy,   \
                             ts_zzzz_xxyyyz,   \
                             ts_zzzz_xxyyzz,   \
                             ts_zzzz_xxyzzz,   \
                             ts_zzzz_xxzzzz,   \
                             ts_zzzz_xyyyyy,   \
                             ts_zzzz_xyyyyz,   \
                             ts_zzzz_xyyyzz,   \
                             ts_zzzz_xyyzzz,   \
                             ts_zzzz_xyzzzz,   \
                             ts_zzzz_xzzzzz,   \
                             ts_zzzz_yyyyyy,   \
                             ts_zzzz_yyyyyz,   \
                             ts_zzzz_yyyyzz,   \
                             ts_zzzz_yyyzzz,   \
                             ts_zzzz_yyzzzz,   \
                             ts_zzzz_yzzzzz,   \
                             ts_zzzz_zzzzzz,   \
                             ts_zzzzzz_xxxxxx, \
                             ts_zzzzzz_xxxxxy, \
                             ts_zzzzzz_xxxxxz, \
                             ts_zzzzzz_xxxxyy, \
                             ts_zzzzzz_xxxxyz, \
                             ts_zzzzzz_xxxxzz, \
                             ts_zzzzzz_xxxyyy, \
                             ts_zzzzzz_xxxyyz, \
                             ts_zzzzzz_xxxyzz, \
                             ts_zzzzzz_xxxzzz, \
                             ts_zzzzzz_xxyyyy, \
                             ts_zzzzzz_xxyyyz, \
                             ts_zzzzzz_xxyyzz, \
                             ts_zzzzzz_xxyzzz, \
                             ts_zzzzzz_xxzzzz, \
                             ts_zzzzzz_xyyyyy, \
                             ts_zzzzzz_xyyyyz, \
                             ts_zzzzzz_xyyyzz, \
                             ts_zzzzzz_xyyzzz, \
                             ts_zzzzzz_xyzzzz, \
                             ts_zzzzzz_xzzzzz, \
                             ts_zzzzzz_yyyyyy, \
                             ts_zzzzzz_yyyyyz, \
                             ts_zzzzzz_yyyyzz, \
                             ts_zzzzzz_yyyzzz, \
                             ts_zzzzzz_yyzzzz, \
                             ts_zzzzzz_yzzzzz, \
                             ts_zzzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzzz_xxxxxx[i] = -10.0 * ts_zzzz_xxxxxx[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxxx[i] * fe_0 + tk_zzzzz_xxxxxx[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_xxxxxx[i] * fz_0;

        tk_zzzzzz_xxxxxy[i] = -10.0 * ts_zzzz_xxxxxy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxxy[i] * fe_0 + tk_zzzzz_xxxxxy[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_xxxxxy[i] * fz_0;

        tk_zzzzzz_xxxxxz[i] = -10.0 * ts_zzzz_xxxxxz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxxz[i] * fe_0 + tk_zzzzz_xxxxx[i] * fe_0 +
                              tk_zzzzz_xxxxxz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxxxz[i] * fz_0;

        tk_zzzzzz_xxxxyy[i] = -10.0 * ts_zzzz_xxxxyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxyy[i] * fe_0 + tk_zzzzz_xxxxyy[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_xxxxyy[i] * fz_0;

        tk_zzzzzz_xxxxyz[i] = -10.0 * ts_zzzz_xxxxyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxyz[i] * fe_0 + tk_zzzzz_xxxxy[i] * fe_0 +
                              tk_zzzzz_xxxxyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxxyz[i] * fz_0;

        tk_zzzzzz_xxxxzz[i] = -10.0 * ts_zzzz_xxxxzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxzz[i] * fe_0 + 2.0 * tk_zzzzz_xxxxz[i] * fe_0 +
                              tk_zzzzz_xxxxzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxxzz[i] * fz_0;

        tk_zzzzzz_xxxyyy[i] = -10.0 * ts_zzzz_xxxyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxyyy[i] * fe_0 + tk_zzzzz_xxxyyy[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_xxxyyy[i] * fz_0;

        tk_zzzzzz_xxxyyz[i] = -10.0 * ts_zzzz_xxxyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxyyz[i] * fe_0 + tk_zzzzz_xxxyy[i] * fe_0 +
                              tk_zzzzz_xxxyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxyyz[i] * fz_0;

        tk_zzzzzz_xxxyzz[i] = -10.0 * ts_zzzz_xxxyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxyzz[i] * fe_0 + 2.0 * tk_zzzzz_xxxyz[i] * fe_0 +
                              tk_zzzzz_xxxyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxyzz[i] * fz_0;

        tk_zzzzzz_xxxzzz[i] = -10.0 * ts_zzzz_xxxzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxzzz[i] * fe_0 + 3.0 * tk_zzzzz_xxxzz[i] * fe_0 +
                              tk_zzzzz_xxxzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxzzz[i] * fz_0;

        tk_zzzzzz_xxyyyy[i] = -10.0 * ts_zzzz_xxyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyyyy[i] * fe_0 + tk_zzzzz_xxyyyy[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_xxyyyy[i] * fz_0;

        tk_zzzzzz_xxyyyz[i] = -10.0 * ts_zzzz_xxyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyyyz[i] * fe_0 + tk_zzzzz_xxyyy[i] * fe_0 +
                              tk_zzzzz_xxyyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyyyz[i] * fz_0;

        tk_zzzzzz_xxyyzz[i] = -10.0 * ts_zzzz_xxyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyyzz[i] * fe_0 + 2.0 * tk_zzzzz_xxyyz[i] * fe_0 +
                              tk_zzzzz_xxyyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyyzz[i] * fz_0;

        tk_zzzzzz_xxyzzz[i] = -10.0 * ts_zzzz_xxyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyzzz[i] * fe_0 + 3.0 * tk_zzzzz_xxyzz[i] * fe_0 +
                              tk_zzzzz_xxyzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyzzz[i] * fz_0;

        tk_zzzzzz_xxzzzz[i] = -10.0 * ts_zzzz_xxzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxzzzz[i] * fe_0 + 4.0 * tk_zzzzz_xxzzz[i] * fe_0 +
                              tk_zzzzz_xxzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxzzzz[i] * fz_0;

        tk_zzzzzz_xyyyyy[i] = -10.0 * ts_zzzz_xyyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyyyy[i] * fe_0 + tk_zzzzz_xyyyyy[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_xyyyyy[i] * fz_0;

        tk_zzzzzz_xyyyyz[i] = -10.0 * ts_zzzz_xyyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyyyz[i] * fe_0 + tk_zzzzz_xyyyy[i] * fe_0 +
                              tk_zzzzz_xyyyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyyyz[i] * fz_0;

        tk_zzzzzz_xyyyzz[i] = -10.0 * ts_zzzz_xyyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyyzz[i] * fe_0 + 2.0 * tk_zzzzz_xyyyz[i] * fe_0 +
                              tk_zzzzz_xyyyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyyzz[i] * fz_0;

        tk_zzzzzz_xyyzzz[i] = -10.0 * ts_zzzz_xyyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyzzz[i] * fe_0 + 3.0 * tk_zzzzz_xyyzz[i] * fe_0 +
                              tk_zzzzz_xyyzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyzzz[i] * fz_0;

        tk_zzzzzz_xyzzzz[i] = -10.0 * ts_zzzz_xyzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyzzzz[i] * fe_0 + 4.0 * tk_zzzzz_xyzzz[i] * fe_0 +
                              tk_zzzzz_xyzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyzzzz[i] * fz_0;

        tk_zzzzzz_xzzzzz[i] = -10.0 * ts_zzzz_xzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xzzzzz[i] * fe_0 + 5.0 * tk_zzzzz_xzzzz[i] * fe_0 +
                              tk_zzzzz_xzzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xzzzzz[i] * fz_0;

        tk_zzzzzz_yyyyyy[i] = -10.0 * ts_zzzz_yyyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyyyy[i] * fe_0 + tk_zzzzz_yyyyyy[i] * pa_z[i] +
                              2.0 * ts_zzzzzz_yyyyyy[i] * fz_0;

        tk_zzzzzz_yyyyyz[i] = -10.0 * ts_zzzz_yyyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyyyz[i] * fe_0 + tk_zzzzz_yyyyy[i] * fe_0 +
                              tk_zzzzz_yyyyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyyyz[i] * fz_0;

        tk_zzzzzz_yyyyzz[i] = -10.0 * ts_zzzz_yyyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyyzz[i] * fe_0 + 2.0 * tk_zzzzz_yyyyz[i] * fe_0 +
                              tk_zzzzz_yyyyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyyzz[i] * fz_0;

        tk_zzzzzz_yyyzzz[i] = -10.0 * ts_zzzz_yyyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyzzz[i] * fe_0 + 3.0 * tk_zzzzz_yyyzz[i] * fe_0 +
                              tk_zzzzz_yyyzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyzzz[i] * fz_0;

        tk_zzzzzz_yyzzzz[i] = -10.0 * ts_zzzz_yyzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyzzzz[i] * fe_0 + 4.0 * tk_zzzzz_yyzzz[i] * fe_0 +
                              tk_zzzzz_yyzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyzzzz[i] * fz_0;

        tk_zzzzzz_yzzzzz[i] = -10.0 * ts_zzzz_yzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yzzzzz[i] * fe_0 + 5.0 * tk_zzzzz_yzzzz[i] * fe_0 +
                              tk_zzzzz_yzzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yzzzzz[i] * fz_0;

        tk_zzzzzz_zzzzzz[i] = -10.0 * ts_zzzz_zzzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_zzzzzz[i] * fe_0 + 6.0 * tk_zzzzz_zzzzz[i] * fe_0 +
                              tk_zzzzz_zzzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_zzzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
