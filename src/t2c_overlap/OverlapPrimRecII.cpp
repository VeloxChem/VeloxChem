#include "OverlapPrimRecII.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_ii(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_ii,
                     const size_t idx_ovl_gi,
                     const size_t idx_ovl_hh,
                     const size_t idx_ovl_hi,
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

    auto ts_xxxy_yyyyyy = pbuffer.data(idx_ovl_gi + 49);

    auto ts_xxxy_yyyyyz = pbuffer.data(idx_ovl_gi + 50);

    auto ts_xxxy_yyyyzz = pbuffer.data(idx_ovl_gi + 51);

    auto ts_xxxy_yyyzzz = pbuffer.data(idx_ovl_gi + 52);

    auto ts_xxxy_yyzzzz = pbuffer.data(idx_ovl_gi + 53);

    auto ts_xxxy_yzzzzz = pbuffer.data(idx_ovl_gi + 54);

    auto ts_xxxz_xxxxxx = pbuffer.data(idx_ovl_gi + 56);

    auto ts_xxxz_xxxxxy = pbuffer.data(idx_ovl_gi + 57);

    auto ts_xxxz_xxxxxz = pbuffer.data(idx_ovl_gi + 58);

    auto ts_xxxz_xxxxyy = pbuffer.data(idx_ovl_gi + 59);

    auto ts_xxxz_xxxxzz = pbuffer.data(idx_ovl_gi + 61);

    auto ts_xxxz_xxxyyy = pbuffer.data(idx_ovl_gi + 62);

    auto ts_xxxz_xxxzzz = pbuffer.data(idx_ovl_gi + 65);

    auto ts_xxxz_xxyyyy = pbuffer.data(idx_ovl_gi + 66);

    auto ts_xxxz_xxzzzz = pbuffer.data(idx_ovl_gi + 70);

    auto ts_xxxz_xyyyyy = pbuffer.data(idx_ovl_gi + 71);

    auto ts_xxxz_xzzzzz = pbuffer.data(idx_ovl_gi + 76);

    auto ts_xxxz_yyyyyz = pbuffer.data(idx_ovl_gi + 78);

    auto ts_xxxz_yyyyzz = pbuffer.data(idx_ovl_gi + 79);

    auto ts_xxxz_yyyzzz = pbuffer.data(idx_ovl_gi + 80);

    auto ts_xxxz_yyzzzz = pbuffer.data(idx_ovl_gi + 81);

    auto ts_xxxz_yzzzzz = pbuffer.data(idx_ovl_gi + 82);

    auto ts_xxxz_zzzzzz = pbuffer.data(idx_ovl_gi + 83);

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

    auto ts_xxyz_xxxxxz = pbuffer.data(idx_ovl_gi + 114);

    auto ts_xxyz_xxxxzz = pbuffer.data(idx_ovl_gi + 117);

    auto ts_xxyz_xxxzzz = pbuffer.data(idx_ovl_gi + 121);

    auto ts_xxyz_xxzzzz = pbuffer.data(idx_ovl_gi + 126);

    auto ts_xxyz_xzzzzz = pbuffer.data(idx_ovl_gi + 132);

    auto ts_xxyz_yyyyyz = pbuffer.data(idx_ovl_gi + 134);

    auto ts_xxyz_yyyyzz = pbuffer.data(idx_ovl_gi + 135);

    auto ts_xxyz_yyyzzz = pbuffer.data(idx_ovl_gi + 136);

    auto ts_xxyz_yyzzzz = pbuffer.data(idx_ovl_gi + 137);

    auto ts_xxyz_yzzzzz = pbuffer.data(idx_ovl_gi + 138);

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

    auto ts_xyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 218);

    auto ts_xyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 219);

    auto ts_xyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 220);

    auto ts_xyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 221);

    auto ts_xyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 222);

    auto ts_xyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 223);

    auto ts_xyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 245);

    auto ts_xyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 246);

    auto ts_xyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 247);

    auto ts_xyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 248);

    auto ts_xyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 249);

    auto ts_xyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 250);

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

    auto ts_yyyz_xxxxxz = pbuffer.data(idx_ovl_gi + 310);

    auto ts_yyyz_xxxxyy = pbuffer.data(idx_ovl_gi + 311);

    auto ts_yyyz_xxxxzz = pbuffer.data(idx_ovl_gi + 313);

    auto ts_yyyz_xxxyyy = pbuffer.data(idx_ovl_gi + 314);

    auto ts_yyyz_xxxzzz = pbuffer.data(idx_ovl_gi + 317);

    auto ts_yyyz_xxyyyy = pbuffer.data(idx_ovl_gi + 318);

    auto ts_yyyz_xxzzzz = pbuffer.data(idx_ovl_gi + 322);

    auto ts_yyyz_xyyyyy = pbuffer.data(idx_ovl_gi + 323);

    auto ts_yyyz_xzzzzz = pbuffer.data(idx_ovl_gi + 328);

    auto ts_yyyz_yyyyyy = pbuffer.data(idx_ovl_gi + 329);

    auto ts_yyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 330);

    auto ts_yyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 331);

    auto ts_yyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 332);

    auto ts_yyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 333);

    auto ts_yyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 334);

    auto ts_yyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 335);

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

    auto ts_yzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 385);

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

    // Set up components of auxiliary buffer : HH

    auto ts_xxxxx_xxxxx = pbuffer.data(idx_ovl_hh);

    auto ts_xxxxx_xxxxy = pbuffer.data(idx_ovl_hh + 1);

    auto ts_xxxxx_xxxxz = pbuffer.data(idx_ovl_hh + 2);

    auto ts_xxxxx_xxxyy = pbuffer.data(idx_ovl_hh + 3);

    auto ts_xxxxx_xxxyz = pbuffer.data(idx_ovl_hh + 4);

    auto ts_xxxxx_xxxzz = pbuffer.data(idx_ovl_hh + 5);

    auto ts_xxxxx_xxyyy = pbuffer.data(idx_ovl_hh + 6);

    auto ts_xxxxx_xxyyz = pbuffer.data(idx_ovl_hh + 7);

    auto ts_xxxxx_xxyzz = pbuffer.data(idx_ovl_hh + 8);

    auto ts_xxxxx_xxzzz = pbuffer.data(idx_ovl_hh + 9);

    auto ts_xxxxx_xyyyy = pbuffer.data(idx_ovl_hh + 10);

    auto ts_xxxxx_xyyyz = pbuffer.data(idx_ovl_hh + 11);

    auto ts_xxxxx_xyyzz = pbuffer.data(idx_ovl_hh + 12);

    auto ts_xxxxx_xyzzz = pbuffer.data(idx_ovl_hh + 13);

    auto ts_xxxxx_xzzzz = pbuffer.data(idx_ovl_hh + 14);

    auto ts_xxxxx_yyyyy = pbuffer.data(idx_ovl_hh + 15);

    auto ts_xxxxx_yyyyz = pbuffer.data(idx_ovl_hh + 16);

    auto ts_xxxxx_yyyzz = pbuffer.data(idx_ovl_hh + 17);

    auto ts_xxxxx_yyzzz = pbuffer.data(idx_ovl_hh + 18);

    auto ts_xxxxx_yzzzz = pbuffer.data(idx_ovl_hh + 19);

    auto ts_xxxxx_zzzzz = pbuffer.data(idx_ovl_hh + 20);

    auto ts_xxxxz_xxxxz = pbuffer.data(idx_ovl_hh + 44);

    auto ts_xxxxz_xxxyz = pbuffer.data(idx_ovl_hh + 46);

    auto ts_xxxxz_xxxzz = pbuffer.data(idx_ovl_hh + 47);

    auto ts_xxxxz_xxyyz = pbuffer.data(idx_ovl_hh + 49);

    auto ts_xxxxz_xxyzz = pbuffer.data(idx_ovl_hh + 50);

    auto ts_xxxxz_xxzzz = pbuffer.data(idx_ovl_hh + 51);

    auto ts_xxxxz_xyyyz = pbuffer.data(idx_ovl_hh + 53);

    auto ts_xxxxz_xyyzz = pbuffer.data(idx_ovl_hh + 54);

    auto ts_xxxxz_xyzzz = pbuffer.data(idx_ovl_hh + 55);

    auto ts_xxxxz_xzzzz = pbuffer.data(idx_ovl_hh + 56);

    auto ts_xxxyy_xxxxy = pbuffer.data(idx_ovl_hh + 64);

    auto ts_xxxyy_xxxyy = pbuffer.data(idx_ovl_hh + 66);

    auto ts_xxxyy_xxxyz = pbuffer.data(idx_ovl_hh + 67);

    auto ts_xxxyy_xxyyy = pbuffer.data(idx_ovl_hh + 69);

    auto ts_xxxyy_xxyyz = pbuffer.data(idx_ovl_hh + 70);

    auto ts_xxxyy_xxyzz = pbuffer.data(idx_ovl_hh + 71);

    auto ts_xxxyy_xyyyy = pbuffer.data(idx_ovl_hh + 73);

    auto ts_xxxyy_xyyyz = pbuffer.data(idx_ovl_hh + 74);

    auto ts_xxxyy_xyyzz = pbuffer.data(idx_ovl_hh + 75);

    auto ts_xxxyy_xyzzz = pbuffer.data(idx_ovl_hh + 76);

    auto ts_xxxyy_yyyyy = pbuffer.data(idx_ovl_hh + 78);

    auto ts_xxxyy_yyyyz = pbuffer.data(idx_ovl_hh + 79);

    auto ts_xxxyy_yyyzz = pbuffer.data(idx_ovl_hh + 80);

    auto ts_xxxyy_yyzzz = pbuffer.data(idx_ovl_hh + 81);

    auto ts_xxxyy_yzzzz = pbuffer.data(idx_ovl_hh + 82);

    auto ts_xxxzz_xxxxx = pbuffer.data(idx_ovl_hh + 105);

    auto ts_xxxzz_xxxxy = pbuffer.data(idx_ovl_hh + 106);

    auto ts_xxxzz_xxxxz = pbuffer.data(idx_ovl_hh + 107);

    auto ts_xxxzz_xxxyy = pbuffer.data(idx_ovl_hh + 108);

    auto ts_xxxzz_xxxyz = pbuffer.data(idx_ovl_hh + 109);

    auto ts_xxxzz_xxxzz = pbuffer.data(idx_ovl_hh + 110);

    auto ts_xxxzz_xxyyy = pbuffer.data(idx_ovl_hh + 111);

    auto ts_xxxzz_xxyyz = pbuffer.data(idx_ovl_hh + 112);

    auto ts_xxxzz_xxyzz = pbuffer.data(idx_ovl_hh + 113);

    auto ts_xxxzz_xxzzz = pbuffer.data(idx_ovl_hh + 114);

    auto ts_xxxzz_xyyyy = pbuffer.data(idx_ovl_hh + 115);

    auto ts_xxxzz_xyyyz = pbuffer.data(idx_ovl_hh + 116);

    auto ts_xxxzz_xyyzz = pbuffer.data(idx_ovl_hh + 117);

    auto ts_xxxzz_xyzzz = pbuffer.data(idx_ovl_hh + 118);

    auto ts_xxxzz_xzzzz = pbuffer.data(idx_ovl_hh + 119);

    auto ts_xxxzz_yyyyz = pbuffer.data(idx_ovl_hh + 121);

    auto ts_xxxzz_yyyzz = pbuffer.data(idx_ovl_hh + 122);

    auto ts_xxxzz_yyzzz = pbuffer.data(idx_ovl_hh + 123);

    auto ts_xxxzz_yzzzz = pbuffer.data(idx_ovl_hh + 124);

    auto ts_xxxzz_zzzzz = pbuffer.data(idx_ovl_hh + 125);

    auto ts_xxyyy_xxxxy = pbuffer.data(idx_ovl_hh + 127);

    auto ts_xxyyy_xxxyy = pbuffer.data(idx_ovl_hh + 129);

    auto ts_xxyyy_xxxyz = pbuffer.data(idx_ovl_hh + 130);

    auto ts_xxyyy_xxyyy = pbuffer.data(idx_ovl_hh + 132);

    auto ts_xxyyy_xxyyz = pbuffer.data(idx_ovl_hh + 133);

    auto ts_xxyyy_xxyzz = pbuffer.data(idx_ovl_hh + 134);

    auto ts_xxyyy_xyyyy = pbuffer.data(idx_ovl_hh + 136);

    auto ts_xxyyy_xyyyz = pbuffer.data(idx_ovl_hh + 137);

    auto ts_xxyyy_xyyzz = pbuffer.data(idx_ovl_hh + 138);

    auto ts_xxyyy_xyzzz = pbuffer.data(idx_ovl_hh + 139);

    auto ts_xxyyy_yyyyy = pbuffer.data(idx_ovl_hh + 141);

    auto ts_xxyyy_yyyyz = pbuffer.data(idx_ovl_hh + 142);

    auto ts_xxyyy_yyyzz = pbuffer.data(idx_ovl_hh + 143);

    auto ts_xxyyy_yyzzz = pbuffer.data(idx_ovl_hh + 144);

    auto ts_xxyyy_yzzzz = pbuffer.data(idx_ovl_hh + 145);

    auto ts_xxzzz_xxxxx = pbuffer.data(idx_ovl_hh + 189);

    auto ts_xxzzz_xxxxy = pbuffer.data(idx_ovl_hh + 190);

    auto ts_xxzzz_xxxxz = pbuffer.data(idx_ovl_hh + 191);

    auto ts_xxzzz_xxxyy = pbuffer.data(idx_ovl_hh + 192);

    auto ts_xxzzz_xxxyz = pbuffer.data(idx_ovl_hh + 193);

    auto ts_xxzzz_xxxzz = pbuffer.data(idx_ovl_hh + 194);

    auto ts_xxzzz_xxyyy = pbuffer.data(idx_ovl_hh + 195);

    auto ts_xxzzz_xxyyz = pbuffer.data(idx_ovl_hh + 196);

    auto ts_xxzzz_xxyzz = pbuffer.data(idx_ovl_hh + 197);

    auto ts_xxzzz_xxzzz = pbuffer.data(idx_ovl_hh + 198);

    auto ts_xxzzz_xyyyy = pbuffer.data(idx_ovl_hh + 199);

    auto ts_xxzzz_xyyyz = pbuffer.data(idx_ovl_hh + 200);

    auto ts_xxzzz_xyyzz = pbuffer.data(idx_ovl_hh + 201);

    auto ts_xxzzz_xyzzz = pbuffer.data(idx_ovl_hh + 202);

    auto ts_xxzzz_xzzzz = pbuffer.data(idx_ovl_hh + 203);

    auto ts_xxzzz_yyyyz = pbuffer.data(idx_ovl_hh + 205);

    auto ts_xxzzz_yyyzz = pbuffer.data(idx_ovl_hh + 206);

    auto ts_xxzzz_yyzzz = pbuffer.data(idx_ovl_hh + 207);

    auto ts_xxzzz_yzzzz = pbuffer.data(idx_ovl_hh + 208);

    auto ts_xxzzz_zzzzz = pbuffer.data(idx_ovl_hh + 209);

    auto ts_xyyyy_xxxxy = pbuffer.data(idx_ovl_hh + 211);

    auto ts_xyyyy_xxxyy = pbuffer.data(idx_ovl_hh + 213);

    auto ts_xyyyy_xxxyz = pbuffer.data(idx_ovl_hh + 214);

    auto ts_xyyyy_xxyyy = pbuffer.data(idx_ovl_hh + 216);

    auto ts_xyyyy_xxyyz = pbuffer.data(idx_ovl_hh + 217);

    auto ts_xyyyy_xxyzz = pbuffer.data(idx_ovl_hh + 218);

    auto ts_xyyyy_xyyyy = pbuffer.data(idx_ovl_hh + 220);

    auto ts_xyyyy_xyyyz = pbuffer.data(idx_ovl_hh + 221);

    auto ts_xyyyy_xyyzz = pbuffer.data(idx_ovl_hh + 222);

    auto ts_xyyyy_xyzzz = pbuffer.data(idx_ovl_hh + 223);

    auto ts_xyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 225);

    auto ts_xyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 226);

    auto ts_xyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 227);

    auto ts_xyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 228);

    auto ts_xyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 229);

    auto ts_xyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 256);

    auto ts_xyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 259);

    auto ts_xyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 260);

    auto ts_xyyzz_xyyyz = pbuffer.data(idx_ovl_hh + 263);

    auto ts_xyyzz_xyyzz = pbuffer.data(idx_ovl_hh + 264);

    auto ts_xyyzz_xyzzz = pbuffer.data(idx_ovl_hh + 265);

    auto ts_xyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 268);

    auto ts_xyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 269);

    auto ts_xyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 270);

    auto ts_xyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 271);

    auto ts_xzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 296);

    auto ts_xzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 298);

    auto ts_xzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 299);

    auto ts_xzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 301);

    auto ts_xzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 302);

    auto ts_xzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 303);

    auto ts_xzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 305);

    auto ts_xzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 306);

    auto ts_xzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 307);

    auto ts_xzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 308);

    auto ts_xzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 310);

    auto ts_xzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 311);

    auto ts_xzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 312);

    auto ts_xzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 313);

    auto ts_xzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 314);

    auto ts_yyyyy_xxxxx = pbuffer.data(idx_ovl_hh + 315);

    auto ts_yyyyy_xxxxy = pbuffer.data(idx_ovl_hh + 316);

    auto ts_yyyyy_xxxxz = pbuffer.data(idx_ovl_hh + 317);

    auto ts_yyyyy_xxxyy = pbuffer.data(idx_ovl_hh + 318);

    auto ts_yyyyy_xxxyz = pbuffer.data(idx_ovl_hh + 319);

    auto ts_yyyyy_xxxzz = pbuffer.data(idx_ovl_hh + 320);

    auto ts_yyyyy_xxyyy = pbuffer.data(idx_ovl_hh + 321);

    auto ts_yyyyy_xxyyz = pbuffer.data(idx_ovl_hh + 322);

    auto ts_yyyyy_xxyzz = pbuffer.data(idx_ovl_hh + 323);

    auto ts_yyyyy_xxzzz = pbuffer.data(idx_ovl_hh + 324);

    auto ts_yyyyy_xyyyy = pbuffer.data(idx_ovl_hh + 325);

    auto ts_yyyyy_xyyyz = pbuffer.data(idx_ovl_hh + 326);

    auto ts_yyyyy_xyyzz = pbuffer.data(idx_ovl_hh + 327);

    auto ts_yyyyy_xyzzz = pbuffer.data(idx_ovl_hh + 328);

    auto ts_yyyyy_xzzzz = pbuffer.data(idx_ovl_hh + 329);

    auto ts_yyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 330);

    auto ts_yyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 331);

    auto ts_yyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 332);

    auto ts_yyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 333);

    auto ts_yyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 334);

    auto ts_yyyyy_zzzzz = pbuffer.data(idx_ovl_hh + 335);

    auto ts_yyyyz_xxxxz = pbuffer.data(idx_ovl_hh + 338);

    auto ts_yyyyz_xxxyz = pbuffer.data(idx_ovl_hh + 340);

    auto ts_yyyyz_xxxzz = pbuffer.data(idx_ovl_hh + 341);

    auto ts_yyyyz_xxyyz = pbuffer.data(idx_ovl_hh + 343);

    auto ts_yyyyz_xxyzz = pbuffer.data(idx_ovl_hh + 344);

    auto ts_yyyyz_xxzzz = pbuffer.data(idx_ovl_hh + 345);

    auto ts_yyyyz_xyyyz = pbuffer.data(idx_ovl_hh + 347);

    auto ts_yyyyz_xyyzz = pbuffer.data(idx_ovl_hh + 348);

    auto ts_yyyyz_xyzzz = pbuffer.data(idx_ovl_hh + 349);

    auto ts_yyyyz_xzzzz = pbuffer.data(idx_ovl_hh + 350);

    auto ts_yyyyz_yyyyz = pbuffer.data(idx_ovl_hh + 352);

    auto ts_yyyyz_yyyzz = pbuffer.data(idx_ovl_hh + 353);

    auto ts_yyyyz_yyzzz = pbuffer.data(idx_ovl_hh + 354);

    auto ts_yyyyz_yzzzz = pbuffer.data(idx_ovl_hh + 355);

    auto ts_yyyyz_zzzzz = pbuffer.data(idx_ovl_hh + 356);

    auto ts_yyyzz_xxxxx = pbuffer.data(idx_ovl_hh + 357);

    auto ts_yyyzz_xxxxy = pbuffer.data(idx_ovl_hh + 358);

    auto ts_yyyzz_xxxxz = pbuffer.data(idx_ovl_hh + 359);

    auto ts_yyyzz_xxxyy = pbuffer.data(idx_ovl_hh + 360);

    auto ts_yyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 361);

    auto ts_yyyzz_xxxzz = pbuffer.data(idx_ovl_hh + 362);

    auto ts_yyyzz_xxyyy = pbuffer.data(idx_ovl_hh + 363);

    auto ts_yyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 364);

    auto ts_yyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 365);

    auto ts_yyyzz_xxzzz = pbuffer.data(idx_ovl_hh + 366);

    auto ts_yyyzz_xyyyy = pbuffer.data(idx_ovl_hh + 367);

    auto ts_yyyzz_xyyyz = pbuffer.data(idx_ovl_hh + 368);

    auto ts_yyyzz_xyyzz = pbuffer.data(idx_ovl_hh + 369);

    auto ts_yyyzz_xyzzz = pbuffer.data(idx_ovl_hh + 370);

    auto ts_yyyzz_xzzzz = pbuffer.data(idx_ovl_hh + 371);

    auto ts_yyyzz_yyyyy = pbuffer.data(idx_ovl_hh + 372);

    auto ts_yyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 373);

    auto ts_yyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 374);

    auto ts_yyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 375);

    auto ts_yyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 376);

    auto ts_yyyzz_zzzzz = pbuffer.data(idx_ovl_hh + 377);

    auto ts_yyzzz_xxxxx = pbuffer.data(idx_ovl_hh + 378);

    auto ts_yyzzz_xxxxy = pbuffer.data(idx_ovl_hh + 379);

    auto ts_yyzzz_xxxxz = pbuffer.data(idx_ovl_hh + 380);

    auto ts_yyzzz_xxxyy = pbuffer.data(idx_ovl_hh + 381);

    auto ts_yyzzz_xxxyz = pbuffer.data(idx_ovl_hh + 382);

    auto ts_yyzzz_xxxzz = pbuffer.data(idx_ovl_hh + 383);

    auto ts_yyzzz_xxyyy = pbuffer.data(idx_ovl_hh + 384);

    auto ts_yyzzz_xxyyz = pbuffer.data(idx_ovl_hh + 385);

    auto ts_yyzzz_xxyzz = pbuffer.data(idx_ovl_hh + 386);

    auto ts_yyzzz_xxzzz = pbuffer.data(idx_ovl_hh + 387);

    auto ts_yyzzz_xyyyy = pbuffer.data(idx_ovl_hh + 388);

    auto ts_yyzzz_xyyyz = pbuffer.data(idx_ovl_hh + 389);

    auto ts_yyzzz_xyyzz = pbuffer.data(idx_ovl_hh + 390);

    auto ts_yyzzz_xyzzz = pbuffer.data(idx_ovl_hh + 391);

    auto ts_yyzzz_xzzzz = pbuffer.data(idx_ovl_hh + 392);

    auto ts_yyzzz_yyyyy = pbuffer.data(idx_ovl_hh + 393);

    auto ts_yyzzz_yyyyz = pbuffer.data(idx_ovl_hh + 394);

    auto ts_yyzzz_yyyzz = pbuffer.data(idx_ovl_hh + 395);

    auto ts_yyzzz_yyzzz = pbuffer.data(idx_ovl_hh + 396);

    auto ts_yyzzz_yzzzz = pbuffer.data(idx_ovl_hh + 397);

    auto ts_yyzzz_zzzzz = pbuffer.data(idx_ovl_hh + 398);

    auto ts_yzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 400);

    auto ts_yzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 401);

    auto ts_yzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 402);

    auto ts_yzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 403);

    auto ts_yzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 404);

    auto ts_yzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 405);

    auto ts_yzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 406);

    auto ts_yzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 407);

    auto ts_yzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 408);

    auto ts_yzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 409);

    auto ts_yzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 410);

    auto ts_yzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 411);

    auto ts_yzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 412);

    auto ts_yzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 413);

    auto ts_yzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 414);

    auto ts_yzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 415);

    auto ts_yzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 416);

    auto ts_yzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 417);

    auto ts_yzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 418);

    auto ts_yzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 419);

    auto ts_zzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 420);

    auto ts_zzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 421);

    auto ts_zzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 422);

    auto ts_zzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 423);

    auto ts_zzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 424);

    auto ts_zzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 425);

    auto ts_zzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 426);

    auto ts_zzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 427);

    auto ts_zzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 428);

    auto ts_zzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 429);

    auto ts_zzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 430);

    auto ts_zzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 431);

    auto ts_zzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 432);

    auto ts_zzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 433);

    auto ts_zzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 434);

    auto ts_zzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 435);

    auto ts_zzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 436);

    auto ts_zzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 437);

    auto ts_zzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 438);

    auto ts_zzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 439);

    auto ts_zzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 440);

    // Set up components of auxiliary buffer : HI

    auto ts_xxxxx_xxxxxx = pbuffer.data(idx_ovl_hi);

    auto ts_xxxxx_xxxxxy = pbuffer.data(idx_ovl_hi + 1);

    auto ts_xxxxx_xxxxxz = pbuffer.data(idx_ovl_hi + 2);

    auto ts_xxxxx_xxxxyy = pbuffer.data(idx_ovl_hi + 3);

    auto ts_xxxxx_xxxxyz = pbuffer.data(idx_ovl_hi + 4);

    auto ts_xxxxx_xxxxzz = pbuffer.data(idx_ovl_hi + 5);

    auto ts_xxxxx_xxxyyy = pbuffer.data(idx_ovl_hi + 6);

    auto ts_xxxxx_xxxyyz = pbuffer.data(idx_ovl_hi + 7);

    auto ts_xxxxx_xxxyzz = pbuffer.data(idx_ovl_hi + 8);

    auto ts_xxxxx_xxxzzz = pbuffer.data(idx_ovl_hi + 9);

    auto ts_xxxxx_xxyyyy = pbuffer.data(idx_ovl_hi + 10);

    auto ts_xxxxx_xxyyyz = pbuffer.data(idx_ovl_hi + 11);

    auto ts_xxxxx_xxyyzz = pbuffer.data(idx_ovl_hi + 12);

    auto ts_xxxxx_xxyzzz = pbuffer.data(idx_ovl_hi + 13);

    auto ts_xxxxx_xxzzzz = pbuffer.data(idx_ovl_hi + 14);

    auto ts_xxxxx_xyyyyy = pbuffer.data(idx_ovl_hi + 15);

    auto ts_xxxxx_xyyyyz = pbuffer.data(idx_ovl_hi + 16);

    auto ts_xxxxx_xyyyzz = pbuffer.data(idx_ovl_hi + 17);

    auto ts_xxxxx_xyyzzz = pbuffer.data(idx_ovl_hi + 18);

    auto ts_xxxxx_xyzzzz = pbuffer.data(idx_ovl_hi + 19);

    auto ts_xxxxx_xzzzzz = pbuffer.data(idx_ovl_hi + 20);

    auto ts_xxxxx_yyyyyy = pbuffer.data(idx_ovl_hi + 21);

    auto ts_xxxxx_yyyyyz = pbuffer.data(idx_ovl_hi + 22);

    auto ts_xxxxx_yyyyzz = pbuffer.data(idx_ovl_hi + 23);

    auto ts_xxxxx_yyyzzz = pbuffer.data(idx_ovl_hi + 24);

    auto ts_xxxxx_yyzzzz = pbuffer.data(idx_ovl_hi + 25);

    auto ts_xxxxx_yzzzzz = pbuffer.data(idx_ovl_hi + 26);

    auto ts_xxxxx_zzzzzz = pbuffer.data(idx_ovl_hi + 27);

    auto ts_xxxxy_xxxxxx = pbuffer.data(idx_ovl_hi + 28);

    auto ts_xxxxy_xxxxxy = pbuffer.data(idx_ovl_hi + 29);

    auto ts_xxxxy_xxxxxz = pbuffer.data(idx_ovl_hi + 30);

    auto ts_xxxxy_xxxxyy = pbuffer.data(idx_ovl_hi + 31);

    auto ts_xxxxy_xxxxzz = pbuffer.data(idx_ovl_hi + 33);

    auto ts_xxxxy_xxxyyy = pbuffer.data(idx_ovl_hi + 34);

    auto ts_xxxxy_xxxzzz = pbuffer.data(idx_ovl_hi + 37);

    auto ts_xxxxy_xxyyyy = pbuffer.data(idx_ovl_hi + 38);

    auto ts_xxxxy_xxzzzz = pbuffer.data(idx_ovl_hi + 42);

    auto ts_xxxxy_xyyyyy = pbuffer.data(idx_ovl_hi + 43);

    auto ts_xxxxy_xzzzzz = pbuffer.data(idx_ovl_hi + 48);

    auto ts_xxxxy_yyyyyy = pbuffer.data(idx_ovl_hi + 49);

    auto ts_xxxxy_yyyyyz = pbuffer.data(idx_ovl_hi + 50);

    auto ts_xxxxy_yyyyzz = pbuffer.data(idx_ovl_hi + 51);

    auto ts_xxxxy_yyyzzz = pbuffer.data(idx_ovl_hi + 52);

    auto ts_xxxxy_yyzzzz = pbuffer.data(idx_ovl_hi + 53);

    auto ts_xxxxy_yzzzzz = pbuffer.data(idx_ovl_hi + 54);

    auto ts_xxxxz_xxxxxx = pbuffer.data(idx_ovl_hi + 56);

    auto ts_xxxxz_xxxxxy = pbuffer.data(idx_ovl_hi + 57);

    auto ts_xxxxz_xxxxxz = pbuffer.data(idx_ovl_hi + 58);

    auto ts_xxxxz_xxxxyy = pbuffer.data(idx_ovl_hi + 59);

    auto ts_xxxxz_xxxxyz = pbuffer.data(idx_ovl_hi + 60);

    auto ts_xxxxz_xxxxzz = pbuffer.data(idx_ovl_hi + 61);

    auto ts_xxxxz_xxxyyy = pbuffer.data(idx_ovl_hi + 62);

    auto ts_xxxxz_xxxyyz = pbuffer.data(idx_ovl_hi + 63);

    auto ts_xxxxz_xxxyzz = pbuffer.data(idx_ovl_hi + 64);

    auto ts_xxxxz_xxxzzz = pbuffer.data(idx_ovl_hi + 65);

    auto ts_xxxxz_xxyyyy = pbuffer.data(idx_ovl_hi + 66);

    auto ts_xxxxz_xxyyyz = pbuffer.data(idx_ovl_hi + 67);

    auto ts_xxxxz_xxyyzz = pbuffer.data(idx_ovl_hi + 68);

    auto ts_xxxxz_xxyzzz = pbuffer.data(idx_ovl_hi + 69);

    auto ts_xxxxz_xxzzzz = pbuffer.data(idx_ovl_hi + 70);

    auto ts_xxxxz_xyyyyy = pbuffer.data(idx_ovl_hi + 71);

    auto ts_xxxxz_xyyyyz = pbuffer.data(idx_ovl_hi + 72);

    auto ts_xxxxz_xyyyzz = pbuffer.data(idx_ovl_hi + 73);

    auto ts_xxxxz_xyyzzz = pbuffer.data(idx_ovl_hi + 74);

    auto ts_xxxxz_xyzzzz = pbuffer.data(idx_ovl_hi + 75);

    auto ts_xxxxz_xzzzzz = pbuffer.data(idx_ovl_hi + 76);

    auto ts_xxxxz_yyyyyz = pbuffer.data(idx_ovl_hi + 78);

    auto ts_xxxxz_yyyyzz = pbuffer.data(idx_ovl_hi + 79);

    auto ts_xxxxz_yyyzzz = pbuffer.data(idx_ovl_hi + 80);

    auto ts_xxxxz_yyzzzz = pbuffer.data(idx_ovl_hi + 81);

    auto ts_xxxxz_yzzzzz = pbuffer.data(idx_ovl_hi + 82);

    auto ts_xxxxz_zzzzzz = pbuffer.data(idx_ovl_hi + 83);

    auto ts_xxxyy_xxxxxx = pbuffer.data(idx_ovl_hi + 84);

    auto ts_xxxyy_xxxxxy = pbuffer.data(idx_ovl_hi + 85);

    auto ts_xxxyy_xxxxxz = pbuffer.data(idx_ovl_hi + 86);

    auto ts_xxxyy_xxxxyy = pbuffer.data(idx_ovl_hi + 87);

    auto ts_xxxyy_xxxxyz = pbuffer.data(idx_ovl_hi + 88);

    auto ts_xxxyy_xxxxzz = pbuffer.data(idx_ovl_hi + 89);

    auto ts_xxxyy_xxxyyy = pbuffer.data(idx_ovl_hi + 90);

    auto ts_xxxyy_xxxyyz = pbuffer.data(idx_ovl_hi + 91);

    auto ts_xxxyy_xxxyzz = pbuffer.data(idx_ovl_hi + 92);

    auto ts_xxxyy_xxxzzz = pbuffer.data(idx_ovl_hi + 93);

    auto ts_xxxyy_xxyyyy = pbuffer.data(idx_ovl_hi + 94);

    auto ts_xxxyy_xxyyyz = pbuffer.data(idx_ovl_hi + 95);

    auto ts_xxxyy_xxyyzz = pbuffer.data(idx_ovl_hi + 96);

    auto ts_xxxyy_xxyzzz = pbuffer.data(idx_ovl_hi + 97);

    auto ts_xxxyy_xxzzzz = pbuffer.data(idx_ovl_hi + 98);

    auto ts_xxxyy_xyyyyy = pbuffer.data(idx_ovl_hi + 99);

    auto ts_xxxyy_xyyyyz = pbuffer.data(idx_ovl_hi + 100);

    auto ts_xxxyy_xyyyzz = pbuffer.data(idx_ovl_hi + 101);

    auto ts_xxxyy_xyyzzz = pbuffer.data(idx_ovl_hi + 102);

    auto ts_xxxyy_xyzzzz = pbuffer.data(idx_ovl_hi + 103);

    auto ts_xxxyy_xzzzzz = pbuffer.data(idx_ovl_hi + 104);

    auto ts_xxxyy_yyyyyy = pbuffer.data(idx_ovl_hi + 105);

    auto ts_xxxyy_yyyyyz = pbuffer.data(idx_ovl_hi + 106);

    auto ts_xxxyy_yyyyzz = pbuffer.data(idx_ovl_hi + 107);

    auto ts_xxxyy_yyyzzz = pbuffer.data(idx_ovl_hi + 108);

    auto ts_xxxyy_yyzzzz = pbuffer.data(idx_ovl_hi + 109);

    auto ts_xxxyy_yzzzzz = pbuffer.data(idx_ovl_hi + 110);

    auto ts_xxxyy_zzzzzz = pbuffer.data(idx_ovl_hi + 111);

    auto ts_xxxyz_xxxxxz = pbuffer.data(idx_ovl_hi + 114);

    auto ts_xxxyz_xxxxzz = pbuffer.data(idx_ovl_hi + 117);

    auto ts_xxxyz_xxxzzz = pbuffer.data(idx_ovl_hi + 121);

    auto ts_xxxyz_xxzzzz = pbuffer.data(idx_ovl_hi + 126);

    auto ts_xxxyz_xzzzzz = pbuffer.data(idx_ovl_hi + 132);

    auto ts_xxxyz_yyyyyz = pbuffer.data(idx_ovl_hi + 134);

    auto ts_xxxyz_yyyyzz = pbuffer.data(idx_ovl_hi + 135);

    auto ts_xxxyz_yyyzzz = pbuffer.data(idx_ovl_hi + 136);

    auto ts_xxxyz_yyzzzz = pbuffer.data(idx_ovl_hi + 137);

    auto ts_xxxyz_yzzzzz = pbuffer.data(idx_ovl_hi + 138);

    auto ts_xxxzz_xxxxxx = pbuffer.data(idx_ovl_hi + 140);

    auto ts_xxxzz_xxxxxy = pbuffer.data(idx_ovl_hi + 141);

    auto ts_xxxzz_xxxxxz = pbuffer.data(idx_ovl_hi + 142);

    auto ts_xxxzz_xxxxyy = pbuffer.data(idx_ovl_hi + 143);

    auto ts_xxxzz_xxxxyz = pbuffer.data(idx_ovl_hi + 144);

    auto ts_xxxzz_xxxxzz = pbuffer.data(idx_ovl_hi + 145);

    auto ts_xxxzz_xxxyyy = pbuffer.data(idx_ovl_hi + 146);

    auto ts_xxxzz_xxxyyz = pbuffer.data(idx_ovl_hi + 147);

    auto ts_xxxzz_xxxyzz = pbuffer.data(idx_ovl_hi + 148);

    auto ts_xxxzz_xxxzzz = pbuffer.data(idx_ovl_hi + 149);

    auto ts_xxxzz_xxyyyy = pbuffer.data(idx_ovl_hi + 150);

    auto ts_xxxzz_xxyyyz = pbuffer.data(idx_ovl_hi + 151);

    auto ts_xxxzz_xxyyzz = pbuffer.data(idx_ovl_hi + 152);

    auto ts_xxxzz_xxyzzz = pbuffer.data(idx_ovl_hi + 153);

    auto ts_xxxzz_xxzzzz = pbuffer.data(idx_ovl_hi + 154);

    auto ts_xxxzz_xyyyyy = pbuffer.data(idx_ovl_hi + 155);

    auto ts_xxxzz_xyyyyz = pbuffer.data(idx_ovl_hi + 156);

    auto ts_xxxzz_xyyyzz = pbuffer.data(idx_ovl_hi + 157);

    auto ts_xxxzz_xyyzzz = pbuffer.data(idx_ovl_hi + 158);

    auto ts_xxxzz_xyzzzz = pbuffer.data(idx_ovl_hi + 159);

    auto ts_xxxzz_xzzzzz = pbuffer.data(idx_ovl_hi + 160);

    auto ts_xxxzz_yyyyyy = pbuffer.data(idx_ovl_hi + 161);

    auto ts_xxxzz_yyyyyz = pbuffer.data(idx_ovl_hi + 162);

    auto ts_xxxzz_yyyyzz = pbuffer.data(idx_ovl_hi + 163);

    auto ts_xxxzz_yyyzzz = pbuffer.data(idx_ovl_hi + 164);

    auto ts_xxxzz_yyzzzz = pbuffer.data(idx_ovl_hi + 165);

    auto ts_xxxzz_yzzzzz = pbuffer.data(idx_ovl_hi + 166);

    auto ts_xxxzz_zzzzzz = pbuffer.data(idx_ovl_hi + 167);

    auto ts_xxyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 168);

    auto ts_xxyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 169);

    auto ts_xxyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 170);

    auto ts_xxyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 171);

    auto ts_xxyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 172);

    auto ts_xxyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 173);

    auto ts_xxyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 174);

    auto ts_xxyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 175);

    auto ts_xxyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 176);

    auto ts_xxyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 177);

    auto ts_xxyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 178);

    auto ts_xxyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 179);

    auto ts_xxyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 180);

    auto ts_xxyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 181);

    auto ts_xxyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 182);

    auto ts_xxyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 183);

    auto ts_xxyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 184);

    auto ts_xxyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 185);

    auto ts_xxyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 186);

    auto ts_xxyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 187);

    auto ts_xxyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 188);

    auto ts_xxyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 189);

    auto ts_xxyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 190);

    auto ts_xxyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 191);

    auto ts_xxyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 192);

    auto ts_xxyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 193);

    auto ts_xxyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 194);

    auto ts_xxyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 195);

    auto ts_xxyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 197);

    auto ts_xxyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 198);

    auto ts_xxyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 199);

    auto ts_xxyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 201);

    auto ts_xxyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 202);

    auto ts_xxyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 205);

    auto ts_xxyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 206);

    auto ts_xxyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 210);

    auto ts_xxyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 211);

    auto ts_xxyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 216);

    auto ts_xxyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 218);

    auto ts_xxyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 219);

    auto ts_xxyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 220);

    auto ts_xxyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 221);

    auto ts_xxyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 222);

    auto ts_xxyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 223);

    auto ts_xxyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 224);

    auto ts_xxyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 226);

    auto ts_xxyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 229);

    auto ts_xxyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 233);

    auto ts_xxyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 238);

    auto ts_xxyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 244);

    auto ts_xxyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 245);

    auto ts_xxyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 246);

    auto ts_xxyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 247);

    auto ts_xxyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 248);

    auto ts_xxyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 249);

    auto ts_xxyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 250);

    auto ts_xxzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 252);

    auto ts_xxzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 253);

    auto ts_xxzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 254);

    auto ts_xxzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 255);

    auto ts_xxzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 256);

    auto ts_xxzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 257);

    auto ts_xxzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 258);

    auto ts_xxzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 259);

    auto ts_xxzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 260);

    auto ts_xxzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 261);

    auto ts_xxzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 262);

    auto ts_xxzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 263);

    auto ts_xxzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 264);

    auto ts_xxzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 265);

    auto ts_xxzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 266);

    auto ts_xxzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 267);

    auto ts_xxzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 268);

    auto ts_xxzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 269);

    auto ts_xxzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 270);

    auto ts_xxzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 271);

    auto ts_xxzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 272);

    auto ts_xxzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 273);

    auto ts_xxzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 274);

    auto ts_xxzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 275);

    auto ts_xxzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 276);

    auto ts_xxzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 277);

    auto ts_xxzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 278);

    auto ts_xxzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 279);

    auto ts_xyyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 280);

    auto ts_xyyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 281);

    auto ts_xyyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 283);

    auto ts_xyyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 284);

    auto ts_xyyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 286);

    auto ts_xyyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 287);

    auto ts_xyyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 288);

    auto ts_xyyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 290);

    auto ts_xyyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 291);

    auto ts_xyyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 292);

    auto ts_xyyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 293);

    auto ts_xyyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 295);

    auto ts_xyyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 296);

    auto ts_xyyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 297);

    auto ts_xyyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 298);

    auto ts_xyyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 299);

    auto ts_xyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 301);

    auto ts_xyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 302);

    auto ts_xyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 303);

    auto ts_xyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 304);

    auto ts_xyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 305);

    auto ts_xyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 306);

    auto ts_xyyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 307);

    auto ts_xyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 330);

    auto ts_xyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 331);

    auto ts_xyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 332);

    auto ts_xyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 333);

    auto ts_xyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 334);

    auto ts_xyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 335);

    auto ts_xyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 340);

    auto ts_xyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 343);

    auto ts_xyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 344);

    auto ts_xyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 347);

    auto ts_xyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 348);

    auto ts_xyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 349);

    auto ts_xyyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 352);

    auto ts_xyyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 353);

    auto ts_xyyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 354);

    auto ts_xyyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 355);

    auto ts_xyyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 357);

    auto ts_xyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 358);

    auto ts_xyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 359);

    auto ts_xyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 360);

    auto ts_xyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 361);

    auto ts_xyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 362);

    auto ts_xyyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 363);

    auto ts_xyzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 385);

    auto ts_xyzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 386);

    auto ts_xyzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 387);

    auto ts_xyzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 388);

    auto ts_xyzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 389);

    auto ts_xyzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 390);

    auto ts_xzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 392);

    auto ts_xzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 394);

    auto ts_xzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 396);

    auto ts_xzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 397);

    auto ts_xzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 399);

    auto ts_xzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 400);

    auto ts_xzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 401);

    auto ts_xzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 403);

    auto ts_xzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 404);

    auto ts_xzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 405);

    auto ts_xzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 406);

    auto ts_xzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 408);

    auto ts_xzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 409);

    auto ts_xzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 410);

    auto ts_xzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 411);

    auto ts_xzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 412);

    auto ts_xzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 413);

    auto ts_xzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 414);

    auto ts_xzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 415);

    auto ts_xzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 416);

    auto ts_xzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 417);

    auto ts_xzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 418);

    auto ts_xzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 419);

    auto ts_yyyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 420);

    auto ts_yyyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 421);

    auto ts_yyyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 422);

    auto ts_yyyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 423);

    auto ts_yyyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 424);

    auto ts_yyyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 425);

    auto ts_yyyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 426);

    auto ts_yyyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 427);

    auto ts_yyyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 428);

    auto ts_yyyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 429);

    auto ts_yyyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 430);

    auto ts_yyyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 431);

    auto ts_yyyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 432);

    auto ts_yyyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 433);

    auto ts_yyyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 434);

    auto ts_yyyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 435);

    auto ts_yyyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 436);

    auto ts_yyyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 437);

    auto ts_yyyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 438);

    auto ts_yyyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 439);

    auto ts_yyyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 440);

    auto ts_yyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 441);

    auto ts_yyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 442);

    auto ts_yyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 443);

    auto ts_yyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 444);

    auto ts_yyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 445);

    auto ts_yyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 446);

    auto ts_yyyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 447);

    auto ts_yyyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 449);

    auto ts_yyyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 450);

    auto ts_yyyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 451);

    auto ts_yyyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 452);

    auto ts_yyyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 453);

    auto ts_yyyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 454);

    auto ts_yyyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 455);

    auto ts_yyyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 456);

    auto ts_yyyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 457);

    auto ts_yyyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 458);

    auto ts_yyyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 459);

    auto ts_yyyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 460);

    auto ts_yyyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 461);

    auto ts_yyyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 462);

    auto ts_yyyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 463);

    auto ts_yyyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 464);

    auto ts_yyyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 465);

    auto ts_yyyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 466);

    auto ts_yyyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 467);

    auto ts_yyyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 468);

    auto ts_yyyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 469);

    auto ts_yyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 470);

    auto ts_yyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 471);

    auto ts_yyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 472);

    auto ts_yyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 473);

    auto ts_yyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 474);

    auto ts_yyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 475);

    auto ts_yyyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 476);

    auto ts_yyyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 477);

    auto ts_yyyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 478);

    auto ts_yyyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 479);

    auto ts_yyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 480);

    auto ts_yyyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 481);

    auto ts_yyyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 482);

    auto ts_yyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 483);

    auto ts_yyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 484);

    auto ts_yyyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 485);

    auto ts_yyyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 486);

    auto ts_yyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 487);

    auto ts_yyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 488);

    auto ts_yyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 489);

    auto ts_yyyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 490);

    auto ts_yyyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 491);

    auto ts_yyyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 492);

    auto ts_yyyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 493);

    auto ts_yyyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 494);

    auto ts_yyyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 495);

    auto ts_yyyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 496);

    auto ts_yyyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 497);

    auto ts_yyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 498);

    auto ts_yyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 499);

    auto ts_yyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 500);

    auto ts_yyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 501);

    auto ts_yyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 502);

    auto ts_yyyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 503);

    auto ts_yyzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 504);

    auto ts_yyzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 505);

    auto ts_yyzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 506);

    auto ts_yyzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 507);

    auto ts_yyzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 508);

    auto ts_yyzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 509);

    auto ts_yyzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 510);

    auto ts_yyzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 511);

    auto ts_yyzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 512);

    auto ts_yyzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 513);

    auto ts_yyzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 514);

    auto ts_yyzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 515);

    auto ts_yyzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 516);

    auto ts_yyzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 517);

    auto ts_yyzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 518);

    auto ts_yyzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 519);

    auto ts_yyzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 520);

    auto ts_yyzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 521);

    auto ts_yyzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 522);

    auto ts_yyzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 523);

    auto ts_yyzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 524);

    auto ts_yyzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 525);

    auto ts_yyzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 526);

    auto ts_yyzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 527);

    auto ts_yyzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 528);

    auto ts_yyzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 529);

    auto ts_yyzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 530);

    auto ts_yyzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 531);

    auto ts_yzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 532);

    auto ts_yzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 533);

    auto ts_yzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 534);

    auto ts_yzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 535);

    auto ts_yzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 536);

    auto ts_yzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 537);

    auto ts_yzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 538);

    auto ts_yzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 539);

    auto ts_yzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 540);

    auto ts_yzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 541);

    auto ts_yzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 542);

    auto ts_yzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 543);

    auto ts_yzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 544);

    auto ts_yzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 545);

    auto ts_yzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 546);

    auto ts_yzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 547);

    auto ts_yzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 548);

    auto ts_yzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 549);

    auto ts_yzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 550);

    auto ts_yzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 551);

    auto ts_yzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 552);

    auto ts_yzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 553);

    auto ts_yzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 554);

    auto ts_yzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 555);

    auto ts_yzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 556);

    auto ts_yzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 557);

    auto ts_yzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 558);

    auto ts_yzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 559);

    auto ts_zzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 560);

    auto ts_zzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 561);

    auto ts_zzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 562);

    auto ts_zzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 563);

    auto ts_zzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 564);

    auto ts_zzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 565);

    auto ts_zzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 566);

    auto ts_zzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 567);

    auto ts_zzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 568);

    auto ts_zzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 569);

    auto ts_zzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 570);

    auto ts_zzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 571);

    auto ts_zzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 572);

    auto ts_zzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 573);

    auto ts_zzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 574);

    auto ts_zzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 575);

    auto ts_zzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 576);

    auto ts_zzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 577);

    auto ts_zzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 578);

    auto ts_zzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 579);

    auto ts_zzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 580);

    auto ts_zzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 581);

    auto ts_zzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 582);

    auto ts_zzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 583);

    auto ts_zzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 584);

    auto ts_zzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 585);

    auto ts_zzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 586);

    auto ts_zzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 587);

    // Set up 0-28 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, ts_xxxx_xxxxxx, ts_xxxx_xxxxxy, ts_xxxx_xxxxxz, ts_xxxx_xxxxyy, ts_xxxx_xxxxyz, ts_xxxx_xxxxzz, ts_xxxx_xxxyyy, ts_xxxx_xxxyyz, ts_xxxx_xxxyzz, ts_xxxx_xxxzzz, ts_xxxx_xxyyyy, ts_xxxx_xxyyyz, ts_xxxx_xxyyzz, ts_xxxx_xxyzzz, ts_xxxx_xxzzzz, ts_xxxx_xyyyyy, ts_xxxx_xyyyyz, ts_xxxx_xyyyzz, ts_xxxx_xyyzzz, ts_xxxx_xyzzzz, ts_xxxx_xzzzzz, ts_xxxx_yyyyyy, ts_xxxx_yyyyyz, ts_xxxx_yyyyzz, ts_xxxx_yyyzzz, ts_xxxx_yyzzzz, ts_xxxx_yzzzzz, ts_xxxx_zzzzzz, ts_xxxxx_xxxxx, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxy, ts_xxxxx_xxxxxz, ts_xxxxx_xxxxy, ts_xxxxx_xxxxyy, ts_xxxxx_xxxxyz, ts_xxxxx_xxxxz, ts_xxxxx_xxxxzz, ts_xxxxx_xxxyy, ts_xxxxx_xxxyyy, ts_xxxxx_xxxyyz, ts_xxxxx_xxxyz, ts_xxxxx_xxxyzz, ts_xxxxx_xxxzz, ts_xxxxx_xxxzzz, ts_xxxxx_xxyyy, ts_xxxxx_xxyyyy, ts_xxxxx_xxyyyz, ts_xxxxx_xxyyz, ts_xxxxx_xxyyzz, ts_xxxxx_xxyzz, ts_xxxxx_xxyzzz, ts_xxxxx_xxzzz, ts_xxxxx_xxzzzz, ts_xxxxx_xyyyy, ts_xxxxx_xyyyyy, ts_xxxxx_xyyyyz, ts_xxxxx_xyyyz, ts_xxxxx_xyyyzz, ts_xxxxx_xyyzz, ts_xxxxx_xyyzzz, ts_xxxxx_xyzzz, ts_xxxxx_xyzzzz, ts_xxxxx_xzzzz, ts_xxxxx_xzzzzz, ts_xxxxx_yyyyy, ts_xxxxx_yyyyyy, ts_xxxxx_yyyyyz, ts_xxxxx_yyyyz, ts_xxxxx_yyyyzz, ts_xxxxx_yyyzz, ts_xxxxx_yyyzzz, ts_xxxxx_yyzzz, ts_xxxxx_yyzzzz, ts_xxxxx_yzzzz, ts_xxxxx_yzzzzz, ts_xxxxx_zzzzz, ts_xxxxx_zzzzzz, ts_xxxxxx_xxxxxx, ts_xxxxxx_xxxxxy, ts_xxxxxx_xxxxxz, ts_xxxxxx_xxxxyy, ts_xxxxxx_xxxxyz, ts_xxxxxx_xxxxzz, ts_xxxxxx_xxxyyy, ts_xxxxxx_xxxyyz, ts_xxxxxx_xxxyzz, ts_xxxxxx_xxxzzz, ts_xxxxxx_xxyyyy, ts_xxxxxx_xxyyyz, ts_xxxxxx_xxyyzz, ts_xxxxxx_xxyzzz, ts_xxxxxx_xxzzzz, ts_xxxxxx_xyyyyy, ts_xxxxxx_xyyyyz, ts_xxxxxx_xyyyzz, ts_xxxxxx_xyyzzz, ts_xxxxxx_xyzzzz, ts_xxxxxx_xzzzzz, ts_xxxxxx_yyyyyy, ts_xxxxxx_yyyyyz, ts_xxxxxx_yyyyzz, ts_xxxxxx_yyyzzz, ts_xxxxxx_yyzzzz, ts_xxxxxx_yzzzzz, ts_xxxxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_xxxxxx[i] = 5.0 * ts_xxxx_xxxxxx[i] * fe_0 + 6.0 * ts_xxxxx_xxxxx[i] * fe_0 + ts_xxxxx_xxxxxx[i] * pa_x[i];

        ts_xxxxxx_xxxxxy[i] = 5.0 * ts_xxxx_xxxxxy[i] * fe_0 + 5.0 * ts_xxxxx_xxxxy[i] * fe_0 + ts_xxxxx_xxxxxy[i] * pa_x[i];

        ts_xxxxxx_xxxxxz[i] = 5.0 * ts_xxxx_xxxxxz[i] * fe_0 + 5.0 * ts_xxxxx_xxxxz[i] * fe_0 + ts_xxxxx_xxxxxz[i] * pa_x[i];

        ts_xxxxxx_xxxxyy[i] = 5.0 * ts_xxxx_xxxxyy[i] * fe_0 + 4.0 * ts_xxxxx_xxxyy[i] * fe_0 + ts_xxxxx_xxxxyy[i] * pa_x[i];

        ts_xxxxxx_xxxxyz[i] = 5.0 * ts_xxxx_xxxxyz[i] * fe_0 + 4.0 * ts_xxxxx_xxxyz[i] * fe_0 + ts_xxxxx_xxxxyz[i] * pa_x[i];

        ts_xxxxxx_xxxxzz[i] = 5.0 * ts_xxxx_xxxxzz[i] * fe_0 + 4.0 * ts_xxxxx_xxxzz[i] * fe_0 + ts_xxxxx_xxxxzz[i] * pa_x[i];

        ts_xxxxxx_xxxyyy[i] = 5.0 * ts_xxxx_xxxyyy[i] * fe_0 + 3.0 * ts_xxxxx_xxyyy[i] * fe_0 + ts_xxxxx_xxxyyy[i] * pa_x[i];

        ts_xxxxxx_xxxyyz[i] = 5.0 * ts_xxxx_xxxyyz[i] * fe_0 + 3.0 * ts_xxxxx_xxyyz[i] * fe_0 + ts_xxxxx_xxxyyz[i] * pa_x[i];

        ts_xxxxxx_xxxyzz[i] = 5.0 * ts_xxxx_xxxyzz[i] * fe_0 + 3.0 * ts_xxxxx_xxyzz[i] * fe_0 + ts_xxxxx_xxxyzz[i] * pa_x[i];

        ts_xxxxxx_xxxzzz[i] = 5.0 * ts_xxxx_xxxzzz[i] * fe_0 + 3.0 * ts_xxxxx_xxzzz[i] * fe_0 + ts_xxxxx_xxxzzz[i] * pa_x[i];

        ts_xxxxxx_xxyyyy[i] = 5.0 * ts_xxxx_xxyyyy[i] * fe_0 + 2.0 * ts_xxxxx_xyyyy[i] * fe_0 + ts_xxxxx_xxyyyy[i] * pa_x[i];

        ts_xxxxxx_xxyyyz[i] = 5.0 * ts_xxxx_xxyyyz[i] * fe_0 + 2.0 * ts_xxxxx_xyyyz[i] * fe_0 + ts_xxxxx_xxyyyz[i] * pa_x[i];

        ts_xxxxxx_xxyyzz[i] = 5.0 * ts_xxxx_xxyyzz[i] * fe_0 + 2.0 * ts_xxxxx_xyyzz[i] * fe_0 + ts_xxxxx_xxyyzz[i] * pa_x[i];

        ts_xxxxxx_xxyzzz[i] = 5.0 * ts_xxxx_xxyzzz[i] * fe_0 + 2.0 * ts_xxxxx_xyzzz[i] * fe_0 + ts_xxxxx_xxyzzz[i] * pa_x[i];

        ts_xxxxxx_xxzzzz[i] = 5.0 * ts_xxxx_xxzzzz[i] * fe_0 + 2.0 * ts_xxxxx_xzzzz[i] * fe_0 + ts_xxxxx_xxzzzz[i] * pa_x[i];

        ts_xxxxxx_xyyyyy[i] = 5.0 * ts_xxxx_xyyyyy[i] * fe_0 + ts_xxxxx_yyyyy[i] * fe_0 + ts_xxxxx_xyyyyy[i] * pa_x[i];

        ts_xxxxxx_xyyyyz[i] = 5.0 * ts_xxxx_xyyyyz[i] * fe_0 + ts_xxxxx_yyyyz[i] * fe_0 + ts_xxxxx_xyyyyz[i] * pa_x[i];

        ts_xxxxxx_xyyyzz[i] = 5.0 * ts_xxxx_xyyyzz[i] * fe_0 + ts_xxxxx_yyyzz[i] * fe_0 + ts_xxxxx_xyyyzz[i] * pa_x[i];

        ts_xxxxxx_xyyzzz[i] = 5.0 * ts_xxxx_xyyzzz[i] * fe_0 + ts_xxxxx_yyzzz[i] * fe_0 + ts_xxxxx_xyyzzz[i] * pa_x[i];

        ts_xxxxxx_xyzzzz[i] = 5.0 * ts_xxxx_xyzzzz[i] * fe_0 + ts_xxxxx_yzzzz[i] * fe_0 + ts_xxxxx_xyzzzz[i] * pa_x[i];

        ts_xxxxxx_xzzzzz[i] = 5.0 * ts_xxxx_xzzzzz[i] * fe_0 + ts_xxxxx_zzzzz[i] * fe_0 + ts_xxxxx_xzzzzz[i] * pa_x[i];

        ts_xxxxxx_yyyyyy[i] = 5.0 * ts_xxxx_yyyyyy[i] * fe_0 + ts_xxxxx_yyyyyy[i] * pa_x[i];

        ts_xxxxxx_yyyyyz[i] = 5.0 * ts_xxxx_yyyyyz[i] * fe_0 + ts_xxxxx_yyyyyz[i] * pa_x[i];

        ts_xxxxxx_yyyyzz[i] = 5.0 * ts_xxxx_yyyyzz[i] * fe_0 + ts_xxxxx_yyyyzz[i] * pa_x[i];

        ts_xxxxxx_yyyzzz[i] = 5.0 * ts_xxxx_yyyzzz[i] * fe_0 + ts_xxxxx_yyyzzz[i] * pa_x[i];

        ts_xxxxxx_yyzzzz[i] = 5.0 * ts_xxxx_yyzzzz[i] * fe_0 + ts_xxxxx_yyzzzz[i] * pa_x[i];

        ts_xxxxxx_yzzzzz[i] = 5.0 * ts_xxxx_yzzzzz[i] * fe_0 + ts_xxxxx_yzzzzz[i] * pa_x[i];

        ts_xxxxxx_zzzzzz[i] = 5.0 * ts_xxxx_zzzzzz[i] * fe_0 + ts_xxxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxxx_xxxxx, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxy, ts_xxxxx_xxxxxz, ts_xxxxx_xxxxy, ts_xxxxx_xxxxyy, ts_xxxxx_xxxxyz, ts_xxxxx_xxxxz, ts_xxxxx_xxxxzz, ts_xxxxx_xxxyy, ts_xxxxx_xxxyyy, ts_xxxxx_xxxyyz, ts_xxxxx_xxxyz, ts_xxxxx_xxxyzz, ts_xxxxx_xxxzz, ts_xxxxx_xxxzzz, ts_xxxxx_xxyyy, ts_xxxxx_xxyyyy, ts_xxxxx_xxyyyz, ts_xxxxx_xxyyz, ts_xxxxx_xxyyzz, ts_xxxxx_xxyzz, ts_xxxxx_xxyzzz, ts_xxxxx_xxzzz, ts_xxxxx_xxzzzz, ts_xxxxx_xyyyy, ts_xxxxx_xyyyyy, ts_xxxxx_xyyyyz, ts_xxxxx_xyyyz, ts_xxxxx_xyyyzz, ts_xxxxx_xyyzz, ts_xxxxx_xyyzzz, ts_xxxxx_xyzzz, ts_xxxxx_xyzzzz, ts_xxxxx_xzzzz, ts_xxxxx_xzzzzz, ts_xxxxx_zzzzzz, ts_xxxxxy_xxxxxx, ts_xxxxxy_xxxxxy, ts_xxxxxy_xxxxxz, ts_xxxxxy_xxxxyy, ts_xxxxxy_xxxxyz, ts_xxxxxy_xxxxzz, ts_xxxxxy_xxxyyy, ts_xxxxxy_xxxyyz, ts_xxxxxy_xxxyzz, ts_xxxxxy_xxxzzz, ts_xxxxxy_xxyyyy, ts_xxxxxy_xxyyyz, ts_xxxxxy_xxyyzz, ts_xxxxxy_xxyzzz, ts_xxxxxy_xxzzzz, ts_xxxxxy_xyyyyy, ts_xxxxxy_xyyyyz, ts_xxxxxy_xyyyzz, ts_xxxxxy_xyyzzz, ts_xxxxxy_xyzzzz, ts_xxxxxy_xzzzzz, ts_xxxxxy_yyyyyy, ts_xxxxxy_yyyyyz, ts_xxxxxy_yyyyzz, ts_xxxxxy_yyyzzz, ts_xxxxxy_yyzzzz, ts_xxxxxy_yzzzzz, ts_xxxxxy_zzzzzz, ts_xxxxy_yyyyyy, ts_xxxxy_yyyyyz, ts_xxxxy_yyyyzz, ts_xxxxy_yyyzzz, ts_xxxxy_yyzzzz, ts_xxxxy_yzzzzz, ts_xxxy_yyyyyy, ts_xxxy_yyyyyz, ts_xxxy_yyyyzz, ts_xxxy_yyyzzz, ts_xxxy_yyzzzz, ts_xxxy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxy_xxxxxx[i] = ts_xxxxx_xxxxxx[i] * pa_y[i];

        ts_xxxxxy_xxxxxy[i] = ts_xxxxx_xxxxx[i] * fe_0 + ts_xxxxx_xxxxxy[i] * pa_y[i];

        ts_xxxxxy_xxxxxz[i] = ts_xxxxx_xxxxxz[i] * pa_y[i];

        ts_xxxxxy_xxxxyy[i] = 2.0 * ts_xxxxx_xxxxy[i] * fe_0 + ts_xxxxx_xxxxyy[i] * pa_y[i];

        ts_xxxxxy_xxxxyz[i] = ts_xxxxx_xxxxz[i] * fe_0 + ts_xxxxx_xxxxyz[i] * pa_y[i];

        ts_xxxxxy_xxxxzz[i] = ts_xxxxx_xxxxzz[i] * pa_y[i];

        ts_xxxxxy_xxxyyy[i] = 3.0 * ts_xxxxx_xxxyy[i] * fe_0 + ts_xxxxx_xxxyyy[i] * pa_y[i];

        ts_xxxxxy_xxxyyz[i] = 2.0 * ts_xxxxx_xxxyz[i] * fe_0 + ts_xxxxx_xxxyyz[i] * pa_y[i];

        ts_xxxxxy_xxxyzz[i] = ts_xxxxx_xxxzz[i] * fe_0 + ts_xxxxx_xxxyzz[i] * pa_y[i];

        ts_xxxxxy_xxxzzz[i] = ts_xxxxx_xxxzzz[i] * pa_y[i];

        ts_xxxxxy_xxyyyy[i] = 4.0 * ts_xxxxx_xxyyy[i] * fe_0 + ts_xxxxx_xxyyyy[i] * pa_y[i];

        ts_xxxxxy_xxyyyz[i] = 3.0 * ts_xxxxx_xxyyz[i] * fe_0 + ts_xxxxx_xxyyyz[i] * pa_y[i];

        ts_xxxxxy_xxyyzz[i] = 2.0 * ts_xxxxx_xxyzz[i] * fe_0 + ts_xxxxx_xxyyzz[i] * pa_y[i];

        ts_xxxxxy_xxyzzz[i] = ts_xxxxx_xxzzz[i] * fe_0 + ts_xxxxx_xxyzzz[i] * pa_y[i];

        ts_xxxxxy_xxzzzz[i] = ts_xxxxx_xxzzzz[i] * pa_y[i];

        ts_xxxxxy_xyyyyy[i] = 5.0 * ts_xxxxx_xyyyy[i] * fe_0 + ts_xxxxx_xyyyyy[i] * pa_y[i];

        ts_xxxxxy_xyyyyz[i] = 4.0 * ts_xxxxx_xyyyz[i] * fe_0 + ts_xxxxx_xyyyyz[i] * pa_y[i];

        ts_xxxxxy_xyyyzz[i] = 3.0 * ts_xxxxx_xyyzz[i] * fe_0 + ts_xxxxx_xyyyzz[i] * pa_y[i];

        ts_xxxxxy_xyyzzz[i] = 2.0 * ts_xxxxx_xyzzz[i] * fe_0 + ts_xxxxx_xyyzzz[i] * pa_y[i];

        ts_xxxxxy_xyzzzz[i] = ts_xxxxx_xzzzz[i] * fe_0 + ts_xxxxx_xyzzzz[i] * pa_y[i];

        ts_xxxxxy_xzzzzz[i] = ts_xxxxx_xzzzzz[i] * pa_y[i];

        ts_xxxxxy_yyyyyy[i] = 4.0 * ts_xxxy_yyyyyy[i] * fe_0 + ts_xxxxy_yyyyyy[i] * pa_x[i];

        ts_xxxxxy_yyyyyz[i] = 4.0 * ts_xxxy_yyyyyz[i] * fe_0 + ts_xxxxy_yyyyyz[i] * pa_x[i];

        ts_xxxxxy_yyyyzz[i] = 4.0 * ts_xxxy_yyyyzz[i] * fe_0 + ts_xxxxy_yyyyzz[i] * pa_x[i];

        ts_xxxxxy_yyyzzz[i] = 4.0 * ts_xxxy_yyyzzz[i] * fe_0 + ts_xxxxy_yyyzzz[i] * pa_x[i];

        ts_xxxxxy_yyzzzz[i] = 4.0 * ts_xxxy_yyzzzz[i] * fe_0 + ts_xxxxy_yyzzzz[i] * pa_x[i];

        ts_xxxxxy_yzzzzz[i] = 4.0 * ts_xxxy_yzzzzz[i] * fe_0 + ts_xxxxy_yzzzzz[i] * pa_x[i];

        ts_xxxxxy_zzzzzz[i] = ts_xxxxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxxx_xxxxx, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxy, ts_xxxxx_xxxxxz, ts_xxxxx_xxxxy, ts_xxxxx_xxxxyy, ts_xxxxx_xxxxyz, ts_xxxxx_xxxxz, ts_xxxxx_xxxxzz, ts_xxxxx_xxxyy, ts_xxxxx_xxxyyy, ts_xxxxx_xxxyyz, ts_xxxxx_xxxyz, ts_xxxxx_xxxyzz, ts_xxxxx_xxxzz, ts_xxxxx_xxxzzz, ts_xxxxx_xxyyy, ts_xxxxx_xxyyyy, ts_xxxxx_xxyyyz, ts_xxxxx_xxyyz, ts_xxxxx_xxyyzz, ts_xxxxx_xxyzz, ts_xxxxx_xxyzzz, ts_xxxxx_xxzzz, ts_xxxxx_xxzzzz, ts_xxxxx_xyyyy, ts_xxxxx_xyyyyy, ts_xxxxx_xyyyyz, ts_xxxxx_xyyyz, ts_xxxxx_xyyyzz, ts_xxxxx_xyyzz, ts_xxxxx_xyyzzz, ts_xxxxx_xyzzz, ts_xxxxx_xyzzzz, ts_xxxxx_xzzzz, ts_xxxxx_xzzzzz, ts_xxxxx_yyyyyy, ts_xxxxxz_xxxxxx, ts_xxxxxz_xxxxxy, ts_xxxxxz_xxxxxz, ts_xxxxxz_xxxxyy, ts_xxxxxz_xxxxyz, ts_xxxxxz_xxxxzz, ts_xxxxxz_xxxyyy, ts_xxxxxz_xxxyyz, ts_xxxxxz_xxxyzz, ts_xxxxxz_xxxzzz, ts_xxxxxz_xxyyyy, ts_xxxxxz_xxyyyz, ts_xxxxxz_xxyyzz, ts_xxxxxz_xxyzzz, ts_xxxxxz_xxzzzz, ts_xxxxxz_xyyyyy, ts_xxxxxz_xyyyyz, ts_xxxxxz_xyyyzz, ts_xxxxxz_xyyzzz, ts_xxxxxz_xyzzzz, ts_xxxxxz_xzzzzz, ts_xxxxxz_yyyyyy, ts_xxxxxz_yyyyyz, ts_xxxxxz_yyyyzz, ts_xxxxxz_yyyzzz, ts_xxxxxz_yyzzzz, ts_xxxxxz_yzzzzz, ts_xxxxxz_zzzzzz, ts_xxxxz_yyyyyz, ts_xxxxz_yyyyzz, ts_xxxxz_yyyzzz, ts_xxxxz_yyzzzz, ts_xxxxz_yzzzzz, ts_xxxxz_zzzzzz, ts_xxxz_yyyyyz, ts_xxxz_yyyyzz, ts_xxxz_yyyzzz, ts_xxxz_yyzzzz, ts_xxxz_yzzzzz, ts_xxxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxz_xxxxxx[i] = ts_xxxxx_xxxxxx[i] * pa_z[i];

        ts_xxxxxz_xxxxxy[i] = ts_xxxxx_xxxxxy[i] * pa_z[i];

        ts_xxxxxz_xxxxxz[i] = ts_xxxxx_xxxxx[i] * fe_0 + ts_xxxxx_xxxxxz[i] * pa_z[i];

        ts_xxxxxz_xxxxyy[i] = ts_xxxxx_xxxxyy[i] * pa_z[i];

        ts_xxxxxz_xxxxyz[i] = ts_xxxxx_xxxxy[i] * fe_0 + ts_xxxxx_xxxxyz[i] * pa_z[i];

        ts_xxxxxz_xxxxzz[i] = 2.0 * ts_xxxxx_xxxxz[i] * fe_0 + ts_xxxxx_xxxxzz[i] * pa_z[i];

        ts_xxxxxz_xxxyyy[i] = ts_xxxxx_xxxyyy[i] * pa_z[i];

        ts_xxxxxz_xxxyyz[i] = ts_xxxxx_xxxyy[i] * fe_0 + ts_xxxxx_xxxyyz[i] * pa_z[i];

        ts_xxxxxz_xxxyzz[i] = 2.0 * ts_xxxxx_xxxyz[i] * fe_0 + ts_xxxxx_xxxyzz[i] * pa_z[i];

        ts_xxxxxz_xxxzzz[i] = 3.0 * ts_xxxxx_xxxzz[i] * fe_0 + ts_xxxxx_xxxzzz[i] * pa_z[i];

        ts_xxxxxz_xxyyyy[i] = ts_xxxxx_xxyyyy[i] * pa_z[i];

        ts_xxxxxz_xxyyyz[i] = ts_xxxxx_xxyyy[i] * fe_0 + ts_xxxxx_xxyyyz[i] * pa_z[i];

        ts_xxxxxz_xxyyzz[i] = 2.0 * ts_xxxxx_xxyyz[i] * fe_0 + ts_xxxxx_xxyyzz[i] * pa_z[i];

        ts_xxxxxz_xxyzzz[i] = 3.0 * ts_xxxxx_xxyzz[i] * fe_0 + ts_xxxxx_xxyzzz[i] * pa_z[i];

        ts_xxxxxz_xxzzzz[i] = 4.0 * ts_xxxxx_xxzzz[i] * fe_0 + ts_xxxxx_xxzzzz[i] * pa_z[i];

        ts_xxxxxz_xyyyyy[i] = ts_xxxxx_xyyyyy[i] * pa_z[i];

        ts_xxxxxz_xyyyyz[i] = ts_xxxxx_xyyyy[i] * fe_0 + ts_xxxxx_xyyyyz[i] * pa_z[i];

        ts_xxxxxz_xyyyzz[i] = 2.0 * ts_xxxxx_xyyyz[i] * fe_0 + ts_xxxxx_xyyyzz[i] * pa_z[i];

        ts_xxxxxz_xyyzzz[i] = 3.0 * ts_xxxxx_xyyzz[i] * fe_0 + ts_xxxxx_xyyzzz[i] * pa_z[i];

        ts_xxxxxz_xyzzzz[i] = 4.0 * ts_xxxxx_xyzzz[i] * fe_0 + ts_xxxxx_xyzzzz[i] * pa_z[i];

        ts_xxxxxz_xzzzzz[i] = 5.0 * ts_xxxxx_xzzzz[i] * fe_0 + ts_xxxxx_xzzzzz[i] * pa_z[i];

        ts_xxxxxz_yyyyyy[i] = ts_xxxxx_yyyyyy[i] * pa_z[i];

        ts_xxxxxz_yyyyyz[i] = 4.0 * ts_xxxz_yyyyyz[i] * fe_0 + ts_xxxxz_yyyyyz[i] * pa_x[i];

        ts_xxxxxz_yyyyzz[i] = 4.0 * ts_xxxz_yyyyzz[i] * fe_0 + ts_xxxxz_yyyyzz[i] * pa_x[i];

        ts_xxxxxz_yyyzzz[i] = 4.0 * ts_xxxz_yyyzzz[i] * fe_0 + ts_xxxxz_yyyzzz[i] * pa_x[i];

        ts_xxxxxz_yyzzzz[i] = 4.0 * ts_xxxz_yyzzzz[i] * fe_0 + ts_xxxxz_yyzzzz[i] * pa_x[i];

        ts_xxxxxz_yzzzzz[i] = 4.0 * ts_xxxz_yzzzzz[i] * fe_0 + ts_xxxxz_yzzzzz[i] * pa_x[i];

        ts_xxxxxz_zzzzzz[i] = 4.0 * ts_xxxz_zzzzzz[i] * fe_0 + ts_xxxxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xxxxxx, ts_xxxx_xxxxxz, ts_xxxx_xxxxzz, ts_xxxx_xxxzzz, ts_xxxx_xxzzzz, ts_xxxx_xzzzzz, ts_xxxxy_xxxxxx, ts_xxxxy_xxxxxz, ts_xxxxy_xxxxzz, ts_xxxxy_xxxzzz, ts_xxxxy_xxzzzz, ts_xxxxy_xzzzzz, ts_xxxxyy_xxxxxx, ts_xxxxyy_xxxxxy, ts_xxxxyy_xxxxxz, ts_xxxxyy_xxxxyy, ts_xxxxyy_xxxxyz, ts_xxxxyy_xxxxzz, ts_xxxxyy_xxxyyy, ts_xxxxyy_xxxyyz, ts_xxxxyy_xxxyzz, ts_xxxxyy_xxxzzz, ts_xxxxyy_xxyyyy, ts_xxxxyy_xxyyyz, ts_xxxxyy_xxyyzz, ts_xxxxyy_xxyzzz, ts_xxxxyy_xxzzzz, ts_xxxxyy_xyyyyy, ts_xxxxyy_xyyyyz, ts_xxxxyy_xyyyzz, ts_xxxxyy_xyyzzz, ts_xxxxyy_xyzzzz, ts_xxxxyy_xzzzzz, ts_xxxxyy_yyyyyy, ts_xxxxyy_yyyyyz, ts_xxxxyy_yyyyzz, ts_xxxxyy_yyyzzz, ts_xxxxyy_yyzzzz, ts_xxxxyy_yzzzzz, ts_xxxxyy_zzzzzz, ts_xxxyy_xxxxxy, ts_xxxyy_xxxxy, ts_xxxyy_xxxxyy, ts_xxxyy_xxxxyz, ts_xxxyy_xxxyy, ts_xxxyy_xxxyyy, ts_xxxyy_xxxyyz, ts_xxxyy_xxxyz, ts_xxxyy_xxxyzz, ts_xxxyy_xxyyy, ts_xxxyy_xxyyyy, ts_xxxyy_xxyyyz, ts_xxxyy_xxyyz, ts_xxxyy_xxyyzz, ts_xxxyy_xxyzz, ts_xxxyy_xxyzzz, ts_xxxyy_xyyyy, ts_xxxyy_xyyyyy, ts_xxxyy_xyyyyz, ts_xxxyy_xyyyz, ts_xxxyy_xyyyzz, ts_xxxyy_xyyzz, ts_xxxyy_xyyzzz, ts_xxxyy_xyzzz, ts_xxxyy_xyzzzz, ts_xxxyy_yyyyy, ts_xxxyy_yyyyyy, ts_xxxyy_yyyyyz, ts_xxxyy_yyyyz, ts_xxxyy_yyyyzz, ts_xxxyy_yyyzz, ts_xxxyy_yyyzzz, ts_xxxyy_yyzzz, ts_xxxyy_yyzzzz, ts_xxxyy_yzzzz, ts_xxxyy_yzzzzz, ts_xxxyy_zzzzzz, ts_xxyy_xxxxxy, ts_xxyy_xxxxyy, ts_xxyy_xxxxyz, ts_xxyy_xxxyyy, ts_xxyy_xxxyyz, ts_xxyy_xxxyzz, ts_xxyy_xxyyyy, ts_xxyy_xxyyyz, ts_xxyy_xxyyzz, ts_xxyy_xxyzzz, ts_xxyy_xyyyyy, ts_xxyy_xyyyyz, ts_xxyy_xyyyzz, ts_xxyy_xyyzzz, ts_xxyy_xyzzzz, ts_xxyy_yyyyyy, ts_xxyy_yyyyyz, ts_xxyy_yyyyzz, ts_xxyy_yyyzzz, ts_xxyy_yyzzzz, ts_xxyy_yzzzzz, ts_xxyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyy_xxxxxx[i] = ts_xxxx_xxxxxx[i] * fe_0 + ts_xxxxy_xxxxxx[i] * pa_y[i];

        ts_xxxxyy_xxxxxy[i] = 3.0 * ts_xxyy_xxxxxy[i] * fe_0 + 5.0 * ts_xxxyy_xxxxy[i] * fe_0 + ts_xxxyy_xxxxxy[i] * pa_x[i];

        ts_xxxxyy_xxxxxz[i] = ts_xxxx_xxxxxz[i] * fe_0 + ts_xxxxy_xxxxxz[i] * pa_y[i];

        ts_xxxxyy_xxxxyy[i] = 3.0 * ts_xxyy_xxxxyy[i] * fe_0 + 4.0 * ts_xxxyy_xxxyy[i] * fe_0 + ts_xxxyy_xxxxyy[i] * pa_x[i];

        ts_xxxxyy_xxxxyz[i] = 3.0 * ts_xxyy_xxxxyz[i] * fe_0 + 4.0 * ts_xxxyy_xxxyz[i] * fe_0 + ts_xxxyy_xxxxyz[i] * pa_x[i];

        ts_xxxxyy_xxxxzz[i] = ts_xxxx_xxxxzz[i] * fe_0 + ts_xxxxy_xxxxzz[i] * pa_y[i];

        ts_xxxxyy_xxxyyy[i] = 3.0 * ts_xxyy_xxxyyy[i] * fe_0 + 3.0 * ts_xxxyy_xxyyy[i] * fe_0 + ts_xxxyy_xxxyyy[i] * pa_x[i];

        ts_xxxxyy_xxxyyz[i] = 3.0 * ts_xxyy_xxxyyz[i] * fe_0 + 3.0 * ts_xxxyy_xxyyz[i] * fe_0 + ts_xxxyy_xxxyyz[i] * pa_x[i];

        ts_xxxxyy_xxxyzz[i] = 3.0 * ts_xxyy_xxxyzz[i] * fe_0 + 3.0 * ts_xxxyy_xxyzz[i] * fe_0 + ts_xxxyy_xxxyzz[i] * pa_x[i];

        ts_xxxxyy_xxxzzz[i] = ts_xxxx_xxxzzz[i] * fe_0 + ts_xxxxy_xxxzzz[i] * pa_y[i];

        ts_xxxxyy_xxyyyy[i] = 3.0 * ts_xxyy_xxyyyy[i] * fe_0 + 2.0 * ts_xxxyy_xyyyy[i] * fe_0 + ts_xxxyy_xxyyyy[i] * pa_x[i];

        ts_xxxxyy_xxyyyz[i] = 3.0 * ts_xxyy_xxyyyz[i] * fe_0 + 2.0 * ts_xxxyy_xyyyz[i] * fe_0 + ts_xxxyy_xxyyyz[i] * pa_x[i];

        ts_xxxxyy_xxyyzz[i] = 3.0 * ts_xxyy_xxyyzz[i] * fe_0 + 2.0 * ts_xxxyy_xyyzz[i] * fe_0 + ts_xxxyy_xxyyzz[i] * pa_x[i];

        ts_xxxxyy_xxyzzz[i] = 3.0 * ts_xxyy_xxyzzz[i] * fe_0 + 2.0 * ts_xxxyy_xyzzz[i] * fe_0 + ts_xxxyy_xxyzzz[i] * pa_x[i];

        ts_xxxxyy_xxzzzz[i] = ts_xxxx_xxzzzz[i] * fe_0 + ts_xxxxy_xxzzzz[i] * pa_y[i];

        ts_xxxxyy_xyyyyy[i] = 3.0 * ts_xxyy_xyyyyy[i] * fe_0 + ts_xxxyy_yyyyy[i] * fe_0 + ts_xxxyy_xyyyyy[i] * pa_x[i];

        ts_xxxxyy_xyyyyz[i] = 3.0 * ts_xxyy_xyyyyz[i] * fe_0 + ts_xxxyy_yyyyz[i] * fe_0 + ts_xxxyy_xyyyyz[i] * pa_x[i];

        ts_xxxxyy_xyyyzz[i] = 3.0 * ts_xxyy_xyyyzz[i] * fe_0 + ts_xxxyy_yyyzz[i] * fe_0 + ts_xxxyy_xyyyzz[i] * pa_x[i];

        ts_xxxxyy_xyyzzz[i] = 3.0 * ts_xxyy_xyyzzz[i] * fe_0 + ts_xxxyy_yyzzz[i] * fe_0 + ts_xxxyy_xyyzzz[i] * pa_x[i];

        ts_xxxxyy_xyzzzz[i] = 3.0 * ts_xxyy_xyzzzz[i] * fe_0 + ts_xxxyy_yzzzz[i] * fe_0 + ts_xxxyy_xyzzzz[i] * pa_x[i];

        ts_xxxxyy_xzzzzz[i] = ts_xxxx_xzzzzz[i] * fe_0 + ts_xxxxy_xzzzzz[i] * pa_y[i];

        ts_xxxxyy_yyyyyy[i] = 3.0 * ts_xxyy_yyyyyy[i] * fe_0 + ts_xxxyy_yyyyyy[i] * pa_x[i];

        ts_xxxxyy_yyyyyz[i] = 3.0 * ts_xxyy_yyyyyz[i] * fe_0 + ts_xxxyy_yyyyyz[i] * pa_x[i];

        ts_xxxxyy_yyyyzz[i] = 3.0 * ts_xxyy_yyyyzz[i] * fe_0 + ts_xxxyy_yyyyzz[i] * pa_x[i];

        ts_xxxxyy_yyyzzz[i] = 3.0 * ts_xxyy_yyyzzz[i] * fe_0 + ts_xxxyy_yyyzzz[i] * pa_x[i];

        ts_xxxxyy_yyzzzz[i] = 3.0 * ts_xxyy_yyzzzz[i] * fe_0 + ts_xxxyy_yyzzzz[i] * pa_x[i];

        ts_xxxxyy_yzzzzz[i] = 3.0 * ts_xxyy_yzzzzz[i] * fe_0 + ts_xxxyy_yzzzzz[i] * pa_x[i];

        ts_xxxxyy_zzzzzz[i] = 3.0 * ts_xxyy_zzzzzz[i] * fe_0 + ts_xxxyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxxy_xxxxxy, ts_xxxxy_xxxxyy, ts_xxxxy_xxxyyy, ts_xxxxy_xxyyyy, ts_xxxxy_xyyyyy, ts_xxxxy_yyyyyy, ts_xxxxyz_xxxxxx, ts_xxxxyz_xxxxxy, ts_xxxxyz_xxxxxz, ts_xxxxyz_xxxxyy, ts_xxxxyz_xxxxyz, ts_xxxxyz_xxxxzz, ts_xxxxyz_xxxyyy, ts_xxxxyz_xxxyyz, ts_xxxxyz_xxxyzz, ts_xxxxyz_xxxzzz, ts_xxxxyz_xxyyyy, ts_xxxxyz_xxyyyz, ts_xxxxyz_xxyyzz, ts_xxxxyz_xxyzzz, ts_xxxxyz_xxzzzz, ts_xxxxyz_xyyyyy, ts_xxxxyz_xyyyyz, ts_xxxxyz_xyyyzz, ts_xxxxyz_xyyzzz, ts_xxxxyz_xyzzzz, ts_xxxxyz_xzzzzz, ts_xxxxyz_yyyyyy, ts_xxxxyz_yyyyyz, ts_xxxxyz_yyyyzz, ts_xxxxyz_yyyzzz, ts_xxxxyz_yyzzzz, ts_xxxxyz_yzzzzz, ts_xxxxyz_zzzzzz, ts_xxxxz_xxxxxx, ts_xxxxz_xxxxxz, ts_xxxxz_xxxxyz, ts_xxxxz_xxxxz, ts_xxxxz_xxxxzz, ts_xxxxz_xxxyyz, ts_xxxxz_xxxyz, ts_xxxxz_xxxyzz, ts_xxxxz_xxxzz, ts_xxxxz_xxxzzz, ts_xxxxz_xxyyyz, ts_xxxxz_xxyyz, ts_xxxxz_xxyyzz, ts_xxxxz_xxyzz, ts_xxxxz_xxyzzz, ts_xxxxz_xxzzz, ts_xxxxz_xxzzzz, ts_xxxxz_xyyyyz, ts_xxxxz_xyyyz, ts_xxxxz_xyyyzz, ts_xxxxz_xyyzz, ts_xxxxz_xyyzzz, ts_xxxxz_xyzzz, ts_xxxxz_xyzzzz, ts_xxxxz_xzzzz, ts_xxxxz_xzzzzz, ts_xxxxz_zzzzzz, ts_xxxyz_yyyyyz, ts_xxxyz_yyyyzz, ts_xxxyz_yyyzzz, ts_xxxyz_yyzzzz, ts_xxxyz_yzzzzz, ts_xxyz_yyyyyz, ts_xxyz_yyyyzz, ts_xxyz_yyyzzz, ts_xxyz_yyzzzz, ts_xxyz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyz_xxxxxx[i] = ts_xxxxz_xxxxxx[i] * pa_y[i];

        ts_xxxxyz_xxxxxy[i] = ts_xxxxy_xxxxxy[i] * pa_z[i];

        ts_xxxxyz_xxxxxz[i] = ts_xxxxz_xxxxxz[i] * pa_y[i];

        ts_xxxxyz_xxxxyy[i] = ts_xxxxy_xxxxyy[i] * pa_z[i];

        ts_xxxxyz_xxxxyz[i] = ts_xxxxz_xxxxz[i] * fe_0 + ts_xxxxz_xxxxyz[i] * pa_y[i];

        ts_xxxxyz_xxxxzz[i] = ts_xxxxz_xxxxzz[i] * pa_y[i];

        ts_xxxxyz_xxxyyy[i] = ts_xxxxy_xxxyyy[i] * pa_z[i];

        ts_xxxxyz_xxxyyz[i] = 2.0 * ts_xxxxz_xxxyz[i] * fe_0 + ts_xxxxz_xxxyyz[i] * pa_y[i];

        ts_xxxxyz_xxxyzz[i] = ts_xxxxz_xxxzz[i] * fe_0 + ts_xxxxz_xxxyzz[i] * pa_y[i];

        ts_xxxxyz_xxxzzz[i] = ts_xxxxz_xxxzzz[i] * pa_y[i];

        ts_xxxxyz_xxyyyy[i] = ts_xxxxy_xxyyyy[i] * pa_z[i];

        ts_xxxxyz_xxyyyz[i] = 3.0 * ts_xxxxz_xxyyz[i] * fe_0 + ts_xxxxz_xxyyyz[i] * pa_y[i];

        ts_xxxxyz_xxyyzz[i] = 2.0 * ts_xxxxz_xxyzz[i] * fe_0 + ts_xxxxz_xxyyzz[i] * pa_y[i];

        ts_xxxxyz_xxyzzz[i] = ts_xxxxz_xxzzz[i] * fe_0 + ts_xxxxz_xxyzzz[i] * pa_y[i];

        ts_xxxxyz_xxzzzz[i] = ts_xxxxz_xxzzzz[i] * pa_y[i];

        ts_xxxxyz_xyyyyy[i] = ts_xxxxy_xyyyyy[i] * pa_z[i];

        ts_xxxxyz_xyyyyz[i] = 4.0 * ts_xxxxz_xyyyz[i] * fe_0 + ts_xxxxz_xyyyyz[i] * pa_y[i];

        ts_xxxxyz_xyyyzz[i] = 3.0 * ts_xxxxz_xyyzz[i] * fe_0 + ts_xxxxz_xyyyzz[i] * pa_y[i];

        ts_xxxxyz_xyyzzz[i] = 2.0 * ts_xxxxz_xyzzz[i] * fe_0 + ts_xxxxz_xyyzzz[i] * pa_y[i];

        ts_xxxxyz_xyzzzz[i] = ts_xxxxz_xzzzz[i] * fe_0 + ts_xxxxz_xyzzzz[i] * pa_y[i];

        ts_xxxxyz_xzzzzz[i] = ts_xxxxz_xzzzzz[i] * pa_y[i];

        ts_xxxxyz_yyyyyy[i] = ts_xxxxy_yyyyyy[i] * pa_z[i];

        ts_xxxxyz_yyyyyz[i] = 3.0 * ts_xxyz_yyyyyz[i] * fe_0 + ts_xxxyz_yyyyyz[i] * pa_x[i];

        ts_xxxxyz_yyyyzz[i] = 3.0 * ts_xxyz_yyyyzz[i] * fe_0 + ts_xxxyz_yyyyzz[i] * pa_x[i];

        ts_xxxxyz_yyyzzz[i] = 3.0 * ts_xxyz_yyyzzz[i] * fe_0 + ts_xxxyz_yyyzzz[i] * pa_x[i];

        ts_xxxxyz_yyzzzz[i] = 3.0 * ts_xxyz_yyzzzz[i] * fe_0 + ts_xxxyz_yyzzzz[i] * pa_x[i];

        ts_xxxxyz_yzzzzz[i] = 3.0 * ts_xxyz_yzzzzz[i] * fe_0 + ts_xxxyz_yzzzzz[i] * pa_x[i];

        ts_xxxxyz_zzzzzz[i] = ts_xxxxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xxxxxx, ts_xxxx_xxxxxy, ts_xxxx_xxxxyy, ts_xxxx_xxxyyy, ts_xxxx_xxyyyy, ts_xxxx_xyyyyy, ts_xxxxz_xxxxxx, ts_xxxxz_xxxxxy, ts_xxxxz_xxxxyy, ts_xxxxz_xxxyyy, ts_xxxxz_xxyyyy, ts_xxxxz_xyyyyy, ts_xxxxzz_xxxxxx, ts_xxxxzz_xxxxxy, ts_xxxxzz_xxxxxz, ts_xxxxzz_xxxxyy, ts_xxxxzz_xxxxyz, ts_xxxxzz_xxxxzz, ts_xxxxzz_xxxyyy, ts_xxxxzz_xxxyyz, ts_xxxxzz_xxxyzz, ts_xxxxzz_xxxzzz, ts_xxxxzz_xxyyyy, ts_xxxxzz_xxyyyz, ts_xxxxzz_xxyyzz, ts_xxxxzz_xxyzzz, ts_xxxxzz_xxzzzz, ts_xxxxzz_xyyyyy, ts_xxxxzz_xyyyyz, ts_xxxxzz_xyyyzz, ts_xxxxzz_xyyzzz, ts_xxxxzz_xyzzzz, ts_xxxxzz_xzzzzz, ts_xxxxzz_yyyyyy, ts_xxxxzz_yyyyyz, ts_xxxxzz_yyyyzz, ts_xxxxzz_yyyzzz, ts_xxxxzz_yyzzzz, ts_xxxxzz_yzzzzz, ts_xxxxzz_zzzzzz, ts_xxxzz_xxxxxz, ts_xxxzz_xxxxyz, ts_xxxzz_xxxxz, ts_xxxzz_xxxxzz, ts_xxxzz_xxxyyz, ts_xxxzz_xxxyz, ts_xxxzz_xxxyzz, ts_xxxzz_xxxzz, ts_xxxzz_xxxzzz, ts_xxxzz_xxyyyz, ts_xxxzz_xxyyz, ts_xxxzz_xxyyzz, ts_xxxzz_xxyzz, ts_xxxzz_xxyzzz, ts_xxxzz_xxzzz, ts_xxxzz_xxzzzz, ts_xxxzz_xyyyyz, ts_xxxzz_xyyyz, ts_xxxzz_xyyyzz, ts_xxxzz_xyyzz, ts_xxxzz_xyyzzz, ts_xxxzz_xyzzz, ts_xxxzz_xyzzzz, ts_xxxzz_xzzzz, ts_xxxzz_xzzzzz, ts_xxxzz_yyyyyy, ts_xxxzz_yyyyyz, ts_xxxzz_yyyyz, ts_xxxzz_yyyyzz, ts_xxxzz_yyyzz, ts_xxxzz_yyyzzz, ts_xxxzz_yyzzz, ts_xxxzz_yyzzzz, ts_xxxzz_yzzzz, ts_xxxzz_yzzzzz, ts_xxxzz_zzzzz, ts_xxxzz_zzzzzz, ts_xxzz_xxxxxz, ts_xxzz_xxxxyz, ts_xxzz_xxxxzz, ts_xxzz_xxxyyz, ts_xxzz_xxxyzz, ts_xxzz_xxxzzz, ts_xxzz_xxyyyz, ts_xxzz_xxyyzz, ts_xxzz_xxyzzz, ts_xxzz_xxzzzz, ts_xxzz_xyyyyz, ts_xxzz_xyyyzz, ts_xxzz_xyyzzz, ts_xxzz_xyzzzz, ts_xxzz_xzzzzz, ts_xxzz_yyyyyy, ts_xxzz_yyyyyz, ts_xxzz_yyyyzz, ts_xxzz_yyyzzz, ts_xxzz_yyzzzz, ts_xxzz_yzzzzz, ts_xxzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxzz_xxxxxx[i] = ts_xxxx_xxxxxx[i] * fe_0 + ts_xxxxz_xxxxxx[i] * pa_z[i];

        ts_xxxxzz_xxxxxy[i] = ts_xxxx_xxxxxy[i] * fe_0 + ts_xxxxz_xxxxxy[i] * pa_z[i];

        ts_xxxxzz_xxxxxz[i] = 3.0 * ts_xxzz_xxxxxz[i] * fe_0 + 5.0 * ts_xxxzz_xxxxz[i] * fe_0 + ts_xxxzz_xxxxxz[i] * pa_x[i];

        ts_xxxxzz_xxxxyy[i] = ts_xxxx_xxxxyy[i] * fe_0 + ts_xxxxz_xxxxyy[i] * pa_z[i];

        ts_xxxxzz_xxxxyz[i] = 3.0 * ts_xxzz_xxxxyz[i] * fe_0 + 4.0 * ts_xxxzz_xxxyz[i] * fe_0 + ts_xxxzz_xxxxyz[i] * pa_x[i];

        ts_xxxxzz_xxxxzz[i] = 3.0 * ts_xxzz_xxxxzz[i] * fe_0 + 4.0 * ts_xxxzz_xxxzz[i] * fe_0 + ts_xxxzz_xxxxzz[i] * pa_x[i];

        ts_xxxxzz_xxxyyy[i] = ts_xxxx_xxxyyy[i] * fe_0 + ts_xxxxz_xxxyyy[i] * pa_z[i];

        ts_xxxxzz_xxxyyz[i] = 3.0 * ts_xxzz_xxxyyz[i] * fe_0 + 3.0 * ts_xxxzz_xxyyz[i] * fe_0 + ts_xxxzz_xxxyyz[i] * pa_x[i];

        ts_xxxxzz_xxxyzz[i] = 3.0 * ts_xxzz_xxxyzz[i] * fe_0 + 3.0 * ts_xxxzz_xxyzz[i] * fe_0 + ts_xxxzz_xxxyzz[i] * pa_x[i];

        ts_xxxxzz_xxxzzz[i] = 3.0 * ts_xxzz_xxxzzz[i] * fe_0 + 3.0 * ts_xxxzz_xxzzz[i] * fe_0 + ts_xxxzz_xxxzzz[i] * pa_x[i];

        ts_xxxxzz_xxyyyy[i] = ts_xxxx_xxyyyy[i] * fe_0 + ts_xxxxz_xxyyyy[i] * pa_z[i];

        ts_xxxxzz_xxyyyz[i] = 3.0 * ts_xxzz_xxyyyz[i] * fe_0 + 2.0 * ts_xxxzz_xyyyz[i] * fe_0 + ts_xxxzz_xxyyyz[i] * pa_x[i];

        ts_xxxxzz_xxyyzz[i] = 3.0 * ts_xxzz_xxyyzz[i] * fe_0 + 2.0 * ts_xxxzz_xyyzz[i] * fe_0 + ts_xxxzz_xxyyzz[i] * pa_x[i];

        ts_xxxxzz_xxyzzz[i] = 3.0 * ts_xxzz_xxyzzz[i] * fe_0 + 2.0 * ts_xxxzz_xyzzz[i] * fe_0 + ts_xxxzz_xxyzzz[i] * pa_x[i];

        ts_xxxxzz_xxzzzz[i] = 3.0 * ts_xxzz_xxzzzz[i] * fe_0 + 2.0 * ts_xxxzz_xzzzz[i] * fe_0 + ts_xxxzz_xxzzzz[i] * pa_x[i];

        ts_xxxxzz_xyyyyy[i] = ts_xxxx_xyyyyy[i] * fe_0 + ts_xxxxz_xyyyyy[i] * pa_z[i];

        ts_xxxxzz_xyyyyz[i] = 3.0 * ts_xxzz_xyyyyz[i] * fe_0 + ts_xxxzz_yyyyz[i] * fe_0 + ts_xxxzz_xyyyyz[i] * pa_x[i];

        ts_xxxxzz_xyyyzz[i] = 3.0 * ts_xxzz_xyyyzz[i] * fe_0 + ts_xxxzz_yyyzz[i] * fe_0 + ts_xxxzz_xyyyzz[i] * pa_x[i];

        ts_xxxxzz_xyyzzz[i] = 3.0 * ts_xxzz_xyyzzz[i] * fe_0 + ts_xxxzz_yyzzz[i] * fe_0 + ts_xxxzz_xyyzzz[i] * pa_x[i];

        ts_xxxxzz_xyzzzz[i] = 3.0 * ts_xxzz_xyzzzz[i] * fe_0 + ts_xxxzz_yzzzz[i] * fe_0 + ts_xxxzz_xyzzzz[i] * pa_x[i];

        ts_xxxxzz_xzzzzz[i] = 3.0 * ts_xxzz_xzzzzz[i] * fe_0 + ts_xxxzz_zzzzz[i] * fe_0 + ts_xxxzz_xzzzzz[i] * pa_x[i];

        ts_xxxxzz_yyyyyy[i] = 3.0 * ts_xxzz_yyyyyy[i] * fe_0 + ts_xxxzz_yyyyyy[i] * pa_x[i];

        ts_xxxxzz_yyyyyz[i] = 3.0 * ts_xxzz_yyyyyz[i] * fe_0 + ts_xxxzz_yyyyyz[i] * pa_x[i];

        ts_xxxxzz_yyyyzz[i] = 3.0 * ts_xxzz_yyyyzz[i] * fe_0 + ts_xxxzz_yyyyzz[i] * pa_x[i];

        ts_xxxxzz_yyyzzz[i] = 3.0 * ts_xxzz_yyyzzz[i] * fe_0 + ts_xxxzz_yyyzzz[i] * pa_x[i];

        ts_xxxxzz_yyzzzz[i] = 3.0 * ts_xxzz_yyzzzz[i] * fe_0 + ts_xxxzz_yyzzzz[i] * pa_x[i];

        ts_xxxxzz_yzzzzz[i] = 3.0 * ts_xxzz_yzzzzz[i] * fe_0 + ts_xxxzz_yzzzzz[i] * pa_x[i];

        ts_xxxxzz_zzzzzz[i] = 3.0 * ts_xxzz_zzzzzz[i] * fe_0 + ts_xxxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxy_xxxxxx, ts_xxxy_xxxxxz, ts_xxxy_xxxxzz, ts_xxxy_xxxzzz, ts_xxxy_xxzzzz, ts_xxxy_xzzzzz, ts_xxxyy_xxxxxx, ts_xxxyy_xxxxxz, ts_xxxyy_xxxxzz, ts_xxxyy_xxxzzz, ts_xxxyy_xxzzzz, ts_xxxyy_xzzzzz, ts_xxxyyy_xxxxxx, ts_xxxyyy_xxxxxy, ts_xxxyyy_xxxxxz, ts_xxxyyy_xxxxyy, ts_xxxyyy_xxxxyz, ts_xxxyyy_xxxxzz, ts_xxxyyy_xxxyyy, ts_xxxyyy_xxxyyz, ts_xxxyyy_xxxyzz, ts_xxxyyy_xxxzzz, ts_xxxyyy_xxyyyy, ts_xxxyyy_xxyyyz, ts_xxxyyy_xxyyzz, ts_xxxyyy_xxyzzz, ts_xxxyyy_xxzzzz, ts_xxxyyy_xyyyyy, ts_xxxyyy_xyyyyz, ts_xxxyyy_xyyyzz, ts_xxxyyy_xyyzzz, ts_xxxyyy_xyzzzz, ts_xxxyyy_xzzzzz, ts_xxxyyy_yyyyyy, ts_xxxyyy_yyyyyz, ts_xxxyyy_yyyyzz, ts_xxxyyy_yyyzzz, ts_xxxyyy_yyzzzz, ts_xxxyyy_yzzzzz, ts_xxxyyy_zzzzzz, ts_xxyyy_xxxxxy, ts_xxyyy_xxxxy, ts_xxyyy_xxxxyy, ts_xxyyy_xxxxyz, ts_xxyyy_xxxyy, ts_xxyyy_xxxyyy, ts_xxyyy_xxxyyz, ts_xxyyy_xxxyz, ts_xxyyy_xxxyzz, ts_xxyyy_xxyyy, ts_xxyyy_xxyyyy, ts_xxyyy_xxyyyz, ts_xxyyy_xxyyz, ts_xxyyy_xxyyzz, ts_xxyyy_xxyzz, ts_xxyyy_xxyzzz, ts_xxyyy_xyyyy, ts_xxyyy_xyyyyy, ts_xxyyy_xyyyyz, ts_xxyyy_xyyyz, ts_xxyyy_xyyyzz, ts_xxyyy_xyyzz, ts_xxyyy_xyyzzz, ts_xxyyy_xyzzz, ts_xxyyy_xyzzzz, ts_xxyyy_yyyyy, ts_xxyyy_yyyyyy, ts_xxyyy_yyyyyz, ts_xxyyy_yyyyz, ts_xxyyy_yyyyzz, ts_xxyyy_yyyzz, ts_xxyyy_yyyzzz, ts_xxyyy_yyzzz, ts_xxyyy_yyzzzz, ts_xxyyy_yzzzz, ts_xxyyy_yzzzzz, ts_xxyyy_zzzzzz, ts_xyyy_xxxxxy, ts_xyyy_xxxxyy, ts_xyyy_xxxxyz, ts_xyyy_xxxyyy, ts_xyyy_xxxyyz, ts_xyyy_xxxyzz, ts_xyyy_xxyyyy, ts_xyyy_xxyyyz, ts_xyyy_xxyyzz, ts_xyyy_xxyzzz, ts_xyyy_xyyyyy, ts_xyyy_xyyyyz, ts_xyyy_xyyyzz, ts_xyyy_xyyzzz, ts_xyyy_xyzzzz, ts_xyyy_yyyyyy, ts_xyyy_yyyyyz, ts_xyyy_yyyyzz, ts_xyyy_yyyzzz, ts_xyyy_yyzzzz, ts_xyyy_yzzzzz, ts_xyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyy_xxxxxx[i] = 2.0 * ts_xxxy_xxxxxx[i] * fe_0 + ts_xxxyy_xxxxxx[i] * pa_y[i];

        ts_xxxyyy_xxxxxy[i] = 2.0 * ts_xyyy_xxxxxy[i] * fe_0 + 5.0 * ts_xxyyy_xxxxy[i] * fe_0 + ts_xxyyy_xxxxxy[i] * pa_x[i];

        ts_xxxyyy_xxxxxz[i] = 2.0 * ts_xxxy_xxxxxz[i] * fe_0 + ts_xxxyy_xxxxxz[i] * pa_y[i];

        ts_xxxyyy_xxxxyy[i] = 2.0 * ts_xyyy_xxxxyy[i] * fe_0 + 4.0 * ts_xxyyy_xxxyy[i] * fe_0 + ts_xxyyy_xxxxyy[i] * pa_x[i];

        ts_xxxyyy_xxxxyz[i] = 2.0 * ts_xyyy_xxxxyz[i] * fe_0 + 4.0 * ts_xxyyy_xxxyz[i] * fe_0 + ts_xxyyy_xxxxyz[i] * pa_x[i];

        ts_xxxyyy_xxxxzz[i] = 2.0 * ts_xxxy_xxxxzz[i] * fe_0 + ts_xxxyy_xxxxzz[i] * pa_y[i];

        ts_xxxyyy_xxxyyy[i] = 2.0 * ts_xyyy_xxxyyy[i] * fe_0 + 3.0 * ts_xxyyy_xxyyy[i] * fe_0 + ts_xxyyy_xxxyyy[i] * pa_x[i];

        ts_xxxyyy_xxxyyz[i] = 2.0 * ts_xyyy_xxxyyz[i] * fe_0 + 3.0 * ts_xxyyy_xxyyz[i] * fe_0 + ts_xxyyy_xxxyyz[i] * pa_x[i];

        ts_xxxyyy_xxxyzz[i] = 2.0 * ts_xyyy_xxxyzz[i] * fe_0 + 3.0 * ts_xxyyy_xxyzz[i] * fe_0 + ts_xxyyy_xxxyzz[i] * pa_x[i];

        ts_xxxyyy_xxxzzz[i] = 2.0 * ts_xxxy_xxxzzz[i] * fe_0 + ts_xxxyy_xxxzzz[i] * pa_y[i];

        ts_xxxyyy_xxyyyy[i] = 2.0 * ts_xyyy_xxyyyy[i] * fe_0 + 2.0 * ts_xxyyy_xyyyy[i] * fe_0 + ts_xxyyy_xxyyyy[i] * pa_x[i];

        ts_xxxyyy_xxyyyz[i] = 2.0 * ts_xyyy_xxyyyz[i] * fe_0 + 2.0 * ts_xxyyy_xyyyz[i] * fe_0 + ts_xxyyy_xxyyyz[i] * pa_x[i];

        ts_xxxyyy_xxyyzz[i] = 2.0 * ts_xyyy_xxyyzz[i] * fe_0 + 2.0 * ts_xxyyy_xyyzz[i] * fe_0 + ts_xxyyy_xxyyzz[i] * pa_x[i];

        ts_xxxyyy_xxyzzz[i] = 2.0 * ts_xyyy_xxyzzz[i] * fe_0 + 2.0 * ts_xxyyy_xyzzz[i] * fe_0 + ts_xxyyy_xxyzzz[i] * pa_x[i];

        ts_xxxyyy_xxzzzz[i] = 2.0 * ts_xxxy_xxzzzz[i] * fe_0 + ts_xxxyy_xxzzzz[i] * pa_y[i];

        ts_xxxyyy_xyyyyy[i] = 2.0 * ts_xyyy_xyyyyy[i] * fe_0 + ts_xxyyy_yyyyy[i] * fe_0 + ts_xxyyy_xyyyyy[i] * pa_x[i];

        ts_xxxyyy_xyyyyz[i] = 2.0 * ts_xyyy_xyyyyz[i] * fe_0 + ts_xxyyy_yyyyz[i] * fe_0 + ts_xxyyy_xyyyyz[i] * pa_x[i];

        ts_xxxyyy_xyyyzz[i] = 2.0 * ts_xyyy_xyyyzz[i] * fe_0 + ts_xxyyy_yyyzz[i] * fe_0 + ts_xxyyy_xyyyzz[i] * pa_x[i];

        ts_xxxyyy_xyyzzz[i] = 2.0 * ts_xyyy_xyyzzz[i] * fe_0 + ts_xxyyy_yyzzz[i] * fe_0 + ts_xxyyy_xyyzzz[i] * pa_x[i];

        ts_xxxyyy_xyzzzz[i] = 2.0 * ts_xyyy_xyzzzz[i] * fe_0 + ts_xxyyy_yzzzz[i] * fe_0 + ts_xxyyy_xyzzzz[i] * pa_x[i];

        ts_xxxyyy_xzzzzz[i] = 2.0 * ts_xxxy_xzzzzz[i] * fe_0 + ts_xxxyy_xzzzzz[i] * pa_y[i];

        ts_xxxyyy_yyyyyy[i] = 2.0 * ts_xyyy_yyyyyy[i] * fe_0 + ts_xxyyy_yyyyyy[i] * pa_x[i];

        ts_xxxyyy_yyyyyz[i] = 2.0 * ts_xyyy_yyyyyz[i] * fe_0 + ts_xxyyy_yyyyyz[i] * pa_x[i];

        ts_xxxyyy_yyyyzz[i] = 2.0 * ts_xyyy_yyyyzz[i] * fe_0 + ts_xxyyy_yyyyzz[i] * pa_x[i];

        ts_xxxyyy_yyyzzz[i] = 2.0 * ts_xyyy_yyyzzz[i] * fe_0 + ts_xxyyy_yyyzzz[i] * pa_x[i];

        ts_xxxyyy_yyzzzz[i] = 2.0 * ts_xyyy_yyzzzz[i] * fe_0 + ts_xxyyy_yyzzzz[i] * pa_x[i];

        ts_xxxyyy_yzzzzz[i] = 2.0 * ts_xyyy_yzzzzz[i] * fe_0 + ts_xxyyy_yzzzzz[i] * pa_x[i];

        ts_xxxyyy_zzzzzz[i] = 2.0 * ts_xyyy_zzzzzz[i] * fe_0 + ts_xxyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxyy_xxxxxx, ts_xxxyy_xxxxxy, ts_xxxyy_xxxxy, ts_xxxyy_xxxxyy, ts_xxxyy_xxxxyz, ts_xxxyy_xxxyy, ts_xxxyy_xxxyyy, ts_xxxyy_xxxyyz, ts_xxxyy_xxxyz, ts_xxxyy_xxxyzz, ts_xxxyy_xxyyy, ts_xxxyy_xxyyyy, ts_xxxyy_xxyyyz, ts_xxxyy_xxyyz, ts_xxxyy_xxyyzz, ts_xxxyy_xxyzz, ts_xxxyy_xxyzzz, ts_xxxyy_xyyyy, ts_xxxyy_xyyyyy, ts_xxxyy_xyyyyz, ts_xxxyy_xyyyz, ts_xxxyy_xyyyzz, ts_xxxyy_xyyzz, ts_xxxyy_xyyzzz, ts_xxxyy_xyzzz, ts_xxxyy_xyzzzz, ts_xxxyy_yyyyyy, ts_xxxyyz_xxxxxx, ts_xxxyyz_xxxxxy, ts_xxxyyz_xxxxxz, ts_xxxyyz_xxxxyy, ts_xxxyyz_xxxxyz, ts_xxxyyz_xxxxzz, ts_xxxyyz_xxxyyy, ts_xxxyyz_xxxyyz, ts_xxxyyz_xxxyzz, ts_xxxyyz_xxxzzz, ts_xxxyyz_xxyyyy, ts_xxxyyz_xxyyyz, ts_xxxyyz_xxyyzz, ts_xxxyyz_xxyzzz, ts_xxxyyz_xxzzzz, ts_xxxyyz_xyyyyy, ts_xxxyyz_xyyyyz, ts_xxxyyz_xyyyzz, ts_xxxyyz_xyyzzz, ts_xxxyyz_xyzzzz, ts_xxxyyz_xzzzzz, ts_xxxyyz_yyyyyy, ts_xxxyyz_yyyyyz, ts_xxxyyz_yyyyzz, ts_xxxyyz_yyyzzz, ts_xxxyyz_yyzzzz, ts_xxxyyz_yzzzzz, ts_xxxyyz_zzzzzz, ts_xxxyz_xxxxxz, ts_xxxyz_xxxxzz, ts_xxxyz_xxxzzz, ts_xxxyz_xxzzzz, ts_xxxyz_xzzzzz, ts_xxxz_xxxxxz, ts_xxxz_xxxxzz, ts_xxxz_xxxzzz, ts_xxxz_xxzzzz, ts_xxxz_xzzzzz, ts_xxyyz_yyyyyz, ts_xxyyz_yyyyzz, ts_xxyyz_yyyzzz, ts_xxyyz_yyzzzz, ts_xxyyz_yzzzzz, ts_xxyyz_zzzzzz, ts_xyyz_yyyyyz, ts_xyyz_yyyyzz, ts_xyyz_yyyzzz, ts_xyyz_yyzzzz, ts_xyyz_yzzzzz, ts_xyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyz_xxxxxx[i] = ts_xxxyy_xxxxxx[i] * pa_z[i];

        ts_xxxyyz_xxxxxy[i] = ts_xxxyy_xxxxxy[i] * pa_z[i];

        ts_xxxyyz_xxxxxz[i] = ts_xxxz_xxxxxz[i] * fe_0 + ts_xxxyz_xxxxxz[i] * pa_y[i];

        ts_xxxyyz_xxxxyy[i] = ts_xxxyy_xxxxyy[i] * pa_z[i];

        ts_xxxyyz_xxxxyz[i] = ts_xxxyy_xxxxy[i] * fe_0 + ts_xxxyy_xxxxyz[i] * pa_z[i];

        ts_xxxyyz_xxxxzz[i] = ts_xxxz_xxxxzz[i] * fe_0 + ts_xxxyz_xxxxzz[i] * pa_y[i];

        ts_xxxyyz_xxxyyy[i] = ts_xxxyy_xxxyyy[i] * pa_z[i];

        ts_xxxyyz_xxxyyz[i] = ts_xxxyy_xxxyy[i] * fe_0 + ts_xxxyy_xxxyyz[i] * pa_z[i];

        ts_xxxyyz_xxxyzz[i] = 2.0 * ts_xxxyy_xxxyz[i] * fe_0 + ts_xxxyy_xxxyzz[i] * pa_z[i];

        ts_xxxyyz_xxxzzz[i] = ts_xxxz_xxxzzz[i] * fe_0 + ts_xxxyz_xxxzzz[i] * pa_y[i];

        ts_xxxyyz_xxyyyy[i] = ts_xxxyy_xxyyyy[i] * pa_z[i];

        ts_xxxyyz_xxyyyz[i] = ts_xxxyy_xxyyy[i] * fe_0 + ts_xxxyy_xxyyyz[i] * pa_z[i];

        ts_xxxyyz_xxyyzz[i] = 2.0 * ts_xxxyy_xxyyz[i] * fe_0 + ts_xxxyy_xxyyzz[i] * pa_z[i];

        ts_xxxyyz_xxyzzz[i] = 3.0 * ts_xxxyy_xxyzz[i] * fe_0 + ts_xxxyy_xxyzzz[i] * pa_z[i];

        ts_xxxyyz_xxzzzz[i] = ts_xxxz_xxzzzz[i] * fe_0 + ts_xxxyz_xxzzzz[i] * pa_y[i];

        ts_xxxyyz_xyyyyy[i] = ts_xxxyy_xyyyyy[i] * pa_z[i];

        ts_xxxyyz_xyyyyz[i] = ts_xxxyy_xyyyy[i] * fe_0 + ts_xxxyy_xyyyyz[i] * pa_z[i];

        ts_xxxyyz_xyyyzz[i] = 2.0 * ts_xxxyy_xyyyz[i] * fe_0 + ts_xxxyy_xyyyzz[i] * pa_z[i];

        ts_xxxyyz_xyyzzz[i] = 3.0 * ts_xxxyy_xyyzz[i] * fe_0 + ts_xxxyy_xyyzzz[i] * pa_z[i];

        ts_xxxyyz_xyzzzz[i] = 4.0 * ts_xxxyy_xyzzz[i] * fe_0 + ts_xxxyy_xyzzzz[i] * pa_z[i];

        ts_xxxyyz_xzzzzz[i] = ts_xxxz_xzzzzz[i] * fe_0 + ts_xxxyz_xzzzzz[i] * pa_y[i];

        ts_xxxyyz_yyyyyy[i] = ts_xxxyy_yyyyyy[i] * pa_z[i];

        ts_xxxyyz_yyyyyz[i] = 2.0 * ts_xyyz_yyyyyz[i] * fe_0 + ts_xxyyz_yyyyyz[i] * pa_x[i];

        ts_xxxyyz_yyyyzz[i] = 2.0 * ts_xyyz_yyyyzz[i] * fe_0 + ts_xxyyz_yyyyzz[i] * pa_x[i];

        ts_xxxyyz_yyyzzz[i] = 2.0 * ts_xyyz_yyyzzz[i] * fe_0 + ts_xxyyz_yyyzzz[i] * pa_x[i];

        ts_xxxyyz_yyzzzz[i] = 2.0 * ts_xyyz_yyzzzz[i] * fe_0 + ts_xxyyz_yyzzzz[i] * pa_x[i];

        ts_xxxyyz_yzzzzz[i] = 2.0 * ts_xyyz_yzzzzz[i] * fe_0 + ts_xxyyz_yzzzzz[i] * pa_x[i];

        ts_xxxyyz_zzzzzz[i] = 2.0 * ts_xyyz_zzzzzz[i] * fe_0 + ts_xxyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 224-252 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxyzz_xxxxxx, ts_xxxyzz_xxxxxy, ts_xxxyzz_xxxxxz, ts_xxxyzz_xxxxyy, ts_xxxyzz_xxxxyz, ts_xxxyzz_xxxxzz, ts_xxxyzz_xxxyyy, ts_xxxyzz_xxxyyz, ts_xxxyzz_xxxyzz, ts_xxxyzz_xxxzzz, ts_xxxyzz_xxyyyy, ts_xxxyzz_xxyyyz, ts_xxxyzz_xxyyzz, ts_xxxyzz_xxyzzz, ts_xxxyzz_xxzzzz, ts_xxxyzz_xyyyyy, ts_xxxyzz_xyyyyz, ts_xxxyzz_xyyyzz, ts_xxxyzz_xyyzzz, ts_xxxyzz_xyzzzz, ts_xxxyzz_xzzzzz, ts_xxxyzz_yyyyyy, ts_xxxyzz_yyyyyz, ts_xxxyzz_yyyyzz, ts_xxxyzz_yyyzzz, ts_xxxyzz_yyzzzz, ts_xxxyzz_yzzzzz, ts_xxxyzz_zzzzzz, ts_xxxzz_xxxxx, ts_xxxzz_xxxxxx, ts_xxxzz_xxxxxy, ts_xxxzz_xxxxxz, ts_xxxzz_xxxxy, ts_xxxzz_xxxxyy, ts_xxxzz_xxxxyz, ts_xxxzz_xxxxz, ts_xxxzz_xxxxzz, ts_xxxzz_xxxyy, ts_xxxzz_xxxyyy, ts_xxxzz_xxxyyz, ts_xxxzz_xxxyz, ts_xxxzz_xxxyzz, ts_xxxzz_xxxzz, ts_xxxzz_xxxzzz, ts_xxxzz_xxyyy, ts_xxxzz_xxyyyy, ts_xxxzz_xxyyyz, ts_xxxzz_xxyyz, ts_xxxzz_xxyyzz, ts_xxxzz_xxyzz, ts_xxxzz_xxyzzz, ts_xxxzz_xxzzz, ts_xxxzz_xxzzzz, ts_xxxzz_xyyyy, ts_xxxzz_xyyyyy, ts_xxxzz_xyyyyz, ts_xxxzz_xyyyz, ts_xxxzz_xyyyzz, ts_xxxzz_xyyzz, ts_xxxzz_xyyzzz, ts_xxxzz_xyzzz, ts_xxxzz_xyzzzz, ts_xxxzz_xzzzz, ts_xxxzz_xzzzzz, ts_xxxzz_zzzzzz, ts_xxyzz_yyyyyy, ts_xxyzz_yyyyyz, ts_xxyzz_yyyyzz, ts_xxyzz_yyyzzz, ts_xxyzz_yyzzzz, ts_xxyzz_yzzzzz, ts_xyzz_yyyyyy, ts_xyzz_yyyyyz, ts_xyzz_yyyyzz, ts_xyzz_yyyzzz, ts_xyzz_yyzzzz, ts_xyzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyzz_xxxxxx[i] = ts_xxxzz_xxxxxx[i] * pa_y[i];

        ts_xxxyzz_xxxxxy[i] = ts_xxxzz_xxxxx[i] * fe_0 + ts_xxxzz_xxxxxy[i] * pa_y[i];

        ts_xxxyzz_xxxxxz[i] = ts_xxxzz_xxxxxz[i] * pa_y[i];

        ts_xxxyzz_xxxxyy[i] = 2.0 * ts_xxxzz_xxxxy[i] * fe_0 + ts_xxxzz_xxxxyy[i] * pa_y[i];

        ts_xxxyzz_xxxxyz[i] = ts_xxxzz_xxxxz[i] * fe_0 + ts_xxxzz_xxxxyz[i] * pa_y[i];

        ts_xxxyzz_xxxxzz[i] = ts_xxxzz_xxxxzz[i] * pa_y[i];

        ts_xxxyzz_xxxyyy[i] = 3.0 * ts_xxxzz_xxxyy[i] * fe_0 + ts_xxxzz_xxxyyy[i] * pa_y[i];

        ts_xxxyzz_xxxyyz[i] = 2.0 * ts_xxxzz_xxxyz[i] * fe_0 + ts_xxxzz_xxxyyz[i] * pa_y[i];

        ts_xxxyzz_xxxyzz[i] = ts_xxxzz_xxxzz[i] * fe_0 + ts_xxxzz_xxxyzz[i] * pa_y[i];

        ts_xxxyzz_xxxzzz[i] = ts_xxxzz_xxxzzz[i] * pa_y[i];

        ts_xxxyzz_xxyyyy[i] = 4.0 * ts_xxxzz_xxyyy[i] * fe_0 + ts_xxxzz_xxyyyy[i] * pa_y[i];

        ts_xxxyzz_xxyyyz[i] = 3.0 * ts_xxxzz_xxyyz[i] * fe_0 + ts_xxxzz_xxyyyz[i] * pa_y[i];

        ts_xxxyzz_xxyyzz[i] = 2.0 * ts_xxxzz_xxyzz[i] * fe_0 + ts_xxxzz_xxyyzz[i] * pa_y[i];

        ts_xxxyzz_xxyzzz[i] = ts_xxxzz_xxzzz[i] * fe_0 + ts_xxxzz_xxyzzz[i] * pa_y[i];

        ts_xxxyzz_xxzzzz[i] = ts_xxxzz_xxzzzz[i] * pa_y[i];

        ts_xxxyzz_xyyyyy[i] = 5.0 * ts_xxxzz_xyyyy[i] * fe_0 + ts_xxxzz_xyyyyy[i] * pa_y[i];

        ts_xxxyzz_xyyyyz[i] = 4.0 * ts_xxxzz_xyyyz[i] * fe_0 + ts_xxxzz_xyyyyz[i] * pa_y[i];

        ts_xxxyzz_xyyyzz[i] = 3.0 * ts_xxxzz_xyyzz[i] * fe_0 + ts_xxxzz_xyyyzz[i] * pa_y[i];

        ts_xxxyzz_xyyzzz[i] = 2.0 * ts_xxxzz_xyzzz[i] * fe_0 + ts_xxxzz_xyyzzz[i] * pa_y[i];

        ts_xxxyzz_xyzzzz[i] = ts_xxxzz_xzzzz[i] * fe_0 + ts_xxxzz_xyzzzz[i] * pa_y[i];

        ts_xxxyzz_xzzzzz[i] = ts_xxxzz_xzzzzz[i] * pa_y[i];

        ts_xxxyzz_yyyyyy[i] = 2.0 * ts_xyzz_yyyyyy[i] * fe_0 + ts_xxyzz_yyyyyy[i] * pa_x[i];

        ts_xxxyzz_yyyyyz[i] = 2.0 * ts_xyzz_yyyyyz[i] * fe_0 + ts_xxyzz_yyyyyz[i] * pa_x[i];

        ts_xxxyzz_yyyyzz[i] = 2.0 * ts_xyzz_yyyyzz[i] * fe_0 + ts_xxyzz_yyyyzz[i] * pa_x[i];

        ts_xxxyzz_yyyzzz[i] = 2.0 * ts_xyzz_yyyzzz[i] * fe_0 + ts_xxyzz_yyyzzz[i] * pa_x[i];

        ts_xxxyzz_yyzzzz[i] = 2.0 * ts_xyzz_yyzzzz[i] * fe_0 + ts_xxyzz_yyzzzz[i] * pa_x[i];

        ts_xxxyzz_yzzzzz[i] = 2.0 * ts_xyzz_yzzzzz[i] * fe_0 + ts_xxyzz_yzzzzz[i] * pa_x[i];

        ts_xxxyzz_zzzzzz[i] = ts_xxxzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxz_xxxxxx, ts_xxxz_xxxxxy, ts_xxxz_xxxxyy, ts_xxxz_xxxyyy, ts_xxxz_xxyyyy, ts_xxxz_xyyyyy, ts_xxxzz_xxxxxx, ts_xxxzz_xxxxxy, ts_xxxzz_xxxxyy, ts_xxxzz_xxxyyy, ts_xxxzz_xxyyyy, ts_xxxzz_xyyyyy, ts_xxxzzz_xxxxxx, ts_xxxzzz_xxxxxy, ts_xxxzzz_xxxxxz, ts_xxxzzz_xxxxyy, ts_xxxzzz_xxxxyz, ts_xxxzzz_xxxxzz, ts_xxxzzz_xxxyyy, ts_xxxzzz_xxxyyz, ts_xxxzzz_xxxyzz, ts_xxxzzz_xxxzzz, ts_xxxzzz_xxyyyy, ts_xxxzzz_xxyyyz, ts_xxxzzz_xxyyzz, ts_xxxzzz_xxyzzz, ts_xxxzzz_xxzzzz, ts_xxxzzz_xyyyyy, ts_xxxzzz_xyyyyz, ts_xxxzzz_xyyyzz, ts_xxxzzz_xyyzzz, ts_xxxzzz_xyzzzz, ts_xxxzzz_xzzzzz, ts_xxxzzz_yyyyyy, ts_xxxzzz_yyyyyz, ts_xxxzzz_yyyyzz, ts_xxxzzz_yyyzzz, ts_xxxzzz_yyzzzz, ts_xxxzzz_yzzzzz, ts_xxxzzz_zzzzzz, ts_xxzzz_xxxxxz, ts_xxzzz_xxxxyz, ts_xxzzz_xxxxz, ts_xxzzz_xxxxzz, ts_xxzzz_xxxyyz, ts_xxzzz_xxxyz, ts_xxzzz_xxxyzz, ts_xxzzz_xxxzz, ts_xxzzz_xxxzzz, ts_xxzzz_xxyyyz, ts_xxzzz_xxyyz, ts_xxzzz_xxyyzz, ts_xxzzz_xxyzz, ts_xxzzz_xxyzzz, ts_xxzzz_xxzzz, ts_xxzzz_xxzzzz, ts_xxzzz_xyyyyz, ts_xxzzz_xyyyz, ts_xxzzz_xyyyzz, ts_xxzzz_xyyzz, ts_xxzzz_xyyzzz, ts_xxzzz_xyzzz, ts_xxzzz_xyzzzz, ts_xxzzz_xzzzz, ts_xxzzz_xzzzzz, ts_xxzzz_yyyyyy, ts_xxzzz_yyyyyz, ts_xxzzz_yyyyz, ts_xxzzz_yyyyzz, ts_xxzzz_yyyzz, ts_xxzzz_yyyzzz, ts_xxzzz_yyzzz, ts_xxzzz_yyzzzz, ts_xxzzz_yzzzz, ts_xxzzz_yzzzzz, ts_xxzzz_zzzzz, ts_xxzzz_zzzzzz, ts_xzzz_xxxxxz, ts_xzzz_xxxxyz, ts_xzzz_xxxxzz, ts_xzzz_xxxyyz, ts_xzzz_xxxyzz, ts_xzzz_xxxzzz, ts_xzzz_xxyyyz, ts_xzzz_xxyyzz, ts_xzzz_xxyzzz, ts_xzzz_xxzzzz, ts_xzzz_xyyyyz, ts_xzzz_xyyyzz, ts_xzzz_xyyzzz, ts_xzzz_xyzzzz, ts_xzzz_xzzzzz, ts_xzzz_yyyyyy, ts_xzzz_yyyyyz, ts_xzzz_yyyyzz, ts_xzzz_yyyzzz, ts_xzzz_yyzzzz, ts_xzzz_yzzzzz, ts_xzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzzz_xxxxxx[i] = 2.0 * ts_xxxz_xxxxxx[i] * fe_0 + ts_xxxzz_xxxxxx[i] * pa_z[i];

        ts_xxxzzz_xxxxxy[i] = 2.0 * ts_xxxz_xxxxxy[i] * fe_0 + ts_xxxzz_xxxxxy[i] * pa_z[i];

        ts_xxxzzz_xxxxxz[i] = 2.0 * ts_xzzz_xxxxxz[i] * fe_0 + 5.0 * ts_xxzzz_xxxxz[i] * fe_0 + ts_xxzzz_xxxxxz[i] * pa_x[i];

        ts_xxxzzz_xxxxyy[i] = 2.0 * ts_xxxz_xxxxyy[i] * fe_0 + ts_xxxzz_xxxxyy[i] * pa_z[i];

        ts_xxxzzz_xxxxyz[i] = 2.0 * ts_xzzz_xxxxyz[i] * fe_0 + 4.0 * ts_xxzzz_xxxyz[i] * fe_0 + ts_xxzzz_xxxxyz[i] * pa_x[i];

        ts_xxxzzz_xxxxzz[i] = 2.0 * ts_xzzz_xxxxzz[i] * fe_0 + 4.0 * ts_xxzzz_xxxzz[i] * fe_0 + ts_xxzzz_xxxxzz[i] * pa_x[i];

        ts_xxxzzz_xxxyyy[i] = 2.0 * ts_xxxz_xxxyyy[i] * fe_0 + ts_xxxzz_xxxyyy[i] * pa_z[i];

        ts_xxxzzz_xxxyyz[i] = 2.0 * ts_xzzz_xxxyyz[i] * fe_0 + 3.0 * ts_xxzzz_xxyyz[i] * fe_0 + ts_xxzzz_xxxyyz[i] * pa_x[i];

        ts_xxxzzz_xxxyzz[i] = 2.0 * ts_xzzz_xxxyzz[i] * fe_0 + 3.0 * ts_xxzzz_xxyzz[i] * fe_0 + ts_xxzzz_xxxyzz[i] * pa_x[i];

        ts_xxxzzz_xxxzzz[i] = 2.0 * ts_xzzz_xxxzzz[i] * fe_0 + 3.0 * ts_xxzzz_xxzzz[i] * fe_0 + ts_xxzzz_xxxzzz[i] * pa_x[i];

        ts_xxxzzz_xxyyyy[i] = 2.0 * ts_xxxz_xxyyyy[i] * fe_0 + ts_xxxzz_xxyyyy[i] * pa_z[i];

        ts_xxxzzz_xxyyyz[i] = 2.0 * ts_xzzz_xxyyyz[i] * fe_0 + 2.0 * ts_xxzzz_xyyyz[i] * fe_0 + ts_xxzzz_xxyyyz[i] * pa_x[i];

        ts_xxxzzz_xxyyzz[i] = 2.0 * ts_xzzz_xxyyzz[i] * fe_0 + 2.0 * ts_xxzzz_xyyzz[i] * fe_0 + ts_xxzzz_xxyyzz[i] * pa_x[i];

        ts_xxxzzz_xxyzzz[i] = 2.0 * ts_xzzz_xxyzzz[i] * fe_0 + 2.0 * ts_xxzzz_xyzzz[i] * fe_0 + ts_xxzzz_xxyzzz[i] * pa_x[i];

        ts_xxxzzz_xxzzzz[i] = 2.0 * ts_xzzz_xxzzzz[i] * fe_0 + 2.0 * ts_xxzzz_xzzzz[i] * fe_0 + ts_xxzzz_xxzzzz[i] * pa_x[i];

        ts_xxxzzz_xyyyyy[i] = 2.0 * ts_xxxz_xyyyyy[i] * fe_0 + ts_xxxzz_xyyyyy[i] * pa_z[i];

        ts_xxxzzz_xyyyyz[i] = 2.0 * ts_xzzz_xyyyyz[i] * fe_0 + ts_xxzzz_yyyyz[i] * fe_0 + ts_xxzzz_xyyyyz[i] * pa_x[i];

        ts_xxxzzz_xyyyzz[i] = 2.0 * ts_xzzz_xyyyzz[i] * fe_0 + ts_xxzzz_yyyzz[i] * fe_0 + ts_xxzzz_xyyyzz[i] * pa_x[i];

        ts_xxxzzz_xyyzzz[i] = 2.0 * ts_xzzz_xyyzzz[i] * fe_0 + ts_xxzzz_yyzzz[i] * fe_0 + ts_xxzzz_xyyzzz[i] * pa_x[i];

        ts_xxxzzz_xyzzzz[i] = 2.0 * ts_xzzz_xyzzzz[i] * fe_0 + ts_xxzzz_yzzzz[i] * fe_0 + ts_xxzzz_xyzzzz[i] * pa_x[i];

        ts_xxxzzz_xzzzzz[i] = 2.0 * ts_xzzz_xzzzzz[i] * fe_0 + ts_xxzzz_zzzzz[i] * fe_0 + ts_xxzzz_xzzzzz[i] * pa_x[i];

        ts_xxxzzz_yyyyyy[i] = 2.0 * ts_xzzz_yyyyyy[i] * fe_0 + ts_xxzzz_yyyyyy[i] * pa_x[i];

        ts_xxxzzz_yyyyyz[i] = 2.0 * ts_xzzz_yyyyyz[i] * fe_0 + ts_xxzzz_yyyyyz[i] * pa_x[i];

        ts_xxxzzz_yyyyzz[i] = 2.0 * ts_xzzz_yyyyzz[i] * fe_0 + ts_xxzzz_yyyyzz[i] * pa_x[i];

        ts_xxxzzz_yyyzzz[i] = 2.0 * ts_xzzz_yyyzzz[i] * fe_0 + ts_xxzzz_yyyzzz[i] * pa_x[i];

        ts_xxxzzz_yyzzzz[i] = 2.0 * ts_xzzz_yyzzzz[i] * fe_0 + ts_xxzzz_yyzzzz[i] * pa_x[i];

        ts_xxxzzz_yzzzzz[i] = 2.0 * ts_xzzz_yzzzzz[i] * fe_0 + ts_xxzzz_yzzzzz[i] * pa_x[i];

        ts_xxxzzz_zzzzzz[i] = 2.0 * ts_xzzz_zzzzzz[i] * fe_0 + ts_xxzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyy_xxxxxx, ts_xxyy_xxxxxz, ts_xxyy_xxxxzz, ts_xxyy_xxxzzz, ts_xxyy_xxzzzz, ts_xxyy_xzzzzz, ts_xxyyy_xxxxxx, ts_xxyyy_xxxxxz, ts_xxyyy_xxxxzz, ts_xxyyy_xxxzzz, ts_xxyyy_xxzzzz, ts_xxyyy_xzzzzz, ts_xxyyyy_xxxxxx, ts_xxyyyy_xxxxxy, ts_xxyyyy_xxxxxz, ts_xxyyyy_xxxxyy, ts_xxyyyy_xxxxyz, ts_xxyyyy_xxxxzz, ts_xxyyyy_xxxyyy, ts_xxyyyy_xxxyyz, ts_xxyyyy_xxxyzz, ts_xxyyyy_xxxzzz, ts_xxyyyy_xxyyyy, ts_xxyyyy_xxyyyz, ts_xxyyyy_xxyyzz, ts_xxyyyy_xxyzzz, ts_xxyyyy_xxzzzz, ts_xxyyyy_xyyyyy, ts_xxyyyy_xyyyyz, ts_xxyyyy_xyyyzz, ts_xxyyyy_xyyzzz, ts_xxyyyy_xyzzzz, ts_xxyyyy_xzzzzz, ts_xxyyyy_yyyyyy, ts_xxyyyy_yyyyyz, ts_xxyyyy_yyyyzz, ts_xxyyyy_yyyzzz, ts_xxyyyy_yyzzzz, ts_xxyyyy_yzzzzz, ts_xxyyyy_zzzzzz, ts_xyyyy_xxxxxy, ts_xyyyy_xxxxy, ts_xyyyy_xxxxyy, ts_xyyyy_xxxxyz, ts_xyyyy_xxxyy, ts_xyyyy_xxxyyy, ts_xyyyy_xxxyyz, ts_xyyyy_xxxyz, ts_xyyyy_xxxyzz, ts_xyyyy_xxyyy, ts_xyyyy_xxyyyy, ts_xyyyy_xxyyyz, ts_xyyyy_xxyyz, ts_xyyyy_xxyyzz, ts_xyyyy_xxyzz, ts_xyyyy_xxyzzz, ts_xyyyy_xyyyy, ts_xyyyy_xyyyyy, ts_xyyyy_xyyyyz, ts_xyyyy_xyyyz, ts_xyyyy_xyyyzz, ts_xyyyy_xyyzz, ts_xyyyy_xyyzzz, ts_xyyyy_xyzzz, ts_xyyyy_xyzzzz, ts_xyyyy_yyyyy, ts_xyyyy_yyyyyy, ts_xyyyy_yyyyyz, ts_xyyyy_yyyyz, ts_xyyyy_yyyyzz, ts_xyyyy_yyyzz, ts_xyyyy_yyyzzz, ts_xyyyy_yyzzz, ts_xyyyy_yyzzzz, ts_xyyyy_yzzzz, ts_xyyyy_yzzzzz, ts_xyyyy_zzzzzz, ts_yyyy_xxxxxy, ts_yyyy_xxxxyy, ts_yyyy_xxxxyz, ts_yyyy_xxxyyy, ts_yyyy_xxxyyz, ts_yyyy_xxxyzz, ts_yyyy_xxyyyy, ts_yyyy_xxyyyz, ts_yyyy_xxyyzz, ts_yyyy_xxyzzz, ts_yyyy_xyyyyy, ts_yyyy_xyyyyz, ts_yyyy_xyyyzz, ts_yyyy_xyyzzz, ts_yyyy_xyzzzz, ts_yyyy_yyyyyy, ts_yyyy_yyyyyz, ts_yyyy_yyyyzz, ts_yyyy_yyyzzz, ts_yyyy_yyzzzz, ts_yyyy_yzzzzz, ts_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyy_xxxxxx[i] = 3.0 * ts_xxyy_xxxxxx[i] * fe_0 + ts_xxyyy_xxxxxx[i] * pa_y[i];

        ts_xxyyyy_xxxxxy[i] = ts_yyyy_xxxxxy[i] * fe_0 + 5.0 * ts_xyyyy_xxxxy[i] * fe_0 + ts_xyyyy_xxxxxy[i] * pa_x[i];

        ts_xxyyyy_xxxxxz[i] = 3.0 * ts_xxyy_xxxxxz[i] * fe_0 + ts_xxyyy_xxxxxz[i] * pa_y[i];

        ts_xxyyyy_xxxxyy[i] = ts_yyyy_xxxxyy[i] * fe_0 + 4.0 * ts_xyyyy_xxxyy[i] * fe_0 + ts_xyyyy_xxxxyy[i] * pa_x[i];

        ts_xxyyyy_xxxxyz[i] = ts_yyyy_xxxxyz[i] * fe_0 + 4.0 * ts_xyyyy_xxxyz[i] * fe_0 + ts_xyyyy_xxxxyz[i] * pa_x[i];

        ts_xxyyyy_xxxxzz[i] = 3.0 * ts_xxyy_xxxxzz[i] * fe_0 + ts_xxyyy_xxxxzz[i] * pa_y[i];

        ts_xxyyyy_xxxyyy[i] = ts_yyyy_xxxyyy[i] * fe_0 + 3.0 * ts_xyyyy_xxyyy[i] * fe_0 + ts_xyyyy_xxxyyy[i] * pa_x[i];

        ts_xxyyyy_xxxyyz[i] = ts_yyyy_xxxyyz[i] * fe_0 + 3.0 * ts_xyyyy_xxyyz[i] * fe_0 + ts_xyyyy_xxxyyz[i] * pa_x[i];

        ts_xxyyyy_xxxyzz[i] = ts_yyyy_xxxyzz[i] * fe_0 + 3.0 * ts_xyyyy_xxyzz[i] * fe_0 + ts_xyyyy_xxxyzz[i] * pa_x[i];

        ts_xxyyyy_xxxzzz[i] = 3.0 * ts_xxyy_xxxzzz[i] * fe_0 + ts_xxyyy_xxxzzz[i] * pa_y[i];

        ts_xxyyyy_xxyyyy[i] = ts_yyyy_xxyyyy[i] * fe_0 + 2.0 * ts_xyyyy_xyyyy[i] * fe_0 + ts_xyyyy_xxyyyy[i] * pa_x[i];

        ts_xxyyyy_xxyyyz[i] = ts_yyyy_xxyyyz[i] * fe_0 + 2.0 * ts_xyyyy_xyyyz[i] * fe_0 + ts_xyyyy_xxyyyz[i] * pa_x[i];

        ts_xxyyyy_xxyyzz[i] = ts_yyyy_xxyyzz[i] * fe_0 + 2.0 * ts_xyyyy_xyyzz[i] * fe_0 + ts_xyyyy_xxyyzz[i] * pa_x[i];

        ts_xxyyyy_xxyzzz[i] = ts_yyyy_xxyzzz[i] * fe_0 + 2.0 * ts_xyyyy_xyzzz[i] * fe_0 + ts_xyyyy_xxyzzz[i] * pa_x[i];

        ts_xxyyyy_xxzzzz[i] = 3.0 * ts_xxyy_xxzzzz[i] * fe_0 + ts_xxyyy_xxzzzz[i] * pa_y[i];

        ts_xxyyyy_xyyyyy[i] = ts_yyyy_xyyyyy[i] * fe_0 + ts_xyyyy_yyyyy[i] * fe_0 + ts_xyyyy_xyyyyy[i] * pa_x[i];

        ts_xxyyyy_xyyyyz[i] = ts_yyyy_xyyyyz[i] * fe_0 + ts_xyyyy_yyyyz[i] * fe_0 + ts_xyyyy_xyyyyz[i] * pa_x[i];

        ts_xxyyyy_xyyyzz[i] = ts_yyyy_xyyyzz[i] * fe_0 + ts_xyyyy_yyyzz[i] * fe_0 + ts_xyyyy_xyyyzz[i] * pa_x[i];

        ts_xxyyyy_xyyzzz[i] = ts_yyyy_xyyzzz[i] * fe_0 + ts_xyyyy_yyzzz[i] * fe_0 + ts_xyyyy_xyyzzz[i] * pa_x[i];

        ts_xxyyyy_xyzzzz[i] = ts_yyyy_xyzzzz[i] * fe_0 + ts_xyyyy_yzzzz[i] * fe_0 + ts_xyyyy_xyzzzz[i] * pa_x[i];

        ts_xxyyyy_xzzzzz[i] = 3.0 * ts_xxyy_xzzzzz[i] * fe_0 + ts_xxyyy_xzzzzz[i] * pa_y[i];

        ts_xxyyyy_yyyyyy[i] = ts_yyyy_yyyyyy[i] * fe_0 + ts_xyyyy_yyyyyy[i] * pa_x[i];

        ts_xxyyyy_yyyyyz[i] = ts_yyyy_yyyyyz[i] * fe_0 + ts_xyyyy_yyyyyz[i] * pa_x[i];

        ts_xxyyyy_yyyyzz[i] = ts_yyyy_yyyyzz[i] * fe_0 + ts_xyyyy_yyyyzz[i] * pa_x[i];

        ts_xxyyyy_yyyzzz[i] = ts_yyyy_yyyzzz[i] * fe_0 + ts_xyyyy_yyyzzz[i] * pa_x[i];

        ts_xxyyyy_yyzzzz[i] = ts_yyyy_yyzzzz[i] * fe_0 + ts_xyyyy_yyzzzz[i] * pa_x[i];

        ts_xxyyyy_yzzzzz[i] = ts_yyyy_yzzzzz[i] * fe_0 + ts_xyyyy_yzzzzz[i] * pa_x[i];

        ts_xxyyyy_zzzzzz[i] = ts_yyyy_zzzzzz[i] * fe_0 + ts_xyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyyy_xxxxxx, ts_xxyyy_xxxxxy, ts_xxyyy_xxxxy, ts_xxyyy_xxxxyy, ts_xxyyy_xxxxyz, ts_xxyyy_xxxyy, ts_xxyyy_xxxyyy, ts_xxyyy_xxxyyz, ts_xxyyy_xxxyz, ts_xxyyy_xxxyzz, ts_xxyyy_xxyyy, ts_xxyyy_xxyyyy, ts_xxyyy_xxyyyz, ts_xxyyy_xxyyz, ts_xxyyy_xxyyzz, ts_xxyyy_xxyzz, ts_xxyyy_xxyzzz, ts_xxyyy_xyyyy, ts_xxyyy_xyyyyy, ts_xxyyy_xyyyyz, ts_xxyyy_xyyyz, ts_xxyyy_xyyyzz, ts_xxyyy_xyyzz, ts_xxyyy_xyyzzz, ts_xxyyy_xyzzz, ts_xxyyy_xyzzzz, ts_xxyyy_yyyyyy, ts_xxyyyz_xxxxxx, ts_xxyyyz_xxxxxy, ts_xxyyyz_xxxxxz, ts_xxyyyz_xxxxyy, ts_xxyyyz_xxxxyz, ts_xxyyyz_xxxxzz, ts_xxyyyz_xxxyyy, ts_xxyyyz_xxxyyz, ts_xxyyyz_xxxyzz, ts_xxyyyz_xxxzzz, ts_xxyyyz_xxyyyy, ts_xxyyyz_xxyyyz, ts_xxyyyz_xxyyzz, ts_xxyyyz_xxyzzz, ts_xxyyyz_xxzzzz, ts_xxyyyz_xyyyyy, ts_xxyyyz_xyyyyz, ts_xxyyyz_xyyyzz, ts_xxyyyz_xyyzzz, ts_xxyyyz_xyzzzz, ts_xxyyyz_xzzzzz, ts_xxyyyz_yyyyyy, ts_xxyyyz_yyyyyz, ts_xxyyyz_yyyyzz, ts_xxyyyz_yyyzzz, ts_xxyyyz_yyzzzz, ts_xxyyyz_yzzzzz, ts_xxyyyz_zzzzzz, ts_xxyyz_xxxxxz, ts_xxyyz_xxxxzz, ts_xxyyz_xxxzzz, ts_xxyyz_xxzzzz, ts_xxyyz_xzzzzz, ts_xxyz_xxxxxz, ts_xxyz_xxxxzz, ts_xxyz_xxxzzz, ts_xxyz_xxzzzz, ts_xxyz_xzzzzz, ts_xyyyz_yyyyyz, ts_xyyyz_yyyyzz, ts_xyyyz_yyyzzz, ts_xyyyz_yyzzzz, ts_xyyyz_yzzzzz, ts_xyyyz_zzzzzz, ts_yyyz_yyyyyz, ts_yyyz_yyyyzz, ts_yyyz_yyyzzz, ts_yyyz_yyzzzz, ts_yyyz_yzzzzz, ts_yyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyz_xxxxxx[i] = ts_xxyyy_xxxxxx[i] * pa_z[i];

        ts_xxyyyz_xxxxxy[i] = ts_xxyyy_xxxxxy[i] * pa_z[i];

        ts_xxyyyz_xxxxxz[i] = 2.0 * ts_xxyz_xxxxxz[i] * fe_0 + ts_xxyyz_xxxxxz[i] * pa_y[i];

        ts_xxyyyz_xxxxyy[i] = ts_xxyyy_xxxxyy[i] * pa_z[i];

        ts_xxyyyz_xxxxyz[i] = ts_xxyyy_xxxxy[i] * fe_0 + ts_xxyyy_xxxxyz[i] * pa_z[i];

        ts_xxyyyz_xxxxzz[i] = 2.0 * ts_xxyz_xxxxzz[i] * fe_0 + ts_xxyyz_xxxxzz[i] * pa_y[i];

        ts_xxyyyz_xxxyyy[i] = ts_xxyyy_xxxyyy[i] * pa_z[i];

        ts_xxyyyz_xxxyyz[i] = ts_xxyyy_xxxyy[i] * fe_0 + ts_xxyyy_xxxyyz[i] * pa_z[i];

        ts_xxyyyz_xxxyzz[i] = 2.0 * ts_xxyyy_xxxyz[i] * fe_0 + ts_xxyyy_xxxyzz[i] * pa_z[i];

        ts_xxyyyz_xxxzzz[i] = 2.0 * ts_xxyz_xxxzzz[i] * fe_0 + ts_xxyyz_xxxzzz[i] * pa_y[i];

        ts_xxyyyz_xxyyyy[i] = ts_xxyyy_xxyyyy[i] * pa_z[i];

        ts_xxyyyz_xxyyyz[i] = ts_xxyyy_xxyyy[i] * fe_0 + ts_xxyyy_xxyyyz[i] * pa_z[i];

        ts_xxyyyz_xxyyzz[i] = 2.0 * ts_xxyyy_xxyyz[i] * fe_0 + ts_xxyyy_xxyyzz[i] * pa_z[i];

        ts_xxyyyz_xxyzzz[i] = 3.0 * ts_xxyyy_xxyzz[i] * fe_0 + ts_xxyyy_xxyzzz[i] * pa_z[i];

        ts_xxyyyz_xxzzzz[i] = 2.0 * ts_xxyz_xxzzzz[i] * fe_0 + ts_xxyyz_xxzzzz[i] * pa_y[i];

        ts_xxyyyz_xyyyyy[i] = ts_xxyyy_xyyyyy[i] * pa_z[i];

        ts_xxyyyz_xyyyyz[i] = ts_xxyyy_xyyyy[i] * fe_0 + ts_xxyyy_xyyyyz[i] * pa_z[i];

        ts_xxyyyz_xyyyzz[i] = 2.0 * ts_xxyyy_xyyyz[i] * fe_0 + ts_xxyyy_xyyyzz[i] * pa_z[i];

        ts_xxyyyz_xyyzzz[i] = 3.0 * ts_xxyyy_xyyzz[i] * fe_0 + ts_xxyyy_xyyzzz[i] * pa_z[i];

        ts_xxyyyz_xyzzzz[i] = 4.0 * ts_xxyyy_xyzzz[i] * fe_0 + ts_xxyyy_xyzzzz[i] * pa_z[i];

        ts_xxyyyz_xzzzzz[i] = 2.0 * ts_xxyz_xzzzzz[i] * fe_0 + ts_xxyyz_xzzzzz[i] * pa_y[i];

        ts_xxyyyz_yyyyyy[i] = ts_xxyyy_yyyyyy[i] * pa_z[i];

        ts_xxyyyz_yyyyyz[i] = ts_yyyz_yyyyyz[i] * fe_0 + ts_xyyyz_yyyyyz[i] * pa_x[i];

        ts_xxyyyz_yyyyzz[i] = ts_yyyz_yyyyzz[i] * fe_0 + ts_xyyyz_yyyyzz[i] * pa_x[i];

        ts_xxyyyz_yyyzzz[i] = ts_yyyz_yyyzzz[i] * fe_0 + ts_xyyyz_yyyzzz[i] * pa_x[i];

        ts_xxyyyz_yyzzzz[i] = ts_yyyz_yyzzzz[i] * fe_0 + ts_xyyyz_yyzzzz[i] * pa_x[i];

        ts_xxyyyz_yzzzzz[i] = ts_yyyz_yzzzzz[i] * fe_0 + ts_xyyyz_yzzzzz[i] * pa_x[i];

        ts_xxyyyz_zzzzzz[i] = ts_yyyz_zzzzzz[i] * fe_0 + ts_xyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 336-364 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xxxxxy, ts_xxyy_xxxxyy, ts_xxyy_xxxyyy, ts_xxyy_xxyyyy, ts_xxyy_xyyyyy, ts_xxyyz_xxxxxy, ts_xxyyz_xxxxyy, ts_xxyyz_xxxyyy, ts_xxyyz_xxyyyy, ts_xxyyz_xyyyyy, ts_xxyyzz_xxxxxx, ts_xxyyzz_xxxxxy, ts_xxyyzz_xxxxxz, ts_xxyyzz_xxxxyy, ts_xxyyzz_xxxxyz, ts_xxyyzz_xxxxzz, ts_xxyyzz_xxxyyy, ts_xxyyzz_xxxyyz, ts_xxyyzz_xxxyzz, ts_xxyyzz_xxxzzz, ts_xxyyzz_xxyyyy, ts_xxyyzz_xxyyyz, ts_xxyyzz_xxyyzz, ts_xxyyzz_xxyzzz, ts_xxyyzz_xxzzzz, ts_xxyyzz_xyyyyy, ts_xxyyzz_xyyyyz, ts_xxyyzz_xyyyzz, ts_xxyyzz_xyyzzz, ts_xxyyzz_xyzzzz, ts_xxyyzz_xzzzzz, ts_xxyyzz_yyyyyy, ts_xxyyzz_yyyyyz, ts_xxyyzz_yyyyzz, ts_xxyyzz_yyyzzz, ts_xxyyzz_yyzzzz, ts_xxyyzz_yzzzzz, ts_xxyyzz_zzzzzz, ts_xxyzz_xxxxxx, ts_xxyzz_xxxxxz, ts_xxyzz_xxxxzz, ts_xxyzz_xxxzzz, ts_xxyzz_xxzzzz, ts_xxyzz_xzzzzz, ts_xxzz_xxxxxx, ts_xxzz_xxxxxz, ts_xxzz_xxxxzz, ts_xxzz_xxxzzz, ts_xxzz_xxzzzz, ts_xxzz_xzzzzz, ts_xyyzz_xxxxyz, ts_xyyzz_xxxyyz, ts_xyyzz_xxxyz, ts_xyyzz_xxxyzz, ts_xyyzz_xxyyyz, ts_xyyzz_xxyyz, ts_xyyzz_xxyyzz, ts_xyyzz_xxyzz, ts_xyyzz_xxyzzz, ts_xyyzz_xyyyyz, ts_xyyzz_xyyyz, ts_xyyzz_xyyyzz, ts_xyyzz_xyyzz, ts_xyyzz_xyyzzz, ts_xyyzz_xyzzz, ts_xyyzz_xyzzzz, ts_xyyzz_yyyyyy, ts_xyyzz_yyyyyz, ts_xyyzz_yyyyz, ts_xyyzz_yyyyzz, ts_xyyzz_yyyzz, ts_xyyzz_yyyzzz, ts_xyyzz_yyzzz, ts_xyyzz_yyzzzz, ts_xyyzz_yzzzz, ts_xyyzz_yzzzzz, ts_xyyzz_zzzzzz, ts_yyzz_xxxxyz, ts_yyzz_xxxyyz, ts_yyzz_xxxyzz, ts_yyzz_xxyyyz, ts_yyzz_xxyyzz, ts_yyzz_xxyzzz, ts_yyzz_xyyyyz, ts_yyzz_xyyyzz, ts_yyzz_xyyzzz, ts_yyzz_xyzzzz, ts_yyzz_yyyyyy, ts_yyzz_yyyyyz, ts_yyzz_yyyyzz, ts_yyzz_yyyzzz, ts_yyzz_yyzzzz, ts_yyzz_yzzzzz, ts_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyzz_xxxxxx[i] = ts_xxzz_xxxxxx[i] * fe_0 + ts_xxyzz_xxxxxx[i] * pa_y[i];

        ts_xxyyzz_xxxxxy[i] = ts_xxyy_xxxxxy[i] * fe_0 + ts_xxyyz_xxxxxy[i] * pa_z[i];

        ts_xxyyzz_xxxxxz[i] = ts_xxzz_xxxxxz[i] * fe_0 + ts_xxyzz_xxxxxz[i] * pa_y[i];

        ts_xxyyzz_xxxxyy[i] = ts_xxyy_xxxxyy[i] * fe_0 + ts_xxyyz_xxxxyy[i] * pa_z[i];

        ts_xxyyzz_xxxxyz[i] = ts_yyzz_xxxxyz[i] * fe_0 + 4.0 * ts_xyyzz_xxxyz[i] * fe_0 + ts_xyyzz_xxxxyz[i] * pa_x[i];

        ts_xxyyzz_xxxxzz[i] = ts_xxzz_xxxxzz[i] * fe_0 + ts_xxyzz_xxxxzz[i] * pa_y[i];

        ts_xxyyzz_xxxyyy[i] = ts_xxyy_xxxyyy[i] * fe_0 + ts_xxyyz_xxxyyy[i] * pa_z[i];

        ts_xxyyzz_xxxyyz[i] = ts_yyzz_xxxyyz[i] * fe_0 + 3.0 * ts_xyyzz_xxyyz[i] * fe_0 + ts_xyyzz_xxxyyz[i] * pa_x[i];

        ts_xxyyzz_xxxyzz[i] = ts_yyzz_xxxyzz[i] * fe_0 + 3.0 * ts_xyyzz_xxyzz[i] * fe_0 + ts_xyyzz_xxxyzz[i] * pa_x[i];

        ts_xxyyzz_xxxzzz[i] = ts_xxzz_xxxzzz[i] * fe_0 + ts_xxyzz_xxxzzz[i] * pa_y[i];

        ts_xxyyzz_xxyyyy[i] = ts_xxyy_xxyyyy[i] * fe_0 + ts_xxyyz_xxyyyy[i] * pa_z[i];

        ts_xxyyzz_xxyyyz[i] = ts_yyzz_xxyyyz[i] * fe_0 + 2.0 * ts_xyyzz_xyyyz[i] * fe_0 + ts_xyyzz_xxyyyz[i] * pa_x[i];

        ts_xxyyzz_xxyyzz[i] = ts_yyzz_xxyyzz[i] * fe_0 + 2.0 * ts_xyyzz_xyyzz[i] * fe_0 + ts_xyyzz_xxyyzz[i] * pa_x[i];

        ts_xxyyzz_xxyzzz[i] = ts_yyzz_xxyzzz[i] * fe_0 + 2.0 * ts_xyyzz_xyzzz[i] * fe_0 + ts_xyyzz_xxyzzz[i] * pa_x[i];

        ts_xxyyzz_xxzzzz[i] = ts_xxzz_xxzzzz[i] * fe_0 + ts_xxyzz_xxzzzz[i] * pa_y[i];

        ts_xxyyzz_xyyyyy[i] = ts_xxyy_xyyyyy[i] * fe_0 + ts_xxyyz_xyyyyy[i] * pa_z[i];

        ts_xxyyzz_xyyyyz[i] = ts_yyzz_xyyyyz[i] * fe_0 + ts_xyyzz_yyyyz[i] * fe_0 + ts_xyyzz_xyyyyz[i] * pa_x[i];

        ts_xxyyzz_xyyyzz[i] = ts_yyzz_xyyyzz[i] * fe_0 + ts_xyyzz_yyyzz[i] * fe_0 + ts_xyyzz_xyyyzz[i] * pa_x[i];

        ts_xxyyzz_xyyzzz[i] = ts_yyzz_xyyzzz[i] * fe_0 + ts_xyyzz_yyzzz[i] * fe_0 + ts_xyyzz_xyyzzz[i] * pa_x[i];

        ts_xxyyzz_xyzzzz[i] = ts_yyzz_xyzzzz[i] * fe_0 + ts_xyyzz_yzzzz[i] * fe_0 + ts_xyyzz_xyzzzz[i] * pa_x[i];

        ts_xxyyzz_xzzzzz[i] = ts_xxzz_xzzzzz[i] * fe_0 + ts_xxyzz_xzzzzz[i] * pa_y[i];

        ts_xxyyzz_yyyyyy[i] = ts_yyzz_yyyyyy[i] * fe_0 + ts_xyyzz_yyyyyy[i] * pa_x[i];

        ts_xxyyzz_yyyyyz[i] = ts_yyzz_yyyyyz[i] * fe_0 + ts_xyyzz_yyyyyz[i] * pa_x[i];

        ts_xxyyzz_yyyyzz[i] = ts_yyzz_yyyyzz[i] * fe_0 + ts_xyyzz_yyyyzz[i] * pa_x[i];

        ts_xxyyzz_yyyzzz[i] = ts_yyzz_yyyzzz[i] * fe_0 + ts_xyyzz_yyyzzz[i] * pa_x[i];

        ts_xxyyzz_yyzzzz[i] = ts_yyzz_yyzzzz[i] * fe_0 + ts_xyyzz_yyzzzz[i] * pa_x[i];

        ts_xxyyzz_yzzzzz[i] = ts_yyzz_yzzzzz[i] * fe_0 + ts_xyyzz_yzzzzz[i] * pa_x[i];

        ts_xxyyzz_zzzzzz[i] = ts_yyzz_zzzzzz[i] * fe_0 + ts_xyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzzz_xxxxxx, ts_xxyzzz_xxxxxy, ts_xxyzzz_xxxxxz, ts_xxyzzz_xxxxyy, ts_xxyzzz_xxxxyz, ts_xxyzzz_xxxxzz, ts_xxyzzz_xxxyyy, ts_xxyzzz_xxxyyz, ts_xxyzzz_xxxyzz, ts_xxyzzz_xxxzzz, ts_xxyzzz_xxyyyy, ts_xxyzzz_xxyyyz, ts_xxyzzz_xxyyzz, ts_xxyzzz_xxyzzz, ts_xxyzzz_xxzzzz, ts_xxyzzz_xyyyyy, ts_xxyzzz_xyyyyz, ts_xxyzzz_xyyyzz, ts_xxyzzz_xyyzzz, ts_xxyzzz_xyzzzz, ts_xxyzzz_xzzzzz, ts_xxyzzz_yyyyyy, ts_xxyzzz_yyyyyz, ts_xxyzzz_yyyyzz, ts_xxyzzz_yyyzzz, ts_xxyzzz_yyzzzz, ts_xxyzzz_yzzzzz, ts_xxyzzz_zzzzzz, ts_xxzzz_xxxxx, ts_xxzzz_xxxxxx, ts_xxzzz_xxxxxy, ts_xxzzz_xxxxxz, ts_xxzzz_xxxxy, ts_xxzzz_xxxxyy, ts_xxzzz_xxxxyz, ts_xxzzz_xxxxz, ts_xxzzz_xxxxzz, ts_xxzzz_xxxyy, ts_xxzzz_xxxyyy, ts_xxzzz_xxxyyz, ts_xxzzz_xxxyz, ts_xxzzz_xxxyzz, ts_xxzzz_xxxzz, ts_xxzzz_xxxzzz, ts_xxzzz_xxyyy, ts_xxzzz_xxyyyy, ts_xxzzz_xxyyyz, ts_xxzzz_xxyyz, ts_xxzzz_xxyyzz, ts_xxzzz_xxyzz, ts_xxzzz_xxyzzz, ts_xxzzz_xxzzz, ts_xxzzz_xxzzzz, ts_xxzzz_xyyyy, ts_xxzzz_xyyyyy, ts_xxzzz_xyyyyz, ts_xxzzz_xyyyz, ts_xxzzz_xyyyzz, ts_xxzzz_xyyzz, ts_xxzzz_xyyzzz, ts_xxzzz_xyzzz, ts_xxzzz_xyzzzz, ts_xxzzz_xzzzz, ts_xxzzz_xzzzzz, ts_xxzzz_zzzzzz, ts_xyzzz_yyyyyy, ts_xyzzz_yyyyyz, ts_xyzzz_yyyyzz, ts_xyzzz_yyyzzz, ts_xyzzz_yyzzzz, ts_xyzzz_yzzzzz, ts_yzzz_yyyyyy, ts_yzzz_yyyyyz, ts_yzzz_yyyyzz, ts_yzzz_yyyzzz, ts_yzzz_yyzzzz, ts_yzzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzzz_xxxxxx[i] = ts_xxzzz_xxxxxx[i] * pa_y[i];

        ts_xxyzzz_xxxxxy[i] = ts_xxzzz_xxxxx[i] * fe_0 + ts_xxzzz_xxxxxy[i] * pa_y[i];

        ts_xxyzzz_xxxxxz[i] = ts_xxzzz_xxxxxz[i] * pa_y[i];

        ts_xxyzzz_xxxxyy[i] = 2.0 * ts_xxzzz_xxxxy[i] * fe_0 + ts_xxzzz_xxxxyy[i] * pa_y[i];

        ts_xxyzzz_xxxxyz[i] = ts_xxzzz_xxxxz[i] * fe_0 + ts_xxzzz_xxxxyz[i] * pa_y[i];

        ts_xxyzzz_xxxxzz[i] = ts_xxzzz_xxxxzz[i] * pa_y[i];

        ts_xxyzzz_xxxyyy[i] = 3.0 * ts_xxzzz_xxxyy[i] * fe_0 + ts_xxzzz_xxxyyy[i] * pa_y[i];

        ts_xxyzzz_xxxyyz[i] = 2.0 * ts_xxzzz_xxxyz[i] * fe_0 + ts_xxzzz_xxxyyz[i] * pa_y[i];

        ts_xxyzzz_xxxyzz[i] = ts_xxzzz_xxxzz[i] * fe_0 + ts_xxzzz_xxxyzz[i] * pa_y[i];

        ts_xxyzzz_xxxzzz[i] = ts_xxzzz_xxxzzz[i] * pa_y[i];

        ts_xxyzzz_xxyyyy[i] = 4.0 * ts_xxzzz_xxyyy[i] * fe_0 + ts_xxzzz_xxyyyy[i] * pa_y[i];

        ts_xxyzzz_xxyyyz[i] = 3.0 * ts_xxzzz_xxyyz[i] * fe_0 + ts_xxzzz_xxyyyz[i] * pa_y[i];

        ts_xxyzzz_xxyyzz[i] = 2.0 * ts_xxzzz_xxyzz[i] * fe_0 + ts_xxzzz_xxyyzz[i] * pa_y[i];

        ts_xxyzzz_xxyzzz[i] = ts_xxzzz_xxzzz[i] * fe_0 + ts_xxzzz_xxyzzz[i] * pa_y[i];

        ts_xxyzzz_xxzzzz[i] = ts_xxzzz_xxzzzz[i] * pa_y[i];

        ts_xxyzzz_xyyyyy[i] = 5.0 * ts_xxzzz_xyyyy[i] * fe_0 + ts_xxzzz_xyyyyy[i] * pa_y[i];

        ts_xxyzzz_xyyyyz[i] = 4.0 * ts_xxzzz_xyyyz[i] * fe_0 + ts_xxzzz_xyyyyz[i] * pa_y[i];

        ts_xxyzzz_xyyyzz[i] = 3.0 * ts_xxzzz_xyyzz[i] * fe_0 + ts_xxzzz_xyyyzz[i] * pa_y[i];

        ts_xxyzzz_xyyzzz[i] = 2.0 * ts_xxzzz_xyzzz[i] * fe_0 + ts_xxzzz_xyyzzz[i] * pa_y[i];

        ts_xxyzzz_xyzzzz[i] = ts_xxzzz_xzzzz[i] * fe_0 + ts_xxzzz_xyzzzz[i] * pa_y[i];

        ts_xxyzzz_xzzzzz[i] = ts_xxzzz_xzzzzz[i] * pa_y[i];

        ts_xxyzzz_yyyyyy[i] = ts_yzzz_yyyyyy[i] * fe_0 + ts_xyzzz_yyyyyy[i] * pa_x[i];

        ts_xxyzzz_yyyyyz[i] = ts_yzzz_yyyyyz[i] * fe_0 + ts_xyzzz_yyyyyz[i] * pa_x[i];

        ts_xxyzzz_yyyyzz[i] = ts_yzzz_yyyyzz[i] * fe_0 + ts_xyzzz_yyyyzz[i] * pa_x[i];

        ts_xxyzzz_yyyzzz[i] = ts_yzzz_yyyzzz[i] * fe_0 + ts_xyzzz_yyyzzz[i] * pa_x[i];

        ts_xxyzzz_yyzzzz[i] = ts_yzzz_yyzzzz[i] * fe_0 + ts_xyzzz_yyzzzz[i] * pa_x[i];

        ts_xxyzzz_yzzzzz[i] = ts_yzzz_yzzzzz[i] * fe_0 + ts_xyzzz_yzzzzz[i] * pa_x[i];

        ts_xxyzzz_zzzzzz[i] = ts_xxzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_z, ts_xxzz_xxxxxx, ts_xxzz_xxxxxy, ts_xxzz_xxxxyy, ts_xxzz_xxxyyy, ts_xxzz_xxyyyy, ts_xxzz_xyyyyy, ts_xxzzz_xxxxxx, ts_xxzzz_xxxxxy, ts_xxzzz_xxxxyy, ts_xxzzz_xxxyyy, ts_xxzzz_xxyyyy, ts_xxzzz_xyyyyy, ts_xxzzzz_xxxxxx, ts_xxzzzz_xxxxxy, ts_xxzzzz_xxxxxz, ts_xxzzzz_xxxxyy, ts_xxzzzz_xxxxyz, ts_xxzzzz_xxxxzz, ts_xxzzzz_xxxyyy, ts_xxzzzz_xxxyyz, ts_xxzzzz_xxxyzz, ts_xxzzzz_xxxzzz, ts_xxzzzz_xxyyyy, ts_xxzzzz_xxyyyz, ts_xxzzzz_xxyyzz, ts_xxzzzz_xxyzzz, ts_xxzzzz_xxzzzz, ts_xxzzzz_xyyyyy, ts_xxzzzz_xyyyyz, ts_xxzzzz_xyyyzz, ts_xxzzzz_xyyzzz, ts_xxzzzz_xyzzzz, ts_xxzzzz_xzzzzz, ts_xxzzzz_yyyyyy, ts_xxzzzz_yyyyyz, ts_xxzzzz_yyyyzz, ts_xxzzzz_yyyzzz, ts_xxzzzz_yyzzzz, ts_xxzzzz_yzzzzz, ts_xxzzzz_zzzzzz, ts_xzzzz_xxxxxz, ts_xzzzz_xxxxyz, ts_xzzzz_xxxxz, ts_xzzzz_xxxxzz, ts_xzzzz_xxxyyz, ts_xzzzz_xxxyz, ts_xzzzz_xxxyzz, ts_xzzzz_xxxzz, ts_xzzzz_xxxzzz, ts_xzzzz_xxyyyz, ts_xzzzz_xxyyz, ts_xzzzz_xxyyzz, ts_xzzzz_xxyzz, ts_xzzzz_xxyzzz, ts_xzzzz_xxzzz, ts_xzzzz_xxzzzz, ts_xzzzz_xyyyyz, ts_xzzzz_xyyyz, ts_xzzzz_xyyyzz, ts_xzzzz_xyyzz, ts_xzzzz_xyyzzz, ts_xzzzz_xyzzz, ts_xzzzz_xyzzzz, ts_xzzzz_xzzzz, ts_xzzzz_xzzzzz, ts_xzzzz_yyyyyy, ts_xzzzz_yyyyyz, ts_xzzzz_yyyyz, ts_xzzzz_yyyyzz, ts_xzzzz_yyyzz, ts_xzzzz_yyyzzz, ts_xzzzz_yyzzz, ts_xzzzz_yyzzzz, ts_xzzzz_yzzzz, ts_xzzzz_yzzzzz, ts_xzzzz_zzzzz, ts_xzzzz_zzzzzz, ts_zzzz_xxxxxz, ts_zzzz_xxxxyz, ts_zzzz_xxxxzz, ts_zzzz_xxxyyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyyy, ts_zzzz_yyyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzzz_xxxxxx[i] = 3.0 * ts_xxzz_xxxxxx[i] * fe_0 + ts_xxzzz_xxxxxx[i] * pa_z[i];

        ts_xxzzzz_xxxxxy[i] = 3.0 * ts_xxzz_xxxxxy[i] * fe_0 + ts_xxzzz_xxxxxy[i] * pa_z[i];

        ts_xxzzzz_xxxxxz[i] = ts_zzzz_xxxxxz[i] * fe_0 + 5.0 * ts_xzzzz_xxxxz[i] * fe_0 + ts_xzzzz_xxxxxz[i] * pa_x[i];

        ts_xxzzzz_xxxxyy[i] = 3.0 * ts_xxzz_xxxxyy[i] * fe_0 + ts_xxzzz_xxxxyy[i] * pa_z[i];

        ts_xxzzzz_xxxxyz[i] = ts_zzzz_xxxxyz[i] * fe_0 + 4.0 * ts_xzzzz_xxxyz[i] * fe_0 + ts_xzzzz_xxxxyz[i] * pa_x[i];

        ts_xxzzzz_xxxxzz[i] = ts_zzzz_xxxxzz[i] * fe_0 + 4.0 * ts_xzzzz_xxxzz[i] * fe_0 + ts_xzzzz_xxxxzz[i] * pa_x[i];

        ts_xxzzzz_xxxyyy[i] = 3.0 * ts_xxzz_xxxyyy[i] * fe_0 + ts_xxzzz_xxxyyy[i] * pa_z[i];

        ts_xxzzzz_xxxyyz[i] = ts_zzzz_xxxyyz[i] * fe_0 + 3.0 * ts_xzzzz_xxyyz[i] * fe_0 + ts_xzzzz_xxxyyz[i] * pa_x[i];

        ts_xxzzzz_xxxyzz[i] = ts_zzzz_xxxyzz[i] * fe_0 + 3.0 * ts_xzzzz_xxyzz[i] * fe_0 + ts_xzzzz_xxxyzz[i] * pa_x[i];

        ts_xxzzzz_xxxzzz[i] = ts_zzzz_xxxzzz[i] * fe_0 + 3.0 * ts_xzzzz_xxzzz[i] * fe_0 + ts_xzzzz_xxxzzz[i] * pa_x[i];

        ts_xxzzzz_xxyyyy[i] = 3.0 * ts_xxzz_xxyyyy[i] * fe_0 + ts_xxzzz_xxyyyy[i] * pa_z[i];

        ts_xxzzzz_xxyyyz[i] = ts_zzzz_xxyyyz[i] * fe_0 + 2.0 * ts_xzzzz_xyyyz[i] * fe_0 + ts_xzzzz_xxyyyz[i] * pa_x[i];

        ts_xxzzzz_xxyyzz[i] = ts_zzzz_xxyyzz[i] * fe_0 + 2.0 * ts_xzzzz_xyyzz[i] * fe_0 + ts_xzzzz_xxyyzz[i] * pa_x[i];

        ts_xxzzzz_xxyzzz[i] = ts_zzzz_xxyzzz[i] * fe_0 + 2.0 * ts_xzzzz_xyzzz[i] * fe_0 + ts_xzzzz_xxyzzz[i] * pa_x[i];

        ts_xxzzzz_xxzzzz[i] = ts_zzzz_xxzzzz[i] * fe_0 + 2.0 * ts_xzzzz_xzzzz[i] * fe_0 + ts_xzzzz_xxzzzz[i] * pa_x[i];

        ts_xxzzzz_xyyyyy[i] = 3.0 * ts_xxzz_xyyyyy[i] * fe_0 + ts_xxzzz_xyyyyy[i] * pa_z[i];

        ts_xxzzzz_xyyyyz[i] = ts_zzzz_xyyyyz[i] * fe_0 + ts_xzzzz_yyyyz[i] * fe_0 + ts_xzzzz_xyyyyz[i] * pa_x[i];

        ts_xxzzzz_xyyyzz[i] = ts_zzzz_xyyyzz[i] * fe_0 + ts_xzzzz_yyyzz[i] * fe_0 + ts_xzzzz_xyyyzz[i] * pa_x[i];

        ts_xxzzzz_xyyzzz[i] = ts_zzzz_xyyzzz[i] * fe_0 + ts_xzzzz_yyzzz[i] * fe_0 + ts_xzzzz_xyyzzz[i] * pa_x[i];

        ts_xxzzzz_xyzzzz[i] = ts_zzzz_xyzzzz[i] * fe_0 + ts_xzzzz_yzzzz[i] * fe_0 + ts_xzzzz_xyzzzz[i] * pa_x[i];

        ts_xxzzzz_xzzzzz[i] = ts_zzzz_xzzzzz[i] * fe_0 + ts_xzzzz_zzzzz[i] * fe_0 + ts_xzzzz_xzzzzz[i] * pa_x[i];

        ts_xxzzzz_yyyyyy[i] = ts_zzzz_yyyyyy[i] * fe_0 + ts_xzzzz_yyyyyy[i] * pa_x[i];

        ts_xxzzzz_yyyyyz[i] = ts_zzzz_yyyyyz[i] * fe_0 + ts_xzzzz_yyyyyz[i] * pa_x[i];

        ts_xxzzzz_yyyyzz[i] = ts_zzzz_yyyyzz[i] * fe_0 + ts_xzzzz_yyyyzz[i] * pa_x[i];

        ts_xxzzzz_yyyzzz[i] = ts_zzzz_yyyzzz[i] * fe_0 + ts_xzzzz_yyyzzz[i] * pa_x[i];

        ts_xxzzzz_yyzzzz[i] = ts_zzzz_yyzzzz[i] * fe_0 + ts_xzzzz_yyzzzz[i] * pa_x[i];

        ts_xxzzzz_yzzzzz[i] = ts_zzzz_yzzzzz[i] * fe_0 + ts_xzzzz_yzzzzz[i] * pa_x[i];

        ts_xxzzzz_zzzzzz[i] = ts_zzzz_zzzzzz[i] * fe_0 + ts_xzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, ts_xyyyyy_xxxxxx, ts_xyyyyy_xxxxxy, ts_xyyyyy_xxxxxz, ts_xyyyyy_xxxxyy, ts_xyyyyy_xxxxyz, ts_xyyyyy_xxxxzz, ts_xyyyyy_xxxyyy, ts_xyyyyy_xxxyyz, ts_xyyyyy_xxxyzz, ts_xyyyyy_xxxzzz, ts_xyyyyy_xxyyyy, ts_xyyyyy_xxyyyz, ts_xyyyyy_xxyyzz, ts_xyyyyy_xxyzzz, ts_xyyyyy_xxzzzz, ts_xyyyyy_xyyyyy, ts_xyyyyy_xyyyyz, ts_xyyyyy_xyyyzz, ts_xyyyyy_xyyzzz, ts_xyyyyy_xyzzzz, ts_xyyyyy_xzzzzz, ts_xyyyyy_yyyyyy, ts_xyyyyy_yyyyyz, ts_xyyyyy_yyyyzz, ts_xyyyyy_yyyzzz, ts_xyyyyy_yyzzzz, ts_xyyyyy_yzzzzz, ts_xyyyyy_zzzzzz, ts_yyyyy_xxxxx, ts_yyyyy_xxxxxx, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxxz, ts_yyyyy_xxxxy, ts_yyyyy_xxxxyy, ts_yyyyy_xxxxyz, ts_yyyyy_xxxxz, ts_yyyyy_xxxxzz, ts_yyyyy_xxxyy, ts_yyyyy_xxxyyy, ts_yyyyy_xxxyyz, ts_yyyyy_xxxyz, ts_yyyyy_xxxyzz, ts_yyyyy_xxxzz, ts_yyyyy_xxxzzz, ts_yyyyy_xxyyy, ts_yyyyy_xxyyyy, ts_yyyyy_xxyyyz, ts_yyyyy_xxyyz, ts_yyyyy_xxyyzz, ts_yyyyy_xxyzz, ts_yyyyy_xxyzzz, ts_yyyyy_xxzzz, ts_yyyyy_xxzzzz, ts_yyyyy_xyyyy, ts_yyyyy_xyyyyy, ts_yyyyy_xyyyyz, ts_yyyyy_xyyyz, ts_yyyyy_xyyyzz, ts_yyyyy_xyyzz, ts_yyyyy_xyyzzz, ts_yyyyy_xyzzz, ts_yyyyy_xyzzzz, ts_yyyyy_xzzzz, ts_yyyyy_xzzzzz, ts_yyyyy_yyyyy, ts_yyyyy_yyyyyy, ts_yyyyy_yyyyyz, ts_yyyyy_yyyyz, ts_yyyyy_yyyyzz, ts_yyyyy_yyyzz, ts_yyyyy_yyyzzz, ts_yyyyy_yyzzz, ts_yyyyy_yyzzzz, ts_yyyyy_yzzzz, ts_yyyyy_yzzzzz, ts_yyyyy_zzzzz, ts_yyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyy_xxxxxx[i] = 6.0 * ts_yyyyy_xxxxx[i] * fe_0 + ts_yyyyy_xxxxxx[i] * pa_x[i];

        ts_xyyyyy_xxxxxy[i] = 5.0 * ts_yyyyy_xxxxy[i] * fe_0 + ts_yyyyy_xxxxxy[i] * pa_x[i];

        ts_xyyyyy_xxxxxz[i] = 5.0 * ts_yyyyy_xxxxz[i] * fe_0 + ts_yyyyy_xxxxxz[i] * pa_x[i];

        ts_xyyyyy_xxxxyy[i] = 4.0 * ts_yyyyy_xxxyy[i] * fe_0 + ts_yyyyy_xxxxyy[i] * pa_x[i];

        ts_xyyyyy_xxxxyz[i] = 4.0 * ts_yyyyy_xxxyz[i] * fe_0 + ts_yyyyy_xxxxyz[i] * pa_x[i];

        ts_xyyyyy_xxxxzz[i] = 4.0 * ts_yyyyy_xxxzz[i] * fe_0 + ts_yyyyy_xxxxzz[i] * pa_x[i];

        ts_xyyyyy_xxxyyy[i] = 3.0 * ts_yyyyy_xxyyy[i] * fe_0 + ts_yyyyy_xxxyyy[i] * pa_x[i];

        ts_xyyyyy_xxxyyz[i] = 3.0 * ts_yyyyy_xxyyz[i] * fe_0 + ts_yyyyy_xxxyyz[i] * pa_x[i];

        ts_xyyyyy_xxxyzz[i] = 3.0 * ts_yyyyy_xxyzz[i] * fe_0 + ts_yyyyy_xxxyzz[i] * pa_x[i];

        ts_xyyyyy_xxxzzz[i] = 3.0 * ts_yyyyy_xxzzz[i] * fe_0 + ts_yyyyy_xxxzzz[i] * pa_x[i];

        ts_xyyyyy_xxyyyy[i] = 2.0 * ts_yyyyy_xyyyy[i] * fe_0 + ts_yyyyy_xxyyyy[i] * pa_x[i];

        ts_xyyyyy_xxyyyz[i] = 2.0 * ts_yyyyy_xyyyz[i] * fe_0 + ts_yyyyy_xxyyyz[i] * pa_x[i];

        ts_xyyyyy_xxyyzz[i] = 2.0 * ts_yyyyy_xyyzz[i] * fe_0 + ts_yyyyy_xxyyzz[i] * pa_x[i];

        ts_xyyyyy_xxyzzz[i] = 2.0 * ts_yyyyy_xyzzz[i] * fe_0 + ts_yyyyy_xxyzzz[i] * pa_x[i];

        ts_xyyyyy_xxzzzz[i] = 2.0 * ts_yyyyy_xzzzz[i] * fe_0 + ts_yyyyy_xxzzzz[i] * pa_x[i];

        ts_xyyyyy_xyyyyy[i] = ts_yyyyy_yyyyy[i] * fe_0 + ts_yyyyy_xyyyyy[i] * pa_x[i];

        ts_xyyyyy_xyyyyz[i] = ts_yyyyy_yyyyz[i] * fe_0 + ts_yyyyy_xyyyyz[i] * pa_x[i];

        ts_xyyyyy_xyyyzz[i] = ts_yyyyy_yyyzz[i] * fe_0 + ts_yyyyy_xyyyzz[i] * pa_x[i];

        ts_xyyyyy_xyyzzz[i] = ts_yyyyy_yyzzz[i] * fe_0 + ts_yyyyy_xyyzzz[i] * pa_x[i];

        ts_xyyyyy_xyzzzz[i] = ts_yyyyy_yzzzz[i] * fe_0 + ts_yyyyy_xyzzzz[i] * pa_x[i];

        ts_xyyyyy_xzzzzz[i] = ts_yyyyy_zzzzz[i] * fe_0 + ts_yyyyy_xzzzzz[i] * pa_x[i];

        ts_xyyyyy_yyyyyy[i] = ts_yyyyy_yyyyyy[i] * pa_x[i];

        ts_xyyyyy_yyyyyz[i] = ts_yyyyy_yyyyyz[i] * pa_x[i];

        ts_xyyyyy_yyyyzz[i] = ts_yyyyy_yyyyzz[i] * pa_x[i];

        ts_xyyyyy_yyyzzz[i] = ts_yyyyy_yyyzzz[i] * pa_x[i];

        ts_xyyyyy_yyzzzz[i] = ts_yyyyy_yyzzzz[i] * pa_x[i];

        ts_xyyyyy_yzzzzz[i] = ts_yyyyy_yzzzzz[i] * pa_x[i];

        ts_xyyyyy_zzzzzz[i] = ts_yyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 448-476 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyyy_xxxxxx, ts_xyyyy_xxxxxy, ts_xyyyy_xxxxyy, ts_xyyyy_xxxyyy, ts_xyyyy_xxyyyy, ts_xyyyy_xyyyyy, ts_xyyyyz_xxxxxx, ts_xyyyyz_xxxxxy, ts_xyyyyz_xxxxxz, ts_xyyyyz_xxxxyy, ts_xyyyyz_xxxxyz, ts_xyyyyz_xxxxzz, ts_xyyyyz_xxxyyy, ts_xyyyyz_xxxyyz, ts_xyyyyz_xxxyzz, ts_xyyyyz_xxxzzz, ts_xyyyyz_xxyyyy, ts_xyyyyz_xxyyyz, ts_xyyyyz_xxyyzz, ts_xyyyyz_xxyzzz, ts_xyyyyz_xxzzzz, ts_xyyyyz_xyyyyy, ts_xyyyyz_xyyyyz, ts_xyyyyz_xyyyzz, ts_xyyyyz_xyyzzz, ts_xyyyyz_xyzzzz, ts_xyyyyz_xzzzzz, ts_xyyyyz_yyyyyy, ts_xyyyyz_yyyyyz, ts_xyyyyz_yyyyzz, ts_xyyyyz_yyyzzz, ts_xyyyyz_yyzzzz, ts_xyyyyz_yzzzzz, ts_xyyyyz_zzzzzz, ts_yyyyz_xxxxxz, ts_yyyyz_xxxxyz, ts_yyyyz_xxxxz, ts_yyyyz_xxxxzz, ts_yyyyz_xxxyyz, ts_yyyyz_xxxyz, ts_yyyyz_xxxyzz, ts_yyyyz_xxxzz, ts_yyyyz_xxxzzz, ts_yyyyz_xxyyyz, ts_yyyyz_xxyyz, ts_yyyyz_xxyyzz, ts_yyyyz_xxyzz, ts_yyyyz_xxyzzz, ts_yyyyz_xxzzz, ts_yyyyz_xxzzzz, ts_yyyyz_xyyyyz, ts_yyyyz_xyyyz, ts_yyyyz_xyyyzz, ts_yyyyz_xyyzz, ts_yyyyz_xyyzzz, ts_yyyyz_xyzzz, ts_yyyyz_xyzzzz, ts_yyyyz_xzzzz, ts_yyyyz_xzzzzz, ts_yyyyz_yyyyyy, ts_yyyyz_yyyyyz, ts_yyyyz_yyyyz, ts_yyyyz_yyyyzz, ts_yyyyz_yyyzz, ts_yyyyz_yyyzzz, ts_yyyyz_yyzzz, ts_yyyyz_yyzzzz, ts_yyyyz_yzzzz, ts_yyyyz_yzzzzz, ts_yyyyz_zzzzz, ts_yyyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyz_xxxxxx[i] = ts_xyyyy_xxxxxx[i] * pa_z[i];

        ts_xyyyyz_xxxxxy[i] = ts_xyyyy_xxxxxy[i] * pa_z[i];

        ts_xyyyyz_xxxxxz[i] = 5.0 * ts_yyyyz_xxxxz[i] * fe_0 + ts_yyyyz_xxxxxz[i] * pa_x[i];

        ts_xyyyyz_xxxxyy[i] = ts_xyyyy_xxxxyy[i] * pa_z[i];

        ts_xyyyyz_xxxxyz[i] = 4.0 * ts_yyyyz_xxxyz[i] * fe_0 + ts_yyyyz_xxxxyz[i] * pa_x[i];

        ts_xyyyyz_xxxxzz[i] = 4.0 * ts_yyyyz_xxxzz[i] * fe_0 + ts_yyyyz_xxxxzz[i] * pa_x[i];

        ts_xyyyyz_xxxyyy[i] = ts_xyyyy_xxxyyy[i] * pa_z[i];

        ts_xyyyyz_xxxyyz[i] = 3.0 * ts_yyyyz_xxyyz[i] * fe_0 + ts_yyyyz_xxxyyz[i] * pa_x[i];

        ts_xyyyyz_xxxyzz[i] = 3.0 * ts_yyyyz_xxyzz[i] * fe_0 + ts_yyyyz_xxxyzz[i] * pa_x[i];

        ts_xyyyyz_xxxzzz[i] = 3.0 * ts_yyyyz_xxzzz[i] * fe_0 + ts_yyyyz_xxxzzz[i] * pa_x[i];

        ts_xyyyyz_xxyyyy[i] = ts_xyyyy_xxyyyy[i] * pa_z[i];

        ts_xyyyyz_xxyyyz[i] = 2.0 * ts_yyyyz_xyyyz[i] * fe_0 + ts_yyyyz_xxyyyz[i] * pa_x[i];

        ts_xyyyyz_xxyyzz[i] = 2.0 * ts_yyyyz_xyyzz[i] * fe_0 + ts_yyyyz_xxyyzz[i] * pa_x[i];

        ts_xyyyyz_xxyzzz[i] = 2.0 * ts_yyyyz_xyzzz[i] * fe_0 + ts_yyyyz_xxyzzz[i] * pa_x[i];

        ts_xyyyyz_xxzzzz[i] = 2.0 * ts_yyyyz_xzzzz[i] * fe_0 + ts_yyyyz_xxzzzz[i] * pa_x[i];

        ts_xyyyyz_xyyyyy[i] = ts_xyyyy_xyyyyy[i] * pa_z[i];

        ts_xyyyyz_xyyyyz[i] = ts_yyyyz_yyyyz[i] * fe_0 + ts_yyyyz_xyyyyz[i] * pa_x[i];

        ts_xyyyyz_xyyyzz[i] = ts_yyyyz_yyyzz[i] * fe_0 + ts_yyyyz_xyyyzz[i] * pa_x[i];

        ts_xyyyyz_xyyzzz[i] = ts_yyyyz_yyzzz[i] * fe_0 + ts_yyyyz_xyyzzz[i] * pa_x[i];

        ts_xyyyyz_xyzzzz[i] = ts_yyyyz_yzzzz[i] * fe_0 + ts_yyyyz_xyzzzz[i] * pa_x[i];

        ts_xyyyyz_xzzzzz[i] = ts_yyyyz_zzzzz[i] * fe_0 + ts_yyyyz_xzzzzz[i] * pa_x[i];

        ts_xyyyyz_yyyyyy[i] = ts_yyyyz_yyyyyy[i] * pa_x[i];

        ts_xyyyyz_yyyyyz[i] = ts_yyyyz_yyyyyz[i] * pa_x[i];

        ts_xyyyyz_yyyyzz[i] = ts_yyyyz_yyyyzz[i] * pa_x[i];

        ts_xyyyyz_yyyzzz[i] = ts_yyyyz_yyyzzz[i] * pa_x[i];

        ts_xyyyyz_yyzzzz[i] = ts_yyyyz_yyzzzz[i] * pa_x[i];

        ts_xyyyyz_yzzzzz[i] = ts_yyyyz_yzzzzz[i] * pa_x[i];

        ts_xyyyyz_zzzzzz[i] = ts_yyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 476-504 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, ts_xyyyzz_xxxxxx, ts_xyyyzz_xxxxxy, ts_xyyyzz_xxxxxz, ts_xyyyzz_xxxxyy, ts_xyyyzz_xxxxyz, ts_xyyyzz_xxxxzz, ts_xyyyzz_xxxyyy, ts_xyyyzz_xxxyyz, ts_xyyyzz_xxxyzz, ts_xyyyzz_xxxzzz, ts_xyyyzz_xxyyyy, ts_xyyyzz_xxyyyz, ts_xyyyzz_xxyyzz, ts_xyyyzz_xxyzzz, ts_xyyyzz_xxzzzz, ts_xyyyzz_xyyyyy, ts_xyyyzz_xyyyyz, ts_xyyyzz_xyyyzz, ts_xyyyzz_xyyzzz, ts_xyyyzz_xyzzzz, ts_xyyyzz_xzzzzz, ts_xyyyzz_yyyyyy, ts_xyyyzz_yyyyyz, ts_xyyyzz_yyyyzz, ts_xyyyzz_yyyzzz, ts_xyyyzz_yyzzzz, ts_xyyyzz_yzzzzz, ts_xyyyzz_zzzzzz, ts_yyyzz_xxxxx, ts_yyyzz_xxxxxx, ts_yyyzz_xxxxxy, ts_yyyzz_xxxxxz, ts_yyyzz_xxxxy, ts_yyyzz_xxxxyy, ts_yyyzz_xxxxyz, ts_yyyzz_xxxxz, ts_yyyzz_xxxxzz, ts_yyyzz_xxxyy, ts_yyyzz_xxxyyy, ts_yyyzz_xxxyyz, ts_yyyzz_xxxyz, ts_yyyzz_xxxyzz, ts_yyyzz_xxxzz, ts_yyyzz_xxxzzz, ts_yyyzz_xxyyy, ts_yyyzz_xxyyyy, ts_yyyzz_xxyyyz, ts_yyyzz_xxyyz, ts_yyyzz_xxyyzz, ts_yyyzz_xxyzz, ts_yyyzz_xxyzzz, ts_yyyzz_xxzzz, ts_yyyzz_xxzzzz, ts_yyyzz_xyyyy, ts_yyyzz_xyyyyy, ts_yyyzz_xyyyyz, ts_yyyzz_xyyyz, ts_yyyzz_xyyyzz, ts_yyyzz_xyyzz, ts_yyyzz_xyyzzz, ts_yyyzz_xyzzz, ts_yyyzz_xyzzzz, ts_yyyzz_xzzzz, ts_yyyzz_xzzzzz, ts_yyyzz_yyyyy, ts_yyyzz_yyyyyy, ts_yyyzz_yyyyyz, ts_yyyzz_yyyyz, ts_yyyzz_yyyyzz, ts_yyyzz_yyyzz, ts_yyyzz_yyyzzz, ts_yyyzz_yyzzz, ts_yyyzz_yyzzzz, ts_yyyzz_yzzzz, ts_yyyzz_yzzzzz, ts_yyyzz_zzzzz, ts_yyyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyzz_xxxxxx[i] = 6.0 * ts_yyyzz_xxxxx[i] * fe_0 + ts_yyyzz_xxxxxx[i] * pa_x[i];

        ts_xyyyzz_xxxxxy[i] = 5.0 * ts_yyyzz_xxxxy[i] * fe_0 + ts_yyyzz_xxxxxy[i] * pa_x[i];

        ts_xyyyzz_xxxxxz[i] = 5.0 * ts_yyyzz_xxxxz[i] * fe_0 + ts_yyyzz_xxxxxz[i] * pa_x[i];

        ts_xyyyzz_xxxxyy[i] = 4.0 * ts_yyyzz_xxxyy[i] * fe_0 + ts_yyyzz_xxxxyy[i] * pa_x[i];

        ts_xyyyzz_xxxxyz[i] = 4.0 * ts_yyyzz_xxxyz[i] * fe_0 + ts_yyyzz_xxxxyz[i] * pa_x[i];

        ts_xyyyzz_xxxxzz[i] = 4.0 * ts_yyyzz_xxxzz[i] * fe_0 + ts_yyyzz_xxxxzz[i] * pa_x[i];

        ts_xyyyzz_xxxyyy[i] = 3.0 * ts_yyyzz_xxyyy[i] * fe_0 + ts_yyyzz_xxxyyy[i] * pa_x[i];

        ts_xyyyzz_xxxyyz[i] = 3.0 * ts_yyyzz_xxyyz[i] * fe_0 + ts_yyyzz_xxxyyz[i] * pa_x[i];

        ts_xyyyzz_xxxyzz[i] = 3.0 * ts_yyyzz_xxyzz[i] * fe_0 + ts_yyyzz_xxxyzz[i] * pa_x[i];

        ts_xyyyzz_xxxzzz[i] = 3.0 * ts_yyyzz_xxzzz[i] * fe_0 + ts_yyyzz_xxxzzz[i] * pa_x[i];

        ts_xyyyzz_xxyyyy[i] = 2.0 * ts_yyyzz_xyyyy[i] * fe_0 + ts_yyyzz_xxyyyy[i] * pa_x[i];

        ts_xyyyzz_xxyyyz[i] = 2.0 * ts_yyyzz_xyyyz[i] * fe_0 + ts_yyyzz_xxyyyz[i] * pa_x[i];

        ts_xyyyzz_xxyyzz[i] = 2.0 * ts_yyyzz_xyyzz[i] * fe_0 + ts_yyyzz_xxyyzz[i] * pa_x[i];

        ts_xyyyzz_xxyzzz[i] = 2.0 * ts_yyyzz_xyzzz[i] * fe_0 + ts_yyyzz_xxyzzz[i] * pa_x[i];

        ts_xyyyzz_xxzzzz[i] = 2.0 * ts_yyyzz_xzzzz[i] * fe_0 + ts_yyyzz_xxzzzz[i] * pa_x[i];

        ts_xyyyzz_xyyyyy[i] = ts_yyyzz_yyyyy[i] * fe_0 + ts_yyyzz_xyyyyy[i] * pa_x[i];

        ts_xyyyzz_xyyyyz[i] = ts_yyyzz_yyyyz[i] * fe_0 + ts_yyyzz_xyyyyz[i] * pa_x[i];

        ts_xyyyzz_xyyyzz[i] = ts_yyyzz_yyyzz[i] * fe_0 + ts_yyyzz_xyyyzz[i] * pa_x[i];

        ts_xyyyzz_xyyzzz[i] = ts_yyyzz_yyzzz[i] * fe_0 + ts_yyyzz_xyyzzz[i] * pa_x[i];

        ts_xyyyzz_xyzzzz[i] = ts_yyyzz_yzzzz[i] * fe_0 + ts_yyyzz_xyzzzz[i] * pa_x[i];

        ts_xyyyzz_xzzzzz[i] = ts_yyyzz_zzzzz[i] * fe_0 + ts_yyyzz_xzzzzz[i] * pa_x[i];

        ts_xyyyzz_yyyyyy[i] = ts_yyyzz_yyyyyy[i] * pa_x[i];

        ts_xyyyzz_yyyyyz[i] = ts_yyyzz_yyyyyz[i] * pa_x[i];

        ts_xyyyzz_yyyyzz[i] = ts_yyyzz_yyyyzz[i] * pa_x[i];

        ts_xyyyzz_yyyzzz[i] = ts_yyyzz_yyyzzz[i] * pa_x[i];

        ts_xyyyzz_yyzzzz[i] = ts_yyyzz_yyzzzz[i] * pa_x[i];

        ts_xyyyzz_yzzzzz[i] = ts_yyyzz_yzzzzz[i] * pa_x[i];

        ts_xyyyzz_zzzzzz[i] = ts_yyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 504-532 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, ts_xyyzzz_xxxxxx, ts_xyyzzz_xxxxxy, ts_xyyzzz_xxxxxz, ts_xyyzzz_xxxxyy, ts_xyyzzz_xxxxyz, ts_xyyzzz_xxxxzz, ts_xyyzzz_xxxyyy, ts_xyyzzz_xxxyyz, ts_xyyzzz_xxxyzz, ts_xyyzzz_xxxzzz, ts_xyyzzz_xxyyyy, ts_xyyzzz_xxyyyz, ts_xyyzzz_xxyyzz, ts_xyyzzz_xxyzzz, ts_xyyzzz_xxzzzz, ts_xyyzzz_xyyyyy, ts_xyyzzz_xyyyyz, ts_xyyzzz_xyyyzz, ts_xyyzzz_xyyzzz, ts_xyyzzz_xyzzzz, ts_xyyzzz_xzzzzz, ts_xyyzzz_yyyyyy, ts_xyyzzz_yyyyyz, ts_xyyzzz_yyyyzz, ts_xyyzzz_yyyzzz, ts_xyyzzz_yyzzzz, ts_xyyzzz_yzzzzz, ts_xyyzzz_zzzzzz, ts_yyzzz_xxxxx, ts_yyzzz_xxxxxx, ts_yyzzz_xxxxxy, ts_yyzzz_xxxxxz, ts_yyzzz_xxxxy, ts_yyzzz_xxxxyy, ts_yyzzz_xxxxyz, ts_yyzzz_xxxxz, ts_yyzzz_xxxxzz, ts_yyzzz_xxxyy, ts_yyzzz_xxxyyy, ts_yyzzz_xxxyyz, ts_yyzzz_xxxyz, ts_yyzzz_xxxyzz, ts_yyzzz_xxxzz, ts_yyzzz_xxxzzz, ts_yyzzz_xxyyy, ts_yyzzz_xxyyyy, ts_yyzzz_xxyyyz, ts_yyzzz_xxyyz, ts_yyzzz_xxyyzz, ts_yyzzz_xxyzz, ts_yyzzz_xxyzzz, ts_yyzzz_xxzzz, ts_yyzzz_xxzzzz, ts_yyzzz_xyyyy, ts_yyzzz_xyyyyy, ts_yyzzz_xyyyyz, ts_yyzzz_xyyyz, ts_yyzzz_xyyyzz, ts_yyzzz_xyyzz, ts_yyzzz_xyyzzz, ts_yyzzz_xyzzz, ts_yyzzz_xyzzzz, ts_yyzzz_xzzzz, ts_yyzzz_xzzzzz, ts_yyzzz_yyyyy, ts_yyzzz_yyyyyy, ts_yyzzz_yyyyyz, ts_yyzzz_yyyyz, ts_yyzzz_yyyyzz, ts_yyzzz_yyyzz, ts_yyzzz_yyyzzz, ts_yyzzz_yyzzz, ts_yyzzz_yyzzzz, ts_yyzzz_yzzzz, ts_yyzzz_yzzzzz, ts_yyzzz_zzzzz, ts_yyzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzzz_xxxxxx[i] = 6.0 * ts_yyzzz_xxxxx[i] * fe_0 + ts_yyzzz_xxxxxx[i] * pa_x[i];

        ts_xyyzzz_xxxxxy[i] = 5.0 * ts_yyzzz_xxxxy[i] * fe_0 + ts_yyzzz_xxxxxy[i] * pa_x[i];

        ts_xyyzzz_xxxxxz[i] = 5.0 * ts_yyzzz_xxxxz[i] * fe_0 + ts_yyzzz_xxxxxz[i] * pa_x[i];

        ts_xyyzzz_xxxxyy[i] = 4.0 * ts_yyzzz_xxxyy[i] * fe_0 + ts_yyzzz_xxxxyy[i] * pa_x[i];

        ts_xyyzzz_xxxxyz[i] = 4.0 * ts_yyzzz_xxxyz[i] * fe_0 + ts_yyzzz_xxxxyz[i] * pa_x[i];

        ts_xyyzzz_xxxxzz[i] = 4.0 * ts_yyzzz_xxxzz[i] * fe_0 + ts_yyzzz_xxxxzz[i] * pa_x[i];

        ts_xyyzzz_xxxyyy[i] = 3.0 * ts_yyzzz_xxyyy[i] * fe_0 + ts_yyzzz_xxxyyy[i] * pa_x[i];

        ts_xyyzzz_xxxyyz[i] = 3.0 * ts_yyzzz_xxyyz[i] * fe_0 + ts_yyzzz_xxxyyz[i] * pa_x[i];

        ts_xyyzzz_xxxyzz[i] = 3.0 * ts_yyzzz_xxyzz[i] * fe_0 + ts_yyzzz_xxxyzz[i] * pa_x[i];

        ts_xyyzzz_xxxzzz[i] = 3.0 * ts_yyzzz_xxzzz[i] * fe_0 + ts_yyzzz_xxxzzz[i] * pa_x[i];

        ts_xyyzzz_xxyyyy[i] = 2.0 * ts_yyzzz_xyyyy[i] * fe_0 + ts_yyzzz_xxyyyy[i] * pa_x[i];

        ts_xyyzzz_xxyyyz[i] = 2.0 * ts_yyzzz_xyyyz[i] * fe_0 + ts_yyzzz_xxyyyz[i] * pa_x[i];

        ts_xyyzzz_xxyyzz[i] = 2.0 * ts_yyzzz_xyyzz[i] * fe_0 + ts_yyzzz_xxyyzz[i] * pa_x[i];

        ts_xyyzzz_xxyzzz[i] = 2.0 * ts_yyzzz_xyzzz[i] * fe_0 + ts_yyzzz_xxyzzz[i] * pa_x[i];

        ts_xyyzzz_xxzzzz[i] = 2.0 * ts_yyzzz_xzzzz[i] * fe_0 + ts_yyzzz_xxzzzz[i] * pa_x[i];

        ts_xyyzzz_xyyyyy[i] = ts_yyzzz_yyyyy[i] * fe_0 + ts_yyzzz_xyyyyy[i] * pa_x[i];

        ts_xyyzzz_xyyyyz[i] = ts_yyzzz_yyyyz[i] * fe_0 + ts_yyzzz_xyyyyz[i] * pa_x[i];

        ts_xyyzzz_xyyyzz[i] = ts_yyzzz_yyyzz[i] * fe_0 + ts_yyzzz_xyyyzz[i] * pa_x[i];

        ts_xyyzzz_xyyzzz[i] = ts_yyzzz_yyzzz[i] * fe_0 + ts_yyzzz_xyyzzz[i] * pa_x[i];

        ts_xyyzzz_xyzzzz[i] = ts_yyzzz_yzzzz[i] * fe_0 + ts_yyzzz_xyzzzz[i] * pa_x[i];

        ts_xyyzzz_xzzzzz[i] = ts_yyzzz_zzzzz[i] * fe_0 + ts_yyzzz_xzzzzz[i] * pa_x[i];

        ts_xyyzzz_yyyyyy[i] = ts_yyzzz_yyyyyy[i] * pa_x[i];

        ts_xyyzzz_yyyyyz[i] = ts_yyzzz_yyyyyz[i] * pa_x[i];

        ts_xyyzzz_yyyyzz[i] = ts_yyzzz_yyyyzz[i] * pa_x[i];

        ts_xyyzzz_yyyzzz[i] = ts_yyzzz_yyyzzz[i] * pa_x[i];

        ts_xyyzzz_yyzzzz[i] = ts_yyzzz_yyzzzz[i] * pa_x[i];

        ts_xyyzzz_yzzzzz[i] = ts_yyzzz_yzzzzz[i] * pa_x[i];

        ts_xyyzzz_zzzzzz[i] = ts_yyzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 532-560 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzzz_xxxxxx, ts_xyzzzz_xxxxxy, ts_xyzzzz_xxxxxz, ts_xyzzzz_xxxxyy, ts_xyzzzz_xxxxyz, ts_xyzzzz_xxxxzz, ts_xyzzzz_xxxyyy, ts_xyzzzz_xxxyyz, ts_xyzzzz_xxxyzz, ts_xyzzzz_xxxzzz, ts_xyzzzz_xxyyyy, ts_xyzzzz_xxyyyz, ts_xyzzzz_xxyyzz, ts_xyzzzz_xxyzzz, ts_xyzzzz_xxzzzz, ts_xyzzzz_xyyyyy, ts_xyzzzz_xyyyyz, ts_xyzzzz_xyyyzz, ts_xyzzzz_xyyzzz, ts_xyzzzz_xyzzzz, ts_xyzzzz_xzzzzz, ts_xyzzzz_yyyyyy, ts_xyzzzz_yyyyyz, ts_xyzzzz_yyyyzz, ts_xyzzzz_yyyzzz, ts_xyzzzz_yyzzzz, ts_xyzzzz_yzzzzz, ts_xyzzzz_zzzzzz, ts_xzzzz_xxxxxx, ts_xzzzz_xxxxxz, ts_xzzzz_xxxxzz, ts_xzzzz_xxxzzz, ts_xzzzz_xxzzzz, ts_xzzzz_xzzzzz, ts_yzzzz_xxxxxy, ts_yzzzz_xxxxy, ts_yzzzz_xxxxyy, ts_yzzzz_xxxxyz, ts_yzzzz_xxxyy, ts_yzzzz_xxxyyy, ts_yzzzz_xxxyyz, ts_yzzzz_xxxyz, ts_yzzzz_xxxyzz, ts_yzzzz_xxyyy, ts_yzzzz_xxyyyy, ts_yzzzz_xxyyyz, ts_yzzzz_xxyyz, ts_yzzzz_xxyyzz, ts_yzzzz_xxyzz, ts_yzzzz_xxyzzz, ts_yzzzz_xyyyy, ts_yzzzz_xyyyyy, ts_yzzzz_xyyyyz, ts_yzzzz_xyyyz, ts_yzzzz_xyyyzz, ts_yzzzz_xyyzz, ts_yzzzz_xyyzzz, ts_yzzzz_xyzzz, ts_yzzzz_xyzzzz, ts_yzzzz_yyyyy, ts_yzzzz_yyyyyy, ts_yzzzz_yyyyyz, ts_yzzzz_yyyyz, ts_yzzzz_yyyyzz, ts_yzzzz_yyyzz, ts_yzzzz_yyyzzz, ts_yzzzz_yyzzz, ts_yzzzz_yyzzzz, ts_yzzzz_yzzzz, ts_yzzzz_yzzzzz, ts_yzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzzz_xxxxxx[i] = ts_xzzzz_xxxxxx[i] * pa_y[i];

        ts_xyzzzz_xxxxxy[i] = 5.0 * ts_yzzzz_xxxxy[i] * fe_0 + ts_yzzzz_xxxxxy[i] * pa_x[i];

        ts_xyzzzz_xxxxxz[i] = ts_xzzzz_xxxxxz[i] * pa_y[i];

        ts_xyzzzz_xxxxyy[i] = 4.0 * ts_yzzzz_xxxyy[i] * fe_0 + ts_yzzzz_xxxxyy[i] * pa_x[i];

        ts_xyzzzz_xxxxyz[i] = 4.0 * ts_yzzzz_xxxyz[i] * fe_0 + ts_yzzzz_xxxxyz[i] * pa_x[i];

        ts_xyzzzz_xxxxzz[i] = ts_xzzzz_xxxxzz[i] * pa_y[i];

        ts_xyzzzz_xxxyyy[i] = 3.0 * ts_yzzzz_xxyyy[i] * fe_0 + ts_yzzzz_xxxyyy[i] * pa_x[i];

        ts_xyzzzz_xxxyyz[i] = 3.0 * ts_yzzzz_xxyyz[i] * fe_0 + ts_yzzzz_xxxyyz[i] * pa_x[i];

        ts_xyzzzz_xxxyzz[i] = 3.0 * ts_yzzzz_xxyzz[i] * fe_0 + ts_yzzzz_xxxyzz[i] * pa_x[i];

        ts_xyzzzz_xxxzzz[i] = ts_xzzzz_xxxzzz[i] * pa_y[i];

        ts_xyzzzz_xxyyyy[i] = 2.0 * ts_yzzzz_xyyyy[i] * fe_0 + ts_yzzzz_xxyyyy[i] * pa_x[i];

        ts_xyzzzz_xxyyyz[i] = 2.0 * ts_yzzzz_xyyyz[i] * fe_0 + ts_yzzzz_xxyyyz[i] * pa_x[i];

        ts_xyzzzz_xxyyzz[i] = 2.0 * ts_yzzzz_xyyzz[i] * fe_0 + ts_yzzzz_xxyyzz[i] * pa_x[i];

        ts_xyzzzz_xxyzzz[i] = 2.0 * ts_yzzzz_xyzzz[i] * fe_0 + ts_yzzzz_xxyzzz[i] * pa_x[i];

        ts_xyzzzz_xxzzzz[i] = ts_xzzzz_xxzzzz[i] * pa_y[i];

        ts_xyzzzz_xyyyyy[i] = ts_yzzzz_yyyyy[i] * fe_0 + ts_yzzzz_xyyyyy[i] * pa_x[i];

        ts_xyzzzz_xyyyyz[i] = ts_yzzzz_yyyyz[i] * fe_0 + ts_yzzzz_xyyyyz[i] * pa_x[i];

        ts_xyzzzz_xyyyzz[i] = ts_yzzzz_yyyzz[i] * fe_0 + ts_yzzzz_xyyyzz[i] * pa_x[i];

        ts_xyzzzz_xyyzzz[i] = ts_yzzzz_yyzzz[i] * fe_0 + ts_yzzzz_xyyzzz[i] * pa_x[i];

        ts_xyzzzz_xyzzzz[i] = ts_yzzzz_yzzzz[i] * fe_0 + ts_yzzzz_xyzzzz[i] * pa_x[i];

        ts_xyzzzz_xzzzzz[i] = ts_xzzzz_xzzzzz[i] * pa_y[i];

        ts_xyzzzz_yyyyyy[i] = ts_yzzzz_yyyyyy[i] * pa_x[i];

        ts_xyzzzz_yyyyyz[i] = ts_yzzzz_yyyyyz[i] * pa_x[i];

        ts_xyzzzz_yyyyzz[i] = ts_yzzzz_yyyyzz[i] * pa_x[i];

        ts_xyzzzz_yyyzzz[i] = ts_yzzzz_yyyzzz[i] * pa_x[i];

        ts_xyzzzz_yyzzzz[i] = ts_yzzzz_yyzzzz[i] * pa_x[i];

        ts_xyzzzz_yzzzzz[i] = ts_yzzzz_yzzzzz[i] * pa_x[i];

        ts_xyzzzz_zzzzzz[i] = ts_yzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 560-588 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_x, ts_xzzzzz_xxxxxx, ts_xzzzzz_xxxxxy, ts_xzzzzz_xxxxxz, ts_xzzzzz_xxxxyy, ts_xzzzzz_xxxxyz, ts_xzzzzz_xxxxzz, ts_xzzzzz_xxxyyy, ts_xzzzzz_xxxyyz, ts_xzzzzz_xxxyzz, ts_xzzzzz_xxxzzz, ts_xzzzzz_xxyyyy, ts_xzzzzz_xxyyyz, ts_xzzzzz_xxyyzz, ts_xzzzzz_xxyzzz, ts_xzzzzz_xxzzzz, ts_xzzzzz_xyyyyy, ts_xzzzzz_xyyyyz, ts_xzzzzz_xyyyzz, ts_xzzzzz_xyyzzz, ts_xzzzzz_xyzzzz, ts_xzzzzz_xzzzzz, ts_xzzzzz_yyyyyy, ts_xzzzzz_yyyyyz, ts_xzzzzz_yyyyzz, ts_xzzzzz_yyyzzz, ts_xzzzzz_yyzzzz, ts_xzzzzz_yzzzzz, ts_xzzzzz_zzzzzz, ts_zzzzz_xxxxx, ts_zzzzz_xxxxxx, ts_zzzzz_xxxxxy, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxy, ts_zzzzz_xxxxyy, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyy, ts_zzzzz_xxxyyy, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyy, ts_zzzzz_xxyyyy, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyy, ts_zzzzz_xyyyyy, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyy, ts_zzzzz_yyyyyy, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzz, ts_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzzz_xxxxxx[i] = 6.0 * ts_zzzzz_xxxxx[i] * fe_0 + ts_zzzzz_xxxxxx[i] * pa_x[i];

        ts_xzzzzz_xxxxxy[i] = 5.0 * ts_zzzzz_xxxxy[i] * fe_0 + ts_zzzzz_xxxxxy[i] * pa_x[i];

        ts_xzzzzz_xxxxxz[i] = 5.0 * ts_zzzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxxxz[i] * pa_x[i];

        ts_xzzzzz_xxxxyy[i] = 4.0 * ts_zzzzz_xxxyy[i] * fe_0 + ts_zzzzz_xxxxyy[i] * pa_x[i];

        ts_xzzzzz_xxxxyz[i] = 4.0 * ts_zzzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxxyz[i] * pa_x[i];

        ts_xzzzzz_xxxxzz[i] = 4.0 * ts_zzzzz_xxxzz[i] * fe_0 + ts_zzzzz_xxxxzz[i] * pa_x[i];

        ts_xzzzzz_xxxyyy[i] = 3.0 * ts_zzzzz_xxyyy[i] * fe_0 + ts_zzzzz_xxxyyy[i] * pa_x[i];

        ts_xzzzzz_xxxyyz[i] = 3.0 * ts_zzzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxxyyz[i] * pa_x[i];

        ts_xzzzzz_xxxyzz[i] = 3.0 * ts_zzzzz_xxyzz[i] * fe_0 + ts_zzzzz_xxxyzz[i] * pa_x[i];

        ts_xzzzzz_xxxzzz[i] = 3.0 * ts_zzzzz_xxzzz[i] * fe_0 + ts_zzzzz_xxxzzz[i] * pa_x[i];

        ts_xzzzzz_xxyyyy[i] = 2.0 * ts_zzzzz_xyyyy[i] * fe_0 + ts_zzzzz_xxyyyy[i] * pa_x[i];

        ts_xzzzzz_xxyyyz[i] = 2.0 * ts_zzzzz_xyyyz[i] * fe_0 + ts_zzzzz_xxyyyz[i] * pa_x[i];

        ts_xzzzzz_xxyyzz[i] = 2.0 * ts_zzzzz_xyyzz[i] * fe_0 + ts_zzzzz_xxyyzz[i] * pa_x[i];

        ts_xzzzzz_xxyzzz[i] = 2.0 * ts_zzzzz_xyzzz[i] * fe_0 + ts_zzzzz_xxyzzz[i] * pa_x[i];

        ts_xzzzzz_xxzzzz[i] = 2.0 * ts_zzzzz_xzzzz[i] * fe_0 + ts_zzzzz_xxzzzz[i] * pa_x[i];

        ts_xzzzzz_xyyyyy[i] = ts_zzzzz_yyyyy[i] * fe_0 + ts_zzzzz_xyyyyy[i] * pa_x[i];

        ts_xzzzzz_xyyyyz[i] = ts_zzzzz_yyyyz[i] * fe_0 + ts_zzzzz_xyyyyz[i] * pa_x[i];

        ts_xzzzzz_xyyyzz[i] = ts_zzzzz_yyyzz[i] * fe_0 + ts_zzzzz_xyyyzz[i] * pa_x[i];

        ts_xzzzzz_xyyzzz[i] = ts_zzzzz_yyzzz[i] * fe_0 + ts_zzzzz_xyyzzz[i] * pa_x[i];

        ts_xzzzzz_xyzzzz[i] = ts_zzzzz_yzzzz[i] * fe_0 + ts_zzzzz_xyzzzz[i] * pa_x[i];

        ts_xzzzzz_xzzzzz[i] = ts_zzzzz_zzzzz[i] * fe_0 + ts_zzzzz_xzzzzz[i] * pa_x[i];

        ts_xzzzzz_yyyyyy[i] = ts_zzzzz_yyyyyy[i] * pa_x[i];

        ts_xzzzzz_yyyyyz[i] = ts_zzzzz_yyyyyz[i] * pa_x[i];

        ts_xzzzzz_yyyyzz[i] = ts_zzzzz_yyyyzz[i] * pa_x[i];

        ts_xzzzzz_yyyzzz[i] = ts_zzzzz_yyyzzz[i] * pa_x[i];

        ts_xzzzzz_yyzzzz[i] = ts_zzzzz_yyzzzz[i] * pa_x[i];

        ts_xzzzzz_yzzzzz[i] = ts_zzzzz_yzzzzz[i] * pa_x[i];

        ts_xzzzzz_zzzzzz[i] = ts_zzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 588-616 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_y, ts_yyyy_xxxxxx, ts_yyyy_xxxxxy, ts_yyyy_xxxxxz, ts_yyyy_xxxxyy, ts_yyyy_xxxxyz, ts_yyyy_xxxxzz, ts_yyyy_xxxyyy, ts_yyyy_xxxyyz, ts_yyyy_xxxyzz, ts_yyyy_xxxzzz, ts_yyyy_xxyyyy, ts_yyyy_xxyyyz, ts_yyyy_xxyyzz, ts_yyyy_xxyzzz, ts_yyyy_xxzzzz, ts_yyyy_xyyyyy, ts_yyyy_xyyyyz, ts_yyyy_xyyyzz, ts_yyyy_xyyzzz, ts_yyyy_xyzzzz, ts_yyyy_xzzzzz, ts_yyyy_yyyyyy, ts_yyyy_yyyyyz, ts_yyyy_yyyyzz, ts_yyyy_yyyzzz, ts_yyyy_yyzzzz, ts_yyyy_yzzzzz, ts_yyyy_zzzzzz, ts_yyyyy_xxxxx, ts_yyyyy_xxxxxx, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxxz, ts_yyyyy_xxxxy, ts_yyyyy_xxxxyy, ts_yyyyy_xxxxyz, ts_yyyyy_xxxxz, ts_yyyyy_xxxxzz, ts_yyyyy_xxxyy, ts_yyyyy_xxxyyy, ts_yyyyy_xxxyyz, ts_yyyyy_xxxyz, ts_yyyyy_xxxyzz, ts_yyyyy_xxxzz, ts_yyyyy_xxxzzz, ts_yyyyy_xxyyy, ts_yyyyy_xxyyyy, ts_yyyyy_xxyyyz, ts_yyyyy_xxyyz, ts_yyyyy_xxyyzz, ts_yyyyy_xxyzz, ts_yyyyy_xxyzzz, ts_yyyyy_xxzzz, ts_yyyyy_xxzzzz, ts_yyyyy_xyyyy, ts_yyyyy_xyyyyy, ts_yyyyy_xyyyyz, ts_yyyyy_xyyyz, ts_yyyyy_xyyyzz, ts_yyyyy_xyyzz, ts_yyyyy_xyyzzz, ts_yyyyy_xyzzz, ts_yyyyy_xyzzzz, ts_yyyyy_xzzzz, ts_yyyyy_xzzzzz, ts_yyyyy_yyyyy, ts_yyyyy_yyyyyy, ts_yyyyy_yyyyyz, ts_yyyyy_yyyyz, ts_yyyyy_yyyyzz, ts_yyyyy_yyyzz, ts_yyyyy_yyyzzz, ts_yyyyy_yyzzz, ts_yyyyy_yyzzzz, ts_yyyyy_yzzzz, ts_yyyyy_yzzzzz, ts_yyyyy_zzzzz, ts_yyyyy_zzzzzz, ts_yyyyyy_xxxxxx, ts_yyyyyy_xxxxxy, ts_yyyyyy_xxxxxz, ts_yyyyyy_xxxxyy, ts_yyyyyy_xxxxyz, ts_yyyyyy_xxxxzz, ts_yyyyyy_xxxyyy, ts_yyyyyy_xxxyyz, ts_yyyyyy_xxxyzz, ts_yyyyyy_xxxzzz, ts_yyyyyy_xxyyyy, ts_yyyyyy_xxyyyz, ts_yyyyyy_xxyyzz, ts_yyyyyy_xxyzzz, ts_yyyyyy_xxzzzz, ts_yyyyyy_xyyyyy, ts_yyyyyy_xyyyyz, ts_yyyyyy_xyyyzz, ts_yyyyyy_xyyzzz, ts_yyyyyy_xyzzzz, ts_yyyyyy_xzzzzz, ts_yyyyyy_yyyyyy, ts_yyyyyy_yyyyyz, ts_yyyyyy_yyyyzz, ts_yyyyyy_yyyzzz, ts_yyyyyy_yyzzzz, ts_yyyyyy_yzzzzz, ts_yyyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyy_xxxxxx[i] = 5.0 * ts_yyyy_xxxxxx[i] * fe_0 + ts_yyyyy_xxxxxx[i] * pa_y[i];

        ts_yyyyyy_xxxxxy[i] = 5.0 * ts_yyyy_xxxxxy[i] * fe_0 + ts_yyyyy_xxxxx[i] * fe_0 + ts_yyyyy_xxxxxy[i] * pa_y[i];

        ts_yyyyyy_xxxxxz[i] = 5.0 * ts_yyyy_xxxxxz[i] * fe_0 + ts_yyyyy_xxxxxz[i] * pa_y[i];

        ts_yyyyyy_xxxxyy[i] = 5.0 * ts_yyyy_xxxxyy[i] * fe_0 + 2.0 * ts_yyyyy_xxxxy[i] * fe_0 + ts_yyyyy_xxxxyy[i] * pa_y[i];

        ts_yyyyyy_xxxxyz[i] = 5.0 * ts_yyyy_xxxxyz[i] * fe_0 + ts_yyyyy_xxxxz[i] * fe_0 + ts_yyyyy_xxxxyz[i] * pa_y[i];

        ts_yyyyyy_xxxxzz[i] = 5.0 * ts_yyyy_xxxxzz[i] * fe_0 + ts_yyyyy_xxxxzz[i] * pa_y[i];

        ts_yyyyyy_xxxyyy[i] = 5.0 * ts_yyyy_xxxyyy[i] * fe_0 + 3.0 * ts_yyyyy_xxxyy[i] * fe_0 + ts_yyyyy_xxxyyy[i] * pa_y[i];

        ts_yyyyyy_xxxyyz[i] = 5.0 * ts_yyyy_xxxyyz[i] * fe_0 + 2.0 * ts_yyyyy_xxxyz[i] * fe_0 + ts_yyyyy_xxxyyz[i] * pa_y[i];

        ts_yyyyyy_xxxyzz[i] = 5.0 * ts_yyyy_xxxyzz[i] * fe_0 + ts_yyyyy_xxxzz[i] * fe_0 + ts_yyyyy_xxxyzz[i] * pa_y[i];

        ts_yyyyyy_xxxzzz[i] = 5.0 * ts_yyyy_xxxzzz[i] * fe_0 + ts_yyyyy_xxxzzz[i] * pa_y[i];

        ts_yyyyyy_xxyyyy[i] = 5.0 * ts_yyyy_xxyyyy[i] * fe_0 + 4.0 * ts_yyyyy_xxyyy[i] * fe_0 + ts_yyyyy_xxyyyy[i] * pa_y[i];

        ts_yyyyyy_xxyyyz[i] = 5.0 * ts_yyyy_xxyyyz[i] * fe_0 + 3.0 * ts_yyyyy_xxyyz[i] * fe_0 + ts_yyyyy_xxyyyz[i] * pa_y[i];

        ts_yyyyyy_xxyyzz[i] = 5.0 * ts_yyyy_xxyyzz[i] * fe_0 + 2.0 * ts_yyyyy_xxyzz[i] * fe_0 + ts_yyyyy_xxyyzz[i] * pa_y[i];

        ts_yyyyyy_xxyzzz[i] = 5.0 * ts_yyyy_xxyzzz[i] * fe_0 + ts_yyyyy_xxzzz[i] * fe_0 + ts_yyyyy_xxyzzz[i] * pa_y[i];

        ts_yyyyyy_xxzzzz[i] = 5.0 * ts_yyyy_xxzzzz[i] * fe_0 + ts_yyyyy_xxzzzz[i] * pa_y[i];

        ts_yyyyyy_xyyyyy[i] = 5.0 * ts_yyyy_xyyyyy[i] * fe_0 + 5.0 * ts_yyyyy_xyyyy[i] * fe_0 + ts_yyyyy_xyyyyy[i] * pa_y[i];

        ts_yyyyyy_xyyyyz[i] = 5.0 * ts_yyyy_xyyyyz[i] * fe_0 + 4.0 * ts_yyyyy_xyyyz[i] * fe_0 + ts_yyyyy_xyyyyz[i] * pa_y[i];

        ts_yyyyyy_xyyyzz[i] = 5.0 * ts_yyyy_xyyyzz[i] * fe_0 + 3.0 * ts_yyyyy_xyyzz[i] * fe_0 + ts_yyyyy_xyyyzz[i] * pa_y[i];

        ts_yyyyyy_xyyzzz[i] = 5.0 * ts_yyyy_xyyzzz[i] * fe_0 + 2.0 * ts_yyyyy_xyzzz[i] * fe_0 + ts_yyyyy_xyyzzz[i] * pa_y[i];

        ts_yyyyyy_xyzzzz[i] = 5.0 * ts_yyyy_xyzzzz[i] * fe_0 + ts_yyyyy_xzzzz[i] * fe_0 + ts_yyyyy_xyzzzz[i] * pa_y[i];

        ts_yyyyyy_xzzzzz[i] = 5.0 * ts_yyyy_xzzzzz[i] * fe_0 + ts_yyyyy_xzzzzz[i] * pa_y[i];

        ts_yyyyyy_yyyyyy[i] = 5.0 * ts_yyyy_yyyyyy[i] * fe_0 + 6.0 * ts_yyyyy_yyyyy[i] * fe_0 + ts_yyyyy_yyyyyy[i] * pa_y[i];

        ts_yyyyyy_yyyyyz[i] = 5.0 * ts_yyyy_yyyyyz[i] * fe_0 + 5.0 * ts_yyyyy_yyyyz[i] * fe_0 + ts_yyyyy_yyyyyz[i] * pa_y[i];

        ts_yyyyyy_yyyyzz[i] = 5.0 * ts_yyyy_yyyyzz[i] * fe_0 + 4.0 * ts_yyyyy_yyyzz[i] * fe_0 + ts_yyyyy_yyyyzz[i] * pa_y[i];

        ts_yyyyyy_yyyzzz[i] = 5.0 * ts_yyyy_yyyzzz[i] * fe_0 + 3.0 * ts_yyyyy_yyzzz[i] * fe_0 + ts_yyyyy_yyyzzz[i] * pa_y[i];

        ts_yyyyyy_yyzzzz[i] = 5.0 * ts_yyyy_yyzzzz[i] * fe_0 + 2.0 * ts_yyyyy_yzzzz[i] * fe_0 + ts_yyyyy_yyzzzz[i] * pa_y[i];

        ts_yyyyyy_yzzzzz[i] = 5.0 * ts_yyyy_yzzzzz[i] * fe_0 + ts_yyyyy_zzzzz[i] * fe_0 + ts_yyyyy_yzzzzz[i] * pa_y[i];

        ts_yyyyyy_zzzzzz[i] = 5.0 * ts_yyyy_zzzzzz[i] * fe_0 + ts_yyyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 616-644 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyyy_xxxxxx, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxy, ts_yyyyy_xxxxyy, ts_yyyyy_xxxxyz, ts_yyyyy_xxxyy, ts_yyyyy_xxxyyy, ts_yyyyy_xxxyyz, ts_yyyyy_xxxyz, ts_yyyyy_xxxyzz, ts_yyyyy_xxyyy, ts_yyyyy_xxyyyy, ts_yyyyy_xxyyyz, ts_yyyyy_xxyyz, ts_yyyyy_xxyyzz, ts_yyyyy_xxyzz, ts_yyyyy_xxyzzz, ts_yyyyy_xyyyy, ts_yyyyy_xyyyyy, ts_yyyyy_xyyyyz, ts_yyyyy_xyyyz, ts_yyyyy_xyyyzz, ts_yyyyy_xyyzz, ts_yyyyy_xyyzzz, ts_yyyyy_xyzzz, ts_yyyyy_xyzzzz, ts_yyyyy_yyyyy, ts_yyyyy_yyyyyy, ts_yyyyy_yyyyyz, ts_yyyyy_yyyyz, ts_yyyyy_yyyyzz, ts_yyyyy_yyyzz, ts_yyyyy_yyyzzz, ts_yyyyy_yyzzz, ts_yyyyy_yyzzzz, ts_yyyyy_yzzzz, ts_yyyyy_yzzzzz, ts_yyyyyz_xxxxxx, ts_yyyyyz_xxxxxy, ts_yyyyyz_xxxxxz, ts_yyyyyz_xxxxyy, ts_yyyyyz_xxxxyz, ts_yyyyyz_xxxxzz, ts_yyyyyz_xxxyyy, ts_yyyyyz_xxxyyz, ts_yyyyyz_xxxyzz, ts_yyyyyz_xxxzzz, ts_yyyyyz_xxyyyy, ts_yyyyyz_xxyyyz, ts_yyyyyz_xxyyzz, ts_yyyyyz_xxyzzz, ts_yyyyyz_xxzzzz, ts_yyyyyz_xyyyyy, ts_yyyyyz_xyyyyz, ts_yyyyyz_xyyyzz, ts_yyyyyz_xyyzzz, ts_yyyyyz_xyzzzz, ts_yyyyyz_xzzzzz, ts_yyyyyz_yyyyyy, ts_yyyyyz_yyyyyz, ts_yyyyyz_yyyyzz, ts_yyyyyz_yyyzzz, ts_yyyyyz_yyzzzz, ts_yyyyyz_yzzzzz, ts_yyyyyz_zzzzzz, ts_yyyyz_xxxxxz, ts_yyyyz_xxxxzz, ts_yyyyz_xxxzzz, ts_yyyyz_xxzzzz, ts_yyyyz_xzzzzz, ts_yyyyz_zzzzzz, ts_yyyz_xxxxxz, ts_yyyz_xxxxzz, ts_yyyz_xxxzzz, ts_yyyz_xxzzzz, ts_yyyz_xzzzzz, ts_yyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyz_xxxxxx[i] = ts_yyyyy_xxxxxx[i] * pa_z[i];

        ts_yyyyyz_xxxxxy[i] = ts_yyyyy_xxxxxy[i] * pa_z[i];

        ts_yyyyyz_xxxxxz[i] = 4.0 * ts_yyyz_xxxxxz[i] * fe_0 + ts_yyyyz_xxxxxz[i] * pa_y[i];

        ts_yyyyyz_xxxxyy[i] = ts_yyyyy_xxxxyy[i] * pa_z[i];

        ts_yyyyyz_xxxxyz[i] = ts_yyyyy_xxxxy[i] * fe_0 + ts_yyyyy_xxxxyz[i] * pa_z[i];

        ts_yyyyyz_xxxxzz[i] = 4.0 * ts_yyyz_xxxxzz[i] * fe_0 + ts_yyyyz_xxxxzz[i] * pa_y[i];

        ts_yyyyyz_xxxyyy[i] = ts_yyyyy_xxxyyy[i] * pa_z[i];

        ts_yyyyyz_xxxyyz[i] = ts_yyyyy_xxxyy[i] * fe_0 + ts_yyyyy_xxxyyz[i] * pa_z[i];

        ts_yyyyyz_xxxyzz[i] = 2.0 * ts_yyyyy_xxxyz[i] * fe_0 + ts_yyyyy_xxxyzz[i] * pa_z[i];

        ts_yyyyyz_xxxzzz[i] = 4.0 * ts_yyyz_xxxzzz[i] * fe_0 + ts_yyyyz_xxxzzz[i] * pa_y[i];

        ts_yyyyyz_xxyyyy[i] = ts_yyyyy_xxyyyy[i] * pa_z[i];

        ts_yyyyyz_xxyyyz[i] = ts_yyyyy_xxyyy[i] * fe_0 + ts_yyyyy_xxyyyz[i] * pa_z[i];

        ts_yyyyyz_xxyyzz[i] = 2.0 * ts_yyyyy_xxyyz[i] * fe_0 + ts_yyyyy_xxyyzz[i] * pa_z[i];

        ts_yyyyyz_xxyzzz[i] = 3.0 * ts_yyyyy_xxyzz[i] * fe_0 + ts_yyyyy_xxyzzz[i] * pa_z[i];

        ts_yyyyyz_xxzzzz[i] = 4.0 * ts_yyyz_xxzzzz[i] * fe_0 + ts_yyyyz_xxzzzz[i] * pa_y[i];

        ts_yyyyyz_xyyyyy[i] = ts_yyyyy_xyyyyy[i] * pa_z[i];

        ts_yyyyyz_xyyyyz[i] = ts_yyyyy_xyyyy[i] * fe_0 + ts_yyyyy_xyyyyz[i] * pa_z[i];

        ts_yyyyyz_xyyyzz[i] = 2.0 * ts_yyyyy_xyyyz[i] * fe_0 + ts_yyyyy_xyyyzz[i] * pa_z[i];

        ts_yyyyyz_xyyzzz[i] = 3.0 * ts_yyyyy_xyyzz[i] * fe_0 + ts_yyyyy_xyyzzz[i] * pa_z[i];

        ts_yyyyyz_xyzzzz[i] = 4.0 * ts_yyyyy_xyzzz[i] * fe_0 + ts_yyyyy_xyzzzz[i] * pa_z[i];

        ts_yyyyyz_xzzzzz[i] = 4.0 * ts_yyyz_xzzzzz[i] * fe_0 + ts_yyyyz_xzzzzz[i] * pa_y[i];

        ts_yyyyyz_yyyyyy[i] = ts_yyyyy_yyyyyy[i] * pa_z[i];

        ts_yyyyyz_yyyyyz[i] = ts_yyyyy_yyyyy[i] * fe_0 + ts_yyyyy_yyyyyz[i] * pa_z[i];

        ts_yyyyyz_yyyyzz[i] = 2.0 * ts_yyyyy_yyyyz[i] * fe_0 + ts_yyyyy_yyyyzz[i] * pa_z[i];

        ts_yyyyyz_yyyzzz[i] = 3.0 * ts_yyyyy_yyyzz[i] * fe_0 + ts_yyyyy_yyyzzz[i] * pa_z[i];

        ts_yyyyyz_yyzzzz[i] = 4.0 * ts_yyyyy_yyzzz[i] * fe_0 + ts_yyyyy_yyzzzz[i] * pa_z[i];

        ts_yyyyyz_yzzzzz[i] = 5.0 * ts_yyyyy_yzzzz[i] * fe_0 + ts_yyyyy_yzzzzz[i] * pa_z[i];

        ts_yyyyyz_zzzzzz[i] = 4.0 * ts_yyyz_zzzzzz[i] * fe_0 + ts_yyyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 644-672 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xxxxxy, ts_yyyy_xxxxyy, ts_yyyy_xxxyyy, ts_yyyy_xxyyyy, ts_yyyy_xyyyyy, ts_yyyy_yyyyyy, ts_yyyyz_xxxxxy, ts_yyyyz_xxxxyy, ts_yyyyz_xxxyyy, ts_yyyyz_xxyyyy, ts_yyyyz_xyyyyy, ts_yyyyz_yyyyyy, ts_yyyyzz_xxxxxx, ts_yyyyzz_xxxxxy, ts_yyyyzz_xxxxxz, ts_yyyyzz_xxxxyy, ts_yyyyzz_xxxxyz, ts_yyyyzz_xxxxzz, ts_yyyyzz_xxxyyy, ts_yyyyzz_xxxyyz, ts_yyyyzz_xxxyzz, ts_yyyyzz_xxxzzz, ts_yyyyzz_xxyyyy, ts_yyyyzz_xxyyyz, ts_yyyyzz_xxyyzz, ts_yyyyzz_xxyzzz, ts_yyyyzz_xxzzzz, ts_yyyyzz_xyyyyy, ts_yyyyzz_xyyyyz, ts_yyyyzz_xyyyzz, ts_yyyyzz_xyyzzz, ts_yyyyzz_xyzzzz, ts_yyyyzz_xzzzzz, ts_yyyyzz_yyyyyy, ts_yyyyzz_yyyyyz, ts_yyyyzz_yyyyzz, ts_yyyyzz_yyyzzz, ts_yyyyzz_yyzzzz, ts_yyyyzz_yzzzzz, ts_yyyyzz_zzzzzz, ts_yyyzz_xxxxxx, ts_yyyzz_xxxxxz, ts_yyyzz_xxxxyz, ts_yyyzz_xxxxz, ts_yyyzz_xxxxzz, ts_yyyzz_xxxyyz, ts_yyyzz_xxxyz, ts_yyyzz_xxxyzz, ts_yyyzz_xxxzz, ts_yyyzz_xxxzzz, ts_yyyzz_xxyyyz, ts_yyyzz_xxyyz, ts_yyyzz_xxyyzz, ts_yyyzz_xxyzz, ts_yyyzz_xxyzzz, ts_yyyzz_xxzzz, ts_yyyzz_xxzzzz, ts_yyyzz_xyyyyz, ts_yyyzz_xyyyz, ts_yyyzz_xyyyzz, ts_yyyzz_xyyzz, ts_yyyzz_xyyzzz, ts_yyyzz_xyzzz, ts_yyyzz_xyzzzz, ts_yyyzz_xzzzz, ts_yyyzz_xzzzzz, ts_yyyzz_yyyyyz, ts_yyyzz_yyyyz, ts_yyyzz_yyyyzz, ts_yyyzz_yyyzz, ts_yyyzz_yyyzzz, ts_yyyzz_yyzzz, ts_yyyzz_yyzzzz, ts_yyyzz_yzzzz, ts_yyyzz_yzzzzz, ts_yyyzz_zzzzz, ts_yyyzz_zzzzzz, ts_yyzz_xxxxxx, ts_yyzz_xxxxxz, ts_yyzz_xxxxyz, ts_yyzz_xxxxzz, ts_yyzz_xxxyyz, ts_yyzz_xxxyzz, ts_yyzz_xxxzzz, ts_yyzz_xxyyyz, ts_yyzz_xxyyzz, ts_yyzz_xxyzzz, ts_yyzz_xxzzzz, ts_yyzz_xyyyyz, ts_yyzz_xyyyzz, ts_yyzz_xyyzzz, ts_yyzz_xyzzzz, ts_yyzz_xzzzzz, ts_yyzz_yyyyyz, ts_yyzz_yyyyzz, ts_yyzz_yyyzzz, ts_yyzz_yyzzzz, ts_yyzz_yzzzzz, ts_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyzz_xxxxxx[i] = 3.0 * ts_yyzz_xxxxxx[i] * fe_0 + ts_yyyzz_xxxxxx[i] * pa_y[i];

        ts_yyyyzz_xxxxxy[i] = ts_yyyy_xxxxxy[i] * fe_0 + ts_yyyyz_xxxxxy[i] * pa_z[i];

        ts_yyyyzz_xxxxxz[i] = 3.0 * ts_yyzz_xxxxxz[i] * fe_0 + ts_yyyzz_xxxxxz[i] * pa_y[i];

        ts_yyyyzz_xxxxyy[i] = ts_yyyy_xxxxyy[i] * fe_0 + ts_yyyyz_xxxxyy[i] * pa_z[i];

        ts_yyyyzz_xxxxyz[i] = 3.0 * ts_yyzz_xxxxyz[i] * fe_0 + ts_yyyzz_xxxxz[i] * fe_0 + ts_yyyzz_xxxxyz[i] * pa_y[i];

        ts_yyyyzz_xxxxzz[i] = 3.0 * ts_yyzz_xxxxzz[i] * fe_0 + ts_yyyzz_xxxxzz[i] * pa_y[i];

        ts_yyyyzz_xxxyyy[i] = ts_yyyy_xxxyyy[i] * fe_0 + ts_yyyyz_xxxyyy[i] * pa_z[i];

        ts_yyyyzz_xxxyyz[i] = 3.0 * ts_yyzz_xxxyyz[i] * fe_0 + 2.0 * ts_yyyzz_xxxyz[i] * fe_0 + ts_yyyzz_xxxyyz[i] * pa_y[i];

        ts_yyyyzz_xxxyzz[i] = 3.0 * ts_yyzz_xxxyzz[i] * fe_0 + ts_yyyzz_xxxzz[i] * fe_0 + ts_yyyzz_xxxyzz[i] * pa_y[i];

        ts_yyyyzz_xxxzzz[i] = 3.0 * ts_yyzz_xxxzzz[i] * fe_0 + ts_yyyzz_xxxzzz[i] * pa_y[i];

        ts_yyyyzz_xxyyyy[i] = ts_yyyy_xxyyyy[i] * fe_0 + ts_yyyyz_xxyyyy[i] * pa_z[i];

        ts_yyyyzz_xxyyyz[i] = 3.0 * ts_yyzz_xxyyyz[i] * fe_0 + 3.0 * ts_yyyzz_xxyyz[i] * fe_0 + ts_yyyzz_xxyyyz[i] * pa_y[i];

        ts_yyyyzz_xxyyzz[i] = 3.0 * ts_yyzz_xxyyzz[i] * fe_0 + 2.0 * ts_yyyzz_xxyzz[i] * fe_0 + ts_yyyzz_xxyyzz[i] * pa_y[i];

        ts_yyyyzz_xxyzzz[i] = 3.0 * ts_yyzz_xxyzzz[i] * fe_0 + ts_yyyzz_xxzzz[i] * fe_0 + ts_yyyzz_xxyzzz[i] * pa_y[i];

        ts_yyyyzz_xxzzzz[i] = 3.0 * ts_yyzz_xxzzzz[i] * fe_0 + ts_yyyzz_xxzzzz[i] * pa_y[i];

        ts_yyyyzz_xyyyyy[i] = ts_yyyy_xyyyyy[i] * fe_0 + ts_yyyyz_xyyyyy[i] * pa_z[i];

        ts_yyyyzz_xyyyyz[i] = 3.0 * ts_yyzz_xyyyyz[i] * fe_0 + 4.0 * ts_yyyzz_xyyyz[i] * fe_0 + ts_yyyzz_xyyyyz[i] * pa_y[i];

        ts_yyyyzz_xyyyzz[i] = 3.0 * ts_yyzz_xyyyzz[i] * fe_0 + 3.0 * ts_yyyzz_xyyzz[i] * fe_0 + ts_yyyzz_xyyyzz[i] * pa_y[i];

        ts_yyyyzz_xyyzzz[i] = 3.0 * ts_yyzz_xyyzzz[i] * fe_0 + 2.0 * ts_yyyzz_xyzzz[i] * fe_0 + ts_yyyzz_xyyzzz[i] * pa_y[i];

        ts_yyyyzz_xyzzzz[i] = 3.0 * ts_yyzz_xyzzzz[i] * fe_0 + ts_yyyzz_xzzzz[i] * fe_0 + ts_yyyzz_xyzzzz[i] * pa_y[i];

        ts_yyyyzz_xzzzzz[i] = 3.0 * ts_yyzz_xzzzzz[i] * fe_0 + ts_yyyzz_xzzzzz[i] * pa_y[i];

        ts_yyyyzz_yyyyyy[i] = ts_yyyy_yyyyyy[i] * fe_0 + ts_yyyyz_yyyyyy[i] * pa_z[i];

        ts_yyyyzz_yyyyyz[i] = 3.0 * ts_yyzz_yyyyyz[i] * fe_0 + 5.0 * ts_yyyzz_yyyyz[i] * fe_0 + ts_yyyzz_yyyyyz[i] * pa_y[i];

        ts_yyyyzz_yyyyzz[i] = 3.0 * ts_yyzz_yyyyzz[i] * fe_0 + 4.0 * ts_yyyzz_yyyzz[i] * fe_0 + ts_yyyzz_yyyyzz[i] * pa_y[i];

        ts_yyyyzz_yyyzzz[i] = 3.0 * ts_yyzz_yyyzzz[i] * fe_0 + 3.0 * ts_yyyzz_yyzzz[i] * fe_0 + ts_yyyzz_yyyzzz[i] * pa_y[i];

        ts_yyyyzz_yyzzzz[i] = 3.0 * ts_yyzz_yyzzzz[i] * fe_0 + 2.0 * ts_yyyzz_yzzzz[i] * fe_0 + ts_yyyzz_yyzzzz[i] * pa_y[i];

        ts_yyyyzz_yzzzzz[i] = 3.0 * ts_yyzz_yzzzzz[i] * fe_0 + ts_yyyzz_zzzzz[i] * fe_0 + ts_yyyzz_yzzzzz[i] * pa_y[i];

        ts_yyyyzz_zzzzzz[i] = 3.0 * ts_yyzz_zzzzzz[i] * fe_0 + ts_yyyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 672-700 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyz_xxxxxy, ts_yyyz_xxxxyy, ts_yyyz_xxxyyy, ts_yyyz_xxyyyy, ts_yyyz_xyyyyy, ts_yyyz_yyyyyy, ts_yyyzz_xxxxxy, ts_yyyzz_xxxxyy, ts_yyyzz_xxxyyy, ts_yyyzz_xxyyyy, ts_yyyzz_xyyyyy, ts_yyyzz_yyyyyy, ts_yyyzzz_xxxxxx, ts_yyyzzz_xxxxxy, ts_yyyzzz_xxxxxz, ts_yyyzzz_xxxxyy, ts_yyyzzz_xxxxyz, ts_yyyzzz_xxxxzz, ts_yyyzzz_xxxyyy, ts_yyyzzz_xxxyyz, ts_yyyzzz_xxxyzz, ts_yyyzzz_xxxzzz, ts_yyyzzz_xxyyyy, ts_yyyzzz_xxyyyz, ts_yyyzzz_xxyyzz, ts_yyyzzz_xxyzzz, ts_yyyzzz_xxzzzz, ts_yyyzzz_xyyyyy, ts_yyyzzz_xyyyyz, ts_yyyzzz_xyyyzz, ts_yyyzzz_xyyzzz, ts_yyyzzz_xyzzzz, ts_yyyzzz_xzzzzz, ts_yyyzzz_yyyyyy, ts_yyyzzz_yyyyyz, ts_yyyzzz_yyyyzz, ts_yyyzzz_yyyzzz, ts_yyyzzz_yyzzzz, ts_yyyzzz_yzzzzz, ts_yyyzzz_zzzzzz, ts_yyzzz_xxxxxx, ts_yyzzz_xxxxxz, ts_yyzzz_xxxxyz, ts_yyzzz_xxxxz, ts_yyzzz_xxxxzz, ts_yyzzz_xxxyyz, ts_yyzzz_xxxyz, ts_yyzzz_xxxyzz, ts_yyzzz_xxxzz, ts_yyzzz_xxxzzz, ts_yyzzz_xxyyyz, ts_yyzzz_xxyyz, ts_yyzzz_xxyyzz, ts_yyzzz_xxyzz, ts_yyzzz_xxyzzz, ts_yyzzz_xxzzz, ts_yyzzz_xxzzzz, ts_yyzzz_xyyyyz, ts_yyzzz_xyyyz, ts_yyzzz_xyyyzz, ts_yyzzz_xyyzz, ts_yyzzz_xyyzzz, ts_yyzzz_xyzzz, ts_yyzzz_xyzzzz, ts_yyzzz_xzzzz, ts_yyzzz_xzzzzz, ts_yyzzz_yyyyyz, ts_yyzzz_yyyyz, ts_yyzzz_yyyyzz, ts_yyzzz_yyyzz, ts_yyzzz_yyyzzz, ts_yyzzz_yyzzz, ts_yyzzz_yyzzzz, ts_yyzzz_yzzzz, ts_yyzzz_yzzzzz, ts_yyzzz_zzzzz, ts_yyzzz_zzzzzz, ts_yzzz_xxxxxx, ts_yzzz_xxxxxz, ts_yzzz_xxxxyz, ts_yzzz_xxxxzz, ts_yzzz_xxxyyz, ts_yzzz_xxxyzz, ts_yzzz_xxxzzz, ts_yzzz_xxyyyz, ts_yzzz_xxyyzz, ts_yzzz_xxyzzz, ts_yzzz_xxzzzz, ts_yzzz_xyyyyz, ts_yzzz_xyyyzz, ts_yzzz_xyyzzz, ts_yzzz_xyzzzz, ts_yzzz_xzzzzz, ts_yzzz_yyyyyz, ts_yzzz_yyyyzz, ts_yzzz_yyyzzz, ts_yzzz_yyzzzz, ts_yzzz_yzzzzz, ts_yzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzzz_xxxxxx[i] = 2.0 * ts_yzzz_xxxxxx[i] * fe_0 + ts_yyzzz_xxxxxx[i] * pa_y[i];

        ts_yyyzzz_xxxxxy[i] = 2.0 * ts_yyyz_xxxxxy[i] * fe_0 + ts_yyyzz_xxxxxy[i] * pa_z[i];

        ts_yyyzzz_xxxxxz[i] = 2.0 * ts_yzzz_xxxxxz[i] * fe_0 + ts_yyzzz_xxxxxz[i] * pa_y[i];

        ts_yyyzzz_xxxxyy[i] = 2.0 * ts_yyyz_xxxxyy[i] * fe_0 + ts_yyyzz_xxxxyy[i] * pa_z[i];

        ts_yyyzzz_xxxxyz[i] = 2.0 * ts_yzzz_xxxxyz[i] * fe_0 + ts_yyzzz_xxxxz[i] * fe_0 + ts_yyzzz_xxxxyz[i] * pa_y[i];

        ts_yyyzzz_xxxxzz[i] = 2.0 * ts_yzzz_xxxxzz[i] * fe_0 + ts_yyzzz_xxxxzz[i] * pa_y[i];

        ts_yyyzzz_xxxyyy[i] = 2.0 * ts_yyyz_xxxyyy[i] * fe_0 + ts_yyyzz_xxxyyy[i] * pa_z[i];

        ts_yyyzzz_xxxyyz[i] = 2.0 * ts_yzzz_xxxyyz[i] * fe_0 + 2.0 * ts_yyzzz_xxxyz[i] * fe_0 + ts_yyzzz_xxxyyz[i] * pa_y[i];

        ts_yyyzzz_xxxyzz[i] = 2.0 * ts_yzzz_xxxyzz[i] * fe_0 + ts_yyzzz_xxxzz[i] * fe_0 + ts_yyzzz_xxxyzz[i] * pa_y[i];

        ts_yyyzzz_xxxzzz[i] = 2.0 * ts_yzzz_xxxzzz[i] * fe_0 + ts_yyzzz_xxxzzz[i] * pa_y[i];

        ts_yyyzzz_xxyyyy[i] = 2.0 * ts_yyyz_xxyyyy[i] * fe_0 + ts_yyyzz_xxyyyy[i] * pa_z[i];

        ts_yyyzzz_xxyyyz[i] = 2.0 * ts_yzzz_xxyyyz[i] * fe_0 + 3.0 * ts_yyzzz_xxyyz[i] * fe_0 + ts_yyzzz_xxyyyz[i] * pa_y[i];

        ts_yyyzzz_xxyyzz[i] = 2.0 * ts_yzzz_xxyyzz[i] * fe_0 + 2.0 * ts_yyzzz_xxyzz[i] * fe_0 + ts_yyzzz_xxyyzz[i] * pa_y[i];

        ts_yyyzzz_xxyzzz[i] = 2.0 * ts_yzzz_xxyzzz[i] * fe_0 + ts_yyzzz_xxzzz[i] * fe_0 + ts_yyzzz_xxyzzz[i] * pa_y[i];

        ts_yyyzzz_xxzzzz[i] = 2.0 * ts_yzzz_xxzzzz[i] * fe_0 + ts_yyzzz_xxzzzz[i] * pa_y[i];

        ts_yyyzzz_xyyyyy[i] = 2.0 * ts_yyyz_xyyyyy[i] * fe_0 + ts_yyyzz_xyyyyy[i] * pa_z[i];

        ts_yyyzzz_xyyyyz[i] = 2.0 * ts_yzzz_xyyyyz[i] * fe_0 + 4.0 * ts_yyzzz_xyyyz[i] * fe_0 + ts_yyzzz_xyyyyz[i] * pa_y[i];

        ts_yyyzzz_xyyyzz[i] = 2.0 * ts_yzzz_xyyyzz[i] * fe_0 + 3.0 * ts_yyzzz_xyyzz[i] * fe_0 + ts_yyzzz_xyyyzz[i] * pa_y[i];

        ts_yyyzzz_xyyzzz[i] = 2.0 * ts_yzzz_xyyzzz[i] * fe_0 + 2.0 * ts_yyzzz_xyzzz[i] * fe_0 + ts_yyzzz_xyyzzz[i] * pa_y[i];

        ts_yyyzzz_xyzzzz[i] = 2.0 * ts_yzzz_xyzzzz[i] * fe_0 + ts_yyzzz_xzzzz[i] * fe_0 + ts_yyzzz_xyzzzz[i] * pa_y[i];

        ts_yyyzzz_xzzzzz[i] = 2.0 * ts_yzzz_xzzzzz[i] * fe_0 + ts_yyzzz_xzzzzz[i] * pa_y[i];

        ts_yyyzzz_yyyyyy[i] = 2.0 * ts_yyyz_yyyyyy[i] * fe_0 + ts_yyyzz_yyyyyy[i] * pa_z[i];

        ts_yyyzzz_yyyyyz[i] = 2.0 * ts_yzzz_yyyyyz[i] * fe_0 + 5.0 * ts_yyzzz_yyyyz[i] * fe_0 + ts_yyzzz_yyyyyz[i] * pa_y[i];

        ts_yyyzzz_yyyyzz[i] = 2.0 * ts_yzzz_yyyyzz[i] * fe_0 + 4.0 * ts_yyzzz_yyyzz[i] * fe_0 + ts_yyzzz_yyyyzz[i] * pa_y[i];

        ts_yyyzzz_yyyzzz[i] = 2.0 * ts_yzzz_yyyzzz[i] * fe_0 + 3.0 * ts_yyzzz_yyzzz[i] * fe_0 + ts_yyzzz_yyyzzz[i] * pa_y[i];

        ts_yyyzzz_yyzzzz[i] = 2.0 * ts_yzzz_yyzzzz[i] * fe_0 + 2.0 * ts_yyzzz_yzzzz[i] * fe_0 + ts_yyzzz_yyzzzz[i] * pa_y[i];

        ts_yyyzzz_yzzzzz[i] = 2.0 * ts_yzzz_yzzzzz[i] * fe_0 + ts_yyzzz_zzzzz[i] * fe_0 + ts_yyzzz_yzzzzz[i] * pa_y[i];

        ts_yyyzzz_zzzzzz[i] = 2.0 * ts_yzzz_zzzzzz[i] * fe_0 + ts_yyzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 700-728 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_y, pa_z, ts_yyzz_xxxxxy, ts_yyzz_xxxxyy, ts_yyzz_xxxyyy, ts_yyzz_xxyyyy, ts_yyzz_xyyyyy, ts_yyzz_yyyyyy, ts_yyzzz_xxxxxy, ts_yyzzz_xxxxyy, ts_yyzzz_xxxyyy, ts_yyzzz_xxyyyy, ts_yyzzz_xyyyyy, ts_yyzzz_yyyyyy, ts_yyzzzz_xxxxxx, ts_yyzzzz_xxxxxy, ts_yyzzzz_xxxxxz, ts_yyzzzz_xxxxyy, ts_yyzzzz_xxxxyz, ts_yyzzzz_xxxxzz, ts_yyzzzz_xxxyyy, ts_yyzzzz_xxxyyz, ts_yyzzzz_xxxyzz, ts_yyzzzz_xxxzzz, ts_yyzzzz_xxyyyy, ts_yyzzzz_xxyyyz, ts_yyzzzz_xxyyzz, ts_yyzzzz_xxyzzz, ts_yyzzzz_xxzzzz, ts_yyzzzz_xyyyyy, ts_yyzzzz_xyyyyz, ts_yyzzzz_xyyyzz, ts_yyzzzz_xyyzzz, ts_yyzzzz_xyzzzz, ts_yyzzzz_xzzzzz, ts_yyzzzz_yyyyyy, ts_yyzzzz_yyyyyz, ts_yyzzzz_yyyyzz, ts_yyzzzz_yyyzzz, ts_yyzzzz_yyzzzz, ts_yyzzzz_yzzzzz, ts_yyzzzz_zzzzzz, ts_yzzzz_xxxxxx, ts_yzzzz_xxxxxz, ts_yzzzz_xxxxyz, ts_yzzzz_xxxxz, ts_yzzzz_xxxxzz, ts_yzzzz_xxxyyz, ts_yzzzz_xxxyz, ts_yzzzz_xxxyzz, ts_yzzzz_xxxzz, ts_yzzzz_xxxzzz, ts_yzzzz_xxyyyz, ts_yzzzz_xxyyz, ts_yzzzz_xxyyzz, ts_yzzzz_xxyzz, ts_yzzzz_xxyzzz, ts_yzzzz_xxzzz, ts_yzzzz_xxzzzz, ts_yzzzz_xyyyyz, ts_yzzzz_xyyyz, ts_yzzzz_xyyyzz, ts_yzzzz_xyyzz, ts_yzzzz_xyyzzz, ts_yzzzz_xyzzz, ts_yzzzz_xyzzzz, ts_yzzzz_xzzzz, ts_yzzzz_xzzzzz, ts_yzzzz_yyyyyz, ts_yzzzz_yyyyz, ts_yzzzz_yyyyzz, ts_yzzzz_yyyzz, ts_yzzzz_yyyzzz, ts_yzzzz_yyzzz, ts_yzzzz_yyzzzz, ts_yzzzz_yzzzz, ts_yzzzz_yzzzzz, ts_yzzzz_zzzzz, ts_yzzzz_zzzzzz, ts_zzzz_xxxxxx, ts_zzzz_xxxxxz, ts_zzzz_xxxxyz, ts_zzzz_xxxxzz, ts_zzzz_xxxyyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzzz_xxxxxx[i] = ts_zzzz_xxxxxx[i] * fe_0 + ts_yzzzz_xxxxxx[i] * pa_y[i];

        ts_yyzzzz_xxxxxy[i] = 3.0 * ts_yyzz_xxxxxy[i] * fe_0 + ts_yyzzz_xxxxxy[i] * pa_z[i];

        ts_yyzzzz_xxxxxz[i] = ts_zzzz_xxxxxz[i] * fe_0 + ts_yzzzz_xxxxxz[i] * pa_y[i];

        ts_yyzzzz_xxxxyy[i] = 3.0 * ts_yyzz_xxxxyy[i] * fe_0 + ts_yyzzz_xxxxyy[i] * pa_z[i];

        ts_yyzzzz_xxxxyz[i] = ts_zzzz_xxxxyz[i] * fe_0 + ts_yzzzz_xxxxz[i] * fe_0 + ts_yzzzz_xxxxyz[i] * pa_y[i];

        ts_yyzzzz_xxxxzz[i] = ts_zzzz_xxxxzz[i] * fe_0 + ts_yzzzz_xxxxzz[i] * pa_y[i];

        ts_yyzzzz_xxxyyy[i] = 3.0 * ts_yyzz_xxxyyy[i] * fe_0 + ts_yyzzz_xxxyyy[i] * pa_z[i];

        ts_yyzzzz_xxxyyz[i] = ts_zzzz_xxxyyz[i] * fe_0 + 2.0 * ts_yzzzz_xxxyz[i] * fe_0 + ts_yzzzz_xxxyyz[i] * pa_y[i];

        ts_yyzzzz_xxxyzz[i] = ts_zzzz_xxxyzz[i] * fe_0 + ts_yzzzz_xxxzz[i] * fe_0 + ts_yzzzz_xxxyzz[i] * pa_y[i];

        ts_yyzzzz_xxxzzz[i] = ts_zzzz_xxxzzz[i] * fe_0 + ts_yzzzz_xxxzzz[i] * pa_y[i];

        ts_yyzzzz_xxyyyy[i] = 3.0 * ts_yyzz_xxyyyy[i] * fe_0 + ts_yyzzz_xxyyyy[i] * pa_z[i];

        ts_yyzzzz_xxyyyz[i] = ts_zzzz_xxyyyz[i] * fe_0 + 3.0 * ts_yzzzz_xxyyz[i] * fe_0 + ts_yzzzz_xxyyyz[i] * pa_y[i];

        ts_yyzzzz_xxyyzz[i] = ts_zzzz_xxyyzz[i] * fe_0 + 2.0 * ts_yzzzz_xxyzz[i] * fe_0 + ts_yzzzz_xxyyzz[i] * pa_y[i];

        ts_yyzzzz_xxyzzz[i] = ts_zzzz_xxyzzz[i] * fe_0 + ts_yzzzz_xxzzz[i] * fe_0 + ts_yzzzz_xxyzzz[i] * pa_y[i];

        ts_yyzzzz_xxzzzz[i] = ts_zzzz_xxzzzz[i] * fe_0 + ts_yzzzz_xxzzzz[i] * pa_y[i];

        ts_yyzzzz_xyyyyy[i] = 3.0 * ts_yyzz_xyyyyy[i] * fe_0 + ts_yyzzz_xyyyyy[i] * pa_z[i];

        ts_yyzzzz_xyyyyz[i] = ts_zzzz_xyyyyz[i] * fe_0 + 4.0 * ts_yzzzz_xyyyz[i] * fe_0 + ts_yzzzz_xyyyyz[i] * pa_y[i];

        ts_yyzzzz_xyyyzz[i] = ts_zzzz_xyyyzz[i] * fe_0 + 3.0 * ts_yzzzz_xyyzz[i] * fe_0 + ts_yzzzz_xyyyzz[i] * pa_y[i];

        ts_yyzzzz_xyyzzz[i] = ts_zzzz_xyyzzz[i] * fe_0 + 2.0 * ts_yzzzz_xyzzz[i] * fe_0 + ts_yzzzz_xyyzzz[i] * pa_y[i];

        ts_yyzzzz_xyzzzz[i] = ts_zzzz_xyzzzz[i] * fe_0 + ts_yzzzz_xzzzz[i] * fe_0 + ts_yzzzz_xyzzzz[i] * pa_y[i];

        ts_yyzzzz_xzzzzz[i] = ts_zzzz_xzzzzz[i] * fe_0 + ts_yzzzz_xzzzzz[i] * pa_y[i];

        ts_yyzzzz_yyyyyy[i] = 3.0 * ts_yyzz_yyyyyy[i] * fe_0 + ts_yyzzz_yyyyyy[i] * pa_z[i];

        ts_yyzzzz_yyyyyz[i] = ts_zzzz_yyyyyz[i] * fe_0 + 5.0 * ts_yzzzz_yyyyz[i] * fe_0 + ts_yzzzz_yyyyyz[i] * pa_y[i];

        ts_yyzzzz_yyyyzz[i] = ts_zzzz_yyyyzz[i] * fe_0 + 4.0 * ts_yzzzz_yyyzz[i] * fe_0 + ts_yzzzz_yyyyzz[i] * pa_y[i];

        ts_yyzzzz_yyyzzz[i] = ts_zzzz_yyyzzz[i] * fe_0 + 3.0 * ts_yzzzz_yyzzz[i] * fe_0 + ts_yzzzz_yyyzzz[i] * pa_y[i];

        ts_yyzzzz_yyzzzz[i] = ts_zzzz_yyzzzz[i] * fe_0 + 2.0 * ts_yzzzz_yzzzz[i] * fe_0 + ts_yzzzz_yyzzzz[i] * pa_y[i];

        ts_yyzzzz_yzzzzz[i] = ts_zzzz_yzzzzz[i] * fe_0 + ts_yzzzz_zzzzz[i] * fe_0 + ts_yzzzz_yzzzzz[i] * pa_y[i];

        ts_yyzzzz_zzzzzz[i] = ts_zzzz_zzzzzz[i] * fe_0 + ts_yzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 728-756 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_y, ts_yzzzzz_xxxxxx, ts_yzzzzz_xxxxxy, ts_yzzzzz_xxxxxz, ts_yzzzzz_xxxxyy, ts_yzzzzz_xxxxyz, ts_yzzzzz_xxxxzz, ts_yzzzzz_xxxyyy, ts_yzzzzz_xxxyyz, ts_yzzzzz_xxxyzz, ts_yzzzzz_xxxzzz, ts_yzzzzz_xxyyyy, ts_yzzzzz_xxyyyz, ts_yzzzzz_xxyyzz, ts_yzzzzz_xxyzzz, ts_yzzzzz_xxzzzz, ts_yzzzzz_xyyyyy, ts_yzzzzz_xyyyyz, ts_yzzzzz_xyyyzz, ts_yzzzzz_xyyzzz, ts_yzzzzz_xyzzzz, ts_yzzzzz_xzzzzz, ts_yzzzzz_yyyyyy, ts_yzzzzz_yyyyyz, ts_yzzzzz_yyyyzz, ts_yzzzzz_yyyzzz, ts_yzzzzz_yyzzzz, ts_yzzzzz_yzzzzz, ts_yzzzzz_zzzzzz, ts_zzzzz_xxxxx, ts_zzzzz_xxxxxx, ts_zzzzz_xxxxxy, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxy, ts_zzzzz_xxxxyy, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyy, ts_zzzzz_xxxyyy, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyy, ts_zzzzz_xxyyyy, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyy, ts_zzzzz_xyyyyy, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyy, ts_zzzzz_yyyyyy, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzz, ts_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzzz_xxxxxx[i] = ts_zzzzz_xxxxxx[i] * pa_y[i];

        ts_yzzzzz_xxxxxy[i] = ts_zzzzz_xxxxx[i] * fe_0 + ts_zzzzz_xxxxxy[i] * pa_y[i];

        ts_yzzzzz_xxxxxz[i] = ts_zzzzz_xxxxxz[i] * pa_y[i];

        ts_yzzzzz_xxxxyy[i] = 2.0 * ts_zzzzz_xxxxy[i] * fe_0 + ts_zzzzz_xxxxyy[i] * pa_y[i];

        ts_yzzzzz_xxxxyz[i] = ts_zzzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxxyz[i] * pa_y[i];

        ts_yzzzzz_xxxxzz[i] = ts_zzzzz_xxxxzz[i] * pa_y[i];

        ts_yzzzzz_xxxyyy[i] = 3.0 * ts_zzzzz_xxxyy[i] * fe_0 + ts_zzzzz_xxxyyy[i] * pa_y[i];

        ts_yzzzzz_xxxyyz[i] = 2.0 * ts_zzzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxyyz[i] * pa_y[i];

        ts_yzzzzz_xxxyzz[i] = ts_zzzzz_xxxzz[i] * fe_0 + ts_zzzzz_xxxyzz[i] * pa_y[i];

        ts_yzzzzz_xxxzzz[i] = ts_zzzzz_xxxzzz[i] * pa_y[i];

        ts_yzzzzz_xxyyyy[i] = 4.0 * ts_zzzzz_xxyyy[i] * fe_0 + ts_zzzzz_xxyyyy[i] * pa_y[i];

        ts_yzzzzz_xxyyyz[i] = 3.0 * ts_zzzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxyyyz[i] * pa_y[i];

        ts_yzzzzz_xxyyzz[i] = 2.0 * ts_zzzzz_xxyzz[i] * fe_0 + ts_zzzzz_xxyyzz[i] * pa_y[i];

        ts_yzzzzz_xxyzzz[i] = ts_zzzzz_xxzzz[i] * fe_0 + ts_zzzzz_xxyzzz[i] * pa_y[i];

        ts_yzzzzz_xxzzzz[i] = ts_zzzzz_xxzzzz[i] * pa_y[i];

        ts_yzzzzz_xyyyyy[i] = 5.0 * ts_zzzzz_xyyyy[i] * fe_0 + ts_zzzzz_xyyyyy[i] * pa_y[i];

        ts_yzzzzz_xyyyyz[i] = 4.0 * ts_zzzzz_xyyyz[i] * fe_0 + ts_zzzzz_xyyyyz[i] * pa_y[i];

        ts_yzzzzz_xyyyzz[i] = 3.0 * ts_zzzzz_xyyzz[i] * fe_0 + ts_zzzzz_xyyyzz[i] * pa_y[i];

        ts_yzzzzz_xyyzzz[i] = 2.0 * ts_zzzzz_xyzzz[i] * fe_0 + ts_zzzzz_xyyzzz[i] * pa_y[i];

        ts_yzzzzz_xyzzzz[i] = ts_zzzzz_xzzzz[i] * fe_0 + ts_zzzzz_xyzzzz[i] * pa_y[i];

        ts_yzzzzz_xzzzzz[i] = ts_zzzzz_xzzzzz[i] * pa_y[i];

        ts_yzzzzz_yyyyyy[i] = 6.0 * ts_zzzzz_yyyyy[i] * fe_0 + ts_zzzzz_yyyyyy[i] * pa_y[i];

        ts_yzzzzz_yyyyyz[i] = 5.0 * ts_zzzzz_yyyyz[i] * fe_0 + ts_zzzzz_yyyyyz[i] * pa_y[i];

        ts_yzzzzz_yyyyzz[i] = 4.0 * ts_zzzzz_yyyzz[i] * fe_0 + ts_zzzzz_yyyyzz[i] * pa_y[i];

        ts_yzzzzz_yyyzzz[i] = 3.0 * ts_zzzzz_yyzzz[i] * fe_0 + ts_zzzzz_yyyzzz[i] * pa_y[i];

        ts_yzzzzz_yyzzzz[i] = 2.0 * ts_zzzzz_yzzzz[i] * fe_0 + ts_zzzzz_yyzzzz[i] * pa_y[i];

        ts_yzzzzz_yzzzzz[i] = ts_zzzzz_zzzzz[i] * fe_0 + ts_zzzzz_yzzzzz[i] * pa_y[i];

        ts_yzzzzz_zzzzzz[i] = ts_zzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 756-784 components of targeted buffer : II

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

    #pragma omp simd aligned(pa_z, ts_zzzz_xxxxxx, ts_zzzz_xxxxxy, ts_zzzz_xxxxxz, ts_zzzz_xxxxyy, ts_zzzz_xxxxyz, ts_zzzz_xxxxzz, ts_zzzz_xxxyyy, ts_zzzz_xxxyyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyyy, ts_zzzz_xxyyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyyy, ts_zzzz_xyyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyyy, ts_zzzz_yyyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzzz, ts_zzzzz_xxxxx, ts_zzzzz_xxxxxx, ts_zzzzz_xxxxxy, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxy, ts_zzzzz_xxxxyy, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyy, ts_zzzzz_xxxyyy, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyy, ts_zzzzz_xxyyyy, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyy, ts_zzzzz_xyyyyy, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyy, ts_zzzzz_yyyyyy, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzz, ts_zzzzz_zzzzzz, ts_zzzzzz_xxxxxx, ts_zzzzzz_xxxxxy, ts_zzzzzz_xxxxxz, ts_zzzzzz_xxxxyy, ts_zzzzzz_xxxxyz, ts_zzzzzz_xxxxzz, ts_zzzzzz_xxxyyy, ts_zzzzzz_xxxyyz, ts_zzzzzz_xxxyzz, ts_zzzzzz_xxxzzz, ts_zzzzzz_xxyyyy, ts_zzzzzz_xxyyyz, ts_zzzzzz_xxyyzz, ts_zzzzzz_xxyzzz, ts_zzzzzz_xxzzzz, ts_zzzzzz_xyyyyy, ts_zzzzzz_xyyyyz, ts_zzzzzz_xyyyzz, ts_zzzzzz_xyyzzz, ts_zzzzzz_xyzzzz, ts_zzzzzz_xzzzzz, ts_zzzzzz_yyyyyy, ts_zzzzzz_yyyyyz, ts_zzzzzz_yyyyzz, ts_zzzzzz_yyyzzz, ts_zzzzzz_yyzzzz, ts_zzzzzz_yzzzzz, ts_zzzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzzz_xxxxxx[i] = 5.0 * ts_zzzz_xxxxxx[i] * fe_0 + ts_zzzzz_xxxxxx[i] * pa_z[i];

        ts_zzzzzz_xxxxxy[i] = 5.0 * ts_zzzz_xxxxxy[i] * fe_0 + ts_zzzzz_xxxxxy[i] * pa_z[i];

        ts_zzzzzz_xxxxxz[i] = 5.0 * ts_zzzz_xxxxxz[i] * fe_0 + ts_zzzzz_xxxxx[i] * fe_0 + ts_zzzzz_xxxxxz[i] * pa_z[i];

        ts_zzzzzz_xxxxyy[i] = 5.0 * ts_zzzz_xxxxyy[i] * fe_0 + ts_zzzzz_xxxxyy[i] * pa_z[i];

        ts_zzzzzz_xxxxyz[i] = 5.0 * ts_zzzz_xxxxyz[i] * fe_0 + ts_zzzzz_xxxxy[i] * fe_0 + ts_zzzzz_xxxxyz[i] * pa_z[i];

        ts_zzzzzz_xxxxzz[i] = 5.0 * ts_zzzz_xxxxzz[i] * fe_0 + 2.0 * ts_zzzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxxzz[i] * pa_z[i];

        ts_zzzzzz_xxxyyy[i] = 5.0 * ts_zzzz_xxxyyy[i] * fe_0 + ts_zzzzz_xxxyyy[i] * pa_z[i];

        ts_zzzzzz_xxxyyz[i] = 5.0 * ts_zzzz_xxxyyz[i] * fe_0 + ts_zzzzz_xxxyy[i] * fe_0 + ts_zzzzz_xxxyyz[i] * pa_z[i];

        ts_zzzzzz_xxxyzz[i] = 5.0 * ts_zzzz_xxxyzz[i] * fe_0 + 2.0 * ts_zzzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxyzz[i] * pa_z[i];

        ts_zzzzzz_xxxzzz[i] = 5.0 * ts_zzzz_xxxzzz[i] * fe_0 + 3.0 * ts_zzzzz_xxxzz[i] * fe_0 + ts_zzzzz_xxxzzz[i] * pa_z[i];

        ts_zzzzzz_xxyyyy[i] = 5.0 * ts_zzzz_xxyyyy[i] * fe_0 + ts_zzzzz_xxyyyy[i] * pa_z[i];

        ts_zzzzzz_xxyyyz[i] = 5.0 * ts_zzzz_xxyyyz[i] * fe_0 + ts_zzzzz_xxyyy[i] * fe_0 + ts_zzzzz_xxyyyz[i] * pa_z[i];

        ts_zzzzzz_xxyyzz[i] = 5.0 * ts_zzzz_xxyyzz[i] * fe_0 + 2.0 * ts_zzzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxyyzz[i] * pa_z[i];

        ts_zzzzzz_xxyzzz[i] = 5.0 * ts_zzzz_xxyzzz[i] * fe_0 + 3.0 * ts_zzzzz_xxyzz[i] * fe_0 + ts_zzzzz_xxyzzz[i] * pa_z[i];

        ts_zzzzzz_xxzzzz[i] = 5.0 * ts_zzzz_xxzzzz[i] * fe_0 + 4.0 * ts_zzzzz_xxzzz[i] * fe_0 + ts_zzzzz_xxzzzz[i] * pa_z[i];

        ts_zzzzzz_xyyyyy[i] = 5.0 * ts_zzzz_xyyyyy[i] * fe_0 + ts_zzzzz_xyyyyy[i] * pa_z[i];

        ts_zzzzzz_xyyyyz[i] = 5.0 * ts_zzzz_xyyyyz[i] * fe_0 + ts_zzzzz_xyyyy[i] * fe_0 + ts_zzzzz_xyyyyz[i] * pa_z[i];

        ts_zzzzzz_xyyyzz[i] = 5.0 * ts_zzzz_xyyyzz[i] * fe_0 + 2.0 * ts_zzzzz_xyyyz[i] * fe_0 + ts_zzzzz_xyyyzz[i] * pa_z[i];

        ts_zzzzzz_xyyzzz[i] = 5.0 * ts_zzzz_xyyzzz[i] * fe_0 + 3.0 * ts_zzzzz_xyyzz[i] * fe_0 + ts_zzzzz_xyyzzz[i] * pa_z[i];

        ts_zzzzzz_xyzzzz[i] = 5.0 * ts_zzzz_xyzzzz[i] * fe_0 + 4.0 * ts_zzzzz_xyzzz[i] * fe_0 + ts_zzzzz_xyzzzz[i] * pa_z[i];

        ts_zzzzzz_xzzzzz[i] = 5.0 * ts_zzzz_xzzzzz[i] * fe_0 + 5.0 * ts_zzzzz_xzzzz[i] * fe_0 + ts_zzzzz_xzzzzz[i] * pa_z[i];

        ts_zzzzzz_yyyyyy[i] = 5.0 * ts_zzzz_yyyyyy[i] * fe_0 + ts_zzzzz_yyyyyy[i] * pa_z[i];

        ts_zzzzzz_yyyyyz[i] = 5.0 * ts_zzzz_yyyyyz[i] * fe_0 + ts_zzzzz_yyyyy[i] * fe_0 + ts_zzzzz_yyyyyz[i] * pa_z[i];

        ts_zzzzzz_yyyyzz[i] = 5.0 * ts_zzzz_yyyyzz[i] * fe_0 + 2.0 * ts_zzzzz_yyyyz[i] * fe_0 + ts_zzzzz_yyyyzz[i] * pa_z[i];

        ts_zzzzzz_yyyzzz[i] = 5.0 * ts_zzzz_yyyzzz[i] * fe_0 + 3.0 * ts_zzzzz_yyyzz[i] * fe_0 + ts_zzzzz_yyyzzz[i] * pa_z[i];

        ts_zzzzzz_yyzzzz[i] = 5.0 * ts_zzzz_yyzzzz[i] * fe_0 + 4.0 * ts_zzzzz_yyzzz[i] * fe_0 + ts_zzzzz_yyzzzz[i] * pa_z[i];

        ts_zzzzzz_yzzzzz[i] = 5.0 * ts_zzzz_yzzzzz[i] * fe_0 + 5.0 * ts_zzzzz_yzzzz[i] * fe_0 + ts_zzzzz_yzzzzz[i] * pa_z[i];

        ts_zzzzzz_zzzzzz[i] = 5.0 * ts_zzzz_zzzzzz[i] * fe_0 + 6.0 * ts_zzzzz_zzzzz[i] * fe_0 + ts_zzzzz_zzzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

