#include "ElectricDipoleMomentumPrimRecII.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_ii(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_ii,
                                      const size_t idx_dip_gi,
                                      const size_t idx_dip_hh,
                                      const size_t idx_ovl_hi,
                                      const size_t idx_dip_hi,
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

    auto tr_x_xxxx_xxxxxx = pbuffer.data(idx_dip_gi);

    auto tr_x_xxxx_xxxxxy = pbuffer.data(idx_dip_gi + 1);

    auto tr_x_xxxx_xxxxxz = pbuffer.data(idx_dip_gi + 2);

    auto tr_x_xxxx_xxxxyy = pbuffer.data(idx_dip_gi + 3);

    auto tr_x_xxxx_xxxxyz = pbuffer.data(idx_dip_gi + 4);

    auto tr_x_xxxx_xxxxzz = pbuffer.data(idx_dip_gi + 5);

    auto tr_x_xxxx_xxxyyy = pbuffer.data(idx_dip_gi + 6);

    auto tr_x_xxxx_xxxyyz = pbuffer.data(idx_dip_gi + 7);

    auto tr_x_xxxx_xxxyzz = pbuffer.data(idx_dip_gi + 8);

    auto tr_x_xxxx_xxxzzz = pbuffer.data(idx_dip_gi + 9);

    auto tr_x_xxxx_xxyyyy = pbuffer.data(idx_dip_gi + 10);

    auto tr_x_xxxx_xxyyyz = pbuffer.data(idx_dip_gi + 11);

    auto tr_x_xxxx_xxyyzz = pbuffer.data(idx_dip_gi + 12);

    auto tr_x_xxxx_xxyzzz = pbuffer.data(idx_dip_gi + 13);

    auto tr_x_xxxx_xxzzzz = pbuffer.data(idx_dip_gi + 14);

    auto tr_x_xxxx_xyyyyy = pbuffer.data(idx_dip_gi + 15);

    auto tr_x_xxxx_xyyyyz = pbuffer.data(idx_dip_gi + 16);

    auto tr_x_xxxx_xyyyzz = pbuffer.data(idx_dip_gi + 17);

    auto tr_x_xxxx_xyyzzz = pbuffer.data(idx_dip_gi + 18);

    auto tr_x_xxxx_xyzzzz = pbuffer.data(idx_dip_gi + 19);

    auto tr_x_xxxx_xzzzzz = pbuffer.data(idx_dip_gi + 20);

    auto tr_x_xxxx_yyyyyy = pbuffer.data(idx_dip_gi + 21);

    auto tr_x_xxxx_yyyyyz = pbuffer.data(idx_dip_gi + 22);

    auto tr_x_xxxx_yyyyzz = pbuffer.data(idx_dip_gi + 23);

    auto tr_x_xxxx_yyyzzz = pbuffer.data(idx_dip_gi + 24);

    auto tr_x_xxxx_yyzzzz = pbuffer.data(idx_dip_gi + 25);

    auto tr_x_xxxx_yzzzzz = pbuffer.data(idx_dip_gi + 26);

    auto tr_x_xxxx_zzzzzz = pbuffer.data(idx_dip_gi + 27);

    auto tr_x_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 28);

    auto tr_x_xxxy_xxxxxy = pbuffer.data(idx_dip_gi + 29);

    auto tr_x_xxxy_xxxxxz = pbuffer.data(idx_dip_gi + 30);

    auto tr_x_xxxy_xxxxyy = pbuffer.data(idx_dip_gi + 31);

    auto tr_x_xxxy_xxxxyz = pbuffer.data(idx_dip_gi + 32);

    auto tr_x_xxxy_xxxxzz = pbuffer.data(idx_dip_gi + 33);

    auto tr_x_xxxy_xxxyyy = pbuffer.data(idx_dip_gi + 34);

    auto tr_x_xxxy_xxxyyz = pbuffer.data(idx_dip_gi + 35);

    auto tr_x_xxxy_xxxyzz = pbuffer.data(idx_dip_gi + 36);

    auto tr_x_xxxy_xxxzzz = pbuffer.data(idx_dip_gi + 37);

    auto tr_x_xxxy_xxyyyy = pbuffer.data(idx_dip_gi + 38);

    auto tr_x_xxxy_xxyyyz = pbuffer.data(idx_dip_gi + 39);

    auto tr_x_xxxy_xxyyzz = pbuffer.data(idx_dip_gi + 40);

    auto tr_x_xxxy_xxyzzz = pbuffer.data(idx_dip_gi + 41);

    auto tr_x_xxxy_xxzzzz = pbuffer.data(idx_dip_gi + 42);

    auto tr_x_xxxy_xyyyyy = pbuffer.data(idx_dip_gi + 43);

    auto tr_x_xxxy_xyyyyz = pbuffer.data(idx_dip_gi + 44);

    auto tr_x_xxxy_xyyyzz = pbuffer.data(idx_dip_gi + 45);

    auto tr_x_xxxy_xyyzzz = pbuffer.data(idx_dip_gi + 46);

    auto tr_x_xxxy_xyzzzz = pbuffer.data(idx_dip_gi + 47);

    auto tr_x_xxxy_xzzzzz = pbuffer.data(idx_dip_gi + 48);

    auto tr_x_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 55);

    auto tr_x_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 56);

    auto tr_x_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 57);

    auto tr_x_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 58);

    auto tr_x_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 59);

    auto tr_x_xxxz_xxxxyz = pbuffer.data(idx_dip_gi + 60);

    auto tr_x_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 61);

    auto tr_x_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 62);

    auto tr_x_xxxz_xxxyyz = pbuffer.data(idx_dip_gi + 63);

    auto tr_x_xxxz_xxxyzz = pbuffer.data(idx_dip_gi + 64);

    auto tr_x_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 65);

    auto tr_x_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 66);

    auto tr_x_xxxz_xxyyyz = pbuffer.data(idx_dip_gi + 67);

    auto tr_x_xxxz_xxyyzz = pbuffer.data(idx_dip_gi + 68);

    auto tr_x_xxxz_xxyzzz = pbuffer.data(idx_dip_gi + 69);

    auto tr_x_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 70);

    auto tr_x_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 71);

    auto tr_x_xxxz_xyyyyz = pbuffer.data(idx_dip_gi + 72);

    auto tr_x_xxxz_xyyyzz = pbuffer.data(idx_dip_gi + 73);

    auto tr_x_xxxz_xyyzzz = pbuffer.data(idx_dip_gi + 74);

    auto tr_x_xxxz_xyzzzz = pbuffer.data(idx_dip_gi + 75);

    auto tr_x_xxxz_xzzzzz = pbuffer.data(idx_dip_gi + 76);

    auto tr_x_xxxz_yyyyyy = pbuffer.data(idx_dip_gi + 77);

    auto tr_x_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 83);

    auto tr_x_xxyy_xxxxxx = pbuffer.data(idx_dip_gi + 84);

    auto tr_x_xxyy_xxxxxy = pbuffer.data(idx_dip_gi + 85);

    auto tr_x_xxyy_xxxxxz = pbuffer.data(idx_dip_gi + 86);

    auto tr_x_xxyy_xxxxyy = pbuffer.data(idx_dip_gi + 87);

    auto tr_x_xxyy_xxxxyz = pbuffer.data(idx_dip_gi + 88);

    auto tr_x_xxyy_xxxxzz = pbuffer.data(idx_dip_gi + 89);

    auto tr_x_xxyy_xxxyyy = pbuffer.data(idx_dip_gi + 90);

    auto tr_x_xxyy_xxxyyz = pbuffer.data(idx_dip_gi + 91);

    auto tr_x_xxyy_xxxyzz = pbuffer.data(idx_dip_gi + 92);

    auto tr_x_xxyy_xxxzzz = pbuffer.data(idx_dip_gi + 93);

    auto tr_x_xxyy_xxyyyy = pbuffer.data(idx_dip_gi + 94);

    auto tr_x_xxyy_xxyyyz = pbuffer.data(idx_dip_gi + 95);

    auto tr_x_xxyy_xxyyzz = pbuffer.data(idx_dip_gi + 96);

    auto tr_x_xxyy_xxyzzz = pbuffer.data(idx_dip_gi + 97);

    auto tr_x_xxyy_xxzzzz = pbuffer.data(idx_dip_gi + 98);

    auto tr_x_xxyy_xyyyyy = pbuffer.data(idx_dip_gi + 99);

    auto tr_x_xxyy_xyyyyz = pbuffer.data(idx_dip_gi + 100);

    auto tr_x_xxyy_xyyyzz = pbuffer.data(idx_dip_gi + 101);

    auto tr_x_xxyy_xyyzzz = pbuffer.data(idx_dip_gi + 102);

    auto tr_x_xxyy_xyzzzz = pbuffer.data(idx_dip_gi + 103);

    auto tr_x_xxyy_xzzzzz = pbuffer.data(idx_dip_gi + 104);

    auto tr_x_xxyy_yyyyyy = pbuffer.data(idx_dip_gi + 105);

    auto tr_x_xxyy_yyyyyz = pbuffer.data(idx_dip_gi + 106);

    auto tr_x_xxyy_yyyyzz = pbuffer.data(idx_dip_gi + 107);

    auto tr_x_xxyy_yyyzzz = pbuffer.data(idx_dip_gi + 108);

    auto tr_x_xxyy_yyzzzz = pbuffer.data(idx_dip_gi + 109);

    auto tr_x_xxyy_yzzzzz = pbuffer.data(idx_dip_gi + 110);

    auto tr_x_xxyy_zzzzzz = pbuffer.data(idx_dip_gi + 111);

    auto tr_x_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 114);

    auto tr_x_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 117);

    auto tr_x_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 121);

    auto tr_x_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 126);

    auto tr_x_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 132);

    auto tr_x_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 139);

    auto tr_x_xxzz_xxxxxx = pbuffer.data(idx_dip_gi + 140);

    auto tr_x_xxzz_xxxxxy = pbuffer.data(idx_dip_gi + 141);

    auto tr_x_xxzz_xxxxxz = pbuffer.data(idx_dip_gi + 142);

    auto tr_x_xxzz_xxxxyy = pbuffer.data(idx_dip_gi + 143);

    auto tr_x_xxzz_xxxxyz = pbuffer.data(idx_dip_gi + 144);

    auto tr_x_xxzz_xxxxzz = pbuffer.data(idx_dip_gi + 145);

    auto tr_x_xxzz_xxxyyy = pbuffer.data(idx_dip_gi + 146);

    auto tr_x_xxzz_xxxyyz = pbuffer.data(idx_dip_gi + 147);

    auto tr_x_xxzz_xxxyzz = pbuffer.data(idx_dip_gi + 148);

    auto tr_x_xxzz_xxxzzz = pbuffer.data(idx_dip_gi + 149);

    auto tr_x_xxzz_xxyyyy = pbuffer.data(idx_dip_gi + 150);

    auto tr_x_xxzz_xxyyyz = pbuffer.data(idx_dip_gi + 151);

    auto tr_x_xxzz_xxyyzz = pbuffer.data(idx_dip_gi + 152);

    auto tr_x_xxzz_xxyzzz = pbuffer.data(idx_dip_gi + 153);

    auto tr_x_xxzz_xxzzzz = pbuffer.data(idx_dip_gi + 154);

    auto tr_x_xxzz_xyyyyy = pbuffer.data(idx_dip_gi + 155);

    auto tr_x_xxzz_xyyyyz = pbuffer.data(idx_dip_gi + 156);

    auto tr_x_xxzz_xyyyzz = pbuffer.data(idx_dip_gi + 157);

    auto tr_x_xxzz_xyyzzz = pbuffer.data(idx_dip_gi + 158);

    auto tr_x_xxzz_xyzzzz = pbuffer.data(idx_dip_gi + 159);

    auto tr_x_xxzz_xzzzzz = pbuffer.data(idx_dip_gi + 160);

    auto tr_x_xxzz_yyyyyy = pbuffer.data(idx_dip_gi + 161);

    auto tr_x_xxzz_yyyyyz = pbuffer.data(idx_dip_gi + 162);

    auto tr_x_xxzz_yyyyzz = pbuffer.data(idx_dip_gi + 163);

    auto tr_x_xxzz_yyyzzz = pbuffer.data(idx_dip_gi + 164);

    auto tr_x_xxzz_yyzzzz = pbuffer.data(idx_dip_gi + 165);

    auto tr_x_xxzz_yzzzzz = pbuffer.data(idx_dip_gi + 166);

    auto tr_x_xxzz_zzzzzz = pbuffer.data(idx_dip_gi + 167);

    auto tr_x_xyyy_xxxxxx = pbuffer.data(idx_dip_gi + 168);

    auto tr_x_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 169);

    auto tr_x_xyyy_xxxxxz = pbuffer.data(idx_dip_gi + 170);

    auto tr_x_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 171);

    auto tr_x_xyyy_xxxxzz = pbuffer.data(idx_dip_gi + 173);

    auto tr_x_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 174);

    auto tr_x_xyyy_xxxzzz = pbuffer.data(idx_dip_gi + 177);

    auto tr_x_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 178);

    auto tr_x_xyyy_xxzzzz = pbuffer.data(idx_dip_gi + 182);

    auto tr_x_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 183);

    auto tr_x_xyyy_xzzzzz = pbuffer.data(idx_dip_gi + 188);

    auto tr_x_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 189);

    auto tr_x_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 190);

    auto tr_x_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 191);

    auto tr_x_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 192);

    auto tr_x_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 193);

    auto tr_x_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 194);

    auto tr_x_xyyz_xxxxxy = pbuffer.data(idx_dip_gi + 197);

    auto tr_x_xyyz_xxxxxz = pbuffer.data(idx_dip_gi + 198);

    auto tr_x_xyyz_xxxxyy = pbuffer.data(idx_dip_gi + 199);

    auto tr_x_xyyz_xxxxzz = pbuffer.data(idx_dip_gi + 201);

    auto tr_x_xyyz_xxxyyy = pbuffer.data(idx_dip_gi + 202);

    auto tr_x_xyyz_xxxzzz = pbuffer.data(idx_dip_gi + 205);

    auto tr_x_xyyz_xxyyyy = pbuffer.data(idx_dip_gi + 206);

    auto tr_x_xyyz_xxzzzz = pbuffer.data(idx_dip_gi + 210);

    auto tr_x_xyyz_xyyyyy = pbuffer.data(idx_dip_gi + 211);

    auto tr_x_xyyz_xzzzzz = pbuffer.data(idx_dip_gi + 216);

    auto tr_x_xyzz_xxxxxx = pbuffer.data(idx_dip_gi + 224);

    auto tr_x_xyzz_xxxxxz = pbuffer.data(idx_dip_gi + 226);

    auto tr_x_xyzz_xxxxzz = pbuffer.data(idx_dip_gi + 229);

    auto tr_x_xyzz_xxxzzz = pbuffer.data(idx_dip_gi + 233);

    auto tr_x_xyzz_xxzzzz = pbuffer.data(idx_dip_gi + 238);

    auto tr_x_xyzz_xzzzzz = pbuffer.data(idx_dip_gi + 244);

    auto tr_x_xzzz_xxxxxx = pbuffer.data(idx_dip_gi + 252);

    auto tr_x_xzzz_xxxxxy = pbuffer.data(idx_dip_gi + 253);

    auto tr_x_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 254);

    auto tr_x_xzzz_xxxxyy = pbuffer.data(idx_dip_gi + 255);

    auto tr_x_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 257);

    auto tr_x_xzzz_xxxyyy = pbuffer.data(idx_dip_gi + 258);

    auto tr_x_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 261);

    auto tr_x_xzzz_xxyyyy = pbuffer.data(idx_dip_gi + 262);

    auto tr_x_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 266);

    auto tr_x_xzzz_xyyyyy = pbuffer.data(idx_dip_gi + 267);

    auto tr_x_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 272);

    auto tr_x_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 274);

    auto tr_x_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 275);

    auto tr_x_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 276);

    auto tr_x_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 277);

    auto tr_x_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 278);

    auto tr_x_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 279);

    auto tr_x_yyyy_xxxxxx = pbuffer.data(idx_dip_gi + 280);

    auto tr_x_yyyy_xxxxxy = pbuffer.data(idx_dip_gi + 281);

    auto tr_x_yyyy_xxxxxz = pbuffer.data(idx_dip_gi + 282);

    auto tr_x_yyyy_xxxxyy = pbuffer.data(idx_dip_gi + 283);

    auto tr_x_yyyy_xxxxyz = pbuffer.data(idx_dip_gi + 284);

    auto tr_x_yyyy_xxxxzz = pbuffer.data(idx_dip_gi + 285);

    auto tr_x_yyyy_xxxyyy = pbuffer.data(idx_dip_gi + 286);

    auto tr_x_yyyy_xxxyyz = pbuffer.data(idx_dip_gi + 287);

    auto tr_x_yyyy_xxxyzz = pbuffer.data(idx_dip_gi + 288);

    auto tr_x_yyyy_xxxzzz = pbuffer.data(idx_dip_gi + 289);

    auto tr_x_yyyy_xxyyyy = pbuffer.data(idx_dip_gi + 290);

    auto tr_x_yyyy_xxyyyz = pbuffer.data(idx_dip_gi + 291);

    auto tr_x_yyyy_xxyyzz = pbuffer.data(idx_dip_gi + 292);

    auto tr_x_yyyy_xxyzzz = pbuffer.data(idx_dip_gi + 293);

    auto tr_x_yyyy_xxzzzz = pbuffer.data(idx_dip_gi + 294);

    auto tr_x_yyyy_xyyyyy = pbuffer.data(idx_dip_gi + 295);

    auto tr_x_yyyy_xyyyyz = pbuffer.data(idx_dip_gi + 296);

    auto tr_x_yyyy_xyyyzz = pbuffer.data(idx_dip_gi + 297);

    auto tr_x_yyyy_xyyzzz = pbuffer.data(idx_dip_gi + 298);

    auto tr_x_yyyy_xyzzzz = pbuffer.data(idx_dip_gi + 299);

    auto tr_x_yyyy_xzzzzz = pbuffer.data(idx_dip_gi + 300);

    auto tr_x_yyyy_yyyyyy = pbuffer.data(idx_dip_gi + 301);

    auto tr_x_yyyy_yyyyyz = pbuffer.data(idx_dip_gi + 302);

    auto tr_x_yyyy_yyyyzz = pbuffer.data(idx_dip_gi + 303);

    auto tr_x_yyyy_yyyzzz = pbuffer.data(idx_dip_gi + 304);

    auto tr_x_yyyy_yyzzzz = pbuffer.data(idx_dip_gi + 305);

    auto tr_x_yyyy_yzzzzz = pbuffer.data(idx_dip_gi + 306);

    auto tr_x_yyyy_zzzzzz = pbuffer.data(idx_dip_gi + 307);

    auto tr_x_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 309);

    auto tr_x_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 310);

    auto tr_x_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 311);

    auto tr_x_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 313);

    auto tr_x_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 314);

    auto tr_x_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 317);

    auto tr_x_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 318);

    auto tr_x_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 322);

    auto tr_x_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 323);

    auto tr_x_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 328);

    auto tr_x_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 329);

    auto tr_x_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 335);

    auto tr_x_yyzz_xxxxxx = pbuffer.data(idx_dip_gi + 336);

    auto tr_x_yyzz_xxxxxy = pbuffer.data(idx_dip_gi + 337);

    auto tr_x_yyzz_xxxxxz = pbuffer.data(idx_dip_gi + 338);

    auto tr_x_yyzz_xxxxyy = pbuffer.data(idx_dip_gi + 339);

    auto tr_x_yyzz_xxxxyz = pbuffer.data(idx_dip_gi + 340);

    auto tr_x_yyzz_xxxxzz = pbuffer.data(idx_dip_gi + 341);

    auto tr_x_yyzz_xxxyyy = pbuffer.data(idx_dip_gi + 342);

    auto tr_x_yyzz_xxxyyz = pbuffer.data(idx_dip_gi + 343);

    auto tr_x_yyzz_xxxyzz = pbuffer.data(idx_dip_gi + 344);

    auto tr_x_yyzz_xxxzzz = pbuffer.data(idx_dip_gi + 345);

    auto tr_x_yyzz_xxyyyy = pbuffer.data(idx_dip_gi + 346);

    auto tr_x_yyzz_xxyyyz = pbuffer.data(idx_dip_gi + 347);

    auto tr_x_yyzz_xxyyzz = pbuffer.data(idx_dip_gi + 348);

    auto tr_x_yyzz_xxyzzz = pbuffer.data(idx_dip_gi + 349);

    auto tr_x_yyzz_xxzzzz = pbuffer.data(idx_dip_gi + 350);

    auto tr_x_yyzz_xyyyyy = pbuffer.data(idx_dip_gi + 351);

    auto tr_x_yyzz_xyyyyz = pbuffer.data(idx_dip_gi + 352);

    auto tr_x_yyzz_xyyyzz = pbuffer.data(idx_dip_gi + 353);

    auto tr_x_yyzz_xyyzzz = pbuffer.data(idx_dip_gi + 354);

    auto tr_x_yyzz_xyzzzz = pbuffer.data(idx_dip_gi + 355);

    auto tr_x_yyzz_xzzzzz = pbuffer.data(idx_dip_gi + 356);

    auto tr_x_yyzz_yyyyyy = pbuffer.data(idx_dip_gi + 357);

    auto tr_x_yyzz_yyyyyz = pbuffer.data(idx_dip_gi + 358);

    auto tr_x_yyzz_yyyyzz = pbuffer.data(idx_dip_gi + 359);

    auto tr_x_yyzz_yyyzzz = pbuffer.data(idx_dip_gi + 360);

    auto tr_x_yyzz_yyzzzz = pbuffer.data(idx_dip_gi + 361);

    auto tr_x_yyzz_yzzzzz = pbuffer.data(idx_dip_gi + 362);

    auto tr_x_yyzz_zzzzzz = pbuffer.data(idx_dip_gi + 363);

    auto tr_x_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 364);

    auto tr_x_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 366);

    auto tr_x_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 368);

    auto tr_x_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 369);

    auto tr_x_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 371);

    auto tr_x_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 372);

    auto tr_x_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 373);

    auto tr_x_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 375);

    auto tr_x_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 376);

    auto tr_x_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 377);

    auto tr_x_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 378);

    auto tr_x_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 380);

    auto tr_x_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 381);

    auto tr_x_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 382);

    auto tr_x_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 383);

    auto tr_x_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 384);

    auto tr_x_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 386);

    auto tr_x_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 387);

    auto tr_x_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 388);

    auto tr_x_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 389);

    auto tr_x_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 390);

    auto tr_x_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 391);

    auto tr_x_zzzz_xxxxxx = pbuffer.data(idx_dip_gi + 392);

    auto tr_x_zzzz_xxxxxy = pbuffer.data(idx_dip_gi + 393);

    auto tr_x_zzzz_xxxxxz = pbuffer.data(idx_dip_gi + 394);

    auto tr_x_zzzz_xxxxyy = pbuffer.data(idx_dip_gi + 395);

    auto tr_x_zzzz_xxxxyz = pbuffer.data(idx_dip_gi + 396);

    auto tr_x_zzzz_xxxxzz = pbuffer.data(idx_dip_gi + 397);

    auto tr_x_zzzz_xxxyyy = pbuffer.data(idx_dip_gi + 398);

    auto tr_x_zzzz_xxxyyz = pbuffer.data(idx_dip_gi + 399);

    auto tr_x_zzzz_xxxyzz = pbuffer.data(idx_dip_gi + 400);

    auto tr_x_zzzz_xxxzzz = pbuffer.data(idx_dip_gi + 401);

    auto tr_x_zzzz_xxyyyy = pbuffer.data(idx_dip_gi + 402);

    auto tr_x_zzzz_xxyyyz = pbuffer.data(idx_dip_gi + 403);

    auto tr_x_zzzz_xxyyzz = pbuffer.data(idx_dip_gi + 404);

    auto tr_x_zzzz_xxyzzz = pbuffer.data(idx_dip_gi + 405);

    auto tr_x_zzzz_xxzzzz = pbuffer.data(idx_dip_gi + 406);

    auto tr_x_zzzz_xyyyyy = pbuffer.data(idx_dip_gi + 407);

    auto tr_x_zzzz_xyyyyz = pbuffer.data(idx_dip_gi + 408);

    auto tr_x_zzzz_xyyyzz = pbuffer.data(idx_dip_gi + 409);

    auto tr_x_zzzz_xyyzzz = pbuffer.data(idx_dip_gi + 410);

    auto tr_x_zzzz_xyzzzz = pbuffer.data(idx_dip_gi + 411);

    auto tr_x_zzzz_xzzzzz = pbuffer.data(idx_dip_gi + 412);

    auto tr_x_zzzz_yyyyyy = pbuffer.data(idx_dip_gi + 413);

    auto tr_x_zzzz_yyyyyz = pbuffer.data(idx_dip_gi + 414);

    auto tr_x_zzzz_yyyyzz = pbuffer.data(idx_dip_gi + 415);

    auto tr_x_zzzz_yyyzzz = pbuffer.data(idx_dip_gi + 416);

    auto tr_x_zzzz_yyzzzz = pbuffer.data(idx_dip_gi + 417);

    auto tr_x_zzzz_yzzzzz = pbuffer.data(idx_dip_gi + 418);

    auto tr_x_zzzz_zzzzzz = pbuffer.data(idx_dip_gi + 419);

    auto tr_y_xxxx_xxxxxx = pbuffer.data(idx_dip_gi + 420);

    auto tr_y_xxxx_xxxxxy = pbuffer.data(idx_dip_gi + 421);

    auto tr_y_xxxx_xxxxxz = pbuffer.data(idx_dip_gi + 422);

    auto tr_y_xxxx_xxxxyy = pbuffer.data(idx_dip_gi + 423);

    auto tr_y_xxxx_xxxxyz = pbuffer.data(idx_dip_gi + 424);

    auto tr_y_xxxx_xxxxzz = pbuffer.data(idx_dip_gi + 425);

    auto tr_y_xxxx_xxxyyy = pbuffer.data(idx_dip_gi + 426);

    auto tr_y_xxxx_xxxyyz = pbuffer.data(idx_dip_gi + 427);

    auto tr_y_xxxx_xxxyzz = pbuffer.data(idx_dip_gi + 428);

    auto tr_y_xxxx_xxxzzz = pbuffer.data(idx_dip_gi + 429);

    auto tr_y_xxxx_xxyyyy = pbuffer.data(idx_dip_gi + 430);

    auto tr_y_xxxx_xxyyyz = pbuffer.data(idx_dip_gi + 431);

    auto tr_y_xxxx_xxyyzz = pbuffer.data(idx_dip_gi + 432);

    auto tr_y_xxxx_xxyzzz = pbuffer.data(idx_dip_gi + 433);

    auto tr_y_xxxx_xxzzzz = pbuffer.data(idx_dip_gi + 434);

    auto tr_y_xxxx_xyyyyy = pbuffer.data(idx_dip_gi + 435);

    auto tr_y_xxxx_xyyyyz = pbuffer.data(idx_dip_gi + 436);

    auto tr_y_xxxx_xyyyzz = pbuffer.data(idx_dip_gi + 437);

    auto tr_y_xxxx_xyyzzz = pbuffer.data(idx_dip_gi + 438);

    auto tr_y_xxxx_xyzzzz = pbuffer.data(idx_dip_gi + 439);

    auto tr_y_xxxx_xzzzzz = pbuffer.data(idx_dip_gi + 440);

    auto tr_y_xxxx_yyyyyy = pbuffer.data(idx_dip_gi + 441);

    auto tr_y_xxxx_yyyyyz = pbuffer.data(idx_dip_gi + 442);

    auto tr_y_xxxx_yyyyzz = pbuffer.data(idx_dip_gi + 443);

    auto tr_y_xxxx_yyyzzz = pbuffer.data(idx_dip_gi + 444);

    auto tr_y_xxxx_yyzzzz = pbuffer.data(idx_dip_gi + 445);

    auto tr_y_xxxx_yzzzzz = pbuffer.data(idx_dip_gi + 446);

    auto tr_y_xxxx_zzzzzz = pbuffer.data(idx_dip_gi + 447);

    auto tr_y_xxxy_xxxxxy = pbuffer.data(idx_dip_gi + 449);

    auto tr_y_xxxy_xxxxyy = pbuffer.data(idx_dip_gi + 451);

    auto tr_y_xxxy_xxxxyz = pbuffer.data(idx_dip_gi + 452);

    auto tr_y_xxxy_xxxyyy = pbuffer.data(idx_dip_gi + 454);

    auto tr_y_xxxy_xxxyyz = pbuffer.data(idx_dip_gi + 455);

    auto tr_y_xxxy_xxxyzz = pbuffer.data(idx_dip_gi + 456);

    auto tr_y_xxxy_xxyyyy = pbuffer.data(idx_dip_gi + 458);

    auto tr_y_xxxy_xxyyyz = pbuffer.data(idx_dip_gi + 459);

    auto tr_y_xxxy_xxyyzz = pbuffer.data(idx_dip_gi + 460);

    auto tr_y_xxxy_xxyzzz = pbuffer.data(idx_dip_gi + 461);

    auto tr_y_xxxy_xyyyyy = pbuffer.data(idx_dip_gi + 463);

    auto tr_y_xxxy_xyyyyz = pbuffer.data(idx_dip_gi + 464);

    auto tr_y_xxxy_xyyyzz = pbuffer.data(idx_dip_gi + 465);

    auto tr_y_xxxy_xyyzzz = pbuffer.data(idx_dip_gi + 466);

    auto tr_y_xxxy_xyzzzz = pbuffer.data(idx_dip_gi + 467);

    auto tr_y_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 469);

    auto tr_y_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 470);

    auto tr_y_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 471);

    auto tr_y_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 472);

    auto tr_y_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 473);

    auto tr_y_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 474);

    auto tr_y_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 475);

    auto tr_y_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 476);

    auto tr_y_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 477);

    auto tr_y_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 479);

    auto tr_y_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 482);

    auto tr_y_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 486);

    auto tr_y_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 491);

    auto tr_y_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 498);

    auto tr_y_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 499);

    auto tr_y_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 500);

    auto tr_y_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 501);

    auto tr_y_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 502);

    auto tr_y_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 503);

    auto tr_y_xxyy_xxxxxx = pbuffer.data(idx_dip_gi + 504);

    auto tr_y_xxyy_xxxxxy = pbuffer.data(idx_dip_gi + 505);

    auto tr_y_xxyy_xxxxxz = pbuffer.data(idx_dip_gi + 506);

    auto tr_y_xxyy_xxxxyy = pbuffer.data(idx_dip_gi + 507);

    auto tr_y_xxyy_xxxxyz = pbuffer.data(idx_dip_gi + 508);

    auto tr_y_xxyy_xxxxzz = pbuffer.data(idx_dip_gi + 509);

    auto tr_y_xxyy_xxxyyy = pbuffer.data(idx_dip_gi + 510);

    auto tr_y_xxyy_xxxyyz = pbuffer.data(idx_dip_gi + 511);

    auto tr_y_xxyy_xxxyzz = pbuffer.data(idx_dip_gi + 512);

    auto tr_y_xxyy_xxxzzz = pbuffer.data(idx_dip_gi + 513);

    auto tr_y_xxyy_xxyyyy = pbuffer.data(idx_dip_gi + 514);

    auto tr_y_xxyy_xxyyyz = pbuffer.data(idx_dip_gi + 515);

    auto tr_y_xxyy_xxyyzz = pbuffer.data(idx_dip_gi + 516);

    auto tr_y_xxyy_xxyzzz = pbuffer.data(idx_dip_gi + 517);

    auto tr_y_xxyy_xxzzzz = pbuffer.data(idx_dip_gi + 518);

    auto tr_y_xxyy_xyyyyy = pbuffer.data(idx_dip_gi + 519);

    auto tr_y_xxyy_xyyyyz = pbuffer.data(idx_dip_gi + 520);

    auto tr_y_xxyy_xyyyzz = pbuffer.data(idx_dip_gi + 521);

    auto tr_y_xxyy_xyyzzz = pbuffer.data(idx_dip_gi + 522);

    auto tr_y_xxyy_xyzzzz = pbuffer.data(idx_dip_gi + 523);

    auto tr_y_xxyy_xzzzzz = pbuffer.data(idx_dip_gi + 524);

    auto tr_y_xxyy_yyyyyy = pbuffer.data(idx_dip_gi + 525);

    auto tr_y_xxyy_yyyyyz = pbuffer.data(idx_dip_gi + 526);

    auto tr_y_xxyy_yyyyzz = pbuffer.data(idx_dip_gi + 527);

    auto tr_y_xxyy_yyyzzz = pbuffer.data(idx_dip_gi + 528);

    auto tr_y_xxyy_yyzzzz = pbuffer.data(idx_dip_gi + 529);

    auto tr_y_xxyy_yzzzzz = pbuffer.data(idx_dip_gi + 530);

    auto tr_y_xxyy_zzzzzz = pbuffer.data(idx_dip_gi + 531);

    auto tr_y_xxyz_xxxxxy = pbuffer.data(idx_dip_gi + 533);

    auto tr_y_xxyz_xxxxyy = pbuffer.data(idx_dip_gi + 535);

    auto tr_y_xxyz_xxxyyy = pbuffer.data(idx_dip_gi + 538);

    auto tr_y_xxyz_xxyyyy = pbuffer.data(idx_dip_gi + 542);

    auto tr_y_xxyz_xyyyyy = pbuffer.data(idx_dip_gi + 547);

    auto tr_y_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 554);

    auto tr_y_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 555);

    auto tr_y_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 556);

    auto tr_y_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 557);

    auto tr_y_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 558);

    auto tr_y_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 559);

    auto tr_y_xxzz_xxxxxx = pbuffer.data(idx_dip_gi + 560);

    auto tr_y_xxzz_xxxxxy = pbuffer.data(idx_dip_gi + 561);

    auto tr_y_xxzz_xxxxxz = pbuffer.data(idx_dip_gi + 562);

    auto tr_y_xxzz_xxxxyy = pbuffer.data(idx_dip_gi + 563);

    auto tr_y_xxzz_xxxxyz = pbuffer.data(idx_dip_gi + 564);

    auto tr_y_xxzz_xxxxzz = pbuffer.data(idx_dip_gi + 565);

    auto tr_y_xxzz_xxxyyy = pbuffer.data(idx_dip_gi + 566);

    auto tr_y_xxzz_xxxyyz = pbuffer.data(idx_dip_gi + 567);

    auto tr_y_xxzz_xxxyzz = pbuffer.data(idx_dip_gi + 568);

    auto tr_y_xxzz_xxxzzz = pbuffer.data(idx_dip_gi + 569);

    auto tr_y_xxzz_xxyyyy = pbuffer.data(idx_dip_gi + 570);

    auto tr_y_xxzz_xxyyyz = pbuffer.data(idx_dip_gi + 571);

    auto tr_y_xxzz_xxyyzz = pbuffer.data(idx_dip_gi + 572);

    auto tr_y_xxzz_xxyzzz = pbuffer.data(idx_dip_gi + 573);

    auto tr_y_xxzz_xxzzzz = pbuffer.data(idx_dip_gi + 574);

    auto tr_y_xxzz_xyyyyy = pbuffer.data(idx_dip_gi + 575);

    auto tr_y_xxzz_xyyyyz = pbuffer.data(idx_dip_gi + 576);

    auto tr_y_xxzz_xyyyzz = pbuffer.data(idx_dip_gi + 577);

    auto tr_y_xxzz_xyyzzz = pbuffer.data(idx_dip_gi + 578);

    auto tr_y_xxzz_xyzzzz = pbuffer.data(idx_dip_gi + 579);

    auto tr_y_xxzz_xzzzzz = pbuffer.data(idx_dip_gi + 580);

    auto tr_y_xxzz_yyyyyy = pbuffer.data(idx_dip_gi + 581);

    auto tr_y_xxzz_yyyyyz = pbuffer.data(idx_dip_gi + 582);

    auto tr_y_xxzz_yyyyzz = pbuffer.data(idx_dip_gi + 583);

    auto tr_y_xxzz_yyyzzz = pbuffer.data(idx_dip_gi + 584);

    auto tr_y_xxzz_yyzzzz = pbuffer.data(idx_dip_gi + 585);

    auto tr_y_xxzz_yzzzzz = pbuffer.data(idx_dip_gi + 586);

    auto tr_y_xxzz_zzzzzz = pbuffer.data(idx_dip_gi + 587);

    auto tr_y_xyyy_xxxxxx = pbuffer.data(idx_dip_gi + 588);

    auto tr_y_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 589);

    auto tr_y_xyyy_xxxxxz = pbuffer.data(idx_dip_gi + 590);

    auto tr_y_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 591);

    auto tr_y_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 592);

    auto tr_y_xyyy_xxxxzz = pbuffer.data(idx_dip_gi + 593);

    auto tr_y_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 594);

    auto tr_y_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 595);

    auto tr_y_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 596);

    auto tr_y_xyyy_xxxzzz = pbuffer.data(idx_dip_gi + 597);

    auto tr_y_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 598);

    auto tr_y_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 599);

    auto tr_y_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 600);

    auto tr_y_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 601);

    auto tr_y_xyyy_xxzzzz = pbuffer.data(idx_dip_gi + 602);

    auto tr_y_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 603);

    auto tr_y_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 604);

    auto tr_y_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 605);

    auto tr_y_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 606);

    auto tr_y_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 607);

    auto tr_y_xyyy_xzzzzz = pbuffer.data(idx_dip_gi + 608);

    auto tr_y_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 609);

    auto tr_y_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 610);

    auto tr_y_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 611);

    auto tr_y_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 612);

    auto tr_y_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 613);

    auto tr_y_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 614);

    auto tr_y_xyyy_zzzzzz = pbuffer.data(idx_dip_gi + 615);

    auto tr_y_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 638);

    auto tr_y_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 639);

    auto tr_y_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 640);

    auto tr_y_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 641);

    auto tr_y_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 642);

    auto tr_y_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 643);

    auto tr_y_xyzz_xxxxyz = pbuffer.data(idx_dip_gi + 648);

    auto tr_y_xyzz_xxxyyz = pbuffer.data(idx_dip_gi + 651);

    auto tr_y_xyzz_xxxyzz = pbuffer.data(idx_dip_gi + 652);

    auto tr_y_xyzz_xxyyyz = pbuffer.data(idx_dip_gi + 655);

    auto tr_y_xyzz_xxyyzz = pbuffer.data(idx_dip_gi + 656);

    auto tr_y_xyzz_xxyzzz = pbuffer.data(idx_dip_gi + 657);

    auto tr_y_xyzz_xyyyyz = pbuffer.data(idx_dip_gi + 660);

    auto tr_y_xyzz_xyyyzz = pbuffer.data(idx_dip_gi + 661);

    auto tr_y_xyzz_xyyzzz = pbuffer.data(idx_dip_gi + 662);

    auto tr_y_xyzz_xyzzzz = pbuffer.data(idx_dip_gi + 663);

    auto tr_y_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 665);

    auto tr_y_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 666);

    auto tr_y_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 667);

    auto tr_y_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 668);

    auto tr_y_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 669);

    auto tr_y_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 670);

    auto tr_y_xyzz_zzzzzz = pbuffer.data(idx_dip_gi + 671);

    auto tr_y_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 674);

    auto tr_y_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 676);

    auto tr_y_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 677);

    auto tr_y_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 679);

    auto tr_y_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 680);

    auto tr_y_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 681);

    auto tr_y_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 683);

    auto tr_y_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 684);

    auto tr_y_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 685);

    auto tr_y_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 686);

    auto tr_y_xzzz_xyyyyz = pbuffer.data(idx_dip_gi + 688);

    auto tr_y_xzzz_xyyyzz = pbuffer.data(idx_dip_gi + 689);

    auto tr_y_xzzz_xyyzzz = pbuffer.data(idx_dip_gi + 690);

    auto tr_y_xzzz_xyzzzz = pbuffer.data(idx_dip_gi + 691);

    auto tr_y_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 692);

    auto tr_y_xzzz_yyyyyy = pbuffer.data(idx_dip_gi + 693);

    auto tr_y_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 694);

    auto tr_y_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 695);

    auto tr_y_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 696);

    auto tr_y_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 697);

    auto tr_y_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 698);

    auto tr_y_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 699);

    auto tr_y_yyyy_xxxxxx = pbuffer.data(idx_dip_gi + 700);

    auto tr_y_yyyy_xxxxxy = pbuffer.data(idx_dip_gi + 701);

    auto tr_y_yyyy_xxxxxz = pbuffer.data(idx_dip_gi + 702);

    auto tr_y_yyyy_xxxxyy = pbuffer.data(idx_dip_gi + 703);

    auto tr_y_yyyy_xxxxyz = pbuffer.data(idx_dip_gi + 704);

    auto tr_y_yyyy_xxxxzz = pbuffer.data(idx_dip_gi + 705);

    auto tr_y_yyyy_xxxyyy = pbuffer.data(idx_dip_gi + 706);

    auto tr_y_yyyy_xxxyyz = pbuffer.data(idx_dip_gi + 707);

    auto tr_y_yyyy_xxxyzz = pbuffer.data(idx_dip_gi + 708);

    auto tr_y_yyyy_xxxzzz = pbuffer.data(idx_dip_gi + 709);

    auto tr_y_yyyy_xxyyyy = pbuffer.data(idx_dip_gi + 710);

    auto tr_y_yyyy_xxyyyz = pbuffer.data(idx_dip_gi + 711);

    auto tr_y_yyyy_xxyyzz = pbuffer.data(idx_dip_gi + 712);

    auto tr_y_yyyy_xxyzzz = pbuffer.data(idx_dip_gi + 713);

    auto tr_y_yyyy_xxzzzz = pbuffer.data(idx_dip_gi + 714);

    auto tr_y_yyyy_xyyyyy = pbuffer.data(idx_dip_gi + 715);

    auto tr_y_yyyy_xyyyyz = pbuffer.data(idx_dip_gi + 716);

    auto tr_y_yyyy_xyyyzz = pbuffer.data(idx_dip_gi + 717);

    auto tr_y_yyyy_xyyzzz = pbuffer.data(idx_dip_gi + 718);

    auto tr_y_yyyy_xyzzzz = pbuffer.data(idx_dip_gi + 719);

    auto tr_y_yyyy_xzzzzz = pbuffer.data(idx_dip_gi + 720);

    auto tr_y_yyyy_yyyyyy = pbuffer.data(idx_dip_gi + 721);

    auto tr_y_yyyy_yyyyyz = pbuffer.data(idx_dip_gi + 722);

    auto tr_y_yyyy_yyyyzz = pbuffer.data(idx_dip_gi + 723);

    auto tr_y_yyyy_yyyzzz = pbuffer.data(idx_dip_gi + 724);

    auto tr_y_yyyy_yyzzzz = pbuffer.data(idx_dip_gi + 725);

    auto tr_y_yyyy_yzzzzz = pbuffer.data(idx_dip_gi + 726);

    auto tr_y_yyyy_zzzzzz = pbuffer.data(idx_dip_gi + 727);

    auto tr_y_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 728);

    auto tr_y_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 729);

    auto tr_y_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 731);

    auto tr_y_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 732);

    auto tr_y_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 734);

    auto tr_y_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 735);

    auto tr_y_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 736);

    auto tr_y_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 738);

    auto tr_y_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 739);

    auto tr_y_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 740);

    auto tr_y_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 741);

    auto tr_y_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 743);

    auto tr_y_yyyz_xyyyyz = pbuffer.data(idx_dip_gi + 744);

    auto tr_y_yyyz_xyyyzz = pbuffer.data(idx_dip_gi + 745);

    auto tr_y_yyyz_xyyzzz = pbuffer.data(idx_dip_gi + 746);

    auto tr_y_yyyz_xyzzzz = pbuffer.data(idx_dip_gi + 747);

    auto tr_y_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 749);

    auto tr_y_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 750);

    auto tr_y_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 751);

    auto tr_y_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 752);

    auto tr_y_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 753);

    auto tr_y_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 754);

    auto tr_y_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 755);

    auto tr_y_yyzz_xxxxxx = pbuffer.data(idx_dip_gi + 756);

    auto tr_y_yyzz_xxxxxy = pbuffer.data(idx_dip_gi + 757);

    auto tr_y_yyzz_xxxxxz = pbuffer.data(idx_dip_gi + 758);

    auto tr_y_yyzz_xxxxyy = pbuffer.data(idx_dip_gi + 759);

    auto tr_y_yyzz_xxxxyz = pbuffer.data(idx_dip_gi + 760);

    auto tr_y_yyzz_xxxxzz = pbuffer.data(idx_dip_gi + 761);

    auto tr_y_yyzz_xxxyyy = pbuffer.data(idx_dip_gi + 762);

    auto tr_y_yyzz_xxxyyz = pbuffer.data(idx_dip_gi + 763);

    auto tr_y_yyzz_xxxyzz = pbuffer.data(idx_dip_gi + 764);

    auto tr_y_yyzz_xxxzzz = pbuffer.data(idx_dip_gi + 765);

    auto tr_y_yyzz_xxyyyy = pbuffer.data(idx_dip_gi + 766);

    auto tr_y_yyzz_xxyyyz = pbuffer.data(idx_dip_gi + 767);

    auto tr_y_yyzz_xxyyzz = pbuffer.data(idx_dip_gi + 768);

    auto tr_y_yyzz_xxyzzz = pbuffer.data(idx_dip_gi + 769);

    auto tr_y_yyzz_xxzzzz = pbuffer.data(idx_dip_gi + 770);

    auto tr_y_yyzz_xyyyyy = pbuffer.data(idx_dip_gi + 771);

    auto tr_y_yyzz_xyyyyz = pbuffer.data(idx_dip_gi + 772);

    auto tr_y_yyzz_xyyyzz = pbuffer.data(idx_dip_gi + 773);

    auto tr_y_yyzz_xyyzzz = pbuffer.data(idx_dip_gi + 774);

    auto tr_y_yyzz_xyzzzz = pbuffer.data(idx_dip_gi + 775);

    auto tr_y_yyzz_xzzzzz = pbuffer.data(idx_dip_gi + 776);

    auto tr_y_yyzz_yyyyyy = pbuffer.data(idx_dip_gi + 777);

    auto tr_y_yyzz_yyyyyz = pbuffer.data(idx_dip_gi + 778);

    auto tr_y_yyzz_yyyyzz = pbuffer.data(idx_dip_gi + 779);

    auto tr_y_yyzz_yyyzzz = pbuffer.data(idx_dip_gi + 780);

    auto tr_y_yyzz_yyzzzz = pbuffer.data(idx_dip_gi + 781);

    auto tr_y_yyzz_yzzzzz = pbuffer.data(idx_dip_gi + 782);

    auto tr_y_yyzz_zzzzzz = pbuffer.data(idx_dip_gi + 783);

    auto tr_y_yzzz_xxxxxy = pbuffer.data(idx_dip_gi + 785);

    auto tr_y_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 786);

    auto tr_y_yzzz_xxxxyy = pbuffer.data(idx_dip_gi + 787);

    auto tr_y_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 788);

    auto tr_y_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 789);

    auto tr_y_yzzz_xxxyyy = pbuffer.data(idx_dip_gi + 790);

    auto tr_y_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 791);

    auto tr_y_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 792);

    auto tr_y_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 793);

    auto tr_y_yzzz_xxyyyy = pbuffer.data(idx_dip_gi + 794);

    auto tr_y_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 795);

    auto tr_y_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 796);

    auto tr_y_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 797);

    auto tr_y_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 798);

    auto tr_y_yzzz_xyyyyy = pbuffer.data(idx_dip_gi + 799);

    auto tr_y_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 800);

    auto tr_y_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 801);

    auto tr_y_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 802);

    auto tr_y_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 803);

    auto tr_y_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 804);

    auto tr_y_yzzz_yyyyyy = pbuffer.data(idx_dip_gi + 805);

    auto tr_y_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 806);

    auto tr_y_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 807);

    auto tr_y_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 808);

    auto tr_y_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 809);

    auto tr_y_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 810);

    auto tr_y_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 811);

    auto tr_y_zzzz_xxxxxx = pbuffer.data(idx_dip_gi + 812);

    auto tr_y_zzzz_xxxxxy = pbuffer.data(idx_dip_gi + 813);

    auto tr_y_zzzz_xxxxxz = pbuffer.data(idx_dip_gi + 814);

    auto tr_y_zzzz_xxxxyy = pbuffer.data(idx_dip_gi + 815);

    auto tr_y_zzzz_xxxxyz = pbuffer.data(idx_dip_gi + 816);

    auto tr_y_zzzz_xxxxzz = pbuffer.data(idx_dip_gi + 817);

    auto tr_y_zzzz_xxxyyy = pbuffer.data(idx_dip_gi + 818);

    auto tr_y_zzzz_xxxyyz = pbuffer.data(idx_dip_gi + 819);

    auto tr_y_zzzz_xxxyzz = pbuffer.data(idx_dip_gi + 820);

    auto tr_y_zzzz_xxxzzz = pbuffer.data(idx_dip_gi + 821);

    auto tr_y_zzzz_xxyyyy = pbuffer.data(idx_dip_gi + 822);

    auto tr_y_zzzz_xxyyyz = pbuffer.data(idx_dip_gi + 823);

    auto tr_y_zzzz_xxyyzz = pbuffer.data(idx_dip_gi + 824);

    auto tr_y_zzzz_xxyzzz = pbuffer.data(idx_dip_gi + 825);

    auto tr_y_zzzz_xxzzzz = pbuffer.data(idx_dip_gi + 826);

    auto tr_y_zzzz_xyyyyy = pbuffer.data(idx_dip_gi + 827);

    auto tr_y_zzzz_xyyyyz = pbuffer.data(idx_dip_gi + 828);

    auto tr_y_zzzz_xyyyzz = pbuffer.data(idx_dip_gi + 829);

    auto tr_y_zzzz_xyyzzz = pbuffer.data(idx_dip_gi + 830);

    auto tr_y_zzzz_xyzzzz = pbuffer.data(idx_dip_gi + 831);

    auto tr_y_zzzz_xzzzzz = pbuffer.data(idx_dip_gi + 832);

    auto tr_y_zzzz_yyyyyy = pbuffer.data(idx_dip_gi + 833);

    auto tr_y_zzzz_yyyyyz = pbuffer.data(idx_dip_gi + 834);

    auto tr_y_zzzz_yyyyzz = pbuffer.data(idx_dip_gi + 835);

    auto tr_y_zzzz_yyyzzz = pbuffer.data(idx_dip_gi + 836);

    auto tr_y_zzzz_yyzzzz = pbuffer.data(idx_dip_gi + 837);

    auto tr_y_zzzz_yzzzzz = pbuffer.data(idx_dip_gi + 838);

    auto tr_y_zzzz_zzzzzz = pbuffer.data(idx_dip_gi + 839);

    auto tr_z_xxxx_xxxxxx = pbuffer.data(idx_dip_gi + 840);

    auto tr_z_xxxx_xxxxxy = pbuffer.data(idx_dip_gi + 841);

    auto tr_z_xxxx_xxxxxz = pbuffer.data(idx_dip_gi + 842);

    auto tr_z_xxxx_xxxxyy = pbuffer.data(idx_dip_gi + 843);

    auto tr_z_xxxx_xxxxyz = pbuffer.data(idx_dip_gi + 844);

    auto tr_z_xxxx_xxxxzz = pbuffer.data(idx_dip_gi + 845);

    auto tr_z_xxxx_xxxyyy = pbuffer.data(idx_dip_gi + 846);

    auto tr_z_xxxx_xxxyyz = pbuffer.data(idx_dip_gi + 847);

    auto tr_z_xxxx_xxxyzz = pbuffer.data(idx_dip_gi + 848);

    auto tr_z_xxxx_xxxzzz = pbuffer.data(idx_dip_gi + 849);

    auto tr_z_xxxx_xxyyyy = pbuffer.data(idx_dip_gi + 850);

    auto tr_z_xxxx_xxyyyz = pbuffer.data(idx_dip_gi + 851);

    auto tr_z_xxxx_xxyyzz = pbuffer.data(idx_dip_gi + 852);

    auto tr_z_xxxx_xxyzzz = pbuffer.data(idx_dip_gi + 853);

    auto tr_z_xxxx_xxzzzz = pbuffer.data(idx_dip_gi + 854);

    auto tr_z_xxxx_xyyyyy = pbuffer.data(idx_dip_gi + 855);

    auto tr_z_xxxx_xyyyyz = pbuffer.data(idx_dip_gi + 856);

    auto tr_z_xxxx_xyyyzz = pbuffer.data(idx_dip_gi + 857);

    auto tr_z_xxxx_xyyzzz = pbuffer.data(idx_dip_gi + 858);

    auto tr_z_xxxx_xyzzzz = pbuffer.data(idx_dip_gi + 859);

    auto tr_z_xxxx_xzzzzz = pbuffer.data(idx_dip_gi + 860);

    auto tr_z_xxxx_yyyyyy = pbuffer.data(idx_dip_gi + 861);

    auto tr_z_xxxx_yyyyyz = pbuffer.data(idx_dip_gi + 862);

    auto tr_z_xxxx_yyyyzz = pbuffer.data(idx_dip_gi + 863);

    auto tr_z_xxxx_yyyzzz = pbuffer.data(idx_dip_gi + 864);

    auto tr_z_xxxx_yyzzzz = pbuffer.data(idx_dip_gi + 865);

    auto tr_z_xxxx_yzzzzz = pbuffer.data(idx_dip_gi + 866);

    auto tr_z_xxxx_zzzzzz = pbuffer.data(idx_dip_gi + 867);

    auto tr_z_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 868);

    auto tr_z_xxxy_xxxxxz = pbuffer.data(idx_dip_gi + 870);

    auto tr_z_xxxy_xxxxzz = pbuffer.data(idx_dip_gi + 873);

    auto tr_z_xxxy_xxxzzz = pbuffer.data(idx_dip_gi + 877);

    auto tr_z_xxxy_xxzzzz = pbuffer.data(idx_dip_gi + 882);

    auto tr_z_xxxy_xzzzzz = pbuffer.data(idx_dip_gi + 888);

    auto tr_z_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 889);

    auto tr_z_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 890);

    auto tr_z_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 891);

    auto tr_z_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 892);

    auto tr_z_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 893);

    auto tr_z_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 894);

    auto tr_z_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 896);

    auto tr_z_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 898);

    auto tr_z_xxxz_xxxxyz = pbuffer.data(idx_dip_gi + 900);

    auto tr_z_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 901);

    auto tr_z_xxxz_xxxyyz = pbuffer.data(idx_dip_gi + 903);

    auto tr_z_xxxz_xxxyzz = pbuffer.data(idx_dip_gi + 904);

    auto tr_z_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 905);

    auto tr_z_xxxz_xxyyyz = pbuffer.data(idx_dip_gi + 907);

    auto tr_z_xxxz_xxyyzz = pbuffer.data(idx_dip_gi + 908);

    auto tr_z_xxxz_xxyzzz = pbuffer.data(idx_dip_gi + 909);

    auto tr_z_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 910);

    auto tr_z_xxxz_xyyyyz = pbuffer.data(idx_dip_gi + 912);

    auto tr_z_xxxz_xyyyzz = pbuffer.data(idx_dip_gi + 913);

    auto tr_z_xxxz_xyyzzz = pbuffer.data(idx_dip_gi + 914);

    auto tr_z_xxxz_xyzzzz = pbuffer.data(idx_dip_gi + 915);

    auto tr_z_xxxz_xzzzzz = pbuffer.data(idx_dip_gi + 916);

    auto tr_z_xxxz_yyyyyy = pbuffer.data(idx_dip_gi + 917);

    auto tr_z_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 918);

    auto tr_z_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 919);

    auto tr_z_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 920);

    auto tr_z_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 921);

    auto tr_z_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 922);

    auto tr_z_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 923);

    auto tr_z_xxyy_xxxxxx = pbuffer.data(idx_dip_gi + 924);

    auto tr_z_xxyy_xxxxxy = pbuffer.data(idx_dip_gi + 925);

    auto tr_z_xxyy_xxxxxz = pbuffer.data(idx_dip_gi + 926);

    auto tr_z_xxyy_xxxxyy = pbuffer.data(idx_dip_gi + 927);

    auto tr_z_xxyy_xxxxyz = pbuffer.data(idx_dip_gi + 928);

    auto tr_z_xxyy_xxxxzz = pbuffer.data(idx_dip_gi + 929);

    auto tr_z_xxyy_xxxyyy = pbuffer.data(idx_dip_gi + 930);

    auto tr_z_xxyy_xxxyyz = pbuffer.data(idx_dip_gi + 931);

    auto tr_z_xxyy_xxxyzz = pbuffer.data(idx_dip_gi + 932);

    auto tr_z_xxyy_xxxzzz = pbuffer.data(idx_dip_gi + 933);

    auto tr_z_xxyy_xxyyyy = pbuffer.data(idx_dip_gi + 934);

    auto tr_z_xxyy_xxyyyz = pbuffer.data(idx_dip_gi + 935);

    auto tr_z_xxyy_xxyyzz = pbuffer.data(idx_dip_gi + 936);

    auto tr_z_xxyy_xxyzzz = pbuffer.data(idx_dip_gi + 937);

    auto tr_z_xxyy_xxzzzz = pbuffer.data(idx_dip_gi + 938);

    auto tr_z_xxyy_xyyyyy = pbuffer.data(idx_dip_gi + 939);

    auto tr_z_xxyy_xyyyyz = pbuffer.data(idx_dip_gi + 940);

    auto tr_z_xxyy_xyyyzz = pbuffer.data(idx_dip_gi + 941);

    auto tr_z_xxyy_xyyzzz = pbuffer.data(idx_dip_gi + 942);

    auto tr_z_xxyy_xyzzzz = pbuffer.data(idx_dip_gi + 943);

    auto tr_z_xxyy_xzzzzz = pbuffer.data(idx_dip_gi + 944);

    auto tr_z_xxyy_yyyyyy = pbuffer.data(idx_dip_gi + 945);

    auto tr_z_xxyy_yyyyyz = pbuffer.data(idx_dip_gi + 946);

    auto tr_z_xxyy_yyyyzz = pbuffer.data(idx_dip_gi + 947);

    auto tr_z_xxyy_yyyzzz = pbuffer.data(idx_dip_gi + 948);

    auto tr_z_xxyy_yyzzzz = pbuffer.data(idx_dip_gi + 949);

    auto tr_z_xxyy_yzzzzz = pbuffer.data(idx_dip_gi + 950);

    auto tr_z_xxyy_zzzzzz = pbuffer.data(idx_dip_gi + 951);

    auto tr_z_xxyz_xxxxxx = pbuffer.data(idx_dip_gi + 952);

    auto tr_z_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 954);

    auto tr_z_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 957);

    auto tr_z_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 961);

    auto tr_z_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 966);

    auto tr_z_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 972);

    auto tr_z_xxyz_yyyyyy = pbuffer.data(idx_dip_gi + 973);

    auto tr_z_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 974);

    auto tr_z_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 975);

    auto tr_z_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 976);

    auto tr_z_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 977);

    auto tr_z_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 978);

    auto tr_z_xxzz_xxxxxx = pbuffer.data(idx_dip_gi + 980);

    auto tr_z_xxzz_xxxxxy = pbuffer.data(idx_dip_gi + 981);

    auto tr_z_xxzz_xxxxxz = pbuffer.data(idx_dip_gi + 982);

    auto tr_z_xxzz_xxxxyy = pbuffer.data(idx_dip_gi + 983);

    auto tr_z_xxzz_xxxxyz = pbuffer.data(idx_dip_gi + 984);

    auto tr_z_xxzz_xxxxzz = pbuffer.data(idx_dip_gi + 985);

    auto tr_z_xxzz_xxxyyy = pbuffer.data(idx_dip_gi + 986);

    auto tr_z_xxzz_xxxyyz = pbuffer.data(idx_dip_gi + 987);

    auto tr_z_xxzz_xxxyzz = pbuffer.data(idx_dip_gi + 988);

    auto tr_z_xxzz_xxxzzz = pbuffer.data(idx_dip_gi + 989);

    auto tr_z_xxzz_xxyyyy = pbuffer.data(idx_dip_gi + 990);

    auto tr_z_xxzz_xxyyyz = pbuffer.data(idx_dip_gi + 991);

    auto tr_z_xxzz_xxyyzz = pbuffer.data(idx_dip_gi + 992);

    auto tr_z_xxzz_xxyzzz = pbuffer.data(idx_dip_gi + 993);

    auto tr_z_xxzz_xxzzzz = pbuffer.data(idx_dip_gi + 994);

    auto tr_z_xxzz_xyyyyy = pbuffer.data(idx_dip_gi + 995);

    auto tr_z_xxzz_xyyyyz = pbuffer.data(idx_dip_gi + 996);

    auto tr_z_xxzz_xyyyzz = pbuffer.data(idx_dip_gi + 997);

    auto tr_z_xxzz_xyyzzz = pbuffer.data(idx_dip_gi + 998);

    auto tr_z_xxzz_xyzzzz = pbuffer.data(idx_dip_gi + 999);

    auto tr_z_xxzz_xzzzzz = pbuffer.data(idx_dip_gi + 1000);

    auto tr_z_xxzz_yyyyyy = pbuffer.data(idx_dip_gi + 1001);

    auto tr_z_xxzz_yyyyyz = pbuffer.data(idx_dip_gi + 1002);

    auto tr_z_xxzz_yyyyzz = pbuffer.data(idx_dip_gi + 1003);

    auto tr_z_xxzz_yyyzzz = pbuffer.data(idx_dip_gi + 1004);

    auto tr_z_xxzz_yyzzzz = pbuffer.data(idx_dip_gi + 1005);

    auto tr_z_xxzz_yzzzzz = pbuffer.data(idx_dip_gi + 1006);

    auto tr_z_xxzz_zzzzzz = pbuffer.data(idx_dip_gi + 1007);

    auto tr_z_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 1009);

    auto tr_z_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 1011);

    auto tr_z_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 1012);

    auto tr_z_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 1014);

    auto tr_z_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 1015);

    auto tr_z_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 1016);

    auto tr_z_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 1018);

    auto tr_z_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 1019);

    auto tr_z_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 1020);

    auto tr_z_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 1021);

    auto tr_z_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 1023);

    auto tr_z_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 1024);

    auto tr_z_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 1025);

    auto tr_z_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 1026);

    auto tr_z_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 1027);

    auto tr_z_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 1029);

    auto tr_z_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 1030);

    auto tr_z_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 1031);

    auto tr_z_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 1032);

    auto tr_z_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 1033);

    auto tr_z_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 1034);

    auto tr_z_xyyy_zzzzzz = pbuffer.data(idx_dip_gi + 1035);

    auto tr_z_xyyz_xxxxyz = pbuffer.data(idx_dip_gi + 1040);

    auto tr_z_xyyz_xxxyyz = pbuffer.data(idx_dip_gi + 1043);

    auto tr_z_xyyz_xxxyzz = pbuffer.data(idx_dip_gi + 1044);

    auto tr_z_xyyz_xxyyyz = pbuffer.data(idx_dip_gi + 1047);

    auto tr_z_xyyz_xxyyzz = pbuffer.data(idx_dip_gi + 1048);

    auto tr_z_xyyz_xxyzzz = pbuffer.data(idx_dip_gi + 1049);

    auto tr_z_xyyz_xyyyyz = pbuffer.data(idx_dip_gi + 1052);

    auto tr_z_xyyz_xyyyzz = pbuffer.data(idx_dip_gi + 1053);

    auto tr_z_xyyz_xyyzzz = pbuffer.data(idx_dip_gi + 1054);

    auto tr_z_xyyz_xyzzzz = pbuffer.data(idx_dip_gi + 1055);

    auto tr_z_xyyz_yyyyyy = pbuffer.data(idx_dip_gi + 1057);

    auto tr_z_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 1058);

    auto tr_z_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 1059);

    auto tr_z_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 1060);

    auto tr_z_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 1061);

    auto tr_z_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 1062);

    auto tr_z_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 1063);

    auto tr_z_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 1085);

    auto tr_z_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 1086);

    auto tr_z_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 1087);

    auto tr_z_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 1088);

    auto tr_z_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 1089);

    auto tr_z_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 1090);

    auto tr_z_xzzz_xxxxxx = pbuffer.data(idx_dip_gi + 1092);

    auto tr_z_xzzz_xxxxxy = pbuffer.data(idx_dip_gi + 1093);

    auto tr_z_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 1094);

    auto tr_z_xzzz_xxxxyy = pbuffer.data(idx_dip_gi + 1095);

    auto tr_z_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 1096);

    auto tr_z_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 1097);

    auto tr_z_xzzz_xxxyyy = pbuffer.data(idx_dip_gi + 1098);

    auto tr_z_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 1099);

    auto tr_z_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 1100);

    auto tr_z_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 1101);

    auto tr_z_xzzz_xxyyyy = pbuffer.data(idx_dip_gi + 1102);

    auto tr_z_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 1103);

    auto tr_z_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 1104);

    auto tr_z_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 1105);

    auto tr_z_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 1106);

    auto tr_z_xzzz_xyyyyy = pbuffer.data(idx_dip_gi + 1107);

    auto tr_z_xzzz_xyyyyz = pbuffer.data(idx_dip_gi + 1108);

    auto tr_z_xzzz_xyyyzz = pbuffer.data(idx_dip_gi + 1109);

    auto tr_z_xzzz_xyyzzz = pbuffer.data(idx_dip_gi + 1110);

    auto tr_z_xzzz_xyzzzz = pbuffer.data(idx_dip_gi + 1111);

    auto tr_z_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 1112);

    auto tr_z_xzzz_yyyyyy = pbuffer.data(idx_dip_gi + 1113);

    auto tr_z_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 1114);

    auto tr_z_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 1115);

    auto tr_z_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 1116);

    auto tr_z_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 1117);

    auto tr_z_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 1118);

    auto tr_z_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 1119);

    auto tr_z_yyyy_xxxxxx = pbuffer.data(idx_dip_gi + 1120);

    auto tr_z_yyyy_xxxxxy = pbuffer.data(idx_dip_gi + 1121);

    auto tr_z_yyyy_xxxxxz = pbuffer.data(idx_dip_gi + 1122);

    auto tr_z_yyyy_xxxxyy = pbuffer.data(idx_dip_gi + 1123);

    auto tr_z_yyyy_xxxxyz = pbuffer.data(idx_dip_gi + 1124);

    auto tr_z_yyyy_xxxxzz = pbuffer.data(idx_dip_gi + 1125);

    auto tr_z_yyyy_xxxyyy = pbuffer.data(idx_dip_gi + 1126);

    auto tr_z_yyyy_xxxyyz = pbuffer.data(idx_dip_gi + 1127);

    auto tr_z_yyyy_xxxyzz = pbuffer.data(idx_dip_gi + 1128);

    auto tr_z_yyyy_xxxzzz = pbuffer.data(idx_dip_gi + 1129);

    auto tr_z_yyyy_xxyyyy = pbuffer.data(idx_dip_gi + 1130);

    auto tr_z_yyyy_xxyyyz = pbuffer.data(idx_dip_gi + 1131);

    auto tr_z_yyyy_xxyyzz = pbuffer.data(idx_dip_gi + 1132);

    auto tr_z_yyyy_xxyzzz = pbuffer.data(idx_dip_gi + 1133);

    auto tr_z_yyyy_xxzzzz = pbuffer.data(idx_dip_gi + 1134);

    auto tr_z_yyyy_xyyyyy = pbuffer.data(idx_dip_gi + 1135);

    auto tr_z_yyyy_xyyyyz = pbuffer.data(idx_dip_gi + 1136);

    auto tr_z_yyyy_xyyyzz = pbuffer.data(idx_dip_gi + 1137);

    auto tr_z_yyyy_xyyzzz = pbuffer.data(idx_dip_gi + 1138);

    auto tr_z_yyyy_xyzzzz = pbuffer.data(idx_dip_gi + 1139);

    auto tr_z_yyyy_xzzzzz = pbuffer.data(idx_dip_gi + 1140);

    auto tr_z_yyyy_yyyyyy = pbuffer.data(idx_dip_gi + 1141);

    auto tr_z_yyyy_yyyyyz = pbuffer.data(idx_dip_gi + 1142);

    auto tr_z_yyyy_yyyyzz = pbuffer.data(idx_dip_gi + 1143);

    auto tr_z_yyyy_yyyzzz = pbuffer.data(idx_dip_gi + 1144);

    auto tr_z_yyyy_yyzzzz = pbuffer.data(idx_dip_gi + 1145);

    auto tr_z_yyyy_yzzzzz = pbuffer.data(idx_dip_gi + 1146);

    auto tr_z_yyyy_zzzzzz = pbuffer.data(idx_dip_gi + 1147);

    auto tr_z_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 1148);

    auto tr_z_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 1150);

    auto tr_z_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 1152);

    auto tr_z_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 1153);

    auto tr_z_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 1155);

    auto tr_z_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 1156);

    auto tr_z_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 1157);

    auto tr_z_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 1159);

    auto tr_z_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 1160);

    auto tr_z_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 1161);

    auto tr_z_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 1162);

    auto tr_z_yyyz_xyyyyz = pbuffer.data(idx_dip_gi + 1164);

    auto tr_z_yyyz_xyyyzz = pbuffer.data(idx_dip_gi + 1165);

    auto tr_z_yyyz_xyyzzz = pbuffer.data(idx_dip_gi + 1166);

    auto tr_z_yyyz_xyzzzz = pbuffer.data(idx_dip_gi + 1167);

    auto tr_z_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 1168);

    auto tr_z_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 1169);

    auto tr_z_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 1170);

    auto tr_z_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 1171);

    auto tr_z_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 1172);

    auto tr_z_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 1173);

    auto tr_z_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 1174);

    auto tr_z_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 1175);

    auto tr_z_yyzz_xxxxxx = pbuffer.data(idx_dip_gi + 1176);

    auto tr_z_yyzz_xxxxxy = pbuffer.data(idx_dip_gi + 1177);

    auto tr_z_yyzz_xxxxxz = pbuffer.data(idx_dip_gi + 1178);

    auto tr_z_yyzz_xxxxyy = pbuffer.data(idx_dip_gi + 1179);

    auto tr_z_yyzz_xxxxyz = pbuffer.data(idx_dip_gi + 1180);

    auto tr_z_yyzz_xxxxzz = pbuffer.data(idx_dip_gi + 1181);

    auto tr_z_yyzz_xxxyyy = pbuffer.data(idx_dip_gi + 1182);

    auto tr_z_yyzz_xxxyyz = pbuffer.data(idx_dip_gi + 1183);

    auto tr_z_yyzz_xxxyzz = pbuffer.data(idx_dip_gi + 1184);

    auto tr_z_yyzz_xxxzzz = pbuffer.data(idx_dip_gi + 1185);

    auto tr_z_yyzz_xxyyyy = pbuffer.data(idx_dip_gi + 1186);

    auto tr_z_yyzz_xxyyyz = pbuffer.data(idx_dip_gi + 1187);

    auto tr_z_yyzz_xxyyzz = pbuffer.data(idx_dip_gi + 1188);

    auto tr_z_yyzz_xxyzzz = pbuffer.data(idx_dip_gi + 1189);

    auto tr_z_yyzz_xxzzzz = pbuffer.data(idx_dip_gi + 1190);

    auto tr_z_yyzz_xyyyyy = pbuffer.data(idx_dip_gi + 1191);

    auto tr_z_yyzz_xyyyyz = pbuffer.data(idx_dip_gi + 1192);

    auto tr_z_yyzz_xyyyzz = pbuffer.data(idx_dip_gi + 1193);

    auto tr_z_yyzz_xyyzzz = pbuffer.data(idx_dip_gi + 1194);

    auto tr_z_yyzz_xyzzzz = pbuffer.data(idx_dip_gi + 1195);

    auto tr_z_yyzz_xzzzzz = pbuffer.data(idx_dip_gi + 1196);

    auto tr_z_yyzz_yyyyyy = pbuffer.data(idx_dip_gi + 1197);

    auto tr_z_yyzz_yyyyyz = pbuffer.data(idx_dip_gi + 1198);

    auto tr_z_yyzz_yyyyzz = pbuffer.data(idx_dip_gi + 1199);

    auto tr_z_yyzz_yyyzzz = pbuffer.data(idx_dip_gi + 1200);

    auto tr_z_yyzz_yyzzzz = pbuffer.data(idx_dip_gi + 1201);

    auto tr_z_yyzz_yzzzzz = pbuffer.data(idx_dip_gi + 1202);

    auto tr_z_yyzz_zzzzzz = pbuffer.data(idx_dip_gi + 1203);

    auto tr_z_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 1204);

    auto tr_z_yzzz_xxxxxy = pbuffer.data(idx_dip_gi + 1205);

    auto tr_z_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 1206);

    auto tr_z_yzzz_xxxxyy = pbuffer.data(idx_dip_gi + 1207);

    auto tr_z_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 1208);

    auto tr_z_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 1209);

    auto tr_z_yzzz_xxxyyy = pbuffer.data(idx_dip_gi + 1210);

    auto tr_z_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 1211);

    auto tr_z_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 1212);

    auto tr_z_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 1213);

    auto tr_z_yzzz_xxyyyy = pbuffer.data(idx_dip_gi + 1214);

    auto tr_z_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 1215);

    auto tr_z_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 1216);

    auto tr_z_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 1217);

    auto tr_z_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 1218);

    auto tr_z_yzzz_xyyyyy = pbuffer.data(idx_dip_gi + 1219);

    auto tr_z_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 1220);

    auto tr_z_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 1221);

    auto tr_z_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 1222);

    auto tr_z_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 1223);

    auto tr_z_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 1224);

    auto tr_z_yzzz_yyyyyy = pbuffer.data(idx_dip_gi + 1225);

    auto tr_z_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 1226);

    auto tr_z_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 1227);

    auto tr_z_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 1228);

    auto tr_z_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 1229);

    auto tr_z_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 1230);

    auto tr_z_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 1231);

    auto tr_z_zzzz_xxxxxx = pbuffer.data(idx_dip_gi + 1232);

    auto tr_z_zzzz_xxxxxy = pbuffer.data(idx_dip_gi + 1233);

    auto tr_z_zzzz_xxxxxz = pbuffer.data(idx_dip_gi + 1234);

    auto tr_z_zzzz_xxxxyy = pbuffer.data(idx_dip_gi + 1235);

    auto tr_z_zzzz_xxxxyz = pbuffer.data(idx_dip_gi + 1236);

    auto tr_z_zzzz_xxxxzz = pbuffer.data(idx_dip_gi + 1237);

    auto tr_z_zzzz_xxxyyy = pbuffer.data(idx_dip_gi + 1238);

    auto tr_z_zzzz_xxxyyz = pbuffer.data(idx_dip_gi + 1239);

    auto tr_z_zzzz_xxxyzz = pbuffer.data(idx_dip_gi + 1240);

    auto tr_z_zzzz_xxxzzz = pbuffer.data(idx_dip_gi + 1241);

    auto tr_z_zzzz_xxyyyy = pbuffer.data(idx_dip_gi + 1242);

    auto tr_z_zzzz_xxyyyz = pbuffer.data(idx_dip_gi + 1243);

    auto tr_z_zzzz_xxyyzz = pbuffer.data(idx_dip_gi + 1244);

    auto tr_z_zzzz_xxyzzz = pbuffer.data(idx_dip_gi + 1245);

    auto tr_z_zzzz_xxzzzz = pbuffer.data(idx_dip_gi + 1246);

    auto tr_z_zzzz_xyyyyy = pbuffer.data(idx_dip_gi + 1247);

    auto tr_z_zzzz_xyyyyz = pbuffer.data(idx_dip_gi + 1248);

    auto tr_z_zzzz_xyyyzz = pbuffer.data(idx_dip_gi + 1249);

    auto tr_z_zzzz_xyyzzz = pbuffer.data(idx_dip_gi + 1250);

    auto tr_z_zzzz_xyzzzz = pbuffer.data(idx_dip_gi + 1251);

    auto tr_z_zzzz_xzzzzz = pbuffer.data(idx_dip_gi + 1252);

    auto tr_z_zzzz_yyyyyy = pbuffer.data(idx_dip_gi + 1253);

    auto tr_z_zzzz_yyyyyz = pbuffer.data(idx_dip_gi + 1254);

    auto tr_z_zzzz_yyyyzz = pbuffer.data(idx_dip_gi + 1255);

    auto tr_z_zzzz_yyyzzz = pbuffer.data(idx_dip_gi + 1256);

    auto tr_z_zzzz_yyzzzz = pbuffer.data(idx_dip_gi + 1257);

    auto tr_z_zzzz_yzzzzz = pbuffer.data(idx_dip_gi + 1258);

    auto tr_z_zzzz_zzzzzz = pbuffer.data(idx_dip_gi + 1259);

    // Set up components of auxiliary buffer : HH

    auto tr_x_xxxxx_xxxxx = pbuffer.data(idx_dip_hh);

    auto tr_x_xxxxx_xxxxy = pbuffer.data(idx_dip_hh + 1);

    auto tr_x_xxxxx_xxxxz = pbuffer.data(idx_dip_hh + 2);

    auto tr_x_xxxxx_xxxyy = pbuffer.data(idx_dip_hh + 3);

    auto tr_x_xxxxx_xxxyz = pbuffer.data(idx_dip_hh + 4);

    auto tr_x_xxxxx_xxxzz = pbuffer.data(idx_dip_hh + 5);

    auto tr_x_xxxxx_xxyyy = pbuffer.data(idx_dip_hh + 6);

    auto tr_x_xxxxx_xxyyz = pbuffer.data(idx_dip_hh + 7);

    auto tr_x_xxxxx_xxyzz = pbuffer.data(idx_dip_hh + 8);

    auto tr_x_xxxxx_xxzzz = pbuffer.data(idx_dip_hh + 9);

    auto tr_x_xxxxx_xyyyy = pbuffer.data(idx_dip_hh + 10);

    auto tr_x_xxxxx_xyyyz = pbuffer.data(idx_dip_hh + 11);

    auto tr_x_xxxxx_xyyzz = pbuffer.data(idx_dip_hh + 12);

    auto tr_x_xxxxx_xyzzz = pbuffer.data(idx_dip_hh + 13);

    auto tr_x_xxxxx_xzzzz = pbuffer.data(idx_dip_hh + 14);

    auto tr_x_xxxxx_yyyyy = pbuffer.data(idx_dip_hh + 15);

    auto tr_x_xxxxx_yyyyz = pbuffer.data(idx_dip_hh + 16);

    auto tr_x_xxxxx_yyyzz = pbuffer.data(idx_dip_hh + 17);

    auto tr_x_xxxxx_yyzzz = pbuffer.data(idx_dip_hh + 18);

    auto tr_x_xxxxx_yzzzz = pbuffer.data(idx_dip_hh + 19);

    auto tr_x_xxxxx_zzzzz = pbuffer.data(idx_dip_hh + 20);

    auto tr_x_xxxxy_xxxxx = pbuffer.data(idx_dip_hh + 21);

    auto tr_x_xxxxy_xxxxy = pbuffer.data(idx_dip_hh + 22);

    auto tr_x_xxxxy_xxxxz = pbuffer.data(idx_dip_hh + 23);

    auto tr_x_xxxxy_xxxyy = pbuffer.data(idx_dip_hh + 24);

    auto tr_x_xxxxy_xxxyz = pbuffer.data(idx_dip_hh + 25);

    auto tr_x_xxxxy_xxxzz = pbuffer.data(idx_dip_hh + 26);

    auto tr_x_xxxxy_xxyyy = pbuffer.data(idx_dip_hh + 27);

    auto tr_x_xxxxy_xxyyz = pbuffer.data(idx_dip_hh + 28);

    auto tr_x_xxxxy_xxyzz = pbuffer.data(idx_dip_hh + 29);

    auto tr_x_xxxxy_xxzzz = pbuffer.data(idx_dip_hh + 30);

    auto tr_x_xxxxy_xyyyy = pbuffer.data(idx_dip_hh + 31);

    auto tr_x_xxxxy_xyyyz = pbuffer.data(idx_dip_hh + 32);

    auto tr_x_xxxxy_xyyzz = pbuffer.data(idx_dip_hh + 33);

    auto tr_x_xxxxy_xyzzz = pbuffer.data(idx_dip_hh + 34);

    auto tr_x_xxxxy_xzzzz = pbuffer.data(idx_dip_hh + 35);

    auto tr_x_xxxxz_xxxxx = pbuffer.data(idx_dip_hh + 42);

    auto tr_x_xxxxz_xxxxy = pbuffer.data(idx_dip_hh + 43);

    auto tr_x_xxxxz_xxxxz = pbuffer.data(idx_dip_hh + 44);

    auto tr_x_xxxxz_xxxyy = pbuffer.data(idx_dip_hh + 45);

    auto tr_x_xxxxz_xxxyz = pbuffer.data(idx_dip_hh + 46);

    auto tr_x_xxxxz_xxxzz = pbuffer.data(idx_dip_hh + 47);

    auto tr_x_xxxxz_xxyyy = pbuffer.data(idx_dip_hh + 48);

    auto tr_x_xxxxz_xxyyz = pbuffer.data(idx_dip_hh + 49);

    auto tr_x_xxxxz_xxyzz = pbuffer.data(idx_dip_hh + 50);

    auto tr_x_xxxxz_xxzzz = pbuffer.data(idx_dip_hh + 51);

    auto tr_x_xxxxz_xyyyy = pbuffer.data(idx_dip_hh + 52);

    auto tr_x_xxxxz_xyyyz = pbuffer.data(idx_dip_hh + 53);

    auto tr_x_xxxxz_xyyzz = pbuffer.data(idx_dip_hh + 54);

    auto tr_x_xxxxz_xyzzz = pbuffer.data(idx_dip_hh + 55);

    auto tr_x_xxxxz_xzzzz = pbuffer.data(idx_dip_hh + 56);

    auto tr_x_xxxxz_yyyyz = pbuffer.data(idx_dip_hh + 58);

    auto tr_x_xxxxz_yyyzz = pbuffer.data(idx_dip_hh + 59);

    auto tr_x_xxxxz_yyzzz = pbuffer.data(idx_dip_hh + 60);

    auto tr_x_xxxxz_yzzzz = pbuffer.data(idx_dip_hh + 61);

    auto tr_x_xxxxz_zzzzz = pbuffer.data(idx_dip_hh + 62);

    auto tr_x_xxxyy_xxxxx = pbuffer.data(idx_dip_hh + 63);

    auto tr_x_xxxyy_xxxxy = pbuffer.data(idx_dip_hh + 64);

    auto tr_x_xxxyy_xxxxz = pbuffer.data(idx_dip_hh + 65);

    auto tr_x_xxxyy_xxxyy = pbuffer.data(idx_dip_hh + 66);

    auto tr_x_xxxyy_xxxyz = pbuffer.data(idx_dip_hh + 67);

    auto tr_x_xxxyy_xxxzz = pbuffer.data(idx_dip_hh + 68);

    auto tr_x_xxxyy_xxyyy = pbuffer.data(idx_dip_hh + 69);

    auto tr_x_xxxyy_xxyyz = pbuffer.data(idx_dip_hh + 70);

    auto tr_x_xxxyy_xxyzz = pbuffer.data(idx_dip_hh + 71);

    auto tr_x_xxxyy_xxzzz = pbuffer.data(idx_dip_hh + 72);

    auto tr_x_xxxyy_xyyyy = pbuffer.data(idx_dip_hh + 73);

    auto tr_x_xxxyy_xyyyz = pbuffer.data(idx_dip_hh + 74);

    auto tr_x_xxxyy_xyyzz = pbuffer.data(idx_dip_hh + 75);

    auto tr_x_xxxyy_xyzzz = pbuffer.data(idx_dip_hh + 76);

    auto tr_x_xxxyy_xzzzz = pbuffer.data(idx_dip_hh + 77);

    auto tr_x_xxxyy_yyyyy = pbuffer.data(idx_dip_hh + 78);

    auto tr_x_xxxyy_yyyyz = pbuffer.data(idx_dip_hh + 79);

    auto tr_x_xxxyy_yyyzz = pbuffer.data(idx_dip_hh + 80);

    auto tr_x_xxxyy_yyzzz = pbuffer.data(idx_dip_hh + 81);

    auto tr_x_xxxyy_yzzzz = pbuffer.data(idx_dip_hh + 82);

    auto tr_x_xxxzz_xxxxx = pbuffer.data(idx_dip_hh + 105);

    auto tr_x_xxxzz_xxxxy = pbuffer.data(idx_dip_hh + 106);

    auto tr_x_xxxzz_xxxxz = pbuffer.data(idx_dip_hh + 107);

    auto tr_x_xxxzz_xxxyy = pbuffer.data(idx_dip_hh + 108);

    auto tr_x_xxxzz_xxxyz = pbuffer.data(idx_dip_hh + 109);

    auto tr_x_xxxzz_xxxzz = pbuffer.data(idx_dip_hh + 110);

    auto tr_x_xxxzz_xxyyy = pbuffer.data(idx_dip_hh + 111);

    auto tr_x_xxxzz_xxyyz = pbuffer.data(idx_dip_hh + 112);

    auto tr_x_xxxzz_xxyzz = pbuffer.data(idx_dip_hh + 113);

    auto tr_x_xxxzz_xxzzz = pbuffer.data(idx_dip_hh + 114);

    auto tr_x_xxxzz_xyyyy = pbuffer.data(idx_dip_hh + 115);

    auto tr_x_xxxzz_xyyyz = pbuffer.data(idx_dip_hh + 116);

    auto tr_x_xxxzz_xyyzz = pbuffer.data(idx_dip_hh + 117);

    auto tr_x_xxxzz_xyzzz = pbuffer.data(idx_dip_hh + 118);

    auto tr_x_xxxzz_xzzzz = pbuffer.data(idx_dip_hh + 119);

    auto tr_x_xxxzz_yyyyy = pbuffer.data(idx_dip_hh + 120);

    auto tr_x_xxxzz_yyyyz = pbuffer.data(idx_dip_hh + 121);

    auto tr_x_xxxzz_yyyzz = pbuffer.data(idx_dip_hh + 122);

    auto tr_x_xxxzz_yyzzz = pbuffer.data(idx_dip_hh + 123);

    auto tr_x_xxxzz_yzzzz = pbuffer.data(idx_dip_hh + 124);

    auto tr_x_xxxzz_zzzzz = pbuffer.data(idx_dip_hh + 125);

    auto tr_x_xxyyy_xxxxx = pbuffer.data(idx_dip_hh + 126);

    auto tr_x_xxyyy_xxxxy = pbuffer.data(idx_dip_hh + 127);

    auto tr_x_xxyyy_xxxxz = pbuffer.data(idx_dip_hh + 128);

    auto tr_x_xxyyy_xxxyy = pbuffer.data(idx_dip_hh + 129);

    auto tr_x_xxyyy_xxxyz = pbuffer.data(idx_dip_hh + 130);

    auto tr_x_xxyyy_xxxzz = pbuffer.data(idx_dip_hh + 131);

    auto tr_x_xxyyy_xxyyy = pbuffer.data(idx_dip_hh + 132);

    auto tr_x_xxyyy_xxyyz = pbuffer.data(idx_dip_hh + 133);

    auto tr_x_xxyyy_xxyzz = pbuffer.data(idx_dip_hh + 134);

    auto tr_x_xxyyy_xxzzz = pbuffer.data(idx_dip_hh + 135);

    auto tr_x_xxyyy_xyyyy = pbuffer.data(idx_dip_hh + 136);

    auto tr_x_xxyyy_xyyyz = pbuffer.data(idx_dip_hh + 137);

    auto tr_x_xxyyy_xyyzz = pbuffer.data(idx_dip_hh + 138);

    auto tr_x_xxyyy_xyzzz = pbuffer.data(idx_dip_hh + 139);

    auto tr_x_xxyyy_xzzzz = pbuffer.data(idx_dip_hh + 140);

    auto tr_x_xxyyy_yyyyy = pbuffer.data(idx_dip_hh + 141);

    auto tr_x_xxyyy_yyyyz = pbuffer.data(idx_dip_hh + 142);

    auto tr_x_xxyyy_yyyzz = pbuffer.data(idx_dip_hh + 143);

    auto tr_x_xxyyy_yyzzz = pbuffer.data(idx_dip_hh + 144);

    auto tr_x_xxyyy_yzzzz = pbuffer.data(idx_dip_hh + 145);

    auto tr_x_xxyzz_xxxxz = pbuffer.data(idx_dip_hh + 170);

    auto tr_x_xxyzz_xxxyz = pbuffer.data(idx_dip_hh + 172);

    auto tr_x_xxyzz_xxxzz = pbuffer.data(idx_dip_hh + 173);

    auto tr_x_xxyzz_xxyyz = pbuffer.data(idx_dip_hh + 175);

    auto tr_x_xxyzz_xxyzz = pbuffer.data(idx_dip_hh + 176);

    auto tr_x_xxyzz_xxzzz = pbuffer.data(idx_dip_hh + 177);

    auto tr_x_xxyzz_xyyyz = pbuffer.data(idx_dip_hh + 179);

    auto tr_x_xxyzz_xyyzz = pbuffer.data(idx_dip_hh + 180);

    auto tr_x_xxyzz_xyzzz = pbuffer.data(idx_dip_hh + 181);

    auto tr_x_xxyzz_xzzzz = pbuffer.data(idx_dip_hh + 182);

    auto tr_x_xxzzz_xxxxx = pbuffer.data(idx_dip_hh + 189);

    auto tr_x_xxzzz_xxxxy = pbuffer.data(idx_dip_hh + 190);

    auto tr_x_xxzzz_xxxxz = pbuffer.data(idx_dip_hh + 191);

    auto tr_x_xxzzz_xxxyy = pbuffer.data(idx_dip_hh + 192);

    auto tr_x_xxzzz_xxxyz = pbuffer.data(idx_dip_hh + 193);

    auto tr_x_xxzzz_xxxzz = pbuffer.data(idx_dip_hh + 194);

    auto tr_x_xxzzz_xxyyy = pbuffer.data(idx_dip_hh + 195);

    auto tr_x_xxzzz_xxyyz = pbuffer.data(idx_dip_hh + 196);

    auto tr_x_xxzzz_xxyzz = pbuffer.data(idx_dip_hh + 197);

    auto tr_x_xxzzz_xxzzz = pbuffer.data(idx_dip_hh + 198);

    auto tr_x_xxzzz_xyyyy = pbuffer.data(idx_dip_hh + 199);

    auto tr_x_xxzzz_xyyyz = pbuffer.data(idx_dip_hh + 200);

    auto tr_x_xxzzz_xyyzz = pbuffer.data(idx_dip_hh + 201);

    auto tr_x_xxzzz_xyzzz = pbuffer.data(idx_dip_hh + 202);

    auto tr_x_xxzzz_xzzzz = pbuffer.data(idx_dip_hh + 203);

    auto tr_x_xxzzz_yyyyy = pbuffer.data(idx_dip_hh + 204);

    auto tr_x_xxzzz_yyyyz = pbuffer.data(idx_dip_hh + 205);

    auto tr_x_xxzzz_yyyzz = pbuffer.data(idx_dip_hh + 206);

    auto tr_x_xxzzz_yyzzz = pbuffer.data(idx_dip_hh + 207);

    auto tr_x_xxzzz_yzzzz = pbuffer.data(idx_dip_hh + 208);

    auto tr_x_xxzzz_zzzzz = pbuffer.data(idx_dip_hh + 209);

    auto tr_x_xyyyy_xxxxy = pbuffer.data(idx_dip_hh + 211);

    auto tr_x_xyyyy_xxxyy = pbuffer.data(idx_dip_hh + 213);

    auto tr_x_xyyyy_xxxyz = pbuffer.data(idx_dip_hh + 214);

    auto tr_x_xyyyy_xxyyy = pbuffer.data(idx_dip_hh + 216);

    auto tr_x_xyyyy_xxyyz = pbuffer.data(idx_dip_hh + 217);

    auto tr_x_xyyyy_xxyzz = pbuffer.data(idx_dip_hh + 218);

    auto tr_x_xyyyy_xyyyy = pbuffer.data(idx_dip_hh + 220);

    auto tr_x_xyyyy_xyyyz = pbuffer.data(idx_dip_hh + 221);

    auto tr_x_xyyyy_xyyzz = pbuffer.data(idx_dip_hh + 222);

    auto tr_x_xyyyy_xyzzz = pbuffer.data(idx_dip_hh + 223);

    auto tr_x_xzzzz_xxxxx = pbuffer.data(idx_dip_hh + 294);

    auto tr_x_xzzzz_xxxxy = pbuffer.data(idx_dip_hh + 295);

    auto tr_x_xzzzz_xxxxz = pbuffer.data(idx_dip_hh + 296);

    auto tr_x_xzzzz_xxxyy = pbuffer.data(idx_dip_hh + 297);

    auto tr_x_xzzzz_xxxyz = pbuffer.data(idx_dip_hh + 298);

    auto tr_x_xzzzz_xxxzz = pbuffer.data(idx_dip_hh + 299);

    auto tr_x_xzzzz_xxyyy = pbuffer.data(idx_dip_hh + 300);

    auto tr_x_xzzzz_xxyyz = pbuffer.data(idx_dip_hh + 301);

    auto tr_x_xzzzz_xxyzz = pbuffer.data(idx_dip_hh + 302);

    auto tr_x_xzzzz_xxzzz = pbuffer.data(idx_dip_hh + 303);

    auto tr_x_xzzzz_xyyyy = pbuffer.data(idx_dip_hh + 304);

    auto tr_x_xzzzz_xyyyz = pbuffer.data(idx_dip_hh + 305);

    auto tr_x_xzzzz_xyyzz = pbuffer.data(idx_dip_hh + 306);

    auto tr_x_xzzzz_xyzzz = pbuffer.data(idx_dip_hh + 307);

    auto tr_x_xzzzz_xzzzz = pbuffer.data(idx_dip_hh + 308);

    auto tr_x_yyyyy_xxxxx = pbuffer.data(idx_dip_hh + 315);

    auto tr_x_yyyyy_xxxxy = pbuffer.data(idx_dip_hh + 316);

    auto tr_x_yyyyy_xxxxz = pbuffer.data(idx_dip_hh + 317);

    auto tr_x_yyyyy_xxxyy = pbuffer.data(idx_dip_hh + 318);

    auto tr_x_yyyyy_xxxyz = pbuffer.data(idx_dip_hh + 319);

    auto tr_x_yyyyy_xxxzz = pbuffer.data(idx_dip_hh + 320);

    auto tr_x_yyyyy_xxyyy = pbuffer.data(idx_dip_hh + 321);

    auto tr_x_yyyyy_xxyyz = pbuffer.data(idx_dip_hh + 322);

    auto tr_x_yyyyy_xxyzz = pbuffer.data(idx_dip_hh + 323);

    auto tr_x_yyyyy_xxzzz = pbuffer.data(idx_dip_hh + 324);

    auto tr_x_yyyyy_xyyyy = pbuffer.data(idx_dip_hh + 325);

    auto tr_x_yyyyy_xyyyz = pbuffer.data(idx_dip_hh + 326);

    auto tr_x_yyyyy_xyyzz = pbuffer.data(idx_dip_hh + 327);

    auto tr_x_yyyyy_xyzzz = pbuffer.data(idx_dip_hh + 328);

    auto tr_x_yyyyy_xzzzz = pbuffer.data(idx_dip_hh + 329);

    auto tr_x_yyyyy_yyyyy = pbuffer.data(idx_dip_hh + 330);

    auto tr_x_yyyyy_yyyyz = pbuffer.data(idx_dip_hh + 331);

    auto tr_x_yyyyy_yyyzz = pbuffer.data(idx_dip_hh + 332);

    auto tr_x_yyyyy_yyzzz = pbuffer.data(idx_dip_hh + 333);

    auto tr_x_yyyyy_yzzzz = pbuffer.data(idx_dip_hh + 334);

    auto tr_x_yyyyy_zzzzz = pbuffer.data(idx_dip_hh + 335);

    auto tr_x_yyyzz_xxxxz = pbuffer.data(idx_dip_hh + 359);

    auto tr_x_yyyzz_xxxyz = pbuffer.data(idx_dip_hh + 361);

    auto tr_x_yyyzz_xxxzz = pbuffer.data(idx_dip_hh + 362);

    auto tr_x_yyyzz_xxyyz = pbuffer.data(idx_dip_hh + 364);

    auto tr_x_yyyzz_xxyzz = pbuffer.data(idx_dip_hh + 365);

    auto tr_x_yyyzz_xxzzz = pbuffer.data(idx_dip_hh + 366);

    auto tr_x_yyyzz_xyyyz = pbuffer.data(idx_dip_hh + 368);

    auto tr_x_yyyzz_xyyzz = pbuffer.data(idx_dip_hh + 369);

    auto tr_x_yyyzz_xyzzz = pbuffer.data(idx_dip_hh + 370);

    auto tr_x_yyyzz_xzzzz = pbuffer.data(idx_dip_hh + 371);

    auto tr_x_yyyzz_yyyyz = pbuffer.data(idx_dip_hh + 373);

    auto tr_x_yyyzz_yyyzz = pbuffer.data(idx_dip_hh + 374);

    auto tr_x_yyyzz_yyzzz = pbuffer.data(idx_dip_hh + 375);

    auto tr_x_yyyzz_yzzzz = pbuffer.data(idx_dip_hh + 376);

    auto tr_x_yyyzz_zzzzz = pbuffer.data(idx_dip_hh + 377);

    auto tr_x_yyzzz_xxxxz = pbuffer.data(idx_dip_hh + 380);

    auto tr_x_yyzzz_xxxyz = pbuffer.data(idx_dip_hh + 382);

    auto tr_x_yyzzz_xxxzz = pbuffer.data(idx_dip_hh + 383);

    auto tr_x_yyzzz_xxyyz = pbuffer.data(idx_dip_hh + 385);

    auto tr_x_yyzzz_xxyzz = pbuffer.data(idx_dip_hh + 386);

    auto tr_x_yyzzz_xxzzz = pbuffer.data(idx_dip_hh + 387);

    auto tr_x_yyzzz_xyyyz = pbuffer.data(idx_dip_hh + 389);

    auto tr_x_yyzzz_xyyzz = pbuffer.data(idx_dip_hh + 390);

    auto tr_x_yyzzz_xyzzz = pbuffer.data(idx_dip_hh + 391);

    auto tr_x_yyzzz_xzzzz = pbuffer.data(idx_dip_hh + 392);

    auto tr_x_yyzzz_yyyyz = pbuffer.data(idx_dip_hh + 394);

    auto tr_x_yyzzz_yyyzz = pbuffer.data(idx_dip_hh + 395);

    auto tr_x_yyzzz_yyzzz = pbuffer.data(idx_dip_hh + 396);

    auto tr_x_yyzzz_yzzzz = pbuffer.data(idx_dip_hh + 397);

    auto tr_x_yyzzz_zzzzz = pbuffer.data(idx_dip_hh + 398);

    auto tr_x_yzzzz_xxxxz = pbuffer.data(idx_dip_hh + 401);

    auto tr_x_yzzzz_xxxyz = pbuffer.data(idx_dip_hh + 403);

    auto tr_x_yzzzz_xxxzz = pbuffer.data(idx_dip_hh + 404);

    auto tr_x_yzzzz_xxyyz = pbuffer.data(idx_dip_hh + 406);

    auto tr_x_yzzzz_xxyzz = pbuffer.data(idx_dip_hh + 407);

    auto tr_x_yzzzz_xxzzz = pbuffer.data(idx_dip_hh + 408);

    auto tr_x_yzzzz_xyyyz = pbuffer.data(idx_dip_hh + 410);

    auto tr_x_yzzzz_xyyzz = pbuffer.data(idx_dip_hh + 411);

    auto tr_x_yzzzz_xyzzz = pbuffer.data(idx_dip_hh + 412);

    auto tr_x_yzzzz_xzzzz = pbuffer.data(idx_dip_hh + 413);

    auto tr_x_yzzzz_yyyyz = pbuffer.data(idx_dip_hh + 415);

    auto tr_x_yzzzz_yyyzz = pbuffer.data(idx_dip_hh + 416);

    auto tr_x_yzzzz_yyzzz = pbuffer.data(idx_dip_hh + 417);

    auto tr_x_yzzzz_yzzzz = pbuffer.data(idx_dip_hh + 418);

    auto tr_x_yzzzz_zzzzz = pbuffer.data(idx_dip_hh + 419);

    auto tr_x_zzzzz_xxxxx = pbuffer.data(idx_dip_hh + 420);

    auto tr_x_zzzzz_xxxxy = pbuffer.data(idx_dip_hh + 421);

    auto tr_x_zzzzz_xxxxz = pbuffer.data(idx_dip_hh + 422);

    auto tr_x_zzzzz_xxxyy = pbuffer.data(idx_dip_hh + 423);

    auto tr_x_zzzzz_xxxyz = pbuffer.data(idx_dip_hh + 424);

    auto tr_x_zzzzz_xxxzz = pbuffer.data(idx_dip_hh + 425);

    auto tr_x_zzzzz_xxyyy = pbuffer.data(idx_dip_hh + 426);

    auto tr_x_zzzzz_xxyyz = pbuffer.data(idx_dip_hh + 427);

    auto tr_x_zzzzz_xxyzz = pbuffer.data(idx_dip_hh + 428);

    auto tr_x_zzzzz_xxzzz = pbuffer.data(idx_dip_hh + 429);

    auto tr_x_zzzzz_xyyyy = pbuffer.data(idx_dip_hh + 430);

    auto tr_x_zzzzz_xyyyz = pbuffer.data(idx_dip_hh + 431);

    auto tr_x_zzzzz_xyyzz = pbuffer.data(idx_dip_hh + 432);

    auto tr_x_zzzzz_xyzzz = pbuffer.data(idx_dip_hh + 433);

    auto tr_x_zzzzz_xzzzz = pbuffer.data(idx_dip_hh + 434);

    auto tr_x_zzzzz_yyyyy = pbuffer.data(idx_dip_hh + 435);

    auto tr_x_zzzzz_yyyyz = pbuffer.data(idx_dip_hh + 436);

    auto tr_x_zzzzz_yyyzz = pbuffer.data(idx_dip_hh + 437);

    auto tr_x_zzzzz_yyzzz = pbuffer.data(idx_dip_hh + 438);

    auto tr_x_zzzzz_yzzzz = pbuffer.data(idx_dip_hh + 439);

    auto tr_x_zzzzz_zzzzz = pbuffer.data(idx_dip_hh + 440);

    auto tr_y_xxxxx_xxxxx = pbuffer.data(idx_dip_hh + 441);

    auto tr_y_xxxxx_xxxxy = pbuffer.data(idx_dip_hh + 442);

    auto tr_y_xxxxx_xxxxz = pbuffer.data(idx_dip_hh + 443);

    auto tr_y_xxxxx_xxxyy = pbuffer.data(idx_dip_hh + 444);

    auto tr_y_xxxxx_xxxyz = pbuffer.data(idx_dip_hh + 445);

    auto tr_y_xxxxx_xxxzz = pbuffer.data(idx_dip_hh + 446);

    auto tr_y_xxxxx_xxyyy = pbuffer.data(idx_dip_hh + 447);

    auto tr_y_xxxxx_xxyyz = pbuffer.data(idx_dip_hh + 448);

    auto tr_y_xxxxx_xxyzz = pbuffer.data(idx_dip_hh + 449);

    auto tr_y_xxxxx_xxzzz = pbuffer.data(idx_dip_hh + 450);

    auto tr_y_xxxxx_xyyyy = pbuffer.data(idx_dip_hh + 451);

    auto tr_y_xxxxx_xyyyz = pbuffer.data(idx_dip_hh + 452);

    auto tr_y_xxxxx_xyyzz = pbuffer.data(idx_dip_hh + 453);

    auto tr_y_xxxxx_xyzzz = pbuffer.data(idx_dip_hh + 454);

    auto tr_y_xxxxx_xzzzz = pbuffer.data(idx_dip_hh + 455);

    auto tr_y_xxxxx_yyyyy = pbuffer.data(idx_dip_hh + 456);

    auto tr_y_xxxxx_yyyyz = pbuffer.data(idx_dip_hh + 457);

    auto tr_y_xxxxx_yyyzz = pbuffer.data(idx_dip_hh + 458);

    auto tr_y_xxxxx_yyzzz = pbuffer.data(idx_dip_hh + 459);

    auto tr_y_xxxxx_yzzzz = pbuffer.data(idx_dip_hh + 460);

    auto tr_y_xxxxx_zzzzz = pbuffer.data(idx_dip_hh + 461);

    auto tr_y_xxxxy_xxxxy = pbuffer.data(idx_dip_hh + 463);

    auto tr_y_xxxxy_xxxyy = pbuffer.data(idx_dip_hh + 465);

    auto tr_y_xxxxy_xxxyz = pbuffer.data(idx_dip_hh + 466);

    auto tr_y_xxxxy_xxyyy = pbuffer.data(idx_dip_hh + 468);

    auto tr_y_xxxxy_xxyyz = pbuffer.data(idx_dip_hh + 469);

    auto tr_y_xxxxy_xxyzz = pbuffer.data(idx_dip_hh + 470);

    auto tr_y_xxxxy_xyyyy = pbuffer.data(idx_dip_hh + 472);

    auto tr_y_xxxxy_xyyyz = pbuffer.data(idx_dip_hh + 473);

    auto tr_y_xxxxy_xyyzz = pbuffer.data(idx_dip_hh + 474);

    auto tr_y_xxxxy_xyzzz = pbuffer.data(idx_dip_hh + 475);

    auto tr_y_xxxxy_yyyyy = pbuffer.data(idx_dip_hh + 477);

    auto tr_y_xxxxy_yyyyz = pbuffer.data(idx_dip_hh + 478);

    auto tr_y_xxxxy_yyyzz = pbuffer.data(idx_dip_hh + 479);

    auto tr_y_xxxxy_yyzzz = pbuffer.data(idx_dip_hh + 480);

    auto tr_y_xxxxy_yzzzz = pbuffer.data(idx_dip_hh + 481);

    auto tr_y_xxxyy_xxxxx = pbuffer.data(idx_dip_hh + 504);

    auto tr_y_xxxyy_xxxxy = pbuffer.data(idx_dip_hh + 505);

    auto tr_y_xxxyy_xxxxz = pbuffer.data(idx_dip_hh + 506);

    auto tr_y_xxxyy_xxxyy = pbuffer.data(idx_dip_hh + 507);

    auto tr_y_xxxyy_xxxyz = pbuffer.data(idx_dip_hh + 508);

    auto tr_y_xxxyy_xxxzz = pbuffer.data(idx_dip_hh + 509);

    auto tr_y_xxxyy_xxyyy = pbuffer.data(idx_dip_hh + 510);

    auto tr_y_xxxyy_xxyyz = pbuffer.data(idx_dip_hh + 511);

    auto tr_y_xxxyy_xxyzz = pbuffer.data(idx_dip_hh + 512);

    auto tr_y_xxxyy_xxzzz = pbuffer.data(idx_dip_hh + 513);

    auto tr_y_xxxyy_xyyyy = pbuffer.data(idx_dip_hh + 514);

    auto tr_y_xxxyy_xyyyz = pbuffer.data(idx_dip_hh + 515);

    auto tr_y_xxxyy_xyyzz = pbuffer.data(idx_dip_hh + 516);

    auto tr_y_xxxyy_xyzzz = pbuffer.data(idx_dip_hh + 517);

    auto tr_y_xxxyy_xzzzz = pbuffer.data(idx_dip_hh + 518);

    auto tr_y_xxxyy_yyyyy = pbuffer.data(idx_dip_hh + 519);

    auto tr_y_xxxyy_yyyyz = pbuffer.data(idx_dip_hh + 520);

    auto tr_y_xxxyy_yyyzz = pbuffer.data(idx_dip_hh + 521);

    auto tr_y_xxxyy_yyzzz = pbuffer.data(idx_dip_hh + 522);

    auto tr_y_xxxyy_yzzzz = pbuffer.data(idx_dip_hh + 523);

    auto tr_y_xxxyy_zzzzz = pbuffer.data(idx_dip_hh + 524);

    auto tr_y_xxxzz_xxxxz = pbuffer.data(idx_dip_hh + 548);

    auto tr_y_xxxzz_xxxyz = pbuffer.data(idx_dip_hh + 550);

    auto tr_y_xxxzz_xxxzz = pbuffer.data(idx_dip_hh + 551);

    auto tr_y_xxxzz_xxyyz = pbuffer.data(idx_dip_hh + 553);

    auto tr_y_xxxzz_xxyzz = pbuffer.data(idx_dip_hh + 554);

    auto tr_y_xxxzz_xxzzz = pbuffer.data(idx_dip_hh + 555);

    auto tr_y_xxxzz_xyyyz = pbuffer.data(idx_dip_hh + 557);

    auto tr_y_xxxzz_xyyzz = pbuffer.data(idx_dip_hh + 558);

    auto tr_y_xxxzz_xyzzz = pbuffer.data(idx_dip_hh + 559);

    auto tr_y_xxxzz_xzzzz = pbuffer.data(idx_dip_hh + 560);

    auto tr_y_xxxzz_yyyyz = pbuffer.data(idx_dip_hh + 562);

    auto tr_y_xxxzz_yyyzz = pbuffer.data(idx_dip_hh + 563);

    auto tr_y_xxxzz_yyzzz = pbuffer.data(idx_dip_hh + 564);

    auto tr_y_xxxzz_yzzzz = pbuffer.data(idx_dip_hh + 565);

    auto tr_y_xxxzz_zzzzz = pbuffer.data(idx_dip_hh + 566);

    auto tr_y_xxyyy_xxxxx = pbuffer.data(idx_dip_hh + 567);

    auto tr_y_xxyyy_xxxxy = pbuffer.data(idx_dip_hh + 568);

    auto tr_y_xxyyy_xxxxz = pbuffer.data(idx_dip_hh + 569);

    auto tr_y_xxyyy_xxxyy = pbuffer.data(idx_dip_hh + 570);

    auto tr_y_xxyyy_xxxyz = pbuffer.data(idx_dip_hh + 571);

    auto tr_y_xxyyy_xxxzz = pbuffer.data(idx_dip_hh + 572);

    auto tr_y_xxyyy_xxyyy = pbuffer.data(idx_dip_hh + 573);

    auto tr_y_xxyyy_xxyyz = pbuffer.data(idx_dip_hh + 574);

    auto tr_y_xxyyy_xxyzz = pbuffer.data(idx_dip_hh + 575);

    auto tr_y_xxyyy_xxzzz = pbuffer.data(idx_dip_hh + 576);

    auto tr_y_xxyyy_xyyyy = pbuffer.data(idx_dip_hh + 577);

    auto tr_y_xxyyy_xyyyz = pbuffer.data(idx_dip_hh + 578);

    auto tr_y_xxyyy_xyyzz = pbuffer.data(idx_dip_hh + 579);

    auto tr_y_xxyyy_xyzzz = pbuffer.data(idx_dip_hh + 580);

    auto tr_y_xxyyy_xzzzz = pbuffer.data(idx_dip_hh + 581);

    auto tr_y_xxyyy_yyyyy = pbuffer.data(idx_dip_hh + 582);

    auto tr_y_xxyyy_yyyyz = pbuffer.data(idx_dip_hh + 583);

    auto tr_y_xxyyy_yyyzz = pbuffer.data(idx_dip_hh + 584);

    auto tr_y_xxyyy_yyzzz = pbuffer.data(idx_dip_hh + 585);

    auto tr_y_xxyyy_yzzzz = pbuffer.data(idx_dip_hh + 586);

    auto tr_y_xxyyy_zzzzz = pbuffer.data(idx_dip_hh + 587);

    auto tr_y_xxyzz_xxxyz = pbuffer.data(idx_dip_hh + 613);

    auto tr_y_xxyzz_xxyyz = pbuffer.data(idx_dip_hh + 616);

    auto tr_y_xxyzz_xxyzz = pbuffer.data(idx_dip_hh + 617);

    auto tr_y_xxyzz_xyyyz = pbuffer.data(idx_dip_hh + 620);

    auto tr_y_xxyzz_xyyzz = pbuffer.data(idx_dip_hh + 621);

    auto tr_y_xxyzz_xyzzz = pbuffer.data(idx_dip_hh + 622);

    auto tr_y_xxyzz_yyyyz = pbuffer.data(idx_dip_hh + 625);

    auto tr_y_xxyzz_yyyzz = pbuffer.data(idx_dip_hh + 626);

    auto tr_y_xxyzz_yyzzz = pbuffer.data(idx_dip_hh + 627);

    auto tr_y_xxyzz_yzzzz = pbuffer.data(idx_dip_hh + 628);

    auto tr_y_xxzzz_xxxxz = pbuffer.data(idx_dip_hh + 632);

    auto tr_y_xxzzz_xxxyz = pbuffer.data(idx_dip_hh + 634);

    auto tr_y_xxzzz_xxxzz = pbuffer.data(idx_dip_hh + 635);

    auto tr_y_xxzzz_xxyyz = pbuffer.data(idx_dip_hh + 637);

    auto tr_y_xxzzz_xxyzz = pbuffer.data(idx_dip_hh + 638);

    auto tr_y_xxzzz_xxzzz = pbuffer.data(idx_dip_hh + 639);

    auto tr_y_xxzzz_xyyyz = pbuffer.data(idx_dip_hh + 641);

    auto tr_y_xxzzz_xyyzz = pbuffer.data(idx_dip_hh + 642);

    auto tr_y_xxzzz_xyzzz = pbuffer.data(idx_dip_hh + 643);

    auto tr_y_xxzzz_xzzzz = pbuffer.data(idx_dip_hh + 644);

    auto tr_y_xxzzz_yyyyz = pbuffer.data(idx_dip_hh + 646);

    auto tr_y_xxzzz_yyyzz = pbuffer.data(idx_dip_hh + 647);

    auto tr_y_xxzzz_yyzzz = pbuffer.data(idx_dip_hh + 648);

    auto tr_y_xxzzz_yzzzz = pbuffer.data(idx_dip_hh + 649);

    auto tr_y_xxzzz_zzzzz = pbuffer.data(idx_dip_hh + 650);

    auto tr_y_xyyyy_xxxxx = pbuffer.data(idx_dip_hh + 651);

    auto tr_y_xyyyy_xxxxy = pbuffer.data(idx_dip_hh + 652);

    auto tr_y_xyyyy_xxxxz = pbuffer.data(idx_dip_hh + 653);

    auto tr_y_xyyyy_xxxyy = pbuffer.data(idx_dip_hh + 654);

    auto tr_y_xyyyy_xxxyz = pbuffer.data(idx_dip_hh + 655);

    auto tr_y_xyyyy_xxxzz = pbuffer.data(idx_dip_hh + 656);

    auto tr_y_xyyyy_xxyyy = pbuffer.data(idx_dip_hh + 657);

    auto tr_y_xyyyy_xxyyz = pbuffer.data(idx_dip_hh + 658);

    auto tr_y_xyyyy_xxyzz = pbuffer.data(idx_dip_hh + 659);

    auto tr_y_xyyyy_xxzzz = pbuffer.data(idx_dip_hh + 660);

    auto tr_y_xyyyy_xyyyy = pbuffer.data(idx_dip_hh + 661);

    auto tr_y_xyyyy_xyyyz = pbuffer.data(idx_dip_hh + 662);

    auto tr_y_xyyyy_xyyzz = pbuffer.data(idx_dip_hh + 663);

    auto tr_y_xyyyy_xyzzz = pbuffer.data(idx_dip_hh + 664);

    auto tr_y_xyyyy_xzzzz = pbuffer.data(idx_dip_hh + 665);

    auto tr_y_xyyyy_yyyyy = pbuffer.data(idx_dip_hh + 666);

    auto tr_y_xyyyy_yyyyz = pbuffer.data(idx_dip_hh + 667);

    auto tr_y_xyyyy_yyyzz = pbuffer.data(idx_dip_hh + 668);

    auto tr_y_xyyyy_yyzzz = pbuffer.data(idx_dip_hh + 669);

    auto tr_y_xyyyy_yzzzz = pbuffer.data(idx_dip_hh + 670);

    auto tr_y_xyyyy_zzzzz = pbuffer.data(idx_dip_hh + 671);

    auto tr_y_xyyzz_xxxxz = pbuffer.data(idx_dip_hh + 695);

    auto tr_y_xyyzz_xxxyz = pbuffer.data(idx_dip_hh + 697);

    auto tr_y_xyyzz_xxxzz = pbuffer.data(idx_dip_hh + 698);

    auto tr_y_xyyzz_xxyyz = pbuffer.data(idx_dip_hh + 700);

    auto tr_y_xyyzz_xxyzz = pbuffer.data(idx_dip_hh + 701);

    auto tr_y_xyyzz_xxzzz = pbuffer.data(idx_dip_hh + 702);

    auto tr_y_xyyzz_xyyyz = pbuffer.data(idx_dip_hh + 704);

    auto tr_y_xyyzz_xyyzz = pbuffer.data(idx_dip_hh + 705);

    auto tr_y_xyyzz_xyzzz = pbuffer.data(idx_dip_hh + 706);

    auto tr_y_xyyzz_xzzzz = pbuffer.data(idx_dip_hh + 707);

    auto tr_y_xyyzz_yyyyz = pbuffer.data(idx_dip_hh + 709);

    auto tr_y_xyyzz_yyyzz = pbuffer.data(idx_dip_hh + 710);

    auto tr_y_xyyzz_yyzzz = pbuffer.data(idx_dip_hh + 711);

    auto tr_y_xyyzz_yzzzz = pbuffer.data(idx_dip_hh + 712);

    auto tr_y_xyyzz_zzzzz = pbuffer.data(idx_dip_hh + 713);

    auto tr_y_xyzzz_xxxyz = pbuffer.data(idx_dip_hh + 718);

    auto tr_y_xyzzz_xxyyz = pbuffer.data(idx_dip_hh + 721);

    auto tr_y_xyzzz_xxyzz = pbuffer.data(idx_dip_hh + 722);

    auto tr_y_xyzzz_xyyyz = pbuffer.data(idx_dip_hh + 725);

    auto tr_y_xyzzz_xyyzz = pbuffer.data(idx_dip_hh + 726);

    auto tr_y_xyzzz_xyzzz = pbuffer.data(idx_dip_hh + 727);

    auto tr_y_xyzzz_yyyyz = pbuffer.data(idx_dip_hh + 730);

    auto tr_y_xyzzz_yyyzz = pbuffer.data(idx_dip_hh + 731);

    auto tr_y_xyzzz_yyzzz = pbuffer.data(idx_dip_hh + 732);

    auto tr_y_xyzzz_yzzzz = pbuffer.data(idx_dip_hh + 733);

    auto tr_y_xzzzz_xxxxz = pbuffer.data(idx_dip_hh + 737);

    auto tr_y_xzzzz_xxxyz = pbuffer.data(idx_dip_hh + 739);

    auto tr_y_xzzzz_xxxzz = pbuffer.data(idx_dip_hh + 740);

    auto tr_y_xzzzz_xxyyz = pbuffer.data(idx_dip_hh + 742);

    auto tr_y_xzzzz_xxyzz = pbuffer.data(idx_dip_hh + 743);

    auto tr_y_xzzzz_xxzzz = pbuffer.data(idx_dip_hh + 744);

    auto tr_y_xzzzz_xyyyz = pbuffer.data(idx_dip_hh + 746);

    auto tr_y_xzzzz_xyyzz = pbuffer.data(idx_dip_hh + 747);

    auto tr_y_xzzzz_xyzzz = pbuffer.data(idx_dip_hh + 748);

    auto tr_y_xzzzz_xzzzz = pbuffer.data(idx_dip_hh + 749);

    auto tr_y_xzzzz_yyyyz = pbuffer.data(idx_dip_hh + 751);

    auto tr_y_xzzzz_yyyzz = pbuffer.data(idx_dip_hh + 752);

    auto tr_y_xzzzz_yyzzz = pbuffer.data(idx_dip_hh + 753);

    auto tr_y_xzzzz_yzzzz = pbuffer.data(idx_dip_hh + 754);

    auto tr_y_xzzzz_zzzzz = pbuffer.data(idx_dip_hh + 755);

    auto tr_y_yyyyy_xxxxx = pbuffer.data(idx_dip_hh + 756);

    auto tr_y_yyyyy_xxxxy = pbuffer.data(idx_dip_hh + 757);

    auto tr_y_yyyyy_xxxxz = pbuffer.data(idx_dip_hh + 758);

    auto tr_y_yyyyy_xxxyy = pbuffer.data(idx_dip_hh + 759);

    auto tr_y_yyyyy_xxxyz = pbuffer.data(idx_dip_hh + 760);

    auto tr_y_yyyyy_xxxzz = pbuffer.data(idx_dip_hh + 761);

    auto tr_y_yyyyy_xxyyy = pbuffer.data(idx_dip_hh + 762);

    auto tr_y_yyyyy_xxyyz = pbuffer.data(idx_dip_hh + 763);

    auto tr_y_yyyyy_xxyzz = pbuffer.data(idx_dip_hh + 764);

    auto tr_y_yyyyy_xxzzz = pbuffer.data(idx_dip_hh + 765);

    auto tr_y_yyyyy_xyyyy = pbuffer.data(idx_dip_hh + 766);

    auto tr_y_yyyyy_xyyyz = pbuffer.data(idx_dip_hh + 767);

    auto tr_y_yyyyy_xyyzz = pbuffer.data(idx_dip_hh + 768);

    auto tr_y_yyyyy_xyzzz = pbuffer.data(idx_dip_hh + 769);

    auto tr_y_yyyyy_xzzzz = pbuffer.data(idx_dip_hh + 770);

    auto tr_y_yyyyy_yyyyy = pbuffer.data(idx_dip_hh + 771);

    auto tr_y_yyyyy_yyyyz = pbuffer.data(idx_dip_hh + 772);

    auto tr_y_yyyyy_yyyzz = pbuffer.data(idx_dip_hh + 773);

    auto tr_y_yyyyy_yyzzz = pbuffer.data(idx_dip_hh + 774);

    auto tr_y_yyyyy_yzzzz = pbuffer.data(idx_dip_hh + 775);

    auto tr_y_yyyyy_zzzzz = pbuffer.data(idx_dip_hh + 776);

    auto tr_y_yyyyz_xxxxy = pbuffer.data(idx_dip_hh + 778);

    auto tr_y_yyyyz_xxxxz = pbuffer.data(idx_dip_hh + 779);

    auto tr_y_yyyyz_xxxyy = pbuffer.data(idx_dip_hh + 780);

    auto tr_y_yyyyz_xxxyz = pbuffer.data(idx_dip_hh + 781);

    auto tr_y_yyyyz_xxxzz = pbuffer.data(idx_dip_hh + 782);

    auto tr_y_yyyyz_xxyyy = pbuffer.data(idx_dip_hh + 783);

    auto tr_y_yyyyz_xxyyz = pbuffer.data(idx_dip_hh + 784);

    auto tr_y_yyyyz_xxyzz = pbuffer.data(idx_dip_hh + 785);

    auto tr_y_yyyyz_xxzzz = pbuffer.data(idx_dip_hh + 786);

    auto tr_y_yyyyz_xyyyy = pbuffer.data(idx_dip_hh + 787);

    auto tr_y_yyyyz_xyyyz = pbuffer.data(idx_dip_hh + 788);

    auto tr_y_yyyyz_xyyzz = pbuffer.data(idx_dip_hh + 789);

    auto tr_y_yyyyz_xyzzz = pbuffer.data(idx_dip_hh + 790);

    auto tr_y_yyyyz_xzzzz = pbuffer.data(idx_dip_hh + 791);

    auto tr_y_yyyyz_yyyyy = pbuffer.data(idx_dip_hh + 792);

    auto tr_y_yyyyz_yyyyz = pbuffer.data(idx_dip_hh + 793);

    auto tr_y_yyyyz_yyyzz = pbuffer.data(idx_dip_hh + 794);

    auto tr_y_yyyyz_yyzzz = pbuffer.data(idx_dip_hh + 795);

    auto tr_y_yyyyz_yzzzz = pbuffer.data(idx_dip_hh + 796);

    auto tr_y_yyyyz_zzzzz = pbuffer.data(idx_dip_hh + 797);

    auto tr_y_yyyzz_xxxxx = pbuffer.data(idx_dip_hh + 798);

    auto tr_y_yyyzz_xxxxy = pbuffer.data(idx_dip_hh + 799);

    auto tr_y_yyyzz_xxxxz = pbuffer.data(idx_dip_hh + 800);

    auto tr_y_yyyzz_xxxyy = pbuffer.data(idx_dip_hh + 801);

    auto tr_y_yyyzz_xxxyz = pbuffer.data(idx_dip_hh + 802);

    auto tr_y_yyyzz_xxxzz = pbuffer.data(idx_dip_hh + 803);

    auto tr_y_yyyzz_xxyyy = pbuffer.data(idx_dip_hh + 804);

    auto tr_y_yyyzz_xxyyz = pbuffer.data(idx_dip_hh + 805);

    auto tr_y_yyyzz_xxyzz = pbuffer.data(idx_dip_hh + 806);

    auto tr_y_yyyzz_xxzzz = pbuffer.data(idx_dip_hh + 807);

    auto tr_y_yyyzz_xyyyy = pbuffer.data(idx_dip_hh + 808);

    auto tr_y_yyyzz_xyyyz = pbuffer.data(idx_dip_hh + 809);

    auto tr_y_yyyzz_xyyzz = pbuffer.data(idx_dip_hh + 810);

    auto tr_y_yyyzz_xyzzz = pbuffer.data(idx_dip_hh + 811);

    auto tr_y_yyyzz_xzzzz = pbuffer.data(idx_dip_hh + 812);

    auto tr_y_yyyzz_yyyyy = pbuffer.data(idx_dip_hh + 813);

    auto tr_y_yyyzz_yyyyz = pbuffer.data(idx_dip_hh + 814);

    auto tr_y_yyyzz_yyyzz = pbuffer.data(idx_dip_hh + 815);

    auto tr_y_yyyzz_yyzzz = pbuffer.data(idx_dip_hh + 816);

    auto tr_y_yyyzz_yzzzz = pbuffer.data(idx_dip_hh + 817);

    auto tr_y_yyyzz_zzzzz = pbuffer.data(idx_dip_hh + 818);

    auto tr_y_yyzzz_xxxxx = pbuffer.data(idx_dip_hh + 819);

    auto tr_y_yyzzz_xxxxy = pbuffer.data(idx_dip_hh + 820);

    auto tr_y_yyzzz_xxxxz = pbuffer.data(idx_dip_hh + 821);

    auto tr_y_yyzzz_xxxyy = pbuffer.data(idx_dip_hh + 822);

    auto tr_y_yyzzz_xxxyz = pbuffer.data(idx_dip_hh + 823);

    auto tr_y_yyzzz_xxxzz = pbuffer.data(idx_dip_hh + 824);

    auto tr_y_yyzzz_xxyyy = pbuffer.data(idx_dip_hh + 825);

    auto tr_y_yyzzz_xxyyz = pbuffer.data(idx_dip_hh + 826);

    auto tr_y_yyzzz_xxyzz = pbuffer.data(idx_dip_hh + 827);

    auto tr_y_yyzzz_xxzzz = pbuffer.data(idx_dip_hh + 828);

    auto tr_y_yyzzz_xyyyy = pbuffer.data(idx_dip_hh + 829);

    auto tr_y_yyzzz_xyyyz = pbuffer.data(idx_dip_hh + 830);

    auto tr_y_yyzzz_xyyzz = pbuffer.data(idx_dip_hh + 831);

    auto tr_y_yyzzz_xyzzz = pbuffer.data(idx_dip_hh + 832);

    auto tr_y_yyzzz_xzzzz = pbuffer.data(idx_dip_hh + 833);

    auto tr_y_yyzzz_yyyyy = pbuffer.data(idx_dip_hh + 834);

    auto tr_y_yyzzz_yyyyz = pbuffer.data(idx_dip_hh + 835);

    auto tr_y_yyzzz_yyyzz = pbuffer.data(idx_dip_hh + 836);

    auto tr_y_yyzzz_yyzzz = pbuffer.data(idx_dip_hh + 837);

    auto tr_y_yyzzz_yzzzz = pbuffer.data(idx_dip_hh + 838);

    auto tr_y_yyzzz_zzzzz = pbuffer.data(idx_dip_hh + 839);

    auto tr_y_yzzzz_xxxxx = pbuffer.data(idx_dip_hh + 840);

    auto tr_y_yzzzz_xxxxy = pbuffer.data(idx_dip_hh + 841);

    auto tr_y_yzzzz_xxxxz = pbuffer.data(idx_dip_hh + 842);

    auto tr_y_yzzzz_xxxyy = pbuffer.data(idx_dip_hh + 843);

    auto tr_y_yzzzz_xxxyz = pbuffer.data(idx_dip_hh + 844);

    auto tr_y_yzzzz_xxxzz = pbuffer.data(idx_dip_hh + 845);

    auto tr_y_yzzzz_xxyyy = pbuffer.data(idx_dip_hh + 846);

    auto tr_y_yzzzz_xxyyz = pbuffer.data(idx_dip_hh + 847);

    auto tr_y_yzzzz_xxyzz = pbuffer.data(idx_dip_hh + 848);

    auto tr_y_yzzzz_xxzzz = pbuffer.data(idx_dip_hh + 849);

    auto tr_y_yzzzz_xyyyy = pbuffer.data(idx_dip_hh + 850);

    auto tr_y_yzzzz_xyyyz = pbuffer.data(idx_dip_hh + 851);

    auto tr_y_yzzzz_xyyzz = pbuffer.data(idx_dip_hh + 852);

    auto tr_y_yzzzz_xyzzz = pbuffer.data(idx_dip_hh + 853);

    auto tr_y_yzzzz_xzzzz = pbuffer.data(idx_dip_hh + 854);

    auto tr_y_yzzzz_yyyyy = pbuffer.data(idx_dip_hh + 855);

    auto tr_y_yzzzz_yyyyz = pbuffer.data(idx_dip_hh + 856);

    auto tr_y_yzzzz_yyyzz = pbuffer.data(idx_dip_hh + 857);

    auto tr_y_yzzzz_yyzzz = pbuffer.data(idx_dip_hh + 858);

    auto tr_y_yzzzz_yzzzz = pbuffer.data(idx_dip_hh + 859);

    auto tr_y_yzzzz_zzzzz = pbuffer.data(idx_dip_hh + 860);

    auto tr_y_zzzzz_xxxxx = pbuffer.data(idx_dip_hh + 861);

    auto tr_y_zzzzz_xxxxy = pbuffer.data(idx_dip_hh + 862);

    auto tr_y_zzzzz_xxxxz = pbuffer.data(idx_dip_hh + 863);

    auto tr_y_zzzzz_xxxyy = pbuffer.data(idx_dip_hh + 864);

    auto tr_y_zzzzz_xxxyz = pbuffer.data(idx_dip_hh + 865);

    auto tr_y_zzzzz_xxxzz = pbuffer.data(idx_dip_hh + 866);

    auto tr_y_zzzzz_xxyyy = pbuffer.data(idx_dip_hh + 867);

    auto tr_y_zzzzz_xxyyz = pbuffer.data(idx_dip_hh + 868);

    auto tr_y_zzzzz_xxyzz = pbuffer.data(idx_dip_hh + 869);

    auto tr_y_zzzzz_xxzzz = pbuffer.data(idx_dip_hh + 870);

    auto tr_y_zzzzz_xyyyy = pbuffer.data(idx_dip_hh + 871);

    auto tr_y_zzzzz_xyyyz = pbuffer.data(idx_dip_hh + 872);

    auto tr_y_zzzzz_xyyzz = pbuffer.data(idx_dip_hh + 873);

    auto tr_y_zzzzz_xyzzz = pbuffer.data(idx_dip_hh + 874);

    auto tr_y_zzzzz_xzzzz = pbuffer.data(idx_dip_hh + 875);

    auto tr_y_zzzzz_yyyyy = pbuffer.data(idx_dip_hh + 876);

    auto tr_y_zzzzz_yyyyz = pbuffer.data(idx_dip_hh + 877);

    auto tr_y_zzzzz_yyyzz = pbuffer.data(idx_dip_hh + 878);

    auto tr_y_zzzzz_yyzzz = pbuffer.data(idx_dip_hh + 879);

    auto tr_y_zzzzz_yzzzz = pbuffer.data(idx_dip_hh + 880);

    auto tr_y_zzzzz_zzzzz = pbuffer.data(idx_dip_hh + 881);

    auto tr_z_xxxxx_xxxxx = pbuffer.data(idx_dip_hh + 882);

    auto tr_z_xxxxx_xxxxy = pbuffer.data(idx_dip_hh + 883);

    auto tr_z_xxxxx_xxxxz = pbuffer.data(idx_dip_hh + 884);

    auto tr_z_xxxxx_xxxyy = pbuffer.data(idx_dip_hh + 885);

    auto tr_z_xxxxx_xxxyz = pbuffer.data(idx_dip_hh + 886);

    auto tr_z_xxxxx_xxxzz = pbuffer.data(idx_dip_hh + 887);

    auto tr_z_xxxxx_xxyyy = pbuffer.data(idx_dip_hh + 888);

    auto tr_z_xxxxx_xxyyz = pbuffer.data(idx_dip_hh + 889);

    auto tr_z_xxxxx_xxyzz = pbuffer.data(idx_dip_hh + 890);

    auto tr_z_xxxxx_xxzzz = pbuffer.data(idx_dip_hh + 891);

    auto tr_z_xxxxx_xyyyy = pbuffer.data(idx_dip_hh + 892);

    auto tr_z_xxxxx_xyyyz = pbuffer.data(idx_dip_hh + 893);

    auto tr_z_xxxxx_xyyzz = pbuffer.data(idx_dip_hh + 894);

    auto tr_z_xxxxx_xyzzz = pbuffer.data(idx_dip_hh + 895);

    auto tr_z_xxxxx_xzzzz = pbuffer.data(idx_dip_hh + 896);

    auto tr_z_xxxxx_yyyyy = pbuffer.data(idx_dip_hh + 897);

    auto tr_z_xxxxx_yyyyz = pbuffer.data(idx_dip_hh + 898);

    auto tr_z_xxxxx_yyyzz = pbuffer.data(idx_dip_hh + 899);

    auto tr_z_xxxxx_yyzzz = pbuffer.data(idx_dip_hh + 900);

    auto tr_z_xxxxx_yzzzz = pbuffer.data(idx_dip_hh + 901);

    auto tr_z_xxxxx_zzzzz = pbuffer.data(idx_dip_hh + 902);

    auto tr_z_xxxxz_xxxxx = pbuffer.data(idx_dip_hh + 924);

    auto tr_z_xxxxz_xxxxy = pbuffer.data(idx_dip_hh + 925);

    auto tr_z_xxxxz_xxxxz = pbuffer.data(idx_dip_hh + 926);

    auto tr_z_xxxxz_xxxyy = pbuffer.data(idx_dip_hh + 927);

    auto tr_z_xxxxz_xxxyz = pbuffer.data(idx_dip_hh + 928);

    auto tr_z_xxxxz_xxxzz = pbuffer.data(idx_dip_hh + 929);

    auto tr_z_xxxxz_xxyyy = pbuffer.data(idx_dip_hh + 930);

    auto tr_z_xxxxz_xxyyz = pbuffer.data(idx_dip_hh + 931);

    auto tr_z_xxxxz_xxyzz = pbuffer.data(idx_dip_hh + 932);

    auto tr_z_xxxxz_xxzzz = pbuffer.data(idx_dip_hh + 933);

    auto tr_z_xxxxz_xyyyy = pbuffer.data(idx_dip_hh + 934);

    auto tr_z_xxxxz_xyyyz = pbuffer.data(idx_dip_hh + 935);

    auto tr_z_xxxxz_xyyzz = pbuffer.data(idx_dip_hh + 936);

    auto tr_z_xxxxz_xyzzz = pbuffer.data(idx_dip_hh + 937);

    auto tr_z_xxxxz_xzzzz = pbuffer.data(idx_dip_hh + 938);

    auto tr_z_xxxxz_yyyyz = pbuffer.data(idx_dip_hh + 940);

    auto tr_z_xxxxz_yyyzz = pbuffer.data(idx_dip_hh + 941);

    auto tr_z_xxxxz_yyzzz = pbuffer.data(idx_dip_hh + 942);

    auto tr_z_xxxxz_yzzzz = pbuffer.data(idx_dip_hh + 943);

    auto tr_z_xxxxz_zzzzz = pbuffer.data(idx_dip_hh + 944);

    auto tr_z_xxxyy_xxxxy = pbuffer.data(idx_dip_hh + 946);

    auto tr_z_xxxyy_xxxyy = pbuffer.data(idx_dip_hh + 948);

    auto tr_z_xxxyy_xxxyz = pbuffer.data(idx_dip_hh + 949);

    auto tr_z_xxxyy_xxyyy = pbuffer.data(idx_dip_hh + 951);

    auto tr_z_xxxyy_xxyyz = pbuffer.data(idx_dip_hh + 952);

    auto tr_z_xxxyy_xxyzz = pbuffer.data(idx_dip_hh + 953);

    auto tr_z_xxxyy_xyyyy = pbuffer.data(idx_dip_hh + 955);

    auto tr_z_xxxyy_xyyyz = pbuffer.data(idx_dip_hh + 956);

    auto tr_z_xxxyy_xyyzz = pbuffer.data(idx_dip_hh + 957);

    auto tr_z_xxxyy_xyzzz = pbuffer.data(idx_dip_hh + 958);

    auto tr_z_xxxyy_yyyyy = pbuffer.data(idx_dip_hh + 960);

    auto tr_z_xxxyy_yyyyz = pbuffer.data(idx_dip_hh + 961);

    auto tr_z_xxxyy_yyyzz = pbuffer.data(idx_dip_hh + 962);

    auto tr_z_xxxyy_yyzzz = pbuffer.data(idx_dip_hh + 963);

    auto tr_z_xxxyy_yzzzz = pbuffer.data(idx_dip_hh + 964);

    auto tr_z_xxxzz_xxxxx = pbuffer.data(idx_dip_hh + 987);

    auto tr_z_xxxzz_xxxxy = pbuffer.data(idx_dip_hh + 988);

    auto tr_z_xxxzz_xxxxz = pbuffer.data(idx_dip_hh + 989);

    auto tr_z_xxxzz_xxxyy = pbuffer.data(idx_dip_hh + 990);

    auto tr_z_xxxzz_xxxyz = pbuffer.data(idx_dip_hh + 991);

    auto tr_z_xxxzz_xxxzz = pbuffer.data(idx_dip_hh + 992);

    auto tr_z_xxxzz_xxyyy = pbuffer.data(idx_dip_hh + 993);

    auto tr_z_xxxzz_xxyyz = pbuffer.data(idx_dip_hh + 994);

    auto tr_z_xxxzz_xxyzz = pbuffer.data(idx_dip_hh + 995);

    auto tr_z_xxxzz_xxzzz = pbuffer.data(idx_dip_hh + 996);

    auto tr_z_xxxzz_xyyyy = pbuffer.data(idx_dip_hh + 997);

    auto tr_z_xxxzz_xyyyz = pbuffer.data(idx_dip_hh + 998);

    auto tr_z_xxxzz_xyyzz = pbuffer.data(idx_dip_hh + 999);

    auto tr_z_xxxzz_xyzzz = pbuffer.data(idx_dip_hh + 1000);

    auto tr_z_xxxzz_xzzzz = pbuffer.data(idx_dip_hh + 1001);

    auto tr_z_xxxzz_yyyyy = pbuffer.data(idx_dip_hh + 1002);

    auto tr_z_xxxzz_yyyyz = pbuffer.data(idx_dip_hh + 1003);

    auto tr_z_xxxzz_yyyzz = pbuffer.data(idx_dip_hh + 1004);

    auto tr_z_xxxzz_yyzzz = pbuffer.data(idx_dip_hh + 1005);

    auto tr_z_xxxzz_yzzzz = pbuffer.data(idx_dip_hh + 1006);

    auto tr_z_xxxzz_zzzzz = pbuffer.data(idx_dip_hh + 1007);

    auto tr_z_xxyyy_xxxxy = pbuffer.data(idx_dip_hh + 1009);

    auto tr_z_xxyyy_xxxyy = pbuffer.data(idx_dip_hh + 1011);

    auto tr_z_xxyyy_xxxyz = pbuffer.data(idx_dip_hh + 1012);

    auto tr_z_xxyyy_xxyyy = pbuffer.data(idx_dip_hh + 1014);

    auto tr_z_xxyyy_xxyyz = pbuffer.data(idx_dip_hh + 1015);

    auto tr_z_xxyyy_xxyzz = pbuffer.data(idx_dip_hh + 1016);

    auto tr_z_xxyyy_xyyyy = pbuffer.data(idx_dip_hh + 1018);

    auto tr_z_xxyyy_xyyyz = pbuffer.data(idx_dip_hh + 1019);

    auto tr_z_xxyyy_xyyzz = pbuffer.data(idx_dip_hh + 1020);

    auto tr_z_xxyyy_xyzzz = pbuffer.data(idx_dip_hh + 1021);

    auto tr_z_xxyyy_yyyyy = pbuffer.data(idx_dip_hh + 1023);

    auto tr_z_xxyyy_yyyyz = pbuffer.data(idx_dip_hh + 1024);

    auto tr_z_xxyyy_yyyzz = pbuffer.data(idx_dip_hh + 1025);

    auto tr_z_xxyyy_yyzzz = pbuffer.data(idx_dip_hh + 1026);

    auto tr_z_xxyyy_yzzzz = pbuffer.data(idx_dip_hh + 1027);

    auto tr_z_xxyyz_xxxyz = pbuffer.data(idx_dip_hh + 1033);

    auto tr_z_xxyyz_xxyyz = pbuffer.data(idx_dip_hh + 1036);

    auto tr_z_xxyyz_xxyzz = pbuffer.data(idx_dip_hh + 1037);

    auto tr_z_xxyyz_xyyyz = pbuffer.data(idx_dip_hh + 1040);

    auto tr_z_xxyyz_xyyzz = pbuffer.data(idx_dip_hh + 1041);

    auto tr_z_xxyyz_xyzzz = pbuffer.data(idx_dip_hh + 1042);

    auto tr_z_xxyyz_yyyyz = pbuffer.data(idx_dip_hh + 1045);

    auto tr_z_xxyyz_yyyzz = pbuffer.data(idx_dip_hh + 1046);

    auto tr_z_xxyyz_yyzzz = pbuffer.data(idx_dip_hh + 1047);

    auto tr_z_xxyyz_yzzzz = pbuffer.data(idx_dip_hh + 1048);

    auto tr_z_xxzzz_xxxxx = pbuffer.data(idx_dip_hh + 1071);

    auto tr_z_xxzzz_xxxxy = pbuffer.data(idx_dip_hh + 1072);

    auto tr_z_xxzzz_xxxxz = pbuffer.data(idx_dip_hh + 1073);

    auto tr_z_xxzzz_xxxyy = pbuffer.data(idx_dip_hh + 1074);

    auto tr_z_xxzzz_xxxyz = pbuffer.data(idx_dip_hh + 1075);

    auto tr_z_xxzzz_xxxzz = pbuffer.data(idx_dip_hh + 1076);

    auto tr_z_xxzzz_xxyyy = pbuffer.data(idx_dip_hh + 1077);

    auto tr_z_xxzzz_xxyyz = pbuffer.data(idx_dip_hh + 1078);

    auto tr_z_xxzzz_xxyzz = pbuffer.data(idx_dip_hh + 1079);

    auto tr_z_xxzzz_xxzzz = pbuffer.data(idx_dip_hh + 1080);

    auto tr_z_xxzzz_xyyyy = pbuffer.data(idx_dip_hh + 1081);

    auto tr_z_xxzzz_xyyyz = pbuffer.data(idx_dip_hh + 1082);

    auto tr_z_xxzzz_xyyzz = pbuffer.data(idx_dip_hh + 1083);

    auto tr_z_xxzzz_xyzzz = pbuffer.data(idx_dip_hh + 1084);

    auto tr_z_xxzzz_xzzzz = pbuffer.data(idx_dip_hh + 1085);

    auto tr_z_xxzzz_yyyyy = pbuffer.data(idx_dip_hh + 1086);

    auto tr_z_xxzzz_yyyyz = pbuffer.data(idx_dip_hh + 1087);

    auto tr_z_xxzzz_yyyzz = pbuffer.data(idx_dip_hh + 1088);

    auto tr_z_xxzzz_yyzzz = pbuffer.data(idx_dip_hh + 1089);

    auto tr_z_xxzzz_yzzzz = pbuffer.data(idx_dip_hh + 1090);

    auto tr_z_xxzzz_zzzzz = pbuffer.data(idx_dip_hh + 1091);

    auto tr_z_xyyyy_xxxxy = pbuffer.data(idx_dip_hh + 1093);

    auto tr_z_xyyyy_xxxyy = pbuffer.data(idx_dip_hh + 1095);

    auto tr_z_xyyyy_xxxyz = pbuffer.data(idx_dip_hh + 1096);

    auto tr_z_xyyyy_xxyyy = pbuffer.data(idx_dip_hh + 1098);

    auto tr_z_xyyyy_xxyyz = pbuffer.data(idx_dip_hh + 1099);

    auto tr_z_xyyyy_xxyzz = pbuffer.data(idx_dip_hh + 1100);

    auto tr_z_xyyyy_xyyyy = pbuffer.data(idx_dip_hh + 1102);

    auto tr_z_xyyyy_xyyyz = pbuffer.data(idx_dip_hh + 1103);

    auto tr_z_xyyyy_xyyzz = pbuffer.data(idx_dip_hh + 1104);

    auto tr_z_xyyyy_xyzzz = pbuffer.data(idx_dip_hh + 1105);

    auto tr_z_xyyyy_yyyyy = pbuffer.data(idx_dip_hh + 1107);

    auto tr_z_xyyyy_yyyyz = pbuffer.data(idx_dip_hh + 1108);

    auto tr_z_xyyyy_yyyzz = pbuffer.data(idx_dip_hh + 1109);

    auto tr_z_xyyyy_yyzzz = pbuffer.data(idx_dip_hh + 1110);

    auto tr_z_xyyyy_yzzzz = pbuffer.data(idx_dip_hh + 1111);

    auto tr_z_xyyyz_xxxyz = pbuffer.data(idx_dip_hh + 1117);

    auto tr_z_xyyyz_xxyyz = pbuffer.data(idx_dip_hh + 1120);

    auto tr_z_xyyyz_xxyzz = pbuffer.data(idx_dip_hh + 1121);

    auto tr_z_xyyyz_xyyyz = pbuffer.data(idx_dip_hh + 1124);

    auto tr_z_xyyyz_xyyzz = pbuffer.data(idx_dip_hh + 1125);

    auto tr_z_xyyyz_xyzzz = pbuffer.data(idx_dip_hh + 1126);

    auto tr_z_xyyyz_yyyyz = pbuffer.data(idx_dip_hh + 1129);

    auto tr_z_xyyyz_yyyzz = pbuffer.data(idx_dip_hh + 1130);

    auto tr_z_xyyyz_yyzzz = pbuffer.data(idx_dip_hh + 1131);

    auto tr_z_xyyyz_yzzzz = pbuffer.data(idx_dip_hh + 1132);

    auto tr_z_xyyzz_xxxxy = pbuffer.data(idx_dip_hh + 1135);

    auto tr_z_xyyzz_xxxyy = pbuffer.data(idx_dip_hh + 1137);

    auto tr_z_xyyzz_xxxyz = pbuffer.data(idx_dip_hh + 1138);

    auto tr_z_xyyzz_xxyyy = pbuffer.data(idx_dip_hh + 1140);

    auto tr_z_xyyzz_xxyyz = pbuffer.data(idx_dip_hh + 1141);

    auto tr_z_xyyzz_xxyzz = pbuffer.data(idx_dip_hh + 1142);

    auto tr_z_xyyzz_xyyyy = pbuffer.data(idx_dip_hh + 1144);

    auto tr_z_xyyzz_xyyyz = pbuffer.data(idx_dip_hh + 1145);

    auto tr_z_xyyzz_xyyzz = pbuffer.data(idx_dip_hh + 1146);

    auto tr_z_xyyzz_xyzzz = pbuffer.data(idx_dip_hh + 1147);

    auto tr_z_xyyzz_yyyyy = pbuffer.data(idx_dip_hh + 1149);

    auto tr_z_xyyzz_yyyyz = pbuffer.data(idx_dip_hh + 1150);

    auto tr_z_xyyzz_yyyzz = pbuffer.data(idx_dip_hh + 1151);

    auto tr_z_xyyzz_yyzzz = pbuffer.data(idx_dip_hh + 1152);

    auto tr_z_xyyzz_yzzzz = pbuffer.data(idx_dip_hh + 1153);

    auto tr_z_xzzzz_xxxxx = pbuffer.data(idx_dip_hh + 1176);

    auto tr_z_xzzzz_xxxxy = pbuffer.data(idx_dip_hh + 1177);

    auto tr_z_xzzzz_xxxxz = pbuffer.data(idx_dip_hh + 1178);

    auto tr_z_xzzzz_xxxyy = pbuffer.data(idx_dip_hh + 1179);

    auto tr_z_xzzzz_xxxyz = pbuffer.data(idx_dip_hh + 1180);

    auto tr_z_xzzzz_xxxzz = pbuffer.data(idx_dip_hh + 1181);

    auto tr_z_xzzzz_xxyyy = pbuffer.data(idx_dip_hh + 1182);

    auto tr_z_xzzzz_xxyyz = pbuffer.data(idx_dip_hh + 1183);

    auto tr_z_xzzzz_xxyzz = pbuffer.data(idx_dip_hh + 1184);

    auto tr_z_xzzzz_xxzzz = pbuffer.data(idx_dip_hh + 1185);

    auto tr_z_xzzzz_xyyyy = pbuffer.data(idx_dip_hh + 1186);

    auto tr_z_xzzzz_xyyyz = pbuffer.data(idx_dip_hh + 1187);

    auto tr_z_xzzzz_xyyzz = pbuffer.data(idx_dip_hh + 1188);

    auto tr_z_xzzzz_xyzzz = pbuffer.data(idx_dip_hh + 1189);

    auto tr_z_xzzzz_xzzzz = pbuffer.data(idx_dip_hh + 1190);

    auto tr_z_xzzzz_yyyyy = pbuffer.data(idx_dip_hh + 1191);

    auto tr_z_xzzzz_yyyyz = pbuffer.data(idx_dip_hh + 1192);

    auto tr_z_xzzzz_yyyzz = pbuffer.data(idx_dip_hh + 1193);

    auto tr_z_xzzzz_yyzzz = pbuffer.data(idx_dip_hh + 1194);

    auto tr_z_xzzzz_yzzzz = pbuffer.data(idx_dip_hh + 1195);

    auto tr_z_xzzzz_zzzzz = pbuffer.data(idx_dip_hh + 1196);

    auto tr_z_yyyyy_xxxxx = pbuffer.data(idx_dip_hh + 1197);

    auto tr_z_yyyyy_xxxxy = pbuffer.data(idx_dip_hh + 1198);

    auto tr_z_yyyyy_xxxxz = pbuffer.data(idx_dip_hh + 1199);

    auto tr_z_yyyyy_xxxyy = pbuffer.data(idx_dip_hh + 1200);

    auto tr_z_yyyyy_xxxyz = pbuffer.data(idx_dip_hh + 1201);

    auto tr_z_yyyyy_xxxzz = pbuffer.data(idx_dip_hh + 1202);

    auto tr_z_yyyyy_xxyyy = pbuffer.data(idx_dip_hh + 1203);

    auto tr_z_yyyyy_xxyyz = pbuffer.data(idx_dip_hh + 1204);

    auto tr_z_yyyyy_xxyzz = pbuffer.data(idx_dip_hh + 1205);

    auto tr_z_yyyyy_xxzzz = pbuffer.data(idx_dip_hh + 1206);

    auto tr_z_yyyyy_xyyyy = pbuffer.data(idx_dip_hh + 1207);

    auto tr_z_yyyyy_xyyyz = pbuffer.data(idx_dip_hh + 1208);

    auto tr_z_yyyyy_xyyzz = pbuffer.data(idx_dip_hh + 1209);

    auto tr_z_yyyyy_xyzzz = pbuffer.data(idx_dip_hh + 1210);

    auto tr_z_yyyyy_xzzzz = pbuffer.data(idx_dip_hh + 1211);

    auto tr_z_yyyyy_yyyyy = pbuffer.data(idx_dip_hh + 1212);

    auto tr_z_yyyyy_yyyyz = pbuffer.data(idx_dip_hh + 1213);

    auto tr_z_yyyyy_yyyzz = pbuffer.data(idx_dip_hh + 1214);

    auto tr_z_yyyyy_yyzzz = pbuffer.data(idx_dip_hh + 1215);

    auto tr_z_yyyyy_yzzzz = pbuffer.data(idx_dip_hh + 1216);

    auto tr_z_yyyyy_zzzzz = pbuffer.data(idx_dip_hh + 1217);

    auto tr_z_yyyyz_xxxxx = pbuffer.data(idx_dip_hh + 1218);

    auto tr_z_yyyyz_xxxxy = pbuffer.data(idx_dip_hh + 1219);

    auto tr_z_yyyyz_xxxxz = pbuffer.data(idx_dip_hh + 1220);

    auto tr_z_yyyyz_xxxyy = pbuffer.data(idx_dip_hh + 1221);

    auto tr_z_yyyyz_xxxyz = pbuffer.data(idx_dip_hh + 1222);

    auto tr_z_yyyyz_xxxzz = pbuffer.data(idx_dip_hh + 1223);

    auto tr_z_yyyyz_xxyyy = pbuffer.data(idx_dip_hh + 1224);

    auto tr_z_yyyyz_xxyyz = pbuffer.data(idx_dip_hh + 1225);

    auto tr_z_yyyyz_xxyzz = pbuffer.data(idx_dip_hh + 1226);

    auto tr_z_yyyyz_xxzzz = pbuffer.data(idx_dip_hh + 1227);

    auto tr_z_yyyyz_xyyyy = pbuffer.data(idx_dip_hh + 1228);

    auto tr_z_yyyyz_xyyyz = pbuffer.data(idx_dip_hh + 1229);

    auto tr_z_yyyyz_xyyzz = pbuffer.data(idx_dip_hh + 1230);

    auto tr_z_yyyyz_xyzzz = pbuffer.data(idx_dip_hh + 1231);

    auto tr_z_yyyyz_xzzzz = pbuffer.data(idx_dip_hh + 1232);

    auto tr_z_yyyyz_yyyyy = pbuffer.data(idx_dip_hh + 1233);

    auto tr_z_yyyyz_yyyyz = pbuffer.data(idx_dip_hh + 1234);

    auto tr_z_yyyyz_yyyzz = pbuffer.data(idx_dip_hh + 1235);

    auto tr_z_yyyyz_yyzzz = pbuffer.data(idx_dip_hh + 1236);

    auto tr_z_yyyyz_yzzzz = pbuffer.data(idx_dip_hh + 1237);

    auto tr_z_yyyyz_zzzzz = pbuffer.data(idx_dip_hh + 1238);

    auto tr_z_yyyzz_xxxxx = pbuffer.data(idx_dip_hh + 1239);

    auto tr_z_yyyzz_xxxxy = pbuffer.data(idx_dip_hh + 1240);

    auto tr_z_yyyzz_xxxxz = pbuffer.data(idx_dip_hh + 1241);

    auto tr_z_yyyzz_xxxyy = pbuffer.data(idx_dip_hh + 1242);

    auto tr_z_yyyzz_xxxyz = pbuffer.data(idx_dip_hh + 1243);

    auto tr_z_yyyzz_xxxzz = pbuffer.data(idx_dip_hh + 1244);

    auto tr_z_yyyzz_xxyyy = pbuffer.data(idx_dip_hh + 1245);

    auto tr_z_yyyzz_xxyyz = pbuffer.data(idx_dip_hh + 1246);

    auto tr_z_yyyzz_xxyzz = pbuffer.data(idx_dip_hh + 1247);

    auto tr_z_yyyzz_xxzzz = pbuffer.data(idx_dip_hh + 1248);

    auto tr_z_yyyzz_xyyyy = pbuffer.data(idx_dip_hh + 1249);

    auto tr_z_yyyzz_xyyyz = pbuffer.data(idx_dip_hh + 1250);

    auto tr_z_yyyzz_xyyzz = pbuffer.data(idx_dip_hh + 1251);

    auto tr_z_yyyzz_xyzzz = pbuffer.data(idx_dip_hh + 1252);

    auto tr_z_yyyzz_xzzzz = pbuffer.data(idx_dip_hh + 1253);

    auto tr_z_yyyzz_yyyyy = pbuffer.data(idx_dip_hh + 1254);

    auto tr_z_yyyzz_yyyyz = pbuffer.data(idx_dip_hh + 1255);

    auto tr_z_yyyzz_yyyzz = pbuffer.data(idx_dip_hh + 1256);

    auto tr_z_yyyzz_yyzzz = pbuffer.data(idx_dip_hh + 1257);

    auto tr_z_yyyzz_yzzzz = pbuffer.data(idx_dip_hh + 1258);

    auto tr_z_yyyzz_zzzzz = pbuffer.data(idx_dip_hh + 1259);

    auto tr_z_yyzzz_xxxxx = pbuffer.data(idx_dip_hh + 1260);

    auto tr_z_yyzzz_xxxxy = pbuffer.data(idx_dip_hh + 1261);

    auto tr_z_yyzzz_xxxxz = pbuffer.data(idx_dip_hh + 1262);

    auto tr_z_yyzzz_xxxyy = pbuffer.data(idx_dip_hh + 1263);

    auto tr_z_yyzzz_xxxyz = pbuffer.data(idx_dip_hh + 1264);

    auto tr_z_yyzzz_xxxzz = pbuffer.data(idx_dip_hh + 1265);

    auto tr_z_yyzzz_xxyyy = pbuffer.data(idx_dip_hh + 1266);

    auto tr_z_yyzzz_xxyyz = pbuffer.data(idx_dip_hh + 1267);

    auto tr_z_yyzzz_xxyzz = pbuffer.data(idx_dip_hh + 1268);

    auto tr_z_yyzzz_xxzzz = pbuffer.data(idx_dip_hh + 1269);

    auto tr_z_yyzzz_xyyyy = pbuffer.data(idx_dip_hh + 1270);

    auto tr_z_yyzzz_xyyyz = pbuffer.data(idx_dip_hh + 1271);

    auto tr_z_yyzzz_xyyzz = pbuffer.data(idx_dip_hh + 1272);

    auto tr_z_yyzzz_xyzzz = pbuffer.data(idx_dip_hh + 1273);

    auto tr_z_yyzzz_xzzzz = pbuffer.data(idx_dip_hh + 1274);

    auto tr_z_yyzzz_yyyyy = pbuffer.data(idx_dip_hh + 1275);

    auto tr_z_yyzzz_yyyyz = pbuffer.data(idx_dip_hh + 1276);

    auto tr_z_yyzzz_yyyzz = pbuffer.data(idx_dip_hh + 1277);

    auto tr_z_yyzzz_yyzzz = pbuffer.data(idx_dip_hh + 1278);

    auto tr_z_yyzzz_yzzzz = pbuffer.data(idx_dip_hh + 1279);

    auto tr_z_yyzzz_zzzzz = pbuffer.data(idx_dip_hh + 1280);

    auto tr_z_yzzzz_xxxxx = pbuffer.data(idx_dip_hh + 1281);

    auto tr_z_yzzzz_xxxxy = pbuffer.data(idx_dip_hh + 1282);

    auto tr_z_yzzzz_xxxxz = pbuffer.data(idx_dip_hh + 1283);

    auto tr_z_yzzzz_xxxyy = pbuffer.data(idx_dip_hh + 1284);

    auto tr_z_yzzzz_xxxyz = pbuffer.data(idx_dip_hh + 1285);

    auto tr_z_yzzzz_xxxzz = pbuffer.data(idx_dip_hh + 1286);

    auto tr_z_yzzzz_xxyyy = pbuffer.data(idx_dip_hh + 1287);

    auto tr_z_yzzzz_xxyyz = pbuffer.data(idx_dip_hh + 1288);

    auto tr_z_yzzzz_xxyzz = pbuffer.data(idx_dip_hh + 1289);

    auto tr_z_yzzzz_xxzzz = pbuffer.data(idx_dip_hh + 1290);

    auto tr_z_yzzzz_xyyyy = pbuffer.data(idx_dip_hh + 1291);

    auto tr_z_yzzzz_xyyyz = pbuffer.data(idx_dip_hh + 1292);

    auto tr_z_yzzzz_xyyzz = pbuffer.data(idx_dip_hh + 1293);

    auto tr_z_yzzzz_xyzzz = pbuffer.data(idx_dip_hh + 1294);

    auto tr_z_yzzzz_xzzzz = pbuffer.data(idx_dip_hh + 1295);

    auto tr_z_yzzzz_yyyyy = pbuffer.data(idx_dip_hh + 1296);

    auto tr_z_yzzzz_yyyyz = pbuffer.data(idx_dip_hh + 1297);

    auto tr_z_yzzzz_yyyzz = pbuffer.data(idx_dip_hh + 1298);

    auto tr_z_yzzzz_yyzzz = pbuffer.data(idx_dip_hh + 1299);

    auto tr_z_yzzzz_yzzzz = pbuffer.data(idx_dip_hh + 1300);

    auto tr_z_yzzzz_zzzzz = pbuffer.data(idx_dip_hh + 1301);

    auto tr_z_zzzzz_xxxxx = pbuffer.data(idx_dip_hh + 1302);

    auto tr_z_zzzzz_xxxxy = pbuffer.data(idx_dip_hh + 1303);

    auto tr_z_zzzzz_xxxxz = pbuffer.data(idx_dip_hh + 1304);

    auto tr_z_zzzzz_xxxyy = pbuffer.data(idx_dip_hh + 1305);

    auto tr_z_zzzzz_xxxyz = pbuffer.data(idx_dip_hh + 1306);

    auto tr_z_zzzzz_xxxzz = pbuffer.data(idx_dip_hh + 1307);

    auto tr_z_zzzzz_xxyyy = pbuffer.data(idx_dip_hh + 1308);

    auto tr_z_zzzzz_xxyyz = pbuffer.data(idx_dip_hh + 1309);

    auto tr_z_zzzzz_xxyzz = pbuffer.data(idx_dip_hh + 1310);

    auto tr_z_zzzzz_xxzzz = pbuffer.data(idx_dip_hh + 1311);

    auto tr_z_zzzzz_xyyyy = pbuffer.data(idx_dip_hh + 1312);

    auto tr_z_zzzzz_xyyyz = pbuffer.data(idx_dip_hh + 1313);

    auto tr_z_zzzzz_xyyzz = pbuffer.data(idx_dip_hh + 1314);

    auto tr_z_zzzzz_xyzzz = pbuffer.data(idx_dip_hh + 1315);

    auto tr_z_zzzzz_xzzzz = pbuffer.data(idx_dip_hh + 1316);

    auto tr_z_zzzzz_yyyyy = pbuffer.data(idx_dip_hh + 1317);

    auto tr_z_zzzzz_yyyyz = pbuffer.data(idx_dip_hh + 1318);

    auto tr_z_zzzzz_yyyzz = pbuffer.data(idx_dip_hh + 1319);

    auto tr_z_zzzzz_yyzzz = pbuffer.data(idx_dip_hh + 1320);

    auto tr_z_zzzzz_yzzzz = pbuffer.data(idx_dip_hh + 1321);

    auto tr_z_zzzzz_zzzzz = pbuffer.data(idx_dip_hh + 1322);

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

    auto ts_xxxxz_xxxxxz = pbuffer.data(idx_ovl_hi + 58);

    auto ts_xxxxz_xxxxzz = pbuffer.data(idx_ovl_hi + 61);

    auto ts_xxxxz_xxxzzz = pbuffer.data(idx_ovl_hi + 65);

    auto ts_xxxxz_xxzzzz = pbuffer.data(idx_ovl_hi + 70);

    auto ts_xxxxz_xzzzzz = pbuffer.data(idx_ovl_hi + 76);

    auto ts_xxxyy_xxxxxy = pbuffer.data(idx_ovl_hi + 85);

    auto ts_xxxyy_xxxxyy = pbuffer.data(idx_ovl_hi + 87);

    auto ts_xxxyy_xxxyyy = pbuffer.data(idx_ovl_hi + 90);

    auto ts_xxxyy_xxyyyy = pbuffer.data(idx_ovl_hi + 94);

    auto ts_xxxyy_xyyyyy = pbuffer.data(idx_ovl_hi + 99);

    auto ts_xxxyy_yyyyyy = pbuffer.data(idx_ovl_hi + 105);

    auto ts_xxxyy_yyyyyz = pbuffer.data(idx_ovl_hi + 106);

    auto ts_xxxyy_yyyyzz = pbuffer.data(idx_ovl_hi + 107);

    auto ts_xxxyy_yyyzzz = pbuffer.data(idx_ovl_hi + 108);

    auto ts_xxxyy_yyzzzz = pbuffer.data(idx_ovl_hi + 109);

    auto ts_xxxyy_yzzzzz = pbuffer.data(idx_ovl_hi + 110);

    auto ts_xxxzz_xxxxxx = pbuffer.data(idx_ovl_hi + 140);

    auto ts_xxxzz_xxxxxz = pbuffer.data(idx_ovl_hi + 142);

    auto ts_xxxzz_xxxxzz = pbuffer.data(idx_ovl_hi + 145);

    auto ts_xxxzz_xxxzzz = pbuffer.data(idx_ovl_hi + 149);

    auto ts_xxxzz_xxzzzz = pbuffer.data(idx_ovl_hi + 154);

    auto ts_xxxzz_xzzzzz = pbuffer.data(idx_ovl_hi + 160);

    auto ts_xxxzz_yyyyyz = pbuffer.data(idx_ovl_hi + 162);

    auto ts_xxxzz_yyyyzz = pbuffer.data(idx_ovl_hi + 163);

    auto ts_xxxzz_yyyzzz = pbuffer.data(idx_ovl_hi + 164);

    auto ts_xxxzz_yyzzzz = pbuffer.data(idx_ovl_hi + 165);

    auto ts_xxxzz_yzzzzz = pbuffer.data(idx_ovl_hi + 166);

    auto ts_xxxzz_zzzzzz = pbuffer.data(idx_ovl_hi + 167);

    auto ts_xxyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 169);

    auto ts_xxyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 171);

    auto ts_xxyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 174);

    auto ts_xxyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 178);

    auto ts_xxyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 183);

    auto ts_xxyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 189);

    auto ts_xxyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 190);

    auto ts_xxyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 191);

    auto ts_xxyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 192);

    auto ts_xxyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 193);

    auto ts_xxyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 194);

    auto ts_xxzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 252);

    auto ts_xxzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 254);

    auto ts_xxzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 257);

    auto ts_xxzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 261);

    auto ts_xxzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 266);

    auto ts_xxzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 272);

    auto ts_xxzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 274);

    auto ts_xxzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 275);

    auto ts_xxzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 276);

    auto ts_xxzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 277);

    auto ts_xxzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 278);

    auto ts_xxzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 279);

    auto ts_xyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 301);

    auto ts_xyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 302);

    auto ts_xyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 303);

    auto ts_xyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 304);

    auto ts_xyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 305);

    auto ts_xyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 306);

    auto ts_xyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 358);

    auto ts_xyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 359);

    auto ts_xyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 360);

    auto ts_xyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 361);

    auto ts_xyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 362);

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

    auto ts_yyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 470);

    auto ts_yyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 471);

    auto ts_yyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 472);

    auto ts_yyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 473);

    auto ts_yyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 474);

    auto ts_yyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 475);

    auto ts_yyyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 478);

    auto ts_yyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 480);

    auto ts_yyyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 481);

    auto ts_yyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 483);

    auto ts_yyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 484);

    auto ts_yyyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 485);

    auto ts_yyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 487);

    auto ts_yyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 488);

    auto ts_yyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 489);

    auto ts_yyyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 490);

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

    auto ts_yyzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 506);

    auto ts_yyzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 508);

    auto ts_yyzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 509);

    auto ts_yyzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 511);

    auto ts_yyzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 512);

    auto ts_yyzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 513);

    auto ts_yyzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 515);

    auto ts_yyzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 516);

    auto ts_yyzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 517);

    auto ts_yyzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 518);

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

    auto ts_yzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 534);

    auto ts_yzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 537);

    auto ts_yzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 541);

    auto ts_yzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 546);

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

    // Set up components of auxiliary buffer : HI

    auto tr_x_xxxxx_xxxxxx = pbuffer.data(idx_dip_hi);

    auto tr_x_xxxxx_xxxxxy = pbuffer.data(idx_dip_hi + 1);

    auto tr_x_xxxxx_xxxxxz = pbuffer.data(idx_dip_hi + 2);

    auto tr_x_xxxxx_xxxxyy = pbuffer.data(idx_dip_hi + 3);

    auto tr_x_xxxxx_xxxxyz = pbuffer.data(idx_dip_hi + 4);

    auto tr_x_xxxxx_xxxxzz = pbuffer.data(idx_dip_hi + 5);

    auto tr_x_xxxxx_xxxyyy = pbuffer.data(idx_dip_hi + 6);

    auto tr_x_xxxxx_xxxyyz = pbuffer.data(idx_dip_hi + 7);

    auto tr_x_xxxxx_xxxyzz = pbuffer.data(idx_dip_hi + 8);

    auto tr_x_xxxxx_xxxzzz = pbuffer.data(idx_dip_hi + 9);

    auto tr_x_xxxxx_xxyyyy = pbuffer.data(idx_dip_hi + 10);

    auto tr_x_xxxxx_xxyyyz = pbuffer.data(idx_dip_hi + 11);

    auto tr_x_xxxxx_xxyyzz = pbuffer.data(idx_dip_hi + 12);

    auto tr_x_xxxxx_xxyzzz = pbuffer.data(idx_dip_hi + 13);

    auto tr_x_xxxxx_xxzzzz = pbuffer.data(idx_dip_hi + 14);

    auto tr_x_xxxxx_xyyyyy = pbuffer.data(idx_dip_hi + 15);

    auto tr_x_xxxxx_xyyyyz = pbuffer.data(idx_dip_hi + 16);

    auto tr_x_xxxxx_xyyyzz = pbuffer.data(idx_dip_hi + 17);

    auto tr_x_xxxxx_xyyzzz = pbuffer.data(idx_dip_hi + 18);

    auto tr_x_xxxxx_xyzzzz = pbuffer.data(idx_dip_hi + 19);

    auto tr_x_xxxxx_xzzzzz = pbuffer.data(idx_dip_hi + 20);

    auto tr_x_xxxxx_yyyyyy = pbuffer.data(idx_dip_hi + 21);

    auto tr_x_xxxxx_yyyyyz = pbuffer.data(idx_dip_hi + 22);

    auto tr_x_xxxxx_yyyyzz = pbuffer.data(idx_dip_hi + 23);

    auto tr_x_xxxxx_yyyzzz = pbuffer.data(idx_dip_hi + 24);

    auto tr_x_xxxxx_yyzzzz = pbuffer.data(idx_dip_hi + 25);

    auto tr_x_xxxxx_yzzzzz = pbuffer.data(idx_dip_hi + 26);

    auto tr_x_xxxxx_zzzzzz = pbuffer.data(idx_dip_hi + 27);

    auto tr_x_xxxxy_xxxxxx = pbuffer.data(idx_dip_hi + 28);

    auto tr_x_xxxxy_xxxxxy = pbuffer.data(idx_dip_hi + 29);

    auto tr_x_xxxxy_xxxxxz = pbuffer.data(idx_dip_hi + 30);

    auto tr_x_xxxxy_xxxxyy = pbuffer.data(idx_dip_hi + 31);

    auto tr_x_xxxxy_xxxxyz = pbuffer.data(idx_dip_hi + 32);

    auto tr_x_xxxxy_xxxxzz = pbuffer.data(idx_dip_hi + 33);

    auto tr_x_xxxxy_xxxyyy = pbuffer.data(idx_dip_hi + 34);

    auto tr_x_xxxxy_xxxyyz = pbuffer.data(idx_dip_hi + 35);

    auto tr_x_xxxxy_xxxyzz = pbuffer.data(idx_dip_hi + 36);

    auto tr_x_xxxxy_xxxzzz = pbuffer.data(idx_dip_hi + 37);

    auto tr_x_xxxxy_xxyyyy = pbuffer.data(idx_dip_hi + 38);

    auto tr_x_xxxxy_xxyyyz = pbuffer.data(idx_dip_hi + 39);

    auto tr_x_xxxxy_xxyyzz = pbuffer.data(idx_dip_hi + 40);

    auto tr_x_xxxxy_xxyzzz = pbuffer.data(idx_dip_hi + 41);

    auto tr_x_xxxxy_xxzzzz = pbuffer.data(idx_dip_hi + 42);

    auto tr_x_xxxxy_xyyyyy = pbuffer.data(idx_dip_hi + 43);

    auto tr_x_xxxxy_xyyyyz = pbuffer.data(idx_dip_hi + 44);

    auto tr_x_xxxxy_xyyyzz = pbuffer.data(idx_dip_hi + 45);

    auto tr_x_xxxxy_xyyzzz = pbuffer.data(idx_dip_hi + 46);

    auto tr_x_xxxxy_xyzzzz = pbuffer.data(idx_dip_hi + 47);

    auto tr_x_xxxxy_xzzzzz = pbuffer.data(idx_dip_hi + 48);

    auto tr_x_xxxxy_yyyyyy = pbuffer.data(idx_dip_hi + 49);

    auto tr_x_xxxxy_zzzzzz = pbuffer.data(idx_dip_hi + 55);

    auto tr_x_xxxxz_xxxxxx = pbuffer.data(idx_dip_hi + 56);

    auto tr_x_xxxxz_xxxxxy = pbuffer.data(idx_dip_hi + 57);

    auto tr_x_xxxxz_xxxxxz = pbuffer.data(idx_dip_hi + 58);

    auto tr_x_xxxxz_xxxxyy = pbuffer.data(idx_dip_hi + 59);

    auto tr_x_xxxxz_xxxxyz = pbuffer.data(idx_dip_hi + 60);

    auto tr_x_xxxxz_xxxxzz = pbuffer.data(idx_dip_hi + 61);

    auto tr_x_xxxxz_xxxyyy = pbuffer.data(idx_dip_hi + 62);

    auto tr_x_xxxxz_xxxyyz = pbuffer.data(idx_dip_hi + 63);

    auto tr_x_xxxxz_xxxyzz = pbuffer.data(idx_dip_hi + 64);

    auto tr_x_xxxxz_xxxzzz = pbuffer.data(idx_dip_hi + 65);

    auto tr_x_xxxxz_xxyyyy = pbuffer.data(idx_dip_hi + 66);

    auto tr_x_xxxxz_xxyyyz = pbuffer.data(idx_dip_hi + 67);

    auto tr_x_xxxxz_xxyyzz = pbuffer.data(idx_dip_hi + 68);

    auto tr_x_xxxxz_xxyzzz = pbuffer.data(idx_dip_hi + 69);

    auto tr_x_xxxxz_xxzzzz = pbuffer.data(idx_dip_hi + 70);

    auto tr_x_xxxxz_xyyyyy = pbuffer.data(idx_dip_hi + 71);

    auto tr_x_xxxxz_xyyyyz = pbuffer.data(idx_dip_hi + 72);

    auto tr_x_xxxxz_xyyyzz = pbuffer.data(idx_dip_hi + 73);

    auto tr_x_xxxxz_xyyzzz = pbuffer.data(idx_dip_hi + 74);

    auto tr_x_xxxxz_xyzzzz = pbuffer.data(idx_dip_hi + 75);

    auto tr_x_xxxxz_xzzzzz = pbuffer.data(idx_dip_hi + 76);

    auto tr_x_xxxxz_yyyyyy = pbuffer.data(idx_dip_hi + 77);

    auto tr_x_xxxxz_yyyyyz = pbuffer.data(idx_dip_hi + 78);

    auto tr_x_xxxxz_yyyyzz = pbuffer.data(idx_dip_hi + 79);

    auto tr_x_xxxxz_yyyzzz = pbuffer.data(idx_dip_hi + 80);

    auto tr_x_xxxxz_yyzzzz = pbuffer.data(idx_dip_hi + 81);

    auto tr_x_xxxxz_yzzzzz = pbuffer.data(idx_dip_hi + 82);

    auto tr_x_xxxxz_zzzzzz = pbuffer.data(idx_dip_hi + 83);

    auto tr_x_xxxyy_xxxxxx = pbuffer.data(idx_dip_hi + 84);

    auto tr_x_xxxyy_xxxxxy = pbuffer.data(idx_dip_hi + 85);

    auto tr_x_xxxyy_xxxxxz = pbuffer.data(idx_dip_hi + 86);

    auto tr_x_xxxyy_xxxxyy = pbuffer.data(idx_dip_hi + 87);

    auto tr_x_xxxyy_xxxxyz = pbuffer.data(idx_dip_hi + 88);

    auto tr_x_xxxyy_xxxxzz = pbuffer.data(idx_dip_hi + 89);

    auto tr_x_xxxyy_xxxyyy = pbuffer.data(idx_dip_hi + 90);

    auto tr_x_xxxyy_xxxyyz = pbuffer.data(idx_dip_hi + 91);

    auto tr_x_xxxyy_xxxyzz = pbuffer.data(idx_dip_hi + 92);

    auto tr_x_xxxyy_xxxzzz = pbuffer.data(idx_dip_hi + 93);

    auto tr_x_xxxyy_xxyyyy = pbuffer.data(idx_dip_hi + 94);

    auto tr_x_xxxyy_xxyyyz = pbuffer.data(idx_dip_hi + 95);

    auto tr_x_xxxyy_xxyyzz = pbuffer.data(idx_dip_hi + 96);

    auto tr_x_xxxyy_xxyzzz = pbuffer.data(idx_dip_hi + 97);

    auto tr_x_xxxyy_xxzzzz = pbuffer.data(idx_dip_hi + 98);

    auto tr_x_xxxyy_xyyyyy = pbuffer.data(idx_dip_hi + 99);

    auto tr_x_xxxyy_xyyyyz = pbuffer.data(idx_dip_hi + 100);

    auto tr_x_xxxyy_xyyyzz = pbuffer.data(idx_dip_hi + 101);

    auto tr_x_xxxyy_xyyzzz = pbuffer.data(idx_dip_hi + 102);

    auto tr_x_xxxyy_xyzzzz = pbuffer.data(idx_dip_hi + 103);

    auto tr_x_xxxyy_xzzzzz = pbuffer.data(idx_dip_hi + 104);

    auto tr_x_xxxyy_yyyyyy = pbuffer.data(idx_dip_hi + 105);

    auto tr_x_xxxyy_yyyyyz = pbuffer.data(idx_dip_hi + 106);

    auto tr_x_xxxyy_yyyyzz = pbuffer.data(idx_dip_hi + 107);

    auto tr_x_xxxyy_yyyzzz = pbuffer.data(idx_dip_hi + 108);

    auto tr_x_xxxyy_yyzzzz = pbuffer.data(idx_dip_hi + 109);

    auto tr_x_xxxyy_yzzzzz = pbuffer.data(idx_dip_hi + 110);

    auto tr_x_xxxyy_zzzzzz = pbuffer.data(idx_dip_hi + 111);

    auto tr_x_xxxyz_xxxxxz = pbuffer.data(idx_dip_hi + 114);

    auto tr_x_xxxyz_xxxxzz = pbuffer.data(idx_dip_hi + 117);

    auto tr_x_xxxyz_xxxzzz = pbuffer.data(idx_dip_hi + 121);

    auto tr_x_xxxyz_xxzzzz = pbuffer.data(idx_dip_hi + 126);

    auto tr_x_xxxyz_xzzzzz = pbuffer.data(idx_dip_hi + 132);

    auto tr_x_xxxyz_zzzzzz = pbuffer.data(idx_dip_hi + 139);

    auto tr_x_xxxzz_xxxxxx = pbuffer.data(idx_dip_hi + 140);

    auto tr_x_xxxzz_xxxxxy = pbuffer.data(idx_dip_hi + 141);

    auto tr_x_xxxzz_xxxxxz = pbuffer.data(idx_dip_hi + 142);

    auto tr_x_xxxzz_xxxxyy = pbuffer.data(idx_dip_hi + 143);

    auto tr_x_xxxzz_xxxxyz = pbuffer.data(idx_dip_hi + 144);

    auto tr_x_xxxzz_xxxxzz = pbuffer.data(idx_dip_hi + 145);

    auto tr_x_xxxzz_xxxyyy = pbuffer.data(idx_dip_hi + 146);

    auto tr_x_xxxzz_xxxyyz = pbuffer.data(idx_dip_hi + 147);

    auto tr_x_xxxzz_xxxyzz = pbuffer.data(idx_dip_hi + 148);

    auto tr_x_xxxzz_xxxzzz = pbuffer.data(idx_dip_hi + 149);

    auto tr_x_xxxzz_xxyyyy = pbuffer.data(idx_dip_hi + 150);

    auto tr_x_xxxzz_xxyyyz = pbuffer.data(idx_dip_hi + 151);

    auto tr_x_xxxzz_xxyyzz = pbuffer.data(idx_dip_hi + 152);

    auto tr_x_xxxzz_xxyzzz = pbuffer.data(idx_dip_hi + 153);

    auto tr_x_xxxzz_xxzzzz = pbuffer.data(idx_dip_hi + 154);

    auto tr_x_xxxzz_xyyyyy = pbuffer.data(idx_dip_hi + 155);

    auto tr_x_xxxzz_xyyyyz = pbuffer.data(idx_dip_hi + 156);

    auto tr_x_xxxzz_xyyyzz = pbuffer.data(idx_dip_hi + 157);

    auto tr_x_xxxzz_xyyzzz = pbuffer.data(idx_dip_hi + 158);

    auto tr_x_xxxzz_xyzzzz = pbuffer.data(idx_dip_hi + 159);

    auto tr_x_xxxzz_xzzzzz = pbuffer.data(idx_dip_hi + 160);

    auto tr_x_xxxzz_yyyyyy = pbuffer.data(idx_dip_hi + 161);

    auto tr_x_xxxzz_yyyyyz = pbuffer.data(idx_dip_hi + 162);

    auto tr_x_xxxzz_yyyyzz = pbuffer.data(idx_dip_hi + 163);

    auto tr_x_xxxzz_yyyzzz = pbuffer.data(idx_dip_hi + 164);

    auto tr_x_xxxzz_yyzzzz = pbuffer.data(idx_dip_hi + 165);

    auto tr_x_xxxzz_yzzzzz = pbuffer.data(idx_dip_hi + 166);

    auto tr_x_xxxzz_zzzzzz = pbuffer.data(idx_dip_hi + 167);

    auto tr_x_xxyyy_xxxxxx = pbuffer.data(idx_dip_hi + 168);

    auto tr_x_xxyyy_xxxxxy = pbuffer.data(idx_dip_hi + 169);

    auto tr_x_xxyyy_xxxxxz = pbuffer.data(idx_dip_hi + 170);

    auto tr_x_xxyyy_xxxxyy = pbuffer.data(idx_dip_hi + 171);

    auto tr_x_xxyyy_xxxxyz = pbuffer.data(idx_dip_hi + 172);

    auto tr_x_xxyyy_xxxxzz = pbuffer.data(idx_dip_hi + 173);

    auto tr_x_xxyyy_xxxyyy = pbuffer.data(idx_dip_hi + 174);

    auto tr_x_xxyyy_xxxyyz = pbuffer.data(idx_dip_hi + 175);

    auto tr_x_xxyyy_xxxyzz = pbuffer.data(idx_dip_hi + 176);

    auto tr_x_xxyyy_xxxzzz = pbuffer.data(idx_dip_hi + 177);

    auto tr_x_xxyyy_xxyyyy = pbuffer.data(idx_dip_hi + 178);

    auto tr_x_xxyyy_xxyyyz = pbuffer.data(idx_dip_hi + 179);

    auto tr_x_xxyyy_xxyyzz = pbuffer.data(idx_dip_hi + 180);

    auto tr_x_xxyyy_xxyzzz = pbuffer.data(idx_dip_hi + 181);

    auto tr_x_xxyyy_xxzzzz = pbuffer.data(idx_dip_hi + 182);

    auto tr_x_xxyyy_xyyyyy = pbuffer.data(idx_dip_hi + 183);

    auto tr_x_xxyyy_xyyyyz = pbuffer.data(idx_dip_hi + 184);

    auto tr_x_xxyyy_xyyyzz = pbuffer.data(idx_dip_hi + 185);

    auto tr_x_xxyyy_xyyzzz = pbuffer.data(idx_dip_hi + 186);

    auto tr_x_xxyyy_xyzzzz = pbuffer.data(idx_dip_hi + 187);

    auto tr_x_xxyyy_xzzzzz = pbuffer.data(idx_dip_hi + 188);

    auto tr_x_xxyyy_yyyyyy = pbuffer.data(idx_dip_hi + 189);

    auto tr_x_xxyyy_yyyyyz = pbuffer.data(idx_dip_hi + 190);

    auto tr_x_xxyyy_yyyyzz = pbuffer.data(idx_dip_hi + 191);

    auto tr_x_xxyyy_yyyzzz = pbuffer.data(idx_dip_hi + 192);

    auto tr_x_xxyyy_yyzzzz = pbuffer.data(idx_dip_hi + 193);

    auto tr_x_xxyyy_yzzzzz = pbuffer.data(idx_dip_hi + 194);

    auto tr_x_xxyyy_zzzzzz = pbuffer.data(idx_dip_hi + 195);

    auto tr_x_xxyyz_xxxxxy = pbuffer.data(idx_dip_hi + 197);

    auto tr_x_xxyyz_xxxxxz = pbuffer.data(idx_dip_hi + 198);

    auto tr_x_xxyyz_xxxxyy = pbuffer.data(idx_dip_hi + 199);

    auto tr_x_xxyyz_xxxxzz = pbuffer.data(idx_dip_hi + 201);

    auto tr_x_xxyyz_xxxyyy = pbuffer.data(idx_dip_hi + 202);

    auto tr_x_xxyyz_xxxzzz = pbuffer.data(idx_dip_hi + 205);

    auto tr_x_xxyyz_xxyyyy = pbuffer.data(idx_dip_hi + 206);

    auto tr_x_xxyyz_xxzzzz = pbuffer.data(idx_dip_hi + 210);

    auto tr_x_xxyyz_xyyyyy = pbuffer.data(idx_dip_hi + 211);

    auto tr_x_xxyyz_xzzzzz = pbuffer.data(idx_dip_hi + 216);

    auto tr_x_xxyyz_yyyyyy = pbuffer.data(idx_dip_hi + 217);

    auto tr_x_xxyyz_zzzzzz = pbuffer.data(idx_dip_hi + 223);

    auto tr_x_xxyzz_xxxxxx = pbuffer.data(idx_dip_hi + 224);

    auto tr_x_xxyzz_xxxxxz = pbuffer.data(idx_dip_hi + 226);

    auto tr_x_xxyzz_xxxxyz = pbuffer.data(idx_dip_hi + 228);

    auto tr_x_xxyzz_xxxxzz = pbuffer.data(idx_dip_hi + 229);

    auto tr_x_xxyzz_xxxyyz = pbuffer.data(idx_dip_hi + 231);

    auto tr_x_xxyzz_xxxyzz = pbuffer.data(idx_dip_hi + 232);

    auto tr_x_xxyzz_xxxzzz = pbuffer.data(idx_dip_hi + 233);

    auto tr_x_xxyzz_xxyyyz = pbuffer.data(idx_dip_hi + 235);

    auto tr_x_xxyzz_xxyyzz = pbuffer.data(idx_dip_hi + 236);

    auto tr_x_xxyzz_xxyzzz = pbuffer.data(idx_dip_hi + 237);

    auto tr_x_xxyzz_xxzzzz = pbuffer.data(idx_dip_hi + 238);

    auto tr_x_xxyzz_xyyyyz = pbuffer.data(idx_dip_hi + 240);

    auto tr_x_xxyzz_xyyyzz = pbuffer.data(idx_dip_hi + 241);

    auto tr_x_xxyzz_xyyzzz = pbuffer.data(idx_dip_hi + 242);

    auto tr_x_xxyzz_xyzzzz = pbuffer.data(idx_dip_hi + 243);

    auto tr_x_xxyzz_xzzzzz = pbuffer.data(idx_dip_hi + 244);

    auto tr_x_xxyzz_zzzzzz = pbuffer.data(idx_dip_hi + 251);

    auto tr_x_xxzzz_xxxxxx = pbuffer.data(idx_dip_hi + 252);

    auto tr_x_xxzzz_xxxxxy = pbuffer.data(idx_dip_hi + 253);

    auto tr_x_xxzzz_xxxxxz = pbuffer.data(idx_dip_hi + 254);

    auto tr_x_xxzzz_xxxxyy = pbuffer.data(idx_dip_hi + 255);

    auto tr_x_xxzzz_xxxxyz = pbuffer.data(idx_dip_hi + 256);

    auto tr_x_xxzzz_xxxxzz = pbuffer.data(idx_dip_hi + 257);

    auto tr_x_xxzzz_xxxyyy = pbuffer.data(idx_dip_hi + 258);

    auto tr_x_xxzzz_xxxyyz = pbuffer.data(idx_dip_hi + 259);

    auto tr_x_xxzzz_xxxyzz = pbuffer.data(idx_dip_hi + 260);

    auto tr_x_xxzzz_xxxzzz = pbuffer.data(idx_dip_hi + 261);

    auto tr_x_xxzzz_xxyyyy = pbuffer.data(idx_dip_hi + 262);

    auto tr_x_xxzzz_xxyyyz = pbuffer.data(idx_dip_hi + 263);

    auto tr_x_xxzzz_xxyyzz = pbuffer.data(idx_dip_hi + 264);

    auto tr_x_xxzzz_xxyzzz = pbuffer.data(idx_dip_hi + 265);

    auto tr_x_xxzzz_xxzzzz = pbuffer.data(idx_dip_hi + 266);

    auto tr_x_xxzzz_xyyyyy = pbuffer.data(idx_dip_hi + 267);

    auto tr_x_xxzzz_xyyyyz = pbuffer.data(idx_dip_hi + 268);

    auto tr_x_xxzzz_xyyyzz = pbuffer.data(idx_dip_hi + 269);

    auto tr_x_xxzzz_xyyzzz = pbuffer.data(idx_dip_hi + 270);

    auto tr_x_xxzzz_xyzzzz = pbuffer.data(idx_dip_hi + 271);

    auto tr_x_xxzzz_xzzzzz = pbuffer.data(idx_dip_hi + 272);

    auto tr_x_xxzzz_yyyyyy = pbuffer.data(idx_dip_hi + 273);

    auto tr_x_xxzzz_yyyyyz = pbuffer.data(idx_dip_hi + 274);

    auto tr_x_xxzzz_yyyyzz = pbuffer.data(idx_dip_hi + 275);

    auto tr_x_xxzzz_yyyzzz = pbuffer.data(idx_dip_hi + 276);

    auto tr_x_xxzzz_yyzzzz = pbuffer.data(idx_dip_hi + 277);

    auto tr_x_xxzzz_yzzzzz = pbuffer.data(idx_dip_hi + 278);

    auto tr_x_xxzzz_zzzzzz = pbuffer.data(idx_dip_hi + 279);

    auto tr_x_xyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 280);

    auto tr_x_xyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 281);

    auto tr_x_xyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 282);

    auto tr_x_xyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 283);

    auto tr_x_xyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 284);

    auto tr_x_xyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 285);

    auto tr_x_xyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 286);

    auto tr_x_xyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 287);

    auto tr_x_xyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 288);

    auto tr_x_xyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 289);

    auto tr_x_xyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 290);

    auto tr_x_xyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 291);

    auto tr_x_xyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 292);

    auto tr_x_xyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 293);

    auto tr_x_xyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 294);

    auto tr_x_xyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 295);

    auto tr_x_xyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 296);

    auto tr_x_xyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 297);

    auto tr_x_xyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 298);

    auto tr_x_xyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 299);

    auto tr_x_xyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 300);

    auto tr_x_xyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 301);

    auto tr_x_xyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 302);

    auto tr_x_xyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 303);

    auto tr_x_xyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 304);

    auto tr_x_xyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 305);

    auto tr_x_xyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 306);

    auto tr_x_xyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 309);

    auto tr_x_xyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 310);

    auto tr_x_xyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 311);

    auto tr_x_xyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 313);

    auto tr_x_xyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 314);

    auto tr_x_xyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 317);

    auto tr_x_xyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 318);

    auto tr_x_xyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 322);

    auto tr_x_xyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 323);

    auto tr_x_xyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 328);

    auto tr_x_xyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 336);

    auto tr_x_xyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 337);

    auto tr_x_xyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 338);

    auto tr_x_xyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 339);

    auto tr_x_xyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 341);

    auto tr_x_xyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 342);

    auto tr_x_xyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 345);

    auto tr_x_xyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 346);

    auto tr_x_xyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 350);

    auto tr_x_xyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 351);

    auto tr_x_xyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 356);

    auto tr_x_xyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 358);

    auto tr_x_xyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 359);

    auto tr_x_xyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 360);

    auto tr_x_xyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 361);

    auto tr_x_xyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 362);

    auto tr_x_xyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 364);

    auto tr_x_xyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 366);

    auto tr_x_xyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 369);

    auto tr_x_xyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 373);

    auto tr_x_xyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 378);

    auto tr_x_xyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 384);

    auto tr_x_xzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 392);

    auto tr_x_xzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 393);

    auto tr_x_xzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 394);

    auto tr_x_xzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 395);

    auto tr_x_xzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 396);

    auto tr_x_xzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 397);

    auto tr_x_xzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 398);

    auto tr_x_xzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 399);

    auto tr_x_xzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 400);

    auto tr_x_xzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 401);

    auto tr_x_xzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 402);

    auto tr_x_xzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 403);

    auto tr_x_xzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 404);

    auto tr_x_xzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 405);

    auto tr_x_xzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 406);

    auto tr_x_xzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 407);

    auto tr_x_xzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 408);

    auto tr_x_xzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 409);

    auto tr_x_xzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 410);

    auto tr_x_xzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 411);

    auto tr_x_xzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 412);

    auto tr_x_xzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 414);

    auto tr_x_xzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 415);

    auto tr_x_xzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 416);

    auto tr_x_xzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 417);

    auto tr_x_xzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 418);

    auto tr_x_xzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 419);

    auto tr_x_yyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 420);

    auto tr_x_yyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 421);

    auto tr_x_yyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 422);

    auto tr_x_yyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 423);

    auto tr_x_yyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 424);

    auto tr_x_yyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 425);

    auto tr_x_yyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 426);

    auto tr_x_yyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 427);

    auto tr_x_yyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 428);

    auto tr_x_yyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 429);

    auto tr_x_yyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 430);

    auto tr_x_yyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 431);

    auto tr_x_yyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 432);

    auto tr_x_yyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 433);

    auto tr_x_yyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 434);

    auto tr_x_yyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 435);

    auto tr_x_yyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 436);

    auto tr_x_yyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 437);

    auto tr_x_yyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 438);

    auto tr_x_yyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 439);

    auto tr_x_yyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 440);

    auto tr_x_yyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 441);

    auto tr_x_yyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 442);

    auto tr_x_yyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 443);

    auto tr_x_yyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 444);

    auto tr_x_yyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 445);

    auto tr_x_yyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 446);

    auto tr_x_yyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 447);

    auto tr_x_yyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 449);

    auto tr_x_yyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 450);

    auto tr_x_yyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 451);

    auto tr_x_yyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 453);

    auto tr_x_yyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 454);

    auto tr_x_yyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 457);

    auto tr_x_yyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 458);

    auto tr_x_yyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 462);

    auto tr_x_yyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 463);

    auto tr_x_yyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 468);

    auto tr_x_yyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 469);

    auto tr_x_yyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 470);

    auto tr_x_yyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 471);

    auto tr_x_yyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 472);

    auto tr_x_yyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 473);

    auto tr_x_yyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 474);

    auto tr_x_yyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 475);

    auto tr_x_yyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 476);

    auto tr_x_yyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 477);

    auto tr_x_yyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 478);

    auto tr_x_yyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 479);

    auto tr_x_yyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 480);

    auto tr_x_yyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 481);

    auto tr_x_yyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 482);

    auto tr_x_yyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 483);

    auto tr_x_yyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 484);

    auto tr_x_yyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 485);

    auto tr_x_yyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 486);

    auto tr_x_yyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 487);

    auto tr_x_yyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 488);

    auto tr_x_yyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 489);

    auto tr_x_yyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 490);

    auto tr_x_yyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 491);

    auto tr_x_yyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 492);

    auto tr_x_yyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 493);

    auto tr_x_yyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 494);

    auto tr_x_yyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 495);

    auto tr_x_yyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 496);

    auto tr_x_yyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 497);

    auto tr_x_yyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 498);

    auto tr_x_yyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 499);

    auto tr_x_yyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 500);

    auto tr_x_yyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 501);

    auto tr_x_yyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 502);

    auto tr_x_yyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 503);

    auto tr_x_yyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 504);

    auto tr_x_yyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 505);

    auto tr_x_yyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 506);

    auto tr_x_yyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 507);

    auto tr_x_yyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 508);

    auto tr_x_yyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 509);

    auto tr_x_yyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 510);

    auto tr_x_yyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 511);

    auto tr_x_yyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 512);

    auto tr_x_yyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 513);

    auto tr_x_yyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 514);

    auto tr_x_yyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 515);

    auto tr_x_yyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 516);

    auto tr_x_yyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 517);

    auto tr_x_yyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 518);

    auto tr_x_yyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 519);

    auto tr_x_yyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 520);

    auto tr_x_yyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 521);

    auto tr_x_yyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 522);

    auto tr_x_yyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 523);

    auto tr_x_yyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 524);

    auto tr_x_yyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 525);

    auto tr_x_yyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 526);

    auto tr_x_yyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 527);

    auto tr_x_yyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 528);

    auto tr_x_yyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 529);

    auto tr_x_yyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 530);

    auto tr_x_yyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 531);

    auto tr_x_yzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 532);

    auto tr_x_yzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 534);

    auto tr_x_yzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 536);

    auto tr_x_yzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 537);

    auto tr_x_yzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 539);

    auto tr_x_yzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 540);

    auto tr_x_yzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 541);

    auto tr_x_yzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 543);

    auto tr_x_yzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 544);

    auto tr_x_yzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 545);

    auto tr_x_yzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 546);

    auto tr_x_yzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 548);

    auto tr_x_yzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 549);

    auto tr_x_yzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 550);

    auto tr_x_yzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 551);

    auto tr_x_yzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 552);

    auto tr_x_yzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 553);

    auto tr_x_yzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 554);

    auto tr_x_yzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 555);

    auto tr_x_yzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 556);

    auto tr_x_yzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 557);

    auto tr_x_yzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 558);

    auto tr_x_yzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 559);

    auto tr_x_zzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 560);

    auto tr_x_zzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 561);

    auto tr_x_zzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 562);

    auto tr_x_zzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 563);

    auto tr_x_zzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 564);

    auto tr_x_zzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 565);

    auto tr_x_zzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 566);

    auto tr_x_zzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 567);

    auto tr_x_zzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 568);

    auto tr_x_zzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 569);

    auto tr_x_zzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 570);

    auto tr_x_zzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 571);

    auto tr_x_zzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 572);

    auto tr_x_zzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 573);

    auto tr_x_zzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 574);

    auto tr_x_zzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 575);

    auto tr_x_zzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 576);

    auto tr_x_zzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 577);

    auto tr_x_zzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 578);

    auto tr_x_zzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 579);

    auto tr_x_zzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 580);

    auto tr_x_zzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 581);

    auto tr_x_zzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 582);

    auto tr_x_zzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 583);

    auto tr_x_zzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 584);

    auto tr_x_zzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 585);

    auto tr_x_zzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 586);

    auto tr_x_zzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 587);

    auto tr_y_xxxxx_xxxxxx = pbuffer.data(idx_dip_hi + 588);

    auto tr_y_xxxxx_xxxxxy = pbuffer.data(idx_dip_hi + 589);

    auto tr_y_xxxxx_xxxxxz = pbuffer.data(idx_dip_hi + 590);

    auto tr_y_xxxxx_xxxxyy = pbuffer.data(idx_dip_hi + 591);

    auto tr_y_xxxxx_xxxxyz = pbuffer.data(idx_dip_hi + 592);

    auto tr_y_xxxxx_xxxxzz = pbuffer.data(idx_dip_hi + 593);

    auto tr_y_xxxxx_xxxyyy = pbuffer.data(idx_dip_hi + 594);

    auto tr_y_xxxxx_xxxyyz = pbuffer.data(idx_dip_hi + 595);

    auto tr_y_xxxxx_xxxyzz = pbuffer.data(idx_dip_hi + 596);

    auto tr_y_xxxxx_xxxzzz = pbuffer.data(idx_dip_hi + 597);

    auto tr_y_xxxxx_xxyyyy = pbuffer.data(idx_dip_hi + 598);

    auto tr_y_xxxxx_xxyyyz = pbuffer.data(idx_dip_hi + 599);

    auto tr_y_xxxxx_xxyyzz = pbuffer.data(idx_dip_hi + 600);

    auto tr_y_xxxxx_xxyzzz = pbuffer.data(idx_dip_hi + 601);

    auto tr_y_xxxxx_xxzzzz = pbuffer.data(idx_dip_hi + 602);

    auto tr_y_xxxxx_xyyyyy = pbuffer.data(idx_dip_hi + 603);

    auto tr_y_xxxxx_xyyyyz = pbuffer.data(idx_dip_hi + 604);

    auto tr_y_xxxxx_xyyyzz = pbuffer.data(idx_dip_hi + 605);

    auto tr_y_xxxxx_xyyzzz = pbuffer.data(idx_dip_hi + 606);

    auto tr_y_xxxxx_xyzzzz = pbuffer.data(idx_dip_hi + 607);

    auto tr_y_xxxxx_xzzzzz = pbuffer.data(idx_dip_hi + 608);

    auto tr_y_xxxxx_yyyyyy = pbuffer.data(idx_dip_hi + 609);

    auto tr_y_xxxxx_yyyyyz = pbuffer.data(idx_dip_hi + 610);

    auto tr_y_xxxxx_yyyyzz = pbuffer.data(idx_dip_hi + 611);

    auto tr_y_xxxxx_yyyzzz = pbuffer.data(idx_dip_hi + 612);

    auto tr_y_xxxxx_yyzzzz = pbuffer.data(idx_dip_hi + 613);

    auto tr_y_xxxxx_yzzzzz = pbuffer.data(idx_dip_hi + 614);

    auto tr_y_xxxxx_zzzzzz = pbuffer.data(idx_dip_hi + 615);

    auto tr_y_xxxxy_xxxxxx = pbuffer.data(idx_dip_hi + 616);

    auto tr_y_xxxxy_xxxxxy = pbuffer.data(idx_dip_hi + 617);

    auto tr_y_xxxxy_xxxxyy = pbuffer.data(idx_dip_hi + 619);

    auto tr_y_xxxxy_xxxxyz = pbuffer.data(idx_dip_hi + 620);

    auto tr_y_xxxxy_xxxyyy = pbuffer.data(idx_dip_hi + 622);

    auto tr_y_xxxxy_xxxyyz = pbuffer.data(idx_dip_hi + 623);

    auto tr_y_xxxxy_xxxyzz = pbuffer.data(idx_dip_hi + 624);

    auto tr_y_xxxxy_xxyyyy = pbuffer.data(idx_dip_hi + 626);

    auto tr_y_xxxxy_xxyyyz = pbuffer.data(idx_dip_hi + 627);

    auto tr_y_xxxxy_xxyyzz = pbuffer.data(idx_dip_hi + 628);

    auto tr_y_xxxxy_xxyzzz = pbuffer.data(idx_dip_hi + 629);

    auto tr_y_xxxxy_xyyyyy = pbuffer.data(idx_dip_hi + 631);

    auto tr_y_xxxxy_xyyyyz = pbuffer.data(idx_dip_hi + 632);

    auto tr_y_xxxxy_xyyyzz = pbuffer.data(idx_dip_hi + 633);

    auto tr_y_xxxxy_xyyzzz = pbuffer.data(idx_dip_hi + 634);

    auto tr_y_xxxxy_xyzzzz = pbuffer.data(idx_dip_hi + 635);

    auto tr_y_xxxxy_yyyyyy = pbuffer.data(idx_dip_hi + 637);

    auto tr_y_xxxxy_yyyyyz = pbuffer.data(idx_dip_hi + 638);

    auto tr_y_xxxxy_yyyyzz = pbuffer.data(idx_dip_hi + 639);

    auto tr_y_xxxxy_yyyzzz = pbuffer.data(idx_dip_hi + 640);

    auto tr_y_xxxxy_yyzzzz = pbuffer.data(idx_dip_hi + 641);

    auto tr_y_xxxxy_yzzzzz = pbuffer.data(idx_dip_hi + 642);

    auto tr_y_xxxxy_zzzzzz = pbuffer.data(idx_dip_hi + 643);

    auto tr_y_xxxxz_xxxxxx = pbuffer.data(idx_dip_hi + 644);

    auto tr_y_xxxxz_xxxxxy = pbuffer.data(idx_dip_hi + 645);

    auto tr_y_xxxxz_xxxxxz = pbuffer.data(idx_dip_hi + 646);

    auto tr_y_xxxxz_xxxxyy = pbuffer.data(idx_dip_hi + 647);

    auto tr_y_xxxxz_xxxxzz = pbuffer.data(idx_dip_hi + 649);

    auto tr_y_xxxxz_xxxyyy = pbuffer.data(idx_dip_hi + 650);

    auto tr_y_xxxxz_xxxzzz = pbuffer.data(idx_dip_hi + 653);

    auto tr_y_xxxxz_xxyyyy = pbuffer.data(idx_dip_hi + 654);

    auto tr_y_xxxxz_xxzzzz = pbuffer.data(idx_dip_hi + 658);

    auto tr_y_xxxxz_xyyyyy = pbuffer.data(idx_dip_hi + 659);

    auto tr_y_xxxxz_xzzzzz = pbuffer.data(idx_dip_hi + 664);

    auto tr_y_xxxxz_yyyyyz = pbuffer.data(idx_dip_hi + 666);

    auto tr_y_xxxxz_yyyyzz = pbuffer.data(idx_dip_hi + 667);

    auto tr_y_xxxxz_yyyzzz = pbuffer.data(idx_dip_hi + 668);

    auto tr_y_xxxxz_yyzzzz = pbuffer.data(idx_dip_hi + 669);

    auto tr_y_xxxxz_yzzzzz = pbuffer.data(idx_dip_hi + 670);

    auto tr_y_xxxxz_zzzzzz = pbuffer.data(idx_dip_hi + 671);

    auto tr_y_xxxyy_xxxxxx = pbuffer.data(idx_dip_hi + 672);

    auto tr_y_xxxyy_xxxxxy = pbuffer.data(idx_dip_hi + 673);

    auto tr_y_xxxyy_xxxxxz = pbuffer.data(idx_dip_hi + 674);

    auto tr_y_xxxyy_xxxxyy = pbuffer.data(idx_dip_hi + 675);

    auto tr_y_xxxyy_xxxxyz = pbuffer.data(idx_dip_hi + 676);

    auto tr_y_xxxyy_xxxxzz = pbuffer.data(idx_dip_hi + 677);

    auto tr_y_xxxyy_xxxyyy = pbuffer.data(idx_dip_hi + 678);

    auto tr_y_xxxyy_xxxyyz = pbuffer.data(idx_dip_hi + 679);

    auto tr_y_xxxyy_xxxyzz = pbuffer.data(idx_dip_hi + 680);

    auto tr_y_xxxyy_xxxzzz = pbuffer.data(idx_dip_hi + 681);

    auto tr_y_xxxyy_xxyyyy = pbuffer.data(idx_dip_hi + 682);

    auto tr_y_xxxyy_xxyyyz = pbuffer.data(idx_dip_hi + 683);

    auto tr_y_xxxyy_xxyyzz = pbuffer.data(idx_dip_hi + 684);

    auto tr_y_xxxyy_xxyzzz = pbuffer.data(idx_dip_hi + 685);

    auto tr_y_xxxyy_xxzzzz = pbuffer.data(idx_dip_hi + 686);

    auto tr_y_xxxyy_xyyyyy = pbuffer.data(idx_dip_hi + 687);

    auto tr_y_xxxyy_xyyyyz = pbuffer.data(idx_dip_hi + 688);

    auto tr_y_xxxyy_xyyyzz = pbuffer.data(idx_dip_hi + 689);

    auto tr_y_xxxyy_xyyzzz = pbuffer.data(idx_dip_hi + 690);

    auto tr_y_xxxyy_xyzzzz = pbuffer.data(idx_dip_hi + 691);

    auto tr_y_xxxyy_xzzzzz = pbuffer.data(idx_dip_hi + 692);

    auto tr_y_xxxyy_yyyyyy = pbuffer.data(idx_dip_hi + 693);

    auto tr_y_xxxyy_yyyyyz = pbuffer.data(idx_dip_hi + 694);

    auto tr_y_xxxyy_yyyyzz = pbuffer.data(idx_dip_hi + 695);

    auto tr_y_xxxyy_yyyzzz = pbuffer.data(idx_dip_hi + 696);

    auto tr_y_xxxyy_yyzzzz = pbuffer.data(idx_dip_hi + 697);

    auto tr_y_xxxyy_yzzzzz = pbuffer.data(idx_dip_hi + 698);

    auto tr_y_xxxyy_zzzzzz = pbuffer.data(idx_dip_hi + 699);

    auto tr_y_xxxyz_xxxxxy = pbuffer.data(idx_dip_hi + 701);

    auto tr_y_xxxyz_xxxxyy = pbuffer.data(idx_dip_hi + 703);

    auto tr_y_xxxyz_xxxyyy = pbuffer.data(idx_dip_hi + 706);

    auto tr_y_xxxyz_xxyyyy = pbuffer.data(idx_dip_hi + 710);

    auto tr_y_xxxyz_xyyyyy = pbuffer.data(idx_dip_hi + 715);

    auto tr_y_xxxyz_yyyyyz = pbuffer.data(idx_dip_hi + 722);

    auto tr_y_xxxyz_yyyyzz = pbuffer.data(idx_dip_hi + 723);

    auto tr_y_xxxyz_yyyzzz = pbuffer.data(idx_dip_hi + 724);

    auto tr_y_xxxyz_yyzzzz = pbuffer.data(idx_dip_hi + 725);

    auto tr_y_xxxyz_yzzzzz = pbuffer.data(idx_dip_hi + 726);

    auto tr_y_xxxyz_zzzzzz = pbuffer.data(idx_dip_hi + 727);

    auto tr_y_xxxzz_xxxxxx = pbuffer.data(idx_dip_hi + 728);

    auto tr_y_xxxzz_xxxxxy = pbuffer.data(idx_dip_hi + 729);

    auto tr_y_xxxzz_xxxxxz = pbuffer.data(idx_dip_hi + 730);

    auto tr_y_xxxzz_xxxxyy = pbuffer.data(idx_dip_hi + 731);

    auto tr_y_xxxzz_xxxxyz = pbuffer.data(idx_dip_hi + 732);

    auto tr_y_xxxzz_xxxxzz = pbuffer.data(idx_dip_hi + 733);

    auto tr_y_xxxzz_xxxyyy = pbuffer.data(idx_dip_hi + 734);

    auto tr_y_xxxzz_xxxyyz = pbuffer.data(idx_dip_hi + 735);

    auto tr_y_xxxzz_xxxyzz = pbuffer.data(idx_dip_hi + 736);

    auto tr_y_xxxzz_xxxzzz = pbuffer.data(idx_dip_hi + 737);

    auto tr_y_xxxzz_xxyyyy = pbuffer.data(idx_dip_hi + 738);

    auto tr_y_xxxzz_xxyyyz = pbuffer.data(idx_dip_hi + 739);

    auto tr_y_xxxzz_xxyyzz = pbuffer.data(idx_dip_hi + 740);

    auto tr_y_xxxzz_xxyzzz = pbuffer.data(idx_dip_hi + 741);

    auto tr_y_xxxzz_xxzzzz = pbuffer.data(idx_dip_hi + 742);

    auto tr_y_xxxzz_xyyyyy = pbuffer.data(idx_dip_hi + 743);

    auto tr_y_xxxzz_xyyyyz = pbuffer.data(idx_dip_hi + 744);

    auto tr_y_xxxzz_xyyyzz = pbuffer.data(idx_dip_hi + 745);

    auto tr_y_xxxzz_xyyzzz = pbuffer.data(idx_dip_hi + 746);

    auto tr_y_xxxzz_xyzzzz = pbuffer.data(idx_dip_hi + 747);

    auto tr_y_xxxzz_xzzzzz = pbuffer.data(idx_dip_hi + 748);

    auto tr_y_xxxzz_yyyyyy = pbuffer.data(idx_dip_hi + 749);

    auto tr_y_xxxzz_yyyyyz = pbuffer.data(idx_dip_hi + 750);

    auto tr_y_xxxzz_yyyyzz = pbuffer.data(idx_dip_hi + 751);

    auto tr_y_xxxzz_yyyzzz = pbuffer.data(idx_dip_hi + 752);

    auto tr_y_xxxzz_yyzzzz = pbuffer.data(idx_dip_hi + 753);

    auto tr_y_xxxzz_yzzzzz = pbuffer.data(idx_dip_hi + 754);

    auto tr_y_xxxzz_zzzzzz = pbuffer.data(idx_dip_hi + 755);

    auto tr_y_xxyyy_xxxxxx = pbuffer.data(idx_dip_hi + 756);

    auto tr_y_xxyyy_xxxxxy = pbuffer.data(idx_dip_hi + 757);

    auto tr_y_xxyyy_xxxxxz = pbuffer.data(idx_dip_hi + 758);

    auto tr_y_xxyyy_xxxxyy = pbuffer.data(idx_dip_hi + 759);

    auto tr_y_xxyyy_xxxxyz = pbuffer.data(idx_dip_hi + 760);

    auto tr_y_xxyyy_xxxxzz = pbuffer.data(idx_dip_hi + 761);

    auto tr_y_xxyyy_xxxyyy = pbuffer.data(idx_dip_hi + 762);

    auto tr_y_xxyyy_xxxyyz = pbuffer.data(idx_dip_hi + 763);

    auto tr_y_xxyyy_xxxyzz = pbuffer.data(idx_dip_hi + 764);

    auto tr_y_xxyyy_xxxzzz = pbuffer.data(idx_dip_hi + 765);

    auto tr_y_xxyyy_xxyyyy = pbuffer.data(idx_dip_hi + 766);

    auto tr_y_xxyyy_xxyyyz = pbuffer.data(idx_dip_hi + 767);

    auto tr_y_xxyyy_xxyyzz = pbuffer.data(idx_dip_hi + 768);

    auto tr_y_xxyyy_xxyzzz = pbuffer.data(idx_dip_hi + 769);

    auto tr_y_xxyyy_xxzzzz = pbuffer.data(idx_dip_hi + 770);

    auto tr_y_xxyyy_xyyyyy = pbuffer.data(idx_dip_hi + 771);

    auto tr_y_xxyyy_xyyyyz = pbuffer.data(idx_dip_hi + 772);

    auto tr_y_xxyyy_xyyyzz = pbuffer.data(idx_dip_hi + 773);

    auto tr_y_xxyyy_xyyzzz = pbuffer.data(idx_dip_hi + 774);

    auto tr_y_xxyyy_xyzzzz = pbuffer.data(idx_dip_hi + 775);

    auto tr_y_xxyyy_xzzzzz = pbuffer.data(idx_dip_hi + 776);

    auto tr_y_xxyyy_yyyyyy = pbuffer.data(idx_dip_hi + 777);

    auto tr_y_xxyyy_yyyyyz = pbuffer.data(idx_dip_hi + 778);

    auto tr_y_xxyyy_yyyyzz = pbuffer.data(idx_dip_hi + 779);

    auto tr_y_xxyyy_yyyzzz = pbuffer.data(idx_dip_hi + 780);

    auto tr_y_xxyyy_yyzzzz = pbuffer.data(idx_dip_hi + 781);

    auto tr_y_xxyyy_yzzzzz = pbuffer.data(idx_dip_hi + 782);

    auto tr_y_xxyyy_zzzzzz = pbuffer.data(idx_dip_hi + 783);

    auto tr_y_xxyyz_xxxxxx = pbuffer.data(idx_dip_hi + 784);

    auto tr_y_xxyyz_xxxxxy = pbuffer.data(idx_dip_hi + 785);

    auto tr_y_xxyyz_xxxxyy = pbuffer.data(idx_dip_hi + 787);

    auto tr_y_xxyyz_xxxyyy = pbuffer.data(idx_dip_hi + 790);

    auto tr_y_xxyyz_xxyyyy = pbuffer.data(idx_dip_hi + 794);

    auto tr_y_xxyyz_xyyyyy = pbuffer.data(idx_dip_hi + 799);

    auto tr_y_xxyyz_yyyyyz = pbuffer.data(idx_dip_hi + 806);

    auto tr_y_xxyyz_yyyyzz = pbuffer.data(idx_dip_hi + 807);

    auto tr_y_xxyyz_yyyzzz = pbuffer.data(idx_dip_hi + 808);

    auto tr_y_xxyyz_yyzzzz = pbuffer.data(idx_dip_hi + 809);

    auto tr_y_xxyyz_yzzzzz = pbuffer.data(idx_dip_hi + 810);

    auto tr_y_xxyyz_zzzzzz = pbuffer.data(idx_dip_hi + 811);

    auto tr_y_xxyzz_xxxxxy = pbuffer.data(idx_dip_hi + 813);

    auto tr_y_xxyzz_xxxxyy = pbuffer.data(idx_dip_hi + 815);

    auto tr_y_xxyzz_xxxxyz = pbuffer.data(idx_dip_hi + 816);

    auto tr_y_xxyzz_xxxyyy = pbuffer.data(idx_dip_hi + 818);

    auto tr_y_xxyzz_xxxyyz = pbuffer.data(idx_dip_hi + 819);

    auto tr_y_xxyzz_xxxyzz = pbuffer.data(idx_dip_hi + 820);

    auto tr_y_xxyzz_xxyyyy = pbuffer.data(idx_dip_hi + 822);

    auto tr_y_xxyzz_xxyyyz = pbuffer.data(idx_dip_hi + 823);

    auto tr_y_xxyzz_xxyyzz = pbuffer.data(idx_dip_hi + 824);

    auto tr_y_xxyzz_xxyzzz = pbuffer.data(idx_dip_hi + 825);

    auto tr_y_xxyzz_xyyyyy = pbuffer.data(idx_dip_hi + 827);

    auto tr_y_xxyzz_xyyyyz = pbuffer.data(idx_dip_hi + 828);

    auto tr_y_xxyzz_xyyyzz = pbuffer.data(idx_dip_hi + 829);

    auto tr_y_xxyzz_xyyzzz = pbuffer.data(idx_dip_hi + 830);

    auto tr_y_xxyzz_xyzzzz = pbuffer.data(idx_dip_hi + 831);

    auto tr_y_xxyzz_yyyyyy = pbuffer.data(idx_dip_hi + 833);

    auto tr_y_xxyzz_yyyyyz = pbuffer.data(idx_dip_hi + 834);

    auto tr_y_xxyzz_yyyyzz = pbuffer.data(idx_dip_hi + 835);

    auto tr_y_xxyzz_yyyzzz = pbuffer.data(idx_dip_hi + 836);

    auto tr_y_xxyzz_yyzzzz = pbuffer.data(idx_dip_hi + 837);

    auto tr_y_xxyzz_yzzzzz = pbuffer.data(idx_dip_hi + 838);

    auto tr_y_xxyzz_zzzzzz = pbuffer.data(idx_dip_hi + 839);

    auto tr_y_xxzzz_xxxxxx = pbuffer.data(idx_dip_hi + 840);

    auto tr_y_xxzzz_xxxxxy = pbuffer.data(idx_dip_hi + 841);

    auto tr_y_xxzzz_xxxxxz = pbuffer.data(idx_dip_hi + 842);

    auto tr_y_xxzzz_xxxxyy = pbuffer.data(idx_dip_hi + 843);

    auto tr_y_xxzzz_xxxxyz = pbuffer.data(idx_dip_hi + 844);

    auto tr_y_xxzzz_xxxxzz = pbuffer.data(idx_dip_hi + 845);

    auto tr_y_xxzzz_xxxyyy = pbuffer.data(idx_dip_hi + 846);

    auto tr_y_xxzzz_xxxyyz = pbuffer.data(idx_dip_hi + 847);

    auto tr_y_xxzzz_xxxyzz = pbuffer.data(idx_dip_hi + 848);

    auto tr_y_xxzzz_xxxzzz = pbuffer.data(idx_dip_hi + 849);

    auto tr_y_xxzzz_xxyyyy = pbuffer.data(idx_dip_hi + 850);

    auto tr_y_xxzzz_xxyyyz = pbuffer.data(idx_dip_hi + 851);

    auto tr_y_xxzzz_xxyyzz = pbuffer.data(idx_dip_hi + 852);

    auto tr_y_xxzzz_xxyzzz = pbuffer.data(idx_dip_hi + 853);

    auto tr_y_xxzzz_xxzzzz = pbuffer.data(idx_dip_hi + 854);

    auto tr_y_xxzzz_xyyyyy = pbuffer.data(idx_dip_hi + 855);

    auto tr_y_xxzzz_xyyyyz = pbuffer.data(idx_dip_hi + 856);

    auto tr_y_xxzzz_xyyyzz = pbuffer.data(idx_dip_hi + 857);

    auto tr_y_xxzzz_xyyzzz = pbuffer.data(idx_dip_hi + 858);

    auto tr_y_xxzzz_xyzzzz = pbuffer.data(idx_dip_hi + 859);

    auto tr_y_xxzzz_xzzzzz = pbuffer.data(idx_dip_hi + 860);

    auto tr_y_xxzzz_yyyyyy = pbuffer.data(idx_dip_hi + 861);

    auto tr_y_xxzzz_yyyyyz = pbuffer.data(idx_dip_hi + 862);

    auto tr_y_xxzzz_yyyyzz = pbuffer.data(idx_dip_hi + 863);

    auto tr_y_xxzzz_yyyzzz = pbuffer.data(idx_dip_hi + 864);

    auto tr_y_xxzzz_yyzzzz = pbuffer.data(idx_dip_hi + 865);

    auto tr_y_xxzzz_yzzzzz = pbuffer.data(idx_dip_hi + 866);

    auto tr_y_xxzzz_zzzzzz = pbuffer.data(idx_dip_hi + 867);

    auto tr_y_xyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 868);

    auto tr_y_xyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 869);

    auto tr_y_xyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 870);

    auto tr_y_xyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 871);

    auto tr_y_xyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 872);

    auto tr_y_xyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 873);

    auto tr_y_xyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 874);

    auto tr_y_xyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 875);

    auto tr_y_xyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 876);

    auto tr_y_xyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 877);

    auto tr_y_xyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 878);

    auto tr_y_xyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 879);

    auto tr_y_xyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 880);

    auto tr_y_xyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 881);

    auto tr_y_xyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 882);

    auto tr_y_xyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 883);

    auto tr_y_xyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 884);

    auto tr_y_xyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 885);

    auto tr_y_xyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 886);

    auto tr_y_xyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 887);

    auto tr_y_xyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 888);

    auto tr_y_xyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 889);

    auto tr_y_xyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 890);

    auto tr_y_xyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 891);

    auto tr_y_xyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 892);

    auto tr_y_xyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 893);

    auto tr_y_xyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 894);

    auto tr_y_xyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 895);

    auto tr_y_xyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 918);

    auto tr_y_xyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 919);

    auto tr_y_xyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 920);

    auto tr_y_xyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 921);

    auto tr_y_xyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 922);

    auto tr_y_xyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 923);

    auto tr_y_xyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 926);

    auto tr_y_xyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 928);

    auto tr_y_xyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 929);

    auto tr_y_xyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 931);

    auto tr_y_xyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 932);

    auto tr_y_xyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 933);

    auto tr_y_xyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 935);

    auto tr_y_xyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 936);

    auto tr_y_xyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 937);

    auto tr_y_xyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 938);

    auto tr_y_xyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 940);

    auto tr_y_xyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 941);

    auto tr_y_xyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 942);

    auto tr_y_xyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 943);

    auto tr_y_xyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 944);

    auto tr_y_xyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 945);

    auto tr_y_xyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 946);

    auto tr_y_xyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 947);

    auto tr_y_xyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 948);

    auto tr_y_xyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 949);

    auto tr_y_xyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 950);

    auto tr_y_xyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 951);

    auto tr_y_xyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 956);

    auto tr_y_xyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 959);

    auto tr_y_xyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 960);

    auto tr_y_xyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 963);

    auto tr_y_xyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 964);

    auto tr_y_xyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 965);

    auto tr_y_xyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 968);

    auto tr_y_xyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 969);

    auto tr_y_xyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 970);

    auto tr_y_xyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 971);

    auto tr_y_xyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 973);

    auto tr_y_xyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 974);

    auto tr_y_xyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 975);

    auto tr_y_xyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 976);

    auto tr_y_xyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 977);

    auto tr_y_xyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 978);

    auto tr_y_xyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 979);

    auto tr_y_xzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 982);

    auto tr_y_xzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 984);

    auto tr_y_xzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 985);

    auto tr_y_xzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 987);

    auto tr_y_xzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 988);

    auto tr_y_xzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 989);

    auto tr_y_xzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 991);

    auto tr_y_xzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 992);

    auto tr_y_xzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 993);

    auto tr_y_xzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 994);

    auto tr_y_xzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 996);

    auto tr_y_xzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 997);

    auto tr_y_xzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 998);

    auto tr_y_xzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 999);

    auto tr_y_xzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1000);

    auto tr_y_xzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1001);

    auto tr_y_xzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1002);

    auto tr_y_xzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1003);

    auto tr_y_xzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1004);

    auto tr_y_xzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1005);

    auto tr_y_xzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1006);

    auto tr_y_xzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1007);

    auto tr_y_yyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1008);

    auto tr_y_yyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1009);

    auto tr_y_yyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1010);

    auto tr_y_yyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1011);

    auto tr_y_yyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1012);

    auto tr_y_yyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1013);

    auto tr_y_yyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1014);

    auto tr_y_yyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1015);

    auto tr_y_yyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1016);

    auto tr_y_yyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1017);

    auto tr_y_yyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1018);

    auto tr_y_yyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1019);

    auto tr_y_yyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1020);

    auto tr_y_yyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1021);

    auto tr_y_yyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1022);

    auto tr_y_yyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1023);

    auto tr_y_yyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1024);

    auto tr_y_yyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1025);

    auto tr_y_yyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1026);

    auto tr_y_yyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1027);

    auto tr_y_yyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1028);

    auto tr_y_yyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1029);

    auto tr_y_yyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1030);

    auto tr_y_yyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1031);

    auto tr_y_yyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1032);

    auto tr_y_yyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1033);

    auto tr_y_yyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1034);

    auto tr_y_yyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1035);

    auto tr_y_yyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1036);

    auto tr_y_yyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 1037);

    auto tr_y_yyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1038);

    auto tr_y_yyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 1039);

    auto tr_y_yyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1040);

    auto tr_y_yyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1041);

    auto tr_y_yyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 1042);

    auto tr_y_yyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1043);

    auto tr_y_yyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1044);

    auto tr_y_yyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1045);

    auto tr_y_yyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 1046);

    auto tr_y_yyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1047);

    auto tr_y_yyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1048);

    auto tr_y_yyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1049);

    auto tr_y_yyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1050);

    auto tr_y_yyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 1051);

    auto tr_y_yyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1052);

    auto tr_y_yyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1053);

    auto tr_y_yyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1054);

    auto tr_y_yyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1055);

    auto tr_y_yyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1056);

    auto tr_y_yyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1057);

    auto tr_y_yyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1058);

    auto tr_y_yyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1059);

    auto tr_y_yyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1060);

    auto tr_y_yyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1061);

    auto tr_y_yyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1062);

    auto tr_y_yyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1063);

    auto tr_y_yyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1064);

    auto tr_y_yyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1065);

    auto tr_y_yyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1066);

    auto tr_y_yyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1067);

    auto tr_y_yyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1068);

    auto tr_y_yyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1069);

    auto tr_y_yyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1070);

    auto tr_y_yyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1071);

    auto tr_y_yyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1072);

    auto tr_y_yyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1073);

    auto tr_y_yyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1074);

    auto tr_y_yyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1075);

    auto tr_y_yyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1076);

    auto tr_y_yyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1077);

    auto tr_y_yyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1078);

    auto tr_y_yyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1079);

    auto tr_y_yyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1080);

    auto tr_y_yyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1081);

    auto tr_y_yyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1082);

    auto tr_y_yyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1083);

    auto tr_y_yyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1084);

    auto tr_y_yyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1085);

    auto tr_y_yyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1086);

    auto tr_y_yyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1087);

    auto tr_y_yyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1088);

    auto tr_y_yyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1089);

    auto tr_y_yyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1090);

    auto tr_y_yyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1091);

    auto tr_y_yyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1092);

    auto tr_y_yyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1093);

    auto tr_y_yyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1094);

    auto tr_y_yyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1095);

    auto tr_y_yyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1096);

    auto tr_y_yyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1097);

    auto tr_y_yyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1098);

    auto tr_y_yyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1099);

    auto tr_y_yyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1100);

    auto tr_y_yyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1101);

    auto tr_y_yyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1102);

    auto tr_y_yyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1103);

    auto tr_y_yyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1104);

    auto tr_y_yyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1105);

    auto tr_y_yyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1106);

    auto tr_y_yyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1107);

    auto tr_y_yyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1108);

    auto tr_y_yyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1109);

    auto tr_y_yyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1110);

    auto tr_y_yyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1111);

    auto tr_y_yyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1112);

    auto tr_y_yyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1113);

    auto tr_y_yyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1114);

    auto tr_y_yyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1115);

    auto tr_y_yyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1116);

    auto tr_y_yyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1117);

    auto tr_y_yyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1118);

    auto tr_y_yyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1119);

    auto tr_y_yzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1120);

    auto tr_y_yzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1121);

    auto tr_y_yzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1122);

    auto tr_y_yzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1123);

    auto tr_y_yzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1124);

    auto tr_y_yzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1125);

    auto tr_y_yzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1126);

    auto tr_y_yzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1127);

    auto tr_y_yzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1128);

    auto tr_y_yzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1129);

    auto tr_y_yzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1130);

    auto tr_y_yzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1131);

    auto tr_y_yzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1132);

    auto tr_y_yzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1133);

    auto tr_y_yzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1134);

    auto tr_y_yzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1135);

    auto tr_y_yzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1136);

    auto tr_y_yzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1137);

    auto tr_y_yzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1138);

    auto tr_y_yzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1139);

    auto tr_y_yzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1140);

    auto tr_y_yzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1141);

    auto tr_y_yzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1142);

    auto tr_y_yzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1143);

    auto tr_y_yzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1144);

    auto tr_y_yzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1145);

    auto tr_y_yzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1146);

    auto tr_y_yzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1147);

    auto tr_y_zzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1148);

    auto tr_y_zzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1149);

    auto tr_y_zzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1150);

    auto tr_y_zzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1151);

    auto tr_y_zzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1152);

    auto tr_y_zzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1153);

    auto tr_y_zzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1154);

    auto tr_y_zzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1155);

    auto tr_y_zzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1156);

    auto tr_y_zzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1157);

    auto tr_y_zzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1158);

    auto tr_y_zzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1159);

    auto tr_y_zzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1160);

    auto tr_y_zzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1161);

    auto tr_y_zzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1162);

    auto tr_y_zzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1163);

    auto tr_y_zzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1164);

    auto tr_y_zzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1165);

    auto tr_y_zzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1166);

    auto tr_y_zzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1167);

    auto tr_y_zzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1168);

    auto tr_y_zzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1169);

    auto tr_y_zzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1170);

    auto tr_y_zzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1171);

    auto tr_y_zzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1172);

    auto tr_y_zzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1173);

    auto tr_y_zzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1174);

    auto tr_y_zzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1175);

    auto tr_z_xxxxx_xxxxxx = pbuffer.data(idx_dip_hi + 1176);

    auto tr_z_xxxxx_xxxxxy = pbuffer.data(idx_dip_hi + 1177);

    auto tr_z_xxxxx_xxxxxz = pbuffer.data(idx_dip_hi + 1178);

    auto tr_z_xxxxx_xxxxyy = pbuffer.data(idx_dip_hi + 1179);

    auto tr_z_xxxxx_xxxxyz = pbuffer.data(idx_dip_hi + 1180);

    auto tr_z_xxxxx_xxxxzz = pbuffer.data(idx_dip_hi + 1181);

    auto tr_z_xxxxx_xxxyyy = pbuffer.data(idx_dip_hi + 1182);

    auto tr_z_xxxxx_xxxyyz = pbuffer.data(idx_dip_hi + 1183);

    auto tr_z_xxxxx_xxxyzz = pbuffer.data(idx_dip_hi + 1184);

    auto tr_z_xxxxx_xxxzzz = pbuffer.data(idx_dip_hi + 1185);

    auto tr_z_xxxxx_xxyyyy = pbuffer.data(idx_dip_hi + 1186);

    auto tr_z_xxxxx_xxyyyz = pbuffer.data(idx_dip_hi + 1187);

    auto tr_z_xxxxx_xxyyzz = pbuffer.data(idx_dip_hi + 1188);

    auto tr_z_xxxxx_xxyzzz = pbuffer.data(idx_dip_hi + 1189);

    auto tr_z_xxxxx_xxzzzz = pbuffer.data(idx_dip_hi + 1190);

    auto tr_z_xxxxx_xyyyyy = pbuffer.data(idx_dip_hi + 1191);

    auto tr_z_xxxxx_xyyyyz = pbuffer.data(idx_dip_hi + 1192);

    auto tr_z_xxxxx_xyyyzz = pbuffer.data(idx_dip_hi + 1193);

    auto tr_z_xxxxx_xyyzzz = pbuffer.data(idx_dip_hi + 1194);

    auto tr_z_xxxxx_xyzzzz = pbuffer.data(idx_dip_hi + 1195);

    auto tr_z_xxxxx_xzzzzz = pbuffer.data(idx_dip_hi + 1196);

    auto tr_z_xxxxx_yyyyyy = pbuffer.data(idx_dip_hi + 1197);

    auto tr_z_xxxxx_yyyyyz = pbuffer.data(idx_dip_hi + 1198);

    auto tr_z_xxxxx_yyyyzz = pbuffer.data(idx_dip_hi + 1199);

    auto tr_z_xxxxx_yyyzzz = pbuffer.data(idx_dip_hi + 1200);

    auto tr_z_xxxxx_yyzzzz = pbuffer.data(idx_dip_hi + 1201);

    auto tr_z_xxxxx_yzzzzz = pbuffer.data(idx_dip_hi + 1202);

    auto tr_z_xxxxx_zzzzzz = pbuffer.data(idx_dip_hi + 1203);

    auto tr_z_xxxxy_xxxxxx = pbuffer.data(idx_dip_hi + 1204);

    auto tr_z_xxxxy_xxxxxz = pbuffer.data(idx_dip_hi + 1206);

    auto tr_z_xxxxy_xxxxzz = pbuffer.data(idx_dip_hi + 1209);

    auto tr_z_xxxxy_xxxzzz = pbuffer.data(idx_dip_hi + 1213);

    auto tr_z_xxxxy_xxzzzz = pbuffer.data(idx_dip_hi + 1218);

    auto tr_z_xxxxy_xzzzzz = pbuffer.data(idx_dip_hi + 1224);

    auto tr_z_xxxxy_yyyyyy = pbuffer.data(idx_dip_hi + 1225);

    auto tr_z_xxxxy_yyyyyz = pbuffer.data(idx_dip_hi + 1226);

    auto tr_z_xxxxy_yyyyzz = pbuffer.data(idx_dip_hi + 1227);

    auto tr_z_xxxxy_yyyzzz = pbuffer.data(idx_dip_hi + 1228);

    auto tr_z_xxxxy_yyzzzz = pbuffer.data(idx_dip_hi + 1229);

    auto tr_z_xxxxy_yzzzzz = pbuffer.data(idx_dip_hi + 1230);

    auto tr_z_xxxxz_xxxxxx = pbuffer.data(idx_dip_hi + 1232);

    auto tr_z_xxxxz_xxxxxy = pbuffer.data(idx_dip_hi + 1233);

    auto tr_z_xxxxz_xxxxxz = pbuffer.data(idx_dip_hi + 1234);

    auto tr_z_xxxxz_xxxxyy = pbuffer.data(idx_dip_hi + 1235);

    auto tr_z_xxxxz_xxxxyz = pbuffer.data(idx_dip_hi + 1236);

    auto tr_z_xxxxz_xxxxzz = pbuffer.data(idx_dip_hi + 1237);

    auto tr_z_xxxxz_xxxyyy = pbuffer.data(idx_dip_hi + 1238);

    auto tr_z_xxxxz_xxxyyz = pbuffer.data(idx_dip_hi + 1239);

    auto tr_z_xxxxz_xxxyzz = pbuffer.data(idx_dip_hi + 1240);

    auto tr_z_xxxxz_xxxzzz = pbuffer.data(idx_dip_hi + 1241);

    auto tr_z_xxxxz_xxyyyy = pbuffer.data(idx_dip_hi + 1242);

    auto tr_z_xxxxz_xxyyyz = pbuffer.data(idx_dip_hi + 1243);

    auto tr_z_xxxxz_xxyyzz = pbuffer.data(idx_dip_hi + 1244);

    auto tr_z_xxxxz_xxyzzz = pbuffer.data(idx_dip_hi + 1245);

    auto tr_z_xxxxz_xxzzzz = pbuffer.data(idx_dip_hi + 1246);

    auto tr_z_xxxxz_xyyyyy = pbuffer.data(idx_dip_hi + 1247);

    auto tr_z_xxxxz_xyyyyz = pbuffer.data(idx_dip_hi + 1248);

    auto tr_z_xxxxz_xyyyzz = pbuffer.data(idx_dip_hi + 1249);

    auto tr_z_xxxxz_xyyzzz = pbuffer.data(idx_dip_hi + 1250);

    auto tr_z_xxxxz_xyzzzz = pbuffer.data(idx_dip_hi + 1251);

    auto tr_z_xxxxz_xzzzzz = pbuffer.data(idx_dip_hi + 1252);

    auto tr_z_xxxxz_yyyyyy = pbuffer.data(idx_dip_hi + 1253);

    auto tr_z_xxxxz_yyyyyz = pbuffer.data(idx_dip_hi + 1254);

    auto tr_z_xxxxz_yyyyzz = pbuffer.data(idx_dip_hi + 1255);

    auto tr_z_xxxxz_yyyzzz = pbuffer.data(idx_dip_hi + 1256);

    auto tr_z_xxxxz_yyzzzz = pbuffer.data(idx_dip_hi + 1257);

    auto tr_z_xxxxz_yzzzzz = pbuffer.data(idx_dip_hi + 1258);

    auto tr_z_xxxxz_zzzzzz = pbuffer.data(idx_dip_hi + 1259);

    auto tr_z_xxxyy_xxxxxx = pbuffer.data(idx_dip_hi + 1260);

    auto tr_z_xxxyy_xxxxxy = pbuffer.data(idx_dip_hi + 1261);

    auto tr_z_xxxyy_xxxxxz = pbuffer.data(idx_dip_hi + 1262);

    auto tr_z_xxxyy_xxxxyy = pbuffer.data(idx_dip_hi + 1263);

    auto tr_z_xxxyy_xxxxyz = pbuffer.data(idx_dip_hi + 1264);

    auto tr_z_xxxyy_xxxxzz = pbuffer.data(idx_dip_hi + 1265);

    auto tr_z_xxxyy_xxxyyy = pbuffer.data(idx_dip_hi + 1266);

    auto tr_z_xxxyy_xxxyyz = pbuffer.data(idx_dip_hi + 1267);

    auto tr_z_xxxyy_xxxyzz = pbuffer.data(idx_dip_hi + 1268);

    auto tr_z_xxxyy_xxxzzz = pbuffer.data(idx_dip_hi + 1269);

    auto tr_z_xxxyy_xxyyyy = pbuffer.data(idx_dip_hi + 1270);

    auto tr_z_xxxyy_xxyyyz = pbuffer.data(idx_dip_hi + 1271);

    auto tr_z_xxxyy_xxyyzz = pbuffer.data(idx_dip_hi + 1272);

    auto tr_z_xxxyy_xxyzzz = pbuffer.data(idx_dip_hi + 1273);

    auto tr_z_xxxyy_xxzzzz = pbuffer.data(idx_dip_hi + 1274);

    auto tr_z_xxxyy_xyyyyy = pbuffer.data(idx_dip_hi + 1275);

    auto tr_z_xxxyy_xyyyyz = pbuffer.data(idx_dip_hi + 1276);

    auto tr_z_xxxyy_xyyyzz = pbuffer.data(idx_dip_hi + 1277);

    auto tr_z_xxxyy_xyyzzz = pbuffer.data(idx_dip_hi + 1278);

    auto tr_z_xxxyy_xyzzzz = pbuffer.data(idx_dip_hi + 1279);

    auto tr_z_xxxyy_xzzzzz = pbuffer.data(idx_dip_hi + 1280);

    auto tr_z_xxxyy_yyyyyy = pbuffer.data(idx_dip_hi + 1281);

    auto tr_z_xxxyy_yyyyyz = pbuffer.data(idx_dip_hi + 1282);

    auto tr_z_xxxyy_yyyyzz = pbuffer.data(idx_dip_hi + 1283);

    auto tr_z_xxxyy_yyyzzz = pbuffer.data(idx_dip_hi + 1284);

    auto tr_z_xxxyy_yyzzzz = pbuffer.data(idx_dip_hi + 1285);

    auto tr_z_xxxyy_yzzzzz = pbuffer.data(idx_dip_hi + 1286);

    auto tr_z_xxxyy_zzzzzz = pbuffer.data(idx_dip_hi + 1287);

    auto tr_z_xxxyz_xxxxxx = pbuffer.data(idx_dip_hi + 1288);

    auto tr_z_xxxyz_xxxxxz = pbuffer.data(idx_dip_hi + 1290);

    auto tr_z_xxxyz_xxxxzz = pbuffer.data(idx_dip_hi + 1293);

    auto tr_z_xxxyz_xxxzzz = pbuffer.data(idx_dip_hi + 1297);

    auto tr_z_xxxyz_xxzzzz = pbuffer.data(idx_dip_hi + 1302);

    auto tr_z_xxxyz_xzzzzz = pbuffer.data(idx_dip_hi + 1308);

    auto tr_z_xxxyz_yyyyyy = pbuffer.data(idx_dip_hi + 1309);

    auto tr_z_xxxyz_yyyyyz = pbuffer.data(idx_dip_hi + 1310);

    auto tr_z_xxxyz_yyyyzz = pbuffer.data(idx_dip_hi + 1311);

    auto tr_z_xxxyz_yyyzzz = pbuffer.data(idx_dip_hi + 1312);

    auto tr_z_xxxyz_yyzzzz = pbuffer.data(idx_dip_hi + 1313);

    auto tr_z_xxxyz_yzzzzz = pbuffer.data(idx_dip_hi + 1314);

    auto tr_z_xxxzz_xxxxxx = pbuffer.data(idx_dip_hi + 1316);

    auto tr_z_xxxzz_xxxxxy = pbuffer.data(idx_dip_hi + 1317);

    auto tr_z_xxxzz_xxxxxz = pbuffer.data(idx_dip_hi + 1318);

    auto tr_z_xxxzz_xxxxyy = pbuffer.data(idx_dip_hi + 1319);

    auto tr_z_xxxzz_xxxxyz = pbuffer.data(idx_dip_hi + 1320);

    auto tr_z_xxxzz_xxxxzz = pbuffer.data(idx_dip_hi + 1321);

    auto tr_z_xxxzz_xxxyyy = pbuffer.data(idx_dip_hi + 1322);

    auto tr_z_xxxzz_xxxyyz = pbuffer.data(idx_dip_hi + 1323);

    auto tr_z_xxxzz_xxxyzz = pbuffer.data(idx_dip_hi + 1324);

    auto tr_z_xxxzz_xxxzzz = pbuffer.data(idx_dip_hi + 1325);

    auto tr_z_xxxzz_xxyyyy = pbuffer.data(idx_dip_hi + 1326);

    auto tr_z_xxxzz_xxyyyz = pbuffer.data(idx_dip_hi + 1327);

    auto tr_z_xxxzz_xxyyzz = pbuffer.data(idx_dip_hi + 1328);

    auto tr_z_xxxzz_xxyzzz = pbuffer.data(idx_dip_hi + 1329);

    auto tr_z_xxxzz_xxzzzz = pbuffer.data(idx_dip_hi + 1330);

    auto tr_z_xxxzz_xyyyyy = pbuffer.data(idx_dip_hi + 1331);

    auto tr_z_xxxzz_xyyyyz = pbuffer.data(idx_dip_hi + 1332);

    auto tr_z_xxxzz_xyyyzz = pbuffer.data(idx_dip_hi + 1333);

    auto tr_z_xxxzz_xyyzzz = pbuffer.data(idx_dip_hi + 1334);

    auto tr_z_xxxzz_xyzzzz = pbuffer.data(idx_dip_hi + 1335);

    auto tr_z_xxxzz_xzzzzz = pbuffer.data(idx_dip_hi + 1336);

    auto tr_z_xxxzz_yyyyyy = pbuffer.data(idx_dip_hi + 1337);

    auto tr_z_xxxzz_yyyyyz = pbuffer.data(idx_dip_hi + 1338);

    auto tr_z_xxxzz_yyyyzz = pbuffer.data(idx_dip_hi + 1339);

    auto tr_z_xxxzz_yyyzzz = pbuffer.data(idx_dip_hi + 1340);

    auto tr_z_xxxzz_yyzzzz = pbuffer.data(idx_dip_hi + 1341);

    auto tr_z_xxxzz_yzzzzz = pbuffer.data(idx_dip_hi + 1342);

    auto tr_z_xxxzz_zzzzzz = pbuffer.data(idx_dip_hi + 1343);

    auto tr_z_xxyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1344);

    auto tr_z_xxyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1345);

    auto tr_z_xxyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1346);

    auto tr_z_xxyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1347);

    auto tr_z_xxyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1348);

    auto tr_z_xxyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1349);

    auto tr_z_xxyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1350);

    auto tr_z_xxyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1351);

    auto tr_z_xxyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1352);

    auto tr_z_xxyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1353);

    auto tr_z_xxyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1354);

    auto tr_z_xxyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1355);

    auto tr_z_xxyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1356);

    auto tr_z_xxyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1357);

    auto tr_z_xxyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1358);

    auto tr_z_xxyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1359);

    auto tr_z_xxyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1360);

    auto tr_z_xxyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1361);

    auto tr_z_xxyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1362);

    auto tr_z_xxyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1363);

    auto tr_z_xxyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1364);

    auto tr_z_xxyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1365);

    auto tr_z_xxyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1366);

    auto tr_z_xxyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1367);

    auto tr_z_xxyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1368);

    auto tr_z_xxyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1369);

    auto tr_z_xxyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1370);

    auto tr_z_xxyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1371);

    auto tr_z_xxyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1372);

    auto tr_z_xxyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1374);

    auto tr_z_xxyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1376);

    auto tr_z_xxyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1377);

    auto tr_z_xxyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1379);

    auto tr_z_xxyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1380);

    auto tr_z_xxyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1381);

    auto tr_z_xxyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1383);

    auto tr_z_xxyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1384);

    auto tr_z_xxyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1385);

    auto tr_z_xxyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1386);

    auto tr_z_xxyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1388);

    auto tr_z_xxyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1389);

    auto tr_z_xxyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1390);

    auto tr_z_xxyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1391);

    auto tr_z_xxyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1392);

    auto tr_z_xxyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1393);

    auto tr_z_xxyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1394);

    auto tr_z_xxyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1395);

    auto tr_z_xxyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1396);

    auto tr_z_xxyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1397);

    auto tr_z_xxyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1398);

    auto tr_z_xxyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1399);

    auto tr_z_xxyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1400);

    auto tr_z_xxyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1402);

    auto tr_z_xxyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1405);

    auto tr_z_xxyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1409);

    auto tr_z_xxyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1414);

    auto tr_z_xxyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1420);

    auto tr_z_xxyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1421);

    auto tr_z_xxyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1422);

    auto tr_z_xxyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1423);

    auto tr_z_xxyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1424);

    auto tr_z_xxyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1425);

    auto tr_z_xxyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1426);

    auto tr_z_xxzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1428);

    auto tr_z_xxzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1429);

    auto tr_z_xxzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1430);

    auto tr_z_xxzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1431);

    auto tr_z_xxzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1432);

    auto tr_z_xxzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1433);

    auto tr_z_xxzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1434);

    auto tr_z_xxzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1435);

    auto tr_z_xxzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1436);

    auto tr_z_xxzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1437);

    auto tr_z_xxzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1438);

    auto tr_z_xxzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1439);

    auto tr_z_xxzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1440);

    auto tr_z_xxzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1441);

    auto tr_z_xxzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1442);

    auto tr_z_xxzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1443);

    auto tr_z_xxzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1444);

    auto tr_z_xxzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1445);

    auto tr_z_xxzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1446);

    auto tr_z_xxzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1447);

    auto tr_z_xxzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1448);

    auto tr_z_xxzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1449);

    auto tr_z_xxzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1450);

    auto tr_z_xxzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1451);

    auto tr_z_xxzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1452);

    auto tr_z_xxzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1453);

    auto tr_z_xxzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1454);

    auto tr_z_xxzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1455);

    auto tr_z_xyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1457);

    auto tr_z_xyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1459);

    auto tr_z_xyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1460);

    auto tr_z_xyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1462);

    auto tr_z_xyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1463);

    auto tr_z_xyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1464);

    auto tr_z_xyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1466);

    auto tr_z_xyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1467);

    auto tr_z_xyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1468);

    auto tr_z_xyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1469);

    auto tr_z_xyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1471);

    auto tr_z_xyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1472);

    auto tr_z_xyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1473);

    auto tr_z_xyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1474);

    auto tr_z_xyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1475);

    auto tr_z_xyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1477);

    auto tr_z_xyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1478);

    auto tr_z_xyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1479);

    auto tr_z_xyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1480);

    auto tr_z_xyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1481);

    auto tr_z_xyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1482);

    auto tr_z_xyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1483);

    auto tr_z_xyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1488);

    auto tr_z_xyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1491);

    auto tr_z_xyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1492);

    auto tr_z_xyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1495);

    auto tr_z_xyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1496);

    auto tr_z_xyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1497);

    auto tr_z_xyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1500);

    auto tr_z_xyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1501);

    auto tr_z_xyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1502);

    auto tr_z_xyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1503);

    auto tr_z_xyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1505);

    auto tr_z_xyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1506);

    auto tr_z_xyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1507);

    auto tr_z_xyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1508);

    auto tr_z_xyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1509);

    auto tr_z_xyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1510);

    auto tr_z_xyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1511);

    auto tr_z_xyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1513);

    auto tr_z_xyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1515);

    auto tr_z_xyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1516);

    auto tr_z_xyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1518);

    auto tr_z_xyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1519);

    auto tr_z_xyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1520);

    auto tr_z_xyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1522);

    auto tr_z_xyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1523);

    auto tr_z_xyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1524);

    auto tr_z_xyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1525);

    auto tr_z_xyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1527);

    auto tr_z_xyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1528);

    auto tr_z_xyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1529);

    auto tr_z_xyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1530);

    auto tr_z_xyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1531);

    auto tr_z_xyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1533);

    auto tr_z_xyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1534);

    auto tr_z_xyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1535);

    auto tr_z_xyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1536);

    auto tr_z_xyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1537);

    auto tr_z_xyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1538);

    auto tr_z_xyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1539);

    auto tr_z_xyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1561);

    auto tr_z_xyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1562);

    auto tr_z_xyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1563);

    auto tr_z_xyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1564);

    auto tr_z_xyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1565);

    auto tr_z_xyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1566);

    auto tr_z_xzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1568);

    auto tr_z_xzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1569);

    auto tr_z_xzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1570);

    auto tr_z_xzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1571);

    auto tr_z_xzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1572);

    auto tr_z_xzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1573);

    auto tr_z_xzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1574);

    auto tr_z_xzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1575);

    auto tr_z_xzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1576);

    auto tr_z_xzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1577);

    auto tr_z_xzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1578);

    auto tr_z_xzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1579);

    auto tr_z_xzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1580);

    auto tr_z_xzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1581);

    auto tr_z_xzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1582);

    auto tr_z_xzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1583);

    auto tr_z_xzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1584);

    auto tr_z_xzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1585);

    auto tr_z_xzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1586);

    auto tr_z_xzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1587);

    auto tr_z_xzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1588);

    auto tr_z_xzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1589);

    auto tr_z_xzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1590);

    auto tr_z_xzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1591);

    auto tr_z_xzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1592);

    auto tr_z_xzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1593);

    auto tr_z_xzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1594);

    auto tr_z_xzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1595);

    auto tr_z_yyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1596);

    auto tr_z_yyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1597);

    auto tr_z_yyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1598);

    auto tr_z_yyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1599);

    auto tr_z_yyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1600);

    auto tr_z_yyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1601);

    auto tr_z_yyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1602);

    auto tr_z_yyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1603);

    auto tr_z_yyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1604);

    auto tr_z_yyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1605);

    auto tr_z_yyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1606);

    auto tr_z_yyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1607);

    auto tr_z_yyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1608);

    auto tr_z_yyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1609);

    auto tr_z_yyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1610);

    auto tr_z_yyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1611);

    auto tr_z_yyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1612);

    auto tr_z_yyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1613);

    auto tr_z_yyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1614);

    auto tr_z_yyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1615);

    auto tr_z_yyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1616);

    auto tr_z_yyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1617);

    auto tr_z_yyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1618);

    auto tr_z_yyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1619);

    auto tr_z_yyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1620);

    auto tr_z_yyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1621);

    auto tr_z_yyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1622);

    auto tr_z_yyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1623);

    auto tr_z_yyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1624);

    auto tr_z_yyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 1625);

    auto tr_z_yyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1626);

    auto tr_z_yyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 1627);

    auto tr_z_yyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1628);

    auto tr_z_yyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1629);

    auto tr_z_yyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 1630);

    auto tr_z_yyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1631);

    auto tr_z_yyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1632);

    auto tr_z_yyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1633);

    auto tr_z_yyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 1634);

    auto tr_z_yyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1635);

    auto tr_z_yyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1636);

    auto tr_z_yyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1637);

    auto tr_z_yyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1638);

    auto tr_z_yyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 1639);

    auto tr_z_yyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1640);

    auto tr_z_yyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1641);

    auto tr_z_yyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1642);

    auto tr_z_yyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1643);

    auto tr_z_yyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1644);

    auto tr_z_yyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1645);

    auto tr_z_yyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1646);

    auto tr_z_yyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1647);

    auto tr_z_yyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1648);

    auto tr_z_yyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1649);

    auto tr_z_yyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1650);

    auto tr_z_yyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1651);

    auto tr_z_yyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1652);

    auto tr_z_yyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1653);

    auto tr_z_yyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1654);

    auto tr_z_yyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1655);

    auto tr_z_yyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1656);

    auto tr_z_yyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1657);

    auto tr_z_yyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1658);

    auto tr_z_yyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1659);

    auto tr_z_yyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1660);

    auto tr_z_yyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1661);

    auto tr_z_yyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1662);

    auto tr_z_yyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1663);

    auto tr_z_yyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1664);

    auto tr_z_yyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1665);

    auto tr_z_yyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1666);

    auto tr_z_yyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1667);

    auto tr_z_yyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1668);

    auto tr_z_yyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1669);

    auto tr_z_yyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1670);

    auto tr_z_yyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1671);

    auto tr_z_yyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1672);

    auto tr_z_yyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1673);

    auto tr_z_yyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1674);

    auto tr_z_yyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1675);

    auto tr_z_yyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1676);

    auto tr_z_yyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1677);

    auto tr_z_yyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1678);

    auto tr_z_yyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1679);

    auto tr_z_yyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1680);

    auto tr_z_yyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1681);

    auto tr_z_yyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1682);

    auto tr_z_yyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1683);

    auto tr_z_yyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1684);

    auto tr_z_yyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1685);

    auto tr_z_yyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1686);

    auto tr_z_yyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1687);

    auto tr_z_yyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1688);

    auto tr_z_yyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1689);

    auto tr_z_yyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1690);

    auto tr_z_yyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1691);

    auto tr_z_yyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1692);

    auto tr_z_yyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1693);

    auto tr_z_yyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1694);

    auto tr_z_yyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1695);

    auto tr_z_yyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1696);

    auto tr_z_yyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1697);

    auto tr_z_yyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1698);

    auto tr_z_yyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1699);

    auto tr_z_yyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1700);

    auto tr_z_yyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1701);

    auto tr_z_yyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1702);

    auto tr_z_yyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1703);

    auto tr_z_yyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1704);

    auto tr_z_yyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1705);

    auto tr_z_yyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1706);

    auto tr_z_yyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1707);

    auto tr_z_yzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1708);

    auto tr_z_yzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1709);

    auto tr_z_yzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1710);

    auto tr_z_yzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1711);

    auto tr_z_yzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1712);

    auto tr_z_yzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1713);

    auto tr_z_yzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1714);

    auto tr_z_yzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1715);

    auto tr_z_yzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1716);

    auto tr_z_yzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1717);

    auto tr_z_yzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1718);

    auto tr_z_yzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1719);

    auto tr_z_yzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1720);

    auto tr_z_yzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1721);

    auto tr_z_yzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1722);

    auto tr_z_yzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1723);

    auto tr_z_yzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1724);

    auto tr_z_yzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1725);

    auto tr_z_yzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1726);

    auto tr_z_yzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1727);

    auto tr_z_yzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1728);

    auto tr_z_yzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1729);

    auto tr_z_yzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1730);

    auto tr_z_yzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1731);

    auto tr_z_yzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1732);

    auto tr_z_yzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1733);

    auto tr_z_yzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1734);

    auto tr_z_yzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1735);

    auto tr_z_zzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1736);

    auto tr_z_zzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1737);

    auto tr_z_zzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1738);

    auto tr_z_zzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1739);

    auto tr_z_zzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1740);

    auto tr_z_zzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1741);

    auto tr_z_zzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1742);

    auto tr_z_zzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1743);

    auto tr_z_zzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1744);

    auto tr_z_zzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1745);

    auto tr_z_zzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1746);

    auto tr_z_zzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1747);

    auto tr_z_zzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1748);

    auto tr_z_zzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1749);

    auto tr_z_zzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1750);

    auto tr_z_zzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1751);

    auto tr_z_zzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1752);

    auto tr_z_zzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1753);

    auto tr_z_zzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1754);

    auto tr_z_zzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1755);

    auto tr_z_zzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1756);

    auto tr_z_zzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1757);

    auto tr_z_zzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1758);

    auto tr_z_zzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1759);

    auto tr_z_zzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1760);

    auto tr_z_zzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1761);

    auto tr_z_zzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1762);

    auto tr_z_zzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1763);

    // Set up 0-28 components of targeted buffer : II

    auto tr_x_xxxxxx_xxxxxx = pbuffer.data(idx_dip_ii);

    auto tr_x_xxxxxx_xxxxxy = pbuffer.data(idx_dip_ii + 1);

    auto tr_x_xxxxxx_xxxxxz = pbuffer.data(idx_dip_ii + 2);

    auto tr_x_xxxxxx_xxxxyy = pbuffer.data(idx_dip_ii + 3);

    auto tr_x_xxxxxx_xxxxyz = pbuffer.data(idx_dip_ii + 4);

    auto tr_x_xxxxxx_xxxxzz = pbuffer.data(idx_dip_ii + 5);

    auto tr_x_xxxxxx_xxxyyy = pbuffer.data(idx_dip_ii + 6);

    auto tr_x_xxxxxx_xxxyyz = pbuffer.data(idx_dip_ii + 7);

    auto tr_x_xxxxxx_xxxyzz = pbuffer.data(idx_dip_ii + 8);

    auto tr_x_xxxxxx_xxxzzz = pbuffer.data(idx_dip_ii + 9);

    auto tr_x_xxxxxx_xxyyyy = pbuffer.data(idx_dip_ii + 10);

    auto tr_x_xxxxxx_xxyyyz = pbuffer.data(idx_dip_ii + 11);

    auto tr_x_xxxxxx_xxyyzz = pbuffer.data(idx_dip_ii + 12);

    auto tr_x_xxxxxx_xxyzzz = pbuffer.data(idx_dip_ii + 13);

    auto tr_x_xxxxxx_xxzzzz = pbuffer.data(idx_dip_ii + 14);

    auto tr_x_xxxxxx_xyyyyy = pbuffer.data(idx_dip_ii + 15);

    auto tr_x_xxxxxx_xyyyyz = pbuffer.data(idx_dip_ii + 16);

    auto tr_x_xxxxxx_xyyyzz = pbuffer.data(idx_dip_ii + 17);

    auto tr_x_xxxxxx_xyyzzz = pbuffer.data(idx_dip_ii + 18);

    auto tr_x_xxxxxx_xyzzzz = pbuffer.data(idx_dip_ii + 19);

    auto tr_x_xxxxxx_xzzzzz = pbuffer.data(idx_dip_ii + 20);

    auto tr_x_xxxxxx_yyyyyy = pbuffer.data(idx_dip_ii + 21);

    auto tr_x_xxxxxx_yyyyyz = pbuffer.data(idx_dip_ii + 22);

    auto tr_x_xxxxxx_yyyyzz = pbuffer.data(idx_dip_ii + 23);

    auto tr_x_xxxxxx_yyyzzz = pbuffer.data(idx_dip_ii + 24);

    auto tr_x_xxxxxx_yyzzzz = pbuffer.data(idx_dip_ii + 25);

    auto tr_x_xxxxxx_yzzzzz = pbuffer.data(idx_dip_ii + 26);

    auto tr_x_xxxxxx_zzzzzz = pbuffer.data(idx_dip_ii + 27);

    #pragma omp simd aligned(pa_x, tr_x_xxxx_xxxxxx, tr_x_xxxx_xxxxxy, tr_x_xxxx_xxxxxz, tr_x_xxxx_xxxxyy, tr_x_xxxx_xxxxyz, tr_x_xxxx_xxxxzz, tr_x_xxxx_xxxyyy, tr_x_xxxx_xxxyyz, tr_x_xxxx_xxxyzz, tr_x_xxxx_xxxzzz, tr_x_xxxx_xxyyyy, tr_x_xxxx_xxyyyz, tr_x_xxxx_xxyyzz, tr_x_xxxx_xxyzzz, tr_x_xxxx_xxzzzz, tr_x_xxxx_xyyyyy, tr_x_xxxx_xyyyyz, tr_x_xxxx_xyyyzz, tr_x_xxxx_xyyzzz, tr_x_xxxx_xyzzzz, tr_x_xxxx_xzzzzz, tr_x_xxxx_yyyyyy, tr_x_xxxx_yyyyyz, tr_x_xxxx_yyyyzz, tr_x_xxxx_yyyzzz, tr_x_xxxx_yyzzzz, tr_x_xxxx_yzzzzz, tr_x_xxxx_zzzzzz, tr_x_xxxxx_xxxxx, tr_x_xxxxx_xxxxxx, tr_x_xxxxx_xxxxxy, tr_x_xxxxx_xxxxxz, tr_x_xxxxx_xxxxy, tr_x_xxxxx_xxxxyy, tr_x_xxxxx_xxxxyz, tr_x_xxxxx_xxxxz, tr_x_xxxxx_xxxxzz, tr_x_xxxxx_xxxyy, tr_x_xxxxx_xxxyyy, tr_x_xxxxx_xxxyyz, tr_x_xxxxx_xxxyz, tr_x_xxxxx_xxxyzz, tr_x_xxxxx_xxxzz, tr_x_xxxxx_xxxzzz, tr_x_xxxxx_xxyyy, tr_x_xxxxx_xxyyyy, tr_x_xxxxx_xxyyyz, tr_x_xxxxx_xxyyz, tr_x_xxxxx_xxyyzz, tr_x_xxxxx_xxyzz, tr_x_xxxxx_xxyzzz, tr_x_xxxxx_xxzzz, tr_x_xxxxx_xxzzzz, tr_x_xxxxx_xyyyy, tr_x_xxxxx_xyyyyy, tr_x_xxxxx_xyyyyz, tr_x_xxxxx_xyyyz, tr_x_xxxxx_xyyyzz, tr_x_xxxxx_xyyzz, tr_x_xxxxx_xyyzzz, tr_x_xxxxx_xyzzz, tr_x_xxxxx_xyzzzz, tr_x_xxxxx_xzzzz, tr_x_xxxxx_xzzzzz, tr_x_xxxxx_yyyyy, tr_x_xxxxx_yyyyyy, tr_x_xxxxx_yyyyyz, tr_x_xxxxx_yyyyz, tr_x_xxxxx_yyyyzz, tr_x_xxxxx_yyyzz, tr_x_xxxxx_yyyzzz, tr_x_xxxxx_yyzzz, tr_x_xxxxx_yyzzzz, tr_x_xxxxx_yzzzz, tr_x_xxxxx_yzzzzz, tr_x_xxxxx_zzzzz, tr_x_xxxxx_zzzzzz, tr_x_xxxxxx_xxxxxx, tr_x_xxxxxx_xxxxxy, tr_x_xxxxxx_xxxxxz, tr_x_xxxxxx_xxxxyy, tr_x_xxxxxx_xxxxyz, tr_x_xxxxxx_xxxxzz, tr_x_xxxxxx_xxxyyy, tr_x_xxxxxx_xxxyyz, tr_x_xxxxxx_xxxyzz, tr_x_xxxxxx_xxxzzz, tr_x_xxxxxx_xxyyyy, tr_x_xxxxxx_xxyyyz, tr_x_xxxxxx_xxyyzz, tr_x_xxxxxx_xxyzzz, tr_x_xxxxxx_xxzzzz, tr_x_xxxxxx_xyyyyy, tr_x_xxxxxx_xyyyyz, tr_x_xxxxxx_xyyyzz, tr_x_xxxxxx_xyyzzz, tr_x_xxxxxx_xyzzzz, tr_x_xxxxxx_xzzzzz, tr_x_xxxxxx_yyyyyy, tr_x_xxxxxx_yyyyyz, tr_x_xxxxxx_yyyyzz, tr_x_xxxxxx_yyyzzz, tr_x_xxxxxx_yyzzzz, tr_x_xxxxxx_yzzzzz, tr_x_xxxxxx_zzzzzz, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxy, ts_xxxxx_xxxxxz, ts_xxxxx_xxxxyy, ts_xxxxx_xxxxyz, ts_xxxxx_xxxxzz, ts_xxxxx_xxxyyy, ts_xxxxx_xxxyyz, ts_xxxxx_xxxyzz, ts_xxxxx_xxxzzz, ts_xxxxx_xxyyyy, ts_xxxxx_xxyyyz, ts_xxxxx_xxyyzz, ts_xxxxx_xxyzzz, ts_xxxxx_xxzzzz, ts_xxxxx_xyyyyy, ts_xxxxx_xyyyyz, ts_xxxxx_xyyyzz, ts_xxxxx_xyyzzz, ts_xxxxx_xyzzzz, ts_xxxxx_xzzzzz, ts_xxxxx_yyyyyy, ts_xxxxx_yyyyyz, ts_xxxxx_yyyyzz, ts_xxxxx_yyyzzz, ts_xxxxx_yyzzzz, ts_xxxxx_yzzzzz, ts_xxxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_xxxxxx[i] = 5.0 * tr_x_xxxx_xxxxxx[i] * fe_0 + 6.0 * tr_x_xxxxx_xxxxx[i] * fe_0 + ts_xxxxx_xxxxxx[i] * fe_0 + tr_x_xxxxx_xxxxxx[i] * pa_x[i];

        tr_x_xxxxxx_xxxxxy[i] = 5.0 * tr_x_xxxx_xxxxxy[i] * fe_0 + 5.0 * tr_x_xxxxx_xxxxy[i] * fe_0 + ts_xxxxx_xxxxxy[i] * fe_0 + tr_x_xxxxx_xxxxxy[i] * pa_x[i];

        tr_x_xxxxxx_xxxxxz[i] = 5.0 * tr_x_xxxx_xxxxxz[i] * fe_0 + 5.0 * tr_x_xxxxx_xxxxz[i] * fe_0 + ts_xxxxx_xxxxxz[i] * fe_0 + tr_x_xxxxx_xxxxxz[i] * pa_x[i];

        tr_x_xxxxxx_xxxxyy[i] = 5.0 * tr_x_xxxx_xxxxyy[i] * fe_0 + 4.0 * tr_x_xxxxx_xxxyy[i] * fe_0 + ts_xxxxx_xxxxyy[i] * fe_0 + tr_x_xxxxx_xxxxyy[i] * pa_x[i];

        tr_x_xxxxxx_xxxxyz[i] = 5.0 * tr_x_xxxx_xxxxyz[i] * fe_0 + 4.0 * tr_x_xxxxx_xxxyz[i] * fe_0 + ts_xxxxx_xxxxyz[i] * fe_0 + tr_x_xxxxx_xxxxyz[i] * pa_x[i];

        tr_x_xxxxxx_xxxxzz[i] = 5.0 * tr_x_xxxx_xxxxzz[i] * fe_0 + 4.0 * tr_x_xxxxx_xxxzz[i] * fe_0 + ts_xxxxx_xxxxzz[i] * fe_0 + tr_x_xxxxx_xxxxzz[i] * pa_x[i];

        tr_x_xxxxxx_xxxyyy[i] = 5.0 * tr_x_xxxx_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxxxx_xxyyy[i] * fe_0 + ts_xxxxx_xxxyyy[i] * fe_0 + tr_x_xxxxx_xxxyyy[i] * pa_x[i];

        tr_x_xxxxxx_xxxyyz[i] = 5.0 * tr_x_xxxx_xxxyyz[i] * fe_0 + 3.0 * tr_x_xxxxx_xxyyz[i] * fe_0 + ts_xxxxx_xxxyyz[i] * fe_0 + tr_x_xxxxx_xxxyyz[i] * pa_x[i];

        tr_x_xxxxxx_xxxyzz[i] = 5.0 * tr_x_xxxx_xxxyzz[i] * fe_0 + 3.0 * tr_x_xxxxx_xxyzz[i] * fe_0 + ts_xxxxx_xxxyzz[i] * fe_0 + tr_x_xxxxx_xxxyzz[i] * pa_x[i];

        tr_x_xxxxxx_xxxzzz[i] = 5.0 * tr_x_xxxx_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxxxx_xxzzz[i] * fe_0 + ts_xxxxx_xxxzzz[i] * fe_0 + tr_x_xxxxx_xxxzzz[i] * pa_x[i];

        tr_x_xxxxxx_xxyyyy[i] = 5.0 * tr_x_xxxx_xxyyyy[i] * fe_0 + 2.0 * tr_x_xxxxx_xyyyy[i] * fe_0 + ts_xxxxx_xxyyyy[i] * fe_0 + tr_x_xxxxx_xxyyyy[i] * pa_x[i];

        tr_x_xxxxxx_xxyyyz[i] = 5.0 * tr_x_xxxx_xxyyyz[i] * fe_0 + 2.0 * tr_x_xxxxx_xyyyz[i] * fe_0 + ts_xxxxx_xxyyyz[i] * fe_0 + tr_x_xxxxx_xxyyyz[i] * pa_x[i];

        tr_x_xxxxxx_xxyyzz[i] = 5.0 * tr_x_xxxx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxxx_xyyzz[i] * fe_0 + ts_xxxxx_xxyyzz[i] * fe_0 + tr_x_xxxxx_xxyyzz[i] * pa_x[i];

        tr_x_xxxxxx_xxyzzz[i] = 5.0 * tr_x_xxxx_xxyzzz[i] * fe_0 + 2.0 * tr_x_xxxxx_xyzzz[i] * fe_0 + ts_xxxxx_xxyzzz[i] * fe_0 + tr_x_xxxxx_xxyzzz[i] * pa_x[i];

        tr_x_xxxxxx_xxzzzz[i] = 5.0 * tr_x_xxxx_xxzzzz[i] * fe_0 + 2.0 * tr_x_xxxxx_xzzzz[i] * fe_0 + ts_xxxxx_xxzzzz[i] * fe_0 + tr_x_xxxxx_xxzzzz[i] * pa_x[i];

        tr_x_xxxxxx_xyyyyy[i] = 5.0 * tr_x_xxxx_xyyyyy[i] * fe_0 + tr_x_xxxxx_yyyyy[i] * fe_0 + ts_xxxxx_xyyyyy[i] * fe_0 + tr_x_xxxxx_xyyyyy[i] * pa_x[i];

        tr_x_xxxxxx_xyyyyz[i] = 5.0 * tr_x_xxxx_xyyyyz[i] * fe_0 + tr_x_xxxxx_yyyyz[i] * fe_0 + ts_xxxxx_xyyyyz[i] * fe_0 + tr_x_xxxxx_xyyyyz[i] * pa_x[i];

        tr_x_xxxxxx_xyyyzz[i] = 5.0 * tr_x_xxxx_xyyyzz[i] * fe_0 + tr_x_xxxxx_yyyzz[i] * fe_0 + ts_xxxxx_xyyyzz[i] * fe_0 + tr_x_xxxxx_xyyyzz[i] * pa_x[i];

        tr_x_xxxxxx_xyyzzz[i] = 5.0 * tr_x_xxxx_xyyzzz[i] * fe_0 + tr_x_xxxxx_yyzzz[i] * fe_0 + ts_xxxxx_xyyzzz[i] * fe_0 + tr_x_xxxxx_xyyzzz[i] * pa_x[i];

        tr_x_xxxxxx_xyzzzz[i] = 5.0 * tr_x_xxxx_xyzzzz[i] * fe_0 + tr_x_xxxxx_yzzzz[i] * fe_0 + ts_xxxxx_xyzzzz[i] * fe_0 + tr_x_xxxxx_xyzzzz[i] * pa_x[i];

        tr_x_xxxxxx_xzzzzz[i] = 5.0 * tr_x_xxxx_xzzzzz[i] * fe_0 + tr_x_xxxxx_zzzzz[i] * fe_0 + ts_xxxxx_xzzzzz[i] * fe_0 + tr_x_xxxxx_xzzzzz[i] * pa_x[i];

        tr_x_xxxxxx_yyyyyy[i] = 5.0 * tr_x_xxxx_yyyyyy[i] * fe_0 + ts_xxxxx_yyyyyy[i] * fe_0 + tr_x_xxxxx_yyyyyy[i] * pa_x[i];

        tr_x_xxxxxx_yyyyyz[i] = 5.0 * tr_x_xxxx_yyyyyz[i] * fe_0 + ts_xxxxx_yyyyyz[i] * fe_0 + tr_x_xxxxx_yyyyyz[i] * pa_x[i];

        tr_x_xxxxxx_yyyyzz[i] = 5.0 * tr_x_xxxx_yyyyzz[i] * fe_0 + ts_xxxxx_yyyyzz[i] * fe_0 + tr_x_xxxxx_yyyyzz[i] * pa_x[i];

        tr_x_xxxxxx_yyyzzz[i] = 5.0 * tr_x_xxxx_yyyzzz[i] * fe_0 + ts_xxxxx_yyyzzz[i] * fe_0 + tr_x_xxxxx_yyyzzz[i] * pa_x[i];

        tr_x_xxxxxx_yyzzzz[i] = 5.0 * tr_x_xxxx_yyzzzz[i] * fe_0 + ts_xxxxx_yyzzzz[i] * fe_0 + tr_x_xxxxx_yyzzzz[i] * pa_x[i];

        tr_x_xxxxxx_yzzzzz[i] = 5.0 * tr_x_xxxx_yzzzzz[i] * fe_0 + ts_xxxxx_yzzzzz[i] * fe_0 + tr_x_xxxxx_yzzzzz[i] * pa_x[i];

        tr_x_xxxxxx_zzzzzz[i] = 5.0 * tr_x_xxxx_zzzzzz[i] * fe_0 + ts_xxxxx_zzzzzz[i] * fe_0 + tr_x_xxxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : II

    auto tr_x_xxxxxy_xxxxxx = pbuffer.data(idx_dip_ii + 28);

    auto tr_x_xxxxxy_xxxxxy = pbuffer.data(idx_dip_ii + 29);

    auto tr_x_xxxxxy_xxxxxz = pbuffer.data(idx_dip_ii + 30);

    auto tr_x_xxxxxy_xxxxyy = pbuffer.data(idx_dip_ii + 31);

    auto tr_x_xxxxxy_xxxxyz = pbuffer.data(idx_dip_ii + 32);

    auto tr_x_xxxxxy_xxxxzz = pbuffer.data(idx_dip_ii + 33);

    auto tr_x_xxxxxy_xxxyyy = pbuffer.data(idx_dip_ii + 34);

    auto tr_x_xxxxxy_xxxyyz = pbuffer.data(idx_dip_ii + 35);

    auto tr_x_xxxxxy_xxxyzz = pbuffer.data(idx_dip_ii + 36);

    auto tr_x_xxxxxy_xxxzzz = pbuffer.data(idx_dip_ii + 37);

    auto tr_x_xxxxxy_xxyyyy = pbuffer.data(idx_dip_ii + 38);

    auto tr_x_xxxxxy_xxyyyz = pbuffer.data(idx_dip_ii + 39);

    auto tr_x_xxxxxy_xxyyzz = pbuffer.data(idx_dip_ii + 40);

    auto tr_x_xxxxxy_xxyzzz = pbuffer.data(idx_dip_ii + 41);

    auto tr_x_xxxxxy_xxzzzz = pbuffer.data(idx_dip_ii + 42);

    auto tr_x_xxxxxy_xyyyyy = pbuffer.data(idx_dip_ii + 43);

    auto tr_x_xxxxxy_xyyyyz = pbuffer.data(idx_dip_ii + 44);

    auto tr_x_xxxxxy_xyyyzz = pbuffer.data(idx_dip_ii + 45);

    auto tr_x_xxxxxy_xyyzzz = pbuffer.data(idx_dip_ii + 46);

    auto tr_x_xxxxxy_xyzzzz = pbuffer.data(idx_dip_ii + 47);

    auto tr_x_xxxxxy_xzzzzz = pbuffer.data(idx_dip_ii + 48);

    auto tr_x_xxxxxy_yyyyyy = pbuffer.data(idx_dip_ii + 49);

    auto tr_x_xxxxxy_yyyyyz = pbuffer.data(idx_dip_ii + 50);

    auto tr_x_xxxxxy_yyyyzz = pbuffer.data(idx_dip_ii + 51);

    auto tr_x_xxxxxy_yyyzzz = pbuffer.data(idx_dip_ii + 52);

    auto tr_x_xxxxxy_yyzzzz = pbuffer.data(idx_dip_ii + 53);

    auto tr_x_xxxxxy_yzzzzz = pbuffer.data(idx_dip_ii + 54);

    auto tr_x_xxxxxy_zzzzzz = pbuffer.data(idx_dip_ii + 55);

    #pragma omp simd aligned(pa_y, tr_x_xxxxx_xxxxx, tr_x_xxxxx_xxxxxx, tr_x_xxxxx_xxxxxy, tr_x_xxxxx_xxxxxz, tr_x_xxxxx_xxxxy, tr_x_xxxxx_xxxxyy, tr_x_xxxxx_xxxxyz, tr_x_xxxxx_xxxxz, tr_x_xxxxx_xxxxzz, tr_x_xxxxx_xxxyy, tr_x_xxxxx_xxxyyy, tr_x_xxxxx_xxxyyz, tr_x_xxxxx_xxxyz, tr_x_xxxxx_xxxyzz, tr_x_xxxxx_xxxzz, tr_x_xxxxx_xxxzzz, tr_x_xxxxx_xxyyy, tr_x_xxxxx_xxyyyy, tr_x_xxxxx_xxyyyz, tr_x_xxxxx_xxyyz, tr_x_xxxxx_xxyyzz, tr_x_xxxxx_xxyzz, tr_x_xxxxx_xxyzzz, tr_x_xxxxx_xxzzz, tr_x_xxxxx_xxzzzz, tr_x_xxxxx_xyyyy, tr_x_xxxxx_xyyyyy, tr_x_xxxxx_xyyyyz, tr_x_xxxxx_xyyyz, tr_x_xxxxx_xyyyzz, tr_x_xxxxx_xyyzz, tr_x_xxxxx_xyyzzz, tr_x_xxxxx_xyzzz, tr_x_xxxxx_xyzzzz, tr_x_xxxxx_xzzzz, tr_x_xxxxx_xzzzzz, tr_x_xxxxx_yyyyy, tr_x_xxxxx_yyyyyy, tr_x_xxxxx_yyyyyz, tr_x_xxxxx_yyyyz, tr_x_xxxxx_yyyyzz, tr_x_xxxxx_yyyzz, tr_x_xxxxx_yyyzzz, tr_x_xxxxx_yyzzz, tr_x_xxxxx_yyzzzz, tr_x_xxxxx_yzzzz, tr_x_xxxxx_yzzzzz, tr_x_xxxxx_zzzzz, tr_x_xxxxx_zzzzzz, tr_x_xxxxxy_xxxxxx, tr_x_xxxxxy_xxxxxy, tr_x_xxxxxy_xxxxxz, tr_x_xxxxxy_xxxxyy, tr_x_xxxxxy_xxxxyz, tr_x_xxxxxy_xxxxzz, tr_x_xxxxxy_xxxyyy, tr_x_xxxxxy_xxxyyz, tr_x_xxxxxy_xxxyzz, tr_x_xxxxxy_xxxzzz, tr_x_xxxxxy_xxyyyy, tr_x_xxxxxy_xxyyyz, tr_x_xxxxxy_xxyyzz, tr_x_xxxxxy_xxyzzz, tr_x_xxxxxy_xxzzzz, tr_x_xxxxxy_xyyyyy, tr_x_xxxxxy_xyyyyz, tr_x_xxxxxy_xyyyzz, tr_x_xxxxxy_xyyzzz, tr_x_xxxxxy_xyzzzz, tr_x_xxxxxy_xzzzzz, tr_x_xxxxxy_yyyyyy, tr_x_xxxxxy_yyyyyz, tr_x_xxxxxy_yyyyzz, tr_x_xxxxxy_yyyzzz, tr_x_xxxxxy_yyzzzz, tr_x_xxxxxy_yzzzzz, tr_x_xxxxxy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxy_xxxxxx[i] = tr_x_xxxxx_xxxxxx[i] * pa_y[i];

        tr_x_xxxxxy_xxxxxy[i] = tr_x_xxxxx_xxxxx[i] * fe_0 + tr_x_xxxxx_xxxxxy[i] * pa_y[i];

        tr_x_xxxxxy_xxxxxz[i] = tr_x_xxxxx_xxxxxz[i] * pa_y[i];

        tr_x_xxxxxy_xxxxyy[i] = 2.0 * tr_x_xxxxx_xxxxy[i] * fe_0 + tr_x_xxxxx_xxxxyy[i] * pa_y[i];

        tr_x_xxxxxy_xxxxyz[i] = tr_x_xxxxx_xxxxz[i] * fe_0 + tr_x_xxxxx_xxxxyz[i] * pa_y[i];

        tr_x_xxxxxy_xxxxzz[i] = tr_x_xxxxx_xxxxzz[i] * pa_y[i];

        tr_x_xxxxxy_xxxyyy[i] = 3.0 * tr_x_xxxxx_xxxyy[i] * fe_0 + tr_x_xxxxx_xxxyyy[i] * pa_y[i];

        tr_x_xxxxxy_xxxyyz[i] = 2.0 * tr_x_xxxxx_xxxyz[i] * fe_0 + tr_x_xxxxx_xxxyyz[i] * pa_y[i];

        tr_x_xxxxxy_xxxyzz[i] = tr_x_xxxxx_xxxzz[i] * fe_0 + tr_x_xxxxx_xxxyzz[i] * pa_y[i];

        tr_x_xxxxxy_xxxzzz[i] = tr_x_xxxxx_xxxzzz[i] * pa_y[i];

        tr_x_xxxxxy_xxyyyy[i] = 4.0 * tr_x_xxxxx_xxyyy[i] * fe_0 + tr_x_xxxxx_xxyyyy[i] * pa_y[i];

        tr_x_xxxxxy_xxyyyz[i] = 3.0 * tr_x_xxxxx_xxyyz[i] * fe_0 + tr_x_xxxxx_xxyyyz[i] * pa_y[i];

        tr_x_xxxxxy_xxyyzz[i] = 2.0 * tr_x_xxxxx_xxyzz[i] * fe_0 + tr_x_xxxxx_xxyyzz[i] * pa_y[i];

        tr_x_xxxxxy_xxyzzz[i] = tr_x_xxxxx_xxzzz[i] * fe_0 + tr_x_xxxxx_xxyzzz[i] * pa_y[i];

        tr_x_xxxxxy_xxzzzz[i] = tr_x_xxxxx_xxzzzz[i] * pa_y[i];

        tr_x_xxxxxy_xyyyyy[i] = 5.0 * tr_x_xxxxx_xyyyy[i] * fe_0 + tr_x_xxxxx_xyyyyy[i] * pa_y[i];

        tr_x_xxxxxy_xyyyyz[i] = 4.0 * tr_x_xxxxx_xyyyz[i] * fe_0 + tr_x_xxxxx_xyyyyz[i] * pa_y[i];

        tr_x_xxxxxy_xyyyzz[i] = 3.0 * tr_x_xxxxx_xyyzz[i] * fe_0 + tr_x_xxxxx_xyyyzz[i] * pa_y[i];

        tr_x_xxxxxy_xyyzzz[i] = 2.0 * tr_x_xxxxx_xyzzz[i] * fe_0 + tr_x_xxxxx_xyyzzz[i] * pa_y[i];

        tr_x_xxxxxy_xyzzzz[i] = tr_x_xxxxx_xzzzz[i] * fe_0 + tr_x_xxxxx_xyzzzz[i] * pa_y[i];

        tr_x_xxxxxy_xzzzzz[i] = tr_x_xxxxx_xzzzzz[i] * pa_y[i];

        tr_x_xxxxxy_yyyyyy[i] = 6.0 * tr_x_xxxxx_yyyyy[i] * fe_0 + tr_x_xxxxx_yyyyyy[i] * pa_y[i];

        tr_x_xxxxxy_yyyyyz[i] = 5.0 * tr_x_xxxxx_yyyyz[i] * fe_0 + tr_x_xxxxx_yyyyyz[i] * pa_y[i];

        tr_x_xxxxxy_yyyyzz[i] = 4.0 * tr_x_xxxxx_yyyzz[i] * fe_0 + tr_x_xxxxx_yyyyzz[i] * pa_y[i];

        tr_x_xxxxxy_yyyzzz[i] = 3.0 * tr_x_xxxxx_yyzzz[i] * fe_0 + tr_x_xxxxx_yyyzzz[i] * pa_y[i];

        tr_x_xxxxxy_yyzzzz[i] = 2.0 * tr_x_xxxxx_yzzzz[i] * fe_0 + tr_x_xxxxx_yyzzzz[i] * pa_y[i];

        tr_x_xxxxxy_yzzzzz[i] = tr_x_xxxxx_zzzzz[i] * fe_0 + tr_x_xxxxx_yzzzzz[i] * pa_y[i];

        tr_x_xxxxxy_zzzzzz[i] = tr_x_xxxxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : II

    auto tr_x_xxxxxz_xxxxxx = pbuffer.data(idx_dip_ii + 56);

    auto tr_x_xxxxxz_xxxxxy = pbuffer.data(idx_dip_ii + 57);

    auto tr_x_xxxxxz_xxxxxz = pbuffer.data(idx_dip_ii + 58);

    auto tr_x_xxxxxz_xxxxyy = pbuffer.data(idx_dip_ii + 59);

    auto tr_x_xxxxxz_xxxxyz = pbuffer.data(idx_dip_ii + 60);

    auto tr_x_xxxxxz_xxxxzz = pbuffer.data(idx_dip_ii + 61);

    auto tr_x_xxxxxz_xxxyyy = pbuffer.data(idx_dip_ii + 62);

    auto tr_x_xxxxxz_xxxyyz = pbuffer.data(idx_dip_ii + 63);

    auto tr_x_xxxxxz_xxxyzz = pbuffer.data(idx_dip_ii + 64);

    auto tr_x_xxxxxz_xxxzzz = pbuffer.data(idx_dip_ii + 65);

    auto tr_x_xxxxxz_xxyyyy = pbuffer.data(idx_dip_ii + 66);

    auto tr_x_xxxxxz_xxyyyz = pbuffer.data(idx_dip_ii + 67);

    auto tr_x_xxxxxz_xxyyzz = pbuffer.data(idx_dip_ii + 68);

    auto tr_x_xxxxxz_xxyzzz = pbuffer.data(idx_dip_ii + 69);

    auto tr_x_xxxxxz_xxzzzz = pbuffer.data(idx_dip_ii + 70);

    auto tr_x_xxxxxz_xyyyyy = pbuffer.data(idx_dip_ii + 71);

    auto tr_x_xxxxxz_xyyyyz = pbuffer.data(idx_dip_ii + 72);

    auto tr_x_xxxxxz_xyyyzz = pbuffer.data(idx_dip_ii + 73);

    auto tr_x_xxxxxz_xyyzzz = pbuffer.data(idx_dip_ii + 74);

    auto tr_x_xxxxxz_xyzzzz = pbuffer.data(idx_dip_ii + 75);

    auto tr_x_xxxxxz_xzzzzz = pbuffer.data(idx_dip_ii + 76);

    auto tr_x_xxxxxz_yyyyyy = pbuffer.data(idx_dip_ii + 77);

    auto tr_x_xxxxxz_yyyyyz = pbuffer.data(idx_dip_ii + 78);

    auto tr_x_xxxxxz_yyyyzz = pbuffer.data(idx_dip_ii + 79);

    auto tr_x_xxxxxz_yyyzzz = pbuffer.data(idx_dip_ii + 80);

    auto tr_x_xxxxxz_yyzzzz = pbuffer.data(idx_dip_ii + 81);

    auto tr_x_xxxxxz_yzzzzz = pbuffer.data(idx_dip_ii + 82);

    auto tr_x_xxxxxz_zzzzzz = pbuffer.data(idx_dip_ii + 83);

    #pragma omp simd aligned(pa_z, tr_x_xxxxx_xxxxx, tr_x_xxxxx_xxxxxx, tr_x_xxxxx_xxxxxy, tr_x_xxxxx_xxxxxz, tr_x_xxxxx_xxxxy, tr_x_xxxxx_xxxxyy, tr_x_xxxxx_xxxxyz, tr_x_xxxxx_xxxxz, tr_x_xxxxx_xxxxzz, tr_x_xxxxx_xxxyy, tr_x_xxxxx_xxxyyy, tr_x_xxxxx_xxxyyz, tr_x_xxxxx_xxxyz, tr_x_xxxxx_xxxyzz, tr_x_xxxxx_xxxzz, tr_x_xxxxx_xxxzzz, tr_x_xxxxx_xxyyy, tr_x_xxxxx_xxyyyy, tr_x_xxxxx_xxyyyz, tr_x_xxxxx_xxyyz, tr_x_xxxxx_xxyyzz, tr_x_xxxxx_xxyzz, tr_x_xxxxx_xxyzzz, tr_x_xxxxx_xxzzz, tr_x_xxxxx_xxzzzz, tr_x_xxxxx_xyyyy, tr_x_xxxxx_xyyyyy, tr_x_xxxxx_xyyyyz, tr_x_xxxxx_xyyyz, tr_x_xxxxx_xyyyzz, tr_x_xxxxx_xyyzz, tr_x_xxxxx_xyyzzz, tr_x_xxxxx_xyzzz, tr_x_xxxxx_xyzzzz, tr_x_xxxxx_xzzzz, tr_x_xxxxx_xzzzzz, tr_x_xxxxx_yyyyy, tr_x_xxxxx_yyyyyy, tr_x_xxxxx_yyyyyz, tr_x_xxxxx_yyyyz, tr_x_xxxxx_yyyyzz, tr_x_xxxxx_yyyzz, tr_x_xxxxx_yyyzzz, tr_x_xxxxx_yyzzz, tr_x_xxxxx_yyzzzz, tr_x_xxxxx_yzzzz, tr_x_xxxxx_yzzzzz, tr_x_xxxxx_zzzzz, tr_x_xxxxx_zzzzzz, tr_x_xxxxxz_xxxxxx, tr_x_xxxxxz_xxxxxy, tr_x_xxxxxz_xxxxxz, tr_x_xxxxxz_xxxxyy, tr_x_xxxxxz_xxxxyz, tr_x_xxxxxz_xxxxzz, tr_x_xxxxxz_xxxyyy, tr_x_xxxxxz_xxxyyz, tr_x_xxxxxz_xxxyzz, tr_x_xxxxxz_xxxzzz, tr_x_xxxxxz_xxyyyy, tr_x_xxxxxz_xxyyyz, tr_x_xxxxxz_xxyyzz, tr_x_xxxxxz_xxyzzz, tr_x_xxxxxz_xxzzzz, tr_x_xxxxxz_xyyyyy, tr_x_xxxxxz_xyyyyz, tr_x_xxxxxz_xyyyzz, tr_x_xxxxxz_xyyzzz, tr_x_xxxxxz_xyzzzz, tr_x_xxxxxz_xzzzzz, tr_x_xxxxxz_yyyyyy, tr_x_xxxxxz_yyyyyz, tr_x_xxxxxz_yyyyzz, tr_x_xxxxxz_yyyzzz, tr_x_xxxxxz_yyzzzz, tr_x_xxxxxz_yzzzzz, tr_x_xxxxxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxz_xxxxxx[i] = tr_x_xxxxx_xxxxxx[i] * pa_z[i];

        tr_x_xxxxxz_xxxxxy[i] = tr_x_xxxxx_xxxxxy[i] * pa_z[i];

        tr_x_xxxxxz_xxxxxz[i] = tr_x_xxxxx_xxxxx[i] * fe_0 + tr_x_xxxxx_xxxxxz[i] * pa_z[i];

        tr_x_xxxxxz_xxxxyy[i] = tr_x_xxxxx_xxxxyy[i] * pa_z[i];

        tr_x_xxxxxz_xxxxyz[i] = tr_x_xxxxx_xxxxy[i] * fe_0 + tr_x_xxxxx_xxxxyz[i] * pa_z[i];

        tr_x_xxxxxz_xxxxzz[i] = 2.0 * tr_x_xxxxx_xxxxz[i] * fe_0 + tr_x_xxxxx_xxxxzz[i] * pa_z[i];

        tr_x_xxxxxz_xxxyyy[i] = tr_x_xxxxx_xxxyyy[i] * pa_z[i];

        tr_x_xxxxxz_xxxyyz[i] = tr_x_xxxxx_xxxyy[i] * fe_0 + tr_x_xxxxx_xxxyyz[i] * pa_z[i];

        tr_x_xxxxxz_xxxyzz[i] = 2.0 * tr_x_xxxxx_xxxyz[i] * fe_0 + tr_x_xxxxx_xxxyzz[i] * pa_z[i];

        tr_x_xxxxxz_xxxzzz[i] = 3.0 * tr_x_xxxxx_xxxzz[i] * fe_0 + tr_x_xxxxx_xxxzzz[i] * pa_z[i];

        tr_x_xxxxxz_xxyyyy[i] = tr_x_xxxxx_xxyyyy[i] * pa_z[i];

        tr_x_xxxxxz_xxyyyz[i] = tr_x_xxxxx_xxyyy[i] * fe_0 + tr_x_xxxxx_xxyyyz[i] * pa_z[i];

        tr_x_xxxxxz_xxyyzz[i] = 2.0 * tr_x_xxxxx_xxyyz[i] * fe_0 + tr_x_xxxxx_xxyyzz[i] * pa_z[i];

        tr_x_xxxxxz_xxyzzz[i] = 3.0 * tr_x_xxxxx_xxyzz[i] * fe_0 + tr_x_xxxxx_xxyzzz[i] * pa_z[i];

        tr_x_xxxxxz_xxzzzz[i] = 4.0 * tr_x_xxxxx_xxzzz[i] * fe_0 + tr_x_xxxxx_xxzzzz[i] * pa_z[i];

        tr_x_xxxxxz_xyyyyy[i] = tr_x_xxxxx_xyyyyy[i] * pa_z[i];

        tr_x_xxxxxz_xyyyyz[i] = tr_x_xxxxx_xyyyy[i] * fe_0 + tr_x_xxxxx_xyyyyz[i] * pa_z[i];

        tr_x_xxxxxz_xyyyzz[i] = 2.0 * tr_x_xxxxx_xyyyz[i] * fe_0 + tr_x_xxxxx_xyyyzz[i] * pa_z[i];

        tr_x_xxxxxz_xyyzzz[i] = 3.0 * tr_x_xxxxx_xyyzz[i] * fe_0 + tr_x_xxxxx_xyyzzz[i] * pa_z[i];

        tr_x_xxxxxz_xyzzzz[i] = 4.0 * tr_x_xxxxx_xyzzz[i] * fe_0 + tr_x_xxxxx_xyzzzz[i] * pa_z[i];

        tr_x_xxxxxz_xzzzzz[i] = 5.0 * tr_x_xxxxx_xzzzz[i] * fe_0 + tr_x_xxxxx_xzzzzz[i] * pa_z[i];

        tr_x_xxxxxz_yyyyyy[i] = tr_x_xxxxx_yyyyyy[i] * pa_z[i];

        tr_x_xxxxxz_yyyyyz[i] = tr_x_xxxxx_yyyyy[i] * fe_0 + tr_x_xxxxx_yyyyyz[i] * pa_z[i];

        tr_x_xxxxxz_yyyyzz[i] = 2.0 * tr_x_xxxxx_yyyyz[i] * fe_0 + tr_x_xxxxx_yyyyzz[i] * pa_z[i];

        tr_x_xxxxxz_yyyzzz[i] = 3.0 * tr_x_xxxxx_yyyzz[i] * fe_0 + tr_x_xxxxx_yyyzzz[i] * pa_z[i];

        tr_x_xxxxxz_yyzzzz[i] = 4.0 * tr_x_xxxxx_yyzzz[i] * fe_0 + tr_x_xxxxx_yyzzzz[i] * pa_z[i];

        tr_x_xxxxxz_yzzzzz[i] = 5.0 * tr_x_xxxxx_yzzzz[i] * fe_0 + tr_x_xxxxx_yzzzzz[i] * pa_z[i];

        tr_x_xxxxxz_zzzzzz[i] = 6.0 * tr_x_xxxxx_zzzzz[i] * fe_0 + tr_x_xxxxx_zzzzzz[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : II

    auto tr_x_xxxxyy_xxxxxx = pbuffer.data(idx_dip_ii + 84);

    auto tr_x_xxxxyy_xxxxxy = pbuffer.data(idx_dip_ii + 85);

    auto tr_x_xxxxyy_xxxxxz = pbuffer.data(idx_dip_ii + 86);

    auto tr_x_xxxxyy_xxxxyy = pbuffer.data(idx_dip_ii + 87);

    auto tr_x_xxxxyy_xxxxyz = pbuffer.data(idx_dip_ii + 88);

    auto tr_x_xxxxyy_xxxxzz = pbuffer.data(idx_dip_ii + 89);

    auto tr_x_xxxxyy_xxxyyy = pbuffer.data(idx_dip_ii + 90);

    auto tr_x_xxxxyy_xxxyyz = pbuffer.data(idx_dip_ii + 91);

    auto tr_x_xxxxyy_xxxyzz = pbuffer.data(idx_dip_ii + 92);

    auto tr_x_xxxxyy_xxxzzz = pbuffer.data(idx_dip_ii + 93);

    auto tr_x_xxxxyy_xxyyyy = pbuffer.data(idx_dip_ii + 94);

    auto tr_x_xxxxyy_xxyyyz = pbuffer.data(idx_dip_ii + 95);

    auto tr_x_xxxxyy_xxyyzz = pbuffer.data(idx_dip_ii + 96);

    auto tr_x_xxxxyy_xxyzzz = pbuffer.data(idx_dip_ii + 97);

    auto tr_x_xxxxyy_xxzzzz = pbuffer.data(idx_dip_ii + 98);

    auto tr_x_xxxxyy_xyyyyy = pbuffer.data(idx_dip_ii + 99);

    auto tr_x_xxxxyy_xyyyyz = pbuffer.data(idx_dip_ii + 100);

    auto tr_x_xxxxyy_xyyyzz = pbuffer.data(idx_dip_ii + 101);

    auto tr_x_xxxxyy_xyyzzz = pbuffer.data(idx_dip_ii + 102);

    auto tr_x_xxxxyy_xyzzzz = pbuffer.data(idx_dip_ii + 103);

    auto tr_x_xxxxyy_xzzzzz = pbuffer.data(idx_dip_ii + 104);

    auto tr_x_xxxxyy_yyyyyy = pbuffer.data(idx_dip_ii + 105);

    auto tr_x_xxxxyy_yyyyyz = pbuffer.data(idx_dip_ii + 106);

    auto tr_x_xxxxyy_yyyyzz = pbuffer.data(idx_dip_ii + 107);

    auto tr_x_xxxxyy_yyyzzz = pbuffer.data(idx_dip_ii + 108);

    auto tr_x_xxxxyy_yyzzzz = pbuffer.data(idx_dip_ii + 109);

    auto tr_x_xxxxyy_yzzzzz = pbuffer.data(idx_dip_ii + 110);

    auto tr_x_xxxxyy_zzzzzz = pbuffer.data(idx_dip_ii + 111);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxxx_xxxxxx, tr_x_xxxx_xxxxxy, tr_x_xxxx_xxxxxz, tr_x_xxxx_xxxxyy, tr_x_xxxx_xxxxyz, tr_x_xxxx_xxxxzz, tr_x_xxxx_xxxyyy, tr_x_xxxx_xxxyyz, tr_x_xxxx_xxxyzz, tr_x_xxxx_xxxzzz, tr_x_xxxx_xxyyyy, tr_x_xxxx_xxyyyz, tr_x_xxxx_xxyyzz, tr_x_xxxx_xxyzzz, tr_x_xxxx_xxzzzz, tr_x_xxxx_xyyyyy, tr_x_xxxx_xyyyyz, tr_x_xxxx_xyyyzz, tr_x_xxxx_xyyzzz, tr_x_xxxx_xyzzzz, tr_x_xxxx_xzzzzz, tr_x_xxxx_zzzzzz, tr_x_xxxxy_xxxxx, tr_x_xxxxy_xxxxxx, tr_x_xxxxy_xxxxxy, tr_x_xxxxy_xxxxxz, tr_x_xxxxy_xxxxy, tr_x_xxxxy_xxxxyy, tr_x_xxxxy_xxxxyz, tr_x_xxxxy_xxxxz, tr_x_xxxxy_xxxxzz, tr_x_xxxxy_xxxyy, tr_x_xxxxy_xxxyyy, tr_x_xxxxy_xxxyyz, tr_x_xxxxy_xxxyz, tr_x_xxxxy_xxxyzz, tr_x_xxxxy_xxxzz, tr_x_xxxxy_xxxzzz, tr_x_xxxxy_xxyyy, tr_x_xxxxy_xxyyyy, tr_x_xxxxy_xxyyyz, tr_x_xxxxy_xxyyz, tr_x_xxxxy_xxyyzz, tr_x_xxxxy_xxyzz, tr_x_xxxxy_xxyzzz, tr_x_xxxxy_xxzzz, tr_x_xxxxy_xxzzzz, tr_x_xxxxy_xyyyy, tr_x_xxxxy_xyyyyy, tr_x_xxxxy_xyyyyz, tr_x_xxxxy_xyyyz, tr_x_xxxxy_xyyyzz, tr_x_xxxxy_xyyzz, tr_x_xxxxy_xyyzzz, tr_x_xxxxy_xyzzz, tr_x_xxxxy_xyzzzz, tr_x_xxxxy_xzzzz, tr_x_xxxxy_xzzzzz, tr_x_xxxxy_zzzzzz, tr_x_xxxxyy_xxxxxx, tr_x_xxxxyy_xxxxxy, tr_x_xxxxyy_xxxxxz, tr_x_xxxxyy_xxxxyy, tr_x_xxxxyy_xxxxyz, tr_x_xxxxyy_xxxxzz, tr_x_xxxxyy_xxxyyy, tr_x_xxxxyy_xxxyyz, tr_x_xxxxyy_xxxyzz, tr_x_xxxxyy_xxxzzz, tr_x_xxxxyy_xxyyyy, tr_x_xxxxyy_xxyyyz, tr_x_xxxxyy_xxyyzz, tr_x_xxxxyy_xxyzzz, tr_x_xxxxyy_xxzzzz, tr_x_xxxxyy_xyyyyy, tr_x_xxxxyy_xyyyyz, tr_x_xxxxyy_xyyyzz, tr_x_xxxxyy_xyyzzz, tr_x_xxxxyy_xyzzzz, tr_x_xxxxyy_xzzzzz, tr_x_xxxxyy_yyyyyy, tr_x_xxxxyy_yyyyyz, tr_x_xxxxyy_yyyyzz, tr_x_xxxxyy_yyyzzz, tr_x_xxxxyy_yyzzzz, tr_x_xxxxyy_yzzzzz, tr_x_xxxxyy_zzzzzz, tr_x_xxxyy_yyyyyy, tr_x_xxxyy_yyyyyz, tr_x_xxxyy_yyyyzz, tr_x_xxxyy_yyyzzz, tr_x_xxxyy_yyzzzz, tr_x_xxxyy_yzzzzz, tr_x_xxyy_yyyyyy, tr_x_xxyy_yyyyyz, tr_x_xxyy_yyyyzz, tr_x_xxyy_yyyzzz, tr_x_xxyy_yyzzzz, tr_x_xxyy_yzzzzz, ts_xxxyy_yyyyyy, ts_xxxyy_yyyyyz, ts_xxxyy_yyyyzz, ts_xxxyy_yyyzzz, ts_xxxyy_yyzzzz, ts_xxxyy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyy_xxxxxx[i] = tr_x_xxxx_xxxxxx[i] * fe_0 + tr_x_xxxxy_xxxxxx[i] * pa_y[i];

        tr_x_xxxxyy_xxxxxy[i] = tr_x_xxxx_xxxxxy[i] * fe_0 + tr_x_xxxxy_xxxxx[i] * fe_0 + tr_x_xxxxy_xxxxxy[i] * pa_y[i];

        tr_x_xxxxyy_xxxxxz[i] = tr_x_xxxx_xxxxxz[i] * fe_0 + tr_x_xxxxy_xxxxxz[i] * pa_y[i];

        tr_x_xxxxyy_xxxxyy[i] = tr_x_xxxx_xxxxyy[i] * fe_0 + 2.0 * tr_x_xxxxy_xxxxy[i] * fe_0 + tr_x_xxxxy_xxxxyy[i] * pa_y[i];

        tr_x_xxxxyy_xxxxyz[i] = tr_x_xxxx_xxxxyz[i] * fe_0 + tr_x_xxxxy_xxxxz[i] * fe_0 + tr_x_xxxxy_xxxxyz[i] * pa_y[i];

        tr_x_xxxxyy_xxxxzz[i] = tr_x_xxxx_xxxxzz[i] * fe_0 + tr_x_xxxxy_xxxxzz[i] * pa_y[i];

        tr_x_xxxxyy_xxxyyy[i] = tr_x_xxxx_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxxxy_xxxyy[i] * fe_0 + tr_x_xxxxy_xxxyyy[i] * pa_y[i];

        tr_x_xxxxyy_xxxyyz[i] = tr_x_xxxx_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxxxy_xxxyz[i] * fe_0 + tr_x_xxxxy_xxxyyz[i] * pa_y[i];

        tr_x_xxxxyy_xxxyzz[i] = tr_x_xxxx_xxxyzz[i] * fe_0 + tr_x_xxxxy_xxxzz[i] * fe_0 + tr_x_xxxxy_xxxyzz[i] * pa_y[i];

        tr_x_xxxxyy_xxxzzz[i] = tr_x_xxxx_xxxzzz[i] * fe_0 + tr_x_xxxxy_xxxzzz[i] * pa_y[i];

        tr_x_xxxxyy_xxyyyy[i] = tr_x_xxxx_xxyyyy[i] * fe_0 + 4.0 * tr_x_xxxxy_xxyyy[i] * fe_0 + tr_x_xxxxy_xxyyyy[i] * pa_y[i];

        tr_x_xxxxyy_xxyyyz[i] = tr_x_xxxx_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxxxy_xxyyz[i] * fe_0 + tr_x_xxxxy_xxyyyz[i] * pa_y[i];

        tr_x_xxxxyy_xxyyzz[i] = tr_x_xxxx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxxy_xxyzz[i] * fe_0 + tr_x_xxxxy_xxyyzz[i] * pa_y[i];

        tr_x_xxxxyy_xxyzzz[i] = tr_x_xxxx_xxyzzz[i] * fe_0 + tr_x_xxxxy_xxzzz[i] * fe_0 + tr_x_xxxxy_xxyzzz[i] * pa_y[i];

        tr_x_xxxxyy_xxzzzz[i] = tr_x_xxxx_xxzzzz[i] * fe_0 + tr_x_xxxxy_xxzzzz[i] * pa_y[i];

        tr_x_xxxxyy_xyyyyy[i] = tr_x_xxxx_xyyyyy[i] * fe_0 + 5.0 * tr_x_xxxxy_xyyyy[i] * fe_0 + tr_x_xxxxy_xyyyyy[i] * pa_y[i];

        tr_x_xxxxyy_xyyyyz[i] = tr_x_xxxx_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxxxy_xyyyz[i] * fe_0 + tr_x_xxxxy_xyyyyz[i] * pa_y[i];

        tr_x_xxxxyy_xyyyzz[i] = tr_x_xxxx_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxxxy_xyyzz[i] * fe_0 + tr_x_xxxxy_xyyyzz[i] * pa_y[i];

        tr_x_xxxxyy_xyyzzz[i] = tr_x_xxxx_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxxxy_xyzzz[i] * fe_0 + tr_x_xxxxy_xyyzzz[i] * pa_y[i];

        tr_x_xxxxyy_xyzzzz[i] = tr_x_xxxx_xyzzzz[i] * fe_0 + tr_x_xxxxy_xzzzz[i] * fe_0 + tr_x_xxxxy_xyzzzz[i] * pa_y[i];

        tr_x_xxxxyy_xzzzzz[i] = tr_x_xxxx_xzzzzz[i] * fe_0 + tr_x_xxxxy_xzzzzz[i] * pa_y[i];

        tr_x_xxxxyy_yyyyyy[i] = 3.0 * tr_x_xxyy_yyyyyy[i] * fe_0 + ts_xxxyy_yyyyyy[i] * fe_0 + tr_x_xxxyy_yyyyyy[i] * pa_x[i];

        tr_x_xxxxyy_yyyyyz[i] = 3.0 * tr_x_xxyy_yyyyyz[i] * fe_0 + ts_xxxyy_yyyyyz[i] * fe_0 + tr_x_xxxyy_yyyyyz[i] * pa_x[i];

        tr_x_xxxxyy_yyyyzz[i] = 3.0 * tr_x_xxyy_yyyyzz[i] * fe_0 + ts_xxxyy_yyyyzz[i] * fe_0 + tr_x_xxxyy_yyyyzz[i] * pa_x[i];

        tr_x_xxxxyy_yyyzzz[i] = 3.0 * tr_x_xxyy_yyyzzz[i] * fe_0 + ts_xxxyy_yyyzzz[i] * fe_0 + tr_x_xxxyy_yyyzzz[i] * pa_x[i];

        tr_x_xxxxyy_yyzzzz[i] = 3.0 * tr_x_xxyy_yyzzzz[i] * fe_0 + ts_xxxyy_yyzzzz[i] * fe_0 + tr_x_xxxyy_yyzzzz[i] * pa_x[i];

        tr_x_xxxxyy_yzzzzz[i] = 3.0 * tr_x_xxyy_yzzzzz[i] * fe_0 + ts_xxxyy_yzzzzz[i] * fe_0 + tr_x_xxxyy_yzzzzz[i] * pa_x[i];

        tr_x_xxxxyy_zzzzzz[i] = tr_x_xxxx_zzzzzz[i] * fe_0 + tr_x_xxxxy_zzzzzz[i] * pa_y[i];
    }

    // Set up 112-140 components of targeted buffer : II

    auto tr_x_xxxxyz_xxxxxx = pbuffer.data(idx_dip_ii + 112);

    auto tr_x_xxxxyz_xxxxxy = pbuffer.data(idx_dip_ii + 113);

    auto tr_x_xxxxyz_xxxxxz = pbuffer.data(idx_dip_ii + 114);

    auto tr_x_xxxxyz_xxxxyy = pbuffer.data(idx_dip_ii + 115);

    auto tr_x_xxxxyz_xxxxyz = pbuffer.data(idx_dip_ii + 116);

    auto tr_x_xxxxyz_xxxxzz = pbuffer.data(idx_dip_ii + 117);

    auto tr_x_xxxxyz_xxxyyy = pbuffer.data(idx_dip_ii + 118);

    auto tr_x_xxxxyz_xxxyyz = pbuffer.data(idx_dip_ii + 119);

    auto tr_x_xxxxyz_xxxyzz = pbuffer.data(idx_dip_ii + 120);

    auto tr_x_xxxxyz_xxxzzz = pbuffer.data(idx_dip_ii + 121);

    auto tr_x_xxxxyz_xxyyyy = pbuffer.data(idx_dip_ii + 122);

    auto tr_x_xxxxyz_xxyyyz = pbuffer.data(idx_dip_ii + 123);

    auto tr_x_xxxxyz_xxyyzz = pbuffer.data(idx_dip_ii + 124);

    auto tr_x_xxxxyz_xxyzzz = pbuffer.data(idx_dip_ii + 125);

    auto tr_x_xxxxyz_xxzzzz = pbuffer.data(idx_dip_ii + 126);

    auto tr_x_xxxxyz_xyyyyy = pbuffer.data(idx_dip_ii + 127);

    auto tr_x_xxxxyz_xyyyyz = pbuffer.data(idx_dip_ii + 128);

    auto tr_x_xxxxyz_xyyyzz = pbuffer.data(idx_dip_ii + 129);

    auto tr_x_xxxxyz_xyyzzz = pbuffer.data(idx_dip_ii + 130);

    auto tr_x_xxxxyz_xyzzzz = pbuffer.data(idx_dip_ii + 131);

    auto tr_x_xxxxyz_xzzzzz = pbuffer.data(idx_dip_ii + 132);

    auto tr_x_xxxxyz_yyyyyy = pbuffer.data(idx_dip_ii + 133);

    auto tr_x_xxxxyz_yyyyyz = pbuffer.data(idx_dip_ii + 134);

    auto tr_x_xxxxyz_yyyyzz = pbuffer.data(idx_dip_ii + 135);

    auto tr_x_xxxxyz_yyyzzz = pbuffer.data(idx_dip_ii + 136);

    auto tr_x_xxxxyz_yyzzzz = pbuffer.data(idx_dip_ii + 137);

    auto tr_x_xxxxyz_yzzzzz = pbuffer.data(idx_dip_ii + 138);

    auto tr_x_xxxxyz_zzzzzz = pbuffer.data(idx_dip_ii + 139);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxxy_xxxxxy, tr_x_xxxxy_xxxxyy, tr_x_xxxxy_xxxyyy, tr_x_xxxxy_xxyyyy, tr_x_xxxxy_xyyyyy, tr_x_xxxxy_yyyyyy, tr_x_xxxxyz_xxxxxx, tr_x_xxxxyz_xxxxxy, tr_x_xxxxyz_xxxxxz, tr_x_xxxxyz_xxxxyy, tr_x_xxxxyz_xxxxyz, tr_x_xxxxyz_xxxxzz, tr_x_xxxxyz_xxxyyy, tr_x_xxxxyz_xxxyyz, tr_x_xxxxyz_xxxyzz, tr_x_xxxxyz_xxxzzz, tr_x_xxxxyz_xxyyyy, tr_x_xxxxyz_xxyyyz, tr_x_xxxxyz_xxyyzz, tr_x_xxxxyz_xxyzzz, tr_x_xxxxyz_xxzzzz, tr_x_xxxxyz_xyyyyy, tr_x_xxxxyz_xyyyyz, tr_x_xxxxyz_xyyyzz, tr_x_xxxxyz_xyyzzz, tr_x_xxxxyz_xyzzzz, tr_x_xxxxyz_xzzzzz, tr_x_xxxxyz_yyyyyy, tr_x_xxxxyz_yyyyyz, tr_x_xxxxyz_yyyyzz, tr_x_xxxxyz_yyyzzz, tr_x_xxxxyz_yyzzzz, tr_x_xxxxyz_yzzzzz, tr_x_xxxxyz_zzzzzz, tr_x_xxxxz_xxxxxx, tr_x_xxxxz_xxxxxz, tr_x_xxxxz_xxxxyz, tr_x_xxxxz_xxxxz, tr_x_xxxxz_xxxxzz, tr_x_xxxxz_xxxyyz, tr_x_xxxxz_xxxyz, tr_x_xxxxz_xxxyzz, tr_x_xxxxz_xxxzz, tr_x_xxxxz_xxxzzz, tr_x_xxxxz_xxyyyz, tr_x_xxxxz_xxyyz, tr_x_xxxxz_xxyyzz, tr_x_xxxxz_xxyzz, tr_x_xxxxz_xxyzzz, tr_x_xxxxz_xxzzz, tr_x_xxxxz_xxzzzz, tr_x_xxxxz_xyyyyz, tr_x_xxxxz_xyyyz, tr_x_xxxxz_xyyyzz, tr_x_xxxxz_xyyzz, tr_x_xxxxz_xyyzzz, tr_x_xxxxz_xyzzz, tr_x_xxxxz_xyzzzz, tr_x_xxxxz_xzzzz, tr_x_xxxxz_xzzzzz, tr_x_xxxxz_yyyyyz, tr_x_xxxxz_yyyyz, tr_x_xxxxz_yyyyzz, tr_x_xxxxz_yyyzz, tr_x_xxxxz_yyyzzz, tr_x_xxxxz_yyzzz, tr_x_xxxxz_yyzzzz, tr_x_xxxxz_yzzzz, tr_x_xxxxz_yzzzzz, tr_x_xxxxz_zzzzz, tr_x_xxxxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyz_xxxxxx[i] = tr_x_xxxxz_xxxxxx[i] * pa_y[i];

        tr_x_xxxxyz_xxxxxy[i] = tr_x_xxxxy_xxxxxy[i] * pa_z[i];

        tr_x_xxxxyz_xxxxxz[i] = tr_x_xxxxz_xxxxxz[i] * pa_y[i];

        tr_x_xxxxyz_xxxxyy[i] = tr_x_xxxxy_xxxxyy[i] * pa_z[i];

        tr_x_xxxxyz_xxxxyz[i] = tr_x_xxxxz_xxxxz[i] * fe_0 + tr_x_xxxxz_xxxxyz[i] * pa_y[i];

        tr_x_xxxxyz_xxxxzz[i] = tr_x_xxxxz_xxxxzz[i] * pa_y[i];

        tr_x_xxxxyz_xxxyyy[i] = tr_x_xxxxy_xxxyyy[i] * pa_z[i];

        tr_x_xxxxyz_xxxyyz[i] = 2.0 * tr_x_xxxxz_xxxyz[i] * fe_0 + tr_x_xxxxz_xxxyyz[i] * pa_y[i];

        tr_x_xxxxyz_xxxyzz[i] = tr_x_xxxxz_xxxzz[i] * fe_0 + tr_x_xxxxz_xxxyzz[i] * pa_y[i];

        tr_x_xxxxyz_xxxzzz[i] = tr_x_xxxxz_xxxzzz[i] * pa_y[i];

        tr_x_xxxxyz_xxyyyy[i] = tr_x_xxxxy_xxyyyy[i] * pa_z[i];

        tr_x_xxxxyz_xxyyyz[i] = 3.0 * tr_x_xxxxz_xxyyz[i] * fe_0 + tr_x_xxxxz_xxyyyz[i] * pa_y[i];

        tr_x_xxxxyz_xxyyzz[i] = 2.0 * tr_x_xxxxz_xxyzz[i] * fe_0 + tr_x_xxxxz_xxyyzz[i] * pa_y[i];

        tr_x_xxxxyz_xxyzzz[i] = tr_x_xxxxz_xxzzz[i] * fe_0 + tr_x_xxxxz_xxyzzz[i] * pa_y[i];

        tr_x_xxxxyz_xxzzzz[i] = tr_x_xxxxz_xxzzzz[i] * pa_y[i];

        tr_x_xxxxyz_xyyyyy[i] = tr_x_xxxxy_xyyyyy[i] * pa_z[i];

        tr_x_xxxxyz_xyyyyz[i] = 4.0 * tr_x_xxxxz_xyyyz[i] * fe_0 + tr_x_xxxxz_xyyyyz[i] * pa_y[i];

        tr_x_xxxxyz_xyyyzz[i] = 3.0 * tr_x_xxxxz_xyyzz[i] * fe_0 + tr_x_xxxxz_xyyyzz[i] * pa_y[i];

        tr_x_xxxxyz_xyyzzz[i] = 2.0 * tr_x_xxxxz_xyzzz[i] * fe_0 + tr_x_xxxxz_xyyzzz[i] * pa_y[i];

        tr_x_xxxxyz_xyzzzz[i] = tr_x_xxxxz_xzzzz[i] * fe_0 + tr_x_xxxxz_xyzzzz[i] * pa_y[i];

        tr_x_xxxxyz_xzzzzz[i] = tr_x_xxxxz_xzzzzz[i] * pa_y[i];

        tr_x_xxxxyz_yyyyyy[i] = tr_x_xxxxy_yyyyyy[i] * pa_z[i];

        tr_x_xxxxyz_yyyyyz[i] = 5.0 * tr_x_xxxxz_yyyyz[i] * fe_0 + tr_x_xxxxz_yyyyyz[i] * pa_y[i];

        tr_x_xxxxyz_yyyyzz[i] = 4.0 * tr_x_xxxxz_yyyzz[i] * fe_0 + tr_x_xxxxz_yyyyzz[i] * pa_y[i];

        tr_x_xxxxyz_yyyzzz[i] = 3.0 * tr_x_xxxxz_yyzzz[i] * fe_0 + tr_x_xxxxz_yyyzzz[i] * pa_y[i];

        tr_x_xxxxyz_yyzzzz[i] = 2.0 * tr_x_xxxxz_yzzzz[i] * fe_0 + tr_x_xxxxz_yyzzzz[i] * pa_y[i];

        tr_x_xxxxyz_yzzzzz[i] = tr_x_xxxxz_zzzzz[i] * fe_0 + tr_x_xxxxz_yzzzzz[i] * pa_y[i];

        tr_x_xxxxyz_zzzzzz[i] = tr_x_xxxxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : II

    auto tr_x_xxxxzz_xxxxxx = pbuffer.data(idx_dip_ii + 140);

    auto tr_x_xxxxzz_xxxxxy = pbuffer.data(idx_dip_ii + 141);

    auto tr_x_xxxxzz_xxxxxz = pbuffer.data(idx_dip_ii + 142);

    auto tr_x_xxxxzz_xxxxyy = pbuffer.data(idx_dip_ii + 143);

    auto tr_x_xxxxzz_xxxxyz = pbuffer.data(idx_dip_ii + 144);

    auto tr_x_xxxxzz_xxxxzz = pbuffer.data(idx_dip_ii + 145);

    auto tr_x_xxxxzz_xxxyyy = pbuffer.data(idx_dip_ii + 146);

    auto tr_x_xxxxzz_xxxyyz = pbuffer.data(idx_dip_ii + 147);

    auto tr_x_xxxxzz_xxxyzz = pbuffer.data(idx_dip_ii + 148);

    auto tr_x_xxxxzz_xxxzzz = pbuffer.data(idx_dip_ii + 149);

    auto tr_x_xxxxzz_xxyyyy = pbuffer.data(idx_dip_ii + 150);

    auto tr_x_xxxxzz_xxyyyz = pbuffer.data(idx_dip_ii + 151);

    auto tr_x_xxxxzz_xxyyzz = pbuffer.data(idx_dip_ii + 152);

    auto tr_x_xxxxzz_xxyzzz = pbuffer.data(idx_dip_ii + 153);

    auto tr_x_xxxxzz_xxzzzz = pbuffer.data(idx_dip_ii + 154);

    auto tr_x_xxxxzz_xyyyyy = pbuffer.data(idx_dip_ii + 155);

    auto tr_x_xxxxzz_xyyyyz = pbuffer.data(idx_dip_ii + 156);

    auto tr_x_xxxxzz_xyyyzz = pbuffer.data(idx_dip_ii + 157);

    auto tr_x_xxxxzz_xyyzzz = pbuffer.data(idx_dip_ii + 158);

    auto tr_x_xxxxzz_xyzzzz = pbuffer.data(idx_dip_ii + 159);

    auto tr_x_xxxxzz_xzzzzz = pbuffer.data(idx_dip_ii + 160);

    auto tr_x_xxxxzz_yyyyyy = pbuffer.data(idx_dip_ii + 161);

    auto tr_x_xxxxzz_yyyyyz = pbuffer.data(idx_dip_ii + 162);

    auto tr_x_xxxxzz_yyyyzz = pbuffer.data(idx_dip_ii + 163);

    auto tr_x_xxxxzz_yyyzzz = pbuffer.data(idx_dip_ii + 164);

    auto tr_x_xxxxzz_yyzzzz = pbuffer.data(idx_dip_ii + 165);

    auto tr_x_xxxxzz_yzzzzz = pbuffer.data(idx_dip_ii + 166);

    auto tr_x_xxxxzz_zzzzzz = pbuffer.data(idx_dip_ii + 167);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxxx_xxxxxx, tr_x_xxxx_xxxxxy, tr_x_xxxx_xxxxxz, tr_x_xxxx_xxxxyy, tr_x_xxxx_xxxxyz, tr_x_xxxx_xxxxzz, tr_x_xxxx_xxxyyy, tr_x_xxxx_xxxyyz, tr_x_xxxx_xxxyzz, tr_x_xxxx_xxxzzz, tr_x_xxxx_xxyyyy, tr_x_xxxx_xxyyyz, tr_x_xxxx_xxyyzz, tr_x_xxxx_xxyzzz, tr_x_xxxx_xxzzzz, tr_x_xxxx_xyyyyy, tr_x_xxxx_xyyyyz, tr_x_xxxx_xyyyzz, tr_x_xxxx_xyyzzz, tr_x_xxxx_xyzzzz, tr_x_xxxx_xzzzzz, tr_x_xxxx_yyyyyy, tr_x_xxxxz_xxxxx, tr_x_xxxxz_xxxxxx, tr_x_xxxxz_xxxxxy, tr_x_xxxxz_xxxxxz, tr_x_xxxxz_xxxxy, tr_x_xxxxz_xxxxyy, tr_x_xxxxz_xxxxyz, tr_x_xxxxz_xxxxz, tr_x_xxxxz_xxxxzz, tr_x_xxxxz_xxxyy, tr_x_xxxxz_xxxyyy, tr_x_xxxxz_xxxyyz, tr_x_xxxxz_xxxyz, tr_x_xxxxz_xxxyzz, tr_x_xxxxz_xxxzz, tr_x_xxxxz_xxxzzz, tr_x_xxxxz_xxyyy, tr_x_xxxxz_xxyyyy, tr_x_xxxxz_xxyyyz, tr_x_xxxxz_xxyyz, tr_x_xxxxz_xxyyzz, tr_x_xxxxz_xxyzz, tr_x_xxxxz_xxyzzz, tr_x_xxxxz_xxzzz, tr_x_xxxxz_xxzzzz, tr_x_xxxxz_xyyyy, tr_x_xxxxz_xyyyyy, tr_x_xxxxz_xyyyyz, tr_x_xxxxz_xyyyz, tr_x_xxxxz_xyyyzz, tr_x_xxxxz_xyyzz, tr_x_xxxxz_xyyzzz, tr_x_xxxxz_xyzzz, tr_x_xxxxz_xyzzzz, tr_x_xxxxz_xzzzz, tr_x_xxxxz_xzzzzz, tr_x_xxxxz_yyyyyy, tr_x_xxxxzz_xxxxxx, tr_x_xxxxzz_xxxxxy, tr_x_xxxxzz_xxxxxz, tr_x_xxxxzz_xxxxyy, tr_x_xxxxzz_xxxxyz, tr_x_xxxxzz_xxxxzz, tr_x_xxxxzz_xxxyyy, tr_x_xxxxzz_xxxyyz, tr_x_xxxxzz_xxxyzz, tr_x_xxxxzz_xxxzzz, tr_x_xxxxzz_xxyyyy, tr_x_xxxxzz_xxyyyz, tr_x_xxxxzz_xxyyzz, tr_x_xxxxzz_xxyzzz, tr_x_xxxxzz_xxzzzz, tr_x_xxxxzz_xyyyyy, tr_x_xxxxzz_xyyyyz, tr_x_xxxxzz_xyyyzz, tr_x_xxxxzz_xyyzzz, tr_x_xxxxzz_xyzzzz, tr_x_xxxxzz_xzzzzz, tr_x_xxxxzz_yyyyyy, tr_x_xxxxzz_yyyyyz, tr_x_xxxxzz_yyyyzz, tr_x_xxxxzz_yyyzzz, tr_x_xxxxzz_yyzzzz, tr_x_xxxxzz_yzzzzz, tr_x_xxxxzz_zzzzzz, tr_x_xxxzz_yyyyyz, tr_x_xxxzz_yyyyzz, tr_x_xxxzz_yyyzzz, tr_x_xxxzz_yyzzzz, tr_x_xxxzz_yzzzzz, tr_x_xxxzz_zzzzzz, tr_x_xxzz_yyyyyz, tr_x_xxzz_yyyyzz, tr_x_xxzz_yyyzzz, tr_x_xxzz_yyzzzz, tr_x_xxzz_yzzzzz, tr_x_xxzz_zzzzzz, ts_xxxzz_yyyyyz, ts_xxxzz_yyyyzz, ts_xxxzz_yyyzzz, ts_xxxzz_yyzzzz, ts_xxxzz_yzzzzz, ts_xxxzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxzz_xxxxxx[i] = tr_x_xxxx_xxxxxx[i] * fe_0 + tr_x_xxxxz_xxxxxx[i] * pa_z[i];

        tr_x_xxxxzz_xxxxxy[i] = tr_x_xxxx_xxxxxy[i] * fe_0 + tr_x_xxxxz_xxxxxy[i] * pa_z[i];

        tr_x_xxxxzz_xxxxxz[i] = tr_x_xxxx_xxxxxz[i] * fe_0 + tr_x_xxxxz_xxxxx[i] * fe_0 + tr_x_xxxxz_xxxxxz[i] * pa_z[i];

        tr_x_xxxxzz_xxxxyy[i] = tr_x_xxxx_xxxxyy[i] * fe_0 + tr_x_xxxxz_xxxxyy[i] * pa_z[i];

        tr_x_xxxxzz_xxxxyz[i] = tr_x_xxxx_xxxxyz[i] * fe_0 + tr_x_xxxxz_xxxxy[i] * fe_0 + tr_x_xxxxz_xxxxyz[i] * pa_z[i];

        tr_x_xxxxzz_xxxxzz[i] = tr_x_xxxx_xxxxzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xxxxz[i] * fe_0 + tr_x_xxxxz_xxxxzz[i] * pa_z[i];

        tr_x_xxxxzz_xxxyyy[i] = tr_x_xxxx_xxxyyy[i] * fe_0 + tr_x_xxxxz_xxxyyy[i] * pa_z[i];

        tr_x_xxxxzz_xxxyyz[i] = tr_x_xxxx_xxxyyz[i] * fe_0 + tr_x_xxxxz_xxxyy[i] * fe_0 + tr_x_xxxxz_xxxyyz[i] * pa_z[i];

        tr_x_xxxxzz_xxxyzz[i] = tr_x_xxxx_xxxyzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xxxyz[i] * fe_0 + tr_x_xxxxz_xxxyzz[i] * pa_z[i];

        tr_x_xxxxzz_xxxzzz[i] = tr_x_xxxx_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxxxz_xxxzz[i] * fe_0 + tr_x_xxxxz_xxxzzz[i] * pa_z[i];

        tr_x_xxxxzz_xxyyyy[i] = tr_x_xxxx_xxyyyy[i] * fe_0 + tr_x_xxxxz_xxyyyy[i] * pa_z[i];

        tr_x_xxxxzz_xxyyyz[i] = tr_x_xxxx_xxyyyz[i] * fe_0 + tr_x_xxxxz_xxyyy[i] * fe_0 + tr_x_xxxxz_xxyyyz[i] * pa_z[i];

        tr_x_xxxxzz_xxyyzz[i] = tr_x_xxxx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xxyyz[i] * fe_0 + tr_x_xxxxz_xxyyzz[i] * pa_z[i];

        tr_x_xxxxzz_xxyzzz[i] = tr_x_xxxx_xxyzzz[i] * fe_0 + 3.0 * tr_x_xxxxz_xxyzz[i] * fe_0 + tr_x_xxxxz_xxyzzz[i] * pa_z[i];

        tr_x_xxxxzz_xxzzzz[i] = tr_x_xxxx_xxzzzz[i] * fe_0 + 4.0 * tr_x_xxxxz_xxzzz[i] * fe_0 + tr_x_xxxxz_xxzzzz[i] * pa_z[i];

        tr_x_xxxxzz_xyyyyy[i] = tr_x_xxxx_xyyyyy[i] * fe_0 + tr_x_xxxxz_xyyyyy[i] * pa_z[i];

        tr_x_xxxxzz_xyyyyz[i] = tr_x_xxxx_xyyyyz[i] * fe_0 + tr_x_xxxxz_xyyyy[i] * fe_0 + tr_x_xxxxz_xyyyyz[i] * pa_z[i];

        tr_x_xxxxzz_xyyyzz[i] = tr_x_xxxx_xyyyzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xyyyz[i] * fe_0 + tr_x_xxxxz_xyyyzz[i] * pa_z[i];

        tr_x_xxxxzz_xyyzzz[i] = tr_x_xxxx_xyyzzz[i] * fe_0 + 3.0 * tr_x_xxxxz_xyyzz[i] * fe_0 + tr_x_xxxxz_xyyzzz[i] * pa_z[i];

        tr_x_xxxxzz_xyzzzz[i] = tr_x_xxxx_xyzzzz[i] * fe_0 + 4.0 * tr_x_xxxxz_xyzzz[i] * fe_0 + tr_x_xxxxz_xyzzzz[i] * pa_z[i];

        tr_x_xxxxzz_xzzzzz[i] = tr_x_xxxx_xzzzzz[i] * fe_0 + 5.0 * tr_x_xxxxz_xzzzz[i] * fe_0 + tr_x_xxxxz_xzzzzz[i] * pa_z[i];

        tr_x_xxxxzz_yyyyyy[i] = tr_x_xxxx_yyyyyy[i] * fe_0 + tr_x_xxxxz_yyyyyy[i] * pa_z[i];

        tr_x_xxxxzz_yyyyyz[i] = 3.0 * tr_x_xxzz_yyyyyz[i] * fe_0 + ts_xxxzz_yyyyyz[i] * fe_0 + tr_x_xxxzz_yyyyyz[i] * pa_x[i];

        tr_x_xxxxzz_yyyyzz[i] = 3.0 * tr_x_xxzz_yyyyzz[i] * fe_0 + ts_xxxzz_yyyyzz[i] * fe_0 + tr_x_xxxzz_yyyyzz[i] * pa_x[i];

        tr_x_xxxxzz_yyyzzz[i] = 3.0 * tr_x_xxzz_yyyzzz[i] * fe_0 + ts_xxxzz_yyyzzz[i] * fe_0 + tr_x_xxxzz_yyyzzz[i] * pa_x[i];

        tr_x_xxxxzz_yyzzzz[i] = 3.0 * tr_x_xxzz_yyzzzz[i] * fe_0 + ts_xxxzz_yyzzzz[i] * fe_0 + tr_x_xxxzz_yyzzzz[i] * pa_x[i];

        tr_x_xxxxzz_yzzzzz[i] = 3.0 * tr_x_xxzz_yzzzzz[i] * fe_0 + ts_xxxzz_yzzzzz[i] * fe_0 + tr_x_xxxzz_yzzzzz[i] * pa_x[i];

        tr_x_xxxxzz_zzzzzz[i] = 3.0 * tr_x_xxzz_zzzzzz[i] * fe_0 + ts_xxxzz_zzzzzz[i] * fe_0 + tr_x_xxxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : II

    auto tr_x_xxxyyy_xxxxxx = pbuffer.data(idx_dip_ii + 168);

    auto tr_x_xxxyyy_xxxxxy = pbuffer.data(idx_dip_ii + 169);

    auto tr_x_xxxyyy_xxxxxz = pbuffer.data(idx_dip_ii + 170);

    auto tr_x_xxxyyy_xxxxyy = pbuffer.data(idx_dip_ii + 171);

    auto tr_x_xxxyyy_xxxxyz = pbuffer.data(idx_dip_ii + 172);

    auto tr_x_xxxyyy_xxxxzz = pbuffer.data(idx_dip_ii + 173);

    auto tr_x_xxxyyy_xxxyyy = pbuffer.data(idx_dip_ii + 174);

    auto tr_x_xxxyyy_xxxyyz = pbuffer.data(idx_dip_ii + 175);

    auto tr_x_xxxyyy_xxxyzz = pbuffer.data(idx_dip_ii + 176);

    auto tr_x_xxxyyy_xxxzzz = pbuffer.data(idx_dip_ii + 177);

    auto tr_x_xxxyyy_xxyyyy = pbuffer.data(idx_dip_ii + 178);

    auto tr_x_xxxyyy_xxyyyz = pbuffer.data(idx_dip_ii + 179);

    auto tr_x_xxxyyy_xxyyzz = pbuffer.data(idx_dip_ii + 180);

    auto tr_x_xxxyyy_xxyzzz = pbuffer.data(idx_dip_ii + 181);

    auto tr_x_xxxyyy_xxzzzz = pbuffer.data(idx_dip_ii + 182);

    auto tr_x_xxxyyy_xyyyyy = pbuffer.data(idx_dip_ii + 183);

    auto tr_x_xxxyyy_xyyyyz = pbuffer.data(idx_dip_ii + 184);

    auto tr_x_xxxyyy_xyyyzz = pbuffer.data(idx_dip_ii + 185);

    auto tr_x_xxxyyy_xyyzzz = pbuffer.data(idx_dip_ii + 186);

    auto tr_x_xxxyyy_xyzzzz = pbuffer.data(idx_dip_ii + 187);

    auto tr_x_xxxyyy_xzzzzz = pbuffer.data(idx_dip_ii + 188);

    auto tr_x_xxxyyy_yyyyyy = pbuffer.data(idx_dip_ii + 189);

    auto tr_x_xxxyyy_yyyyyz = pbuffer.data(idx_dip_ii + 190);

    auto tr_x_xxxyyy_yyyyzz = pbuffer.data(idx_dip_ii + 191);

    auto tr_x_xxxyyy_yyyzzz = pbuffer.data(idx_dip_ii + 192);

    auto tr_x_xxxyyy_yyzzzz = pbuffer.data(idx_dip_ii + 193);

    auto tr_x_xxxyyy_yzzzzz = pbuffer.data(idx_dip_ii + 194);

    auto tr_x_xxxyyy_zzzzzz = pbuffer.data(idx_dip_ii + 195);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxxy_xxxxxx, tr_x_xxxy_xxxxxy, tr_x_xxxy_xxxxxz, tr_x_xxxy_xxxxyy, tr_x_xxxy_xxxxyz, tr_x_xxxy_xxxxzz, tr_x_xxxy_xxxyyy, tr_x_xxxy_xxxyyz, tr_x_xxxy_xxxyzz, tr_x_xxxy_xxxzzz, tr_x_xxxy_xxyyyy, tr_x_xxxy_xxyyyz, tr_x_xxxy_xxyyzz, tr_x_xxxy_xxyzzz, tr_x_xxxy_xxzzzz, tr_x_xxxy_xyyyyy, tr_x_xxxy_xyyyyz, tr_x_xxxy_xyyyzz, tr_x_xxxy_xyyzzz, tr_x_xxxy_xyzzzz, tr_x_xxxy_xzzzzz, tr_x_xxxy_zzzzzz, tr_x_xxxyy_xxxxx, tr_x_xxxyy_xxxxxx, tr_x_xxxyy_xxxxxy, tr_x_xxxyy_xxxxxz, tr_x_xxxyy_xxxxy, tr_x_xxxyy_xxxxyy, tr_x_xxxyy_xxxxyz, tr_x_xxxyy_xxxxz, tr_x_xxxyy_xxxxzz, tr_x_xxxyy_xxxyy, tr_x_xxxyy_xxxyyy, tr_x_xxxyy_xxxyyz, tr_x_xxxyy_xxxyz, tr_x_xxxyy_xxxyzz, tr_x_xxxyy_xxxzz, tr_x_xxxyy_xxxzzz, tr_x_xxxyy_xxyyy, tr_x_xxxyy_xxyyyy, tr_x_xxxyy_xxyyyz, tr_x_xxxyy_xxyyz, tr_x_xxxyy_xxyyzz, tr_x_xxxyy_xxyzz, tr_x_xxxyy_xxyzzz, tr_x_xxxyy_xxzzz, tr_x_xxxyy_xxzzzz, tr_x_xxxyy_xyyyy, tr_x_xxxyy_xyyyyy, tr_x_xxxyy_xyyyyz, tr_x_xxxyy_xyyyz, tr_x_xxxyy_xyyyzz, tr_x_xxxyy_xyyzz, tr_x_xxxyy_xyyzzz, tr_x_xxxyy_xyzzz, tr_x_xxxyy_xyzzzz, tr_x_xxxyy_xzzzz, tr_x_xxxyy_xzzzzz, tr_x_xxxyy_zzzzzz, tr_x_xxxyyy_xxxxxx, tr_x_xxxyyy_xxxxxy, tr_x_xxxyyy_xxxxxz, tr_x_xxxyyy_xxxxyy, tr_x_xxxyyy_xxxxyz, tr_x_xxxyyy_xxxxzz, tr_x_xxxyyy_xxxyyy, tr_x_xxxyyy_xxxyyz, tr_x_xxxyyy_xxxyzz, tr_x_xxxyyy_xxxzzz, tr_x_xxxyyy_xxyyyy, tr_x_xxxyyy_xxyyyz, tr_x_xxxyyy_xxyyzz, tr_x_xxxyyy_xxyzzz, tr_x_xxxyyy_xxzzzz, tr_x_xxxyyy_xyyyyy, tr_x_xxxyyy_xyyyyz, tr_x_xxxyyy_xyyyzz, tr_x_xxxyyy_xyyzzz, tr_x_xxxyyy_xyzzzz, tr_x_xxxyyy_xzzzzz, tr_x_xxxyyy_yyyyyy, tr_x_xxxyyy_yyyyyz, tr_x_xxxyyy_yyyyzz, tr_x_xxxyyy_yyyzzz, tr_x_xxxyyy_yyzzzz, tr_x_xxxyyy_yzzzzz, tr_x_xxxyyy_zzzzzz, tr_x_xxyyy_yyyyyy, tr_x_xxyyy_yyyyyz, tr_x_xxyyy_yyyyzz, tr_x_xxyyy_yyyzzz, tr_x_xxyyy_yyzzzz, tr_x_xxyyy_yzzzzz, tr_x_xyyy_yyyyyy, tr_x_xyyy_yyyyyz, tr_x_xyyy_yyyyzz, tr_x_xyyy_yyyzzz, tr_x_xyyy_yyzzzz, tr_x_xyyy_yzzzzz, ts_xxyyy_yyyyyy, ts_xxyyy_yyyyyz, ts_xxyyy_yyyyzz, ts_xxyyy_yyyzzz, ts_xxyyy_yyzzzz, ts_xxyyy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyy_xxxxxx[i] = 2.0 * tr_x_xxxy_xxxxxx[i] * fe_0 + tr_x_xxxyy_xxxxxx[i] * pa_y[i];

        tr_x_xxxyyy_xxxxxy[i] = 2.0 * tr_x_xxxy_xxxxxy[i] * fe_0 + tr_x_xxxyy_xxxxx[i] * fe_0 + tr_x_xxxyy_xxxxxy[i] * pa_y[i];

        tr_x_xxxyyy_xxxxxz[i] = 2.0 * tr_x_xxxy_xxxxxz[i] * fe_0 + tr_x_xxxyy_xxxxxz[i] * pa_y[i];

        tr_x_xxxyyy_xxxxyy[i] = 2.0 * tr_x_xxxy_xxxxyy[i] * fe_0 + 2.0 * tr_x_xxxyy_xxxxy[i] * fe_0 + tr_x_xxxyy_xxxxyy[i] * pa_y[i];

        tr_x_xxxyyy_xxxxyz[i] = 2.0 * tr_x_xxxy_xxxxyz[i] * fe_0 + tr_x_xxxyy_xxxxz[i] * fe_0 + tr_x_xxxyy_xxxxyz[i] * pa_y[i];

        tr_x_xxxyyy_xxxxzz[i] = 2.0 * tr_x_xxxy_xxxxzz[i] * fe_0 + tr_x_xxxyy_xxxxzz[i] * pa_y[i];

        tr_x_xxxyyy_xxxyyy[i] = 2.0 * tr_x_xxxy_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxxyy_xxxyy[i] * fe_0 + tr_x_xxxyy_xxxyyy[i] * pa_y[i];

        tr_x_xxxyyy_xxxyyz[i] = 2.0 * tr_x_xxxy_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxxyy_xxxyz[i] * fe_0 + tr_x_xxxyy_xxxyyz[i] * pa_y[i];

        tr_x_xxxyyy_xxxyzz[i] = 2.0 * tr_x_xxxy_xxxyzz[i] * fe_0 + tr_x_xxxyy_xxxzz[i] * fe_0 + tr_x_xxxyy_xxxyzz[i] * pa_y[i];

        tr_x_xxxyyy_xxxzzz[i] = 2.0 * tr_x_xxxy_xxxzzz[i] * fe_0 + tr_x_xxxyy_xxxzzz[i] * pa_y[i];

        tr_x_xxxyyy_xxyyyy[i] = 2.0 * tr_x_xxxy_xxyyyy[i] * fe_0 + 4.0 * tr_x_xxxyy_xxyyy[i] * fe_0 + tr_x_xxxyy_xxyyyy[i] * pa_y[i];

        tr_x_xxxyyy_xxyyyz[i] = 2.0 * tr_x_xxxy_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxxyy_xxyyz[i] * fe_0 + tr_x_xxxyy_xxyyyz[i] * pa_y[i];

        tr_x_xxxyyy_xxyyzz[i] = 2.0 * tr_x_xxxy_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxyy_xxyzz[i] * fe_0 + tr_x_xxxyy_xxyyzz[i] * pa_y[i];

        tr_x_xxxyyy_xxyzzz[i] = 2.0 * tr_x_xxxy_xxyzzz[i] * fe_0 + tr_x_xxxyy_xxzzz[i] * fe_0 + tr_x_xxxyy_xxyzzz[i] * pa_y[i];

        tr_x_xxxyyy_xxzzzz[i] = 2.0 * tr_x_xxxy_xxzzzz[i] * fe_0 + tr_x_xxxyy_xxzzzz[i] * pa_y[i];

        tr_x_xxxyyy_xyyyyy[i] = 2.0 * tr_x_xxxy_xyyyyy[i] * fe_0 + 5.0 * tr_x_xxxyy_xyyyy[i] * fe_0 + tr_x_xxxyy_xyyyyy[i] * pa_y[i];

        tr_x_xxxyyy_xyyyyz[i] = 2.0 * tr_x_xxxy_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxxyy_xyyyz[i] * fe_0 + tr_x_xxxyy_xyyyyz[i] * pa_y[i];

        tr_x_xxxyyy_xyyyzz[i] = 2.0 * tr_x_xxxy_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxxyy_xyyzz[i] * fe_0 + tr_x_xxxyy_xyyyzz[i] * pa_y[i];

        tr_x_xxxyyy_xyyzzz[i] = 2.0 * tr_x_xxxy_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxxyy_xyzzz[i] * fe_0 + tr_x_xxxyy_xyyzzz[i] * pa_y[i];

        tr_x_xxxyyy_xyzzzz[i] = 2.0 * tr_x_xxxy_xyzzzz[i] * fe_0 + tr_x_xxxyy_xzzzz[i] * fe_0 + tr_x_xxxyy_xyzzzz[i] * pa_y[i];

        tr_x_xxxyyy_xzzzzz[i] = 2.0 * tr_x_xxxy_xzzzzz[i] * fe_0 + tr_x_xxxyy_xzzzzz[i] * pa_y[i];

        tr_x_xxxyyy_yyyyyy[i] = 2.0 * tr_x_xyyy_yyyyyy[i] * fe_0 + ts_xxyyy_yyyyyy[i] * fe_0 + tr_x_xxyyy_yyyyyy[i] * pa_x[i];

        tr_x_xxxyyy_yyyyyz[i] = 2.0 * tr_x_xyyy_yyyyyz[i] * fe_0 + ts_xxyyy_yyyyyz[i] * fe_0 + tr_x_xxyyy_yyyyyz[i] * pa_x[i];

        tr_x_xxxyyy_yyyyzz[i] = 2.0 * tr_x_xyyy_yyyyzz[i] * fe_0 + ts_xxyyy_yyyyzz[i] * fe_0 + tr_x_xxyyy_yyyyzz[i] * pa_x[i];

        tr_x_xxxyyy_yyyzzz[i] = 2.0 * tr_x_xyyy_yyyzzz[i] * fe_0 + ts_xxyyy_yyyzzz[i] * fe_0 + tr_x_xxyyy_yyyzzz[i] * pa_x[i];

        tr_x_xxxyyy_yyzzzz[i] = 2.0 * tr_x_xyyy_yyzzzz[i] * fe_0 + ts_xxyyy_yyzzzz[i] * fe_0 + tr_x_xxyyy_yyzzzz[i] * pa_x[i];

        tr_x_xxxyyy_yzzzzz[i] = 2.0 * tr_x_xyyy_yzzzzz[i] * fe_0 + ts_xxyyy_yzzzzz[i] * fe_0 + tr_x_xxyyy_yzzzzz[i] * pa_x[i];

        tr_x_xxxyyy_zzzzzz[i] = 2.0 * tr_x_xxxy_zzzzzz[i] * fe_0 + tr_x_xxxyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 196-224 components of targeted buffer : II

    auto tr_x_xxxyyz_xxxxxx = pbuffer.data(idx_dip_ii + 196);

    auto tr_x_xxxyyz_xxxxxy = pbuffer.data(idx_dip_ii + 197);

    auto tr_x_xxxyyz_xxxxxz = pbuffer.data(idx_dip_ii + 198);

    auto tr_x_xxxyyz_xxxxyy = pbuffer.data(idx_dip_ii + 199);

    auto tr_x_xxxyyz_xxxxyz = pbuffer.data(idx_dip_ii + 200);

    auto tr_x_xxxyyz_xxxxzz = pbuffer.data(idx_dip_ii + 201);

    auto tr_x_xxxyyz_xxxyyy = pbuffer.data(idx_dip_ii + 202);

    auto tr_x_xxxyyz_xxxyyz = pbuffer.data(idx_dip_ii + 203);

    auto tr_x_xxxyyz_xxxyzz = pbuffer.data(idx_dip_ii + 204);

    auto tr_x_xxxyyz_xxxzzz = pbuffer.data(idx_dip_ii + 205);

    auto tr_x_xxxyyz_xxyyyy = pbuffer.data(idx_dip_ii + 206);

    auto tr_x_xxxyyz_xxyyyz = pbuffer.data(idx_dip_ii + 207);

    auto tr_x_xxxyyz_xxyyzz = pbuffer.data(idx_dip_ii + 208);

    auto tr_x_xxxyyz_xxyzzz = pbuffer.data(idx_dip_ii + 209);

    auto tr_x_xxxyyz_xxzzzz = pbuffer.data(idx_dip_ii + 210);

    auto tr_x_xxxyyz_xyyyyy = pbuffer.data(idx_dip_ii + 211);

    auto tr_x_xxxyyz_xyyyyz = pbuffer.data(idx_dip_ii + 212);

    auto tr_x_xxxyyz_xyyyzz = pbuffer.data(idx_dip_ii + 213);

    auto tr_x_xxxyyz_xyyzzz = pbuffer.data(idx_dip_ii + 214);

    auto tr_x_xxxyyz_xyzzzz = pbuffer.data(idx_dip_ii + 215);

    auto tr_x_xxxyyz_xzzzzz = pbuffer.data(idx_dip_ii + 216);

    auto tr_x_xxxyyz_yyyyyy = pbuffer.data(idx_dip_ii + 217);

    auto tr_x_xxxyyz_yyyyyz = pbuffer.data(idx_dip_ii + 218);

    auto tr_x_xxxyyz_yyyyzz = pbuffer.data(idx_dip_ii + 219);

    auto tr_x_xxxyyz_yyyzzz = pbuffer.data(idx_dip_ii + 220);

    auto tr_x_xxxyyz_yyzzzz = pbuffer.data(idx_dip_ii + 221);

    auto tr_x_xxxyyz_yzzzzz = pbuffer.data(idx_dip_ii + 222);

    auto tr_x_xxxyyz_zzzzzz = pbuffer.data(idx_dip_ii + 223);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxyy_xxxxxx, tr_x_xxxyy_xxxxxy, tr_x_xxxyy_xxxxy, tr_x_xxxyy_xxxxyy, tr_x_xxxyy_xxxxyz, tr_x_xxxyy_xxxyy, tr_x_xxxyy_xxxyyy, tr_x_xxxyy_xxxyyz, tr_x_xxxyy_xxxyz, tr_x_xxxyy_xxxyzz, tr_x_xxxyy_xxyyy, tr_x_xxxyy_xxyyyy, tr_x_xxxyy_xxyyyz, tr_x_xxxyy_xxyyz, tr_x_xxxyy_xxyyzz, tr_x_xxxyy_xxyzz, tr_x_xxxyy_xxyzzz, tr_x_xxxyy_xyyyy, tr_x_xxxyy_xyyyyy, tr_x_xxxyy_xyyyyz, tr_x_xxxyy_xyyyz, tr_x_xxxyy_xyyyzz, tr_x_xxxyy_xyyzz, tr_x_xxxyy_xyyzzz, tr_x_xxxyy_xyzzz, tr_x_xxxyy_xyzzzz, tr_x_xxxyy_yyyyy, tr_x_xxxyy_yyyyyy, tr_x_xxxyy_yyyyyz, tr_x_xxxyy_yyyyz, tr_x_xxxyy_yyyyzz, tr_x_xxxyy_yyyzz, tr_x_xxxyy_yyyzzz, tr_x_xxxyy_yyzzz, tr_x_xxxyy_yyzzzz, tr_x_xxxyy_yzzzz, tr_x_xxxyy_yzzzzz, tr_x_xxxyyz_xxxxxx, tr_x_xxxyyz_xxxxxy, tr_x_xxxyyz_xxxxxz, tr_x_xxxyyz_xxxxyy, tr_x_xxxyyz_xxxxyz, tr_x_xxxyyz_xxxxzz, tr_x_xxxyyz_xxxyyy, tr_x_xxxyyz_xxxyyz, tr_x_xxxyyz_xxxyzz, tr_x_xxxyyz_xxxzzz, tr_x_xxxyyz_xxyyyy, tr_x_xxxyyz_xxyyyz, tr_x_xxxyyz_xxyyzz, tr_x_xxxyyz_xxyzzz, tr_x_xxxyyz_xxzzzz, tr_x_xxxyyz_xyyyyy, tr_x_xxxyyz_xyyyyz, tr_x_xxxyyz_xyyyzz, tr_x_xxxyyz_xyyzzz, tr_x_xxxyyz_xyzzzz, tr_x_xxxyyz_xzzzzz, tr_x_xxxyyz_yyyyyy, tr_x_xxxyyz_yyyyyz, tr_x_xxxyyz_yyyyzz, tr_x_xxxyyz_yyyzzz, tr_x_xxxyyz_yyzzzz, tr_x_xxxyyz_yzzzzz, tr_x_xxxyyz_zzzzzz, tr_x_xxxyz_xxxxxz, tr_x_xxxyz_xxxxzz, tr_x_xxxyz_xxxzzz, tr_x_xxxyz_xxzzzz, tr_x_xxxyz_xzzzzz, tr_x_xxxyz_zzzzzz, tr_x_xxxz_xxxxxz, tr_x_xxxz_xxxxzz, tr_x_xxxz_xxxzzz, tr_x_xxxz_xxzzzz, tr_x_xxxz_xzzzzz, tr_x_xxxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyz_xxxxxx[i] = tr_x_xxxyy_xxxxxx[i] * pa_z[i];

        tr_x_xxxyyz_xxxxxy[i] = tr_x_xxxyy_xxxxxy[i] * pa_z[i];

        tr_x_xxxyyz_xxxxxz[i] = tr_x_xxxz_xxxxxz[i] * fe_0 + tr_x_xxxyz_xxxxxz[i] * pa_y[i];

        tr_x_xxxyyz_xxxxyy[i] = tr_x_xxxyy_xxxxyy[i] * pa_z[i];

        tr_x_xxxyyz_xxxxyz[i] = tr_x_xxxyy_xxxxy[i] * fe_0 + tr_x_xxxyy_xxxxyz[i] * pa_z[i];

        tr_x_xxxyyz_xxxxzz[i] = tr_x_xxxz_xxxxzz[i] * fe_0 + tr_x_xxxyz_xxxxzz[i] * pa_y[i];

        tr_x_xxxyyz_xxxyyy[i] = tr_x_xxxyy_xxxyyy[i] * pa_z[i];

        tr_x_xxxyyz_xxxyyz[i] = tr_x_xxxyy_xxxyy[i] * fe_0 + tr_x_xxxyy_xxxyyz[i] * pa_z[i];

        tr_x_xxxyyz_xxxyzz[i] = 2.0 * tr_x_xxxyy_xxxyz[i] * fe_0 + tr_x_xxxyy_xxxyzz[i] * pa_z[i];

        tr_x_xxxyyz_xxxzzz[i] = tr_x_xxxz_xxxzzz[i] * fe_0 + tr_x_xxxyz_xxxzzz[i] * pa_y[i];

        tr_x_xxxyyz_xxyyyy[i] = tr_x_xxxyy_xxyyyy[i] * pa_z[i];

        tr_x_xxxyyz_xxyyyz[i] = tr_x_xxxyy_xxyyy[i] * fe_0 + tr_x_xxxyy_xxyyyz[i] * pa_z[i];

        tr_x_xxxyyz_xxyyzz[i] = 2.0 * tr_x_xxxyy_xxyyz[i] * fe_0 + tr_x_xxxyy_xxyyzz[i] * pa_z[i];

        tr_x_xxxyyz_xxyzzz[i] = 3.0 * tr_x_xxxyy_xxyzz[i] * fe_0 + tr_x_xxxyy_xxyzzz[i] * pa_z[i];

        tr_x_xxxyyz_xxzzzz[i] = tr_x_xxxz_xxzzzz[i] * fe_0 + tr_x_xxxyz_xxzzzz[i] * pa_y[i];

        tr_x_xxxyyz_xyyyyy[i] = tr_x_xxxyy_xyyyyy[i] * pa_z[i];

        tr_x_xxxyyz_xyyyyz[i] = tr_x_xxxyy_xyyyy[i] * fe_0 + tr_x_xxxyy_xyyyyz[i] * pa_z[i];

        tr_x_xxxyyz_xyyyzz[i] = 2.0 * tr_x_xxxyy_xyyyz[i] * fe_0 + tr_x_xxxyy_xyyyzz[i] * pa_z[i];

        tr_x_xxxyyz_xyyzzz[i] = 3.0 * tr_x_xxxyy_xyyzz[i] * fe_0 + tr_x_xxxyy_xyyzzz[i] * pa_z[i];

        tr_x_xxxyyz_xyzzzz[i] = 4.0 * tr_x_xxxyy_xyzzz[i] * fe_0 + tr_x_xxxyy_xyzzzz[i] * pa_z[i];

        tr_x_xxxyyz_xzzzzz[i] = tr_x_xxxz_xzzzzz[i] * fe_0 + tr_x_xxxyz_xzzzzz[i] * pa_y[i];

        tr_x_xxxyyz_yyyyyy[i] = tr_x_xxxyy_yyyyyy[i] * pa_z[i];

        tr_x_xxxyyz_yyyyyz[i] = tr_x_xxxyy_yyyyy[i] * fe_0 + tr_x_xxxyy_yyyyyz[i] * pa_z[i];

        tr_x_xxxyyz_yyyyzz[i] = 2.0 * tr_x_xxxyy_yyyyz[i] * fe_0 + tr_x_xxxyy_yyyyzz[i] * pa_z[i];

        tr_x_xxxyyz_yyyzzz[i] = 3.0 * tr_x_xxxyy_yyyzz[i] * fe_0 + tr_x_xxxyy_yyyzzz[i] * pa_z[i];

        tr_x_xxxyyz_yyzzzz[i] = 4.0 * tr_x_xxxyy_yyzzz[i] * fe_0 + tr_x_xxxyy_yyzzzz[i] * pa_z[i];

        tr_x_xxxyyz_yzzzzz[i] = 5.0 * tr_x_xxxyy_yzzzz[i] * fe_0 + tr_x_xxxyy_yzzzzz[i] * pa_z[i];

        tr_x_xxxyyz_zzzzzz[i] = tr_x_xxxz_zzzzzz[i] * fe_0 + tr_x_xxxyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 224-252 components of targeted buffer : II

    auto tr_x_xxxyzz_xxxxxx = pbuffer.data(idx_dip_ii + 224);

    auto tr_x_xxxyzz_xxxxxy = pbuffer.data(idx_dip_ii + 225);

    auto tr_x_xxxyzz_xxxxxz = pbuffer.data(idx_dip_ii + 226);

    auto tr_x_xxxyzz_xxxxyy = pbuffer.data(idx_dip_ii + 227);

    auto tr_x_xxxyzz_xxxxyz = pbuffer.data(idx_dip_ii + 228);

    auto tr_x_xxxyzz_xxxxzz = pbuffer.data(idx_dip_ii + 229);

    auto tr_x_xxxyzz_xxxyyy = pbuffer.data(idx_dip_ii + 230);

    auto tr_x_xxxyzz_xxxyyz = pbuffer.data(idx_dip_ii + 231);

    auto tr_x_xxxyzz_xxxyzz = pbuffer.data(idx_dip_ii + 232);

    auto tr_x_xxxyzz_xxxzzz = pbuffer.data(idx_dip_ii + 233);

    auto tr_x_xxxyzz_xxyyyy = pbuffer.data(idx_dip_ii + 234);

    auto tr_x_xxxyzz_xxyyyz = pbuffer.data(idx_dip_ii + 235);

    auto tr_x_xxxyzz_xxyyzz = pbuffer.data(idx_dip_ii + 236);

    auto tr_x_xxxyzz_xxyzzz = pbuffer.data(idx_dip_ii + 237);

    auto tr_x_xxxyzz_xxzzzz = pbuffer.data(idx_dip_ii + 238);

    auto tr_x_xxxyzz_xyyyyy = pbuffer.data(idx_dip_ii + 239);

    auto tr_x_xxxyzz_xyyyyz = pbuffer.data(idx_dip_ii + 240);

    auto tr_x_xxxyzz_xyyyzz = pbuffer.data(idx_dip_ii + 241);

    auto tr_x_xxxyzz_xyyzzz = pbuffer.data(idx_dip_ii + 242);

    auto tr_x_xxxyzz_xyzzzz = pbuffer.data(idx_dip_ii + 243);

    auto tr_x_xxxyzz_xzzzzz = pbuffer.data(idx_dip_ii + 244);

    auto tr_x_xxxyzz_yyyyyy = pbuffer.data(idx_dip_ii + 245);

    auto tr_x_xxxyzz_yyyyyz = pbuffer.data(idx_dip_ii + 246);

    auto tr_x_xxxyzz_yyyyzz = pbuffer.data(idx_dip_ii + 247);

    auto tr_x_xxxyzz_yyyzzz = pbuffer.data(idx_dip_ii + 248);

    auto tr_x_xxxyzz_yyzzzz = pbuffer.data(idx_dip_ii + 249);

    auto tr_x_xxxyzz_yzzzzz = pbuffer.data(idx_dip_ii + 250);

    auto tr_x_xxxyzz_zzzzzz = pbuffer.data(idx_dip_ii + 251);

    #pragma omp simd aligned(pa_y, tr_x_xxxyzz_xxxxxx, tr_x_xxxyzz_xxxxxy, tr_x_xxxyzz_xxxxxz, tr_x_xxxyzz_xxxxyy, tr_x_xxxyzz_xxxxyz, tr_x_xxxyzz_xxxxzz, tr_x_xxxyzz_xxxyyy, tr_x_xxxyzz_xxxyyz, tr_x_xxxyzz_xxxyzz, tr_x_xxxyzz_xxxzzz, tr_x_xxxyzz_xxyyyy, tr_x_xxxyzz_xxyyyz, tr_x_xxxyzz_xxyyzz, tr_x_xxxyzz_xxyzzz, tr_x_xxxyzz_xxzzzz, tr_x_xxxyzz_xyyyyy, tr_x_xxxyzz_xyyyyz, tr_x_xxxyzz_xyyyzz, tr_x_xxxyzz_xyyzzz, tr_x_xxxyzz_xyzzzz, tr_x_xxxyzz_xzzzzz, tr_x_xxxyzz_yyyyyy, tr_x_xxxyzz_yyyyyz, tr_x_xxxyzz_yyyyzz, tr_x_xxxyzz_yyyzzz, tr_x_xxxyzz_yyzzzz, tr_x_xxxyzz_yzzzzz, tr_x_xxxyzz_zzzzzz, tr_x_xxxzz_xxxxx, tr_x_xxxzz_xxxxxx, tr_x_xxxzz_xxxxxy, tr_x_xxxzz_xxxxxz, tr_x_xxxzz_xxxxy, tr_x_xxxzz_xxxxyy, tr_x_xxxzz_xxxxyz, tr_x_xxxzz_xxxxz, tr_x_xxxzz_xxxxzz, tr_x_xxxzz_xxxyy, tr_x_xxxzz_xxxyyy, tr_x_xxxzz_xxxyyz, tr_x_xxxzz_xxxyz, tr_x_xxxzz_xxxyzz, tr_x_xxxzz_xxxzz, tr_x_xxxzz_xxxzzz, tr_x_xxxzz_xxyyy, tr_x_xxxzz_xxyyyy, tr_x_xxxzz_xxyyyz, tr_x_xxxzz_xxyyz, tr_x_xxxzz_xxyyzz, tr_x_xxxzz_xxyzz, tr_x_xxxzz_xxyzzz, tr_x_xxxzz_xxzzz, tr_x_xxxzz_xxzzzz, tr_x_xxxzz_xyyyy, tr_x_xxxzz_xyyyyy, tr_x_xxxzz_xyyyyz, tr_x_xxxzz_xyyyz, tr_x_xxxzz_xyyyzz, tr_x_xxxzz_xyyzz, tr_x_xxxzz_xyyzzz, tr_x_xxxzz_xyzzz, tr_x_xxxzz_xyzzzz, tr_x_xxxzz_xzzzz, tr_x_xxxzz_xzzzzz, tr_x_xxxzz_yyyyy, tr_x_xxxzz_yyyyyy, tr_x_xxxzz_yyyyyz, tr_x_xxxzz_yyyyz, tr_x_xxxzz_yyyyzz, tr_x_xxxzz_yyyzz, tr_x_xxxzz_yyyzzz, tr_x_xxxzz_yyzzz, tr_x_xxxzz_yyzzzz, tr_x_xxxzz_yzzzz, tr_x_xxxzz_yzzzzz, tr_x_xxxzz_zzzzz, tr_x_xxxzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyzz_xxxxxx[i] = tr_x_xxxzz_xxxxxx[i] * pa_y[i];

        tr_x_xxxyzz_xxxxxy[i] = tr_x_xxxzz_xxxxx[i] * fe_0 + tr_x_xxxzz_xxxxxy[i] * pa_y[i];

        tr_x_xxxyzz_xxxxxz[i] = tr_x_xxxzz_xxxxxz[i] * pa_y[i];

        tr_x_xxxyzz_xxxxyy[i] = 2.0 * tr_x_xxxzz_xxxxy[i] * fe_0 + tr_x_xxxzz_xxxxyy[i] * pa_y[i];

        tr_x_xxxyzz_xxxxyz[i] = tr_x_xxxzz_xxxxz[i] * fe_0 + tr_x_xxxzz_xxxxyz[i] * pa_y[i];

        tr_x_xxxyzz_xxxxzz[i] = tr_x_xxxzz_xxxxzz[i] * pa_y[i];

        tr_x_xxxyzz_xxxyyy[i] = 3.0 * tr_x_xxxzz_xxxyy[i] * fe_0 + tr_x_xxxzz_xxxyyy[i] * pa_y[i];

        tr_x_xxxyzz_xxxyyz[i] = 2.0 * tr_x_xxxzz_xxxyz[i] * fe_0 + tr_x_xxxzz_xxxyyz[i] * pa_y[i];

        tr_x_xxxyzz_xxxyzz[i] = tr_x_xxxzz_xxxzz[i] * fe_0 + tr_x_xxxzz_xxxyzz[i] * pa_y[i];

        tr_x_xxxyzz_xxxzzz[i] = tr_x_xxxzz_xxxzzz[i] * pa_y[i];

        tr_x_xxxyzz_xxyyyy[i] = 4.0 * tr_x_xxxzz_xxyyy[i] * fe_0 + tr_x_xxxzz_xxyyyy[i] * pa_y[i];

        tr_x_xxxyzz_xxyyyz[i] = 3.0 * tr_x_xxxzz_xxyyz[i] * fe_0 + tr_x_xxxzz_xxyyyz[i] * pa_y[i];

        tr_x_xxxyzz_xxyyzz[i] = 2.0 * tr_x_xxxzz_xxyzz[i] * fe_0 + tr_x_xxxzz_xxyyzz[i] * pa_y[i];

        tr_x_xxxyzz_xxyzzz[i] = tr_x_xxxzz_xxzzz[i] * fe_0 + tr_x_xxxzz_xxyzzz[i] * pa_y[i];

        tr_x_xxxyzz_xxzzzz[i] = tr_x_xxxzz_xxzzzz[i] * pa_y[i];

        tr_x_xxxyzz_xyyyyy[i] = 5.0 * tr_x_xxxzz_xyyyy[i] * fe_0 + tr_x_xxxzz_xyyyyy[i] * pa_y[i];

        tr_x_xxxyzz_xyyyyz[i] = 4.0 * tr_x_xxxzz_xyyyz[i] * fe_0 + tr_x_xxxzz_xyyyyz[i] * pa_y[i];

        tr_x_xxxyzz_xyyyzz[i] = 3.0 * tr_x_xxxzz_xyyzz[i] * fe_0 + tr_x_xxxzz_xyyyzz[i] * pa_y[i];

        tr_x_xxxyzz_xyyzzz[i] = 2.0 * tr_x_xxxzz_xyzzz[i] * fe_0 + tr_x_xxxzz_xyyzzz[i] * pa_y[i];

        tr_x_xxxyzz_xyzzzz[i] = tr_x_xxxzz_xzzzz[i] * fe_0 + tr_x_xxxzz_xyzzzz[i] * pa_y[i];

        tr_x_xxxyzz_xzzzzz[i] = tr_x_xxxzz_xzzzzz[i] * pa_y[i];

        tr_x_xxxyzz_yyyyyy[i] = 6.0 * tr_x_xxxzz_yyyyy[i] * fe_0 + tr_x_xxxzz_yyyyyy[i] * pa_y[i];

        tr_x_xxxyzz_yyyyyz[i] = 5.0 * tr_x_xxxzz_yyyyz[i] * fe_0 + tr_x_xxxzz_yyyyyz[i] * pa_y[i];

        tr_x_xxxyzz_yyyyzz[i] = 4.0 * tr_x_xxxzz_yyyzz[i] * fe_0 + tr_x_xxxzz_yyyyzz[i] * pa_y[i];

        tr_x_xxxyzz_yyyzzz[i] = 3.0 * tr_x_xxxzz_yyzzz[i] * fe_0 + tr_x_xxxzz_yyyzzz[i] * pa_y[i];

        tr_x_xxxyzz_yyzzzz[i] = 2.0 * tr_x_xxxzz_yzzzz[i] * fe_0 + tr_x_xxxzz_yyzzzz[i] * pa_y[i];

        tr_x_xxxyzz_yzzzzz[i] = tr_x_xxxzz_zzzzz[i] * fe_0 + tr_x_xxxzz_yzzzzz[i] * pa_y[i];

        tr_x_xxxyzz_zzzzzz[i] = tr_x_xxxzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : II

    auto tr_x_xxxzzz_xxxxxx = pbuffer.data(idx_dip_ii + 252);

    auto tr_x_xxxzzz_xxxxxy = pbuffer.data(idx_dip_ii + 253);

    auto tr_x_xxxzzz_xxxxxz = pbuffer.data(idx_dip_ii + 254);

    auto tr_x_xxxzzz_xxxxyy = pbuffer.data(idx_dip_ii + 255);

    auto tr_x_xxxzzz_xxxxyz = pbuffer.data(idx_dip_ii + 256);

    auto tr_x_xxxzzz_xxxxzz = pbuffer.data(idx_dip_ii + 257);

    auto tr_x_xxxzzz_xxxyyy = pbuffer.data(idx_dip_ii + 258);

    auto tr_x_xxxzzz_xxxyyz = pbuffer.data(idx_dip_ii + 259);

    auto tr_x_xxxzzz_xxxyzz = pbuffer.data(idx_dip_ii + 260);

    auto tr_x_xxxzzz_xxxzzz = pbuffer.data(idx_dip_ii + 261);

    auto tr_x_xxxzzz_xxyyyy = pbuffer.data(idx_dip_ii + 262);

    auto tr_x_xxxzzz_xxyyyz = pbuffer.data(idx_dip_ii + 263);

    auto tr_x_xxxzzz_xxyyzz = pbuffer.data(idx_dip_ii + 264);

    auto tr_x_xxxzzz_xxyzzz = pbuffer.data(idx_dip_ii + 265);

    auto tr_x_xxxzzz_xxzzzz = pbuffer.data(idx_dip_ii + 266);

    auto tr_x_xxxzzz_xyyyyy = pbuffer.data(idx_dip_ii + 267);

    auto tr_x_xxxzzz_xyyyyz = pbuffer.data(idx_dip_ii + 268);

    auto tr_x_xxxzzz_xyyyzz = pbuffer.data(idx_dip_ii + 269);

    auto tr_x_xxxzzz_xyyzzz = pbuffer.data(idx_dip_ii + 270);

    auto tr_x_xxxzzz_xyzzzz = pbuffer.data(idx_dip_ii + 271);

    auto tr_x_xxxzzz_xzzzzz = pbuffer.data(idx_dip_ii + 272);

    auto tr_x_xxxzzz_yyyyyy = pbuffer.data(idx_dip_ii + 273);

    auto tr_x_xxxzzz_yyyyyz = pbuffer.data(idx_dip_ii + 274);

    auto tr_x_xxxzzz_yyyyzz = pbuffer.data(idx_dip_ii + 275);

    auto tr_x_xxxzzz_yyyzzz = pbuffer.data(idx_dip_ii + 276);

    auto tr_x_xxxzzz_yyzzzz = pbuffer.data(idx_dip_ii + 277);

    auto tr_x_xxxzzz_yzzzzz = pbuffer.data(idx_dip_ii + 278);

    auto tr_x_xxxzzz_zzzzzz = pbuffer.data(idx_dip_ii + 279);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxxz_xxxxxx, tr_x_xxxz_xxxxxy, tr_x_xxxz_xxxxxz, tr_x_xxxz_xxxxyy, tr_x_xxxz_xxxxyz, tr_x_xxxz_xxxxzz, tr_x_xxxz_xxxyyy, tr_x_xxxz_xxxyyz, tr_x_xxxz_xxxyzz, tr_x_xxxz_xxxzzz, tr_x_xxxz_xxyyyy, tr_x_xxxz_xxyyyz, tr_x_xxxz_xxyyzz, tr_x_xxxz_xxyzzz, tr_x_xxxz_xxzzzz, tr_x_xxxz_xyyyyy, tr_x_xxxz_xyyyyz, tr_x_xxxz_xyyyzz, tr_x_xxxz_xyyzzz, tr_x_xxxz_xyzzzz, tr_x_xxxz_xzzzzz, tr_x_xxxz_yyyyyy, tr_x_xxxzz_xxxxx, tr_x_xxxzz_xxxxxx, tr_x_xxxzz_xxxxxy, tr_x_xxxzz_xxxxxz, tr_x_xxxzz_xxxxy, tr_x_xxxzz_xxxxyy, tr_x_xxxzz_xxxxyz, tr_x_xxxzz_xxxxz, tr_x_xxxzz_xxxxzz, tr_x_xxxzz_xxxyy, tr_x_xxxzz_xxxyyy, tr_x_xxxzz_xxxyyz, tr_x_xxxzz_xxxyz, tr_x_xxxzz_xxxyzz, tr_x_xxxzz_xxxzz, tr_x_xxxzz_xxxzzz, tr_x_xxxzz_xxyyy, tr_x_xxxzz_xxyyyy, tr_x_xxxzz_xxyyyz, tr_x_xxxzz_xxyyz, tr_x_xxxzz_xxyyzz, tr_x_xxxzz_xxyzz, tr_x_xxxzz_xxyzzz, tr_x_xxxzz_xxzzz, tr_x_xxxzz_xxzzzz, tr_x_xxxzz_xyyyy, tr_x_xxxzz_xyyyyy, tr_x_xxxzz_xyyyyz, tr_x_xxxzz_xyyyz, tr_x_xxxzz_xyyyzz, tr_x_xxxzz_xyyzz, tr_x_xxxzz_xyyzzz, tr_x_xxxzz_xyzzz, tr_x_xxxzz_xyzzzz, tr_x_xxxzz_xzzzz, tr_x_xxxzz_xzzzzz, tr_x_xxxzz_yyyyyy, tr_x_xxxzzz_xxxxxx, tr_x_xxxzzz_xxxxxy, tr_x_xxxzzz_xxxxxz, tr_x_xxxzzz_xxxxyy, tr_x_xxxzzz_xxxxyz, tr_x_xxxzzz_xxxxzz, tr_x_xxxzzz_xxxyyy, tr_x_xxxzzz_xxxyyz, tr_x_xxxzzz_xxxyzz, tr_x_xxxzzz_xxxzzz, tr_x_xxxzzz_xxyyyy, tr_x_xxxzzz_xxyyyz, tr_x_xxxzzz_xxyyzz, tr_x_xxxzzz_xxyzzz, tr_x_xxxzzz_xxzzzz, tr_x_xxxzzz_xyyyyy, tr_x_xxxzzz_xyyyyz, tr_x_xxxzzz_xyyyzz, tr_x_xxxzzz_xyyzzz, tr_x_xxxzzz_xyzzzz, tr_x_xxxzzz_xzzzzz, tr_x_xxxzzz_yyyyyy, tr_x_xxxzzz_yyyyyz, tr_x_xxxzzz_yyyyzz, tr_x_xxxzzz_yyyzzz, tr_x_xxxzzz_yyzzzz, tr_x_xxxzzz_yzzzzz, tr_x_xxxzzz_zzzzzz, tr_x_xxzzz_yyyyyz, tr_x_xxzzz_yyyyzz, tr_x_xxzzz_yyyzzz, tr_x_xxzzz_yyzzzz, tr_x_xxzzz_yzzzzz, tr_x_xxzzz_zzzzzz, tr_x_xzzz_yyyyyz, tr_x_xzzz_yyyyzz, tr_x_xzzz_yyyzzz, tr_x_xzzz_yyzzzz, tr_x_xzzz_yzzzzz, tr_x_xzzz_zzzzzz, ts_xxzzz_yyyyyz, ts_xxzzz_yyyyzz, ts_xxzzz_yyyzzz, ts_xxzzz_yyzzzz, ts_xxzzz_yzzzzz, ts_xxzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzzz_xxxxxx[i] = 2.0 * tr_x_xxxz_xxxxxx[i] * fe_0 + tr_x_xxxzz_xxxxxx[i] * pa_z[i];

        tr_x_xxxzzz_xxxxxy[i] = 2.0 * tr_x_xxxz_xxxxxy[i] * fe_0 + tr_x_xxxzz_xxxxxy[i] * pa_z[i];

        tr_x_xxxzzz_xxxxxz[i] = 2.0 * tr_x_xxxz_xxxxxz[i] * fe_0 + tr_x_xxxzz_xxxxx[i] * fe_0 + tr_x_xxxzz_xxxxxz[i] * pa_z[i];

        tr_x_xxxzzz_xxxxyy[i] = 2.0 * tr_x_xxxz_xxxxyy[i] * fe_0 + tr_x_xxxzz_xxxxyy[i] * pa_z[i];

        tr_x_xxxzzz_xxxxyz[i] = 2.0 * tr_x_xxxz_xxxxyz[i] * fe_0 + tr_x_xxxzz_xxxxy[i] * fe_0 + tr_x_xxxzz_xxxxyz[i] * pa_z[i];

        tr_x_xxxzzz_xxxxzz[i] = 2.0 * tr_x_xxxz_xxxxzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xxxxz[i] * fe_0 + tr_x_xxxzz_xxxxzz[i] * pa_z[i];

        tr_x_xxxzzz_xxxyyy[i] = 2.0 * tr_x_xxxz_xxxyyy[i] * fe_0 + tr_x_xxxzz_xxxyyy[i] * pa_z[i];

        tr_x_xxxzzz_xxxyyz[i] = 2.0 * tr_x_xxxz_xxxyyz[i] * fe_0 + tr_x_xxxzz_xxxyy[i] * fe_0 + tr_x_xxxzz_xxxyyz[i] * pa_z[i];

        tr_x_xxxzzz_xxxyzz[i] = 2.0 * tr_x_xxxz_xxxyzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xxxyz[i] * fe_0 + tr_x_xxxzz_xxxyzz[i] * pa_z[i];

        tr_x_xxxzzz_xxxzzz[i] = 2.0 * tr_x_xxxz_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxxzz_xxxzz[i] * fe_0 + tr_x_xxxzz_xxxzzz[i] * pa_z[i];

        tr_x_xxxzzz_xxyyyy[i] = 2.0 * tr_x_xxxz_xxyyyy[i] * fe_0 + tr_x_xxxzz_xxyyyy[i] * pa_z[i];

        tr_x_xxxzzz_xxyyyz[i] = 2.0 * tr_x_xxxz_xxyyyz[i] * fe_0 + tr_x_xxxzz_xxyyy[i] * fe_0 + tr_x_xxxzz_xxyyyz[i] * pa_z[i];

        tr_x_xxxzzz_xxyyzz[i] = 2.0 * tr_x_xxxz_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xxyyz[i] * fe_0 + tr_x_xxxzz_xxyyzz[i] * pa_z[i];

        tr_x_xxxzzz_xxyzzz[i] = 2.0 * tr_x_xxxz_xxyzzz[i] * fe_0 + 3.0 * tr_x_xxxzz_xxyzz[i] * fe_0 + tr_x_xxxzz_xxyzzz[i] * pa_z[i];

        tr_x_xxxzzz_xxzzzz[i] = 2.0 * tr_x_xxxz_xxzzzz[i] * fe_0 + 4.0 * tr_x_xxxzz_xxzzz[i] * fe_0 + tr_x_xxxzz_xxzzzz[i] * pa_z[i];

        tr_x_xxxzzz_xyyyyy[i] = 2.0 * tr_x_xxxz_xyyyyy[i] * fe_0 + tr_x_xxxzz_xyyyyy[i] * pa_z[i];

        tr_x_xxxzzz_xyyyyz[i] = 2.0 * tr_x_xxxz_xyyyyz[i] * fe_0 + tr_x_xxxzz_xyyyy[i] * fe_0 + tr_x_xxxzz_xyyyyz[i] * pa_z[i];

        tr_x_xxxzzz_xyyyzz[i] = 2.0 * tr_x_xxxz_xyyyzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xyyyz[i] * fe_0 + tr_x_xxxzz_xyyyzz[i] * pa_z[i];

        tr_x_xxxzzz_xyyzzz[i] = 2.0 * tr_x_xxxz_xyyzzz[i] * fe_0 + 3.0 * tr_x_xxxzz_xyyzz[i] * fe_0 + tr_x_xxxzz_xyyzzz[i] * pa_z[i];

        tr_x_xxxzzz_xyzzzz[i] = 2.0 * tr_x_xxxz_xyzzzz[i] * fe_0 + 4.0 * tr_x_xxxzz_xyzzz[i] * fe_0 + tr_x_xxxzz_xyzzzz[i] * pa_z[i];

        tr_x_xxxzzz_xzzzzz[i] = 2.0 * tr_x_xxxz_xzzzzz[i] * fe_0 + 5.0 * tr_x_xxxzz_xzzzz[i] * fe_0 + tr_x_xxxzz_xzzzzz[i] * pa_z[i];

        tr_x_xxxzzz_yyyyyy[i] = 2.0 * tr_x_xxxz_yyyyyy[i] * fe_0 + tr_x_xxxzz_yyyyyy[i] * pa_z[i];

        tr_x_xxxzzz_yyyyyz[i] = 2.0 * tr_x_xzzz_yyyyyz[i] * fe_0 + ts_xxzzz_yyyyyz[i] * fe_0 + tr_x_xxzzz_yyyyyz[i] * pa_x[i];

        tr_x_xxxzzz_yyyyzz[i] = 2.0 * tr_x_xzzz_yyyyzz[i] * fe_0 + ts_xxzzz_yyyyzz[i] * fe_0 + tr_x_xxzzz_yyyyzz[i] * pa_x[i];

        tr_x_xxxzzz_yyyzzz[i] = 2.0 * tr_x_xzzz_yyyzzz[i] * fe_0 + ts_xxzzz_yyyzzz[i] * fe_0 + tr_x_xxzzz_yyyzzz[i] * pa_x[i];

        tr_x_xxxzzz_yyzzzz[i] = 2.0 * tr_x_xzzz_yyzzzz[i] * fe_0 + ts_xxzzz_yyzzzz[i] * fe_0 + tr_x_xxzzz_yyzzzz[i] * pa_x[i];

        tr_x_xxxzzz_yzzzzz[i] = 2.0 * tr_x_xzzz_yzzzzz[i] * fe_0 + ts_xxzzz_yzzzzz[i] * fe_0 + tr_x_xxzzz_yzzzzz[i] * pa_x[i];

        tr_x_xxxzzz_zzzzzz[i] = 2.0 * tr_x_xzzz_zzzzzz[i] * fe_0 + ts_xxzzz_zzzzzz[i] * fe_0 + tr_x_xxzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : II

    auto tr_x_xxyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 280);

    auto tr_x_xxyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 281);

    auto tr_x_xxyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 282);

    auto tr_x_xxyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 283);

    auto tr_x_xxyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 284);

    auto tr_x_xxyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 285);

    auto tr_x_xxyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 286);

    auto tr_x_xxyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 287);

    auto tr_x_xxyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 288);

    auto tr_x_xxyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 289);

    auto tr_x_xxyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 290);

    auto tr_x_xxyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 291);

    auto tr_x_xxyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 292);

    auto tr_x_xxyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 293);

    auto tr_x_xxyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 294);

    auto tr_x_xxyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 295);

    auto tr_x_xxyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 296);

    auto tr_x_xxyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 297);

    auto tr_x_xxyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 298);

    auto tr_x_xxyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 299);

    auto tr_x_xxyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 300);

    auto tr_x_xxyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 301);

    auto tr_x_xxyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 302);

    auto tr_x_xxyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 303);

    auto tr_x_xxyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 304);

    auto tr_x_xxyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 305);

    auto tr_x_xxyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 306);

    auto tr_x_xxyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 307);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxyy_xxxxxx, tr_x_xxyy_xxxxxy, tr_x_xxyy_xxxxxz, tr_x_xxyy_xxxxyy, tr_x_xxyy_xxxxyz, tr_x_xxyy_xxxxzz, tr_x_xxyy_xxxyyy, tr_x_xxyy_xxxyyz, tr_x_xxyy_xxxyzz, tr_x_xxyy_xxxzzz, tr_x_xxyy_xxyyyy, tr_x_xxyy_xxyyyz, tr_x_xxyy_xxyyzz, tr_x_xxyy_xxyzzz, tr_x_xxyy_xxzzzz, tr_x_xxyy_xyyyyy, tr_x_xxyy_xyyyyz, tr_x_xxyy_xyyyzz, tr_x_xxyy_xyyzzz, tr_x_xxyy_xyzzzz, tr_x_xxyy_xzzzzz, tr_x_xxyy_zzzzzz, tr_x_xxyyy_xxxxx, tr_x_xxyyy_xxxxxx, tr_x_xxyyy_xxxxxy, tr_x_xxyyy_xxxxxz, tr_x_xxyyy_xxxxy, tr_x_xxyyy_xxxxyy, tr_x_xxyyy_xxxxyz, tr_x_xxyyy_xxxxz, tr_x_xxyyy_xxxxzz, tr_x_xxyyy_xxxyy, tr_x_xxyyy_xxxyyy, tr_x_xxyyy_xxxyyz, tr_x_xxyyy_xxxyz, tr_x_xxyyy_xxxyzz, tr_x_xxyyy_xxxzz, tr_x_xxyyy_xxxzzz, tr_x_xxyyy_xxyyy, tr_x_xxyyy_xxyyyy, tr_x_xxyyy_xxyyyz, tr_x_xxyyy_xxyyz, tr_x_xxyyy_xxyyzz, tr_x_xxyyy_xxyzz, tr_x_xxyyy_xxyzzz, tr_x_xxyyy_xxzzz, tr_x_xxyyy_xxzzzz, tr_x_xxyyy_xyyyy, tr_x_xxyyy_xyyyyy, tr_x_xxyyy_xyyyyz, tr_x_xxyyy_xyyyz, tr_x_xxyyy_xyyyzz, tr_x_xxyyy_xyyzz, tr_x_xxyyy_xyyzzz, tr_x_xxyyy_xyzzz, tr_x_xxyyy_xyzzzz, tr_x_xxyyy_xzzzz, tr_x_xxyyy_xzzzzz, tr_x_xxyyy_zzzzzz, tr_x_xxyyyy_xxxxxx, tr_x_xxyyyy_xxxxxy, tr_x_xxyyyy_xxxxxz, tr_x_xxyyyy_xxxxyy, tr_x_xxyyyy_xxxxyz, tr_x_xxyyyy_xxxxzz, tr_x_xxyyyy_xxxyyy, tr_x_xxyyyy_xxxyyz, tr_x_xxyyyy_xxxyzz, tr_x_xxyyyy_xxxzzz, tr_x_xxyyyy_xxyyyy, tr_x_xxyyyy_xxyyyz, tr_x_xxyyyy_xxyyzz, tr_x_xxyyyy_xxyzzz, tr_x_xxyyyy_xxzzzz, tr_x_xxyyyy_xyyyyy, tr_x_xxyyyy_xyyyyz, tr_x_xxyyyy_xyyyzz, tr_x_xxyyyy_xyyzzz, tr_x_xxyyyy_xyzzzz, tr_x_xxyyyy_xzzzzz, tr_x_xxyyyy_yyyyyy, tr_x_xxyyyy_yyyyyz, tr_x_xxyyyy_yyyyzz, tr_x_xxyyyy_yyyzzz, tr_x_xxyyyy_yyzzzz, tr_x_xxyyyy_yzzzzz, tr_x_xxyyyy_zzzzzz, tr_x_xyyyy_yyyyyy, tr_x_xyyyy_yyyyyz, tr_x_xyyyy_yyyyzz, tr_x_xyyyy_yyyzzz, tr_x_xyyyy_yyzzzz, tr_x_xyyyy_yzzzzz, tr_x_yyyy_yyyyyy, tr_x_yyyy_yyyyyz, tr_x_yyyy_yyyyzz, tr_x_yyyy_yyyzzz, tr_x_yyyy_yyzzzz, tr_x_yyyy_yzzzzz, ts_xyyyy_yyyyyy, ts_xyyyy_yyyyyz, ts_xyyyy_yyyyzz, ts_xyyyy_yyyzzz, ts_xyyyy_yyzzzz, ts_xyyyy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyy_xxxxxx[i] = 3.0 * tr_x_xxyy_xxxxxx[i] * fe_0 + tr_x_xxyyy_xxxxxx[i] * pa_y[i];

        tr_x_xxyyyy_xxxxxy[i] = 3.0 * tr_x_xxyy_xxxxxy[i] * fe_0 + tr_x_xxyyy_xxxxx[i] * fe_0 + tr_x_xxyyy_xxxxxy[i] * pa_y[i];

        tr_x_xxyyyy_xxxxxz[i] = 3.0 * tr_x_xxyy_xxxxxz[i] * fe_0 + tr_x_xxyyy_xxxxxz[i] * pa_y[i];

        tr_x_xxyyyy_xxxxyy[i] = 3.0 * tr_x_xxyy_xxxxyy[i] * fe_0 + 2.0 * tr_x_xxyyy_xxxxy[i] * fe_0 + tr_x_xxyyy_xxxxyy[i] * pa_y[i];

        tr_x_xxyyyy_xxxxyz[i] = 3.0 * tr_x_xxyy_xxxxyz[i] * fe_0 + tr_x_xxyyy_xxxxz[i] * fe_0 + tr_x_xxyyy_xxxxyz[i] * pa_y[i];

        tr_x_xxyyyy_xxxxzz[i] = 3.0 * tr_x_xxyy_xxxxzz[i] * fe_0 + tr_x_xxyyy_xxxxzz[i] * pa_y[i];

        tr_x_xxyyyy_xxxyyy[i] = 3.0 * tr_x_xxyy_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxyyy_xxxyy[i] * fe_0 + tr_x_xxyyy_xxxyyy[i] * pa_y[i];

        tr_x_xxyyyy_xxxyyz[i] = 3.0 * tr_x_xxyy_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxyyy_xxxyz[i] * fe_0 + tr_x_xxyyy_xxxyyz[i] * pa_y[i];

        tr_x_xxyyyy_xxxyzz[i] = 3.0 * tr_x_xxyy_xxxyzz[i] * fe_0 + tr_x_xxyyy_xxxzz[i] * fe_0 + tr_x_xxyyy_xxxyzz[i] * pa_y[i];

        tr_x_xxyyyy_xxxzzz[i] = 3.0 * tr_x_xxyy_xxxzzz[i] * fe_0 + tr_x_xxyyy_xxxzzz[i] * pa_y[i];

        tr_x_xxyyyy_xxyyyy[i] = 3.0 * tr_x_xxyy_xxyyyy[i] * fe_0 + 4.0 * tr_x_xxyyy_xxyyy[i] * fe_0 + tr_x_xxyyy_xxyyyy[i] * pa_y[i];

        tr_x_xxyyyy_xxyyyz[i] = 3.0 * tr_x_xxyy_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxyyy_xxyyz[i] * fe_0 + tr_x_xxyyy_xxyyyz[i] * pa_y[i];

        tr_x_xxyyyy_xxyyzz[i] = 3.0 * tr_x_xxyy_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxyyy_xxyzz[i] * fe_0 + tr_x_xxyyy_xxyyzz[i] * pa_y[i];

        tr_x_xxyyyy_xxyzzz[i] = 3.0 * tr_x_xxyy_xxyzzz[i] * fe_0 + tr_x_xxyyy_xxzzz[i] * fe_0 + tr_x_xxyyy_xxyzzz[i] * pa_y[i];

        tr_x_xxyyyy_xxzzzz[i] = 3.0 * tr_x_xxyy_xxzzzz[i] * fe_0 + tr_x_xxyyy_xxzzzz[i] * pa_y[i];

        tr_x_xxyyyy_xyyyyy[i] = 3.0 * tr_x_xxyy_xyyyyy[i] * fe_0 + 5.0 * tr_x_xxyyy_xyyyy[i] * fe_0 + tr_x_xxyyy_xyyyyy[i] * pa_y[i];

        tr_x_xxyyyy_xyyyyz[i] = 3.0 * tr_x_xxyy_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxyyy_xyyyz[i] * fe_0 + tr_x_xxyyy_xyyyyz[i] * pa_y[i];

        tr_x_xxyyyy_xyyyzz[i] = 3.0 * tr_x_xxyy_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxyyy_xyyzz[i] * fe_0 + tr_x_xxyyy_xyyyzz[i] * pa_y[i];

        tr_x_xxyyyy_xyyzzz[i] = 3.0 * tr_x_xxyy_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxyyy_xyzzz[i] * fe_0 + tr_x_xxyyy_xyyzzz[i] * pa_y[i];

        tr_x_xxyyyy_xyzzzz[i] = 3.0 * tr_x_xxyy_xyzzzz[i] * fe_0 + tr_x_xxyyy_xzzzz[i] * fe_0 + tr_x_xxyyy_xyzzzz[i] * pa_y[i];

        tr_x_xxyyyy_xzzzzz[i] = 3.0 * tr_x_xxyy_xzzzzz[i] * fe_0 + tr_x_xxyyy_xzzzzz[i] * pa_y[i];

        tr_x_xxyyyy_yyyyyy[i] = tr_x_yyyy_yyyyyy[i] * fe_0 + ts_xyyyy_yyyyyy[i] * fe_0 + tr_x_xyyyy_yyyyyy[i] * pa_x[i];

        tr_x_xxyyyy_yyyyyz[i] = tr_x_yyyy_yyyyyz[i] * fe_0 + ts_xyyyy_yyyyyz[i] * fe_0 + tr_x_xyyyy_yyyyyz[i] * pa_x[i];

        tr_x_xxyyyy_yyyyzz[i] = tr_x_yyyy_yyyyzz[i] * fe_0 + ts_xyyyy_yyyyzz[i] * fe_0 + tr_x_xyyyy_yyyyzz[i] * pa_x[i];

        tr_x_xxyyyy_yyyzzz[i] = tr_x_yyyy_yyyzzz[i] * fe_0 + ts_xyyyy_yyyzzz[i] * fe_0 + tr_x_xyyyy_yyyzzz[i] * pa_x[i];

        tr_x_xxyyyy_yyzzzz[i] = tr_x_yyyy_yyzzzz[i] * fe_0 + ts_xyyyy_yyzzzz[i] * fe_0 + tr_x_xyyyy_yyzzzz[i] * pa_x[i];

        tr_x_xxyyyy_yzzzzz[i] = tr_x_yyyy_yzzzzz[i] * fe_0 + ts_xyyyy_yzzzzz[i] * fe_0 + tr_x_xyyyy_yzzzzz[i] * pa_x[i];

        tr_x_xxyyyy_zzzzzz[i] = 3.0 * tr_x_xxyy_zzzzzz[i] * fe_0 + tr_x_xxyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 308-336 components of targeted buffer : II

    auto tr_x_xxyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 308);

    auto tr_x_xxyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 309);

    auto tr_x_xxyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 310);

    auto tr_x_xxyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 311);

    auto tr_x_xxyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 312);

    auto tr_x_xxyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 313);

    auto tr_x_xxyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 314);

    auto tr_x_xxyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 315);

    auto tr_x_xxyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 316);

    auto tr_x_xxyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 317);

    auto tr_x_xxyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 318);

    auto tr_x_xxyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 319);

    auto tr_x_xxyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 320);

    auto tr_x_xxyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 321);

    auto tr_x_xxyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 322);

    auto tr_x_xxyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 323);

    auto tr_x_xxyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 324);

    auto tr_x_xxyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 325);

    auto tr_x_xxyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 326);

    auto tr_x_xxyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 327);

    auto tr_x_xxyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 328);

    auto tr_x_xxyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 329);

    auto tr_x_xxyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 330);

    auto tr_x_xxyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 331);

    auto tr_x_xxyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 332);

    auto tr_x_xxyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 333);

    auto tr_x_xxyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 334);

    auto tr_x_xxyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 335);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxyyy_xxxxxx, tr_x_xxyyy_xxxxxy, tr_x_xxyyy_xxxxy, tr_x_xxyyy_xxxxyy, tr_x_xxyyy_xxxxyz, tr_x_xxyyy_xxxyy, tr_x_xxyyy_xxxyyy, tr_x_xxyyy_xxxyyz, tr_x_xxyyy_xxxyz, tr_x_xxyyy_xxxyzz, tr_x_xxyyy_xxyyy, tr_x_xxyyy_xxyyyy, tr_x_xxyyy_xxyyyz, tr_x_xxyyy_xxyyz, tr_x_xxyyy_xxyyzz, tr_x_xxyyy_xxyzz, tr_x_xxyyy_xxyzzz, tr_x_xxyyy_xyyyy, tr_x_xxyyy_xyyyyy, tr_x_xxyyy_xyyyyz, tr_x_xxyyy_xyyyz, tr_x_xxyyy_xyyyzz, tr_x_xxyyy_xyyzz, tr_x_xxyyy_xyyzzz, tr_x_xxyyy_xyzzz, tr_x_xxyyy_xyzzzz, tr_x_xxyyy_yyyyy, tr_x_xxyyy_yyyyyy, tr_x_xxyyy_yyyyyz, tr_x_xxyyy_yyyyz, tr_x_xxyyy_yyyyzz, tr_x_xxyyy_yyyzz, tr_x_xxyyy_yyyzzz, tr_x_xxyyy_yyzzz, tr_x_xxyyy_yyzzzz, tr_x_xxyyy_yzzzz, tr_x_xxyyy_yzzzzz, tr_x_xxyyyz_xxxxxx, tr_x_xxyyyz_xxxxxy, tr_x_xxyyyz_xxxxxz, tr_x_xxyyyz_xxxxyy, tr_x_xxyyyz_xxxxyz, tr_x_xxyyyz_xxxxzz, tr_x_xxyyyz_xxxyyy, tr_x_xxyyyz_xxxyyz, tr_x_xxyyyz_xxxyzz, tr_x_xxyyyz_xxxzzz, tr_x_xxyyyz_xxyyyy, tr_x_xxyyyz_xxyyyz, tr_x_xxyyyz_xxyyzz, tr_x_xxyyyz_xxyzzz, tr_x_xxyyyz_xxzzzz, tr_x_xxyyyz_xyyyyy, tr_x_xxyyyz_xyyyyz, tr_x_xxyyyz_xyyyzz, tr_x_xxyyyz_xyyzzz, tr_x_xxyyyz_xyzzzz, tr_x_xxyyyz_xzzzzz, tr_x_xxyyyz_yyyyyy, tr_x_xxyyyz_yyyyyz, tr_x_xxyyyz_yyyyzz, tr_x_xxyyyz_yyyzzz, tr_x_xxyyyz_yyzzzz, tr_x_xxyyyz_yzzzzz, tr_x_xxyyyz_zzzzzz, tr_x_xxyyz_xxxxxz, tr_x_xxyyz_xxxxzz, tr_x_xxyyz_xxxzzz, tr_x_xxyyz_xxzzzz, tr_x_xxyyz_xzzzzz, tr_x_xxyyz_zzzzzz, tr_x_xxyz_xxxxxz, tr_x_xxyz_xxxxzz, tr_x_xxyz_xxxzzz, tr_x_xxyz_xxzzzz, tr_x_xxyz_xzzzzz, tr_x_xxyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyz_xxxxxx[i] = tr_x_xxyyy_xxxxxx[i] * pa_z[i];

        tr_x_xxyyyz_xxxxxy[i] = tr_x_xxyyy_xxxxxy[i] * pa_z[i];

        tr_x_xxyyyz_xxxxxz[i] = 2.0 * tr_x_xxyz_xxxxxz[i] * fe_0 + tr_x_xxyyz_xxxxxz[i] * pa_y[i];

        tr_x_xxyyyz_xxxxyy[i] = tr_x_xxyyy_xxxxyy[i] * pa_z[i];

        tr_x_xxyyyz_xxxxyz[i] = tr_x_xxyyy_xxxxy[i] * fe_0 + tr_x_xxyyy_xxxxyz[i] * pa_z[i];

        tr_x_xxyyyz_xxxxzz[i] = 2.0 * tr_x_xxyz_xxxxzz[i] * fe_0 + tr_x_xxyyz_xxxxzz[i] * pa_y[i];

        tr_x_xxyyyz_xxxyyy[i] = tr_x_xxyyy_xxxyyy[i] * pa_z[i];

        tr_x_xxyyyz_xxxyyz[i] = tr_x_xxyyy_xxxyy[i] * fe_0 + tr_x_xxyyy_xxxyyz[i] * pa_z[i];

        tr_x_xxyyyz_xxxyzz[i] = 2.0 * tr_x_xxyyy_xxxyz[i] * fe_0 + tr_x_xxyyy_xxxyzz[i] * pa_z[i];

        tr_x_xxyyyz_xxxzzz[i] = 2.0 * tr_x_xxyz_xxxzzz[i] * fe_0 + tr_x_xxyyz_xxxzzz[i] * pa_y[i];

        tr_x_xxyyyz_xxyyyy[i] = tr_x_xxyyy_xxyyyy[i] * pa_z[i];

        tr_x_xxyyyz_xxyyyz[i] = tr_x_xxyyy_xxyyy[i] * fe_0 + tr_x_xxyyy_xxyyyz[i] * pa_z[i];

        tr_x_xxyyyz_xxyyzz[i] = 2.0 * tr_x_xxyyy_xxyyz[i] * fe_0 + tr_x_xxyyy_xxyyzz[i] * pa_z[i];

        tr_x_xxyyyz_xxyzzz[i] = 3.0 * tr_x_xxyyy_xxyzz[i] * fe_0 + tr_x_xxyyy_xxyzzz[i] * pa_z[i];

        tr_x_xxyyyz_xxzzzz[i] = 2.0 * tr_x_xxyz_xxzzzz[i] * fe_0 + tr_x_xxyyz_xxzzzz[i] * pa_y[i];

        tr_x_xxyyyz_xyyyyy[i] = tr_x_xxyyy_xyyyyy[i] * pa_z[i];

        tr_x_xxyyyz_xyyyyz[i] = tr_x_xxyyy_xyyyy[i] * fe_0 + tr_x_xxyyy_xyyyyz[i] * pa_z[i];

        tr_x_xxyyyz_xyyyzz[i] = 2.0 * tr_x_xxyyy_xyyyz[i] * fe_0 + tr_x_xxyyy_xyyyzz[i] * pa_z[i];

        tr_x_xxyyyz_xyyzzz[i] = 3.0 * tr_x_xxyyy_xyyzz[i] * fe_0 + tr_x_xxyyy_xyyzzz[i] * pa_z[i];

        tr_x_xxyyyz_xyzzzz[i] = 4.0 * tr_x_xxyyy_xyzzz[i] * fe_0 + tr_x_xxyyy_xyzzzz[i] * pa_z[i];

        tr_x_xxyyyz_xzzzzz[i] = 2.0 * tr_x_xxyz_xzzzzz[i] * fe_0 + tr_x_xxyyz_xzzzzz[i] * pa_y[i];

        tr_x_xxyyyz_yyyyyy[i] = tr_x_xxyyy_yyyyyy[i] * pa_z[i];

        tr_x_xxyyyz_yyyyyz[i] = tr_x_xxyyy_yyyyy[i] * fe_0 + tr_x_xxyyy_yyyyyz[i] * pa_z[i];

        tr_x_xxyyyz_yyyyzz[i] = 2.0 * tr_x_xxyyy_yyyyz[i] * fe_0 + tr_x_xxyyy_yyyyzz[i] * pa_z[i];

        tr_x_xxyyyz_yyyzzz[i] = 3.0 * tr_x_xxyyy_yyyzz[i] * fe_0 + tr_x_xxyyy_yyyzzz[i] * pa_z[i];

        tr_x_xxyyyz_yyzzzz[i] = 4.0 * tr_x_xxyyy_yyzzz[i] * fe_0 + tr_x_xxyyy_yyzzzz[i] * pa_z[i];

        tr_x_xxyyyz_yzzzzz[i] = 5.0 * tr_x_xxyyy_yzzzz[i] * fe_0 + tr_x_xxyyy_yzzzzz[i] * pa_z[i];

        tr_x_xxyyyz_zzzzzz[i] = 2.0 * tr_x_xxyz_zzzzzz[i] * fe_0 + tr_x_xxyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 336-364 components of targeted buffer : II

    auto tr_x_xxyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 336);

    auto tr_x_xxyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 337);

    auto tr_x_xxyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 338);

    auto tr_x_xxyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 339);

    auto tr_x_xxyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 340);

    auto tr_x_xxyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 341);

    auto tr_x_xxyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 342);

    auto tr_x_xxyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 343);

    auto tr_x_xxyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 344);

    auto tr_x_xxyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 345);

    auto tr_x_xxyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 346);

    auto tr_x_xxyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 347);

    auto tr_x_xxyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 348);

    auto tr_x_xxyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 349);

    auto tr_x_xxyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 350);

    auto tr_x_xxyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 351);

    auto tr_x_xxyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 352);

    auto tr_x_xxyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 353);

    auto tr_x_xxyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 354);

    auto tr_x_xxyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 355);

    auto tr_x_xxyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 356);

    auto tr_x_xxyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 357);

    auto tr_x_xxyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 358);

    auto tr_x_xxyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 359);

    auto tr_x_xxyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 360);

    auto tr_x_xxyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 361);

    auto tr_x_xxyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 362);

    auto tr_x_xxyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 363);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xxyy_xxxxxy, tr_x_xxyy_xxxxyy, tr_x_xxyy_xxxyyy, tr_x_xxyy_xxyyyy, tr_x_xxyy_xyyyyy, tr_x_xxyy_yyyyyy, tr_x_xxyyz_xxxxxy, tr_x_xxyyz_xxxxyy, tr_x_xxyyz_xxxyyy, tr_x_xxyyz_xxyyyy, tr_x_xxyyz_xyyyyy, tr_x_xxyyz_yyyyyy, tr_x_xxyyzz_xxxxxx, tr_x_xxyyzz_xxxxxy, tr_x_xxyyzz_xxxxxz, tr_x_xxyyzz_xxxxyy, tr_x_xxyyzz_xxxxyz, tr_x_xxyyzz_xxxxzz, tr_x_xxyyzz_xxxyyy, tr_x_xxyyzz_xxxyyz, tr_x_xxyyzz_xxxyzz, tr_x_xxyyzz_xxxzzz, tr_x_xxyyzz_xxyyyy, tr_x_xxyyzz_xxyyyz, tr_x_xxyyzz_xxyyzz, tr_x_xxyyzz_xxyzzz, tr_x_xxyyzz_xxzzzz, tr_x_xxyyzz_xyyyyy, tr_x_xxyyzz_xyyyyz, tr_x_xxyyzz_xyyyzz, tr_x_xxyyzz_xyyzzz, tr_x_xxyyzz_xyzzzz, tr_x_xxyyzz_xzzzzz, tr_x_xxyyzz_yyyyyy, tr_x_xxyyzz_yyyyyz, tr_x_xxyyzz_yyyyzz, tr_x_xxyyzz_yyyzzz, tr_x_xxyyzz_yyzzzz, tr_x_xxyyzz_yzzzzz, tr_x_xxyyzz_zzzzzz, tr_x_xxyzz_xxxxxx, tr_x_xxyzz_xxxxxz, tr_x_xxyzz_xxxxyz, tr_x_xxyzz_xxxxz, tr_x_xxyzz_xxxxzz, tr_x_xxyzz_xxxyyz, tr_x_xxyzz_xxxyz, tr_x_xxyzz_xxxyzz, tr_x_xxyzz_xxxzz, tr_x_xxyzz_xxxzzz, tr_x_xxyzz_xxyyyz, tr_x_xxyzz_xxyyz, tr_x_xxyzz_xxyyzz, tr_x_xxyzz_xxyzz, tr_x_xxyzz_xxyzzz, tr_x_xxyzz_xxzzz, tr_x_xxyzz_xxzzzz, tr_x_xxyzz_xyyyyz, tr_x_xxyzz_xyyyz, tr_x_xxyzz_xyyyzz, tr_x_xxyzz_xyyzz, tr_x_xxyzz_xyyzzz, tr_x_xxyzz_xyzzz, tr_x_xxyzz_xyzzzz, tr_x_xxyzz_xzzzz, tr_x_xxyzz_xzzzzz, tr_x_xxyzz_zzzzzz, tr_x_xxzz_xxxxxx, tr_x_xxzz_xxxxxz, tr_x_xxzz_xxxxyz, tr_x_xxzz_xxxxzz, tr_x_xxzz_xxxyyz, tr_x_xxzz_xxxyzz, tr_x_xxzz_xxxzzz, tr_x_xxzz_xxyyyz, tr_x_xxzz_xxyyzz, tr_x_xxzz_xxyzzz, tr_x_xxzz_xxzzzz, tr_x_xxzz_xyyyyz, tr_x_xxzz_xyyyzz, tr_x_xxzz_xyyzzz, tr_x_xxzz_xyzzzz, tr_x_xxzz_xzzzzz, tr_x_xxzz_zzzzzz, tr_x_xyyzz_yyyyyz, tr_x_xyyzz_yyyyzz, tr_x_xyyzz_yyyzzz, tr_x_xyyzz_yyzzzz, tr_x_xyyzz_yzzzzz, tr_x_yyzz_yyyyyz, tr_x_yyzz_yyyyzz, tr_x_yyzz_yyyzzz, tr_x_yyzz_yyzzzz, tr_x_yyzz_yzzzzz, ts_xyyzz_yyyyyz, ts_xyyzz_yyyyzz, ts_xyyzz_yyyzzz, ts_xyyzz_yyzzzz, ts_xyyzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyzz_xxxxxx[i] = tr_x_xxzz_xxxxxx[i] * fe_0 + tr_x_xxyzz_xxxxxx[i] * pa_y[i];

        tr_x_xxyyzz_xxxxxy[i] = tr_x_xxyy_xxxxxy[i] * fe_0 + tr_x_xxyyz_xxxxxy[i] * pa_z[i];

        tr_x_xxyyzz_xxxxxz[i] = tr_x_xxzz_xxxxxz[i] * fe_0 + tr_x_xxyzz_xxxxxz[i] * pa_y[i];

        tr_x_xxyyzz_xxxxyy[i] = tr_x_xxyy_xxxxyy[i] * fe_0 + tr_x_xxyyz_xxxxyy[i] * pa_z[i];

        tr_x_xxyyzz_xxxxyz[i] = tr_x_xxzz_xxxxyz[i] * fe_0 + tr_x_xxyzz_xxxxz[i] * fe_0 + tr_x_xxyzz_xxxxyz[i] * pa_y[i];

        tr_x_xxyyzz_xxxxzz[i] = tr_x_xxzz_xxxxzz[i] * fe_0 + tr_x_xxyzz_xxxxzz[i] * pa_y[i];

        tr_x_xxyyzz_xxxyyy[i] = tr_x_xxyy_xxxyyy[i] * fe_0 + tr_x_xxyyz_xxxyyy[i] * pa_z[i];

        tr_x_xxyyzz_xxxyyz[i] = tr_x_xxzz_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxyzz_xxxyz[i] * fe_0 + tr_x_xxyzz_xxxyyz[i] * pa_y[i];

        tr_x_xxyyzz_xxxyzz[i] = tr_x_xxzz_xxxyzz[i] * fe_0 + tr_x_xxyzz_xxxzz[i] * fe_0 + tr_x_xxyzz_xxxyzz[i] * pa_y[i];

        tr_x_xxyyzz_xxxzzz[i] = tr_x_xxzz_xxxzzz[i] * fe_0 + tr_x_xxyzz_xxxzzz[i] * pa_y[i];

        tr_x_xxyyzz_xxyyyy[i] = tr_x_xxyy_xxyyyy[i] * fe_0 + tr_x_xxyyz_xxyyyy[i] * pa_z[i];

        tr_x_xxyyzz_xxyyyz[i] = tr_x_xxzz_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxyzz_xxyyz[i] * fe_0 + tr_x_xxyzz_xxyyyz[i] * pa_y[i];

        tr_x_xxyyzz_xxyyzz[i] = tr_x_xxzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxyzz_xxyzz[i] * fe_0 + tr_x_xxyzz_xxyyzz[i] * pa_y[i];

        tr_x_xxyyzz_xxyzzz[i] = tr_x_xxzz_xxyzzz[i] * fe_0 + tr_x_xxyzz_xxzzz[i] * fe_0 + tr_x_xxyzz_xxyzzz[i] * pa_y[i];

        tr_x_xxyyzz_xxzzzz[i] = tr_x_xxzz_xxzzzz[i] * fe_0 + tr_x_xxyzz_xxzzzz[i] * pa_y[i];

        tr_x_xxyyzz_xyyyyy[i] = tr_x_xxyy_xyyyyy[i] * fe_0 + tr_x_xxyyz_xyyyyy[i] * pa_z[i];

        tr_x_xxyyzz_xyyyyz[i] = tr_x_xxzz_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxyzz_xyyyz[i] * fe_0 + tr_x_xxyzz_xyyyyz[i] * pa_y[i];

        tr_x_xxyyzz_xyyyzz[i] = tr_x_xxzz_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxyzz_xyyzz[i] * fe_0 + tr_x_xxyzz_xyyyzz[i] * pa_y[i];

        tr_x_xxyyzz_xyyzzz[i] = tr_x_xxzz_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxyzz_xyzzz[i] * fe_0 + tr_x_xxyzz_xyyzzz[i] * pa_y[i];

        tr_x_xxyyzz_xyzzzz[i] = tr_x_xxzz_xyzzzz[i] * fe_0 + tr_x_xxyzz_xzzzz[i] * fe_0 + tr_x_xxyzz_xyzzzz[i] * pa_y[i];

        tr_x_xxyyzz_xzzzzz[i] = tr_x_xxzz_xzzzzz[i] * fe_0 + tr_x_xxyzz_xzzzzz[i] * pa_y[i];

        tr_x_xxyyzz_yyyyyy[i] = tr_x_xxyy_yyyyyy[i] * fe_0 + tr_x_xxyyz_yyyyyy[i] * pa_z[i];

        tr_x_xxyyzz_yyyyyz[i] = tr_x_yyzz_yyyyyz[i] * fe_0 + ts_xyyzz_yyyyyz[i] * fe_0 + tr_x_xyyzz_yyyyyz[i] * pa_x[i];

        tr_x_xxyyzz_yyyyzz[i] = tr_x_yyzz_yyyyzz[i] * fe_0 + ts_xyyzz_yyyyzz[i] * fe_0 + tr_x_xyyzz_yyyyzz[i] * pa_x[i];

        tr_x_xxyyzz_yyyzzz[i] = tr_x_yyzz_yyyzzz[i] * fe_0 + ts_xyyzz_yyyzzz[i] * fe_0 + tr_x_xyyzz_yyyzzz[i] * pa_x[i];

        tr_x_xxyyzz_yyzzzz[i] = tr_x_yyzz_yyzzzz[i] * fe_0 + ts_xyyzz_yyzzzz[i] * fe_0 + tr_x_xyyzz_yyzzzz[i] * pa_x[i];

        tr_x_xxyyzz_yzzzzz[i] = tr_x_yyzz_yzzzzz[i] * fe_0 + ts_xyyzz_yzzzzz[i] * fe_0 + tr_x_xyyzz_yzzzzz[i] * pa_x[i];

        tr_x_xxyyzz_zzzzzz[i] = tr_x_xxzz_zzzzzz[i] * fe_0 + tr_x_xxyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 364-392 components of targeted buffer : II

    auto tr_x_xxyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 364);

    auto tr_x_xxyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 365);

    auto tr_x_xxyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 366);

    auto tr_x_xxyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 367);

    auto tr_x_xxyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 368);

    auto tr_x_xxyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 369);

    auto tr_x_xxyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 370);

    auto tr_x_xxyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 371);

    auto tr_x_xxyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 372);

    auto tr_x_xxyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 373);

    auto tr_x_xxyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 374);

    auto tr_x_xxyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 375);

    auto tr_x_xxyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 376);

    auto tr_x_xxyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 377);

    auto tr_x_xxyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 378);

    auto tr_x_xxyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 379);

    auto tr_x_xxyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 380);

    auto tr_x_xxyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 381);

    auto tr_x_xxyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 382);

    auto tr_x_xxyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 383);

    auto tr_x_xxyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 384);

    auto tr_x_xxyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 385);

    auto tr_x_xxyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 386);

    auto tr_x_xxyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 387);

    auto tr_x_xxyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 388);

    auto tr_x_xxyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 389);

    auto tr_x_xxyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 390);

    auto tr_x_xxyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 391);

    #pragma omp simd aligned(pa_y, tr_x_xxyzzz_xxxxxx, tr_x_xxyzzz_xxxxxy, tr_x_xxyzzz_xxxxxz, tr_x_xxyzzz_xxxxyy, tr_x_xxyzzz_xxxxyz, tr_x_xxyzzz_xxxxzz, tr_x_xxyzzz_xxxyyy, tr_x_xxyzzz_xxxyyz, tr_x_xxyzzz_xxxyzz, tr_x_xxyzzz_xxxzzz, tr_x_xxyzzz_xxyyyy, tr_x_xxyzzz_xxyyyz, tr_x_xxyzzz_xxyyzz, tr_x_xxyzzz_xxyzzz, tr_x_xxyzzz_xxzzzz, tr_x_xxyzzz_xyyyyy, tr_x_xxyzzz_xyyyyz, tr_x_xxyzzz_xyyyzz, tr_x_xxyzzz_xyyzzz, tr_x_xxyzzz_xyzzzz, tr_x_xxyzzz_xzzzzz, tr_x_xxyzzz_yyyyyy, tr_x_xxyzzz_yyyyyz, tr_x_xxyzzz_yyyyzz, tr_x_xxyzzz_yyyzzz, tr_x_xxyzzz_yyzzzz, tr_x_xxyzzz_yzzzzz, tr_x_xxyzzz_zzzzzz, tr_x_xxzzz_xxxxx, tr_x_xxzzz_xxxxxx, tr_x_xxzzz_xxxxxy, tr_x_xxzzz_xxxxxz, tr_x_xxzzz_xxxxy, tr_x_xxzzz_xxxxyy, tr_x_xxzzz_xxxxyz, tr_x_xxzzz_xxxxz, tr_x_xxzzz_xxxxzz, tr_x_xxzzz_xxxyy, tr_x_xxzzz_xxxyyy, tr_x_xxzzz_xxxyyz, tr_x_xxzzz_xxxyz, tr_x_xxzzz_xxxyzz, tr_x_xxzzz_xxxzz, tr_x_xxzzz_xxxzzz, tr_x_xxzzz_xxyyy, tr_x_xxzzz_xxyyyy, tr_x_xxzzz_xxyyyz, tr_x_xxzzz_xxyyz, tr_x_xxzzz_xxyyzz, tr_x_xxzzz_xxyzz, tr_x_xxzzz_xxyzzz, tr_x_xxzzz_xxzzz, tr_x_xxzzz_xxzzzz, tr_x_xxzzz_xyyyy, tr_x_xxzzz_xyyyyy, tr_x_xxzzz_xyyyyz, tr_x_xxzzz_xyyyz, tr_x_xxzzz_xyyyzz, tr_x_xxzzz_xyyzz, tr_x_xxzzz_xyyzzz, tr_x_xxzzz_xyzzz, tr_x_xxzzz_xyzzzz, tr_x_xxzzz_xzzzz, tr_x_xxzzz_xzzzzz, tr_x_xxzzz_yyyyy, tr_x_xxzzz_yyyyyy, tr_x_xxzzz_yyyyyz, tr_x_xxzzz_yyyyz, tr_x_xxzzz_yyyyzz, tr_x_xxzzz_yyyzz, tr_x_xxzzz_yyyzzz, tr_x_xxzzz_yyzzz, tr_x_xxzzz_yyzzzz, tr_x_xxzzz_yzzzz, tr_x_xxzzz_yzzzzz, tr_x_xxzzz_zzzzz, tr_x_xxzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzzz_xxxxxx[i] = tr_x_xxzzz_xxxxxx[i] * pa_y[i];

        tr_x_xxyzzz_xxxxxy[i] = tr_x_xxzzz_xxxxx[i] * fe_0 + tr_x_xxzzz_xxxxxy[i] * pa_y[i];

        tr_x_xxyzzz_xxxxxz[i] = tr_x_xxzzz_xxxxxz[i] * pa_y[i];

        tr_x_xxyzzz_xxxxyy[i] = 2.0 * tr_x_xxzzz_xxxxy[i] * fe_0 + tr_x_xxzzz_xxxxyy[i] * pa_y[i];

        tr_x_xxyzzz_xxxxyz[i] = tr_x_xxzzz_xxxxz[i] * fe_0 + tr_x_xxzzz_xxxxyz[i] * pa_y[i];

        tr_x_xxyzzz_xxxxzz[i] = tr_x_xxzzz_xxxxzz[i] * pa_y[i];

        tr_x_xxyzzz_xxxyyy[i] = 3.0 * tr_x_xxzzz_xxxyy[i] * fe_0 + tr_x_xxzzz_xxxyyy[i] * pa_y[i];

        tr_x_xxyzzz_xxxyyz[i] = 2.0 * tr_x_xxzzz_xxxyz[i] * fe_0 + tr_x_xxzzz_xxxyyz[i] * pa_y[i];

        tr_x_xxyzzz_xxxyzz[i] = tr_x_xxzzz_xxxzz[i] * fe_0 + tr_x_xxzzz_xxxyzz[i] * pa_y[i];

        tr_x_xxyzzz_xxxzzz[i] = tr_x_xxzzz_xxxzzz[i] * pa_y[i];

        tr_x_xxyzzz_xxyyyy[i] = 4.0 * tr_x_xxzzz_xxyyy[i] * fe_0 + tr_x_xxzzz_xxyyyy[i] * pa_y[i];

        tr_x_xxyzzz_xxyyyz[i] = 3.0 * tr_x_xxzzz_xxyyz[i] * fe_0 + tr_x_xxzzz_xxyyyz[i] * pa_y[i];

        tr_x_xxyzzz_xxyyzz[i] = 2.0 * tr_x_xxzzz_xxyzz[i] * fe_0 + tr_x_xxzzz_xxyyzz[i] * pa_y[i];

        tr_x_xxyzzz_xxyzzz[i] = tr_x_xxzzz_xxzzz[i] * fe_0 + tr_x_xxzzz_xxyzzz[i] * pa_y[i];

        tr_x_xxyzzz_xxzzzz[i] = tr_x_xxzzz_xxzzzz[i] * pa_y[i];

        tr_x_xxyzzz_xyyyyy[i] = 5.0 * tr_x_xxzzz_xyyyy[i] * fe_0 + tr_x_xxzzz_xyyyyy[i] * pa_y[i];

        tr_x_xxyzzz_xyyyyz[i] = 4.0 * tr_x_xxzzz_xyyyz[i] * fe_0 + tr_x_xxzzz_xyyyyz[i] * pa_y[i];

        tr_x_xxyzzz_xyyyzz[i] = 3.0 * tr_x_xxzzz_xyyzz[i] * fe_0 + tr_x_xxzzz_xyyyzz[i] * pa_y[i];

        tr_x_xxyzzz_xyyzzz[i] = 2.0 * tr_x_xxzzz_xyzzz[i] * fe_0 + tr_x_xxzzz_xyyzzz[i] * pa_y[i];

        tr_x_xxyzzz_xyzzzz[i] = tr_x_xxzzz_xzzzz[i] * fe_0 + tr_x_xxzzz_xyzzzz[i] * pa_y[i];

        tr_x_xxyzzz_xzzzzz[i] = tr_x_xxzzz_xzzzzz[i] * pa_y[i];

        tr_x_xxyzzz_yyyyyy[i] = 6.0 * tr_x_xxzzz_yyyyy[i] * fe_0 + tr_x_xxzzz_yyyyyy[i] * pa_y[i];

        tr_x_xxyzzz_yyyyyz[i] = 5.0 * tr_x_xxzzz_yyyyz[i] * fe_0 + tr_x_xxzzz_yyyyyz[i] * pa_y[i];

        tr_x_xxyzzz_yyyyzz[i] = 4.0 * tr_x_xxzzz_yyyzz[i] * fe_0 + tr_x_xxzzz_yyyyzz[i] * pa_y[i];

        tr_x_xxyzzz_yyyzzz[i] = 3.0 * tr_x_xxzzz_yyzzz[i] * fe_0 + tr_x_xxzzz_yyyzzz[i] * pa_y[i];

        tr_x_xxyzzz_yyzzzz[i] = 2.0 * tr_x_xxzzz_yzzzz[i] * fe_0 + tr_x_xxzzz_yyzzzz[i] * pa_y[i];

        tr_x_xxyzzz_yzzzzz[i] = tr_x_xxzzz_zzzzz[i] * fe_0 + tr_x_xxzzz_yzzzzz[i] * pa_y[i];

        tr_x_xxyzzz_zzzzzz[i] = tr_x_xxzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : II

    auto tr_x_xxzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 392);

    auto tr_x_xxzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 393);

    auto tr_x_xxzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 394);

    auto tr_x_xxzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 395);

    auto tr_x_xxzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 396);

    auto tr_x_xxzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 397);

    auto tr_x_xxzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 398);

    auto tr_x_xxzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 399);

    auto tr_x_xxzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 400);

    auto tr_x_xxzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 401);

    auto tr_x_xxzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 402);

    auto tr_x_xxzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 403);

    auto tr_x_xxzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 404);

    auto tr_x_xxzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 405);

    auto tr_x_xxzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 406);

    auto tr_x_xxzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 407);

    auto tr_x_xxzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 408);

    auto tr_x_xxzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 409);

    auto tr_x_xxzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 410);

    auto tr_x_xxzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 411);

    auto tr_x_xxzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 412);

    auto tr_x_xxzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 413);

    auto tr_x_xxzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 414);

    auto tr_x_xxzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 415);

    auto tr_x_xxzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 416);

    auto tr_x_xxzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 417);

    auto tr_x_xxzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 418);

    auto tr_x_xxzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 419);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxzz_xxxxxx, tr_x_xxzz_xxxxxy, tr_x_xxzz_xxxxxz, tr_x_xxzz_xxxxyy, tr_x_xxzz_xxxxyz, tr_x_xxzz_xxxxzz, tr_x_xxzz_xxxyyy, tr_x_xxzz_xxxyyz, tr_x_xxzz_xxxyzz, tr_x_xxzz_xxxzzz, tr_x_xxzz_xxyyyy, tr_x_xxzz_xxyyyz, tr_x_xxzz_xxyyzz, tr_x_xxzz_xxyzzz, tr_x_xxzz_xxzzzz, tr_x_xxzz_xyyyyy, tr_x_xxzz_xyyyyz, tr_x_xxzz_xyyyzz, tr_x_xxzz_xyyzzz, tr_x_xxzz_xyzzzz, tr_x_xxzz_xzzzzz, tr_x_xxzz_yyyyyy, tr_x_xxzzz_xxxxx, tr_x_xxzzz_xxxxxx, tr_x_xxzzz_xxxxxy, tr_x_xxzzz_xxxxxz, tr_x_xxzzz_xxxxy, tr_x_xxzzz_xxxxyy, tr_x_xxzzz_xxxxyz, tr_x_xxzzz_xxxxz, tr_x_xxzzz_xxxxzz, tr_x_xxzzz_xxxyy, tr_x_xxzzz_xxxyyy, tr_x_xxzzz_xxxyyz, tr_x_xxzzz_xxxyz, tr_x_xxzzz_xxxyzz, tr_x_xxzzz_xxxzz, tr_x_xxzzz_xxxzzz, tr_x_xxzzz_xxyyy, tr_x_xxzzz_xxyyyy, tr_x_xxzzz_xxyyyz, tr_x_xxzzz_xxyyz, tr_x_xxzzz_xxyyzz, tr_x_xxzzz_xxyzz, tr_x_xxzzz_xxyzzz, tr_x_xxzzz_xxzzz, tr_x_xxzzz_xxzzzz, tr_x_xxzzz_xyyyy, tr_x_xxzzz_xyyyyy, tr_x_xxzzz_xyyyyz, tr_x_xxzzz_xyyyz, tr_x_xxzzz_xyyyzz, tr_x_xxzzz_xyyzz, tr_x_xxzzz_xyyzzz, tr_x_xxzzz_xyzzz, tr_x_xxzzz_xyzzzz, tr_x_xxzzz_xzzzz, tr_x_xxzzz_xzzzzz, tr_x_xxzzz_yyyyyy, tr_x_xxzzzz_xxxxxx, tr_x_xxzzzz_xxxxxy, tr_x_xxzzzz_xxxxxz, tr_x_xxzzzz_xxxxyy, tr_x_xxzzzz_xxxxyz, tr_x_xxzzzz_xxxxzz, tr_x_xxzzzz_xxxyyy, tr_x_xxzzzz_xxxyyz, tr_x_xxzzzz_xxxyzz, tr_x_xxzzzz_xxxzzz, tr_x_xxzzzz_xxyyyy, tr_x_xxzzzz_xxyyyz, tr_x_xxzzzz_xxyyzz, tr_x_xxzzzz_xxyzzz, tr_x_xxzzzz_xxzzzz, tr_x_xxzzzz_xyyyyy, tr_x_xxzzzz_xyyyyz, tr_x_xxzzzz_xyyyzz, tr_x_xxzzzz_xyyzzz, tr_x_xxzzzz_xyzzzz, tr_x_xxzzzz_xzzzzz, tr_x_xxzzzz_yyyyyy, tr_x_xxzzzz_yyyyyz, tr_x_xxzzzz_yyyyzz, tr_x_xxzzzz_yyyzzz, tr_x_xxzzzz_yyzzzz, tr_x_xxzzzz_yzzzzz, tr_x_xxzzzz_zzzzzz, tr_x_xzzzz_yyyyyz, tr_x_xzzzz_yyyyzz, tr_x_xzzzz_yyyzzz, tr_x_xzzzz_yyzzzz, tr_x_xzzzz_yzzzzz, tr_x_xzzzz_zzzzzz, tr_x_zzzz_yyyyyz, tr_x_zzzz_yyyyzz, tr_x_zzzz_yyyzzz, tr_x_zzzz_yyzzzz, tr_x_zzzz_yzzzzz, tr_x_zzzz_zzzzzz, ts_xzzzz_yyyyyz, ts_xzzzz_yyyyzz, ts_xzzzz_yyyzzz, ts_xzzzz_yyzzzz, ts_xzzzz_yzzzzz, ts_xzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzzz_xxxxxx[i] = 3.0 * tr_x_xxzz_xxxxxx[i] * fe_0 + tr_x_xxzzz_xxxxxx[i] * pa_z[i];

        tr_x_xxzzzz_xxxxxy[i] = 3.0 * tr_x_xxzz_xxxxxy[i] * fe_0 + tr_x_xxzzz_xxxxxy[i] * pa_z[i];

        tr_x_xxzzzz_xxxxxz[i] = 3.0 * tr_x_xxzz_xxxxxz[i] * fe_0 + tr_x_xxzzz_xxxxx[i] * fe_0 + tr_x_xxzzz_xxxxxz[i] * pa_z[i];

        tr_x_xxzzzz_xxxxyy[i] = 3.0 * tr_x_xxzz_xxxxyy[i] * fe_0 + tr_x_xxzzz_xxxxyy[i] * pa_z[i];

        tr_x_xxzzzz_xxxxyz[i] = 3.0 * tr_x_xxzz_xxxxyz[i] * fe_0 + tr_x_xxzzz_xxxxy[i] * fe_0 + tr_x_xxzzz_xxxxyz[i] * pa_z[i];

        tr_x_xxzzzz_xxxxzz[i] = 3.0 * tr_x_xxzz_xxxxzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xxxxz[i] * fe_0 + tr_x_xxzzz_xxxxzz[i] * pa_z[i];

        tr_x_xxzzzz_xxxyyy[i] = 3.0 * tr_x_xxzz_xxxyyy[i] * fe_0 + tr_x_xxzzz_xxxyyy[i] * pa_z[i];

        tr_x_xxzzzz_xxxyyz[i] = 3.0 * tr_x_xxzz_xxxyyz[i] * fe_0 + tr_x_xxzzz_xxxyy[i] * fe_0 + tr_x_xxzzz_xxxyyz[i] * pa_z[i];

        tr_x_xxzzzz_xxxyzz[i] = 3.0 * tr_x_xxzz_xxxyzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xxxyz[i] * fe_0 + tr_x_xxzzz_xxxyzz[i] * pa_z[i];

        tr_x_xxzzzz_xxxzzz[i] = 3.0 * tr_x_xxzz_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxzzz_xxxzz[i] * fe_0 + tr_x_xxzzz_xxxzzz[i] * pa_z[i];

        tr_x_xxzzzz_xxyyyy[i] = 3.0 * tr_x_xxzz_xxyyyy[i] * fe_0 + tr_x_xxzzz_xxyyyy[i] * pa_z[i];

        tr_x_xxzzzz_xxyyyz[i] = 3.0 * tr_x_xxzz_xxyyyz[i] * fe_0 + tr_x_xxzzz_xxyyy[i] * fe_0 + tr_x_xxzzz_xxyyyz[i] * pa_z[i];

        tr_x_xxzzzz_xxyyzz[i] = 3.0 * tr_x_xxzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xxyyz[i] * fe_0 + tr_x_xxzzz_xxyyzz[i] * pa_z[i];

        tr_x_xxzzzz_xxyzzz[i] = 3.0 * tr_x_xxzz_xxyzzz[i] * fe_0 + 3.0 * tr_x_xxzzz_xxyzz[i] * fe_0 + tr_x_xxzzz_xxyzzz[i] * pa_z[i];

        tr_x_xxzzzz_xxzzzz[i] = 3.0 * tr_x_xxzz_xxzzzz[i] * fe_0 + 4.0 * tr_x_xxzzz_xxzzz[i] * fe_0 + tr_x_xxzzz_xxzzzz[i] * pa_z[i];

        tr_x_xxzzzz_xyyyyy[i] = 3.0 * tr_x_xxzz_xyyyyy[i] * fe_0 + tr_x_xxzzz_xyyyyy[i] * pa_z[i];

        tr_x_xxzzzz_xyyyyz[i] = 3.0 * tr_x_xxzz_xyyyyz[i] * fe_0 + tr_x_xxzzz_xyyyy[i] * fe_0 + tr_x_xxzzz_xyyyyz[i] * pa_z[i];

        tr_x_xxzzzz_xyyyzz[i] = 3.0 * tr_x_xxzz_xyyyzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xyyyz[i] * fe_0 + tr_x_xxzzz_xyyyzz[i] * pa_z[i];

        tr_x_xxzzzz_xyyzzz[i] = 3.0 * tr_x_xxzz_xyyzzz[i] * fe_0 + 3.0 * tr_x_xxzzz_xyyzz[i] * fe_0 + tr_x_xxzzz_xyyzzz[i] * pa_z[i];

        tr_x_xxzzzz_xyzzzz[i] = 3.0 * tr_x_xxzz_xyzzzz[i] * fe_0 + 4.0 * tr_x_xxzzz_xyzzz[i] * fe_0 + tr_x_xxzzz_xyzzzz[i] * pa_z[i];

        tr_x_xxzzzz_xzzzzz[i] = 3.0 * tr_x_xxzz_xzzzzz[i] * fe_0 + 5.0 * tr_x_xxzzz_xzzzz[i] * fe_0 + tr_x_xxzzz_xzzzzz[i] * pa_z[i];

        tr_x_xxzzzz_yyyyyy[i] = 3.0 * tr_x_xxzz_yyyyyy[i] * fe_0 + tr_x_xxzzz_yyyyyy[i] * pa_z[i];

        tr_x_xxzzzz_yyyyyz[i] = tr_x_zzzz_yyyyyz[i] * fe_0 + ts_xzzzz_yyyyyz[i] * fe_0 + tr_x_xzzzz_yyyyyz[i] * pa_x[i];

        tr_x_xxzzzz_yyyyzz[i] = tr_x_zzzz_yyyyzz[i] * fe_0 + ts_xzzzz_yyyyzz[i] * fe_0 + tr_x_xzzzz_yyyyzz[i] * pa_x[i];

        tr_x_xxzzzz_yyyzzz[i] = tr_x_zzzz_yyyzzz[i] * fe_0 + ts_xzzzz_yyyzzz[i] * fe_0 + tr_x_xzzzz_yyyzzz[i] * pa_x[i];

        tr_x_xxzzzz_yyzzzz[i] = tr_x_zzzz_yyzzzz[i] * fe_0 + ts_xzzzz_yyzzzz[i] * fe_0 + tr_x_xzzzz_yyzzzz[i] * pa_x[i];

        tr_x_xxzzzz_yzzzzz[i] = tr_x_zzzz_yzzzzz[i] * fe_0 + ts_xzzzz_yzzzzz[i] * fe_0 + tr_x_xzzzz_yzzzzz[i] * pa_x[i];

        tr_x_xxzzzz_zzzzzz[i] = tr_x_zzzz_zzzzzz[i] * fe_0 + ts_xzzzz_zzzzzz[i] * fe_0 + tr_x_xzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : II

    auto tr_x_xyyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 420);

    auto tr_x_xyyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 421);

    auto tr_x_xyyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 422);

    auto tr_x_xyyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 423);

    auto tr_x_xyyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 424);

    auto tr_x_xyyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 425);

    auto tr_x_xyyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 426);

    auto tr_x_xyyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 427);

    auto tr_x_xyyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 428);

    auto tr_x_xyyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 429);

    auto tr_x_xyyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 430);

    auto tr_x_xyyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 431);

    auto tr_x_xyyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 432);

    auto tr_x_xyyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 433);

    auto tr_x_xyyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 434);

    auto tr_x_xyyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 435);

    auto tr_x_xyyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 436);

    auto tr_x_xyyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 437);

    auto tr_x_xyyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 438);

    auto tr_x_xyyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 439);

    auto tr_x_xyyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 440);

    auto tr_x_xyyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 441);

    auto tr_x_xyyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 442);

    auto tr_x_xyyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 443);

    auto tr_x_xyyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 444);

    auto tr_x_xyyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 445);

    auto tr_x_xyyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 446);

    auto tr_x_xyyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 447);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyyy_xxxxxx, tr_x_xyyy_xxxxxz, tr_x_xyyy_xxxxzz, tr_x_xyyy_xxxzzz, tr_x_xyyy_xxzzzz, tr_x_xyyy_xzzzzz, tr_x_xyyyy_xxxxxx, tr_x_xyyyy_xxxxxz, tr_x_xyyyy_xxxxzz, tr_x_xyyyy_xxxzzz, tr_x_xyyyy_xxzzzz, tr_x_xyyyy_xzzzzz, tr_x_xyyyyy_xxxxxx, tr_x_xyyyyy_xxxxxy, tr_x_xyyyyy_xxxxxz, tr_x_xyyyyy_xxxxyy, tr_x_xyyyyy_xxxxyz, tr_x_xyyyyy_xxxxzz, tr_x_xyyyyy_xxxyyy, tr_x_xyyyyy_xxxyyz, tr_x_xyyyyy_xxxyzz, tr_x_xyyyyy_xxxzzz, tr_x_xyyyyy_xxyyyy, tr_x_xyyyyy_xxyyyz, tr_x_xyyyyy_xxyyzz, tr_x_xyyyyy_xxyzzz, tr_x_xyyyyy_xxzzzz, tr_x_xyyyyy_xyyyyy, tr_x_xyyyyy_xyyyyz, tr_x_xyyyyy_xyyyzz, tr_x_xyyyyy_xyyzzz, tr_x_xyyyyy_xyzzzz, tr_x_xyyyyy_xzzzzz, tr_x_xyyyyy_yyyyyy, tr_x_xyyyyy_yyyyyz, tr_x_xyyyyy_yyyyzz, tr_x_xyyyyy_yyyzzz, tr_x_xyyyyy_yyzzzz, tr_x_xyyyyy_yzzzzz, tr_x_xyyyyy_zzzzzz, tr_x_yyyyy_xxxxxy, tr_x_yyyyy_xxxxy, tr_x_yyyyy_xxxxyy, tr_x_yyyyy_xxxxyz, tr_x_yyyyy_xxxyy, tr_x_yyyyy_xxxyyy, tr_x_yyyyy_xxxyyz, tr_x_yyyyy_xxxyz, tr_x_yyyyy_xxxyzz, tr_x_yyyyy_xxyyy, tr_x_yyyyy_xxyyyy, tr_x_yyyyy_xxyyyz, tr_x_yyyyy_xxyyz, tr_x_yyyyy_xxyyzz, tr_x_yyyyy_xxyzz, tr_x_yyyyy_xxyzzz, tr_x_yyyyy_xyyyy, tr_x_yyyyy_xyyyyy, tr_x_yyyyy_xyyyyz, tr_x_yyyyy_xyyyz, tr_x_yyyyy_xyyyzz, tr_x_yyyyy_xyyzz, tr_x_yyyyy_xyyzzz, tr_x_yyyyy_xyzzz, tr_x_yyyyy_xyzzzz, tr_x_yyyyy_yyyyy, tr_x_yyyyy_yyyyyy, tr_x_yyyyy_yyyyyz, tr_x_yyyyy_yyyyz, tr_x_yyyyy_yyyyzz, tr_x_yyyyy_yyyzz, tr_x_yyyyy_yyyzzz, tr_x_yyyyy_yyzzz, tr_x_yyyyy_yyzzzz, tr_x_yyyyy_yzzzz, tr_x_yyyyy_yzzzzz, tr_x_yyyyy_zzzzzz, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxyy, ts_yyyyy_xxxxyz, ts_yyyyy_xxxyyy, ts_yyyyy_xxxyyz, ts_yyyyy_xxxyzz, ts_yyyyy_xxyyyy, ts_yyyyy_xxyyyz, ts_yyyyy_xxyyzz, ts_yyyyy_xxyzzz, ts_yyyyy_xyyyyy, ts_yyyyy_xyyyyz, ts_yyyyy_xyyyzz, ts_yyyyy_xyyzzz, ts_yyyyy_xyzzzz, ts_yyyyy_yyyyyy, ts_yyyyy_yyyyyz, ts_yyyyy_yyyyzz, ts_yyyyy_yyyzzz, ts_yyyyy_yyzzzz, ts_yyyyy_yzzzzz, ts_yyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyy_xxxxxx[i] = 4.0 * tr_x_xyyy_xxxxxx[i] * fe_0 + tr_x_xyyyy_xxxxxx[i] * pa_y[i];

        tr_x_xyyyyy_xxxxxy[i] = 5.0 * tr_x_yyyyy_xxxxy[i] * fe_0 + ts_yyyyy_xxxxxy[i] * fe_0 + tr_x_yyyyy_xxxxxy[i] * pa_x[i];

        tr_x_xyyyyy_xxxxxz[i] = 4.0 * tr_x_xyyy_xxxxxz[i] * fe_0 + tr_x_xyyyy_xxxxxz[i] * pa_y[i];

        tr_x_xyyyyy_xxxxyy[i] = 4.0 * tr_x_yyyyy_xxxyy[i] * fe_0 + ts_yyyyy_xxxxyy[i] * fe_0 + tr_x_yyyyy_xxxxyy[i] * pa_x[i];

        tr_x_xyyyyy_xxxxyz[i] = 4.0 * tr_x_yyyyy_xxxyz[i] * fe_0 + ts_yyyyy_xxxxyz[i] * fe_0 + tr_x_yyyyy_xxxxyz[i] * pa_x[i];

        tr_x_xyyyyy_xxxxzz[i] = 4.0 * tr_x_xyyy_xxxxzz[i] * fe_0 + tr_x_xyyyy_xxxxzz[i] * pa_y[i];

        tr_x_xyyyyy_xxxyyy[i] = 3.0 * tr_x_yyyyy_xxyyy[i] * fe_0 + ts_yyyyy_xxxyyy[i] * fe_0 + tr_x_yyyyy_xxxyyy[i] * pa_x[i];

        tr_x_xyyyyy_xxxyyz[i] = 3.0 * tr_x_yyyyy_xxyyz[i] * fe_0 + ts_yyyyy_xxxyyz[i] * fe_0 + tr_x_yyyyy_xxxyyz[i] * pa_x[i];

        tr_x_xyyyyy_xxxyzz[i] = 3.0 * tr_x_yyyyy_xxyzz[i] * fe_0 + ts_yyyyy_xxxyzz[i] * fe_0 + tr_x_yyyyy_xxxyzz[i] * pa_x[i];

        tr_x_xyyyyy_xxxzzz[i] = 4.0 * tr_x_xyyy_xxxzzz[i] * fe_0 + tr_x_xyyyy_xxxzzz[i] * pa_y[i];

        tr_x_xyyyyy_xxyyyy[i] = 2.0 * tr_x_yyyyy_xyyyy[i] * fe_0 + ts_yyyyy_xxyyyy[i] * fe_0 + tr_x_yyyyy_xxyyyy[i] * pa_x[i];

        tr_x_xyyyyy_xxyyyz[i] = 2.0 * tr_x_yyyyy_xyyyz[i] * fe_0 + ts_yyyyy_xxyyyz[i] * fe_0 + tr_x_yyyyy_xxyyyz[i] * pa_x[i];

        tr_x_xyyyyy_xxyyzz[i] = 2.0 * tr_x_yyyyy_xyyzz[i] * fe_0 + ts_yyyyy_xxyyzz[i] * fe_0 + tr_x_yyyyy_xxyyzz[i] * pa_x[i];

        tr_x_xyyyyy_xxyzzz[i] = 2.0 * tr_x_yyyyy_xyzzz[i] * fe_0 + ts_yyyyy_xxyzzz[i] * fe_0 + tr_x_yyyyy_xxyzzz[i] * pa_x[i];

        tr_x_xyyyyy_xxzzzz[i] = 4.0 * tr_x_xyyy_xxzzzz[i] * fe_0 + tr_x_xyyyy_xxzzzz[i] * pa_y[i];

        tr_x_xyyyyy_xyyyyy[i] = tr_x_yyyyy_yyyyy[i] * fe_0 + ts_yyyyy_xyyyyy[i] * fe_0 + tr_x_yyyyy_xyyyyy[i] * pa_x[i];

        tr_x_xyyyyy_xyyyyz[i] = tr_x_yyyyy_yyyyz[i] * fe_0 + ts_yyyyy_xyyyyz[i] * fe_0 + tr_x_yyyyy_xyyyyz[i] * pa_x[i];

        tr_x_xyyyyy_xyyyzz[i] = tr_x_yyyyy_yyyzz[i] * fe_0 + ts_yyyyy_xyyyzz[i] * fe_0 + tr_x_yyyyy_xyyyzz[i] * pa_x[i];

        tr_x_xyyyyy_xyyzzz[i] = tr_x_yyyyy_yyzzz[i] * fe_0 + ts_yyyyy_xyyzzz[i] * fe_0 + tr_x_yyyyy_xyyzzz[i] * pa_x[i];

        tr_x_xyyyyy_xyzzzz[i] = tr_x_yyyyy_yzzzz[i] * fe_0 + ts_yyyyy_xyzzzz[i] * fe_0 + tr_x_yyyyy_xyzzzz[i] * pa_x[i];

        tr_x_xyyyyy_xzzzzz[i] = 4.0 * tr_x_xyyy_xzzzzz[i] * fe_0 + tr_x_xyyyy_xzzzzz[i] * pa_y[i];

        tr_x_xyyyyy_yyyyyy[i] = ts_yyyyy_yyyyyy[i] * fe_0 + tr_x_yyyyy_yyyyyy[i] * pa_x[i];

        tr_x_xyyyyy_yyyyyz[i] = ts_yyyyy_yyyyyz[i] * fe_0 + tr_x_yyyyy_yyyyyz[i] * pa_x[i];

        tr_x_xyyyyy_yyyyzz[i] = ts_yyyyy_yyyyzz[i] * fe_0 + tr_x_yyyyy_yyyyzz[i] * pa_x[i];

        tr_x_xyyyyy_yyyzzz[i] = ts_yyyyy_yyyzzz[i] * fe_0 + tr_x_yyyyy_yyyzzz[i] * pa_x[i];

        tr_x_xyyyyy_yyzzzz[i] = ts_yyyyy_yyzzzz[i] * fe_0 + tr_x_yyyyy_yyzzzz[i] * pa_x[i];

        tr_x_xyyyyy_yzzzzz[i] = ts_yyyyy_yzzzzz[i] * fe_0 + tr_x_yyyyy_yzzzzz[i] * pa_x[i];

        tr_x_xyyyyy_zzzzzz[i] = ts_yyyyy_zzzzzz[i] * fe_0 + tr_x_yyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 448-476 components of targeted buffer : II

    auto tr_x_xyyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 448);

    auto tr_x_xyyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 449);

    auto tr_x_xyyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 450);

    auto tr_x_xyyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 451);

    auto tr_x_xyyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 452);

    auto tr_x_xyyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 453);

    auto tr_x_xyyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 454);

    auto tr_x_xyyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 455);

    auto tr_x_xyyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 456);

    auto tr_x_xyyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 457);

    auto tr_x_xyyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 458);

    auto tr_x_xyyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 459);

    auto tr_x_xyyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 460);

    auto tr_x_xyyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 461);

    auto tr_x_xyyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 462);

    auto tr_x_xyyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 463);

    auto tr_x_xyyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 464);

    auto tr_x_xyyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 465);

    auto tr_x_xyyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 466);

    auto tr_x_xyyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 467);

    auto tr_x_xyyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 468);

    auto tr_x_xyyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 469);

    auto tr_x_xyyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 470);

    auto tr_x_xyyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 471);

    auto tr_x_xyyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 472);

    auto tr_x_xyyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 473);

    auto tr_x_xyyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 474);

    auto tr_x_xyyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 475);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyyy_xxxxxx, tr_x_xyyyy_xxxxxy, tr_x_xyyyy_xxxxy, tr_x_xyyyy_xxxxyy, tr_x_xyyyy_xxxxyz, tr_x_xyyyy_xxxyy, tr_x_xyyyy_xxxyyy, tr_x_xyyyy_xxxyyz, tr_x_xyyyy_xxxyz, tr_x_xyyyy_xxxyzz, tr_x_xyyyy_xxyyy, tr_x_xyyyy_xxyyyy, tr_x_xyyyy_xxyyyz, tr_x_xyyyy_xxyyz, tr_x_xyyyy_xxyyzz, tr_x_xyyyy_xxyzz, tr_x_xyyyy_xxyzzz, tr_x_xyyyy_xyyyy, tr_x_xyyyy_xyyyyy, tr_x_xyyyy_xyyyyz, tr_x_xyyyy_xyyyz, tr_x_xyyyy_xyyyzz, tr_x_xyyyy_xyyzz, tr_x_xyyyy_xyyzzz, tr_x_xyyyy_xyzzz, tr_x_xyyyy_xyzzzz, tr_x_xyyyy_yyyyyy, tr_x_xyyyyz_xxxxxx, tr_x_xyyyyz_xxxxxy, tr_x_xyyyyz_xxxxxz, tr_x_xyyyyz_xxxxyy, tr_x_xyyyyz_xxxxyz, tr_x_xyyyyz_xxxxzz, tr_x_xyyyyz_xxxyyy, tr_x_xyyyyz_xxxyyz, tr_x_xyyyyz_xxxyzz, tr_x_xyyyyz_xxxzzz, tr_x_xyyyyz_xxyyyy, tr_x_xyyyyz_xxyyyz, tr_x_xyyyyz_xxyyzz, tr_x_xyyyyz_xxyzzz, tr_x_xyyyyz_xxzzzz, tr_x_xyyyyz_xyyyyy, tr_x_xyyyyz_xyyyyz, tr_x_xyyyyz_xyyyzz, tr_x_xyyyyz_xyyzzz, tr_x_xyyyyz_xyzzzz, tr_x_xyyyyz_xzzzzz, tr_x_xyyyyz_yyyyyy, tr_x_xyyyyz_yyyyyz, tr_x_xyyyyz_yyyyzz, tr_x_xyyyyz_yyyzzz, tr_x_xyyyyz_yyzzzz, tr_x_xyyyyz_yzzzzz, tr_x_xyyyyz_zzzzzz, tr_x_xyyyz_xxxxxz, tr_x_xyyyz_xxxxzz, tr_x_xyyyz_xxxzzz, tr_x_xyyyz_xxzzzz, tr_x_xyyyz_xzzzzz, tr_x_xyyz_xxxxxz, tr_x_xyyz_xxxxzz, tr_x_xyyz_xxxzzz, tr_x_xyyz_xxzzzz, tr_x_xyyz_xzzzzz, tr_x_yyyyz_yyyyyz, tr_x_yyyyz_yyyyzz, tr_x_yyyyz_yyyzzz, tr_x_yyyyz_yyzzzz, tr_x_yyyyz_yzzzzz, tr_x_yyyyz_zzzzzz, ts_yyyyz_yyyyyz, ts_yyyyz_yyyyzz, ts_yyyyz_yyyzzz, ts_yyyyz_yyzzzz, ts_yyyyz_yzzzzz, ts_yyyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyz_xxxxxx[i] = tr_x_xyyyy_xxxxxx[i] * pa_z[i];

        tr_x_xyyyyz_xxxxxy[i] = tr_x_xyyyy_xxxxxy[i] * pa_z[i];

        tr_x_xyyyyz_xxxxxz[i] = 3.0 * tr_x_xyyz_xxxxxz[i] * fe_0 + tr_x_xyyyz_xxxxxz[i] * pa_y[i];

        tr_x_xyyyyz_xxxxyy[i] = tr_x_xyyyy_xxxxyy[i] * pa_z[i];

        tr_x_xyyyyz_xxxxyz[i] = tr_x_xyyyy_xxxxy[i] * fe_0 + tr_x_xyyyy_xxxxyz[i] * pa_z[i];

        tr_x_xyyyyz_xxxxzz[i] = 3.0 * tr_x_xyyz_xxxxzz[i] * fe_0 + tr_x_xyyyz_xxxxzz[i] * pa_y[i];

        tr_x_xyyyyz_xxxyyy[i] = tr_x_xyyyy_xxxyyy[i] * pa_z[i];

        tr_x_xyyyyz_xxxyyz[i] = tr_x_xyyyy_xxxyy[i] * fe_0 + tr_x_xyyyy_xxxyyz[i] * pa_z[i];

        tr_x_xyyyyz_xxxyzz[i] = 2.0 * tr_x_xyyyy_xxxyz[i] * fe_0 + tr_x_xyyyy_xxxyzz[i] * pa_z[i];

        tr_x_xyyyyz_xxxzzz[i] = 3.0 * tr_x_xyyz_xxxzzz[i] * fe_0 + tr_x_xyyyz_xxxzzz[i] * pa_y[i];

        tr_x_xyyyyz_xxyyyy[i] = tr_x_xyyyy_xxyyyy[i] * pa_z[i];

        tr_x_xyyyyz_xxyyyz[i] = tr_x_xyyyy_xxyyy[i] * fe_0 + tr_x_xyyyy_xxyyyz[i] * pa_z[i];

        tr_x_xyyyyz_xxyyzz[i] = 2.0 * tr_x_xyyyy_xxyyz[i] * fe_0 + tr_x_xyyyy_xxyyzz[i] * pa_z[i];

        tr_x_xyyyyz_xxyzzz[i] = 3.0 * tr_x_xyyyy_xxyzz[i] * fe_0 + tr_x_xyyyy_xxyzzz[i] * pa_z[i];

        tr_x_xyyyyz_xxzzzz[i] = 3.0 * tr_x_xyyz_xxzzzz[i] * fe_0 + tr_x_xyyyz_xxzzzz[i] * pa_y[i];

        tr_x_xyyyyz_xyyyyy[i] = tr_x_xyyyy_xyyyyy[i] * pa_z[i];

        tr_x_xyyyyz_xyyyyz[i] = tr_x_xyyyy_xyyyy[i] * fe_0 + tr_x_xyyyy_xyyyyz[i] * pa_z[i];

        tr_x_xyyyyz_xyyyzz[i] = 2.0 * tr_x_xyyyy_xyyyz[i] * fe_0 + tr_x_xyyyy_xyyyzz[i] * pa_z[i];

        tr_x_xyyyyz_xyyzzz[i] = 3.0 * tr_x_xyyyy_xyyzz[i] * fe_0 + tr_x_xyyyy_xyyzzz[i] * pa_z[i];

        tr_x_xyyyyz_xyzzzz[i] = 4.0 * tr_x_xyyyy_xyzzz[i] * fe_0 + tr_x_xyyyy_xyzzzz[i] * pa_z[i];

        tr_x_xyyyyz_xzzzzz[i] = 3.0 * tr_x_xyyz_xzzzzz[i] * fe_0 + tr_x_xyyyz_xzzzzz[i] * pa_y[i];

        tr_x_xyyyyz_yyyyyy[i] = tr_x_xyyyy_yyyyyy[i] * pa_z[i];

        tr_x_xyyyyz_yyyyyz[i] = ts_yyyyz_yyyyyz[i] * fe_0 + tr_x_yyyyz_yyyyyz[i] * pa_x[i];

        tr_x_xyyyyz_yyyyzz[i] = ts_yyyyz_yyyyzz[i] * fe_0 + tr_x_yyyyz_yyyyzz[i] * pa_x[i];

        tr_x_xyyyyz_yyyzzz[i] = ts_yyyyz_yyyzzz[i] * fe_0 + tr_x_yyyyz_yyyzzz[i] * pa_x[i];

        tr_x_xyyyyz_yyzzzz[i] = ts_yyyyz_yyzzzz[i] * fe_0 + tr_x_yyyyz_yyzzzz[i] * pa_x[i];

        tr_x_xyyyyz_yzzzzz[i] = ts_yyyyz_yzzzzz[i] * fe_0 + tr_x_yyyyz_yzzzzz[i] * pa_x[i];

        tr_x_xyyyyz_zzzzzz[i] = ts_yyyyz_zzzzzz[i] * fe_0 + tr_x_yyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 476-504 components of targeted buffer : II

    auto tr_x_xyyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 476);

    auto tr_x_xyyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 477);

    auto tr_x_xyyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 478);

    auto tr_x_xyyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 479);

    auto tr_x_xyyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 480);

    auto tr_x_xyyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 481);

    auto tr_x_xyyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 482);

    auto tr_x_xyyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 483);

    auto tr_x_xyyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 484);

    auto tr_x_xyyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 485);

    auto tr_x_xyyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 486);

    auto tr_x_xyyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 487);

    auto tr_x_xyyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 488);

    auto tr_x_xyyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 489);

    auto tr_x_xyyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 490);

    auto tr_x_xyyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 491);

    auto tr_x_xyyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 492);

    auto tr_x_xyyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 493);

    auto tr_x_xyyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 494);

    auto tr_x_xyyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 495);

    auto tr_x_xyyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 496);

    auto tr_x_xyyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 497);

    auto tr_x_xyyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 498);

    auto tr_x_xyyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 499);

    auto tr_x_xyyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 500);

    auto tr_x_xyyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 501);

    auto tr_x_xyyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 502);

    auto tr_x_xyyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 503);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyy_xxxxxy, tr_x_xyyy_xxxxyy, tr_x_xyyy_xxxyyy, tr_x_xyyy_xxyyyy, tr_x_xyyy_xyyyyy, tr_x_xyyyz_xxxxxy, tr_x_xyyyz_xxxxyy, tr_x_xyyyz_xxxyyy, tr_x_xyyyz_xxyyyy, tr_x_xyyyz_xyyyyy, tr_x_xyyyzz_xxxxxx, tr_x_xyyyzz_xxxxxy, tr_x_xyyyzz_xxxxxz, tr_x_xyyyzz_xxxxyy, tr_x_xyyyzz_xxxxyz, tr_x_xyyyzz_xxxxzz, tr_x_xyyyzz_xxxyyy, tr_x_xyyyzz_xxxyyz, tr_x_xyyyzz_xxxyzz, tr_x_xyyyzz_xxxzzz, tr_x_xyyyzz_xxyyyy, tr_x_xyyyzz_xxyyyz, tr_x_xyyyzz_xxyyzz, tr_x_xyyyzz_xxyzzz, tr_x_xyyyzz_xxzzzz, tr_x_xyyyzz_xyyyyy, tr_x_xyyyzz_xyyyyz, tr_x_xyyyzz_xyyyzz, tr_x_xyyyzz_xyyzzz, tr_x_xyyyzz_xyzzzz, tr_x_xyyyzz_xzzzzz, tr_x_xyyyzz_yyyyyy, tr_x_xyyyzz_yyyyyz, tr_x_xyyyzz_yyyyzz, tr_x_xyyyzz_yyyzzz, tr_x_xyyyzz_yyzzzz, tr_x_xyyyzz_yzzzzz, tr_x_xyyyzz_zzzzzz, tr_x_xyyzz_xxxxxx, tr_x_xyyzz_xxxxxz, tr_x_xyyzz_xxxxzz, tr_x_xyyzz_xxxzzz, tr_x_xyyzz_xxzzzz, tr_x_xyyzz_xzzzzz, tr_x_xyzz_xxxxxx, tr_x_xyzz_xxxxxz, tr_x_xyzz_xxxxzz, tr_x_xyzz_xxxzzz, tr_x_xyzz_xxzzzz, tr_x_xyzz_xzzzzz, tr_x_yyyzz_xxxxyz, tr_x_yyyzz_xxxyyz, tr_x_yyyzz_xxxyz, tr_x_yyyzz_xxxyzz, tr_x_yyyzz_xxyyyz, tr_x_yyyzz_xxyyz, tr_x_yyyzz_xxyyzz, tr_x_yyyzz_xxyzz, tr_x_yyyzz_xxyzzz, tr_x_yyyzz_xyyyyz, tr_x_yyyzz_xyyyz, tr_x_yyyzz_xyyyzz, tr_x_yyyzz_xyyzz, tr_x_yyyzz_xyyzzz, tr_x_yyyzz_xyzzz, tr_x_yyyzz_xyzzzz, tr_x_yyyzz_yyyyyy, tr_x_yyyzz_yyyyyz, tr_x_yyyzz_yyyyz, tr_x_yyyzz_yyyyzz, tr_x_yyyzz_yyyzz, tr_x_yyyzz_yyyzzz, tr_x_yyyzz_yyzzz, tr_x_yyyzz_yyzzzz, tr_x_yyyzz_yzzzz, tr_x_yyyzz_yzzzzz, tr_x_yyyzz_zzzzzz, ts_yyyzz_xxxxyz, ts_yyyzz_xxxyyz, ts_yyyzz_xxxyzz, ts_yyyzz_xxyyyz, ts_yyyzz_xxyyzz, ts_yyyzz_xxyzzz, ts_yyyzz_xyyyyz, ts_yyyzz_xyyyzz, ts_yyyzz_xyyzzz, ts_yyyzz_xyzzzz, ts_yyyzz_yyyyyy, ts_yyyzz_yyyyyz, ts_yyyzz_yyyyzz, ts_yyyzz_yyyzzz, ts_yyyzz_yyzzzz, ts_yyyzz_yzzzzz, ts_yyyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyzz_xxxxxx[i] = 2.0 * tr_x_xyzz_xxxxxx[i] * fe_0 + tr_x_xyyzz_xxxxxx[i] * pa_y[i];

        tr_x_xyyyzz_xxxxxy[i] = tr_x_xyyy_xxxxxy[i] * fe_0 + tr_x_xyyyz_xxxxxy[i] * pa_z[i];

        tr_x_xyyyzz_xxxxxz[i] = 2.0 * tr_x_xyzz_xxxxxz[i] * fe_0 + tr_x_xyyzz_xxxxxz[i] * pa_y[i];

        tr_x_xyyyzz_xxxxyy[i] = tr_x_xyyy_xxxxyy[i] * fe_0 + tr_x_xyyyz_xxxxyy[i] * pa_z[i];

        tr_x_xyyyzz_xxxxyz[i] = 4.0 * tr_x_yyyzz_xxxyz[i] * fe_0 + ts_yyyzz_xxxxyz[i] * fe_0 + tr_x_yyyzz_xxxxyz[i] * pa_x[i];

        tr_x_xyyyzz_xxxxzz[i] = 2.0 * tr_x_xyzz_xxxxzz[i] * fe_0 + tr_x_xyyzz_xxxxzz[i] * pa_y[i];

        tr_x_xyyyzz_xxxyyy[i] = tr_x_xyyy_xxxyyy[i] * fe_0 + tr_x_xyyyz_xxxyyy[i] * pa_z[i];

        tr_x_xyyyzz_xxxyyz[i] = 3.0 * tr_x_yyyzz_xxyyz[i] * fe_0 + ts_yyyzz_xxxyyz[i] * fe_0 + tr_x_yyyzz_xxxyyz[i] * pa_x[i];

        tr_x_xyyyzz_xxxyzz[i] = 3.0 * tr_x_yyyzz_xxyzz[i] * fe_0 + ts_yyyzz_xxxyzz[i] * fe_0 + tr_x_yyyzz_xxxyzz[i] * pa_x[i];

        tr_x_xyyyzz_xxxzzz[i] = 2.0 * tr_x_xyzz_xxxzzz[i] * fe_0 + tr_x_xyyzz_xxxzzz[i] * pa_y[i];

        tr_x_xyyyzz_xxyyyy[i] = tr_x_xyyy_xxyyyy[i] * fe_0 + tr_x_xyyyz_xxyyyy[i] * pa_z[i];

        tr_x_xyyyzz_xxyyyz[i] = 2.0 * tr_x_yyyzz_xyyyz[i] * fe_0 + ts_yyyzz_xxyyyz[i] * fe_0 + tr_x_yyyzz_xxyyyz[i] * pa_x[i];

        tr_x_xyyyzz_xxyyzz[i] = 2.0 * tr_x_yyyzz_xyyzz[i] * fe_0 + ts_yyyzz_xxyyzz[i] * fe_0 + tr_x_yyyzz_xxyyzz[i] * pa_x[i];

        tr_x_xyyyzz_xxyzzz[i] = 2.0 * tr_x_yyyzz_xyzzz[i] * fe_0 + ts_yyyzz_xxyzzz[i] * fe_0 + tr_x_yyyzz_xxyzzz[i] * pa_x[i];

        tr_x_xyyyzz_xxzzzz[i] = 2.0 * tr_x_xyzz_xxzzzz[i] * fe_0 + tr_x_xyyzz_xxzzzz[i] * pa_y[i];

        tr_x_xyyyzz_xyyyyy[i] = tr_x_xyyy_xyyyyy[i] * fe_0 + tr_x_xyyyz_xyyyyy[i] * pa_z[i];

        tr_x_xyyyzz_xyyyyz[i] = tr_x_yyyzz_yyyyz[i] * fe_0 + ts_yyyzz_xyyyyz[i] * fe_0 + tr_x_yyyzz_xyyyyz[i] * pa_x[i];

        tr_x_xyyyzz_xyyyzz[i] = tr_x_yyyzz_yyyzz[i] * fe_0 + ts_yyyzz_xyyyzz[i] * fe_0 + tr_x_yyyzz_xyyyzz[i] * pa_x[i];

        tr_x_xyyyzz_xyyzzz[i] = tr_x_yyyzz_yyzzz[i] * fe_0 + ts_yyyzz_xyyzzz[i] * fe_0 + tr_x_yyyzz_xyyzzz[i] * pa_x[i];

        tr_x_xyyyzz_xyzzzz[i] = tr_x_yyyzz_yzzzz[i] * fe_0 + ts_yyyzz_xyzzzz[i] * fe_0 + tr_x_yyyzz_xyzzzz[i] * pa_x[i];

        tr_x_xyyyzz_xzzzzz[i] = 2.0 * tr_x_xyzz_xzzzzz[i] * fe_0 + tr_x_xyyzz_xzzzzz[i] * pa_y[i];

        tr_x_xyyyzz_yyyyyy[i] = ts_yyyzz_yyyyyy[i] * fe_0 + tr_x_yyyzz_yyyyyy[i] * pa_x[i];

        tr_x_xyyyzz_yyyyyz[i] = ts_yyyzz_yyyyyz[i] * fe_0 + tr_x_yyyzz_yyyyyz[i] * pa_x[i];

        tr_x_xyyyzz_yyyyzz[i] = ts_yyyzz_yyyyzz[i] * fe_0 + tr_x_yyyzz_yyyyzz[i] * pa_x[i];

        tr_x_xyyyzz_yyyzzz[i] = ts_yyyzz_yyyzzz[i] * fe_0 + tr_x_yyyzz_yyyzzz[i] * pa_x[i];

        tr_x_xyyyzz_yyzzzz[i] = ts_yyyzz_yyzzzz[i] * fe_0 + tr_x_yyyzz_yyzzzz[i] * pa_x[i];

        tr_x_xyyyzz_yzzzzz[i] = ts_yyyzz_yzzzzz[i] * fe_0 + tr_x_yyyzz_yzzzzz[i] * pa_x[i];

        tr_x_xyyyzz_zzzzzz[i] = ts_yyyzz_zzzzzz[i] * fe_0 + tr_x_yyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 504-532 components of targeted buffer : II

    auto tr_x_xyyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 504);

    auto tr_x_xyyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 505);

    auto tr_x_xyyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 506);

    auto tr_x_xyyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 507);

    auto tr_x_xyyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 508);

    auto tr_x_xyyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 509);

    auto tr_x_xyyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 510);

    auto tr_x_xyyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 511);

    auto tr_x_xyyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 512);

    auto tr_x_xyyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 513);

    auto tr_x_xyyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 514);

    auto tr_x_xyyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 515);

    auto tr_x_xyyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 516);

    auto tr_x_xyyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 517);

    auto tr_x_xyyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 518);

    auto tr_x_xyyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 519);

    auto tr_x_xyyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 520);

    auto tr_x_xyyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 521);

    auto tr_x_xyyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 522);

    auto tr_x_xyyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 523);

    auto tr_x_xyyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 524);

    auto tr_x_xyyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 525);

    auto tr_x_xyyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 526);

    auto tr_x_xyyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 527);

    auto tr_x_xyyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 528);

    auto tr_x_xyyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 529);

    auto tr_x_xyyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 530);

    auto tr_x_xyyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 531);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyz_xxxxxy, tr_x_xyyz_xxxxyy, tr_x_xyyz_xxxyyy, tr_x_xyyz_xxyyyy, tr_x_xyyz_xyyyyy, tr_x_xyyzz_xxxxxy, tr_x_xyyzz_xxxxyy, tr_x_xyyzz_xxxyyy, tr_x_xyyzz_xxyyyy, tr_x_xyyzz_xyyyyy, tr_x_xyyzzz_xxxxxx, tr_x_xyyzzz_xxxxxy, tr_x_xyyzzz_xxxxxz, tr_x_xyyzzz_xxxxyy, tr_x_xyyzzz_xxxxyz, tr_x_xyyzzz_xxxxzz, tr_x_xyyzzz_xxxyyy, tr_x_xyyzzz_xxxyyz, tr_x_xyyzzz_xxxyzz, tr_x_xyyzzz_xxxzzz, tr_x_xyyzzz_xxyyyy, tr_x_xyyzzz_xxyyyz, tr_x_xyyzzz_xxyyzz, tr_x_xyyzzz_xxyzzz, tr_x_xyyzzz_xxzzzz, tr_x_xyyzzz_xyyyyy, tr_x_xyyzzz_xyyyyz, tr_x_xyyzzz_xyyyzz, tr_x_xyyzzz_xyyzzz, tr_x_xyyzzz_xyzzzz, tr_x_xyyzzz_xzzzzz, tr_x_xyyzzz_yyyyyy, tr_x_xyyzzz_yyyyyz, tr_x_xyyzzz_yyyyzz, tr_x_xyyzzz_yyyzzz, tr_x_xyyzzz_yyzzzz, tr_x_xyyzzz_yzzzzz, tr_x_xyyzzz_zzzzzz, tr_x_xyzzz_xxxxxx, tr_x_xyzzz_xxxxxz, tr_x_xyzzz_xxxxzz, tr_x_xyzzz_xxxzzz, tr_x_xyzzz_xxzzzz, tr_x_xyzzz_xzzzzz, tr_x_xzzz_xxxxxx, tr_x_xzzz_xxxxxz, tr_x_xzzz_xxxxzz, tr_x_xzzz_xxxzzz, tr_x_xzzz_xxzzzz, tr_x_xzzz_xzzzzz, tr_x_yyzzz_xxxxyz, tr_x_yyzzz_xxxyyz, tr_x_yyzzz_xxxyz, tr_x_yyzzz_xxxyzz, tr_x_yyzzz_xxyyyz, tr_x_yyzzz_xxyyz, tr_x_yyzzz_xxyyzz, tr_x_yyzzz_xxyzz, tr_x_yyzzz_xxyzzz, tr_x_yyzzz_xyyyyz, tr_x_yyzzz_xyyyz, tr_x_yyzzz_xyyyzz, tr_x_yyzzz_xyyzz, tr_x_yyzzz_xyyzzz, tr_x_yyzzz_xyzzz, tr_x_yyzzz_xyzzzz, tr_x_yyzzz_yyyyyy, tr_x_yyzzz_yyyyyz, tr_x_yyzzz_yyyyz, tr_x_yyzzz_yyyyzz, tr_x_yyzzz_yyyzz, tr_x_yyzzz_yyyzzz, tr_x_yyzzz_yyzzz, tr_x_yyzzz_yyzzzz, tr_x_yyzzz_yzzzz, tr_x_yyzzz_yzzzzz, tr_x_yyzzz_zzzzzz, ts_yyzzz_xxxxyz, ts_yyzzz_xxxyyz, ts_yyzzz_xxxyzz, ts_yyzzz_xxyyyz, ts_yyzzz_xxyyzz, ts_yyzzz_xxyzzz, ts_yyzzz_xyyyyz, ts_yyzzz_xyyyzz, ts_yyzzz_xyyzzz, ts_yyzzz_xyzzzz, ts_yyzzz_yyyyyy, ts_yyzzz_yyyyyz, ts_yyzzz_yyyyzz, ts_yyzzz_yyyzzz, ts_yyzzz_yyzzzz, ts_yyzzz_yzzzzz, ts_yyzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzzz_xxxxxx[i] = tr_x_xzzz_xxxxxx[i] * fe_0 + tr_x_xyzzz_xxxxxx[i] * pa_y[i];

        tr_x_xyyzzz_xxxxxy[i] = 2.0 * tr_x_xyyz_xxxxxy[i] * fe_0 + tr_x_xyyzz_xxxxxy[i] * pa_z[i];

        tr_x_xyyzzz_xxxxxz[i] = tr_x_xzzz_xxxxxz[i] * fe_0 + tr_x_xyzzz_xxxxxz[i] * pa_y[i];

        tr_x_xyyzzz_xxxxyy[i] = 2.0 * tr_x_xyyz_xxxxyy[i] * fe_0 + tr_x_xyyzz_xxxxyy[i] * pa_z[i];

        tr_x_xyyzzz_xxxxyz[i] = 4.0 * tr_x_yyzzz_xxxyz[i] * fe_0 + ts_yyzzz_xxxxyz[i] * fe_0 + tr_x_yyzzz_xxxxyz[i] * pa_x[i];

        tr_x_xyyzzz_xxxxzz[i] = tr_x_xzzz_xxxxzz[i] * fe_0 + tr_x_xyzzz_xxxxzz[i] * pa_y[i];

        tr_x_xyyzzz_xxxyyy[i] = 2.0 * tr_x_xyyz_xxxyyy[i] * fe_0 + tr_x_xyyzz_xxxyyy[i] * pa_z[i];

        tr_x_xyyzzz_xxxyyz[i] = 3.0 * tr_x_yyzzz_xxyyz[i] * fe_0 + ts_yyzzz_xxxyyz[i] * fe_0 + tr_x_yyzzz_xxxyyz[i] * pa_x[i];

        tr_x_xyyzzz_xxxyzz[i] = 3.0 * tr_x_yyzzz_xxyzz[i] * fe_0 + ts_yyzzz_xxxyzz[i] * fe_0 + tr_x_yyzzz_xxxyzz[i] * pa_x[i];

        tr_x_xyyzzz_xxxzzz[i] = tr_x_xzzz_xxxzzz[i] * fe_0 + tr_x_xyzzz_xxxzzz[i] * pa_y[i];

        tr_x_xyyzzz_xxyyyy[i] = 2.0 * tr_x_xyyz_xxyyyy[i] * fe_0 + tr_x_xyyzz_xxyyyy[i] * pa_z[i];

        tr_x_xyyzzz_xxyyyz[i] = 2.0 * tr_x_yyzzz_xyyyz[i] * fe_0 + ts_yyzzz_xxyyyz[i] * fe_0 + tr_x_yyzzz_xxyyyz[i] * pa_x[i];

        tr_x_xyyzzz_xxyyzz[i] = 2.0 * tr_x_yyzzz_xyyzz[i] * fe_0 + ts_yyzzz_xxyyzz[i] * fe_0 + tr_x_yyzzz_xxyyzz[i] * pa_x[i];

        tr_x_xyyzzz_xxyzzz[i] = 2.0 * tr_x_yyzzz_xyzzz[i] * fe_0 + ts_yyzzz_xxyzzz[i] * fe_0 + tr_x_yyzzz_xxyzzz[i] * pa_x[i];

        tr_x_xyyzzz_xxzzzz[i] = tr_x_xzzz_xxzzzz[i] * fe_0 + tr_x_xyzzz_xxzzzz[i] * pa_y[i];

        tr_x_xyyzzz_xyyyyy[i] = 2.0 * tr_x_xyyz_xyyyyy[i] * fe_0 + tr_x_xyyzz_xyyyyy[i] * pa_z[i];

        tr_x_xyyzzz_xyyyyz[i] = tr_x_yyzzz_yyyyz[i] * fe_0 + ts_yyzzz_xyyyyz[i] * fe_0 + tr_x_yyzzz_xyyyyz[i] * pa_x[i];

        tr_x_xyyzzz_xyyyzz[i] = tr_x_yyzzz_yyyzz[i] * fe_0 + ts_yyzzz_xyyyzz[i] * fe_0 + tr_x_yyzzz_xyyyzz[i] * pa_x[i];

        tr_x_xyyzzz_xyyzzz[i] = tr_x_yyzzz_yyzzz[i] * fe_0 + ts_yyzzz_xyyzzz[i] * fe_0 + tr_x_yyzzz_xyyzzz[i] * pa_x[i];

        tr_x_xyyzzz_xyzzzz[i] = tr_x_yyzzz_yzzzz[i] * fe_0 + ts_yyzzz_xyzzzz[i] * fe_0 + tr_x_yyzzz_xyzzzz[i] * pa_x[i];

        tr_x_xyyzzz_xzzzzz[i] = tr_x_xzzz_xzzzzz[i] * fe_0 + tr_x_xyzzz_xzzzzz[i] * pa_y[i];

        tr_x_xyyzzz_yyyyyy[i] = ts_yyzzz_yyyyyy[i] * fe_0 + tr_x_yyzzz_yyyyyy[i] * pa_x[i];

        tr_x_xyyzzz_yyyyyz[i] = ts_yyzzz_yyyyyz[i] * fe_0 + tr_x_yyzzz_yyyyyz[i] * pa_x[i];

        tr_x_xyyzzz_yyyyzz[i] = ts_yyzzz_yyyyzz[i] * fe_0 + tr_x_yyzzz_yyyyzz[i] * pa_x[i];

        tr_x_xyyzzz_yyyzzz[i] = ts_yyzzz_yyyzzz[i] * fe_0 + tr_x_yyzzz_yyyzzz[i] * pa_x[i];

        tr_x_xyyzzz_yyzzzz[i] = ts_yyzzz_yyzzzz[i] * fe_0 + tr_x_yyzzz_yyzzzz[i] * pa_x[i];

        tr_x_xyyzzz_yzzzzz[i] = ts_yyzzz_yzzzzz[i] * fe_0 + tr_x_yyzzz_yzzzzz[i] * pa_x[i];

        tr_x_xyyzzz_zzzzzz[i] = ts_yyzzz_zzzzzz[i] * fe_0 + tr_x_yyzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 532-560 components of targeted buffer : II

    auto tr_x_xyzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 532);

    auto tr_x_xyzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 533);

    auto tr_x_xyzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 534);

    auto tr_x_xyzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 535);

    auto tr_x_xyzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 536);

    auto tr_x_xyzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 537);

    auto tr_x_xyzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 538);

    auto tr_x_xyzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 539);

    auto tr_x_xyzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 540);

    auto tr_x_xyzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 541);

    auto tr_x_xyzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 542);

    auto tr_x_xyzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 543);

    auto tr_x_xyzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 544);

    auto tr_x_xyzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 545);

    auto tr_x_xyzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 546);

    auto tr_x_xyzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 547);

    auto tr_x_xyzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 548);

    auto tr_x_xyzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 549);

    auto tr_x_xyzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 550);

    auto tr_x_xyzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 551);

    auto tr_x_xyzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 552);

    auto tr_x_xyzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 553);

    auto tr_x_xyzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 554);

    auto tr_x_xyzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 555);

    auto tr_x_xyzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 556);

    auto tr_x_xyzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 557);

    auto tr_x_xyzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 558);

    auto tr_x_xyzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 559);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyzzzz_xxxxxx, tr_x_xyzzzz_xxxxxy, tr_x_xyzzzz_xxxxxz, tr_x_xyzzzz_xxxxyy, tr_x_xyzzzz_xxxxyz, tr_x_xyzzzz_xxxxzz, tr_x_xyzzzz_xxxyyy, tr_x_xyzzzz_xxxyyz, tr_x_xyzzzz_xxxyzz, tr_x_xyzzzz_xxxzzz, tr_x_xyzzzz_xxyyyy, tr_x_xyzzzz_xxyyyz, tr_x_xyzzzz_xxyyzz, tr_x_xyzzzz_xxyzzz, tr_x_xyzzzz_xxzzzz, tr_x_xyzzzz_xyyyyy, tr_x_xyzzzz_xyyyyz, tr_x_xyzzzz_xyyyzz, tr_x_xyzzzz_xyyzzz, tr_x_xyzzzz_xyzzzz, tr_x_xyzzzz_xzzzzz, tr_x_xyzzzz_yyyyyy, tr_x_xyzzzz_yyyyyz, tr_x_xyzzzz_yyyyzz, tr_x_xyzzzz_yyyzzz, tr_x_xyzzzz_yyzzzz, tr_x_xyzzzz_yzzzzz, tr_x_xyzzzz_zzzzzz, tr_x_xzzzz_xxxxx, tr_x_xzzzz_xxxxxx, tr_x_xzzzz_xxxxxy, tr_x_xzzzz_xxxxxz, tr_x_xzzzz_xxxxy, tr_x_xzzzz_xxxxyy, tr_x_xzzzz_xxxxyz, tr_x_xzzzz_xxxxz, tr_x_xzzzz_xxxxzz, tr_x_xzzzz_xxxyy, tr_x_xzzzz_xxxyyy, tr_x_xzzzz_xxxyyz, tr_x_xzzzz_xxxyz, tr_x_xzzzz_xxxyzz, tr_x_xzzzz_xxxzz, tr_x_xzzzz_xxxzzz, tr_x_xzzzz_xxyyy, tr_x_xzzzz_xxyyyy, tr_x_xzzzz_xxyyyz, tr_x_xzzzz_xxyyz, tr_x_xzzzz_xxyyzz, tr_x_xzzzz_xxyzz, tr_x_xzzzz_xxyzzz, tr_x_xzzzz_xxzzz, tr_x_xzzzz_xxzzzz, tr_x_xzzzz_xyyyy, tr_x_xzzzz_xyyyyy, tr_x_xzzzz_xyyyyz, tr_x_xzzzz_xyyyz, tr_x_xzzzz_xyyyzz, tr_x_xzzzz_xyyzz, tr_x_xzzzz_xyyzzz, tr_x_xzzzz_xyzzz, tr_x_xzzzz_xyzzzz, tr_x_xzzzz_xzzzz, tr_x_xzzzz_xzzzzz, tr_x_xzzzz_zzzzzz, tr_x_yzzzz_yyyyyy, tr_x_yzzzz_yyyyyz, tr_x_yzzzz_yyyyzz, tr_x_yzzzz_yyyzzz, tr_x_yzzzz_yyzzzz, tr_x_yzzzz_yzzzzz, ts_yzzzz_yyyyyy, ts_yzzzz_yyyyyz, ts_yzzzz_yyyyzz, ts_yzzzz_yyyzzz, ts_yzzzz_yyzzzz, ts_yzzzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzzz_xxxxxx[i] = tr_x_xzzzz_xxxxxx[i] * pa_y[i];

        tr_x_xyzzzz_xxxxxy[i] = tr_x_xzzzz_xxxxx[i] * fe_0 + tr_x_xzzzz_xxxxxy[i] * pa_y[i];

        tr_x_xyzzzz_xxxxxz[i] = tr_x_xzzzz_xxxxxz[i] * pa_y[i];

        tr_x_xyzzzz_xxxxyy[i] = 2.0 * tr_x_xzzzz_xxxxy[i] * fe_0 + tr_x_xzzzz_xxxxyy[i] * pa_y[i];

        tr_x_xyzzzz_xxxxyz[i] = tr_x_xzzzz_xxxxz[i] * fe_0 + tr_x_xzzzz_xxxxyz[i] * pa_y[i];

        tr_x_xyzzzz_xxxxzz[i] = tr_x_xzzzz_xxxxzz[i] * pa_y[i];

        tr_x_xyzzzz_xxxyyy[i] = 3.0 * tr_x_xzzzz_xxxyy[i] * fe_0 + tr_x_xzzzz_xxxyyy[i] * pa_y[i];

        tr_x_xyzzzz_xxxyyz[i] = 2.0 * tr_x_xzzzz_xxxyz[i] * fe_0 + tr_x_xzzzz_xxxyyz[i] * pa_y[i];

        tr_x_xyzzzz_xxxyzz[i] = tr_x_xzzzz_xxxzz[i] * fe_0 + tr_x_xzzzz_xxxyzz[i] * pa_y[i];

        tr_x_xyzzzz_xxxzzz[i] = tr_x_xzzzz_xxxzzz[i] * pa_y[i];

        tr_x_xyzzzz_xxyyyy[i] = 4.0 * tr_x_xzzzz_xxyyy[i] * fe_0 + tr_x_xzzzz_xxyyyy[i] * pa_y[i];

        tr_x_xyzzzz_xxyyyz[i] = 3.0 * tr_x_xzzzz_xxyyz[i] * fe_0 + tr_x_xzzzz_xxyyyz[i] * pa_y[i];

        tr_x_xyzzzz_xxyyzz[i] = 2.0 * tr_x_xzzzz_xxyzz[i] * fe_0 + tr_x_xzzzz_xxyyzz[i] * pa_y[i];

        tr_x_xyzzzz_xxyzzz[i] = tr_x_xzzzz_xxzzz[i] * fe_0 + tr_x_xzzzz_xxyzzz[i] * pa_y[i];

        tr_x_xyzzzz_xxzzzz[i] = tr_x_xzzzz_xxzzzz[i] * pa_y[i];

        tr_x_xyzzzz_xyyyyy[i] = 5.0 * tr_x_xzzzz_xyyyy[i] * fe_0 + tr_x_xzzzz_xyyyyy[i] * pa_y[i];

        tr_x_xyzzzz_xyyyyz[i] = 4.0 * tr_x_xzzzz_xyyyz[i] * fe_0 + tr_x_xzzzz_xyyyyz[i] * pa_y[i];

        tr_x_xyzzzz_xyyyzz[i] = 3.0 * tr_x_xzzzz_xyyzz[i] * fe_0 + tr_x_xzzzz_xyyyzz[i] * pa_y[i];

        tr_x_xyzzzz_xyyzzz[i] = 2.0 * tr_x_xzzzz_xyzzz[i] * fe_0 + tr_x_xzzzz_xyyzzz[i] * pa_y[i];

        tr_x_xyzzzz_xyzzzz[i] = tr_x_xzzzz_xzzzz[i] * fe_0 + tr_x_xzzzz_xyzzzz[i] * pa_y[i];

        tr_x_xyzzzz_xzzzzz[i] = tr_x_xzzzz_xzzzzz[i] * pa_y[i];

        tr_x_xyzzzz_yyyyyy[i] = ts_yzzzz_yyyyyy[i] * fe_0 + tr_x_yzzzz_yyyyyy[i] * pa_x[i];

        tr_x_xyzzzz_yyyyyz[i] = ts_yzzzz_yyyyyz[i] * fe_0 + tr_x_yzzzz_yyyyyz[i] * pa_x[i];

        tr_x_xyzzzz_yyyyzz[i] = ts_yzzzz_yyyyzz[i] * fe_0 + tr_x_yzzzz_yyyyzz[i] * pa_x[i];

        tr_x_xyzzzz_yyyzzz[i] = ts_yzzzz_yyyzzz[i] * fe_0 + tr_x_yzzzz_yyyzzz[i] * pa_x[i];

        tr_x_xyzzzz_yyzzzz[i] = ts_yzzzz_yyzzzz[i] * fe_0 + tr_x_yzzzz_yyzzzz[i] * pa_x[i];

        tr_x_xyzzzz_yzzzzz[i] = ts_yzzzz_yzzzzz[i] * fe_0 + tr_x_yzzzz_yzzzzz[i] * pa_x[i];

        tr_x_xyzzzz_zzzzzz[i] = tr_x_xzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 560-588 components of targeted buffer : II

    auto tr_x_xzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 560);

    auto tr_x_xzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 561);

    auto tr_x_xzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 562);

    auto tr_x_xzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 563);

    auto tr_x_xzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 564);

    auto tr_x_xzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 565);

    auto tr_x_xzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 566);

    auto tr_x_xzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 567);

    auto tr_x_xzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 568);

    auto tr_x_xzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 569);

    auto tr_x_xzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 570);

    auto tr_x_xzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 571);

    auto tr_x_xzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 572);

    auto tr_x_xzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 573);

    auto tr_x_xzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 574);

    auto tr_x_xzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 575);

    auto tr_x_xzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 576);

    auto tr_x_xzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 577);

    auto tr_x_xzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 578);

    auto tr_x_xzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 579);

    auto tr_x_xzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 580);

    auto tr_x_xzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 581);

    auto tr_x_xzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 582);

    auto tr_x_xzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 583);

    auto tr_x_xzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 584);

    auto tr_x_xzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 585);

    auto tr_x_xzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 586);

    auto tr_x_xzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 587);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xzzz_xxxxxx, tr_x_xzzz_xxxxxy, tr_x_xzzz_xxxxyy, tr_x_xzzz_xxxyyy, tr_x_xzzz_xxyyyy, tr_x_xzzz_xyyyyy, tr_x_xzzzz_xxxxxx, tr_x_xzzzz_xxxxxy, tr_x_xzzzz_xxxxyy, tr_x_xzzzz_xxxyyy, tr_x_xzzzz_xxyyyy, tr_x_xzzzz_xyyyyy, tr_x_xzzzzz_xxxxxx, tr_x_xzzzzz_xxxxxy, tr_x_xzzzzz_xxxxxz, tr_x_xzzzzz_xxxxyy, tr_x_xzzzzz_xxxxyz, tr_x_xzzzzz_xxxxzz, tr_x_xzzzzz_xxxyyy, tr_x_xzzzzz_xxxyyz, tr_x_xzzzzz_xxxyzz, tr_x_xzzzzz_xxxzzz, tr_x_xzzzzz_xxyyyy, tr_x_xzzzzz_xxyyyz, tr_x_xzzzzz_xxyyzz, tr_x_xzzzzz_xxyzzz, tr_x_xzzzzz_xxzzzz, tr_x_xzzzzz_xyyyyy, tr_x_xzzzzz_xyyyyz, tr_x_xzzzzz_xyyyzz, tr_x_xzzzzz_xyyzzz, tr_x_xzzzzz_xyzzzz, tr_x_xzzzzz_xzzzzz, tr_x_xzzzzz_yyyyyy, tr_x_xzzzzz_yyyyyz, tr_x_xzzzzz_yyyyzz, tr_x_xzzzzz_yyyzzz, tr_x_xzzzzz_yyzzzz, tr_x_xzzzzz_yzzzzz, tr_x_xzzzzz_zzzzzz, tr_x_zzzzz_xxxxxz, tr_x_zzzzz_xxxxyz, tr_x_zzzzz_xxxxz, tr_x_zzzzz_xxxxzz, tr_x_zzzzz_xxxyyz, tr_x_zzzzz_xxxyz, tr_x_zzzzz_xxxyzz, tr_x_zzzzz_xxxzz, tr_x_zzzzz_xxxzzz, tr_x_zzzzz_xxyyyz, tr_x_zzzzz_xxyyz, tr_x_zzzzz_xxyyzz, tr_x_zzzzz_xxyzz, tr_x_zzzzz_xxyzzz, tr_x_zzzzz_xxzzz, tr_x_zzzzz_xxzzzz, tr_x_zzzzz_xyyyyz, tr_x_zzzzz_xyyyz, tr_x_zzzzz_xyyyzz, tr_x_zzzzz_xyyzz, tr_x_zzzzz_xyyzzz, tr_x_zzzzz_xyzzz, tr_x_zzzzz_xyzzzz, tr_x_zzzzz_xzzzz, tr_x_zzzzz_xzzzzz, tr_x_zzzzz_yyyyyy, tr_x_zzzzz_yyyyyz, tr_x_zzzzz_yyyyz, tr_x_zzzzz_yyyyzz, tr_x_zzzzz_yyyzz, tr_x_zzzzz_yyyzzz, tr_x_zzzzz_yyzzz, tr_x_zzzzz_yyzzzz, tr_x_zzzzz_yzzzz, tr_x_zzzzz_yzzzzz, tr_x_zzzzz_zzzzz, tr_x_zzzzz_zzzzzz, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyyy, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzzz_xxxxxx[i] = 4.0 * tr_x_xzzz_xxxxxx[i] * fe_0 + tr_x_xzzzz_xxxxxx[i] * pa_z[i];

        tr_x_xzzzzz_xxxxxy[i] = 4.0 * tr_x_xzzz_xxxxxy[i] * fe_0 + tr_x_xzzzz_xxxxxy[i] * pa_z[i];

        tr_x_xzzzzz_xxxxxz[i] = 5.0 * tr_x_zzzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxxxz[i] * fe_0 + tr_x_zzzzz_xxxxxz[i] * pa_x[i];

        tr_x_xzzzzz_xxxxyy[i] = 4.0 * tr_x_xzzz_xxxxyy[i] * fe_0 + tr_x_xzzzz_xxxxyy[i] * pa_z[i];

        tr_x_xzzzzz_xxxxyz[i] = 4.0 * tr_x_zzzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxxyz[i] * fe_0 + tr_x_zzzzz_xxxxyz[i] * pa_x[i];

        tr_x_xzzzzz_xxxxzz[i] = 4.0 * tr_x_zzzzz_xxxzz[i] * fe_0 + ts_zzzzz_xxxxzz[i] * fe_0 + tr_x_zzzzz_xxxxzz[i] * pa_x[i];

        tr_x_xzzzzz_xxxyyy[i] = 4.0 * tr_x_xzzz_xxxyyy[i] * fe_0 + tr_x_xzzzz_xxxyyy[i] * pa_z[i];

        tr_x_xzzzzz_xxxyyz[i] = 3.0 * tr_x_zzzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxxyyz[i] * fe_0 + tr_x_zzzzz_xxxyyz[i] * pa_x[i];

        tr_x_xzzzzz_xxxyzz[i] = 3.0 * tr_x_zzzzz_xxyzz[i] * fe_0 + ts_zzzzz_xxxyzz[i] * fe_0 + tr_x_zzzzz_xxxyzz[i] * pa_x[i];

        tr_x_xzzzzz_xxxzzz[i] = 3.0 * tr_x_zzzzz_xxzzz[i] * fe_0 + ts_zzzzz_xxxzzz[i] * fe_0 + tr_x_zzzzz_xxxzzz[i] * pa_x[i];

        tr_x_xzzzzz_xxyyyy[i] = 4.0 * tr_x_xzzz_xxyyyy[i] * fe_0 + tr_x_xzzzz_xxyyyy[i] * pa_z[i];

        tr_x_xzzzzz_xxyyyz[i] = 2.0 * tr_x_zzzzz_xyyyz[i] * fe_0 + ts_zzzzz_xxyyyz[i] * fe_0 + tr_x_zzzzz_xxyyyz[i] * pa_x[i];

        tr_x_xzzzzz_xxyyzz[i] = 2.0 * tr_x_zzzzz_xyyzz[i] * fe_0 + ts_zzzzz_xxyyzz[i] * fe_0 + tr_x_zzzzz_xxyyzz[i] * pa_x[i];

        tr_x_xzzzzz_xxyzzz[i] = 2.0 * tr_x_zzzzz_xyzzz[i] * fe_0 + ts_zzzzz_xxyzzz[i] * fe_0 + tr_x_zzzzz_xxyzzz[i] * pa_x[i];

        tr_x_xzzzzz_xxzzzz[i] = 2.0 * tr_x_zzzzz_xzzzz[i] * fe_0 + ts_zzzzz_xxzzzz[i] * fe_0 + tr_x_zzzzz_xxzzzz[i] * pa_x[i];

        tr_x_xzzzzz_xyyyyy[i] = 4.0 * tr_x_xzzz_xyyyyy[i] * fe_0 + tr_x_xzzzz_xyyyyy[i] * pa_z[i];

        tr_x_xzzzzz_xyyyyz[i] = tr_x_zzzzz_yyyyz[i] * fe_0 + ts_zzzzz_xyyyyz[i] * fe_0 + tr_x_zzzzz_xyyyyz[i] * pa_x[i];

        tr_x_xzzzzz_xyyyzz[i] = tr_x_zzzzz_yyyzz[i] * fe_0 + ts_zzzzz_xyyyzz[i] * fe_0 + tr_x_zzzzz_xyyyzz[i] * pa_x[i];

        tr_x_xzzzzz_xyyzzz[i] = tr_x_zzzzz_yyzzz[i] * fe_0 + ts_zzzzz_xyyzzz[i] * fe_0 + tr_x_zzzzz_xyyzzz[i] * pa_x[i];

        tr_x_xzzzzz_xyzzzz[i] = tr_x_zzzzz_yzzzz[i] * fe_0 + ts_zzzzz_xyzzzz[i] * fe_0 + tr_x_zzzzz_xyzzzz[i] * pa_x[i];

        tr_x_xzzzzz_xzzzzz[i] = tr_x_zzzzz_zzzzz[i] * fe_0 + ts_zzzzz_xzzzzz[i] * fe_0 + tr_x_zzzzz_xzzzzz[i] * pa_x[i];

        tr_x_xzzzzz_yyyyyy[i] = ts_zzzzz_yyyyyy[i] * fe_0 + tr_x_zzzzz_yyyyyy[i] * pa_x[i];

        tr_x_xzzzzz_yyyyyz[i] = ts_zzzzz_yyyyyz[i] * fe_0 + tr_x_zzzzz_yyyyyz[i] * pa_x[i];

        tr_x_xzzzzz_yyyyzz[i] = ts_zzzzz_yyyyzz[i] * fe_0 + tr_x_zzzzz_yyyyzz[i] * pa_x[i];

        tr_x_xzzzzz_yyyzzz[i] = ts_zzzzz_yyyzzz[i] * fe_0 + tr_x_zzzzz_yyyzzz[i] * pa_x[i];

        tr_x_xzzzzz_yyzzzz[i] = ts_zzzzz_yyzzzz[i] * fe_0 + tr_x_zzzzz_yyzzzz[i] * pa_x[i];

        tr_x_xzzzzz_yzzzzz[i] = ts_zzzzz_yzzzzz[i] * fe_0 + tr_x_zzzzz_yzzzzz[i] * pa_x[i];

        tr_x_xzzzzz_zzzzzz[i] = ts_zzzzz_zzzzzz[i] * fe_0 + tr_x_zzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 588-616 components of targeted buffer : II

    auto tr_x_yyyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 588);

    auto tr_x_yyyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 589);

    auto tr_x_yyyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 590);

    auto tr_x_yyyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 591);

    auto tr_x_yyyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 592);

    auto tr_x_yyyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 593);

    auto tr_x_yyyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 594);

    auto tr_x_yyyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 595);

    auto tr_x_yyyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 596);

    auto tr_x_yyyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 597);

    auto tr_x_yyyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 598);

    auto tr_x_yyyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 599);

    auto tr_x_yyyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 600);

    auto tr_x_yyyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 601);

    auto tr_x_yyyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 602);

    auto tr_x_yyyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 603);

    auto tr_x_yyyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 604);

    auto tr_x_yyyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 605);

    auto tr_x_yyyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 606);

    auto tr_x_yyyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 607);

    auto tr_x_yyyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 608);

    auto tr_x_yyyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 609);

    auto tr_x_yyyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 610);

    auto tr_x_yyyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 611);

    auto tr_x_yyyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 612);

    auto tr_x_yyyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 613);

    auto tr_x_yyyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 614);

    auto tr_x_yyyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 615);

    #pragma omp simd aligned(pa_y, tr_x_yyyy_xxxxxx, tr_x_yyyy_xxxxxy, tr_x_yyyy_xxxxxz, tr_x_yyyy_xxxxyy, tr_x_yyyy_xxxxyz, tr_x_yyyy_xxxxzz, tr_x_yyyy_xxxyyy, tr_x_yyyy_xxxyyz, tr_x_yyyy_xxxyzz, tr_x_yyyy_xxxzzz, tr_x_yyyy_xxyyyy, tr_x_yyyy_xxyyyz, tr_x_yyyy_xxyyzz, tr_x_yyyy_xxyzzz, tr_x_yyyy_xxzzzz, tr_x_yyyy_xyyyyy, tr_x_yyyy_xyyyyz, tr_x_yyyy_xyyyzz, tr_x_yyyy_xyyzzz, tr_x_yyyy_xyzzzz, tr_x_yyyy_xzzzzz, tr_x_yyyy_yyyyyy, tr_x_yyyy_yyyyyz, tr_x_yyyy_yyyyzz, tr_x_yyyy_yyyzzz, tr_x_yyyy_yyzzzz, tr_x_yyyy_yzzzzz, tr_x_yyyy_zzzzzz, tr_x_yyyyy_xxxxx, tr_x_yyyyy_xxxxxx, tr_x_yyyyy_xxxxxy, tr_x_yyyyy_xxxxxz, tr_x_yyyyy_xxxxy, tr_x_yyyyy_xxxxyy, tr_x_yyyyy_xxxxyz, tr_x_yyyyy_xxxxz, tr_x_yyyyy_xxxxzz, tr_x_yyyyy_xxxyy, tr_x_yyyyy_xxxyyy, tr_x_yyyyy_xxxyyz, tr_x_yyyyy_xxxyz, tr_x_yyyyy_xxxyzz, tr_x_yyyyy_xxxzz, tr_x_yyyyy_xxxzzz, tr_x_yyyyy_xxyyy, tr_x_yyyyy_xxyyyy, tr_x_yyyyy_xxyyyz, tr_x_yyyyy_xxyyz, tr_x_yyyyy_xxyyzz, tr_x_yyyyy_xxyzz, tr_x_yyyyy_xxyzzz, tr_x_yyyyy_xxzzz, tr_x_yyyyy_xxzzzz, tr_x_yyyyy_xyyyy, tr_x_yyyyy_xyyyyy, tr_x_yyyyy_xyyyyz, tr_x_yyyyy_xyyyz, tr_x_yyyyy_xyyyzz, tr_x_yyyyy_xyyzz, tr_x_yyyyy_xyyzzz, tr_x_yyyyy_xyzzz, tr_x_yyyyy_xyzzzz, tr_x_yyyyy_xzzzz, tr_x_yyyyy_xzzzzz, tr_x_yyyyy_yyyyy, tr_x_yyyyy_yyyyyy, tr_x_yyyyy_yyyyyz, tr_x_yyyyy_yyyyz, tr_x_yyyyy_yyyyzz, tr_x_yyyyy_yyyzz, tr_x_yyyyy_yyyzzz, tr_x_yyyyy_yyzzz, tr_x_yyyyy_yyzzzz, tr_x_yyyyy_yzzzz, tr_x_yyyyy_yzzzzz, tr_x_yyyyy_zzzzz, tr_x_yyyyy_zzzzzz, tr_x_yyyyyy_xxxxxx, tr_x_yyyyyy_xxxxxy, tr_x_yyyyyy_xxxxxz, tr_x_yyyyyy_xxxxyy, tr_x_yyyyyy_xxxxyz, tr_x_yyyyyy_xxxxzz, tr_x_yyyyyy_xxxyyy, tr_x_yyyyyy_xxxyyz, tr_x_yyyyyy_xxxyzz, tr_x_yyyyyy_xxxzzz, tr_x_yyyyyy_xxyyyy, tr_x_yyyyyy_xxyyyz, tr_x_yyyyyy_xxyyzz, tr_x_yyyyyy_xxyzzz, tr_x_yyyyyy_xxzzzz, tr_x_yyyyyy_xyyyyy, tr_x_yyyyyy_xyyyyz, tr_x_yyyyyy_xyyyzz, tr_x_yyyyyy_xyyzzz, tr_x_yyyyyy_xyzzzz, tr_x_yyyyyy_xzzzzz, tr_x_yyyyyy_yyyyyy, tr_x_yyyyyy_yyyyyz, tr_x_yyyyyy_yyyyzz, tr_x_yyyyyy_yyyzzz, tr_x_yyyyyy_yyzzzz, tr_x_yyyyyy_yzzzzz, tr_x_yyyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyy_xxxxxx[i] = 5.0 * tr_x_yyyy_xxxxxx[i] * fe_0 + tr_x_yyyyy_xxxxxx[i] * pa_y[i];

        tr_x_yyyyyy_xxxxxy[i] = 5.0 * tr_x_yyyy_xxxxxy[i] * fe_0 + tr_x_yyyyy_xxxxx[i] * fe_0 + tr_x_yyyyy_xxxxxy[i] * pa_y[i];

        tr_x_yyyyyy_xxxxxz[i] = 5.0 * tr_x_yyyy_xxxxxz[i] * fe_0 + tr_x_yyyyy_xxxxxz[i] * pa_y[i];

        tr_x_yyyyyy_xxxxyy[i] = 5.0 * tr_x_yyyy_xxxxyy[i] * fe_0 + 2.0 * tr_x_yyyyy_xxxxy[i] * fe_0 + tr_x_yyyyy_xxxxyy[i] * pa_y[i];

        tr_x_yyyyyy_xxxxyz[i] = 5.0 * tr_x_yyyy_xxxxyz[i] * fe_0 + tr_x_yyyyy_xxxxz[i] * fe_0 + tr_x_yyyyy_xxxxyz[i] * pa_y[i];

        tr_x_yyyyyy_xxxxzz[i] = 5.0 * tr_x_yyyy_xxxxzz[i] * fe_0 + tr_x_yyyyy_xxxxzz[i] * pa_y[i];

        tr_x_yyyyyy_xxxyyy[i] = 5.0 * tr_x_yyyy_xxxyyy[i] * fe_0 + 3.0 * tr_x_yyyyy_xxxyy[i] * fe_0 + tr_x_yyyyy_xxxyyy[i] * pa_y[i];

        tr_x_yyyyyy_xxxyyz[i] = 5.0 * tr_x_yyyy_xxxyyz[i] * fe_0 + 2.0 * tr_x_yyyyy_xxxyz[i] * fe_0 + tr_x_yyyyy_xxxyyz[i] * pa_y[i];

        tr_x_yyyyyy_xxxyzz[i] = 5.0 * tr_x_yyyy_xxxyzz[i] * fe_0 + tr_x_yyyyy_xxxzz[i] * fe_0 + tr_x_yyyyy_xxxyzz[i] * pa_y[i];

        tr_x_yyyyyy_xxxzzz[i] = 5.0 * tr_x_yyyy_xxxzzz[i] * fe_0 + tr_x_yyyyy_xxxzzz[i] * pa_y[i];

        tr_x_yyyyyy_xxyyyy[i] = 5.0 * tr_x_yyyy_xxyyyy[i] * fe_0 + 4.0 * tr_x_yyyyy_xxyyy[i] * fe_0 + tr_x_yyyyy_xxyyyy[i] * pa_y[i];

        tr_x_yyyyyy_xxyyyz[i] = 5.0 * tr_x_yyyy_xxyyyz[i] * fe_0 + 3.0 * tr_x_yyyyy_xxyyz[i] * fe_0 + tr_x_yyyyy_xxyyyz[i] * pa_y[i];

        tr_x_yyyyyy_xxyyzz[i] = 5.0 * tr_x_yyyy_xxyyzz[i] * fe_0 + 2.0 * tr_x_yyyyy_xxyzz[i] * fe_0 + tr_x_yyyyy_xxyyzz[i] * pa_y[i];

        tr_x_yyyyyy_xxyzzz[i] = 5.0 * tr_x_yyyy_xxyzzz[i] * fe_0 + tr_x_yyyyy_xxzzz[i] * fe_0 + tr_x_yyyyy_xxyzzz[i] * pa_y[i];

        tr_x_yyyyyy_xxzzzz[i] = 5.0 * tr_x_yyyy_xxzzzz[i] * fe_0 + tr_x_yyyyy_xxzzzz[i] * pa_y[i];

        tr_x_yyyyyy_xyyyyy[i] = 5.0 * tr_x_yyyy_xyyyyy[i] * fe_0 + 5.0 * tr_x_yyyyy_xyyyy[i] * fe_0 + tr_x_yyyyy_xyyyyy[i] * pa_y[i];

        tr_x_yyyyyy_xyyyyz[i] = 5.0 * tr_x_yyyy_xyyyyz[i] * fe_0 + 4.0 * tr_x_yyyyy_xyyyz[i] * fe_0 + tr_x_yyyyy_xyyyyz[i] * pa_y[i];

        tr_x_yyyyyy_xyyyzz[i] = 5.0 * tr_x_yyyy_xyyyzz[i] * fe_0 + 3.0 * tr_x_yyyyy_xyyzz[i] * fe_0 + tr_x_yyyyy_xyyyzz[i] * pa_y[i];

        tr_x_yyyyyy_xyyzzz[i] = 5.0 * tr_x_yyyy_xyyzzz[i] * fe_0 + 2.0 * tr_x_yyyyy_xyzzz[i] * fe_0 + tr_x_yyyyy_xyyzzz[i] * pa_y[i];

        tr_x_yyyyyy_xyzzzz[i] = 5.0 * tr_x_yyyy_xyzzzz[i] * fe_0 + tr_x_yyyyy_xzzzz[i] * fe_0 + tr_x_yyyyy_xyzzzz[i] * pa_y[i];

        tr_x_yyyyyy_xzzzzz[i] = 5.0 * tr_x_yyyy_xzzzzz[i] * fe_0 + tr_x_yyyyy_xzzzzz[i] * pa_y[i];

        tr_x_yyyyyy_yyyyyy[i] = 5.0 * tr_x_yyyy_yyyyyy[i] * fe_0 + 6.0 * tr_x_yyyyy_yyyyy[i] * fe_0 + tr_x_yyyyy_yyyyyy[i] * pa_y[i];

        tr_x_yyyyyy_yyyyyz[i] = 5.0 * tr_x_yyyy_yyyyyz[i] * fe_0 + 5.0 * tr_x_yyyyy_yyyyz[i] * fe_0 + tr_x_yyyyy_yyyyyz[i] * pa_y[i];

        tr_x_yyyyyy_yyyyzz[i] = 5.0 * tr_x_yyyy_yyyyzz[i] * fe_0 + 4.0 * tr_x_yyyyy_yyyzz[i] * fe_0 + tr_x_yyyyy_yyyyzz[i] * pa_y[i];

        tr_x_yyyyyy_yyyzzz[i] = 5.0 * tr_x_yyyy_yyyzzz[i] * fe_0 + 3.0 * tr_x_yyyyy_yyzzz[i] * fe_0 + tr_x_yyyyy_yyyzzz[i] * pa_y[i];

        tr_x_yyyyyy_yyzzzz[i] = 5.0 * tr_x_yyyy_yyzzzz[i] * fe_0 + 2.0 * tr_x_yyyyy_yzzzz[i] * fe_0 + tr_x_yyyyy_yyzzzz[i] * pa_y[i];

        tr_x_yyyyyy_yzzzzz[i] = 5.0 * tr_x_yyyy_yzzzzz[i] * fe_0 + tr_x_yyyyy_zzzzz[i] * fe_0 + tr_x_yyyyy_yzzzzz[i] * pa_y[i];

        tr_x_yyyyyy_zzzzzz[i] = 5.0 * tr_x_yyyy_zzzzzz[i] * fe_0 + tr_x_yyyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 616-644 components of targeted buffer : II

    auto tr_x_yyyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 616);

    auto tr_x_yyyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 617);

    auto tr_x_yyyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 618);

    auto tr_x_yyyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 619);

    auto tr_x_yyyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 620);

    auto tr_x_yyyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 621);

    auto tr_x_yyyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 622);

    auto tr_x_yyyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 623);

    auto tr_x_yyyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 624);

    auto tr_x_yyyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 625);

    auto tr_x_yyyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 626);

    auto tr_x_yyyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 627);

    auto tr_x_yyyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 628);

    auto tr_x_yyyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 629);

    auto tr_x_yyyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 630);

    auto tr_x_yyyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 631);

    auto tr_x_yyyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 632);

    auto tr_x_yyyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 633);

    auto tr_x_yyyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 634);

    auto tr_x_yyyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 635);

    auto tr_x_yyyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 636);

    auto tr_x_yyyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 637);

    auto tr_x_yyyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 638);

    auto tr_x_yyyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 639);

    auto tr_x_yyyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 640);

    auto tr_x_yyyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 641);

    auto tr_x_yyyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 642);

    auto tr_x_yyyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 643);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyyy_xxxxxx, tr_x_yyyyy_xxxxxy, tr_x_yyyyy_xxxxy, tr_x_yyyyy_xxxxyy, tr_x_yyyyy_xxxxyz, tr_x_yyyyy_xxxyy, tr_x_yyyyy_xxxyyy, tr_x_yyyyy_xxxyyz, tr_x_yyyyy_xxxyz, tr_x_yyyyy_xxxyzz, tr_x_yyyyy_xxyyy, tr_x_yyyyy_xxyyyy, tr_x_yyyyy_xxyyyz, tr_x_yyyyy_xxyyz, tr_x_yyyyy_xxyyzz, tr_x_yyyyy_xxyzz, tr_x_yyyyy_xxyzzz, tr_x_yyyyy_xyyyy, tr_x_yyyyy_xyyyyy, tr_x_yyyyy_xyyyyz, tr_x_yyyyy_xyyyz, tr_x_yyyyy_xyyyzz, tr_x_yyyyy_xyyzz, tr_x_yyyyy_xyyzzz, tr_x_yyyyy_xyzzz, tr_x_yyyyy_xyzzzz, tr_x_yyyyy_yyyyy, tr_x_yyyyy_yyyyyy, tr_x_yyyyy_yyyyyz, tr_x_yyyyy_yyyyz, tr_x_yyyyy_yyyyzz, tr_x_yyyyy_yyyzz, tr_x_yyyyy_yyyzzz, tr_x_yyyyy_yyzzz, tr_x_yyyyy_yyzzzz, tr_x_yyyyy_yzzzz, tr_x_yyyyy_yzzzzz, tr_x_yyyyyz_xxxxxx, tr_x_yyyyyz_xxxxxy, tr_x_yyyyyz_xxxxxz, tr_x_yyyyyz_xxxxyy, tr_x_yyyyyz_xxxxyz, tr_x_yyyyyz_xxxxzz, tr_x_yyyyyz_xxxyyy, tr_x_yyyyyz_xxxyyz, tr_x_yyyyyz_xxxyzz, tr_x_yyyyyz_xxxzzz, tr_x_yyyyyz_xxyyyy, tr_x_yyyyyz_xxyyyz, tr_x_yyyyyz_xxyyzz, tr_x_yyyyyz_xxyzzz, tr_x_yyyyyz_xxzzzz, tr_x_yyyyyz_xyyyyy, tr_x_yyyyyz_xyyyyz, tr_x_yyyyyz_xyyyzz, tr_x_yyyyyz_xyyzzz, tr_x_yyyyyz_xyzzzz, tr_x_yyyyyz_xzzzzz, tr_x_yyyyyz_yyyyyy, tr_x_yyyyyz_yyyyyz, tr_x_yyyyyz_yyyyzz, tr_x_yyyyyz_yyyzzz, tr_x_yyyyyz_yyzzzz, tr_x_yyyyyz_yzzzzz, tr_x_yyyyyz_zzzzzz, tr_x_yyyyz_xxxxxz, tr_x_yyyyz_xxxxzz, tr_x_yyyyz_xxxzzz, tr_x_yyyyz_xxzzzz, tr_x_yyyyz_xzzzzz, tr_x_yyyyz_zzzzzz, tr_x_yyyz_xxxxxz, tr_x_yyyz_xxxxzz, tr_x_yyyz_xxxzzz, tr_x_yyyz_xxzzzz, tr_x_yyyz_xzzzzz, tr_x_yyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyz_xxxxxx[i] = tr_x_yyyyy_xxxxxx[i] * pa_z[i];

        tr_x_yyyyyz_xxxxxy[i] = tr_x_yyyyy_xxxxxy[i] * pa_z[i];

        tr_x_yyyyyz_xxxxxz[i] = 4.0 * tr_x_yyyz_xxxxxz[i] * fe_0 + tr_x_yyyyz_xxxxxz[i] * pa_y[i];

        tr_x_yyyyyz_xxxxyy[i] = tr_x_yyyyy_xxxxyy[i] * pa_z[i];

        tr_x_yyyyyz_xxxxyz[i] = tr_x_yyyyy_xxxxy[i] * fe_0 + tr_x_yyyyy_xxxxyz[i] * pa_z[i];

        tr_x_yyyyyz_xxxxzz[i] = 4.0 * tr_x_yyyz_xxxxzz[i] * fe_0 + tr_x_yyyyz_xxxxzz[i] * pa_y[i];

        tr_x_yyyyyz_xxxyyy[i] = tr_x_yyyyy_xxxyyy[i] * pa_z[i];

        tr_x_yyyyyz_xxxyyz[i] = tr_x_yyyyy_xxxyy[i] * fe_0 + tr_x_yyyyy_xxxyyz[i] * pa_z[i];

        tr_x_yyyyyz_xxxyzz[i] = 2.0 * tr_x_yyyyy_xxxyz[i] * fe_0 + tr_x_yyyyy_xxxyzz[i] * pa_z[i];

        tr_x_yyyyyz_xxxzzz[i] = 4.0 * tr_x_yyyz_xxxzzz[i] * fe_0 + tr_x_yyyyz_xxxzzz[i] * pa_y[i];

        tr_x_yyyyyz_xxyyyy[i] = tr_x_yyyyy_xxyyyy[i] * pa_z[i];

        tr_x_yyyyyz_xxyyyz[i] = tr_x_yyyyy_xxyyy[i] * fe_0 + tr_x_yyyyy_xxyyyz[i] * pa_z[i];

        tr_x_yyyyyz_xxyyzz[i] = 2.0 * tr_x_yyyyy_xxyyz[i] * fe_0 + tr_x_yyyyy_xxyyzz[i] * pa_z[i];

        tr_x_yyyyyz_xxyzzz[i] = 3.0 * tr_x_yyyyy_xxyzz[i] * fe_0 + tr_x_yyyyy_xxyzzz[i] * pa_z[i];

        tr_x_yyyyyz_xxzzzz[i] = 4.0 * tr_x_yyyz_xxzzzz[i] * fe_0 + tr_x_yyyyz_xxzzzz[i] * pa_y[i];

        tr_x_yyyyyz_xyyyyy[i] = tr_x_yyyyy_xyyyyy[i] * pa_z[i];

        tr_x_yyyyyz_xyyyyz[i] = tr_x_yyyyy_xyyyy[i] * fe_0 + tr_x_yyyyy_xyyyyz[i] * pa_z[i];

        tr_x_yyyyyz_xyyyzz[i] = 2.0 * tr_x_yyyyy_xyyyz[i] * fe_0 + tr_x_yyyyy_xyyyzz[i] * pa_z[i];

        tr_x_yyyyyz_xyyzzz[i] = 3.0 * tr_x_yyyyy_xyyzz[i] * fe_0 + tr_x_yyyyy_xyyzzz[i] * pa_z[i];

        tr_x_yyyyyz_xyzzzz[i] = 4.0 * tr_x_yyyyy_xyzzz[i] * fe_0 + tr_x_yyyyy_xyzzzz[i] * pa_z[i];

        tr_x_yyyyyz_xzzzzz[i] = 4.0 * tr_x_yyyz_xzzzzz[i] * fe_0 + tr_x_yyyyz_xzzzzz[i] * pa_y[i];

        tr_x_yyyyyz_yyyyyy[i] = tr_x_yyyyy_yyyyyy[i] * pa_z[i];

        tr_x_yyyyyz_yyyyyz[i] = tr_x_yyyyy_yyyyy[i] * fe_0 + tr_x_yyyyy_yyyyyz[i] * pa_z[i];

        tr_x_yyyyyz_yyyyzz[i] = 2.0 * tr_x_yyyyy_yyyyz[i] * fe_0 + tr_x_yyyyy_yyyyzz[i] * pa_z[i];

        tr_x_yyyyyz_yyyzzz[i] = 3.0 * tr_x_yyyyy_yyyzz[i] * fe_0 + tr_x_yyyyy_yyyzzz[i] * pa_z[i];

        tr_x_yyyyyz_yyzzzz[i] = 4.0 * tr_x_yyyyy_yyzzz[i] * fe_0 + tr_x_yyyyy_yyzzzz[i] * pa_z[i];

        tr_x_yyyyyz_yzzzzz[i] = 5.0 * tr_x_yyyyy_yzzzz[i] * fe_0 + tr_x_yyyyy_yzzzzz[i] * pa_z[i];

        tr_x_yyyyyz_zzzzzz[i] = 4.0 * tr_x_yyyz_zzzzzz[i] * fe_0 + tr_x_yyyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 644-672 components of targeted buffer : II

    auto tr_x_yyyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 644);

    auto tr_x_yyyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 645);

    auto tr_x_yyyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 646);

    auto tr_x_yyyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 647);

    auto tr_x_yyyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 648);

    auto tr_x_yyyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 649);

    auto tr_x_yyyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 650);

    auto tr_x_yyyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 651);

    auto tr_x_yyyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 652);

    auto tr_x_yyyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 653);

    auto tr_x_yyyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 654);

    auto tr_x_yyyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 655);

    auto tr_x_yyyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 656);

    auto tr_x_yyyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 657);

    auto tr_x_yyyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 658);

    auto tr_x_yyyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 659);

    auto tr_x_yyyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 660);

    auto tr_x_yyyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 661);

    auto tr_x_yyyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 662);

    auto tr_x_yyyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 663);

    auto tr_x_yyyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 664);

    auto tr_x_yyyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 665);

    auto tr_x_yyyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 666);

    auto tr_x_yyyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 667);

    auto tr_x_yyyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 668);

    auto tr_x_yyyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 669);

    auto tr_x_yyyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 670);

    auto tr_x_yyyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 671);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyy_xxxxxy, tr_x_yyyy_xxxxyy, tr_x_yyyy_xxxyyy, tr_x_yyyy_xxyyyy, tr_x_yyyy_xyyyyy, tr_x_yyyy_yyyyyy, tr_x_yyyyz_xxxxxy, tr_x_yyyyz_xxxxyy, tr_x_yyyyz_xxxyyy, tr_x_yyyyz_xxyyyy, tr_x_yyyyz_xyyyyy, tr_x_yyyyz_yyyyyy, tr_x_yyyyzz_xxxxxx, tr_x_yyyyzz_xxxxxy, tr_x_yyyyzz_xxxxxz, tr_x_yyyyzz_xxxxyy, tr_x_yyyyzz_xxxxyz, tr_x_yyyyzz_xxxxzz, tr_x_yyyyzz_xxxyyy, tr_x_yyyyzz_xxxyyz, tr_x_yyyyzz_xxxyzz, tr_x_yyyyzz_xxxzzz, tr_x_yyyyzz_xxyyyy, tr_x_yyyyzz_xxyyyz, tr_x_yyyyzz_xxyyzz, tr_x_yyyyzz_xxyzzz, tr_x_yyyyzz_xxzzzz, tr_x_yyyyzz_xyyyyy, tr_x_yyyyzz_xyyyyz, tr_x_yyyyzz_xyyyzz, tr_x_yyyyzz_xyyzzz, tr_x_yyyyzz_xyzzzz, tr_x_yyyyzz_xzzzzz, tr_x_yyyyzz_yyyyyy, tr_x_yyyyzz_yyyyyz, tr_x_yyyyzz_yyyyzz, tr_x_yyyyzz_yyyzzz, tr_x_yyyyzz_yyzzzz, tr_x_yyyyzz_yzzzzz, tr_x_yyyyzz_zzzzzz, tr_x_yyyzz_xxxxxx, tr_x_yyyzz_xxxxxz, tr_x_yyyzz_xxxxyz, tr_x_yyyzz_xxxxz, tr_x_yyyzz_xxxxzz, tr_x_yyyzz_xxxyyz, tr_x_yyyzz_xxxyz, tr_x_yyyzz_xxxyzz, tr_x_yyyzz_xxxzz, tr_x_yyyzz_xxxzzz, tr_x_yyyzz_xxyyyz, tr_x_yyyzz_xxyyz, tr_x_yyyzz_xxyyzz, tr_x_yyyzz_xxyzz, tr_x_yyyzz_xxyzzz, tr_x_yyyzz_xxzzz, tr_x_yyyzz_xxzzzz, tr_x_yyyzz_xyyyyz, tr_x_yyyzz_xyyyz, tr_x_yyyzz_xyyyzz, tr_x_yyyzz_xyyzz, tr_x_yyyzz_xyyzzz, tr_x_yyyzz_xyzzz, tr_x_yyyzz_xyzzzz, tr_x_yyyzz_xzzzz, tr_x_yyyzz_xzzzzz, tr_x_yyyzz_yyyyyz, tr_x_yyyzz_yyyyz, tr_x_yyyzz_yyyyzz, tr_x_yyyzz_yyyzz, tr_x_yyyzz_yyyzzz, tr_x_yyyzz_yyzzz, tr_x_yyyzz_yyzzzz, tr_x_yyyzz_yzzzz, tr_x_yyyzz_yzzzzz, tr_x_yyyzz_zzzzz, tr_x_yyyzz_zzzzzz, tr_x_yyzz_xxxxxx, tr_x_yyzz_xxxxxz, tr_x_yyzz_xxxxyz, tr_x_yyzz_xxxxzz, tr_x_yyzz_xxxyyz, tr_x_yyzz_xxxyzz, tr_x_yyzz_xxxzzz, tr_x_yyzz_xxyyyz, tr_x_yyzz_xxyyzz, tr_x_yyzz_xxyzzz, tr_x_yyzz_xxzzzz, tr_x_yyzz_xyyyyz, tr_x_yyzz_xyyyzz, tr_x_yyzz_xyyzzz, tr_x_yyzz_xyzzzz, tr_x_yyzz_xzzzzz, tr_x_yyzz_yyyyyz, tr_x_yyzz_yyyyzz, tr_x_yyzz_yyyzzz, tr_x_yyzz_yyzzzz, tr_x_yyzz_yzzzzz, tr_x_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyzz_xxxxxx[i] = 3.0 * tr_x_yyzz_xxxxxx[i] * fe_0 + tr_x_yyyzz_xxxxxx[i] * pa_y[i];

        tr_x_yyyyzz_xxxxxy[i] = tr_x_yyyy_xxxxxy[i] * fe_0 + tr_x_yyyyz_xxxxxy[i] * pa_z[i];

        tr_x_yyyyzz_xxxxxz[i] = 3.0 * tr_x_yyzz_xxxxxz[i] * fe_0 + tr_x_yyyzz_xxxxxz[i] * pa_y[i];

        tr_x_yyyyzz_xxxxyy[i] = tr_x_yyyy_xxxxyy[i] * fe_0 + tr_x_yyyyz_xxxxyy[i] * pa_z[i];

        tr_x_yyyyzz_xxxxyz[i] = 3.0 * tr_x_yyzz_xxxxyz[i] * fe_0 + tr_x_yyyzz_xxxxz[i] * fe_0 + tr_x_yyyzz_xxxxyz[i] * pa_y[i];

        tr_x_yyyyzz_xxxxzz[i] = 3.0 * tr_x_yyzz_xxxxzz[i] * fe_0 + tr_x_yyyzz_xxxxzz[i] * pa_y[i];

        tr_x_yyyyzz_xxxyyy[i] = tr_x_yyyy_xxxyyy[i] * fe_0 + tr_x_yyyyz_xxxyyy[i] * pa_z[i];

        tr_x_yyyyzz_xxxyyz[i] = 3.0 * tr_x_yyzz_xxxyyz[i] * fe_0 + 2.0 * tr_x_yyyzz_xxxyz[i] * fe_0 + tr_x_yyyzz_xxxyyz[i] * pa_y[i];

        tr_x_yyyyzz_xxxyzz[i] = 3.0 * tr_x_yyzz_xxxyzz[i] * fe_0 + tr_x_yyyzz_xxxzz[i] * fe_0 + tr_x_yyyzz_xxxyzz[i] * pa_y[i];

        tr_x_yyyyzz_xxxzzz[i] = 3.0 * tr_x_yyzz_xxxzzz[i] * fe_0 + tr_x_yyyzz_xxxzzz[i] * pa_y[i];

        tr_x_yyyyzz_xxyyyy[i] = tr_x_yyyy_xxyyyy[i] * fe_0 + tr_x_yyyyz_xxyyyy[i] * pa_z[i];

        tr_x_yyyyzz_xxyyyz[i] = 3.0 * tr_x_yyzz_xxyyyz[i] * fe_0 + 3.0 * tr_x_yyyzz_xxyyz[i] * fe_0 + tr_x_yyyzz_xxyyyz[i] * pa_y[i];

        tr_x_yyyyzz_xxyyzz[i] = 3.0 * tr_x_yyzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_yyyzz_xxyzz[i] * fe_0 + tr_x_yyyzz_xxyyzz[i] * pa_y[i];

        tr_x_yyyyzz_xxyzzz[i] = 3.0 * tr_x_yyzz_xxyzzz[i] * fe_0 + tr_x_yyyzz_xxzzz[i] * fe_0 + tr_x_yyyzz_xxyzzz[i] * pa_y[i];

        tr_x_yyyyzz_xxzzzz[i] = 3.0 * tr_x_yyzz_xxzzzz[i] * fe_0 + tr_x_yyyzz_xxzzzz[i] * pa_y[i];

        tr_x_yyyyzz_xyyyyy[i] = tr_x_yyyy_xyyyyy[i] * fe_0 + tr_x_yyyyz_xyyyyy[i] * pa_z[i];

        tr_x_yyyyzz_xyyyyz[i] = 3.0 * tr_x_yyzz_xyyyyz[i] * fe_0 + 4.0 * tr_x_yyyzz_xyyyz[i] * fe_0 + tr_x_yyyzz_xyyyyz[i] * pa_y[i];

        tr_x_yyyyzz_xyyyzz[i] = 3.0 * tr_x_yyzz_xyyyzz[i] * fe_0 + 3.0 * tr_x_yyyzz_xyyzz[i] * fe_0 + tr_x_yyyzz_xyyyzz[i] * pa_y[i];

        tr_x_yyyyzz_xyyzzz[i] = 3.0 * tr_x_yyzz_xyyzzz[i] * fe_0 + 2.0 * tr_x_yyyzz_xyzzz[i] * fe_0 + tr_x_yyyzz_xyyzzz[i] * pa_y[i];

        tr_x_yyyyzz_xyzzzz[i] = 3.0 * tr_x_yyzz_xyzzzz[i] * fe_0 + tr_x_yyyzz_xzzzz[i] * fe_0 + tr_x_yyyzz_xyzzzz[i] * pa_y[i];

        tr_x_yyyyzz_xzzzzz[i] = 3.0 * tr_x_yyzz_xzzzzz[i] * fe_0 + tr_x_yyyzz_xzzzzz[i] * pa_y[i];

        tr_x_yyyyzz_yyyyyy[i] = tr_x_yyyy_yyyyyy[i] * fe_0 + tr_x_yyyyz_yyyyyy[i] * pa_z[i];

        tr_x_yyyyzz_yyyyyz[i] = 3.0 * tr_x_yyzz_yyyyyz[i] * fe_0 + 5.0 * tr_x_yyyzz_yyyyz[i] * fe_0 + tr_x_yyyzz_yyyyyz[i] * pa_y[i];

        tr_x_yyyyzz_yyyyzz[i] = 3.0 * tr_x_yyzz_yyyyzz[i] * fe_0 + 4.0 * tr_x_yyyzz_yyyzz[i] * fe_0 + tr_x_yyyzz_yyyyzz[i] * pa_y[i];

        tr_x_yyyyzz_yyyzzz[i] = 3.0 * tr_x_yyzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_yyyzz_yyzzz[i] * fe_0 + tr_x_yyyzz_yyyzzz[i] * pa_y[i];

        tr_x_yyyyzz_yyzzzz[i] = 3.0 * tr_x_yyzz_yyzzzz[i] * fe_0 + 2.0 * tr_x_yyyzz_yzzzz[i] * fe_0 + tr_x_yyyzz_yyzzzz[i] * pa_y[i];

        tr_x_yyyyzz_yzzzzz[i] = 3.0 * tr_x_yyzz_yzzzzz[i] * fe_0 + tr_x_yyyzz_zzzzz[i] * fe_0 + tr_x_yyyzz_yzzzzz[i] * pa_y[i];

        tr_x_yyyyzz_zzzzzz[i] = 3.0 * tr_x_yyzz_zzzzzz[i] * fe_0 + tr_x_yyyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 672-700 components of targeted buffer : II

    auto tr_x_yyyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 672);

    auto tr_x_yyyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 673);

    auto tr_x_yyyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 674);

    auto tr_x_yyyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 675);

    auto tr_x_yyyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 676);

    auto tr_x_yyyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 677);

    auto tr_x_yyyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 678);

    auto tr_x_yyyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 679);

    auto tr_x_yyyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 680);

    auto tr_x_yyyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 681);

    auto tr_x_yyyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 682);

    auto tr_x_yyyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 683);

    auto tr_x_yyyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 684);

    auto tr_x_yyyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 685);

    auto tr_x_yyyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 686);

    auto tr_x_yyyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 687);

    auto tr_x_yyyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 688);

    auto tr_x_yyyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 689);

    auto tr_x_yyyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 690);

    auto tr_x_yyyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 691);

    auto tr_x_yyyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 692);

    auto tr_x_yyyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 693);

    auto tr_x_yyyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 694);

    auto tr_x_yyyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 695);

    auto tr_x_yyyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 696);

    auto tr_x_yyyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 697);

    auto tr_x_yyyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 698);

    auto tr_x_yyyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 699);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyz_xxxxxy, tr_x_yyyz_xxxxyy, tr_x_yyyz_xxxyyy, tr_x_yyyz_xxyyyy, tr_x_yyyz_xyyyyy, tr_x_yyyz_yyyyyy, tr_x_yyyzz_xxxxxy, tr_x_yyyzz_xxxxyy, tr_x_yyyzz_xxxyyy, tr_x_yyyzz_xxyyyy, tr_x_yyyzz_xyyyyy, tr_x_yyyzz_yyyyyy, tr_x_yyyzzz_xxxxxx, tr_x_yyyzzz_xxxxxy, tr_x_yyyzzz_xxxxxz, tr_x_yyyzzz_xxxxyy, tr_x_yyyzzz_xxxxyz, tr_x_yyyzzz_xxxxzz, tr_x_yyyzzz_xxxyyy, tr_x_yyyzzz_xxxyyz, tr_x_yyyzzz_xxxyzz, tr_x_yyyzzz_xxxzzz, tr_x_yyyzzz_xxyyyy, tr_x_yyyzzz_xxyyyz, tr_x_yyyzzz_xxyyzz, tr_x_yyyzzz_xxyzzz, tr_x_yyyzzz_xxzzzz, tr_x_yyyzzz_xyyyyy, tr_x_yyyzzz_xyyyyz, tr_x_yyyzzz_xyyyzz, tr_x_yyyzzz_xyyzzz, tr_x_yyyzzz_xyzzzz, tr_x_yyyzzz_xzzzzz, tr_x_yyyzzz_yyyyyy, tr_x_yyyzzz_yyyyyz, tr_x_yyyzzz_yyyyzz, tr_x_yyyzzz_yyyzzz, tr_x_yyyzzz_yyzzzz, tr_x_yyyzzz_yzzzzz, tr_x_yyyzzz_zzzzzz, tr_x_yyzzz_xxxxxx, tr_x_yyzzz_xxxxxz, tr_x_yyzzz_xxxxyz, tr_x_yyzzz_xxxxz, tr_x_yyzzz_xxxxzz, tr_x_yyzzz_xxxyyz, tr_x_yyzzz_xxxyz, tr_x_yyzzz_xxxyzz, tr_x_yyzzz_xxxzz, tr_x_yyzzz_xxxzzz, tr_x_yyzzz_xxyyyz, tr_x_yyzzz_xxyyz, tr_x_yyzzz_xxyyzz, tr_x_yyzzz_xxyzz, tr_x_yyzzz_xxyzzz, tr_x_yyzzz_xxzzz, tr_x_yyzzz_xxzzzz, tr_x_yyzzz_xyyyyz, tr_x_yyzzz_xyyyz, tr_x_yyzzz_xyyyzz, tr_x_yyzzz_xyyzz, tr_x_yyzzz_xyyzzz, tr_x_yyzzz_xyzzz, tr_x_yyzzz_xyzzzz, tr_x_yyzzz_xzzzz, tr_x_yyzzz_xzzzzz, tr_x_yyzzz_yyyyyz, tr_x_yyzzz_yyyyz, tr_x_yyzzz_yyyyzz, tr_x_yyzzz_yyyzz, tr_x_yyzzz_yyyzzz, tr_x_yyzzz_yyzzz, tr_x_yyzzz_yyzzzz, tr_x_yyzzz_yzzzz, tr_x_yyzzz_yzzzzz, tr_x_yyzzz_zzzzz, tr_x_yyzzz_zzzzzz, tr_x_yzzz_xxxxxx, tr_x_yzzz_xxxxxz, tr_x_yzzz_xxxxyz, tr_x_yzzz_xxxxzz, tr_x_yzzz_xxxyyz, tr_x_yzzz_xxxyzz, tr_x_yzzz_xxxzzz, tr_x_yzzz_xxyyyz, tr_x_yzzz_xxyyzz, tr_x_yzzz_xxyzzz, tr_x_yzzz_xxzzzz, tr_x_yzzz_xyyyyz, tr_x_yzzz_xyyyzz, tr_x_yzzz_xyyzzz, tr_x_yzzz_xyzzzz, tr_x_yzzz_xzzzzz, tr_x_yzzz_yyyyyz, tr_x_yzzz_yyyyzz, tr_x_yzzz_yyyzzz, tr_x_yzzz_yyzzzz, tr_x_yzzz_yzzzzz, tr_x_yzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzzz_xxxxxx[i] = 2.0 * tr_x_yzzz_xxxxxx[i] * fe_0 + tr_x_yyzzz_xxxxxx[i] * pa_y[i];

        tr_x_yyyzzz_xxxxxy[i] = 2.0 * tr_x_yyyz_xxxxxy[i] * fe_0 + tr_x_yyyzz_xxxxxy[i] * pa_z[i];

        tr_x_yyyzzz_xxxxxz[i] = 2.0 * tr_x_yzzz_xxxxxz[i] * fe_0 + tr_x_yyzzz_xxxxxz[i] * pa_y[i];

        tr_x_yyyzzz_xxxxyy[i] = 2.0 * tr_x_yyyz_xxxxyy[i] * fe_0 + tr_x_yyyzz_xxxxyy[i] * pa_z[i];

        tr_x_yyyzzz_xxxxyz[i] = 2.0 * tr_x_yzzz_xxxxyz[i] * fe_0 + tr_x_yyzzz_xxxxz[i] * fe_0 + tr_x_yyzzz_xxxxyz[i] * pa_y[i];

        tr_x_yyyzzz_xxxxzz[i] = 2.0 * tr_x_yzzz_xxxxzz[i] * fe_0 + tr_x_yyzzz_xxxxzz[i] * pa_y[i];

        tr_x_yyyzzz_xxxyyy[i] = 2.0 * tr_x_yyyz_xxxyyy[i] * fe_0 + tr_x_yyyzz_xxxyyy[i] * pa_z[i];

        tr_x_yyyzzz_xxxyyz[i] = 2.0 * tr_x_yzzz_xxxyyz[i] * fe_0 + 2.0 * tr_x_yyzzz_xxxyz[i] * fe_0 + tr_x_yyzzz_xxxyyz[i] * pa_y[i];

        tr_x_yyyzzz_xxxyzz[i] = 2.0 * tr_x_yzzz_xxxyzz[i] * fe_0 + tr_x_yyzzz_xxxzz[i] * fe_0 + tr_x_yyzzz_xxxyzz[i] * pa_y[i];

        tr_x_yyyzzz_xxxzzz[i] = 2.0 * tr_x_yzzz_xxxzzz[i] * fe_0 + tr_x_yyzzz_xxxzzz[i] * pa_y[i];

        tr_x_yyyzzz_xxyyyy[i] = 2.0 * tr_x_yyyz_xxyyyy[i] * fe_0 + tr_x_yyyzz_xxyyyy[i] * pa_z[i];

        tr_x_yyyzzz_xxyyyz[i] = 2.0 * tr_x_yzzz_xxyyyz[i] * fe_0 + 3.0 * tr_x_yyzzz_xxyyz[i] * fe_0 + tr_x_yyzzz_xxyyyz[i] * pa_y[i];

        tr_x_yyyzzz_xxyyzz[i] = 2.0 * tr_x_yzzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_yyzzz_xxyzz[i] * fe_0 + tr_x_yyzzz_xxyyzz[i] * pa_y[i];

        tr_x_yyyzzz_xxyzzz[i] = 2.0 * tr_x_yzzz_xxyzzz[i] * fe_0 + tr_x_yyzzz_xxzzz[i] * fe_0 + tr_x_yyzzz_xxyzzz[i] * pa_y[i];

        tr_x_yyyzzz_xxzzzz[i] = 2.0 * tr_x_yzzz_xxzzzz[i] * fe_0 + tr_x_yyzzz_xxzzzz[i] * pa_y[i];

        tr_x_yyyzzz_xyyyyy[i] = 2.0 * tr_x_yyyz_xyyyyy[i] * fe_0 + tr_x_yyyzz_xyyyyy[i] * pa_z[i];

        tr_x_yyyzzz_xyyyyz[i] = 2.0 * tr_x_yzzz_xyyyyz[i] * fe_0 + 4.0 * tr_x_yyzzz_xyyyz[i] * fe_0 + tr_x_yyzzz_xyyyyz[i] * pa_y[i];

        tr_x_yyyzzz_xyyyzz[i] = 2.0 * tr_x_yzzz_xyyyzz[i] * fe_0 + 3.0 * tr_x_yyzzz_xyyzz[i] * fe_0 + tr_x_yyzzz_xyyyzz[i] * pa_y[i];

        tr_x_yyyzzz_xyyzzz[i] = 2.0 * tr_x_yzzz_xyyzzz[i] * fe_0 + 2.0 * tr_x_yyzzz_xyzzz[i] * fe_0 + tr_x_yyzzz_xyyzzz[i] * pa_y[i];

        tr_x_yyyzzz_xyzzzz[i] = 2.0 * tr_x_yzzz_xyzzzz[i] * fe_0 + tr_x_yyzzz_xzzzz[i] * fe_0 + tr_x_yyzzz_xyzzzz[i] * pa_y[i];

        tr_x_yyyzzz_xzzzzz[i] = 2.0 * tr_x_yzzz_xzzzzz[i] * fe_0 + tr_x_yyzzz_xzzzzz[i] * pa_y[i];

        tr_x_yyyzzz_yyyyyy[i] = 2.0 * tr_x_yyyz_yyyyyy[i] * fe_0 + tr_x_yyyzz_yyyyyy[i] * pa_z[i];

        tr_x_yyyzzz_yyyyyz[i] = 2.0 * tr_x_yzzz_yyyyyz[i] * fe_0 + 5.0 * tr_x_yyzzz_yyyyz[i] * fe_0 + tr_x_yyzzz_yyyyyz[i] * pa_y[i];

        tr_x_yyyzzz_yyyyzz[i] = 2.0 * tr_x_yzzz_yyyyzz[i] * fe_0 + 4.0 * tr_x_yyzzz_yyyzz[i] * fe_0 + tr_x_yyzzz_yyyyzz[i] * pa_y[i];

        tr_x_yyyzzz_yyyzzz[i] = 2.0 * tr_x_yzzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_yyzzz_yyzzz[i] * fe_0 + tr_x_yyzzz_yyyzzz[i] * pa_y[i];

        tr_x_yyyzzz_yyzzzz[i] = 2.0 * tr_x_yzzz_yyzzzz[i] * fe_0 + 2.0 * tr_x_yyzzz_yzzzz[i] * fe_0 + tr_x_yyzzz_yyzzzz[i] * pa_y[i];

        tr_x_yyyzzz_yzzzzz[i] = 2.0 * tr_x_yzzz_yzzzzz[i] * fe_0 + tr_x_yyzzz_zzzzz[i] * fe_0 + tr_x_yyzzz_yzzzzz[i] * pa_y[i];

        tr_x_yyyzzz_zzzzzz[i] = 2.0 * tr_x_yzzz_zzzzzz[i] * fe_0 + tr_x_yyzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 700-728 components of targeted buffer : II

    auto tr_x_yyzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 700);

    auto tr_x_yyzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 701);

    auto tr_x_yyzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 702);

    auto tr_x_yyzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 703);

    auto tr_x_yyzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 704);

    auto tr_x_yyzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 705);

    auto tr_x_yyzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 706);

    auto tr_x_yyzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 707);

    auto tr_x_yyzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 708);

    auto tr_x_yyzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 709);

    auto tr_x_yyzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 710);

    auto tr_x_yyzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 711);

    auto tr_x_yyzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 712);

    auto tr_x_yyzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 713);

    auto tr_x_yyzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 714);

    auto tr_x_yyzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 715);

    auto tr_x_yyzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 716);

    auto tr_x_yyzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 717);

    auto tr_x_yyzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 718);

    auto tr_x_yyzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 719);

    auto tr_x_yyzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 720);

    auto tr_x_yyzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 721);

    auto tr_x_yyzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 722);

    auto tr_x_yyzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 723);

    auto tr_x_yyzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 724);

    auto tr_x_yyzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 725);

    auto tr_x_yyzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 726);

    auto tr_x_yyzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 727);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyzz_xxxxxy, tr_x_yyzz_xxxxyy, tr_x_yyzz_xxxyyy, tr_x_yyzz_xxyyyy, tr_x_yyzz_xyyyyy, tr_x_yyzz_yyyyyy, tr_x_yyzzz_xxxxxy, tr_x_yyzzz_xxxxyy, tr_x_yyzzz_xxxyyy, tr_x_yyzzz_xxyyyy, tr_x_yyzzz_xyyyyy, tr_x_yyzzz_yyyyyy, tr_x_yyzzzz_xxxxxx, tr_x_yyzzzz_xxxxxy, tr_x_yyzzzz_xxxxxz, tr_x_yyzzzz_xxxxyy, tr_x_yyzzzz_xxxxyz, tr_x_yyzzzz_xxxxzz, tr_x_yyzzzz_xxxyyy, tr_x_yyzzzz_xxxyyz, tr_x_yyzzzz_xxxyzz, tr_x_yyzzzz_xxxzzz, tr_x_yyzzzz_xxyyyy, tr_x_yyzzzz_xxyyyz, tr_x_yyzzzz_xxyyzz, tr_x_yyzzzz_xxyzzz, tr_x_yyzzzz_xxzzzz, tr_x_yyzzzz_xyyyyy, tr_x_yyzzzz_xyyyyz, tr_x_yyzzzz_xyyyzz, tr_x_yyzzzz_xyyzzz, tr_x_yyzzzz_xyzzzz, tr_x_yyzzzz_xzzzzz, tr_x_yyzzzz_yyyyyy, tr_x_yyzzzz_yyyyyz, tr_x_yyzzzz_yyyyzz, tr_x_yyzzzz_yyyzzz, tr_x_yyzzzz_yyzzzz, tr_x_yyzzzz_yzzzzz, tr_x_yyzzzz_zzzzzz, tr_x_yzzzz_xxxxxx, tr_x_yzzzz_xxxxxz, tr_x_yzzzz_xxxxyz, tr_x_yzzzz_xxxxz, tr_x_yzzzz_xxxxzz, tr_x_yzzzz_xxxyyz, tr_x_yzzzz_xxxyz, tr_x_yzzzz_xxxyzz, tr_x_yzzzz_xxxzz, tr_x_yzzzz_xxxzzz, tr_x_yzzzz_xxyyyz, tr_x_yzzzz_xxyyz, tr_x_yzzzz_xxyyzz, tr_x_yzzzz_xxyzz, tr_x_yzzzz_xxyzzz, tr_x_yzzzz_xxzzz, tr_x_yzzzz_xxzzzz, tr_x_yzzzz_xyyyyz, tr_x_yzzzz_xyyyz, tr_x_yzzzz_xyyyzz, tr_x_yzzzz_xyyzz, tr_x_yzzzz_xyyzzz, tr_x_yzzzz_xyzzz, tr_x_yzzzz_xyzzzz, tr_x_yzzzz_xzzzz, tr_x_yzzzz_xzzzzz, tr_x_yzzzz_yyyyyz, tr_x_yzzzz_yyyyz, tr_x_yzzzz_yyyyzz, tr_x_yzzzz_yyyzz, tr_x_yzzzz_yyyzzz, tr_x_yzzzz_yyzzz, tr_x_yzzzz_yyzzzz, tr_x_yzzzz_yzzzz, tr_x_yzzzz_yzzzzz, tr_x_yzzzz_zzzzz, tr_x_yzzzz_zzzzzz, tr_x_zzzz_xxxxxx, tr_x_zzzz_xxxxxz, tr_x_zzzz_xxxxyz, tr_x_zzzz_xxxxzz, tr_x_zzzz_xxxyyz, tr_x_zzzz_xxxyzz, tr_x_zzzz_xxxzzz, tr_x_zzzz_xxyyyz, tr_x_zzzz_xxyyzz, tr_x_zzzz_xxyzzz, tr_x_zzzz_xxzzzz, tr_x_zzzz_xyyyyz, tr_x_zzzz_xyyyzz, tr_x_zzzz_xyyzzz, tr_x_zzzz_xyzzzz, tr_x_zzzz_xzzzzz, tr_x_zzzz_yyyyyz, tr_x_zzzz_yyyyzz, tr_x_zzzz_yyyzzz, tr_x_zzzz_yyzzzz, tr_x_zzzz_yzzzzz, tr_x_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzzz_xxxxxx[i] = tr_x_zzzz_xxxxxx[i] * fe_0 + tr_x_yzzzz_xxxxxx[i] * pa_y[i];

        tr_x_yyzzzz_xxxxxy[i] = 3.0 * tr_x_yyzz_xxxxxy[i] * fe_0 + tr_x_yyzzz_xxxxxy[i] * pa_z[i];

        tr_x_yyzzzz_xxxxxz[i] = tr_x_zzzz_xxxxxz[i] * fe_0 + tr_x_yzzzz_xxxxxz[i] * pa_y[i];

        tr_x_yyzzzz_xxxxyy[i] = 3.0 * tr_x_yyzz_xxxxyy[i] * fe_0 + tr_x_yyzzz_xxxxyy[i] * pa_z[i];

        tr_x_yyzzzz_xxxxyz[i] = tr_x_zzzz_xxxxyz[i] * fe_0 + tr_x_yzzzz_xxxxz[i] * fe_0 + tr_x_yzzzz_xxxxyz[i] * pa_y[i];

        tr_x_yyzzzz_xxxxzz[i] = tr_x_zzzz_xxxxzz[i] * fe_0 + tr_x_yzzzz_xxxxzz[i] * pa_y[i];

        tr_x_yyzzzz_xxxyyy[i] = 3.0 * tr_x_yyzz_xxxyyy[i] * fe_0 + tr_x_yyzzz_xxxyyy[i] * pa_z[i];

        tr_x_yyzzzz_xxxyyz[i] = tr_x_zzzz_xxxyyz[i] * fe_0 + 2.0 * tr_x_yzzzz_xxxyz[i] * fe_0 + tr_x_yzzzz_xxxyyz[i] * pa_y[i];

        tr_x_yyzzzz_xxxyzz[i] = tr_x_zzzz_xxxyzz[i] * fe_0 + tr_x_yzzzz_xxxzz[i] * fe_0 + tr_x_yzzzz_xxxyzz[i] * pa_y[i];

        tr_x_yyzzzz_xxxzzz[i] = tr_x_zzzz_xxxzzz[i] * fe_0 + tr_x_yzzzz_xxxzzz[i] * pa_y[i];

        tr_x_yyzzzz_xxyyyy[i] = 3.0 * tr_x_yyzz_xxyyyy[i] * fe_0 + tr_x_yyzzz_xxyyyy[i] * pa_z[i];

        tr_x_yyzzzz_xxyyyz[i] = tr_x_zzzz_xxyyyz[i] * fe_0 + 3.0 * tr_x_yzzzz_xxyyz[i] * fe_0 + tr_x_yzzzz_xxyyyz[i] * pa_y[i];

        tr_x_yyzzzz_xxyyzz[i] = tr_x_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_yzzzz_xxyzz[i] * fe_0 + tr_x_yzzzz_xxyyzz[i] * pa_y[i];

        tr_x_yyzzzz_xxyzzz[i] = tr_x_zzzz_xxyzzz[i] * fe_0 + tr_x_yzzzz_xxzzz[i] * fe_0 + tr_x_yzzzz_xxyzzz[i] * pa_y[i];

        tr_x_yyzzzz_xxzzzz[i] = tr_x_zzzz_xxzzzz[i] * fe_0 + tr_x_yzzzz_xxzzzz[i] * pa_y[i];

        tr_x_yyzzzz_xyyyyy[i] = 3.0 * tr_x_yyzz_xyyyyy[i] * fe_0 + tr_x_yyzzz_xyyyyy[i] * pa_z[i];

        tr_x_yyzzzz_xyyyyz[i] = tr_x_zzzz_xyyyyz[i] * fe_0 + 4.0 * tr_x_yzzzz_xyyyz[i] * fe_0 + tr_x_yzzzz_xyyyyz[i] * pa_y[i];

        tr_x_yyzzzz_xyyyzz[i] = tr_x_zzzz_xyyyzz[i] * fe_0 + 3.0 * tr_x_yzzzz_xyyzz[i] * fe_0 + tr_x_yzzzz_xyyyzz[i] * pa_y[i];

        tr_x_yyzzzz_xyyzzz[i] = tr_x_zzzz_xyyzzz[i] * fe_0 + 2.0 * tr_x_yzzzz_xyzzz[i] * fe_0 + tr_x_yzzzz_xyyzzz[i] * pa_y[i];

        tr_x_yyzzzz_xyzzzz[i] = tr_x_zzzz_xyzzzz[i] * fe_0 + tr_x_yzzzz_xzzzz[i] * fe_0 + tr_x_yzzzz_xyzzzz[i] * pa_y[i];

        tr_x_yyzzzz_xzzzzz[i] = tr_x_zzzz_xzzzzz[i] * fe_0 + tr_x_yzzzz_xzzzzz[i] * pa_y[i];

        tr_x_yyzzzz_yyyyyy[i] = 3.0 * tr_x_yyzz_yyyyyy[i] * fe_0 + tr_x_yyzzz_yyyyyy[i] * pa_z[i];

        tr_x_yyzzzz_yyyyyz[i] = tr_x_zzzz_yyyyyz[i] * fe_0 + 5.0 * tr_x_yzzzz_yyyyz[i] * fe_0 + tr_x_yzzzz_yyyyyz[i] * pa_y[i];

        tr_x_yyzzzz_yyyyzz[i] = tr_x_zzzz_yyyyzz[i] * fe_0 + 4.0 * tr_x_yzzzz_yyyzz[i] * fe_0 + tr_x_yzzzz_yyyyzz[i] * pa_y[i];

        tr_x_yyzzzz_yyyzzz[i] = tr_x_zzzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_yzzzz_yyzzz[i] * fe_0 + tr_x_yzzzz_yyyzzz[i] * pa_y[i];

        tr_x_yyzzzz_yyzzzz[i] = tr_x_zzzz_yyzzzz[i] * fe_0 + 2.0 * tr_x_yzzzz_yzzzz[i] * fe_0 + tr_x_yzzzz_yyzzzz[i] * pa_y[i];

        tr_x_yyzzzz_yzzzzz[i] = tr_x_zzzz_yzzzzz[i] * fe_0 + tr_x_yzzzz_zzzzz[i] * fe_0 + tr_x_yzzzz_yzzzzz[i] * pa_y[i];

        tr_x_yyzzzz_zzzzzz[i] = tr_x_zzzz_zzzzzz[i] * fe_0 + tr_x_yzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 728-756 components of targeted buffer : II

    auto tr_x_yzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 728);

    auto tr_x_yzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 729);

    auto tr_x_yzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 730);

    auto tr_x_yzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 731);

    auto tr_x_yzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 732);

    auto tr_x_yzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 733);

    auto tr_x_yzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 734);

    auto tr_x_yzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 735);

    auto tr_x_yzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 736);

    auto tr_x_yzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 737);

    auto tr_x_yzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 738);

    auto tr_x_yzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 739);

    auto tr_x_yzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 740);

    auto tr_x_yzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 741);

    auto tr_x_yzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 742);

    auto tr_x_yzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 743);

    auto tr_x_yzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 744);

    auto tr_x_yzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 745);

    auto tr_x_yzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 746);

    auto tr_x_yzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 747);

    auto tr_x_yzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 748);

    auto tr_x_yzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 749);

    auto tr_x_yzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 750);

    auto tr_x_yzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 751);

    auto tr_x_yzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 752);

    auto tr_x_yzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 753);

    auto tr_x_yzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 754);

    auto tr_x_yzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 755);

    #pragma omp simd aligned(pa_y, tr_x_yzzzzz_xxxxxx, tr_x_yzzzzz_xxxxxy, tr_x_yzzzzz_xxxxxz, tr_x_yzzzzz_xxxxyy, tr_x_yzzzzz_xxxxyz, tr_x_yzzzzz_xxxxzz, tr_x_yzzzzz_xxxyyy, tr_x_yzzzzz_xxxyyz, tr_x_yzzzzz_xxxyzz, tr_x_yzzzzz_xxxzzz, tr_x_yzzzzz_xxyyyy, tr_x_yzzzzz_xxyyyz, tr_x_yzzzzz_xxyyzz, tr_x_yzzzzz_xxyzzz, tr_x_yzzzzz_xxzzzz, tr_x_yzzzzz_xyyyyy, tr_x_yzzzzz_xyyyyz, tr_x_yzzzzz_xyyyzz, tr_x_yzzzzz_xyyzzz, tr_x_yzzzzz_xyzzzz, tr_x_yzzzzz_xzzzzz, tr_x_yzzzzz_yyyyyy, tr_x_yzzzzz_yyyyyz, tr_x_yzzzzz_yyyyzz, tr_x_yzzzzz_yyyzzz, tr_x_yzzzzz_yyzzzz, tr_x_yzzzzz_yzzzzz, tr_x_yzzzzz_zzzzzz, tr_x_zzzzz_xxxxx, tr_x_zzzzz_xxxxxx, tr_x_zzzzz_xxxxxy, tr_x_zzzzz_xxxxxz, tr_x_zzzzz_xxxxy, tr_x_zzzzz_xxxxyy, tr_x_zzzzz_xxxxyz, tr_x_zzzzz_xxxxz, tr_x_zzzzz_xxxxzz, tr_x_zzzzz_xxxyy, tr_x_zzzzz_xxxyyy, tr_x_zzzzz_xxxyyz, tr_x_zzzzz_xxxyz, tr_x_zzzzz_xxxyzz, tr_x_zzzzz_xxxzz, tr_x_zzzzz_xxxzzz, tr_x_zzzzz_xxyyy, tr_x_zzzzz_xxyyyy, tr_x_zzzzz_xxyyyz, tr_x_zzzzz_xxyyz, tr_x_zzzzz_xxyyzz, tr_x_zzzzz_xxyzz, tr_x_zzzzz_xxyzzz, tr_x_zzzzz_xxzzz, tr_x_zzzzz_xxzzzz, tr_x_zzzzz_xyyyy, tr_x_zzzzz_xyyyyy, tr_x_zzzzz_xyyyyz, tr_x_zzzzz_xyyyz, tr_x_zzzzz_xyyyzz, tr_x_zzzzz_xyyzz, tr_x_zzzzz_xyyzzz, tr_x_zzzzz_xyzzz, tr_x_zzzzz_xyzzzz, tr_x_zzzzz_xzzzz, tr_x_zzzzz_xzzzzz, tr_x_zzzzz_yyyyy, tr_x_zzzzz_yyyyyy, tr_x_zzzzz_yyyyyz, tr_x_zzzzz_yyyyz, tr_x_zzzzz_yyyyzz, tr_x_zzzzz_yyyzz, tr_x_zzzzz_yyyzzz, tr_x_zzzzz_yyzzz, tr_x_zzzzz_yyzzzz, tr_x_zzzzz_yzzzz, tr_x_zzzzz_yzzzzz, tr_x_zzzzz_zzzzz, tr_x_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzzz_xxxxxx[i] = tr_x_zzzzz_xxxxxx[i] * pa_y[i];

        tr_x_yzzzzz_xxxxxy[i] = tr_x_zzzzz_xxxxx[i] * fe_0 + tr_x_zzzzz_xxxxxy[i] * pa_y[i];

        tr_x_yzzzzz_xxxxxz[i] = tr_x_zzzzz_xxxxxz[i] * pa_y[i];

        tr_x_yzzzzz_xxxxyy[i] = 2.0 * tr_x_zzzzz_xxxxy[i] * fe_0 + tr_x_zzzzz_xxxxyy[i] * pa_y[i];

        tr_x_yzzzzz_xxxxyz[i] = tr_x_zzzzz_xxxxz[i] * fe_0 + tr_x_zzzzz_xxxxyz[i] * pa_y[i];

        tr_x_yzzzzz_xxxxzz[i] = tr_x_zzzzz_xxxxzz[i] * pa_y[i];

        tr_x_yzzzzz_xxxyyy[i] = 3.0 * tr_x_zzzzz_xxxyy[i] * fe_0 + tr_x_zzzzz_xxxyyy[i] * pa_y[i];

        tr_x_yzzzzz_xxxyyz[i] = 2.0 * tr_x_zzzzz_xxxyz[i] * fe_0 + tr_x_zzzzz_xxxyyz[i] * pa_y[i];

        tr_x_yzzzzz_xxxyzz[i] = tr_x_zzzzz_xxxzz[i] * fe_0 + tr_x_zzzzz_xxxyzz[i] * pa_y[i];

        tr_x_yzzzzz_xxxzzz[i] = tr_x_zzzzz_xxxzzz[i] * pa_y[i];

        tr_x_yzzzzz_xxyyyy[i] = 4.0 * tr_x_zzzzz_xxyyy[i] * fe_0 + tr_x_zzzzz_xxyyyy[i] * pa_y[i];

        tr_x_yzzzzz_xxyyyz[i] = 3.0 * tr_x_zzzzz_xxyyz[i] * fe_0 + tr_x_zzzzz_xxyyyz[i] * pa_y[i];

        tr_x_yzzzzz_xxyyzz[i] = 2.0 * tr_x_zzzzz_xxyzz[i] * fe_0 + tr_x_zzzzz_xxyyzz[i] * pa_y[i];

        tr_x_yzzzzz_xxyzzz[i] = tr_x_zzzzz_xxzzz[i] * fe_0 + tr_x_zzzzz_xxyzzz[i] * pa_y[i];

        tr_x_yzzzzz_xxzzzz[i] = tr_x_zzzzz_xxzzzz[i] * pa_y[i];

        tr_x_yzzzzz_xyyyyy[i] = 5.0 * tr_x_zzzzz_xyyyy[i] * fe_0 + tr_x_zzzzz_xyyyyy[i] * pa_y[i];

        tr_x_yzzzzz_xyyyyz[i] = 4.0 * tr_x_zzzzz_xyyyz[i] * fe_0 + tr_x_zzzzz_xyyyyz[i] * pa_y[i];

        tr_x_yzzzzz_xyyyzz[i] = 3.0 * tr_x_zzzzz_xyyzz[i] * fe_0 + tr_x_zzzzz_xyyyzz[i] * pa_y[i];

        tr_x_yzzzzz_xyyzzz[i] = 2.0 * tr_x_zzzzz_xyzzz[i] * fe_0 + tr_x_zzzzz_xyyzzz[i] * pa_y[i];

        tr_x_yzzzzz_xyzzzz[i] = tr_x_zzzzz_xzzzz[i] * fe_0 + tr_x_zzzzz_xyzzzz[i] * pa_y[i];

        tr_x_yzzzzz_xzzzzz[i] = tr_x_zzzzz_xzzzzz[i] * pa_y[i];

        tr_x_yzzzzz_yyyyyy[i] = 6.0 * tr_x_zzzzz_yyyyy[i] * fe_0 + tr_x_zzzzz_yyyyyy[i] * pa_y[i];

        tr_x_yzzzzz_yyyyyz[i] = 5.0 * tr_x_zzzzz_yyyyz[i] * fe_0 + tr_x_zzzzz_yyyyyz[i] * pa_y[i];

        tr_x_yzzzzz_yyyyzz[i] = 4.0 * tr_x_zzzzz_yyyzz[i] * fe_0 + tr_x_zzzzz_yyyyzz[i] * pa_y[i];

        tr_x_yzzzzz_yyyzzz[i] = 3.0 * tr_x_zzzzz_yyzzz[i] * fe_0 + tr_x_zzzzz_yyyzzz[i] * pa_y[i];

        tr_x_yzzzzz_yyzzzz[i] = 2.0 * tr_x_zzzzz_yzzzz[i] * fe_0 + tr_x_zzzzz_yyzzzz[i] * pa_y[i];

        tr_x_yzzzzz_yzzzzz[i] = tr_x_zzzzz_zzzzz[i] * fe_0 + tr_x_zzzzz_yzzzzz[i] * pa_y[i];

        tr_x_yzzzzz_zzzzzz[i] = tr_x_zzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 756-784 components of targeted buffer : II

    auto tr_x_zzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 756);

    auto tr_x_zzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 757);

    auto tr_x_zzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 758);

    auto tr_x_zzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 759);

    auto tr_x_zzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 760);

    auto tr_x_zzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 761);

    auto tr_x_zzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 762);

    auto tr_x_zzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 763);

    auto tr_x_zzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 764);

    auto tr_x_zzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 765);

    auto tr_x_zzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 766);

    auto tr_x_zzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 767);

    auto tr_x_zzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 768);

    auto tr_x_zzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 769);

    auto tr_x_zzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 770);

    auto tr_x_zzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 771);

    auto tr_x_zzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 772);

    auto tr_x_zzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 773);

    auto tr_x_zzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 774);

    auto tr_x_zzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 775);

    auto tr_x_zzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 776);

    auto tr_x_zzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 777);

    auto tr_x_zzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 778);

    auto tr_x_zzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 779);

    auto tr_x_zzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 780);

    auto tr_x_zzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 781);

    auto tr_x_zzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 782);

    auto tr_x_zzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 783);

    #pragma omp simd aligned(pa_z, tr_x_zzzz_xxxxxx, tr_x_zzzz_xxxxxy, tr_x_zzzz_xxxxxz, tr_x_zzzz_xxxxyy, tr_x_zzzz_xxxxyz, tr_x_zzzz_xxxxzz, tr_x_zzzz_xxxyyy, tr_x_zzzz_xxxyyz, tr_x_zzzz_xxxyzz, tr_x_zzzz_xxxzzz, tr_x_zzzz_xxyyyy, tr_x_zzzz_xxyyyz, tr_x_zzzz_xxyyzz, tr_x_zzzz_xxyzzz, tr_x_zzzz_xxzzzz, tr_x_zzzz_xyyyyy, tr_x_zzzz_xyyyyz, tr_x_zzzz_xyyyzz, tr_x_zzzz_xyyzzz, tr_x_zzzz_xyzzzz, tr_x_zzzz_xzzzzz, tr_x_zzzz_yyyyyy, tr_x_zzzz_yyyyyz, tr_x_zzzz_yyyyzz, tr_x_zzzz_yyyzzz, tr_x_zzzz_yyzzzz, tr_x_zzzz_yzzzzz, tr_x_zzzz_zzzzzz, tr_x_zzzzz_xxxxx, tr_x_zzzzz_xxxxxx, tr_x_zzzzz_xxxxxy, tr_x_zzzzz_xxxxxz, tr_x_zzzzz_xxxxy, tr_x_zzzzz_xxxxyy, tr_x_zzzzz_xxxxyz, tr_x_zzzzz_xxxxz, tr_x_zzzzz_xxxxzz, tr_x_zzzzz_xxxyy, tr_x_zzzzz_xxxyyy, tr_x_zzzzz_xxxyyz, tr_x_zzzzz_xxxyz, tr_x_zzzzz_xxxyzz, tr_x_zzzzz_xxxzz, tr_x_zzzzz_xxxzzz, tr_x_zzzzz_xxyyy, tr_x_zzzzz_xxyyyy, tr_x_zzzzz_xxyyyz, tr_x_zzzzz_xxyyz, tr_x_zzzzz_xxyyzz, tr_x_zzzzz_xxyzz, tr_x_zzzzz_xxyzzz, tr_x_zzzzz_xxzzz, tr_x_zzzzz_xxzzzz, tr_x_zzzzz_xyyyy, tr_x_zzzzz_xyyyyy, tr_x_zzzzz_xyyyyz, tr_x_zzzzz_xyyyz, tr_x_zzzzz_xyyyzz, tr_x_zzzzz_xyyzz, tr_x_zzzzz_xyyzzz, tr_x_zzzzz_xyzzz, tr_x_zzzzz_xyzzzz, tr_x_zzzzz_xzzzz, tr_x_zzzzz_xzzzzz, tr_x_zzzzz_yyyyy, tr_x_zzzzz_yyyyyy, tr_x_zzzzz_yyyyyz, tr_x_zzzzz_yyyyz, tr_x_zzzzz_yyyyzz, tr_x_zzzzz_yyyzz, tr_x_zzzzz_yyyzzz, tr_x_zzzzz_yyzzz, tr_x_zzzzz_yyzzzz, tr_x_zzzzz_yzzzz, tr_x_zzzzz_yzzzzz, tr_x_zzzzz_zzzzz, tr_x_zzzzz_zzzzzz, tr_x_zzzzzz_xxxxxx, tr_x_zzzzzz_xxxxxy, tr_x_zzzzzz_xxxxxz, tr_x_zzzzzz_xxxxyy, tr_x_zzzzzz_xxxxyz, tr_x_zzzzzz_xxxxzz, tr_x_zzzzzz_xxxyyy, tr_x_zzzzzz_xxxyyz, tr_x_zzzzzz_xxxyzz, tr_x_zzzzzz_xxxzzz, tr_x_zzzzzz_xxyyyy, tr_x_zzzzzz_xxyyyz, tr_x_zzzzzz_xxyyzz, tr_x_zzzzzz_xxyzzz, tr_x_zzzzzz_xxzzzz, tr_x_zzzzzz_xyyyyy, tr_x_zzzzzz_xyyyyz, tr_x_zzzzzz_xyyyzz, tr_x_zzzzzz_xyyzzz, tr_x_zzzzzz_xyzzzz, tr_x_zzzzzz_xzzzzz, tr_x_zzzzzz_yyyyyy, tr_x_zzzzzz_yyyyyz, tr_x_zzzzzz_yyyyzz, tr_x_zzzzzz_yyyzzz, tr_x_zzzzzz_yyzzzz, tr_x_zzzzzz_yzzzzz, tr_x_zzzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzzz_xxxxxx[i] = 5.0 * tr_x_zzzz_xxxxxx[i] * fe_0 + tr_x_zzzzz_xxxxxx[i] * pa_z[i];

        tr_x_zzzzzz_xxxxxy[i] = 5.0 * tr_x_zzzz_xxxxxy[i] * fe_0 + tr_x_zzzzz_xxxxxy[i] * pa_z[i];

        tr_x_zzzzzz_xxxxxz[i] = 5.0 * tr_x_zzzz_xxxxxz[i] * fe_0 + tr_x_zzzzz_xxxxx[i] * fe_0 + tr_x_zzzzz_xxxxxz[i] * pa_z[i];

        tr_x_zzzzzz_xxxxyy[i] = 5.0 * tr_x_zzzz_xxxxyy[i] * fe_0 + tr_x_zzzzz_xxxxyy[i] * pa_z[i];

        tr_x_zzzzzz_xxxxyz[i] = 5.0 * tr_x_zzzz_xxxxyz[i] * fe_0 + tr_x_zzzzz_xxxxy[i] * fe_0 + tr_x_zzzzz_xxxxyz[i] * pa_z[i];

        tr_x_zzzzzz_xxxxzz[i] = 5.0 * tr_x_zzzz_xxxxzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xxxxz[i] * fe_0 + tr_x_zzzzz_xxxxzz[i] * pa_z[i];

        tr_x_zzzzzz_xxxyyy[i] = 5.0 * tr_x_zzzz_xxxyyy[i] * fe_0 + tr_x_zzzzz_xxxyyy[i] * pa_z[i];

        tr_x_zzzzzz_xxxyyz[i] = 5.0 * tr_x_zzzz_xxxyyz[i] * fe_0 + tr_x_zzzzz_xxxyy[i] * fe_0 + tr_x_zzzzz_xxxyyz[i] * pa_z[i];

        tr_x_zzzzzz_xxxyzz[i] = 5.0 * tr_x_zzzz_xxxyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xxxyz[i] * fe_0 + tr_x_zzzzz_xxxyzz[i] * pa_z[i];

        tr_x_zzzzzz_xxxzzz[i] = 5.0 * tr_x_zzzz_xxxzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_xxxzz[i] * fe_0 + tr_x_zzzzz_xxxzzz[i] * pa_z[i];

        tr_x_zzzzzz_xxyyyy[i] = 5.0 * tr_x_zzzz_xxyyyy[i] * fe_0 + tr_x_zzzzz_xxyyyy[i] * pa_z[i];

        tr_x_zzzzzz_xxyyyz[i] = 5.0 * tr_x_zzzz_xxyyyz[i] * fe_0 + tr_x_zzzzz_xxyyy[i] * fe_0 + tr_x_zzzzz_xxyyyz[i] * pa_z[i];

        tr_x_zzzzzz_xxyyzz[i] = 5.0 * tr_x_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xxyyz[i] * fe_0 + tr_x_zzzzz_xxyyzz[i] * pa_z[i];

        tr_x_zzzzzz_xxyzzz[i] = 5.0 * tr_x_zzzz_xxyzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_xxyzz[i] * fe_0 + tr_x_zzzzz_xxyzzz[i] * pa_z[i];

        tr_x_zzzzzz_xxzzzz[i] = 5.0 * tr_x_zzzz_xxzzzz[i] * fe_0 + 4.0 * tr_x_zzzzz_xxzzz[i] * fe_0 + tr_x_zzzzz_xxzzzz[i] * pa_z[i];

        tr_x_zzzzzz_xyyyyy[i] = 5.0 * tr_x_zzzz_xyyyyy[i] * fe_0 + tr_x_zzzzz_xyyyyy[i] * pa_z[i];

        tr_x_zzzzzz_xyyyyz[i] = 5.0 * tr_x_zzzz_xyyyyz[i] * fe_0 + tr_x_zzzzz_xyyyy[i] * fe_0 + tr_x_zzzzz_xyyyyz[i] * pa_z[i];

        tr_x_zzzzzz_xyyyzz[i] = 5.0 * tr_x_zzzz_xyyyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xyyyz[i] * fe_0 + tr_x_zzzzz_xyyyzz[i] * pa_z[i];

        tr_x_zzzzzz_xyyzzz[i] = 5.0 * tr_x_zzzz_xyyzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_xyyzz[i] * fe_0 + tr_x_zzzzz_xyyzzz[i] * pa_z[i];

        tr_x_zzzzzz_xyzzzz[i] = 5.0 * tr_x_zzzz_xyzzzz[i] * fe_0 + 4.0 * tr_x_zzzzz_xyzzz[i] * fe_0 + tr_x_zzzzz_xyzzzz[i] * pa_z[i];

        tr_x_zzzzzz_xzzzzz[i] = 5.0 * tr_x_zzzz_xzzzzz[i] * fe_0 + 5.0 * tr_x_zzzzz_xzzzz[i] * fe_0 + tr_x_zzzzz_xzzzzz[i] * pa_z[i];

        tr_x_zzzzzz_yyyyyy[i] = 5.0 * tr_x_zzzz_yyyyyy[i] * fe_0 + tr_x_zzzzz_yyyyyy[i] * pa_z[i];

        tr_x_zzzzzz_yyyyyz[i] = 5.0 * tr_x_zzzz_yyyyyz[i] * fe_0 + tr_x_zzzzz_yyyyy[i] * fe_0 + tr_x_zzzzz_yyyyyz[i] * pa_z[i];

        tr_x_zzzzzz_yyyyzz[i] = 5.0 * tr_x_zzzz_yyyyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_yyyyz[i] * fe_0 + tr_x_zzzzz_yyyyzz[i] * pa_z[i];

        tr_x_zzzzzz_yyyzzz[i] = 5.0 * tr_x_zzzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_yyyzz[i] * fe_0 + tr_x_zzzzz_yyyzzz[i] * pa_z[i];

        tr_x_zzzzzz_yyzzzz[i] = 5.0 * tr_x_zzzz_yyzzzz[i] * fe_0 + 4.0 * tr_x_zzzzz_yyzzz[i] * fe_0 + tr_x_zzzzz_yyzzzz[i] * pa_z[i];

        tr_x_zzzzzz_yzzzzz[i] = 5.0 * tr_x_zzzz_yzzzzz[i] * fe_0 + 5.0 * tr_x_zzzzz_yzzzz[i] * fe_0 + tr_x_zzzzz_yzzzzz[i] * pa_z[i];

        tr_x_zzzzzz_zzzzzz[i] = 5.0 * tr_x_zzzz_zzzzzz[i] * fe_0 + 6.0 * tr_x_zzzzz_zzzzz[i] * fe_0 + tr_x_zzzzz_zzzzzz[i] * pa_z[i];
    }

    // Set up 784-812 components of targeted buffer : II

    auto tr_y_xxxxxx_xxxxxx = pbuffer.data(idx_dip_ii + 784);

    auto tr_y_xxxxxx_xxxxxy = pbuffer.data(idx_dip_ii + 785);

    auto tr_y_xxxxxx_xxxxxz = pbuffer.data(idx_dip_ii + 786);

    auto tr_y_xxxxxx_xxxxyy = pbuffer.data(idx_dip_ii + 787);

    auto tr_y_xxxxxx_xxxxyz = pbuffer.data(idx_dip_ii + 788);

    auto tr_y_xxxxxx_xxxxzz = pbuffer.data(idx_dip_ii + 789);

    auto tr_y_xxxxxx_xxxyyy = pbuffer.data(idx_dip_ii + 790);

    auto tr_y_xxxxxx_xxxyyz = pbuffer.data(idx_dip_ii + 791);

    auto tr_y_xxxxxx_xxxyzz = pbuffer.data(idx_dip_ii + 792);

    auto tr_y_xxxxxx_xxxzzz = pbuffer.data(idx_dip_ii + 793);

    auto tr_y_xxxxxx_xxyyyy = pbuffer.data(idx_dip_ii + 794);

    auto tr_y_xxxxxx_xxyyyz = pbuffer.data(idx_dip_ii + 795);

    auto tr_y_xxxxxx_xxyyzz = pbuffer.data(idx_dip_ii + 796);

    auto tr_y_xxxxxx_xxyzzz = pbuffer.data(idx_dip_ii + 797);

    auto tr_y_xxxxxx_xxzzzz = pbuffer.data(idx_dip_ii + 798);

    auto tr_y_xxxxxx_xyyyyy = pbuffer.data(idx_dip_ii + 799);

    auto tr_y_xxxxxx_xyyyyz = pbuffer.data(idx_dip_ii + 800);

    auto tr_y_xxxxxx_xyyyzz = pbuffer.data(idx_dip_ii + 801);

    auto tr_y_xxxxxx_xyyzzz = pbuffer.data(idx_dip_ii + 802);

    auto tr_y_xxxxxx_xyzzzz = pbuffer.data(idx_dip_ii + 803);

    auto tr_y_xxxxxx_xzzzzz = pbuffer.data(idx_dip_ii + 804);

    auto tr_y_xxxxxx_yyyyyy = pbuffer.data(idx_dip_ii + 805);

    auto tr_y_xxxxxx_yyyyyz = pbuffer.data(idx_dip_ii + 806);

    auto tr_y_xxxxxx_yyyyzz = pbuffer.data(idx_dip_ii + 807);

    auto tr_y_xxxxxx_yyyzzz = pbuffer.data(idx_dip_ii + 808);

    auto tr_y_xxxxxx_yyzzzz = pbuffer.data(idx_dip_ii + 809);

    auto tr_y_xxxxxx_yzzzzz = pbuffer.data(idx_dip_ii + 810);

    auto tr_y_xxxxxx_zzzzzz = pbuffer.data(idx_dip_ii + 811);

    #pragma omp simd aligned(pa_x, tr_y_xxxx_xxxxxx, tr_y_xxxx_xxxxxy, tr_y_xxxx_xxxxxz, tr_y_xxxx_xxxxyy, tr_y_xxxx_xxxxyz, tr_y_xxxx_xxxxzz, tr_y_xxxx_xxxyyy, tr_y_xxxx_xxxyyz, tr_y_xxxx_xxxyzz, tr_y_xxxx_xxxzzz, tr_y_xxxx_xxyyyy, tr_y_xxxx_xxyyyz, tr_y_xxxx_xxyyzz, tr_y_xxxx_xxyzzz, tr_y_xxxx_xxzzzz, tr_y_xxxx_xyyyyy, tr_y_xxxx_xyyyyz, tr_y_xxxx_xyyyzz, tr_y_xxxx_xyyzzz, tr_y_xxxx_xyzzzz, tr_y_xxxx_xzzzzz, tr_y_xxxx_yyyyyy, tr_y_xxxx_yyyyyz, tr_y_xxxx_yyyyzz, tr_y_xxxx_yyyzzz, tr_y_xxxx_yyzzzz, tr_y_xxxx_yzzzzz, tr_y_xxxx_zzzzzz, tr_y_xxxxx_xxxxx, tr_y_xxxxx_xxxxxx, tr_y_xxxxx_xxxxxy, tr_y_xxxxx_xxxxxz, tr_y_xxxxx_xxxxy, tr_y_xxxxx_xxxxyy, tr_y_xxxxx_xxxxyz, tr_y_xxxxx_xxxxz, tr_y_xxxxx_xxxxzz, tr_y_xxxxx_xxxyy, tr_y_xxxxx_xxxyyy, tr_y_xxxxx_xxxyyz, tr_y_xxxxx_xxxyz, tr_y_xxxxx_xxxyzz, tr_y_xxxxx_xxxzz, tr_y_xxxxx_xxxzzz, tr_y_xxxxx_xxyyy, tr_y_xxxxx_xxyyyy, tr_y_xxxxx_xxyyyz, tr_y_xxxxx_xxyyz, tr_y_xxxxx_xxyyzz, tr_y_xxxxx_xxyzz, tr_y_xxxxx_xxyzzz, tr_y_xxxxx_xxzzz, tr_y_xxxxx_xxzzzz, tr_y_xxxxx_xyyyy, tr_y_xxxxx_xyyyyy, tr_y_xxxxx_xyyyyz, tr_y_xxxxx_xyyyz, tr_y_xxxxx_xyyyzz, tr_y_xxxxx_xyyzz, tr_y_xxxxx_xyyzzz, tr_y_xxxxx_xyzzz, tr_y_xxxxx_xyzzzz, tr_y_xxxxx_xzzzz, tr_y_xxxxx_xzzzzz, tr_y_xxxxx_yyyyy, tr_y_xxxxx_yyyyyy, tr_y_xxxxx_yyyyyz, tr_y_xxxxx_yyyyz, tr_y_xxxxx_yyyyzz, tr_y_xxxxx_yyyzz, tr_y_xxxxx_yyyzzz, tr_y_xxxxx_yyzzz, tr_y_xxxxx_yyzzzz, tr_y_xxxxx_yzzzz, tr_y_xxxxx_yzzzzz, tr_y_xxxxx_zzzzz, tr_y_xxxxx_zzzzzz, tr_y_xxxxxx_xxxxxx, tr_y_xxxxxx_xxxxxy, tr_y_xxxxxx_xxxxxz, tr_y_xxxxxx_xxxxyy, tr_y_xxxxxx_xxxxyz, tr_y_xxxxxx_xxxxzz, tr_y_xxxxxx_xxxyyy, tr_y_xxxxxx_xxxyyz, tr_y_xxxxxx_xxxyzz, tr_y_xxxxxx_xxxzzz, tr_y_xxxxxx_xxyyyy, tr_y_xxxxxx_xxyyyz, tr_y_xxxxxx_xxyyzz, tr_y_xxxxxx_xxyzzz, tr_y_xxxxxx_xxzzzz, tr_y_xxxxxx_xyyyyy, tr_y_xxxxxx_xyyyyz, tr_y_xxxxxx_xyyyzz, tr_y_xxxxxx_xyyzzz, tr_y_xxxxxx_xyzzzz, tr_y_xxxxxx_xzzzzz, tr_y_xxxxxx_yyyyyy, tr_y_xxxxxx_yyyyyz, tr_y_xxxxxx_yyyyzz, tr_y_xxxxxx_yyyzzz, tr_y_xxxxxx_yyzzzz, tr_y_xxxxxx_yzzzzz, tr_y_xxxxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxx_xxxxxx[i] = 5.0 * tr_y_xxxx_xxxxxx[i] * fe_0 + 6.0 * tr_y_xxxxx_xxxxx[i] * fe_0 + tr_y_xxxxx_xxxxxx[i] * pa_x[i];

        tr_y_xxxxxx_xxxxxy[i] = 5.0 * tr_y_xxxx_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxxxx_xxxxy[i] * fe_0 + tr_y_xxxxx_xxxxxy[i] * pa_x[i];

        tr_y_xxxxxx_xxxxxz[i] = 5.0 * tr_y_xxxx_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxxxx_xxxxz[i] * fe_0 + tr_y_xxxxx_xxxxxz[i] * pa_x[i];

        tr_y_xxxxxx_xxxxyy[i] = 5.0 * tr_y_xxxx_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxxxx_xxxyy[i] * fe_0 + tr_y_xxxxx_xxxxyy[i] * pa_x[i];

        tr_y_xxxxxx_xxxxyz[i] = 5.0 * tr_y_xxxx_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxxxx_xxxyz[i] * fe_0 + tr_y_xxxxx_xxxxyz[i] * pa_x[i];

        tr_y_xxxxxx_xxxxzz[i] = 5.0 * tr_y_xxxx_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxxxx_xxxzz[i] * fe_0 + tr_y_xxxxx_xxxxzz[i] * pa_x[i];

        tr_y_xxxxxx_xxxyyy[i] = 5.0 * tr_y_xxxx_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxxxx_xxyyy[i] * fe_0 + tr_y_xxxxx_xxxyyy[i] * pa_x[i];

        tr_y_xxxxxx_xxxyyz[i] = 5.0 * tr_y_xxxx_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxxxx_xxyyz[i] * fe_0 + tr_y_xxxxx_xxxyyz[i] * pa_x[i];

        tr_y_xxxxxx_xxxyzz[i] = 5.0 * tr_y_xxxx_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxxxx_xxyzz[i] * fe_0 + tr_y_xxxxx_xxxyzz[i] * pa_x[i];

        tr_y_xxxxxx_xxxzzz[i] = 5.0 * tr_y_xxxx_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxxxx_xxzzz[i] * fe_0 + tr_y_xxxxx_xxxzzz[i] * pa_x[i];

        tr_y_xxxxxx_xxyyyy[i] = 5.0 * tr_y_xxxx_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxxxx_xyyyy[i] * fe_0 + tr_y_xxxxx_xxyyyy[i] * pa_x[i];

        tr_y_xxxxxx_xxyyyz[i] = 5.0 * tr_y_xxxx_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxxxx_xyyyz[i] * fe_0 + tr_y_xxxxx_xxyyyz[i] * pa_x[i];

        tr_y_xxxxxx_xxyyzz[i] = 5.0 * tr_y_xxxx_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxxxx_xyyzz[i] * fe_0 + tr_y_xxxxx_xxyyzz[i] * pa_x[i];

        tr_y_xxxxxx_xxyzzz[i] = 5.0 * tr_y_xxxx_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxxxx_xyzzz[i] * fe_0 + tr_y_xxxxx_xxyzzz[i] * pa_x[i];

        tr_y_xxxxxx_xxzzzz[i] = 5.0 * tr_y_xxxx_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxxxx_xzzzz[i] * fe_0 + tr_y_xxxxx_xxzzzz[i] * pa_x[i];

        tr_y_xxxxxx_xyyyyy[i] = 5.0 * tr_y_xxxx_xyyyyy[i] * fe_0 + tr_y_xxxxx_yyyyy[i] * fe_0 + tr_y_xxxxx_xyyyyy[i] * pa_x[i];

        tr_y_xxxxxx_xyyyyz[i] = 5.0 * tr_y_xxxx_xyyyyz[i] * fe_0 + tr_y_xxxxx_yyyyz[i] * fe_0 + tr_y_xxxxx_xyyyyz[i] * pa_x[i];

        tr_y_xxxxxx_xyyyzz[i] = 5.0 * tr_y_xxxx_xyyyzz[i] * fe_0 + tr_y_xxxxx_yyyzz[i] * fe_0 + tr_y_xxxxx_xyyyzz[i] * pa_x[i];

        tr_y_xxxxxx_xyyzzz[i] = 5.0 * tr_y_xxxx_xyyzzz[i] * fe_0 + tr_y_xxxxx_yyzzz[i] * fe_0 + tr_y_xxxxx_xyyzzz[i] * pa_x[i];

        tr_y_xxxxxx_xyzzzz[i] = 5.0 * tr_y_xxxx_xyzzzz[i] * fe_0 + tr_y_xxxxx_yzzzz[i] * fe_0 + tr_y_xxxxx_xyzzzz[i] * pa_x[i];

        tr_y_xxxxxx_xzzzzz[i] = 5.0 * tr_y_xxxx_xzzzzz[i] * fe_0 + tr_y_xxxxx_zzzzz[i] * fe_0 + tr_y_xxxxx_xzzzzz[i] * pa_x[i];

        tr_y_xxxxxx_yyyyyy[i] = 5.0 * tr_y_xxxx_yyyyyy[i] * fe_0 + tr_y_xxxxx_yyyyyy[i] * pa_x[i];

        tr_y_xxxxxx_yyyyyz[i] = 5.0 * tr_y_xxxx_yyyyyz[i] * fe_0 + tr_y_xxxxx_yyyyyz[i] * pa_x[i];

        tr_y_xxxxxx_yyyyzz[i] = 5.0 * tr_y_xxxx_yyyyzz[i] * fe_0 + tr_y_xxxxx_yyyyzz[i] * pa_x[i];

        tr_y_xxxxxx_yyyzzz[i] = 5.0 * tr_y_xxxx_yyyzzz[i] * fe_0 + tr_y_xxxxx_yyyzzz[i] * pa_x[i];

        tr_y_xxxxxx_yyzzzz[i] = 5.0 * tr_y_xxxx_yyzzzz[i] * fe_0 + tr_y_xxxxx_yyzzzz[i] * pa_x[i];

        tr_y_xxxxxx_yzzzzz[i] = 5.0 * tr_y_xxxx_yzzzzz[i] * fe_0 + tr_y_xxxxx_yzzzzz[i] * pa_x[i];

        tr_y_xxxxxx_zzzzzz[i] = 5.0 * tr_y_xxxx_zzzzzz[i] * fe_0 + tr_y_xxxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 812-840 components of targeted buffer : II

    auto tr_y_xxxxxy_xxxxxx = pbuffer.data(idx_dip_ii + 812);

    auto tr_y_xxxxxy_xxxxxy = pbuffer.data(idx_dip_ii + 813);

    auto tr_y_xxxxxy_xxxxxz = pbuffer.data(idx_dip_ii + 814);

    auto tr_y_xxxxxy_xxxxyy = pbuffer.data(idx_dip_ii + 815);

    auto tr_y_xxxxxy_xxxxyz = pbuffer.data(idx_dip_ii + 816);

    auto tr_y_xxxxxy_xxxxzz = pbuffer.data(idx_dip_ii + 817);

    auto tr_y_xxxxxy_xxxyyy = pbuffer.data(idx_dip_ii + 818);

    auto tr_y_xxxxxy_xxxyyz = pbuffer.data(idx_dip_ii + 819);

    auto tr_y_xxxxxy_xxxyzz = pbuffer.data(idx_dip_ii + 820);

    auto tr_y_xxxxxy_xxxzzz = pbuffer.data(idx_dip_ii + 821);

    auto tr_y_xxxxxy_xxyyyy = pbuffer.data(idx_dip_ii + 822);

    auto tr_y_xxxxxy_xxyyyz = pbuffer.data(idx_dip_ii + 823);

    auto tr_y_xxxxxy_xxyyzz = pbuffer.data(idx_dip_ii + 824);

    auto tr_y_xxxxxy_xxyzzz = pbuffer.data(idx_dip_ii + 825);

    auto tr_y_xxxxxy_xxzzzz = pbuffer.data(idx_dip_ii + 826);

    auto tr_y_xxxxxy_xyyyyy = pbuffer.data(idx_dip_ii + 827);

    auto tr_y_xxxxxy_xyyyyz = pbuffer.data(idx_dip_ii + 828);

    auto tr_y_xxxxxy_xyyyzz = pbuffer.data(idx_dip_ii + 829);

    auto tr_y_xxxxxy_xyyzzz = pbuffer.data(idx_dip_ii + 830);

    auto tr_y_xxxxxy_xyzzzz = pbuffer.data(idx_dip_ii + 831);

    auto tr_y_xxxxxy_xzzzzz = pbuffer.data(idx_dip_ii + 832);

    auto tr_y_xxxxxy_yyyyyy = pbuffer.data(idx_dip_ii + 833);

    auto tr_y_xxxxxy_yyyyyz = pbuffer.data(idx_dip_ii + 834);

    auto tr_y_xxxxxy_yyyyzz = pbuffer.data(idx_dip_ii + 835);

    auto tr_y_xxxxxy_yyyzzz = pbuffer.data(idx_dip_ii + 836);

    auto tr_y_xxxxxy_yyzzzz = pbuffer.data(idx_dip_ii + 837);

    auto tr_y_xxxxxy_yzzzzz = pbuffer.data(idx_dip_ii + 838);

    auto tr_y_xxxxxy_zzzzzz = pbuffer.data(idx_dip_ii + 839);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xxxxx_xxxxxx, tr_y_xxxxx_xxxxxz, tr_y_xxxxx_xxxxzz, tr_y_xxxxx_xxxzzz, tr_y_xxxxx_xxzzzz, tr_y_xxxxx_xzzzzz, tr_y_xxxxxy_xxxxxx, tr_y_xxxxxy_xxxxxy, tr_y_xxxxxy_xxxxxz, tr_y_xxxxxy_xxxxyy, tr_y_xxxxxy_xxxxyz, tr_y_xxxxxy_xxxxzz, tr_y_xxxxxy_xxxyyy, tr_y_xxxxxy_xxxyyz, tr_y_xxxxxy_xxxyzz, tr_y_xxxxxy_xxxzzz, tr_y_xxxxxy_xxyyyy, tr_y_xxxxxy_xxyyyz, tr_y_xxxxxy_xxyyzz, tr_y_xxxxxy_xxyzzz, tr_y_xxxxxy_xxzzzz, tr_y_xxxxxy_xyyyyy, tr_y_xxxxxy_xyyyyz, tr_y_xxxxxy_xyyyzz, tr_y_xxxxxy_xyyzzz, tr_y_xxxxxy_xyzzzz, tr_y_xxxxxy_xzzzzz, tr_y_xxxxxy_yyyyyy, tr_y_xxxxxy_yyyyyz, tr_y_xxxxxy_yyyyzz, tr_y_xxxxxy_yyyzzz, tr_y_xxxxxy_yyzzzz, tr_y_xxxxxy_yzzzzz, tr_y_xxxxxy_zzzzzz, tr_y_xxxxy_xxxxxy, tr_y_xxxxy_xxxxy, tr_y_xxxxy_xxxxyy, tr_y_xxxxy_xxxxyz, tr_y_xxxxy_xxxyy, tr_y_xxxxy_xxxyyy, tr_y_xxxxy_xxxyyz, tr_y_xxxxy_xxxyz, tr_y_xxxxy_xxxyzz, tr_y_xxxxy_xxyyy, tr_y_xxxxy_xxyyyy, tr_y_xxxxy_xxyyyz, tr_y_xxxxy_xxyyz, tr_y_xxxxy_xxyyzz, tr_y_xxxxy_xxyzz, tr_y_xxxxy_xxyzzz, tr_y_xxxxy_xyyyy, tr_y_xxxxy_xyyyyy, tr_y_xxxxy_xyyyyz, tr_y_xxxxy_xyyyz, tr_y_xxxxy_xyyyzz, tr_y_xxxxy_xyyzz, tr_y_xxxxy_xyyzzz, tr_y_xxxxy_xyzzz, tr_y_xxxxy_xyzzzz, tr_y_xxxxy_yyyyy, tr_y_xxxxy_yyyyyy, tr_y_xxxxy_yyyyyz, tr_y_xxxxy_yyyyz, tr_y_xxxxy_yyyyzz, tr_y_xxxxy_yyyzz, tr_y_xxxxy_yyyzzz, tr_y_xxxxy_yyzzz, tr_y_xxxxy_yyzzzz, tr_y_xxxxy_yzzzz, tr_y_xxxxy_yzzzzz, tr_y_xxxxy_zzzzzz, tr_y_xxxy_xxxxxy, tr_y_xxxy_xxxxyy, tr_y_xxxy_xxxxyz, tr_y_xxxy_xxxyyy, tr_y_xxxy_xxxyyz, tr_y_xxxy_xxxyzz, tr_y_xxxy_xxyyyy, tr_y_xxxy_xxyyyz, tr_y_xxxy_xxyyzz, tr_y_xxxy_xxyzzz, tr_y_xxxy_xyyyyy, tr_y_xxxy_xyyyyz, tr_y_xxxy_xyyyzz, tr_y_xxxy_xyyzzz, tr_y_xxxy_xyzzzz, tr_y_xxxy_yyyyyy, tr_y_xxxy_yyyyyz, tr_y_xxxy_yyyyzz, tr_y_xxxy_yyyzzz, tr_y_xxxy_yyzzzz, tr_y_xxxy_yzzzzz, tr_y_xxxy_zzzzzz, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxz, ts_xxxxx_xxxxzz, ts_xxxxx_xxxzzz, ts_xxxxx_xxzzzz, ts_xxxxx_xzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxy_xxxxxx[i] = ts_xxxxx_xxxxxx[i] * fe_0 + tr_y_xxxxx_xxxxxx[i] * pa_y[i];

        tr_y_xxxxxy_xxxxxy[i] = 4.0 * tr_y_xxxy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxxxy_xxxxy[i] * fe_0 + tr_y_xxxxy_xxxxxy[i] * pa_x[i];

        tr_y_xxxxxy_xxxxxz[i] = ts_xxxxx_xxxxxz[i] * fe_0 + tr_y_xxxxx_xxxxxz[i] * pa_y[i];

        tr_y_xxxxxy_xxxxyy[i] = 4.0 * tr_y_xxxy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxxxy_xxxyy[i] * fe_0 + tr_y_xxxxy_xxxxyy[i] * pa_x[i];

        tr_y_xxxxxy_xxxxyz[i] = 4.0 * tr_y_xxxy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxxxy_xxxyz[i] * fe_0 + tr_y_xxxxy_xxxxyz[i] * pa_x[i];

        tr_y_xxxxxy_xxxxzz[i] = ts_xxxxx_xxxxzz[i] * fe_0 + tr_y_xxxxx_xxxxzz[i] * pa_y[i];

        tr_y_xxxxxy_xxxyyy[i] = 4.0 * tr_y_xxxy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxxxy_xxyyy[i] * fe_0 + tr_y_xxxxy_xxxyyy[i] * pa_x[i];

        tr_y_xxxxxy_xxxyyz[i] = 4.0 * tr_y_xxxy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxxxy_xxyyz[i] * fe_0 + tr_y_xxxxy_xxxyyz[i] * pa_x[i];

        tr_y_xxxxxy_xxxyzz[i] = 4.0 * tr_y_xxxy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxxxy_xxyzz[i] * fe_0 + tr_y_xxxxy_xxxyzz[i] * pa_x[i];

        tr_y_xxxxxy_xxxzzz[i] = ts_xxxxx_xxxzzz[i] * fe_0 + tr_y_xxxxx_xxxzzz[i] * pa_y[i];

        tr_y_xxxxxy_xxyyyy[i] = 4.0 * tr_y_xxxy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxxxy_xyyyy[i] * fe_0 + tr_y_xxxxy_xxyyyy[i] * pa_x[i];

        tr_y_xxxxxy_xxyyyz[i] = 4.0 * tr_y_xxxy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxxxy_xyyyz[i] * fe_0 + tr_y_xxxxy_xxyyyz[i] * pa_x[i];

        tr_y_xxxxxy_xxyyzz[i] = 4.0 * tr_y_xxxy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxxxy_xyyzz[i] * fe_0 + tr_y_xxxxy_xxyyzz[i] * pa_x[i];

        tr_y_xxxxxy_xxyzzz[i] = 4.0 * tr_y_xxxy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxxxy_xyzzz[i] * fe_0 + tr_y_xxxxy_xxyzzz[i] * pa_x[i];

        tr_y_xxxxxy_xxzzzz[i] = ts_xxxxx_xxzzzz[i] * fe_0 + tr_y_xxxxx_xxzzzz[i] * pa_y[i];

        tr_y_xxxxxy_xyyyyy[i] = 4.0 * tr_y_xxxy_xyyyyy[i] * fe_0 + tr_y_xxxxy_yyyyy[i] * fe_0 + tr_y_xxxxy_xyyyyy[i] * pa_x[i];

        tr_y_xxxxxy_xyyyyz[i] = 4.0 * tr_y_xxxy_xyyyyz[i] * fe_0 + tr_y_xxxxy_yyyyz[i] * fe_0 + tr_y_xxxxy_xyyyyz[i] * pa_x[i];

        tr_y_xxxxxy_xyyyzz[i] = 4.0 * tr_y_xxxy_xyyyzz[i] * fe_0 + tr_y_xxxxy_yyyzz[i] * fe_0 + tr_y_xxxxy_xyyyzz[i] * pa_x[i];

        tr_y_xxxxxy_xyyzzz[i] = 4.0 * tr_y_xxxy_xyyzzz[i] * fe_0 + tr_y_xxxxy_yyzzz[i] * fe_0 + tr_y_xxxxy_xyyzzz[i] * pa_x[i];

        tr_y_xxxxxy_xyzzzz[i] = 4.0 * tr_y_xxxy_xyzzzz[i] * fe_0 + tr_y_xxxxy_yzzzz[i] * fe_0 + tr_y_xxxxy_xyzzzz[i] * pa_x[i];

        tr_y_xxxxxy_xzzzzz[i] = ts_xxxxx_xzzzzz[i] * fe_0 + tr_y_xxxxx_xzzzzz[i] * pa_y[i];

        tr_y_xxxxxy_yyyyyy[i] = 4.0 * tr_y_xxxy_yyyyyy[i] * fe_0 + tr_y_xxxxy_yyyyyy[i] * pa_x[i];

        tr_y_xxxxxy_yyyyyz[i] = 4.0 * tr_y_xxxy_yyyyyz[i] * fe_0 + tr_y_xxxxy_yyyyyz[i] * pa_x[i];

        tr_y_xxxxxy_yyyyzz[i] = 4.0 * tr_y_xxxy_yyyyzz[i] * fe_0 + tr_y_xxxxy_yyyyzz[i] * pa_x[i];

        tr_y_xxxxxy_yyyzzz[i] = 4.0 * tr_y_xxxy_yyyzzz[i] * fe_0 + tr_y_xxxxy_yyyzzz[i] * pa_x[i];

        tr_y_xxxxxy_yyzzzz[i] = 4.0 * tr_y_xxxy_yyzzzz[i] * fe_0 + tr_y_xxxxy_yyzzzz[i] * pa_x[i];

        tr_y_xxxxxy_yzzzzz[i] = 4.0 * tr_y_xxxy_yzzzzz[i] * fe_0 + tr_y_xxxxy_yzzzzz[i] * pa_x[i];

        tr_y_xxxxxy_zzzzzz[i] = 4.0 * tr_y_xxxy_zzzzzz[i] * fe_0 + tr_y_xxxxy_zzzzzz[i] * pa_x[i];
    }

    // Set up 840-868 components of targeted buffer : II

    auto tr_y_xxxxxz_xxxxxx = pbuffer.data(idx_dip_ii + 840);

    auto tr_y_xxxxxz_xxxxxy = pbuffer.data(idx_dip_ii + 841);

    auto tr_y_xxxxxz_xxxxxz = pbuffer.data(idx_dip_ii + 842);

    auto tr_y_xxxxxz_xxxxyy = pbuffer.data(idx_dip_ii + 843);

    auto tr_y_xxxxxz_xxxxyz = pbuffer.data(idx_dip_ii + 844);

    auto tr_y_xxxxxz_xxxxzz = pbuffer.data(idx_dip_ii + 845);

    auto tr_y_xxxxxz_xxxyyy = pbuffer.data(idx_dip_ii + 846);

    auto tr_y_xxxxxz_xxxyyz = pbuffer.data(idx_dip_ii + 847);

    auto tr_y_xxxxxz_xxxyzz = pbuffer.data(idx_dip_ii + 848);

    auto tr_y_xxxxxz_xxxzzz = pbuffer.data(idx_dip_ii + 849);

    auto tr_y_xxxxxz_xxyyyy = pbuffer.data(idx_dip_ii + 850);

    auto tr_y_xxxxxz_xxyyyz = pbuffer.data(idx_dip_ii + 851);

    auto tr_y_xxxxxz_xxyyzz = pbuffer.data(idx_dip_ii + 852);

    auto tr_y_xxxxxz_xxyzzz = pbuffer.data(idx_dip_ii + 853);

    auto tr_y_xxxxxz_xxzzzz = pbuffer.data(idx_dip_ii + 854);

    auto tr_y_xxxxxz_xyyyyy = pbuffer.data(idx_dip_ii + 855);

    auto tr_y_xxxxxz_xyyyyz = pbuffer.data(idx_dip_ii + 856);

    auto tr_y_xxxxxz_xyyyzz = pbuffer.data(idx_dip_ii + 857);

    auto tr_y_xxxxxz_xyyzzz = pbuffer.data(idx_dip_ii + 858);

    auto tr_y_xxxxxz_xyzzzz = pbuffer.data(idx_dip_ii + 859);

    auto tr_y_xxxxxz_xzzzzz = pbuffer.data(idx_dip_ii + 860);

    auto tr_y_xxxxxz_yyyyyy = pbuffer.data(idx_dip_ii + 861);

    auto tr_y_xxxxxz_yyyyyz = pbuffer.data(idx_dip_ii + 862);

    auto tr_y_xxxxxz_yyyyzz = pbuffer.data(idx_dip_ii + 863);

    auto tr_y_xxxxxz_yyyzzz = pbuffer.data(idx_dip_ii + 864);

    auto tr_y_xxxxxz_yyzzzz = pbuffer.data(idx_dip_ii + 865);

    auto tr_y_xxxxxz_yzzzzz = pbuffer.data(idx_dip_ii + 866);

    auto tr_y_xxxxxz_zzzzzz = pbuffer.data(idx_dip_ii + 867);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxxx_xxxxx, tr_y_xxxxx_xxxxxx, tr_y_xxxxx_xxxxxy, tr_y_xxxxx_xxxxxz, tr_y_xxxxx_xxxxy, tr_y_xxxxx_xxxxyy, tr_y_xxxxx_xxxxyz, tr_y_xxxxx_xxxxz, tr_y_xxxxx_xxxxzz, tr_y_xxxxx_xxxyy, tr_y_xxxxx_xxxyyy, tr_y_xxxxx_xxxyyz, tr_y_xxxxx_xxxyz, tr_y_xxxxx_xxxyzz, tr_y_xxxxx_xxxzz, tr_y_xxxxx_xxxzzz, tr_y_xxxxx_xxyyy, tr_y_xxxxx_xxyyyy, tr_y_xxxxx_xxyyyz, tr_y_xxxxx_xxyyz, tr_y_xxxxx_xxyyzz, tr_y_xxxxx_xxyzz, tr_y_xxxxx_xxyzzz, tr_y_xxxxx_xxzzz, tr_y_xxxxx_xxzzzz, tr_y_xxxxx_xyyyy, tr_y_xxxxx_xyyyyy, tr_y_xxxxx_xyyyyz, tr_y_xxxxx_xyyyz, tr_y_xxxxx_xyyyzz, tr_y_xxxxx_xyyzz, tr_y_xxxxx_xyyzzz, tr_y_xxxxx_xyzzz, tr_y_xxxxx_xyzzzz, tr_y_xxxxx_xzzzz, tr_y_xxxxx_xzzzzz, tr_y_xxxxx_yyyyyy, tr_y_xxxxxz_xxxxxx, tr_y_xxxxxz_xxxxxy, tr_y_xxxxxz_xxxxxz, tr_y_xxxxxz_xxxxyy, tr_y_xxxxxz_xxxxyz, tr_y_xxxxxz_xxxxzz, tr_y_xxxxxz_xxxyyy, tr_y_xxxxxz_xxxyyz, tr_y_xxxxxz_xxxyzz, tr_y_xxxxxz_xxxzzz, tr_y_xxxxxz_xxyyyy, tr_y_xxxxxz_xxyyyz, tr_y_xxxxxz_xxyyzz, tr_y_xxxxxz_xxyzzz, tr_y_xxxxxz_xxzzzz, tr_y_xxxxxz_xyyyyy, tr_y_xxxxxz_xyyyyz, tr_y_xxxxxz_xyyyzz, tr_y_xxxxxz_xyyzzz, tr_y_xxxxxz_xyzzzz, tr_y_xxxxxz_xzzzzz, tr_y_xxxxxz_yyyyyy, tr_y_xxxxxz_yyyyyz, tr_y_xxxxxz_yyyyzz, tr_y_xxxxxz_yyyzzz, tr_y_xxxxxz_yyzzzz, tr_y_xxxxxz_yzzzzz, tr_y_xxxxxz_zzzzzz, tr_y_xxxxz_yyyyyz, tr_y_xxxxz_yyyyzz, tr_y_xxxxz_yyyzzz, tr_y_xxxxz_yyzzzz, tr_y_xxxxz_yzzzzz, tr_y_xxxxz_zzzzzz, tr_y_xxxz_yyyyyz, tr_y_xxxz_yyyyzz, tr_y_xxxz_yyyzzz, tr_y_xxxz_yyzzzz, tr_y_xxxz_yzzzzz, tr_y_xxxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxz_xxxxxx[i] = tr_y_xxxxx_xxxxxx[i] * pa_z[i];

        tr_y_xxxxxz_xxxxxy[i] = tr_y_xxxxx_xxxxxy[i] * pa_z[i];

        tr_y_xxxxxz_xxxxxz[i] = tr_y_xxxxx_xxxxx[i] * fe_0 + tr_y_xxxxx_xxxxxz[i] * pa_z[i];

        tr_y_xxxxxz_xxxxyy[i] = tr_y_xxxxx_xxxxyy[i] * pa_z[i];

        tr_y_xxxxxz_xxxxyz[i] = tr_y_xxxxx_xxxxy[i] * fe_0 + tr_y_xxxxx_xxxxyz[i] * pa_z[i];

        tr_y_xxxxxz_xxxxzz[i] = 2.0 * tr_y_xxxxx_xxxxz[i] * fe_0 + tr_y_xxxxx_xxxxzz[i] * pa_z[i];

        tr_y_xxxxxz_xxxyyy[i] = tr_y_xxxxx_xxxyyy[i] * pa_z[i];

        tr_y_xxxxxz_xxxyyz[i] = tr_y_xxxxx_xxxyy[i] * fe_0 + tr_y_xxxxx_xxxyyz[i] * pa_z[i];

        tr_y_xxxxxz_xxxyzz[i] = 2.0 * tr_y_xxxxx_xxxyz[i] * fe_0 + tr_y_xxxxx_xxxyzz[i] * pa_z[i];

        tr_y_xxxxxz_xxxzzz[i] = 3.0 * tr_y_xxxxx_xxxzz[i] * fe_0 + tr_y_xxxxx_xxxzzz[i] * pa_z[i];

        tr_y_xxxxxz_xxyyyy[i] = tr_y_xxxxx_xxyyyy[i] * pa_z[i];

        tr_y_xxxxxz_xxyyyz[i] = tr_y_xxxxx_xxyyy[i] * fe_0 + tr_y_xxxxx_xxyyyz[i] * pa_z[i];

        tr_y_xxxxxz_xxyyzz[i] = 2.0 * tr_y_xxxxx_xxyyz[i] * fe_0 + tr_y_xxxxx_xxyyzz[i] * pa_z[i];

        tr_y_xxxxxz_xxyzzz[i] = 3.0 * tr_y_xxxxx_xxyzz[i] * fe_0 + tr_y_xxxxx_xxyzzz[i] * pa_z[i];

        tr_y_xxxxxz_xxzzzz[i] = 4.0 * tr_y_xxxxx_xxzzz[i] * fe_0 + tr_y_xxxxx_xxzzzz[i] * pa_z[i];

        tr_y_xxxxxz_xyyyyy[i] = tr_y_xxxxx_xyyyyy[i] * pa_z[i];

        tr_y_xxxxxz_xyyyyz[i] = tr_y_xxxxx_xyyyy[i] * fe_0 + tr_y_xxxxx_xyyyyz[i] * pa_z[i];

        tr_y_xxxxxz_xyyyzz[i] = 2.0 * tr_y_xxxxx_xyyyz[i] * fe_0 + tr_y_xxxxx_xyyyzz[i] * pa_z[i];

        tr_y_xxxxxz_xyyzzz[i] = 3.0 * tr_y_xxxxx_xyyzz[i] * fe_0 + tr_y_xxxxx_xyyzzz[i] * pa_z[i];

        tr_y_xxxxxz_xyzzzz[i] = 4.0 * tr_y_xxxxx_xyzzz[i] * fe_0 + tr_y_xxxxx_xyzzzz[i] * pa_z[i];

        tr_y_xxxxxz_xzzzzz[i] = 5.0 * tr_y_xxxxx_xzzzz[i] * fe_0 + tr_y_xxxxx_xzzzzz[i] * pa_z[i];

        tr_y_xxxxxz_yyyyyy[i] = tr_y_xxxxx_yyyyyy[i] * pa_z[i];

        tr_y_xxxxxz_yyyyyz[i] = 4.0 * tr_y_xxxz_yyyyyz[i] * fe_0 + tr_y_xxxxz_yyyyyz[i] * pa_x[i];

        tr_y_xxxxxz_yyyyzz[i] = 4.0 * tr_y_xxxz_yyyyzz[i] * fe_0 + tr_y_xxxxz_yyyyzz[i] * pa_x[i];

        tr_y_xxxxxz_yyyzzz[i] = 4.0 * tr_y_xxxz_yyyzzz[i] * fe_0 + tr_y_xxxxz_yyyzzz[i] * pa_x[i];

        tr_y_xxxxxz_yyzzzz[i] = 4.0 * tr_y_xxxz_yyzzzz[i] * fe_0 + tr_y_xxxxz_yyzzzz[i] * pa_x[i];

        tr_y_xxxxxz_yzzzzz[i] = 4.0 * tr_y_xxxz_yzzzzz[i] * fe_0 + tr_y_xxxxz_yzzzzz[i] * pa_x[i];

        tr_y_xxxxxz_zzzzzz[i] = 4.0 * tr_y_xxxz_zzzzzz[i] * fe_0 + tr_y_xxxxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 868-896 components of targeted buffer : II

    auto tr_y_xxxxyy_xxxxxx = pbuffer.data(idx_dip_ii + 868);

    auto tr_y_xxxxyy_xxxxxy = pbuffer.data(idx_dip_ii + 869);

    auto tr_y_xxxxyy_xxxxxz = pbuffer.data(idx_dip_ii + 870);

    auto tr_y_xxxxyy_xxxxyy = pbuffer.data(idx_dip_ii + 871);

    auto tr_y_xxxxyy_xxxxyz = pbuffer.data(idx_dip_ii + 872);

    auto tr_y_xxxxyy_xxxxzz = pbuffer.data(idx_dip_ii + 873);

    auto tr_y_xxxxyy_xxxyyy = pbuffer.data(idx_dip_ii + 874);

    auto tr_y_xxxxyy_xxxyyz = pbuffer.data(idx_dip_ii + 875);

    auto tr_y_xxxxyy_xxxyzz = pbuffer.data(idx_dip_ii + 876);

    auto tr_y_xxxxyy_xxxzzz = pbuffer.data(idx_dip_ii + 877);

    auto tr_y_xxxxyy_xxyyyy = pbuffer.data(idx_dip_ii + 878);

    auto tr_y_xxxxyy_xxyyyz = pbuffer.data(idx_dip_ii + 879);

    auto tr_y_xxxxyy_xxyyzz = pbuffer.data(idx_dip_ii + 880);

    auto tr_y_xxxxyy_xxyzzz = pbuffer.data(idx_dip_ii + 881);

    auto tr_y_xxxxyy_xxzzzz = pbuffer.data(idx_dip_ii + 882);

    auto tr_y_xxxxyy_xyyyyy = pbuffer.data(idx_dip_ii + 883);

    auto tr_y_xxxxyy_xyyyyz = pbuffer.data(idx_dip_ii + 884);

    auto tr_y_xxxxyy_xyyyzz = pbuffer.data(idx_dip_ii + 885);

    auto tr_y_xxxxyy_xyyzzz = pbuffer.data(idx_dip_ii + 886);

    auto tr_y_xxxxyy_xyzzzz = pbuffer.data(idx_dip_ii + 887);

    auto tr_y_xxxxyy_xzzzzz = pbuffer.data(idx_dip_ii + 888);

    auto tr_y_xxxxyy_yyyyyy = pbuffer.data(idx_dip_ii + 889);

    auto tr_y_xxxxyy_yyyyyz = pbuffer.data(idx_dip_ii + 890);

    auto tr_y_xxxxyy_yyyyzz = pbuffer.data(idx_dip_ii + 891);

    auto tr_y_xxxxyy_yyyzzz = pbuffer.data(idx_dip_ii + 892);

    auto tr_y_xxxxyy_yyzzzz = pbuffer.data(idx_dip_ii + 893);

    auto tr_y_xxxxyy_yzzzzz = pbuffer.data(idx_dip_ii + 894);

    auto tr_y_xxxxyy_zzzzzz = pbuffer.data(idx_dip_ii + 895);

    #pragma omp simd aligned(pa_x, tr_y_xxxxyy_xxxxxx, tr_y_xxxxyy_xxxxxy, tr_y_xxxxyy_xxxxxz, tr_y_xxxxyy_xxxxyy, tr_y_xxxxyy_xxxxyz, tr_y_xxxxyy_xxxxzz, tr_y_xxxxyy_xxxyyy, tr_y_xxxxyy_xxxyyz, tr_y_xxxxyy_xxxyzz, tr_y_xxxxyy_xxxzzz, tr_y_xxxxyy_xxyyyy, tr_y_xxxxyy_xxyyyz, tr_y_xxxxyy_xxyyzz, tr_y_xxxxyy_xxyzzz, tr_y_xxxxyy_xxzzzz, tr_y_xxxxyy_xyyyyy, tr_y_xxxxyy_xyyyyz, tr_y_xxxxyy_xyyyzz, tr_y_xxxxyy_xyyzzz, tr_y_xxxxyy_xyzzzz, tr_y_xxxxyy_xzzzzz, tr_y_xxxxyy_yyyyyy, tr_y_xxxxyy_yyyyyz, tr_y_xxxxyy_yyyyzz, tr_y_xxxxyy_yyyzzz, tr_y_xxxxyy_yyzzzz, tr_y_xxxxyy_yzzzzz, tr_y_xxxxyy_zzzzzz, tr_y_xxxyy_xxxxx, tr_y_xxxyy_xxxxxx, tr_y_xxxyy_xxxxxy, tr_y_xxxyy_xxxxxz, tr_y_xxxyy_xxxxy, tr_y_xxxyy_xxxxyy, tr_y_xxxyy_xxxxyz, tr_y_xxxyy_xxxxz, tr_y_xxxyy_xxxxzz, tr_y_xxxyy_xxxyy, tr_y_xxxyy_xxxyyy, tr_y_xxxyy_xxxyyz, tr_y_xxxyy_xxxyz, tr_y_xxxyy_xxxyzz, tr_y_xxxyy_xxxzz, tr_y_xxxyy_xxxzzz, tr_y_xxxyy_xxyyy, tr_y_xxxyy_xxyyyy, tr_y_xxxyy_xxyyyz, tr_y_xxxyy_xxyyz, tr_y_xxxyy_xxyyzz, tr_y_xxxyy_xxyzz, tr_y_xxxyy_xxyzzz, tr_y_xxxyy_xxzzz, tr_y_xxxyy_xxzzzz, tr_y_xxxyy_xyyyy, tr_y_xxxyy_xyyyyy, tr_y_xxxyy_xyyyyz, tr_y_xxxyy_xyyyz, tr_y_xxxyy_xyyyzz, tr_y_xxxyy_xyyzz, tr_y_xxxyy_xyyzzz, tr_y_xxxyy_xyzzz, tr_y_xxxyy_xyzzzz, tr_y_xxxyy_xzzzz, tr_y_xxxyy_xzzzzz, tr_y_xxxyy_yyyyy, tr_y_xxxyy_yyyyyy, tr_y_xxxyy_yyyyyz, tr_y_xxxyy_yyyyz, tr_y_xxxyy_yyyyzz, tr_y_xxxyy_yyyzz, tr_y_xxxyy_yyyzzz, tr_y_xxxyy_yyzzz, tr_y_xxxyy_yyzzzz, tr_y_xxxyy_yzzzz, tr_y_xxxyy_yzzzzz, tr_y_xxxyy_zzzzz, tr_y_xxxyy_zzzzzz, tr_y_xxyy_xxxxxx, tr_y_xxyy_xxxxxy, tr_y_xxyy_xxxxxz, tr_y_xxyy_xxxxyy, tr_y_xxyy_xxxxyz, tr_y_xxyy_xxxxzz, tr_y_xxyy_xxxyyy, tr_y_xxyy_xxxyyz, tr_y_xxyy_xxxyzz, tr_y_xxyy_xxxzzz, tr_y_xxyy_xxyyyy, tr_y_xxyy_xxyyyz, tr_y_xxyy_xxyyzz, tr_y_xxyy_xxyzzz, tr_y_xxyy_xxzzzz, tr_y_xxyy_xyyyyy, tr_y_xxyy_xyyyyz, tr_y_xxyy_xyyyzz, tr_y_xxyy_xyyzzz, tr_y_xxyy_xyzzzz, tr_y_xxyy_xzzzzz, tr_y_xxyy_yyyyyy, tr_y_xxyy_yyyyyz, tr_y_xxyy_yyyyzz, tr_y_xxyy_yyyzzz, tr_y_xxyy_yyzzzz, tr_y_xxyy_yzzzzz, tr_y_xxyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyy_xxxxxx[i] = 3.0 * tr_y_xxyy_xxxxxx[i] * fe_0 + 6.0 * tr_y_xxxyy_xxxxx[i] * fe_0 + tr_y_xxxyy_xxxxxx[i] * pa_x[i];

        tr_y_xxxxyy_xxxxxy[i] = 3.0 * tr_y_xxyy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxxyy_xxxxy[i] * fe_0 + tr_y_xxxyy_xxxxxy[i] * pa_x[i];

        tr_y_xxxxyy_xxxxxz[i] = 3.0 * tr_y_xxyy_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxxyy_xxxxz[i] * fe_0 + tr_y_xxxyy_xxxxxz[i] * pa_x[i];

        tr_y_xxxxyy_xxxxyy[i] = 3.0 * tr_y_xxyy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxxyy_xxxyy[i] * fe_0 + tr_y_xxxyy_xxxxyy[i] * pa_x[i];

        tr_y_xxxxyy_xxxxyz[i] = 3.0 * tr_y_xxyy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxxyy_xxxyz[i] * fe_0 + tr_y_xxxyy_xxxxyz[i] * pa_x[i];

        tr_y_xxxxyy_xxxxzz[i] = 3.0 * tr_y_xxyy_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxxyy_xxxzz[i] * fe_0 + tr_y_xxxyy_xxxxzz[i] * pa_x[i];

        tr_y_xxxxyy_xxxyyy[i] = 3.0 * tr_y_xxyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxxyy_xxyyy[i] * fe_0 + tr_y_xxxyy_xxxyyy[i] * pa_x[i];

        tr_y_xxxxyy_xxxyyz[i] = 3.0 * tr_y_xxyy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxxyy_xxyyz[i] * fe_0 + tr_y_xxxyy_xxxyyz[i] * pa_x[i];

        tr_y_xxxxyy_xxxyzz[i] = 3.0 * tr_y_xxyy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxxyy_xxyzz[i] * fe_0 + tr_y_xxxyy_xxxyzz[i] * pa_x[i];

        tr_y_xxxxyy_xxxzzz[i] = 3.0 * tr_y_xxyy_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxxyy_xxzzz[i] * fe_0 + tr_y_xxxyy_xxxzzz[i] * pa_x[i];

        tr_y_xxxxyy_xxyyyy[i] = 3.0 * tr_y_xxyy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxxyy_xyyyy[i] * fe_0 + tr_y_xxxyy_xxyyyy[i] * pa_x[i];

        tr_y_xxxxyy_xxyyyz[i] = 3.0 * tr_y_xxyy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxxyy_xyyyz[i] * fe_0 + tr_y_xxxyy_xxyyyz[i] * pa_x[i];

        tr_y_xxxxyy_xxyyzz[i] = 3.0 * tr_y_xxyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxxyy_xyyzz[i] * fe_0 + tr_y_xxxyy_xxyyzz[i] * pa_x[i];

        tr_y_xxxxyy_xxyzzz[i] = 3.0 * tr_y_xxyy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxxyy_xyzzz[i] * fe_0 + tr_y_xxxyy_xxyzzz[i] * pa_x[i];

        tr_y_xxxxyy_xxzzzz[i] = 3.0 * tr_y_xxyy_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxxyy_xzzzz[i] * fe_0 + tr_y_xxxyy_xxzzzz[i] * pa_x[i];

        tr_y_xxxxyy_xyyyyy[i] = 3.0 * tr_y_xxyy_xyyyyy[i] * fe_0 + tr_y_xxxyy_yyyyy[i] * fe_0 + tr_y_xxxyy_xyyyyy[i] * pa_x[i];

        tr_y_xxxxyy_xyyyyz[i] = 3.0 * tr_y_xxyy_xyyyyz[i] * fe_0 + tr_y_xxxyy_yyyyz[i] * fe_0 + tr_y_xxxyy_xyyyyz[i] * pa_x[i];

        tr_y_xxxxyy_xyyyzz[i] = 3.0 * tr_y_xxyy_xyyyzz[i] * fe_0 + tr_y_xxxyy_yyyzz[i] * fe_0 + tr_y_xxxyy_xyyyzz[i] * pa_x[i];

        tr_y_xxxxyy_xyyzzz[i] = 3.0 * tr_y_xxyy_xyyzzz[i] * fe_0 + tr_y_xxxyy_yyzzz[i] * fe_0 + tr_y_xxxyy_xyyzzz[i] * pa_x[i];

        tr_y_xxxxyy_xyzzzz[i] = 3.0 * tr_y_xxyy_xyzzzz[i] * fe_0 + tr_y_xxxyy_yzzzz[i] * fe_0 + tr_y_xxxyy_xyzzzz[i] * pa_x[i];

        tr_y_xxxxyy_xzzzzz[i] = 3.0 * tr_y_xxyy_xzzzzz[i] * fe_0 + tr_y_xxxyy_zzzzz[i] * fe_0 + tr_y_xxxyy_xzzzzz[i] * pa_x[i];

        tr_y_xxxxyy_yyyyyy[i] = 3.0 * tr_y_xxyy_yyyyyy[i] * fe_0 + tr_y_xxxyy_yyyyyy[i] * pa_x[i];

        tr_y_xxxxyy_yyyyyz[i] = 3.0 * tr_y_xxyy_yyyyyz[i] * fe_0 + tr_y_xxxyy_yyyyyz[i] * pa_x[i];

        tr_y_xxxxyy_yyyyzz[i] = 3.0 * tr_y_xxyy_yyyyzz[i] * fe_0 + tr_y_xxxyy_yyyyzz[i] * pa_x[i];

        tr_y_xxxxyy_yyyzzz[i] = 3.0 * tr_y_xxyy_yyyzzz[i] * fe_0 + tr_y_xxxyy_yyyzzz[i] * pa_x[i];

        tr_y_xxxxyy_yyzzzz[i] = 3.0 * tr_y_xxyy_yyzzzz[i] * fe_0 + tr_y_xxxyy_yyzzzz[i] * pa_x[i];

        tr_y_xxxxyy_yzzzzz[i] = 3.0 * tr_y_xxyy_yzzzzz[i] * fe_0 + tr_y_xxxyy_yzzzzz[i] * pa_x[i];

        tr_y_xxxxyy_zzzzzz[i] = 3.0 * tr_y_xxyy_zzzzzz[i] * fe_0 + tr_y_xxxyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 896-924 components of targeted buffer : II

    auto tr_y_xxxxyz_xxxxxx = pbuffer.data(idx_dip_ii + 896);

    auto tr_y_xxxxyz_xxxxxy = pbuffer.data(idx_dip_ii + 897);

    auto tr_y_xxxxyz_xxxxxz = pbuffer.data(idx_dip_ii + 898);

    auto tr_y_xxxxyz_xxxxyy = pbuffer.data(idx_dip_ii + 899);

    auto tr_y_xxxxyz_xxxxyz = pbuffer.data(idx_dip_ii + 900);

    auto tr_y_xxxxyz_xxxxzz = pbuffer.data(idx_dip_ii + 901);

    auto tr_y_xxxxyz_xxxyyy = pbuffer.data(idx_dip_ii + 902);

    auto tr_y_xxxxyz_xxxyyz = pbuffer.data(idx_dip_ii + 903);

    auto tr_y_xxxxyz_xxxyzz = pbuffer.data(idx_dip_ii + 904);

    auto tr_y_xxxxyz_xxxzzz = pbuffer.data(idx_dip_ii + 905);

    auto tr_y_xxxxyz_xxyyyy = pbuffer.data(idx_dip_ii + 906);

    auto tr_y_xxxxyz_xxyyyz = pbuffer.data(idx_dip_ii + 907);

    auto tr_y_xxxxyz_xxyyzz = pbuffer.data(idx_dip_ii + 908);

    auto tr_y_xxxxyz_xxyzzz = pbuffer.data(idx_dip_ii + 909);

    auto tr_y_xxxxyz_xxzzzz = pbuffer.data(idx_dip_ii + 910);

    auto tr_y_xxxxyz_xyyyyy = pbuffer.data(idx_dip_ii + 911);

    auto tr_y_xxxxyz_xyyyyz = pbuffer.data(idx_dip_ii + 912);

    auto tr_y_xxxxyz_xyyyzz = pbuffer.data(idx_dip_ii + 913);

    auto tr_y_xxxxyz_xyyzzz = pbuffer.data(idx_dip_ii + 914);

    auto tr_y_xxxxyz_xyzzzz = pbuffer.data(idx_dip_ii + 915);

    auto tr_y_xxxxyz_xzzzzz = pbuffer.data(idx_dip_ii + 916);

    auto tr_y_xxxxyz_yyyyyy = pbuffer.data(idx_dip_ii + 917);

    auto tr_y_xxxxyz_yyyyyz = pbuffer.data(idx_dip_ii + 918);

    auto tr_y_xxxxyz_yyyyzz = pbuffer.data(idx_dip_ii + 919);

    auto tr_y_xxxxyz_yyyzzz = pbuffer.data(idx_dip_ii + 920);

    auto tr_y_xxxxyz_yyzzzz = pbuffer.data(idx_dip_ii + 921);

    auto tr_y_xxxxyz_yzzzzz = pbuffer.data(idx_dip_ii + 922);

    auto tr_y_xxxxyz_zzzzzz = pbuffer.data(idx_dip_ii + 923);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxxy_xxxxxx, tr_y_xxxxy_xxxxxy, tr_y_xxxxy_xxxxy, tr_y_xxxxy_xxxxyy, tr_y_xxxxy_xxxxyz, tr_y_xxxxy_xxxyy, tr_y_xxxxy_xxxyyy, tr_y_xxxxy_xxxyyz, tr_y_xxxxy_xxxyz, tr_y_xxxxy_xxxyzz, tr_y_xxxxy_xxyyy, tr_y_xxxxy_xxyyyy, tr_y_xxxxy_xxyyyz, tr_y_xxxxy_xxyyz, tr_y_xxxxy_xxyyzz, tr_y_xxxxy_xxyzz, tr_y_xxxxy_xxyzzz, tr_y_xxxxy_xyyyy, tr_y_xxxxy_xyyyyy, tr_y_xxxxy_xyyyyz, tr_y_xxxxy_xyyyz, tr_y_xxxxy_xyyyzz, tr_y_xxxxy_xyyzz, tr_y_xxxxy_xyyzzz, tr_y_xxxxy_xyzzz, tr_y_xxxxy_xyzzzz, tr_y_xxxxy_yyyyyy, tr_y_xxxxyz_xxxxxx, tr_y_xxxxyz_xxxxxy, tr_y_xxxxyz_xxxxxz, tr_y_xxxxyz_xxxxyy, tr_y_xxxxyz_xxxxyz, tr_y_xxxxyz_xxxxzz, tr_y_xxxxyz_xxxyyy, tr_y_xxxxyz_xxxyyz, tr_y_xxxxyz_xxxyzz, tr_y_xxxxyz_xxxzzz, tr_y_xxxxyz_xxyyyy, tr_y_xxxxyz_xxyyyz, tr_y_xxxxyz_xxyyzz, tr_y_xxxxyz_xxyzzz, tr_y_xxxxyz_xxzzzz, tr_y_xxxxyz_xyyyyy, tr_y_xxxxyz_xyyyyz, tr_y_xxxxyz_xyyyzz, tr_y_xxxxyz_xyyzzz, tr_y_xxxxyz_xyzzzz, tr_y_xxxxyz_xzzzzz, tr_y_xxxxyz_yyyyyy, tr_y_xxxxyz_yyyyyz, tr_y_xxxxyz_yyyyzz, tr_y_xxxxyz_yyyzzz, tr_y_xxxxyz_yyzzzz, tr_y_xxxxyz_yzzzzz, tr_y_xxxxyz_zzzzzz, tr_y_xxxxz_xxxxxz, tr_y_xxxxz_xxxxzz, tr_y_xxxxz_xxxzzz, tr_y_xxxxz_xxzzzz, tr_y_xxxxz_xzzzzz, tr_y_xxxyz_yyyyyz, tr_y_xxxyz_yyyyzz, tr_y_xxxyz_yyyzzz, tr_y_xxxyz_yyzzzz, tr_y_xxxyz_yzzzzz, tr_y_xxxyz_zzzzzz, tr_y_xxyz_yyyyyz, tr_y_xxyz_yyyyzz, tr_y_xxyz_yyyzzz, tr_y_xxyz_yyzzzz, tr_y_xxyz_yzzzzz, tr_y_xxyz_zzzzzz, ts_xxxxz_xxxxxz, ts_xxxxz_xxxxzz, ts_xxxxz_xxxzzz, ts_xxxxz_xxzzzz, ts_xxxxz_xzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyz_xxxxxx[i] = tr_y_xxxxy_xxxxxx[i] * pa_z[i];

        tr_y_xxxxyz_xxxxxy[i] = tr_y_xxxxy_xxxxxy[i] * pa_z[i];

        tr_y_xxxxyz_xxxxxz[i] = ts_xxxxz_xxxxxz[i] * fe_0 + tr_y_xxxxz_xxxxxz[i] * pa_y[i];

        tr_y_xxxxyz_xxxxyy[i] = tr_y_xxxxy_xxxxyy[i] * pa_z[i];

        tr_y_xxxxyz_xxxxyz[i] = tr_y_xxxxy_xxxxy[i] * fe_0 + tr_y_xxxxy_xxxxyz[i] * pa_z[i];

        tr_y_xxxxyz_xxxxzz[i] = ts_xxxxz_xxxxzz[i] * fe_0 + tr_y_xxxxz_xxxxzz[i] * pa_y[i];

        tr_y_xxxxyz_xxxyyy[i] = tr_y_xxxxy_xxxyyy[i] * pa_z[i];

        tr_y_xxxxyz_xxxyyz[i] = tr_y_xxxxy_xxxyy[i] * fe_0 + tr_y_xxxxy_xxxyyz[i] * pa_z[i];

        tr_y_xxxxyz_xxxyzz[i] = 2.0 * tr_y_xxxxy_xxxyz[i] * fe_0 + tr_y_xxxxy_xxxyzz[i] * pa_z[i];

        tr_y_xxxxyz_xxxzzz[i] = ts_xxxxz_xxxzzz[i] * fe_0 + tr_y_xxxxz_xxxzzz[i] * pa_y[i];

        tr_y_xxxxyz_xxyyyy[i] = tr_y_xxxxy_xxyyyy[i] * pa_z[i];

        tr_y_xxxxyz_xxyyyz[i] = tr_y_xxxxy_xxyyy[i] * fe_0 + tr_y_xxxxy_xxyyyz[i] * pa_z[i];

        tr_y_xxxxyz_xxyyzz[i] = 2.0 * tr_y_xxxxy_xxyyz[i] * fe_0 + tr_y_xxxxy_xxyyzz[i] * pa_z[i];

        tr_y_xxxxyz_xxyzzz[i] = 3.0 * tr_y_xxxxy_xxyzz[i] * fe_0 + tr_y_xxxxy_xxyzzz[i] * pa_z[i];

        tr_y_xxxxyz_xxzzzz[i] = ts_xxxxz_xxzzzz[i] * fe_0 + tr_y_xxxxz_xxzzzz[i] * pa_y[i];

        tr_y_xxxxyz_xyyyyy[i] = tr_y_xxxxy_xyyyyy[i] * pa_z[i];

        tr_y_xxxxyz_xyyyyz[i] = tr_y_xxxxy_xyyyy[i] * fe_0 + tr_y_xxxxy_xyyyyz[i] * pa_z[i];

        tr_y_xxxxyz_xyyyzz[i] = 2.0 * tr_y_xxxxy_xyyyz[i] * fe_0 + tr_y_xxxxy_xyyyzz[i] * pa_z[i];

        tr_y_xxxxyz_xyyzzz[i] = 3.0 * tr_y_xxxxy_xyyzz[i] * fe_0 + tr_y_xxxxy_xyyzzz[i] * pa_z[i];

        tr_y_xxxxyz_xyzzzz[i] = 4.0 * tr_y_xxxxy_xyzzz[i] * fe_0 + tr_y_xxxxy_xyzzzz[i] * pa_z[i];

        tr_y_xxxxyz_xzzzzz[i] = ts_xxxxz_xzzzzz[i] * fe_0 + tr_y_xxxxz_xzzzzz[i] * pa_y[i];

        tr_y_xxxxyz_yyyyyy[i] = tr_y_xxxxy_yyyyyy[i] * pa_z[i];

        tr_y_xxxxyz_yyyyyz[i] = 3.0 * tr_y_xxyz_yyyyyz[i] * fe_0 + tr_y_xxxyz_yyyyyz[i] * pa_x[i];

        tr_y_xxxxyz_yyyyzz[i] = 3.0 * tr_y_xxyz_yyyyzz[i] * fe_0 + tr_y_xxxyz_yyyyzz[i] * pa_x[i];

        tr_y_xxxxyz_yyyzzz[i] = 3.0 * tr_y_xxyz_yyyzzz[i] * fe_0 + tr_y_xxxyz_yyyzzz[i] * pa_x[i];

        tr_y_xxxxyz_yyzzzz[i] = 3.0 * tr_y_xxyz_yyzzzz[i] * fe_0 + tr_y_xxxyz_yyzzzz[i] * pa_x[i];

        tr_y_xxxxyz_yzzzzz[i] = 3.0 * tr_y_xxyz_yzzzzz[i] * fe_0 + tr_y_xxxyz_yzzzzz[i] * pa_x[i];

        tr_y_xxxxyz_zzzzzz[i] = 3.0 * tr_y_xxyz_zzzzzz[i] * fe_0 + tr_y_xxxyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 924-952 components of targeted buffer : II

    auto tr_y_xxxxzz_xxxxxx = pbuffer.data(idx_dip_ii + 924);

    auto tr_y_xxxxzz_xxxxxy = pbuffer.data(idx_dip_ii + 925);

    auto tr_y_xxxxzz_xxxxxz = pbuffer.data(idx_dip_ii + 926);

    auto tr_y_xxxxzz_xxxxyy = pbuffer.data(idx_dip_ii + 927);

    auto tr_y_xxxxzz_xxxxyz = pbuffer.data(idx_dip_ii + 928);

    auto tr_y_xxxxzz_xxxxzz = pbuffer.data(idx_dip_ii + 929);

    auto tr_y_xxxxzz_xxxyyy = pbuffer.data(idx_dip_ii + 930);

    auto tr_y_xxxxzz_xxxyyz = pbuffer.data(idx_dip_ii + 931);

    auto tr_y_xxxxzz_xxxyzz = pbuffer.data(idx_dip_ii + 932);

    auto tr_y_xxxxzz_xxxzzz = pbuffer.data(idx_dip_ii + 933);

    auto tr_y_xxxxzz_xxyyyy = pbuffer.data(idx_dip_ii + 934);

    auto tr_y_xxxxzz_xxyyyz = pbuffer.data(idx_dip_ii + 935);

    auto tr_y_xxxxzz_xxyyzz = pbuffer.data(idx_dip_ii + 936);

    auto tr_y_xxxxzz_xxyzzz = pbuffer.data(idx_dip_ii + 937);

    auto tr_y_xxxxzz_xxzzzz = pbuffer.data(idx_dip_ii + 938);

    auto tr_y_xxxxzz_xyyyyy = pbuffer.data(idx_dip_ii + 939);

    auto tr_y_xxxxzz_xyyyyz = pbuffer.data(idx_dip_ii + 940);

    auto tr_y_xxxxzz_xyyyzz = pbuffer.data(idx_dip_ii + 941);

    auto tr_y_xxxxzz_xyyzzz = pbuffer.data(idx_dip_ii + 942);

    auto tr_y_xxxxzz_xyzzzz = pbuffer.data(idx_dip_ii + 943);

    auto tr_y_xxxxzz_xzzzzz = pbuffer.data(idx_dip_ii + 944);

    auto tr_y_xxxxzz_yyyyyy = pbuffer.data(idx_dip_ii + 945);

    auto tr_y_xxxxzz_yyyyyz = pbuffer.data(idx_dip_ii + 946);

    auto tr_y_xxxxzz_yyyyzz = pbuffer.data(idx_dip_ii + 947);

    auto tr_y_xxxxzz_yyyzzz = pbuffer.data(idx_dip_ii + 948);

    auto tr_y_xxxxzz_yyzzzz = pbuffer.data(idx_dip_ii + 949);

    auto tr_y_xxxxzz_yzzzzz = pbuffer.data(idx_dip_ii + 950);

    auto tr_y_xxxxzz_zzzzzz = pbuffer.data(idx_dip_ii + 951);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxx_xxxxxx, tr_y_xxxx_xxxxxy, tr_y_xxxx_xxxxyy, tr_y_xxxx_xxxyyy, tr_y_xxxx_xxyyyy, tr_y_xxxx_xyyyyy, tr_y_xxxxz_xxxxxx, tr_y_xxxxz_xxxxxy, tr_y_xxxxz_xxxxyy, tr_y_xxxxz_xxxyyy, tr_y_xxxxz_xxyyyy, tr_y_xxxxz_xyyyyy, tr_y_xxxxzz_xxxxxx, tr_y_xxxxzz_xxxxxy, tr_y_xxxxzz_xxxxxz, tr_y_xxxxzz_xxxxyy, tr_y_xxxxzz_xxxxyz, tr_y_xxxxzz_xxxxzz, tr_y_xxxxzz_xxxyyy, tr_y_xxxxzz_xxxyyz, tr_y_xxxxzz_xxxyzz, tr_y_xxxxzz_xxxzzz, tr_y_xxxxzz_xxyyyy, tr_y_xxxxzz_xxyyyz, tr_y_xxxxzz_xxyyzz, tr_y_xxxxzz_xxyzzz, tr_y_xxxxzz_xxzzzz, tr_y_xxxxzz_xyyyyy, tr_y_xxxxzz_xyyyyz, tr_y_xxxxzz_xyyyzz, tr_y_xxxxzz_xyyzzz, tr_y_xxxxzz_xyzzzz, tr_y_xxxxzz_xzzzzz, tr_y_xxxxzz_yyyyyy, tr_y_xxxxzz_yyyyyz, tr_y_xxxxzz_yyyyzz, tr_y_xxxxzz_yyyzzz, tr_y_xxxxzz_yyzzzz, tr_y_xxxxzz_yzzzzz, tr_y_xxxxzz_zzzzzz, tr_y_xxxzz_xxxxxz, tr_y_xxxzz_xxxxyz, tr_y_xxxzz_xxxxz, tr_y_xxxzz_xxxxzz, tr_y_xxxzz_xxxyyz, tr_y_xxxzz_xxxyz, tr_y_xxxzz_xxxyzz, tr_y_xxxzz_xxxzz, tr_y_xxxzz_xxxzzz, tr_y_xxxzz_xxyyyz, tr_y_xxxzz_xxyyz, tr_y_xxxzz_xxyyzz, tr_y_xxxzz_xxyzz, tr_y_xxxzz_xxyzzz, tr_y_xxxzz_xxzzz, tr_y_xxxzz_xxzzzz, tr_y_xxxzz_xyyyyz, tr_y_xxxzz_xyyyz, tr_y_xxxzz_xyyyzz, tr_y_xxxzz_xyyzz, tr_y_xxxzz_xyyzzz, tr_y_xxxzz_xyzzz, tr_y_xxxzz_xyzzzz, tr_y_xxxzz_xzzzz, tr_y_xxxzz_xzzzzz, tr_y_xxxzz_yyyyyy, tr_y_xxxzz_yyyyyz, tr_y_xxxzz_yyyyz, tr_y_xxxzz_yyyyzz, tr_y_xxxzz_yyyzz, tr_y_xxxzz_yyyzzz, tr_y_xxxzz_yyzzz, tr_y_xxxzz_yyzzzz, tr_y_xxxzz_yzzzz, tr_y_xxxzz_yzzzzz, tr_y_xxxzz_zzzzz, tr_y_xxxzz_zzzzzz, tr_y_xxzz_xxxxxz, tr_y_xxzz_xxxxyz, tr_y_xxzz_xxxxzz, tr_y_xxzz_xxxyyz, tr_y_xxzz_xxxyzz, tr_y_xxzz_xxxzzz, tr_y_xxzz_xxyyyz, tr_y_xxzz_xxyyzz, tr_y_xxzz_xxyzzz, tr_y_xxzz_xxzzzz, tr_y_xxzz_xyyyyz, tr_y_xxzz_xyyyzz, tr_y_xxzz_xyyzzz, tr_y_xxzz_xyzzzz, tr_y_xxzz_xzzzzz, tr_y_xxzz_yyyyyy, tr_y_xxzz_yyyyyz, tr_y_xxzz_yyyyzz, tr_y_xxzz_yyyzzz, tr_y_xxzz_yyzzzz, tr_y_xxzz_yzzzzz, tr_y_xxzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxzz_xxxxxx[i] = tr_y_xxxx_xxxxxx[i] * fe_0 + tr_y_xxxxz_xxxxxx[i] * pa_z[i];

        tr_y_xxxxzz_xxxxxy[i] = tr_y_xxxx_xxxxxy[i] * fe_0 + tr_y_xxxxz_xxxxxy[i] * pa_z[i];

        tr_y_xxxxzz_xxxxxz[i] = 3.0 * tr_y_xxzz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxxzz_xxxxz[i] * fe_0 + tr_y_xxxzz_xxxxxz[i] * pa_x[i];

        tr_y_xxxxzz_xxxxyy[i] = tr_y_xxxx_xxxxyy[i] * fe_0 + tr_y_xxxxz_xxxxyy[i] * pa_z[i];

        tr_y_xxxxzz_xxxxyz[i] = 3.0 * tr_y_xxzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxxzz_xxxyz[i] * fe_0 + tr_y_xxxzz_xxxxyz[i] * pa_x[i];

        tr_y_xxxxzz_xxxxzz[i] = 3.0 * tr_y_xxzz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxxzz_xxxzz[i] * fe_0 + tr_y_xxxzz_xxxxzz[i] * pa_x[i];

        tr_y_xxxxzz_xxxyyy[i] = tr_y_xxxx_xxxyyy[i] * fe_0 + tr_y_xxxxz_xxxyyy[i] * pa_z[i];

        tr_y_xxxxzz_xxxyyz[i] = 3.0 * tr_y_xxzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxxzz_xxyyz[i] * fe_0 + tr_y_xxxzz_xxxyyz[i] * pa_x[i];

        tr_y_xxxxzz_xxxyzz[i] = 3.0 * tr_y_xxzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxxzz_xxyzz[i] * fe_0 + tr_y_xxxzz_xxxyzz[i] * pa_x[i];

        tr_y_xxxxzz_xxxzzz[i] = 3.0 * tr_y_xxzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxxzz_xxzzz[i] * fe_0 + tr_y_xxxzz_xxxzzz[i] * pa_x[i];

        tr_y_xxxxzz_xxyyyy[i] = tr_y_xxxx_xxyyyy[i] * fe_0 + tr_y_xxxxz_xxyyyy[i] * pa_z[i];

        tr_y_xxxxzz_xxyyyz[i] = 3.0 * tr_y_xxzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxxzz_xyyyz[i] * fe_0 + tr_y_xxxzz_xxyyyz[i] * pa_x[i];

        tr_y_xxxxzz_xxyyzz[i] = 3.0 * tr_y_xxzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxxzz_xyyzz[i] * fe_0 + tr_y_xxxzz_xxyyzz[i] * pa_x[i];

        tr_y_xxxxzz_xxyzzz[i] = 3.0 * tr_y_xxzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxxzz_xyzzz[i] * fe_0 + tr_y_xxxzz_xxyzzz[i] * pa_x[i];

        tr_y_xxxxzz_xxzzzz[i] = 3.0 * tr_y_xxzz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxxzz_xzzzz[i] * fe_0 + tr_y_xxxzz_xxzzzz[i] * pa_x[i];

        tr_y_xxxxzz_xyyyyy[i] = tr_y_xxxx_xyyyyy[i] * fe_0 + tr_y_xxxxz_xyyyyy[i] * pa_z[i];

        tr_y_xxxxzz_xyyyyz[i] = 3.0 * tr_y_xxzz_xyyyyz[i] * fe_0 + tr_y_xxxzz_yyyyz[i] * fe_0 + tr_y_xxxzz_xyyyyz[i] * pa_x[i];

        tr_y_xxxxzz_xyyyzz[i] = 3.0 * tr_y_xxzz_xyyyzz[i] * fe_0 + tr_y_xxxzz_yyyzz[i] * fe_0 + tr_y_xxxzz_xyyyzz[i] * pa_x[i];

        tr_y_xxxxzz_xyyzzz[i] = 3.0 * tr_y_xxzz_xyyzzz[i] * fe_0 + tr_y_xxxzz_yyzzz[i] * fe_0 + tr_y_xxxzz_xyyzzz[i] * pa_x[i];

        tr_y_xxxxzz_xyzzzz[i] = 3.0 * tr_y_xxzz_xyzzzz[i] * fe_0 + tr_y_xxxzz_yzzzz[i] * fe_0 + tr_y_xxxzz_xyzzzz[i] * pa_x[i];

        tr_y_xxxxzz_xzzzzz[i] = 3.0 * tr_y_xxzz_xzzzzz[i] * fe_0 + tr_y_xxxzz_zzzzz[i] * fe_0 + tr_y_xxxzz_xzzzzz[i] * pa_x[i];

        tr_y_xxxxzz_yyyyyy[i] = 3.0 * tr_y_xxzz_yyyyyy[i] * fe_0 + tr_y_xxxzz_yyyyyy[i] * pa_x[i];

        tr_y_xxxxzz_yyyyyz[i] = 3.0 * tr_y_xxzz_yyyyyz[i] * fe_0 + tr_y_xxxzz_yyyyyz[i] * pa_x[i];

        tr_y_xxxxzz_yyyyzz[i] = 3.0 * tr_y_xxzz_yyyyzz[i] * fe_0 + tr_y_xxxzz_yyyyzz[i] * pa_x[i];

        tr_y_xxxxzz_yyyzzz[i] = 3.0 * tr_y_xxzz_yyyzzz[i] * fe_0 + tr_y_xxxzz_yyyzzz[i] * pa_x[i];

        tr_y_xxxxzz_yyzzzz[i] = 3.0 * tr_y_xxzz_yyzzzz[i] * fe_0 + tr_y_xxxzz_yyzzzz[i] * pa_x[i];

        tr_y_xxxxzz_yzzzzz[i] = 3.0 * tr_y_xxzz_yzzzzz[i] * fe_0 + tr_y_xxxzz_yzzzzz[i] * pa_x[i];

        tr_y_xxxxzz_zzzzzz[i] = 3.0 * tr_y_xxzz_zzzzzz[i] * fe_0 + tr_y_xxxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 952-980 components of targeted buffer : II

    auto tr_y_xxxyyy_xxxxxx = pbuffer.data(idx_dip_ii + 952);

    auto tr_y_xxxyyy_xxxxxy = pbuffer.data(idx_dip_ii + 953);

    auto tr_y_xxxyyy_xxxxxz = pbuffer.data(idx_dip_ii + 954);

    auto tr_y_xxxyyy_xxxxyy = pbuffer.data(idx_dip_ii + 955);

    auto tr_y_xxxyyy_xxxxyz = pbuffer.data(idx_dip_ii + 956);

    auto tr_y_xxxyyy_xxxxzz = pbuffer.data(idx_dip_ii + 957);

    auto tr_y_xxxyyy_xxxyyy = pbuffer.data(idx_dip_ii + 958);

    auto tr_y_xxxyyy_xxxyyz = pbuffer.data(idx_dip_ii + 959);

    auto tr_y_xxxyyy_xxxyzz = pbuffer.data(idx_dip_ii + 960);

    auto tr_y_xxxyyy_xxxzzz = pbuffer.data(idx_dip_ii + 961);

    auto tr_y_xxxyyy_xxyyyy = pbuffer.data(idx_dip_ii + 962);

    auto tr_y_xxxyyy_xxyyyz = pbuffer.data(idx_dip_ii + 963);

    auto tr_y_xxxyyy_xxyyzz = pbuffer.data(idx_dip_ii + 964);

    auto tr_y_xxxyyy_xxyzzz = pbuffer.data(idx_dip_ii + 965);

    auto tr_y_xxxyyy_xxzzzz = pbuffer.data(idx_dip_ii + 966);

    auto tr_y_xxxyyy_xyyyyy = pbuffer.data(idx_dip_ii + 967);

    auto tr_y_xxxyyy_xyyyyz = pbuffer.data(idx_dip_ii + 968);

    auto tr_y_xxxyyy_xyyyzz = pbuffer.data(idx_dip_ii + 969);

    auto tr_y_xxxyyy_xyyzzz = pbuffer.data(idx_dip_ii + 970);

    auto tr_y_xxxyyy_xyzzzz = pbuffer.data(idx_dip_ii + 971);

    auto tr_y_xxxyyy_xzzzzz = pbuffer.data(idx_dip_ii + 972);

    auto tr_y_xxxyyy_yyyyyy = pbuffer.data(idx_dip_ii + 973);

    auto tr_y_xxxyyy_yyyyyz = pbuffer.data(idx_dip_ii + 974);

    auto tr_y_xxxyyy_yyyyzz = pbuffer.data(idx_dip_ii + 975);

    auto tr_y_xxxyyy_yyyzzz = pbuffer.data(idx_dip_ii + 976);

    auto tr_y_xxxyyy_yyzzzz = pbuffer.data(idx_dip_ii + 977);

    auto tr_y_xxxyyy_yzzzzz = pbuffer.data(idx_dip_ii + 978);

    auto tr_y_xxxyyy_zzzzzz = pbuffer.data(idx_dip_ii + 979);

    #pragma omp simd aligned(pa_x, tr_y_xxxyyy_xxxxxx, tr_y_xxxyyy_xxxxxy, tr_y_xxxyyy_xxxxxz, tr_y_xxxyyy_xxxxyy, tr_y_xxxyyy_xxxxyz, tr_y_xxxyyy_xxxxzz, tr_y_xxxyyy_xxxyyy, tr_y_xxxyyy_xxxyyz, tr_y_xxxyyy_xxxyzz, tr_y_xxxyyy_xxxzzz, tr_y_xxxyyy_xxyyyy, tr_y_xxxyyy_xxyyyz, tr_y_xxxyyy_xxyyzz, tr_y_xxxyyy_xxyzzz, tr_y_xxxyyy_xxzzzz, tr_y_xxxyyy_xyyyyy, tr_y_xxxyyy_xyyyyz, tr_y_xxxyyy_xyyyzz, tr_y_xxxyyy_xyyzzz, tr_y_xxxyyy_xyzzzz, tr_y_xxxyyy_xzzzzz, tr_y_xxxyyy_yyyyyy, tr_y_xxxyyy_yyyyyz, tr_y_xxxyyy_yyyyzz, tr_y_xxxyyy_yyyzzz, tr_y_xxxyyy_yyzzzz, tr_y_xxxyyy_yzzzzz, tr_y_xxxyyy_zzzzzz, tr_y_xxyyy_xxxxx, tr_y_xxyyy_xxxxxx, tr_y_xxyyy_xxxxxy, tr_y_xxyyy_xxxxxz, tr_y_xxyyy_xxxxy, tr_y_xxyyy_xxxxyy, tr_y_xxyyy_xxxxyz, tr_y_xxyyy_xxxxz, tr_y_xxyyy_xxxxzz, tr_y_xxyyy_xxxyy, tr_y_xxyyy_xxxyyy, tr_y_xxyyy_xxxyyz, tr_y_xxyyy_xxxyz, tr_y_xxyyy_xxxyzz, tr_y_xxyyy_xxxzz, tr_y_xxyyy_xxxzzz, tr_y_xxyyy_xxyyy, tr_y_xxyyy_xxyyyy, tr_y_xxyyy_xxyyyz, tr_y_xxyyy_xxyyz, tr_y_xxyyy_xxyyzz, tr_y_xxyyy_xxyzz, tr_y_xxyyy_xxyzzz, tr_y_xxyyy_xxzzz, tr_y_xxyyy_xxzzzz, tr_y_xxyyy_xyyyy, tr_y_xxyyy_xyyyyy, tr_y_xxyyy_xyyyyz, tr_y_xxyyy_xyyyz, tr_y_xxyyy_xyyyzz, tr_y_xxyyy_xyyzz, tr_y_xxyyy_xyyzzz, tr_y_xxyyy_xyzzz, tr_y_xxyyy_xyzzzz, tr_y_xxyyy_xzzzz, tr_y_xxyyy_xzzzzz, tr_y_xxyyy_yyyyy, tr_y_xxyyy_yyyyyy, tr_y_xxyyy_yyyyyz, tr_y_xxyyy_yyyyz, tr_y_xxyyy_yyyyzz, tr_y_xxyyy_yyyzz, tr_y_xxyyy_yyyzzz, tr_y_xxyyy_yyzzz, tr_y_xxyyy_yyzzzz, tr_y_xxyyy_yzzzz, tr_y_xxyyy_yzzzzz, tr_y_xxyyy_zzzzz, tr_y_xxyyy_zzzzzz, tr_y_xyyy_xxxxxx, tr_y_xyyy_xxxxxy, tr_y_xyyy_xxxxxz, tr_y_xyyy_xxxxyy, tr_y_xyyy_xxxxyz, tr_y_xyyy_xxxxzz, tr_y_xyyy_xxxyyy, tr_y_xyyy_xxxyyz, tr_y_xyyy_xxxyzz, tr_y_xyyy_xxxzzz, tr_y_xyyy_xxyyyy, tr_y_xyyy_xxyyyz, tr_y_xyyy_xxyyzz, tr_y_xyyy_xxyzzz, tr_y_xyyy_xxzzzz, tr_y_xyyy_xyyyyy, tr_y_xyyy_xyyyyz, tr_y_xyyy_xyyyzz, tr_y_xyyy_xyyzzz, tr_y_xyyy_xyzzzz, tr_y_xyyy_xzzzzz, tr_y_xyyy_yyyyyy, tr_y_xyyy_yyyyyz, tr_y_xyyy_yyyyzz, tr_y_xyyy_yyyzzz, tr_y_xyyy_yyzzzz, tr_y_xyyy_yzzzzz, tr_y_xyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyy_xxxxxx[i] = 2.0 * tr_y_xyyy_xxxxxx[i] * fe_0 + 6.0 * tr_y_xxyyy_xxxxx[i] * fe_0 + tr_y_xxyyy_xxxxxx[i] * pa_x[i];

        tr_y_xxxyyy_xxxxxy[i] = 2.0 * tr_y_xyyy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxyyy_xxxxy[i] * fe_0 + tr_y_xxyyy_xxxxxy[i] * pa_x[i];

        tr_y_xxxyyy_xxxxxz[i] = 2.0 * tr_y_xyyy_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxyyy_xxxxz[i] * fe_0 + tr_y_xxyyy_xxxxxz[i] * pa_x[i];

        tr_y_xxxyyy_xxxxyy[i] = 2.0 * tr_y_xyyy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxyyy_xxxyy[i] * fe_0 + tr_y_xxyyy_xxxxyy[i] * pa_x[i];

        tr_y_xxxyyy_xxxxyz[i] = 2.0 * tr_y_xyyy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxyyy_xxxyz[i] * fe_0 + tr_y_xxyyy_xxxxyz[i] * pa_x[i];

        tr_y_xxxyyy_xxxxzz[i] = 2.0 * tr_y_xyyy_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxyyy_xxxzz[i] * fe_0 + tr_y_xxyyy_xxxxzz[i] * pa_x[i];

        tr_y_xxxyyy_xxxyyy[i] = 2.0 * tr_y_xyyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxyyy_xxyyy[i] * fe_0 + tr_y_xxyyy_xxxyyy[i] * pa_x[i];

        tr_y_xxxyyy_xxxyyz[i] = 2.0 * tr_y_xyyy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxyyy_xxyyz[i] * fe_0 + tr_y_xxyyy_xxxyyz[i] * pa_x[i];

        tr_y_xxxyyy_xxxyzz[i] = 2.0 * tr_y_xyyy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxyyy_xxyzz[i] * fe_0 + tr_y_xxyyy_xxxyzz[i] * pa_x[i];

        tr_y_xxxyyy_xxxzzz[i] = 2.0 * tr_y_xyyy_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxyyy_xxzzz[i] * fe_0 + tr_y_xxyyy_xxxzzz[i] * pa_x[i];

        tr_y_xxxyyy_xxyyyy[i] = 2.0 * tr_y_xyyy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxyyy_xyyyy[i] * fe_0 + tr_y_xxyyy_xxyyyy[i] * pa_x[i];

        tr_y_xxxyyy_xxyyyz[i] = 2.0 * tr_y_xyyy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxyyy_xyyyz[i] * fe_0 + tr_y_xxyyy_xxyyyz[i] * pa_x[i];

        tr_y_xxxyyy_xxyyzz[i] = 2.0 * tr_y_xyyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxyyy_xyyzz[i] * fe_0 + tr_y_xxyyy_xxyyzz[i] * pa_x[i];

        tr_y_xxxyyy_xxyzzz[i] = 2.0 * tr_y_xyyy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxyyy_xyzzz[i] * fe_0 + tr_y_xxyyy_xxyzzz[i] * pa_x[i];

        tr_y_xxxyyy_xxzzzz[i] = 2.0 * tr_y_xyyy_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxyyy_xzzzz[i] * fe_0 + tr_y_xxyyy_xxzzzz[i] * pa_x[i];

        tr_y_xxxyyy_xyyyyy[i] = 2.0 * tr_y_xyyy_xyyyyy[i] * fe_0 + tr_y_xxyyy_yyyyy[i] * fe_0 + tr_y_xxyyy_xyyyyy[i] * pa_x[i];

        tr_y_xxxyyy_xyyyyz[i] = 2.0 * tr_y_xyyy_xyyyyz[i] * fe_0 + tr_y_xxyyy_yyyyz[i] * fe_0 + tr_y_xxyyy_xyyyyz[i] * pa_x[i];

        tr_y_xxxyyy_xyyyzz[i] = 2.0 * tr_y_xyyy_xyyyzz[i] * fe_0 + tr_y_xxyyy_yyyzz[i] * fe_0 + tr_y_xxyyy_xyyyzz[i] * pa_x[i];

        tr_y_xxxyyy_xyyzzz[i] = 2.0 * tr_y_xyyy_xyyzzz[i] * fe_0 + tr_y_xxyyy_yyzzz[i] * fe_0 + tr_y_xxyyy_xyyzzz[i] * pa_x[i];

        tr_y_xxxyyy_xyzzzz[i] = 2.0 * tr_y_xyyy_xyzzzz[i] * fe_0 + tr_y_xxyyy_yzzzz[i] * fe_0 + tr_y_xxyyy_xyzzzz[i] * pa_x[i];

        tr_y_xxxyyy_xzzzzz[i] = 2.0 * tr_y_xyyy_xzzzzz[i] * fe_0 + tr_y_xxyyy_zzzzz[i] * fe_0 + tr_y_xxyyy_xzzzzz[i] * pa_x[i];

        tr_y_xxxyyy_yyyyyy[i] = 2.0 * tr_y_xyyy_yyyyyy[i] * fe_0 + tr_y_xxyyy_yyyyyy[i] * pa_x[i];

        tr_y_xxxyyy_yyyyyz[i] = 2.0 * tr_y_xyyy_yyyyyz[i] * fe_0 + tr_y_xxyyy_yyyyyz[i] * pa_x[i];

        tr_y_xxxyyy_yyyyzz[i] = 2.0 * tr_y_xyyy_yyyyzz[i] * fe_0 + tr_y_xxyyy_yyyyzz[i] * pa_x[i];

        tr_y_xxxyyy_yyyzzz[i] = 2.0 * tr_y_xyyy_yyyzzz[i] * fe_0 + tr_y_xxyyy_yyyzzz[i] * pa_x[i];

        tr_y_xxxyyy_yyzzzz[i] = 2.0 * tr_y_xyyy_yyzzzz[i] * fe_0 + tr_y_xxyyy_yyzzzz[i] * pa_x[i];

        tr_y_xxxyyy_yzzzzz[i] = 2.0 * tr_y_xyyy_yzzzzz[i] * fe_0 + tr_y_xxyyy_yzzzzz[i] * pa_x[i];

        tr_y_xxxyyy_zzzzzz[i] = 2.0 * tr_y_xyyy_zzzzzz[i] * fe_0 + tr_y_xxyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 980-1008 components of targeted buffer : II

    auto tr_y_xxxyyz_xxxxxx = pbuffer.data(idx_dip_ii + 980);

    auto tr_y_xxxyyz_xxxxxy = pbuffer.data(idx_dip_ii + 981);

    auto tr_y_xxxyyz_xxxxxz = pbuffer.data(idx_dip_ii + 982);

    auto tr_y_xxxyyz_xxxxyy = pbuffer.data(idx_dip_ii + 983);

    auto tr_y_xxxyyz_xxxxyz = pbuffer.data(idx_dip_ii + 984);

    auto tr_y_xxxyyz_xxxxzz = pbuffer.data(idx_dip_ii + 985);

    auto tr_y_xxxyyz_xxxyyy = pbuffer.data(idx_dip_ii + 986);

    auto tr_y_xxxyyz_xxxyyz = pbuffer.data(idx_dip_ii + 987);

    auto tr_y_xxxyyz_xxxyzz = pbuffer.data(idx_dip_ii + 988);

    auto tr_y_xxxyyz_xxxzzz = pbuffer.data(idx_dip_ii + 989);

    auto tr_y_xxxyyz_xxyyyy = pbuffer.data(idx_dip_ii + 990);

    auto tr_y_xxxyyz_xxyyyz = pbuffer.data(idx_dip_ii + 991);

    auto tr_y_xxxyyz_xxyyzz = pbuffer.data(idx_dip_ii + 992);

    auto tr_y_xxxyyz_xxyzzz = pbuffer.data(idx_dip_ii + 993);

    auto tr_y_xxxyyz_xxzzzz = pbuffer.data(idx_dip_ii + 994);

    auto tr_y_xxxyyz_xyyyyy = pbuffer.data(idx_dip_ii + 995);

    auto tr_y_xxxyyz_xyyyyz = pbuffer.data(idx_dip_ii + 996);

    auto tr_y_xxxyyz_xyyyzz = pbuffer.data(idx_dip_ii + 997);

    auto tr_y_xxxyyz_xyyzzz = pbuffer.data(idx_dip_ii + 998);

    auto tr_y_xxxyyz_xyzzzz = pbuffer.data(idx_dip_ii + 999);

    auto tr_y_xxxyyz_xzzzzz = pbuffer.data(idx_dip_ii + 1000);

    auto tr_y_xxxyyz_yyyyyy = pbuffer.data(idx_dip_ii + 1001);

    auto tr_y_xxxyyz_yyyyyz = pbuffer.data(idx_dip_ii + 1002);

    auto tr_y_xxxyyz_yyyyzz = pbuffer.data(idx_dip_ii + 1003);

    auto tr_y_xxxyyz_yyyzzz = pbuffer.data(idx_dip_ii + 1004);

    auto tr_y_xxxyyz_yyzzzz = pbuffer.data(idx_dip_ii + 1005);

    auto tr_y_xxxyyz_yzzzzz = pbuffer.data(idx_dip_ii + 1006);

    auto tr_y_xxxyyz_zzzzzz = pbuffer.data(idx_dip_ii + 1007);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxyy_xxxxx, tr_y_xxxyy_xxxxxx, tr_y_xxxyy_xxxxxy, tr_y_xxxyy_xxxxxz, tr_y_xxxyy_xxxxy, tr_y_xxxyy_xxxxyy, tr_y_xxxyy_xxxxyz, tr_y_xxxyy_xxxxz, tr_y_xxxyy_xxxxzz, tr_y_xxxyy_xxxyy, tr_y_xxxyy_xxxyyy, tr_y_xxxyy_xxxyyz, tr_y_xxxyy_xxxyz, tr_y_xxxyy_xxxyzz, tr_y_xxxyy_xxxzz, tr_y_xxxyy_xxxzzz, tr_y_xxxyy_xxyyy, tr_y_xxxyy_xxyyyy, tr_y_xxxyy_xxyyyz, tr_y_xxxyy_xxyyz, tr_y_xxxyy_xxyyzz, tr_y_xxxyy_xxyzz, tr_y_xxxyy_xxyzzz, tr_y_xxxyy_xxzzz, tr_y_xxxyy_xxzzzz, tr_y_xxxyy_xyyyy, tr_y_xxxyy_xyyyyy, tr_y_xxxyy_xyyyyz, tr_y_xxxyy_xyyyz, tr_y_xxxyy_xyyyzz, tr_y_xxxyy_xyyzz, tr_y_xxxyy_xyyzzz, tr_y_xxxyy_xyzzz, tr_y_xxxyy_xyzzzz, tr_y_xxxyy_xzzzz, tr_y_xxxyy_xzzzzz, tr_y_xxxyy_yyyyyy, tr_y_xxxyyz_xxxxxx, tr_y_xxxyyz_xxxxxy, tr_y_xxxyyz_xxxxxz, tr_y_xxxyyz_xxxxyy, tr_y_xxxyyz_xxxxyz, tr_y_xxxyyz_xxxxzz, tr_y_xxxyyz_xxxyyy, tr_y_xxxyyz_xxxyyz, tr_y_xxxyyz_xxxyzz, tr_y_xxxyyz_xxxzzz, tr_y_xxxyyz_xxyyyy, tr_y_xxxyyz_xxyyyz, tr_y_xxxyyz_xxyyzz, tr_y_xxxyyz_xxyzzz, tr_y_xxxyyz_xxzzzz, tr_y_xxxyyz_xyyyyy, tr_y_xxxyyz_xyyyyz, tr_y_xxxyyz_xyyyzz, tr_y_xxxyyz_xyyzzz, tr_y_xxxyyz_xyzzzz, tr_y_xxxyyz_xzzzzz, tr_y_xxxyyz_yyyyyy, tr_y_xxxyyz_yyyyyz, tr_y_xxxyyz_yyyyzz, tr_y_xxxyyz_yyyzzz, tr_y_xxxyyz_yyzzzz, tr_y_xxxyyz_yzzzzz, tr_y_xxxyyz_zzzzzz, tr_y_xxyyz_yyyyyz, tr_y_xxyyz_yyyyzz, tr_y_xxyyz_yyyzzz, tr_y_xxyyz_yyzzzz, tr_y_xxyyz_yzzzzz, tr_y_xxyyz_zzzzzz, tr_y_xyyz_yyyyyz, tr_y_xyyz_yyyyzz, tr_y_xyyz_yyyzzz, tr_y_xyyz_yyzzzz, tr_y_xyyz_yzzzzz, tr_y_xyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyz_xxxxxx[i] = tr_y_xxxyy_xxxxxx[i] * pa_z[i];

        tr_y_xxxyyz_xxxxxy[i] = tr_y_xxxyy_xxxxxy[i] * pa_z[i];

        tr_y_xxxyyz_xxxxxz[i] = tr_y_xxxyy_xxxxx[i] * fe_0 + tr_y_xxxyy_xxxxxz[i] * pa_z[i];

        tr_y_xxxyyz_xxxxyy[i] = tr_y_xxxyy_xxxxyy[i] * pa_z[i];

        tr_y_xxxyyz_xxxxyz[i] = tr_y_xxxyy_xxxxy[i] * fe_0 + tr_y_xxxyy_xxxxyz[i] * pa_z[i];

        tr_y_xxxyyz_xxxxzz[i] = 2.0 * tr_y_xxxyy_xxxxz[i] * fe_0 + tr_y_xxxyy_xxxxzz[i] * pa_z[i];

        tr_y_xxxyyz_xxxyyy[i] = tr_y_xxxyy_xxxyyy[i] * pa_z[i];

        tr_y_xxxyyz_xxxyyz[i] = tr_y_xxxyy_xxxyy[i] * fe_0 + tr_y_xxxyy_xxxyyz[i] * pa_z[i];

        tr_y_xxxyyz_xxxyzz[i] = 2.0 * tr_y_xxxyy_xxxyz[i] * fe_0 + tr_y_xxxyy_xxxyzz[i] * pa_z[i];

        tr_y_xxxyyz_xxxzzz[i] = 3.0 * tr_y_xxxyy_xxxzz[i] * fe_0 + tr_y_xxxyy_xxxzzz[i] * pa_z[i];

        tr_y_xxxyyz_xxyyyy[i] = tr_y_xxxyy_xxyyyy[i] * pa_z[i];

        tr_y_xxxyyz_xxyyyz[i] = tr_y_xxxyy_xxyyy[i] * fe_0 + tr_y_xxxyy_xxyyyz[i] * pa_z[i];

        tr_y_xxxyyz_xxyyzz[i] = 2.0 * tr_y_xxxyy_xxyyz[i] * fe_0 + tr_y_xxxyy_xxyyzz[i] * pa_z[i];

        tr_y_xxxyyz_xxyzzz[i] = 3.0 * tr_y_xxxyy_xxyzz[i] * fe_0 + tr_y_xxxyy_xxyzzz[i] * pa_z[i];

        tr_y_xxxyyz_xxzzzz[i] = 4.0 * tr_y_xxxyy_xxzzz[i] * fe_0 + tr_y_xxxyy_xxzzzz[i] * pa_z[i];

        tr_y_xxxyyz_xyyyyy[i] = tr_y_xxxyy_xyyyyy[i] * pa_z[i];

        tr_y_xxxyyz_xyyyyz[i] = tr_y_xxxyy_xyyyy[i] * fe_0 + tr_y_xxxyy_xyyyyz[i] * pa_z[i];

        tr_y_xxxyyz_xyyyzz[i] = 2.0 * tr_y_xxxyy_xyyyz[i] * fe_0 + tr_y_xxxyy_xyyyzz[i] * pa_z[i];

        tr_y_xxxyyz_xyyzzz[i] = 3.0 * tr_y_xxxyy_xyyzz[i] * fe_0 + tr_y_xxxyy_xyyzzz[i] * pa_z[i];

        tr_y_xxxyyz_xyzzzz[i] = 4.0 * tr_y_xxxyy_xyzzz[i] * fe_0 + tr_y_xxxyy_xyzzzz[i] * pa_z[i];

        tr_y_xxxyyz_xzzzzz[i] = 5.0 * tr_y_xxxyy_xzzzz[i] * fe_0 + tr_y_xxxyy_xzzzzz[i] * pa_z[i];

        tr_y_xxxyyz_yyyyyy[i] = tr_y_xxxyy_yyyyyy[i] * pa_z[i];

        tr_y_xxxyyz_yyyyyz[i] = 2.0 * tr_y_xyyz_yyyyyz[i] * fe_0 + tr_y_xxyyz_yyyyyz[i] * pa_x[i];

        tr_y_xxxyyz_yyyyzz[i] = 2.0 * tr_y_xyyz_yyyyzz[i] * fe_0 + tr_y_xxyyz_yyyyzz[i] * pa_x[i];

        tr_y_xxxyyz_yyyzzz[i] = 2.0 * tr_y_xyyz_yyyzzz[i] * fe_0 + tr_y_xxyyz_yyyzzz[i] * pa_x[i];

        tr_y_xxxyyz_yyzzzz[i] = 2.0 * tr_y_xyyz_yyzzzz[i] * fe_0 + tr_y_xxyyz_yyzzzz[i] * pa_x[i];

        tr_y_xxxyyz_yzzzzz[i] = 2.0 * tr_y_xyyz_yzzzzz[i] * fe_0 + tr_y_xxyyz_yzzzzz[i] * pa_x[i];

        tr_y_xxxyyz_zzzzzz[i] = 2.0 * tr_y_xyyz_zzzzzz[i] * fe_0 + tr_y_xxyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1008-1036 components of targeted buffer : II

    auto tr_y_xxxyzz_xxxxxx = pbuffer.data(idx_dip_ii + 1008);

    auto tr_y_xxxyzz_xxxxxy = pbuffer.data(idx_dip_ii + 1009);

    auto tr_y_xxxyzz_xxxxxz = pbuffer.data(idx_dip_ii + 1010);

    auto tr_y_xxxyzz_xxxxyy = pbuffer.data(idx_dip_ii + 1011);

    auto tr_y_xxxyzz_xxxxyz = pbuffer.data(idx_dip_ii + 1012);

    auto tr_y_xxxyzz_xxxxzz = pbuffer.data(idx_dip_ii + 1013);

    auto tr_y_xxxyzz_xxxyyy = pbuffer.data(idx_dip_ii + 1014);

    auto tr_y_xxxyzz_xxxyyz = pbuffer.data(idx_dip_ii + 1015);

    auto tr_y_xxxyzz_xxxyzz = pbuffer.data(idx_dip_ii + 1016);

    auto tr_y_xxxyzz_xxxzzz = pbuffer.data(idx_dip_ii + 1017);

    auto tr_y_xxxyzz_xxyyyy = pbuffer.data(idx_dip_ii + 1018);

    auto tr_y_xxxyzz_xxyyyz = pbuffer.data(idx_dip_ii + 1019);

    auto tr_y_xxxyzz_xxyyzz = pbuffer.data(idx_dip_ii + 1020);

    auto tr_y_xxxyzz_xxyzzz = pbuffer.data(idx_dip_ii + 1021);

    auto tr_y_xxxyzz_xxzzzz = pbuffer.data(idx_dip_ii + 1022);

    auto tr_y_xxxyzz_xyyyyy = pbuffer.data(idx_dip_ii + 1023);

    auto tr_y_xxxyzz_xyyyyz = pbuffer.data(idx_dip_ii + 1024);

    auto tr_y_xxxyzz_xyyyzz = pbuffer.data(idx_dip_ii + 1025);

    auto tr_y_xxxyzz_xyyzzz = pbuffer.data(idx_dip_ii + 1026);

    auto tr_y_xxxyzz_xyzzzz = pbuffer.data(idx_dip_ii + 1027);

    auto tr_y_xxxyzz_xzzzzz = pbuffer.data(idx_dip_ii + 1028);

    auto tr_y_xxxyzz_yyyyyy = pbuffer.data(idx_dip_ii + 1029);

    auto tr_y_xxxyzz_yyyyyz = pbuffer.data(idx_dip_ii + 1030);

    auto tr_y_xxxyzz_yyyyzz = pbuffer.data(idx_dip_ii + 1031);

    auto tr_y_xxxyzz_yyyzzz = pbuffer.data(idx_dip_ii + 1032);

    auto tr_y_xxxyzz_yyzzzz = pbuffer.data(idx_dip_ii + 1033);

    auto tr_y_xxxyzz_yzzzzz = pbuffer.data(idx_dip_ii + 1034);

    auto tr_y_xxxyzz_zzzzzz = pbuffer.data(idx_dip_ii + 1035);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxy_xxxxxy, tr_y_xxxy_xxxxyy, tr_y_xxxy_xxxyyy, tr_y_xxxy_xxyyyy, tr_y_xxxy_xyyyyy, tr_y_xxxyz_xxxxxy, tr_y_xxxyz_xxxxyy, tr_y_xxxyz_xxxyyy, tr_y_xxxyz_xxyyyy, tr_y_xxxyz_xyyyyy, tr_y_xxxyzz_xxxxxx, tr_y_xxxyzz_xxxxxy, tr_y_xxxyzz_xxxxxz, tr_y_xxxyzz_xxxxyy, tr_y_xxxyzz_xxxxyz, tr_y_xxxyzz_xxxxzz, tr_y_xxxyzz_xxxyyy, tr_y_xxxyzz_xxxyyz, tr_y_xxxyzz_xxxyzz, tr_y_xxxyzz_xxxzzz, tr_y_xxxyzz_xxyyyy, tr_y_xxxyzz_xxyyyz, tr_y_xxxyzz_xxyyzz, tr_y_xxxyzz_xxyzzz, tr_y_xxxyzz_xxzzzz, tr_y_xxxyzz_xyyyyy, tr_y_xxxyzz_xyyyyz, tr_y_xxxyzz_xyyyzz, tr_y_xxxyzz_xyyzzz, tr_y_xxxyzz_xyzzzz, tr_y_xxxyzz_xzzzzz, tr_y_xxxyzz_yyyyyy, tr_y_xxxyzz_yyyyyz, tr_y_xxxyzz_yyyyzz, tr_y_xxxyzz_yyyzzz, tr_y_xxxyzz_yyzzzz, tr_y_xxxyzz_yzzzzz, tr_y_xxxyzz_zzzzzz, tr_y_xxxzz_xxxxxx, tr_y_xxxzz_xxxxxz, tr_y_xxxzz_xxxxzz, tr_y_xxxzz_xxxzzz, tr_y_xxxzz_xxzzzz, tr_y_xxxzz_xzzzzz, tr_y_xxyzz_xxxxyz, tr_y_xxyzz_xxxyyz, tr_y_xxyzz_xxxyz, tr_y_xxyzz_xxxyzz, tr_y_xxyzz_xxyyyz, tr_y_xxyzz_xxyyz, tr_y_xxyzz_xxyyzz, tr_y_xxyzz_xxyzz, tr_y_xxyzz_xxyzzz, tr_y_xxyzz_xyyyyz, tr_y_xxyzz_xyyyz, tr_y_xxyzz_xyyyzz, tr_y_xxyzz_xyyzz, tr_y_xxyzz_xyyzzz, tr_y_xxyzz_xyzzz, tr_y_xxyzz_xyzzzz, tr_y_xxyzz_yyyyyy, tr_y_xxyzz_yyyyyz, tr_y_xxyzz_yyyyz, tr_y_xxyzz_yyyyzz, tr_y_xxyzz_yyyzz, tr_y_xxyzz_yyyzzz, tr_y_xxyzz_yyzzz, tr_y_xxyzz_yyzzzz, tr_y_xxyzz_yzzzz, tr_y_xxyzz_yzzzzz, tr_y_xxyzz_zzzzzz, tr_y_xyzz_xxxxyz, tr_y_xyzz_xxxyyz, tr_y_xyzz_xxxyzz, tr_y_xyzz_xxyyyz, tr_y_xyzz_xxyyzz, tr_y_xyzz_xxyzzz, tr_y_xyzz_xyyyyz, tr_y_xyzz_xyyyzz, tr_y_xyzz_xyyzzz, tr_y_xyzz_xyzzzz, tr_y_xyzz_yyyyyy, tr_y_xyzz_yyyyyz, tr_y_xyzz_yyyyzz, tr_y_xyzz_yyyzzz, tr_y_xyzz_yyzzzz, tr_y_xyzz_yzzzzz, tr_y_xyzz_zzzzzz, ts_xxxzz_xxxxxx, ts_xxxzz_xxxxxz, ts_xxxzz_xxxxzz, ts_xxxzz_xxxzzz, ts_xxxzz_xxzzzz, ts_xxxzz_xzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyzz_xxxxxx[i] = ts_xxxzz_xxxxxx[i] * fe_0 + tr_y_xxxzz_xxxxxx[i] * pa_y[i];

        tr_y_xxxyzz_xxxxxy[i] = tr_y_xxxy_xxxxxy[i] * fe_0 + tr_y_xxxyz_xxxxxy[i] * pa_z[i];

        tr_y_xxxyzz_xxxxxz[i] = ts_xxxzz_xxxxxz[i] * fe_0 + tr_y_xxxzz_xxxxxz[i] * pa_y[i];

        tr_y_xxxyzz_xxxxyy[i] = tr_y_xxxy_xxxxyy[i] * fe_0 + tr_y_xxxyz_xxxxyy[i] * pa_z[i];

        tr_y_xxxyzz_xxxxyz[i] = 2.0 * tr_y_xyzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxyzz_xxxyz[i] * fe_0 + tr_y_xxyzz_xxxxyz[i] * pa_x[i];

        tr_y_xxxyzz_xxxxzz[i] = ts_xxxzz_xxxxzz[i] * fe_0 + tr_y_xxxzz_xxxxzz[i] * pa_y[i];

        tr_y_xxxyzz_xxxyyy[i] = tr_y_xxxy_xxxyyy[i] * fe_0 + tr_y_xxxyz_xxxyyy[i] * pa_z[i];

        tr_y_xxxyzz_xxxyyz[i] = 2.0 * tr_y_xyzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxyzz_xxyyz[i] * fe_0 + tr_y_xxyzz_xxxyyz[i] * pa_x[i];

        tr_y_xxxyzz_xxxyzz[i] = 2.0 * tr_y_xyzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxyzz_xxyzz[i] * fe_0 + tr_y_xxyzz_xxxyzz[i] * pa_x[i];

        tr_y_xxxyzz_xxxzzz[i] = ts_xxxzz_xxxzzz[i] * fe_0 + tr_y_xxxzz_xxxzzz[i] * pa_y[i];

        tr_y_xxxyzz_xxyyyy[i] = tr_y_xxxy_xxyyyy[i] * fe_0 + tr_y_xxxyz_xxyyyy[i] * pa_z[i];

        tr_y_xxxyzz_xxyyyz[i] = 2.0 * tr_y_xyzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxyzz_xyyyz[i] * fe_0 + tr_y_xxyzz_xxyyyz[i] * pa_x[i];

        tr_y_xxxyzz_xxyyzz[i] = 2.0 * tr_y_xyzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxyzz_xyyzz[i] * fe_0 + tr_y_xxyzz_xxyyzz[i] * pa_x[i];

        tr_y_xxxyzz_xxyzzz[i] = 2.0 * tr_y_xyzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxyzz_xyzzz[i] * fe_0 + tr_y_xxyzz_xxyzzz[i] * pa_x[i];

        tr_y_xxxyzz_xxzzzz[i] = ts_xxxzz_xxzzzz[i] * fe_0 + tr_y_xxxzz_xxzzzz[i] * pa_y[i];

        tr_y_xxxyzz_xyyyyy[i] = tr_y_xxxy_xyyyyy[i] * fe_0 + tr_y_xxxyz_xyyyyy[i] * pa_z[i];

        tr_y_xxxyzz_xyyyyz[i] = 2.0 * tr_y_xyzz_xyyyyz[i] * fe_0 + tr_y_xxyzz_yyyyz[i] * fe_0 + tr_y_xxyzz_xyyyyz[i] * pa_x[i];

        tr_y_xxxyzz_xyyyzz[i] = 2.0 * tr_y_xyzz_xyyyzz[i] * fe_0 + tr_y_xxyzz_yyyzz[i] * fe_0 + tr_y_xxyzz_xyyyzz[i] * pa_x[i];

        tr_y_xxxyzz_xyyzzz[i] = 2.0 * tr_y_xyzz_xyyzzz[i] * fe_0 + tr_y_xxyzz_yyzzz[i] * fe_0 + tr_y_xxyzz_xyyzzz[i] * pa_x[i];

        tr_y_xxxyzz_xyzzzz[i] = 2.0 * tr_y_xyzz_xyzzzz[i] * fe_0 + tr_y_xxyzz_yzzzz[i] * fe_0 + tr_y_xxyzz_xyzzzz[i] * pa_x[i];

        tr_y_xxxyzz_xzzzzz[i] = ts_xxxzz_xzzzzz[i] * fe_0 + tr_y_xxxzz_xzzzzz[i] * pa_y[i];

        tr_y_xxxyzz_yyyyyy[i] = 2.0 * tr_y_xyzz_yyyyyy[i] * fe_0 + tr_y_xxyzz_yyyyyy[i] * pa_x[i];

        tr_y_xxxyzz_yyyyyz[i] = 2.0 * tr_y_xyzz_yyyyyz[i] * fe_0 + tr_y_xxyzz_yyyyyz[i] * pa_x[i];

        tr_y_xxxyzz_yyyyzz[i] = 2.0 * tr_y_xyzz_yyyyzz[i] * fe_0 + tr_y_xxyzz_yyyyzz[i] * pa_x[i];

        tr_y_xxxyzz_yyyzzz[i] = 2.0 * tr_y_xyzz_yyyzzz[i] * fe_0 + tr_y_xxyzz_yyyzzz[i] * pa_x[i];

        tr_y_xxxyzz_yyzzzz[i] = 2.0 * tr_y_xyzz_yyzzzz[i] * fe_0 + tr_y_xxyzz_yyzzzz[i] * pa_x[i];

        tr_y_xxxyzz_yzzzzz[i] = 2.0 * tr_y_xyzz_yzzzzz[i] * fe_0 + tr_y_xxyzz_yzzzzz[i] * pa_x[i];

        tr_y_xxxyzz_zzzzzz[i] = 2.0 * tr_y_xyzz_zzzzzz[i] * fe_0 + tr_y_xxyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1036-1064 components of targeted buffer : II

    auto tr_y_xxxzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1036);

    auto tr_y_xxxzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1037);

    auto tr_y_xxxzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1038);

    auto tr_y_xxxzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1039);

    auto tr_y_xxxzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1040);

    auto tr_y_xxxzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1041);

    auto tr_y_xxxzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1042);

    auto tr_y_xxxzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1043);

    auto tr_y_xxxzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1044);

    auto tr_y_xxxzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1045);

    auto tr_y_xxxzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1046);

    auto tr_y_xxxzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1047);

    auto tr_y_xxxzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1048);

    auto tr_y_xxxzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1049);

    auto tr_y_xxxzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1050);

    auto tr_y_xxxzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1051);

    auto tr_y_xxxzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1052);

    auto tr_y_xxxzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1053);

    auto tr_y_xxxzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1054);

    auto tr_y_xxxzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1055);

    auto tr_y_xxxzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1056);

    auto tr_y_xxxzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1057);

    auto tr_y_xxxzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1058);

    auto tr_y_xxxzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1059);

    auto tr_y_xxxzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1060);

    auto tr_y_xxxzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1061);

    auto tr_y_xxxzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1062);

    auto tr_y_xxxzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1063);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxz_xxxxxx, tr_y_xxxz_xxxxxy, tr_y_xxxz_xxxxyy, tr_y_xxxz_xxxyyy, tr_y_xxxz_xxyyyy, tr_y_xxxz_xyyyyy, tr_y_xxxzz_xxxxxx, tr_y_xxxzz_xxxxxy, tr_y_xxxzz_xxxxyy, tr_y_xxxzz_xxxyyy, tr_y_xxxzz_xxyyyy, tr_y_xxxzz_xyyyyy, tr_y_xxxzzz_xxxxxx, tr_y_xxxzzz_xxxxxy, tr_y_xxxzzz_xxxxxz, tr_y_xxxzzz_xxxxyy, tr_y_xxxzzz_xxxxyz, tr_y_xxxzzz_xxxxzz, tr_y_xxxzzz_xxxyyy, tr_y_xxxzzz_xxxyyz, tr_y_xxxzzz_xxxyzz, tr_y_xxxzzz_xxxzzz, tr_y_xxxzzz_xxyyyy, tr_y_xxxzzz_xxyyyz, tr_y_xxxzzz_xxyyzz, tr_y_xxxzzz_xxyzzz, tr_y_xxxzzz_xxzzzz, tr_y_xxxzzz_xyyyyy, tr_y_xxxzzz_xyyyyz, tr_y_xxxzzz_xyyyzz, tr_y_xxxzzz_xyyzzz, tr_y_xxxzzz_xyzzzz, tr_y_xxxzzz_xzzzzz, tr_y_xxxzzz_yyyyyy, tr_y_xxxzzz_yyyyyz, tr_y_xxxzzz_yyyyzz, tr_y_xxxzzz_yyyzzz, tr_y_xxxzzz_yyzzzz, tr_y_xxxzzz_yzzzzz, tr_y_xxxzzz_zzzzzz, tr_y_xxzzz_xxxxxz, tr_y_xxzzz_xxxxyz, tr_y_xxzzz_xxxxz, tr_y_xxzzz_xxxxzz, tr_y_xxzzz_xxxyyz, tr_y_xxzzz_xxxyz, tr_y_xxzzz_xxxyzz, tr_y_xxzzz_xxxzz, tr_y_xxzzz_xxxzzz, tr_y_xxzzz_xxyyyz, tr_y_xxzzz_xxyyz, tr_y_xxzzz_xxyyzz, tr_y_xxzzz_xxyzz, tr_y_xxzzz_xxyzzz, tr_y_xxzzz_xxzzz, tr_y_xxzzz_xxzzzz, tr_y_xxzzz_xyyyyz, tr_y_xxzzz_xyyyz, tr_y_xxzzz_xyyyzz, tr_y_xxzzz_xyyzz, tr_y_xxzzz_xyyzzz, tr_y_xxzzz_xyzzz, tr_y_xxzzz_xyzzzz, tr_y_xxzzz_xzzzz, tr_y_xxzzz_xzzzzz, tr_y_xxzzz_yyyyyy, tr_y_xxzzz_yyyyyz, tr_y_xxzzz_yyyyz, tr_y_xxzzz_yyyyzz, tr_y_xxzzz_yyyzz, tr_y_xxzzz_yyyzzz, tr_y_xxzzz_yyzzz, tr_y_xxzzz_yyzzzz, tr_y_xxzzz_yzzzz, tr_y_xxzzz_yzzzzz, tr_y_xxzzz_zzzzz, tr_y_xxzzz_zzzzzz, tr_y_xzzz_xxxxxz, tr_y_xzzz_xxxxyz, tr_y_xzzz_xxxxzz, tr_y_xzzz_xxxyyz, tr_y_xzzz_xxxyzz, tr_y_xzzz_xxxzzz, tr_y_xzzz_xxyyyz, tr_y_xzzz_xxyyzz, tr_y_xzzz_xxyzzz, tr_y_xzzz_xxzzzz, tr_y_xzzz_xyyyyz, tr_y_xzzz_xyyyzz, tr_y_xzzz_xyyzzz, tr_y_xzzz_xyzzzz, tr_y_xzzz_xzzzzz, tr_y_xzzz_yyyyyy, tr_y_xzzz_yyyyyz, tr_y_xzzz_yyyyzz, tr_y_xzzz_yyyzzz, tr_y_xzzz_yyzzzz, tr_y_xzzz_yzzzzz, tr_y_xzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzzz_xxxxxx[i] = 2.0 * tr_y_xxxz_xxxxxx[i] * fe_0 + tr_y_xxxzz_xxxxxx[i] * pa_z[i];

        tr_y_xxxzzz_xxxxxy[i] = 2.0 * tr_y_xxxz_xxxxxy[i] * fe_0 + tr_y_xxxzz_xxxxxy[i] * pa_z[i];

        tr_y_xxxzzz_xxxxxz[i] = 2.0 * tr_y_xzzz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxzzz_xxxxz[i] * fe_0 + tr_y_xxzzz_xxxxxz[i] * pa_x[i];

        tr_y_xxxzzz_xxxxyy[i] = 2.0 * tr_y_xxxz_xxxxyy[i] * fe_0 + tr_y_xxxzz_xxxxyy[i] * pa_z[i];

        tr_y_xxxzzz_xxxxyz[i] = 2.0 * tr_y_xzzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxzzz_xxxyz[i] * fe_0 + tr_y_xxzzz_xxxxyz[i] * pa_x[i];

        tr_y_xxxzzz_xxxxzz[i] = 2.0 * tr_y_xzzz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxzzz_xxxzz[i] * fe_0 + tr_y_xxzzz_xxxxzz[i] * pa_x[i];

        tr_y_xxxzzz_xxxyyy[i] = 2.0 * tr_y_xxxz_xxxyyy[i] * fe_0 + tr_y_xxxzz_xxxyyy[i] * pa_z[i];

        tr_y_xxxzzz_xxxyyz[i] = 2.0 * tr_y_xzzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxzzz_xxyyz[i] * fe_0 + tr_y_xxzzz_xxxyyz[i] * pa_x[i];

        tr_y_xxxzzz_xxxyzz[i] = 2.0 * tr_y_xzzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxzzz_xxyzz[i] * fe_0 + tr_y_xxzzz_xxxyzz[i] * pa_x[i];

        tr_y_xxxzzz_xxxzzz[i] = 2.0 * tr_y_xzzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxzzz_xxzzz[i] * fe_0 + tr_y_xxzzz_xxxzzz[i] * pa_x[i];

        tr_y_xxxzzz_xxyyyy[i] = 2.0 * tr_y_xxxz_xxyyyy[i] * fe_0 + tr_y_xxxzz_xxyyyy[i] * pa_z[i];

        tr_y_xxxzzz_xxyyyz[i] = 2.0 * tr_y_xzzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxzzz_xyyyz[i] * fe_0 + tr_y_xxzzz_xxyyyz[i] * pa_x[i];

        tr_y_xxxzzz_xxyyzz[i] = 2.0 * tr_y_xzzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxzzz_xyyzz[i] * fe_0 + tr_y_xxzzz_xxyyzz[i] * pa_x[i];

        tr_y_xxxzzz_xxyzzz[i] = 2.0 * tr_y_xzzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxzzz_xyzzz[i] * fe_0 + tr_y_xxzzz_xxyzzz[i] * pa_x[i];

        tr_y_xxxzzz_xxzzzz[i] = 2.0 * tr_y_xzzz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxzzz_xzzzz[i] * fe_0 + tr_y_xxzzz_xxzzzz[i] * pa_x[i];

        tr_y_xxxzzz_xyyyyy[i] = 2.0 * tr_y_xxxz_xyyyyy[i] * fe_0 + tr_y_xxxzz_xyyyyy[i] * pa_z[i];

        tr_y_xxxzzz_xyyyyz[i] = 2.0 * tr_y_xzzz_xyyyyz[i] * fe_0 + tr_y_xxzzz_yyyyz[i] * fe_0 + tr_y_xxzzz_xyyyyz[i] * pa_x[i];

        tr_y_xxxzzz_xyyyzz[i] = 2.0 * tr_y_xzzz_xyyyzz[i] * fe_0 + tr_y_xxzzz_yyyzz[i] * fe_0 + tr_y_xxzzz_xyyyzz[i] * pa_x[i];

        tr_y_xxxzzz_xyyzzz[i] = 2.0 * tr_y_xzzz_xyyzzz[i] * fe_0 + tr_y_xxzzz_yyzzz[i] * fe_0 + tr_y_xxzzz_xyyzzz[i] * pa_x[i];

        tr_y_xxxzzz_xyzzzz[i] = 2.0 * tr_y_xzzz_xyzzzz[i] * fe_0 + tr_y_xxzzz_yzzzz[i] * fe_0 + tr_y_xxzzz_xyzzzz[i] * pa_x[i];

        tr_y_xxxzzz_xzzzzz[i] = 2.0 * tr_y_xzzz_xzzzzz[i] * fe_0 + tr_y_xxzzz_zzzzz[i] * fe_0 + tr_y_xxzzz_xzzzzz[i] * pa_x[i];

        tr_y_xxxzzz_yyyyyy[i] = 2.0 * tr_y_xzzz_yyyyyy[i] * fe_0 + tr_y_xxzzz_yyyyyy[i] * pa_x[i];

        tr_y_xxxzzz_yyyyyz[i] = 2.0 * tr_y_xzzz_yyyyyz[i] * fe_0 + tr_y_xxzzz_yyyyyz[i] * pa_x[i];

        tr_y_xxxzzz_yyyyzz[i] = 2.0 * tr_y_xzzz_yyyyzz[i] * fe_0 + tr_y_xxzzz_yyyyzz[i] * pa_x[i];

        tr_y_xxxzzz_yyyzzz[i] = 2.0 * tr_y_xzzz_yyyzzz[i] * fe_0 + tr_y_xxzzz_yyyzzz[i] * pa_x[i];

        tr_y_xxxzzz_yyzzzz[i] = 2.0 * tr_y_xzzz_yyzzzz[i] * fe_0 + tr_y_xxzzz_yyzzzz[i] * pa_x[i];

        tr_y_xxxzzz_yzzzzz[i] = 2.0 * tr_y_xzzz_yzzzzz[i] * fe_0 + tr_y_xxzzz_yzzzzz[i] * pa_x[i];

        tr_y_xxxzzz_zzzzzz[i] = 2.0 * tr_y_xzzz_zzzzzz[i] * fe_0 + tr_y_xxzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1064-1092 components of targeted buffer : II

    auto tr_y_xxyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 1064);

    auto tr_y_xxyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 1065);

    auto tr_y_xxyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 1066);

    auto tr_y_xxyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 1067);

    auto tr_y_xxyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 1068);

    auto tr_y_xxyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 1069);

    auto tr_y_xxyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 1070);

    auto tr_y_xxyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 1071);

    auto tr_y_xxyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 1072);

    auto tr_y_xxyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 1073);

    auto tr_y_xxyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 1074);

    auto tr_y_xxyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 1075);

    auto tr_y_xxyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 1076);

    auto tr_y_xxyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 1077);

    auto tr_y_xxyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 1078);

    auto tr_y_xxyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 1079);

    auto tr_y_xxyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 1080);

    auto tr_y_xxyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 1081);

    auto tr_y_xxyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 1082);

    auto tr_y_xxyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 1083);

    auto tr_y_xxyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 1084);

    auto tr_y_xxyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 1085);

    auto tr_y_xxyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 1086);

    auto tr_y_xxyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 1087);

    auto tr_y_xxyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 1088);

    auto tr_y_xxyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 1089);

    auto tr_y_xxyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 1090);

    auto tr_y_xxyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 1091);

    #pragma omp simd aligned(pa_x, tr_y_xxyyyy_xxxxxx, tr_y_xxyyyy_xxxxxy, tr_y_xxyyyy_xxxxxz, tr_y_xxyyyy_xxxxyy, tr_y_xxyyyy_xxxxyz, tr_y_xxyyyy_xxxxzz, tr_y_xxyyyy_xxxyyy, tr_y_xxyyyy_xxxyyz, tr_y_xxyyyy_xxxyzz, tr_y_xxyyyy_xxxzzz, tr_y_xxyyyy_xxyyyy, tr_y_xxyyyy_xxyyyz, tr_y_xxyyyy_xxyyzz, tr_y_xxyyyy_xxyzzz, tr_y_xxyyyy_xxzzzz, tr_y_xxyyyy_xyyyyy, tr_y_xxyyyy_xyyyyz, tr_y_xxyyyy_xyyyzz, tr_y_xxyyyy_xyyzzz, tr_y_xxyyyy_xyzzzz, tr_y_xxyyyy_xzzzzz, tr_y_xxyyyy_yyyyyy, tr_y_xxyyyy_yyyyyz, tr_y_xxyyyy_yyyyzz, tr_y_xxyyyy_yyyzzz, tr_y_xxyyyy_yyzzzz, tr_y_xxyyyy_yzzzzz, tr_y_xxyyyy_zzzzzz, tr_y_xyyyy_xxxxx, tr_y_xyyyy_xxxxxx, tr_y_xyyyy_xxxxxy, tr_y_xyyyy_xxxxxz, tr_y_xyyyy_xxxxy, tr_y_xyyyy_xxxxyy, tr_y_xyyyy_xxxxyz, tr_y_xyyyy_xxxxz, tr_y_xyyyy_xxxxzz, tr_y_xyyyy_xxxyy, tr_y_xyyyy_xxxyyy, tr_y_xyyyy_xxxyyz, tr_y_xyyyy_xxxyz, tr_y_xyyyy_xxxyzz, tr_y_xyyyy_xxxzz, tr_y_xyyyy_xxxzzz, tr_y_xyyyy_xxyyy, tr_y_xyyyy_xxyyyy, tr_y_xyyyy_xxyyyz, tr_y_xyyyy_xxyyz, tr_y_xyyyy_xxyyzz, tr_y_xyyyy_xxyzz, tr_y_xyyyy_xxyzzz, tr_y_xyyyy_xxzzz, tr_y_xyyyy_xxzzzz, tr_y_xyyyy_xyyyy, tr_y_xyyyy_xyyyyy, tr_y_xyyyy_xyyyyz, tr_y_xyyyy_xyyyz, tr_y_xyyyy_xyyyzz, tr_y_xyyyy_xyyzz, tr_y_xyyyy_xyyzzz, tr_y_xyyyy_xyzzz, tr_y_xyyyy_xyzzzz, tr_y_xyyyy_xzzzz, tr_y_xyyyy_xzzzzz, tr_y_xyyyy_yyyyy, tr_y_xyyyy_yyyyyy, tr_y_xyyyy_yyyyyz, tr_y_xyyyy_yyyyz, tr_y_xyyyy_yyyyzz, tr_y_xyyyy_yyyzz, tr_y_xyyyy_yyyzzz, tr_y_xyyyy_yyzzz, tr_y_xyyyy_yyzzzz, tr_y_xyyyy_yzzzz, tr_y_xyyyy_yzzzzz, tr_y_xyyyy_zzzzz, tr_y_xyyyy_zzzzzz, tr_y_yyyy_xxxxxx, tr_y_yyyy_xxxxxy, tr_y_yyyy_xxxxxz, tr_y_yyyy_xxxxyy, tr_y_yyyy_xxxxyz, tr_y_yyyy_xxxxzz, tr_y_yyyy_xxxyyy, tr_y_yyyy_xxxyyz, tr_y_yyyy_xxxyzz, tr_y_yyyy_xxxzzz, tr_y_yyyy_xxyyyy, tr_y_yyyy_xxyyyz, tr_y_yyyy_xxyyzz, tr_y_yyyy_xxyzzz, tr_y_yyyy_xxzzzz, tr_y_yyyy_xyyyyy, tr_y_yyyy_xyyyyz, tr_y_yyyy_xyyyzz, tr_y_yyyy_xyyzzz, tr_y_yyyy_xyzzzz, tr_y_yyyy_xzzzzz, tr_y_yyyy_yyyyyy, tr_y_yyyy_yyyyyz, tr_y_yyyy_yyyyzz, tr_y_yyyy_yyyzzz, tr_y_yyyy_yyzzzz, tr_y_yyyy_yzzzzz, tr_y_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyy_xxxxxx[i] = tr_y_yyyy_xxxxxx[i] * fe_0 + 6.0 * tr_y_xyyyy_xxxxx[i] * fe_0 + tr_y_xyyyy_xxxxxx[i] * pa_x[i];

        tr_y_xxyyyy_xxxxxy[i] = tr_y_yyyy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xyyyy_xxxxy[i] * fe_0 + tr_y_xyyyy_xxxxxy[i] * pa_x[i];

        tr_y_xxyyyy_xxxxxz[i] = tr_y_yyyy_xxxxxz[i] * fe_0 + 5.0 * tr_y_xyyyy_xxxxz[i] * fe_0 + tr_y_xyyyy_xxxxxz[i] * pa_x[i];

        tr_y_xxyyyy_xxxxyy[i] = tr_y_yyyy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xyyyy_xxxyy[i] * fe_0 + tr_y_xyyyy_xxxxyy[i] * pa_x[i];

        tr_y_xxyyyy_xxxxyz[i] = tr_y_yyyy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xyyyy_xxxyz[i] * fe_0 + tr_y_xyyyy_xxxxyz[i] * pa_x[i];

        tr_y_xxyyyy_xxxxzz[i] = tr_y_yyyy_xxxxzz[i] * fe_0 + 4.0 * tr_y_xyyyy_xxxzz[i] * fe_0 + tr_y_xyyyy_xxxxzz[i] * pa_x[i];

        tr_y_xxyyyy_xxxyyy[i] = tr_y_yyyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xyyyy_xxyyy[i] * fe_0 + tr_y_xyyyy_xxxyyy[i] * pa_x[i];

        tr_y_xxyyyy_xxxyyz[i] = tr_y_yyyy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xyyyy_xxyyz[i] * fe_0 + tr_y_xyyyy_xxxyyz[i] * pa_x[i];

        tr_y_xxyyyy_xxxyzz[i] = tr_y_yyyy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xyyyy_xxyzz[i] * fe_0 + tr_y_xyyyy_xxxyzz[i] * pa_x[i];

        tr_y_xxyyyy_xxxzzz[i] = tr_y_yyyy_xxxzzz[i] * fe_0 + 3.0 * tr_y_xyyyy_xxzzz[i] * fe_0 + tr_y_xyyyy_xxxzzz[i] * pa_x[i];

        tr_y_xxyyyy_xxyyyy[i] = tr_y_yyyy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xyyyy_xyyyy[i] * fe_0 + tr_y_xyyyy_xxyyyy[i] * pa_x[i];

        tr_y_xxyyyy_xxyyyz[i] = tr_y_yyyy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xyyyy_xyyyz[i] * fe_0 + tr_y_xyyyy_xxyyyz[i] * pa_x[i];

        tr_y_xxyyyy_xxyyzz[i] = tr_y_yyyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xyyyy_xyyzz[i] * fe_0 + tr_y_xyyyy_xxyyzz[i] * pa_x[i];

        tr_y_xxyyyy_xxyzzz[i] = tr_y_yyyy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xyyyy_xyzzz[i] * fe_0 + tr_y_xyyyy_xxyzzz[i] * pa_x[i];

        tr_y_xxyyyy_xxzzzz[i] = tr_y_yyyy_xxzzzz[i] * fe_0 + 2.0 * tr_y_xyyyy_xzzzz[i] * fe_0 + tr_y_xyyyy_xxzzzz[i] * pa_x[i];

        tr_y_xxyyyy_xyyyyy[i] = tr_y_yyyy_xyyyyy[i] * fe_0 + tr_y_xyyyy_yyyyy[i] * fe_0 + tr_y_xyyyy_xyyyyy[i] * pa_x[i];

        tr_y_xxyyyy_xyyyyz[i] = tr_y_yyyy_xyyyyz[i] * fe_0 + tr_y_xyyyy_yyyyz[i] * fe_0 + tr_y_xyyyy_xyyyyz[i] * pa_x[i];

        tr_y_xxyyyy_xyyyzz[i] = tr_y_yyyy_xyyyzz[i] * fe_0 + tr_y_xyyyy_yyyzz[i] * fe_0 + tr_y_xyyyy_xyyyzz[i] * pa_x[i];

        tr_y_xxyyyy_xyyzzz[i] = tr_y_yyyy_xyyzzz[i] * fe_0 + tr_y_xyyyy_yyzzz[i] * fe_0 + tr_y_xyyyy_xyyzzz[i] * pa_x[i];

        tr_y_xxyyyy_xyzzzz[i] = tr_y_yyyy_xyzzzz[i] * fe_0 + tr_y_xyyyy_yzzzz[i] * fe_0 + tr_y_xyyyy_xyzzzz[i] * pa_x[i];

        tr_y_xxyyyy_xzzzzz[i] = tr_y_yyyy_xzzzzz[i] * fe_0 + tr_y_xyyyy_zzzzz[i] * fe_0 + tr_y_xyyyy_xzzzzz[i] * pa_x[i];

        tr_y_xxyyyy_yyyyyy[i] = tr_y_yyyy_yyyyyy[i] * fe_0 + tr_y_xyyyy_yyyyyy[i] * pa_x[i];

        tr_y_xxyyyy_yyyyyz[i] = tr_y_yyyy_yyyyyz[i] * fe_0 + tr_y_xyyyy_yyyyyz[i] * pa_x[i];

        tr_y_xxyyyy_yyyyzz[i] = tr_y_yyyy_yyyyzz[i] * fe_0 + tr_y_xyyyy_yyyyzz[i] * pa_x[i];

        tr_y_xxyyyy_yyyzzz[i] = tr_y_yyyy_yyyzzz[i] * fe_0 + tr_y_xyyyy_yyyzzz[i] * pa_x[i];

        tr_y_xxyyyy_yyzzzz[i] = tr_y_yyyy_yyzzzz[i] * fe_0 + tr_y_xyyyy_yyzzzz[i] * pa_x[i];

        tr_y_xxyyyy_yzzzzz[i] = tr_y_yyyy_yzzzzz[i] * fe_0 + tr_y_xyyyy_yzzzzz[i] * pa_x[i];

        tr_y_xxyyyy_zzzzzz[i] = tr_y_yyyy_zzzzzz[i] * fe_0 + tr_y_xyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1092-1120 components of targeted buffer : II

    auto tr_y_xxyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 1092);

    auto tr_y_xxyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 1093);

    auto tr_y_xxyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 1094);

    auto tr_y_xxyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 1095);

    auto tr_y_xxyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 1096);

    auto tr_y_xxyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 1097);

    auto tr_y_xxyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 1098);

    auto tr_y_xxyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 1099);

    auto tr_y_xxyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 1100);

    auto tr_y_xxyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 1101);

    auto tr_y_xxyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 1102);

    auto tr_y_xxyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 1103);

    auto tr_y_xxyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 1104);

    auto tr_y_xxyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 1105);

    auto tr_y_xxyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 1106);

    auto tr_y_xxyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 1107);

    auto tr_y_xxyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 1108);

    auto tr_y_xxyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 1109);

    auto tr_y_xxyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 1110);

    auto tr_y_xxyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 1111);

    auto tr_y_xxyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 1112);

    auto tr_y_xxyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 1113);

    auto tr_y_xxyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 1114);

    auto tr_y_xxyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 1115);

    auto tr_y_xxyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 1116);

    auto tr_y_xxyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 1117);

    auto tr_y_xxyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 1118);

    auto tr_y_xxyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 1119);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyyy_xxxxx, tr_y_xxyyy_xxxxxx, tr_y_xxyyy_xxxxxy, tr_y_xxyyy_xxxxxz, tr_y_xxyyy_xxxxy, tr_y_xxyyy_xxxxyy, tr_y_xxyyy_xxxxyz, tr_y_xxyyy_xxxxz, tr_y_xxyyy_xxxxzz, tr_y_xxyyy_xxxyy, tr_y_xxyyy_xxxyyy, tr_y_xxyyy_xxxyyz, tr_y_xxyyy_xxxyz, tr_y_xxyyy_xxxyzz, tr_y_xxyyy_xxxzz, tr_y_xxyyy_xxxzzz, tr_y_xxyyy_xxyyy, tr_y_xxyyy_xxyyyy, tr_y_xxyyy_xxyyyz, tr_y_xxyyy_xxyyz, tr_y_xxyyy_xxyyzz, tr_y_xxyyy_xxyzz, tr_y_xxyyy_xxyzzz, tr_y_xxyyy_xxzzz, tr_y_xxyyy_xxzzzz, tr_y_xxyyy_xyyyy, tr_y_xxyyy_xyyyyy, tr_y_xxyyy_xyyyyz, tr_y_xxyyy_xyyyz, tr_y_xxyyy_xyyyzz, tr_y_xxyyy_xyyzz, tr_y_xxyyy_xyyzzz, tr_y_xxyyy_xyzzz, tr_y_xxyyy_xyzzzz, tr_y_xxyyy_xzzzz, tr_y_xxyyy_xzzzzz, tr_y_xxyyy_yyyyyy, tr_y_xxyyyz_xxxxxx, tr_y_xxyyyz_xxxxxy, tr_y_xxyyyz_xxxxxz, tr_y_xxyyyz_xxxxyy, tr_y_xxyyyz_xxxxyz, tr_y_xxyyyz_xxxxzz, tr_y_xxyyyz_xxxyyy, tr_y_xxyyyz_xxxyyz, tr_y_xxyyyz_xxxyzz, tr_y_xxyyyz_xxxzzz, tr_y_xxyyyz_xxyyyy, tr_y_xxyyyz_xxyyyz, tr_y_xxyyyz_xxyyzz, tr_y_xxyyyz_xxyzzz, tr_y_xxyyyz_xxzzzz, tr_y_xxyyyz_xyyyyy, tr_y_xxyyyz_xyyyyz, tr_y_xxyyyz_xyyyzz, tr_y_xxyyyz_xyyzzz, tr_y_xxyyyz_xyzzzz, tr_y_xxyyyz_xzzzzz, tr_y_xxyyyz_yyyyyy, tr_y_xxyyyz_yyyyyz, tr_y_xxyyyz_yyyyzz, tr_y_xxyyyz_yyyzzz, tr_y_xxyyyz_yyzzzz, tr_y_xxyyyz_yzzzzz, tr_y_xxyyyz_zzzzzz, tr_y_xyyyz_yyyyyz, tr_y_xyyyz_yyyyzz, tr_y_xyyyz_yyyzzz, tr_y_xyyyz_yyzzzz, tr_y_xyyyz_yzzzzz, tr_y_xyyyz_zzzzzz, tr_y_yyyz_yyyyyz, tr_y_yyyz_yyyyzz, tr_y_yyyz_yyyzzz, tr_y_yyyz_yyzzzz, tr_y_yyyz_yzzzzz, tr_y_yyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyz_xxxxxx[i] = tr_y_xxyyy_xxxxxx[i] * pa_z[i];

        tr_y_xxyyyz_xxxxxy[i] = tr_y_xxyyy_xxxxxy[i] * pa_z[i];

        tr_y_xxyyyz_xxxxxz[i] = tr_y_xxyyy_xxxxx[i] * fe_0 + tr_y_xxyyy_xxxxxz[i] * pa_z[i];

        tr_y_xxyyyz_xxxxyy[i] = tr_y_xxyyy_xxxxyy[i] * pa_z[i];

        tr_y_xxyyyz_xxxxyz[i] = tr_y_xxyyy_xxxxy[i] * fe_0 + tr_y_xxyyy_xxxxyz[i] * pa_z[i];

        tr_y_xxyyyz_xxxxzz[i] = 2.0 * tr_y_xxyyy_xxxxz[i] * fe_0 + tr_y_xxyyy_xxxxzz[i] * pa_z[i];

        tr_y_xxyyyz_xxxyyy[i] = tr_y_xxyyy_xxxyyy[i] * pa_z[i];

        tr_y_xxyyyz_xxxyyz[i] = tr_y_xxyyy_xxxyy[i] * fe_0 + tr_y_xxyyy_xxxyyz[i] * pa_z[i];

        tr_y_xxyyyz_xxxyzz[i] = 2.0 * tr_y_xxyyy_xxxyz[i] * fe_0 + tr_y_xxyyy_xxxyzz[i] * pa_z[i];

        tr_y_xxyyyz_xxxzzz[i] = 3.0 * tr_y_xxyyy_xxxzz[i] * fe_0 + tr_y_xxyyy_xxxzzz[i] * pa_z[i];

        tr_y_xxyyyz_xxyyyy[i] = tr_y_xxyyy_xxyyyy[i] * pa_z[i];

        tr_y_xxyyyz_xxyyyz[i] = tr_y_xxyyy_xxyyy[i] * fe_0 + tr_y_xxyyy_xxyyyz[i] * pa_z[i];

        tr_y_xxyyyz_xxyyzz[i] = 2.0 * tr_y_xxyyy_xxyyz[i] * fe_0 + tr_y_xxyyy_xxyyzz[i] * pa_z[i];

        tr_y_xxyyyz_xxyzzz[i] = 3.0 * tr_y_xxyyy_xxyzz[i] * fe_0 + tr_y_xxyyy_xxyzzz[i] * pa_z[i];

        tr_y_xxyyyz_xxzzzz[i] = 4.0 * tr_y_xxyyy_xxzzz[i] * fe_0 + tr_y_xxyyy_xxzzzz[i] * pa_z[i];

        tr_y_xxyyyz_xyyyyy[i] = tr_y_xxyyy_xyyyyy[i] * pa_z[i];

        tr_y_xxyyyz_xyyyyz[i] = tr_y_xxyyy_xyyyy[i] * fe_0 + tr_y_xxyyy_xyyyyz[i] * pa_z[i];

        tr_y_xxyyyz_xyyyzz[i] = 2.0 * tr_y_xxyyy_xyyyz[i] * fe_0 + tr_y_xxyyy_xyyyzz[i] * pa_z[i];

        tr_y_xxyyyz_xyyzzz[i] = 3.0 * tr_y_xxyyy_xyyzz[i] * fe_0 + tr_y_xxyyy_xyyzzz[i] * pa_z[i];

        tr_y_xxyyyz_xyzzzz[i] = 4.0 * tr_y_xxyyy_xyzzz[i] * fe_0 + tr_y_xxyyy_xyzzzz[i] * pa_z[i];

        tr_y_xxyyyz_xzzzzz[i] = 5.0 * tr_y_xxyyy_xzzzz[i] * fe_0 + tr_y_xxyyy_xzzzzz[i] * pa_z[i];

        tr_y_xxyyyz_yyyyyy[i] = tr_y_xxyyy_yyyyyy[i] * pa_z[i];

        tr_y_xxyyyz_yyyyyz[i] = tr_y_yyyz_yyyyyz[i] * fe_0 + tr_y_xyyyz_yyyyyz[i] * pa_x[i];

        tr_y_xxyyyz_yyyyzz[i] = tr_y_yyyz_yyyyzz[i] * fe_0 + tr_y_xyyyz_yyyyzz[i] * pa_x[i];

        tr_y_xxyyyz_yyyzzz[i] = tr_y_yyyz_yyyzzz[i] * fe_0 + tr_y_xyyyz_yyyzzz[i] * pa_x[i];

        tr_y_xxyyyz_yyzzzz[i] = tr_y_yyyz_yyzzzz[i] * fe_0 + tr_y_xyyyz_yyzzzz[i] * pa_x[i];

        tr_y_xxyyyz_yzzzzz[i] = tr_y_yyyz_yzzzzz[i] * fe_0 + tr_y_xyyyz_yzzzzz[i] * pa_x[i];

        tr_y_xxyyyz_zzzzzz[i] = tr_y_yyyz_zzzzzz[i] * fe_0 + tr_y_xyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1120-1148 components of targeted buffer : II

    auto tr_y_xxyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 1120);

    auto tr_y_xxyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 1121);

    auto tr_y_xxyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 1122);

    auto tr_y_xxyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 1123);

    auto tr_y_xxyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 1124);

    auto tr_y_xxyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 1125);

    auto tr_y_xxyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 1126);

    auto tr_y_xxyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 1127);

    auto tr_y_xxyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 1128);

    auto tr_y_xxyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 1129);

    auto tr_y_xxyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 1130);

    auto tr_y_xxyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 1131);

    auto tr_y_xxyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 1132);

    auto tr_y_xxyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 1133);

    auto tr_y_xxyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 1134);

    auto tr_y_xxyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 1135);

    auto tr_y_xxyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 1136);

    auto tr_y_xxyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 1137);

    auto tr_y_xxyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 1138);

    auto tr_y_xxyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 1139);

    auto tr_y_xxyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 1140);

    auto tr_y_xxyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 1141);

    auto tr_y_xxyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 1142);

    auto tr_y_xxyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 1143);

    auto tr_y_xxyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 1144);

    auto tr_y_xxyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 1145);

    auto tr_y_xxyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 1146);

    auto tr_y_xxyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 1147);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyy_xxxxxx, tr_y_xxyy_xxxxxy, tr_y_xxyy_xxxxyy, tr_y_xxyy_xxxyyy, tr_y_xxyy_xxyyyy, tr_y_xxyy_xyyyyy, tr_y_xxyyz_xxxxxx, tr_y_xxyyz_xxxxxy, tr_y_xxyyz_xxxxyy, tr_y_xxyyz_xxxyyy, tr_y_xxyyz_xxyyyy, tr_y_xxyyz_xyyyyy, tr_y_xxyyzz_xxxxxx, tr_y_xxyyzz_xxxxxy, tr_y_xxyyzz_xxxxxz, tr_y_xxyyzz_xxxxyy, tr_y_xxyyzz_xxxxyz, tr_y_xxyyzz_xxxxzz, tr_y_xxyyzz_xxxyyy, tr_y_xxyyzz_xxxyyz, tr_y_xxyyzz_xxxyzz, tr_y_xxyyzz_xxxzzz, tr_y_xxyyzz_xxyyyy, tr_y_xxyyzz_xxyyyz, tr_y_xxyyzz_xxyyzz, tr_y_xxyyzz_xxyzzz, tr_y_xxyyzz_xxzzzz, tr_y_xxyyzz_xyyyyy, tr_y_xxyyzz_xyyyyz, tr_y_xxyyzz_xyyyzz, tr_y_xxyyzz_xyyzzz, tr_y_xxyyzz_xyzzzz, tr_y_xxyyzz_xzzzzz, tr_y_xxyyzz_yyyyyy, tr_y_xxyyzz_yyyyyz, tr_y_xxyyzz_yyyyzz, tr_y_xxyyzz_yyyzzz, tr_y_xxyyzz_yyzzzz, tr_y_xxyyzz_yzzzzz, tr_y_xxyyzz_zzzzzz, tr_y_xyyzz_xxxxxz, tr_y_xyyzz_xxxxyz, tr_y_xyyzz_xxxxz, tr_y_xyyzz_xxxxzz, tr_y_xyyzz_xxxyyz, tr_y_xyyzz_xxxyz, tr_y_xyyzz_xxxyzz, tr_y_xyyzz_xxxzz, tr_y_xyyzz_xxxzzz, tr_y_xyyzz_xxyyyz, tr_y_xyyzz_xxyyz, tr_y_xyyzz_xxyyzz, tr_y_xyyzz_xxyzz, tr_y_xyyzz_xxyzzz, tr_y_xyyzz_xxzzz, tr_y_xyyzz_xxzzzz, tr_y_xyyzz_xyyyyz, tr_y_xyyzz_xyyyz, tr_y_xyyzz_xyyyzz, tr_y_xyyzz_xyyzz, tr_y_xyyzz_xyyzzz, tr_y_xyyzz_xyzzz, tr_y_xyyzz_xyzzzz, tr_y_xyyzz_xzzzz, tr_y_xyyzz_xzzzzz, tr_y_xyyzz_yyyyyy, tr_y_xyyzz_yyyyyz, tr_y_xyyzz_yyyyz, tr_y_xyyzz_yyyyzz, tr_y_xyyzz_yyyzz, tr_y_xyyzz_yyyzzz, tr_y_xyyzz_yyzzz, tr_y_xyyzz_yyzzzz, tr_y_xyyzz_yzzzz, tr_y_xyyzz_yzzzzz, tr_y_xyyzz_zzzzz, tr_y_xyyzz_zzzzzz, tr_y_yyzz_xxxxxz, tr_y_yyzz_xxxxyz, tr_y_yyzz_xxxxzz, tr_y_yyzz_xxxyyz, tr_y_yyzz_xxxyzz, tr_y_yyzz_xxxzzz, tr_y_yyzz_xxyyyz, tr_y_yyzz_xxyyzz, tr_y_yyzz_xxyzzz, tr_y_yyzz_xxzzzz, tr_y_yyzz_xyyyyz, tr_y_yyzz_xyyyzz, tr_y_yyzz_xyyzzz, tr_y_yyzz_xyzzzz, tr_y_yyzz_xzzzzz, tr_y_yyzz_yyyyyy, tr_y_yyzz_yyyyyz, tr_y_yyzz_yyyyzz, tr_y_yyzz_yyyzzz, tr_y_yyzz_yyzzzz, tr_y_yyzz_yzzzzz, tr_y_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyzz_xxxxxx[i] = tr_y_xxyy_xxxxxx[i] * fe_0 + tr_y_xxyyz_xxxxxx[i] * pa_z[i];

        tr_y_xxyyzz_xxxxxy[i] = tr_y_xxyy_xxxxxy[i] * fe_0 + tr_y_xxyyz_xxxxxy[i] * pa_z[i];

        tr_y_xxyyzz_xxxxxz[i] = tr_y_yyzz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xyyzz_xxxxz[i] * fe_0 + tr_y_xyyzz_xxxxxz[i] * pa_x[i];

        tr_y_xxyyzz_xxxxyy[i] = tr_y_xxyy_xxxxyy[i] * fe_0 + tr_y_xxyyz_xxxxyy[i] * pa_z[i];

        tr_y_xxyyzz_xxxxyz[i] = tr_y_yyzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xyyzz_xxxyz[i] * fe_0 + tr_y_xyyzz_xxxxyz[i] * pa_x[i];

        tr_y_xxyyzz_xxxxzz[i] = tr_y_yyzz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xyyzz_xxxzz[i] * fe_0 + tr_y_xyyzz_xxxxzz[i] * pa_x[i];

        tr_y_xxyyzz_xxxyyy[i] = tr_y_xxyy_xxxyyy[i] * fe_0 + tr_y_xxyyz_xxxyyy[i] * pa_z[i];

        tr_y_xxyyzz_xxxyyz[i] = tr_y_yyzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xyyzz_xxyyz[i] * fe_0 + tr_y_xyyzz_xxxyyz[i] * pa_x[i];

        tr_y_xxyyzz_xxxyzz[i] = tr_y_yyzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xyyzz_xxyzz[i] * fe_0 + tr_y_xyyzz_xxxyzz[i] * pa_x[i];

        tr_y_xxyyzz_xxxzzz[i] = tr_y_yyzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xyyzz_xxzzz[i] * fe_0 + tr_y_xyyzz_xxxzzz[i] * pa_x[i];

        tr_y_xxyyzz_xxyyyy[i] = tr_y_xxyy_xxyyyy[i] * fe_0 + tr_y_xxyyz_xxyyyy[i] * pa_z[i];

        tr_y_xxyyzz_xxyyyz[i] = tr_y_yyzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xyyzz_xyyyz[i] * fe_0 + tr_y_xyyzz_xxyyyz[i] * pa_x[i];

        tr_y_xxyyzz_xxyyzz[i] = tr_y_yyzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xyyzz_xyyzz[i] * fe_0 + tr_y_xyyzz_xxyyzz[i] * pa_x[i];

        tr_y_xxyyzz_xxyzzz[i] = tr_y_yyzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xyyzz_xyzzz[i] * fe_0 + tr_y_xyyzz_xxyzzz[i] * pa_x[i];

        tr_y_xxyyzz_xxzzzz[i] = tr_y_yyzz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xyyzz_xzzzz[i] * fe_0 + tr_y_xyyzz_xxzzzz[i] * pa_x[i];

        tr_y_xxyyzz_xyyyyy[i] = tr_y_xxyy_xyyyyy[i] * fe_0 + tr_y_xxyyz_xyyyyy[i] * pa_z[i];

        tr_y_xxyyzz_xyyyyz[i] = tr_y_yyzz_xyyyyz[i] * fe_0 + tr_y_xyyzz_yyyyz[i] * fe_0 + tr_y_xyyzz_xyyyyz[i] * pa_x[i];

        tr_y_xxyyzz_xyyyzz[i] = tr_y_yyzz_xyyyzz[i] * fe_0 + tr_y_xyyzz_yyyzz[i] * fe_0 + tr_y_xyyzz_xyyyzz[i] * pa_x[i];

        tr_y_xxyyzz_xyyzzz[i] = tr_y_yyzz_xyyzzz[i] * fe_0 + tr_y_xyyzz_yyzzz[i] * fe_0 + tr_y_xyyzz_xyyzzz[i] * pa_x[i];

        tr_y_xxyyzz_xyzzzz[i] = tr_y_yyzz_xyzzzz[i] * fe_0 + tr_y_xyyzz_yzzzz[i] * fe_0 + tr_y_xyyzz_xyzzzz[i] * pa_x[i];

        tr_y_xxyyzz_xzzzzz[i] = tr_y_yyzz_xzzzzz[i] * fe_0 + tr_y_xyyzz_zzzzz[i] * fe_0 + tr_y_xyyzz_xzzzzz[i] * pa_x[i];

        tr_y_xxyyzz_yyyyyy[i] = tr_y_yyzz_yyyyyy[i] * fe_0 + tr_y_xyyzz_yyyyyy[i] * pa_x[i];

        tr_y_xxyyzz_yyyyyz[i] = tr_y_yyzz_yyyyyz[i] * fe_0 + tr_y_xyyzz_yyyyyz[i] * pa_x[i];

        tr_y_xxyyzz_yyyyzz[i] = tr_y_yyzz_yyyyzz[i] * fe_0 + tr_y_xyyzz_yyyyzz[i] * pa_x[i];

        tr_y_xxyyzz_yyyzzz[i] = tr_y_yyzz_yyyzzz[i] * fe_0 + tr_y_xyyzz_yyyzzz[i] * pa_x[i];

        tr_y_xxyyzz_yyzzzz[i] = tr_y_yyzz_yyzzzz[i] * fe_0 + tr_y_xyyzz_yyzzzz[i] * pa_x[i];

        tr_y_xxyyzz_yzzzzz[i] = tr_y_yyzz_yzzzzz[i] * fe_0 + tr_y_xyyzz_yzzzzz[i] * pa_x[i];

        tr_y_xxyyzz_zzzzzz[i] = tr_y_yyzz_zzzzzz[i] * fe_0 + tr_y_xyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1148-1176 components of targeted buffer : II

    auto tr_y_xxyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1148);

    auto tr_y_xxyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1149);

    auto tr_y_xxyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1150);

    auto tr_y_xxyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1151);

    auto tr_y_xxyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1152);

    auto tr_y_xxyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1153);

    auto tr_y_xxyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1154);

    auto tr_y_xxyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1155);

    auto tr_y_xxyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1156);

    auto tr_y_xxyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1157);

    auto tr_y_xxyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1158);

    auto tr_y_xxyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1159);

    auto tr_y_xxyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1160);

    auto tr_y_xxyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1161);

    auto tr_y_xxyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1162);

    auto tr_y_xxyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1163);

    auto tr_y_xxyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1164);

    auto tr_y_xxyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1165);

    auto tr_y_xxyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1166);

    auto tr_y_xxyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1167);

    auto tr_y_xxyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1168);

    auto tr_y_xxyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1169);

    auto tr_y_xxyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1170);

    auto tr_y_xxyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1171);

    auto tr_y_xxyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1172);

    auto tr_y_xxyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1173);

    auto tr_y_xxyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1174);

    auto tr_y_xxyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1175);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxyz_xxxxxy, tr_y_xxyz_xxxxyy, tr_y_xxyz_xxxyyy, tr_y_xxyz_xxyyyy, tr_y_xxyz_xyyyyy, tr_y_xxyzz_xxxxxy, tr_y_xxyzz_xxxxyy, tr_y_xxyzz_xxxyyy, tr_y_xxyzz_xxyyyy, tr_y_xxyzz_xyyyyy, tr_y_xxyzzz_xxxxxx, tr_y_xxyzzz_xxxxxy, tr_y_xxyzzz_xxxxxz, tr_y_xxyzzz_xxxxyy, tr_y_xxyzzz_xxxxyz, tr_y_xxyzzz_xxxxzz, tr_y_xxyzzz_xxxyyy, tr_y_xxyzzz_xxxyyz, tr_y_xxyzzz_xxxyzz, tr_y_xxyzzz_xxxzzz, tr_y_xxyzzz_xxyyyy, tr_y_xxyzzz_xxyyyz, tr_y_xxyzzz_xxyyzz, tr_y_xxyzzz_xxyzzz, tr_y_xxyzzz_xxzzzz, tr_y_xxyzzz_xyyyyy, tr_y_xxyzzz_xyyyyz, tr_y_xxyzzz_xyyyzz, tr_y_xxyzzz_xyyzzz, tr_y_xxyzzz_xyzzzz, tr_y_xxyzzz_xzzzzz, tr_y_xxyzzz_yyyyyy, tr_y_xxyzzz_yyyyyz, tr_y_xxyzzz_yyyyzz, tr_y_xxyzzz_yyyzzz, tr_y_xxyzzz_yyzzzz, tr_y_xxyzzz_yzzzzz, tr_y_xxyzzz_zzzzzz, tr_y_xxzzz_xxxxxx, tr_y_xxzzz_xxxxxz, tr_y_xxzzz_xxxxzz, tr_y_xxzzz_xxxzzz, tr_y_xxzzz_xxzzzz, tr_y_xxzzz_xzzzzz, tr_y_xyzzz_xxxxyz, tr_y_xyzzz_xxxyyz, tr_y_xyzzz_xxxyz, tr_y_xyzzz_xxxyzz, tr_y_xyzzz_xxyyyz, tr_y_xyzzz_xxyyz, tr_y_xyzzz_xxyyzz, tr_y_xyzzz_xxyzz, tr_y_xyzzz_xxyzzz, tr_y_xyzzz_xyyyyz, tr_y_xyzzz_xyyyz, tr_y_xyzzz_xyyyzz, tr_y_xyzzz_xyyzz, tr_y_xyzzz_xyyzzz, tr_y_xyzzz_xyzzz, tr_y_xyzzz_xyzzzz, tr_y_xyzzz_yyyyyy, tr_y_xyzzz_yyyyyz, tr_y_xyzzz_yyyyz, tr_y_xyzzz_yyyyzz, tr_y_xyzzz_yyyzz, tr_y_xyzzz_yyyzzz, tr_y_xyzzz_yyzzz, tr_y_xyzzz_yyzzzz, tr_y_xyzzz_yzzzz, tr_y_xyzzz_yzzzzz, tr_y_xyzzz_zzzzzz, tr_y_yzzz_xxxxyz, tr_y_yzzz_xxxyyz, tr_y_yzzz_xxxyzz, tr_y_yzzz_xxyyyz, tr_y_yzzz_xxyyzz, tr_y_yzzz_xxyzzz, tr_y_yzzz_xyyyyz, tr_y_yzzz_xyyyzz, tr_y_yzzz_xyyzzz, tr_y_yzzz_xyzzzz, tr_y_yzzz_yyyyyy, tr_y_yzzz_yyyyyz, tr_y_yzzz_yyyyzz, tr_y_yzzz_yyyzzz, tr_y_yzzz_yyzzzz, tr_y_yzzz_yzzzzz, tr_y_yzzz_zzzzzz, ts_xxzzz_xxxxxx, ts_xxzzz_xxxxxz, ts_xxzzz_xxxxzz, ts_xxzzz_xxxzzz, ts_xxzzz_xxzzzz, ts_xxzzz_xzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzzz_xxxxxx[i] = ts_xxzzz_xxxxxx[i] * fe_0 + tr_y_xxzzz_xxxxxx[i] * pa_y[i];

        tr_y_xxyzzz_xxxxxy[i] = 2.0 * tr_y_xxyz_xxxxxy[i] * fe_0 + tr_y_xxyzz_xxxxxy[i] * pa_z[i];

        tr_y_xxyzzz_xxxxxz[i] = ts_xxzzz_xxxxxz[i] * fe_0 + tr_y_xxzzz_xxxxxz[i] * pa_y[i];

        tr_y_xxyzzz_xxxxyy[i] = 2.0 * tr_y_xxyz_xxxxyy[i] * fe_0 + tr_y_xxyzz_xxxxyy[i] * pa_z[i];

        tr_y_xxyzzz_xxxxyz[i] = tr_y_yzzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xyzzz_xxxyz[i] * fe_0 + tr_y_xyzzz_xxxxyz[i] * pa_x[i];

        tr_y_xxyzzz_xxxxzz[i] = ts_xxzzz_xxxxzz[i] * fe_0 + tr_y_xxzzz_xxxxzz[i] * pa_y[i];

        tr_y_xxyzzz_xxxyyy[i] = 2.0 * tr_y_xxyz_xxxyyy[i] * fe_0 + tr_y_xxyzz_xxxyyy[i] * pa_z[i];

        tr_y_xxyzzz_xxxyyz[i] = tr_y_yzzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xyzzz_xxyyz[i] * fe_0 + tr_y_xyzzz_xxxyyz[i] * pa_x[i];

        tr_y_xxyzzz_xxxyzz[i] = tr_y_yzzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xyzzz_xxyzz[i] * fe_0 + tr_y_xyzzz_xxxyzz[i] * pa_x[i];

        tr_y_xxyzzz_xxxzzz[i] = ts_xxzzz_xxxzzz[i] * fe_0 + tr_y_xxzzz_xxxzzz[i] * pa_y[i];

        tr_y_xxyzzz_xxyyyy[i] = 2.0 * tr_y_xxyz_xxyyyy[i] * fe_0 + tr_y_xxyzz_xxyyyy[i] * pa_z[i];

        tr_y_xxyzzz_xxyyyz[i] = tr_y_yzzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xyzzz_xyyyz[i] * fe_0 + tr_y_xyzzz_xxyyyz[i] * pa_x[i];

        tr_y_xxyzzz_xxyyzz[i] = tr_y_yzzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xyzzz_xyyzz[i] * fe_0 + tr_y_xyzzz_xxyyzz[i] * pa_x[i];

        tr_y_xxyzzz_xxyzzz[i] = tr_y_yzzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xyzzz_xyzzz[i] * fe_0 + tr_y_xyzzz_xxyzzz[i] * pa_x[i];

        tr_y_xxyzzz_xxzzzz[i] = ts_xxzzz_xxzzzz[i] * fe_0 + tr_y_xxzzz_xxzzzz[i] * pa_y[i];

        tr_y_xxyzzz_xyyyyy[i] = 2.0 * tr_y_xxyz_xyyyyy[i] * fe_0 + tr_y_xxyzz_xyyyyy[i] * pa_z[i];

        tr_y_xxyzzz_xyyyyz[i] = tr_y_yzzz_xyyyyz[i] * fe_0 + tr_y_xyzzz_yyyyz[i] * fe_0 + tr_y_xyzzz_xyyyyz[i] * pa_x[i];

        tr_y_xxyzzz_xyyyzz[i] = tr_y_yzzz_xyyyzz[i] * fe_0 + tr_y_xyzzz_yyyzz[i] * fe_0 + tr_y_xyzzz_xyyyzz[i] * pa_x[i];

        tr_y_xxyzzz_xyyzzz[i] = tr_y_yzzz_xyyzzz[i] * fe_0 + tr_y_xyzzz_yyzzz[i] * fe_0 + tr_y_xyzzz_xyyzzz[i] * pa_x[i];

        tr_y_xxyzzz_xyzzzz[i] = tr_y_yzzz_xyzzzz[i] * fe_0 + tr_y_xyzzz_yzzzz[i] * fe_0 + tr_y_xyzzz_xyzzzz[i] * pa_x[i];

        tr_y_xxyzzz_xzzzzz[i] = ts_xxzzz_xzzzzz[i] * fe_0 + tr_y_xxzzz_xzzzzz[i] * pa_y[i];

        tr_y_xxyzzz_yyyyyy[i] = tr_y_yzzz_yyyyyy[i] * fe_0 + tr_y_xyzzz_yyyyyy[i] * pa_x[i];

        tr_y_xxyzzz_yyyyyz[i] = tr_y_yzzz_yyyyyz[i] * fe_0 + tr_y_xyzzz_yyyyyz[i] * pa_x[i];

        tr_y_xxyzzz_yyyyzz[i] = tr_y_yzzz_yyyyzz[i] * fe_0 + tr_y_xyzzz_yyyyzz[i] * pa_x[i];

        tr_y_xxyzzz_yyyzzz[i] = tr_y_yzzz_yyyzzz[i] * fe_0 + tr_y_xyzzz_yyyzzz[i] * pa_x[i];

        tr_y_xxyzzz_yyzzzz[i] = tr_y_yzzz_yyzzzz[i] * fe_0 + tr_y_xyzzz_yyzzzz[i] * pa_x[i];

        tr_y_xxyzzz_yzzzzz[i] = tr_y_yzzz_yzzzzz[i] * fe_0 + tr_y_xyzzz_yzzzzz[i] * pa_x[i];

        tr_y_xxyzzz_zzzzzz[i] = tr_y_yzzz_zzzzzz[i] * fe_0 + tr_y_xyzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1176-1204 components of targeted buffer : II

    auto tr_y_xxzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1176);

    auto tr_y_xxzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1177);

    auto tr_y_xxzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1178);

    auto tr_y_xxzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1179);

    auto tr_y_xxzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1180);

    auto tr_y_xxzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1181);

    auto tr_y_xxzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1182);

    auto tr_y_xxzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1183);

    auto tr_y_xxzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1184);

    auto tr_y_xxzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1185);

    auto tr_y_xxzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1186);

    auto tr_y_xxzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1187);

    auto tr_y_xxzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1188);

    auto tr_y_xxzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1189);

    auto tr_y_xxzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1190);

    auto tr_y_xxzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1191);

    auto tr_y_xxzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1192);

    auto tr_y_xxzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1193);

    auto tr_y_xxzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1194);

    auto tr_y_xxzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1195);

    auto tr_y_xxzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1196);

    auto tr_y_xxzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1197);

    auto tr_y_xxzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1198);

    auto tr_y_xxzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1199);

    auto tr_y_xxzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1200);

    auto tr_y_xxzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1201);

    auto tr_y_xxzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1202);

    auto tr_y_xxzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1203);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxzz_xxxxxx, tr_y_xxzz_xxxxxy, tr_y_xxzz_xxxxyy, tr_y_xxzz_xxxyyy, tr_y_xxzz_xxyyyy, tr_y_xxzz_xyyyyy, tr_y_xxzzz_xxxxxx, tr_y_xxzzz_xxxxxy, tr_y_xxzzz_xxxxyy, tr_y_xxzzz_xxxyyy, tr_y_xxzzz_xxyyyy, tr_y_xxzzz_xyyyyy, tr_y_xxzzzz_xxxxxx, tr_y_xxzzzz_xxxxxy, tr_y_xxzzzz_xxxxxz, tr_y_xxzzzz_xxxxyy, tr_y_xxzzzz_xxxxyz, tr_y_xxzzzz_xxxxzz, tr_y_xxzzzz_xxxyyy, tr_y_xxzzzz_xxxyyz, tr_y_xxzzzz_xxxyzz, tr_y_xxzzzz_xxxzzz, tr_y_xxzzzz_xxyyyy, tr_y_xxzzzz_xxyyyz, tr_y_xxzzzz_xxyyzz, tr_y_xxzzzz_xxyzzz, tr_y_xxzzzz_xxzzzz, tr_y_xxzzzz_xyyyyy, tr_y_xxzzzz_xyyyyz, tr_y_xxzzzz_xyyyzz, tr_y_xxzzzz_xyyzzz, tr_y_xxzzzz_xyzzzz, tr_y_xxzzzz_xzzzzz, tr_y_xxzzzz_yyyyyy, tr_y_xxzzzz_yyyyyz, tr_y_xxzzzz_yyyyzz, tr_y_xxzzzz_yyyzzz, tr_y_xxzzzz_yyzzzz, tr_y_xxzzzz_yzzzzz, tr_y_xxzzzz_zzzzzz, tr_y_xzzzz_xxxxxz, tr_y_xzzzz_xxxxyz, tr_y_xzzzz_xxxxz, tr_y_xzzzz_xxxxzz, tr_y_xzzzz_xxxyyz, tr_y_xzzzz_xxxyz, tr_y_xzzzz_xxxyzz, tr_y_xzzzz_xxxzz, tr_y_xzzzz_xxxzzz, tr_y_xzzzz_xxyyyz, tr_y_xzzzz_xxyyz, tr_y_xzzzz_xxyyzz, tr_y_xzzzz_xxyzz, tr_y_xzzzz_xxyzzz, tr_y_xzzzz_xxzzz, tr_y_xzzzz_xxzzzz, tr_y_xzzzz_xyyyyz, tr_y_xzzzz_xyyyz, tr_y_xzzzz_xyyyzz, tr_y_xzzzz_xyyzz, tr_y_xzzzz_xyyzzz, tr_y_xzzzz_xyzzz, tr_y_xzzzz_xyzzzz, tr_y_xzzzz_xzzzz, tr_y_xzzzz_xzzzzz, tr_y_xzzzz_yyyyyy, tr_y_xzzzz_yyyyyz, tr_y_xzzzz_yyyyz, tr_y_xzzzz_yyyyzz, tr_y_xzzzz_yyyzz, tr_y_xzzzz_yyyzzz, tr_y_xzzzz_yyzzz, tr_y_xzzzz_yyzzzz, tr_y_xzzzz_yzzzz, tr_y_xzzzz_yzzzzz, tr_y_xzzzz_zzzzz, tr_y_xzzzz_zzzzzz, tr_y_zzzz_xxxxxz, tr_y_zzzz_xxxxyz, tr_y_zzzz_xxxxzz, tr_y_zzzz_xxxyyz, tr_y_zzzz_xxxyzz, tr_y_zzzz_xxxzzz, tr_y_zzzz_xxyyyz, tr_y_zzzz_xxyyzz, tr_y_zzzz_xxyzzz, tr_y_zzzz_xxzzzz, tr_y_zzzz_xyyyyz, tr_y_zzzz_xyyyzz, tr_y_zzzz_xyyzzz, tr_y_zzzz_xyzzzz, tr_y_zzzz_xzzzzz, tr_y_zzzz_yyyyyy, tr_y_zzzz_yyyyyz, tr_y_zzzz_yyyyzz, tr_y_zzzz_yyyzzz, tr_y_zzzz_yyzzzz, tr_y_zzzz_yzzzzz, tr_y_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzzz_xxxxxx[i] = 3.0 * tr_y_xxzz_xxxxxx[i] * fe_0 + tr_y_xxzzz_xxxxxx[i] * pa_z[i];

        tr_y_xxzzzz_xxxxxy[i] = 3.0 * tr_y_xxzz_xxxxxy[i] * fe_0 + tr_y_xxzzz_xxxxxy[i] * pa_z[i];

        tr_y_xxzzzz_xxxxxz[i] = tr_y_zzzz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xzzzz_xxxxz[i] * fe_0 + tr_y_xzzzz_xxxxxz[i] * pa_x[i];

        tr_y_xxzzzz_xxxxyy[i] = 3.0 * tr_y_xxzz_xxxxyy[i] * fe_0 + tr_y_xxzzz_xxxxyy[i] * pa_z[i];

        tr_y_xxzzzz_xxxxyz[i] = tr_y_zzzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xzzzz_xxxyz[i] * fe_0 + tr_y_xzzzz_xxxxyz[i] * pa_x[i];

        tr_y_xxzzzz_xxxxzz[i] = tr_y_zzzz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xzzzz_xxxzz[i] * fe_0 + tr_y_xzzzz_xxxxzz[i] * pa_x[i];

        tr_y_xxzzzz_xxxyyy[i] = 3.0 * tr_y_xxzz_xxxyyy[i] * fe_0 + tr_y_xxzzz_xxxyyy[i] * pa_z[i];

        tr_y_xxzzzz_xxxyyz[i] = tr_y_zzzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xzzzz_xxyyz[i] * fe_0 + tr_y_xzzzz_xxxyyz[i] * pa_x[i];

        tr_y_xxzzzz_xxxyzz[i] = tr_y_zzzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xzzzz_xxyzz[i] * fe_0 + tr_y_xzzzz_xxxyzz[i] * pa_x[i];

        tr_y_xxzzzz_xxxzzz[i] = tr_y_zzzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xzzzz_xxzzz[i] * fe_0 + tr_y_xzzzz_xxxzzz[i] * pa_x[i];

        tr_y_xxzzzz_xxyyyy[i] = 3.0 * tr_y_xxzz_xxyyyy[i] * fe_0 + tr_y_xxzzz_xxyyyy[i] * pa_z[i];

        tr_y_xxzzzz_xxyyyz[i] = tr_y_zzzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xzzzz_xyyyz[i] * fe_0 + tr_y_xzzzz_xxyyyz[i] * pa_x[i];

        tr_y_xxzzzz_xxyyzz[i] = tr_y_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xzzzz_xyyzz[i] * fe_0 + tr_y_xzzzz_xxyyzz[i] * pa_x[i];

        tr_y_xxzzzz_xxyzzz[i] = tr_y_zzzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xzzzz_xyzzz[i] * fe_0 + tr_y_xzzzz_xxyzzz[i] * pa_x[i];

        tr_y_xxzzzz_xxzzzz[i] = tr_y_zzzz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xzzzz_xzzzz[i] * fe_0 + tr_y_xzzzz_xxzzzz[i] * pa_x[i];

        tr_y_xxzzzz_xyyyyy[i] = 3.0 * tr_y_xxzz_xyyyyy[i] * fe_0 + tr_y_xxzzz_xyyyyy[i] * pa_z[i];

        tr_y_xxzzzz_xyyyyz[i] = tr_y_zzzz_xyyyyz[i] * fe_0 + tr_y_xzzzz_yyyyz[i] * fe_0 + tr_y_xzzzz_xyyyyz[i] * pa_x[i];

        tr_y_xxzzzz_xyyyzz[i] = tr_y_zzzz_xyyyzz[i] * fe_0 + tr_y_xzzzz_yyyzz[i] * fe_0 + tr_y_xzzzz_xyyyzz[i] * pa_x[i];

        tr_y_xxzzzz_xyyzzz[i] = tr_y_zzzz_xyyzzz[i] * fe_0 + tr_y_xzzzz_yyzzz[i] * fe_0 + tr_y_xzzzz_xyyzzz[i] * pa_x[i];

        tr_y_xxzzzz_xyzzzz[i] = tr_y_zzzz_xyzzzz[i] * fe_0 + tr_y_xzzzz_yzzzz[i] * fe_0 + tr_y_xzzzz_xyzzzz[i] * pa_x[i];

        tr_y_xxzzzz_xzzzzz[i] = tr_y_zzzz_xzzzzz[i] * fe_0 + tr_y_xzzzz_zzzzz[i] * fe_0 + tr_y_xzzzz_xzzzzz[i] * pa_x[i];

        tr_y_xxzzzz_yyyyyy[i] = tr_y_zzzz_yyyyyy[i] * fe_0 + tr_y_xzzzz_yyyyyy[i] * pa_x[i];

        tr_y_xxzzzz_yyyyyz[i] = tr_y_zzzz_yyyyyz[i] * fe_0 + tr_y_xzzzz_yyyyyz[i] * pa_x[i];

        tr_y_xxzzzz_yyyyzz[i] = tr_y_zzzz_yyyyzz[i] * fe_0 + tr_y_xzzzz_yyyyzz[i] * pa_x[i];

        tr_y_xxzzzz_yyyzzz[i] = tr_y_zzzz_yyyzzz[i] * fe_0 + tr_y_xzzzz_yyyzzz[i] * pa_x[i];

        tr_y_xxzzzz_yyzzzz[i] = tr_y_zzzz_yyzzzz[i] * fe_0 + tr_y_xzzzz_yyzzzz[i] * pa_x[i];

        tr_y_xxzzzz_yzzzzz[i] = tr_y_zzzz_yzzzzz[i] * fe_0 + tr_y_xzzzz_yzzzzz[i] * pa_x[i];

        tr_y_xxzzzz_zzzzzz[i] = tr_y_zzzz_zzzzzz[i] * fe_0 + tr_y_xzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1204-1232 components of targeted buffer : II

    auto tr_y_xyyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 1204);

    auto tr_y_xyyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 1205);

    auto tr_y_xyyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 1206);

    auto tr_y_xyyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 1207);

    auto tr_y_xyyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 1208);

    auto tr_y_xyyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 1209);

    auto tr_y_xyyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 1210);

    auto tr_y_xyyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 1211);

    auto tr_y_xyyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 1212);

    auto tr_y_xyyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 1213);

    auto tr_y_xyyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 1214);

    auto tr_y_xyyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 1215);

    auto tr_y_xyyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 1216);

    auto tr_y_xyyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 1217);

    auto tr_y_xyyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 1218);

    auto tr_y_xyyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 1219);

    auto tr_y_xyyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 1220);

    auto tr_y_xyyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 1221);

    auto tr_y_xyyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 1222);

    auto tr_y_xyyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 1223);

    auto tr_y_xyyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 1224);

    auto tr_y_xyyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 1225);

    auto tr_y_xyyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 1226);

    auto tr_y_xyyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 1227);

    auto tr_y_xyyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 1228);

    auto tr_y_xyyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 1229);

    auto tr_y_xyyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 1230);

    auto tr_y_xyyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 1231);

    #pragma omp simd aligned(pa_x, tr_y_xyyyyy_xxxxxx, tr_y_xyyyyy_xxxxxy, tr_y_xyyyyy_xxxxxz, tr_y_xyyyyy_xxxxyy, tr_y_xyyyyy_xxxxyz, tr_y_xyyyyy_xxxxzz, tr_y_xyyyyy_xxxyyy, tr_y_xyyyyy_xxxyyz, tr_y_xyyyyy_xxxyzz, tr_y_xyyyyy_xxxzzz, tr_y_xyyyyy_xxyyyy, tr_y_xyyyyy_xxyyyz, tr_y_xyyyyy_xxyyzz, tr_y_xyyyyy_xxyzzz, tr_y_xyyyyy_xxzzzz, tr_y_xyyyyy_xyyyyy, tr_y_xyyyyy_xyyyyz, tr_y_xyyyyy_xyyyzz, tr_y_xyyyyy_xyyzzz, tr_y_xyyyyy_xyzzzz, tr_y_xyyyyy_xzzzzz, tr_y_xyyyyy_yyyyyy, tr_y_xyyyyy_yyyyyz, tr_y_xyyyyy_yyyyzz, tr_y_xyyyyy_yyyzzz, tr_y_xyyyyy_yyzzzz, tr_y_xyyyyy_yzzzzz, tr_y_xyyyyy_zzzzzz, tr_y_yyyyy_xxxxx, tr_y_yyyyy_xxxxxx, tr_y_yyyyy_xxxxxy, tr_y_yyyyy_xxxxxz, tr_y_yyyyy_xxxxy, tr_y_yyyyy_xxxxyy, tr_y_yyyyy_xxxxyz, tr_y_yyyyy_xxxxz, tr_y_yyyyy_xxxxzz, tr_y_yyyyy_xxxyy, tr_y_yyyyy_xxxyyy, tr_y_yyyyy_xxxyyz, tr_y_yyyyy_xxxyz, tr_y_yyyyy_xxxyzz, tr_y_yyyyy_xxxzz, tr_y_yyyyy_xxxzzz, tr_y_yyyyy_xxyyy, tr_y_yyyyy_xxyyyy, tr_y_yyyyy_xxyyyz, tr_y_yyyyy_xxyyz, tr_y_yyyyy_xxyyzz, tr_y_yyyyy_xxyzz, tr_y_yyyyy_xxyzzz, tr_y_yyyyy_xxzzz, tr_y_yyyyy_xxzzzz, tr_y_yyyyy_xyyyy, tr_y_yyyyy_xyyyyy, tr_y_yyyyy_xyyyyz, tr_y_yyyyy_xyyyz, tr_y_yyyyy_xyyyzz, tr_y_yyyyy_xyyzz, tr_y_yyyyy_xyyzzz, tr_y_yyyyy_xyzzz, tr_y_yyyyy_xyzzzz, tr_y_yyyyy_xzzzz, tr_y_yyyyy_xzzzzz, tr_y_yyyyy_yyyyy, tr_y_yyyyy_yyyyyy, tr_y_yyyyy_yyyyyz, tr_y_yyyyy_yyyyz, tr_y_yyyyy_yyyyzz, tr_y_yyyyy_yyyzz, tr_y_yyyyy_yyyzzz, tr_y_yyyyy_yyzzz, tr_y_yyyyy_yyzzzz, tr_y_yyyyy_yzzzz, tr_y_yyyyy_yzzzzz, tr_y_yyyyy_zzzzz, tr_y_yyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyy_xxxxxx[i] = 6.0 * tr_y_yyyyy_xxxxx[i] * fe_0 + tr_y_yyyyy_xxxxxx[i] * pa_x[i];

        tr_y_xyyyyy_xxxxxy[i] = 5.0 * tr_y_yyyyy_xxxxy[i] * fe_0 + tr_y_yyyyy_xxxxxy[i] * pa_x[i];

        tr_y_xyyyyy_xxxxxz[i] = 5.0 * tr_y_yyyyy_xxxxz[i] * fe_0 + tr_y_yyyyy_xxxxxz[i] * pa_x[i];

        tr_y_xyyyyy_xxxxyy[i] = 4.0 * tr_y_yyyyy_xxxyy[i] * fe_0 + tr_y_yyyyy_xxxxyy[i] * pa_x[i];

        tr_y_xyyyyy_xxxxyz[i] = 4.0 * tr_y_yyyyy_xxxyz[i] * fe_0 + tr_y_yyyyy_xxxxyz[i] * pa_x[i];

        tr_y_xyyyyy_xxxxzz[i] = 4.0 * tr_y_yyyyy_xxxzz[i] * fe_0 + tr_y_yyyyy_xxxxzz[i] * pa_x[i];

        tr_y_xyyyyy_xxxyyy[i] = 3.0 * tr_y_yyyyy_xxyyy[i] * fe_0 + tr_y_yyyyy_xxxyyy[i] * pa_x[i];

        tr_y_xyyyyy_xxxyyz[i] = 3.0 * tr_y_yyyyy_xxyyz[i] * fe_0 + tr_y_yyyyy_xxxyyz[i] * pa_x[i];

        tr_y_xyyyyy_xxxyzz[i] = 3.0 * tr_y_yyyyy_xxyzz[i] * fe_0 + tr_y_yyyyy_xxxyzz[i] * pa_x[i];

        tr_y_xyyyyy_xxxzzz[i] = 3.0 * tr_y_yyyyy_xxzzz[i] * fe_0 + tr_y_yyyyy_xxxzzz[i] * pa_x[i];

        tr_y_xyyyyy_xxyyyy[i] = 2.0 * tr_y_yyyyy_xyyyy[i] * fe_0 + tr_y_yyyyy_xxyyyy[i] * pa_x[i];

        tr_y_xyyyyy_xxyyyz[i] = 2.0 * tr_y_yyyyy_xyyyz[i] * fe_0 + tr_y_yyyyy_xxyyyz[i] * pa_x[i];

        tr_y_xyyyyy_xxyyzz[i] = 2.0 * tr_y_yyyyy_xyyzz[i] * fe_0 + tr_y_yyyyy_xxyyzz[i] * pa_x[i];

        tr_y_xyyyyy_xxyzzz[i] = 2.0 * tr_y_yyyyy_xyzzz[i] * fe_0 + tr_y_yyyyy_xxyzzz[i] * pa_x[i];

        tr_y_xyyyyy_xxzzzz[i] = 2.0 * tr_y_yyyyy_xzzzz[i] * fe_0 + tr_y_yyyyy_xxzzzz[i] * pa_x[i];

        tr_y_xyyyyy_xyyyyy[i] = tr_y_yyyyy_yyyyy[i] * fe_0 + tr_y_yyyyy_xyyyyy[i] * pa_x[i];

        tr_y_xyyyyy_xyyyyz[i] = tr_y_yyyyy_yyyyz[i] * fe_0 + tr_y_yyyyy_xyyyyz[i] * pa_x[i];

        tr_y_xyyyyy_xyyyzz[i] = tr_y_yyyyy_yyyzz[i] * fe_0 + tr_y_yyyyy_xyyyzz[i] * pa_x[i];

        tr_y_xyyyyy_xyyzzz[i] = tr_y_yyyyy_yyzzz[i] * fe_0 + tr_y_yyyyy_xyyzzz[i] * pa_x[i];

        tr_y_xyyyyy_xyzzzz[i] = tr_y_yyyyy_yzzzz[i] * fe_0 + tr_y_yyyyy_xyzzzz[i] * pa_x[i];

        tr_y_xyyyyy_xzzzzz[i] = tr_y_yyyyy_zzzzz[i] * fe_0 + tr_y_yyyyy_xzzzzz[i] * pa_x[i];

        tr_y_xyyyyy_yyyyyy[i] = tr_y_yyyyy_yyyyyy[i] * pa_x[i];

        tr_y_xyyyyy_yyyyyz[i] = tr_y_yyyyy_yyyyyz[i] * pa_x[i];

        tr_y_xyyyyy_yyyyzz[i] = tr_y_yyyyy_yyyyzz[i] * pa_x[i];

        tr_y_xyyyyy_yyyzzz[i] = tr_y_yyyyy_yyyzzz[i] * pa_x[i];

        tr_y_xyyyyy_yyzzzz[i] = tr_y_yyyyy_yyzzzz[i] * pa_x[i];

        tr_y_xyyyyy_yzzzzz[i] = tr_y_yyyyy_yzzzzz[i] * pa_x[i];

        tr_y_xyyyyy_zzzzzz[i] = tr_y_yyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1232-1260 components of targeted buffer : II

    auto tr_y_xyyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 1232);

    auto tr_y_xyyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 1233);

    auto tr_y_xyyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 1234);

    auto tr_y_xyyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 1235);

    auto tr_y_xyyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 1236);

    auto tr_y_xyyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 1237);

    auto tr_y_xyyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 1238);

    auto tr_y_xyyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 1239);

    auto tr_y_xyyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 1240);

    auto tr_y_xyyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 1241);

    auto tr_y_xyyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 1242);

    auto tr_y_xyyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 1243);

    auto tr_y_xyyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 1244);

    auto tr_y_xyyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 1245);

    auto tr_y_xyyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 1246);

    auto tr_y_xyyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 1247);

    auto tr_y_xyyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 1248);

    auto tr_y_xyyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 1249);

    auto tr_y_xyyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 1250);

    auto tr_y_xyyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 1251);

    auto tr_y_xyyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 1252);

    auto tr_y_xyyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 1253);

    auto tr_y_xyyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 1254);

    auto tr_y_xyyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 1255);

    auto tr_y_xyyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 1256);

    auto tr_y_xyyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 1257);

    auto tr_y_xyyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 1258);

    auto tr_y_xyyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 1259);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xyyyy_xxxxxx, tr_y_xyyyy_xxxxxy, tr_y_xyyyy_xxxxyy, tr_y_xyyyy_xxxyyy, tr_y_xyyyy_xxyyyy, tr_y_xyyyy_xyyyyy, tr_y_xyyyyz_xxxxxx, tr_y_xyyyyz_xxxxxy, tr_y_xyyyyz_xxxxxz, tr_y_xyyyyz_xxxxyy, tr_y_xyyyyz_xxxxyz, tr_y_xyyyyz_xxxxzz, tr_y_xyyyyz_xxxyyy, tr_y_xyyyyz_xxxyyz, tr_y_xyyyyz_xxxyzz, tr_y_xyyyyz_xxxzzz, tr_y_xyyyyz_xxyyyy, tr_y_xyyyyz_xxyyyz, tr_y_xyyyyz_xxyyzz, tr_y_xyyyyz_xxyzzz, tr_y_xyyyyz_xxzzzz, tr_y_xyyyyz_xyyyyy, tr_y_xyyyyz_xyyyyz, tr_y_xyyyyz_xyyyzz, tr_y_xyyyyz_xyyzzz, tr_y_xyyyyz_xyzzzz, tr_y_xyyyyz_xzzzzz, tr_y_xyyyyz_yyyyyy, tr_y_xyyyyz_yyyyyz, tr_y_xyyyyz_yyyyzz, tr_y_xyyyyz_yyyzzz, tr_y_xyyyyz_yyzzzz, tr_y_xyyyyz_yzzzzz, tr_y_xyyyyz_zzzzzz, tr_y_yyyyz_xxxxxz, tr_y_yyyyz_xxxxyz, tr_y_yyyyz_xxxxz, tr_y_yyyyz_xxxxzz, tr_y_yyyyz_xxxyyz, tr_y_yyyyz_xxxyz, tr_y_yyyyz_xxxyzz, tr_y_yyyyz_xxxzz, tr_y_yyyyz_xxxzzz, tr_y_yyyyz_xxyyyz, tr_y_yyyyz_xxyyz, tr_y_yyyyz_xxyyzz, tr_y_yyyyz_xxyzz, tr_y_yyyyz_xxyzzz, tr_y_yyyyz_xxzzz, tr_y_yyyyz_xxzzzz, tr_y_yyyyz_xyyyyz, tr_y_yyyyz_xyyyz, tr_y_yyyyz_xyyyzz, tr_y_yyyyz_xyyzz, tr_y_yyyyz_xyyzzz, tr_y_yyyyz_xyzzz, tr_y_yyyyz_xyzzzz, tr_y_yyyyz_xzzzz, tr_y_yyyyz_xzzzzz, tr_y_yyyyz_yyyyyy, tr_y_yyyyz_yyyyyz, tr_y_yyyyz_yyyyz, tr_y_yyyyz_yyyyzz, tr_y_yyyyz_yyyzz, tr_y_yyyyz_yyyzzz, tr_y_yyyyz_yyzzz, tr_y_yyyyz_yyzzzz, tr_y_yyyyz_yzzzz, tr_y_yyyyz_yzzzzz, tr_y_yyyyz_zzzzz, tr_y_yyyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyz_xxxxxx[i] = tr_y_xyyyy_xxxxxx[i] * pa_z[i];

        tr_y_xyyyyz_xxxxxy[i] = tr_y_xyyyy_xxxxxy[i] * pa_z[i];

        tr_y_xyyyyz_xxxxxz[i] = 5.0 * tr_y_yyyyz_xxxxz[i] * fe_0 + tr_y_yyyyz_xxxxxz[i] * pa_x[i];

        tr_y_xyyyyz_xxxxyy[i] = tr_y_xyyyy_xxxxyy[i] * pa_z[i];

        tr_y_xyyyyz_xxxxyz[i] = 4.0 * tr_y_yyyyz_xxxyz[i] * fe_0 + tr_y_yyyyz_xxxxyz[i] * pa_x[i];

        tr_y_xyyyyz_xxxxzz[i] = 4.0 * tr_y_yyyyz_xxxzz[i] * fe_0 + tr_y_yyyyz_xxxxzz[i] * pa_x[i];

        tr_y_xyyyyz_xxxyyy[i] = tr_y_xyyyy_xxxyyy[i] * pa_z[i];

        tr_y_xyyyyz_xxxyyz[i] = 3.0 * tr_y_yyyyz_xxyyz[i] * fe_0 + tr_y_yyyyz_xxxyyz[i] * pa_x[i];

        tr_y_xyyyyz_xxxyzz[i] = 3.0 * tr_y_yyyyz_xxyzz[i] * fe_0 + tr_y_yyyyz_xxxyzz[i] * pa_x[i];

        tr_y_xyyyyz_xxxzzz[i] = 3.0 * tr_y_yyyyz_xxzzz[i] * fe_0 + tr_y_yyyyz_xxxzzz[i] * pa_x[i];

        tr_y_xyyyyz_xxyyyy[i] = tr_y_xyyyy_xxyyyy[i] * pa_z[i];

        tr_y_xyyyyz_xxyyyz[i] = 2.0 * tr_y_yyyyz_xyyyz[i] * fe_0 + tr_y_yyyyz_xxyyyz[i] * pa_x[i];

        tr_y_xyyyyz_xxyyzz[i] = 2.0 * tr_y_yyyyz_xyyzz[i] * fe_0 + tr_y_yyyyz_xxyyzz[i] * pa_x[i];

        tr_y_xyyyyz_xxyzzz[i] = 2.0 * tr_y_yyyyz_xyzzz[i] * fe_0 + tr_y_yyyyz_xxyzzz[i] * pa_x[i];

        tr_y_xyyyyz_xxzzzz[i] = 2.0 * tr_y_yyyyz_xzzzz[i] * fe_0 + tr_y_yyyyz_xxzzzz[i] * pa_x[i];

        tr_y_xyyyyz_xyyyyy[i] = tr_y_xyyyy_xyyyyy[i] * pa_z[i];

        tr_y_xyyyyz_xyyyyz[i] = tr_y_yyyyz_yyyyz[i] * fe_0 + tr_y_yyyyz_xyyyyz[i] * pa_x[i];

        tr_y_xyyyyz_xyyyzz[i] = tr_y_yyyyz_yyyzz[i] * fe_0 + tr_y_yyyyz_xyyyzz[i] * pa_x[i];

        tr_y_xyyyyz_xyyzzz[i] = tr_y_yyyyz_yyzzz[i] * fe_0 + tr_y_yyyyz_xyyzzz[i] * pa_x[i];

        tr_y_xyyyyz_xyzzzz[i] = tr_y_yyyyz_yzzzz[i] * fe_0 + tr_y_yyyyz_xyzzzz[i] * pa_x[i];

        tr_y_xyyyyz_xzzzzz[i] = tr_y_yyyyz_zzzzz[i] * fe_0 + tr_y_yyyyz_xzzzzz[i] * pa_x[i];

        tr_y_xyyyyz_yyyyyy[i] = tr_y_yyyyz_yyyyyy[i] * pa_x[i];

        tr_y_xyyyyz_yyyyyz[i] = tr_y_yyyyz_yyyyyz[i] * pa_x[i];

        tr_y_xyyyyz_yyyyzz[i] = tr_y_yyyyz_yyyyzz[i] * pa_x[i];

        tr_y_xyyyyz_yyyzzz[i] = tr_y_yyyyz_yyyzzz[i] * pa_x[i];

        tr_y_xyyyyz_yyzzzz[i] = tr_y_yyyyz_yyzzzz[i] * pa_x[i];

        tr_y_xyyyyz_yzzzzz[i] = tr_y_yyyyz_yzzzzz[i] * pa_x[i];

        tr_y_xyyyyz_zzzzzz[i] = tr_y_yyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1260-1288 components of targeted buffer : II

    auto tr_y_xyyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 1260);

    auto tr_y_xyyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 1261);

    auto tr_y_xyyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 1262);

    auto tr_y_xyyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 1263);

    auto tr_y_xyyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 1264);

    auto tr_y_xyyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 1265);

    auto tr_y_xyyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 1266);

    auto tr_y_xyyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 1267);

    auto tr_y_xyyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 1268);

    auto tr_y_xyyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 1269);

    auto tr_y_xyyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 1270);

    auto tr_y_xyyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 1271);

    auto tr_y_xyyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 1272);

    auto tr_y_xyyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 1273);

    auto tr_y_xyyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 1274);

    auto tr_y_xyyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 1275);

    auto tr_y_xyyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 1276);

    auto tr_y_xyyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 1277);

    auto tr_y_xyyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 1278);

    auto tr_y_xyyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 1279);

    auto tr_y_xyyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 1280);

    auto tr_y_xyyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 1281);

    auto tr_y_xyyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 1282);

    auto tr_y_xyyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 1283);

    auto tr_y_xyyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 1284);

    auto tr_y_xyyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 1285);

    auto tr_y_xyyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 1286);

    auto tr_y_xyyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 1287);

    #pragma omp simd aligned(pa_x, tr_y_xyyyzz_xxxxxx, tr_y_xyyyzz_xxxxxy, tr_y_xyyyzz_xxxxxz, tr_y_xyyyzz_xxxxyy, tr_y_xyyyzz_xxxxyz, tr_y_xyyyzz_xxxxzz, tr_y_xyyyzz_xxxyyy, tr_y_xyyyzz_xxxyyz, tr_y_xyyyzz_xxxyzz, tr_y_xyyyzz_xxxzzz, tr_y_xyyyzz_xxyyyy, tr_y_xyyyzz_xxyyyz, tr_y_xyyyzz_xxyyzz, tr_y_xyyyzz_xxyzzz, tr_y_xyyyzz_xxzzzz, tr_y_xyyyzz_xyyyyy, tr_y_xyyyzz_xyyyyz, tr_y_xyyyzz_xyyyzz, tr_y_xyyyzz_xyyzzz, tr_y_xyyyzz_xyzzzz, tr_y_xyyyzz_xzzzzz, tr_y_xyyyzz_yyyyyy, tr_y_xyyyzz_yyyyyz, tr_y_xyyyzz_yyyyzz, tr_y_xyyyzz_yyyzzz, tr_y_xyyyzz_yyzzzz, tr_y_xyyyzz_yzzzzz, tr_y_xyyyzz_zzzzzz, tr_y_yyyzz_xxxxx, tr_y_yyyzz_xxxxxx, tr_y_yyyzz_xxxxxy, tr_y_yyyzz_xxxxxz, tr_y_yyyzz_xxxxy, tr_y_yyyzz_xxxxyy, tr_y_yyyzz_xxxxyz, tr_y_yyyzz_xxxxz, tr_y_yyyzz_xxxxzz, tr_y_yyyzz_xxxyy, tr_y_yyyzz_xxxyyy, tr_y_yyyzz_xxxyyz, tr_y_yyyzz_xxxyz, tr_y_yyyzz_xxxyzz, tr_y_yyyzz_xxxzz, tr_y_yyyzz_xxxzzz, tr_y_yyyzz_xxyyy, tr_y_yyyzz_xxyyyy, tr_y_yyyzz_xxyyyz, tr_y_yyyzz_xxyyz, tr_y_yyyzz_xxyyzz, tr_y_yyyzz_xxyzz, tr_y_yyyzz_xxyzzz, tr_y_yyyzz_xxzzz, tr_y_yyyzz_xxzzzz, tr_y_yyyzz_xyyyy, tr_y_yyyzz_xyyyyy, tr_y_yyyzz_xyyyyz, tr_y_yyyzz_xyyyz, tr_y_yyyzz_xyyyzz, tr_y_yyyzz_xyyzz, tr_y_yyyzz_xyyzzz, tr_y_yyyzz_xyzzz, tr_y_yyyzz_xyzzzz, tr_y_yyyzz_xzzzz, tr_y_yyyzz_xzzzzz, tr_y_yyyzz_yyyyy, tr_y_yyyzz_yyyyyy, tr_y_yyyzz_yyyyyz, tr_y_yyyzz_yyyyz, tr_y_yyyzz_yyyyzz, tr_y_yyyzz_yyyzz, tr_y_yyyzz_yyyzzz, tr_y_yyyzz_yyzzz, tr_y_yyyzz_yyzzzz, tr_y_yyyzz_yzzzz, tr_y_yyyzz_yzzzzz, tr_y_yyyzz_zzzzz, tr_y_yyyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyzz_xxxxxx[i] = 6.0 * tr_y_yyyzz_xxxxx[i] * fe_0 + tr_y_yyyzz_xxxxxx[i] * pa_x[i];

        tr_y_xyyyzz_xxxxxy[i] = 5.0 * tr_y_yyyzz_xxxxy[i] * fe_0 + tr_y_yyyzz_xxxxxy[i] * pa_x[i];

        tr_y_xyyyzz_xxxxxz[i] = 5.0 * tr_y_yyyzz_xxxxz[i] * fe_0 + tr_y_yyyzz_xxxxxz[i] * pa_x[i];

        tr_y_xyyyzz_xxxxyy[i] = 4.0 * tr_y_yyyzz_xxxyy[i] * fe_0 + tr_y_yyyzz_xxxxyy[i] * pa_x[i];

        tr_y_xyyyzz_xxxxyz[i] = 4.0 * tr_y_yyyzz_xxxyz[i] * fe_0 + tr_y_yyyzz_xxxxyz[i] * pa_x[i];

        tr_y_xyyyzz_xxxxzz[i] = 4.0 * tr_y_yyyzz_xxxzz[i] * fe_0 + tr_y_yyyzz_xxxxzz[i] * pa_x[i];

        tr_y_xyyyzz_xxxyyy[i] = 3.0 * tr_y_yyyzz_xxyyy[i] * fe_0 + tr_y_yyyzz_xxxyyy[i] * pa_x[i];

        tr_y_xyyyzz_xxxyyz[i] = 3.0 * tr_y_yyyzz_xxyyz[i] * fe_0 + tr_y_yyyzz_xxxyyz[i] * pa_x[i];

        tr_y_xyyyzz_xxxyzz[i] = 3.0 * tr_y_yyyzz_xxyzz[i] * fe_0 + tr_y_yyyzz_xxxyzz[i] * pa_x[i];

        tr_y_xyyyzz_xxxzzz[i] = 3.0 * tr_y_yyyzz_xxzzz[i] * fe_0 + tr_y_yyyzz_xxxzzz[i] * pa_x[i];

        tr_y_xyyyzz_xxyyyy[i] = 2.0 * tr_y_yyyzz_xyyyy[i] * fe_0 + tr_y_yyyzz_xxyyyy[i] * pa_x[i];

        tr_y_xyyyzz_xxyyyz[i] = 2.0 * tr_y_yyyzz_xyyyz[i] * fe_0 + tr_y_yyyzz_xxyyyz[i] * pa_x[i];

        tr_y_xyyyzz_xxyyzz[i] = 2.0 * tr_y_yyyzz_xyyzz[i] * fe_0 + tr_y_yyyzz_xxyyzz[i] * pa_x[i];

        tr_y_xyyyzz_xxyzzz[i] = 2.0 * tr_y_yyyzz_xyzzz[i] * fe_0 + tr_y_yyyzz_xxyzzz[i] * pa_x[i];

        tr_y_xyyyzz_xxzzzz[i] = 2.0 * tr_y_yyyzz_xzzzz[i] * fe_0 + tr_y_yyyzz_xxzzzz[i] * pa_x[i];

        tr_y_xyyyzz_xyyyyy[i] = tr_y_yyyzz_yyyyy[i] * fe_0 + tr_y_yyyzz_xyyyyy[i] * pa_x[i];

        tr_y_xyyyzz_xyyyyz[i] = tr_y_yyyzz_yyyyz[i] * fe_0 + tr_y_yyyzz_xyyyyz[i] * pa_x[i];

        tr_y_xyyyzz_xyyyzz[i] = tr_y_yyyzz_yyyzz[i] * fe_0 + tr_y_yyyzz_xyyyzz[i] * pa_x[i];

        tr_y_xyyyzz_xyyzzz[i] = tr_y_yyyzz_yyzzz[i] * fe_0 + tr_y_yyyzz_xyyzzz[i] * pa_x[i];

        tr_y_xyyyzz_xyzzzz[i] = tr_y_yyyzz_yzzzz[i] * fe_0 + tr_y_yyyzz_xyzzzz[i] * pa_x[i];

        tr_y_xyyyzz_xzzzzz[i] = tr_y_yyyzz_zzzzz[i] * fe_0 + tr_y_yyyzz_xzzzzz[i] * pa_x[i];

        tr_y_xyyyzz_yyyyyy[i] = tr_y_yyyzz_yyyyyy[i] * pa_x[i];

        tr_y_xyyyzz_yyyyyz[i] = tr_y_yyyzz_yyyyyz[i] * pa_x[i];

        tr_y_xyyyzz_yyyyzz[i] = tr_y_yyyzz_yyyyzz[i] * pa_x[i];

        tr_y_xyyyzz_yyyzzz[i] = tr_y_yyyzz_yyyzzz[i] * pa_x[i];

        tr_y_xyyyzz_yyzzzz[i] = tr_y_yyyzz_yyzzzz[i] * pa_x[i];

        tr_y_xyyyzz_yzzzzz[i] = tr_y_yyyzz_yzzzzz[i] * pa_x[i];

        tr_y_xyyyzz_zzzzzz[i] = tr_y_yyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1288-1316 components of targeted buffer : II

    auto tr_y_xyyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1288);

    auto tr_y_xyyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1289);

    auto tr_y_xyyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1290);

    auto tr_y_xyyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1291);

    auto tr_y_xyyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1292);

    auto tr_y_xyyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1293);

    auto tr_y_xyyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1294);

    auto tr_y_xyyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1295);

    auto tr_y_xyyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1296);

    auto tr_y_xyyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1297);

    auto tr_y_xyyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1298);

    auto tr_y_xyyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1299);

    auto tr_y_xyyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1300);

    auto tr_y_xyyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1301);

    auto tr_y_xyyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1302);

    auto tr_y_xyyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1303);

    auto tr_y_xyyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1304);

    auto tr_y_xyyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1305);

    auto tr_y_xyyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1306);

    auto tr_y_xyyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1307);

    auto tr_y_xyyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1308);

    auto tr_y_xyyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1309);

    auto tr_y_xyyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1310);

    auto tr_y_xyyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1311);

    auto tr_y_xyyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1312);

    auto tr_y_xyyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1313);

    auto tr_y_xyyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1314);

    auto tr_y_xyyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1315);

    #pragma omp simd aligned(pa_x, tr_y_xyyzzz_xxxxxx, tr_y_xyyzzz_xxxxxy, tr_y_xyyzzz_xxxxxz, tr_y_xyyzzz_xxxxyy, tr_y_xyyzzz_xxxxyz, tr_y_xyyzzz_xxxxzz, tr_y_xyyzzz_xxxyyy, tr_y_xyyzzz_xxxyyz, tr_y_xyyzzz_xxxyzz, tr_y_xyyzzz_xxxzzz, tr_y_xyyzzz_xxyyyy, tr_y_xyyzzz_xxyyyz, tr_y_xyyzzz_xxyyzz, tr_y_xyyzzz_xxyzzz, tr_y_xyyzzz_xxzzzz, tr_y_xyyzzz_xyyyyy, tr_y_xyyzzz_xyyyyz, tr_y_xyyzzz_xyyyzz, tr_y_xyyzzz_xyyzzz, tr_y_xyyzzz_xyzzzz, tr_y_xyyzzz_xzzzzz, tr_y_xyyzzz_yyyyyy, tr_y_xyyzzz_yyyyyz, tr_y_xyyzzz_yyyyzz, tr_y_xyyzzz_yyyzzz, tr_y_xyyzzz_yyzzzz, tr_y_xyyzzz_yzzzzz, tr_y_xyyzzz_zzzzzz, tr_y_yyzzz_xxxxx, tr_y_yyzzz_xxxxxx, tr_y_yyzzz_xxxxxy, tr_y_yyzzz_xxxxxz, tr_y_yyzzz_xxxxy, tr_y_yyzzz_xxxxyy, tr_y_yyzzz_xxxxyz, tr_y_yyzzz_xxxxz, tr_y_yyzzz_xxxxzz, tr_y_yyzzz_xxxyy, tr_y_yyzzz_xxxyyy, tr_y_yyzzz_xxxyyz, tr_y_yyzzz_xxxyz, tr_y_yyzzz_xxxyzz, tr_y_yyzzz_xxxzz, tr_y_yyzzz_xxxzzz, tr_y_yyzzz_xxyyy, tr_y_yyzzz_xxyyyy, tr_y_yyzzz_xxyyyz, tr_y_yyzzz_xxyyz, tr_y_yyzzz_xxyyzz, tr_y_yyzzz_xxyzz, tr_y_yyzzz_xxyzzz, tr_y_yyzzz_xxzzz, tr_y_yyzzz_xxzzzz, tr_y_yyzzz_xyyyy, tr_y_yyzzz_xyyyyy, tr_y_yyzzz_xyyyyz, tr_y_yyzzz_xyyyz, tr_y_yyzzz_xyyyzz, tr_y_yyzzz_xyyzz, tr_y_yyzzz_xyyzzz, tr_y_yyzzz_xyzzz, tr_y_yyzzz_xyzzzz, tr_y_yyzzz_xzzzz, tr_y_yyzzz_xzzzzz, tr_y_yyzzz_yyyyy, tr_y_yyzzz_yyyyyy, tr_y_yyzzz_yyyyyz, tr_y_yyzzz_yyyyz, tr_y_yyzzz_yyyyzz, tr_y_yyzzz_yyyzz, tr_y_yyzzz_yyyzzz, tr_y_yyzzz_yyzzz, tr_y_yyzzz_yyzzzz, tr_y_yyzzz_yzzzz, tr_y_yyzzz_yzzzzz, tr_y_yyzzz_zzzzz, tr_y_yyzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzzz_xxxxxx[i] = 6.0 * tr_y_yyzzz_xxxxx[i] * fe_0 + tr_y_yyzzz_xxxxxx[i] * pa_x[i];

        tr_y_xyyzzz_xxxxxy[i] = 5.0 * tr_y_yyzzz_xxxxy[i] * fe_0 + tr_y_yyzzz_xxxxxy[i] * pa_x[i];

        tr_y_xyyzzz_xxxxxz[i] = 5.0 * tr_y_yyzzz_xxxxz[i] * fe_0 + tr_y_yyzzz_xxxxxz[i] * pa_x[i];

        tr_y_xyyzzz_xxxxyy[i] = 4.0 * tr_y_yyzzz_xxxyy[i] * fe_0 + tr_y_yyzzz_xxxxyy[i] * pa_x[i];

        tr_y_xyyzzz_xxxxyz[i] = 4.0 * tr_y_yyzzz_xxxyz[i] * fe_0 + tr_y_yyzzz_xxxxyz[i] * pa_x[i];

        tr_y_xyyzzz_xxxxzz[i] = 4.0 * tr_y_yyzzz_xxxzz[i] * fe_0 + tr_y_yyzzz_xxxxzz[i] * pa_x[i];

        tr_y_xyyzzz_xxxyyy[i] = 3.0 * tr_y_yyzzz_xxyyy[i] * fe_0 + tr_y_yyzzz_xxxyyy[i] * pa_x[i];

        tr_y_xyyzzz_xxxyyz[i] = 3.0 * tr_y_yyzzz_xxyyz[i] * fe_0 + tr_y_yyzzz_xxxyyz[i] * pa_x[i];

        tr_y_xyyzzz_xxxyzz[i] = 3.0 * tr_y_yyzzz_xxyzz[i] * fe_0 + tr_y_yyzzz_xxxyzz[i] * pa_x[i];

        tr_y_xyyzzz_xxxzzz[i] = 3.0 * tr_y_yyzzz_xxzzz[i] * fe_0 + tr_y_yyzzz_xxxzzz[i] * pa_x[i];

        tr_y_xyyzzz_xxyyyy[i] = 2.0 * tr_y_yyzzz_xyyyy[i] * fe_0 + tr_y_yyzzz_xxyyyy[i] * pa_x[i];

        tr_y_xyyzzz_xxyyyz[i] = 2.0 * tr_y_yyzzz_xyyyz[i] * fe_0 + tr_y_yyzzz_xxyyyz[i] * pa_x[i];

        tr_y_xyyzzz_xxyyzz[i] = 2.0 * tr_y_yyzzz_xyyzz[i] * fe_0 + tr_y_yyzzz_xxyyzz[i] * pa_x[i];

        tr_y_xyyzzz_xxyzzz[i] = 2.0 * tr_y_yyzzz_xyzzz[i] * fe_0 + tr_y_yyzzz_xxyzzz[i] * pa_x[i];

        tr_y_xyyzzz_xxzzzz[i] = 2.0 * tr_y_yyzzz_xzzzz[i] * fe_0 + tr_y_yyzzz_xxzzzz[i] * pa_x[i];

        tr_y_xyyzzz_xyyyyy[i] = tr_y_yyzzz_yyyyy[i] * fe_0 + tr_y_yyzzz_xyyyyy[i] * pa_x[i];

        tr_y_xyyzzz_xyyyyz[i] = tr_y_yyzzz_yyyyz[i] * fe_0 + tr_y_yyzzz_xyyyyz[i] * pa_x[i];

        tr_y_xyyzzz_xyyyzz[i] = tr_y_yyzzz_yyyzz[i] * fe_0 + tr_y_yyzzz_xyyyzz[i] * pa_x[i];

        tr_y_xyyzzz_xyyzzz[i] = tr_y_yyzzz_yyzzz[i] * fe_0 + tr_y_yyzzz_xyyzzz[i] * pa_x[i];

        tr_y_xyyzzz_xyzzzz[i] = tr_y_yyzzz_yzzzz[i] * fe_0 + tr_y_yyzzz_xyzzzz[i] * pa_x[i];

        tr_y_xyyzzz_xzzzzz[i] = tr_y_yyzzz_zzzzz[i] * fe_0 + tr_y_yyzzz_xzzzzz[i] * pa_x[i];

        tr_y_xyyzzz_yyyyyy[i] = tr_y_yyzzz_yyyyyy[i] * pa_x[i];

        tr_y_xyyzzz_yyyyyz[i] = tr_y_yyzzz_yyyyyz[i] * pa_x[i];

        tr_y_xyyzzz_yyyyzz[i] = tr_y_yyzzz_yyyyzz[i] * pa_x[i];

        tr_y_xyyzzz_yyyzzz[i] = tr_y_yyzzz_yyyzzz[i] * pa_x[i];

        tr_y_xyyzzz_yyzzzz[i] = tr_y_yyzzz_yyzzzz[i] * pa_x[i];

        tr_y_xyyzzz_yzzzzz[i] = tr_y_yyzzz_yzzzzz[i] * pa_x[i];

        tr_y_xyyzzz_zzzzzz[i] = tr_y_yyzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1316-1344 components of targeted buffer : II

    auto tr_y_xyzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1316);

    auto tr_y_xyzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1317);

    auto tr_y_xyzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1318);

    auto tr_y_xyzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1319);

    auto tr_y_xyzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1320);

    auto tr_y_xyzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1321);

    auto tr_y_xyzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1322);

    auto tr_y_xyzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1323);

    auto tr_y_xyzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1324);

    auto tr_y_xyzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1325);

    auto tr_y_xyzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1326);

    auto tr_y_xyzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1327);

    auto tr_y_xyzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1328);

    auto tr_y_xyzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1329);

    auto tr_y_xyzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1330);

    auto tr_y_xyzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1331);

    auto tr_y_xyzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1332);

    auto tr_y_xyzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1333);

    auto tr_y_xyzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1334);

    auto tr_y_xyzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1335);

    auto tr_y_xyzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1336);

    auto tr_y_xyzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1337);

    auto tr_y_xyzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1338);

    auto tr_y_xyzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1339);

    auto tr_y_xyzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1340);

    auto tr_y_xyzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1341);

    auto tr_y_xyzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1342);

    auto tr_y_xyzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1343);

    #pragma omp simd aligned(pa_x, tr_y_xyzzzz_xxxxxx, tr_y_xyzzzz_xxxxxy, tr_y_xyzzzz_xxxxxz, tr_y_xyzzzz_xxxxyy, tr_y_xyzzzz_xxxxyz, tr_y_xyzzzz_xxxxzz, tr_y_xyzzzz_xxxyyy, tr_y_xyzzzz_xxxyyz, tr_y_xyzzzz_xxxyzz, tr_y_xyzzzz_xxxzzz, tr_y_xyzzzz_xxyyyy, tr_y_xyzzzz_xxyyyz, tr_y_xyzzzz_xxyyzz, tr_y_xyzzzz_xxyzzz, tr_y_xyzzzz_xxzzzz, tr_y_xyzzzz_xyyyyy, tr_y_xyzzzz_xyyyyz, tr_y_xyzzzz_xyyyzz, tr_y_xyzzzz_xyyzzz, tr_y_xyzzzz_xyzzzz, tr_y_xyzzzz_xzzzzz, tr_y_xyzzzz_yyyyyy, tr_y_xyzzzz_yyyyyz, tr_y_xyzzzz_yyyyzz, tr_y_xyzzzz_yyyzzz, tr_y_xyzzzz_yyzzzz, tr_y_xyzzzz_yzzzzz, tr_y_xyzzzz_zzzzzz, tr_y_yzzzz_xxxxx, tr_y_yzzzz_xxxxxx, tr_y_yzzzz_xxxxxy, tr_y_yzzzz_xxxxxz, tr_y_yzzzz_xxxxy, tr_y_yzzzz_xxxxyy, tr_y_yzzzz_xxxxyz, tr_y_yzzzz_xxxxz, tr_y_yzzzz_xxxxzz, tr_y_yzzzz_xxxyy, tr_y_yzzzz_xxxyyy, tr_y_yzzzz_xxxyyz, tr_y_yzzzz_xxxyz, tr_y_yzzzz_xxxyzz, tr_y_yzzzz_xxxzz, tr_y_yzzzz_xxxzzz, tr_y_yzzzz_xxyyy, tr_y_yzzzz_xxyyyy, tr_y_yzzzz_xxyyyz, tr_y_yzzzz_xxyyz, tr_y_yzzzz_xxyyzz, tr_y_yzzzz_xxyzz, tr_y_yzzzz_xxyzzz, tr_y_yzzzz_xxzzz, tr_y_yzzzz_xxzzzz, tr_y_yzzzz_xyyyy, tr_y_yzzzz_xyyyyy, tr_y_yzzzz_xyyyyz, tr_y_yzzzz_xyyyz, tr_y_yzzzz_xyyyzz, tr_y_yzzzz_xyyzz, tr_y_yzzzz_xyyzzz, tr_y_yzzzz_xyzzz, tr_y_yzzzz_xyzzzz, tr_y_yzzzz_xzzzz, tr_y_yzzzz_xzzzzz, tr_y_yzzzz_yyyyy, tr_y_yzzzz_yyyyyy, tr_y_yzzzz_yyyyyz, tr_y_yzzzz_yyyyz, tr_y_yzzzz_yyyyzz, tr_y_yzzzz_yyyzz, tr_y_yzzzz_yyyzzz, tr_y_yzzzz_yyzzz, tr_y_yzzzz_yyzzzz, tr_y_yzzzz_yzzzz, tr_y_yzzzz_yzzzzz, tr_y_yzzzz_zzzzz, tr_y_yzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzzz_xxxxxx[i] = 6.0 * tr_y_yzzzz_xxxxx[i] * fe_0 + tr_y_yzzzz_xxxxxx[i] * pa_x[i];

        tr_y_xyzzzz_xxxxxy[i] = 5.0 * tr_y_yzzzz_xxxxy[i] * fe_0 + tr_y_yzzzz_xxxxxy[i] * pa_x[i];

        tr_y_xyzzzz_xxxxxz[i] = 5.0 * tr_y_yzzzz_xxxxz[i] * fe_0 + tr_y_yzzzz_xxxxxz[i] * pa_x[i];

        tr_y_xyzzzz_xxxxyy[i] = 4.0 * tr_y_yzzzz_xxxyy[i] * fe_0 + tr_y_yzzzz_xxxxyy[i] * pa_x[i];

        tr_y_xyzzzz_xxxxyz[i] = 4.0 * tr_y_yzzzz_xxxyz[i] * fe_0 + tr_y_yzzzz_xxxxyz[i] * pa_x[i];

        tr_y_xyzzzz_xxxxzz[i] = 4.0 * tr_y_yzzzz_xxxzz[i] * fe_0 + tr_y_yzzzz_xxxxzz[i] * pa_x[i];

        tr_y_xyzzzz_xxxyyy[i] = 3.0 * tr_y_yzzzz_xxyyy[i] * fe_0 + tr_y_yzzzz_xxxyyy[i] * pa_x[i];

        tr_y_xyzzzz_xxxyyz[i] = 3.0 * tr_y_yzzzz_xxyyz[i] * fe_0 + tr_y_yzzzz_xxxyyz[i] * pa_x[i];

        tr_y_xyzzzz_xxxyzz[i] = 3.0 * tr_y_yzzzz_xxyzz[i] * fe_0 + tr_y_yzzzz_xxxyzz[i] * pa_x[i];

        tr_y_xyzzzz_xxxzzz[i] = 3.0 * tr_y_yzzzz_xxzzz[i] * fe_0 + tr_y_yzzzz_xxxzzz[i] * pa_x[i];

        tr_y_xyzzzz_xxyyyy[i] = 2.0 * tr_y_yzzzz_xyyyy[i] * fe_0 + tr_y_yzzzz_xxyyyy[i] * pa_x[i];

        tr_y_xyzzzz_xxyyyz[i] = 2.0 * tr_y_yzzzz_xyyyz[i] * fe_0 + tr_y_yzzzz_xxyyyz[i] * pa_x[i];

        tr_y_xyzzzz_xxyyzz[i] = 2.0 * tr_y_yzzzz_xyyzz[i] * fe_0 + tr_y_yzzzz_xxyyzz[i] * pa_x[i];

        tr_y_xyzzzz_xxyzzz[i] = 2.0 * tr_y_yzzzz_xyzzz[i] * fe_0 + tr_y_yzzzz_xxyzzz[i] * pa_x[i];

        tr_y_xyzzzz_xxzzzz[i] = 2.0 * tr_y_yzzzz_xzzzz[i] * fe_0 + tr_y_yzzzz_xxzzzz[i] * pa_x[i];

        tr_y_xyzzzz_xyyyyy[i] = tr_y_yzzzz_yyyyy[i] * fe_0 + tr_y_yzzzz_xyyyyy[i] * pa_x[i];

        tr_y_xyzzzz_xyyyyz[i] = tr_y_yzzzz_yyyyz[i] * fe_0 + tr_y_yzzzz_xyyyyz[i] * pa_x[i];

        tr_y_xyzzzz_xyyyzz[i] = tr_y_yzzzz_yyyzz[i] * fe_0 + tr_y_yzzzz_xyyyzz[i] * pa_x[i];

        tr_y_xyzzzz_xyyzzz[i] = tr_y_yzzzz_yyzzz[i] * fe_0 + tr_y_yzzzz_xyyzzz[i] * pa_x[i];

        tr_y_xyzzzz_xyzzzz[i] = tr_y_yzzzz_yzzzz[i] * fe_0 + tr_y_yzzzz_xyzzzz[i] * pa_x[i];

        tr_y_xyzzzz_xzzzzz[i] = tr_y_yzzzz_zzzzz[i] * fe_0 + tr_y_yzzzz_xzzzzz[i] * pa_x[i];

        tr_y_xyzzzz_yyyyyy[i] = tr_y_yzzzz_yyyyyy[i] * pa_x[i];

        tr_y_xyzzzz_yyyyyz[i] = tr_y_yzzzz_yyyyyz[i] * pa_x[i];

        tr_y_xyzzzz_yyyyzz[i] = tr_y_yzzzz_yyyyzz[i] * pa_x[i];

        tr_y_xyzzzz_yyyzzz[i] = tr_y_yzzzz_yyyzzz[i] * pa_x[i];

        tr_y_xyzzzz_yyzzzz[i] = tr_y_yzzzz_yyzzzz[i] * pa_x[i];

        tr_y_xyzzzz_yzzzzz[i] = tr_y_yzzzz_yzzzzz[i] * pa_x[i];

        tr_y_xyzzzz_zzzzzz[i] = tr_y_yzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1344-1372 components of targeted buffer : II

    auto tr_y_xzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1344);

    auto tr_y_xzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1345);

    auto tr_y_xzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1346);

    auto tr_y_xzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1347);

    auto tr_y_xzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1348);

    auto tr_y_xzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1349);

    auto tr_y_xzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1350);

    auto tr_y_xzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1351);

    auto tr_y_xzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1352);

    auto tr_y_xzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1353);

    auto tr_y_xzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1354);

    auto tr_y_xzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1355);

    auto tr_y_xzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1356);

    auto tr_y_xzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1357);

    auto tr_y_xzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1358);

    auto tr_y_xzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1359);

    auto tr_y_xzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1360);

    auto tr_y_xzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1361);

    auto tr_y_xzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1362);

    auto tr_y_xzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1363);

    auto tr_y_xzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1364);

    auto tr_y_xzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1365);

    auto tr_y_xzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1366);

    auto tr_y_xzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1367);

    auto tr_y_xzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1368);

    auto tr_y_xzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1369);

    auto tr_y_xzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1370);

    auto tr_y_xzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1371);

    #pragma omp simd aligned(pa_x, tr_y_xzzzzz_xxxxxx, tr_y_xzzzzz_xxxxxy, tr_y_xzzzzz_xxxxxz, tr_y_xzzzzz_xxxxyy, tr_y_xzzzzz_xxxxyz, tr_y_xzzzzz_xxxxzz, tr_y_xzzzzz_xxxyyy, tr_y_xzzzzz_xxxyyz, tr_y_xzzzzz_xxxyzz, tr_y_xzzzzz_xxxzzz, tr_y_xzzzzz_xxyyyy, tr_y_xzzzzz_xxyyyz, tr_y_xzzzzz_xxyyzz, tr_y_xzzzzz_xxyzzz, tr_y_xzzzzz_xxzzzz, tr_y_xzzzzz_xyyyyy, tr_y_xzzzzz_xyyyyz, tr_y_xzzzzz_xyyyzz, tr_y_xzzzzz_xyyzzz, tr_y_xzzzzz_xyzzzz, tr_y_xzzzzz_xzzzzz, tr_y_xzzzzz_yyyyyy, tr_y_xzzzzz_yyyyyz, tr_y_xzzzzz_yyyyzz, tr_y_xzzzzz_yyyzzz, tr_y_xzzzzz_yyzzzz, tr_y_xzzzzz_yzzzzz, tr_y_xzzzzz_zzzzzz, tr_y_zzzzz_xxxxx, tr_y_zzzzz_xxxxxx, tr_y_zzzzz_xxxxxy, tr_y_zzzzz_xxxxxz, tr_y_zzzzz_xxxxy, tr_y_zzzzz_xxxxyy, tr_y_zzzzz_xxxxyz, tr_y_zzzzz_xxxxz, tr_y_zzzzz_xxxxzz, tr_y_zzzzz_xxxyy, tr_y_zzzzz_xxxyyy, tr_y_zzzzz_xxxyyz, tr_y_zzzzz_xxxyz, tr_y_zzzzz_xxxyzz, tr_y_zzzzz_xxxzz, tr_y_zzzzz_xxxzzz, tr_y_zzzzz_xxyyy, tr_y_zzzzz_xxyyyy, tr_y_zzzzz_xxyyyz, tr_y_zzzzz_xxyyz, tr_y_zzzzz_xxyyzz, tr_y_zzzzz_xxyzz, tr_y_zzzzz_xxyzzz, tr_y_zzzzz_xxzzz, tr_y_zzzzz_xxzzzz, tr_y_zzzzz_xyyyy, tr_y_zzzzz_xyyyyy, tr_y_zzzzz_xyyyyz, tr_y_zzzzz_xyyyz, tr_y_zzzzz_xyyyzz, tr_y_zzzzz_xyyzz, tr_y_zzzzz_xyyzzz, tr_y_zzzzz_xyzzz, tr_y_zzzzz_xyzzzz, tr_y_zzzzz_xzzzz, tr_y_zzzzz_xzzzzz, tr_y_zzzzz_yyyyy, tr_y_zzzzz_yyyyyy, tr_y_zzzzz_yyyyyz, tr_y_zzzzz_yyyyz, tr_y_zzzzz_yyyyzz, tr_y_zzzzz_yyyzz, tr_y_zzzzz_yyyzzz, tr_y_zzzzz_yyzzz, tr_y_zzzzz_yyzzzz, tr_y_zzzzz_yzzzz, tr_y_zzzzz_yzzzzz, tr_y_zzzzz_zzzzz, tr_y_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzzz_xxxxxx[i] = 6.0 * tr_y_zzzzz_xxxxx[i] * fe_0 + tr_y_zzzzz_xxxxxx[i] * pa_x[i];

        tr_y_xzzzzz_xxxxxy[i] = 5.0 * tr_y_zzzzz_xxxxy[i] * fe_0 + tr_y_zzzzz_xxxxxy[i] * pa_x[i];

        tr_y_xzzzzz_xxxxxz[i] = 5.0 * tr_y_zzzzz_xxxxz[i] * fe_0 + tr_y_zzzzz_xxxxxz[i] * pa_x[i];

        tr_y_xzzzzz_xxxxyy[i] = 4.0 * tr_y_zzzzz_xxxyy[i] * fe_0 + tr_y_zzzzz_xxxxyy[i] * pa_x[i];

        tr_y_xzzzzz_xxxxyz[i] = 4.0 * tr_y_zzzzz_xxxyz[i] * fe_0 + tr_y_zzzzz_xxxxyz[i] * pa_x[i];

        tr_y_xzzzzz_xxxxzz[i] = 4.0 * tr_y_zzzzz_xxxzz[i] * fe_0 + tr_y_zzzzz_xxxxzz[i] * pa_x[i];

        tr_y_xzzzzz_xxxyyy[i] = 3.0 * tr_y_zzzzz_xxyyy[i] * fe_0 + tr_y_zzzzz_xxxyyy[i] * pa_x[i];

        tr_y_xzzzzz_xxxyyz[i] = 3.0 * tr_y_zzzzz_xxyyz[i] * fe_0 + tr_y_zzzzz_xxxyyz[i] * pa_x[i];

        tr_y_xzzzzz_xxxyzz[i] = 3.0 * tr_y_zzzzz_xxyzz[i] * fe_0 + tr_y_zzzzz_xxxyzz[i] * pa_x[i];

        tr_y_xzzzzz_xxxzzz[i] = 3.0 * tr_y_zzzzz_xxzzz[i] * fe_0 + tr_y_zzzzz_xxxzzz[i] * pa_x[i];

        tr_y_xzzzzz_xxyyyy[i] = 2.0 * tr_y_zzzzz_xyyyy[i] * fe_0 + tr_y_zzzzz_xxyyyy[i] * pa_x[i];

        tr_y_xzzzzz_xxyyyz[i] = 2.0 * tr_y_zzzzz_xyyyz[i] * fe_0 + tr_y_zzzzz_xxyyyz[i] * pa_x[i];

        tr_y_xzzzzz_xxyyzz[i] = 2.0 * tr_y_zzzzz_xyyzz[i] * fe_0 + tr_y_zzzzz_xxyyzz[i] * pa_x[i];

        tr_y_xzzzzz_xxyzzz[i] = 2.0 * tr_y_zzzzz_xyzzz[i] * fe_0 + tr_y_zzzzz_xxyzzz[i] * pa_x[i];

        tr_y_xzzzzz_xxzzzz[i] = 2.0 * tr_y_zzzzz_xzzzz[i] * fe_0 + tr_y_zzzzz_xxzzzz[i] * pa_x[i];

        tr_y_xzzzzz_xyyyyy[i] = tr_y_zzzzz_yyyyy[i] * fe_0 + tr_y_zzzzz_xyyyyy[i] * pa_x[i];

        tr_y_xzzzzz_xyyyyz[i] = tr_y_zzzzz_yyyyz[i] * fe_0 + tr_y_zzzzz_xyyyyz[i] * pa_x[i];

        tr_y_xzzzzz_xyyyzz[i] = tr_y_zzzzz_yyyzz[i] * fe_0 + tr_y_zzzzz_xyyyzz[i] * pa_x[i];

        tr_y_xzzzzz_xyyzzz[i] = tr_y_zzzzz_yyzzz[i] * fe_0 + tr_y_zzzzz_xyyzzz[i] * pa_x[i];

        tr_y_xzzzzz_xyzzzz[i] = tr_y_zzzzz_yzzzz[i] * fe_0 + tr_y_zzzzz_xyzzzz[i] * pa_x[i];

        tr_y_xzzzzz_xzzzzz[i] = tr_y_zzzzz_zzzzz[i] * fe_0 + tr_y_zzzzz_xzzzzz[i] * pa_x[i];

        tr_y_xzzzzz_yyyyyy[i] = tr_y_zzzzz_yyyyyy[i] * pa_x[i];

        tr_y_xzzzzz_yyyyyz[i] = tr_y_zzzzz_yyyyyz[i] * pa_x[i];

        tr_y_xzzzzz_yyyyzz[i] = tr_y_zzzzz_yyyyzz[i] * pa_x[i];

        tr_y_xzzzzz_yyyzzz[i] = tr_y_zzzzz_yyyzzz[i] * pa_x[i];

        tr_y_xzzzzz_yyzzzz[i] = tr_y_zzzzz_yyzzzz[i] * pa_x[i];

        tr_y_xzzzzz_yzzzzz[i] = tr_y_zzzzz_yzzzzz[i] * pa_x[i];

        tr_y_xzzzzz_zzzzzz[i] = tr_y_zzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1372-1400 components of targeted buffer : II

    auto tr_y_yyyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 1372);

    auto tr_y_yyyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 1373);

    auto tr_y_yyyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 1374);

    auto tr_y_yyyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 1375);

    auto tr_y_yyyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 1376);

    auto tr_y_yyyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 1377);

    auto tr_y_yyyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 1378);

    auto tr_y_yyyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 1379);

    auto tr_y_yyyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 1380);

    auto tr_y_yyyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 1381);

    auto tr_y_yyyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 1382);

    auto tr_y_yyyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 1383);

    auto tr_y_yyyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 1384);

    auto tr_y_yyyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 1385);

    auto tr_y_yyyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 1386);

    auto tr_y_yyyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 1387);

    auto tr_y_yyyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 1388);

    auto tr_y_yyyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 1389);

    auto tr_y_yyyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 1390);

    auto tr_y_yyyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 1391);

    auto tr_y_yyyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 1392);

    auto tr_y_yyyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 1393);

    auto tr_y_yyyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 1394);

    auto tr_y_yyyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 1395);

    auto tr_y_yyyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 1396);

    auto tr_y_yyyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 1397);

    auto tr_y_yyyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 1398);

    auto tr_y_yyyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 1399);

    #pragma omp simd aligned(pa_y, tr_y_yyyy_xxxxxx, tr_y_yyyy_xxxxxy, tr_y_yyyy_xxxxxz, tr_y_yyyy_xxxxyy, tr_y_yyyy_xxxxyz, tr_y_yyyy_xxxxzz, tr_y_yyyy_xxxyyy, tr_y_yyyy_xxxyyz, tr_y_yyyy_xxxyzz, tr_y_yyyy_xxxzzz, tr_y_yyyy_xxyyyy, tr_y_yyyy_xxyyyz, tr_y_yyyy_xxyyzz, tr_y_yyyy_xxyzzz, tr_y_yyyy_xxzzzz, tr_y_yyyy_xyyyyy, tr_y_yyyy_xyyyyz, tr_y_yyyy_xyyyzz, tr_y_yyyy_xyyzzz, tr_y_yyyy_xyzzzz, tr_y_yyyy_xzzzzz, tr_y_yyyy_yyyyyy, tr_y_yyyy_yyyyyz, tr_y_yyyy_yyyyzz, tr_y_yyyy_yyyzzz, tr_y_yyyy_yyzzzz, tr_y_yyyy_yzzzzz, tr_y_yyyy_zzzzzz, tr_y_yyyyy_xxxxx, tr_y_yyyyy_xxxxxx, tr_y_yyyyy_xxxxxy, tr_y_yyyyy_xxxxxz, tr_y_yyyyy_xxxxy, tr_y_yyyyy_xxxxyy, tr_y_yyyyy_xxxxyz, tr_y_yyyyy_xxxxz, tr_y_yyyyy_xxxxzz, tr_y_yyyyy_xxxyy, tr_y_yyyyy_xxxyyy, tr_y_yyyyy_xxxyyz, tr_y_yyyyy_xxxyz, tr_y_yyyyy_xxxyzz, tr_y_yyyyy_xxxzz, tr_y_yyyyy_xxxzzz, tr_y_yyyyy_xxyyy, tr_y_yyyyy_xxyyyy, tr_y_yyyyy_xxyyyz, tr_y_yyyyy_xxyyz, tr_y_yyyyy_xxyyzz, tr_y_yyyyy_xxyzz, tr_y_yyyyy_xxyzzz, tr_y_yyyyy_xxzzz, tr_y_yyyyy_xxzzzz, tr_y_yyyyy_xyyyy, tr_y_yyyyy_xyyyyy, tr_y_yyyyy_xyyyyz, tr_y_yyyyy_xyyyz, tr_y_yyyyy_xyyyzz, tr_y_yyyyy_xyyzz, tr_y_yyyyy_xyyzzz, tr_y_yyyyy_xyzzz, tr_y_yyyyy_xyzzzz, tr_y_yyyyy_xzzzz, tr_y_yyyyy_xzzzzz, tr_y_yyyyy_yyyyy, tr_y_yyyyy_yyyyyy, tr_y_yyyyy_yyyyyz, tr_y_yyyyy_yyyyz, tr_y_yyyyy_yyyyzz, tr_y_yyyyy_yyyzz, tr_y_yyyyy_yyyzzz, tr_y_yyyyy_yyzzz, tr_y_yyyyy_yyzzzz, tr_y_yyyyy_yzzzz, tr_y_yyyyy_yzzzzz, tr_y_yyyyy_zzzzz, tr_y_yyyyy_zzzzzz, tr_y_yyyyyy_xxxxxx, tr_y_yyyyyy_xxxxxy, tr_y_yyyyyy_xxxxxz, tr_y_yyyyyy_xxxxyy, tr_y_yyyyyy_xxxxyz, tr_y_yyyyyy_xxxxzz, tr_y_yyyyyy_xxxyyy, tr_y_yyyyyy_xxxyyz, tr_y_yyyyyy_xxxyzz, tr_y_yyyyyy_xxxzzz, tr_y_yyyyyy_xxyyyy, tr_y_yyyyyy_xxyyyz, tr_y_yyyyyy_xxyyzz, tr_y_yyyyyy_xxyzzz, tr_y_yyyyyy_xxzzzz, tr_y_yyyyyy_xyyyyy, tr_y_yyyyyy_xyyyyz, tr_y_yyyyyy_xyyyzz, tr_y_yyyyyy_xyyzzz, tr_y_yyyyyy_xyzzzz, tr_y_yyyyyy_xzzzzz, tr_y_yyyyyy_yyyyyy, tr_y_yyyyyy_yyyyyz, tr_y_yyyyyy_yyyyzz, tr_y_yyyyyy_yyyzzz, tr_y_yyyyyy_yyzzzz, tr_y_yyyyyy_yzzzzz, tr_y_yyyyyy_zzzzzz, ts_yyyyy_xxxxxx, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxxz, ts_yyyyy_xxxxyy, ts_yyyyy_xxxxyz, ts_yyyyy_xxxxzz, ts_yyyyy_xxxyyy, ts_yyyyy_xxxyyz, ts_yyyyy_xxxyzz, ts_yyyyy_xxxzzz, ts_yyyyy_xxyyyy, ts_yyyyy_xxyyyz, ts_yyyyy_xxyyzz, ts_yyyyy_xxyzzz, ts_yyyyy_xxzzzz, ts_yyyyy_xyyyyy, ts_yyyyy_xyyyyz, ts_yyyyy_xyyyzz, ts_yyyyy_xyyzzz, ts_yyyyy_xyzzzz, ts_yyyyy_xzzzzz, ts_yyyyy_yyyyyy, ts_yyyyy_yyyyyz, ts_yyyyy_yyyyzz, ts_yyyyy_yyyzzz, ts_yyyyy_yyzzzz, ts_yyyyy_yzzzzz, ts_yyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyy_xxxxxx[i] = 5.0 * tr_y_yyyy_xxxxxx[i] * fe_0 + ts_yyyyy_xxxxxx[i] * fe_0 + tr_y_yyyyy_xxxxxx[i] * pa_y[i];

        tr_y_yyyyyy_xxxxxy[i] = 5.0 * tr_y_yyyy_xxxxxy[i] * fe_0 + tr_y_yyyyy_xxxxx[i] * fe_0 + ts_yyyyy_xxxxxy[i] * fe_0 + tr_y_yyyyy_xxxxxy[i] * pa_y[i];

        tr_y_yyyyyy_xxxxxz[i] = 5.0 * tr_y_yyyy_xxxxxz[i] * fe_0 + ts_yyyyy_xxxxxz[i] * fe_0 + tr_y_yyyyy_xxxxxz[i] * pa_y[i];

        tr_y_yyyyyy_xxxxyy[i] = 5.0 * tr_y_yyyy_xxxxyy[i] * fe_0 + 2.0 * tr_y_yyyyy_xxxxy[i] * fe_0 + ts_yyyyy_xxxxyy[i] * fe_0 + tr_y_yyyyy_xxxxyy[i] * pa_y[i];

        tr_y_yyyyyy_xxxxyz[i] = 5.0 * tr_y_yyyy_xxxxyz[i] * fe_0 + tr_y_yyyyy_xxxxz[i] * fe_0 + ts_yyyyy_xxxxyz[i] * fe_0 + tr_y_yyyyy_xxxxyz[i] * pa_y[i];

        tr_y_yyyyyy_xxxxzz[i] = 5.0 * tr_y_yyyy_xxxxzz[i] * fe_0 + ts_yyyyy_xxxxzz[i] * fe_0 + tr_y_yyyyy_xxxxzz[i] * pa_y[i];

        tr_y_yyyyyy_xxxyyy[i] = 5.0 * tr_y_yyyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_yyyyy_xxxyy[i] * fe_0 + ts_yyyyy_xxxyyy[i] * fe_0 + tr_y_yyyyy_xxxyyy[i] * pa_y[i];

        tr_y_yyyyyy_xxxyyz[i] = 5.0 * tr_y_yyyy_xxxyyz[i] * fe_0 + 2.0 * tr_y_yyyyy_xxxyz[i] * fe_0 + ts_yyyyy_xxxyyz[i] * fe_0 + tr_y_yyyyy_xxxyyz[i] * pa_y[i];

        tr_y_yyyyyy_xxxyzz[i] = 5.0 * tr_y_yyyy_xxxyzz[i] * fe_0 + tr_y_yyyyy_xxxzz[i] * fe_0 + ts_yyyyy_xxxyzz[i] * fe_0 + tr_y_yyyyy_xxxyzz[i] * pa_y[i];

        tr_y_yyyyyy_xxxzzz[i] = 5.0 * tr_y_yyyy_xxxzzz[i] * fe_0 + ts_yyyyy_xxxzzz[i] * fe_0 + tr_y_yyyyy_xxxzzz[i] * pa_y[i];

        tr_y_yyyyyy_xxyyyy[i] = 5.0 * tr_y_yyyy_xxyyyy[i] * fe_0 + 4.0 * tr_y_yyyyy_xxyyy[i] * fe_0 + ts_yyyyy_xxyyyy[i] * fe_0 + tr_y_yyyyy_xxyyyy[i] * pa_y[i];

        tr_y_yyyyyy_xxyyyz[i] = 5.0 * tr_y_yyyy_xxyyyz[i] * fe_0 + 3.0 * tr_y_yyyyy_xxyyz[i] * fe_0 + ts_yyyyy_xxyyyz[i] * fe_0 + tr_y_yyyyy_xxyyyz[i] * pa_y[i];

        tr_y_yyyyyy_xxyyzz[i] = 5.0 * tr_y_yyyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyyyy_xxyzz[i] * fe_0 + ts_yyyyy_xxyyzz[i] * fe_0 + tr_y_yyyyy_xxyyzz[i] * pa_y[i];

        tr_y_yyyyyy_xxyzzz[i] = 5.0 * tr_y_yyyy_xxyzzz[i] * fe_0 + tr_y_yyyyy_xxzzz[i] * fe_0 + ts_yyyyy_xxyzzz[i] * fe_0 + tr_y_yyyyy_xxyzzz[i] * pa_y[i];

        tr_y_yyyyyy_xxzzzz[i] = 5.0 * tr_y_yyyy_xxzzzz[i] * fe_0 + ts_yyyyy_xxzzzz[i] * fe_0 + tr_y_yyyyy_xxzzzz[i] * pa_y[i];

        tr_y_yyyyyy_xyyyyy[i] = 5.0 * tr_y_yyyy_xyyyyy[i] * fe_0 + 5.0 * tr_y_yyyyy_xyyyy[i] * fe_0 + ts_yyyyy_xyyyyy[i] * fe_0 + tr_y_yyyyy_xyyyyy[i] * pa_y[i];

        tr_y_yyyyyy_xyyyyz[i] = 5.0 * tr_y_yyyy_xyyyyz[i] * fe_0 + 4.0 * tr_y_yyyyy_xyyyz[i] * fe_0 + ts_yyyyy_xyyyyz[i] * fe_0 + tr_y_yyyyy_xyyyyz[i] * pa_y[i];

        tr_y_yyyyyy_xyyyzz[i] = 5.0 * tr_y_yyyy_xyyyzz[i] * fe_0 + 3.0 * tr_y_yyyyy_xyyzz[i] * fe_0 + ts_yyyyy_xyyyzz[i] * fe_0 + tr_y_yyyyy_xyyyzz[i] * pa_y[i];

        tr_y_yyyyyy_xyyzzz[i] = 5.0 * tr_y_yyyy_xyyzzz[i] * fe_0 + 2.0 * tr_y_yyyyy_xyzzz[i] * fe_0 + ts_yyyyy_xyyzzz[i] * fe_0 + tr_y_yyyyy_xyyzzz[i] * pa_y[i];

        tr_y_yyyyyy_xyzzzz[i] = 5.0 * tr_y_yyyy_xyzzzz[i] * fe_0 + tr_y_yyyyy_xzzzz[i] * fe_0 + ts_yyyyy_xyzzzz[i] * fe_0 + tr_y_yyyyy_xyzzzz[i] * pa_y[i];

        tr_y_yyyyyy_xzzzzz[i] = 5.0 * tr_y_yyyy_xzzzzz[i] * fe_0 + ts_yyyyy_xzzzzz[i] * fe_0 + tr_y_yyyyy_xzzzzz[i] * pa_y[i];

        tr_y_yyyyyy_yyyyyy[i] = 5.0 * tr_y_yyyy_yyyyyy[i] * fe_0 + 6.0 * tr_y_yyyyy_yyyyy[i] * fe_0 + ts_yyyyy_yyyyyy[i] * fe_0 + tr_y_yyyyy_yyyyyy[i] * pa_y[i];

        tr_y_yyyyyy_yyyyyz[i] = 5.0 * tr_y_yyyy_yyyyyz[i] * fe_0 + 5.0 * tr_y_yyyyy_yyyyz[i] * fe_0 + ts_yyyyy_yyyyyz[i] * fe_0 + tr_y_yyyyy_yyyyyz[i] * pa_y[i];

        tr_y_yyyyyy_yyyyzz[i] = 5.0 * tr_y_yyyy_yyyyzz[i] * fe_0 + 4.0 * tr_y_yyyyy_yyyzz[i] * fe_0 + ts_yyyyy_yyyyzz[i] * fe_0 + tr_y_yyyyy_yyyyzz[i] * pa_y[i];

        tr_y_yyyyyy_yyyzzz[i] = 5.0 * tr_y_yyyy_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyyyy_yyzzz[i] * fe_0 + ts_yyyyy_yyyzzz[i] * fe_0 + tr_y_yyyyy_yyyzzz[i] * pa_y[i];

        tr_y_yyyyyy_yyzzzz[i] = 5.0 * tr_y_yyyy_yyzzzz[i] * fe_0 + 2.0 * tr_y_yyyyy_yzzzz[i] * fe_0 + ts_yyyyy_yyzzzz[i] * fe_0 + tr_y_yyyyy_yyzzzz[i] * pa_y[i];

        tr_y_yyyyyy_yzzzzz[i] = 5.0 * tr_y_yyyy_yzzzzz[i] * fe_0 + tr_y_yyyyy_zzzzz[i] * fe_0 + ts_yyyyy_yzzzzz[i] * fe_0 + tr_y_yyyyy_yzzzzz[i] * pa_y[i];

        tr_y_yyyyyy_zzzzzz[i] = 5.0 * tr_y_yyyy_zzzzzz[i] * fe_0 + ts_yyyyy_zzzzzz[i] * fe_0 + tr_y_yyyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 1400-1428 components of targeted buffer : II

    auto tr_y_yyyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 1400);

    auto tr_y_yyyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 1401);

    auto tr_y_yyyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 1402);

    auto tr_y_yyyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 1403);

    auto tr_y_yyyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 1404);

    auto tr_y_yyyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 1405);

    auto tr_y_yyyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 1406);

    auto tr_y_yyyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 1407);

    auto tr_y_yyyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 1408);

    auto tr_y_yyyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 1409);

    auto tr_y_yyyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 1410);

    auto tr_y_yyyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 1411);

    auto tr_y_yyyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 1412);

    auto tr_y_yyyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 1413);

    auto tr_y_yyyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 1414);

    auto tr_y_yyyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 1415);

    auto tr_y_yyyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 1416);

    auto tr_y_yyyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 1417);

    auto tr_y_yyyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 1418);

    auto tr_y_yyyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 1419);

    auto tr_y_yyyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 1420);

    auto tr_y_yyyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 1421);

    auto tr_y_yyyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 1422);

    auto tr_y_yyyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 1423);

    auto tr_y_yyyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 1424);

    auto tr_y_yyyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 1425);

    auto tr_y_yyyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 1426);

    auto tr_y_yyyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 1427);

    #pragma omp simd aligned(pa_z, tr_y_yyyyy_xxxxx, tr_y_yyyyy_xxxxxx, tr_y_yyyyy_xxxxxy, tr_y_yyyyy_xxxxxz, tr_y_yyyyy_xxxxy, tr_y_yyyyy_xxxxyy, tr_y_yyyyy_xxxxyz, tr_y_yyyyy_xxxxz, tr_y_yyyyy_xxxxzz, tr_y_yyyyy_xxxyy, tr_y_yyyyy_xxxyyy, tr_y_yyyyy_xxxyyz, tr_y_yyyyy_xxxyz, tr_y_yyyyy_xxxyzz, tr_y_yyyyy_xxxzz, tr_y_yyyyy_xxxzzz, tr_y_yyyyy_xxyyy, tr_y_yyyyy_xxyyyy, tr_y_yyyyy_xxyyyz, tr_y_yyyyy_xxyyz, tr_y_yyyyy_xxyyzz, tr_y_yyyyy_xxyzz, tr_y_yyyyy_xxyzzz, tr_y_yyyyy_xxzzz, tr_y_yyyyy_xxzzzz, tr_y_yyyyy_xyyyy, tr_y_yyyyy_xyyyyy, tr_y_yyyyy_xyyyyz, tr_y_yyyyy_xyyyz, tr_y_yyyyy_xyyyzz, tr_y_yyyyy_xyyzz, tr_y_yyyyy_xyyzzz, tr_y_yyyyy_xyzzz, tr_y_yyyyy_xyzzzz, tr_y_yyyyy_xzzzz, tr_y_yyyyy_xzzzzz, tr_y_yyyyy_yyyyy, tr_y_yyyyy_yyyyyy, tr_y_yyyyy_yyyyyz, tr_y_yyyyy_yyyyz, tr_y_yyyyy_yyyyzz, tr_y_yyyyy_yyyzz, tr_y_yyyyy_yyyzzz, tr_y_yyyyy_yyzzz, tr_y_yyyyy_yyzzzz, tr_y_yyyyy_yzzzz, tr_y_yyyyy_yzzzzz, tr_y_yyyyy_zzzzz, tr_y_yyyyy_zzzzzz, tr_y_yyyyyz_xxxxxx, tr_y_yyyyyz_xxxxxy, tr_y_yyyyyz_xxxxxz, tr_y_yyyyyz_xxxxyy, tr_y_yyyyyz_xxxxyz, tr_y_yyyyyz_xxxxzz, tr_y_yyyyyz_xxxyyy, tr_y_yyyyyz_xxxyyz, tr_y_yyyyyz_xxxyzz, tr_y_yyyyyz_xxxzzz, tr_y_yyyyyz_xxyyyy, tr_y_yyyyyz_xxyyyz, tr_y_yyyyyz_xxyyzz, tr_y_yyyyyz_xxyzzz, tr_y_yyyyyz_xxzzzz, tr_y_yyyyyz_xyyyyy, tr_y_yyyyyz_xyyyyz, tr_y_yyyyyz_xyyyzz, tr_y_yyyyyz_xyyzzz, tr_y_yyyyyz_xyzzzz, tr_y_yyyyyz_xzzzzz, tr_y_yyyyyz_yyyyyy, tr_y_yyyyyz_yyyyyz, tr_y_yyyyyz_yyyyzz, tr_y_yyyyyz_yyyzzz, tr_y_yyyyyz_yyzzzz, tr_y_yyyyyz_yzzzzz, tr_y_yyyyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyz_xxxxxx[i] = tr_y_yyyyy_xxxxxx[i] * pa_z[i];

        tr_y_yyyyyz_xxxxxy[i] = tr_y_yyyyy_xxxxxy[i] * pa_z[i];

        tr_y_yyyyyz_xxxxxz[i] = tr_y_yyyyy_xxxxx[i] * fe_0 + tr_y_yyyyy_xxxxxz[i] * pa_z[i];

        tr_y_yyyyyz_xxxxyy[i] = tr_y_yyyyy_xxxxyy[i] * pa_z[i];

        tr_y_yyyyyz_xxxxyz[i] = tr_y_yyyyy_xxxxy[i] * fe_0 + tr_y_yyyyy_xxxxyz[i] * pa_z[i];

        tr_y_yyyyyz_xxxxzz[i] = 2.0 * tr_y_yyyyy_xxxxz[i] * fe_0 + tr_y_yyyyy_xxxxzz[i] * pa_z[i];

        tr_y_yyyyyz_xxxyyy[i] = tr_y_yyyyy_xxxyyy[i] * pa_z[i];

        tr_y_yyyyyz_xxxyyz[i] = tr_y_yyyyy_xxxyy[i] * fe_0 + tr_y_yyyyy_xxxyyz[i] * pa_z[i];

        tr_y_yyyyyz_xxxyzz[i] = 2.0 * tr_y_yyyyy_xxxyz[i] * fe_0 + tr_y_yyyyy_xxxyzz[i] * pa_z[i];

        tr_y_yyyyyz_xxxzzz[i] = 3.0 * tr_y_yyyyy_xxxzz[i] * fe_0 + tr_y_yyyyy_xxxzzz[i] * pa_z[i];

        tr_y_yyyyyz_xxyyyy[i] = tr_y_yyyyy_xxyyyy[i] * pa_z[i];

        tr_y_yyyyyz_xxyyyz[i] = tr_y_yyyyy_xxyyy[i] * fe_0 + tr_y_yyyyy_xxyyyz[i] * pa_z[i];

        tr_y_yyyyyz_xxyyzz[i] = 2.0 * tr_y_yyyyy_xxyyz[i] * fe_0 + tr_y_yyyyy_xxyyzz[i] * pa_z[i];

        tr_y_yyyyyz_xxyzzz[i] = 3.0 * tr_y_yyyyy_xxyzz[i] * fe_0 + tr_y_yyyyy_xxyzzz[i] * pa_z[i];

        tr_y_yyyyyz_xxzzzz[i] = 4.0 * tr_y_yyyyy_xxzzz[i] * fe_0 + tr_y_yyyyy_xxzzzz[i] * pa_z[i];

        tr_y_yyyyyz_xyyyyy[i] = tr_y_yyyyy_xyyyyy[i] * pa_z[i];

        tr_y_yyyyyz_xyyyyz[i] = tr_y_yyyyy_xyyyy[i] * fe_0 + tr_y_yyyyy_xyyyyz[i] * pa_z[i];

        tr_y_yyyyyz_xyyyzz[i] = 2.0 * tr_y_yyyyy_xyyyz[i] * fe_0 + tr_y_yyyyy_xyyyzz[i] * pa_z[i];

        tr_y_yyyyyz_xyyzzz[i] = 3.0 * tr_y_yyyyy_xyyzz[i] * fe_0 + tr_y_yyyyy_xyyzzz[i] * pa_z[i];

        tr_y_yyyyyz_xyzzzz[i] = 4.0 * tr_y_yyyyy_xyzzz[i] * fe_0 + tr_y_yyyyy_xyzzzz[i] * pa_z[i];

        tr_y_yyyyyz_xzzzzz[i] = 5.0 * tr_y_yyyyy_xzzzz[i] * fe_0 + tr_y_yyyyy_xzzzzz[i] * pa_z[i];

        tr_y_yyyyyz_yyyyyy[i] = tr_y_yyyyy_yyyyyy[i] * pa_z[i];

        tr_y_yyyyyz_yyyyyz[i] = tr_y_yyyyy_yyyyy[i] * fe_0 + tr_y_yyyyy_yyyyyz[i] * pa_z[i];

        tr_y_yyyyyz_yyyyzz[i] = 2.0 * tr_y_yyyyy_yyyyz[i] * fe_0 + tr_y_yyyyy_yyyyzz[i] * pa_z[i];

        tr_y_yyyyyz_yyyzzz[i] = 3.0 * tr_y_yyyyy_yyyzz[i] * fe_0 + tr_y_yyyyy_yyyzzz[i] * pa_z[i];

        tr_y_yyyyyz_yyzzzz[i] = 4.0 * tr_y_yyyyy_yyzzz[i] * fe_0 + tr_y_yyyyy_yyzzzz[i] * pa_z[i];

        tr_y_yyyyyz_yzzzzz[i] = 5.0 * tr_y_yyyyy_yzzzz[i] * fe_0 + tr_y_yyyyy_yzzzzz[i] * pa_z[i];

        tr_y_yyyyyz_zzzzzz[i] = 6.0 * tr_y_yyyyy_zzzzz[i] * fe_0 + tr_y_yyyyy_zzzzzz[i] * pa_z[i];
    }

    // Set up 1428-1456 components of targeted buffer : II

    auto tr_y_yyyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 1428);

    auto tr_y_yyyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 1429);

    auto tr_y_yyyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 1430);

    auto tr_y_yyyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 1431);

    auto tr_y_yyyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 1432);

    auto tr_y_yyyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 1433);

    auto tr_y_yyyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 1434);

    auto tr_y_yyyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 1435);

    auto tr_y_yyyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 1436);

    auto tr_y_yyyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 1437);

    auto tr_y_yyyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 1438);

    auto tr_y_yyyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 1439);

    auto tr_y_yyyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 1440);

    auto tr_y_yyyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 1441);

    auto tr_y_yyyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 1442);

    auto tr_y_yyyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 1443);

    auto tr_y_yyyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 1444);

    auto tr_y_yyyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 1445);

    auto tr_y_yyyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 1446);

    auto tr_y_yyyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 1447);

    auto tr_y_yyyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 1448);

    auto tr_y_yyyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 1449);

    auto tr_y_yyyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 1450);

    auto tr_y_yyyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 1451);

    auto tr_y_yyyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 1452);

    auto tr_y_yyyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 1453);

    auto tr_y_yyyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 1454);

    auto tr_y_yyyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 1455);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyyy_xxxxxx, tr_y_yyyy_xxxxxy, tr_y_yyyy_xxxxyy, tr_y_yyyy_xxxxyz, tr_y_yyyy_xxxyyy, tr_y_yyyy_xxxyyz, tr_y_yyyy_xxxyzz, tr_y_yyyy_xxyyyy, tr_y_yyyy_xxyyyz, tr_y_yyyy_xxyyzz, tr_y_yyyy_xxyzzz, tr_y_yyyy_xyyyyy, tr_y_yyyy_xyyyyz, tr_y_yyyy_xyyyzz, tr_y_yyyy_xyyzzz, tr_y_yyyy_xyzzzz, tr_y_yyyy_yyyyyy, tr_y_yyyy_yyyyyz, tr_y_yyyy_yyyyzz, tr_y_yyyy_yyyzzz, tr_y_yyyy_yyzzzz, tr_y_yyyy_yzzzzz, tr_y_yyyyz_xxxxxx, tr_y_yyyyz_xxxxxy, tr_y_yyyyz_xxxxy, tr_y_yyyyz_xxxxyy, tr_y_yyyyz_xxxxyz, tr_y_yyyyz_xxxyy, tr_y_yyyyz_xxxyyy, tr_y_yyyyz_xxxyyz, tr_y_yyyyz_xxxyz, tr_y_yyyyz_xxxyzz, tr_y_yyyyz_xxyyy, tr_y_yyyyz_xxyyyy, tr_y_yyyyz_xxyyyz, tr_y_yyyyz_xxyyz, tr_y_yyyyz_xxyyzz, tr_y_yyyyz_xxyzz, tr_y_yyyyz_xxyzzz, tr_y_yyyyz_xyyyy, tr_y_yyyyz_xyyyyy, tr_y_yyyyz_xyyyyz, tr_y_yyyyz_xyyyz, tr_y_yyyyz_xyyyzz, tr_y_yyyyz_xyyzz, tr_y_yyyyz_xyyzzz, tr_y_yyyyz_xyzzz, tr_y_yyyyz_xyzzzz, tr_y_yyyyz_yyyyy, tr_y_yyyyz_yyyyyy, tr_y_yyyyz_yyyyyz, tr_y_yyyyz_yyyyz, tr_y_yyyyz_yyyyzz, tr_y_yyyyz_yyyzz, tr_y_yyyyz_yyyzzz, tr_y_yyyyz_yyzzz, tr_y_yyyyz_yyzzzz, tr_y_yyyyz_yzzzz, tr_y_yyyyz_yzzzzz, tr_y_yyyyzz_xxxxxx, tr_y_yyyyzz_xxxxxy, tr_y_yyyyzz_xxxxxz, tr_y_yyyyzz_xxxxyy, tr_y_yyyyzz_xxxxyz, tr_y_yyyyzz_xxxxzz, tr_y_yyyyzz_xxxyyy, tr_y_yyyyzz_xxxyyz, tr_y_yyyyzz_xxxyzz, tr_y_yyyyzz_xxxzzz, tr_y_yyyyzz_xxyyyy, tr_y_yyyyzz_xxyyyz, tr_y_yyyyzz_xxyyzz, tr_y_yyyyzz_xxyzzz, tr_y_yyyyzz_xxzzzz, tr_y_yyyyzz_xyyyyy, tr_y_yyyyzz_xyyyyz, tr_y_yyyyzz_xyyyzz, tr_y_yyyyzz_xyyzzz, tr_y_yyyyzz_xyzzzz, tr_y_yyyyzz_xzzzzz, tr_y_yyyyzz_yyyyyy, tr_y_yyyyzz_yyyyyz, tr_y_yyyyzz_yyyyzz, tr_y_yyyyzz_yyyzzz, tr_y_yyyyzz_yyzzzz, tr_y_yyyyzz_yzzzzz, tr_y_yyyyzz_zzzzzz, tr_y_yyyzz_xxxxxz, tr_y_yyyzz_xxxxzz, tr_y_yyyzz_xxxzzz, tr_y_yyyzz_xxzzzz, tr_y_yyyzz_xzzzzz, tr_y_yyyzz_zzzzzz, tr_y_yyzz_xxxxxz, tr_y_yyzz_xxxxzz, tr_y_yyzz_xxxzzz, tr_y_yyzz_xxzzzz, tr_y_yyzz_xzzzzz, tr_y_yyzz_zzzzzz, ts_yyyzz_xxxxxz, ts_yyyzz_xxxxzz, ts_yyyzz_xxxzzz, ts_yyyzz_xxzzzz, ts_yyyzz_xzzzzz, ts_yyyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyzz_xxxxxx[i] = tr_y_yyyy_xxxxxx[i] * fe_0 + tr_y_yyyyz_xxxxxx[i] * pa_z[i];

        tr_y_yyyyzz_xxxxxy[i] = tr_y_yyyy_xxxxxy[i] * fe_0 + tr_y_yyyyz_xxxxxy[i] * pa_z[i];

        tr_y_yyyyzz_xxxxxz[i] = 3.0 * tr_y_yyzz_xxxxxz[i] * fe_0 + ts_yyyzz_xxxxxz[i] * fe_0 + tr_y_yyyzz_xxxxxz[i] * pa_y[i];

        tr_y_yyyyzz_xxxxyy[i] = tr_y_yyyy_xxxxyy[i] * fe_0 + tr_y_yyyyz_xxxxyy[i] * pa_z[i];

        tr_y_yyyyzz_xxxxyz[i] = tr_y_yyyy_xxxxyz[i] * fe_0 + tr_y_yyyyz_xxxxy[i] * fe_0 + tr_y_yyyyz_xxxxyz[i] * pa_z[i];

        tr_y_yyyyzz_xxxxzz[i] = 3.0 * tr_y_yyzz_xxxxzz[i] * fe_0 + ts_yyyzz_xxxxzz[i] * fe_0 + tr_y_yyyzz_xxxxzz[i] * pa_y[i];

        tr_y_yyyyzz_xxxyyy[i] = tr_y_yyyy_xxxyyy[i] * fe_0 + tr_y_yyyyz_xxxyyy[i] * pa_z[i];

        tr_y_yyyyzz_xxxyyz[i] = tr_y_yyyy_xxxyyz[i] * fe_0 + tr_y_yyyyz_xxxyy[i] * fe_0 + tr_y_yyyyz_xxxyyz[i] * pa_z[i];

        tr_y_yyyyzz_xxxyzz[i] = tr_y_yyyy_xxxyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_xxxyz[i] * fe_0 + tr_y_yyyyz_xxxyzz[i] * pa_z[i];

        tr_y_yyyyzz_xxxzzz[i] = 3.0 * tr_y_yyzz_xxxzzz[i] * fe_0 + ts_yyyzz_xxxzzz[i] * fe_0 + tr_y_yyyzz_xxxzzz[i] * pa_y[i];

        tr_y_yyyyzz_xxyyyy[i] = tr_y_yyyy_xxyyyy[i] * fe_0 + tr_y_yyyyz_xxyyyy[i] * pa_z[i];

        tr_y_yyyyzz_xxyyyz[i] = tr_y_yyyy_xxyyyz[i] * fe_0 + tr_y_yyyyz_xxyyy[i] * fe_0 + tr_y_yyyyz_xxyyyz[i] * pa_z[i];

        tr_y_yyyyzz_xxyyzz[i] = tr_y_yyyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_xxyyz[i] * fe_0 + tr_y_yyyyz_xxyyzz[i] * pa_z[i];

        tr_y_yyyyzz_xxyzzz[i] = tr_y_yyyy_xxyzzz[i] * fe_0 + 3.0 * tr_y_yyyyz_xxyzz[i] * fe_0 + tr_y_yyyyz_xxyzzz[i] * pa_z[i];

        tr_y_yyyyzz_xxzzzz[i] = 3.0 * tr_y_yyzz_xxzzzz[i] * fe_0 + ts_yyyzz_xxzzzz[i] * fe_0 + tr_y_yyyzz_xxzzzz[i] * pa_y[i];

        tr_y_yyyyzz_xyyyyy[i] = tr_y_yyyy_xyyyyy[i] * fe_0 + tr_y_yyyyz_xyyyyy[i] * pa_z[i];

        tr_y_yyyyzz_xyyyyz[i] = tr_y_yyyy_xyyyyz[i] * fe_0 + tr_y_yyyyz_xyyyy[i] * fe_0 + tr_y_yyyyz_xyyyyz[i] * pa_z[i];

        tr_y_yyyyzz_xyyyzz[i] = tr_y_yyyy_xyyyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_xyyyz[i] * fe_0 + tr_y_yyyyz_xyyyzz[i] * pa_z[i];

        tr_y_yyyyzz_xyyzzz[i] = tr_y_yyyy_xyyzzz[i] * fe_0 + 3.0 * tr_y_yyyyz_xyyzz[i] * fe_0 + tr_y_yyyyz_xyyzzz[i] * pa_z[i];

        tr_y_yyyyzz_xyzzzz[i] = tr_y_yyyy_xyzzzz[i] * fe_0 + 4.0 * tr_y_yyyyz_xyzzz[i] * fe_0 + tr_y_yyyyz_xyzzzz[i] * pa_z[i];

        tr_y_yyyyzz_xzzzzz[i] = 3.0 * tr_y_yyzz_xzzzzz[i] * fe_0 + ts_yyyzz_xzzzzz[i] * fe_0 + tr_y_yyyzz_xzzzzz[i] * pa_y[i];

        tr_y_yyyyzz_yyyyyy[i] = tr_y_yyyy_yyyyyy[i] * fe_0 + tr_y_yyyyz_yyyyyy[i] * pa_z[i];

        tr_y_yyyyzz_yyyyyz[i] = tr_y_yyyy_yyyyyz[i] * fe_0 + tr_y_yyyyz_yyyyy[i] * fe_0 + tr_y_yyyyz_yyyyyz[i] * pa_z[i];

        tr_y_yyyyzz_yyyyzz[i] = tr_y_yyyy_yyyyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_yyyyz[i] * fe_0 + tr_y_yyyyz_yyyyzz[i] * pa_z[i];

        tr_y_yyyyzz_yyyzzz[i] = tr_y_yyyy_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyyyz_yyyzz[i] * fe_0 + tr_y_yyyyz_yyyzzz[i] * pa_z[i];

        tr_y_yyyyzz_yyzzzz[i] = tr_y_yyyy_yyzzzz[i] * fe_0 + 4.0 * tr_y_yyyyz_yyzzz[i] * fe_0 + tr_y_yyyyz_yyzzzz[i] * pa_z[i];

        tr_y_yyyyzz_yzzzzz[i] = tr_y_yyyy_yzzzzz[i] * fe_0 + 5.0 * tr_y_yyyyz_yzzzz[i] * fe_0 + tr_y_yyyyz_yzzzzz[i] * pa_z[i];

        tr_y_yyyyzz_zzzzzz[i] = 3.0 * tr_y_yyzz_zzzzzz[i] * fe_0 + ts_yyyzz_zzzzzz[i] * fe_0 + tr_y_yyyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1456-1484 components of targeted buffer : II

    auto tr_y_yyyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1456);

    auto tr_y_yyyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1457);

    auto tr_y_yyyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1458);

    auto tr_y_yyyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1459);

    auto tr_y_yyyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1460);

    auto tr_y_yyyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1461);

    auto tr_y_yyyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1462);

    auto tr_y_yyyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1463);

    auto tr_y_yyyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1464);

    auto tr_y_yyyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1465);

    auto tr_y_yyyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1466);

    auto tr_y_yyyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1467);

    auto tr_y_yyyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1468);

    auto tr_y_yyyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1469);

    auto tr_y_yyyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1470);

    auto tr_y_yyyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1471);

    auto tr_y_yyyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1472);

    auto tr_y_yyyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1473);

    auto tr_y_yyyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1474);

    auto tr_y_yyyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1475);

    auto tr_y_yyyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1476);

    auto tr_y_yyyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1477);

    auto tr_y_yyyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1478);

    auto tr_y_yyyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1479);

    auto tr_y_yyyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1480);

    auto tr_y_yyyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1481);

    auto tr_y_yyyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1482);

    auto tr_y_yyyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1483);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyyz_xxxxxx, tr_y_yyyz_xxxxxy, tr_y_yyyz_xxxxyy, tr_y_yyyz_xxxxyz, tr_y_yyyz_xxxyyy, tr_y_yyyz_xxxyyz, tr_y_yyyz_xxxyzz, tr_y_yyyz_xxyyyy, tr_y_yyyz_xxyyyz, tr_y_yyyz_xxyyzz, tr_y_yyyz_xxyzzz, tr_y_yyyz_xyyyyy, tr_y_yyyz_xyyyyz, tr_y_yyyz_xyyyzz, tr_y_yyyz_xyyzzz, tr_y_yyyz_xyzzzz, tr_y_yyyz_yyyyyy, tr_y_yyyz_yyyyyz, tr_y_yyyz_yyyyzz, tr_y_yyyz_yyyzzz, tr_y_yyyz_yyzzzz, tr_y_yyyz_yzzzzz, tr_y_yyyzz_xxxxxx, tr_y_yyyzz_xxxxxy, tr_y_yyyzz_xxxxy, tr_y_yyyzz_xxxxyy, tr_y_yyyzz_xxxxyz, tr_y_yyyzz_xxxyy, tr_y_yyyzz_xxxyyy, tr_y_yyyzz_xxxyyz, tr_y_yyyzz_xxxyz, tr_y_yyyzz_xxxyzz, tr_y_yyyzz_xxyyy, tr_y_yyyzz_xxyyyy, tr_y_yyyzz_xxyyyz, tr_y_yyyzz_xxyyz, tr_y_yyyzz_xxyyzz, tr_y_yyyzz_xxyzz, tr_y_yyyzz_xxyzzz, tr_y_yyyzz_xyyyy, tr_y_yyyzz_xyyyyy, tr_y_yyyzz_xyyyyz, tr_y_yyyzz_xyyyz, tr_y_yyyzz_xyyyzz, tr_y_yyyzz_xyyzz, tr_y_yyyzz_xyyzzz, tr_y_yyyzz_xyzzz, tr_y_yyyzz_xyzzzz, tr_y_yyyzz_yyyyy, tr_y_yyyzz_yyyyyy, tr_y_yyyzz_yyyyyz, tr_y_yyyzz_yyyyz, tr_y_yyyzz_yyyyzz, tr_y_yyyzz_yyyzz, tr_y_yyyzz_yyyzzz, tr_y_yyyzz_yyzzz, tr_y_yyyzz_yyzzzz, tr_y_yyyzz_yzzzz, tr_y_yyyzz_yzzzzz, tr_y_yyyzzz_xxxxxx, tr_y_yyyzzz_xxxxxy, tr_y_yyyzzz_xxxxxz, tr_y_yyyzzz_xxxxyy, tr_y_yyyzzz_xxxxyz, tr_y_yyyzzz_xxxxzz, tr_y_yyyzzz_xxxyyy, tr_y_yyyzzz_xxxyyz, tr_y_yyyzzz_xxxyzz, tr_y_yyyzzz_xxxzzz, tr_y_yyyzzz_xxyyyy, tr_y_yyyzzz_xxyyyz, tr_y_yyyzzz_xxyyzz, tr_y_yyyzzz_xxyzzz, tr_y_yyyzzz_xxzzzz, tr_y_yyyzzz_xyyyyy, tr_y_yyyzzz_xyyyyz, tr_y_yyyzzz_xyyyzz, tr_y_yyyzzz_xyyzzz, tr_y_yyyzzz_xyzzzz, tr_y_yyyzzz_xzzzzz, tr_y_yyyzzz_yyyyyy, tr_y_yyyzzz_yyyyyz, tr_y_yyyzzz_yyyyzz, tr_y_yyyzzz_yyyzzz, tr_y_yyyzzz_yyzzzz, tr_y_yyyzzz_yzzzzz, tr_y_yyyzzz_zzzzzz, tr_y_yyzzz_xxxxxz, tr_y_yyzzz_xxxxzz, tr_y_yyzzz_xxxzzz, tr_y_yyzzz_xxzzzz, tr_y_yyzzz_xzzzzz, tr_y_yyzzz_zzzzzz, tr_y_yzzz_xxxxxz, tr_y_yzzz_xxxxzz, tr_y_yzzz_xxxzzz, tr_y_yzzz_xxzzzz, tr_y_yzzz_xzzzzz, tr_y_yzzz_zzzzzz, ts_yyzzz_xxxxxz, ts_yyzzz_xxxxzz, ts_yyzzz_xxxzzz, ts_yyzzz_xxzzzz, ts_yyzzz_xzzzzz, ts_yyzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzzz_xxxxxx[i] = 2.0 * tr_y_yyyz_xxxxxx[i] * fe_0 + tr_y_yyyzz_xxxxxx[i] * pa_z[i];

        tr_y_yyyzzz_xxxxxy[i] = 2.0 * tr_y_yyyz_xxxxxy[i] * fe_0 + tr_y_yyyzz_xxxxxy[i] * pa_z[i];

        tr_y_yyyzzz_xxxxxz[i] = 2.0 * tr_y_yzzz_xxxxxz[i] * fe_0 + ts_yyzzz_xxxxxz[i] * fe_0 + tr_y_yyzzz_xxxxxz[i] * pa_y[i];

        tr_y_yyyzzz_xxxxyy[i] = 2.0 * tr_y_yyyz_xxxxyy[i] * fe_0 + tr_y_yyyzz_xxxxyy[i] * pa_z[i];

        tr_y_yyyzzz_xxxxyz[i] = 2.0 * tr_y_yyyz_xxxxyz[i] * fe_0 + tr_y_yyyzz_xxxxy[i] * fe_0 + tr_y_yyyzz_xxxxyz[i] * pa_z[i];

        tr_y_yyyzzz_xxxxzz[i] = 2.0 * tr_y_yzzz_xxxxzz[i] * fe_0 + ts_yyzzz_xxxxzz[i] * fe_0 + tr_y_yyzzz_xxxxzz[i] * pa_y[i];

        tr_y_yyyzzz_xxxyyy[i] = 2.0 * tr_y_yyyz_xxxyyy[i] * fe_0 + tr_y_yyyzz_xxxyyy[i] * pa_z[i];

        tr_y_yyyzzz_xxxyyz[i] = 2.0 * tr_y_yyyz_xxxyyz[i] * fe_0 + tr_y_yyyzz_xxxyy[i] * fe_0 + tr_y_yyyzz_xxxyyz[i] * pa_z[i];

        tr_y_yyyzzz_xxxyzz[i] = 2.0 * tr_y_yyyz_xxxyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_xxxyz[i] * fe_0 + tr_y_yyyzz_xxxyzz[i] * pa_z[i];

        tr_y_yyyzzz_xxxzzz[i] = 2.0 * tr_y_yzzz_xxxzzz[i] * fe_0 + ts_yyzzz_xxxzzz[i] * fe_0 + tr_y_yyzzz_xxxzzz[i] * pa_y[i];

        tr_y_yyyzzz_xxyyyy[i] = 2.0 * tr_y_yyyz_xxyyyy[i] * fe_0 + tr_y_yyyzz_xxyyyy[i] * pa_z[i];

        tr_y_yyyzzz_xxyyyz[i] = 2.0 * tr_y_yyyz_xxyyyz[i] * fe_0 + tr_y_yyyzz_xxyyy[i] * fe_0 + tr_y_yyyzz_xxyyyz[i] * pa_z[i];

        tr_y_yyyzzz_xxyyzz[i] = 2.0 * tr_y_yyyz_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_xxyyz[i] * fe_0 + tr_y_yyyzz_xxyyzz[i] * pa_z[i];

        tr_y_yyyzzz_xxyzzz[i] = 2.0 * tr_y_yyyz_xxyzzz[i] * fe_0 + 3.0 * tr_y_yyyzz_xxyzz[i] * fe_0 + tr_y_yyyzz_xxyzzz[i] * pa_z[i];

        tr_y_yyyzzz_xxzzzz[i] = 2.0 * tr_y_yzzz_xxzzzz[i] * fe_0 + ts_yyzzz_xxzzzz[i] * fe_0 + tr_y_yyzzz_xxzzzz[i] * pa_y[i];

        tr_y_yyyzzz_xyyyyy[i] = 2.0 * tr_y_yyyz_xyyyyy[i] * fe_0 + tr_y_yyyzz_xyyyyy[i] * pa_z[i];

        tr_y_yyyzzz_xyyyyz[i] = 2.0 * tr_y_yyyz_xyyyyz[i] * fe_0 + tr_y_yyyzz_xyyyy[i] * fe_0 + tr_y_yyyzz_xyyyyz[i] * pa_z[i];

        tr_y_yyyzzz_xyyyzz[i] = 2.0 * tr_y_yyyz_xyyyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_xyyyz[i] * fe_0 + tr_y_yyyzz_xyyyzz[i] * pa_z[i];

        tr_y_yyyzzz_xyyzzz[i] = 2.0 * tr_y_yyyz_xyyzzz[i] * fe_0 + 3.0 * tr_y_yyyzz_xyyzz[i] * fe_0 + tr_y_yyyzz_xyyzzz[i] * pa_z[i];

        tr_y_yyyzzz_xyzzzz[i] = 2.0 * tr_y_yyyz_xyzzzz[i] * fe_0 + 4.0 * tr_y_yyyzz_xyzzz[i] * fe_0 + tr_y_yyyzz_xyzzzz[i] * pa_z[i];

        tr_y_yyyzzz_xzzzzz[i] = 2.0 * tr_y_yzzz_xzzzzz[i] * fe_0 + ts_yyzzz_xzzzzz[i] * fe_0 + tr_y_yyzzz_xzzzzz[i] * pa_y[i];

        tr_y_yyyzzz_yyyyyy[i] = 2.0 * tr_y_yyyz_yyyyyy[i] * fe_0 + tr_y_yyyzz_yyyyyy[i] * pa_z[i];

        tr_y_yyyzzz_yyyyyz[i] = 2.0 * tr_y_yyyz_yyyyyz[i] * fe_0 + tr_y_yyyzz_yyyyy[i] * fe_0 + tr_y_yyyzz_yyyyyz[i] * pa_z[i];

        tr_y_yyyzzz_yyyyzz[i] = 2.0 * tr_y_yyyz_yyyyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_yyyyz[i] * fe_0 + tr_y_yyyzz_yyyyzz[i] * pa_z[i];

        tr_y_yyyzzz_yyyzzz[i] = 2.0 * tr_y_yyyz_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyyzz_yyyzz[i] * fe_0 + tr_y_yyyzz_yyyzzz[i] * pa_z[i];

        tr_y_yyyzzz_yyzzzz[i] = 2.0 * tr_y_yyyz_yyzzzz[i] * fe_0 + 4.0 * tr_y_yyyzz_yyzzz[i] * fe_0 + tr_y_yyyzz_yyzzzz[i] * pa_z[i];

        tr_y_yyyzzz_yzzzzz[i] = 2.0 * tr_y_yyyz_yzzzzz[i] * fe_0 + 5.0 * tr_y_yyyzz_yzzzz[i] * fe_0 + tr_y_yyyzz_yzzzzz[i] * pa_z[i];

        tr_y_yyyzzz_zzzzzz[i] = 2.0 * tr_y_yzzz_zzzzzz[i] * fe_0 + ts_yyzzz_zzzzzz[i] * fe_0 + tr_y_yyzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1484-1512 components of targeted buffer : II

    auto tr_y_yyzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1484);

    auto tr_y_yyzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1485);

    auto tr_y_yyzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1486);

    auto tr_y_yyzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1487);

    auto tr_y_yyzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1488);

    auto tr_y_yyzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1489);

    auto tr_y_yyzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1490);

    auto tr_y_yyzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1491);

    auto tr_y_yyzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1492);

    auto tr_y_yyzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1493);

    auto tr_y_yyzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1494);

    auto tr_y_yyzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1495);

    auto tr_y_yyzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1496);

    auto tr_y_yyzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1497);

    auto tr_y_yyzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1498);

    auto tr_y_yyzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1499);

    auto tr_y_yyzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1500);

    auto tr_y_yyzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1501);

    auto tr_y_yyzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1502);

    auto tr_y_yyzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1503);

    auto tr_y_yyzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1504);

    auto tr_y_yyzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1505);

    auto tr_y_yyzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1506);

    auto tr_y_yyzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1507);

    auto tr_y_yyzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1508);

    auto tr_y_yyzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1509);

    auto tr_y_yyzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1510);

    auto tr_y_yyzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1511);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyzz_xxxxxx, tr_y_yyzz_xxxxxy, tr_y_yyzz_xxxxyy, tr_y_yyzz_xxxxyz, tr_y_yyzz_xxxyyy, tr_y_yyzz_xxxyyz, tr_y_yyzz_xxxyzz, tr_y_yyzz_xxyyyy, tr_y_yyzz_xxyyyz, tr_y_yyzz_xxyyzz, tr_y_yyzz_xxyzzz, tr_y_yyzz_xyyyyy, tr_y_yyzz_xyyyyz, tr_y_yyzz_xyyyzz, tr_y_yyzz_xyyzzz, tr_y_yyzz_xyzzzz, tr_y_yyzz_yyyyyy, tr_y_yyzz_yyyyyz, tr_y_yyzz_yyyyzz, tr_y_yyzz_yyyzzz, tr_y_yyzz_yyzzzz, tr_y_yyzz_yzzzzz, tr_y_yyzzz_xxxxxx, tr_y_yyzzz_xxxxxy, tr_y_yyzzz_xxxxy, tr_y_yyzzz_xxxxyy, tr_y_yyzzz_xxxxyz, tr_y_yyzzz_xxxyy, tr_y_yyzzz_xxxyyy, tr_y_yyzzz_xxxyyz, tr_y_yyzzz_xxxyz, tr_y_yyzzz_xxxyzz, tr_y_yyzzz_xxyyy, tr_y_yyzzz_xxyyyy, tr_y_yyzzz_xxyyyz, tr_y_yyzzz_xxyyz, tr_y_yyzzz_xxyyzz, tr_y_yyzzz_xxyzz, tr_y_yyzzz_xxyzzz, tr_y_yyzzz_xyyyy, tr_y_yyzzz_xyyyyy, tr_y_yyzzz_xyyyyz, tr_y_yyzzz_xyyyz, tr_y_yyzzz_xyyyzz, tr_y_yyzzz_xyyzz, tr_y_yyzzz_xyyzzz, tr_y_yyzzz_xyzzz, tr_y_yyzzz_xyzzzz, tr_y_yyzzz_yyyyy, tr_y_yyzzz_yyyyyy, tr_y_yyzzz_yyyyyz, tr_y_yyzzz_yyyyz, tr_y_yyzzz_yyyyzz, tr_y_yyzzz_yyyzz, tr_y_yyzzz_yyyzzz, tr_y_yyzzz_yyzzz, tr_y_yyzzz_yyzzzz, tr_y_yyzzz_yzzzz, tr_y_yyzzz_yzzzzz, tr_y_yyzzzz_xxxxxx, tr_y_yyzzzz_xxxxxy, tr_y_yyzzzz_xxxxxz, tr_y_yyzzzz_xxxxyy, tr_y_yyzzzz_xxxxyz, tr_y_yyzzzz_xxxxzz, tr_y_yyzzzz_xxxyyy, tr_y_yyzzzz_xxxyyz, tr_y_yyzzzz_xxxyzz, tr_y_yyzzzz_xxxzzz, tr_y_yyzzzz_xxyyyy, tr_y_yyzzzz_xxyyyz, tr_y_yyzzzz_xxyyzz, tr_y_yyzzzz_xxyzzz, tr_y_yyzzzz_xxzzzz, tr_y_yyzzzz_xyyyyy, tr_y_yyzzzz_xyyyyz, tr_y_yyzzzz_xyyyzz, tr_y_yyzzzz_xyyzzz, tr_y_yyzzzz_xyzzzz, tr_y_yyzzzz_xzzzzz, tr_y_yyzzzz_yyyyyy, tr_y_yyzzzz_yyyyyz, tr_y_yyzzzz_yyyyzz, tr_y_yyzzzz_yyyzzz, tr_y_yyzzzz_yyzzzz, tr_y_yyzzzz_yzzzzz, tr_y_yyzzzz_zzzzzz, tr_y_yzzzz_xxxxxz, tr_y_yzzzz_xxxxzz, tr_y_yzzzz_xxxzzz, tr_y_yzzzz_xxzzzz, tr_y_yzzzz_xzzzzz, tr_y_yzzzz_zzzzzz, tr_y_zzzz_xxxxxz, tr_y_zzzz_xxxxzz, tr_y_zzzz_xxxzzz, tr_y_zzzz_xxzzzz, tr_y_zzzz_xzzzzz, tr_y_zzzz_zzzzzz, ts_yzzzz_xxxxxz, ts_yzzzz_xxxxzz, ts_yzzzz_xxxzzz, ts_yzzzz_xxzzzz, ts_yzzzz_xzzzzz, ts_yzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzzz_xxxxxx[i] = 3.0 * tr_y_yyzz_xxxxxx[i] * fe_0 + tr_y_yyzzz_xxxxxx[i] * pa_z[i];

        tr_y_yyzzzz_xxxxxy[i] = 3.0 * tr_y_yyzz_xxxxxy[i] * fe_0 + tr_y_yyzzz_xxxxxy[i] * pa_z[i];

        tr_y_yyzzzz_xxxxxz[i] = tr_y_zzzz_xxxxxz[i] * fe_0 + ts_yzzzz_xxxxxz[i] * fe_0 + tr_y_yzzzz_xxxxxz[i] * pa_y[i];

        tr_y_yyzzzz_xxxxyy[i] = 3.0 * tr_y_yyzz_xxxxyy[i] * fe_0 + tr_y_yyzzz_xxxxyy[i] * pa_z[i];

        tr_y_yyzzzz_xxxxyz[i] = 3.0 * tr_y_yyzz_xxxxyz[i] * fe_0 + tr_y_yyzzz_xxxxy[i] * fe_0 + tr_y_yyzzz_xxxxyz[i] * pa_z[i];

        tr_y_yyzzzz_xxxxzz[i] = tr_y_zzzz_xxxxzz[i] * fe_0 + ts_yzzzz_xxxxzz[i] * fe_0 + tr_y_yzzzz_xxxxzz[i] * pa_y[i];

        tr_y_yyzzzz_xxxyyy[i] = 3.0 * tr_y_yyzz_xxxyyy[i] * fe_0 + tr_y_yyzzz_xxxyyy[i] * pa_z[i];

        tr_y_yyzzzz_xxxyyz[i] = 3.0 * tr_y_yyzz_xxxyyz[i] * fe_0 + tr_y_yyzzz_xxxyy[i] * fe_0 + tr_y_yyzzz_xxxyyz[i] * pa_z[i];

        tr_y_yyzzzz_xxxyzz[i] = 3.0 * tr_y_yyzz_xxxyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_xxxyz[i] * fe_0 + tr_y_yyzzz_xxxyzz[i] * pa_z[i];

        tr_y_yyzzzz_xxxzzz[i] = tr_y_zzzz_xxxzzz[i] * fe_0 + ts_yzzzz_xxxzzz[i] * fe_0 + tr_y_yzzzz_xxxzzz[i] * pa_y[i];

        tr_y_yyzzzz_xxyyyy[i] = 3.0 * tr_y_yyzz_xxyyyy[i] * fe_0 + tr_y_yyzzz_xxyyyy[i] * pa_z[i];

        tr_y_yyzzzz_xxyyyz[i] = 3.0 * tr_y_yyzz_xxyyyz[i] * fe_0 + tr_y_yyzzz_xxyyy[i] * fe_0 + tr_y_yyzzz_xxyyyz[i] * pa_z[i];

        tr_y_yyzzzz_xxyyzz[i] = 3.0 * tr_y_yyzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_xxyyz[i] * fe_0 + tr_y_yyzzz_xxyyzz[i] * pa_z[i];

        tr_y_yyzzzz_xxyzzz[i] = 3.0 * tr_y_yyzz_xxyzzz[i] * fe_0 + 3.0 * tr_y_yyzzz_xxyzz[i] * fe_0 + tr_y_yyzzz_xxyzzz[i] * pa_z[i];

        tr_y_yyzzzz_xxzzzz[i] = tr_y_zzzz_xxzzzz[i] * fe_0 + ts_yzzzz_xxzzzz[i] * fe_0 + tr_y_yzzzz_xxzzzz[i] * pa_y[i];

        tr_y_yyzzzz_xyyyyy[i] = 3.0 * tr_y_yyzz_xyyyyy[i] * fe_0 + tr_y_yyzzz_xyyyyy[i] * pa_z[i];

        tr_y_yyzzzz_xyyyyz[i] = 3.0 * tr_y_yyzz_xyyyyz[i] * fe_0 + tr_y_yyzzz_xyyyy[i] * fe_0 + tr_y_yyzzz_xyyyyz[i] * pa_z[i];

        tr_y_yyzzzz_xyyyzz[i] = 3.0 * tr_y_yyzz_xyyyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_xyyyz[i] * fe_0 + tr_y_yyzzz_xyyyzz[i] * pa_z[i];

        tr_y_yyzzzz_xyyzzz[i] = 3.0 * tr_y_yyzz_xyyzzz[i] * fe_0 + 3.0 * tr_y_yyzzz_xyyzz[i] * fe_0 + tr_y_yyzzz_xyyzzz[i] * pa_z[i];

        tr_y_yyzzzz_xyzzzz[i] = 3.0 * tr_y_yyzz_xyzzzz[i] * fe_0 + 4.0 * tr_y_yyzzz_xyzzz[i] * fe_0 + tr_y_yyzzz_xyzzzz[i] * pa_z[i];

        tr_y_yyzzzz_xzzzzz[i] = tr_y_zzzz_xzzzzz[i] * fe_0 + ts_yzzzz_xzzzzz[i] * fe_0 + tr_y_yzzzz_xzzzzz[i] * pa_y[i];

        tr_y_yyzzzz_yyyyyy[i] = 3.0 * tr_y_yyzz_yyyyyy[i] * fe_0 + tr_y_yyzzz_yyyyyy[i] * pa_z[i];

        tr_y_yyzzzz_yyyyyz[i] = 3.0 * tr_y_yyzz_yyyyyz[i] * fe_0 + tr_y_yyzzz_yyyyy[i] * fe_0 + tr_y_yyzzz_yyyyyz[i] * pa_z[i];

        tr_y_yyzzzz_yyyyzz[i] = 3.0 * tr_y_yyzz_yyyyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_yyyyz[i] * fe_0 + tr_y_yyzzz_yyyyzz[i] * pa_z[i];

        tr_y_yyzzzz_yyyzzz[i] = 3.0 * tr_y_yyzz_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyzzz_yyyzz[i] * fe_0 + tr_y_yyzzz_yyyzzz[i] * pa_z[i];

        tr_y_yyzzzz_yyzzzz[i] = 3.0 * tr_y_yyzz_yyzzzz[i] * fe_0 + 4.0 * tr_y_yyzzz_yyzzz[i] * fe_0 + tr_y_yyzzz_yyzzzz[i] * pa_z[i];

        tr_y_yyzzzz_yzzzzz[i] = 3.0 * tr_y_yyzz_yzzzzz[i] * fe_0 + 5.0 * tr_y_yyzzz_yzzzz[i] * fe_0 + tr_y_yyzzz_yzzzzz[i] * pa_z[i];

        tr_y_yyzzzz_zzzzzz[i] = tr_y_zzzz_zzzzzz[i] * fe_0 + ts_yzzzz_zzzzzz[i] * fe_0 + tr_y_yzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1512-1540 components of targeted buffer : II

    auto tr_y_yzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1512);

    auto tr_y_yzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1513);

    auto tr_y_yzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1514);

    auto tr_y_yzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1515);

    auto tr_y_yzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1516);

    auto tr_y_yzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1517);

    auto tr_y_yzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1518);

    auto tr_y_yzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1519);

    auto tr_y_yzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1520);

    auto tr_y_yzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1521);

    auto tr_y_yzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1522);

    auto tr_y_yzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1523);

    auto tr_y_yzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1524);

    auto tr_y_yzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1525);

    auto tr_y_yzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1526);

    auto tr_y_yzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1527);

    auto tr_y_yzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1528);

    auto tr_y_yzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1529);

    auto tr_y_yzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1530);

    auto tr_y_yzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1531);

    auto tr_y_yzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1532);

    auto tr_y_yzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1533);

    auto tr_y_yzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1534);

    auto tr_y_yzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1535);

    auto tr_y_yzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1536);

    auto tr_y_yzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1537);

    auto tr_y_yzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1538);

    auto tr_y_yzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1539);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yzzz_xxxxxy, tr_y_yzzz_xxxxyy, tr_y_yzzz_xxxyyy, tr_y_yzzz_xxyyyy, tr_y_yzzz_xyyyyy, tr_y_yzzz_yyyyyy, tr_y_yzzzz_xxxxxy, tr_y_yzzzz_xxxxyy, tr_y_yzzzz_xxxyyy, tr_y_yzzzz_xxyyyy, tr_y_yzzzz_xyyyyy, tr_y_yzzzz_yyyyyy, tr_y_yzzzzz_xxxxxx, tr_y_yzzzzz_xxxxxy, tr_y_yzzzzz_xxxxxz, tr_y_yzzzzz_xxxxyy, tr_y_yzzzzz_xxxxyz, tr_y_yzzzzz_xxxxzz, tr_y_yzzzzz_xxxyyy, tr_y_yzzzzz_xxxyyz, tr_y_yzzzzz_xxxyzz, tr_y_yzzzzz_xxxzzz, tr_y_yzzzzz_xxyyyy, tr_y_yzzzzz_xxyyyz, tr_y_yzzzzz_xxyyzz, tr_y_yzzzzz_xxyzzz, tr_y_yzzzzz_xxzzzz, tr_y_yzzzzz_xyyyyy, tr_y_yzzzzz_xyyyyz, tr_y_yzzzzz_xyyyzz, tr_y_yzzzzz_xyyzzz, tr_y_yzzzzz_xyzzzz, tr_y_yzzzzz_xzzzzz, tr_y_yzzzzz_yyyyyy, tr_y_yzzzzz_yyyyyz, tr_y_yzzzzz_yyyyzz, tr_y_yzzzzz_yyyzzz, tr_y_yzzzzz_yyzzzz, tr_y_yzzzzz_yzzzzz, tr_y_yzzzzz_zzzzzz, tr_y_zzzzz_xxxxxx, tr_y_zzzzz_xxxxxz, tr_y_zzzzz_xxxxyz, tr_y_zzzzz_xxxxz, tr_y_zzzzz_xxxxzz, tr_y_zzzzz_xxxyyz, tr_y_zzzzz_xxxyz, tr_y_zzzzz_xxxyzz, tr_y_zzzzz_xxxzz, tr_y_zzzzz_xxxzzz, tr_y_zzzzz_xxyyyz, tr_y_zzzzz_xxyyz, tr_y_zzzzz_xxyyzz, tr_y_zzzzz_xxyzz, tr_y_zzzzz_xxyzzz, tr_y_zzzzz_xxzzz, tr_y_zzzzz_xxzzzz, tr_y_zzzzz_xyyyyz, tr_y_zzzzz_xyyyz, tr_y_zzzzz_xyyyzz, tr_y_zzzzz_xyyzz, tr_y_zzzzz_xyyzzz, tr_y_zzzzz_xyzzz, tr_y_zzzzz_xyzzzz, tr_y_zzzzz_xzzzz, tr_y_zzzzz_xzzzzz, tr_y_zzzzz_yyyyyz, tr_y_zzzzz_yyyyz, tr_y_zzzzz_yyyyzz, tr_y_zzzzz_yyyzz, tr_y_zzzzz_yyyzzz, tr_y_zzzzz_yyzzz, tr_y_zzzzz_yyzzzz, tr_y_zzzzz_yzzzz, tr_y_zzzzz_yzzzzz, tr_y_zzzzz_zzzzz, tr_y_zzzzz_zzzzzz, ts_zzzzz_xxxxxx, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzzz_xxxxxx[i] = ts_zzzzz_xxxxxx[i] * fe_0 + tr_y_zzzzz_xxxxxx[i] * pa_y[i];

        tr_y_yzzzzz_xxxxxy[i] = 4.0 * tr_y_yzzz_xxxxxy[i] * fe_0 + tr_y_yzzzz_xxxxxy[i] * pa_z[i];

        tr_y_yzzzzz_xxxxxz[i] = ts_zzzzz_xxxxxz[i] * fe_0 + tr_y_zzzzz_xxxxxz[i] * pa_y[i];

        tr_y_yzzzzz_xxxxyy[i] = 4.0 * tr_y_yzzz_xxxxyy[i] * fe_0 + tr_y_yzzzz_xxxxyy[i] * pa_z[i];

        tr_y_yzzzzz_xxxxyz[i] = tr_y_zzzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxxyz[i] * fe_0 + tr_y_zzzzz_xxxxyz[i] * pa_y[i];

        tr_y_yzzzzz_xxxxzz[i] = ts_zzzzz_xxxxzz[i] * fe_0 + tr_y_zzzzz_xxxxzz[i] * pa_y[i];

        tr_y_yzzzzz_xxxyyy[i] = 4.0 * tr_y_yzzz_xxxyyy[i] * fe_0 + tr_y_yzzzz_xxxyyy[i] * pa_z[i];

        tr_y_yzzzzz_xxxyyz[i] = 2.0 * tr_y_zzzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxyyz[i] * fe_0 + tr_y_zzzzz_xxxyyz[i] * pa_y[i];

        tr_y_yzzzzz_xxxyzz[i] = tr_y_zzzzz_xxxzz[i] * fe_0 + ts_zzzzz_xxxyzz[i] * fe_0 + tr_y_zzzzz_xxxyzz[i] * pa_y[i];

        tr_y_yzzzzz_xxxzzz[i] = ts_zzzzz_xxxzzz[i] * fe_0 + tr_y_zzzzz_xxxzzz[i] * pa_y[i];

        tr_y_yzzzzz_xxyyyy[i] = 4.0 * tr_y_yzzz_xxyyyy[i] * fe_0 + tr_y_yzzzz_xxyyyy[i] * pa_z[i];

        tr_y_yzzzzz_xxyyyz[i] = 3.0 * tr_y_zzzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxyyyz[i] * fe_0 + tr_y_zzzzz_xxyyyz[i] * pa_y[i];

        tr_y_yzzzzz_xxyyzz[i] = 2.0 * tr_y_zzzzz_xxyzz[i] * fe_0 + ts_zzzzz_xxyyzz[i] * fe_0 + tr_y_zzzzz_xxyyzz[i] * pa_y[i];

        tr_y_yzzzzz_xxyzzz[i] = tr_y_zzzzz_xxzzz[i] * fe_0 + ts_zzzzz_xxyzzz[i] * fe_0 + tr_y_zzzzz_xxyzzz[i] * pa_y[i];

        tr_y_yzzzzz_xxzzzz[i] = ts_zzzzz_xxzzzz[i] * fe_0 + tr_y_zzzzz_xxzzzz[i] * pa_y[i];

        tr_y_yzzzzz_xyyyyy[i] = 4.0 * tr_y_yzzz_xyyyyy[i] * fe_0 + tr_y_yzzzz_xyyyyy[i] * pa_z[i];

        tr_y_yzzzzz_xyyyyz[i] = 4.0 * tr_y_zzzzz_xyyyz[i] * fe_0 + ts_zzzzz_xyyyyz[i] * fe_0 + tr_y_zzzzz_xyyyyz[i] * pa_y[i];

        tr_y_yzzzzz_xyyyzz[i] = 3.0 * tr_y_zzzzz_xyyzz[i] * fe_0 + ts_zzzzz_xyyyzz[i] * fe_0 + tr_y_zzzzz_xyyyzz[i] * pa_y[i];

        tr_y_yzzzzz_xyyzzz[i] = 2.0 * tr_y_zzzzz_xyzzz[i] * fe_0 + ts_zzzzz_xyyzzz[i] * fe_0 + tr_y_zzzzz_xyyzzz[i] * pa_y[i];

        tr_y_yzzzzz_xyzzzz[i] = tr_y_zzzzz_xzzzz[i] * fe_0 + ts_zzzzz_xyzzzz[i] * fe_0 + tr_y_zzzzz_xyzzzz[i] * pa_y[i];

        tr_y_yzzzzz_xzzzzz[i] = ts_zzzzz_xzzzzz[i] * fe_0 + tr_y_zzzzz_xzzzzz[i] * pa_y[i];

        tr_y_yzzzzz_yyyyyy[i] = 4.0 * tr_y_yzzz_yyyyyy[i] * fe_0 + tr_y_yzzzz_yyyyyy[i] * pa_z[i];

        tr_y_yzzzzz_yyyyyz[i] = 5.0 * tr_y_zzzzz_yyyyz[i] * fe_0 + ts_zzzzz_yyyyyz[i] * fe_0 + tr_y_zzzzz_yyyyyz[i] * pa_y[i];

        tr_y_yzzzzz_yyyyzz[i] = 4.0 * tr_y_zzzzz_yyyzz[i] * fe_0 + ts_zzzzz_yyyyzz[i] * fe_0 + tr_y_zzzzz_yyyyzz[i] * pa_y[i];

        tr_y_yzzzzz_yyyzzz[i] = 3.0 * tr_y_zzzzz_yyzzz[i] * fe_0 + ts_zzzzz_yyyzzz[i] * fe_0 + tr_y_zzzzz_yyyzzz[i] * pa_y[i];

        tr_y_yzzzzz_yyzzzz[i] = 2.0 * tr_y_zzzzz_yzzzz[i] * fe_0 + ts_zzzzz_yyzzzz[i] * fe_0 + tr_y_zzzzz_yyzzzz[i] * pa_y[i];

        tr_y_yzzzzz_yzzzzz[i] = tr_y_zzzzz_zzzzz[i] * fe_0 + ts_zzzzz_yzzzzz[i] * fe_0 + tr_y_zzzzz_yzzzzz[i] * pa_y[i];

        tr_y_yzzzzz_zzzzzz[i] = ts_zzzzz_zzzzzz[i] * fe_0 + tr_y_zzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1540-1568 components of targeted buffer : II

    auto tr_y_zzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1540);

    auto tr_y_zzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1541);

    auto tr_y_zzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1542);

    auto tr_y_zzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1543);

    auto tr_y_zzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1544);

    auto tr_y_zzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1545);

    auto tr_y_zzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1546);

    auto tr_y_zzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1547);

    auto tr_y_zzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1548);

    auto tr_y_zzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1549);

    auto tr_y_zzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1550);

    auto tr_y_zzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1551);

    auto tr_y_zzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1552);

    auto tr_y_zzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1553);

    auto tr_y_zzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1554);

    auto tr_y_zzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1555);

    auto tr_y_zzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1556);

    auto tr_y_zzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1557);

    auto tr_y_zzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1558);

    auto tr_y_zzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1559);

    auto tr_y_zzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1560);

    auto tr_y_zzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1561);

    auto tr_y_zzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1562);

    auto tr_y_zzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1563);

    auto tr_y_zzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1564);

    auto tr_y_zzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1565);

    auto tr_y_zzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1566);

    auto tr_y_zzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1567);

    #pragma omp simd aligned(pa_z, tr_y_zzzz_xxxxxx, tr_y_zzzz_xxxxxy, tr_y_zzzz_xxxxxz, tr_y_zzzz_xxxxyy, tr_y_zzzz_xxxxyz, tr_y_zzzz_xxxxzz, tr_y_zzzz_xxxyyy, tr_y_zzzz_xxxyyz, tr_y_zzzz_xxxyzz, tr_y_zzzz_xxxzzz, tr_y_zzzz_xxyyyy, tr_y_zzzz_xxyyyz, tr_y_zzzz_xxyyzz, tr_y_zzzz_xxyzzz, tr_y_zzzz_xxzzzz, tr_y_zzzz_xyyyyy, tr_y_zzzz_xyyyyz, tr_y_zzzz_xyyyzz, tr_y_zzzz_xyyzzz, tr_y_zzzz_xyzzzz, tr_y_zzzz_xzzzzz, tr_y_zzzz_yyyyyy, tr_y_zzzz_yyyyyz, tr_y_zzzz_yyyyzz, tr_y_zzzz_yyyzzz, tr_y_zzzz_yyzzzz, tr_y_zzzz_yzzzzz, tr_y_zzzz_zzzzzz, tr_y_zzzzz_xxxxx, tr_y_zzzzz_xxxxxx, tr_y_zzzzz_xxxxxy, tr_y_zzzzz_xxxxxz, tr_y_zzzzz_xxxxy, tr_y_zzzzz_xxxxyy, tr_y_zzzzz_xxxxyz, tr_y_zzzzz_xxxxz, tr_y_zzzzz_xxxxzz, tr_y_zzzzz_xxxyy, tr_y_zzzzz_xxxyyy, tr_y_zzzzz_xxxyyz, tr_y_zzzzz_xxxyz, tr_y_zzzzz_xxxyzz, tr_y_zzzzz_xxxzz, tr_y_zzzzz_xxxzzz, tr_y_zzzzz_xxyyy, tr_y_zzzzz_xxyyyy, tr_y_zzzzz_xxyyyz, tr_y_zzzzz_xxyyz, tr_y_zzzzz_xxyyzz, tr_y_zzzzz_xxyzz, tr_y_zzzzz_xxyzzz, tr_y_zzzzz_xxzzz, tr_y_zzzzz_xxzzzz, tr_y_zzzzz_xyyyy, tr_y_zzzzz_xyyyyy, tr_y_zzzzz_xyyyyz, tr_y_zzzzz_xyyyz, tr_y_zzzzz_xyyyzz, tr_y_zzzzz_xyyzz, tr_y_zzzzz_xyyzzz, tr_y_zzzzz_xyzzz, tr_y_zzzzz_xyzzzz, tr_y_zzzzz_xzzzz, tr_y_zzzzz_xzzzzz, tr_y_zzzzz_yyyyy, tr_y_zzzzz_yyyyyy, tr_y_zzzzz_yyyyyz, tr_y_zzzzz_yyyyz, tr_y_zzzzz_yyyyzz, tr_y_zzzzz_yyyzz, tr_y_zzzzz_yyyzzz, tr_y_zzzzz_yyzzz, tr_y_zzzzz_yyzzzz, tr_y_zzzzz_yzzzz, tr_y_zzzzz_yzzzzz, tr_y_zzzzz_zzzzz, tr_y_zzzzz_zzzzzz, tr_y_zzzzzz_xxxxxx, tr_y_zzzzzz_xxxxxy, tr_y_zzzzzz_xxxxxz, tr_y_zzzzzz_xxxxyy, tr_y_zzzzzz_xxxxyz, tr_y_zzzzzz_xxxxzz, tr_y_zzzzzz_xxxyyy, tr_y_zzzzzz_xxxyyz, tr_y_zzzzzz_xxxyzz, tr_y_zzzzzz_xxxzzz, tr_y_zzzzzz_xxyyyy, tr_y_zzzzzz_xxyyyz, tr_y_zzzzzz_xxyyzz, tr_y_zzzzzz_xxyzzz, tr_y_zzzzzz_xxzzzz, tr_y_zzzzzz_xyyyyy, tr_y_zzzzzz_xyyyyz, tr_y_zzzzzz_xyyyzz, tr_y_zzzzzz_xyyzzz, tr_y_zzzzzz_xyzzzz, tr_y_zzzzzz_xzzzzz, tr_y_zzzzzz_yyyyyy, tr_y_zzzzzz_yyyyyz, tr_y_zzzzzz_yyyyzz, tr_y_zzzzzz_yyyzzz, tr_y_zzzzzz_yyzzzz, tr_y_zzzzzz_yzzzzz, tr_y_zzzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzzz_xxxxxx[i] = 5.0 * tr_y_zzzz_xxxxxx[i] * fe_0 + tr_y_zzzzz_xxxxxx[i] * pa_z[i];

        tr_y_zzzzzz_xxxxxy[i] = 5.0 * tr_y_zzzz_xxxxxy[i] * fe_0 + tr_y_zzzzz_xxxxxy[i] * pa_z[i];

        tr_y_zzzzzz_xxxxxz[i] = 5.0 * tr_y_zzzz_xxxxxz[i] * fe_0 + tr_y_zzzzz_xxxxx[i] * fe_0 + tr_y_zzzzz_xxxxxz[i] * pa_z[i];

        tr_y_zzzzzz_xxxxyy[i] = 5.0 * tr_y_zzzz_xxxxyy[i] * fe_0 + tr_y_zzzzz_xxxxyy[i] * pa_z[i];

        tr_y_zzzzzz_xxxxyz[i] = 5.0 * tr_y_zzzz_xxxxyz[i] * fe_0 + tr_y_zzzzz_xxxxy[i] * fe_0 + tr_y_zzzzz_xxxxyz[i] * pa_z[i];

        tr_y_zzzzzz_xxxxzz[i] = 5.0 * tr_y_zzzz_xxxxzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xxxxz[i] * fe_0 + tr_y_zzzzz_xxxxzz[i] * pa_z[i];

        tr_y_zzzzzz_xxxyyy[i] = 5.0 * tr_y_zzzz_xxxyyy[i] * fe_0 + tr_y_zzzzz_xxxyyy[i] * pa_z[i];

        tr_y_zzzzzz_xxxyyz[i] = 5.0 * tr_y_zzzz_xxxyyz[i] * fe_0 + tr_y_zzzzz_xxxyy[i] * fe_0 + tr_y_zzzzz_xxxyyz[i] * pa_z[i];

        tr_y_zzzzzz_xxxyzz[i] = 5.0 * tr_y_zzzz_xxxyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xxxyz[i] * fe_0 + tr_y_zzzzz_xxxyzz[i] * pa_z[i];

        tr_y_zzzzzz_xxxzzz[i] = 5.0 * tr_y_zzzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_xxxzz[i] * fe_0 + tr_y_zzzzz_xxxzzz[i] * pa_z[i];

        tr_y_zzzzzz_xxyyyy[i] = 5.0 * tr_y_zzzz_xxyyyy[i] * fe_0 + tr_y_zzzzz_xxyyyy[i] * pa_z[i];

        tr_y_zzzzzz_xxyyyz[i] = 5.0 * tr_y_zzzz_xxyyyz[i] * fe_0 + tr_y_zzzzz_xxyyy[i] * fe_0 + tr_y_zzzzz_xxyyyz[i] * pa_z[i];

        tr_y_zzzzzz_xxyyzz[i] = 5.0 * tr_y_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xxyyz[i] * fe_0 + tr_y_zzzzz_xxyyzz[i] * pa_z[i];

        tr_y_zzzzzz_xxyzzz[i] = 5.0 * tr_y_zzzz_xxyzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_xxyzz[i] * fe_0 + tr_y_zzzzz_xxyzzz[i] * pa_z[i];

        tr_y_zzzzzz_xxzzzz[i] = 5.0 * tr_y_zzzz_xxzzzz[i] * fe_0 + 4.0 * tr_y_zzzzz_xxzzz[i] * fe_0 + tr_y_zzzzz_xxzzzz[i] * pa_z[i];

        tr_y_zzzzzz_xyyyyy[i] = 5.0 * tr_y_zzzz_xyyyyy[i] * fe_0 + tr_y_zzzzz_xyyyyy[i] * pa_z[i];

        tr_y_zzzzzz_xyyyyz[i] = 5.0 * tr_y_zzzz_xyyyyz[i] * fe_0 + tr_y_zzzzz_xyyyy[i] * fe_0 + tr_y_zzzzz_xyyyyz[i] * pa_z[i];

        tr_y_zzzzzz_xyyyzz[i] = 5.0 * tr_y_zzzz_xyyyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xyyyz[i] * fe_0 + tr_y_zzzzz_xyyyzz[i] * pa_z[i];

        tr_y_zzzzzz_xyyzzz[i] = 5.0 * tr_y_zzzz_xyyzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_xyyzz[i] * fe_0 + tr_y_zzzzz_xyyzzz[i] * pa_z[i];

        tr_y_zzzzzz_xyzzzz[i] = 5.0 * tr_y_zzzz_xyzzzz[i] * fe_0 + 4.0 * tr_y_zzzzz_xyzzz[i] * fe_0 + tr_y_zzzzz_xyzzzz[i] * pa_z[i];

        tr_y_zzzzzz_xzzzzz[i] = 5.0 * tr_y_zzzz_xzzzzz[i] * fe_0 + 5.0 * tr_y_zzzzz_xzzzz[i] * fe_0 + tr_y_zzzzz_xzzzzz[i] * pa_z[i];

        tr_y_zzzzzz_yyyyyy[i] = 5.0 * tr_y_zzzz_yyyyyy[i] * fe_0 + tr_y_zzzzz_yyyyyy[i] * pa_z[i];

        tr_y_zzzzzz_yyyyyz[i] = 5.0 * tr_y_zzzz_yyyyyz[i] * fe_0 + tr_y_zzzzz_yyyyy[i] * fe_0 + tr_y_zzzzz_yyyyyz[i] * pa_z[i];

        tr_y_zzzzzz_yyyyzz[i] = 5.0 * tr_y_zzzz_yyyyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_yyyyz[i] * fe_0 + tr_y_zzzzz_yyyyzz[i] * pa_z[i];

        tr_y_zzzzzz_yyyzzz[i] = 5.0 * tr_y_zzzz_yyyzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_yyyzz[i] * fe_0 + tr_y_zzzzz_yyyzzz[i] * pa_z[i];

        tr_y_zzzzzz_yyzzzz[i] = 5.0 * tr_y_zzzz_yyzzzz[i] * fe_0 + 4.0 * tr_y_zzzzz_yyzzz[i] * fe_0 + tr_y_zzzzz_yyzzzz[i] * pa_z[i];

        tr_y_zzzzzz_yzzzzz[i] = 5.0 * tr_y_zzzz_yzzzzz[i] * fe_0 + 5.0 * tr_y_zzzzz_yzzzz[i] * fe_0 + tr_y_zzzzz_yzzzzz[i] * pa_z[i];

        tr_y_zzzzzz_zzzzzz[i] = 5.0 * tr_y_zzzz_zzzzzz[i] * fe_0 + 6.0 * tr_y_zzzzz_zzzzz[i] * fe_0 + tr_y_zzzzz_zzzzzz[i] * pa_z[i];
    }

    // Set up 1568-1596 components of targeted buffer : II

    auto tr_z_xxxxxx_xxxxxx = pbuffer.data(idx_dip_ii + 1568);

    auto tr_z_xxxxxx_xxxxxy = pbuffer.data(idx_dip_ii + 1569);

    auto tr_z_xxxxxx_xxxxxz = pbuffer.data(idx_dip_ii + 1570);

    auto tr_z_xxxxxx_xxxxyy = pbuffer.data(idx_dip_ii + 1571);

    auto tr_z_xxxxxx_xxxxyz = pbuffer.data(idx_dip_ii + 1572);

    auto tr_z_xxxxxx_xxxxzz = pbuffer.data(idx_dip_ii + 1573);

    auto tr_z_xxxxxx_xxxyyy = pbuffer.data(idx_dip_ii + 1574);

    auto tr_z_xxxxxx_xxxyyz = pbuffer.data(idx_dip_ii + 1575);

    auto tr_z_xxxxxx_xxxyzz = pbuffer.data(idx_dip_ii + 1576);

    auto tr_z_xxxxxx_xxxzzz = pbuffer.data(idx_dip_ii + 1577);

    auto tr_z_xxxxxx_xxyyyy = pbuffer.data(idx_dip_ii + 1578);

    auto tr_z_xxxxxx_xxyyyz = pbuffer.data(idx_dip_ii + 1579);

    auto tr_z_xxxxxx_xxyyzz = pbuffer.data(idx_dip_ii + 1580);

    auto tr_z_xxxxxx_xxyzzz = pbuffer.data(idx_dip_ii + 1581);

    auto tr_z_xxxxxx_xxzzzz = pbuffer.data(idx_dip_ii + 1582);

    auto tr_z_xxxxxx_xyyyyy = pbuffer.data(idx_dip_ii + 1583);

    auto tr_z_xxxxxx_xyyyyz = pbuffer.data(idx_dip_ii + 1584);

    auto tr_z_xxxxxx_xyyyzz = pbuffer.data(idx_dip_ii + 1585);

    auto tr_z_xxxxxx_xyyzzz = pbuffer.data(idx_dip_ii + 1586);

    auto tr_z_xxxxxx_xyzzzz = pbuffer.data(idx_dip_ii + 1587);

    auto tr_z_xxxxxx_xzzzzz = pbuffer.data(idx_dip_ii + 1588);

    auto tr_z_xxxxxx_yyyyyy = pbuffer.data(idx_dip_ii + 1589);

    auto tr_z_xxxxxx_yyyyyz = pbuffer.data(idx_dip_ii + 1590);

    auto tr_z_xxxxxx_yyyyzz = pbuffer.data(idx_dip_ii + 1591);

    auto tr_z_xxxxxx_yyyzzz = pbuffer.data(idx_dip_ii + 1592);

    auto tr_z_xxxxxx_yyzzzz = pbuffer.data(idx_dip_ii + 1593);

    auto tr_z_xxxxxx_yzzzzz = pbuffer.data(idx_dip_ii + 1594);

    auto tr_z_xxxxxx_zzzzzz = pbuffer.data(idx_dip_ii + 1595);

    #pragma omp simd aligned(pa_x, tr_z_xxxx_xxxxxx, tr_z_xxxx_xxxxxy, tr_z_xxxx_xxxxxz, tr_z_xxxx_xxxxyy, tr_z_xxxx_xxxxyz, tr_z_xxxx_xxxxzz, tr_z_xxxx_xxxyyy, tr_z_xxxx_xxxyyz, tr_z_xxxx_xxxyzz, tr_z_xxxx_xxxzzz, tr_z_xxxx_xxyyyy, tr_z_xxxx_xxyyyz, tr_z_xxxx_xxyyzz, tr_z_xxxx_xxyzzz, tr_z_xxxx_xxzzzz, tr_z_xxxx_xyyyyy, tr_z_xxxx_xyyyyz, tr_z_xxxx_xyyyzz, tr_z_xxxx_xyyzzz, tr_z_xxxx_xyzzzz, tr_z_xxxx_xzzzzz, tr_z_xxxx_yyyyyy, tr_z_xxxx_yyyyyz, tr_z_xxxx_yyyyzz, tr_z_xxxx_yyyzzz, tr_z_xxxx_yyzzzz, tr_z_xxxx_yzzzzz, tr_z_xxxx_zzzzzz, tr_z_xxxxx_xxxxx, tr_z_xxxxx_xxxxxx, tr_z_xxxxx_xxxxxy, tr_z_xxxxx_xxxxxz, tr_z_xxxxx_xxxxy, tr_z_xxxxx_xxxxyy, tr_z_xxxxx_xxxxyz, tr_z_xxxxx_xxxxz, tr_z_xxxxx_xxxxzz, tr_z_xxxxx_xxxyy, tr_z_xxxxx_xxxyyy, tr_z_xxxxx_xxxyyz, tr_z_xxxxx_xxxyz, tr_z_xxxxx_xxxyzz, tr_z_xxxxx_xxxzz, tr_z_xxxxx_xxxzzz, tr_z_xxxxx_xxyyy, tr_z_xxxxx_xxyyyy, tr_z_xxxxx_xxyyyz, tr_z_xxxxx_xxyyz, tr_z_xxxxx_xxyyzz, tr_z_xxxxx_xxyzz, tr_z_xxxxx_xxyzzz, tr_z_xxxxx_xxzzz, tr_z_xxxxx_xxzzzz, tr_z_xxxxx_xyyyy, tr_z_xxxxx_xyyyyy, tr_z_xxxxx_xyyyyz, tr_z_xxxxx_xyyyz, tr_z_xxxxx_xyyyzz, tr_z_xxxxx_xyyzz, tr_z_xxxxx_xyyzzz, tr_z_xxxxx_xyzzz, tr_z_xxxxx_xyzzzz, tr_z_xxxxx_xzzzz, tr_z_xxxxx_xzzzzz, tr_z_xxxxx_yyyyy, tr_z_xxxxx_yyyyyy, tr_z_xxxxx_yyyyyz, tr_z_xxxxx_yyyyz, tr_z_xxxxx_yyyyzz, tr_z_xxxxx_yyyzz, tr_z_xxxxx_yyyzzz, tr_z_xxxxx_yyzzz, tr_z_xxxxx_yyzzzz, tr_z_xxxxx_yzzzz, tr_z_xxxxx_yzzzzz, tr_z_xxxxx_zzzzz, tr_z_xxxxx_zzzzzz, tr_z_xxxxxx_xxxxxx, tr_z_xxxxxx_xxxxxy, tr_z_xxxxxx_xxxxxz, tr_z_xxxxxx_xxxxyy, tr_z_xxxxxx_xxxxyz, tr_z_xxxxxx_xxxxzz, tr_z_xxxxxx_xxxyyy, tr_z_xxxxxx_xxxyyz, tr_z_xxxxxx_xxxyzz, tr_z_xxxxxx_xxxzzz, tr_z_xxxxxx_xxyyyy, tr_z_xxxxxx_xxyyyz, tr_z_xxxxxx_xxyyzz, tr_z_xxxxxx_xxyzzz, tr_z_xxxxxx_xxzzzz, tr_z_xxxxxx_xyyyyy, tr_z_xxxxxx_xyyyyz, tr_z_xxxxxx_xyyyzz, tr_z_xxxxxx_xyyzzz, tr_z_xxxxxx_xyzzzz, tr_z_xxxxxx_xzzzzz, tr_z_xxxxxx_yyyyyy, tr_z_xxxxxx_yyyyyz, tr_z_xxxxxx_yyyyzz, tr_z_xxxxxx_yyyzzz, tr_z_xxxxxx_yyzzzz, tr_z_xxxxxx_yzzzzz, tr_z_xxxxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxx_xxxxxx[i] = 5.0 * tr_z_xxxx_xxxxxx[i] * fe_0 + 6.0 * tr_z_xxxxx_xxxxx[i] * fe_0 + tr_z_xxxxx_xxxxxx[i] * pa_x[i];

        tr_z_xxxxxx_xxxxxy[i] = 5.0 * tr_z_xxxx_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxxxx_xxxxy[i] * fe_0 + tr_z_xxxxx_xxxxxy[i] * pa_x[i];

        tr_z_xxxxxx_xxxxxz[i] = 5.0 * tr_z_xxxx_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxxxx_xxxxz[i] * fe_0 + tr_z_xxxxx_xxxxxz[i] * pa_x[i];

        tr_z_xxxxxx_xxxxyy[i] = 5.0 * tr_z_xxxx_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxxxx_xxxyy[i] * fe_0 + tr_z_xxxxx_xxxxyy[i] * pa_x[i];

        tr_z_xxxxxx_xxxxyz[i] = 5.0 * tr_z_xxxx_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxxxx_xxxyz[i] * fe_0 + tr_z_xxxxx_xxxxyz[i] * pa_x[i];

        tr_z_xxxxxx_xxxxzz[i] = 5.0 * tr_z_xxxx_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxxxx_xxxzz[i] * fe_0 + tr_z_xxxxx_xxxxzz[i] * pa_x[i];

        tr_z_xxxxxx_xxxyyy[i] = 5.0 * tr_z_xxxx_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxxxx_xxyyy[i] * fe_0 + tr_z_xxxxx_xxxyyy[i] * pa_x[i];

        tr_z_xxxxxx_xxxyyz[i] = 5.0 * tr_z_xxxx_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxxxx_xxyyz[i] * fe_0 + tr_z_xxxxx_xxxyyz[i] * pa_x[i];

        tr_z_xxxxxx_xxxyzz[i] = 5.0 * tr_z_xxxx_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxxxx_xxyzz[i] * fe_0 + tr_z_xxxxx_xxxyzz[i] * pa_x[i];

        tr_z_xxxxxx_xxxzzz[i] = 5.0 * tr_z_xxxx_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxxxx_xxzzz[i] * fe_0 + tr_z_xxxxx_xxxzzz[i] * pa_x[i];

        tr_z_xxxxxx_xxyyyy[i] = 5.0 * tr_z_xxxx_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxxxx_xyyyy[i] * fe_0 + tr_z_xxxxx_xxyyyy[i] * pa_x[i];

        tr_z_xxxxxx_xxyyyz[i] = 5.0 * tr_z_xxxx_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxxxx_xyyyz[i] * fe_0 + tr_z_xxxxx_xxyyyz[i] * pa_x[i];

        tr_z_xxxxxx_xxyyzz[i] = 5.0 * tr_z_xxxx_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxxxx_xyyzz[i] * fe_0 + tr_z_xxxxx_xxyyzz[i] * pa_x[i];

        tr_z_xxxxxx_xxyzzz[i] = 5.0 * tr_z_xxxx_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxxxx_xyzzz[i] * fe_0 + tr_z_xxxxx_xxyzzz[i] * pa_x[i];

        tr_z_xxxxxx_xxzzzz[i] = 5.0 * tr_z_xxxx_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxxxx_xzzzz[i] * fe_0 + tr_z_xxxxx_xxzzzz[i] * pa_x[i];

        tr_z_xxxxxx_xyyyyy[i] = 5.0 * tr_z_xxxx_xyyyyy[i] * fe_0 + tr_z_xxxxx_yyyyy[i] * fe_0 + tr_z_xxxxx_xyyyyy[i] * pa_x[i];

        tr_z_xxxxxx_xyyyyz[i] = 5.0 * tr_z_xxxx_xyyyyz[i] * fe_0 + tr_z_xxxxx_yyyyz[i] * fe_0 + tr_z_xxxxx_xyyyyz[i] * pa_x[i];

        tr_z_xxxxxx_xyyyzz[i] = 5.0 * tr_z_xxxx_xyyyzz[i] * fe_0 + tr_z_xxxxx_yyyzz[i] * fe_0 + tr_z_xxxxx_xyyyzz[i] * pa_x[i];

        tr_z_xxxxxx_xyyzzz[i] = 5.0 * tr_z_xxxx_xyyzzz[i] * fe_0 + tr_z_xxxxx_yyzzz[i] * fe_0 + tr_z_xxxxx_xyyzzz[i] * pa_x[i];

        tr_z_xxxxxx_xyzzzz[i] = 5.0 * tr_z_xxxx_xyzzzz[i] * fe_0 + tr_z_xxxxx_yzzzz[i] * fe_0 + tr_z_xxxxx_xyzzzz[i] * pa_x[i];

        tr_z_xxxxxx_xzzzzz[i] = 5.0 * tr_z_xxxx_xzzzzz[i] * fe_0 + tr_z_xxxxx_zzzzz[i] * fe_0 + tr_z_xxxxx_xzzzzz[i] * pa_x[i];

        tr_z_xxxxxx_yyyyyy[i] = 5.0 * tr_z_xxxx_yyyyyy[i] * fe_0 + tr_z_xxxxx_yyyyyy[i] * pa_x[i];

        tr_z_xxxxxx_yyyyyz[i] = 5.0 * tr_z_xxxx_yyyyyz[i] * fe_0 + tr_z_xxxxx_yyyyyz[i] * pa_x[i];

        tr_z_xxxxxx_yyyyzz[i] = 5.0 * tr_z_xxxx_yyyyzz[i] * fe_0 + tr_z_xxxxx_yyyyzz[i] * pa_x[i];

        tr_z_xxxxxx_yyyzzz[i] = 5.0 * tr_z_xxxx_yyyzzz[i] * fe_0 + tr_z_xxxxx_yyyzzz[i] * pa_x[i];

        tr_z_xxxxxx_yyzzzz[i] = 5.0 * tr_z_xxxx_yyzzzz[i] * fe_0 + tr_z_xxxxx_yyzzzz[i] * pa_x[i];

        tr_z_xxxxxx_yzzzzz[i] = 5.0 * tr_z_xxxx_yzzzzz[i] * fe_0 + tr_z_xxxxx_yzzzzz[i] * pa_x[i];

        tr_z_xxxxxx_zzzzzz[i] = 5.0 * tr_z_xxxx_zzzzzz[i] * fe_0 + tr_z_xxxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 1596-1624 components of targeted buffer : II

    auto tr_z_xxxxxy_xxxxxx = pbuffer.data(idx_dip_ii + 1596);

    auto tr_z_xxxxxy_xxxxxy = pbuffer.data(idx_dip_ii + 1597);

    auto tr_z_xxxxxy_xxxxxz = pbuffer.data(idx_dip_ii + 1598);

    auto tr_z_xxxxxy_xxxxyy = pbuffer.data(idx_dip_ii + 1599);

    auto tr_z_xxxxxy_xxxxyz = pbuffer.data(idx_dip_ii + 1600);

    auto tr_z_xxxxxy_xxxxzz = pbuffer.data(idx_dip_ii + 1601);

    auto tr_z_xxxxxy_xxxyyy = pbuffer.data(idx_dip_ii + 1602);

    auto tr_z_xxxxxy_xxxyyz = pbuffer.data(idx_dip_ii + 1603);

    auto tr_z_xxxxxy_xxxyzz = pbuffer.data(idx_dip_ii + 1604);

    auto tr_z_xxxxxy_xxxzzz = pbuffer.data(idx_dip_ii + 1605);

    auto tr_z_xxxxxy_xxyyyy = pbuffer.data(idx_dip_ii + 1606);

    auto tr_z_xxxxxy_xxyyyz = pbuffer.data(idx_dip_ii + 1607);

    auto tr_z_xxxxxy_xxyyzz = pbuffer.data(idx_dip_ii + 1608);

    auto tr_z_xxxxxy_xxyzzz = pbuffer.data(idx_dip_ii + 1609);

    auto tr_z_xxxxxy_xxzzzz = pbuffer.data(idx_dip_ii + 1610);

    auto tr_z_xxxxxy_xyyyyy = pbuffer.data(idx_dip_ii + 1611);

    auto tr_z_xxxxxy_xyyyyz = pbuffer.data(idx_dip_ii + 1612);

    auto tr_z_xxxxxy_xyyyzz = pbuffer.data(idx_dip_ii + 1613);

    auto tr_z_xxxxxy_xyyzzz = pbuffer.data(idx_dip_ii + 1614);

    auto tr_z_xxxxxy_xyzzzz = pbuffer.data(idx_dip_ii + 1615);

    auto tr_z_xxxxxy_xzzzzz = pbuffer.data(idx_dip_ii + 1616);

    auto tr_z_xxxxxy_yyyyyy = pbuffer.data(idx_dip_ii + 1617);

    auto tr_z_xxxxxy_yyyyyz = pbuffer.data(idx_dip_ii + 1618);

    auto tr_z_xxxxxy_yyyyzz = pbuffer.data(idx_dip_ii + 1619);

    auto tr_z_xxxxxy_yyyzzz = pbuffer.data(idx_dip_ii + 1620);

    auto tr_z_xxxxxy_yyzzzz = pbuffer.data(idx_dip_ii + 1621);

    auto tr_z_xxxxxy_yzzzzz = pbuffer.data(idx_dip_ii + 1622);

    auto tr_z_xxxxxy_zzzzzz = pbuffer.data(idx_dip_ii + 1623);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxx_xxxxx, tr_z_xxxxx_xxxxxx, tr_z_xxxxx_xxxxxy, tr_z_xxxxx_xxxxxz, tr_z_xxxxx_xxxxy, tr_z_xxxxx_xxxxyy, tr_z_xxxxx_xxxxyz, tr_z_xxxxx_xxxxz, tr_z_xxxxx_xxxxzz, tr_z_xxxxx_xxxyy, tr_z_xxxxx_xxxyyy, tr_z_xxxxx_xxxyyz, tr_z_xxxxx_xxxyz, tr_z_xxxxx_xxxyzz, tr_z_xxxxx_xxxzz, tr_z_xxxxx_xxxzzz, tr_z_xxxxx_xxyyy, tr_z_xxxxx_xxyyyy, tr_z_xxxxx_xxyyyz, tr_z_xxxxx_xxyyz, tr_z_xxxxx_xxyyzz, tr_z_xxxxx_xxyzz, tr_z_xxxxx_xxyzzz, tr_z_xxxxx_xxzzz, tr_z_xxxxx_xxzzzz, tr_z_xxxxx_xyyyy, tr_z_xxxxx_xyyyyy, tr_z_xxxxx_xyyyyz, tr_z_xxxxx_xyyyz, tr_z_xxxxx_xyyyzz, tr_z_xxxxx_xyyzz, tr_z_xxxxx_xyyzzz, tr_z_xxxxx_xyzzz, tr_z_xxxxx_xyzzzz, tr_z_xxxxx_xzzzz, tr_z_xxxxx_xzzzzz, tr_z_xxxxx_zzzzzz, tr_z_xxxxxy_xxxxxx, tr_z_xxxxxy_xxxxxy, tr_z_xxxxxy_xxxxxz, tr_z_xxxxxy_xxxxyy, tr_z_xxxxxy_xxxxyz, tr_z_xxxxxy_xxxxzz, tr_z_xxxxxy_xxxyyy, tr_z_xxxxxy_xxxyyz, tr_z_xxxxxy_xxxyzz, tr_z_xxxxxy_xxxzzz, tr_z_xxxxxy_xxyyyy, tr_z_xxxxxy_xxyyyz, tr_z_xxxxxy_xxyyzz, tr_z_xxxxxy_xxyzzz, tr_z_xxxxxy_xxzzzz, tr_z_xxxxxy_xyyyyy, tr_z_xxxxxy_xyyyyz, tr_z_xxxxxy_xyyyzz, tr_z_xxxxxy_xyyzzz, tr_z_xxxxxy_xyzzzz, tr_z_xxxxxy_xzzzzz, tr_z_xxxxxy_yyyyyy, tr_z_xxxxxy_yyyyyz, tr_z_xxxxxy_yyyyzz, tr_z_xxxxxy_yyyzzz, tr_z_xxxxxy_yyzzzz, tr_z_xxxxxy_yzzzzz, tr_z_xxxxxy_zzzzzz, tr_z_xxxxy_yyyyyy, tr_z_xxxxy_yyyyyz, tr_z_xxxxy_yyyyzz, tr_z_xxxxy_yyyzzz, tr_z_xxxxy_yyzzzz, tr_z_xxxxy_yzzzzz, tr_z_xxxy_yyyyyy, tr_z_xxxy_yyyyyz, tr_z_xxxy_yyyyzz, tr_z_xxxy_yyyzzz, tr_z_xxxy_yyzzzz, tr_z_xxxy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxy_xxxxxx[i] = tr_z_xxxxx_xxxxxx[i] * pa_y[i];

        tr_z_xxxxxy_xxxxxy[i] = tr_z_xxxxx_xxxxx[i] * fe_0 + tr_z_xxxxx_xxxxxy[i] * pa_y[i];

        tr_z_xxxxxy_xxxxxz[i] = tr_z_xxxxx_xxxxxz[i] * pa_y[i];

        tr_z_xxxxxy_xxxxyy[i] = 2.0 * tr_z_xxxxx_xxxxy[i] * fe_0 + tr_z_xxxxx_xxxxyy[i] * pa_y[i];

        tr_z_xxxxxy_xxxxyz[i] = tr_z_xxxxx_xxxxz[i] * fe_0 + tr_z_xxxxx_xxxxyz[i] * pa_y[i];

        tr_z_xxxxxy_xxxxzz[i] = tr_z_xxxxx_xxxxzz[i] * pa_y[i];

        tr_z_xxxxxy_xxxyyy[i] = 3.0 * tr_z_xxxxx_xxxyy[i] * fe_0 + tr_z_xxxxx_xxxyyy[i] * pa_y[i];

        tr_z_xxxxxy_xxxyyz[i] = 2.0 * tr_z_xxxxx_xxxyz[i] * fe_0 + tr_z_xxxxx_xxxyyz[i] * pa_y[i];

        tr_z_xxxxxy_xxxyzz[i] = tr_z_xxxxx_xxxzz[i] * fe_0 + tr_z_xxxxx_xxxyzz[i] * pa_y[i];

        tr_z_xxxxxy_xxxzzz[i] = tr_z_xxxxx_xxxzzz[i] * pa_y[i];

        tr_z_xxxxxy_xxyyyy[i] = 4.0 * tr_z_xxxxx_xxyyy[i] * fe_0 + tr_z_xxxxx_xxyyyy[i] * pa_y[i];

        tr_z_xxxxxy_xxyyyz[i] = 3.0 * tr_z_xxxxx_xxyyz[i] * fe_0 + tr_z_xxxxx_xxyyyz[i] * pa_y[i];

        tr_z_xxxxxy_xxyyzz[i] = 2.0 * tr_z_xxxxx_xxyzz[i] * fe_0 + tr_z_xxxxx_xxyyzz[i] * pa_y[i];

        tr_z_xxxxxy_xxyzzz[i] = tr_z_xxxxx_xxzzz[i] * fe_0 + tr_z_xxxxx_xxyzzz[i] * pa_y[i];

        tr_z_xxxxxy_xxzzzz[i] = tr_z_xxxxx_xxzzzz[i] * pa_y[i];

        tr_z_xxxxxy_xyyyyy[i] = 5.0 * tr_z_xxxxx_xyyyy[i] * fe_0 + tr_z_xxxxx_xyyyyy[i] * pa_y[i];

        tr_z_xxxxxy_xyyyyz[i] = 4.0 * tr_z_xxxxx_xyyyz[i] * fe_0 + tr_z_xxxxx_xyyyyz[i] * pa_y[i];

        tr_z_xxxxxy_xyyyzz[i] = 3.0 * tr_z_xxxxx_xyyzz[i] * fe_0 + tr_z_xxxxx_xyyyzz[i] * pa_y[i];

        tr_z_xxxxxy_xyyzzz[i] = 2.0 * tr_z_xxxxx_xyzzz[i] * fe_0 + tr_z_xxxxx_xyyzzz[i] * pa_y[i];

        tr_z_xxxxxy_xyzzzz[i] = tr_z_xxxxx_xzzzz[i] * fe_0 + tr_z_xxxxx_xyzzzz[i] * pa_y[i];

        tr_z_xxxxxy_xzzzzz[i] = tr_z_xxxxx_xzzzzz[i] * pa_y[i];

        tr_z_xxxxxy_yyyyyy[i] = 4.0 * tr_z_xxxy_yyyyyy[i] * fe_0 + tr_z_xxxxy_yyyyyy[i] * pa_x[i];

        tr_z_xxxxxy_yyyyyz[i] = 4.0 * tr_z_xxxy_yyyyyz[i] * fe_0 + tr_z_xxxxy_yyyyyz[i] * pa_x[i];

        tr_z_xxxxxy_yyyyzz[i] = 4.0 * tr_z_xxxy_yyyyzz[i] * fe_0 + tr_z_xxxxy_yyyyzz[i] * pa_x[i];

        tr_z_xxxxxy_yyyzzz[i] = 4.0 * tr_z_xxxy_yyyzzz[i] * fe_0 + tr_z_xxxxy_yyyzzz[i] * pa_x[i];

        tr_z_xxxxxy_yyzzzz[i] = 4.0 * tr_z_xxxy_yyzzzz[i] * fe_0 + tr_z_xxxxy_yyzzzz[i] * pa_x[i];

        tr_z_xxxxxy_yzzzzz[i] = 4.0 * tr_z_xxxy_yzzzzz[i] * fe_0 + tr_z_xxxxy_yzzzzz[i] * pa_x[i];

        tr_z_xxxxxy_zzzzzz[i] = tr_z_xxxxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 1624-1652 components of targeted buffer : II

    auto tr_z_xxxxxz_xxxxxx = pbuffer.data(idx_dip_ii + 1624);

    auto tr_z_xxxxxz_xxxxxy = pbuffer.data(idx_dip_ii + 1625);

    auto tr_z_xxxxxz_xxxxxz = pbuffer.data(idx_dip_ii + 1626);

    auto tr_z_xxxxxz_xxxxyy = pbuffer.data(idx_dip_ii + 1627);

    auto tr_z_xxxxxz_xxxxyz = pbuffer.data(idx_dip_ii + 1628);

    auto tr_z_xxxxxz_xxxxzz = pbuffer.data(idx_dip_ii + 1629);

    auto tr_z_xxxxxz_xxxyyy = pbuffer.data(idx_dip_ii + 1630);

    auto tr_z_xxxxxz_xxxyyz = pbuffer.data(idx_dip_ii + 1631);

    auto tr_z_xxxxxz_xxxyzz = pbuffer.data(idx_dip_ii + 1632);

    auto tr_z_xxxxxz_xxxzzz = pbuffer.data(idx_dip_ii + 1633);

    auto tr_z_xxxxxz_xxyyyy = pbuffer.data(idx_dip_ii + 1634);

    auto tr_z_xxxxxz_xxyyyz = pbuffer.data(idx_dip_ii + 1635);

    auto tr_z_xxxxxz_xxyyzz = pbuffer.data(idx_dip_ii + 1636);

    auto tr_z_xxxxxz_xxyzzz = pbuffer.data(idx_dip_ii + 1637);

    auto tr_z_xxxxxz_xxzzzz = pbuffer.data(idx_dip_ii + 1638);

    auto tr_z_xxxxxz_xyyyyy = pbuffer.data(idx_dip_ii + 1639);

    auto tr_z_xxxxxz_xyyyyz = pbuffer.data(idx_dip_ii + 1640);

    auto tr_z_xxxxxz_xyyyzz = pbuffer.data(idx_dip_ii + 1641);

    auto tr_z_xxxxxz_xyyzzz = pbuffer.data(idx_dip_ii + 1642);

    auto tr_z_xxxxxz_xyzzzz = pbuffer.data(idx_dip_ii + 1643);

    auto tr_z_xxxxxz_xzzzzz = pbuffer.data(idx_dip_ii + 1644);

    auto tr_z_xxxxxz_yyyyyy = pbuffer.data(idx_dip_ii + 1645);

    auto tr_z_xxxxxz_yyyyyz = pbuffer.data(idx_dip_ii + 1646);

    auto tr_z_xxxxxz_yyyyzz = pbuffer.data(idx_dip_ii + 1647);

    auto tr_z_xxxxxz_yyyzzz = pbuffer.data(idx_dip_ii + 1648);

    auto tr_z_xxxxxz_yyzzzz = pbuffer.data(idx_dip_ii + 1649);

    auto tr_z_xxxxxz_yzzzzz = pbuffer.data(idx_dip_ii + 1650);

    auto tr_z_xxxxxz_zzzzzz = pbuffer.data(idx_dip_ii + 1651);

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xxxxx_xxxxxx, tr_z_xxxxx_xxxxxy, tr_z_xxxxx_xxxxyy, tr_z_xxxxx_xxxyyy, tr_z_xxxxx_xxyyyy, tr_z_xxxxx_xyyyyy, tr_z_xxxxxz_xxxxxx, tr_z_xxxxxz_xxxxxy, tr_z_xxxxxz_xxxxxz, tr_z_xxxxxz_xxxxyy, tr_z_xxxxxz_xxxxyz, tr_z_xxxxxz_xxxxzz, tr_z_xxxxxz_xxxyyy, tr_z_xxxxxz_xxxyyz, tr_z_xxxxxz_xxxyzz, tr_z_xxxxxz_xxxzzz, tr_z_xxxxxz_xxyyyy, tr_z_xxxxxz_xxyyyz, tr_z_xxxxxz_xxyyzz, tr_z_xxxxxz_xxyzzz, tr_z_xxxxxz_xxzzzz, tr_z_xxxxxz_xyyyyy, tr_z_xxxxxz_xyyyyz, tr_z_xxxxxz_xyyyzz, tr_z_xxxxxz_xyyzzz, tr_z_xxxxxz_xyzzzz, tr_z_xxxxxz_xzzzzz, tr_z_xxxxxz_yyyyyy, tr_z_xxxxxz_yyyyyz, tr_z_xxxxxz_yyyyzz, tr_z_xxxxxz_yyyzzz, tr_z_xxxxxz_yyzzzz, tr_z_xxxxxz_yzzzzz, tr_z_xxxxxz_zzzzzz, tr_z_xxxxz_xxxxxz, tr_z_xxxxz_xxxxyz, tr_z_xxxxz_xxxxz, tr_z_xxxxz_xxxxzz, tr_z_xxxxz_xxxyyz, tr_z_xxxxz_xxxyz, tr_z_xxxxz_xxxyzz, tr_z_xxxxz_xxxzz, tr_z_xxxxz_xxxzzz, tr_z_xxxxz_xxyyyz, tr_z_xxxxz_xxyyz, tr_z_xxxxz_xxyyzz, tr_z_xxxxz_xxyzz, tr_z_xxxxz_xxyzzz, tr_z_xxxxz_xxzzz, tr_z_xxxxz_xxzzzz, tr_z_xxxxz_xyyyyz, tr_z_xxxxz_xyyyz, tr_z_xxxxz_xyyyzz, tr_z_xxxxz_xyyzz, tr_z_xxxxz_xyyzzz, tr_z_xxxxz_xyzzz, tr_z_xxxxz_xyzzzz, tr_z_xxxxz_xzzzz, tr_z_xxxxz_xzzzzz, tr_z_xxxxz_yyyyyy, tr_z_xxxxz_yyyyyz, tr_z_xxxxz_yyyyz, tr_z_xxxxz_yyyyzz, tr_z_xxxxz_yyyzz, tr_z_xxxxz_yyyzzz, tr_z_xxxxz_yyzzz, tr_z_xxxxz_yyzzzz, tr_z_xxxxz_yzzzz, tr_z_xxxxz_yzzzzz, tr_z_xxxxz_zzzzz, tr_z_xxxxz_zzzzzz, tr_z_xxxz_xxxxxz, tr_z_xxxz_xxxxyz, tr_z_xxxz_xxxxzz, tr_z_xxxz_xxxyyz, tr_z_xxxz_xxxyzz, tr_z_xxxz_xxxzzz, tr_z_xxxz_xxyyyz, tr_z_xxxz_xxyyzz, tr_z_xxxz_xxyzzz, tr_z_xxxz_xxzzzz, tr_z_xxxz_xyyyyz, tr_z_xxxz_xyyyzz, tr_z_xxxz_xyyzzz, tr_z_xxxz_xyzzzz, tr_z_xxxz_xzzzzz, tr_z_xxxz_yyyyyy, tr_z_xxxz_yyyyyz, tr_z_xxxz_yyyyzz, tr_z_xxxz_yyyzzz, tr_z_xxxz_yyzzzz, tr_z_xxxz_yzzzzz, tr_z_xxxz_zzzzzz, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxy, ts_xxxxx_xxxxyy, ts_xxxxx_xxxyyy, ts_xxxxx_xxyyyy, ts_xxxxx_xyyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxz_xxxxxx[i] = ts_xxxxx_xxxxxx[i] * fe_0 + tr_z_xxxxx_xxxxxx[i] * pa_z[i];

        tr_z_xxxxxz_xxxxxy[i] = ts_xxxxx_xxxxxy[i] * fe_0 + tr_z_xxxxx_xxxxxy[i] * pa_z[i];

        tr_z_xxxxxz_xxxxxz[i] = 4.0 * tr_z_xxxz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxxxz_xxxxz[i] * fe_0 + tr_z_xxxxz_xxxxxz[i] * pa_x[i];

        tr_z_xxxxxz_xxxxyy[i] = ts_xxxxx_xxxxyy[i] * fe_0 + tr_z_xxxxx_xxxxyy[i] * pa_z[i];

        tr_z_xxxxxz_xxxxyz[i] = 4.0 * tr_z_xxxz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxxxz_xxxyz[i] * fe_0 + tr_z_xxxxz_xxxxyz[i] * pa_x[i];

        tr_z_xxxxxz_xxxxzz[i] = 4.0 * tr_z_xxxz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxxxz_xxxzz[i] * fe_0 + tr_z_xxxxz_xxxxzz[i] * pa_x[i];

        tr_z_xxxxxz_xxxyyy[i] = ts_xxxxx_xxxyyy[i] * fe_0 + tr_z_xxxxx_xxxyyy[i] * pa_z[i];

        tr_z_xxxxxz_xxxyyz[i] = 4.0 * tr_z_xxxz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxxxz_xxyyz[i] * fe_0 + tr_z_xxxxz_xxxyyz[i] * pa_x[i];

        tr_z_xxxxxz_xxxyzz[i] = 4.0 * tr_z_xxxz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxxxz_xxyzz[i] * fe_0 + tr_z_xxxxz_xxxyzz[i] * pa_x[i];

        tr_z_xxxxxz_xxxzzz[i] = 4.0 * tr_z_xxxz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxxxz_xxzzz[i] * fe_0 + tr_z_xxxxz_xxxzzz[i] * pa_x[i];

        tr_z_xxxxxz_xxyyyy[i] = ts_xxxxx_xxyyyy[i] * fe_0 + tr_z_xxxxx_xxyyyy[i] * pa_z[i];

        tr_z_xxxxxz_xxyyyz[i] = 4.0 * tr_z_xxxz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxxxz_xyyyz[i] * fe_0 + tr_z_xxxxz_xxyyyz[i] * pa_x[i];

        tr_z_xxxxxz_xxyyzz[i] = 4.0 * tr_z_xxxz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxxxz_xyyzz[i] * fe_0 + tr_z_xxxxz_xxyyzz[i] * pa_x[i];

        tr_z_xxxxxz_xxyzzz[i] = 4.0 * tr_z_xxxz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxxxz_xyzzz[i] * fe_0 + tr_z_xxxxz_xxyzzz[i] * pa_x[i];

        tr_z_xxxxxz_xxzzzz[i] = 4.0 * tr_z_xxxz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxxxz_xzzzz[i] * fe_0 + tr_z_xxxxz_xxzzzz[i] * pa_x[i];

        tr_z_xxxxxz_xyyyyy[i] = ts_xxxxx_xyyyyy[i] * fe_0 + tr_z_xxxxx_xyyyyy[i] * pa_z[i];

        tr_z_xxxxxz_xyyyyz[i] = 4.0 * tr_z_xxxz_xyyyyz[i] * fe_0 + tr_z_xxxxz_yyyyz[i] * fe_0 + tr_z_xxxxz_xyyyyz[i] * pa_x[i];

        tr_z_xxxxxz_xyyyzz[i] = 4.0 * tr_z_xxxz_xyyyzz[i] * fe_0 + tr_z_xxxxz_yyyzz[i] * fe_0 + tr_z_xxxxz_xyyyzz[i] * pa_x[i];

        tr_z_xxxxxz_xyyzzz[i] = 4.0 * tr_z_xxxz_xyyzzz[i] * fe_0 + tr_z_xxxxz_yyzzz[i] * fe_0 + tr_z_xxxxz_xyyzzz[i] * pa_x[i];

        tr_z_xxxxxz_xyzzzz[i] = 4.0 * tr_z_xxxz_xyzzzz[i] * fe_0 + tr_z_xxxxz_yzzzz[i] * fe_0 + tr_z_xxxxz_xyzzzz[i] * pa_x[i];

        tr_z_xxxxxz_xzzzzz[i] = 4.0 * tr_z_xxxz_xzzzzz[i] * fe_0 + tr_z_xxxxz_zzzzz[i] * fe_0 + tr_z_xxxxz_xzzzzz[i] * pa_x[i];

        tr_z_xxxxxz_yyyyyy[i] = 4.0 * tr_z_xxxz_yyyyyy[i] * fe_0 + tr_z_xxxxz_yyyyyy[i] * pa_x[i];

        tr_z_xxxxxz_yyyyyz[i] = 4.0 * tr_z_xxxz_yyyyyz[i] * fe_0 + tr_z_xxxxz_yyyyyz[i] * pa_x[i];

        tr_z_xxxxxz_yyyyzz[i] = 4.0 * tr_z_xxxz_yyyyzz[i] * fe_0 + tr_z_xxxxz_yyyyzz[i] * pa_x[i];

        tr_z_xxxxxz_yyyzzz[i] = 4.0 * tr_z_xxxz_yyyzzz[i] * fe_0 + tr_z_xxxxz_yyyzzz[i] * pa_x[i];

        tr_z_xxxxxz_yyzzzz[i] = 4.0 * tr_z_xxxz_yyzzzz[i] * fe_0 + tr_z_xxxxz_yyzzzz[i] * pa_x[i];

        tr_z_xxxxxz_yzzzzz[i] = 4.0 * tr_z_xxxz_yzzzzz[i] * fe_0 + tr_z_xxxxz_yzzzzz[i] * pa_x[i];

        tr_z_xxxxxz_zzzzzz[i] = 4.0 * tr_z_xxxz_zzzzzz[i] * fe_0 + tr_z_xxxxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1652-1680 components of targeted buffer : II

    auto tr_z_xxxxyy_xxxxxx = pbuffer.data(idx_dip_ii + 1652);

    auto tr_z_xxxxyy_xxxxxy = pbuffer.data(idx_dip_ii + 1653);

    auto tr_z_xxxxyy_xxxxxz = pbuffer.data(idx_dip_ii + 1654);

    auto tr_z_xxxxyy_xxxxyy = pbuffer.data(idx_dip_ii + 1655);

    auto tr_z_xxxxyy_xxxxyz = pbuffer.data(idx_dip_ii + 1656);

    auto tr_z_xxxxyy_xxxxzz = pbuffer.data(idx_dip_ii + 1657);

    auto tr_z_xxxxyy_xxxyyy = pbuffer.data(idx_dip_ii + 1658);

    auto tr_z_xxxxyy_xxxyyz = pbuffer.data(idx_dip_ii + 1659);

    auto tr_z_xxxxyy_xxxyzz = pbuffer.data(idx_dip_ii + 1660);

    auto tr_z_xxxxyy_xxxzzz = pbuffer.data(idx_dip_ii + 1661);

    auto tr_z_xxxxyy_xxyyyy = pbuffer.data(idx_dip_ii + 1662);

    auto tr_z_xxxxyy_xxyyyz = pbuffer.data(idx_dip_ii + 1663);

    auto tr_z_xxxxyy_xxyyzz = pbuffer.data(idx_dip_ii + 1664);

    auto tr_z_xxxxyy_xxyzzz = pbuffer.data(idx_dip_ii + 1665);

    auto tr_z_xxxxyy_xxzzzz = pbuffer.data(idx_dip_ii + 1666);

    auto tr_z_xxxxyy_xyyyyy = pbuffer.data(idx_dip_ii + 1667);

    auto tr_z_xxxxyy_xyyyyz = pbuffer.data(idx_dip_ii + 1668);

    auto tr_z_xxxxyy_xyyyzz = pbuffer.data(idx_dip_ii + 1669);

    auto tr_z_xxxxyy_xyyzzz = pbuffer.data(idx_dip_ii + 1670);

    auto tr_z_xxxxyy_xyzzzz = pbuffer.data(idx_dip_ii + 1671);

    auto tr_z_xxxxyy_xzzzzz = pbuffer.data(idx_dip_ii + 1672);

    auto tr_z_xxxxyy_yyyyyy = pbuffer.data(idx_dip_ii + 1673);

    auto tr_z_xxxxyy_yyyyyz = pbuffer.data(idx_dip_ii + 1674);

    auto tr_z_xxxxyy_yyyyzz = pbuffer.data(idx_dip_ii + 1675);

    auto tr_z_xxxxyy_yyyzzz = pbuffer.data(idx_dip_ii + 1676);

    auto tr_z_xxxxyy_yyzzzz = pbuffer.data(idx_dip_ii + 1677);

    auto tr_z_xxxxyy_yzzzzz = pbuffer.data(idx_dip_ii + 1678);

    auto tr_z_xxxxyy_zzzzzz = pbuffer.data(idx_dip_ii + 1679);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxx_xxxxxx, tr_z_xxxx_xxxxxz, tr_z_xxxx_xxxxzz, tr_z_xxxx_xxxzzz, tr_z_xxxx_xxzzzz, tr_z_xxxx_xzzzzz, tr_z_xxxxy_xxxxxx, tr_z_xxxxy_xxxxxz, tr_z_xxxxy_xxxxzz, tr_z_xxxxy_xxxzzz, tr_z_xxxxy_xxzzzz, tr_z_xxxxy_xzzzzz, tr_z_xxxxyy_xxxxxx, tr_z_xxxxyy_xxxxxy, tr_z_xxxxyy_xxxxxz, tr_z_xxxxyy_xxxxyy, tr_z_xxxxyy_xxxxyz, tr_z_xxxxyy_xxxxzz, tr_z_xxxxyy_xxxyyy, tr_z_xxxxyy_xxxyyz, tr_z_xxxxyy_xxxyzz, tr_z_xxxxyy_xxxzzz, tr_z_xxxxyy_xxyyyy, tr_z_xxxxyy_xxyyyz, tr_z_xxxxyy_xxyyzz, tr_z_xxxxyy_xxyzzz, tr_z_xxxxyy_xxzzzz, tr_z_xxxxyy_xyyyyy, tr_z_xxxxyy_xyyyyz, tr_z_xxxxyy_xyyyzz, tr_z_xxxxyy_xyyzzz, tr_z_xxxxyy_xyzzzz, tr_z_xxxxyy_xzzzzz, tr_z_xxxxyy_yyyyyy, tr_z_xxxxyy_yyyyyz, tr_z_xxxxyy_yyyyzz, tr_z_xxxxyy_yyyzzz, tr_z_xxxxyy_yyzzzz, tr_z_xxxxyy_yzzzzz, tr_z_xxxxyy_zzzzzz, tr_z_xxxyy_xxxxxy, tr_z_xxxyy_xxxxy, tr_z_xxxyy_xxxxyy, tr_z_xxxyy_xxxxyz, tr_z_xxxyy_xxxyy, tr_z_xxxyy_xxxyyy, tr_z_xxxyy_xxxyyz, tr_z_xxxyy_xxxyz, tr_z_xxxyy_xxxyzz, tr_z_xxxyy_xxyyy, tr_z_xxxyy_xxyyyy, tr_z_xxxyy_xxyyyz, tr_z_xxxyy_xxyyz, tr_z_xxxyy_xxyyzz, tr_z_xxxyy_xxyzz, tr_z_xxxyy_xxyzzz, tr_z_xxxyy_xyyyy, tr_z_xxxyy_xyyyyy, tr_z_xxxyy_xyyyyz, tr_z_xxxyy_xyyyz, tr_z_xxxyy_xyyyzz, tr_z_xxxyy_xyyzz, tr_z_xxxyy_xyyzzz, tr_z_xxxyy_xyzzz, tr_z_xxxyy_xyzzzz, tr_z_xxxyy_yyyyy, tr_z_xxxyy_yyyyyy, tr_z_xxxyy_yyyyyz, tr_z_xxxyy_yyyyz, tr_z_xxxyy_yyyyzz, tr_z_xxxyy_yyyzz, tr_z_xxxyy_yyyzzz, tr_z_xxxyy_yyzzz, tr_z_xxxyy_yyzzzz, tr_z_xxxyy_yzzzz, tr_z_xxxyy_yzzzzz, tr_z_xxxyy_zzzzzz, tr_z_xxyy_xxxxxy, tr_z_xxyy_xxxxyy, tr_z_xxyy_xxxxyz, tr_z_xxyy_xxxyyy, tr_z_xxyy_xxxyyz, tr_z_xxyy_xxxyzz, tr_z_xxyy_xxyyyy, tr_z_xxyy_xxyyyz, tr_z_xxyy_xxyyzz, tr_z_xxyy_xxyzzz, tr_z_xxyy_xyyyyy, tr_z_xxyy_xyyyyz, tr_z_xxyy_xyyyzz, tr_z_xxyy_xyyzzz, tr_z_xxyy_xyzzzz, tr_z_xxyy_yyyyyy, tr_z_xxyy_yyyyyz, tr_z_xxyy_yyyyzz, tr_z_xxyy_yyyzzz, tr_z_xxyy_yyzzzz, tr_z_xxyy_yzzzzz, tr_z_xxyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyy_xxxxxx[i] = tr_z_xxxx_xxxxxx[i] * fe_0 + tr_z_xxxxy_xxxxxx[i] * pa_y[i];

        tr_z_xxxxyy_xxxxxy[i] = 3.0 * tr_z_xxyy_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxxyy_xxxxy[i] * fe_0 + tr_z_xxxyy_xxxxxy[i] * pa_x[i];

        tr_z_xxxxyy_xxxxxz[i] = tr_z_xxxx_xxxxxz[i] * fe_0 + tr_z_xxxxy_xxxxxz[i] * pa_y[i];

        tr_z_xxxxyy_xxxxyy[i] = 3.0 * tr_z_xxyy_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxxyy_xxxyy[i] * fe_0 + tr_z_xxxyy_xxxxyy[i] * pa_x[i];

        tr_z_xxxxyy_xxxxyz[i] = 3.0 * tr_z_xxyy_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxxyy_xxxyz[i] * fe_0 + tr_z_xxxyy_xxxxyz[i] * pa_x[i];

        tr_z_xxxxyy_xxxxzz[i] = tr_z_xxxx_xxxxzz[i] * fe_0 + tr_z_xxxxy_xxxxzz[i] * pa_y[i];

        tr_z_xxxxyy_xxxyyy[i] = 3.0 * tr_z_xxyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxxyy_xxyyy[i] * fe_0 + tr_z_xxxyy_xxxyyy[i] * pa_x[i];

        tr_z_xxxxyy_xxxyyz[i] = 3.0 * tr_z_xxyy_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxxyy_xxyyz[i] * fe_0 + tr_z_xxxyy_xxxyyz[i] * pa_x[i];

        tr_z_xxxxyy_xxxyzz[i] = 3.0 * tr_z_xxyy_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxxyy_xxyzz[i] * fe_0 + tr_z_xxxyy_xxxyzz[i] * pa_x[i];

        tr_z_xxxxyy_xxxzzz[i] = tr_z_xxxx_xxxzzz[i] * fe_0 + tr_z_xxxxy_xxxzzz[i] * pa_y[i];

        tr_z_xxxxyy_xxyyyy[i] = 3.0 * tr_z_xxyy_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxxyy_xyyyy[i] * fe_0 + tr_z_xxxyy_xxyyyy[i] * pa_x[i];

        tr_z_xxxxyy_xxyyyz[i] = 3.0 * tr_z_xxyy_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxxyy_xyyyz[i] * fe_0 + tr_z_xxxyy_xxyyyz[i] * pa_x[i];

        tr_z_xxxxyy_xxyyzz[i] = 3.0 * tr_z_xxyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxxyy_xyyzz[i] * fe_0 + tr_z_xxxyy_xxyyzz[i] * pa_x[i];

        tr_z_xxxxyy_xxyzzz[i] = 3.0 * tr_z_xxyy_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxxyy_xyzzz[i] * fe_0 + tr_z_xxxyy_xxyzzz[i] * pa_x[i];

        tr_z_xxxxyy_xxzzzz[i] = tr_z_xxxx_xxzzzz[i] * fe_0 + tr_z_xxxxy_xxzzzz[i] * pa_y[i];

        tr_z_xxxxyy_xyyyyy[i] = 3.0 * tr_z_xxyy_xyyyyy[i] * fe_0 + tr_z_xxxyy_yyyyy[i] * fe_0 + tr_z_xxxyy_xyyyyy[i] * pa_x[i];

        tr_z_xxxxyy_xyyyyz[i] = 3.0 * tr_z_xxyy_xyyyyz[i] * fe_0 + tr_z_xxxyy_yyyyz[i] * fe_0 + tr_z_xxxyy_xyyyyz[i] * pa_x[i];

        tr_z_xxxxyy_xyyyzz[i] = 3.0 * tr_z_xxyy_xyyyzz[i] * fe_0 + tr_z_xxxyy_yyyzz[i] * fe_0 + tr_z_xxxyy_xyyyzz[i] * pa_x[i];

        tr_z_xxxxyy_xyyzzz[i] = 3.0 * tr_z_xxyy_xyyzzz[i] * fe_0 + tr_z_xxxyy_yyzzz[i] * fe_0 + tr_z_xxxyy_xyyzzz[i] * pa_x[i];

        tr_z_xxxxyy_xyzzzz[i] = 3.0 * tr_z_xxyy_xyzzzz[i] * fe_0 + tr_z_xxxyy_yzzzz[i] * fe_0 + tr_z_xxxyy_xyzzzz[i] * pa_x[i];

        tr_z_xxxxyy_xzzzzz[i] = tr_z_xxxx_xzzzzz[i] * fe_0 + tr_z_xxxxy_xzzzzz[i] * pa_y[i];

        tr_z_xxxxyy_yyyyyy[i] = 3.0 * tr_z_xxyy_yyyyyy[i] * fe_0 + tr_z_xxxyy_yyyyyy[i] * pa_x[i];

        tr_z_xxxxyy_yyyyyz[i] = 3.0 * tr_z_xxyy_yyyyyz[i] * fe_0 + tr_z_xxxyy_yyyyyz[i] * pa_x[i];

        tr_z_xxxxyy_yyyyzz[i] = 3.0 * tr_z_xxyy_yyyyzz[i] * fe_0 + tr_z_xxxyy_yyyyzz[i] * pa_x[i];

        tr_z_xxxxyy_yyyzzz[i] = 3.0 * tr_z_xxyy_yyyzzz[i] * fe_0 + tr_z_xxxyy_yyyzzz[i] * pa_x[i];

        tr_z_xxxxyy_yyzzzz[i] = 3.0 * tr_z_xxyy_yyzzzz[i] * fe_0 + tr_z_xxxyy_yyzzzz[i] * pa_x[i];

        tr_z_xxxxyy_yzzzzz[i] = 3.0 * tr_z_xxyy_yzzzzz[i] * fe_0 + tr_z_xxxyy_yzzzzz[i] * pa_x[i];

        tr_z_xxxxyy_zzzzzz[i] = 3.0 * tr_z_xxyy_zzzzzz[i] * fe_0 + tr_z_xxxyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1680-1708 components of targeted buffer : II

    auto tr_z_xxxxyz_xxxxxx = pbuffer.data(idx_dip_ii + 1680);

    auto tr_z_xxxxyz_xxxxxy = pbuffer.data(idx_dip_ii + 1681);

    auto tr_z_xxxxyz_xxxxxz = pbuffer.data(idx_dip_ii + 1682);

    auto tr_z_xxxxyz_xxxxyy = pbuffer.data(idx_dip_ii + 1683);

    auto tr_z_xxxxyz_xxxxyz = pbuffer.data(idx_dip_ii + 1684);

    auto tr_z_xxxxyz_xxxxzz = pbuffer.data(idx_dip_ii + 1685);

    auto tr_z_xxxxyz_xxxyyy = pbuffer.data(idx_dip_ii + 1686);

    auto tr_z_xxxxyz_xxxyyz = pbuffer.data(idx_dip_ii + 1687);

    auto tr_z_xxxxyz_xxxyzz = pbuffer.data(idx_dip_ii + 1688);

    auto tr_z_xxxxyz_xxxzzz = pbuffer.data(idx_dip_ii + 1689);

    auto tr_z_xxxxyz_xxyyyy = pbuffer.data(idx_dip_ii + 1690);

    auto tr_z_xxxxyz_xxyyyz = pbuffer.data(idx_dip_ii + 1691);

    auto tr_z_xxxxyz_xxyyzz = pbuffer.data(idx_dip_ii + 1692);

    auto tr_z_xxxxyz_xxyzzz = pbuffer.data(idx_dip_ii + 1693);

    auto tr_z_xxxxyz_xxzzzz = pbuffer.data(idx_dip_ii + 1694);

    auto tr_z_xxxxyz_xyyyyy = pbuffer.data(idx_dip_ii + 1695);

    auto tr_z_xxxxyz_xyyyyz = pbuffer.data(idx_dip_ii + 1696);

    auto tr_z_xxxxyz_xyyyzz = pbuffer.data(idx_dip_ii + 1697);

    auto tr_z_xxxxyz_xyyzzz = pbuffer.data(idx_dip_ii + 1698);

    auto tr_z_xxxxyz_xyzzzz = pbuffer.data(idx_dip_ii + 1699);

    auto tr_z_xxxxyz_xzzzzz = pbuffer.data(idx_dip_ii + 1700);

    auto tr_z_xxxxyz_yyyyyy = pbuffer.data(idx_dip_ii + 1701);

    auto tr_z_xxxxyz_yyyyyz = pbuffer.data(idx_dip_ii + 1702);

    auto tr_z_xxxxyz_yyyyzz = pbuffer.data(idx_dip_ii + 1703);

    auto tr_z_xxxxyz_yyyzzz = pbuffer.data(idx_dip_ii + 1704);

    auto tr_z_xxxxyz_yyzzzz = pbuffer.data(idx_dip_ii + 1705);

    auto tr_z_xxxxyz_yzzzzz = pbuffer.data(idx_dip_ii + 1706);

    auto tr_z_xxxxyz_zzzzzz = pbuffer.data(idx_dip_ii + 1707);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxyz_xxxxxx, tr_z_xxxxyz_xxxxxy, tr_z_xxxxyz_xxxxxz, tr_z_xxxxyz_xxxxyy, tr_z_xxxxyz_xxxxyz, tr_z_xxxxyz_xxxxzz, tr_z_xxxxyz_xxxyyy, tr_z_xxxxyz_xxxyyz, tr_z_xxxxyz_xxxyzz, tr_z_xxxxyz_xxxzzz, tr_z_xxxxyz_xxyyyy, tr_z_xxxxyz_xxyyyz, tr_z_xxxxyz_xxyyzz, tr_z_xxxxyz_xxyzzz, tr_z_xxxxyz_xxzzzz, tr_z_xxxxyz_xyyyyy, tr_z_xxxxyz_xyyyyz, tr_z_xxxxyz_xyyyzz, tr_z_xxxxyz_xyyzzz, tr_z_xxxxyz_xyzzzz, tr_z_xxxxyz_xzzzzz, tr_z_xxxxyz_yyyyyy, tr_z_xxxxyz_yyyyyz, tr_z_xxxxyz_yyyyzz, tr_z_xxxxyz_yyyzzz, tr_z_xxxxyz_yyzzzz, tr_z_xxxxyz_yzzzzz, tr_z_xxxxyz_zzzzzz, tr_z_xxxxz_xxxxx, tr_z_xxxxz_xxxxxx, tr_z_xxxxz_xxxxxy, tr_z_xxxxz_xxxxxz, tr_z_xxxxz_xxxxy, tr_z_xxxxz_xxxxyy, tr_z_xxxxz_xxxxyz, tr_z_xxxxz_xxxxz, tr_z_xxxxz_xxxxzz, tr_z_xxxxz_xxxyy, tr_z_xxxxz_xxxyyy, tr_z_xxxxz_xxxyyz, tr_z_xxxxz_xxxyz, tr_z_xxxxz_xxxyzz, tr_z_xxxxz_xxxzz, tr_z_xxxxz_xxxzzz, tr_z_xxxxz_xxyyy, tr_z_xxxxz_xxyyyy, tr_z_xxxxz_xxyyyz, tr_z_xxxxz_xxyyz, tr_z_xxxxz_xxyyzz, tr_z_xxxxz_xxyzz, tr_z_xxxxz_xxyzzz, tr_z_xxxxz_xxzzz, tr_z_xxxxz_xxzzzz, tr_z_xxxxz_xyyyy, tr_z_xxxxz_xyyyyy, tr_z_xxxxz_xyyyyz, tr_z_xxxxz_xyyyz, tr_z_xxxxz_xyyyzz, tr_z_xxxxz_xyyzz, tr_z_xxxxz_xyyzzz, tr_z_xxxxz_xyzzz, tr_z_xxxxz_xyzzzz, tr_z_xxxxz_xzzzz, tr_z_xxxxz_xzzzzz, tr_z_xxxxz_zzzzzz, tr_z_xxxyz_yyyyyy, tr_z_xxxyz_yyyyyz, tr_z_xxxyz_yyyyzz, tr_z_xxxyz_yyyzzz, tr_z_xxxyz_yyzzzz, tr_z_xxxyz_yzzzzz, tr_z_xxyz_yyyyyy, tr_z_xxyz_yyyyyz, tr_z_xxyz_yyyyzz, tr_z_xxyz_yyyzzz, tr_z_xxyz_yyzzzz, tr_z_xxyz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyz_xxxxxx[i] = tr_z_xxxxz_xxxxxx[i] * pa_y[i];

        tr_z_xxxxyz_xxxxxy[i] = tr_z_xxxxz_xxxxx[i] * fe_0 + tr_z_xxxxz_xxxxxy[i] * pa_y[i];

        tr_z_xxxxyz_xxxxxz[i] = tr_z_xxxxz_xxxxxz[i] * pa_y[i];

        tr_z_xxxxyz_xxxxyy[i] = 2.0 * tr_z_xxxxz_xxxxy[i] * fe_0 + tr_z_xxxxz_xxxxyy[i] * pa_y[i];

        tr_z_xxxxyz_xxxxyz[i] = tr_z_xxxxz_xxxxz[i] * fe_0 + tr_z_xxxxz_xxxxyz[i] * pa_y[i];

        tr_z_xxxxyz_xxxxzz[i] = tr_z_xxxxz_xxxxzz[i] * pa_y[i];

        tr_z_xxxxyz_xxxyyy[i] = 3.0 * tr_z_xxxxz_xxxyy[i] * fe_0 + tr_z_xxxxz_xxxyyy[i] * pa_y[i];

        tr_z_xxxxyz_xxxyyz[i] = 2.0 * tr_z_xxxxz_xxxyz[i] * fe_0 + tr_z_xxxxz_xxxyyz[i] * pa_y[i];

        tr_z_xxxxyz_xxxyzz[i] = tr_z_xxxxz_xxxzz[i] * fe_0 + tr_z_xxxxz_xxxyzz[i] * pa_y[i];

        tr_z_xxxxyz_xxxzzz[i] = tr_z_xxxxz_xxxzzz[i] * pa_y[i];

        tr_z_xxxxyz_xxyyyy[i] = 4.0 * tr_z_xxxxz_xxyyy[i] * fe_0 + tr_z_xxxxz_xxyyyy[i] * pa_y[i];

        tr_z_xxxxyz_xxyyyz[i] = 3.0 * tr_z_xxxxz_xxyyz[i] * fe_0 + tr_z_xxxxz_xxyyyz[i] * pa_y[i];

        tr_z_xxxxyz_xxyyzz[i] = 2.0 * tr_z_xxxxz_xxyzz[i] * fe_0 + tr_z_xxxxz_xxyyzz[i] * pa_y[i];

        tr_z_xxxxyz_xxyzzz[i] = tr_z_xxxxz_xxzzz[i] * fe_0 + tr_z_xxxxz_xxyzzz[i] * pa_y[i];

        tr_z_xxxxyz_xxzzzz[i] = tr_z_xxxxz_xxzzzz[i] * pa_y[i];

        tr_z_xxxxyz_xyyyyy[i] = 5.0 * tr_z_xxxxz_xyyyy[i] * fe_0 + tr_z_xxxxz_xyyyyy[i] * pa_y[i];

        tr_z_xxxxyz_xyyyyz[i] = 4.0 * tr_z_xxxxz_xyyyz[i] * fe_0 + tr_z_xxxxz_xyyyyz[i] * pa_y[i];

        tr_z_xxxxyz_xyyyzz[i] = 3.0 * tr_z_xxxxz_xyyzz[i] * fe_0 + tr_z_xxxxz_xyyyzz[i] * pa_y[i];

        tr_z_xxxxyz_xyyzzz[i] = 2.0 * tr_z_xxxxz_xyzzz[i] * fe_0 + tr_z_xxxxz_xyyzzz[i] * pa_y[i];

        tr_z_xxxxyz_xyzzzz[i] = tr_z_xxxxz_xzzzz[i] * fe_0 + tr_z_xxxxz_xyzzzz[i] * pa_y[i];

        tr_z_xxxxyz_xzzzzz[i] = tr_z_xxxxz_xzzzzz[i] * pa_y[i];

        tr_z_xxxxyz_yyyyyy[i] = 3.0 * tr_z_xxyz_yyyyyy[i] * fe_0 + tr_z_xxxyz_yyyyyy[i] * pa_x[i];

        tr_z_xxxxyz_yyyyyz[i] = 3.0 * tr_z_xxyz_yyyyyz[i] * fe_0 + tr_z_xxxyz_yyyyyz[i] * pa_x[i];

        tr_z_xxxxyz_yyyyzz[i] = 3.0 * tr_z_xxyz_yyyyzz[i] * fe_0 + tr_z_xxxyz_yyyyzz[i] * pa_x[i];

        tr_z_xxxxyz_yyyzzz[i] = 3.0 * tr_z_xxyz_yyyzzz[i] * fe_0 + tr_z_xxxyz_yyyzzz[i] * pa_x[i];

        tr_z_xxxxyz_yyzzzz[i] = 3.0 * tr_z_xxyz_yyzzzz[i] * fe_0 + tr_z_xxxyz_yyzzzz[i] * pa_x[i];

        tr_z_xxxxyz_yzzzzz[i] = 3.0 * tr_z_xxyz_yzzzzz[i] * fe_0 + tr_z_xxxyz_yzzzzz[i] * pa_x[i];

        tr_z_xxxxyz_zzzzzz[i] = tr_z_xxxxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1708-1736 components of targeted buffer : II

    auto tr_z_xxxxzz_xxxxxx = pbuffer.data(idx_dip_ii + 1708);

    auto tr_z_xxxxzz_xxxxxy = pbuffer.data(idx_dip_ii + 1709);

    auto tr_z_xxxxzz_xxxxxz = pbuffer.data(idx_dip_ii + 1710);

    auto tr_z_xxxxzz_xxxxyy = pbuffer.data(idx_dip_ii + 1711);

    auto tr_z_xxxxzz_xxxxyz = pbuffer.data(idx_dip_ii + 1712);

    auto tr_z_xxxxzz_xxxxzz = pbuffer.data(idx_dip_ii + 1713);

    auto tr_z_xxxxzz_xxxyyy = pbuffer.data(idx_dip_ii + 1714);

    auto tr_z_xxxxzz_xxxyyz = pbuffer.data(idx_dip_ii + 1715);

    auto tr_z_xxxxzz_xxxyzz = pbuffer.data(idx_dip_ii + 1716);

    auto tr_z_xxxxzz_xxxzzz = pbuffer.data(idx_dip_ii + 1717);

    auto tr_z_xxxxzz_xxyyyy = pbuffer.data(idx_dip_ii + 1718);

    auto tr_z_xxxxzz_xxyyyz = pbuffer.data(idx_dip_ii + 1719);

    auto tr_z_xxxxzz_xxyyzz = pbuffer.data(idx_dip_ii + 1720);

    auto tr_z_xxxxzz_xxyzzz = pbuffer.data(idx_dip_ii + 1721);

    auto tr_z_xxxxzz_xxzzzz = pbuffer.data(idx_dip_ii + 1722);

    auto tr_z_xxxxzz_xyyyyy = pbuffer.data(idx_dip_ii + 1723);

    auto tr_z_xxxxzz_xyyyyz = pbuffer.data(idx_dip_ii + 1724);

    auto tr_z_xxxxzz_xyyyzz = pbuffer.data(idx_dip_ii + 1725);

    auto tr_z_xxxxzz_xyyzzz = pbuffer.data(idx_dip_ii + 1726);

    auto tr_z_xxxxzz_xyzzzz = pbuffer.data(idx_dip_ii + 1727);

    auto tr_z_xxxxzz_xzzzzz = pbuffer.data(idx_dip_ii + 1728);

    auto tr_z_xxxxzz_yyyyyy = pbuffer.data(idx_dip_ii + 1729);

    auto tr_z_xxxxzz_yyyyyz = pbuffer.data(idx_dip_ii + 1730);

    auto tr_z_xxxxzz_yyyyzz = pbuffer.data(idx_dip_ii + 1731);

    auto tr_z_xxxxzz_yyyzzz = pbuffer.data(idx_dip_ii + 1732);

    auto tr_z_xxxxzz_yyzzzz = pbuffer.data(idx_dip_ii + 1733);

    auto tr_z_xxxxzz_yzzzzz = pbuffer.data(idx_dip_ii + 1734);

    auto tr_z_xxxxzz_zzzzzz = pbuffer.data(idx_dip_ii + 1735);

    #pragma omp simd aligned(pa_x, tr_z_xxxxzz_xxxxxx, tr_z_xxxxzz_xxxxxy, tr_z_xxxxzz_xxxxxz, tr_z_xxxxzz_xxxxyy, tr_z_xxxxzz_xxxxyz, tr_z_xxxxzz_xxxxzz, tr_z_xxxxzz_xxxyyy, tr_z_xxxxzz_xxxyyz, tr_z_xxxxzz_xxxyzz, tr_z_xxxxzz_xxxzzz, tr_z_xxxxzz_xxyyyy, tr_z_xxxxzz_xxyyyz, tr_z_xxxxzz_xxyyzz, tr_z_xxxxzz_xxyzzz, tr_z_xxxxzz_xxzzzz, tr_z_xxxxzz_xyyyyy, tr_z_xxxxzz_xyyyyz, tr_z_xxxxzz_xyyyzz, tr_z_xxxxzz_xyyzzz, tr_z_xxxxzz_xyzzzz, tr_z_xxxxzz_xzzzzz, tr_z_xxxxzz_yyyyyy, tr_z_xxxxzz_yyyyyz, tr_z_xxxxzz_yyyyzz, tr_z_xxxxzz_yyyzzz, tr_z_xxxxzz_yyzzzz, tr_z_xxxxzz_yzzzzz, tr_z_xxxxzz_zzzzzz, tr_z_xxxzz_xxxxx, tr_z_xxxzz_xxxxxx, tr_z_xxxzz_xxxxxy, tr_z_xxxzz_xxxxxz, tr_z_xxxzz_xxxxy, tr_z_xxxzz_xxxxyy, tr_z_xxxzz_xxxxyz, tr_z_xxxzz_xxxxz, tr_z_xxxzz_xxxxzz, tr_z_xxxzz_xxxyy, tr_z_xxxzz_xxxyyy, tr_z_xxxzz_xxxyyz, tr_z_xxxzz_xxxyz, tr_z_xxxzz_xxxyzz, tr_z_xxxzz_xxxzz, tr_z_xxxzz_xxxzzz, tr_z_xxxzz_xxyyy, tr_z_xxxzz_xxyyyy, tr_z_xxxzz_xxyyyz, tr_z_xxxzz_xxyyz, tr_z_xxxzz_xxyyzz, tr_z_xxxzz_xxyzz, tr_z_xxxzz_xxyzzz, tr_z_xxxzz_xxzzz, tr_z_xxxzz_xxzzzz, tr_z_xxxzz_xyyyy, tr_z_xxxzz_xyyyyy, tr_z_xxxzz_xyyyyz, tr_z_xxxzz_xyyyz, tr_z_xxxzz_xyyyzz, tr_z_xxxzz_xyyzz, tr_z_xxxzz_xyyzzz, tr_z_xxxzz_xyzzz, tr_z_xxxzz_xyzzzz, tr_z_xxxzz_xzzzz, tr_z_xxxzz_xzzzzz, tr_z_xxxzz_yyyyy, tr_z_xxxzz_yyyyyy, tr_z_xxxzz_yyyyyz, tr_z_xxxzz_yyyyz, tr_z_xxxzz_yyyyzz, tr_z_xxxzz_yyyzz, tr_z_xxxzz_yyyzzz, tr_z_xxxzz_yyzzz, tr_z_xxxzz_yyzzzz, tr_z_xxxzz_yzzzz, tr_z_xxxzz_yzzzzz, tr_z_xxxzz_zzzzz, tr_z_xxxzz_zzzzzz, tr_z_xxzz_xxxxxx, tr_z_xxzz_xxxxxy, tr_z_xxzz_xxxxxz, tr_z_xxzz_xxxxyy, tr_z_xxzz_xxxxyz, tr_z_xxzz_xxxxzz, tr_z_xxzz_xxxyyy, tr_z_xxzz_xxxyyz, tr_z_xxzz_xxxyzz, tr_z_xxzz_xxxzzz, tr_z_xxzz_xxyyyy, tr_z_xxzz_xxyyyz, tr_z_xxzz_xxyyzz, tr_z_xxzz_xxyzzz, tr_z_xxzz_xxzzzz, tr_z_xxzz_xyyyyy, tr_z_xxzz_xyyyyz, tr_z_xxzz_xyyyzz, tr_z_xxzz_xyyzzz, tr_z_xxzz_xyzzzz, tr_z_xxzz_xzzzzz, tr_z_xxzz_yyyyyy, tr_z_xxzz_yyyyyz, tr_z_xxzz_yyyyzz, tr_z_xxzz_yyyzzz, tr_z_xxzz_yyzzzz, tr_z_xxzz_yzzzzz, tr_z_xxzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxzz_xxxxxx[i] = 3.0 * tr_z_xxzz_xxxxxx[i] * fe_0 + 6.0 * tr_z_xxxzz_xxxxx[i] * fe_0 + tr_z_xxxzz_xxxxxx[i] * pa_x[i];

        tr_z_xxxxzz_xxxxxy[i] = 3.0 * tr_z_xxzz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxxzz_xxxxy[i] * fe_0 + tr_z_xxxzz_xxxxxy[i] * pa_x[i];

        tr_z_xxxxzz_xxxxxz[i] = 3.0 * tr_z_xxzz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxxzz_xxxxz[i] * fe_0 + tr_z_xxxzz_xxxxxz[i] * pa_x[i];

        tr_z_xxxxzz_xxxxyy[i] = 3.0 * tr_z_xxzz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxxzz_xxxyy[i] * fe_0 + tr_z_xxxzz_xxxxyy[i] * pa_x[i];

        tr_z_xxxxzz_xxxxyz[i] = 3.0 * tr_z_xxzz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxxzz_xxxyz[i] * fe_0 + tr_z_xxxzz_xxxxyz[i] * pa_x[i];

        tr_z_xxxxzz_xxxxzz[i] = 3.0 * tr_z_xxzz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxxzz_xxxzz[i] * fe_0 + tr_z_xxxzz_xxxxzz[i] * pa_x[i];

        tr_z_xxxxzz_xxxyyy[i] = 3.0 * tr_z_xxzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxxzz_xxyyy[i] * fe_0 + tr_z_xxxzz_xxxyyy[i] * pa_x[i];

        tr_z_xxxxzz_xxxyyz[i] = 3.0 * tr_z_xxzz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxxzz_xxyyz[i] * fe_0 + tr_z_xxxzz_xxxyyz[i] * pa_x[i];

        tr_z_xxxxzz_xxxyzz[i] = 3.0 * tr_z_xxzz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxxzz_xxyzz[i] * fe_0 + tr_z_xxxzz_xxxyzz[i] * pa_x[i];

        tr_z_xxxxzz_xxxzzz[i] = 3.0 * tr_z_xxzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxxzz_xxzzz[i] * fe_0 + tr_z_xxxzz_xxxzzz[i] * pa_x[i];

        tr_z_xxxxzz_xxyyyy[i] = 3.0 * tr_z_xxzz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxxzz_xyyyy[i] * fe_0 + tr_z_xxxzz_xxyyyy[i] * pa_x[i];

        tr_z_xxxxzz_xxyyyz[i] = 3.0 * tr_z_xxzz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxxzz_xyyyz[i] * fe_0 + tr_z_xxxzz_xxyyyz[i] * pa_x[i];

        tr_z_xxxxzz_xxyyzz[i] = 3.0 * tr_z_xxzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxxzz_xyyzz[i] * fe_0 + tr_z_xxxzz_xxyyzz[i] * pa_x[i];

        tr_z_xxxxzz_xxyzzz[i] = 3.0 * tr_z_xxzz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxxzz_xyzzz[i] * fe_0 + tr_z_xxxzz_xxyzzz[i] * pa_x[i];

        tr_z_xxxxzz_xxzzzz[i] = 3.0 * tr_z_xxzz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxxzz_xzzzz[i] * fe_0 + tr_z_xxxzz_xxzzzz[i] * pa_x[i];

        tr_z_xxxxzz_xyyyyy[i] = 3.0 * tr_z_xxzz_xyyyyy[i] * fe_0 + tr_z_xxxzz_yyyyy[i] * fe_0 + tr_z_xxxzz_xyyyyy[i] * pa_x[i];

        tr_z_xxxxzz_xyyyyz[i] = 3.0 * tr_z_xxzz_xyyyyz[i] * fe_0 + tr_z_xxxzz_yyyyz[i] * fe_0 + tr_z_xxxzz_xyyyyz[i] * pa_x[i];

        tr_z_xxxxzz_xyyyzz[i] = 3.0 * tr_z_xxzz_xyyyzz[i] * fe_0 + tr_z_xxxzz_yyyzz[i] * fe_0 + tr_z_xxxzz_xyyyzz[i] * pa_x[i];

        tr_z_xxxxzz_xyyzzz[i] = 3.0 * tr_z_xxzz_xyyzzz[i] * fe_0 + tr_z_xxxzz_yyzzz[i] * fe_0 + tr_z_xxxzz_xyyzzz[i] * pa_x[i];

        tr_z_xxxxzz_xyzzzz[i] = 3.0 * tr_z_xxzz_xyzzzz[i] * fe_0 + tr_z_xxxzz_yzzzz[i] * fe_0 + tr_z_xxxzz_xyzzzz[i] * pa_x[i];

        tr_z_xxxxzz_xzzzzz[i] = 3.0 * tr_z_xxzz_xzzzzz[i] * fe_0 + tr_z_xxxzz_zzzzz[i] * fe_0 + tr_z_xxxzz_xzzzzz[i] * pa_x[i];

        tr_z_xxxxzz_yyyyyy[i] = 3.0 * tr_z_xxzz_yyyyyy[i] * fe_0 + tr_z_xxxzz_yyyyyy[i] * pa_x[i];

        tr_z_xxxxzz_yyyyyz[i] = 3.0 * tr_z_xxzz_yyyyyz[i] * fe_0 + tr_z_xxxzz_yyyyyz[i] * pa_x[i];

        tr_z_xxxxzz_yyyyzz[i] = 3.0 * tr_z_xxzz_yyyyzz[i] * fe_0 + tr_z_xxxzz_yyyyzz[i] * pa_x[i];

        tr_z_xxxxzz_yyyzzz[i] = 3.0 * tr_z_xxzz_yyyzzz[i] * fe_0 + tr_z_xxxzz_yyyzzz[i] * pa_x[i];

        tr_z_xxxxzz_yyzzzz[i] = 3.0 * tr_z_xxzz_yyzzzz[i] * fe_0 + tr_z_xxxzz_yyzzzz[i] * pa_x[i];

        tr_z_xxxxzz_yzzzzz[i] = 3.0 * tr_z_xxzz_yzzzzz[i] * fe_0 + tr_z_xxxzz_yzzzzz[i] * pa_x[i];

        tr_z_xxxxzz_zzzzzz[i] = 3.0 * tr_z_xxzz_zzzzzz[i] * fe_0 + tr_z_xxxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1736-1764 components of targeted buffer : II

    auto tr_z_xxxyyy_xxxxxx = pbuffer.data(idx_dip_ii + 1736);

    auto tr_z_xxxyyy_xxxxxy = pbuffer.data(idx_dip_ii + 1737);

    auto tr_z_xxxyyy_xxxxxz = pbuffer.data(idx_dip_ii + 1738);

    auto tr_z_xxxyyy_xxxxyy = pbuffer.data(idx_dip_ii + 1739);

    auto tr_z_xxxyyy_xxxxyz = pbuffer.data(idx_dip_ii + 1740);

    auto tr_z_xxxyyy_xxxxzz = pbuffer.data(idx_dip_ii + 1741);

    auto tr_z_xxxyyy_xxxyyy = pbuffer.data(idx_dip_ii + 1742);

    auto tr_z_xxxyyy_xxxyyz = pbuffer.data(idx_dip_ii + 1743);

    auto tr_z_xxxyyy_xxxyzz = pbuffer.data(idx_dip_ii + 1744);

    auto tr_z_xxxyyy_xxxzzz = pbuffer.data(idx_dip_ii + 1745);

    auto tr_z_xxxyyy_xxyyyy = pbuffer.data(idx_dip_ii + 1746);

    auto tr_z_xxxyyy_xxyyyz = pbuffer.data(idx_dip_ii + 1747);

    auto tr_z_xxxyyy_xxyyzz = pbuffer.data(idx_dip_ii + 1748);

    auto tr_z_xxxyyy_xxyzzz = pbuffer.data(idx_dip_ii + 1749);

    auto tr_z_xxxyyy_xxzzzz = pbuffer.data(idx_dip_ii + 1750);

    auto tr_z_xxxyyy_xyyyyy = pbuffer.data(idx_dip_ii + 1751);

    auto tr_z_xxxyyy_xyyyyz = pbuffer.data(idx_dip_ii + 1752);

    auto tr_z_xxxyyy_xyyyzz = pbuffer.data(idx_dip_ii + 1753);

    auto tr_z_xxxyyy_xyyzzz = pbuffer.data(idx_dip_ii + 1754);

    auto tr_z_xxxyyy_xyzzzz = pbuffer.data(idx_dip_ii + 1755);

    auto tr_z_xxxyyy_xzzzzz = pbuffer.data(idx_dip_ii + 1756);

    auto tr_z_xxxyyy_yyyyyy = pbuffer.data(idx_dip_ii + 1757);

    auto tr_z_xxxyyy_yyyyyz = pbuffer.data(idx_dip_ii + 1758);

    auto tr_z_xxxyyy_yyyyzz = pbuffer.data(idx_dip_ii + 1759);

    auto tr_z_xxxyyy_yyyzzz = pbuffer.data(idx_dip_ii + 1760);

    auto tr_z_xxxyyy_yyzzzz = pbuffer.data(idx_dip_ii + 1761);

    auto tr_z_xxxyyy_yzzzzz = pbuffer.data(idx_dip_ii + 1762);

    auto tr_z_xxxyyy_zzzzzz = pbuffer.data(idx_dip_ii + 1763);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxy_xxxxxx, tr_z_xxxy_xxxxxz, tr_z_xxxy_xxxxzz, tr_z_xxxy_xxxzzz, tr_z_xxxy_xxzzzz, tr_z_xxxy_xzzzzz, tr_z_xxxyy_xxxxxx, tr_z_xxxyy_xxxxxz, tr_z_xxxyy_xxxxzz, tr_z_xxxyy_xxxzzz, tr_z_xxxyy_xxzzzz, tr_z_xxxyy_xzzzzz, tr_z_xxxyyy_xxxxxx, tr_z_xxxyyy_xxxxxy, tr_z_xxxyyy_xxxxxz, tr_z_xxxyyy_xxxxyy, tr_z_xxxyyy_xxxxyz, tr_z_xxxyyy_xxxxzz, tr_z_xxxyyy_xxxyyy, tr_z_xxxyyy_xxxyyz, tr_z_xxxyyy_xxxyzz, tr_z_xxxyyy_xxxzzz, tr_z_xxxyyy_xxyyyy, tr_z_xxxyyy_xxyyyz, tr_z_xxxyyy_xxyyzz, tr_z_xxxyyy_xxyzzz, tr_z_xxxyyy_xxzzzz, tr_z_xxxyyy_xyyyyy, tr_z_xxxyyy_xyyyyz, tr_z_xxxyyy_xyyyzz, tr_z_xxxyyy_xyyzzz, tr_z_xxxyyy_xyzzzz, tr_z_xxxyyy_xzzzzz, tr_z_xxxyyy_yyyyyy, tr_z_xxxyyy_yyyyyz, tr_z_xxxyyy_yyyyzz, tr_z_xxxyyy_yyyzzz, tr_z_xxxyyy_yyzzzz, tr_z_xxxyyy_yzzzzz, tr_z_xxxyyy_zzzzzz, tr_z_xxyyy_xxxxxy, tr_z_xxyyy_xxxxy, tr_z_xxyyy_xxxxyy, tr_z_xxyyy_xxxxyz, tr_z_xxyyy_xxxyy, tr_z_xxyyy_xxxyyy, tr_z_xxyyy_xxxyyz, tr_z_xxyyy_xxxyz, tr_z_xxyyy_xxxyzz, tr_z_xxyyy_xxyyy, tr_z_xxyyy_xxyyyy, tr_z_xxyyy_xxyyyz, tr_z_xxyyy_xxyyz, tr_z_xxyyy_xxyyzz, tr_z_xxyyy_xxyzz, tr_z_xxyyy_xxyzzz, tr_z_xxyyy_xyyyy, tr_z_xxyyy_xyyyyy, tr_z_xxyyy_xyyyyz, tr_z_xxyyy_xyyyz, tr_z_xxyyy_xyyyzz, tr_z_xxyyy_xyyzz, tr_z_xxyyy_xyyzzz, tr_z_xxyyy_xyzzz, tr_z_xxyyy_xyzzzz, tr_z_xxyyy_yyyyy, tr_z_xxyyy_yyyyyy, tr_z_xxyyy_yyyyyz, tr_z_xxyyy_yyyyz, tr_z_xxyyy_yyyyzz, tr_z_xxyyy_yyyzz, tr_z_xxyyy_yyyzzz, tr_z_xxyyy_yyzzz, tr_z_xxyyy_yyzzzz, tr_z_xxyyy_yzzzz, tr_z_xxyyy_yzzzzz, tr_z_xxyyy_zzzzzz, tr_z_xyyy_xxxxxy, tr_z_xyyy_xxxxyy, tr_z_xyyy_xxxxyz, tr_z_xyyy_xxxyyy, tr_z_xyyy_xxxyyz, tr_z_xyyy_xxxyzz, tr_z_xyyy_xxyyyy, tr_z_xyyy_xxyyyz, tr_z_xyyy_xxyyzz, tr_z_xyyy_xxyzzz, tr_z_xyyy_xyyyyy, tr_z_xyyy_xyyyyz, tr_z_xyyy_xyyyzz, tr_z_xyyy_xyyzzz, tr_z_xyyy_xyzzzz, tr_z_xyyy_yyyyyy, tr_z_xyyy_yyyyyz, tr_z_xyyy_yyyyzz, tr_z_xyyy_yyyzzz, tr_z_xyyy_yyzzzz, tr_z_xyyy_yzzzzz, tr_z_xyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyy_xxxxxx[i] = 2.0 * tr_z_xxxy_xxxxxx[i] * fe_0 + tr_z_xxxyy_xxxxxx[i] * pa_y[i];

        tr_z_xxxyyy_xxxxxy[i] = 2.0 * tr_z_xyyy_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxyyy_xxxxy[i] * fe_0 + tr_z_xxyyy_xxxxxy[i] * pa_x[i];

        tr_z_xxxyyy_xxxxxz[i] = 2.0 * tr_z_xxxy_xxxxxz[i] * fe_0 + tr_z_xxxyy_xxxxxz[i] * pa_y[i];

        tr_z_xxxyyy_xxxxyy[i] = 2.0 * tr_z_xyyy_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxyyy_xxxyy[i] * fe_0 + tr_z_xxyyy_xxxxyy[i] * pa_x[i];

        tr_z_xxxyyy_xxxxyz[i] = 2.0 * tr_z_xyyy_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxyyy_xxxyz[i] * fe_0 + tr_z_xxyyy_xxxxyz[i] * pa_x[i];

        tr_z_xxxyyy_xxxxzz[i] = 2.0 * tr_z_xxxy_xxxxzz[i] * fe_0 + tr_z_xxxyy_xxxxzz[i] * pa_y[i];

        tr_z_xxxyyy_xxxyyy[i] = 2.0 * tr_z_xyyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxyyy_xxyyy[i] * fe_0 + tr_z_xxyyy_xxxyyy[i] * pa_x[i];

        tr_z_xxxyyy_xxxyyz[i] = 2.0 * tr_z_xyyy_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxyyy_xxyyz[i] * fe_0 + tr_z_xxyyy_xxxyyz[i] * pa_x[i];

        tr_z_xxxyyy_xxxyzz[i] = 2.0 * tr_z_xyyy_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxyyy_xxyzz[i] * fe_0 + tr_z_xxyyy_xxxyzz[i] * pa_x[i];

        tr_z_xxxyyy_xxxzzz[i] = 2.0 * tr_z_xxxy_xxxzzz[i] * fe_0 + tr_z_xxxyy_xxxzzz[i] * pa_y[i];

        tr_z_xxxyyy_xxyyyy[i] = 2.0 * tr_z_xyyy_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxyyy_xyyyy[i] * fe_0 + tr_z_xxyyy_xxyyyy[i] * pa_x[i];

        tr_z_xxxyyy_xxyyyz[i] = 2.0 * tr_z_xyyy_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxyyy_xyyyz[i] * fe_0 + tr_z_xxyyy_xxyyyz[i] * pa_x[i];

        tr_z_xxxyyy_xxyyzz[i] = 2.0 * tr_z_xyyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxyyy_xyyzz[i] * fe_0 + tr_z_xxyyy_xxyyzz[i] * pa_x[i];

        tr_z_xxxyyy_xxyzzz[i] = 2.0 * tr_z_xyyy_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxyyy_xyzzz[i] * fe_0 + tr_z_xxyyy_xxyzzz[i] * pa_x[i];

        tr_z_xxxyyy_xxzzzz[i] = 2.0 * tr_z_xxxy_xxzzzz[i] * fe_0 + tr_z_xxxyy_xxzzzz[i] * pa_y[i];

        tr_z_xxxyyy_xyyyyy[i] = 2.0 * tr_z_xyyy_xyyyyy[i] * fe_0 + tr_z_xxyyy_yyyyy[i] * fe_0 + tr_z_xxyyy_xyyyyy[i] * pa_x[i];

        tr_z_xxxyyy_xyyyyz[i] = 2.0 * tr_z_xyyy_xyyyyz[i] * fe_0 + tr_z_xxyyy_yyyyz[i] * fe_0 + tr_z_xxyyy_xyyyyz[i] * pa_x[i];

        tr_z_xxxyyy_xyyyzz[i] = 2.0 * tr_z_xyyy_xyyyzz[i] * fe_0 + tr_z_xxyyy_yyyzz[i] * fe_0 + tr_z_xxyyy_xyyyzz[i] * pa_x[i];

        tr_z_xxxyyy_xyyzzz[i] = 2.0 * tr_z_xyyy_xyyzzz[i] * fe_0 + tr_z_xxyyy_yyzzz[i] * fe_0 + tr_z_xxyyy_xyyzzz[i] * pa_x[i];

        tr_z_xxxyyy_xyzzzz[i] = 2.0 * tr_z_xyyy_xyzzzz[i] * fe_0 + tr_z_xxyyy_yzzzz[i] * fe_0 + tr_z_xxyyy_xyzzzz[i] * pa_x[i];

        tr_z_xxxyyy_xzzzzz[i] = 2.0 * tr_z_xxxy_xzzzzz[i] * fe_0 + tr_z_xxxyy_xzzzzz[i] * pa_y[i];

        tr_z_xxxyyy_yyyyyy[i] = 2.0 * tr_z_xyyy_yyyyyy[i] * fe_0 + tr_z_xxyyy_yyyyyy[i] * pa_x[i];

        tr_z_xxxyyy_yyyyyz[i] = 2.0 * tr_z_xyyy_yyyyyz[i] * fe_0 + tr_z_xxyyy_yyyyyz[i] * pa_x[i];

        tr_z_xxxyyy_yyyyzz[i] = 2.0 * tr_z_xyyy_yyyyzz[i] * fe_0 + tr_z_xxyyy_yyyyzz[i] * pa_x[i];

        tr_z_xxxyyy_yyyzzz[i] = 2.0 * tr_z_xyyy_yyyzzz[i] * fe_0 + tr_z_xxyyy_yyyzzz[i] * pa_x[i];

        tr_z_xxxyyy_yyzzzz[i] = 2.0 * tr_z_xyyy_yyzzzz[i] * fe_0 + tr_z_xxyyy_yyzzzz[i] * pa_x[i];

        tr_z_xxxyyy_yzzzzz[i] = 2.0 * tr_z_xyyy_yzzzzz[i] * fe_0 + tr_z_xxyyy_yzzzzz[i] * pa_x[i];

        tr_z_xxxyyy_zzzzzz[i] = 2.0 * tr_z_xyyy_zzzzzz[i] * fe_0 + tr_z_xxyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1764-1792 components of targeted buffer : II

    auto tr_z_xxxyyz_xxxxxx = pbuffer.data(idx_dip_ii + 1764);

    auto tr_z_xxxyyz_xxxxxy = pbuffer.data(idx_dip_ii + 1765);

    auto tr_z_xxxyyz_xxxxxz = pbuffer.data(idx_dip_ii + 1766);

    auto tr_z_xxxyyz_xxxxyy = pbuffer.data(idx_dip_ii + 1767);

    auto tr_z_xxxyyz_xxxxyz = pbuffer.data(idx_dip_ii + 1768);

    auto tr_z_xxxyyz_xxxxzz = pbuffer.data(idx_dip_ii + 1769);

    auto tr_z_xxxyyz_xxxyyy = pbuffer.data(idx_dip_ii + 1770);

    auto tr_z_xxxyyz_xxxyyz = pbuffer.data(idx_dip_ii + 1771);

    auto tr_z_xxxyyz_xxxyzz = pbuffer.data(idx_dip_ii + 1772);

    auto tr_z_xxxyyz_xxxzzz = pbuffer.data(idx_dip_ii + 1773);

    auto tr_z_xxxyyz_xxyyyy = pbuffer.data(idx_dip_ii + 1774);

    auto tr_z_xxxyyz_xxyyyz = pbuffer.data(idx_dip_ii + 1775);

    auto tr_z_xxxyyz_xxyyzz = pbuffer.data(idx_dip_ii + 1776);

    auto tr_z_xxxyyz_xxyzzz = pbuffer.data(idx_dip_ii + 1777);

    auto tr_z_xxxyyz_xxzzzz = pbuffer.data(idx_dip_ii + 1778);

    auto tr_z_xxxyyz_xyyyyy = pbuffer.data(idx_dip_ii + 1779);

    auto tr_z_xxxyyz_xyyyyz = pbuffer.data(idx_dip_ii + 1780);

    auto tr_z_xxxyyz_xyyyzz = pbuffer.data(idx_dip_ii + 1781);

    auto tr_z_xxxyyz_xyyzzz = pbuffer.data(idx_dip_ii + 1782);

    auto tr_z_xxxyyz_xyzzzz = pbuffer.data(idx_dip_ii + 1783);

    auto tr_z_xxxyyz_xzzzzz = pbuffer.data(idx_dip_ii + 1784);

    auto tr_z_xxxyyz_yyyyyy = pbuffer.data(idx_dip_ii + 1785);

    auto tr_z_xxxyyz_yyyyyz = pbuffer.data(idx_dip_ii + 1786);

    auto tr_z_xxxyyz_yyyyzz = pbuffer.data(idx_dip_ii + 1787);

    auto tr_z_xxxyyz_yyyzzz = pbuffer.data(idx_dip_ii + 1788);

    auto tr_z_xxxyyz_yyzzzz = pbuffer.data(idx_dip_ii + 1789);

    auto tr_z_xxxyyz_yzzzzz = pbuffer.data(idx_dip_ii + 1790);

    auto tr_z_xxxyyz_zzzzzz = pbuffer.data(idx_dip_ii + 1791);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxxyy_xxxxxy, tr_z_xxxyy_xxxxyy, tr_z_xxxyy_xxxyyy, tr_z_xxxyy_xxyyyy, tr_z_xxxyy_xyyyyy, tr_z_xxxyyz_xxxxxx, tr_z_xxxyyz_xxxxxy, tr_z_xxxyyz_xxxxxz, tr_z_xxxyyz_xxxxyy, tr_z_xxxyyz_xxxxyz, tr_z_xxxyyz_xxxxzz, tr_z_xxxyyz_xxxyyy, tr_z_xxxyyz_xxxyyz, tr_z_xxxyyz_xxxyzz, tr_z_xxxyyz_xxxzzz, tr_z_xxxyyz_xxyyyy, tr_z_xxxyyz_xxyyyz, tr_z_xxxyyz_xxyyzz, tr_z_xxxyyz_xxyzzz, tr_z_xxxyyz_xxzzzz, tr_z_xxxyyz_xyyyyy, tr_z_xxxyyz_xyyyyz, tr_z_xxxyyz_xyyyzz, tr_z_xxxyyz_xyyzzz, tr_z_xxxyyz_xyzzzz, tr_z_xxxyyz_xzzzzz, tr_z_xxxyyz_yyyyyy, tr_z_xxxyyz_yyyyyz, tr_z_xxxyyz_yyyyzz, tr_z_xxxyyz_yyyzzz, tr_z_xxxyyz_yyzzzz, tr_z_xxxyyz_yzzzzz, tr_z_xxxyyz_zzzzzz, tr_z_xxxyz_xxxxxx, tr_z_xxxyz_xxxxxz, tr_z_xxxyz_xxxxzz, tr_z_xxxyz_xxxzzz, tr_z_xxxyz_xxzzzz, tr_z_xxxyz_xzzzzz, tr_z_xxxz_xxxxxx, tr_z_xxxz_xxxxxz, tr_z_xxxz_xxxxzz, tr_z_xxxz_xxxzzz, tr_z_xxxz_xxzzzz, tr_z_xxxz_xzzzzz, tr_z_xxyyz_xxxxyz, tr_z_xxyyz_xxxyyz, tr_z_xxyyz_xxxyz, tr_z_xxyyz_xxxyzz, tr_z_xxyyz_xxyyyz, tr_z_xxyyz_xxyyz, tr_z_xxyyz_xxyyzz, tr_z_xxyyz_xxyzz, tr_z_xxyyz_xxyzzz, tr_z_xxyyz_xyyyyz, tr_z_xxyyz_xyyyz, tr_z_xxyyz_xyyyzz, tr_z_xxyyz_xyyzz, tr_z_xxyyz_xyyzzz, tr_z_xxyyz_xyzzz, tr_z_xxyyz_xyzzzz, tr_z_xxyyz_yyyyyy, tr_z_xxyyz_yyyyyz, tr_z_xxyyz_yyyyz, tr_z_xxyyz_yyyyzz, tr_z_xxyyz_yyyzz, tr_z_xxyyz_yyyzzz, tr_z_xxyyz_yyzzz, tr_z_xxyyz_yyzzzz, tr_z_xxyyz_yzzzz, tr_z_xxyyz_yzzzzz, tr_z_xxyyz_zzzzzz, tr_z_xyyz_xxxxyz, tr_z_xyyz_xxxyyz, tr_z_xyyz_xxxyzz, tr_z_xyyz_xxyyyz, tr_z_xyyz_xxyyzz, tr_z_xyyz_xxyzzz, tr_z_xyyz_xyyyyz, tr_z_xyyz_xyyyzz, tr_z_xyyz_xyyzzz, tr_z_xyyz_xyzzzz, tr_z_xyyz_yyyyyy, tr_z_xyyz_yyyyyz, tr_z_xyyz_yyyyzz, tr_z_xyyz_yyyzzz, tr_z_xyyz_yyzzzz, tr_z_xyyz_yzzzzz, tr_z_xyyz_zzzzzz, ts_xxxyy_xxxxxy, ts_xxxyy_xxxxyy, ts_xxxyy_xxxyyy, ts_xxxyy_xxyyyy, ts_xxxyy_xyyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyz_xxxxxx[i] = tr_z_xxxz_xxxxxx[i] * fe_0 + tr_z_xxxyz_xxxxxx[i] * pa_y[i];

        tr_z_xxxyyz_xxxxxy[i] = ts_xxxyy_xxxxxy[i] * fe_0 + tr_z_xxxyy_xxxxxy[i] * pa_z[i];

        tr_z_xxxyyz_xxxxxz[i] = tr_z_xxxz_xxxxxz[i] * fe_0 + tr_z_xxxyz_xxxxxz[i] * pa_y[i];

        tr_z_xxxyyz_xxxxyy[i] = ts_xxxyy_xxxxyy[i] * fe_0 + tr_z_xxxyy_xxxxyy[i] * pa_z[i];

        tr_z_xxxyyz_xxxxyz[i] = 2.0 * tr_z_xyyz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxyyz_xxxyz[i] * fe_0 + tr_z_xxyyz_xxxxyz[i] * pa_x[i];

        tr_z_xxxyyz_xxxxzz[i] = tr_z_xxxz_xxxxzz[i] * fe_0 + tr_z_xxxyz_xxxxzz[i] * pa_y[i];

        tr_z_xxxyyz_xxxyyy[i] = ts_xxxyy_xxxyyy[i] * fe_0 + tr_z_xxxyy_xxxyyy[i] * pa_z[i];

        tr_z_xxxyyz_xxxyyz[i] = 2.0 * tr_z_xyyz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxyyz_xxyyz[i] * fe_0 + tr_z_xxyyz_xxxyyz[i] * pa_x[i];

        tr_z_xxxyyz_xxxyzz[i] = 2.0 * tr_z_xyyz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxyyz_xxyzz[i] * fe_0 + tr_z_xxyyz_xxxyzz[i] * pa_x[i];

        tr_z_xxxyyz_xxxzzz[i] = tr_z_xxxz_xxxzzz[i] * fe_0 + tr_z_xxxyz_xxxzzz[i] * pa_y[i];

        tr_z_xxxyyz_xxyyyy[i] = ts_xxxyy_xxyyyy[i] * fe_0 + tr_z_xxxyy_xxyyyy[i] * pa_z[i];

        tr_z_xxxyyz_xxyyyz[i] = 2.0 * tr_z_xyyz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxyyz_xyyyz[i] * fe_0 + tr_z_xxyyz_xxyyyz[i] * pa_x[i];

        tr_z_xxxyyz_xxyyzz[i] = 2.0 * tr_z_xyyz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxyyz_xyyzz[i] * fe_0 + tr_z_xxyyz_xxyyzz[i] * pa_x[i];

        tr_z_xxxyyz_xxyzzz[i] = 2.0 * tr_z_xyyz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxyyz_xyzzz[i] * fe_0 + tr_z_xxyyz_xxyzzz[i] * pa_x[i];

        tr_z_xxxyyz_xxzzzz[i] = tr_z_xxxz_xxzzzz[i] * fe_0 + tr_z_xxxyz_xxzzzz[i] * pa_y[i];

        tr_z_xxxyyz_xyyyyy[i] = ts_xxxyy_xyyyyy[i] * fe_0 + tr_z_xxxyy_xyyyyy[i] * pa_z[i];

        tr_z_xxxyyz_xyyyyz[i] = 2.0 * tr_z_xyyz_xyyyyz[i] * fe_0 + tr_z_xxyyz_yyyyz[i] * fe_0 + tr_z_xxyyz_xyyyyz[i] * pa_x[i];

        tr_z_xxxyyz_xyyyzz[i] = 2.0 * tr_z_xyyz_xyyyzz[i] * fe_0 + tr_z_xxyyz_yyyzz[i] * fe_0 + tr_z_xxyyz_xyyyzz[i] * pa_x[i];

        tr_z_xxxyyz_xyyzzz[i] = 2.0 * tr_z_xyyz_xyyzzz[i] * fe_0 + tr_z_xxyyz_yyzzz[i] * fe_0 + tr_z_xxyyz_xyyzzz[i] * pa_x[i];

        tr_z_xxxyyz_xyzzzz[i] = 2.0 * tr_z_xyyz_xyzzzz[i] * fe_0 + tr_z_xxyyz_yzzzz[i] * fe_0 + tr_z_xxyyz_xyzzzz[i] * pa_x[i];

        tr_z_xxxyyz_xzzzzz[i] = tr_z_xxxz_xzzzzz[i] * fe_0 + tr_z_xxxyz_xzzzzz[i] * pa_y[i];

        tr_z_xxxyyz_yyyyyy[i] = 2.0 * tr_z_xyyz_yyyyyy[i] * fe_0 + tr_z_xxyyz_yyyyyy[i] * pa_x[i];

        tr_z_xxxyyz_yyyyyz[i] = 2.0 * tr_z_xyyz_yyyyyz[i] * fe_0 + tr_z_xxyyz_yyyyyz[i] * pa_x[i];

        tr_z_xxxyyz_yyyyzz[i] = 2.0 * tr_z_xyyz_yyyyzz[i] * fe_0 + tr_z_xxyyz_yyyyzz[i] * pa_x[i];

        tr_z_xxxyyz_yyyzzz[i] = 2.0 * tr_z_xyyz_yyyzzz[i] * fe_0 + tr_z_xxyyz_yyyzzz[i] * pa_x[i];

        tr_z_xxxyyz_yyzzzz[i] = 2.0 * tr_z_xyyz_yyzzzz[i] * fe_0 + tr_z_xxyyz_yyzzzz[i] * pa_x[i];

        tr_z_xxxyyz_yzzzzz[i] = 2.0 * tr_z_xyyz_yzzzzz[i] * fe_0 + tr_z_xxyyz_yzzzzz[i] * pa_x[i];

        tr_z_xxxyyz_zzzzzz[i] = 2.0 * tr_z_xyyz_zzzzzz[i] * fe_0 + tr_z_xxyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1792-1820 components of targeted buffer : II

    auto tr_z_xxxyzz_xxxxxx = pbuffer.data(idx_dip_ii + 1792);

    auto tr_z_xxxyzz_xxxxxy = pbuffer.data(idx_dip_ii + 1793);

    auto tr_z_xxxyzz_xxxxxz = pbuffer.data(idx_dip_ii + 1794);

    auto tr_z_xxxyzz_xxxxyy = pbuffer.data(idx_dip_ii + 1795);

    auto tr_z_xxxyzz_xxxxyz = pbuffer.data(idx_dip_ii + 1796);

    auto tr_z_xxxyzz_xxxxzz = pbuffer.data(idx_dip_ii + 1797);

    auto tr_z_xxxyzz_xxxyyy = pbuffer.data(idx_dip_ii + 1798);

    auto tr_z_xxxyzz_xxxyyz = pbuffer.data(idx_dip_ii + 1799);

    auto tr_z_xxxyzz_xxxyzz = pbuffer.data(idx_dip_ii + 1800);

    auto tr_z_xxxyzz_xxxzzz = pbuffer.data(idx_dip_ii + 1801);

    auto tr_z_xxxyzz_xxyyyy = pbuffer.data(idx_dip_ii + 1802);

    auto tr_z_xxxyzz_xxyyyz = pbuffer.data(idx_dip_ii + 1803);

    auto tr_z_xxxyzz_xxyyzz = pbuffer.data(idx_dip_ii + 1804);

    auto tr_z_xxxyzz_xxyzzz = pbuffer.data(idx_dip_ii + 1805);

    auto tr_z_xxxyzz_xxzzzz = pbuffer.data(idx_dip_ii + 1806);

    auto tr_z_xxxyzz_xyyyyy = pbuffer.data(idx_dip_ii + 1807);

    auto tr_z_xxxyzz_xyyyyz = pbuffer.data(idx_dip_ii + 1808);

    auto tr_z_xxxyzz_xyyyzz = pbuffer.data(idx_dip_ii + 1809);

    auto tr_z_xxxyzz_xyyzzz = pbuffer.data(idx_dip_ii + 1810);

    auto tr_z_xxxyzz_xyzzzz = pbuffer.data(idx_dip_ii + 1811);

    auto tr_z_xxxyzz_xzzzzz = pbuffer.data(idx_dip_ii + 1812);

    auto tr_z_xxxyzz_yyyyyy = pbuffer.data(idx_dip_ii + 1813);

    auto tr_z_xxxyzz_yyyyyz = pbuffer.data(idx_dip_ii + 1814);

    auto tr_z_xxxyzz_yyyyzz = pbuffer.data(idx_dip_ii + 1815);

    auto tr_z_xxxyzz_yyyzzz = pbuffer.data(idx_dip_ii + 1816);

    auto tr_z_xxxyzz_yyzzzz = pbuffer.data(idx_dip_ii + 1817);

    auto tr_z_xxxyzz_yzzzzz = pbuffer.data(idx_dip_ii + 1818);

    auto tr_z_xxxyzz_zzzzzz = pbuffer.data(idx_dip_ii + 1819);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxyzz_xxxxxx, tr_z_xxxyzz_xxxxxy, tr_z_xxxyzz_xxxxxz, tr_z_xxxyzz_xxxxyy, tr_z_xxxyzz_xxxxyz, tr_z_xxxyzz_xxxxzz, tr_z_xxxyzz_xxxyyy, tr_z_xxxyzz_xxxyyz, tr_z_xxxyzz_xxxyzz, tr_z_xxxyzz_xxxzzz, tr_z_xxxyzz_xxyyyy, tr_z_xxxyzz_xxyyyz, tr_z_xxxyzz_xxyyzz, tr_z_xxxyzz_xxyzzz, tr_z_xxxyzz_xxzzzz, tr_z_xxxyzz_xyyyyy, tr_z_xxxyzz_xyyyyz, tr_z_xxxyzz_xyyyzz, tr_z_xxxyzz_xyyzzz, tr_z_xxxyzz_xyzzzz, tr_z_xxxyzz_xzzzzz, tr_z_xxxyzz_yyyyyy, tr_z_xxxyzz_yyyyyz, tr_z_xxxyzz_yyyyzz, tr_z_xxxyzz_yyyzzz, tr_z_xxxyzz_yyzzzz, tr_z_xxxyzz_yzzzzz, tr_z_xxxyzz_zzzzzz, tr_z_xxxzz_xxxxx, tr_z_xxxzz_xxxxxx, tr_z_xxxzz_xxxxxy, tr_z_xxxzz_xxxxxz, tr_z_xxxzz_xxxxy, tr_z_xxxzz_xxxxyy, tr_z_xxxzz_xxxxyz, tr_z_xxxzz_xxxxz, tr_z_xxxzz_xxxxzz, tr_z_xxxzz_xxxyy, tr_z_xxxzz_xxxyyy, tr_z_xxxzz_xxxyyz, tr_z_xxxzz_xxxyz, tr_z_xxxzz_xxxyzz, tr_z_xxxzz_xxxzz, tr_z_xxxzz_xxxzzz, tr_z_xxxzz_xxyyy, tr_z_xxxzz_xxyyyy, tr_z_xxxzz_xxyyyz, tr_z_xxxzz_xxyyz, tr_z_xxxzz_xxyyzz, tr_z_xxxzz_xxyzz, tr_z_xxxzz_xxyzzz, tr_z_xxxzz_xxzzz, tr_z_xxxzz_xxzzzz, tr_z_xxxzz_xyyyy, tr_z_xxxzz_xyyyyy, tr_z_xxxzz_xyyyyz, tr_z_xxxzz_xyyyz, tr_z_xxxzz_xyyyzz, tr_z_xxxzz_xyyzz, tr_z_xxxzz_xyyzzz, tr_z_xxxzz_xyzzz, tr_z_xxxzz_xyzzzz, tr_z_xxxzz_xzzzz, tr_z_xxxzz_xzzzzz, tr_z_xxxzz_zzzzzz, tr_z_xxyzz_yyyyyy, tr_z_xxyzz_yyyyyz, tr_z_xxyzz_yyyyzz, tr_z_xxyzz_yyyzzz, tr_z_xxyzz_yyzzzz, tr_z_xxyzz_yzzzzz, tr_z_xyzz_yyyyyy, tr_z_xyzz_yyyyyz, tr_z_xyzz_yyyyzz, tr_z_xyzz_yyyzzz, tr_z_xyzz_yyzzzz, tr_z_xyzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyzz_xxxxxx[i] = tr_z_xxxzz_xxxxxx[i] * pa_y[i];

        tr_z_xxxyzz_xxxxxy[i] = tr_z_xxxzz_xxxxx[i] * fe_0 + tr_z_xxxzz_xxxxxy[i] * pa_y[i];

        tr_z_xxxyzz_xxxxxz[i] = tr_z_xxxzz_xxxxxz[i] * pa_y[i];

        tr_z_xxxyzz_xxxxyy[i] = 2.0 * tr_z_xxxzz_xxxxy[i] * fe_0 + tr_z_xxxzz_xxxxyy[i] * pa_y[i];

        tr_z_xxxyzz_xxxxyz[i] = tr_z_xxxzz_xxxxz[i] * fe_0 + tr_z_xxxzz_xxxxyz[i] * pa_y[i];

        tr_z_xxxyzz_xxxxzz[i] = tr_z_xxxzz_xxxxzz[i] * pa_y[i];

        tr_z_xxxyzz_xxxyyy[i] = 3.0 * tr_z_xxxzz_xxxyy[i] * fe_0 + tr_z_xxxzz_xxxyyy[i] * pa_y[i];

        tr_z_xxxyzz_xxxyyz[i] = 2.0 * tr_z_xxxzz_xxxyz[i] * fe_0 + tr_z_xxxzz_xxxyyz[i] * pa_y[i];

        tr_z_xxxyzz_xxxyzz[i] = tr_z_xxxzz_xxxzz[i] * fe_0 + tr_z_xxxzz_xxxyzz[i] * pa_y[i];

        tr_z_xxxyzz_xxxzzz[i] = tr_z_xxxzz_xxxzzz[i] * pa_y[i];

        tr_z_xxxyzz_xxyyyy[i] = 4.0 * tr_z_xxxzz_xxyyy[i] * fe_0 + tr_z_xxxzz_xxyyyy[i] * pa_y[i];

        tr_z_xxxyzz_xxyyyz[i] = 3.0 * tr_z_xxxzz_xxyyz[i] * fe_0 + tr_z_xxxzz_xxyyyz[i] * pa_y[i];

        tr_z_xxxyzz_xxyyzz[i] = 2.0 * tr_z_xxxzz_xxyzz[i] * fe_0 + tr_z_xxxzz_xxyyzz[i] * pa_y[i];

        tr_z_xxxyzz_xxyzzz[i] = tr_z_xxxzz_xxzzz[i] * fe_0 + tr_z_xxxzz_xxyzzz[i] * pa_y[i];

        tr_z_xxxyzz_xxzzzz[i] = tr_z_xxxzz_xxzzzz[i] * pa_y[i];

        tr_z_xxxyzz_xyyyyy[i] = 5.0 * tr_z_xxxzz_xyyyy[i] * fe_0 + tr_z_xxxzz_xyyyyy[i] * pa_y[i];

        tr_z_xxxyzz_xyyyyz[i] = 4.0 * tr_z_xxxzz_xyyyz[i] * fe_0 + tr_z_xxxzz_xyyyyz[i] * pa_y[i];

        tr_z_xxxyzz_xyyyzz[i] = 3.0 * tr_z_xxxzz_xyyzz[i] * fe_0 + tr_z_xxxzz_xyyyzz[i] * pa_y[i];

        tr_z_xxxyzz_xyyzzz[i] = 2.0 * tr_z_xxxzz_xyzzz[i] * fe_0 + tr_z_xxxzz_xyyzzz[i] * pa_y[i];

        tr_z_xxxyzz_xyzzzz[i] = tr_z_xxxzz_xzzzz[i] * fe_0 + tr_z_xxxzz_xyzzzz[i] * pa_y[i];

        tr_z_xxxyzz_xzzzzz[i] = tr_z_xxxzz_xzzzzz[i] * pa_y[i];

        tr_z_xxxyzz_yyyyyy[i] = 2.0 * tr_z_xyzz_yyyyyy[i] * fe_0 + tr_z_xxyzz_yyyyyy[i] * pa_x[i];

        tr_z_xxxyzz_yyyyyz[i] = 2.0 * tr_z_xyzz_yyyyyz[i] * fe_0 + tr_z_xxyzz_yyyyyz[i] * pa_x[i];

        tr_z_xxxyzz_yyyyzz[i] = 2.0 * tr_z_xyzz_yyyyzz[i] * fe_0 + tr_z_xxyzz_yyyyzz[i] * pa_x[i];

        tr_z_xxxyzz_yyyzzz[i] = 2.0 * tr_z_xyzz_yyyzzz[i] * fe_0 + tr_z_xxyzz_yyyzzz[i] * pa_x[i];

        tr_z_xxxyzz_yyzzzz[i] = 2.0 * tr_z_xyzz_yyzzzz[i] * fe_0 + tr_z_xxyzz_yyzzzz[i] * pa_x[i];

        tr_z_xxxyzz_yzzzzz[i] = 2.0 * tr_z_xyzz_yzzzzz[i] * fe_0 + tr_z_xxyzz_yzzzzz[i] * pa_x[i];

        tr_z_xxxyzz_zzzzzz[i] = tr_z_xxxzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1820-1848 components of targeted buffer : II

    auto tr_z_xxxzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1820);

    auto tr_z_xxxzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1821);

    auto tr_z_xxxzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1822);

    auto tr_z_xxxzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1823);

    auto tr_z_xxxzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1824);

    auto tr_z_xxxzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1825);

    auto tr_z_xxxzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1826);

    auto tr_z_xxxzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1827);

    auto tr_z_xxxzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1828);

    auto tr_z_xxxzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1829);

    auto tr_z_xxxzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1830);

    auto tr_z_xxxzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1831);

    auto tr_z_xxxzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1832);

    auto tr_z_xxxzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1833);

    auto tr_z_xxxzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1834);

    auto tr_z_xxxzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1835);

    auto tr_z_xxxzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1836);

    auto tr_z_xxxzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1837);

    auto tr_z_xxxzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1838);

    auto tr_z_xxxzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1839);

    auto tr_z_xxxzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1840);

    auto tr_z_xxxzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1841);

    auto tr_z_xxxzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1842);

    auto tr_z_xxxzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1843);

    auto tr_z_xxxzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1844);

    auto tr_z_xxxzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1845);

    auto tr_z_xxxzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1846);

    auto tr_z_xxxzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1847);

    #pragma omp simd aligned(pa_x, tr_z_xxxzzz_xxxxxx, tr_z_xxxzzz_xxxxxy, tr_z_xxxzzz_xxxxxz, tr_z_xxxzzz_xxxxyy, tr_z_xxxzzz_xxxxyz, tr_z_xxxzzz_xxxxzz, tr_z_xxxzzz_xxxyyy, tr_z_xxxzzz_xxxyyz, tr_z_xxxzzz_xxxyzz, tr_z_xxxzzz_xxxzzz, tr_z_xxxzzz_xxyyyy, tr_z_xxxzzz_xxyyyz, tr_z_xxxzzz_xxyyzz, tr_z_xxxzzz_xxyzzz, tr_z_xxxzzz_xxzzzz, tr_z_xxxzzz_xyyyyy, tr_z_xxxzzz_xyyyyz, tr_z_xxxzzz_xyyyzz, tr_z_xxxzzz_xyyzzz, tr_z_xxxzzz_xyzzzz, tr_z_xxxzzz_xzzzzz, tr_z_xxxzzz_yyyyyy, tr_z_xxxzzz_yyyyyz, tr_z_xxxzzz_yyyyzz, tr_z_xxxzzz_yyyzzz, tr_z_xxxzzz_yyzzzz, tr_z_xxxzzz_yzzzzz, tr_z_xxxzzz_zzzzzz, tr_z_xxzzz_xxxxx, tr_z_xxzzz_xxxxxx, tr_z_xxzzz_xxxxxy, tr_z_xxzzz_xxxxxz, tr_z_xxzzz_xxxxy, tr_z_xxzzz_xxxxyy, tr_z_xxzzz_xxxxyz, tr_z_xxzzz_xxxxz, tr_z_xxzzz_xxxxzz, tr_z_xxzzz_xxxyy, tr_z_xxzzz_xxxyyy, tr_z_xxzzz_xxxyyz, tr_z_xxzzz_xxxyz, tr_z_xxzzz_xxxyzz, tr_z_xxzzz_xxxzz, tr_z_xxzzz_xxxzzz, tr_z_xxzzz_xxyyy, tr_z_xxzzz_xxyyyy, tr_z_xxzzz_xxyyyz, tr_z_xxzzz_xxyyz, tr_z_xxzzz_xxyyzz, tr_z_xxzzz_xxyzz, tr_z_xxzzz_xxyzzz, tr_z_xxzzz_xxzzz, tr_z_xxzzz_xxzzzz, tr_z_xxzzz_xyyyy, tr_z_xxzzz_xyyyyy, tr_z_xxzzz_xyyyyz, tr_z_xxzzz_xyyyz, tr_z_xxzzz_xyyyzz, tr_z_xxzzz_xyyzz, tr_z_xxzzz_xyyzzz, tr_z_xxzzz_xyzzz, tr_z_xxzzz_xyzzzz, tr_z_xxzzz_xzzzz, tr_z_xxzzz_xzzzzz, tr_z_xxzzz_yyyyy, tr_z_xxzzz_yyyyyy, tr_z_xxzzz_yyyyyz, tr_z_xxzzz_yyyyz, tr_z_xxzzz_yyyyzz, tr_z_xxzzz_yyyzz, tr_z_xxzzz_yyyzzz, tr_z_xxzzz_yyzzz, tr_z_xxzzz_yyzzzz, tr_z_xxzzz_yzzzz, tr_z_xxzzz_yzzzzz, tr_z_xxzzz_zzzzz, tr_z_xxzzz_zzzzzz, tr_z_xzzz_xxxxxx, tr_z_xzzz_xxxxxy, tr_z_xzzz_xxxxxz, tr_z_xzzz_xxxxyy, tr_z_xzzz_xxxxyz, tr_z_xzzz_xxxxzz, tr_z_xzzz_xxxyyy, tr_z_xzzz_xxxyyz, tr_z_xzzz_xxxyzz, tr_z_xzzz_xxxzzz, tr_z_xzzz_xxyyyy, tr_z_xzzz_xxyyyz, tr_z_xzzz_xxyyzz, tr_z_xzzz_xxyzzz, tr_z_xzzz_xxzzzz, tr_z_xzzz_xyyyyy, tr_z_xzzz_xyyyyz, tr_z_xzzz_xyyyzz, tr_z_xzzz_xyyzzz, tr_z_xzzz_xyzzzz, tr_z_xzzz_xzzzzz, tr_z_xzzz_yyyyyy, tr_z_xzzz_yyyyyz, tr_z_xzzz_yyyyzz, tr_z_xzzz_yyyzzz, tr_z_xzzz_yyzzzz, tr_z_xzzz_yzzzzz, tr_z_xzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzzz_xxxxxx[i] = 2.0 * tr_z_xzzz_xxxxxx[i] * fe_0 + 6.0 * tr_z_xxzzz_xxxxx[i] * fe_0 + tr_z_xxzzz_xxxxxx[i] * pa_x[i];

        tr_z_xxxzzz_xxxxxy[i] = 2.0 * tr_z_xzzz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxzzz_xxxxy[i] * fe_0 + tr_z_xxzzz_xxxxxy[i] * pa_x[i];

        tr_z_xxxzzz_xxxxxz[i] = 2.0 * tr_z_xzzz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxzzz_xxxxz[i] * fe_0 + tr_z_xxzzz_xxxxxz[i] * pa_x[i];

        tr_z_xxxzzz_xxxxyy[i] = 2.0 * tr_z_xzzz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxzzz_xxxyy[i] * fe_0 + tr_z_xxzzz_xxxxyy[i] * pa_x[i];

        tr_z_xxxzzz_xxxxyz[i] = 2.0 * tr_z_xzzz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxzzz_xxxyz[i] * fe_0 + tr_z_xxzzz_xxxxyz[i] * pa_x[i];

        tr_z_xxxzzz_xxxxzz[i] = 2.0 * tr_z_xzzz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxzzz_xxxzz[i] * fe_0 + tr_z_xxzzz_xxxxzz[i] * pa_x[i];

        tr_z_xxxzzz_xxxyyy[i] = 2.0 * tr_z_xzzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxzzz_xxyyy[i] * fe_0 + tr_z_xxzzz_xxxyyy[i] * pa_x[i];

        tr_z_xxxzzz_xxxyyz[i] = 2.0 * tr_z_xzzz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxzzz_xxyyz[i] * fe_0 + tr_z_xxzzz_xxxyyz[i] * pa_x[i];

        tr_z_xxxzzz_xxxyzz[i] = 2.0 * tr_z_xzzz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxzzz_xxyzz[i] * fe_0 + tr_z_xxzzz_xxxyzz[i] * pa_x[i];

        tr_z_xxxzzz_xxxzzz[i] = 2.0 * tr_z_xzzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxzzz_xxzzz[i] * fe_0 + tr_z_xxzzz_xxxzzz[i] * pa_x[i];

        tr_z_xxxzzz_xxyyyy[i] = 2.0 * tr_z_xzzz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxzzz_xyyyy[i] * fe_0 + tr_z_xxzzz_xxyyyy[i] * pa_x[i];

        tr_z_xxxzzz_xxyyyz[i] = 2.0 * tr_z_xzzz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxzzz_xyyyz[i] * fe_0 + tr_z_xxzzz_xxyyyz[i] * pa_x[i];

        tr_z_xxxzzz_xxyyzz[i] = 2.0 * tr_z_xzzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxzzz_xyyzz[i] * fe_0 + tr_z_xxzzz_xxyyzz[i] * pa_x[i];

        tr_z_xxxzzz_xxyzzz[i] = 2.0 * tr_z_xzzz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxzzz_xyzzz[i] * fe_0 + tr_z_xxzzz_xxyzzz[i] * pa_x[i];

        tr_z_xxxzzz_xxzzzz[i] = 2.0 * tr_z_xzzz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxzzz_xzzzz[i] * fe_0 + tr_z_xxzzz_xxzzzz[i] * pa_x[i];

        tr_z_xxxzzz_xyyyyy[i] = 2.0 * tr_z_xzzz_xyyyyy[i] * fe_0 + tr_z_xxzzz_yyyyy[i] * fe_0 + tr_z_xxzzz_xyyyyy[i] * pa_x[i];

        tr_z_xxxzzz_xyyyyz[i] = 2.0 * tr_z_xzzz_xyyyyz[i] * fe_0 + tr_z_xxzzz_yyyyz[i] * fe_0 + tr_z_xxzzz_xyyyyz[i] * pa_x[i];

        tr_z_xxxzzz_xyyyzz[i] = 2.0 * tr_z_xzzz_xyyyzz[i] * fe_0 + tr_z_xxzzz_yyyzz[i] * fe_0 + tr_z_xxzzz_xyyyzz[i] * pa_x[i];

        tr_z_xxxzzz_xyyzzz[i] = 2.0 * tr_z_xzzz_xyyzzz[i] * fe_0 + tr_z_xxzzz_yyzzz[i] * fe_0 + tr_z_xxzzz_xyyzzz[i] * pa_x[i];

        tr_z_xxxzzz_xyzzzz[i] = 2.0 * tr_z_xzzz_xyzzzz[i] * fe_0 + tr_z_xxzzz_yzzzz[i] * fe_0 + tr_z_xxzzz_xyzzzz[i] * pa_x[i];

        tr_z_xxxzzz_xzzzzz[i] = 2.0 * tr_z_xzzz_xzzzzz[i] * fe_0 + tr_z_xxzzz_zzzzz[i] * fe_0 + tr_z_xxzzz_xzzzzz[i] * pa_x[i];

        tr_z_xxxzzz_yyyyyy[i] = 2.0 * tr_z_xzzz_yyyyyy[i] * fe_0 + tr_z_xxzzz_yyyyyy[i] * pa_x[i];

        tr_z_xxxzzz_yyyyyz[i] = 2.0 * tr_z_xzzz_yyyyyz[i] * fe_0 + tr_z_xxzzz_yyyyyz[i] * pa_x[i];

        tr_z_xxxzzz_yyyyzz[i] = 2.0 * tr_z_xzzz_yyyyzz[i] * fe_0 + tr_z_xxzzz_yyyyzz[i] * pa_x[i];

        tr_z_xxxzzz_yyyzzz[i] = 2.0 * tr_z_xzzz_yyyzzz[i] * fe_0 + tr_z_xxzzz_yyyzzz[i] * pa_x[i];

        tr_z_xxxzzz_yyzzzz[i] = 2.0 * tr_z_xzzz_yyzzzz[i] * fe_0 + tr_z_xxzzz_yyzzzz[i] * pa_x[i];

        tr_z_xxxzzz_yzzzzz[i] = 2.0 * tr_z_xzzz_yzzzzz[i] * fe_0 + tr_z_xxzzz_yzzzzz[i] * pa_x[i];

        tr_z_xxxzzz_zzzzzz[i] = 2.0 * tr_z_xzzz_zzzzzz[i] * fe_0 + tr_z_xxzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1848-1876 components of targeted buffer : II

    auto tr_z_xxyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 1848);

    auto tr_z_xxyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 1849);

    auto tr_z_xxyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 1850);

    auto tr_z_xxyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 1851);

    auto tr_z_xxyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 1852);

    auto tr_z_xxyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 1853);

    auto tr_z_xxyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 1854);

    auto tr_z_xxyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 1855);

    auto tr_z_xxyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 1856);

    auto tr_z_xxyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 1857);

    auto tr_z_xxyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 1858);

    auto tr_z_xxyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 1859);

    auto tr_z_xxyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 1860);

    auto tr_z_xxyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 1861);

    auto tr_z_xxyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 1862);

    auto tr_z_xxyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 1863);

    auto tr_z_xxyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 1864);

    auto tr_z_xxyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 1865);

    auto tr_z_xxyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 1866);

    auto tr_z_xxyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 1867);

    auto tr_z_xxyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 1868);

    auto tr_z_xxyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 1869);

    auto tr_z_xxyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 1870);

    auto tr_z_xxyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 1871);

    auto tr_z_xxyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 1872);

    auto tr_z_xxyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 1873);

    auto tr_z_xxyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 1874);

    auto tr_z_xxyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 1875);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyy_xxxxxx, tr_z_xxyy_xxxxxz, tr_z_xxyy_xxxxzz, tr_z_xxyy_xxxzzz, tr_z_xxyy_xxzzzz, tr_z_xxyy_xzzzzz, tr_z_xxyyy_xxxxxx, tr_z_xxyyy_xxxxxz, tr_z_xxyyy_xxxxzz, tr_z_xxyyy_xxxzzz, tr_z_xxyyy_xxzzzz, tr_z_xxyyy_xzzzzz, tr_z_xxyyyy_xxxxxx, tr_z_xxyyyy_xxxxxy, tr_z_xxyyyy_xxxxxz, tr_z_xxyyyy_xxxxyy, tr_z_xxyyyy_xxxxyz, tr_z_xxyyyy_xxxxzz, tr_z_xxyyyy_xxxyyy, tr_z_xxyyyy_xxxyyz, tr_z_xxyyyy_xxxyzz, tr_z_xxyyyy_xxxzzz, tr_z_xxyyyy_xxyyyy, tr_z_xxyyyy_xxyyyz, tr_z_xxyyyy_xxyyzz, tr_z_xxyyyy_xxyzzz, tr_z_xxyyyy_xxzzzz, tr_z_xxyyyy_xyyyyy, tr_z_xxyyyy_xyyyyz, tr_z_xxyyyy_xyyyzz, tr_z_xxyyyy_xyyzzz, tr_z_xxyyyy_xyzzzz, tr_z_xxyyyy_xzzzzz, tr_z_xxyyyy_yyyyyy, tr_z_xxyyyy_yyyyyz, tr_z_xxyyyy_yyyyzz, tr_z_xxyyyy_yyyzzz, tr_z_xxyyyy_yyzzzz, tr_z_xxyyyy_yzzzzz, tr_z_xxyyyy_zzzzzz, tr_z_xyyyy_xxxxxy, tr_z_xyyyy_xxxxy, tr_z_xyyyy_xxxxyy, tr_z_xyyyy_xxxxyz, tr_z_xyyyy_xxxyy, tr_z_xyyyy_xxxyyy, tr_z_xyyyy_xxxyyz, tr_z_xyyyy_xxxyz, tr_z_xyyyy_xxxyzz, tr_z_xyyyy_xxyyy, tr_z_xyyyy_xxyyyy, tr_z_xyyyy_xxyyyz, tr_z_xyyyy_xxyyz, tr_z_xyyyy_xxyyzz, tr_z_xyyyy_xxyzz, tr_z_xyyyy_xxyzzz, tr_z_xyyyy_xyyyy, tr_z_xyyyy_xyyyyy, tr_z_xyyyy_xyyyyz, tr_z_xyyyy_xyyyz, tr_z_xyyyy_xyyyzz, tr_z_xyyyy_xyyzz, tr_z_xyyyy_xyyzzz, tr_z_xyyyy_xyzzz, tr_z_xyyyy_xyzzzz, tr_z_xyyyy_yyyyy, tr_z_xyyyy_yyyyyy, tr_z_xyyyy_yyyyyz, tr_z_xyyyy_yyyyz, tr_z_xyyyy_yyyyzz, tr_z_xyyyy_yyyzz, tr_z_xyyyy_yyyzzz, tr_z_xyyyy_yyzzz, tr_z_xyyyy_yyzzzz, tr_z_xyyyy_yzzzz, tr_z_xyyyy_yzzzzz, tr_z_xyyyy_zzzzzz, tr_z_yyyy_xxxxxy, tr_z_yyyy_xxxxyy, tr_z_yyyy_xxxxyz, tr_z_yyyy_xxxyyy, tr_z_yyyy_xxxyyz, tr_z_yyyy_xxxyzz, tr_z_yyyy_xxyyyy, tr_z_yyyy_xxyyyz, tr_z_yyyy_xxyyzz, tr_z_yyyy_xxyzzz, tr_z_yyyy_xyyyyy, tr_z_yyyy_xyyyyz, tr_z_yyyy_xyyyzz, tr_z_yyyy_xyyzzz, tr_z_yyyy_xyzzzz, tr_z_yyyy_yyyyyy, tr_z_yyyy_yyyyyz, tr_z_yyyy_yyyyzz, tr_z_yyyy_yyyzzz, tr_z_yyyy_yyzzzz, tr_z_yyyy_yzzzzz, tr_z_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyy_xxxxxx[i] = 3.0 * tr_z_xxyy_xxxxxx[i] * fe_0 + tr_z_xxyyy_xxxxxx[i] * pa_y[i];

        tr_z_xxyyyy_xxxxxy[i] = tr_z_yyyy_xxxxxy[i] * fe_0 + 5.0 * tr_z_xyyyy_xxxxy[i] * fe_0 + tr_z_xyyyy_xxxxxy[i] * pa_x[i];

        tr_z_xxyyyy_xxxxxz[i] = 3.0 * tr_z_xxyy_xxxxxz[i] * fe_0 + tr_z_xxyyy_xxxxxz[i] * pa_y[i];

        tr_z_xxyyyy_xxxxyy[i] = tr_z_yyyy_xxxxyy[i] * fe_0 + 4.0 * tr_z_xyyyy_xxxyy[i] * fe_0 + tr_z_xyyyy_xxxxyy[i] * pa_x[i];

        tr_z_xxyyyy_xxxxyz[i] = tr_z_yyyy_xxxxyz[i] * fe_0 + 4.0 * tr_z_xyyyy_xxxyz[i] * fe_0 + tr_z_xyyyy_xxxxyz[i] * pa_x[i];

        tr_z_xxyyyy_xxxxzz[i] = 3.0 * tr_z_xxyy_xxxxzz[i] * fe_0 + tr_z_xxyyy_xxxxzz[i] * pa_y[i];

        tr_z_xxyyyy_xxxyyy[i] = tr_z_yyyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_xyyyy_xxyyy[i] * fe_0 + tr_z_xyyyy_xxxyyy[i] * pa_x[i];

        tr_z_xxyyyy_xxxyyz[i] = tr_z_yyyy_xxxyyz[i] * fe_0 + 3.0 * tr_z_xyyyy_xxyyz[i] * fe_0 + tr_z_xyyyy_xxxyyz[i] * pa_x[i];

        tr_z_xxyyyy_xxxyzz[i] = tr_z_yyyy_xxxyzz[i] * fe_0 + 3.0 * tr_z_xyyyy_xxyzz[i] * fe_0 + tr_z_xyyyy_xxxyzz[i] * pa_x[i];

        tr_z_xxyyyy_xxxzzz[i] = 3.0 * tr_z_xxyy_xxxzzz[i] * fe_0 + tr_z_xxyyy_xxxzzz[i] * pa_y[i];

        tr_z_xxyyyy_xxyyyy[i] = tr_z_yyyy_xxyyyy[i] * fe_0 + 2.0 * tr_z_xyyyy_xyyyy[i] * fe_0 + tr_z_xyyyy_xxyyyy[i] * pa_x[i];

        tr_z_xxyyyy_xxyyyz[i] = tr_z_yyyy_xxyyyz[i] * fe_0 + 2.0 * tr_z_xyyyy_xyyyz[i] * fe_0 + tr_z_xyyyy_xxyyyz[i] * pa_x[i];

        tr_z_xxyyyy_xxyyzz[i] = tr_z_yyyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_xyyyy_xyyzz[i] * fe_0 + tr_z_xyyyy_xxyyzz[i] * pa_x[i];

        tr_z_xxyyyy_xxyzzz[i] = tr_z_yyyy_xxyzzz[i] * fe_0 + 2.0 * tr_z_xyyyy_xyzzz[i] * fe_0 + tr_z_xyyyy_xxyzzz[i] * pa_x[i];

        tr_z_xxyyyy_xxzzzz[i] = 3.0 * tr_z_xxyy_xxzzzz[i] * fe_0 + tr_z_xxyyy_xxzzzz[i] * pa_y[i];

        tr_z_xxyyyy_xyyyyy[i] = tr_z_yyyy_xyyyyy[i] * fe_0 + tr_z_xyyyy_yyyyy[i] * fe_0 + tr_z_xyyyy_xyyyyy[i] * pa_x[i];

        tr_z_xxyyyy_xyyyyz[i] = tr_z_yyyy_xyyyyz[i] * fe_0 + tr_z_xyyyy_yyyyz[i] * fe_0 + tr_z_xyyyy_xyyyyz[i] * pa_x[i];

        tr_z_xxyyyy_xyyyzz[i] = tr_z_yyyy_xyyyzz[i] * fe_0 + tr_z_xyyyy_yyyzz[i] * fe_0 + tr_z_xyyyy_xyyyzz[i] * pa_x[i];

        tr_z_xxyyyy_xyyzzz[i] = tr_z_yyyy_xyyzzz[i] * fe_0 + tr_z_xyyyy_yyzzz[i] * fe_0 + tr_z_xyyyy_xyyzzz[i] * pa_x[i];

        tr_z_xxyyyy_xyzzzz[i] = tr_z_yyyy_xyzzzz[i] * fe_0 + tr_z_xyyyy_yzzzz[i] * fe_0 + tr_z_xyyyy_xyzzzz[i] * pa_x[i];

        tr_z_xxyyyy_xzzzzz[i] = 3.0 * tr_z_xxyy_xzzzzz[i] * fe_0 + tr_z_xxyyy_xzzzzz[i] * pa_y[i];

        tr_z_xxyyyy_yyyyyy[i] = tr_z_yyyy_yyyyyy[i] * fe_0 + tr_z_xyyyy_yyyyyy[i] * pa_x[i];

        tr_z_xxyyyy_yyyyyz[i] = tr_z_yyyy_yyyyyz[i] * fe_0 + tr_z_xyyyy_yyyyyz[i] * pa_x[i];

        tr_z_xxyyyy_yyyyzz[i] = tr_z_yyyy_yyyyzz[i] * fe_0 + tr_z_xyyyy_yyyyzz[i] * pa_x[i];

        tr_z_xxyyyy_yyyzzz[i] = tr_z_yyyy_yyyzzz[i] * fe_0 + tr_z_xyyyy_yyyzzz[i] * pa_x[i];

        tr_z_xxyyyy_yyzzzz[i] = tr_z_yyyy_yyzzzz[i] * fe_0 + tr_z_xyyyy_yyzzzz[i] * pa_x[i];

        tr_z_xxyyyy_yzzzzz[i] = tr_z_yyyy_yzzzzz[i] * fe_0 + tr_z_xyyyy_yzzzzz[i] * pa_x[i];

        tr_z_xxyyyy_zzzzzz[i] = tr_z_yyyy_zzzzzz[i] * fe_0 + tr_z_xyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1876-1904 components of targeted buffer : II

    auto tr_z_xxyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 1876);

    auto tr_z_xxyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 1877);

    auto tr_z_xxyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 1878);

    auto tr_z_xxyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 1879);

    auto tr_z_xxyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 1880);

    auto tr_z_xxyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 1881);

    auto tr_z_xxyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 1882);

    auto tr_z_xxyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 1883);

    auto tr_z_xxyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 1884);

    auto tr_z_xxyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 1885);

    auto tr_z_xxyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 1886);

    auto tr_z_xxyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 1887);

    auto tr_z_xxyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 1888);

    auto tr_z_xxyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 1889);

    auto tr_z_xxyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 1890);

    auto tr_z_xxyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 1891);

    auto tr_z_xxyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 1892);

    auto tr_z_xxyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 1893);

    auto tr_z_xxyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 1894);

    auto tr_z_xxyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 1895);

    auto tr_z_xxyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 1896);

    auto tr_z_xxyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 1897);

    auto tr_z_xxyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 1898);

    auto tr_z_xxyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 1899);

    auto tr_z_xxyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 1900);

    auto tr_z_xxyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 1901);

    auto tr_z_xxyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 1902);

    auto tr_z_xxyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 1903);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxyyy_xxxxxy, tr_z_xxyyy_xxxxyy, tr_z_xxyyy_xxxyyy, tr_z_xxyyy_xxyyyy, tr_z_xxyyy_xyyyyy, tr_z_xxyyyz_xxxxxx, tr_z_xxyyyz_xxxxxy, tr_z_xxyyyz_xxxxxz, tr_z_xxyyyz_xxxxyy, tr_z_xxyyyz_xxxxyz, tr_z_xxyyyz_xxxxzz, tr_z_xxyyyz_xxxyyy, tr_z_xxyyyz_xxxyyz, tr_z_xxyyyz_xxxyzz, tr_z_xxyyyz_xxxzzz, tr_z_xxyyyz_xxyyyy, tr_z_xxyyyz_xxyyyz, tr_z_xxyyyz_xxyyzz, tr_z_xxyyyz_xxyzzz, tr_z_xxyyyz_xxzzzz, tr_z_xxyyyz_xyyyyy, tr_z_xxyyyz_xyyyyz, tr_z_xxyyyz_xyyyzz, tr_z_xxyyyz_xyyzzz, tr_z_xxyyyz_xyzzzz, tr_z_xxyyyz_xzzzzz, tr_z_xxyyyz_yyyyyy, tr_z_xxyyyz_yyyyyz, tr_z_xxyyyz_yyyyzz, tr_z_xxyyyz_yyyzzz, tr_z_xxyyyz_yyzzzz, tr_z_xxyyyz_yzzzzz, tr_z_xxyyyz_zzzzzz, tr_z_xxyyz_xxxxxx, tr_z_xxyyz_xxxxxz, tr_z_xxyyz_xxxxzz, tr_z_xxyyz_xxxzzz, tr_z_xxyyz_xxzzzz, tr_z_xxyyz_xzzzzz, tr_z_xxyz_xxxxxx, tr_z_xxyz_xxxxxz, tr_z_xxyz_xxxxzz, tr_z_xxyz_xxxzzz, tr_z_xxyz_xxzzzz, tr_z_xxyz_xzzzzz, tr_z_xyyyz_xxxxyz, tr_z_xyyyz_xxxyyz, tr_z_xyyyz_xxxyz, tr_z_xyyyz_xxxyzz, tr_z_xyyyz_xxyyyz, tr_z_xyyyz_xxyyz, tr_z_xyyyz_xxyyzz, tr_z_xyyyz_xxyzz, tr_z_xyyyz_xxyzzz, tr_z_xyyyz_xyyyyz, tr_z_xyyyz_xyyyz, tr_z_xyyyz_xyyyzz, tr_z_xyyyz_xyyzz, tr_z_xyyyz_xyyzzz, tr_z_xyyyz_xyzzz, tr_z_xyyyz_xyzzzz, tr_z_xyyyz_yyyyyy, tr_z_xyyyz_yyyyyz, tr_z_xyyyz_yyyyz, tr_z_xyyyz_yyyyzz, tr_z_xyyyz_yyyzz, tr_z_xyyyz_yyyzzz, tr_z_xyyyz_yyzzz, tr_z_xyyyz_yyzzzz, tr_z_xyyyz_yzzzz, tr_z_xyyyz_yzzzzz, tr_z_xyyyz_zzzzzz, tr_z_yyyz_xxxxyz, tr_z_yyyz_xxxyyz, tr_z_yyyz_xxxyzz, tr_z_yyyz_xxyyyz, tr_z_yyyz_xxyyzz, tr_z_yyyz_xxyzzz, tr_z_yyyz_xyyyyz, tr_z_yyyz_xyyyzz, tr_z_yyyz_xyyzzz, tr_z_yyyz_xyzzzz, tr_z_yyyz_yyyyyy, tr_z_yyyz_yyyyyz, tr_z_yyyz_yyyyzz, tr_z_yyyz_yyyzzz, tr_z_yyyz_yyzzzz, tr_z_yyyz_yzzzzz, tr_z_yyyz_zzzzzz, ts_xxyyy_xxxxxy, ts_xxyyy_xxxxyy, ts_xxyyy_xxxyyy, ts_xxyyy_xxyyyy, ts_xxyyy_xyyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyz_xxxxxx[i] = 2.0 * tr_z_xxyz_xxxxxx[i] * fe_0 + tr_z_xxyyz_xxxxxx[i] * pa_y[i];

        tr_z_xxyyyz_xxxxxy[i] = ts_xxyyy_xxxxxy[i] * fe_0 + tr_z_xxyyy_xxxxxy[i] * pa_z[i];

        tr_z_xxyyyz_xxxxxz[i] = 2.0 * tr_z_xxyz_xxxxxz[i] * fe_0 + tr_z_xxyyz_xxxxxz[i] * pa_y[i];

        tr_z_xxyyyz_xxxxyy[i] = ts_xxyyy_xxxxyy[i] * fe_0 + tr_z_xxyyy_xxxxyy[i] * pa_z[i];

        tr_z_xxyyyz_xxxxyz[i] = tr_z_yyyz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xyyyz_xxxyz[i] * fe_0 + tr_z_xyyyz_xxxxyz[i] * pa_x[i];

        tr_z_xxyyyz_xxxxzz[i] = 2.0 * tr_z_xxyz_xxxxzz[i] * fe_0 + tr_z_xxyyz_xxxxzz[i] * pa_y[i];

        tr_z_xxyyyz_xxxyyy[i] = ts_xxyyy_xxxyyy[i] * fe_0 + tr_z_xxyyy_xxxyyy[i] * pa_z[i];

        tr_z_xxyyyz_xxxyyz[i] = tr_z_yyyz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xyyyz_xxyyz[i] * fe_0 + tr_z_xyyyz_xxxyyz[i] * pa_x[i];

        tr_z_xxyyyz_xxxyzz[i] = tr_z_yyyz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xyyyz_xxyzz[i] * fe_0 + tr_z_xyyyz_xxxyzz[i] * pa_x[i];

        tr_z_xxyyyz_xxxzzz[i] = 2.0 * tr_z_xxyz_xxxzzz[i] * fe_0 + tr_z_xxyyz_xxxzzz[i] * pa_y[i];

        tr_z_xxyyyz_xxyyyy[i] = ts_xxyyy_xxyyyy[i] * fe_0 + tr_z_xxyyy_xxyyyy[i] * pa_z[i];

        tr_z_xxyyyz_xxyyyz[i] = tr_z_yyyz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xyyyz_xyyyz[i] * fe_0 + tr_z_xyyyz_xxyyyz[i] * pa_x[i];

        tr_z_xxyyyz_xxyyzz[i] = tr_z_yyyz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xyyyz_xyyzz[i] * fe_0 + tr_z_xyyyz_xxyyzz[i] * pa_x[i];

        tr_z_xxyyyz_xxyzzz[i] = tr_z_yyyz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xyyyz_xyzzz[i] * fe_0 + tr_z_xyyyz_xxyzzz[i] * pa_x[i];

        tr_z_xxyyyz_xxzzzz[i] = 2.0 * tr_z_xxyz_xxzzzz[i] * fe_0 + tr_z_xxyyz_xxzzzz[i] * pa_y[i];

        tr_z_xxyyyz_xyyyyy[i] = ts_xxyyy_xyyyyy[i] * fe_0 + tr_z_xxyyy_xyyyyy[i] * pa_z[i];

        tr_z_xxyyyz_xyyyyz[i] = tr_z_yyyz_xyyyyz[i] * fe_0 + tr_z_xyyyz_yyyyz[i] * fe_0 + tr_z_xyyyz_xyyyyz[i] * pa_x[i];

        tr_z_xxyyyz_xyyyzz[i] = tr_z_yyyz_xyyyzz[i] * fe_0 + tr_z_xyyyz_yyyzz[i] * fe_0 + tr_z_xyyyz_xyyyzz[i] * pa_x[i];

        tr_z_xxyyyz_xyyzzz[i] = tr_z_yyyz_xyyzzz[i] * fe_0 + tr_z_xyyyz_yyzzz[i] * fe_0 + tr_z_xyyyz_xyyzzz[i] * pa_x[i];

        tr_z_xxyyyz_xyzzzz[i] = tr_z_yyyz_xyzzzz[i] * fe_0 + tr_z_xyyyz_yzzzz[i] * fe_0 + tr_z_xyyyz_xyzzzz[i] * pa_x[i];

        tr_z_xxyyyz_xzzzzz[i] = 2.0 * tr_z_xxyz_xzzzzz[i] * fe_0 + tr_z_xxyyz_xzzzzz[i] * pa_y[i];

        tr_z_xxyyyz_yyyyyy[i] = tr_z_yyyz_yyyyyy[i] * fe_0 + tr_z_xyyyz_yyyyyy[i] * pa_x[i];

        tr_z_xxyyyz_yyyyyz[i] = tr_z_yyyz_yyyyyz[i] * fe_0 + tr_z_xyyyz_yyyyyz[i] * pa_x[i];

        tr_z_xxyyyz_yyyyzz[i] = tr_z_yyyz_yyyyzz[i] * fe_0 + tr_z_xyyyz_yyyyzz[i] * pa_x[i];

        tr_z_xxyyyz_yyyzzz[i] = tr_z_yyyz_yyyzzz[i] * fe_0 + tr_z_xyyyz_yyyzzz[i] * pa_x[i];

        tr_z_xxyyyz_yyzzzz[i] = tr_z_yyyz_yyzzzz[i] * fe_0 + tr_z_xyyyz_yyzzzz[i] * pa_x[i];

        tr_z_xxyyyz_yzzzzz[i] = tr_z_yyyz_yzzzzz[i] * fe_0 + tr_z_xyyyz_yzzzzz[i] * pa_x[i];

        tr_z_xxyyyz_zzzzzz[i] = tr_z_yyyz_zzzzzz[i] * fe_0 + tr_z_xyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1904-1932 components of targeted buffer : II

    auto tr_z_xxyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 1904);

    auto tr_z_xxyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 1905);

    auto tr_z_xxyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 1906);

    auto tr_z_xxyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 1907);

    auto tr_z_xxyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 1908);

    auto tr_z_xxyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 1909);

    auto tr_z_xxyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 1910);

    auto tr_z_xxyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 1911);

    auto tr_z_xxyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 1912);

    auto tr_z_xxyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 1913);

    auto tr_z_xxyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 1914);

    auto tr_z_xxyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 1915);

    auto tr_z_xxyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 1916);

    auto tr_z_xxyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 1917);

    auto tr_z_xxyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 1918);

    auto tr_z_xxyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 1919);

    auto tr_z_xxyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 1920);

    auto tr_z_xxyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 1921);

    auto tr_z_xxyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 1922);

    auto tr_z_xxyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 1923);

    auto tr_z_xxyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 1924);

    auto tr_z_xxyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 1925);

    auto tr_z_xxyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 1926);

    auto tr_z_xxyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 1927);

    auto tr_z_xxyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 1928);

    auto tr_z_xxyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 1929);

    auto tr_z_xxyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 1930);

    auto tr_z_xxyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 1931);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyyzz_xxxxxx, tr_z_xxyyzz_xxxxxy, tr_z_xxyyzz_xxxxxz, tr_z_xxyyzz_xxxxyy, tr_z_xxyyzz_xxxxyz, tr_z_xxyyzz_xxxxzz, tr_z_xxyyzz_xxxyyy, tr_z_xxyyzz_xxxyyz, tr_z_xxyyzz_xxxyzz, tr_z_xxyyzz_xxxzzz, tr_z_xxyyzz_xxyyyy, tr_z_xxyyzz_xxyyyz, tr_z_xxyyzz_xxyyzz, tr_z_xxyyzz_xxyzzz, tr_z_xxyyzz_xxzzzz, tr_z_xxyyzz_xyyyyy, tr_z_xxyyzz_xyyyyz, tr_z_xxyyzz_xyyyzz, tr_z_xxyyzz_xyyzzz, tr_z_xxyyzz_xyzzzz, tr_z_xxyyzz_xzzzzz, tr_z_xxyyzz_yyyyyy, tr_z_xxyyzz_yyyyyz, tr_z_xxyyzz_yyyyzz, tr_z_xxyyzz_yyyzzz, tr_z_xxyyzz_yyzzzz, tr_z_xxyyzz_yzzzzz, tr_z_xxyyzz_zzzzzz, tr_z_xxyzz_xxxxxx, tr_z_xxyzz_xxxxxz, tr_z_xxyzz_xxxxzz, tr_z_xxyzz_xxxzzz, tr_z_xxyzz_xxzzzz, tr_z_xxyzz_xzzzzz, tr_z_xxzz_xxxxxx, tr_z_xxzz_xxxxxz, tr_z_xxzz_xxxxzz, tr_z_xxzz_xxxzzz, tr_z_xxzz_xxzzzz, tr_z_xxzz_xzzzzz, tr_z_xyyzz_xxxxxy, tr_z_xyyzz_xxxxy, tr_z_xyyzz_xxxxyy, tr_z_xyyzz_xxxxyz, tr_z_xyyzz_xxxyy, tr_z_xyyzz_xxxyyy, tr_z_xyyzz_xxxyyz, tr_z_xyyzz_xxxyz, tr_z_xyyzz_xxxyzz, tr_z_xyyzz_xxyyy, tr_z_xyyzz_xxyyyy, tr_z_xyyzz_xxyyyz, tr_z_xyyzz_xxyyz, tr_z_xyyzz_xxyyzz, tr_z_xyyzz_xxyzz, tr_z_xyyzz_xxyzzz, tr_z_xyyzz_xyyyy, tr_z_xyyzz_xyyyyy, tr_z_xyyzz_xyyyyz, tr_z_xyyzz_xyyyz, tr_z_xyyzz_xyyyzz, tr_z_xyyzz_xyyzz, tr_z_xyyzz_xyyzzz, tr_z_xyyzz_xyzzz, tr_z_xyyzz_xyzzzz, tr_z_xyyzz_yyyyy, tr_z_xyyzz_yyyyyy, tr_z_xyyzz_yyyyyz, tr_z_xyyzz_yyyyz, tr_z_xyyzz_yyyyzz, tr_z_xyyzz_yyyzz, tr_z_xyyzz_yyyzzz, tr_z_xyyzz_yyzzz, tr_z_xyyzz_yyzzzz, tr_z_xyyzz_yzzzz, tr_z_xyyzz_yzzzzz, tr_z_xyyzz_zzzzzz, tr_z_yyzz_xxxxxy, tr_z_yyzz_xxxxyy, tr_z_yyzz_xxxxyz, tr_z_yyzz_xxxyyy, tr_z_yyzz_xxxyyz, tr_z_yyzz_xxxyzz, tr_z_yyzz_xxyyyy, tr_z_yyzz_xxyyyz, tr_z_yyzz_xxyyzz, tr_z_yyzz_xxyzzz, tr_z_yyzz_xyyyyy, tr_z_yyzz_xyyyyz, tr_z_yyzz_xyyyzz, tr_z_yyzz_xyyzzz, tr_z_yyzz_xyzzzz, tr_z_yyzz_yyyyyy, tr_z_yyzz_yyyyyz, tr_z_yyzz_yyyyzz, tr_z_yyzz_yyyzzz, tr_z_yyzz_yyzzzz, tr_z_yyzz_yzzzzz, tr_z_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyzz_xxxxxx[i] = tr_z_xxzz_xxxxxx[i] * fe_0 + tr_z_xxyzz_xxxxxx[i] * pa_y[i];

        tr_z_xxyyzz_xxxxxy[i] = tr_z_yyzz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xyyzz_xxxxy[i] * fe_0 + tr_z_xyyzz_xxxxxy[i] * pa_x[i];

        tr_z_xxyyzz_xxxxxz[i] = tr_z_xxzz_xxxxxz[i] * fe_0 + tr_z_xxyzz_xxxxxz[i] * pa_y[i];

        tr_z_xxyyzz_xxxxyy[i] = tr_z_yyzz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xyyzz_xxxyy[i] * fe_0 + tr_z_xyyzz_xxxxyy[i] * pa_x[i];

        tr_z_xxyyzz_xxxxyz[i] = tr_z_yyzz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xyyzz_xxxyz[i] * fe_0 + tr_z_xyyzz_xxxxyz[i] * pa_x[i];

        tr_z_xxyyzz_xxxxzz[i] = tr_z_xxzz_xxxxzz[i] * fe_0 + tr_z_xxyzz_xxxxzz[i] * pa_y[i];

        tr_z_xxyyzz_xxxyyy[i] = tr_z_yyzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xyyzz_xxyyy[i] * fe_0 + tr_z_xyyzz_xxxyyy[i] * pa_x[i];

        tr_z_xxyyzz_xxxyyz[i] = tr_z_yyzz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xyyzz_xxyyz[i] * fe_0 + tr_z_xyyzz_xxxyyz[i] * pa_x[i];

        tr_z_xxyyzz_xxxyzz[i] = tr_z_yyzz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xyyzz_xxyzz[i] * fe_0 + tr_z_xyyzz_xxxyzz[i] * pa_x[i];

        tr_z_xxyyzz_xxxzzz[i] = tr_z_xxzz_xxxzzz[i] * fe_0 + tr_z_xxyzz_xxxzzz[i] * pa_y[i];

        tr_z_xxyyzz_xxyyyy[i] = tr_z_yyzz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xyyzz_xyyyy[i] * fe_0 + tr_z_xyyzz_xxyyyy[i] * pa_x[i];

        tr_z_xxyyzz_xxyyyz[i] = tr_z_yyzz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xyyzz_xyyyz[i] * fe_0 + tr_z_xyyzz_xxyyyz[i] * pa_x[i];

        tr_z_xxyyzz_xxyyzz[i] = tr_z_yyzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xyyzz_xyyzz[i] * fe_0 + tr_z_xyyzz_xxyyzz[i] * pa_x[i];

        tr_z_xxyyzz_xxyzzz[i] = tr_z_yyzz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xyyzz_xyzzz[i] * fe_0 + tr_z_xyyzz_xxyzzz[i] * pa_x[i];

        tr_z_xxyyzz_xxzzzz[i] = tr_z_xxzz_xxzzzz[i] * fe_0 + tr_z_xxyzz_xxzzzz[i] * pa_y[i];

        tr_z_xxyyzz_xyyyyy[i] = tr_z_yyzz_xyyyyy[i] * fe_0 + tr_z_xyyzz_yyyyy[i] * fe_0 + tr_z_xyyzz_xyyyyy[i] * pa_x[i];

        tr_z_xxyyzz_xyyyyz[i] = tr_z_yyzz_xyyyyz[i] * fe_0 + tr_z_xyyzz_yyyyz[i] * fe_0 + tr_z_xyyzz_xyyyyz[i] * pa_x[i];

        tr_z_xxyyzz_xyyyzz[i] = tr_z_yyzz_xyyyzz[i] * fe_0 + tr_z_xyyzz_yyyzz[i] * fe_0 + tr_z_xyyzz_xyyyzz[i] * pa_x[i];

        tr_z_xxyyzz_xyyzzz[i] = tr_z_yyzz_xyyzzz[i] * fe_0 + tr_z_xyyzz_yyzzz[i] * fe_0 + tr_z_xyyzz_xyyzzz[i] * pa_x[i];

        tr_z_xxyyzz_xyzzzz[i] = tr_z_yyzz_xyzzzz[i] * fe_0 + tr_z_xyyzz_yzzzz[i] * fe_0 + tr_z_xyyzz_xyzzzz[i] * pa_x[i];

        tr_z_xxyyzz_xzzzzz[i] = tr_z_xxzz_xzzzzz[i] * fe_0 + tr_z_xxyzz_xzzzzz[i] * pa_y[i];

        tr_z_xxyyzz_yyyyyy[i] = tr_z_yyzz_yyyyyy[i] * fe_0 + tr_z_xyyzz_yyyyyy[i] * pa_x[i];

        tr_z_xxyyzz_yyyyyz[i] = tr_z_yyzz_yyyyyz[i] * fe_0 + tr_z_xyyzz_yyyyyz[i] * pa_x[i];

        tr_z_xxyyzz_yyyyzz[i] = tr_z_yyzz_yyyyzz[i] * fe_0 + tr_z_xyyzz_yyyyzz[i] * pa_x[i];

        tr_z_xxyyzz_yyyzzz[i] = tr_z_yyzz_yyyzzz[i] * fe_0 + tr_z_xyyzz_yyyzzz[i] * pa_x[i];

        tr_z_xxyyzz_yyzzzz[i] = tr_z_yyzz_yyzzzz[i] * fe_0 + tr_z_xyyzz_yyzzzz[i] * pa_x[i];

        tr_z_xxyyzz_yzzzzz[i] = tr_z_yyzz_yzzzzz[i] * fe_0 + tr_z_xyyzz_yzzzzz[i] * pa_x[i];

        tr_z_xxyyzz_zzzzzz[i] = tr_z_yyzz_zzzzzz[i] * fe_0 + tr_z_xyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1932-1960 components of targeted buffer : II

    auto tr_z_xxyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1932);

    auto tr_z_xxyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1933);

    auto tr_z_xxyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1934);

    auto tr_z_xxyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1935);

    auto tr_z_xxyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1936);

    auto tr_z_xxyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1937);

    auto tr_z_xxyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1938);

    auto tr_z_xxyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1939);

    auto tr_z_xxyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1940);

    auto tr_z_xxyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1941);

    auto tr_z_xxyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1942);

    auto tr_z_xxyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1943);

    auto tr_z_xxyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1944);

    auto tr_z_xxyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1945);

    auto tr_z_xxyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1946);

    auto tr_z_xxyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1947);

    auto tr_z_xxyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1948);

    auto tr_z_xxyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1949);

    auto tr_z_xxyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1950);

    auto tr_z_xxyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1951);

    auto tr_z_xxyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1952);

    auto tr_z_xxyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1953);

    auto tr_z_xxyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1954);

    auto tr_z_xxyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1955);

    auto tr_z_xxyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1956);

    auto tr_z_xxyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1957);

    auto tr_z_xxyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1958);

    auto tr_z_xxyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1959);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyzzz_xxxxxx, tr_z_xxyzzz_xxxxxy, tr_z_xxyzzz_xxxxxz, tr_z_xxyzzz_xxxxyy, tr_z_xxyzzz_xxxxyz, tr_z_xxyzzz_xxxxzz, tr_z_xxyzzz_xxxyyy, tr_z_xxyzzz_xxxyyz, tr_z_xxyzzz_xxxyzz, tr_z_xxyzzz_xxxzzz, tr_z_xxyzzz_xxyyyy, tr_z_xxyzzz_xxyyyz, tr_z_xxyzzz_xxyyzz, tr_z_xxyzzz_xxyzzz, tr_z_xxyzzz_xxzzzz, tr_z_xxyzzz_xyyyyy, tr_z_xxyzzz_xyyyyz, tr_z_xxyzzz_xyyyzz, tr_z_xxyzzz_xyyzzz, tr_z_xxyzzz_xyzzzz, tr_z_xxyzzz_xzzzzz, tr_z_xxyzzz_yyyyyy, tr_z_xxyzzz_yyyyyz, tr_z_xxyzzz_yyyyzz, tr_z_xxyzzz_yyyzzz, tr_z_xxyzzz_yyzzzz, tr_z_xxyzzz_yzzzzz, tr_z_xxyzzz_zzzzzz, tr_z_xxzzz_xxxxx, tr_z_xxzzz_xxxxxx, tr_z_xxzzz_xxxxxy, tr_z_xxzzz_xxxxxz, tr_z_xxzzz_xxxxy, tr_z_xxzzz_xxxxyy, tr_z_xxzzz_xxxxyz, tr_z_xxzzz_xxxxz, tr_z_xxzzz_xxxxzz, tr_z_xxzzz_xxxyy, tr_z_xxzzz_xxxyyy, tr_z_xxzzz_xxxyyz, tr_z_xxzzz_xxxyz, tr_z_xxzzz_xxxyzz, tr_z_xxzzz_xxxzz, tr_z_xxzzz_xxxzzz, tr_z_xxzzz_xxyyy, tr_z_xxzzz_xxyyyy, tr_z_xxzzz_xxyyyz, tr_z_xxzzz_xxyyz, tr_z_xxzzz_xxyyzz, tr_z_xxzzz_xxyzz, tr_z_xxzzz_xxyzzz, tr_z_xxzzz_xxzzz, tr_z_xxzzz_xxzzzz, tr_z_xxzzz_xyyyy, tr_z_xxzzz_xyyyyy, tr_z_xxzzz_xyyyyz, tr_z_xxzzz_xyyyz, tr_z_xxzzz_xyyyzz, tr_z_xxzzz_xyyzz, tr_z_xxzzz_xyyzzz, tr_z_xxzzz_xyzzz, tr_z_xxzzz_xyzzzz, tr_z_xxzzz_xzzzz, tr_z_xxzzz_xzzzzz, tr_z_xxzzz_zzzzzz, tr_z_xyzzz_yyyyyy, tr_z_xyzzz_yyyyyz, tr_z_xyzzz_yyyyzz, tr_z_xyzzz_yyyzzz, tr_z_xyzzz_yyzzzz, tr_z_xyzzz_yzzzzz, tr_z_yzzz_yyyyyy, tr_z_yzzz_yyyyyz, tr_z_yzzz_yyyyzz, tr_z_yzzz_yyyzzz, tr_z_yzzz_yyzzzz, tr_z_yzzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzzz_xxxxxx[i] = tr_z_xxzzz_xxxxxx[i] * pa_y[i];

        tr_z_xxyzzz_xxxxxy[i] = tr_z_xxzzz_xxxxx[i] * fe_0 + tr_z_xxzzz_xxxxxy[i] * pa_y[i];

        tr_z_xxyzzz_xxxxxz[i] = tr_z_xxzzz_xxxxxz[i] * pa_y[i];

        tr_z_xxyzzz_xxxxyy[i] = 2.0 * tr_z_xxzzz_xxxxy[i] * fe_0 + tr_z_xxzzz_xxxxyy[i] * pa_y[i];

        tr_z_xxyzzz_xxxxyz[i] = tr_z_xxzzz_xxxxz[i] * fe_0 + tr_z_xxzzz_xxxxyz[i] * pa_y[i];

        tr_z_xxyzzz_xxxxzz[i] = tr_z_xxzzz_xxxxzz[i] * pa_y[i];

        tr_z_xxyzzz_xxxyyy[i] = 3.0 * tr_z_xxzzz_xxxyy[i] * fe_0 + tr_z_xxzzz_xxxyyy[i] * pa_y[i];

        tr_z_xxyzzz_xxxyyz[i] = 2.0 * tr_z_xxzzz_xxxyz[i] * fe_0 + tr_z_xxzzz_xxxyyz[i] * pa_y[i];

        tr_z_xxyzzz_xxxyzz[i] = tr_z_xxzzz_xxxzz[i] * fe_0 + tr_z_xxzzz_xxxyzz[i] * pa_y[i];

        tr_z_xxyzzz_xxxzzz[i] = tr_z_xxzzz_xxxzzz[i] * pa_y[i];

        tr_z_xxyzzz_xxyyyy[i] = 4.0 * tr_z_xxzzz_xxyyy[i] * fe_0 + tr_z_xxzzz_xxyyyy[i] * pa_y[i];

        tr_z_xxyzzz_xxyyyz[i] = 3.0 * tr_z_xxzzz_xxyyz[i] * fe_0 + tr_z_xxzzz_xxyyyz[i] * pa_y[i];

        tr_z_xxyzzz_xxyyzz[i] = 2.0 * tr_z_xxzzz_xxyzz[i] * fe_0 + tr_z_xxzzz_xxyyzz[i] * pa_y[i];

        tr_z_xxyzzz_xxyzzz[i] = tr_z_xxzzz_xxzzz[i] * fe_0 + tr_z_xxzzz_xxyzzz[i] * pa_y[i];

        tr_z_xxyzzz_xxzzzz[i] = tr_z_xxzzz_xxzzzz[i] * pa_y[i];

        tr_z_xxyzzz_xyyyyy[i] = 5.0 * tr_z_xxzzz_xyyyy[i] * fe_0 + tr_z_xxzzz_xyyyyy[i] * pa_y[i];

        tr_z_xxyzzz_xyyyyz[i] = 4.0 * tr_z_xxzzz_xyyyz[i] * fe_0 + tr_z_xxzzz_xyyyyz[i] * pa_y[i];

        tr_z_xxyzzz_xyyyzz[i] = 3.0 * tr_z_xxzzz_xyyzz[i] * fe_0 + tr_z_xxzzz_xyyyzz[i] * pa_y[i];

        tr_z_xxyzzz_xyyzzz[i] = 2.0 * tr_z_xxzzz_xyzzz[i] * fe_0 + tr_z_xxzzz_xyyzzz[i] * pa_y[i];

        tr_z_xxyzzz_xyzzzz[i] = tr_z_xxzzz_xzzzz[i] * fe_0 + tr_z_xxzzz_xyzzzz[i] * pa_y[i];

        tr_z_xxyzzz_xzzzzz[i] = tr_z_xxzzz_xzzzzz[i] * pa_y[i];

        tr_z_xxyzzz_yyyyyy[i] = tr_z_yzzz_yyyyyy[i] * fe_0 + tr_z_xyzzz_yyyyyy[i] * pa_x[i];

        tr_z_xxyzzz_yyyyyz[i] = tr_z_yzzz_yyyyyz[i] * fe_0 + tr_z_xyzzz_yyyyyz[i] * pa_x[i];

        tr_z_xxyzzz_yyyyzz[i] = tr_z_yzzz_yyyyzz[i] * fe_0 + tr_z_xyzzz_yyyyzz[i] * pa_x[i];

        tr_z_xxyzzz_yyyzzz[i] = tr_z_yzzz_yyyzzz[i] * fe_0 + tr_z_xyzzz_yyyzzz[i] * pa_x[i];

        tr_z_xxyzzz_yyzzzz[i] = tr_z_yzzz_yyzzzz[i] * fe_0 + tr_z_xyzzz_yyzzzz[i] * pa_x[i];

        tr_z_xxyzzz_yzzzzz[i] = tr_z_yzzz_yzzzzz[i] * fe_0 + tr_z_xyzzz_yzzzzz[i] * pa_x[i];

        tr_z_xxyzzz_zzzzzz[i] = tr_z_xxzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1960-1988 components of targeted buffer : II

    auto tr_z_xxzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 1960);

    auto tr_z_xxzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 1961);

    auto tr_z_xxzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 1962);

    auto tr_z_xxzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 1963);

    auto tr_z_xxzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 1964);

    auto tr_z_xxzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 1965);

    auto tr_z_xxzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 1966);

    auto tr_z_xxzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 1967);

    auto tr_z_xxzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 1968);

    auto tr_z_xxzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 1969);

    auto tr_z_xxzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 1970);

    auto tr_z_xxzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 1971);

    auto tr_z_xxzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 1972);

    auto tr_z_xxzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 1973);

    auto tr_z_xxzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 1974);

    auto tr_z_xxzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 1975);

    auto tr_z_xxzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 1976);

    auto tr_z_xxzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 1977);

    auto tr_z_xxzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 1978);

    auto tr_z_xxzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 1979);

    auto tr_z_xxzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 1980);

    auto tr_z_xxzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 1981);

    auto tr_z_xxzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 1982);

    auto tr_z_xxzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 1983);

    auto tr_z_xxzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 1984);

    auto tr_z_xxzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 1985);

    auto tr_z_xxzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 1986);

    auto tr_z_xxzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 1987);

    #pragma omp simd aligned(pa_x, tr_z_xxzzzz_xxxxxx, tr_z_xxzzzz_xxxxxy, tr_z_xxzzzz_xxxxxz, tr_z_xxzzzz_xxxxyy, tr_z_xxzzzz_xxxxyz, tr_z_xxzzzz_xxxxzz, tr_z_xxzzzz_xxxyyy, tr_z_xxzzzz_xxxyyz, tr_z_xxzzzz_xxxyzz, tr_z_xxzzzz_xxxzzz, tr_z_xxzzzz_xxyyyy, tr_z_xxzzzz_xxyyyz, tr_z_xxzzzz_xxyyzz, tr_z_xxzzzz_xxyzzz, tr_z_xxzzzz_xxzzzz, tr_z_xxzzzz_xyyyyy, tr_z_xxzzzz_xyyyyz, tr_z_xxzzzz_xyyyzz, tr_z_xxzzzz_xyyzzz, tr_z_xxzzzz_xyzzzz, tr_z_xxzzzz_xzzzzz, tr_z_xxzzzz_yyyyyy, tr_z_xxzzzz_yyyyyz, tr_z_xxzzzz_yyyyzz, tr_z_xxzzzz_yyyzzz, tr_z_xxzzzz_yyzzzz, tr_z_xxzzzz_yzzzzz, tr_z_xxzzzz_zzzzzz, tr_z_xzzzz_xxxxx, tr_z_xzzzz_xxxxxx, tr_z_xzzzz_xxxxxy, tr_z_xzzzz_xxxxxz, tr_z_xzzzz_xxxxy, tr_z_xzzzz_xxxxyy, tr_z_xzzzz_xxxxyz, tr_z_xzzzz_xxxxz, tr_z_xzzzz_xxxxzz, tr_z_xzzzz_xxxyy, tr_z_xzzzz_xxxyyy, tr_z_xzzzz_xxxyyz, tr_z_xzzzz_xxxyz, tr_z_xzzzz_xxxyzz, tr_z_xzzzz_xxxzz, tr_z_xzzzz_xxxzzz, tr_z_xzzzz_xxyyy, tr_z_xzzzz_xxyyyy, tr_z_xzzzz_xxyyyz, tr_z_xzzzz_xxyyz, tr_z_xzzzz_xxyyzz, tr_z_xzzzz_xxyzz, tr_z_xzzzz_xxyzzz, tr_z_xzzzz_xxzzz, tr_z_xzzzz_xxzzzz, tr_z_xzzzz_xyyyy, tr_z_xzzzz_xyyyyy, tr_z_xzzzz_xyyyyz, tr_z_xzzzz_xyyyz, tr_z_xzzzz_xyyyzz, tr_z_xzzzz_xyyzz, tr_z_xzzzz_xyyzzz, tr_z_xzzzz_xyzzz, tr_z_xzzzz_xyzzzz, tr_z_xzzzz_xzzzz, tr_z_xzzzz_xzzzzz, tr_z_xzzzz_yyyyy, tr_z_xzzzz_yyyyyy, tr_z_xzzzz_yyyyyz, tr_z_xzzzz_yyyyz, tr_z_xzzzz_yyyyzz, tr_z_xzzzz_yyyzz, tr_z_xzzzz_yyyzzz, tr_z_xzzzz_yyzzz, tr_z_xzzzz_yyzzzz, tr_z_xzzzz_yzzzz, tr_z_xzzzz_yzzzzz, tr_z_xzzzz_zzzzz, tr_z_xzzzz_zzzzzz, tr_z_zzzz_xxxxxx, tr_z_zzzz_xxxxxy, tr_z_zzzz_xxxxxz, tr_z_zzzz_xxxxyy, tr_z_zzzz_xxxxyz, tr_z_zzzz_xxxxzz, tr_z_zzzz_xxxyyy, tr_z_zzzz_xxxyyz, tr_z_zzzz_xxxyzz, tr_z_zzzz_xxxzzz, tr_z_zzzz_xxyyyy, tr_z_zzzz_xxyyyz, tr_z_zzzz_xxyyzz, tr_z_zzzz_xxyzzz, tr_z_zzzz_xxzzzz, tr_z_zzzz_xyyyyy, tr_z_zzzz_xyyyyz, tr_z_zzzz_xyyyzz, tr_z_zzzz_xyyzzz, tr_z_zzzz_xyzzzz, tr_z_zzzz_xzzzzz, tr_z_zzzz_yyyyyy, tr_z_zzzz_yyyyyz, tr_z_zzzz_yyyyzz, tr_z_zzzz_yyyzzz, tr_z_zzzz_yyzzzz, tr_z_zzzz_yzzzzz, tr_z_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzzz_xxxxxx[i] = tr_z_zzzz_xxxxxx[i] * fe_0 + 6.0 * tr_z_xzzzz_xxxxx[i] * fe_0 + tr_z_xzzzz_xxxxxx[i] * pa_x[i];

        tr_z_xxzzzz_xxxxxy[i] = tr_z_zzzz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xzzzz_xxxxy[i] * fe_0 + tr_z_xzzzz_xxxxxy[i] * pa_x[i];

        tr_z_xxzzzz_xxxxxz[i] = tr_z_zzzz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xzzzz_xxxxz[i] * fe_0 + tr_z_xzzzz_xxxxxz[i] * pa_x[i];

        tr_z_xxzzzz_xxxxyy[i] = tr_z_zzzz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xzzzz_xxxyy[i] * fe_0 + tr_z_xzzzz_xxxxyy[i] * pa_x[i];

        tr_z_xxzzzz_xxxxyz[i] = tr_z_zzzz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xzzzz_xxxyz[i] * fe_0 + tr_z_xzzzz_xxxxyz[i] * pa_x[i];

        tr_z_xxzzzz_xxxxzz[i] = tr_z_zzzz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xzzzz_xxxzz[i] * fe_0 + tr_z_xzzzz_xxxxzz[i] * pa_x[i];

        tr_z_xxzzzz_xxxyyy[i] = tr_z_zzzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xzzzz_xxyyy[i] * fe_0 + tr_z_xzzzz_xxxyyy[i] * pa_x[i];

        tr_z_xxzzzz_xxxyyz[i] = tr_z_zzzz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xzzzz_xxyyz[i] * fe_0 + tr_z_xzzzz_xxxyyz[i] * pa_x[i];

        tr_z_xxzzzz_xxxyzz[i] = tr_z_zzzz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xzzzz_xxyzz[i] * fe_0 + tr_z_xzzzz_xxxyzz[i] * pa_x[i];

        tr_z_xxzzzz_xxxzzz[i] = tr_z_zzzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xzzzz_xxzzz[i] * fe_0 + tr_z_xzzzz_xxxzzz[i] * pa_x[i];

        tr_z_xxzzzz_xxyyyy[i] = tr_z_zzzz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xzzzz_xyyyy[i] * fe_0 + tr_z_xzzzz_xxyyyy[i] * pa_x[i];

        tr_z_xxzzzz_xxyyyz[i] = tr_z_zzzz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xzzzz_xyyyz[i] * fe_0 + tr_z_xzzzz_xxyyyz[i] * pa_x[i];

        tr_z_xxzzzz_xxyyzz[i] = tr_z_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xzzzz_xyyzz[i] * fe_0 + tr_z_xzzzz_xxyyzz[i] * pa_x[i];

        tr_z_xxzzzz_xxyzzz[i] = tr_z_zzzz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xzzzz_xyzzz[i] * fe_0 + tr_z_xzzzz_xxyzzz[i] * pa_x[i];

        tr_z_xxzzzz_xxzzzz[i] = tr_z_zzzz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xzzzz_xzzzz[i] * fe_0 + tr_z_xzzzz_xxzzzz[i] * pa_x[i];

        tr_z_xxzzzz_xyyyyy[i] = tr_z_zzzz_xyyyyy[i] * fe_0 + tr_z_xzzzz_yyyyy[i] * fe_0 + tr_z_xzzzz_xyyyyy[i] * pa_x[i];

        tr_z_xxzzzz_xyyyyz[i] = tr_z_zzzz_xyyyyz[i] * fe_0 + tr_z_xzzzz_yyyyz[i] * fe_0 + tr_z_xzzzz_xyyyyz[i] * pa_x[i];

        tr_z_xxzzzz_xyyyzz[i] = tr_z_zzzz_xyyyzz[i] * fe_0 + tr_z_xzzzz_yyyzz[i] * fe_0 + tr_z_xzzzz_xyyyzz[i] * pa_x[i];

        tr_z_xxzzzz_xyyzzz[i] = tr_z_zzzz_xyyzzz[i] * fe_0 + tr_z_xzzzz_yyzzz[i] * fe_0 + tr_z_xzzzz_xyyzzz[i] * pa_x[i];

        tr_z_xxzzzz_xyzzzz[i] = tr_z_zzzz_xyzzzz[i] * fe_0 + tr_z_xzzzz_yzzzz[i] * fe_0 + tr_z_xzzzz_xyzzzz[i] * pa_x[i];

        tr_z_xxzzzz_xzzzzz[i] = tr_z_zzzz_xzzzzz[i] * fe_0 + tr_z_xzzzz_zzzzz[i] * fe_0 + tr_z_xzzzz_xzzzzz[i] * pa_x[i];

        tr_z_xxzzzz_yyyyyy[i] = tr_z_zzzz_yyyyyy[i] * fe_0 + tr_z_xzzzz_yyyyyy[i] * pa_x[i];

        tr_z_xxzzzz_yyyyyz[i] = tr_z_zzzz_yyyyyz[i] * fe_0 + tr_z_xzzzz_yyyyyz[i] * pa_x[i];

        tr_z_xxzzzz_yyyyzz[i] = tr_z_zzzz_yyyyzz[i] * fe_0 + tr_z_xzzzz_yyyyzz[i] * pa_x[i];

        tr_z_xxzzzz_yyyzzz[i] = tr_z_zzzz_yyyzzz[i] * fe_0 + tr_z_xzzzz_yyyzzz[i] * pa_x[i];

        tr_z_xxzzzz_yyzzzz[i] = tr_z_zzzz_yyzzzz[i] * fe_0 + tr_z_xzzzz_yyzzzz[i] * pa_x[i];

        tr_z_xxzzzz_yzzzzz[i] = tr_z_zzzz_yzzzzz[i] * fe_0 + tr_z_xzzzz_yzzzzz[i] * pa_x[i];

        tr_z_xxzzzz_zzzzzz[i] = tr_z_zzzz_zzzzzz[i] * fe_0 + tr_z_xzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1988-2016 components of targeted buffer : II

    auto tr_z_xyyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 1988);

    auto tr_z_xyyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 1989);

    auto tr_z_xyyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 1990);

    auto tr_z_xyyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 1991);

    auto tr_z_xyyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 1992);

    auto tr_z_xyyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 1993);

    auto tr_z_xyyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 1994);

    auto tr_z_xyyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 1995);

    auto tr_z_xyyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 1996);

    auto tr_z_xyyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 1997);

    auto tr_z_xyyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 1998);

    auto tr_z_xyyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 1999);

    auto tr_z_xyyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 2000);

    auto tr_z_xyyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 2001);

    auto tr_z_xyyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 2002);

    auto tr_z_xyyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 2003);

    auto tr_z_xyyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 2004);

    auto tr_z_xyyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 2005);

    auto tr_z_xyyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 2006);

    auto tr_z_xyyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 2007);

    auto tr_z_xyyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 2008);

    auto tr_z_xyyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 2009);

    auto tr_z_xyyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 2010);

    auto tr_z_xyyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 2011);

    auto tr_z_xyyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 2012);

    auto tr_z_xyyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 2013);

    auto tr_z_xyyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 2014);

    auto tr_z_xyyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 2015);

    #pragma omp simd aligned(pa_x, tr_z_xyyyyy_xxxxxx, tr_z_xyyyyy_xxxxxy, tr_z_xyyyyy_xxxxxz, tr_z_xyyyyy_xxxxyy, tr_z_xyyyyy_xxxxyz, tr_z_xyyyyy_xxxxzz, tr_z_xyyyyy_xxxyyy, tr_z_xyyyyy_xxxyyz, tr_z_xyyyyy_xxxyzz, tr_z_xyyyyy_xxxzzz, tr_z_xyyyyy_xxyyyy, tr_z_xyyyyy_xxyyyz, tr_z_xyyyyy_xxyyzz, tr_z_xyyyyy_xxyzzz, tr_z_xyyyyy_xxzzzz, tr_z_xyyyyy_xyyyyy, tr_z_xyyyyy_xyyyyz, tr_z_xyyyyy_xyyyzz, tr_z_xyyyyy_xyyzzz, tr_z_xyyyyy_xyzzzz, tr_z_xyyyyy_xzzzzz, tr_z_xyyyyy_yyyyyy, tr_z_xyyyyy_yyyyyz, tr_z_xyyyyy_yyyyzz, tr_z_xyyyyy_yyyzzz, tr_z_xyyyyy_yyzzzz, tr_z_xyyyyy_yzzzzz, tr_z_xyyyyy_zzzzzz, tr_z_yyyyy_xxxxx, tr_z_yyyyy_xxxxxx, tr_z_yyyyy_xxxxxy, tr_z_yyyyy_xxxxxz, tr_z_yyyyy_xxxxy, tr_z_yyyyy_xxxxyy, tr_z_yyyyy_xxxxyz, tr_z_yyyyy_xxxxz, tr_z_yyyyy_xxxxzz, tr_z_yyyyy_xxxyy, tr_z_yyyyy_xxxyyy, tr_z_yyyyy_xxxyyz, tr_z_yyyyy_xxxyz, tr_z_yyyyy_xxxyzz, tr_z_yyyyy_xxxzz, tr_z_yyyyy_xxxzzz, tr_z_yyyyy_xxyyy, tr_z_yyyyy_xxyyyy, tr_z_yyyyy_xxyyyz, tr_z_yyyyy_xxyyz, tr_z_yyyyy_xxyyzz, tr_z_yyyyy_xxyzz, tr_z_yyyyy_xxyzzz, tr_z_yyyyy_xxzzz, tr_z_yyyyy_xxzzzz, tr_z_yyyyy_xyyyy, tr_z_yyyyy_xyyyyy, tr_z_yyyyy_xyyyyz, tr_z_yyyyy_xyyyz, tr_z_yyyyy_xyyyzz, tr_z_yyyyy_xyyzz, tr_z_yyyyy_xyyzzz, tr_z_yyyyy_xyzzz, tr_z_yyyyy_xyzzzz, tr_z_yyyyy_xzzzz, tr_z_yyyyy_xzzzzz, tr_z_yyyyy_yyyyy, tr_z_yyyyy_yyyyyy, tr_z_yyyyy_yyyyyz, tr_z_yyyyy_yyyyz, tr_z_yyyyy_yyyyzz, tr_z_yyyyy_yyyzz, tr_z_yyyyy_yyyzzz, tr_z_yyyyy_yyzzz, tr_z_yyyyy_yyzzzz, tr_z_yyyyy_yzzzz, tr_z_yyyyy_yzzzzz, tr_z_yyyyy_zzzzz, tr_z_yyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyy_xxxxxx[i] = 6.0 * tr_z_yyyyy_xxxxx[i] * fe_0 + tr_z_yyyyy_xxxxxx[i] * pa_x[i];

        tr_z_xyyyyy_xxxxxy[i] = 5.0 * tr_z_yyyyy_xxxxy[i] * fe_0 + tr_z_yyyyy_xxxxxy[i] * pa_x[i];

        tr_z_xyyyyy_xxxxxz[i] = 5.0 * tr_z_yyyyy_xxxxz[i] * fe_0 + tr_z_yyyyy_xxxxxz[i] * pa_x[i];

        tr_z_xyyyyy_xxxxyy[i] = 4.0 * tr_z_yyyyy_xxxyy[i] * fe_0 + tr_z_yyyyy_xxxxyy[i] * pa_x[i];

        tr_z_xyyyyy_xxxxyz[i] = 4.0 * tr_z_yyyyy_xxxyz[i] * fe_0 + tr_z_yyyyy_xxxxyz[i] * pa_x[i];

        tr_z_xyyyyy_xxxxzz[i] = 4.0 * tr_z_yyyyy_xxxzz[i] * fe_0 + tr_z_yyyyy_xxxxzz[i] * pa_x[i];

        tr_z_xyyyyy_xxxyyy[i] = 3.0 * tr_z_yyyyy_xxyyy[i] * fe_0 + tr_z_yyyyy_xxxyyy[i] * pa_x[i];

        tr_z_xyyyyy_xxxyyz[i] = 3.0 * tr_z_yyyyy_xxyyz[i] * fe_0 + tr_z_yyyyy_xxxyyz[i] * pa_x[i];

        tr_z_xyyyyy_xxxyzz[i] = 3.0 * tr_z_yyyyy_xxyzz[i] * fe_0 + tr_z_yyyyy_xxxyzz[i] * pa_x[i];

        tr_z_xyyyyy_xxxzzz[i] = 3.0 * tr_z_yyyyy_xxzzz[i] * fe_0 + tr_z_yyyyy_xxxzzz[i] * pa_x[i];

        tr_z_xyyyyy_xxyyyy[i] = 2.0 * tr_z_yyyyy_xyyyy[i] * fe_0 + tr_z_yyyyy_xxyyyy[i] * pa_x[i];

        tr_z_xyyyyy_xxyyyz[i] = 2.0 * tr_z_yyyyy_xyyyz[i] * fe_0 + tr_z_yyyyy_xxyyyz[i] * pa_x[i];

        tr_z_xyyyyy_xxyyzz[i] = 2.0 * tr_z_yyyyy_xyyzz[i] * fe_0 + tr_z_yyyyy_xxyyzz[i] * pa_x[i];

        tr_z_xyyyyy_xxyzzz[i] = 2.0 * tr_z_yyyyy_xyzzz[i] * fe_0 + tr_z_yyyyy_xxyzzz[i] * pa_x[i];

        tr_z_xyyyyy_xxzzzz[i] = 2.0 * tr_z_yyyyy_xzzzz[i] * fe_0 + tr_z_yyyyy_xxzzzz[i] * pa_x[i];

        tr_z_xyyyyy_xyyyyy[i] = tr_z_yyyyy_yyyyy[i] * fe_0 + tr_z_yyyyy_xyyyyy[i] * pa_x[i];

        tr_z_xyyyyy_xyyyyz[i] = tr_z_yyyyy_yyyyz[i] * fe_0 + tr_z_yyyyy_xyyyyz[i] * pa_x[i];

        tr_z_xyyyyy_xyyyzz[i] = tr_z_yyyyy_yyyzz[i] * fe_0 + tr_z_yyyyy_xyyyzz[i] * pa_x[i];

        tr_z_xyyyyy_xyyzzz[i] = tr_z_yyyyy_yyzzz[i] * fe_0 + tr_z_yyyyy_xyyzzz[i] * pa_x[i];

        tr_z_xyyyyy_xyzzzz[i] = tr_z_yyyyy_yzzzz[i] * fe_0 + tr_z_yyyyy_xyzzzz[i] * pa_x[i];

        tr_z_xyyyyy_xzzzzz[i] = tr_z_yyyyy_zzzzz[i] * fe_0 + tr_z_yyyyy_xzzzzz[i] * pa_x[i];

        tr_z_xyyyyy_yyyyyy[i] = tr_z_yyyyy_yyyyyy[i] * pa_x[i];

        tr_z_xyyyyy_yyyyyz[i] = tr_z_yyyyy_yyyyyz[i] * pa_x[i];

        tr_z_xyyyyy_yyyyzz[i] = tr_z_yyyyy_yyyyzz[i] * pa_x[i];

        tr_z_xyyyyy_yyyzzz[i] = tr_z_yyyyy_yyyzzz[i] * pa_x[i];

        tr_z_xyyyyy_yyzzzz[i] = tr_z_yyyyy_yyzzzz[i] * pa_x[i];

        tr_z_xyyyyy_yzzzzz[i] = tr_z_yyyyy_yzzzzz[i] * pa_x[i];

        tr_z_xyyyyy_zzzzzz[i] = tr_z_yyyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 2016-2044 components of targeted buffer : II

    auto tr_z_xyyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 2016);

    auto tr_z_xyyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 2017);

    auto tr_z_xyyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 2018);

    auto tr_z_xyyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 2019);

    auto tr_z_xyyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 2020);

    auto tr_z_xyyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 2021);

    auto tr_z_xyyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 2022);

    auto tr_z_xyyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 2023);

    auto tr_z_xyyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 2024);

    auto tr_z_xyyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 2025);

    auto tr_z_xyyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 2026);

    auto tr_z_xyyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 2027);

    auto tr_z_xyyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 2028);

    auto tr_z_xyyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 2029);

    auto tr_z_xyyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 2030);

    auto tr_z_xyyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 2031);

    auto tr_z_xyyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 2032);

    auto tr_z_xyyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 2033);

    auto tr_z_xyyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 2034);

    auto tr_z_xyyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 2035);

    auto tr_z_xyyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 2036);

    auto tr_z_xyyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 2037);

    auto tr_z_xyyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 2038);

    auto tr_z_xyyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 2039);

    auto tr_z_xyyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 2040);

    auto tr_z_xyyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 2041);

    auto tr_z_xyyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 2042);

    auto tr_z_xyyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 2043);

    #pragma omp simd aligned(pa_x, tr_z_xyyyyz_xxxxxx, tr_z_xyyyyz_xxxxxy, tr_z_xyyyyz_xxxxxz, tr_z_xyyyyz_xxxxyy, tr_z_xyyyyz_xxxxyz, tr_z_xyyyyz_xxxxzz, tr_z_xyyyyz_xxxyyy, tr_z_xyyyyz_xxxyyz, tr_z_xyyyyz_xxxyzz, tr_z_xyyyyz_xxxzzz, tr_z_xyyyyz_xxyyyy, tr_z_xyyyyz_xxyyyz, tr_z_xyyyyz_xxyyzz, tr_z_xyyyyz_xxyzzz, tr_z_xyyyyz_xxzzzz, tr_z_xyyyyz_xyyyyy, tr_z_xyyyyz_xyyyyz, tr_z_xyyyyz_xyyyzz, tr_z_xyyyyz_xyyzzz, tr_z_xyyyyz_xyzzzz, tr_z_xyyyyz_xzzzzz, tr_z_xyyyyz_yyyyyy, tr_z_xyyyyz_yyyyyz, tr_z_xyyyyz_yyyyzz, tr_z_xyyyyz_yyyzzz, tr_z_xyyyyz_yyzzzz, tr_z_xyyyyz_yzzzzz, tr_z_xyyyyz_zzzzzz, tr_z_yyyyz_xxxxx, tr_z_yyyyz_xxxxxx, tr_z_yyyyz_xxxxxy, tr_z_yyyyz_xxxxxz, tr_z_yyyyz_xxxxy, tr_z_yyyyz_xxxxyy, tr_z_yyyyz_xxxxyz, tr_z_yyyyz_xxxxz, tr_z_yyyyz_xxxxzz, tr_z_yyyyz_xxxyy, tr_z_yyyyz_xxxyyy, tr_z_yyyyz_xxxyyz, tr_z_yyyyz_xxxyz, tr_z_yyyyz_xxxyzz, tr_z_yyyyz_xxxzz, tr_z_yyyyz_xxxzzz, tr_z_yyyyz_xxyyy, tr_z_yyyyz_xxyyyy, tr_z_yyyyz_xxyyyz, tr_z_yyyyz_xxyyz, tr_z_yyyyz_xxyyzz, tr_z_yyyyz_xxyzz, tr_z_yyyyz_xxyzzz, tr_z_yyyyz_xxzzz, tr_z_yyyyz_xxzzzz, tr_z_yyyyz_xyyyy, tr_z_yyyyz_xyyyyy, tr_z_yyyyz_xyyyyz, tr_z_yyyyz_xyyyz, tr_z_yyyyz_xyyyzz, tr_z_yyyyz_xyyzz, tr_z_yyyyz_xyyzzz, tr_z_yyyyz_xyzzz, tr_z_yyyyz_xyzzzz, tr_z_yyyyz_xzzzz, tr_z_yyyyz_xzzzzz, tr_z_yyyyz_yyyyy, tr_z_yyyyz_yyyyyy, tr_z_yyyyz_yyyyyz, tr_z_yyyyz_yyyyz, tr_z_yyyyz_yyyyzz, tr_z_yyyyz_yyyzz, tr_z_yyyyz_yyyzzz, tr_z_yyyyz_yyzzz, tr_z_yyyyz_yyzzzz, tr_z_yyyyz_yzzzz, tr_z_yyyyz_yzzzzz, tr_z_yyyyz_zzzzz, tr_z_yyyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyz_xxxxxx[i] = 6.0 * tr_z_yyyyz_xxxxx[i] * fe_0 + tr_z_yyyyz_xxxxxx[i] * pa_x[i];

        tr_z_xyyyyz_xxxxxy[i] = 5.0 * tr_z_yyyyz_xxxxy[i] * fe_0 + tr_z_yyyyz_xxxxxy[i] * pa_x[i];

        tr_z_xyyyyz_xxxxxz[i] = 5.0 * tr_z_yyyyz_xxxxz[i] * fe_0 + tr_z_yyyyz_xxxxxz[i] * pa_x[i];

        tr_z_xyyyyz_xxxxyy[i] = 4.0 * tr_z_yyyyz_xxxyy[i] * fe_0 + tr_z_yyyyz_xxxxyy[i] * pa_x[i];

        tr_z_xyyyyz_xxxxyz[i] = 4.0 * tr_z_yyyyz_xxxyz[i] * fe_0 + tr_z_yyyyz_xxxxyz[i] * pa_x[i];

        tr_z_xyyyyz_xxxxzz[i] = 4.0 * tr_z_yyyyz_xxxzz[i] * fe_0 + tr_z_yyyyz_xxxxzz[i] * pa_x[i];

        tr_z_xyyyyz_xxxyyy[i] = 3.0 * tr_z_yyyyz_xxyyy[i] * fe_0 + tr_z_yyyyz_xxxyyy[i] * pa_x[i];

        tr_z_xyyyyz_xxxyyz[i] = 3.0 * tr_z_yyyyz_xxyyz[i] * fe_0 + tr_z_yyyyz_xxxyyz[i] * pa_x[i];

        tr_z_xyyyyz_xxxyzz[i] = 3.0 * tr_z_yyyyz_xxyzz[i] * fe_0 + tr_z_yyyyz_xxxyzz[i] * pa_x[i];

        tr_z_xyyyyz_xxxzzz[i] = 3.0 * tr_z_yyyyz_xxzzz[i] * fe_0 + tr_z_yyyyz_xxxzzz[i] * pa_x[i];

        tr_z_xyyyyz_xxyyyy[i] = 2.0 * tr_z_yyyyz_xyyyy[i] * fe_0 + tr_z_yyyyz_xxyyyy[i] * pa_x[i];

        tr_z_xyyyyz_xxyyyz[i] = 2.0 * tr_z_yyyyz_xyyyz[i] * fe_0 + tr_z_yyyyz_xxyyyz[i] * pa_x[i];

        tr_z_xyyyyz_xxyyzz[i] = 2.0 * tr_z_yyyyz_xyyzz[i] * fe_0 + tr_z_yyyyz_xxyyzz[i] * pa_x[i];

        tr_z_xyyyyz_xxyzzz[i] = 2.0 * tr_z_yyyyz_xyzzz[i] * fe_0 + tr_z_yyyyz_xxyzzz[i] * pa_x[i];

        tr_z_xyyyyz_xxzzzz[i] = 2.0 * tr_z_yyyyz_xzzzz[i] * fe_0 + tr_z_yyyyz_xxzzzz[i] * pa_x[i];

        tr_z_xyyyyz_xyyyyy[i] = tr_z_yyyyz_yyyyy[i] * fe_0 + tr_z_yyyyz_xyyyyy[i] * pa_x[i];

        tr_z_xyyyyz_xyyyyz[i] = tr_z_yyyyz_yyyyz[i] * fe_0 + tr_z_yyyyz_xyyyyz[i] * pa_x[i];

        tr_z_xyyyyz_xyyyzz[i] = tr_z_yyyyz_yyyzz[i] * fe_0 + tr_z_yyyyz_xyyyzz[i] * pa_x[i];

        tr_z_xyyyyz_xyyzzz[i] = tr_z_yyyyz_yyzzz[i] * fe_0 + tr_z_yyyyz_xyyzzz[i] * pa_x[i];

        tr_z_xyyyyz_xyzzzz[i] = tr_z_yyyyz_yzzzz[i] * fe_0 + tr_z_yyyyz_xyzzzz[i] * pa_x[i];

        tr_z_xyyyyz_xzzzzz[i] = tr_z_yyyyz_zzzzz[i] * fe_0 + tr_z_yyyyz_xzzzzz[i] * pa_x[i];

        tr_z_xyyyyz_yyyyyy[i] = tr_z_yyyyz_yyyyyy[i] * pa_x[i];

        tr_z_xyyyyz_yyyyyz[i] = tr_z_yyyyz_yyyyyz[i] * pa_x[i];

        tr_z_xyyyyz_yyyyzz[i] = tr_z_yyyyz_yyyyzz[i] * pa_x[i];

        tr_z_xyyyyz_yyyzzz[i] = tr_z_yyyyz_yyyzzz[i] * pa_x[i];

        tr_z_xyyyyz_yyzzzz[i] = tr_z_yyyyz_yyzzzz[i] * pa_x[i];

        tr_z_xyyyyz_yzzzzz[i] = tr_z_yyyyz_yzzzzz[i] * pa_x[i];

        tr_z_xyyyyz_zzzzzz[i] = tr_z_yyyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 2044-2072 components of targeted buffer : II

    auto tr_z_xyyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 2044);

    auto tr_z_xyyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 2045);

    auto tr_z_xyyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 2046);

    auto tr_z_xyyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 2047);

    auto tr_z_xyyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 2048);

    auto tr_z_xyyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 2049);

    auto tr_z_xyyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 2050);

    auto tr_z_xyyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 2051);

    auto tr_z_xyyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 2052);

    auto tr_z_xyyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 2053);

    auto tr_z_xyyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 2054);

    auto tr_z_xyyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 2055);

    auto tr_z_xyyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 2056);

    auto tr_z_xyyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 2057);

    auto tr_z_xyyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 2058);

    auto tr_z_xyyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 2059);

    auto tr_z_xyyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 2060);

    auto tr_z_xyyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 2061);

    auto tr_z_xyyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 2062);

    auto tr_z_xyyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 2063);

    auto tr_z_xyyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 2064);

    auto tr_z_xyyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 2065);

    auto tr_z_xyyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 2066);

    auto tr_z_xyyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 2067);

    auto tr_z_xyyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 2068);

    auto tr_z_xyyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 2069);

    auto tr_z_xyyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 2070);

    auto tr_z_xyyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 2071);

    #pragma omp simd aligned(pa_x, tr_z_xyyyzz_xxxxxx, tr_z_xyyyzz_xxxxxy, tr_z_xyyyzz_xxxxxz, tr_z_xyyyzz_xxxxyy, tr_z_xyyyzz_xxxxyz, tr_z_xyyyzz_xxxxzz, tr_z_xyyyzz_xxxyyy, tr_z_xyyyzz_xxxyyz, tr_z_xyyyzz_xxxyzz, tr_z_xyyyzz_xxxzzz, tr_z_xyyyzz_xxyyyy, tr_z_xyyyzz_xxyyyz, tr_z_xyyyzz_xxyyzz, tr_z_xyyyzz_xxyzzz, tr_z_xyyyzz_xxzzzz, tr_z_xyyyzz_xyyyyy, tr_z_xyyyzz_xyyyyz, tr_z_xyyyzz_xyyyzz, tr_z_xyyyzz_xyyzzz, tr_z_xyyyzz_xyzzzz, tr_z_xyyyzz_xzzzzz, tr_z_xyyyzz_yyyyyy, tr_z_xyyyzz_yyyyyz, tr_z_xyyyzz_yyyyzz, tr_z_xyyyzz_yyyzzz, tr_z_xyyyzz_yyzzzz, tr_z_xyyyzz_yzzzzz, tr_z_xyyyzz_zzzzzz, tr_z_yyyzz_xxxxx, tr_z_yyyzz_xxxxxx, tr_z_yyyzz_xxxxxy, tr_z_yyyzz_xxxxxz, tr_z_yyyzz_xxxxy, tr_z_yyyzz_xxxxyy, tr_z_yyyzz_xxxxyz, tr_z_yyyzz_xxxxz, tr_z_yyyzz_xxxxzz, tr_z_yyyzz_xxxyy, tr_z_yyyzz_xxxyyy, tr_z_yyyzz_xxxyyz, tr_z_yyyzz_xxxyz, tr_z_yyyzz_xxxyzz, tr_z_yyyzz_xxxzz, tr_z_yyyzz_xxxzzz, tr_z_yyyzz_xxyyy, tr_z_yyyzz_xxyyyy, tr_z_yyyzz_xxyyyz, tr_z_yyyzz_xxyyz, tr_z_yyyzz_xxyyzz, tr_z_yyyzz_xxyzz, tr_z_yyyzz_xxyzzz, tr_z_yyyzz_xxzzz, tr_z_yyyzz_xxzzzz, tr_z_yyyzz_xyyyy, tr_z_yyyzz_xyyyyy, tr_z_yyyzz_xyyyyz, tr_z_yyyzz_xyyyz, tr_z_yyyzz_xyyyzz, tr_z_yyyzz_xyyzz, tr_z_yyyzz_xyyzzz, tr_z_yyyzz_xyzzz, tr_z_yyyzz_xyzzzz, tr_z_yyyzz_xzzzz, tr_z_yyyzz_xzzzzz, tr_z_yyyzz_yyyyy, tr_z_yyyzz_yyyyyy, tr_z_yyyzz_yyyyyz, tr_z_yyyzz_yyyyz, tr_z_yyyzz_yyyyzz, tr_z_yyyzz_yyyzz, tr_z_yyyzz_yyyzzz, tr_z_yyyzz_yyzzz, tr_z_yyyzz_yyzzzz, tr_z_yyyzz_yzzzz, tr_z_yyyzz_yzzzzz, tr_z_yyyzz_zzzzz, tr_z_yyyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyzz_xxxxxx[i] = 6.0 * tr_z_yyyzz_xxxxx[i] * fe_0 + tr_z_yyyzz_xxxxxx[i] * pa_x[i];

        tr_z_xyyyzz_xxxxxy[i] = 5.0 * tr_z_yyyzz_xxxxy[i] * fe_0 + tr_z_yyyzz_xxxxxy[i] * pa_x[i];

        tr_z_xyyyzz_xxxxxz[i] = 5.0 * tr_z_yyyzz_xxxxz[i] * fe_0 + tr_z_yyyzz_xxxxxz[i] * pa_x[i];

        tr_z_xyyyzz_xxxxyy[i] = 4.0 * tr_z_yyyzz_xxxyy[i] * fe_0 + tr_z_yyyzz_xxxxyy[i] * pa_x[i];

        tr_z_xyyyzz_xxxxyz[i] = 4.0 * tr_z_yyyzz_xxxyz[i] * fe_0 + tr_z_yyyzz_xxxxyz[i] * pa_x[i];

        tr_z_xyyyzz_xxxxzz[i] = 4.0 * tr_z_yyyzz_xxxzz[i] * fe_0 + tr_z_yyyzz_xxxxzz[i] * pa_x[i];

        tr_z_xyyyzz_xxxyyy[i] = 3.0 * tr_z_yyyzz_xxyyy[i] * fe_0 + tr_z_yyyzz_xxxyyy[i] * pa_x[i];

        tr_z_xyyyzz_xxxyyz[i] = 3.0 * tr_z_yyyzz_xxyyz[i] * fe_0 + tr_z_yyyzz_xxxyyz[i] * pa_x[i];

        tr_z_xyyyzz_xxxyzz[i] = 3.0 * tr_z_yyyzz_xxyzz[i] * fe_0 + tr_z_yyyzz_xxxyzz[i] * pa_x[i];

        tr_z_xyyyzz_xxxzzz[i] = 3.0 * tr_z_yyyzz_xxzzz[i] * fe_0 + tr_z_yyyzz_xxxzzz[i] * pa_x[i];

        tr_z_xyyyzz_xxyyyy[i] = 2.0 * tr_z_yyyzz_xyyyy[i] * fe_0 + tr_z_yyyzz_xxyyyy[i] * pa_x[i];

        tr_z_xyyyzz_xxyyyz[i] = 2.0 * tr_z_yyyzz_xyyyz[i] * fe_0 + tr_z_yyyzz_xxyyyz[i] * pa_x[i];

        tr_z_xyyyzz_xxyyzz[i] = 2.0 * tr_z_yyyzz_xyyzz[i] * fe_0 + tr_z_yyyzz_xxyyzz[i] * pa_x[i];

        tr_z_xyyyzz_xxyzzz[i] = 2.0 * tr_z_yyyzz_xyzzz[i] * fe_0 + tr_z_yyyzz_xxyzzz[i] * pa_x[i];

        tr_z_xyyyzz_xxzzzz[i] = 2.0 * tr_z_yyyzz_xzzzz[i] * fe_0 + tr_z_yyyzz_xxzzzz[i] * pa_x[i];

        tr_z_xyyyzz_xyyyyy[i] = tr_z_yyyzz_yyyyy[i] * fe_0 + tr_z_yyyzz_xyyyyy[i] * pa_x[i];

        tr_z_xyyyzz_xyyyyz[i] = tr_z_yyyzz_yyyyz[i] * fe_0 + tr_z_yyyzz_xyyyyz[i] * pa_x[i];

        tr_z_xyyyzz_xyyyzz[i] = tr_z_yyyzz_yyyzz[i] * fe_0 + tr_z_yyyzz_xyyyzz[i] * pa_x[i];

        tr_z_xyyyzz_xyyzzz[i] = tr_z_yyyzz_yyzzz[i] * fe_0 + tr_z_yyyzz_xyyzzz[i] * pa_x[i];

        tr_z_xyyyzz_xyzzzz[i] = tr_z_yyyzz_yzzzz[i] * fe_0 + tr_z_yyyzz_xyzzzz[i] * pa_x[i];

        tr_z_xyyyzz_xzzzzz[i] = tr_z_yyyzz_zzzzz[i] * fe_0 + tr_z_yyyzz_xzzzzz[i] * pa_x[i];

        tr_z_xyyyzz_yyyyyy[i] = tr_z_yyyzz_yyyyyy[i] * pa_x[i];

        tr_z_xyyyzz_yyyyyz[i] = tr_z_yyyzz_yyyyyz[i] * pa_x[i];

        tr_z_xyyyzz_yyyyzz[i] = tr_z_yyyzz_yyyyzz[i] * pa_x[i];

        tr_z_xyyyzz_yyyzzz[i] = tr_z_yyyzz_yyyzzz[i] * pa_x[i];

        tr_z_xyyyzz_yyzzzz[i] = tr_z_yyyzz_yyzzzz[i] * pa_x[i];

        tr_z_xyyyzz_yzzzzz[i] = tr_z_yyyzz_yzzzzz[i] * pa_x[i];

        tr_z_xyyyzz_zzzzzz[i] = tr_z_yyyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 2072-2100 components of targeted buffer : II

    auto tr_z_xyyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2072);

    auto tr_z_xyyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2073);

    auto tr_z_xyyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2074);

    auto tr_z_xyyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2075);

    auto tr_z_xyyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2076);

    auto tr_z_xyyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2077);

    auto tr_z_xyyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2078);

    auto tr_z_xyyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2079);

    auto tr_z_xyyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2080);

    auto tr_z_xyyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2081);

    auto tr_z_xyyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2082);

    auto tr_z_xyyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2083);

    auto tr_z_xyyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2084);

    auto tr_z_xyyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2085);

    auto tr_z_xyyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2086);

    auto tr_z_xyyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2087);

    auto tr_z_xyyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2088);

    auto tr_z_xyyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2089);

    auto tr_z_xyyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2090);

    auto tr_z_xyyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2091);

    auto tr_z_xyyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2092);

    auto tr_z_xyyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2093);

    auto tr_z_xyyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2094);

    auto tr_z_xyyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2095);

    auto tr_z_xyyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2096);

    auto tr_z_xyyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2097);

    auto tr_z_xyyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2098);

    auto tr_z_xyyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2099);

    #pragma omp simd aligned(pa_x, tr_z_xyyzzz_xxxxxx, tr_z_xyyzzz_xxxxxy, tr_z_xyyzzz_xxxxxz, tr_z_xyyzzz_xxxxyy, tr_z_xyyzzz_xxxxyz, tr_z_xyyzzz_xxxxzz, tr_z_xyyzzz_xxxyyy, tr_z_xyyzzz_xxxyyz, tr_z_xyyzzz_xxxyzz, tr_z_xyyzzz_xxxzzz, tr_z_xyyzzz_xxyyyy, tr_z_xyyzzz_xxyyyz, tr_z_xyyzzz_xxyyzz, tr_z_xyyzzz_xxyzzz, tr_z_xyyzzz_xxzzzz, tr_z_xyyzzz_xyyyyy, tr_z_xyyzzz_xyyyyz, tr_z_xyyzzz_xyyyzz, tr_z_xyyzzz_xyyzzz, tr_z_xyyzzz_xyzzzz, tr_z_xyyzzz_xzzzzz, tr_z_xyyzzz_yyyyyy, tr_z_xyyzzz_yyyyyz, tr_z_xyyzzz_yyyyzz, tr_z_xyyzzz_yyyzzz, tr_z_xyyzzz_yyzzzz, tr_z_xyyzzz_yzzzzz, tr_z_xyyzzz_zzzzzz, tr_z_yyzzz_xxxxx, tr_z_yyzzz_xxxxxx, tr_z_yyzzz_xxxxxy, tr_z_yyzzz_xxxxxz, tr_z_yyzzz_xxxxy, tr_z_yyzzz_xxxxyy, tr_z_yyzzz_xxxxyz, tr_z_yyzzz_xxxxz, tr_z_yyzzz_xxxxzz, tr_z_yyzzz_xxxyy, tr_z_yyzzz_xxxyyy, tr_z_yyzzz_xxxyyz, tr_z_yyzzz_xxxyz, tr_z_yyzzz_xxxyzz, tr_z_yyzzz_xxxzz, tr_z_yyzzz_xxxzzz, tr_z_yyzzz_xxyyy, tr_z_yyzzz_xxyyyy, tr_z_yyzzz_xxyyyz, tr_z_yyzzz_xxyyz, tr_z_yyzzz_xxyyzz, tr_z_yyzzz_xxyzz, tr_z_yyzzz_xxyzzz, tr_z_yyzzz_xxzzz, tr_z_yyzzz_xxzzzz, tr_z_yyzzz_xyyyy, tr_z_yyzzz_xyyyyy, tr_z_yyzzz_xyyyyz, tr_z_yyzzz_xyyyz, tr_z_yyzzz_xyyyzz, tr_z_yyzzz_xyyzz, tr_z_yyzzz_xyyzzz, tr_z_yyzzz_xyzzz, tr_z_yyzzz_xyzzzz, tr_z_yyzzz_xzzzz, tr_z_yyzzz_xzzzzz, tr_z_yyzzz_yyyyy, tr_z_yyzzz_yyyyyy, tr_z_yyzzz_yyyyyz, tr_z_yyzzz_yyyyz, tr_z_yyzzz_yyyyzz, tr_z_yyzzz_yyyzz, tr_z_yyzzz_yyyzzz, tr_z_yyzzz_yyzzz, tr_z_yyzzz_yyzzzz, tr_z_yyzzz_yzzzz, tr_z_yyzzz_yzzzzz, tr_z_yyzzz_zzzzz, tr_z_yyzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzzz_xxxxxx[i] = 6.0 * tr_z_yyzzz_xxxxx[i] * fe_0 + tr_z_yyzzz_xxxxxx[i] * pa_x[i];

        tr_z_xyyzzz_xxxxxy[i] = 5.0 * tr_z_yyzzz_xxxxy[i] * fe_0 + tr_z_yyzzz_xxxxxy[i] * pa_x[i];

        tr_z_xyyzzz_xxxxxz[i] = 5.0 * tr_z_yyzzz_xxxxz[i] * fe_0 + tr_z_yyzzz_xxxxxz[i] * pa_x[i];

        tr_z_xyyzzz_xxxxyy[i] = 4.0 * tr_z_yyzzz_xxxyy[i] * fe_0 + tr_z_yyzzz_xxxxyy[i] * pa_x[i];

        tr_z_xyyzzz_xxxxyz[i] = 4.0 * tr_z_yyzzz_xxxyz[i] * fe_0 + tr_z_yyzzz_xxxxyz[i] * pa_x[i];

        tr_z_xyyzzz_xxxxzz[i] = 4.0 * tr_z_yyzzz_xxxzz[i] * fe_0 + tr_z_yyzzz_xxxxzz[i] * pa_x[i];

        tr_z_xyyzzz_xxxyyy[i] = 3.0 * tr_z_yyzzz_xxyyy[i] * fe_0 + tr_z_yyzzz_xxxyyy[i] * pa_x[i];

        tr_z_xyyzzz_xxxyyz[i] = 3.0 * tr_z_yyzzz_xxyyz[i] * fe_0 + tr_z_yyzzz_xxxyyz[i] * pa_x[i];

        tr_z_xyyzzz_xxxyzz[i] = 3.0 * tr_z_yyzzz_xxyzz[i] * fe_0 + tr_z_yyzzz_xxxyzz[i] * pa_x[i];

        tr_z_xyyzzz_xxxzzz[i] = 3.0 * tr_z_yyzzz_xxzzz[i] * fe_0 + tr_z_yyzzz_xxxzzz[i] * pa_x[i];

        tr_z_xyyzzz_xxyyyy[i] = 2.0 * tr_z_yyzzz_xyyyy[i] * fe_0 + tr_z_yyzzz_xxyyyy[i] * pa_x[i];

        tr_z_xyyzzz_xxyyyz[i] = 2.0 * tr_z_yyzzz_xyyyz[i] * fe_0 + tr_z_yyzzz_xxyyyz[i] * pa_x[i];

        tr_z_xyyzzz_xxyyzz[i] = 2.0 * tr_z_yyzzz_xyyzz[i] * fe_0 + tr_z_yyzzz_xxyyzz[i] * pa_x[i];

        tr_z_xyyzzz_xxyzzz[i] = 2.0 * tr_z_yyzzz_xyzzz[i] * fe_0 + tr_z_yyzzz_xxyzzz[i] * pa_x[i];

        tr_z_xyyzzz_xxzzzz[i] = 2.0 * tr_z_yyzzz_xzzzz[i] * fe_0 + tr_z_yyzzz_xxzzzz[i] * pa_x[i];

        tr_z_xyyzzz_xyyyyy[i] = tr_z_yyzzz_yyyyy[i] * fe_0 + tr_z_yyzzz_xyyyyy[i] * pa_x[i];

        tr_z_xyyzzz_xyyyyz[i] = tr_z_yyzzz_yyyyz[i] * fe_0 + tr_z_yyzzz_xyyyyz[i] * pa_x[i];

        tr_z_xyyzzz_xyyyzz[i] = tr_z_yyzzz_yyyzz[i] * fe_0 + tr_z_yyzzz_xyyyzz[i] * pa_x[i];

        tr_z_xyyzzz_xyyzzz[i] = tr_z_yyzzz_yyzzz[i] * fe_0 + tr_z_yyzzz_xyyzzz[i] * pa_x[i];

        tr_z_xyyzzz_xyzzzz[i] = tr_z_yyzzz_yzzzz[i] * fe_0 + tr_z_yyzzz_xyzzzz[i] * pa_x[i];

        tr_z_xyyzzz_xzzzzz[i] = tr_z_yyzzz_zzzzz[i] * fe_0 + tr_z_yyzzz_xzzzzz[i] * pa_x[i];

        tr_z_xyyzzz_yyyyyy[i] = tr_z_yyzzz_yyyyyy[i] * pa_x[i];

        tr_z_xyyzzz_yyyyyz[i] = tr_z_yyzzz_yyyyyz[i] * pa_x[i];

        tr_z_xyyzzz_yyyyzz[i] = tr_z_yyzzz_yyyyzz[i] * pa_x[i];

        tr_z_xyyzzz_yyyzzz[i] = tr_z_yyzzz_yyyzzz[i] * pa_x[i];

        tr_z_xyyzzz_yyzzzz[i] = tr_z_yyzzz_yyzzzz[i] * pa_x[i];

        tr_z_xyyzzz_yzzzzz[i] = tr_z_yyzzz_yzzzzz[i] * pa_x[i];

        tr_z_xyyzzz_zzzzzz[i] = tr_z_yyzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 2100-2128 components of targeted buffer : II

    auto tr_z_xyzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2100);

    auto tr_z_xyzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2101);

    auto tr_z_xyzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2102);

    auto tr_z_xyzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2103);

    auto tr_z_xyzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2104);

    auto tr_z_xyzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2105);

    auto tr_z_xyzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2106);

    auto tr_z_xyzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2107);

    auto tr_z_xyzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2108);

    auto tr_z_xyzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2109);

    auto tr_z_xyzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2110);

    auto tr_z_xyzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2111);

    auto tr_z_xyzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2112);

    auto tr_z_xyzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2113);

    auto tr_z_xyzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2114);

    auto tr_z_xyzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2115);

    auto tr_z_xyzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2116);

    auto tr_z_xyzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2117);

    auto tr_z_xyzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2118);

    auto tr_z_xyzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2119);

    auto tr_z_xyzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2120);

    auto tr_z_xyzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2121);

    auto tr_z_xyzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2122);

    auto tr_z_xyzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2123);

    auto tr_z_xyzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2124);

    auto tr_z_xyzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2125);

    auto tr_z_xyzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2126);

    auto tr_z_xyzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2127);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyzzzz_xxxxxx, tr_z_xyzzzz_xxxxxy, tr_z_xyzzzz_xxxxxz, tr_z_xyzzzz_xxxxyy, tr_z_xyzzzz_xxxxyz, tr_z_xyzzzz_xxxxzz, tr_z_xyzzzz_xxxyyy, tr_z_xyzzzz_xxxyyz, tr_z_xyzzzz_xxxyzz, tr_z_xyzzzz_xxxzzz, tr_z_xyzzzz_xxyyyy, tr_z_xyzzzz_xxyyyz, tr_z_xyzzzz_xxyyzz, tr_z_xyzzzz_xxyzzz, tr_z_xyzzzz_xxzzzz, tr_z_xyzzzz_xyyyyy, tr_z_xyzzzz_xyyyyz, tr_z_xyzzzz_xyyyzz, tr_z_xyzzzz_xyyzzz, tr_z_xyzzzz_xyzzzz, tr_z_xyzzzz_xzzzzz, tr_z_xyzzzz_yyyyyy, tr_z_xyzzzz_yyyyyz, tr_z_xyzzzz_yyyyzz, tr_z_xyzzzz_yyyzzz, tr_z_xyzzzz_yyzzzz, tr_z_xyzzzz_yzzzzz, tr_z_xyzzzz_zzzzzz, tr_z_xzzzz_xxxxxx, tr_z_xzzzz_xxxxxz, tr_z_xzzzz_xxxxzz, tr_z_xzzzz_xxxzzz, tr_z_xzzzz_xxzzzz, tr_z_xzzzz_xzzzzz, tr_z_yzzzz_xxxxxy, tr_z_yzzzz_xxxxy, tr_z_yzzzz_xxxxyy, tr_z_yzzzz_xxxxyz, tr_z_yzzzz_xxxyy, tr_z_yzzzz_xxxyyy, tr_z_yzzzz_xxxyyz, tr_z_yzzzz_xxxyz, tr_z_yzzzz_xxxyzz, tr_z_yzzzz_xxyyy, tr_z_yzzzz_xxyyyy, tr_z_yzzzz_xxyyyz, tr_z_yzzzz_xxyyz, tr_z_yzzzz_xxyyzz, tr_z_yzzzz_xxyzz, tr_z_yzzzz_xxyzzz, tr_z_yzzzz_xyyyy, tr_z_yzzzz_xyyyyy, tr_z_yzzzz_xyyyyz, tr_z_yzzzz_xyyyz, tr_z_yzzzz_xyyyzz, tr_z_yzzzz_xyyzz, tr_z_yzzzz_xyyzzz, tr_z_yzzzz_xyzzz, tr_z_yzzzz_xyzzzz, tr_z_yzzzz_yyyyy, tr_z_yzzzz_yyyyyy, tr_z_yzzzz_yyyyyz, tr_z_yzzzz_yyyyz, tr_z_yzzzz_yyyyzz, tr_z_yzzzz_yyyzz, tr_z_yzzzz_yyyzzz, tr_z_yzzzz_yyzzz, tr_z_yzzzz_yyzzzz, tr_z_yzzzz_yzzzz, tr_z_yzzzz_yzzzzz, tr_z_yzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzzz_xxxxxx[i] = tr_z_xzzzz_xxxxxx[i] * pa_y[i];

        tr_z_xyzzzz_xxxxxy[i] = 5.0 * tr_z_yzzzz_xxxxy[i] * fe_0 + tr_z_yzzzz_xxxxxy[i] * pa_x[i];

        tr_z_xyzzzz_xxxxxz[i] = tr_z_xzzzz_xxxxxz[i] * pa_y[i];

        tr_z_xyzzzz_xxxxyy[i] = 4.0 * tr_z_yzzzz_xxxyy[i] * fe_0 + tr_z_yzzzz_xxxxyy[i] * pa_x[i];

        tr_z_xyzzzz_xxxxyz[i] = 4.0 * tr_z_yzzzz_xxxyz[i] * fe_0 + tr_z_yzzzz_xxxxyz[i] * pa_x[i];

        tr_z_xyzzzz_xxxxzz[i] = tr_z_xzzzz_xxxxzz[i] * pa_y[i];

        tr_z_xyzzzz_xxxyyy[i] = 3.0 * tr_z_yzzzz_xxyyy[i] * fe_0 + tr_z_yzzzz_xxxyyy[i] * pa_x[i];

        tr_z_xyzzzz_xxxyyz[i] = 3.0 * tr_z_yzzzz_xxyyz[i] * fe_0 + tr_z_yzzzz_xxxyyz[i] * pa_x[i];

        tr_z_xyzzzz_xxxyzz[i] = 3.0 * tr_z_yzzzz_xxyzz[i] * fe_0 + tr_z_yzzzz_xxxyzz[i] * pa_x[i];

        tr_z_xyzzzz_xxxzzz[i] = tr_z_xzzzz_xxxzzz[i] * pa_y[i];

        tr_z_xyzzzz_xxyyyy[i] = 2.0 * tr_z_yzzzz_xyyyy[i] * fe_0 + tr_z_yzzzz_xxyyyy[i] * pa_x[i];

        tr_z_xyzzzz_xxyyyz[i] = 2.0 * tr_z_yzzzz_xyyyz[i] * fe_0 + tr_z_yzzzz_xxyyyz[i] * pa_x[i];

        tr_z_xyzzzz_xxyyzz[i] = 2.0 * tr_z_yzzzz_xyyzz[i] * fe_0 + tr_z_yzzzz_xxyyzz[i] * pa_x[i];

        tr_z_xyzzzz_xxyzzz[i] = 2.0 * tr_z_yzzzz_xyzzz[i] * fe_0 + tr_z_yzzzz_xxyzzz[i] * pa_x[i];

        tr_z_xyzzzz_xxzzzz[i] = tr_z_xzzzz_xxzzzz[i] * pa_y[i];

        tr_z_xyzzzz_xyyyyy[i] = tr_z_yzzzz_yyyyy[i] * fe_0 + tr_z_yzzzz_xyyyyy[i] * pa_x[i];

        tr_z_xyzzzz_xyyyyz[i] = tr_z_yzzzz_yyyyz[i] * fe_0 + tr_z_yzzzz_xyyyyz[i] * pa_x[i];

        tr_z_xyzzzz_xyyyzz[i] = tr_z_yzzzz_yyyzz[i] * fe_0 + tr_z_yzzzz_xyyyzz[i] * pa_x[i];

        tr_z_xyzzzz_xyyzzz[i] = tr_z_yzzzz_yyzzz[i] * fe_0 + tr_z_yzzzz_xyyzzz[i] * pa_x[i];

        tr_z_xyzzzz_xyzzzz[i] = tr_z_yzzzz_yzzzz[i] * fe_0 + tr_z_yzzzz_xyzzzz[i] * pa_x[i];

        tr_z_xyzzzz_xzzzzz[i] = tr_z_xzzzz_xzzzzz[i] * pa_y[i];

        tr_z_xyzzzz_yyyyyy[i] = tr_z_yzzzz_yyyyyy[i] * pa_x[i];

        tr_z_xyzzzz_yyyyyz[i] = tr_z_yzzzz_yyyyyz[i] * pa_x[i];

        tr_z_xyzzzz_yyyyzz[i] = tr_z_yzzzz_yyyyzz[i] * pa_x[i];

        tr_z_xyzzzz_yyyzzz[i] = tr_z_yzzzz_yyyzzz[i] * pa_x[i];

        tr_z_xyzzzz_yyzzzz[i] = tr_z_yzzzz_yyzzzz[i] * pa_x[i];

        tr_z_xyzzzz_yzzzzz[i] = tr_z_yzzzz_yzzzzz[i] * pa_x[i];

        tr_z_xyzzzz_zzzzzz[i] = tr_z_yzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 2128-2156 components of targeted buffer : II

    auto tr_z_xzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2128);

    auto tr_z_xzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2129);

    auto tr_z_xzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2130);

    auto tr_z_xzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2131);

    auto tr_z_xzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2132);

    auto tr_z_xzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2133);

    auto tr_z_xzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2134);

    auto tr_z_xzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2135);

    auto tr_z_xzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2136);

    auto tr_z_xzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2137);

    auto tr_z_xzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2138);

    auto tr_z_xzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2139);

    auto tr_z_xzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2140);

    auto tr_z_xzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2141);

    auto tr_z_xzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2142);

    auto tr_z_xzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2143);

    auto tr_z_xzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2144);

    auto tr_z_xzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2145);

    auto tr_z_xzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2146);

    auto tr_z_xzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2147);

    auto tr_z_xzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2148);

    auto tr_z_xzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2149);

    auto tr_z_xzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2150);

    auto tr_z_xzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2151);

    auto tr_z_xzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2152);

    auto tr_z_xzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2153);

    auto tr_z_xzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2154);

    auto tr_z_xzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2155);

    #pragma omp simd aligned(pa_x, tr_z_xzzzzz_xxxxxx, tr_z_xzzzzz_xxxxxy, tr_z_xzzzzz_xxxxxz, tr_z_xzzzzz_xxxxyy, tr_z_xzzzzz_xxxxyz, tr_z_xzzzzz_xxxxzz, tr_z_xzzzzz_xxxyyy, tr_z_xzzzzz_xxxyyz, tr_z_xzzzzz_xxxyzz, tr_z_xzzzzz_xxxzzz, tr_z_xzzzzz_xxyyyy, tr_z_xzzzzz_xxyyyz, tr_z_xzzzzz_xxyyzz, tr_z_xzzzzz_xxyzzz, tr_z_xzzzzz_xxzzzz, tr_z_xzzzzz_xyyyyy, tr_z_xzzzzz_xyyyyz, tr_z_xzzzzz_xyyyzz, tr_z_xzzzzz_xyyzzz, tr_z_xzzzzz_xyzzzz, tr_z_xzzzzz_xzzzzz, tr_z_xzzzzz_yyyyyy, tr_z_xzzzzz_yyyyyz, tr_z_xzzzzz_yyyyzz, tr_z_xzzzzz_yyyzzz, tr_z_xzzzzz_yyzzzz, tr_z_xzzzzz_yzzzzz, tr_z_xzzzzz_zzzzzz, tr_z_zzzzz_xxxxx, tr_z_zzzzz_xxxxxx, tr_z_zzzzz_xxxxxy, tr_z_zzzzz_xxxxxz, tr_z_zzzzz_xxxxy, tr_z_zzzzz_xxxxyy, tr_z_zzzzz_xxxxyz, tr_z_zzzzz_xxxxz, tr_z_zzzzz_xxxxzz, tr_z_zzzzz_xxxyy, tr_z_zzzzz_xxxyyy, tr_z_zzzzz_xxxyyz, tr_z_zzzzz_xxxyz, tr_z_zzzzz_xxxyzz, tr_z_zzzzz_xxxzz, tr_z_zzzzz_xxxzzz, tr_z_zzzzz_xxyyy, tr_z_zzzzz_xxyyyy, tr_z_zzzzz_xxyyyz, tr_z_zzzzz_xxyyz, tr_z_zzzzz_xxyyzz, tr_z_zzzzz_xxyzz, tr_z_zzzzz_xxyzzz, tr_z_zzzzz_xxzzz, tr_z_zzzzz_xxzzzz, tr_z_zzzzz_xyyyy, tr_z_zzzzz_xyyyyy, tr_z_zzzzz_xyyyyz, tr_z_zzzzz_xyyyz, tr_z_zzzzz_xyyyzz, tr_z_zzzzz_xyyzz, tr_z_zzzzz_xyyzzz, tr_z_zzzzz_xyzzz, tr_z_zzzzz_xyzzzz, tr_z_zzzzz_xzzzz, tr_z_zzzzz_xzzzzz, tr_z_zzzzz_yyyyy, tr_z_zzzzz_yyyyyy, tr_z_zzzzz_yyyyyz, tr_z_zzzzz_yyyyz, tr_z_zzzzz_yyyyzz, tr_z_zzzzz_yyyzz, tr_z_zzzzz_yyyzzz, tr_z_zzzzz_yyzzz, tr_z_zzzzz_yyzzzz, tr_z_zzzzz_yzzzz, tr_z_zzzzz_yzzzzz, tr_z_zzzzz_zzzzz, tr_z_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzzz_xxxxxx[i] = 6.0 * tr_z_zzzzz_xxxxx[i] * fe_0 + tr_z_zzzzz_xxxxxx[i] * pa_x[i];

        tr_z_xzzzzz_xxxxxy[i] = 5.0 * tr_z_zzzzz_xxxxy[i] * fe_0 + tr_z_zzzzz_xxxxxy[i] * pa_x[i];

        tr_z_xzzzzz_xxxxxz[i] = 5.0 * tr_z_zzzzz_xxxxz[i] * fe_0 + tr_z_zzzzz_xxxxxz[i] * pa_x[i];

        tr_z_xzzzzz_xxxxyy[i] = 4.0 * tr_z_zzzzz_xxxyy[i] * fe_0 + tr_z_zzzzz_xxxxyy[i] * pa_x[i];

        tr_z_xzzzzz_xxxxyz[i] = 4.0 * tr_z_zzzzz_xxxyz[i] * fe_0 + tr_z_zzzzz_xxxxyz[i] * pa_x[i];

        tr_z_xzzzzz_xxxxzz[i] = 4.0 * tr_z_zzzzz_xxxzz[i] * fe_0 + tr_z_zzzzz_xxxxzz[i] * pa_x[i];

        tr_z_xzzzzz_xxxyyy[i] = 3.0 * tr_z_zzzzz_xxyyy[i] * fe_0 + tr_z_zzzzz_xxxyyy[i] * pa_x[i];

        tr_z_xzzzzz_xxxyyz[i] = 3.0 * tr_z_zzzzz_xxyyz[i] * fe_0 + tr_z_zzzzz_xxxyyz[i] * pa_x[i];

        tr_z_xzzzzz_xxxyzz[i] = 3.0 * tr_z_zzzzz_xxyzz[i] * fe_0 + tr_z_zzzzz_xxxyzz[i] * pa_x[i];

        tr_z_xzzzzz_xxxzzz[i] = 3.0 * tr_z_zzzzz_xxzzz[i] * fe_0 + tr_z_zzzzz_xxxzzz[i] * pa_x[i];

        tr_z_xzzzzz_xxyyyy[i] = 2.0 * tr_z_zzzzz_xyyyy[i] * fe_0 + tr_z_zzzzz_xxyyyy[i] * pa_x[i];

        tr_z_xzzzzz_xxyyyz[i] = 2.0 * tr_z_zzzzz_xyyyz[i] * fe_0 + tr_z_zzzzz_xxyyyz[i] * pa_x[i];

        tr_z_xzzzzz_xxyyzz[i] = 2.0 * tr_z_zzzzz_xyyzz[i] * fe_0 + tr_z_zzzzz_xxyyzz[i] * pa_x[i];

        tr_z_xzzzzz_xxyzzz[i] = 2.0 * tr_z_zzzzz_xyzzz[i] * fe_0 + tr_z_zzzzz_xxyzzz[i] * pa_x[i];

        tr_z_xzzzzz_xxzzzz[i] = 2.0 * tr_z_zzzzz_xzzzz[i] * fe_0 + tr_z_zzzzz_xxzzzz[i] * pa_x[i];

        tr_z_xzzzzz_xyyyyy[i] = tr_z_zzzzz_yyyyy[i] * fe_0 + tr_z_zzzzz_xyyyyy[i] * pa_x[i];

        tr_z_xzzzzz_xyyyyz[i] = tr_z_zzzzz_yyyyz[i] * fe_0 + tr_z_zzzzz_xyyyyz[i] * pa_x[i];

        tr_z_xzzzzz_xyyyzz[i] = tr_z_zzzzz_yyyzz[i] * fe_0 + tr_z_zzzzz_xyyyzz[i] * pa_x[i];

        tr_z_xzzzzz_xyyzzz[i] = tr_z_zzzzz_yyzzz[i] * fe_0 + tr_z_zzzzz_xyyzzz[i] * pa_x[i];

        tr_z_xzzzzz_xyzzzz[i] = tr_z_zzzzz_yzzzz[i] * fe_0 + tr_z_zzzzz_xyzzzz[i] * pa_x[i];

        tr_z_xzzzzz_xzzzzz[i] = tr_z_zzzzz_zzzzz[i] * fe_0 + tr_z_zzzzz_xzzzzz[i] * pa_x[i];

        tr_z_xzzzzz_yyyyyy[i] = tr_z_zzzzz_yyyyyy[i] * pa_x[i];

        tr_z_xzzzzz_yyyyyz[i] = tr_z_zzzzz_yyyyyz[i] * pa_x[i];

        tr_z_xzzzzz_yyyyzz[i] = tr_z_zzzzz_yyyyzz[i] * pa_x[i];

        tr_z_xzzzzz_yyyzzz[i] = tr_z_zzzzz_yyyzzz[i] * pa_x[i];

        tr_z_xzzzzz_yyzzzz[i] = tr_z_zzzzz_yyzzzz[i] * pa_x[i];

        tr_z_xzzzzz_yzzzzz[i] = tr_z_zzzzz_yzzzzz[i] * pa_x[i];

        tr_z_xzzzzz_zzzzzz[i] = tr_z_zzzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 2156-2184 components of targeted buffer : II

    auto tr_z_yyyyyy_xxxxxx = pbuffer.data(idx_dip_ii + 2156);

    auto tr_z_yyyyyy_xxxxxy = pbuffer.data(idx_dip_ii + 2157);

    auto tr_z_yyyyyy_xxxxxz = pbuffer.data(idx_dip_ii + 2158);

    auto tr_z_yyyyyy_xxxxyy = pbuffer.data(idx_dip_ii + 2159);

    auto tr_z_yyyyyy_xxxxyz = pbuffer.data(idx_dip_ii + 2160);

    auto tr_z_yyyyyy_xxxxzz = pbuffer.data(idx_dip_ii + 2161);

    auto tr_z_yyyyyy_xxxyyy = pbuffer.data(idx_dip_ii + 2162);

    auto tr_z_yyyyyy_xxxyyz = pbuffer.data(idx_dip_ii + 2163);

    auto tr_z_yyyyyy_xxxyzz = pbuffer.data(idx_dip_ii + 2164);

    auto tr_z_yyyyyy_xxxzzz = pbuffer.data(idx_dip_ii + 2165);

    auto tr_z_yyyyyy_xxyyyy = pbuffer.data(idx_dip_ii + 2166);

    auto tr_z_yyyyyy_xxyyyz = pbuffer.data(idx_dip_ii + 2167);

    auto tr_z_yyyyyy_xxyyzz = pbuffer.data(idx_dip_ii + 2168);

    auto tr_z_yyyyyy_xxyzzz = pbuffer.data(idx_dip_ii + 2169);

    auto tr_z_yyyyyy_xxzzzz = pbuffer.data(idx_dip_ii + 2170);

    auto tr_z_yyyyyy_xyyyyy = pbuffer.data(idx_dip_ii + 2171);

    auto tr_z_yyyyyy_xyyyyz = pbuffer.data(idx_dip_ii + 2172);

    auto tr_z_yyyyyy_xyyyzz = pbuffer.data(idx_dip_ii + 2173);

    auto tr_z_yyyyyy_xyyzzz = pbuffer.data(idx_dip_ii + 2174);

    auto tr_z_yyyyyy_xyzzzz = pbuffer.data(idx_dip_ii + 2175);

    auto tr_z_yyyyyy_xzzzzz = pbuffer.data(idx_dip_ii + 2176);

    auto tr_z_yyyyyy_yyyyyy = pbuffer.data(idx_dip_ii + 2177);

    auto tr_z_yyyyyy_yyyyyz = pbuffer.data(idx_dip_ii + 2178);

    auto tr_z_yyyyyy_yyyyzz = pbuffer.data(idx_dip_ii + 2179);

    auto tr_z_yyyyyy_yyyzzz = pbuffer.data(idx_dip_ii + 2180);

    auto tr_z_yyyyyy_yyzzzz = pbuffer.data(idx_dip_ii + 2181);

    auto tr_z_yyyyyy_yzzzzz = pbuffer.data(idx_dip_ii + 2182);

    auto tr_z_yyyyyy_zzzzzz = pbuffer.data(idx_dip_ii + 2183);

    #pragma omp simd aligned(pa_y, tr_z_yyyy_xxxxxx, tr_z_yyyy_xxxxxy, tr_z_yyyy_xxxxxz, tr_z_yyyy_xxxxyy, tr_z_yyyy_xxxxyz, tr_z_yyyy_xxxxzz, tr_z_yyyy_xxxyyy, tr_z_yyyy_xxxyyz, tr_z_yyyy_xxxyzz, tr_z_yyyy_xxxzzz, tr_z_yyyy_xxyyyy, tr_z_yyyy_xxyyyz, tr_z_yyyy_xxyyzz, tr_z_yyyy_xxyzzz, tr_z_yyyy_xxzzzz, tr_z_yyyy_xyyyyy, tr_z_yyyy_xyyyyz, tr_z_yyyy_xyyyzz, tr_z_yyyy_xyyzzz, tr_z_yyyy_xyzzzz, tr_z_yyyy_xzzzzz, tr_z_yyyy_yyyyyy, tr_z_yyyy_yyyyyz, tr_z_yyyy_yyyyzz, tr_z_yyyy_yyyzzz, tr_z_yyyy_yyzzzz, tr_z_yyyy_yzzzzz, tr_z_yyyy_zzzzzz, tr_z_yyyyy_xxxxx, tr_z_yyyyy_xxxxxx, tr_z_yyyyy_xxxxxy, tr_z_yyyyy_xxxxxz, tr_z_yyyyy_xxxxy, tr_z_yyyyy_xxxxyy, tr_z_yyyyy_xxxxyz, tr_z_yyyyy_xxxxz, tr_z_yyyyy_xxxxzz, tr_z_yyyyy_xxxyy, tr_z_yyyyy_xxxyyy, tr_z_yyyyy_xxxyyz, tr_z_yyyyy_xxxyz, tr_z_yyyyy_xxxyzz, tr_z_yyyyy_xxxzz, tr_z_yyyyy_xxxzzz, tr_z_yyyyy_xxyyy, tr_z_yyyyy_xxyyyy, tr_z_yyyyy_xxyyyz, tr_z_yyyyy_xxyyz, tr_z_yyyyy_xxyyzz, tr_z_yyyyy_xxyzz, tr_z_yyyyy_xxyzzz, tr_z_yyyyy_xxzzz, tr_z_yyyyy_xxzzzz, tr_z_yyyyy_xyyyy, tr_z_yyyyy_xyyyyy, tr_z_yyyyy_xyyyyz, tr_z_yyyyy_xyyyz, tr_z_yyyyy_xyyyzz, tr_z_yyyyy_xyyzz, tr_z_yyyyy_xyyzzz, tr_z_yyyyy_xyzzz, tr_z_yyyyy_xyzzzz, tr_z_yyyyy_xzzzz, tr_z_yyyyy_xzzzzz, tr_z_yyyyy_yyyyy, tr_z_yyyyy_yyyyyy, tr_z_yyyyy_yyyyyz, tr_z_yyyyy_yyyyz, tr_z_yyyyy_yyyyzz, tr_z_yyyyy_yyyzz, tr_z_yyyyy_yyyzzz, tr_z_yyyyy_yyzzz, tr_z_yyyyy_yyzzzz, tr_z_yyyyy_yzzzz, tr_z_yyyyy_yzzzzz, tr_z_yyyyy_zzzzz, tr_z_yyyyy_zzzzzz, tr_z_yyyyyy_xxxxxx, tr_z_yyyyyy_xxxxxy, tr_z_yyyyyy_xxxxxz, tr_z_yyyyyy_xxxxyy, tr_z_yyyyyy_xxxxyz, tr_z_yyyyyy_xxxxzz, tr_z_yyyyyy_xxxyyy, tr_z_yyyyyy_xxxyyz, tr_z_yyyyyy_xxxyzz, tr_z_yyyyyy_xxxzzz, tr_z_yyyyyy_xxyyyy, tr_z_yyyyyy_xxyyyz, tr_z_yyyyyy_xxyyzz, tr_z_yyyyyy_xxyzzz, tr_z_yyyyyy_xxzzzz, tr_z_yyyyyy_xyyyyy, tr_z_yyyyyy_xyyyyz, tr_z_yyyyyy_xyyyzz, tr_z_yyyyyy_xyyzzz, tr_z_yyyyyy_xyzzzz, tr_z_yyyyyy_xzzzzz, tr_z_yyyyyy_yyyyyy, tr_z_yyyyyy_yyyyyz, tr_z_yyyyyy_yyyyzz, tr_z_yyyyyy_yyyzzz, tr_z_yyyyyy_yyzzzz, tr_z_yyyyyy_yzzzzz, tr_z_yyyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyy_xxxxxx[i] = 5.0 * tr_z_yyyy_xxxxxx[i] * fe_0 + tr_z_yyyyy_xxxxxx[i] * pa_y[i];

        tr_z_yyyyyy_xxxxxy[i] = 5.0 * tr_z_yyyy_xxxxxy[i] * fe_0 + tr_z_yyyyy_xxxxx[i] * fe_0 + tr_z_yyyyy_xxxxxy[i] * pa_y[i];

        tr_z_yyyyyy_xxxxxz[i] = 5.0 * tr_z_yyyy_xxxxxz[i] * fe_0 + tr_z_yyyyy_xxxxxz[i] * pa_y[i];

        tr_z_yyyyyy_xxxxyy[i] = 5.0 * tr_z_yyyy_xxxxyy[i] * fe_0 + 2.0 * tr_z_yyyyy_xxxxy[i] * fe_0 + tr_z_yyyyy_xxxxyy[i] * pa_y[i];

        tr_z_yyyyyy_xxxxyz[i] = 5.0 * tr_z_yyyy_xxxxyz[i] * fe_0 + tr_z_yyyyy_xxxxz[i] * fe_0 + tr_z_yyyyy_xxxxyz[i] * pa_y[i];

        tr_z_yyyyyy_xxxxzz[i] = 5.0 * tr_z_yyyy_xxxxzz[i] * fe_0 + tr_z_yyyyy_xxxxzz[i] * pa_y[i];

        tr_z_yyyyyy_xxxyyy[i] = 5.0 * tr_z_yyyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_yyyyy_xxxyy[i] * fe_0 + tr_z_yyyyy_xxxyyy[i] * pa_y[i];

        tr_z_yyyyyy_xxxyyz[i] = 5.0 * tr_z_yyyy_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyyyy_xxxyz[i] * fe_0 + tr_z_yyyyy_xxxyyz[i] * pa_y[i];

        tr_z_yyyyyy_xxxyzz[i] = 5.0 * tr_z_yyyy_xxxyzz[i] * fe_0 + tr_z_yyyyy_xxxzz[i] * fe_0 + tr_z_yyyyy_xxxyzz[i] * pa_y[i];

        tr_z_yyyyyy_xxxzzz[i] = 5.0 * tr_z_yyyy_xxxzzz[i] * fe_0 + tr_z_yyyyy_xxxzzz[i] * pa_y[i];

        tr_z_yyyyyy_xxyyyy[i] = 5.0 * tr_z_yyyy_xxyyyy[i] * fe_0 + 4.0 * tr_z_yyyyy_xxyyy[i] * fe_0 + tr_z_yyyyy_xxyyyy[i] * pa_y[i];

        tr_z_yyyyyy_xxyyyz[i] = 5.0 * tr_z_yyyy_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyyyy_xxyyz[i] * fe_0 + tr_z_yyyyy_xxyyyz[i] * pa_y[i];

        tr_z_yyyyyy_xxyyzz[i] = 5.0 * tr_z_yyyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyyyy_xxyzz[i] * fe_0 + tr_z_yyyyy_xxyyzz[i] * pa_y[i];

        tr_z_yyyyyy_xxyzzz[i] = 5.0 * tr_z_yyyy_xxyzzz[i] * fe_0 + tr_z_yyyyy_xxzzz[i] * fe_0 + tr_z_yyyyy_xxyzzz[i] * pa_y[i];

        tr_z_yyyyyy_xxzzzz[i] = 5.0 * tr_z_yyyy_xxzzzz[i] * fe_0 + tr_z_yyyyy_xxzzzz[i] * pa_y[i];

        tr_z_yyyyyy_xyyyyy[i] = 5.0 * tr_z_yyyy_xyyyyy[i] * fe_0 + 5.0 * tr_z_yyyyy_xyyyy[i] * fe_0 + tr_z_yyyyy_xyyyyy[i] * pa_y[i];

        tr_z_yyyyyy_xyyyyz[i] = 5.0 * tr_z_yyyy_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyyyy_xyyyz[i] * fe_0 + tr_z_yyyyy_xyyyyz[i] * pa_y[i];

        tr_z_yyyyyy_xyyyzz[i] = 5.0 * tr_z_yyyy_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyyyy_xyyzz[i] * fe_0 + tr_z_yyyyy_xyyyzz[i] * pa_y[i];

        tr_z_yyyyyy_xyyzzz[i] = 5.0 * tr_z_yyyy_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyyyy_xyzzz[i] * fe_0 + tr_z_yyyyy_xyyzzz[i] * pa_y[i];

        tr_z_yyyyyy_xyzzzz[i] = 5.0 * tr_z_yyyy_xyzzzz[i] * fe_0 + tr_z_yyyyy_xzzzz[i] * fe_0 + tr_z_yyyyy_xyzzzz[i] * pa_y[i];

        tr_z_yyyyyy_xzzzzz[i] = 5.0 * tr_z_yyyy_xzzzzz[i] * fe_0 + tr_z_yyyyy_xzzzzz[i] * pa_y[i];

        tr_z_yyyyyy_yyyyyy[i] = 5.0 * tr_z_yyyy_yyyyyy[i] * fe_0 + 6.0 * tr_z_yyyyy_yyyyy[i] * fe_0 + tr_z_yyyyy_yyyyyy[i] * pa_y[i];

        tr_z_yyyyyy_yyyyyz[i] = 5.0 * tr_z_yyyy_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyyyy_yyyyz[i] * fe_0 + tr_z_yyyyy_yyyyyz[i] * pa_y[i];

        tr_z_yyyyyy_yyyyzz[i] = 5.0 * tr_z_yyyy_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyyyy_yyyzz[i] * fe_0 + tr_z_yyyyy_yyyyzz[i] * pa_y[i];

        tr_z_yyyyyy_yyyzzz[i] = 5.0 * tr_z_yyyy_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyyyy_yyzzz[i] * fe_0 + tr_z_yyyyy_yyyzzz[i] * pa_y[i];

        tr_z_yyyyyy_yyzzzz[i] = 5.0 * tr_z_yyyy_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyyyy_yzzzz[i] * fe_0 + tr_z_yyyyy_yyzzzz[i] * pa_y[i];

        tr_z_yyyyyy_yzzzzz[i] = 5.0 * tr_z_yyyy_yzzzzz[i] * fe_0 + tr_z_yyyyy_zzzzz[i] * fe_0 + tr_z_yyyyy_yzzzzz[i] * pa_y[i];

        tr_z_yyyyyy_zzzzzz[i] = 5.0 * tr_z_yyyy_zzzzzz[i] * fe_0 + tr_z_yyyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 2184-2212 components of targeted buffer : II

    auto tr_z_yyyyyz_xxxxxx = pbuffer.data(idx_dip_ii + 2184);

    auto tr_z_yyyyyz_xxxxxy = pbuffer.data(idx_dip_ii + 2185);

    auto tr_z_yyyyyz_xxxxxz = pbuffer.data(idx_dip_ii + 2186);

    auto tr_z_yyyyyz_xxxxyy = pbuffer.data(idx_dip_ii + 2187);

    auto tr_z_yyyyyz_xxxxyz = pbuffer.data(idx_dip_ii + 2188);

    auto tr_z_yyyyyz_xxxxzz = pbuffer.data(idx_dip_ii + 2189);

    auto tr_z_yyyyyz_xxxyyy = pbuffer.data(idx_dip_ii + 2190);

    auto tr_z_yyyyyz_xxxyyz = pbuffer.data(idx_dip_ii + 2191);

    auto tr_z_yyyyyz_xxxyzz = pbuffer.data(idx_dip_ii + 2192);

    auto tr_z_yyyyyz_xxxzzz = pbuffer.data(idx_dip_ii + 2193);

    auto tr_z_yyyyyz_xxyyyy = pbuffer.data(idx_dip_ii + 2194);

    auto tr_z_yyyyyz_xxyyyz = pbuffer.data(idx_dip_ii + 2195);

    auto tr_z_yyyyyz_xxyyzz = pbuffer.data(idx_dip_ii + 2196);

    auto tr_z_yyyyyz_xxyzzz = pbuffer.data(idx_dip_ii + 2197);

    auto tr_z_yyyyyz_xxzzzz = pbuffer.data(idx_dip_ii + 2198);

    auto tr_z_yyyyyz_xyyyyy = pbuffer.data(idx_dip_ii + 2199);

    auto tr_z_yyyyyz_xyyyyz = pbuffer.data(idx_dip_ii + 2200);

    auto tr_z_yyyyyz_xyyyzz = pbuffer.data(idx_dip_ii + 2201);

    auto tr_z_yyyyyz_xyyzzz = pbuffer.data(idx_dip_ii + 2202);

    auto tr_z_yyyyyz_xyzzzz = pbuffer.data(idx_dip_ii + 2203);

    auto tr_z_yyyyyz_xzzzzz = pbuffer.data(idx_dip_ii + 2204);

    auto tr_z_yyyyyz_yyyyyy = pbuffer.data(idx_dip_ii + 2205);

    auto tr_z_yyyyyz_yyyyyz = pbuffer.data(idx_dip_ii + 2206);

    auto tr_z_yyyyyz_yyyyzz = pbuffer.data(idx_dip_ii + 2207);

    auto tr_z_yyyyyz_yyyzzz = pbuffer.data(idx_dip_ii + 2208);

    auto tr_z_yyyyyz_yyzzzz = pbuffer.data(idx_dip_ii + 2209);

    auto tr_z_yyyyyz_yzzzzz = pbuffer.data(idx_dip_ii + 2210);

    auto tr_z_yyyyyz_zzzzzz = pbuffer.data(idx_dip_ii + 2211);

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yyyyy_xxxxxy, tr_z_yyyyy_xxxxyy, tr_z_yyyyy_xxxyyy, tr_z_yyyyy_xxyyyy, tr_z_yyyyy_xyyyyy, tr_z_yyyyy_yyyyyy, tr_z_yyyyyz_xxxxxx, tr_z_yyyyyz_xxxxxy, tr_z_yyyyyz_xxxxxz, tr_z_yyyyyz_xxxxyy, tr_z_yyyyyz_xxxxyz, tr_z_yyyyyz_xxxxzz, tr_z_yyyyyz_xxxyyy, tr_z_yyyyyz_xxxyyz, tr_z_yyyyyz_xxxyzz, tr_z_yyyyyz_xxxzzz, tr_z_yyyyyz_xxyyyy, tr_z_yyyyyz_xxyyyz, tr_z_yyyyyz_xxyyzz, tr_z_yyyyyz_xxyzzz, tr_z_yyyyyz_xxzzzz, tr_z_yyyyyz_xyyyyy, tr_z_yyyyyz_xyyyyz, tr_z_yyyyyz_xyyyzz, tr_z_yyyyyz_xyyzzz, tr_z_yyyyyz_xyzzzz, tr_z_yyyyyz_xzzzzz, tr_z_yyyyyz_yyyyyy, tr_z_yyyyyz_yyyyyz, tr_z_yyyyyz_yyyyzz, tr_z_yyyyyz_yyyzzz, tr_z_yyyyyz_yyzzzz, tr_z_yyyyyz_yzzzzz, tr_z_yyyyyz_zzzzzz, tr_z_yyyyz_xxxxxx, tr_z_yyyyz_xxxxxz, tr_z_yyyyz_xxxxyz, tr_z_yyyyz_xxxxz, tr_z_yyyyz_xxxxzz, tr_z_yyyyz_xxxyyz, tr_z_yyyyz_xxxyz, tr_z_yyyyz_xxxyzz, tr_z_yyyyz_xxxzz, tr_z_yyyyz_xxxzzz, tr_z_yyyyz_xxyyyz, tr_z_yyyyz_xxyyz, tr_z_yyyyz_xxyyzz, tr_z_yyyyz_xxyzz, tr_z_yyyyz_xxyzzz, tr_z_yyyyz_xxzzz, tr_z_yyyyz_xxzzzz, tr_z_yyyyz_xyyyyz, tr_z_yyyyz_xyyyz, tr_z_yyyyz_xyyyzz, tr_z_yyyyz_xyyzz, tr_z_yyyyz_xyyzzz, tr_z_yyyyz_xyzzz, tr_z_yyyyz_xyzzzz, tr_z_yyyyz_xzzzz, tr_z_yyyyz_xzzzzz, tr_z_yyyyz_yyyyyz, tr_z_yyyyz_yyyyz, tr_z_yyyyz_yyyyzz, tr_z_yyyyz_yyyzz, tr_z_yyyyz_yyyzzz, tr_z_yyyyz_yyzzz, tr_z_yyyyz_yyzzzz, tr_z_yyyyz_yzzzz, tr_z_yyyyz_yzzzzz, tr_z_yyyyz_zzzzz, tr_z_yyyyz_zzzzzz, tr_z_yyyz_xxxxxx, tr_z_yyyz_xxxxxz, tr_z_yyyz_xxxxyz, tr_z_yyyz_xxxxzz, tr_z_yyyz_xxxyyz, tr_z_yyyz_xxxyzz, tr_z_yyyz_xxxzzz, tr_z_yyyz_xxyyyz, tr_z_yyyz_xxyyzz, tr_z_yyyz_xxyzzz, tr_z_yyyz_xxzzzz, tr_z_yyyz_xyyyyz, tr_z_yyyz_xyyyzz, tr_z_yyyz_xyyzzz, tr_z_yyyz_xyzzzz, tr_z_yyyz_xzzzzz, tr_z_yyyz_yyyyyz, tr_z_yyyz_yyyyzz, tr_z_yyyz_yyyzzz, tr_z_yyyz_yyzzzz, tr_z_yyyz_yzzzzz, tr_z_yyyz_zzzzzz, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxyy, ts_yyyyy_xxxyyy, ts_yyyyy_xxyyyy, ts_yyyyy_xyyyyy, ts_yyyyy_yyyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyz_xxxxxx[i] = 4.0 * tr_z_yyyz_xxxxxx[i] * fe_0 + tr_z_yyyyz_xxxxxx[i] * pa_y[i];

        tr_z_yyyyyz_xxxxxy[i] = ts_yyyyy_xxxxxy[i] * fe_0 + tr_z_yyyyy_xxxxxy[i] * pa_z[i];

        tr_z_yyyyyz_xxxxxz[i] = 4.0 * tr_z_yyyz_xxxxxz[i] * fe_0 + tr_z_yyyyz_xxxxxz[i] * pa_y[i];

        tr_z_yyyyyz_xxxxyy[i] = ts_yyyyy_xxxxyy[i] * fe_0 + tr_z_yyyyy_xxxxyy[i] * pa_z[i];

        tr_z_yyyyyz_xxxxyz[i] = 4.0 * tr_z_yyyz_xxxxyz[i] * fe_0 + tr_z_yyyyz_xxxxz[i] * fe_0 + tr_z_yyyyz_xxxxyz[i] * pa_y[i];

        tr_z_yyyyyz_xxxxzz[i] = 4.0 * tr_z_yyyz_xxxxzz[i] * fe_0 + tr_z_yyyyz_xxxxzz[i] * pa_y[i];

        tr_z_yyyyyz_xxxyyy[i] = ts_yyyyy_xxxyyy[i] * fe_0 + tr_z_yyyyy_xxxyyy[i] * pa_z[i];

        tr_z_yyyyyz_xxxyyz[i] = 4.0 * tr_z_yyyz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyyyz_xxxyz[i] * fe_0 + tr_z_yyyyz_xxxyyz[i] * pa_y[i];

        tr_z_yyyyyz_xxxyzz[i] = 4.0 * tr_z_yyyz_xxxyzz[i] * fe_0 + tr_z_yyyyz_xxxzz[i] * fe_0 + tr_z_yyyyz_xxxyzz[i] * pa_y[i];

        tr_z_yyyyyz_xxxzzz[i] = 4.0 * tr_z_yyyz_xxxzzz[i] * fe_0 + tr_z_yyyyz_xxxzzz[i] * pa_y[i];

        tr_z_yyyyyz_xxyyyy[i] = ts_yyyyy_xxyyyy[i] * fe_0 + tr_z_yyyyy_xxyyyy[i] * pa_z[i];

        tr_z_yyyyyz_xxyyyz[i] = 4.0 * tr_z_yyyz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyyyz_xxyyz[i] * fe_0 + tr_z_yyyyz_xxyyyz[i] * pa_y[i];

        tr_z_yyyyyz_xxyyzz[i] = 4.0 * tr_z_yyyz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyyyz_xxyzz[i] * fe_0 + tr_z_yyyyz_xxyyzz[i] * pa_y[i];

        tr_z_yyyyyz_xxyzzz[i] = 4.0 * tr_z_yyyz_xxyzzz[i] * fe_0 + tr_z_yyyyz_xxzzz[i] * fe_0 + tr_z_yyyyz_xxyzzz[i] * pa_y[i];

        tr_z_yyyyyz_xxzzzz[i] = 4.0 * tr_z_yyyz_xxzzzz[i] * fe_0 + tr_z_yyyyz_xxzzzz[i] * pa_y[i];

        tr_z_yyyyyz_xyyyyy[i] = ts_yyyyy_xyyyyy[i] * fe_0 + tr_z_yyyyy_xyyyyy[i] * pa_z[i];

        tr_z_yyyyyz_xyyyyz[i] = 4.0 * tr_z_yyyz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyyyz_xyyyz[i] * fe_0 + tr_z_yyyyz_xyyyyz[i] * pa_y[i];

        tr_z_yyyyyz_xyyyzz[i] = 4.0 * tr_z_yyyz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyyyz_xyyzz[i] * fe_0 + tr_z_yyyyz_xyyyzz[i] * pa_y[i];

        tr_z_yyyyyz_xyyzzz[i] = 4.0 * tr_z_yyyz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyyyz_xyzzz[i] * fe_0 + tr_z_yyyyz_xyyzzz[i] * pa_y[i];

        tr_z_yyyyyz_xyzzzz[i] = 4.0 * tr_z_yyyz_xyzzzz[i] * fe_0 + tr_z_yyyyz_xzzzz[i] * fe_0 + tr_z_yyyyz_xyzzzz[i] * pa_y[i];

        tr_z_yyyyyz_xzzzzz[i] = 4.0 * tr_z_yyyz_xzzzzz[i] * fe_0 + tr_z_yyyyz_xzzzzz[i] * pa_y[i];

        tr_z_yyyyyz_yyyyyy[i] = ts_yyyyy_yyyyyy[i] * fe_0 + tr_z_yyyyy_yyyyyy[i] * pa_z[i];

        tr_z_yyyyyz_yyyyyz[i] = 4.0 * tr_z_yyyz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyyyz_yyyyz[i] * fe_0 + tr_z_yyyyz_yyyyyz[i] * pa_y[i];

        tr_z_yyyyyz_yyyyzz[i] = 4.0 * tr_z_yyyz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyyyz_yyyzz[i] * fe_0 + tr_z_yyyyz_yyyyzz[i] * pa_y[i];

        tr_z_yyyyyz_yyyzzz[i] = 4.0 * tr_z_yyyz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyyyz_yyzzz[i] * fe_0 + tr_z_yyyyz_yyyzzz[i] * pa_y[i];

        tr_z_yyyyyz_yyzzzz[i] = 4.0 * tr_z_yyyz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyyyz_yzzzz[i] * fe_0 + tr_z_yyyyz_yyzzzz[i] * pa_y[i];

        tr_z_yyyyyz_yzzzzz[i] = 4.0 * tr_z_yyyz_yzzzzz[i] * fe_0 + tr_z_yyyyz_zzzzz[i] * fe_0 + tr_z_yyyyz_yzzzzz[i] * pa_y[i];

        tr_z_yyyyyz_zzzzzz[i] = 4.0 * tr_z_yyyz_zzzzzz[i] * fe_0 + tr_z_yyyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 2212-2240 components of targeted buffer : II

    auto tr_z_yyyyzz_xxxxxx = pbuffer.data(idx_dip_ii + 2212);

    auto tr_z_yyyyzz_xxxxxy = pbuffer.data(idx_dip_ii + 2213);

    auto tr_z_yyyyzz_xxxxxz = pbuffer.data(idx_dip_ii + 2214);

    auto tr_z_yyyyzz_xxxxyy = pbuffer.data(idx_dip_ii + 2215);

    auto tr_z_yyyyzz_xxxxyz = pbuffer.data(idx_dip_ii + 2216);

    auto tr_z_yyyyzz_xxxxzz = pbuffer.data(idx_dip_ii + 2217);

    auto tr_z_yyyyzz_xxxyyy = pbuffer.data(idx_dip_ii + 2218);

    auto tr_z_yyyyzz_xxxyyz = pbuffer.data(idx_dip_ii + 2219);

    auto tr_z_yyyyzz_xxxyzz = pbuffer.data(idx_dip_ii + 2220);

    auto tr_z_yyyyzz_xxxzzz = pbuffer.data(idx_dip_ii + 2221);

    auto tr_z_yyyyzz_xxyyyy = pbuffer.data(idx_dip_ii + 2222);

    auto tr_z_yyyyzz_xxyyyz = pbuffer.data(idx_dip_ii + 2223);

    auto tr_z_yyyyzz_xxyyzz = pbuffer.data(idx_dip_ii + 2224);

    auto tr_z_yyyyzz_xxyzzz = pbuffer.data(idx_dip_ii + 2225);

    auto tr_z_yyyyzz_xxzzzz = pbuffer.data(idx_dip_ii + 2226);

    auto tr_z_yyyyzz_xyyyyy = pbuffer.data(idx_dip_ii + 2227);

    auto tr_z_yyyyzz_xyyyyz = pbuffer.data(idx_dip_ii + 2228);

    auto tr_z_yyyyzz_xyyyzz = pbuffer.data(idx_dip_ii + 2229);

    auto tr_z_yyyyzz_xyyzzz = pbuffer.data(idx_dip_ii + 2230);

    auto tr_z_yyyyzz_xyzzzz = pbuffer.data(idx_dip_ii + 2231);

    auto tr_z_yyyyzz_xzzzzz = pbuffer.data(idx_dip_ii + 2232);

    auto tr_z_yyyyzz_yyyyyy = pbuffer.data(idx_dip_ii + 2233);

    auto tr_z_yyyyzz_yyyyyz = pbuffer.data(idx_dip_ii + 2234);

    auto tr_z_yyyyzz_yyyyzz = pbuffer.data(idx_dip_ii + 2235);

    auto tr_z_yyyyzz_yyyzzz = pbuffer.data(idx_dip_ii + 2236);

    auto tr_z_yyyyzz_yyzzzz = pbuffer.data(idx_dip_ii + 2237);

    auto tr_z_yyyyzz_yzzzzz = pbuffer.data(idx_dip_ii + 2238);

    auto tr_z_yyyyzz_zzzzzz = pbuffer.data(idx_dip_ii + 2239);

    #pragma omp simd aligned(pa_y, tr_z_yyyyzz_xxxxxx, tr_z_yyyyzz_xxxxxy, tr_z_yyyyzz_xxxxxz, tr_z_yyyyzz_xxxxyy, tr_z_yyyyzz_xxxxyz, tr_z_yyyyzz_xxxxzz, tr_z_yyyyzz_xxxyyy, tr_z_yyyyzz_xxxyyz, tr_z_yyyyzz_xxxyzz, tr_z_yyyyzz_xxxzzz, tr_z_yyyyzz_xxyyyy, tr_z_yyyyzz_xxyyyz, tr_z_yyyyzz_xxyyzz, tr_z_yyyyzz_xxyzzz, tr_z_yyyyzz_xxzzzz, tr_z_yyyyzz_xyyyyy, tr_z_yyyyzz_xyyyyz, tr_z_yyyyzz_xyyyzz, tr_z_yyyyzz_xyyzzz, tr_z_yyyyzz_xyzzzz, tr_z_yyyyzz_xzzzzz, tr_z_yyyyzz_yyyyyy, tr_z_yyyyzz_yyyyyz, tr_z_yyyyzz_yyyyzz, tr_z_yyyyzz_yyyzzz, tr_z_yyyyzz_yyzzzz, tr_z_yyyyzz_yzzzzz, tr_z_yyyyzz_zzzzzz, tr_z_yyyzz_xxxxx, tr_z_yyyzz_xxxxxx, tr_z_yyyzz_xxxxxy, tr_z_yyyzz_xxxxxz, tr_z_yyyzz_xxxxy, tr_z_yyyzz_xxxxyy, tr_z_yyyzz_xxxxyz, tr_z_yyyzz_xxxxz, tr_z_yyyzz_xxxxzz, tr_z_yyyzz_xxxyy, tr_z_yyyzz_xxxyyy, tr_z_yyyzz_xxxyyz, tr_z_yyyzz_xxxyz, tr_z_yyyzz_xxxyzz, tr_z_yyyzz_xxxzz, tr_z_yyyzz_xxxzzz, tr_z_yyyzz_xxyyy, tr_z_yyyzz_xxyyyy, tr_z_yyyzz_xxyyyz, tr_z_yyyzz_xxyyz, tr_z_yyyzz_xxyyzz, tr_z_yyyzz_xxyzz, tr_z_yyyzz_xxyzzz, tr_z_yyyzz_xxzzz, tr_z_yyyzz_xxzzzz, tr_z_yyyzz_xyyyy, tr_z_yyyzz_xyyyyy, tr_z_yyyzz_xyyyyz, tr_z_yyyzz_xyyyz, tr_z_yyyzz_xyyyzz, tr_z_yyyzz_xyyzz, tr_z_yyyzz_xyyzzz, tr_z_yyyzz_xyzzz, tr_z_yyyzz_xyzzzz, tr_z_yyyzz_xzzzz, tr_z_yyyzz_xzzzzz, tr_z_yyyzz_yyyyy, tr_z_yyyzz_yyyyyy, tr_z_yyyzz_yyyyyz, tr_z_yyyzz_yyyyz, tr_z_yyyzz_yyyyzz, tr_z_yyyzz_yyyzz, tr_z_yyyzz_yyyzzz, tr_z_yyyzz_yyzzz, tr_z_yyyzz_yyzzzz, tr_z_yyyzz_yzzzz, tr_z_yyyzz_yzzzzz, tr_z_yyyzz_zzzzz, tr_z_yyyzz_zzzzzz, tr_z_yyzz_xxxxxx, tr_z_yyzz_xxxxxy, tr_z_yyzz_xxxxxz, tr_z_yyzz_xxxxyy, tr_z_yyzz_xxxxyz, tr_z_yyzz_xxxxzz, tr_z_yyzz_xxxyyy, tr_z_yyzz_xxxyyz, tr_z_yyzz_xxxyzz, tr_z_yyzz_xxxzzz, tr_z_yyzz_xxyyyy, tr_z_yyzz_xxyyyz, tr_z_yyzz_xxyyzz, tr_z_yyzz_xxyzzz, tr_z_yyzz_xxzzzz, tr_z_yyzz_xyyyyy, tr_z_yyzz_xyyyyz, tr_z_yyzz_xyyyzz, tr_z_yyzz_xyyzzz, tr_z_yyzz_xyzzzz, tr_z_yyzz_xzzzzz, tr_z_yyzz_yyyyyy, tr_z_yyzz_yyyyyz, tr_z_yyzz_yyyyzz, tr_z_yyzz_yyyzzz, tr_z_yyzz_yyzzzz, tr_z_yyzz_yzzzzz, tr_z_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyzz_xxxxxx[i] = 3.0 * tr_z_yyzz_xxxxxx[i] * fe_0 + tr_z_yyyzz_xxxxxx[i] * pa_y[i];

        tr_z_yyyyzz_xxxxxy[i] = 3.0 * tr_z_yyzz_xxxxxy[i] * fe_0 + tr_z_yyyzz_xxxxx[i] * fe_0 + tr_z_yyyzz_xxxxxy[i] * pa_y[i];

        tr_z_yyyyzz_xxxxxz[i] = 3.0 * tr_z_yyzz_xxxxxz[i] * fe_0 + tr_z_yyyzz_xxxxxz[i] * pa_y[i];

        tr_z_yyyyzz_xxxxyy[i] = 3.0 * tr_z_yyzz_xxxxyy[i] * fe_0 + 2.0 * tr_z_yyyzz_xxxxy[i] * fe_0 + tr_z_yyyzz_xxxxyy[i] * pa_y[i];

        tr_z_yyyyzz_xxxxyz[i] = 3.0 * tr_z_yyzz_xxxxyz[i] * fe_0 + tr_z_yyyzz_xxxxz[i] * fe_0 + tr_z_yyyzz_xxxxyz[i] * pa_y[i];

        tr_z_yyyyzz_xxxxzz[i] = 3.0 * tr_z_yyzz_xxxxzz[i] * fe_0 + tr_z_yyyzz_xxxxzz[i] * pa_y[i];

        tr_z_yyyyzz_xxxyyy[i] = 3.0 * tr_z_yyzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_yyyzz_xxxyy[i] * fe_0 + tr_z_yyyzz_xxxyyy[i] * pa_y[i];

        tr_z_yyyyzz_xxxyyz[i] = 3.0 * tr_z_yyzz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyyzz_xxxyz[i] * fe_0 + tr_z_yyyzz_xxxyyz[i] * pa_y[i];

        tr_z_yyyyzz_xxxyzz[i] = 3.0 * tr_z_yyzz_xxxyzz[i] * fe_0 + tr_z_yyyzz_xxxzz[i] * fe_0 + tr_z_yyyzz_xxxyzz[i] * pa_y[i];

        tr_z_yyyyzz_xxxzzz[i] = 3.0 * tr_z_yyzz_xxxzzz[i] * fe_0 + tr_z_yyyzz_xxxzzz[i] * pa_y[i];

        tr_z_yyyyzz_xxyyyy[i] = 3.0 * tr_z_yyzz_xxyyyy[i] * fe_0 + 4.0 * tr_z_yyyzz_xxyyy[i] * fe_0 + tr_z_yyyzz_xxyyyy[i] * pa_y[i];

        tr_z_yyyyzz_xxyyyz[i] = 3.0 * tr_z_yyzz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyyzz_xxyyz[i] * fe_0 + tr_z_yyyzz_xxyyyz[i] * pa_y[i];

        tr_z_yyyyzz_xxyyzz[i] = 3.0 * tr_z_yyzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyyzz_xxyzz[i] * fe_0 + tr_z_yyyzz_xxyyzz[i] * pa_y[i];

        tr_z_yyyyzz_xxyzzz[i] = 3.0 * tr_z_yyzz_xxyzzz[i] * fe_0 + tr_z_yyyzz_xxzzz[i] * fe_0 + tr_z_yyyzz_xxyzzz[i] * pa_y[i];

        tr_z_yyyyzz_xxzzzz[i] = 3.0 * tr_z_yyzz_xxzzzz[i] * fe_0 + tr_z_yyyzz_xxzzzz[i] * pa_y[i];

        tr_z_yyyyzz_xyyyyy[i] = 3.0 * tr_z_yyzz_xyyyyy[i] * fe_0 + 5.0 * tr_z_yyyzz_xyyyy[i] * fe_0 + tr_z_yyyzz_xyyyyy[i] * pa_y[i];

        tr_z_yyyyzz_xyyyyz[i] = 3.0 * tr_z_yyzz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyyzz_xyyyz[i] * fe_0 + tr_z_yyyzz_xyyyyz[i] * pa_y[i];

        tr_z_yyyyzz_xyyyzz[i] = 3.0 * tr_z_yyzz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyyzz_xyyzz[i] * fe_0 + tr_z_yyyzz_xyyyzz[i] * pa_y[i];

        tr_z_yyyyzz_xyyzzz[i] = 3.0 * tr_z_yyzz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyyzz_xyzzz[i] * fe_0 + tr_z_yyyzz_xyyzzz[i] * pa_y[i];

        tr_z_yyyyzz_xyzzzz[i] = 3.0 * tr_z_yyzz_xyzzzz[i] * fe_0 + tr_z_yyyzz_xzzzz[i] * fe_0 + tr_z_yyyzz_xyzzzz[i] * pa_y[i];

        tr_z_yyyyzz_xzzzzz[i] = 3.0 * tr_z_yyzz_xzzzzz[i] * fe_0 + tr_z_yyyzz_xzzzzz[i] * pa_y[i];

        tr_z_yyyyzz_yyyyyy[i] = 3.0 * tr_z_yyzz_yyyyyy[i] * fe_0 + 6.0 * tr_z_yyyzz_yyyyy[i] * fe_0 + tr_z_yyyzz_yyyyyy[i] * pa_y[i];

        tr_z_yyyyzz_yyyyyz[i] = 3.0 * tr_z_yyzz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyyzz_yyyyz[i] * fe_0 + tr_z_yyyzz_yyyyyz[i] * pa_y[i];

        tr_z_yyyyzz_yyyyzz[i] = 3.0 * tr_z_yyzz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyyzz_yyyzz[i] * fe_0 + tr_z_yyyzz_yyyyzz[i] * pa_y[i];

        tr_z_yyyyzz_yyyzzz[i] = 3.0 * tr_z_yyzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyyzz_yyzzz[i] * fe_0 + tr_z_yyyzz_yyyzzz[i] * pa_y[i];

        tr_z_yyyyzz_yyzzzz[i] = 3.0 * tr_z_yyzz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyyzz_yzzzz[i] * fe_0 + tr_z_yyyzz_yyzzzz[i] * pa_y[i];

        tr_z_yyyyzz_yzzzzz[i] = 3.0 * tr_z_yyzz_yzzzzz[i] * fe_0 + tr_z_yyyzz_zzzzz[i] * fe_0 + tr_z_yyyzz_yzzzzz[i] * pa_y[i];

        tr_z_yyyyzz_zzzzzz[i] = 3.0 * tr_z_yyzz_zzzzzz[i] * fe_0 + tr_z_yyyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 2240-2268 components of targeted buffer : II

    auto tr_z_yyyzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2240);

    auto tr_z_yyyzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2241);

    auto tr_z_yyyzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2242);

    auto tr_z_yyyzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2243);

    auto tr_z_yyyzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2244);

    auto tr_z_yyyzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2245);

    auto tr_z_yyyzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2246);

    auto tr_z_yyyzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2247);

    auto tr_z_yyyzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2248);

    auto tr_z_yyyzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2249);

    auto tr_z_yyyzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2250);

    auto tr_z_yyyzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2251);

    auto tr_z_yyyzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2252);

    auto tr_z_yyyzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2253);

    auto tr_z_yyyzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2254);

    auto tr_z_yyyzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2255);

    auto tr_z_yyyzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2256);

    auto tr_z_yyyzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2257);

    auto tr_z_yyyzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2258);

    auto tr_z_yyyzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2259);

    auto tr_z_yyyzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2260);

    auto tr_z_yyyzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2261);

    auto tr_z_yyyzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2262);

    auto tr_z_yyyzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2263);

    auto tr_z_yyyzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2264);

    auto tr_z_yyyzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2265);

    auto tr_z_yyyzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2266);

    auto tr_z_yyyzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2267);

    #pragma omp simd aligned(pa_y, tr_z_yyyzzz_xxxxxx, tr_z_yyyzzz_xxxxxy, tr_z_yyyzzz_xxxxxz, tr_z_yyyzzz_xxxxyy, tr_z_yyyzzz_xxxxyz, tr_z_yyyzzz_xxxxzz, tr_z_yyyzzz_xxxyyy, tr_z_yyyzzz_xxxyyz, tr_z_yyyzzz_xxxyzz, tr_z_yyyzzz_xxxzzz, tr_z_yyyzzz_xxyyyy, tr_z_yyyzzz_xxyyyz, tr_z_yyyzzz_xxyyzz, tr_z_yyyzzz_xxyzzz, tr_z_yyyzzz_xxzzzz, tr_z_yyyzzz_xyyyyy, tr_z_yyyzzz_xyyyyz, tr_z_yyyzzz_xyyyzz, tr_z_yyyzzz_xyyzzz, tr_z_yyyzzz_xyzzzz, tr_z_yyyzzz_xzzzzz, tr_z_yyyzzz_yyyyyy, tr_z_yyyzzz_yyyyyz, tr_z_yyyzzz_yyyyzz, tr_z_yyyzzz_yyyzzz, tr_z_yyyzzz_yyzzzz, tr_z_yyyzzz_yzzzzz, tr_z_yyyzzz_zzzzzz, tr_z_yyzzz_xxxxx, tr_z_yyzzz_xxxxxx, tr_z_yyzzz_xxxxxy, tr_z_yyzzz_xxxxxz, tr_z_yyzzz_xxxxy, tr_z_yyzzz_xxxxyy, tr_z_yyzzz_xxxxyz, tr_z_yyzzz_xxxxz, tr_z_yyzzz_xxxxzz, tr_z_yyzzz_xxxyy, tr_z_yyzzz_xxxyyy, tr_z_yyzzz_xxxyyz, tr_z_yyzzz_xxxyz, tr_z_yyzzz_xxxyzz, tr_z_yyzzz_xxxzz, tr_z_yyzzz_xxxzzz, tr_z_yyzzz_xxyyy, tr_z_yyzzz_xxyyyy, tr_z_yyzzz_xxyyyz, tr_z_yyzzz_xxyyz, tr_z_yyzzz_xxyyzz, tr_z_yyzzz_xxyzz, tr_z_yyzzz_xxyzzz, tr_z_yyzzz_xxzzz, tr_z_yyzzz_xxzzzz, tr_z_yyzzz_xyyyy, tr_z_yyzzz_xyyyyy, tr_z_yyzzz_xyyyyz, tr_z_yyzzz_xyyyz, tr_z_yyzzz_xyyyzz, tr_z_yyzzz_xyyzz, tr_z_yyzzz_xyyzzz, tr_z_yyzzz_xyzzz, tr_z_yyzzz_xyzzzz, tr_z_yyzzz_xzzzz, tr_z_yyzzz_xzzzzz, tr_z_yyzzz_yyyyy, tr_z_yyzzz_yyyyyy, tr_z_yyzzz_yyyyyz, tr_z_yyzzz_yyyyz, tr_z_yyzzz_yyyyzz, tr_z_yyzzz_yyyzz, tr_z_yyzzz_yyyzzz, tr_z_yyzzz_yyzzz, tr_z_yyzzz_yyzzzz, tr_z_yyzzz_yzzzz, tr_z_yyzzz_yzzzzz, tr_z_yyzzz_zzzzz, tr_z_yyzzz_zzzzzz, tr_z_yzzz_xxxxxx, tr_z_yzzz_xxxxxy, tr_z_yzzz_xxxxxz, tr_z_yzzz_xxxxyy, tr_z_yzzz_xxxxyz, tr_z_yzzz_xxxxzz, tr_z_yzzz_xxxyyy, tr_z_yzzz_xxxyyz, tr_z_yzzz_xxxyzz, tr_z_yzzz_xxxzzz, tr_z_yzzz_xxyyyy, tr_z_yzzz_xxyyyz, tr_z_yzzz_xxyyzz, tr_z_yzzz_xxyzzz, tr_z_yzzz_xxzzzz, tr_z_yzzz_xyyyyy, tr_z_yzzz_xyyyyz, tr_z_yzzz_xyyyzz, tr_z_yzzz_xyyzzz, tr_z_yzzz_xyzzzz, tr_z_yzzz_xzzzzz, tr_z_yzzz_yyyyyy, tr_z_yzzz_yyyyyz, tr_z_yzzz_yyyyzz, tr_z_yzzz_yyyzzz, tr_z_yzzz_yyzzzz, tr_z_yzzz_yzzzzz, tr_z_yzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzzz_xxxxxx[i] = 2.0 * tr_z_yzzz_xxxxxx[i] * fe_0 + tr_z_yyzzz_xxxxxx[i] * pa_y[i];

        tr_z_yyyzzz_xxxxxy[i] = 2.0 * tr_z_yzzz_xxxxxy[i] * fe_0 + tr_z_yyzzz_xxxxx[i] * fe_0 + tr_z_yyzzz_xxxxxy[i] * pa_y[i];

        tr_z_yyyzzz_xxxxxz[i] = 2.0 * tr_z_yzzz_xxxxxz[i] * fe_0 + tr_z_yyzzz_xxxxxz[i] * pa_y[i];

        tr_z_yyyzzz_xxxxyy[i] = 2.0 * tr_z_yzzz_xxxxyy[i] * fe_0 + 2.0 * tr_z_yyzzz_xxxxy[i] * fe_0 + tr_z_yyzzz_xxxxyy[i] * pa_y[i];

        tr_z_yyyzzz_xxxxyz[i] = 2.0 * tr_z_yzzz_xxxxyz[i] * fe_0 + tr_z_yyzzz_xxxxz[i] * fe_0 + tr_z_yyzzz_xxxxyz[i] * pa_y[i];

        tr_z_yyyzzz_xxxxzz[i] = 2.0 * tr_z_yzzz_xxxxzz[i] * fe_0 + tr_z_yyzzz_xxxxzz[i] * pa_y[i];

        tr_z_yyyzzz_xxxyyy[i] = 2.0 * tr_z_yzzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_yyzzz_xxxyy[i] * fe_0 + tr_z_yyzzz_xxxyyy[i] * pa_y[i];

        tr_z_yyyzzz_xxxyyz[i] = 2.0 * tr_z_yzzz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyzzz_xxxyz[i] * fe_0 + tr_z_yyzzz_xxxyyz[i] * pa_y[i];

        tr_z_yyyzzz_xxxyzz[i] = 2.0 * tr_z_yzzz_xxxyzz[i] * fe_0 + tr_z_yyzzz_xxxzz[i] * fe_0 + tr_z_yyzzz_xxxyzz[i] * pa_y[i];

        tr_z_yyyzzz_xxxzzz[i] = 2.0 * tr_z_yzzz_xxxzzz[i] * fe_0 + tr_z_yyzzz_xxxzzz[i] * pa_y[i];

        tr_z_yyyzzz_xxyyyy[i] = 2.0 * tr_z_yzzz_xxyyyy[i] * fe_0 + 4.0 * tr_z_yyzzz_xxyyy[i] * fe_0 + tr_z_yyzzz_xxyyyy[i] * pa_y[i];

        tr_z_yyyzzz_xxyyyz[i] = 2.0 * tr_z_yzzz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyzzz_xxyyz[i] * fe_0 + tr_z_yyzzz_xxyyyz[i] * pa_y[i];

        tr_z_yyyzzz_xxyyzz[i] = 2.0 * tr_z_yzzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyzzz_xxyzz[i] * fe_0 + tr_z_yyzzz_xxyyzz[i] * pa_y[i];

        tr_z_yyyzzz_xxyzzz[i] = 2.0 * tr_z_yzzz_xxyzzz[i] * fe_0 + tr_z_yyzzz_xxzzz[i] * fe_0 + tr_z_yyzzz_xxyzzz[i] * pa_y[i];

        tr_z_yyyzzz_xxzzzz[i] = 2.0 * tr_z_yzzz_xxzzzz[i] * fe_0 + tr_z_yyzzz_xxzzzz[i] * pa_y[i];

        tr_z_yyyzzz_xyyyyy[i] = 2.0 * tr_z_yzzz_xyyyyy[i] * fe_0 + 5.0 * tr_z_yyzzz_xyyyy[i] * fe_0 + tr_z_yyzzz_xyyyyy[i] * pa_y[i];

        tr_z_yyyzzz_xyyyyz[i] = 2.0 * tr_z_yzzz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyzzz_xyyyz[i] * fe_0 + tr_z_yyzzz_xyyyyz[i] * pa_y[i];

        tr_z_yyyzzz_xyyyzz[i] = 2.0 * tr_z_yzzz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyzzz_xyyzz[i] * fe_0 + tr_z_yyzzz_xyyyzz[i] * pa_y[i];

        tr_z_yyyzzz_xyyzzz[i] = 2.0 * tr_z_yzzz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyzzz_xyzzz[i] * fe_0 + tr_z_yyzzz_xyyzzz[i] * pa_y[i];

        tr_z_yyyzzz_xyzzzz[i] = 2.0 * tr_z_yzzz_xyzzzz[i] * fe_0 + tr_z_yyzzz_xzzzz[i] * fe_0 + tr_z_yyzzz_xyzzzz[i] * pa_y[i];

        tr_z_yyyzzz_xzzzzz[i] = 2.0 * tr_z_yzzz_xzzzzz[i] * fe_0 + tr_z_yyzzz_xzzzzz[i] * pa_y[i];

        tr_z_yyyzzz_yyyyyy[i] = 2.0 * tr_z_yzzz_yyyyyy[i] * fe_0 + 6.0 * tr_z_yyzzz_yyyyy[i] * fe_0 + tr_z_yyzzz_yyyyyy[i] * pa_y[i];

        tr_z_yyyzzz_yyyyyz[i] = 2.0 * tr_z_yzzz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyzzz_yyyyz[i] * fe_0 + tr_z_yyzzz_yyyyyz[i] * pa_y[i];

        tr_z_yyyzzz_yyyyzz[i] = 2.0 * tr_z_yzzz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyzzz_yyyzz[i] * fe_0 + tr_z_yyzzz_yyyyzz[i] * pa_y[i];

        tr_z_yyyzzz_yyyzzz[i] = 2.0 * tr_z_yzzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyzzz_yyzzz[i] * fe_0 + tr_z_yyzzz_yyyzzz[i] * pa_y[i];

        tr_z_yyyzzz_yyzzzz[i] = 2.0 * tr_z_yzzz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyzzz_yzzzz[i] * fe_0 + tr_z_yyzzz_yyzzzz[i] * pa_y[i];

        tr_z_yyyzzz_yzzzzz[i] = 2.0 * tr_z_yzzz_yzzzzz[i] * fe_0 + tr_z_yyzzz_zzzzz[i] * fe_0 + tr_z_yyzzz_yzzzzz[i] * pa_y[i];

        tr_z_yyyzzz_zzzzzz[i] = 2.0 * tr_z_yzzz_zzzzzz[i] * fe_0 + tr_z_yyzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 2268-2296 components of targeted buffer : II

    auto tr_z_yyzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2268);

    auto tr_z_yyzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2269);

    auto tr_z_yyzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2270);

    auto tr_z_yyzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2271);

    auto tr_z_yyzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2272);

    auto tr_z_yyzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2273);

    auto tr_z_yyzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2274);

    auto tr_z_yyzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2275);

    auto tr_z_yyzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2276);

    auto tr_z_yyzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2277);

    auto tr_z_yyzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2278);

    auto tr_z_yyzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2279);

    auto tr_z_yyzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2280);

    auto tr_z_yyzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2281);

    auto tr_z_yyzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2282);

    auto tr_z_yyzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2283);

    auto tr_z_yyzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2284);

    auto tr_z_yyzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2285);

    auto tr_z_yyzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2286);

    auto tr_z_yyzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2287);

    auto tr_z_yyzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2288);

    auto tr_z_yyzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2289);

    auto tr_z_yyzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2290);

    auto tr_z_yyzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2291);

    auto tr_z_yyzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2292);

    auto tr_z_yyzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2293);

    auto tr_z_yyzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2294);

    auto tr_z_yyzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2295);

    #pragma omp simd aligned(pa_y, tr_z_yyzzzz_xxxxxx, tr_z_yyzzzz_xxxxxy, tr_z_yyzzzz_xxxxxz, tr_z_yyzzzz_xxxxyy, tr_z_yyzzzz_xxxxyz, tr_z_yyzzzz_xxxxzz, tr_z_yyzzzz_xxxyyy, tr_z_yyzzzz_xxxyyz, tr_z_yyzzzz_xxxyzz, tr_z_yyzzzz_xxxzzz, tr_z_yyzzzz_xxyyyy, tr_z_yyzzzz_xxyyyz, tr_z_yyzzzz_xxyyzz, tr_z_yyzzzz_xxyzzz, tr_z_yyzzzz_xxzzzz, tr_z_yyzzzz_xyyyyy, tr_z_yyzzzz_xyyyyz, tr_z_yyzzzz_xyyyzz, tr_z_yyzzzz_xyyzzz, tr_z_yyzzzz_xyzzzz, tr_z_yyzzzz_xzzzzz, tr_z_yyzzzz_yyyyyy, tr_z_yyzzzz_yyyyyz, tr_z_yyzzzz_yyyyzz, tr_z_yyzzzz_yyyzzz, tr_z_yyzzzz_yyzzzz, tr_z_yyzzzz_yzzzzz, tr_z_yyzzzz_zzzzzz, tr_z_yzzzz_xxxxx, tr_z_yzzzz_xxxxxx, tr_z_yzzzz_xxxxxy, tr_z_yzzzz_xxxxxz, tr_z_yzzzz_xxxxy, tr_z_yzzzz_xxxxyy, tr_z_yzzzz_xxxxyz, tr_z_yzzzz_xxxxz, tr_z_yzzzz_xxxxzz, tr_z_yzzzz_xxxyy, tr_z_yzzzz_xxxyyy, tr_z_yzzzz_xxxyyz, tr_z_yzzzz_xxxyz, tr_z_yzzzz_xxxyzz, tr_z_yzzzz_xxxzz, tr_z_yzzzz_xxxzzz, tr_z_yzzzz_xxyyy, tr_z_yzzzz_xxyyyy, tr_z_yzzzz_xxyyyz, tr_z_yzzzz_xxyyz, tr_z_yzzzz_xxyyzz, tr_z_yzzzz_xxyzz, tr_z_yzzzz_xxyzzz, tr_z_yzzzz_xxzzz, tr_z_yzzzz_xxzzzz, tr_z_yzzzz_xyyyy, tr_z_yzzzz_xyyyyy, tr_z_yzzzz_xyyyyz, tr_z_yzzzz_xyyyz, tr_z_yzzzz_xyyyzz, tr_z_yzzzz_xyyzz, tr_z_yzzzz_xyyzzz, tr_z_yzzzz_xyzzz, tr_z_yzzzz_xyzzzz, tr_z_yzzzz_xzzzz, tr_z_yzzzz_xzzzzz, tr_z_yzzzz_yyyyy, tr_z_yzzzz_yyyyyy, tr_z_yzzzz_yyyyyz, tr_z_yzzzz_yyyyz, tr_z_yzzzz_yyyyzz, tr_z_yzzzz_yyyzz, tr_z_yzzzz_yyyzzz, tr_z_yzzzz_yyzzz, tr_z_yzzzz_yyzzzz, tr_z_yzzzz_yzzzz, tr_z_yzzzz_yzzzzz, tr_z_yzzzz_zzzzz, tr_z_yzzzz_zzzzzz, tr_z_zzzz_xxxxxx, tr_z_zzzz_xxxxxy, tr_z_zzzz_xxxxxz, tr_z_zzzz_xxxxyy, tr_z_zzzz_xxxxyz, tr_z_zzzz_xxxxzz, tr_z_zzzz_xxxyyy, tr_z_zzzz_xxxyyz, tr_z_zzzz_xxxyzz, tr_z_zzzz_xxxzzz, tr_z_zzzz_xxyyyy, tr_z_zzzz_xxyyyz, tr_z_zzzz_xxyyzz, tr_z_zzzz_xxyzzz, tr_z_zzzz_xxzzzz, tr_z_zzzz_xyyyyy, tr_z_zzzz_xyyyyz, tr_z_zzzz_xyyyzz, tr_z_zzzz_xyyzzz, tr_z_zzzz_xyzzzz, tr_z_zzzz_xzzzzz, tr_z_zzzz_yyyyyy, tr_z_zzzz_yyyyyz, tr_z_zzzz_yyyyzz, tr_z_zzzz_yyyzzz, tr_z_zzzz_yyzzzz, tr_z_zzzz_yzzzzz, tr_z_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzzz_xxxxxx[i] = tr_z_zzzz_xxxxxx[i] * fe_0 + tr_z_yzzzz_xxxxxx[i] * pa_y[i];

        tr_z_yyzzzz_xxxxxy[i] = tr_z_zzzz_xxxxxy[i] * fe_0 + tr_z_yzzzz_xxxxx[i] * fe_0 + tr_z_yzzzz_xxxxxy[i] * pa_y[i];

        tr_z_yyzzzz_xxxxxz[i] = tr_z_zzzz_xxxxxz[i] * fe_0 + tr_z_yzzzz_xxxxxz[i] * pa_y[i];

        tr_z_yyzzzz_xxxxyy[i] = tr_z_zzzz_xxxxyy[i] * fe_0 + 2.0 * tr_z_yzzzz_xxxxy[i] * fe_0 + tr_z_yzzzz_xxxxyy[i] * pa_y[i];

        tr_z_yyzzzz_xxxxyz[i] = tr_z_zzzz_xxxxyz[i] * fe_0 + tr_z_yzzzz_xxxxz[i] * fe_0 + tr_z_yzzzz_xxxxyz[i] * pa_y[i];

        tr_z_yyzzzz_xxxxzz[i] = tr_z_zzzz_xxxxzz[i] * fe_0 + tr_z_yzzzz_xxxxzz[i] * pa_y[i];

        tr_z_yyzzzz_xxxyyy[i] = tr_z_zzzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_yzzzz_xxxyy[i] * fe_0 + tr_z_yzzzz_xxxyyy[i] * pa_y[i];

        tr_z_yyzzzz_xxxyyz[i] = tr_z_zzzz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yzzzz_xxxyz[i] * fe_0 + tr_z_yzzzz_xxxyyz[i] * pa_y[i];

        tr_z_yyzzzz_xxxyzz[i] = tr_z_zzzz_xxxyzz[i] * fe_0 + tr_z_yzzzz_xxxzz[i] * fe_0 + tr_z_yzzzz_xxxyzz[i] * pa_y[i];

        tr_z_yyzzzz_xxxzzz[i] = tr_z_zzzz_xxxzzz[i] * fe_0 + tr_z_yzzzz_xxxzzz[i] * pa_y[i];

        tr_z_yyzzzz_xxyyyy[i] = tr_z_zzzz_xxyyyy[i] * fe_0 + 4.0 * tr_z_yzzzz_xxyyy[i] * fe_0 + tr_z_yzzzz_xxyyyy[i] * pa_y[i];

        tr_z_yyzzzz_xxyyyz[i] = tr_z_zzzz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yzzzz_xxyyz[i] * fe_0 + tr_z_yzzzz_xxyyyz[i] * pa_y[i];

        tr_z_yyzzzz_xxyyzz[i] = tr_z_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yzzzz_xxyzz[i] * fe_0 + tr_z_yzzzz_xxyyzz[i] * pa_y[i];

        tr_z_yyzzzz_xxyzzz[i] = tr_z_zzzz_xxyzzz[i] * fe_0 + tr_z_yzzzz_xxzzz[i] * fe_0 + tr_z_yzzzz_xxyzzz[i] * pa_y[i];

        tr_z_yyzzzz_xxzzzz[i] = tr_z_zzzz_xxzzzz[i] * fe_0 + tr_z_yzzzz_xxzzzz[i] * pa_y[i];

        tr_z_yyzzzz_xyyyyy[i] = tr_z_zzzz_xyyyyy[i] * fe_0 + 5.0 * tr_z_yzzzz_xyyyy[i] * fe_0 + tr_z_yzzzz_xyyyyy[i] * pa_y[i];

        tr_z_yyzzzz_xyyyyz[i] = tr_z_zzzz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yzzzz_xyyyz[i] * fe_0 + tr_z_yzzzz_xyyyyz[i] * pa_y[i];

        tr_z_yyzzzz_xyyyzz[i] = tr_z_zzzz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yzzzz_xyyzz[i] * fe_0 + tr_z_yzzzz_xyyyzz[i] * pa_y[i];

        tr_z_yyzzzz_xyyzzz[i] = tr_z_zzzz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yzzzz_xyzzz[i] * fe_0 + tr_z_yzzzz_xyyzzz[i] * pa_y[i];

        tr_z_yyzzzz_xyzzzz[i] = tr_z_zzzz_xyzzzz[i] * fe_0 + tr_z_yzzzz_xzzzz[i] * fe_0 + tr_z_yzzzz_xyzzzz[i] * pa_y[i];

        tr_z_yyzzzz_xzzzzz[i] = tr_z_zzzz_xzzzzz[i] * fe_0 + tr_z_yzzzz_xzzzzz[i] * pa_y[i];

        tr_z_yyzzzz_yyyyyy[i] = tr_z_zzzz_yyyyyy[i] * fe_0 + 6.0 * tr_z_yzzzz_yyyyy[i] * fe_0 + tr_z_yzzzz_yyyyyy[i] * pa_y[i];

        tr_z_yyzzzz_yyyyyz[i] = tr_z_zzzz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yzzzz_yyyyz[i] * fe_0 + tr_z_yzzzz_yyyyyz[i] * pa_y[i];

        tr_z_yyzzzz_yyyyzz[i] = tr_z_zzzz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yzzzz_yyyzz[i] * fe_0 + tr_z_yzzzz_yyyyzz[i] * pa_y[i];

        tr_z_yyzzzz_yyyzzz[i] = tr_z_zzzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yzzzz_yyzzz[i] * fe_0 + tr_z_yzzzz_yyyzzz[i] * pa_y[i];

        tr_z_yyzzzz_yyzzzz[i] = tr_z_zzzz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yzzzz_yzzzz[i] * fe_0 + tr_z_yzzzz_yyzzzz[i] * pa_y[i];

        tr_z_yyzzzz_yzzzzz[i] = tr_z_zzzz_yzzzzz[i] * fe_0 + tr_z_yzzzz_zzzzz[i] * fe_0 + tr_z_yzzzz_yzzzzz[i] * pa_y[i];

        tr_z_yyzzzz_zzzzzz[i] = tr_z_zzzz_zzzzzz[i] * fe_0 + tr_z_yzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 2296-2324 components of targeted buffer : II

    auto tr_z_yzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2296);

    auto tr_z_yzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2297);

    auto tr_z_yzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2298);

    auto tr_z_yzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2299);

    auto tr_z_yzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2300);

    auto tr_z_yzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2301);

    auto tr_z_yzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2302);

    auto tr_z_yzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2303);

    auto tr_z_yzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2304);

    auto tr_z_yzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2305);

    auto tr_z_yzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2306);

    auto tr_z_yzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2307);

    auto tr_z_yzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2308);

    auto tr_z_yzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2309);

    auto tr_z_yzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2310);

    auto tr_z_yzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2311);

    auto tr_z_yzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2312);

    auto tr_z_yzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2313);

    auto tr_z_yzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2314);

    auto tr_z_yzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2315);

    auto tr_z_yzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2316);

    auto tr_z_yzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2317);

    auto tr_z_yzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2318);

    auto tr_z_yzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2319);

    auto tr_z_yzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2320);

    auto tr_z_yzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2321);

    auto tr_z_yzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2322);

    auto tr_z_yzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2323);

    #pragma omp simd aligned(pa_y, tr_z_yzzzzz_xxxxxx, tr_z_yzzzzz_xxxxxy, tr_z_yzzzzz_xxxxxz, tr_z_yzzzzz_xxxxyy, tr_z_yzzzzz_xxxxyz, tr_z_yzzzzz_xxxxzz, tr_z_yzzzzz_xxxyyy, tr_z_yzzzzz_xxxyyz, tr_z_yzzzzz_xxxyzz, tr_z_yzzzzz_xxxzzz, tr_z_yzzzzz_xxyyyy, tr_z_yzzzzz_xxyyyz, tr_z_yzzzzz_xxyyzz, tr_z_yzzzzz_xxyzzz, tr_z_yzzzzz_xxzzzz, tr_z_yzzzzz_xyyyyy, tr_z_yzzzzz_xyyyyz, tr_z_yzzzzz_xyyyzz, tr_z_yzzzzz_xyyzzz, tr_z_yzzzzz_xyzzzz, tr_z_yzzzzz_xzzzzz, tr_z_yzzzzz_yyyyyy, tr_z_yzzzzz_yyyyyz, tr_z_yzzzzz_yyyyzz, tr_z_yzzzzz_yyyzzz, tr_z_yzzzzz_yyzzzz, tr_z_yzzzzz_yzzzzz, tr_z_yzzzzz_zzzzzz, tr_z_zzzzz_xxxxx, tr_z_zzzzz_xxxxxx, tr_z_zzzzz_xxxxxy, tr_z_zzzzz_xxxxxz, tr_z_zzzzz_xxxxy, tr_z_zzzzz_xxxxyy, tr_z_zzzzz_xxxxyz, tr_z_zzzzz_xxxxz, tr_z_zzzzz_xxxxzz, tr_z_zzzzz_xxxyy, tr_z_zzzzz_xxxyyy, tr_z_zzzzz_xxxyyz, tr_z_zzzzz_xxxyz, tr_z_zzzzz_xxxyzz, tr_z_zzzzz_xxxzz, tr_z_zzzzz_xxxzzz, tr_z_zzzzz_xxyyy, tr_z_zzzzz_xxyyyy, tr_z_zzzzz_xxyyyz, tr_z_zzzzz_xxyyz, tr_z_zzzzz_xxyyzz, tr_z_zzzzz_xxyzz, tr_z_zzzzz_xxyzzz, tr_z_zzzzz_xxzzz, tr_z_zzzzz_xxzzzz, tr_z_zzzzz_xyyyy, tr_z_zzzzz_xyyyyy, tr_z_zzzzz_xyyyyz, tr_z_zzzzz_xyyyz, tr_z_zzzzz_xyyyzz, tr_z_zzzzz_xyyzz, tr_z_zzzzz_xyyzzz, tr_z_zzzzz_xyzzz, tr_z_zzzzz_xyzzzz, tr_z_zzzzz_xzzzz, tr_z_zzzzz_xzzzzz, tr_z_zzzzz_yyyyy, tr_z_zzzzz_yyyyyy, tr_z_zzzzz_yyyyyz, tr_z_zzzzz_yyyyz, tr_z_zzzzz_yyyyzz, tr_z_zzzzz_yyyzz, tr_z_zzzzz_yyyzzz, tr_z_zzzzz_yyzzz, tr_z_zzzzz_yyzzzz, tr_z_zzzzz_yzzzz, tr_z_zzzzz_yzzzzz, tr_z_zzzzz_zzzzz, tr_z_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzzz_xxxxxx[i] = tr_z_zzzzz_xxxxxx[i] * pa_y[i];

        tr_z_yzzzzz_xxxxxy[i] = tr_z_zzzzz_xxxxx[i] * fe_0 + tr_z_zzzzz_xxxxxy[i] * pa_y[i];

        tr_z_yzzzzz_xxxxxz[i] = tr_z_zzzzz_xxxxxz[i] * pa_y[i];

        tr_z_yzzzzz_xxxxyy[i] = 2.0 * tr_z_zzzzz_xxxxy[i] * fe_0 + tr_z_zzzzz_xxxxyy[i] * pa_y[i];

        tr_z_yzzzzz_xxxxyz[i] = tr_z_zzzzz_xxxxz[i] * fe_0 + tr_z_zzzzz_xxxxyz[i] * pa_y[i];

        tr_z_yzzzzz_xxxxzz[i] = tr_z_zzzzz_xxxxzz[i] * pa_y[i];

        tr_z_yzzzzz_xxxyyy[i] = 3.0 * tr_z_zzzzz_xxxyy[i] * fe_0 + tr_z_zzzzz_xxxyyy[i] * pa_y[i];

        tr_z_yzzzzz_xxxyyz[i] = 2.0 * tr_z_zzzzz_xxxyz[i] * fe_0 + tr_z_zzzzz_xxxyyz[i] * pa_y[i];

        tr_z_yzzzzz_xxxyzz[i] = tr_z_zzzzz_xxxzz[i] * fe_0 + tr_z_zzzzz_xxxyzz[i] * pa_y[i];

        tr_z_yzzzzz_xxxzzz[i] = tr_z_zzzzz_xxxzzz[i] * pa_y[i];

        tr_z_yzzzzz_xxyyyy[i] = 4.0 * tr_z_zzzzz_xxyyy[i] * fe_0 + tr_z_zzzzz_xxyyyy[i] * pa_y[i];

        tr_z_yzzzzz_xxyyyz[i] = 3.0 * tr_z_zzzzz_xxyyz[i] * fe_0 + tr_z_zzzzz_xxyyyz[i] * pa_y[i];

        tr_z_yzzzzz_xxyyzz[i] = 2.0 * tr_z_zzzzz_xxyzz[i] * fe_0 + tr_z_zzzzz_xxyyzz[i] * pa_y[i];

        tr_z_yzzzzz_xxyzzz[i] = tr_z_zzzzz_xxzzz[i] * fe_0 + tr_z_zzzzz_xxyzzz[i] * pa_y[i];

        tr_z_yzzzzz_xxzzzz[i] = tr_z_zzzzz_xxzzzz[i] * pa_y[i];

        tr_z_yzzzzz_xyyyyy[i] = 5.0 * tr_z_zzzzz_xyyyy[i] * fe_0 + tr_z_zzzzz_xyyyyy[i] * pa_y[i];

        tr_z_yzzzzz_xyyyyz[i] = 4.0 * tr_z_zzzzz_xyyyz[i] * fe_0 + tr_z_zzzzz_xyyyyz[i] * pa_y[i];

        tr_z_yzzzzz_xyyyzz[i] = 3.0 * tr_z_zzzzz_xyyzz[i] * fe_0 + tr_z_zzzzz_xyyyzz[i] * pa_y[i];

        tr_z_yzzzzz_xyyzzz[i] = 2.0 * tr_z_zzzzz_xyzzz[i] * fe_0 + tr_z_zzzzz_xyyzzz[i] * pa_y[i];

        tr_z_yzzzzz_xyzzzz[i] = tr_z_zzzzz_xzzzz[i] * fe_0 + tr_z_zzzzz_xyzzzz[i] * pa_y[i];

        tr_z_yzzzzz_xzzzzz[i] = tr_z_zzzzz_xzzzzz[i] * pa_y[i];

        tr_z_yzzzzz_yyyyyy[i] = 6.0 * tr_z_zzzzz_yyyyy[i] * fe_0 + tr_z_zzzzz_yyyyyy[i] * pa_y[i];

        tr_z_yzzzzz_yyyyyz[i] = 5.0 * tr_z_zzzzz_yyyyz[i] * fe_0 + tr_z_zzzzz_yyyyyz[i] * pa_y[i];

        tr_z_yzzzzz_yyyyzz[i] = 4.0 * tr_z_zzzzz_yyyzz[i] * fe_0 + tr_z_zzzzz_yyyyzz[i] * pa_y[i];

        tr_z_yzzzzz_yyyzzz[i] = 3.0 * tr_z_zzzzz_yyzzz[i] * fe_0 + tr_z_zzzzz_yyyzzz[i] * pa_y[i];

        tr_z_yzzzzz_yyzzzz[i] = 2.0 * tr_z_zzzzz_yzzzz[i] * fe_0 + tr_z_zzzzz_yyzzzz[i] * pa_y[i];

        tr_z_yzzzzz_yzzzzz[i] = tr_z_zzzzz_zzzzz[i] * fe_0 + tr_z_zzzzz_yzzzzz[i] * pa_y[i];

        tr_z_yzzzzz_zzzzzz[i] = tr_z_zzzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 2324-2352 components of targeted buffer : II

    auto tr_z_zzzzzz_xxxxxx = pbuffer.data(idx_dip_ii + 2324);

    auto tr_z_zzzzzz_xxxxxy = pbuffer.data(idx_dip_ii + 2325);

    auto tr_z_zzzzzz_xxxxxz = pbuffer.data(idx_dip_ii + 2326);

    auto tr_z_zzzzzz_xxxxyy = pbuffer.data(idx_dip_ii + 2327);

    auto tr_z_zzzzzz_xxxxyz = pbuffer.data(idx_dip_ii + 2328);

    auto tr_z_zzzzzz_xxxxzz = pbuffer.data(idx_dip_ii + 2329);

    auto tr_z_zzzzzz_xxxyyy = pbuffer.data(idx_dip_ii + 2330);

    auto tr_z_zzzzzz_xxxyyz = pbuffer.data(idx_dip_ii + 2331);

    auto tr_z_zzzzzz_xxxyzz = pbuffer.data(idx_dip_ii + 2332);

    auto tr_z_zzzzzz_xxxzzz = pbuffer.data(idx_dip_ii + 2333);

    auto tr_z_zzzzzz_xxyyyy = pbuffer.data(idx_dip_ii + 2334);

    auto tr_z_zzzzzz_xxyyyz = pbuffer.data(idx_dip_ii + 2335);

    auto tr_z_zzzzzz_xxyyzz = pbuffer.data(idx_dip_ii + 2336);

    auto tr_z_zzzzzz_xxyzzz = pbuffer.data(idx_dip_ii + 2337);

    auto tr_z_zzzzzz_xxzzzz = pbuffer.data(idx_dip_ii + 2338);

    auto tr_z_zzzzzz_xyyyyy = pbuffer.data(idx_dip_ii + 2339);

    auto tr_z_zzzzzz_xyyyyz = pbuffer.data(idx_dip_ii + 2340);

    auto tr_z_zzzzzz_xyyyzz = pbuffer.data(idx_dip_ii + 2341);

    auto tr_z_zzzzzz_xyyzzz = pbuffer.data(idx_dip_ii + 2342);

    auto tr_z_zzzzzz_xyzzzz = pbuffer.data(idx_dip_ii + 2343);

    auto tr_z_zzzzzz_xzzzzz = pbuffer.data(idx_dip_ii + 2344);

    auto tr_z_zzzzzz_yyyyyy = pbuffer.data(idx_dip_ii + 2345);

    auto tr_z_zzzzzz_yyyyyz = pbuffer.data(idx_dip_ii + 2346);

    auto tr_z_zzzzzz_yyyyzz = pbuffer.data(idx_dip_ii + 2347);

    auto tr_z_zzzzzz_yyyzzz = pbuffer.data(idx_dip_ii + 2348);

    auto tr_z_zzzzzz_yyzzzz = pbuffer.data(idx_dip_ii + 2349);

    auto tr_z_zzzzzz_yzzzzz = pbuffer.data(idx_dip_ii + 2350);

    auto tr_z_zzzzzz_zzzzzz = pbuffer.data(idx_dip_ii + 2351);

    #pragma omp simd aligned(pa_z, tr_z_zzzz_xxxxxx, tr_z_zzzz_xxxxxy, tr_z_zzzz_xxxxxz, tr_z_zzzz_xxxxyy, tr_z_zzzz_xxxxyz, tr_z_zzzz_xxxxzz, tr_z_zzzz_xxxyyy, tr_z_zzzz_xxxyyz, tr_z_zzzz_xxxyzz, tr_z_zzzz_xxxzzz, tr_z_zzzz_xxyyyy, tr_z_zzzz_xxyyyz, tr_z_zzzz_xxyyzz, tr_z_zzzz_xxyzzz, tr_z_zzzz_xxzzzz, tr_z_zzzz_xyyyyy, tr_z_zzzz_xyyyyz, tr_z_zzzz_xyyyzz, tr_z_zzzz_xyyzzz, tr_z_zzzz_xyzzzz, tr_z_zzzz_xzzzzz, tr_z_zzzz_yyyyyy, tr_z_zzzz_yyyyyz, tr_z_zzzz_yyyyzz, tr_z_zzzz_yyyzzz, tr_z_zzzz_yyzzzz, tr_z_zzzz_yzzzzz, tr_z_zzzz_zzzzzz, tr_z_zzzzz_xxxxx, tr_z_zzzzz_xxxxxx, tr_z_zzzzz_xxxxxy, tr_z_zzzzz_xxxxxz, tr_z_zzzzz_xxxxy, tr_z_zzzzz_xxxxyy, tr_z_zzzzz_xxxxyz, tr_z_zzzzz_xxxxz, tr_z_zzzzz_xxxxzz, tr_z_zzzzz_xxxyy, tr_z_zzzzz_xxxyyy, tr_z_zzzzz_xxxyyz, tr_z_zzzzz_xxxyz, tr_z_zzzzz_xxxyzz, tr_z_zzzzz_xxxzz, tr_z_zzzzz_xxxzzz, tr_z_zzzzz_xxyyy, tr_z_zzzzz_xxyyyy, tr_z_zzzzz_xxyyyz, tr_z_zzzzz_xxyyz, tr_z_zzzzz_xxyyzz, tr_z_zzzzz_xxyzz, tr_z_zzzzz_xxyzzz, tr_z_zzzzz_xxzzz, tr_z_zzzzz_xxzzzz, tr_z_zzzzz_xyyyy, tr_z_zzzzz_xyyyyy, tr_z_zzzzz_xyyyyz, tr_z_zzzzz_xyyyz, tr_z_zzzzz_xyyyzz, tr_z_zzzzz_xyyzz, tr_z_zzzzz_xyyzzz, tr_z_zzzzz_xyzzz, tr_z_zzzzz_xyzzzz, tr_z_zzzzz_xzzzz, tr_z_zzzzz_xzzzzz, tr_z_zzzzz_yyyyy, tr_z_zzzzz_yyyyyy, tr_z_zzzzz_yyyyyz, tr_z_zzzzz_yyyyz, tr_z_zzzzz_yyyyzz, tr_z_zzzzz_yyyzz, tr_z_zzzzz_yyyzzz, tr_z_zzzzz_yyzzz, tr_z_zzzzz_yyzzzz, tr_z_zzzzz_yzzzz, tr_z_zzzzz_yzzzzz, tr_z_zzzzz_zzzzz, tr_z_zzzzz_zzzzzz, tr_z_zzzzzz_xxxxxx, tr_z_zzzzzz_xxxxxy, tr_z_zzzzzz_xxxxxz, tr_z_zzzzzz_xxxxyy, tr_z_zzzzzz_xxxxyz, tr_z_zzzzzz_xxxxzz, tr_z_zzzzzz_xxxyyy, tr_z_zzzzzz_xxxyyz, tr_z_zzzzzz_xxxyzz, tr_z_zzzzzz_xxxzzz, tr_z_zzzzzz_xxyyyy, tr_z_zzzzzz_xxyyyz, tr_z_zzzzzz_xxyyzz, tr_z_zzzzzz_xxyzzz, tr_z_zzzzzz_xxzzzz, tr_z_zzzzzz_xyyyyy, tr_z_zzzzzz_xyyyyz, tr_z_zzzzzz_xyyyzz, tr_z_zzzzzz_xyyzzz, tr_z_zzzzzz_xyzzzz, tr_z_zzzzzz_xzzzzz, tr_z_zzzzzz_yyyyyy, tr_z_zzzzzz_yyyyyz, tr_z_zzzzzz_yyyyzz, tr_z_zzzzzz_yyyzzz, tr_z_zzzzzz_yyzzzz, tr_z_zzzzzz_yzzzzz, tr_z_zzzzzz_zzzzzz, ts_zzzzz_xxxxxx, ts_zzzzz_xxxxxy, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxyy, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyyy, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyyy, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyyy, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyyy, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzzz_xxxxxx[i] = 5.0 * tr_z_zzzz_xxxxxx[i] * fe_0 + ts_zzzzz_xxxxxx[i] * fe_0 + tr_z_zzzzz_xxxxxx[i] * pa_z[i];

        tr_z_zzzzzz_xxxxxy[i] = 5.0 * tr_z_zzzz_xxxxxy[i] * fe_0 + ts_zzzzz_xxxxxy[i] * fe_0 + tr_z_zzzzz_xxxxxy[i] * pa_z[i];

        tr_z_zzzzzz_xxxxxz[i] = 5.0 * tr_z_zzzz_xxxxxz[i] * fe_0 + tr_z_zzzzz_xxxxx[i] * fe_0 + ts_zzzzz_xxxxxz[i] * fe_0 + tr_z_zzzzz_xxxxxz[i] * pa_z[i];

        tr_z_zzzzzz_xxxxyy[i] = 5.0 * tr_z_zzzz_xxxxyy[i] * fe_0 + ts_zzzzz_xxxxyy[i] * fe_0 + tr_z_zzzzz_xxxxyy[i] * pa_z[i];

        tr_z_zzzzzz_xxxxyz[i] = 5.0 * tr_z_zzzz_xxxxyz[i] * fe_0 + tr_z_zzzzz_xxxxy[i] * fe_0 + ts_zzzzz_xxxxyz[i] * fe_0 + tr_z_zzzzz_xxxxyz[i] * pa_z[i];

        tr_z_zzzzzz_xxxxzz[i] = 5.0 * tr_z_zzzz_xxxxzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxxzz[i] * fe_0 + tr_z_zzzzz_xxxxzz[i] * pa_z[i];

        tr_z_zzzzzz_xxxyyy[i] = 5.0 * tr_z_zzzz_xxxyyy[i] * fe_0 + ts_zzzzz_xxxyyy[i] * fe_0 + tr_z_zzzzz_xxxyyy[i] * pa_z[i];

        tr_z_zzzzzz_xxxyyz[i] = 5.0 * tr_z_zzzz_xxxyyz[i] * fe_0 + tr_z_zzzzz_xxxyy[i] * fe_0 + ts_zzzzz_xxxyyz[i] * fe_0 + tr_z_zzzzz_xxxyyz[i] * pa_z[i];

        tr_z_zzzzzz_xxxyzz[i] = 5.0 * tr_z_zzzz_xxxyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxyzz[i] * fe_0 + tr_z_zzzzz_xxxyzz[i] * pa_z[i];

        tr_z_zzzzzz_xxxzzz[i] = 5.0 * tr_z_zzzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_xxxzz[i] * fe_0 + ts_zzzzz_xxxzzz[i] * fe_0 + tr_z_zzzzz_xxxzzz[i] * pa_z[i];

        tr_z_zzzzzz_xxyyyy[i] = 5.0 * tr_z_zzzz_xxyyyy[i] * fe_0 + ts_zzzzz_xxyyyy[i] * fe_0 + tr_z_zzzzz_xxyyyy[i] * pa_z[i];

        tr_z_zzzzzz_xxyyyz[i] = 5.0 * tr_z_zzzz_xxyyyz[i] * fe_0 + tr_z_zzzzz_xxyyy[i] * fe_0 + ts_zzzzz_xxyyyz[i] * fe_0 + tr_z_zzzzz_xxyyyz[i] * pa_z[i];

        tr_z_zzzzzz_xxyyzz[i] = 5.0 * tr_z_zzzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxyyzz[i] * fe_0 + tr_z_zzzzz_xxyyzz[i] * pa_z[i];

        tr_z_zzzzzz_xxyzzz[i] = 5.0 * tr_z_zzzz_xxyzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_xxyzz[i] * fe_0 + ts_zzzzz_xxyzzz[i] * fe_0 + tr_z_zzzzz_xxyzzz[i] * pa_z[i];

        tr_z_zzzzzz_xxzzzz[i] = 5.0 * tr_z_zzzz_xxzzzz[i] * fe_0 + 4.0 * tr_z_zzzzz_xxzzz[i] * fe_0 + ts_zzzzz_xxzzzz[i] * fe_0 + tr_z_zzzzz_xxzzzz[i] * pa_z[i];

        tr_z_zzzzzz_xyyyyy[i] = 5.0 * tr_z_zzzz_xyyyyy[i] * fe_0 + ts_zzzzz_xyyyyy[i] * fe_0 + tr_z_zzzzz_xyyyyy[i] * pa_z[i];

        tr_z_zzzzzz_xyyyyz[i] = 5.0 * tr_z_zzzz_xyyyyz[i] * fe_0 + tr_z_zzzzz_xyyyy[i] * fe_0 + ts_zzzzz_xyyyyz[i] * fe_0 + tr_z_zzzzz_xyyyyz[i] * pa_z[i];

        tr_z_zzzzzz_xyyyzz[i] = 5.0 * tr_z_zzzz_xyyyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xyyyz[i] * fe_0 + ts_zzzzz_xyyyzz[i] * fe_0 + tr_z_zzzzz_xyyyzz[i] * pa_z[i];

        tr_z_zzzzzz_xyyzzz[i] = 5.0 * tr_z_zzzz_xyyzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_xyyzz[i] * fe_0 + ts_zzzzz_xyyzzz[i] * fe_0 + tr_z_zzzzz_xyyzzz[i] * pa_z[i];

        tr_z_zzzzzz_xyzzzz[i] = 5.0 * tr_z_zzzz_xyzzzz[i] * fe_0 + 4.0 * tr_z_zzzzz_xyzzz[i] * fe_0 + ts_zzzzz_xyzzzz[i] * fe_0 + tr_z_zzzzz_xyzzzz[i] * pa_z[i];

        tr_z_zzzzzz_xzzzzz[i] = 5.0 * tr_z_zzzz_xzzzzz[i] * fe_0 + 5.0 * tr_z_zzzzz_xzzzz[i] * fe_0 + ts_zzzzz_xzzzzz[i] * fe_0 + tr_z_zzzzz_xzzzzz[i] * pa_z[i];

        tr_z_zzzzzz_yyyyyy[i] = 5.0 * tr_z_zzzz_yyyyyy[i] * fe_0 + ts_zzzzz_yyyyyy[i] * fe_0 + tr_z_zzzzz_yyyyyy[i] * pa_z[i];

        tr_z_zzzzzz_yyyyyz[i] = 5.0 * tr_z_zzzz_yyyyyz[i] * fe_0 + tr_z_zzzzz_yyyyy[i] * fe_0 + ts_zzzzz_yyyyyz[i] * fe_0 + tr_z_zzzzz_yyyyyz[i] * pa_z[i];

        tr_z_zzzzzz_yyyyzz[i] = 5.0 * tr_z_zzzz_yyyyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_yyyyz[i] * fe_0 + ts_zzzzz_yyyyzz[i] * fe_0 + tr_z_zzzzz_yyyyzz[i] * pa_z[i];

        tr_z_zzzzzz_yyyzzz[i] = 5.0 * tr_z_zzzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_yyyzz[i] * fe_0 + ts_zzzzz_yyyzzz[i] * fe_0 + tr_z_zzzzz_yyyzzz[i] * pa_z[i];

        tr_z_zzzzzz_yyzzzz[i] = 5.0 * tr_z_zzzz_yyzzzz[i] * fe_0 + 4.0 * tr_z_zzzzz_yyzzz[i] * fe_0 + ts_zzzzz_yyzzzz[i] * fe_0 + tr_z_zzzzz_yyzzzz[i] * pa_z[i];

        tr_z_zzzzzz_yzzzzz[i] = 5.0 * tr_z_zzzz_yzzzzz[i] * fe_0 + 5.0 * tr_z_zzzzz_yzzzz[i] * fe_0 + ts_zzzzz_yzzzzz[i] * fe_0 + tr_z_zzzzz_yzzzzz[i] * pa_z[i];

        tr_z_zzzzzz_zzzzzz[i] = 5.0 * tr_z_zzzz_zzzzzz[i] * fe_0 + 6.0 * tr_z_zzzzz_zzzzz[i] * fe_0 + ts_zzzzz_zzzzzz[i] * fe_0 + tr_z_zzzzz_zzzzzz[i] * pa_z[i];
    }

}

} // diprec namespace

