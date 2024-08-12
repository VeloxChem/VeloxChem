#include "NuclearPotentialPrimRecII.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_ii(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_ii,
                               const size_t idx_npot_0_gi,
                               const size_t idx_npot_1_gi,
                               const size_t idx_npot_0_hh,
                               const size_t idx_npot_1_hh,
                               const size_t idx_npot_0_hi,
                               const size_t idx_npot_1_hi,
                               const CSimdArray<double>& factors,
                               const size_t idx_rpa,
                               const size_t idx_rpc,
                               const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : GI

    auto ta_xxxx_xxxxxx_0 = pbuffer.data(idx_npot_0_gi);

    auto ta_xxxx_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 1);

    auto ta_xxxx_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 2);

    auto ta_xxxx_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 3);

    auto ta_xxxx_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 4);

    auto ta_xxxx_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 5);

    auto ta_xxxx_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 6);

    auto ta_xxxx_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 7);

    auto ta_xxxx_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 8);

    auto ta_xxxx_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 9);

    auto ta_xxxx_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 10);

    auto ta_xxxx_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 11);

    auto ta_xxxx_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 12);

    auto ta_xxxx_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 13);

    auto ta_xxxx_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 14);

    auto ta_xxxx_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 15);

    auto ta_xxxx_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 16);

    auto ta_xxxx_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 17);

    auto ta_xxxx_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 18);

    auto ta_xxxx_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 19);

    auto ta_xxxx_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 20);

    auto ta_xxxx_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 21);

    auto ta_xxxx_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 22);

    auto ta_xxxx_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 23);

    auto ta_xxxx_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 24);

    auto ta_xxxx_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 25);

    auto ta_xxxx_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 26);

    auto ta_xxxx_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 27);

    auto ta_xxxy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 28);

    auto ta_xxxy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 30);

    auto ta_xxxy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 33);

    auto ta_xxxy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 37);

    auto ta_xxxy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 42);

    auto ta_xxxy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 48);

    auto ta_xxxy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 49);

    auto ta_xxxy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 50);

    auto ta_xxxy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 51);

    auto ta_xxxy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 52);

    auto ta_xxxy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 53);

    auto ta_xxxy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 54);

    auto ta_xxxz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 56);

    auto ta_xxxz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 57);

    auto ta_xxxz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 58);

    auto ta_xxxz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 59);

    auto ta_xxxz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 61);

    auto ta_xxxz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 62);

    auto ta_xxxz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 65);

    auto ta_xxxz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 66);

    auto ta_xxxz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 70);

    auto ta_xxxz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 71);

    auto ta_xxxz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 76);

    auto ta_xxxz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 78);

    auto ta_xxxz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 79);

    auto ta_xxxz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 80);

    auto ta_xxxz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 81);

    auto ta_xxxz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 82);

    auto ta_xxxz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 83);

    auto ta_xxyy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 84);

    auto ta_xxyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 85);

    auto ta_xxyy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 86);

    auto ta_xxyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 87);

    auto ta_xxyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 88);

    auto ta_xxyy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 89);

    auto ta_xxyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 90);

    auto ta_xxyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 91);

    auto ta_xxyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 92);

    auto ta_xxyy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 93);

    auto ta_xxyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 94);

    auto ta_xxyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 95);

    auto ta_xxyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 96);

    auto ta_xxyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 97);

    auto ta_xxyy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 98);

    auto ta_xxyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 99);

    auto ta_xxyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 100);

    auto ta_xxyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 101);

    auto ta_xxyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 102);

    auto ta_xxyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 103);

    auto ta_xxyy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 104);

    auto ta_xxyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 105);

    auto ta_xxyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 106);

    auto ta_xxyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 107);

    auto ta_xxyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 108);

    auto ta_xxyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 109);

    auto ta_xxyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 110);

    auto ta_xxyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 111);

    auto ta_xxyz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 114);

    auto ta_xxyz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 117);

    auto ta_xxyz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 121);

    auto ta_xxyz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 126);

    auto ta_xxyz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 132);

    auto ta_xxyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 134);

    auto ta_xxyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 135);

    auto ta_xxyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 136);

    auto ta_xxyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 137);

    auto ta_xxyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 138);

    auto ta_xxzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 140);

    auto ta_xxzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 141);

    auto ta_xxzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 142);

    auto ta_xxzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 143);

    auto ta_xxzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 144);

    auto ta_xxzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 145);

    auto ta_xxzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 146);

    auto ta_xxzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 147);

    auto ta_xxzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 148);

    auto ta_xxzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 149);

    auto ta_xxzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 150);

    auto ta_xxzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 151);

    auto ta_xxzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 152);

    auto ta_xxzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 153);

    auto ta_xxzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 154);

    auto ta_xxzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 155);

    auto ta_xxzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 156);

    auto ta_xxzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 157);

    auto ta_xxzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 158);

    auto ta_xxzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 159);

    auto ta_xxzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 160);

    auto ta_xxzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 161);

    auto ta_xxzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 162);

    auto ta_xxzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 163);

    auto ta_xxzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 164);

    auto ta_xxzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 165);

    auto ta_xxzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 166);

    auto ta_xxzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 167);

    auto ta_xyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 169);

    auto ta_xyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 171);

    auto ta_xyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 172);

    auto ta_xyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 174);

    auto ta_xyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 175);

    auto ta_xyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 176);

    auto ta_xyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 178);

    auto ta_xyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 179);

    auto ta_xyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 180);

    auto ta_xyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 181);

    auto ta_xyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 183);

    auto ta_xyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 184);

    auto ta_xyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 185);

    auto ta_xyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 186);

    auto ta_xyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 187);

    auto ta_xyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 189);

    auto ta_xyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 190);

    auto ta_xyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 191);

    auto ta_xyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 192);

    auto ta_xyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 193);

    auto ta_xyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 194);

    auto ta_xyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 195);

    auto ta_xyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 218);

    auto ta_xyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 219);

    auto ta_xyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 220);

    auto ta_xyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 221);

    auto ta_xyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 222);

    auto ta_xyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 223);

    auto ta_xyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 245);

    auto ta_xyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 246);

    auto ta_xyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 247);

    auto ta_xyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 248);

    auto ta_xyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 249);

    auto ta_xyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 250);

    auto ta_xzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 254);

    auto ta_xzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 256);

    auto ta_xzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 257);

    auto ta_xzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 259);

    auto ta_xzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 260);

    auto ta_xzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 261);

    auto ta_xzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 263);

    auto ta_xzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 264);

    auto ta_xzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 265);

    auto ta_xzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 266);

    auto ta_xzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 268);

    auto ta_xzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 269);

    auto ta_xzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 270);

    auto ta_xzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 271);

    auto ta_xzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 272);

    auto ta_xzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 273);

    auto ta_xzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 274);

    auto ta_xzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 275);

    auto ta_xzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 276);

    auto ta_xzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 277);

    auto ta_xzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 278);

    auto ta_xzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 279);

    auto ta_yyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 280);

    auto ta_yyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 281);

    auto ta_yyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 282);

    auto ta_yyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 283);

    auto ta_yyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 284);

    auto ta_yyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 285);

    auto ta_yyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 286);

    auto ta_yyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 287);

    auto ta_yyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 288);

    auto ta_yyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 289);

    auto ta_yyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 290);

    auto ta_yyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 291);

    auto ta_yyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 292);

    auto ta_yyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 293);

    auto ta_yyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 294);

    auto ta_yyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 295);

    auto ta_yyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 296);

    auto ta_yyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 297);

    auto ta_yyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 298);

    auto ta_yyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 299);

    auto ta_yyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 300);

    auto ta_yyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 301);

    auto ta_yyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 302);

    auto ta_yyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 303);

    auto ta_yyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 304);

    auto ta_yyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 305);

    auto ta_yyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 306);

    auto ta_yyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 307);

    auto ta_yyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 309);

    auto ta_yyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 310);

    auto ta_yyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 311);

    auto ta_yyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 313);

    auto ta_yyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 314);

    auto ta_yyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 317);

    auto ta_yyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 318);

    auto ta_yyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 322);

    auto ta_yyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 323);

    auto ta_yyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 328);

    auto ta_yyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 329);

    auto ta_yyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 330);

    auto ta_yyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 331);

    auto ta_yyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 332);

    auto ta_yyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 333);

    auto ta_yyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 334);

    auto ta_yyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 335);

    auto ta_yyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 336);

    auto ta_yyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 337);

    auto ta_yyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 338);

    auto ta_yyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 339);

    auto ta_yyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 340);

    auto ta_yyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 341);

    auto ta_yyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 342);

    auto ta_yyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 343);

    auto ta_yyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 344);

    auto ta_yyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 345);

    auto ta_yyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 346);

    auto ta_yyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 347);

    auto ta_yyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 348);

    auto ta_yyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 349);

    auto ta_yyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 350);

    auto ta_yyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 351);

    auto ta_yyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 352);

    auto ta_yyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 353);

    auto ta_yyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 354);

    auto ta_yyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 355);

    auto ta_yyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 356);

    auto ta_yyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 357);

    auto ta_yyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 358);

    auto ta_yyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 359);

    auto ta_yyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 360);

    auto ta_yyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 361);

    auto ta_yyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 362);

    auto ta_yyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 363);

    auto ta_yzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 364);

    auto ta_yzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 366);

    auto ta_yzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 368);

    auto ta_yzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 369);

    auto ta_yzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 371);

    auto ta_yzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 372);

    auto ta_yzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 373);

    auto ta_yzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 375);

    auto ta_yzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 376);

    auto ta_yzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 377);

    auto ta_yzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 378);

    auto ta_yzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 380);

    auto ta_yzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 381);

    auto ta_yzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 382);

    auto ta_yzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 383);

    auto ta_yzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 384);

    auto ta_yzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 385);

    auto ta_yzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 386);

    auto ta_yzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 387);

    auto ta_yzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 388);

    auto ta_yzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 389);

    auto ta_yzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 390);

    auto ta_yzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 391);

    auto ta_zzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 392);

    auto ta_zzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 393);

    auto ta_zzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 394);

    auto ta_zzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 395);

    auto ta_zzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 396);

    auto ta_zzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 397);

    auto ta_zzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 398);

    auto ta_zzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 399);

    auto ta_zzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 400);

    auto ta_zzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 401);

    auto ta_zzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 402);

    auto ta_zzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 403);

    auto ta_zzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 404);

    auto ta_zzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 405);

    auto ta_zzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 406);

    auto ta_zzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 407);

    auto ta_zzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 408);

    auto ta_zzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 409);

    auto ta_zzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 410);

    auto ta_zzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 411);

    auto ta_zzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 412);

    auto ta_zzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 413);

    auto ta_zzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 414);

    auto ta_zzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 415);

    auto ta_zzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 416);

    auto ta_zzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 417);

    auto ta_zzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 418);

    auto ta_zzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 419);

    // Set up components of auxiliary buffer : GI

    auto ta_xxxx_xxxxxx_1 = pbuffer.data(idx_npot_1_gi);

    auto ta_xxxx_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 1);

    auto ta_xxxx_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 2);

    auto ta_xxxx_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 3);

    auto ta_xxxx_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 4);

    auto ta_xxxx_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 5);

    auto ta_xxxx_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 6);

    auto ta_xxxx_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 7);

    auto ta_xxxx_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 8);

    auto ta_xxxx_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 9);

    auto ta_xxxx_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 10);

    auto ta_xxxx_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 11);

    auto ta_xxxx_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 12);

    auto ta_xxxx_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 13);

    auto ta_xxxx_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 14);

    auto ta_xxxx_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 15);

    auto ta_xxxx_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 16);

    auto ta_xxxx_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 17);

    auto ta_xxxx_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 18);

    auto ta_xxxx_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 19);

    auto ta_xxxx_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 20);

    auto ta_xxxx_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 21);

    auto ta_xxxx_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 22);

    auto ta_xxxx_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 23);

    auto ta_xxxx_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 24);

    auto ta_xxxx_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 25);

    auto ta_xxxx_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 26);

    auto ta_xxxx_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 27);

    auto ta_xxxy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 28);

    auto ta_xxxy_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 30);

    auto ta_xxxy_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 33);

    auto ta_xxxy_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 37);

    auto ta_xxxy_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 42);

    auto ta_xxxy_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 48);

    auto ta_xxxy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 49);

    auto ta_xxxy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 50);

    auto ta_xxxy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 51);

    auto ta_xxxy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 52);

    auto ta_xxxy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 53);

    auto ta_xxxy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 54);

    auto ta_xxxz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 56);

    auto ta_xxxz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 57);

    auto ta_xxxz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 58);

    auto ta_xxxz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 59);

    auto ta_xxxz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 61);

    auto ta_xxxz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 62);

    auto ta_xxxz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 65);

    auto ta_xxxz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 66);

    auto ta_xxxz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 70);

    auto ta_xxxz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 71);

    auto ta_xxxz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 76);

    auto ta_xxxz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 78);

    auto ta_xxxz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 79);

    auto ta_xxxz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 80);

    auto ta_xxxz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 81);

    auto ta_xxxz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 82);

    auto ta_xxxz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 83);

    auto ta_xxyy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 84);

    auto ta_xxyy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 85);

    auto ta_xxyy_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 86);

    auto ta_xxyy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 87);

    auto ta_xxyy_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 88);

    auto ta_xxyy_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 89);

    auto ta_xxyy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 90);

    auto ta_xxyy_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 91);

    auto ta_xxyy_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 92);

    auto ta_xxyy_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 93);

    auto ta_xxyy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 94);

    auto ta_xxyy_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 95);

    auto ta_xxyy_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 96);

    auto ta_xxyy_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 97);

    auto ta_xxyy_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 98);

    auto ta_xxyy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 99);

    auto ta_xxyy_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 100);

    auto ta_xxyy_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 101);

    auto ta_xxyy_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 102);

    auto ta_xxyy_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 103);

    auto ta_xxyy_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 104);

    auto ta_xxyy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 105);

    auto ta_xxyy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 106);

    auto ta_xxyy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 107);

    auto ta_xxyy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 108);

    auto ta_xxyy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 109);

    auto ta_xxyy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 110);

    auto ta_xxyy_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 111);

    auto ta_xxyz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 114);

    auto ta_xxyz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 117);

    auto ta_xxyz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 121);

    auto ta_xxyz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 126);

    auto ta_xxyz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 132);

    auto ta_xxyz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 134);

    auto ta_xxyz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 135);

    auto ta_xxyz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 136);

    auto ta_xxyz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 137);

    auto ta_xxyz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 138);

    auto ta_xxzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 140);

    auto ta_xxzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 141);

    auto ta_xxzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 142);

    auto ta_xxzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 143);

    auto ta_xxzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 144);

    auto ta_xxzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 145);

    auto ta_xxzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 146);

    auto ta_xxzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 147);

    auto ta_xxzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 148);

    auto ta_xxzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 149);

    auto ta_xxzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 150);

    auto ta_xxzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 151);

    auto ta_xxzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 152);

    auto ta_xxzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 153);

    auto ta_xxzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 154);

    auto ta_xxzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 155);

    auto ta_xxzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 156);

    auto ta_xxzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 157);

    auto ta_xxzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 158);

    auto ta_xxzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 159);

    auto ta_xxzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 160);

    auto ta_xxzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 161);

    auto ta_xxzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 162);

    auto ta_xxzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 163);

    auto ta_xxzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 164);

    auto ta_xxzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 165);

    auto ta_xxzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 166);

    auto ta_xxzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 167);

    auto ta_xyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 169);

    auto ta_xyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 171);

    auto ta_xyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 172);

    auto ta_xyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 174);

    auto ta_xyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 175);

    auto ta_xyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 176);

    auto ta_xyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 178);

    auto ta_xyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 179);

    auto ta_xyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 180);

    auto ta_xyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 181);

    auto ta_xyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 183);

    auto ta_xyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 184);

    auto ta_xyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 185);

    auto ta_xyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 186);

    auto ta_xyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 187);

    auto ta_xyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 189);

    auto ta_xyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 190);

    auto ta_xyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 191);

    auto ta_xyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 192);

    auto ta_xyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 193);

    auto ta_xyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 194);

    auto ta_xyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 195);

    auto ta_xyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 218);

    auto ta_xyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 219);

    auto ta_xyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 220);

    auto ta_xyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 221);

    auto ta_xyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 222);

    auto ta_xyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 223);

    auto ta_xyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 245);

    auto ta_xyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 246);

    auto ta_xyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 247);

    auto ta_xyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 248);

    auto ta_xyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 249);

    auto ta_xyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 250);

    auto ta_xzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 254);

    auto ta_xzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 256);

    auto ta_xzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 257);

    auto ta_xzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 259);

    auto ta_xzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 260);

    auto ta_xzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 261);

    auto ta_xzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 263);

    auto ta_xzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 264);

    auto ta_xzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 265);

    auto ta_xzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 266);

    auto ta_xzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 268);

    auto ta_xzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 269);

    auto ta_xzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 270);

    auto ta_xzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 271);

    auto ta_xzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 272);

    auto ta_xzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 273);

    auto ta_xzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 274);

    auto ta_xzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 275);

    auto ta_xzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 276);

    auto ta_xzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 277);

    auto ta_xzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 278);

    auto ta_xzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 279);

    auto ta_yyyy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 280);

    auto ta_yyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 281);

    auto ta_yyyy_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 282);

    auto ta_yyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 283);

    auto ta_yyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 284);

    auto ta_yyyy_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 285);

    auto ta_yyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 286);

    auto ta_yyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 287);

    auto ta_yyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 288);

    auto ta_yyyy_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 289);

    auto ta_yyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 290);

    auto ta_yyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 291);

    auto ta_yyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 292);

    auto ta_yyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 293);

    auto ta_yyyy_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 294);

    auto ta_yyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 295);

    auto ta_yyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 296);

    auto ta_yyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 297);

    auto ta_yyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 298);

    auto ta_yyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 299);

    auto ta_yyyy_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 300);

    auto ta_yyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 301);

    auto ta_yyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 302);

    auto ta_yyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 303);

    auto ta_yyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 304);

    auto ta_yyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 305);

    auto ta_yyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 306);

    auto ta_yyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 307);

    auto ta_yyyz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 309);

    auto ta_yyyz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 310);

    auto ta_yyyz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 311);

    auto ta_yyyz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 313);

    auto ta_yyyz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 314);

    auto ta_yyyz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 317);

    auto ta_yyyz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 318);

    auto ta_yyyz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 322);

    auto ta_yyyz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 323);

    auto ta_yyyz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 328);

    auto ta_yyyz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 329);

    auto ta_yyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 330);

    auto ta_yyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 331);

    auto ta_yyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 332);

    auto ta_yyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 333);

    auto ta_yyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 334);

    auto ta_yyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 335);

    auto ta_yyzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 336);

    auto ta_yyzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 337);

    auto ta_yyzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 338);

    auto ta_yyzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 339);

    auto ta_yyzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 340);

    auto ta_yyzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 341);

    auto ta_yyzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 342);

    auto ta_yyzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 343);

    auto ta_yyzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 344);

    auto ta_yyzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 345);

    auto ta_yyzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 346);

    auto ta_yyzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 347);

    auto ta_yyzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 348);

    auto ta_yyzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 349);

    auto ta_yyzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 350);

    auto ta_yyzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 351);

    auto ta_yyzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 352);

    auto ta_yyzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 353);

    auto ta_yyzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 354);

    auto ta_yyzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 355);

    auto ta_yyzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 356);

    auto ta_yyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 357);

    auto ta_yyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 358);

    auto ta_yyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 359);

    auto ta_yyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 360);

    auto ta_yyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 361);

    auto ta_yyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 362);

    auto ta_yyzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 363);

    auto ta_yzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 364);

    auto ta_yzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 366);

    auto ta_yzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 368);

    auto ta_yzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 369);

    auto ta_yzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 371);

    auto ta_yzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 372);

    auto ta_yzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 373);

    auto ta_yzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 375);

    auto ta_yzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 376);

    auto ta_yzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 377);

    auto ta_yzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 378);

    auto ta_yzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 380);

    auto ta_yzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 381);

    auto ta_yzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 382);

    auto ta_yzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 383);

    auto ta_yzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 384);

    auto ta_yzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 385);

    auto ta_yzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 386);

    auto ta_yzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 387);

    auto ta_yzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 388);

    auto ta_yzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 389);

    auto ta_yzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 390);

    auto ta_yzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 391);

    auto ta_zzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 392);

    auto ta_zzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 393);

    auto ta_zzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 394);

    auto ta_zzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 395);

    auto ta_zzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 396);

    auto ta_zzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 397);

    auto ta_zzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 398);

    auto ta_zzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 399);

    auto ta_zzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 400);

    auto ta_zzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 401);

    auto ta_zzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 402);

    auto ta_zzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 403);

    auto ta_zzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 404);

    auto ta_zzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 405);

    auto ta_zzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 406);

    auto ta_zzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 407);

    auto ta_zzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 408);

    auto ta_zzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 409);

    auto ta_zzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 410);

    auto ta_zzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 411);

    auto ta_zzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 412);

    auto ta_zzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 413);

    auto ta_zzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 414);

    auto ta_zzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 415);

    auto ta_zzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 416);

    auto ta_zzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 417);

    auto ta_zzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 418);

    auto ta_zzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 419);

    // Set up components of auxiliary buffer : HH

    auto ta_xxxxx_xxxxx_0 = pbuffer.data(idx_npot_0_hh);

    auto ta_xxxxx_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 1);

    auto ta_xxxxx_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 2);

    auto ta_xxxxx_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 3);

    auto ta_xxxxx_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 4);

    auto ta_xxxxx_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 5);

    auto ta_xxxxx_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 6);

    auto ta_xxxxx_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 7);

    auto ta_xxxxx_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 8);

    auto ta_xxxxx_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 9);

    auto ta_xxxxx_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 10);

    auto ta_xxxxx_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 11);

    auto ta_xxxxx_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 12);

    auto ta_xxxxx_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 13);

    auto ta_xxxxx_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 14);

    auto ta_xxxxx_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 15);

    auto ta_xxxxx_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 16);

    auto ta_xxxxx_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 17);

    auto ta_xxxxx_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 18);

    auto ta_xxxxx_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 19);

    auto ta_xxxxx_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 20);

    auto ta_xxxxz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 44);

    auto ta_xxxxz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 46);

    auto ta_xxxxz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 47);

    auto ta_xxxxz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 49);

    auto ta_xxxxz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 50);

    auto ta_xxxxz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 51);

    auto ta_xxxxz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 53);

    auto ta_xxxxz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 54);

    auto ta_xxxxz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 55);

    auto ta_xxxxz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 56);

    auto ta_xxxyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 64);

    auto ta_xxxyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 66);

    auto ta_xxxyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 67);

    auto ta_xxxyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 69);

    auto ta_xxxyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 70);

    auto ta_xxxyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 71);

    auto ta_xxxyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 73);

    auto ta_xxxyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 74);

    auto ta_xxxyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 75);

    auto ta_xxxyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 76);

    auto ta_xxxyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 78);

    auto ta_xxxyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 79);

    auto ta_xxxyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 80);

    auto ta_xxxyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 81);

    auto ta_xxxyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 82);

    auto ta_xxxzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 105);

    auto ta_xxxzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 106);

    auto ta_xxxzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 107);

    auto ta_xxxzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 108);

    auto ta_xxxzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 109);

    auto ta_xxxzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 110);

    auto ta_xxxzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 111);

    auto ta_xxxzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 112);

    auto ta_xxxzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 113);

    auto ta_xxxzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 114);

    auto ta_xxxzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 115);

    auto ta_xxxzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 116);

    auto ta_xxxzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 117);

    auto ta_xxxzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 118);

    auto ta_xxxzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 119);

    auto ta_xxxzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 121);

    auto ta_xxxzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 122);

    auto ta_xxxzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 123);

    auto ta_xxxzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 124);

    auto ta_xxxzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 125);

    auto ta_xxyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 127);

    auto ta_xxyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 129);

    auto ta_xxyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 130);

    auto ta_xxyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 132);

    auto ta_xxyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 133);

    auto ta_xxyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 134);

    auto ta_xxyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 136);

    auto ta_xxyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 137);

    auto ta_xxyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 138);

    auto ta_xxyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 139);

    auto ta_xxyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 141);

    auto ta_xxyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 142);

    auto ta_xxyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 143);

    auto ta_xxyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 144);

    auto ta_xxyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 145);

    auto ta_xxzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 189);

    auto ta_xxzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 190);

    auto ta_xxzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 191);

    auto ta_xxzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 192);

    auto ta_xxzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 193);

    auto ta_xxzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 194);

    auto ta_xxzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 195);

    auto ta_xxzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 196);

    auto ta_xxzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 197);

    auto ta_xxzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 198);

    auto ta_xxzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 199);

    auto ta_xxzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 200);

    auto ta_xxzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 201);

    auto ta_xxzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 202);

    auto ta_xxzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 203);

    auto ta_xxzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 205);

    auto ta_xxzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 206);

    auto ta_xxzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 207);

    auto ta_xxzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 208);

    auto ta_xxzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 209);

    auto ta_xyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 211);

    auto ta_xyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 213);

    auto ta_xyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 214);

    auto ta_xyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 216);

    auto ta_xyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 217);

    auto ta_xyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 218);

    auto ta_xyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 220);

    auto ta_xyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 221);

    auto ta_xyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 222);

    auto ta_xyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 223);

    auto ta_xyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 225);

    auto ta_xyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 226);

    auto ta_xyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 227);

    auto ta_xyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 228);

    auto ta_xyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 229);

    auto ta_xyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 256);

    auto ta_xyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 259);

    auto ta_xyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 260);

    auto ta_xyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 263);

    auto ta_xyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 264);

    auto ta_xyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 265);

    auto ta_xyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 268);

    auto ta_xyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 269);

    auto ta_xyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 270);

    auto ta_xyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 271);

    auto ta_xzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 296);

    auto ta_xzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 298);

    auto ta_xzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 299);

    auto ta_xzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 301);

    auto ta_xzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 302);

    auto ta_xzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 303);

    auto ta_xzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 305);

    auto ta_xzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 306);

    auto ta_xzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 307);

    auto ta_xzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 308);

    auto ta_xzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 310);

    auto ta_xzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 311);

    auto ta_xzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 312);

    auto ta_xzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 313);

    auto ta_xzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 314);

    auto ta_yyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 315);

    auto ta_yyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 316);

    auto ta_yyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 317);

    auto ta_yyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 318);

    auto ta_yyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 319);

    auto ta_yyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 320);

    auto ta_yyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 321);

    auto ta_yyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 322);

    auto ta_yyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 323);

    auto ta_yyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 324);

    auto ta_yyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 325);

    auto ta_yyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 326);

    auto ta_yyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 327);

    auto ta_yyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 328);

    auto ta_yyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 329);

    auto ta_yyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 330);

    auto ta_yyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 331);

    auto ta_yyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 332);

    auto ta_yyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 333);

    auto ta_yyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 334);

    auto ta_yyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 335);

    auto ta_yyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 338);

    auto ta_yyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 340);

    auto ta_yyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 341);

    auto ta_yyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 343);

    auto ta_yyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 344);

    auto ta_yyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 345);

    auto ta_yyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 347);

    auto ta_yyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 348);

    auto ta_yyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 349);

    auto ta_yyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 350);

    auto ta_yyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 352);

    auto ta_yyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 353);

    auto ta_yyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 354);

    auto ta_yyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 355);

    auto ta_yyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 356);

    auto ta_yyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 357);

    auto ta_yyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 358);

    auto ta_yyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 359);

    auto ta_yyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 360);

    auto ta_yyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 361);

    auto ta_yyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 362);

    auto ta_yyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 363);

    auto ta_yyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 364);

    auto ta_yyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 365);

    auto ta_yyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 366);

    auto ta_yyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 367);

    auto ta_yyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 368);

    auto ta_yyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 369);

    auto ta_yyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 370);

    auto ta_yyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 371);

    auto ta_yyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 372);

    auto ta_yyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 373);

    auto ta_yyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 374);

    auto ta_yyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 375);

    auto ta_yyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 376);

    auto ta_yyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 377);

    auto ta_yyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 378);

    auto ta_yyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 379);

    auto ta_yyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 380);

    auto ta_yyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 381);

    auto ta_yyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 382);

    auto ta_yyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 383);

    auto ta_yyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 384);

    auto ta_yyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 385);

    auto ta_yyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 386);

    auto ta_yyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 387);

    auto ta_yyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 388);

    auto ta_yyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 389);

    auto ta_yyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 390);

    auto ta_yyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 391);

    auto ta_yyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 392);

    auto ta_yyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 393);

    auto ta_yyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 394);

    auto ta_yyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 395);

    auto ta_yyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 396);

    auto ta_yyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 397);

    auto ta_yyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 398);

    auto ta_yzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 400);

    auto ta_yzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 401);

    auto ta_yzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 402);

    auto ta_yzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 403);

    auto ta_yzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 404);

    auto ta_yzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 405);

    auto ta_yzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 406);

    auto ta_yzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 407);

    auto ta_yzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 408);

    auto ta_yzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 409);

    auto ta_yzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 410);

    auto ta_yzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 411);

    auto ta_yzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 412);

    auto ta_yzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 413);

    auto ta_yzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 414);

    auto ta_yzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 415);

    auto ta_yzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 416);

    auto ta_yzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 417);

    auto ta_yzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 418);

    auto ta_yzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 419);

    auto ta_zzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 420);

    auto ta_zzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 421);

    auto ta_zzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 422);

    auto ta_zzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 423);

    auto ta_zzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 424);

    auto ta_zzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 425);

    auto ta_zzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 426);

    auto ta_zzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 427);

    auto ta_zzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 428);

    auto ta_zzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 429);

    auto ta_zzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 430);

    auto ta_zzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 431);

    auto ta_zzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 432);

    auto ta_zzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 433);

    auto ta_zzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 434);

    auto ta_zzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 435);

    auto ta_zzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 436);

    auto ta_zzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 437);

    auto ta_zzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 438);

    auto ta_zzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 439);

    auto ta_zzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 440);

    // Set up components of auxiliary buffer : HH

    auto ta_xxxxx_xxxxx_1 = pbuffer.data(idx_npot_1_hh);

    auto ta_xxxxx_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 1);

    auto ta_xxxxx_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 2);

    auto ta_xxxxx_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 3);

    auto ta_xxxxx_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 4);

    auto ta_xxxxx_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 5);

    auto ta_xxxxx_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 6);

    auto ta_xxxxx_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 7);

    auto ta_xxxxx_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 8);

    auto ta_xxxxx_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 9);

    auto ta_xxxxx_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 10);

    auto ta_xxxxx_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 11);

    auto ta_xxxxx_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 12);

    auto ta_xxxxx_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 13);

    auto ta_xxxxx_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 14);

    auto ta_xxxxx_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 15);

    auto ta_xxxxx_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 16);

    auto ta_xxxxx_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 17);

    auto ta_xxxxx_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 18);

    auto ta_xxxxx_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 19);

    auto ta_xxxxx_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 20);

    auto ta_xxxxz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 44);

    auto ta_xxxxz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 46);

    auto ta_xxxxz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 47);

    auto ta_xxxxz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 49);

    auto ta_xxxxz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 50);

    auto ta_xxxxz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 51);

    auto ta_xxxxz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 53);

    auto ta_xxxxz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 54);

    auto ta_xxxxz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 55);

    auto ta_xxxxz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 56);

    auto ta_xxxyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 64);

    auto ta_xxxyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 66);

    auto ta_xxxyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 67);

    auto ta_xxxyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 69);

    auto ta_xxxyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 70);

    auto ta_xxxyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 71);

    auto ta_xxxyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 73);

    auto ta_xxxyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 74);

    auto ta_xxxyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 75);

    auto ta_xxxyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 76);

    auto ta_xxxyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 78);

    auto ta_xxxyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 79);

    auto ta_xxxyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 80);

    auto ta_xxxyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 81);

    auto ta_xxxyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 82);

    auto ta_xxxzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 105);

    auto ta_xxxzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 106);

    auto ta_xxxzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 107);

    auto ta_xxxzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 108);

    auto ta_xxxzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 109);

    auto ta_xxxzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 110);

    auto ta_xxxzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 111);

    auto ta_xxxzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 112);

    auto ta_xxxzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 113);

    auto ta_xxxzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 114);

    auto ta_xxxzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 115);

    auto ta_xxxzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 116);

    auto ta_xxxzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 117);

    auto ta_xxxzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 118);

    auto ta_xxxzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 119);

    auto ta_xxxzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 121);

    auto ta_xxxzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 122);

    auto ta_xxxzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 123);

    auto ta_xxxzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 124);

    auto ta_xxxzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 125);

    auto ta_xxyyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 127);

    auto ta_xxyyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 129);

    auto ta_xxyyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 130);

    auto ta_xxyyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 132);

    auto ta_xxyyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 133);

    auto ta_xxyyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 134);

    auto ta_xxyyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 136);

    auto ta_xxyyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 137);

    auto ta_xxyyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 138);

    auto ta_xxyyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 139);

    auto ta_xxyyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 141);

    auto ta_xxyyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 142);

    auto ta_xxyyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 143);

    auto ta_xxyyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 144);

    auto ta_xxyyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 145);

    auto ta_xxzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 189);

    auto ta_xxzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 190);

    auto ta_xxzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 191);

    auto ta_xxzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 192);

    auto ta_xxzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 193);

    auto ta_xxzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 194);

    auto ta_xxzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 195);

    auto ta_xxzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 196);

    auto ta_xxzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 197);

    auto ta_xxzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 198);

    auto ta_xxzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 199);

    auto ta_xxzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 200);

    auto ta_xxzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 201);

    auto ta_xxzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 202);

    auto ta_xxzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 203);

    auto ta_xxzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 205);

    auto ta_xxzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 206);

    auto ta_xxzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 207);

    auto ta_xxzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 208);

    auto ta_xxzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 209);

    auto ta_xyyyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 211);

    auto ta_xyyyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 213);

    auto ta_xyyyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 214);

    auto ta_xyyyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 216);

    auto ta_xyyyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 217);

    auto ta_xyyyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 218);

    auto ta_xyyyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 220);

    auto ta_xyyyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 221);

    auto ta_xyyyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 222);

    auto ta_xyyyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 223);

    auto ta_xyyyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 225);

    auto ta_xyyyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 226);

    auto ta_xyyyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 227);

    auto ta_xyyyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 228);

    auto ta_xyyyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 229);

    auto ta_xyyzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 256);

    auto ta_xyyzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 259);

    auto ta_xyyzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 260);

    auto ta_xyyzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 263);

    auto ta_xyyzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 264);

    auto ta_xyyzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 265);

    auto ta_xyyzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 268);

    auto ta_xyyzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 269);

    auto ta_xyyzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 270);

    auto ta_xyyzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 271);

    auto ta_xzzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 296);

    auto ta_xzzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 298);

    auto ta_xzzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 299);

    auto ta_xzzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 301);

    auto ta_xzzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 302);

    auto ta_xzzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 303);

    auto ta_xzzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 305);

    auto ta_xzzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 306);

    auto ta_xzzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 307);

    auto ta_xzzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 308);

    auto ta_xzzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 310);

    auto ta_xzzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 311);

    auto ta_xzzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 312);

    auto ta_xzzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 313);

    auto ta_xzzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 314);

    auto ta_yyyyy_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 315);

    auto ta_yyyyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 316);

    auto ta_yyyyy_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 317);

    auto ta_yyyyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 318);

    auto ta_yyyyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 319);

    auto ta_yyyyy_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 320);

    auto ta_yyyyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 321);

    auto ta_yyyyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 322);

    auto ta_yyyyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 323);

    auto ta_yyyyy_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 324);

    auto ta_yyyyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 325);

    auto ta_yyyyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 326);

    auto ta_yyyyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 327);

    auto ta_yyyyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 328);

    auto ta_yyyyy_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 329);

    auto ta_yyyyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 330);

    auto ta_yyyyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 331);

    auto ta_yyyyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 332);

    auto ta_yyyyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 333);

    auto ta_yyyyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 334);

    auto ta_yyyyy_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 335);

    auto ta_yyyyz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 338);

    auto ta_yyyyz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 340);

    auto ta_yyyyz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 341);

    auto ta_yyyyz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 343);

    auto ta_yyyyz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 344);

    auto ta_yyyyz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 345);

    auto ta_yyyyz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 347);

    auto ta_yyyyz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 348);

    auto ta_yyyyz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 349);

    auto ta_yyyyz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 350);

    auto ta_yyyyz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 352);

    auto ta_yyyyz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 353);

    auto ta_yyyyz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 354);

    auto ta_yyyyz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 355);

    auto ta_yyyyz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 356);

    auto ta_yyyzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 357);

    auto ta_yyyzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 358);

    auto ta_yyyzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 359);

    auto ta_yyyzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 360);

    auto ta_yyyzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 361);

    auto ta_yyyzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 362);

    auto ta_yyyzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 363);

    auto ta_yyyzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 364);

    auto ta_yyyzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 365);

    auto ta_yyyzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 366);

    auto ta_yyyzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 367);

    auto ta_yyyzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 368);

    auto ta_yyyzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 369);

    auto ta_yyyzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 370);

    auto ta_yyyzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 371);

    auto ta_yyyzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 372);

    auto ta_yyyzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 373);

    auto ta_yyyzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 374);

    auto ta_yyyzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 375);

    auto ta_yyyzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 376);

    auto ta_yyyzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 377);

    auto ta_yyzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 378);

    auto ta_yyzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 379);

    auto ta_yyzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 380);

    auto ta_yyzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 381);

    auto ta_yyzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 382);

    auto ta_yyzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 383);

    auto ta_yyzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 384);

    auto ta_yyzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 385);

    auto ta_yyzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 386);

    auto ta_yyzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 387);

    auto ta_yyzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 388);

    auto ta_yyzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 389);

    auto ta_yyzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 390);

    auto ta_yyzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 391);

    auto ta_yyzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 392);

    auto ta_yyzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 393);

    auto ta_yyzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 394);

    auto ta_yyzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 395);

    auto ta_yyzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 396);

    auto ta_yyzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 397);

    auto ta_yyzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 398);

    auto ta_yzzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 400);

    auto ta_yzzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 401);

    auto ta_yzzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 402);

    auto ta_yzzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 403);

    auto ta_yzzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 404);

    auto ta_yzzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 405);

    auto ta_yzzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 406);

    auto ta_yzzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 407);

    auto ta_yzzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 408);

    auto ta_yzzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 409);

    auto ta_yzzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 410);

    auto ta_yzzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 411);

    auto ta_yzzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 412);

    auto ta_yzzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 413);

    auto ta_yzzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 414);

    auto ta_yzzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 415);

    auto ta_yzzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 416);

    auto ta_yzzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 417);

    auto ta_yzzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 418);

    auto ta_yzzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 419);

    auto ta_zzzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 420);

    auto ta_zzzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 421);

    auto ta_zzzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 422);

    auto ta_zzzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 423);

    auto ta_zzzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 424);

    auto ta_zzzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 425);

    auto ta_zzzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 426);

    auto ta_zzzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 427);

    auto ta_zzzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 428);

    auto ta_zzzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 429);

    auto ta_zzzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 430);

    auto ta_zzzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 431);

    auto ta_zzzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 432);

    auto ta_zzzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 433);

    auto ta_zzzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 434);

    auto ta_zzzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 435);

    auto ta_zzzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 436);

    auto ta_zzzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 437);

    auto ta_zzzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 438);

    auto ta_zzzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 439);

    auto ta_zzzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 440);

    // Set up components of auxiliary buffer : HI

    auto ta_xxxxx_xxxxxx_0 = pbuffer.data(idx_npot_0_hi);

    auto ta_xxxxx_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 1);

    auto ta_xxxxx_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 2);

    auto ta_xxxxx_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 3);

    auto ta_xxxxx_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 4);

    auto ta_xxxxx_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 5);

    auto ta_xxxxx_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 6);

    auto ta_xxxxx_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 7);

    auto ta_xxxxx_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 8);

    auto ta_xxxxx_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 9);

    auto ta_xxxxx_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 10);

    auto ta_xxxxx_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 11);

    auto ta_xxxxx_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 12);

    auto ta_xxxxx_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 13);

    auto ta_xxxxx_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 14);

    auto ta_xxxxx_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 15);

    auto ta_xxxxx_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 16);

    auto ta_xxxxx_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 17);

    auto ta_xxxxx_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 18);

    auto ta_xxxxx_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 19);

    auto ta_xxxxx_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 20);

    auto ta_xxxxx_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 21);

    auto ta_xxxxx_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 22);

    auto ta_xxxxx_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 23);

    auto ta_xxxxx_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 24);

    auto ta_xxxxx_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 25);

    auto ta_xxxxx_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 26);

    auto ta_xxxxx_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 27);

    auto ta_xxxxy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 28);

    auto ta_xxxxy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 29);

    auto ta_xxxxy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 30);

    auto ta_xxxxy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 31);

    auto ta_xxxxy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 33);

    auto ta_xxxxy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 34);

    auto ta_xxxxy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 37);

    auto ta_xxxxy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 38);

    auto ta_xxxxy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 42);

    auto ta_xxxxy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 43);

    auto ta_xxxxy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 48);

    auto ta_xxxxy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 49);

    auto ta_xxxxy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 50);

    auto ta_xxxxy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 51);

    auto ta_xxxxy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 52);

    auto ta_xxxxy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 53);

    auto ta_xxxxy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 54);

    auto ta_xxxxz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 56);

    auto ta_xxxxz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 57);

    auto ta_xxxxz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 58);

    auto ta_xxxxz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 59);

    auto ta_xxxxz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 60);

    auto ta_xxxxz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 61);

    auto ta_xxxxz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 62);

    auto ta_xxxxz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 63);

    auto ta_xxxxz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 64);

    auto ta_xxxxz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 65);

    auto ta_xxxxz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 66);

    auto ta_xxxxz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 67);

    auto ta_xxxxz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 68);

    auto ta_xxxxz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 69);

    auto ta_xxxxz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 70);

    auto ta_xxxxz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 71);

    auto ta_xxxxz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 72);

    auto ta_xxxxz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 73);

    auto ta_xxxxz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 74);

    auto ta_xxxxz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 75);

    auto ta_xxxxz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 76);

    auto ta_xxxxz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 78);

    auto ta_xxxxz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 79);

    auto ta_xxxxz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 80);

    auto ta_xxxxz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 81);

    auto ta_xxxxz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 82);

    auto ta_xxxxz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 83);

    auto ta_xxxyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 84);

    auto ta_xxxyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 85);

    auto ta_xxxyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 86);

    auto ta_xxxyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 87);

    auto ta_xxxyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 88);

    auto ta_xxxyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 89);

    auto ta_xxxyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 90);

    auto ta_xxxyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 91);

    auto ta_xxxyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 92);

    auto ta_xxxyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 93);

    auto ta_xxxyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 94);

    auto ta_xxxyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 95);

    auto ta_xxxyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 96);

    auto ta_xxxyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 97);

    auto ta_xxxyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 98);

    auto ta_xxxyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 99);

    auto ta_xxxyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 100);

    auto ta_xxxyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 101);

    auto ta_xxxyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 102);

    auto ta_xxxyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 103);

    auto ta_xxxyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 104);

    auto ta_xxxyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 105);

    auto ta_xxxyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 106);

    auto ta_xxxyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 107);

    auto ta_xxxyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 108);

    auto ta_xxxyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 109);

    auto ta_xxxyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 110);

    auto ta_xxxyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 111);

    auto ta_xxxyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 114);

    auto ta_xxxyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 117);

    auto ta_xxxyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 121);

    auto ta_xxxyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 126);

    auto ta_xxxyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 132);

    auto ta_xxxyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 134);

    auto ta_xxxyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 135);

    auto ta_xxxyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 136);

    auto ta_xxxyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 137);

    auto ta_xxxyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 138);

    auto ta_xxxzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 140);

    auto ta_xxxzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 141);

    auto ta_xxxzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 142);

    auto ta_xxxzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 143);

    auto ta_xxxzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 144);

    auto ta_xxxzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 145);

    auto ta_xxxzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 146);

    auto ta_xxxzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 147);

    auto ta_xxxzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 148);

    auto ta_xxxzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 149);

    auto ta_xxxzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 150);

    auto ta_xxxzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 151);

    auto ta_xxxzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 152);

    auto ta_xxxzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 153);

    auto ta_xxxzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 154);

    auto ta_xxxzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 155);

    auto ta_xxxzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 156);

    auto ta_xxxzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 157);

    auto ta_xxxzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 158);

    auto ta_xxxzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 159);

    auto ta_xxxzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 160);

    auto ta_xxxzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 161);

    auto ta_xxxzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 162);

    auto ta_xxxzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 163);

    auto ta_xxxzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 164);

    auto ta_xxxzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 165);

    auto ta_xxxzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 166);

    auto ta_xxxzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 167);

    auto ta_xxyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 168);

    auto ta_xxyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 169);

    auto ta_xxyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 170);

    auto ta_xxyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 171);

    auto ta_xxyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 172);

    auto ta_xxyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 173);

    auto ta_xxyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 174);

    auto ta_xxyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 175);

    auto ta_xxyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 176);

    auto ta_xxyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 177);

    auto ta_xxyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 178);

    auto ta_xxyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 179);

    auto ta_xxyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 180);

    auto ta_xxyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 181);

    auto ta_xxyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 182);

    auto ta_xxyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 183);

    auto ta_xxyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 184);

    auto ta_xxyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 185);

    auto ta_xxyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 186);

    auto ta_xxyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 187);

    auto ta_xxyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 188);

    auto ta_xxyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 189);

    auto ta_xxyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 190);

    auto ta_xxyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 191);

    auto ta_xxyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 192);

    auto ta_xxyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 193);

    auto ta_xxyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 194);

    auto ta_xxyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 195);

    auto ta_xxyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 197);

    auto ta_xxyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 198);

    auto ta_xxyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 199);

    auto ta_xxyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 201);

    auto ta_xxyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 202);

    auto ta_xxyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 205);

    auto ta_xxyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 206);

    auto ta_xxyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 210);

    auto ta_xxyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 211);

    auto ta_xxyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 216);

    auto ta_xxyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 218);

    auto ta_xxyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 219);

    auto ta_xxyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 220);

    auto ta_xxyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 221);

    auto ta_xxyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 222);

    auto ta_xxyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 223);

    auto ta_xxyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 224);

    auto ta_xxyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 226);

    auto ta_xxyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 229);

    auto ta_xxyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 233);

    auto ta_xxyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 238);

    auto ta_xxyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 244);

    auto ta_xxyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 245);

    auto ta_xxyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 246);

    auto ta_xxyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 247);

    auto ta_xxyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 248);

    auto ta_xxyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 249);

    auto ta_xxyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 250);

    auto ta_xxzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 252);

    auto ta_xxzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 253);

    auto ta_xxzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 254);

    auto ta_xxzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 255);

    auto ta_xxzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 256);

    auto ta_xxzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 257);

    auto ta_xxzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 258);

    auto ta_xxzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 259);

    auto ta_xxzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 260);

    auto ta_xxzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 261);

    auto ta_xxzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 262);

    auto ta_xxzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 263);

    auto ta_xxzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 264);

    auto ta_xxzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 265);

    auto ta_xxzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 266);

    auto ta_xxzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 267);

    auto ta_xxzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 268);

    auto ta_xxzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 269);

    auto ta_xxzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 270);

    auto ta_xxzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 271);

    auto ta_xxzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 272);

    auto ta_xxzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 273);

    auto ta_xxzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 274);

    auto ta_xxzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 275);

    auto ta_xxzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 276);

    auto ta_xxzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 277);

    auto ta_xxzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 278);

    auto ta_xxzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 279);

    auto ta_xyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 280);

    auto ta_xyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 281);

    auto ta_xyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 283);

    auto ta_xyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 284);

    auto ta_xyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 286);

    auto ta_xyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 287);

    auto ta_xyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 288);

    auto ta_xyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 290);

    auto ta_xyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 291);

    auto ta_xyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 292);

    auto ta_xyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 293);

    auto ta_xyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 295);

    auto ta_xyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 296);

    auto ta_xyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 297);

    auto ta_xyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 298);

    auto ta_xyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 299);

    auto ta_xyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 301);

    auto ta_xyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 302);

    auto ta_xyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 303);

    auto ta_xyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 304);

    auto ta_xyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 305);

    auto ta_xyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 306);

    auto ta_xyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 307);

    auto ta_xyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 330);

    auto ta_xyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 331);

    auto ta_xyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 332);

    auto ta_xyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 333);

    auto ta_xyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 334);

    auto ta_xyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 335);

    auto ta_xyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 340);

    auto ta_xyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 343);

    auto ta_xyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 344);

    auto ta_xyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 347);

    auto ta_xyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 348);

    auto ta_xyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 349);

    auto ta_xyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 352);

    auto ta_xyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 353);

    auto ta_xyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 354);

    auto ta_xyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 355);

    auto ta_xyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 357);

    auto ta_xyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 358);

    auto ta_xyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 359);

    auto ta_xyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 360);

    auto ta_xyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 361);

    auto ta_xyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 362);

    auto ta_xyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 363);

    auto ta_xyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 385);

    auto ta_xyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 386);

    auto ta_xyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 387);

    auto ta_xyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 388);

    auto ta_xyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 389);

    auto ta_xyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 390);

    auto ta_xzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 392);

    auto ta_xzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 394);

    auto ta_xzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 396);

    auto ta_xzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 397);

    auto ta_xzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 399);

    auto ta_xzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 400);

    auto ta_xzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 401);

    auto ta_xzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 403);

    auto ta_xzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 404);

    auto ta_xzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 405);

    auto ta_xzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 406);

    auto ta_xzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 408);

    auto ta_xzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 409);

    auto ta_xzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 410);

    auto ta_xzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 411);

    auto ta_xzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 412);

    auto ta_xzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 413);

    auto ta_xzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 414);

    auto ta_xzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 415);

    auto ta_xzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 416);

    auto ta_xzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 417);

    auto ta_xzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 418);

    auto ta_xzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 419);

    auto ta_yyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 420);

    auto ta_yyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 421);

    auto ta_yyyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 422);

    auto ta_yyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 423);

    auto ta_yyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 424);

    auto ta_yyyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 425);

    auto ta_yyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 426);

    auto ta_yyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 427);

    auto ta_yyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 428);

    auto ta_yyyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 429);

    auto ta_yyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 430);

    auto ta_yyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 431);

    auto ta_yyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 432);

    auto ta_yyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 433);

    auto ta_yyyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 434);

    auto ta_yyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 435);

    auto ta_yyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 436);

    auto ta_yyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 437);

    auto ta_yyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 438);

    auto ta_yyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 439);

    auto ta_yyyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 440);

    auto ta_yyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 441);

    auto ta_yyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 442);

    auto ta_yyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 443);

    auto ta_yyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 444);

    auto ta_yyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 445);

    auto ta_yyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 446);

    auto ta_yyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 447);

    auto ta_yyyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 449);

    auto ta_yyyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 450);

    auto ta_yyyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 451);

    auto ta_yyyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 452);

    auto ta_yyyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 453);

    auto ta_yyyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 454);

    auto ta_yyyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 455);

    auto ta_yyyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 456);

    auto ta_yyyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 457);

    auto ta_yyyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 458);

    auto ta_yyyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 459);

    auto ta_yyyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 460);

    auto ta_yyyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 461);

    auto ta_yyyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 462);

    auto ta_yyyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 463);

    auto ta_yyyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 464);

    auto ta_yyyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 465);

    auto ta_yyyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 466);

    auto ta_yyyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 467);

    auto ta_yyyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 468);

    auto ta_yyyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 469);

    auto ta_yyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 470);

    auto ta_yyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 471);

    auto ta_yyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 472);

    auto ta_yyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 473);

    auto ta_yyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 474);

    auto ta_yyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 475);

    auto ta_yyyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 476);

    auto ta_yyyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 477);

    auto ta_yyyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 478);

    auto ta_yyyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 479);

    auto ta_yyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 480);

    auto ta_yyyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 481);

    auto ta_yyyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 482);

    auto ta_yyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 483);

    auto ta_yyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 484);

    auto ta_yyyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 485);

    auto ta_yyyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 486);

    auto ta_yyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 487);

    auto ta_yyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 488);

    auto ta_yyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 489);

    auto ta_yyyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 490);

    auto ta_yyyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 491);

    auto ta_yyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 492);

    auto ta_yyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 493);

    auto ta_yyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 494);

    auto ta_yyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 495);

    auto ta_yyyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 496);

    auto ta_yyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 497);

    auto ta_yyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 498);

    auto ta_yyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 499);

    auto ta_yyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 500);

    auto ta_yyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 501);

    auto ta_yyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 502);

    auto ta_yyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 503);

    auto ta_yyzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 504);

    auto ta_yyzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 505);

    auto ta_yyzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 506);

    auto ta_yyzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 507);

    auto ta_yyzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 508);

    auto ta_yyzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 509);

    auto ta_yyzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 510);

    auto ta_yyzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 511);

    auto ta_yyzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 512);

    auto ta_yyzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 513);

    auto ta_yyzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 514);

    auto ta_yyzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 515);

    auto ta_yyzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 516);

    auto ta_yyzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 517);

    auto ta_yyzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 518);

    auto ta_yyzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 519);

    auto ta_yyzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 520);

    auto ta_yyzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 521);

    auto ta_yyzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 522);

    auto ta_yyzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 523);

    auto ta_yyzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 524);

    auto ta_yyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 525);

    auto ta_yyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 526);

    auto ta_yyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 527);

    auto ta_yyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 528);

    auto ta_yyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 529);

    auto ta_yyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 530);

    auto ta_yyzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 531);

    auto ta_yzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 532);

    auto ta_yzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 533);

    auto ta_yzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 534);

    auto ta_yzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 535);

    auto ta_yzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 536);

    auto ta_yzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 537);

    auto ta_yzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 538);

    auto ta_yzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 539);

    auto ta_yzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 540);

    auto ta_yzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 541);

    auto ta_yzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 542);

    auto ta_yzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 543);

    auto ta_yzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 544);

    auto ta_yzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 545);

    auto ta_yzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 546);

    auto ta_yzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 547);

    auto ta_yzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 548);

    auto ta_yzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 549);

    auto ta_yzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 550);

    auto ta_yzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 551);

    auto ta_yzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 552);

    auto ta_yzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 553);

    auto ta_yzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 554);

    auto ta_yzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 555);

    auto ta_yzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 556);

    auto ta_yzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 557);

    auto ta_yzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 558);

    auto ta_yzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 559);

    auto ta_zzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 560);

    auto ta_zzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 561);

    auto ta_zzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 562);

    auto ta_zzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 563);

    auto ta_zzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 564);

    auto ta_zzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 565);

    auto ta_zzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 566);

    auto ta_zzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 567);

    auto ta_zzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 568);

    auto ta_zzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 569);

    auto ta_zzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 570);

    auto ta_zzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 571);

    auto ta_zzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 572);

    auto ta_zzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 573);

    auto ta_zzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 574);

    auto ta_zzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 575);

    auto ta_zzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 576);

    auto ta_zzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 577);

    auto ta_zzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 578);

    auto ta_zzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 579);

    auto ta_zzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 580);

    auto ta_zzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 581);

    auto ta_zzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 582);

    auto ta_zzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 583);

    auto ta_zzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 584);

    auto ta_zzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 585);

    auto ta_zzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 586);

    auto ta_zzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 587);

    // Set up components of auxiliary buffer : HI

    auto ta_xxxxx_xxxxxx_1 = pbuffer.data(idx_npot_1_hi);

    auto ta_xxxxx_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 1);

    auto ta_xxxxx_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 2);

    auto ta_xxxxx_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 3);

    auto ta_xxxxx_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 4);

    auto ta_xxxxx_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 5);

    auto ta_xxxxx_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 6);

    auto ta_xxxxx_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 7);

    auto ta_xxxxx_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 8);

    auto ta_xxxxx_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 9);

    auto ta_xxxxx_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 10);

    auto ta_xxxxx_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 11);

    auto ta_xxxxx_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 12);

    auto ta_xxxxx_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 13);

    auto ta_xxxxx_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 14);

    auto ta_xxxxx_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 15);

    auto ta_xxxxx_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 16);

    auto ta_xxxxx_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 17);

    auto ta_xxxxx_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 18);

    auto ta_xxxxx_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 19);

    auto ta_xxxxx_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 20);

    auto ta_xxxxx_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 21);

    auto ta_xxxxx_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 22);

    auto ta_xxxxx_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 23);

    auto ta_xxxxx_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 24);

    auto ta_xxxxx_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 25);

    auto ta_xxxxx_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 26);

    auto ta_xxxxx_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 27);

    auto ta_xxxxy_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 28);

    auto ta_xxxxy_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 29);

    auto ta_xxxxy_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 30);

    auto ta_xxxxy_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 31);

    auto ta_xxxxy_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 33);

    auto ta_xxxxy_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 34);

    auto ta_xxxxy_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 37);

    auto ta_xxxxy_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 38);

    auto ta_xxxxy_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 42);

    auto ta_xxxxy_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 43);

    auto ta_xxxxy_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 48);

    auto ta_xxxxy_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 49);

    auto ta_xxxxy_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 50);

    auto ta_xxxxy_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 51);

    auto ta_xxxxy_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 52);

    auto ta_xxxxy_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 53);

    auto ta_xxxxy_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 54);

    auto ta_xxxxz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 56);

    auto ta_xxxxz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 57);

    auto ta_xxxxz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 58);

    auto ta_xxxxz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 59);

    auto ta_xxxxz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 60);

    auto ta_xxxxz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 61);

    auto ta_xxxxz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 62);

    auto ta_xxxxz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 63);

    auto ta_xxxxz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 64);

    auto ta_xxxxz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 65);

    auto ta_xxxxz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 66);

    auto ta_xxxxz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 67);

    auto ta_xxxxz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 68);

    auto ta_xxxxz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 69);

    auto ta_xxxxz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 70);

    auto ta_xxxxz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 71);

    auto ta_xxxxz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 72);

    auto ta_xxxxz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 73);

    auto ta_xxxxz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 74);

    auto ta_xxxxz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 75);

    auto ta_xxxxz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 76);

    auto ta_xxxxz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 78);

    auto ta_xxxxz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 79);

    auto ta_xxxxz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 80);

    auto ta_xxxxz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 81);

    auto ta_xxxxz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 82);

    auto ta_xxxxz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 83);

    auto ta_xxxyy_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 84);

    auto ta_xxxyy_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 85);

    auto ta_xxxyy_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 86);

    auto ta_xxxyy_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 87);

    auto ta_xxxyy_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 88);

    auto ta_xxxyy_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 89);

    auto ta_xxxyy_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 90);

    auto ta_xxxyy_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 91);

    auto ta_xxxyy_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 92);

    auto ta_xxxyy_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 93);

    auto ta_xxxyy_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 94);

    auto ta_xxxyy_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 95);

    auto ta_xxxyy_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 96);

    auto ta_xxxyy_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 97);

    auto ta_xxxyy_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 98);

    auto ta_xxxyy_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 99);

    auto ta_xxxyy_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 100);

    auto ta_xxxyy_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 101);

    auto ta_xxxyy_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 102);

    auto ta_xxxyy_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 103);

    auto ta_xxxyy_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 104);

    auto ta_xxxyy_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 105);

    auto ta_xxxyy_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 106);

    auto ta_xxxyy_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 107);

    auto ta_xxxyy_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 108);

    auto ta_xxxyy_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 109);

    auto ta_xxxyy_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 110);

    auto ta_xxxyy_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 111);

    auto ta_xxxyz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 114);

    auto ta_xxxyz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 117);

    auto ta_xxxyz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 121);

    auto ta_xxxyz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 126);

    auto ta_xxxyz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 132);

    auto ta_xxxyz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 134);

    auto ta_xxxyz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 135);

    auto ta_xxxyz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 136);

    auto ta_xxxyz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 137);

    auto ta_xxxyz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 138);

    auto ta_xxxzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 140);

    auto ta_xxxzz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 141);

    auto ta_xxxzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 142);

    auto ta_xxxzz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 143);

    auto ta_xxxzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 144);

    auto ta_xxxzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 145);

    auto ta_xxxzz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 146);

    auto ta_xxxzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 147);

    auto ta_xxxzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 148);

    auto ta_xxxzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 149);

    auto ta_xxxzz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 150);

    auto ta_xxxzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 151);

    auto ta_xxxzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 152);

    auto ta_xxxzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 153);

    auto ta_xxxzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 154);

    auto ta_xxxzz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 155);

    auto ta_xxxzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 156);

    auto ta_xxxzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 157);

    auto ta_xxxzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 158);

    auto ta_xxxzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 159);

    auto ta_xxxzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 160);

    auto ta_xxxzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 161);

    auto ta_xxxzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 162);

    auto ta_xxxzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 163);

    auto ta_xxxzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 164);

    auto ta_xxxzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 165);

    auto ta_xxxzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 166);

    auto ta_xxxzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 167);

    auto ta_xxyyy_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 168);

    auto ta_xxyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 169);

    auto ta_xxyyy_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 170);

    auto ta_xxyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 171);

    auto ta_xxyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 172);

    auto ta_xxyyy_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 173);

    auto ta_xxyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 174);

    auto ta_xxyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 175);

    auto ta_xxyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 176);

    auto ta_xxyyy_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 177);

    auto ta_xxyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 178);

    auto ta_xxyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 179);

    auto ta_xxyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 180);

    auto ta_xxyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 181);

    auto ta_xxyyy_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 182);

    auto ta_xxyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 183);

    auto ta_xxyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 184);

    auto ta_xxyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 185);

    auto ta_xxyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 186);

    auto ta_xxyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 187);

    auto ta_xxyyy_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 188);

    auto ta_xxyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 189);

    auto ta_xxyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 190);

    auto ta_xxyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 191);

    auto ta_xxyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 192);

    auto ta_xxyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 193);

    auto ta_xxyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 194);

    auto ta_xxyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 195);

    auto ta_xxyyz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 197);

    auto ta_xxyyz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 198);

    auto ta_xxyyz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 199);

    auto ta_xxyyz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 201);

    auto ta_xxyyz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 202);

    auto ta_xxyyz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 205);

    auto ta_xxyyz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 206);

    auto ta_xxyyz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 210);

    auto ta_xxyyz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 211);

    auto ta_xxyyz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 216);

    auto ta_xxyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 218);

    auto ta_xxyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 219);

    auto ta_xxyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 220);

    auto ta_xxyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 221);

    auto ta_xxyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 222);

    auto ta_xxyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 223);

    auto ta_xxyzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 224);

    auto ta_xxyzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 226);

    auto ta_xxyzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 229);

    auto ta_xxyzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 233);

    auto ta_xxyzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 238);

    auto ta_xxyzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 244);

    auto ta_xxyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 245);

    auto ta_xxyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 246);

    auto ta_xxyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 247);

    auto ta_xxyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 248);

    auto ta_xxyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 249);

    auto ta_xxyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 250);

    auto ta_xxzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 252);

    auto ta_xxzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 253);

    auto ta_xxzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 254);

    auto ta_xxzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 255);

    auto ta_xxzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 256);

    auto ta_xxzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 257);

    auto ta_xxzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 258);

    auto ta_xxzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 259);

    auto ta_xxzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 260);

    auto ta_xxzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 261);

    auto ta_xxzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 262);

    auto ta_xxzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 263);

    auto ta_xxzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 264);

    auto ta_xxzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 265);

    auto ta_xxzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 266);

    auto ta_xxzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 267);

    auto ta_xxzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 268);

    auto ta_xxzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 269);

    auto ta_xxzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 270);

    auto ta_xxzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 271);

    auto ta_xxzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 272);

    auto ta_xxzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 273);

    auto ta_xxzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 274);

    auto ta_xxzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 275);

    auto ta_xxzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 276);

    auto ta_xxzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 277);

    auto ta_xxzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 278);

    auto ta_xxzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 279);

    auto ta_xyyyy_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 280);

    auto ta_xyyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 281);

    auto ta_xyyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 283);

    auto ta_xyyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 284);

    auto ta_xyyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 286);

    auto ta_xyyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 287);

    auto ta_xyyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 288);

    auto ta_xyyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 290);

    auto ta_xyyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 291);

    auto ta_xyyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 292);

    auto ta_xyyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 293);

    auto ta_xyyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 295);

    auto ta_xyyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 296);

    auto ta_xyyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 297);

    auto ta_xyyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 298);

    auto ta_xyyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 299);

    auto ta_xyyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 301);

    auto ta_xyyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 302);

    auto ta_xyyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 303);

    auto ta_xyyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 304);

    auto ta_xyyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 305);

    auto ta_xyyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 306);

    auto ta_xyyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 307);

    auto ta_xyyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 330);

    auto ta_xyyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 331);

    auto ta_xyyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 332);

    auto ta_xyyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 333);

    auto ta_xyyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 334);

    auto ta_xyyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 335);

    auto ta_xyyzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 340);

    auto ta_xyyzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 343);

    auto ta_xyyzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 344);

    auto ta_xyyzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 347);

    auto ta_xyyzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 348);

    auto ta_xyyzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 349);

    auto ta_xyyzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 352);

    auto ta_xyyzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 353);

    auto ta_xyyzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 354);

    auto ta_xyyzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 355);

    auto ta_xyyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 357);

    auto ta_xyyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 358);

    auto ta_xyyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 359);

    auto ta_xyyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 360);

    auto ta_xyyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 361);

    auto ta_xyyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 362);

    auto ta_xyyzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 363);

    auto ta_xyzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 385);

    auto ta_xyzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 386);

    auto ta_xyzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 387);

    auto ta_xyzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 388);

    auto ta_xyzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 389);

    auto ta_xyzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 390);

    auto ta_xzzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 392);

    auto ta_xzzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 394);

    auto ta_xzzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 396);

    auto ta_xzzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 397);

    auto ta_xzzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 399);

    auto ta_xzzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 400);

    auto ta_xzzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 401);

    auto ta_xzzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 403);

    auto ta_xzzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 404);

    auto ta_xzzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 405);

    auto ta_xzzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 406);

    auto ta_xzzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 408);

    auto ta_xzzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 409);

    auto ta_xzzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 410);

    auto ta_xzzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 411);

    auto ta_xzzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 412);

    auto ta_xzzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 413);

    auto ta_xzzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 414);

    auto ta_xzzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 415);

    auto ta_xzzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 416);

    auto ta_xzzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 417);

    auto ta_xzzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 418);

    auto ta_xzzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 419);

    auto ta_yyyyy_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 420);

    auto ta_yyyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 421);

    auto ta_yyyyy_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 422);

    auto ta_yyyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 423);

    auto ta_yyyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 424);

    auto ta_yyyyy_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 425);

    auto ta_yyyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 426);

    auto ta_yyyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 427);

    auto ta_yyyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 428);

    auto ta_yyyyy_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 429);

    auto ta_yyyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 430);

    auto ta_yyyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 431);

    auto ta_yyyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 432);

    auto ta_yyyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 433);

    auto ta_yyyyy_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 434);

    auto ta_yyyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 435);

    auto ta_yyyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 436);

    auto ta_yyyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 437);

    auto ta_yyyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 438);

    auto ta_yyyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 439);

    auto ta_yyyyy_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 440);

    auto ta_yyyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 441);

    auto ta_yyyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 442);

    auto ta_yyyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 443);

    auto ta_yyyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 444);

    auto ta_yyyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 445);

    auto ta_yyyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 446);

    auto ta_yyyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 447);

    auto ta_yyyyz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 449);

    auto ta_yyyyz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 450);

    auto ta_yyyyz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 451);

    auto ta_yyyyz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 452);

    auto ta_yyyyz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 453);

    auto ta_yyyyz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 454);

    auto ta_yyyyz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 455);

    auto ta_yyyyz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 456);

    auto ta_yyyyz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 457);

    auto ta_yyyyz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 458);

    auto ta_yyyyz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 459);

    auto ta_yyyyz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 460);

    auto ta_yyyyz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 461);

    auto ta_yyyyz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 462);

    auto ta_yyyyz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 463);

    auto ta_yyyyz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 464);

    auto ta_yyyyz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 465);

    auto ta_yyyyz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 466);

    auto ta_yyyyz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 467);

    auto ta_yyyyz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 468);

    auto ta_yyyyz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 469);

    auto ta_yyyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 470);

    auto ta_yyyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 471);

    auto ta_yyyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 472);

    auto ta_yyyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 473);

    auto ta_yyyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 474);

    auto ta_yyyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 475);

    auto ta_yyyzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 476);

    auto ta_yyyzz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 477);

    auto ta_yyyzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 478);

    auto ta_yyyzz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 479);

    auto ta_yyyzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 480);

    auto ta_yyyzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 481);

    auto ta_yyyzz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 482);

    auto ta_yyyzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 483);

    auto ta_yyyzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 484);

    auto ta_yyyzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 485);

    auto ta_yyyzz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 486);

    auto ta_yyyzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 487);

    auto ta_yyyzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 488);

    auto ta_yyyzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 489);

    auto ta_yyyzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 490);

    auto ta_yyyzz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 491);

    auto ta_yyyzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 492);

    auto ta_yyyzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 493);

    auto ta_yyyzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 494);

    auto ta_yyyzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 495);

    auto ta_yyyzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 496);

    auto ta_yyyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 497);

    auto ta_yyyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 498);

    auto ta_yyyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 499);

    auto ta_yyyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 500);

    auto ta_yyyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 501);

    auto ta_yyyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 502);

    auto ta_yyyzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 503);

    auto ta_yyzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 504);

    auto ta_yyzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 505);

    auto ta_yyzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 506);

    auto ta_yyzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 507);

    auto ta_yyzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 508);

    auto ta_yyzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 509);

    auto ta_yyzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 510);

    auto ta_yyzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 511);

    auto ta_yyzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 512);

    auto ta_yyzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 513);

    auto ta_yyzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 514);

    auto ta_yyzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 515);

    auto ta_yyzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 516);

    auto ta_yyzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 517);

    auto ta_yyzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 518);

    auto ta_yyzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 519);

    auto ta_yyzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 520);

    auto ta_yyzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 521);

    auto ta_yyzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 522);

    auto ta_yyzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 523);

    auto ta_yyzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 524);

    auto ta_yyzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 525);

    auto ta_yyzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 526);

    auto ta_yyzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 527);

    auto ta_yyzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 528);

    auto ta_yyzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 529);

    auto ta_yyzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 530);

    auto ta_yyzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 531);

    auto ta_yzzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 532);

    auto ta_yzzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 533);

    auto ta_yzzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 534);

    auto ta_yzzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 535);

    auto ta_yzzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 536);

    auto ta_yzzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 537);

    auto ta_yzzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 538);

    auto ta_yzzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 539);

    auto ta_yzzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 540);

    auto ta_yzzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 541);

    auto ta_yzzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 542);

    auto ta_yzzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 543);

    auto ta_yzzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 544);

    auto ta_yzzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 545);

    auto ta_yzzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 546);

    auto ta_yzzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 547);

    auto ta_yzzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 548);

    auto ta_yzzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 549);

    auto ta_yzzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 550);

    auto ta_yzzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 551);

    auto ta_yzzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 552);

    auto ta_yzzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 553);

    auto ta_yzzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 554);

    auto ta_yzzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 555);

    auto ta_yzzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 556);

    auto ta_yzzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 557);

    auto ta_yzzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 558);

    auto ta_yzzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 559);

    auto ta_zzzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_hi + 560);

    auto ta_zzzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_hi + 561);

    auto ta_zzzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_hi + 562);

    auto ta_zzzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_hi + 563);

    auto ta_zzzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_hi + 564);

    auto ta_zzzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_hi + 565);

    auto ta_zzzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_hi + 566);

    auto ta_zzzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_hi + 567);

    auto ta_zzzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_hi + 568);

    auto ta_zzzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_hi + 569);

    auto ta_zzzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_hi + 570);

    auto ta_zzzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_hi + 571);

    auto ta_zzzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_hi + 572);

    auto ta_zzzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_hi + 573);

    auto ta_zzzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_hi + 574);

    auto ta_zzzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_hi + 575);

    auto ta_zzzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_hi + 576);

    auto ta_zzzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_hi + 577);

    auto ta_zzzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_hi + 578);

    auto ta_zzzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_hi + 579);

    auto ta_zzzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_hi + 580);

    auto ta_zzzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_hi + 581);

    auto ta_zzzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_hi + 582);

    auto ta_zzzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_hi + 583);

    auto ta_zzzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_hi + 584);

    auto ta_zzzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_hi + 585);

    auto ta_zzzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_hi + 586);

    auto ta_zzzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_hi + 587);

    // Set up 0-28 components of targeted buffer : II

    auto ta_xxxxxx_xxxxxx_0 = pbuffer.data(idx_npot_0_ii);

    auto ta_xxxxxx_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 1);

    auto ta_xxxxxx_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 2);

    auto ta_xxxxxx_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 3);

    auto ta_xxxxxx_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 4);

    auto ta_xxxxxx_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 5);

    auto ta_xxxxxx_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 6);

    auto ta_xxxxxx_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 7);

    auto ta_xxxxxx_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 8);

    auto ta_xxxxxx_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 9);

    auto ta_xxxxxx_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 10);

    auto ta_xxxxxx_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 11);

    auto ta_xxxxxx_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 12);

    auto ta_xxxxxx_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 13);

    auto ta_xxxxxx_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 14);

    auto ta_xxxxxx_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 15);

    auto ta_xxxxxx_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 16);

    auto ta_xxxxxx_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 17);

    auto ta_xxxxxx_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 18);

    auto ta_xxxxxx_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 19);

    auto ta_xxxxxx_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 20);

    auto ta_xxxxxx_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 21);

    auto ta_xxxxxx_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 22);

    auto ta_xxxxxx_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 23);

    auto ta_xxxxxx_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 24);

    auto ta_xxxxxx_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 25);

    auto ta_xxxxxx_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 26);

    auto ta_xxxxxx_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 27);

    #pragma omp simd aligned(pa_x, pc_x, ta_xxxx_xxxxxx_0, ta_xxxx_xxxxxx_1, ta_xxxx_xxxxxy_0, ta_xxxx_xxxxxy_1, ta_xxxx_xxxxxz_0, ta_xxxx_xxxxxz_1, ta_xxxx_xxxxyy_0, ta_xxxx_xxxxyy_1, ta_xxxx_xxxxyz_0, ta_xxxx_xxxxyz_1, ta_xxxx_xxxxzz_0, ta_xxxx_xxxxzz_1, ta_xxxx_xxxyyy_0, ta_xxxx_xxxyyy_1, ta_xxxx_xxxyyz_0, ta_xxxx_xxxyyz_1, ta_xxxx_xxxyzz_0, ta_xxxx_xxxyzz_1, ta_xxxx_xxxzzz_0, ta_xxxx_xxxzzz_1, ta_xxxx_xxyyyy_0, ta_xxxx_xxyyyy_1, ta_xxxx_xxyyyz_0, ta_xxxx_xxyyyz_1, ta_xxxx_xxyyzz_0, ta_xxxx_xxyyzz_1, ta_xxxx_xxyzzz_0, ta_xxxx_xxyzzz_1, ta_xxxx_xxzzzz_0, ta_xxxx_xxzzzz_1, ta_xxxx_xyyyyy_0, ta_xxxx_xyyyyy_1, ta_xxxx_xyyyyz_0, ta_xxxx_xyyyyz_1, ta_xxxx_xyyyzz_0, ta_xxxx_xyyyzz_1, ta_xxxx_xyyzzz_0, ta_xxxx_xyyzzz_1, ta_xxxx_xyzzzz_0, ta_xxxx_xyzzzz_1, ta_xxxx_xzzzzz_0, ta_xxxx_xzzzzz_1, ta_xxxx_yyyyyy_0, ta_xxxx_yyyyyy_1, ta_xxxx_yyyyyz_0, ta_xxxx_yyyyyz_1, ta_xxxx_yyyyzz_0, ta_xxxx_yyyyzz_1, ta_xxxx_yyyzzz_0, ta_xxxx_yyyzzz_1, ta_xxxx_yyzzzz_0, ta_xxxx_yyzzzz_1, ta_xxxx_yzzzzz_0, ta_xxxx_yzzzzz_1, ta_xxxx_zzzzzz_0, ta_xxxx_zzzzzz_1, ta_xxxxx_xxxxx_0, ta_xxxxx_xxxxx_1, ta_xxxxx_xxxxxx_0, ta_xxxxx_xxxxxx_1, ta_xxxxx_xxxxxy_0, ta_xxxxx_xxxxxy_1, ta_xxxxx_xxxxxz_0, ta_xxxxx_xxxxxz_1, ta_xxxxx_xxxxy_0, ta_xxxxx_xxxxy_1, ta_xxxxx_xxxxyy_0, ta_xxxxx_xxxxyy_1, ta_xxxxx_xxxxyz_0, ta_xxxxx_xxxxyz_1, ta_xxxxx_xxxxz_0, ta_xxxxx_xxxxz_1, ta_xxxxx_xxxxzz_0, ta_xxxxx_xxxxzz_1, ta_xxxxx_xxxyy_0, ta_xxxxx_xxxyy_1, ta_xxxxx_xxxyyy_0, ta_xxxxx_xxxyyy_1, ta_xxxxx_xxxyyz_0, ta_xxxxx_xxxyyz_1, ta_xxxxx_xxxyz_0, ta_xxxxx_xxxyz_1, ta_xxxxx_xxxyzz_0, ta_xxxxx_xxxyzz_1, ta_xxxxx_xxxzz_0, ta_xxxxx_xxxzz_1, ta_xxxxx_xxxzzz_0, ta_xxxxx_xxxzzz_1, ta_xxxxx_xxyyy_0, ta_xxxxx_xxyyy_1, ta_xxxxx_xxyyyy_0, ta_xxxxx_xxyyyy_1, ta_xxxxx_xxyyyz_0, ta_xxxxx_xxyyyz_1, ta_xxxxx_xxyyz_0, ta_xxxxx_xxyyz_1, ta_xxxxx_xxyyzz_0, ta_xxxxx_xxyyzz_1, ta_xxxxx_xxyzz_0, ta_xxxxx_xxyzz_1, ta_xxxxx_xxyzzz_0, ta_xxxxx_xxyzzz_1, ta_xxxxx_xxzzz_0, ta_xxxxx_xxzzz_1, ta_xxxxx_xxzzzz_0, ta_xxxxx_xxzzzz_1, ta_xxxxx_xyyyy_0, ta_xxxxx_xyyyy_1, ta_xxxxx_xyyyyy_0, ta_xxxxx_xyyyyy_1, ta_xxxxx_xyyyyz_0, ta_xxxxx_xyyyyz_1, ta_xxxxx_xyyyz_0, ta_xxxxx_xyyyz_1, ta_xxxxx_xyyyzz_0, ta_xxxxx_xyyyzz_1, ta_xxxxx_xyyzz_0, ta_xxxxx_xyyzz_1, ta_xxxxx_xyyzzz_0, ta_xxxxx_xyyzzz_1, ta_xxxxx_xyzzz_0, ta_xxxxx_xyzzz_1, ta_xxxxx_xyzzzz_0, ta_xxxxx_xyzzzz_1, ta_xxxxx_xzzzz_0, ta_xxxxx_xzzzz_1, ta_xxxxx_xzzzzz_0, ta_xxxxx_xzzzzz_1, ta_xxxxx_yyyyy_0, ta_xxxxx_yyyyy_1, ta_xxxxx_yyyyyy_0, ta_xxxxx_yyyyyy_1, ta_xxxxx_yyyyyz_0, ta_xxxxx_yyyyyz_1, ta_xxxxx_yyyyz_0, ta_xxxxx_yyyyz_1, ta_xxxxx_yyyyzz_0, ta_xxxxx_yyyyzz_1, ta_xxxxx_yyyzz_0, ta_xxxxx_yyyzz_1, ta_xxxxx_yyyzzz_0, ta_xxxxx_yyyzzz_1, ta_xxxxx_yyzzz_0, ta_xxxxx_yyzzz_1, ta_xxxxx_yyzzzz_0, ta_xxxxx_yyzzzz_1, ta_xxxxx_yzzzz_0, ta_xxxxx_yzzzz_1, ta_xxxxx_yzzzzz_0, ta_xxxxx_yzzzzz_1, ta_xxxxx_zzzzz_0, ta_xxxxx_zzzzz_1, ta_xxxxx_zzzzzz_0, ta_xxxxx_zzzzzz_1, ta_xxxxxx_xxxxxx_0, ta_xxxxxx_xxxxxy_0, ta_xxxxxx_xxxxxz_0, ta_xxxxxx_xxxxyy_0, ta_xxxxxx_xxxxyz_0, ta_xxxxxx_xxxxzz_0, ta_xxxxxx_xxxyyy_0, ta_xxxxxx_xxxyyz_0, ta_xxxxxx_xxxyzz_0, ta_xxxxxx_xxxzzz_0, ta_xxxxxx_xxyyyy_0, ta_xxxxxx_xxyyyz_0, ta_xxxxxx_xxyyzz_0, ta_xxxxxx_xxyzzz_0, ta_xxxxxx_xxzzzz_0, ta_xxxxxx_xyyyyy_0, ta_xxxxxx_xyyyyz_0, ta_xxxxxx_xyyyzz_0, ta_xxxxxx_xyyzzz_0, ta_xxxxxx_xyzzzz_0, ta_xxxxxx_xzzzzz_0, ta_xxxxxx_yyyyyy_0, ta_xxxxxx_yyyyyz_0, ta_xxxxxx_yyyyzz_0, ta_xxxxxx_yyyzzz_0, ta_xxxxxx_yyzzzz_0, ta_xxxxxx_yzzzzz_0, ta_xxxxxx_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_xxxxxx_0[i] = 5.0 * ta_xxxx_xxxxxx_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxxx_1[i] * fe_0 + 6.0 * ta_xxxxx_xxxxx_0[i] * fe_0 - 6.0 * ta_xxxxx_xxxxx_1[i] * fe_0 + ta_xxxxx_xxxxxx_0[i] * pa_x[i] - ta_xxxxx_xxxxxx_1[i] * pc_x[i];

        ta_xxxxxx_xxxxxy_0[i] = 5.0 * ta_xxxx_xxxxxy_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxxy_1[i] * fe_0 + 5.0 * ta_xxxxx_xxxxy_0[i] * fe_0 - 5.0 * ta_xxxxx_xxxxy_1[i] * fe_0 + ta_xxxxx_xxxxxy_0[i] * pa_x[i] - ta_xxxxx_xxxxxy_1[i] * pc_x[i];

        ta_xxxxxx_xxxxxz_0[i] = 5.0 * ta_xxxx_xxxxxz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxxz_1[i] * fe_0 + 5.0 * ta_xxxxx_xxxxz_0[i] * fe_0 - 5.0 * ta_xxxxx_xxxxz_1[i] * fe_0 + ta_xxxxx_xxxxxz_0[i] * pa_x[i] - ta_xxxxx_xxxxxz_1[i] * pc_x[i];

        ta_xxxxxx_xxxxyy_0[i] = 5.0 * ta_xxxx_xxxxyy_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxyy_1[i] * fe_0 + 4.0 * ta_xxxxx_xxxyy_0[i] * fe_0 - 4.0 * ta_xxxxx_xxxyy_1[i] * fe_0 + ta_xxxxx_xxxxyy_0[i] * pa_x[i] - ta_xxxxx_xxxxyy_1[i] * pc_x[i];

        ta_xxxxxx_xxxxyz_0[i] = 5.0 * ta_xxxx_xxxxyz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxxxx_xxxyz_0[i] * fe_0 - 4.0 * ta_xxxxx_xxxyz_1[i] * fe_0 + ta_xxxxx_xxxxyz_0[i] * pa_x[i] - ta_xxxxx_xxxxyz_1[i] * pc_x[i];

        ta_xxxxxx_xxxxzz_0[i] = 5.0 * ta_xxxx_xxxxzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxzz_1[i] * fe_0 + 4.0 * ta_xxxxx_xxxzz_0[i] * fe_0 - 4.0 * ta_xxxxx_xxxzz_1[i] * fe_0 + ta_xxxxx_xxxxzz_0[i] * pa_x[i] - ta_xxxxx_xxxxzz_1[i] * pc_x[i];

        ta_xxxxxx_xxxyyy_0[i] = 5.0 * ta_xxxx_xxxyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xxxyyy_1[i] * fe_0 + 3.0 * ta_xxxxx_xxyyy_0[i] * fe_0 - 3.0 * ta_xxxxx_xxyyy_1[i] * fe_0 + ta_xxxxx_xxxyyy_0[i] * pa_x[i] - ta_xxxxx_xxxyyy_1[i] * pc_x[i];

        ta_xxxxxx_xxxyyz_0[i] = 5.0 * ta_xxxx_xxxyyz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxxxx_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxyyz_1[i] * fe_0 + ta_xxxxx_xxxyyz_0[i] * pa_x[i] - ta_xxxxx_xxxyyz_1[i] * pc_x[i];

        ta_xxxxxx_xxxyzz_0[i] = 5.0 * ta_xxxx_xxxyzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxxxx_xxyzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxyzz_1[i] * fe_0 + ta_xxxxx_xxxyzz_0[i] * pa_x[i] - ta_xxxxx_xxxyzz_1[i] * pc_x[i];

        ta_xxxxxx_xxxzzz_0[i] = 5.0 * ta_xxxx_xxxzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxzzz_1[i] * fe_0 + 3.0 * ta_xxxxx_xxzzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxzzz_1[i] * fe_0 + ta_xxxxx_xxxzzz_0[i] * pa_x[i] - ta_xxxxx_xxxzzz_1[i] * pc_x[i];

        ta_xxxxxx_xxyyyy_0[i] = 5.0 * ta_xxxx_xxyyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xxyyyy_1[i] * fe_0 + 2.0 * ta_xxxxx_xyyyy_0[i] * fe_0 - 2.0 * ta_xxxxx_xyyyy_1[i] * fe_0 + ta_xxxxx_xxyyyy_0[i] * pa_x[i] - ta_xxxxx_xxyyyy_1[i] * pc_x[i];

        ta_xxxxxx_xxyyyz_0[i] = 5.0 * ta_xxxx_xxyyyz_0[i] * fe_0 - 5.0 * ta_xxxx_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxxxx_xyyyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyyyz_1[i] * fe_0 + ta_xxxxx_xxyyyz_0[i] * pa_x[i] - ta_xxxxx_xxyyyz_1[i] * pc_x[i];

        ta_xxxxxx_xxyyzz_0[i] = 5.0 * ta_xxxx_xxyyzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxxxx_xyyzz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyyzz_1[i] * fe_0 + ta_xxxxx_xxyyzz_0[i] * pa_x[i] - ta_xxxxx_xxyyzz_1[i] * pc_x[i];

        ta_xxxxxx_xxyzzz_0[i] = 5.0 * ta_xxxx_xxyzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxxxx_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyzzz_1[i] * fe_0 + ta_xxxxx_xxyzzz_0[i] * pa_x[i] - ta_xxxxx_xxyzzz_1[i] * pc_x[i];

        ta_xxxxxx_xxzzzz_0[i] = 5.0 * ta_xxxx_xxzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxzzzz_1[i] * fe_0 + 2.0 * ta_xxxxx_xzzzz_0[i] * fe_0 - 2.0 * ta_xxxxx_xzzzz_1[i] * fe_0 + ta_xxxxx_xxzzzz_0[i] * pa_x[i] - ta_xxxxx_xxzzzz_1[i] * pc_x[i];

        ta_xxxxxx_xyyyyy_0[i] = 5.0 * ta_xxxx_xyyyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xyyyyy_1[i] * fe_0 + ta_xxxxx_yyyyy_0[i] * fe_0 - ta_xxxxx_yyyyy_1[i] * fe_0 + ta_xxxxx_xyyyyy_0[i] * pa_x[i] - ta_xxxxx_xyyyyy_1[i] * pc_x[i];

        ta_xxxxxx_xyyyyz_0[i] = 5.0 * ta_xxxx_xyyyyz_0[i] * fe_0 - 5.0 * ta_xxxx_xyyyyz_1[i] * fe_0 + ta_xxxxx_yyyyz_0[i] * fe_0 - ta_xxxxx_yyyyz_1[i] * fe_0 + ta_xxxxx_xyyyyz_0[i] * pa_x[i] - ta_xxxxx_xyyyyz_1[i] * pc_x[i];

        ta_xxxxxx_xyyyzz_0[i] = 5.0 * ta_xxxx_xyyyzz_0[i] * fe_0 - 5.0 * ta_xxxx_xyyyzz_1[i] * fe_0 + ta_xxxxx_yyyzz_0[i] * fe_0 - ta_xxxxx_yyyzz_1[i] * fe_0 + ta_xxxxx_xyyyzz_0[i] * pa_x[i] - ta_xxxxx_xyyyzz_1[i] * pc_x[i];

        ta_xxxxxx_xyyzzz_0[i] = 5.0 * ta_xxxx_xyyzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xyyzzz_1[i] * fe_0 + ta_xxxxx_yyzzz_0[i] * fe_0 - ta_xxxxx_yyzzz_1[i] * fe_0 + ta_xxxxx_xyyzzz_0[i] * pa_x[i] - ta_xxxxx_xyyzzz_1[i] * pc_x[i];

        ta_xxxxxx_xyzzzz_0[i] = 5.0 * ta_xxxx_xyzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xyzzzz_1[i] * fe_0 + ta_xxxxx_yzzzz_0[i] * fe_0 - ta_xxxxx_yzzzz_1[i] * fe_0 + ta_xxxxx_xyzzzz_0[i] * pa_x[i] - ta_xxxxx_xyzzzz_1[i] * pc_x[i];

        ta_xxxxxx_xzzzzz_0[i] = 5.0 * ta_xxxx_xzzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xzzzzz_1[i] * fe_0 + ta_xxxxx_zzzzz_0[i] * fe_0 - ta_xxxxx_zzzzz_1[i] * fe_0 + ta_xxxxx_xzzzzz_0[i] * pa_x[i] - ta_xxxxx_xzzzzz_1[i] * pc_x[i];

        ta_xxxxxx_yyyyyy_0[i] = 5.0 * ta_xxxx_yyyyyy_0[i] * fe_0 - 5.0 * ta_xxxx_yyyyyy_1[i] * fe_0 + ta_xxxxx_yyyyyy_0[i] * pa_x[i] - ta_xxxxx_yyyyyy_1[i] * pc_x[i];

        ta_xxxxxx_yyyyyz_0[i] = 5.0 * ta_xxxx_yyyyyz_0[i] * fe_0 - 5.0 * ta_xxxx_yyyyyz_1[i] * fe_0 + ta_xxxxx_yyyyyz_0[i] * pa_x[i] - ta_xxxxx_yyyyyz_1[i] * pc_x[i];

        ta_xxxxxx_yyyyzz_0[i] = 5.0 * ta_xxxx_yyyyzz_0[i] * fe_0 - 5.0 * ta_xxxx_yyyyzz_1[i] * fe_0 + ta_xxxxx_yyyyzz_0[i] * pa_x[i] - ta_xxxxx_yyyyzz_1[i] * pc_x[i];

        ta_xxxxxx_yyyzzz_0[i] = 5.0 * ta_xxxx_yyyzzz_0[i] * fe_0 - 5.0 * ta_xxxx_yyyzzz_1[i] * fe_0 + ta_xxxxx_yyyzzz_0[i] * pa_x[i] - ta_xxxxx_yyyzzz_1[i] * pc_x[i];

        ta_xxxxxx_yyzzzz_0[i] = 5.0 * ta_xxxx_yyzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_yyzzzz_1[i] * fe_0 + ta_xxxxx_yyzzzz_0[i] * pa_x[i] - ta_xxxxx_yyzzzz_1[i] * pc_x[i];

        ta_xxxxxx_yzzzzz_0[i] = 5.0 * ta_xxxx_yzzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_yzzzzz_1[i] * fe_0 + ta_xxxxx_yzzzzz_0[i] * pa_x[i] - ta_xxxxx_yzzzzz_1[i] * pc_x[i];

        ta_xxxxxx_zzzzzz_0[i] = 5.0 * ta_xxxx_zzzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_zzzzzz_1[i] * fe_0 + ta_xxxxx_zzzzzz_0[i] * pa_x[i] - ta_xxxxx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : II

    auto ta_xxxxxy_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 28);

    auto ta_xxxxxy_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 29);

    auto ta_xxxxxy_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 30);

    auto ta_xxxxxy_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 31);

    auto ta_xxxxxy_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 32);

    auto ta_xxxxxy_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 33);

    auto ta_xxxxxy_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 34);

    auto ta_xxxxxy_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 35);

    auto ta_xxxxxy_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 36);

    auto ta_xxxxxy_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 37);

    auto ta_xxxxxy_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 38);

    auto ta_xxxxxy_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 39);

    auto ta_xxxxxy_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 40);

    auto ta_xxxxxy_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 41);

    auto ta_xxxxxy_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 42);

    auto ta_xxxxxy_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 43);

    auto ta_xxxxxy_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 44);

    auto ta_xxxxxy_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 45);

    auto ta_xxxxxy_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 46);

    auto ta_xxxxxy_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 47);

    auto ta_xxxxxy_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 48);

    auto ta_xxxxxy_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 49);

    auto ta_xxxxxy_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 50);

    auto ta_xxxxxy_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 51);

    auto ta_xxxxxy_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 52);

    auto ta_xxxxxy_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 53);

    auto ta_xxxxxy_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 54);

    auto ta_xxxxxy_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 55);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxxxx_xxxxx_0, ta_xxxxx_xxxxx_1, ta_xxxxx_xxxxxx_0, ta_xxxxx_xxxxxx_1, ta_xxxxx_xxxxxy_0, ta_xxxxx_xxxxxy_1, ta_xxxxx_xxxxxz_0, ta_xxxxx_xxxxxz_1, ta_xxxxx_xxxxy_0, ta_xxxxx_xxxxy_1, ta_xxxxx_xxxxyy_0, ta_xxxxx_xxxxyy_1, ta_xxxxx_xxxxyz_0, ta_xxxxx_xxxxyz_1, ta_xxxxx_xxxxz_0, ta_xxxxx_xxxxz_1, ta_xxxxx_xxxxzz_0, ta_xxxxx_xxxxzz_1, ta_xxxxx_xxxyy_0, ta_xxxxx_xxxyy_1, ta_xxxxx_xxxyyy_0, ta_xxxxx_xxxyyy_1, ta_xxxxx_xxxyyz_0, ta_xxxxx_xxxyyz_1, ta_xxxxx_xxxyz_0, ta_xxxxx_xxxyz_1, ta_xxxxx_xxxyzz_0, ta_xxxxx_xxxyzz_1, ta_xxxxx_xxxzz_0, ta_xxxxx_xxxzz_1, ta_xxxxx_xxxzzz_0, ta_xxxxx_xxxzzz_1, ta_xxxxx_xxyyy_0, ta_xxxxx_xxyyy_1, ta_xxxxx_xxyyyy_0, ta_xxxxx_xxyyyy_1, ta_xxxxx_xxyyyz_0, ta_xxxxx_xxyyyz_1, ta_xxxxx_xxyyz_0, ta_xxxxx_xxyyz_1, ta_xxxxx_xxyyzz_0, ta_xxxxx_xxyyzz_1, ta_xxxxx_xxyzz_0, ta_xxxxx_xxyzz_1, ta_xxxxx_xxyzzz_0, ta_xxxxx_xxyzzz_1, ta_xxxxx_xxzzz_0, ta_xxxxx_xxzzz_1, ta_xxxxx_xxzzzz_0, ta_xxxxx_xxzzzz_1, ta_xxxxx_xyyyy_0, ta_xxxxx_xyyyy_1, ta_xxxxx_xyyyyy_0, ta_xxxxx_xyyyyy_1, ta_xxxxx_xyyyyz_0, ta_xxxxx_xyyyyz_1, ta_xxxxx_xyyyz_0, ta_xxxxx_xyyyz_1, ta_xxxxx_xyyyzz_0, ta_xxxxx_xyyyzz_1, ta_xxxxx_xyyzz_0, ta_xxxxx_xyyzz_1, ta_xxxxx_xyyzzz_0, ta_xxxxx_xyyzzz_1, ta_xxxxx_xyzzz_0, ta_xxxxx_xyzzz_1, ta_xxxxx_xyzzzz_0, ta_xxxxx_xyzzzz_1, ta_xxxxx_xzzzz_0, ta_xxxxx_xzzzz_1, ta_xxxxx_xzzzzz_0, ta_xxxxx_xzzzzz_1, ta_xxxxx_zzzzzz_0, ta_xxxxx_zzzzzz_1, ta_xxxxxy_xxxxxx_0, ta_xxxxxy_xxxxxy_0, ta_xxxxxy_xxxxxz_0, ta_xxxxxy_xxxxyy_0, ta_xxxxxy_xxxxyz_0, ta_xxxxxy_xxxxzz_0, ta_xxxxxy_xxxyyy_0, ta_xxxxxy_xxxyyz_0, ta_xxxxxy_xxxyzz_0, ta_xxxxxy_xxxzzz_0, ta_xxxxxy_xxyyyy_0, ta_xxxxxy_xxyyyz_0, ta_xxxxxy_xxyyzz_0, ta_xxxxxy_xxyzzz_0, ta_xxxxxy_xxzzzz_0, ta_xxxxxy_xyyyyy_0, ta_xxxxxy_xyyyyz_0, ta_xxxxxy_xyyyzz_0, ta_xxxxxy_xyyzzz_0, ta_xxxxxy_xyzzzz_0, ta_xxxxxy_xzzzzz_0, ta_xxxxxy_yyyyyy_0, ta_xxxxxy_yyyyyz_0, ta_xxxxxy_yyyyzz_0, ta_xxxxxy_yyyzzz_0, ta_xxxxxy_yyzzzz_0, ta_xxxxxy_yzzzzz_0, ta_xxxxxy_zzzzzz_0, ta_xxxxy_yyyyyy_0, ta_xxxxy_yyyyyy_1, ta_xxxxy_yyyyyz_0, ta_xxxxy_yyyyyz_1, ta_xxxxy_yyyyzz_0, ta_xxxxy_yyyyzz_1, ta_xxxxy_yyyzzz_0, ta_xxxxy_yyyzzz_1, ta_xxxxy_yyzzzz_0, ta_xxxxy_yyzzzz_1, ta_xxxxy_yzzzzz_0, ta_xxxxy_yzzzzz_1, ta_xxxy_yyyyyy_0, ta_xxxy_yyyyyy_1, ta_xxxy_yyyyyz_0, ta_xxxy_yyyyyz_1, ta_xxxy_yyyyzz_0, ta_xxxy_yyyyzz_1, ta_xxxy_yyyzzz_0, ta_xxxy_yyyzzz_1, ta_xxxy_yyzzzz_0, ta_xxxy_yyzzzz_1, ta_xxxy_yzzzzz_0, ta_xxxy_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxy_xxxxxx_0[i] = ta_xxxxx_xxxxxx_0[i] * pa_y[i] - ta_xxxxx_xxxxxx_1[i] * pc_y[i];

        ta_xxxxxy_xxxxxy_0[i] = ta_xxxxx_xxxxx_0[i] * fe_0 - ta_xxxxx_xxxxx_1[i] * fe_0 + ta_xxxxx_xxxxxy_0[i] * pa_y[i] - ta_xxxxx_xxxxxy_1[i] * pc_y[i];

        ta_xxxxxy_xxxxxz_0[i] = ta_xxxxx_xxxxxz_0[i] * pa_y[i] - ta_xxxxx_xxxxxz_1[i] * pc_y[i];

        ta_xxxxxy_xxxxyy_0[i] = 2.0 * ta_xxxxx_xxxxy_0[i] * fe_0 - 2.0 * ta_xxxxx_xxxxy_1[i] * fe_0 + ta_xxxxx_xxxxyy_0[i] * pa_y[i] - ta_xxxxx_xxxxyy_1[i] * pc_y[i];

        ta_xxxxxy_xxxxyz_0[i] = ta_xxxxx_xxxxz_0[i] * fe_0 - ta_xxxxx_xxxxz_1[i] * fe_0 + ta_xxxxx_xxxxyz_0[i] * pa_y[i] - ta_xxxxx_xxxxyz_1[i] * pc_y[i];

        ta_xxxxxy_xxxxzz_0[i] = ta_xxxxx_xxxxzz_0[i] * pa_y[i] - ta_xxxxx_xxxxzz_1[i] * pc_y[i];

        ta_xxxxxy_xxxyyy_0[i] = 3.0 * ta_xxxxx_xxxyy_0[i] * fe_0 - 3.0 * ta_xxxxx_xxxyy_1[i] * fe_0 + ta_xxxxx_xxxyyy_0[i] * pa_y[i] - ta_xxxxx_xxxyyy_1[i] * pc_y[i];

        ta_xxxxxy_xxxyyz_0[i] = 2.0 * ta_xxxxx_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxxyz_1[i] * fe_0 + ta_xxxxx_xxxyyz_0[i] * pa_y[i] - ta_xxxxx_xxxyyz_1[i] * pc_y[i];

        ta_xxxxxy_xxxyzz_0[i] = ta_xxxxx_xxxzz_0[i] * fe_0 - ta_xxxxx_xxxzz_1[i] * fe_0 + ta_xxxxx_xxxyzz_0[i] * pa_y[i] - ta_xxxxx_xxxyzz_1[i] * pc_y[i];

        ta_xxxxxy_xxxzzz_0[i] = ta_xxxxx_xxxzzz_0[i] * pa_y[i] - ta_xxxxx_xxxzzz_1[i] * pc_y[i];

        ta_xxxxxy_xxyyyy_0[i] = 4.0 * ta_xxxxx_xxyyy_0[i] * fe_0 - 4.0 * ta_xxxxx_xxyyy_1[i] * fe_0 + ta_xxxxx_xxyyyy_0[i] * pa_y[i] - ta_xxxxx_xxyyyy_1[i] * pc_y[i];

        ta_xxxxxy_xxyyyz_0[i] = 3.0 * ta_xxxxx_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxyyz_1[i] * fe_0 + ta_xxxxx_xxyyyz_0[i] * pa_y[i] - ta_xxxxx_xxyyyz_1[i] * pc_y[i];

        ta_xxxxxy_xxyyzz_0[i] = 2.0 * ta_xxxxx_xxyzz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxyzz_1[i] * fe_0 + ta_xxxxx_xxyyzz_0[i] * pa_y[i] - ta_xxxxx_xxyyzz_1[i] * pc_y[i];

        ta_xxxxxy_xxyzzz_0[i] = ta_xxxxx_xxzzz_0[i] * fe_0 - ta_xxxxx_xxzzz_1[i] * fe_0 + ta_xxxxx_xxyzzz_0[i] * pa_y[i] - ta_xxxxx_xxyzzz_1[i] * pc_y[i];

        ta_xxxxxy_xxzzzz_0[i] = ta_xxxxx_xxzzzz_0[i] * pa_y[i] - ta_xxxxx_xxzzzz_1[i] * pc_y[i];

        ta_xxxxxy_xyyyyy_0[i] = 5.0 * ta_xxxxx_xyyyy_0[i] * fe_0 - 5.0 * ta_xxxxx_xyyyy_1[i] * fe_0 + ta_xxxxx_xyyyyy_0[i] * pa_y[i] - ta_xxxxx_xyyyyy_1[i] * pc_y[i];

        ta_xxxxxy_xyyyyz_0[i] = 4.0 * ta_xxxxx_xyyyz_0[i] * fe_0 - 4.0 * ta_xxxxx_xyyyz_1[i] * fe_0 + ta_xxxxx_xyyyyz_0[i] * pa_y[i] - ta_xxxxx_xyyyyz_1[i] * pc_y[i];

        ta_xxxxxy_xyyyzz_0[i] = 3.0 * ta_xxxxx_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xyyzz_1[i] * fe_0 + ta_xxxxx_xyyyzz_0[i] * pa_y[i] - ta_xxxxx_xyyyzz_1[i] * pc_y[i];

        ta_xxxxxy_xyyzzz_0[i] = 2.0 * ta_xxxxx_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyzzz_1[i] * fe_0 + ta_xxxxx_xyyzzz_0[i] * pa_y[i] - ta_xxxxx_xyyzzz_1[i] * pc_y[i];

        ta_xxxxxy_xyzzzz_0[i] = ta_xxxxx_xzzzz_0[i] * fe_0 - ta_xxxxx_xzzzz_1[i] * fe_0 + ta_xxxxx_xyzzzz_0[i] * pa_y[i] - ta_xxxxx_xyzzzz_1[i] * pc_y[i];

        ta_xxxxxy_xzzzzz_0[i] = ta_xxxxx_xzzzzz_0[i] * pa_y[i] - ta_xxxxx_xzzzzz_1[i] * pc_y[i];

        ta_xxxxxy_yyyyyy_0[i] = 4.0 * ta_xxxy_yyyyyy_0[i] * fe_0 - 4.0 * ta_xxxy_yyyyyy_1[i] * fe_0 + ta_xxxxy_yyyyyy_0[i] * pa_x[i] - ta_xxxxy_yyyyyy_1[i] * pc_x[i];

        ta_xxxxxy_yyyyyz_0[i] = 4.0 * ta_xxxy_yyyyyz_0[i] * fe_0 - 4.0 * ta_xxxy_yyyyyz_1[i] * fe_0 + ta_xxxxy_yyyyyz_0[i] * pa_x[i] - ta_xxxxy_yyyyyz_1[i] * pc_x[i];

        ta_xxxxxy_yyyyzz_0[i] = 4.0 * ta_xxxy_yyyyzz_0[i] * fe_0 - 4.0 * ta_xxxy_yyyyzz_1[i] * fe_0 + ta_xxxxy_yyyyzz_0[i] * pa_x[i] - ta_xxxxy_yyyyzz_1[i] * pc_x[i];

        ta_xxxxxy_yyyzzz_0[i] = 4.0 * ta_xxxy_yyyzzz_0[i] * fe_0 - 4.0 * ta_xxxy_yyyzzz_1[i] * fe_0 + ta_xxxxy_yyyzzz_0[i] * pa_x[i] - ta_xxxxy_yyyzzz_1[i] * pc_x[i];

        ta_xxxxxy_yyzzzz_0[i] = 4.0 * ta_xxxy_yyzzzz_0[i] * fe_0 - 4.0 * ta_xxxy_yyzzzz_1[i] * fe_0 + ta_xxxxy_yyzzzz_0[i] * pa_x[i] - ta_xxxxy_yyzzzz_1[i] * pc_x[i];

        ta_xxxxxy_yzzzzz_0[i] = 4.0 * ta_xxxy_yzzzzz_0[i] * fe_0 - 4.0 * ta_xxxy_yzzzzz_1[i] * fe_0 + ta_xxxxy_yzzzzz_0[i] * pa_x[i] - ta_xxxxy_yzzzzz_1[i] * pc_x[i];

        ta_xxxxxy_zzzzzz_0[i] = ta_xxxxx_zzzzzz_0[i] * pa_y[i] - ta_xxxxx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : II

    auto ta_xxxxxz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 56);

    auto ta_xxxxxz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 57);

    auto ta_xxxxxz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 58);

    auto ta_xxxxxz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 59);

    auto ta_xxxxxz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 60);

    auto ta_xxxxxz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 61);

    auto ta_xxxxxz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 62);

    auto ta_xxxxxz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 63);

    auto ta_xxxxxz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 64);

    auto ta_xxxxxz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 65);

    auto ta_xxxxxz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 66);

    auto ta_xxxxxz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 67);

    auto ta_xxxxxz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 68);

    auto ta_xxxxxz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 69);

    auto ta_xxxxxz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 70);

    auto ta_xxxxxz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 71);

    auto ta_xxxxxz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 72);

    auto ta_xxxxxz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 73);

    auto ta_xxxxxz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 74);

    auto ta_xxxxxz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 75);

    auto ta_xxxxxz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 76);

    auto ta_xxxxxz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 77);

    auto ta_xxxxxz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 78);

    auto ta_xxxxxz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 79);

    auto ta_xxxxxz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 80);

    auto ta_xxxxxz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 81);

    auto ta_xxxxxz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 82);

    auto ta_xxxxxz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 83);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xxxxx_xxxxx_0, ta_xxxxx_xxxxx_1, ta_xxxxx_xxxxxx_0, ta_xxxxx_xxxxxx_1, ta_xxxxx_xxxxxy_0, ta_xxxxx_xxxxxy_1, ta_xxxxx_xxxxxz_0, ta_xxxxx_xxxxxz_1, ta_xxxxx_xxxxy_0, ta_xxxxx_xxxxy_1, ta_xxxxx_xxxxyy_0, ta_xxxxx_xxxxyy_1, ta_xxxxx_xxxxyz_0, ta_xxxxx_xxxxyz_1, ta_xxxxx_xxxxz_0, ta_xxxxx_xxxxz_1, ta_xxxxx_xxxxzz_0, ta_xxxxx_xxxxzz_1, ta_xxxxx_xxxyy_0, ta_xxxxx_xxxyy_1, ta_xxxxx_xxxyyy_0, ta_xxxxx_xxxyyy_1, ta_xxxxx_xxxyyz_0, ta_xxxxx_xxxyyz_1, ta_xxxxx_xxxyz_0, ta_xxxxx_xxxyz_1, ta_xxxxx_xxxyzz_0, ta_xxxxx_xxxyzz_1, ta_xxxxx_xxxzz_0, ta_xxxxx_xxxzz_1, ta_xxxxx_xxxzzz_0, ta_xxxxx_xxxzzz_1, ta_xxxxx_xxyyy_0, ta_xxxxx_xxyyy_1, ta_xxxxx_xxyyyy_0, ta_xxxxx_xxyyyy_1, ta_xxxxx_xxyyyz_0, ta_xxxxx_xxyyyz_1, ta_xxxxx_xxyyz_0, ta_xxxxx_xxyyz_1, ta_xxxxx_xxyyzz_0, ta_xxxxx_xxyyzz_1, ta_xxxxx_xxyzz_0, ta_xxxxx_xxyzz_1, ta_xxxxx_xxyzzz_0, ta_xxxxx_xxyzzz_1, ta_xxxxx_xxzzz_0, ta_xxxxx_xxzzz_1, ta_xxxxx_xxzzzz_0, ta_xxxxx_xxzzzz_1, ta_xxxxx_xyyyy_0, ta_xxxxx_xyyyy_1, ta_xxxxx_xyyyyy_0, ta_xxxxx_xyyyyy_1, ta_xxxxx_xyyyyz_0, ta_xxxxx_xyyyyz_1, ta_xxxxx_xyyyz_0, ta_xxxxx_xyyyz_1, ta_xxxxx_xyyyzz_0, ta_xxxxx_xyyyzz_1, ta_xxxxx_xyyzz_0, ta_xxxxx_xyyzz_1, ta_xxxxx_xyyzzz_0, ta_xxxxx_xyyzzz_1, ta_xxxxx_xyzzz_0, ta_xxxxx_xyzzz_1, ta_xxxxx_xyzzzz_0, ta_xxxxx_xyzzzz_1, ta_xxxxx_xzzzz_0, ta_xxxxx_xzzzz_1, ta_xxxxx_xzzzzz_0, ta_xxxxx_xzzzzz_1, ta_xxxxx_yyyyyy_0, ta_xxxxx_yyyyyy_1, ta_xxxxxz_xxxxxx_0, ta_xxxxxz_xxxxxy_0, ta_xxxxxz_xxxxxz_0, ta_xxxxxz_xxxxyy_0, ta_xxxxxz_xxxxyz_0, ta_xxxxxz_xxxxzz_0, ta_xxxxxz_xxxyyy_0, ta_xxxxxz_xxxyyz_0, ta_xxxxxz_xxxyzz_0, ta_xxxxxz_xxxzzz_0, ta_xxxxxz_xxyyyy_0, ta_xxxxxz_xxyyyz_0, ta_xxxxxz_xxyyzz_0, ta_xxxxxz_xxyzzz_0, ta_xxxxxz_xxzzzz_0, ta_xxxxxz_xyyyyy_0, ta_xxxxxz_xyyyyz_0, ta_xxxxxz_xyyyzz_0, ta_xxxxxz_xyyzzz_0, ta_xxxxxz_xyzzzz_0, ta_xxxxxz_xzzzzz_0, ta_xxxxxz_yyyyyy_0, ta_xxxxxz_yyyyyz_0, ta_xxxxxz_yyyyzz_0, ta_xxxxxz_yyyzzz_0, ta_xxxxxz_yyzzzz_0, ta_xxxxxz_yzzzzz_0, ta_xxxxxz_zzzzzz_0, ta_xxxxz_yyyyyz_0, ta_xxxxz_yyyyyz_1, ta_xxxxz_yyyyzz_0, ta_xxxxz_yyyyzz_1, ta_xxxxz_yyyzzz_0, ta_xxxxz_yyyzzz_1, ta_xxxxz_yyzzzz_0, ta_xxxxz_yyzzzz_1, ta_xxxxz_yzzzzz_0, ta_xxxxz_yzzzzz_1, ta_xxxxz_zzzzzz_0, ta_xxxxz_zzzzzz_1, ta_xxxz_yyyyyz_0, ta_xxxz_yyyyyz_1, ta_xxxz_yyyyzz_0, ta_xxxz_yyyyzz_1, ta_xxxz_yyyzzz_0, ta_xxxz_yyyzzz_1, ta_xxxz_yyzzzz_0, ta_xxxz_yyzzzz_1, ta_xxxz_yzzzzz_0, ta_xxxz_yzzzzz_1, ta_xxxz_zzzzzz_0, ta_xxxz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxz_xxxxxx_0[i] = ta_xxxxx_xxxxxx_0[i] * pa_z[i] - ta_xxxxx_xxxxxx_1[i] * pc_z[i];

        ta_xxxxxz_xxxxxy_0[i] = ta_xxxxx_xxxxxy_0[i] * pa_z[i] - ta_xxxxx_xxxxxy_1[i] * pc_z[i];

        ta_xxxxxz_xxxxxz_0[i] = ta_xxxxx_xxxxx_0[i] * fe_0 - ta_xxxxx_xxxxx_1[i] * fe_0 + ta_xxxxx_xxxxxz_0[i] * pa_z[i] - ta_xxxxx_xxxxxz_1[i] * pc_z[i];

        ta_xxxxxz_xxxxyy_0[i] = ta_xxxxx_xxxxyy_0[i] * pa_z[i] - ta_xxxxx_xxxxyy_1[i] * pc_z[i];

        ta_xxxxxz_xxxxyz_0[i] = ta_xxxxx_xxxxy_0[i] * fe_0 - ta_xxxxx_xxxxy_1[i] * fe_0 + ta_xxxxx_xxxxyz_0[i] * pa_z[i] - ta_xxxxx_xxxxyz_1[i] * pc_z[i];

        ta_xxxxxz_xxxxzz_0[i] = 2.0 * ta_xxxxx_xxxxz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxxxz_1[i] * fe_0 + ta_xxxxx_xxxxzz_0[i] * pa_z[i] - ta_xxxxx_xxxxzz_1[i] * pc_z[i];

        ta_xxxxxz_xxxyyy_0[i] = ta_xxxxx_xxxyyy_0[i] * pa_z[i] - ta_xxxxx_xxxyyy_1[i] * pc_z[i];

        ta_xxxxxz_xxxyyz_0[i] = ta_xxxxx_xxxyy_0[i] * fe_0 - ta_xxxxx_xxxyy_1[i] * fe_0 + ta_xxxxx_xxxyyz_0[i] * pa_z[i] - ta_xxxxx_xxxyyz_1[i] * pc_z[i];

        ta_xxxxxz_xxxyzz_0[i] = 2.0 * ta_xxxxx_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxxyz_1[i] * fe_0 + ta_xxxxx_xxxyzz_0[i] * pa_z[i] - ta_xxxxx_xxxyzz_1[i] * pc_z[i];

        ta_xxxxxz_xxxzzz_0[i] = 3.0 * ta_xxxxx_xxxzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxxzz_1[i] * fe_0 + ta_xxxxx_xxxzzz_0[i] * pa_z[i] - ta_xxxxx_xxxzzz_1[i] * pc_z[i];

        ta_xxxxxz_xxyyyy_0[i] = ta_xxxxx_xxyyyy_0[i] * pa_z[i] - ta_xxxxx_xxyyyy_1[i] * pc_z[i];

        ta_xxxxxz_xxyyyz_0[i] = ta_xxxxx_xxyyy_0[i] * fe_0 - ta_xxxxx_xxyyy_1[i] * fe_0 + ta_xxxxx_xxyyyz_0[i] * pa_z[i] - ta_xxxxx_xxyyyz_1[i] * pc_z[i];

        ta_xxxxxz_xxyyzz_0[i] = 2.0 * ta_xxxxx_xxyyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxyyz_1[i] * fe_0 + ta_xxxxx_xxyyzz_0[i] * pa_z[i] - ta_xxxxx_xxyyzz_1[i] * pc_z[i];

        ta_xxxxxz_xxyzzz_0[i] = 3.0 * ta_xxxxx_xxyzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxyzz_1[i] * fe_0 + ta_xxxxx_xxyzzz_0[i] * pa_z[i] - ta_xxxxx_xxyzzz_1[i] * pc_z[i];

        ta_xxxxxz_xxzzzz_0[i] = 4.0 * ta_xxxxx_xxzzz_0[i] * fe_0 - 4.0 * ta_xxxxx_xxzzz_1[i] * fe_0 + ta_xxxxx_xxzzzz_0[i] * pa_z[i] - ta_xxxxx_xxzzzz_1[i] * pc_z[i];

        ta_xxxxxz_xyyyyy_0[i] = ta_xxxxx_xyyyyy_0[i] * pa_z[i] - ta_xxxxx_xyyyyy_1[i] * pc_z[i];

        ta_xxxxxz_xyyyyz_0[i] = ta_xxxxx_xyyyy_0[i] * fe_0 - ta_xxxxx_xyyyy_1[i] * fe_0 + ta_xxxxx_xyyyyz_0[i] * pa_z[i] - ta_xxxxx_xyyyyz_1[i] * pc_z[i];

        ta_xxxxxz_xyyyzz_0[i] = 2.0 * ta_xxxxx_xyyyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyyyz_1[i] * fe_0 + ta_xxxxx_xyyyzz_0[i] * pa_z[i] - ta_xxxxx_xyyyzz_1[i] * pc_z[i];

        ta_xxxxxz_xyyzzz_0[i] = 3.0 * ta_xxxxx_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xyyzz_1[i] * fe_0 + ta_xxxxx_xyyzzz_0[i] * pa_z[i] - ta_xxxxx_xyyzzz_1[i] * pc_z[i];

        ta_xxxxxz_xyzzzz_0[i] = 4.0 * ta_xxxxx_xyzzz_0[i] * fe_0 - 4.0 * ta_xxxxx_xyzzz_1[i] * fe_0 + ta_xxxxx_xyzzzz_0[i] * pa_z[i] - ta_xxxxx_xyzzzz_1[i] * pc_z[i];

        ta_xxxxxz_xzzzzz_0[i] = 5.0 * ta_xxxxx_xzzzz_0[i] * fe_0 - 5.0 * ta_xxxxx_xzzzz_1[i] * fe_0 + ta_xxxxx_xzzzzz_0[i] * pa_z[i] - ta_xxxxx_xzzzzz_1[i] * pc_z[i];

        ta_xxxxxz_yyyyyy_0[i] = ta_xxxxx_yyyyyy_0[i] * pa_z[i] - ta_xxxxx_yyyyyy_1[i] * pc_z[i];

        ta_xxxxxz_yyyyyz_0[i] = 4.0 * ta_xxxz_yyyyyz_0[i] * fe_0 - 4.0 * ta_xxxz_yyyyyz_1[i] * fe_0 + ta_xxxxz_yyyyyz_0[i] * pa_x[i] - ta_xxxxz_yyyyyz_1[i] * pc_x[i];

        ta_xxxxxz_yyyyzz_0[i] = 4.0 * ta_xxxz_yyyyzz_0[i] * fe_0 - 4.0 * ta_xxxz_yyyyzz_1[i] * fe_0 + ta_xxxxz_yyyyzz_0[i] * pa_x[i] - ta_xxxxz_yyyyzz_1[i] * pc_x[i];

        ta_xxxxxz_yyyzzz_0[i] = 4.0 * ta_xxxz_yyyzzz_0[i] * fe_0 - 4.0 * ta_xxxz_yyyzzz_1[i] * fe_0 + ta_xxxxz_yyyzzz_0[i] * pa_x[i] - ta_xxxxz_yyyzzz_1[i] * pc_x[i];

        ta_xxxxxz_yyzzzz_0[i] = 4.0 * ta_xxxz_yyzzzz_0[i] * fe_0 - 4.0 * ta_xxxz_yyzzzz_1[i] * fe_0 + ta_xxxxz_yyzzzz_0[i] * pa_x[i] - ta_xxxxz_yyzzzz_1[i] * pc_x[i];

        ta_xxxxxz_yzzzzz_0[i] = 4.0 * ta_xxxz_yzzzzz_0[i] * fe_0 - 4.0 * ta_xxxz_yzzzzz_1[i] * fe_0 + ta_xxxxz_yzzzzz_0[i] * pa_x[i] - ta_xxxxz_yzzzzz_1[i] * pc_x[i];

        ta_xxxxxz_zzzzzz_0[i] = 4.0 * ta_xxxz_zzzzzz_0[i] * fe_0 - 4.0 * ta_xxxz_zzzzzz_1[i] * fe_0 + ta_xxxxz_zzzzzz_0[i] * pa_x[i] - ta_xxxxz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 84-112 components of targeted buffer : II

    auto ta_xxxxyy_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 84);

    auto ta_xxxxyy_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 85);

    auto ta_xxxxyy_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 86);

    auto ta_xxxxyy_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 87);

    auto ta_xxxxyy_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 88);

    auto ta_xxxxyy_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 89);

    auto ta_xxxxyy_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 90);

    auto ta_xxxxyy_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 91);

    auto ta_xxxxyy_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 92);

    auto ta_xxxxyy_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 93);

    auto ta_xxxxyy_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 94);

    auto ta_xxxxyy_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 95);

    auto ta_xxxxyy_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 96);

    auto ta_xxxxyy_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 97);

    auto ta_xxxxyy_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 98);

    auto ta_xxxxyy_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 99);

    auto ta_xxxxyy_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 100);

    auto ta_xxxxyy_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 101);

    auto ta_xxxxyy_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 102);

    auto ta_xxxxyy_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 103);

    auto ta_xxxxyy_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 104);

    auto ta_xxxxyy_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 105);

    auto ta_xxxxyy_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 106);

    auto ta_xxxxyy_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 107);

    auto ta_xxxxyy_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 108);

    auto ta_xxxxyy_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 109);

    auto ta_xxxxyy_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 110);

    auto ta_xxxxyy_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 111);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxxx_xxxxxx_0, ta_xxxx_xxxxxx_1, ta_xxxx_xxxxxz_0, ta_xxxx_xxxxxz_1, ta_xxxx_xxxxzz_0, ta_xxxx_xxxxzz_1, ta_xxxx_xxxzzz_0, ta_xxxx_xxxzzz_1, ta_xxxx_xxzzzz_0, ta_xxxx_xxzzzz_1, ta_xxxx_xzzzzz_0, ta_xxxx_xzzzzz_1, ta_xxxxy_xxxxxx_0, ta_xxxxy_xxxxxx_1, ta_xxxxy_xxxxxz_0, ta_xxxxy_xxxxxz_1, ta_xxxxy_xxxxzz_0, ta_xxxxy_xxxxzz_1, ta_xxxxy_xxxzzz_0, ta_xxxxy_xxxzzz_1, ta_xxxxy_xxzzzz_0, ta_xxxxy_xxzzzz_1, ta_xxxxy_xzzzzz_0, ta_xxxxy_xzzzzz_1, ta_xxxxyy_xxxxxx_0, ta_xxxxyy_xxxxxy_0, ta_xxxxyy_xxxxxz_0, ta_xxxxyy_xxxxyy_0, ta_xxxxyy_xxxxyz_0, ta_xxxxyy_xxxxzz_0, ta_xxxxyy_xxxyyy_0, ta_xxxxyy_xxxyyz_0, ta_xxxxyy_xxxyzz_0, ta_xxxxyy_xxxzzz_0, ta_xxxxyy_xxyyyy_0, ta_xxxxyy_xxyyyz_0, ta_xxxxyy_xxyyzz_0, ta_xxxxyy_xxyzzz_0, ta_xxxxyy_xxzzzz_0, ta_xxxxyy_xyyyyy_0, ta_xxxxyy_xyyyyz_0, ta_xxxxyy_xyyyzz_0, ta_xxxxyy_xyyzzz_0, ta_xxxxyy_xyzzzz_0, ta_xxxxyy_xzzzzz_0, ta_xxxxyy_yyyyyy_0, ta_xxxxyy_yyyyyz_0, ta_xxxxyy_yyyyzz_0, ta_xxxxyy_yyyzzz_0, ta_xxxxyy_yyzzzz_0, ta_xxxxyy_yzzzzz_0, ta_xxxxyy_zzzzzz_0, ta_xxxyy_xxxxxy_0, ta_xxxyy_xxxxxy_1, ta_xxxyy_xxxxy_0, ta_xxxyy_xxxxy_1, ta_xxxyy_xxxxyy_0, ta_xxxyy_xxxxyy_1, ta_xxxyy_xxxxyz_0, ta_xxxyy_xxxxyz_1, ta_xxxyy_xxxyy_0, ta_xxxyy_xxxyy_1, ta_xxxyy_xxxyyy_0, ta_xxxyy_xxxyyy_1, ta_xxxyy_xxxyyz_0, ta_xxxyy_xxxyyz_1, ta_xxxyy_xxxyz_0, ta_xxxyy_xxxyz_1, ta_xxxyy_xxxyzz_0, ta_xxxyy_xxxyzz_1, ta_xxxyy_xxyyy_0, ta_xxxyy_xxyyy_1, ta_xxxyy_xxyyyy_0, ta_xxxyy_xxyyyy_1, ta_xxxyy_xxyyyz_0, ta_xxxyy_xxyyyz_1, ta_xxxyy_xxyyz_0, ta_xxxyy_xxyyz_1, ta_xxxyy_xxyyzz_0, ta_xxxyy_xxyyzz_1, ta_xxxyy_xxyzz_0, ta_xxxyy_xxyzz_1, ta_xxxyy_xxyzzz_0, ta_xxxyy_xxyzzz_1, ta_xxxyy_xyyyy_0, ta_xxxyy_xyyyy_1, ta_xxxyy_xyyyyy_0, ta_xxxyy_xyyyyy_1, ta_xxxyy_xyyyyz_0, ta_xxxyy_xyyyyz_1, ta_xxxyy_xyyyz_0, ta_xxxyy_xyyyz_1, ta_xxxyy_xyyyzz_0, ta_xxxyy_xyyyzz_1, ta_xxxyy_xyyzz_0, ta_xxxyy_xyyzz_1, ta_xxxyy_xyyzzz_0, ta_xxxyy_xyyzzz_1, ta_xxxyy_xyzzz_0, ta_xxxyy_xyzzz_1, ta_xxxyy_xyzzzz_0, ta_xxxyy_xyzzzz_1, ta_xxxyy_yyyyy_0, ta_xxxyy_yyyyy_1, ta_xxxyy_yyyyyy_0, ta_xxxyy_yyyyyy_1, ta_xxxyy_yyyyyz_0, ta_xxxyy_yyyyyz_1, ta_xxxyy_yyyyz_0, ta_xxxyy_yyyyz_1, ta_xxxyy_yyyyzz_0, ta_xxxyy_yyyyzz_1, ta_xxxyy_yyyzz_0, ta_xxxyy_yyyzz_1, ta_xxxyy_yyyzzz_0, ta_xxxyy_yyyzzz_1, ta_xxxyy_yyzzz_0, ta_xxxyy_yyzzz_1, ta_xxxyy_yyzzzz_0, ta_xxxyy_yyzzzz_1, ta_xxxyy_yzzzz_0, ta_xxxyy_yzzzz_1, ta_xxxyy_yzzzzz_0, ta_xxxyy_yzzzzz_1, ta_xxxyy_zzzzzz_0, ta_xxxyy_zzzzzz_1, ta_xxyy_xxxxxy_0, ta_xxyy_xxxxxy_1, ta_xxyy_xxxxyy_0, ta_xxyy_xxxxyy_1, ta_xxyy_xxxxyz_0, ta_xxyy_xxxxyz_1, ta_xxyy_xxxyyy_0, ta_xxyy_xxxyyy_1, ta_xxyy_xxxyyz_0, ta_xxyy_xxxyyz_1, ta_xxyy_xxxyzz_0, ta_xxyy_xxxyzz_1, ta_xxyy_xxyyyy_0, ta_xxyy_xxyyyy_1, ta_xxyy_xxyyyz_0, ta_xxyy_xxyyyz_1, ta_xxyy_xxyyzz_0, ta_xxyy_xxyyzz_1, ta_xxyy_xxyzzz_0, ta_xxyy_xxyzzz_1, ta_xxyy_xyyyyy_0, ta_xxyy_xyyyyy_1, ta_xxyy_xyyyyz_0, ta_xxyy_xyyyyz_1, ta_xxyy_xyyyzz_0, ta_xxyy_xyyyzz_1, ta_xxyy_xyyzzz_0, ta_xxyy_xyyzzz_1, ta_xxyy_xyzzzz_0, ta_xxyy_xyzzzz_1, ta_xxyy_yyyyyy_0, ta_xxyy_yyyyyy_1, ta_xxyy_yyyyyz_0, ta_xxyy_yyyyyz_1, ta_xxyy_yyyyzz_0, ta_xxyy_yyyyzz_1, ta_xxyy_yyyzzz_0, ta_xxyy_yyyzzz_1, ta_xxyy_yyzzzz_0, ta_xxyy_yyzzzz_1, ta_xxyy_yzzzzz_0, ta_xxyy_yzzzzz_1, ta_xxyy_zzzzzz_0, ta_xxyy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyy_xxxxxx_0[i] = ta_xxxx_xxxxxx_0[i] * fe_0 - ta_xxxx_xxxxxx_1[i] * fe_0 + ta_xxxxy_xxxxxx_0[i] * pa_y[i] - ta_xxxxy_xxxxxx_1[i] * pc_y[i];

        ta_xxxxyy_xxxxxy_0[i] = 3.0 * ta_xxyy_xxxxxy_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxxy_1[i] * fe_0 + 5.0 * ta_xxxyy_xxxxy_0[i] * fe_0 - 5.0 * ta_xxxyy_xxxxy_1[i] * fe_0 + ta_xxxyy_xxxxxy_0[i] * pa_x[i] - ta_xxxyy_xxxxxy_1[i] * pc_x[i];

        ta_xxxxyy_xxxxxz_0[i] = ta_xxxx_xxxxxz_0[i] * fe_0 - ta_xxxx_xxxxxz_1[i] * fe_0 + ta_xxxxy_xxxxxz_0[i] * pa_y[i] - ta_xxxxy_xxxxxz_1[i] * pc_y[i];

        ta_xxxxyy_xxxxyy_0[i] = 3.0 * ta_xxyy_xxxxyy_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxyy_1[i] * fe_0 + 4.0 * ta_xxxyy_xxxyy_0[i] * fe_0 - 4.0 * ta_xxxyy_xxxyy_1[i] * fe_0 + ta_xxxyy_xxxxyy_0[i] * pa_x[i] - ta_xxxyy_xxxxyy_1[i] * pc_x[i];

        ta_xxxxyy_xxxxyz_0[i] = 3.0 * ta_xxyy_xxxxyz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxxyy_xxxyz_0[i] * fe_0 - 4.0 * ta_xxxyy_xxxyz_1[i] * fe_0 + ta_xxxyy_xxxxyz_0[i] * pa_x[i] - ta_xxxyy_xxxxyz_1[i] * pc_x[i];

        ta_xxxxyy_xxxxzz_0[i] = ta_xxxx_xxxxzz_0[i] * fe_0 - ta_xxxx_xxxxzz_1[i] * fe_0 + ta_xxxxy_xxxxzz_0[i] * pa_y[i] - ta_xxxxy_xxxxzz_1[i] * pc_y[i];

        ta_xxxxyy_xxxyyy_0[i] = 3.0 * ta_xxyy_xxxyyy_0[i] * fe_0 - 3.0 * ta_xxyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_xxxyy_xxyyy_0[i] * fe_0 - 3.0 * ta_xxxyy_xxyyy_1[i] * fe_0 + ta_xxxyy_xxxyyy_0[i] * pa_x[i] - ta_xxxyy_xxxyyy_1[i] * pc_x[i];

        ta_xxxxyy_xxxyyz_0[i] = 3.0 * ta_xxyy_xxxyyz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxxyy_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxyy_xxyyz_1[i] * fe_0 + ta_xxxyy_xxxyyz_0[i] * pa_x[i] - ta_xxxyy_xxxyyz_1[i] * pc_x[i];

        ta_xxxxyy_xxxyzz_0[i] = 3.0 * ta_xxyy_xxxyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxxyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xxxyy_xxyzz_1[i] * fe_0 + ta_xxxyy_xxxyzz_0[i] * pa_x[i] - ta_xxxyy_xxxyzz_1[i] * pc_x[i];

        ta_xxxxyy_xxxzzz_0[i] = ta_xxxx_xxxzzz_0[i] * fe_0 - ta_xxxx_xxxzzz_1[i] * fe_0 + ta_xxxxy_xxxzzz_0[i] * pa_y[i] - ta_xxxxy_xxxzzz_1[i] * pc_y[i];

        ta_xxxxyy_xxyyyy_0[i] = 3.0 * ta_xxyy_xxyyyy_0[i] * fe_0 - 3.0 * ta_xxyy_xxyyyy_1[i] * fe_0 + 2.0 * ta_xxxyy_xyyyy_0[i] * fe_0 - 2.0 * ta_xxxyy_xyyyy_1[i] * fe_0 + ta_xxxyy_xxyyyy_0[i] * pa_x[i] - ta_xxxyy_xxyyyy_1[i] * pc_x[i];

        ta_xxxxyy_xxyyyz_0[i] = 3.0 * ta_xxyy_xxyyyz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxxyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xyyyz_1[i] * fe_0 + ta_xxxyy_xxyyyz_0[i] * pa_x[i] - ta_xxxyy_xxyyyz_1[i] * pc_x[i];

        ta_xxxxyy_xxyyzz_0[i] = 3.0 * ta_xxyy_xxyyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxxyy_xyyzz_0[i] * fe_0 - 2.0 * ta_xxxyy_xyyzz_1[i] * fe_0 + ta_xxxyy_xxyyzz_0[i] * pa_x[i] - ta_xxxyy_xxyyzz_1[i] * pc_x[i];

        ta_xxxxyy_xxyzzz_0[i] = 3.0 * ta_xxyy_xxyzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxxyy_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxyy_xyzzz_1[i] * fe_0 + ta_xxxyy_xxyzzz_0[i] * pa_x[i] - ta_xxxyy_xxyzzz_1[i] * pc_x[i];

        ta_xxxxyy_xxzzzz_0[i] = ta_xxxx_xxzzzz_0[i] * fe_0 - ta_xxxx_xxzzzz_1[i] * fe_0 + ta_xxxxy_xxzzzz_0[i] * pa_y[i] - ta_xxxxy_xxzzzz_1[i] * pc_y[i];

        ta_xxxxyy_xyyyyy_0[i] = 3.0 * ta_xxyy_xyyyyy_0[i] * fe_0 - 3.0 * ta_xxyy_xyyyyy_1[i] * fe_0 + ta_xxxyy_yyyyy_0[i] * fe_0 - ta_xxxyy_yyyyy_1[i] * fe_0 + ta_xxxyy_xyyyyy_0[i] * pa_x[i] - ta_xxxyy_xyyyyy_1[i] * pc_x[i];

        ta_xxxxyy_xyyyyz_0[i] = 3.0 * ta_xxyy_xyyyyz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyyyz_1[i] * fe_0 + ta_xxxyy_yyyyz_0[i] * fe_0 - ta_xxxyy_yyyyz_1[i] * fe_0 + ta_xxxyy_xyyyyz_0[i] * pa_x[i] - ta_xxxyy_xyyyyz_1[i] * pc_x[i];

        ta_xxxxyy_xyyyzz_0[i] = 3.0 * ta_xxyy_xyyyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyyzz_1[i] * fe_0 + ta_xxxyy_yyyzz_0[i] * fe_0 - ta_xxxyy_yyyzz_1[i] * fe_0 + ta_xxxyy_xyyyzz_0[i] * pa_x[i] - ta_xxxyy_xyyyzz_1[i] * pc_x[i];

        ta_xxxxyy_xyyzzz_0[i] = 3.0 * ta_xxyy_xyyzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyzzz_1[i] * fe_0 + ta_xxxyy_yyzzz_0[i] * fe_0 - ta_xxxyy_yyzzz_1[i] * fe_0 + ta_xxxyy_xyyzzz_0[i] * pa_x[i] - ta_xxxyy_xyyzzz_1[i] * pc_x[i];

        ta_xxxxyy_xyzzzz_0[i] = 3.0 * ta_xxyy_xyzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyzzzz_1[i] * fe_0 + ta_xxxyy_yzzzz_0[i] * fe_0 - ta_xxxyy_yzzzz_1[i] * fe_0 + ta_xxxyy_xyzzzz_0[i] * pa_x[i] - ta_xxxyy_xyzzzz_1[i] * pc_x[i];

        ta_xxxxyy_xzzzzz_0[i] = ta_xxxx_xzzzzz_0[i] * fe_0 - ta_xxxx_xzzzzz_1[i] * fe_0 + ta_xxxxy_xzzzzz_0[i] * pa_y[i] - ta_xxxxy_xzzzzz_1[i] * pc_y[i];

        ta_xxxxyy_yyyyyy_0[i] = 3.0 * ta_xxyy_yyyyyy_0[i] * fe_0 - 3.0 * ta_xxyy_yyyyyy_1[i] * fe_0 + ta_xxxyy_yyyyyy_0[i] * pa_x[i] - ta_xxxyy_yyyyyy_1[i] * pc_x[i];

        ta_xxxxyy_yyyyyz_0[i] = 3.0 * ta_xxyy_yyyyyz_0[i] * fe_0 - 3.0 * ta_xxyy_yyyyyz_1[i] * fe_0 + ta_xxxyy_yyyyyz_0[i] * pa_x[i] - ta_xxxyy_yyyyyz_1[i] * pc_x[i];

        ta_xxxxyy_yyyyzz_0[i] = 3.0 * ta_xxyy_yyyyzz_0[i] * fe_0 - 3.0 * ta_xxyy_yyyyzz_1[i] * fe_0 + ta_xxxyy_yyyyzz_0[i] * pa_x[i] - ta_xxxyy_yyyyzz_1[i] * pc_x[i];

        ta_xxxxyy_yyyzzz_0[i] = 3.0 * ta_xxyy_yyyzzz_0[i] * fe_0 - 3.0 * ta_xxyy_yyyzzz_1[i] * fe_0 + ta_xxxyy_yyyzzz_0[i] * pa_x[i] - ta_xxxyy_yyyzzz_1[i] * pc_x[i];

        ta_xxxxyy_yyzzzz_0[i] = 3.0 * ta_xxyy_yyzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_yyzzzz_1[i] * fe_0 + ta_xxxyy_yyzzzz_0[i] * pa_x[i] - ta_xxxyy_yyzzzz_1[i] * pc_x[i];

        ta_xxxxyy_yzzzzz_0[i] = 3.0 * ta_xxyy_yzzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_yzzzzz_1[i] * fe_0 + ta_xxxyy_yzzzzz_0[i] * pa_x[i] - ta_xxxyy_yzzzzz_1[i] * pc_x[i];

        ta_xxxxyy_zzzzzz_0[i] = 3.0 * ta_xxyy_zzzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_zzzzzz_1[i] * fe_0 + ta_xxxyy_zzzzzz_0[i] * pa_x[i] - ta_xxxyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 112-140 components of targeted buffer : II

    auto ta_xxxxyz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 112);

    auto ta_xxxxyz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 113);

    auto ta_xxxxyz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 114);

    auto ta_xxxxyz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 115);

    auto ta_xxxxyz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 116);

    auto ta_xxxxyz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 117);

    auto ta_xxxxyz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 118);

    auto ta_xxxxyz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 119);

    auto ta_xxxxyz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 120);

    auto ta_xxxxyz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 121);

    auto ta_xxxxyz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 122);

    auto ta_xxxxyz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 123);

    auto ta_xxxxyz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 124);

    auto ta_xxxxyz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 125);

    auto ta_xxxxyz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 126);

    auto ta_xxxxyz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 127);

    auto ta_xxxxyz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 128);

    auto ta_xxxxyz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 129);

    auto ta_xxxxyz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 130);

    auto ta_xxxxyz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 131);

    auto ta_xxxxyz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 132);

    auto ta_xxxxyz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 133);

    auto ta_xxxxyz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 134);

    auto ta_xxxxyz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 135);

    auto ta_xxxxyz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 136);

    auto ta_xxxxyz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 137);

    auto ta_xxxxyz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 138);

    auto ta_xxxxyz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 139);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxxxy_xxxxxy_0, ta_xxxxy_xxxxxy_1, ta_xxxxy_xxxxyy_0, ta_xxxxy_xxxxyy_1, ta_xxxxy_xxxyyy_0, ta_xxxxy_xxxyyy_1, ta_xxxxy_xxyyyy_0, ta_xxxxy_xxyyyy_1, ta_xxxxy_xyyyyy_0, ta_xxxxy_xyyyyy_1, ta_xxxxy_yyyyyy_0, ta_xxxxy_yyyyyy_1, ta_xxxxyz_xxxxxx_0, ta_xxxxyz_xxxxxy_0, ta_xxxxyz_xxxxxz_0, ta_xxxxyz_xxxxyy_0, ta_xxxxyz_xxxxyz_0, ta_xxxxyz_xxxxzz_0, ta_xxxxyz_xxxyyy_0, ta_xxxxyz_xxxyyz_0, ta_xxxxyz_xxxyzz_0, ta_xxxxyz_xxxzzz_0, ta_xxxxyz_xxyyyy_0, ta_xxxxyz_xxyyyz_0, ta_xxxxyz_xxyyzz_0, ta_xxxxyz_xxyzzz_0, ta_xxxxyz_xxzzzz_0, ta_xxxxyz_xyyyyy_0, ta_xxxxyz_xyyyyz_0, ta_xxxxyz_xyyyzz_0, ta_xxxxyz_xyyzzz_0, ta_xxxxyz_xyzzzz_0, ta_xxxxyz_xzzzzz_0, ta_xxxxyz_yyyyyy_0, ta_xxxxyz_yyyyyz_0, ta_xxxxyz_yyyyzz_0, ta_xxxxyz_yyyzzz_0, ta_xxxxyz_yyzzzz_0, ta_xxxxyz_yzzzzz_0, ta_xxxxyz_zzzzzz_0, ta_xxxxz_xxxxxx_0, ta_xxxxz_xxxxxx_1, ta_xxxxz_xxxxxz_0, ta_xxxxz_xxxxxz_1, ta_xxxxz_xxxxyz_0, ta_xxxxz_xxxxyz_1, ta_xxxxz_xxxxz_0, ta_xxxxz_xxxxz_1, ta_xxxxz_xxxxzz_0, ta_xxxxz_xxxxzz_1, ta_xxxxz_xxxyyz_0, ta_xxxxz_xxxyyz_1, ta_xxxxz_xxxyz_0, ta_xxxxz_xxxyz_1, ta_xxxxz_xxxyzz_0, ta_xxxxz_xxxyzz_1, ta_xxxxz_xxxzz_0, ta_xxxxz_xxxzz_1, ta_xxxxz_xxxzzz_0, ta_xxxxz_xxxzzz_1, ta_xxxxz_xxyyyz_0, ta_xxxxz_xxyyyz_1, ta_xxxxz_xxyyz_0, ta_xxxxz_xxyyz_1, ta_xxxxz_xxyyzz_0, ta_xxxxz_xxyyzz_1, ta_xxxxz_xxyzz_0, ta_xxxxz_xxyzz_1, ta_xxxxz_xxyzzz_0, ta_xxxxz_xxyzzz_1, ta_xxxxz_xxzzz_0, ta_xxxxz_xxzzz_1, ta_xxxxz_xxzzzz_0, ta_xxxxz_xxzzzz_1, ta_xxxxz_xyyyyz_0, ta_xxxxz_xyyyyz_1, ta_xxxxz_xyyyz_0, ta_xxxxz_xyyyz_1, ta_xxxxz_xyyyzz_0, ta_xxxxz_xyyyzz_1, ta_xxxxz_xyyzz_0, ta_xxxxz_xyyzz_1, ta_xxxxz_xyyzzz_0, ta_xxxxz_xyyzzz_1, ta_xxxxz_xyzzz_0, ta_xxxxz_xyzzz_1, ta_xxxxz_xyzzzz_0, ta_xxxxz_xyzzzz_1, ta_xxxxz_xzzzz_0, ta_xxxxz_xzzzz_1, ta_xxxxz_xzzzzz_0, ta_xxxxz_xzzzzz_1, ta_xxxxz_zzzzzz_0, ta_xxxxz_zzzzzz_1, ta_xxxyz_yyyyyz_0, ta_xxxyz_yyyyyz_1, ta_xxxyz_yyyyzz_0, ta_xxxyz_yyyyzz_1, ta_xxxyz_yyyzzz_0, ta_xxxyz_yyyzzz_1, ta_xxxyz_yyzzzz_0, ta_xxxyz_yyzzzz_1, ta_xxxyz_yzzzzz_0, ta_xxxyz_yzzzzz_1, ta_xxyz_yyyyyz_0, ta_xxyz_yyyyyz_1, ta_xxyz_yyyyzz_0, ta_xxyz_yyyyzz_1, ta_xxyz_yyyzzz_0, ta_xxyz_yyyzzz_1, ta_xxyz_yyzzzz_0, ta_xxyz_yyzzzz_1, ta_xxyz_yzzzzz_0, ta_xxyz_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyz_xxxxxx_0[i] = ta_xxxxz_xxxxxx_0[i] * pa_y[i] - ta_xxxxz_xxxxxx_1[i] * pc_y[i];

        ta_xxxxyz_xxxxxy_0[i] = ta_xxxxy_xxxxxy_0[i] * pa_z[i] - ta_xxxxy_xxxxxy_1[i] * pc_z[i];

        ta_xxxxyz_xxxxxz_0[i] = ta_xxxxz_xxxxxz_0[i] * pa_y[i] - ta_xxxxz_xxxxxz_1[i] * pc_y[i];

        ta_xxxxyz_xxxxyy_0[i] = ta_xxxxy_xxxxyy_0[i] * pa_z[i] - ta_xxxxy_xxxxyy_1[i] * pc_z[i];

        ta_xxxxyz_xxxxyz_0[i] = ta_xxxxz_xxxxz_0[i] * fe_0 - ta_xxxxz_xxxxz_1[i] * fe_0 + ta_xxxxz_xxxxyz_0[i] * pa_y[i] - ta_xxxxz_xxxxyz_1[i] * pc_y[i];

        ta_xxxxyz_xxxxzz_0[i] = ta_xxxxz_xxxxzz_0[i] * pa_y[i] - ta_xxxxz_xxxxzz_1[i] * pc_y[i];

        ta_xxxxyz_xxxyyy_0[i] = ta_xxxxy_xxxyyy_0[i] * pa_z[i] - ta_xxxxy_xxxyyy_1[i] * pc_z[i];

        ta_xxxxyz_xxxyyz_0[i] = 2.0 * ta_xxxxz_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxxz_xxxyz_1[i] * fe_0 + ta_xxxxz_xxxyyz_0[i] * pa_y[i] - ta_xxxxz_xxxyyz_1[i] * pc_y[i];

        ta_xxxxyz_xxxyzz_0[i] = ta_xxxxz_xxxzz_0[i] * fe_0 - ta_xxxxz_xxxzz_1[i] * fe_0 + ta_xxxxz_xxxyzz_0[i] * pa_y[i] - ta_xxxxz_xxxyzz_1[i] * pc_y[i];

        ta_xxxxyz_xxxzzz_0[i] = ta_xxxxz_xxxzzz_0[i] * pa_y[i] - ta_xxxxz_xxxzzz_1[i] * pc_y[i];

        ta_xxxxyz_xxyyyy_0[i] = ta_xxxxy_xxyyyy_0[i] * pa_z[i] - ta_xxxxy_xxyyyy_1[i] * pc_z[i];

        ta_xxxxyz_xxyyyz_0[i] = 3.0 * ta_xxxxz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxxz_xxyyz_1[i] * fe_0 + ta_xxxxz_xxyyyz_0[i] * pa_y[i] - ta_xxxxz_xxyyyz_1[i] * pc_y[i];

        ta_xxxxyz_xxyyzz_0[i] = 2.0 * ta_xxxxz_xxyzz_0[i] * fe_0 - 2.0 * ta_xxxxz_xxyzz_1[i] * fe_0 + ta_xxxxz_xxyyzz_0[i] * pa_y[i] - ta_xxxxz_xxyyzz_1[i] * pc_y[i];

        ta_xxxxyz_xxyzzz_0[i] = ta_xxxxz_xxzzz_0[i] * fe_0 - ta_xxxxz_xxzzz_1[i] * fe_0 + ta_xxxxz_xxyzzz_0[i] * pa_y[i] - ta_xxxxz_xxyzzz_1[i] * pc_y[i];

        ta_xxxxyz_xxzzzz_0[i] = ta_xxxxz_xxzzzz_0[i] * pa_y[i] - ta_xxxxz_xxzzzz_1[i] * pc_y[i];

        ta_xxxxyz_xyyyyy_0[i] = ta_xxxxy_xyyyyy_0[i] * pa_z[i] - ta_xxxxy_xyyyyy_1[i] * pc_z[i];

        ta_xxxxyz_xyyyyz_0[i] = 4.0 * ta_xxxxz_xyyyz_0[i] * fe_0 - 4.0 * ta_xxxxz_xyyyz_1[i] * fe_0 + ta_xxxxz_xyyyyz_0[i] * pa_y[i] - ta_xxxxz_xyyyyz_1[i] * pc_y[i];

        ta_xxxxyz_xyyyzz_0[i] = 3.0 * ta_xxxxz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxxz_xyyzz_1[i] * fe_0 + ta_xxxxz_xyyyzz_0[i] * pa_y[i] - ta_xxxxz_xyyyzz_1[i] * pc_y[i];

        ta_xxxxyz_xyyzzz_0[i] = 2.0 * ta_xxxxz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxxz_xyzzz_1[i] * fe_0 + ta_xxxxz_xyyzzz_0[i] * pa_y[i] - ta_xxxxz_xyyzzz_1[i] * pc_y[i];

        ta_xxxxyz_xyzzzz_0[i] = ta_xxxxz_xzzzz_0[i] * fe_0 - ta_xxxxz_xzzzz_1[i] * fe_0 + ta_xxxxz_xyzzzz_0[i] * pa_y[i] - ta_xxxxz_xyzzzz_1[i] * pc_y[i];

        ta_xxxxyz_xzzzzz_0[i] = ta_xxxxz_xzzzzz_0[i] * pa_y[i] - ta_xxxxz_xzzzzz_1[i] * pc_y[i];

        ta_xxxxyz_yyyyyy_0[i] = ta_xxxxy_yyyyyy_0[i] * pa_z[i] - ta_xxxxy_yyyyyy_1[i] * pc_z[i];

        ta_xxxxyz_yyyyyz_0[i] = 3.0 * ta_xxyz_yyyyyz_0[i] * fe_0 - 3.0 * ta_xxyz_yyyyyz_1[i] * fe_0 + ta_xxxyz_yyyyyz_0[i] * pa_x[i] - ta_xxxyz_yyyyyz_1[i] * pc_x[i];

        ta_xxxxyz_yyyyzz_0[i] = 3.0 * ta_xxyz_yyyyzz_0[i] * fe_0 - 3.0 * ta_xxyz_yyyyzz_1[i] * fe_0 + ta_xxxyz_yyyyzz_0[i] * pa_x[i] - ta_xxxyz_yyyyzz_1[i] * pc_x[i];

        ta_xxxxyz_yyyzzz_0[i] = 3.0 * ta_xxyz_yyyzzz_0[i] * fe_0 - 3.0 * ta_xxyz_yyyzzz_1[i] * fe_0 + ta_xxxyz_yyyzzz_0[i] * pa_x[i] - ta_xxxyz_yyyzzz_1[i] * pc_x[i];

        ta_xxxxyz_yyzzzz_0[i] = 3.0 * ta_xxyz_yyzzzz_0[i] * fe_0 - 3.0 * ta_xxyz_yyzzzz_1[i] * fe_0 + ta_xxxyz_yyzzzz_0[i] * pa_x[i] - ta_xxxyz_yyzzzz_1[i] * pc_x[i];

        ta_xxxxyz_yzzzzz_0[i] = 3.0 * ta_xxyz_yzzzzz_0[i] * fe_0 - 3.0 * ta_xxyz_yzzzzz_1[i] * fe_0 + ta_xxxyz_yzzzzz_0[i] * pa_x[i] - ta_xxxyz_yzzzzz_1[i] * pc_x[i];

        ta_xxxxyz_zzzzzz_0[i] = ta_xxxxz_zzzzzz_0[i] * pa_y[i] - ta_xxxxz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : II

    auto ta_xxxxzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 140);

    auto ta_xxxxzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 141);

    auto ta_xxxxzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 142);

    auto ta_xxxxzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 143);

    auto ta_xxxxzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 144);

    auto ta_xxxxzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 145);

    auto ta_xxxxzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 146);

    auto ta_xxxxzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 147);

    auto ta_xxxxzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 148);

    auto ta_xxxxzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 149);

    auto ta_xxxxzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 150);

    auto ta_xxxxzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 151);

    auto ta_xxxxzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 152);

    auto ta_xxxxzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 153);

    auto ta_xxxxzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 154);

    auto ta_xxxxzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 155);

    auto ta_xxxxzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 156);

    auto ta_xxxxzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 157);

    auto ta_xxxxzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 158);

    auto ta_xxxxzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 159);

    auto ta_xxxxzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 160);

    auto ta_xxxxzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 161);

    auto ta_xxxxzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 162);

    auto ta_xxxxzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 163);

    auto ta_xxxxzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 164);

    auto ta_xxxxzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 165);

    auto ta_xxxxzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 166);

    auto ta_xxxxzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 167);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xxxx_xxxxxx_0, ta_xxxx_xxxxxx_1, ta_xxxx_xxxxxy_0, ta_xxxx_xxxxxy_1, ta_xxxx_xxxxyy_0, ta_xxxx_xxxxyy_1, ta_xxxx_xxxyyy_0, ta_xxxx_xxxyyy_1, ta_xxxx_xxyyyy_0, ta_xxxx_xxyyyy_1, ta_xxxx_xyyyyy_0, ta_xxxx_xyyyyy_1, ta_xxxxz_xxxxxx_0, ta_xxxxz_xxxxxx_1, ta_xxxxz_xxxxxy_0, ta_xxxxz_xxxxxy_1, ta_xxxxz_xxxxyy_0, ta_xxxxz_xxxxyy_1, ta_xxxxz_xxxyyy_0, ta_xxxxz_xxxyyy_1, ta_xxxxz_xxyyyy_0, ta_xxxxz_xxyyyy_1, ta_xxxxz_xyyyyy_0, ta_xxxxz_xyyyyy_1, ta_xxxxzz_xxxxxx_0, ta_xxxxzz_xxxxxy_0, ta_xxxxzz_xxxxxz_0, ta_xxxxzz_xxxxyy_0, ta_xxxxzz_xxxxyz_0, ta_xxxxzz_xxxxzz_0, ta_xxxxzz_xxxyyy_0, ta_xxxxzz_xxxyyz_0, ta_xxxxzz_xxxyzz_0, ta_xxxxzz_xxxzzz_0, ta_xxxxzz_xxyyyy_0, ta_xxxxzz_xxyyyz_0, ta_xxxxzz_xxyyzz_0, ta_xxxxzz_xxyzzz_0, ta_xxxxzz_xxzzzz_0, ta_xxxxzz_xyyyyy_0, ta_xxxxzz_xyyyyz_0, ta_xxxxzz_xyyyzz_0, ta_xxxxzz_xyyzzz_0, ta_xxxxzz_xyzzzz_0, ta_xxxxzz_xzzzzz_0, ta_xxxxzz_yyyyyy_0, ta_xxxxzz_yyyyyz_0, ta_xxxxzz_yyyyzz_0, ta_xxxxzz_yyyzzz_0, ta_xxxxzz_yyzzzz_0, ta_xxxxzz_yzzzzz_0, ta_xxxxzz_zzzzzz_0, ta_xxxzz_xxxxxz_0, ta_xxxzz_xxxxxz_1, ta_xxxzz_xxxxyz_0, ta_xxxzz_xxxxyz_1, ta_xxxzz_xxxxz_0, ta_xxxzz_xxxxz_1, ta_xxxzz_xxxxzz_0, ta_xxxzz_xxxxzz_1, ta_xxxzz_xxxyyz_0, ta_xxxzz_xxxyyz_1, ta_xxxzz_xxxyz_0, ta_xxxzz_xxxyz_1, ta_xxxzz_xxxyzz_0, ta_xxxzz_xxxyzz_1, ta_xxxzz_xxxzz_0, ta_xxxzz_xxxzz_1, ta_xxxzz_xxxzzz_0, ta_xxxzz_xxxzzz_1, ta_xxxzz_xxyyyz_0, ta_xxxzz_xxyyyz_1, ta_xxxzz_xxyyz_0, ta_xxxzz_xxyyz_1, ta_xxxzz_xxyyzz_0, ta_xxxzz_xxyyzz_1, ta_xxxzz_xxyzz_0, ta_xxxzz_xxyzz_1, ta_xxxzz_xxyzzz_0, ta_xxxzz_xxyzzz_1, ta_xxxzz_xxzzz_0, ta_xxxzz_xxzzz_1, ta_xxxzz_xxzzzz_0, ta_xxxzz_xxzzzz_1, ta_xxxzz_xyyyyz_0, ta_xxxzz_xyyyyz_1, ta_xxxzz_xyyyz_0, ta_xxxzz_xyyyz_1, ta_xxxzz_xyyyzz_0, ta_xxxzz_xyyyzz_1, ta_xxxzz_xyyzz_0, ta_xxxzz_xyyzz_1, ta_xxxzz_xyyzzz_0, ta_xxxzz_xyyzzz_1, ta_xxxzz_xyzzz_0, ta_xxxzz_xyzzz_1, ta_xxxzz_xyzzzz_0, ta_xxxzz_xyzzzz_1, ta_xxxzz_xzzzz_0, ta_xxxzz_xzzzz_1, ta_xxxzz_xzzzzz_0, ta_xxxzz_xzzzzz_1, ta_xxxzz_yyyyyy_0, ta_xxxzz_yyyyyy_1, ta_xxxzz_yyyyyz_0, ta_xxxzz_yyyyyz_1, ta_xxxzz_yyyyz_0, ta_xxxzz_yyyyz_1, ta_xxxzz_yyyyzz_0, ta_xxxzz_yyyyzz_1, ta_xxxzz_yyyzz_0, ta_xxxzz_yyyzz_1, ta_xxxzz_yyyzzz_0, ta_xxxzz_yyyzzz_1, ta_xxxzz_yyzzz_0, ta_xxxzz_yyzzz_1, ta_xxxzz_yyzzzz_0, ta_xxxzz_yyzzzz_1, ta_xxxzz_yzzzz_0, ta_xxxzz_yzzzz_1, ta_xxxzz_yzzzzz_0, ta_xxxzz_yzzzzz_1, ta_xxxzz_zzzzz_0, ta_xxxzz_zzzzz_1, ta_xxxzz_zzzzzz_0, ta_xxxzz_zzzzzz_1, ta_xxzz_xxxxxz_0, ta_xxzz_xxxxxz_1, ta_xxzz_xxxxyz_0, ta_xxzz_xxxxyz_1, ta_xxzz_xxxxzz_0, ta_xxzz_xxxxzz_1, ta_xxzz_xxxyyz_0, ta_xxzz_xxxyyz_1, ta_xxzz_xxxyzz_0, ta_xxzz_xxxyzz_1, ta_xxzz_xxxzzz_0, ta_xxzz_xxxzzz_1, ta_xxzz_xxyyyz_0, ta_xxzz_xxyyyz_1, ta_xxzz_xxyyzz_0, ta_xxzz_xxyyzz_1, ta_xxzz_xxyzzz_0, ta_xxzz_xxyzzz_1, ta_xxzz_xxzzzz_0, ta_xxzz_xxzzzz_1, ta_xxzz_xyyyyz_0, ta_xxzz_xyyyyz_1, ta_xxzz_xyyyzz_0, ta_xxzz_xyyyzz_1, ta_xxzz_xyyzzz_0, ta_xxzz_xyyzzz_1, ta_xxzz_xyzzzz_0, ta_xxzz_xyzzzz_1, ta_xxzz_xzzzzz_0, ta_xxzz_xzzzzz_1, ta_xxzz_yyyyyy_0, ta_xxzz_yyyyyy_1, ta_xxzz_yyyyyz_0, ta_xxzz_yyyyyz_1, ta_xxzz_yyyyzz_0, ta_xxzz_yyyyzz_1, ta_xxzz_yyyzzz_0, ta_xxzz_yyyzzz_1, ta_xxzz_yyzzzz_0, ta_xxzz_yyzzzz_1, ta_xxzz_yzzzzz_0, ta_xxzz_yzzzzz_1, ta_xxzz_zzzzzz_0, ta_xxzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxzz_xxxxxx_0[i] = ta_xxxx_xxxxxx_0[i] * fe_0 - ta_xxxx_xxxxxx_1[i] * fe_0 + ta_xxxxz_xxxxxx_0[i] * pa_z[i] - ta_xxxxz_xxxxxx_1[i] * pc_z[i];

        ta_xxxxzz_xxxxxy_0[i] = ta_xxxx_xxxxxy_0[i] * fe_0 - ta_xxxx_xxxxxy_1[i] * fe_0 + ta_xxxxz_xxxxxy_0[i] * pa_z[i] - ta_xxxxz_xxxxxy_1[i] * pc_z[i];

        ta_xxxxzz_xxxxxz_0[i] = 3.0 * ta_xxzz_xxxxxz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxxz_1[i] * fe_0 + 5.0 * ta_xxxzz_xxxxz_0[i] * fe_0 - 5.0 * ta_xxxzz_xxxxz_1[i] * fe_0 + ta_xxxzz_xxxxxz_0[i] * pa_x[i] - ta_xxxzz_xxxxxz_1[i] * pc_x[i];

        ta_xxxxzz_xxxxyy_0[i] = ta_xxxx_xxxxyy_0[i] * fe_0 - ta_xxxx_xxxxyy_1[i] * fe_0 + ta_xxxxz_xxxxyy_0[i] * pa_z[i] - ta_xxxxz_xxxxyy_1[i] * pc_z[i];

        ta_xxxxzz_xxxxyz_0[i] = 3.0 * ta_xxzz_xxxxyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxxzz_xxxyz_0[i] * fe_0 - 4.0 * ta_xxxzz_xxxyz_1[i] * fe_0 + ta_xxxzz_xxxxyz_0[i] * pa_x[i] - ta_xxxzz_xxxxyz_1[i] * pc_x[i];

        ta_xxxxzz_xxxxzz_0[i] = 3.0 * ta_xxzz_xxxxzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxzz_1[i] * fe_0 + 4.0 * ta_xxxzz_xxxzz_0[i] * fe_0 - 4.0 * ta_xxxzz_xxxzz_1[i] * fe_0 + ta_xxxzz_xxxxzz_0[i] * pa_x[i] - ta_xxxzz_xxxxzz_1[i] * pc_x[i];

        ta_xxxxzz_xxxyyy_0[i] = ta_xxxx_xxxyyy_0[i] * fe_0 - ta_xxxx_xxxyyy_1[i] * fe_0 + ta_xxxxz_xxxyyy_0[i] * pa_z[i] - ta_xxxxz_xxxyyy_1[i] * pc_z[i];

        ta_xxxxzz_xxxyyz_0[i] = 3.0 * ta_xxzz_xxxyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxxzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxzz_xxyyz_1[i] * fe_0 + ta_xxxzz_xxxyyz_0[i] * pa_x[i] - ta_xxxzz_xxxyyz_1[i] * pc_x[i];

        ta_xxxxzz_xxxyzz_0[i] = 3.0 * ta_xxzz_xxxyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxxzz_xxyzz_0[i] * fe_0 - 3.0 * ta_xxxzz_xxyzz_1[i] * fe_0 + ta_xxxzz_xxxyzz_0[i] * pa_x[i] - ta_xxxzz_xxxyzz_1[i] * pc_x[i];

        ta_xxxxzz_xxxzzz_0[i] = 3.0 * ta_xxzz_xxxzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_xxxzz_xxzzz_0[i] * fe_0 - 3.0 * ta_xxxzz_xxzzz_1[i] * fe_0 + ta_xxxzz_xxxzzz_0[i] * pa_x[i] - ta_xxxzz_xxxzzz_1[i] * pc_x[i];

        ta_xxxxzz_xxyyyy_0[i] = ta_xxxx_xxyyyy_0[i] * fe_0 - ta_xxxx_xxyyyy_1[i] * fe_0 + ta_xxxxz_xxyyyy_0[i] * pa_z[i] - ta_xxxxz_xxyyyy_1[i] * pc_z[i];

        ta_xxxxzz_xxyyyz_0[i] = 3.0 * ta_xxzz_xxyyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxxzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xxxzz_xyyyz_1[i] * fe_0 + ta_xxxzz_xxyyyz_0[i] * pa_x[i] - ta_xxxzz_xxyyyz_1[i] * pc_x[i];

        ta_xxxxzz_xxyyzz_0[i] = 3.0 * ta_xxzz_xxyyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxxzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xxxzz_xyyzz_1[i] * fe_0 + ta_xxxzz_xxyyzz_0[i] * pa_x[i] - ta_xxxzz_xxyyzz_1[i] * pc_x[i];

        ta_xxxxzz_xxyzzz_0[i] = 3.0 * ta_xxzz_xxyzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxxzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxzz_xyzzz_1[i] * fe_0 + ta_xxxzz_xxyzzz_0[i] * pa_x[i] - ta_xxxzz_xxyzzz_1[i] * pc_x[i];

        ta_xxxxzz_xxzzzz_0[i] = 3.0 * ta_xxzz_xxzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxzzzz_1[i] * fe_0 + 2.0 * ta_xxxzz_xzzzz_0[i] * fe_0 - 2.0 * ta_xxxzz_xzzzz_1[i] * fe_0 + ta_xxxzz_xxzzzz_0[i] * pa_x[i] - ta_xxxzz_xxzzzz_1[i] * pc_x[i];

        ta_xxxxzz_xyyyyy_0[i] = ta_xxxx_xyyyyy_0[i] * fe_0 - ta_xxxx_xyyyyy_1[i] * fe_0 + ta_xxxxz_xyyyyy_0[i] * pa_z[i] - ta_xxxxz_xyyyyy_1[i] * pc_z[i];

        ta_xxxxzz_xyyyyz_0[i] = 3.0 * ta_xxzz_xyyyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyyyz_1[i] * fe_0 + ta_xxxzz_yyyyz_0[i] * fe_0 - ta_xxxzz_yyyyz_1[i] * fe_0 + ta_xxxzz_xyyyyz_0[i] * pa_x[i] - ta_xxxzz_xyyyyz_1[i] * pc_x[i];

        ta_xxxxzz_xyyyzz_0[i] = 3.0 * ta_xxzz_xyyyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyyzz_1[i] * fe_0 + ta_xxxzz_yyyzz_0[i] * fe_0 - ta_xxxzz_yyyzz_1[i] * fe_0 + ta_xxxzz_xyyyzz_0[i] * pa_x[i] - ta_xxxzz_xyyyzz_1[i] * pc_x[i];

        ta_xxxxzz_xyyzzz_0[i] = 3.0 * ta_xxzz_xyyzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyzzz_1[i] * fe_0 + ta_xxxzz_yyzzz_0[i] * fe_0 - ta_xxxzz_yyzzz_1[i] * fe_0 + ta_xxxzz_xyyzzz_0[i] * pa_x[i] - ta_xxxzz_xyyzzz_1[i] * pc_x[i];

        ta_xxxxzz_xyzzzz_0[i] = 3.0 * ta_xxzz_xyzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyzzzz_1[i] * fe_0 + ta_xxxzz_yzzzz_0[i] * fe_0 - ta_xxxzz_yzzzz_1[i] * fe_0 + ta_xxxzz_xyzzzz_0[i] * pa_x[i] - ta_xxxzz_xyzzzz_1[i] * pc_x[i];

        ta_xxxxzz_xzzzzz_0[i] = 3.0 * ta_xxzz_xzzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xzzzzz_1[i] * fe_0 + ta_xxxzz_zzzzz_0[i] * fe_0 - ta_xxxzz_zzzzz_1[i] * fe_0 + ta_xxxzz_xzzzzz_0[i] * pa_x[i] - ta_xxxzz_xzzzzz_1[i] * pc_x[i];

        ta_xxxxzz_yyyyyy_0[i] = 3.0 * ta_xxzz_yyyyyy_0[i] * fe_0 - 3.0 * ta_xxzz_yyyyyy_1[i] * fe_0 + ta_xxxzz_yyyyyy_0[i] * pa_x[i] - ta_xxxzz_yyyyyy_1[i] * pc_x[i];

        ta_xxxxzz_yyyyyz_0[i] = 3.0 * ta_xxzz_yyyyyz_0[i] * fe_0 - 3.0 * ta_xxzz_yyyyyz_1[i] * fe_0 + ta_xxxzz_yyyyyz_0[i] * pa_x[i] - ta_xxxzz_yyyyyz_1[i] * pc_x[i];

        ta_xxxxzz_yyyyzz_0[i] = 3.0 * ta_xxzz_yyyyzz_0[i] * fe_0 - 3.0 * ta_xxzz_yyyyzz_1[i] * fe_0 + ta_xxxzz_yyyyzz_0[i] * pa_x[i] - ta_xxxzz_yyyyzz_1[i] * pc_x[i];

        ta_xxxxzz_yyyzzz_0[i] = 3.0 * ta_xxzz_yyyzzz_0[i] * fe_0 - 3.0 * ta_xxzz_yyyzzz_1[i] * fe_0 + ta_xxxzz_yyyzzz_0[i] * pa_x[i] - ta_xxxzz_yyyzzz_1[i] * pc_x[i];

        ta_xxxxzz_yyzzzz_0[i] = 3.0 * ta_xxzz_yyzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_yyzzzz_1[i] * fe_0 + ta_xxxzz_yyzzzz_0[i] * pa_x[i] - ta_xxxzz_yyzzzz_1[i] * pc_x[i];

        ta_xxxxzz_yzzzzz_0[i] = 3.0 * ta_xxzz_yzzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_yzzzzz_1[i] * fe_0 + ta_xxxzz_yzzzzz_0[i] * pa_x[i] - ta_xxxzz_yzzzzz_1[i] * pc_x[i];

        ta_xxxxzz_zzzzzz_0[i] = 3.0 * ta_xxzz_zzzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_zzzzzz_1[i] * fe_0 + ta_xxxzz_zzzzzz_0[i] * pa_x[i] - ta_xxxzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 168-196 components of targeted buffer : II

    auto ta_xxxyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 168);

    auto ta_xxxyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 169);

    auto ta_xxxyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 170);

    auto ta_xxxyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 171);

    auto ta_xxxyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 172);

    auto ta_xxxyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 173);

    auto ta_xxxyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 174);

    auto ta_xxxyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 175);

    auto ta_xxxyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 176);

    auto ta_xxxyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 177);

    auto ta_xxxyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 178);

    auto ta_xxxyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 179);

    auto ta_xxxyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 180);

    auto ta_xxxyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 181);

    auto ta_xxxyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 182);

    auto ta_xxxyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 183);

    auto ta_xxxyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 184);

    auto ta_xxxyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 185);

    auto ta_xxxyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 186);

    auto ta_xxxyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 187);

    auto ta_xxxyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 188);

    auto ta_xxxyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 189);

    auto ta_xxxyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 190);

    auto ta_xxxyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 191);

    auto ta_xxxyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 192);

    auto ta_xxxyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 193);

    auto ta_xxxyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 194);

    auto ta_xxxyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 195);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxxy_xxxxxx_0, ta_xxxy_xxxxxx_1, ta_xxxy_xxxxxz_0, ta_xxxy_xxxxxz_1, ta_xxxy_xxxxzz_0, ta_xxxy_xxxxzz_1, ta_xxxy_xxxzzz_0, ta_xxxy_xxxzzz_1, ta_xxxy_xxzzzz_0, ta_xxxy_xxzzzz_1, ta_xxxy_xzzzzz_0, ta_xxxy_xzzzzz_1, ta_xxxyy_xxxxxx_0, ta_xxxyy_xxxxxx_1, ta_xxxyy_xxxxxz_0, ta_xxxyy_xxxxxz_1, ta_xxxyy_xxxxzz_0, ta_xxxyy_xxxxzz_1, ta_xxxyy_xxxzzz_0, ta_xxxyy_xxxzzz_1, ta_xxxyy_xxzzzz_0, ta_xxxyy_xxzzzz_1, ta_xxxyy_xzzzzz_0, ta_xxxyy_xzzzzz_1, ta_xxxyyy_xxxxxx_0, ta_xxxyyy_xxxxxy_0, ta_xxxyyy_xxxxxz_0, ta_xxxyyy_xxxxyy_0, ta_xxxyyy_xxxxyz_0, ta_xxxyyy_xxxxzz_0, ta_xxxyyy_xxxyyy_0, ta_xxxyyy_xxxyyz_0, ta_xxxyyy_xxxyzz_0, ta_xxxyyy_xxxzzz_0, ta_xxxyyy_xxyyyy_0, ta_xxxyyy_xxyyyz_0, ta_xxxyyy_xxyyzz_0, ta_xxxyyy_xxyzzz_0, ta_xxxyyy_xxzzzz_0, ta_xxxyyy_xyyyyy_0, ta_xxxyyy_xyyyyz_0, ta_xxxyyy_xyyyzz_0, ta_xxxyyy_xyyzzz_0, ta_xxxyyy_xyzzzz_0, ta_xxxyyy_xzzzzz_0, ta_xxxyyy_yyyyyy_0, ta_xxxyyy_yyyyyz_0, ta_xxxyyy_yyyyzz_0, ta_xxxyyy_yyyzzz_0, ta_xxxyyy_yyzzzz_0, ta_xxxyyy_yzzzzz_0, ta_xxxyyy_zzzzzz_0, ta_xxyyy_xxxxxy_0, ta_xxyyy_xxxxxy_1, ta_xxyyy_xxxxy_0, ta_xxyyy_xxxxy_1, ta_xxyyy_xxxxyy_0, ta_xxyyy_xxxxyy_1, ta_xxyyy_xxxxyz_0, ta_xxyyy_xxxxyz_1, ta_xxyyy_xxxyy_0, ta_xxyyy_xxxyy_1, ta_xxyyy_xxxyyy_0, ta_xxyyy_xxxyyy_1, ta_xxyyy_xxxyyz_0, ta_xxyyy_xxxyyz_1, ta_xxyyy_xxxyz_0, ta_xxyyy_xxxyz_1, ta_xxyyy_xxxyzz_0, ta_xxyyy_xxxyzz_1, ta_xxyyy_xxyyy_0, ta_xxyyy_xxyyy_1, ta_xxyyy_xxyyyy_0, ta_xxyyy_xxyyyy_1, ta_xxyyy_xxyyyz_0, ta_xxyyy_xxyyyz_1, ta_xxyyy_xxyyz_0, ta_xxyyy_xxyyz_1, ta_xxyyy_xxyyzz_0, ta_xxyyy_xxyyzz_1, ta_xxyyy_xxyzz_0, ta_xxyyy_xxyzz_1, ta_xxyyy_xxyzzz_0, ta_xxyyy_xxyzzz_1, ta_xxyyy_xyyyy_0, ta_xxyyy_xyyyy_1, ta_xxyyy_xyyyyy_0, ta_xxyyy_xyyyyy_1, ta_xxyyy_xyyyyz_0, ta_xxyyy_xyyyyz_1, ta_xxyyy_xyyyz_0, ta_xxyyy_xyyyz_1, ta_xxyyy_xyyyzz_0, ta_xxyyy_xyyyzz_1, ta_xxyyy_xyyzz_0, ta_xxyyy_xyyzz_1, ta_xxyyy_xyyzzz_0, ta_xxyyy_xyyzzz_1, ta_xxyyy_xyzzz_0, ta_xxyyy_xyzzz_1, ta_xxyyy_xyzzzz_0, ta_xxyyy_xyzzzz_1, ta_xxyyy_yyyyy_0, ta_xxyyy_yyyyy_1, ta_xxyyy_yyyyyy_0, ta_xxyyy_yyyyyy_1, ta_xxyyy_yyyyyz_0, ta_xxyyy_yyyyyz_1, ta_xxyyy_yyyyz_0, ta_xxyyy_yyyyz_1, ta_xxyyy_yyyyzz_0, ta_xxyyy_yyyyzz_1, ta_xxyyy_yyyzz_0, ta_xxyyy_yyyzz_1, ta_xxyyy_yyyzzz_0, ta_xxyyy_yyyzzz_1, ta_xxyyy_yyzzz_0, ta_xxyyy_yyzzz_1, ta_xxyyy_yyzzzz_0, ta_xxyyy_yyzzzz_1, ta_xxyyy_yzzzz_0, ta_xxyyy_yzzzz_1, ta_xxyyy_yzzzzz_0, ta_xxyyy_yzzzzz_1, ta_xxyyy_zzzzzz_0, ta_xxyyy_zzzzzz_1, ta_xyyy_xxxxxy_0, ta_xyyy_xxxxxy_1, ta_xyyy_xxxxyy_0, ta_xyyy_xxxxyy_1, ta_xyyy_xxxxyz_0, ta_xyyy_xxxxyz_1, ta_xyyy_xxxyyy_0, ta_xyyy_xxxyyy_1, ta_xyyy_xxxyyz_0, ta_xyyy_xxxyyz_1, ta_xyyy_xxxyzz_0, ta_xyyy_xxxyzz_1, ta_xyyy_xxyyyy_0, ta_xyyy_xxyyyy_1, ta_xyyy_xxyyyz_0, ta_xyyy_xxyyyz_1, ta_xyyy_xxyyzz_0, ta_xyyy_xxyyzz_1, ta_xyyy_xxyzzz_0, ta_xyyy_xxyzzz_1, ta_xyyy_xyyyyy_0, ta_xyyy_xyyyyy_1, ta_xyyy_xyyyyz_0, ta_xyyy_xyyyyz_1, ta_xyyy_xyyyzz_0, ta_xyyy_xyyyzz_1, ta_xyyy_xyyzzz_0, ta_xyyy_xyyzzz_1, ta_xyyy_xyzzzz_0, ta_xyyy_xyzzzz_1, ta_xyyy_yyyyyy_0, ta_xyyy_yyyyyy_1, ta_xyyy_yyyyyz_0, ta_xyyy_yyyyyz_1, ta_xyyy_yyyyzz_0, ta_xyyy_yyyyzz_1, ta_xyyy_yyyzzz_0, ta_xyyy_yyyzzz_1, ta_xyyy_yyzzzz_0, ta_xyyy_yyzzzz_1, ta_xyyy_yzzzzz_0, ta_xyyy_yzzzzz_1, ta_xyyy_zzzzzz_0, ta_xyyy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyy_xxxxxx_0[i] = 2.0 * ta_xxxy_xxxxxx_0[i] * fe_0 - 2.0 * ta_xxxy_xxxxxx_1[i] * fe_0 + ta_xxxyy_xxxxxx_0[i] * pa_y[i] - ta_xxxyy_xxxxxx_1[i] * pc_y[i];

        ta_xxxyyy_xxxxxy_0[i] = 2.0 * ta_xyyy_xxxxxy_0[i] * fe_0 - 2.0 * ta_xyyy_xxxxxy_1[i] * fe_0 + 5.0 * ta_xxyyy_xxxxy_0[i] * fe_0 - 5.0 * ta_xxyyy_xxxxy_1[i] * fe_0 + ta_xxyyy_xxxxxy_0[i] * pa_x[i] - ta_xxyyy_xxxxxy_1[i] * pc_x[i];

        ta_xxxyyy_xxxxxz_0[i] = 2.0 * ta_xxxy_xxxxxz_0[i] * fe_0 - 2.0 * ta_xxxy_xxxxxz_1[i] * fe_0 + ta_xxxyy_xxxxxz_0[i] * pa_y[i] - ta_xxxyy_xxxxxz_1[i] * pc_y[i];

        ta_xxxyyy_xxxxyy_0[i] = 2.0 * ta_xyyy_xxxxyy_0[i] * fe_0 - 2.0 * ta_xyyy_xxxxyy_1[i] * fe_0 + 4.0 * ta_xxyyy_xxxyy_0[i] * fe_0 - 4.0 * ta_xxyyy_xxxyy_1[i] * fe_0 + ta_xxyyy_xxxxyy_0[i] * pa_x[i] - ta_xxyyy_xxxxyy_1[i] * pc_x[i];

        ta_xxxyyy_xxxxyz_0[i] = 2.0 * ta_xyyy_xxxxyz_0[i] * fe_0 - 2.0 * ta_xyyy_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxyyy_xxxyz_0[i] * fe_0 - 4.0 * ta_xxyyy_xxxyz_1[i] * fe_0 + ta_xxyyy_xxxxyz_0[i] * pa_x[i] - ta_xxyyy_xxxxyz_1[i] * pc_x[i];

        ta_xxxyyy_xxxxzz_0[i] = 2.0 * ta_xxxy_xxxxzz_0[i] * fe_0 - 2.0 * ta_xxxy_xxxxzz_1[i] * fe_0 + ta_xxxyy_xxxxzz_0[i] * pa_y[i] - ta_xxxyy_xxxxzz_1[i] * pc_y[i];

        ta_xxxyyy_xxxyyy_0[i] = 2.0 * ta_xyyy_xxxyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_xxyyy_xxyyy_0[i] * fe_0 - 3.0 * ta_xxyyy_xxyyy_1[i] * fe_0 + ta_xxyyy_xxxyyy_0[i] * pa_x[i] - ta_xxyyy_xxxyyy_1[i] * pc_x[i];

        ta_xxxyyy_xxxyyz_0[i] = 2.0 * ta_xyyy_xxxyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxyyy_xxyyz_0[i] * fe_0 - 3.0 * ta_xxyyy_xxyyz_1[i] * fe_0 + ta_xxyyy_xxxyyz_0[i] * pa_x[i] - ta_xxyyy_xxxyyz_1[i] * pc_x[i];

        ta_xxxyyy_xxxyzz_0[i] = 2.0 * ta_xyyy_xxxyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xxyyy_xxyzz_1[i] * fe_0 + ta_xxyyy_xxxyzz_0[i] * pa_x[i] - ta_xxyyy_xxxyzz_1[i] * pc_x[i];

        ta_xxxyyy_xxxzzz_0[i] = 2.0 * ta_xxxy_xxxzzz_0[i] * fe_0 - 2.0 * ta_xxxy_xxxzzz_1[i] * fe_0 + ta_xxxyy_xxxzzz_0[i] * pa_y[i] - ta_xxxyy_xxxzzz_1[i] * pc_y[i];

        ta_xxxyyy_xxyyyy_0[i] = 2.0 * ta_xyyy_xxyyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xxyyyy_1[i] * fe_0 + 2.0 * ta_xxyyy_xyyyy_0[i] * fe_0 - 2.0 * ta_xxyyy_xyyyy_1[i] * fe_0 + ta_xxyyy_xxyyyy_0[i] * pa_x[i] - ta_xxyyy_xxyyyy_1[i] * pc_x[i];

        ta_xxxyyy_xxyyyz_0[i] = 2.0 * ta_xyyy_xxyyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xyyyz_1[i] * fe_0 + ta_xxyyy_xxyyyz_0[i] * pa_x[i] - ta_xxyyy_xxyyyz_1[i] * pc_x[i];

        ta_xxxyyy_xxyyzz_0[i] = 2.0 * ta_xyyy_xxyyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxyyy_xyyzz_0[i] * fe_0 - 2.0 * ta_xxyyy_xyyzz_1[i] * fe_0 + ta_xxyyy_xxyyzz_0[i] * pa_x[i] - ta_xxyyy_xxyyzz_1[i] * pc_x[i];

        ta_xxxyyy_xxyzzz_0[i] = 2.0 * ta_xyyy_xxyzzz_0[i] * fe_0 - 2.0 * ta_xyyy_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxyyy_xyzzz_0[i] * fe_0 - 2.0 * ta_xxyyy_xyzzz_1[i] * fe_0 + ta_xxyyy_xxyzzz_0[i] * pa_x[i] - ta_xxyyy_xxyzzz_1[i] * pc_x[i];

        ta_xxxyyy_xxzzzz_0[i] = 2.0 * ta_xxxy_xxzzzz_0[i] * fe_0 - 2.0 * ta_xxxy_xxzzzz_1[i] * fe_0 + ta_xxxyy_xxzzzz_0[i] * pa_y[i] - ta_xxxyy_xxzzzz_1[i] * pc_y[i];

        ta_xxxyyy_xyyyyy_0[i] = 2.0 * ta_xyyy_xyyyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xyyyyy_1[i] * fe_0 + ta_xxyyy_yyyyy_0[i] * fe_0 - ta_xxyyy_yyyyy_1[i] * fe_0 + ta_xxyyy_xyyyyy_0[i] * pa_x[i] - ta_xxyyy_xyyyyy_1[i] * pc_x[i];

        ta_xxxyyy_xyyyyz_0[i] = 2.0 * ta_xyyy_xyyyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyyyz_1[i] * fe_0 + ta_xxyyy_yyyyz_0[i] * fe_0 - ta_xxyyy_yyyyz_1[i] * fe_0 + ta_xxyyy_xyyyyz_0[i] * pa_x[i] - ta_xxyyy_xyyyyz_1[i] * pc_x[i];

        ta_xxxyyy_xyyyzz_0[i] = 2.0 * ta_xyyy_xyyyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyyzz_1[i] * fe_0 + ta_xxyyy_yyyzz_0[i] * fe_0 - ta_xxyyy_yyyzz_1[i] * fe_0 + ta_xxyyy_xyyyzz_0[i] * pa_x[i] - ta_xxyyy_xyyyzz_1[i] * pc_x[i];

        ta_xxxyyy_xyyzzz_0[i] = 2.0 * ta_xyyy_xyyzzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyzzz_1[i] * fe_0 + ta_xxyyy_yyzzz_0[i] * fe_0 - ta_xxyyy_yyzzz_1[i] * fe_0 + ta_xxyyy_xyyzzz_0[i] * pa_x[i] - ta_xxyyy_xyyzzz_1[i] * pc_x[i];

        ta_xxxyyy_xyzzzz_0[i] = 2.0 * ta_xyyy_xyzzzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyzzzz_1[i] * fe_0 + ta_xxyyy_yzzzz_0[i] * fe_0 - ta_xxyyy_yzzzz_1[i] * fe_0 + ta_xxyyy_xyzzzz_0[i] * pa_x[i] - ta_xxyyy_xyzzzz_1[i] * pc_x[i];

        ta_xxxyyy_xzzzzz_0[i] = 2.0 * ta_xxxy_xzzzzz_0[i] * fe_0 - 2.0 * ta_xxxy_xzzzzz_1[i] * fe_0 + ta_xxxyy_xzzzzz_0[i] * pa_y[i] - ta_xxxyy_xzzzzz_1[i] * pc_y[i];

        ta_xxxyyy_yyyyyy_0[i] = 2.0 * ta_xyyy_yyyyyy_0[i] * fe_0 - 2.0 * ta_xyyy_yyyyyy_1[i] * fe_0 + ta_xxyyy_yyyyyy_0[i] * pa_x[i] - ta_xxyyy_yyyyyy_1[i] * pc_x[i];

        ta_xxxyyy_yyyyyz_0[i] = 2.0 * ta_xyyy_yyyyyz_0[i] * fe_0 - 2.0 * ta_xyyy_yyyyyz_1[i] * fe_0 + ta_xxyyy_yyyyyz_0[i] * pa_x[i] - ta_xxyyy_yyyyyz_1[i] * pc_x[i];

        ta_xxxyyy_yyyyzz_0[i] = 2.0 * ta_xyyy_yyyyzz_0[i] * fe_0 - 2.0 * ta_xyyy_yyyyzz_1[i] * fe_0 + ta_xxyyy_yyyyzz_0[i] * pa_x[i] - ta_xxyyy_yyyyzz_1[i] * pc_x[i];

        ta_xxxyyy_yyyzzz_0[i] = 2.0 * ta_xyyy_yyyzzz_0[i] * fe_0 - 2.0 * ta_xyyy_yyyzzz_1[i] * fe_0 + ta_xxyyy_yyyzzz_0[i] * pa_x[i] - ta_xxyyy_yyyzzz_1[i] * pc_x[i];

        ta_xxxyyy_yyzzzz_0[i] = 2.0 * ta_xyyy_yyzzzz_0[i] * fe_0 - 2.0 * ta_xyyy_yyzzzz_1[i] * fe_0 + ta_xxyyy_yyzzzz_0[i] * pa_x[i] - ta_xxyyy_yyzzzz_1[i] * pc_x[i];

        ta_xxxyyy_yzzzzz_0[i] = 2.0 * ta_xyyy_yzzzzz_0[i] * fe_0 - 2.0 * ta_xyyy_yzzzzz_1[i] * fe_0 + ta_xxyyy_yzzzzz_0[i] * pa_x[i] - ta_xxyyy_yzzzzz_1[i] * pc_x[i];

        ta_xxxyyy_zzzzzz_0[i] = 2.0 * ta_xyyy_zzzzzz_0[i] * fe_0 - 2.0 * ta_xyyy_zzzzzz_1[i] * fe_0 + ta_xxyyy_zzzzzz_0[i] * pa_x[i] - ta_xxyyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 196-224 components of targeted buffer : II

    auto ta_xxxyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 196);

    auto ta_xxxyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 197);

    auto ta_xxxyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 198);

    auto ta_xxxyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 199);

    auto ta_xxxyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 200);

    auto ta_xxxyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 201);

    auto ta_xxxyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 202);

    auto ta_xxxyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 203);

    auto ta_xxxyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 204);

    auto ta_xxxyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 205);

    auto ta_xxxyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 206);

    auto ta_xxxyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 207);

    auto ta_xxxyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 208);

    auto ta_xxxyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 209);

    auto ta_xxxyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 210);

    auto ta_xxxyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 211);

    auto ta_xxxyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 212);

    auto ta_xxxyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 213);

    auto ta_xxxyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 214);

    auto ta_xxxyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 215);

    auto ta_xxxyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 216);

    auto ta_xxxyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 217);

    auto ta_xxxyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 218);

    auto ta_xxxyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 219);

    auto ta_xxxyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 220);

    auto ta_xxxyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 221);

    auto ta_xxxyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 222);

    auto ta_xxxyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 223);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxxyy_xxxxxx_0, ta_xxxyy_xxxxxx_1, ta_xxxyy_xxxxxy_0, ta_xxxyy_xxxxxy_1, ta_xxxyy_xxxxy_0, ta_xxxyy_xxxxy_1, ta_xxxyy_xxxxyy_0, ta_xxxyy_xxxxyy_1, ta_xxxyy_xxxxyz_0, ta_xxxyy_xxxxyz_1, ta_xxxyy_xxxyy_0, ta_xxxyy_xxxyy_1, ta_xxxyy_xxxyyy_0, ta_xxxyy_xxxyyy_1, ta_xxxyy_xxxyyz_0, ta_xxxyy_xxxyyz_1, ta_xxxyy_xxxyz_0, ta_xxxyy_xxxyz_1, ta_xxxyy_xxxyzz_0, ta_xxxyy_xxxyzz_1, ta_xxxyy_xxyyy_0, ta_xxxyy_xxyyy_1, ta_xxxyy_xxyyyy_0, ta_xxxyy_xxyyyy_1, ta_xxxyy_xxyyyz_0, ta_xxxyy_xxyyyz_1, ta_xxxyy_xxyyz_0, ta_xxxyy_xxyyz_1, ta_xxxyy_xxyyzz_0, ta_xxxyy_xxyyzz_1, ta_xxxyy_xxyzz_0, ta_xxxyy_xxyzz_1, ta_xxxyy_xxyzzz_0, ta_xxxyy_xxyzzz_1, ta_xxxyy_xyyyy_0, ta_xxxyy_xyyyy_1, ta_xxxyy_xyyyyy_0, ta_xxxyy_xyyyyy_1, ta_xxxyy_xyyyyz_0, ta_xxxyy_xyyyyz_1, ta_xxxyy_xyyyz_0, ta_xxxyy_xyyyz_1, ta_xxxyy_xyyyzz_0, ta_xxxyy_xyyyzz_1, ta_xxxyy_xyyzz_0, ta_xxxyy_xyyzz_1, ta_xxxyy_xyyzzz_0, ta_xxxyy_xyyzzz_1, ta_xxxyy_xyzzz_0, ta_xxxyy_xyzzz_1, ta_xxxyy_xyzzzz_0, ta_xxxyy_xyzzzz_1, ta_xxxyy_yyyyyy_0, ta_xxxyy_yyyyyy_1, ta_xxxyyz_xxxxxx_0, ta_xxxyyz_xxxxxy_0, ta_xxxyyz_xxxxxz_0, ta_xxxyyz_xxxxyy_0, ta_xxxyyz_xxxxyz_0, ta_xxxyyz_xxxxzz_0, ta_xxxyyz_xxxyyy_0, ta_xxxyyz_xxxyyz_0, ta_xxxyyz_xxxyzz_0, ta_xxxyyz_xxxzzz_0, ta_xxxyyz_xxyyyy_0, ta_xxxyyz_xxyyyz_0, ta_xxxyyz_xxyyzz_0, ta_xxxyyz_xxyzzz_0, ta_xxxyyz_xxzzzz_0, ta_xxxyyz_xyyyyy_0, ta_xxxyyz_xyyyyz_0, ta_xxxyyz_xyyyzz_0, ta_xxxyyz_xyyzzz_0, ta_xxxyyz_xyzzzz_0, ta_xxxyyz_xzzzzz_0, ta_xxxyyz_yyyyyy_0, ta_xxxyyz_yyyyyz_0, ta_xxxyyz_yyyyzz_0, ta_xxxyyz_yyyzzz_0, ta_xxxyyz_yyzzzz_0, ta_xxxyyz_yzzzzz_0, ta_xxxyyz_zzzzzz_0, ta_xxxyz_xxxxxz_0, ta_xxxyz_xxxxxz_1, ta_xxxyz_xxxxzz_0, ta_xxxyz_xxxxzz_1, ta_xxxyz_xxxzzz_0, ta_xxxyz_xxxzzz_1, ta_xxxyz_xxzzzz_0, ta_xxxyz_xxzzzz_1, ta_xxxyz_xzzzzz_0, ta_xxxyz_xzzzzz_1, ta_xxxz_xxxxxz_0, ta_xxxz_xxxxxz_1, ta_xxxz_xxxxzz_0, ta_xxxz_xxxxzz_1, ta_xxxz_xxxzzz_0, ta_xxxz_xxxzzz_1, ta_xxxz_xxzzzz_0, ta_xxxz_xxzzzz_1, ta_xxxz_xzzzzz_0, ta_xxxz_xzzzzz_1, ta_xxyyz_yyyyyz_0, ta_xxyyz_yyyyyz_1, ta_xxyyz_yyyyzz_0, ta_xxyyz_yyyyzz_1, ta_xxyyz_yyyzzz_0, ta_xxyyz_yyyzzz_1, ta_xxyyz_yyzzzz_0, ta_xxyyz_yyzzzz_1, ta_xxyyz_yzzzzz_0, ta_xxyyz_yzzzzz_1, ta_xxyyz_zzzzzz_0, ta_xxyyz_zzzzzz_1, ta_xyyz_yyyyyz_0, ta_xyyz_yyyyyz_1, ta_xyyz_yyyyzz_0, ta_xyyz_yyyyzz_1, ta_xyyz_yyyzzz_0, ta_xyyz_yyyzzz_1, ta_xyyz_yyzzzz_0, ta_xyyz_yyzzzz_1, ta_xyyz_yzzzzz_0, ta_xyyz_yzzzzz_1, ta_xyyz_zzzzzz_0, ta_xyyz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyz_xxxxxx_0[i] = ta_xxxyy_xxxxxx_0[i] * pa_z[i] - ta_xxxyy_xxxxxx_1[i] * pc_z[i];

        ta_xxxyyz_xxxxxy_0[i] = ta_xxxyy_xxxxxy_0[i] * pa_z[i] - ta_xxxyy_xxxxxy_1[i] * pc_z[i];

        ta_xxxyyz_xxxxxz_0[i] = ta_xxxz_xxxxxz_0[i] * fe_0 - ta_xxxz_xxxxxz_1[i] * fe_0 + ta_xxxyz_xxxxxz_0[i] * pa_y[i] - ta_xxxyz_xxxxxz_1[i] * pc_y[i];

        ta_xxxyyz_xxxxyy_0[i] = ta_xxxyy_xxxxyy_0[i] * pa_z[i] - ta_xxxyy_xxxxyy_1[i] * pc_z[i];

        ta_xxxyyz_xxxxyz_0[i] = ta_xxxyy_xxxxy_0[i] * fe_0 - ta_xxxyy_xxxxy_1[i] * fe_0 + ta_xxxyy_xxxxyz_0[i] * pa_z[i] - ta_xxxyy_xxxxyz_1[i] * pc_z[i];

        ta_xxxyyz_xxxxzz_0[i] = ta_xxxz_xxxxzz_0[i] * fe_0 - ta_xxxz_xxxxzz_1[i] * fe_0 + ta_xxxyz_xxxxzz_0[i] * pa_y[i] - ta_xxxyz_xxxxzz_1[i] * pc_y[i];

        ta_xxxyyz_xxxyyy_0[i] = ta_xxxyy_xxxyyy_0[i] * pa_z[i] - ta_xxxyy_xxxyyy_1[i] * pc_z[i];

        ta_xxxyyz_xxxyyz_0[i] = ta_xxxyy_xxxyy_0[i] * fe_0 - ta_xxxyy_xxxyy_1[i] * fe_0 + ta_xxxyy_xxxyyz_0[i] * pa_z[i] - ta_xxxyy_xxxyyz_1[i] * pc_z[i];

        ta_xxxyyz_xxxyzz_0[i] = 2.0 * ta_xxxyy_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xxxyz_1[i] * fe_0 + ta_xxxyy_xxxyzz_0[i] * pa_z[i] - ta_xxxyy_xxxyzz_1[i] * pc_z[i];

        ta_xxxyyz_xxxzzz_0[i] = ta_xxxz_xxxzzz_0[i] * fe_0 - ta_xxxz_xxxzzz_1[i] * fe_0 + ta_xxxyz_xxxzzz_0[i] * pa_y[i] - ta_xxxyz_xxxzzz_1[i] * pc_y[i];

        ta_xxxyyz_xxyyyy_0[i] = ta_xxxyy_xxyyyy_0[i] * pa_z[i] - ta_xxxyy_xxyyyy_1[i] * pc_z[i];

        ta_xxxyyz_xxyyyz_0[i] = ta_xxxyy_xxyyy_0[i] * fe_0 - ta_xxxyy_xxyyy_1[i] * fe_0 + ta_xxxyy_xxyyyz_0[i] * pa_z[i] - ta_xxxyy_xxyyyz_1[i] * pc_z[i];

        ta_xxxyyz_xxyyzz_0[i] = 2.0 * ta_xxxyy_xxyyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xxyyz_1[i] * fe_0 + ta_xxxyy_xxyyzz_0[i] * pa_z[i] - ta_xxxyy_xxyyzz_1[i] * pc_z[i];

        ta_xxxyyz_xxyzzz_0[i] = 3.0 * ta_xxxyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xxxyy_xxyzz_1[i] * fe_0 + ta_xxxyy_xxyzzz_0[i] * pa_z[i] - ta_xxxyy_xxyzzz_1[i] * pc_z[i];

        ta_xxxyyz_xxzzzz_0[i] = ta_xxxz_xxzzzz_0[i] * fe_0 - ta_xxxz_xxzzzz_1[i] * fe_0 + ta_xxxyz_xxzzzz_0[i] * pa_y[i] - ta_xxxyz_xxzzzz_1[i] * pc_y[i];

        ta_xxxyyz_xyyyyy_0[i] = ta_xxxyy_xyyyyy_0[i] * pa_z[i] - ta_xxxyy_xyyyyy_1[i] * pc_z[i];

        ta_xxxyyz_xyyyyz_0[i] = ta_xxxyy_xyyyy_0[i] * fe_0 - ta_xxxyy_xyyyy_1[i] * fe_0 + ta_xxxyy_xyyyyz_0[i] * pa_z[i] - ta_xxxyy_xyyyyz_1[i] * pc_z[i];

        ta_xxxyyz_xyyyzz_0[i] = 2.0 * ta_xxxyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xyyyz_1[i] * fe_0 + ta_xxxyy_xyyyzz_0[i] * pa_z[i] - ta_xxxyy_xyyyzz_1[i] * pc_z[i];

        ta_xxxyyz_xyyzzz_0[i] = 3.0 * ta_xxxyy_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxyy_xyyzz_1[i] * fe_0 + ta_xxxyy_xyyzzz_0[i] * pa_z[i] - ta_xxxyy_xyyzzz_1[i] * pc_z[i];

        ta_xxxyyz_xyzzzz_0[i] = 4.0 * ta_xxxyy_xyzzz_0[i] * fe_0 - 4.0 * ta_xxxyy_xyzzz_1[i] * fe_0 + ta_xxxyy_xyzzzz_0[i] * pa_z[i] - ta_xxxyy_xyzzzz_1[i] * pc_z[i];

        ta_xxxyyz_xzzzzz_0[i] = ta_xxxz_xzzzzz_0[i] * fe_0 - ta_xxxz_xzzzzz_1[i] * fe_0 + ta_xxxyz_xzzzzz_0[i] * pa_y[i] - ta_xxxyz_xzzzzz_1[i] * pc_y[i];

        ta_xxxyyz_yyyyyy_0[i] = ta_xxxyy_yyyyyy_0[i] * pa_z[i] - ta_xxxyy_yyyyyy_1[i] * pc_z[i];

        ta_xxxyyz_yyyyyz_0[i] = 2.0 * ta_xyyz_yyyyyz_0[i] * fe_0 - 2.0 * ta_xyyz_yyyyyz_1[i] * fe_0 + ta_xxyyz_yyyyyz_0[i] * pa_x[i] - ta_xxyyz_yyyyyz_1[i] * pc_x[i];

        ta_xxxyyz_yyyyzz_0[i] = 2.0 * ta_xyyz_yyyyzz_0[i] * fe_0 - 2.0 * ta_xyyz_yyyyzz_1[i] * fe_0 + ta_xxyyz_yyyyzz_0[i] * pa_x[i] - ta_xxyyz_yyyyzz_1[i] * pc_x[i];

        ta_xxxyyz_yyyzzz_0[i] = 2.0 * ta_xyyz_yyyzzz_0[i] * fe_0 - 2.0 * ta_xyyz_yyyzzz_1[i] * fe_0 + ta_xxyyz_yyyzzz_0[i] * pa_x[i] - ta_xxyyz_yyyzzz_1[i] * pc_x[i];

        ta_xxxyyz_yyzzzz_0[i] = 2.0 * ta_xyyz_yyzzzz_0[i] * fe_0 - 2.0 * ta_xyyz_yyzzzz_1[i] * fe_0 + ta_xxyyz_yyzzzz_0[i] * pa_x[i] - ta_xxyyz_yyzzzz_1[i] * pc_x[i];

        ta_xxxyyz_yzzzzz_0[i] = 2.0 * ta_xyyz_yzzzzz_0[i] * fe_0 - 2.0 * ta_xyyz_yzzzzz_1[i] * fe_0 + ta_xxyyz_yzzzzz_0[i] * pa_x[i] - ta_xxyyz_yzzzzz_1[i] * pc_x[i];

        ta_xxxyyz_zzzzzz_0[i] = 2.0 * ta_xyyz_zzzzzz_0[i] * fe_0 - 2.0 * ta_xyyz_zzzzzz_1[i] * fe_0 + ta_xxyyz_zzzzzz_0[i] * pa_x[i] - ta_xxyyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 224-252 components of targeted buffer : II

    auto ta_xxxyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 224);

    auto ta_xxxyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 225);

    auto ta_xxxyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 226);

    auto ta_xxxyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 227);

    auto ta_xxxyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 228);

    auto ta_xxxyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 229);

    auto ta_xxxyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 230);

    auto ta_xxxyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 231);

    auto ta_xxxyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 232);

    auto ta_xxxyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 233);

    auto ta_xxxyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 234);

    auto ta_xxxyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 235);

    auto ta_xxxyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 236);

    auto ta_xxxyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 237);

    auto ta_xxxyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 238);

    auto ta_xxxyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 239);

    auto ta_xxxyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 240);

    auto ta_xxxyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 241);

    auto ta_xxxyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 242);

    auto ta_xxxyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 243);

    auto ta_xxxyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 244);

    auto ta_xxxyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 245);

    auto ta_xxxyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 246);

    auto ta_xxxyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 247);

    auto ta_xxxyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 248);

    auto ta_xxxyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 249);

    auto ta_xxxyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 250);

    auto ta_xxxyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 251);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxxyzz_xxxxxx_0, ta_xxxyzz_xxxxxy_0, ta_xxxyzz_xxxxxz_0, ta_xxxyzz_xxxxyy_0, ta_xxxyzz_xxxxyz_0, ta_xxxyzz_xxxxzz_0, ta_xxxyzz_xxxyyy_0, ta_xxxyzz_xxxyyz_0, ta_xxxyzz_xxxyzz_0, ta_xxxyzz_xxxzzz_0, ta_xxxyzz_xxyyyy_0, ta_xxxyzz_xxyyyz_0, ta_xxxyzz_xxyyzz_0, ta_xxxyzz_xxyzzz_0, ta_xxxyzz_xxzzzz_0, ta_xxxyzz_xyyyyy_0, ta_xxxyzz_xyyyyz_0, ta_xxxyzz_xyyyzz_0, ta_xxxyzz_xyyzzz_0, ta_xxxyzz_xyzzzz_0, ta_xxxyzz_xzzzzz_0, ta_xxxyzz_yyyyyy_0, ta_xxxyzz_yyyyyz_0, ta_xxxyzz_yyyyzz_0, ta_xxxyzz_yyyzzz_0, ta_xxxyzz_yyzzzz_0, ta_xxxyzz_yzzzzz_0, ta_xxxyzz_zzzzzz_0, ta_xxxzz_xxxxx_0, ta_xxxzz_xxxxx_1, ta_xxxzz_xxxxxx_0, ta_xxxzz_xxxxxx_1, ta_xxxzz_xxxxxy_0, ta_xxxzz_xxxxxy_1, ta_xxxzz_xxxxxz_0, ta_xxxzz_xxxxxz_1, ta_xxxzz_xxxxy_0, ta_xxxzz_xxxxy_1, ta_xxxzz_xxxxyy_0, ta_xxxzz_xxxxyy_1, ta_xxxzz_xxxxyz_0, ta_xxxzz_xxxxyz_1, ta_xxxzz_xxxxz_0, ta_xxxzz_xxxxz_1, ta_xxxzz_xxxxzz_0, ta_xxxzz_xxxxzz_1, ta_xxxzz_xxxyy_0, ta_xxxzz_xxxyy_1, ta_xxxzz_xxxyyy_0, ta_xxxzz_xxxyyy_1, ta_xxxzz_xxxyyz_0, ta_xxxzz_xxxyyz_1, ta_xxxzz_xxxyz_0, ta_xxxzz_xxxyz_1, ta_xxxzz_xxxyzz_0, ta_xxxzz_xxxyzz_1, ta_xxxzz_xxxzz_0, ta_xxxzz_xxxzz_1, ta_xxxzz_xxxzzz_0, ta_xxxzz_xxxzzz_1, ta_xxxzz_xxyyy_0, ta_xxxzz_xxyyy_1, ta_xxxzz_xxyyyy_0, ta_xxxzz_xxyyyy_1, ta_xxxzz_xxyyyz_0, ta_xxxzz_xxyyyz_1, ta_xxxzz_xxyyz_0, ta_xxxzz_xxyyz_1, ta_xxxzz_xxyyzz_0, ta_xxxzz_xxyyzz_1, ta_xxxzz_xxyzz_0, ta_xxxzz_xxyzz_1, ta_xxxzz_xxyzzz_0, ta_xxxzz_xxyzzz_1, ta_xxxzz_xxzzz_0, ta_xxxzz_xxzzz_1, ta_xxxzz_xxzzzz_0, ta_xxxzz_xxzzzz_1, ta_xxxzz_xyyyy_0, ta_xxxzz_xyyyy_1, ta_xxxzz_xyyyyy_0, ta_xxxzz_xyyyyy_1, ta_xxxzz_xyyyyz_0, ta_xxxzz_xyyyyz_1, ta_xxxzz_xyyyz_0, ta_xxxzz_xyyyz_1, ta_xxxzz_xyyyzz_0, ta_xxxzz_xyyyzz_1, ta_xxxzz_xyyzz_0, ta_xxxzz_xyyzz_1, ta_xxxzz_xyyzzz_0, ta_xxxzz_xyyzzz_1, ta_xxxzz_xyzzz_0, ta_xxxzz_xyzzz_1, ta_xxxzz_xyzzzz_0, ta_xxxzz_xyzzzz_1, ta_xxxzz_xzzzz_0, ta_xxxzz_xzzzz_1, ta_xxxzz_xzzzzz_0, ta_xxxzz_xzzzzz_1, ta_xxxzz_zzzzzz_0, ta_xxxzz_zzzzzz_1, ta_xxyzz_yyyyyy_0, ta_xxyzz_yyyyyy_1, ta_xxyzz_yyyyyz_0, ta_xxyzz_yyyyyz_1, ta_xxyzz_yyyyzz_0, ta_xxyzz_yyyyzz_1, ta_xxyzz_yyyzzz_0, ta_xxyzz_yyyzzz_1, ta_xxyzz_yyzzzz_0, ta_xxyzz_yyzzzz_1, ta_xxyzz_yzzzzz_0, ta_xxyzz_yzzzzz_1, ta_xyzz_yyyyyy_0, ta_xyzz_yyyyyy_1, ta_xyzz_yyyyyz_0, ta_xyzz_yyyyyz_1, ta_xyzz_yyyyzz_0, ta_xyzz_yyyyzz_1, ta_xyzz_yyyzzz_0, ta_xyzz_yyyzzz_1, ta_xyzz_yyzzzz_0, ta_xyzz_yyzzzz_1, ta_xyzz_yzzzzz_0, ta_xyzz_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyzz_xxxxxx_0[i] = ta_xxxzz_xxxxxx_0[i] * pa_y[i] - ta_xxxzz_xxxxxx_1[i] * pc_y[i];

        ta_xxxyzz_xxxxxy_0[i] = ta_xxxzz_xxxxx_0[i] * fe_0 - ta_xxxzz_xxxxx_1[i] * fe_0 + ta_xxxzz_xxxxxy_0[i] * pa_y[i] - ta_xxxzz_xxxxxy_1[i] * pc_y[i];

        ta_xxxyzz_xxxxxz_0[i] = ta_xxxzz_xxxxxz_0[i] * pa_y[i] - ta_xxxzz_xxxxxz_1[i] * pc_y[i];

        ta_xxxyzz_xxxxyy_0[i] = 2.0 * ta_xxxzz_xxxxy_0[i] * fe_0 - 2.0 * ta_xxxzz_xxxxy_1[i] * fe_0 + ta_xxxzz_xxxxyy_0[i] * pa_y[i] - ta_xxxzz_xxxxyy_1[i] * pc_y[i];

        ta_xxxyzz_xxxxyz_0[i] = ta_xxxzz_xxxxz_0[i] * fe_0 - ta_xxxzz_xxxxz_1[i] * fe_0 + ta_xxxzz_xxxxyz_0[i] * pa_y[i] - ta_xxxzz_xxxxyz_1[i] * pc_y[i];

        ta_xxxyzz_xxxxzz_0[i] = ta_xxxzz_xxxxzz_0[i] * pa_y[i] - ta_xxxzz_xxxxzz_1[i] * pc_y[i];

        ta_xxxyzz_xxxyyy_0[i] = 3.0 * ta_xxxzz_xxxyy_0[i] * fe_0 - 3.0 * ta_xxxzz_xxxyy_1[i] * fe_0 + ta_xxxzz_xxxyyy_0[i] * pa_y[i] - ta_xxxzz_xxxyyy_1[i] * pc_y[i];

        ta_xxxyzz_xxxyyz_0[i] = 2.0 * ta_xxxzz_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxzz_xxxyz_1[i] * fe_0 + ta_xxxzz_xxxyyz_0[i] * pa_y[i] - ta_xxxzz_xxxyyz_1[i] * pc_y[i];

        ta_xxxyzz_xxxyzz_0[i] = ta_xxxzz_xxxzz_0[i] * fe_0 - ta_xxxzz_xxxzz_1[i] * fe_0 + ta_xxxzz_xxxyzz_0[i] * pa_y[i] - ta_xxxzz_xxxyzz_1[i] * pc_y[i];

        ta_xxxyzz_xxxzzz_0[i] = ta_xxxzz_xxxzzz_0[i] * pa_y[i] - ta_xxxzz_xxxzzz_1[i] * pc_y[i];

        ta_xxxyzz_xxyyyy_0[i] = 4.0 * ta_xxxzz_xxyyy_0[i] * fe_0 - 4.0 * ta_xxxzz_xxyyy_1[i] * fe_0 + ta_xxxzz_xxyyyy_0[i] * pa_y[i] - ta_xxxzz_xxyyyy_1[i] * pc_y[i];

        ta_xxxyzz_xxyyyz_0[i] = 3.0 * ta_xxxzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxzz_xxyyz_1[i] * fe_0 + ta_xxxzz_xxyyyz_0[i] * pa_y[i] - ta_xxxzz_xxyyyz_1[i] * pc_y[i];

        ta_xxxyzz_xxyyzz_0[i] = 2.0 * ta_xxxzz_xxyzz_0[i] * fe_0 - 2.0 * ta_xxxzz_xxyzz_1[i] * fe_0 + ta_xxxzz_xxyyzz_0[i] * pa_y[i] - ta_xxxzz_xxyyzz_1[i] * pc_y[i];

        ta_xxxyzz_xxyzzz_0[i] = ta_xxxzz_xxzzz_0[i] * fe_0 - ta_xxxzz_xxzzz_1[i] * fe_0 + ta_xxxzz_xxyzzz_0[i] * pa_y[i] - ta_xxxzz_xxyzzz_1[i] * pc_y[i];

        ta_xxxyzz_xxzzzz_0[i] = ta_xxxzz_xxzzzz_0[i] * pa_y[i] - ta_xxxzz_xxzzzz_1[i] * pc_y[i];

        ta_xxxyzz_xyyyyy_0[i] = 5.0 * ta_xxxzz_xyyyy_0[i] * fe_0 - 5.0 * ta_xxxzz_xyyyy_1[i] * fe_0 + ta_xxxzz_xyyyyy_0[i] * pa_y[i] - ta_xxxzz_xyyyyy_1[i] * pc_y[i];

        ta_xxxyzz_xyyyyz_0[i] = 4.0 * ta_xxxzz_xyyyz_0[i] * fe_0 - 4.0 * ta_xxxzz_xyyyz_1[i] * fe_0 + ta_xxxzz_xyyyyz_0[i] * pa_y[i] - ta_xxxzz_xyyyyz_1[i] * pc_y[i];

        ta_xxxyzz_xyyyzz_0[i] = 3.0 * ta_xxxzz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxzz_xyyzz_1[i] * fe_0 + ta_xxxzz_xyyyzz_0[i] * pa_y[i] - ta_xxxzz_xyyyzz_1[i] * pc_y[i];

        ta_xxxyzz_xyyzzz_0[i] = 2.0 * ta_xxxzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxzz_xyzzz_1[i] * fe_0 + ta_xxxzz_xyyzzz_0[i] * pa_y[i] - ta_xxxzz_xyyzzz_1[i] * pc_y[i];

        ta_xxxyzz_xyzzzz_0[i] = ta_xxxzz_xzzzz_0[i] * fe_0 - ta_xxxzz_xzzzz_1[i] * fe_0 + ta_xxxzz_xyzzzz_0[i] * pa_y[i] - ta_xxxzz_xyzzzz_1[i] * pc_y[i];

        ta_xxxyzz_xzzzzz_0[i] = ta_xxxzz_xzzzzz_0[i] * pa_y[i] - ta_xxxzz_xzzzzz_1[i] * pc_y[i];

        ta_xxxyzz_yyyyyy_0[i] = 2.0 * ta_xyzz_yyyyyy_0[i] * fe_0 - 2.0 * ta_xyzz_yyyyyy_1[i] * fe_0 + ta_xxyzz_yyyyyy_0[i] * pa_x[i] - ta_xxyzz_yyyyyy_1[i] * pc_x[i];

        ta_xxxyzz_yyyyyz_0[i] = 2.0 * ta_xyzz_yyyyyz_0[i] * fe_0 - 2.0 * ta_xyzz_yyyyyz_1[i] * fe_0 + ta_xxyzz_yyyyyz_0[i] * pa_x[i] - ta_xxyzz_yyyyyz_1[i] * pc_x[i];

        ta_xxxyzz_yyyyzz_0[i] = 2.0 * ta_xyzz_yyyyzz_0[i] * fe_0 - 2.0 * ta_xyzz_yyyyzz_1[i] * fe_0 + ta_xxyzz_yyyyzz_0[i] * pa_x[i] - ta_xxyzz_yyyyzz_1[i] * pc_x[i];

        ta_xxxyzz_yyyzzz_0[i] = 2.0 * ta_xyzz_yyyzzz_0[i] * fe_0 - 2.0 * ta_xyzz_yyyzzz_1[i] * fe_0 + ta_xxyzz_yyyzzz_0[i] * pa_x[i] - ta_xxyzz_yyyzzz_1[i] * pc_x[i];

        ta_xxxyzz_yyzzzz_0[i] = 2.0 * ta_xyzz_yyzzzz_0[i] * fe_0 - 2.0 * ta_xyzz_yyzzzz_1[i] * fe_0 + ta_xxyzz_yyzzzz_0[i] * pa_x[i] - ta_xxyzz_yyzzzz_1[i] * pc_x[i];

        ta_xxxyzz_yzzzzz_0[i] = 2.0 * ta_xyzz_yzzzzz_0[i] * fe_0 - 2.0 * ta_xyzz_yzzzzz_1[i] * fe_0 + ta_xxyzz_yzzzzz_0[i] * pa_x[i] - ta_xxyzz_yzzzzz_1[i] * pc_x[i];

        ta_xxxyzz_zzzzzz_0[i] = ta_xxxzz_zzzzzz_0[i] * pa_y[i] - ta_xxxzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 252-280 components of targeted buffer : II

    auto ta_xxxzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 252);

    auto ta_xxxzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 253);

    auto ta_xxxzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 254);

    auto ta_xxxzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 255);

    auto ta_xxxzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 256);

    auto ta_xxxzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 257);

    auto ta_xxxzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 258);

    auto ta_xxxzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 259);

    auto ta_xxxzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 260);

    auto ta_xxxzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 261);

    auto ta_xxxzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 262);

    auto ta_xxxzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 263);

    auto ta_xxxzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 264);

    auto ta_xxxzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 265);

    auto ta_xxxzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 266);

    auto ta_xxxzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 267);

    auto ta_xxxzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 268);

    auto ta_xxxzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 269);

    auto ta_xxxzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 270);

    auto ta_xxxzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 271);

    auto ta_xxxzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 272);

    auto ta_xxxzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 273);

    auto ta_xxxzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 274);

    auto ta_xxxzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 275);

    auto ta_xxxzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 276);

    auto ta_xxxzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 277);

    auto ta_xxxzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 278);

    auto ta_xxxzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 279);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xxxz_xxxxxx_0, ta_xxxz_xxxxxx_1, ta_xxxz_xxxxxy_0, ta_xxxz_xxxxxy_1, ta_xxxz_xxxxyy_0, ta_xxxz_xxxxyy_1, ta_xxxz_xxxyyy_0, ta_xxxz_xxxyyy_1, ta_xxxz_xxyyyy_0, ta_xxxz_xxyyyy_1, ta_xxxz_xyyyyy_0, ta_xxxz_xyyyyy_1, ta_xxxzz_xxxxxx_0, ta_xxxzz_xxxxxx_1, ta_xxxzz_xxxxxy_0, ta_xxxzz_xxxxxy_1, ta_xxxzz_xxxxyy_0, ta_xxxzz_xxxxyy_1, ta_xxxzz_xxxyyy_0, ta_xxxzz_xxxyyy_1, ta_xxxzz_xxyyyy_0, ta_xxxzz_xxyyyy_1, ta_xxxzz_xyyyyy_0, ta_xxxzz_xyyyyy_1, ta_xxxzzz_xxxxxx_0, ta_xxxzzz_xxxxxy_0, ta_xxxzzz_xxxxxz_0, ta_xxxzzz_xxxxyy_0, ta_xxxzzz_xxxxyz_0, ta_xxxzzz_xxxxzz_0, ta_xxxzzz_xxxyyy_0, ta_xxxzzz_xxxyyz_0, ta_xxxzzz_xxxyzz_0, ta_xxxzzz_xxxzzz_0, ta_xxxzzz_xxyyyy_0, ta_xxxzzz_xxyyyz_0, ta_xxxzzz_xxyyzz_0, ta_xxxzzz_xxyzzz_0, ta_xxxzzz_xxzzzz_0, ta_xxxzzz_xyyyyy_0, ta_xxxzzz_xyyyyz_0, ta_xxxzzz_xyyyzz_0, ta_xxxzzz_xyyzzz_0, ta_xxxzzz_xyzzzz_0, ta_xxxzzz_xzzzzz_0, ta_xxxzzz_yyyyyy_0, ta_xxxzzz_yyyyyz_0, ta_xxxzzz_yyyyzz_0, ta_xxxzzz_yyyzzz_0, ta_xxxzzz_yyzzzz_0, ta_xxxzzz_yzzzzz_0, ta_xxxzzz_zzzzzz_0, ta_xxzzz_xxxxxz_0, ta_xxzzz_xxxxxz_1, ta_xxzzz_xxxxyz_0, ta_xxzzz_xxxxyz_1, ta_xxzzz_xxxxz_0, ta_xxzzz_xxxxz_1, ta_xxzzz_xxxxzz_0, ta_xxzzz_xxxxzz_1, ta_xxzzz_xxxyyz_0, ta_xxzzz_xxxyyz_1, ta_xxzzz_xxxyz_0, ta_xxzzz_xxxyz_1, ta_xxzzz_xxxyzz_0, ta_xxzzz_xxxyzz_1, ta_xxzzz_xxxzz_0, ta_xxzzz_xxxzz_1, ta_xxzzz_xxxzzz_0, ta_xxzzz_xxxzzz_1, ta_xxzzz_xxyyyz_0, ta_xxzzz_xxyyyz_1, ta_xxzzz_xxyyz_0, ta_xxzzz_xxyyz_1, ta_xxzzz_xxyyzz_0, ta_xxzzz_xxyyzz_1, ta_xxzzz_xxyzz_0, ta_xxzzz_xxyzz_1, ta_xxzzz_xxyzzz_0, ta_xxzzz_xxyzzz_1, ta_xxzzz_xxzzz_0, ta_xxzzz_xxzzz_1, ta_xxzzz_xxzzzz_0, ta_xxzzz_xxzzzz_1, ta_xxzzz_xyyyyz_0, ta_xxzzz_xyyyyz_1, ta_xxzzz_xyyyz_0, ta_xxzzz_xyyyz_1, ta_xxzzz_xyyyzz_0, ta_xxzzz_xyyyzz_1, ta_xxzzz_xyyzz_0, ta_xxzzz_xyyzz_1, ta_xxzzz_xyyzzz_0, ta_xxzzz_xyyzzz_1, ta_xxzzz_xyzzz_0, ta_xxzzz_xyzzz_1, ta_xxzzz_xyzzzz_0, ta_xxzzz_xyzzzz_1, ta_xxzzz_xzzzz_0, ta_xxzzz_xzzzz_1, ta_xxzzz_xzzzzz_0, ta_xxzzz_xzzzzz_1, ta_xxzzz_yyyyyy_0, ta_xxzzz_yyyyyy_1, ta_xxzzz_yyyyyz_0, ta_xxzzz_yyyyyz_1, ta_xxzzz_yyyyz_0, ta_xxzzz_yyyyz_1, ta_xxzzz_yyyyzz_0, ta_xxzzz_yyyyzz_1, ta_xxzzz_yyyzz_0, ta_xxzzz_yyyzz_1, ta_xxzzz_yyyzzz_0, ta_xxzzz_yyyzzz_1, ta_xxzzz_yyzzz_0, ta_xxzzz_yyzzz_1, ta_xxzzz_yyzzzz_0, ta_xxzzz_yyzzzz_1, ta_xxzzz_yzzzz_0, ta_xxzzz_yzzzz_1, ta_xxzzz_yzzzzz_0, ta_xxzzz_yzzzzz_1, ta_xxzzz_zzzzz_0, ta_xxzzz_zzzzz_1, ta_xxzzz_zzzzzz_0, ta_xxzzz_zzzzzz_1, ta_xzzz_xxxxxz_0, ta_xzzz_xxxxxz_1, ta_xzzz_xxxxyz_0, ta_xzzz_xxxxyz_1, ta_xzzz_xxxxzz_0, ta_xzzz_xxxxzz_1, ta_xzzz_xxxyyz_0, ta_xzzz_xxxyyz_1, ta_xzzz_xxxyzz_0, ta_xzzz_xxxyzz_1, ta_xzzz_xxxzzz_0, ta_xzzz_xxxzzz_1, ta_xzzz_xxyyyz_0, ta_xzzz_xxyyyz_1, ta_xzzz_xxyyzz_0, ta_xzzz_xxyyzz_1, ta_xzzz_xxyzzz_0, ta_xzzz_xxyzzz_1, ta_xzzz_xxzzzz_0, ta_xzzz_xxzzzz_1, ta_xzzz_xyyyyz_0, ta_xzzz_xyyyyz_1, ta_xzzz_xyyyzz_0, ta_xzzz_xyyyzz_1, ta_xzzz_xyyzzz_0, ta_xzzz_xyyzzz_1, ta_xzzz_xyzzzz_0, ta_xzzz_xyzzzz_1, ta_xzzz_xzzzzz_0, ta_xzzz_xzzzzz_1, ta_xzzz_yyyyyy_0, ta_xzzz_yyyyyy_1, ta_xzzz_yyyyyz_0, ta_xzzz_yyyyyz_1, ta_xzzz_yyyyzz_0, ta_xzzz_yyyyzz_1, ta_xzzz_yyyzzz_0, ta_xzzz_yyyzzz_1, ta_xzzz_yyzzzz_0, ta_xzzz_yyzzzz_1, ta_xzzz_yzzzzz_0, ta_xzzz_yzzzzz_1, ta_xzzz_zzzzzz_0, ta_xzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzzz_xxxxxx_0[i] = 2.0 * ta_xxxz_xxxxxx_0[i] * fe_0 - 2.0 * ta_xxxz_xxxxxx_1[i] * fe_0 + ta_xxxzz_xxxxxx_0[i] * pa_z[i] - ta_xxxzz_xxxxxx_1[i] * pc_z[i];

        ta_xxxzzz_xxxxxy_0[i] = 2.0 * ta_xxxz_xxxxxy_0[i] * fe_0 - 2.0 * ta_xxxz_xxxxxy_1[i] * fe_0 + ta_xxxzz_xxxxxy_0[i] * pa_z[i] - ta_xxxzz_xxxxxy_1[i] * pc_z[i];

        ta_xxxzzz_xxxxxz_0[i] = 2.0 * ta_xzzz_xxxxxz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxxxz_1[i] * fe_0 + 5.0 * ta_xxzzz_xxxxz_0[i] * fe_0 - 5.0 * ta_xxzzz_xxxxz_1[i] * fe_0 + ta_xxzzz_xxxxxz_0[i] * pa_x[i] - ta_xxzzz_xxxxxz_1[i] * pc_x[i];

        ta_xxxzzz_xxxxyy_0[i] = 2.0 * ta_xxxz_xxxxyy_0[i] * fe_0 - 2.0 * ta_xxxz_xxxxyy_1[i] * fe_0 + ta_xxxzz_xxxxyy_0[i] * pa_z[i] - ta_xxxzz_xxxxyy_1[i] * pc_z[i];

        ta_xxxzzz_xxxxyz_0[i] = 2.0 * ta_xzzz_xxxxyz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_xxzzz_xxxyz_1[i] * fe_0 + ta_xxzzz_xxxxyz_0[i] * pa_x[i] - ta_xxzzz_xxxxyz_1[i] * pc_x[i];

        ta_xxxzzz_xxxxzz_0[i] = 2.0 * ta_xzzz_xxxxzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxxzz_1[i] * fe_0 + 4.0 * ta_xxzzz_xxxzz_0[i] * fe_0 - 4.0 * ta_xxzzz_xxxzz_1[i] * fe_0 + ta_xxzzz_xxxxzz_0[i] * pa_x[i] - ta_xxzzz_xxxxzz_1[i] * pc_x[i];

        ta_xxxzzz_xxxyyy_0[i] = 2.0 * ta_xxxz_xxxyyy_0[i] * fe_0 - 2.0 * ta_xxxz_xxxyyy_1[i] * fe_0 + ta_xxxzz_xxxyyy_0[i] * pa_z[i] - ta_xxxzz_xxxyyy_1[i] * pc_z[i];

        ta_xxxzzz_xxxyyz_0[i] = 2.0 * ta_xzzz_xxxyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxzzz_xxyyz_1[i] * fe_0 + ta_xxzzz_xxxyyz_0[i] * pa_x[i] - ta_xxzzz_xxxyyz_1[i] * pc_x[i];

        ta_xxxzzz_xxxyzz_0[i] = 2.0 * ta_xzzz_xxxyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_xxzzz_xxyzz_1[i] * fe_0 + ta_xxzzz_xxxyzz_0[i] * pa_x[i] - ta_xxzzz_xxxyzz_1[i] * pc_x[i];

        ta_xxxzzz_xxxzzz_0[i] = 2.0 * ta_xzzz_xxxzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_xxzzz_xxzzz_0[i] * fe_0 - 3.0 * ta_xxzzz_xxzzz_1[i] * fe_0 + ta_xxzzz_xxxzzz_0[i] * pa_x[i] - ta_xxzzz_xxxzzz_1[i] * pc_x[i];

        ta_xxxzzz_xxyyyy_0[i] = 2.0 * ta_xxxz_xxyyyy_0[i] * fe_0 - 2.0 * ta_xxxz_xxyyyy_1[i] * fe_0 + ta_xxxzz_xxyyyy_0[i] * pa_z[i] - ta_xxxzz_xxyyyy_1[i] * pc_z[i];

        ta_xxxzzz_xxyyyz_0[i] = 2.0 * ta_xzzz_xxyyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xxzzz_xyyyz_1[i] * fe_0 + ta_xxzzz_xxyyyz_0[i] * pa_x[i] - ta_xxzzz_xxyyyz_1[i] * pc_x[i];

        ta_xxxzzz_xxyyzz_0[i] = 2.0 * ta_xzzz_xxyyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xxzzz_xyyzz_1[i] * fe_0 + ta_xxzzz_xxyyzz_0[i] * pa_x[i] - ta_xxzzz_xxyyzz_1[i] * pc_x[i];

        ta_xxxzzz_xxyzzz_0[i] = 2.0 * ta_xzzz_xxyzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxzzz_xyzzz_1[i] * fe_0 + ta_xxzzz_xxyzzz_0[i] * pa_x[i] - ta_xxzzz_xxyzzz_1[i] * pc_x[i];

        ta_xxxzzz_xxzzzz_0[i] = 2.0 * ta_xzzz_xxzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxzzzz_1[i] * fe_0 + 2.0 * ta_xxzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_xxzzz_xzzzz_1[i] * fe_0 + ta_xxzzz_xxzzzz_0[i] * pa_x[i] - ta_xxzzz_xxzzzz_1[i] * pc_x[i];

        ta_xxxzzz_xyyyyy_0[i] = 2.0 * ta_xxxz_xyyyyy_0[i] * fe_0 - 2.0 * ta_xxxz_xyyyyy_1[i] * fe_0 + ta_xxxzz_xyyyyy_0[i] * pa_z[i] - ta_xxxzz_xyyyyy_1[i] * pc_z[i];

        ta_xxxzzz_xyyyyz_0[i] = 2.0 * ta_xzzz_xyyyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyyyz_1[i] * fe_0 + ta_xxzzz_yyyyz_0[i] * fe_0 - ta_xxzzz_yyyyz_1[i] * fe_0 + ta_xxzzz_xyyyyz_0[i] * pa_x[i] - ta_xxzzz_xyyyyz_1[i] * pc_x[i];

        ta_xxxzzz_xyyyzz_0[i] = 2.0 * ta_xzzz_xyyyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyyzz_1[i] * fe_0 + ta_xxzzz_yyyzz_0[i] * fe_0 - ta_xxzzz_yyyzz_1[i] * fe_0 + ta_xxzzz_xyyyzz_0[i] * pa_x[i] - ta_xxzzz_xyyyzz_1[i] * pc_x[i];

        ta_xxxzzz_xyyzzz_0[i] = 2.0 * ta_xzzz_xyyzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyzzz_1[i] * fe_0 + ta_xxzzz_yyzzz_0[i] * fe_0 - ta_xxzzz_yyzzz_1[i] * fe_0 + ta_xxzzz_xyyzzz_0[i] * pa_x[i] - ta_xxzzz_xyyzzz_1[i] * pc_x[i];

        ta_xxxzzz_xyzzzz_0[i] = 2.0 * ta_xzzz_xyzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyzzzz_1[i] * fe_0 + ta_xxzzz_yzzzz_0[i] * fe_0 - ta_xxzzz_yzzzz_1[i] * fe_0 + ta_xxzzz_xyzzzz_0[i] * pa_x[i] - ta_xxzzz_xyzzzz_1[i] * pc_x[i];

        ta_xxxzzz_xzzzzz_0[i] = 2.0 * ta_xzzz_xzzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xzzzzz_1[i] * fe_0 + ta_xxzzz_zzzzz_0[i] * fe_0 - ta_xxzzz_zzzzz_1[i] * fe_0 + ta_xxzzz_xzzzzz_0[i] * pa_x[i] - ta_xxzzz_xzzzzz_1[i] * pc_x[i];

        ta_xxxzzz_yyyyyy_0[i] = 2.0 * ta_xzzz_yyyyyy_0[i] * fe_0 - 2.0 * ta_xzzz_yyyyyy_1[i] * fe_0 + ta_xxzzz_yyyyyy_0[i] * pa_x[i] - ta_xxzzz_yyyyyy_1[i] * pc_x[i];

        ta_xxxzzz_yyyyyz_0[i] = 2.0 * ta_xzzz_yyyyyz_0[i] * fe_0 - 2.0 * ta_xzzz_yyyyyz_1[i] * fe_0 + ta_xxzzz_yyyyyz_0[i] * pa_x[i] - ta_xxzzz_yyyyyz_1[i] * pc_x[i];

        ta_xxxzzz_yyyyzz_0[i] = 2.0 * ta_xzzz_yyyyzz_0[i] * fe_0 - 2.0 * ta_xzzz_yyyyzz_1[i] * fe_0 + ta_xxzzz_yyyyzz_0[i] * pa_x[i] - ta_xxzzz_yyyyzz_1[i] * pc_x[i];

        ta_xxxzzz_yyyzzz_0[i] = 2.0 * ta_xzzz_yyyzzz_0[i] * fe_0 - 2.0 * ta_xzzz_yyyzzz_1[i] * fe_0 + ta_xxzzz_yyyzzz_0[i] * pa_x[i] - ta_xxzzz_yyyzzz_1[i] * pc_x[i];

        ta_xxxzzz_yyzzzz_0[i] = 2.0 * ta_xzzz_yyzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_yyzzzz_1[i] * fe_0 + ta_xxzzz_yyzzzz_0[i] * pa_x[i] - ta_xxzzz_yyzzzz_1[i] * pc_x[i];

        ta_xxxzzz_yzzzzz_0[i] = 2.0 * ta_xzzz_yzzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_yzzzzz_1[i] * fe_0 + ta_xxzzz_yzzzzz_0[i] * pa_x[i] - ta_xxzzz_yzzzzz_1[i] * pc_x[i];

        ta_xxxzzz_zzzzzz_0[i] = 2.0 * ta_xzzz_zzzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_zzzzzz_1[i] * fe_0 + ta_xxzzz_zzzzzz_0[i] * pa_x[i] - ta_xxzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 280-308 components of targeted buffer : II

    auto ta_xxyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 280);

    auto ta_xxyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 281);

    auto ta_xxyyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 282);

    auto ta_xxyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 283);

    auto ta_xxyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 284);

    auto ta_xxyyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 285);

    auto ta_xxyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 286);

    auto ta_xxyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 287);

    auto ta_xxyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 288);

    auto ta_xxyyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 289);

    auto ta_xxyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 290);

    auto ta_xxyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 291);

    auto ta_xxyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 292);

    auto ta_xxyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 293);

    auto ta_xxyyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 294);

    auto ta_xxyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 295);

    auto ta_xxyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 296);

    auto ta_xxyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 297);

    auto ta_xxyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 298);

    auto ta_xxyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 299);

    auto ta_xxyyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 300);

    auto ta_xxyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 301);

    auto ta_xxyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 302);

    auto ta_xxyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 303);

    auto ta_xxyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 304);

    auto ta_xxyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 305);

    auto ta_xxyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 306);

    auto ta_xxyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 307);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxyy_xxxxxx_0, ta_xxyy_xxxxxx_1, ta_xxyy_xxxxxz_0, ta_xxyy_xxxxxz_1, ta_xxyy_xxxxzz_0, ta_xxyy_xxxxzz_1, ta_xxyy_xxxzzz_0, ta_xxyy_xxxzzz_1, ta_xxyy_xxzzzz_0, ta_xxyy_xxzzzz_1, ta_xxyy_xzzzzz_0, ta_xxyy_xzzzzz_1, ta_xxyyy_xxxxxx_0, ta_xxyyy_xxxxxx_1, ta_xxyyy_xxxxxz_0, ta_xxyyy_xxxxxz_1, ta_xxyyy_xxxxzz_0, ta_xxyyy_xxxxzz_1, ta_xxyyy_xxxzzz_0, ta_xxyyy_xxxzzz_1, ta_xxyyy_xxzzzz_0, ta_xxyyy_xxzzzz_1, ta_xxyyy_xzzzzz_0, ta_xxyyy_xzzzzz_1, ta_xxyyyy_xxxxxx_0, ta_xxyyyy_xxxxxy_0, ta_xxyyyy_xxxxxz_0, ta_xxyyyy_xxxxyy_0, ta_xxyyyy_xxxxyz_0, ta_xxyyyy_xxxxzz_0, ta_xxyyyy_xxxyyy_0, ta_xxyyyy_xxxyyz_0, ta_xxyyyy_xxxyzz_0, ta_xxyyyy_xxxzzz_0, ta_xxyyyy_xxyyyy_0, ta_xxyyyy_xxyyyz_0, ta_xxyyyy_xxyyzz_0, ta_xxyyyy_xxyzzz_0, ta_xxyyyy_xxzzzz_0, ta_xxyyyy_xyyyyy_0, ta_xxyyyy_xyyyyz_0, ta_xxyyyy_xyyyzz_0, ta_xxyyyy_xyyzzz_0, ta_xxyyyy_xyzzzz_0, ta_xxyyyy_xzzzzz_0, ta_xxyyyy_yyyyyy_0, ta_xxyyyy_yyyyyz_0, ta_xxyyyy_yyyyzz_0, ta_xxyyyy_yyyzzz_0, ta_xxyyyy_yyzzzz_0, ta_xxyyyy_yzzzzz_0, ta_xxyyyy_zzzzzz_0, ta_xyyyy_xxxxxy_0, ta_xyyyy_xxxxxy_1, ta_xyyyy_xxxxy_0, ta_xyyyy_xxxxy_1, ta_xyyyy_xxxxyy_0, ta_xyyyy_xxxxyy_1, ta_xyyyy_xxxxyz_0, ta_xyyyy_xxxxyz_1, ta_xyyyy_xxxyy_0, ta_xyyyy_xxxyy_1, ta_xyyyy_xxxyyy_0, ta_xyyyy_xxxyyy_1, ta_xyyyy_xxxyyz_0, ta_xyyyy_xxxyyz_1, ta_xyyyy_xxxyz_0, ta_xyyyy_xxxyz_1, ta_xyyyy_xxxyzz_0, ta_xyyyy_xxxyzz_1, ta_xyyyy_xxyyy_0, ta_xyyyy_xxyyy_1, ta_xyyyy_xxyyyy_0, ta_xyyyy_xxyyyy_1, ta_xyyyy_xxyyyz_0, ta_xyyyy_xxyyyz_1, ta_xyyyy_xxyyz_0, ta_xyyyy_xxyyz_1, ta_xyyyy_xxyyzz_0, ta_xyyyy_xxyyzz_1, ta_xyyyy_xxyzz_0, ta_xyyyy_xxyzz_1, ta_xyyyy_xxyzzz_0, ta_xyyyy_xxyzzz_1, ta_xyyyy_xyyyy_0, ta_xyyyy_xyyyy_1, ta_xyyyy_xyyyyy_0, ta_xyyyy_xyyyyy_1, ta_xyyyy_xyyyyz_0, ta_xyyyy_xyyyyz_1, ta_xyyyy_xyyyz_0, ta_xyyyy_xyyyz_1, ta_xyyyy_xyyyzz_0, ta_xyyyy_xyyyzz_1, ta_xyyyy_xyyzz_0, ta_xyyyy_xyyzz_1, ta_xyyyy_xyyzzz_0, ta_xyyyy_xyyzzz_1, ta_xyyyy_xyzzz_0, ta_xyyyy_xyzzz_1, ta_xyyyy_xyzzzz_0, ta_xyyyy_xyzzzz_1, ta_xyyyy_yyyyy_0, ta_xyyyy_yyyyy_1, ta_xyyyy_yyyyyy_0, ta_xyyyy_yyyyyy_1, ta_xyyyy_yyyyyz_0, ta_xyyyy_yyyyyz_1, ta_xyyyy_yyyyz_0, ta_xyyyy_yyyyz_1, ta_xyyyy_yyyyzz_0, ta_xyyyy_yyyyzz_1, ta_xyyyy_yyyzz_0, ta_xyyyy_yyyzz_1, ta_xyyyy_yyyzzz_0, ta_xyyyy_yyyzzz_1, ta_xyyyy_yyzzz_0, ta_xyyyy_yyzzz_1, ta_xyyyy_yyzzzz_0, ta_xyyyy_yyzzzz_1, ta_xyyyy_yzzzz_0, ta_xyyyy_yzzzz_1, ta_xyyyy_yzzzzz_0, ta_xyyyy_yzzzzz_1, ta_xyyyy_zzzzzz_0, ta_xyyyy_zzzzzz_1, ta_yyyy_xxxxxy_0, ta_yyyy_xxxxxy_1, ta_yyyy_xxxxyy_0, ta_yyyy_xxxxyy_1, ta_yyyy_xxxxyz_0, ta_yyyy_xxxxyz_1, ta_yyyy_xxxyyy_0, ta_yyyy_xxxyyy_1, ta_yyyy_xxxyyz_0, ta_yyyy_xxxyyz_1, ta_yyyy_xxxyzz_0, ta_yyyy_xxxyzz_1, ta_yyyy_xxyyyy_0, ta_yyyy_xxyyyy_1, ta_yyyy_xxyyyz_0, ta_yyyy_xxyyyz_1, ta_yyyy_xxyyzz_0, ta_yyyy_xxyyzz_1, ta_yyyy_xxyzzz_0, ta_yyyy_xxyzzz_1, ta_yyyy_xyyyyy_0, ta_yyyy_xyyyyy_1, ta_yyyy_xyyyyz_0, ta_yyyy_xyyyyz_1, ta_yyyy_xyyyzz_0, ta_yyyy_xyyyzz_1, ta_yyyy_xyyzzz_0, ta_yyyy_xyyzzz_1, ta_yyyy_xyzzzz_0, ta_yyyy_xyzzzz_1, ta_yyyy_yyyyyy_0, ta_yyyy_yyyyyy_1, ta_yyyy_yyyyyz_0, ta_yyyy_yyyyyz_1, ta_yyyy_yyyyzz_0, ta_yyyy_yyyyzz_1, ta_yyyy_yyyzzz_0, ta_yyyy_yyyzzz_1, ta_yyyy_yyzzzz_0, ta_yyyy_yyzzzz_1, ta_yyyy_yzzzzz_0, ta_yyyy_yzzzzz_1, ta_yyyy_zzzzzz_0, ta_yyyy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyy_xxxxxx_0[i] = 3.0 * ta_xxyy_xxxxxx_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxxx_1[i] * fe_0 + ta_xxyyy_xxxxxx_0[i] * pa_y[i] - ta_xxyyy_xxxxxx_1[i] * pc_y[i];

        ta_xxyyyy_xxxxxy_0[i] = ta_yyyy_xxxxxy_0[i] * fe_0 - ta_yyyy_xxxxxy_1[i] * fe_0 + 5.0 * ta_xyyyy_xxxxy_0[i] * fe_0 - 5.0 * ta_xyyyy_xxxxy_1[i] * fe_0 + ta_xyyyy_xxxxxy_0[i] * pa_x[i] - ta_xyyyy_xxxxxy_1[i] * pc_x[i];

        ta_xxyyyy_xxxxxz_0[i] = 3.0 * ta_xxyy_xxxxxz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxxz_1[i] * fe_0 + ta_xxyyy_xxxxxz_0[i] * pa_y[i] - ta_xxyyy_xxxxxz_1[i] * pc_y[i];

        ta_xxyyyy_xxxxyy_0[i] = ta_yyyy_xxxxyy_0[i] * fe_0 - ta_yyyy_xxxxyy_1[i] * fe_0 + 4.0 * ta_xyyyy_xxxyy_0[i] * fe_0 - 4.0 * ta_xyyyy_xxxyy_1[i] * fe_0 + ta_xyyyy_xxxxyy_0[i] * pa_x[i] - ta_xyyyy_xxxxyy_1[i] * pc_x[i];

        ta_xxyyyy_xxxxyz_0[i] = ta_yyyy_xxxxyz_0[i] * fe_0 - ta_yyyy_xxxxyz_1[i] * fe_0 + 4.0 * ta_xyyyy_xxxyz_0[i] * fe_0 - 4.0 * ta_xyyyy_xxxyz_1[i] * fe_0 + ta_xyyyy_xxxxyz_0[i] * pa_x[i] - ta_xyyyy_xxxxyz_1[i] * pc_x[i];

        ta_xxyyyy_xxxxzz_0[i] = 3.0 * ta_xxyy_xxxxzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxzz_1[i] * fe_0 + ta_xxyyy_xxxxzz_0[i] * pa_y[i] - ta_xxyyy_xxxxzz_1[i] * pc_y[i];

        ta_xxyyyy_xxxyyy_0[i] = ta_yyyy_xxxyyy_0[i] * fe_0 - ta_yyyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_xyyyy_xxyyy_0[i] * fe_0 - 3.0 * ta_xyyyy_xxyyy_1[i] * fe_0 + ta_xyyyy_xxxyyy_0[i] * pa_x[i] - ta_xyyyy_xxxyyy_1[i] * pc_x[i];

        ta_xxyyyy_xxxyyz_0[i] = ta_yyyy_xxxyyz_0[i] * fe_0 - ta_yyyy_xxxyyz_1[i] * fe_0 + 3.0 * ta_xyyyy_xxyyz_0[i] * fe_0 - 3.0 * ta_xyyyy_xxyyz_1[i] * fe_0 + ta_xyyyy_xxxyyz_0[i] * pa_x[i] - ta_xyyyy_xxxyyz_1[i] * pc_x[i];

        ta_xxyyyy_xxxyzz_0[i] = ta_yyyy_xxxyzz_0[i] * fe_0 - ta_yyyy_xxxyzz_1[i] * fe_0 + 3.0 * ta_xyyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xyyyy_xxyzz_1[i] * fe_0 + ta_xyyyy_xxxyzz_0[i] * pa_x[i] - ta_xyyyy_xxxyzz_1[i] * pc_x[i];

        ta_xxyyyy_xxxzzz_0[i] = 3.0 * ta_xxyy_xxxzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxzzz_1[i] * fe_0 + ta_xxyyy_xxxzzz_0[i] * pa_y[i] - ta_xxyyy_xxxzzz_1[i] * pc_y[i];

        ta_xxyyyy_xxyyyy_0[i] = ta_yyyy_xxyyyy_0[i] * fe_0 - ta_yyyy_xxyyyy_1[i] * fe_0 + 2.0 * ta_xyyyy_xyyyy_0[i] * fe_0 - 2.0 * ta_xyyyy_xyyyy_1[i] * fe_0 + ta_xyyyy_xxyyyy_0[i] * pa_x[i] - ta_xyyyy_xxyyyy_1[i] * pc_x[i];

        ta_xxyyyy_xxyyyz_0[i] = ta_yyyy_xxyyyz_0[i] * fe_0 - ta_yyyy_xxyyyz_1[i] * fe_0 + 2.0 * ta_xyyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xyyyy_xyyyz_1[i] * fe_0 + ta_xyyyy_xxyyyz_0[i] * pa_x[i] - ta_xyyyy_xxyyyz_1[i] * pc_x[i];

        ta_xxyyyy_xxyyzz_0[i] = ta_yyyy_xxyyzz_0[i] * fe_0 - ta_yyyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_xyyyy_xyyzz_0[i] * fe_0 - 2.0 * ta_xyyyy_xyyzz_1[i] * fe_0 + ta_xyyyy_xxyyzz_0[i] * pa_x[i] - ta_xyyyy_xxyyzz_1[i] * pc_x[i];

        ta_xxyyyy_xxyzzz_0[i] = ta_yyyy_xxyzzz_0[i] * fe_0 - ta_yyyy_xxyzzz_1[i] * fe_0 + 2.0 * ta_xyyyy_xyzzz_0[i] * fe_0 - 2.0 * ta_xyyyy_xyzzz_1[i] * fe_0 + ta_xyyyy_xxyzzz_0[i] * pa_x[i] - ta_xyyyy_xxyzzz_1[i] * pc_x[i];

        ta_xxyyyy_xxzzzz_0[i] = 3.0 * ta_xxyy_xxzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxzzzz_1[i] * fe_0 + ta_xxyyy_xxzzzz_0[i] * pa_y[i] - ta_xxyyy_xxzzzz_1[i] * pc_y[i];

        ta_xxyyyy_xyyyyy_0[i] = ta_yyyy_xyyyyy_0[i] * fe_0 - ta_yyyy_xyyyyy_1[i] * fe_0 + ta_xyyyy_yyyyy_0[i] * fe_0 - ta_xyyyy_yyyyy_1[i] * fe_0 + ta_xyyyy_xyyyyy_0[i] * pa_x[i] - ta_xyyyy_xyyyyy_1[i] * pc_x[i];

        ta_xxyyyy_xyyyyz_0[i] = ta_yyyy_xyyyyz_0[i] * fe_0 - ta_yyyy_xyyyyz_1[i] * fe_0 + ta_xyyyy_yyyyz_0[i] * fe_0 - ta_xyyyy_yyyyz_1[i] * fe_0 + ta_xyyyy_xyyyyz_0[i] * pa_x[i] - ta_xyyyy_xyyyyz_1[i] * pc_x[i];

        ta_xxyyyy_xyyyzz_0[i] = ta_yyyy_xyyyzz_0[i] * fe_0 - ta_yyyy_xyyyzz_1[i] * fe_0 + ta_xyyyy_yyyzz_0[i] * fe_0 - ta_xyyyy_yyyzz_1[i] * fe_0 + ta_xyyyy_xyyyzz_0[i] * pa_x[i] - ta_xyyyy_xyyyzz_1[i] * pc_x[i];

        ta_xxyyyy_xyyzzz_0[i] = ta_yyyy_xyyzzz_0[i] * fe_0 - ta_yyyy_xyyzzz_1[i] * fe_0 + ta_xyyyy_yyzzz_0[i] * fe_0 - ta_xyyyy_yyzzz_1[i] * fe_0 + ta_xyyyy_xyyzzz_0[i] * pa_x[i] - ta_xyyyy_xyyzzz_1[i] * pc_x[i];

        ta_xxyyyy_xyzzzz_0[i] = ta_yyyy_xyzzzz_0[i] * fe_0 - ta_yyyy_xyzzzz_1[i] * fe_0 + ta_xyyyy_yzzzz_0[i] * fe_0 - ta_xyyyy_yzzzz_1[i] * fe_0 + ta_xyyyy_xyzzzz_0[i] * pa_x[i] - ta_xyyyy_xyzzzz_1[i] * pc_x[i];

        ta_xxyyyy_xzzzzz_0[i] = 3.0 * ta_xxyy_xzzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xzzzzz_1[i] * fe_0 + ta_xxyyy_xzzzzz_0[i] * pa_y[i] - ta_xxyyy_xzzzzz_1[i] * pc_y[i];

        ta_xxyyyy_yyyyyy_0[i] = ta_yyyy_yyyyyy_0[i] * fe_0 - ta_yyyy_yyyyyy_1[i] * fe_0 + ta_xyyyy_yyyyyy_0[i] * pa_x[i] - ta_xyyyy_yyyyyy_1[i] * pc_x[i];

        ta_xxyyyy_yyyyyz_0[i] = ta_yyyy_yyyyyz_0[i] * fe_0 - ta_yyyy_yyyyyz_1[i] * fe_0 + ta_xyyyy_yyyyyz_0[i] * pa_x[i] - ta_xyyyy_yyyyyz_1[i] * pc_x[i];

        ta_xxyyyy_yyyyzz_0[i] = ta_yyyy_yyyyzz_0[i] * fe_0 - ta_yyyy_yyyyzz_1[i] * fe_0 + ta_xyyyy_yyyyzz_0[i] * pa_x[i] - ta_xyyyy_yyyyzz_1[i] * pc_x[i];

        ta_xxyyyy_yyyzzz_0[i] = ta_yyyy_yyyzzz_0[i] * fe_0 - ta_yyyy_yyyzzz_1[i] * fe_0 + ta_xyyyy_yyyzzz_0[i] * pa_x[i] - ta_xyyyy_yyyzzz_1[i] * pc_x[i];

        ta_xxyyyy_yyzzzz_0[i] = ta_yyyy_yyzzzz_0[i] * fe_0 - ta_yyyy_yyzzzz_1[i] * fe_0 + ta_xyyyy_yyzzzz_0[i] * pa_x[i] - ta_xyyyy_yyzzzz_1[i] * pc_x[i];

        ta_xxyyyy_yzzzzz_0[i] = ta_yyyy_yzzzzz_0[i] * fe_0 - ta_yyyy_yzzzzz_1[i] * fe_0 + ta_xyyyy_yzzzzz_0[i] * pa_x[i] - ta_xyyyy_yzzzzz_1[i] * pc_x[i];

        ta_xxyyyy_zzzzzz_0[i] = ta_yyyy_zzzzzz_0[i] * fe_0 - ta_yyyy_zzzzzz_1[i] * fe_0 + ta_xyyyy_zzzzzz_0[i] * pa_x[i] - ta_xyyyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 308-336 components of targeted buffer : II

    auto ta_xxyyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 308);

    auto ta_xxyyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 309);

    auto ta_xxyyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 310);

    auto ta_xxyyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 311);

    auto ta_xxyyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 312);

    auto ta_xxyyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 313);

    auto ta_xxyyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 314);

    auto ta_xxyyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 315);

    auto ta_xxyyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 316);

    auto ta_xxyyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 317);

    auto ta_xxyyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 318);

    auto ta_xxyyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 319);

    auto ta_xxyyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 320);

    auto ta_xxyyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 321);

    auto ta_xxyyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 322);

    auto ta_xxyyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 323);

    auto ta_xxyyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 324);

    auto ta_xxyyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 325);

    auto ta_xxyyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 326);

    auto ta_xxyyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 327);

    auto ta_xxyyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 328);

    auto ta_xxyyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 329);

    auto ta_xxyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 330);

    auto ta_xxyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 331);

    auto ta_xxyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 332);

    auto ta_xxyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 333);

    auto ta_xxyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 334);

    auto ta_xxyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 335);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxyyy_xxxxxx_0, ta_xxyyy_xxxxxx_1, ta_xxyyy_xxxxxy_0, ta_xxyyy_xxxxxy_1, ta_xxyyy_xxxxy_0, ta_xxyyy_xxxxy_1, ta_xxyyy_xxxxyy_0, ta_xxyyy_xxxxyy_1, ta_xxyyy_xxxxyz_0, ta_xxyyy_xxxxyz_1, ta_xxyyy_xxxyy_0, ta_xxyyy_xxxyy_1, ta_xxyyy_xxxyyy_0, ta_xxyyy_xxxyyy_1, ta_xxyyy_xxxyyz_0, ta_xxyyy_xxxyyz_1, ta_xxyyy_xxxyz_0, ta_xxyyy_xxxyz_1, ta_xxyyy_xxxyzz_0, ta_xxyyy_xxxyzz_1, ta_xxyyy_xxyyy_0, ta_xxyyy_xxyyy_1, ta_xxyyy_xxyyyy_0, ta_xxyyy_xxyyyy_1, ta_xxyyy_xxyyyz_0, ta_xxyyy_xxyyyz_1, ta_xxyyy_xxyyz_0, ta_xxyyy_xxyyz_1, ta_xxyyy_xxyyzz_0, ta_xxyyy_xxyyzz_1, ta_xxyyy_xxyzz_0, ta_xxyyy_xxyzz_1, ta_xxyyy_xxyzzz_0, ta_xxyyy_xxyzzz_1, ta_xxyyy_xyyyy_0, ta_xxyyy_xyyyy_1, ta_xxyyy_xyyyyy_0, ta_xxyyy_xyyyyy_1, ta_xxyyy_xyyyyz_0, ta_xxyyy_xyyyyz_1, ta_xxyyy_xyyyz_0, ta_xxyyy_xyyyz_1, ta_xxyyy_xyyyzz_0, ta_xxyyy_xyyyzz_1, ta_xxyyy_xyyzz_0, ta_xxyyy_xyyzz_1, ta_xxyyy_xyyzzz_0, ta_xxyyy_xyyzzz_1, ta_xxyyy_xyzzz_0, ta_xxyyy_xyzzz_1, ta_xxyyy_xyzzzz_0, ta_xxyyy_xyzzzz_1, ta_xxyyy_yyyyyy_0, ta_xxyyy_yyyyyy_1, ta_xxyyyz_xxxxxx_0, ta_xxyyyz_xxxxxy_0, ta_xxyyyz_xxxxxz_0, ta_xxyyyz_xxxxyy_0, ta_xxyyyz_xxxxyz_0, ta_xxyyyz_xxxxzz_0, ta_xxyyyz_xxxyyy_0, ta_xxyyyz_xxxyyz_0, ta_xxyyyz_xxxyzz_0, ta_xxyyyz_xxxzzz_0, ta_xxyyyz_xxyyyy_0, ta_xxyyyz_xxyyyz_0, ta_xxyyyz_xxyyzz_0, ta_xxyyyz_xxyzzz_0, ta_xxyyyz_xxzzzz_0, ta_xxyyyz_xyyyyy_0, ta_xxyyyz_xyyyyz_0, ta_xxyyyz_xyyyzz_0, ta_xxyyyz_xyyzzz_0, ta_xxyyyz_xyzzzz_0, ta_xxyyyz_xzzzzz_0, ta_xxyyyz_yyyyyy_0, ta_xxyyyz_yyyyyz_0, ta_xxyyyz_yyyyzz_0, ta_xxyyyz_yyyzzz_0, ta_xxyyyz_yyzzzz_0, ta_xxyyyz_yzzzzz_0, ta_xxyyyz_zzzzzz_0, ta_xxyyz_xxxxxz_0, ta_xxyyz_xxxxxz_1, ta_xxyyz_xxxxzz_0, ta_xxyyz_xxxxzz_1, ta_xxyyz_xxxzzz_0, ta_xxyyz_xxxzzz_1, ta_xxyyz_xxzzzz_0, ta_xxyyz_xxzzzz_1, ta_xxyyz_xzzzzz_0, ta_xxyyz_xzzzzz_1, ta_xxyz_xxxxxz_0, ta_xxyz_xxxxxz_1, ta_xxyz_xxxxzz_0, ta_xxyz_xxxxzz_1, ta_xxyz_xxxzzz_0, ta_xxyz_xxxzzz_1, ta_xxyz_xxzzzz_0, ta_xxyz_xxzzzz_1, ta_xxyz_xzzzzz_0, ta_xxyz_xzzzzz_1, ta_xyyyz_yyyyyz_0, ta_xyyyz_yyyyyz_1, ta_xyyyz_yyyyzz_0, ta_xyyyz_yyyyzz_1, ta_xyyyz_yyyzzz_0, ta_xyyyz_yyyzzz_1, ta_xyyyz_yyzzzz_0, ta_xyyyz_yyzzzz_1, ta_xyyyz_yzzzzz_0, ta_xyyyz_yzzzzz_1, ta_xyyyz_zzzzzz_0, ta_xyyyz_zzzzzz_1, ta_yyyz_yyyyyz_0, ta_yyyz_yyyyyz_1, ta_yyyz_yyyyzz_0, ta_yyyz_yyyyzz_1, ta_yyyz_yyyzzz_0, ta_yyyz_yyyzzz_1, ta_yyyz_yyzzzz_0, ta_yyyz_yyzzzz_1, ta_yyyz_yzzzzz_0, ta_yyyz_yzzzzz_1, ta_yyyz_zzzzzz_0, ta_yyyz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyz_xxxxxx_0[i] = ta_xxyyy_xxxxxx_0[i] * pa_z[i] - ta_xxyyy_xxxxxx_1[i] * pc_z[i];

        ta_xxyyyz_xxxxxy_0[i] = ta_xxyyy_xxxxxy_0[i] * pa_z[i] - ta_xxyyy_xxxxxy_1[i] * pc_z[i];

        ta_xxyyyz_xxxxxz_0[i] = 2.0 * ta_xxyz_xxxxxz_0[i] * fe_0 - 2.0 * ta_xxyz_xxxxxz_1[i] * fe_0 + ta_xxyyz_xxxxxz_0[i] * pa_y[i] - ta_xxyyz_xxxxxz_1[i] * pc_y[i];

        ta_xxyyyz_xxxxyy_0[i] = ta_xxyyy_xxxxyy_0[i] * pa_z[i] - ta_xxyyy_xxxxyy_1[i] * pc_z[i];

        ta_xxyyyz_xxxxyz_0[i] = ta_xxyyy_xxxxy_0[i] * fe_0 - ta_xxyyy_xxxxy_1[i] * fe_0 + ta_xxyyy_xxxxyz_0[i] * pa_z[i] - ta_xxyyy_xxxxyz_1[i] * pc_z[i];

        ta_xxyyyz_xxxxzz_0[i] = 2.0 * ta_xxyz_xxxxzz_0[i] * fe_0 - 2.0 * ta_xxyz_xxxxzz_1[i] * fe_0 + ta_xxyyz_xxxxzz_0[i] * pa_y[i] - ta_xxyyz_xxxxzz_1[i] * pc_y[i];

        ta_xxyyyz_xxxyyy_0[i] = ta_xxyyy_xxxyyy_0[i] * pa_z[i] - ta_xxyyy_xxxyyy_1[i] * pc_z[i];

        ta_xxyyyz_xxxyyz_0[i] = ta_xxyyy_xxxyy_0[i] * fe_0 - ta_xxyyy_xxxyy_1[i] * fe_0 + ta_xxyyy_xxxyyz_0[i] * pa_z[i] - ta_xxyyy_xxxyyz_1[i] * pc_z[i];

        ta_xxyyyz_xxxyzz_0[i] = 2.0 * ta_xxyyy_xxxyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xxxyz_1[i] * fe_0 + ta_xxyyy_xxxyzz_0[i] * pa_z[i] - ta_xxyyy_xxxyzz_1[i] * pc_z[i];

        ta_xxyyyz_xxxzzz_0[i] = 2.0 * ta_xxyz_xxxzzz_0[i] * fe_0 - 2.0 * ta_xxyz_xxxzzz_1[i] * fe_0 + ta_xxyyz_xxxzzz_0[i] * pa_y[i] - ta_xxyyz_xxxzzz_1[i] * pc_y[i];

        ta_xxyyyz_xxyyyy_0[i] = ta_xxyyy_xxyyyy_0[i] * pa_z[i] - ta_xxyyy_xxyyyy_1[i] * pc_z[i];

        ta_xxyyyz_xxyyyz_0[i] = ta_xxyyy_xxyyy_0[i] * fe_0 - ta_xxyyy_xxyyy_1[i] * fe_0 + ta_xxyyy_xxyyyz_0[i] * pa_z[i] - ta_xxyyy_xxyyyz_1[i] * pc_z[i];

        ta_xxyyyz_xxyyzz_0[i] = 2.0 * ta_xxyyy_xxyyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xxyyz_1[i] * fe_0 + ta_xxyyy_xxyyzz_0[i] * pa_z[i] - ta_xxyyy_xxyyzz_1[i] * pc_z[i];

        ta_xxyyyz_xxyzzz_0[i] = 3.0 * ta_xxyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xxyyy_xxyzz_1[i] * fe_0 + ta_xxyyy_xxyzzz_0[i] * pa_z[i] - ta_xxyyy_xxyzzz_1[i] * pc_z[i];

        ta_xxyyyz_xxzzzz_0[i] = 2.0 * ta_xxyz_xxzzzz_0[i] * fe_0 - 2.0 * ta_xxyz_xxzzzz_1[i] * fe_0 + ta_xxyyz_xxzzzz_0[i] * pa_y[i] - ta_xxyyz_xxzzzz_1[i] * pc_y[i];

        ta_xxyyyz_xyyyyy_0[i] = ta_xxyyy_xyyyyy_0[i] * pa_z[i] - ta_xxyyy_xyyyyy_1[i] * pc_z[i];

        ta_xxyyyz_xyyyyz_0[i] = ta_xxyyy_xyyyy_0[i] * fe_0 - ta_xxyyy_xyyyy_1[i] * fe_0 + ta_xxyyy_xyyyyz_0[i] * pa_z[i] - ta_xxyyy_xyyyyz_1[i] * pc_z[i];

        ta_xxyyyz_xyyyzz_0[i] = 2.0 * ta_xxyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xyyyz_1[i] * fe_0 + ta_xxyyy_xyyyzz_0[i] * pa_z[i] - ta_xxyyy_xyyyzz_1[i] * pc_z[i];

        ta_xxyyyz_xyyzzz_0[i] = 3.0 * ta_xxyyy_xyyzz_0[i] * fe_0 - 3.0 * ta_xxyyy_xyyzz_1[i] * fe_0 + ta_xxyyy_xyyzzz_0[i] * pa_z[i] - ta_xxyyy_xyyzzz_1[i] * pc_z[i];

        ta_xxyyyz_xyzzzz_0[i] = 4.0 * ta_xxyyy_xyzzz_0[i] * fe_0 - 4.0 * ta_xxyyy_xyzzz_1[i] * fe_0 + ta_xxyyy_xyzzzz_0[i] * pa_z[i] - ta_xxyyy_xyzzzz_1[i] * pc_z[i];

        ta_xxyyyz_xzzzzz_0[i] = 2.0 * ta_xxyz_xzzzzz_0[i] * fe_0 - 2.0 * ta_xxyz_xzzzzz_1[i] * fe_0 + ta_xxyyz_xzzzzz_0[i] * pa_y[i] - ta_xxyyz_xzzzzz_1[i] * pc_y[i];

        ta_xxyyyz_yyyyyy_0[i] = ta_xxyyy_yyyyyy_0[i] * pa_z[i] - ta_xxyyy_yyyyyy_1[i] * pc_z[i];

        ta_xxyyyz_yyyyyz_0[i] = ta_yyyz_yyyyyz_0[i] * fe_0 - ta_yyyz_yyyyyz_1[i] * fe_0 + ta_xyyyz_yyyyyz_0[i] * pa_x[i] - ta_xyyyz_yyyyyz_1[i] * pc_x[i];

        ta_xxyyyz_yyyyzz_0[i] = ta_yyyz_yyyyzz_0[i] * fe_0 - ta_yyyz_yyyyzz_1[i] * fe_0 + ta_xyyyz_yyyyzz_0[i] * pa_x[i] - ta_xyyyz_yyyyzz_1[i] * pc_x[i];

        ta_xxyyyz_yyyzzz_0[i] = ta_yyyz_yyyzzz_0[i] * fe_0 - ta_yyyz_yyyzzz_1[i] * fe_0 + ta_xyyyz_yyyzzz_0[i] * pa_x[i] - ta_xyyyz_yyyzzz_1[i] * pc_x[i];

        ta_xxyyyz_yyzzzz_0[i] = ta_yyyz_yyzzzz_0[i] * fe_0 - ta_yyyz_yyzzzz_1[i] * fe_0 + ta_xyyyz_yyzzzz_0[i] * pa_x[i] - ta_xyyyz_yyzzzz_1[i] * pc_x[i];

        ta_xxyyyz_yzzzzz_0[i] = ta_yyyz_yzzzzz_0[i] * fe_0 - ta_yyyz_yzzzzz_1[i] * fe_0 + ta_xyyyz_yzzzzz_0[i] * pa_x[i] - ta_xyyyz_yzzzzz_1[i] * pc_x[i];

        ta_xxyyyz_zzzzzz_0[i] = ta_yyyz_zzzzzz_0[i] * fe_0 - ta_yyyz_zzzzzz_1[i] * fe_0 + ta_xyyyz_zzzzzz_0[i] * pa_x[i] - ta_xyyyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 336-364 components of targeted buffer : II

    auto ta_xxyyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 336);

    auto ta_xxyyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 337);

    auto ta_xxyyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 338);

    auto ta_xxyyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 339);

    auto ta_xxyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 340);

    auto ta_xxyyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 341);

    auto ta_xxyyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 342);

    auto ta_xxyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 343);

    auto ta_xxyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 344);

    auto ta_xxyyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 345);

    auto ta_xxyyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 346);

    auto ta_xxyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 347);

    auto ta_xxyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 348);

    auto ta_xxyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 349);

    auto ta_xxyyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 350);

    auto ta_xxyyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 351);

    auto ta_xxyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 352);

    auto ta_xxyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 353);

    auto ta_xxyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 354);

    auto ta_xxyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 355);

    auto ta_xxyyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 356);

    auto ta_xxyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 357);

    auto ta_xxyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 358);

    auto ta_xxyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 359);

    auto ta_xxyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 360);

    auto ta_xxyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 361);

    auto ta_xxyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 362);

    auto ta_xxyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 363);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxyy_xxxxxy_0, ta_xxyy_xxxxxy_1, ta_xxyy_xxxxyy_0, ta_xxyy_xxxxyy_1, ta_xxyy_xxxyyy_0, ta_xxyy_xxxyyy_1, ta_xxyy_xxyyyy_0, ta_xxyy_xxyyyy_1, ta_xxyy_xyyyyy_0, ta_xxyy_xyyyyy_1, ta_xxyyz_xxxxxy_0, ta_xxyyz_xxxxxy_1, ta_xxyyz_xxxxyy_0, ta_xxyyz_xxxxyy_1, ta_xxyyz_xxxyyy_0, ta_xxyyz_xxxyyy_1, ta_xxyyz_xxyyyy_0, ta_xxyyz_xxyyyy_1, ta_xxyyz_xyyyyy_0, ta_xxyyz_xyyyyy_1, ta_xxyyzz_xxxxxx_0, ta_xxyyzz_xxxxxy_0, ta_xxyyzz_xxxxxz_0, ta_xxyyzz_xxxxyy_0, ta_xxyyzz_xxxxyz_0, ta_xxyyzz_xxxxzz_0, ta_xxyyzz_xxxyyy_0, ta_xxyyzz_xxxyyz_0, ta_xxyyzz_xxxyzz_0, ta_xxyyzz_xxxzzz_0, ta_xxyyzz_xxyyyy_0, ta_xxyyzz_xxyyyz_0, ta_xxyyzz_xxyyzz_0, ta_xxyyzz_xxyzzz_0, ta_xxyyzz_xxzzzz_0, ta_xxyyzz_xyyyyy_0, ta_xxyyzz_xyyyyz_0, ta_xxyyzz_xyyyzz_0, ta_xxyyzz_xyyzzz_0, ta_xxyyzz_xyzzzz_0, ta_xxyyzz_xzzzzz_0, ta_xxyyzz_yyyyyy_0, ta_xxyyzz_yyyyyz_0, ta_xxyyzz_yyyyzz_0, ta_xxyyzz_yyyzzz_0, ta_xxyyzz_yyzzzz_0, ta_xxyyzz_yzzzzz_0, ta_xxyyzz_zzzzzz_0, ta_xxyzz_xxxxxx_0, ta_xxyzz_xxxxxx_1, ta_xxyzz_xxxxxz_0, ta_xxyzz_xxxxxz_1, ta_xxyzz_xxxxzz_0, ta_xxyzz_xxxxzz_1, ta_xxyzz_xxxzzz_0, ta_xxyzz_xxxzzz_1, ta_xxyzz_xxzzzz_0, ta_xxyzz_xxzzzz_1, ta_xxyzz_xzzzzz_0, ta_xxyzz_xzzzzz_1, ta_xxzz_xxxxxx_0, ta_xxzz_xxxxxx_1, ta_xxzz_xxxxxz_0, ta_xxzz_xxxxxz_1, ta_xxzz_xxxxzz_0, ta_xxzz_xxxxzz_1, ta_xxzz_xxxzzz_0, ta_xxzz_xxxzzz_1, ta_xxzz_xxzzzz_0, ta_xxzz_xxzzzz_1, ta_xxzz_xzzzzz_0, ta_xxzz_xzzzzz_1, ta_xyyzz_xxxxyz_0, ta_xyyzz_xxxxyz_1, ta_xyyzz_xxxyyz_0, ta_xyyzz_xxxyyz_1, ta_xyyzz_xxxyz_0, ta_xyyzz_xxxyz_1, ta_xyyzz_xxxyzz_0, ta_xyyzz_xxxyzz_1, ta_xyyzz_xxyyyz_0, ta_xyyzz_xxyyyz_1, ta_xyyzz_xxyyz_0, ta_xyyzz_xxyyz_1, ta_xyyzz_xxyyzz_0, ta_xyyzz_xxyyzz_1, ta_xyyzz_xxyzz_0, ta_xyyzz_xxyzz_1, ta_xyyzz_xxyzzz_0, ta_xyyzz_xxyzzz_1, ta_xyyzz_xyyyyz_0, ta_xyyzz_xyyyyz_1, ta_xyyzz_xyyyz_0, ta_xyyzz_xyyyz_1, ta_xyyzz_xyyyzz_0, ta_xyyzz_xyyyzz_1, ta_xyyzz_xyyzz_0, ta_xyyzz_xyyzz_1, ta_xyyzz_xyyzzz_0, ta_xyyzz_xyyzzz_1, ta_xyyzz_xyzzz_0, ta_xyyzz_xyzzz_1, ta_xyyzz_xyzzzz_0, ta_xyyzz_xyzzzz_1, ta_xyyzz_yyyyyy_0, ta_xyyzz_yyyyyy_1, ta_xyyzz_yyyyyz_0, ta_xyyzz_yyyyyz_1, ta_xyyzz_yyyyz_0, ta_xyyzz_yyyyz_1, ta_xyyzz_yyyyzz_0, ta_xyyzz_yyyyzz_1, ta_xyyzz_yyyzz_0, ta_xyyzz_yyyzz_1, ta_xyyzz_yyyzzz_0, ta_xyyzz_yyyzzz_1, ta_xyyzz_yyzzz_0, ta_xyyzz_yyzzz_1, ta_xyyzz_yyzzzz_0, ta_xyyzz_yyzzzz_1, ta_xyyzz_yzzzz_0, ta_xyyzz_yzzzz_1, ta_xyyzz_yzzzzz_0, ta_xyyzz_yzzzzz_1, ta_xyyzz_zzzzzz_0, ta_xyyzz_zzzzzz_1, ta_yyzz_xxxxyz_0, ta_yyzz_xxxxyz_1, ta_yyzz_xxxyyz_0, ta_yyzz_xxxyyz_1, ta_yyzz_xxxyzz_0, ta_yyzz_xxxyzz_1, ta_yyzz_xxyyyz_0, ta_yyzz_xxyyyz_1, ta_yyzz_xxyyzz_0, ta_yyzz_xxyyzz_1, ta_yyzz_xxyzzz_0, ta_yyzz_xxyzzz_1, ta_yyzz_xyyyyz_0, ta_yyzz_xyyyyz_1, ta_yyzz_xyyyzz_0, ta_yyzz_xyyyzz_1, ta_yyzz_xyyzzz_0, ta_yyzz_xyyzzz_1, ta_yyzz_xyzzzz_0, ta_yyzz_xyzzzz_1, ta_yyzz_yyyyyy_0, ta_yyzz_yyyyyy_1, ta_yyzz_yyyyyz_0, ta_yyzz_yyyyyz_1, ta_yyzz_yyyyzz_0, ta_yyzz_yyyyzz_1, ta_yyzz_yyyzzz_0, ta_yyzz_yyyzzz_1, ta_yyzz_yyzzzz_0, ta_yyzz_yyzzzz_1, ta_yyzz_yzzzzz_0, ta_yyzz_yzzzzz_1, ta_yyzz_zzzzzz_0, ta_yyzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyzz_xxxxxx_0[i] = ta_xxzz_xxxxxx_0[i] * fe_0 - ta_xxzz_xxxxxx_1[i] * fe_0 + ta_xxyzz_xxxxxx_0[i] * pa_y[i] - ta_xxyzz_xxxxxx_1[i] * pc_y[i];

        ta_xxyyzz_xxxxxy_0[i] = ta_xxyy_xxxxxy_0[i] * fe_0 - ta_xxyy_xxxxxy_1[i] * fe_0 + ta_xxyyz_xxxxxy_0[i] * pa_z[i] - ta_xxyyz_xxxxxy_1[i] * pc_z[i];

        ta_xxyyzz_xxxxxz_0[i] = ta_xxzz_xxxxxz_0[i] * fe_0 - ta_xxzz_xxxxxz_1[i] * fe_0 + ta_xxyzz_xxxxxz_0[i] * pa_y[i] - ta_xxyzz_xxxxxz_1[i] * pc_y[i];

        ta_xxyyzz_xxxxyy_0[i] = ta_xxyy_xxxxyy_0[i] * fe_0 - ta_xxyy_xxxxyy_1[i] * fe_0 + ta_xxyyz_xxxxyy_0[i] * pa_z[i] - ta_xxyyz_xxxxyy_1[i] * pc_z[i];

        ta_xxyyzz_xxxxyz_0[i] = ta_yyzz_xxxxyz_0[i] * fe_0 - ta_yyzz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xyyzz_xxxyz_0[i] * fe_0 - 4.0 * ta_xyyzz_xxxyz_1[i] * fe_0 + ta_xyyzz_xxxxyz_0[i] * pa_x[i] - ta_xyyzz_xxxxyz_1[i] * pc_x[i];

        ta_xxyyzz_xxxxzz_0[i] = ta_xxzz_xxxxzz_0[i] * fe_0 - ta_xxzz_xxxxzz_1[i] * fe_0 + ta_xxyzz_xxxxzz_0[i] * pa_y[i] - ta_xxyzz_xxxxzz_1[i] * pc_y[i];

        ta_xxyyzz_xxxyyy_0[i] = ta_xxyy_xxxyyy_0[i] * fe_0 - ta_xxyy_xxxyyy_1[i] * fe_0 + ta_xxyyz_xxxyyy_0[i] * pa_z[i] - ta_xxyyz_xxxyyy_1[i] * pc_z[i];

        ta_xxyyzz_xxxyyz_0[i] = ta_yyzz_xxxyyz_0[i] * fe_0 - ta_yyzz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xyyzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xyyzz_xxyyz_1[i] * fe_0 + ta_xyyzz_xxxyyz_0[i] * pa_x[i] - ta_xyyzz_xxxyyz_1[i] * pc_x[i];

        ta_xxyyzz_xxxyzz_0[i] = ta_yyzz_xxxyzz_0[i] * fe_0 - ta_yyzz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xyyzz_xxyzz_0[i] * fe_0 - 3.0 * ta_xyyzz_xxyzz_1[i] * fe_0 + ta_xyyzz_xxxyzz_0[i] * pa_x[i] - ta_xyyzz_xxxyzz_1[i] * pc_x[i];

        ta_xxyyzz_xxxzzz_0[i] = ta_xxzz_xxxzzz_0[i] * fe_0 - ta_xxzz_xxxzzz_1[i] * fe_0 + ta_xxyzz_xxxzzz_0[i] * pa_y[i] - ta_xxyzz_xxxzzz_1[i] * pc_y[i];

        ta_xxyyzz_xxyyyy_0[i] = ta_xxyy_xxyyyy_0[i] * fe_0 - ta_xxyy_xxyyyy_1[i] * fe_0 + ta_xxyyz_xxyyyy_0[i] * pa_z[i] - ta_xxyyz_xxyyyy_1[i] * pc_z[i];

        ta_xxyyzz_xxyyyz_0[i] = ta_yyzz_xxyyyz_0[i] * fe_0 - ta_yyzz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xyyzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xyyzz_xyyyz_1[i] * fe_0 + ta_xyyzz_xxyyyz_0[i] * pa_x[i] - ta_xyyzz_xxyyyz_1[i] * pc_x[i];

        ta_xxyyzz_xxyyzz_0[i] = ta_yyzz_xxyyzz_0[i] * fe_0 - ta_yyzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xyyzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xyyzz_xyyzz_1[i] * fe_0 + ta_xyyzz_xxyyzz_0[i] * pa_x[i] - ta_xyyzz_xxyyzz_1[i] * pc_x[i];

        ta_xxyyzz_xxyzzz_0[i] = ta_yyzz_xxyzzz_0[i] * fe_0 - ta_yyzz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xyyzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xyyzz_xyzzz_1[i] * fe_0 + ta_xyyzz_xxyzzz_0[i] * pa_x[i] - ta_xyyzz_xxyzzz_1[i] * pc_x[i];

        ta_xxyyzz_xxzzzz_0[i] = ta_xxzz_xxzzzz_0[i] * fe_0 - ta_xxzz_xxzzzz_1[i] * fe_0 + ta_xxyzz_xxzzzz_0[i] * pa_y[i] - ta_xxyzz_xxzzzz_1[i] * pc_y[i];

        ta_xxyyzz_xyyyyy_0[i] = ta_xxyy_xyyyyy_0[i] * fe_0 - ta_xxyy_xyyyyy_1[i] * fe_0 + ta_xxyyz_xyyyyy_0[i] * pa_z[i] - ta_xxyyz_xyyyyy_1[i] * pc_z[i];

        ta_xxyyzz_xyyyyz_0[i] = ta_yyzz_xyyyyz_0[i] * fe_0 - ta_yyzz_xyyyyz_1[i] * fe_0 + ta_xyyzz_yyyyz_0[i] * fe_0 - ta_xyyzz_yyyyz_1[i] * fe_0 + ta_xyyzz_xyyyyz_0[i] * pa_x[i] - ta_xyyzz_xyyyyz_1[i] * pc_x[i];

        ta_xxyyzz_xyyyzz_0[i] = ta_yyzz_xyyyzz_0[i] * fe_0 - ta_yyzz_xyyyzz_1[i] * fe_0 + ta_xyyzz_yyyzz_0[i] * fe_0 - ta_xyyzz_yyyzz_1[i] * fe_0 + ta_xyyzz_xyyyzz_0[i] * pa_x[i] - ta_xyyzz_xyyyzz_1[i] * pc_x[i];

        ta_xxyyzz_xyyzzz_0[i] = ta_yyzz_xyyzzz_0[i] * fe_0 - ta_yyzz_xyyzzz_1[i] * fe_0 + ta_xyyzz_yyzzz_0[i] * fe_0 - ta_xyyzz_yyzzz_1[i] * fe_0 + ta_xyyzz_xyyzzz_0[i] * pa_x[i] - ta_xyyzz_xyyzzz_1[i] * pc_x[i];

        ta_xxyyzz_xyzzzz_0[i] = ta_yyzz_xyzzzz_0[i] * fe_0 - ta_yyzz_xyzzzz_1[i] * fe_0 + ta_xyyzz_yzzzz_0[i] * fe_0 - ta_xyyzz_yzzzz_1[i] * fe_0 + ta_xyyzz_xyzzzz_0[i] * pa_x[i] - ta_xyyzz_xyzzzz_1[i] * pc_x[i];

        ta_xxyyzz_xzzzzz_0[i] = ta_xxzz_xzzzzz_0[i] * fe_0 - ta_xxzz_xzzzzz_1[i] * fe_0 + ta_xxyzz_xzzzzz_0[i] * pa_y[i] - ta_xxyzz_xzzzzz_1[i] * pc_y[i];

        ta_xxyyzz_yyyyyy_0[i] = ta_yyzz_yyyyyy_0[i] * fe_0 - ta_yyzz_yyyyyy_1[i] * fe_0 + ta_xyyzz_yyyyyy_0[i] * pa_x[i] - ta_xyyzz_yyyyyy_1[i] * pc_x[i];

        ta_xxyyzz_yyyyyz_0[i] = ta_yyzz_yyyyyz_0[i] * fe_0 - ta_yyzz_yyyyyz_1[i] * fe_0 + ta_xyyzz_yyyyyz_0[i] * pa_x[i] - ta_xyyzz_yyyyyz_1[i] * pc_x[i];

        ta_xxyyzz_yyyyzz_0[i] = ta_yyzz_yyyyzz_0[i] * fe_0 - ta_yyzz_yyyyzz_1[i] * fe_0 + ta_xyyzz_yyyyzz_0[i] * pa_x[i] - ta_xyyzz_yyyyzz_1[i] * pc_x[i];

        ta_xxyyzz_yyyzzz_0[i] = ta_yyzz_yyyzzz_0[i] * fe_0 - ta_yyzz_yyyzzz_1[i] * fe_0 + ta_xyyzz_yyyzzz_0[i] * pa_x[i] - ta_xyyzz_yyyzzz_1[i] * pc_x[i];

        ta_xxyyzz_yyzzzz_0[i] = ta_yyzz_yyzzzz_0[i] * fe_0 - ta_yyzz_yyzzzz_1[i] * fe_0 + ta_xyyzz_yyzzzz_0[i] * pa_x[i] - ta_xyyzz_yyzzzz_1[i] * pc_x[i];

        ta_xxyyzz_yzzzzz_0[i] = ta_yyzz_yzzzzz_0[i] * fe_0 - ta_yyzz_yzzzzz_1[i] * fe_0 + ta_xyyzz_yzzzzz_0[i] * pa_x[i] - ta_xyyzz_yzzzzz_1[i] * pc_x[i];

        ta_xxyyzz_zzzzzz_0[i] = ta_yyzz_zzzzzz_0[i] * fe_0 - ta_yyzz_zzzzzz_1[i] * fe_0 + ta_xyyzz_zzzzzz_0[i] * pa_x[i] - ta_xyyzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 364-392 components of targeted buffer : II

    auto ta_xxyzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 364);

    auto ta_xxyzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 365);

    auto ta_xxyzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 366);

    auto ta_xxyzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 367);

    auto ta_xxyzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 368);

    auto ta_xxyzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 369);

    auto ta_xxyzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 370);

    auto ta_xxyzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 371);

    auto ta_xxyzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 372);

    auto ta_xxyzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 373);

    auto ta_xxyzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 374);

    auto ta_xxyzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 375);

    auto ta_xxyzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 376);

    auto ta_xxyzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 377);

    auto ta_xxyzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 378);

    auto ta_xxyzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 379);

    auto ta_xxyzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 380);

    auto ta_xxyzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 381);

    auto ta_xxyzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 382);

    auto ta_xxyzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 383);

    auto ta_xxyzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 384);

    auto ta_xxyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 385);

    auto ta_xxyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 386);

    auto ta_xxyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 387);

    auto ta_xxyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 388);

    auto ta_xxyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 389);

    auto ta_xxyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 390);

    auto ta_xxyzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 391);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxyzzz_xxxxxx_0, ta_xxyzzz_xxxxxy_0, ta_xxyzzz_xxxxxz_0, ta_xxyzzz_xxxxyy_0, ta_xxyzzz_xxxxyz_0, ta_xxyzzz_xxxxzz_0, ta_xxyzzz_xxxyyy_0, ta_xxyzzz_xxxyyz_0, ta_xxyzzz_xxxyzz_0, ta_xxyzzz_xxxzzz_0, ta_xxyzzz_xxyyyy_0, ta_xxyzzz_xxyyyz_0, ta_xxyzzz_xxyyzz_0, ta_xxyzzz_xxyzzz_0, ta_xxyzzz_xxzzzz_0, ta_xxyzzz_xyyyyy_0, ta_xxyzzz_xyyyyz_0, ta_xxyzzz_xyyyzz_0, ta_xxyzzz_xyyzzz_0, ta_xxyzzz_xyzzzz_0, ta_xxyzzz_xzzzzz_0, ta_xxyzzz_yyyyyy_0, ta_xxyzzz_yyyyyz_0, ta_xxyzzz_yyyyzz_0, ta_xxyzzz_yyyzzz_0, ta_xxyzzz_yyzzzz_0, ta_xxyzzz_yzzzzz_0, ta_xxyzzz_zzzzzz_0, ta_xxzzz_xxxxx_0, ta_xxzzz_xxxxx_1, ta_xxzzz_xxxxxx_0, ta_xxzzz_xxxxxx_1, ta_xxzzz_xxxxxy_0, ta_xxzzz_xxxxxy_1, ta_xxzzz_xxxxxz_0, ta_xxzzz_xxxxxz_1, ta_xxzzz_xxxxy_0, ta_xxzzz_xxxxy_1, ta_xxzzz_xxxxyy_0, ta_xxzzz_xxxxyy_1, ta_xxzzz_xxxxyz_0, ta_xxzzz_xxxxyz_1, ta_xxzzz_xxxxz_0, ta_xxzzz_xxxxz_1, ta_xxzzz_xxxxzz_0, ta_xxzzz_xxxxzz_1, ta_xxzzz_xxxyy_0, ta_xxzzz_xxxyy_1, ta_xxzzz_xxxyyy_0, ta_xxzzz_xxxyyy_1, ta_xxzzz_xxxyyz_0, ta_xxzzz_xxxyyz_1, ta_xxzzz_xxxyz_0, ta_xxzzz_xxxyz_1, ta_xxzzz_xxxyzz_0, ta_xxzzz_xxxyzz_1, ta_xxzzz_xxxzz_0, ta_xxzzz_xxxzz_1, ta_xxzzz_xxxzzz_0, ta_xxzzz_xxxzzz_1, ta_xxzzz_xxyyy_0, ta_xxzzz_xxyyy_1, ta_xxzzz_xxyyyy_0, ta_xxzzz_xxyyyy_1, ta_xxzzz_xxyyyz_0, ta_xxzzz_xxyyyz_1, ta_xxzzz_xxyyz_0, ta_xxzzz_xxyyz_1, ta_xxzzz_xxyyzz_0, ta_xxzzz_xxyyzz_1, ta_xxzzz_xxyzz_0, ta_xxzzz_xxyzz_1, ta_xxzzz_xxyzzz_0, ta_xxzzz_xxyzzz_1, ta_xxzzz_xxzzz_0, ta_xxzzz_xxzzz_1, ta_xxzzz_xxzzzz_0, ta_xxzzz_xxzzzz_1, ta_xxzzz_xyyyy_0, ta_xxzzz_xyyyy_1, ta_xxzzz_xyyyyy_0, ta_xxzzz_xyyyyy_1, ta_xxzzz_xyyyyz_0, ta_xxzzz_xyyyyz_1, ta_xxzzz_xyyyz_0, ta_xxzzz_xyyyz_1, ta_xxzzz_xyyyzz_0, ta_xxzzz_xyyyzz_1, ta_xxzzz_xyyzz_0, ta_xxzzz_xyyzz_1, ta_xxzzz_xyyzzz_0, ta_xxzzz_xyyzzz_1, ta_xxzzz_xyzzz_0, ta_xxzzz_xyzzz_1, ta_xxzzz_xyzzzz_0, ta_xxzzz_xyzzzz_1, ta_xxzzz_xzzzz_0, ta_xxzzz_xzzzz_1, ta_xxzzz_xzzzzz_0, ta_xxzzz_xzzzzz_1, ta_xxzzz_zzzzzz_0, ta_xxzzz_zzzzzz_1, ta_xyzzz_yyyyyy_0, ta_xyzzz_yyyyyy_1, ta_xyzzz_yyyyyz_0, ta_xyzzz_yyyyyz_1, ta_xyzzz_yyyyzz_0, ta_xyzzz_yyyyzz_1, ta_xyzzz_yyyzzz_0, ta_xyzzz_yyyzzz_1, ta_xyzzz_yyzzzz_0, ta_xyzzz_yyzzzz_1, ta_xyzzz_yzzzzz_0, ta_xyzzz_yzzzzz_1, ta_yzzz_yyyyyy_0, ta_yzzz_yyyyyy_1, ta_yzzz_yyyyyz_0, ta_yzzz_yyyyyz_1, ta_yzzz_yyyyzz_0, ta_yzzz_yyyyzz_1, ta_yzzz_yyyzzz_0, ta_yzzz_yyyzzz_1, ta_yzzz_yyzzzz_0, ta_yzzz_yyzzzz_1, ta_yzzz_yzzzzz_0, ta_yzzz_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzzz_xxxxxx_0[i] = ta_xxzzz_xxxxxx_0[i] * pa_y[i] - ta_xxzzz_xxxxxx_1[i] * pc_y[i];

        ta_xxyzzz_xxxxxy_0[i] = ta_xxzzz_xxxxx_0[i] * fe_0 - ta_xxzzz_xxxxx_1[i] * fe_0 + ta_xxzzz_xxxxxy_0[i] * pa_y[i] - ta_xxzzz_xxxxxy_1[i] * pc_y[i];

        ta_xxyzzz_xxxxxz_0[i] = ta_xxzzz_xxxxxz_0[i] * pa_y[i] - ta_xxzzz_xxxxxz_1[i] * pc_y[i];

        ta_xxyzzz_xxxxyy_0[i] = 2.0 * ta_xxzzz_xxxxy_0[i] * fe_0 - 2.0 * ta_xxzzz_xxxxy_1[i] * fe_0 + ta_xxzzz_xxxxyy_0[i] * pa_y[i] - ta_xxzzz_xxxxyy_1[i] * pc_y[i];

        ta_xxyzzz_xxxxyz_0[i] = ta_xxzzz_xxxxz_0[i] * fe_0 - ta_xxzzz_xxxxz_1[i] * fe_0 + ta_xxzzz_xxxxyz_0[i] * pa_y[i] - ta_xxzzz_xxxxyz_1[i] * pc_y[i];

        ta_xxyzzz_xxxxzz_0[i] = ta_xxzzz_xxxxzz_0[i] * pa_y[i] - ta_xxzzz_xxxxzz_1[i] * pc_y[i];

        ta_xxyzzz_xxxyyy_0[i] = 3.0 * ta_xxzzz_xxxyy_0[i] * fe_0 - 3.0 * ta_xxzzz_xxxyy_1[i] * fe_0 + ta_xxzzz_xxxyyy_0[i] * pa_y[i] - ta_xxzzz_xxxyyy_1[i] * pc_y[i];

        ta_xxyzzz_xxxyyz_0[i] = 2.0 * ta_xxzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_xxzzz_xxxyz_1[i] * fe_0 + ta_xxzzz_xxxyyz_0[i] * pa_y[i] - ta_xxzzz_xxxyyz_1[i] * pc_y[i];

        ta_xxyzzz_xxxyzz_0[i] = ta_xxzzz_xxxzz_0[i] * fe_0 - ta_xxzzz_xxxzz_1[i] * fe_0 + ta_xxzzz_xxxyzz_0[i] * pa_y[i] - ta_xxzzz_xxxyzz_1[i] * pc_y[i];

        ta_xxyzzz_xxxzzz_0[i] = ta_xxzzz_xxxzzz_0[i] * pa_y[i] - ta_xxzzz_xxxzzz_1[i] * pc_y[i];

        ta_xxyzzz_xxyyyy_0[i] = 4.0 * ta_xxzzz_xxyyy_0[i] * fe_0 - 4.0 * ta_xxzzz_xxyyy_1[i] * fe_0 + ta_xxzzz_xxyyyy_0[i] * pa_y[i] - ta_xxzzz_xxyyyy_1[i] * pc_y[i];

        ta_xxyzzz_xxyyyz_0[i] = 3.0 * ta_xxzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxzzz_xxyyz_1[i] * fe_0 + ta_xxzzz_xxyyyz_0[i] * pa_y[i] - ta_xxzzz_xxyyyz_1[i] * pc_y[i];

        ta_xxyzzz_xxyyzz_0[i] = 2.0 * ta_xxzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_xxzzz_xxyzz_1[i] * fe_0 + ta_xxzzz_xxyyzz_0[i] * pa_y[i] - ta_xxzzz_xxyyzz_1[i] * pc_y[i];

        ta_xxyzzz_xxyzzz_0[i] = ta_xxzzz_xxzzz_0[i] * fe_0 - ta_xxzzz_xxzzz_1[i] * fe_0 + ta_xxzzz_xxyzzz_0[i] * pa_y[i] - ta_xxzzz_xxyzzz_1[i] * pc_y[i];

        ta_xxyzzz_xxzzzz_0[i] = ta_xxzzz_xxzzzz_0[i] * pa_y[i] - ta_xxzzz_xxzzzz_1[i] * pc_y[i];

        ta_xxyzzz_xyyyyy_0[i] = 5.0 * ta_xxzzz_xyyyy_0[i] * fe_0 - 5.0 * ta_xxzzz_xyyyy_1[i] * fe_0 + ta_xxzzz_xyyyyy_0[i] * pa_y[i] - ta_xxzzz_xyyyyy_1[i] * pc_y[i];

        ta_xxyzzz_xyyyyz_0[i] = 4.0 * ta_xxzzz_xyyyz_0[i] * fe_0 - 4.0 * ta_xxzzz_xyyyz_1[i] * fe_0 + ta_xxzzz_xyyyyz_0[i] * pa_y[i] - ta_xxzzz_xyyyyz_1[i] * pc_y[i];

        ta_xxyzzz_xyyyzz_0[i] = 3.0 * ta_xxzzz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxzzz_xyyzz_1[i] * fe_0 + ta_xxzzz_xyyyzz_0[i] * pa_y[i] - ta_xxzzz_xyyyzz_1[i] * pc_y[i];

        ta_xxyzzz_xyyzzz_0[i] = 2.0 * ta_xxzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxzzz_xyzzz_1[i] * fe_0 + ta_xxzzz_xyyzzz_0[i] * pa_y[i] - ta_xxzzz_xyyzzz_1[i] * pc_y[i];

        ta_xxyzzz_xyzzzz_0[i] = ta_xxzzz_xzzzz_0[i] * fe_0 - ta_xxzzz_xzzzz_1[i] * fe_0 + ta_xxzzz_xyzzzz_0[i] * pa_y[i] - ta_xxzzz_xyzzzz_1[i] * pc_y[i];

        ta_xxyzzz_xzzzzz_0[i] = ta_xxzzz_xzzzzz_0[i] * pa_y[i] - ta_xxzzz_xzzzzz_1[i] * pc_y[i];

        ta_xxyzzz_yyyyyy_0[i] = ta_yzzz_yyyyyy_0[i] * fe_0 - ta_yzzz_yyyyyy_1[i] * fe_0 + ta_xyzzz_yyyyyy_0[i] * pa_x[i] - ta_xyzzz_yyyyyy_1[i] * pc_x[i];

        ta_xxyzzz_yyyyyz_0[i] = ta_yzzz_yyyyyz_0[i] * fe_0 - ta_yzzz_yyyyyz_1[i] * fe_0 + ta_xyzzz_yyyyyz_0[i] * pa_x[i] - ta_xyzzz_yyyyyz_1[i] * pc_x[i];

        ta_xxyzzz_yyyyzz_0[i] = ta_yzzz_yyyyzz_0[i] * fe_0 - ta_yzzz_yyyyzz_1[i] * fe_0 + ta_xyzzz_yyyyzz_0[i] * pa_x[i] - ta_xyzzz_yyyyzz_1[i] * pc_x[i];

        ta_xxyzzz_yyyzzz_0[i] = ta_yzzz_yyyzzz_0[i] * fe_0 - ta_yzzz_yyyzzz_1[i] * fe_0 + ta_xyzzz_yyyzzz_0[i] * pa_x[i] - ta_xyzzz_yyyzzz_1[i] * pc_x[i];

        ta_xxyzzz_yyzzzz_0[i] = ta_yzzz_yyzzzz_0[i] * fe_0 - ta_yzzz_yyzzzz_1[i] * fe_0 + ta_xyzzz_yyzzzz_0[i] * pa_x[i] - ta_xyzzz_yyzzzz_1[i] * pc_x[i];

        ta_xxyzzz_yzzzzz_0[i] = ta_yzzz_yzzzzz_0[i] * fe_0 - ta_yzzz_yzzzzz_1[i] * fe_0 + ta_xyzzz_yzzzzz_0[i] * pa_x[i] - ta_xyzzz_yzzzzz_1[i] * pc_x[i];

        ta_xxyzzz_zzzzzz_0[i] = ta_xxzzz_zzzzzz_0[i] * pa_y[i] - ta_xxzzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 392-420 components of targeted buffer : II

    auto ta_xxzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 392);

    auto ta_xxzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 393);

    auto ta_xxzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 394);

    auto ta_xxzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 395);

    auto ta_xxzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 396);

    auto ta_xxzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 397);

    auto ta_xxzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 398);

    auto ta_xxzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 399);

    auto ta_xxzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 400);

    auto ta_xxzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 401);

    auto ta_xxzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 402);

    auto ta_xxzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 403);

    auto ta_xxzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 404);

    auto ta_xxzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 405);

    auto ta_xxzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 406);

    auto ta_xxzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 407);

    auto ta_xxzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 408);

    auto ta_xxzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 409);

    auto ta_xxzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 410);

    auto ta_xxzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 411);

    auto ta_xxzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 412);

    auto ta_xxzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 413);

    auto ta_xxzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 414);

    auto ta_xxzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 415);

    auto ta_xxzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 416);

    auto ta_xxzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 417);

    auto ta_xxzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 418);

    auto ta_xxzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 419);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xxzz_xxxxxx_0, ta_xxzz_xxxxxx_1, ta_xxzz_xxxxxy_0, ta_xxzz_xxxxxy_1, ta_xxzz_xxxxyy_0, ta_xxzz_xxxxyy_1, ta_xxzz_xxxyyy_0, ta_xxzz_xxxyyy_1, ta_xxzz_xxyyyy_0, ta_xxzz_xxyyyy_1, ta_xxzz_xyyyyy_0, ta_xxzz_xyyyyy_1, ta_xxzzz_xxxxxx_0, ta_xxzzz_xxxxxx_1, ta_xxzzz_xxxxxy_0, ta_xxzzz_xxxxxy_1, ta_xxzzz_xxxxyy_0, ta_xxzzz_xxxxyy_1, ta_xxzzz_xxxyyy_0, ta_xxzzz_xxxyyy_1, ta_xxzzz_xxyyyy_0, ta_xxzzz_xxyyyy_1, ta_xxzzz_xyyyyy_0, ta_xxzzz_xyyyyy_1, ta_xxzzzz_xxxxxx_0, ta_xxzzzz_xxxxxy_0, ta_xxzzzz_xxxxxz_0, ta_xxzzzz_xxxxyy_0, ta_xxzzzz_xxxxyz_0, ta_xxzzzz_xxxxzz_0, ta_xxzzzz_xxxyyy_0, ta_xxzzzz_xxxyyz_0, ta_xxzzzz_xxxyzz_0, ta_xxzzzz_xxxzzz_0, ta_xxzzzz_xxyyyy_0, ta_xxzzzz_xxyyyz_0, ta_xxzzzz_xxyyzz_0, ta_xxzzzz_xxyzzz_0, ta_xxzzzz_xxzzzz_0, ta_xxzzzz_xyyyyy_0, ta_xxzzzz_xyyyyz_0, ta_xxzzzz_xyyyzz_0, ta_xxzzzz_xyyzzz_0, ta_xxzzzz_xyzzzz_0, ta_xxzzzz_xzzzzz_0, ta_xxzzzz_yyyyyy_0, ta_xxzzzz_yyyyyz_0, ta_xxzzzz_yyyyzz_0, ta_xxzzzz_yyyzzz_0, ta_xxzzzz_yyzzzz_0, ta_xxzzzz_yzzzzz_0, ta_xxzzzz_zzzzzz_0, ta_xzzzz_xxxxxz_0, ta_xzzzz_xxxxxz_1, ta_xzzzz_xxxxyz_0, ta_xzzzz_xxxxyz_1, ta_xzzzz_xxxxz_0, ta_xzzzz_xxxxz_1, ta_xzzzz_xxxxzz_0, ta_xzzzz_xxxxzz_1, ta_xzzzz_xxxyyz_0, ta_xzzzz_xxxyyz_1, ta_xzzzz_xxxyz_0, ta_xzzzz_xxxyz_1, ta_xzzzz_xxxyzz_0, ta_xzzzz_xxxyzz_1, ta_xzzzz_xxxzz_0, ta_xzzzz_xxxzz_1, ta_xzzzz_xxxzzz_0, ta_xzzzz_xxxzzz_1, ta_xzzzz_xxyyyz_0, ta_xzzzz_xxyyyz_1, ta_xzzzz_xxyyz_0, ta_xzzzz_xxyyz_1, ta_xzzzz_xxyyzz_0, ta_xzzzz_xxyyzz_1, ta_xzzzz_xxyzz_0, ta_xzzzz_xxyzz_1, ta_xzzzz_xxyzzz_0, ta_xzzzz_xxyzzz_1, ta_xzzzz_xxzzz_0, ta_xzzzz_xxzzz_1, ta_xzzzz_xxzzzz_0, ta_xzzzz_xxzzzz_1, ta_xzzzz_xyyyyz_0, ta_xzzzz_xyyyyz_1, ta_xzzzz_xyyyz_0, ta_xzzzz_xyyyz_1, ta_xzzzz_xyyyzz_0, ta_xzzzz_xyyyzz_1, ta_xzzzz_xyyzz_0, ta_xzzzz_xyyzz_1, ta_xzzzz_xyyzzz_0, ta_xzzzz_xyyzzz_1, ta_xzzzz_xyzzz_0, ta_xzzzz_xyzzz_1, ta_xzzzz_xyzzzz_0, ta_xzzzz_xyzzzz_1, ta_xzzzz_xzzzz_0, ta_xzzzz_xzzzz_1, ta_xzzzz_xzzzzz_0, ta_xzzzz_xzzzzz_1, ta_xzzzz_yyyyyy_0, ta_xzzzz_yyyyyy_1, ta_xzzzz_yyyyyz_0, ta_xzzzz_yyyyyz_1, ta_xzzzz_yyyyz_0, ta_xzzzz_yyyyz_1, ta_xzzzz_yyyyzz_0, ta_xzzzz_yyyyzz_1, ta_xzzzz_yyyzz_0, ta_xzzzz_yyyzz_1, ta_xzzzz_yyyzzz_0, ta_xzzzz_yyyzzz_1, ta_xzzzz_yyzzz_0, ta_xzzzz_yyzzz_1, ta_xzzzz_yyzzzz_0, ta_xzzzz_yyzzzz_1, ta_xzzzz_yzzzz_0, ta_xzzzz_yzzzz_1, ta_xzzzz_yzzzzz_0, ta_xzzzz_yzzzzz_1, ta_xzzzz_zzzzz_0, ta_xzzzz_zzzzz_1, ta_xzzzz_zzzzzz_0, ta_xzzzz_zzzzzz_1, ta_zzzz_xxxxxz_0, ta_zzzz_xxxxxz_1, ta_zzzz_xxxxyz_0, ta_zzzz_xxxxyz_1, ta_zzzz_xxxxzz_0, ta_zzzz_xxxxzz_1, ta_zzzz_xxxyyz_0, ta_zzzz_xxxyyz_1, ta_zzzz_xxxyzz_0, ta_zzzz_xxxyzz_1, ta_zzzz_xxxzzz_0, ta_zzzz_xxxzzz_1, ta_zzzz_xxyyyz_0, ta_zzzz_xxyyyz_1, ta_zzzz_xxyyzz_0, ta_zzzz_xxyyzz_1, ta_zzzz_xxyzzz_0, ta_zzzz_xxyzzz_1, ta_zzzz_xxzzzz_0, ta_zzzz_xxzzzz_1, ta_zzzz_xyyyyz_0, ta_zzzz_xyyyyz_1, ta_zzzz_xyyyzz_0, ta_zzzz_xyyyzz_1, ta_zzzz_xyyzzz_0, ta_zzzz_xyyzzz_1, ta_zzzz_xyzzzz_0, ta_zzzz_xyzzzz_1, ta_zzzz_xzzzzz_0, ta_zzzz_xzzzzz_1, ta_zzzz_yyyyyy_0, ta_zzzz_yyyyyy_1, ta_zzzz_yyyyyz_0, ta_zzzz_yyyyyz_1, ta_zzzz_yyyyzz_0, ta_zzzz_yyyyzz_1, ta_zzzz_yyyzzz_0, ta_zzzz_yyyzzz_1, ta_zzzz_yyzzzz_0, ta_zzzz_yyzzzz_1, ta_zzzz_yzzzzz_0, ta_zzzz_yzzzzz_1, ta_zzzz_zzzzzz_0, ta_zzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzzz_xxxxxx_0[i] = 3.0 * ta_xxzz_xxxxxx_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxxx_1[i] * fe_0 + ta_xxzzz_xxxxxx_0[i] * pa_z[i] - ta_xxzzz_xxxxxx_1[i] * pc_z[i];

        ta_xxzzzz_xxxxxy_0[i] = 3.0 * ta_xxzz_xxxxxy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxxy_1[i] * fe_0 + ta_xxzzz_xxxxxy_0[i] * pa_z[i] - ta_xxzzz_xxxxxy_1[i] * pc_z[i];

        ta_xxzzzz_xxxxxz_0[i] = ta_zzzz_xxxxxz_0[i] * fe_0 - ta_zzzz_xxxxxz_1[i] * fe_0 + 5.0 * ta_xzzzz_xxxxz_0[i] * fe_0 - 5.0 * ta_xzzzz_xxxxz_1[i] * fe_0 + ta_xzzzz_xxxxxz_0[i] * pa_x[i] - ta_xzzzz_xxxxxz_1[i] * pc_x[i];

        ta_xxzzzz_xxxxyy_0[i] = 3.0 * ta_xxzz_xxxxyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxyy_1[i] * fe_0 + ta_xxzzz_xxxxyy_0[i] * pa_z[i] - ta_xxzzz_xxxxyy_1[i] * pc_z[i];

        ta_xxzzzz_xxxxyz_0[i] = ta_zzzz_xxxxyz_0[i] * fe_0 - ta_zzzz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xzzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_xzzzz_xxxyz_1[i] * fe_0 + ta_xzzzz_xxxxyz_0[i] * pa_x[i] - ta_xzzzz_xxxxyz_1[i] * pc_x[i];

        ta_xxzzzz_xxxxzz_0[i] = ta_zzzz_xxxxzz_0[i] * fe_0 - ta_zzzz_xxxxzz_1[i] * fe_0 + 4.0 * ta_xzzzz_xxxzz_0[i] * fe_0 - 4.0 * ta_xzzzz_xxxzz_1[i] * fe_0 + ta_xzzzz_xxxxzz_0[i] * pa_x[i] - ta_xzzzz_xxxxzz_1[i] * pc_x[i];

        ta_xxzzzz_xxxyyy_0[i] = 3.0 * ta_xxzz_xxxyyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxyyy_1[i] * fe_0 + ta_xxzzz_xxxyyy_0[i] * pa_z[i] - ta_xxzzz_xxxyyy_1[i] * pc_z[i];

        ta_xxzzzz_xxxyyz_0[i] = ta_zzzz_xxxyyz_0[i] * fe_0 - ta_zzzz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xzzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xzzzz_xxyyz_1[i] * fe_0 + ta_xzzzz_xxxyyz_0[i] * pa_x[i] - ta_xzzzz_xxxyyz_1[i] * pc_x[i];

        ta_xxzzzz_xxxyzz_0[i] = ta_zzzz_xxxyzz_0[i] * fe_0 - ta_zzzz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xzzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_xzzzz_xxyzz_1[i] * fe_0 + ta_xzzzz_xxxyzz_0[i] * pa_x[i] - ta_xzzzz_xxxyzz_1[i] * pc_x[i];

        ta_xxzzzz_xxxzzz_0[i] = ta_zzzz_xxxzzz_0[i] * fe_0 - ta_zzzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_xzzzz_xxzzz_0[i] * fe_0 - 3.0 * ta_xzzzz_xxzzz_1[i] * fe_0 + ta_xzzzz_xxxzzz_0[i] * pa_x[i] - ta_xzzzz_xxxzzz_1[i] * pc_x[i];

        ta_xxzzzz_xxyyyy_0[i] = 3.0 * ta_xxzz_xxyyyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxyyyy_1[i] * fe_0 + ta_xxzzz_xxyyyy_0[i] * pa_z[i] - ta_xxzzz_xxyyyy_1[i] * pc_z[i];

        ta_xxzzzz_xxyyyz_0[i] = ta_zzzz_xxyyyz_0[i] * fe_0 - ta_zzzz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xzzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xzzzz_xyyyz_1[i] * fe_0 + ta_xzzzz_xxyyyz_0[i] * pa_x[i] - ta_xzzzz_xxyyyz_1[i] * pc_x[i];

        ta_xxzzzz_xxyyzz_0[i] = ta_zzzz_xxyyzz_0[i] * fe_0 - ta_zzzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xzzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xzzzz_xyyzz_1[i] * fe_0 + ta_xzzzz_xxyyzz_0[i] * pa_x[i] - ta_xzzzz_xxyyzz_1[i] * pc_x[i];

        ta_xxzzzz_xxyzzz_0[i] = ta_zzzz_xxyzzz_0[i] * fe_0 - ta_zzzz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xzzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xzzzz_xyzzz_1[i] * fe_0 + ta_xzzzz_xxyzzz_0[i] * pa_x[i] - ta_xzzzz_xxyzzz_1[i] * pc_x[i];

        ta_xxzzzz_xxzzzz_0[i] = ta_zzzz_xxzzzz_0[i] * fe_0 - ta_zzzz_xxzzzz_1[i] * fe_0 + 2.0 * ta_xzzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_xzzzz_xzzzz_1[i] * fe_0 + ta_xzzzz_xxzzzz_0[i] * pa_x[i] - ta_xzzzz_xxzzzz_1[i] * pc_x[i];

        ta_xxzzzz_xyyyyy_0[i] = 3.0 * ta_xxzz_xyyyyy_0[i] * fe_0 - 3.0 * ta_xxzz_xyyyyy_1[i] * fe_0 + ta_xxzzz_xyyyyy_0[i] * pa_z[i] - ta_xxzzz_xyyyyy_1[i] * pc_z[i];

        ta_xxzzzz_xyyyyz_0[i] = ta_zzzz_xyyyyz_0[i] * fe_0 - ta_zzzz_xyyyyz_1[i] * fe_0 + ta_xzzzz_yyyyz_0[i] * fe_0 - ta_xzzzz_yyyyz_1[i] * fe_0 + ta_xzzzz_xyyyyz_0[i] * pa_x[i] - ta_xzzzz_xyyyyz_1[i] * pc_x[i];

        ta_xxzzzz_xyyyzz_0[i] = ta_zzzz_xyyyzz_0[i] * fe_0 - ta_zzzz_xyyyzz_1[i] * fe_0 + ta_xzzzz_yyyzz_0[i] * fe_0 - ta_xzzzz_yyyzz_1[i] * fe_0 + ta_xzzzz_xyyyzz_0[i] * pa_x[i] - ta_xzzzz_xyyyzz_1[i] * pc_x[i];

        ta_xxzzzz_xyyzzz_0[i] = ta_zzzz_xyyzzz_0[i] * fe_0 - ta_zzzz_xyyzzz_1[i] * fe_0 + ta_xzzzz_yyzzz_0[i] * fe_0 - ta_xzzzz_yyzzz_1[i] * fe_0 + ta_xzzzz_xyyzzz_0[i] * pa_x[i] - ta_xzzzz_xyyzzz_1[i] * pc_x[i];

        ta_xxzzzz_xyzzzz_0[i] = ta_zzzz_xyzzzz_0[i] * fe_0 - ta_zzzz_xyzzzz_1[i] * fe_0 + ta_xzzzz_yzzzz_0[i] * fe_0 - ta_xzzzz_yzzzz_1[i] * fe_0 + ta_xzzzz_xyzzzz_0[i] * pa_x[i] - ta_xzzzz_xyzzzz_1[i] * pc_x[i];

        ta_xxzzzz_xzzzzz_0[i] = ta_zzzz_xzzzzz_0[i] * fe_0 - ta_zzzz_xzzzzz_1[i] * fe_0 + ta_xzzzz_zzzzz_0[i] * fe_0 - ta_xzzzz_zzzzz_1[i] * fe_0 + ta_xzzzz_xzzzzz_0[i] * pa_x[i] - ta_xzzzz_xzzzzz_1[i] * pc_x[i];

        ta_xxzzzz_yyyyyy_0[i] = ta_zzzz_yyyyyy_0[i] * fe_0 - ta_zzzz_yyyyyy_1[i] * fe_0 + ta_xzzzz_yyyyyy_0[i] * pa_x[i] - ta_xzzzz_yyyyyy_1[i] * pc_x[i];

        ta_xxzzzz_yyyyyz_0[i] = ta_zzzz_yyyyyz_0[i] * fe_0 - ta_zzzz_yyyyyz_1[i] * fe_0 + ta_xzzzz_yyyyyz_0[i] * pa_x[i] - ta_xzzzz_yyyyyz_1[i] * pc_x[i];

        ta_xxzzzz_yyyyzz_0[i] = ta_zzzz_yyyyzz_0[i] * fe_0 - ta_zzzz_yyyyzz_1[i] * fe_0 + ta_xzzzz_yyyyzz_0[i] * pa_x[i] - ta_xzzzz_yyyyzz_1[i] * pc_x[i];

        ta_xxzzzz_yyyzzz_0[i] = ta_zzzz_yyyzzz_0[i] * fe_0 - ta_zzzz_yyyzzz_1[i] * fe_0 + ta_xzzzz_yyyzzz_0[i] * pa_x[i] - ta_xzzzz_yyyzzz_1[i] * pc_x[i];

        ta_xxzzzz_yyzzzz_0[i] = ta_zzzz_yyzzzz_0[i] * fe_0 - ta_zzzz_yyzzzz_1[i] * fe_0 + ta_xzzzz_yyzzzz_0[i] * pa_x[i] - ta_xzzzz_yyzzzz_1[i] * pc_x[i];

        ta_xxzzzz_yzzzzz_0[i] = ta_zzzz_yzzzzz_0[i] * fe_0 - ta_zzzz_yzzzzz_1[i] * fe_0 + ta_xzzzz_yzzzzz_0[i] * pa_x[i] - ta_xzzzz_yzzzzz_1[i] * pc_x[i];

        ta_xxzzzz_zzzzzz_0[i] = ta_zzzz_zzzzzz_0[i] * fe_0 - ta_zzzz_zzzzzz_1[i] * fe_0 + ta_xzzzz_zzzzzz_0[i] * pa_x[i] - ta_xzzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 420-448 components of targeted buffer : II

    auto ta_xyyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 420);

    auto ta_xyyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 421);

    auto ta_xyyyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 422);

    auto ta_xyyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 423);

    auto ta_xyyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 424);

    auto ta_xyyyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 425);

    auto ta_xyyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 426);

    auto ta_xyyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 427);

    auto ta_xyyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 428);

    auto ta_xyyyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 429);

    auto ta_xyyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 430);

    auto ta_xyyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 431);

    auto ta_xyyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 432);

    auto ta_xyyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 433);

    auto ta_xyyyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 434);

    auto ta_xyyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 435);

    auto ta_xyyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 436);

    auto ta_xyyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 437);

    auto ta_xyyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 438);

    auto ta_xyyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 439);

    auto ta_xyyyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 440);

    auto ta_xyyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 441);

    auto ta_xyyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 442);

    auto ta_xyyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 443);

    auto ta_xyyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 444);

    auto ta_xyyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 445);

    auto ta_xyyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 446);

    auto ta_xyyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 447);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyyyyy_xxxxxx_0, ta_xyyyyy_xxxxxy_0, ta_xyyyyy_xxxxxz_0, ta_xyyyyy_xxxxyy_0, ta_xyyyyy_xxxxyz_0, ta_xyyyyy_xxxxzz_0, ta_xyyyyy_xxxyyy_0, ta_xyyyyy_xxxyyz_0, ta_xyyyyy_xxxyzz_0, ta_xyyyyy_xxxzzz_0, ta_xyyyyy_xxyyyy_0, ta_xyyyyy_xxyyyz_0, ta_xyyyyy_xxyyzz_0, ta_xyyyyy_xxyzzz_0, ta_xyyyyy_xxzzzz_0, ta_xyyyyy_xyyyyy_0, ta_xyyyyy_xyyyyz_0, ta_xyyyyy_xyyyzz_0, ta_xyyyyy_xyyzzz_0, ta_xyyyyy_xyzzzz_0, ta_xyyyyy_xzzzzz_0, ta_xyyyyy_yyyyyy_0, ta_xyyyyy_yyyyyz_0, ta_xyyyyy_yyyyzz_0, ta_xyyyyy_yyyzzz_0, ta_xyyyyy_yyzzzz_0, ta_xyyyyy_yzzzzz_0, ta_xyyyyy_zzzzzz_0, ta_yyyyy_xxxxx_0, ta_yyyyy_xxxxx_1, ta_yyyyy_xxxxxx_0, ta_yyyyy_xxxxxx_1, ta_yyyyy_xxxxxy_0, ta_yyyyy_xxxxxy_1, ta_yyyyy_xxxxxz_0, ta_yyyyy_xxxxxz_1, ta_yyyyy_xxxxy_0, ta_yyyyy_xxxxy_1, ta_yyyyy_xxxxyy_0, ta_yyyyy_xxxxyy_1, ta_yyyyy_xxxxyz_0, ta_yyyyy_xxxxyz_1, ta_yyyyy_xxxxz_0, ta_yyyyy_xxxxz_1, ta_yyyyy_xxxxzz_0, ta_yyyyy_xxxxzz_1, ta_yyyyy_xxxyy_0, ta_yyyyy_xxxyy_1, ta_yyyyy_xxxyyy_0, ta_yyyyy_xxxyyy_1, ta_yyyyy_xxxyyz_0, ta_yyyyy_xxxyyz_1, ta_yyyyy_xxxyz_0, ta_yyyyy_xxxyz_1, ta_yyyyy_xxxyzz_0, ta_yyyyy_xxxyzz_1, ta_yyyyy_xxxzz_0, ta_yyyyy_xxxzz_1, ta_yyyyy_xxxzzz_0, ta_yyyyy_xxxzzz_1, ta_yyyyy_xxyyy_0, ta_yyyyy_xxyyy_1, ta_yyyyy_xxyyyy_0, ta_yyyyy_xxyyyy_1, ta_yyyyy_xxyyyz_0, ta_yyyyy_xxyyyz_1, ta_yyyyy_xxyyz_0, ta_yyyyy_xxyyz_1, ta_yyyyy_xxyyzz_0, ta_yyyyy_xxyyzz_1, ta_yyyyy_xxyzz_0, ta_yyyyy_xxyzz_1, ta_yyyyy_xxyzzz_0, ta_yyyyy_xxyzzz_1, ta_yyyyy_xxzzz_0, ta_yyyyy_xxzzz_1, ta_yyyyy_xxzzzz_0, ta_yyyyy_xxzzzz_1, ta_yyyyy_xyyyy_0, ta_yyyyy_xyyyy_1, ta_yyyyy_xyyyyy_0, ta_yyyyy_xyyyyy_1, ta_yyyyy_xyyyyz_0, ta_yyyyy_xyyyyz_1, ta_yyyyy_xyyyz_0, ta_yyyyy_xyyyz_1, ta_yyyyy_xyyyzz_0, ta_yyyyy_xyyyzz_1, ta_yyyyy_xyyzz_0, ta_yyyyy_xyyzz_1, ta_yyyyy_xyyzzz_0, ta_yyyyy_xyyzzz_1, ta_yyyyy_xyzzz_0, ta_yyyyy_xyzzz_1, ta_yyyyy_xyzzzz_0, ta_yyyyy_xyzzzz_1, ta_yyyyy_xzzzz_0, ta_yyyyy_xzzzz_1, ta_yyyyy_xzzzzz_0, ta_yyyyy_xzzzzz_1, ta_yyyyy_yyyyy_0, ta_yyyyy_yyyyy_1, ta_yyyyy_yyyyyy_0, ta_yyyyy_yyyyyy_1, ta_yyyyy_yyyyyz_0, ta_yyyyy_yyyyyz_1, ta_yyyyy_yyyyz_0, ta_yyyyy_yyyyz_1, ta_yyyyy_yyyyzz_0, ta_yyyyy_yyyyzz_1, ta_yyyyy_yyyzz_0, ta_yyyyy_yyyzz_1, ta_yyyyy_yyyzzz_0, ta_yyyyy_yyyzzz_1, ta_yyyyy_yyzzz_0, ta_yyyyy_yyzzz_1, ta_yyyyy_yyzzzz_0, ta_yyyyy_yyzzzz_1, ta_yyyyy_yzzzz_0, ta_yyyyy_yzzzz_1, ta_yyyyy_yzzzzz_0, ta_yyyyy_yzzzzz_1, ta_yyyyy_zzzzz_0, ta_yyyyy_zzzzz_1, ta_yyyyy_zzzzzz_0, ta_yyyyy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyy_xxxxxx_0[i] = 6.0 * ta_yyyyy_xxxxx_0[i] * fe_0 - 6.0 * ta_yyyyy_xxxxx_1[i] * fe_0 + ta_yyyyy_xxxxxx_0[i] * pa_x[i] - ta_yyyyy_xxxxxx_1[i] * pc_x[i];

        ta_xyyyyy_xxxxxy_0[i] = 5.0 * ta_yyyyy_xxxxy_0[i] * fe_0 - 5.0 * ta_yyyyy_xxxxy_1[i] * fe_0 + ta_yyyyy_xxxxxy_0[i] * pa_x[i] - ta_yyyyy_xxxxxy_1[i] * pc_x[i];

        ta_xyyyyy_xxxxxz_0[i] = 5.0 * ta_yyyyy_xxxxz_0[i] * fe_0 - 5.0 * ta_yyyyy_xxxxz_1[i] * fe_0 + ta_yyyyy_xxxxxz_0[i] * pa_x[i] - ta_yyyyy_xxxxxz_1[i] * pc_x[i];

        ta_xyyyyy_xxxxyy_0[i] = 4.0 * ta_yyyyy_xxxyy_0[i] * fe_0 - 4.0 * ta_yyyyy_xxxyy_1[i] * fe_0 + ta_yyyyy_xxxxyy_0[i] * pa_x[i] - ta_yyyyy_xxxxyy_1[i] * pc_x[i];

        ta_xyyyyy_xxxxyz_0[i] = 4.0 * ta_yyyyy_xxxyz_0[i] * fe_0 - 4.0 * ta_yyyyy_xxxyz_1[i] * fe_0 + ta_yyyyy_xxxxyz_0[i] * pa_x[i] - ta_yyyyy_xxxxyz_1[i] * pc_x[i];

        ta_xyyyyy_xxxxzz_0[i] = 4.0 * ta_yyyyy_xxxzz_0[i] * fe_0 - 4.0 * ta_yyyyy_xxxzz_1[i] * fe_0 + ta_yyyyy_xxxxzz_0[i] * pa_x[i] - ta_yyyyy_xxxxzz_1[i] * pc_x[i];

        ta_xyyyyy_xxxyyy_0[i] = 3.0 * ta_yyyyy_xxyyy_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyyy_1[i] * fe_0 + ta_yyyyy_xxxyyy_0[i] * pa_x[i] - ta_yyyyy_xxxyyy_1[i] * pc_x[i];

        ta_xyyyyy_xxxyyz_0[i] = 3.0 * ta_yyyyy_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyyz_1[i] * fe_0 + ta_yyyyy_xxxyyz_0[i] * pa_x[i] - ta_yyyyy_xxxyyz_1[i] * pc_x[i];

        ta_xyyyyy_xxxyzz_0[i] = 3.0 * ta_yyyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyzz_1[i] * fe_0 + ta_yyyyy_xxxyzz_0[i] * pa_x[i] - ta_yyyyy_xxxyzz_1[i] * pc_x[i];

        ta_xyyyyy_xxxzzz_0[i] = 3.0 * ta_yyyyy_xxzzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxzzz_1[i] * fe_0 + ta_yyyyy_xxxzzz_0[i] * pa_x[i] - ta_yyyyy_xxxzzz_1[i] * pc_x[i];

        ta_xyyyyy_xxyyyy_0[i] = 2.0 * ta_yyyyy_xyyyy_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyyy_1[i] * fe_0 + ta_yyyyy_xxyyyy_0[i] * pa_x[i] - ta_yyyyy_xxyyyy_1[i] * pc_x[i];

        ta_xyyyyy_xxyyyz_0[i] = 2.0 * ta_yyyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyyz_1[i] * fe_0 + ta_yyyyy_xxyyyz_0[i] * pa_x[i] - ta_yyyyy_xxyyyz_1[i] * pc_x[i];

        ta_xyyyyy_xxyyzz_0[i] = 2.0 * ta_yyyyy_xyyzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyzz_1[i] * fe_0 + ta_yyyyy_xxyyzz_0[i] * pa_x[i] - ta_yyyyy_xxyyzz_1[i] * pc_x[i];

        ta_xyyyyy_xxyzzz_0[i] = 2.0 * ta_yyyyy_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyzzz_1[i] * fe_0 + ta_yyyyy_xxyzzz_0[i] * pa_x[i] - ta_yyyyy_xxyzzz_1[i] * pc_x[i];

        ta_xyyyyy_xxzzzz_0[i] = 2.0 * ta_yyyyy_xzzzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xzzzz_1[i] * fe_0 + ta_yyyyy_xxzzzz_0[i] * pa_x[i] - ta_yyyyy_xxzzzz_1[i] * pc_x[i];

        ta_xyyyyy_xyyyyy_0[i] = ta_yyyyy_yyyyy_0[i] * fe_0 - ta_yyyyy_yyyyy_1[i] * fe_0 + ta_yyyyy_xyyyyy_0[i] * pa_x[i] - ta_yyyyy_xyyyyy_1[i] * pc_x[i];

        ta_xyyyyy_xyyyyz_0[i] = ta_yyyyy_yyyyz_0[i] * fe_0 - ta_yyyyy_yyyyz_1[i] * fe_0 + ta_yyyyy_xyyyyz_0[i] * pa_x[i] - ta_yyyyy_xyyyyz_1[i] * pc_x[i];

        ta_xyyyyy_xyyyzz_0[i] = ta_yyyyy_yyyzz_0[i] * fe_0 - ta_yyyyy_yyyzz_1[i] * fe_0 + ta_yyyyy_xyyyzz_0[i] * pa_x[i] - ta_yyyyy_xyyyzz_1[i] * pc_x[i];

        ta_xyyyyy_xyyzzz_0[i] = ta_yyyyy_yyzzz_0[i] * fe_0 - ta_yyyyy_yyzzz_1[i] * fe_0 + ta_yyyyy_xyyzzz_0[i] * pa_x[i] - ta_yyyyy_xyyzzz_1[i] * pc_x[i];

        ta_xyyyyy_xyzzzz_0[i] = ta_yyyyy_yzzzz_0[i] * fe_0 - ta_yyyyy_yzzzz_1[i] * fe_0 + ta_yyyyy_xyzzzz_0[i] * pa_x[i] - ta_yyyyy_xyzzzz_1[i] * pc_x[i];

        ta_xyyyyy_xzzzzz_0[i] = ta_yyyyy_zzzzz_0[i] * fe_0 - ta_yyyyy_zzzzz_1[i] * fe_0 + ta_yyyyy_xzzzzz_0[i] * pa_x[i] - ta_yyyyy_xzzzzz_1[i] * pc_x[i];

        ta_xyyyyy_yyyyyy_0[i] = ta_yyyyy_yyyyyy_0[i] * pa_x[i] - ta_yyyyy_yyyyyy_1[i] * pc_x[i];

        ta_xyyyyy_yyyyyz_0[i] = ta_yyyyy_yyyyyz_0[i] * pa_x[i] - ta_yyyyy_yyyyyz_1[i] * pc_x[i];

        ta_xyyyyy_yyyyzz_0[i] = ta_yyyyy_yyyyzz_0[i] * pa_x[i] - ta_yyyyy_yyyyzz_1[i] * pc_x[i];

        ta_xyyyyy_yyyzzz_0[i] = ta_yyyyy_yyyzzz_0[i] * pa_x[i] - ta_yyyyy_yyyzzz_1[i] * pc_x[i];

        ta_xyyyyy_yyzzzz_0[i] = ta_yyyyy_yyzzzz_0[i] * pa_x[i] - ta_yyyyy_yyzzzz_1[i] * pc_x[i];

        ta_xyyyyy_yzzzzz_0[i] = ta_yyyyy_yzzzzz_0[i] * pa_x[i] - ta_yyyyy_yzzzzz_1[i] * pc_x[i];

        ta_xyyyyy_zzzzzz_0[i] = ta_yyyyy_zzzzzz_0[i] * pa_x[i] - ta_yyyyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 448-476 components of targeted buffer : II

    auto ta_xyyyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 448);

    auto ta_xyyyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 449);

    auto ta_xyyyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 450);

    auto ta_xyyyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 451);

    auto ta_xyyyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 452);

    auto ta_xyyyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 453);

    auto ta_xyyyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 454);

    auto ta_xyyyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 455);

    auto ta_xyyyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 456);

    auto ta_xyyyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 457);

    auto ta_xyyyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 458);

    auto ta_xyyyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 459);

    auto ta_xyyyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 460);

    auto ta_xyyyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 461);

    auto ta_xyyyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 462);

    auto ta_xyyyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 463);

    auto ta_xyyyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 464);

    auto ta_xyyyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 465);

    auto ta_xyyyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 466);

    auto ta_xyyyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 467);

    auto ta_xyyyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 468);

    auto ta_xyyyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 469);

    auto ta_xyyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 470);

    auto ta_xyyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 471);

    auto ta_xyyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 472);

    auto ta_xyyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 473);

    auto ta_xyyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 474);

    auto ta_xyyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 475);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xyyyy_xxxxxx_0, ta_xyyyy_xxxxxx_1, ta_xyyyy_xxxxxy_0, ta_xyyyy_xxxxxy_1, ta_xyyyy_xxxxyy_0, ta_xyyyy_xxxxyy_1, ta_xyyyy_xxxyyy_0, ta_xyyyy_xxxyyy_1, ta_xyyyy_xxyyyy_0, ta_xyyyy_xxyyyy_1, ta_xyyyy_xyyyyy_0, ta_xyyyy_xyyyyy_1, ta_xyyyyz_xxxxxx_0, ta_xyyyyz_xxxxxy_0, ta_xyyyyz_xxxxxz_0, ta_xyyyyz_xxxxyy_0, ta_xyyyyz_xxxxyz_0, ta_xyyyyz_xxxxzz_0, ta_xyyyyz_xxxyyy_0, ta_xyyyyz_xxxyyz_0, ta_xyyyyz_xxxyzz_0, ta_xyyyyz_xxxzzz_0, ta_xyyyyz_xxyyyy_0, ta_xyyyyz_xxyyyz_0, ta_xyyyyz_xxyyzz_0, ta_xyyyyz_xxyzzz_0, ta_xyyyyz_xxzzzz_0, ta_xyyyyz_xyyyyy_0, ta_xyyyyz_xyyyyz_0, ta_xyyyyz_xyyyzz_0, ta_xyyyyz_xyyzzz_0, ta_xyyyyz_xyzzzz_0, ta_xyyyyz_xzzzzz_0, ta_xyyyyz_yyyyyy_0, ta_xyyyyz_yyyyyz_0, ta_xyyyyz_yyyyzz_0, ta_xyyyyz_yyyzzz_0, ta_xyyyyz_yyzzzz_0, ta_xyyyyz_yzzzzz_0, ta_xyyyyz_zzzzzz_0, ta_yyyyz_xxxxxz_0, ta_yyyyz_xxxxxz_1, ta_yyyyz_xxxxyz_0, ta_yyyyz_xxxxyz_1, ta_yyyyz_xxxxz_0, ta_yyyyz_xxxxz_1, ta_yyyyz_xxxxzz_0, ta_yyyyz_xxxxzz_1, ta_yyyyz_xxxyyz_0, ta_yyyyz_xxxyyz_1, ta_yyyyz_xxxyz_0, ta_yyyyz_xxxyz_1, ta_yyyyz_xxxyzz_0, ta_yyyyz_xxxyzz_1, ta_yyyyz_xxxzz_0, ta_yyyyz_xxxzz_1, ta_yyyyz_xxxzzz_0, ta_yyyyz_xxxzzz_1, ta_yyyyz_xxyyyz_0, ta_yyyyz_xxyyyz_1, ta_yyyyz_xxyyz_0, ta_yyyyz_xxyyz_1, ta_yyyyz_xxyyzz_0, ta_yyyyz_xxyyzz_1, ta_yyyyz_xxyzz_0, ta_yyyyz_xxyzz_1, ta_yyyyz_xxyzzz_0, ta_yyyyz_xxyzzz_1, ta_yyyyz_xxzzz_0, ta_yyyyz_xxzzz_1, ta_yyyyz_xxzzzz_0, ta_yyyyz_xxzzzz_1, ta_yyyyz_xyyyyz_0, ta_yyyyz_xyyyyz_1, ta_yyyyz_xyyyz_0, ta_yyyyz_xyyyz_1, ta_yyyyz_xyyyzz_0, ta_yyyyz_xyyyzz_1, ta_yyyyz_xyyzz_0, ta_yyyyz_xyyzz_1, ta_yyyyz_xyyzzz_0, ta_yyyyz_xyyzzz_1, ta_yyyyz_xyzzz_0, ta_yyyyz_xyzzz_1, ta_yyyyz_xyzzzz_0, ta_yyyyz_xyzzzz_1, ta_yyyyz_xzzzz_0, ta_yyyyz_xzzzz_1, ta_yyyyz_xzzzzz_0, ta_yyyyz_xzzzzz_1, ta_yyyyz_yyyyyy_0, ta_yyyyz_yyyyyy_1, ta_yyyyz_yyyyyz_0, ta_yyyyz_yyyyyz_1, ta_yyyyz_yyyyz_0, ta_yyyyz_yyyyz_1, ta_yyyyz_yyyyzz_0, ta_yyyyz_yyyyzz_1, ta_yyyyz_yyyzz_0, ta_yyyyz_yyyzz_1, ta_yyyyz_yyyzzz_0, ta_yyyyz_yyyzzz_1, ta_yyyyz_yyzzz_0, ta_yyyyz_yyzzz_1, ta_yyyyz_yyzzzz_0, ta_yyyyz_yyzzzz_1, ta_yyyyz_yzzzz_0, ta_yyyyz_yzzzz_1, ta_yyyyz_yzzzzz_0, ta_yyyyz_yzzzzz_1, ta_yyyyz_zzzzz_0, ta_yyyyz_zzzzz_1, ta_yyyyz_zzzzzz_0, ta_yyyyz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyz_xxxxxx_0[i] = ta_xyyyy_xxxxxx_0[i] * pa_z[i] - ta_xyyyy_xxxxxx_1[i] * pc_z[i];

        ta_xyyyyz_xxxxxy_0[i] = ta_xyyyy_xxxxxy_0[i] * pa_z[i] - ta_xyyyy_xxxxxy_1[i] * pc_z[i];

        ta_xyyyyz_xxxxxz_0[i] = 5.0 * ta_yyyyz_xxxxz_0[i] * fe_0 - 5.0 * ta_yyyyz_xxxxz_1[i] * fe_0 + ta_yyyyz_xxxxxz_0[i] * pa_x[i] - ta_yyyyz_xxxxxz_1[i] * pc_x[i];

        ta_xyyyyz_xxxxyy_0[i] = ta_xyyyy_xxxxyy_0[i] * pa_z[i] - ta_xyyyy_xxxxyy_1[i] * pc_z[i];

        ta_xyyyyz_xxxxyz_0[i] = 4.0 * ta_yyyyz_xxxyz_0[i] * fe_0 - 4.0 * ta_yyyyz_xxxyz_1[i] * fe_0 + ta_yyyyz_xxxxyz_0[i] * pa_x[i] - ta_yyyyz_xxxxyz_1[i] * pc_x[i];

        ta_xyyyyz_xxxxzz_0[i] = 4.0 * ta_yyyyz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyyyz_xxxzz_1[i] * fe_0 + ta_yyyyz_xxxxzz_0[i] * pa_x[i] - ta_yyyyz_xxxxzz_1[i] * pc_x[i];

        ta_xyyyyz_xxxyyy_0[i] = ta_xyyyy_xxxyyy_0[i] * pa_z[i] - ta_xyyyy_xxxyyy_1[i] * pc_z[i];

        ta_xyyyyz_xxxyyz_0[i] = 3.0 * ta_yyyyz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyyz_xxyyz_1[i] * fe_0 + ta_yyyyz_xxxyyz_0[i] * pa_x[i] - ta_yyyyz_xxxyyz_1[i] * pc_x[i];

        ta_xyyyyz_xxxyzz_0[i] = 3.0 * ta_yyyyz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyyz_xxyzz_1[i] * fe_0 + ta_yyyyz_xxxyzz_0[i] * pa_x[i] - ta_yyyyz_xxxyzz_1[i] * pc_x[i];

        ta_xyyyyz_xxxzzz_0[i] = 3.0 * ta_yyyyz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyyyz_xxzzz_1[i] * fe_0 + ta_yyyyz_xxxzzz_0[i] * pa_x[i] - ta_yyyyz_xxxzzz_1[i] * pc_x[i];

        ta_xyyyyz_xxyyyy_0[i] = ta_xyyyy_xxyyyy_0[i] * pa_z[i] - ta_xyyyy_xxyyyy_1[i] * pc_z[i];

        ta_xyyyyz_xxyyyz_0[i] = 2.0 * ta_yyyyz_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyyz_xyyyz_1[i] * fe_0 + ta_yyyyz_xxyyyz_0[i] * pa_x[i] - ta_yyyyz_xxyyyz_1[i] * pc_x[i];

        ta_xyyyyz_xxyyzz_0[i] = 2.0 * ta_yyyyz_xyyzz_0[i] * fe_0 - 2.0 * ta_yyyyz_xyyzz_1[i] * fe_0 + ta_yyyyz_xxyyzz_0[i] * pa_x[i] - ta_yyyyz_xxyyzz_1[i] * pc_x[i];

        ta_xyyyyz_xxyzzz_0[i] = 2.0 * ta_yyyyz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyyz_xyzzz_1[i] * fe_0 + ta_yyyyz_xxyzzz_0[i] * pa_x[i] - ta_yyyyz_xxyzzz_1[i] * pc_x[i];

        ta_xyyyyz_xxzzzz_0[i] = 2.0 * ta_yyyyz_xzzzz_0[i] * fe_0 - 2.0 * ta_yyyyz_xzzzz_1[i] * fe_0 + ta_yyyyz_xxzzzz_0[i] * pa_x[i] - ta_yyyyz_xxzzzz_1[i] * pc_x[i];

        ta_xyyyyz_xyyyyy_0[i] = ta_xyyyy_xyyyyy_0[i] * pa_z[i] - ta_xyyyy_xyyyyy_1[i] * pc_z[i];

        ta_xyyyyz_xyyyyz_0[i] = ta_yyyyz_yyyyz_0[i] * fe_0 - ta_yyyyz_yyyyz_1[i] * fe_0 + ta_yyyyz_xyyyyz_0[i] * pa_x[i] - ta_yyyyz_xyyyyz_1[i] * pc_x[i];

        ta_xyyyyz_xyyyzz_0[i] = ta_yyyyz_yyyzz_0[i] * fe_0 - ta_yyyyz_yyyzz_1[i] * fe_0 + ta_yyyyz_xyyyzz_0[i] * pa_x[i] - ta_yyyyz_xyyyzz_1[i] * pc_x[i];

        ta_xyyyyz_xyyzzz_0[i] = ta_yyyyz_yyzzz_0[i] * fe_0 - ta_yyyyz_yyzzz_1[i] * fe_0 + ta_yyyyz_xyyzzz_0[i] * pa_x[i] - ta_yyyyz_xyyzzz_1[i] * pc_x[i];

        ta_xyyyyz_xyzzzz_0[i] = ta_yyyyz_yzzzz_0[i] * fe_0 - ta_yyyyz_yzzzz_1[i] * fe_0 + ta_yyyyz_xyzzzz_0[i] * pa_x[i] - ta_yyyyz_xyzzzz_1[i] * pc_x[i];

        ta_xyyyyz_xzzzzz_0[i] = ta_yyyyz_zzzzz_0[i] * fe_0 - ta_yyyyz_zzzzz_1[i] * fe_0 + ta_yyyyz_xzzzzz_0[i] * pa_x[i] - ta_yyyyz_xzzzzz_1[i] * pc_x[i];

        ta_xyyyyz_yyyyyy_0[i] = ta_yyyyz_yyyyyy_0[i] * pa_x[i] - ta_yyyyz_yyyyyy_1[i] * pc_x[i];

        ta_xyyyyz_yyyyyz_0[i] = ta_yyyyz_yyyyyz_0[i] * pa_x[i] - ta_yyyyz_yyyyyz_1[i] * pc_x[i];

        ta_xyyyyz_yyyyzz_0[i] = ta_yyyyz_yyyyzz_0[i] * pa_x[i] - ta_yyyyz_yyyyzz_1[i] * pc_x[i];

        ta_xyyyyz_yyyzzz_0[i] = ta_yyyyz_yyyzzz_0[i] * pa_x[i] - ta_yyyyz_yyyzzz_1[i] * pc_x[i];

        ta_xyyyyz_yyzzzz_0[i] = ta_yyyyz_yyzzzz_0[i] * pa_x[i] - ta_yyyyz_yyzzzz_1[i] * pc_x[i];

        ta_xyyyyz_yzzzzz_0[i] = ta_yyyyz_yzzzzz_0[i] * pa_x[i] - ta_yyyyz_yzzzzz_1[i] * pc_x[i];

        ta_xyyyyz_zzzzzz_0[i] = ta_yyyyz_zzzzzz_0[i] * pa_x[i] - ta_yyyyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 476-504 components of targeted buffer : II

    auto ta_xyyyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 476);

    auto ta_xyyyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 477);

    auto ta_xyyyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 478);

    auto ta_xyyyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 479);

    auto ta_xyyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 480);

    auto ta_xyyyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 481);

    auto ta_xyyyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 482);

    auto ta_xyyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 483);

    auto ta_xyyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 484);

    auto ta_xyyyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 485);

    auto ta_xyyyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 486);

    auto ta_xyyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 487);

    auto ta_xyyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 488);

    auto ta_xyyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 489);

    auto ta_xyyyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 490);

    auto ta_xyyyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 491);

    auto ta_xyyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 492);

    auto ta_xyyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 493);

    auto ta_xyyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 494);

    auto ta_xyyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 495);

    auto ta_xyyyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 496);

    auto ta_xyyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 497);

    auto ta_xyyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 498);

    auto ta_xyyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 499);

    auto ta_xyyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 500);

    auto ta_xyyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 501);

    auto ta_xyyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 502);

    auto ta_xyyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 503);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyyyzz_xxxxxx_0, ta_xyyyzz_xxxxxy_0, ta_xyyyzz_xxxxxz_0, ta_xyyyzz_xxxxyy_0, ta_xyyyzz_xxxxyz_0, ta_xyyyzz_xxxxzz_0, ta_xyyyzz_xxxyyy_0, ta_xyyyzz_xxxyyz_0, ta_xyyyzz_xxxyzz_0, ta_xyyyzz_xxxzzz_0, ta_xyyyzz_xxyyyy_0, ta_xyyyzz_xxyyyz_0, ta_xyyyzz_xxyyzz_0, ta_xyyyzz_xxyzzz_0, ta_xyyyzz_xxzzzz_0, ta_xyyyzz_xyyyyy_0, ta_xyyyzz_xyyyyz_0, ta_xyyyzz_xyyyzz_0, ta_xyyyzz_xyyzzz_0, ta_xyyyzz_xyzzzz_0, ta_xyyyzz_xzzzzz_0, ta_xyyyzz_yyyyyy_0, ta_xyyyzz_yyyyyz_0, ta_xyyyzz_yyyyzz_0, ta_xyyyzz_yyyzzz_0, ta_xyyyzz_yyzzzz_0, ta_xyyyzz_yzzzzz_0, ta_xyyyzz_zzzzzz_0, ta_yyyzz_xxxxx_0, ta_yyyzz_xxxxx_1, ta_yyyzz_xxxxxx_0, ta_yyyzz_xxxxxx_1, ta_yyyzz_xxxxxy_0, ta_yyyzz_xxxxxy_1, ta_yyyzz_xxxxxz_0, ta_yyyzz_xxxxxz_1, ta_yyyzz_xxxxy_0, ta_yyyzz_xxxxy_1, ta_yyyzz_xxxxyy_0, ta_yyyzz_xxxxyy_1, ta_yyyzz_xxxxyz_0, ta_yyyzz_xxxxyz_1, ta_yyyzz_xxxxz_0, ta_yyyzz_xxxxz_1, ta_yyyzz_xxxxzz_0, ta_yyyzz_xxxxzz_1, ta_yyyzz_xxxyy_0, ta_yyyzz_xxxyy_1, ta_yyyzz_xxxyyy_0, ta_yyyzz_xxxyyy_1, ta_yyyzz_xxxyyz_0, ta_yyyzz_xxxyyz_1, ta_yyyzz_xxxyz_0, ta_yyyzz_xxxyz_1, ta_yyyzz_xxxyzz_0, ta_yyyzz_xxxyzz_1, ta_yyyzz_xxxzz_0, ta_yyyzz_xxxzz_1, ta_yyyzz_xxxzzz_0, ta_yyyzz_xxxzzz_1, ta_yyyzz_xxyyy_0, ta_yyyzz_xxyyy_1, ta_yyyzz_xxyyyy_0, ta_yyyzz_xxyyyy_1, ta_yyyzz_xxyyyz_0, ta_yyyzz_xxyyyz_1, ta_yyyzz_xxyyz_0, ta_yyyzz_xxyyz_1, ta_yyyzz_xxyyzz_0, ta_yyyzz_xxyyzz_1, ta_yyyzz_xxyzz_0, ta_yyyzz_xxyzz_1, ta_yyyzz_xxyzzz_0, ta_yyyzz_xxyzzz_1, ta_yyyzz_xxzzz_0, ta_yyyzz_xxzzz_1, ta_yyyzz_xxzzzz_0, ta_yyyzz_xxzzzz_1, ta_yyyzz_xyyyy_0, ta_yyyzz_xyyyy_1, ta_yyyzz_xyyyyy_0, ta_yyyzz_xyyyyy_1, ta_yyyzz_xyyyyz_0, ta_yyyzz_xyyyyz_1, ta_yyyzz_xyyyz_0, ta_yyyzz_xyyyz_1, ta_yyyzz_xyyyzz_0, ta_yyyzz_xyyyzz_1, ta_yyyzz_xyyzz_0, ta_yyyzz_xyyzz_1, ta_yyyzz_xyyzzz_0, ta_yyyzz_xyyzzz_1, ta_yyyzz_xyzzz_0, ta_yyyzz_xyzzz_1, ta_yyyzz_xyzzzz_0, ta_yyyzz_xyzzzz_1, ta_yyyzz_xzzzz_0, ta_yyyzz_xzzzz_1, ta_yyyzz_xzzzzz_0, ta_yyyzz_xzzzzz_1, ta_yyyzz_yyyyy_0, ta_yyyzz_yyyyy_1, ta_yyyzz_yyyyyy_0, ta_yyyzz_yyyyyy_1, ta_yyyzz_yyyyyz_0, ta_yyyzz_yyyyyz_1, ta_yyyzz_yyyyz_0, ta_yyyzz_yyyyz_1, ta_yyyzz_yyyyzz_0, ta_yyyzz_yyyyzz_1, ta_yyyzz_yyyzz_0, ta_yyyzz_yyyzz_1, ta_yyyzz_yyyzzz_0, ta_yyyzz_yyyzzz_1, ta_yyyzz_yyzzz_0, ta_yyyzz_yyzzz_1, ta_yyyzz_yyzzzz_0, ta_yyyzz_yyzzzz_1, ta_yyyzz_yzzzz_0, ta_yyyzz_yzzzz_1, ta_yyyzz_yzzzzz_0, ta_yyyzz_yzzzzz_1, ta_yyyzz_zzzzz_0, ta_yyyzz_zzzzz_1, ta_yyyzz_zzzzzz_0, ta_yyyzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyzz_xxxxxx_0[i] = 6.0 * ta_yyyzz_xxxxx_0[i] * fe_0 - 6.0 * ta_yyyzz_xxxxx_1[i] * fe_0 + ta_yyyzz_xxxxxx_0[i] * pa_x[i] - ta_yyyzz_xxxxxx_1[i] * pc_x[i];

        ta_xyyyzz_xxxxxy_0[i] = 5.0 * ta_yyyzz_xxxxy_0[i] * fe_0 - 5.0 * ta_yyyzz_xxxxy_1[i] * fe_0 + ta_yyyzz_xxxxxy_0[i] * pa_x[i] - ta_yyyzz_xxxxxy_1[i] * pc_x[i];

        ta_xyyyzz_xxxxxz_0[i] = 5.0 * ta_yyyzz_xxxxz_0[i] * fe_0 - 5.0 * ta_yyyzz_xxxxz_1[i] * fe_0 + ta_yyyzz_xxxxxz_0[i] * pa_x[i] - ta_yyyzz_xxxxxz_1[i] * pc_x[i];

        ta_xyyyzz_xxxxyy_0[i] = 4.0 * ta_yyyzz_xxxyy_0[i] * fe_0 - 4.0 * ta_yyyzz_xxxyy_1[i] * fe_0 + ta_yyyzz_xxxxyy_0[i] * pa_x[i] - ta_yyyzz_xxxxyy_1[i] * pc_x[i];

        ta_xyyyzz_xxxxyz_0[i] = 4.0 * ta_yyyzz_xxxyz_0[i] * fe_0 - 4.0 * ta_yyyzz_xxxyz_1[i] * fe_0 + ta_yyyzz_xxxxyz_0[i] * pa_x[i] - ta_yyyzz_xxxxyz_1[i] * pc_x[i];

        ta_xyyyzz_xxxxzz_0[i] = 4.0 * ta_yyyzz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyyzz_xxxzz_1[i] * fe_0 + ta_yyyzz_xxxxzz_0[i] * pa_x[i] - ta_yyyzz_xxxxzz_1[i] * pc_x[i];

        ta_xyyyzz_xxxyyy_0[i] = 3.0 * ta_yyyzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yyyzz_xxyyy_1[i] * fe_0 + ta_yyyzz_xxxyyy_0[i] * pa_x[i] - ta_yyyzz_xxxyyy_1[i] * pc_x[i];

        ta_xyyyzz_xxxyyz_0[i] = 3.0 * ta_yyyzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxyyz_1[i] * fe_0 + ta_yyyzz_xxxyyz_0[i] * pa_x[i] - ta_yyyzz_xxxyyz_1[i] * pc_x[i];

        ta_xyyyzz_xxxyzz_0[i] = 3.0 * ta_yyyzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxyzz_1[i] * fe_0 + ta_yyyzz_xxxyzz_0[i] * pa_x[i] - ta_yyyzz_xxxyzz_1[i] * pc_x[i];

        ta_xyyyzz_xxxzzz_0[i] = 3.0 * ta_yyyzz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxzzz_1[i] * fe_0 + ta_yyyzz_xxxzzz_0[i] * pa_x[i] - ta_yyyzz_xxxzzz_1[i] * pc_x[i];

        ta_xyyyzz_xxyyyy_0[i] = 2.0 * ta_yyyzz_xyyyy_0[i] * fe_0 - 2.0 * ta_yyyzz_xyyyy_1[i] * fe_0 + ta_yyyzz_xxyyyy_0[i] * pa_x[i] - ta_yyyzz_xxyyyy_1[i] * pc_x[i];

        ta_xyyyzz_xxyyyz_0[i] = 2.0 * ta_yyyzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyyyz_1[i] * fe_0 + ta_yyyzz_xxyyyz_0[i] * pa_x[i] - ta_yyyzz_xxyyyz_1[i] * pc_x[i];

        ta_xyyyzz_xxyyzz_0[i] = 2.0 * ta_yyyzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyyzz_1[i] * fe_0 + ta_yyyzz_xxyyzz_0[i] * pa_x[i] - ta_yyyzz_xxyyzz_1[i] * pc_x[i];

        ta_xyyyzz_xxyzzz_0[i] = 2.0 * ta_yyyzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyzzz_1[i] * fe_0 + ta_yyyzz_xxyzzz_0[i] * pa_x[i] - ta_yyyzz_xxyzzz_1[i] * pc_x[i];

        ta_xyyyzz_xxzzzz_0[i] = 2.0 * ta_yyyzz_xzzzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xzzzz_1[i] * fe_0 + ta_yyyzz_xxzzzz_0[i] * pa_x[i] - ta_yyyzz_xxzzzz_1[i] * pc_x[i];

        ta_xyyyzz_xyyyyy_0[i] = ta_yyyzz_yyyyy_0[i] * fe_0 - ta_yyyzz_yyyyy_1[i] * fe_0 + ta_yyyzz_xyyyyy_0[i] * pa_x[i] - ta_yyyzz_xyyyyy_1[i] * pc_x[i];

        ta_xyyyzz_xyyyyz_0[i] = ta_yyyzz_yyyyz_0[i] * fe_0 - ta_yyyzz_yyyyz_1[i] * fe_0 + ta_yyyzz_xyyyyz_0[i] * pa_x[i] - ta_yyyzz_xyyyyz_1[i] * pc_x[i];

        ta_xyyyzz_xyyyzz_0[i] = ta_yyyzz_yyyzz_0[i] * fe_0 - ta_yyyzz_yyyzz_1[i] * fe_0 + ta_yyyzz_xyyyzz_0[i] * pa_x[i] - ta_yyyzz_xyyyzz_1[i] * pc_x[i];

        ta_xyyyzz_xyyzzz_0[i] = ta_yyyzz_yyzzz_0[i] * fe_0 - ta_yyyzz_yyzzz_1[i] * fe_0 + ta_yyyzz_xyyzzz_0[i] * pa_x[i] - ta_yyyzz_xyyzzz_1[i] * pc_x[i];

        ta_xyyyzz_xyzzzz_0[i] = ta_yyyzz_yzzzz_0[i] * fe_0 - ta_yyyzz_yzzzz_1[i] * fe_0 + ta_yyyzz_xyzzzz_0[i] * pa_x[i] - ta_yyyzz_xyzzzz_1[i] * pc_x[i];

        ta_xyyyzz_xzzzzz_0[i] = ta_yyyzz_zzzzz_0[i] * fe_0 - ta_yyyzz_zzzzz_1[i] * fe_0 + ta_yyyzz_xzzzzz_0[i] * pa_x[i] - ta_yyyzz_xzzzzz_1[i] * pc_x[i];

        ta_xyyyzz_yyyyyy_0[i] = ta_yyyzz_yyyyyy_0[i] * pa_x[i] - ta_yyyzz_yyyyyy_1[i] * pc_x[i];

        ta_xyyyzz_yyyyyz_0[i] = ta_yyyzz_yyyyyz_0[i] * pa_x[i] - ta_yyyzz_yyyyyz_1[i] * pc_x[i];

        ta_xyyyzz_yyyyzz_0[i] = ta_yyyzz_yyyyzz_0[i] * pa_x[i] - ta_yyyzz_yyyyzz_1[i] * pc_x[i];

        ta_xyyyzz_yyyzzz_0[i] = ta_yyyzz_yyyzzz_0[i] * pa_x[i] - ta_yyyzz_yyyzzz_1[i] * pc_x[i];

        ta_xyyyzz_yyzzzz_0[i] = ta_yyyzz_yyzzzz_0[i] * pa_x[i] - ta_yyyzz_yyzzzz_1[i] * pc_x[i];

        ta_xyyyzz_yzzzzz_0[i] = ta_yyyzz_yzzzzz_0[i] * pa_x[i] - ta_yyyzz_yzzzzz_1[i] * pc_x[i];

        ta_xyyyzz_zzzzzz_0[i] = ta_yyyzz_zzzzzz_0[i] * pa_x[i] - ta_yyyzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 504-532 components of targeted buffer : II

    auto ta_xyyzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 504);

    auto ta_xyyzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 505);

    auto ta_xyyzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 506);

    auto ta_xyyzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 507);

    auto ta_xyyzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 508);

    auto ta_xyyzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 509);

    auto ta_xyyzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 510);

    auto ta_xyyzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 511);

    auto ta_xyyzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 512);

    auto ta_xyyzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 513);

    auto ta_xyyzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 514);

    auto ta_xyyzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 515);

    auto ta_xyyzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 516);

    auto ta_xyyzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 517);

    auto ta_xyyzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 518);

    auto ta_xyyzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 519);

    auto ta_xyyzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 520);

    auto ta_xyyzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 521);

    auto ta_xyyzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 522);

    auto ta_xyyzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 523);

    auto ta_xyyzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 524);

    auto ta_xyyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 525);

    auto ta_xyyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 526);

    auto ta_xyyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 527);

    auto ta_xyyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 528);

    auto ta_xyyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 529);

    auto ta_xyyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 530);

    auto ta_xyyzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 531);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyyzzz_xxxxxx_0, ta_xyyzzz_xxxxxy_0, ta_xyyzzz_xxxxxz_0, ta_xyyzzz_xxxxyy_0, ta_xyyzzz_xxxxyz_0, ta_xyyzzz_xxxxzz_0, ta_xyyzzz_xxxyyy_0, ta_xyyzzz_xxxyyz_0, ta_xyyzzz_xxxyzz_0, ta_xyyzzz_xxxzzz_0, ta_xyyzzz_xxyyyy_0, ta_xyyzzz_xxyyyz_0, ta_xyyzzz_xxyyzz_0, ta_xyyzzz_xxyzzz_0, ta_xyyzzz_xxzzzz_0, ta_xyyzzz_xyyyyy_0, ta_xyyzzz_xyyyyz_0, ta_xyyzzz_xyyyzz_0, ta_xyyzzz_xyyzzz_0, ta_xyyzzz_xyzzzz_0, ta_xyyzzz_xzzzzz_0, ta_xyyzzz_yyyyyy_0, ta_xyyzzz_yyyyyz_0, ta_xyyzzz_yyyyzz_0, ta_xyyzzz_yyyzzz_0, ta_xyyzzz_yyzzzz_0, ta_xyyzzz_yzzzzz_0, ta_xyyzzz_zzzzzz_0, ta_yyzzz_xxxxx_0, ta_yyzzz_xxxxx_1, ta_yyzzz_xxxxxx_0, ta_yyzzz_xxxxxx_1, ta_yyzzz_xxxxxy_0, ta_yyzzz_xxxxxy_1, ta_yyzzz_xxxxxz_0, ta_yyzzz_xxxxxz_1, ta_yyzzz_xxxxy_0, ta_yyzzz_xxxxy_1, ta_yyzzz_xxxxyy_0, ta_yyzzz_xxxxyy_1, ta_yyzzz_xxxxyz_0, ta_yyzzz_xxxxyz_1, ta_yyzzz_xxxxz_0, ta_yyzzz_xxxxz_1, ta_yyzzz_xxxxzz_0, ta_yyzzz_xxxxzz_1, ta_yyzzz_xxxyy_0, ta_yyzzz_xxxyy_1, ta_yyzzz_xxxyyy_0, ta_yyzzz_xxxyyy_1, ta_yyzzz_xxxyyz_0, ta_yyzzz_xxxyyz_1, ta_yyzzz_xxxyz_0, ta_yyzzz_xxxyz_1, ta_yyzzz_xxxyzz_0, ta_yyzzz_xxxyzz_1, ta_yyzzz_xxxzz_0, ta_yyzzz_xxxzz_1, ta_yyzzz_xxxzzz_0, ta_yyzzz_xxxzzz_1, ta_yyzzz_xxyyy_0, ta_yyzzz_xxyyy_1, ta_yyzzz_xxyyyy_0, ta_yyzzz_xxyyyy_1, ta_yyzzz_xxyyyz_0, ta_yyzzz_xxyyyz_1, ta_yyzzz_xxyyz_0, ta_yyzzz_xxyyz_1, ta_yyzzz_xxyyzz_0, ta_yyzzz_xxyyzz_1, ta_yyzzz_xxyzz_0, ta_yyzzz_xxyzz_1, ta_yyzzz_xxyzzz_0, ta_yyzzz_xxyzzz_1, ta_yyzzz_xxzzz_0, ta_yyzzz_xxzzz_1, ta_yyzzz_xxzzzz_0, ta_yyzzz_xxzzzz_1, ta_yyzzz_xyyyy_0, ta_yyzzz_xyyyy_1, ta_yyzzz_xyyyyy_0, ta_yyzzz_xyyyyy_1, ta_yyzzz_xyyyyz_0, ta_yyzzz_xyyyyz_1, ta_yyzzz_xyyyz_0, ta_yyzzz_xyyyz_1, ta_yyzzz_xyyyzz_0, ta_yyzzz_xyyyzz_1, ta_yyzzz_xyyzz_0, ta_yyzzz_xyyzz_1, ta_yyzzz_xyyzzz_0, ta_yyzzz_xyyzzz_1, ta_yyzzz_xyzzz_0, ta_yyzzz_xyzzz_1, ta_yyzzz_xyzzzz_0, ta_yyzzz_xyzzzz_1, ta_yyzzz_xzzzz_0, ta_yyzzz_xzzzz_1, ta_yyzzz_xzzzzz_0, ta_yyzzz_xzzzzz_1, ta_yyzzz_yyyyy_0, ta_yyzzz_yyyyy_1, ta_yyzzz_yyyyyy_0, ta_yyzzz_yyyyyy_1, ta_yyzzz_yyyyyz_0, ta_yyzzz_yyyyyz_1, ta_yyzzz_yyyyz_0, ta_yyzzz_yyyyz_1, ta_yyzzz_yyyyzz_0, ta_yyzzz_yyyyzz_1, ta_yyzzz_yyyzz_0, ta_yyzzz_yyyzz_1, ta_yyzzz_yyyzzz_0, ta_yyzzz_yyyzzz_1, ta_yyzzz_yyzzz_0, ta_yyzzz_yyzzz_1, ta_yyzzz_yyzzzz_0, ta_yyzzz_yyzzzz_1, ta_yyzzz_yzzzz_0, ta_yyzzz_yzzzz_1, ta_yyzzz_yzzzzz_0, ta_yyzzz_yzzzzz_1, ta_yyzzz_zzzzz_0, ta_yyzzz_zzzzz_1, ta_yyzzz_zzzzzz_0, ta_yyzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzzz_xxxxxx_0[i] = 6.0 * ta_yyzzz_xxxxx_0[i] * fe_0 - 6.0 * ta_yyzzz_xxxxx_1[i] * fe_0 + ta_yyzzz_xxxxxx_0[i] * pa_x[i] - ta_yyzzz_xxxxxx_1[i] * pc_x[i];

        ta_xyyzzz_xxxxxy_0[i] = 5.0 * ta_yyzzz_xxxxy_0[i] * fe_0 - 5.0 * ta_yyzzz_xxxxy_1[i] * fe_0 + ta_yyzzz_xxxxxy_0[i] * pa_x[i] - ta_yyzzz_xxxxxy_1[i] * pc_x[i];

        ta_xyyzzz_xxxxxz_0[i] = 5.0 * ta_yyzzz_xxxxz_0[i] * fe_0 - 5.0 * ta_yyzzz_xxxxz_1[i] * fe_0 + ta_yyzzz_xxxxxz_0[i] * pa_x[i] - ta_yyzzz_xxxxxz_1[i] * pc_x[i];

        ta_xyyzzz_xxxxyy_0[i] = 4.0 * ta_yyzzz_xxxyy_0[i] * fe_0 - 4.0 * ta_yyzzz_xxxyy_1[i] * fe_0 + ta_yyzzz_xxxxyy_0[i] * pa_x[i] - ta_yyzzz_xxxxyy_1[i] * pc_x[i];

        ta_xyyzzz_xxxxyz_0[i] = 4.0 * ta_yyzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_yyzzz_xxxyz_1[i] * fe_0 + ta_yyzzz_xxxxyz_0[i] * pa_x[i] - ta_yyzzz_xxxxyz_1[i] * pc_x[i];

        ta_xyyzzz_xxxxzz_0[i] = 4.0 * ta_yyzzz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyzzz_xxxzz_1[i] * fe_0 + ta_yyzzz_xxxxzz_0[i] * pa_x[i] - ta_yyzzz_xxxxzz_1[i] * pc_x[i];

        ta_xyyzzz_xxxyyy_0[i] = 3.0 * ta_yyzzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yyzzz_xxyyy_1[i] * fe_0 + ta_yyzzz_xxxyyy_0[i] * pa_x[i] - ta_yyzzz_xxxyyy_1[i] * pc_x[i];

        ta_xyyzzz_xxxyyz_0[i] = 3.0 * ta_yyzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxyyz_1[i] * fe_0 + ta_yyzzz_xxxyyz_0[i] * pa_x[i] - ta_yyzzz_xxxyyz_1[i] * pc_x[i];

        ta_xyyzzz_xxxyzz_0[i] = 3.0 * ta_yyzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxyzz_1[i] * fe_0 + ta_yyzzz_xxxyzz_0[i] * pa_x[i] - ta_yyzzz_xxxyzz_1[i] * pc_x[i];

        ta_xyyzzz_xxxzzz_0[i] = 3.0 * ta_yyzzz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxzzz_1[i] * fe_0 + ta_yyzzz_xxxzzz_0[i] * pa_x[i] - ta_yyzzz_xxxzzz_1[i] * pc_x[i];

        ta_xyyzzz_xxyyyy_0[i] = 2.0 * ta_yyzzz_xyyyy_0[i] * fe_0 - 2.0 * ta_yyzzz_xyyyy_1[i] * fe_0 + ta_yyzzz_xxyyyy_0[i] * pa_x[i] - ta_yyzzz_xxyyyy_1[i] * pc_x[i];

        ta_xyyzzz_xxyyyz_0[i] = 2.0 * ta_yyzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyyyz_1[i] * fe_0 + ta_yyzzz_xxyyyz_0[i] * pa_x[i] - ta_yyzzz_xxyyyz_1[i] * pc_x[i];

        ta_xyyzzz_xxyyzz_0[i] = 2.0 * ta_yyzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyyzz_1[i] * fe_0 + ta_yyzzz_xxyyzz_0[i] * pa_x[i] - ta_yyzzz_xxyyzz_1[i] * pc_x[i];

        ta_xyyzzz_xxyzzz_0[i] = 2.0 * ta_yyzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyzzz_1[i] * fe_0 + ta_yyzzz_xxyzzz_0[i] * pa_x[i] - ta_yyzzz_xxyzzz_1[i] * pc_x[i];

        ta_xyyzzz_xxzzzz_0[i] = 2.0 * ta_yyzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xzzzz_1[i] * fe_0 + ta_yyzzz_xxzzzz_0[i] * pa_x[i] - ta_yyzzz_xxzzzz_1[i] * pc_x[i];

        ta_xyyzzz_xyyyyy_0[i] = ta_yyzzz_yyyyy_0[i] * fe_0 - ta_yyzzz_yyyyy_1[i] * fe_0 + ta_yyzzz_xyyyyy_0[i] * pa_x[i] - ta_yyzzz_xyyyyy_1[i] * pc_x[i];

        ta_xyyzzz_xyyyyz_0[i] = ta_yyzzz_yyyyz_0[i] * fe_0 - ta_yyzzz_yyyyz_1[i] * fe_0 + ta_yyzzz_xyyyyz_0[i] * pa_x[i] - ta_yyzzz_xyyyyz_1[i] * pc_x[i];

        ta_xyyzzz_xyyyzz_0[i] = ta_yyzzz_yyyzz_0[i] * fe_0 - ta_yyzzz_yyyzz_1[i] * fe_0 + ta_yyzzz_xyyyzz_0[i] * pa_x[i] - ta_yyzzz_xyyyzz_1[i] * pc_x[i];

        ta_xyyzzz_xyyzzz_0[i] = ta_yyzzz_yyzzz_0[i] * fe_0 - ta_yyzzz_yyzzz_1[i] * fe_0 + ta_yyzzz_xyyzzz_0[i] * pa_x[i] - ta_yyzzz_xyyzzz_1[i] * pc_x[i];

        ta_xyyzzz_xyzzzz_0[i] = ta_yyzzz_yzzzz_0[i] * fe_0 - ta_yyzzz_yzzzz_1[i] * fe_0 + ta_yyzzz_xyzzzz_0[i] * pa_x[i] - ta_yyzzz_xyzzzz_1[i] * pc_x[i];

        ta_xyyzzz_xzzzzz_0[i] = ta_yyzzz_zzzzz_0[i] * fe_0 - ta_yyzzz_zzzzz_1[i] * fe_0 + ta_yyzzz_xzzzzz_0[i] * pa_x[i] - ta_yyzzz_xzzzzz_1[i] * pc_x[i];

        ta_xyyzzz_yyyyyy_0[i] = ta_yyzzz_yyyyyy_0[i] * pa_x[i] - ta_yyzzz_yyyyyy_1[i] * pc_x[i];

        ta_xyyzzz_yyyyyz_0[i] = ta_yyzzz_yyyyyz_0[i] * pa_x[i] - ta_yyzzz_yyyyyz_1[i] * pc_x[i];

        ta_xyyzzz_yyyyzz_0[i] = ta_yyzzz_yyyyzz_0[i] * pa_x[i] - ta_yyzzz_yyyyzz_1[i] * pc_x[i];

        ta_xyyzzz_yyyzzz_0[i] = ta_yyzzz_yyyzzz_0[i] * pa_x[i] - ta_yyzzz_yyyzzz_1[i] * pc_x[i];

        ta_xyyzzz_yyzzzz_0[i] = ta_yyzzz_yyzzzz_0[i] * pa_x[i] - ta_yyzzz_yyzzzz_1[i] * pc_x[i];

        ta_xyyzzz_yzzzzz_0[i] = ta_yyzzz_yzzzzz_0[i] * pa_x[i] - ta_yyzzz_yzzzzz_1[i] * pc_x[i];

        ta_xyyzzz_zzzzzz_0[i] = ta_yyzzz_zzzzzz_0[i] * pa_x[i] - ta_yyzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 532-560 components of targeted buffer : II

    auto ta_xyzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 532);

    auto ta_xyzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 533);

    auto ta_xyzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 534);

    auto ta_xyzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 535);

    auto ta_xyzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 536);

    auto ta_xyzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 537);

    auto ta_xyzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 538);

    auto ta_xyzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 539);

    auto ta_xyzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 540);

    auto ta_xyzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 541);

    auto ta_xyzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 542);

    auto ta_xyzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 543);

    auto ta_xyzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 544);

    auto ta_xyzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 545);

    auto ta_xyzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 546);

    auto ta_xyzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 547);

    auto ta_xyzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 548);

    auto ta_xyzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 549);

    auto ta_xyzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 550);

    auto ta_xyzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 551);

    auto ta_xyzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 552);

    auto ta_xyzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 553);

    auto ta_xyzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 554);

    auto ta_xyzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 555);

    auto ta_xyzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 556);

    auto ta_xyzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 557);

    auto ta_xyzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 558);

    auto ta_xyzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 559);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xyzzzz_xxxxxx_0, ta_xyzzzz_xxxxxy_0, ta_xyzzzz_xxxxxz_0, ta_xyzzzz_xxxxyy_0, ta_xyzzzz_xxxxyz_0, ta_xyzzzz_xxxxzz_0, ta_xyzzzz_xxxyyy_0, ta_xyzzzz_xxxyyz_0, ta_xyzzzz_xxxyzz_0, ta_xyzzzz_xxxzzz_0, ta_xyzzzz_xxyyyy_0, ta_xyzzzz_xxyyyz_0, ta_xyzzzz_xxyyzz_0, ta_xyzzzz_xxyzzz_0, ta_xyzzzz_xxzzzz_0, ta_xyzzzz_xyyyyy_0, ta_xyzzzz_xyyyyz_0, ta_xyzzzz_xyyyzz_0, ta_xyzzzz_xyyzzz_0, ta_xyzzzz_xyzzzz_0, ta_xyzzzz_xzzzzz_0, ta_xyzzzz_yyyyyy_0, ta_xyzzzz_yyyyyz_0, ta_xyzzzz_yyyyzz_0, ta_xyzzzz_yyyzzz_0, ta_xyzzzz_yyzzzz_0, ta_xyzzzz_yzzzzz_0, ta_xyzzzz_zzzzzz_0, ta_xzzzz_xxxxxx_0, ta_xzzzz_xxxxxx_1, ta_xzzzz_xxxxxz_0, ta_xzzzz_xxxxxz_1, ta_xzzzz_xxxxzz_0, ta_xzzzz_xxxxzz_1, ta_xzzzz_xxxzzz_0, ta_xzzzz_xxxzzz_1, ta_xzzzz_xxzzzz_0, ta_xzzzz_xxzzzz_1, ta_xzzzz_xzzzzz_0, ta_xzzzz_xzzzzz_1, ta_yzzzz_xxxxxy_0, ta_yzzzz_xxxxxy_1, ta_yzzzz_xxxxy_0, ta_yzzzz_xxxxy_1, ta_yzzzz_xxxxyy_0, ta_yzzzz_xxxxyy_1, ta_yzzzz_xxxxyz_0, ta_yzzzz_xxxxyz_1, ta_yzzzz_xxxyy_0, ta_yzzzz_xxxyy_1, ta_yzzzz_xxxyyy_0, ta_yzzzz_xxxyyy_1, ta_yzzzz_xxxyyz_0, ta_yzzzz_xxxyyz_1, ta_yzzzz_xxxyz_0, ta_yzzzz_xxxyz_1, ta_yzzzz_xxxyzz_0, ta_yzzzz_xxxyzz_1, ta_yzzzz_xxyyy_0, ta_yzzzz_xxyyy_1, ta_yzzzz_xxyyyy_0, ta_yzzzz_xxyyyy_1, ta_yzzzz_xxyyyz_0, ta_yzzzz_xxyyyz_1, ta_yzzzz_xxyyz_0, ta_yzzzz_xxyyz_1, ta_yzzzz_xxyyzz_0, ta_yzzzz_xxyyzz_1, ta_yzzzz_xxyzz_0, ta_yzzzz_xxyzz_1, ta_yzzzz_xxyzzz_0, ta_yzzzz_xxyzzz_1, ta_yzzzz_xyyyy_0, ta_yzzzz_xyyyy_1, ta_yzzzz_xyyyyy_0, ta_yzzzz_xyyyyy_1, ta_yzzzz_xyyyyz_0, ta_yzzzz_xyyyyz_1, ta_yzzzz_xyyyz_0, ta_yzzzz_xyyyz_1, ta_yzzzz_xyyyzz_0, ta_yzzzz_xyyyzz_1, ta_yzzzz_xyyzz_0, ta_yzzzz_xyyzz_1, ta_yzzzz_xyyzzz_0, ta_yzzzz_xyyzzz_1, ta_yzzzz_xyzzz_0, ta_yzzzz_xyzzz_1, ta_yzzzz_xyzzzz_0, ta_yzzzz_xyzzzz_1, ta_yzzzz_yyyyy_0, ta_yzzzz_yyyyy_1, ta_yzzzz_yyyyyy_0, ta_yzzzz_yyyyyy_1, ta_yzzzz_yyyyyz_0, ta_yzzzz_yyyyyz_1, ta_yzzzz_yyyyz_0, ta_yzzzz_yyyyz_1, ta_yzzzz_yyyyzz_0, ta_yzzzz_yyyyzz_1, ta_yzzzz_yyyzz_0, ta_yzzzz_yyyzz_1, ta_yzzzz_yyyzzz_0, ta_yzzzz_yyyzzz_1, ta_yzzzz_yyzzz_0, ta_yzzzz_yyzzz_1, ta_yzzzz_yyzzzz_0, ta_yzzzz_yyzzzz_1, ta_yzzzz_yzzzz_0, ta_yzzzz_yzzzz_1, ta_yzzzz_yzzzzz_0, ta_yzzzz_yzzzzz_1, ta_yzzzz_zzzzzz_0, ta_yzzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzzz_xxxxxx_0[i] = ta_xzzzz_xxxxxx_0[i] * pa_y[i] - ta_xzzzz_xxxxxx_1[i] * pc_y[i];

        ta_xyzzzz_xxxxxy_0[i] = 5.0 * ta_yzzzz_xxxxy_0[i] * fe_0 - 5.0 * ta_yzzzz_xxxxy_1[i] * fe_0 + ta_yzzzz_xxxxxy_0[i] * pa_x[i] - ta_yzzzz_xxxxxy_1[i] * pc_x[i];

        ta_xyzzzz_xxxxxz_0[i] = ta_xzzzz_xxxxxz_0[i] * pa_y[i] - ta_xzzzz_xxxxxz_1[i] * pc_y[i];

        ta_xyzzzz_xxxxyy_0[i] = 4.0 * ta_yzzzz_xxxyy_0[i] * fe_0 - 4.0 * ta_yzzzz_xxxyy_1[i] * fe_0 + ta_yzzzz_xxxxyy_0[i] * pa_x[i] - ta_yzzzz_xxxxyy_1[i] * pc_x[i];

        ta_xyzzzz_xxxxyz_0[i] = 4.0 * ta_yzzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_yzzzz_xxxyz_1[i] * fe_0 + ta_yzzzz_xxxxyz_0[i] * pa_x[i] - ta_yzzzz_xxxxyz_1[i] * pc_x[i];

        ta_xyzzzz_xxxxzz_0[i] = ta_xzzzz_xxxxzz_0[i] * pa_y[i] - ta_xzzzz_xxxxzz_1[i] * pc_y[i];

        ta_xyzzzz_xxxyyy_0[i] = 3.0 * ta_yzzzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yzzzz_xxyyy_1[i] * fe_0 + ta_yzzzz_xxxyyy_0[i] * pa_x[i] - ta_yzzzz_xxxyyy_1[i] * pc_x[i];

        ta_xyzzzz_xxxyyz_0[i] = 3.0 * ta_yzzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yzzzz_xxyyz_1[i] * fe_0 + ta_yzzzz_xxxyyz_0[i] * pa_x[i] - ta_yzzzz_xxxyyz_1[i] * pc_x[i];

        ta_xyzzzz_xxxyzz_0[i] = 3.0 * ta_yzzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yzzzz_xxyzz_1[i] * fe_0 + ta_yzzzz_xxxyzz_0[i] * pa_x[i] - ta_yzzzz_xxxyzz_1[i] * pc_x[i];

        ta_xyzzzz_xxxzzz_0[i] = ta_xzzzz_xxxzzz_0[i] * pa_y[i] - ta_xzzzz_xxxzzz_1[i] * pc_y[i];

        ta_xyzzzz_xxyyyy_0[i] = 2.0 * ta_yzzzz_xyyyy_0[i] * fe_0 - 2.0 * ta_yzzzz_xyyyy_1[i] * fe_0 + ta_yzzzz_xxyyyy_0[i] * pa_x[i] - ta_yzzzz_xxyyyy_1[i] * pc_x[i];

        ta_xyzzzz_xxyyyz_0[i] = 2.0 * ta_yzzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyyyz_1[i] * fe_0 + ta_yzzzz_xxyyyz_0[i] * pa_x[i] - ta_yzzzz_xxyyyz_1[i] * pc_x[i];

        ta_xyzzzz_xxyyzz_0[i] = 2.0 * ta_yzzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyyzz_1[i] * fe_0 + ta_yzzzz_xxyyzz_0[i] * pa_x[i] - ta_yzzzz_xxyyzz_1[i] * pc_x[i];

        ta_xyzzzz_xxyzzz_0[i] = 2.0 * ta_yzzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyzzz_1[i] * fe_0 + ta_yzzzz_xxyzzz_0[i] * pa_x[i] - ta_yzzzz_xxyzzz_1[i] * pc_x[i];

        ta_xyzzzz_xxzzzz_0[i] = ta_xzzzz_xxzzzz_0[i] * pa_y[i] - ta_xzzzz_xxzzzz_1[i] * pc_y[i];

        ta_xyzzzz_xyyyyy_0[i] = ta_yzzzz_yyyyy_0[i] * fe_0 - ta_yzzzz_yyyyy_1[i] * fe_0 + ta_yzzzz_xyyyyy_0[i] * pa_x[i] - ta_yzzzz_xyyyyy_1[i] * pc_x[i];

        ta_xyzzzz_xyyyyz_0[i] = ta_yzzzz_yyyyz_0[i] * fe_0 - ta_yzzzz_yyyyz_1[i] * fe_0 + ta_yzzzz_xyyyyz_0[i] * pa_x[i] - ta_yzzzz_xyyyyz_1[i] * pc_x[i];

        ta_xyzzzz_xyyyzz_0[i] = ta_yzzzz_yyyzz_0[i] * fe_0 - ta_yzzzz_yyyzz_1[i] * fe_0 + ta_yzzzz_xyyyzz_0[i] * pa_x[i] - ta_yzzzz_xyyyzz_1[i] * pc_x[i];

        ta_xyzzzz_xyyzzz_0[i] = ta_yzzzz_yyzzz_0[i] * fe_0 - ta_yzzzz_yyzzz_1[i] * fe_0 + ta_yzzzz_xyyzzz_0[i] * pa_x[i] - ta_yzzzz_xyyzzz_1[i] * pc_x[i];

        ta_xyzzzz_xyzzzz_0[i] = ta_yzzzz_yzzzz_0[i] * fe_0 - ta_yzzzz_yzzzz_1[i] * fe_0 + ta_yzzzz_xyzzzz_0[i] * pa_x[i] - ta_yzzzz_xyzzzz_1[i] * pc_x[i];

        ta_xyzzzz_xzzzzz_0[i] = ta_xzzzz_xzzzzz_0[i] * pa_y[i] - ta_xzzzz_xzzzzz_1[i] * pc_y[i];

        ta_xyzzzz_yyyyyy_0[i] = ta_yzzzz_yyyyyy_0[i] * pa_x[i] - ta_yzzzz_yyyyyy_1[i] * pc_x[i];

        ta_xyzzzz_yyyyyz_0[i] = ta_yzzzz_yyyyyz_0[i] * pa_x[i] - ta_yzzzz_yyyyyz_1[i] * pc_x[i];

        ta_xyzzzz_yyyyzz_0[i] = ta_yzzzz_yyyyzz_0[i] * pa_x[i] - ta_yzzzz_yyyyzz_1[i] * pc_x[i];

        ta_xyzzzz_yyyzzz_0[i] = ta_yzzzz_yyyzzz_0[i] * pa_x[i] - ta_yzzzz_yyyzzz_1[i] * pc_x[i];

        ta_xyzzzz_yyzzzz_0[i] = ta_yzzzz_yyzzzz_0[i] * pa_x[i] - ta_yzzzz_yyzzzz_1[i] * pc_x[i];

        ta_xyzzzz_yzzzzz_0[i] = ta_yzzzz_yzzzzz_0[i] * pa_x[i] - ta_yzzzz_yzzzzz_1[i] * pc_x[i];

        ta_xyzzzz_zzzzzz_0[i] = ta_yzzzz_zzzzzz_0[i] * pa_x[i] - ta_yzzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 560-588 components of targeted buffer : II

    auto ta_xzzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 560);

    auto ta_xzzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 561);

    auto ta_xzzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 562);

    auto ta_xzzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 563);

    auto ta_xzzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 564);

    auto ta_xzzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 565);

    auto ta_xzzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 566);

    auto ta_xzzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 567);

    auto ta_xzzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 568);

    auto ta_xzzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 569);

    auto ta_xzzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 570);

    auto ta_xzzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 571);

    auto ta_xzzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 572);

    auto ta_xzzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 573);

    auto ta_xzzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 574);

    auto ta_xzzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 575);

    auto ta_xzzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 576);

    auto ta_xzzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 577);

    auto ta_xzzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 578);

    auto ta_xzzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 579);

    auto ta_xzzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 580);

    auto ta_xzzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 581);

    auto ta_xzzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 582);

    auto ta_xzzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 583);

    auto ta_xzzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 584);

    auto ta_xzzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 585);

    auto ta_xzzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 586);

    auto ta_xzzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 587);

    #pragma omp simd aligned(pa_x, pc_x, ta_xzzzzz_xxxxxx_0, ta_xzzzzz_xxxxxy_0, ta_xzzzzz_xxxxxz_0, ta_xzzzzz_xxxxyy_0, ta_xzzzzz_xxxxyz_0, ta_xzzzzz_xxxxzz_0, ta_xzzzzz_xxxyyy_0, ta_xzzzzz_xxxyyz_0, ta_xzzzzz_xxxyzz_0, ta_xzzzzz_xxxzzz_0, ta_xzzzzz_xxyyyy_0, ta_xzzzzz_xxyyyz_0, ta_xzzzzz_xxyyzz_0, ta_xzzzzz_xxyzzz_0, ta_xzzzzz_xxzzzz_0, ta_xzzzzz_xyyyyy_0, ta_xzzzzz_xyyyyz_0, ta_xzzzzz_xyyyzz_0, ta_xzzzzz_xyyzzz_0, ta_xzzzzz_xyzzzz_0, ta_xzzzzz_xzzzzz_0, ta_xzzzzz_yyyyyy_0, ta_xzzzzz_yyyyyz_0, ta_xzzzzz_yyyyzz_0, ta_xzzzzz_yyyzzz_0, ta_xzzzzz_yyzzzz_0, ta_xzzzzz_yzzzzz_0, ta_xzzzzz_zzzzzz_0, ta_zzzzz_xxxxx_0, ta_zzzzz_xxxxx_1, ta_zzzzz_xxxxxx_0, ta_zzzzz_xxxxxx_1, ta_zzzzz_xxxxxy_0, ta_zzzzz_xxxxxy_1, ta_zzzzz_xxxxxz_0, ta_zzzzz_xxxxxz_1, ta_zzzzz_xxxxy_0, ta_zzzzz_xxxxy_1, ta_zzzzz_xxxxyy_0, ta_zzzzz_xxxxyy_1, ta_zzzzz_xxxxyz_0, ta_zzzzz_xxxxyz_1, ta_zzzzz_xxxxz_0, ta_zzzzz_xxxxz_1, ta_zzzzz_xxxxzz_0, ta_zzzzz_xxxxzz_1, ta_zzzzz_xxxyy_0, ta_zzzzz_xxxyy_1, ta_zzzzz_xxxyyy_0, ta_zzzzz_xxxyyy_1, ta_zzzzz_xxxyyz_0, ta_zzzzz_xxxyyz_1, ta_zzzzz_xxxyz_0, ta_zzzzz_xxxyz_1, ta_zzzzz_xxxyzz_0, ta_zzzzz_xxxyzz_1, ta_zzzzz_xxxzz_0, ta_zzzzz_xxxzz_1, ta_zzzzz_xxxzzz_0, ta_zzzzz_xxxzzz_1, ta_zzzzz_xxyyy_0, ta_zzzzz_xxyyy_1, ta_zzzzz_xxyyyy_0, ta_zzzzz_xxyyyy_1, ta_zzzzz_xxyyyz_0, ta_zzzzz_xxyyyz_1, ta_zzzzz_xxyyz_0, ta_zzzzz_xxyyz_1, ta_zzzzz_xxyyzz_0, ta_zzzzz_xxyyzz_1, ta_zzzzz_xxyzz_0, ta_zzzzz_xxyzz_1, ta_zzzzz_xxyzzz_0, ta_zzzzz_xxyzzz_1, ta_zzzzz_xxzzz_0, ta_zzzzz_xxzzz_1, ta_zzzzz_xxzzzz_0, ta_zzzzz_xxzzzz_1, ta_zzzzz_xyyyy_0, ta_zzzzz_xyyyy_1, ta_zzzzz_xyyyyy_0, ta_zzzzz_xyyyyy_1, ta_zzzzz_xyyyyz_0, ta_zzzzz_xyyyyz_1, ta_zzzzz_xyyyz_0, ta_zzzzz_xyyyz_1, ta_zzzzz_xyyyzz_0, ta_zzzzz_xyyyzz_1, ta_zzzzz_xyyzz_0, ta_zzzzz_xyyzz_1, ta_zzzzz_xyyzzz_0, ta_zzzzz_xyyzzz_1, ta_zzzzz_xyzzz_0, ta_zzzzz_xyzzz_1, ta_zzzzz_xyzzzz_0, ta_zzzzz_xyzzzz_1, ta_zzzzz_xzzzz_0, ta_zzzzz_xzzzz_1, ta_zzzzz_xzzzzz_0, ta_zzzzz_xzzzzz_1, ta_zzzzz_yyyyy_0, ta_zzzzz_yyyyy_1, ta_zzzzz_yyyyyy_0, ta_zzzzz_yyyyyy_1, ta_zzzzz_yyyyyz_0, ta_zzzzz_yyyyyz_1, ta_zzzzz_yyyyz_0, ta_zzzzz_yyyyz_1, ta_zzzzz_yyyyzz_0, ta_zzzzz_yyyyzz_1, ta_zzzzz_yyyzz_0, ta_zzzzz_yyyzz_1, ta_zzzzz_yyyzzz_0, ta_zzzzz_yyyzzz_1, ta_zzzzz_yyzzz_0, ta_zzzzz_yyzzz_1, ta_zzzzz_yyzzzz_0, ta_zzzzz_yyzzzz_1, ta_zzzzz_yzzzz_0, ta_zzzzz_yzzzz_1, ta_zzzzz_yzzzzz_0, ta_zzzzz_yzzzzz_1, ta_zzzzz_zzzzz_0, ta_zzzzz_zzzzz_1, ta_zzzzz_zzzzzz_0, ta_zzzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzzz_xxxxxx_0[i] = 6.0 * ta_zzzzz_xxxxx_0[i] * fe_0 - 6.0 * ta_zzzzz_xxxxx_1[i] * fe_0 + ta_zzzzz_xxxxxx_0[i] * pa_x[i] - ta_zzzzz_xxxxxx_1[i] * pc_x[i];

        ta_xzzzzz_xxxxxy_0[i] = 5.0 * ta_zzzzz_xxxxy_0[i] * fe_0 - 5.0 * ta_zzzzz_xxxxy_1[i] * fe_0 + ta_zzzzz_xxxxxy_0[i] * pa_x[i] - ta_zzzzz_xxxxxy_1[i] * pc_x[i];

        ta_xzzzzz_xxxxxz_0[i] = 5.0 * ta_zzzzz_xxxxz_0[i] * fe_0 - 5.0 * ta_zzzzz_xxxxz_1[i] * fe_0 + ta_zzzzz_xxxxxz_0[i] * pa_x[i] - ta_zzzzz_xxxxxz_1[i] * pc_x[i];

        ta_xzzzzz_xxxxyy_0[i] = 4.0 * ta_zzzzz_xxxyy_0[i] * fe_0 - 4.0 * ta_zzzzz_xxxyy_1[i] * fe_0 + ta_zzzzz_xxxxyy_0[i] * pa_x[i] - ta_zzzzz_xxxxyy_1[i] * pc_x[i];

        ta_xzzzzz_xxxxyz_0[i] = 4.0 * ta_zzzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_zzzzz_xxxyz_1[i] * fe_0 + ta_zzzzz_xxxxyz_0[i] * pa_x[i] - ta_zzzzz_xxxxyz_1[i] * pc_x[i];

        ta_xzzzzz_xxxxzz_0[i] = 4.0 * ta_zzzzz_xxxzz_0[i] * fe_0 - 4.0 * ta_zzzzz_xxxzz_1[i] * fe_0 + ta_zzzzz_xxxxzz_0[i] * pa_x[i] - ta_zzzzz_xxxxzz_1[i] * pc_x[i];

        ta_xzzzzz_xxxyyy_0[i] = 3.0 * ta_zzzzz_xxyyy_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyyy_1[i] * fe_0 + ta_zzzzz_xxxyyy_0[i] * pa_x[i] - ta_zzzzz_xxxyyy_1[i] * pc_x[i];

        ta_xzzzzz_xxxyyz_0[i] = 3.0 * ta_zzzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyyz_1[i] * fe_0 + ta_zzzzz_xxxyyz_0[i] * pa_x[i] - ta_zzzzz_xxxyyz_1[i] * pc_x[i];

        ta_xzzzzz_xxxyzz_0[i] = 3.0 * ta_zzzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyzz_1[i] * fe_0 + ta_zzzzz_xxxyzz_0[i] * pa_x[i] - ta_zzzzz_xxxyzz_1[i] * pc_x[i];

        ta_xzzzzz_xxxzzz_0[i] = 3.0 * ta_zzzzz_xxzzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxzzz_1[i] * fe_0 + ta_zzzzz_xxxzzz_0[i] * pa_x[i] - ta_zzzzz_xxxzzz_1[i] * pc_x[i];

        ta_xzzzzz_xxyyyy_0[i] = 2.0 * ta_zzzzz_xyyyy_0[i] * fe_0 - 2.0 * ta_zzzzz_xyyyy_1[i] * fe_0 + ta_zzzzz_xxyyyy_0[i] * pa_x[i] - ta_zzzzz_xxyyyy_1[i] * pc_x[i];

        ta_xzzzzz_xxyyyz_0[i] = 2.0 * ta_zzzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyyyz_1[i] * fe_0 + ta_zzzzz_xxyyyz_0[i] * pa_x[i] - ta_zzzzz_xxyyyz_1[i] * pc_x[i];

        ta_xzzzzz_xxyyzz_0[i] = 2.0 * ta_zzzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyyzz_1[i] * fe_0 + ta_zzzzz_xxyyzz_0[i] * pa_x[i] - ta_zzzzz_xxyyzz_1[i] * pc_x[i];

        ta_xzzzzz_xxyzzz_0[i] = 2.0 * ta_zzzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyzzz_1[i] * fe_0 + ta_zzzzz_xxyzzz_0[i] * pa_x[i] - ta_zzzzz_xxyzzz_1[i] * pc_x[i];

        ta_xzzzzz_xxzzzz_0[i] = 2.0 * ta_zzzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xzzzz_1[i] * fe_0 + ta_zzzzz_xxzzzz_0[i] * pa_x[i] - ta_zzzzz_xxzzzz_1[i] * pc_x[i];

        ta_xzzzzz_xyyyyy_0[i] = ta_zzzzz_yyyyy_0[i] * fe_0 - ta_zzzzz_yyyyy_1[i] * fe_0 + ta_zzzzz_xyyyyy_0[i] * pa_x[i] - ta_zzzzz_xyyyyy_1[i] * pc_x[i];

        ta_xzzzzz_xyyyyz_0[i] = ta_zzzzz_yyyyz_0[i] * fe_0 - ta_zzzzz_yyyyz_1[i] * fe_0 + ta_zzzzz_xyyyyz_0[i] * pa_x[i] - ta_zzzzz_xyyyyz_1[i] * pc_x[i];

        ta_xzzzzz_xyyyzz_0[i] = ta_zzzzz_yyyzz_0[i] * fe_0 - ta_zzzzz_yyyzz_1[i] * fe_0 + ta_zzzzz_xyyyzz_0[i] * pa_x[i] - ta_zzzzz_xyyyzz_1[i] * pc_x[i];

        ta_xzzzzz_xyyzzz_0[i] = ta_zzzzz_yyzzz_0[i] * fe_0 - ta_zzzzz_yyzzz_1[i] * fe_0 + ta_zzzzz_xyyzzz_0[i] * pa_x[i] - ta_zzzzz_xyyzzz_1[i] * pc_x[i];

        ta_xzzzzz_xyzzzz_0[i] = ta_zzzzz_yzzzz_0[i] * fe_0 - ta_zzzzz_yzzzz_1[i] * fe_0 + ta_zzzzz_xyzzzz_0[i] * pa_x[i] - ta_zzzzz_xyzzzz_1[i] * pc_x[i];

        ta_xzzzzz_xzzzzz_0[i] = ta_zzzzz_zzzzz_0[i] * fe_0 - ta_zzzzz_zzzzz_1[i] * fe_0 + ta_zzzzz_xzzzzz_0[i] * pa_x[i] - ta_zzzzz_xzzzzz_1[i] * pc_x[i];

        ta_xzzzzz_yyyyyy_0[i] = ta_zzzzz_yyyyyy_0[i] * pa_x[i] - ta_zzzzz_yyyyyy_1[i] * pc_x[i];

        ta_xzzzzz_yyyyyz_0[i] = ta_zzzzz_yyyyyz_0[i] * pa_x[i] - ta_zzzzz_yyyyyz_1[i] * pc_x[i];

        ta_xzzzzz_yyyyzz_0[i] = ta_zzzzz_yyyyzz_0[i] * pa_x[i] - ta_zzzzz_yyyyzz_1[i] * pc_x[i];

        ta_xzzzzz_yyyzzz_0[i] = ta_zzzzz_yyyzzz_0[i] * pa_x[i] - ta_zzzzz_yyyzzz_1[i] * pc_x[i];

        ta_xzzzzz_yyzzzz_0[i] = ta_zzzzz_yyzzzz_0[i] * pa_x[i] - ta_zzzzz_yyzzzz_1[i] * pc_x[i];

        ta_xzzzzz_yzzzzz_0[i] = ta_zzzzz_yzzzzz_0[i] * pa_x[i] - ta_zzzzz_yzzzzz_1[i] * pc_x[i];

        ta_xzzzzz_zzzzzz_0[i] = ta_zzzzz_zzzzzz_0[i] * pa_x[i] - ta_zzzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 588-616 components of targeted buffer : II

    auto ta_yyyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 588);

    auto ta_yyyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 589);

    auto ta_yyyyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 590);

    auto ta_yyyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 591);

    auto ta_yyyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 592);

    auto ta_yyyyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 593);

    auto ta_yyyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 594);

    auto ta_yyyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 595);

    auto ta_yyyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 596);

    auto ta_yyyyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 597);

    auto ta_yyyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 598);

    auto ta_yyyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 599);

    auto ta_yyyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 600);

    auto ta_yyyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 601);

    auto ta_yyyyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 602);

    auto ta_yyyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 603);

    auto ta_yyyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 604);

    auto ta_yyyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 605);

    auto ta_yyyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 606);

    auto ta_yyyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 607);

    auto ta_yyyyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 608);

    auto ta_yyyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 609);

    auto ta_yyyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 610);

    auto ta_yyyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 611);

    auto ta_yyyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 612);

    auto ta_yyyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 613);

    auto ta_yyyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 614);

    auto ta_yyyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 615);

    #pragma omp simd aligned(pa_y, pc_y, ta_yyyy_xxxxxx_0, ta_yyyy_xxxxxx_1, ta_yyyy_xxxxxy_0, ta_yyyy_xxxxxy_1, ta_yyyy_xxxxxz_0, ta_yyyy_xxxxxz_1, ta_yyyy_xxxxyy_0, ta_yyyy_xxxxyy_1, ta_yyyy_xxxxyz_0, ta_yyyy_xxxxyz_1, ta_yyyy_xxxxzz_0, ta_yyyy_xxxxzz_1, ta_yyyy_xxxyyy_0, ta_yyyy_xxxyyy_1, ta_yyyy_xxxyyz_0, ta_yyyy_xxxyyz_1, ta_yyyy_xxxyzz_0, ta_yyyy_xxxyzz_1, ta_yyyy_xxxzzz_0, ta_yyyy_xxxzzz_1, ta_yyyy_xxyyyy_0, ta_yyyy_xxyyyy_1, ta_yyyy_xxyyyz_0, ta_yyyy_xxyyyz_1, ta_yyyy_xxyyzz_0, ta_yyyy_xxyyzz_1, ta_yyyy_xxyzzz_0, ta_yyyy_xxyzzz_1, ta_yyyy_xxzzzz_0, ta_yyyy_xxzzzz_1, ta_yyyy_xyyyyy_0, ta_yyyy_xyyyyy_1, ta_yyyy_xyyyyz_0, ta_yyyy_xyyyyz_1, ta_yyyy_xyyyzz_0, ta_yyyy_xyyyzz_1, ta_yyyy_xyyzzz_0, ta_yyyy_xyyzzz_1, ta_yyyy_xyzzzz_0, ta_yyyy_xyzzzz_1, ta_yyyy_xzzzzz_0, ta_yyyy_xzzzzz_1, ta_yyyy_yyyyyy_0, ta_yyyy_yyyyyy_1, ta_yyyy_yyyyyz_0, ta_yyyy_yyyyyz_1, ta_yyyy_yyyyzz_0, ta_yyyy_yyyyzz_1, ta_yyyy_yyyzzz_0, ta_yyyy_yyyzzz_1, ta_yyyy_yyzzzz_0, ta_yyyy_yyzzzz_1, ta_yyyy_yzzzzz_0, ta_yyyy_yzzzzz_1, ta_yyyy_zzzzzz_0, ta_yyyy_zzzzzz_1, ta_yyyyy_xxxxx_0, ta_yyyyy_xxxxx_1, ta_yyyyy_xxxxxx_0, ta_yyyyy_xxxxxx_1, ta_yyyyy_xxxxxy_0, ta_yyyyy_xxxxxy_1, ta_yyyyy_xxxxxz_0, ta_yyyyy_xxxxxz_1, ta_yyyyy_xxxxy_0, ta_yyyyy_xxxxy_1, ta_yyyyy_xxxxyy_0, ta_yyyyy_xxxxyy_1, ta_yyyyy_xxxxyz_0, ta_yyyyy_xxxxyz_1, ta_yyyyy_xxxxz_0, ta_yyyyy_xxxxz_1, ta_yyyyy_xxxxzz_0, ta_yyyyy_xxxxzz_1, ta_yyyyy_xxxyy_0, ta_yyyyy_xxxyy_1, ta_yyyyy_xxxyyy_0, ta_yyyyy_xxxyyy_1, ta_yyyyy_xxxyyz_0, ta_yyyyy_xxxyyz_1, ta_yyyyy_xxxyz_0, ta_yyyyy_xxxyz_1, ta_yyyyy_xxxyzz_0, ta_yyyyy_xxxyzz_1, ta_yyyyy_xxxzz_0, ta_yyyyy_xxxzz_1, ta_yyyyy_xxxzzz_0, ta_yyyyy_xxxzzz_1, ta_yyyyy_xxyyy_0, ta_yyyyy_xxyyy_1, ta_yyyyy_xxyyyy_0, ta_yyyyy_xxyyyy_1, ta_yyyyy_xxyyyz_0, ta_yyyyy_xxyyyz_1, ta_yyyyy_xxyyz_0, ta_yyyyy_xxyyz_1, ta_yyyyy_xxyyzz_0, ta_yyyyy_xxyyzz_1, ta_yyyyy_xxyzz_0, ta_yyyyy_xxyzz_1, ta_yyyyy_xxyzzz_0, ta_yyyyy_xxyzzz_1, ta_yyyyy_xxzzz_0, ta_yyyyy_xxzzz_1, ta_yyyyy_xxzzzz_0, ta_yyyyy_xxzzzz_1, ta_yyyyy_xyyyy_0, ta_yyyyy_xyyyy_1, ta_yyyyy_xyyyyy_0, ta_yyyyy_xyyyyy_1, ta_yyyyy_xyyyyz_0, ta_yyyyy_xyyyyz_1, ta_yyyyy_xyyyz_0, ta_yyyyy_xyyyz_1, ta_yyyyy_xyyyzz_0, ta_yyyyy_xyyyzz_1, ta_yyyyy_xyyzz_0, ta_yyyyy_xyyzz_1, ta_yyyyy_xyyzzz_0, ta_yyyyy_xyyzzz_1, ta_yyyyy_xyzzz_0, ta_yyyyy_xyzzz_1, ta_yyyyy_xyzzzz_0, ta_yyyyy_xyzzzz_1, ta_yyyyy_xzzzz_0, ta_yyyyy_xzzzz_1, ta_yyyyy_xzzzzz_0, ta_yyyyy_xzzzzz_1, ta_yyyyy_yyyyy_0, ta_yyyyy_yyyyy_1, ta_yyyyy_yyyyyy_0, ta_yyyyy_yyyyyy_1, ta_yyyyy_yyyyyz_0, ta_yyyyy_yyyyyz_1, ta_yyyyy_yyyyz_0, ta_yyyyy_yyyyz_1, ta_yyyyy_yyyyzz_0, ta_yyyyy_yyyyzz_1, ta_yyyyy_yyyzz_0, ta_yyyyy_yyyzz_1, ta_yyyyy_yyyzzz_0, ta_yyyyy_yyyzzz_1, ta_yyyyy_yyzzz_0, ta_yyyyy_yyzzz_1, ta_yyyyy_yyzzzz_0, ta_yyyyy_yyzzzz_1, ta_yyyyy_yzzzz_0, ta_yyyyy_yzzzz_1, ta_yyyyy_yzzzzz_0, ta_yyyyy_yzzzzz_1, ta_yyyyy_zzzzz_0, ta_yyyyy_zzzzz_1, ta_yyyyy_zzzzzz_0, ta_yyyyy_zzzzzz_1, ta_yyyyyy_xxxxxx_0, ta_yyyyyy_xxxxxy_0, ta_yyyyyy_xxxxxz_0, ta_yyyyyy_xxxxyy_0, ta_yyyyyy_xxxxyz_0, ta_yyyyyy_xxxxzz_0, ta_yyyyyy_xxxyyy_0, ta_yyyyyy_xxxyyz_0, ta_yyyyyy_xxxyzz_0, ta_yyyyyy_xxxzzz_0, ta_yyyyyy_xxyyyy_0, ta_yyyyyy_xxyyyz_0, ta_yyyyyy_xxyyzz_0, ta_yyyyyy_xxyzzz_0, ta_yyyyyy_xxzzzz_0, ta_yyyyyy_xyyyyy_0, ta_yyyyyy_xyyyyz_0, ta_yyyyyy_xyyyzz_0, ta_yyyyyy_xyyzzz_0, ta_yyyyyy_xyzzzz_0, ta_yyyyyy_xzzzzz_0, ta_yyyyyy_yyyyyy_0, ta_yyyyyy_yyyyyz_0, ta_yyyyyy_yyyyzz_0, ta_yyyyyy_yyyzzz_0, ta_yyyyyy_yyzzzz_0, ta_yyyyyy_yzzzzz_0, ta_yyyyyy_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyy_xxxxxx_0[i] = 5.0 * ta_yyyy_xxxxxx_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxxx_1[i] * fe_0 + ta_yyyyy_xxxxxx_0[i] * pa_y[i] - ta_yyyyy_xxxxxx_1[i] * pc_y[i];

        ta_yyyyyy_xxxxxy_0[i] = 5.0 * ta_yyyy_xxxxxy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxxy_1[i] * fe_0 + ta_yyyyy_xxxxx_0[i] * fe_0 - ta_yyyyy_xxxxx_1[i] * fe_0 + ta_yyyyy_xxxxxy_0[i] * pa_y[i] - ta_yyyyy_xxxxxy_1[i] * pc_y[i];

        ta_yyyyyy_xxxxxz_0[i] = 5.0 * ta_yyyy_xxxxxz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxxz_1[i] * fe_0 + ta_yyyyy_xxxxxz_0[i] * pa_y[i] - ta_yyyyy_xxxxxz_1[i] * pc_y[i];

        ta_yyyyyy_xxxxyy_0[i] = 5.0 * ta_yyyy_xxxxyy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxyy_1[i] * fe_0 + 2.0 * ta_yyyyy_xxxxy_0[i] * fe_0 - 2.0 * ta_yyyyy_xxxxy_1[i] * fe_0 + ta_yyyyy_xxxxyy_0[i] * pa_y[i] - ta_yyyyy_xxxxyy_1[i] * pc_y[i];

        ta_yyyyyy_xxxxyz_0[i] = 5.0 * ta_yyyy_xxxxyz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxyz_1[i] * fe_0 + ta_yyyyy_xxxxz_0[i] * fe_0 - ta_yyyyy_xxxxz_1[i] * fe_0 + ta_yyyyy_xxxxyz_0[i] * pa_y[i] - ta_yyyyy_xxxxyz_1[i] * pc_y[i];

        ta_yyyyyy_xxxxzz_0[i] = 5.0 * ta_yyyy_xxxxzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxzz_1[i] * fe_0 + ta_yyyyy_xxxxzz_0[i] * pa_y[i] - ta_yyyyy_xxxxzz_1[i] * pc_y[i];

        ta_yyyyyy_xxxyyy_0[i] = 5.0 * ta_yyyy_xxxyyy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_yyyyy_xxxyy_0[i] * fe_0 - 3.0 * ta_yyyyy_xxxyy_1[i] * fe_0 + ta_yyyyy_xxxyyy_0[i] * pa_y[i] - ta_yyyyy_xxxyyy_1[i] * pc_y[i];

        ta_yyyyyy_xxxyyz_0[i] = 5.0 * ta_yyyy_xxxyyz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxyyz_1[i] * fe_0 + 2.0 * ta_yyyyy_xxxyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xxxyz_1[i] * fe_0 + ta_yyyyy_xxxyyz_0[i] * pa_y[i] - ta_yyyyy_xxxyyz_1[i] * pc_y[i];

        ta_yyyyyy_xxxyzz_0[i] = 5.0 * ta_yyyy_xxxyzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxyzz_1[i] * fe_0 + ta_yyyyy_xxxzz_0[i] * fe_0 - ta_yyyyy_xxxzz_1[i] * fe_0 + ta_yyyyy_xxxyzz_0[i] * pa_y[i] - ta_yyyyy_xxxyzz_1[i] * pc_y[i];

        ta_yyyyyy_xxxzzz_0[i] = 5.0 * ta_yyyy_xxxzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxzzz_1[i] * fe_0 + ta_yyyyy_xxxzzz_0[i] * pa_y[i] - ta_yyyyy_xxxzzz_1[i] * pc_y[i];

        ta_yyyyyy_xxyyyy_0[i] = 5.0 * ta_yyyy_xxyyyy_0[i] * fe_0 - 5.0 * ta_yyyy_xxyyyy_1[i] * fe_0 + 4.0 * ta_yyyyy_xxyyy_0[i] * fe_0 - 4.0 * ta_yyyyy_xxyyy_1[i] * fe_0 + ta_yyyyy_xxyyyy_0[i] * pa_y[i] - ta_yyyyy_xxyyyy_1[i] * pc_y[i];

        ta_yyyyyy_xxyyyz_0[i] = 5.0 * ta_yyyy_xxyyyz_0[i] * fe_0 - 5.0 * ta_yyyy_xxyyyz_1[i] * fe_0 + 3.0 * ta_yyyyy_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyyz_1[i] * fe_0 + ta_yyyyy_xxyyyz_0[i] * pa_y[i] - ta_yyyyy_xxyyyz_1[i] * pc_y[i];

        ta_yyyyyy_xxyyzz_0[i] = 5.0 * ta_yyyy_xxyyzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_yyyyy_xxyzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xxyzz_1[i] * fe_0 + ta_yyyyy_xxyyzz_0[i] * pa_y[i] - ta_yyyyy_xxyyzz_1[i] * pc_y[i];

        ta_yyyyyy_xxyzzz_0[i] = 5.0 * ta_yyyy_xxyzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxyzzz_1[i] * fe_0 + ta_yyyyy_xxzzz_0[i] * fe_0 - ta_yyyyy_xxzzz_1[i] * fe_0 + ta_yyyyy_xxyzzz_0[i] * pa_y[i] - ta_yyyyy_xxyzzz_1[i] * pc_y[i];

        ta_yyyyyy_xxzzzz_0[i] = 5.0 * ta_yyyy_xxzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxzzzz_1[i] * fe_0 + ta_yyyyy_xxzzzz_0[i] * pa_y[i] - ta_yyyyy_xxzzzz_1[i] * pc_y[i];

        ta_yyyyyy_xyyyyy_0[i] = 5.0 * ta_yyyy_xyyyyy_0[i] * fe_0 - 5.0 * ta_yyyy_xyyyyy_1[i] * fe_0 + 5.0 * ta_yyyyy_xyyyy_0[i] * fe_0 - 5.0 * ta_yyyyy_xyyyy_1[i] * fe_0 + ta_yyyyy_xyyyyy_0[i] * pa_y[i] - ta_yyyyy_xyyyyy_1[i] * pc_y[i];

        ta_yyyyyy_xyyyyz_0[i] = 5.0 * ta_yyyy_xyyyyz_0[i] * fe_0 - 5.0 * ta_yyyy_xyyyyz_1[i] * fe_0 + 4.0 * ta_yyyyy_xyyyz_0[i] * fe_0 - 4.0 * ta_yyyyy_xyyyz_1[i] * fe_0 + ta_yyyyy_xyyyyz_0[i] * pa_y[i] - ta_yyyyy_xyyyyz_1[i] * pc_y[i];

        ta_yyyyyy_xyyyzz_0[i] = 5.0 * ta_yyyy_xyyyzz_0[i] * fe_0 - 5.0 * ta_yyyy_xyyyzz_1[i] * fe_0 + 3.0 * ta_yyyyy_xyyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xyyzz_1[i] * fe_0 + ta_yyyyy_xyyyzz_0[i] * pa_y[i] - ta_yyyyy_xyyyzz_1[i] * pc_y[i];

        ta_yyyyyy_xyyzzz_0[i] = 5.0 * ta_yyyy_xyyzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xyyzzz_1[i] * fe_0 + 2.0 * ta_yyyyy_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyzzz_1[i] * fe_0 + ta_yyyyy_xyyzzz_0[i] * pa_y[i] - ta_yyyyy_xyyzzz_1[i] * pc_y[i];

        ta_yyyyyy_xyzzzz_0[i] = 5.0 * ta_yyyy_xyzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xyzzzz_1[i] * fe_0 + ta_yyyyy_xzzzz_0[i] * fe_0 - ta_yyyyy_xzzzz_1[i] * fe_0 + ta_yyyyy_xyzzzz_0[i] * pa_y[i] - ta_yyyyy_xyzzzz_1[i] * pc_y[i];

        ta_yyyyyy_xzzzzz_0[i] = 5.0 * ta_yyyy_xzzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xzzzzz_1[i] * fe_0 + ta_yyyyy_xzzzzz_0[i] * pa_y[i] - ta_yyyyy_xzzzzz_1[i] * pc_y[i];

        ta_yyyyyy_yyyyyy_0[i] = 5.0 * ta_yyyy_yyyyyy_0[i] * fe_0 - 5.0 * ta_yyyy_yyyyyy_1[i] * fe_0 + 6.0 * ta_yyyyy_yyyyy_0[i] * fe_0 - 6.0 * ta_yyyyy_yyyyy_1[i] * fe_0 + ta_yyyyy_yyyyyy_0[i] * pa_y[i] - ta_yyyyy_yyyyyy_1[i] * pc_y[i];

        ta_yyyyyy_yyyyyz_0[i] = 5.0 * ta_yyyy_yyyyyz_0[i] * fe_0 - 5.0 * ta_yyyy_yyyyyz_1[i] * fe_0 + 5.0 * ta_yyyyy_yyyyz_0[i] * fe_0 - 5.0 * ta_yyyyy_yyyyz_1[i] * fe_0 + ta_yyyyy_yyyyyz_0[i] * pa_y[i] - ta_yyyyy_yyyyyz_1[i] * pc_y[i];

        ta_yyyyyy_yyyyzz_0[i] = 5.0 * ta_yyyy_yyyyzz_0[i] * fe_0 - 5.0 * ta_yyyy_yyyyzz_1[i] * fe_0 + 4.0 * ta_yyyyy_yyyzz_0[i] * fe_0 - 4.0 * ta_yyyyy_yyyzz_1[i] * fe_0 + ta_yyyyy_yyyyzz_0[i] * pa_y[i] - ta_yyyyy_yyyyzz_1[i] * pc_y[i];

        ta_yyyyyy_yyyzzz_0[i] = 5.0 * ta_yyyy_yyyzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yyyzzz_1[i] * fe_0 + 3.0 * ta_yyyyy_yyzzz_0[i] * fe_0 - 3.0 * ta_yyyyy_yyzzz_1[i] * fe_0 + ta_yyyyy_yyyzzz_0[i] * pa_y[i] - ta_yyyyy_yyyzzz_1[i] * pc_y[i];

        ta_yyyyyy_yyzzzz_0[i] = 5.0 * ta_yyyy_yyzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yyzzzz_1[i] * fe_0 + 2.0 * ta_yyyyy_yzzzz_0[i] * fe_0 - 2.0 * ta_yyyyy_yzzzz_1[i] * fe_0 + ta_yyyyy_yyzzzz_0[i] * pa_y[i] - ta_yyyyy_yyzzzz_1[i] * pc_y[i];

        ta_yyyyyy_yzzzzz_0[i] = 5.0 * ta_yyyy_yzzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yzzzzz_1[i] * fe_0 + ta_yyyyy_zzzzz_0[i] * fe_0 - ta_yyyyy_zzzzz_1[i] * fe_0 + ta_yyyyy_yzzzzz_0[i] * pa_y[i] - ta_yyyyy_yzzzzz_1[i] * pc_y[i];

        ta_yyyyyy_zzzzzz_0[i] = 5.0 * ta_yyyy_zzzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_zzzzzz_1[i] * fe_0 + ta_yyyyy_zzzzzz_0[i] * pa_y[i] - ta_yyyyy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 616-644 components of targeted buffer : II

    auto ta_yyyyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 616);

    auto ta_yyyyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 617);

    auto ta_yyyyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 618);

    auto ta_yyyyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 619);

    auto ta_yyyyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 620);

    auto ta_yyyyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 621);

    auto ta_yyyyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 622);

    auto ta_yyyyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 623);

    auto ta_yyyyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 624);

    auto ta_yyyyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 625);

    auto ta_yyyyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 626);

    auto ta_yyyyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 627);

    auto ta_yyyyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 628);

    auto ta_yyyyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 629);

    auto ta_yyyyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 630);

    auto ta_yyyyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 631);

    auto ta_yyyyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 632);

    auto ta_yyyyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 633);

    auto ta_yyyyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 634);

    auto ta_yyyyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 635);

    auto ta_yyyyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 636);

    auto ta_yyyyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 637);

    auto ta_yyyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 638);

    auto ta_yyyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 639);

    auto ta_yyyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 640);

    auto ta_yyyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 641);

    auto ta_yyyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 642);

    auto ta_yyyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 643);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yyyyy_xxxxxx_0, ta_yyyyy_xxxxxx_1, ta_yyyyy_xxxxxy_0, ta_yyyyy_xxxxxy_1, ta_yyyyy_xxxxy_0, ta_yyyyy_xxxxy_1, ta_yyyyy_xxxxyy_0, ta_yyyyy_xxxxyy_1, ta_yyyyy_xxxxyz_0, ta_yyyyy_xxxxyz_1, ta_yyyyy_xxxyy_0, ta_yyyyy_xxxyy_1, ta_yyyyy_xxxyyy_0, ta_yyyyy_xxxyyy_1, ta_yyyyy_xxxyyz_0, ta_yyyyy_xxxyyz_1, ta_yyyyy_xxxyz_0, ta_yyyyy_xxxyz_1, ta_yyyyy_xxxyzz_0, ta_yyyyy_xxxyzz_1, ta_yyyyy_xxyyy_0, ta_yyyyy_xxyyy_1, ta_yyyyy_xxyyyy_0, ta_yyyyy_xxyyyy_1, ta_yyyyy_xxyyyz_0, ta_yyyyy_xxyyyz_1, ta_yyyyy_xxyyz_0, ta_yyyyy_xxyyz_1, ta_yyyyy_xxyyzz_0, ta_yyyyy_xxyyzz_1, ta_yyyyy_xxyzz_0, ta_yyyyy_xxyzz_1, ta_yyyyy_xxyzzz_0, ta_yyyyy_xxyzzz_1, ta_yyyyy_xyyyy_0, ta_yyyyy_xyyyy_1, ta_yyyyy_xyyyyy_0, ta_yyyyy_xyyyyy_1, ta_yyyyy_xyyyyz_0, ta_yyyyy_xyyyyz_1, ta_yyyyy_xyyyz_0, ta_yyyyy_xyyyz_1, ta_yyyyy_xyyyzz_0, ta_yyyyy_xyyyzz_1, ta_yyyyy_xyyzz_0, ta_yyyyy_xyyzz_1, ta_yyyyy_xyyzzz_0, ta_yyyyy_xyyzzz_1, ta_yyyyy_xyzzz_0, ta_yyyyy_xyzzz_1, ta_yyyyy_xyzzzz_0, ta_yyyyy_xyzzzz_1, ta_yyyyy_yyyyy_0, ta_yyyyy_yyyyy_1, ta_yyyyy_yyyyyy_0, ta_yyyyy_yyyyyy_1, ta_yyyyy_yyyyyz_0, ta_yyyyy_yyyyyz_1, ta_yyyyy_yyyyz_0, ta_yyyyy_yyyyz_1, ta_yyyyy_yyyyzz_0, ta_yyyyy_yyyyzz_1, ta_yyyyy_yyyzz_0, ta_yyyyy_yyyzz_1, ta_yyyyy_yyyzzz_0, ta_yyyyy_yyyzzz_1, ta_yyyyy_yyzzz_0, ta_yyyyy_yyzzz_1, ta_yyyyy_yyzzzz_0, ta_yyyyy_yyzzzz_1, ta_yyyyy_yzzzz_0, ta_yyyyy_yzzzz_1, ta_yyyyy_yzzzzz_0, ta_yyyyy_yzzzzz_1, ta_yyyyyz_xxxxxx_0, ta_yyyyyz_xxxxxy_0, ta_yyyyyz_xxxxxz_0, ta_yyyyyz_xxxxyy_0, ta_yyyyyz_xxxxyz_0, ta_yyyyyz_xxxxzz_0, ta_yyyyyz_xxxyyy_0, ta_yyyyyz_xxxyyz_0, ta_yyyyyz_xxxyzz_0, ta_yyyyyz_xxxzzz_0, ta_yyyyyz_xxyyyy_0, ta_yyyyyz_xxyyyz_0, ta_yyyyyz_xxyyzz_0, ta_yyyyyz_xxyzzz_0, ta_yyyyyz_xxzzzz_0, ta_yyyyyz_xyyyyy_0, ta_yyyyyz_xyyyyz_0, ta_yyyyyz_xyyyzz_0, ta_yyyyyz_xyyzzz_0, ta_yyyyyz_xyzzzz_0, ta_yyyyyz_xzzzzz_0, ta_yyyyyz_yyyyyy_0, ta_yyyyyz_yyyyyz_0, ta_yyyyyz_yyyyzz_0, ta_yyyyyz_yyyzzz_0, ta_yyyyyz_yyzzzz_0, ta_yyyyyz_yzzzzz_0, ta_yyyyyz_zzzzzz_0, ta_yyyyz_xxxxxz_0, ta_yyyyz_xxxxxz_1, ta_yyyyz_xxxxzz_0, ta_yyyyz_xxxxzz_1, ta_yyyyz_xxxzzz_0, ta_yyyyz_xxxzzz_1, ta_yyyyz_xxzzzz_0, ta_yyyyz_xxzzzz_1, ta_yyyyz_xzzzzz_0, ta_yyyyz_xzzzzz_1, ta_yyyyz_zzzzzz_0, ta_yyyyz_zzzzzz_1, ta_yyyz_xxxxxz_0, ta_yyyz_xxxxxz_1, ta_yyyz_xxxxzz_0, ta_yyyz_xxxxzz_1, ta_yyyz_xxxzzz_0, ta_yyyz_xxxzzz_1, ta_yyyz_xxzzzz_0, ta_yyyz_xxzzzz_1, ta_yyyz_xzzzzz_0, ta_yyyz_xzzzzz_1, ta_yyyz_zzzzzz_0, ta_yyyz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyz_xxxxxx_0[i] = ta_yyyyy_xxxxxx_0[i] * pa_z[i] - ta_yyyyy_xxxxxx_1[i] * pc_z[i];

        ta_yyyyyz_xxxxxy_0[i] = ta_yyyyy_xxxxxy_0[i] * pa_z[i] - ta_yyyyy_xxxxxy_1[i] * pc_z[i];

        ta_yyyyyz_xxxxxz_0[i] = 4.0 * ta_yyyz_xxxxxz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxxxz_1[i] * fe_0 + ta_yyyyz_xxxxxz_0[i] * pa_y[i] - ta_yyyyz_xxxxxz_1[i] * pc_y[i];

        ta_yyyyyz_xxxxyy_0[i] = ta_yyyyy_xxxxyy_0[i] * pa_z[i] - ta_yyyyy_xxxxyy_1[i] * pc_z[i];

        ta_yyyyyz_xxxxyz_0[i] = ta_yyyyy_xxxxy_0[i] * fe_0 - ta_yyyyy_xxxxy_1[i] * fe_0 + ta_yyyyy_xxxxyz_0[i] * pa_z[i] - ta_yyyyy_xxxxyz_1[i] * pc_z[i];

        ta_yyyyyz_xxxxzz_0[i] = 4.0 * ta_yyyz_xxxxzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxxzz_1[i] * fe_0 + ta_yyyyz_xxxxzz_0[i] * pa_y[i] - ta_yyyyz_xxxxzz_1[i] * pc_y[i];

        ta_yyyyyz_xxxyyy_0[i] = ta_yyyyy_xxxyyy_0[i] * pa_z[i] - ta_yyyyy_xxxyyy_1[i] * pc_z[i];

        ta_yyyyyz_xxxyyz_0[i] = ta_yyyyy_xxxyy_0[i] * fe_0 - ta_yyyyy_xxxyy_1[i] * fe_0 + ta_yyyyy_xxxyyz_0[i] * pa_z[i] - ta_yyyyy_xxxyyz_1[i] * pc_z[i];

        ta_yyyyyz_xxxyzz_0[i] = 2.0 * ta_yyyyy_xxxyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xxxyz_1[i] * fe_0 + ta_yyyyy_xxxyzz_0[i] * pa_z[i] - ta_yyyyy_xxxyzz_1[i] * pc_z[i];

        ta_yyyyyz_xxxzzz_0[i] = 4.0 * ta_yyyz_xxxzzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxzzz_1[i] * fe_0 + ta_yyyyz_xxxzzz_0[i] * pa_y[i] - ta_yyyyz_xxxzzz_1[i] * pc_y[i];

        ta_yyyyyz_xxyyyy_0[i] = ta_yyyyy_xxyyyy_0[i] * pa_z[i] - ta_yyyyy_xxyyyy_1[i] * pc_z[i];

        ta_yyyyyz_xxyyyz_0[i] = ta_yyyyy_xxyyy_0[i] * fe_0 - ta_yyyyy_xxyyy_1[i] * fe_0 + ta_yyyyy_xxyyyz_0[i] * pa_z[i] - ta_yyyyy_xxyyyz_1[i] * pc_z[i];

        ta_yyyyyz_xxyyzz_0[i] = 2.0 * ta_yyyyy_xxyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xxyyz_1[i] * fe_0 + ta_yyyyy_xxyyzz_0[i] * pa_z[i] - ta_yyyyy_xxyyzz_1[i] * pc_z[i];

        ta_yyyyyz_xxyzzz_0[i] = 3.0 * ta_yyyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyzz_1[i] * fe_0 + ta_yyyyy_xxyzzz_0[i] * pa_z[i] - ta_yyyyy_xxyzzz_1[i] * pc_z[i];

        ta_yyyyyz_xxzzzz_0[i] = 4.0 * ta_yyyz_xxzzzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxzzzz_1[i] * fe_0 + ta_yyyyz_xxzzzz_0[i] * pa_y[i] - ta_yyyyz_xxzzzz_1[i] * pc_y[i];

        ta_yyyyyz_xyyyyy_0[i] = ta_yyyyy_xyyyyy_0[i] * pa_z[i] - ta_yyyyy_xyyyyy_1[i] * pc_z[i];

        ta_yyyyyz_xyyyyz_0[i] = ta_yyyyy_xyyyy_0[i] * fe_0 - ta_yyyyy_xyyyy_1[i] * fe_0 + ta_yyyyy_xyyyyz_0[i] * pa_z[i] - ta_yyyyy_xyyyyz_1[i] * pc_z[i];

        ta_yyyyyz_xyyyzz_0[i] = 2.0 * ta_yyyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyyz_1[i] * fe_0 + ta_yyyyy_xyyyzz_0[i] * pa_z[i] - ta_yyyyy_xyyyzz_1[i] * pc_z[i];

        ta_yyyyyz_xyyzzz_0[i] = 3.0 * ta_yyyyy_xyyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xyyzz_1[i] * fe_0 + ta_yyyyy_xyyzzz_0[i] * pa_z[i] - ta_yyyyy_xyyzzz_1[i] * pc_z[i];

        ta_yyyyyz_xyzzzz_0[i] = 4.0 * ta_yyyyy_xyzzz_0[i] * fe_0 - 4.0 * ta_yyyyy_xyzzz_1[i] * fe_0 + ta_yyyyy_xyzzzz_0[i] * pa_z[i] - ta_yyyyy_xyzzzz_1[i] * pc_z[i];

        ta_yyyyyz_xzzzzz_0[i] = 4.0 * ta_yyyz_xzzzzz_0[i] * fe_0 - 4.0 * ta_yyyz_xzzzzz_1[i] * fe_0 + ta_yyyyz_xzzzzz_0[i] * pa_y[i] - ta_yyyyz_xzzzzz_1[i] * pc_y[i];

        ta_yyyyyz_yyyyyy_0[i] = ta_yyyyy_yyyyyy_0[i] * pa_z[i] - ta_yyyyy_yyyyyy_1[i] * pc_z[i];

        ta_yyyyyz_yyyyyz_0[i] = ta_yyyyy_yyyyy_0[i] * fe_0 - ta_yyyyy_yyyyy_1[i] * fe_0 + ta_yyyyy_yyyyyz_0[i] * pa_z[i] - ta_yyyyy_yyyyyz_1[i] * pc_z[i];

        ta_yyyyyz_yyyyzz_0[i] = 2.0 * ta_yyyyy_yyyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_yyyyz_1[i] * fe_0 + ta_yyyyy_yyyyzz_0[i] * pa_z[i] - ta_yyyyy_yyyyzz_1[i] * pc_z[i];

        ta_yyyyyz_yyyzzz_0[i] = 3.0 * ta_yyyyy_yyyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_yyyzz_1[i] * fe_0 + ta_yyyyy_yyyzzz_0[i] * pa_z[i] - ta_yyyyy_yyyzzz_1[i] * pc_z[i];

        ta_yyyyyz_yyzzzz_0[i] = 4.0 * ta_yyyyy_yyzzz_0[i] * fe_0 - 4.0 * ta_yyyyy_yyzzz_1[i] * fe_0 + ta_yyyyy_yyzzzz_0[i] * pa_z[i] - ta_yyyyy_yyzzzz_1[i] * pc_z[i];

        ta_yyyyyz_yzzzzz_0[i] = 5.0 * ta_yyyyy_yzzzz_0[i] * fe_0 - 5.0 * ta_yyyyy_yzzzz_1[i] * fe_0 + ta_yyyyy_yzzzzz_0[i] * pa_z[i] - ta_yyyyy_yzzzzz_1[i] * pc_z[i];

        ta_yyyyyz_zzzzzz_0[i] = 4.0 * ta_yyyz_zzzzzz_0[i] * fe_0 - 4.0 * ta_yyyz_zzzzzz_1[i] * fe_0 + ta_yyyyz_zzzzzz_0[i] * pa_y[i] - ta_yyyyz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 644-672 components of targeted buffer : II

    auto ta_yyyyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 644);

    auto ta_yyyyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 645);

    auto ta_yyyyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 646);

    auto ta_yyyyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 647);

    auto ta_yyyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 648);

    auto ta_yyyyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 649);

    auto ta_yyyyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 650);

    auto ta_yyyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 651);

    auto ta_yyyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 652);

    auto ta_yyyyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 653);

    auto ta_yyyyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 654);

    auto ta_yyyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 655);

    auto ta_yyyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 656);

    auto ta_yyyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 657);

    auto ta_yyyyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 658);

    auto ta_yyyyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 659);

    auto ta_yyyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 660);

    auto ta_yyyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 661);

    auto ta_yyyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 662);

    auto ta_yyyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 663);

    auto ta_yyyyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 664);

    auto ta_yyyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 665);

    auto ta_yyyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 666);

    auto ta_yyyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 667);

    auto ta_yyyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 668);

    auto ta_yyyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 669);

    auto ta_yyyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 670);

    auto ta_yyyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 671);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yyyy_xxxxxy_0, ta_yyyy_xxxxxy_1, ta_yyyy_xxxxyy_0, ta_yyyy_xxxxyy_1, ta_yyyy_xxxyyy_0, ta_yyyy_xxxyyy_1, ta_yyyy_xxyyyy_0, ta_yyyy_xxyyyy_1, ta_yyyy_xyyyyy_0, ta_yyyy_xyyyyy_1, ta_yyyy_yyyyyy_0, ta_yyyy_yyyyyy_1, ta_yyyyz_xxxxxy_0, ta_yyyyz_xxxxxy_1, ta_yyyyz_xxxxyy_0, ta_yyyyz_xxxxyy_1, ta_yyyyz_xxxyyy_0, ta_yyyyz_xxxyyy_1, ta_yyyyz_xxyyyy_0, ta_yyyyz_xxyyyy_1, ta_yyyyz_xyyyyy_0, ta_yyyyz_xyyyyy_1, ta_yyyyz_yyyyyy_0, ta_yyyyz_yyyyyy_1, ta_yyyyzz_xxxxxx_0, ta_yyyyzz_xxxxxy_0, ta_yyyyzz_xxxxxz_0, ta_yyyyzz_xxxxyy_0, ta_yyyyzz_xxxxyz_0, ta_yyyyzz_xxxxzz_0, ta_yyyyzz_xxxyyy_0, ta_yyyyzz_xxxyyz_0, ta_yyyyzz_xxxyzz_0, ta_yyyyzz_xxxzzz_0, ta_yyyyzz_xxyyyy_0, ta_yyyyzz_xxyyyz_0, ta_yyyyzz_xxyyzz_0, ta_yyyyzz_xxyzzz_0, ta_yyyyzz_xxzzzz_0, ta_yyyyzz_xyyyyy_0, ta_yyyyzz_xyyyyz_0, ta_yyyyzz_xyyyzz_0, ta_yyyyzz_xyyzzz_0, ta_yyyyzz_xyzzzz_0, ta_yyyyzz_xzzzzz_0, ta_yyyyzz_yyyyyy_0, ta_yyyyzz_yyyyyz_0, ta_yyyyzz_yyyyzz_0, ta_yyyyzz_yyyzzz_0, ta_yyyyzz_yyzzzz_0, ta_yyyyzz_yzzzzz_0, ta_yyyyzz_zzzzzz_0, ta_yyyzz_xxxxxx_0, ta_yyyzz_xxxxxx_1, ta_yyyzz_xxxxxz_0, ta_yyyzz_xxxxxz_1, ta_yyyzz_xxxxyz_0, ta_yyyzz_xxxxyz_1, ta_yyyzz_xxxxz_0, ta_yyyzz_xxxxz_1, ta_yyyzz_xxxxzz_0, ta_yyyzz_xxxxzz_1, ta_yyyzz_xxxyyz_0, ta_yyyzz_xxxyyz_1, ta_yyyzz_xxxyz_0, ta_yyyzz_xxxyz_1, ta_yyyzz_xxxyzz_0, ta_yyyzz_xxxyzz_1, ta_yyyzz_xxxzz_0, ta_yyyzz_xxxzz_1, ta_yyyzz_xxxzzz_0, ta_yyyzz_xxxzzz_1, ta_yyyzz_xxyyyz_0, ta_yyyzz_xxyyyz_1, ta_yyyzz_xxyyz_0, ta_yyyzz_xxyyz_1, ta_yyyzz_xxyyzz_0, ta_yyyzz_xxyyzz_1, ta_yyyzz_xxyzz_0, ta_yyyzz_xxyzz_1, ta_yyyzz_xxyzzz_0, ta_yyyzz_xxyzzz_1, ta_yyyzz_xxzzz_0, ta_yyyzz_xxzzz_1, ta_yyyzz_xxzzzz_0, ta_yyyzz_xxzzzz_1, ta_yyyzz_xyyyyz_0, ta_yyyzz_xyyyyz_1, ta_yyyzz_xyyyz_0, ta_yyyzz_xyyyz_1, ta_yyyzz_xyyyzz_0, ta_yyyzz_xyyyzz_1, ta_yyyzz_xyyzz_0, ta_yyyzz_xyyzz_1, ta_yyyzz_xyyzzz_0, ta_yyyzz_xyyzzz_1, ta_yyyzz_xyzzz_0, ta_yyyzz_xyzzz_1, ta_yyyzz_xyzzzz_0, ta_yyyzz_xyzzzz_1, ta_yyyzz_xzzzz_0, ta_yyyzz_xzzzz_1, ta_yyyzz_xzzzzz_0, ta_yyyzz_xzzzzz_1, ta_yyyzz_yyyyyz_0, ta_yyyzz_yyyyyz_1, ta_yyyzz_yyyyz_0, ta_yyyzz_yyyyz_1, ta_yyyzz_yyyyzz_0, ta_yyyzz_yyyyzz_1, ta_yyyzz_yyyzz_0, ta_yyyzz_yyyzz_1, ta_yyyzz_yyyzzz_0, ta_yyyzz_yyyzzz_1, ta_yyyzz_yyzzz_0, ta_yyyzz_yyzzz_1, ta_yyyzz_yyzzzz_0, ta_yyyzz_yyzzzz_1, ta_yyyzz_yzzzz_0, ta_yyyzz_yzzzz_1, ta_yyyzz_yzzzzz_0, ta_yyyzz_yzzzzz_1, ta_yyyzz_zzzzz_0, ta_yyyzz_zzzzz_1, ta_yyyzz_zzzzzz_0, ta_yyyzz_zzzzzz_1, ta_yyzz_xxxxxx_0, ta_yyzz_xxxxxx_1, ta_yyzz_xxxxxz_0, ta_yyzz_xxxxxz_1, ta_yyzz_xxxxyz_0, ta_yyzz_xxxxyz_1, ta_yyzz_xxxxzz_0, ta_yyzz_xxxxzz_1, ta_yyzz_xxxyyz_0, ta_yyzz_xxxyyz_1, ta_yyzz_xxxyzz_0, ta_yyzz_xxxyzz_1, ta_yyzz_xxxzzz_0, ta_yyzz_xxxzzz_1, ta_yyzz_xxyyyz_0, ta_yyzz_xxyyyz_1, ta_yyzz_xxyyzz_0, ta_yyzz_xxyyzz_1, ta_yyzz_xxyzzz_0, ta_yyzz_xxyzzz_1, ta_yyzz_xxzzzz_0, ta_yyzz_xxzzzz_1, ta_yyzz_xyyyyz_0, ta_yyzz_xyyyyz_1, ta_yyzz_xyyyzz_0, ta_yyzz_xyyyzz_1, ta_yyzz_xyyzzz_0, ta_yyzz_xyyzzz_1, ta_yyzz_xyzzzz_0, ta_yyzz_xyzzzz_1, ta_yyzz_xzzzzz_0, ta_yyzz_xzzzzz_1, ta_yyzz_yyyyyz_0, ta_yyzz_yyyyyz_1, ta_yyzz_yyyyzz_0, ta_yyzz_yyyyzz_1, ta_yyzz_yyyzzz_0, ta_yyzz_yyyzzz_1, ta_yyzz_yyzzzz_0, ta_yyzz_yyzzzz_1, ta_yyzz_yzzzzz_0, ta_yyzz_yzzzzz_1, ta_yyzz_zzzzzz_0, ta_yyzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyzz_xxxxxx_0[i] = 3.0 * ta_yyzz_xxxxxx_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxxx_1[i] * fe_0 + ta_yyyzz_xxxxxx_0[i] * pa_y[i] - ta_yyyzz_xxxxxx_1[i] * pc_y[i];

        ta_yyyyzz_xxxxxy_0[i] = ta_yyyy_xxxxxy_0[i] * fe_0 - ta_yyyy_xxxxxy_1[i] * fe_0 + ta_yyyyz_xxxxxy_0[i] * pa_z[i] - ta_yyyyz_xxxxxy_1[i] * pc_z[i];

        ta_yyyyzz_xxxxxz_0[i] = 3.0 * ta_yyzz_xxxxxz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxxz_1[i] * fe_0 + ta_yyyzz_xxxxxz_0[i] * pa_y[i] - ta_yyyzz_xxxxxz_1[i] * pc_y[i];

        ta_yyyyzz_xxxxyy_0[i] = ta_yyyy_xxxxyy_0[i] * fe_0 - ta_yyyy_xxxxyy_1[i] * fe_0 + ta_yyyyz_xxxxyy_0[i] * pa_z[i] - ta_yyyyz_xxxxyy_1[i] * pc_z[i];

        ta_yyyyzz_xxxxyz_0[i] = 3.0 * ta_yyzz_xxxxyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxyz_1[i] * fe_0 + ta_yyyzz_xxxxz_0[i] * fe_0 - ta_yyyzz_xxxxz_1[i] * fe_0 + ta_yyyzz_xxxxyz_0[i] * pa_y[i] - ta_yyyzz_xxxxyz_1[i] * pc_y[i];

        ta_yyyyzz_xxxxzz_0[i] = 3.0 * ta_yyzz_xxxxzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxzz_1[i] * fe_0 + ta_yyyzz_xxxxzz_0[i] * pa_y[i] - ta_yyyzz_xxxxzz_1[i] * pc_y[i];

        ta_yyyyzz_xxxyyy_0[i] = ta_yyyy_xxxyyy_0[i] * fe_0 - ta_yyyy_xxxyyy_1[i] * fe_0 + ta_yyyyz_xxxyyy_0[i] * pa_z[i] - ta_yyyyz_xxxyyy_1[i] * pc_z[i];

        ta_yyyyzz_xxxyyz_0[i] = 3.0 * ta_yyzz_xxxyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxyyz_1[i] * fe_0 + 2.0 * ta_yyyzz_xxxyz_0[i] * fe_0 - 2.0 * ta_yyyzz_xxxyz_1[i] * fe_0 + ta_yyyzz_xxxyyz_0[i] * pa_y[i] - ta_yyyzz_xxxyyz_1[i] * pc_y[i];

        ta_yyyyzz_xxxyzz_0[i] = 3.0 * ta_yyzz_xxxyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxyzz_1[i] * fe_0 + ta_yyyzz_xxxzz_0[i] * fe_0 - ta_yyyzz_xxxzz_1[i] * fe_0 + ta_yyyzz_xxxyzz_0[i] * pa_y[i] - ta_yyyzz_xxxyzz_1[i] * pc_y[i];

        ta_yyyyzz_xxxzzz_0[i] = 3.0 * ta_yyzz_xxxzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxzzz_1[i] * fe_0 + ta_yyyzz_xxxzzz_0[i] * pa_y[i] - ta_yyyzz_xxxzzz_1[i] * pc_y[i];

        ta_yyyyzz_xxyyyy_0[i] = ta_yyyy_xxyyyy_0[i] * fe_0 - ta_yyyy_xxyyyy_1[i] * fe_0 + ta_yyyyz_xxyyyy_0[i] * pa_z[i] - ta_yyyyz_xxyyyy_1[i] * pc_z[i];

        ta_yyyyzz_xxyyyz_0[i] = 3.0 * ta_yyzz_xxyyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyyz_1[i] * fe_0 + 3.0 * ta_yyyzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxyyz_1[i] * fe_0 + ta_yyyzz_xxyyyz_0[i] * pa_y[i] - ta_yyyzz_xxyyyz_1[i] * pc_y[i];

        ta_yyyyzz_xxyyzz_0[i] = 3.0 * ta_yyzz_xxyyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_yyyzz_xxyzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xxyzz_1[i] * fe_0 + ta_yyyzz_xxyyzz_0[i] * pa_y[i] - ta_yyyzz_xxyyzz_1[i] * pc_y[i];

        ta_yyyyzz_xxyzzz_0[i] = 3.0 * ta_yyzz_xxyzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyzzz_1[i] * fe_0 + ta_yyyzz_xxzzz_0[i] * fe_0 - ta_yyyzz_xxzzz_1[i] * fe_0 + ta_yyyzz_xxyzzz_0[i] * pa_y[i] - ta_yyyzz_xxyzzz_1[i] * pc_y[i];

        ta_yyyyzz_xxzzzz_0[i] = 3.0 * ta_yyzz_xxzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxzzzz_1[i] * fe_0 + ta_yyyzz_xxzzzz_0[i] * pa_y[i] - ta_yyyzz_xxzzzz_1[i] * pc_y[i];

        ta_yyyyzz_xyyyyy_0[i] = ta_yyyy_xyyyyy_0[i] * fe_0 - ta_yyyy_xyyyyy_1[i] * fe_0 + ta_yyyyz_xyyyyy_0[i] * pa_z[i] - ta_yyyyz_xyyyyy_1[i] * pc_z[i];

        ta_yyyyzz_xyyyyz_0[i] = 3.0 * ta_yyzz_xyyyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xyyyyz_1[i] * fe_0 + 4.0 * ta_yyyzz_xyyyz_0[i] * fe_0 - 4.0 * ta_yyyzz_xyyyz_1[i] * fe_0 + ta_yyyzz_xyyyyz_0[i] * pa_y[i] - ta_yyyzz_xyyyyz_1[i] * pc_y[i];

        ta_yyyyzz_xyyyzz_0[i] = 3.0 * ta_yyzz_xyyyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xyyyzz_1[i] * fe_0 + 3.0 * ta_yyyzz_xyyzz_0[i] * fe_0 - 3.0 * ta_yyyzz_xyyzz_1[i] * fe_0 + ta_yyyzz_xyyyzz_0[i] * pa_y[i] - ta_yyyzz_xyyyzz_1[i] * pc_y[i];

        ta_yyyyzz_xyyzzz_0[i] = 3.0 * ta_yyzz_xyyzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xyyzzz_1[i] * fe_0 + 2.0 * ta_yyyzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyzzz_1[i] * fe_0 + ta_yyyzz_xyyzzz_0[i] * pa_y[i] - ta_yyyzz_xyyzzz_1[i] * pc_y[i];

        ta_yyyyzz_xyzzzz_0[i] = 3.0 * ta_yyzz_xyzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xyzzzz_1[i] * fe_0 + ta_yyyzz_xzzzz_0[i] * fe_0 - ta_yyyzz_xzzzz_1[i] * fe_0 + ta_yyyzz_xyzzzz_0[i] * pa_y[i] - ta_yyyzz_xyzzzz_1[i] * pc_y[i];

        ta_yyyyzz_xzzzzz_0[i] = 3.0 * ta_yyzz_xzzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xzzzzz_1[i] * fe_0 + ta_yyyzz_xzzzzz_0[i] * pa_y[i] - ta_yyyzz_xzzzzz_1[i] * pc_y[i];

        ta_yyyyzz_yyyyyy_0[i] = ta_yyyy_yyyyyy_0[i] * fe_0 - ta_yyyy_yyyyyy_1[i] * fe_0 + ta_yyyyz_yyyyyy_0[i] * pa_z[i] - ta_yyyyz_yyyyyy_1[i] * pc_z[i];

        ta_yyyyzz_yyyyyz_0[i] = 3.0 * ta_yyzz_yyyyyz_0[i] * fe_0 - 3.0 * ta_yyzz_yyyyyz_1[i] * fe_0 + 5.0 * ta_yyyzz_yyyyz_0[i] * fe_0 - 5.0 * ta_yyyzz_yyyyz_1[i] * fe_0 + ta_yyyzz_yyyyyz_0[i] * pa_y[i] - ta_yyyzz_yyyyyz_1[i] * pc_y[i];

        ta_yyyyzz_yyyyzz_0[i] = 3.0 * ta_yyzz_yyyyzz_0[i] * fe_0 - 3.0 * ta_yyzz_yyyyzz_1[i] * fe_0 + 4.0 * ta_yyyzz_yyyzz_0[i] * fe_0 - 4.0 * ta_yyyzz_yyyzz_1[i] * fe_0 + ta_yyyzz_yyyyzz_0[i] * pa_y[i] - ta_yyyzz_yyyyzz_1[i] * pc_y[i];

        ta_yyyyzz_yyyzzz_0[i] = 3.0 * ta_yyzz_yyyzzz_0[i] * fe_0 - 3.0 * ta_yyzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_yyyzz_yyzzz_0[i] * fe_0 - 3.0 * ta_yyyzz_yyzzz_1[i] * fe_0 + ta_yyyzz_yyyzzz_0[i] * pa_y[i] - ta_yyyzz_yyyzzz_1[i] * pc_y[i];

        ta_yyyyzz_yyzzzz_0[i] = 3.0 * ta_yyzz_yyzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_yyzzzz_1[i] * fe_0 + 2.0 * ta_yyyzz_yzzzz_0[i] * fe_0 - 2.0 * ta_yyyzz_yzzzz_1[i] * fe_0 + ta_yyyzz_yyzzzz_0[i] * pa_y[i] - ta_yyyzz_yyzzzz_1[i] * pc_y[i];

        ta_yyyyzz_yzzzzz_0[i] = 3.0 * ta_yyzz_yzzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_yzzzzz_1[i] * fe_0 + ta_yyyzz_zzzzz_0[i] * fe_0 - ta_yyyzz_zzzzz_1[i] * fe_0 + ta_yyyzz_yzzzzz_0[i] * pa_y[i] - ta_yyyzz_yzzzzz_1[i] * pc_y[i];

        ta_yyyyzz_zzzzzz_0[i] = 3.0 * ta_yyzz_zzzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_zzzzzz_1[i] * fe_0 + ta_yyyzz_zzzzzz_0[i] * pa_y[i] - ta_yyyzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 672-700 components of targeted buffer : II

    auto ta_yyyzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 672);

    auto ta_yyyzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 673);

    auto ta_yyyzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 674);

    auto ta_yyyzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 675);

    auto ta_yyyzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 676);

    auto ta_yyyzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 677);

    auto ta_yyyzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 678);

    auto ta_yyyzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 679);

    auto ta_yyyzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 680);

    auto ta_yyyzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 681);

    auto ta_yyyzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 682);

    auto ta_yyyzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 683);

    auto ta_yyyzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 684);

    auto ta_yyyzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 685);

    auto ta_yyyzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 686);

    auto ta_yyyzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 687);

    auto ta_yyyzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 688);

    auto ta_yyyzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 689);

    auto ta_yyyzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 690);

    auto ta_yyyzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 691);

    auto ta_yyyzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 692);

    auto ta_yyyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 693);

    auto ta_yyyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 694);

    auto ta_yyyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 695);

    auto ta_yyyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 696);

    auto ta_yyyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 697);

    auto ta_yyyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 698);

    auto ta_yyyzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 699);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yyyz_xxxxxy_0, ta_yyyz_xxxxxy_1, ta_yyyz_xxxxyy_0, ta_yyyz_xxxxyy_1, ta_yyyz_xxxyyy_0, ta_yyyz_xxxyyy_1, ta_yyyz_xxyyyy_0, ta_yyyz_xxyyyy_1, ta_yyyz_xyyyyy_0, ta_yyyz_xyyyyy_1, ta_yyyz_yyyyyy_0, ta_yyyz_yyyyyy_1, ta_yyyzz_xxxxxy_0, ta_yyyzz_xxxxxy_1, ta_yyyzz_xxxxyy_0, ta_yyyzz_xxxxyy_1, ta_yyyzz_xxxyyy_0, ta_yyyzz_xxxyyy_1, ta_yyyzz_xxyyyy_0, ta_yyyzz_xxyyyy_1, ta_yyyzz_xyyyyy_0, ta_yyyzz_xyyyyy_1, ta_yyyzz_yyyyyy_0, ta_yyyzz_yyyyyy_1, ta_yyyzzz_xxxxxx_0, ta_yyyzzz_xxxxxy_0, ta_yyyzzz_xxxxxz_0, ta_yyyzzz_xxxxyy_0, ta_yyyzzz_xxxxyz_0, ta_yyyzzz_xxxxzz_0, ta_yyyzzz_xxxyyy_0, ta_yyyzzz_xxxyyz_0, ta_yyyzzz_xxxyzz_0, ta_yyyzzz_xxxzzz_0, ta_yyyzzz_xxyyyy_0, ta_yyyzzz_xxyyyz_0, ta_yyyzzz_xxyyzz_0, ta_yyyzzz_xxyzzz_0, ta_yyyzzz_xxzzzz_0, ta_yyyzzz_xyyyyy_0, ta_yyyzzz_xyyyyz_0, ta_yyyzzz_xyyyzz_0, ta_yyyzzz_xyyzzz_0, ta_yyyzzz_xyzzzz_0, ta_yyyzzz_xzzzzz_0, ta_yyyzzz_yyyyyy_0, ta_yyyzzz_yyyyyz_0, ta_yyyzzz_yyyyzz_0, ta_yyyzzz_yyyzzz_0, ta_yyyzzz_yyzzzz_0, ta_yyyzzz_yzzzzz_0, ta_yyyzzz_zzzzzz_0, ta_yyzzz_xxxxxx_0, ta_yyzzz_xxxxxx_1, ta_yyzzz_xxxxxz_0, ta_yyzzz_xxxxxz_1, ta_yyzzz_xxxxyz_0, ta_yyzzz_xxxxyz_1, ta_yyzzz_xxxxz_0, ta_yyzzz_xxxxz_1, ta_yyzzz_xxxxzz_0, ta_yyzzz_xxxxzz_1, ta_yyzzz_xxxyyz_0, ta_yyzzz_xxxyyz_1, ta_yyzzz_xxxyz_0, ta_yyzzz_xxxyz_1, ta_yyzzz_xxxyzz_0, ta_yyzzz_xxxyzz_1, ta_yyzzz_xxxzz_0, ta_yyzzz_xxxzz_1, ta_yyzzz_xxxzzz_0, ta_yyzzz_xxxzzz_1, ta_yyzzz_xxyyyz_0, ta_yyzzz_xxyyyz_1, ta_yyzzz_xxyyz_0, ta_yyzzz_xxyyz_1, ta_yyzzz_xxyyzz_0, ta_yyzzz_xxyyzz_1, ta_yyzzz_xxyzz_0, ta_yyzzz_xxyzz_1, ta_yyzzz_xxyzzz_0, ta_yyzzz_xxyzzz_1, ta_yyzzz_xxzzz_0, ta_yyzzz_xxzzz_1, ta_yyzzz_xxzzzz_0, ta_yyzzz_xxzzzz_1, ta_yyzzz_xyyyyz_0, ta_yyzzz_xyyyyz_1, ta_yyzzz_xyyyz_0, ta_yyzzz_xyyyz_1, ta_yyzzz_xyyyzz_0, ta_yyzzz_xyyyzz_1, ta_yyzzz_xyyzz_0, ta_yyzzz_xyyzz_1, ta_yyzzz_xyyzzz_0, ta_yyzzz_xyyzzz_1, ta_yyzzz_xyzzz_0, ta_yyzzz_xyzzz_1, ta_yyzzz_xyzzzz_0, ta_yyzzz_xyzzzz_1, ta_yyzzz_xzzzz_0, ta_yyzzz_xzzzz_1, ta_yyzzz_xzzzzz_0, ta_yyzzz_xzzzzz_1, ta_yyzzz_yyyyyz_0, ta_yyzzz_yyyyyz_1, ta_yyzzz_yyyyz_0, ta_yyzzz_yyyyz_1, ta_yyzzz_yyyyzz_0, ta_yyzzz_yyyyzz_1, ta_yyzzz_yyyzz_0, ta_yyzzz_yyyzz_1, ta_yyzzz_yyyzzz_0, ta_yyzzz_yyyzzz_1, ta_yyzzz_yyzzz_0, ta_yyzzz_yyzzz_1, ta_yyzzz_yyzzzz_0, ta_yyzzz_yyzzzz_1, ta_yyzzz_yzzzz_0, ta_yyzzz_yzzzz_1, ta_yyzzz_yzzzzz_0, ta_yyzzz_yzzzzz_1, ta_yyzzz_zzzzz_0, ta_yyzzz_zzzzz_1, ta_yyzzz_zzzzzz_0, ta_yyzzz_zzzzzz_1, ta_yzzz_xxxxxx_0, ta_yzzz_xxxxxx_1, ta_yzzz_xxxxxz_0, ta_yzzz_xxxxxz_1, ta_yzzz_xxxxyz_0, ta_yzzz_xxxxyz_1, ta_yzzz_xxxxzz_0, ta_yzzz_xxxxzz_1, ta_yzzz_xxxyyz_0, ta_yzzz_xxxyyz_1, ta_yzzz_xxxyzz_0, ta_yzzz_xxxyzz_1, ta_yzzz_xxxzzz_0, ta_yzzz_xxxzzz_1, ta_yzzz_xxyyyz_0, ta_yzzz_xxyyyz_1, ta_yzzz_xxyyzz_0, ta_yzzz_xxyyzz_1, ta_yzzz_xxyzzz_0, ta_yzzz_xxyzzz_1, ta_yzzz_xxzzzz_0, ta_yzzz_xxzzzz_1, ta_yzzz_xyyyyz_0, ta_yzzz_xyyyyz_1, ta_yzzz_xyyyzz_0, ta_yzzz_xyyyzz_1, ta_yzzz_xyyzzz_0, ta_yzzz_xyyzzz_1, ta_yzzz_xyzzzz_0, ta_yzzz_xyzzzz_1, ta_yzzz_xzzzzz_0, ta_yzzz_xzzzzz_1, ta_yzzz_yyyyyz_0, ta_yzzz_yyyyyz_1, ta_yzzz_yyyyzz_0, ta_yzzz_yyyyzz_1, ta_yzzz_yyyzzz_0, ta_yzzz_yyyzzz_1, ta_yzzz_yyzzzz_0, ta_yzzz_yyzzzz_1, ta_yzzz_yzzzzz_0, ta_yzzz_yzzzzz_1, ta_yzzz_zzzzzz_0, ta_yzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzzz_xxxxxx_0[i] = 2.0 * ta_yzzz_xxxxxx_0[i] * fe_0 - 2.0 * ta_yzzz_xxxxxx_1[i] * fe_0 + ta_yyzzz_xxxxxx_0[i] * pa_y[i] - ta_yyzzz_xxxxxx_1[i] * pc_y[i];

        ta_yyyzzz_xxxxxy_0[i] = 2.0 * ta_yyyz_xxxxxy_0[i] * fe_0 - 2.0 * ta_yyyz_xxxxxy_1[i] * fe_0 + ta_yyyzz_xxxxxy_0[i] * pa_z[i] - ta_yyyzz_xxxxxy_1[i] * pc_z[i];

        ta_yyyzzz_xxxxxz_0[i] = 2.0 * ta_yzzz_xxxxxz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxxxz_1[i] * fe_0 + ta_yyzzz_xxxxxz_0[i] * pa_y[i] - ta_yyzzz_xxxxxz_1[i] * pc_y[i];

        ta_yyyzzz_xxxxyy_0[i] = 2.0 * ta_yyyz_xxxxyy_0[i] * fe_0 - 2.0 * ta_yyyz_xxxxyy_1[i] * fe_0 + ta_yyyzz_xxxxyy_0[i] * pa_z[i] - ta_yyyzz_xxxxyy_1[i] * pc_z[i];

        ta_yyyzzz_xxxxyz_0[i] = 2.0 * ta_yzzz_xxxxyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxxyz_1[i] * fe_0 + ta_yyzzz_xxxxz_0[i] * fe_0 - ta_yyzzz_xxxxz_1[i] * fe_0 + ta_yyzzz_xxxxyz_0[i] * pa_y[i] - ta_yyzzz_xxxxyz_1[i] * pc_y[i];

        ta_yyyzzz_xxxxzz_0[i] = 2.0 * ta_yzzz_xxxxzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxxzz_1[i] * fe_0 + ta_yyzzz_xxxxzz_0[i] * pa_y[i] - ta_yyzzz_xxxxzz_1[i] * pc_y[i];

        ta_yyyzzz_xxxyyy_0[i] = 2.0 * ta_yyyz_xxxyyy_0[i] * fe_0 - 2.0 * ta_yyyz_xxxyyy_1[i] * fe_0 + ta_yyyzz_xxxyyy_0[i] * pa_z[i] - ta_yyyzz_xxxyyy_1[i] * pc_z[i];

        ta_yyyzzz_xxxyyz_0[i] = 2.0 * ta_yzzz_xxxyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxyyz_1[i] * fe_0 + 2.0 * ta_yyzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_yyzzz_xxxyz_1[i] * fe_0 + ta_yyzzz_xxxyyz_0[i] * pa_y[i] - ta_yyzzz_xxxyyz_1[i] * pc_y[i];

        ta_yyyzzz_xxxyzz_0[i] = 2.0 * ta_yzzz_xxxyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxyzz_1[i] * fe_0 + ta_yyzzz_xxxzz_0[i] * fe_0 - ta_yyzzz_xxxzz_1[i] * fe_0 + ta_yyzzz_xxxyzz_0[i] * pa_y[i] - ta_yyzzz_xxxyzz_1[i] * pc_y[i];

        ta_yyyzzz_xxxzzz_0[i] = 2.0 * ta_yzzz_xxxzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxzzz_1[i] * fe_0 + ta_yyzzz_xxxzzz_0[i] * pa_y[i] - ta_yyzzz_xxxzzz_1[i] * pc_y[i];

        ta_yyyzzz_xxyyyy_0[i] = 2.0 * ta_yyyz_xxyyyy_0[i] * fe_0 - 2.0 * ta_yyyz_xxyyyy_1[i] * fe_0 + ta_yyyzz_xxyyyy_0[i] * pa_z[i] - ta_yyyzz_xxyyyy_1[i] * pc_z[i];

        ta_yyyzzz_xxyyyz_0[i] = 2.0 * ta_yzzz_xxyyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyyyz_1[i] * fe_0 + 3.0 * ta_yyzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxyyz_1[i] * fe_0 + ta_yyzzz_xxyyyz_0[i] * pa_y[i] - ta_yyzzz_xxyyyz_1[i] * pc_y[i];

        ta_yyyzzz_xxyyzz_0[i] = 2.0 * ta_yzzz_xxyyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_yyzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xxyzz_1[i] * fe_0 + ta_yyzzz_xxyyzz_0[i] * pa_y[i] - ta_yyzzz_xxyyzz_1[i] * pc_y[i];

        ta_yyyzzz_xxyzzz_0[i] = 2.0 * ta_yzzz_xxyzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyzzz_1[i] * fe_0 + ta_yyzzz_xxzzz_0[i] * fe_0 - ta_yyzzz_xxzzz_1[i] * fe_0 + ta_yyzzz_xxyzzz_0[i] * pa_y[i] - ta_yyzzz_xxyzzz_1[i] * pc_y[i];

        ta_yyyzzz_xxzzzz_0[i] = 2.0 * ta_yzzz_xxzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxzzzz_1[i] * fe_0 + ta_yyzzz_xxzzzz_0[i] * pa_y[i] - ta_yyzzz_xxzzzz_1[i] * pc_y[i];

        ta_yyyzzz_xyyyyy_0[i] = 2.0 * ta_yyyz_xyyyyy_0[i] * fe_0 - 2.0 * ta_yyyz_xyyyyy_1[i] * fe_0 + ta_yyyzz_xyyyyy_0[i] * pa_z[i] - ta_yyyzz_xyyyyy_1[i] * pc_z[i];

        ta_yyyzzz_xyyyyz_0[i] = 2.0 * ta_yzzz_xyyyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyyyz_1[i] * fe_0 + 4.0 * ta_yyzzz_xyyyz_0[i] * fe_0 - 4.0 * ta_yyzzz_xyyyz_1[i] * fe_0 + ta_yyzzz_xyyyyz_0[i] * pa_y[i] - ta_yyzzz_xyyyyz_1[i] * pc_y[i];

        ta_yyyzzz_xyyyzz_0[i] = 2.0 * ta_yzzz_xyyyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyyzz_1[i] * fe_0 + 3.0 * ta_yyzzz_xyyzz_0[i] * fe_0 - 3.0 * ta_yyzzz_xyyzz_1[i] * fe_0 + ta_yyzzz_xyyyzz_0[i] * pa_y[i] - ta_yyzzz_xyyyzz_1[i] * pc_y[i];

        ta_yyyzzz_xyyzzz_0[i] = 2.0 * ta_yzzz_xyyzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyzzz_1[i] * fe_0 + 2.0 * ta_yyzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyzzz_1[i] * fe_0 + ta_yyzzz_xyyzzz_0[i] * pa_y[i] - ta_yyzzz_xyyzzz_1[i] * pc_y[i];

        ta_yyyzzz_xyzzzz_0[i] = 2.0 * ta_yzzz_xyzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyzzzz_1[i] * fe_0 + ta_yyzzz_xzzzz_0[i] * fe_0 - ta_yyzzz_xzzzz_1[i] * fe_0 + ta_yyzzz_xyzzzz_0[i] * pa_y[i] - ta_yyzzz_xyzzzz_1[i] * pc_y[i];

        ta_yyyzzz_xzzzzz_0[i] = 2.0 * ta_yzzz_xzzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xzzzzz_1[i] * fe_0 + ta_yyzzz_xzzzzz_0[i] * pa_y[i] - ta_yyzzz_xzzzzz_1[i] * pc_y[i];

        ta_yyyzzz_yyyyyy_0[i] = 2.0 * ta_yyyz_yyyyyy_0[i] * fe_0 - 2.0 * ta_yyyz_yyyyyy_1[i] * fe_0 + ta_yyyzz_yyyyyy_0[i] * pa_z[i] - ta_yyyzz_yyyyyy_1[i] * pc_z[i];

        ta_yyyzzz_yyyyyz_0[i] = 2.0 * ta_yzzz_yyyyyz_0[i] * fe_0 - 2.0 * ta_yzzz_yyyyyz_1[i] * fe_0 + 5.0 * ta_yyzzz_yyyyz_0[i] * fe_0 - 5.0 * ta_yyzzz_yyyyz_1[i] * fe_0 + ta_yyzzz_yyyyyz_0[i] * pa_y[i] - ta_yyzzz_yyyyyz_1[i] * pc_y[i];

        ta_yyyzzz_yyyyzz_0[i] = 2.0 * ta_yzzz_yyyyzz_0[i] * fe_0 - 2.0 * ta_yzzz_yyyyzz_1[i] * fe_0 + 4.0 * ta_yyzzz_yyyzz_0[i] * fe_0 - 4.0 * ta_yyzzz_yyyzz_1[i] * fe_0 + ta_yyzzz_yyyyzz_0[i] * pa_y[i] - ta_yyzzz_yyyyzz_1[i] * pc_y[i];

        ta_yyyzzz_yyyzzz_0[i] = 2.0 * ta_yzzz_yyyzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_yyzzz_yyzzz_0[i] * fe_0 - 3.0 * ta_yyzzz_yyzzz_1[i] * fe_0 + ta_yyzzz_yyyzzz_0[i] * pa_y[i] - ta_yyzzz_yyyzzz_1[i] * pc_y[i];

        ta_yyyzzz_yyzzzz_0[i] = 2.0 * ta_yzzz_yyzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yyzzzz_1[i] * fe_0 + 2.0 * ta_yyzzz_yzzzz_0[i] * fe_0 - 2.0 * ta_yyzzz_yzzzz_1[i] * fe_0 + ta_yyzzz_yyzzzz_0[i] * pa_y[i] - ta_yyzzz_yyzzzz_1[i] * pc_y[i];

        ta_yyyzzz_yzzzzz_0[i] = 2.0 * ta_yzzz_yzzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yzzzzz_1[i] * fe_0 + ta_yyzzz_zzzzz_0[i] * fe_0 - ta_yyzzz_zzzzz_1[i] * fe_0 + ta_yyzzz_yzzzzz_0[i] * pa_y[i] - ta_yyzzz_yzzzzz_1[i] * pc_y[i];

        ta_yyyzzz_zzzzzz_0[i] = 2.0 * ta_yzzz_zzzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_zzzzzz_1[i] * fe_0 + ta_yyzzz_zzzzzz_0[i] * pa_y[i] - ta_yyzzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 700-728 components of targeted buffer : II

    auto ta_yyzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 700);

    auto ta_yyzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 701);

    auto ta_yyzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 702);

    auto ta_yyzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 703);

    auto ta_yyzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 704);

    auto ta_yyzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 705);

    auto ta_yyzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 706);

    auto ta_yyzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 707);

    auto ta_yyzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 708);

    auto ta_yyzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 709);

    auto ta_yyzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 710);

    auto ta_yyzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 711);

    auto ta_yyzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 712);

    auto ta_yyzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 713);

    auto ta_yyzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 714);

    auto ta_yyzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 715);

    auto ta_yyzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 716);

    auto ta_yyzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 717);

    auto ta_yyzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 718);

    auto ta_yyzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 719);

    auto ta_yyzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 720);

    auto ta_yyzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 721);

    auto ta_yyzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 722);

    auto ta_yyzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 723);

    auto ta_yyzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 724);

    auto ta_yyzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 725);

    auto ta_yyzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 726);

    auto ta_yyzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 727);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yyzz_xxxxxy_0, ta_yyzz_xxxxxy_1, ta_yyzz_xxxxyy_0, ta_yyzz_xxxxyy_1, ta_yyzz_xxxyyy_0, ta_yyzz_xxxyyy_1, ta_yyzz_xxyyyy_0, ta_yyzz_xxyyyy_1, ta_yyzz_xyyyyy_0, ta_yyzz_xyyyyy_1, ta_yyzz_yyyyyy_0, ta_yyzz_yyyyyy_1, ta_yyzzz_xxxxxy_0, ta_yyzzz_xxxxxy_1, ta_yyzzz_xxxxyy_0, ta_yyzzz_xxxxyy_1, ta_yyzzz_xxxyyy_0, ta_yyzzz_xxxyyy_1, ta_yyzzz_xxyyyy_0, ta_yyzzz_xxyyyy_1, ta_yyzzz_xyyyyy_0, ta_yyzzz_xyyyyy_1, ta_yyzzz_yyyyyy_0, ta_yyzzz_yyyyyy_1, ta_yyzzzz_xxxxxx_0, ta_yyzzzz_xxxxxy_0, ta_yyzzzz_xxxxxz_0, ta_yyzzzz_xxxxyy_0, ta_yyzzzz_xxxxyz_0, ta_yyzzzz_xxxxzz_0, ta_yyzzzz_xxxyyy_0, ta_yyzzzz_xxxyyz_0, ta_yyzzzz_xxxyzz_0, ta_yyzzzz_xxxzzz_0, ta_yyzzzz_xxyyyy_0, ta_yyzzzz_xxyyyz_0, ta_yyzzzz_xxyyzz_0, ta_yyzzzz_xxyzzz_0, ta_yyzzzz_xxzzzz_0, ta_yyzzzz_xyyyyy_0, ta_yyzzzz_xyyyyz_0, ta_yyzzzz_xyyyzz_0, ta_yyzzzz_xyyzzz_0, ta_yyzzzz_xyzzzz_0, ta_yyzzzz_xzzzzz_0, ta_yyzzzz_yyyyyy_0, ta_yyzzzz_yyyyyz_0, ta_yyzzzz_yyyyzz_0, ta_yyzzzz_yyyzzz_0, ta_yyzzzz_yyzzzz_0, ta_yyzzzz_yzzzzz_0, ta_yyzzzz_zzzzzz_0, ta_yzzzz_xxxxxx_0, ta_yzzzz_xxxxxx_1, ta_yzzzz_xxxxxz_0, ta_yzzzz_xxxxxz_1, ta_yzzzz_xxxxyz_0, ta_yzzzz_xxxxyz_1, ta_yzzzz_xxxxz_0, ta_yzzzz_xxxxz_1, ta_yzzzz_xxxxzz_0, ta_yzzzz_xxxxzz_1, ta_yzzzz_xxxyyz_0, ta_yzzzz_xxxyyz_1, ta_yzzzz_xxxyz_0, ta_yzzzz_xxxyz_1, ta_yzzzz_xxxyzz_0, ta_yzzzz_xxxyzz_1, ta_yzzzz_xxxzz_0, ta_yzzzz_xxxzz_1, ta_yzzzz_xxxzzz_0, ta_yzzzz_xxxzzz_1, ta_yzzzz_xxyyyz_0, ta_yzzzz_xxyyyz_1, ta_yzzzz_xxyyz_0, ta_yzzzz_xxyyz_1, ta_yzzzz_xxyyzz_0, ta_yzzzz_xxyyzz_1, ta_yzzzz_xxyzz_0, ta_yzzzz_xxyzz_1, ta_yzzzz_xxyzzz_0, ta_yzzzz_xxyzzz_1, ta_yzzzz_xxzzz_0, ta_yzzzz_xxzzz_1, ta_yzzzz_xxzzzz_0, ta_yzzzz_xxzzzz_1, ta_yzzzz_xyyyyz_0, ta_yzzzz_xyyyyz_1, ta_yzzzz_xyyyz_0, ta_yzzzz_xyyyz_1, ta_yzzzz_xyyyzz_0, ta_yzzzz_xyyyzz_1, ta_yzzzz_xyyzz_0, ta_yzzzz_xyyzz_1, ta_yzzzz_xyyzzz_0, ta_yzzzz_xyyzzz_1, ta_yzzzz_xyzzz_0, ta_yzzzz_xyzzz_1, ta_yzzzz_xyzzzz_0, ta_yzzzz_xyzzzz_1, ta_yzzzz_xzzzz_0, ta_yzzzz_xzzzz_1, ta_yzzzz_xzzzzz_0, ta_yzzzz_xzzzzz_1, ta_yzzzz_yyyyyz_0, ta_yzzzz_yyyyyz_1, ta_yzzzz_yyyyz_0, ta_yzzzz_yyyyz_1, ta_yzzzz_yyyyzz_0, ta_yzzzz_yyyyzz_1, ta_yzzzz_yyyzz_0, ta_yzzzz_yyyzz_1, ta_yzzzz_yyyzzz_0, ta_yzzzz_yyyzzz_1, ta_yzzzz_yyzzz_0, ta_yzzzz_yyzzz_1, ta_yzzzz_yyzzzz_0, ta_yzzzz_yyzzzz_1, ta_yzzzz_yzzzz_0, ta_yzzzz_yzzzz_1, ta_yzzzz_yzzzzz_0, ta_yzzzz_yzzzzz_1, ta_yzzzz_zzzzz_0, ta_yzzzz_zzzzz_1, ta_yzzzz_zzzzzz_0, ta_yzzzz_zzzzzz_1, ta_zzzz_xxxxxx_0, ta_zzzz_xxxxxx_1, ta_zzzz_xxxxxz_0, ta_zzzz_xxxxxz_1, ta_zzzz_xxxxyz_0, ta_zzzz_xxxxyz_1, ta_zzzz_xxxxzz_0, ta_zzzz_xxxxzz_1, ta_zzzz_xxxyyz_0, ta_zzzz_xxxyyz_1, ta_zzzz_xxxyzz_0, ta_zzzz_xxxyzz_1, ta_zzzz_xxxzzz_0, ta_zzzz_xxxzzz_1, ta_zzzz_xxyyyz_0, ta_zzzz_xxyyyz_1, ta_zzzz_xxyyzz_0, ta_zzzz_xxyyzz_1, ta_zzzz_xxyzzz_0, ta_zzzz_xxyzzz_1, ta_zzzz_xxzzzz_0, ta_zzzz_xxzzzz_1, ta_zzzz_xyyyyz_0, ta_zzzz_xyyyyz_1, ta_zzzz_xyyyzz_0, ta_zzzz_xyyyzz_1, ta_zzzz_xyyzzz_0, ta_zzzz_xyyzzz_1, ta_zzzz_xyzzzz_0, ta_zzzz_xyzzzz_1, ta_zzzz_xzzzzz_0, ta_zzzz_xzzzzz_1, ta_zzzz_yyyyyz_0, ta_zzzz_yyyyyz_1, ta_zzzz_yyyyzz_0, ta_zzzz_yyyyzz_1, ta_zzzz_yyyzzz_0, ta_zzzz_yyyzzz_1, ta_zzzz_yyzzzz_0, ta_zzzz_yyzzzz_1, ta_zzzz_yzzzzz_0, ta_zzzz_yzzzzz_1, ta_zzzz_zzzzzz_0, ta_zzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzzz_xxxxxx_0[i] = ta_zzzz_xxxxxx_0[i] * fe_0 - ta_zzzz_xxxxxx_1[i] * fe_0 + ta_yzzzz_xxxxxx_0[i] * pa_y[i] - ta_yzzzz_xxxxxx_1[i] * pc_y[i];

        ta_yyzzzz_xxxxxy_0[i] = 3.0 * ta_yyzz_xxxxxy_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxxy_1[i] * fe_0 + ta_yyzzz_xxxxxy_0[i] * pa_z[i] - ta_yyzzz_xxxxxy_1[i] * pc_z[i];

        ta_yyzzzz_xxxxxz_0[i] = ta_zzzz_xxxxxz_0[i] * fe_0 - ta_zzzz_xxxxxz_1[i] * fe_0 + ta_yzzzz_xxxxxz_0[i] * pa_y[i] - ta_yzzzz_xxxxxz_1[i] * pc_y[i];

        ta_yyzzzz_xxxxyy_0[i] = 3.0 * ta_yyzz_xxxxyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxyy_1[i] * fe_0 + ta_yyzzz_xxxxyy_0[i] * pa_z[i] - ta_yyzzz_xxxxyy_1[i] * pc_z[i];

        ta_yyzzzz_xxxxyz_0[i] = ta_zzzz_xxxxyz_0[i] * fe_0 - ta_zzzz_xxxxyz_1[i] * fe_0 + ta_yzzzz_xxxxz_0[i] * fe_0 - ta_yzzzz_xxxxz_1[i] * fe_0 + ta_yzzzz_xxxxyz_0[i] * pa_y[i] - ta_yzzzz_xxxxyz_1[i] * pc_y[i];

        ta_yyzzzz_xxxxzz_0[i] = ta_zzzz_xxxxzz_0[i] * fe_0 - ta_zzzz_xxxxzz_1[i] * fe_0 + ta_yzzzz_xxxxzz_0[i] * pa_y[i] - ta_yzzzz_xxxxzz_1[i] * pc_y[i];

        ta_yyzzzz_xxxyyy_0[i] = 3.0 * ta_yyzz_xxxyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxxyyy_1[i] * fe_0 + ta_yyzzz_xxxyyy_0[i] * pa_z[i] - ta_yyzzz_xxxyyy_1[i] * pc_z[i];

        ta_yyzzzz_xxxyyz_0[i] = ta_zzzz_xxxyyz_0[i] * fe_0 - ta_zzzz_xxxyyz_1[i] * fe_0 + 2.0 * ta_yzzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_yzzzz_xxxyz_1[i] * fe_0 + ta_yzzzz_xxxyyz_0[i] * pa_y[i] - ta_yzzzz_xxxyyz_1[i] * pc_y[i];

        ta_yyzzzz_xxxyzz_0[i] = ta_zzzz_xxxyzz_0[i] * fe_0 - ta_zzzz_xxxyzz_1[i] * fe_0 + ta_yzzzz_xxxzz_0[i] * fe_0 - ta_yzzzz_xxxzz_1[i] * fe_0 + ta_yzzzz_xxxyzz_0[i] * pa_y[i] - ta_yzzzz_xxxyzz_1[i] * pc_y[i];

        ta_yyzzzz_xxxzzz_0[i] = ta_zzzz_xxxzzz_0[i] * fe_0 - ta_zzzz_xxxzzz_1[i] * fe_0 + ta_yzzzz_xxxzzz_0[i] * pa_y[i] - ta_yzzzz_xxxzzz_1[i] * pc_y[i];

        ta_yyzzzz_xxyyyy_0[i] = 3.0 * ta_yyzz_xxyyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyyy_1[i] * fe_0 + ta_yyzzz_xxyyyy_0[i] * pa_z[i] - ta_yyzzz_xxyyyy_1[i] * pc_z[i];

        ta_yyzzzz_xxyyyz_0[i] = ta_zzzz_xxyyyz_0[i] * fe_0 - ta_zzzz_xxyyyz_1[i] * fe_0 + 3.0 * ta_yzzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yzzzz_xxyyz_1[i] * fe_0 + ta_yzzzz_xxyyyz_0[i] * pa_y[i] - ta_yzzzz_xxyyyz_1[i] * pc_y[i];

        ta_yyzzzz_xxyyzz_0[i] = ta_zzzz_xxyyzz_0[i] * fe_0 - ta_zzzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_yzzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_yzzzz_xxyzz_1[i] * fe_0 + ta_yzzzz_xxyyzz_0[i] * pa_y[i] - ta_yzzzz_xxyyzz_1[i] * pc_y[i];

        ta_yyzzzz_xxyzzz_0[i] = ta_zzzz_xxyzzz_0[i] * fe_0 - ta_zzzz_xxyzzz_1[i] * fe_0 + ta_yzzzz_xxzzz_0[i] * fe_0 - ta_yzzzz_xxzzz_1[i] * fe_0 + ta_yzzzz_xxyzzz_0[i] * pa_y[i] - ta_yzzzz_xxyzzz_1[i] * pc_y[i];

        ta_yyzzzz_xxzzzz_0[i] = ta_zzzz_xxzzzz_0[i] * fe_0 - ta_zzzz_xxzzzz_1[i] * fe_0 + ta_yzzzz_xxzzzz_0[i] * pa_y[i] - ta_yzzzz_xxzzzz_1[i] * pc_y[i];

        ta_yyzzzz_xyyyyy_0[i] = 3.0 * ta_yyzz_xyyyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xyyyyy_1[i] * fe_0 + ta_yyzzz_xyyyyy_0[i] * pa_z[i] - ta_yyzzz_xyyyyy_1[i] * pc_z[i];

        ta_yyzzzz_xyyyyz_0[i] = ta_zzzz_xyyyyz_0[i] * fe_0 - ta_zzzz_xyyyyz_1[i] * fe_0 + 4.0 * ta_yzzzz_xyyyz_0[i] * fe_0 - 4.0 * ta_yzzzz_xyyyz_1[i] * fe_0 + ta_yzzzz_xyyyyz_0[i] * pa_y[i] - ta_yzzzz_xyyyyz_1[i] * pc_y[i];

        ta_yyzzzz_xyyyzz_0[i] = ta_zzzz_xyyyzz_0[i] * fe_0 - ta_zzzz_xyyyzz_1[i] * fe_0 + 3.0 * ta_yzzzz_xyyzz_0[i] * fe_0 - 3.0 * ta_yzzzz_xyyzz_1[i] * fe_0 + ta_yzzzz_xyyyzz_0[i] * pa_y[i] - ta_yzzzz_xyyyzz_1[i] * pc_y[i];

        ta_yyzzzz_xyyzzz_0[i] = ta_zzzz_xyyzzz_0[i] * fe_0 - ta_zzzz_xyyzzz_1[i] * fe_0 + 2.0 * ta_yzzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyzzz_1[i] * fe_0 + ta_yzzzz_xyyzzz_0[i] * pa_y[i] - ta_yzzzz_xyyzzz_1[i] * pc_y[i];

        ta_yyzzzz_xyzzzz_0[i] = ta_zzzz_xyzzzz_0[i] * fe_0 - ta_zzzz_xyzzzz_1[i] * fe_0 + ta_yzzzz_xzzzz_0[i] * fe_0 - ta_yzzzz_xzzzz_1[i] * fe_0 + ta_yzzzz_xyzzzz_0[i] * pa_y[i] - ta_yzzzz_xyzzzz_1[i] * pc_y[i];

        ta_yyzzzz_xzzzzz_0[i] = ta_zzzz_xzzzzz_0[i] * fe_0 - ta_zzzz_xzzzzz_1[i] * fe_0 + ta_yzzzz_xzzzzz_0[i] * pa_y[i] - ta_yzzzz_xzzzzz_1[i] * pc_y[i];

        ta_yyzzzz_yyyyyy_0[i] = 3.0 * ta_yyzz_yyyyyy_0[i] * fe_0 - 3.0 * ta_yyzz_yyyyyy_1[i] * fe_0 + ta_yyzzz_yyyyyy_0[i] * pa_z[i] - ta_yyzzz_yyyyyy_1[i] * pc_z[i];

        ta_yyzzzz_yyyyyz_0[i] = ta_zzzz_yyyyyz_0[i] * fe_0 - ta_zzzz_yyyyyz_1[i] * fe_0 + 5.0 * ta_yzzzz_yyyyz_0[i] * fe_0 - 5.0 * ta_yzzzz_yyyyz_1[i] * fe_0 + ta_yzzzz_yyyyyz_0[i] * pa_y[i] - ta_yzzzz_yyyyyz_1[i] * pc_y[i];

        ta_yyzzzz_yyyyzz_0[i] = ta_zzzz_yyyyzz_0[i] * fe_0 - ta_zzzz_yyyyzz_1[i] * fe_0 + 4.0 * ta_yzzzz_yyyzz_0[i] * fe_0 - 4.0 * ta_yzzzz_yyyzz_1[i] * fe_0 + ta_yzzzz_yyyyzz_0[i] * pa_y[i] - ta_yzzzz_yyyyzz_1[i] * pc_y[i];

        ta_yyzzzz_yyyzzz_0[i] = ta_zzzz_yyyzzz_0[i] * fe_0 - ta_zzzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_yzzzz_yyzzz_0[i] * fe_0 - 3.0 * ta_yzzzz_yyzzz_1[i] * fe_0 + ta_yzzzz_yyyzzz_0[i] * pa_y[i] - ta_yzzzz_yyyzzz_1[i] * pc_y[i];

        ta_yyzzzz_yyzzzz_0[i] = ta_zzzz_yyzzzz_0[i] * fe_0 - ta_zzzz_yyzzzz_1[i] * fe_0 + 2.0 * ta_yzzzz_yzzzz_0[i] * fe_0 - 2.0 * ta_yzzzz_yzzzz_1[i] * fe_0 + ta_yzzzz_yyzzzz_0[i] * pa_y[i] - ta_yzzzz_yyzzzz_1[i] * pc_y[i];

        ta_yyzzzz_yzzzzz_0[i] = ta_zzzz_yzzzzz_0[i] * fe_0 - ta_zzzz_yzzzzz_1[i] * fe_0 + ta_yzzzz_zzzzz_0[i] * fe_0 - ta_yzzzz_zzzzz_1[i] * fe_0 + ta_yzzzz_yzzzzz_0[i] * pa_y[i] - ta_yzzzz_yzzzzz_1[i] * pc_y[i];

        ta_yyzzzz_zzzzzz_0[i] = ta_zzzz_zzzzzz_0[i] * fe_0 - ta_zzzz_zzzzzz_1[i] * fe_0 + ta_yzzzz_zzzzzz_0[i] * pa_y[i] - ta_yzzzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 728-756 components of targeted buffer : II

    auto ta_yzzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 728);

    auto ta_yzzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 729);

    auto ta_yzzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 730);

    auto ta_yzzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 731);

    auto ta_yzzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 732);

    auto ta_yzzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 733);

    auto ta_yzzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 734);

    auto ta_yzzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 735);

    auto ta_yzzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 736);

    auto ta_yzzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 737);

    auto ta_yzzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 738);

    auto ta_yzzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 739);

    auto ta_yzzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 740);

    auto ta_yzzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 741);

    auto ta_yzzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 742);

    auto ta_yzzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 743);

    auto ta_yzzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 744);

    auto ta_yzzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 745);

    auto ta_yzzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 746);

    auto ta_yzzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 747);

    auto ta_yzzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 748);

    auto ta_yzzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 749);

    auto ta_yzzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 750);

    auto ta_yzzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 751);

    auto ta_yzzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 752);

    auto ta_yzzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 753);

    auto ta_yzzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 754);

    auto ta_yzzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 755);

    #pragma omp simd aligned(pa_y, pc_y, ta_yzzzzz_xxxxxx_0, ta_yzzzzz_xxxxxy_0, ta_yzzzzz_xxxxxz_0, ta_yzzzzz_xxxxyy_0, ta_yzzzzz_xxxxyz_0, ta_yzzzzz_xxxxzz_0, ta_yzzzzz_xxxyyy_0, ta_yzzzzz_xxxyyz_0, ta_yzzzzz_xxxyzz_0, ta_yzzzzz_xxxzzz_0, ta_yzzzzz_xxyyyy_0, ta_yzzzzz_xxyyyz_0, ta_yzzzzz_xxyyzz_0, ta_yzzzzz_xxyzzz_0, ta_yzzzzz_xxzzzz_0, ta_yzzzzz_xyyyyy_0, ta_yzzzzz_xyyyyz_0, ta_yzzzzz_xyyyzz_0, ta_yzzzzz_xyyzzz_0, ta_yzzzzz_xyzzzz_0, ta_yzzzzz_xzzzzz_0, ta_yzzzzz_yyyyyy_0, ta_yzzzzz_yyyyyz_0, ta_yzzzzz_yyyyzz_0, ta_yzzzzz_yyyzzz_0, ta_yzzzzz_yyzzzz_0, ta_yzzzzz_yzzzzz_0, ta_yzzzzz_zzzzzz_0, ta_zzzzz_xxxxx_0, ta_zzzzz_xxxxx_1, ta_zzzzz_xxxxxx_0, ta_zzzzz_xxxxxx_1, ta_zzzzz_xxxxxy_0, ta_zzzzz_xxxxxy_1, ta_zzzzz_xxxxxz_0, ta_zzzzz_xxxxxz_1, ta_zzzzz_xxxxy_0, ta_zzzzz_xxxxy_1, ta_zzzzz_xxxxyy_0, ta_zzzzz_xxxxyy_1, ta_zzzzz_xxxxyz_0, ta_zzzzz_xxxxyz_1, ta_zzzzz_xxxxz_0, ta_zzzzz_xxxxz_1, ta_zzzzz_xxxxzz_0, ta_zzzzz_xxxxzz_1, ta_zzzzz_xxxyy_0, ta_zzzzz_xxxyy_1, ta_zzzzz_xxxyyy_0, ta_zzzzz_xxxyyy_1, ta_zzzzz_xxxyyz_0, ta_zzzzz_xxxyyz_1, ta_zzzzz_xxxyz_0, ta_zzzzz_xxxyz_1, ta_zzzzz_xxxyzz_0, ta_zzzzz_xxxyzz_1, ta_zzzzz_xxxzz_0, ta_zzzzz_xxxzz_1, ta_zzzzz_xxxzzz_0, ta_zzzzz_xxxzzz_1, ta_zzzzz_xxyyy_0, ta_zzzzz_xxyyy_1, ta_zzzzz_xxyyyy_0, ta_zzzzz_xxyyyy_1, ta_zzzzz_xxyyyz_0, ta_zzzzz_xxyyyz_1, ta_zzzzz_xxyyz_0, ta_zzzzz_xxyyz_1, ta_zzzzz_xxyyzz_0, ta_zzzzz_xxyyzz_1, ta_zzzzz_xxyzz_0, ta_zzzzz_xxyzz_1, ta_zzzzz_xxyzzz_0, ta_zzzzz_xxyzzz_1, ta_zzzzz_xxzzz_0, ta_zzzzz_xxzzz_1, ta_zzzzz_xxzzzz_0, ta_zzzzz_xxzzzz_1, ta_zzzzz_xyyyy_0, ta_zzzzz_xyyyy_1, ta_zzzzz_xyyyyy_0, ta_zzzzz_xyyyyy_1, ta_zzzzz_xyyyyz_0, ta_zzzzz_xyyyyz_1, ta_zzzzz_xyyyz_0, ta_zzzzz_xyyyz_1, ta_zzzzz_xyyyzz_0, ta_zzzzz_xyyyzz_1, ta_zzzzz_xyyzz_0, ta_zzzzz_xyyzz_1, ta_zzzzz_xyyzzz_0, ta_zzzzz_xyyzzz_1, ta_zzzzz_xyzzz_0, ta_zzzzz_xyzzz_1, ta_zzzzz_xyzzzz_0, ta_zzzzz_xyzzzz_1, ta_zzzzz_xzzzz_0, ta_zzzzz_xzzzz_1, ta_zzzzz_xzzzzz_0, ta_zzzzz_xzzzzz_1, ta_zzzzz_yyyyy_0, ta_zzzzz_yyyyy_1, ta_zzzzz_yyyyyy_0, ta_zzzzz_yyyyyy_1, ta_zzzzz_yyyyyz_0, ta_zzzzz_yyyyyz_1, ta_zzzzz_yyyyz_0, ta_zzzzz_yyyyz_1, ta_zzzzz_yyyyzz_0, ta_zzzzz_yyyyzz_1, ta_zzzzz_yyyzz_0, ta_zzzzz_yyyzz_1, ta_zzzzz_yyyzzz_0, ta_zzzzz_yyyzzz_1, ta_zzzzz_yyzzz_0, ta_zzzzz_yyzzz_1, ta_zzzzz_yyzzzz_0, ta_zzzzz_yyzzzz_1, ta_zzzzz_yzzzz_0, ta_zzzzz_yzzzz_1, ta_zzzzz_yzzzzz_0, ta_zzzzz_yzzzzz_1, ta_zzzzz_zzzzz_0, ta_zzzzz_zzzzz_1, ta_zzzzz_zzzzzz_0, ta_zzzzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzzz_xxxxxx_0[i] = ta_zzzzz_xxxxxx_0[i] * pa_y[i] - ta_zzzzz_xxxxxx_1[i] * pc_y[i];

        ta_yzzzzz_xxxxxy_0[i] = ta_zzzzz_xxxxx_0[i] * fe_0 - ta_zzzzz_xxxxx_1[i] * fe_0 + ta_zzzzz_xxxxxy_0[i] * pa_y[i] - ta_zzzzz_xxxxxy_1[i] * pc_y[i];

        ta_yzzzzz_xxxxxz_0[i] = ta_zzzzz_xxxxxz_0[i] * pa_y[i] - ta_zzzzz_xxxxxz_1[i] * pc_y[i];

        ta_yzzzzz_xxxxyy_0[i] = 2.0 * ta_zzzzz_xxxxy_0[i] * fe_0 - 2.0 * ta_zzzzz_xxxxy_1[i] * fe_0 + ta_zzzzz_xxxxyy_0[i] * pa_y[i] - ta_zzzzz_xxxxyy_1[i] * pc_y[i];

        ta_yzzzzz_xxxxyz_0[i] = ta_zzzzz_xxxxz_0[i] * fe_0 - ta_zzzzz_xxxxz_1[i] * fe_0 + ta_zzzzz_xxxxyz_0[i] * pa_y[i] - ta_zzzzz_xxxxyz_1[i] * pc_y[i];

        ta_yzzzzz_xxxxzz_0[i] = ta_zzzzz_xxxxzz_0[i] * pa_y[i] - ta_zzzzz_xxxxzz_1[i] * pc_y[i];

        ta_yzzzzz_xxxyyy_0[i] = 3.0 * ta_zzzzz_xxxyy_0[i] * fe_0 - 3.0 * ta_zzzzz_xxxyy_1[i] * fe_0 + ta_zzzzz_xxxyyy_0[i] * pa_y[i] - ta_zzzzz_xxxyyy_1[i] * pc_y[i];

        ta_yzzzzz_xxxyyz_0[i] = 2.0 * ta_zzzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xxxyz_1[i] * fe_0 + ta_zzzzz_xxxyyz_0[i] * pa_y[i] - ta_zzzzz_xxxyyz_1[i] * pc_y[i];

        ta_yzzzzz_xxxyzz_0[i] = ta_zzzzz_xxxzz_0[i] * fe_0 - ta_zzzzz_xxxzz_1[i] * fe_0 + ta_zzzzz_xxxyzz_0[i] * pa_y[i] - ta_zzzzz_xxxyzz_1[i] * pc_y[i];

        ta_yzzzzz_xxxzzz_0[i] = ta_zzzzz_xxxzzz_0[i] * pa_y[i] - ta_zzzzz_xxxzzz_1[i] * pc_y[i];

        ta_yzzzzz_xxyyyy_0[i] = 4.0 * ta_zzzzz_xxyyy_0[i] * fe_0 - 4.0 * ta_zzzzz_xxyyy_1[i] * fe_0 + ta_zzzzz_xxyyyy_0[i] * pa_y[i] - ta_zzzzz_xxyyyy_1[i] * pc_y[i];

        ta_yzzzzz_xxyyyz_0[i] = 3.0 * ta_zzzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyyz_1[i] * fe_0 + ta_zzzzz_xxyyyz_0[i] * pa_y[i] - ta_zzzzz_xxyyyz_1[i] * pc_y[i];

        ta_yzzzzz_xxyyzz_0[i] = 2.0 * ta_zzzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xxyzz_1[i] * fe_0 + ta_zzzzz_xxyyzz_0[i] * pa_y[i] - ta_zzzzz_xxyyzz_1[i] * pc_y[i];

        ta_yzzzzz_xxyzzz_0[i] = ta_zzzzz_xxzzz_0[i] * fe_0 - ta_zzzzz_xxzzz_1[i] * fe_0 + ta_zzzzz_xxyzzz_0[i] * pa_y[i] - ta_zzzzz_xxyzzz_1[i] * pc_y[i];

        ta_yzzzzz_xxzzzz_0[i] = ta_zzzzz_xxzzzz_0[i] * pa_y[i] - ta_zzzzz_xxzzzz_1[i] * pc_y[i];

        ta_yzzzzz_xyyyyy_0[i] = 5.0 * ta_zzzzz_xyyyy_0[i] * fe_0 - 5.0 * ta_zzzzz_xyyyy_1[i] * fe_0 + ta_zzzzz_xyyyyy_0[i] * pa_y[i] - ta_zzzzz_xyyyyy_1[i] * pc_y[i];

        ta_yzzzzz_xyyyyz_0[i] = 4.0 * ta_zzzzz_xyyyz_0[i] * fe_0 - 4.0 * ta_zzzzz_xyyyz_1[i] * fe_0 + ta_zzzzz_xyyyyz_0[i] * pa_y[i] - ta_zzzzz_xyyyyz_1[i] * pc_y[i];

        ta_yzzzzz_xyyyzz_0[i] = 3.0 * ta_zzzzz_xyyzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xyyzz_1[i] * fe_0 + ta_zzzzz_xyyyzz_0[i] * pa_y[i] - ta_zzzzz_xyyyzz_1[i] * pc_y[i];

        ta_yzzzzz_xyyzzz_0[i] = 2.0 * ta_zzzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyzzz_1[i] * fe_0 + ta_zzzzz_xyyzzz_0[i] * pa_y[i] - ta_zzzzz_xyyzzz_1[i] * pc_y[i];

        ta_yzzzzz_xyzzzz_0[i] = ta_zzzzz_xzzzz_0[i] * fe_0 - ta_zzzzz_xzzzz_1[i] * fe_0 + ta_zzzzz_xyzzzz_0[i] * pa_y[i] - ta_zzzzz_xyzzzz_1[i] * pc_y[i];

        ta_yzzzzz_xzzzzz_0[i] = ta_zzzzz_xzzzzz_0[i] * pa_y[i] - ta_zzzzz_xzzzzz_1[i] * pc_y[i];

        ta_yzzzzz_yyyyyy_0[i] = 6.0 * ta_zzzzz_yyyyy_0[i] * fe_0 - 6.0 * ta_zzzzz_yyyyy_1[i] * fe_0 + ta_zzzzz_yyyyyy_0[i] * pa_y[i] - ta_zzzzz_yyyyyy_1[i] * pc_y[i];

        ta_yzzzzz_yyyyyz_0[i] = 5.0 * ta_zzzzz_yyyyz_0[i] * fe_0 - 5.0 * ta_zzzzz_yyyyz_1[i] * fe_0 + ta_zzzzz_yyyyyz_0[i] * pa_y[i] - ta_zzzzz_yyyyyz_1[i] * pc_y[i];

        ta_yzzzzz_yyyyzz_0[i] = 4.0 * ta_zzzzz_yyyzz_0[i] * fe_0 - 4.0 * ta_zzzzz_yyyzz_1[i] * fe_0 + ta_zzzzz_yyyyzz_0[i] * pa_y[i] - ta_zzzzz_yyyyzz_1[i] * pc_y[i];

        ta_yzzzzz_yyyzzz_0[i] = 3.0 * ta_zzzzz_yyzzz_0[i] * fe_0 - 3.0 * ta_zzzzz_yyzzz_1[i] * fe_0 + ta_zzzzz_yyyzzz_0[i] * pa_y[i] - ta_zzzzz_yyyzzz_1[i] * pc_y[i];

        ta_yzzzzz_yyzzzz_0[i] = 2.0 * ta_zzzzz_yzzzz_0[i] * fe_0 - 2.0 * ta_zzzzz_yzzzz_1[i] * fe_0 + ta_zzzzz_yyzzzz_0[i] * pa_y[i] - ta_zzzzz_yyzzzz_1[i] * pc_y[i];

        ta_yzzzzz_yzzzzz_0[i] = ta_zzzzz_zzzzz_0[i] * fe_0 - ta_zzzzz_zzzzz_1[i] * fe_0 + ta_zzzzz_yzzzzz_0[i] * pa_y[i] - ta_zzzzz_yzzzzz_1[i] * pc_y[i];

        ta_yzzzzz_zzzzzz_0[i] = ta_zzzzz_zzzzzz_0[i] * pa_y[i] - ta_zzzzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 756-784 components of targeted buffer : II

    auto ta_zzzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_ii + 756);

    auto ta_zzzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_ii + 757);

    auto ta_zzzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_ii + 758);

    auto ta_zzzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_ii + 759);

    auto ta_zzzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_ii + 760);

    auto ta_zzzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_ii + 761);

    auto ta_zzzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_ii + 762);

    auto ta_zzzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_ii + 763);

    auto ta_zzzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_ii + 764);

    auto ta_zzzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_ii + 765);

    auto ta_zzzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_ii + 766);

    auto ta_zzzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_ii + 767);

    auto ta_zzzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_ii + 768);

    auto ta_zzzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_ii + 769);

    auto ta_zzzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_ii + 770);

    auto ta_zzzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_ii + 771);

    auto ta_zzzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_ii + 772);

    auto ta_zzzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_ii + 773);

    auto ta_zzzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_ii + 774);

    auto ta_zzzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_ii + 775);

    auto ta_zzzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_ii + 776);

    auto ta_zzzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_ii + 777);

    auto ta_zzzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_ii + 778);

    auto ta_zzzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_ii + 779);

    auto ta_zzzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_ii + 780);

    auto ta_zzzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_ii + 781);

    auto ta_zzzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_ii + 782);

    auto ta_zzzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_ii + 783);

    #pragma omp simd aligned(pa_z, pc_z, ta_zzzz_xxxxxx_0, ta_zzzz_xxxxxx_1, ta_zzzz_xxxxxy_0, ta_zzzz_xxxxxy_1, ta_zzzz_xxxxxz_0, ta_zzzz_xxxxxz_1, ta_zzzz_xxxxyy_0, ta_zzzz_xxxxyy_1, ta_zzzz_xxxxyz_0, ta_zzzz_xxxxyz_1, ta_zzzz_xxxxzz_0, ta_zzzz_xxxxzz_1, ta_zzzz_xxxyyy_0, ta_zzzz_xxxyyy_1, ta_zzzz_xxxyyz_0, ta_zzzz_xxxyyz_1, ta_zzzz_xxxyzz_0, ta_zzzz_xxxyzz_1, ta_zzzz_xxxzzz_0, ta_zzzz_xxxzzz_1, ta_zzzz_xxyyyy_0, ta_zzzz_xxyyyy_1, ta_zzzz_xxyyyz_0, ta_zzzz_xxyyyz_1, ta_zzzz_xxyyzz_0, ta_zzzz_xxyyzz_1, ta_zzzz_xxyzzz_0, ta_zzzz_xxyzzz_1, ta_zzzz_xxzzzz_0, ta_zzzz_xxzzzz_1, ta_zzzz_xyyyyy_0, ta_zzzz_xyyyyy_1, ta_zzzz_xyyyyz_0, ta_zzzz_xyyyyz_1, ta_zzzz_xyyyzz_0, ta_zzzz_xyyyzz_1, ta_zzzz_xyyzzz_0, ta_zzzz_xyyzzz_1, ta_zzzz_xyzzzz_0, ta_zzzz_xyzzzz_1, ta_zzzz_xzzzzz_0, ta_zzzz_xzzzzz_1, ta_zzzz_yyyyyy_0, ta_zzzz_yyyyyy_1, ta_zzzz_yyyyyz_0, ta_zzzz_yyyyyz_1, ta_zzzz_yyyyzz_0, ta_zzzz_yyyyzz_1, ta_zzzz_yyyzzz_0, ta_zzzz_yyyzzz_1, ta_zzzz_yyzzzz_0, ta_zzzz_yyzzzz_1, ta_zzzz_yzzzzz_0, ta_zzzz_yzzzzz_1, ta_zzzz_zzzzzz_0, ta_zzzz_zzzzzz_1, ta_zzzzz_xxxxx_0, ta_zzzzz_xxxxx_1, ta_zzzzz_xxxxxx_0, ta_zzzzz_xxxxxx_1, ta_zzzzz_xxxxxy_0, ta_zzzzz_xxxxxy_1, ta_zzzzz_xxxxxz_0, ta_zzzzz_xxxxxz_1, ta_zzzzz_xxxxy_0, ta_zzzzz_xxxxy_1, ta_zzzzz_xxxxyy_0, ta_zzzzz_xxxxyy_1, ta_zzzzz_xxxxyz_0, ta_zzzzz_xxxxyz_1, ta_zzzzz_xxxxz_0, ta_zzzzz_xxxxz_1, ta_zzzzz_xxxxzz_0, ta_zzzzz_xxxxzz_1, ta_zzzzz_xxxyy_0, ta_zzzzz_xxxyy_1, ta_zzzzz_xxxyyy_0, ta_zzzzz_xxxyyy_1, ta_zzzzz_xxxyyz_0, ta_zzzzz_xxxyyz_1, ta_zzzzz_xxxyz_0, ta_zzzzz_xxxyz_1, ta_zzzzz_xxxyzz_0, ta_zzzzz_xxxyzz_1, ta_zzzzz_xxxzz_0, ta_zzzzz_xxxzz_1, ta_zzzzz_xxxzzz_0, ta_zzzzz_xxxzzz_1, ta_zzzzz_xxyyy_0, ta_zzzzz_xxyyy_1, ta_zzzzz_xxyyyy_0, ta_zzzzz_xxyyyy_1, ta_zzzzz_xxyyyz_0, ta_zzzzz_xxyyyz_1, ta_zzzzz_xxyyz_0, ta_zzzzz_xxyyz_1, ta_zzzzz_xxyyzz_0, ta_zzzzz_xxyyzz_1, ta_zzzzz_xxyzz_0, ta_zzzzz_xxyzz_1, ta_zzzzz_xxyzzz_0, ta_zzzzz_xxyzzz_1, ta_zzzzz_xxzzz_0, ta_zzzzz_xxzzz_1, ta_zzzzz_xxzzzz_0, ta_zzzzz_xxzzzz_1, ta_zzzzz_xyyyy_0, ta_zzzzz_xyyyy_1, ta_zzzzz_xyyyyy_0, ta_zzzzz_xyyyyy_1, ta_zzzzz_xyyyyz_0, ta_zzzzz_xyyyyz_1, ta_zzzzz_xyyyz_0, ta_zzzzz_xyyyz_1, ta_zzzzz_xyyyzz_0, ta_zzzzz_xyyyzz_1, ta_zzzzz_xyyzz_0, ta_zzzzz_xyyzz_1, ta_zzzzz_xyyzzz_0, ta_zzzzz_xyyzzz_1, ta_zzzzz_xyzzz_0, ta_zzzzz_xyzzz_1, ta_zzzzz_xyzzzz_0, ta_zzzzz_xyzzzz_1, ta_zzzzz_xzzzz_0, ta_zzzzz_xzzzz_1, ta_zzzzz_xzzzzz_0, ta_zzzzz_xzzzzz_1, ta_zzzzz_yyyyy_0, ta_zzzzz_yyyyy_1, ta_zzzzz_yyyyyy_0, ta_zzzzz_yyyyyy_1, ta_zzzzz_yyyyyz_0, ta_zzzzz_yyyyyz_1, ta_zzzzz_yyyyz_0, ta_zzzzz_yyyyz_1, ta_zzzzz_yyyyzz_0, ta_zzzzz_yyyyzz_1, ta_zzzzz_yyyzz_0, ta_zzzzz_yyyzz_1, ta_zzzzz_yyyzzz_0, ta_zzzzz_yyyzzz_1, ta_zzzzz_yyzzz_0, ta_zzzzz_yyzzz_1, ta_zzzzz_yyzzzz_0, ta_zzzzz_yyzzzz_1, ta_zzzzz_yzzzz_0, ta_zzzzz_yzzzz_1, ta_zzzzz_yzzzzz_0, ta_zzzzz_yzzzzz_1, ta_zzzzz_zzzzz_0, ta_zzzzz_zzzzz_1, ta_zzzzz_zzzzzz_0, ta_zzzzz_zzzzzz_1, ta_zzzzzz_xxxxxx_0, ta_zzzzzz_xxxxxy_0, ta_zzzzzz_xxxxxz_0, ta_zzzzzz_xxxxyy_0, ta_zzzzzz_xxxxyz_0, ta_zzzzzz_xxxxzz_0, ta_zzzzzz_xxxyyy_0, ta_zzzzzz_xxxyyz_0, ta_zzzzzz_xxxyzz_0, ta_zzzzzz_xxxzzz_0, ta_zzzzzz_xxyyyy_0, ta_zzzzzz_xxyyyz_0, ta_zzzzzz_xxyyzz_0, ta_zzzzzz_xxyzzz_0, ta_zzzzzz_xxzzzz_0, ta_zzzzzz_xyyyyy_0, ta_zzzzzz_xyyyyz_0, ta_zzzzzz_xyyyzz_0, ta_zzzzzz_xyyzzz_0, ta_zzzzzz_xyzzzz_0, ta_zzzzzz_xzzzzz_0, ta_zzzzzz_yyyyyy_0, ta_zzzzzz_yyyyyz_0, ta_zzzzzz_yyyyzz_0, ta_zzzzzz_yyyzzz_0, ta_zzzzzz_yyzzzz_0, ta_zzzzzz_yzzzzz_0, ta_zzzzzz_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzzz_xxxxxx_0[i] = 5.0 * ta_zzzz_xxxxxx_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxxx_1[i] * fe_0 + ta_zzzzz_xxxxxx_0[i] * pa_z[i] - ta_zzzzz_xxxxxx_1[i] * pc_z[i];

        ta_zzzzzz_xxxxxy_0[i] = 5.0 * ta_zzzz_xxxxxy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxxy_1[i] * fe_0 + ta_zzzzz_xxxxxy_0[i] * pa_z[i] - ta_zzzzz_xxxxxy_1[i] * pc_z[i];

        ta_zzzzzz_xxxxxz_0[i] = 5.0 * ta_zzzz_xxxxxz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxxz_1[i] * fe_0 + ta_zzzzz_xxxxx_0[i] * fe_0 - ta_zzzzz_xxxxx_1[i] * fe_0 + ta_zzzzz_xxxxxz_0[i] * pa_z[i] - ta_zzzzz_xxxxxz_1[i] * pc_z[i];

        ta_zzzzzz_xxxxyy_0[i] = 5.0 * ta_zzzz_xxxxyy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxyy_1[i] * fe_0 + ta_zzzzz_xxxxyy_0[i] * pa_z[i] - ta_zzzzz_xxxxyy_1[i] * pc_z[i];

        ta_zzzzzz_xxxxyz_0[i] = 5.0 * ta_zzzz_xxxxyz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxyz_1[i] * fe_0 + ta_zzzzz_xxxxy_0[i] * fe_0 - ta_zzzzz_xxxxy_1[i] * fe_0 + ta_zzzzz_xxxxyz_0[i] * pa_z[i] - ta_zzzzz_xxxxyz_1[i] * pc_z[i];

        ta_zzzzzz_xxxxzz_0[i] = 5.0 * ta_zzzz_xxxxzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xxxxz_0[i] * fe_0 - 2.0 * ta_zzzzz_xxxxz_1[i] * fe_0 + ta_zzzzz_xxxxzz_0[i] * pa_z[i] - ta_zzzzz_xxxxzz_1[i] * pc_z[i];

        ta_zzzzzz_xxxyyy_0[i] = 5.0 * ta_zzzz_xxxyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxyyy_1[i] * fe_0 + ta_zzzzz_xxxyyy_0[i] * pa_z[i] - ta_zzzzz_xxxyyy_1[i] * pc_z[i];

        ta_zzzzzz_xxxyyz_0[i] = 5.0 * ta_zzzz_xxxyyz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxyyz_1[i] * fe_0 + ta_zzzzz_xxxyy_0[i] * fe_0 - ta_zzzzz_xxxyy_1[i] * fe_0 + ta_zzzzz_xxxyyz_0[i] * pa_z[i] - ta_zzzzz_xxxyyz_1[i] * pc_z[i];

        ta_zzzzzz_xxxyzz_0[i] = 5.0 * ta_zzzz_xxxyzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xxxyz_1[i] * fe_0 + ta_zzzzz_xxxyzz_0[i] * pa_z[i] - ta_zzzzz_xxxyzz_1[i] * pc_z[i];

        ta_zzzzzz_xxxzzz_0[i] = 5.0 * ta_zzzz_xxxzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_xxxzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxxzz_1[i] * fe_0 + ta_zzzzz_xxxzzz_0[i] * pa_z[i] - ta_zzzzz_xxxzzz_1[i] * pc_z[i];

        ta_zzzzzz_xxyyyy_0[i] = 5.0 * ta_zzzz_xxyyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xxyyyy_1[i] * fe_0 + ta_zzzzz_xxyyyy_0[i] * pa_z[i] - ta_zzzzz_xxyyyy_1[i] * pc_z[i];

        ta_zzzzzz_xxyyyz_0[i] = 5.0 * ta_zzzz_xxyyyz_0[i] * fe_0 - 5.0 * ta_zzzz_xxyyyz_1[i] * fe_0 + ta_zzzzz_xxyyy_0[i] * fe_0 - ta_zzzzz_xxyyy_1[i] * fe_0 + ta_zzzzz_xxyyyz_0[i] * pa_z[i] - ta_zzzzz_xxyyyz_1[i] * pc_z[i];

        ta_zzzzzz_xxyyzz_0[i] = 5.0 * ta_zzzz_xxyyzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xxyyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xxyyz_1[i] * fe_0 + ta_zzzzz_xxyyzz_0[i] * pa_z[i] - ta_zzzzz_xxyyzz_1[i] * pc_z[i];

        ta_zzzzzz_xxyzzz_0[i] = 5.0 * ta_zzzz_xxyzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxyzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyzz_1[i] * fe_0 + ta_zzzzz_xxyzzz_0[i] * pa_z[i] - ta_zzzzz_xxyzzz_1[i] * pc_z[i];

        ta_zzzzzz_xxzzzz_0[i] = 5.0 * ta_zzzz_xxzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxzzzz_1[i] * fe_0 + 4.0 * ta_zzzzz_xxzzz_0[i] * fe_0 - 4.0 * ta_zzzzz_xxzzz_1[i] * fe_0 + ta_zzzzz_xxzzzz_0[i] * pa_z[i] - ta_zzzzz_xxzzzz_1[i] * pc_z[i];

        ta_zzzzzz_xyyyyy_0[i] = 5.0 * ta_zzzz_xyyyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xyyyyy_1[i] * fe_0 + ta_zzzzz_xyyyyy_0[i] * pa_z[i] - ta_zzzzz_xyyyyy_1[i] * pc_z[i];

        ta_zzzzzz_xyyyyz_0[i] = 5.0 * ta_zzzz_xyyyyz_0[i] * fe_0 - 5.0 * ta_zzzz_xyyyyz_1[i] * fe_0 + ta_zzzzz_xyyyy_0[i] * fe_0 - ta_zzzzz_xyyyy_1[i] * fe_0 + ta_zzzzz_xyyyyz_0[i] * pa_z[i] - ta_zzzzz_xyyyyz_1[i] * pc_z[i];

        ta_zzzzzz_xyyyzz_0[i] = 5.0 * ta_zzzz_xyyyzz_0[i] * fe_0 - 5.0 * ta_zzzz_xyyyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyyyz_1[i] * fe_0 + ta_zzzzz_xyyyzz_0[i] * pa_z[i] - ta_zzzzz_xyyyzz_1[i] * pc_z[i];

        ta_zzzzzz_xyyzzz_0[i] = 5.0 * ta_zzzz_xyyzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xyyzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_xyyzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xyyzz_1[i] * fe_0 + ta_zzzzz_xyyzzz_0[i] * pa_z[i] - ta_zzzzz_xyyzzz_1[i] * pc_z[i];

        ta_zzzzzz_xyzzzz_0[i] = 5.0 * ta_zzzz_xyzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xyzzzz_1[i] * fe_0 + 4.0 * ta_zzzzz_xyzzz_0[i] * fe_0 - 4.0 * ta_zzzzz_xyzzz_1[i] * fe_0 + ta_zzzzz_xyzzzz_0[i] * pa_z[i] - ta_zzzzz_xyzzzz_1[i] * pc_z[i];

        ta_zzzzzz_xzzzzz_0[i] = 5.0 * ta_zzzz_xzzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xzzzzz_1[i] * fe_0 + 5.0 * ta_zzzzz_xzzzz_0[i] * fe_0 - 5.0 * ta_zzzzz_xzzzz_1[i] * fe_0 + ta_zzzzz_xzzzzz_0[i] * pa_z[i] - ta_zzzzz_xzzzzz_1[i] * pc_z[i];

        ta_zzzzzz_yyyyyy_0[i] = 5.0 * ta_zzzz_yyyyyy_0[i] * fe_0 - 5.0 * ta_zzzz_yyyyyy_1[i] * fe_0 + ta_zzzzz_yyyyyy_0[i] * pa_z[i] - ta_zzzzz_yyyyyy_1[i] * pc_z[i];

        ta_zzzzzz_yyyyyz_0[i] = 5.0 * ta_zzzz_yyyyyz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyyyz_1[i] * fe_0 + ta_zzzzz_yyyyy_0[i] * fe_0 - ta_zzzzz_yyyyy_1[i] * fe_0 + ta_zzzzz_yyyyyz_0[i] * pa_z[i] - ta_zzzzz_yyyyyz_1[i] * pc_z[i];

        ta_zzzzzz_yyyyzz_0[i] = 5.0 * ta_zzzz_yyyyzz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_yyyyz_0[i] * fe_0 - 2.0 * ta_zzzzz_yyyyz_1[i] * fe_0 + ta_zzzzz_yyyyzz_0[i] * pa_z[i] - ta_zzzzz_yyyyzz_1[i] * pc_z[i];

        ta_zzzzzz_yyyzzz_0[i] = 5.0 * ta_zzzz_yyyzzz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_yyyzz_0[i] * fe_0 - 3.0 * ta_zzzzz_yyyzz_1[i] * fe_0 + ta_zzzzz_yyyzzz_0[i] * pa_z[i] - ta_zzzzz_yyyzzz_1[i] * pc_z[i];

        ta_zzzzzz_yyzzzz_0[i] = 5.0 * ta_zzzz_yyzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_yyzzzz_1[i] * fe_0 + 4.0 * ta_zzzzz_yyzzz_0[i] * fe_0 - 4.0 * ta_zzzzz_yyzzz_1[i] * fe_0 + ta_zzzzz_yyzzzz_0[i] * pa_z[i] - ta_zzzzz_yyzzzz_1[i] * pc_z[i];

        ta_zzzzzz_yzzzzz_0[i] = 5.0 * ta_zzzz_yzzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_yzzzzz_1[i] * fe_0 + 5.0 * ta_zzzzz_yzzzz_0[i] * fe_0 - 5.0 * ta_zzzzz_yzzzz_1[i] * fe_0 + ta_zzzzz_yzzzzz_0[i] * pa_z[i] - ta_zzzzz_yzzzzz_1[i] * pc_z[i];

        ta_zzzzzz_zzzzzz_0[i] = 5.0 * ta_zzzz_zzzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_zzzzzz_1[i] * fe_0 + 6.0 * ta_zzzzz_zzzzz_0[i] * fe_0 - 6.0 * ta_zzzzz_zzzzz_1[i] * fe_0 + ta_zzzzz_zzzzzz_0[i] * pa_z[i] - ta_zzzzz_zzzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

