#include "ElectricDipoleMomentumPrimRecGI.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_gi(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_gi,
                                      const size_t idx_dip_di,
                                      const size_t idx_dip_fh,
                                      const size_t idx_ovl_fi,
                                      const size_t idx_dip_fi,
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

    // Set up components of auxiliary buffer : DI

    auto tr_x_xx_xxxxxx = pbuffer.data(idx_dip_di);

    auto tr_x_xx_xxxxxy = pbuffer.data(idx_dip_di + 1);

    auto tr_x_xx_xxxxxz = pbuffer.data(idx_dip_di + 2);

    auto tr_x_xx_xxxxyy = pbuffer.data(idx_dip_di + 3);

    auto tr_x_xx_xxxxyz = pbuffer.data(idx_dip_di + 4);

    auto tr_x_xx_xxxxzz = pbuffer.data(idx_dip_di + 5);

    auto tr_x_xx_xxxyyy = pbuffer.data(idx_dip_di + 6);

    auto tr_x_xx_xxxyyz = pbuffer.data(idx_dip_di + 7);

    auto tr_x_xx_xxxyzz = pbuffer.data(idx_dip_di + 8);

    auto tr_x_xx_xxxzzz = pbuffer.data(idx_dip_di + 9);

    auto tr_x_xx_xxyyyy = pbuffer.data(idx_dip_di + 10);

    auto tr_x_xx_xxyyyz = pbuffer.data(idx_dip_di + 11);

    auto tr_x_xx_xxyyzz = pbuffer.data(idx_dip_di + 12);

    auto tr_x_xx_xxyzzz = pbuffer.data(idx_dip_di + 13);

    auto tr_x_xx_xxzzzz = pbuffer.data(idx_dip_di + 14);

    auto tr_x_xx_xyyyyy = pbuffer.data(idx_dip_di + 15);

    auto tr_x_xx_xyyyyz = pbuffer.data(idx_dip_di + 16);

    auto tr_x_xx_xyyyzz = pbuffer.data(idx_dip_di + 17);

    auto tr_x_xx_xyyzzz = pbuffer.data(idx_dip_di + 18);

    auto tr_x_xx_xyzzzz = pbuffer.data(idx_dip_di + 19);

    auto tr_x_xx_xzzzzz = pbuffer.data(idx_dip_di + 20);

    auto tr_x_xx_yyyyyy = pbuffer.data(idx_dip_di + 21);

    auto tr_x_xx_yyyyyz = pbuffer.data(idx_dip_di + 22);

    auto tr_x_xx_yyyyzz = pbuffer.data(idx_dip_di + 23);

    auto tr_x_xx_yyyzzz = pbuffer.data(idx_dip_di + 24);

    auto tr_x_xx_yyzzzz = pbuffer.data(idx_dip_di + 25);

    auto tr_x_xx_yzzzzz = pbuffer.data(idx_dip_di + 26);

    auto tr_x_xx_zzzzzz = pbuffer.data(idx_dip_di + 27);

    auto tr_x_xy_xxxxxx = pbuffer.data(idx_dip_di + 28);

    auto tr_x_xy_xxxxxz = pbuffer.data(idx_dip_di + 30);

    auto tr_x_xy_xxxxzz = pbuffer.data(idx_dip_di + 33);

    auto tr_x_xy_xxxzzz = pbuffer.data(idx_dip_di + 37);

    auto tr_x_xy_xxzzzz = pbuffer.data(idx_dip_di + 42);

    auto tr_x_xy_xzzzzz = pbuffer.data(idx_dip_di + 48);

    auto tr_x_xz_xxxxxx = pbuffer.data(idx_dip_di + 56);

    auto tr_x_xz_xxxxxy = pbuffer.data(idx_dip_di + 57);

    auto tr_x_xz_xxxxxz = pbuffer.data(idx_dip_di + 58);

    auto tr_x_xz_xxxxyy = pbuffer.data(idx_dip_di + 59);

    auto tr_x_xz_xxxxzz = pbuffer.data(idx_dip_di + 61);

    auto tr_x_xz_xxxyyy = pbuffer.data(idx_dip_di + 62);

    auto tr_x_xz_xxxzzz = pbuffer.data(idx_dip_di + 65);

    auto tr_x_xz_xxyyyy = pbuffer.data(idx_dip_di + 66);

    auto tr_x_xz_xxzzzz = pbuffer.data(idx_dip_di + 70);

    auto tr_x_xz_xyyyyy = pbuffer.data(idx_dip_di + 71);

    auto tr_x_xz_xzzzzz = pbuffer.data(idx_dip_di + 76);

    auto tr_x_yy_xxxxxx = pbuffer.data(idx_dip_di + 84);

    auto tr_x_yy_xxxxxy = pbuffer.data(idx_dip_di + 85);

    auto tr_x_yy_xxxxxz = pbuffer.data(idx_dip_di + 86);

    auto tr_x_yy_xxxxyy = pbuffer.data(idx_dip_di + 87);

    auto tr_x_yy_xxxxyz = pbuffer.data(idx_dip_di + 88);

    auto tr_x_yy_xxxxzz = pbuffer.data(idx_dip_di + 89);

    auto tr_x_yy_xxxyyy = pbuffer.data(idx_dip_di + 90);

    auto tr_x_yy_xxxyyz = pbuffer.data(idx_dip_di + 91);

    auto tr_x_yy_xxxyzz = pbuffer.data(idx_dip_di + 92);

    auto tr_x_yy_xxxzzz = pbuffer.data(idx_dip_di + 93);

    auto tr_x_yy_xxyyyy = pbuffer.data(idx_dip_di + 94);

    auto tr_x_yy_xxyyyz = pbuffer.data(idx_dip_di + 95);

    auto tr_x_yy_xxyyzz = pbuffer.data(idx_dip_di + 96);

    auto tr_x_yy_xxyzzz = pbuffer.data(idx_dip_di + 97);

    auto tr_x_yy_xxzzzz = pbuffer.data(idx_dip_di + 98);

    auto tr_x_yy_xyyyyy = pbuffer.data(idx_dip_di + 99);

    auto tr_x_yy_xyyyyz = pbuffer.data(idx_dip_di + 100);

    auto tr_x_yy_xyyyzz = pbuffer.data(idx_dip_di + 101);

    auto tr_x_yy_xyyzzz = pbuffer.data(idx_dip_di + 102);

    auto tr_x_yy_xyzzzz = pbuffer.data(idx_dip_di + 103);

    auto tr_x_yy_xzzzzz = pbuffer.data(idx_dip_di + 104);

    auto tr_x_yy_yyyyyy = pbuffer.data(idx_dip_di + 105);

    auto tr_x_yy_yyyyyz = pbuffer.data(idx_dip_di + 106);

    auto tr_x_yy_yyyyzz = pbuffer.data(idx_dip_di + 107);

    auto tr_x_yy_yyyzzz = pbuffer.data(idx_dip_di + 108);

    auto tr_x_yy_yyzzzz = pbuffer.data(idx_dip_di + 109);

    auto tr_x_yy_yzzzzz = pbuffer.data(idx_dip_di + 110);

    auto tr_x_yy_zzzzzz = pbuffer.data(idx_dip_di + 111);

    auto tr_x_yz_xxxxxz = pbuffer.data(idx_dip_di + 114);

    auto tr_x_yz_xxxxzz = pbuffer.data(idx_dip_di + 117);

    auto tr_x_yz_xxxzzz = pbuffer.data(idx_dip_di + 121);

    auto tr_x_yz_xxzzzz = pbuffer.data(idx_dip_di + 126);

    auto tr_x_yz_xzzzzz = pbuffer.data(idx_dip_di + 132);

    auto tr_x_yz_zzzzzz = pbuffer.data(idx_dip_di + 139);

    auto tr_x_zz_xxxxxx = pbuffer.data(idx_dip_di + 140);

    auto tr_x_zz_xxxxxy = pbuffer.data(idx_dip_di + 141);

    auto tr_x_zz_xxxxxz = pbuffer.data(idx_dip_di + 142);

    auto tr_x_zz_xxxxyy = pbuffer.data(idx_dip_di + 143);

    auto tr_x_zz_xxxxyz = pbuffer.data(idx_dip_di + 144);

    auto tr_x_zz_xxxxzz = pbuffer.data(idx_dip_di + 145);

    auto tr_x_zz_xxxyyy = pbuffer.data(idx_dip_di + 146);

    auto tr_x_zz_xxxyyz = pbuffer.data(idx_dip_di + 147);

    auto tr_x_zz_xxxyzz = pbuffer.data(idx_dip_di + 148);

    auto tr_x_zz_xxxzzz = pbuffer.data(idx_dip_di + 149);

    auto tr_x_zz_xxyyyy = pbuffer.data(idx_dip_di + 150);

    auto tr_x_zz_xxyyyz = pbuffer.data(idx_dip_di + 151);

    auto tr_x_zz_xxyyzz = pbuffer.data(idx_dip_di + 152);

    auto tr_x_zz_xxyzzz = pbuffer.data(idx_dip_di + 153);

    auto tr_x_zz_xxzzzz = pbuffer.data(idx_dip_di + 154);

    auto tr_x_zz_xyyyyy = pbuffer.data(idx_dip_di + 155);

    auto tr_x_zz_xyyyyz = pbuffer.data(idx_dip_di + 156);

    auto tr_x_zz_xyyyzz = pbuffer.data(idx_dip_di + 157);

    auto tr_x_zz_xyyzzz = pbuffer.data(idx_dip_di + 158);

    auto tr_x_zz_xyzzzz = pbuffer.data(idx_dip_di + 159);

    auto tr_x_zz_xzzzzz = pbuffer.data(idx_dip_di + 160);

    auto tr_x_zz_yyyyyy = pbuffer.data(idx_dip_di + 161);

    auto tr_x_zz_yyyyyz = pbuffer.data(idx_dip_di + 162);

    auto tr_x_zz_yyyyzz = pbuffer.data(idx_dip_di + 163);

    auto tr_x_zz_yyyzzz = pbuffer.data(idx_dip_di + 164);

    auto tr_x_zz_yyzzzz = pbuffer.data(idx_dip_di + 165);

    auto tr_x_zz_yzzzzz = pbuffer.data(idx_dip_di + 166);

    auto tr_x_zz_zzzzzz = pbuffer.data(idx_dip_di + 167);

    auto tr_y_xx_xxxxxx = pbuffer.data(idx_dip_di + 168);

    auto tr_y_xx_xxxxxy = pbuffer.data(idx_dip_di + 169);

    auto tr_y_xx_xxxxxz = pbuffer.data(idx_dip_di + 170);

    auto tr_y_xx_xxxxyy = pbuffer.data(idx_dip_di + 171);

    auto tr_y_xx_xxxxyz = pbuffer.data(idx_dip_di + 172);

    auto tr_y_xx_xxxxzz = pbuffer.data(idx_dip_di + 173);

    auto tr_y_xx_xxxyyy = pbuffer.data(idx_dip_di + 174);

    auto tr_y_xx_xxxyyz = pbuffer.data(idx_dip_di + 175);

    auto tr_y_xx_xxxyzz = pbuffer.data(idx_dip_di + 176);

    auto tr_y_xx_xxxzzz = pbuffer.data(idx_dip_di + 177);

    auto tr_y_xx_xxyyyy = pbuffer.data(idx_dip_di + 178);

    auto tr_y_xx_xxyyyz = pbuffer.data(idx_dip_di + 179);

    auto tr_y_xx_xxyyzz = pbuffer.data(idx_dip_di + 180);

    auto tr_y_xx_xxyzzz = pbuffer.data(idx_dip_di + 181);

    auto tr_y_xx_xxzzzz = pbuffer.data(idx_dip_di + 182);

    auto tr_y_xx_xyyyyy = pbuffer.data(idx_dip_di + 183);

    auto tr_y_xx_xyyyyz = pbuffer.data(idx_dip_di + 184);

    auto tr_y_xx_xyyyzz = pbuffer.data(idx_dip_di + 185);

    auto tr_y_xx_xyyzzz = pbuffer.data(idx_dip_di + 186);

    auto tr_y_xx_xyzzzz = pbuffer.data(idx_dip_di + 187);

    auto tr_y_xx_xzzzzz = pbuffer.data(idx_dip_di + 188);

    auto tr_y_xx_yyyyyy = pbuffer.data(idx_dip_di + 189);

    auto tr_y_xx_yyyyyz = pbuffer.data(idx_dip_di + 190);

    auto tr_y_xx_yyyyzz = pbuffer.data(idx_dip_di + 191);

    auto tr_y_xx_yyyzzz = pbuffer.data(idx_dip_di + 192);

    auto tr_y_xx_yyzzzz = pbuffer.data(idx_dip_di + 193);

    auto tr_y_xx_yzzzzz = pbuffer.data(idx_dip_di + 194);

    auto tr_y_xx_zzzzzz = pbuffer.data(idx_dip_di + 195);

    auto tr_y_xy_xxxxxy = pbuffer.data(idx_dip_di + 197);

    auto tr_y_xy_xxxxyy = pbuffer.data(idx_dip_di + 199);

    auto tr_y_xy_xxxxyz = pbuffer.data(idx_dip_di + 200);

    auto tr_y_xy_xxxyyy = pbuffer.data(idx_dip_di + 202);

    auto tr_y_xy_xxxyyz = pbuffer.data(idx_dip_di + 203);

    auto tr_y_xy_xxxyzz = pbuffer.data(idx_dip_di + 204);

    auto tr_y_xy_xxyyyy = pbuffer.data(idx_dip_di + 206);

    auto tr_y_xy_xxyyyz = pbuffer.data(idx_dip_di + 207);

    auto tr_y_xy_xxyyzz = pbuffer.data(idx_dip_di + 208);

    auto tr_y_xy_xxyzzz = pbuffer.data(idx_dip_di + 209);

    auto tr_y_xy_xyyyyy = pbuffer.data(idx_dip_di + 211);

    auto tr_y_xy_xyyyyz = pbuffer.data(idx_dip_di + 212);

    auto tr_y_xy_xyyyzz = pbuffer.data(idx_dip_di + 213);

    auto tr_y_xy_xyyzzz = pbuffer.data(idx_dip_di + 214);

    auto tr_y_xy_xyzzzz = pbuffer.data(idx_dip_di + 215);

    auto tr_y_xy_yyyyyy = pbuffer.data(idx_dip_di + 217);

    auto tr_y_xy_yyyyyz = pbuffer.data(idx_dip_di + 218);

    auto tr_y_xy_yyyyzz = pbuffer.data(idx_dip_di + 219);

    auto tr_y_xy_yyyzzz = pbuffer.data(idx_dip_di + 220);

    auto tr_y_xy_yyzzzz = pbuffer.data(idx_dip_di + 221);

    auto tr_y_xy_yzzzzz = pbuffer.data(idx_dip_di + 222);

    auto tr_y_xy_zzzzzz = pbuffer.data(idx_dip_di + 223);

    auto tr_y_xz_yyyyyz = pbuffer.data(idx_dip_di + 246);

    auto tr_y_xz_yyyyzz = pbuffer.data(idx_dip_di + 247);

    auto tr_y_xz_yyyzzz = pbuffer.data(idx_dip_di + 248);

    auto tr_y_xz_yyzzzz = pbuffer.data(idx_dip_di + 249);

    auto tr_y_xz_yzzzzz = pbuffer.data(idx_dip_di + 250);

    auto tr_y_xz_zzzzzz = pbuffer.data(idx_dip_di + 251);

    auto tr_y_yy_xxxxxx = pbuffer.data(idx_dip_di + 252);

    auto tr_y_yy_xxxxxy = pbuffer.data(idx_dip_di + 253);

    auto tr_y_yy_xxxxxz = pbuffer.data(idx_dip_di + 254);

    auto tr_y_yy_xxxxyy = pbuffer.data(idx_dip_di + 255);

    auto tr_y_yy_xxxxyz = pbuffer.data(idx_dip_di + 256);

    auto tr_y_yy_xxxxzz = pbuffer.data(idx_dip_di + 257);

    auto tr_y_yy_xxxyyy = pbuffer.data(idx_dip_di + 258);

    auto tr_y_yy_xxxyyz = pbuffer.data(idx_dip_di + 259);

    auto tr_y_yy_xxxyzz = pbuffer.data(idx_dip_di + 260);

    auto tr_y_yy_xxxzzz = pbuffer.data(idx_dip_di + 261);

    auto tr_y_yy_xxyyyy = pbuffer.data(idx_dip_di + 262);

    auto tr_y_yy_xxyyyz = pbuffer.data(idx_dip_di + 263);

    auto tr_y_yy_xxyyzz = pbuffer.data(idx_dip_di + 264);

    auto tr_y_yy_xxyzzz = pbuffer.data(idx_dip_di + 265);

    auto tr_y_yy_xxzzzz = pbuffer.data(idx_dip_di + 266);

    auto tr_y_yy_xyyyyy = pbuffer.data(idx_dip_di + 267);

    auto tr_y_yy_xyyyyz = pbuffer.data(idx_dip_di + 268);

    auto tr_y_yy_xyyyzz = pbuffer.data(idx_dip_di + 269);

    auto tr_y_yy_xyyzzz = pbuffer.data(idx_dip_di + 270);

    auto tr_y_yy_xyzzzz = pbuffer.data(idx_dip_di + 271);

    auto tr_y_yy_xzzzzz = pbuffer.data(idx_dip_di + 272);

    auto tr_y_yy_yyyyyy = pbuffer.data(idx_dip_di + 273);

    auto tr_y_yy_yyyyyz = pbuffer.data(idx_dip_di + 274);

    auto tr_y_yy_yyyyzz = pbuffer.data(idx_dip_di + 275);

    auto tr_y_yy_yyyzzz = pbuffer.data(idx_dip_di + 276);

    auto tr_y_yy_yyzzzz = pbuffer.data(idx_dip_di + 277);

    auto tr_y_yy_yzzzzz = pbuffer.data(idx_dip_di + 278);

    auto tr_y_yy_zzzzzz = pbuffer.data(idx_dip_di + 279);

    auto tr_y_yz_xxxxxy = pbuffer.data(idx_dip_di + 281);

    auto tr_y_yz_xxxxyy = pbuffer.data(idx_dip_di + 283);

    auto tr_y_yz_xxxyyy = pbuffer.data(idx_dip_di + 286);

    auto tr_y_yz_xxyyyy = pbuffer.data(idx_dip_di + 290);

    auto tr_y_yz_xyyyyy = pbuffer.data(idx_dip_di + 295);

    auto tr_y_yz_yyyyyy = pbuffer.data(idx_dip_di + 301);

    auto tr_y_yz_yyyyyz = pbuffer.data(idx_dip_di + 302);

    auto tr_y_yz_yyyyzz = pbuffer.data(idx_dip_di + 303);

    auto tr_y_yz_yyyzzz = pbuffer.data(idx_dip_di + 304);

    auto tr_y_yz_yyzzzz = pbuffer.data(idx_dip_di + 305);

    auto tr_y_yz_yzzzzz = pbuffer.data(idx_dip_di + 306);

    auto tr_y_yz_zzzzzz = pbuffer.data(idx_dip_di + 307);

    auto tr_y_zz_xxxxxx = pbuffer.data(idx_dip_di + 308);

    auto tr_y_zz_xxxxxy = pbuffer.data(idx_dip_di + 309);

    auto tr_y_zz_xxxxxz = pbuffer.data(idx_dip_di + 310);

    auto tr_y_zz_xxxxyy = pbuffer.data(idx_dip_di + 311);

    auto tr_y_zz_xxxxyz = pbuffer.data(idx_dip_di + 312);

    auto tr_y_zz_xxxxzz = pbuffer.data(idx_dip_di + 313);

    auto tr_y_zz_xxxyyy = pbuffer.data(idx_dip_di + 314);

    auto tr_y_zz_xxxyyz = pbuffer.data(idx_dip_di + 315);

    auto tr_y_zz_xxxyzz = pbuffer.data(idx_dip_di + 316);

    auto tr_y_zz_xxxzzz = pbuffer.data(idx_dip_di + 317);

    auto tr_y_zz_xxyyyy = pbuffer.data(idx_dip_di + 318);

    auto tr_y_zz_xxyyyz = pbuffer.data(idx_dip_di + 319);

    auto tr_y_zz_xxyyzz = pbuffer.data(idx_dip_di + 320);

    auto tr_y_zz_xxyzzz = pbuffer.data(idx_dip_di + 321);

    auto tr_y_zz_xxzzzz = pbuffer.data(idx_dip_di + 322);

    auto tr_y_zz_xyyyyy = pbuffer.data(idx_dip_di + 323);

    auto tr_y_zz_xyyyyz = pbuffer.data(idx_dip_di + 324);

    auto tr_y_zz_xyyyzz = pbuffer.data(idx_dip_di + 325);

    auto tr_y_zz_xyyzzz = pbuffer.data(idx_dip_di + 326);

    auto tr_y_zz_xyzzzz = pbuffer.data(idx_dip_di + 327);

    auto tr_y_zz_xzzzzz = pbuffer.data(idx_dip_di + 328);

    auto tr_y_zz_yyyyyy = pbuffer.data(idx_dip_di + 329);

    auto tr_y_zz_yyyyyz = pbuffer.data(idx_dip_di + 330);

    auto tr_y_zz_yyyyzz = pbuffer.data(idx_dip_di + 331);

    auto tr_y_zz_yyyzzz = pbuffer.data(idx_dip_di + 332);

    auto tr_y_zz_yyzzzz = pbuffer.data(idx_dip_di + 333);

    auto tr_y_zz_yzzzzz = pbuffer.data(idx_dip_di + 334);

    auto tr_y_zz_zzzzzz = pbuffer.data(idx_dip_di + 335);

    auto tr_z_xx_xxxxxx = pbuffer.data(idx_dip_di + 336);

    auto tr_z_xx_xxxxxy = pbuffer.data(idx_dip_di + 337);

    auto tr_z_xx_xxxxxz = pbuffer.data(idx_dip_di + 338);

    auto tr_z_xx_xxxxyy = pbuffer.data(idx_dip_di + 339);

    auto tr_z_xx_xxxxyz = pbuffer.data(idx_dip_di + 340);

    auto tr_z_xx_xxxxzz = pbuffer.data(idx_dip_di + 341);

    auto tr_z_xx_xxxyyy = pbuffer.data(idx_dip_di + 342);

    auto tr_z_xx_xxxyyz = pbuffer.data(idx_dip_di + 343);

    auto tr_z_xx_xxxyzz = pbuffer.data(idx_dip_di + 344);

    auto tr_z_xx_xxxzzz = pbuffer.data(idx_dip_di + 345);

    auto tr_z_xx_xxyyyy = pbuffer.data(idx_dip_di + 346);

    auto tr_z_xx_xxyyyz = pbuffer.data(idx_dip_di + 347);

    auto tr_z_xx_xxyyzz = pbuffer.data(idx_dip_di + 348);

    auto tr_z_xx_xxyzzz = pbuffer.data(idx_dip_di + 349);

    auto tr_z_xx_xxzzzz = pbuffer.data(idx_dip_di + 350);

    auto tr_z_xx_xyyyyy = pbuffer.data(idx_dip_di + 351);

    auto tr_z_xx_xyyyyz = pbuffer.data(idx_dip_di + 352);

    auto tr_z_xx_xyyyzz = pbuffer.data(idx_dip_di + 353);

    auto tr_z_xx_xyyzzz = pbuffer.data(idx_dip_di + 354);

    auto tr_z_xx_xyzzzz = pbuffer.data(idx_dip_di + 355);

    auto tr_z_xx_xzzzzz = pbuffer.data(idx_dip_di + 356);

    auto tr_z_xx_yyyyyy = pbuffer.data(idx_dip_di + 357);

    auto tr_z_xx_yyyyyz = pbuffer.data(idx_dip_di + 358);

    auto tr_z_xx_yyyyzz = pbuffer.data(idx_dip_di + 359);

    auto tr_z_xx_yyyzzz = pbuffer.data(idx_dip_di + 360);

    auto tr_z_xx_yyzzzz = pbuffer.data(idx_dip_di + 361);

    auto tr_z_xx_yzzzzz = pbuffer.data(idx_dip_di + 362);

    auto tr_z_xx_zzzzzz = pbuffer.data(idx_dip_di + 363);

    auto tr_z_xy_yyyyyy = pbuffer.data(idx_dip_di + 385);

    auto tr_z_xy_yyyyyz = pbuffer.data(idx_dip_di + 386);

    auto tr_z_xy_yyyyzz = pbuffer.data(idx_dip_di + 387);

    auto tr_z_xy_yyyzzz = pbuffer.data(idx_dip_di + 388);

    auto tr_z_xy_yyzzzz = pbuffer.data(idx_dip_di + 389);

    auto tr_z_xy_yzzzzz = pbuffer.data(idx_dip_di + 390);

    auto tr_z_xz_xxxxxz = pbuffer.data(idx_dip_di + 394);

    auto tr_z_xz_xxxxyz = pbuffer.data(idx_dip_di + 396);

    auto tr_z_xz_xxxxzz = pbuffer.data(idx_dip_di + 397);

    auto tr_z_xz_xxxyyz = pbuffer.data(idx_dip_di + 399);

    auto tr_z_xz_xxxyzz = pbuffer.data(idx_dip_di + 400);

    auto tr_z_xz_xxxzzz = pbuffer.data(idx_dip_di + 401);

    auto tr_z_xz_xxyyyz = pbuffer.data(idx_dip_di + 403);

    auto tr_z_xz_xxyyzz = pbuffer.data(idx_dip_di + 404);

    auto tr_z_xz_xxyzzz = pbuffer.data(idx_dip_di + 405);

    auto tr_z_xz_xxzzzz = pbuffer.data(idx_dip_di + 406);

    auto tr_z_xz_xyyyyz = pbuffer.data(idx_dip_di + 408);

    auto tr_z_xz_xyyyzz = pbuffer.data(idx_dip_di + 409);

    auto tr_z_xz_xyyzzz = pbuffer.data(idx_dip_di + 410);

    auto tr_z_xz_xyzzzz = pbuffer.data(idx_dip_di + 411);

    auto tr_z_xz_xzzzzz = pbuffer.data(idx_dip_di + 412);

    auto tr_z_xz_yyyyyy = pbuffer.data(idx_dip_di + 413);

    auto tr_z_xz_yyyyyz = pbuffer.data(idx_dip_di + 414);

    auto tr_z_xz_yyyyzz = pbuffer.data(idx_dip_di + 415);

    auto tr_z_xz_yyyzzz = pbuffer.data(idx_dip_di + 416);

    auto tr_z_xz_yyzzzz = pbuffer.data(idx_dip_di + 417);

    auto tr_z_xz_yzzzzz = pbuffer.data(idx_dip_di + 418);

    auto tr_z_xz_zzzzzz = pbuffer.data(idx_dip_di + 419);

    auto tr_z_yy_xxxxxx = pbuffer.data(idx_dip_di + 420);

    auto tr_z_yy_xxxxxy = pbuffer.data(idx_dip_di + 421);

    auto tr_z_yy_xxxxxz = pbuffer.data(idx_dip_di + 422);

    auto tr_z_yy_xxxxyy = pbuffer.data(idx_dip_di + 423);

    auto tr_z_yy_xxxxyz = pbuffer.data(idx_dip_di + 424);

    auto tr_z_yy_xxxxzz = pbuffer.data(idx_dip_di + 425);

    auto tr_z_yy_xxxyyy = pbuffer.data(idx_dip_di + 426);

    auto tr_z_yy_xxxyyz = pbuffer.data(idx_dip_di + 427);

    auto tr_z_yy_xxxyzz = pbuffer.data(idx_dip_di + 428);

    auto tr_z_yy_xxxzzz = pbuffer.data(idx_dip_di + 429);

    auto tr_z_yy_xxyyyy = pbuffer.data(idx_dip_di + 430);

    auto tr_z_yy_xxyyyz = pbuffer.data(idx_dip_di + 431);

    auto tr_z_yy_xxyyzz = pbuffer.data(idx_dip_di + 432);

    auto tr_z_yy_xxyzzz = pbuffer.data(idx_dip_di + 433);

    auto tr_z_yy_xxzzzz = pbuffer.data(idx_dip_di + 434);

    auto tr_z_yy_xyyyyy = pbuffer.data(idx_dip_di + 435);

    auto tr_z_yy_xyyyyz = pbuffer.data(idx_dip_di + 436);

    auto tr_z_yy_xyyyzz = pbuffer.data(idx_dip_di + 437);

    auto tr_z_yy_xyyzzz = pbuffer.data(idx_dip_di + 438);

    auto tr_z_yy_xyzzzz = pbuffer.data(idx_dip_di + 439);

    auto tr_z_yy_xzzzzz = pbuffer.data(idx_dip_di + 440);

    auto tr_z_yy_yyyyyy = pbuffer.data(idx_dip_di + 441);

    auto tr_z_yy_yyyyyz = pbuffer.data(idx_dip_di + 442);

    auto tr_z_yy_yyyyzz = pbuffer.data(idx_dip_di + 443);

    auto tr_z_yy_yyyzzz = pbuffer.data(idx_dip_di + 444);

    auto tr_z_yy_yyzzzz = pbuffer.data(idx_dip_di + 445);

    auto tr_z_yy_yzzzzz = pbuffer.data(idx_dip_di + 446);

    auto tr_z_yy_zzzzzz = pbuffer.data(idx_dip_di + 447);

    auto tr_z_yz_xxxxxx = pbuffer.data(idx_dip_di + 448);

    auto tr_z_yz_xxxxxz = pbuffer.data(idx_dip_di + 450);

    auto tr_z_yz_xxxxyz = pbuffer.data(idx_dip_di + 452);

    auto tr_z_yz_xxxxzz = pbuffer.data(idx_dip_di + 453);

    auto tr_z_yz_xxxyyz = pbuffer.data(idx_dip_di + 455);

    auto tr_z_yz_xxxyzz = pbuffer.data(idx_dip_di + 456);

    auto tr_z_yz_xxxzzz = pbuffer.data(idx_dip_di + 457);

    auto tr_z_yz_xxyyyz = pbuffer.data(idx_dip_di + 459);

    auto tr_z_yz_xxyyzz = pbuffer.data(idx_dip_di + 460);

    auto tr_z_yz_xxyzzz = pbuffer.data(idx_dip_di + 461);

    auto tr_z_yz_xxzzzz = pbuffer.data(idx_dip_di + 462);

    auto tr_z_yz_xyyyyz = pbuffer.data(idx_dip_di + 464);

    auto tr_z_yz_xyyyzz = pbuffer.data(idx_dip_di + 465);

    auto tr_z_yz_xyyzzz = pbuffer.data(idx_dip_di + 466);

    auto tr_z_yz_xyzzzz = pbuffer.data(idx_dip_di + 467);

    auto tr_z_yz_xzzzzz = pbuffer.data(idx_dip_di + 468);

    auto tr_z_yz_yyyyyy = pbuffer.data(idx_dip_di + 469);

    auto tr_z_yz_yyyyyz = pbuffer.data(idx_dip_di + 470);

    auto tr_z_yz_yyyyzz = pbuffer.data(idx_dip_di + 471);

    auto tr_z_yz_yyyzzz = pbuffer.data(idx_dip_di + 472);

    auto tr_z_yz_yyzzzz = pbuffer.data(idx_dip_di + 473);

    auto tr_z_yz_yzzzzz = pbuffer.data(idx_dip_di + 474);

    auto tr_z_yz_zzzzzz = pbuffer.data(idx_dip_di + 475);

    auto tr_z_zz_xxxxxx = pbuffer.data(idx_dip_di + 476);

    auto tr_z_zz_xxxxxy = pbuffer.data(idx_dip_di + 477);

    auto tr_z_zz_xxxxxz = pbuffer.data(idx_dip_di + 478);

    auto tr_z_zz_xxxxyy = pbuffer.data(idx_dip_di + 479);

    auto tr_z_zz_xxxxyz = pbuffer.data(idx_dip_di + 480);

    auto tr_z_zz_xxxxzz = pbuffer.data(idx_dip_di + 481);

    auto tr_z_zz_xxxyyy = pbuffer.data(idx_dip_di + 482);

    auto tr_z_zz_xxxyyz = pbuffer.data(idx_dip_di + 483);

    auto tr_z_zz_xxxyzz = pbuffer.data(idx_dip_di + 484);

    auto tr_z_zz_xxxzzz = pbuffer.data(idx_dip_di + 485);

    auto tr_z_zz_xxyyyy = pbuffer.data(idx_dip_di + 486);

    auto tr_z_zz_xxyyyz = pbuffer.data(idx_dip_di + 487);

    auto tr_z_zz_xxyyzz = pbuffer.data(idx_dip_di + 488);

    auto tr_z_zz_xxyzzz = pbuffer.data(idx_dip_di + 489);

    auto tr_z_zz_xxzzzz = pbuffer.data(idx_dip_di + 490);

    auto tr_z_zz_xyyyyy = pbuffer.data(idx_dip_di + 491);

    auto tr_z_zz_xyyyyz = pbuffer.data(idx_dip_di + 492);

    auto tr_z_zz_xyyyzz = pbuffer.data(idx_dip_di + 493);

    auto tr_z_zz_xyyzzz = pbuffer.data(idx_dip_di + 494);

    auto tr_z_zz_xyzzzz = pbuffer.data(idx_dip_di + 495);

    auto tr_z_zz_xzzzzz = pbuffer.data(idx_dip_di + 496);

    auto tr_z_zz_yyyyyy = pbuffer.data(idx_dip_di + 497);

    auto tr_z_zz_yyyyyz = pbuffer.data(idx_dip_di + 498);

    auto tr_z_zz_yyyyzz = pbuffer.data(idx_dip_di + 499);

    auto tr_z_zz_yyyzzz = pbuffer.data(idx_dip_di + 500);

    auto tr_z_zz_yyzzzz = pbuffer.data(idx_dip_di + 501);

    auto tr_z_zz_yzzzzz = pbuffer.data(idx_dip_di + 502);

    auto tr_z_zz_zzzzzz = pbuffer.data(idx_dip_di + 503);

    // Set up components of auxiliary buffer : FH

    auto tr_x_xxx_xxxxx = pbuffer.data(idx_dip_fh);

    auto tr_x_xxx_xxxxy = pbuffer.data(idx_dip_fh + 1);

    auto tr_x_xxx_xxxxz = pbuffer.data(idx_dip_fh + 2);

    auto tr_x_xxx_xxxyy = pbuffer.data(idx_dip_fh + 3);

    auto tr_x_xxx_xxxyz = pbuffer.data(idx_dip_fh + 4);

    auto tr_x_xxx_xxxzz = pbuffer.data(idx_dip_fh + 5);

    auto tr_x_xxx_xxyyy = pbuffer.data(idx_dip_fh + 6);

    auto tr_x_xxx_xxyyz = pbuffer.data(idx_dip_fh + 7);

    auto tr_x_xxx_xxyzz = pbuffer.data(idx_dip_fh + 8);

    auto tr_x_xxx_xxzzz = pbuffer.data(idx_dip_fh + 9);

    auto tr_x_xxx_xyyyy = pbuffer.data(idx_dip_fh + 10);

    auto tr_x_xxx_xyyyz = pbuffer.data(idx_dip_fh + 11);

    auto tr_x_xxx_xyyzz = pbuffer.data(idx_dip_fh + 12);

    auto tr_x_xxx_xyzzz = pbuffer.data(idx_dip_fh + 13);

    auto tr_x_xxx_xzzzz = pbuffer.data(idx_dip_fh + 14);

    auto tr_x_xxx_yyyyy = pbuffer.data(idx_dip_fh + 15);

    auto tr_x_xxx_yyyyz = pbuffer.data(idx_dip_fh + 16);

    auto tr_x_xxx_yyyzz = pbuffer.data(idx_dip_fh + 17);

    auto tr_x_xxx_yyzzz = pbuffer.data(idx_dip_fh + 18);

    auto tr_x_xxx_yzzzz = pbuffer.data(idx_dip_fh + 19);

    auto tr_x_xxx_zzzzz = pbuffer.data(idx_dip_fh + 20);

    auto tr_x_xxy_xxxxx = pbuffer.data(idx_dip_fh + 21);

    auto tr_x_xxy_xxxxy = pbuffer.data(idx_dip_fh + 22);

    auto tr_x_xxy_xxxxz = pbuffer.data(idx_dip_fh + 23);

    auto tr_x_xxy_xxxyy = pbuffer.data(idx_dip_fh + 24);

    auto tr_x_xxy_xxxyz = pbuffer.data(idx_dip_fh + 25);

    auto tr_x_xxy_xxxzz = pbuffer.data(idx_dip_fh + 26);

    auto tr_x_xxy_xxyyy = pbuffer.data(idx_dip_fh + 27);

    auto tr_x_xxy_xxyyz = pbuffer.data(idx_dip_fh + 28);

    auto tr_x_xxy_xxyzz = pbuffer.data(idx_dip_fh + 29);

    auto tr_x_xxy_xxzzz = pbuffer.data(idx_dip_fh + 30);

    auto tr_x_xxy_xyyyy = pbuffer.data(idx_dip_fh + 31);

    auto tr_x_xxy_xyyyz = pbuffer.data(idx_dip_fh + 32);

    auto tr_x_xxy_xyyzz = pbuffer.data(idx_dip_fh + 33);

    auto tr_x_xxy_xyzzz = pbuffer.data(idx_dip_fh + 34);

    auto tr_x_xxy_xzzzz = pbuffer.data(idx_dip_fh + 35);

    auto tr_x_xxz_xxxxx = pbuffer.data(idx_dip_fh + 42);

    auto tr_x_xxz_xxxxy = pbuffer.data(idx_dip_fh + 43);

    auto tr_x_xxz_xxxxz = pbuffer.data(idx_dip_fh + 44);

    auto tr_x_xxz_xxxyy = pbuffer.data(idx_dip_fh + 45);

    auto tr_x_xxz_xxxyz = pbuffer.data(idx_dip_fh + 46);

    auto tr_x_xxz_xxxzz = pbuffer.data(idx_dip_fh + 47);

    auto tr_x_xxz_xxyyy = pbuffer.data(idx_dip_fh + 48);

    auto tr_x_xxz_xxyyz = pbuffer.data(idx_dip_fh + 49);

    auto tr_x_xxz_xxyzz = pbuffer.data(idx_dip_fh + 50);

    auto tr_x_xxz_xxzzz = pbuffer.data(idx_dip_fh + 51);

    auto tr_x_xxz_xyyyy = pbuffer.data(idx_dip_fh + 52);

    auto tr_x_xxz_xyyyz = pbuffer.data(idx_dip_fh + 53);

    auto tr_x_xxz_xyyzz = pbuffer.data(idx_dip_fh + 54);

    auto tr_x_xxz_xyzzz = pbuffer.data(idx_dip_fh + 55);

    auto tr_x_xxz_xzzzz = pbuffer.data(idx_dip_fh + 56);

    auto tr_x_xxz_yyyyz = pbuffer.data(idx_dip_fh + 58);

    auto tr_x_xxz_yyyzz = pbuffer.data(idx_dip_fh + 59);

    auto tr_x_xxz_yyzzz = pbuffer.data(idx_dip_fh + 60);

    auto tr_x_xxz_yzzzz = pbuffer.data(idx_dip_fh + 61);

    auto tr_x_xxz_zzzzz = pbuffer.data(idx_dip_fh + 62);

    auto tr_x_xyy_xxxxy = pbuffer.data(idx_dip_fh + 64);

    auto tr_x_xyy_xxxyy = pbuffer.data(idx_dip_fh + 66);

    auto tr_x_xyy_xxxyz = pbuffer.data(idx_dip_fh + 67);

    auto tr_x_xyy_xxyyy = pbuffer.data(idx_dip_fh + 69);

    auto tr_x_xyy_xxyyz = pbuffer.data(idx_dip_fh + 70);

    auto tr_x_xyy_xxyzz = pbuffer.data(idx_dip_fh + 71);

    auto tr_x_xyy_xyyyy = pbuffer.data(idx_dip_fh + 73);

    auto tr_x_xyy_xyyyz = pbuffer.data(idx_dip_fh + 74);

    auto tr_x_xyy_xyyzz = pbuffer.data(idx_dip_fh + 75);

    auto tr_x_xyy_xyzzz = pbuffer.data(idx_dip_fh + 76);

    auto tr_x_xzz_xxxxx = pbuffer.data(idx_dip_fh + 105);

    auto tr_x_xzz_xxxxy = pbuffer.data(idx_dip_fh + 106);

    auto tr_x_xzz_xxxxz = pbuffer.data(idx_dip_fh + 107);

    auto tr_x_xzz_xxxyy = pbuffer.data(idx_dip_fh + 108);

    auto tr_x_xzz_xxxyz = pbuffer.data(idx_dip_fh + 109);

    auto tr_x_xzz_xxxzz = pbuffer.data(idx_dip_fh + 110);

    auto tr_x_xzz_xxyyy = pbuffer.data(idx_dip_fh + 111);

    auto tr_x_xzz_xxyyz = pbuffer.data(idx_dip_fh + 112);

    auto tr_x_xzz_xxyzz = pbuffer.data(idx_dip_fh + 113);

    auto tr_x_xzz_xxzzz = pbuffer.data(idx_dip_fh + 114);

    auto tr_x_xzz_xyyyy = pbuffer.data(idx_dip_fh + 115);

    auto tr_x_xzz_xyyyz = pbuffer.data(idx_dip_fh + 116);

    auto tr_x_xzz_xyyzz = pbuffer.data(idx_dip_fh + 117);

    auto tr_x_xzz_xyzzz = pbuffer.data(idx_dip_fh + 118);

    auto tr_x_xzz_xzzzz = pbuffer.data(idx_dip_fh + 119);

    auto tr_x_yyy_xxxxx = pbuffer.data(idx_dip_fh + 126);

    auto tr_x_yyy_xxxxy = pbuffer.data(idx_dip_fh + 127);

    auto tr_x_yyy_xxxxz = pbuffer.data(idx_dip_fh + 128);

    auto tr_x_yyy_xxxyy = pbuffer.data(idx_dip_fh + 129);

    auto tr_x_yyy_xxxyz = pbuffer.data(idx_dip_fh + 130);

    auto tr_x_yyy_xxxzz = pbuffer.data(idx_dip_fh + 131);

    auto tr_x_yyy_xxyyy = pbuffer.data(idx_dip_fh + 132);

    auto tr_x_yyy_xxyyz = pbuffer.data(idx_dip_fh + 133);

    auto tr_x_yyy_xxyzz = pbuffer.data(idx_dip_fh + 134);

    auto tr_x_yyy_xxzzz = pbuffer.data(idx_dip_fh + 135);

    auto tr_x_yyy_xyyyy = pbuffer.data(idx_dip_fh + 136);

    auto tr_x_yyy_xyyyz = pbuffer.data(idx_dip_fh + 137);

    auto tr_x_yyy_xyyzz = pbuffer.data(idx_dip_fh + 138);

    auto tr_x_yyy_xyzzz = pbuffer.data(idx_dip_fh + 139);

    auto tr_x_yyy_xzzzz = pbuffer.data(idx_dip_fh + 140);

    auto tr_x_yyy_yyyyy = pbuffer.data(idx_dip_fh + 141);

    auto tr_x_yyy_yyyyz = pbuffer.data(idx_dip_fh + 142);

    auto tr_x_yyy_yyyzz = pbuffer.data(idx_dip_fh + 143);

    auto tr_x_yyy_yyzzz = pbuffer.data(idx_dip_fh + 144);

    auto tr_x_yyy_yzzzz = pbuffer.data(idx_dip_fh + 145);

    auto tr_x_yyy_zzzzz = pbuffer.data(idx_dip_fh + 146);

    auto tr_x_yzz_xxxxz = pbuffer.data(idx_dip_fh + 170);

    auto tr_x_yzz_xxxyz = pbuffer.data(idx_dip_fh + 172);

    auto tr_x_yzz_xxxzz = pbuffer.data(idx_dip_fh + 173);

    auto tr_x_yzz_xxyyz = pbuffer.data(idx_dip_fh + 175);

    auto tr_x_yzz_xxyzz = pbuffer.data(idx_dip_fh + 176);

    auto tr_x_yzz_xxzzz = pbuffer.data(idx_dip_fh + 177);

    auto tr_x_yzz_xyyyz = pbuffer.data(idx_dip_fh + 179);

    auto tr_x_yzz_xyyzz = pbuffer.data(idx_dip_fh + 180);

    auto tr_x_yzz_xyzzz = pbuffer.data(idx_dip_fh + 181);

    auto tr_x_yzz_xzzzz = pbuffer.data(idx_dip_fh + 182);

    auto tr_x_yzz_yyyyz = pbuffer.data(idx_dip_fh + 184);

    auto tr_x_yzz_yyyzz = pbuffer.data(idx_dip_fh + 185);

    auto tr_x_yzz_yyzzz = pbuffer.data(idx_dip_fh + 186);

    auto tr_x_yzz_yzzzz = pbuffer.data(idx_dip_fh + 187);

    auto tr_x_yzz_zzzzz = pbuffer.data(idx_dip_fh + 188);

    auto tr_x_zzz_xxxxx = pbuffer.data(idx_dip_fh + 189);

    auto tr_x_zzz_xxxxy = pbuffer.data(idx_dip_fh + 190);

    auto tr_x_zzz_xxxxz = pbuffer.data(idx_dip_fh + 191);

    auto tr_x_zzz_xxxyy = pbuffer.data(idx_dip_fh + 192);

    auto tr_x_zzz_xxxyz = pbuffer.data(idx_dip_fh + 193);

    auto tr_x_zzz_xxxzz = pbuffer.data(idx_dip_fh + 194);

    auto tr_x_zzz_xxyyy = pbuffer.data(idx_dip_fh + 195);

    auto tr_x_zzz_xxyyz = pbuffer.data(idx_dip_fh + 196);

    auto tr_x_zzz_xxyzz = pbuffer.data(idx_dip_fh + 197);

    auto tr_x_zzz_xxzzz = pbuffer.data(idx_dip_fh + 198);

    auto tr_x_zzz_xyyyy = pbuffer.data(idx_dip_fh + 199);

    auto tr_x_zzz_xyyyz = pbuffer.data(idx_dip_fh + 200);

    auto tr_x_zzz_xyyzz = pbuffer.data(idx_dip_fh + 201);

    auto tr_x_zzz_xyzzz = pbuffer.data(idx_dip_fh + 202);

    auto tr_x_zzz_xzzzz = pbuffer.data(idx_dip_fh + 203);

    auto tr_x_zzz_yyyyy = pbuffer.data(idx_dip_fh + 204);

    auto tr_x_zzz_yyyyz = pbuffer.data(idx_dip_fh + 205);

    auto tr_x_zzz_yyyzz = pbuffer.data(idx_dip_fh + 206);

    auto tr_x_zzz_yyzzz = pbuffer.data(idx_dip_fh + 207);

    auto tr_x_zzz_yzzzz = pbuffer.data(idx_dip_fh + 208);

    auto tr_x_zzz_zzzzz = pbuffer.data(idx_dip_fh + 209);

    auto tr_y_xxx_xxxxx = pbuffer.data(idx_dip_fh + 210);

    auto tr_y_xxx_xxxxy = pbuffer.data(idx_dip_fh + 211);

    auto tr_y_xxx_xxxxz = pbuffer.data(idx_dip_fh + 212);

    auto tr_y_xxx_xxxyy = pbuffer.data(idx_dip_fh + 213);

    auto tr_y_xxx_xxxyz = pbuffer.data(idx_dip_fh + 214);

    auto tr_y_xxx_xxxzz = pbuffer.data(idx_dip_fh + 215);

    auto tr_y_xxx_xxyyy = pbuffer.data(idx_dip_fh + 216);

    auto tr_y_xxx_xxyyz = pbuffer.data(idx_dip_fh + 217);

    auto tr_y_xxx_xxyzz = pbuffer.data(idx_dip_fh + 218);

    auto tr_y_xxx_xxzzz = pbuffer.data(idx_dip_fh + 219);

    auto tr_y_xxx_xyyyy = pbuffer.data(idx_dip_fh + 220);

    auto tr_y_xxx_xyyyz = pbuffer.data(idx_dip_fh + 221);

    auto tr_y_xxx_xyyzz = pbuffer.data(idx_dip_fh + 222);

    auto tr_y_xxx_xyzzz = pbuffer.data(idx_dip_fh + 223);

    auto tr_y_xxx_xzzzz = pbuffer.data(idx_dip_fh + 224);

    auto tr_y_xxx_yyyyy = pbuffer.data(idx_dip_fh + 225);

    auto tr_y_xxx_yyyyz = pbuffer.data(idx_dip_fh + 226);

    auto tr_y_xxx_yyyzz = pbuffer.data(idx_dip_fh + 227);

    auto tr_y_xxx_yyzzz = pbuffer.data(idx_dip_fh + 228);

    auto tr_y_xxx_yzzzz = pbuffer.data(idx_dip_fh + 229);

    auto tr_y_xxx_zzzzz = pbuffer.data(idx_dip_fh + 230);

    auto tr_y_xxy_xxxxy = pbuffer.data(idx_dip_fh + 232);

    auto tr_y_xxy_xxxyy = pbuffer.data(idx_dip_fh + 234);

    auto tr_y_xxy_xxxyz = pbuffer.data(idx_dip_fh + 235);

    auto tr_y_xxy_xxyyy = pbuffer.data(idx_dip_fh + 237);

    auto tr_y_xxy_xxyyz = pbuffer.data(idx_dip_fh + 238);

    auto tr_y_xxy_xxyzz = pbuffer.data(idx_dip_fh + 239);

    auto tr_y_xxy_xyyyy = pbuffer.data(idx_dip_fh + 241);

    auto tr_y_xxy_xyyyz = pbuffer.data(idx_dip_fh + 242);

    auto tr_y_xxy_xyyzz = pbuffer.data(idx_dip_fh + 243);

    auto tr_y_xxy_xyzzz = pbuffer.data(idx_dip_fh + 244);

    auto tr_y_xxy_yyyyy = pbuffer.data(idx_dip_fh + 246);

    auto tr_y_xxy_yyyyz = pbuffer.data(idx_dip_fh + 247);

    auto tr_y_xxy_yyyzz = pbuffer.data(idx_dip_fh + 248);

    auto tr_y_xxy_yyzzz = pbuffer.data(idx_dip_fh + 249);

    auto tr_y_xxy_yzzzz = pbuffer.data(idx_dip_fh + 250);

    auto tr_y_xyy_xxxxx = pbuffer.data(idx_dip_fh + 273);

    auto tr_y_xyy_xxxxy = pbuffer.data(idx_dip_fh + 274);

    auto tr_y_xyy_xxxxz = pbuffer.data(idx_dip_fh + 275);

    auto tr_y_xyy_xxxyy = pbuffer.data(idx_dip_fh + 276);

    auto tr_y_xyy_xxxyz = pbuffer.data(idx_dip_fh + 277);

    auto tr_y_xyy_xxxzz = pbuffer.data(idx_dip_fh + 278);

    auto tr_y_xyy_xxyyy = pbuffer.data(idx_dip_fh + 279);

    auto tr_y_xyy_xxyyz = pbuffer.data(idx_dip_fh + 280);

    auto tr_y_xyy_xxyzz = pbuffer.data(idx_dip_fh + 281);

    auto tr_y_xyy_xxzzz = pbuffer.data(idx_dip_fh + 282);

    auto tr_y_xyy_xyyyy = pbuffer.data(idx_dip_fh + 283);

    auto tr_y_xyy_xyyyz = pbuffer.data(idx_dip_fh + 284);

    auto tr_y_xyy_xyyzz = pbuffer.data(idx_dip_fh + 285);

    auto tr_y_xyy_xyzzz = pbuffer.data(idx_dip_fh + 286);

    auto tr_y_xyy_xzzzz = pbuffer.data(idx_dip_fh + 287);

    auto tr_y_xyy_yyyyy = pbuffer.data(idx_dip_fh + 288);

    auto tr_y_xyy_yyyyz = pbuffer.data(idx_dip_fh + 289);

    auto tr_y_xyy_yyyzz = pbuffer.data(idx_dip_fh + 290);

    auto tr_y_xyy_yyzzz = pbuffer.data(idx_dip_fh + 291);

    auto tr_y_xyy_yzzzz = pbuffer.data(idx_dip_fh + 292);

    auto tr_y_xyy_zzzzz = pbuffer.data(idx_dip_fh + 293);

    auto tr_y_xzz_xxxxz = pbuffer.data(idx_dip_fh + 317);

    auto tr_y_xzz_xxxyz = pbuffer.data(idx_dip_fh + 319);

    auto tr_y_xzz_xxxzz = pbuffer.data(idx_dip_fh + 320);

    auto tr_y_xzz_xxyyz = pbuffer.data(idx_dip_fh + 322);

    auto tr_y_xzz_xxyzz = pbuffer.data(idx_dip_fh + 323);

    auto tr_y_xzz_xxzzz = pbuffer.data(idx_dip_fh + 324);

    auto tr_y_xzz_xyyyz = pbuffer.data(idx_dip_fh + 326);

    auto tr_y_xzz_xyyzz = pbuffer.data(idx_dip_fh + 327);

    auto tr_y_xzz_xyzzz = pbuffer.data(idx_dip_fh + 328);

    auto tr_y_xzz_xzzzz = pbuffer.data(idx_dip_fh + 329);

    auto tr_y_xzz_yyyyz = pbuffer.data(idx_dip_fh + 331);

    auto tr_y_xzz_yyyzz = pbuffer.data(idx_dip_fh + 332);

    auto tr_y_xzz_yyzzz = pbuffer.data(idx_dip_fh + 333);

    auto tr_y_xzz_yzzzz = pbuffer.data(idx_dip_fh + 334);

    auto tr_y_xzz_zzzzz = pbuffer.data(idx_dip_fh + 335);

    auto tr_y_yyy_xxxxx = pbuffer.data(idx_dip_fh + 336);

    auto tr_y_yyy_xxxxy = pbuffer.data(idx_dip_fh + 337);

    auto tr_y_yyy_xxxxz = pbuffer.data(idx_dip_fh + 338);

    auto tr_y_yyy_xxxyy = pbuffer.data(idx_dip_fh + 339);

    auto tr_y_yyy_xxxyz = pbuffer.data(idx_dip_fh + 340);

    auto tr_y_yyy_xxxzz = pbuffer.data(idx_dip_fh + 341);

    auto tr_y_yyy_xxyyy = pbuffer.data(idx_dip_fh + 342);

    auto tr_y_yyy_xxyyz = pbuffer.data(idx_dip_fh + 343);

    auto tr_y_yyy_xxyzz = pbuffer.data(idx_dip_fh + 344);

    auto tr_y_yyy_xxzzz = pbuffer.data(idx_dip_fh + 345);

    auto tr_y_yyy_xyyyy = pbuffer.data(idx_dip_fh + 346);

    auto tr_y_yyy_xyyyz = pbuffer.data(idx_dip_fh + 347);

    auto tr_y_yyy_xyyzz = pbuffer.data(idx_dip_fh + 348);

    auto tr_y_yyy_xyzzz = pbuffer.data(idx_dip_fh + 349);

    auto tr_y_yyy_xzzzz = pbuffer.data(idx_dip_fh + 350);

    auto tr_y_yyy_yyyyy = pbuffer.data(idx_dip_fh + 351);

    auto tr_y_yyy_yyyyz = pbuffer.data(idx_dip_fh + 352);

    auto tr_y_yyy_yyyzz = pbuffer.data(idx_dip_fh + 353);

    auto tr_y_yyy_yyzzz = pbuffer.data(idx_dip_fh + 354);

    auto tr_y_yyy_yzzzz = pbuffer.data(idx_dip_fh + 355);

    auto tr_y_yyy_zzzzz = pbuffer.data(idx_dip_fh + 356);

    auto tr_y_yyz_xxxxy = pbuffer.data(idx_dip_fh + 358);

    auto tr_y_yyz_xxxxz = pbuffer.data(idx_dip_fh + 359);

    auto tr_y_yyz_xxxyy = pbuffer.data(idx_dip_fh + 360);

    auto tr_y_yyz_xxxyz = pbuffer.data(idx_dip_fh + 361);

    auto tr_y_yyz_xxxzz = pbuffer.data(idx_dip_fh + 362);

    auto tr_y_yyz_xxyyy = pbuffer.data(idx_dip_fh + 363);

    auto tr_y_yyz_xxyyz = pbuffer.data(idx_dip_fh + 364);

    auto tr_y_yyz_xxyzz = pbuffer.data(idx_dip_fh + 365);

    auto tr_y_yyz_xxzzz = pbuffer.data(idx_dip_fh + 366);

    auto tr_y_yyz_xyyyy = pbuffer.data(idx_dip_fh + 367);

    auto tr_y_yyz_xyyyz = pbuffer.data(idx_dip_fh + 368);

    auto tr_y_yyz_xyyzz = pbuffer.data(idx_dip_fh + 369);

    auto tr_y_yyz_xyzzz = pbuffer.data(idx_dip_fh + 370);

    auto tr_y_yyz_xzzzz = pbuffer.data(idx_dip_fh + 371);

    auto tr_y_yyz_yyyyy = pbuffer.data(idx_dip_fh + 372);

    auto tr_y_yyz_yyyyz = pbuffer.data(idx_dip_fh + 373);

    auto tr_y_yyz_yyyzz = pbuffer.data(idx_dip_fh + 374);

    auto tr_y_yyz_yyzzz = pbuffer.data(idx_dip_fh + 375);

    auto tr_y_yyz_yzzzz = pbuffer.data(idx_dip_fh + 376);

    auto tr_y_yyz_zzzzz = pbuffer.data(idx_dip_fh + 377);

    auto tr_y_yzz_xxxxx = pbuffer.data(idx_dip_fh + 378);

    auto tr_y_yzz_xxxxy = pbuffer.data(idx_dip_fh + 379);

    auto tr_y_yzz_xxxxz = pbuffer.data(idx_dip_fh + 380);

    auto tr_y_yzz_xxxyy = pbuffer.data(idx_dip_fh + 381);

    auto tr_y_yzz_xxxyz = pbuffer.data(idx_dip_fh + 382);

    auto tr_y_yzz_xxxzz = pbuffer.data(idx_dip_fh + 383);

    auto tr_y_yzz_xxyyy = pbuffer.data(idx_dip_fh + 384);

    auto tr_y_yzz_xxyyz = pbuffer.data(idx_dip_fh + 385);

    auto tr_y_yzz_xxyzz = pbuffer.data(idx_dip_fh + 386);

    auto tr_y_yzz_xxzzz = pbuffer.data(idx_dip_fh + 387);

    auto tr_y_yzz_xyyyy = pbuffer.data(idx_dip_fh + 388);

    auto tr_y_yzz_xyyyz = pbuffer.data(idx_dip_fh + 389);

    auto tr_y_yzz_xyyzz = pbuffer.data(idx_dip_fh + 390);

    auto tr_y_yzz_xyzzz = pbuffer.data(idx_dip_fh + 391);

    auto tr_y_yzz_xzzzz = pbuffer.data(idx_dip_fh + 392);

    auto tr_y_yzz_yyyyy = pbuffer.data(idx_dip_fh + 393);

    auto tr_y_yzz_yyyyz = pbuffer.data(idx_dip_fh + 394);

    auto tr_y_yzz_yyyzz = pbuffer.data(idx_dip_fh + 395);

    auto tr_y_yzz_yyzzz = pbuffer.data(idx_dip_fh + 396);

    auto tr_y_yzz_yzzzz = pbuffer.data(idx_dip_fh + 397);

    auto tr_y_yzz_zzzzz = pbuffer.data(idx_dip_fh + 398);

    auto tr_y_zzz_xxxxx = pbuffer.data(idx_dip_fh + 399);

    auto tr_y_zzz_xxxxy = pbuffer.data(idx_dip_fh + 400);

    auto tr_y_zzz_xxxxz = pbuffer.data(idx_dip_fh + 401);

    auto tr_y_zzz_xxxyy = pbuffer.data(idx_dip_fh + 402);

    auto tr_y_zzz_xxxyz = pbuffer.data(idx_dip_fh + 403);

    auto tr_y_zzz_xxxzz = pbuffer.data(idx_dip_fh + 404);

    auto tr_y_zzz_xxyyy = pbuffer.data(idx_dip_fh + 405);

    auto tr_y_zzz_xxyyz = pbuffer.data(idx_dip_fh + 406);

    auto tr_y_zzz_xxyzz = pbuffer.data(idx_dip_fh + 407);

    auto tr_y_zzz_xxzzz = pbuffer.data(idx_dip_fh + 408);

    auto tr_y_zzz_xyyyy = pbuffer.data(idx_dip_fh + 409);

    auto tr_y_zzz_xyyyz = pbuffer.data(idx_dip_fh + 410);

    auto tr_y_zzz_xyyzz = pbuffer.data(idx_dip_fh + 411);

    auto tr_y_zzz_xyzzz = pbuffer.data(idx_dip_fh + 412);

    auto tr_y_zzz_xzzzz = pbuffer.data(idx_dip_fh + 413);

    auto tr_y_zzz_yyyyy = pbuffer.data(idx_dip_fh + 414);

    auto tr_y_zzz_yyyyz = pbuffer.data(idx_dip_fh + 415);

    auto tr_y_zzz_yyyzz = pbuffer.data(idx_dip_fh + 416);

    auto tr_y_zzz_yyzzz = pbuffer.data(idx_dip_fh + 417);

    auto tr_y_zzz_yzzzz = pbuffer.data(idx_dip_fh + 418);

    auto tr_y_zzz_zzzzz = pbuffer.data(idx_dip_fh + 419);

    auto tr_z_xxx_xxxxx = pbuffer.data(idx_dip_fh + 420);

    auto tr_z_xxx_xxxxy = pbuffer.data(idx_dip_fh + 421);

    auto tr_z_xxx_xxxxz = pbuffer.data(idx_dip_fh + 422);

    auto tr_z_xxx_xxxyy = pbuffer.data(idx_dip_fh + 423);

    auto tr_z_xxx_xxxyz = pbuffer.data(idx_dip_fh + 424);

    auto tr_z_xxx_xxxzz = pbuffer.data(idx_dip_fh + 425);

    auto tr_z_xxx_xxyyy = pbuffer.data(idx_dip_fh + 426);

    auto tr_z_xxx_xxyyz = pbuffer.data(idx_dip_fh + 427);

    auto tr_z_xxx_xxyzz = pbuffer.data(idx_dip_fh + 428);

    auto tr_z_xxx_xxzzz = pbuffer.data(idx_dip_fh + 429);

    auto tr_z_xxx_xyyyy = pbuffer.data(idx_dip_fh + 430);

    auto tr_z_xxx_xyyyz = pbuffer.data(idx_dip_fh + 431);

    auto tr_z_xxx_xyyzz = pbuffer.data(idx_dip_fh + 432);

    auto tr_z_xxx_xyzzz = pbuffer.data(idx_dip_fh + 433);

    auto tr_z_xxx_xzzzz = pbuffer.data(idx_dip_fh + 434);

    auto tr_z_xxx_yyyyy = pbuffer.data(idx_dip_fh + 435);

    auto tr_z_xxx_yyyyz = pbuffer.data(idx_dip_fh + 436);

    auto tr_z_xxx_yyyzz = pbuffer.data(idx_dip_fh + 437);

    auto tr_z_xxx_yyzzz = pbuffer.data(idx_dip_fh + 438);

    auto tr_z_xxx_yzzzz = pbuffer.data(idx_dip_fh + 439);

    auto tr_z_xxx_zzzzz = pbuffer.data(idx_dip_fh + 440);

    auto tr_z_xxz_xxxxx = pbuffer.data(idx_dip_fh + 462);

    auto tr_z_xxz_xxxxy = pbuffer.data(idx_dip_fh + 463);

    auto tr_z_xxz_xxxxz = pbuffer.data(idx_dip_fh + 464);

    auto tr_z_xxz_xxxyy = pbuffer.data(idx_dip_fh + 465);

    auto tr_z_xxz_xxxyz = pbuffer.data(idx_dip_fh + 466);

    auto tr_z_xxz_xxxzz = pbuffer.data(idx_dip_fh + 467);

    auto tr_z_xxz_xxyyy = pbuffer.data(idx_dip_fh + 468);

    auto tr_z_xxz_xxyyz = pbuffer.data(idx_dip_fh + 469);

    auto tr_z_xxz_xxyzz = pbuffer.data(idx_dip_fh + 470);

    auto tr_z_xxz_xxzzz = pbuffer.data(idx_dip_fh + 471);

    auto tr_z_xxz_xyyyy = pbuffer.data(idx_dip_fh + 472);

    auto tr_z_xxz_xyyyz = pbuffer.data(idx_dip_fh + 473);

    auto tr_z_xxz_xyyzz = pbuffer.data(idx_dip_fh + 474);

    auto tr_z_xxz_xyzzz = pbuffer.data(idx_dip_fh + 475);

    auto tr_z_xxz_xzzzz = pbuffer.data(idx_dip_fh + 476);

    auto tr_z_xxz_yyyyz = pbuffer.data(idx_dip_fh + 478);

    auto tr_z_xxz_yyyzz = pbuffer.data(idx_dip_fh + 479);

    auto tr_z_xxz_yyzzz = pbuffer.data(idx_dip_fh + 480);

    auto tr_z_xxz_yzzzz = pbuffer.data(idx_dip_fh + 481);

    auto tr_z_xxz_zzzzz = pbuffer.data(idx_dip_fh + 482);

    auto tr_z_xyy_xxxxy = pbuffer.data(idx_dip_fh + 484);

    auto tr_z_xyy_xxxyy = pbuffer.data(idx_dip_fh + 486);

    auto tr_z_xyy_xxxyz = pbuffer.data(idx_dip_fh + 487);

    auto tr_z_xyy_xxyyy = pbuffer.data(idx_dip_fh + 489);

    auto tr_z_xyy_xxyyz = pbuffer.data(idx_dip_fh + 490);

    auto tr_z_xyy_xxyzz = pbuffer.data(idx_dip_fh + 491);

    auto tr_z_xyy_xyyyy = pbuffer.data(idx_dip_fh + 493);

    auto tr_z_xyy_xyyyz = pbuffer.data(idx_dip_fh + 494);

    auto tr_z_xyy_xyyzz = pbuffer.data(idx_dip_fh + 495);

    auto tr_z_xyy_xyzzz = pbuffer.data(idx_dip_fh + 496);

    auto tr_z_xyy_yyyyy = pbuffer.data(idx_dip_fh + 498);

    auto tr_z_xyy_yyyyz = pbuffer.data(idx_dip_fh + 499);

    auto tr_z_xyy_yyyzz = pbuffer.data(idx_dip_fh + 500);

    auto tr_z_xyy_yyzzz = pbuffer.data(idx_dip_fh + 501);

    auto tr_z_xyy_yzzzz = pbuffer.data(idx_dip_fh + 502);

    auto tr_z_xzz_xxxxx = pbuffer.data(idx_dip_fh + 525);

    auto tr_z_xzz_xxxxy = pbuffer.data(idx_dip_fh + 526);

    auto tr_z_xzz_xxxxz = pbuffer.data(idx_dip_fh + 527);

    auto tr_z_xzz_xxxyy = pbuffer.data(idx_dip_fh + 528);

    auto tr_z_xzz_xxxyz = pbuffer.data(idx_dip_fh + 529);

    auto tr_z_xzz_xxxzz = pbuffer.data(idx_dip_fh + 530);

    auto tr_z_xzz_xxyyy = pbuffer.data(idx_dip_fh + 531);

    auto tr_z_xzz_xxyyz = pbuffer.data(idx_dip_fh + 532);

    auto tr_z_xzz_xxyzz = pbuffer.data(idx_dip_fh + 533);

    auto tr_z_xzz_xxzzz = pbuffer.data(idx_dip_fh + 534);

    auto tr_z_xzz_xyyyy = pbuffer.data(idx_dip_fh + 535);

    auto tr_z_xzz_xyyyz = pbuffer.data(idx_dip_fh + 536);

    auto tr_z_xzz_xyyzz = pbuffer.data(idx_dip_fh + 537);

    auto tr_z_xzz_xyzzz = pbuffer.data(idx_dip_fh + 538);

    auto tr_z_xzz_xzzzz = pbuffer.data(idx_dip_fh + 539);

    auto tr_z_xzz_yyyyy = pbuffer.data(idx_dip_fh + 540);

    auto tr_z_xzz_yyyyz = pbuffer.data(idx_dip_fh + 541);

    auto tr_z_xzz_yyyzz = pbuffer.data(idx_dip_fh + 542);

    auto tr_z_xzz_yyzzz = pbuffer.data(idx_dip_fh + 543);

    auto tr_z_xzz_yzzzz = pbuffer.data(idx_dip_fh + 544);

    auto tr_z_xzz_zzzzz = pbuffer.data(idx_dip_fh + 545);

    auto tr_z_yyy_xxxxx = pbuffer.data(idx_dip_fh + 546);

    auto tr_z_yyy_xxxxy = pbuffer.data(idx_dip_fh + 547);

    auto tr_z_yyy_xxxxz = pbuffer.data(idx_dip_fh + 548);

    auto tr_z_yyy_xxxyy = pbuffer.data(idx_dip_fh + 549);

    auto tr_z_yyy_xxxyz = pbuffer.data(idx_dip_fh + 550);

    auto tr_z_yyy_xxxzz = pbuffer.data(idx_dip_fh + 551);

    auto tr_z_yyy_xxyyy = pbuffer.data(idx_dip_fh + 552);

    auto tr_z_yyy_xxyyz = pbuffer.data(idx_dip_fh + 553);

    auto tr_z_yyy_xxyzz = pbuffer.data(idx_dip_fh + 554);

    auto tr_z_yyy_xxzzz = pbuffer.data(idx_dip_fh + 555);

    auto tr_z_yyy_xyyyy = pbuffer.data(idx_dip_fh + 556);

    auto tr_z_yyy_xyyyz = pbuffer.data(idx_dip_fh + 557);

    auto tr_z_yyy_xyyzz = pbuffer.data(idx_dip_fh + 558);

    auto tr_z_yyy_xyzzz = pbuffer.data(idx_dip_fh + 559);

    auto tr_z_yyy_xzzzz = pbuffer.data(idx_dip_fh + 560);

    auto tr_z_yyy_yyyyy = pbuffer.data(idx_dip_fh + 561);

    auto tr_z_yyy_yyyyz = pbuffer.data(idx_dip_fh + 562);

    auto tr_z_yyy_yyyzz = pbuffer.data(idx_dip_fh + 563);

    auto tr_z_yyy_yyzzz = pbuffer.data(idx_dip_fh + 564);

    auto tr_z_yyy_yzzzz = pbuffer.data(idx_dip_fh + 565);

    auto tr_z_yyy_zzzzz = pbuffer.data(idx_dip_fh + 566);

    auto tr_z_yyz_xxxxx = pbuffer.data(idx_dip_fh + 567);

    auto tr_z_yyz_xxxxy = pbuffer.data(idx_dip_fh + 568);

    auto tr_z_yyz_xxxxz = pbuffer.data(idx_dip_fh + 569);

    auto tr_z_yyz_xxxyy = pbuffer.data(idx_dip_fh + 570);

    auto tr_z_yyz_xxxyz = pbuffer.data(idx_dip_fh + 571);

    auto tr_z_yyz_xxxzz = pbuffer.data(idx_dip_fh + 572);

    auto tr_z_yyz_xxyyy = pbuffer.data(idx_dip_fh + 573);

    auto tr_z_yyz_xxyyz = pbuffer.data(idx_dip_fh + 574);

    auto tr_z_yyz_xxyzz = pbuffer.data(idx_dip_fh + 575);

    auto tr_z_yyz_xxzzz = pbuffer.data(idx_dip_fh + 576);

    auto tr_z_yyz_xyyyy = pbuffer.data(idx_dip_fh + 577);

    auto tr_z_yyz_xyyyz = pbuffer.data(idx_dip_fh + 578);

    auto tr_z_yyz_xyyzz = pbuffer.data(idx_dip_fh + 579);

    auto tr_z_yyz_xyzzz = pbuffer.data(idx_dip_fh + 580);

    auto tr_z_yyz_xzzzz = pbuffer.data(idx_dip_fh + 581);

    auto tr_z_yyz_yyyyy = pbuffer.data(idx_dip_fh + 582);

    auto tr_z_yyz_yyyyz = pbuffer.data(idx_dip_fh + 583);

    auto tr_z_yyz_yyyzz = pbuffer.data(idx_dip_fh + 584);

    auto tr_z_yyz_yyzzz = pbuffer.data(idx_dip_fh + 585);

    auto tr_z_yyz_yzzzz = pbuffer.data(idx_dip_fh + 586);

    auto tr_z_yyz_zzzzz = pbuffer.data(idx_dip_fh + 587);

    auto tr_z_yzz_xxxxx = pbuffer.data(idx_dip_fh + 588);

    auto tr_z_yzz_xxxxy = pbuffer.data(idx_dip_fh + 589);

    auto tr_z_yzz_xxxxz = pbuffer.data(idx_dip_fh + 590);

    auto tr_z_yzz_xxxyy = pbuffer.data(idx_dip_fh + 591);

    auto tr_z_yzz_xxxyz = pbuffer.data(idx_dip_fh + 592);

    auto tr_z_yzz_xxxzz = pbuffer.data(idx_dip_fh + 593);

    auto tr_z_yzz_xxyyy = pbuffer.data(idx_dip_fh + 594);

    auto tr_z_yzz_xxyyz = pbuffer.data(idx_dip_fh + 595);

    auto tr_z_yzz_xxyzz = pbuffer.data(idx_dip_fh + 596);

    auto tr_z_yzz_xxzzz = pbuffer.data(idx_dip_fh + 597);

    auto tr_z_yzz_xyyyy = pbuffer.data(idx_dip_fh + 598);

    auto tr_z_yzz_xyyyz = pbuffer.data(idx_dip_fh + 599);

    auto tr_z_yzz_xyyzz = pbuffer.data(idx_dip_fh + 600);

    auto tr_z_yzz_xyzzz = pbuffer.data(idx_dip_fh + 601);

    auto tr_z_yzz_xzzzz = pbuffer.data(idx_dip_fh + 602);

    auto tr_z_yzz_yyyyy = pbuffer.data(idx_dip_fh + 603);

    auto tr_z_yzz_yyyyz = pbuffer.data(idx_dip_fh + 604);

    auto tr_z_yzz_yyyzz = pbuffer.data(idx_dip_fh + 605);

    auto tr_z_yzz_yyzzz = pbuffer.data(idx_dip_fh + 606);

    auto tr_z_yzz_yzzzz = pbuffer.data(idx_dip_fh + 607);

    auto tr_z_yzz_zzzzz = pbuffer.data(idx_dip_fh + 608);

    auto tr_z_zzz_xxxxx = pbuffer.data(idx_dip_fh + 609);

    auto tr_z_zzz_xxxxy = pbuffer.data(idx_dip_fh + 610);

    auto tr_z_zzz_xxxxz = pbuffer.data(idx_dip_fh + 611);

    auto tr_z_zzz_xxxyy = pbuffer.data(idx_dip_fh + 612);

    auto tr_z_zzz_xxxyz = pbuffer.data(idx_dip_fh + 613);

    auto tr_z_zzz_xxxzz = pbuffer.data(idx_dip_fh + 614);

    auto tr_z_zzz_xxyyy = pbuffer.data(idx_dip_fh + 615);

    auto tr_z_zzz_xxyyz = pbuffer.data(idx_dip_fh + 616);

    auto tr_z_zzz_xxyzz = pbuffer.data(idx_dip_fh + 617);

    auto tr_z_zzz_xxzzz = pbuffer.data(idx_dip_fh + 618);

    auto tr_z_zzz_xyyyy = pbuffer.data(idx_dip_fh + 619);

    auto tr_z_zzz_xyyyz = pbuffer.data(idx_dip_fh + 620);

    auto tr_z_zzz_xyyzz = pbuffer.data(idx_dip_fh + 621);

    auto tr_z_zzz_xyzzz = pbuffer.data(idx_dip_fh + 622);

    auto tr_z_zzz_xzzzz = pbuffer.data(idx_dip_fh + 623);

    auto tr_z_zzz_yyyyy = pbuffer.data(idx_dip_fh + 624);

    auto tr_z_zzz_yyyyz = pbuffer.data(idx_dip_fh + 625);

    auto tr_z_zzz_yyyzz = pbuffer.data(idx_dip_fh + 626);

    auto tr_z_zzz_yyzzz = pbuffer.data(idx_dip_fh + 627);

    auto tr_z_zzz_yzzzz = pbuffer.data(idx_dip_fh + 628);

    auto tr_z_zzz_zzzzz = pbuffer.data(idx_dip_fh + 629);

    // Set up components of auxiliary buffer : FI

    auto ts_xxx_xxxxxx = pbuffer.data(idx_ovl_fi);

    auto ts_xxx_xxxxxy = pbuffer.data(idx_ovl_fi + 1);

    auto ts_xxx_xxxxxz = pbuffer.data(idx_ovl_fi + 2);

    auto ts_xxx_xxxxyy = pbuffer.data(idx_ovl_fi + 3);

    auto ts_xxx_xxxxyz = pbuffer.data(idx_ovl_fi + 4);

    auto ts_xxx_xxxxzz = pbuffer.data(idx_ovl_fi + 5);

    auto ts_xxx_xxxyyy = pbuffer.data(idx_ovl_fi + 6);

    auto ts_xxx_xxxyyz = pbuffer.data(idx_ovl_fi + 7);

    auto ts_xxx_xxxyzz = pbuffer.data(idx_ovl_fi + 8);

    auto ts_xxx_xxxzzz = pbuffer.data(idx_ovl_fi + 9);

    auto ts_xxx_xxyyyy = pbuffer.data(idx_ovl_fi + 10);

    auto ts_xxx_xxyyyz = pbuffer.data(idx_ovl_fi + 11);

    auto ts_xxx_xxyyzz = pbuffer.data(idx_ovl_fi + 12);

    auto ts_xxx_xxyzzz = pbuffer.data(idx_ovl_fi + 13);

    auto ts_xxx_xxzzzz = pbuffer.data(idx_ovl_fi + 14);

    auto ts_xxx_xyyyyy = pbuffer.data(idx_ovl_fi + 15);

    auto ts_xxx_xyyyyz = pbuffer.data(idx_ovl_fi + 16);

    auto ts_xxx_xyyyzz = pbuffer.data(idx_ovl_fi + 17);

    auto ts_xxx_xyyzzz = pbuffer.data(idx_ovl_fi + 18);

    auto ts_xxx_xyzzzz = pbuffer.data(idx_ovl_fi + 19);

    auto ts_xxx_xzzzzz = pbuffer.data(idx_ovl_fi + 20);

    auto ts_xxx_yyyyyy = pbuffer.data(idx_ovl_fi + 21);

    auto ts_xxx_yyyyyz = pbuffer.data(idx_ovl_fi + 22);

    auto ts_xxx_yyyyzz = pbuffer.data(idx_ovl_fi + 23);

    auto ts_xxx_yyyzzz = pbuffer.data(idx_ovl_fi + 24);

    auto ts_xxx_yyzzzz = pbuffer.data(idx_ovl_fi + 25);

    auto ts_xxx_yzzzzz = pbuffer.data(idx_ovl_fi + 26);

    auto ts_xxx_zzzzzz = pbuffer.data(idx_ovl_fi + 27);

    auto ts_xxz_xxxxxz = pbuffer.data(idx_ovl_fi + 58);

    auto ts_xxz_xxxxzz = pbuffer.data(idx_ovl_fi + 61);

    auto ts_xxz_xxxzzz = pbuffer.data(idx_ovl_fi + 65);

    auto ts_xxz_xxzzzz = pbuffer.data(idx_ovl_fi + 70);

    auto ts_xxz_xzzzzz = pbuffer.data(idx_ovl_fi + 76);

    auto ts_xyy_yyyyyy = pbuffer.data(idx_ovl_fi + 105);

    auto ts_xyy_yyyyyz = pbuffer.data(idx_ovl_fi + 106);

    auto ts_xyy_yyyyzz = pbuffer.data(idx_ovl_fi + 107);

    auto ts_xyy_yyyzzz = pbuffer.data(idx_ovl_fi + 108);

    auto ts_xyy_yyzzzz = pbuffer.data(idx_ovl_fi + 109);

    auto ts_xyy_yzzzzz = pbuffer.data(idx_ovl_fi + 110);

    auto ts_xzz_yyyyyz = pbuffer.data(idx_ovl_fi + 162);

    auto ts_xzz_yyyyzz = pbuffer.data(idx_ovl_fi + 163);

    auto ts_xzz_yyyzzz = pbuffer.data(idx_ovl_fi + 164);

    auto ts_xzz_yyzzzz = pbuffer.data(idx_ovl_fi + 165);

    auto ts_xzz_yzzzzz = pbuffer.data(idx_ovl_fi + 166);

    auto ts_xzz_zzzzzz = pbuffer.data(idx_ovl_fi + 167);

    auto ts_yyy_xxxxxx = pbuffer.data(idx_ovl_fi + 168);

    auto ts_yyy_xxxxxy = pbuffer.data(idx_ovl_fi + 169);

    auto ts_yyy_xxxxxz = pbuffer.data(idx_ovl_fi + 170);

    auto ts_yyy_xxxxyy = pbuffer.data(idx_ovl_fi + 171);

    auto ts_yyy_xxxxyz = pbuffer.data(idx_ovl_fi + 172);

    auto ts_yyy_xxxxzz = pbuffer.data(idx_ovl_fi + 173);

    auto ts_yyy_xxxyyy = pbuffer.data(idx_ovl_fi + 174);

    auto ts_yyy_xxxyyz = pbuffer.data(idx_ovl_fi + 175);

    auto ts_yyy_xxxyzz = pbuffer.data(idx_ovl_fi + 176);

    auto ts_yyy_xxxzzz = pbuffer.data(idx_ovl_fi + 177);

    auto ts_yyy_xxyyyy = pbuffer.data(idx_ovl_fi + 178);

    auto ts_yyy_xxyyyz = pbuffer.data(idx_ovl_fi + 179);

    auto ts_yyy_xxyyzz = pbuffer.data(idx_ovl_fi + 180);

    auto ts_yyy_xxyzzz = pbuffer.data(idx_ovl_fi + 181);

    auto ts_yyy_xxzzzz = pbuffer.data(idx_ovl_fi + 182);

    auto ts_yyy_xyyyyy = pbuffer.data(idx_ovl_fi + 183);

    auto ts_yyy_xyyyyz = pbuffer.data(idx_ovl_fi + 184);

    auto ts_yyy_xyyyzz = pbuffer.data(idx_ovl_fi + 185);

    auto ts_yyy_xyyzzz = pbuffer.data(idx_ovl_fi + 186);

    auto ts_yyy_xyzzzz = pbuffer.data(idx_ovl_fi + 187);

    auto ts_yyy_xzzzzz = pbuffer.data(idx_ovl_fi + 188);

    auto ts_yyy_yyyyyy = pbuffer.data(idx_ovl_fi + 189);

    auto ts_yyy_yyyyyz = pbuffer.data(idx_ovl_fi + 190);

    auto ts_yyy_yyyyzz = pbuffer.data(idx_ovl_fi + 191);

    auto ts_yyy_yyyzzz = pbuffer.data(idx_ovl_fi + 192);

    auto ts_yyy_yyzzzz = pbuffer.data(idx_ovl_fi + 193);

    auto ts_yyy_yzzzzz = pbuffer.data(idx_ovl_fi + 194);

    auto ts_yyy_zzzzzz = pbuffer.data(idx_ovl_fi + 195);

    auto ts_yyz_yyyyyz = pbuffer.data(idx_ovl_fi + 218);

    auto ts_yyz_yyyyzz = pbuffer.data(idx_ovl_fi + 219);

    auto ts_yyz_yyyzzz = pbuffer.data(idx_ovl_fi + 220);

    auto ts_yyz_yyzzzz = pbuffer.data(idx_ovl_fi + 221);

    auto ts_yyz_yzzzzz = pbuffer.data(idx_ovl_fi + 222);

    auto ts_yyz_zzzzzz = pbuffer.data(idx_ovl_fi + 223);

    auto ts_yzz_xxxxxz = pbuffer.data(idx_ovl_fi + 226);

    auto ts_yzz_xxxxzz = pbuffer.data(idx_ovl_fi + 229);

    auto ts_yzz_xxxzzz = pbuffer.data(idx_ovl_fi + 233);

    auto ts_yzz_xxzzzz = pbuffer.data(idx_ovl_fi + 238);

    auto ts_yzz_xzzzzz = pbuffer.data(idx_ovl_fi + 244);

    auto ts_yzz_yyyyyy = pbuffer.data(idx_ovl_fi + 245);

    auto ts_yzz_yyyyyz = pbuffer.data(idx_ovl_fi + 246);

    auto ts_yzz_yyyyzz = pbuffer.data(idx_ovl_fi + 247);

    auto ts_yzz_yyyzzz = pbuffer.data(idx_ovl_fi + 248);

    auto ts_yzz_yyzzzz = pbuffer.data(idx_ovl_fi + 249);

    auto ts_yzz_yzzzzz = pbuffer.data(idx_ovl_fi + 250);

    auto ts_yzz_zzzzzz = pbuffer.data(idx_ovl_fi + 251);

    auto ts_zzz_xxxxxx = pbuffer.data(idx_ovl_fi + 252);

    auto ts_zzz_xxxxxy = pbuffer.data(idx_ovl_fi + 253);

    auto ts_zzz_xxxxxz = pbuffer.data(idx_ovl_fi + 254);

    auto ts_zzz_xxxxyy = pbuffer.data(idx_ovl_fi + 255);

    auto ts_zzz_xxxxyz = pbuffer.data(idx_ovl_fi + 256);

    auto ts_zzz_xxxxzz = pbuffer.data(idx_ovl_fi + 257);

    auto ts_zzz_xxxyyy = pbuffer.data(idx_ovl_fi + 258);

    auto ts_zzz_xxxyyz = pbuffer.data(idx_ovl_fi + 259);

    auto ts_zzz_xxxyzz = pbuffer.data(idx_ovl_fi + 260);

    auto ts_zzz_xxxzzz = pbuffer.data(idx_ovl_fi + 261);

    auto ts_zzz_xxyyyy = pbuffer.data(idx_ovl_fi + 262);

    auto ts_zzz_xxyyyz = pbuffer.data(idx_ovl_fi + 263);

    auto ts_zzz_xxyyzz = pbuffer.data(idx_ovl_fi + 264);

    auto ts_zzz_xxyzzz = pbuffer.data(idx_ovl_fi + 265);

    auto ts_zzz_xxzzzz = pbuffer.data(idx_ovl_fi + 266);

    auto ts_zzz_xyyyyy = pbuffer.data(idx_ovl_fi + 267);

    auto ts_zzz_xyyyyz = pbuffer.data(idx_ovl_fi + 268);

    auto ts_zzz_xyyyzz = pbuffer.data(idx_ovl_fi + 269);

    auto ts_zzz_xyyzzz = pbuffer.data(idx_ovl_fi + 270);

    auto ts_zzz_xyzzzz = pbuffer.data(idx_ovl_fi + 271);

    auto ts_zzz_xzzzzz = pbuffer.data(idx_ovl_fi + 272);

    auto ts_zzz_yyyyyy = pbuffer.data(idx_ovl_fi + 273);

    auto ts_zzz_yyyyyz = pbuffer.data(idx_ovl_fi + 274);

    auto ts_zzz_yyyyzz = pbuffer.data(idx_ovl_fi + 275);

    auto ts_zzz_yyyzzz = pbuffer.data(idx_ovl_fi + 276);

    auto ts_zzz_yyzzzz = pbuffer.data(idx_ovl_fi + 277);

    auto ts_zzz_yzzzzz = pbuffer.data(idx_ovl_fi + 278);

    auto ts_zzz_zzzzzz = pbuffer.data(idx_ovl_fi + 279);

    // Set up components of auxiliary buffer : FI

    auto tr_x_xxx_xxxxxx = pbuffer.data(idx_dip_fi);

    auto tr_x_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 1);

    auto tr_x_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 2);

    auto tr_x_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 3);

    auto tr_x_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 4);

    auto tr_x_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 5);

    auto tr_x_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 6);

    auto tr_x_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 7);

    auto tr_x_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 8);

    auto tr_x_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 9);

    auto tr_x_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 10);

    auto tr_x_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 11);

    auto tr_x_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 12);

    auto tr_x_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 13);

    auto tr_x_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 14);

    auto tr_x_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 15);

    auto tr_x_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 16);

    auto tr_x_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 17);

    auto tr_x_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 18);

    auto tr_x_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 19);

    auto tr_x_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 20);

    auto tr_x_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 21);

    auto tr_x_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 22);

    auto tr_x_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 23);

    auto tr_x_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 24);

    auto tr_x_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 25);

    auto tr_x_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 26);

    auto tr_x_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 27);

    auto tr_x_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 28);

    auto tr_x_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 29);

    auto tr_x_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 30);

    auto tr_x_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 31);

    auto tr_x_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 32);

    auto tr_x_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 33);

    auto tr_x_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 34);

    auto tr_x_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 35);

    auto tr_x_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 36);

    auto tr_x_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 37);

    auto tr_x_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 38);

    auto tr_x_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 39);

    auto tr_x_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 40);

    auto tr_x_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 41);

    auto tr_x_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 42);

    auto tr_x_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 43);

    auto tr_x_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 44);

    auto tr_x_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 45);

    auto tr_x_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 46);

    auto tr_x_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 47);

    auto tr_x_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 48);

    auto tr_x_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 49);

    auto tr_x_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 55);

    auto tr_x_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 56);

    auto tr_x_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 57);

    auto tr_x_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 58);

    auto tr_x_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 59);

    auto tr_x_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 60);

    auto tr_x_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 61);

    auto tr_x_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 62);

    auto tr_x_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 63);

    auto tr_x_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 64);

    auto tr_x_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 65);

    auto tr_x_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 66);

    auto tr_x_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 67);

    auto tr_x_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 68);

    auto tr_x_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 69);

    auto tr_x_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 70);

    auto tr_x_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 71);

    auto tr_x_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 72);

    auto tr_x_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 73);

    auto tr_x_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 74);

    auto tr_x_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 75);

    auto tr_x_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 76);

    auto tr_x_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 77);

    auto tr_x_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 78);

    auto tr_x_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 79);

    auto tr_x_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 80);

    auto tr_x_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 81);

    auto tr_x_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 82);

    auto tr_x_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 83);

    auto tr_x_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 84);

    auto tr_x_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 85);

    auto tr_x_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 86);

    auto tr_x_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 87);

    auto tr_x_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 88);

    auto tr_x_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 89);

    auto tr_x_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 90);

    auto tr_x_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 91);

    auto tr_x_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 92);

    auto tr_x_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 93);

    auto tr_x_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 94);

    auto tr_x_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 95);

    auto tr_x_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 96);

    auto tr_x_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 97);

    auto tr_x_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 98);

    auto tr_x_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 99);

    auto tr_x_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 100);

    auto tr_x_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 101);

    auto tr_x_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 102);

    auto tr_x_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 103);

    auto tr_x_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 104);

    auto tr_x_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 105);

    auto tr_x_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 106);

    auto tr_x_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 107);

    auto tr_x_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 108);

    auto tr_x_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 109);

    auto tr_x_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 110);

    auto tr_x_xyz_xxxxxz = pbuffer.data(idx_dip_fi + 114);

    auto tr_x_xyz_xxxxzz = pbuffer.data(idx_dip_fi + 117);

    auto tr_x_xyz_xxxzzz = pbuffer.data(idx_dip_fi + 121);

    auto tr_x_xyz_xxzzzz = pbuffer.data(idx_dip_fi + 126);

    auto tr_x_xyz_xzzzzz = pbuffer.data(idx_dip_fi + 132);

    auto tr_x_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 140);

    auto tr_x_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 141);

    auto tr_x_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 142);

    auto tr_x_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 143);

    auto tr_x_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 144);

    auto tr_x_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 145);

    auto tr_x_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 146);

    auto tr_x_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 147);

    auto tr_x_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 148);

    auto tr_x_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 149);

    auto tr_x_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 150);

    auto tr_x_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 151);

    auto tr_x_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 152);

    auto tr_x_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 153);

    auto tr_x_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 154);

    auto tr_x_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 155);

    auto tr_x_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 156);

    auto tr_x_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 157);

    auto tr_x_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 158);

    auto tr_x_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 159);

    auto tr_x_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 160);

    auto tr_x_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 162);

    auto tr_x_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 163);

    auto tr_x_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 164);

    auto tr_x_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 165);

    auto tr_x_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 166);

    auto tr_x_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 167);

    auto tr_x_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 168);

    auto tr_x_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 169);

    auto tr_x_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 170);

    auto tr_x_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 171);

    auto tr_x_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 172);

    auto tr_x_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 173);

    auto tr_x_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 174);

    auto tr_x_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 175);

    auto tr_x_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 176);

    auto tr_x_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 177);

    auto tr_x_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 178);

    auto tr_x_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 179);

    auto tr_x_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 180);

    auto tr_x_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 181);

    auto tr_x_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 182);

    auto tr_x_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 183);

    auto tr_x_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 184);

    auto tr_x_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 185);

    auto tr_x_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 186);

    auto tr_x_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 187);

    auto tr_x_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 188);

    auto tr_x_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 189);

    auto tr_x_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 190);

    auto tr_x_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 191);

    auto tr_x_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 192);

    auto tr_x_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 193);

    auto tr_x_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 194);

    auto tr_x_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 195);

    auto tr_x_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 197);

    auto tr_x_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 198);

    auto tr_x_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 199);

    auto tr_x_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 201);

    auto tr_x_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 202);

    auto tr_x_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 205);

    auto tr_x_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 206);

    auto tr_x_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 210);

    auto tr_x_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 211);

    auto tr_x_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 216);

    auto tr_x_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 217);

    auto tr_x_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 218);

    auto tr_x_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 219);

    auto tr_x_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 220);

    auto tr_x_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 221);

    auto tr_x_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 222);

    auto tr_x_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 223);

    auto tr_x_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 224);

    auto tr_x_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 226);

    auto tr_x_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 228);

    auto tr_x_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 229);

    auto tr_x_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 231);

    auto tr_x_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 232);

    auto tr_x_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 233);

    auto tr_x_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 235);

    auto tr_x_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 236);

    auto tr_x_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 237);

    auto tr_x_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 238);

    auto tr_x_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 240);

    auto tr_x_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 241);

    auto tr_x_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 242);

    auto tr_x_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 243);

    auto tr_x_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 244);

    auto tr_x_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 245);

    auto tr_x_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 246);

    auto tr_x_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 247);

    auto tr_x_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 248);

    auto tr_x_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 249);

    auto tr_x_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 250);

    auto tr_x_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 251);

    auto tr_x_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 252);

    auto tr_x_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 253);

    auto tr_x_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 254);

    auto tr_x_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 255);

    auto tr_x_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 256);

    auto tr_x_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 257);

    auto tr_x_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 258);

    auto tr_x_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 259);

    auto tr_x_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 260);

    auto tr_x_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 261);

    auto tr_x_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 262);

    auto tr_x_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 263);

    auto tr_x_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 264);

    auto tr_x_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 265);

    auto tr_x_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 266);

    auto tr_x_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 267);

    auto tr_x_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 268);

    auto tr_x_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 269);

    auto tr_x_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 270);

    auto tr_x_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 271);

    auto tr_x_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 272);

    auto tr_x_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 273);

    auto tr_x_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 274);

    auto tr_x_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 275);

    auto tr_x_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 276);

    auto tr_x_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 277);

    auto tr_x_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 278);

    auto tr_x_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 279);

    auto tr_y_xxx_xxxxxx = pbuffer.data(idx_dip_fi + 280);

    auto tr_y_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 281);

    auto tr_y_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 282);

    auto tr_y_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 283);

    auto tr_y_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 284);

    auto tr_y_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 285);

    auto tr_y_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 286);

    auto tr_y_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 287);

    auto tr_y_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 288);

    auto tr_y_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 289);

    auto tr_y_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 290);

    auto tr_y_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 291);

    auto tr_y_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 292);

    auto tr_y_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 293);

    auto tr_y_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 294);

    auto tr_y_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 295);

    auto tr_y_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 296);

    auto tr_y_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 297);

    auto tr_y_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 298);

    auto tr_y_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 299);

    auto tr_y_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 300);

    auto tr_y_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 301);

    auto tr_y_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 302);

    auto tr_y_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 303);

    auto tr_y_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 304);

    auto tr_y_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 305);

    auto tr_y_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 306);

    auto tr_y_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 307);

    auto tr_y_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 308);

    auto tr_y_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 309);

    auto tr_y_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 311);

    auto tr_y_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 312);

    auto tr_y_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 314);

    auto tr_y_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 315);

    auto tr_y_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 316);

    auto tr_y_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 318);

    auto tr_y_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 319);

    auto tr_y_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 320);

    auto tr_y_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 321);

    auto tr_y_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 323);

    auto tr_y_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 324);

    auto tr_y_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 325);

    auto tr_y_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 326);

    auto tr_y_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 327);

    auto tr_y_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 329);

    auto tr_y_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 330);

    auto tr_y_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 331);

    auto tr_y_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 332);

    auto tr_y_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 333);

    auto tr_y_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 334);

    auto tr_y_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 335);

    auto tr_y_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 336);

    auto tr_y_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 337);

    auto tr_y_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 338);

    auto tr_y_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 339);

    auto tr_y_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 341);

    auto tr_y_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 342);

    auto tr_y_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 345);

    auto tr_y_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 346);

    auto tr_y_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 350);

    auto tr_y_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 351);

    auto tr_y_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 356);

    auto tr_y_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 358);

    auto tr_y_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 359);

    auto tr_y_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 360);

    auto tr_y_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 361);

    auto tr_y_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 362);

    auto tr_y_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 363);

    auto tr_y_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 364);

    auto tr_y_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 365);

    auto tr_y_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 366);

    auto tr_y_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 367);

    auto tr_y_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 368);

    auto tr_y_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 369);

    auto tr_y_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 370);

    auto tr_y_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 371);

    auto tr_y_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 372);

    auto tr_y_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 373);

    auto tr_y_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 374);

    auto tr_y_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 375);

    auto tr_y_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 376);

    auto tr_y_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 377);

    auto tr_y_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 378);

    auto tr_y_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 379);

    auto tr_y_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 380);

    auto tr_y_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 381);

    auto tr_y_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 382);

    auto tr_y_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 383);

    auto tr_y_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 384);

    auto tr_y_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 385);

    auto tr_y_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 386);

    auto tr_y_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 387);

    auto tr_y_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 388);

    auto tr_y_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 389);

    auto tr_y_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 390);

    auto tr_y_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 391);

    auto tr_y_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 414);

    auto tr_y_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 415);

    auto tr_y_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 416);

    auto tr_y_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 417);

    auto tr_y_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 418);

    auto tr_y_xyz_zzzzzz = pbuffer.data(idx_dip_fi + 419);

    auto tr_y_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 422);

    auto tr_y_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 424);

    auto tr_y_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 425);

    auto tr_y_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 427);

    auto tr_y_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 428);

    auto tr_y_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 429);

    auto tr_y_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 431);

    auto tr_y_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 432);

    auto tr_y_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 433);

    auto tr_y_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 434);

    auto tr_y_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 436);

    auto tr_y_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 437);

    auto tr_y_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 438);

    auto tr_y_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 439);

    auto tr_y_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 440);

    auto tr_y_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 441);

    auto tr_y_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 442);

    auto tr_y_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 443);

    auto tr_y_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 444);

    auto tr_y_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 445);

    auto tr_y_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 446);

    auto tr_y_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 447);

    auto tr_y_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 448);

    auto tr_y_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 449);

    auto tr_y_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 450);

    auto tr_y_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 451);

    auto tr_y_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 452);

    auto tr_y_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 453);

    auto tr_y_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 454);

    auto tr_y_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 455);

    auto tr_y_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 456);

    auto tr_y_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 457);

    auto tr_y_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 458);

    auto tr_y_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 459);

    auto tr_y_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 460);

    auto tr_y_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 461);

    auto tr_y_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 462);

    auto tr_y_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 463);

    auto tr_y_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 464);

    auto tr_y_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 465);

    auto tr_y_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 466);

    auto tr_y_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 467);

    auto tr_y_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 468);

    auto tr_y_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 469);

    auto tr_y_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 470);

    auto tr_y_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 471);

    auto tr_y_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 472);

    auto tr_y_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 473);

    auto tr_y_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 474);

    auto tr_y_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 475);

    auto tr_y_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 476);

    auto tr_y_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 477);

    auto tr_y_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 478);

    auto tr_y_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 479);

    auto tr_y_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 480);

    auto tr_y_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 481);

    auto tr_y_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 482);

    auto tr_y_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 483);

    auto tr_y_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 484);

    auto tr_y_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 485);

    auto tr_y_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 486);

    auto tr_y_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 487);

    auto tr_y_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 488);

    auto tr_y_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 489);

    auto tr_y_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 490);

    auto tr_y_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 491);

    auto tr_y_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 492);

    auto tr_y_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 493);

    auto tr_y_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 494);

    auto tr_y_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 495);

    auto tr_y_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 496);

    auto tr_y_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 497);

    auto tr_y_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 498);

    auto tr_y_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 499);

    auto tr_y_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 500);

    auto tr_y_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 501);

    auto tr_y_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 502);

    auto tr_y_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 503);

    auto tr_y_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 504);

    auto tr_y_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 505);

    auto tr_y_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 506);

    auto tr_y_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 507);

    auto tr_y_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 508);

    auto tr_y_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 509);

    auto tr_y_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 510);

    auto tr_y_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 511);

    auto tr_y_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 512);

    auto tr_y_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 513);

    auto tr_y_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 514);

    auto tr_y_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 515);

    auto tr_y_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 516);

    auto tr_y_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 517);

    auto tr_y_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 518);

    auto tr_y_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 519);

    auto tr_y_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 520);

    auto tr_y_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 521);

    auto tr_y_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 522);

    auto tr_y_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 523);

    auto tr_y_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 524);

    auto tr_y_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 525);

    auto tr_y_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 526);

    auto tr_y_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 527);

    auto tr_y_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 528);

    auto tr_y_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 529);

    auto tr_y_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 530);

    auto tr_y_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 531);

    auto tr_y_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 532);

    auto tr_y_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 533);

    auto tr_y_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 534);

    auto tr_y_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 535);

    auto tr_y_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 536);

    auto tr_y_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 537);

    auto tr_y_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 538);

    auto tr_y_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 539);

    auto tr_y_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 540);

    auto tr_y_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 541);

    auto tr_y_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 542);

    auto tr_y_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 543);

    auto tr_y_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 544);

    auto tr_y_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 545);

    auto tr_y_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 546);

    auto tr_y_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 547);

    auto tr_y_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 548);

    auto tr_y_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 549);

    auto tr_y_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 550);

    auto tr_y_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 551);

    auto tr_y_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 552);

    auto tr_y_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 553);

    auto tr_y_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 554);

    auto tr_y_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 555);

    auto tr_y_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 556);

    auto tr_y_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 557);

    auto tr_y_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 558);

    auto tr_y_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 559);

    auto tr_z_xxx_xxxxxx = pbuffer.data(idx_dip_fi + 560);

    auto tr_z_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 561);

    auto tr_z_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 562);

    auto tr_z_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 563);

    auto tr_z_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 564);

    auto tr_z_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 565);

    auto tr_z_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 566);

    auto tr_z_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 567);

    auto tr_z_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 568);

    auto tr_z_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 569);

    auto tr_z_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 570);

    auto tr_z_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 571);

    auto tr_z_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 572);

    auto tr_z_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 573);

    auto tr_z_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 574);

    auto tr_z_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 575);

    auto tr_z_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 576);

    auto tr_z_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 577);

    auto tr_z_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 578);

    auto tr_z_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 579);

    auto tr_z_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 580);

    auto tr_z_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 581);

    auto tr_z_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 582);

    auto tr_z_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 583);

    auto tr_z_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 584);

    auto tr_z_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 585);

    auto tr_z_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 586);

    auto tr_z_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 587);

    auto tr_z_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 588);

    auto tr_z_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 590);

    auto tr_z_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 593);

    auto tr_z_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 597);

    auto tr_z_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 602);

    auto tr_z_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 608);

    auto tr_z_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 609);

    auto tr_z_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 610);

    auto tr_z_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 611);

    auto tr_z_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 612);

    auto tr_z_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 613);

    auto tr_z_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 614);

    auto tr_z_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 616);

    auto tr_z_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 617);

    auto tr_z_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 618);

    auto tr_z_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 619);

    auto tr_z_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 620);

    auto tr_z_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 621);

    auto tr_z_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 622);

    auto tr_z_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 623);

    auto tr_z_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 624);

    auto tr_z_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 625);

    auto tr_z_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 626);

    auto tr_z_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 627);

    auto tr_z_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 628);

    auto tr_z_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 629);

    auto tr_z_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 630);

    auto tr_z_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 631);

    auto tr_z_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 632);

    auto tr_z_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 633);

    auto tr_z_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 634);

    auto tr_z_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 635);

    auto tr_z_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 636);

    auto tr_z_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 637);

    auto tr_z_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 638);

    auto tr_z_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 639);

    auto tr_z_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 640);

    auto tr_z_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 641);

    auto tr_z_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 642);

    auto tr_z_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 643);

    auto tr_z_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 645);

    auto tr_z_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 647);

    auto tr_z_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 648);

    auto tr_z_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 650);

    auto tr_z_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 651);

    auto tr_z_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 652);

    auto tr_z_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 654);

    auto tr_z_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 655);

    auto tr_z_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 656);

    auto tr_z_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 657);

    auto tr_z_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 659);

    auto tr_z_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 660);

    auto tr_z_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 661);

    auto tr_z_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 662);

    auto tr_z_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 663);

    auto tr_z_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 665);

    auto tr_z_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 666);

    auto tr_z_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 667);

    auto tr_z_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 668);

    auto tr_z_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 669);

    auto tr_z_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 670);

    auto tr_z_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 671);

    auto tr_z_xyz_yyyyyy = pbuffer.data(idx_dip_fi + 693);

    auto tr_z_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 694);

    auto tr_z_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 695);

    auto tr_z_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 696);

    auto tr_z_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 697);

    auto tr_z_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 698);

    auto tr_z_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 700);

    auto tr_z_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 701);

    auto tr_z_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 702);

    auto tr_z_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 703);

    auto tr_z_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 704);

    auto tr_z_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 705);

    auto tr_z_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 706);

    auto tr_z_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 707);

    auto tr_z_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 708);

    auto tr_z_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 709);

    auto tr_z_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 710);

    auto tr_z_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 711);

    auto tr_z_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 712);

    auto tr_z_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 713);

    auto tr_z_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 714);

    auto tr_z_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 715);

    auto tr_z_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 716);

    auto tr_z_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 717);

    auto tr_z_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 718);

    auto tr_z_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 719);

    auto tr_z_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 720);

    auto tr_z_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 721);

    auto tr_z_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 722);

    auto tr_z_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 723);

    auto tr_z_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 724);

    auto tr_z_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 725);

    auto tr_z_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 726);

    auto tr_z_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 727);

    auto tr_z_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 728);

    auto tr_z_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 729);

    auto tr_z_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 730);

    auto tr_z_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 731);

    auto tr_z_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 732);

    auto tr_z_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 733);

    auto tr_z_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 734);

    auto tr_z_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 735);

    auto tr_z_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 736);

    auto tr_z_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 737);

    auto tr_z_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 738);

    auto tr_z_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 739);

    auto tr_z_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 740);

    auto tr_z_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 741);

    auto tr_z_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 742);

    auto tr_z_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 743);

    auto tr_z_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 744);

    auto tr_z_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 745);

    auto tr_z_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 746);

    auto tr_z_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 747);

    auto tr_z_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 748);

    auto tr_z_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 749);

    auto tr_z_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 750);

    auto tr_z_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 751);

    auto tr_z_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 752);

    auto tr_z_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 753);

    auto tr_z_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 754);

    auto tr_z_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 755);

    auto tr_z_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 756);

    auto tr_z_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 757);

    auto tr_z_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 758);

    auto tr_z_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 759);

    auto tr_z_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 760);

    auto tr_z_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 761);

    auto tr_z_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 762);

    auto tr_z_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 763);

    auto tr_z_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 764);

    auto tr_z_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 765);

    auto tr_z_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 766);

    auto tr_z_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 767);

    auto tr_z_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 768);

    auto tr_z_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 769);

    auto tr_z_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 770);

    auto tr_z_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 771);

    auto tr_z_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 772);

    auto tr_z_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 773);

    auto tr_z_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 774);

    auto tr_z_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 775);

    auto tr_z_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 776);

    auto tr_z_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 777);

    auto tr_z_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 778);

    auto tr_z_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 779);

    auto tr_z_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 780);

    auto tr_z_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 781);

    auto tr_z_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 782);

    auto tr_z_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 783);

    auto tr_z_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 784);

    auto tr_z_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 785);

    auto tr_z_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 786);

    auto tr_z_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 787);

    auto tr_z_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 788);

    auto tr_z_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 789);

    auto tr_z_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 790);

    auto tr_z_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 791);

    auto tr_z_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 792);

    auto tr_z_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 793);

    auto tr_z_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 794);

    auto tr_z_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 795);

    auto tr_z_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 796);

    auto tr_z_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 797);

    auto tr_z_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 798);

    auto tr_z_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 799);

    auto tr_z_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 800);

    auto tr_z_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 801);

    auto tr_z_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 802);

    auto tr_z_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 803);

    auto tr_z_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 804);

    auto tr_z_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 805);

    auto tr_z_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 806);

    auto tr_z_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 807);

    auto tr_z_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 808);

    auto tr_z_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 809);

    auto tr_z_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 810);

    auto tr_z_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 811);

    auto tr_z_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 812);

    auto tr_z_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 813);

    auto tr_z_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 814);

    auto tr_z_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 815);

    auto tr_z_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 816);

    auto tr_z_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 817);

    auto tr_z_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 818);

    auto tr_z_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 819);

    auto tr_z_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 820);

    auto tr_z_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 821);

    auto tr_z_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 822);

    auto tr_z_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 823);

    auto tr_z_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 824);

    auto tr_z_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 825);

    auto tr_z_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 826);

    auto tr_z_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 827);

    auto tr_z_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 828);

    auto tr_z_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 829);

    auto tr_z_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 830);

    auto tr_z_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 831);

    auto tr_z_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 832);

    auto tr_z_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 833);

    auto tr_z_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 834);

    auto tr_z_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 835);

    auto tr_z_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 836);

    auto tr_z_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 837);

    auto tr_z_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 838);

    auto tr_z_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 839);

    // Set up 0-28 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_x_xx_xxxxxx, tr_x_xx_xxxxxy, tr_x_xx_xxxxxz, tr_x_xx_xxxxyy, tr_x_xx_xxxxyz, tr_x_xx_xxxxzz, tr_x_xx_xxxyyy, tr_x_xx_xxxyyz, tr_x_xx_xxxyzz, tr_x_xx_xxxzzz, tr_x_xx_xxyyyy, tr_x_xx_xxyyyz, tr_x_xx_xxyyzz, tr_x_xx_xxyzzz, tr_x_xx_xxzzzz, tr_x_xx_xyyyyy, tr_x_xx_xyyyyz, tr_x_xx_xyyyzz, tr_x_xx_xyyzzz, tr_x_xx_xyzzzz, tr_x_xx_xzzzzz, tr_x_xx_yyyyyy, tr_x_xx_yyyyyz, tr_x_xx_yyyyzz, tr_x_xx_yyyzzz, tr_x_xx_yyzzzz, tr_x_xx_yzzzzz, tr_x_xx_zzzzzz, tr_x_xxx_xxxxx, tr_x_xxx_xxxxxx, tr_x_xxx_xxxxxy, tr_x_xxx_xxxxxz, tr_x_xxx_xxxxy, tr_x_xxx_xxxxyy, tr_x_xxx_xxxxyz, tr_x_xxx_xxxxz, tr_x_xxx_xxxxzz, tr_x_xxx_xxxyy, tr_x_xxx_xxxyyy, tr_x_xxx_xxxyyz, tr_x_xxx_xxxyz, tr_x_xxx_xxxyzz, tr_x_xxx_xxxzz, tr_x_xxx_xxxzzz, tr_x_xxx_xxyyy, tr_x_xxx_xxyyyy, tr_x_xxx_xxyyyz, tr_x_xxx_xxyyz, tr_x_xxx_xxyyzz, tr_x_xxx_xxyzz, tr_x_xxx_xxyzzz, tr_x_xxx_xxzzz, tr_x_xxx_xxzzzz, tr_x_xxx_xyyyy, tr_x_xxx_xyyyyy, tr_x_xxx_xyyyyz, tr_x_xxx_xyyyz, tr_x_xxx_xyyyzz, tr_x_xxx_xyyzz, tr_x_xxx_xyyzzz, tr_x_xxx_xyzzz, tr_x_xxx_xyzzzz, tr_x_xxx_xzzzz, tr_x_xxx_xzzzzz, tr_x_xxx_yyyyy, tr_x_xxx_yyyyyy, tr_x_xxx_yyyyyz, tr_x_xxx_yyyyz, tr_x_xxx_yyyyzz, tr_x_xxx_yyyzz, tr_x_xxx_yyyzzz, tr_x_xxx_yyzzz, tr_x_xxx_yyzzzz, tr_x_xxx_yzzzz, tr_x_xxx_yzzzzz, tr_x_xxx_zzzzz, tr_x_xxx_zzzzzz, tr_x_xxxx_xxxxxx, tr_x_xxxx_xxxxxy, tr_x_xxxx_xxxxxz, tr_x_xxxx_xxxxyy, tr_x_xxxx_xxxxyz, tr_x_xxxx_xxxxzz, tr_x_xxxx_xxxyyy, tr_x_xxxx_xxxyyz, tr_x_xxxx_xxxyzz, tr_x_xxxx_xxxzzz, tr_x_xxxx_xxyyyy, tr_x_xxxx_xxyyyz, tr_x_xxxx_xxyyzz, tr_x_xxxx_xxyzzz, tr_x_xxxx_xxzzzz, tr_x_xxxx_xyyyyy, tr_x_xxxx_xyyyyz, tr_x_xxxx_xyyyzz, tr_x_xxxx_xyyzzz, tr_x_xxxx_xyzzzz, tr_x_xxxx_xzzzzz, tr_x_xxxx_yyyyyy, tr_x_xxxx_yyyyyz, tr_x_xxxx_yyyyzz, tr_x_xxxx_yyyzzz, tr_x_xxxx_yyzzzz, tr_x_xxxx_yzzzzz, tr_x_xxxx_zzzzzz, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxxz, ts_xxx_xxxxyy, ts_xxx_xxxxyz, ts_xxx_xxxxzz, ts_xxx_xxxyyy, ts_xxx_xxxyyz, ts_xxx_xxxyzz, ts_xxx_xxxzzz, ts_xxx_xxyyyy, ts_xxx_xxyyyz, ts_xxx_xxyyzz, ts_xxx_xxyzzz, ts_xxx_xxzzzz, ts_xxx_xyyyyy, ts_xxx_xyyyyz, ts_xxx_xyyyzz, ts_xxx_xyyzzz, ts_xxx_xyzzzz, ts_xxx_xzzzzz, ts_xxx_yyyyyy, ts_xxx_yyyyyz, ts_xxx_yyyyzz, ts_xxx_yyyzzz, ts_xxx_yyzzzz, ts_xxx_yzzzzz, ts_xxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_xxxxxx[i] = 3.0 * tr_x_xx_xxxxxx[i] * fe_0 + 6.0 * tr_x_xxx_xxxxx[i] * fe_0 + ts_xxx_xxxxxx[i] * fe_0 + tr_x_xxx_xxxxxx[i] * pa_x[i];

        tr_x_xxxx_xxxxxy[i] = 3.0 * tr_x_xx_xxxxxy[i] * fe_0 + 5.0 * tr_x_xxx_xxxxy[i] * fe_0 + ts_xxx_xxxxxy[i] * fe_0 + tr_x_xxx_xxxxxy[i] * pa_x[i];

        tr_x_xxxx_xxxxxz[i] = 3.0 * tr_x_xx_xxxxxz[i] * fe_0 + 5.0 * tr_x_xxx_xxxxz[i] * fe_0 + ts_xxx_xxxxxz[i] * fe_0 + tr_x_xxx_xxxxxz[i] * pa_x[i];

        tr_x_xxxx_xxxxyy[i] = 3.0 * tr_x_xx_xxxxyy[i] * fe_0 + 4.0 * tr_x_xxx_xxxyy[i] * fe_0 + ts_xxx_xxxxyy[i] * fe_0 + tr_x_xxx_xxxxyy[i] * pa_x[i];

        tr_x_xxxx_xxxxyz[i] = 3.0 * tr_x_xx_xxxxyz[i] * fe_0 + 4.0 * tr_x_xxx_xxxyz[i] * fe_0 + ts_xxx_xxxxyz[i] * fe_0 + tr_x_xxx_xxxxyz[i] * pa_x[i];

        tr_x_xxxx_xxxxzz[i] = 3.0 * tr_x_xx_xxxxzz[i] * fe_0 + 4.0 * tr_x_xxx_xxxzz[i] * fe_0 + ts_xxx_xxxxzz[i] * fe_0 + tr_x_xxx_xxxxzz[i] * pa_x[i];

        tr_x_xxxx_xxxyyy[i] = 3.0 * tr_x_xx_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxx_xxyyy[i] * fe_0 + ts_xxx_xxxyyy[i] * fe_0 + tr_x_xxx_xxxyyy[i] * pa_x[i];

        tr_x_xxxx_xxxyyz[i] = 3.0 * tr_x_xx_xxxyyz[i] * fe_0 + 3.0 * tr_x_xxx_xxyyz[i] * fe_0 + ts_xxx_xxxyyz[i] * fe_0 + tr_x_xxx_xxxyyz[i] * pa_x[i];

        tr_x_xxxx_xxxyzz[i] = 3.0 * tr_x_xx_xxxyzz[i] * fe_0 + 3.0 * tr_x_xxx_xxyzz[i] * fe_0 + ts_xxx_xxxyzz[i] * fe_0 + tr_x_xxx_xxxyzz[i] * pa_x[i];

        tr_x_xxxx_xxxzzz[i] = 3.0 * tr_x_xx_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxx_xxzzz[i] * fe_0 + ts_xxx_xxxzzz[i] * fe_0 + tr_x_xxx_xxxzzz[i] * pa_x[i];

        tr_x_xxxx_xxyyyy[i] = 3.0 * tr_x_xx_xxyyyy[i] * fe_0 + 2.0 * tr_x_xxx_xyyyy[i] * fe_0 + ts_xxx_xxyyyy[i] * fe_0 + tr_x_xxx_xxyyyy[i] * pa_x[i];

        tr_x_xxxx_xxyyyz[i] = 3.0 * tr_x_xx_xxyyyz[i] * fe_0 + 2.0 * tr_x_xxx_xyyyz[i] * fe_0 + ts_xxx_xxyyyz[i] * fe_0 + tr_x_xxx_xxyyyz[i] * pa_x[i];

        tr_x_xxxx_xxyyzz[i] = 3.0 * tr_x_xx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxx_xyyzz[i] * fe_0 + ts_xxx_xxyyzz[i] * fe_0 + tr_x_xxx_xxyyzz[i] * pa_x[i];

        tr_x_xxxx_xxyzzz[i] = 3.0 * tr_x_xx_xxyzzz[i] * fe_0 + 2.0 * tr_x_xxx_xyzzz[i] * fe_0 + ts_xxx_xxyzzz[i] * fe_0 + tr_x_xxx_xxyzzz[i] * pa_x[i];

        tr_x_xxxx_xxzzzz[i] = 3.0 * tr_x_xx_xxzzzz[i] * fe_0 + 2.0 * tr_x_xxx_xzzzz[i] * fe_0 + ts_xxx_xxzzzz[i] * fe_0 + tr_x_xxx_xxzzzz[i] * pa_x[i];

        tr_x_xxxx_xyyyyy[i] = 3.0 * tr_x_xx_xyyyyy[i] * fe_0 + tr_x_xxx_yyyyy[i] * fe_0 + ts_xxx_xyyyyy[i] * fe_0 + tr_x_xxx_xyyyyy[i] * pa_x[i];

        tr_x_xxxx_xyyyyz[i] = 3.0 * tr_x_xx_xyyyyz[i] * fe_0 + tr_x_xxx_yyyyz[i] * fe_0 + ts_xxx_xyyyyz[i] * fe_0 + tr_x_xxx_xyyyyz[i] * pa_x[i];

        tr_x_xxxx_xyyyzz[i] = 3.0 * tr_x_xx_xyyyzz[i] * fe_0 + tr_x_xxx_yyyzz[i] * fe_0 + ts_xxx_xyyyzz[i] * fe_0 + tr_x_xxx_xyyyzz[i] * pa_x[i];

        tr_x_xxxx_xyyzzz[i] = 3.0 * tr_x_xx_xyyzzz[i] * fe_0 + tr_x_xxx_yyzzz[i] * fe_0 + ts_xxx_xyyzzz[i] * fe_0 + tr_x_xxx_xyyzzz[i] * pa_x[i];

        tr_x_xxxx_xyzzzz[i] = 3.0 * tr_x_xx_xyzzzz[i] * fe_0 + tr_x_xxx_yzzzz[i] * fe_0 + ts_xxx_xyzzzz[i] * fe_0 + tr_x_xxx_xyzzzz[i] * pa_x[i];

        tr_x_xxxx_xzzzzz[i] = 3.0 * tr_x_xx_xzzzzz[i] * fe_0 + tr_x_xxx_zzzzz[i] * fe_0 + ts_xxx_xzzzzz[i] * fe_0 + tr_x_xxx_xzzzzz[i] * pa_x[i];

        tr_x_xxxx_yyyyyy[i] = 3.0 * tr_x_xx_yyyyyy[i] * fe_0 + ts_xxx_yyyyyy[i] * fe_0 + tr_x_xxx_yyyyyy[i] * pa_x[i];

        tr_x_xxxx_yyyyyz[i] = 3.0 * tr_x_xx_yyyyyz[i] * fe_0 + ts_xxx_yyyyyz[i] * fe_0 + tr_x_xxx_yyyyyz[i] * pa_x[i];

        tr_x_xxxx_yyyyzz[i] = 3.0 * tr_x_xx_yyyyzz[i] * fe_0 + ts_xxx_yyyyzz[i] * fe_0 + tr_x_xxx_yyyyzz[i] * pa_x[i];

        tr_x_xxxx_yyyzzz[i] = 3.0 * tr_x_xx_yyyzzz[i] * fe_0 + ts_xxx_yyyzzz[i] * fe_0 + tr_x_xxx_yyyzzz[i] * pa_x[i];

        tr_x_xxxx_yyzzzz[i] = 3.0 * tr_x_xx_yyzzzz[i] * fe_0 + ts_xxx_yyzzzz[i] * fe_0 + tr_x_xxx_yyzzzz[i] * pa_x[i];

        tr_x_xxxx_yzzzzz[i] = 3.0 * tr_x_xx_yzzzzz[i] * fe_0 + ts_xxx_yzzzzz[i] * fe_0 + tr_x_xxx_yzzzzz[i] * pa_x[i];

        tr_x_xxxx_zzzzzz[i] = 3.0 * tr_x_xx_zzzzzz[i] * fe_0 + ts_xxx_zzzzzz[i] * fe_0 + tr_x_xxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : GI

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

    auto tr_x_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 49);

    auto tr_x_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 50);

    auto tr_x_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 51);

    auto tr_x_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 52);

    auto tr_x_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 53);

    auto tr_x_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 54);

    auto tr_x_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 55);

    #pragma omp simd aligned(pa_y, tr_x_xxx_xxxxx, tr_x_xxx_xxxxxx, tr_x_xxx_xxxxxy, tr_x_xxx_xxxxxz, tr_x_xxx_xxxxy, tr_x_xxx_xxxxyy, tr_x_xxx_xxxxyz, tr_x_xxx_xxxxz, tr_x_xxx_xxxxzz, tr_x_xxx_xxxyy, tr_x_xxx_xxxyyy, tr_x_xxx_xxxyyz, tr_x_xxx_xxxyz, tr_x_xxx_xxxyzz, tr_x_xxx_xxxzz, tr_x_xxx_xxxzzz, tr_x_xxx_xxyyy, tr_x_xxx_xxyyyy, tr_x_xxx_xxyyyz, tr_x_xxx_xxyyz, tr_x_xxx_xxyyzz, tr_x_xxx_xxyzz, tr_x_xxx_xxyzzz, tr_x_xxx_xxzzz, tr_x_xxx_xxzzzz, tr_x_xxx_xyyyy, tr_x_xxx_xyyyyy, tr_x_xxx_xyyyyz, tr_x_xxx_xyyyz, tr_x_xxx_xyyyzz, tr_x_xxx_xyyzz, tr_x_xxx_xyyzzz, tr_x_xxx_xyzzz, tr_x_xxx_xyzzzz, tr_x_xxx_xzzzz, tr_x_xxx_xzzzzz, tr_x_xxx_yyyyy, tr_x_xxx_yyyyyy, tr_x_xxx_yyyyyz, tr_x_xxx_yyyyz, tr_x_xxx_yyyyzz, tr_x_xxx_yyyzz, tr_x_xxx_yyyzzz, tr_x_xxx_yyzzz, tr_x_xxx_yyzzzz, tr_x_xxx_yzzzz, tr_x_xxx_yzzzzz, tr_x_xxx_zzzzz, tr_x_xxx_zzzzzz, tr_x_xxxy_xxxxxx, tr_x_xxxy_xxxxxy, tr_x_xxxy_xxxxxz, tr_x_xxxy_xxxxyy, tr_x_xxxy_xxxxyz, tr_x_xxxy_xxxxzz, tr_x_xxxy_xxxyyy, tr_x_xxxy_xxxyyz, tr_x_xxxy_xxxyzz, tr_x_xxxy_xxxzzz, tr_x_xxxy_xxyyyy, tr_x_xxxy_xxyyyz, tr_x_xxxy_xxyyzz, tr_x_xxxy_xxyzzz, tr_x_xxxy_xxzzzz, tr_x_xxxy_xyyyyy, tr_x_xxxy_xyyyyz, tr_x_xxxy_xyyyzz, tr_x_xxxy_xyyzzz, tr_x_xxxy_xyzzzz, tr_x_xxxy_xzzzzz, tr_x_xxxy_yyyyyy, tr_x_xxxy_yyyyyz, tr_x_xxxy_yyyyzz, tr_x_xxxy_yyyzzz, tr_x_xxxy_yyzzzz, tr_x_xxxy_yzzzzz, tr_x_xxxy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxy_xxxxxx[i] = tr_x_xxx_xxxxxx[i] * pa_y[i];

        tr_x_xxxy_xxxxxy[i] = tr_x_xxx_xxxxx[i] * fe_0 + tr_x_xxx_xxxxxy[i] * pa_y[i];

        tr_x_xxxy_xxxxxz[i] = tr_x_xxx_xxxxxz[i] * pa_y[i];

        tr_x_xxxy_xxxxyy[i] = 2.0 * tr_x_xxx_xxxxy[i] * fe_0 + tr_x_xxx_xxxxyy[i] * pa_y[i];

        tr_x_xxxy_xxxxyz[i] = tr_x_xxx_xxxxz[i] * fe_0 + tr_x_xxx_xxxxyz[i] * pa_y[i];

        tr_x_xxxy_xxxxzz[i] = tr_x_xxx_xxxxzz[i] * pa_y[i];

        tr_x_xxxy_xxxyyy[i] = 3.0 * tr_x_xxx_xxxyy[i] * fe_0 + tr_x_xxx_xxxyyy[i] * pa_y[i];

        tr_x_xxxy_xxxyyz[i] = 2.0 * tr_x_xxx_xxxyz[i] * fe_0 + tr_x_xxx_xxxyyz[i] * pa_y[i];

        tr_x_xxxy_xxxyzz[i] = tr_x_xxx_xxxzz[i] * fe_0 + tr_x_xxx_xxxyzz[i] * pa_y[i];

        tr_x_xxxy_xxxzzz[i] = tr_x_xxx_xxxzzz[i] * pa_y[i];

        tr_x_xxxy_xxyyyy[i] = 4.0 * tr_x_xxx_xxyyy[i] * fe_0 + tr_x_xxx_xxyyyy[i] * pa_y[i];

        tr_x_xxxy_xxyyyz[i] = 3.0 * tr_x_xxx_xxyyz[i] * fe_0 + tr_x_xxx_xxyyyz[i] * pa_y[i];

        tr_x_xxxy_xxyyzz[i] = 2.0 * tr_x_xxx_xxyzz[i] * fe_0 + tr_x_xxx_xxyyzz[i] * pa_y[i];

        tr_x_xxxy_xxyzzz[i] = tr_x_xxx_xxzzz[i] * fe_0 + tr_x_xxx_xxyzzz[i] * pa_y[i];

        tr_x_xxxy_xxzzzz[i] = tr_x_xxx_xxzzzz[i] * pa_y[i];

        tr_x_xxxy_xyyyyy[i] = 5.0 * tr_x_xxx_xyyyy[i] * fe_0 + tr_x_xxx_xyyyyy[i] * pa_y[i];

        tr_x_xxxy_xyyyyz[i] = 4.0 * tr_x_xxx_xyyyz[i] * fe_0 + tr_x_xxx_xyyyyz[i] * pa_y[i];

        tr_x_xxxy_xyyyzz[i] = 3.0 * tr_x_xxx_xyyzz[i] * fe_0 + tr_x_xxx_xyyyzz[i] * pa_y[i];

        tr_x_xxxy_xyyzzz[i] = 2.0 * tr_x_xxx_xyzzz[i] * fe_0 + tr_x_xxx_xyyzzz[i] * pa_y[i];

        tr_x_xxxy_xyzzzz[i] = tr_x_xxx_xzzzz[i] * fe_0 + tr_x_xxx_xyzzzz[i] * pa_y[i];

        tr_x_xxxy_xzzzzz[i] = tr_x_xxx_xzzzzz[i] * pa_y[i];

        tr_x_xxxy_yyyyyy[i] = 6.0 * tr_x_xxx_yyyyy[i] * fe_0 + tr_x_xxx_yyyyyy[i] * pa_y[i];

        tr_x_xxxy_yyyyyz[i] = 5.0 * tr_x_xxx_yyyyz[i] * fe_0 + tr_x_xxx_yyyyyz[i] * pa_y[i];

        tr_x_xxxy_yyyyzz[i] = 4.0 * tr_x_xxx_yyyzz[i] * fe_0 + tr_x_xxx_yyyyzz[i] * pa_y[i];

        tr_x_xxxy_yyyzzz[i] = 3.0 * tr_x_xxx_yyzzz[i] * fe_0 + tr_x_xxx_yyyzzz[i] * pa_y[i];

        tr_x_xxxy_yyzzzz[i] = 2.0 * tr_x_xxx_yzzzz[i] * fe_0 + tr_x_xxx_yyzzzz[i] * pa_y[i];

        tr_x_xxxy_yzzzzz[i] = tr_x_xxx_zzzzz[i] * fe_0 + tr_x_xxx_yzzzzz[i] * pa_y[i];

        tr_x_xxxy_zzzzzz[i] = tr_x_xxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : GI

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

    auto tr_x_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 78);

    auto tr_x_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 79);

    auto tr_x_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 80);

    auto tr_x_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 81);

    auto tr_x_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 82);

    auto tr_x_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 83);

    #pragma omp simd aligned(pa_z, tr_x_xxx_xxxxx, tr_x_xxx_xxxxxx, tr_x_xxx_xxxxxy, tr_x_xxx_xxxxxz, tr_x_xxx_xxxxy, tr_x_xxx_xxxxyy, tr_x_xxx_xxxxyz, tr_x_xxx_xxxxz, tr_x_xxx_xxxxzz, tr_x_xxx_xxxyy, tr_x_xxx_xxxyyy, tr_x_xxx_xxxyyz, tr_x_xxx_xxxyz, tr_x_xxx_xxxyzz, tr_x_xxx_xxxzz, tr_x_xxx_xxxzzz, tr_x_xxx_xxyyy, tr_x_xxx_xxyyyy, tr_x_xxx_xxyyyz, tr_x_xxx_xxyyz, tr_x_xxx_xxyyzz, tr_x_xxx_xxyzz, tr_x_xxx_xxyzzz, tr_x_xxx_xxzzz, tr_x_xxx_xxzzzz, tr_x_xxx_xyyyy, tr_x_xxx_xyyyyy, tr_x_xxx_xyyyyz, tr_x_xxx_xyyyz, tr_x_xxx_xyyyzz, tr_x_xxx_xyyzz, tr_x_xxx_xyyzzz, tr_x_xxx_xyzzz, tr_x_xxx_xyzzzz, tr_x_xxx_xzzzz, tr_x_xxx_xzzzzz, tr_x_xxx_yyyyy, tr_x_xxx_yyyyyy, tr_x_xxx_yyyyyz, tr_x_xxx_yyyyz, tr_x_xxx_yyyyzz, tr_x_xxx_yyyzz, tr_x_xxx_yyyzzz, tr_x_xxx_yyzzz, tr_x_xxx_yyzzzz, tr_x_xxx_yzzzz, tr_x_xxx_yzzzzz, tr_x_xxx_zzzzz, tr_x_xxx_zzzzzz, tr_x_xxxz_xxxxxx, tr_x_xxxz_xxxxxy, tr_x_xxxz_xxxxxz, tr_x_xxxz_xxxxyy, tr_x_xxxz_xxxxyz, tr_x_xxxz_xxxxzz, tr_x_xxxz_xxxyyy, tr_x_xxxz_xxxyyz, tr_x_xxxz_xxxyzz, tr_x_xxxz_xxxzzz, tr_x_xxxz_xxyyyy, tr_x_xxxz_xxyyyz, tr_x_xxxz_xxyyzz, tr_x_xxxz_xxyzzz, tr_x_xxxz_xxzzzz, tr_x_xxxz_xyyyyy, tr_x_xxxz_xyyyyz, tr_x_xxxz_xyyyzz, tr_x_xxxz_xyyzzz, tr_x_xxxz_xyzzzz, tr_x_xxxz_xzzzzz, tr_x_xxxz_yyyyyy, tr_x_xxxz_yyyyyz, tr_x_xxxz_yyyyzz, tr_x_xxxz_yyyzzz, tr_x_xxxz_yyzzzz, tr_x_xxxz_yzzzzz, tr_x_xxxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxz_xxxxxx[i] = tr_x_xxx_xxxxxx[i] * pa_z[i];

        tr_x_xxxz_xxxxxy[i] = tr_x_xxx_xxxxxy[i] * pa_z[i];

        tr_x_xxxz_xxxxxz[i] = tr_x_xxx_xxxxx[i] * fe_0 + tr_x_xxx_xxxxxz[i] * pa_z[i];

        tr_x_xxxz_xxxxyy[i] = tr_x_xxx_xxxxyy[i] * pa_z[i];

        tr_x_xxxz_xxxxyz[i] = tr_x_xxx_xxxxy[i] * fe_0 + tr_x_xxx_xxxxyz[i] * pa_z[i];

        tr_x_xxxz_xxxxzz[i] = 2.0 * tr_x_xxx_xxxxz[i] * fe_0 + tr_x_xxx_xxxxzz[i] * pa_z[i];

        tr_x_xxxz_xxxyyy[i] = tr_x_xxx_xxxyyy[i] * pa_z[i];

        tr_x_xxxz_xxxyyz[i] = tr_x_xxx_xxxyy[i] * fe_0 + tr_x_xxx_xxxyyz[i] * pa_z[i];

        tr_x_xxxz_xxxyzz[i] = 2.0 * tr_x_xxx_xxxyz[i] * fe_0 + tr_x_xxx_xxxyzz[i] * pa_z[i];

        tr_x_xxxz_xxxzzz[i] = 3.0 * tr_x_xxx_xxxzz[i] * fe_0 + tr_x_xxx_xxxzzz[i] * pa_z[i];

        tr_x_xxxz_xxyyyy[i] = tr_x_xxx_xxyyyy[i] * pa_z[i];

        tr_x_xxxz_xxyyyz[i] = tr_x_xxx_xxyyy[i] * fe_0 + tr_x_xxx_xxyyyz[i] * pa_z[i];

        tr_x_xxxz_xxyyzz[i] = 2.0 * tr_x_xxx_xxyyz[i] * fe_0 + tr_x_xxx_xxyyzz[i] * pa_z[i];

        tr_x_xxxz_xxyzzz[i] = 3.0 * tr_x_xxx_xxyzz[i] * fe_0 + tr_x_xxx_xxyzzz[i] * pa_z[i];

        tr_x_xxxz_xxzzzz[i] = 4.0 * tr_x_xxx_xxzzz[i] * fe_0 + tr_x_xxx_xxzzzz[i] * pa_z[i];

        tr_x_xxxz_xyyyyy[i] = tr_x_xxx_xyyyyy[i] * pa_z[i];

        tr_x_xxxz_xyyyyz[i] = tr_x_xxx_xyyyy[i] * fe_0 + tr_x_xxx_xyyyyz[i] * pa_z[i];

        tr_x_xxxz_xyyyzz[i] = 2.0 * tr_x_xxx_xyyyz[i] * fe_0 + tr_x_xxx_xyyyzz[i] * pa_z[i];

        tr_x_xxxz_xyyzzz[i] = 3.0 * tr_x_xxx_xyyzz[i] * fe_0 + tr_x_xxx_xyyzzz[i] * pa_z[i];

        tr_x_xxxz_xyzzzz[i] = 4.0 * tr_x_xxx_xyzzz[i] * fe_0 + tr_x_xxx_xyzzzz[i] * pa_z[i];

        tr_x_xxxz_xzzzzz[i] = 5.0 * tr_x_xxx_xzzzz[i] * fe_0 + tr_x_xxx_xzzzzz[i] * pa_z[i];

        tr_x_xxxz_yyyyyy[i] = tr_x_xxx_yyyyyy[i] * pa_z[i];

        tr_x_xxxz_yyyyyz[i] = tr_x_xxx_yyyyy[i] * fe_0 + tr_x_xxx_yyyyyz[i] * pa_z[i];

        tr_x_xxxz_yyyyzz[i] = 2.0 * tr_x_xxx_yyyyz[i] * fe_0 + tr_x_xxx_yyyyzz[i] * pa_z[i];

        tr_x_xxxz_yyyzzz[i] = 3.0 * tr_x_xxx_yyyzz[i] * fe_0 + tr_x_xxx_yyyzzz[i] * pa_z[i];

        tr_x_xxxz_yyzzzz[i] = 4.0 * tr_x_xxx_yyzzz[i] * fe_0 + tr_x_xxx_yyzzzz[i] * pa_z[i];

        tr_x_xxxz_yzzzzz[i] = 5.0 * tr_x_xxx_yzzzz[i] * fe_0 + tr_x_xxx_yzzzzz[i] * pa_z[i];

        tr_x_xxxz_zzzzzz[i] = 6.0 * tr_x_xxx_zzzzz[i] * fe_0 + tr_x_xxx_zzzzzz[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xx_xxxxxx, tr_x_xx_xxxxxy, tr_x_xx_xxxxxz, tr_x_xx_xxxxyy, tr_x_xx_xxxxyz, tr_x_xx_xxxxzz, tr_x_xx_xxxyyy, tr_x_xx_xxxyyz, tr_x_xx_xxxyzz, tr_x_xx_xxxzzz, tr_x_xx_xxyyyy, tr_x_xx_xxyyyz, tr_x_xx_xxyyzz, tr_x_xx_xxyzzz, tr_x_xx_xxzzzz, tr_x_xx_xyyyyy, tr_x_xx_xyyyyz, tr_x_xx_xyyyzz, tr_x_xx_xyyzzz, tr_x_xx_xyzzzz, tr_x_xx_xzzzzz, tr_x_xx_zzzzzz, tr_x_xxy_xxxxx, tr_x_xxy_xxxxxx, tr_x_xxy_xxxxxy, tr_x_xxy_xxxxxz, tr_x_xxy_xxxxy, tr_x_xxy_xxxxyy, tr_x_xxy_xxxxyz, tr_x_xxy_xxxxz, tr_x_xxy_xxxxzz, tr_x_xxy_xxxyy, tr_x_xxy_xxxyyy, tr_x_xxy_xxxyyz, tr_x_xxy_xxxyz, tr_x_xxy_xxxyzz, tr_x_xxy_xxxzz, tr_x_xxy_xxxzzz, tr_x_xxy_xxyyy, tr_x_xxy_xxyyyy, tr_x_xxy_xxyyyz, tr_x_xxy_xxyyz, tr_x_xxy_xxyyzz, tr_x_xxy_xxyzz, tr_x_xxy_xxyzzz, tr_x_xxy_xxzzz, tr_x_xxy_xxzzzz, tr_x_xxy_xyyyy, tr_x_xxy_xyyyyy, tr_x_xxy_xyyyyz, tr_x_xxy_xyyyz, tr_x_xxy_xyyyzz, tr_x_xxy_xyyzz, tr_x_xxy_xyyzzz, tr_x_xxy_xyzzz, tr_x_xxy_xyzzzz, tr_x_xxy_xzzzz, tr_x_xxy_xzzzzz, tr_x_xxy_zzzzzz, tr_x_xxyy_xxxxxx, tr_x_xxyy_xxxxxy, tr_x_xxyy_xxxxxz, tr_x_xxyy_xxxxyy, tr_x_xxyy_xxxxyz, tr_x_xxyy_xxxxzz, tr_x_xxyy_xxxyyy, tr_x_xxyy_xxxyyz, tr_x_xxyy_xxxyzz, tr_x_xxyy_xxxzzz, tr_x_xxyy_xxyyyy, tr_x_xxyy_xxyyyz, tr_x_xxyy_xxyyzz, tr_x_xxyy_xxyzzz, tr_x_xxyy_xxzzzz, tr_x_xxyy_xyyyyy, tr_x_xxyy_xyyyyz, tr_x_xxyy_xyyyzz, tr_x_xxyy_xyyzzz, tr_x_xxyy_xyzzzz, tr_x_xxyy_xzzzzz, tr_x_xxyy_yyyyyy, tr_x_xxyy_yyyyyz, tr_x_xxyy_yyyyzz, tr_x_xxyy_yyyzzz, tr_x_xxyy_yyzzzz, tr_x_xxyy_yzzzzz, tr_x_xxyy_zzzzzz, tr_x_xyy_yyyyyy, tr_x_xyy_yyyyyz, tr_x_xyy_yyyyzz, tr_x_xyy_yyyzzz, tr_x_xyy_yyzzzz, tr_x_xyy_yzzzzz, tr_x_yy_yyyyyy, tr_x_yy_yyyyyz, tr_x_yy_yyyyzz, tr_x_yy_yyyzzz, tr_x_yy_yyzzzz, tr_x_yy_yzzzzz, ts_xyy_yyyyyy, ts_xyy_yyyyyz, ts_xyy_yyyyzz, ts_xyy_yyyzzz, ts_xyy_yyzzzz, ts_xyy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyy_xxxxxx[i] = tr_x_xx_xxxxxx[i] * fe_0 + tr_x_xxy_xxxxxx[i] * pa_y[i];

        tr_x_xxyy_xxxxxy[i] = tr_x_xx_xxxxxy[i] * fe_0 + tr_x_xxy_xxxxx[i] * fe_0 + tr_x_xxy_xxxxxy[i] * pa_y[i];

        tr_x_xxyy_xxxxxz[i] = tr_x_xx_xxxxxz[i] * fe_0 + tr_x_xxy_xxxxxz[i] * pa_y[i];

        tr_x_xxyy_xxxxyy[i] = tr_x_xx_xxxxyy[i] * fe_0 + 2.0 * tr_x_xxy_xxxxy[i] * fe_0 + tr_x_xxy_xxxxyy[i] * pa_y[i];

        tr_x_xxyy_xxxxyz[i] = tr_x_xx_xxxxyz[i] * fe_0 + tr_x_xxy_xxxxz[i] * fe_0 + tr_x_xxy_xxxxyz[i] * pa_y[i];

        tr_x_xxyy_xxxxzz[i] = tr_x_xx_xxxxzz[i] * fe_0 + tr_x_xxy_xxxxzz[i] * pa_y[i];

        tr_x_xxyy_xxxyyy[i] = tr_x_xx_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxy_xxxyy[i] * fe_0 + tr_x_xxy_xxxyyy[i] * pa_y[i];

        tr_x_xxyy_xxxyyz[i] = tr_x_xx_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxy_xxxyz[i] * fe_0 + tr_x_xxy_xxxyyz[i] * pa_y[i];

        tr_x_xxyy_xxxyzz[i] = tr_x_xx_xxxyzz[i] * fe_0 + tr_x_xxy_xxxzz[i] * fe_0 + tr_x_xxy_xxxyzz[i] * pa_y[i];

        tr_x_xxyy_xxxzzz[i] = tr_x_xx_xxxzzz[i] * fe_0 + tr_x_xxy_xxxzzz[i] * pa_y[i];

        tr_x_xxyy_xxyyyy[i] = tr_x_xx_xxyyyy[i] * fe_0 + 4.0 * tr_x_xxy_xxyyy[i] * fe_0 + tr_x_xxy_xxyyyy[i] * pa_y[i];

        tr_x_xxyy_xxyyyz[i] = tr_x_xx_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxy_xxyyz[i] * fe_0 + tr_x_xxy_xxyyyz[i] * pa_y[i];

        tr_x_xxyy_xxyyzz[i] = tr_x_xx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxy_xxyzz[i] * fe_0 + tr_x_xxy_xxyyzz[i] * pa_y[i];

        tr_x_xxyy_xxyzzz[i] = tr_x_xx_xxyzzz[i] * fe_0 + tr_x_xxy_xxzzz[i] * fe_0 + tr_x_xxy_xxyzzz[i] * pa_y[i];

        tr_x_xxyy_xxzzzz[i] = tr_x_xx_xxzzzz[i] * fe_0 + tr_x_xxy_xxzzzz[i] * pa_y[i];

        tr_x_xxyy_xyyyyy[i] = tr_x_xx_xyyyyy[i] * fe_0 + 5.0 * tr_x_xxy_xyyyy[i] * fe_0 + tr_x_xxy_xyyyyy[i] * pa_y[i];

        tr_x_xxyy_xyyyyz[i] = tr_x_xx_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxy_xyyyz[i] * fe_0 + tr_x_xxy_xyyyyz[i] * pa_y[i];

        tr_x_xxyy_xyyyzz[i] = tr_x_xx_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxy_xyyzz[i] * fe_0 + tr_x_xxy_xyyyzz[i] * pa_y[i];

        tr_x_xxyy_xyyzzz[i] = tr_x_xx_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxy_xyzzz[i] * fe_0 + tr_x_xxy_xyyzzz[i] * pa_y[i];

        tr_x_xxyy_xyzzzz[i] = tr_x_xx_xyzzzz[i] * fe_0 + tr_x_xxy_xzzzz[i] * fe_0 + tr_x_xxy_xyzzzz[i] * pa_y[i];

        tr_x_xxyy_xzzzzz[i] = tr_x_xx_xzzzzz[i] * fe_0 + tr_x_xxy_xzzzzz[i] * pa_y[i];

        tr_x_xxyy_yyyyyy[i] = tr_x_yy_yyyyyy[i] * fe_0 + ts_xyy_yyyyyy[i] * fe_0 + tr_x_xyy_yyyyyy[i] * pa_x[i];

        tr_x_xxyy_yyyyyz[i] = tr_x_yy_yyyyyz[i] * fe_0 + ts_xyy_yyyyyz[i] * fe_0 + tr_x_xyy_yyyyyz[i] * pa_x[i];

        tr_x_xxyy_yyyyzz[i] = tr_x_yy_yyyyzz[i] * fe_0 + ts_xyy_yyyyzz[i] * fe_0 + tr_x_xyy_yyyyzz[i] * pa_x[i];

        tr_x_xxyy_yyyzzz[i] = tr_x_yy_yyyzzz[i] * fe_0 + ts_xyy_yyyzzz[i] * fe_0 + tr_x_xyy_yyyzzz[i] * pa_x[i];

        tr_x_xxyy_yyzzzz[i] = tr_x_yy_yyzzzz[i] * fe_0 + ts_xyy_yyzzzz[i] * fe_0 + tr_x_xyy_yyzzzz[i] * pa_x[i];

        tr_x_xxyy_yzzzzz[i] = tr_x_yy_yzzzzz[i] * fe_0 + ts_xyy_yzzzzz[i] * fe_0 + tr_x_xyy_yzzzzz[i] * pa_x[i];

        tr_x_xxyy_zzzzzz[i] = tr_x_xx_zzzzzz[i] * fe_0 + tr_x_xxy_zzzzzz[i] * pa_y[i];
    }

    // Set up 112-140 components of targeted buffer : GI

    auto tr_x_xxyz_xxxxxx = pbuffer.data(idx_dip_gi + 112);

    auto tr_x_xxyz_xxxxxy = pbuffer.data(idx_dip_gi + 113);

    auto tr_x_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 114);

    auto tr_x_xxyz_xxxxyy = pbuffer.data(idx_dip_gi + 115);

    auto tr_x_xxyz_xxxxyz = pbuffer.data(idx_dip_gi + 116);

    auto tr_x_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 117);

    auto tr_x_xxyz_xxxyyy = pbuffer.data(idx_dip_gi + 118);

    auto tr_x_xxyz_xxxyyz = pbuffer.data(idx_dip_gi + 119);

    auto tr_x_xxyz_xxxyzz = pbuffer.data(idx_dip_gi + 120);

    auto tr_x_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 121);

    auto tr_x_xxyz_xxyyyy = pbuffer.data(idx_dip_gi + 122);

    auto tr_x_xxyz_xxyyyz = pbuffer.data(idx_dip_gi + 123);

    auto tr_x_xxyz_xxyyzz = pbuffer.data(idx_dip_gi + 124);

    auto tr_x_xxyz_xxyzzz = pbuffer.data(idx_dip_gi + 125);

    auto tr_x_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 126);

    auto tr_x_xxyz_xyyyyy = pbuffer.data(idx_dip_gi + 127);

    auto tr_x_xxyz_xyyyyz = pbuffer.data(idx_dip_gi + 128);

    auto tr_x_xxyz_xyyyzz = pbuffer.data(idx_dip_gi + 129);

    auto tr_x_xxyz_xyyzzz = pbuffer.data(idx_dip_gi + 130);

    auto tr_x_xxyz_xyzzzz = pbuffer.data(idx_dip_gi + 131);

    auto tr_x_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 132);

    auto tr_x_xxyz_yyyyyy = pbuffer.data(idx_dip_gi + 133);

    auto tr_x_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 134);

    auto tr_x_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 135);

    auto tr_x_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 136);

    auto tr_x_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 137);

    auto tr_x_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 138);

    auto tr_x_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 139);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxy_xxxxxy, tr_x_xxy_xxxxyy, tr_x_xxy_xxxyyy, tr_x_xxy_xxyyyy, tr_x_xxy_xyyyyy, tr_x_xxy_yyyyyy, tr_x_xxyz_xxxxxx, tr_x_xxyz_xxxxxy, tr_x_xxyz_xxxxxz, tr_x_xxyz_xxxxyy, tr_x_xxyz_xxxxyz, tr_x_xxyz_xxxxzz, tr_x_xxyz_xxxyyy, tr_x_xxyz_xxxyyz, tr_x_xxyz_xxxyzz, tr_x_xxyz_xxxzzz, tr_x_xxyz_xxyyyy, tr_x_xxyz_xxyyyz, tr_x_xxyz_xxyyzz, tr_x_xxyz_xxyzzz, tr_x_xxyz_xxzzzz, tr_x_xxyz_xyyyyy, tr_x_xxyz_xyyyyz, tr_x_xxyz_xyyyzz, tr_x_xxyz_xyyzzz, tr_x_xxyz_xyzzzz, tr_x_xxyz_xzzzzz, tr_x_xxyz_yyyyyy, tr_x_xxyz_yyyyyz, tr_x_xxyz_yyyyzz, tr_x_xxyz_yyyzzz, tr_x_xxyz_yyzzzz, tr_x_xxyz_yzzzzz, tr_x_xxyz_zzzzzz, tr_x_xxz_xxxxxx, tr_x_xxz_xxxxxz, tr_x_xxz_xxxxyz, tr_x_xxz_xxxxz, tr_x_xxz_xxxxzz, tr_x_xxz_xxxyyz, tr_x_xxz_xxxyz, tr_x_xxz_xxxyzz, tr_x_xxz_xxxzz, tr_x_xxz_xxxzzz, tr_x_xxz_xxyyyz, tr_x_xxz_xxyyz, tr_x_xxz_xxyyzz, tr_x_xxz_xxyzz, tr_x_xxz_xxyzzz, tr_x_xxz_xxzzz, tr_x_xxz_xxzzzz, tr_x_xxz_xyyyyz, tr_x_xxz_xyyyz, tr_x_xxz_xyyyzz, tr_x_xxz_xyyzz, tr_x_xxz_xyyzzz, tr_x_xxz_xyzzz, tr_x_xxz_xyzzzz, tr_x_xxz_xzzzz, tr_x_xxz_xzzzzz, tr_x_xxz_yyyyyz, tr_x_xxz_yyyyz, tr_x_xxz_yyyyzz, tr_x_xxz_yyyzz, tr_x_xxz_yyyzzz, tr_x_xxz_yyzzz, tr_x_xxz_yyzzzz, tr_x_xxz_yzzzz, tr_x_xxz_yzzzzz, tr_x_xxz_zzzzz, tr_x_xxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyz_xxxxxx[i] = tr_x_xxz_xxxxxx[i] * pa_y[i];

        tr_x_xxyz_xxxxxy[i] = tr_x_xxy_xxxxxy[i] * pa_z[i];

        tr_x_xxyz_xxxxxz[i] = tr_x_xxz_xxxxxz[i] * pa_y[i];

        tr_x_xxyz_xxxxyy[i] = tr_x_xxy_xxxxyy[i] * pa_z[i];

        tr_x_xxyz_xxxxyz[i] = tr_x_xxz_xxxxz[i] * fe_0 + tr_x_xxz_xxxxyz[i] * pa_y[i];

        tr_x_xxyz_xxxxzz[i] = tr_x_xxz_xxxxzz[i] * pa_y[i];

        tr_x_xxyz_xxxyyy[i] = tr_x_xxy_xxxyyy[i] * pa_z[i];

        tr_x_xxyz_xxxyyz[i] = 2.0 * tr_x_xxz_xxxyz[i] * fe_0 + tr_x_xxz_xxxyyz[i] * pa_y[i];

        tr_x_xxyz_xxxyzz[i] = tr_x_xxz_xxxzz[i] * fe_0 + tr_x_xxz_xxxyzz[i] * pa_y[i];

        tr_x_xxyz_xxxzzz[i] = tr_x_xxz_xxxzzz[i] * pa_y[i];

        tr_x_xxyz_xxyyyy[i] = tr_x_xxy_xxyyyy[i] * pa_z[i];

        tr_x_xxyz_xxyyyz[i] = 3.0 * tr_x_xxz_xxyyz[i] * fe_0 + tr_x_xxz_xxyyyz[i] * pa_y[i];

        tr_x_xxyz_xxyyzz[i] = 2.0 * tr_x_xxz_xxyzz[i] * fe_0 + tr_x_xxz_xxyyzz[i] * pa_y[i];

        tr_x_xxyz_xxyzzz[i] = tr_x_xxz_xxzzz[i] * fe_0 + tr_x_xxz_xxyzzz[i] * pa_y[i];

        tr_x_xxyz_xxzzzz[i] = tr_x_xxz_xxzzzz[i] * pa_y[i];

        tr_x_xxyz_xyyyyy[i] = tr_x_xxy_xyyyyy[i] * pa_z[i];

        tr_x_xxyz_xyyyyz[i] = 4.0 * tr_x_xxz_xyyyz[i] * fe_0 + tr_x_xxz_xyyyyz[i] * pa_y[i];

        tr_x_xxyz_xyyyzz[i] = 3.0 * tr_x_xxz_xyyzz[i] * fe_0 + tr_x_xxz_xyyyzz[i] * pa_y[i];

        tr_x_xxyz_xyyzzz[i] = 2.0 * tr_x_xxz_xyzzz[i] * fe_0 + tr_x_xxz_xyyzzz[i] * pa_y[i];

        tr_x_xxyz_xyzzzz[i] = tr_x_xxz_xzzzz[i] * fe_0 + tr_x_xxz_xyzzzz[i] * pa_y[i];

        tr_x_xxyz_xzzzzz[i] = tr_x_xxz_xzzzzz[i] * pa_y[i];

        tr_x_xxyz_yyyyyy[i] = tr_x_xxy_yyyyyy[i] * pa_z[i];

        tr_x_xxyz_yyyyyz[i] = 5.0 * tr_x_xxz_yyyyz[i] * fe_0 + tr_x_xxz_yyyyyz[i] * pa_y[i];

        tr_x_xxyz_yyyyzz[i] = 4.0 * tr_x_xxz_yyyzz[i] * fe_0 + tr_x_xxz_yyyyzz[i] * pa_y[i];

        tr_x_xxyz_yyyzzz[i] = 3.0 * tr_x_xxz_yyzzz[i] * fe_0 + tr_x_xxz_yyyzzz[i] * pa_y[i];

        tr_x_xxyz_yyzzzz[i] = 2.0 * tr_x_xxz_yzzzz[i] * fe_0 + tr_x_xxz_yyzzzz[i] * pa_y[i];

        tr_x_xxyz_yzzzzz[i] = tr_x_xxz_zzzzz[i] * fe_0 + tr_x_xxz_yzzzzz[i] * pa_y[i];

        tr_x_xxyz_zzzzzz[i] = tr_x_xxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xx_xxxxxx, tr_x_xx_xxxxxy, tr_x_xx_xxxxxz, tr_x_xx_xxxxyy, tr_x_xx_xxxxyz, tr_x_xx_xxxxzz, tr_x_xx_xxxyyy, tr_x_xx_xxxyyz, tr_x_xx_xxxyzz, tr_x_xx_xxxzzz, tr_x_xx_xxyyyy, tr_x_xx_xxyyyz, tr_x_xx_xxyyzz, tr_x_xx_xxyzzz, tr_x_xx_xxzzzz, tr_x_xx_xyyyyy, tr_x_xx_xyyyyz, tr_x_xx_xyyyzz, tr_x_xx_xyyzzz, tr_x_xx_xyzzzz, tr_x_xx_xzzzzz, tr_x_xx_yyyyyy, tr_x_xxz_xxxxx, tr_x_xxz_xxxxxx, tr_x_xxz_xxxxxy, tr_x_xxz_xxxxxz, tr_x_xxz_xxxxy, tr_x_xxz_xxxxyy, tr_x_xxz_xxxxyz, tr_x_xxz_xxxxz, tr_x_xxz_xxxxzz, tr_x_xxz_xxxyy, tr_x_xxz_xxxyyy, tr_x_xxz_xxxyyz, tr_x_xxz_xxxyz, tr_x_xxz_xxxyzz, tr_x_xxz_xxxzz, tr_x_xxz_xxxzzz, tr_x_xxz_xxyyy, tr_x_xxz_xxyyyy, tr_x_xxz_xxyyyz, tr_x_xxz_xxyyz, tr_x_xxz_xxyyzz, tr_x_xxz_xxyzz, tr_x_xxz_xxyzzz, tr_x_xxz_xxzzz, tr_x_xxz_xxzzzz, tr_x_xxz_xyyyy, tr_x_xxz_xyyyyy, tr_x_xxz_xyyyyz, tr_x_xxz_xyyyz, tr_x_xxz_xyyyzz, tr_x_xxz_xyyzz, tr_x_xxz_xyyzzz, tr_x_xxz_xyzzz, tr_x_xxz_xyzzzz, tr_x_xxz_xzzzz, tr_x_xxz_xzzzzz, tr_x_xxz_yyyyyy, tr_x_xxzz_xxxxxx, tr_x_xxzz_xxxxxy, tr_x_xxzz_xxxxxz, tr_x_xxzz_xxxxyy, tr_x_xxzz_xxxxyz, tr_x_xxzz_xxxxzz, tr_x_xxzz_xxxyyy, tr_x_xxzz_xxxyyz, tr_x_xxzz_xxxyzz, tr_x_xxzz_xxxzzz, tr_x_xxzz_xxyyyy, tr_x_xxzz_xxyyyz, tr_x_xxzz_xxyyzz, tr_x_xxzz_xxyzzz, tr_x_xxzz_xxzzzz, tr_x_xxzz_xyyyyy, tr_x_xxzz_xyyyyz, tr_x_xxzz_xyyyzz, tr_x_xxzz_xyyzzz, tr_x_xxzz_xyzzzz, tr_x_xxzz_xzzzzz, tr_x_xxzz_yyyyyy, tr_x_xxzz_yyyyyz, tr_x_xxzz_yyyyzz, tr_x_xxzz_yyyzzz, tr_x_xxzz_yyzzzz, tr_x_xxzz_yzzzzz, tr_x_xxzz_zzzzzz, tr_x_xzz_yyyyyz, tr_x_xzz_yyyyzz, tr_x_xzz_yyyzzz, tr_x_xzz_yyzzzz, tr_x_xzz_yzzzzz, tr_x_xzz_zzzzzz, tr_x_zz_yyyyyz, tr_x_zz_yyyyzz, tr_x_zz_yyyzzz, tr_x_zz_yyzzzz, tr_x_zz_yzzzzz, tr_x_zz_zzzzzz, ts_xzz_yyyyyz, ts_xzz_yyyyzz, ts_xzz_yyyzzz, ts_xzz_yyzzzz, ts_xzz_yzzzzz, ts_xzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzz_xxxxxx[i] = tr_x_xx_xxxxxx[i] * fe_0 + tr_x_xxz_xxxxxx[i] * pa_z[i];

        tr_x_xxzz_xxxxxy[i] = tr_x_xx_xxxxxy[i] * fe_0 + tr_x_xxz_xxxxxy[i] * pa_z[i];

        tr_x_xxzz_xxxxxz[i] = tr_x_xx_xxxxxz[i] * fe_0 + tr_x_xxz_xxxxx[i] * fe_0 + tr_x_xxz_xxxxxz[i] * pa_z[i];

        tr_x_xxzz_xxxxyy[i] = tr_x_xx_xxxxyy[i] * fe_0 + tr_x_xxz_xxxxyy[i] * pa_z[i];

        tr_x_xxzz_xxxxyz[i] = tr_x_xx_xxxxyz[i] * fe_0 + tr_x_xxz_xxxxy[i] * fe_0 + tr_x_xxz_xxxxyz[i] * pa_z[i];

        tr_x_xxzz_xxxxzz[i] = tr_x_xx_xxxxzz[i] * fe_0 + 2.0 * tr_x_xxz_xxxxz[i] * fe_0 + tr_x_xxz_xxxxzz[i] * pa_z[i];

        tr_x_xxzz_xxxyyy[i] = tr_x_xx_xxxyyy[i] * fe_0 + tr_x_xxz_xxxyyy[i] * pa_z[i];

        tr_x_xxzz_xxxyyz[i] = tr_x_xx_xxxyyz[i] * fe_0 + tr_x_xxz_xxxyy[i] * fe_0 + tr_x_xxz_xxxyyz[i] * pa_z[i];

        tr_x_xxzz_xxxyzz[i] = tr_x_xx_xxxyzz[i] * fe_0 + 2.0 * tr_x_xxz_xxxyz[i] * fe_0 + tr_x_xxz_xxxyzz[i] * pa_z[i];

        tr_x_xxzz_xxxzzz[i] = tr_x_xx_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxz_xxxzz[i] * fe_0 + tr_x_xxz_xxxzzz[i] * pa_z[i];

        tr_x_xxzz_xxyyyy[i] = tr_x_xx_xxyyyy[i] * fe_0 + tr_x_xxz_xxyyyy[i] * pa_z[i];

        tr_x_xxzz_xxyyyz[i] = tr_x_xx_xxyyyz[i] * fe_0 + tr_x_xxz_xxyyy[i] * fe_0 + tr_x_xxz_xxyyyz[i] * pa_z[i];

        tr_x_xxzz_xxyyzz[i] = tr_x_xx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxz_xxyyz[i] * fe_0 + tr_x_xxz_xxyyzz[i] * pa_z[i];

        tr_x_xxzz_xxyzzz[i] = tr_x_xx_xxyzzz[i] * fe_0 + 3.0 * tr_x_xxz_xxyzz[i] * fe_0 + tr_x_xxz_xxyzzz[i] * pa_z[i];

        tr_x_xxzz_xxzzzz[i] = tr_x_xx_xxzzzz[i] * fe_0 + 4.0 * tr_x_xxz_xxzzz[i] * fe_0 + tr_x_xxz_xxzzzz[i] * pa_z[i];

        tr_x_xxzz_xyyyyy[i] = tr_x_xx_xyyyyy[i] * fe_0 + tr_x_xxz_xyyyyy[i] * pa_z[i];

        tr_x_xxzz_xyyyyz[i] = tr_x_xx_xyyyyz[i] * fe_0 + tr_x_xxz_xyyyy[i] * fe_0 + tr_x_xxz_xyyyyz[i] * pa_z[i];

        tr_x_xxzz_xyyyzz[i] = tr_x_xx_xyyyzz[i] * fe_0 + 2.0 * tr_x_xxz_xyyyz[i] * fe_0 + tr_x_xxz_xyyyzz[i] * pa_z[i];

        tr_x_xxzz_xyyzzz[i] = tr_x_xx_xyyzzz[i] * fe_0 + 3.0 * tr_x_xxz_xyyzz[i] * fe_0 + tr_x_xxz_xyyzzz[i] * pa_z[i];

        tr_x_xxzz_xyzzzz[i] = tr_x_xx_xyzzzz[i] * fe_0 + 4.0 * tr_x_xxz_xyzzz[i] * fe_0 + tr_x_xxz_xyzzzz[i] * pa_z[i];

        tr_x_xxzz_xzzzzz[i] = tr_x_xx_xzzzzz[i] * fe_0 + 5.0 * tr_x_xxz_xzzzz[i] * fe_0 + tr_x_xxz_xzzzzz[i] * pa_z[i];

        tr_x_xxzz_yyyyyy[i] = tr_x_xx_yyyyyy[i] * fe_0 + tr_x_xxz_yyyyyy[i] * pa_z[i];

        tr_x_xxzz_yyyyyz[i] = tr_x_zz_yyyyyz[i] * fe_0 + ts_xzz_yyyyyz[i] * fe_0 + tr_x_xzz_yyyyyz[i] * pa_x[i];

        tr_x_xxzz_yyyyzz[i] = tr_x_zz_yyyyzz[i] * fe_0 + ts_xzz_yyyyzz[i] * fe_0 + tr_x_xzz_yyyyzz[i] * pa_x[i];

        tr_x_xxzz_yyyzzz[i] = tr_x_zz_yyyzzz[i] * fe_0 + ts_xzz_yyyzzz[i] * fe_0 + tr_x_xzz_yyyzzz[i] * pa_x[i];

        tr_x_xxzz_yyzzzz[i] = tr_x_zz_yyzzzz[i] * fe_0 + ts_xzz_yyzzzz[i] * fe_0 + tr_x_xzz_yyzzzz[i] * pa_x[i];

        tr_x_xxzz_yzzzzz[i] = tr_x_zz_yzzzzz[i] * fe_0 + ts_xzz_yzzzzz[i] * fe_0 + tr_x_xzz_yzzzzz[i] * pa_x[i];

        tr_x_xxzz_zzzzzz[i] = tr_x_zz_zzzzzz[i] * fe_0 + ts_xzz_zzzzzz[i] * fe_0 + tr_x_xzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : GI

    auto tr_x_xyyy_xxxxxx = pbuffer.data(idx_dip_gi + 168);

    auto tr_x_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 169);

    auto tr_x_xyyy_xxxxxz = pbuffer.data(idx_dip_gi + 170);

    auto tr_x_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 171);

    auto tr_x_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 172);

    auto tr_x_xyyy_xxxxzz = pbuffer.data(idx_dip_gi + 173);

    auto tr_x_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 174);

    auto tr_x_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 175);

    auto tr_x_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 176);

    auto tr_x_xyyy_xxxzzz = pbuffer.data(idx_dip_gi + 177);

    auto tr_x_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 178);

    auto tr_x_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 179);

    auto tr_x_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 180);

    auto tr_x_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 181);

    auto tr_x_xyyy_xxzzzz = pbuffer.data(idx_dip_gi + 182);

    auto tr_x_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 183);

    auto tr_x_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 184);

    auto tr_x_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 185);

    auto tr_x_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 186);

    auto tr_x_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 187);

    auto tr_x_xyyy_xzzzzz = pbuffer.data(idx_dip_gi + 188);

    auto tr_x_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 189);

    auto tr_x_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 190);

    auto tr_x_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 191);

    auto tr_x_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 192);

    auto tr_x_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 193);

    auto tr_x_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 194);

    auto tr_x_xyyy_zzzzzz = pbuffer.data(idx_dip_gi + 195);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xy_xxxxxx, tr_x_xy_xxxxxz, tr_x_xy_xxxxzz, tr_x_xy_xxxzzz, tr_x_xy_xxzzzz, tr_x_xy_xzzzzz, tr_x_xyy_xxxxxx, tr_x_xyy_xxxxxz, tr_x_xyy_xxxxzz, tr_x_xyy_xxxzzz, tr_x_xyy_xxzzzz, tr_x_xyy_xzzzzz, tr_x_xyyy_xxxxxx, tr_x_xyyy_xxxxxy, tr_x_xyyy_xxxxxz, tr_x_xyyy_xxxxyy, tr_x_xyyy_xxxxyz, tr_x_xyyy_xxxxzz, tr_x_xyyy_xxxyyy, tr_x_xyyy_xxxyyz, tr_x_xyyy_xxxyzz, tr_x_xyyy_xxxzzz, tr_x_xyyy_xxyyyy, tr_x_xyyy_xxyyyz, tr_x_xyyy_xxyyzz, tr_x_xyyy_xxyzzz, tr_x_xyyy_xxzzzz, tr_x_xyyy_xyyyyy, tr_x_xyyy_xyyyyz, tr_x_xyyy_xyyyzz, tr_x_xyyy_xyyzzz, tr_x_xyyy_xyzzzz, tr_x_xyyy_xzzzzz, tr_x_xyyy_yyyyyy, tr_x_xyyy_yyyyyz, tr_x_xyyy_yyyyzz, tr_x_xyyy_yyyzzz, tr_x_xyyy_yyzzzz, tr_x_xyyy_yzzzzz, tr_x_xyyy_zzzzzz, tr_x_yyy_xxxxxy, tr_x_yyy_xxxxy, tr_x_yyy_xxxxyy, tr_x_yyy_xxxxyz, tr_x_yyy_xxxyy, tr_x_yyy_xxxyyy, tr_x_yyy_xxxyyz, tr_x_yyy_xxxyz, tr_x_yyy_xxxyzz, tr_x_yyy_xxyyy, tr_x_yyy_xxyyyy, tr_x_yyy_xxyyyz, tr_x_yyy_xxyyz, tr_x_yyy_xxyyzz, tr_x_yyy_xxyzz, tr_x_yyy_xxyzzz, tr_x_yyy_xyyyy, tr_x_yyy_xyyyyy, tr_x_yyy_xyyyyz, tr_x_yyy_xyyyz, tr_x_yyy_xyyyzz, tr_x_yyy_xyyzz, tr_x_yyy_xyyzzz, tr_x_yyy_xyzzz, tr_x_yyy_xyzzzz, tr_x_yyy_yyyyy, tr_x_yyy_yyyyyy, tr_x_yyy_yyyyyz, tr_x_yyy_yyyyz, tr_x_yyy_yyyyzz, tr_x_yyy_yyyzz, tr_x_yyy_yyyzzz, tr_x_yyy_yyzzz, tr_x_yyy_yyzzzz, tr_x_yyy_yzzzz, tr_x_yyy_yzzzzz, tr_x_yyy_zzzzzz, ts_yyy_xxxxxy, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyzz, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyzz, ts_yyy_xxyzzz, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzzz, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzzz, ts_yyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyy_xxxxxx[i] = 2.0 * tr_x_xy_xxxxxx[i] * fe_0 + tr_x_xyy_xxxxxx[i] * pa_y[i];

        tr_x_xyyy_xxxxxy[i] = 5.0 * tr_x_yyy_xxxxy[i] * fe_0 + ts_yyy_xxxxxy[i] * fe_0 + tr_x_yyy_xxxxxy[i] * pa_x[i];

        tr_x_xyyy_xxxxxz[i] = 2.0 * tr_x_xy_xxxxxz[i] * fe_0 + tr_x_xyy_xxxxxz[i] * pa_y[i];

        tr_x_xyyy_xxxxyy[i] = 4.0 * tr_x_yyy_xxxyy[i] * fe_0 + ts_yyy_xxxxyy[i] * fe_0 + tr_x_yyy_xxxxyy[i] * pa_x[i];

        tr_x_xyyy_xxxxyz[i] = 4.0 * tr_x_yyy_xxxyz[i] * fe_0 + ts_yyy_xxxxyz[i] * fe_0 + tr_x_yyy_xxxxyz[i] * pa_x[i];

        tr_x_xyyy_xxxxzz[i] = 2.0 * tr_x_xy_xxxxzz[i] * fe_0 + tr_x_xyy_xxxxzz[i] * pa_y[i];

        tr_x_xyyy_xxxyyy[i] = 3.0 * tr_x_yyy_xxyyy[i] * fe_0 + ts_yyy_xxxyyy[i] * fe_0 + tr_x_yyy_xxxyyy[i] * pa_x[i];

        tr_x_xyyy_xxxyyz[i] = 3.0 * tr_x_yyy_xxyyz[i] * fe_0 + ts_yyy_xxxyyz[i] * fe_0 + tr_x_yyy_xxxyyz[i] * pa_x[i];

        tr_x_xyyy_xxxyzz[i] = 3.0 * tr_x_yyy_xxyzz[i] * fe_0 + ts_yyy_xxxyzz[i] * fe_0 + tr_x_yyy_xxxyzz[i] * pa_x[i];

        tr_x_xyyy_xxxzzz[i] = 2.0 * tr_x_xy_xxxzzz[i] * fe_0 + tr_x_xyy_xxxzzz[i] * pa_y[i];

        tr_x_xyyy_xxyyyy[i] = 2.0 * tr_x_yyy_xyyyy[i] * fe_0 + ts_yyy_xxyyyy[i] * fe_0 + tr_x_yyy_xxyyyy[i] * pa_x[i];

        tr_x_xyyy_xxyyyz[i] = 2.0 * tr_x_yyy_xyyyz[i] * fe_0 + ts_yyy_xxyyyz[i] * fe_0 + tr_x_yyy_xxyyyz[i] * pa_x[i];

        tr_x_xyyy_xxyyzz[i] = 2.0 * tr_x_yyy_xyyzz[i] * fe_0 + ts_yyy_xxyyzz[i] * fe_0 + tr_x_yyy_xxyyzz[i] * pa_x[i];

        tr_x_xyyy_xxyzzz[i] = 2.0 * tr_x_yyy_xyzzz[i] * fe_0 + ts_yyy_xxyzzz[i] * fe_0 + tr_x_yyy_xxyzzz[i] * pa_x[i];

        tr_x_xyyy_xxzzzz[i] = 2.0 * tr_x_xy_xxzzzz[i] * fe_0 + tr_x_xyy_xxzzzz[i] * pa_y[i];

        tr_x_xyyy_xyyyyy[i] = tr_x_yyy_yyyyy[i] * fe_0 + ts_yyy_xyyyyy[i] * fe_0 + tr_x_yyy_xyyyyy[i] * pa_x[i];

        tr_x_xyyy_xyyyyz[i] = tr_x_yyy_yyyyz[i] * fe_0 + ts_yyy_xyyyyz[i] * fe_0 + tr_x_yyy_xyyyyz[i] * pa_x[i];

        tr_x_xyyy_xyyyzz[i] = tr_x_yyy_yyyzz[i] * fe_0 + ts_yyy_xyyyzz[i] * fe_0 + tr_x_yyy_xyyyzz[i] * pa_x[i];

        tr_x_xyyy_xyyzzz[i] = tr_x_yyy_yyzzz[i] * fe_0 + ts_yyy_xyyzzz[i] * fe_0 + tr_x_yyy_xyyzzz[i] * pa_x[i];

        tr_x_xyyy_xyzzzz[i] = tr_x_yyy_yzzzz[i] * fe_0 + ts_yyy_xyzzzz[i] * fe_0 + tr_x_yyy_xyzzzz[i] * pa_x[i];

        tr_x_xyyy_xzzzzz[i] = 2.0 * tr_x_xy_xzzzzz[i] * fe_0 + tr_x_xyy_xzzzzz[i] * pa_y[i];

        tr_x_xyyy_yyyyyy[i] = ts_yyy_yyyyyy[i] * fe_0 + tr_x_yyy_yyyyyy[i] * pa_x[i];

        tr_x_xyyy_yyyyyz[i] = ts_yyy_yyyyyz[i] * fe_0 + tr_x_yyy_yyyyyz[i] * pa_x[i];

        tr_x_xyyy_yyyyzz[i] = ts_yyy_yyyyzz[i] * fe_0 + tr_x_yyy_yyyyzz[i] * pa_x[i];

        tr_x_xyyy_yyyzzz[i] = ts_yyy_yyyzzz[i] * fe_0 + tr_x_yyy_yyyzzz[i] * pa_x[i];

        tr_x_xyyy_yyzzzz[i] = ts_yyy_yyzzzz[i] * fe_0 + tr_x_yyy_yyzzzz[i] * pa_x[i];

        tr_x_xyyy_yzzzzz[i] = ts_yyy_yzzzzz[i] * fe_0 + tr_x_yyy_yzzzzz[i] * pa_x[i];

        tr_x_xyyy_zzzzzz[i] = ts_yyy_zzzzzz[i] * fe_0 + tr_x_yyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : GI

    auto tr_x_xyyz_xxxxxx = pbuffer.data(idx_dip_gi + 196);

    auto tr_x_xyyz_xxxxxy = pbuffer.data(idx_dip_gi + 197);

    auto tr_x_xyyz_xxxxxz = pbuffer.data(idx_dip_gi + 198);

    auto tr_x_xyyz_xxxxyy = pbuffer.data(idx_dip_gi + 199);

    auto tr_x_xyyz_xxxxyz = pbuffer.data(idx_dip_gi + 200);

    auto tr_x_xyyz_xxxxzz = pbuffer.data(idx_dip_gi + 201);

    auto tr_x_xyyz_xxxyyy = pbuffer.data(idx_dip_gi + 202);

    auto tr_x_xyyz_xxxyyz = pbuffer.data(idx_dip_gi + 203);

    auto tr_x_xyyz_xxxyzz = pbuffer.data(idx_dip_gi + 204);

    auto tr_x_xyyz_xxxzzz = pbuffer.data(idx_dip_gi + 205);

    auto tr_x_xyyz_xxyyyy = pbuffer.data(idx_dip_gi + 206);

    auto tr_x_xyyz_xxyyyz = pbuffer.data(idx_dip_gi + 207);

    auto tr_x_xyyz_xxyyzz = pbuffer.data(idx_dip_gi + 208);

    auto tr_x_xyyz_xxyzzz = pbuffer.data(idx_dip_gi + 209);

    auto tr_x_xyyz_xxzzzz = pbuffer.data(idx_dip_gi + 210);

    auto tr_x_xyyz_xyyyyy = pbuffer.data(idx_dip_gi + 211);

    auto tr_x_xyyz_xyyyyz = pbuffer.data(idx_dip_gi + 212);

    auto tr_x_xyyz_xyyyzz = pbuffer.data(idx_dip_gi + 213);

    auto tr_x_xyyz_xyyzzz = pbuffer.data(idx_dip_gi + 214);

    auto tr_x_xyyz_xyzzzz = pbuffer.data(idx_dip_gi + 215);

    auto tr_x_xyyz_xzzzzz = pbuffer.data(idx_dip_gi + 216);

    auto tr_x_xyyz_yyyyyy = pbuffer.data(idx_dip_gi + 217);

    auto tr_x_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 218);

    auto tr_x_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 219);

    auto tr_x_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 220);

    auto tr_x_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 221);

    auto tr_x_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 222);

    auto tr_x_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 223);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyy_xxxxxx, tr_x_xyy_xxxxxy, tr_x_xyy_xxxxy, tr_x_xyy_xxxxyy, tr_x_xyy_xxxxyz, tr_x_xyy_xxxyy, tr_x_xyy_xxxyyy, tr_x_xyy_xxxyyz, tr_x_xyy_xxxyz, tr_x_xyy_xxxyzz, tr_x_xyy_xxyyy, tr_x_xyy_xxyyyy, tr_x_xyy_xxyyyz, tr_x_xyy_xxyyz, tr_x_xyy_xxyyzz, tr_x_xyy_xxyzz, tr_x_xyy_xxyzzz, tr_x_xyy_xyyyy, tr_x_xyy_xyyyyy, tr_x_xyy_xyyyyz, tr_x_xyy_xyyyz, tr_x_xyy_xyyyzz, tr_x_xyy_xyyzz, tr_x_xyy_xyyzzz, tr_x_xyy_xyzzz, tr_x_xyy_xyzzzz, tr_x_xyy_yyyyyy, tr_x_xyyz_xxxxxx, tr_x_xyyz_xxxxxy, tr_x_xyyz_xxxxxz, tr_x_xyyz_xxxxyy, tr_x_xyyz_xxxxyz, tr_x_xyyz_xxxxzz, tr_x_xyyz_xxxyyy, tr_x_xyyz_xxxyyz, tr_x_xyyz_xxxyzz, tr_x_xyyz_xxxzzz, tr_x_xyyz_xxyyyy, tr_x_xyyz_xxyyyz, tr_x_xyyz_xxyyzz, tr_x_xyyz_xxyzzz, tr_x_xyyz_xxzzzz, tr_x_xyyz_xyyyyy, tr_x_xyyz_xyyyyz, tr_x_xyyz_xyyyzz, tr_x_xyyz_xyyzzz, tr_x_xyyz_xyzzzz, tr_x_xyyz_xzzzzz, tr_x_xyyz_yyyyyy, tr_x_xyyz_yyyyyz, tr_x_xyyz_yyyyzz, tr_x_xyyz_yyyzzz, tr_x_xyyz_yyzzzz, tr_x_xyyz_yzzzzz, tr_x_xyyz_zzzzzz, tr_x_xyz_xxxxxz, tr_x_xyz_xxxxzz, tr_x_xyz_xxxzzz, tr_x_xyz_xxzzzz, tr_x_xyz_xzzzzz, tr_x_xz_xxxxxz, tr_x_xz_xxxxzz, tr_x_xz_xxxzzz, tr_x_xz_xxzzzz, tr_x_xz_xzzzzz, tr_x_yyz_yyyyyz, tr_x_yyz_yyyyzz, tr_x_yyz_yyyzzz, tr_x_yyz_yyzzzz, tr_x_yyz_yzzzzz, tr_x_yyz_zzzzzz, ts_yyz_yyyyyz, ts_yyz_yyyyzz, ts_yyz_yyyzzz, ts_yyz_yyzzzz, ts_yyz_yzzzzz, ts_yyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyz_xxxxxx[i] = tr_x_xyy_xxxxxx[i] * pa_z[i];

        tr_x_xyyz_xxxxxy[i] = tr_x_xyy_xxxxxy[i] * pa_z[i];

        tr_x_xyyz_xxxxxz[i] = tr_x_xz_xxxxxz[i] * fe_0 + tr_x_xyz_xxxxxz[i] * pa_y[i];

        tr_x_xyyz_xxxxyy[i] = tr_x_xyy_xxxxyy[i] * pa_z[i];

        tr_x_xyyz_xxxxyz[i] = tr_x_xyy_xxxxy[i] * fe_0 + tr_x_xyy_xxxxyz[i] * pa_z[i];

        tr_x_xyyz_xxxxzz[i] = tr_x_xz_xxxxzz[i] * fe_0 + tr_x_xyz_xxxxzz[i] * pa_y[i];

        tr_x_xyyz_xxxyyy[i] = tr_x_xyy_xxxyyy[i] * pa_z[i];

        tr_x_xyyz_xxxyyz[i] = tr_x_xyy_xxxyy[i] * fe_0 + tr_x_xyy_xxxyyz[i] * pa_z[i];

        tr_x_xyyz_xxxyzz[i] = 2.0 * tr_x_xyy_xxxyz[i] * fe_0 + tr_x_xyy_xxxyzz[i] * pa_z[i];

        tr_x_xyyz_xxxzzz[i] = tr_x_xz_xxxzzz[i] * fe_0 + tr_x_xyz_xxxzzz[i] * pa_y[i];

        tr_x_xyyz_xxyyyy[i] = tr_x_xyy_xxyyyy[i] * pa_z[i];

        tr_x_xyyz_xxyyyz[i] = tr_x_xyy_xxyyy[i] * fe_0 + tr_x_xyy_xxyyyz[i] * pa_z[i];

        tr_x_xyyz_xxyyzz[i] = 2.0 * tr_x_xyy_xxyyz[i] * fe_0 + tr_x_xyy_xxyyzz[i] * pa_z[i];

        tr_x_xyyz_xxyzzz[i] = 3.0 * tr_x_xyy_xxyzz[i] * fe_0 + tr_x_xyy_xxyzzz[i] * pa_z[i];

        tr_x_xyyz_xxzzzz[i] = tr_x_xz_xxzzzz[i] * fe_0 + tr_x_xyz_xxzzzz[i] * pa_y[i];

        tr_x_xyyz_xyyyyy[i] = tr_x_xyy_xyyyyy[i] * pa_z[i];

        tr_x_xyyz_xyyyyz[i] = tr_x_xyy_xyyyy[i] * fe_0 + tr_x_xyy_xyyyyz[i] * pa_z[i];

        tr_x_xyyz_xyyyzz[i] = 2.0 * tr_x_xyy_xyyyz[i] * fe_0 + tr_x_xyy_xyyyzz[i] * pa_z[i];

        tr_x_xyyz_xyyzzz[i] = 3.0 * tr_x_xyy_xyyzz[i] * fe_0 + tr_x_xyy_xyyzzz[i] * pa_z[i];

        tr_x_xyyz_xyzzzz[i] = 4.0 * tr_x_xyy_xyzzz[i] * fe_0 + tr_x_xyy_xyzzzz[i] * pa_z[i];

        tr_x_xyyz_xzzzzz[i] = tr_x_xz_xzzzzz[i] * fe_0 + tr_x_xyz_xzzzzz[i] * pa_y[i];

        tr_x_xyyz_yyyyyy[i] = tr_x_xyy_yyyyyy[i] * pa_z[i];

        tr_x_xyyz_yyyyyz[i] = ts_yyz_yyyyyz[i] * fe_0 + tr_x_yyz_yyyyyz[i] * pa_x[i];

        tr_x_xyyz_yyyyzz[i] = ts_yyz_yyyyzz[i] * fe_0 + tr_x_yyz_yyyyzz[i] * pa_x[i];

        tr_x_xyyz_yyyzzz[i] = ts_yyz_yyyzzz[i] * fe_0 + tr_x_yyz_yyyzzz[i] * pa_x[i];

        tr_x_xyyz_yyzzzz[i] = ts_yyz_yyzzzz[i] * fe_0 + tr_x_yyz_yyzzzz[i] * pa_x[i];

        tr_x_xyyz_yzzzzz[i] = ts_yyz_yzzzzz[i] * fe_0 + tr_x_yyz_yzzzzz[i] * pa_x[i];

        tr_x_xyyz_zzzzzz[i] = ts_yyz_zzzzzz[i] * fe_0 + tr_x_yyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 224-252 components of targeted buffer : GI

    auto tr_x_xyzz_xxxxxx = pbuffer.data(idx_dip_gi + 224);

    auto tr_x_xyzz_xxxxxy = pbuffer.data(idx_dip_gi + 225);

    auto tr_x_xyzz_xxxxxz = pbuffer.data(idx_dip_gi + 226);

    auto tr_x_xyzz_xxxxyy = pbuffer.data(idx_dip_gi + 227);

    auto tr_x_xyzz_xxxxyz = pbuffer.data(idx_dip_gi + 228);

    auto tr_x_xyzz_xxxxzz = pbuffer.data(idx_dip_gi + 229);

    auto tr_x_xyzz_xxxyyy = pbuffer.data(idx_dip_gi + 230);

    auto tr_x_xyzz_xxxyyz = pbuffer.data(idx_dip_gi + 231);

    auto tr_x_xyzz_xxxyzz = pbuffer.data(idx_dip_gi + 232);

    auto tr_x_xyzz_xxxzzz = pbuffer.data(idx_dip_gi + 233);

    auto tr_x_xyzz_xxyyyy = pbuffer.data(idx_dip_gi + 234);

    auto tr_x_xyzz_xxyyyz = pbuffer.data(idx_dip_gi + 235);

    auto tr_x_xyzz_xxyyzz = pbuffer.data(idx_dip_gi + 236);

    auto tr_x_xyzz_xxyzzz = pbuffer.data(idx_dip_gi + 237);

    auto tr_x_xyzz_xxzzzz = pbuffer.data(idx_dip_gi + 238);

    auto tr_x_xyzz_xyyyyy = pbuffer.data(idx_dip_gi + 239);

    auto tr_x_xyzz_xyyyyz = pbuffer.data(idx_dip_gi + 240);

    auto tr_x_xyzz_xyyyzz = pbuffer.data(idx_dip_gi + 241);

    auto tr_x_xyzz_xyyzzz = pbuffer.data(idx_dip_gi + 242);

    auto tr_x_xyzz_xyzzzz = pbuffer.data(idx_dip_gi + 243);

    auto tr_x_xyzz_xzzzzz = pbuffer.data(idx_dip_gi + 244);

    auto tr_x_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 245);

    auto tr_x_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 246);

    auto tr_x_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 247);

    auto tr_x_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 248);

    auto tr_x_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 249);

    auto tr_x_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 250);

    auto tr_x_xyzz_zzzzzz = pbuffer.data(idx_dip_gi + 251);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyzz_xxxxxx, tr_x_xyzz_xxxxxy, tr_x_xyzz_xxxxxz, tr_x_xyzz_xxxxyy, tr_x_xyzz_xxxxyz, tr_x_xyzz_xxxxzz, tr_x_xyzz_xxxyyy, tr_x_xyzz_xxxyyz, tr_x_xyzz_xxxyzz, tr_x_xyzz_xxxzzz, tr_x_xyzz_xxyyyy, tr_x_xyzz_xxyyyz, tr_x_xyzz_xxyyzz, tr_x_xyzz_xxyzzz, tr_x_xyzz_xxzzzz, tr_x_xyzz_xyyyyy, tr_x_xyzz_xyyyyz, tr_x_xyzz_xyyyzz, tr_x_xyzz_xyyzzz, tr_x_xyzz_xyzzzz, tr_x_xyzz_xzzzzz, tr_x_xyzz_yyyyyy, tr_x_xyzz_yyyyyz, tr_x_xyzz_yyyyzz, tr_x_xyzz_yyyzzz, tr_x_xyzz_yyzzzz, tr_x_xyzz_yzzzzz, tr_x_xyzz_zzzzzz, tr_x_xzz_xxxxx, tr_x_xzz_xxxxxx, tr_x_xzz_xxxxxy, tr_x_xzz_xxxxxz, tr_x_xzz_xxxxy, tr_x_xzz_xxxxyy, tr_x_xzz_xxxxyz, tr_x_xzz_xxxxz, tr_x_xzz_xxxxzz, tr_x_xzz_xxxyy, tr_x_xzz_xxxyyy, tr_x_xzz_xxxyyz, tr_x_xzz_xxxyz, tr_x_xzz_xxxyzz, tr_x_xzz_xxxzz, tr_x_xzz_xxxzzz, tr_x_xzz_xxyyy, tr_x_xzz_xxyyyy, tr_x_xzz_xxyyyz, tr_x_xzz_xxyyz, tr_x_xzz_xxyyzz, tr_x_xzz_xxyzz, tr_x_xzz_xxyzzz, tr_x_xzz_xxzzz, tr_x_xzz_xxzzzz, tr_x_xzz_xyyyy, tr_x_xzz_xyyyyy, tr_x_xzz_xyyyyz, tr_x_xzz_xyyyz, tr_x_xzz_xyyyzz, tr_x_xzz_xyyzz, tr_x_xzz_xyyzzz, tr_x_xzz_xyzzz, tr_x_xzz_xyzzzz, tr_x_xzz_xzzzz, tr_x_xzz_xzzzzz, tr_x_xzz_zzzzzz, tr_x_yzz_yyyyyy, tr_x_yzz_yyyyyz, tr_x_yzz_yyyyzz, tr_x_yzz_yyyzzz, tr_x_yzz_yyzzzz, tr_x_yzz_yzzzzz, ts_yzz_yyyyyy, ts_yzz_yyyyyz, ts_yzz_yyyyzz, ts_yzz_yyyzzz, ts_yzz_yyzzzz, ts_yzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzz_xxxxxx[i] = tr_x_xzz_xxxxxx[i] * pa_y[i];

        tr_x_xyzz_xxxxxy[i] = tr_x_xzz_xxxxx[i] * fe_0 + tr_x_xzz_xxxxxy[i] * pa_y[i];

        tr_x_xyzz_xxxxxz[i] = tr_x_xzz_xxxxxz[i] * pa_y[i];

        tr_x_xyzz_xxxxyy[i] = 2.0 * tr_x_xzz_xxxxy[i] * fe_0 + tr_x_xzz_xxxxyy[i] * pa_y[i];

        tr_x_xyzz_xxxxyz[i] = tr_x_xzz_xxxxz[i] * fe_0 + tr_x_xzz_xxxxyz[i] * pa_y[i];

        tr_x_xyzz_xxxxzz[i] = tr_x_xzz_xxxxzz[i] * pa_y[i];

        tr_x_xyzz_xxxyyy[i] = 3.0 * tr_x_xzz_xxxyy[i] * fe_0 + tr_x_xzz_xxxyyy[i] * pa_y[i];

        tr_x_xyzz_xxxyyz[i] = 2.0 * tr_x_xzz_xxxyz[i] * fe_0 + tr_x_xzz_xxxyyz[i] * pa_y[i];

        tr_x_xyzz_xxxyzz[i] = tr_x_xzz_xxxzz[i] * fe_0 + tr_x_xzz_xxxyzz[i] * pa_y[i];

        tr_x_xyzz_xxxzzz[i] = tr_x_xzz_xxxzzz[i] * pa_y[i];

        tr_x_xyzz_xxyyyy[i] = 4.0 * tr_x_xzz_xxyyy[i] * fe_0 + tr_x_xzz_xxyyyy[i] * pa_y[i];

        tr_x_xyzz_xxyyyz[i] = 3.0 * tr_x_xzz_xxyyz[i] * fe_0 + tr_x_xzz_xxyyyz[i] * pa_y[i];

        tr_x_xyzz_xxyyzz[i] = 2.0 * tr_x_xzz_xxyzz[i] * fe_0 + tr_x_xzz_xxyyzz[i] * pa_y[i];

        tr_x_xyzz_xxyzzz[i] = tr_x_xzz_xxzzz[i] * fe_0 + tr_x_xzz_xxyzzz[i] * pa_y[i];

        tr_x_xyzz_xxzzzz[i] = tr_x_xzz_xxzzzz[i] * pa_y[i];

        tr_x_xyzz_xyyyyy[i] = 5.0 * tr_x_xzz_xyyyy[i] * fe_0 + tr_x_xzz_xyyyyy[i] * pa_y[i];

        tr_x_xyzz_xyyyyz[i] = 4.0 * tr_x_xzz_xyyyz[i] * fe_0 + tr_x_xzz_xyyyyz[i] * pa_y[i];

        tr_x_xyzz_xyyyzz[i] = 3.0 * tr_x_xzz_xyyzz[i] * fe_0 + tr_x_xzz_xyyyzz[i] * pa_y[i];

        tr_x_xyzz_xyyzzz[i] = 2.0 * tr_x_xzz_xyzzz[i] * fe_0 + tr_x_xzz_xyyzzz[i] * pa_y[i];

        tr_x_xyzz_xyzzzz[i] = tr_x_xzz_xzzzz[i] * fe_0 + tr_x_xzz_xyzzzz[i] * pa_y[i];

        tr_x_xyzz_xzzzzz[i] = tr_x_xzz_xzzzzz[i] * pa_y[i];

        tr_x_xyzz_yyyyyy[i] = ts_yzz_yyyyyy[i] * fe_0 + tr_x_yzz_yyyyyy[i] * pa_x[i];

        tr_x_xyzz_yyyyyz[i] = ts_yzz_yyyyyz[i] * fe_0 + tr_x_yzz_yyyyyz[i] * pa_x[i];

        tr_x_xyzz_yyyyzz[i] = ts_yzz_yyyyzz[i] * fe_0 + tr_x_yzz_yyyyzz[i] * pa_x[i];

        tr_x_xyzz_yyyzzz[i] = ts_yzz_yyyzzz[i] * fe_0 + tr_x_yzz_yyyzzz[i] * pa_x[i];

        tr_x_xyzz_yyzzzz[i] = ts_yzz_yyzzzz[i] * fe_0 + tr_x_yzz_yyzzzz[i] * pa_x[i];

        tr_x_xyzz_yzzzzz[i] = ts_yzz_yzzzzz[i] * fe_0 + tr_x_yzz_yzzzzz[i] * pa_x[i];

        tr_x_xyzz_zzzzzz[i] = tr_x_xzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : GI

    auto tr_x_xzzz_xxxxxx = pbuffer.data(idx_dip_gi + 252);

    auto tr_x_xzzz_xxxxxy = pbuffer.data(idx_dip_gi + 253);

    auto tr_x_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 254);

    auto tr_x_xzzz_xxxxyy = pbuffer.data(idx_dip_gi + 255);

    auto tr_x_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 256);

    auto tr_x_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 257);

    auto tr_x_xzzz_xxxyyy = pbuffer.data(idx_dip_gi + 258);

    auto tr_x_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 259);

    auto tr_x_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 260);

    auto tr_x_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 261);

    auto tr_x_xzzz_xxyyyy = pbuffer.data(idx_dip_gi + 262);

    auto tr_x_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 263);

    auto tr_x_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 264);

    auto tr_x_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 265);

    auto tr_x_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 266);

    auto tr_x_xzzz_xyyyyy = pbuffer.data(idx_dip_gi + 267);

    auto tr_x_xzzz_xyyyyz = pbuffer.data(idx_dip_gi + 268);

    auto tr_x_xzzz_xyyyzz = pbuffer.data(idx_dip_gi + 269);

    auto tr_x_xzzz_xyyzzz = pbuffer.data(idx_dip_gi + 270);

    auto tr_x_xzzz_xyzzzz = pbuffer.data(idx_dip_gi + 271);

    auto tr_x_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 272);

    auto tr_x_xzzz_yyyyyy = pbuffer.data(idx_dip_gi + 273);

    auto tr_x_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 274);

    auto tr_x_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 275);

    auto tr_x_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 276);

    auto tr_x_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 277);

    auto tr_x_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 278);

    auto tr_x_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 279);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xz_xxxxxx, tr_x_xz_xxxxxy, tr_x_xz_xxxxyy, tr_x_xz_xxxyyy, tr_x_xz_xxyyyy, tr_x_xz_xyyyyy, tr_x_xzz_xxxxxx, tr_x_xzz_xxxxxy, tr_x_xzz_xxxxyy, tr_x_xzz_xxxyyy, tr_x_xzz_xxyyyy, tr_x_xzz_xyyyyy, tr_x_xzzz_xxxxxx, tr_x_xzzz_xxxxxy, tr_x_xzzz_xxxxxz, tr_x_xzzz_xxxxyy, tr_x_xzzz_xxxxyz, tr_x_xzzz_xxxxzz, tr_x_xzzz_xxxyyy, tr_x_xzzz_xxxyyz, tr_x_xzzz_xxxyzz, tr_x_xzzz_xxxzzz, tr_x_xzzz_xxyyyy, tr_x_xzzz_xxyyyz, tr_x_xzzz_xxyyzz, tr_x_xzzz_xxyzzz, tr_x_xzzz_xxzzzz, tr_x_xzzz_xyyyyy, tr_x_xzzz_xyyyyz, tr_x_xzzz_xyyyzz, tr_x_xzzz_xyyzzz, tr_x_xzzz_xyzzzz, tr_x_xzzz_xzzzzz, tr_x_xzzz_yyyyyy, tr_x_xzzz_yyyyyz, tr_x_xzzz_yyyyzz, tr_x_xzzz_yyyzzz, tr_x_xzzz_yyzzzz, tr_x_xzzz_yzzzzz, tr_x_xzzz_zzzzzz, tr_x_zzz_xxxxxz, tr_x_zzz_xxxxyz, tr_x_zzz_xxxxz, tr_x_zzz_xxxxzz, tr_x_zzz_xxxyyz, tr_x_zzz_xxxyz, tr_x_zzz_xxxyzz, tr_x_zzz_xxxzz, tr_x_zzz_xxxzzz, tr_x_zzz_xxyyyz, tr_x_zzz_xxyyz, tr_x_zzz_xxyyzz, tr_x_zzz_xxyzz, tr_x_zzz_xxyzzz, tr_x_zzz_xxzzz, tr_x_zzz_xxzzzz, tr_x_zzz_xyyyyz, tr_x_zzz_xyyyz, tr_x_zzz_xyyyzz, tr_x_zzz_xyyzz, tr_x_zzz_xyyzzz, tr_x_zzz_xyzzz, tr_x_zzz_xyzzzz, tr_x_zzz_xzzzz, tr_x_zzz_xzzzzz, tr_x_zzz_yyyyyy, tr_x_zzz_yyyyyz, tr_x_zzz_yyyyz, tr_x_zzz_yyyyzz, tr_x_zzz_yyyzz, tr_x_zzz_yyyzzz, tr_x_zzz_yyzzz, tr_x_zzz_yyzzzz, tr_x_zzz_yzzzz, tr_x_zzz_yzzzzz, tr_x_zzz_zzzzz, tr_x_zzz_zzzzzz, ts_zzz_xxxxxz, ts_zzz_xxxxyz, ts_zzz_xxxxzz, ts_zzz_xxxyyz, ts_zzz_xxxyzz, ts_zzz_xxxzzz, ts_zzz_xxyyyz, ts_zzz_xxyyzz, ts_zzz_xxyzzz, ts_zzz_xxzzzz, ts_zzz_xyyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzz_xxxxxx[i] = 2.0 * tr_x_xz_xxxxxx[i] * fe_0 + tr_x_xzz_xxxxxx[i] * pa_z[i];

        tr_x_xzzz_xxxxxy[i] = 2.0 * tr_x_xz_xxxxxy[i] * fe_0 + tr_x_xzz_xxxxxy[i] * pa_z[i];

        tr_x_xzzz_xxxxxz[i] = 5.0 * tr_x_zzz_xxxxz[i] * fe_0 + ts_zzz_xxxxxz[i] * fe_0 + tr_x_zzz_xxxxxz[i] * pa_x[i];

        tr_x_xzzz_xxxxyy[i] = 2.0 * tr_x_xz_xxxxyy[i] * fe_0 + tr_x_xzz_xxxxyy[i] * pa_z[i];

        tr_x_xzzz_xxxxyz[i] = 4.0 * tr_x_zzz_xxxyz[i] * fe_0 + ts_zzz_xxxxyz[i] * fe_0 + tr_x_zzz_xxxxyz[i] * pa_x[i];

        tr_x_xzzz_xxxxzz[i] = 4.0 * tr_x_zzz_xxxzz[i] * fe_0 + ts_zzz_xxxxzz[i] * fe_0 + tr_x_zzz_xxxxzz[i] * pa_x[i];

        tr_x_xzzz_xxxyyy[i] = 2.0 * tr_x_xz_xxxyyy[i] * fe_0 + tr_x_xzz_xxxyyy[i] * pa_z[i];

        tr_x_xzzz_xxxyyz[i] = 3.0 * tr_x_zzz_xxyyz[i] * fe_0 + ts_zzz_xxxyyz[i] * fe_0 + tr_x_zzz_xxxyyz[i] * pa_x[i];

        tr_x_xzzz_xxxyzz[i] = 3.0 * tr_x_zzz_xxyzz[i] * fe_0 + ts_zzz_xxxyzz[i] * fe_0 + tr_x_zzz_xxxyzz[i] * pa_x[i];

        tr_x_xzzz_xxxzzz[i] = 3.0 * tr_x_zzz_xxzzz[i] * fe_0 + ts_zzz_xxxzzz[i] * fe_0 + tr_x_zzz_xxxzzz[i] * pa_x[i];

        tr_x_xzzz_xxyyyy[i] = 2.0 * tr_x_xz_xxyyyy[i] * fe_0 + tr_x_xzz_xxyyyy[i] * pa_z[i];

        tr_x_xzzz_xxyyyz[i] = 2.0 * tr_x_zzz_xyyyz[i] * fe_0 + ts_zzz_xxyyyz[i] * fe_0 + tr_x_zzz_xxyyyz[i] * pa_x[i];

        tr_x_xzzz_xxyyzz[i] = 2.0 * tr_x_zzz_xyyzz[i] * fe_0 + ts_zzz_xxyyzz[i] * fe_0 + tr_x_zzz_xxyyzz[i] * pa_x[i];

        tr_x_xzzz_xxyzzz[i] = 2.0 * tr_x_zzz_xyzzz[i] * fe_0 + ts_zzz_xxyzzz[i] * fe_0 + tr_x_zzz_xxyzzz[i] * pa_x[i];

        tr_x_xzzz_xxzzzz[i] = 2.0 * tr_x_zzz_xzzzz[i] * fe_0 + ts_zzz_xxzzzz[i] * fe_0 + tr_x_zzz_xxzzzz[i] * pa_x[i];

        tr_x_xzzz_xyyyyy[i] = 2.0 * tr_x_xz_xyyyyy[i] * fe_0 + tr_x_xzz_xyyyyy[i] * pa_z[i];

        tr_x_xzzz_xyyyyz[i] = tr_x_zzz_yyyyz[i] * fe_0 + ts_zzz_xyyyyz[i] * fe_0 + tr_x_zzz_xyyyyz[i] * pa_x[i];

        tr_x_xzzz_xyyyzz[i] = tr_x_zzz_yyyzz[i] * fe_0 + ts_zzz_xyyyzz[i] * fe_0 + tr_x_zzz_xyyyzz[i] * pa_x[i];

        tr_x_xzzz_xyyzzz[i] = tr_x_zzz_yyzzz[i] * fe_0 + ts_zzz_xyyzzz[i] * fe_0 + tr_x_zzz_xyyzzz[i] * pa_x[i];

        tr_x_xzzz_xyzzzz[i] = tr_x_zzz_yzzzz[i] * fe_0 + ts_zzz_xyzzzz[i] * fe_0 + tr_x_zzz_xyzzzz[i] * pa_x[i];

        tr_x_xzzz_xzzzzz[i] = tr_x_zzz_zzzzz[i] * fe_0 + ts_zzz_xzzzzz[i] * fe_0 + tr_x_zzz_xzzzzz[i] * pa_x[i];

        tr_x_xzzz_yyyyyy[i] = ts_zzz_yyyyyy[i] * fe_0 + tr_x_zzz_yyyyyy[i] * pa_x[i];

        tr_x_xzzz_yyyyyz[i] = ts_zzz_yyyyyz[i] * fe_0 + tr_x_zzz_yyyyyz[i] * pa_x[i];

        tr_x_xzzz_yyyyzz[i] = ts_zzz_yyyyzz[i] * fe_0 + tr_x_zzz_yyyyzz[i] * pa_x[i];

        tr_x_xzzz_yyyzzz[i] = ts_zzz_yyyzzz[i] * fe_0 + tr_x_zzz_yyyzzz[i] * pa_x[i];

        tr_x_xzzz_yyzzzz[i] = ts_zzz_yyzzzz[i] * fe_0 + tr_x_zzz_yyzzzz[i] * pa_x[i];

        tr_x_xzzz_yzzzzz[i] = ts_zzz_yzzzzz[i] * fe_0 + tr_x_zzz_yzzzzz[i] * pa_x[i];

        tr_x_xzzz_zzzzzz[i] = ts_zzz_zzzzzz[i] * fe_0 + tr_x_zzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, tr_x_yy_xxxxxx, tr_x_yy_xxxxxy, tr_x_yy_xxxxxz, tr_x_yy_xxxxyy, tr_x_yy_xxxxyz, tr_x_yy_xxxxzz, tr_x_yy_xxxyyy, tr_x_yy_xxxyyz, tr_x_yy_xxxyzz, tr_x_yy_xxxzzz, tr_x_yy_xxyyyy, tr_x_yy_xxyyyz, tr_x_yy_xxyyzz, tr_x_yy_xxyzzz, tr_x_yy_xxzzzz, tr_x_yy_xyyyyy, tr_x_yy_xyyyyz, tr_x_yy_xyyyzz, tr_x_yy_xyyzzz, tr_x_yy_xyzzzz, tr_x_yy_xzzzzz, tr_x_yy_yyyyyy, tr_x_yy_yyyyyz, tr_x_yy_yyyyzz, tr_x_yy_yyyzzz, tr_x_yy_yyzzzz, tr_x_yy_yzzzzz, tr_x_yy_zzzzzz, tr_x_yyy_xxxxx, tr_x_yyy_xxxxxx, tr_x_yyy_xxxxxy, tr_x_yyy_xxxxxz, tr_x_yyy_xxxxy, tr_x_yyy_xxxxyy, tr_x_yyy_xxxxyz, tr_x_yyy_xxxxz, tr_x_yyy_xxxxzz, tr_x_yyy_xxxyy, tr_x_yyy_xxxyyy, tr_x_yyy_xxxyyz, tr_x_yyy_xxxyz, tr_x_yyy_xxxyzz, tr_x_yyy_xxxzz, tr_x_yyy_xxxzzz, tr_x_yyy_xxyyy, tr_x_yyy_xxyyyy, tr_x_yyy_xxyyyz, tr_x_yyy_xxyyz, tr_x_yyy_xxyyzz, tr_x_yyy_xxyzz, tr_x_yyy_xxyzzz, tr_x_yyy_xxzzz, tr_x_yyy_xxzzzz, tr_x_yyy_xyyyy, tr_x_yyy_xyyyyy, tr_x_yyy_xyyyyz, tr_x_yyy_xyyyz, tr_x_yyy_xyyyzz, tr_x_yyy_xyyzz, tr_x_yyy_xyyzzz, tr_x_yyy_xyzzz, tr_x_yyy_xyzzzz, tr_x_yyy_xzzzz, tr_x_yyy_xzzzzz, tr_x_yyy_yyyyy, tr_x_yyy_yyyyyy, tr_x_yyy_yyyyyz, tr_x_yyy_yyyyz, tr_x_yyy_yyyyzz, tr_x_yyy_yyyzz, tr_x_yyy_yyyzzz, tr_x_yyy_yyzzz, tr_x_yyy_yyzzzz, tr_x_yyy_yzzzz, tr_x_yyy_yzzzzz, tr_x_yyy_zzzzz, tr_x_yyy_zzzzzz, tr_x_yyyy_xxxxxx, tr_x_yyyy_xxxxxy, tr_x_yyyy_xxxxxz, tr_x_yyyy_xxxxyy, tr_x_yyyy_xxxxyz, tr_x_yyyy_xxxxzz, tr_x_yyyy_xxxyyy, tr_x_yyyy_xxxyyz, tr_x_yyyy_xxxyzz, tr_x_yyyy_xxxzzz, tr_x_yyyy_xxyyyy, tr_x_yyyy_xxyyyz, tr_x_yyyy_xxyyzz, tr_x_yyyy_xxyzzz, tr_x_yyyy_xxzzzz, tr_x_yyyy_xyyyyy, tr_x_yyyy_xyyyyz, tr_x_yyyy_xyyyzz, tr_x_yyyy_xyyzzz, tr_x_yyyy_xyzzzz, tr_x_yyyy_xzzzzz, tr_x_yyyy_yyyyyy, tr_x_yyyy_yyyyyz, tr_x_yyyy_yyyyzz, tr_x_yyyy_yyyzzz, tr_x_yyyy_yyzzzz, tr_x_yyyy_yzzzzz, tr_x_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyy_xxxxxx[i] = 3.0 * tr_x_yy_xxxxxx[i] * fe_0 + tr_x_yyy_xxxxxx[i] * pa_y[i];

        tr_x_yyyy_xxxxxy[i] = 3.0 * tr_x_yy_xxxxxy[i] * fe_0 + tr_x_yyy_xxxxx[i] * fe_0 + tr_x_yyy_xxxxxy[i] * pa_y[i];

        tr_x_yyyy_xxxxxz[i] = 3.0 * tr_x_yy_xxxxxz[i] * fe_0 + tr_x_yyy_xxxxxz[i] * pa_y[i];

        tr_x_yyyy_xxxxyy[i] = 3.0 * tr_x_yy_xxxxyy[i] * fe_0 + 2.0 * tr_x_yyy_xxxxy[i] * fe_0 + tr_x_yyy_xxxxyy[i] * pa_y[i];

        tr_x_yyyy_xxxxyz[i] = 3.0 * tr_x_yy_xxxxyz[i] * fe_0 + tr_x_yyy_xxxxz[i] * fe_0 + tr_x_yyy_xxxxyz[i] * pa_y[i];

        tr_x_yyyy_xxxxzz[i] = 3.0 * tr_x_yy_xxxxzz[i] * fe_0 + tr_x_yyy_xxxxzz[i] * pa_y[i];

        tr_x_yyyy_xxxyyy[i] = 3.0 * tr_x_yy_xxxyyy[i] * fe_0 + 3.0 * tr_x_yyy_xxxyy[i] * fe_0 + tr_x_yyy_xxxyyy[i] * pa_y[i];

        tr_x_yyyy_xxxyyz[i] = 3.0 * tr_x_yy_xxxyyz[i] * fe_0 + 2.0 * tr_x_yyy_xxxyz[i] * fe_0 + tr_x_yyy_xxxyyz[i] * pa_y[i];

        tr_x_yyyy_xxxyzz[i] = 3.0 * tr_x_yy_xxxyzz[i] * fe_0 + tr_x_yyy_xxxzz[i] * fe_0 + tr_x_yyy_xxxyzz[i] * pa_y[i];

        tr_x_yyyy_xxxzzz[i] = 3.0 * tr_x_yy_xxxzzz[i] * fe_0 + tr_x_yyy_xxxzzz[i] * pa_y[i];

        tr_x_yyyy_xxyyyy[i] = 3.0 * tr_x_yy_xxyyyy[i] * fe_0 + 4.0 * tr_x_yyy_xxyyy[i] * fe_0 + tr_x_yyy_xxyyyy[i] * pa_y[i];

        tr_x_yyyy_xxyyyz[i] = 3.0 * tr_x_yy_xxyyyz[i] * fe_0 + 3.0 * tr_x_yyy_xxyyz[i] * fe_0 + tr_x_yyy_xxyyyz[i] * pa_y[i];

        tr_x_yyyy_xxyyzz[i] = 3.0 * tr_x_yy_xxyyzz[i] * fe_0 + 2.0 * tr_x_yyy_xxyzz[i] * fe_0 + tr_x_yyy_xxyyzz[i] * pa_y[i];

        tr_x_yyyy_xxyzzz[i] = 3.0 * tr_x_yy_xxyzzz[i] * fe_0 + tr_x_yyy_xxzzz[i] * fe_0 + tr_x_yyy_xxyzzz[i] * pa_y[i];

        tr_x_yyyy_xxzzzz[i] = 3.0 * tr_x_yy_xxzzzz[i] * fe_0 + tr_x_yyy_xxzzzz[i] * pa_y[i];

        tr_x_yyyy_xyyyyy[i] = 3.0 * tr_x_yy_xyyyyy[i] * fe_0 + 5.0 * tr_x_yyy_xyyyy[i] * fe_0 + tr_x_yyy_xyyyyy[i] * pa_y[i];

        tr_x_yyyy_xyyyyz[i] = 3.0 * tr_x_yy_xyyyyz[i] * fe_0 + 4.0 * tr_x_yyy_xyyyz[i] * fe_0 + tr_x_yyy_xyyyyz[i] * pa_y[i];

        tr_x_yyyy_xyyyzz[i] = 3.0 * tr_x_yy_xyyyzz[i] * fe_0 + 3.0 * tr_x_yyy_xyyzz[i] * fe_0 + tr_x_yyy_xyyyzz[i] * pa_y[i];

        tr_x_yyyy_xyyzzz[i] = 3.0 * tr_x_yy_xyyzzz[i] * fe_0 + 2.0 * tr_x_yyy_xyzzz[i] * fe_0 + tr_x_yyy_xyyzzz[i] * pa_y[i];

        tr_x_yyyy_xyzzzz[i] = 3.0 * tr_x_yy_xyzzzz[i] * fe_0 + tr_x_yyy_xzzzz[i] * fe_0 + tr_x_yyy_xyzzzz[i] * pa_y[i];

        tr_x_yyyy_xzzzzz[i] = 3.0 * tr_x_yy_xzzzzz[i] * fe_0 + tr_x_yyy_xzzzzz[i] * pa_y[i];

        tr_x_yyyy_yyyyyy[i] = 3.0 * tr_x_yy_yyyyyy[i] * fe_0 + 6.0 * tr_x_yyy_yyyyy[i] * fe_0 + tr_x_yyy_yyyyyy[i] * pa_y[i];

        tr_x_yyyy_yyyyyz[i] = 3.0 * tr_x_yy_yyyyyz[i] * fe_0 + 5.0 * tr_x_yyy_yyyyz[i] * fe_0 + tr_x_yyy_yyyyyz[i] * pa_y[i];

        tr_x_yyyy_yyyyzz[i] = 3.0 * tr_x_yy_yyyyzz[i] * fe_0 + 4.0 * tr_x_yyy_yyyzz[i] * fe_0 + tr_x_yyy_yyyyzz[i] * pa_y[i];

        tr_x_yyyy_yyyzzz[i] = 3.0 * tr_x_yy_yyyzzz[i] * fe_0 + 3.0 * tr_x_yyy_yyzzz[i] * fe_0 + tr_x_yyy_yyyzzz[i] * pa_y[i];

        tr_x_yyyy_yyzzzz[i] = 3.0 * tr_x_yy_yyzzzz[i] * fe_0 + 2.0 * tr_x_yyy_yzzzz[i] * fe_0 + tr_x_yyy_yyzzzz[i] * pa_y[i];

        tr_x_yyyy_yzzzzz[i] = 3.0 * tr_x_yy_yzzzzz[i] * fe_0 + tr_x_yyy_zzzzz[i] * fe_0 + tr_x_yyy_yzzzzz[i] * pa_y[i];

        tr_x_yyyy_zzzzzz[i] = 3.0 * tr_x_yy_zzzzzz[i] * fe_0 + tr_x_yyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 308-336 components of targeted buffer : GI

    auto tr_x_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 308);

    auto tr_x_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 309);

    auto tr_x_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 310);

    auto tr_x_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 311);

    auto tr_x_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 312);

    auto tr_x_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 313);

    auto tr_x_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 314);

    auto tr_x_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 315);

    auto tr_x_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 316);

    auto tr_x_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 317);

    auto tr_x_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 318);

    auto tr_x_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 319);

    auto tr_x_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 320);

    auto tr_x_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 321);

    auto tr_x_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 322);

    auto tr_x_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 323);

    auto tr_x_yyyz_xyyyyz = pbuffer.data(idx_dip_gi + 324);

    auto tr_x_yyyz_xyyyzz = pbuffer.data(idx_dip_gi + 325);

    auto tr_x_yyyz_xyyzzz = pbuffer.data(idx_dip_gi + 326);

    auto tr_x_yyyz_xyzzzz = pbuffer.data(idx_dip_gi + 327);

    auto tr_x_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 328);

    auto tr_x_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 329);

    auto tr_x_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 330);

    auto tr_x_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 331);

    auto tr_x_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 332);

    auto tr_x_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 333);

    auto tr_x_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 334);

    auto tr_x_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 335);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyy_xxxxxx, tr_x_yyy_xxxxxy, tr_x_yyy_xxxxy, tr_x_yyy_xxxxyy, tr_x_yyy_xxxxyz, tr_x_yyy_xxxyy, tr_x_yyy_xxxyyy, tr_x_yyy_xxxyyz, tr_x_yyy_xxxyz, tr_x_yyy_xxxyzz, tr_x_yyy_xxyyy, tr_x_yyy_xxyyyy, tr_x_yyy_xxyyyz, tr_x_yyy_xxyyz, tr_x_yyy_xxyyzz, tr_x_yyy_xxyzz, tr_x_yyy_xxyzzz, tr_x_yyy_xyyyy, tr_x_yyy_xyyyyy, tr_x_yyy_xyyyyz, tr_x_yyy_xyyyz, tr_x_yyy_xyyyzz, tr_x_yyy_xyyzz, tr_x_yyy_xyyzzz, tr_x_yyy_xyzzz, tr_x_yyy_xyzzzz, tr_x_yyy_yyyyy, tr_x_yyy_yyyyyy, tr_x_yyy_yyyyyz, tr_x_yyy_yyyyz, tr_x_yyy_yyyyzz, tr_x_yyy_yyyzz, tr_x_yyy_yyyzzz, tr_x_yyy_yyzzz, tr_x_yyy_yyzzzz, tr_x_yyy_yzzzz, tr_x_yyy_yzzzzz, tr_x_yyyz_xxxxxx, tr_x_yyyz_xxxxxy, tr_x_yyyz_xxxxxz, tr_x_yyyz_xxxxyy, tr_x_yyyz_xxxxyz, tr_x_yyyz_xxxxzz, tr_x_yyyz_xxxyyy, tr_x_yyyz_xxxyyz, tr_x_yyyz_xxxyzz, tr_x_yyyz_xxxzzz, tr_x_yyyz_xxyyyy, tr_x_yyyz_xxyyyz, tr_x_yyyz_xxyyzz, tr_x_yyyz_xxyzzz, tr_x_yyyz_xxzzzz, tr_x_yyyz_xyyyyy, tr_x_yyyz_xyyyyz, tr_x_yyyz_xyyyzz, tr_x_yyyz_xyyzzz, tr_x_yyyz_xyzzzz, tr_x_yyyz_xzzzzz, tr_x_yyyz_yyyyyy, tr_x_yyyz_yyyyyz, tr_x_yyyz_yyyyzz, tr_x_yyyz_yyyzzz, tr_x_yyyz_yyzzzz, tr_x_yyyz_yzzzzz, tr_x_yyyz_zzzzzz, tr_x_yyz_xxxxxz, tr_x_yyz_xxxxzz, tr_x_yyz_xxxzzz, tr_x_yyz_xxzzzz, tr_x_yyz_xzzzzz, tr_x_yyz_zzzzzz, tr_x_yz_xxxxxz, tr_x_yz_xxxxzz, tr_x_yz_xxxzzz, tr_x_yz_xxzzzz, tr_x_yz_xzzzzz, tr_x_yz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyz_xxxxxx[i] = tr_x_yyy_xxxxxx[i] * pa_z[i];

        tr_x_yyyz_xxxxxy[i] = tr_x_yyy_xxxxxy[i] * pa_z[i];

        tr_x_yyyz_xxxxxz[i] = 2.0 * tr_x_yz_xxxxxz[i] * fe_0 + tr_x_yyz_xxxxxz[i] * pa_y[i];

        tr_x_yyyz_xxxxyy[i] = tr_x_yyy_xxxxyy[i] * pa_z[i];

        tr_x_yyyz_xxxxyz[i] = tr_x_yyy_xxxxy[i] * fe_0 + tr_x_yyy_xxxxyz[i] * pa_z[i];

        tr_x_yyyz_xxxxzz[i] = 2.0 * tr_x_yz_xxxxzz[i] * fe_0 + tr_x_yyz_xxxxzz[i] * pa_y[i];

        tr_x_yyyz_xxxyyy[i] = tr_x_yyy_xxxyyy[i] * pa_z[i];

        tr_x_yyyz_xxxyyz[i] = tr_x_yyy_xxxyy[i] * fe_0 + tr_x_yyy_xxxyyz[i] * pa_z[i];

        tr_x_yyyz_xxxyzz[i] = 2.0 * tr_x_yyy_xxxyz[i] * fe_0 + tr_x_yyy_xxxyzz[i] * pa_z[i];

        tr_x_yyyz_xxxzzz[i] = 2.0 * tr_x_yz_xxxzzz[i] * fe_0 + tr_x_yyz_xxxzzz[i] * pa_y[i];

        tr_x_yyyz_xxyyyy[i] = tr_x_yyy_xxyyyy[i] * pa_z[i];

        tr_x_yyyz_xxyyyz[i] = tr_x_yyy_xxyyy[i] * fe_0 + tr_x_yyy_xxyyyz[i] * pa_z[i];

        tr_x_yyyz_xxyyzz[i] = 2.0 * tr_x_yyy_xxyyz[i] * fe_0 + tr_x_yyy_xxyyzz[i] * pa_z[i];

        tr_x_yyyz_xxyzzz[i] = 3.0 * tr_x_yyy_xxyzz[i] * fe_0 + tr_x_yyy_xxyzzz[i] * pa_z[i];

        tr_x_yyyz_xxzzzz[i] = 2.0 * tr_x_yz_xxzzzz[i] * fe_0 + tr_x_yyz_xxzzzz[i] * pa_y[i];

        tr_x_yyyz_xyyyyy[i] = tr_x_yyy_xyyyyy[i] * pa_z[i];

        tr_x_yyyz_xyyyyz[i] = tr_x_yyy_xyyyy[i] * fe_0 + tr_x_yyy_xyyyyz[i] * pa_z[i];

        tr_x_yyyz_xyyyzz[i] = 2.0 * tr_x_yyy_xyyyz[i] * fe_0 + tr_x_yyy_xyyyzz[i] * pa_z[i];

        tr_x_yyyz_xyyzzz[i] = 3.0 * tr_x_yyy_xyyzz[i] * fe_0 + tr_x_yyy_xyyzzz[i] * pa_z[i];

        tr_x_yyyz_xyzzzz[i] = 4.0 * tr_x_yyy_xyzzz[i] * fe_0 + tr_x_yyy_xyzzzz[i] * pa_z[i];

        tr_x_yyyz_xzzzzz[i] = 2.0 * tr_x_yz_xzzzzz[i] * fe_0 + tr_x_yyz_xzzzzz[i] * pa_y[i];

        tr_x_yyyz_yyyyyy[i] = tr_x_yyy_yyyyyy[i] * pa_z[i];

        tr_x_yyyz_yyyyyz[i] = tr_x_yyy_yyyyy[i] * fe_0 + tr_x_yyy_yyyyyz[i] * pa_z[i];

        tr_x_yyyz_yyyyzz[i] = 2.0 * tr_x_yyy_yyyyz[i] * fe_0 + tr_x_yyy_yyyyzz[i] * pa_z[i];

        tr_x_yyyz_yyyzzz[i] = 3.0 * tr_x_yyy_yyyzz[i] * fe_0 + tr_x_yyy_yyyzzz[i] * pa_z[i];

        tr_x_yyyz_yyzzzz[i] = 4.0 * tr_x_yyy_yyzzz[i] * fe_0 + tr_x_yyy_yyzzzz[i] * pa_z[i];

        tr_x_yyyz_yzzzzz[i] = 5.0 * tr_x_yyy_yzzzz[i] * fe_0 + tr_x_yyy_yzzzzz[i] * pa_z[i];

        tr_x_yyyz_zzzzzz[i] = 2.0 * tr_x_yz_zzzzzz[i] * fe_0 + tr_x_yyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 336-364 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yy_xxxxxy, tr_x_yy_xxxxyy, tr_x_yy_xxxyyy, tr_x_yy_xxyyyy, tr_x_yy_xyyyyy, tr_x_yy_yyyyyy, tr_x_yyz_xxxxxy, tr_x_yyz_xxxxyy, tr_x_yyz_xxxyyy, tr_x_yyz_xxyyyy, tr_x_yyz_xyyyyy, tr_x_yyz_yyyyyy, tr_x_yyzz_xxxxxx, tr_x_yyzz_xxxxxy, tr_x_yyzz_xxxxxz, tr_x_yyzz_xxxxyy, tr_x_yyzz_xxxxyz, tr_x_yyzz_xxxxzz, tr_x_yyzz_xxxyyy, tr_x_yyzz_xxxyyz, tr_x_yyzz_xxxyzz, tr_x_yyzz_xxxzzz, tr_x_yyzz_xxyyyy, tr_x_yyzz_xxyyyz, tr_x_yyzz_xxyyzz, tr_x_yyzz_xxyzzz, tr_x_yyzz_xxzzzz, tr_x_yyzz_xyyyyy, tr_x_yyzz_xyyyyz, tr_x_yyzz_xyyyzz, tr_x_yyzz_xyyzzz, tr_x_yyzz_xyzzzz, tr_x_yyzz_xzzzzz, tr_x_yyzz_yyyyyy, tr_x_yyzz_yyyyyz, tr_x_yyzz_yyyyzz, tr_x_yyzz_yyyzzz, tr_x_yyzz_yyzzzz, tr_x_yyzz_yzzzzz, tr_x_yyzz_zzzzzz, tr_x_yzz_xxxxxx, tr_x_yzz_xxxxxz, tr_x_yzz_xxxxyz, tr_x_yzz_xxxxz, tr_x_yzz_xxxxzz, tr_x_yzz_xxxyyz, tr_x_yzz_xxxyz, tr_x_yzz_xxxyzz, tr_x_yzz_xxxzz, tr_x_yzz_xxxzzz, tr_x_yzz_xxyyyz, tr_x_yzz_xxyyz, tr_x_yzz_xxyyzz, tr_x_yzz_xxyzz, tr_x_yzz_xxyzzz, tr_x_yzz_xxzzz, tr_x_yzz_xxzzzz, tr_x_yzz_xyyyyz, tr_x_yzz_xyyyz, tr_x_yzz_xyyyzz, tr_x_yzz_xyyzz, tr_x_yzz_xyyzzz, tr_x_yzz_xyzzz, tr_x_yzz_xyzzzz, tr_x_yzz_xzzzz, tr_x_yzz_xzzzzz, tr_x_yzz_yyyyyz, tr_x_yzz_yyyyz, tr_x_yzz_yyyyzz, tr_x_yzz_yyyzz, tr_x_yzz_yyyzzz, tr_x_yzz_yyzzz, tr_x_yzz_yyzzzz, tr_x_yzz_yzzzz, tr_x_yzz_yzzzzz, tr_x_yzz_zzzzz, tr_x_yzz_zzzzzz, tr_x_zz_xxxxxx, tr_x_zz_xxxxxz, tr_x_zz_xxxxyz, tr_x_zz_xxxxzz, tr_x_zz_xxxyyz, tr_x_zz_xxxyzz, tr_x_zz_xxxzzz, tr_x_zz_xxyyyz, tr_x_zz_xxyyzz, tr_x_zz_xxyzzz, tr_x_zz_xxzzzz, tr_x_zz_xyyyyz, tr_x_zz_xyyyzz, tr_x_zz_xyyzzz, tr_x_zz_xyzzzz, tr_x_zz_xzzzzz, tr_x_zz_yyyyyz, tr_x_zz_yyyyzz, tr_x_zz_yyyzzz, tr_x_zz_yyzzzz, tr_x_zz_yzzzzz, tr_x_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzz_xxxxxx[i] = tr_x_zz_xxxxxx[i] * fe_0 + tr_x_yzz_xxxxxx[i] * pa_y[i];

        tr_x_yyzz_xxxxxy[i] = tr_x_yy_xxxxxy[i] * fe_0 + tr_x_yyz_xxxxxy[i] * pa_z[i];

        tr_x_yyzz_xxxxxz[i] = tr_x_zz_xxxxxz[i] * fe_0 + tr_x_yzz_xxxxxz[i] * pa_y[i];

        tr_x_yyzz_xxxxyy[i] = tr_x_yy_xxxxyy[i] * fe_0 + tr_x_yyz_xxxxyy[i] * pa_z[i];

        tr_x_yyzz_xxxxyz[i] = tr_x_zz_xxxxyz[i] * fe_0 + tr_x_yzz_xxxxz[i] * fe_0 + tr_x_yzz_xxxxyz[i] * pa_y[i];

        tr_x_yyzz_xxxxzz[i] = tr_x_zz_xxxxzz[i] * fe_0 + tr_x_yzz_xxxxzz[i] * pa_y[i];

        tr_x_yyzz_xxxyyy[i] = tr_x_yy_xxxyyy[i] * fe_0 + tr_x_yyz_xxxyyy[i] * pa_z[i];

        tr_x_yyzz_xxxyyz[i] = tr_x_zz_xxxyyz[i] * fe_0 + 2.0 * tr_x_yzz_xxxyz[i] * fe_0 + tr_x_yzz_xxxyyz[i] * pa_y[i];

        tr_x_yyzz_xxxyzz[i] = tr_x_zz_xxxyzz[i] * fe_0 + tr_x_yzz_xxxzz[i] * fe_0 + tr_x_yzz_xxxyzz[i] * pa_y[i];

        tr_x_yyzz_xxxzzz[i] = tr_x_zz_xxxzzz[i] * fe_0 + tr_x_yzz_xxxzzz[i] * pa_y[i];

        tr_x_yyzz_xxyyyy[i] = tr_x_yy_xxyyyy[i] * fe_0 + tr_x_yyz_xxyyyy[i] * pa_z[i];

        tr_x_yyzz_xxyyyz[i] = tr_x_zz_xxyyyz[i] * fe_0 + 3.0 * tr_x_yzz_xxyyz[i] * fe_0 + tr_x_yzz_xxyyyz[i] * pa_y[i];

        tr_x_yyzz_xxyyzz[i] = tr_x_zz_xxyyzz[i] * fe_0 + 2.0 * tr_x_yzz_xxyzz[i] * fe_0 + tr_x_yzz_xxyyzz[i] * pa_y[i];

        tr_x_yyzz_xxyzzz[i] = tr_x_zz_xxyzzz[i] * fe_0 + tr_x_yzz_xxzzz[i] * fe_0 + tr_x_yzz_xxyzzz[i] * pa_y[i];

        tr_x_yyzz_xxzzzz[i] = tr_x_zz_xxzzzz[i] * fe_0 + tr_x_yzz_xxzzzz[i] * pa_y[i];

        tr_x_yyzz_xyyyyy[i] = tr_x_yy_xyyyyy[i] * fe_0 + tr_x_yyz_xyyyyy[i] * pa_z[i];

        tr_x_yyzz_xyyyyz[i] = tr_x_zz_xyyyyz[i] * fe_0 + 4.0 * tr_x_yzz_xyyyz[i] * fe_0 + tr_x_yzz_xyyyyz[i] * pa_y[i];

        tr_x_yyzz_xyyyzz[i] = tr_x_zz_xyyyzz[i] * fe_0 + 3.0 * tr_x_yzz_xyyzz[i] * fe_0 + tr_x_yzz_xyyyzz[i] * pa_y[i];

        tr_x_yyzz_xyyzzz[i] = tr_x_zz_xyyzzz[i] * fe_0 + 2.0 * tr_x_yzz_xyzzz[i] * fe_0 + tr_x_yzz_xyyzzz[i] * pa_y[i];

        tr_x_yyzz_xyzzzz[i] = tr_x_zz_xyzzzz[i] * fe_0 + tr_x_yzz_xzzzz[i] * fe_0 + tr_x_yzz_xyzzzz[i] * pa_y[i];

        tr_x_yyzz_xzzzzz[i] = tr_x_zz_xzzzzz[i] * fe_0 + tr_x_yzz_xzzzzz[i] * pa_y[i];

        tr_x_yyzz_yyyyyy[i] = tr_x_yy_yyyyyy[i] * fe_0 + tr_x_yyz_yyyyyy[i] * pa_z[i];

        tr_x_yyzz_yyyyyz[i] = tr_x_zz_yyyyyz[i] * fe_0 + 5.0 * tr_x_yzz_yyyyz[i] * fe_0 + tr_x_yzz_yyyyyz[i] * pa_y[i];

        tr_x_yyzz_yyyyzz[i] = tr_x_zz_yyyyzz[i] * fe_0 + 4.0 * tr_x_yzz_yyyzz[i] * fe_0 + tr_x_yzz_yyyyzz[i] * pa_y[i];

        tr_x_yyzz_yyyzzz[i] = tr_x_zz_yyyzzz[i] * fe_0 + 3.0 * tr_x_yzz_yyzzz[i] * fe_0 + tr_x_yzz_yyyzzz[i] * pa_y[i];

        tr_x_yyzz_yyzzzz[i] = tr_x_zz_yyzzzz[i] * fe_0 + 2.0 * tr_x_yzz_yzzzz[i] * fe_0 + tr_x_yzz_yyzzzz[i] * pa_y[i];

        tr_x_yyzz_yzzzzz[i] = tr_x_zz_yzzzzz[i] * fe_0 + tr_x_yzz_zzzzz[i] * fe_0 + tr_x_yzz_yzzzzz[i] * pa_y[i];

        tr_x_yyzz_zzzzzz[i] = tr_x_zz_zzzzzz[i] * fe_0 + tr_x_yzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 364-392 components of targeted buffer : GI

    auto tr_x_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 364);

    auto tr_x_yzzz_xxxxxy = pbuffer.data(idx_dip_gi + 365);

    auto tr_x_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 366);

    auto tr_x_yzzz_xxxxyy = pbuffer.data(idx_dip_gi + 367);

    auto tr_x_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 368);

    auto tr_x_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 369);

    auto tr_x_yzzz_xxxyyy = pbuffer.data(idx_dip_gi + 370);

    auto tr_x_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 371);

    auto tr_x_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 372);

    auto tr_x_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 373);

    auto tr_x_yzzz_xxyyyy = pbuffer.data(idx_dip_gi + 374);

    auto tr_x_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 375);

    auto tr_x_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 376);

    auto tr_x_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 377);

    auto tr_x_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 378);

    auto tr_x_yzzz_xyyyyy = pbuffer.data(idx_dip_gi + 379);

    auto tr_x_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 380);

    auto tr_x_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 381);

    auto tr_x_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 382);

    auto tr_x_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 383);

    auto tr_x_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 384);

    auto tr_x_yzzz_yyyyyy = pbuffer.data(idx_dip_gi + 385);

    auto tr_x_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 386);

    auto tr_x_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 387);

    auto tr_x_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 388);

    auto tr_x_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 389);

    auto tr_x_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 390);

    auto tr_x_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 391);

    #pragma omp simd aligned(pa_y, tr_x_yzzz_xxxxxx, tr_x_yzzz_xxxxxy, tr_x_yzzz_xxxxxz, tr_x_yzzz_xxxxyy, tr_x_yzzz_xxxxyz, tr_x_yzzz_xxxxzz, tr_x_yzzz_xxxyyy, tr_x_yzzz_xxxyyz, tr_x_yzzz_xxxyzz, tr_x_yzzz_xxxzzz, tr_x_yzzz_xxyyyy, tr_x_yzzz_xxyyyz, tr_x_yzzz_xxyyzz, tr_x_yzzz_xxyzzz, tr_x_yzzz_xxzzzz, tr_x_yzzz_xyyyyy, tr_x_yzzz_xyyyyz, tr_x_yzzz_xyyyzz, tr_x_yzzz_xyyzzz, tr_x_yzzz_xyzzzz, tr_x_yzzz_xzzzzz, tr_x_yzzz_yyyyyy, tr_x_yzzz_yyyyyz, tr_x_yzzz_yyyyzz, tr_x_yzzz_yyyzzz, tr_x_yzzz_yyzzzz, tr_x_yzzz_yzzzzz, tr_x_yzzz_zzzzzz, tr_x_zzz_xxxxx, tr_x_zzz_xxxxxx, tr_x_zzz_xxxxxy, tr_x_zzz_xxxxxz, tr_x_zzz_xxxxy, tr_x_zzz_xxxxyy, tr_x_zzz_xxxxyz, tr_x_zzz_xxxxz, tr_x_zzz_xxxxzz, tr_x_zzz_xxxyy, tr_x_zzz_xxxyyy, tr_x_zzz_xxxyyz, tr_x_zzz_xxxyz, tr_x_zzz_xxxyzz, tr_x_zzz_xxxzz, tr_x_zzz_xxxzzz, tr_x_zzz_xxyyy, tr_x_zzz_xxyyyy, tr_x_zzz_xxyyyz, tr_x_zzz_xxyyz, tr_x_zzz_xxyyzz, tr_x_zzz_xxyzz, tr_x_zzz_xxyzzz, tr_x_zzz_xxzzz, tr_x_zzz_xxzzzz, tr_x_zzz_xyyyy, tr_x_zzz_xyyyyy, tr_x_zzz_xyyyyz, tr_x_zzz_xyyyz, tr_x_zzz_xyyyzz, tr_x_zzz_xyyzz, tr_x_zzz_xyyzzz, tr_x_zzz_xyzzz, tr_x_zzz_xyzzzz, tr_x_zzz_xzzzz, tr_x_zzz_xzzzzz, tr_x_zzz_yyyyy, tr_x_zzz_yyyyyy, tr_x_zzz_yyyyyz, tr_x_zzz_yyyyz, tr_x_zzz_yyyyzz, tr_x_zzz_yyyzz, tr_x_zzz_yyyzzz, tr_x_zzz_yyzzz, tr_x_zzz_yyzzzz, tr_x_zzz_yzzzz, tr_x_zzz_yzzzzz, tr_x_zzz_zzzzz, tr_x_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzz_xxxxxx[i] = tr_x_zzz_xxxxxx[i] * pa_y[i];

        tr_x_yzzz_xxxxxy[i] = tr_x_zzz_xxxxx[i] * fe_0 + tr_x_zzz_xxxxxy[i] * pa_y[i];

        tr_x_yzzz_xxxxxz[i] = tr_x_zzz_xxxxxz[i] * pa_y[i];

        tr_x_yzzz_xxxxyy[i] = 2.0 * tr_x_zzz_xxxxy[i] * fe_0 + tr_x_zzz_xxxxyy[i] * pa_y[i];

        tr_x_yzzz_xxxxyz[i] = tr_x_zzz_xxxxz[i] * fe_0 + tr_x_zzz_xxxxyz[i] * pa_y[i];

        tr_x_yzzz_xxxxzz[i] = tr_x_zzz_xxxxzz[i] * pa_y[i];

        tr_x_yzzz_xxxyyy[i] = 3.0 * tr_x_zzz_xxxyy[i] * fe_0 + tr_x_zzz_xxxyyy[i] * pa_y[i];

        tr_x_yzzz_xxxyyz[i] = 2.0 * tr_x_zzz_xxxyz[i] * fe_0 + tr_x_zzz_xxxyyz[i] * pa_y[i];

        tr_x_yzzz_xxxyzz[i] = tr_x_zzz_xxxzz[i] * fe_0 + tr_x_zzz_xxxyzz[i] * pa_y[i];

        tr_x_yzzz_xxxzzz[i] = tr_x_zzz_xxxzzz[i] * pa_y[i];

        tr_x_yzzz_xxyyyy[i] = 4.0 * tr_x_zzz_xxyyy[i] * fe_0 + tr_x_zzz_xxyyyy[i] * pa_y[i];

        tr_x_yzzz_xxyyyz[i] = 3.0 * tr_x_zzz_xxyyz[i] * fe_0 + tr_x_zzz_xxyyyz[i] * pa_y[i];

        tr_x_yzzz_xxyyzz[i] = 2.0 * tr_x_zzz_xxyzz[i] * fe_0 + tr_x_zzz_xxyyzz[i] * pa_y[i];

        tr_x_yzzz_xxyzzz[i] = tr_x_zzz_xxzzz[i] * fe_0 + tr_x_zzz_xxyzzz[i] * pa_y[i];

        tr_x_yzzz_xxzzzz[i] = tr_x_zzz_xxzzzz[i] * pa_y[i];

        tr_x_yzzz_xyyyyy[i] = 5.0 * tr_x_zzz_xyyyy[i] * fe_0 + tr_x_zzz_xyyyyy[i] * pa_y[i];

        tr_x_yzzz_xyyyyz[i] = 4.0 * tr_x_zzz_xyyyz[i] * fe_0 + tr_x_zzz_xyyyyz[i] * pa_y[i];

        tr_x_yzzz_xyyyzz[i] = 3.0 * tr_x_zzz_xyyzz[i] * fe_0 + tr_x_zzz_xyyyzz[i] * pa_y[i];

        tr_x_yzzz_xyyzzz[i] = 2.0 * tr_x_zzz_xyzzz[i] * fe_0 + tr_x_zzz_xyyzzz[i] * pa_y[i];

        tr_x_yzzz_xyzzzz[i] = tr_x_zzz_xzzzz[i] * fe_0 + tr_x_zzz_xyzzzz[i] * pa_y[i];

        tr_x_yzzz_xzzzzz[i] = tr_x_zzz_xzzzzz[i] * pa_y[i];

        tr_x_yzzz_yyyyyy[i] = 6.0 * tr_x_zzz_yyyyy[i] * fe_0 + tr_x_zzz_yyyyyy[i] * pa_y[i];

        tr_x_yzzz_yyyyyz[i] = 5.0 * tr_x_zzz_yyyyz[i] * fe_0 + tr_x_zzz_yyyyyz[i] * pa_y[i];

        tr_x_yzzz_yyyyzz[i] = 4.0 * tr_x_zzz_yyyzz[i] * fe_0 + tr_x_zzz_yyyyzz[i] * pa_y[i];

        tr_x_yzzz_yyyzzz[i] = 3.0 * tr_x_zzz_yyzzz[i] * fe_0 + tr_x_zzz_yyyzzz[i] * pa_y[i];

        tr_x_yzzz_yyzzzz[i] = 2.0 * tr_x_zzz_yzzzz[i] * fe_0 + tr_x_zzz_yyzzzz[i] * pa_y[i];

        tr_x_yzzz_yzzzzz[i] = tr_x_zzz_zzzzz[i] * fe_0 + tr_x_zzz_yzzzzz[i] * pa_y[i];

        tr_x_yzzz_zzzzzz[i] = tr_x_zzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_z, tr_x_zz_xxxxxx, tr_x_zz_xxxxxy, tr_x_zz_xxxxxz, tr_x_zz_xxxxyy, tr_x_zz_xxxxyz, tr_x_zz_xxxxzz, tr_x_zz_xxxyyy, tr_x_zz_xxxyyz, tr_x_zz_xxxyzz, tr_x_zz_xxxzzz, tr_x_zz_xxyyyy, tr_x_zz_xxyyyz, tr_x_zz_xxyyzz, tr_x_zz_xxyzzz, tr_x_zz_xxzzzz, tr_x_zz_xyyyyy, tr_x_zz_xyyyyz, tr_x_zz_xyyyzz, tr_x_zz_xyyzzz, tr_x_zz_xyzzzz, tr_x_zz_xzzzzz, tr_x_zz_yyyyyy, tr_x_zz_yyyyyz, tr_x_zz_yyyyzz, tr_x_zz_yyyzzz, tr_x_zz_yyzzzz, tr_x_zz_yzzzzz, tr_x_zz_zzzzzz, tr_x_zzz_xxxxx, tr_x_zzz_xxxxxx, tr_x_zzz_xxxxxy, tr_x_zzz_xxxxxz, tr_x_zzz_xxxxy, tr_x_zzz_xxxxyy, tr_x_zzz_xxxxyz, tr_x_zzz_xxxxz, tr_x_zzz_xxxxzz, tr_x_zzz_xxxyy, tr_x_zzz_xxxyyy, tr_x_zzz_xxxyyz, tr_x_zzz_xxxyz, tr_x_zzz_xxxyzz, tr_x_zzz_xxxzz, tr_x_zzz_xxxzzz, tr_x_zzz_xxyyy, tr_x_zzz_xxyyyy, tr_x_zzz_xxyyyz, tr_x_zzz_xxyyz, tr_x_zzz_xxyyzz, tr_x_zzz_xxyzz, tr_x_zzz_xxyzzz, tr_x_zzz_xxzzz, tr_x_zzz_xxzzzz, tr_x_zzz_xyyyy, tr_x_zzz_xyyyyy, tr_x_zzz_xyyyyz, tr_x_zzz_xyyyz, tr_x_zzz_xyyyzz, tr_x_zzz_xyyzz, tr_x_zzz_xyyzzz, tr_x_zzz_xyzzz, tr_x_zzz_xyzzzz, tr_x_zzz_xzzzz, tr_x_zzz_xzzzzz, tr_x_zzz_yyyyy, tr_x_zzz_yyyyyy, tr_x_zzz_yyyyyz, tr_x_zzz_yyyyz, tr_x_zzz_yyyyzz, tr_x_zzz_yyyzz, tr_x_zzz_yyyzzz, tr_x_zzz_yyzzz, tr_x_zzz_yyzzzz, tr_x_zzz_yzzzz, tr_x_zzz_yzzzzz, tr_x_zzz_zzzzz, tr_x_zzz_zzzzzz, tr_x_zzzz_xxxxxx, tr_x_zzzz_xxxxxy, tr_x_zzzz_xxxxxz, tr_x_zzzz_xxxxyy, tr_x_zzzz_xxxxyz, tr_x_zzzz_xxxxzz, tr_x_zzzz_xxxyyy, tr_x_zzzz_xxxyyz, tr_x_zzzz_xxxyzz, tr_x_zzzz_xxxzzz, tr_x_zzzz_xxyyyy, tr_x_zzzz_xxyyyz, tr_x_zzzz_xxyyzz, tr_x_zzzz_xxyzzz, tr_x_zzzz_xxzzzz, tr_x_zzzz_xyyyyy, tr_x_zzzz_xyyyyz, tr_x_zzzz_xyyyzz, tr_x_zzzz_xyyzzz, tr_x_zzzz_xyzzzz, tr_x_zzzz_xzzzzz, tr_x_zzzz_yyyyyy, tr_x_zzzz_yyyyyz, tr_x_zzzz_yyyyzz, tr_x_zzzz_yyyzzz, tr_x_zzzz_yyzzzz, tr_x_zzzz_yzzzzz, tr_x_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzz_xxxxxx[i] = 3.0 * tr_x_zz_xxxxxx[i] * fe_0 + tr_x_zzz_xxxxxx[i] * pa_z[i];

        tr_x_zzzz_xxxxxy[i] = 3.0 * tr_x_zz_xxxxxy[i] * fe_0 + tr_x_zzz_xxxxxy[i] * pa_z[i];

        tr_x_zzzz_xxxxxz[i] = 3.0 * tr_x_zz_xxxxxz[i] * fe_0 + tr_x_zzz_xxxxx[i] * fe_0 + tr_x_zzz_xxxxxz[i] * pa_z[i];

        tr_x_zzzz_xxxxyy[i] = 3.0 * tr_x_zz_xxxxyy[i] * fe_0 + tr_x_zzz_xxxxyy[i] * pa_z[i];

        tr_x_zzzz_xxxxyz[i] = 3.0 * tr_x_zz_xxxxyz[i] * fe_0 + tr_x_zzz_xxxxy[i] * fe_0 + tr_x_zzz_xxxxyz[i] * pa_z[i];

        tr_x_zzzz_xxxxzz[i] = 3.0 * tr_x_zz_xxxxzz[i] * fe_0 + 2.0 * tr_x_zzz_xxxxz[i] * fe_0 + tr_x_zzz_xxxxzz[i] * pa_z[i];

        tr_x_zzzz_xxxyyy[i] = 3.0 * tr_x_zz_xxxyyy[i] * fe_0 + tr_x_zzz_xxxyyy[i] * pa_z[i];

        tr_x_zzzz_xxxyyz[i] = 3.0 * tr_x_zz_xxxyyz[i] * fe_0 + tr_x_zzz_xxxyy[i] * fe_0 + tr_x_zzz_xxxyyz[i] * pa_z[i];

        tr_x_zzzz_xxxyzz[i] = 3.0 * tr_x_zz_xxxyzz[i] * fe_0 + 2.0 * tr_x_zzz_xxxyz[i] * fe_0 + tr_x_zzz_xxxyzz[i] * pa_z[i];

        tr_x_zzzz_xxxzzz[i] = 3.0 * tr_x_zz_xxxzzz[i] * fe_0 + 3.0 * tr_x_zzz_xxxzz[i] * fe_0 + tr_x_zzz_xxxzzz[i] * pa_z[i];

        tr_x_zzzz_xxyyyy[i] = 3.0 * tr_x_zz_xxyyyy[i] * fe_0 + tr_x_zzz_xxyyyy[i] * pa_z[i];

        tr_x_zzzz_xxyyyz[i] = 3.0 * tr_x_zz_xxyyyz[i] * fe_0 + tr_x_zzz_xxyyy[i] * fe_0 + tr_x_zzz_xxyyyz[i] * pa_z[i];

        tr_x_zzzz_xxyyzz[i] = 3.0 * tr_x_zz_xxyyzz[i] * fe_0 + 2.0 * tr_x_zzz_xxyyz[i] * fe_0 + tr_x_zzz_xxyyzz[i] * pa_z[i];

        tr_x_zzzz_xxyzzz[i] = 3.0 * tr_x_zz_xxyzzz[i] * fe_0 + 3.0 * tr_x_zzz_xxyzz[i] * fe_0 + tr_x_zzz_xxyzzz[i] * pa_z[i];

        tr_x_zzzz_xxzzzz[i] = 3.0 * tr_x_zz_xxzzzz[i] * fe_0 + 4.0 * tr_x_zzz_xxzzz[i] * fe_0 + tr_x_zzz_xxzzzz[i] * pa_z[i];

        tr_x_zzzz_xyyyyy[i] = 3.0 * tr_x_zz_xyyyyy[i] * fe_0 + tr_x_zzz_xyyyyy[i] * pa_z[i];

        tr_x_zzzz_xyyyyz[i] = 3.0 * tr_x_zz_xyyyyz[i] * fe_0 + tr_x_zzz_xyyyy[i] * fe_0 + tr_x_zzz_xyyyyz[i] * pa_z[i];

        tr_x_zzzz_xyyyzz[i] = 3.0 * tr_x_zz_xyyyzz[i] * fe_0 + 2.0 * tr_x_zzz_xyyyz[i] * fe_0 + tr_x_zzz_xyyyzz[i] * pa_z[i];

        tr_x_zzzz_xyyzzz[i] = 3.0 * tr_x_zz_xyyzzz[i] * fe_0 + 3.0 * tr_x_zzz_xyyzz[i] * fe_0 + tr_x_zzz_xyyzzz[i] * pa_z[i];

        tr_x_zzzz_xyzzzz[i] = 3.0 * tr_x_zz_xyzzzz[i] * fe_0 + 4.0 * tr_x_zzz_xyzzz[i] * fe_0 + tr_x_zzz_xyzzzz[i] * pa_z[i];

        tr_x_zzzz_xzzzzz[i] = 3.0 * tr_x_zz_xzzzzz[i] * fe_0 + 5.0 * tr_x_zzz_xzzzz[i] * fe_0 + tr_x_zzz_xzzzzz[i] * pa_z[i];

        tr_x_zzzz_yyyyyy[i] = 3.0 * tr_x_zz_yyyyyy[i] * fe_0 + tr_x_zzz_yyyyyy[i] * pa_z[i];

        tr_x_zzzz_yyyyyz[i] = 3.0 * tr_x_zz_yyyyyz[i] * fe_0 + tr_x_zzz_yyyyy[i] * fe_0 + tr_x_zzz_yyyyyz[i] * pa_z[i];

        tr_x_zzzz_yyyyzz[i] = 3.0 * tr_x_zz_yyyyzz[i] * fe_0 + 2.0 * tr_x_zzz_yyyyz[i] * fe_0 + tr_x_zzz_yyyyzz[i] * pa_z[i];

        tr_x_zzzz_yyyzzz[i] = 3.0 * tr_x_zz_yyyzzz[i] * fe_0 + 3.0 * tr_x_zzz_yyyzz[i] * fe_0 + tr_x_zzz_yyyzzz[i] * pa_z[i];

        tr_x_zzzz_yyzzzz[i] = 3.0 * tr_x_zz_yyzzzz[i] * fe_0 + 4.0 * tr_x_zzz_yyzzz[i] * fe_0 + tr_x_zzz_yyzzzz[i] * pa_z[i];

        tr_x_zzzz_yzzzzz[i] = 3.0 * tr_x_zz_yzzzzz[i] * fe_0 + 5.0 * tr_x_zzz_yzzzz[i] * fe_0 + tr_x_zzz_yzzzzz[i] * pa_z[i];

        tr_x_zzzz_zzzzzz[i] = 3.0 * tr_x_zz_zzzzzz[i] * fe_0 + 6.0 * tr_x_zzz_zzzzz[i] * fe_0 + tr_x_zzz_zzzzzz[i] * pa_z[i];
    }

    // Set up 420-448 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_y_xx_xxxxxx, tr_y_xx_xxxxxy, tr_y_xx_xxxxxz, tr_y_xx_xxxxyy, tr_y_xx_xxxxyz, tr_y_xx_xxxxzz, tr_y_xx_xxxyyy, tr_y_xx_xxxyyz, tr_y_xx_xxxyzz, tr_y_xx_xxxzzz, tr_y_xx_xxyyyy, tr_y_xx_xxyyyz, tr_y_xx_xxyyzz, tr_y_xx_xxyzzz, tr_y_xx_xxzzzz, tr_y_xx_xyyyyy, tr_y_xx_xyyyyz, tr_y_xx_xyyyzz, tr_y_xx_xyyzzz, tr_y_xx_xyzzzz, tr_y_xx_xzzzzz, tr_y_xx_yyyyyy, tr_y_xx_yyyyyz, tr_y_xx_yyyyzz, tr_y_xx_yyyzzz, tr_y_xx_yyzzzz, tr_y_xx_yzzzzz, tr_y_xx_zzzzzz, tr_y_xxx_xxxxx, tr_y_xxx_xxxxxx, tr_y_xxx_xxxxxy, tr_y_xxx_xxxxxz, tr_y_xxx_xxxxy, tr_y_xxx_xxxxyy, tr_y_xxx_xxxxyz, tr_y_xxx_xxxxz, tr_y_xxx_xxxxzz, tr_y_xxx_xxxyy, tr_y_xxx_xxxyyy, tr_y_xxx_xxxyyz, tr_y_xxx_xxxyz, tr_y_xxx_xxxyzz, tr_y_xxx_xxxzz, tr_y_xxx_xxxzzz, tr_y_xxx_xxyyy, tr_y_xxx_xxyyyy, tr_y_xxx_xxyyyz, tr_y_xxx_xxyyz, tr_y_xxx_xxyyzz, tr_y_xxx_xxyzz, tr_y_xxx_xxyzzz, tr_y_xxx_xxzzz, tr_y_xxx_xxzzzz, tr_y_xxx_xyyyy, tr_y_xxx_xyyyyy, tr_y_xxx_xyyyyz, tr_y_xxx_xyyyz, tr_y_xxx_xyyyzz, tr_y_xxx_xyyzz, tr_y_xxx_xyyzzz, tr_y_xxx_xyzzz, tr_y_xxx_xyzzzz, tr_y_xxx_xzzzz, tr_y_xxx_xzzzzz, tr_y_xxx_yyyyy, tr_y_xxx_yyyyyy, tr_y_xxx_yyyyyz, tr_y_xxx_yyyyz, tr_y_xxx_yyyyzz, tr_y_xxx_yyyzz, tr_y_xxx_yyyzzz, tr_y_xxx_yyzzz, tr_y_xxx_yyzzzz, tr_y_xxx_yzzzz, tr_y_xxx_yzzzzz, tr_y_xxx_zzzzz, tr_y_xxx_zzzzzz, tr_y_xxxx_xxxxxx, tr_y_xxxx_xxxxxy, tr_y_xxxx_xxxxxz, tr_y_xxxx_xxxxyy, tr_y_xxxx_xxxxyz, tr_y_xxxx_xxxxzz, tr_y_xxxx_xxxyyy, tr_y_xxxx_xxxyyz, tr_y_xxxx_xxxyzz, tr_y_xxxx_xxxzzz, tr_y_xxxx_xxyyyy, tr_y_xxxx_xxyyyz, tr_y_xxxx_xxyyzz, tr_y_xxxx_xxyzzz, tr_y_xxxx_xxzzzz, tr_y_xxxx_xyyyyy, tr_y_xxxx_xyyyyz, tr_y_xxxx_xyyyzz, tr_y_xxxx_xyyzzz, tr_y_xxxx_xyzzzz, tr_y_xxxx_xzzzzz, tr_y_xxxx_yyyyyy, tr_y_xxxx_yyyyyz, tr_y_xxxx_yyyyzz, tr_y_xxxx_yyyzzz, tr_y_xxxx_yyzzzz, tr_y_xxxx_yzzzzz, tr_y_xxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxx_xxxxxx[i] = 3.0 * tr_y_xx_xxxxxx[i] * fe_0 + 6.0 * tr_y_xxx_xxxxx[i] * fe_0 + tr_y_xxx_xxxxxx[i] * pa_x[i];

        tr_y_xxxx_xxxxxy[i] = 3.0 * tr_y_xx_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxx_xxxxy[i] * fe_0 + tr_y_xxx_xxxxxy[i] * pa_x[i];

        tr_y_xxxx_xxxxxz[i] = 3.0 * tr_y_xx_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxx_xxxxz[i] * fe_0 + tr_y_xxx_xxxxxz[i] * pa_x[i];

        tr_y_xxxx_xxxxyy[i] = 3.0 * tr_y_xx_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxx_xxxyy[i] * fe_0 + tr_y_xxx_xxxxyy[i] * pa_x[i];

        tr_y_xxxx_xxxxyz[i] = 3.0 * tr_y_xx_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxx_xxxyz[i] * fe_0 + tr_y_xxx_xxxxyz[i] * pa_x[i];

        tr_y_xxxx_xxxxzz[i] = 3.0 * tr_y_xx_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxx_xxxzz[i] * fe_0 + tr_y_xxx_xxxxzz[i] * pa_x[i];

        tr_y_xxxx_xxxyyy[i] = 3.0 * tr_y_xx_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxx_xxyyy[i] * fe_0 + tr_y_xxx_xxxyyy[i] * pa_x[i];

        tr_y_xxxx_xxxyyz[i] = 3.0 * tr_y_xx_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxx_xxyyz[i] * fe_0 + tr_y_xxx_xxxyyz[i] * pa_x[i];

        tr_y_xxxx_xxxyzz[i] = 3.0 * tr_y_xx_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxx_xxyzz[i] * fe_0 + tr_y_xxx_xxxyzz[i] * pa_x[i];

        tr_y_xxxx_xxxzzz[i] = 3.0 * tr_y_xx_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxx_xxzzz[i] * fe_0 + tr_y_xxx_xxxzzz[i] * pa_x[i];

        tr_y_xxxx_xxyyyy[i] = 3.0 * tr_y_xx_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxx_xyyyy[i] * fe_0 + tr_y_xxx_xxyyyy[i] * pa_x[i];

        tr_y_xxxx_xxyyyz[i] = 3.0 * tr_y_xx_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxx_xyyyz[i] * fe_0 + tr_y_xxx_xxyyyz[i] * pa_x[i];

        tr_y_xxxx_xxyyzz[i] = 3.0 * tr_y_xx_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxx_xyyzz[i] * fe_0 + tr_y_xxx_xxyyzz[i] * pa_x[i];

        tr_y_xxxx_xxyzzz[i] = 3.0 * tr_y_xx_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxx_xyzzz[i] * fe_0 + tr_y_xxx_xxyzzz[i] * pa_x[i];

        tr_y_xxxx_xxzzzz[i] = 3.0 * tr_y_xx_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxx_xzzzz[i] * fe_0 + tr_y_xxx_xxzzzz[i] * pa_x[i];

        tr_y_xxxx_xyyyyy[i] = 3.0 * tr_y_xx_xyyyyy[i] * fe_0 + tr_y_xxx_yyyyy[i] * fe_0 + tr_y_xxx_xyyyyy[i] * pa_x[i];

        tr_y_xxxx_xyyyyz[i] = 3.0 * tr_y_xx_xyyyyz[i] * fe_0 + tr_y_xxx_yyyyz[i] * fe_0 + tr_y_xxx_xyyyyz[i] * pa_x[i];

        tr_y_xxxx_xyyyzz[i] = 3.0 * tr_y_xx_xyyyzz[i] * fe_0 + tr_y_xxx_yyyzz[i] * fe_0 + tr_y_xxx_xyyyzz[i] * pa_x[i];

        tr_y_xxxx_xyyzzz[i] = 3.0 * tr_y_xx_xyyzzz[i] * fe_0 + tr_y_xxx_yyzzz[i] * fe_0 + tr_y_xxx_xyyzzz[i] * pa_x[i];

        tr_y_xxxx_xyzzzz[i] = 3.0 * tr_y_xx_xyzzzz[i] * fe_0 + tr_y_xxx_yzzzz[i] * fe_0 + tr_y_xxx_xyzzzz[i] * pa_x[i];

        tr_y_xxxx_xzzzzz[i] = 3.0 * tr_y_xx_xzzzzz[i] * fe_0 + tr_y_xxx_zzzzz[i] * fe_0 + tr_y_xxx_xzzzzz[i] * pa_x[i];

        tr_y_xxxx_yyyyyy[i] = 3.0 * tr_y_xx_yyyyyy[i] * fe_0 + tr_y_xxx_yyyyyy[i] * pa_x[i];

        tr_y_xxxx_yyyyyz[i] = 3.0 * tr_y_xx_yyyyyz[i] * fe_0 + tr_y_xxx_yyyyyz[i] * pa_x[i];

        tr_y_xxxx_yyyyzz[i] = 3.0 * tr_y_xx_yyyyzz[i] * fe_0 + tr_y_xxx_yyyyzz[i] * pa_x[i];

        tr_y_xxxx_yyyzzz[i] = 3.0 * tr_y_xx_yyyzzz[i] * fe_0 + tr_y_xxx_yyyzzz[i] * pa_x[i];

        tr_y_xxxx_yyzzzz[i] = 3.0 * tr_y_xx_yyzzzz[i] * fe_0 + tr_y_xxx_yyzzzz[i] * pa_x[i];

        tr_y_xxxx_yzzzzz[i] = 3.0 * tr_y_xx_yzzzzz[i] * fe_0 + tr_y_xxx_yzzzzz[i] * pa_x[i];

        tr_y_xxxx_zzzzzz[i] = 3.0 * tr_y_xx_zzzzzz[i] * fe_0 + tr_y_xxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 448-476 components of targeted buffer : GI

    auto tr_y_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 448);

    auto tr_y_xxxy_xxxxxy = pbuffer.data(idx_dip_gi + 449);

    auto tr_y_xxxy_xxxxxz = pbuffer.data(idx_dip_gi + 450);

    auto tr_y_xxxy_xxxxyy = pbuffer.data(idx_dip_gi + 451);

    auto tr_y_xxxy_xxxxyz = pbuffer.data(idx_dip_gi + 452);

    auto tr_y_xxxy_xxxxzz = pbuffer.data(idx_dip_gi + 453);

    auto tr_y_xxxy_xxxyyy = pbuffer.data(idx_dip_gi + 454);

    auto tr_y_xxxy_xxxyyz = pbuffer.data(idx_dip_gi + 455);

    auto tr_y_xxxy_xxxyzz = pbuffer.data(idx_dip_gi + 456);

    auto tr_y_xxxy_xxxzzz = pbuffer.data(idx_dip_gi + 457);

    auto tr_y_xxxy_xxyyyy = pbuffer.data(idx_dip_gi + 458);

    auto tr_y_xxxy_xxyyyz = pbuffer.data(idx_dip_gi + 459);

    auto tr_y_xxxy_xxyyzz = pbuffer.data(idx_dip_gi + 460);

    auto tr_y_xxxy_xxyzzz = pbuffer.data(idx_dip_gi + 461);

    auto tr_y_xxxy_xxzzzz = pbuffer.data(idx_dip_gi + 462);

    auto tr_y_xxxy_xyyyyy = pbuffer.data(idx_dip_gi + 463);

    auto tr_y_xxxy_xyyyyz = pbuffer.data(idx_dip_gi + 464);

    auto tr_y_xxxy_xyyyzz = pbuffer.data(idx_dip_gi + 465);

    auto tr_y_xxxy_xyyzzz = pbuffer.data(idx_dip_gi + 466);

    auto tr_y_xxxy_xyzzzz = pbuffer.data(idx_dip_gi + 467);

    auto tr_y_xxxy_xzzzzz = pbuffer.data(idx_dip_gi + 468);

    auto tr_y_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 469);

    auto tr_y_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 470);

    auto tr_y_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 471);

    auto tr_y_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 472);

    auto tr_y_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 473);

    auto tr_y_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 474);

    auto tr_y_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 475);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xxx_xxxxxx, tr_y_xxx_xxxxxz, tr_y_xxx_xxxxzz, tr_y_xxx_xxxzzz, tr_y_xxx_xxzzzz, tr_y_xxx_xzzzzz, tr_y_xxxy_xxxxxx, tr_y_xxxy_xxxxxy, tr_y_xxxy_xxxxxz, tr_y_xxxy_xxxxyy, tr_y_xxxy_xxxxyz, tr_y_xxxy_xxxxzz, tr_y_xxxy_xxxyyy, tr_y_xxxy_xxxyyz, tr_y_xxxy_xxxyzz, tr_y_xxxy_xxxzzz, tr_y_xxxy_xxyyyy, tr_y_xxxy_xxyyyz, tr_y_xxxy_xxyyzz, tr_y_xxxy_xxyzzz, tr_y_xxxy_xxzzzz, tr_y_xxxy_xyyyyy, tr_y_xxxy_xyyyyz, tr_y_xxxy_xyyyzz, tr_y_xxxy_xyyzzz, tr_y_xxxy_xyzzzz, tr_y_xxxy_xzzzzz, tr_y_xxxy_yyyyyy, tr_y_xxxy_yyyyyz, tr_y_xxxy_yyyyzz, tr_y_xxxy_yyyzzz, tr_y_xxxy_yyzzzz, tr_y_xxxy_yzzzzz, tr_y_xxxy_zzzzzz, tr_y_xxy_xxxxxy, tr_y_xxy_xxxxy, tr_y_xxy_xxxxyy, tr_y_xxy_xxxxyz, tr_y_xxy_xxxyy, tr_y_xxy_xxxyyy, tr_y_xxy_xxxyyz, tr_y_xxy_xxxyz, tr_y_xxy_xxxyzz, tr_y_xxy_xxyyy, tr_y_xxy_xxyyyy, tr_y_xxy_xxyyyz, tr_y_xxy_xxyyz, tr_y_xxy_xxyyzz, tr_y_xxy_xxyzz, tr_y_xxy_xxyzzz, tr_y_xxy_xyyyy, tr_y_xxy_xyyyyy, tr_y_xxy_xyyyyz, tr_y_xxy_xyyyz, tr_y_xxy_xyyyzz, tr_y_xxy_xyyzz, tr_y_xxy_xyyzzz, tr_y_xxy_xyzzz, tr_y_xxy_xyzzzz, tr_y_xxy_yyyyy, tr_y_xxy_yyyyyy, tr_y_xxy_yyyyyz, tr_y_xxy_yyyyz, tr_y_xxy_yyyyzz, tr_y_xxy_yyyzz, tr_y_xxy_yyyzzz, tr_y_xxy_yyzzz, tr_y_xxy_yyzzzz, tr_y_xxy_yzzzz, tr_y_xxy_yzzzzz, tr_y_xxy_zzzzzz, tr_y_xy_xxxxxy, tr_y_xy_xxxxyy, tr_y_xy_xxxxyz, tr_y_xy_xxxyyy, tr_y_xy_xxxyyz, tr_y_xy_xxxyzz, tr_y_xy_xxyyyy, tr_y_xy_xxyyyz, tr_y_xy_xxyyzz, tr_y_xy_xxyzzz, tr_y_xy_xyyyyy, tr_y_xy_xyyyyz, tr_y_xy_xyyyzz, tr_y_xy_xyyzzz, tr_y_xy_xyzzzz, tr_y_xy_yyyyyy, tr_y_xy_yyyyyz, tr_y_xy_yyyyzz, tr_y_xy_yyyzzz, tr_y_xy_yyzzzz, tr_y_xy_yzzzzz, tr_y_xy_zzzzzz, ts_xxx_xxxxxx, ts_xxx_xxxxxz, ts_xxx_xxxxzz, ts_xxx_xxxzzz, ts_xxx_xxzzzz, ts_xxx_xzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxy_xxxxxx[i] = ts_xxx_xxxxxx[i] * fe_0 + tr_y_xxx_xxxxxx[i] * pa_y[i];

        tr_y_xxxy_xxxxxy[i] = 2.0 * tr_y_xy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxy_xxxxy[i] * fe_0 + tr_y_xxy_xxxxxy[i] * pa_x[i];

        tr_y_xxxy_xxxxxz[i] = ts_xxx_xxxxxz[i] * fe_0 + tr_y_xxx_xxxxxz[i] * pa_y[i];

        tr_y_xxxy_xxxxyy[i] = 2.0 * tr_y_xy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxy_xxxyy[i] * fe_0 + tr_y_xxy_xxxxyy[i] * pa_x[i];

        tr_y_xxxy_xxxxyz[i] = 2.0 * tr_y_xy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxy_xxxyz[i] * fe_0 + tr_y_xxy_xxxxyz[i] * pa_x[i];

        tr_y_xxxy_xxxxzz[i] = ts_xxx_xxxxzz[i] * fe_0 + tr_y_xxx_xxxxzz[i] * pa_y[i];

        tr_y_xxxy_xxxyyy[i] = 2.0 * tr_y_xy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxy_xxyyy[i] * fe_0 + tr_y_xxy_xxxyyy[i] * pa_x[i];

        tr_y_xxxy_xxxyyz[i] = 2.0 * tr_y_xy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxy_xxyyz[i] * fe_0 + tr_y_xxy_xxxyyz[i] * pa_x[i];

        tr_y_xxxy_xxxyzz[i] = 2.0 * tr_y_xy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxy_xxyzz[i] * fe_0 + tr_y_xxy_xxxyzz[i] * pa_x[i];

        tr_y_xxxy_xxxzzz[i] = ts_xxx_xxxzzz[i] * fe_0 + tr_y_xxx_xxxzzz[i] * pa_y[i];

        tr_y_xxxy_xxyyyy[i] = 2.0 * tr_y_xy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxy_xyyyy[i] * fe_0 + tr_y_xxy_xxyyyy[i] * pa_x[i];

        tr_y_xxxy_xxyyyz[i] = 2.0 * tr_y_xy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxy_xyyyz[i] * fe_0 + tr_y_xxy_xxyyyz[i] * pa_x[i];

        tr_y_xxxy_xxyyzz[i] = 2.0 * tr_y_xy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxy_xyyzz[i] * fe_0 + tr_y_xxy_xxyyzz[i] * pa_x[i];

        tr_y_xxxy_xxyzzz[i] = 2.0 * tr_y_xy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxy_xyzzz[i] * fe_0 + tr_y_xxy_xxyzzz[i] * pa_x[i];

        tr_y_xxxy_xxzzzz[i] = ts_xxx_xxzzzz[i] * fe_0 + tr_y_xxx_xxzzzz[i] * pa_y[i];

        tr_y_xxxy_xyyyyy[i] = 2.0 * tr_y_xy_xyyyyy[i] * fe_0 + tr_y_xxy_yyyyy[i] * fe_0 + tr_y_xxy_xyyyyy[i] * pa_x[i];

        tr_y_xxxy_xyyyyz[i] = 2.0 * tr_y_xy_xyyyyz[i] * fe_0 + tr_y_xxy_yyyyz[i] * fe_0 + tr_y_xxy_xyyyyz[i] * pa_x[i];

        tr_y_xxxy_xyyyzz[i] = 2.0 * tr_y_xy_xyyyzz[i] * fe_0 + tr_y_xxy_yyyzz[i] * fe_0 + tr_y_xxy_xyyyzz[i] * pa_x[i];

        tr_y_xxxy_xyyzzz[i] = 2.0 * tr_y_xy_xyyzzz[i] * fe_0 + tr_y_xxy_yyzzz[i] * fe_0 + tr_y_xxy_xyyzzz[i] * pa_x[i];

        tr_y_xxxy_xyzzzz[i] = 2.0 * tr_y_xy_xyzzzz[i] * fe_0 + tr_y_xxy_yzzzz[i] * fe_0 + tr_y_xxy_xyzzzz[i] * pa_x[i];

        tr_y_xxxy_xzzzzz[i] = ts_xxx_xzzzzz[i] * fe_0 + tr_y_xxx_xzzzzz[i] * pa_y[i];

        tr_y_xxxy_yyyyyy[i] = 2.0 * tr_y_xy_yyyyyy[i] * fe_0 + tr_y_xxy_yyyyyy[i] * pa_x[i];

        tr_y_xxxy_yyyyyz[i] = 2.0 * tr_y_xy_yyyyyz[i] * fe_0 + tr_y_xxy_yyyyyz[i] * pa_x[i];

        tr_y_xxxy_yyyyzz[i] = 2.0 * tr_y_xy_yyyyzz[i] * fe_0 + tr_y_xxy_yyyyzz[i] * pa_x[i];

        tr_y_xxxy_yyyzzz[i] = 2.0 * tr_y_xy_yyyzzz[i] * fe_0 + tr_y_xxy_yyyzzz[i] * pa_x[i];

        tr_y_xxxy_yyzzzz[i] = 2.0 * tr_y_xy_yyzzzz[i] * fe_0 + tr_y_xxy_yyzzzz[i] * pa_x[i];

        tr_y_xxxy_yzzzzz[i] = 2.0 * tr_y_xy_yzzzzz[i] * fe_0 + tr_y_xxy_yzzzzz[i] * pa_x[i];

        tr_y_xxxy_zzzzzz[i] = 2.0 * tr_y_xy_zzzzzz[i] * fe_0 + tr_y_xxy_zzzzzz[i] * pa_x[i];
    }

    // Set up 476-504 components of targeted buffer : GI

    auto tr_y_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 476);

    auto tr_y_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 477);

    auto tr_y_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 478);

    auto tr_y_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 479);

    auto tr_y_xxxz_xxxxyz = pbuffer.data(idx_dip_gi + 480);

    auto tr_y_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 481);

    auto tr_y_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 482);

    auto tr_y_xxxz_xxxyyz = pbuffer.data(idx_dip_gi + 483);

    auto tr_y_xxxz_xxxyzz = pbuffer.data(idx_dip_gi + 484);

    auto tr_y_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 485);

    auto tr_y_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 486);

    auto tr_y_xxxz_xxyyyz = pbuffer.data(idx_dip_gi + 487);

    auto tr_y_xxxz_xxyyzz = pbuffer.data(idx_dip_gi + 488);

    auto tr_y_xxxz_xxyzzz = pbuffer.data(idx_dip_gi + 489);

    auto tr_y_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 490);

    auto tr_y_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 491);

    auto tr_y_xxxz_xyyyyz = pbuffer.data(idx_dip_gi + 492);

    auto tr_y_xxxz_xyyyzz = pbuffer.data(idx_dip_gi + 493);

    auto tr_y_xxxz_xyyzzz = pbuffer.data(idx_dip_gi + 494);

    auto tr_y_xxxz_xyzzzz = pbuffer.data(idx_dip_gi + 495);

    auto tr_y_xxxz_xzzzzz = pbuffer.data(idx_dip_gi + 496);

    auto tr_y_xxxz_yyyyyy = pbuffer.data(idx_dip_gi + 497);

    auto tr_y_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 498);

    auto tr_y_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 499);

    auto tr_y_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 500);

    auto tr_y_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 501);

    auto tr_y_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 502);

    auto tr_y_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 503);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxx_xxxxx, tr_y_xxx_xxxxxx, tr_y_xxx_xxxxxy, tr_y_xxx_xxxxxz, tr_y_xxx_xxxxy, tr_y_xxx_xxxxyy, tr_y_xxx_xxxxyz, tr_y_xxx_xxxxz, tr_y_xxx_xxxxzz, tr_y_xxx_xxxyy, tr_y_xxx_xxxyyy, tr_y_xxx_xxxyyz, tr_y_xxx_xxxyz, tr_y_xxx_xxxyzz, tr_y_xxx_xxxzz, tr_y_xxx_xxxzzz, tr_y_xxx_xxyyy, tr_y_xxx_xxyyyy, tr_y_xxx_xxyyyz, tr_y_xxx_xxyyz, tr_y_xxx_xxyyzz, tr_y_xxx_xxyzz, tr_y_xxx_xxyzzz, tr_y_xxx_xxzzz, tr_y_xxx_xxzzzz, tr_y_xxx_xyyyy, tr_y_xxx_xyyyyy, tr_y_xxx_xyyyyz, tr_y_xxx_xyyyz, tr_y_xxx_xyyyzz, tr_y_xxx_xyyzz, tr_y_xxx_xyyzzz, tr_y_xxx_xyzzz, tr_y_xxx_xyzzzz, tr_y_xxx_xzzzz, tr_y_xxx_xzzzzz, tr_y_xxx_yyyyyy, tr_y_xxxz_xxxxxx, tr_y_xxxz_xxxxxy, tr_y_xxxz_xxxxxz, tr_y_xxxz_xxxxyy, tr_y_xxxz_xxxxyz, tr_y_xxxz_xxxxzz, tr_y_xxxz_xxxyyy, tr_y_xxxz_xxxyyz, tr_y_xxxz_xxxyzz, tr_y_xxxz_xxxzzz, tr_y_xxxz_xxyyyy, tr_y_xxxz_xxyyyz, tr_y_xxxz_xxyyzz, tr_y_xxxz_xxyzzz, tr_y_xxxz_xxzzzz, tr_y_xxxz_xyyyyy, tr_y_xxxz_xyyyyz, tr_y_xxxz_xyyyzz, tr_y_xxxz_xyyzzz, tr_y_xxxz_xyzzzz, tr_y_xxxz_xzzzzz, tr_y_xxxz_yyyyyy, tr_y_xxxz_yyyyyz, tr_y_xxxz_yyyyzz, tr_y_xxxz_yyyzzz, tr_y_xxxz_yyzzzz, tr_y_xxxz_yzzzzz, tr_y_xxxz_zzzzzz, tr_y_xxz_yyyyyz, tr_y_xxz_yyyyzz, tr_y_xxz_yyyzzz, tr_y_xxz_yyzzzz, tr_y_xxz_yzzzzz, tr_y_xxz_zzzzzz, tr_y_xz_yyyyyz, tr_y_xz_yyyyzz, tr_y_xz_yyyzzz, tr_y_xz_yyzzzz, tr_y_xz_yzzzzz, tr_y_xz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxz_xxxxxx[i] = tr_y_xxx_xxxxxx[i] * pa_z[i];

        tr_y_xxxz_xxxxxy[i] = tr_y_xxx_xxxxxy[i] * pa_z[i];

        tr_y_xxxz_xxxxxz[i] = tr_y_xxx_xxxxx[i] * fe_0 + tr_y_xxx_xxxxxz[i] * pa_z[i];

        tr_y_xxxz_xxxxyy[i] = tr_y_xxx_xxxxyy[i] * pa_z[i];

        tr_y_xxxz_xxxxyz[i] = tr_y_xxx_xxxxy[i] * fe_0 + tr_y_xxx_xxxxyz[i] * pa_z[i];

        tr_y_xxxz_xxxxzz[i] = 2.0 * tr_y_xxx_xxxxz[i] * fe_0 + tr_y_xxx_xxxxzz[i] * pa_z[i];

        tr_y_xxxz_xxxyyy[i] = tr_y_xxx_xxxyyy[i] * pa_z[i];

        tr_y_xxxz_xxxyyz[i] = tr_y_xxx_xxxyy[i] * fe_0 + tr_y_xxx_xxxyyz[i] * pa_z[i];

        tr_y_xxxz_xxxyzz[i] = 2.0 * tr_y_xxx_xxxyz[i] * fe_0 + tr_y_xxx_xxxyzz[i] * pa_z[i];

        tr_y_xxxz_xxxzzz[i] = 3.0 * tr_y_xxx_xxxzz[i] * fe_0 + tr_y_xxx_xxxzzz[i] * pa_z[i];

        tr_y_xxxz_xxyyyy[i] = tr_y_xxx_xxyyyy[i] * pa_z[i];

        tr_y_xxxz_xxyyyz[i] = tr_y_xxx_xxyyy[i] * fe_0 + tr_y_xxx_xxyyyz[i] * pa_z[i];

        tr_y_xxxz_xxyyzz[i] = 2.0 * tr_y_xxx_xxyyz[i] * fe_0 + tr_y_xxx_xxyyzz[i] * pa_z[i];

        tr_y_xxxz_xxyzzz[i] = 3.0 * tr_y_xxx_xxyzz[i] * fe_0 + tr_y_xxx_xxyzzz[i] * pa_z[i];

        tr_y_xxxz_xxzzzz[i] = 4.0 * tr_y_xxx_xxzzz[i] * fe_0 + tr_y_xxx_xxzzzz[i] * pa_z[i];

        tr_y_xxxz_xyyyyy[i] = tr_y_xxx_xyyyyy[i] * pa_z[i];

        tr_y_xxxz_xyyyyz[i] = tr_y_xxx_xyyyy[i] * fe_0 + tr_y_xxx_xyyyyz[i] * pa_z[i];

        tr_y_xxxz_xyyyzz[i] = 2.0 * tr_y_xxx_xyyyz[i] * fe_0 + tr_y_xxx_xyyyzz[i] * pa_z[i];

        tr_y_xxxz_xyyzzz[i] = 3.0 * tr_y_xxx_xyyzz[i] * fe_0 + tr_y_xxx_xyyzzz[i] * pa_z[i];

        tr_y_xxxz_xyzzzz[i] = 4.0 * tr_y_xxx_xyzzz[i] * fe_0 + tr_y_xxx_xyzzzz[i] * pa_z[i];

        tr_y_xxxz_xzzzzz[i] = 5.0 * tr_y_xxx_xzzzz[i] * fe_0 + tr_y_xxx_xzzzzz[i] * pa_z[i];

        tr_y_xxxz_yyyyyy[i] = tr_y_xxx_yyyyyy[i] * pa_z[i];

        tr_y_xxxz_yyyyyz[i] = 2.0 * tr_y_xz_yyyyyz[i] * fe_0 + tr_y_xxz_yyyyyz[i] * pa_x[i];

        tr_y_xxxz_yyyyzz[i] = 2.0 * tr_y_xz_yyyyzz[i] * fe_0 + tr_y_xxz_yyyyzz[i] * pa_x[i];

        tr_y_xxxz_yyyzzz[i] = 2.0 * tr_y_xz_yyyzzz[i] * fe_0 + tr_y_xxz_yyyzzz[i] * pa_x[i];

        tr_y_xxxz_yyzzzz[i] = 2.0 * tr_y_xz_yyzzzz[i] * fe_0 + tr_y_xxz_yyzzzz[i] * pa_x[i];

        tr_y_xxxz_yzzzzz[i] = 2.0 * tr_y_xz_yzzzzz[i] * fe_0 + tr_y_xxz_yzzzzz[i] * pa_x[i];

        tr_y_xxxz_zzzzzz[i] = 2.0 * tr_y_xz_zzzzzz[i] * fe_0 + tr_y_xxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 504-532 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_y_xxyy_xxxxxx, tr_y_xxyy_xxxxxy, tr_y_xxyy_xxxxxz, tr_y_xxyy_xxxxyy, tr_y_xxyy_xxxxyz, tr_y_xxyy_xxxxzz, tr_y_xxyy_xxxyyy, tr_y_xxyy_xxxyyz, tr_y_xxyy_xxxyzz, tr_y_xxyy_xxxzzz, tr_y_xxyy_xxyyyy, tr_y_xxyy_xxyyyz, tr_y_xxyy_xxyyzz, tr_y_xxyy_xxyzzz, tr_y_xxyy_xxzzzz, tr_y_xxyy_xyyyyy, tr_y_xxyy_xyyyyz, tr_y_xxyy_xyyyzz, tr_y_xxyy_xyyzzz, tr_y_xxyy_xyzzzz, tr_y_xxyy_xzzzzz, tr_y_xxyy_yyyyyy, tr_y_xxyy_yyyyyz, tr_y_xxyy_yyyyzz, tr_y_xxyy_yyyzzz, tr_y_xxyy_yyzzzz, tr_y_xxyy_yzzzzz, tr_y_xxyy_zzzzzz, tr_y_xyy_xxxxx, tr_y_xyy_xxxxxx, tr_y_xyy_xxxxxy, tr_y_xyy_xxxxxz, tr_y_xyy_xxxxy, tr_y_xyy_xxxxyy, tr_y_xyy_xxxxyz, tr_y_xyy_xxxxz, tr_y_xyy_xxxxzz, tr_y_xyy_xxxyy, tr_y_xyy_xxxyyy, tr_y_xyy_xxxyyz, tr_y_xyy_xxxyz, tr_y_xyy_xxxyzz, tr_y_xyy_xxxzz, tr_y_xyy_xxxzzz, tr_y_xyy_xxyyy, tr_y_xyy_xxyyyy, tr_y_xyy_xxyyyz, tr_y_xyy_xxyyz, tr_y_xyy_xxyyzz, tr_y_xyy_xxyzz, tr_y_xyy_xxyzzz, tr_y_xyy_xxzzz, tr_y_xyy_xxzzzz, tr_y_xyy_xyyyy, tr_y_xyy_xyyyyy, tr_y_xyy_xyyyyz, tr_y_xyy_xyyyz, tr_y_xyy_xyyyzz, tr_y_xyy_xyyzz, tr_y_xyy_xyyzzz, tr_y_xyy_xyzzz, tr_y_xyy_xyzzzz, tr_y_xyy_xzzzz, tr_y_xyy_xzzzzz, tr_y_xyy_yyyyy, tr_y_xyy_yyyyyy, tr_y_xyy_yyyyyz, tr_y_xyy_yyyyz, tr_y_xyy_yyyyzz, tr_y_xyy_yyyzz, tr_y_xyy_yyyzzz, tr_y_xyy_yyzzz, tr_y_xyy_yyzzzz, tr_y_xyy_yzzzz, tr_y_xyy_yzzzzz, tr_y_xyy_zzzzz, tr_y_xyy_zzzzzz, tr_y_yy_xxxxxx, tr_y_yy_xxxxxy, tr_y_yy_xxxxxz, tr_y_yy_xxxxyy, tr_y_yy_xxxxyz, tr_y_yy_xxxxzz, tr_y_yy_xxxyyy, tr_y_yy_xxxyyz, tr_y_yy_xxxyzz, tr_y_yy_xxxzzz, tr_y_yy_xxyyyy, tr_y_yy_xxyyyz, tr_y_yy_xxyyzz, tr_y_yy_xxyzzz, tr_y_yy_xxzzzz, tr_y_yy_xyyyyy, tr_y_yy_xyyyyz, tr_y_yy_xyyyzz, tr_y_yy_xyyzzz, tr_y_yy_xyzzzz, tr_y_yy_xzzzzz, tr_y_yy_yyyyyy, tr_y_yy_yyyyyz, tr_y_yy_yyyyzz, tr_y_yy_yyyzzz, tr_y_yy_yyzzzz, tr_y_yy_yzzzzz, tr_y_yy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyy_xxxxxx[i] = tr_y_yy_xxxxxx[i] * fe_0 + 6.0 * tr_y_xyy_xxxxx[i] * fe_0 + tr_y_xyy_xxxxxx[i] * pa_x[i];

        tr_y_xxyy_xxxxxy[i] = tr_y_yy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xyy_xxxxy[i] * fe_0 + tr_y_xyy_xxxxxy[i] * pa_x[i];

        tr_y_xxyy_xxxxxz[i] = tr_y_yy_xxxxxz[i] * fe_0 + 5.0 * tr_y_xyy_xxxxz[i] * fe_0 + tr_y_xyy_xxxxxz[i] * pa_x[i];

        tr_y_xxyy_xxxxyy[i] = tr_y_yy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xyy_xxxyy[i] * fe_0 + tr_y_xyy_xxxxyy[i] * pa_x[i];

        tr_y_xxyy_xxxxyz[i] = tr_y_yy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xyy_xxxyz[i] * fe_0 + tr_y_xyy_xxxxyz[i] * pa_x[i];

        tr_y_xxyy_xxxxzz[i] = tr_y_yy_xxxxzz[i] * fe_0 + 4.0 * tr_y_xyy_xxxzz[i] * fe_0 + tr_y_xyy_xxxxzz[i] * pa_x[i];

        tr_y_xxyy_xxxyyy[i] = tr_y_yy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xyy_xxyyy[i] * fe_0 + tr_y_xyy_xxxyyy[i] * pa_x[i];

        tr_y_xxyy_xxxyyz[i] = tr_y_yy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xyy_xxyyz[i] * fe_0 + tr_y_xyy_xxxyyz[i] * pa_x[i];

        tr_y_xxyy_xxxyzz[i] = tr_y_yy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xyy_xxyzz[i] * fe_0 + tr_y_xyy_xxxyzz[i] * pa_x[i];

        tr_y_xxyy_xxxzzz[i] = tr_y_yy_xxxzzz[i] * fe_0 + 3.0 * tr_y_xyy_xxzzz[i] * fe_0 + tr_y_xyy_xxxzzz[i] * pa_x[i];

        tr_y_xxyy_xxyyyy[i] = tr_y_yy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xyy_xyyyy[i] * fe_0 + tr_y_xyy_xxyyyy[i] * pa_x[i];

        tr_y_xxyy_xxyyyz[i] = tr_y_yy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xyy_xyyyz[i] * fe_0 + tr_y_xyy_xxyyyz[i] * pa_x[i];

        tr_y_xxyy_xxyyzz[i] = tr_y_yy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xyy_xyyzz[i] * fe_0 + tr_y_xyy_xxyyzz[i] * pa_x[i];

        tr_y_xxyy_xxyzzz[i] = tr_y_yy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xyy_xyzzz[i] * fe_0 + tr_y_xyy_xxyzzz[i] * pa_x[i];

        tr_y_xxyy_xxzzzz[i] = tr_y_yy_xxzzzz[i] * fe_0 + 2.0 * tr_y_xyy_xzzzz[i] * fe_0 + tr_y_xyy_xxzzzz[i] * pa_x[i];

        tr_y_xxyy_xyyyyy[i] = tr_y_yy_xyyyyy[i] * fe_0 + tr_y_xyy_yyyyy[i] * fe_0 + tr_y_xyy_xyyyyy[i] * pa_x[i];

        tr_y_xxyy_xyyyyz[i] = tr_y_yy_xyyyyz[i] * fe_0 + tr_y_xyy_yyyyz[i] * fe_0 + tr_y_xyy_xyyyyz[i] * pa_x[i];

        tr_y_xxyy_xyyyzz[i] = tr_y_yy_xyyyzz[i] * fe_0 + tr_y_xyy_yyyzz[i] * fe_0 + tr_y_xyy_xyyyzz[i] * pa_x[i];

        tr_y_xxyy_xyyzzz[i] = tr_y_yy_xyyzzz[i] * fe_0 + tr_y_xyy_yyzzz[i] * fe_0 + tr_y_xyy_xyyzzz[i] * pa_x[i];

        tr_y_xxyy_xyzzzz[i] = tr_y_yy_xyzzzz[i] * fe_0 + tr_y_xyy_yzzzz[i] * fe_0 + tr_y_xyy_xyzzzz[i] * pa_x[i];

        tr_y_xxyy_xzzzzz[i] = tr_y_yy_xzzzzz[i] * fe_0 + tr_y_xyy_zzzzz[i] * fe_0 + tr_y_xyy_xzzzzz[i] * pa_x[i];

        tr_y_xxyy_yyyyyy[i] = tr_y_yy_yyyyyy[i] * fe_0 + tr_y_xyy_yyyyyy[i] * pa_x[i];

        tr_y_xxyy_yyyyyz[i] = tr_y_yy_yyyyyz[i] * fe_0 + tr_y_xyy_yyyyyz[i] * pa_x[i];

        tr_y_xxyy_yyyyzz[i] = tr_y_yy_yyyyzz[i] * fe_0 + tr_y_xyy_yyyyzz[i] * pa_x[i];

        tr_y_xxyy_yyyzzz[i] = tr_y_yy_yyyzzz[i] * fe_0 + tr_y_xyy_yyyzzz[i] * pa_x[i];

        tr_y_xxyy_yyzzzz[i] = tr_y_yy_yyzzzz[i] * fe_0 + tr_y_xyy_yyzzzz[i] * pa_x[i];

        tr_y_xxyy_yzzzzz[i] = tr_y_yy_yzzzzz[i] * fe_0 + tr_y_xyy_yzzzzz[i] * pa_x[i];

        tr_y_xxyy_zzzzzz[i] = tr_y_yy_zzzzzz[i] * fe_0 + tr_y_xyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 532-560 components of targeted buffer : GI

    auto tr_y_xxyz_xxxxxx = pbuffer.data(idx_dip_gi + 532);

    auto tr_y_xxyz_xxxxxy = pbuffer.data(idx_dip_gi + 533);

    auto tr_y_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 534);

    auto tr_y_xxyz_xxxxyy = pbuffer.data(idx_dip_gi + 535);

    auto tr_y_xxyz_xxxxyz = pbuffer.data(idx_dip_gi + 536);

    auto tr_y_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 537);

    auto tr_y_xxyz_xxxyyy = pbuffer.data(idx_dip_gi + 538);

    auto tr_y_xxyz_xxxyyz = pbuffer.data(idx_dip_gi + 539);

    auto tr_y_xxyz_xxxyzz = pbuffer.data(idx_dip_gi + 540);

    auto tr_y_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 541);

    auto tr_y_xxyz_xxyyyy = pbuffer.data(idx_dip_gi + 542);

    auto tr_y_xxyz_xxyyyz = pbuffer.data(idx_dip_gi + 543);

    auto tr_y_xxyz_xxyyzz = pbuffer.data(idx_dip_gi + 544);

    auto tr_y_xxyz_xxyzzz = pbuffer.data(idx_dip_gi + 545);

    auto tr_y_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 546);

    auto tr_y_xxyz_xyyyyy = pbuffer.data(idx_dip_gi + 547);

    auto tr_y_xxyz_xyyyyz = pbuffer.data(idx_dip_gi + 548);

    auto tr_y_xxyz_xyyyzz = pbuffer.data(idx_dip_gi + 549);

    auto tr_y_xxyz_xyyzzz = pbuffer.data(idx_dip_gi + 550);

    auto tr_y_xxyz_xyzzzz = pbuffer.data(idx_dip_gi + 551);

    auto tr_y_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 552);

    auto tr_y_xxyz_yyyyyy = pbuffer.data(idx_dip_gi + 553);

    auto tr_y_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 554);

    auto tr_y_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 555);

    auto tr_y_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 556);

    auto tr_y_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 557);

    auto tr_y_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 558);

    auto tr_y_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 559);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxy_xxxxxx, tr_y_xxy_xxxxxy, tr_y_xxy_xxxxy, tr_y_xxy_xxxxyy, tr_y_xxy_xxxxyz, tr_y_xxy_xxxyy, tr_y_xxy_xxxyyy, tr_y_xxy_xxxyyz, tr_y_xxy_xxxyz, tr_y_xxy_xxxyzz, tr_y_xxy_xxyyy, tr_y_xxy_xxyyyy, tr_y_xxy_xxyyyz, tr_y_xxy_xxyyz, tr_y_xxy_xxyyzz, tr_y_xxy_xxyzz, tr_y_xxy_xxyzzz, tr_y_xxy_xyyyy, tr_y_xxy_xyyyyy, tr_y_xxy_xyyyyz, tr_y_xxy_xyyyz, tr_y_xxy_xyyyzz, tr_y_xxy_xyyzz, tr_y_xxy_xyyzzz, tr_y_xxy_xyzzz, tr_y_xxy_xyzzzz, tr_y_xxy_yyyyyy, tr_y_xxyz_xxxxxx, tr_y_xxyz_xxxxxy, tr_y_xxyz_xxxxxz, tr_y_xxyz_xxxxyy, tr_y_xxyz_xxxxyz, tr_y_xxyz_xxxxzz, tr_y_xxyz_xxxyyy, tr_y_xxyz_xxxyyz, tr_y_xxyz_xxxyzz, tr_y_xxyz_xxxzzz, tr_y_xxyz_xxyyyy, tr_y_xxyz_xxyyyz, tr_y_xxyz_xxyyzz, tr_y_xxyz_xxyzzz, tr_y_xxyz_xxzzzz, tr_y_xxyz_xyyyyy, tr_y_xxyz_xyyyyz, tr_y_xxyz_xyyyzz, tr_y_xxyz_xyyzzz, tr_y_xxyz_xyzzzz, tr_y_xxyz_xzzzzz, tr_y_xxyz_yyyyyy, tr_y_xxyz_yyyyyz, tr_y_xxyz_yyyyzz, tr_y_xxyz_yyyzzz, tr_y_xxyz_yyzzzz, tr_y_xxyz_yzzzzz, tr_y_xxyz_zzzzzz, tr_y_xxz_xxxxxz, tr_y_xxz_xxxxzz, tr_y_xxz_xxxzzz, tr_y_xxz_xxzzzz, tr_y_xxz_xzzzzz, tr_y_xyz_yyyyyz, tr_y_xyz_yyyyzz, tr_y_xyz_yyyzzz, tr_y_xyz_yyzzzz, tr_y_xyz_yzzzzz, tr_y_xyz_zzzzzz, tr_y_yz_yyyyyz, tr_y_yz_yyyyzz, tr_y_yz_yyyzzz, tr_y_yz_yyzzzz, tr_y_yz_yzzzzz, tr_y_yz_zzzzzz, ts_xxz_xxxxxz, ts_xxz_xxxxzz, ts_xxz_xxxzzz, ts_xxz_xxzzzz, ts_xxz_xzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyz_xxxxxx[i] = tr_y_xxy_xxxxxx[i] * pa_z[i];

        tr_y_xxyz_xxxxxy[i] = tr_y_xxy_xxxxxy[i] * pa_z[i];

        tr_y_xxyz_xxxxxz[i] = ts_xxz_xxxxxz[i] * fe_0 + tr_y_xxz_xxxxxz[i] * pa_y[i];

        tr_y_xxyz_xxxxyy[i] = tr_y_xxy_xxxxyy[i] * pa_z[i];

        tr_y_xxyz_xxxxyz[i] = tr_y_xxy_xxxxy[i] * fe_0 + tr_y_xxy_xxxxyz[i] * pa_z[i];

        tr_y_xxyz_xxxxzz[i] = ts_xxz_xxxxzz[i] * fe_0 + tr_y_xxz_xxxxzz[i] * pa_y[i];

        tr_y_xxyz_xxxyyy[i] = tr_y_xxy_xxxyyy[i] * pa_z[i];

        tr_y_xxyz_xxxyyz[i] = tr_y_xxy_xxxyy[i] * fe_0 + tr_y_xxy_xxxyyz[i] * pa_z[i];

        tr_y_xxyz_xxxyzz[i] = 2.0 * tr_y_xxy_xxxyz[i] * fe_0 + tr_y_xxy_xxxyzz[i] * pa_z[i];

        tr_y_xxyz_xxxzzz[i] = ts_xxz_xxxzzz[i] * fe_0 + tr_y_xxz_xxxzzz[i] * pa_y[i];

        tr_y_xxyz_xxyyyy[i] = tr_y_xxy_xxyyyy[i] * pa_z[i];

        tr_y_xxyz_xxyyyz[i] = tr_y_xxy_xxyyy[i] * fe_0 + tr_y_xxy_xxyyyz[i] * pa_z[i];

        tr_y_xxyz_xxyyzz[i] = 2.0 * tr_y_xxy_xxyyz[i] * fe_0 + tr_y_xxy_xxyyzz[i] * pa_z[i];

        tr_y_xxyz_xxyzzz[i] = 3.0 * tr_y_xxy_xxyzz[i] * fe_0 + tr_y_xxy_xxyzzz[i] * pa_z[i];

        tr_y_xxyz_xxzzzz[i] = ts_xxz_xxzzzz[i] * fe_0 + tr_y_xxz_xxzzzz[i] * pa_y[i];

        tr_y_xxyz_xyyyyy[i] = tr_y_xxy_xyyyyy[i] * pa_z[i];

        tr_y_xxyz_xyyyyz[i] = tr_y_xxy_xyyyy[i] * fe_0 + tr_y_xxy_xyyyyz[i] * pa_z[i];

        tr_y_xxyz_xyyyzz[i] = 2.0 * tr_y_xxy_xyyyz[i] * fe_0 + tr_y_xxy_xyyyzz[i] * pa_z[i];

        tr_y_xxyz_xyyzzz[i] = 3.0 * tr_y_xxy_xyyzz[i] * fe_0 + tr_y_xxy_xyyzzz[i] * pa_z[i];

        tr_y_xxyz_xyzzzz[i] = 4.0 * tr_y_xxy_xyzzz[i] * fe_0 + tr_y_xxy_xyzzzz[i] * pa_z[i];

        tr_y_xxyz_xzzzzz[i] = ts_xxz_xzzzzz[i] * fe_0 + tr_y_xxz_xzzzzz[i] * pa_y[i];

        tr_y_xxyz_yyyyyy[i] = tr_y_xxy_yyyyyy[i] * pa_z[i];

        tr_y_xxyz_yyyyyz[i] = tr_y_yz_yyyyyz[i] * fe_0 + tr_y_xyz_yyyyyz[i] * pa_x[i];

        tr_y_xxyz_yyyyzz[i] = tr_y_yz_yyyyzz[i] * fe_0 + tr_y_xyz_yyyyzz[i] * pa_x[i];

        tr_y_xxyz_yyyzzz[i] = tr_y_yz_yyyzzz[i] * fe_0 + tr_y_xyz_yyyzzz[i] * pa_x[i];

        tr_y_xxyz_yyzzzz[i] = tr_y_yz_yyzzzz[i] * fe_0 + tr_y_xyz_yyzzzz[i] * pa_x[i];

        tr_y_xxyz_yzzzzz[i] = tr_y_yz_yzzzzz[i] * fe_0 + tr_y_xyz_yzzzzz[i] * pa_x[i];

        tr_y_xxyz_zzzzzz[i] = tr_y_yz_zzzzzz[i] * fe_0 + tr_y_xyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 560-588 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xx_xxxxxx, tr_y_xx_xxxxxy, tr_y_xx_xxxxyy, tr_y_xx_xxxyyy, tr_y_xx_xxyyyy, tr_y_xx_xyyyyy, tr_y_xxz_xxxxxx, tr_y_xxz_xxxxxy, tr_y_xxz_xxxxyy, tr_y_xxz_xxxyyy, tr_y_xxz_xxyyyy, tr_y_xxz_xyyyyy, tr_y_xxzz_xxxxxx, tr_y_xxzz_xxxxxy, tr_y_xxzz_xxxxxz, tr_y_xxzz_xxxxyy, tr_y_xxzz_xxxxyz, tr_y_xxzz_xxxxzz, tr_y_xxzz_xxxyyy, tr_y_xxzz_xxxyyz, tr_y_xxzz_xxxyzz, tr_y_xxzz_xxxzzz, tr_y_xxzz_xxyyyy, tr_y_xxzz_xxyyyz, tr_y_xxzz_xxyyzz, tr_y_xxzz_xxyzzz, tr_y_xxzz_xxzzzz, tr_y_xxzz_xyyyyy, tr_y_xxzz_xyyyyz, tr_y_xxzz_xyyyzz, tr_y_xxzz_xyyzzz, tr_y_xxzz_xyzzzz, tr_y_xxzz_xzzzzz, tr_y_xxzz_yyyyyy, tr_y_xxzz_yyyyyz, tr_y_xxzz_yyyyzz, tr_y_xxzz_yyyzzz, tr_y_xxzz_yyzzzz, tr_y_xxzz_yzzzzz, tr_y_xxzz_zzzzzz, tr_y_xzz_xxxxxz, tr_y_xzz_xxxxyz, tr_y_xzz_xxxxz, tr_y_xzz_xxxxzz, tr_y_xzz_xxxyyz, tr_y_xzz_xxxyz, tr_y_xzz_xxxyzz, tr_y_xzz_xxxzz, tr_y_xzz_xxxzzz, tr_y_xzz_xxyyyz, tr_y_xzz_xxyyz, tr_y_xzz_xxyyzz, tr_y_xzz_xxyzz, tr_y_xzz_xxyzzz, tr_y_xzz_xxzzz, tr_y_xzz_xxzzzz, tr_y_xzz_xyyyyz, tr_y_xzz_xyyyz, tr_y_xzz_xyyyzz, tr_y_xzz_xyyzz, tr_y_xzz_xyyzzz, tr_y_xzz_xyzzz, tr_y_xzz_xyzzzz, tr_y_xzz_xzzzz, tr_y_xzz_xzzzzz, tr_y_xzz_yyyyyy, tr_y_xzz_yyyyyz, tr_y_xzz_yyyyz, tr_y_xzz_yyyyzz, tr_y_xzz_yyyzz, tr_y_xzz_yyyzzz, tr_y_xzz_yyzzz, tr_y_xzz_yyzzzz, tr_y_xzz_yzzzz, tr_y_xzz_yzzzzz, tr_y_xzz_zzzzz, tr_y_xzz_zzzzzz, tr_y_zz_xxxxxz, tr_y_zz_xxxxyz, tr_y_zz_xxxxzz, tr_y_zz_xxxyyz, tr_y_zz_xxxyzz, tr_y_zz_xxxzzz, tr_y_zz_xxyyyz, tr_y_zz_xxyyzz, tr_y_zz_xxyzzz, tr_y_zz_xxzzzz, tr_y_zz_xyyyyz, tr_y_zz_xyyyzz, tr_y_zz_xyyzzz, tr_y_zz_xyzzzz, tr_y_zz_xzzzzz, tr_y_zz_yyyyyy, tr_y_zz_yyyyyz, tr_y_zz_yyyyzz, tr_y_zz_yyyzzz, tr_y_zz_yyzzzz, tr_y_zz_yzzzzz, tr_y_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzz_xxxxxx[i] = tr_y_xx_xxxxxx[i] * fe_0 + tr_y_xxz_xxxxxx[i] * pa_z[i];

        tr_y_xxzz_xxxxxy[i] = tr_y_xx_xxxxxy[i] * fe_0 + tr_y_xxz_xxxxxy[i] * pa_z[i];

        tr_y_xxzz_xxxxxz[i] = tr_y_zz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xzz_xxxxz[i] * fe_0 + tr_y_xzz_xxxxxz[i] * pa_x[i];

        tr_y_xxzz_xxxxyy[i] = tr_y_xx_xxxxyy[i] * fe_0 + tr_y_xxz_xxxxyy[i] * pa_z[i];

        tr_y_xxzz_xxxxyz[i] = tr_y_zz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xzz_xxxyz[i] * fe_0 + tr_y_xzz_xxxxyz[i] * pa_x[i];

        tr_y_xxzz_xxxxzz[i] = tr_y_zz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xzz_xxxzz[i] * fe_0 + tr_y_xzz_xxxxzz[i] * pa_x[i];

        tr_y_xxzz_xxxyyy[i] = tr_y_xx_xxxyyy[i] * fe_0 + tr_y_xxz_xxxyyy[i] * pa_z[i];

        tr_y_xxzz_xxxyyz[i] = tr_y_zz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xzz_xxyyz[i] * fe_0 + tr_y_xzz_xxxyyz[i] * pa_x[i];

        tr_y_xxzz_xxxyzz[i] = tr_y_zz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xzz_xxyzz[i] * fe_0 + tr_y_xzz_xxxyzz[i] * pa_x[i];

        tr_y_xxzz_xxxzzz[i] = tr_y_zz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xzz_xxzzz[i] * fe_0 + tr_y_xzz_xxxzzz[i] * pa_x[i];

        tr_y_xxzz_xxyyyy[i] = tr_y_xx_xxyyyy[i] * fe_0 + tr_y_xxz_xxyyyy[i] * pa_z[i];

        tr_y_xxzz_xxyyyz[i] = tr_y_zz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xzz_xyyyz[i] * fe_0 + tr_y_xzz_xxyyyz[i] * pa_x[i];

        tr_y_xxzz_xxyyzz[i] = tr_y_zz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xzz_xyyzz[i] * fe_0 + tr_y_xzz_xxyyzz[i] * pa_x[i];

        tr_y_xxzz_xxyzzz[i] = tr_y_zz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xzz_xyzzz[i] * fe_0 + tr_y_xzz_xxyzzz[i] * pa_x[i];

        tr_y_xxzz_xxzzzz[i] = tr_y_zz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xzz_xzzzz[i] * fe_0 + tr_y_xzz_xxzzzz[i] * pa_x[i];

        tr_y_xxzz_xyyyyy[i] = tr_y_xx_xyyyyy[i] * fe_0 + tr_y_xxz_xyyyyy[i] * pa_z[i];

        tr_y_xxzz_xyyyyz[i] = tr_y_zz_xyyyyz[i] * fe_0 + tr_y_xzz_yyyyz[i] * fe_0 + tr_y_xzz_xyyyyz[i] * pa_x[i];

        tr_y_xxzz_xyyyzz[i] = tr_y_zz_xyyyzz[i] * fe_0 + tr_y_xzz_yyyzz[i] * fe_0 + tr_y_xzz_xyyyzz[i] * pa_x[i];

        tr_y_xxzz_xyyzzz[i] = tr_y_zz_xyyzzz[i] * fe_0 + tr_y_xzz_yyzzz[i] * fe_0 + tr_y_xzz_xyyzzz[i] * pa_x[i];

        tr_y_xxzz_xyzzzz[i] = tr_y_zz_xyzzzz[i] * fe_0 + tr_y_xzz_yzzzz[i] * fe_0 + tr_y_xzz_xyzzzz[i] * pa_x[i];

        tr_y_xxzz_xzzzzz[i] = tr_y_zz_xzzzzz[i] * fe_0 + tr_y_xzz_zzzzz[i] * fe_0 + tr_y_xzz_xzzzzz[i] * pa_x[i];

        tr_y_xxzz_yyyyyy[i] = tr_y_zz_yyyyyy[i] * fe_0 + tr_y_xzz_yyyyyy[i] * pa_x[i];

        tr_y_xxzz_yyyyyz[i] = tr_y_zz_yyyyyz[i] * fe_0 + tr_y_xzz_yyyyyz[i] * pa_x[i];

        tr_y_xxzz_yyyyzz[i] = tr_y_zz_yyyyzz[i] * fe_0 + tr_y_xzz_yyyyzz[i] * pa_x[i];

        tr_y_xxzz_yyyzzz[i] = tr_y_zz_yyyzzz[i] * fe_0 + tr_y_xzz_yyyzzz[i] * pa_x[i];

        tr_y_xxzz_yyzzzz[i] = tr_y_zz_yyzzzz[i] * fe_0 + tr_y_xzz_yyzzzz[i] * pa_x[i];

        tr_y_xxzz_yzzzzz[i] = tr_y_zz_yzzzzz[i] * fe_0 + tr_y_xzz_yzzzzz[i] * pa_x[i];

        tr_y_xxzz_zzzzzz[i] = tr_y_zz_zzzzzz[i] * fe_0 + tr_y_xzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 588-616 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_y_xyyy_xxxxxx, tr_y_xyyy_xxxxxy, tr_y_xyyy_xxxxxz, tr_y_xyyy_xxxxyy, tr_y_xyyy_xxxxyz, tr_y_xyyy_xxxxzz, tr_y_xyyy_xxxyyy, tr_y_xyyy_xxxyyz, tr_y_xyyy_xxxyzz, tr_y_xyyy_xxxzzz, tr_y_xyyy_xxyyyy, tr_y_xyyy_xxyyyz, tr_y_xyyy_xxyyzz, tr_y_xyyy_xxyzzz, tr_y_xyyy_xxzzzz, tr_y_xyyy_xyyyyy, tr_y_xyyy_xyyyyz, tr_y_xyyy_xyyyzz, tr_y_xyyy_xyyzzz, tr_y_xyyy_xyzzzz, tr_y_xyyy_xzzzzz, tr_y_xyyy_yyyyyy, tr_y_xyyy_yyyyyz, tr_y_xyyy_yyyyzz, tr_y_xyyy_yyyzzz, tr_y_xyyy_yyzzzz, tr_y_xyyy_yzzzzz, tr_y_xyyy_zzzzzz, tr_y_yyy_xxxxx, tr_y_yyy_xxxxxx, tr_y_yyy_xxxxxy, tr_y_yyy_xxxxxz, tr_y_yyy_xxxxy, tr_y_yyy_xxxxyy, tr_y_yyy_xxxxyz, tr_y_yyy_xxxxz, tr_y_yyy_xxxxzz, tr_y_yyy_xxxyy, tr_y_yyy_xxxyyy, tr_y_yyy_xxxyyz, tr_y_yyy_xxxyz, tr_y_yyy_xxxyzz, tr_y_yyy_xxxzz, tr_y_yyy_xxxzzz, tr_y_yyy_xxyyy, tr_y_yyy_xxyyyy, tr_y_yyy_xxyyyz, tr_y_yyy_xxyyz, tr_y_yyy_xxyyzz, tr_y_yyy_xxyzz, tr_y_yyy_xxyzzz, tr_y_yyy_xxzzz, tr_y_yyy_xxzzzz, tr_y_yyy_xyyyy, tr_y_yyy_xyyyyy, tr_y_yyy_xyyyyz, tr_y_yyy_xyyyz, tr_y_yyy_xyyyzz, tr_y_yyy_xyyzz, tr_y_yyy_xyyzzz, tr_y_yyy_xyzzz, tr_y_yyy_xyzzzz, tr_y_yyy_xzzzz, tr_y_yyy_xzzzzz, tr_y_yyy_yyyyy, tr_y_yyy_yyyyyy, tr_y_yyy_yyyyyz, tr_y_yyy_yyyyz, tr_y_yyy_yyyyzz, tr_y_yyy_yyyzz, tr_y_yyy_yyyzzz, tr_y_yyy_yyzzz, tr_y_yyy_yyzzzz, tr_y_yyy_yzzzz, tr_y_yyy_yzzzzz, tr_y_yyy_zzzzz, tr_y_yyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyy_xxxxxx[i] = 6.0 * tr_y_yyy_xxxxx[i] * fe_0 + tr_y_yyy_xxxxxx[i] * pa_x[i];

        tr_y_xyyy_xxxxxy[i] = 5.0 * tr_y_yyy_xxxxy[i] * fe_0 + tr_y_yyy_xxxxxy[i] * pa_x[i];

        tr_y_xyyy_xxxxxz[i] = 5.0 * tr_y_yyy_xxxxz[i] * fe_0 + tr_y_yyy_xxxxxz[i] * pa_x[i];

        tr_y_xyyy_xxxxyy[i] = 4.0 * tr_y_yyy_xxxyy[i] * fe_0 + tr_y_yyy_xxxxyy[i] * pa_x[i];

        tr_y_xyyy_xxxxyz[i] = 4.0 * tr_y_yyy_xxxyz[i] * fe_0 + tr_y_yyy_xxxxyz[i] * pa_x[i];

        tr_y_xyyy_xxxxzz[i] = 4.0 * tr_y_yyy_xxxzz[i] * fe_0 + tr_y_yyy_xxxxzz[i] * pa_x[i];

        tr_y_xyyy_xxxyyy[i] = 3.0 * tr_y_yyy_xxyyy[i] * fe_0 + tr_y_yyy_xxxyyy[i] * pa_x[i];

        tr_y_xyyy_xxxyyz[i] = 3.0 * tr_y_yyy_xxyyz[i] * fe_0 + tr_y_yyy_xxxyyz[i] * pa_x[i];

        tr_y_xyyy_xxxyzz[i] = 3.0 * tr_y_yyy_xxyzz[i] * fe_0 + tr_y_yyy_xxxyzz[i] * pa_x[i];

        tr_y_xyyy_xxxzzz[i] = 3.0 * tr_y_yyy_xxzzz[i] * fe_0 + tr_y_yyy_xxxzzz[i] * pa_x[i];

        tr_y_xyyy_xxyyyy[i] = 2.0 * tr_y_yyy_xyyyy[i] * fe_0 + tr_y_yyy_xxyyyy[i] * pa_x[i];

        tr_y_xyyy_xxyyyz[i] = 2.0 * tr_y_yyy_xyyyz[i] * fe_0 + tr_y_yyy_xxyyyz[i] * pa_x[i];

        tr_y_xyyy_xxyyzz[i] = 2.0 * tr_y_yyy_xyyzz[i] * fe_0 + tr_y_yyy_xxyyzz[i] * pa_x[i];

        tr_y_xyyy_xxyzzz[i] = 2.0 * tr_y_yyy_xyzzz[i] * fe_0 + tr_y_yyy_xxyzzz[i] * pa_x[i];

        tr_y_xyyy_xxzzzz[i] = 2.0 * tr_y_yyy_xzzzz[i] * fe_0 + tr_y_yyy_xxzzzz[i] * pa_x[i];

        tr_y_xyyy_xyyyyy[i] = tr_y_yyy_yyyyy[i] * fe_0 + tr_y_yyy_xyyyyy[i] * pa_x[i];

        tr_y_xyyy_xyyyyz[i] = tr_y_yyy_yyyyz[i] * fe_0 + tr_y_yyy_xyyyyz[i] * pa_x[i];

        tr_y_xyyy_xyyyzz[i] = tr_y_yyy_yyyzz[i] * fe_0 + tr_y_yyy_xyyyzz[i] * pa_x[i];

        tr_y_xyyy_xyyzzz[i] = tr_y_yyy_yyzzz[i] * fe_0 + tr_y_yyy_xyyzzz[i] * pa_x[i];

        tr_y_xyyy_xyzzzz[i] = tr_y_yyy_yzzzz[i] * fe_0 + tr_y_yyy_xyzzzz[i] * pa_x[i];

        tr_y_xyyy_xzzzzz[i] = tr_y_yyy_zzzzz[i] * fe_0 + tr_y_yyy_xzzzzz[i] * pa_x[i];

        tr_y_xyyy_yyyyyy[i] = tr_y_yyy_yyyyyy[i] * pa_x[i];

        tr_y_xyyy_yyyyyz[i] = tr_y_yyy_yyyyyz[i] * pa_x[i];

        tr_y_xyyy_yyyyzz[i] = tr_y_yyy_yyyyzz[i] * pa_x[i];

        tr_y_xyyy_yyyzzz[i] = tr_y_yyy_yyyzzz[i] * pa_x[i];

        tr_y_xyyy_yyzzzz[i] = tr_y_yyy_yyzzzz[i] * pa_x[i];

        tr_y_xyyy_yzzzzz[i] = tr_y_yyy_yzzzzz[i] * pa_x[i];

        tr_y_xyyy_zzzzzz[i] = tr_y_yyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 616-644 components of targeted buffer : GI

    auto tr_y_xyyz_xxxxxx = pbuffer.data(idx_dip_gi + 616);

    auto tr_y_xyyz_xxxxxy = pbuffer.data(idx_dip_gi + 617);

    auto tr_y_xyyz_xxxxxz = pbuffer.data(idx_dip_gi + 618);

    auto tr_y_xyyz_xxxxyy = pbuffer.data(idx_dip_gi + 619);

    auto tr_y_xyyz_xxxxyz = pbuffer.data(idx_dip_gi + 620);

    auto tr_y_xyyz_xxxxzz = pbuffer.data(idx_dip_gi + 621);

    auto tr_y_xyyz_xxxyyy = pbuffer.data(idx_dip_gi + 622);

    auto tr_y_xyyz_xxxyyz = pbuffer.data(idx_dip_gi + 623);

    auto tr_y_xyyz_xxxyzz = pbuffer.data(idx_dip_gi + 624);

    auto tr_y_xyyz_xxxzzz = pbuffer.data(idx_dip_gi + 625);

    auto tr_y_xyyz_xxyyyy = pbuffer.data(idx_dip_gi + 626);

    auto tr_y_xyyz_xxyyyz = pbuffer.data(idx_dip_gi + 627);

    auto tr_y_xyyz_xxyyzz = pbuffer.data(idx_dip_gi + 628);

    auto tr_y_xyyz_xxyzzz = pbuffer.data(idx_dip_gi + 629);

    auto tr_y_xyyz_xxzzzz = pbuffer.data(idx_dip_gi + 630);

    auto tr_y_xyyz_xyyyyy = pbuffer.data(idx_dip_gi + 631);

    auto tr_y_xyyz_xyyyyz = pbuffer.data(idx_dip_gi + 632);

    auto tr_y_xyyz_xyyyzz = pbuffer.data(idx_dip_gi + 633);

    auto tr_y_xyyz_xyyzzz = pbuffer.data(idx_dip_gi + 634);

    auto tr_y_xyyz_xyzzzz = pbuffer.data(idx_dip_gi + 635);

    auto tr_y_xyyz_xzzzzz = pbuffer.data(idx_dip_gi + 636);

    auto tr_y_xyyz_yyyyyy = pbuffer.data(idx_dip_gi + 637);

    auto tr_y_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 638);

    auto tr_y_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 639);

    auto tr_y_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 640);

    auto tr_y_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 641);

    auto tr_y_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 642);

    auto tr_y_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 643);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xyy_xxxxxx, tr_y_xyy_xxxxxy, tr_y_xyy_xxxxyy, tr_y_xyy_xxxyyy, tr_y_xyy_xxyyyy, tr_y_xyy_xyyyyy, tr_y_xyyz_xxxxxx, tr_y_xyyz_xxxxxy, tr_y_xyyz_xxxxxz, tr_y_xyyz_xxxxyy, tr_y_xyyz_xxxxyz, tr_y_xyyz_xxxxzz, tr_y_xyyz_xxxyyy, tr_y_xyyz_xxxyyz, tr_y_xyyz_xxxyzz, tr_y_xyyz_xxxzzz, tr_y_xyyz_xxyyyy, tr_y_xyyz_xxyyyz, tr_y_xyyz_xxyyzz, tr_y_xyyz_xxyzzz, tr_y_xyyz_xxzzzz, tr_y_xyyz_xyyyyy, tr_y_xyyz_xyyyyz, tr_y_xyyz_xyyyzz, tr_y_xyyz_xyyzzz, tr_y_xyyz_xyzzzz, tr_y_xyyz_xzzzzz, tr_y_xyyz_yyyyyy, tr_y_xyyz_yyyyyz, tr_y_xyyz_yyyyzz, tr_y_xyyz_yyyzzz, tr_y_xyyz_yyzzzz, tr_y_xyyz_yzzzzz, tr_y_xyyz_zzzzzz, tr_y_yyz_xxxxxz, tr_y_yyz_xxxxyz, tr_y_yyz_xxxxz, tr_y_yyz_xxxxzz, tr_y_yyz_xxxyyz, tr_y_yyz_xxxyz, tr_y_yyz_xxxyzz, tr_y_yyz_xxxzz, tr_y_yyz_xxxzzz, tr_y_yyz_xxyyyz, tr_y_yyz_xxyyz, tr_y_yyz_xxyyzz, tr_y_yyz_xxyzz, tr_y_yyz_xxyzzz, tr_y_yyz_xxzzz, tr_y_yyz_xxzzzz, tr_y_yyz_xyyyyz, tr_y_yyz_xyyyz, tr_y_yyz_xyyyzz, tr_y_yyz_xyyzz, tr_y_yyz_xyyzzz, tr_y_yyz_xyzzz, tr_y_yyz_xyzzzz, tr_y_yyz_xzzzz, tr_y_yyz_xzzzzz, tr_y_yyz_yyyyyy, tr_y_yyz_yyyyyz, tr_y_yyz_yyyyz, tr_y_yyz_yyyyzz, tr_y_yyz_yyyzz, tr_y_yyz_yyyzzz, tr_y_yyz_yyzzz, tr_y_yyz_yyzzzz, tr_y_yyz_yzzzz, tr_y_yyz_yzzzzz, tr_y_yyz_zzzzz, tr_y_yyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyz_xxxxxx[i] = tr_y_xyy_xxxxxx[i] * pa_z[i];

        tr_y_xyyz_xxxxxy[i] = tr_y_xyy_xxxxxy[i] * pa_z[i];

        tr_y_xyyz_xxxxxz[i] = 5.0 * tr_y_yyz_xxxxz[i] * fe_0 + tr_y_yyz_xxxxxz[i] * pa_x[i];

        tr_y_xyyz_xxxxyy[i] = tr_y_xyy_xxxxyy[i] * pa_z[i];

        tr_y_xyyz_xxxxyz[i] = 4.0 * tr_y_yyz_xxxyz[i] * fe_0 + tr_y_yyz_xxxxyz[i] * pa_x[i];

        tr_y_xyyz_xxxxzz[i] = 4.0 * tr_y_yyz_xxxzz[i] * fe_0 + tr_y_yyz_xxxxzz[i] * pa_x[i];

        tr_y_xyyz_xxxyyy[i] = tr_y_xyy_xxxyyy[i] * pa_z[i];

        tr_y_xyyz_xxxyyz[i] = 3.0 * tr_y_yyz_xxyyz[i] * fe_0 + tr_y_yyz_xxxyyz[i] * pa_x[i];

        tr_y_xyyz_xxxyzz[i] = 3.0 * tr_y_yyz_xxyzz[i] * fe_0 + tr_y_yyz_xxxyzz[i] * pa_x[i];

        tr_y_xyyz_xxxzzz[i] = 3.0 * tr_y_yyz_xxzzz[i] * fe_0 + tr_y_yyz_xxxzzz[i] * pa_x[i];

        tr_y_xyyz_xxyyyy[i] = tr_y_xyy_xxyyyy[i] * pa_z[i];

        tr_y_xyyz_xxyyyz[i] = 2.0 * tr_y_yyz_xyyyz[i] * fe_0 + tr_y_yyz_xxyyyz[i] * pa_x[i];

        tr_y_xyyz_xxyyzz[i] = 2.0 * tr_y_yyz_xyyzz[i] * fe_0 + tr_y_yyz_xxyyzz[i] * pa_x[i];

        tr_y_xyyz_xxyzzz[i] = 2.0 * tr_y_yyz_xyzzz[i] * fe_0 + tr_y_yyz_xxyzzz[i] * pa_x[i];

        tr_y_xyyz_xxzzzz[i] = 2.0 * tr_y_yyz_xzzzz[i] * fe_0 + tr_y_yyz_xxzzzz[i] * pa_x[i];

        tr_y_xyyz_xyyyyy[i] = tr_y_xyy_xyyyyy[i] * pa_z[i];

        tr_y_xyyz_xyyyyz[i] = tr_y_yyz_yyyyz[i] * fe_0 + tr_y_yyz_xyyyyz[i] * pa_x[i];

        tr_y_xyyz_xyyyzz[i] = tr_y_yyz_yyyzz[i] * fe_0 + tr_y_yyz_xyyyzz[i] * pa_x[i];

        tr_y_xyyz_xyyzzz[i] = tr_y_yyz_yyzzz[i] * fe_0 + tr_y_yyz_xyyzzz[i] * pa_x[i];

        tr_y_xyyz_xyzzzz[i] = tr_y_yyz_yzzzz[i] * fe_0 + tr_y_yyz_xyzzzz[i] * pa_x[i];

        tr_y_xyyz_xzzzzz[i] = tr_y_yyz_zzzzz[i] * fe_0 + tr_y_yyz_xzzzzz[i] * pa_x[i];

        tr_y_xyyz_yyyyyy[i] = tr_y_yyz_yyyyyy[i] * pa_x[i];

        tr_y_xyyz_yyyyyz[i] = tr_y_yyz_yyyyyz[i] * pa_x[i];

        tr_y_xyyz_yyyyzz[i] = tr_y_yyz_yyyyzz[i] * pa_x[i];

        tr_y_xyyz_yyyzzz[i] = tr_y_yyz_yyyzzz[i] * pa_x[i];

        tr_y_xyyz_yyzzzz[i] = tr_y_yyz_yyzzzz[i] * pa_x[i];

        tr_y_xyyz_yzzzzz[i] = tr_y_yyz_yzzzzz[i] * pa_x[i];

        tr_y_xyyz_zzzzzz[i] = tr_y_yyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 644-672 components of targeted buffer : GI

    auto tr_y_xyzz_xxxxxx = pbuffer.data(idx_dip_gi + 644);

    auto tr_y_xyzz_xxxxxy = pbuffer.data(idx_dip_gi + 645);

    auto tr_y_xyzz_xxxxxz = pbuffer.data(idx_dip_gi + 646);

    auto tr_y_xyzz_xxxxyy = pbuffer.data(idx_dip_gi + 647);

    auto tr_y_xyzz_xxxxyz = pbuffer.data(idx_dip_gi + 648);

    auto tr_y_xyzz_xxxxzz = pbuffer.data(idx_dip_gi + 649);

    auto tr_y_xyzz_xxxyyy = pbuffer.data(idx_dip_gi + 650);

    auto tr_y_xyzz_xxxyyz = pbuffer.data(idx_dip_gi + 651);

    auto tr_y_xyzz_xxxyzz = pbuffer.data(idx_dip_gi + 652);

    auto tr_y_xyzz_xxxzzz = pbuffer.data(idx_dip_gi + 653);

    auto tr_y_xyzz_xxyyyy = pbuffer.data(idx_dip_gi + 654);

    auto tr_y_xyzz_xxyyyz = pbuffer.data(idx_dip_gi + 655);

    auto tr_y_xyzz_xxyyzz = pbuffer.data(idx_dip_gi + 656);

    auto tr_y_xyzz_xxyzzz = pbuffer.data(idx_dip_gi + 657);

    auto tr_y_xyzz_xxzzzz = pbuffer.data(idx_dip_gi + 658);

    auto tr_y_xyzz_xyyyyy = pbuffer.data(idx_dip_gi + 659);

    auto tr_y_xyzz_xyyyyz = pbuffer.data(idx_dip_gi + 660);

    auto tr_y_xyzz_xyyyzz = pbuffer.data(idx_dip_gi + 661);

    auto tr_y_xyzz_xyyzzz = pbuffer.data(idx_dip_gi + 662);

    auto tr_y_xyzz_xyzzzz = pbuffer.data(idx_dip_gi + 663);

    auto tr_y_xyzz_xzzzzz = pbuffer.data(idx_dip_gi + 664);

    auto tr_y_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 665);

    auto tr_y_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 666);

    auto tr_y_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 667);

    auto tr_y_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 668);

    auto tr_y_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 669);

    auto tr_y_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 670);

    auto tr_y_xyzz_zzzzzz = pbuffer.data(idx_dip_gi + 671);

    #pragma omp simd aligned(pa_x, tr_y_xyzz_xxxxxx, tr_y_xyzz_xxxxxy, tr_y_xyzz_xxxxxz, tr_y_xyzz_xxxxyy, tr_y_xyzz_xxxxyz, tr_y_xyzz_xxxxzz, tr_y_xyzz_xxxyyy, tr_y_xyzz_xxxyyz, tr_y_xyzz_xxxyzz, tr_y_xyzz_xxxzzz, tr_y_xyzz_xxyyyy, tr_y_xyzz_xxyyyz, tr_y_xyzz_xxyyzz, tr_y_xyzz_xxyzzz, tr_y_xyzz_xxzzzz, tr_y_xyzz_xyyyyy, tr_y_xyzz_xyyyyz, tr_y_xyzz_xyyyzz, tr_y_xyzz_xyyzzz, tr_y_xyzz_xyzzzz, tr_y_xyzz_xzzzzz, tr_y_xyzz_yyyyyy, tr_y_xyzz_yyyyyz, tr_y_xyzz_yyyyzz, tr_y_xyzz_yyyzzz, tr_y_xyzz_yyzzzz, tr_y_xyzz_yzzzzz, tr_y_xyzz_zzzzzz, tr_y_yzz_xxxxx, tr_y_yzz_xxxxxx, tr_y_yzz_xxxxxy, tr_y_yzz_xxxxxz, tr_y_yzz_xxxxy, tr_y_yzz_xxxxyy, tr_y_yzz_xxxxyz, tr_y_yzz_xxxxz, tr_y_yzz_xxxxzz, tr_y_yzz_xxxyy, tr_y_yzz_xxxyyy, tr_y_yzz_xxxyyz, tr_y_yzz_xxxyz, tr_y_yzz_xxxyzz, tr_y_yzz_xxxzz, tr_y_yzz_xxxzzz, tr_y_yzz_xxyyy, tr_y_yzz_xxyyyy, tr_y_yzz_xxyyyz, tr_y_yzz_xxyyz, tr_y_yzz_xxyyzz, tr_y_yzz_xxyzz, tr_y_yzz_xxyzzz, tr_y_yzz_xxzzz, tr_y_yzz_xxzzzz, tr_y_yzz_xyyyy, tr_y_yzz_xyyyyy, tr_y_yzz_xyyyyz, tr_y_yzz_xyyyz, tr_y_yzz_xyyyzz, tr_y_yzz_xyyzz, tr_y_yzz_xyyzzz, tr_y_yzz_xyzzz, tr_y_yzz_xyzzzz, tr_y_yzz_xzzzz, tr_y_yzz_xzzzzz, tr_y_yzz_yyyyy, tr_y_yzz_yyyyyy, tr_y_yzz_yyyyyz, tr_y_yzz_yyyyz, tr_y_yzz_yyyyzz, tr_y_yzz_yyyzz, tr_y_yzz_yyyzzz, tr_y_yzz_yyzzz, tr_y_yzz_yyzzzz, tr_y_yzz_yzzzz, tr_y_yzz_yzzzzz, tr_y_yzz_zzzzz, tr_y_yzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzz_xxxxxx[i] = 6.0 * tr_y_yzz_xxxxx[i] * fe_0 + tr_y_yzz_xxxxxx[i] * pa_x[i];

        tr_y_xyzz_xxxxxy[i] = 5.0 * tr_y_yzz_xxxxy[i] * fe_0 + tr_y_yzz_xxxxxy[i] * pa_x[i];

        tr_y_xyzz_xxxxxz[i] = 5.0 * tr_y_yzz_xxxxz[i] * fe_0 + tr_y_yzz_xxxxxz[i] * pa_x[i];

        tr_y_xyzz_xxxxyy[i] = 4.0 * tr_y_yzz_xxxyy[i] * fe_0 + tr_y_yzz_xxxxyy[i] * pa_x[i];

        tr_y_xyzz_xxxxyz[i] = 4.0 * tr_y_yzz_xxxyz[i] * fe_0 + tr_y_yzz_xxxxyz[i] * pa_x[i];

        tr_y_xyzz_xxxxzz[i] = 4.0 * tr_y_yzz_xxxzz[i] * fe_0 + tr_y_yzz_xxxxzz[i] * pa_x[i];

        tr_y_xyzz_xxxyyy[i] = 3.0 * tr_y_yzz_xxyyy[i] * fe_0 + tr_y_yzz_xxxyyy[i] * pa_x[i];

        tr_y_xyzz_xxxyyz[i] = 3.0 * tr_y_yzz_xxyyz[i] * fe_0 + tr_y_yzz_xxxyyz[i] * pa_x[i];

        tr_y_xyzz_xxxyzz[i] = 3.0 * tr_y_yzz_xxyzz[i] * fe_0 + tr_y_yzz_xxxyzz[i] * pa_x[i];

        tr_y_xyzz_xxxzzz[i] = 3.0 * tr_y_yzz_xxzzz[i] * fe_0 + tr_y_yzz_xxxzzz[i] * pa_x[i];

        tr_y_xyzz_xxyyyy[i] = 2.0 * tr_y_yzz_xyyyy[i] * fe_0 + tr_y_yzz_xxyyyy[i] * pa_x[i];

        tr_y_xyzz_xxyyyz[i] = 2.0 * tr_y_yzz_xyyyz[i] * fe_0 + tr_y_yzz_xxyyyz[i] * pa_x[i];

        tr_y_xyzz_xxyyzz[i] = 2.0 * tr_y_yzz_xyyzz[i] * fe_0 + tr_y_yzz_xxyyzz[i] * pa_x[i];

        tr_y_xyzz_xxyzzz[i] = 2.0 * tr_y_yzz_xyzzz[i] * fe_0 + tr_y_yzz_xxyzzz[i] * pa_x[i];

        tr_y_xyzz_xxzzzz[i] = 2.0 * tr_y_yzz_xzzzz[i] * fe_0 + tr_y_yzz_xxzzzz[i] * pa_x[i];

        tr_y_xyzz_xyyyyy[i] = tr_y_yzz_yyyyy[i] * fe_0 + tr_y_yzz_xyyyyy[i] * pa_x[i];

        tr_y_xyzz_xyyyyz[i] = tr_y_yzz_yyyyz[i] * fe_0 + tr_y_yzz_xyyyyz[i] * pa_x[i];

        tr_y_xyzz_xyyyzz[i] = tr_y_yzz_yyyzz[i] * fe_0 + tr_y_yzz_xyyyzz[i] * pa_x[i];

        tr_y_xyzz_xyyzzz[i] = tr_y_yzz_yyzzz[i] * fe_0 + tr_y_yzz_xyyzzz[i] * pa_x[i];

        tr_y_xyzz_xyzzzz[i] = tr_y_yzz_yzzzz[i] * fe_0 + tr_y_yzz_xyzzzz[i] * pa_x[i];

        tr_y_xyzz_xzzzzz[i] = tr_y_yzz_zzzzz[i] * fe_0 + tr_y_yzz_xzzzzz[i] * pa_x[i];

        tr_y_xyzz_yyyyyy[i] = tr_y_yzz_yyyyyy[i] * pa_x[i];

        tr_y_xyzz_yyyyyz[i] = tr_y_yzz_yyyyyz[i] * pa_x[i];

        tr_y_xyzz_yyyyzz[i] = tr_y_yzz_yyyyzz[i] * pa_x[i];

        tr_y_xyzz_yyyzzz[i] = tr_y_yzz_yyyzzz[i] * pa_x[i];

        tr_y_xyzz_yyzzzz[i] = tr_y_yzz_yyzzzz[i] * pa_x[i];

        tr_y_xyzz_yzzzzz[i] = tr_y_yzz_yzzzzz[i] * pa_x[i];

        tr_y_xyzz_zzzzzz[i] = tr_y_yzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 672-700 components of targeted buffer : GI

    auto tr_y_xzzz_xxxxxx = pbuffer.data(idx_dip_gi + 672);

    auto tr_y_xzzz_xxxxxy = pbuffer.data(idx_dip_gi + 673);

    auto tr_y_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 674);

    auto tr_y_xzzz_xxxxyy = pbuffer.data(idx_dip_gi + 675);

    auto tr_y_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 676);

    auto tr_y_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 677);

    auto tr_y_xzzz_xxxyyy = pbuffer.data(idx_dip_gi + 678);

    auto tr_y_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 679);

    auto tr_y_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 680);

    auto tr_y_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 681);

    auto tr_y_xzzz_xxyyyy = pbuffer.data(idx_dip_gi + 682);

    auto tr_y_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 683);

    auto tr_y_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 684);

    auto tr_y_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 685);

    auto tr_y_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 686);

    auto tr_y_xzzz_xyyyyy = pbuffer.data(idx_dip_gi + 687);

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

    #pragma omp simd aligned(pa_x, tr_y_xzzz_xxxxxx, tr_y_xzzz_xxxxxy, tr_y_xzzz_xxxxxz, tr_y_xzzz_xxxxyy, tr_y_xzzz_xxxxyz, tr_y_xzzz_xxxxzz, tr_y_xzzz_xxxyyy, tr_y_xzzz_xxxyyz, tr_y_xzzz_xxxyzz, tr_y_xzzz_xxxzzz, tr_y_xzzz_xxyyyy, tr_y_xzzz_xxyyyz, tr_y_xzzz_xxyyzz, tr_y_xzzz_xxyzzz, tr_y_xzzz_xxzzzz, tr_y_xzzz_xyyyyy, tr_y_xzzz_xyyyyz, tr_y_xzzz_xyyyzz, tr_y_xzzz_xyyzzz, tr_y_xzzz_xyzzzz, tr_y_xzzz_xzzzzz, tr_y_xzzz_yyyyyy, tr_y_xzzz_yyyyyz, tr_y_xzzz_yyyyzz, tr_y_xzzz_yyyzzz, tr_y_xzzz_yyzzzz, tr_y_xzzz_yzzzzz, tr_y_xzzz_zzzzzz, tr_y_zzz_xxxxx, tr_y_zzz_xxxxxx, tr_y_zzz_xxxxxy, tr_y_zzz_xxxxxz, tr_y_zzz_xxxxy, tr_y_zzz_xxxxyy, tr_y_zzz_xxxxyz, tr_y_zzz_xxxxz, tr_y_zzz_xxxxzz, tr_y_zzz_xxxyy, tr_y_zzz_xxxyyy, tr_y_zzz_xxxyyz, tr_y_zzz_xxxyz, tr_y_zzz_xxxyzz, tr_y_zzz_xxxzz, tr_y_zzz_xxxzzz, tr_y_zzz_xxyyy, tr_y_zzz_xxyyyy, tr_y_zzz_xxyyyz, tr_y_zzz_xxyyz, tr_y_zzz_xxyyzz, tr_y_zzz_xxyzz, tr_y_zzz_xxyzzz, tr_y_zzz_xxzzz, tr_y_zzz_xxzzzz, tr_y_zzz_xyyyy, tr_y_zzz_xyyyyy, tr_y_zzz_xyyyyz, tr_y_zzz_xyyyz, tr_y_zzz_xyyyzz, tr_y_zzz_xyyzz, tr_y_zzz_xyyzzz, tr_y_zzz_xyzzz, tr_y_zzz_xyzzzz, tr_y_zzz_xzzzz, tr_y_zzz_xzzzzz, tr_y_zzz_yyyyy, tr_y_zzz_yyyyyy, tr_y_zzz_yyyyyz, tr_y_zzz_yyyyz, tr_y_zzz_yyyyzz, tr_y_zzz_yyyzz, tr_y_zzz_yyyzzz, tr_y_zzz_yyzzz, tr_y_zzz_yyzzzz, tr_y_zzz_yzzzz, tr_y_zzz_yzzzzz, tr_y_zzz_zzzzz, tr_y_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzz_xxxxxx[i] = 6.0 * tr_y_zzz_xxxxx[i] * fe_0 + tr_y_zzz_xxxxxx[i] * pa_x[i];

        tr_y_xzzz_xxxxxy[i] = 5.0 * tr_y_zzz_xxxxy[i] * fe_0 + tr_y_zzz_xxxxxy[i] * pa_x[i];

        tr_y_xzzz_xxxxxz[i] = 5.0 * tr_y_zzz_xxxxz[i] * fe_0 + tr_y_zzz_xxxxxz[i] * pa_x[i];

        tr_y_xzzz_xxxxyy[i] = 4.0 * tr_y_zzz_xxxyy[i] * fe_0 + tr_y_zzz_xxxxyy[i] * pa_x[i];

        tr_y_xzzz_xxxxyz[i] = 4.0 * tr_y_zzz_xxxyz[i] * fe_0 + tr_y_zzz_xxxxyz[i] * pa_x[i];

        tr_y_xzzz_xxxxzz[i] = 4.0 * tr_y_zzz_xxxzz[i] * fe_0 + tr_y_zzz_xxxxzz[i] * pa_x[i];

        tr_y_xzzz_xxxyyy[i] = 3.0 * tr_y_zzz_xxyyy[i] * fe_0 + tr_y_zzz_xxxyyy[i] * pa_x[i];

        tr_y_xzzz_xxxyyz[i] = 3.0 * tr_y_zzz_xxyyz[i] * fe_0 + tr_y_zzz_xxxyyz[i] * pa_x[i];

        tr_y_xzzz_xxxyzz[i] = 3.0 * tr_y_zzz_xxyzz[i] * fe_0 + tr_y_zzz_xxxyzz[i] * pa_x[i];

        tr_y_xzzz_xxxzzz[i] = 3.0 * tr_y_zzz_xxzzz[i] * fe_0 + tr_y_zzz_xxxzzz[i] * pa_x[i];

        tr_y_xzzz_xxyyyy[i] = 2.0 * tr_y_zzz_xyyyy[i] * fe_0 + tr_y_zzz_xxyyyy[i] * pa_x[i];

        tr_y_xzzz_xxyyyz[i] = 2.0 * tr_y_zzz_xyyyz[i] * fe_0 + tr_y_zzz_xxyyyz[i] * pa_x[i];

        tr_y_xzzz_xxyyzz[i] = 2.0 * tr_y_zzz_xyyzz[i] * fe_0 + tr_y_zzz_xxyyzz[i] * pa_x[i];

        tr_y_xzzz_xxyzzz[i] = 2.0 * tr_y_zzz_xyzzz[i] * fe_0 + tr_y_zzz_xxyzzz[i] * pa_x[i];

        tr_y_xzzz_xxzzzz[i] = 2.0 * tr_y_zzz_xzzzz[i] * fe_0 + tr_y_zzz_xxzzzz[i] * pa_x[i];

        tr_y_xzzz_xyyyyy[i] = tr_y_zzz_yyyyy[i] * fe_0 + tr_y_zzz_xyyyyy[i] * pa_x[i];

        tr_y_xzzz_xyyyyz[i] = tr_y_zzz_yyyyz[i] * fe_0 + tr_y_zzz_xyyyyz[i] * pa_x[i];

        tr_y_xzzz_xyyyzz[i] = tr_y_zzz_yyyzz[i] * fe_0 + tr_y_zzz_xyyyzz[i] * pa_x[i];

        tr_y_xzzz_xyyzzz[i] = tr_y_zzz_yyzzz[i] * fe_0 + tr_y_zzz_xyyzzz[i] * pa_x[i];

        tr_y_xzzz_xyzzzz[i] = tr_y_zzz_yzzzz[i] * fe_0 + tr_y_zzz_xyzzzz[i] * pa_x[i];

        tr_y_xzzz_xzzzzz[i] = tr_y_zzz_zzzzz[i] * fe_0 + tr_y_zzz_xzzzzz[i] * pa_x[i];

        tr_y_xzzz_yyyyyy[i] = tr_y_zzz_yyyyyy[i] * pa_x[i];

        tr_y_xzzz_yyyyyz[i] = tr_y_zzz_yyyyyz[i] * pa_x[i];

        tr_y_xzzz_yyyyzz[i] = tr_y_zzz_yyyyzz[i] * pa_x[i];

        tr_y_xzzz_yyyzzz[i] = tr_y_zzz_yyyzzz[i] * pa_x[i];

        tr_y_xzzz_yyzzzz[i] = tr_y_zzz_yyzzzz[i] * pa_x[i];

        tr_y_xzzz_yzzzzz[i] = tr_y_zzz_yzzzzz[i] * pa_x[i];

        tr_y_xzzz_zzzzzz[i] = tr_y_zzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 700-728 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, tr_y_yy_xxxxxx, tr_y_yy_xxxxxy, tr_y_yy_xxxxxz, tr_y_yy_xxxxyy, tr_y_yy_xxxxyz, tr_y_yy_xxxxzz, tr_y_yy_xxxyyy, tr_y_yy_xxxyyz, tr_y_yy_xxxyzz, tr_y_yy_xxxzzz, tr_y_yy_xxyyyy, tr_y_yy_xxyyyz, tr_y_yy_xxyyzz, tr_y_yy_xxyzzz, tr_y_yy_xxzzzz, tr_y_yy_xyyyyy, tr_y_yy_xyyyyz, tr_y_yy_xyyyzz, tr_y_yy_xyyzzz, tr_y_yy_xyzzzz, tr_y_yy_xzzzzz, tr_y_yy_yyyyyy, tr_y_yy_yyyyyz, tr_y_yy_yyyyzz, tr_y_yy_yyyzzz, tr_y_yy_yyzzzz, tr_y_yy_yzzzzz, tr_y_yy_zzzzzz, tr_y_yyy_xxxxx, tr_y_yyy_xxxxxx, tr_y_yyy_xxxxxy, tr_y_yyy_xxxxxz, tr_y_yyy_xxxxy, tr_y_yyy_xxxxyy, tr_y_yyy_xxxxyz, tr_y_yyy_xxxxz, tr_y_yyy_xxxxzz, tr_y_yyy_xxxyy, tr_y_yyy_xxxyyy, tr_y_yyy_xxxyyz, tr_y_yyy_xxxyz, tr_y_yyy_xxxyzz, tr_y_yyy_xxxzz, tr_y_yyy_xxxzzz, tr_y_yyy_xxyyy, tr_y_yyy_xxyyyy, tr_y_yyy_xxyyyz, tr_y_yyy_xxyyz, tr_y_yyy_xxyyzz, tr_y_yyy_xxyzz, tr_y_yyy_xxyzzz, tr_y_yyy_xxzzz, tr_y_yyy_xxzzzz, tr_y_yyy_xyyyy, tr_y_yyy_xyyyyy, tr_y_yyy_xyyyyz, tr_y_yyy_xyyyz, tr_y_yyy_xyyyzz, tr_y_yyy_xyyzz, tr_y_yyy_xyyzzz, tr_y_yyy_xyzzz, tr_y_yyy_xyzzzz, tr_y_yyy_xzzzz, tr_y_yyy_xzzzzz, tr_y_yyy_yyyyy, tr_y_yyy_yyyyyy, tr_y_yyy_yyyyyz, tr_y_yyy_yyyyz, tr_y_yyy_yyyyzz, tr_y_yyy_yyyzz, tr_y_yyy_yyyzzz, tr_y_yyy_yyzzz, tr_y_yyy_yyzzzz, tr_y_yyy_yzzzz, tr_y_yyy_yzzzzz, tr_y_yyy_zzzzz, tr_y_yyy_zzzzzz, tr_y_yyyy_xxxxxx, tr_y_yyyy_xxxxxy, tr_y_yyyy_xxxxxz, tr_y_yyyy_xxxxyy, tr_y_yyyy_xxxxyz, tr_y_yyyy_xxxxzz, tr_y_yyyy_xxxyyy, tr_y_yyyy_xxxyyz, tr_y_yyyy_xxxyzz, tr_y_yyyy_xxxzzz, tr_y_yyyy_xxyyyy, tr_y_yyyy_xxyyyz, tr_y_yyyy_xxyyzz, tr_y_yyyy_xxyzzz, tr_y_yyyy_xxzzzz, tr_y_yyyy_xyyyyy, tr_y_yyyy_xyyyyz, tr_y_yyyy_xyyyzz, tr_y_yyyy_xyyzzz, tr_y_yyyy_xyzzzz, tr_y_yyyy_xzzzzz, tr_y_yyyy_yyyyyy, tr_y_yyyy_yyyyyz, tr_y_yyyy_yyyyzz, tr_y_yyyy_yyyzzz, tr_y_yyyy_yyzzzz, tr_y_yyyy_yzzzzz, tr_y_yyyy_zzzzzz, ts_yyy_xxxxxx, ts_yyy_xxxxxy, ts_yyy_xxxxxz, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxxzz, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyzz, ts_yyy_xxxzzz, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyzz, ts_yyy_xxyzzz, ts_yyy_xxzzzz, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzzz, ts_yyy_xzzzzz, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzzz, ts_yyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyy_xxxxxx[i] = 3.0 * tr_y_yy_xxxxxx[i] * fe_0 + ts_yyy_xxxxxx[i] * fe_0 + tr_y_yyy_xxxxxx[i] * pa_y[i];

        tr_y_yyyy_xxxxxy[i] = 3.0 * tr_y_yy_xxxxxy[i] * fe_0 + tr_y_yyy_xxxxx[i] * fe_0 + ts_yyy_xxxxxy[i] * fe_0 + tr_y_yyy_xxxxxy[i] * pa_y[i];

        tr_y_yyyy_xxxxxz[i] = 3.0 * tr_y_yy_xxxxxz[i] * fe_0 + ts_yyy_xxxxxz[i] * fe_0 + tr_y_yyy_xxxxxz[i] * pa_y[i];

        tr_y_yyyy_xxxxyy[i] = 3.0 * tr_y_yy_xxxxyy[i] * fe_0 + 2.0 * tr_y_yyy_xxxxy[i] * fe_0 + ts_yyy_xxxxyy[i] * fe_0 + tr_y_yyy_xxxxyy[i] * pa_y[i];

        tr_y_yyyy_xxxxyz[i] = 3.0 * tr_y_yy_xxxxyz[i] * fe_0 + tr_y_yyy_xxxxz[i] * fe_0 + ts_yyy_xxxxyz[i] * fe_0 + tr_y_yyy_xxxxyz[i] * pa_y[i];

        tr_y_yyyy_xxxxzz[i] = 3.0 * tr_y_yy_xxxxzz[i] * fe_0 + ts_yyy_xxxxzz[i] * fe_0 + tr_y_yyy_xxxxzz[i] * pa_y[i];

        tr_y_yyyy_xxxyyy[i] = 3.0 * tr_y_yy_xxxyyy[i] * fe_0 + 3.0 * tr_y_yyy_xxxyy[i] * fe_0 + ts_yyy_xxxyyy[i] * fe_0 + tr_y_yyy_xxxyyy[i] * pa_y[i];

        tr_y_yyyy_xxxyyz[i] = 3.0 * tr_y_yy_xxxyyz[i] * fe_0 + 2.0 * tr_y_yyy_xxxyz[i] * fe_0 + ts_yyy_xxxyyz[i] * fe_0 + tr_y_yyy_xxxyyz[i] * pa_y[i];

        tr_y_yyyy_xxxyzz[i] = 3.0 * tr_y_yy_xxxyzz[i] * fe_0 + tr_y_yyy_xxxzz[i] * fe_0 + ts_yyy_xxxyzz[i] * fe_0 + tr_y_yyy_xxxyzz[i] * pa_y[i];

        tr_y_yyyy_xxxzzz[i] = 3.0 * tr_y_yy_xxxzzz[i] * fe_0 + ts_yyy_xxxzzz[i] * fe_0 + tr_y_yyy_xxxzzz[i] * pa_y[i];

        tr_y_yyyy_xxyyyy[i] = 3.0 * tr_y_yy_xxyyyy[i] * fe_0 + 4.0 * tr_y_yyy_xxyyy[i] * fe_0 + ts_yyy_xxyyyy[i] * fe_0 + tr_y_yyy_xxyyyy[i] * pa_y[i];

        tr_y_yyyy_xxyyyz[i] = 3.0 * tr_y_yy_xxyyyz[i] * fe_0 + 3.0 * tr_y_yyy_xxyyz[i] * fe_0 + ts_yyy_xxyyyz[i] * fe_0 + tr_y_yyy_xxyyyz[i] * pa_y[i];

        tr_y_yyyy_xxyyzz[i] = 3.0 * tr_y_yy_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyy_xxyzz[i] * fe_0 + ts_yyy_xxyyzz[i] * fe_0 + tr_y_yyy_xxyyzz[i] * pa_y[i];

        tr_y_yyyy_xxyzzz[i] = 3.0 * tr_y_yy_xxyzzz[i] * fe_0 + tr_y_yyy_xxzzz[i] * fe_0 + ts_yyy_xxyzzz[i] * fe_0 + tr_y_yyy_xxyzzz[i] * pa_y[i];

        tr_y_yyyy_xxzzzz[i] = 3.0 * tr_y_yy_xxzzzz[i] * fe_0 + ts_yyy_xxzzzz[i] * fe_0 + tr_y_yyy_xxzzzz[i] * pa_y[i];

        tr_y_yyyy_xyyyyy[i] = 3.0 * tr_y_yy_xyyyyy[i] * fe_0 + 5.0 * tr_y_yyy_xyyyy[i] * fe_0 + ts_yyy_xyyyyy[i] * fe_0 + tr_y_yyy_xyyyyy[i] * pa_y[i];

        tr_y_yyyy_xyyyyz[i] = 3.0 * tr_y_yy_xyyyyz[i] * fe_0 + 4.0 * tr_y_yyy_xyyyz[i] * fe_0 + ts_yyy_xyyyyz[i] * fe_0 + tr_y_yyy_xyyyyz[i] * pa_y[i];

        tr_y_yyyy_xyyyzz[i] = 3.0 * tr_y_yy_xyyyzz[i] * fe_0 + 3.0 * tr_y_yyy_xyyzz[i] * fe_0 + ts_yyy_xyyyzz[i] * fe_0 + tr_y_yyy_xyyyzz[i] * pa_y[i];

        tr_y_yyyy_xyyzzz[i] = 3.0 * tr_y_yy_xyyzzz[i] * fe_0 + 2.0 * tr_y_yyy_xyzzz[i] * fe_0 + ts_yyy_xyyzzz[i] * fe_0 + tr_y_yyy_xyyzzz[i] * pa_y[i];

        tr_y_yyyy_xyzzzz[i] = 3.0 * tr_y_yy_xyzzzz[i] * fe_0 + tr_y_yyy_xzzzz[i] * fe_0 + ts_yyy_xyzzzz[i] * fe_0 + tr_y_yyy_xyzzzz[i] * pa_y[i];

        tr_y_yyyy_xzzzzz[i] = 3.0 * tr_y_yy_xzzzzz[i] * fe_0 + ts_yyy_xzzzzz[i] * fe_0 + tr_y_yyy_xzzzzz[i] * pa_y[i];

        tr_y_yyyy_yyyyyy[i] = 3.0 * tr_y_yy_yyyyyy[i] * fe_0 + 6.0 * tr_y_yyy_yyyyy[i] * fe_0 + ts_yyy_yyyyyy[i] * fe_0 + tr_y_yyy_yyyyyy[i] * pa_y[i];

        tr_y_yyyy_yyyyyz[i] = 3.0 * tr_y_yy_yyyyyz[i] * fe_0 + 5.0 * tr_y_yyy_yyyyz[i] * fe_0 + ts_yyy_yyyyyz[i] * fe_0 + tr_y_yyy_yyyyyz[i] * pa_y[i];

        tr_y_yyyy_yyyyzz[i] = 3.0 * tr_y_yy_yyyyzz[i] * fe_0 + 4.0 * tr_y_yyy_yyyzz[i] * fe_0 + ts_yyy_yyyyzz[i] * fe_0 + tr_y_yyy_yyyyzz[i] * pa_y[i];

        tr_y_yyyy_yyyzzz[i] = 3.0 * tr_y_yy_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyy_yyzzz[i] * fe_0 + ts_yyy_yyyzzz[i] * fe_0 + tr_y_yyy_yyyzzz[i] * pa_y[i];

        tr_y_yyyy_yyzzzz[i] = 3.0 * tr_y_yy_yyzzzz[i] * fe_0 + 2.0 * tr_y_yyy_yzzzz[i] * fe_0 + ts_yyy_yyzzzz[i] * fe_0 + tr_y_yyy_yyzzzz[i] * pa_y[i];

        tr_y_yyyy_yzzzzz[i] = 3.0 * tr_y_yy_yzzzzz[i] * fe_0 + tr_y_yyy_zzzzz[i] * fe_0 + ts_yyy_yzzzzz[i] * fe_0 + tr_y_yyy_yzzzzz[i] * pa_y[i];

        tr_y_yyyy_zzzzzz[i] = 3.0 * tr_y_yy_zzzzzz[i] * fe_0 + ts_yyy_zzzzzz[i] * fe_0 + tr_y_yyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 728-756 components of targeted buffer : GI

    auto tr_y_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 728);

    auto tr_y_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 729);

    auto tr_y_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 730);

    auto tr_y_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 731);

    auto tr_y_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 732);

    auto tr_y_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 733);

    auto tr_y_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 734);

    auto tr_y_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 735);

    auto tr_y_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 736);

    auto tr_y_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 737);

    auto tr_y_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 738);

    auto tr_y_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 739);

    auto tr_y_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 740);

    auto tr_y_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 741);

    auto tr_y_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 742);

    auto tr_y_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 743);

    auto tr_y_yyyz_xyyyyz = pbuffer.data(idx_dip_gi + 744);

    auto tr_y_yyyz_xyyyzz = pbuffer.data(idx_dip_gi + 745);

    auto tr_y_yyyz_xyyzzz = pbuffer.data(idx_dip_gi + 746);

    auto tr_y_yyyz_xyzzzz = pbuffer.data(idx_dip_gi + 747);

    auto tr_y_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 748);

    auto tr_y_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 749);

    auto tr_y_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 750);

    auto tr_y_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 751);

    auto tr_y_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 752);

    auto tr_y_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 753);

    auto tr_y_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 754);

    auto tr_y_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 755);

    #pragma omp simd aligned(pa_z, tr_y_yyy_xxxxx, tr_y_yyy_xxxxxx, tr_y_yyy_xxxxxy, tr_y_yyy_xxxxxz, tr_y_yyy_xxxxy, tr_y_yyy_xxxxyy, tr_y_yyy_xxxxyz, tr_y_yyy_xxxxz, tr_y_yyy_xxxxzz, tr_y_yyy_xxxyy, tr_y_yyy_xxxyyy, tr_y_yyy_xxxyyz, tr_y_yyy_xxxyz, tr_y_yyy_xxxyzz, tr_y_yyy_xxxzz, tr_y_yyy_xxxzzz, tr_y_yyy_xxyyy, tr_y_yyy_xxyyyy, tr_y_yyy_xxyyyz, tr_y_yyy_xxyyz, tr_y_yyy_xxyyzz, tr_y_yyy_xxyzz, tr_y_yyy_xxyzzz, tr_y_yyy_xxzzz, tr_y_yyy_xxzzzz, tr_y_yyy_xyyyy, tr_y_yyy_xyyyyy, tr_y_yyy_xyyyyz, tr_y_yyy_xyyyz, tr_y_yyy_xyyyzz, tr_y_yyy_xyyzz, tr_y_yyy_xyyzzz, tr_y_yyy_xyzzz, tr_y_yyy_xyzzzz, tr_y_yyy_xzzzz, tr_y_yyy_xzzzzz, tr_y_yyy_yyyyy, tr_y_yyy_yyyyyy, tr_y_yyy_yyyyyz, tr_y_yyy_yyyyz, tr_y_yyy_yyyyzz, tr_y_yyy_yyyzz, tr_y_yyy_yyyzzz, tr_y_yyy_yyzzz, tr_y_yyy_yyzzzz, tr_y_yyy_yzzzz, tr_y_yyy_yzzzzz, tr_y_yyy_zzzzz, tr_y_yyy_zzzzzz, tr_y_yyyz_xxxxxx, tr_y_yyyz_xxxxxy, tr_y_yyyz_xxxxxz, tr_y_yyyz_xxxxyy, tr_y_yyyz_xxxxyz, tr_y_yyyz_xxxxzz, tr_y_yyyz_xxxyyy, tr_y_yyyz_xxxyyz, tr_y_yyyz_xxxyzz, tr_y_yyyz_xxxzzz, tr_y_yyyz_xxyyyy, tr_y_yyyz_xxyyyz, tr_y_yyyz_xxyyzz, tr_y_yyyz_xxyzzz, tr_y_yyyz_xxzzzz, tr_y_yyyz_xyyyyy, tr_y_yyyz_xyyyyz, tr_y_yyyz_xyyyzz, tr_y_yyyz_xyyzzz, tr_y_yyyz_xyzzzz, tr_y_yyyz_xzzzzz, tr_y_yyyz_yyyyyy, tr_y_yyyz_yyyyyz, tr_y_yyyz_yyyyzz, tr_y_yyyz_yyyzzz, tr_y_yyyz_yyzzzz, tr_y_yyyz_yzzzzz, tr_y_yyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyz_xxxxxx[i] = tr_y_yyy_xxxxxx[i] * pa_z[i];

        tr_y_yyyz_xxxxxy[i] = tr_y_yyy_xxxxxy[i] * pa_z[i];

        tr_y_yyyz_xxxxxz[i] = tr_y_yyy_xxxxx[i] * fe_0 + tr_y_yyy_xxxxxz[i] * pa_z[i];

        tr_y_yyyz_xxxxyy[i] = tr_y_yyy_xxxxyy[i] * pa_z[i];

        tr_y_yyyz_xxxxyz[i] = tr_y_yyy_xxxxy[i] * fe_0 + tr_y_yyy_xxxxyz[i] * pa_z[i];

        tr_y_yyyz_xxxxzz[i] = 2.0 * tr_y_yyy_xxxxz[i] * fe_0 + tr_y_yyy_xxxxzz[i] * pa_z[i];

        tr_y_yyyz_xxxyyy[i] = tr_y_yyy_xxxyyy[i] * pa_z[i];

        tr_y_yyyz_xxxyyz[i] = tr_y_yyy_xxxyy[i] * fe_0 + tr_y_yyy_xxxyyz[i] * pa_z[i];

        tr_y_yyyz_xxxyzz[i] = 2.0 * tr_y_yyy_xxxyz[i] * fe_0 + tr_y_yyy_xxxyzz[i] * pa_z[i];

        tr_y_yyyz_xxxzzz[i] = 3.0 * tr_y_yyy_xxxzz[i] * fe_0 + tr_y_yyy_xxxzzz[i] * pa_z[i];

        tr_y_yyyz_xxyyyy[i] = tr_y_yyy_xxyyyy[i] * pa_z[i];

        tr_y_yyyz_xxyyyz[i] = tr_y_yyy_xxyyy[i] * fe_0 + tr_y_yyy_xxyyyz[i] * pa_z[i];

        tr_y_yyyz_xxyyzz[i] = 2.0 * tr_y_yyy_xxyyz[i] * fe_0 + tr_y_yyy_xxyyzz[i] * pa_z[i];

        tr_y_yyyz_xxyzzz[i] = 3.0 * tr_y_yyy_xxyzz[i] * fe_0 + tr_y_yyy_xxyzzz[i] * pa_z[i];

        tr_y_yyyz_xxzzzz[i] = 4.0 * tr_y_yyy_xxzzz[i] * fe_0 + tr_y_yyy_xxzzzz[i] * pa_z[i];

        tr_y_yyyz_xyyyyy[i] = tr_y_yyy_xyyyyy[i] * pa_z[i];

        tr_y_yyyz_xyyyyz[i] = tr_y_yyy_xyyyy[i] * fe_0 + tr_y_yyy_xyyyyz[i] * pa_z[i];

        tr_y_yyyz_xyyyzz[i] = 2.0 * tr_y_yyy_xyyyz[i] * fe_0 + tr_y_yyy_xyyyzz[i] * pa_z[i];

        tr_y_yyyz_xyyzzz[i] = 3.0 * tr_y_yyy_xyyzz[i] * fe_0 + tr_y_yyy_xyyzzz[i] * pa_z[i];

        tr_y_yyyz_xyzzzz[i] = 4.0 * tr_y_yyy_xyzzz[i] * fe_0 + tr_y_yyy_xyzzzz[i] * pa_z[i];

        tr_y_yyyz_xzzzzz[i] = 5.0 * tr_y_yyy_xzzzz[i] * fe_0 + tr_y_yyy_xzzzzz[i] * pa_z[i];

        tr_y_yyyz_yyyyyy[i] = tr_y_yyy_yyyyyy[i] * pa_z[i];

        tr_y_yyyz_yyyyyz[i] = tr_y_yyy_yyyyy[i] * fe_0 + tr_y_yyy_yyyyyz[i] * pa_z[i];

        tr_y_yyyz_yyyyzz[i] = 2.0 * tr_y_yyy_yyyyz[i] * fe_0 + tr_y_yyy_yyyyzz[i] * pa_z[i];

        tr_y_yyyz_yyyzzz[i] = 3.0 * tr_y_yyy_yyyzz[i] * fe_0 + tr_y_yyy_yyyzzz[i] * pa_z[i];

        tr_y_yyyz_yyzzzz[i] = 4.0 * tr_y_yyy_yyzzz[i] * fe_0 + tr_y_yyy_yyzzzz[i] * pa_z[i];

        tr_y_yyyz_yzzzzz[i] = 5.0 * tr_y_yyy_yzzzz[i] * fe_0 + tr_y_yyy_yzzzzz[i] * pa_z[i];

        tr_y_yyyz_zzzzzz[i] = 6.0 * tr_y_yyy_zzzzz[i] * fe_0 + tr_y_yyy_zzzzzz[i] * pa_z[i];
    }

    // Set up 756-784 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yy_xxxxxx, tr_y_yy_xxxxxy, tr_y_yy_xxxxyy, tr_y_yy_xxxxyz, tr_y_yy_xxxyyy, tr_y_yy_xxxyyz, tr_y_yy_xxxyzz, tr_y_yy_xxyyyy, tr_y_yy_xxyyyz, tr_y_yy_xxyyzz, tr_y_yy_xxyzzz, tr_y_yy_xyyyyy, tr_y_yy_xyyyyz, tr_y_yy_xyyyzz, tr_y_yy_xyyzzz, tr_y_yy_xyzzzz, tr_y_yy_yyyyyy, tr_y_yy_yyyyyz, tr_y_yy_yyyyzz, tr_y_yy_yyyzzz, tr_y_yy_yyzzzz, tr_y_yy_yzzzzz, tr_y_yyz_xxxxxx, tr_y_yyz_xxxxxy, tr_y_yyz_xxxxy, tr_y_yyz_xxxxyy, tr_y_yyz_xxxxyz, tr_y_yyz_xxxyy, tr_y_yyz_xxxyyy, tr_y_yyz_xxxyyz, tr_y_yyz_xxxyz, tr_y_yyz_xxxyzz, tr_y_yyz_xxyyy, tr_y_yyz_xxyyyy, tr_y_yyz_xxyyyz, tr_y_yyz_xxyyz, tr_y_yyz_xxyyzz, tr_y_yyz_xxyzz, tr_y_yyz_xxyzzz, tr_y_yyz_xyyyy, tr_y_yyz_xyyyyy, tr_y_yyz_xyyyyz, tr_y_yyz_xyyyz, tr_y_yyz_xyyyzz, tr_y_yyz_xyyzz, tr_y_yyz_xyyzzz, tr_y_yyz_xyzzz, tr_y_yyz_xyzzzz, tr_y_yyz_yyyyy, tr_y_yyz_yyyyyy, tr_y_yyz_yyyyyz, tr_y_yyz_yyyyz, tr_y_yyz_yyyyzz, tr_y_yyz_yyyzz, tr_y_yyz_yyyzzz, tr_y_yyz_yyzzz, tr_y_yyz_yyzzzz, tr_y_yyz_yzzzz, tr_y_yyz_yzzzzz, tr_y_yyzz_xxxxxx, tr_y_yyzz_xxxxxy, tr_y_yyzz_xxxxxz, tr_y_yyzz_xxxxyy, tr_y_yyzz_xxxxyz, tr_y_yyzz_xxxxzz, tr_y_yyzz_xxxyyy, tr_y_yyzz_xxxyyz, tr_y_yyzz_xxxyzz, tr_y_yyzz_xxxzzz, tr_y_yyzz_xxyyyy, tr_y_yyzz_xxyyyz, tr_y_yyzz_xxyyzz, tr_y_yyzz_xxyzzz, tr_y_yyzz_xxzzzz, tr_y_yyzz_xyyyyy, tr_y_yyzz_xyyyyz, tr_y_yyzz_xyyyzz, tr_y_yyzz_xyyzzz, tr_y_yyzz_xyzzzz, tr_y_yyzz_xzzzzz, tr_y_yyzz_yyyyyy, tr_y_yyzz_yyyyyz, tr_y_yyzz_yyyyzz, tr_y_yyzz_yyyzzz, tr_y_yyzz_yyzzzz, tr_y_yyzz_yzzzzz, tr_y_yyzz_zzzzzz, tr_y_yzz_xxxxxz, tr_y_yzz_xxxxzz, tr_y_yzz_xxxzzz, tr_y_yzz_xxzzzz, tr_y_yzz_xzzzzz, tr_y_yzz_zzzzzz, tr_y_zz_xxxxxz, tr_y_zz_xxxxzz, tr_y_zz_xxxzzz, tr_y_zz_xxzzzz, tr_y_zz_xzzzzz, tr_y_zz_zzzzzz, ts_yzz_xxxxxz, ts_yzz_xxxxzz, ts_yzz_xxxzzz, ts_yzz_xxzzzz, ts_yzz_xzzzzz, ts_yzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzz_xxxxxx[i] = tr_y_yy_xxxxxx[i] * fe_0 + tr_y_yyz_xxxxxx[i] * pa_z[i];

        tr_y_yyzz_xxxxxy[i] = tr_y_yy_xxxxxy[i] * fe_0 + tr_y_yyz_xxxxxy[i] * pa_z[i];

        tr_y_yyzz_xxxxxz[i] = tr_y_zz_xxxxxz[i] * fe_0 + ts_yzz_xxxxxz[i] * fe_0 + tr_y_yzz_xxxxxz[i] * pa_y[i];

        tr_y_yyzz_xxxxyy[i] = tr_y_yy_xxxxyy[i] * fe_0 + tr_y_yyz_xxxxyy[i] * pa_z[i];

        tr_y_yyzz_xxxxyz[i] = tr_y_yy_xxxxyz[i] * fe_0 + tr_y_yyz_xxxxy[i] * fe_0 + tr_y_yyz_xxxxyz[i] * pa_z[i];

        tr_y_yyzz_xxxxzz[i] = tr_y_zz_xxxxzz[i] * fe_0 + ts_yzz_xxxxzz[i] * fe_0 + tr_y_yzz_xxxxzz[i] * pa_y[i];

        tr_y_yyzz_xxxyyy[i] = tr_y_yy_xxxyyy[i] * fe_0 + tr_y_yyz_xxxyyy[i] * pa_z[i];

        tr_y_yyzz_xxxyyz[i] = tr_y_yy_xxxyyz[i] * fe_0 + tr_y_yyz_xxxyy[i] * fe_0 + tr_y_yyz_xxxyyz[i] * pa_z[i];

        tr_y_yyzz_xxxyzz[i] = tr_y_yy_xxxyzz[i] * fe_0 + 2.0 * tr_y_yyz_xxxyz[i] * fe_0 + tr_y_yyz_xxxyzz[i] * pa_z[i];

        tr_y_yyzz_xxxzzz[i] = tr_y_zz_xxxzzz[i] * fe_0 + ts_yzz_xxxzzz[i] * fe_0 + tr_y_yzz_xxxzzz[i] * pa_y[i];

        tr_y_yyzz_xxyyyy[i] = tr_y_yy_xxyyyy[i] * fe_0 + tr_y_yyz_xxyyyy[i] * pa_z[i];

        tr_y_yyzz_xxyyyz[i] = tr_y_yy_xxyyyz[i] * fe_0 + tr_y_yyz_xxyyy[i] * fe_0 + tr_y_yyz_xxyyyz[i] * pa_z[i];

        tr_y_yyzz_xxyyzz[i] = tr_y_yy_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyz_xxyyz[i] * fe_0 + tr_y_yyz_xxyyzz[i] * pa_z[i];

        tr_y_yyzz_xxyzzz[i] = tr_y_yy_xxyzzz[i] * fe_0 + 3.0 * tr_y_yyz_xxyzz[i] * fe_0 + tr_y_yyz_xxyzzz[i] * pa_z[i];

        tr_y_yyzz_xxzzzz[i] = tr_y_zz_xxzzzz[i] * fe_0 + ts_yzz_xxzzzz[i] * fe_0 + tr_y_yzz_xxzzzz[i] * pa_y[i];

        tr_y_yyzz_xyyyyy[i] = tr_y_yy_xyyyyy[i] * fe_0 + tr_y_yyz_xyyyyy[i] * pa_z[i];

        tr_y_yyzz_xyyyyz[i] = tr_y_yy_xyyyyz[i] * fe_0 + tr_y_yyz_xyyyy[i] * fe_0 + tr_y_yyz_xyyyyz[i] * pa_z[i];

        tr_y_yyzz_xyyyzz[i] = tr_y_yy_xyyyzz[i] * fe_0 + 2.0 * tr_y_yyz_xyyyz[i] * fe_0 + tr_y_yyz_xyyyzz[i] * pa_z[i];

        tr_y_yyzz_xyyzzz[i] = tr_y_yy_xyyzzz[i] * fe_0 + 3.0 * tr_y_yyz_xyyzz[i] * fe_0 + tr_y_yyz_xyyzzz[i] * pa_z[i];

        tr_y_yyzz_xyzzzz[i] = tr_y_yy_xyzzzz[i] * fe_0 + 4.0 * tr_y_yyz_xyzzz[i] * fe_0 + tr_y_yyz_xyzzzz[i] * pa_z[i];

        tr_y_yyzz_xzzzzz[i] = tr_y_zz_xzzzzz[i] * fe_0 + ts_yzz_xzzzzz[i] * fe_0 + tr_y_yzz_xzzzzz[i] * pa_y[i];

        tr_y_yyzz_yyyyyy[i] = tr_y_yy_yyyyyy[i] * fe_0 + tr_y_yyz_yyyyyy[i] * pa_z[i];

        tr_y_yyzz_yyyyyz[i] = tr_y_yy_yyyyyz[i] * fe_0 + tr_y_yyz_yyyyy[i] * fe_0 + tr_y_yyz_yyyyyz[i] * pa_z[i];

        tr_y_yyzz_yyyyzz[i] = tr_y_yy_yyyyzz[i] * fe_0 + 2.0 * tr_y_yyz_yyyyz[i] * fe_0 + tr_y_yyz_yyyyzz[i] * pa_z[i];

        tr_y_yyzz_yyyzzz[i] = tr_y_yy_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyz_yyyzz[i] * fe_0 + tr_y_yyz_yyyzzz[i] * pa_z[i];

        tr_y_yyzz_yyzzzz[i] = tr_y_yy_yyzzzz[i] * fe_0 + 4.0 * tr_y_yyz_yyzzz[i] * fe_0 + tr_y_yyz_yyzzzz[i] * pa_z[i];

        tr_y_yyzz_yzzzzz[i] = tr_y_yy_yzzzzz[i] * fe_0 + 5.0 * tr_y_yyz_yzzzz[i] * fe_0 + tr_y_yyz_yzzzzz[i] * pa_z[i];

        tr_y_yyzz_zzzzzz[i] = tr_y_zz_zzzzzz[i] * fe_0 + ts_yzz_zzzzzz[i] * fe_0 + tr_y_yzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 784-812 components of targeted buffer : GI

    auto tr_y_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 784);

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

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yz_xxxxxy, tr_y_yz_xxxxyy, tr_y_yz_xxxyyy, tr_y_yz_xxyyyy, tr_y_yz_xyyyyy, tr_y_yz_yyyyyy, tr_y_yzz_xxxxxy, tr_y_yzz_xxxxyy, tr_y_yzz_xxxyyy, tr_y_yzz_xxyyyy, tr_y_yzz_xyyyyy, tr_y_yzz_yyyyyy, tr_y_yzzz_xxxxxx, tr_y_yzzz_xxxxxy, tr_y_yzzz_xxxxxz, tr_y_yzzz_xxxxyy, tr_y_yzzz_xxxxyz, tr_y_yzzz_xxxxzz, tr_y_yzzz_xxxyyy, tr_y_yzzz_xxxyyz, tr_y_yzzz_xxxyzz, tr_y_yzzz_xxxzzz, tr_y_yzzz_xxyyyy, tr_y_yzzz_xxyyyz, tr_y_yzzz_xxyyzz, tr_y_yzzz_xxyzzz, tr_y_yzzz_xxzzzz, tr_y_yzzz_xyyyyy, tr_y_yzzz_xyyyyz, tr_y_yzzz_xyyyzz, tr_y_yzzz_xyyzzz, tr_y_yzzz_xyzzzz, tr_y_yzzz_xzzzzz, tr_y_yzzz_yyyyyy, tr_y_yzzz_yyyyyz, tr_y_yzzz_yyyyzz, tr_y_yzzz_yyyzzz, tr_y_yzzz_yyzzzz, tr_y_yzzz_yzzzzz, tr_y_yzzz_zzzzzz, tr_y_zzz_xxxxxx, tr_y_zzz_xxxxxz, tr_y_zzz_xxxxyz, tr_y_zzz_xxxxz, tr_y_zzz_xxxxzz, tr_y_zzz_xxxyyz, tr_y_zzz_xxxyz, tr_y_zzz_xxxyzz, tr_y_zzz_xxxzz, tr_y_zzz_xxxzzz, tr_y_zzz_xxyyyz, tr_y_zzz_xxyyz, tr_y_zzz_xxyyzz, tr_y_zzz_xxyzz, tr_y_zzz_xxyzzz, tr_y_zzz_xxzzz, tr_y_zzz_xxzzzz, tr_y_zzz_xyyyyz, tr_y_zzz_xyyyz, tr_y_zzz_xyyyzz, tr_y_zzz_xyyzz, tr_y_zzz_xyyzzz, tr_y_zzz_xyzzz, tr_y_zzz_xyzzzz, tr_y_zzz_xzzzz, tr_y_zzz_xzzzzz, tr_y_zzz_yyyyyz, tr_y_zzz_yyyyz, tr_y_zzz_yyyyzz, tr_y_zzz_yyyzz, tr_y_zzz_yyyzzz, tr_y_zzz_yyzzz, tr_y_zzz_yyzzzz, tr_y_zzz_yzzzz, tr_y_zzz_yzzzzz, tr_y_zzz_zzzzz, tr_y_zzz_zzzzzz, ts_zzz_xxxxxx, ts_zzz_xxxxxz, ts_zzz_xxxxyz, ts_zzz_xxxxzz, ts_zzz_xxxyyz, ts_zzz_xxxyzz, ts_zzz_xxxzzz, ts_zzz_xxyyyz, ts_zzz_xxyyzz, ts_zzz_xxyzzz, ts_zzz_xxzzzz, ts_zzz_xyyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzz_xxxxxx[i] = ts_zzz_xxxxxx[i] * fe_0 + tr_y_zzz_xxxxxx[i] * pa_y[i];

        tr_y_yzzz_xxxxxy[i] = 2.0 * tr_y_yz_xxxxxy[i] * fe_0 + tr_y_yzz_xxxxxy[i] * pa_z[i];

        tr_y_yzzz_xxxxxz[i] = ts_zzz_xxxxxz[i] * fe_0 + tr_y_zzz_xxxxxz[i] * pa_y[i];

        tr_y_yzzz_xxxxyy[i] = 2.0 * tr_y_yz_xxxxyy[i] * fe_0 + tr_y_yzz_xxxxyy[i] * pa_z[i];

        tr_y_yzzz_xxxxyz[i] = tr_y_zzz_xxxxz[i] * fe_0 + ts_zzz_xxxxyz[i] * fe_0 + tr_y_zzz_xxxxyz[i] * pa_y[i];

        tr_y_yzzz_xxxxzz[i] = ts_zzz_xxxxzz[i] * fe_0 + tr_y_zzz_xxxxzz[i] * pa_y[i];

        tr_y_yzzz_xxxyyy[i] = 2.0 * tr_y_yz_xxxyyy[i] * fe_0 + tr_y_yzz_xxxyyy[i] * pa_z[i];

        tr_y_yzzz_xxxyyz[i] = 2.0 * tr_y_zzz_xxxyz[i] * fe_0 + ts_zzz_xxxyyz[i] * fe_0 + tr_y_zzz_xxxyyz[i] * pa_y[i];

        tr_y_yzzz_xxxyzz[i] = tr_y_zzz_xxxzz[i] * fe_0 + ts_zzz_xxxyzz[i] * fe_0 + tr_y_zzz_xxxyzz[i] * pa_y[i];

        tr_y_yzzz_xxxzzz[i] = ts_zzz_xxxzzz[i] * fe_0 + tr_y_zzz_xxxzzz[i] * pa_y[i];

        tr_y_yzzz_xxyyyy[i] = 2.0 * tr_y_yz_xxyyyy[i] * fe_0 + tr_y_yzz_xxyyyy[i] * pa_z[i];

        tr_y_yzzz_xxyyyz[i] = 3.0 * tr_y_zzz_xxyyz[i] * fe_0 + ts_zzz_xxyyyz[i] * fe_0 + tr_y_zzz_xxyyyz[i] * pa_y[i];

        tr_y_yzzz_xxyyzz[i] = 2.0 * tr_y_zzz_xxyzz[i] * fe_0 + ts_zzz_xxyyzz[i] * fe_0 + tr_y_zzz_xxyyzz[i] * pa_y[i];

        tr_y_yzzz_xxyzzz[i] = tr_y_zzz_xxzzz[i] * fe_0 + ts_zzz_xxyzzz[i] * fe_0 + tr_y_zzz_xxyzzz[i] * pa_y[i];

        tr_y_yzzz_xxzzzz[i] = ts_zzz_xxzzzz[i] * fe_0 + tr_y_zzz_xxzzzz[i] * pa_y[i];

        tr_y_yzzz_xyyyyy[i] = 2.0 * tr_y_yz_xyyyyy[i] * fe_0 + tr_y_yzz_xyyyyy[i] * pa_z[i];

        tr_y_yzzz_xyyyyz[i] = 4.0 * tr_y_zzz_xyyyz[i] * fe_0 + ts_zzz_xyyyyz[i] * fe_0 + tr_y_zzz_xyyyyz[i] * pa_y[i];

        tr_y_yzzz_xyyyzz[i] = 3.0 * tr_y_zzz_xyyzz[i] * fe_0 + ts_zzz_xyyyzz[i] * fe_0 + tr_y_zzz_xyyyzz[i] * pa_y[i];

        tr_y_yzzz_xyyzzz[i] = 2.0 * tr_y_zzz_xyzzz[i] * fe_0 + ts_zzz_xyyzzz[i] * fe_0 + tr_y_zzz_xyyzzz[i] * pa_y[i];

        tr_y_yzzz_xyzzzz[i] = tr_y_zzz_xzzzz[i] * fe_0 + ts_zzz_xyzzzz[i] * fe_0 + tr_y_zzz_xyzzzz[i] * pa_y[i];

        tr_y_yzzz_xzzzzz[i] = ts_zzz_xzzzzz[i] * fe_0 + tr_y_zzz_xzzzzz[i] * pa_y[i];

        tr_y_yzzz_yyyyyy[i] = 2.0 * tr_y_yz_yyyyyy[i] * fe_0 + tr_y_yzz_yyyyyy[i] * pa_z[i];

        tr_y_yzzz_yyyyyz[i] = 5.0 * tr_y_zzz_yyyyz[i] * fe_0 + ts_zzz_yyyyyz[i] * fe_0 + tr_y_zzz_yyyyyz[i] * pa_y[i];

        tr_y_yzzz_yyyyzz[i] = 4.0 * tr_y_zzz_yyyzz[i] * fe_0 + ts_zzz_yyyyzz[i] * fe_0 + tr_y_zzz_yyyyzz[i] * pa_y[i];

        tr_y_yzzz_yyyzzz[i] = 3.0 * tr_y_zzz_yyzzz[i] * fe_0 + ts_zzz_yyyzzz[i] * fe_0 + tr_y_zzz_yyyzzz[i] * pa_y[i];

        tr_y_yzzz_yyzzzz[i] = 2.0 * tr_y_zzz_yzzzz[i] * fe_0 + ts_zzz_yyzzzz[i] * fe_0 + tr_y_zzz_yyzzzz[i] * pa_y[i];

        tr_y_yzzz_yzzzzz[i] = tr_y_zzz_zzzzz[i] * fe_0 + ts_zzz_yzzzzz[i] * fe_0 + tr_y_zzz_yzzzzz[i] * pa_y[i];

        tr_y_yzzz_zzzzzz[i] = ts_zzz_zzzzzz[i] * fe_0 + tr_y_zzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 812-840 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_z, tr_y_zz_xxxxxx, tr_y_zz_xxxxxy, tr_y_zz_xxxxxz, tr_y_zz_xxxxyy, tr_y_zz_xxxxyz, tr_y_zz_xxxxzz, tr_y_zz_xxxyyy, tr_y_zz_xxxyyz, tr_y_zz_xxxyzz, tr_y_zz_xxxzzz, tr_y_zz_xxyyyy, tr_y_zz_xxyyyz, tr_y_zz_xxyyzz, tr_y_zz_xxyzzz, tr_y_zz_xxzzzz, tr_y_zz_xyyyyy, tr_y_zz_xyyyyz, tr_y_zz_xyyyzz, tr_y_zz_xyyzzz, tr_y_zz_xyzzzz, tr_y_zz_xzzzzz, tr_y_zz_yyyyyy, tr_y_zz_yyyyyz, tr_y_zz_yyyyzz, tr_y_zz_yyyzzz, tr_y_zz_yyzzzz, tr_y_zz_yzzzzz, tr_y_zz_zzzzzz, tr_y_zzz_xxxxx, tr_y_zzz_xxxxxx, tr_y_zzz_xxxxxy, tr_y_zzz_xxxxxz, tr_y_zzz_xxxxy, tr_y_zzz_xxxxyy, tr_y_zzz_xxxxyz, tr_y_zzz_xxxxz, tr_y_zzz_xxxxzz, tr_y_zzz_xxxyy, tr_y_zzz_xxxyyy, tr_y_zzz_xxxyyz, tr_y_zzz_xxxyz, tr_y_zzz_xxxyzz, tr_y_zzz_xxxzz, tr_y_zzz_xxxzzz, tr_y_zzz_xxyyy, tr_y_zzz_xxyyyy, tr_y_zzz_xxyyyz, tr_y_zzz_xxyyz, tr_y_zzz_xxyyzz, tr_y_zzz_xxyzz, tr_y_zzz_xxyzzz, tr_y_zzz_xxzzz, tr_y_zzz_xxzzzz, tr_y_zzz_xyyyy, tr_y_zzz_xyyyyy, tr_y_zzz_xyyyyz, tr_y_zzz_xyyyz, tr_y_zzz_xyyyzz, tr_y_zzz_xyyzz, tr_y_zzz_xyyzzz, tr_y_zzz_xyzzz, tr_y_zzz_xyzzzz, tr_y_zzz_xzzzz, tr_y_zzz_xzzzzz, tr_y_zzz_yyyyy, tr_y_zzz_yyyyyy, tr_y_zzz_yyyyyz, tr_y_zzz_yyyyz, tr_y_zzz_yyyyzz, tr_y_zzz_yyyzz, tr_y_zzz_yyyzzz, tr_y_zzz_yyzzz, tr_y_zzz_yyzzzz, tr_y_zzz_yzzzz, tr_y_zzz_yzzzzz, tr_y_zzz_zzzzz, tr_y_zzz_zzzzzz, tr_y_zzzz_xxxxxx, tr_y_zzzz_xxxxxy, tr_y_zzzz_xxxxxz, tr_y_zzzz_xxxxyy, tr_y_zzzz_xxxxyz, tr_y_zzzz_xxxxzz, tr_y_zzzz_xxxyyy, tr_y_zzzz_xxxyyz, tr_y_zzzz_xxxyzz, tr_y_zzzz_xxxzzz, tr_y_zzzz_xxyyyy, tr_y_zzzz_xxyyyz, tr_y_zzzz_xxyyzz, tr_y_zzzz_xxyzzz, tr_y_zzzz_xxzzzz, tr_y_zzzz_xyyyyy, tr_y_zzzz_xyyyyz, tr_y_zzzz_xyyyzz, tr_y_zzzz_xyyzzz, tr_y_zzzz_xyzzzz, tr_y_zzzz_xzzzzz, tr_y_zzzz_yyyyyy, tr_y_zzzz_yyyyyz, tr_y_zzzz_yyyyzz, tr_y_zzzz_yyyzzz, tr_y_zzzz_yyzzzz, tr_y_zzzz_yzzzzz, tr_y_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzz_xxxxxx[i] = 3.0 * tr_y_zz_xxxxxx[i] * fe_0 + tr_y_zzz_xxxxxx[i] * pa_z[i];

        tr_y_zzzz_xxxxxy[i] = 3.0 * tr_y_zz_xxxxxy[i] * fe_0 + tr_y_zzz_xxxxxy[i] * pa_z[i];

        tr_y_zzzz_xxxxxz[i] = 3.0 * tr_y_zz_xxxxxz[i] * fe_0 + tr_y_zzz_xxxxx[i] * fe_0 + tr_y_zzz_xxxxxz[i] * pa_z[i];

        tr_y_zzzz_xxxxyy[i] = 3.0 * tr_y_zz_xxxxyy[i] * fe_0 + tr_y_zzz_xxxxyy[i] * pa_z[i];

        tr_y_zzzz_xxxxyz[i] = 3.0 * tr_y_zz_xxxxyz[i] * fe_0 + tr_y_zzz_xxxxy[i] * fe_0 + tr_y_zzz_xxxxyz[i] * pa_z[i];

        tr_y_zzzz_xxxxzz[i] = 3.0 * tr_y_zz_xxxxzz[i] * fe_0 + 2.0 * tr_y_zzz_xxxxz[i] * fe_0 + tr_y_zzz_xxxxzz[i] * pa_z[i];

        tr_y_zzzz_xxxyyy[i] = 3.0 * tr_y_zz_xxxyyy[i] * fe_0 + tr_y_zzz_xxxyyy[i] * pa_z[i];

        tr_y_zzzz_xxxyyz[i] = 3.0 * tr_y_zz_xxxyyz[i] * fe_0 + tr_y_zzz_xxxyy[i] * fe_0 + tr_y_zzz_xxxyyz[i] * pa_z[i];

        tr_y_zzzz_xxxyzz[i] = 3.0 * tr_y_zz_xxxyzz[i] * fe_0 + 2.0 * tr_y_zzz_xxxyz[i] * fe_0 + tr_y_zzz_xxxyzz[i] * pa_z[i];

        tr_y_zzzz_xxxzzz[i] = 3.0 * tr_y_zz_xxxzzz[i] * fe_0 + 3.0 * tr_y_zzz_xxxzz[i] * fe_0 + tr_y_zzz_xxxzzz[i] * pa_z[i];

        tr_y_zzzz_xxyyyy[i] = 3.0 * tr_y_zz_xxyyyy[i] * fe_0 + tr_y_zzz_xxyyyy[i] * pa_z[i];

        tr_y_zzzz_xxyyyz[i] = 3.0 * tr_y_zz_xxyyyz[i] * fe_0 + tr_y_zzz_xxyyy[i] * fe_0 + tr_y_zzz_xxyyyz[i] * pa_z[i];

        tr_y_zzzz_xxyyzz[i] = 3.0 * tr_y_zz_xxyyzz[i] * fe_0 + 2.0 * tr_y_zzz_xxyyz[i] * fe_0 + tr_y_zzz_xxyyzz[i] * pa_z[i];

        tr_y_zzzz_xxyzzz[i] = 3.0 * tr_y_zz_xxyzzz[i] * fe_0 + 3.0 * tr_y_zzz_xxyzz[i] * fe_0 + tr_y_zzz_xxyzzz[i] * pa_z[i];

        tr_y_zzzz_xxzzzz[i] = 3.0 * tr_y_zz_xxzzzz[i] * fe_0 + 4.0 * tr_y_zzz_xxzzz[i] * fe_0 + tr_y_zzz_xxzzzz[i] * pa_z[i];

        tr_y_zzzz_xyyyyy[i] = 3.0 * tr_y_zz_xyyyyy[i] * fe_0 + tr_y_zzz_xyyyyy[i] * pa_z[i];

        tr_y_zzzz_xyyyyz[i] = 3.0 * tr_y_zz_xyyyyz[i] * fe_0 + tr_y_zzz_xyyyy[i] * fe_0 + tr_y_zzz_xyyyyz[i] * pa_z[i];

        tr_y_zzzz_xyyyzz[i] = 3.0 * tr_y_zz_xyyyzz[i] * fe_0 + 2.0 * tr_y_zzz_xyyyz[i] * fe_0 + tr_y_zzz_xyyyzz[i] * pa_z[i];

        tr_y_zzzz_xyyzzz[i] = 3.0 * tr_y_zz_xyyzzz[i] * fe_0 + 3.0 * tr_y_zzz_xyyzz[i] * fe_0 + tr_y_zzz_xyyzzz[i] * pa_z[i];

        tr_y_zzzz_xyzzzz[i] = 3.0 * tr_y_zz_xyzzzz[i] * fe_0 + 4.0 * tr_y_zzz_xyzzz[i] * fe_0 + tr_y_zzz_xyzzzz[i] * pa_z[i];

        tr_y_zzzz_xzzzzz[i] = 3.0 * tr_y_zz_xzzzzz[i] * fe_0 + 5.0 * tr_y_zzz_xzzzz[i] * fe_0 + tr_y_zzz_xzzzzz[i] * pa_z[i];

        tr_y_zzzz_yyyyyy[i] = 3.0 * tr_y_zz_yyyyyy[i] * fe_0 + tr_y_zzz_yyyyyy[i] * pa_z[i];

        tr_y_zzzz_yyyyyz[i] = 3.0 * tr_y_zz_yyyyyz[i] * fe_0 + tr_y_zzz_yyyyy[i] * fe_0 + tr_y_zzz_yyyyyz[i] * pa_z[i];

        tr_y_zzzz_yyyyzz[i] = 3.0 * tr_y_zz_yyyyzz[i] * fe_0 + 2.0 * tr_y_zzz_yyyyz[i] * fe_0 + tr_y_zzz_yyyyzz[i] * pa_z[i];

        tr_y_zzzz_yyyzzz[i] = 3.0 * tr_y_zz_yyyzzz[i] * fe_0 + 3.0 * tr_y_zzz_yyyzz[i] * fe_0 + tr_y_zzz_yyyzzz[i] * pa_z[i];

        tr_y_zzzz_yyzzzz[i] = 3.0 * tr_y_zz_yyzzzz[i] * fe_0 + 4.0 * tr_y_zzz_yyzzz[i] * fe_0 + tr_y_zzz_yyzzzz[i] * pa_z[i];

        tr_y_zzzz_yzzzzz[i] = 3.0 * tr_y_zz_yzzzzz[i] * fe_0 + 5.0 * tr_y_zzz_yzzzz[i] * fe_0 + tr_y_zzz_yzzzzz[i] * pa_z[i];

        tr_y_zzzz_zzzzzz[i] = 3.0 * tr_y_zz_zzzzzz[i] * fe_0 + 6.0 * tr_y_zzz_zzzzz[i] * fe_0 + tr_y_zzz_zzzzzz[i] * pa_z[i];
    }

    // Set up 840-868 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_z_xx_xxxxxx, tr_z_xx_xxxxxy, tr_z_xx_xxxxxz, tr_z_xx_xxxxyy, tr_z_xx_xxxxyz, tr_z_xx_xxxxzz, tr_z_xx_xxxyyy, tr_z_xx_xxxyyz, tr_z_xx_xxxyzz, tr_z_xx_xxxzzz, tr_z_xx_xxyyyy, tr_z_xx_xxyyyz, tr_z_xx_xxyyzz, tr_z_xx_xxyzzz, tr_z_xx_xxzzzz, tr_z_xx_xyyyyy, tr_z_xx_xyyyyz, tr_z_xx_xyyyzz, tr_z_xx_xyyzzz, tr_z_xx_xyzzzz, tr_z_xx_xzzzzz, tr_z_xx_yyyyyy, tr_z_xx_yyyyyz, tr_z_xx_yyyyzz, tr_z_xx_yyyzzz, tr_z_xx_yyzzzz, tr_z_xx_yzzzzz, tr_z_xx_zzzzzz, tr_z_xxx_xxxxx, tr_z_xxx_xxxxxx, tr_z_xxx_xxxxxy, tr_z_xxx_xxxxxz, tr_z_xxx_xxxxy, tr_z_xxx_xxxxyy, tr_z_xxx_xxxxyz, tr_z_xxx_xxxxz, tr_z_xxx_xxxxzz, tr_z_xxx_xxxyy, tr_z_xxx_xxxyyy, tr_z_xxx_xxxyyz, tr_z_xxx_xxxyz, tr_z_xxx_xxxyzz, tr_z_xxx_xxxzz, tr_z_xxx_xxxzzz, tr_z_xxx_xxyyy, tr_z_xxx_xxyyyy, tr_z_xxx_xxyyyz, tr_z_xxx_xxyyz, tr_z_xxx_xxyyzz, tr_z_xxx_xxyzz, tr_z_xxx_xxyzzz, tr_z_xxx_xxzzz, tr_z_xxx_xxzzzz, tr_z_xxx_xyyyy, tr_z_xxx_xyyyyy, tr_z_xxx_xyyyyz, tr_z_xxx_xyyyz, tr_z_xxx_xyyyzz, tr_z_xxx_xyyzz, tr_z_xxx_xyyzzz, tr_z_xxx_xyzzz, tr_z_xxx_xyzzzz, tr_z_xxx_xzzzz, tr_z_xxx_xzzzzz, tr_z_xxx_yyyyy, tr_z_xxx_yyyyyy, tr_z_xxx_yyyyyz, tr_z_xxx_yyyyz, tr_z_xxx_yyyyzz, tr_z_xxx_yyyzz, tr_z_xxx_yyyzzz, tr_z_xxx_yyzzz, tr_z_xxx_yyzzzz, tr_z_xxx_yzzzz, tr_z_xxx_yzzzzz, tr_z_xxx_zzzzz, tr_z_xxx_zzzzzz, tr_z_xxxx_xxxxxx, tr_z_xxxx_xxxxxy, tr_z_xxxx_xxxxxz, tr_z_xxxx_xxxxyy, tr_z_xxxx_xxxxyz, tr_z_xxxx_xxxxzz, tr_z_xxxx_xxxyyy, tr_z_xxxx_xxxyyz, tr_z_xxxx_xxxyzz, tr_z_xxxx_xxxzzz, tr_z_xxxx_xxyyyy, tr_z_xxxx_xxyyyz, tr_z_xxxx_xxyyzz, tr_z_xxxx_xxyzzz, tr_z_xxxx_xxzzzz, tr_z_xxxx_xyyyyy, tr_z_xxxx_xyyyyz, tr_z_xxxx_xyyyzz, tr_z_xxxx_xyyzzz, tr_z_xxxx_xyzzzz, tr_z_xxxx_xzzzzz, tr_z_xxxx_yyyyyy, tr_z_xxxx_yyyyyz, tr_z_xxxx_yyyyzz, tr_z_xxxx_yyyzzz, tr_z_xxxx_yyzzzz, tr_z_xxxx_yzzzzz, tr_z_xxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxx_xxxxxx[i] = 3.0 * tr_z_xx_xxxxxx[i] * fe_0 + 6.0 * tr_z_xxx_xxxxx[i] * fe_0 + tr_z_xxx_xxxxxx[i] * pa_x[i];

        tr_z_xxxx_xxxxxy[i] = 3.0 * tr_z_xx_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxx_xxxxy[i] * fe_0 + tr_z_xxx_xxxxxy[i] * pa_x[i];

        tr_z_xxxx_xxxxxz[i] = 3.0 * tr_z_xx_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxx_xxxxz[i] * fe_0 + tr_z_xxx_xxxxxz[i] * pa_x[i];

        tr_z_xxxx_xxxxyy[i] = 3.0 * tr_z_xx_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxx_xxxyy[i] * fe_0 + tr_z_xxx_xxxxyy[i] * pa_x[i];

        tr_z_xxxx_xxxxyz[i] = 3.0 * tr_z_xx_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxx_xxxyz[i] * fe_0 + tr_z_xxx_xxxxyz[i] * pa_x[i];

        tr_z_xxxx_xxxxzz[i] = 3.0 * tr_z_xx_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxx_xxxzz[i] * fe_0 + tr_z_xxx_xxxxzz[i] * pa_x[i];

        tr_z_xxxx_xxxyyy[i] = 3.0 * tr_z_xx_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxx_xxyyy[i] * fe_0 + tr_z_xxx_xxxyyy[i] * pa_x[i];

        tr_z_xxxx_xxxyyz[i] = 3.0 * tr_z_xx_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxx_xxyyz[i] * fe_0 + tr_z_xxx_xxxyyz[i] * pa_x[i];

        tr_z_xxxx_xxxyzz[i] = 3.0 * tr_z_xx_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxx_xxyzz[i] * fe_0 + tr_z_xxx_xxxyzz[i] * pa_x[i];

        tr_z_xxxx_xxxzzz[i] = 3.0 * tr_z_xx_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxx_xxzzz[i] * fe_0 + tr_z_xxx_xxxzzz[i] * pa_x[i];

        tr_z_xxxx_xxyyyy[i] = 3.0 * tr_z_xx_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxx_xyyyy[i] * fe_0 + tr_z_xxx_xxyyyy[i] * pa_x[i];

        tr_z_xxxx_xxyyyz[i] = 3.0 * tr_z_xx_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxx_xyyyz[i] * fe_0 + tr_z_xxx_xxyyyz[i] * pa_x[i];

        tr_z_xxxx_xxyyzz[i] = 3.0 * tr_z_xx_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxx_xyyzz[i] * fe_0 + tr_z_xxx_xxyyzz[i] * pa_x[i];

        tr_z_xxxx_xxyzzz[i] = 3.0 * tr_z_xx_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxx_xyzzz[i] * fe_0 + tr_z_xxx_xxyzzz[i] * pa_x[i];

        tr_z_xxxx_xxzzzz[i] = 3.0 * tr_z_xx_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxx_xzzzz[i] * fe_0 + tr_z_xxx_xxzzzz[i] * pa_x[i];

        tr_z_xxxx_xyyyyy[i] = 3.0 * tr_z_xx_xyyyyy[i] * fe_0 + tr_z_xxx_yyyyy[i] * fe_0 + tr_z_xxx_xyyyyy[i] * pa_x[i];

        tr_z_xxxx_xyyyyz[i] = 3.0 * tr_z_xx_xyyyyz[i] * fe_0 + tr_z_xxx_yyyyz[i] * fe_0 + tr_z_xxx_xyyyyz[i] * pa_x[i];

        tr_z_xxxx_xyyyzz[i] = 3.0 * tr_z_xx_xyyyzz[i] * fe_0 + tr_z_xxx_yyyzz[i] * fe_0 + tr_z_xxx_xyyyzz[i] * pa_x[i];

        tr_z_xxxx_xyyzzz[i] = 3.0 * tr_z_xx_xyyzzz[i] * fe_0 + tr_z_xxx_yyzzz[i] * fe_0 + tr_z_xxx_xyyzzz[i] * pa_x[i];

        tr_z_xxxx_xyzzzz[i] = 3.0 * tr_z_xx_xyzzzz[i] * fe_0 + tr_z_xxx_yzzzz[i] * fe_0 + tr_z_xxx_xyzzzz[i] * pa_x[i];

        tr_z_xxxx_xzzzzz[i] = 3.0 * tr_z_xx_xzzzzz[i] * fe_0 + tr_z_xxx_zzzzz[i] * fe_0 + tr_z_xxx_xzzzzz[i] * pa_x[i];

        tr_z_xxxx_yyyyyy[i] = 3.0 * tr_z_xx_yyyyyy[i] * fe_0 + tr_z_xxx_yyyyyy[i] * pa_x[i];

        tr_z_xxxx_yyyyyz[i] = 3.0 * tr_z_xx_yyyyyz[i] * fe_0 + tr_z_xxx_yyyyyz[i] * pa_x[i];

        tr_z_xxxx_yyyyzz[i] = 3.0 * tr_z_xx_yyyyzz[i] * fe_0 + tr_z_xxx_yyyyzz[i] * pa_x[i];

        tr_z_xxxx_yyyzzz[i] = 3.0 * tr_z_xx_yyyzzz[i] * fe_0 + tr_z_xxx_yyyzzz[i] * pa_x[i];

        tr_z_xxxx_yyzzzz[i] = 3.0 * tr_z_xx_yyzzzz[i] * fe_0 + tr_z_xxx_yyzzzz[i] * pa_x[i];

        tr_z_xxxx_yzzzzz[i] = 3.0 * tr_z_xx_yzzzzz[i] * fe_0 + tr_z_xxx_yzzzzz[i] * pa_x[i];

        tr_z_xxxx_zzzzzz[i] = 3.0 * tr_z_xx_zzzzzz[i] * fe_0 + tr_z_xxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 868-896 components of targeted buffer : GI

    auto tr_z_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 868);

    auto tr_z_xxxy_xxxxxy = pbuffer.data(idx_dip_gi + 869);

    auto tr_z_xxxy_xxxxxz = pbuffer.data(idx_dip_gi + 870);

    auto tr_z_xxxy_xxxxyy = pbuffer.data(idx_dip_gi + 871);

    auto tr_z_xxxy_xxxxyz = pbuffer.data(idx_dip_gi + 872);

    auto tr_z_xxxy_xxxxzz = pbuffer.data(idx_dip_gi + 873);

    auto tr_z_xxxy_xxxyyy = pbuffer.data(idx_dip_gi + 874);

    auto tr_z_xxxy_xxxyyz = pbuffer.data(idx_dip_gi + 875);

    auto tr_z_xxxy_xxxyzz = pbuffer.data(idx_dip_gi + 876);

    auto tr_z_xxxy_xxxzzz = pbuffer.data(idx_dip_gi + 877);

    auto tr_z_xxxy_xxyyyy = pbuffer.data(idx_dip_gi + 878);

    auto tr_z_xxxy_xxyyyz = pbuffer.data(idx_dip_gi + 879);

    auto tr_z_xxxy_xxyyzz = pbuffer.data(idx_dip_gi + 880);

    auto tr_z_xxxy_xxyzzz = pbuffer.data(idx_dip_gi + 881);

    auto tr_z_xxxy_xxzzzz = pbuffer.data(idx_dip_gi + 882);

    auto tr_z_xxxy_xyyyyy = pbuffer.data(idx_dip_gi + 883);

    auto tr_z_xxxy_xyyyyz = pbuffer.data(idx_dip_gi + 884);

    auto tr_z_xxxy_xyyyzz = pbuffer.data(idx_dip_gi + 885);

    auto tr_z_xxxy_xyyzzz = pbuffer.data(idx_dip_gi + 886);

    auto tr_z_xxxy_xyzzzz = pbuffer.data(idx_dip_gi + 887);

    auto tr_z_xxxy_xzzzzz = pbuffer.data(idx_dip_gi + 888);

    auto tr_z_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 889);

    auto tr_z_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 890);

    auto tr_z_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 891);

    auto tr_z_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 892);

    auto tr_z_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 893);

    auto tr_z_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 894);

    auto tr_z_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 895);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxx_xxxxx, tr_z_xxx_xxxxxx, tr_z_xxx_xxxxxy, tr_z_xxx_xxxxxz, tr_z_xxx_xxxxy, tr_z_xxx_xxxxyy, tr_z_xxx_xxxxyz, tr_z_xxx_xxxxz, tr_z_xxx_xxxxzz, tr_z_xxx_xxxyy, tr_z_xxx_xxxyyy, tr_z_xxx_xxxyyz, tr_z_xxx_xxxyz, tr_z_xxx_xxxyzz, tr_z_xxx_xxxzz, tr_z_xxx_xxxzzz, tr_z_xxx_xxyyy, tr_z_xxx_xxyyyy, tr_z_xxx_xxyyyz, tr_z_xxx_xxyyz, tr_z_xxx_xxyyzz, tr_z_xxx_xxyzz, tr_z_xxx_xxyzzz, tr_z_xxx_xxzzz, tr_z_xxx_xxzzzz, tr_z_xxx_xyyyy, tr_z_xxx_xyyyyy, tr_z_xxx_xyyyyz, tr_z_xxx_xyyyz, tr_z_xxx_xyyyzz, tr_z_xxx_xyyzz, tr_z_xxx_xyyzzz, tr_z_xxx_xyzzz, tr_z_xxx_xyzzzz, tr_z_xxx_xzzzz, tr_z_xxx_xzzzzz, tr_z_xxx_zzzzzz, tr_z_xxxy_xxxxxx, tr_z_xxxy_xxxxxy, tr_z_xxxy_xxxxxz, tr_z_xxxy_xxxxyy, tr_z_xxxy_xxxxyz, tr_z_xxxy_xxxxzz, tr_z_xxxy_xxxyyy, tr_z_xxxy_xxxyyz, tr_z_xxxy_xxxyzz, tr_z_xxxy_xxxzzz, tr_z_xxxy_xxyyyy, tr_z_xxxy_xxyyyz, tr_z_xxxy_xxyyzz, tr_z_xxxy_xxyzzz, tr_z_xxxy_xxzzzz, tr_z_xxxy_xyyyyy, tr_z_xxxy_xyyyyz, tr_z_xxxy_xyyyzz, tr_z_xxxy_xyyzzz, tr_z_xxxy_xyzzzz, tr_z_xxxy_xzzzzz, tr_z_xxxy_yyyyyy, tr_z_xxxy_yyyyyz, tr_z_xxxy_yyyyzz, tr_z_xxxy_yyyzzz, tr_z_xxxy_yyzzzz, tr_z_xxxy_yzzzzz, tr_z_xxxy_zzzzzz, tr_z_xxy_yyyyyy, tr_z_xxy_yyyyyz, tr_z_xxy_yyyyzz, tr_z_xxy_yyyzzz, tr_z_xxy_yyzzzz, tr_z_xxy_yzzzzz, tr_z_xy_yyyyyy, tr_z_xy_yyyyyz, tr_z_xy_yyyyzz, tr_z_xy_yyyzzz, tr_z_xy_yyzzzz, tr_z_xy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxy_xxxxxx[i] = tr_z_xxx_xxxxxx[i] * pa_y[i];

        tr_z_xxxy_xxxxxy[i] = tr_z_xxx_xxxxx[i] * fe_0 + tr_z_xxx_xxxxxy[i] * pa_y[i];

        tr_z_xxxy_xxxxxz[i] = tr_z_xxx_xxxxxz[i] * pa_y[i];

        tr_z_xxxy_xxxxyy[i] = 2.0 * tr_z_xxx_xxxxy[i] * fe_0 + tr_z_xxx_xxxxyy[i] * pa_y[i];

        tr_z_xxxy_xxxxyz[i] = tr_z_xxx_xxxxz[i] * fe_0 + tr_z_xxx_xxxxyz[i] * pa_y[i];

        tr_z_xxxy_xxxxzz[i] = tr_z_xxx_xxxxzz[i] * pa_y[i];

        tr_z_xxxy_xxxyyy[i] = 3.0 * tr_z_xxx_xxxyy[i] * fe_0 + tr_z_xxx_xxxyyy[i] * pa_y[i];

        tr_z_xxxy_xxxyyz[i] = 2.0 * tr_z_xxx_xxxyz[i] * fe_0 + tr_z_xxx_xxxyyz[i] * pa_y[i];

        tr_z_xxxy_xxxyzz[i] = tr_z_xxx_xxxzz[i] * fe_0 + tr_z_xxx_xxxyzz[i] * pa_y[i];

        tr_z_xxxy_xxxzzz[i] = tr_z_xxx_xxxzzz[i] * pa_y[i];

        tr_z_xxxy_xxyyyy[i] = 4.0 * tr_z_xxx_xxyyy[i] * fe_0 + tr_z_xxx_xxyyyy[i] * pa_y[i];

        tr_z_xxxy_xxyyyz[i] = 3.0 * tr_z_xxx_xxyyz[i] * fe_0 + tr_z_xxx_xxyyyz[i] * pa_y[i];

        tr_z_xxxy_xxyyzz[i] = 2.0 * tr_z_xxx_xxyzz[i] * fe_0 + tr_z_xxx_xxyyzz[i] * pa_y[i];

        tr_z_xxxy_xxyzzz[i] = tr_z_xxx_xxzzz[i] * fe_0 + tr_z_xxx_xxyzzz[i] * pa_y[i];

        tr_z_xxxy_xxzzzz[i] = tr_z_xxx_xxzzzz[i] * pa_y[i];

        tr_z_xxxy_xyyyyy[i] = 5.0 * tr_z_xxx_xyyyy[i] * fe_0 + tr_z_xxx_xyyyyy[i] * pa_y[i];

        tr_z_xxxy_xyyyyz[i] = 4.0 * tr_z_xxx_xyyyz[i] * fe_0 + tr_z_xxx_xyyyyz[i] * pa_y[i];

        tr_z_xxxy_xyyyzz[i] = 3.0 * tr_z_xxx_xyyzz[i] * fe_0 + tr_z_xxx_xyyyzz[i] * pa_y[i];

        tr_z_xxxy_xyyzzz[i] = 2.0 * tr_z_xxx_xyzzz[i] * fe_0 + tr_z_xxx_xyyzzz[i] * pa_y[i];

        tr_z_xxxy_xyzzzz[i] = tr_z_xxx_xzzzz[i] * fe_0 + tr_z_xxx_xyzzzz[i] * pa_y[i];

        tr_z_xxxy_xzzzzz[i] = tr_z_xxx_xzzzzz[i] * pa_y[i];

        tr_z_xxxy_yyyyyy[i] = 2.0 * tr_z_xy_yyyyyy[i] * fe_0 + tr_z_xxy_yyyyyy[i] * pa_x[i];

        tr_z_xxxy_yyyyyz[i] = 2.0 * tr_z_xy_yyyyyz[i] * fe_0 + tr_z_xxy_yyyyyz[i] * pa_x[i];

        tr_z_xxxy_yyyyzz[i] = 2.0 * tr_z_xy_yyyyzz[i] * fe_0 + tr_z_xxy_yyyyzz[i] * pa_x[i];

        tr_z_xxxy_yyyzzz[i] = 2.0 * tr_z_xy_yyyzzz[i] * fe_0 + tr_z_xxy_yyyzzz[i] * pa_x[i];

        tr_z_xxxy_yyzzzz[i] = 2.0 * tr_z_xy_yyzzzz[i] * fe_0 + tr_z_xxy_yyzzzz[i] * pa_x[i];

        tr_z_xxxy_yzzzzz[i] = 2.0 * tr_z_xy_yzzzzz[i] * fe_0 + tr_z_xxy_yzzzzz[i] * pa_x[i];

        tr_z_xxxy_zzzzzz[i] = tr_z_xxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 896-924 components of targeted buffer : GI

    auto tr_z_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 896);

    auto tr_z_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 897);

    auto tr_z_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 898);

    auto tr_z_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 899);

    auto tr_z_xxxz_xxxxyz = pbuffer.data(idx_dip_gi + 900);

    auto tr_z_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 901);

    auto tr_z_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 902);

    auto tr_z_xxxz_xxxyyz = pbuffer.data(idx_dip_gi + 903);

    auto tr_z_xxxz_xxxyzz = pbuffer.data(idx_dip_gi + 904);

    auto tr_z_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 905);

    auto tr_z_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 906);

    auto tr_z_xxxz_xxyyyz = pbuffer.data(idx_dip_gi + 907);

    auto tr_z_xxxz_xxyyzz = pbuffer.data(idx_dip_gi + 908);

    auto tr_z_xxxz_xxyzzz = pbuffer.data(idx_dip_gi + 909);

    auto tr_z_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 910);

    auto tr_z_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 911);

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

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xxx_xxxxxx, tr_z_xxx_xxxxxy, tr_z_xxx_xxxxyy, tr_z_xxx_xxxyyy, tr_z_xxx_xxyyyy, tr_z_xxx_xyyyyy, tr_z_xxxz_xxxxxx, tr_z_xxxz_xxxxxy, tr_z_xxxz_xxxxxz, tr_z_xxxz_xxxxyy, tr_z_xxxz_xxxxyz, tr_z_xxxz_xxxxzz, tr_z_xxxz_xxxyyy, tr_z_xxxz_xxxyyz, tr_z_xxxz_xxxyzz, tr_z_xxxz_xxxzzz, tr_z_xxxz_xxyyyy, tr_z_xxxz_xxyyyz, tr_z_xxxz_xxyyzz, tr_z_xxxz_xxyzzz, tr_z_xxxz_xxzzzz, tr_z_xxxz_xyyyyy, tr_z_xxxz_xyyyyz, tr_z_xxxz_xyyyzz, tr_z_xxxz_xyyzzz, tr_z_xxxz_xyzzzz, tr_z_xxxz_xzzzzz, tr_z_xxxz_yyyyyy, tr_z_xxxz_yyyyyz, tr_z_xxxz_yyyyzz, tr_z_xxxz_yyyzzz, tr_z_xxxz_yyzzzz, tr_z_xxxz_yzzzzz, tr_z_xxxz_zzzzzz, tr_z_xxz_xxxxxz, tr_z_xxz_xxxxyz, tr_z_xxz_xxxxz, tr_z_xxz_xxxxzz, tr_z_xxz_xxxyyz, tr_z_xxz_xxxyz, tr_z_xxz_xxxyzz, tr_z_xxz_xxxzz, tr_z_xxz_xxxzzz, tr_z_xxz_xxyyyz, tr_z_xxz_xxyyz, tr_z_xxz_xxyyzz, tr_z_xxz_xxyzz, tr_z_xxz_xxyzzz, tr_z_xxz_xxzzz, tr_z_xxz_xxzzzz, tr_z_xxz_xyyyyz, tr_z_xxz_xyyyz, tr_z_xxz_xyyyzz, tr_z_xxz_xyyzz, tr_z_xxz_xyyzzz, tr_z_xxz_xyzzz, tr_z_xxz_xyzzzz, tr_z_xxz_xzzzz, tr_z_xxz_xzzzzz, tr_z_xxz_yyyyyy, tr_z_xxz_yyyyyz, tr_z_xxz_yyyyz, tr_z_xxz_yyyyzz, tr_z_xxz_yyyzz, tr_z_xxz_yyyzzz, tr_z_xxz_yyzzz, tr_z_xxz_yyzzzz, tr_z_xxz_yzzzz, tr_z_xxz_yzzzzz, tr_z_xxz_zzzzz, tr_z_xxz_zzzzzz, tr_z_xz_xxxxxz, tr_z_xz_xxxxyz, tr_z_xz_xxxxzz, tr_z_xz_xxxyyz, tr_z_xz_xxxyzz, tr_z_xz_xxxzzz, tr_z_xz_xxyyyz, tr_z_xz_xxyyzz, tr_z_xz_xxyzzz, tr_z_xz_xxzzzz, tr_z_xz_xyyyyz, tr_z_xz_xyyyzz, tr_z_xz_xyyzzz, tr_z_xz_xyzzzz, tr_z_xz_xzzzzz, tr_z_xz_yyyyyy, tr_z_xz_yyyyyz, tr_z_xz_yyyyzz, tr_z_xz_yyyzzz, tr_z_xz_yyzzzz, tr_z_xz_yzzzzz, tr_z_xz_zzzzzz, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxyy, ts_xxx_xxxyyy, ts_xxx_xxyyyy, ts_xxx_xyyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxz_xxxxxx[i] = ts_xxx_xxxxxx[i] * fe_0 + tr_z_xxx_xxxxxx[i] * pa_z[i];

        tr_z_xxxz_xxxxxy[i] = ts_xxx_xxxxxy[i] * fe_0 + tr_z_xxx_xxxxxy[i] * pa_z[i];

        tr_z_xxxz_xxxxxz[i] = 2.0 * tr_z_xz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxz_xxxxz[i] * fe_0 + tr_z_xxz_xxxxxz[i] * pa_x[i];

        tr_z_xxxz_xxxxyy[i] = ts_xxx_xxxxyy[i] * fe_0 + tr_z_xxx_xxxxyy[i] * pa_z[i];

        tr_z_xxxz_xxxxyz[i] = 2.0 * tr_z_xz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxz_xxxyz[i] * fe_0 + tr_z_xxz_xxxxyz[i] * pa_x[i];

        tr_z_xxxz_xxxxzz[i] = 2.0 * tr_z_xz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxz_xxxzz[i] * fe_0 + tr_z_xxz_xxxxzz[i] * pa_x[i];

        tr_z_xxxz_xxxyyy[i] = ts_xxx_xxxyyy[i] * fe_0 + tr_z_xxx_xxxyyy[i] * pa_z[i];

        tr_z_xxxz_xxxyyz[i] = 2.0 * tr_z_xz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxz_xxyyz[i] * fe_0 + tr_z_xxz_xxxyyz[i] * pa_x[i];

        tr_z_xxxz_xxxyzz[i] = 2.0 * tr_z_xz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxz_xxyzz[i] * fe_0 + tr_z_xxz_xxxyzz[i] * pa_x[i];

        tr_z_xxxz_xxxzzz[i] = 2.0 * tr_z_xz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxz_xxzzz[i] * fe_0 + tr_z_xxz_xxxzzz[i] * pa_x[i];

        tr_z_xxxz_xxyyyy[i] = ts_xxx_xxyyyy[i] * fe_0 + tr_z_xxx_xxyyyy[i] * pa_z[i];

        tr_z_xxxz_xxyyyz[i] = 2.0 * tr_z_xz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxz_xyyyz[i] * fe_0 + tr_z_xxz_xxyyyz[i] * pa_x[i];

        tr_z_xxxz_xxyyzz[i] = 2.0 * tr_z_xz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxz_xyyzz[i] * fe_0 + tr_z_xxz_xxyyzz[i] * pa_x[i];

        tr_z_xxxz_xxyzzz[i] = 2.0 * tr_z_xz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxz_xyzzz[i] * fe_0 + tr_z_xxz_xxyzzz[i] * pa_x[i];

        tr_z_xxxz_xxzzzz[i] = 2.0 * tr_z_xz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxz_xzzzz[i] * fe_0 + tr_z_xxz_xxzzzz[i] * pa_x[i];

        tr_z_xxxz_xyyyyy[i] = ts_xxx_xyyyyy[i] * fe_0 + tr_z_xxx_xyyyyy[i] * pa_z[i];

        tr_z_xxxz_xyyyyz[i] = 2.0 * tr_z_xz_xyyyyz[i] * fe_0 + tr_z_xxz_yyyyz[i] * fe_0 + tr_z_xxz_xyyyyz[i] * pa_x[i];

        tr_z_xxxz_xyyyzz[i] = 2.0 * tr_z_xz_xyyyzz[i] * fe_0 + tr_z_xxz_yyyzz[i] * fe_0 + tr_z_xxz_xyyyzz[i] * pa_x[i];

        tr_z_xxxz_xyyzzz[i] = 2.0 * tr_z_xz_xyyzzz[i] * fe_0 + tr_z_xxz_yyzzz[i] * fe_0 + tr_z_xxz_xyyzzz[i] * pa_x[i];

        tr_z_xxxz_xyzzzz[i] = 2.0 * tr_z_xz_xyzzzz[i] * fe_0 + tr_z_xxz_yzzzz[i] * fe_0 + tr_z_xxz_xyzzzz[i] * pa_x[i];

        tr_z_xxxz_xzzzzz[i] = 2.0 * tr_z_xz_xzzzzz[i] * fe_0 + tr_z_xxz_zzzzz[i] * fe_0 + tr_z_xxz_xzzzzz[i] * pa_x[i];

        tr_z_xxxz_yyyyyy[i] = 2.0 * tr_z_xz_yyyyyy[i] * fe_0 + tr_z_xxz_yyyyyy[i] * pa_x[i];

        tr_z_xxxz_yyyyyz[i] = 2.0 * tr_z_xz_yyyyyz[i] * fe_0 + tr_z_xxz_yyyyyz[i] * pa_x[i];

        tr_z_xxxz_yyyyzz[i] = 2.0 * tr_z_xz_yyyyzz[i] * fe_0 + tr_z_xxz_yyyyzz[i] * pa_x[i];

        tr_z_xxxz_yyyzzz[i] = 2.0 * tr_z_xz_yyyzzz[i] * fe_0 + tr_z_xxz_yyyzzz[i] * pa_x[i];

        tr_z_xxxz_yyzzzz[i] = 2.0 * tr_z_xz_yyzzzz[i] * fe_0 + tr_z_xxz_yyzzzz[i] * pa_x[i];

        tr_z_xxxz_yzzzzz[i] = 2.0 * tr_z_xz_yzzzzz[i] * fe_0 + tr_z_xxz_yzzzzz[i] * pa_x[i];

        tr_z_xxxz_zzzzzz[i] = 2.0 * tr_z_xz_zzzzzz[i] * fe_0 + tr_z_xxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 924-952 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xx_xxxxxx, tr_z_xx_xxxxxz, tr_z_xx_xxxxzz, tr_z_xx_xxxzzz, tr_z_xx_xxzzzz, tr_z_xx_xzzzzz, tr_z_xxy_xxxxxx, tr_z_xxy_xxxxxz, tr_z_xxy_xxxxzz, tr_z_xxy_xxxzzz, tr_z_xxy_xxzzzz, tr_z_xxy_xzzzzz, tr_z_xxyy_xxxxxx, tr_z_xxyy_xxxxxy, tr_z_xxyy_xxxxxz, tr_z_xxyy_xxxxyy, tr_z_xxyy_xxxxyz, tr_z_xxyy_xxxxzz, tr_z_xxyy_xxxyyy, tr_z_xxyy_xxxyyz, tr_z_xxyy_xxxyzz, tr_z_xxyy_xxxzzz, tr_z_xxyy_xxyyyy, tr_z_xxyy_xxyyyz, tr_z_xxyy_xxyyzz, tr_z_xxyy_xxyzzz, tr_z_xxyy_xxzzzz, tr_z_xxyy_xyyyyy, tr_z_xxyy_xyyyyz, tr_z_xxyy_xyyyzz, tr_z_xxyy_xyyzzz, tr_z_xxyy_xyzzzz, tr_z_xxyy_xzzzzz, tr_z_xxyy_yyyyyy, tr_z_xxyy_yyyyyz, tr_z_xxyy_yyyyzz, tr_z_xxyy_yyyzzz, tr_z_xxyy_yyzzzz, tr_z_xxyy_yzzzzz, tr_z_xxyy_zzzzzz, tr_z_xyy_xxxxxy, tr_z_xyy_xxxxy, tr_z_xyy_xxxxyy, tr_z_xyy_xxxxyz, tr_z_xyy_xxxyy, tr_z_xyy_xxxyyy, tr_z_xyy_xxxyyz, tr_z_xyy_xxxyz, tr_z_xyy_xxxyzz, tr_z_xyy_xxyyy, tr_z_xyy_xxyyyy, tr_z_xyy_xxyyyz, tr_z_xyy_xxyyz, tr_z_xyy_xxyyzz, tr_z_xyy_xxyzz, tr_z_xyy_xxyzzz, tr_z_xyy_xyyyy, tr_z_xyy_xyyyyy, tr_z_xyy_xyyyyz, tr_z_xyy_xyyyz, tr_z_xyy_xyyyzz, tr_z_xyy_xyyzz, tr_z_xyy_xyyzzz, tr_z_xyy_xyzzz, tr_z_xyy_xyzzzz, tr_z_xyy_yyyyy, tr_z_xyy_yyyyyy, tr_z_xyy_yyyyyz, tr_z_xyy_yyyyz, tr_z_xyy_yyyyzz, tr_z_xyy_yyyzz, tr_z_xyy_yyyzzz, tr_z_xyy_yyzzz, tr_z_xyy_yyzzzz, tr_z_xyy_yzzzz, tr_z_xyy_yzzzzz, tr_z_xyy_zzzzzz, tr_z_yy_xxxxxy, tr_z_yy_xxxxyy, tr_z_yy_xxxxyz, tr_z_yy_xxxyyy, tr_z_yy_xxxyyz, tr_z_yy_xxxyzz, tr_z_yy_xxyyyy, tr_z_yy_xxyyyz, tr_z_yy_xxyyzz, tr_z_yy_xxyzzz, tr_z_yy_xyyyyy, tr_z_yy_xyyyyz, tr_z_yy_xyyyzz, tr_z_yy_xyyzzz, tr_z_yy_xyzzzz, tr_z_yy_yyyyyy, tr_z_yy_yyyyyz, tr_z_yy_yyyyzz, tr_z_yy_yyyzzz, tr_z_yy_yyzzzz, tr_z_yy_yzzzzz, tr_z_yy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyy_xxxxxx[i] = tr_z_xx_xxxxxx[i] * fe_0 + tr_z_xxy_xxxxxx[i] * pa_y[i];

        tr_z_xxyy_xxxxxy[i] = tr_z_yy_xxxxxy[i] * fe_0 + 5.0 * tr_z_xyy_xxxxy[i] * fe_0 + tr_z_xyy_xxxxxy[i] * pa_x[i];

        tr_z_xxyy_xxxxxz[i] = tr_z_xx_xxxxxz[i] * fe_0 + tr_z_xxy_xxxxxz[i] * pa_y[i];

        tr_z_xxyy_xxxxyy[i] = tr_z_yy_xxxxyy[i] * fe_0 + 4.0 * tr_z_xyy_xxxyy[i] * fe_0 + tr_z_xyy_xxxxyy[i] * pa_x[i];

        tr_z_xxyy_xxxxyz[i] = tr_z_yy_xxxxyz[i] * fe_0 + 4.0 * tr_z_xyy_xxxyz[i] * fe_0 + tr_z_xyy_xxxxyz[i] * pa_x[i];

        tr_z_xxyy_xxxxzz[i] = tr_z_xx_xxxxzz[i] * fe_0 + tr_z_xxy_xxxxzz[i] * pa_y[i];

        tr_z_xxyy_xxxyyy[i] = tr_z_yy_xxxyyy[i] * fe_0 + 3.0 * tr_z_xyy_xxyyy[i] * fe_0 + tr_z_xyy_xxxyyy[i] * pa_x[i];

        tr_z_xxyy_xxxyyz[i] = tr_z_yy_xxxyyz[i] * fe_0 + 3.0 * tr_z_xyy_xxyyz[i] * fe_0 + tr_z_xyy_xxxyyz[i] * pa_x[i];

        tr_z_xxyy_xxxyzz[i] = tr_z_yy_xxxyzz[i] * fe_0 + 3.0 * tr_z_xyy_xxyzz[i] * fe_0 + tr_z_xyy_xxxyzz[i] * pa_x[i];

        tr_z_xxyy_xxxzzz[i] = tr_z_xx_xxxzzz[i] * fe_0 + tr_z_xxy_xxxzzz[i] * pa_y[i];

        tr_z_xxyy_xxyyyy[i] = tr_z_yy_xxyyyy[i] * fe_0 + 2.0 * tr_z_xyy_xyyyy[i] * fe_0 + tr_z_xyy_xxyyyy[i] * pa_x[i];

        tr_z_xxyy_xxyyyz[i] = tr_z_yy_xxyyyz[i] * fe_0 + 2.0 * tr_z_xyy_xyyyz[i] * fe_0 + tr_z_xyy_xxyyyz[i] * pa_x[i];

        tr_z_xxyy_xxyyzz[i] = tr_z_yy_xxyyzz[i] * fe_0 + 2.0 * tr_z_xyy_xyyzz[i] * fe_0 + tr_z_xyy_xxyyzz[i] * pa_x[i];

        tr_z_xxyy_xxyzzz[i] = tr_z_yy_xxyzzz[i] * fe_0 + 2.0 * tr_z_xyy_xyzzz[i] * fe_0 + tr_z_xyy_xxyzzz[i] * pa_x[i];

        tr_z_xxyy_xxzzzz[i] = tr_z_xx_xxzzzz[i] * fe_0 + tr_z_xxy_xxzzzz[i] * pa_y[i];

        tr_z_xxyy_xyyyyy[i] = tr_z_yy_xyyyyy[i] * fe_0 + tr_z_xyy_yyyyy[i] * fe_0 + tr_z_xyy_xyyyyy[i] * pa_x[i];

        tr_z_xxyy_xyyyyz[i] = tr_z_yy_xyyyyz[i] * fe_0 + tr_z_xyy_yyyyz[i] * fe_0 + tr_z_xyy_xyyyyz[i] * pa_x[i];

        tr_z_xxyy_xyyyzz[i] = tr_z_yy_xyyyzz[i] * fe_0 + tr_z_xyy_yyyzz[i] * fe_0 + tr_z_xyy_xyyyzz[i] * pa_x[i];

        tr_z_xxyy_xyyzzz[i] = tr_z_yy_xyyzzz[i] * fe_0 + tr_z_xyy_yyzzz[i] * fe_0 + tr_z_xyy_xyyzzz[i] * pa_x[i];

        tr_z_xxyy_xyzzzz[i] = tr_z_yy_xyzzzz[i] * fe_0 + tr_z_xyy_yzzzz[i] * fe_0 + tr_z_xyy_xyzzzz[i] * pa_x[i];

        tr_z_xxyy_xzzzzz[i] = tr_z_xx_xzzzzz[i] * fe_0 + tr_z_xxy_xzzzzz[i] * pa_y[i];

        tr_z_xxyy_yyyyyy[i] = tr_z_yy_yyyyyy[i] * fe_0 + tr_z_xyy_yyyyyy[i] * pa_x[i];

        tr_z_xxyy_yyyyyz[i] = tr_z_yy_yyyyyz[i] * fe_0 + tr_z_xyy_yyyyyz[i] * pa_x[i];

        tr_z_xxyy_yyyyzz[i] = tr_z_yy_yyyyzz[i] * fe_0 + tr_z_xyy_yyyyzz[i] * pa_x[i];

        tr_z_xxyy_yyyzzz[i] = tr_z_yy_yyyzzz[i] * fe_0 + tr_z_xyy_yyyzzz[i] * pa_x[i];

        tr_z_xxyy_yyzzzz[i] = tr_z_yy_yyzzzz[i] * fe_0 + tr_z_xyy_yyzzzz[i] * pa_x[i];

        tr_z_xxyy_yzzzzz[i] = tr_z_yy_yzzzzz[i] * fe_0 + tr_z_xyy_yzzzzz[i] * pa_x[i];

        tr_z_xxyy_zzzzzz[i] = tr_z_yy_zzzzzz[i] * fe_0 + tr_z_xyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 952-980 components of targeted buffer : GI

    auto tr_z_xxyz_xxxxxx = pbuffer.data(idx_dip_gi + 952);

    auto tr_z_xxyz_xxxxxy = pbuffer.data(idx_dip_gi + 953);

    auto tr_z_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 954);

    auto tr_z_xxyz_xxxxyy = pbuffer.data(idx_dip_gi + 955);

    auto tr_z_xxyz_xxxxyz = pbuffer.data(idx_dip_gi + 956);

    auto tr_z_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 957);

    auto tr_z_xxyz_xxxyyy = pbuffer.data(idx_dip_gi + 958);

    auto tr_z_xxyz_xxxyyz = pbuffer.data(idx_dip_gi + 959);

    auto tr_z_xxyz_xxxyzz = pbuffer.data(idx_dip_gi + 960);

    auto tr_z_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 961);

    auto tr_z_xxyz_xxyyyy = pbuffer.data(idx_dip_gi + 962);

    auto tr_z_xxyz_xxyyyz = pbuffer.data(idx_dip_gi + 963);

    auto tr_z_xxyz_xxyyzz = pbuffer.data(idx_dip_gi + 964);

    auto tr_z_xxyz_xxyzzz = pbuffer.data(idx_dip_gi + 965);

    auto tr_z_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 966);

    auto tr_z_xxyz_xyyyyy = pbuffer.data(idx_dip_gi + 967);

    auto tr_z_xxyz_xyyyyz = pbuffer.data(idx_dip_gi + 968);

    auto tr_z_xxyz_xyyyzz = pbuffer.data(idx_dip_gi + 969);

    auto tr_z_xxyz_xyyzzz = pbuffer.data(idx_dip_gi + 970);

    auto tr_z_xxyz_xyzzzz = pbuffer.data(idx_dip_gi + 971);

    auto tr_z_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 972);

    auto tr_z_xxyz_yyyyyy = pbuffer.data(idx_dip_gi + 973);

    auto tr_z_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 974);

    auto tr_z_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 975);

    auto tr_z_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 976);

    auto tr_z_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 977);

    auto tr_z_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 978);

    auto tr_z_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 979);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyz_xxxxxx, tr_z_xxyz_xxxxxy, tr_z_xxyz_xxxxxz, tr_z_xxyz_xxxxyy, tr_z_xxyz_xxxxyz, tr_z_xxyz_xxxxzz, tr_z_xxyz_xxxyyy, tr_z_xxyz_xxxyyz, tr_z_xxyz_xxxyzz, tr_z_xxyz_xxxzzz, tr_z_xxyz_xxyyyy, tr_z_xxyz_xxyyyz, tr_z_xxyz_xxyyzz, tr_z_xxyz_xxyzzz, tr_z_xxyz_xxzzzz, tr_z_xxyz_xyyyyy, tr_z_xxyz_xyyyyz, tr_z_xxyz_xyyyzz, tr_z_xxyz_xyyzzz, tr_z_xxyz_xyzzzz, tr_z_xxyz_xzzzzz, tr_z_xxyz_yyyyyy, tr_z_xxyz_yyyyyz, tr_z_xxyz_yyyyzz, tr_z_xxyz_yyyzzz, tr_z_xxyz_yyzzzz, tr_z_xxyz_yzzzzz, tr_z_xxyz_zzzzzz, tr_z_xxz_xxxxx, tr_z_xxz_xxxxxx, tr_z_xxz_xxxxxy, tr_z_xxz_xxxxxz, tr_z_xxz_xxxxy, tr_z_xxz_xxxxyy, tr_z_xxz_xxxxyz, tr_z_xxz_xxxxz, tr_z_xxz_xxxxzz, tr_z_xxz_xxxyy, tr_z_xxz_xxxyyy, tr_z_xxz_xxxyyz, tr_z_xxz_xxxyz, tr_z_xxz_xxxyzz, tr_z_xxz_xxxzz, tr_z_xxz_xxxzzz, tr_z_xxz_xxyyy, tr_z_xxz_xxyyyy, tr_z_xxz_xxyyyz, tr_z_xxz_xxyyz, tr_z_xxz_xxyyzz, tr_z_xxz_xxyzz, tr_z_xxz_xxyzzz, tr_z_xxz_xxzzz, tr_z_xxz_xxzzzz, tr_z_xxz_xyyyy, tr_z_xxz_xyyyyy, tr_z_xxz_xyyyyz, tr_z_xxz_xyyyz, tr_z_xxz_xyyyzz, tr_z_xxz_xyyzz, tr_z_xxz_xyyzzz, tr_z_xxz_xyzzz, tr_z_xxz_xyzzzz, tr_z_xxz_xzzzz, tr_z_xxz_xzzzzz, tr_z_xxz_zzzzzz, tr_z_xyz_yyyyyy, tr_z_xyz_yyyyyz, tr_z_xyz_yyyyzz, tr_z_xyz_yyyzzz, tr_z_xyz_yyzzzz, tr_z_xyz_yzzzzz, tr_z_yz_yyyyyy, tr_z_yz_yyyyyz, tr_z_yz_yyyyzz, tr_z_yz_yyyzzz, tr_z_yz_yyzzzz, tr_z_yz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyz_xxxxxx[i] = tr_z_xxz_xxxxxx[i] * pa_y[i];

        tr_z_xxyz_xxxxxy[i] = tr_z_xxz_xxxxx[i] * fe_0 + tr_z_xxz_xxxxxy[i] * pa_y[i];

        tr_z_xxyz_xxxxxz[i] = tr_z_xxz_xxxxxz[i] * pa_y[i];

        tr_z_xxyz_xxxxyy[i] = 2.0 * tr_z_xxz_xxxxy[i] * fe_0 + tr_z_xxz_xxxxyy[i] * pa_y[i];

        tr_z_xxyz_xxxxyz[i] = tr_z_xxz_xxxxz[i] * fe_0 + tr_z_xxz_xxxxyz[i] * pa_y[i];

        tr_z_xxyz_xxxxzz[i] = tr_z_xxz_xxxxzz[i] * pa_y[i];

        tr_z_xxyz_xxxyyy[i] = 3.0 * tr_z_xxz_xxxyy[i] * fe_0 + tr_z_xxz_xxxyyy[i] * pa_y[i];

        tr_z_xxyz_xxxyyz[i] = 2.0 * tr_z_xxz_xxxyz[i] * fe_0 + tr_z_xxz_xxxyyz[i] * pa_y[i];

        tr_z_xxyz_xxxyzz[i] = tr_z_xxz_xxxzz[i] * fe_0 + tr_z_xxz_xxxyzz[i] * pa_y[i];

        tr_z_xxyz_xxxzzz[i] = tr_z_xxz_xxxzzz[i] * pa_y[i];

        tr_z_xxyz_xxyyyy[i] = 4.0 * tr_z_xxz_xxyyy[i] * fe_0 + tr_z_xxz_xxyyyy[i] * pa_y[i];

        tr_z_xxyz_xxyyyz[i] = 3.0 * tr_z_xxz_xxyyz[i] * fe_0 + tr_z_xxz_xxyyyz[i] * pa_y[i];

        tr_z_xxyz_xxyyzz[i] = 2.0 * tr_z_xxz_xxyzz[i] * fe_0 + tr_z_xxz_xxyyzz[i] * pa_y[i];

        tr_z_xxyz_xxyzzz[i] = tr_z_xxz_xxzzz[i] * fe_0 + tr_z_xxz_xxyzzz[i] * pa_y[i];

        tr_z_xxyz_xxzzzz[i] = tr_z_xxz_xxzzzz[i] * pa_y[i];

        tr_z_xxyz_xyyyyy[i] = 5.0 * tr_z_xxz_xyyyy[i] * fe_0 + tr_z_xxz_xyyyyy[i] * pa_y[i];

        tr_z_xxyz_xyyyyz[i] = 4.0 * tr_z_xxz_xyyyz[i] * fe_0 + tr_z_xxz_xyyyyz[i] * pa_y[i];

        tr_z_xxyz_xyyyzz[i] = 3.0 * tr_z_xxz_xyyzz[i] * fe_0 + tr_z_xxz_xyyyzz[i] * pa_y[i];

        tr_z_xxyz_xyyzzz[i] = 2.0 * tr_z_xxz_xyzzz[i] * fe_0 + tr_z_xxz_xyyzzz[i] * pa_y[i];

        tr_z_xxyz_xyzzzz[i] = tr_z_xxz_xzzzz[i] * fe_0 + tr_z_xxz_xyzzzz[i] * pa_y[i];

        tr_z_xxyz_xzzzzz[i] = tr_z_xxz_xzzzzz[i] * pa_y[i];

        tr_z_xxyz_yyyyyy[i] = tr_z_yz_yyyyyy[i] * fe_0 + tr_z_xyz_yyyyyy[i] * pa_x[i];

        tr_z_xxyz_yyyyyz[i] = tr_z_yz_yyyyyz[i] * fe_0 + tr_z_xyz_yyyyyz[i] * pa_x[i];

        tr_z_xxyz_yyyyzz[i] = tr_z_yz_yyyyzz[i] * fe_0 + tr_z_xyz_yyyyzz[i] * pa_x[i];

        tr_z_xxyz_yyyzzz[i] = tr_z_yz_yyyzzz[i] * fe_0 + tr_z_xyz_yyyzzz[i] * pa_x[i];

        tr_z_xxyz_yyzzzz[i] = tr_z_yz_yyzzzz[i] * fe_0 + tr_z_xyz_yyzzzz[i] * pa_x[i];

        tr_z_xxyz_yzzzzz[i] = tr_z_yz_yzzzzz[i] * fe_0 + tr_z_xyz_yzzzzz[i] * pa_x[i];

        tr_z_xxyz_zzzzzz[i] = tr_z_xxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 980-1008 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_z_xxzz_xxxxxx, tr_z_xxzz_xxxxxy, tr_z_xxzz_xxxxxz, tr_z_xxzz_xxxxyy, tr_z_xxzz_xxxxyz, tr_z_xxzz_xxxxzz, tr_z_xxzz_xxxyyy, tr_z_xxzz_xxxyyz, tr_z_xxzz_xxxyzz, tr_z_xxzz_xxxzzz, tr_z_xxzz_xxyyyy, tr_z_xxzz_xxyyyz, tr_z_xxzz_xxyyzz, tr_z_xxzz_xxyzzz, tr_z_xxzz_xxzzzz, tr_z_xxzz_xyyyyy, tr_z_xxzz_xyyyyz, tr_z_xxzz_xyyyzz, tr_z_xxzz_xyyzzz, tr_z_xxzz_xyzzzz, tr_z_xxzz_xzzzzz, tr_z_xxzz_yyyyyy, tr_z_xxzz_yyyyyz, tr_z_xxzz_yyyyzz, tr_z_xxzz_yyyzzz, tr_z_xxzz_yyzzzz, tr_z_xxzz_yzzzzz, tr_z_xxzz_zzzzzz, tr_z_xzz_xxxxx, tr_z_xzz_xxxxxx, tr_z_xzz_xxxxxy, tr_z_xzz_xxxxxz, tr_z_xzz_xxxxy, tr_z_xzz_xxxxyy, tr_z_xzz_xxxxyz, tr_z_xzz_xxxxz, tr_z_xzz_xxxxzz, tr_z_xzz_xxxyy, tr_z_xzz_xxxyyy, tr_z_xzz_xxxyyz, tr_z_xzz_xxxyz, tr_z_xzz_xxxyzz, tr_z_xzz_xxxzz, tr_z_xzz_xxxzzz, tr_z_xzz_xxyyy, tr_z_xzz_xxyyyy, tr_z_xzz_xxyyyz, tr_z_xzz_xxyyz, tr_z_xzz_xxyyzz, tr_z_xzz_xxyzz, tr_z_xzz_xxyzzz, tr_z_xzz_xxzzz, tr_z_xzz_xxzzzz, tr_z_xzz_xyyyy, tr_z_xzz_xyyyyy, tr_z_xzz_xyyyyz, tr_z_xzz_xyyyz, tr_z_xzz_xyyyzz, tr_z_xzz_xyyzz, tr_z_xzz_xyyzzz, tr_z_xzz_xyzzz, tr_z_xzz_xyzzzz, tr_z_xzz_xzzzz, tr_z_xzz_xzzzzz, tr_z_xzz_yyyyy, tr_z_xzz_yyyyyy, tr_z_xzz_yyyyyz, tr_z_xzz_yyyyz, tr_z_xzz_yyyyzz, tr_z_xzz_yyyzz, tr_z_xzz_yyyzzz, tr_z_xzz_yyzzz, tr_z_xzz_yyzzzz, tr_z_xzz_yzzzz, tr_z_xzz_yzzzzz, tr_z_xzz_zzzzz, tr_z_xzz_zzzzzz, tr_z_zz_xxxxxx, tr_z_zz_xxxxxy, tr_z_zz_xxxxxz, tr_z_zz_xxxxyy, tr_z_zz_xxxxyz, tr_z_zz_xxxxzz, tr_z_zz_xxxyyy, tr_z_zz_xxxyyz, tr_z_zz_xxxyzz, tr_z_zz_xxxzzz, tr_z_zz_xxyyyy, tr_z_zz_xxyyyz, tr_z_zz_xxyyzz, tr_z_zz_xxyzzz, tr_z_zz_xxzzzz, tr_z_zz_xyyyyy, tr_z_zz_xyyyyz, tr_z_zz_xyyyzz, tr_z_zz_xyyzzz, tr_z_zz_xyzzzz, tr_z_zz_xzzzzz, tr_z_zz_yyyyyy, tr_z_zz_yyyyyz, tr_z_zz_yyyyzz, tr_z_zz_yyyzzz, tr_z_zz_yyzzzz, tr_z_zz_yzzzzz, tr_z_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzz_xxxxxx[i] = tr_z_zz_xxxxxx[i] * fe_0 + 6.0 * tr_z_xzz_xxxxx[i] * fe_0 + tr_z_xzz_xxxxxx[i] * pa_x[i];

        tr_z_xxzz_xxxxxy[i] = tr_z_zz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xzz_xxxxy[i] * fe_0 + tr_z_xzz_xxxxxy[i] * pa_x[i];

        tr_z_xxzz_xxxxxz[i] = tr_z_zz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xzz_xxxxz[i] * fe_0 + tr_z_xzz_xxxxxz[i] * pa_x[i];

        tr_z_xxzz_xxxxyy[i] = tr_z_zz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xzz_xxxyy[i] * fe_0 + tr_z_xzz_xxxxyy[i] * pa_x[i];

        tr_z_xxzz_xxxxyz[i] = tr_z_zz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xzz_xxxyz[i] * fe_0 + tr_z_xzz_xxxxyz[i] * pa_x[i];

        tr_z_xxzz_xxxxzz[i] = tr_z_zz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xzz_xxxzz[i] * fe_0 + tr_z_xzz_xxxxzz[i] * pa_x[i];

        tr_z_xxzz_xxxyyy[i] = tr_z_zz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xzz_xxyyy[i] * fe_0 + tr_z_xzz_xxxyyy[i] * pa_x[i];

        tr_z_xxzz_xxxyyz[i] = tr_z_zz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xzz_xxyyz[i] * fe_0 + tr_z_xzz_xxxyyz[i] * pa_x[i];

        tr_z_xxzz_xxxyzz[i] = tr_z_zz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xzz_xxyzz[i] * fe_0 + tr_z_xzz_xxxyzz[i] * pa_x[i];

        tr_z_xxzz_xxxzzz[i] = tr_z_zz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xzz_xxzzz[i] * fe_0 + tr_z_xzz_xxxzzz[i] * pa_x[i];

        tr_z_xxzz_xxyyyy[i] = tr_z_zz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xzz_xyyyy[i] * fe_0 + tr_z_xzz_xxyyyy[i] * pa_x[i];

        tr_z_xxzz_xxyyyz[i] = tr_z_zz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xzz_xyyyz[i] * fe_0 + tr_z_xzz_xxyyyz[i] * pa_x[i];

        tr_z_xxzz_xxyyzz[i] = tr_z_zz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xzz_xyyzz[i] * fe_0 + tr_z_xzz_xxyyzz[i] * pa_x[i];

        tr_z_xxzz_xxyzzz[i] = tr_z_zz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xzz_xyzzz[i] * fe_0 + tr_z_xzz_xxyzzz[i] * pa_x[i];

        tr_z_xxzz_xxzzzz[i] = tr_z_zz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xzz_xzzzz[i] * fe_0 + tr_z_xzz_xxzzzz[i] * pa_x[i];

        tr_z_xxzz_xyyyyy[i] = tr_z_zz_xyyyyy[i] * fe_0 + tr_z_xzz_yyyyy[i] * fe_0 + tr_z_xzz_xyyyyy[i] * pa_x[i];

        tr_z_xxzz_xyyyyz[i] = tr_z_zz_xyyyyz[i] * fe_0 + tr_z_xzz_yyyyz[i] * fe_0 + tr_z_xzz_xyyyyz[i] * pa_x[i];

        tr_z_xxzz_xyyyzz[i] = tr_z_zz_xyyyzz[i] * fe_0 + tr_z_xzz_yyyzz[i] * fe_0 + tr_z_xzz_xyyyzz[i] * pa_x[i];

        tr_z_xxzz_xyyzzz[i] = tr_z_zz_xyyzzz[i] * fe_0 + tr_z_xzz_yyzzz[i] * fe_0 + tr_z_xzz_xyyzzz[i] * pa_x[i];

        tr_z_xxzz_xyzzzz[i] = tr_z_zz_xyzzzz[i] * fe_0 + tr_z_xzz_yzzzz[i] * fe_0 + tr_z_xzz_xyzzzz[i] * pa_x[i];

        tr_z_xxzz_xzzzzz[i] = tr_z_zz_xzzzzz[i] * fe_0 + tr_z_xzz_zzzzz[i] * fe_0 + tr_z_xzz_xzzzzz[i] * pa_x[i];

        tr_z_xxzz_yyyyyy[i] = tr_z_zz_yyyyyy[i] * fe_0 + tr_z_xzz_yyyyyy[i] * pa_x[i];

        tr_z_xxzz_yyyyyz[i] = tr_z_zz_yyyyyz[i] * fe_0 + tr_z_xzz_yyyyyz[i] * pa_x[i];

        tr_z_xxzz_yyyyzz[i] = tr_z_zz_yyyyzz[i] * fe_0 + tr_z_xzz_yyyyzz[i] * pa_x[i];

        tr_z_xxzz_yyyzzz[i] = tr_z_zz_yyyzzz[i] * fe_0 + tr_z_xzz_yyyzzz[i] * pa_x[i];

        tr_z_xxzz_yyzzzz[i] = tr_z_zz_yyzzzz[i] * fe_0 + tr_z_xzz_yyzzzz[i] * pa_x[i];

        tr_z_xxzz_yzzzzz[i] = tr_z_zz_yzzzzz[i] * fe_0 + tr_z_xzz_yzzzzz[i] * pa_x[i];

        tr_z_xxzz_zzzzzz[i] = tr_z_zz_zzzzzz[i] * fe_0 + tr_z_xzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1008-1036 components of targeted buffer : GI

    auto tr_z_xyyy_xxxxxx = pbuffer.data(idx_dip_gi + 1008);

    auto tr_z_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 1009);

    auto tr_z_xyyy_xxxxxz = pbuffer.data(idx_dip_gi + 1010);

    auto tr_z_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 1011);

    auto tr_z_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 1012);

    auto tr_z_xyyy_xxxxzz = pbuffer.data(idx_dip_gi + 1013);

    auto tr_z_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 1014);

    auto tr_z_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 1015);

    auto tr_z_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 1016);

    auto tr_z_xyyy_xxxzzz = pbuffer.data(idx_dip_gi + 1017);

    auto tr_z_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 1018);

    auto tr_z_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 1019);

    auto tr_z_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 1020);

    auto tr_z_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 1021);

    auto tr_z_xyyy_xxzzzz = pbuffer.data(idx_dip_gi + 1022);

    auto tr_z_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 1023);

    auto tr_z_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 1024);

    auto tr_z_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 1025);

    auto tr_z_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 1026);

    auto tr_z_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 1027);

    auto tr_z_xyyy_xzzzzz = pbuffer.data(idx_dip_gi + 1028);

    auto tr_z_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 1029);

    auto tr_z_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 1030);

    auto tr_z_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 1031);

    auto tr_z_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 1032);

    auto tr_z_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 1033);

    auto tr_z_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 1034);

    auto tr_z_xyyy_zzzzzz = pbuffer.data(idx_dip_gi + 1035);

    #pragma omp simd aligned(pa_x, tr_z_xyyy_xxxxxx, tr_z_xyyy_xxxxxy, tr_z_xyyy_xxxxxz, tr_z_xyyy_xxxxyy, tr_z_xyyy_xxxxyz, tr_z_xyyy_xxxxzz, tr_z_xyyy_xxxyyy, tr_z_xyyy_xxxyyz, tr_z_xyyy_xxxyzz, tr_z_xyyy_xxxzzz, tr_z_xyyy_xxyyyy, tr_z_xyyy_xxyyyz, tr_z_xyyy_xxyyzz, tr_z_xyyy_xxyzzz, tr_z_xyyy_xxzzzz, tr_z_xyyy_xyyyyy, tr_z_xyyy_xyyyyz, tr_z_xyyy_xyyyzz, tr_z_xyyy_xyyzzz, tr_z_xyyy_xyzzzz, tr_z_xyyy_xzzzzz, tr_z_xyyy_yyyyyy, tr_z_xyyy_yyyyyz, tr_z_xyyy_yyyyzz, tr_z_xyyy_yyyzzz, tr_z_xyyy_yyzzzz, tr_z_xyyy_yzzzzz, tr_z_xyyy_zzzzzz, tr_z_yyy_xxxxx, tr_z_yyy_xxxxxx, tr_z_yyy_xxxxxy, tr_z_yyy_xxxxxz, tr_z_yyy_xxxxy, tr_z_yyy_xxxxyy, tr_z_yyy_xxxxyz, tr_z_yyy_xxxxz, tr_z_yyy_xxxxzz, tr_z_yyy_xxxyy, tr_z_yyy_xxxyyy, tr_z_yyy_xxxyyz, tr_z_yyy_xxxyz, tr_z_yyy_xxxyzz, tr_z_yyy_xxxzz, tr_z_yyy_xxxzzz, tr_z_yyy_xxyyy, tr_z_yyy_xxyyyy, tr_z_yyy_xxyyyz, tr_z_yyy_xxyyz, tr_z_yyy_xxyyzz, tr_z_yyy_xxyzz, tr_z_yyy_xxyzzz, tr_z_yyy_xxzzz, tr_z_yyy_xxzzzz, tr_z_yyy_xyyyy, tr_z_yyy_xyyyyy, tr_z_yyy_xyyyyz, tr_z_yyy_xyyyz, tr_z_yyy_xyyyzz, tr_z_yyy_xyyzz, tr_z_yyy_xyyzzz, tr_z_yyy_xyzzz, tr_z_yyy_xyzzzz, tr_z_yyy_xzzzz, tr_z_yyy_xzzzzz, tr_z_yyy_yyyyy, tr_z_yyy_yyyyyy, tr_z_yyy_yyyyyz, tr_z_yyy_yyyyz, tr_z_yyy_yyyyzz, tr_z_yyy_yyyzz, tr_z_yyy_yyyzzz, tr_z_yyy_yyzzz, tr_z_yyy_yyzzzz, tr_z_yyy_yzzzz, tr_z_yyy_yzzzzz, tr_z_yyy_zzzzz, tr_z_yyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyy_xxxxxx[i] = 6.0 * tr_z_yyy_xxxxx[i] * fe_0 + tr_z_yyy_xxxxxx[i] * pa_x[i];

        tr_z_xyyy_xxxxxy[i] = 5.0 * tr_z_yyy_xxxxy[i] * fe_0 + tr_z_yyy_xxxxxy[i] * pa_x[i];

        tr_z_xyyy_xxxxxz[i] = 5.0 * tr_z_yyy_xxxxz[i] * fe_0 + tr_z_yyy_xxxxxz[i] * pa_x[i];

        tr_z_xyyy_xxxxyy[i] = 4.0 * tr_z_yyy_xxxyy[i] * fe_0 + tr_z_yyy_xxxxyy[i] * pa_x[i];

        tr_z_xyyy_xxxxyz[i] = 4.0 * tr_z_yyy_xxxyz[i] * fe_0 + tr_z_yyy_xxxxyz[i] * pa_x[i];

        tr_z_xyyy_xxxxzz[i] = 4.0 * tr_z_yyy_xxxzz[i] * fe_0 + tr_z_yyy_xxxxzz[i] * pa_x[i];

        tr_z_xyyy_xxxyyy[i] = 3.0 * tr_z_yyy_xxyyy[i] * fe_0 + tr_z_yyy_xxxyyy[i] * pa_x[i];

        tr_z_xyyy_xxxyyz[i] = 3.0 * tr_z_yyy_xxyyz[i] * fe_0 + tr_z_yyy_xxxyyz[i] * pa_x[i];

        tr_z_xyyy_xxxyzz[i] = 3.0 * tr_z_yyy_xxyzz[i] * fe_0 + tr_z_yyy_xxxyzz[i] * pa_x[i];

        tr_z_xyyy_xxxzzz[i] = 3.0 * tr_z_yyy_xxzzz[i] * fe_0 + tr_z_yyy_xxxzzz[i] * pa_x[i];

        tr_z_xyyy_xxyyyy[i] = 2.0 * tr_z_yyy_xyyyy[i] * fe_0 + tr_z_yyy_xxyyyy[i] * pa_x[i];

        tr_z_xyyy_xxyyyz[i] = 2.0 * tr_z_yyy_xyyyz[i] * fe_0 + tr_z_yyy_xxyyyz[i] * pa_x[i];

        tr_z_xyyy_xxyyzz[i] = 2.0 * tr_z_yyy_xyyzz[i] * fe_0 + tr_z_yyy_xxyyzz[i] * pa_x[i];

        tr_z_xyyy_xxyzzz[i] = 2.0 * tr_z_yyy_xyzzz[i] * fe_0 + tr_z_yyy_xxyzzz[i] * pa_x[i];

        tr_z_xyyy_xxzzzz[i] = 2.0 * tr_z_yyy_xzzzz[i] * fe_0 + tr_z_yyy_xxzzzz[i] * pa_x[i];

        tr_z_xyyy_xyyyyy[i] = tr_z_yyy_yyyyy[i] * fe_0 + tr_z_yyy_xyyyyy[i] * pa_x[i];

        tr_z_xyyy_xyyyyz[i] = tr_z_yyy_yyyyz[i] * fe_0 + tr_z_yyy_xyyyyz[i] * pa_x[i];

        tr_z_xyyy_xyyyzz[i] = tr_z_yyy_yyyzz[i] * fe_0 + tr_z_yyy_xyyyzz[i] * pa_x[i];

        tr_z_xyyy_xyyzzz[i] = tr_z_yyy_yyzzz[i] * fe_0 + tr_z_yyy_xyyzzz[i] * pa_x[i];

        tr_z_xyyy_xyzzzz[i] = tr_z_yyy_yzzzz[i] * fe_0 + tr_z_yyy_xyzzzz[i] * pa_x[i];

        tr_z_xyyy_xzzzzz[i] = tr_z_yyy_zzzzz[i] * fe_0 + tr_z_yyy_xzzzzz[i] * pa_x[i];

        tr_z_xyyy_yyyyyy[i] = tr_z_yyy_yyyyyy[i] * pa_x[i];

        tr_z_xyyy_yyyyyz[i] = tr_z_yyy_yyyyyz[i] * pa_x[i];

        tr_z_xyyy_yyyyzz[i] = tr_z_yyy_yyyyzz[i] * pa_x[i];

        tr_z_xyyy_yyyzzz[i] = tr_z_yyy_yyyzzz[i] * pa_x[i];

        tr_z_xyyy_yyzzzz[i] = tr_z_yyy_yyzzzz[i] * pa_x[i];

        tr_z_xyyy_yzzzzz[i] = tr_z_yyy_yzzzzz[i] * pa_x[i];

        tr_z_xyyy_zzzzzz[i] = tr_z_yyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1036-1064 components of targeted buffer : GI

    auto tr_z_xyyz_xxxxxx = pbuffer.data(idx_dip_gi + 1036);

    auto tr_z_xyyz_xxxxxy = pbuffer.data(idx_dip_gi + 1037);

    auto tr_z_xyyz_xxxxxz = pbuffer.data(idx_dip_gi + 1038);

    auto tr_z_xyyz_xxxxyy = pbuffer.data(idx_dip_gi + 1039);

    auto tr_z_xyyz_xxxxyz = pbuffer.data(idx_dip_gi + 1040);

    auto tr_z_xyyz_xxxxzz = pbuffer.data(idx_dip_gi + 1041);

    auto tr_z_xyyz_xxxyyy = pbuffer.data(idx_dip_gi + 1042);

    auto tr_z_xyyz_xxxyyz = pbuffer.data(idx_dip_gi + 1043);

    auto tr_z_xyyz_xxxyzz = pbuffer.data(idx_dip_gi + 1044);

    auto tr_z_xyyz_xxxzzz = pbuffer.data(idx_dip_gi + 1045);

    auto tr_z_xyyz_xxyyyy = pbuffer.data(idx_dip_gi + 1046);

    auto tr_z_xyyz_xxyyyz = pbuffer.data(idx_dip_gi + 1047);

    auto tr_z_xyyz_xxyyzz = pbuffer.data(idx_dip_gi + 1048);

    auto tr_z_xyyz_xxyzzz = pbuffer.data(idx_dip_gi + 1049);

    auto tr_z_xyyz_xxzzzz = pbuffer.data(idx_dip_gi + 1050);

    auto tr_z_xyyz_xyyyyy = pbuffer.data(idx_dip_gi + 1051);

    auto tr_z_xyyz_xyyyyz = pbuffer.data(idx_dip_gi + 1052);

    auto tr_z_xyyz_xyyyzz = pbuffer.data(idx_dip_gi + 1053);

    auto tr_z_xyyz_xyyzzz = pbuffer.data(idx_dip_gi + 1054);

    auto tr_z_xyyz_xyzzzz = pbuffer.data(idx_dip_gi + 1055);

    auto tr_z_xyyz_xzzzzz = pbuffer.data(idx_dip_gi + 1056);

    auto tr_z_xyyz_yyyyyy = pbuffer.data(idx_dip_gi + 1057);

    auto tr_z_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 1058);

    auto tr_z_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 1059);

    auto tr_z_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 1060);

    auto tr_z_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 1061);

    auto tr_z_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 1062);

    auto tr_z_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 1063);

    #pragma omp simd aligned(pa_x, tr_z_xyyz_xxxxxx, tr_z_xyyz_xxxxxy, tr_z_xyyz_xxxxxz, tr_z_xyyz_xxxxyy, tr_z_xyyz_xxxxyz, tr_z_xyyz_xxxxzz, tr_z_xyyz_xxxyyy, tr_z_xyyz_xxxyyz, tr_z_xyyz_xxxyzz, tr_z_xyyz_xxxzzz, tr_z_xyyz_xxyyyy, tr_z_xyyz_xxyyyz, tr_z_xyyz_xxyyzz, tr_z_xyyz_xxyzzz, tr_z_xyyz_xxzzzz, tr_z_xyyz_xyyyyy, tr_z_xyyz_xyyyyz, tr_z_xyyz_xyyyzz, tr_z_xyyz_xyyzzz, tr_z_xyyz_xyzzzz, tr_z_xyyz_xzzzzz, tr_z_xyyz_yyyyyy, tr_z_xyyz_yyyyyz, tr_z_xyyz_yyyyzz, tr_z_xyyz_yyyzzz, tr_z_xyyz_yyzzzz, tr_z_xyyz_yzzzzz, tr_z_xyyz_zzzzzz, tr_z_yyz_xxxxx, tr_z_yyz_xxxxxx, tr_z_yyz_xxxxxy, tr_z_yyz_xxxxxz, tr_z_yyz_xxxxy, tr_z_yyz_xxxxyy, tr_z_yyz_xxxxyz, tr_z_yyz_xxxxz, tr_z_yyz_xxxxzz, tr_z_yyz_xxxyy, tr_z_yyz_xxxyyy, tr_z_yyz_xxxyyz, tr_z_yyz_xxxyz, tr_z_yyz_xxxyzz, tr_z_yyz_xxxzz, tr_z_yyz_xxxzzz, tr_z_yyz_xxyyy, tr_z_yyz_xxyyyy, tr_z_yyz_xxyyyz, tr_z_yyz_xxyyz, tr_z_yyz_xxyyzz, tr_z_yyz_xxyzz, tr_z_yyz_xxyzzz, tr_z_yyz_xxzzz, tr_z_yyz_xxzzzz, tr_z_yyz_xyyyy, tr_z_yyz_xyyyyy, tr_z_yyz_xyyyyz, tr_z_yyz_xyyyz, tr_z_yyz_xyyyzz, tr_z_yyz_xyyzz, tr_z_yyz_xyyzzz, tr_z_yyz_xyzzz, tr_z_yyz_xyzzzz, tr_z_yyz_xzzzz, tr_z_yyz_xzzzzz, tr_z_yyz_yyyyy, tr_z_yyz_yyyyyy, tr_z_yyz_yyyyyz, tr_z_yyz_yyyyz, tr_z_yyz_yyyyzz, tr_z_yyz_yyyzz, tr_z_yyz_yyyzzz, tr_z_yyz_yyzzz, tr_z_yyz_yyzzzz, tr_z_yyz_yzzzz, tr_z_yyz_yzzzzz, tr_z_yyz_zzzzz, tr_z_yyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyz_xxxxxx[i] = 6.0 * tr_z_yyz_xxxxx[i] * fe_0 + tr_z_yyz_xxxxxx[i] * pa_x[i];

        tr_z_xyyz_xxxxxy[i] = 5.0 * tr_z_yyz_xxxxy[i] * fe_0 + tr_z_yyz_xxxxxy[i] * pa_x[i];

        tr_z_xyyz_xxxxxz[i] = 5.0 * tr_z_yyz_xxxxz[i] * fe_0 + tr_z_yyz_xxxxxz[i] * pa_x[i];

        tr_z_xyyz_xxxxyy[i] = 4.0 * tr_z_yyz_xxxyy[i] * fe_0 + tr_z_yyz_xxxxyy[i] * pa_x[i];

        tr_z_xyyz_xxxxyz[i] = 4.0 * tr_z_yyz_xxxyz[i] * fe_0 + tr_z_yyz_xxxxyz[i] * pa_x[i];

        tr_z_xyyz_xxxxzz[i] = 4.0 * tr_z_yyz_xxxzz[i] * fe_0 + tr_z_yyz_xxxxzz[i] * pa_x[i];

        tr_z_xyyz_xxxyyy[i] = 3.0 * tr_z_yyz_xxyyy[i] * fe_0 + tr_z_yyz_xxxyyy[i] * pa_x[i];

        tr_z_xyyz_xxxyyz[i] = 3.0 * tr_z_yyz_xxyyz[i] * fe_0 + tr_z_yyz_xxxyyz[i] * pa_x[i];

        tr_z_xyyz_xxxyzz[i] = 3.0 * tr_z_yyz_xxyzz[i] * fe_0 + tr_z_yyz_xxxyzz[i] * pa_x[i];

        tr_z_xyyz_xxxzzz[i] = 3.0 * tr_z_yyz_xxzzz[i] * fe_0 + tr_z_yyz_xxxzzz[i] * pa_x[i];

        tr_z_xyyz_xxyyyy[i] = 2.0 * tr_z_yyz_xyyyy[i] * fe_0 + tr_z_yyz_xxyyyy[i] * pa_x[i];

        tr_z_xyyz_xxyyyz[i] = 2.0 * tr_z_yyz_xyyyz[i] * fe_0 + tr_z_yyz_xxyyyz[i] * pa_x[i];

        tr_z_xyyz_xxyyzz[i] = 2.0 * tr_z_yyz_xyyzz[i] * fe_0 + tr_z_yyz_xxyyzz[i] * pa_x[i];

        tr_z_xyyz_xxyzzz[i] = 2.0 * tr_z_yyz_xyzzz[i] * fe_0 + tr_z_yyz_xxyzzz[i] * pa_x[i];

        tr_z_xyyz_xxzzzz[i] = 2.0 * tr_z_yyz_xzzzz[i] * fe_0 + tr_z_yyz_xxzzzz[i] * pa_x[i];

        tr_z_xyyz_xyyyyy[i] = tr_z_yyz_yyyyy[i] * fe_0 + tr_z_yyz_xyyyyy[i] * pa_x[i];

        tr_z_xyyz_xyyyyz[i] = tr_z_yyz_yyyyz[i] * fe_0 + tr_z_yyz_xyyyyz[i] * pa_x[i];

        tr_z_xyyz_xyyyzz[i] = tr_z_yyz_yyyzz[i] * fe_0 + tr_z_yyz_xyyyzz[i] * pa_x[i];

        tr_z_xyyz_xyyzzz[i] = tr_z_yyz_yyzzz[i] * fe_0 + tr_z_yyz_xyyzzz[i] * pa_x[i];

        tr_z_xyyz_xyzzzz[i] = tr_z_yyz_yzzzz[i] * fe_0 + tr_z_yyz_xyzzzz[i] * pa_x[i];

        tr_z_xyyz_xzzzzz[i] = tr_z_yyz_zzzzz[i] * fe_0 + tr_z_yyz_xzzzzz[i] * pa_x[i];

        tr_z_xyyz_yyyyyy[i] = tr_z_yyz_yyyyyy[i] * pa_x[i];

        tr_z_xyyz_yyyyyz[i] = tr_z_yyz_yyyyyz[i] * pa_x[i];

        tr_z_xyyz_yyyyzz[i] = tr_z_yyz_yyyyzz[i] * pa_x[i];

        tr_z_xyyz_yyyzzz[i] = tr_z_yyz_yyyzzz[i] * pa_x[i];

        tr_z_xyyz_yyzzzz[i] = tr_z_yyz_yyzzzz[i] * pa_x[i];

        tr_z_xyyz_yzzzzz[i] = tr_z_yyz_yzzzzz[i] * pa_x[i];

        tr_z_xyyz_zzzzzz[i] = tr_z_yyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1064-1092 components of targeted buffer : GI

    auto tr_z_xyzz_xxxxxx = pbuffer.data(idx_dip_gi + 1064);

    auto tr_z_xyzz_xxxxxy = pbuffer.data(idx_dip_gi + 1065);

    auto tr_z_xyzz_xxxxxz = pbuffer.data(idx_dip_gi + 1066);

    auto tr_z_xyzz_xxxxyy = pbuffer.data(idx_dip_gi + 1067);

    auto tr_z_xyzz_xxxxyz = pbuffer.data(idx_dip_gi + 1068);

    auto tr_z_xyzz_xxxxzz = pbuffer.data(idx_dip_gi + 1069);

    auto tr_z_xyzz_xxxyyy = pbuffer.data(idx_dip_gi + 1070);

    auto tr_z_xyzz_xxxyyz = pbuffer.data(idx_dip_gi + 1071);

    auto tr_z_xyzz_xxxyzz = pbuffer.data(idx_dip_gi + 1072);

    auto tr_z_xyzz_xxxzzz = pbuffer.data(idx_dip_gi + 1073);

    auto tr_z_xyzz_xxyyyy = pbuffer.data(idx_dip_gi + 1074);

    auto tr_z_xyzz_xxyyyz = pbuffer.data(idx_dip_gi + 1075);

    auto tr_z_xyzz_xxyyzz = pbuffer.data(idx_dip_gi + 1076);

    auto tr_z_xyzz_xxyzzz = pbuffer.data(idx_dip_gi + 1077);

    auto tr_z_xyzz_xxzzzz = pbuffer.data(idx_dip_gi + 1078);

    auto tr_z_xyzz_xyyyyy = pbuffer.data(idx_dip_gi + 1079);

    auto tr_z_xyzz_xyyyyz = pbuffer.data(idx_dip_gi + 1080);

    auto tr_z_xyzz_xyyyzz = pbuffer.data(idx_dip_gi + 1081);

    auto tr_z_xyzz_xyyzzz = pbuffer.data(idx_dip_gi + 1082);

    auto tr_z_xyzz_xyzzzz = pbuffer.data(idx_dip_gi + 1083);

    auto tr_z_xyzz_xzzzzz = pbuffer.data(idx_dip_gi + 1084);

    auto tr_z_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 1085);

    auto tr_z_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 1086);

    auto tr_z_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 1087);

    auto tr_z_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 1088);

    auto tr_z_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 1089);

    auto tr_z_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 1090);

    auto tr_z_xyzz_zzzzzz = pbuffer.data(idx_dip_gi + 1091);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyzz_xxxxxx, tr_z_xyzz_xxxxxy, tr_z_xyzz_xxxxxz, tr_z_xyzz_xxxxyy, tr_z_xyzz_xxxxyz, tr_z_xyzz_xxxxzz, tr_z_xyzz_xxxyyy, tr_z_xyzz_xxxyyz, tr_z_xyzz_xxxyzz, tr_z_xyzz_xxxzzz, tr_z_xyzz_xxyyyy, tr_z_xyzz_xxyyyz, tr_z_xyzz_xxyyzz, tr_z_xyzz_xxyzzz, tr_z_xyzz_xxzzzz, tr_z_xyzz_xyyyyy, tr_z_xyzz_xyyyyz, tr_z_xyzz_xyyyzz, tr_z_xyzz_xyyzzz, tr_z_xyzz_xyzzzz, tr_z_xyzz_xzzzzz, tr_z_xyzz_yyyyyy, tr_z_xyzz_yyyyyz, tr_z_xyzz_yyyyzz, tr_z_xyzz_yyyzzz, tr_z_xyzz_yyzzzz, tr_z_xyzz_yzzzzz, tr_z_xyzz_zzzzzz, tr_z_xzz_xxxxxx, tr_z_xzz_xxxxxz, tr_z_xzz_xxxxzz, tr_z_xzz_xxxzzz, tr_z_xzz_xxzzzz, tr_z_xzz_xzzzzz, tr_z_yzz_xxxxxy, tr_z_yzz_xxxxy, tr_z_yzz_xxxxyy, tr_z_yzz_xxxxyz, tr_z_yzz_xxxyy, tr_z_yzz_xxxyyy, tr_z_yzz_xxxyyz, tr_z_yzz_xxxyz, tr_z_yzz_xxxyzz, tr_z_yzz_xxyyy, tr_z_yzz_xxyyyy, tr_z_yzz_xxyyyz, tr_z_yzz_xxyyz, tr_z_yzz_xxyyzz, tr_z_yzz_xxyzz, tr_z_yzz_xxyzzz, tr_z_yzz_xyyyy, tr_z_yzz_xyyyyy, tr_z_yzz_xyyyyz, tr_z_yzz_xyyyz, tr_z_yzz_xyyyzz, tr_z_yzz_xyyzz, tr_z_yzz_xyyzzz, tr_z_yzz_xyzzz, tr_z_yzz_xyzzzz, tr_z_yzz_yyyyy, tr_z_yzz_yyyyyy, tr_z_yzz_yyyyyz, tr_z_yzz_yyyyz, tr_z_yzz_yyyyzz, tr_z_yzz_yyyzz, tr_z_yzz_yyyzzz, tr_z_yzz_yyzzz, tr_z_yzz_yyzzzz, tr_z_yzz_yzzzz, tr_z_yzz_yzzzzz, tr_z_yzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzz_xxxxxx[i] = tr_z_xzz_xxxxxx[i] * pa_y[i];

        tr_z_xyzz_xxxxxy[i] = 5.0 * tr_z_yzz_xxxxy[i] * fe_0 + tr_z_yzz_xxxxxy[i] * pa_x[i];

        tr_z_xyzz_xxxxxz[i] = tr_z_xzz_xxxxxz[i] * pa_y[i];

        tr_z_xyzz_xxxxyy[i] = 4.0 * tr_z_yzz_xxxyy[i] * fe_0 + tr_z_yzz_xxxxyy[i] * pa_x[i];

        tr_z_xyzz_xxxxyz[i] = 4.0 * tr_z_yzz_xxxyz[i] * fe_0 + tr_z_yzz_xxxxyz[i] * pa_x[i];

        tr_z_xyzz_xxxxzz[i] = tr_z_xzz_xxxxzz[i] * pa_y[i];

        tr_z_xyzz_xxxyyy[i] = 3.0 * tr_z_yzz_xxyyy[i] * fe_0 + tr_z_yzz_xxxyyy[i] * pa_x[i];

        tr_z_xyzz_xxxyyz[i] = 3.0 * tr_z_yzz_xxyyz[i] * fe_0 + tr_z_yzz_xxxyyz[i] * pa_x[i];

        tr_z_xyzz_xxxyzz[i] = 3.0 * tr_z_yzz_xxyzz[i] * fe_0 + tr_z_yzz_xxxyzz[i] * pa_x[i];

        tr_z_xyzz_xxxzzz[i] = tr_z_xzz_xxxzzz[i] * pa_y[i];

        tr_z_xyzz_xxyyyy[i] = 2.0 * tr_z_yzz_xyyyy[i] * fe_0 + tr_z_yzz_xxyyyy[i] * pa_x[i];

        tr_z_xyzz_xxyyyz[i] = 2.0 * tr_z_yzz_xyyyz[i] * fe_0 + tr_z_yzz_xxyyyz[i] * pa_x[i];

        tr_z_xyzz_xxyyzz[i] = 2.0 * tr_z_yzz_xyyzz[i] * fe_0 + tr_z_yzz_xxyyzz[i] * pa_x[i];

        tr_z_xyzz_xxyzzz[i] = 2.0 * tr_z_yzz_xyzzz[i] * fe_0 + tr_z_yzz_xxyzzz[i] * pa_x[i];

        tr_z_xyzz_xxzzzz[i] = tr_z_xzz_xxzzzz[i] * pa_y[i];

        tr_z_xyzz_xyyyyy[i] = tr_z_yzz_yyyyy[i] * fe_0 + tr_z_yzz_xyyyyy[i] * pa_x[i];

        tr_z_xyzz_xyyyyz[i] = tr_z_yzz_yyyyz[i] * fe_0 + tr_z_yzz_xyyyyz[i] * pa_x[i];

        tr_z_xyzz_xyyyzz[i] = tr_z_yzz_yyyzz[i] * fe_0 + tr_z_yzz_xyyyzz[i] * pa_x[i];

        tr_z_xyzz_xyyzzz[i] = tr_z_yzz_yyzzz[i] * fe_0 + tr_z_yzz_xyyzzz[i] * pa_x[i];

        tr_z_xyzz_xyzzzz[i] = tr_z_yzz_yzzzz[i] * fe_0 + tr_z_yzz_xyzzzz[i] * pa_x[i];

        tr_z_xyzz_xzzzzz[i] = tr_z_xzz_xzzzzz[i] * pa_y[i];

        tr_z_xyzz_yyyyyy[i] = tr_z_yzz_yyyyyy[i] * pa_x[i];

        tr_z_xyzz_yyyyyz[i] = tr_z_yzz_yyyyyz[i] * pa_x[i];

        tr_z_xyzz_yyyyzz[i] = tr_z_yzz_yyyyzz[i] * pa_x[i];

        tr_z_xyzz_yyyzzz[i] = tr_z_yzz_yyyzzz[i] * pa_x[i];

        tr_z_xyzz_yyzzzz[i] = tr_z_yzz_yyzzzz[i] * pa_x[i];

        tr_z_xyzz_yzzzzz[i] = tr_z_yzz_yzzzzz[i] * pa_x[i];

        tr_z_xyzz_zzzzzz[i] = tr_z_yzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1092-1120 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_x, tr_z_xzzz_xxxxxx, tr_z_xzzz_xxxxxy, tr_z_xzzz_xxxxxz, tr_z_xzzz_xxxxyy, tr_z_xzzz_xxxxyz, tr_z_xzzz_xxxxzz, tr_z_xzzz_xxxyyy, tr_z_xzzz_xxxyyz, tr_z_xzzz_xxxyzz, tr_z_xzzz_xxxzzz, tr_z_xzzz_xxyyyy, tr_z_xzzz_xxyyyz, tr_z_xzzz_xxyyzz, tr_z_xzzz_xxyzzz, tr_z_xzzz_xxzzzz, tr_z_xzzz_xyyyyy, tr_z_xzzz_xyyyyz, tr_z_xzzz_xyyyzz, tr_z_xzzz_xyyzzz, tr_z_xzzz_xyzzzz, tr_z_xzzz_xzzzzz, tr_z_xzzz_yyyyyy, tr_z_xzzz_yyyyyz, tr_z_xzzz_yyyyzz, tr_z_xzzz_yyyzzz, tr_z_xzzz_yyzzzz, tr_z_xzzz_yzzzzz, tr_z_xzzz_zzzzzz, tr_z_zzz_xxxxx, tr_z_zzz_xxxxxx, tr_z_zzz_xxxxxy, tr_z_zzz_xxxxxz, tr_z_zzz_xxxxy, tr_z_zzz_xxxxyy, tr_z_zzz_xxxxyz, tr_z_zzz_xxxxz, tr_z_zzz_xxxxzz, tr_z_zzz_xxxyy, tr_z_zzz_xxxyyy, tr_z_zzz_xxxyyz, tr_z_zzz_xxxyz, tr_z_zzz_xxxyzz, tr_z_zzz_xxxzz, tr_z_zzz_xxxzzz, tr_z_zzz_xxyyy, tr_z_zzz_xxyyyy, tr_z_zzz_xxyyyz, tr_z_zzz_xxyyz, tr_z_zzz_xxyyzz, tr_z_zzz_xxyzz, tr_z_zzz_xxyzzz, tr_z_zzz_xxzzz, tr_z_zzz_xxzzzz, tr_z_zzz_xyyyy, tr_z_zzz_xyyyyy, tr_z_zzz_xyyyyz, tr_z_zzz_xyyyz, tr_z_zzz_xyyyzz, tr_z_zzz_xyyzz, tr_z_zzz_xyyzzz, tr_z_zzz_xyzzz, tr_z_zzz_xyzzzz, tr_z_zzz_xzzzz, tr_z_zzz_xzzzzz, tr_z_zzz_yyyyy, tr_z_zzz_yyyyyy, tr_z_zzz_yyyyyz, tr_z_zzz_yyyyz, tr_z_zzz_yyyyzz, tr_z_zzz_yyyzz, tr_z_zzz_yyyzzz, tr_z_zzz_yyzzz, tr_z_zzz_yyzzzz, tr_z_zzz_yzzzz, tr_z_zzz_yzzzzz, tr_z_zzz_zzzzz, tr_z_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzz_xxxxxx[i] = 6.0 * tr_z_zzz_xxxxx[i] * fe_0 + tr_z_zzz_xxxxxx[i] * pa_x[i];

        tr_z_xzzz_xxxxxy[i] = 5.0 * tr_z_zzz_xxxxy[i] * fe_0 + tr_z_zzz_xxxxxy[i] * pa_x[i];

        tr_z_xzzz_xxxxxz[i] = 5.0 * tr_z_zzz_xxxxz[i] * fe_0 + tr_z_zzz_xxxxxz[i] * pa_x[i];

        tr_z_xzzz_xxxxyy[i] = 4.0 * tr_z_zzz_xxxyy[i] * fe_0 + tr_z_zzz_xxxxyy[i] * pa_x[i];

        tr_z_xzzz_xxxxyz[i] = 4.0 * tr_z_zzz_xxxyz[i] * fe_0 + tr_z_zzz_xxxxyz[i] * pa_x[i];

        tr_z_xzzz_xxxxzz[i] = 4.0 * tr_z_zzz_xxxzz[i] * fe_0 + tr_z_zzz_xxxxzz[i] * pa_x[i];

        tr_z_xzzz_xxxyyy[i] = 3.0 * tr_z_zzz_xxyyy[i] * fe_0 + tr_z_zzz_xxxyyy[i] * pa_x[i];

        tr_z_xzzz_xxxyyz[i] = 3.0 * tr_z_zzz_xxyyz[i] * fe_0 + tr_z_zzz_xxxyyz[i] * pa_x[i];

        tr_z_xzzz_xxxyzz[i] = 3.0 * tr_z_zzz_xxyzz[i] * fe_0 + tr_z_zzz_xxxyzz[i] * pa_x[i];

        tr_z_xzzz_xxxzzz[i] = 3.0 * tr_z_zzz_xxzzz[i] * fe_0 + tr_z_zzz_xxxzzz[i] * pa_x[i];

        tr_z_xzzz_xxyyyy[i] = 2.0 * tr_z_zzz_xyyyy[i] * fe_0 + tr_z_zzz_xxyyyy[i] * pa_x[i];

        tr_z_xzzz_xxyyyz[i] = 2.0 * tr_z_zzz_xyyyz[i] * fe_0 + tr_z_zzz_xxyyyz[i] * pa_x[i];

        tr_z_xzzz_xxyyzz[i] = 2.0 * tr_z_zzz_xyyzz[i] * fe_0 + tr_z_zzz_xxyyzz[i] * pa_x[i];

        tr_z_xzzz_xxyzzz[i] = 2.0 * tr_z_zzz_xyzzz[i] * fe_0 + tr_z_zzz_xxyzzz[i] * pa_x[i];

        tr_z_xzzz_xxzzzz[i] = 2.0 * tr_z_zzz_xzzzz[i] * fe_0 + tr_z_zzz_xxzzzz[i] * pa_x[i];

        tr_z_xzzz_xyyyyy[i] = tr_z_zzz_yyyyy[i] * fe_0 + tr_z_zzz_xyyyyy[i] * pa_x[i];

        tr_z_xzzz_xyyyyz[i] = tr_z_zzz_yyyyz[i] * fe_0 + tr_z_zzz_xyyyyz[i] * pa_x[i];

        tr_z_xzzz_xyyyzz[i] = tr_z_zzz_yyyzz[i] * fe_0 + tr_z_zzz_xyyyzz[i] * pa_x[i];

        tr_z_xzzz_xyyzzz[i] = tr_z_zzz_yyzzz[i] * fe_0 + tr_z_zzz_xyyzzz[i] * pa_x[i];

        tr_z_xzzz_xyzzzz[i] = tr_z_zzz_yzzzz[i] * fe_0 + tr_z_zzz_xyzzzz[i] * pa_x[i];

        tr_z_xzzz_xzzzzz[i] = tr_z_zzz_zzzzz[i] * fe_0 + tr_z_zzz_xzzzzz[i] * pa_x[i];

        tr_z_xzzz_yyyyyy[i] = tr_z_zzz_yyyyyy[i] * pa_x[i];

        tr_z_xzzz_yyyyyz[i] = tr_z_zzz_yyyyyz[i] * pa_x[i];

        tr_z_xzzz_yyyyzz[i] = tr_z_zzz_yyyyzz[i] * pa_x[i];

        tr_z_xzzz_yyyzzz[i] = tr_z_zzz_yyyzzz[i] * pa_x[i];

        tr_z_xzzz_yyzzzz[i] = tr_z_zzz_yyzzzz[i] * pa_x[i];

        tr_z_xzzz_yzzzzz[i] = tr_z_zzz_yzzzzz[i] * pa_x[i];

        tr_z_xzzz_zzzzzz[i] = tr_z_zzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1120-1148 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, tr_z_yy_xxxxxx, tr_z_yy_xxxxxy, tr_z_yy_xxxxxz, tr_z_yy_xxxxyy, tr_z_yy_xxxxyz, tr_z_yy_xxxxzz, tr_z_yy_xxxyyy, tr_z_yy_xxxyyz, tr_z_yy_xxxyzz, tr_z_yy_xxxzzz, tr_z_yy_xxyyyy, tr_z_yy_xxyyyz, tr_z_yy_xxyyzz, tr_z_yy_xxyzzz, tr_z_yy_xxzzzz, tr_z_yy_xyyyyy, tr_z_yy_xyyyyz, tr_z_yy_xyyyzz, tr_z_yy_xyyzzz, tr_z_yy_xyzzzz, tr_z_yy_xzzzzz, tr_z_yy_yyyyyy, tr_z_yy_yyyyyz, tr_z_yy_yyyyzz, tr_z_yy_yyyzzz, tr_z_yy_yyzzzz, tr_z_yy_yzzzzz, tr_z_yy_zzzzzz, tr_z_yyy_xxxxx, tr_z_yyy_xxxxxx, tr_z_yyy_xxxxxy, tr_z_yyy_xxxxxz, tr_z_yyy_xxxxy, tr_z_yyy_xxxxyy, tr_z_yyy_xxxxyz, tr_z_yyy_xxxxz, tr_z_yyy_xxxxzz, tr_z_yyy_xxxyy, tr_z_yyy_xxxyyy, tr_z_yyy_xxxyyz, tr_z_yyy_xxxyz, tr_z_yyy_xxxyzz, tr_z_yyy_xxxzz, tr_z_yyy_xxxzzz, tr_z_yyy_xxyyy, tr_z_yyy_xxyyyy, tr_z_yyy_xxyyyz, tr_z_yyy_xxyyz, tr_z_yyy_xxyyzz, tr_z_yyy_xxyzz, tr_z_yyy_xxyzzz, tr_z_yyy_xxzzz, tr_z_yyy_xxzzzz, tr_z_yyy_xyyyy, tr_z_yyy_xyyyyy, tr_z_yyy_xyyyyz, tr_z_yyy_xyyyz, tr_z_yyy_xyyyzz, tr_z_yyy_xyyzz, tr_z_yyy_xyyzzz, tr_z_yyy_xyzzz, tr_z_yyy_xyzzzz, tr_z_yyy_xzzzz, tr_z_yyy_xzzzzz, tr_z_yyy_yyyyy, tr_z_yyy_yyyyyy, tr_z_yyy_yyyyyz, tr_z_yyy_yyyyz, tr_z_yyy_yyyyzz, tr_z_yyy_yyyzz, tr_z_yyy_yyyzzz, tr_z_yyy_yyzzz, tr_z_yyy_yyzzzz, tr_z_yyy_yzzzz, tr_z_yyy_yzzzzz, tr_z_yyy_zzzzz, tr_z_yyy_zzzzzz, tr_z_yyyy_xxxxxx, tr_z_yyyy_xxxxxy, tr_z_yyyy_xxxxxz, tr_z_yyyy_xxxxyy, tr_z_yyyy_xxxxyz, tr_z_yyyy_xxxxzz, tr_z_yyyy_xxxyyy, tr_z_yyyy_xxxyyz, tr_z_yyyy_xxxyzz, tr_z_yyyy_xxxzzz, tr_z_yyyy_xxyyyy, tr_z_yyyy_xxyyyz, tr_z_yyyy_xxyyzz, tr_z_yyyy_xxyzzz, tr_z_yyyy_xxzzzz, tr_z_yyyy_xyyyyy, tr_z_yyyy_xyyyyz, tr_z_yyyy_xyyyzz, tr_z_yyyy_xyyzzz, tr_z_yyyy_xyzzzz, tr_z_yyyy_xzzzzz, tr_z_yyyy_yyyyyy, tr_z_yyyy_yyyyyz, tr_z_yyyy_yyyyzz, tr_z_yyyy_yyyzzz, tr_z_yyyy_yyzzzz, tr_z_yyyy_yzzzzz, tr_z_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyy_xxxxxx[i] = 3.0 * tr_z_yy_xxxxxx[i] * fe_0 + tr_z_yyy_xxxxxx[i] * pa_y[i];

        tr_z_yyyy_xxxxxy[i] = 3.0 * tr_z_yy_xxxxxy[i] * fe_0 + tr_z_yyy_xxxxx[i] * fe_0 + tr_z_yyy_xxxxxy[i] * pa_y[i];

        tr_z_yyyy_xxxxxz[i] = 3.0 * tr_z_yy_xxxxxz[i] * fe_0 + tr_z_yyy_xxxxxz[i] * pa_y[i];

        tr_z_yyyy_xxxxyy[i] = 3.0 * tr_z_yy_xxxxyy[i] * fe_0 + 2.0 * tr_z_yyy_xxxxy[i] * fe_0 + tr_z_yyy_xxxxyy[i] * pa_y[i];

        tr_z_yyyy_xxxxyz[i] = 3.0 * tr_z_yy_xxxxyz[i] * fe_0 + tr_z_yyy_xxxxz[i] * fe_0 + tr_z_yyy_xxxxyz[i] * pa_y[i];

        tr_z_yyyy_xxxxzz[i] = 3.0 * tr_z_yy_xxxxzz[i] * fe_0 + tr_z_yyy_xxxxzz[i] * pa_y[i];

        tr_z_yyyy_xxxyyy[i] = 3.0 * tr_z_yy_xxxyyy[i] * fe_0 + 3.0 * tr_z_yyy_xxxyy[i] * fe_0 + tr_z_yyy_xxxyyy[i] * pa_y[i];

        tr_z_yyyy_xxxyyz[i] = 3.0 * tr_z_yy_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyy_xxxyz[i] * fe_0 + tr_z_yyy_xxxyyz[i] * pa_y[i];

        tr_z_yyyy_xxxyzz[i] = 3.0 * tr_z_yy_xxxyzz[i] * fe_0 + tr_z_yyy_xxxzz[i] * fe_0 + tr_z_yyy_xxxyzz[i] * pa_y[i];

        tr_z_yyyy_xxxzzz[i] = 3.0 * tr_z_yy_xxxzzz[i] * fe_0 + tr_z_yyy_xxxzzz[i] * pa_y[i];

        tr_z_yyyy_xxyyyy[i] = 3.0 * tr_z_yy_xxyyyy[i] * fe_0 + 4.0 * tr_z_yyy_xxyyy[i] * fe_0 + tr_z_yyy_xxyyyy[i] * pa_y[i];

        tr_z_yyyy_xxyyyz[i] = 3.0 * tr_z_yy_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyy_xxyyz[i] * fe_0 + tr_z_yyy_xxyyyz[i] * pa_y[i];

        tr_z_yyyy_xxyyzz[i] = 3.0 * tr_z_yy_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyy_xxyzz[i] * fe_0 + tr_z_yyy_xxyyzz[i] * pa_y[i];

        tr_z_yyyy_xxyzzz[i] = 3.0 * tr_z_yy_xxyzzz[i] * fe_0 + tr_z_yyy_xxzzz[i] * fe_0 + tr_z_yyy_xxyzzz[i] * pa_y[i];

        tr_z_yyyy_xxzzzz[i] = 3.0 * tr_z_yy_xxzzzz[i] * fe_0 + tr_z_yyy_xxzzzz[i] * pa_y[i];

        tr_z_yyyy_xyyyyy[i] = 3.0 * tr_z_yy_xyyyyy[i] * fe_0 + 5.0 * tr_z_yyy_xyyyy[i] * fe_0 + tr_z_yyy_xyyyyy[i] * pa_y[i];

        tr_z_yyyy_xyyyyz[i] = 3.0 * tr_z_yy_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyy_xyyyz[i] * fe_0 + tr_z_yyy_xyyyyz[i] * pa_y[i];

        tr_z_yyyy_xyyyzz[i] = 3.0 * tr_z_yy_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyy_xyyzz[i] * fe_0 + tr_z_yyy_xyyyzz[i] * pa_y[i];

        tr_z_yyyy_xyyzzz[i] = 3.0 * tr_z_yy_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyy_xyzzz[i] * fe_0 + tr_z_yyy_xyyzzz[i] * pa_y[i];

        tr_z_yyyy_xyzzzz[i] = 3.0 * tr_z_yy_xyzzzz[i] * fe_0 + tr_z_yyy_xzzzz[i] * fe_0 + tr_z_yyy_xyzzzz[i] * pa_y[i];

        tr_z_yyyy_xzzzzz[i] = 3.0 * tr_z_yy_xzzzzz[i] * fe_0 + tr_z_yyy_xzzzzz[i] * pa_y[i];

        tr_z_yyyy_yyyyyy[i] = 3.0 * tr_z_yy_yyyyyy[i] * fe_0 + 6.0 * tr_z_yyy_yyyyy[i] * fe_0 + tr_z_yyy_yyyyyy[i] * pa_y[i];

        tr_z_yyyy_yyyyyz[i] = 3.0 * tr_z_yy_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyy_yyyyz[i] * fe_0 + tr_z_yyy_yyyyyz[i] * pa_y[i];

        tr_z_yyyy_yyyyzz[i] = 3.0 * tr_z_yy_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyy_yyyzz[i] * fe_0 + tr_z_yyy_yyyyzz[i] * pa_y[i];

        tr_z_yyyy_yyyzzz[i] = 3.0 * tr_z_yy_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyy_yyzzz[i] * fe_0 + tr_z_yyy_yyyzzz[i] * pa_y[i];

        tr_z_yyyy_yyzzzz[i] = 3.0 * tr_z_yy_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyy_yzzzz[i] * fe_0 + tr_z_yyy_yyzzzz[i] * pa_y[i];

        tr_z_yyyy_yzzzzz[i] = 3.0 * tr_z_yy_yzzzzz[i] * fe_0 + tr_z_yyy_zzzzz[i] * fe_0 + tr_z_yyy_yzzzzz[i] * pa_y[i];

        tr_z_yyyy_zzzzzz[i] = 3.0 * tr_z_yy_zzzzzz[i] * fe_0 + tr_z_yyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 1148-1176 components of targeted buffer : GI

    auto tr_z_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 1148);

    auto tr_z_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 1149);

    auto tr_z_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 1150);

    auto tr_z_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 1151);

    auto tr_z_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 1152);

    auto tr_z_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 1153);

    auto tr_z_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 1154);

    auto tr_z_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 1155);

    auto tr_z_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 1156);

    auto tr_z_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 1157);

    auto tr_z_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 1158);

    auto tr_z_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 1159);

    auto tr_z_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 1160);

    auto tr_z_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 1161);

    auto tr_z_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 1162);

    auto tr_z_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 1163);

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

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yyy_xxxxxy, tr_z_yyy_xxxxyy, tr_z_yyy_xxxyyy, tr_z_yyy_xxyyyy, tr_z_yyy_xyyyyy, tr_z_yyy_yyyyyy, tr_z_yyyz_xxxxxx, tr_z_yyyz_xxxxxy, tr_z_yyyz_xxxxxz, tr_z_yyyz_xxxxyy, tr_z_yyyz_xxxxyz, tr_z_yyyz_xxxxzz, tr_z_yyyz_xxxyyy, tr_z_yyyz_xxxyyz, tr_z_yyyz_xxxyzz, tr_z_yyyz_xxxzzz, tr_z_yyyz_xxyyyy, tr_z_yyyz_xxyyyz, tr_z_yyyz_xxyyzz, tr_z_yyyz_xxyzzz, tr_z_yyyz_xxzzzz, tr_z_yyyz_xyyyyy, tr_z_yyyz_xyyyyz, tr_z_yyyz_xyyyzz, tr_z_yyyz_xyyzzz, tr_z_yyyz_xyzzzz, tr_z_yyyz_xzzzzz, tr_z_yyyz_yyyyyy, tr_z_yyyz_yyyyyz, tr_z_yyyz_yyyyzz, tr_z_yyyz_yyyzzz, tr_z_yyyz_yyzzzz, tr_z_yyyz_yzzzzz, tr_z_yyyz_zzzzzz, tr_z_yyz_xxxxxx, tr_z_yyz_xxxxxz, tr_z_yyz_xxxxyz, tr_z_yyz_xxxxz, tr_z_yyz_xxxxzz, tr_z_yyz_xxxyyz, tr_z_yyz_xxxyz, tr_z_yyz_xxxyzz, tr_z_yyz_xxxzz, tr_z_yyz_xxxzzz, tr_z_yyz_xxyyyz, tr_z_yyz_xxyyz, tr_z_yyz_xxyyzz, tr_z_yyz_xxyzz, tr_z_yyz_xxyzzz, tr_z_yyz_xxzzz, tr_z_yyz_xxzzzz, tr_z_yyz_xyyyyz, tr_z_yyz_xyyyz, tr_z_yyz_xyyyzz, tr_z_yyz_xyyzz, tr_z_yyz_xyyzzz, tr_z_yyz_xyzzz, tr_z_yyz_xyzzzz, tr_z_yyz_xzzzz, tr_z_yyz_xzzzzz, tr_z_yyz_yyyyyz, tr_z_yyz_yyyyz, tr_z_yyz_yyyyzz, tr_z_yyz_yyyzz, tr_z_yyz_yyyzzz, tr_z_yyz_yyzzz, tr_z_yyz_yyzzzz, tr_z_yyz_yzzzz, tr_z_yyz_yzzzzz, tr_z_yyz_zzzzz, tr_z_yyz_zzzzzz, tr_z_yz_xxxxxx, tr_z_yz_xxxxxz, tr_z_yz_xxxxyz, tr_z_yz_xxxxzz, tr_z_yz_xxxyyz, tr_z_yz_xxxyzz, tr_z_yz_xxxzzz, tr_z_yz_xxyyyz, tr_z_yz_xxyyzz, tr_z_yz_xxyzzz, tr_z_yz_xxzzzz, tr_z_yz_xyyyyz, tr_z_yz_xyyyzz, tr_z_yz_xyyzzz, tr_z_yz_xyzzzz, tr_z_yz_xzzzzz, tr_z_yz_yyyyyz, tr_z_yz_yyyyzz, tr_z_yz_yyyzzz, tr_z_yz_yyzzzz, tr_z_yz_yzzzzz, tr_z_yz_zzzzzz, ts_yyy_xxxxxy, ts_yyy_xxxxyy, ts_yyy_xxxyyy, ts_yyy_xxyyyy, ts_yyy_xyyyyy, ts_yyy_yyyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyz_xxxxxx[i] = 2.0 * tr_z_yz_xxxxxx[i] * fe_0 + tr_z_yyz_xxxxxx[i] * pa_y[i];

        tr_z_yyyz_xxxxxy[i] = ts_yyy_xxxxxy[i] * fe_0 + tr_z_yyy_xxxxxy[i] * pa_z[i];

        tr_z_yyyz_xxxxxz[i] = 2.0 * tr_z_yz_xxxxxz[i] * fe_0 + tr_z_yyz_xxxxxz[i] * pa_y[i];

        tr_z_yyyz_xxxxyy[i] = ts_yyy_xxxxyy[i] * fe_0 + tr_z_yyy_xxxxyy[i] * pa_z[i];

        tr_z_yyyz_xxxxyz[i] = 2.0 * tr_z_yz_xxxxyz[i] * fe_0 + tr_z_yyz_xxxxz[i] * fe_0 + tr_z_yyz_xxxxyz[i] * pa_y[i];

        tr_z_yyyz_xxxxzz[i] = 2.0 * tr_z_yz_xxxxzz[i] * fe_0 + tr_z_yyz_xxxxzz[i] * pa_y[i];

        tr_z_yyyz_xxxyyy[i] = ts_yyy_xxxyyy[i] * fe_0 + tr_z_yyy_xxxyyy[i] * pa_z[i];

        tr_z_yyyz_xxxyyz[i] = 2.0 * tr_z_yz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyz_xxxyz[i] * fe_0 + tr_z_yyz_xxxyyz[i] * pa_y[i];

        tr_z_yyyz_xxxyzz[i] = 2.0 * tr_z_yz_xxxyzz[i] * fe_0 + tr_z_yyz_xxxzz[i] * fe_0 + tr_z_yyz_xxxyzz[i] * pa_y[i];

        tr_z_yyyz_xxxzzz[i] = 2.0 * tr_z_yz_xxxzzz[i] * fe_0 + tr_z_yyz_xxxzzz[i] * pa_y[i];

        tr_z_yyyz_xxyyyy[i] = ts_yyy_xxyyyy[i] * fe_0 + tr_z_yyy_xxyyyy[i] * pa_z[i];

        tr_z_yyyz_xxyyyz[i] = 2.0 * tr_z_yz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyz_xxyyz[i] * fe_0 + tr_z_yyz_xxyyyz[i] * pa_y[i];

        tr_z_yyyz_xxyyzz[i] = 2.0 * tr_z_yz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyz_xxyzz[i] * fe_0 + tr_z_yyz_xxyyzz[i] * pa_y[i];

        tr_z_yyyz_xxyzzz[i] = 2.0 * tr_z_yz_xxyzzz[i] * fe_0 + tr_z_yyz_xxzzz[i] * fe_0 + tr_z_yyz_xxyzzz[i] * pa_y[i];

        tr_z_yyyz_xxzzzz[i] = 2.0 * tr_z_yz_xxzzzz[i] * fe_0 + tr_z_yyz_xxzzzz[i] * pa_y[i];

        tr_z_yyyz_xyyyyy[i] = ts_yyy_xyyyyy[i] * fe_0 + tr_z_yyy_xyyyyy[i] * pa_z[i];

        tr_z_yyyz_xyyyyz[i] = 2.0 * tr_z_yz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyz_xyyyz[i] * fe_0 + tr_z_yyz_xyyyyz[i] * pa_y[i];

        tr_z_yyyz_xyyyzz[i] = 2.0 * tr_z_yz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyz_xyyzz[i] * fe_0 + tr_z_yyz_xyyyzz[i] * pa_y[i];

        tr_z_yyyz_xyyzzz[i] = 2.0 * tr_z_yz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyz_xyzzz[i] * fe_0 + tr_z_yyz_xyyzzz[i] * pa_y[i];

        tr_z_yyyz_xyzzzz[i] = 2.0 * tr_z_yz_xyzzzz[i] * fe_0 + tr_z_yyz_xzzzz[i] * fe_0 + tr_z_yyz_xyzzzz[i] * pa_y[i];

        tr_z_yyyz_xzzzzz[i] = 2.0 * tr_z_yz_xzzzzz[i] * fe_0 + tr_z_yyz_xzzzzz[i] * pa_y[i];

        tr_z_yyyz_yyyyyy[i] = ts_yyy_yyyyyy[i] * fe_0 + tr_z_yyy_yyyyyy[i] * pa_z[i];

        tr_z_yyyz_yyyyyz[i] = 2.0 * tr_z_yz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyz_yyyyz[i] * fe_0 + tr_z_yyz_yyyyyz[i] * pa_y[i];

        tr_z_yyyz_yyyyzz[i] = 2.0 * tr_z_yz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyz_yyyzz[i] * fe_0 + tr_z_yyz_yyyyzz[i] * pa_y[i];

        tr_z_yyyz_yyyzzz[i] = 2.0 * tr_z_yz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyz_yyzzz[i] * fe_0 + tr_z_yyz_yyyzzz[i] * pa_y[i];

        tr_z_yyyz_yyzzzz[i] = 2.0 * tr_z_yz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyz_yzzzz[i] * fe_0 + tr_z_yyz_yyzzzz[i] * pa_y[i];

        tr_z_yyyz_yzzzzz[i] = 2.0 * tr_z_yz_yzzzzz[i] * fe_0 + tr_z_yyz_zzzzz[i] * fe_0 + tr_z_yyz_yzzzzz[i] * pa_y[i];

        tr_z_yyyz_zzzzzz[i] = 2.0 * tr_z_yz_zzzzzz[i] * fe_0 + tr_z_yyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1176-1204 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, tr_z_yyzz_xxxxxx, tr_z_yyzz_xxxxxy, tr_z_yyzz_xxxxxz, tr_z_yyzz_xxxxyy, tr_z_yyzz_xxxxyz, tr_z_yyzz_xxxxzz, tr_z_yyzz_xxxyyy, tr_z_yyzz_xxxyyz, tr_z_yyzz_xxxyzz, tr_z_yyzz_xxxzzz, tr_z_yyzz_xxyyyy, tr_z_yyzz_xxyyyz, tr_z_yyzz_xxyyzz, tr_z_yyzz_xxyzzz, tr_z_yyzz_xxzzzz, tr_z_yyzz_xyyyyy, tr_z_yyzz_xyyyyz, tr_z_yyzz_xyyyzz, tr_z_yyzz_xyyzzz, tr_z_yyzz_xyzzzz, tr_z_yyzz_xzzzzz, tr_z_yyzz_yyyyyy, tr_z_yyzz_yyyyyz, tr_z_yyzz_yyyyzz, tr_z_yyzz_yyyzzz, tr_z_yyzz_yyzzzz, tr_z_yyzz_yzzzzz, tr_z_yyzz_zzzzzz, tr_z_yzz_xxxxx, tr_z_yzz_xxxxxx, tr_z_yzz_xxxxxy, tr_z_yzz_xxxxxz, tr_z_yzz_xxxxy, tr_z_yzz_xxxxyy, tr_z_yzz_xxxxyz, tr_z_yzz_xxxxz, tr_z_yzz_xxxxzz, tr_z_yzz_xxxyy, tr_z_yzz_xxxyyy, tr_z_yzz_xxxyyz, tr_z_yzz_xxxyz, tr_z_yzz_xxxyzz, tr_z_yzz_xxxzz, tr_z_yzz_xxxzzz, tr_z_yzz_xxyyy, tr_z_yzz_xxyyyy, tr_z_yzz_xxyyyz, tr_z_yzz_xxyyz, tr_z_yzz_xxyyzz, tr_z_yzz_xxyzz, tr_z_yzz_xxyzzz, tr_z_yzz_xxzzz, tr_z_yzz_xxzzzz, tr_z_yzz_xyyyy, tr_z_yzz_xyyyyy, tr_z_yzz_xyyyyz, tr_z_yzz_xyyyz, tr_z_yzz_xyyyzz, tr_z_yzz_xyyzz, tr_z_yzz_xyyzzz, tr_z_yzz_xyzzz, tr_z_yzz_xyzzzz, tr_z_yzz_xzzzz, tr_z_yzz_xzzzzz, tr_z_yzz_yyyyy, tr_z_yzz_yyyyyy, tr_z_yzz_yyyyyz, tr_z_yzz_yyyyz, tr_z_yzz_yyyyzz, tr_z_yzz_yyyzz, tr_z_yzz_yyyzzz, tr_z_yzz_yyzzz, tr_z_yzz_yyzzzz, tr_z_yzz_yzzzz, tr_z_yzz_yzzzzz, tr_z_yzz_zzzzz, tr_z_yzz_zzzzzz, tr_z_zz_xxxxxx, tr_z_zz_xxxxxy, tr_z_zz_xxxxxz, tr_z_zz_xxxxyy, tr_z_zz_xxxxyz, tr_z_zz_xxxxzz, tr_z_zz_xxxyyy, tr_z_zz_xxxyyz, tr_z_zz_xxxyzz, tr_z_zz_xxxzzz, tr_z_zz_xxyyyy, tr_z_zz_xxyyyz, tr_z_zz_xxyyzz, tr_z_zz_xxyzzz, tr_z_zz_xxzzzz, tr_z_zz_xyyyyy, tr_z_zz_xyyyyz, tr_z_zz_xyyyzz, tr_z_zz_xyyzzz, tr_z_zz_xyzzzz, tr_z_zz_xzzzzz, tr_z_zz_yyyyyy, tr_z_zz_yyyyyz, tr_z_zz_yyyyzz, tr_z_zz_yyyzzz, tr_z_zz_yyzzzz, tr_z_zz_yzzzzz, tr_z_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzz_xxxxxx[i] = tr_z_zz_xxxxxx[i] * fe_0 + tr_z_yzz_xxxxxx[i] * pa_y[i];

        tr_z_yyzz_xxxxxy[i] = tr_z_zz_xxxxxy[i] * fe_0 + tr_z_yzz_xxxxx[i] * fe_0 + tr_z_yzz_xxxxxy[i] * pa_y[i];

        tr_z_yyzz_xxxxxz[i] = tr_z_zz_xxxxxz[i] * fe_0 + tr_z_yzz_xxxxxz[i] * pa_y[i];

        tr_z_yyzz_xxxxyy[i] = tr_z_zz_xxxxyy[i] * fe_0 + 2.0 * tr_z_yzz_xxxxy[i] * fe_0 + tr_z_yzz_xxxxyy[i] * pa_y[i];

        tr_z_yyzz_xxxxyz[i] = tr_z_zz_xxxxyz[i] * fe_0 + tr_z_yzz_xxxxz[i] * fe_0 + tr_z_yzz_xxxxyz[i] * pa_y[i];

        tr_z_yyzz_xxxxzz[i] = tr_z_zz_xxxxzz[i] * fe_0 + tr_z_yzz_xxxxzz[i] * pa_y[i];

        tr_z_yyzz_xxxyyy[i] = tr_z_zz_xxxyyy[i] * fe_0 + 3.0 * tr_z_yzz_xxxyy[i] * fe_0 + tr_z_yzz_xxxyyy[i] * pa_y[i];

        tr_z_yyzz_xxxyyz[i] = tr_z_zz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yzz_xxxyz[i] * fe_0 + tr_z_yzz_xxxyyz[i] * pa_y[i];

        tr_z_yyzz_xxxyzz[i] = tr_z_zz_xxxyzz[i] * fe_0 + tr_z_yzz_xxxzz[i] * fe_0 + tr_z_yzz_xxxyzz[i] * pa_y[i];

        tr_z_yyzz_xxxzzz[i] = tr_z_zz_xxxzzz[i] * fe_0 + tr_z_yzz_xxxzzz[i] * pa_y[i];

        tr_z_yyzz_xxyyyy[i] = tr_z_zz_xxyyyy[i] * fe_0 + 4.0 * tr_z_yzz_xxyyy[i] * fe_0 + tr_z_yzz_xxyyyy[i] * pa_y[i];

        tr_z_yyzz_xxyyyz[i] = tr_z_zz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yzz_xxyyz[i] * fe_0 + tr_z_yzz_xxyyyz[i] * pa_y[i];

        tr_z_yyzz_xxyyzz[i] = tr_z_zz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yzz_xxyzz[i] * fe_0 + tr_z_yzz_xxyyzz[i] * pa_y[i];

        tr_z_yyzz_xxyzzz[i] = tr_z_zz_xxyzzz[i] * fe_0 + tr_z_yzz_xxzzz[i] * fe_0 + tr_z_yzz_xxyzzz[i] * pa_y[i];

        tr_z_yyzz_xxzzzz[i] = tr_z_zz_xxzzzz[i] * fe_0 + tr_z_yzz_xxzzzz[i] * pa_y[i];

        tr_z_yyzz_xyyyyy[i] = tr_z_zz_xyyyyy[i] * fe_0 + 5.0 * tr_z_yzz_xyyyy[i] * fe_0 + tr_z_yzz_xyyyyy[i] * pa_y[i];

        tr_z_yyzz_xyyyyz[i] = tr_z_zz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yzz_xyyyz[i] * fe_0 + tr_z_yzz_xyyyyz[i] * pa_y[i];

        tr_z_yyzz_xyyyzz[i] = tr_z_zz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yzz_xyyzz[i] * fe_0 + tr_z_yzz_xyyyzz[i] * pa_y[i];

        tr_z_yyzz_xyyzzz[i] = tr_z_zz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yzz_xyzzz[i] * fe_0 + tr_z_yzz_xyyzzz[i] * pa_y[i];

        tr_z_yyzz_xyzzzz[i] = tr_z_zz_xyzzzz[i] * fe_0 + tr_z_yzz_xzzzz[i] * fe_0 + tr_z_yzz_xyzzzz[i] * pa_y[i];

        tr_z_yyzz_xzzzzz[i] = tr_z_zz_xzzzzz[i] * fe_0 + tr_z_yzz_xzzzzz[i] * pa_y[i];

        tr_z_yyzz_yyyyyy[i] = tr_z_zz_yyyyyy[i] * fe_0 + 6.0 * tr_z_yzz_yyyyy[i] * fe_0 + tr_z_yzz_yyyyyy[i] * pa_y[i];

        tr_z_yyzz_yyyyyz[i] = tr_z_zz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yzz_yyyyz[i] * fe_0 + tr_z_yzz_yyyyyz[i] * pa_y[i];

        tr_z_yyzz_yyyyzz[i] = tr_z_zz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yzz_yyyzz[i] * fe_0 + tr_z_yzz_yyyyzz[i] * pa_y[i];

        tr_z_yyzz_yyyzzz[i] = tr_z_zz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yzz_yyzzz[i] * fe_0 + tr_z_yzz_yyyzzz[i] * pa_y[i];

        tr_z_yyzz_yyzzzz[i] = tr_z_zz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yzz_yzzzz[i] * fe_0 + tr_z_yzz_yyzzzz[i] * pa_y[i];

        tr_z_yyzz_yzzzzz[i] = tr_z_zz_yzzzzz[i] * fe_0 + tr_z_yzz_zzzzz[i] * fe_0 + tr_z_yzz_yzzzzz[i] * pa_y[i];

        tr_z_yyzz_zzzzzz[i] = tr_z_zz_zzzzzz[i] * fe_0 + tr_z_yzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1204-1232 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_y, tr_z_yzzz_xxxxxx, tr_z_yzzz_xxxxxy, tr_z_yzzz_xxxxxz, tr_z_yzzz_xxxxyy, tr_z_yzzz_xxxxyz, tr_z_yzzz_xxxxzz, tr_z_yzzz_xxxyyy, tr_z_yzzz_xxxyyz, tr_z_yzzz_xxxyzz, tr_z_yzzz_xxxzzz, tr_z_yzzz_xxyyyy, tr_z_yzzz_xxyyyz, tr_z_yzzz_xxyyzz, tr_z_yzzz_xxyzzz, tr_z_yzzz_xxzzzz, tr_z_yzzz_xyyyyy, tr_z_yzzz_xyyyyz, tr_z_yzzz_xyyyzz, tr_z_yzzz_xyyzzz, tr_z_yzzz_xyzzzz, tr_z_yzzz_xzzzzz, tr_z_yzzz_yyyyyy, tr_z_yzzz_yyyyyz, tr_z_yzzz_yyyyzz, tr_z_yzzz_yyyzzz, tr_z_yzzz_yyzzzz, tr_z_yzzz_yzzzzz, tr_z_yzzz_zzzzzz, tr_z_zzz_xxxxx, tr_z_zzz_xxxxxx, tr_z_zzz_xxxxxy, tr_z_zzz_xxxxxz, tr_z_zzz_xxxxy, tr_z_zzz_xxxxyy, tr_z_zzz_xxxxyz, tr_z_zzz_xxxxz, tr_z_zzz_xxxxzz, tr_z_zzz_xxxyy, tr_z_zzz_xxxyyy, tr_z_zzz_xxxyyz, tr_z_zzz_xxxyz, tr_z_zzz_xxxyzz, tr_z_zzz_xxxzz, tr_z_zzz_xxxzzz, tr_z_zzz_xxyyy, tr_z_zzz_xxyyyy, tr_z_zzz_xxyyyz, tr_z_zzz_xxyyz, tr_z_zzz_xxyyzz, tr_z_zzz_xxyzz, tr_z_zzz_xxyzzz, tr_z_zzz_xxzzz, tr_z_zzz_xxzzzz, tr_z_zzz_xyyyy, tr_z_zzz_xyyyyy, tr_z_zzz_xyyyyz, tr_z_zzz_xyyyz, tr_z_zzz_xyyyzz, tr_z_zzz_xyyzz, tr_z_zzz_xyyzzz, tr_z_zzz_xyzzz, tr_z_zzz_xyzzzz, tr_z_zzz_xzzzz, tr_z_zzz_xzzzzz, tr_z_zzz_yyyyy, tr_z_zzz_yyyyyy, tr_z_zzz_yyyyyz, tr_z_zzz_yyyyz, tr_z_zzz_yyyyzz, tr_z_zzz_yyyzz, tr_z_zzz_yyyzzz, tr_z_zzz_yyzzz, tr_z_zzz_yyzzzz, tr_z_zzz_yzzzz, tr_z_zzz_yzzzzz, tr_z_zzz_zzzzz, tr_z_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzz_xxxxxx[i] = tr_z_zzz_xxxxxx[i] * pa_y[i];

        tr_z_yzzz_xxxxxy[i] = tr_z_zzz_xxxxx[i] * fe_0 + tr_z_zzz_xxxxxy[i] * pa_y[i];

        tr_z_yzzz_xxxxxz[i] = tr_z_zzz_xxxxxz[i] * pa_y[i];

        tr_z_yzzz_xxxxyy[i] = 2.0 * tr_z_zzz_xxxxy[i] * fe_0 + tr_z_zzz_xxxxyy[i] * pa_y[i];

        tr_z_yzzz_xxxxyz[i] = tr_z_zzz_xxxxz[i] * fe_0 + tr_z_zzz_xxxxyz[i] * pa_y[i];

        tr_z_yzzz_xxxxzz[i] = tr_z_zzz_xxxxzz[i] * pa_y[i];

        tr_z_yzzz_xxxyyy[i] = 3.0 * tr_z_zzz_xxxyy[i] * fe_0 + tr_z_zzz_xxxyyy[i] * pa_y[i];

        tr_z_yzzz_xxxyyz[i] = 2.0 * tr_z_zzz_xxxyz[i] * fe_0 + tr_z_zzz_xxxyyz[i] * pa_y[i];

        tr_z_yzzz_xxxyzz[i] = tr_z_zzz_xxxzz[i] * fe_0 + tr_z_zzz_xxxyzz[i] * pa_y[i];

        tr_z_yzzz_xxxzzz[i] = tr_z_zzz_xxxzzz[i] * pa_y[i];

        tr_z_yzzz_xxyyyy[i] = 4.0 * tr_z_zzz_xxyyy[i] * fe_0 + tr_z_zzz_xxyyyy[i] * pa_y[i];

        tr_z_yzzz_xxyyyz[i] = 3.0 * tr_z_zzz_xxyyz[i] * fe_0 + tr_z_zzz_xxyyyz[i] * pa_y[i];

        tr_z_yzzz_xxyyzz[i] = 2.0 * tr_z_zzz_xxyzz[i] * fe_0 + tr_z_zzz_xxyyzz[i] * pa_y[i];

        tr_z_yzzz_xxyzzz[i] = tr_z_zzz_xxzzz[i] * fe_0 + tr_z_zzz_xxyzzz[i] * pa_y[i];

        tr_z_yzzz_xxzzzz[i] = tr_z_zzz_xxzzzz[i] * pa_y[i];

        tr_z_yzzz_xyyyyy[i] = 5.0 * tr_z_zzz_xyyyy[i] * fe_0 + tr_z_zzz_xyyyyy[i] * pa_y[i];

        tr_z_yzzz_xyyyyz[i] = 4.0 * tr_z_zzz_xyyyz[i] * fe_0 + tr_z_zzz_xyyyyz[i] * pa_y[i];

        tr_z_yzzz_xyyyzz[i] = 3.0 * tr_z_zzz_xyyzz[i] * fe_0 + tr_z_zzz_xyyyzz[i] * pa_y[i];

        tr_z_yzzz_xyyzzz[i] = 2.0 * tr_z_zzz_xyzzz[i] * fe_0 + tr_z_zzz_xyyzzz[i] * pa_y[i];

        tr_z_yzzz_xyzzzz[i] = tr_z_zzz_xzzzz[i] * fe_0 + tr_z_zzz_xyzzzz[i] * pa_y[i];

        tr_z_yzzz_xzzzzz[i] = tr_z_zzz_xzzzzz[i] * pa_y[i];

        tr_z_yzzz_yyyyyy[i] = 6.0 * tr_z_zzz_yyyyy[i] * fe_0 + tr_z_zzz_yyyyyy[i] * pa_y[i];

        tr_z_yzzz_yyyyyz[i] = 5.0 * tr_z_zzz_yyyyz[i] * fe_0 + tr_z_zzz_yyyyyz[i] * pa_y[i];

        tr_z_yzzz_yyyyzz[i] = 4.0 * tr_z_zzz_yyyzz[i] * fe_0 + tr_z_zzz_yyyyzz[i] * pa_y[i];

        tr_z_yzzz_yyyzzz[i] = 3.0 * tr_z_zzz_yyzzz[i] * fe_0 + tr_z_zzz_yyyzzz[i] * pa_y[i];

        tr_z_yzzz_yyzzzz[i] = 2.0 * tr_z_zzz_yzzzz[i] * fe_0 + tr_z_zzz_yyzzzz[i] * pa_y[i];

        tr_z_yzzz_yzzzzz[i] = tr_z_zzz_zzzzz[i] * fe_0 + tr_z_zzz_yzzzzz[i] * pa_y[i];

        tr_z_yzzz_zzzzzz[i] = tr_z_zzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1232-1260 components of targeted buffer : GI

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

    #pragma omp simd aligned(pa_z, tr_z_zz_xxxxxx, tr_z_zz_xxxxxy, tr_z_zz_xxxxxz, tr_z_zz_xxxxyy, tr_z_zz_xxxxyz, tr_z_zz_xxxxzz, tr_z_zz_xxxyyy, tr_z_zz_xxxyyz, tr_z_zz_xxxyzz, tr_z_zz_xxxzzz, tr_z_zz_xxyyyy, tr_z_zz_xxyyyz, tr_z_zz_xxyyzz, tr_z_zz_xxyzzz, tr_z_zz_xxzzzz, tr_z_zz_xyyyyy, tr_z_zz_xyyyyz, tr_z_zz_xyyyzz, tr_z_zz_xyyzzz, tr_z_zz_xyzzzz, tr_z_zz_xzzzzz, tr_z_zz_yyyyyy, tr_z_zz_yyyyyz, tr_z_zz_yyyyzz, tr_z_zz_yyyzzz, tr_z_zz_yyzzzz, tr_z_zz_yzzzzz, tr_z_zz_zzzzzz, tr_z_zzz_xxxxx, tr_z_zzz_xxxxxx, tr_z_zzz_xxxxxy, tr_z_zzz_xxxxxz, tr_z_zzz_xxxxy, tr_z_zzz_xxxxyy, tr_z_zzz_xxxxyz, tr_z_zzz_xxxxz, tr_z_zzz_xxxxzz, tr_z_zzz_xxxyy, tr_z_zzz_xxxyyy, tr_z_zzz_xxxyyz, tr_z_zzz_xxxyz, tr_z_zzz_xxxyzz, tr_z_zzz_xxxzz, tr_z_zzz_xxxzzz, tr_z_zzz_xxyyy, tr_z_zzz_xxyyyy, tr_z_zzz_xxyyyz, tr_z_zzz_xxyyz, tr_z_zzz_xxyyzz, tr_z_zzz_xxyzz, tr_z_zzz_xxyzzz, tr_z_zzz_xxzzz, tr_z_zzz_xxzzzz, tr_z_zzz_xyyyy, tr_z_zzz_xyyyyy, tr_z_zzz_xyyyyz, tr_z_zzz_xyyyz, tr_z_zzz_xyyyzz, tr_z_zzz_xyyzz, tr_z_zzz_xyyzzz, tr_z_zzz_xyzzz, tr_z_zzz_xyzzzz, tr_z_zzz_xzzzz, tr_z_zzz_xzzzzz, tr_z_zzz_yyyyy, tr_z_zzz_yyyyyy, tr_z_zzz_yyyyyz, tr_z_zzz_yyyyz, tr_z_zzz_yyyyzz, tr_z_zzz_yyyzz, tr_z_zzz_yyyzzz, tr_z_zzz_yyzzz, tr_z_zzz_yyzzzz, tr_z_zzz_yzzzz, tr_z_zzz_yzzzzz, tr_z_zzz_zzzzz, tr_z_zzz_zzzzzz, tr_z_zzzz_xxxxxx, tr_z_zzzz_xxxxxy, tr_z_zzzz_xxxxxz, tr_z_zzzz_xxxxyy, tr_z_zzzz_xxxxyz, tr_z_zzzz_xxxxzz, tr_z_zzzz_xxxyyy, tr_z_zzzz_xxxyyz, tr_z_zzzz_xxxyzz, tr_z_zzzz_xxxzzz, tr_z_zzzz_xxyyyy, tr_z_zzzz_xxyyyz, tr_z_zzzz_xxyyzz, tr_z_zzzz_xxyzzz, tr_z_zzzz_xxzzzz, tr_z_zzzz_xyyyyy, tr_z_zzzz_xyyyyz, tr_z_zzzz_xyyyzz, tr_z_zzzz_xyyzzz, tr_z_zzzz_xyzzzz, tr_z_zzzz_xzzzzz, tr_z_zzzz_yyyyyy, tr_z_zzzz_yyyyyz, tr_z_zzzz_yyyyzz, tr_z_zzzz_yyyzzz, tr_z_zzzz_yyzzzz, tr_z_zzzz_yzzzzz, tr_z_zzzz_zzzzzz, ts_zzz_xxxxxx, ts_zzz_xxxxxy, ts_zzz_xxxxxz, ts_zzz_xxxxyy, ts_zzz_xxxxyz, ts_zzz_xxxxzz, ts_zzz_xxxyyy, ts_zzz_xxxyyz, ts_zzz_xxxyzz, ts_zzz_xxxzzz, ts_zzz_xxyyyy, ts_zzz_xxyyyz, ts_zzz_xxyyzz, ts_zzz_xxyzzz, ts_zzz_xxzzzz, ts_zzz_xyyyyy, ts_zzz_xyyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzz_xxxxxx[i] = 3.0 * tr_z_zz_xxxxxx[i] * fe_0 + ts_zzz_xxxxxx[i] * fe_0 + tr_z_zzz_xxxxxx[i] * pa_z[i];

        tr_z_zzzz_xxxxxy[i] = 3.0 * tr_z_zz_xxxxxy[i] * fe_0 + ts_zzz_xxxxxy[i] * fe_0 + tr_z_zzz_xxxxxy[i] * pa_z[i];

        tr_z_zzzz_xxxxxz[i] = 3.0 * tr_z_zz_xxxxxz[i] * fe_0 + tr_z_zzz_xxxxx[i] * fe_0 + ts_zzz_xxxxxz[i] * fe_0 + tr_z_zzz_xxxxxz[i] * pa_z[i];

        tr_z_zzzz_xxxxyy[i] = 3.0 * tr_z_zz_xxxxyy[i] * fe_0 + ts_zzz_xxxxyy[i] * fe_0 + tr_z_zzz_xxxxyy[i] * pa_z[i];

        tr_z_zzzz_xxxxyz[i] = 3.0 * tr_z_zz_xxxxyz[i] * fe_0 + tr_z_zzz_xxxxy[i] * fe_0 + ts_zzz_xxxxyz[i] * fe_0 + tr_z_zzz_xxxxyz[i] * pa_z[i];

        tr_z_zzzz_xxxxzz[i] = 3.0 * tr_z_zz_xxxxzz[i] * fe_0 + 2.0 * tr_z_zzz_xxxxz[i] * fe_0 + ts_zzz_xxxxzz[i] * fe_0 + tr_z_zzz_xxxxzz[i] * pa_z[i];

        tr_z_zzzz_xxxyyy[i] = 3.0 * tr_z_zz_xxxyyy[i] * fe_0 + ts_zzz_xxxyyy[i] * fe_0 + tr_z_zzz_xxxyyy[i] * pa_z[i];

        tr_z_zzzz_xxxyyz[i] = 3.0 * tr_z_zz_xxxyyz[i] * fe_0 + tr_z_zzz_xxxyy[i] * fe_0 + ts_zzz_xxxyyz[i] * fe_0 + tr_z_zzz_xxxyyz[i] * pa_z[i];

        tr_z_zzzz_xxxyzz[i] = 3.0 * tr_z_zz_xxxyzz[i] * fe_0 + 2.0 * tr_z_zzz_xxxyz[i] * fe_0 + ts_zzz_xxxyzz[i] * fe_0 + tr_z_zzz_xxxyzz[i] * pa_z[i];

        tr_z_zzzz_xxxzzz[i] = 3.0 * tr_z_zz_xxxzzz[i] * fe_0 + 3.0 * tr_z_zzz_xxxzz[i] * fe_0 + ts_zzz_xxxzzz[i] * fe_0 + tr_z_zzz_xxxzzz[i] * pa_z[i];

        tr_z_zzzz_xxyyyy[i] = 3.0 * tr_z_zz_xxyyyy[i] * fe_0 + ts_zzz_xxyyyy[i] * fe_0 + tr_z_zzz_xxyyyy[i] * pa_z[i];

        tr_z_zzzz_xxyyyz[i] = 3.0 * tr_z_zz_xxyyyz[i] * fe_0 + tr_z_zzz_xxyyy[i] * fe_0 + ts_zzz_xxyyyz[i] * fe_0 + tr_z_zzz_xxyyyz[i] * pa_z[i];

        tr_z_zzzz_xxyyzz[i] = 3.0 * tr_z_zz_xxyyzz[i] * fe_0 + 2.0 * tr_z_zzz_xxyyz[i] * fe_0 + ts_zzz_xxyyzz[i] * fe_0 + tr_z_zzz_xxyyzz[i] * pa_z[i];

        tr_z_zzzz_xxyzzz[i] = 3.0 * tr_z_zz_xxyzzz[i] * fe_0 + 3.0 * tr_z_zzz_xxyzz[i] * fe_0 + ts_zzz_xxyzzz[i] * fe_0 + tr_z_zzz_xxyzzz[i] * pa_z[i];

        tr_z_zzzz_xxzzzz[i] = 3.0 * tr_z_zz_xxzzzz[i] * fe_0 + 4.0 * tr_z_zzz_xxzzz[i] * fe_0 + ts_zzz_xxzzzz[i] * fe_0 + tr_z_zzz_xxzzzz[i] * pa_z[i];

        tr_z_zzzz_xyyyyy[i] = 3.0 * tr_z_zz_xyyyyy[i] * fe_0 + ts_zzz_xyyyyy[i] * fe_0 + tr_z_zzz_xyyyyy[i] * pa_z[i];

        tr_z_zzzz_xyyyyz[i] = 3.0 * tr_z_zz_xyyyyz[i] * fe_0 + tr_z_zzz_xyyyy[i] * fe_0 + ts_zzz_xyyyyz[i] * fe_0 + tr_z_zzz_xyyyyz[i] * pa_z[i];

        tr_z_zzzz_xyyyzz[i] = 3.0 * tr_z_zz_xyyyzz[i] * fe_0 + 2.0 * tr_z_zzz_xyyyz[i] * fe_0 + ts_zzz_xyyyzz[i] * fe_0 + tr_z_zzz_xyyyzz[i] * pa_z[i];

        tr_z_zzzz_xyyzzz[i] = 3.0 * tr_z_zz_xyyzzz[i] * fe_0 + 3.0 * tr_z_zzz_xyyzz[i] * fe_0 + ts_zzz_xyyzzz[i] * fe_0 + tr_z_zzz_xyyzzz[i] * pa_z[i];

        tr_z_zzzz_xyzzzz[i] = 3.0 * tr_z_zz_xyzzzz[i] * fe_0 + 4.0 * tr_z_zzz_xyzzz[i] * fe_0 + ts_zzz_xyzzzz[i] * fe_0 + tr_z_zzz_xyzzzz[i] * pa_z[i];

        tr_z_zzzz_xzzzzz[i] = 3.0 * tr_z_zz_xzzzzz[i] * fe_0 + 5.0 * tr_z_zzz_xzzzz[i] * fe_0 + ts_zzz_xzzzzz[i] * fe_0 + tr_z_zzz_xzzzzz[i] * pa_z[i];

        tr_z_zzzz_yyyyyy[i] = 3.0 * tr_z_zz_yyyyyy[i] * fe_0 + ts_zzz_yyyyyy[i] * fe_0 + tr_z_zzz_yyyyyy[i] * pa_z[i];

        tr_z_zzzz_yyyyyz[i] = 3.0 * tr_z_zz_yyyyyz[i] * fe_0 + tr_z_zzz_yyyyy[i] * fe_0 + ts_zzz_yyyyyz[i] * fe_0 + tr_z_zzz_yyyyyz[i] * pa_z[i];

        tr_z_zzzz_yyyyzz[i] = 3.0 * tr_z_zz_yyyyzz[i] * fe_0 + 2.0 * tr_z_zzz_yyyyz[i] * fe_0 + ts_zzz_yyyyzz[i] * fe_0 + tr_z_zzz_yyyyzz[i] * pa_z[i];

        tr_z_zzzz_yyyzzz[i] = 3.0 * tr_z_zz_yyyzzz[i] * fe_0 + 3.0 * tr_z_zzz_yyyzz[i] * fe_0 + ts_zzz_yyyzzz[i] * fe_0 + tr_z_zzz_yyyzzz[i] * pa_z[i];

        tr_z_zzzz_yyzzzz[i] = 3.0 * tr_z_zz_yyzzzz[i] * fe_0 + 4.0 * tr_z_zzz_yyzzz[i] * fe_0 + ts_zzz_yyzzzz[i] * fe_0 + tr_z_zzz_yyzzzz[i] * pa_z[i];

        tr_z_zzzz_yzzzzz[i] = 3.0 * tr_z_zz_yzzzzz[i] * fe_0 + 5.0 * tr_z_zzz_yzzzz[i] * fe_0 + ts_zzz_yzzzzz[i] * fe_0 + tr_z_zzz_yzzzzz[i] * pa_z[i];

        tr_z_zzzz_zzzzzz[i] = 3.0 * tr_z_zz_zzzzzz[i] * fe_0 + 6.0 * tr_z_zzz_zzzzz[i] * fe_0 + ts_zzz_zzzzzz[i] * fe_0 + tr_z_zzz_zzzzzz[i] * pa_z[i];
    }

}

} // diprec namespace

