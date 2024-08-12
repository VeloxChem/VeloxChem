#include "ElectricDipoleMomentumPrimRecIH.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_ih(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_ih,
                                      const size_t idx_dip_gh,
                                      const size_t idx_dip_hg,
                                      const size_t idx_ovl_hh,
                                      const size_t idx_dip_hh,
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

    // Set up components of auxiliary buffer : GH

    auto tr_x_xxxx_xxxxx = pbuffer.data(idx_dip_gh);

    auto tr_x_xxxx_xxxxy = pbuffer.data(idx_dip_gh + 1);

    auto tr_x_xxxx_xxxxz = pbuffer.data(idx_dip_gh + 2);

    auto tr_x_xxxx_xxxyy = pbuffer.data(idx_dip_gh + 3);

    auto tr_x_xxxx_xxxyz = pbuffer.data(idx_dip_gh + 4);

    auto tr_x_xxxx_xxxzz = pbuffer.data(idx_dip_gh + 5);

    auto tr_x_xxxx_xxyyy = pbuffer.data(idx_dip_gh + 6);

    auto tr_x_xxxx_xxyyz = pbuffer.data(idx_dip_gh + 7);

    auto tr_x_xxxx_xxyzz = pbuffer.data(idx_dip_gh + 8);

    auto tr_x_xxxx_xxzzz = pbuffer.data(idx_dip_gh + 9);

    auto tr_x_xxxx_xyyyy = pbuffer.data(idx_dip_gh + 10);

    auto tr_x_xxxx_xyyyz = pbuffer.data(idx_dip_gh + 11);

    auto tr_x_xxxx_xyyzz = pbuffer.data(idx_dip_gh + 12);

    auto tr_x_xxxx_xyzzz = pbuffer.data(idx_dip_gh + 13);

    auto tr_x_xxxx_xzzzz = pbuffer.data(idx_dip_gh + 14);

    auto tr_x_xxxx_yyyyy = pbuffer.data(idx_dip_gh + 15);

    auto tr_x_xxxx_yyyyz = pbuffer.data(idx_dip_gh + 16);

    auto tr_x_xxxx_yyyzz = pbuffer.data(idx_dip_gh + 17);

    auto tr_x_xxxx_yyzzz = pbuffer.data(idx_dip_gh + 18);

    auto tr_x_xxxx_yzzzz = pbuffer.data(idx_dip_gh + 19);

    auto tr_x_xxxx_zzzzz = pbuffer.data(idx_dip_gh + 20);

    auto tr_x_xxxy_xxxxx = pbuffer.data(idx_dip_gh + 21);

    auto tr_x_xxxy_xxxxy = pbuffer.data(idx_dip_gh + 22);

    auto tr_x_xxxy_xxxxz = pbuffer.data(idx_dip_gh + 23);

    auto tr_x_xxxy_xxxyy = pbuffer.data(idx_dip_gh + 24);

    auto tr_x_xxxy_xxxyz = pbuffer.data(idx_dip_gh + 25);

    auto tr_x_xxxy_xxxzz = pbuffer.data(idx_dip_gh + 26);

    auto tr_x_xxxy_xxyyy = pbuffer.data(idx_dip_gh + 27);

    auto tr_x_xxxy_xxyyz = pbuffer.data(idx_dip_gh + 28);

    auto tr_x_xxxy_xxyzz = pbuffer.data(idx_dip_gh + 29);

    auto tr_x_xxxy_xxzzz = pbuffer.data(idx_dip_gh + 30);

    auto tr_x_xxxy_xyyyy = pbuffer.data(idx_dip_gh + 31);

    auto tr_x_xxxy_xyyyz = pbuffer.data(idx_dip_gh + 32);

    auto tr_x_xxxy_xyyzz = pbuffer.data(idx_dip_gh + 33);

    auto tr_x_xxxy_xyzzz = pbuffer.data(idx_dip_gh + 34);

    auto tr_x_xxxy_xzzzz = pbuffer.data(idx_dip_gh + 35);

    auto tr_x_xxxy_zzzzz = pbuffer.data(idx_dip_gh + 41);

    auto tr_x_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 42);

    auto tr_x_xxxz_xxxxy = pbuffer.data(idx_dip_gh + 43);

    auto tr_x_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 44);

    auto tr_x_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 45);

    auto tr_x_xxxz_xxxyz = pbuffer.data(idx_dip_gh + 46);

    auto tr_x_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 47);

    auto tr_x_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 48);

    auto tr_x_xxxz_xxyyz = pbuffer.data(idx_dip_gh + 49);

    auto tr_x_xxxz_xxyzz = pbuffer.data(idx_dip_gh + 50);

    auto tr_x_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 51);

    auto tr_x_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 52);

    auto tr_x_xxxz_xyyyz = pbuffer.data(idx_dip_gh + 53);

    auto tr_x_xxxz_xyyzz = pbuffer.data(idx_dip_gh + 54);

    auto tr_x_xxxz_xyzzz = pbuffer.data(idx_dip_gh + 55);

    auto tr_x_xxxz_xzzzz = pbuffer.data(idx_dip_gh + 56);

    auto tr_x_xxxz_yyyyy = pbuffer.data(idx_dip_gh + 57);

    auto tr_x_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 62);

    auto tr_x_xxyy_xxxxx = pbuffer.data(idx_dip_gh + 63);

    auto tr_x_xxyy_xxxxy = pbuffer.data(idx_dip_gh + 64);

    auto tr_x_xxyy_xxxxz = pbuffer.data(idx_dip_gh + 65);

    auto tr_x_xxyy_xxxyy = pbuffer.data(idx_dip_gh + 66);

    auto tr_x_xxyy_xxxyz = pbuffer.data(idx_dip_gh + 67);

    auto tr_x_xxyy_xxxzz = pbuffer.data(idx_dip_gh + 68);

    auto tr_x_xxyy_xxyyy = pbuffer.data(idx_dip_gh + 69);

    auto tr_x_xxyy_xxyyz = pbuffer.data(idx_dip_gh + 70);

    auto tr_x_xxyy_xxyzz = pbuffer.data(idx_dip_gh + 71);

    auto tr_x_xxyy_xxzzz = pbuffer.data(idx_dip_gh + 72);

    auto tr_x_xxyy_xyyyy = pbuffer.data(idx_dip_gh + 73);

    auto tr_x_xxyy_xyyyz = pbuffer.data(idx_dip_gh + 74);

    auto tr_x_xxyy_xyyzz = pbuffer.data(idx_dip_gh + 75);

    auto tr_x_xxyy_xyzzz = pbuffer.data(idx_dip_gh + 76);

    auto tr_x_xxyy_xzzzz = pbuffer.data(idx_dip_gh + 77);

    auto tr_x_xxyy_yyyyy = pbuffer.data(idx_dip_gh + 78);

    auto tr_x_xxyy_yyyyz = pbuffer.data(idx_dip_gh + 79);

    auto tr_x_xxyy_yyyzz = pbuffer.data(idx_dip_gh + 80);

    auto tr_x_xxyy_yyzzz = pbuffer.data(idx_dip_gh + 81);

    auto tr_x_xxyy_yzzzz = pbuffer.data(idx_dip_gh + 82);

    auto tr_x_xxyy_zzzzz = pbuffer.data(idx_dip_gh + 83);

    auto tr_x_xxyz_xxxxz = pbuffer.data(idx_dip_gh + 86);

    auto tr_x_xxyz_xxxzz = pbuffer.data(idx_dip_gh + 89);

    auto tr_x_xxyz_xxzzz = pbuffer.data(idx_dip_gh + 93);

    auto tr_x_xxyz_xzzzz = pbuffer.data(idx_dip_gh + 98);

    auto tr_x_xxyz_zzzzz = pbuffer.data(idx_dip_gh + 104);

    auto tr_x_xxzz_xxxxx = pbuffer.data(idx_dip_gh + 105);

    auto tr_x_xxzz_xxxxy = pbuffer.data(idx_dip_gh + 106);

    auto tr_x_xxzz_xxxxz = pbuffer.data(idx_dip_gh + 107);

    auto tr_x_xxzz_xxxyy = pbuffer.data(idx_dip_gh + 108);

    auto tr_x_xxzz_xxxyz = pbuffer.data(idx_dip_gh + 109);

    auto tr_x_xxzz_xxxzz = pbuffer.data(idx_dip_gh + 110);

    auto tr_x_xxzz_xxyyy = pbuffer.data(idx_dip_gh + 111);

    auto tr_x_xxzz_xxyyz = pbuffer.data(idx_dip_gh + 112);

    auto tr_x_xxzz_xxyzz = pbuffer.data(idx_dip_gh + 113);

    auto tr_x_xxzz_xxzzz = pbuffer.data(idx_dip_gh + 114);

    auto tr_x_xxzz_xyyyy = pbuffer.data(idx_dip_gh + 115);

    auto tr_x_xxzz_xyyyz = pbuffer.data(idx_dip_gh + 116);

    auto tr_x_xxzz_xyyzz = pbuffer.data(idx_dip_gh + 117);

    auto tr_x_xxzz_xyzzz = pbuffer.data(idx_dip_gh + 118);

    auto tr_x_xxzz_xzzzz = pbuffer.data(idx_dip_gh + 119);

    auto tr_x_xxzz_yyyyy = pbuffer.data(idx_dip_gh + 120);

    auto tr_x_xxzz_yyyyz = pbuffer.data(idx_dip_gh + 121);

    auto tr_x_xxzz_yyyzz = pbuffer.data(idx_dip_gh + 122);

    auto tr_x_xxzz_yyzzz = pbuffer.data(idx_dip_gh + 123);

    auto tr_x_xxzz_yzzzz = pbuffer.data(idx_dip_gh + 124);

    auto tr_x_xxzz_zzzzz = pbuffer.data(idx_dip_gh + 125);

    auto tr_x_xyyy_xxxxx = pbuffer.data(idx_dip_gh + 126);

    auto tr_x_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 127);

    auto tr_x_xyyy_xxxxz = pbuffer.data(idx_dip_gh + 128);

    auto tr_x_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 129);

    auto tr_x_xyyy_xxxzz = pbuffer.data(idx_dip_gh + 131);

    auto tr_x_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 132);

    auto tr_x_xyyy_xxzzz = pbuffer.data(idx_dip_gh + 135);

    auto tr_x_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 136);

    auto tr_x_xyyy_xzzzz = pbuffer.data(idx_dip_gh + 140);

    auto tr_x_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 141);

    auto tr_x_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 142);

    auto tr_x_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 143);

    auto tr_x_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 144);

    auto tr_x_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 145);

    auto tr_x_xyyz_xxxxy = pbuffer.data(idx_dip_gh + 148);

    auto tr_x_xyyz_xxxxz = pbuffer.data(idx_dip_gh + 149);

    auto tr_x_xyyz_xxxyy = pbuffer.data(idx_dip_gh + 150);

    auto tr_x_xyyz_xxxzz = pbuffer.data(idx_dip_gh + 152);

    auto tr_x_xyyz_xxyyy = pbuffer.data(idx_dip_gh + 153);

    auto tr_x_xyyz_xxzzz = pbuffer.data(idx_dip_gh + 156);

    auto tr_x_xyyz_xyyyy = pbuffer.data(idx_dip_gh + 157);

    auto tr_x_xyyz_xzzzz = pbuffer.data(idx_dip_gh + 161);

    auto tr_x_xyzz_xxxxx = pbuffer.data(idx_dip_gh + 168);

    auto tr_x_xyzz_xxxxz = pbuffer.data(idx_dip_gh + 170);

    auto tr_x_xyzz_xxxzz = pbuffer.data(idx_dip_gh + 173);

    auto tr_x_xyzz_xxzzz = pbuffer.data(idx_dip_gh + 177);

    auto tr_x_xyzz_xzzzz = pbuffer.data(idx_dip_gh + 182);

    auto tr_x_xzzz_xxxxx = pbuffer.data(idx_dip_gh + 189);

    auto tr_x_xzzz_xxxxy = pbuffer.data(idx_dip_gh + 190);

    auto tr_x_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 191);

    auto tr_x_xzzz_xxxyy = pbuffer.data(idx_dip_gh + 192);

    auto tr_x_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 194);

    auto tr_x_xzzz_xxyyy = pbuffer.data(idx_dip_gh + 195);

    auto tr_x_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 198);

    auto tr_x_xzzz_xyyyy = pbuffer.data(idx_dip_gh + 199);

    auto tr_x_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 203);

    auto tr_x_xzzz_yyyyz = pbuffer.data(idx_dip_gh + 205);

    auto tr_x_xzzz_yyyzz = pbuffer.data(idx_dip_gh + 206);

    auto tr_x_xzzz_yyzzz = pbuffer.data(idx_dip_gh + 207);

    auto tr_x_xzzz_yzzzz = pbuffer.data(idx_dip_gh + 208);

    auto tr_x_xzzz_zzzzz = pbuffer.data(idx_dip_gh + 209);

    auto tr_x_yyyy_xxxxx = pbuffer.data(idx_dip_gh + 210);

    auto tr_x_yyyy_xxxxy = pbuffer.data(idx_dip_gh + 211);

    auto tr_x_yyyy_xxxxz = pbuffer.data(idx_dip_gh + 212);

    auto tr_x_yyyy_xxxyy = pbuffer.data(idx_dip_gh + 213);

    auto tr_x_yyyy_xxxyz = pbuffer.data(idx_dip_gh + 214);

    auto tr_x_yyyy_xxxzz = pbuffer.data(idx_dip_gh + 215);

    auto tr_x_yyyy_xxyyy = pbuffer.data(idx_dip_gh + 216);

    auto tr_x_yyyy_xxyyz = pbuffer.data(idx_dip_gh + 217);

    auto tr_x_yyyy_xxyzz = pbuffer.data(idx_dip_gh + 218);

    auto tr_x_yyyy_xxzzz = pbuffer.data(idx_dip_gh + 219);

    auto tr_x_yyyy_xyyyy = pbuffer.data(idx_dip_gh + 220);

    auto tr_x_yyyy_xyyyz = pbuffer.data(idx_dip_gh + 221);

    auto tr_x_yyyy_xyyzz = pbuffer.data(idx_dip_gh + 222);

    auto tr_x_yyyy_xyzzz = pbuffer.data(idx_dip_gh + 223);

    auto tr_x_yyyy_xzzzz = pbuffer.data(idx_dip_gh + 224);

    auto tr_x_yyyy_yyyyy = pbuffer.data(idx_dip_gh + 225);

    auto tr_x_yyyy_yyyyz = pbuffer.data(idx_dip_gh + 226);

    auto tr_x_yyyy_yyyzz = pbuffer.data(idx_dip_gh + 227);

    auto tr_x_yyyy_yyzzz = pbuffer.data(idx_dip_gh + 228);

    auto tr_x_yyyy_yzzzz = pbuffer.data(idx_dip_gh + 229);

    auto tr_x_yyyy_zzzzz = pbuffer.data(idx_dip_gh + 230);

    auto tr_x_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 232);

    auto tr_x_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 233);

    auto tr_x_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 234);

    auto tr_x_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 236);

    auto tr_x_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 237);

    auto tr_x_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 240);

    auto tr_x_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 241);

    auto tr_x_yyyz_xzzzz = pbuffer.data(idx_dip_gh + 245);

    auto tr_x_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 246);

    auto tr_x_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 251);

    auto tr_x_yyzz_xxxxx = pbuffer.data(idx_dip_gh + 252);

    auto tr_x_yyzz_xxxxy = pbuffer.data(idx_dip_gh + 253);

    auto tr_x_yyzz_xxxxz = pbuffer.data(idx_dip_gh + 254);

    auto tr_x_yyzz_xxxyy = pbuffer.data(idx_dip_gh + 255);

    auto tr_x_yyzz_xxxyz = pbuffer.data(idx_dip_gh + 256);

    auto tr_x_yyzz_xxxzz = pbuffer.data(idx_dip_gh + 257);

    auto tr_x_yyzz_xxyyy = pbuffer.data(idx_dip_gh + 258);

    auto tr_x_yyzz_xxyyz = pbuffer.data(idx_dip_gh + 259);

    auto tr_x_yyzz_xxyzz = pbuffer.data(idx_dip_gh + 260);

    auto tr_x_yyzz_xxzzz = pbuffer.data(idx_dip_gh + 261);

    auto tr_x_yyzz_xyyyy = pbuffer.data(idx_dip_gh + 262);

    auto tr_x_yyzz_xyyyz = pbuffer.data(idx_dip_gh + 263);

    auto tr_x_yyzz_xyyzz = pbuffer.data(idx_dip_gh + 264);

    auto tr_x_yyzz_xyzzz = pbuffer.data(idx_dip_gh + 265);

    auto tr_x_yyzz_xzzzz = pbuffer.data(idx_dip_gh + 266);

    auto tr_x_yyzz_yyyyy = pbuffer.data(idx_dip_gh + 267);

    auto tr_x_yyzz_yyyyz = pbuffer.data(idx_dip_gh + 268);

    auto tr_x_yyzz_yyyzz = pbuffer.data(idx_dip_gh + 269);

    auto tr_x_yyzz_yyzzz = pbuffer.data(idx_dip_gh + 270);

    auto tr_x_yyzz_yzzzz = pbuffer.data(idx_dip_gh + 271);

    auto tr_x_yyzz_zzzzz = pbuffer.data(idx_dip_gh + 272);

    auto tr_x_yzzz_xxxxx = pbuffer.data(idx_dip_gh + 273);

    auto tr_x_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 275);

    auto tr_x_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 277);

    auto tr_x_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 278);

    auto tr_x_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 280);

    auto tr_x_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 281);

    auto tr_x_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 282);

    auto tr_x_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 284);

    auto tr_x_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 285);

    auto tr_x_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 286);

    auto tr_x_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 287);

    auto tr_x_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 289);

    auto tr_x_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 290);

    auto tr_x_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 291);

    auto tr_x_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 292);

    auto tr_x_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 293);

    auto tr_x_zzzz_xxxxx = pbuffer.data(idx_dip_gh + 294);

    auto tr_x_zzzz_xxxxy = pbuffer.data(idx_dip_gh + 295);

    auto tr_x_zzzz_xxxxz = pbuffer.data(idx_dip_gh + 296);

    auto tr_x_zzzz_xxxyy = pbuffer.data(idx_dip_gh + 297);

    auto tr_x_zzzz_xxxyz = pbuffer.data(idx_dip_gh + 298);

    auto tr_x_zzzz_xxxzz = pbuffer.data(idx_dip_gh + 299);

    auto tr_x_zzzz_xxyyy = pbuffer.data(idx_dip_gh + 300);

    auto tr_x_zzzz_xxyyz = pbuffer.data(idx_dip_gh + 301);

    auto tr_x_zzzz_xxyzz = pbuffer.data(idx_dip_gh + 302);

    auto tr_x_zzzz_xxzzz = pbuffer.data(idx_dip_gh + 303);

    auto tr_x_zzzz_xyyyy = pbuffer.data(idx_dip_gh + 304);

    auto tr_x_zzzz_xyyyz = pbuffer.data(idx_dip_gh + 305);

    auto tr_x_zzzz_xyyzz = pbuffer.data(idx_dip_gh + 306);

    auto tr_x_zzzz_xyzzz = pbuffer.data(idx_dip_gh + 307);

    auto tr_x_zzzz_xzzzz = pbuffer.data(idx_dip_gh + 308);

    auto tr_x_zzzz_yyyyy = pbuffer.data(idx_dip_gh + 309);

    auto tr_x_zzzz_yyyyz = pbuffer.data(idx_dip_gh + 310);

    auto tr_x_zzzz_yyyzz = pbuffer.data(idx_dip_gh + 311);

    auto tr_x_zzzz_yyzzz = pbuffer.data(idx_dip_gh + 312);

    auto tr_x_zzzz_yzzzz = pbuffer.data(idx_dip_gh + 313);

    auto tr_x_zzzz_zzzzz = pbuffer.data(idx_dip_gh + 314);

    auto tr_y_xxxx_xxxxx = pbuffer.data(idx_dip_gh + 315);

    auto tr_y_xxxx_xxxxy = pbuffer.data(idx_dip_gh + 316);

    auto tr_y_xxxx_xxxxz = pbuffer.data(idx_dip_gh + 317);

    auto tr_y_xxxx_xxxyy = pbuffer.data(idx_dip_gh + 318);

    auto tr_y_xxxx_xxxyz = pbuffer.data(idx_dip_gh + 319);

    auto tr_y_xxxx_xxxzz = pbuffer.data(idx_dip_gh + 320);

    auto tr_y_xxxx_xxyyy = pbuffer.data(idx_dip_gh + 321);

    auto tr_y_xxxx_xxyyz = pbuffer.data(idx_dip_gh + 322);

    auto tr_y_xxxx_xxyzz = pbuffer.data(idx_dip_gh + 323);

    auto tr_y_xxxx_xxzzz = pbuffer.data(idx_dip_gh + 324);

    auto tr_y_xxxx_xyyyy = pbuffer.data(idx_dip_gh + 325);

    auto tr_y_xxxx_xyyyz = pbuffer.data(idx_dip_gh + 326);

    auto tr_y_xxxx_xyyzz = pbuffer.data(idx_dip_gh + 327);

    auto tr_y_xxxx_xyzzz = pbuffer.data(idx_dip_gh + 328);

    auto tr_y_xxxx_xzzzz = pbuffer.data(idx_dip_gh + 329);

    auto tr_y_xxxx_yyyyy = pbuffer.data(idx_dip_gh + 330);

    auto tr_y_xxxx_yyyyz = pbuffer.data(idx_dip_gh + 331);

    auto tr_y_xxxx_yyyzz = pbuffer.data(idx_dip_gh + 332);

    auto tr_y_xxxx_yyzzz = pbuffer.data(idx_dip_gh + 333);

    auto tr_y_xxxx_yzzzz = pbuffer.data(idx_dip_gh + 334);

    auto tr_y_xxxx_zzzzz = pbuffer.data(idx_dip_gh + 335);

    auto tr_y_xxxy_xxxxy = pbuffer.data(idx_dip_gh + 337);

    auto tr_y_xxxy_xxxyy = pbuffer.data(idx_dip_gh + 339);

    auto tr_y_xxxy_xxxyz = pbuffer.data(idx_dip_gh + 340);

    auto tr_y_xxxy_xxyyy = pbuffer.data(idx_dip_gh + 342);

    auto tr_y_xxxy_xxyyz = pbuffer.data(idx_dip_gh + 343);

    auto tr_y_xxxy_xxyzz = pbuffer.data(idx_dip_gh + 344);

    auto tr_y_xxxy_xyyyy = pbuffer.data(idx_dip_gh + 346);

    auto tr_y_xxxy_xyyyz = pbuffer.data(idx_dip_gh + 347);

    auto tr_y_xxxy_xyyzz = pbuffer.data(idx_dip_gh + 348);

    auto tr_y_xxxy_xyzzz = pbuffer.data(idx_dip_gh + 349);

    auto tr_y_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 351);

    auto tr_y_xxxy_yyyyz = pbuffer.data(idx_dip_gh + 352);

    auto tr_y_xxxy_yyyzz = pbuffer.data(idx_dip_gh + 353);

    auto tr_y_xxxy_yyzzz = pbuffer.data(idx_dip_gh + 354);

    auto tr_y_xxxy_yzzzz = pbuffer.data(idx_dip_gh + 355);

    auto tr_y_xxxy_zzzzz = pbuffer.data(idx_dip_gh + 356);

    auto tr_y_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 357);

    auto tr_y_xxxz_xxxxy = pbuffer.data(idx_dip_gh + 358);

    auto tr_y_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 360);

    auto tr_y_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 363);

    auto tr_y_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 367);

    auto tr_y_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 373);

    auto tr_y_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 374);

    auto tr_y_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 375);

    auto tr_y_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 376);

    auto tr_y_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 377);

    auto tr_y_xxyy_xxxxx = pbuffer.data(idx_dip_gh + 378);

    auto tr_y_xxyy_xxxxy = pbuffer.data(idx_dip_gh + 379);

    auto tr_y_xxyy_xxxxz = pbuffer.data(idx_dip_gh + 380);

    auto tr_y_xxyy_xxxyy = pbuffer.data(idx_dip_gh + 381);

    auto tr_y_xxyy_xxxyz = pbuffer.data(idx_dip_gh + 382);

    auto tr_y_xxyy_xxxzz = pbuffer.data(idx_dip_gh + 383);

    auto tr_y_xxyy_xxyyy = pbuffer.data(idx_dip_gh + 384);

    auto tr_y_xxyy_xxyyz = pbuffer.data(idx_dip_gh + 385);

    auto tr_y_xxyy_xxyzz = pbuffer.data(idx_dip_gh + 386);

    auto tr_y_xxyy_xxzzz = pbuffer.data(idx_dip_gh + 387);

    auto tr_y_xxyy_xyyyy = pbuffer.data(idx_dip_gh + 388);

    auto tr_y_xxyy_xyyyz = pbuffer.data(idx_dip_gh + 389);

    auto tr_y_xxyy_xyyzz = pbuffer.data(idx_dip_gh + 390);

    auto tr_y_xxyy_xyzzz = pbuffer.data(idx_dip_gh + 391);

    auto tr_y_xxyy_xzzzz = pbuffer.data(idx_dip_gh + 392);

    auto tr_y_xxyy_yyyyy = pbuffer.data(idx_dip_gh + 393);

    auto tr_y_xxyy_yyyyz = pbuffer.data(idx_dip_gh + 394);

    auto tr_y_xxyy_yyyzz = pbuffer.data(idx_dip_gh + 395);

    auto tr_y_xxyy_yyzzz = pbuffer.data(idx_dip_gh + 396);

    auto tr_y_xxyy_yzzzz = pbuffer.data(idx_dip_gh + 397);

    auto tr_y_xxyy_zzzzz = pbuffer.data(idx_dip_gh + 398);

    auto tr_y_xxyz_xxxxy = pbuffer.data(idx_dip_gh + 400);

    auto tr_y_xxyz_xxxyy = pbuffer.data(idx_dip_gh + 402);

    auto tr_y_xxyz_xxyyy = pbuffer.data(idx_dip_gh + 405);

    auto tr_y_xxyz_xyyyy = pbuffer.data(idx_dip_gh + 409);

    auto tr_y_xxyz_yyyyz = pbuffer.data(idx_dip_gh + 415);

    auto tr_y_xxyz_yyyzz = pbuffer.data(idx_dip_gh + 416);

    auto tr_y_xxyz_yyzzz = pbuffer.data(idx_dip_gh + 417);

    auto tr_y_xxyz_yzzzz = pbuffer.data(idx_dip_gh + 418);

    auto tr_y_xxyz_zzzzz = pbuffer.data(idx_dip_gh + 419);

    auto tr_y_xxzz_xxxxx = pbuffer.data(idx_dip_gh + 420);

    auto tr_y_xxzz_xxxxy = pbuffer.data(idx_dip_gh + 421);

    auto tr_y_xxzz_xxxxz = pbuffer.data(idx_dip_gh + 422);

    auto tr_y_xxzz_xxxyy = pbuffer.data(idx_dip_gh + 423);

    auto tr_y_xxzz_xxxyz = pbuffer.data(idx_dip_gh + 424);

    auto tr_y_xxzz_xxxzz = pbuffer.data(idx_dip_gh + 425);

    auto tr_y_xxzz_xxyyy = pbuffer.data(idx_dip_gh + 426);

    auto tr_y_xxzz_xxyyz = pbuffer.data(idx_dip_gh + 427);

    auto tr_y_xxzz_xxyzz = pbuffer.data(idx_dip_gh + 428);

    auto tr_y_xxzz_xxzzz = pbuffer.data(idx_dip_gh + 429);

    auto tr_y_xxzz_xyyyy = pbuffer.data(idx_dip_gh + 430);

    auto tr_y_xxzz_xyyyz = pbuffer.data(idx_dip_gh + 431);

    auto tr_y_xxzz_xyyzz = pbuffer.data(idx_dip_gh + 432);

    auto tr_y_xxzz_xyzzz = pbuffer.data(idx_dip_gh + 433);

    auto tr_y_xxzz_xzzzz = pbuffer.data(idx_dip_gh + 434);

    auto tr_y_xxzz_yyyyy = pbuffer.data(idx_dip_gh + 435);

    auto tr_y_xxzz_yyyyz = pbuffer.data(idx_dip_gh + 436);

    auto tr_y_xxzz_yyyzz = pbuffer.data(idx_dip_gh + 437);

    auto tr_y_xxzz_yyzzz = pbuffer.data(idx_dip_gh + 438);

    auto tr_y_xxzz_yzzzz = pbuffer.data(idx_dip_gh + 439);

    auto tr_y_xxzz_zzzzz = pbuffer.data(idx_dip_gh + 440);

    auto tr_y_xyyy_xxxxx = pbuffer.data(idx_dip_gh + 441);

    auto tr_y_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 442);

    auto tr_y_xyyy_xxxxz = pbuffer.data(idx_dip_gh + 443);

    auto tr_y_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 444);

    auto tr_y_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 445);

    auto tr_y_xyyy_xxxzz = pbuffer.data(idx_dip_gh + 446);

    auto tr_y_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 447);

    auto tr_y_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 448);

    auto tr_y_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 449);

    auto tr_y_xyyy_xxzzz = pbuffer.data(idx_dip_gh + 450);

    auto tr_y_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 451);

    auto tr_y_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 452);

    auto tr_y_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 453);

    auto tr_y_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 454);

    auto tr_y_xyyy_xzzzz = pbuffer.data(idx_dip_gh + 455);

    auto tr_y_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 456);

    auto tr_y_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 457);

    auto tr_y_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 458);

    auto tr_y_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 459);

    auto tr_y_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 460);

    auto tr_y_xyyy_zzzzz = pbuffer.data(idx_dip_gh + 461);

    auto tr_y_xyyz_yyyyz = pbuffer.data(idx_dip_gh + 478);

    auto tr_y_xyyz_yyyzz = pbuffer.data(idx_dip_gh + 479);

    auto tr_y_xyyz_yyzzz = pbuffer.data(idx_dip_gh + 480);

    auto tr_y_xyyz_yzzzz = pbuffer.data(idx_dip_gh + 481);

    auto tr_y_xyyz_zzzzz = pbuffer.data(idx_dip_gh + 482);

    auto tr_y_xyzz_xxxyz = pbuffer.data(idx_dip_gh + 487);

    auto tr_y_xyzz_xxyyz = pbuffer.data(idx_dip_gh + 490);

    auto tr_y_xyzz_xxyzz = pbuffer.data(idx_dip_gh + 491);

    auto tr_y_xyzz_xyyyz = pbuffer.data(idx_dip_gh + 494);

    auto tr_y_xyzz_xyyzz = pbuffer.data(idx_dip_gh + 495);

    auto tr_y_xyzz_xyzzz = pbuffer.data(idx_dip_gh + 496);

    auto tr_y_xyzz_yyyyy = pbuffer.data(idx_dip_gh + 498);

    auto tr_y_xyzz_yyyyz = pbuffer.data(idx_dip_gh + 499);

    auto tr_y_xyzz_yyyzz = pbuffer.data(idx_dip_gh + 500);

    auto tr_y_xyzz_yyzzz = pbuffer.data(idx_dip_gh + 501);

    auto tr_y_xyzz_yzzzz = pbuffer.data(idx_dip_gh + 502);

    auto tr_y_xyzz_zzzzz = pbuffer.data(idx_dip_gh + 503);

    auto tr_y_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 506);

    auto tr_y_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 508);

    auto tr_y_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 509);

    auto tr_y_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 511);

    auto tr_y_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 512);

    auto tr_y_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 513);

    auto tr_y_xzzz_xyyyz = pbuffer.data(idx_dip_gh + 515);

    auto tr_y_xzzz_xyyzz = pbuffer.data(idx_dip_gh + 516);

    auto tr_y_xzzz_xyzzz = pbuffer.data(idx_dip_gh + 517);

    auto tr_y_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 518);

    auto tr_y_xzzz_yyyyy = pbuffer.data(idx_dip_gh + 519);

    auto tr_y_xzzz_yyyyz = pbuffer.data(idx_dip_gh + 520);

    auto tr_y_xzzz_yyyzz = pbuffer.data(idx_dip_gh + 521);

    auto tr_y_xzzz_yyzzz = pbuffer.data(idx_dip_gh + 522);

    auto tr_y_xzzz_yzzzz = pbuffer.data(idx_dip_gh + 523);

    auto tr_y_xzzz_zzzzz = pbuffer.data(idx_dip_gh + 524);

    auto tr_y_yyyy_xxxxx = pbuffer.data(idx_dip_gh + 525);

    auto tr_y_yyyy_xxxxy = pbuffer.data(idx_dip_gh + 526);

    auto tr_y_yyyy_xxxxz = pbuffer.data(idx_dip_gh + 527);

    auto tr_y_yyyy_xxxyy = pbuffer.data(idx_dip_gh + 528);

    auto tr_y_yyyy_xxxyz = pbuffer.data(idx_dip_gh + 529);

    auto tr_y_yyyy_xxxzz = pbuffer.data(idx_dip_gh + 530);

    auto tr_y_yyyy_xxyyy = pbuffer.data(idx_dip_gh + 531);

    auto tr_y_yyyy_xxyyz = pbuffer.data(idx_dip_gh + 532);

    auto tr_y_yyyy_xxyzz = pbuffer.data(idx_dip_gh + 533);

    auto tr_y_yyyy_xxzzz = pbuffer.data(idx_dip_gh + 534);

    auto tr_y_yyyy_xyyyy = pbuffer.data(idx_dip_gh + 535);

    auto tr_y_yyyy_xyyyz = pbuffer.data(idx_dip_gh + 536);

    auto tr_y_yyyy_xyyzz = pbuffer.data(idx_dip_gh + 537);

    auto tr_y_yyyy_xyzzz = pbuffer.data(idx_dip_gh + 538);

    auto tr_y_yyyy_xzzzz = pbuffer.data(idx_dip_gh + 539);

    auto tr_y_yyyy_yyyyy = pbuffer.data(idx_dip_gh + 540);

    auto tr_y_yyyy_yyyyz = pbuffer.data(idx_dip_gh + 541);

    auto tr_y_yyyy_yyyzz = pbuffer.data(idx_dip_gh + 542);

    auto tr_y_yyyy_yyzzz = pbuffer.data(idx_dip_gh + 543);

    auto tr_y_yyyy_yzzzz = pbuffer.data(idx_dip_gh + 544);

    auto tr_y_yyyy_zzzzz = pbuffer.data(idx_dip_gh + 545);

    auto tr_y_yyyz_xxxxx = pbuffer.data(idx_dip_gh + 546);

    auto tr_y_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 547);

    auto tr_y_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 549);

    auto tr_y_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 550);

    auto tr_y_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 552);

    auto tr_y_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 553);

    auto tr_y_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 554);

    auto tr_y_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 556);

    auto tr_y_yyyz_xyyyz = pbuffer.data(idx_dip_gh + 557);

    auto tr_y_yyyz_xyyzz = pbuffer.data(idx_dip_gh + 558);

    auto tr_y_yyyz_xyzzz = pbuffer.data(idx_dip_gh + 559);

    auto tr_y_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 561);

    auto tr_y_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 562);

    auto tr_y_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 563);

    auto tr_y_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 564);

    auto tr_y_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 565);

    auto tr_y_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 566);

    auto tr_y_yyzz_xxxxx = pbuffer.data(idx_dip_gh + 567);

    auto tr_y_yyzz_xxxxy = pbuffer.data(idx_dip_gh + 568);

    auto tr_y_yyzz_xxxxz = pbuffer.data(idx_dip_gh + 569);

    auto tr_y_yyzz_xxxyy = pbuffer.data(idx_dip_gh + 570);

    auto tr_y_yyzz_xxxyz = pbuffer.data(idx_dip_gh + 571);

    auto tr_y_yyzz_xxxzz = pbuffer.data(idx_dip_gh + 572);

    auto tr_y_yyzz_xxyyy = pbuffer.data(idx_dip_gh + 573);

    auto tr_y_yyzz_xxyyz = pbuffer.data(idx_dip_gh + 574);

    auto tr_y_yyzz_xxyzz = pbuffer.data(idx_dip_gh + 575);

    auto tr_y_yyzz_xxzzz = pbuffer.data(idx_dip_gh + 576);

    auto tr_y_yyzz_xyyyy = pbuffer.data(idx_dip_gh + 577);

    auto tr_y_yyzz_xyyyz = pbuffer.data(idx_dip_gh + 578);

    auto tr_y_yyzz_xyyzz = pbuffer.data(idx_dip_gh + 579);

    auto tr_y_yyzz_xyzzz = pbuffer.data(idx_dip_gh + 580);

    auto tr_y_yyzz_xzzzz = pbuffer.data(idx_dip_gh + 581);

    auto tr_y_yyzz_yyyyy = pbuffer.data(idx_dip_gh + 582);

    auto tr_y_yyzz_yyyyz = pbuffer.data(idx_dip_gh + 583);

    auto tr_y_yyzz_yyyzz = pbuffer.data(idx_dip_gh + 584);

    auto tr_y_yyzz_yyzzz = pbuffer.data(idx_dip_gh + 585);

    auto tr_y_yyzz_yzzzz = pbuffer.data(idx_dip_gh + 586);

    auto tr_y_yyzz_zzzzz = pbuffer.data(idx_dip_gh + 587);

    auto tr_y_yzzz_xxxxy = pbuffer.data(idx_dip_gh + 589);

    auto tr_y_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 590);

    auto tr_y_yzzz_xxxyy = pbuffer.data(idx_dip_gh + 591);

    auto tr_y_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 592);

    auto tr_y_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 593);

    auto tr_y_yzzz_xxyyy = pbuffer.data(idx_dip_gh + 594);

    auto tr_y_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 595);

    auto tr_y_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 596);

    auto tr_y_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 597);

    auto tr_y_yzzz_xyyyy = pbuffer.data(idx_dip_gh + 598);

    auto tr_y_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 599);

    auto tr_y_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 600);

    auto tr_y_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 601);

    auto tr_y_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 602);

    auto tr_y_yzzz_yyyyy = pbuffer.data(idx_dip_gh + 603);

    auto tr_y_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 604);

    auto tr_y_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 605);

    auto tr_y_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 606);

    auto tr_y_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 607);

    auto tr_y_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 608);

    auto tr_y_zzzz_xxxxx = pbuffer.data(idx_dip_gh + 609);

    auto tr_y_zzzz_xxxxy = pbuffer.data(idx_dip_gh + 610);

    auto tr_y_zzzz_xxxxz = pbuffer.data(idx_dip_gh + 611);

    auto tr_y_zzzz_xxxyy = pbuffer.data(idx_dip_gh + 612);

    auto tr_y_zzzz_xxxyz = pbuffer.data(idx_dip_gh + 613);

    auto tr_y_zzzz_xxxzz = pbuffer.data(idx_dip_gh + 614);

    auto tr_y_zzzz_xxyyy = pbuffer.data(idx_dip_gh + 615);

    auto tr_y_zzzz_xxyyz = pbuffer.data(idx_dip_gh + 616);

    auto tr_y_zzzz_xxyzz = pbuffer.data(idx_dip_gh + 617);

    auto tr_y_zzzz_xxzzz = pbuffer.data(idx_dip_gh + 618);

    auto tr_y_zzzz_xyyyy = pbuffer.data(idx_dip_gh + 619);

    auto tr_y_zzzz_xyyyz = pbuffer.data(idx_dip_gh + 620);

    auto tr_y_zzzz_xyyzz = pbuffer.data(idx_dip_gh + 621);

    auto tr_y_zzzz_xyzzz = pbuffer.data(idx_dip_gh + 622);

    auto tr_y_zzzz_xzzzz = pbuffer.data(idx_dip_gh + 623);

    auto tr_y_zzzz_yyyyy = pbuffer.data(idx_dip_gh + 624);

    auto tr_y_zzzz_yyyyz = pbuffer.data(idx_dip_gh + 625);

    auto tr_y_zzzz_yyyzz = pbuffer.data(idx_dip_gh + 626);

    auto tr_y_zzzz_yyzzz = pbuffer.data(idx_dip_gh + 627);

    auto tr_y_zzzz_yzzzz = pbuffer.data(idx_dip_gh + 628);

    auto tr_y_zzzz_zzzzz = pbuffer.data(idx_dip_gh + 629);

    auto tr_z_xxxx_xxxxx = pbuffer.data(idx_dip_gh + 630);

    auto tr_z_xxxx_xxxxy = pbuffer.data(idx_dip_gh + 631);

    auto tr_z_xxxx_xxxxz = pbuffer.data(idx_dip_gh + 632);

    auto tr_z_xxxx_xxxyy = pbuffer.data(idx_dip_gh + 633);

    auto tr_z_xxxx_xxxyz = pbuffer.data(idx_dip_gh + 634);

    auto tr_z_xxxx_xxxzz = pbuffer.data(idx_dip_gh + 635);

    auto tr_z_xxxx_xxyyy = pbuffer.data(idx_dip_gh + 636);

    auto tr_z_xxxx_xxyyz = pbuffer.data(idx_dip_gh + 637);

    auto tr_z_xxxx_xxyzz = pbuffer.data(idx_dip_gh + 638);

    auto tr_z_xxxx_xxzzz = pbuffer.data(idx_dip_gh + 639);

    auto tr_z_xxxx_xyyyy = pbuffer.data(idx_dip_gh + 640);

    auto tr_z_xxxx_xyyyz = pbuffer.data(idx_dip_gh + 641);

    auto tr_z_xxxx_xyyzz = pbuffer.data(idx_dip_gh + 642);

    auto tr_z_xxxx_xyzzz = pbuffer.data(idx_dip_gh + 643);

    auto tr_z_xxxx_xzzzz = pbuffer.data(idx_dip_gh + 644);

    auto tr_z_xxxx_yyyyy = pbuffer.data(idx_dip_gh + 645);

    auto tr_z_xxxx_yyyyz = pbuffer.data(idx_dip_gh + 646);

    auto tr_z_xxxx_yyyzz = pbuffer.data(idx_dip_gh + 647);

    auto tr_z_xxxx_yyzzz = pbuffer.data(idx_dip_gh + 648);

    auto tr_z_xxxx_yzzzz = pbuffer.data(idx_dip_gh + 649);

    auto tr_z_xxxx_zzzzz = pbuffer.data(idx_dip_gh + 650);

    auto tr_z_xxxy_xxxxx = pbuffer.data(idx_dip_gh + 651);

    auto tr_z_xxxy_xxxxz = pbuffer.data(idx_dip_gh + 653);

    auto tr_z_xxxy_xxxzz = pbuffer.data(idx_dip_gh + 656);

    auto tr_z_xxxy_xxzzz = pbuffer.data(idx_dip_gh + 660);

    auto tr_z_xxxy_xzzzz = pbuffer.data(idx_dip_gh + 665);

    auto tr_z_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 666);

    auto tr_z_xxxy_yyyyz = pbuffer.data(idx_dip_gh + 667);

    auto tr_z_xxxy_yyyzz = pbuffer.data(idx_dip_gh + 668);

    auto tr_z_xxxy_yyzzz = pbuffer.data(idx_dip_gh + 669);

    auto tr_z_xxxy_yzzzz = pbuffer.data(idx_dip_gh + 670);

    auto tr_z_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 672);

    auto tr_z_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 674);

    auto tr_z_xxxz_xxxyz = pbuffer.data(idx_dip_gh + 676);

    auto tr_z_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 677);

    auto tr_z_xxxz_xxyyz = pbuffer.data(idx_dip_gh + 679);

    auto tr_z_xxxz_xxyzz = pbuffer.data(idx_dip_gh + 680);

    auto tr_z_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 681);

    auto tr_z_xxxz_xyyyz = pbuffer.data(idx_dip_gh + 683);

    auto tr_z_xxxz_xyyzz = pbuffer.data(idx_dip_gh + 684);

    auto tr_z_xxxz_xyzzz = pbuffer.data(idx_dip_gh + 685);

    auto tr_z_xxxz_xzzzz = pbuffer.data(idx_dip_gh + 686);

    auto tr_z_xxxz_yyyyy = pbuffer.data(idx_dip_gh + 687);

    auto tr_z_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 688);

    auto tr_z_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 689);

    auto tr_z_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 690);

    auto tr_z_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 691);

    auto tr_z_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 692);

    auto tr_z_xxyy_xxxxx = pbuffer.data(idx_dip_gh + 693);

    auto tr_z_xxyy_xxxxy = pbuffer.data(idx_dip_gh + 694);

    auto tr_z_xxyy_xxxxz = pbuffer.data(idx_dip_gh + 695);

    auto tr_z_xxyy_xxxyy = pbuffer.data(idx_dip_gh + 696);

    auto tr_z_xxyy_xxxyz = pbuffer.data(idx_dip_gh + 697);

    auto tr_z_xxyy_xxxzz = pbuffer.data(idx_dip_gh + 698);

    auto tr_z_xxyy_xxyyy = pbuffer.data(idx_dip_gh + 699);

    auto tr_z_xxyy_xxyyz = pbuffer.data(idx_dip_gh + 700);

    auto tr_z_xxyy_xxyzz = pbuffer.data(idx_dip_gh + 701);

    auto tr_z_xxyy_xxzzz = pbuffer.data(idx_dip_gh + 702);

    auto tr_z_xxyy_xyyyy = pbuffer.data(idx_dip_gh + 703);

    auto tr_z_xxyy_xyyyz = pbuffer.data(idx_dip_gh + 704);

    auto tr_z_xxyy_xyyzz = pbuffer.data(idx_dip_gh + 705);

    auto tr_z_xxyy_xyzzz = pbuffer.data(idx_dip_gh + 706);

    auto tr_z_xxyy_xzzzz = pbuffer.data(idx_dip_gh + 707);

    auto tr_z_xxyy_yyyyy = pbuffer.data(idx_dip_gh + 708);

    auto tr_z_xxyy_yyyyz = pbuffer.data(idx_dip_gh + 709);

    auto tr_z_xxyy_yyyzz = pbuffer.data(idx_dip_gh + 710);

    auto tr_z_xxyy_yyzzz = pbuffer.data(idx_dip_gh + 711);

    auto tr_z_xxyy_yzzzz = pbuffer.data(idx_dip_gh + 712);

    auto tr_z_xxyy_zzzzz = pbuffer.data(idx_dip_gh + 713);

    auto tr_z_xxyz_xxxxx = pbuffer.data(idx_dip_gh + 714);

    auto tr_z_xxyz_xxxxz = pbuffer.data(idx_dip_gh + 716);

    auto tr_z_xxyz_xxxzz = pbuffer.data(idx_dip_gh + 719);

    auto tr_z_xxyz_xxzzz = pbuffer.data(idx_dip_gh + 723);

    auto tr_z_xxyz_xzzzz = pbuffer.data(idx_dip_gh + 728);

    auto tr_z_xxyz_yyyyy = pbuffer.data(idx_dip_gh + 729);

    auto tr_z_xxyz_yyyyz = pbuffer.data(idx_dip_gh + 730);

    auto tr_z_xxyz_yyyzz = pbuffer.data(idx_dip_gh + 731);

    auto tr_z_xxyz_yyzzz = pbuffer.data(idx_dip_gh + 732);

    auto tr_z_xxyz_yzzzz = pbuffer.data(idx_dip_gh + 733);

    auto tr_z_xxzz_xxxxx = pbuffer.data(idx_dip_gh + 735);

    auto tr_z_xxzz_xxxxy = pbuffer.data(idx_dip_gh + 736);

    auto tr_z_xxzz_xxxxz = pbuffer.data(idx_dip_gh + 737);

    auto tr_z_xxzz_xxxyy = pbuffer.data(idx_dip_gh + 738);

    auto tr_z_xxzz_xxxyz = pbuffer.data(idx_dip_gh + 739);

    auto tr_z_xxzz_xxxzz = pbuffer.data(idx_dip_gh + 740);

    auto tr_z_xxzz_xxyyy = pbuffer.data(idx_dip_gh + 741);

    auto tr_z_xxzz_xxyyz = pbuffer.data(idx_dip_gh + 742);

    auto tr_z_xxzz_xxyzz = pbuffer.data(idx_dip_gh + 743);

    auto tr_z_xxzz_xxzzz = pbuffer.data(idx_dip_gh + 744);

    auto tr_z_xxzz_xyyyy = pbuffer.data(idx_dip_gh + 745);

    auto tr_z_xxzz_xyyyz = pbuffer.data(idx_dip_gh + 746);

    auto tr_z_xxzz_xyyzz = pbuffer.data(idx_dip_gh + 747);

    auto tr_z_xxzz_xyzzz = pbuffer.data(idx_dip_gh + 748);

    auto tr_z_xxzz_xzzzz = pbuffer.data(idx_dip_gh + 749);

    auto tr_z_xxzz_yyyyy = pbuffer.data(idx_dip_gh + 750);

    auto tr_z_xxzz_yyyyz = pbuffer.data(idx_dip_gh + 751);

    auto tr_z_xxzz_yyyzz = pbuffer.data(idx_dip_gh + 752);

    auto tr_z_xxzz_yyzzz = pbuffer.data(idx_dip_gh + 753);

    auto tr_z_xxzz_yzzzz = pbuffer.data(idx_dip_gh + 754);

    auto tr_z_xxzz_zzzzz = pbuffer.data(idx_dip_gh + 755);

    auto tr_z_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 757);

    auto tr_z_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 759);

    auto tr_z_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 760);

    auto tr_z_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 762);

    auto tr_z_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 763);

    auto tr_z_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 764);

    auto tr_z_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 766);

    auto tr_z_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 767);

    auto tr_z_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 768);

    auto tr_z_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 769);

    auto tr_z_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 771);

    auto tr_z_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 772);

    auto tr_z_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 773);

    auto tr_z_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 774);

    auto tr_z_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 775);

    auto tr_z_xyyy_zzzzz = pbuffer.data(idx_dip_gh + 776);

    auto tr_z_xyyz_xxxyz = pbuffer.data(idx_dip_gh + 781);

    auto tr_z_xyyz_xxyyz = pbuffer.data(idx_dip_gh + 784);

    auto tr_z_xyyz_xxyzz = pbuffer.data(idx_dip_gh + 785);

    auto tr_z_xyyz_xyyyz = pbuffer.data(idx_dip_gh + 788);

    auto tr_z_xyyz_xyyzz = pbuffer.data(idx_dip_gh + 789);

    auto tr_z_xyyz_xyzzz = pbuffer.data(idx_dip_gh + 790);

    auto tr_z_xyyz_yyyyy = pbuffer.data(idx_dip_gh + 792);

    auto tr_z_xyyz_yyyyz = pbuffer.data(idx_dip_gh + 793);

    auto tr_z_xyyz_yyyzz = pbuffer.data(idx_dip_gh + 794);

    auto tr_z_xyyz_yyzzz = pbuffer.data(idx_dip_gh + 795);

    auto tr_z_xyyz_yzzzz = pbuffer.data(idx_dip_gh + 796);

    auto tr_z_xyyz_zzzzz = pbuffer.data(idx_dip_gh + 797);

    auto tr_z_xyzz_yyyyy = pbuffer.data(idx_dip_gh + 813);

    auto tr_z_xyzz_yyyyz = pbuffer.data(idx_dip_gh + 814);

    auto tr_z_xyzz_yyyzz = pbuffer.data(idx_dip_gh + 815);

    auto tr_z_xyzz_yyzzz = pbuffer.data(idx_dip_gh + 816);

    auto tr_z_xyzz_yzzzz = pbuffer.data(idx_dip_gh + 817);

    auto tr_z_xzzz_xxxxx = pbuffer.data(idx_dip_gh + 819);

    auto tr_z_xzzz_xxxxy = pbuffer.data(idx_dip_gh + 820);

    auto tr_z_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 821);

    auto tr_z_xzzz_xxxyy = pbuffer.data(idx_dip_gh + 822);

    auto tr_z_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 823);

    auto tr_z_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 824);

    auto tr_z_xzzz_xxyyy = pbuffer.data(idx_dip_gh + 825);

    auto tr_z_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 826);

    auto tr_z_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 827);

    auto tr_z_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 828);

    auto tr_z_xzzz_xyyyy = pbuffer.data(idx_dip_gh + 829);

    auto tr_z_xzzz_xyyyz = pbuffer.data(idx_dip_gh + 830);

    auto tr_z_xzzz_xyyzz = pbuffer.data(idx_dip_gh + 831);

    auto tr_z_xzzz_xyzzz = pbuffer.data(idx_dip_gh + 832);

    auto tr_z_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 833);

    auto tr_z_xzzz_yyyyy = pbuffer.data(idx_dip_gh + 834);

    auto tr_z_xzzz_yyyyz = pbuffer.data(idx_dip_gh + 835);

    auto tr_z_xzzz_yyyzz = pbuffer.data(idx_dip_gh + 836);

    auto tr_z_xzzz_yyzzz = pbuffer.data(idx_dip_gh + 837);

    auto tr_z_xzzz_yzzzz = pbuffer.data(idx_dip_gh + 838);

    auto tr_z_xzzz_zzzzz = pbuffer.data(idx_dip_gh + 839);

    auto tr_z_yyyy_xxxxx = pbuffer.data(idx_dip_gh + 840);

    auto tr_z_yyyy_xxxxy = pbuffer.data(idx_dip_gh + 841);

    auto tr_z_yyyy_xxxxz = pbuffer.data(idx_dip_gh + 842);

    auto tr_z_yyyy_xxxyy = pbuffer.data(idx_dip_gh + 843);

    auto tr_z_yyyy_xxxyz = pbuffer.data(idx_dip_gh + 844);

    auto tr_z_yyyy_xxxzz = pbuffer.data(idx_dip_gh + 845);

    auto tr_z_yyyy_xxyyy = pbuffer.data(idx_dip_gh + 846);

    auto tr_z_yyyy_xxyyz = pbuffer.data(idx_dip_gh + 847);

    auto tr_z_yyyy_xxyzz = pbuffer.data(idx_dip_gh + 848);

    auto tr_z_yyyy_xxzzz = pbuffer.data(idx_dip_gh + 849);

    auto tr_z_yyyy_xyyyy = pbuffer.data(idx_dip_gh + 850);

    auto tr_z_yyyy_xyyyz = pbuffer.data(idx_dip_gh + 851);

    auto tr_z_yyyy_xyyzz = pbuffer.data(idx_dip_gh + 852);

    auto tr_z_yyyy_xyzzz = pbuffer.data(idx_dip_gh + 853);

    auto tr_z_yyyy_xzzzz = pbuffer.data(idx_dip_gh + 854);

    auto tr_z_yyyy_yyyyy = pbuffer.data(idx_dip_gh + 855);

    auto tr_z_yyyy_yyyyz = pbuffer.data(idx_dip_gh + 856);

    auto tr_z_yyyy_yyyzz = pbuffer.data(idx_dip_gh + 857);

    auto tr_z_yyyy_yyzzz = pbuffer.data(idx_dip_gh + 858);

    auto tr_z_yyyy_yzzzz = pbuffer.data(idx_dip_gh + 859);

    auto tr_z_yyyy_zzzzz = pbuffer.data(idx_dip_gh + 860);

    auto tr_z_yyyz_xxxxx = pbuffer.data(idx_dip_gh + 861);

    auto tr_z_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 863);

    auto tr_z_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 865);

    auto tr_z_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 866);

    auto tr_z_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 868);

    auto tr_z_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 869);

    auto tr_z_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 870);

    auto tr_z_yyyz_xyyyz = pbuffer.data(idx_dip_gh + 872);

    auto tr_z_yyyz_xyyzz = pbuffer.data(idx_dip_gh + 873);

    auto tr_z_yyyz_xyzzz = pbuffer.data(idx_dip_gh + 874);

    auto tr_z_yyyz_xzzzz = pbuffer.data(idx_dip_gh + 875);

    auto tr_z_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 876);

    auto tr_z_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 877);

    auto tr_z_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 878);

    auto tr_z_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 879);

    auto tr_z_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 880);

    auto tr_z_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 881);

    auto tr_z_yyzz_xxxxx = pbuffer.data(idx_dip_gh + 882);

    auto tr_z_yyzz_xxxxy = pbuffer.data(idx_dip_gh + 883);

    auto tr_z_yyzz_xxxxz = pbuffer.data(idx_dip_gh + 884);

    auto tr_z_yyzz_xxxyy = pbuffer.data(idx_dip_gh + 885);

    auto tr_z_yyzz_xxxyz = pbuffer.data(idx_dip_gh + 886);

    auto tr_z_yyzz_xxxzz = pbuffer.data(idx_dip_gh + 887);

    auto tr_z_yyzz_xxyyy = pbuffer.data(idx_dip_gh + 888);

    auto tr_z_yyzz_xxyyz = pbuffer.data(idx_dip_gh + 889);

    auto tr_z_yyzz_xxyzz = pbuffer.data(idx_dip_gh + 890);

    auto tr_z_yyzz_xxzzz = pbuffer.data(idx_dip_gh + 891);

    auto tr_z_yyzz_xyyyy = pbuffer.data(idx_dip_gh + 892);

    auto tr_z_yyzz_xyyyz = pbuffer.data(idx_dip_gh + 893);

    auto tr_z_yyzz_xyyzz = pbuffer.data(idx_dip_gh + 894);

    auto tr_z_yyzz_xyzzz = pbuffer.data(idx_dip_gh + 895);

    auto tr_z_yyzz_xzzzz = pbuffer.data(idx_dip_gh + 896);

    auto tr_z_yyzz_yyyyy = pbuffer.data(idx_dip_gh + 897);

    auto tr_z_yyzz_yyyyz = pbuffer.data(idx_dip_gh + 898);

    auto tr_z_yyzz_yyyzz = pbuffer.data(idx_dip_gh + 899);

    auto tr_z_yyzz_yyzzz = pbuffer.data(idx_dip_gh + 900);

    auto tr_z_yyzz_yzzzz = pbuffer.data(idx_dip_gh + 901);

    auto tr_z_yyzz_zzzzz = pbuffer.data(idx_dip_gh + 902);

    auto tr_z_yzzz_xxxxx = pbuffer.data(idx_dip_gh + 903);

    auto tr_z_yzzz_xxxxy = pbuffer.data(idx_dip_gh + 904);

    auto tr_z_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 905);

    auto tr_z_yzzz_xxxyy = pbuffer.data(idx_dip_gh + 906);

    auto tr_z_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 907);

    auto tr_z_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 908);

    auto tr_z_yzzz_xxyyy = pbuffer.data(idx_dip_gh + 909);

    auto tr_z_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 910);

    auto tr_z_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 911);

    auto tr_z_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 912);

    auto tr_z_yzzz_xyyyy = pbuffer.data(idx_dip_gh + 913);

    auto tr_z_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 914);

    auto tr_z_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 915);

    auto tr_z_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 916);

    auto tr_z_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 917);

    auto tr_z_yzzz_yyyyy = pbuffer.data(idx_dip_gh + 918);

    auto tr_z_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 919);

    auto tr_z_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 920);

    auto tr_z_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 921);

    auto tr_z_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 922);

    auto tr_z_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 923);

    auto tr_z_zzzz_xxxxx = pbuffer.data(idx_dip_gh + 924);

    auto tr_z_zzzz_xxxxy = pbuffer.data(idx_dip_gh + 925);

    auto tr_z_zzzz_xxxxz = pbuffer.data(idx_dip_gh + 926);

    auto tr_z_zzzz_xxxyy = pbuffer.data(idx_dip_gh + 927);

    auto tr_z_zzzz_xxxyz = pbuffer.data(idx_dip_gh + 928);

    auto tr_z_zzzz_xxxzz = pbuffer.data(idx_dip_gh + 929);

    auto tr_z_zzzz_xxyyy = pbuffer.data(idx_dip_gh + 930);

    auto tr_z_zzzz_xxyyz = pbuffer.data(idx_dip_gh + 931);

    auto tr_z_zzzz_xxyzz = pbuffer.data(idx_dip_gh + 932);

    auto tr_z_zzzz_xxzzz = pbuffer.data(idx_dip_gh + 933);

    auto tr_z_zzzz_xyyyy = pbuffer.data(idx_dip_gh + 934);

    auto tr_z_zzzz_xyyyz = pbuffer.data(idx_dip_gh + 935);

    auto tr_z_zzzz_xyyzz = pbuffer.data(idx_dip_gh + 936);

    auto tr_z_zzzz_xyzzz = pbuffer.data(idx_dip_gh + 937);

    auto tr_z_zzzz_xzzzz = pbuffer.data(idx_dip_gh + 938);

    auto tr_z_zzzz_yyyyy = pbuffer.data(idx_dip_gh + 939);

    auto tr_z_zzzz_yyyyz = pbuffer.data(idx_dip_gh + 940);

    auto tr_z_zzzz_yyyzz = pbuffer.data(idx_dip_gh + 941);

    auto tr_z_zzzz_yyzzz = pbuffer.data(idx_dip_gh + 942);

    auto tr_z_zzzz_yzzzz = pbuffer.data(idx_dip_gh + 943);

    auto tr_z_zzzz_zzzzz = pbuffer.data(idx_dip_gh + 944);

    // Set up components of auxiliary buffer : HG

    auto tr_x_xxxxx_xxxx = pbuffer.data(idx_dip_hg);

    auto tr_x_xxxxx_xxxy = pbuffer.data(idx_dip_hg + 1);

    auto tr_x_xxxxx_xxxz = pbuffer.data(idx_dip_hg + 2);

    auto tr_x_xxxxx_xxyy = pbuffer.data(idx_dip_hg + 3);

    auto tr_x_xxxxx_xxyz = pbuffer.data(idx_dip_hg + 4);

    auto tr_x_xxxxx_xxzz = pbuffer.data(idx_dip_hg + 5);

    auto tr_x_xxxxx_xyyy = pbuffer.data(idx_dip_hg + 6);

    auto tr_x_xxxxx_xyyz = pbuffer.data(idx_dip_hg + 7);

    auto tr_x_xxxxx_xyzz = pbuffer.data(idx_dip_hg + 8);

    auto tr_x_xxxxx_xzzz = pbuffer.data(idx_dip_hg + 9);

    auto tr_x_xxxxx_yyyy = pbuffer.data(idx_dip_hg + 10);

    auto tr_x_xxxxx_yyyz = pbuffer.data(idx_dip_hg + 11);

    auto tr_x_xxxxx_yyzz = pbuffer.data(idx_dip_hg + 12);

    auto tr_x_xxxxx_yzzz = pbuffer.data(idx_dip_hg + 13);

    auto tr_x_xxxxx_zzzz = pbuffer.data(idx_dip_hg + 14);

    auto tr_x_xxxxy_xxxx = pbuffer.data(idx_dip_hg + 15);

    auto tr_x_xxxxy_xxxy = pbuffer.data(idx_dip_hg + 16);

    auto tr_x_xxxxy_xxxz = pbuffer.data(idx_dip_hg + 17);

    auto tr_x_xxxxy_xxyy = pbuffer.data(idx_dip_hg + 18);

    auto tr_x_xxxxy_xxyz = pbuffer.data(idx_dip_hg + 19);

    auto tr_x_xxxxy_xxzz = pbuffer.data(idx_dip_hg + 20);

    auto tr_x_xxxxy_xyyy = pbuffer.data(idx_dip_hg + 21);

    auto tr_x_xxxxy_xyyz = pbuffer.data(idx_dip_hg + 22);

    auto tr_x_xxxxy_xyzz = pbuffer.data(idx_dip_hg + 23);

    auto tr_x_xxxxy_xzzz = pbuffer.data(idx_dip_hg + 24);

    auto tr_x_xxxxz_xxxx = pbuffer.data(idx_dip_hg + 30);

    auto tr_x_xxxxz_xxxy = pbuffer.data(idx_dip_hg + 31);

    auto tr_x_xxxxz_xxxz = pbuffer.data(idx_dip_hg + 32);

    auto tr_x_xxxxz_xxyy = pbuffer.data(idx_dip_hg + 33);

    auto tr_x_xxxxz_xxyz = pbuffer.data(idx_dip_hg + 34);

    auto tr_x_xxxxz_xxzz = pbuffer.data(idx_dip_hg + 35);

    auto tr_x_xxxxz_xyyy = pbuffer.data(idx_dip_hg + 36);

    auto tr_x_xxxxz_xyyz = pbuffer.data(idx_dip_hg + 37);

    auto tr_x_xxxxz_xyzz = pbuffer.data(idx_dip_hg + 38);

    auto tr_x_xxxxz_xzzz = pbuffer.data(idx_dip_hg + 39);

    auto tr_x_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 41);

    auto tr_x_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 42);

    auto tr_x_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 43);

    auto tr_x_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 44);

    auto tr_x_xxxyy_xxxx = pbuffer.data(idx_dip_hg + 45);

    auto tr_x_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 46);

    auto tr_x_xxxyy_xxxz = pbuffer.data(idx_dip_hg + 47);

    auto tr_x_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 48);

    auto tr_x_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 49);

    auto tr_x_xxxyy_xxzz = pbuffer.data(idx_dip_hg + 50);

    auto tr_x_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 51);

    auto tr_x_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 52);

    auto tr_x_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 53);

    auto tr_x_xxxyy_xzzz = pbuffer.data(idx_dip_hg + 54);

    auto tr_x_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 55);

    auto tr_x_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 56);

    auto tr_x_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 57);

    auto tr_x_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 58);

    auto tr_x_xxxzz_xxxx = pbuffer.data(idx_dip_hg + 75);

    auto tr_x_xxxzz_xxxy = pbuffer.data(idx_dip_hg + 76);

    auto tr_x_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 77);

    auto tr_x_xxxzz_xxyy = pbuffer.data(idx_dip_hg + 78);

    auto tr_x_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 79);

    auto tr_x_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 80);

    auto tr_x_xxxzz_xyyy = pbuffer.data(idx_dip_hg + 81);

    auto tr_x_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 82);

    auto tr_x_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 83);

    auto tr_x_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 84);

    auto tr_x_xxxzz_yyyy = pbuffer.data(idx_dip_hg + 85);

    auto tr_x_xxxzz_yyyz = pbuffer.data(idx_dip_hg + 86);

    auto tr_x_xxxzz_yyzz = pbuffer.data(idx_dip_hg + 87);

    auto tr_x_xxxzz_yzzz = pbuffer.data(idx_dip_hg + 88);

    auto tr_x_xxxzz_zzzz = pbuffer.data(idx_dip_hg + 89);

    auto tr_x_xxyyy_xxxx = pbuffer.data(idx_dip_hg + 90);

    auto tr_x_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 91);

    auto tr_x_xxyyy_xxxz = pbuffer.data(idx_dip_hg + 92);

    auto tr_x_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 93);

    auto tr_x_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 94);

    auto tr_x_xxyyy_xxzz = pbuffer.data(idx_dip_hg + 95);

    auto tr_x_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 96);

    auto tr_x_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 97);

    auto tr_x_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 98);

    auto tr_x_xxyyy_xzzz = pbuffer.data(idx_dip_hg + 99);

    auto tr_x_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 100);

    auto tr_x_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 101);

    auto tr_x_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 102);

    auto tr_x_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 103);

    auto tr_x_xxyzz_xxxz = pbuffer.data(idx_dip_hg + 122);

    auto tr_x_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 124);

    auto tr_x_xxyzz_xxzz = pbuffer.data(idx_dip_hg + 125);

    auto tr_x_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 127);

    auto tr_x_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 128);

    auto tr_x_xxyzz_xzzz = pbuffer.data(idx_dip_hg + 129);

    auto tr_x_xxzzz_xxxx = pbuffer.data(idx_dip_hg + 135);

    auto tr_x_xxzzz_xxxy = pbuffer.data(idx_dip_hg + 136);

    auto tr_x_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 137);

    auto tr_x_xxzzz_xxyy = pbuffer.data(idx_dip_hg + 138);

    auto tr_x_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 139);

    auto tr_x_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 140);

    auto tr_x_xxzzz_xyyy = pbuffer.data(idx_dip_hg + 141);

    auto tr_x_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 142);

    auto tr_x_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 143);

    auto tr_x_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 144);

    auto tr_x_xxzzz_yyyy = pbuffer.data(idx_dip_hg + 145);

    auto tr_x_xxzzz_yyyz = pbuffer.data(idx_dip_hg + 146);

    auto tr_x_xxzzz_yyzz = pbuffer.data(idx_dip_hg + 147);

    auto tr_x_xxzzz_yzzz = pbuffer.data(idx_dip_hg + 148);

    auto tr_x_xxzzz_zzzz = pbuffer.data(idx_dip_hg + 149);

    auto tr_x_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 151);

    auto tr_x_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 153);

    auto tr_x_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 154);

    auto tr_x_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 156);

    auto tr_x_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 157);

    auto tr_x_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 158);

    auto tr_x_xzzzz_xxxx = pbuffer.data(idx_dip_hg + 210);

    auto tr_x_xzzzz_xxxy = pbuffer.data(idx_dip_hg + 211);

    auto tr_x_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 212);

    auto tr_x_xzzzz_xxyy = pbuffer.data(idx_dip_hg + 213);

    auto tr_x_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 214);

    auto tr_x_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 215);

    auto tr_x_xzzzz_xyyy = pbuffer.data(idx_dip_hg + 216);

    auto tr_x_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 217);

    auto tr_x_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 218);

    auto tr_x_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 219);

    auto tr_x_yyyyy_xxxx = pbuffer.data(idx_dip_hg + 225);

    auto tr_x_yyyyy_xxxy = pbuffer.data(idx_dip_hg + 226);

    auto tr_x_yyyyy_xxxz = pbuffer.data(idx_dip_hg + 227);

    auto tr_x_yyyyy_xxyy = pbuffer.data(idx_dip_hg + 228);

    auto tr_x_yyyyy_xxyz = pbuffer.data(idx_dip_hg + 229);

    auto tr_x_yyyyy_xxzz = pbuffer.data(idx_dip_hg + 230);

    auto tr_x_yyyyy_xyyy = pbuffer.data(idx_dip_hg + 231);

    auto tr_x_yyyyy_xyyz = pbuffer.data(idx_dip_hg + 232);

    auto tr_x_yyyyy_xyzz = pbuffer.data(idx_dip_hg + 233);

    auto tr_x_yyyyy_xzzz = pbuffer.data(idx_dip_hg + 234);

    auto tr_x_yyyyy_yyyy = pbuffer.data(idx_dip_hg + 235);

    auto tr_x_yyyyy_yyyz = pbuffer.data(idx_dip_hg + 236);

    auto tr_x_yyyyy_yyzz = pbuffer.data(idx_dip_hg + 237);

    auto tr_x_yyyyy_yzzz = pbuffer.data(idx_dip_hg + 238);

    auto tr_x_yyyyy_zzzz = pbuffer.data(idx_dip_hg + 239);

    auto tr_x_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 257);

    auto tr_x_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 259);

    auto tr_x_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 260);

    auto tr_x_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 262);

    auto tr_x_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 263);

    auto tr_x_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 264);

    auto tr_x_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 266);

    auto tr_x_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 267);

    auto tr_x_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 268);

    auto tr_x_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 269);

    auto tr_x_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 272);

    auto tr_x_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 274);

    auto tr_x_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 275);

    auto tr_x_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 277);

    auto tr_x_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 278);

    auto tr_x_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 279);

    auto tr_x_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 281);

    auto tr_x_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 282);

    auto tr_x_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 283);

    auto tr_x_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 284);

    auto tr_x_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 287);

    auto tr_x_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 289);

    auto tr_x_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 290);

    auto tr_x_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 292);

    auto tr_x_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 293);

    auto tr_x_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 294);

    auto tr_x_yzzzz_yyyz = pbuffer.data(idx_dip_hg + 296);

    auto tr_x_yzzzz_yyzz = pbuffer.data(idx_dip_hg + 297);

    auto tr_x_yzzzz_yzzz = pbuffer.data(idx_dip_hg + 298);

    auto tr_x_yzzzz_zzzz = pbuffer.data(idx_dip_hg + 299);

    auto tr_x_zzzzz_xxxx = pbuffer.data(idx_dip_hg + 300);

    auto tr_x_zzzzz_xxxy = pbuffer.data(idx_dip_hg + 301);

    auto tr_x_zzzzz_xxxz = pbuffer.data(idx_dip_hg + 302);

    auto tr_x_zzzzz_xxyy = pbuffer.data(idx_dip_hg + 303);

    auto tr_x_zzzzz_xxyz = pbuffer.data(idx_dip_hg + 304);

    auto tr_x_zzzzz_xxzz = pbuffer.data(idx_dip_hg + 305);

    auto tr_x_zzzzz_xyyy = pbuffer.data(idx_dip_hg + 306);

    auto tr_x_zzzzz_xyyz = pbuffer.data(idx_dip_hg + 307);

    auto tr_x_zzzzz_xyzz = pbuffer.data(idx_dip_hg + 308);

    auto tr_x_zzzzz_xzzz = pbuffer.data(idx_dip_hg + 309);

    auto tr_x_zzzzz_yyyy = pbuffer.data(idx_dip_hg + 310);

    auto tr_x_zzzzz_yyyz = pbuffer.data(idx_dip_hg + 311);

    auto tr_x_zzzzz_yyzz = pbuffer.data(idx_dip_hg + 312);

    auto tr_x_zzzzz_yzzz = pbuffer.data(idx_dip_hg + 313);

    auto tr_x_zzzzz_zzzz = pbuffer.data(idx_dip_hg + 314);

    auto tr_y_xxxxx_xxxx = pbuffer.data(idx_dip_hg + 315);

    auto tr_y_xxxxx_xxxy = pbuffer.data(idx_dip_hg + 316);

    auto tr_y_xxxxx_xxxz = pbuffer.data(idx_dip_hg + 317);

    auto tr_y_xxxxx_xxyy = pbuffer.data(idx_dip_hg + 318);

    auto tr_y_xxxxx_xxyz = pbuffer.data(idx_dip_hg + 319);

    auto tr_y_xxxxx_xxzz = pbuffer.data(idx_dip_hg + 320);

    auto tr_y_xxxxx_xyyy = pbuffer.data(idx_dip_hg + 321);

    auto tr_y_xxxxx_xyyz = pbuffer.data(idx_dip_hg + 322);

    auto tr_y_xxxxx_xyzz = pbuffer.data(idx_dip_hg + 323);

    auto tr_y_xxxxx_xzzz = pbuffer.data(idx_dip_hg + 324);

    auto tr_y_xxxxx_yyyy = pbuffer.data(idx_dip_hg + 325);

    auto tr_y_xxxxx_yyyz = pbuffer.data(idx_dip_hg + 326);

    auto tr_y_xxxxx_yyzz = pbuffer.data(idx_dip_hg + 327);

    auto tr_y_xxxxx_yzzz = pbuffer.data(idx_dip_hg + 328);

    auto tr_y_xxxxx_zzzz = pbuffer.data(idx_dip_hg + 329);

    auto tr_y_xxxxy_xxxy = pbuffer.data(idx_dip_hg + 331);

    auto tr_y_xxxxy_xxyy = pbuffer.data(idx_dip_hg + 333);

    auto tr_y_xxxxy_xxyz = pbuffer.data(idx_dip_hg + 334);

    auto tr_y_xxxxy_xyyy = pbuffer.data(idx_dip_hg + 336);

    auto tr_y_xxxxy_xyyz = pbuffer.data(idx_dip_hg + 337);

    auto tr_y_xxxxy_xyzz = pbuffer.data(idx_dip_hg + 338);

    auto tr_y_xxxxy_yyyy = pbuffer.data(idx_dip_hg + 340);

    auto tr_y_xxxxy_yyyz = pbuffer.data(idx_dip_hg + 341);

    auto tr_y_xxxxy_yyzz = pbuffer.data(idx_dip_hg + 342);

    auto tr_y_xxxxy_yzzz = pbuffer.data(idx_dip_hg + 343);

    auto tr_y_xxxyy_xxxx = pbuffer.data(idx_dip_hg + 360);

    auto tr_y_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 361);

    auto tr_y_xxxyy_xxxz = pbuffer.data(idx_dip_hg + 362);

    auto tr_y_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 363);

    auto tr_y_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 364);

    auto tr_y_xxxyy_xxzz = pbuffer.data(idx_dip_hg + 365);

    auto tr_y_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 366);

    auto tr_y_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 367);

    auto tr_y_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 368);

    auto tr_y_xxxyy_xzzz = pbuffer.data(idx_dip_hg + 369);

    auto tr_y_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 370);

    auto tr_y_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 371);

    auto tr_y_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 372);

    auto tr_y_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 373);

    auto tr_y_xxxyy_zzzz = pbuffer.data(idx_dip_hg + 374);

    auto tr_y_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 392);

    auto tr_y_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 394);

    auto tr_y_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 395);

    auto tr_y_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 397);

    auto tr_y_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 398);

    auto tr_y_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 399);

    auto tr_y_xxxzz_yyyz = pbuffer.data(idx_dip_hg + 401);

    auto tr_y_xxxzz_yyzz = pbuffer.data(idx_dip_hg + 402);

    auto tr_y_xxxzz_yzzz = pbuffer.data(idx_dip_hg + 403);

    auto tr_y_xxxzz_zzzz = pbuffer.data(idx_dip_hg + 404);

    auto tr_y_xxyyy_xxxx = pbuffer.data(idx_dip_hg + 405);

    auto tr_y_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 406);

    auto tr_y_xxyyy_xxxz = pbuffer.data(idx_dip_hg + 407);

    auto tr_y_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 408);

    auto tr_y_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 409);

    auto tr_y_xxyyy_xxzz = pbuffer.data(idx_dip_hg + 410);

    auto tr_y_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 411);

    auto tr_y_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 412);

    auto tr_y_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 413);

    auto tr_y_xxyyy_xzzz = pbuffer.data(idx_dip_hg + 414);

    auto tr_y_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 415);

    auto tr_y_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 416);

    auto tr_y_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 417);

    auto tr_y_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 418);

    auto tr_y_xxyyy_zzzz = pbuffer.data(idx_dip_hg + 419);

    auto tr_y_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 439);

    auto tr_y_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 442);

    auto tr_y_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 443);

    auto tr_y_xxyzz_yyyz = pbuffer.data(idx_dip_hg + 446);

    auto tr_y_xxyzz_yyzz = pbuffer.data(idx_dip_hg + 447);

    auto tr_y_xxyzz_yzzz = pbuffer.data(idx_dip_hg + 448);

    auto tr_y_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 452);

    auto tr_y_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 454);

    auto tr_y_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 455);

    auto tr_y_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 457);

    auto tr_y_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 458);

    auto tr_y_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 459);

    auto tr_y_xxzzz_yyyz = pbuffer.data(idx_dip_hg + 461);

    auto tr_y_xxzzz_yyzz = pbuffer.data(idx_dip_hg + 462);

    auto tr_y_xxzzz_yzzz = pbuffer.data(idx_dip_hg + 463);

    auto tr_y_xxzzz_zzzz = pbuffer.data(idx_dip_hg + 464);

    auto tr_y_xyyyy_xxxx = pbuffer.data(idx_dip_hg + 465);

    auto tr_y_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 466);

    auto tr_y_xyyyy_xxxz = pbuffer.data(idx_dip_hg + 467);

    auto tr_y_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 468);

    auto tr_y_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 469);

    auto tr_y_xyyyy_xxzz = pbuffer.data(idx_dip_hg + 470);

    auto tr_y_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 471);

    auto tr_y_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 472);

    auto tr_y_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 473);

    auto tr_y_xyyyy_xzzz = pbuffer.data(idx_dip_hg + 474);

    auto tr_y_xyyyy_yyyy = pbuffer.data(idx_dip_hg + 475);

    auto tr_y_xyyyy_yyyz = pbuffer.data(idx_dip_hg + 476);

    auto tr_y_xyyyy_yyzz = pbuffer.data(idx_dip_hg + 477);

    auto tr_y_xyyyy_yzzz = pbuffer.data(idx_dip_hg + 478);

    auto tr_y_xyyyy_zzzz = pbuffer.data(idx_dip_hg + 479);

    auto tr_y_xyyzz_xxxz = pbuffer.data(idx_dip_hg + 497);

    auto tr_y_xyyzz_xxyz = pbuffer.data(idx_dip_hg + 499);

    auto tr_y_xyyzz_xxzz = pbuffer.data(idx_dip_hg + 500);

    auto tr_y_xyyzz_xyyz = pbuffer.data(idx_dip_hg + 502);

    auto tr_y_xyyzz_xyzz = pbuffer.data(idx_dip_hg + 503);

    auto tr_y_xyyzz_xzzz = pbuffer.data(idx_dip_hg + 504);

    auto tr_y_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 506);

    auto tr_y_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 507);

    auto tr_y_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 508);

    auto tr_y_xyyzz_zzzz = pbuffer.data(idx_dip_hg + 509);

    auto tr_y_xyzzz_xxyz = pbuffer.data(idx_dip_hg + 514);

    auto tr_y_xyzzz_xyyz = pbuffer.data(idx_dip_hg + 517);

    auto tr_y_xyzzz_xyzz = pbuffer.data(idx_dip_hg + 518);

    auto tr_y_xyzzz_yyyz = pbuffer.data(idx_dip_hg + 521);

    auto tr_y_xyzzz_yyzz = pbuffer.data(idx_dip_hg + 522);

    auto tr_y_xyzzz_yzzz = pbuffer.data(idx_dip_hg + 523);

    auto tr_y_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 527);

    auto tr_y_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 529);

    auto tr_y_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 530);

    auto tr_y_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 532);

    auto tr_y_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 533);

    auto tr_y_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 534);

    auto tr_y_xzzzz_yyyz = pbuffer.data(idx_dip_hg + 536);

    auto tr_y_xzzzz_yyzz = pbuffer.data(idx_dip_hg + 537);

    auto tr_y_xzzzz_yzzz = pbuffer.data(idx_dip_hg + 538);

    auto tr_y_xzzzz_zzzz = pbuffer.data(idx_dip_hg + 539);

    auto tr_y_yyyyy_xxxx = pbuffer.data(idx_dip_hg + 540);

    auto tr_y_yyyyy_xxxy = pbuffer.data(idx_dip_hg + 541);

    auto tr_y_yyyyy_xxxz = pbuffer.data(idx_dip_hg + 542);

    auto tr_y_yyyyy_xxyy = pbuffer.data(idx_dip_hg + 543);

    auto tr_y_yyyyy_xxyz = pbuffer.data(idx_dip_hg + 544);

    auto tr_y_yyyyy_xxzz = pbuffer.data(idx_dip_hg + 545);

    auto tr_y_yyyyy_xyyy = pbuffer.data(idx_dip_hg + 546);

    auto tr_y_yyyyy_xyyz = pbuffer.data(idx_dip_hg + 547);

    auto tr_y_yyyyy_xyzz = pbuffer.data(idx_dip_hg + 548);

    auto tr_y_yyyyy_xzzz = pbuffer.data(idx_dip_hg + 549);

    auto tr_y_yyyyy_yyyy = pbuffer.data(idx_dip_hg + 550);

    auto tr_y_yyyyy_yyyz = pbuffer.data(idx_dip_hg + 551);

    auto tr_y_yyyyy_yyzz = pbuffer.data(idx_dip_hg + 552);

    auto tr_y_yyyyy_yzzz = pbuffer.data(idx_dip_hg + 553);

    auto tr_y_yyyyy_zzzz = pbuffer.data(idx_dip_hg + 554);

    auto tr_y_yyyyz_xxxy = pbuffer.data(idx_dip_hg + 556);

    auto tr_y_yyyyz_xxxz = pbuffer.data(idx_dip_hg + 557);

    auto tr_y_yyyyz_xxyy = pbuffer.data(idx_dip_hg + 558);

    auto tr_y_yyyyz_xxyz = pbuffer.data(idx_dip_hg + 559);

    auto tr_y_yyyyz_xxzz = pbuffer.data(idx_dip_hg + 560);

    auto tr_y_yyyyz_xyyy = pbuffer.data(idx_dip_hg + 561);

    auto tr_y_yyyyz_xyyz = pbuffer.data(idx_dip_hg + 562);

    auto tr_y_yyyyz_xyzz = pbuffer.data(idx_dip_hg + 563);

    auto tr_y_yyyyz_xzzz = pbuffer.data(idx_dip_hg + 564);

    auto tr_y_yyyyz_yyyy = pbuffer.data(idx_dip_hg + 565);

    auto tr_y_yyyyz_yyyz = pbuffer.data(idx_dip_hg + 566);

    auto tr_y_yyyyz_yyzz = pbuffer.data(idx_dip_hg + 567);

    auto tr_y_yyyyz_yzzz = pbuffer.data(idx_dip_hg + 568);

    auto tr_y_yyyyz_zzzz = pbuffer.data(idx_dip_hg + 569);

    auto tr_y_yyyzz_xxxx = pbuffer.data(idx_dip_hg + 570);

    auto tr_y_yyyzz_xxxy = pbuffer.data(idx_dip_hg + 571);

    auto tr_y_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 572);

    auto tr_y_yyyzz_xxyy = pbuffer.data(idx_dip_hg + 573);

    auto tr_y_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 574);

    auto tr_y_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 575);

    auto tr_y_yyyzz_xyyy = pbuffer.data(idx_dip_hg + 576);

    auto tr_y_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 577);

    auto tr_y_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 578);

    auto tr_y_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 579);

    auto tr_y_yyyzz_yyyy = pbuffer.data(idx_dip_hg + 580);

    auto tr_y_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 581);

    auto tr_y_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 582);

    auto tr_y_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 583);

    auto tr_y_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 584);

    auto tr_y_yyzzz_xxxx = pbuffer.data(idx_dip_hg + 585);

    auto tr_y_yyzzz_xxxy = pbuffer.data(idx_dip_hg + 586);

    auto tr_y_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 587);

    auto tr_y_yyzzz_xxyy = pbuffer.data(idx_dip_hg + 588);

    auto tr_y_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 589);

    auto tr_y_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 590);

    auto tr_y_yyzzz_xyyy = pbuffer.data(idx_dip_hg + 591);

    auto tr_y_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 592);

    auto tr_y_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 593);

    auto tr_y_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 594);

    auto tr_y_yyzzz_yyyy = pbuffer.data(idx_dip_hg + 595);

    auto tr_y_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 596);

    auto tr_y_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 597);

    auto tr_y_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 598);

    auto tr_y_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 599);

    auto tr_y_yzzzz_xxxx = pbuffer.data(idx_dip_hg + 600);

    auto tr_y_yzzzz_xxxy = pbuffer.data(idx_dip_hg + 601);

    auto tr_y_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 602);

    auto tr_y_yzzzz_xxyy = pbuffer.data(idx_dip_hg + 603);

    auto tr_y_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 604);

    auto tr_y_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 605);

    auto tr_y_yzzzz_xyyy = pbuffer.data(idx_dip_hg + 606);

    auto tr_y_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 607);

    auto tr_y_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 608);

    auto tr_y_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 609);

    auto tr_y_yzzzz_yyyy = pbuffer.data(idx_dip_hg + 610);

    auto tr_y_yzzzz_yyyz = pbuffer.data(idx_dip_hg + 611);

    auto tr_y_yzzzz_yyzz = pbuffer.data(idx_dip_hg + 612);

    auto tr_y_yzzzz_yzzz = pbuffer.data(idx_dip_hg + 613);

    auto tr_y_yzzzz_zzzz = pbuffer.data(idx_dip_hg + 614);

    auto tr_y_zzzzz_xxxx = pbuffer.data(idx_dip_hg + 615);

    auto tr_y_zzzzz_xxxy = pbuffer.data(idx_dip_hg + 616);

    auto tr_y_zzzzz_xxxz = pbuffer.data(idx_dip_hg + 617);

    auto tr_y_zzzzz_xxyy = pbuffer.data(idx_dip_hg + 618);

    auto tr_y_zzzzz_xxyz = pbuffer.data(idx_dip_hg + 619);

    auto tr_y_zzzzz_xxzz = pbuffer.data(idx_dip_hg + 620);

    auto tr_y_zzzzz_xyyy = pbuffer.data(idx_dip_hg + 621);

    auto tr_y_zzzzz_xyyz = pbuffer.data(idx_dip_hg + 622);

    auto tr_y_zzzzz_xyzz = pbuffer.data(idx_dip_hg + 623);

    auto tr_y_zzzzz_xzzz = pbuffer.data(idx_dip_hg + 624);

    auto tr_y_zzzzz_yyyy = pbuffer.data(idx_dip_hg + 625);

    auto tr_y_zzzzz_yyyz = pbuffer.data(idx_dip_hg + 626);

    auto tr_y_zzzzz_yyzz = pbuffer.data(idx_dip_hg + 627);

    auto tr_y_zzzzz_yzzz = pbuffer.data(idx_dip_hg + 628);

    auto tr_y_zzzzz_zzzz = pbuffer.data(idx_dip_hg + 629);

    auto tr_z_xxxxx_xxxx = pbuffer.data(idx_dip_hg + 630);

    auto tr_z_xxxxx_xxxy = pbuffer.data(idx_dip_hg + 631);

    auto tr_z_xxxxx_xxxz = pbuffer.data(idx_dip_hg + 632);

    auto tr_z_xxxxx_xxyy = pbuffer.data(idx_dip_hg + 633);

    auto tr_z_xxxxx_xxyz = pbuffer.data(idx_dip_hg + 634);

    auto tr_z_xxxxx_xxzz = pbuffer.data(idx_dip_hg + 635);

    auto tr_z_xxxxx_xyyy = pbuffer.data(idx_dip_hg + 636);

    auto tr_z_xxxxx_xyyz = pbuffer.data(idx_dip_hg + 637);

    auto tr_z_xxxxx_xyzz = pbuffer.data(idx_dip_hg + 638);

    auto tr_z_xxxxx_xzzz = pbuffer.data(idx_dip_hg + 639);

    auto tr_z_xxxxx_yyyy = pbuffer.data(idx_dip_hg + 640);

    auto tr_z_xxxxx_yyyz = pbuffer.data(idx_dip_hg + 641);

    auto tr_z_xxxxx_yyzz = pbuffer.data(idx_dip_hg + 642);

    auto tr_z_xxxxx_yzzz = pbuffer.data(idx_dip_hg + 643);

    auto tr_z_xxxxx_zzzz = pbuffer.data(idx_dip_hg + 644);

    auto tr_z_xxxxz_xxxx = pbuffer.data(idx_dip_hg + 660);

    auto tr_z_xxxxz_xxxy = pbuffer.data(idx_dip_hg + 661);

    auto tr_z_xxxxz_xxxz = pbuffer.data(idx_dip_hg + 662);

    auto tr_z_xxxxz_xxyy = pbuffer.data(idx_dip_hg + 663);

    auto tr_z_xxxxz_xxyz = pbuffer.data(idx_dip_hg + 664);

    auto tr_z_xxxxz_xxzz = pbuffer.data(idx_dip_hg + 665);

    auto tr_z_xxxxz_xyyy = pbuffer.data(idx_dip_hg + 666);

    auto tr_z_xxxxz_xyyz = pbuffer.data(idx_dip_hg + 667);

    auto tr_z_xxxxz_xyzz = pbuffer.data(idx_dip_hg + 668);

    auto tr_z_xxxxz_xzzz = pbuffer.data(idx_dip_hg + 669);

    auto tr_z_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 671);

    auto tr_z_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 672);

    auto tr_z_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 673);

    auto tr_z_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 674);

    auto tr_z_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 676);

    auto tr_z_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 678);

    auto tr_z_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 679);

    auto tr_z_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 681);

    auto tr_z_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 682);

    auto tr_z_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 683);

    auto tr_z_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 685);

    auto tr_z_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 686);

    auto tr_z_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 687);

    auto tr_z_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 688);

    auto tr_z_xxxzz_xxxx = pbuffer.data(idx_dip_hg + 705);

    auto tr_z_xxxzz_xxxy = pbuffer.data(idx_dip_hg + 706);

    auto tr_z_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 707);

    auto tr_z_xxxzz_xxyy = pbuffer.data(idx_dip_hg + 708);

    auto tr_z_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 709);

    auto tr_z_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 710);

    auto tr_z_xxxzz_xyyy = pbuffer.data(idx_dip_hg + 711);

    auto tr_z_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 712);

    auto tr_z_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 713);

    auto tr_z_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 714);

    auto tr_z_xxxzz_yyyy = pbuffer.data(idx_dip_hg + 715);

    auto tr_z_xxxzz_yyyz = pbuffer.data(idx_dip_hg + 716);

    auto tr_z_xxxzz_yyzz = pbuffer.data(idx_dip_hg + 717);

    auto tr_z_xxxzz_yzzz = pbuffer.data(idx_dip_hg + 718);

    auto tr_z_xxxzz_zzzz = pbuffer.data(idx_dip_hg + 719);

    auto tr_z_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 721);

    auto tr_z_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 723);

    auto tr_z_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 724);

    auto tr_z_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 726);

    auto tr_z_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 727);

    auto tr_z_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 728);

    auto tr_z_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 730);

    auto tr_z_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 731);

    auto tr_z_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 732);

    auto tr_z_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 733);

    auto tr_z_xxyyz_xxyz = pbuffer.data(idx_dip_hg + 739);

    auto tr_z_xxyyz_xyyz = pbuffer.data(idx_dip_hg + 742);

    auto tr_z_xxyyz_xyzz = pbuffer.data(idx_dip_hg + 743);

    auto tr_z_xxyyz_yyyz = pbuffer.data(idx_dip_hg + 746);

    auto tr_z_xxyyz_yyzz = pbuffer.data(idx_dip_hg + 747);

    auto tr_z_xxyyz_yzzz = pbuffer.data(idx_dip_hg + 748);

    auto tr_z_xxzzz_xxxx = pbuffer.data(idx_dip_hg + 765);

    auto tr_z_xxzzz_xxxy = pbuffer.data(idx_dip_hg + 766);

    auto tr_z_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 767);

    auto tr_z_xxzzz_xxyy = pbuffer.data(idx_dip_hg + 768);

    auto tr_z_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 769);

    auto tr_z_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 770);

    auto tr_z_xxzzz_xyyy = pbuffer.data(idx_dip_hg + 771);

    auto tr_z_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 772);

    auto tr_z_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 773);

    auto tr_z_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 774);

    auto tr_z_xxzzz_yyyy = pbuffer.data(idx_dip_hg + 775);

    auto tr_z_xxzzz_yyyz = pbuffer.data(idx_dip_hg + 776);

    auto tr_z_xxzzz_yyzz = pbuffer.data(idx_dip_hg + 777);

    auto tr_z_xxzzz_yzzz = pbuffer.data(idx_dip_hg + 778);

    auto tr_z_xxzzz_zzzz = pbuffer.data(idx_dip_hg + 779);

    auto tr_z_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 781);

    auto tr_z_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 783);

    auto tr_z_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 784);

    auto tr_z_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 786);

    auto tr_z_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 787);

    auto tr_z_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 788);

    auto tr_z_xyyyy_yyyy = pbuffer.data(idx_dip_hg + 790);

    auto tr_z_xyyyy_yyyz = pbuffer.data(idx_dip_hg + 791);

    auto tr_z_xyyyy_yyzz = pbuffer.data(idx_dip_hg + 792);

    auto tr_z_xyyyy_yzzz = pbuffer.data(idx_dip_hg + 793);

    auto tr_z_xyyyz_xxyz = pbuffer.data(idx_dip_hg + 799);

    auto tr_z_xyyyz_xyyz = pbuffer.data(idx_dip_hg + 802);

    auto tr_z_xyyyz_xyzz = pbuffer.data(idx_dip_hg + 803);

    auto tr_z_xyyyz_yyyz = pbuffer.data(idx_dip_hg + 806);

    auto tr_z_xyyyz_yyzz = pbuffer.data(idx_dip_hg + 807);

    auto tr_z_xyyyz_yzzz = pbuffer.data(idx_dip_hg + 808);

    auto tr_z_xyyzz_xxxy = pbuffer.data(idx_dip_hg + 811);

    auto tr_z_xyyzz_xxyy = pbuffer.data(idx_dip_hg + 813);

    auto tr_z_xyyzz_xxyz = pbuffer.data(idx_dip_hg + 814);

    auto tr_z_xyyzz_xyyy = pbuffer.data(idx_dip_hg + 816);

    auto tr_z_xyyzz_xyyz = pbuffer.data(idx_dip_hg + 817);

    auto tr_z_xyyzz_xyzz = pbuffer.data(idx_dip_hg + 818);

    auto tr_z_xyyzz_yyyy = pbuffer.data(idx_dip_hg + 820);

    auto tr_z_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 821);

    auto tr_z_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 822);

    auto tr_z_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 823);

    auto tr_z_xzzzz_xxxx = pbuffer.data(idx_dip_hg + 840);

    auto tr_z_xzzzz_xxxy = pbuffer.data(idx_dip_hg + 841);

    auto tr_z_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 842);

    auto tr_z_xzzzz_xxyy = pbuffer.data(idx_dip_hg + 843);

    auto tr_z_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 844);

    auto tr_z_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 845);

    auto tr_z_xzzzz_xyyy = pbuffer.data(idx_dip_hg + 846);

    auto tr_z_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 847);

    auto tr_z_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 848);

    auto tr_z_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 849);

    auto tr_z_xzzzz_yyyy = pbuffer.data(idx_dip_hg + 850);

    auto tr_z_xzzzz_yyyz = pbuffer.data(idx_dip_hg + 851);

    auto tr_z_xzzzz_yyzz = pbuffer.data(idx_dip_hg + 852);

    auto tr_z_xzzzz_yzzz = pbuffer.data(idx_dip_hg + 853);

    auto tr_z_xzzzz_zzzz = pbuffer.data(idx_dip_hg + 854);

    auto tr_z_yyyyy_xxxx = pbuffer.data(idx_dip_hg + 855);

    auto tr_z_yyyyy_xxxy = pbuffer.data(idx_dip_hg + 856);

    auto tr_z_yyyyy_xxxz = pbuffer.data(idx_dip_hg + 857);

    auto tr_z_yyyyy_xxyy = pbuffer.data(idx_dip_hg + 858);

    auto tr_z_yyyyy_xxyz = pbuffer.data(idx_dip_hg + 859);

    auto tr_z_yyyyy_xxzz = pbuffer.data(idx_dip_hg + 860);

    auto tr_z_yyyyy_xyyy = pbuffer.data(idx_dip_hg + 861);

    auto tr_z_yyyyy_xyyz = pbuffer.data(idx_dip_hg + 862);

    auto tr_z_yyyyy_xyzz = pbuffer.data(idx_dip_hg + 863);

    auto tr_z_yyyyy_xzzz = pbuffer.data(idx_dip_hg + 864);

    auto tr_z_yyyyy_yyyy = pbuffer.data(idx_dip_hg + 865);

    auto tr_z_yyyyy_yyyz = pbuffer.data(idx_dip_hg + 866);

    auto tr_z_yyyyy_yyzz = pbuffer.data(idx_dip_hg + 867);

    auto tr_z_yyyyy_yzzz = pbuffer.data(idx_dip_hg + 868);

    auto tr_z_yyyyy_zzzz = pbuffer.data(idx_dip_hg + 869);

    auto tr_z_yyyyz_xxxx = pbuffer.data(idx_dip_hg + 870);

    auto tr_z_yyyyz_xxxy = pbuffer.data(idx_dip_hg + 871);

    auto tr_z_yyyyz_xxxz = pbuffer.data(idx_dip_hg + 872);

    auto tr_z_yyyyz_xxyy = pbuffer.data(idx_dip_hg + 873);

    auto tr_z_yyyyz_xxyz = pbuffer.data(idx_dip_hg + 874);

    auto tr_z_yyyyz_xxzz = pbuffer.data(idx_dip_hg + 875);

    auto tr_z_yyyyz_xyyy = pbuffer.data(idx_dip_hg + 876);

    auto tr_z_yyyyz_xyyz = pbuffer.data(idx_dip_hg + 877);

    auto tr_z_yyyyz_xyzz = pbuffer.data(idx_dip_hg + 878);

    auto tr_z_yyyyz_xzzz = pbuffer.data(idx_dip_hg + 879);

    auto tr_z_yyyyz_yyyy = pbuffer.data(idx_dip_hg + 880);

    auto tr_z_yyyyz_yyyz = pbuffer.data(idx_dip_hg + 881);

    auto tr_z_yyyyz_yyzz = pbuffer.data(idx_dip_hg + 882);

    auto tr_z_yyyyz_yzzz = pbuffer.data(idx_dip_hg + 883);

    auto tr_z_yyyyz_zzzz = pbuffer.data(idx_dip_hg + 884);

    auto tr_z_yyyzz_xxxx = pbuffer.data(idx_dip_hg + 885);

    auto tr_z_yyyzz_xxxy = pbuffer.data(idx_dip_hg + 886);

    auto tr_z_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 887);

    auto tr_z_yyyzz_xxyy = pbuffer.data(idx_dip_hg + 888);

    auto tr_z_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 889);

    auto tr_z_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 890);

    auto tr_z_yyyzz_xyyy = pbuffer.data(idx_dip_hg + 891);

    auto tr_z_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 892);

    auto tr_z_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 893);

    auto tr_z_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 894);

    auto tr_z_yyyzz_yyyy = pbuffer.data(idx_dip_hg + 895);

    auto tr_z_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 896);

    auto tr_z_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 897);

    auto tr_z_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 898);

    auto tr_z_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 899);

    auto tr_z_yyzzz_xxxx = pbuffer.data(idx_dip_hg + 900);

    auto tr_z_yyzzz_xxxy = pbuffer.data(idx_dip_hg + 901);

    auto tr_z_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 902);

    auto tr_z_yyzzz_xxyy = pbuffer.data(idx_dip_hg + 903);

    auto tr_z_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 904);

    auto tr_z_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 905);

    auto tr_z_yyzzz_xyyy = pbuffer.data(idx_dip_hg + 906);

    auto tr_z_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 907);

    auto tr_z_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 908);

    auto tr_z_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 909);

    auto tr_z_yyzzz_yyyy = pbuffer.data(idx_dip_hg + 910);

    auto tr_z_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 911);

    auto tr_z_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 912);

    auto tr_z_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 913);

    auto tr_z_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 914);

    auto tr_z_yzzzz_xxxx = pbuffer.data(idx_dip_hg + 915);

    auto tr_z_yzzzz_xxxy = pbuffer.data(idx_dip_hg + 916);

    auto tr_z_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 917);

    auto tr_z_yzzzz_xxyy = pbuffer.data(idx_dip_hg + 918);

    auto tr_z_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 919);

    auto tr_z_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 920);

    auto tr_z_yzzzz_xyyy = pbuffer.data(idx_dip_hg + 921);

    auto tr_z_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 922);

    auto tr_z_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 923);

    auto tr_z_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 924);

    auto tr_z_yzzzz_yyyy = pbuffer.data(idx_dip_hg + 925);

    auto tr_z_yzzzz_yyyz = pbuffer.data(idx_dip_hg + 926);

    auto tr_z_yzzzz_yyzz = pbuffer.data(idx_dip_hg + 927);

    auto tr_z_yzzzz_yzzz = pbuffer.data(idx_dip_hg + 928);

    auto tr_z_yzzzz_zzzz = pbuffer.data(idx_dip_hg + 929);

    auto tr_z_zzzzz_xxxx = pbuffer.data(idx_dip_hg + 930);

    auto tr_z_zzzzz_xxxy = pbuffer.data(idx_dip_hg + 931);

    auto tr_z_zzzzz_xxxz = pbuffer.data(idx_dip_hg + 932);

    auto tr_z_zzzzz_xxyy = pbuffer.data(idx_dip_hg + 933);

    auto tr_z_zzzzz_xxyz = pbuffer.data(idx_dip_hg + 934);

    auto tr_z_zzzzz_xxzz = pbuffer.data(idx_dip_hg + 935);

    auto tr_z_zzzzz_xyyy = pbuffer.data(idx_dip_hg + 936);

    auto tr_z_zzzzz_xyyz = pbuffer.data(idx_dip_hg + 937);

    auto tr_z_zzzzz_xyzz = pbuffer.data(idx_dip_hg + 938);

    auto tr_z_zzzzz_xzzz = pbuffer.data(idx_dip_hg + 939);

    auto tr_z_zzzzz_yyyy = pbuffer.data(idx_dip_hg + 940);

    auto tr_z_zzzzz_yyyz = pbuffer.data(idx_dip_hg + 941);

    auto tr_z_zzzzz_yyzz = pbuffer.data(idx_dip_hg + 942);

    auto tr_z_zzzzz_yzzz = pbuffer.data(idx_dip_hg + 943);

    auto tr_z_zzzzz_zzzz = pbuffer.data(idx_dip_hg + 944);

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

    auto ts_xxxxz_xxxzz = pbuffer.data(idx_ovl_hh + 47);

    auto ts_xxxxz_xxzzz = pbuffer.data(idx_ovl_hh + 51);

    auto ts_xxxxz_xzzzz = pbuffer.data(idx_ovl_hh + 56);

    auto ts_xxxyy_xxxxy = pbuffer.data(idx_ovl_hh + 64);

    auto ts_xxxyy_xxxyy = pbuffer.data(idx_ovl_hh + 66);

    auto ts_xxxyy_xxyyy = pbuffer.data(idx_ovl_hh + 69);

    auto ts_xxxyy_xyyyy = pbuffer.data(idx_ovl_hh + 73);

    auto ts_xxxyy_yyyyy = pbuffer.data(idx_ovl_hh + 78);

    auto ts_xxxyy_yyyyz = pbuffer.data(idx_ovl_hh + 79);

    auto ts_xxxyy_yyyzz = pbuffer.data(idx_ovl_hh + 80);

    auto ts_xxxyy_yyzzz = pbuffer.data(idx_ovl_hh + 81);

    auto ts_xxxyy_yzzzz = pbuffer.data(idx_ovl_hh + 82);

    auto ts_xxxzz_xxxxx = pbuffer.data(idx_ovl_hh + 105);

    auto ts_xxxzz_xxxxz = pbuffer.data(idx_ovl_hh + 107);

    auto ts_xxxzz_xxxzz = pbuffer.data(idx_ovl_hh + 110);

    auto ts_xxxzz_xxzzz = pbuffer.data(idx_ovl_hh + 114);

    auto ts_xxxzz_xzzzz = pbuffer.data(idx_ovl_hh + 119);

    auto ts_xxxzz_yyyyz = pbuffer.data(idx_ovl_hh + 121);

    auto ts_xxxzz_yyyzz = pbuffer.data(idx_ovl_hh + 122);

    auto ts_xxxzz_yyzzz = pbuffer.data(idx_ovl_hh + 123);

    auto ts_xxxzz_yzzzz = pbuffer.data(idx_ovl_hh + 124);

    auto ts_xxxzz_zzzzz = pbuffer.data(idx_ovl_hh + 125);

    auto ts_xxyyy_xxxxy = pbuffer.data(idx_ovl_hh + 127);

    auto ts_xxyyy_xxxyy = pbuffer.data(idx_ovl_hh + 129);

    auto ts_xxyyy_xxyyy = pbuffer.data(idx_ovl_hh + 132);

    auto ts_xxyyy_xyyyy = pbuffer.data(idx_ovl_hh + 136);

    auto ts_xxyyy_yyyyy = pbuffer.data(idx_ovl_hh + 141);

    auto ts_xxyyy_yyyyz = pbuffer.data(idx_ovl_hh + 142);

    auto ts_xxyyy_yyyzz = pbuffer.data(idx_ovl_hh + 143);

    auto ts_xxyyy_yyzzz = pbuffer.data(idx_ovl_hh + 144);

    auto ts_xxyyy_yzzzz = pbuffer.data(idx_ovl_hh + 145);

    auto ts_xxzzz_xxxxx = pbuffer.data(idx_ovl_hh + 189);

    auto ts_xxzzz_xxxxz = pbuffer.data(idx_ovl_hh + 191);

    auto ts_xxzzz_xxxzz = pbuffer.data(idx_ovl_hh + 194);

    auto ts_xxzzz_xxzzz = pbuffer.data(idx_ovl_hh + 198);

    auto ts_xxzzz_xzzzz = pbuffer.data(idx_ovl_hh + 203);

    auto ts_xxzzz_yyyyz = pbuffer.data(idx_ovl_hh + 205);

    auto ts_xxzzz_yyyzz = pbuffer.data(idx_ovl_hh + 206);

    auto ts_xxzzz_yyzzz = pbuffer.data(idx_ovl_hh + 207);

    auto ts_xxzzz_yzzzz = pbuffer.data(idx_ovl_hh + 208);

    auto ts_xxzzz_zzzzz = pbuffer.data(idx_ovl_hh + 209);

    auto ts_xyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 225);

    auto ts_xyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 226);

    auto ts_xyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 227);

    auto ts_xyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 228);

    auto ts_xyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 229);

    auto ts_xyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 268);

    auto ts_xyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 269);

    auto ts_xyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 270);

    auto ts_xyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 271);

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

    auto ts_yyyyz_yyyyz = pbuffer.data(idx_ovl_hh + 352);

    auto ts_yyyyz_yyyzz = pbuffer.data(idx_ovl_hh + 353);

    auto ts_yyyyz_yyzzz = pbuffer.data(idx_ovl_hh + 354);

    auto ts_yyyyz_yzzzz = pbuffer.data(idx_ovl_hh + 355);

    auto ts_yyyyz_zzzzz = pbuffer.data(idx_ovl_hh + 356);

    auto ts_yyyzz_xxxxz = pbuffer.data(idx_ovl_hh + 359);

    auto ts_yyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 361);

    auto ts_yyyzz_xxxzz = pbuffer.data(idx_ovl_hh + 362);

    auto ts_yyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 364);

    auto ts_yyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 365);

    auto ts_yyyzz_xxzzz = pbuffer.data(idx_ovl_hh + 366);

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

    auto ts_yyzzz_xxxxz = pbuffer.data(idx_ovl_hh + 380);

    auto ts_yyzzz_xxxyz = pbuffer.data(idx_ovl_hh + 382);

    auto ts_yyzzz_xxxzz = pbuffer.data(idx_ovl_hh + 383);

    auto ts_yyzzz_xxyyz = pbuffer.data(idx_ovl_hh + 385);

    auto ts_yyzzz_xxyzz = pbuffer.data(idx_ovl_hh + 386);

    auto ts_yyzzz_xxzzz = pbuffer.data(idx_ovl_hh + 387);

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

    auto ts_yzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 401);

    auto ts_yzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 404);

    auto ts_yzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 408);

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

    auto tr_x_xxxxy_yyyyy = pbuffer.data(idx_dip_hh + 36);

    auto tr_x_xxxxy_zzzzz = pbuffer.data(idx_dip_hh + 41);

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

    auto tr_x_xxxxz_yyyyy = pbuffer.data(idx_dip_hh + 57);

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

    auto tr_x_xxxyy_zzzzz = pbuffer.data(idx_dip_hh + 83);

    auto tr_x_xxxyz_xxxxz = pbuffer.data(idx_dip_hh + 86);

    auto tr_x_xxxyz_xxxzz = pbuffer.data(idx_dip_hh + 89);

    auto tr_x_xxxyz_xxzzz = pbuffer.data(idx_dip_hh + 93);

    auto tr_x_xxxyz_xzzzz = pbuffer.data(idx_dip_hh + 98);

    auto tr_x_xxxyz_zzzzz = pbuffer.data(idx_dip_hh + 104);

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

    auto tr_x_xxyyy_zzzzz = pbuffer.data(idx_dip_hh + 146);

    auto tr_x_xxyyz_xxxxy = pbuffer.data(idx_dip_hh + 148);

    auto tr_x_xxyyz_xxxxz = pbuffer.data(idx_dip_hh + 149);

    auto tr_x_xxyyz_xxxyy = pbuffer.data(idx_dip_hh + 150);

    auto tr_x_xxyyz_xxxzz = pbuffer.data(idx_dip_hh + 152);

    auto tr_x_xxyyz_xxyyy = pbuffer.data(idx_dip_hh + 153);

    auto tr_x_xxyyz_xxzzz = pbuffer.data(idx_dip_hh + 156);

    auto tr_x_xxyyz_xyyyy = pbuffer.data(idx_dip_hh + 157);

    auto tr_x_xxyyz_xzzzz = pbuffer.data(idx_dip_hh + 161);

    auto tr_x_xxyyz_yyyyy = pbuffer.data(idx_dip_hh + 162);

    auto tr_x_xxyyz_zzzzz = pbuffer.data(idx_dip_hh + 167);

    auto tr_x_xxyzz_xxxxx = pbuffer.data(idx_dip_hh + 168);

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

    auto tr_x_xxyzz_zzzzz = pbuffer.data(idx_dip_hh + 188);

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

    auto tr_x_xyyyy_xxxxx = pbuffer.data(idx_dip_hh + 210);

    auto tr_x_xyyyy_xxxxy = pbuffer.data(idx_dip_hh + 211);

    auto tr_x_xyyyy_xxxxz = pbuffer.data(idx_dip_hh + 212);

    auto tr_x_xyyyy_xxxyy = pbuffer.data(idx_dip_hh + 213);

    auto tr_x_xyyyy_xxxyz = pbuffer.data(idx_dip_hh + 214);

    auto tr_x_xyyyy_xxxzz = pbuffer.data(idx_dip_hh + 215);

    auto tr_x_xyyyy_xxyyy = pbuffer.data(idx_dip_hh + 216);

    auto tr_x_xyyyy_xxyyz = pbuffer.data(idx_dip_hh + 217);

    auto tr_x_xyyyy_xxyzz = pbuffer.data(idx_dip_hh + 218);

    auto tr_x_xyyyy_xxzzz = pbuffer.data(idx_dip_hh + 219);

    auto tr_x_xyyyy_xyyyy = pbuffer.data(idx_dip_hh + 220);

    auto tr_x_xyyyy_xyyyz = pbuffer.data(idx_dip_hh + 221);

    auto tr_x_xyyyy_xyyzz = pbuffer.data(idx_dip_hh + 222);

    auto tr_x_xyyyy_xyzzz = pbuffer.data(idx_dip_hh + 223);

    auto tr_x_xyyyy_xzzzz = pbuffer.data(idx_dip_hh + 224);

    auto tr_x_xyyyy_yyyyy = pbuffer.data(idx_dip_hh + 225);

    auto tr_x_xyyyy_yyyyz = pbuffer.data(idx_dip_hh + 226);

    auto tr_x_xyyyy_yyyzz = pbuffer.data(idx_dip_hh + 227);

    auto tr_x_xyyyy_yyzzz = pbuffer.data(idx_dip_hh + 228);

    auto tr_x_xyyyy_yzzzz = pbuffer.data(idx_dip_hh + 229);

    auto tr_x_xyyyz_xxxxy = pbuffer.data(idx_dip_hh + 232);

    auto tr_x_xyyyz_xxxxz = pbuffer.data(idx_dip_hh + 233);

    auto tr_x_xyyyz_xxxyy = pbuffer.data(idx_dip_hh + 234);

    auto tr_x_xyyyz_xxxzz = pbuffer.data(idx_dip_hh + 236);

    auto tr_x_xyyyz_xxyyy = pbuffer.data(idx_dip_hh + 237);

    auto tr_x_xyyyz_xxzzz = pbuffer.data(idx_dip_hh + 240);

    auto tr_x_xyyyz_xyyyy = pbuffer.data(idx_dip_hh + 241);

    auto tr_x_xyyyz_xzzzz = pbuffer.data(idx_dip_hh + 245);

    auto tr_x_xyyzz_xxxxx = pbuffer.data(idx_dip_hh + 252);

    auto tr_x_xyyzz_xxxxy = pbuffer.data(idx_dip_hh + 253);

    auto tr_x_xyyzz_xxxxz = pbuffer.data(idx_dip_hh + 254);

    auto tr_x_xyyzz_xxxyy = pbuffer.data(idx_dip_hh + 255);

    auto tr_x_xyyzz_xxxzz = pbuffer.data(idx_dip_hh + 257);

    auto tr_x_xyyzz_xxyyy = pbuffer.data(idx_dip_hh + 258);

    auto tr_x_xyyzz_xxzzz = pbuffer.data(idx_dip_hh + 261);

    auto tr_x_xyyzz_xyyyy = pbuffer.data(idx_dip_hh + 262);

    auto tr_x_xyyzz_xzzzz = pbuffer.data(idx_dip_hh + 266);

    auto tr_x_xyyzz_yyyyz = pbuffer.data(idx_dip_hh + 268);

    auto tr_x_xyyzz_yyyzz = pbuffer.data(idx_dip_hh + 269);

    auto tr_x_xyyzz_yyzzz = pbuffer.data(idx_dip_hh + 270);

    auto tr_x_xyyzz_yzzzz = pbuffer.data(idx_dip_hh + 271);

    auto tr_x_xyzzz_xxxxx = pbuffer.data(idx_dip_hh + 273);

    auto tr_x_xyzzz_xxxxz = pbuffer.data(idx_dip_hh + 275);

    auto tr_x_xyzzz_xxxzz = pbuffer.data(idx_dip_hh + 278);

    auto tr_x_xyzzz_xxzzz = pbuffer.data(idx_dip_hh + 282);

    auto tr_x_xyzzz_xzzzz = pbuffer.data(idx_dip_hh + 287);

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

    auto tr_x_xzzzz_yyyyz = pbuffer.data(idx_dip_hh + 310);

    auto tr_x_xzzzz_yyyzz = pbuffer.data(idx_dip_hh + 311);

    auto tr_x_xzzzz_yyzzz = pbuffer.data(idx_dip_hh + 312);

    auto tr_x_xzzzz_yzzzz = pbuffer.data(idx_dip_hh + 313);

    auto tr_x_xzzzz_zzzzz = pbuffer.data(idx_dip_hh + 314);

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

    auto tr_x_yyyyz_xxxxy = pbuffer.data(idx_dip_hh + 337);

    auto tr_x_yyyyz_xxxxz = pbuffer.data(idx_dip_hh + 338);

    auto tr_x_yyyyz_xxxyy = pbuffer.data(idx_dip_hh + 339);

    auto tr_x_yyyyz_xxxzz = pbuffer.data(idx_dip_hh + 341);

    auto tr_x_yyyyz_xxyyy = pbuffer.data(idx_dip_hh + 342);

    auto tr_x_yyyyz_xxzzz = pbuffer.data(idx_dip_hh + 345);

    auto tr_x_yyyyz_xyyyy = pbuffer.data(idx_dip_hh + 346);

    auto tr_x_yyyyz_xzzzz = pbuffer.data(idx_dip_hh + 350);

    auto tr_x_yyyyz_yyyyy = pbuffer.data(idx_dip_hh + 351);

    auto tr_x_yyyyz_yyyyz = pbuffer.data(idx_dip_hh + 352);

    auto tr_x_yyyyz_yyyzz = pbuffer.data(idx_dip_hh + 353);

    auto tr_x_yyyyz_yyzzz = pbuffer.data(idx_dip_hh + 354);

    auto tr_x_yyyyz_yzzzz = pbuffer.data(idx_dip_hh + 355);

    auto tr_x_yyyyz_zzzzz = pbuffer.data(idx_dip_hh + 356);

    auto tr_x_yyyzz_xxxxx = pbuffer.data(idx_dip_hh + 357);

    auto tr_x_yyyzz_xxxxy = pbuffer.data(idx_dip_hh + 358);

    auto tr_x_yyyzz_xxxxz = pbuffer.data(idx_dip_hh + 359);

    auto tr_x_yyyzz_xxxyy = pbuffer.data(idx_dip_hh + 360);

    auto tr_x_yyyzz_xxxyz = pbuffer.data(idx_dip_hh + 361);

    auto tr_x_yyyzz_xxxzz = pbuffer.data(idx_dip_hh + 362);

    auto tr_x_yyyzz_xxyyy = pbuffer.data(idx_dip_hh + 363);

    auto tr_x_yyyzz_xxyyz = pbuffer.data(idx_dip_hh + 364);

    auto tr_x_yyyzz_xxyzz = pbuffer.data(idx_dip_hh + 365);

    auto tr_x_yyyzz_xxzzz = pbuffer.data(idx_dip_hh + 366);

    auto tr_x_yyyzz_xyyyy = pbuffer.data(idx_dip_hh + 367);

    auto tr_x_yyyzz_xyyyz = pbuffer.data(idx_dip_hh + 368);

    auto tr_x_yyyzz_xyyzz = pbuffer.data(idx_dip_hh + 369);

    auto tr_x_yyyzz_xyzzz = pbuffer.data(idx_dip_hh + 370);

    auto tr_x_yyyzz_xzzzz = pbuffer.data(idx_dip_hh + 371);

    auto tr_x_yyyzz_yyyyy = pbuffer.data(idx_dip_hh + 372);

    auto tr_x_yyyzz_yyyyz = pbuffer.data(idx_dip_hh + 373);

    auto tr_x_yyyzz_yyyzz = pbuffer.data(idx_dip_hh + 374);

    auto tr_x_yyyzz_yyzzz = pbuffer.data(idx_dip_hh + 375);

    auto tr_x_yyyzz_yzzzz = pbuffer.data(idx_dip_hh + 376);

    auto tr_x_yyyzz_zzzzz = pbuffer.data(idx_dip_hh + 377);

    auto tr_x_yyzzz_xxxxx = pbuffer.data(idx_dip_hh + 378);

    auto tr_x_yyzzz_xxxxy = pbuffer.data(idx_dip_hh + 379);

    auto tr_x_yyzzz_xxxxz = pbuffer.data(idx_dip_hh + 380);

    auto tr_x_yyzzz_xxxyy = pbuffer.data(idx_dip_hh + 381);

    auto tr_x_yyzzz_xxxyz = pbuffer.data(idx_dip_hh + 382);

    auto tr_x_yyzzz_xxxzz = pbuffer.data(idx_dip_hh + 383);

    auto tr_x_yyzzz_xxyyy = pbuffer.data(idx_dip_hh + 384);

    auto tr_x_yyzzz_xxyyz = pbuffer.data(idx_dip_hh + 385);

    auto tr_x_yyzzz_xxyzz = pbuffer.data(idx_dip_hh + 386);

    auto tr_x_yyzzz_xxzzz = pbuffer.data(idx_dip_hh + 387);

    auto tr_x_yyzzz_xyyyy = pbuffer.data(idx_dip_hh + 388);

    auto tr_x_yyzzz_xyyyz = pbuffer.data(idx_dip_hh + 389);

    auto tr_x_yyzzz_xyyzz = pbuffer.data(idx_dip_hh + 390);

    auto tr_x_yyzzz_xyzzz = pbuffer.data(idx_dip_hh + 391);

    auto tr_x_yyzzz_xzzzz = pbuffer.data(idx_dip_hh + 392);

    auto tr_x_yyzzz_yyyyy = pbuffer.data(idx_dip_hh + 393);

    auto tr_x_yyzzz_yyyyz = pbuffer.data(idx_dip_hh + 394);

    auto tr_x_yyzzz_yyyzz = pbuffer.data(idx_dip_hh + 395);

    auto tr_x_yyzzz_yyzzz = pbuffer.data(idx_dip_hh + 396);

    auto tr_x_yyzzz_yzzzz = pbuffer.data(idx_dip_hh + 397);

    auto tr_x_yyzzz_zzzzz = pbuffer.data(idx_dip_hh + 398);

    auto tr_x_yzzzz_xxxxx = pbuffer.data(idx_dip_hh + 399);

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

    auto tr_x_yzzzz_yyyyy = pbuffer.data(idx_dip_hh + 414);

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

    auto tr_y_xxxxy_xxxxx = pbuffer.data(idx_dip_hh + 462);

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

    auto tr_y_xxxxy_zzzzz = pbuffer.data(idx_dip_hh + 482);

    auto tr_y_xxxxz_xxxxx = pbuffer.data(idx_dip_hh + 483);

    auto tr_y_xxxxz_xxxxy = pbuffer.data(idx_dip_hh + 484);

    auto tr_y_xxxxz_xxxxz = pbuffer.data(idx_dip_hh + 485);

    auto tr_y_xxxxz_xxxyy = pbuffer.data(idx_dip_hh + 486);

    auto tr_y_xxxxz_xxxzz = pbuffer.data(idx_dip_hh + 488);

    auto tr_y_xxxxz_xxyyy = pbuffer.data(idx_dip_hh + 489);

    auto tr_y_xxxxz_xxzzz = pbuffer.data(idx_dip_hh + 492);

    auto tr_y_xxxxz_xyyyy = pbuffer.data(idx_dip_hh + 493);

    auto tr_y_xxxxz_xzzzz = pbuffer.data(idx_dip_hh + 497);

    auto tr_y_xxxxz_yyyyz = pbuffer.data(idx_dip_hh + 499);

    auto tr_y_xxxxz_yyyzz = pbuffer.data(idx_dip_hh + 500);

    auto tr_y_xxxxz_yyzzz = pbuffer.data(idx_dip_hh + 501);

    auto tr_y_xxxxz_yzzzz = pbuffer.data(idx_dip_hh + 502);

    auto tr_y_xxxxz_zzzzz = pbuffer.data(idx_dip_hh + 503);

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

    auto tr_y_xxxyz_xxxxy = pbuffer.data(idx_dip_hh + 526);

    auto tr_y_xxxyz_xxxyy = pbuffer.data(idx_dip_hh + 528);

    auto tr_y_xxxyz_xxyyy = pbuffer.data(idx_dip_hh + 531);

    auto tr_y_xxxyz_xyyyy = pbuffer.data(idx_dip_hh + 535);

    auto tr_y_xxxyz_yyyyz = pbuffer.data(idx_dip_hh + 541);

    auto tr_y_xxxyz_yyyzz = pbuffer.data(idx_dip_hh + 542);

    auto tr_y_xxxyz_yyzzz = pbuffer.data(idx_dip_hh + 543);

    auto tr_y_xxxyz_yzzzz = pbuffer.data(idx_dip_hh + 544);

    auto tr_y_xxxyz_zzzzz = pbuffer.data(idx_dip_hh + 545);

    auto tr_y_xxxzz_xxxxx = pbuffer.data(idx_dip_hh + 546);

    auto tr_y_xxxzz_xxxxy = pbuffer.data(idx_dip_hh + 547);

    auto tr_y_xxxzz_xxxxz = pbuffer.data(idx_dip_hh + 548);

    auto tr_y_xxxzz_xxxyy = pbuffer.data(idx_dip_hh + 549);

    auto tr_y_xxxzz_xxxyz = pbuffer.data(idx_dip_hh + 550);

    auto tr_y_xxxzz_xxxzz = pbuffer.data(idx_dip_hh + 551);

    auto tr_y_xxxzz_xxyyy = pbuffer.data(idx_dip_hh + 552);

    auto tr_y_xxxzz_xxyyz = pbuffer.data(idx_dip_hh + 553);

    auto tr_y_xxxzz_xxyzz = pbuffer.data(idx_dip_hh + 554);

    auto tr_y_xxxzz_xxzzz = pbuffer.data(idx_dip_hh + 555);

    auto tr_y_xxxzz_xyyyy = pbuffer.data(idx_dip_hh + 556);

    auto tr_y_xxxzz_xyyyz = pbuffer.data(idx_dip_hh + 557);

    auto tr_y_xxxzz_xyyzz = pbuffer.data(idx_dip_hh + 558);

    auto tr_y_xxxzz_xyzzz = pbuffer.data(idx_dip_hh + 559);

    auto tr_y_xxxzz_xzzzz = pbuffer.data(idx_dip_hh + 560);

    auto tr_y_xxxzz_yyyyy = pbuffer.data(idx_dip_hh + 561);

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

    auto tr_y_xxyyz_xxxxx = pbuffer.data(idx_dip_hh + 588);

    auto tr_y_xxyyz_xxxxy = pbuffer.data(idx_dip_hh + 589);

    auto tr_y_xxyyz_xxxyy = pbuffer.data(idx_dip_hh + 591);

    auto tr_y_xxyyz_xxyyy = pbuffer.data(idx_dip_hh + 594);

    auto tr_y_xxyyz_xyyyy = pbuffer.data(idx_dip_hh + 598);

    auto tr_y_xxyyz_yyyyz = pbuffer.data(idx_dip_hh + 604);

    auto tr_y_xxyyz_yyyzz = pbuffer.data(idx_dip_hh + 605);

    auto tr_y_xxyyz_yyzzz = pbuffer.data(idx_dip_hh + 606);

    auto tr_y_xxyyz_yzzzz = pbuffer.data(idx_dip_hh + 607);

    auto tr_y_xxyyz_zzzzz = pbuffer.data(idx_dip_hh + 608);

    auto tr_y_xxyzz_xxxxy = pbuffer.data(idx_dip_hh + 610);

    auto tr_y_xxyzz_xxxyy = pbuffer.data(idx_dip_hh + 612);

    auto tr_y_xxyzz_xxxyz = pbuffer.data(idx_dip_hh + 613);

    auto tr_y_xxyzz_xxyyy = pbuffer.data(idx_dip_hh + 615);

    auto tr_y_xxyzz_xxyyz = pbuffer.data(idx_dip_hh + 616);

    auto tr_y_xxyzz_xxyzz = pbuffer.data(idx_dip_hh + 617);

    auto tr_y_xxyzz_xyyyy = pbuffer.data(idx_dip_hh + 619);

    auto tr_y_xxyzz_xyyyz = pbuffer.data(idx_dip_hh + 620);

    auto tr_y_xxyzz_xyyzz = pbuffer.data(idx_dip_hh + 621);

    auto tr_y_xxyzz_xyzzz = pbuffer.data(idx_dip_hh + 622);

    auto tr_y_xxyzz_yyyyy = pbuffer.data(idx_dip_hh + 624);

    auto tr_y_xxyzz_yyyyz = pbuffer.data(idx_dip_hh + 625);

    auto tr_y_xxyzz_yyyzz = pbuffer.data(idx_dip_hh + 626);

    auto tr_y_xxyzz_yyzzz = pbuffer.data(idx_dip_hh + 627);

    auto tr_y_xxyzz_yzzzz = pbuffer.data(idx_dip_hh + 628);

    auto tr_y_xxyzz_zzzzz = pbuffer.data(idx_dip_hh + 629);

    auto tr_y_xxzzz_xxxxx = pbuffer.data(idx_dip_hh + 630);

    auto tr_y_xxzzz_xxxxy = pbuffer.data(idx_dip_hh + 631);

    auto tr_y_xxzzz_xxxxz = pbuffer.data(idx_dip_hh + 632);

    auto tr_y_xxzzz_xxxyy = pbuffer.data(idx_dip_hh + 633);

    auto tr_y_xxzzz_xxxyz = pbuffer.data(idx_dip_hh + 634);

    auto tr_y_xxzzz_xxxzz = pbuffer.data(idx_dip_hh + 635);

    auto tr_y_xxzzz_xxyyy = pbuffer.data(idx_dip_hh + 636);

    auto tr_y_xxzzz_xxyyz = pbuffer.data(idx_dip_hh + 637);

    auto tr_y_xxzzz_xxyzz = pbuffer.data(idx_dip_hh + 638);

    auto tr_y_xxzzz_xxzzz = pbuffer.data(idx_dip_hh + 639);

    auto tr_y_xxzzz_xyyyy = pbuffer.data(idx_dip_hh + 640);

    auto tr_y_xxzzz_xyyyz = pbuffer.data(idx_dip_hh + 641);

    auto tr_y_xxzzz_xyyzz = pbuffer.data(idx_dip_hh + 642);

    auto tr_y_xxzzz_xyzzz = pbuffer.data(idx_dip_hh + 643);

    auto tr_y_xxzzz_xzzzz = pbuffer.data(idx_dip_hh + 644);

    auto tr_y_xxzzz_yyyyy = pbuffer.data(idx_dip_hh + 645);

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

    auto tr_y_xyyyz_yyyyz = pbuffer.data(idx_dip_hh + 688);

    auto tr_y_xyyyz_yyyzz = pbuffer.data(idx_dip_hh + 689);

    auto tr_y_xyyyz_yyzzz = pbuffer.data(idx_dip_hh + 690);

    auto tr_y_xyyyz_yzzzz = pbuffer.data(idx_dip_hh + 691);

    auto tr_y_xyyyz_zzzzz = pbuffer.data(idx_dip_hh + 692);

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

    auto tr_y_xyyzz_yyyyy = pbuffer.data(idx_dip_hh + 708);

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

    auto tr_y_xyzzz_yyyyy = pbuffer.data(idx_dip_hh + 729);

    auto tr_y_xyzzz_yyyyz = pbuffer.data(idx_dip_hh + 730);

    auto tr_y_xyzzz_yyyzz = pbuffer.data(idx_dip_hh + 731);

    auto tr_y_xyzzz_yyzzz = pbuffer.data(idx_dip_hh + 732);

    auto tr_y_xyzzz_yzzzz = pbuffer.data(idx_dip_hh + 733);

    auto tr_y_xyzzz_zzzzz = pbuffer.data(idx_dip_hh + 734);

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

    auto tr_y_xzzzz_yyyyy = pbuffer.data(idx_dip_hh + 750);

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

    auto tr_y_yyyyz_xxxxx = pbuffer.data(idx_dip_hh + 777);

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

    auto tr_z_xxxxy_xxxxx = pbuffer.data(idx_dip_hh + 903);

    auto tr_z_xxxxy_xxxxz = pbuffer.data(idx_dip_hh + 905);

    auto tr_z_xxxxy_xxxzz = pbuffer.data(idx_dip_hh + 908);

    auto tr_z_xxxxy_xxzzz = pbuffer.data(idx_dip_hh + 912);

    auto tr_z_xxxxy_xzzzz = pbuffer.data(idx_dip_hh + 917);

    auto tr_z_xxxxy_yyyyy = pbuffer.data(idx_dip_hh + 918);

    auto tr_z_xxxxy_yyyyz = pbuffer.data(idx_dip_hh + 919);

    auto tr_z_xxxxy_yyyzz = pbuffer.data(idx_dip_hh + 920);

    auto tr_z_xxxxy_yyzzz = pbuffer.data(idx_dip_hh + 921);

    auto tr_z_xxxxy_yzzzz = pbuffer.data(idx_dip_hh + 922);

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

    auto tr_z_xxxxz_yyyyy = pbuffer.data(idx_dip_hh + 939);

    auto tr_z_xxxxz_yyyyz = pbuffer.data(idx_dip_hh + 940);

    auto tr_z_xxxxz_yyyzz = pbuffer.data(idx_dip_hh + 941);

    auto tr_z_xxxxz_yyzzz = pbuffer.data(idx_dip_hh + 942);

    auto tr_z_xxxxz_yzzzz = pbuffer.data(idx_dip_hh + 943);

    auto tr_z_xxxxz_zzzzz = pbuffer.data(idx_dip_hh + 944);

    auto tr_z_xxxyy_xxxxx = pbuffer.data(idx_dip_hh + 945);

    auto tr_z_xxxyy_xxxxy = pbuffer.data(idx_dip_hh + 946);

    auto tr_z_xxxyy_xxxxz = pbuffer.data(idx_dip_hh + 947);

    auto tr_z_xxxyy_xxxyy = pbuffer.data(idx_dip_hh + 948);

    auto tr_z_xxxyy_xxxyz = pbuffer.data(idx_dip_hh + 949);

    auto tr_z_xxxyy_xxxzz = pbuffer.data(idx_dip_hh + 950);

    auto tr_z_xxxyy_xxyyy = pbuffer.data(idx_dip_hh + 951);

    auto tr_z_xxxyy_xxyyz = pbuffer.data(idx_dip_hh + 952);

    auto tr_z_xxxyy_xxyzz = pbuffer.data(idx_dip_hh + 953);

    auto tr_z_xxxyy_xxzzz = pbuffer.data(idx_dip_hh + 954);

    auto tr_z_xxxyy_xyyyy = pbuffer.data(idx_dip_hh + 955);

    auto tr_z_xxxyy_xyyyz = pbuffer.data(idx_dip_hh + 956);

    auto tr_z_xxxyy_xyyzz = pbuffer.data(idx_dip_hh + 957);

    auto tr_z_xxxyy_xyzzz = pbuffer.data(idx_dip_hh + 958);

    auto tr_z_xxxyy_xzzzz = pbuffer.data(idx_dip_hh + 959);

    auto tr_z_xxxyy_yyyyy = pbuffer.data(idx_dip_hh + 960);

    auto tr_z_xxxyy_yyyyz = pbuffer.data(idx_dip_hh + 961);

    auto tr_z_xxxyy_yyyzz = pbuffer.data(idx_dip_hh + 962);

    auto tr_z_xxxyy_yyzzz = pbuffer.data(idx_dip_hh + 963);

    auto tr_z_xxxyy_yzzzz = pbuffer.data(idx_dip_hh + 964);

    auto tr_z_xxxyy_zzzzz = pbuffer.data(idx_dip_hh + 965);

    auto tr_z_xxxyz_xxxxx = pbuffer.data(idx_dip_hh + 966);

    auto tr_z_xxxyz_xxxxz = pbuffer.data(idx_dip_hh + 968);

    auto tr_z_xxxyz_xxxzz = pbuffer.data(idx_dip_hh + 971);

    auto tr_z_xxxyz_xxzzz = pbuffer.data(idx_dip_hh + 975);

    auto tr_z_xxxyz_xzzzz = pbuffer.data(idx_dip_hh + 980);

    auto tr_z_xxxyz_yyyyy = pbuffer.data(idx_dip_hh + 981);

    auto tr_z_xxxyz_yyyyz = pbuffer.data(idx_dip_hh + 982);

    auto tr_z_xxxyz_yyyzz = pbuffer.data(idx_dip_hh + 983);

    auto tr_z_xxxyz_yyzzz = pbuffer.data(idx_dip_hh + 984);

    auto tr_z_xxxyz_yzzzz = pbuffer.data(idx_dip_hh + 985);

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

    auto tr_z_xxyyy_xxxxx = pbuffer.data(idx_dip_hh + 1008);

    auto tr_z_xxyyy_xxxxy = pbuffer.data(idx_dip_hh + 1009);

    auto tr_z_xxyyy_xxxxz = pbuffer.data(idx_dip_hh + 1010);

    auto tr_z_xxyyy_xxxyy = pbuffer.data(idx_dip_hh + 1011);

    auto tr_z_xxyyy_xxxyz = pbuffer.data(idx_dip_hh + 1012);

    auto tr_z_xxyyy_xxxzz = pbuffer.data(idx_dip_hh + 1013);

    auto tr_z_xxyyy_xxyyy = pbuffer.data(idx_dip_hh + 1014);

    auto tr_z_xxyyy_xxyyz = pbuffer.data(idx_dip_hh + 1015);

    auto tr_z_xxyyy_xxyzz = pbuffer.data(idx_dip_hh + 1016);

    auto tr_z_xxyyy_xxzzz = pbuffer.data(idx_dip_hh + 1017);

    auto tr_z_xxyyy_xyyyy = pbuffer.data(idx_dip_hh + 1018);

    auto tr_z_xxyyy_xyyyz = pbuffer.data(idx_dip_hh + 1019);

    auto tr_z_xxyyy_xyyzz = pbuffer.data(idx_dip_hh + 1020);

    auto tr_z_xxyyy_xyzzz = pbuffer.data(idx_dip_hh + 1021);

    auto tr_z_xxyyy_xzzzz = pbuffer.data(idx_dip_hh + 1022);

    auto tr_z_xxyyy_yyyyy = pbuffer.data(idx_dip_hh + 1023);

    auto tr_z_xxyyy_yyyyz = pbuffer.data(idx_dip_hh + 1024);

    auto tr_z_xxyyy_yyyzz = pbuffer.data(idx_dip_hh + 1025);

    auto tr_z_xxyyy_yyzzz = pbuffer.data(idx_dip_hh + 1026);

    auto tr_z_xxyyy_yzzzz = pbuffer.data(idx_dip_hh + 1027);

    auto tr_z_xxyyy_zzzzz = pbuffer.data(idx_dip_hh + 1028);

    auto tr_z_xxyyz_xxxxx = pbuffer.data(idx_dip_hh + 1029);

    auto tr_z_xxyyz_xxxxz = pbuffer.data(idx_dip_hh + 1031);

    auto tr_z_xxyyz_xxxyz = pbuffer.data(idx_dip_hh + 1033);

    auto tr_z_xxyyz_xxxzz = pbuffer.data(idx_dip_hh + 1034);

    auto tr_z_xxyyz_xxyyz = pbuffer.data(idx_dip_hh + 1036);

    auto tr_z_xxyyz_xxyzz = pbuffer.data(idx_dip_hh + 1037);

    auto tr_z_xxyyz_xxzzz = pbuffer.data(idx_dip_hh + 1038);

    auto tr_z_xxyyz_xyyyz = pbuffer.data(idx_dip_hh + 1040);

    auto tr_z_xxyyz_xyyzz = pbuffer.data(idx_dip_hh + 1041);

    auto tr_z_xxyyz_xyzzz = pbuffer.data(idx_dip_hh + 1042);

    auto tr_z_xxyyz_xzzzz = pbuffer.data(idx_dip_hh + 1043);

    auto tr_z_xxyyz_yyyyy = pbuffer.data(idx_dip_hh + 1044);

    auto tr_z_xxyyz_yyyyz = pbuffer.data(idx_dip_hh + 1045);

    auto tr_z_xxyyz_yyyzz = pbuffer.data(idx_dip_hh + 1046);

    auto tr_z_xxyyz_yyzzz = pbuffer.data(idx_dip_hh + 1047);

    auto tr_z_xxyyz_yzzzz = pbuffer.data(idx_dip_hh + 1048);

    auto tr_z_xxyyz_zzzzz = pbuffer.data(idx_dip_hh + 1049);

    auto tr_z_xxyzz_xxxxx = pbuffer.data(idx_dip_hh + 1050);

    auto tr_z_xxyzz_xxxxz = pbuffer.data(idx_dip_hh + 1052);

    auto tr_z_xxyzz_xxxzz = pbuffer.data(idx_dip_hh + 1055);

    auto tr_z_xxyzz_xxzzz = pbuffer.data(idx_dip_hh + 1059);

    auto tr_z_xxyzz_xzzzz = pbuffer.data(idx_dip_hh + 1064);

    auto tr_z_xxyzz_yyyyy = pbuffer.data(idx_dip_hh + 1065);

    auto tr_z_xxyzz_yyyyz = pbuffer.data(idx_dip_hh + 1066);

    auto tr_z_xxyzz_yyyzz = pbuffer.data(idx_dip_hh + 1067);

    auto tr_z_xxyzz_yyzzz = pbuffer.data(idx_dip_hh + 1068);

    auto tr_z_xxyzz_yzzzz = pbuffer.data(idx_dip_hh + 1069);

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

    auto tr_z_xyyyy_zzzzz = pbuffer.data(idx_dip_hh + 1112);

    auto tr_z_xyyyz_xxxyz = pbuffer.data(idx_dip_hh + 1117);

    auto tr_z_xyyyz_xxyyz = pbuffer.data(idx_dip_hh + 1120);

    auto tr_z_xyyyz_xxyzz = pbuffer.data(idx_dip_hh + 1121);

    auto tr_z_xyyyz_xyyyz = pbuffer.data(idx_dip_hh + 1124);

    auto tr_z_xyyyz_xyyzz = pbuffer.data(idx_dip_hh + 1125);

    auto tr_z_xyyyz_xyzzz = pbuffer.data(idx_dip_hh + 1126);

    auto tr_z_xyyyz_yyyyy = pbuffer.data(idx_dip_hh + 1128);

    auto tr_z_xyyyz_yyyyz = pbuffer.data(idx_dip_hh + 1129);

    auto tr_z_xyyyz_yyyzz = pbuffer.data(idx_dip_hh + 1130);

    auto tr_z_xyyyz_yyzzz = pbuffer.data(idx_dip_hh + 1131);

    auto tr_z_xyyyz_yzzzz = pbuffer.data(idx_dip_hh + 1132);

    auto tr_z_xyyyz_zzzzz = pbuffer.data(idx_dip_hh + 1133);

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

    auto tr_z_xyyzz_zzzzz = pbuffer.data(idx_dip_hh + 1154);

    auto tr_z_xyzzz_yyyyy = pbuffer.data(idx_dip_hh + 1170);

    auto tr_z_xyzzz_yyyyz = pbuffer.data(idx_dip_hh + 1171);

    auto tr_z_xyzzz_yyyzz = pbuffer.data(idx_dip_hh + 1172);

    auto tr_z_xyzzz_yyzzz = pbuffer.data(idx_dip_hh + 1173);

    auto tr_z_xyzzz_yzzzz = pbuffer.data(idx_dip_hh + 1174);

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

    // Set up 0-21 components of targeted buffer : IH

    auto tr_x_xxxxxx_xxxxx = pbuffer.data(idx_dip_ih);

    auto tr_x_xxxxxx_xxxxy = pbuffer.data(idx_dip_ih + 1);

    auto tr_x_xxxxxx_xxxxz = pbuffer.data(idx_dip_ih + 2);

    auto tr_x_xxxxxx_xxxyy = pbuffer.data(idx_dip_ih + 3);

    auto tr_x_xxxxxx_xxxyz = pbuffer.data(idx_dip_ih + 4);

    auto tr_x_xxxxxx_xxxzz = pbuffer.data(idx_dip_ih + 5);

    auto tr_x_xxxxxx_xxyyy = pbuffer.data(idx_dip_ih + 6);

    auto tr_x_xxxxxx_xxyyz = pbuffer.data(idx_dip_ih + 7);

    auto tr_x_xxxxxx_xxyzz = pbuffer.data(idx_dip_ih + 8);

    auto tr_x_xxxxxx_xxzzz = pbuffer.data(idx_dip_ih + 9);

    auto tr_x_xxxxxx_xyyyy = pbuffer.data(idx_dip_ih + 10);

    auto tr_x_xxxxxx_xyyyz = pbuffer.data(idx_dip_ih + 11);

    auto tr_x_xxxxxx_xyyzz = pbuffer.data(idx_dip_ih + 12);

    auto tr_x_xxxxxx_xyzzz = pbuffer.data(idx_dip_ih + 13);

    auto tr_x_xxxxxx_xzzzz = pbuffer.data(idx_dip_ih + 14);

    auto tr_x_xxxxxx_yyyyy = pbuffer.data(idx_dip_ih + 15);

    auto tr_x_xxxxxx_yyyyz = pbuffer.data(idx_dip_ih + 16);

    auto tr_x_xxxxxx_yyyzz = pbuffer.data(idx_dip_ih + 17);

    auto tr_x_xxxxxx_yyzzz = pbuffer.data(idx_dip_ih + 18);

    auto tr_x_xxxxxx_yzzzz = pbuffer.data(idx_dip_ih + 19);

    auto tr_x_xxxxxx_zzzzz = pbuffer.data(idx_dip_ih + 20);

    #pragma omp simd aligned(pa_x, tr_x_xxxx_xxxxx, tr_x_xxxx_xxxxy, tr_x_xxxx_xxxxz, tr_x_xxxx_xxxyy, tr_x_xxxx_xxxyz, tr_x_xxxx_xxxzz, tr_x_xxxx_xxyyy, tr_x_xxxx_xxyyz, tr_x_xxxx_xxyzz, tr_x_xxxx_xxzzz, tr_x_xxxx_xyyyy, tr_x_xxxx_xyyyz, tr_x_xxxx_xyyzz, tr_x_xxxx_xyzzz, tr_x_xxxx_xzzzz, tr_x_xxxx_yyyyy, tr_x_xxxx_yyyyz, tr_x_xxxx_yyyzz, tr_x_xxxx_yyzzz, tr_x_xxxx_yzzzz, tr_x_xxxx_zzzzz, tr_x_xxxxx_xxxx, tr_x_xxxxx_xxxxx, tr_x_xxxxx_xxxxy, tr_x_xxxxx_xxxxz, tr_x_xxxxx_xxxy, tr_x_xxxxx_xxxyy, tr_x_xxxxx_xxxyz, tr_x_xxxxx_xxxz, tr_x_xxxxx_xxxzz, tr_x_xxxxx_xxyy, tr_x_xxxxx_xxyyy, tr_x_xxxxx_xxyyz, tr_x_xxxxx_xxyz, tr_x_xxxxx_xxyzz, tr_x_xxxxx_xxzz, tr_x_xxxxx_xxzzz, tr_x_xxxxx_xyyy, tr_x_xxxxx_xyyyy, tr_x_xxxxx_xyyyz, tr_x_xxxxx_xyyz, tr_x_xxxxx_xyyzz, tr_x_xxxxx_xyzz, tr_x_xxxxx_xyzzz, tr_x_xxxxx_xzzz, tr_x_xxxxx_xzzzz, tr_x_xxxxx_yyyy, tr_x_xxxxx_yyyyy, tr_x_xxxxx_yyyyz, tr_x_xxxxx_yyyz, tr_x_xxxxx_yyyzz, tr_x_xxxxx_yyzz, tr_x_xxxxx_yyzzz, tr_x_xxxxx_yzzz, tr_x_xxxxx_yzzzz, tr_x_xxxxx_zzzz, tr_x_xxxxx_zzzzz, tr_x_xxxxxx_xxxxx, tr_x_xxxxxx_xxxxy, tr_x_xxxxxx_xxxxz, tr_x_xxxxxx_xxxyy, tr_x_xxxxxx_xxxyz, tr_x_xxxxxx_xxxzz, tr_x_xxxxxx_xxyyy, tr_x_xxxxxx_xxyyz, tr_x_xxxxxx_xxyzz, tr_x_xxxxxx_xxzzz, tr_x_xxxxxx_xyyyy, tr_x_xxxxxx_xyyyz, tr_x_xxxxxx_xyyzz, tr_x_xxxxxx_xyzzz, tr_x_xxxxxx_xzzzz, tr_x_xxxxxx_yyyyy, tr_x_xxxxxx_yyyyz, tr_x_xxxxxx_yyyzz, tr_x_xxxxxx_yyzzz, tr_x_xxxxxx_yzzzz, tr_x_xxxxxx_zzzzz, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxzz, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyyy, ts_xxxxx_yyyyz, ts_xxxxx_yyyzz, ts_xxxxx_yyzzz, ts_xxxxx_yzzzz, ts_xxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_xxxxx[i] = 5.0 * tr_x_xxxx_xxxxx[i] * fe_0 + 5.0 * tr_x_xxxxx_xxxx[i] * fe_0 + ts_xxxxx_xxxxx[i] * fe_0 + tr_x_xxxxx_xxxxx[i] * pa_x[i];

        tr_x_xxxxxx_xxxxy[i] = 5.0 * tr_x_xxxx_xxxxy[i] * fe_0 + 4.0 * tr_x_xxxxx_xxxy[i] * fe_0 + ts_xxxxx_xxxxy[i] * fe_0 + tr_x_xxxxx_xxxxy[i] * pa_x[i];

        tr_x_xxxxxx_xxxxz[i] = 5.0 * tr_x_xxxx_xxxxz[i] * fe_0 + 4.0 * tr_x_xxxxx_xxxz[i] * fe_0 + ts_xxxxx_xxxxz[i] * fe_0 + tr_x_xxxxx_xxxxz[i] * pa_x[i];

        tr_x_xxxxxx_xxxyy[i] = 5.0 * tr_x_xxxx_xxxyy[i] * fe_0 + 3.0 * tr_x_xxxxx_xxyy[i] * fe_0 + ts_xxxxx_xxxyy[i] * fe_0 + tr_x_xxxxx_xxxyy[i] * pa_x[i];

        tr_x_xxxxxx_xxxyz[i] = 5.0 * tr_x_xxxx_xxxyz[i] * fe_0 + 3.0 * tr_x_xxxxx_xxyz[i] * fe_0 + ts_xxxxx_xxxyz[i] * fe_0 + tr_x_xxxxx_xxxyz[i] * pa_x[i];

        tr_x_xxxxxx_xxxzz[i] = 5.0 * tr_x_xxxx_xxxzz[i] * fe_0 + 3.0 * tr_x_xxxxx_xxzz[i] * fe_0 + ts_xxxxx_xxxzz[i] * fe_0 + tr_x_xxxxx_xxxzz[i] * pa_x[i];

        tr_x_xxxxxx_xxyyy[i] = 5.0 * tr_x_xxxx_xxyyy[i] * fe_0 + 2.0 * tr_x_xxxxx_xyyy[i] * fe_0 + ts_xxxxx_xxyyy[i] * fe_0 + tr_x_xxxxx_xxyyy[i] * pa_x[i];

        tr_x_xxxxxx_xxyyz[i] = 5.0 * tr_x_xxxx_xxyyz[i] * fe_0 + 2.0 * tr_x_xxxxx_xyyz[i] * fe_0 + ts_xxxxx_xxyyz[i] * fe_0 + tr_x_xxxxx_xxyyz[i] * pa_x[i];

        tr_x_xxxxxx_xxyzz[i] = 5.0 * tr_x_xxxx_xxyzz[i] * fe_0 + 2.0 * tr_x_xxxxx_xyzz[i] * fe_0 + ts_xxxxx_xxyzz[i] * fe_0 + tr_x_xxxxx_xxyzz[i] * pa_x[i];

        tr_x_xxxxxx_xxzzz[i] = 5.0 * tr_x_xxxx_xxzzz[i] * fe_0 + 2.0 * tr_x_xxxxx_xzzz[i] * fe_0 + ts_xxxxx_xxzzz[i] * fe_0 + tr_x_xxxxx_xxzzz[i] * pa_x[i];

        tr_x_xxxxxx_xyyyy[i] = 5.0 * tr_x_xxxx_xyyyy[i] * fe_0 + tr_x_xxxxx_yyyy[i] * fe_0 + ts_xxxxx_xyyyy[i] * fe_0 + tr_x_xxxxx_xyyyy[i] * pa_x[i];

        tr_x_xxxxxx_xyyyz[i] = 5.0 * tr_x_xxxx_xyyyz[i] * fe_0 + tr_x_xxxxx_yyyz[i] * fe_0 + ts_xxxxx_xyyyz[i] * fe_0 + tr_x_xxxxx_xyyyz[i] * pa_x[i];

        tr_x_xxxxxx_xyyzz[i] = 5.0 * tr_x_xxxx_xyyzz[i] * fe_0 + tr_x_xxxxx_yyzz[i] * fe_0 + ts_xxxxx_xyyzz[i] * fe_0 + tr_x_xxxxx_xyyzz[i] * pa_x[i];

        tr_x_xxxxxx_xyzzz[i] = 5.0 * tr_x_xxxx_xyzzz[i] * fe_0 + tr_x_xxxxx_yzzz[i] * fe_0 + ts_xxxxx_xyzzz[i] * fe_0 + tr_x_xxxxx_xyzzz[i] * pa_x[i];

        tr_x_xxxxxx_xzzzz[i] = 5.0 * tr_x_xxxx_xzzzz[i] * fe_0 + tr_x_xxxxx_zzzz[i] * fe_0 + ts_xxxxx_xzzzz[i] * fe_0 + tr_x_xxxxx_xzzzz[i] * pa_x[i];

        tr_x_xxxxxx_yyyyy[i] = 5.0 * tr_x_xxxx_yyyyy[i] * fe_0 + ts_xxxxx_yyyyy[i] * fe_0 + tr_x_xxxxx_yyyyy[i] * pa_x[i];

        tr_x_xxxxxx_yyyyz[i] = 5.0 * tr_x_xxxx_yyyyz[i] * fe_0 + ts_xxxxx_yyyyz[i] * fe_0 + tr_x_xxxxx_yyyyz[i] * pa_x[i];

        tr_x_xxxxxx_yyyzz[i] = 5.0 * tr_x_xxxx_yyyzz[i] * fe_0 + ts_xxxxx_yyyzz[i] * fe_0 + tr_x_xxxxx_yyyzz[i] * pa_x[i];

        tr_x_xxxxxx_yyzzz[i] = 5.0 * tr_x_xxxx_yyzzz[i] * fe_0 + ts_xxxxx_yyzzz[i] * fe_0 + tr_x_xxxxx_yyzzz[i] * pa_x[i];

        tr_x_xxxxxx_yzzzz[i] = 5.0 * tr_x_xxxx_yzzzz[i] * fe_0 + ts_xxxxx_yzzzz[i] * fe_0 + tr_x_xxxxx_yzzzz[i] * pa_x[i];

        tr_x_xxxxxx_zzzzz[i] = 5.0 * tr_x_xxxx_zzzzz[i] * fe_0 + ts_xxxxx_zzzzz[i] * fe_0 + tr_x_xxxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : IH

    auto tr_x_xxxxxy_xxxxx = pbuffer.data(idx_dip_ih + 21);

    auto tr_x_xxxxxy_xxxxy = pbuffer.data(idx_dip_ih + 22);

    auto tr_x_xxxxxy_xxxxz = pbuffer.data(idx_dip_ih + 23);

    auto tr_x_xxxxxy_xxxyy = pbuffer.data(idx_dip_ih + 24);

    auto tr_x_xxxxxy_xxxyz = pbuffer.data(idx_dip_ih + 25);

    auto tr_x_xxxxxy_xxxzz = pbuffer.data(idx_dip_ih + 26);

    auto tr_x_xxxxxy_xxyyy = pbuffer.data(idx_dip_ih + 27);

    auto tr_x_xxxxxy_xxyyz = pbuffer.data(idx_dip_ih + 28);

    auto tr_x_xxxxxy_xxyzz = pbuffer.data(idx_dip_ih + 29);

    auto tr_x_xxxxxy_xxzzz = pbuffer.data(idx_dip_ih + 30);

    auto tr_x_xxxxxy_xyyyy = pbuffer.data(idx_dip_ih + 31);

    auto tr_x_xxxxxy_xyyyz = pbuffer.data(idx_dip_ih + 32);

    auto tr_x_xxxxxy_xyyzz = pbuffer.data(idx_dip_ih + 33);

    auto tr_x_xxxxxy_xyzzz = pbuffer.data(idx_dip_ih + 34);

    auto tr_x_xxxxxy_xzzzz = pbuffer.data(idx_dip_ih + 35);

    auto tr_x_xxxxxy_yyyyy = pbuffer.data(idx_dip_ih + 36);

    auto tr_x_xxxxxy_yyyyz = pbuffer.data(idx_dip_ih + 37);

    auto tr_x_xxxxxy_yyyzz = pbuffer.data(idx_dip_ih + 38);

    auto tr_x_xxxxxy_yyzzz = pbuffer.data(idx_dip_ih + 39);

    auto tr_x_xxxxxy_yzzzz = pbuffer.data(idx_dip_ih + 40);

    auto tr_x_xxxxxy_zzzzz = pbuffer.data(idx_dip_ih + 41);

    #pragma omp simd aligned(pa_y, tr_x_xxxxx_xxxx, tr_x_xxxxx_xxxxx, tr_x_xxxxx_xxxxy, tr_x_xxxxx_xxxxz, tr_x_xxxxx_xxxy, tr_x_xxxxx_xxxyy, tr_x_xxxxx_xxxyz, tr_x_xxxxx_xxxz, tr_x_xxxxx_xxxzz, tr_x_xxxxx_xxyy, tr_x_xxxxx_xxyyy, tr_x_xxxxx_xxyyz, tr_x_xxxxx_xxyz, tr_x_xxxxx_xxyzz, tr_x_xxxxx_xxzz, tr_x_xxxxx_xxzzz, tr_x_xxxxx_xyyy, tr_x_xxxxx_xyyyy, tr_x_xxxxx_xyyyz, tr_x_xxxxx_xyyz, tr_x_xxxxx_xyyzz, tr_x_xxxxx_xyzz, tr_x_xxxxx_xyzzz, tr_x_xxxxx_xzzz, tr_x_xxxxx_xzzzz, tr_x_xxxxx_yyyy, tr_x_xxxxx_yyyyy, tr_x_xxxxx_yyyyz, tr_x_xxxxx_yyyz, tr_x_xxxxx_yyyzz, tr_x_xxxxx_yyzz, tr_x_xxxxx_yyzzz, tr_x_xxxxx_yzzz, tr_x_xxxxx_yzzzz, tr_x_xxxxx_zzzz, tr_x_xxxxx_zzzzz, tr_x_xxxxxy_xxxxx, tr_x_xxxxxy_xxxxy, tr_x_xxxxxy_xxxxz, tr_x_xxxxxy_xxxyy, tr_x_xxxxxy_xxxyz, tr_x_xxxxxy_xxxzz, tr_x_xxxxxy_xxyyy, tr_x_xxxxxy_xxyyz, tr_x_xxxxxy_xxyzz, tr_x_xxxxxy_xxzzz, tr_x_xxxxxy_xyyyy, tr_x_xxxxxy_xyyyz, tr_x_xxxxxy_xyyzz, tr_x_xxxxxy_xyzzz, tr_x_xxxxxy_xzzzz, tr_x_xxxxxy_yyyyy, tr_x_xxxxxy_yyyyz, tr_x_xxxxxy_yyyzz, tr_x_xxxxxy_yyzzz, tr_x_xxxxxy_yzzzz, tr_x_xxxxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxy_xxxxx[i] = tr_x_xxxxx_xxxxx[i] * pa_y[i];

        tr_x_xxxxxy_xxxxy[i] = tr_x_xxxxx_xxxx[i] * fe_0 + tr_x_xxxxx_xxxxy[i] * pa_y[i];

        tr_x_xxxxxy_xxxxz[i] = tr_x_xxxxx_xxxxz[i] * pa_y[i];

        tr_x_xxxxxy_xxxyy[i] = 2.0 * tr_x_xxxxx_xxxy[i] * fe_0 + tr_x_xxxxx_xxxyy[i] * pa_y[i];

        tr_x_xxxxxy_xxxyz[i] = tr_x_xxxxx_xxxz[i] * fe_0 + tr_x_xxxxx_xxxyz[i] * pa_y[i];

        tr_x_xxxxxy_xxxzz[i] = tr_x_xxxxx_xxxzz[i] * pa_y[i];

        tr_x_xxxxxy_xxyyy[i] = 3.0 * tr_x_xxxxx_xxyy[i] * fe_0 + tr_x_xxxxx_xxyyy[i] * pa_y[i];

        tr_x_xxxxxy_xxyyz[i] = 2.0 * tr_x_xxxxx_xxyz[i] * fe_0 + tr_x_xxxxx_xxyyz[i] * pa_y[i];

        tr_x_xxxxxy_xxyzz[i] = tr_x_xxxxx_xxzz[i] * fe_0 + tr_x_xxxxx_xxyzz[i] * pa_y[i];

        tr_x_xxxxxy_xxzzz[i] = tr_x_xxxxx_xxzzz[i] * pa_y[i];

        tr_x_xxxxxy_xyyyy[i] = 4.0 * tr_x_xxxxx_xyyy[i] * fe_0 + tr_x_xxxxx_xyyyy[i] * pa_y[i];

        tr_x_xxxxxy_xyyyz[i] = 3.0 * tr_x_xxxxx_xyyz[i] * fe_0 + tr_x_xxxxx_xyyyz[i] * pa_y[i];

        tr_x_xxxxxy_xyyzz[i] = 2.0 * tr_x_xxxxx_xyzz[i] * fe_0 + tr_x_xxxxx_xyyzz[i] * pa_y[i];

        tr_x_xxxxxy_xyzzz[i] = tr_x_xxxxx_xzzz[i] * fe_0 + tr_x_xxxxx_xyzzz[i] * pa_y[i];

        tr_x_xxxxxy_xzzzz[i] = tr_x_xxxxx_xzzzz[i] * pa_y[i];

        tr_x_xxxxxy_yyyyy[i] = 5.0 * tr_x_xxxxx_yyyy[i] * fe_0 + tr_x_xxxxx_yyyyy[i] * pa_y[i];

        tr_x_xxxxxy_yyyyz[i] = 4.0 * tr_x_xxxxx_yyyz[i] * fe_0 + tr_x_xxxxx_yyyyz[i] * pa_y[i];

        tr_x_xxxxxy_yyyzz[i] = 3.0 * tr_x_xxxxx_yyzz[i] * fe_0 + tr_x_xxxxx_yyyzz[i] * pa_y[i];

        tr_x_xxxxxy_yyzzz[i] = 2.0 * tr_x_xxxxx_yzzz[i] * fe_0 + tr_x_xxxxx_yyzzz[i] * pa_y[i];

        tr_x_xxxxxy_yzzzz[i] = tr_x_xxxxx_zzzz[i] * fe_0 + tr_x_xxxxx_yzzzz[i] * pa_y[i];

        tr_x_xxxxxy_zzzzz[i] = tr_x_xxxxx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : IH

    auto tr_x_xxxxxz_xxxxx = pbuffer.data(idx_dip_ih + 42);

    auto tr_x_xxxxxz_xxxxy = pbuffer.data(idx_dip_ih + 43);

    auto tr_x_xxxxxz_xxxxz = pbuffer.data(idx_dip_ih + 44);

    auto tr_x_xxxxxz_xxxyy = pbuffer.data(idx_dip_ih + 45);

    auto tr_x_xxxxxz_xxxyz = pbuffer.data(idx_dip_ih + 46);

    auto tr_x_xxxxxz_xxxzz = pbuffer.data(idx_dip_ih + 47);

    auto tr_x_xxxxxz_xxyyy = pbuffer.data(idx_dip_ih + 48);

    auto tr_x_xxxxxz_xxyyz = pbuffer.data(idx_dip_ih + 49);

    auto tr_x_xxxxxz_xxyzz = pbuffer.data(idx_dip_ih + 50);

    auto tr_x_xxxxxz_xxzzz = pbuffer.data(idx_dip_ih + 51);

    auto tr_x_xxxxxz_xyyyy = pbuffer.data(idx_dip_ih + 52);

    auto tr_x_xxxxxz_xyyyz = pbuffer.data(idx_dip_ih + 53);

    auto tr_x_xxxxxz_xyyzz = pbuffer.data(idx_dip_ih + 54);

    auto tr_x_xxxxxz_xyzzz = pbuffer.data(idx_dip_ih + 55);

    auto tr_x_xxxxxz_xzzzz = pbuffer.data(idx_dip_ih + 56);

    auto tr_x_xxxxxz_yyyyy = pbuffer.data(idx_dip_ih + 57);

    auto tr_x_xxxxxz_yyyyz = pbuffer.data(idx_dip_ih + 58);

    auto tr_x_xxxxxz_yyyzz = pbuffer.data(idx_dip_ih + 59);

    auto tr_x_xxxxxz_yyzzz = pbuffer.data(idx_dip_ih + 60);

    auto tr_x_xxxxxz_yzzzz = pbuffer.data(idx_dip_ih + 61);

    auto tr_x_xxxxxz_zzzzz = pbuffer.data(idx_dip_ih + 62);

    #pragma omp simd aligned(pa_z, tr_x_xxxxx_xxxx, tr_x_xxxxx_xxxxx, tr_x_xxxxx_xxxxy, tr_x_xxxxx_xxxxz, tr_x_xxxxx_xxxy, tr_x_xxxxx_xxxyy, tr_x_xxxxx_xxxyz, tr_x_xxxxx_xxxz, tr_x_xxxxx_xxxzz, tr_x_xxxxx_xxyy, tr_x_xxxxx_xxyyy, tr_x_xxxxx_xxyyz, tr_x_xxxxx_xxyz, tr_x_xxxxx_xxyzz, tr_x_xxxxx_xxzz, tr_x_xxxxx_xxzzz, tr_x_xxxxx_xyyy, tr_x_xxxxx_xyyyy, tr_x_xxxxx_xyyyz, tr_x_xxxxx_xyyz, tr_x_xxxxx_xyyzz, tr_x_xxxxx_xyzz, tr_x_xxxxx_xyzzz, tr_x_xxxxx_xzzz, tr_x_xxxxx_xzzzz, tr_x_xxxxx_yyyy, tr_x_xxxxx_yyyyy, tr_x_xxxxx_yyyyz, tr_x_xxxxx_yyyz, tr_x_xxxxx_yyyzz, tr_x_xxxxx_yyzz, tr_x_xxxxx_yyzzz, tr_x_xxxxx_yzzz, tr_x_xxxxx_yzzzz, tr_x_xxxxx_zzzz, tr_x_xxxxx_zzzzz, tr_x_xxxxxz_xxxxx, tr_x_xxxxxz_xxxxy, tr_x_xxxxxz_xxxxz, tr_x_xxxxxz_xxxyy, tr_x_xxxxxz_xxxyz, tr_x_xxxxxz_xxxzz, tr_x_xxxxxz_xxyyy, tr_x_xxxxxz_xxyyz, tr_x_xxxxxz_xxyzz, tr_x_xxxxxz_xxzzz, tr_x_xxxxxz_xyyyy, tr_x_xxxxxz_xyyyz, tr_x_xxxxxz_xyyzz, tr_x_xxxxxz_xyzzz, tr_x_xxxxxz_xzzzz, tr_x_xxxxxz_yyyyy, tr_x_xxxxxz_yyyyz, tr_x_xxxxxz_yyyzz, tr_x_xxxxxz_yyzzz, tr_x_xxxxxz_yzzzz, tr_x_xxxxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxz_xxxxx[i] = tr_x_xxxxx_xxxxx[i] * pa_z[i];

        tr_x_xxxxxz_xxxxy[i] = tr_x_xxxxx_xxxxy[i] * pa_z[i];

        tr_x_xxxxxz_xxxxz[i] = tr_x_xxxxx_xxxx[i] * fe_0 + tr_x_xxxxx_xxxxz[i] * pa_z[i];

        tr_x_xxxxxz_xxxyy[i] = tr_x_xxxxx_xxxyy[i] * pa_z[i];

        tr_x_xxxxxz_xxxyz[i] = tr_x_xxxxx_xxxy[i] * fe_0 + tr_x_xxxxx_xxxyz[i] * pa_z[i];

        tr_x_xxxxxz_xxxzz[i] = 2.0 * tr_x_xxxxx_xxxz[i] * fe_0 + tr_x_xxxxx_xxxzz[i] * pa_z[i];

        tr_x_xxxxxz_xxyyy[i] = tr_x_xxxxx_xxyyy[i] * pa_z[i];

        tr_x_xxxxxz_xxyyz[i] = tr_x_xxxxx_xxyy[i] * fe_0 + tr_x_xxxxx_xxyyz[i] * pa_z[i];

        tr_x_xxxxxz_xxyzz[i] = 2.0 * tr_x_xxxxx_xxyz[i] * fe_0 + tr_x_xxxxx_xxyzz[i] * pa_z[i];

        tr_x_xxxxxz_xxzzz[i] = 3.0 * tr_x_xxxxx_xxzz[i] * fe_0 + tr_x_xxxxx_xxzzz[i] * pa_z[i];

        tr_x_xxxxxz_xyyyy[i] = tr_x_xxxxx_xyyyy[i] * pa_z[i];

        tr_x_xxxxxz_xyyyz[i] = tr_x_xxxxx_xyyy[i] * fe_0 + tr_x_xxxxx_xyyyz[i] * pa_z[i];

        tr_x_xxxxxz_xyyzz[i] = 2.0 * tr_x_xxxxx_xyyz[i] * fe_0 + tr_x_xxxxx_xyyzz[i] * pa_z[i];

        tr_x_xxxxxz_xyzzz[i] = 3.0 * tr_x_xxxxx_xyzz[i] * fe_0 + tr_x_xxxxx_xyzzz[i] * pa_z[i];

        tr_x_xxxxxz_xzzzz[i] = 4.0 * tr_x_xxxxx_xzzz[i] * fe_0 + tr_x_xxxxx_xzzzz[i] * pa_z[i];

        tr_x_xxxxxz_yyyyy[i] = tr_x_xxxxx_yyyyy[i] * pa_z[i];

        tr_x_xxxxxz_yyyyz[i] = tr_x_xxxxx_yyyy[i] * fe_0 + tr_x_xxxxx_yyyyz[i] * pa_z[i];

        tr_x_xxxxxz_yyyzz[i] = 2.0 * tr_x_xxxxx_yyyz[i] * fe_0 + tr_x_xxxxx_yyyzz[i] * pa_z[i];

        tr_x_xxxxxz_yyzzz[i] = 3.0 * tr_x_xxxxx_yyzz[i] * fe_0 + tr_x_xxxxx_yyzzz[i] * pa_z[i];

        tr_x_xxxxxz_yzzzz[i] = 4.0 * tr_x_xxxxx_yzzz[i] * fe_0 + tr_x_xxxxx_yzzzz[i] * pa_z[i];

        tr_x_xxxxxz_zzzzz[i] = 5.0 * tr_x_xxxxx_zzzz[i] * fe_0 + tr_x_xxxxx_zzzzz[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : IH

    auto tr_x_xxxxyy_xxxxx = pbuffer.data(idx_dip_ih + 63);

    auto tr_x_xxxxyy_xxxxy = pbuffer.data(idx_dip_ih + 64);

    auto tr_x_xxxxyy_xxxxz = pbuffer.data(idx_dip_ih + 65);

    auto tr_x_xxxxyy_xxxyy = pbuffer.data(idx_dip_ih + 66);

    auto tr_x_xxxxyy_xxxyz = pbuffer.data(idx_dip_ih + 67);

    auto tr_x_xxxxyy_xxxzz = pbuffer.data(idx_dip_ih + 68);

    auto tr_x_xxxxyy_xxyyy = pbuffer.data(idx_dip_ih + 69);

    auto tr_x_xxxxyy_xxyyz = pbuffer.data(idx_dip_ih + 70);

    auto tr_x_xxxxyy_xxyzz = pbuffer.data(idx_dip_ih + 71);

    auto tr_x_xxxxyy_xxzzz = pbuffer.data(idx_dip_ih + 72);

    auto tr_x_xxxxyy_xyyyy = pbuffer.data(idx_dip_ih + 73);

    auto tr_x_xxxxyy_xyyyz = pbuffer.data(idx_dip_ih + 74);

    auto tr_x_xxxxyy_xyyzz = pbuffer.data(idx_dip_ih + 75);

    auto tr_x_xxxxyy_xyzzz = pbuffer.data(idx_dip_ih + 76);

    auto tr_x_xxxxyy_xzzzz = pbuffer.data(idx_dip_ih + 77);

    auto tr_x_xxxxyy_yyyyy = pbuffer.data(idx_dip_ih + 78);

    auto tr_x_xxxxyy_yyyyz = pbuffer.data(idx_dip_ih + 79);

    auto tr_x_xxxxyy_yyyzz = pbuffer.data(idx_dip_ih + 80);

    auto tr_x_xxxxyy_yyzzz = pbuffer.data(idx_dip_ih + 81);

    auto tr_x_xxxxyy_yzzzz = pbuffer.data(idx_dip_ih + 82);

    auto tr_x_xxxxyy_zzzzz = pbuffer.data(idx_dip_ih + 83);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxxx_xxxxx, tr_x_xxxx_xxxxy, tr_x_xxxx_xxxxz, tr_x_xxxx_xxxyy, tr_x_xxxx_xxxyz, tr_x_xxxx_xxxzz, tr_x_xxxx_xxyyy, tr_x_xxxx_xxyyz, tr_x_xxxx_xxyzz, tr_x_xxxx_xxzzz, tr_x_xxxx_xyyyy, tr_x_xxxx_xyyyz, tr_x_xxxx_xyyzz, tr_x_xxxx_xyzzz, tr_x_xxxx_xzzzz, tr_x_xxxx_zzzzz, tr_x_xxxxy_xxxx, tr_x_xxxxy_xxxxx, tr_x_xxxxy_xxxxy, tr_x_xxxxy_xxxxz, tr_x_xxxxy_xxxy, tr_x_xxxxy_xxxyy, tr_x_xxxxy_xxxyz, tr_x_xxxxy_xxxz, tr_x_xxxxy_xxxzz, tr_x_xxxxy_xxyy, tr_x_xxxxy_xxyyy, tr_x_xxxxy_xxyyz, tr_x_xxxxy_xxyz, tr_x_xxxxy_xxyzz, tr_x_xxxxy_xxzz, tr_x_xxxxy_xxzzz, tr_x_xxxxy_xyyy, tr_x_xxxxy_xyyyy, tr_x_xxxxy_xyyyz, tr_x_xxxxy_xyyz, tr_x_xxxxy_xyyzz, tr_x_xxxxy_xyzz, tr_x_xxxxy_xyzzz, tr_x_xxxxy_xzzz, tr_x_xxxxy_xzzzz, tr_x_xxxxy_zzzzz, tr_x_xxxxyy_xxxxx, tr_x_xxxxyy_xxxxy, tr_x_xxxxyy_xxxxz, tr_x_xxxxyy_xxxyy, tr_x_xxxxyy_xxxyz, tr_x_xxxxyy_xxxzz, tr_x_xxxxyy_xxyyy, tr_x_xxxxyy_xxyyz, tr_x_xxxxyy_xxyzz, tr_x_xxxxyy_xxzzz, tr_x_xxxxyy_xyyyy, tr_x_xxxxyy_xyyyz, tr_x_xxxxyy_xyyzz, tr_x_xxxxyy_xyzzz, tr_x_xxxxyy_xzzzz, tr_x_xxxxyy_yyyyy, tr_x_xxxxyy_yyyyz, tr_x_xxxxyy_yyyzz, tr_x_xxxxyy_yyzzz, tr_x_xxxxyy_yzzzz, tr_x_xxxxyy_zzzzz, tr_x_xxxyy_yyyyy, tr_x_xxxyy_yyyyz, tr_x_xxxyy_yyyzz, tr_x_xxxyy_yyzzz, tr_x_xxxyy_yzzzz, tr_x_xxyy_yyyyy, tr_x_xxyy_yyyyz, tr_x_xxyy_yyyzz, tr_x_xxyy_yyzzz, tr_x_xxyy_yzzzz, ts_xxxyy_yyyyy, ts_xxxyy_yyyyz, ts_xxxyy_yyyzz, ts_xxxyy_yyzzz, ts_xxxyy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyy_xxxxx[i] = tr_x_xxxx_xxxxx[i] * fe_0 + tr_x_xxxxy_xxxxx[i] * pa_y[i];

        tr_x_xxxxyy_xxxxy[i] = tr_x_xxxx_xxxxy[i] * fe_0 + tr_x_xxxxy_xxxx[i] * fe_0 + tr_x_xxxxy_xxxxy[i] * pa_y[i];

        tr_x_xxxxyy_xxxxz[i] = tr_x_xxxx_xxxxz[i] * fe_0 + tr_x_xxxxy_xxxxz[i] * pa_y[i];

        tr_x_xxxxyy_xxxyy[i] = tr_x_xxxx_xxxyy[i] * fe_0 + 2.0 * tr_x_xxxxy_xxxy[i] * fe_0 + tr_x_xxxxy_xxxyy[i] * pa_y[i];

        tr_x_xxxxyy_xxxyz[i] = tr_x_xxxx_xxxyz[i] * fe_0 + tr_x_xxxxy_xxxz[i] * fe_0 + tr_x_xxxxy_xxxyz[i] * pa_y[i];

        tr_x_xxxxyy_xxxzz[i] = tr_x_xxxx_xxxzz[i] * fe_0 + tr_x_xxxxy_xxxzz[i] * pa_y[i];

        tr_x_xxxxyy_xxyyy[i] = tr_x_xxxx_xxyyy[i] * fe_0 + 3.0 * tr_x_xxxxy_xxyy[i] * fe_0 + tr_x_xxxxy_xxyyy[i] * pa_y[i];

        tr_x_xxxxyy_xxyyz[i] = tr_x_xxxx_xxyyz[i] * fe_0 + 2.0 * tr_x_xxxxy_xxyz[i] * fe_0 + tr_x_xxxxy_xxyyz[i] * pa_y[i];

        tr_x_xxxxyy_xxyzz[i] = tr_x_xxxx_xxyzz[i] * fe_0 + tr_x_xxxxy_xxzz[i] * fe_0 + tr_x_xxxxy_xxyzz[i] * pa_y[i];

        tr_x_xxxxyy_xxzzz[i] = tr_x_xxxx_xxzzz[i] * fe_0 + tr_x_xxxxy_xxzzz[i] * pa_y[i];

        tr_x_xxxxyy_xyyyy[i] = tr_x_xxxx_xyyyy[i] * fe_0 + 4.0 * tr_x_xxxxy_xyyy[i] * fe_0 + tr_x_xxxxy_xyyyy[i] * pa_y[i];

        tr_x_xxxxyy_xyyyz[i] = tr_x_xxxx_xyyyz[i] * fe_0 + 3.0 * tr_x_xxxxy_xyyz[i] * fe_0 + tr_x_xxxxy_xyyyz[i] * pa_y[i];

        tr_x_xxxxyy_xyyzz[i] = tr_x_xxxx_xyyzz[i] * fe_0 + 2.0 * tr_x_xxxxy_xyzz[i] * fe_0 + tr_x_xxxxy_xyyzz[i] * pa_y[i];

        tr_x_xxxxyy_xyzzz[i] = tr_x_xxxx_xyzzz[i] * fe_0 + tr_x_xxxxy_xzzz[i] * fe_0 + tr_x_xxxxy_xyzzz[i] * pa_y[i];

        tr_x_xxxxyy_xzzzz[i] = tr_x_xxxx_xzzzz[i] * fe_0 + tr_x_xxxxy_xzzzz[i] * pa_y[i];

        tr_x_xxxxyy_yyyyy[i] = 3.0 * tr_x_xxyy_yyyyy[i] * fe_0 + ts_xxxyy_yyyyy[i] * fe_0 + tr_x_xxxyy_yyyyy[i] * pa_x[i];

        tr_x_xxxxyy_yyyyz[i] = 3.0 * tr_x_xxyy_yyyyz[i] * fe_0 + ts_xxxyy_yyyyz[i] * fe_0 + tr_x_xxxyy_yyyyz[i] * pa_x[i];

        tr_x_xxxxyy_yyyzz[i] = 3.0 * tr_x_xxyy_yyyzz[i] * fe_0 + ts_xxxyy_yyyzz[i] * fe_0 + tr_x_xxxyy_yyyzz[i] * pa_x[i];

        tr_x_xxxxyy_yyzzz[i] = 3.0 * tr_x_xxyy_yyzzz[i] * fe_0 + ts_xxxyy_yyzzz[i] * fe_0 + tr_x_xxxyy_yyzzz[i] * pa_x[i];

        tr_x_xxxxyy_yzzzz[i] = 3.0 * tr_x_xxyy_yzzzz[i] * fe_0 + ts_xxxyy_yzzzz[i] * fe_0 + tr_x_xxxyy_yzzzz[i] * pa_x[i];

        tr_x_xxxxyy_zzzzz[i] = tr_x_xxxx_zzzzz[i] * fe_0 + tr_x_xxxxy_zzzzz[i] * pa_y[i];
    }

    // Set up 84-105 components of targeted buffer : IH

    auto tr_x_xxxxyz_xxxxx = pbuffer.data(idx_dip_ih + 84);

    auto tr_x_xxxxyz_xxxxy = pbuffer.data(idx_dip_ih + 85);

    auto tr_x_xxxxyz_xxxxz = pbuffer.data(idx_dip_ih + 86);

    auto tr_x_xxxxyz_xxxyy = pbuffer.data(idx_dip_ih + 87);

    auto tr_x_xxxxyz_xxxyz = pbuffer.data(idx_dip_ih + 88);

    auto tr_x_xxxxyz_xxxzz = pbuffer.data(idx_dip_ih + 89);

    auto tr_x_xxxxyz_xxyyy = pbuffer.data(idx_dip_ih + 90);

    auto tr_x_xxxxyz_xxyyz = pbuffer.data(idx_dip_ih + 91);

    auto tr_x_xxxxyz_xxyzz = pbuffer.data(idx_dip_ih + 92);

    auto tr_x_xxxxyz_xxzzz = pbuffer.data(idx_dip_ih + 93);

    auto tr_x_xxxxyz_xyyyy = pbuffer.data(idx_dip_ih + 94);

    auto tr_x_xxxxyz_xyyyz = pbuffer.data(idx_dip_ih + 95);

    auto tr_x_xxxxyz_xyyzz = pbuffer.data(idx_dip_ih + 96);

    auto tr_x_xxxxyz_xyzzz = pbuffer.data(idx_dip_ih + 97);

    auto tr_x_xxxxyz_xzzzz = pbuffer.data(idx_dip_ih + 98);

    auto tr_x_xxxxyz_yyyyy = pbuffer.data(idx_dip_ih + 99);

    auto tr_x_xxxxyz_yyyyz = pbuffer.data(idx_dip_ih + 100);

    auto tr_x_xxxxyz_yyyzz = pbuffer.data(idx_dip_ih + 101);

    auto tr_x_xxxxyz_yyzzz = pbuffer.data(idx_dip_ih + 102);

    auto tr_x_xxxxyz_yzzzz = pbuffer.data(idx_dip_ih + 103);

    auto tr_x_xxxxyz_zzzzz = pbuffer.data(idx_dip_ih + 104);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxxy_xxxxy, tr_x_xxxxy_xxxyy, tr_x_xxxxy_xxyyy, tr_x_xxxxy_xyyyy, tr_x_xxxxy_yyyyy, tr_x_xxxxyz_xxxxx, tr_x_xxxxyz_xxxxy, tr_x_xxxxyz_xxxxz, tr_x_xxxxyz_xxxyy, tr_x_xxxxyz_xxxyz, tr_x_xxxxyz_xxxzz, tr_x_xxxxyz_xxyyy, tr_x_xxxxyz_xxyyz, tr_x_xxxxyz_xxyzz, tr_x_xxxxyz_xxzzz, tr_x_xxxxyz_xyyyy, tr_x_xxxxyz_xyyyz, tr_x_xxxxyz_xyyzz, tr_x_xxxxyz_xyzzz, tr_x_xxxxyz_xzzzz, tr_x_xxxxyz_yyyyy, tr_x_xxxxyz_yyyyz, tr_x_xxxxyz_yyyzz, tr_x_xxxxyz_yyzzz, tr_x_xxxxyz_yzzzz, tr_x_xxxxyz_zzzzz, tr_x_xxxxz_xxxxx, tr_x_xxxxz_xxxxz, tr_x_xxxxz_xxxyz, tr_x_xxxxz_xxxz, tr_x_xxxxz_xxxzz, tr_x_xxxxz_xxyyz, tr_x_xxxxz_xxyz, tr_x_xxxxz_xxyzz, tr_x_xxxxz_xxzz, tr_x_xxxxz_xxzzz, tr_x_xxxxz_xyyyz, tr_x_xxxxz_xyyz, tr_x_xxxxz_xyyzz, tr_x_xxxxz_xyzz, tr_x_xxxxz_xyzzz, tr_x_xxxxz_xzzz, tr_x_xxxxz_xzzzz, tr_x_xxxxz_yyyyz, tr_x_xxxxz_yyyz, tr_x_xxxxz_yyyzz, tr_x_xxxxz_yyzz, tr_x_xxxxz_yyzzz, tr_x_xxxxz_yzzz, tr_x_xxxxz_yzzzz, tr_x_xxxxz_zzzz, tr_x_xxxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyz_xxxxx[i] = tr_x_xxxxz_xxxxx[i] * pa_y[i];

        tr_x_xxxxyz_xxxxy[i] = tr_x_xxxxy_xxxxy[i] * pa_z[i];

        tr_x_xxxxyz_xxxxz[i] = tr_x_xxxxz_xxxxz[i] * pa_y[i];

        tr_x_xxxxyz_xxxyy[i] = tr_x_xxxxy_xxxyy[i] * pa_z[i];

        tr_x_xxxxyz_xxxyz[i] = tr_x_xxxxz_xxxz[i] * fe_0 + tr_x_xxxxz_xxxyz[i] * pa_y[i];

        tr_x_xxxxyz_xxxzz[i] = tr_x_xxxxz_xxxzz[i] * pa_y[i];

        tr_x_xxxxyz_xxyyy[i] = tr_x_xxxxy_xxyyy[i] * pa_z[i];

        tr_x_xxxxyz_xxyyz[i] = 2.0 * tr_x_xxxxz_xxyz[i] * fe_0 + tr_x_xxxxz_xxyyz[i] * pa_y[i];

        tr_x_xxxxyz_xxyzz[i] = tr_x_xxxxz_xxzz[i] * fe_0 + tr_x_xxxxz_xxyzz[i] * pa_y[i];

        tr_x_xxxxyz_xxzzz[i] = tr_x_xxxxz_xxzzz[i] * pa_y[i];

        tr_x_xxxxyz_xyyyy[i] = tr_x_xxxxy_xyyyy[i] * pa_z[i];

        tr_x_xxxxyz_xyyyz[i] = 3.0 * tr_x_xxxxz_xyyz[i] * fe_0 + tr_x_xxxxz_xyyyz[i] * pa_y[i];

        tr_x_xxxxyz_xyyzz[i] = 2.0 * tr_x_xxxxz_xyzz[i] * fe_0 + tr_x_xxxxz_xyyzz[i] * pa_y[i];

        tr_x_xxxxyz_xyzzz[i] = tr_x_xxxxz_xzzz[i] * fe_0 + tr_x_xxxxz_xyzzz[i] * pa_y[i];

        tr_x_xxxxyz_xzzzz[i] = tr_x_xxxxz_xzzzz[i] * pa_y[i];

        tr_x_xxxxyz_yyyyy[i] = tr_x_xxxxy_yyyyy[i] * pa_z[i];

        tr_x_xxxxyz_yyyyz[i] = 4.0 * tr_x_xxxxz_yyyz[i] * fe_0 + tr_x_xxxxz_yyyyz[i] * pa_y[i];

        tr_x_xxxxyz_yyyzz[i] = 3.0 * tr_x_xxxxz_yyzz[i] * fe_0 + tr_x_xxxxz_yyyzz[i] * pa_y[i];

        tr_x_xxxxyz_yyzzz[i] = 2.0 * tr_x_xxxxz_yzzz[i] * fe_0 + tr_x_xxxxz_yyzzz[i] * pa_y[i];

        tr_x_xxxxyz_yzzzz[i] = tr_x_xxxxz_zzzz[i] * fe_0 + tr_x_xxxxz_yzzzz[i] * pa_y[i];

        tr_x_xxxxyz_zzzzz[i] = tr_x_xxxxz_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : IH

    auto tr_x_xxxxzz_xxxxx = pbuffer.data(idx_dip_ih + 105);

    auto tr_x_xxxxzz_xxxxy = pbuffer.data(idx_dip_ih + 106);

    auto tr_x_xxxxzz_xxxxz = pbuffer.data(idx_dip_ih + 107);

    auto tr_x_xxxxzz_xxxyy = pbuffer.data(idx_dip_ih + 108);

    auto tr_x_xxxxzz_xxxyz = pbuffer.data(idx_dip_ih + 109);

    auto tr_x_xxxxzz_xxxzz = pbuffer.data(idx_dip_ih + 110);

    auto tr_x_xxxxzz_xxyyy = pbuffer.data(idx_dip_ih + 111);

    auto tr_x_xxxxzz_xxyyz = pbuffer.data(idx_dip_ih + 112);

    auto tr_x_xxxxzz_xxyzz = pbuffer.data(idx_dip_ih + 113);

    auto tr_x_xxxxzz_xxzzz = pbuffer.data(idx_dip_ih + 114);

    auto tr_x_xxxxzz_xyyyy = pbuffer.data(idx_dip_ih + 115);

    auto tr_x_xxxxzz_xyyyz = pbuffer.data(idx_dip_ih + 116);

    auto tr_x_xxxxzz_xyyzz = pbuffer.data(idx_dip_ih + 117);

    auto tr_x_xxxxzz_xyzzz = pbuffer.data(idx_dip_ih + 118);

    auto tr_x_xxxxzz_xzzzz = pbuffer.data(idx_dip_ih + 119);

    auto tr_x_xxxxzz_yyyyy = pbuffer.data(idx_dip_ih + 120);

    auto tr_x_xxxxzz_yyyyz = pbuffer.data(idx_dip_ih + 121);

    auto tr_x_xxxxzz_yyyzz = pbuffer.data(idx_dip_ih + 122);

    auto tr_x_xxxxzz_yyzzz = pbuffer.data(idx_dip_ih + 123);

    auto tr_x_xxxxzz_yzzzz = pbuffer.data(idx_dip_ih + 124);

    auto tr_x_xxxxzz_zzzzz = pbuffer.data(idx_dip_ih + 125);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxxx_xxxxx, tr_x_xxxx_xxxxy, tr_x_xxxx_xxxxz, tr_x_xxxx_xxxyy, tr_x_xxxx_xxxyz, tr_x_xxxx_xxxzz, tr_x_xxxx_xxyyy, tr_x_xxxx_xxyyz, tr_x_xxxx_xxyzz, tr_x_xxxx_xxzzz, tr_x_xxxx_xyyyy, tr_x_xxxx_xyyyz, tr_x_xxxx_xyyzz, tr_x_xxxx_xyzzz, tr_x_xxxx_xzzzz, tr_x_xxxx_yyyyy, tr_x_xxxxz_xxxx, tr_x_xxxxz_xxxxx, tr_x_xxxxz_xxxxy, tr_x_xxxxz_xxxxz, tr_x_xxxxz_xxxy, tr_x_xxxxz_xxxyy, tr_x_xxxxz_xxxyz, tr_x_xxxxz_xxxz, tr_x_xxxxz_xxxzz, tr_x_xxxxz_xxyy, tr_x_xxxxz_xxyyy, tr_x_xxxxz_xxyyz, tr_x_xxxxz_xxyz, tr_x_xxxxz_xxyzz, tr_x_xxxxz_xxzz, tr_x_xxxxz_xxzzz, tr_x_xxxxz_xyyy, tr_x_xxxxz_xyyyy, tr_x_xxxxz_xyyyz, tr_x_xxxxz_xyyz, tr_x_xxxxz_xyyzz, tr_x_xxxxz_xyzz, tr_x_xxxxz_xyzzz, tr_x_xxxxz_xzzz, tr_x_xxxxz_xzzzz, tr_x_xxxxz_yyyyy, tr_x_xxxxzz_xxxxx, tr_x_xxxxzz_xxxxy, tr_x_xxxxzz_xxxxz, tr_x_xxxxzz_xxxyy, tr_x_xxxxzz_xxxyz, tr_x_xxxxzz_xxxzz, tr_x_xxxxzz_xxyyy, tr_x_xxxxzz_xxyyz, tr_x_xxxxzz_xxyzz, tr_x_xxxxzz_xxzzz, tr_x_xxxxzz_xyyyy, tr_x_xxxxzz_xyyyz, tr_x_xxxxzz_xyyzz, tr_x_xxxxzz_xyzzz, tr_x_xxxxzz_xzzzz, tr_x_xxxxzz_yyyyy, tr_x_xxxxzz_yyyyz, tr_x_xxxxzz_yyyzz, tr_x_xxxxzz_yyzzz, tr_x_xxxxzz_yzzzz, tr_x_xxxxzz_zzzzz, tr_x_xxxzz_yyyyz, tr_x_xxxzz_yyyzz, tr_x_xxxzz_yyzzz, tr_x_xxxzz_yzzzz, tr_x_xxxzz_zzzzz, tr_x_xxzz_yyyyz, tr_x_xxzz_yyyzz, tr_x_xxzz_yyzzz, tr_x_xxzz_yzzzz, tr_x_xxzz_zzzzz, ts_xxxzz_yyyyz, ts_xxxzz_yyyzz, ts_xxxzz_yyzzz, ts_xxxzz_yzzzz, ts_xxxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxzz_xxxxx[i] = tr_x_xxxx_xxxxx[i] * fe_0 + tr_x_xxxxz_xxxxx[i] * pa_z[i];

        tr_x_xxxxzz_xxxxy[i] = tr_x_xxxx_xxxxy[i] * fe_0 + tr_x_xxxxz_xxxxy[i] * pa_z[i];

        tr_x_xxxxzz_xxxxz[i] = tr_x_xxxx_xxxxz[i] * fe_0 + tr_x_xxxxz_xxxx[i] * fe_0 + tr_x_xxxxz_xxxxz[i] * pa_z[i];

        tr_x_xxxxzz_xxxyy[i] = tr_x_xxxx_xxxyy[i] * fe_0 + tr_x_xxxxz_xxxyy[i] * pa_z[i];

        tr_x_xxxxzz_xxxyz[i] = tr_x_xxxx_xxxyz[i] * fe_0 + tr_x_xxxxz_xxxy[i] * fe_0 + tr_x_xxxxz_xxxyz[i] * pa_z[i];

        tr_x_xxxxzz_xxxzz[i] = tr_x_xxxx_xxxzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xxxz[i] * fe_0 + tr_x_xxxxz_xxxzz[i] * pa_z[i];

        tr_x_xxxxzz_xxyyy[i] = tr_x_xxxx_xxyyy[i] * fe_0 + tr_x_xxxxz_xxyyy[i] * pa_z[i];

        tr_x_xxxxzz_xxyyz[i] = tr_x_xxxx_xxyyz[i] * fe_0 + tr_x_xxxxz_xxyy[i] * fe_0 + tr_x_xxxxz_xxyyz[i] * pa_z[i];

        tr_x_xxxxzz_xxyzz[i] = tr_x_xxxx_xxyzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xxyz[i] * fe_0 + tr_x_xxxxz_xxyzz[i] * pa_z[i];

        tr_x_xxxxzz_xxzzz[i] = tr_x_xxxx_xxzzz[i] * fe_0 + 3.0 * tr_x_xxxxz_xxzz[i] * fe_0 + tr_x_xxxxz_xxzzz[i] * pa_z[i];

        tr_x_xxxxzz_xyyyy[i] = tr_x_xxxx_xyyyy[i] * fe_0 + tr_x_xxxxz_xyyyy[i] * pa_z[i];

        tr_x_xxxxzz_xyyyz[i] = tr_x_xxxx_xyyyz[i] * fe_0 + tr_x_xxxxz_xyyy[i] * fe_0 + tr_x_xxxxz_xyyyz[i] * pa_z[i];

        tr_x_xxxxzz_xyyzz[i] = tr_x_xxxx_xyyzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xyyz[i] * fe_0 + tr_x_xxxxz_xyyzz[i] * pa_z[i];

        tr_x_xxxxzz_xyzzz[i] = tr_x_xxxx_xyzzz[i] * fe_0 + 3.0 * tr_x_xxxxz_xyzz[i] * fe_0 + tr_x_xxxxz_xyzzz[i] * pa_z[i];

        tr_x_xxxxzz_xzzzz[i] = tr_x_xxxx_xzzzz[i] * fe_0 + 4.0 * tr_x_xxxxz_xzzz[i] * fe_0 + tr_x_xxxxz_xzzzz[i] * pa_z[i];

        tr_x_xxxxzz_yyyyy[i] = tr_x_xxxx_yyyyy[i] * fe_0 + tr_x_xxxxz_yyyyy[i] * pa_z[i];

        tr_x_xxxxzz_yyyyz[i] = 3.0 * tr_x_xxzz_yyyyz[i] * fe_0 + ts_xxxzz_yyyyz[i] * fe_0 + tr_x_xxxzz_yyyyz[i] * pa_x[i];

        tr_x_xxxxzz_yyyzz[i] = 3.0 * tr_x_xxzz_yyyzz[i] * fe_0 + ts_xxxzz_yyyzz[i] * fe_0 + tr_x_xxxzz_yyyzz[i] * pa_x[i];

        tr_x_xxxxzz_yyzzz[i] = 3.0 * tr_x_xxzz_yyzzz[i] * fe_0 + ts_xxxzz_yyzzz[i] * fe_0 + tr_x_xxxzz_yyzzz[i] * pa_x[i];

        tr_x_xxxxzz_yzzzz[i] = 3.0 * tr_x_xxzz_yzzzz[i] * fe_0 + ts_xxxzz_yzzzz[i] * fe_0 + tr_x_xxxzz_yzzzz[i] * pa_x[i];

        tr_x_xxxxzz_zzzzz[i] = 3.0 * tr_x_xxzz_zzzzz[i] * fe_0 + ts_xxxzz_zzzzz[i] * fe_0 + tr_x_xxxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : IH

    auto tr_x_xxxyyy_xxxxx = pbuffer.data(idx_dip_ih + 126);

    auto tr_x_xxxyyy_xxxxy = pbuffer.data(idx_dip_ih + 127);

    auto tr_x_xxxyyy_xxxxz = pbuffer.data(idx_dip_ih + 128);

    auto tr_x_xxxyyy_xxxyy = pbuffer.data(idx_dip_ih + 129);

    auto tr_x_xxxyyy_xxxyz = pbuffer.data(idx_dip_ih + 130);

    auto tr_x_xxxyyy_xxxzz = pbuffer.data(idx_dip_ih + 131);

    auto tr_x_xxxyyy_xxyyy = pbuffer.data(idx_dip_ih + 132);

    auto tr_x_xxxyyy_xxyyz = pbuffer.data(idx_dip_ih + 133);

    auto tr_x_xxxyyy_xxyzz = pbuffer.data(idx_dip_ih + 134);

    auto tr_x_xxxyyy_xxzzz = pbuffer.data(idx_dip_ih + 135);

    auto tr_x_xxxyyy_xyyyy = pbuffer.data(idx_dip_ih + 136);

    auto tr_x_xxxyyy_xyyyz = pbuffer.data(idx_dip_ih + 137);

    auto tr_x_xxxyyy_xyyzz = pbuffer.data(idx_dip_ih + 138);

    auto tr_x_xxxyyy_xyzzz = pbuffer.data(idx_dip_ih + 139);

    auto tr_x_xxxyyy_xzzzz = pbuffer.data(idx_dip_ih + 140);

    auto tr_x_xxxyyy_yyyyy = pbuffer.data(idx_dip_ih + 141);

    auto tr_x_xxxyyy_yyyyz = pbuffer.data(idx_dip_ih + 142);

    auto tr_x_xxxyyy_yyyzz = pbuffer.data(idx_dip_ih + 143);

    auto tr_x_xxxyyy_yyzzz = pbuffer.data(idx_dip_ih + 144);

    auto tr_x_xxxyyy_yzzzz = pbuffer.data(idx_dip_ih + 145);

    auto tr_x_xxxyyy_zzzzz = pbuffer.data(idx_dip_ih + 146);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxxy_xxxxx, tr_x_xxxy_xxxxy, tr_x_xxxy_xxxxz, tr_x_xxxy_xxxyy, tr_x_xxxy_xxxyz, tr_x_xxxy_xxxzz, tr_x_xxxy_xxyyy, tr_x_xxxy_xxyyz, tr_x_xxxy_xxyzz, tr_x_xxxy_xxzzz, tr_x_xxxy_xyyyy, tr_x_xxxy_xyyyz, tr_x_xxxy_xyyzz, tr_x_xxxy_xyzzz, tr_x_xxxy_xzzzz, tr_x_xxxy_zzzzz, tr_x_xxxyy_xxxx, tr_x_xxxyy_xxxxx, tr_x_xxxyy_xxxxy, tr_x_xxxyy_xxxxz, tr_x_xxxyy_xxxy, tr_x_xxxyy_xxxyy, tr_x_xxxyy_xxxyz, tr_x_xxxyy_xxxz, tr_x_xxxyy_xxxzz, tr_x_xxxyy_xxyy, tr_x_xxxyy_xxyyy, tr_x_xxxyy_xxyyz, tr_x_xxxyy_xxyz, tr_x_xxxyy_xxyzz, tr_x_xxxyy_xxzz, tr_x_xxxyy_xxzzz, tr_x_xxxyy_xyyy, tr_x_xxxyy_xyyyy, tr_x_xxxyy_xyyyz, tr_x_xxxyy_xyyz, tr_x_xxxyy_xyyzz, tr_x_xxxyy_xyzz, tr_x_xxxyy_xyzzz, tr_x_xxxyy_xzzz, tr_x_xxxyy_xzzzz, tr_x_xxxyy_zzzzz, tr_x_xxxyyy_xxxxx, tr_x_xxxyyy_xxxxy, tr_x_xxxyyy_xxxxz, tr_x_xxxyyy_xxxyy, tr_x_xxxyyy_xxxyz, tr_x_xxxyyy_xxxzz, tr_x_xxxyyy_xxyyy, tr_x_xxxyyy_xxyyz, tr_x_xxxyyy_xxyzz, tr_x_xxxyyy_xxzzz, tr_x_xxxyyy_xyyyy, tr_x_xxxyyy_xyyyz, tr_x_xxxyyy_xyyzz, tr_x_xxxyyy_xyzzz, tr_x_xxxyyy_xzzzz, tr_x_xxxyyy_yyyyy, tr_x_xxxyyy_yyyyz, tr_x_xxxyyy_yyyzz, tr_x_xxxyyy_yyzzz, tr_x_xxxyyy_yzzzz, tr_x_xxxyyy_zzzzz, tr_x_xxyyy_yyyyy, tr_x_xxyyy_yyyyz, tr_x_xxyyy_yyyzz, tr_x_xxyyy_yyzzz, tr_x_xxyyy_yzzzz, tr_x_xyyy_yyyyy, tr_x_xyyy_yyyyz, tr_x_xyyy_yyyzz, tr_x_xyyy_yyzzz, tr_x_xyyy_yzzzz, ts_xxyyy_yyyyy, ts_xxyyy_yyyyz, ts_xxyyy_yyyzz, ts_xxyyy_yyzzz, ts_xxyyy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyy_xxxxx[i] = 2.0 * tr_x_xxxy_xxxxx[i] * fe_0 + tr_x_xxxyy_xxxxx[i] * pa_y[i];

        tr_x_xxxyyy_xxxxy[i] = 2.0 * tr_x_xxxy_xxxxy[i] * fe_0 + tr_x_xxxyy_xxxx[i] * fe_0 + tr_x_xxxyy_xxxxy[i] * pa_y[i];

        tr_x_xxxyyy_xxxxz[i] = 2.0 * tr_x_xxxy_xxxxz[i] * fe_0 + tr_x_xxxyy_xxxxz[i] * pa_y[i];

        tr_x_xxxyyy_xxxyy[i] = 2.0 * tr_x_xxxy_xxxyy[i] * fe_0 + 2.0 * tr_x_xxxyy_xxxy[i] * fe_0 + tr_x_xxxyy_xxxyy[i] * pa_y[i];

        tr_x_xxxyyy_xxxyz[i] = 2.0 * tr_x_xxxy_xxxyz[i] * fe_0 + tr_x_xxxyy_xxxz[i] * fe_0 + tr_x_xxxyy_xxxyz[i] * pa_y[i];

        tr_x_xxxyyy_xxxzz[i] = 2.0 * tr_x_xxxy_xxxzz[i] * fe_0 + tr_x_xxxyy_xxxzz[i] * pa_y[i];

        tr_x_xxxyyy_xxyyy[i] = 2.0 * tr_x_xxxy_xxyyy[i] * fe_0 + 3.0 * tr_x_xxxyy_xxyy[i] * fe_0 + tr_x_xxxyy_xxyyy[i] * pa_y[i];

        tr_x_xxxyyy_xxyyz[i] = 2.0 * tr_x_xxxy_xxyyz[i] * fe_0 + 2.0 * tr_x_xxxyy_xxyz[i] * fe_0 + tr_x_xxxyy_xxyyz[i] * pa_y[i];

        tr_x_xxxyyy_xxyzz[i] = 2.0 * tr_x_xxxy_xxyzz[i] * fe_0 + tr_x_xxxyy_xxzz[i] * fe_0 + tr_x_xxxyy_xxyzz[i] * pa_y[i];

        tr_x_xxxyyy_xxzzz[i] = 2.0 * tr_x_xxxy_xxzzz[i] * fe_0 + tr_x_xxxyy_xxzzz[i] * pa_y[i];

        tr_x_xxxyyy_xyyyy[i] = 2.0 * tr_x_xxxy_xyyyy[i] * fe_0 + 4.0 * tr_x_xxxyy_xyyy[i] * fe_0 + tr_x_xxxyy_xyyyy[i] * pa_y[i];

        tr_x_xxxyyy_xyyyz[i] = 2.0 * tr_x_xxxy_xyyyz[i] * fe_0 + 3.0 * tr_x_xxxyy_xyyz[i] * fe_0 + tr_x_xxxyy_xyyyz[i] * pa_y[i];

        tr_x_xxxyyy_xyyzz[i] = 2.0 * tr_x_xxxy_xyyzz[i] * fe_0 + 2.0 * tr_x_xxxyy_xyzz[i] * fe_0 + tr_x_xxxyy_xyyzz[i] * pa_y[i];

        tr_x_xxxyyy_xyzzz[i] = 2.0 * tr_x_xxxy_xyzzz[i] * fe_0 + tr_x_xxxyy_xzzz[i] * fe_0 + tr_x_xxxyy_xyzzz[i] * pa_y[i];

        tr_x_xxxyyy_xzzzz[i] = 2.0 * tr_x_xxxy_xzzzz[i] * fe_0 + tr_x_xxxyy_xzzzz[i] * pa_y[i];

        tr_x_xxxyyy_yyyyy[i] = 2.0 * tr_x_xyyy_yyyyy[i] * fe_0 + ts_xxyyy_yyyyy[i] * fe_0 + tr_x_xxyyy_yyyyy[i] * pa_x[i];

        tr_x_xxxyyy_yyyyz[i] = 2.0 * tr_x_xyyy_yyyyz[i] * fe_0 + ts_xxyyy_yyyyz[i] * fe_0 + tr_x_xxyyy_yyyyz[i] * pa_x[i];

        tr_x_xxxyyy_yyyzz[i] = 2.0 * tr_x_xyyy_yyyzz[i] * fe_0 + ts_xxyyy_yyyzz[i] * fe_0 + tr_x_xxyyy_yyyzz[i] * pa_x[i];

        tr_x_xxxyyy_yyzzz[i] = 2.0 * tr_x_xyyy_yyzzz[i] * fe_0 + ts_xxyyy_yyzzz[i] * fe_0 + tr_x_xxyyy_yyzzz[i] * pa_x[i];

        tr_x_xxxyyy_yzzzz[i] = 2.0 * tr_x_xyyy_yzzzz[i] * fe_0 + ts_xxyyy_yzzzz[i] * fe_0 + tr_x_xxyyy_yzzzz[i] * pa_x[i];

        tr_x_xxxyyy_zzzzz[i] = 2.0 * tr_x_xxxy_zzzzz[i] * fe_0 + tr_x_xxxyy_zzzzz[i] * pa_y[i];
    }

    // Set up 147-168 components of targeted buffer : IH

    auto tr_x_xxxyyz_xxxxx = pbuffer.data(idx_dip_ih + 147);

    auto tr_x_xxxyyz_xxxxy = pbuffer.data(idx_dip_ih + 148);

    auto tr_x_xxxyyz_xxxxz = pbuffer.data(idx_dip_ih + 149);

    auto tr_x_xxxyyz_xxxyy = pbuffer.data(idx_dip_ih + 150);

    auto tr_x_xxxyyz_xxxyz = pbuffer.data(idx_dip_ih + 151);

    auto tr_x_xxxyyz_xxxzz = pbuffer.data(idx_dip_ih + 152);

    auto tr_x_xxxyyz_xxyyy = pbuffer.data(idx_dip_ih + 153);

    auto tr_x_xxxyyz_xxyyz = pbuffer.data(idx_dip_ih + 154);

    auto tr_x_xxxyyz_xxyzz = pbuffer.data(idx_dip_ih + 155);

    auto tr_x_xxxyyz_xxzzz = pbuffer.data(idx_dip_ih + 156);

    auto tr_x_xxxyyz_xyyyy = pbuffer.data(idx_dip_ih + 157);

    auto tr_x_xxxyyz_xyyyz = pbuffer.data(idx_dip_ih + 158);

    auto tr_x_xxxyyz_xyyzz = pbuffer.data(idx_dip_ih + 159);

    auto tr_x_xxxyyz_xyzzz = pbuffer.data(idx_dip_ih + 160);

    auto tr_x_xxxyyz_xzzzz = pbuffer.data(idx_dip_ih + 161);

    auto tr_x_xxxyyz_yyyyy = pbuffer.data(idx_dip_ih + 162);

    auto tr_x_xxxyyz_yyyyz = pbuffer.data(idx_dip_ih + 163);

    auto tr_x_xxxyyz_yyyzz = pbuffer.data(idx_dip_ih + 164);

    auto tr_x_xxxyyz_yyzzz = pbuffer.data(idx_dip_ih + 165);

    auto tr_x_xxxyyz_yzzzz = pbuffer.data(idx_dip_ih + 166);

    auto tr_x_xxxyyz_zzzzz = pbuffer.data(idx_dip_ih + 167);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxyy_xxxxx, tr_x_xxxyy_xxxxy, tr_x_xxxyy_xxxy, tr_x_xxxyy_xxxyy, tr_x_xxxyy_xxxyz, tr_x_xxxyy_xxyy, tr_x_xxxyy_xxyyy, tr_x_xxxyy_xxyyz, tr_x_xxxyy_xxyz, tr_x_xxxyy_xxyzz, tr_x_xxxyy_xyyy, tr_x_xxxyy_xyyyy, tr_x_xxxyy_xyyyz, tr_x_xxxyy_xyyz, tr_x_xxxyy_xyyzz, tr_x_xxxyy_xyzz, tr_x_xxxyy_xyzzz, tr_x_xxxyy_yyyy, tr_x_xxxyy_yyyyy, tr_x_xxxyy_yyyyz, tr_x_xxxyy_yyyz, tr_x_xxxyy_yyyzz, tr_x_xxxyy_yyzz, tr_x_xxxyy_yyzzz, tr_x_xxxyy_yzzz, tr_x_xxxyy_yzzzz, tr_x_xxxyyz_xxxxx, tr_x_xxxyyz_xxxxy, tr_x_xxxyyz_xxxxz, tr_x_xxxyyz_xxxyy, tr_x_xxxyyz_xxxyz, tr_x_xxxyyz_xxxzz, tr_x_xxxyyz_xxyyy, tr_x_xxxyyz_xxyyz, tr_x_xxxyyz_xxyzz, tr_x_xxxyyz_xxzzz, tr_x_xxxyyz_xyyyy, tr_x_xxxyyz_xyyyz, tr_x_xxxyyz_xyyzz, tr_x_xxxyyz_xyzzz, tr_x_xxxyyz_xzzzz, tr_x_xxxyyz_yyyyy, tr_x_xxxyyz_yyyyz, tr_x_xxxyyz_yyyzz, tr_x_xxxyyz_yyzzz, tr_x_xxxyyz_yzzzz, tr_x_xxxyyz_zzzzz, tr_x_xxxyz_xxxxz, tr_x_xxxyz_xxxzz, tr_x_xxxyz_xxzzz, tr_x_xxxyz_xzzzz, tr_x_xxxyz_zzzzz, tr_x_xxxz_xxxxz, tr_x_xxxz_xxxzz, tr_x_xxxz_xxzzz, tr_x_xxxz_xzzzz, tr_x_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyz_xxxxx[i] = tr_x_xxxyy_xxxxx[i] * pa_z[i];

        tr_x_xxxyyz_xxxxy[i] = tr_x_xxxyy_xxxxy[i] * pa_z[i];

        tr_x_xxxyyz_xxxxz[i] = tr_x_xxxz_xxxxz[i] * fe_0 + tr_x_xxxyz_xxxxz[i] * pa_y[i];

        tr_x_xxxyyz_xxxyy[i] = tr_x_xxxyy_xxxyy[i] * pa_z[i];

        tr_x_xxxyyz_xxxyz[i] = tr_x_xxxyy_xxxy[i] * fe_0 + tr_x_xxxyy_xxxyz[i] * pa_z[i];

        tr_x_xxxyyz_xxxzz[i] = tr_x_xxxz_xxxzz[i] * fe_0 + tr_x_xxxyz_xxxzz[i] * pa_y[i];

        tr_x_xxxyyz_xxyyy[i] = tr_x_xxxyy_xxyyy[i] * pa_z[i];

        tr_x_xxxyyz_xxyyz[i] = tr_x_xxxyy_xxyy[i] * fe_0 + tr_x_xxxyy_xxyyz[i] * pa_z[i];

        tr_x_xxxyyz_xxyzz[i] = 2.0 * tr_x_xxxyy_xxyz[i] * fe_0 + tr_x_xxxyy_xxyzz[i] * pa_z[i];

        tr_x_xxxyyz_xxzzz[i] = tr_x_xxxz_xxzzz[i] * fe_0 + tr_x_xxxyz_xxzzz[i] * pa_y[i];

        tr_x_xxxyyz_xyyyy[i] = tr_x_xxxyy_xyyyy[i] * pa_z[i];

        tr_x_xxxyyz_xyyyz[i] = tr_x_xxxyy_xyyy[i] * fe_0 + tr_x_xxxyy_xyyyz[i] * pa_z[i];

        tr_x_xxxyyz_xyyzz[i] = 2.0 * tr_x_xxxyy_xyyz[i] * fe_0 + tr_x_xxxyy_xyyzz[i] * pa_z[i];

        tr_x_xxxyyz_xyzzz[i] = 3.0 * tr_x_xxxyy_xyzz[i] * fe_0 + tr_x_xxxyy_xyzzz[i] * pa_z[i];

        tr_x_xxxyyz_xzzzz[i] = tr_x_xxxz_xzzzz[i] * fe_0 + tr_x_xxxyz_xzzzz[i] * pa_y[i];

        tr_x_xxxyyz_yyyyy[i] = tr_x_xxxyy_yyyyy[i] * pa_z[i];

        tr_x_xxxyyz_yyyyz[i] = tr_x_xxxyy_yyyy[i] * fe_0 + tr_x_xxxyy_yyyyz[i] * pa_z[i];

        tr_x_xxxyyz_yyyzz[i] = 2.0 * tr_x_xxxyy_yyyz[i] * fe_0 + tr_x_xxxyy_yyyzz[i] * pa_z[i];

        tr_x_xxxyyz_yyzzz[i] = 3.0 * tr_x_xxxyy_yyzz[i] * fe_0 + tr_x_xxxyy_yyzzz[i] * pa_z[i];

        tr_x_xxxyyz_yzzzz[i] = 4.0 * tr_x_xxxyy_yzzz[i] * fe_0 + tr_x_xxxyy_yzzzz[i] * pa_z[i];

        tr_x_xxxyyz_zzzzz[i] = tr_x_xxxz_zzzzz[i] * fe_0 + tr_x_xxxyz_zzzzz[i] * pa_y[i];
    }

    // Set up 168-189 components of targeted buffer : IH

    auto tr_x_xxxyzz_xxxxx = pbuffer.data(idx_dip_ih + 168);

    auto tr_x_xxxyzz_xxxxy = pbuffer.data(idx_dip_ih + 169);

    auto tr_x_xxxyzz_xxxxz = pbuffer.data(idx_dip_ih + 170);

    auto tr_x_xxxyzz_xxxyy = pbuffer.data(idx_dip_ih + 171);

    auto tr_x_xxxyzz_xxxyz = pbuffer.data(idx_dip_ih + 172);

    auto tr_x_xxxyzz_xxxzz = pbuffer.data(idx_dip_ih + 173);

    auto tr_x_xxxyzz_xxyyy = pbuffer.data(idx_dip_ih + 174);

    auto tr_x_xxxyzz_xxyyz = pbuffer.data(idx_dip_ih + 175);

    auto tr_x_xxxyzz_xxyzz = pbuffer.data(idx_dip_ih + 176);

    auto tr_x_xxxyzz_xxzzz = pbuffer.data(idx_dip_ih + 177);

    auto tr_x_xxxyzz_xyyyy = pbuffer.data(idx_dip_ih + 178);

    auto tr_x_xxxyzz_xyyyz = pbuffer.data(idx_dip_ih + 179);

    auto tr_x_xxxyzz_xyyzz = pbuffer.data(idx_dip_ih + 180);

    auto tr_x_xxxyzz_xyzzz = pbuffer.data(idx_dip_ih + 181);

    auto tr_x_xxxyzz_xzzzz = pbuffer.data(idx_dip_ih + 182);

    auto tr_x_xxxyzz_yyyyy = pbuffer.data(idx_dip_ih + 183);

    auto tr_x_xxxyzz_yyyyz = pbuffer.data(idx_dip_ih + 184);

    auto tr_x_xxxyzz_yyyzz = pbuffer.data(idx_dip_ih + 185);

    auto tr_x_xxxyzz_yyzzz = pbuffer.data(idx_dip_ih + 186);

    auto tr_x_xxxyzz_yzzzz = pbuffer.data(idx_dip_ih + 187);

    auto tr_x_xxxyzz_zzzzz = pbuffer.data(idx_dip_ih + 188);

    #pragma omp simd aligned(pa_y, tr_x_xxxyzz_xxxxx, tr_x_xxxyzz_xxxxy, tr_x_xxxyzz_xxxxz, tr_x_xxxyzz_xxxyy, tr_x_xxxyzz_xxxyz, tr_x_xxxyzz_xxxzz, tr_x_xxxyzz_xxyyy, tr_x_xxxyzz_xxyyz, tr_x_xxxyzz_xxyzz, tr_x_xxxyzz_xxzzz, tr_x_xxxyzz_xyyyy, tr_x_xxxyzz_xyyyz, tr_x_xxxyzz_xyyzz, tr_x_xxxyzz_xyzzz, tr_x_xxxyzz_xzzzz, tr_x_xxxyzz_yyyyy, tr_x_xxxyzz_yyyyz, tr_x_xxxyzz_yyyzz, tr_x_xxxyzz_yyzzz, tr_x_xxxyzz_yzzzz, tr_x_xxxyzz_zzzzz, tr_x_xxxzz_xxxx, tr_x_xxxzz_xxxxx, tr_x_xxxzz_xxxxy, tr_x_xxxzz_xxxxz, tr_x_xxxzz_xxxy, tr_x_xxxzz_xxxyy, tr_x_xxxzz_xxxyz, tr_x_xxxzz_xxxz, tr_x_xxxzz_xxxzz, tr_x_xxxzz_xxyy, tr_x_xxxzz_xxyyy, tr_x_xxxzz_xxyyz, tr_x_xxxzz_xxyz, tr_x_xxxzz_xxyzz, tr_x_xxxzz_xxzz, tr_x_xxxzz_xxzzz, tr_x_xxxzz_xyyy, tr_x_xxxzz_xyyyy, tr_x_xxxzz_xyyyz, tr_x_xxxzz_xyyz, tr_x_xxxzz_xyyzz, tr_x_xxxzz_xyzz, tr_x_xxxzz_xyzzz, tr_x_xxxzz_xzzz, tr_x_xxxzz_xzzzz, tr_x_xxxzz_yyyy, tr_x_xxxzz_yyyyy, tr_x_xxxzz_yyyyz, tr_x_xxxzz_yyyz, tr_x_xxxzz_yyyzz, tr_x_xxxzz_yyzz, tr_x_xxxzz_yyzzz, tr_x_xxxzz_yzzz, tr_x_xxxzz_yzzzz, tr_x_xxxzz_zzzz, tr_x_xxxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyzz_xxxxx[i] = tr_x_xxxzz_xxxxx[i] * pa_y[i];

        tr_x_xxxyzz_xxxxy[i] = tr_x_xxxzz_xxxx[i] * fe_0 + tr_x_xxxzz_xxxxy[i] * pa_y[i];

        tr_x_xxxyzz_xxxxz[i] = tr_x_xxxzz_xxxxz[i] * pa_y[i];

        tr_x_xxxyzz_xxxyy[i] = 2.0 * tr_x_xxxzz_xxxy[i] * fe_0 + tr_x_xxxzz_xxxyy[i] * pa_y[i];

        tr_x_xxxyzz_xxxyz[i] = tr_x_xxxzz_xxxz[i] * fe_0 + tr_x_xxxzz_xxxyz[i] * pa_y[i];

        tr_x_xxxyzz_xxxzz[i] = tr_x_xxxzz_xxxzz[i] * pa_y[i];

        tr_x_xxxyzz_xxyyy[i] = 3.0 * tr_x_xxxzz_xxyy[i] * fe_0 + tr_x_xxxzz_xxyyy[i] * pa_y[i];

        tr_x_xxxyzz_xxyyz[i] = 2.0 * tr_x_xxxzz_xxyz[i] * fe_0 + tr_x_xxxzz_xxyyz[i] * pa_y[i];

        tr_x_xxxyzz_xxyzz[i] = tr_x_xxxzz_xxzz[i] * fe_0 + tr_x_xxxzz_xxyzz[i] * pa_y[i];

        tr_x_xxxyzz_xxzzz[i] = tr_x_xxxzz_xxzzz[i] * pa_y[i];

        tr_x_xxxyzz_xyyyy[i] = 4.0 * tr_x_xxxzz_xyyy[i] * fe_0 + tr_x_xxxzz_xyyyy[i] * pa_y[i];

        tr_x_xxxyzz_xyyyz[i] = 3.0 * tr_x_xxxzz_xyyz[i] * fe_0 + tr_x_xxxzz_xyyyz[i] * pa_y[i];

        tr_x_xxxyzz_xyyzz[i] = 2.0 * tr_x_xxxzz_xyzz[i] * fe_0 + tr_x_xxxzz_xyyzz[i] * pa_y[i];

        tr_x_xxxyzz_xyzzz[i] = tr_x_xxxzz_xzzz[i] * fe_0 + tr_x_xxxzz_xyzzz[i] * pa_y[i];

        tr_x_xxxyzz_xzzzz[i] = tr_x_xxxzz_xzzzz[i] * pa_y[i];

        tr_x_xxxyzz_yyyyy[i] = 5.0 * tr_x_xxxzz_yyyy[i] * fe_0 + tr_x_xxxzz_yyyyy[i] * pa_y[i];

        tr_x_xxxyzz_yyyyz[i] = 4.0 * tr_x_xxxzz_yyyz[i] * fe_0 + tr_x_xxxzz_yyyyz[i] * pa_y[i];

        tr_x_xxxyzz_yyyzz[i] = 3.0 * tr_x_xxxzz_yyzz[i] * fe_0 + tr_x_xxxzz_yyyzz[i] * pa_y[i];

        tr_x_xxxyzz_yyzzz[i] = 2.0 * tr_x_xxxzz_yzzz[i] * fe_0 + tr_x_xxxzz_yyzzz[i] * pa_y[i];

        tr_x_xxxyzz_yzzzz[i] = tr_x_xxxzz_zzzz[i] * fe_0 + tr_x_xxxzz_yzzzz[i] * pa_y[i];

        tr_x_xxxyzz_zzzzz[i] = tr_x_xxxzz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : IH

    auto tr_x_xxxzzz_xxxxx = pbuffer.data(idx_dip_ih + 189);

    auto tr_x_xxxzzz_xxxxy = pbuffer.data(idx_dip_ih + 190);

    auto tr_x_xxxzzz_xxxxz = pbuffer.data(idx_dip_ih + 191);

    auto tr_x_xxxzzz_xxxyy = pbuffer.data(idx_dip_ih + 192);

    auto tr_x_xxxzzz_xxxyz = pbuffer.data(idx_dip_ih + 193);

    auto tr_x_xxxzzz_xxxzz = pbuffer.data(idx_dip_ih + 194);

    auto tr_x_xxxzzz_xxyyy = pbuffer.data(idx_dip_ih + 195);

    auto tr_x_xxxzzz_xxyyz = pbuffer.data(idx_dip_ih + 196);

    auto tr_x_xxxzzz_xxyzz = pbuffer.data(idx_dip_ih + 197);

    auto tr_x_xxxzzz_xxzzz = pbuffer.data(idx_dip_ih + 198);

    auto tr_x_xxxzzz_xyyyy = pbuffer.data(idx_dip_ih + 199);

    auto tr_x_xxxzzz_xyyyz = pbuffer.data(idx_dip_ih + 200);

    auto tr_x_xxxzzz_xyyzz = pbuffer.data(idx_dip_ih + 201);

    auto tr_x_xxxzzz_xyzzz = pbuffer.data(idx_dip_ih + 202);

    auto tr_x_xxxzzz_xzzzz = pbuffer.data(idx_dip_ih + 203);

    auto tr_x_xxxzzz_yyyyy = pbuffer.data(idx_dip_ih + 204);

    auto tr_x_xxxzzz_yyyyz = pbuffer.data(idx_dip_ih + 205);

    auto tr_x_xxxzzz_yyyzz = pbuffer.data(idx_dip_ih + 206);

    auto tr_x_xxxzzz_yyzzz = pbuffer.data(idx_dip_ih + 207);

    auto tr_x_xxxzzz_yzzzz = pbuffer.data(idx_dip_ih + 208);

    auto tr_x_xxxzzz_zzzzz = pbuffer.data(idx_dip_ih + 209);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxxz_xxxxx, tr_x_xxxz_xxxxy, tr_x_xxxz_xxxxz, tr_x_xxxz_xxxyy, tr_x_xxxz_xxxyz, tr_x_xxxz_xxxzz, tr_x_xxxz_xxyyy, tr_x_xxxz_xxyyz, tr_x_xxxz_xxyzz, tr_x_xxxz_xxzzz, tr_x_xxxz_xyyyy, tr_x_xxxz_xyyyz, tr_x_xxxz_xyyzz, tr_x_xxxz_xyzzz, tr_x_xxxz_xzzzz, tr_x_xxxz_yyyyy, tr_x_xxxzz_xxxx, tr_x_xxxzz_xxxxx, tr_x_xxxzz_xxxxy, tr_x_xxxzz_xxxxz, tr_x_xxxzz_xxxy, tr_x_xxxzz_xxxyy, tr_x_xxxzz_xxxyz, tr_x_xxxzz_xxxz, tr_x_xxxzz_xxxzz, tr_x_xxxzz_xxyy, tr_x_xxxzz_xxyyy, tr_x_xxxzz_xxyyz, tr_x_xxxzz_xxyz, tr_x_xxxzz_xxyzz, tr_x_xxxzz_xxzz, tr_x_xxxzz_xxzzz, tr_x_xxxzz_xyyy, tr_x_xxxzz_xyyyy, tr_x_xxxzz_xyyyz, tr_x_xxxzz_xyyz, tr_x_xxxzz_xyyzz, tr_x_xxxzz_xyzz, tr_x_xxxzz_xyzzz, tr_x_xxxzz_xzzz, tr_x_xxxzz_xzzzz, tr_x_xxxzz_yyyyy, tr_x_xxxzzz_xxxxx, tr_x_xxxzzz_xxxxy, tr_x_xxxzzz_xxxxz, tr_x_xxxzzz_xxxyy, tr_x_xxxzzz_xxxyz, tr_x_xxxzzz_xxxzz, tr_x_xxxzzz_xxyyy, tr_x_xxxzzz_xxyyz, tr_x_xxxzzz_xxyzz, tr_x_xxxzzz_xxzzz, tr_x_xxxzzz_xyyyy, tr_x_xxxzzz_xyyyz, tr_x_xxxzzz_xyyzz, tr_x_xxxzzz_xyzzz, tr_x_xxxzzz_xzzzz, tr_x_xxxzzz_yyyyy, tr_x_xxxzzz_yyyyz, tr_x_xxxzzz_yyyzz, tr_x_xxxzzz_yyzzz, tr_x_xxxzzz_yzzzz, tr_x_xxxzzz_zzzzz, tr_x_xxzzz_yyyyz, tr_x_xxzzz_yyyzz, tr_x_xxzzz_yyzzz, tr_x_xxzzz_yzzzz, tr_x_xxzzz_zzzzz, tr_x_xzzz_yyyyz, tr_x_xzzz_yyyzz, tr_x_xzzz_yyzzz, tr_x_xzzz_yzzzz, tr_x_xzzz_zzzzz, ts_xxzzz_yyyyz, ts_xxzzz_yyyzz, ts_xxzzz_yyzzz, ts_xxzzz_yzzzz, ts_xxzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzzz_xxxxx[i] = 2.0 * tr_x_xxxz_xxxxx[i] * fe_0 + tr_x_xxxzz_xxxxx[i] * pa_z[i];

        tr_x_xxxzzz_xxxxy[i] = 2.0 * tr_x_xxxz_xxxxy[i] * fe_0 + tr_x_xxxzz_xxxxy[i] * pa_z[i];

        tr_x_xxxzzz_xxxxz[i] = 2.0 * tr_x_xxxz_xxxxz[i] * fe_0 + tr_x_xxxzz_xxxx[i] * fe_0 + tr_x_xxxzz_xxxxz[i] * pa_z[i];

        tr_x_xxxzzz_xxxyy[i] = 2.0 * tr_x_xxxz_xxxyy[i] * fe_0 + tr_x_xxxzz_xxxyy[i] * pa_z[i];

        tr_x_xxxzzz_xxxyz[i] = 2.0 * tr_x_xxxz_xxxyz[i] * fe_0 + tr_x_xxxzz_xxxy[i] * fe_0 + tr_x_xxxzz_xxxyz[i] * pa_z[i];

        tr_x_xxxzzz_xxxzz[i] = 2.0 * tr_x_xxxz_xxxzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xxxz[i] * fe_0 + tr_x_xxxzz_xxxzz[i] * pa_z[i];

        tr_x_xxxzzz_xxyyy[i] = 2.0 * tr_x_xxxz_xxyyy[i] * fe_0 + tr_x_xxxzz_xxyyy[i] * pa_z[i];

        tr_x_xxxzzz_xxyyz[i] = 2.0 * tr_x_xxxz_xxyyz[i] * fe_0 + tr_x_xxxzz_xxyy[i] * fe_0 + tr_x_xxxzz_xxyyz[i] * pa_z[i];

        tr_x_xxxzzz_xxyzz[i] = 2.0 * tr_x_xxxz_xxyzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xxyz[i] * fe_0 + tr_x_xxxzz_xxyzz[i] * pa_z[i];

        tr_x_xxxzzz_xxzzz[i] = 2.0 * tr_x_xxxz_xxzzz[i] * fe_0 + 3.0 * tr_x_xxxzz_xxzz[i] * fe_0 + tr_x_xxxzz_xxzzz[i] * pa_z[i];

        tr_x_xxxzzz_xyyyy[i] = 2.0 * tr_x_xxxz_xyyyy[i] * fe_0 + tr_x_xxxzz_xyyyy[i] * pa_z[i];

        tr_x_xxxzzz_xyyyz[i] = 2.0 * tr_x_xxxz_xyyyz[i] * fe_0 + tr_x_xxxzz_xyyy[i] * fe_0 + tr_x_xxxzz_xyyyz[i] * pa_z[i];

        tr_x_xxxzzz_xyyzz[i] = 2.0 * tr_x_xxxz_xyyzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xyyz[i] * fe_0 + tr_x_xxxzz_xyyzz[i] * pa_z[i];

        tr_x_xxxzzz_xyzzz[i] = 2.0 * tr_x_xxxz_xyzzz[i] * fe_0 + 3.0 * tr_x_xxxzz_xyzz[i] * fe_0 + tr_x_xxxzz_xyzzz[i] * pa_z[i];

        tr_x_xxxzzz_xzzzz[i] = 2.0 * tr_x_xxxz_xzzzz[i] * fe_0 + 4.0 * tr_x_xxxzz_xzzz[i] * fe_0 + tr_x_xxxzz_xzzzz[i] * pa_z[i];

        tr_x_xxxzzz_yyyyy[i] = 2.0 * tr_x_xxxz_yyyyy[i] * fe_0 + tr_x_xxxzz_yyyyy[i] * pa_z[i];

        tr_x_xxxzzz_yyyyz[i] = 2.0 * tr_x_xzzz_yyyyz[i] * fe_0 + ts_xxzzz_yyyyz[i] * fe_0 + tr_x_xxzzz_yyyyz[i] * pa_x[i];

        tr_x_xxxzzz_yyyzz[i] = 2.0 * tr_x_xzzz_yyyzz[i] * fe_0 + ts_xxzzz_yyyzz[i] * fe_0 + tr_x_xxzzz_yyyzz[i] * pa_x[i];

        tr_x_xxxzzz_yyzzz[i] = 2.0 * tr_x_xzzz_yyzzz[i] * fe_0 + ts_xxzzz_yyzzz[i] * fe_0 + tr_x_xxzzz_yyzzz[i] * pa_x[i];

        tr_x_xxxzzz_yzzzz[i] = 2.0 * tr_x_xzzz_yzzzz[i] * fe_0 + ts_xxzzz_yzzzz[i] * fe_0 + tr_x_xxzzz_yzzzz[i] * pa_x[i];

        tr_x_xxxzzz_zzzzz[i] = 2.0 * tr_x_xzzz_zzzzz[i] * fe_0 + ts_xxzzz_zzzzz[i] * fe_0 + tr_x_xxzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : IH

    auto tr_x_xxyyyy_xxxxx = pbuffer.data(idx_dip_ih + 210);

    auto tr_x_xxyyyy_xxxxy = pbuffer.data(idx_dip_ih + 211);

    auto tr_x_xxyyyy_xxxxz = pbuffer.data(idx_dip_ih + 212);

    auto tr_x_xxyyyy_xxxyy = pbuffer.data(idx_dip_ih + 213);

    auto tr_x_xxyyyy_xxxyz = pbuffer.data(idx_dip_ih + 214);

    auto tr_x_xxyyyy_xxxzz = pbuffer.data(idx_dip_ih + 215);

    auto tr_x_xxyyyy_xxyyy = pbuffer.data(idx_dip_ih + 216);

    auto tr_x_xxyyyy_xxyyz = pbuffer.data(idx_dip_ih + 217);

    auto tr_x_xxyyyy_xxyzz = pbuffer.data(idx_dip_ih + 218);

    auto tr_x_xxyyyy_xxzzz = pbuffer.data(idx_dip_ih + 219);

    auto tr_x_xxyyyy_xyyyy = pbuffer.data(idx_dip_ih + 220);

    auto tr_x_xxyyyy_xyyyz = pbuffer.data(idx_dip_ih + 221);

    auto tr_x_xxyyyy_xyyzz = pbuffer.data(idx_dip_ih + 222);

    auto tr_x_xxyyyy_xyzzz = pbuffer.data(idx_dip_ih + 223);

    auto tr_x_xxyyyy_xzzzz = pbuffer.data(idx_dip_ih + 224);

    auto tr_x_xxyyyy_yyyyy = pbuffer.data(idx_dip_ih + 225);

    auto tr_x_xxyyyy_yyyyz = pbuffer.data(idx_dip_ih + 226);

    auto tr_x_xxyyyy_yyyzz = pbuffer.data(idx_dip_ih + 227);

    auto tr_x_xxyyyy_yyzzz = pbuffer.data(idx_dip_ih + 228);

    auto tr_x_xxyyyy_yzzzz = pbuffer.data(idx_dip_ih + 229);

    auto tr_x_xxyyyy_zzzzz = pbuffer.data(idx_dip_ih + 230);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxyy_xxxxx, tr_x_xxyy_xxxxy, tr_x_xxyy_xxxxz, tr_x_xxyy_xxxyy, tr_x_xxyy_xxxyz, tr_x_xxyy_xxxzz, tr_x_xxyy_xxyyy, tr_x_xxyy_xxyyz, tr_x_xxyy_xxyzz, tr_x_xxyy_xxzzz, tr_x_xxyy_xyyyy, tr_x_xxyy_xyyyz, tr_x_xxyy_xyyzz, tr_x_xxyy_xyzzz, tr_x_xxyy_xzzzz, tr_x_xxyy_zzzzz, tr_x_xxyyy_xxxx, tr_x_xxyyy_xxxxx, tr_x_xxyyy_xxxxy, tr_x_xxyyy_xxxxz, tr_x_xxyyy_xxxy, tr_x_xxyyy_xxxyy, tr_x_xxyyy_xxxyz, tr_x_xxyyy_xxxz, tr_x_xxyyy_xxxzz, tr_x_xxyyy_xxyy, tr_x_xxyyy_xxyyy, tr_x_xxyyy_xxyyz, tr_x_xxyyy_xxyz, tr_x_xxyyy_xxyzz, tr_x_xxyyy_xxzz, tr_x_xxyyy_xxzzz, tr_x_xxyyy_xyyy, tr_x_xxyyy_xyyyy, tr_x_xxyyy_xyyyz, tr_x_xxyyy_xyyz, tr_x_xxyyy_xyyzz, tr_x_xxyyy_xyzz, tr_x_xxyyy_xyzzz, tr_x_xxyyy_xzzz, tr_x_xxyyy_xzzzz, tr_x_xxyyy_zzzzz, tr_x_xxyyyy_xxxxx, tr_x_xxyyyy_xxxxy, tr_x_xxyyyy_xxxxz, tr_x_xxyyyy_xxxyy, tr_x_xxyyyy_xxxyz, tr_x_xxyyyy_xxxzz, tr_x_xxyyyy_xxyyy, tr_x_xxyyyy_xxyyz, tr_x_xxyyyy_xxyzz, tr_x_xxyyyy_xxzzz, tr_x_xxyyyy_xyyyy, tr_x_xxyyyy_xyyyz, tr_x_xxyyyy_xyyzz, tr_x_xxyyyy_xyzzz, tr_x_xxyyyy_xzzzz, tr_x_xxyyyy_yyyyy, tr_x_xxyyyy_yyyyz, tr_x_xxyyyy_yyyzz, tr_x_xxyyyy_yyzzz, tr_x_xxyyyy_yzzzz, tr_x_xxyyyy_zzzzz, tr_x_xyyyy_yyyyy, tr_x_xyyyy_yyyyz, tr_x_xyyyy_yyyzz, tr_x_xyyyy_yyzzz, tr_x_xyyyy_yzzzz, tr_x_yyyy_yyyyy, tr_x_yyyy_yyyyz, tr_x_yyyy_yyyzz, tr_x_yyyy_yyzzz, tr_x_yyyy_yzzzz, ts_xyyyy_yyyyy, ts_xyyyy_yyyyz, ts_xyyyy_yyyzz, ts_xyyyy_yyzzz, ts_xyyyy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyy_xxxxx[i] = 3.0 * tr_x_xxyy_xxxxx[i] * fe_0 + tr_x_xxyyy_xxxxx[i] * pa_y[i];

        tr_x_xxyyyy_xxxxy[i] = 3.0 * tr_x_xxyy_xxxxy[i] * fe_0 + tr_x_xxyyy_xxxx[i] * fe_0 + tr_x_xxyyy_xxxxy[i] * pa_y[i];

        tr_x_xxyyyy_xxxxz[i] = 3.0 * tr_x_xxyy_xxxxz[i] * fe_0 + tr_x_xxyyy_xxxxz[i] * pa_y[i];

        tr_x_xxyyyy_xxxyy[i] = 3.0 * tr_x_xxyy_xxxyy[i] * fe_0 + 2.0 * tr_x_xxyyy_xxxy[i] * fe_0 + tr_x_xxyyy_xxxyy[i] * pa_y[i];

        tr_x_xxyyyy_xxxyz[i] = 3.0 * tr_x_xxyy_xxxyz[i] * fe_0 + tr_x_xxyyy_xxxz[i] * fe_0 + tr_x_xxyyy_xxxyz[i] * pa_y[i];

        tr_x_xxyyyy_xxxzz[i] = 3.0 * tr_x_xxyy_xxxzz[i] * fe_0 + tr_x_xxyyy_xxxzz[i] * pa_y[i];

        tr_x_xxyyyy_xxyyy[i] = 3.0 * tr_x_xxyy_xxyyy[i] * fe_0 + 3.0 * tr_x_xxyyy_xxyy[i] * fe_0 + tr_x_xxyyy_xxyyy[i] * pa_y[i];

        tr_x_xxyyyy_xxyyz[i] = 3.0 * tr_x_xxyy_xxyyz[i] * fe_0 + 2.0 * tr_x_xxyyy_xxyz[i] * fe_0 + tr_x_xxyyy_xxyyz[i] * pa_y[i];

        tr_x_xxyyyy_xxyzz[i] = 3.0 * tr_x_xxyy_xxyzz[i] * fe_0 + tr_x_xxyyy_xxzz[i] * fe_0 + tr_x_xxyyy_xxyzz[i] * pa_y[i];

        tr_x_xxyyyy_xxzzz[i] = 3.0 * tr_x_xxyy_xxzzz[i] * fe_0 + tr_x_xxyyy_xxzzz[i] * pa_y[i];

        tr_x_xxyyyy_xyyyy[i] = 3.0 * tr_x_xxyy_xyyyy[i] * fe_0 + 4.0 * tr_x_xxyyy_xyyy[i] * fe_0 + tr_x_xxyyy_xyyyy[i] * pa_y[i];

        tr_x_xxyyyy_xyyyz[i] = 3.0 * tr_x_xxyy_xyyyz[i] * fe_0 + 3.0 * tr_x_xxyyy_xyyz[i] * fe_0 + tr_x_xxyyy_xyyyz[i] * pa_y[i];

        tr_x_xxyyyy_xyyzz[i] = 3.0 * tr_x_xxyy_xyyzz[i] * fe_0 + 2.0 * tr_x_xxyyy_xyzz[i] * fe_0 + tr_x_xxyyy_xyyzz[i] * pa_y[i];

        tr_x_xxyyyy_xyzzz[i] = 3.0 * tr_x_xxyy_xyzzz[i] * fe_0 + tr_x_xxyyy_xzzz[i] * fe_0 + tr_x_xxyyy_xyzzz[i] * pa_y[i];

        tr_x_xxyyyy_xzzzz[i] = 3.0 * tr_x_xxyy_xzzzz[i] * fe_0 + tr_x_xxyyy_xzzzz[i] * pa_y[i];

        tr_x_xxyyyy_yyyyy[i] = tr_x_yyyy_yyyyy[i] * fe_0 + ts_xyyyy_yyyyy[i] * fe_0 + tr_x_xyyyy_yyyyy[i] * pa_x[i];

        tr_x_xxyyyy_yyyyz[i] = tr_x_yyyy_yyyyz[i] * fe_0 + ts_xyyyy_yyyyz[i] * fe_0 + tr_x_xyyyy_yyyyz[i] * pa_x[i];

        tr_x_xxyyyy_yyyzz[i] = tr_x_yyyy_yyyzz[i] * fe_0 + ts_xyyyy_yyyzz[i] * fe_0 + tr_x_xyyyy_yyyzz[i] * pa_x[i];

        tr_x_xxyyyy_yyzzz[i] = tr_x_yyyy_yyzzz[i] * fe_0 + ts_xyyyy_yyzzz[i] * fe_0 + tr_x_xyyyy_yyzzz[i] * pa_x[i];

        tr_x_xxyyyy_yzzzz[i] = tr_x_yyyy_yzzzz[i] * fe_0 + ts_xyyyy_yzzzz[i] * fe_0 + tr_x_xyyyy_yzzzz[i] * pa_x[i];

        tr_x_xxyyyy_zzzzz[i] = 3.0 * tr_x_xxyy_zzzzz[i] * fe_0 + tr_x_xxyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 231-252 components of targeted buffer : IH

    auto tr_x_xxyyyz_xxxxx = pbuffer.data(idx_dip_ih + 231);

    auto tr_x_xxyyyz_xxxxy = pbuffer.data(idx_dip_ih + 232);

    auto tr_x_xxyyyz_xxxxz = pbuffer.data(idx_dip_ih + 233);

    auto tr_x_xxyyyz_xxxyy = pbuffer.data(idx_dip_ih + 234);

    auto tr_x_xxyyyz_xxxyz = pbuffer.data(idx_dip_ih + 235);

    auto tr_x_xxyyyz_xxxzz = pbuffer.data(idx_dip_ih + 236);

    auto tr_x_xxyyyz_xxyyy = pbuffer.data(idx_dip_ih + 237);

    auto tr_x_xxyyyz_xxyyz = pbuffer.data(idx_dip_ih + 238);

    auto tr_x_xxyyyz_xxyzz = pbuffer.data(idx_dip_ih + 239);

    auto tr_x_xxyyyz_xxzzz = pbuffer.data(idx_dip_ih + 240);

    auto tr_x_xxyyyz_xyyyy = pbuffer.data(idx_dip_ih + 241);

    auto tr_x_xxyyyz_xyyyz = pbuffer.data(idx_dip_ih + 242);

    auto tr_x_xxyyyz_xyyzz = pbuffer.data(idx_dip_ih + 243);

    auto tr_x_xxyyyz_xyzzz = pbuffer.data(idx_dip_ih + 244);

    auto tr_x_xxyyyz_xzzzz = pbuffer.data(idx_dip_ih + 245);

    auto tr_x_xxyyyz_yyyyy = pbuffer.data(idx_dip_ih + 246);

    auto tr_x_xxyyyz_yyyyz = pbuffer.data(idx_dip_ih + 247);

    auto tr_x_xxyyyz_yyyzz = pbuffer.data(idx_dip_ih + 248);

    auto tr_x_xxyyyz_yyzzz = pbuffer.data(idx_dip_ih + 249);

    auto tr_x_xxyyyz_yzzzz = pbuffer.data(idx_dip_ih + 250);

    auto tr_x_xxyyyz_zzzzz = pbuffer.data(idx_dip_ih + 251);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxyyy_xxxxx, tr_x_xxyyy_xxxxy, tr_x_xxyyy_xxxy, tr_x_xxyyy_xxxyy, tr_x_xxyyy_xxxyz, tr_x_xxyyy_xxyy, tr_x_xxyyy_xxyyy, tr_x_xxyyy_xxyyz, tr_x_xxyyy_xxyz, tr_x_xxyyy_xxyzz, tr_x_xxyyy_xyyy, tr_x_xxyyy_xyyyy, tr_x_xxyyy_xyyyz, tr_x_xxyyy_xyyz, tr_x_xxyyy_xyyzz, tr_x_xxyyy_xyzz, tr_x_xxyyy_xyzzz, tr_x_xxyyy_yyyy, tr_x_xxyyy_yyyyy, tr_x_xxyyy_yyyyz, tr_x_xxyyy_yyyz, tr_x_xxyyy_yyyzz, tr_x_xxyyy_yyzz, tr_x_xxyyy_yyzzz, tr_x_xxyyy_yzzz, tr_x_xxyyy_yzzzz, tr_x_xxyyyz_xxxxx, tr_x_xxyyyz_xxxxy, tr_x_xxyyyz_xxxxz, tr_x_xxyyyz_xxxyy, tr_x_xxyyyz_xxxyz, tr_x_xxyyyz_xxxzz, tr_x_xxyyyz_xxyyy, tr_x_xxyyyz_xxyyz, tr_x_xxyyyz_xxyzz, tr_x_xxyyyz_xxzzz, tr_x_xxyyyz_xyyyy, tr_x_xxyyyz_xyyyz, tr_x_xxyyyz_xyyzz, tr_x_xxyyyz_xyzzz, tr_x_xxyyyz_xzzzz, tr_x_xxyyyz_yyyyy, tr_x_xxyyyz_yyyyz, tr_x_xxyyyz_yyyzz, tr_x_xxyyyz_yyzzz, tr_x_xxyyyz_yzzzz, tr_x_xxyyyz_zzzzz, tr_x_xxyyz_xxxxz, tr_x_xxyyz_xxxzz, tr_x_xxyyz_xxzzz, tr_x_xxyyz_xzzzz, tr_x_xxyyz_zzzzz, tr_x_xxyz_xxxxz, tr_x_xxyz_xxxzz, tr_x_xxyz_xxzzz, tr_x_xxyz_xzzzz, tr_x_xxyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyz_xxxxx[i] = tr_x_xxyyy_xxxxx[i] * pa_z[i];

        tr_x_xxyyyz_xxxxy[i] = tr_x_xxyyy_xxxxy[i] * pa_z[i];

        tr_x_xxyyyz_xxxxz[i] = 2.0 * tr_x_xxyz_xxxxz[i] * fe_0 + tr_x_xxyyz_xxxxz[i] * pa_y[i];

        tr_x_xxyyyz_xxxyy[i] = tr_x_xxyyy_xxxyy[i] * pa_z[i];

        tr_x_xxyyyz_xxxyz[i] = tr_x_xxyyy_xxxy[i] * fe_0 + tr_x_xxyyy_xxxyz[i] * pa_z[i];

        tr_x_xxyyyz_xxxzz[i] = 2.0 * tr_x_xxyz_xxxzz[i] * fe_0 + tr_x_xxyyz_xxxzz[i] * pa_y[i];

        tr_x_xxyyyz_xxyyy[i] = tr_x_xxyyy_xxyyy[i] * pa_z[i];

        tr_x_xxyyyz_xxyyz[i] = tr_x_xxyyy_xxyy[i] * fe_0 + tr_x_xxyyy_xxyyz[i] * pa_z[i];

        tr_x_xxyyyz_xxyzz[i] = 2.0 * tr_x_xxyyy_xxyz[i] * fe_0 + tr_x_xxyyy_xxyzz[i] * pa_z[i];

        tr_x_xxyyyz_xxzzz[i] = 2.0 * tr_x_xxyz_xxzzz[i] * fe_0 + tr_x_xxyyz_xxzzz[i] * pa_y[i];

        tr_x_xxyyyz_xyyyy[i] = tr_x_xxyyy_xyyyy[i] * pa_z[i];

        tr_x_xxyyyz_xyyyz[i] = tr_x_xxyyy_xyyy[i] * fe_0 + tr_x_xxyyy_xyyyz[i] * pa_z[i];

        tr_x_xxyyyz_xyyzz[i] = 2.0 * tr_x_xxyyy_xyyz[i] * fe_0 + tr_x_xxyyy_xyyzz[i] * pa_z[i];

        tr_x_xxyyyz_xyzzz[i] = 3.0 * tr_x_xxyyy_xyzz[i] * fe_0 + tr_x_xxyyy_xyzzz[i] * pa_z[i];

        tr_x_xxyyyz_xzzzz[i] = 2.0 * tr_x_xxyz_xzzzz[i] * fe_0 + tr_x_xxyyz_xzzzz[i] * pa_y[i];

        tr_x_xxyyyz_yyyyy[i] = tr_x_xxyyy_yyyyy[i] * pa_z[i];

        tr_x_xxyyyz_yyyyz[i] = tr_x_xxyyy_yyyy[i] * fe_0 + tr_x_xxyyy_yyyyz[i] * pa_z[i];

        tr_x_xxyyyz_yyyzz[i] = 2.0 * tr_x_xxyyy_yyyz[i] * fe_0 + tr_x_xxyyy_yyyzz[i] * pa_z[i];

        tr_x_xxyyyz_yyzzz[i] = 3.0 * tr_x_xxyyy_yyzz[i] * fe_0 + tr_x_xxyyy_yyzzz[i] * pa_z[i];

        tr_x_xxyyyz_yzzzz[i] = 4.0 * tr_x_xxyyy_yzzz[i] * fe_0 + tr_x_xxyyy_yzzzz[i] * pa_z[i];

        tr_x_xxyyyz_zzzzz[i] = 2.0 * tr_x_xxyz_zzzzz[i] * fe_0 + tr_x_xxyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 252-273 components of targeted buffer : IH

    auto tr_x_xxyyzz_xxxxx = pbuffer.data(idx_dip_ih + 252);

    auto tr_x_xxyyzz_xxxxy = pbuffer.data(idx_dip_ih + 253);

    auto tr_x_xxyyzz_xxxxz = pbuffer.data(idx_dip_ih + 254);

    auto tr_x_xxyyzz_xxxyy = pbuffer.data(idx_dip_ih + 255);

    auto tr_x_xxyyzz_xxxyz = pbuffer.data(idx_dip_ih + 256);

    auto tr_x_xxyyzz_xxxzz = pbuffer.data(idx_dip_ih + 257);

    auto tr_x_xxyyzz_xxyyy = pbuffer.data(idx_dip_ih + 258);

    auto tr_x_xxyyzz_xxyyz = pbuffer.data(idx_dip_ih + 259);

    auto tr_x_xxyyzz_xxyzz = pbuffer.data(idx_dip_ih + 260);

    auto tr_x_xxyyzz_xxzzz = pbuffer.data(idx_dip_ih + 261);

    auto tr_x_xxyyzz_xyyyy = pbuffer.data(idx_dip_ih + 262);

    auto tr_x_xxyyzz_xyyyz = pbuffer.data(idx_dip_ih + 263);

    auto tr_x_xxyyzz_xyyzz = pbuffer.data(idx_dip_ih + 264);

    auto tr_x_xxyyzz_xyzzz = pbuffer.data(idx_dip_ih + 265);

    auto tr_x_xxyyzz_xzzzz = pbuffer.data(idx_dip_ih + 266);

    auto tr_x_xxyyzz_yyyyy = pbuffer.data(idx_dip_ih + 267);

    auto tr_x_xxyyzz_yyyyz = pbuffer.data(idx_dip_ih + 268);

    auto tr_x_xxyyzz_yyyzz = pbuffer.data(idx_dip_ih + 269);

    auto tr_x_xxyyzz_yyzzz = pbuffer.data(idx_dip_ih + 270);

    auto tr_x_xxyyzz_yzzzz = pbuffer.data(idx_dip_ih + 271);

    auto tr_x_xxyyzz_zzzzz = pbuffer.data(idx_dip_ih + 272);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xxyy_xxxxy, tr_x_xxyy_xxxyy, tr_x_xxyy_xxyyy, tr_x_xxyy_xyyyy, tr_x_xxyy_yyyyy, tr_x_xxyyz_xxxxy, tr_x_xxyyz_xxxyy, tr_x_xxyyz_xxyyy, tr_x_xxyyz_xyyyy, tr_x_xxyyz_yyyyy, tr_x_xxyyzz_xxxxx, tr_x_xxyyzz_xxxxy, tr_x_xxyyzz_xxxxz, tr_x_xxyyzz_xxxyy, tr_x_xxyyzz_xxxyz, tr_x_xxyyzz_xxxzz, tr_x_xxyyzz_xxyyy, tr_x_xxyyzz_xxyyz, tr_x_xxyyzz_xxyzz, tr_x_xxyyzz_xxzzz, tr_x_xxyyzz_xyyyy, tr_x_xxyyzz_xyyyz, tr_x_xxyyzz_xyyzz, tr_x_xxyyzz_xyzzz, tr_x_xxyyzz_xzzzz, tr_x_xxyyzz_yyyyy, tr_x_xxyyzz_yyyyz, tr_x_xxyyzz_yyyzz, tr_x_xxyyzz_yyzzz, tr_x_xxyyzz_yzzzz, tr_x_xxyyzz_zzzzz, tr_x_xxyzz_xxxxx, tr_x_xxyzz_xxxxz, tr_x_xxyzz_xxxyz, tr_x_xxyzz_xxxz, tr_x_xxyzz_xxxzz, tr_x_xxyzz_xxyyz, tr_x_xxyzz_xxyz, tr_x_xxyzz_xxyzz, tr_x_xxyzz_xxzz, tr_x_xxyzz_xxzzz, tr_x_xxyzz_xyyyz, tr_x_xxyzz_xyyz, tr_x_xxyzz_xyyzz, tr_x_xxyzz_xyzz, tr_x_xxyzz_xyzzz, tr_x_xxyzz_xzzz, tr_x_xxyzz_xzzzz, tr_x_xxyzz_zzzzz, tr_x_xxzz_xxxxx, tr_x_xxzz_xxxxz, tr_x_xxzz_xxxyz, tr_x_xxzz_xxxzz, tr_x_xxzz_xxyyz, tr_x_xxzz_xxyzz, tr_x_xxzz_xxzzz, tr_x_xxzz_xyyyz, tr_x_xxzz_xyyzz, tr_x_xxzz_xyzzz, tr_x_xxzz_xzzzz, tr_x_xxzz_zzzzz, tr_x_xyyzz_yyyyz, tr_x_xyyzz_yyyzz, tr_x_xyyzz_yyzzz, tr_x_xyyzz_yzzzz, tr_x_yyzz_yyyyz, tr_x_yyzz_yyyzz, tr_x_yyzz_yyzzz, tr_x_yyzz_yzzzz, ts_xyyzz_yyyyz, ts_xyyzz_yyyzz, ts_xyyzz_yyzzz, ts_xyyzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyzz_xxxxx[i] = tr_x_xxzz_xxxxx[i] * fe_0 + tr_x_xxyzz_xxxxx[i] * pa_y[i];

        tr_x_xxyyzz_xxxxy[i] = tr_x_xxyy_xxxxy[i] * fe_0 + tr_x_xxyyz_xxxxy[i] * pa_z[i];

        tr_x_xxyyzz_xxxxz[i] = tr_x_xxzz_xxxxz[i] * fe_0 + tr_x_xxyzz_xxxxz[i] * pa_y[i];

        tr_x_xxyyzz_xxxyy[i] = tr_x_xxyy_xxxyy[i] * fe_0 + tr_x_xxyyz_xxxyy[i] * pa_z[i];

        tr_x_xxyyzz_xxxyz[i] = tr_x_xxzz_xxxyz[i] * fe_0 + tr_x_xxyzz_xxxz[i] * fe_0 + tr_x_xxyzz_xxxyz[i] * pa_y[i];

        tr_x_xxyyzz_xxxzz[i] = tr_x_xxzz_xxxzz[i] * fe_0 + tr_x_xxyzz_xxxzz[i] * pa_y[i];

        tr_x_xxyyzz_xxyyy[i] = tr_x_xxyy_xxyyy[i] * fe_0 + tr_x_xxyyz_xxyyy[i] * pa_z[i];

        tr_x_xxyyzz_xxyyz[i] = tr_x_xxzz_xxyyz[i] * fe_0 + 2.0 * tr_x_xxyzz_xxyz[i] * fe_0 + tr_x_xxyzz_xxyyz[i] * pa_y[i];

        tr_x_xxyyzz_xxyzz[i] = tr_x_xxzz_xxyzz[i] * fe_0 + tr_x_xxyzz_xxzz[i] * fe_0 + tr_x_xxyzz_xxyzz[i] * pa_y[i];

        tr_x_xxyyzz_xxzzz[i] = tr_x_xxzz_xxzzz[i] * fe_0 + tr_x_xxyzz_xxzzz[i] * pa_y[i];

        tr_x_xxyyzz_xyyyy[i] = tr_x_xxyy_xyyyy[i] * fe_0 + tr_x_xxyyz_xyyyy[i] * pa_z[i];

        tr_x_xxyyzz_xyyyz[i] = tr_x_xxzz_xyyyz[i] * fe_0 + 3.0 * tr_x_xxyzz_xyyz[i] * fe_0 + tr_x_xxyzz_xyyyz[i] * pa_y[i];

        tr_x_xxyyzz_xyyzz[i] = tr_x_xxzz_xyyzz[i] * fe_0 + 2.0 * tr_x_xxyzz_xyzz[i] * fe_0 + tr_x_xxyzz_xyyzz[i] * pa_y[i];

        tr_x_xxyyzz_xyzzz[i] = tr_x_xxzz_xyzzz[i] * fe_0 + tr_x_xxyzz_xzzz[i] * fe_0 + tr_x_xxyzz_xyzzz[i] * pa_y[i];

        tr_x_xxyyzz_xzzzz[i] = tr_x_xxzz_xzzzz[i] * fe_0 + tr_x_xxyzz_xzzzz[i] * pa_y[i];

        tr_x_xxyyzz_yyyyy[i] = tr_x_xxyy_yyyyy[i] * fe_0 + tr_x_xxyyz_yyyyy[i] * pa_z[i];

        tr_x_xxyyzz_yyyyz[i] = tr_x_yyzz_yyyyz[i] * fe_0 + ts_xyyzz_yyyyz[i] * fe_0 + tr_x_xyyzz_yyyyz[i] * pa_x[i];

        tr_x_xxyyzz_yyyzz[i] = tr_x_yyzz_yyyzz[i] * fe_0 + ts_xyyzz_yyyzz[i] * fe_0 + tr_x_xyyzz_yyyzz[i] * pa_x[i];

        tr_x_xxyyzz_yyzzz[i] = tr_x_yyzz_yyzzz[i] * fe_0 + ts_xyyzz_yyzzz[i] * fe_0 + tr_x_xyyzz_yyzzz[i] * pa_x[i];

        tr_x_xxyyzz_yzzzz[i] = tr_x_yyzz_yzzzz[i] * fe_0 + ts_xyyzz_yzzzz[i] * fe_0 + tr_x_xyyzz_yzzzz[i] * pa_x[i];

        tr_x_xxyyzz_zzzzz[i] = tr_x_xxzz_zzzzz[i] * fe_0 + tr_x_xxyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 273-294 components of targeted buffer : IH

    auto tr_x_xxyzzz_xxxxx = pbuffer.data(idx_dip_ih + 273);

    auto tr_x_xxyzzz_xxxxy = pbuffer.data(idx_dip_ih + 274);

    auto tr_x_xxyzzz_xxxxz = pbuffer.data(idx_dip_ih + 275);

    auto tr_x_xxyzzz_xxxyy = pbuffer.data(idx_dip_ih + 276);

    auto tr_x_xxyzzz_xxxyz = pbuffer.data(idx_dip_ih + 277);

    auto tr_x_xxyzzz_xxxzz = pbuffer.data(idx_dip_ih + 278);

    auto tr_x_xxyzzz_xxyyy = pbuffer.data(idx_dip_ih + 279);

    auto tr_x_xxyzzz_xxyyz = pbuffer.data(idx_dip_ih + 280);

    auto tr_x_xxyzzz_xxyzz = pbuffer.data(idx_dip_ih + 281);

    auto tr_x_xxyzzz_xxzzz = pbuffer.data(idx_dip_ih + 282);

    auto tr_x_xxyzzz_xyyyy = pbuffer.data(idx_dip_ih + 283);

    auto tr_x_xxyzzz_xyyyz = pbuffer.data(idx_dip_ih + 284);

    auto tr_x_xxyzzz_xyyzz = pbuffer.data(idx_dip_ih + 285);

    auto tr_x_xxyzzz_xyzzz = pbuffer.data(idx_dip_ih + 286);

    auto tr_x_xxyzzz_xzzzz = pbuffer.data(idx_dip_ih + 287);

    auto tr_x_xxyzzz_yyyyy = pbuffer.data(idx_dip_ih + 288);

    auto tr_x_xxyzzz_yyyyz = pbuffer.data(idx_dip_ih + 289);

    auto tr_x_xxyzzz_yyyzz = pbuffer.data(idx_dip_ih + 290);

    auto tr_x_xxyzzz_yyzzz = pbuffer.data(idx_dip_ih + 291);

    auto tr_x_xxyzzz_yzzzz = pbuffer.data(idx_dip_ih + 292);

    auto tr_x_xxyzzz_zzzzz = pbuffer.data(idx_dip_ih + 293);

    #pragma omp simd aligned(pa_y, tr_x_xxyzzz_xxxxx, tr_x_xxyzzz_xxxxy, tr_x_xxyzzz_xxxxz, tr_x_xxyzzz_xxxyy, tr_x_xxyzzz_xxxyz, tr_x_xxyzzz_xxxzz, tr_x_xxyzzz_xxyyy, tr_x_xxyzzz_xxyyz, tr_x_xxyzzz_xxyzz, tr_x_xxyzzz_xxzzz, tr_x_xxyzzz_xyyyy, tr_x_xxyzzz_xyyyz, tr_x_xxyzzz_xyyzz, tr_x_xxyzzz_xyzzz, tr_x_xxyzzz_xzzzz, tr_x_xxyzzz_yyyyy, tr_x_xxyzzz_yyyyz, tr_x_xxyzzz_yyyzz, tr_x_xxyzzz_yyzzz, tr_x_xxyzzz_yzzzz, tr_x_xxyzzz_zzzzz, tr_x_xxzzz_xxxx, tr_x_xxzzz_xxxxx, tr_x_xxzzz_xxxxy, tr_x_xxzzz_xxxxz, tr_x_xxzzz_xxxy, tr_x_xxzzz_xxxyy, tr_x_xxzzz_xxxyz, tr_x_xxzzz_xxxz, tr_x_xxzzz_xxxzz, tr_x_xxzzz_xxyy, tr_x_xxzzz_xxyyy, tr_x_xxzzz_xxyyz, tr_x_xxzzz_xxyz, tr_x_xxzzz_xxyzz, tr_x_xxzzz_xxzz, tr_x_xxzzz_xxzzz, tr_x_xxzzz_xyyy, tr_x_xxzzz_xyyyy, tr_x_xxzzz_xyyyz, tr_x_xxzzz_xyyz, tr_x_xxzzz_xyyzz, tr_x_xxzzz_xyzz, tr_x_xxzzz_xyzzz, tr_x_xxzzz_xzzz, tr_x_xxzzz_xzzzz, tr_x_xxzzz_yyyy, tr_x_xxzzz_yyyyy, tr_x_xxzzz_yyyyz, tr_x_xxzzz_yyyz, tr_x_xxzzz_yyyzz, tr_x_xxzzz_yyzz, tr_x_xxzzz_yyzzz, tr_x_xxzzz_yzzz, tr_x_xxzzz_yzzzz, tr_x_xxzzz_zzzz, tr_x_xxzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzzz_xxxxx[i] = tr_x_xxzzz_xxxxx[i] * pa_y[i];

        tr_x_xxyzzz_xxxxy[i] = tr_x_xxzzz_xxxx[i] * fe_0 + tr_x_xxzzz_xxxxy[i] * pa_y[i];

        tr_x_xxyzzz_xxxxz[i] = tr_x_xxzzz_xxxxz[i] * pa_y[i];

        tr_x_xxyzzz_xxxyy[i] = 2.0 * tr_x_xxzzz_xxxy[i] * fe_0 + tr_x_xxzzz_xxxyy[i] * pa_y[i];

        tr_x_xxyzzz_xxxyz[i] = tr_x_xxzzz_xxxz[i] * fe_0 + tr_x_xxzzz_xxxyz[i] * pa_y[i];

        tr_x_xxyzzz_xxxzz[i] = tr_x_xxzzz_xxxzz[i] * pa_y[i];

        tr_x_xxyzzz_xxyyy[i] = 3.0 * tr_x_xxzzz_xxyy[i] * fe_0 + tr_x_xxzzz_xxyyy[i] * pa_y[i];

        tr_x_xxyzzz_xxyyz[i] = 2.0 * tr_x_xxzzz_xxyz[i] * fe_0 + tr_x_xxzzz_xxyyz[i] * pa_y[i];

        tr_x_xxyzzz_xxyzz[i] = tr_x_xxzzz_xxzz[i] * fe_0 + tr_x_xxzzz_xxyzz[i] * pa_y[i];

        tr_x_xxyzzz_xxzzz[i] = tr_x_xxzzz_xxzzz[i] * pa_y[i];

        tr_x_xxyzzz_xyyyy[i] = 4.0 * tr_x_xxzzz_xyyy[i] * fe_0 + tr_x_xxzzz_xyyyy[i] * pa_y[i];

        tr_x_xxyzzz_xyyyz[i] = 3.0 * tr_x_xxzzz_xyyz[i] * fe_0 + tr_x_xxzzz_xyyyz[i] * pa_y[i];

        tr_x_xxyzzz_xyyzz[i] = 2.0 * tr_x_xxzzz_xyzz[i] * fe_0 + tr_x_xxzzz_xyyzz[i] * pa_y[i];

        tr_x_xxyzzz_xyzzz[i] = tr_x_xxzzz_xzzz[i] * fe_0 + tr_x_xxzzz_xyzzz[i] * pa_y[i];

        tr_x_xxyzzz_xzzzz[i] = tr_x_xxzzz_xzzzz[i] * pa_y[i];

        tr_x_xxyzzz_yyyyy[i] = 5.0 * tr_x_xxzzz_yyyy[i] * fe_0 + tr_x_xxzzz_yyyyy[i] * pa_y[i];

        tr_x_xxyzzz_yyyyz[i] = 4.0 * tr_x_xxzzz_yyyz[i] * fe_0 + tr_x_xxzzz_yyyyz[i] * pa_y[i];

        tr_x_xxyzzz_yyyzz[i] = 3.0 * tr_x_xxzzz_yyzz[i] * fe_0 + tr_x_xxzzz_yyyzz[i] * pa_y[i];

        tr_x_xxyzzz_yyzzz[i] = 2.0 * tr_x_xxzzz_yzzz[i] * fe_0 + tr_x_xxzzz_yyzzz[i] * pa_y[i];

        tr_x_xxyzzz_yzzzz[i] = tr_x_xxzzz_zzzz[i] * fe_0 + tr_x_xxzzz_yzzzz[i] * pa_y[i];

        tr_x_xxyzzz_zzzzz[i] = tr_x_xxzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : IH

    auto tr_x_xxzzzz_xxxxx = pbuffer.data(idx_dip_ih + 294);

    auto tr_x_xxzzzz_xxxxy = pbuffer.data(idx_dip_ih + 295);

    auto tr_x_xxzzzz_xxxxz = pbuffer.data(idx_dip_ih + 296);

    auto tr_x_xxzzzz_xxxyy = pbuffer.data(idx_dip_ih + 297);

    auto tr_x_xxzzzz_xxxyz = pbuffer.data(idx_dip_ih + 298);

    auto tr_x_xxzzzz_xxxzz = pbuffer.data(idx_dip_ih + 299);

    auto tr_x_xxzzzz_xxyyy = pbuffer.data(idx_dip_ih + 300);

    auto tr_x_xxzzzz_xxyyz = pbuffer.data(idx_dip_ih + 301);

    auto tr_x_xxzzzz_xxyzz = pbuffer.data(idx_dip_ih + 302);

    auto tr_x_xxzzzz_xxzzz = pbuffer.data(idx_dip_ih + 303);

    auto tr_x_xxzzzz_xyyyy = pbuffer.data(idx_dip_ih + 304);

    auto tr_x_xxzzzz_xyyyz = pbuffer.data(idx_dip_ih + 305);

    auto tr_x_xxzzzz_xyyzz = pbuffer.data(idx_dip_ih + 306);

    auto tr_x_xxzzzz_xyzzz = pbuffer.data(idx_dip_ih + 307);

    auto tr_x_xxzzzz_xzzzz = pbuffer.data(idx_dip_ih + 308);

    auto tr_x_xxzzzz_yyyyy = pbuffer.data(idx_dip_ih + 309);

    auto tr_x_xxzzzz_yyyyz = pbuffer.data(idx_dip_ih + 310);

    auto tr_x_xxzzzz_yyyzz = pbuffer.data(idx_dip_ih + 311);

    auto tr_x_xxzzzz_yyzzz = pbuffer.data(idx_dip_ih + 312);

    auto tr_x_xxzzzz_yzzzz = pbuffer.data(idx_dip_ih + 313);

    auto tr_x_xxzzzz_zzzzz = pbuffer.data(idx_dip_ih + 314);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxzz_xxxxx, tr_x_xxzz_xxxxy, tr_x_xxzz_xxxxz, tr_x_xxzz_xxxyy, tr_x_xxzz_xxxyz, tr_x_xxzz_xxxzz, tr_x_xxzz_xxyyy, tr_x_xxzz_xxyyz, tr_x_xxzz_xxyzz, tr_x_xxzz_xxzzz, tr_x_xxzz_xyyyy, tr_x_xxzz_xyyyz, tr_x_xxzz_xyyzz, tr_x_xxzz_xyzzz, tr_x_xxzz_xzzzz, tr_x_xxzz_yyyyy, tr_x_xxzzz_xxxx, tr_x_xxzzz_xxxxx, tr_x_xxzzz_xxxxy, tr_x_xxzzz_xxxxz, tr_x_xxzzz_xxxy, tr_x_xxzzz_xxxyy, tr_x_xxzzz_xxxyz, tr_x_xxzzz_xxxz, tr_x_xxzzz_xxxzz, tr_x_xxzzz_xxyy, tr_x_xxzzz_xxyyy, tr_x_xxzzz_xxyyz, tr_x_xxzzz_xxyz, tr_x_xxzzz_xxyzz, tr_x_xxzzz_xxzz, tr_x_xxzzz_xxzzz, tr_x_xxzzz_xyyy, tr_x_xxzzz_xyyyy, tr_x_xxzzz_xyyyz, tr_x_xxzzz_xyyz, tr_x_xxzzz_xyyzz, tr_x_xxzzz_xyzz, tr_x_xxzzz_xyzzz, tr_x_xxzzz_xzzz, tr_x_xxzzz_xzzzz, tr_x_xxzzz_yyyyy, tr_x_xxzzzz_xxxxx, tr_x_xxzzzz_xxxxy, tr_x_xxzzzz_xxxxz, tr_x_xxzzzz_xxxyy, tr_x_xxzzzz_xxxyz, tr_x_xxzzzz_xxxzz, tr_x_xxzzzz_xxyyy, tr_x_xxzzzz_xxyyz, tr_x_xxzzzz_xxyzz, tr_x_xxzzzz_xxzzz, tr_x_xxzzzz_xyyyy, tr_x_xxzzzz_xyyyz, tr_x_xxzzzz_xyyzz, tr_x_xxzzzz_xyzzz, tr_x_xxzzzz_xzzzz, tr_x_xxzzzz_yyyyy, tr_x_xxzzzz_yyyyz, tr_x_xxzzzz_yyyzz, tr_x_xxzzzz_yyzzz, tr_x_xxzzzz_yzzzz, tr_x_xxzzzz_zzzzz, tr_x_xzzzz_yyyyz, tr_x_xzzzz_yyyzz, tr_x_xzzzz_yyzzz, tr_x_xzzzz_yzzzz, tr_x_xzzzz_zzzzz, tr_x_zzzz_yyyyz, tr_x_zzzz_yyyzz, tr_x_zzzz_yyzzz, tr_x_zzzz_yzzzz, tr_x_zzzz_zzzzz, ts_xzzzz_yyyyz, ts_xzzzz_yyyzz, ts_xzzzz_yyzzz, ts_xzzzz_yzzzz, ts_xzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzzz_xxxxx[i] = 3.0 * tr_x_xxzz_xxxxx[i] * fe_0 + tr_x_xxzzz_xxxxx[i] * pa_z[i];

        tr_x_xxzzzz_xxxxy[i] = 3.0 * tr_x_xxzz_xxxxy[i] * fe_0 + tr_x_xxzzz_xxxxy[i] * pa_z[i];

        tr_x_xxzzzz_xxxxz[i] = 3.0 * tr_x_xxzz_xxxxz[i] * fe_0 + tr_x_xxzzz_xxxx[i] * fe_0 + tr_x_xxzzz_xxxxz[i] * pa_z[i];

        tr_x_xxzzzz_xxxyy[i] = 3.0 * tr_x_xxzz_xxxyy[i] * fe_0 + tr_x_xxzzz_xxxyy[i] * pa_z[i];

        tr_x_xxzzzz_xxxyz[i] = 3.0 * tr_x_xxzz_xxxyz[i] * fe_0 + tr_x_xxzzz_xxxy[i] * fe_0 + tr_x_xxzzz_xxxyz[i] * pa_z[i];

        tr_x_xxzzzz_xxxzz[i] = 3.0 * tr_x_xxzz_xxxzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xxxz[i] * fe_0 + tr_x_xxzzz_xxxzz[i] * pa_z[i];

        tr_x_xxzzzz_xxyyy[i] = 3.0 * tr_x_xxzz_xxyyy[i] * fe_0 + tr_x_xxzzz_xxyyy[i] * pa_z[i];

        tr_x_xxzzzz_xxyyz[i] = 3.0 * tr_x_xxzz_xxyyz[i] * fe_0 + tr_x_xxzzz_xxyy[i] * fe_0 + tr_x_xxzzz_xxyyz[i] * pa_z[i];

        tr_x_xxzzzz_xxyzz[i] = 3.0 * tr_x_xxzz_xxyzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xxyz[i] * fe_0 + tr_x_xxzzz_xxyzz[i] * pa_z[i];

        tr_x_xxzzzz_xxzzz[i] = 3.0 * tr_x_xxzz_xxzzz[i] * fe_0 + 3.0 * tr_x_xxzzz_xxzz[i] * fe_0 + tr_x_xxzzz_xxzzz[i] * pa_z[i];

        tr_x_xxzzzz_xyyyy[i] = 3.0 * tr_x_xxzz_xyyyy[i] * fe_0 + tr_x_xxzzz_xyyyy[i] * pa_z[i];

        tr_x_xxzzzz_xyyyz[i] = 3.0 * tr_x_xxzz_xyyyz[i] * fe_0 + tr_x_xxzzz_xyyy[i] * fe_0 + tr_x_xxzzz_xyyyz[i] * pa_z[i];

        tr_x_xxzzzz_xyyzz[i] = 3.0 * tr_x_xxzz_xyyzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xyyz[i] * fe_0 + tr_x_xxzzz_xyyzz[i] * pa_z[i];

        tr_x_xxzzzz_xyzzz[i] = 3.0 * tr_x_xxzz_xyzzz[i] * fe_0 + 3.0 * tr_x_xxzzz_xyzz[i] * fe_0 + tr_x_xxzzz_xyzzz[i] * pa_z[i];

        tr_x_xxzzzz_xzzzz[i] = 3.0 * tr_x_xxzz_xzzzz[i] * fe_0 + 4.0 * tr_x_xxzzz_xzzz[i] * fe_0 + tr_x_xxzzz_xzzzz[i] * pa_z[i];

        tr_x_xxzzzz_yyyyy[i] = 3.0 * tr_x_xxzz_yyyyy[i] * fe_0 + tr_x_xxzzz_yyyyy[i] * pa_z[i];

        tr_x_xxzzzz_yyyyz[i] = tr_x_zzzz_yyyyz[i] * fe_0 + ts_xzzzz_yyyyz[i] * fe_0 + tr_x_xzzzz_yyyyz[i] * pa_x[i];

        tr_x_xxzzzz_yyyzz[i] = tr_x_zzzz_yyyzz[i] * fe_0 + ts_xzzzz_yyyzz[i] * fe_0 + tr_x_xzzzz_yyyzz[i] * pa_x[i];

        tr_x_xxzzzz_yyzzz[i] = tr_x_zzzz_yyzzz[i] * fe_0 + ts_xzzzz_yyzzz[i] * fe_0 + tr_x_xzzzz_yyzzz[i] * pa_x[i];

        tr_x_xxzzzz_yzzzz[i] = tr_x_zzzz_yzzzz[i] * fe_0 + ts_xzzzz_yzzzz[i] * fe_0 + tr_x_xzzzz_yzzzz[i] * pa_x[i];

        tr_x_xxzzzz_zzzzz[i] = tr_x_zzzz_zzzzz[i] * fe_0 + ts_xzzzz_zzzzz[i] * fe_0 + tr_x_xzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : IH

    auto tr_x_xyyyyy_xxxxx = pbuffer.data(idx_dip_ih + 315);

    auto tr_x_xyyyyy_xxxxy = pbuffer.data(idx_dip_ih + 316);

    auto tr_x_xyyyyy_xxxxz = pbuffer.data(idx_dip_ih + 317);

    auto tr_x_xyyyyy_xxxyy = pbuffer.data(idx_dip_ih + 318);

    auto tr_x_xyyyyy_xxxyz = pbuffer.data(idx_dip_ih + 319);

    auto tr_x_xyyyyy_xxxzz = pbuffer.data(idx_dip_ih + 320);

    auto tr_x_xyyyyy_xxyyy = pbuffer.data(idx_dip_ih + 321);

    auto tr_x_xyyyyy_xxyyz = pbuffer.data(idx_dip_ih + 322);

    auto tr_x_xyyyyy_xxyzz = pbuffer.data(idx_dip_ih + 323);

    auto tr_x_xyyyyy_xxzzz = pbuffer.data(idx_dip_ih + 324);

    auto tr_x_xyyyyy_xyyyy = pbuffer.data(idx_dip_ih + 325);

    auto tr_x_xyyyyy_xyyyz = pbuffer.data(idx_dip_ih + 326);

    auto tr_x_xyyyyy_xyyzz = pbuffer.data(idx_dip_ih + 327);

    auto tr_x_xyyyyy_xyzzz = pbuffer.data(idx_dip_ih + 328);

    auto tr_x_xyyyyy_xzzzz = pbuffer.data(idx_dip_ih + 329);

    auto tr_x_xyyyyy_yyyyy = pbuffer.data(idx_dip_ih + 330);

    auto tr_x_xyyyyy_yyyyz = pbuffer.data(idx_dip_ih + 331);

    auto tr_x_xyyyyy_yyyzz = pbuffer.data(idx_dip_ih + 332);

    auto tr_x_xyyyyy_yyzzz = pbuffer.data(idx_dip_ih + 333);

    auto tr_x_xyyyyy_yzzzz = pbuffer.data(idx_dip_ih + 334);

    auto tr_x_xyyyyy_zzzzz = pbuffer.data(idx_dip_ih + 335);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyyy_xxxxx, tr_x_xyyy_xxxxz, tr_x_xyyy_xxxzz, tr_x_xyyy_xxzzz, tr_x_xyyy_xzzzz, tr_x_xyyyy_xxxxx, tr_x_xyyyy_xxxxz, tr_x_xyyyy_xxxzz, tr_x_xyyyy_xxzzz, tr_x_xyyyy_xzzzz, tr_x_xyyyyy_xxxxx, tr_x_xyyyyy_xxxxy, tr_x_xyyyyy_xxxxz, tr_x_xyyyyy_xxxyy, tr_x_xyyyyy_xxxyz, tr_x_xyyyyy_xxxzz, tr_x_xyyyyy_xxyyy, tr_x_xyyyyy_xxyyz, tr_x_xyyyyy_xxyzz, tr_x_xyyyyy_xxzzz, tr_x_xyyyyy_xyyyy, tr_x_xyyyyy_xyyyz, tr_x_xyyyyy_xyyzz, tr_x_xyyyyy_xyzzz, tr_x_xyyyyy_xzzzz, tr_x_xyyyyy_yyyyy, tr_x_xyyyyy_yyyyz, tr_x_xyyyyy_yyyzz, tr_x_xyyyyy_yyzzz, tr_x_xyyyyy_yzzzz, tr_x_xyyyyy_zzzzz, tr_x_yyyyy_xxxxy, tr_x_yyyyy_xxxy, tr_x_yyyyy_xxxyy, tr_x_yyyyy_xxxyz, tr_x_yyyyy_xxyy, tr_x_yyyyy_xxyyy, tr_x_yyyyy_xxyyz, tr_x_yyyyy_xxyz, tr_x_yyyyy_xxyzz, tr_x_yyyyy_xyyy, tr_x_yyyyy_xyyyy, tr_x_yyyyy_xyyyz, tr_x_yyyyy_xyyz, tr_x_yyyyy_xyyzz, tr_x_yyyyy_xyzz, tr_x_yyyyy_xyzzz, tr_x_yyyyy_yyyy, tr_x_yyyyy_yyyyy, tr_x_yyyyy_yyyyz, tr_x_yyyyy_yyyz, tr_x_yyyyy_yyyzz, tr_x_yyyyy_yyzz, tr_x_yyyyy_yyzzz, tr_x_yyyyy_yzzz, tr_x_yyyyy_yzzzz, tr_x_yyyyy_zzzzz, ts_yyyyy_xxxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyzz, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzzz, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyy_xxxxx[i] = 4.0 * tr_x_xyyy_xxxxx[i] * fe_0 + tr_x_xyyyy_xxxxx[i] * pa_y[i];

        tr_x_xyyyyy_xxxxy[i] = 4.0 * tr_x_yyyyy_xxxy[i] * fe_0 + ts_yyyyy_xxxxy[i] * fe_0 + tr_x_yyyyy_xxxxy[i] * pa_x[i];

        tr_x_xyyyyy_xxxxz[i] = 4.0 * tr_x_xyyy_xxxxz[i] * fe_0 + tr_x_xyyyy_xxxxz[i] * pa_y[i];

        tr_x_xyyyyy_xxxyy[i] = 3.0 * tr_x_yyyyy_xxyy[i] * fe_0 + ts_yyyyy_xxxyy[i] * fe_0 + tr_x_yyyyy_xxxyy[i] * pa_x[i];

        tr_x_xyyyyy_xxxyz[i] = 3.0 * tr_x_yyyyy_xxyz[i] * fe_0 + ts_yyyyy_xxxyz[i] * fe_0 + tr_x_yyyyy_xxxyz[i] * pa_x[i];

        tr_x_xyyyyy_xxxzz[i] = 4.0 * tr_x_xyyy_xxxzz[i] * fe_0 + tr_x_xyyyy_xxxzz[i] * pa_y[i];

        tr_x_xyyyyy_xxyyy[i] = 2.0 * tr_x_yyyyy_xyyy[i] * fe_0 + ts_yyyyy_xxyyy[i] * fe_0 + tr_x_yyyyy_xxyyy[i] * pa_x[i];

        tr_x_xyyyyy_xxyyz[i] = 2.0 * tr_x_yyyyy_xyyz[i] * fe_0 + ts_yyyyy_xxyyz[i] * fe_0 + tr_x_yyyyy_xxyyz[i] * pa_x[i];

        tr_x_xyyyyy_xxyzz[i] = 2.0 * tr_x_yyyyy_xyzz[i] * fe_0 + ts_yyyyy_xxyzz[i] * fe_0 + tr_x_yyyyy_xxyzz[i] * pa_x[i];

        tr_x_xyyyyy_xxzzz[i] = 4.0 * tr_x_xyyy_xxzzz[i] * fe_0 + tr_x_xyyyy_xxzzz[i] * pa_y[i];

        tr_x_xyyyyy_xyyyy[i] = tr_x_yyyyy_yyyy[i] * fe_0 + ts_yyyyy_xyyyy[i] * fe_0 + tr_x_yyyyy_xyyyy[i] * pa_x[i];

        tr_x_xyyyyy_xyyyz[i] = tr_x_yyyyy_yyyz[i] * fe_0 + ts_yyyyy_xyyyz[i] * fe_0 + tr_x_yyyyy_xyyyz[i] * pa_x[i];

        tr_x_xyyyyy_xyyzz[i] = tr_x_yyyyy_yyzz[i] * fe_0 + ts_yyyyy_xyyzz[i] * fe_0 + tr_x_yyyyy_xyyzz[i] * pa_x[i];

        tr_x_xyyyyy_xyzzz[i] = tr_x_yyyyy_yzzz[i] * fe_0 + ts_yyyyy_xyzzz[i] * fe_0 + tr_x_yyyyy_xyzzz[i] * pa_x[i];

        tr_x_xyyyyy_xzzzz[i] = 4.0 * tr_x_xyyy_xzzzz[i] * fe_0 + tr_x_xyyyy_xzzzz[i] * pa_y[i];

        tr_x_xyyyyy_yyyyy[i] = ts_yyyyy_yyyyy[i] * fe_0 + tr_x_yyyyy_yyyyy[i] * pa_x[i];

        tr_x_xyyyyy_yyyyz[i] = ts_yyyyy_yyyyz[i] * fe_0 + tr_x_yyyyy_yyyyz[i] * pa_x[i];

        tr_x_xyyyyy_yyyzz[i] = ts_yyyyy_yyyzz[i] * fe_0 + tr_x_yyyyy_yyyzz[i] * pa_x[i];

        tr_x_xyyyyy_yyzzz[i] = ts_yyyyy_yyzzz[i] * fe_0 + tr_x_yyyyy_yyzzz[i] * pa_x[i];

        tr_x_xyyyyy_yzzzz[i] = ts_yyyyy_yzzzz[i] * fe_0 + tr_x_yyyyy_yzzzz[i] * pa_x[i];

        tr_x_xyyyyy_zzzzz[i] = ts_yyyyy_zzzzz[i] * fe_0 + tr_x_yyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 336-357 components of targeted buffer : IH

    auto tr_x_xyyyyz_xxxxx = pbuffer.data(idx_dip_ih + 336);

    auto tr_x_xyyyyz_xxxxy = pbuffer.data(idx_dip_ih + 337);

    auto tr_x_xyyyyz_xxxxz = pbuffer.data(idx_dip_ih + 338);

    auto tr_x_xyyyyz_xxxyy = pbuffer.data(idx_dip_ih + 339);

    auto tr_x_xyyyyz_xxxyz = pbuffer.data(idx_dip_ih + 340);

    auto tr_x_xyyyyz_xxxzz = pbuffer.data(idx_dip_ih + 341);

    auto tr_x_xyyyyz_xxyyy = pbuffer.data(idx_dip_ih + 342);

    auto tr_x_xyyyyz_xxyyz = pbuffer.data(idx_dip_ih + 343);

    auto tr_x_xyyyyz_xxyzz = pbuffer.data(idx_dip_ih + 344);

    auto tr_x_xyyyyz_xxzzz = pbuffer.data(idx_dip_ih + 345);

    auto tr_x_xyyyyz_xyyyy = pbuffer.data(idx_dip_ih + 346);

    auto tr_x_xyyyyz_xyyyz = pbuffer.data(idx_dip_ih + 347);

    auto tr_x_xyyyyz_xyyzz = pbuffer.data(idx_dip_ih + 348);

    auto tr_x_xyyyyz_xyzzz = pbuffer.data(idx_dip_ih + 349);

    auto tr_x_xyyyyz_xzzzz = pbuffer.data(idx_dip_ih + 350);

    auto tr_x_xyyyyz_yyyyy = pbuffer.data(idx_dip_ih + 351);

    auto tr_x_xyyyyz_yyyyz = pbuffer.data(idx_dip_ih + 352);

    auto tr_x_xyyyyz_yyyzz = pbuffer.data(idx_dip_ih + 353);

    auto tr_x_xyyyyz_yyzzz = pbuffer.data(idx_dip_ih + 354);

    auto tr_x_xyyyyz_yzzzz = pbuffer.data(idx_dip_ih + 355);

    auto tr_x_xyyyyz_zzzzz = pbuffer.data(idx_dip_ih + 356);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyyy_xxxxx, tr_x_xyyyy_xxxxy, tr_x_xyyyy_xxxy, tr_x_xyyyy_xxxyy, tr_x_xyyyy_xxxyz, tr_x_xyyyy_xxyy, tr_x_xyyyy_xxyyy, tr_x_xyyyy_xxyyz, tr_x_xyyyy_xxyz, tr_x_xyyyy_xxyzz, tr_x_xyyyy_xyyy, tr_x_xyyyy_xyyyy, tr_x_xyyyy_xyyyz, tr_x_xyyyy_xyyz, tr_x_xyyyy_xyyzz, tr_x_xyyyy_xyzz, tr_x_xyyyy_xyzzz, tr_x_xyyyy_yyyyy, tr_x_xyyyyz_xxxxx, tr_x_xyyyyz_xxxxy, tr_x_xyyyyz_xxxxz, tr_x_xyyyyz_xxxyy, tr_x_xyyyyz_xxxyz, tr_x_xyyyyz_xxxzz, tr_x_xyyyyz_xxyyy, tr_x_xyyyyz_xxyyz, tr_x_xyyyyz_xxyzz, tr_x_xyyyyz_xxzzz, tr_x_xyyyyz_xyyyy, tr_x_xyyyyz_xyyyz, tr_x_xyyyyz_xyyzz, tr_x_xyyyyz_xyzzz, tr_x_xyyyyz_xzzzz, tr_x_xyyyyz_yyyyy, tr_x_xyyyyz_yyyyz, tr_x_xyyyyz_yyyzz, tr_x_xyyyyz_yyzzz, tr_x_xyyyyz_yzzzz, tr_x_xyyyyz_zzzzz, tr_x_xyyyz_xxxxz, tr_x_xyyyz_xxxzz, tr_x_xyyyz_xxzzz, tr_x_xyyyz_xzzzz, tr_x_xyyz_xxxxz, tr_x_xyyz_xxxzz, tr_x_xyyz_xxzzz, tr_x_xyyz_xzzzz, tr_x_yyyyz_yyyyz, tr_x_yyyyz_yyyzz, tr_x_yyyyz_yyzzz, tr_x_yyyyz_yzzzz, tr_x_yyyyz_zzzzz, ts_yyyyz_yyyyz, ts_yyyyz_yyyzz, ts_yyyyz_yyzzz, ts_yyyyz_yzzzz, ts_yyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyz_xxxxx[i] = tr_x_xyyyy_xxxxx[i] * pa_z[i];

        tr_x_xyyyyz_xxxxy[i] = tr_x_xyyyy_xxxxy[i] * pa_z[i];

        tr_x_xyyyyz_xxxxz[i] = 3.0 * tr_x_xyyz_xxxxz[i] * fe_0 + tr_x_xyyyz_xxxxz[i] * pa_y[i];

        tr_x_xyyyyz_xxxyy[i] = tr_x_xyyyy_xxxyy[i] * pa_z[i];

        tr_x_xyyyyz_xxxyz[i] = tr_x_xyyyy_xxxy[i] * fe_0 + tr_x_xyyyy_xxxyz[i] * pa_z[i];

        tr_x_xyyyyz_xxxzz[i] = 3.0 * tr_x_xyyz_xxxzz[i] * fe_0 + tr_x_xyyyz_xxxzz[i] * pa_y[i];

        tr_x_xyyyyz_xxyyy[i] = tr_x_xyyyy_xxyyy[i] * pa_z[i];

        tr_x_xyyyyz_xxyyz[i] = tr_x_xyyyy_xxyy[i] * fe_0 + tr_x_xyyyy_xxyyz[i] * pa_z[i];

        tr_x_xyyyyz_xxyzz[i] = 2.0 * tr_x_xyyyy_xxyz[i] * fe_0 + tr_x_xyyyy_xxyzz[i] * pa_z[i];

        tr_x_xyyyyz_xxzzz[i] = 3.0 * tr_x_xyyz_xxzzz[i] * fe_0 + tr_x_xyyyz_xxzzz[i] * pa_y[i];

        tr_x_xyyyyz_xyyyy[i] = tr_x_xyyyy_xyyyy[i] * pa_z[i];

        tr_x_xyyyyz_xyyyz[i] = tr_x_xyyyy_xyyy[i] * fe_0 + tr_x_xyyyy_xyyyz[i] * pa_z[i];

        tr_x_xyyyyz_xyyzz[i] = 2.0 * tr_x_xyyyy_xyyz[i] * fe_0 + tr_x_xyyyy_xyyzz[i] * pa_z[i];

        tr_x_xyyyyz_xyzzz[i] = 3.0 * tr_x_xyyyy_xyzz[i] * fe_0 + tr_x_xyyyy_xyzzz[i] * pa_z[i];

        tr_x_xyyyyz_xzzzz[i] = 3.0 * tr_x_xyyz_xzzzz[i] * fe_0 + tr_x_xyyyz_xzzzz[i] * pa_y[i];

        tr_x_xyyyyz_yyyyy[i] = tr_x_xyyyy_yyyyy[i] * pa_z[i];

        tr_x_xyyyyz_yyyyz[i] = ts_yyyyz_yyyyz[i] * fe_0 + tr_x_yyyyz_yyyyz[i] * pa_x[i];

        tr_x_xyyyyz_yyyzz[i] = ts_yyyyz_yyyzz[i] * fe_0 + tr_x_yyyyz_yyyzz[i] * pa_x[i];

        tr_x_xyyyyz_yyzzz[i] = ts_yyyyz_yyzzz[i] * fe_0 + tr_x_yyyyz_yyzzz[i] * pa_x[i];

        tr_x_xyyyyz_yzzzz[i] = ts_yyyyz_yzzzz[i] * fe_0 + tr_x_yyyyz_yzzzz[i] * pa_x[i];

        tr_x_xyyyyz_zzzzz[i] = ts_yyyyz_zzzzz[i] * fe_0 + tr_x_yyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 357-378 components of targeted buffer : IH

    auto tr_x_xyyyzz_xxxxx = pbuffer.data(idx_dip_ih + 357);

    auto tr_x_xyyyzz_xxxxy = pbuffer.data(idx_dip_ih + 358);

    auto tr_x_xyyyzz_xxxxz = pbuffer.data(idx_dip_ih + 359);

    auto tr_x_xyyyzz_xxxyy = pbuffer.data(idx_dip_ih + 360);

    auto tr_x_xyyyzz_xxxyz = pbuffer.data(idx_dip_ih + 361);

    auto tr_x_xyyyzz_xxxzz = pbuffer.data(idx_dip_ih + 362);

    auto tr_x_xyyyzz_xxyyy = pbuffer.data(idx_dip_ih + 363);

    auto tr_x_xyyyzz_xxyyz = pbuffer.data(idx_dip_ih + 364);

    auto tr_x_xyyyzz_xxyzz = pbuffer.data(idx_dip_ih + 365);

    auto tr_x_xyyyzz_xxzzz = pbuffer.data(idx_dip_ih + 366);

    auto tr_x_xyyyzz_xyyyy = pbuffer.data(idx_dip_ih + 367);

    auto tr_x_xyyyzz_xyyyz = pbuffer.data(idx_dip_ih + 368);

    auto tr_x_xyyyzz_xyyzz = pbuffer.data(idx_dip_ih + 369);

    auto tr_x_xyyyzz_xyzzz = pbuffer.data(idx_dip_ih + 370);

    auto tr_x_xyyyzz_xzzzz = pbuffer.data(idx_dip_ih + 371);

    auto tr_x_xyyyzz_yyyyy = pbuffer.data(idx_dip_ih + 372);

    auto tr_x_xyyyzz_yyyyz = pbuffer.data(idx_dip_ih + 373);

    auto tr_x_xyyyzz_yyyzz = pbuffer.data(idx_dip_ih + 374);

    auto tr_x_xyyyzz_yyzzz = pbuffer.data(idx_dip_ih + 375);

    auto tr_x_xyyyzz_yzzzz = pbuffer.data(idx_dip_ih + 376);

    auto tr_x_xyyyzz_zzzzz = pbuffer.data(idx_dip_ih + 377);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyy_xxxxy, tr_x_xyyy_xxxyy, tr_x_xyyy_xxyyy, tr_x_xyyy_xyyyy, tr_x_xyyyz_xxxxy, tr_x_xyyyz_xxxyy, tr_x_xyyyz_xxyyy, tr_x_xyyyz_xyyyy, tr_x_xyyyzz_xxxxx, tr_x_xyyyzz_xxxxy, tr_x_xyyyzz_xxxxz, tr_x_xyyyzz_xxxyy, tr_x_xyyyzz_xxxyz, tr_x_xyyyzz_xxxzz, tr_x_xyyyzz_xxyyy, tr_x_xyyyzz_xxyyz, tr_x_xyyyzz_xxyzz, tr_x_xyyyzz_xxzzz, tr_x_xyyyzz_xyyyy, tr_x_xyyyzz_xyyyz, tr_x_xyyyzz_xyyzz, tr_x_xyyyzz_xyzzz, tr_x_xyyyzz_xzzzz, tr_x_xyyyzz_yyyyy, tr_x_xyyyzz_yyyyz, tr_x_xyyyzz_yyyzz, tr_x_xyyyzz_yyzzz, tr_x_xyyyzz_yzzzz, tr_x_xyyyzz_zzzzz, tr_x_xyyzz_xxxxx, tr_x_xyyzz_xxxxz, tr_x_xyyzz_xxxzz, tr_x_xyyzz_xxzzz, tr_x_xyyzz_xzzzz, tr_x_xyzz_xxxxx, tr_x_xyzz_xxxxz, tr_x_xyzz_xxxzz, tr_x_xyzz_xxzzz, tr_x_xyzz_xzzzz, tr_x_yyyzz_xxxyz, tr_x_yyyzz_xxyyz, tr_x_yyyzz_xxyz, tr_x_yyyzz_xxyzz, tr_x_yyyzz_xyyyz, tr_x_yyyzz_xyyz, tr_x_yyyzz_xyyzz, tr_x_yyyzz_xyzz, tr_x_yyyzz_xyzzz, tr_x_yyyzz_yyyyy, tr_x_yyyzz_yyyyz, tr_x_yyyzz_yyyz, tr_x_yyyzz_yyyzz, tr_x_yyyzz_yyzz, tr_x_yyyzz_yyzzz, tr_x_yyyzz_yzzz, tr_x_yyyzz_yzzzz, tr_x_yyyzz_zzzzz, ts_yyyzz_xxxyz, ts_yyyzz_xxyyz, ts_yyyzz_xxyzz, ts_yyyzz_xyyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzzz, ts_yyyzz_yyyyy, ts_yyyzz_yyyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyzz_xxxxx[i] = 2.0 * tr_x_xyzz_xxxxx[i] * fe_0 + tr_x_xyyzz_xxxxx[i] * pa_y[i];

        tr_x_xyyyzz_xxxxy[i] = tr_x_xyyy_xxxxy[i] * fe_0 + tr_x_xyyyz_xxxxy[i] * pa_z[i];

        tr_x_xyyyzz_xxxxz[i] = 2.0 * tr_x_xyzz_xxxxz[i] * fe_0 + tr_x_xyyzz_xxxxz[i] * pa_y[i];

        tr_x_xyyyzz_xxxyy[i] = tr_x_xyyy_xxxyy[i] * fe_0 + tr_x_xyyyz_xxxyy[i] * pa_z[i];

        tr_x_xyyyzz_xxxyz[i] = 3.0 * tr_x_yyyzz_xxyz[i] * fe_0 + ts_yyyzz_xxxyz[i] * fe_0 + tr_x_yyyzz_xxxyz[i] * pa_x[i];

        tr_x_xyyyzz_xxxzz[i] = 2.0 * tr_x_xyzz_xxxzz[i] * fe_0 + tr_x_xyyzz_xxxzz[i] * pa_y[i];

        tr_x_xyyyzz_xxyyy[i] = tr_x_xyyy_xxyyy[i] * fe_0 + tr_x_xyyyz_xxyyy[i] * pa_z[i];

        tr_x_xyyyzz_xxyyz[i] = 2.0 * tr_x_yyyzz_xyyz[i] * fe_0 + ts_yyyzz_xxyyz[i] * fe_0 + tr_x_yyyzz_xxyyz[i] * pa_x[i];

        tr_x_xyyyzz_xxyzz[i] = 2.0 * tr_x_yyyzz_xyzz[i] * fe_0 + ts_yyyzz_xxyzz[i] * fe_0 + tr_x_yyyzz_xxyzz[i] * pa_x[i];

        tr_x_xyyyzz_xxzzz[i] = 2.0 * tr_x_xyzz_xxzzz[i] * fe_0 + tr_x_xyyzz_xxzzz[i] * pa_y[i];

        tr_x_xyyyzz_xyyyy[i] = tr_x_xyyy_xyyyy[i] * fe_0 + tr_x_xyyyz_xyyyy[i] * pa_z[i];

        tr_x_xyyyzz_xyyyz[i] = tr_x_yyyzz_yyyz[i] * fe_0 + ts_yyyzz_xyyyz[i] * fe_0 + tr_x_yyyzz_xyyyz[i] * pa_x[i];

        tr_x_xyyyzz_xyyzz[i] = tr_x_yyyzz_yyzz[i] * fe_0 + ts_yyyzz_xyyzz[i] * fe_0 + tr_x_yyyzz_xyyzz[i] * pa_x[i];

        tr_x_xyyyzz_xyzzz[i] = tr_x_yyyzz_yzzz[i] * fe_0 + ts_yyyzz_xyzzz[i] * fe_0 + tr_x_yyyzz_xyzzz[i] * pa_x[i];

        tr_x_xyyyzz_xzzzz[i] = 2.0 * tr_x_xyzz_xzzzz[i] * fe_0 + tr_x_xyyzz_xzzzz[i] * pa_y[i];

        tr_x_xyyyzz_yyyyy[i] = ts_yyyzz_yyyyy[i] * fe_0 + tr_x_yyyzz_yyyyy[i] * pa_x[i];

        tr_x_xyyyzz_yyyyz[i] = ts_yyyzz_yyyyz[i] * fe_0 + tr_x_yyyzz_yyyyz[i] * pa_x[i];

        tr_x_xyyyzz_yyyzz[i] = ts_yyyzz_yyyzz[i] * fe_0 + tr_x_yyyzz_yyyzz[i] * pa_x[i];

        tr_x_xyyyzz_yyzzz[i] = ts_yyyzz_yyzzz[i] * fe_0 + tr_x_yyyzz_yyzzz[i] * pa_x[i];

        tr_x_xyyyzz_yzzzz[i] = ts_yyyzz_yzzzz[i] * fe_0 + tr_x_yyyzz_yzzzz[i] * pa_x[i];

        tr_x_xyyyzz_zzzzz[i] = ts_yyyzz_zzzzz[i] * fe_0 + tr_x_yyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 378-399 components of targeted buffer : IH

    auto tr_x_xyyzzz_xxxxx = pbuffer.data(idx_dip_ih + 378);

    auto tr_x_xyyzzz_xxxxy = pbuffer.data(idx_dip_ih + 379);

    auto tr_x_xyyzzz_xxxxz = pbuffer.data(idx_dip_ih + 380);

    auto tr_x_xyyzzz_xxxyy = pbuffer.data(idx_dip_ih + 381);

    auto tr_x_xyyzzz_xxxyz = pbuffer.data(idx_dip_ih + 382);

    auto tr_x_xyyzzz_xxxzz = pbuffer.data(idx_dip_ih + 383);

    auto tr_x_xyyzzz_xxyyy = pbuffer.data(idx_dip_ih + 384);

    auto tr_x_xyyzzz_xxyyz = pbuffer.data(idx_dip_ih + 385);

    auto tr_x_xyyzzz_xxyzz = pbuffer.data(idx_dip_ih + 386);

    auto tr_x_xyyzzz_xxzzz = pbuffer.data(idx_dip_ih + 387);

    auto tr_x_xyyzzz_xyyyy = pbuffer.data(idx_dip_ih + 388);

    auto tr_x_xyyzzz_xyyyz = pbuffer.data(idx_dip_ih + 389);

    auto tr_x_xyyzzz_xyyzz = pbuffer.data(idx_dip_ih + 390);

    auto tr_x_xyyzzz_xyzzz = pbuffer.data(idx_dip_ih + 391);

    auto tr_x_xyyzzz_xzzzz = pbuffer.data(idx_dip_ih + 392);

    auto tr_x_xyyzzz_yyyyy = pbuffer.data(idx_dip_ih + 393);

    auto tr_x_xyyzzz_yyyyz = pbuffer.data(idx_dip_ih + 394);

    auto tr_x_xyyzzz_yyyzz = pbuffer.data(idx_dip_ih + 395);

    auto tr_x_xyyzzz_yyzzz = pbuffer.data(idx_dip_ih + 396);

    auto tr_x_xyyzzz_yzzzz = pbuffer.data(idx_dip_ih + 397);

    auto tr_x_xyyzzz_zzzzz = pbuffer.data(idx_dip_ih + 398);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyz_xxxxy, tr_x_xyyz_xxxyy, tr_x_xyyz_xxyyy, tr_x_xyyz_xyyyy, tr_x_xyyzz_xxxxy, tr_x_xyyzz_xxxyy, tr_x_xyyzz_xxyyy, tr_x_xyyzz_xyyyy, tr_x_xyyzzz_xxxxx, tr_x_xyyzzz_xxxxy, tr_x_xyyzzz_xxxxz, tr_x_xyyzzz_xxxyy, tr_x_xyyzzz_xxxyz, tr_x_xyyzzz_xxxzz, tr_x_xyyzzz_xxyyy, tr_x_xyyzzz_xxyyz, tr_x_xyyzzz_xxyzz, tr_x_xyyzzz_xxzzz, tr_x_xyyzzz_xyyyy, tr_x_xyyzzz_xyyyz, tr_x_xyyzzz_xyyzz, tr_x_xyyzzz_xyzzz, tr_x_xyyzzz_xzzzz, tr_x_xyyzzz_yyyyy, tr_x_xyyzzz_yyyyz, tr_x_xyyzzz_yyyzz, tr_x_xyyzzz_yyzzz, tr_x_xyyzzz_yzzzz, tr_x_xyyzzz_zzzzz, tr_x_xyzzz_xxxxx, tr_x_xyzzz_xxxxz, tr_x_xyzzz_xxxzz, tr_x_xyzzz_xxzzz, tr_x_xyzzz_xzzzz, tr_x_xzzz_xxxxx, tr_x_xzzz_xxxxz, tr_x_xzzz_xxxzz, tr_x_xzzz_xxzzz, tr_x_xzzz_xzzzz, tr_x_yyzzz_xxxyz, tr_x_yyzzz_xxyyz, tr_x_yyzzz_xxyz, tr_x_yyzzz_xxyzz, tr_x_yyzzz_xyyyz, tr_x_yyzzz_xyyz, tr_x_yyzzz_xyyzz, tr_x_yyzzz_xyzz, tr_x_yyzzz_xyzzz, tr_x_yyzzz_yyyyy, tr_x_yyzzz_yyyyz, tr_x_yyzzz_yyyz, tr_x_yyzzz_yyyzz, tr_x_yyzzz_yyzz, tr_x_yyzzz_yyzzz, tr_x_yyzzz_yzzz, tr_x_yyzzz_yzzzz, tr_x_yyzzz_zzzzz, ts_yyzzz_xxxyz, ts_yyzzz_xxyyz, ts_yyzzz_xxyzz, ts_yyzzz_xyyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzzz, ts_yyzzz_yyyyy, ts_yyzzz_yyyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzzz_xxxxx[i] = tr_x_xzzz_xxxxx[i] * fe_0 + tr_x_xyzzz_xxxxx[i] * pa_y[i];

        tr_x_xyyzzz_xxxxy[i] = 2.0 * tr_x_xyyz_xxxxy[i] * fe_0 + tr_x_xyyzz_xxxxy[i] * pa_z[i];

        tr_x_xyyzzz_xxxxz[i] = tr_x_xzzz_xxxxz[i] * fe_0 + tr_x_xyzzz_xxxxz[i] * pa_y[i];

        tr_x_xyyzzz_xxxyy[i] = 2.0 * tr_x_xyyz_xxxyy[i] * fe_0 + tr_x_xyyzz_xxxyy[i] * pa_z[i];

        tr_x_xyyzzz_xxxyz[i] = 3.0 * tr_x_yyzzz_xxyz[i] * fe_0 + ts_yyzzz_xxxyz[i] * fe_0 + tr_x_yyzzz_xxxyz[i] * pa_x[i];

        tr_x_xyyzzz_xxxzz[i] = tr_x_xzzz_xxxzz[i] * fe_0 + tr_x_xyzzz_xxxzz[i] * pa_y[i];

        tr_x_xyyzzz_xxyyy[i] = 2.0 * tr_x_xyyz_xxyyy[i] * fe_0 + tr_x_xyyzz_xxyyy[i] * pa_z[i];

        tr_x_xyyzzz_xxyyz[i] = 2.0 * tr_x_yyzzz_xyyz[i] * fe_0 + ts_yyzzz_xxyyz[i] * fe_0 + tr_x_yyzzz_xxyyz[i] * pa_x[i];

        tr_x_xyyzzz_xxyzz[i] = 2.0 * tr_x_yyzzz_xyzz[i] * fe_0 + ts_yyzzz_xxyzz[i] * fe_0 + tr_x_yyzzz_xxyzz[i] * pa_x[i];

        tr_x_xyyzzz_xxzzz[i] = tr_x_xzzz_xxzzz[i] * fe_0 + tr_x_xyzzz_xxzzz[i] * pa_y[i];

        tr_x_xyyzzz_xyyyy[i] = 2.0 * tr_x_xyyz_xyyyy[i] * fe_0 + tr_x_xyyzz_xyyyy[i] * pa_z[i];

        tr_x_xyyzzz_xyyyz[i] = tr_x_yyzzz_yyyz[i] * fe_0 + ts_yyzzz_xyyyz[i] * fe_0 + tr_x_yyzzz_xyyyz[i] * pa_x[i];

        tr_x_xyyzzz_xyyzz[i] = tr_x_yyzzz_yyzz[i] * fe_0 + ts_yyzzz_xyyzz[i] * fe_0 + tr_x_yyzzz_xyyzz[i] * pa_x[i];

        tr_x_xyyzzz_xyzzz[i] = tr_x_yyzzz_yzzz[i] * fe_0 + ts_yyzzz_xyzzz[i] * fe_0 + tr_x_yyzzz_xyzzz[i] * pa_x[i];

        tr_x_xyyzzz_xzzzz[i] = tr_x_xzzz_xzzzz[i] * fe_0 + tr_x_xyzzz_xzzzz[i] * pa_y[i];

        tr_x_xyyzzz_yyyyy[i] = ts_yyzzz_yyyyy[i] * fe_0 + tr_x_yyzzz_yyyyy[i] * pa_x[i];

        tr_x_xyyzzz_yyyyz[i] = ts_yyzzz_yyyyz[i] * fe_0 + tr_x_yyzzz_yyyyz[i] * pa_x[i];

        tr_x_xyyzzz_yyyzz[i] = ts_yyzzz_yyyzz[i] * fe_0 + tr_x_yyzzz_yyyzz[i] * pa_x[i];

        tr_x_xyyzzz_yyzzz[i] = ts_yyzzz_yyzzz[i] * fe_0 + tr_x_yyzzz_yyzzz[i] * pa_x[i];

        tr_x_xyyzzz_yzzzz[i] = ts_yyzzz_yzzzz[i] * fe_0 + tr_x_yyzzz_yzzzz[i] * pa_x[i];

        tr_x_xyyzzz_zzzzz[i] = ts_yyzzz_zzzzz[i] * fe_0 + tr_x_yyzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 399-420 components of targeted buffer : IH

    auto tr_x_xyzzzz_xxxxx = pbuffer.data(idx_dip_ih + 399);

    auto tr_x_xyzzzz_xxxxy = pbuffer.data(idx_dip_ih + 400);

    auto tr_x_xyzzzz_xxxxz = pbuffer.data(idx_dip_ih + 401);

    auto tr_x_xyzzzz_xxxyy = pbuffer.data(idx_dip_ih + 402);

    auto tr_x_xyzzzz_xxxyz = pbuffer.data(idx_dip_ih + 403);

    auto tr_x_xyzzzz_xxxzz = pbuffer.data(idx_dip_ih + 404);

    auto tr_x_xyzzzz_xxyyy = pbuffer.data(idx_dip_ih + 405);

    auto tr_x_xyzzzz_xxyyz = pbuffer.data(idx_dip_ih + 406);

    auto tr_x_xyzzzz_xxyzz = pbuffer.data(idx_dip_ih + 407);

    auto tr_x_xyzzzz_xxzzz = pbuffer.data(idx_dip_ih + 408);

    auto tr_x_xyzzzz_xyyyy = pbuffer.data(idx_dip_ih + 409);

    auto tr_x_xyzzzz_xyyyz = pbuffer.data(idx_dip_ih + 410);

    auto tr_x_xyzzzz_xyyzz = pbuffer.data(idx_dip_ih + 411);

    auto tr_x_xyzzzz_xyzzz = pbuffer.data(idx_dip_ih + 412);

    auto tr_x_xyzzzz_xzzzz = pbuffer.data(idx_dip_ih + 413);

    auto tr_x_xyzzzz_yyyyy = pbuffer.data(idx_dip_ih + 414);

    auto tr_x_xyzzzz_yyyyz = pbuffer.data(idx_dip_ih + 415);

    auto tr_x_xyzzzz_yyyzz = pbuffer.data(idx_dip_ih + 416);

    auto tr_x_xyzzzz_yyzzz = pbuffer.data(idx_dip_ih + 417);

    auto tr_x_xyzzzz_yzzzz = pbuffer.data(idx_dip_ih + 418);

    auto tr_x_xyzzzz_zzzzz = pbuffer.data(idx_dip_ih + 419);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyzzzz_xxxxx, tr_x_xyzzzz_xxxxy, tr_x_xyzzzz_xxxxz, tr_x_xyzzzz_xxxyy, tr_x_xyzzzz_xxxyz, tr_x_xyzzzz_xxxzz, tr_x_xyzzzz_xxyyy, tr_x_xyzzzz_xxyyz, tr_x_xyzzzz_xxyzz, tr_x_xyzzzz_xxzzz, tr_x_xyzzzz_xyyyy, tr_x_xyzzzz_xyyyz, tr_x_xyzzzz_xyyzz, tr_x_xyzzzz_xyzzz, tr_x_xyzzzz_xzzzz, tr_x_xyzzzz_yyyyy, tr_x_xyzzzz_yyyyz, tr_x_xyzzzz_yyyzz, tr_x_xyzzzz_yyzzz, tr_x_xyzzzz_yzzzz, tr_x_xyzzzz_zzzzz, tr_x_xzzzz_xxxx, tr_x_xzzzz_xxxxx, tr_x_xzzzz_xxxxy, tr_x_xzzzz_xxxxz, tr_x_xzzzz_xxxy, tr_x_xzzzz_xxxyy, tr_x_xzzzz_xxxyz, tr_x_xzzzz_xxxz, tr_x_xzzzz_xxxzz, tr_x_xzzzz_xxyy, tr_x_xzzzz_xxyyy, tr_x_xzzzz_xxyyz, tr_x_xzzzz_xxyz, tr_x_xzzzz_xxyzz, tr_x_xzzzz_xxzz, tr_x_xzzzz_xxzzz, tr_x_xzzzz_xyyy, tr_x_xzzzz_xyyyy, tr_x_xzzzz_xyyyz, tr_x_xzzzz_xyyz, tr_x_xzzzz_xyyzz, tr_x_xzzzz_xyzz, tr_x_xzzzz_xyzzz, tr_x_xzzzz_xzzz, tr_x_xzzzz_xzzzz, tr_x_xzzzz_zzzzz, tr_x_yzzzz_yyyyy, tr_x_yzzzz_yyyyz, tr_x_yzzzz_yyyzz, tr_x_yzzzz_yyzzz, tr_x_yzzzz_yzzzz, ts_yzzzz_yyyyy, ts_yzzzz_yyyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzzz_xxxxx[i] = tr_x_xzzzz_xxxxx[i] * pa_y[i];

        tr_x_xyzzzz_xxxxy[i] = tr_x_xzzzz_xxxx[i] * fe_0 + tr_x_xzzzz_xxxxy[i] * pa_y[i];

        tr_x_xyzzzz_xxxxz[i] = tr_x_xzzzz_xxxxz[i] * pa_y[i];

        tr_x_xyzzzz_xxxyy[i] = 2.0 * tr_x_xzzzz_xxxy[i] * fe_0 + tr_x_xzzzz_xxxyy[i] * pa_y[i];

        tr_x_xyzzzz_xxxyz[i] = tr_x_xzzzz_xxxz[i] * fe_0 + tr_x_xzzzz_xxxyz[i] * pa_y[i];

        tr_x_xyzzzz_xxxzz[i] = tr_x_xzzzz_xxxzz[i] * pa_y[i];

        tr_x_xyzzzz_xxyyy[i] = 3.0 * tr_x_xzzzz_xxyy[i] * fe_0 + tr_x_xzzzz_xxyyy[i] * pa_y[i];

        tr_x_xyzzzz_xxyyz[i] = 2.0 * tr_x_xzzzz_xxyz[i] * fe_0 + tr_x_xzzzz_xxyyz[i] * pa_y[i];

        tr_x_xyzzzz_xxyzz[i] = tr_x_xzzzz_xxzz[i] * fe_0 + tr_x_xzzzz_xxyzz[i] * pa_y[i];

        tr_x_xyzzzz_xxzzz[i] = tr_x_xzzzz_xxzzz[i] * pa_y[i];

        tr_x_xyzzzz_xyyyy[i] = 4.0 * tr_x_xzzzz_xyyy[i] * fe_0 + tr_x_xzzzz_xyyyy[i] * pa_y[i];

        tr_x_xyzzzz_xyyyz[i] = 3.0 * tr_x_xzzzz_xyyz[i] * fe_0 + tr_x_xzzzz_xyyyz[i] * pa_y[i];

        tr_x_xyzzzz_xyyzz[i] = 2.0 * tr_x_xzzzz_xyzz[i] * fe_0 + tr_x_xzzzz_xyyzz[i] * pa_y[i];

        tr_x_xyzzzz_xyzzz[i] = tr_x_xzzzz_xzzz[i] * fe_0 + tr_x_xzzzz_xyzzz[i] * pa_y[i];

        tr_x_xyzzzz_xzzzz[i] = tr_x_xzzzz_xzzzz[i] * pa_y[i];

        tr_x_xyzzzz_yyyyy[i] = ts_yzzzz_yyyyy[i] * fe_0 + tr_x_yzzzz_yyyyy[i] * pa_x[i];

        tr_x_xyzzzz_yyyyz[i] = ts_yzzzz_yyyyz[i] * fe_0 + tr_x_yzzzz_yyyyz[i] * pa_x[i];

        tr_x_xyzzzz_yyyzz[i] = ts_yzzzz_yyyzz[i] * fe_0 + tr_x_yzzzz_yyyzz[i] * pa_x[i];

        tr_x_xyzzzz_yyzzz[i] = ts_yzzzz_yyzzz[i] * fe_0 + tr_x_yzzzz_yyzzz[i] * pa_x[i];

        tr_x_xyzzzz_yzzzz[i] = ts_yzzzz_yzzzz[i] * fe_0 + tr_x_yzzzz_yzzzz[i] * pa_x[i];

        tr_x_xyzzzz_zzzzz[i] = tr_x_xzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 420-441 components of targeted buffer : IH

    auto tr_x_xzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 420);

    auto tr_x_xzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 421);

    auto tr_x_xzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 422);

    auto tr_x_xzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 423);

    auto tr_x_xzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 424);

    auto tr_x_xzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 425);

    auto tr_x_xzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 426);

    auto tr_x_xzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 427);

    auto tr_x_xzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 428);

    auto tr_x_xzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 429);

    auto tr_x_xzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 430);

    auto tr_x_xzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 431);

    auto tr_x_xzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 432);

    auto tr_x_xzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 433);

    auto tr_x_xzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 434);

    auto tr_x_xzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 435);

    auto tr_x_xzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 436);

    auto tr_x_xzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 437);

    auto tr_x_xzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 438);

    auto tr_x_xzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 439);

    auto tr_x_xzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 440);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xzzz_xxxxx, tr_x_xzzz_xxxxy, tr_x_xzzz_xxxyy, tr_x_xzzz_xxyyy, tr_x_xzzz_xyyyy, tr_x_xzzzz_xxxxx, tr_x_xzzzz_xxxxy, tr_x_xzzzz_xxxyy, tr_x_xzzzz_xxyyy, tr_x_xzzzz_xyyyy, tr_x_xzzzzz_xxxxx, tr_x_xzzzzz_xxxxy, tr_x_xzzzzz_xxxxz, tr_x_xzzzzz_xxxyy, tr_x_xzzzzz_xxxyz, tr_x_xzzzzz_xxxzz, tr_x_xzzzzz_xxyyy, tr_x_xzzzzz_xxyyz, tr_x_xzzzzz_xxyzz, tr_x_xzzzzz_xxzzz, tr_x_xzzzzz_xyyyy, tr_x_xzzzzz_xyyyz, tr_x_xzzzzz_xyyzz, tr_x_xzzzzz_xyzzz, tr_x_xzzzzz_xzzzz, tr_x_xzzzzz_yyyyy, tr_x_xzzzzz_yyyyz, tr_x_xzzzzz_yyyzz, tr_x_xzzzzz_yyzzz, tr_x_xzzzzz_yzzzz, tr_x_xzzzzz_zzzzz, tr_x_zzzzz_xxxxz, tr_x_zzzzz_xxxyz, tr_x_zzzzz_xxxz, tr_x_zzzzz_xxxzz, tr_x_zzzzz_xxyyz, tr_x_zzzzz_xxyz, tr_x_zzzzz_xxyzz, tr_x_zzzzz_xxzz, tr_x_zzzzz_xxzzz, tr_x_zzzzz_xyyyz, tr_x_zzzzz_xyyz, tr_x_zzzzz_xyyzz, tr_x_zzzzz_xyzz, tr_x_zzzzz_xyzzz, tr_x_zzzzz_xzzz, tr_x_zzzzz_xzzzz, tr_x_zzzzz_yyyyy, tr_x_zzzzz_yyyyz, tr_x_zzzzz_yyyz, tr_x_zzzzz_yyyzz, tr_x_zzzzz_yyzz, tr_x_zzzzz_yyzzz, tr_x_zzzzz_yzzz, tr_x_zzzzz_yzzzz, tr_x_zzzzz_zzzz, tr_x_zzzzz_zzzzz, ts_zzzzz_xxxxz, ts_zzzzz_xxxyz, ts_zzzzz_xxxzz, ts_zzzzz_xxyyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzzz_xxxxx[i] = 4.0 * tr_x_xzzz_xxxxx[i] * fe_0 + tr_x_xzzzz_xxxxx[i] * pa_z[i];

        tr_x_xzzzzz_xxxxy[i] = 4.0 * tr_x_xzzz_xxxxy[i] * fe_0 + tr_x_xzzzz_xxxxy[i] * pa_z[i];

        tr_x_xzzzzz_xxxxz[i] = 4.0 * tr_x_zzzzz_xxxz[i] * fe_0 + ts_zzzzz_xxxxz[i] * fe_0 + tr_x_zzzzz_xxxxz[i] * pa_x[i];

        tr_x_xzzzzz_xxxyy[i] = 4.0 * tr_x_xzzz_xxxyy[i] * fe_0 + tr_x_xzzzz_xxxyy[i] * pa_z[i];

        tr_x_xzzzzz_xxxyz[i] = 3.0 * tr_x_zzzzz_xxyz[i] * fe_0 + ts_zzzzz_xxxyz[i] * fe_0 + tr_x_zzzzz_xxxyz[i] * pa_x[i];

        tr_x_xzzzzz_xxxzz[i] = 3.0 * tr_x_zzzzz_xxzz[i] * fe_0 + ts_zzzzz_xxxzz[i] * fe_0 + tr_x_zzzzz_xxxzz[i] * pa_x[i];

        tr_x_xzzzzz_xxyyy[i] = 4.0 * tr_x_xzzz_xxyyy[i] * fe_0 + tr_x_xzzzz_xxyyy[i] * pa_z[i];

        tr_x_xzzzzz_xxyyz[i] = 2.0 * tr_x_zzzzz_xyyz[i] * fe_0 + ts_zzzzz_xxyyz[i] * fe_0 + tr_x_zzzzz_xxyyz[i] * pa_x[i];

        tr_x_xzzzzz_xxyzz[i] = 2.0 * tr_x_zzzzz_xyzz[i] * fe_0 + ts_zzzzz_xxyzz[i] * fe_0 + tr_x_zzzzz_xxyzz[i] * pa_x[i];

        tr_x_xzzzzz_xxzzz[i] = 2.0 * tr_x_zzzzz_xzzz[i] * fe_0 + ts_zzzzz_xxzzz[i] * fe_0 + tr_x_zzzzz_xxzzz[i] * pa_x[i];

        tr_x_xzzzzz_xyyyy[i] = 4.0 * tr_x_xzzz_xyyyy[i] * fe_0 + tr_x_xzzzz_xyyyy[i] * pa_z[i];

        tr_x_xzzzzz_xyyyz[i] = tr_x_zzzzz_yyyz[i] * fe_0 + ts_zzzzz_xyyyz[i] * fe_0 + tr_x_zzzzz_xyyyz[i] * pa_x[i];

        tr_x_xzzzzz_xyyzz[i] = tr_x_zzzzz_yyzz[i] * fe_0 + ts_zzzzz_xyyzz[i] * fe_0 + tr_x_zzzzz_xyyzz[i] * pa_x[i];

        tr_x_xzzzzz_xyzzz[i] = tr_x_zzzzz_yzzz[i] * fe_0 + ts_zzzzz_xyzzz[i] * fe_0 + tr_x_zzzzz_xyzzz[i] * pa_x[i];

        tr_x_xzzzzz_xzzzz[i] = tr_x_zzzzz_zzzz[i] * fe_0 + ts_zzzzz_xzzzz[i] * fe_0 + tr_x_zzzzz_xzzzz[i] * pa_x[i];

        tr_x_xzzzzz_yyyyy[i] = ts_zzzzz_yyyyy[i] * fe_0 + tr_x_zzzzz_yyyyy[i] * pa_x[i];

        tr_x_xzzzzz_yyyyz[i] = ts_zzzzz_yyyyz[i] * fe_0 + tr_x_zzzzz_yyyyz[i] * pa_x[i];

        tr_x_xzzzzz_yyyzz[i] = ts_zzzzz_yyyzz[i] * fe_0 + tr_x_zzzzz_yyyzz[i] * pa_x[i];

        tr_x_xzzzzz_yyzzz[i] = ts_zzzzz_yyzzz[i] * fe_0 + tr_x_zzzzz_yyzzz[i] * pa_x[i];

        tr_x_xzzzzz_yzzzz[i] = ts_zzzzz_yzzzz[i] * fe_0 + tr_x_zzzzz_yzzzz[i] * pa_x[i];

        tr_x_xzzzzz_zzzzz[i] = ts_zzzzz_zzzzz[i] * fe_0 + tr_x_zzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 441-462 components of targeted buffer : IH

    auto tr_x_yyyyyy_xxxxx = pbuffer.data(idx_dip_ih + 441);

    auto tr_x_yyyyyy_xxxxy = pbuffer.data(idx_dip_ih + 442);

    auto tr_x_yyyyyy_xxxxz = pbuffer.data(idx_dip_ih + 443);

    auto tr_x_yyyyyy_xxxyy = pbuffer.data(idx_dip_ih + 444);

    auto tr_x_yyyyyy_xxxyz = pbuffer.data(idx_dip_ih + 445);

    auto tr_x_yyyyyy_xxxzz = pbuffer.data(idx_dip_ih + 446);

    auto tr_x_yyyyyy_xxyyy = pbuffer.data(idx_dip_ih + 447);

    auto tr_x_yyyyyy_xxyyz = pbuffer.data(idx_dip_ih + 448);

    auto tr_x_yyyyyy_xxyzz = pbuffer.data(idx_dip_ih + 449);

    auto tr_x_yyyyyy_xxzzz = pbuffer.data(idx_dip_ih + 450);

    auto tr_x_yyyyyy_xyyyy = pbuffer.data(idx_dip_ih + 451);

    auto tr_x_yyyyyy_xyyyz = pbuffer.data(idx_dip_ih + 452);

    auto tr_x_yyyyyy_xyyzz = pbuffer.data(idx_dip_ih + 453);

    auto tr_x_yyyyyy_xyzzz = pbuffer.data(idx_dip_ih + 454);

    auto tr_x_yyyyyy_xzzzz = pbuffer.data(idx_dip_ih + 455);

    auto tr_x_yyyyyy_yyyyy = pbuffer.data(idx_dip_ih + 456);

    auto tr_x_yyyyyy_yyyyz = pbuffer.data(idx_dip_ih + 457);

    auto tr_x_yyyyyy_yyyzz = pbuffer.data(idx_dip_ih + 458);

    auto tr_x_yyyyyy_yyzzz = pbuffer.data(idx_dip_ih + 459);

    auto tr_x_yyyyyy_yzzzz = pbuffer.data(idx_dip_ih + 460);

    auto tr_x_yyyyyy_zzzzz = pbuffer.data(idx_dip_ih + 461);

    #pragma omp simd aligned(pa_y, tr_x_yyyy_xxxxx, tr_x_yyyy_xxxxy, tr_x_yyyy_xxxxz, tr_x_yyyy_xxxyy, tr_x_yyyy_xxxyz, tr_x_yyyy_xxxzz, tr_x_yyyy_xxyyy, tr_x_yyyy_xxyyz, tr_x_yyyy_xxyzz, tr_x_yyyy_xxzzz, tr_x_yyyy_xyyyy, tr_x_yyyy_xyyyz, tr_x_yyyy_xyyzz, tr_x_yyyy_xyzzz, tr_x_yyyy_xzzzz, tr_x_yyyy_yyyyy, tr_x_yyyy_yyyyz, tr_x_yyyy_yyyzz, tr_x_yyyy_yyzzz, tr_x_yyyy_yzzzz, tr_x_yyyy_zzzzz, tr_x_yyyyy_xxxx, tr_x_yyyyy_xxxxx, tr_x_yyyyy_xxxxy, tr_x_yyyyy_xxxxz, tr_x_yyyyy_xxxy, tr_x_yyyyy_xxxyy, tr_x_yyyyy_xxxyz, tr_x_yyyyy_xxxz, tr_x_yyyyy_xxxzz, tr_x_yyyyy_xxyy, tr_x_yyyyy_xxyyy, tr_x_yyyyy_xxyyz, tr_x_yyyyy_xxyz, tr_x_yyyyy_xxyzz, tr_x_yyyyy_xxzz, tr_x_yyyyy_xxzzz, tr_x_yyyyy_xyyy, tr_x_yyyyy_xyyyy, tr_x_yyyyy_xyyyz, tr_x_yyyyy_xyyz, tr_x_yyyyy_xyyzz, tr_x_yyyyy_xyzz, tr_x_yyyyy_xyzzz, tr_x_yyyyy_xzzz, tr_x_yyyyy_xzzzz, tr_x_yyyyy_yyyy, tr_x_yyyyy_yyyyy, tr_x_yyyyy_yyyyz, tr_x_yyyyy_yyyz, tr_x_yyyyy_yyyzz, tr_x_yyyyy_yyzz, tr_x_yyyyy_yyzzz, tr_x_yyyyy_yzzz, tr_x_yyyyy_yzzzz, tr_x_yyyyy_zzzz, tr_x_yyyyy_zzzzz, tr_x_yyyyyy_xxxxx, tr_x_yyyyyy_xxxxy, tr_x_yyyyyy_xxxxz, tr_x_yyyyyy_xxxyy, tr_x_yyyyyy_xxxyz, tr_x_yyyyyy_xxxzz, tr_x_yyyyyy_xxyyy, tr_x_yyyyyy_xxyyz, tr_x_yyyyyy_xxyzz, tr_x_yyyyyy_xxzzz, tr_x_yyyyyy_xyyyy, tr_x_yyyyyy_xyyyz, tr_x_yyyyyy_xyyzz, tr_x_yyyyyy_xyzzz, tr_x_yyyyyy_xzzzz, tr_x_yyyyyy_yyyyy, tr_x_yyyyyy_yyyyz, tr_x_yyyyyy_yyyzz, tr_x_yyyyyy_yyzzz, tr_x_yyyyyy_yzzzz, tr_x_yyyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyy_xxxxx[i] = 5.0 * tr_x_yyyy_xxxxx[i] * fe_0 + tr_x_yyyyy_xxxxx[i] * pa_y[i];

        tr_x_yyyyyy_xxxxy[i] = 5.0 * tr_x_yyyy_xxxxy[i] * fe_0 + tr_x_yyyyy_xxxx[i] * fe_0 + tr_x_yyyyy_xxxxy[i] * pa_y[i];

        tr_x_yyyyyy_xxxxz[i] = 5.0 * tr_x_yyyy_xxxxz[i] * fe_0 + tr_x_yyyyy_xxxxz[i] * pa_y[i];

        tr_x_yyyyyy_xxxyy[i] = 5.0 * tr_x_yyyy_xxxyy[i] * fe_0 + 2.0 * tr_x_yyyyy_xxxy[i] * fe_0 + tr_x_yyyyy_xxxyy[i] * pa_y[i];

        tr_x_yyyyyy_xxxyz[i] = 5.0 * tr_x_yyyy_xxxyz[i] * fe_0 + tr_x_yyyyy_xxxz[i] * fe_0 + tr_x_yyyyy_xxxyz[i] * pa_y[i];

        tr_x_yyyyyy_xxxzz[i] = 5.0 * tr_x_yyyy_xxxzz[i] * fe_0 + tr_x_yyyyy_xxxzz[i] * pa_y[i];

        tr_x_yyyyyy_xxyyy[i] = 5.0 * tr_x_yyyy_xxyyy[i] * fe_0 + 3.0 * tr_x_yyyyy_xxyy[i] * fe_0 + tr_x_yyyyy_xxyyy[i] * pa_y[i];

        tr_x_yyyyyy_xxyyz[i] = 5.0 * tr_x_yyyy_xxyyz[i] * fe_0 + 2.0 * tr_x_yyyyy_xxyz[i] * fe_0 + tr_x_yyyyy_xxyyz[i] * pa_y[i];

        tr_x_yyyyyy_xxyzz[i] = 5.0 * tr_x_yyyy_xxyzz[i] * fe_0 + tr_x_yyyyy_xxzz[i] * fe_0 + tr_x_yyyyy_xxyzz[i] * pa_y[i];

        tr_x_yyyyyy_xxzzz[i] = 5.0 * tr_x_yyyy_xxzzz[i] * fe_0 + tr_x_yyyyy_xxzzz[i] * pa_y[i];

        tr_x_yyyyyy_xyyyy[i] = 5.0 * tr_x_yyyy_xyyyy[i] * fe_0 + 4.0 * tr_x_yyyyy_xyyy[i] * fe_0 + tr_x_yyyyy_xyyyy[i] * pa_y[i];

        tr_x_yyyyyy_xyyyz[i] = 5.0 * tr_x_yyyy_xyyyz[i] * fe_0 + 3.0 * tr_x_yyyyy_xyyz[i] * fe_0 + tr_x_yyyyy_xyyyz[i] * pa_y[i];

        tr_x_yyyyyy_xyyzz[i] = 5.0 * tr_x_yyyy_xyyzz[i] * fe_0 + 2.0 * tr_x_yyyyy_xyzz[i] * fe_0 + tr_x_yyyyy_xyyzz[i] * pa_y[i];

        tr_x_yyyyyy_xyzzz[i] = 5.0 * tr_x_yyyy_xyzzz[i] * fe_0 + tr_x_yyyyy_xzzz[i] * fe_0 + tr_x_yyyyy_xyzzz[i] * pa_y[i];

        tr_x_yyyyyy_xzzzz[i] = 5.0 * tr_x_yyyy_xzzzz[i] * fe_0 + tr_x_yyyyy_xzzzz[i] * pa_y[i];

        tr_x_yyyyyy_yyyyy[i] = 5.0 * tr_x_yyyy_yyyyy[i] * fe_0 + 5.0 * tr_x_yyyyy_yyyy[i] * fe_0 + tr_x_yyyyy_yyyyy[i] * pa_y[i];

        tr_x_yyyyyy_yyyyz[i] = 5.0 * tr_x_yyyy_yyyyz[i] * fe_0 + 4.0 * tr_x_yyyyy_yyyz[i] * fe_0 + tr_x_yyyyy_yyyyz[i] * pa_y[i];

        tr_x_yyyyyy_yyyzz[i] = 5.0 * tr_x_yyyy_yyyzz[i] * fe_0 + 3.0 * tr_x_yyyyy_yyzz[i] * fe_0 + tr_x_yyyyy_yyyzz[i] * pa_y[i];

        tr_x_yyyyyy_yyzzz[i] = 5.0 * tr_x_yyyy_yyzzz[i] * fe_0 + 2.0 * tr_x_yyyyy_yzzz[i] * fe_0 + tr_x_yyyyy_yyzzz[i] * pa_y[i];

        tr_x_yyyyyy_yzzzz[i] = 5.0 * tr_x_yyyy_yzzzz[i] * fe_0 + tr_x_yyyyy_zzzz[i] * fe_0 + tr_x_yyyyy_yzzzz[i] * pa_y[i];

        tr_x_yyyyyy_zzzzz[i] = 5.0 * tr_x_yyyy_zzzzz[i] * fe_0 + tr_x_yyyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 462-483 components of targeted buffer : IH

    auto tr_x_yyyyyz_xxxxx = pbuffer.data(idx_dip_ih + 462);

    auto tr_x_yyyyyz_xxxxy = pbuffer.data(idx_dip_ih + 463);

    auto tr_x_yyyyyz_xxxxz = pbuffer.data(idx_dip_ih + 464);

    auto tr_x_yyyyyz_xxxyy = pbuffer.data(idx_dip_ih + 465);

    auto tr_x_yyyyyz_xxxyz = pbuffer.data(idx_dip_ih + 466);

    auto tr_x_yyyyyz_xxxzz = pbuffer.data(idx_dip_ih + 467);

    auto tr_x_yyyyyz_xxyyy = pbuffer.data(idx_dip_ih + 468);

    auto tr_x_yyyyyz_xxyyz = pbuffer.data(idx_dip_ih + 469);

    auto tr_x_yyyyyz_xxyzz = pbuffer.data(idx_dip_ih + 470);

    auto tr_x_yyyyyz_xxzzz = pbuffer.data(idx_dip_ih + 471);

    auto tr_x_yyyyyz_xyyyy = pbuffer.data(idx_dip_ih + 472);

    auto tr_x_yyyyyz_xyyyz = pbuffer.data(idx_dip_ih + 473);

    auto tr_x_yyyyyz_xyyzz = pbuffer.data(idx_dip_ih + 474);

    auto tr_x_yyyyyz_xyzzz = pbuffer.data(idx_dip_ih + 475);

    auto tr_x_yyyyyz_xzzzz = pbuffer.data(idx_dip_ih + 476);

    auto tr_x_yyyyyz_yyyyy = pbuffer.data(idx_dip_ih + 477);

    auto tr_x_yyyyyz_yyyyz = pbuffer.data(idx_dip_ih + 478);

    auto tr_x_yyyyyz_yyyzz = pbuffer.data(idx_dip_ih + 479);

    auto tr_x_yyyyyz_yyzzz = pbuffer.data(idx_dip_ih + 480);

    auto tr_x_yyyyyz_yzzzz = pbuffer.data(idx_dip_ih + 481);

    auto tr_x_yyyyyz_zzzzz = pbuffer.data(idx_dip_ih + 482);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyyy_xxxxx, tr_x_yyyyy_xxxxy, tr_x_yyyyy_xxxy, tr_x_yyyyy_xxxyy, tr_x_yyyyy_xxxyz, tr_x_yyyyy_xxyy, tr_x_yyyyy_xxyyy, tr_x_yyyyy_xxyyz, tr_x_yyyyy_xxyz, tr_x_yyyyy_xxyzz, tr_x_yyyyy_xyyy, tr_x_yyyyy_xyyyy, tr_x_yyyyy_xyyyz, tr_x_yyyyy_xyyz, tr_x_yyyyy_xyyzz, tr_x_yyyyy_xyzz, tr_x_yyyyy_xyzzz, tr_x_yyyyy_yyyy, tr_x_yyyyy_yyyyy, tr_x_yyyyy_yyyyz, tr_x_yyyyy_yyyz, tr_x_yyyyy_yyyzz, tr_x_yyyyy_yyzz, tr_x_yyyyy_yyzzz, tr_x_yyyyy_yzzz, tr_x_yyyyy_yzzzz, tr_x_yyyyyz_xxxxx, tr_x_yyyyyz_xxxxy, tr_x_yyyyyz_xxxxz, tr_x_yyyyyz_xxxyy, tr_x_yyyyyz_xxxyz, tr_x_yyyyyz_xxxzz, tr_x_yyyyyz_xxyyy, tr_x_yyyyyz_xxyyz, tr_x_yyyyyz_xxyzz, tr_x_yyyyyz_xxzzz, tr_x_yyyyyz_xyyyy, tr_x_yyyyyz_xyyyz, tr_x_yyyyyz_xyyzz, tr_x_yyyyyz_xyzzz, tr_x_yyyyyz_xzzzz, tr_x_yyyyyz_yyyyy, tr_x_yyyyyz_yyyyz, tr_x_yyyyyz_yyyzz, tr_x_yyyyyz_yyzzz, tr_x_yyyyyz_yzzzz, tr_x_yyyyyz_zzzzz, tr_x_yyyyz_xxxxz, tr_x_yyyyz_xxxzz, tr_x_yyyyz_xxzzz, tr_x_yyyyz_xzzzz, tr_x_yyyyz_zzzzz, tr_x_yyyz_xxxxz, tr_x_yyyz_xxxzz, tr_x_yyyz_xxzzz, tr_x_yyyz_xzzzz, tr_x_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyz_xxxxx[i] = tr_x_yyyyy_xxxxx[i] * pa_z[i];

        tr_x_yyyyyz_xxxxy[i] = tr_x_yyyyy_xxxxy[i] * pa_z[i];

        tr_x_yyyyyz_xxxxz[i] = 4.0 * tr_x_yyyz_xxxxz[i] * fe_0 + tr_x_yyyyz_xxxxz[i] * pa_y[i];

        tr_x_yyyyyz_xxxyy[i] = tr_x_yyyyy_xxxyy[i] * pa_z[i];

        tr_x_yyyyyz_xxxyz[i] = tr_x_yyyyy_xxxy[i] * fe_0 + tr_x_yyyyy_xxxyz[i] * pa_z[i];

        tr_x_yyyyyz_xxxzz[i] = 4.0 * tr_x_yyyz_xxxzz[i] * fe_0 + tr_x_yyyyz_xxxzz[i] * pa_y[i];

        tr_x_yyyyyz_xxyyy[i] = tr_x_yyyyy_xxyyy[i] * pa_z[i];

        tr_x_yyyyyz_xxyyz[i] = tr_x_yyyyy_xxyy[i] * fe_0 + tr_x_yyyyy_xxyyz[i] * pa_z[i];

        tr_x_yyyyyz_xxyzz[i] = 2.0 * tr_x_yyyyy_xxyz[i] * fe_0 + tr_x_yyyyy_xxyzz[i] * pa_z[i];

        tr_x_yyyyyz_xxzzz[i] = 4.0 * tr_x_yyyz_xxzzz[i] * fe_0 + tr_x_yyyyz_xxzzz[i] * pa_y[i];

        tr_x_yyyyyz_xyyyy[i] = tr_x_yyyyy_xyyyy[i] * pa_z[i];

        tr_x_yyyyyz_xyyyz[i] = tr_x_yyyyy_xyyy[i] * fe_0 + tr_x_yyyyy_xyyyz[i] * pa_z[i];

        tr_x_yyyyyz_xyyzz[i] = 2.0 * tr_x_yyyyy_xyyz[i] * fe_0 + tr_x_yyyyy_xyyzz[i] * pa_z[i];

        tr_x_yyyyyz_xyzzz[i] = 3.0 * tr_x_yyyyy_xyzz[i] * fe_0 + tr_x_yyyyy_xyzzz[i] * pa_z[i];

        tr_x_yyyyyz_xzzzz[i] = 4.0 * tr_x_yyyz_xzzzz[i] * fe_0 + tr_x_yyyyz_xzzzz[i] * pa_y[i];

        tr_x_yyyyyz_yyyyy[i] = tr_x_yyyyy_yyyyy[i] * pa_z[i];

        tr_x_yyyyyz_yyyyz[i] = tr_x_yyyyy_yyyy[i] * fe_0 + tr_x_yyyyy_yyyyz[i] * pa_z[i];

        tr_x_yyyyyz_yyyzz[i] = 2.0 * tr_x_yyyyy_yyyz[i] * fe_0 + tr_x_yyyyy_yyyzz[i] * pa_z[i];

        tr_x_yyyyyz_yyzzz[i] = 3.0 * tr_x_yyyyy_yyzz[i] * fe_0 + tr_x_yyyyy_yyzzz[i] * pa_z[i];

        tr_x_yyyyyz_yzzzz[i] = 4.0 * tr_x_yyyyy_yzzz[i] * fe_0 + tr_x_yyyyy_yzzzz[i] * pa_z[i];

        tr_x_yyyyyz_zzzzz[i] = 4.0 * tr_x_yyyz_zzzzz[i] * fe_0 + tr_x_yyyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 483-504 components of targeted buffer : IH

    auto tr_x_yyyyzz_xxxxx = pbuffer.data(idx_dip_ih + 483);

    auto tr_x_yyyyzz_xxxxy = pbuffer.data(idx_dip_ih + 484);

    auto tr_x_yyyyzz_xxxxz = pbuffer.data(idx_dip_ih + 485);

    auto tr_x_yyyyzz_xxxyy = pbuffer.data(idx_dip_ih + 486);

    auto tr_x_yyyyzz_xxxyz = pbuffer.data(idx_dip_ih + 487);

    auto tr_x_yyyyzz_xxxzz = pbuffer.data(idx_dip_ih + 488);

    auto tr_x_yyyyzz_xxyyy = pbuffer.data(idx_dip_ih + 489);

    auto tr_x_yyyyzz_xxyyz = pbuffer.data(idx_dip_ih + 490);

    auto tr_x_yyyyzz_xxyzz = pbuffer.data(idx_dip_ih + 491);

    auto tr_x_yyyyzz_xxzzz = pbuffer.data(idx_dip_ih + 492);

    auto tr_x_yyyyzz_xyyyy = pbuffer.data(idx_dip_ih + 493);

    auto tr_x_yyyyzz_xyyyz = pbuffer.data(idx_dip_ih + 494);

    auto tr_x_yyyyzz_xyyzz = pbuffer.data(idx_dip_ih + 495);

    auto tr_x_yyyyzz_xyzzz = pbuffer.data(idx_dip_ih + 496);

    auto tr_x_yyyyzz_xzzzz = pbuffer.data(idx_dip_ih + 497);

    auto tr_x_yyyyzz_yyyyy = pbuffer.data(idx_dip_ih + 498);

    auto tr_x_yyyyzz_yyyyz = pbuffer.data(idx_dip_ih + 499);

    auto tr_x_yyyyzz_yyyzz = pbuffer.data(idx_dip_ih + 500);

    auto tr_x_yyyyzz_yyzzz = pbuffer.data(idx_dip_ih + 501);

    auto tr_x_yyyyzz_yzzzz = pbuffer.data(idx_dip_ih + 502);

    auto tr_x_yyyyzz_zzzzz = pbuffer.data(idx_dip_ih + 503);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyy_xxxxy, tr_x_yyyy_xxxyy, tr_x_yyyy_xxyyy, tr_x_yyyy_xyyyy, tr_x_yyyy_yyyyy, tr_x_yyyyz_xxxxy, tr_x_yyyyz_xxxyy, tr_x_yyyyz_xxyyy, tr_x_yyyyz_xyyyy, tr_x_yyyyz_yyyyy, tr_x_yyyyzz_xxxxx, tr_x_yyyyzz_xxxxy, tr_x_yyyyzz_xxxxz, tr_x_yyyyzz_xxxyy, tr_x_yyyyzz_xxxyz, tr_x_yyyyzz_xxxzz, tr_x_yyyyzz_xxyyy, tr_x_yyyyzz_xxyyz, tr_x_yyyyzz_xxyzz, tr_x_yyyyzz_xxzzz, tr_x_yyyyzz_xyyyy, tr_x_yyyyzz_xyyyz, tr_x_yyyyzz_xyyzz, tr_x_yyyyzz_xyzzz, tr_x_yyyyzz_xzzzz, tr_x_yyyyzz_yyyyy, tr_x_yyyyzz_yyyyz, tr_x_yyyyzz_yyyzz, tr_x_yyyyzz_yyzzz, tr_x_yyyyzz_yzzzz, tr_x_yyyyzz_zzzzz, tr_x_yyyzz_xxxxx, tr_x_yyyzz_xxxxz, tr_x_yyyzz_xxxyz, tr_x_yyyzz_xxxz, tr_x_yyyzz_xxxzz, tr_x_yyyzz_xxyyz, tr_x_yyyzz_xxyz, tr_x_yyyzz_xxyzz, tr_x_yyyzz_xxzz, tr_x_yyyzz_xxzzz, tr_x_yyyzz_xyyyz, tr_x_yyyzz_xyyz, tr_x_yyyzz_xyyzz, tr_x_yyyzz_xyzz, tr_x_yyyzz_xyzzz, tr_x_yyyzz_xzzz, tr_x_yyyzz_xzzzz, tr_x_yyyzz_yyyyz, tr_x_yyyzz_yyyz, tr_x_yyyzz_yyyzz, tr_x_yyyzz_yyzz, tr_x_yyyzz_yyzzz, tr_x_yyyzz_yzzz, tr_x_yyyzz_yzzzz, tr_x_yyyzz_zzzz, tr_x_yyyzz_zzzzz, tr_x_yyzz_xxxxx, tr_x_yyzz_xxxxz, tr_x_yyzz_xxxyz, tr_x_yyzz_xxxzz, tr_x_yyzz_xxyyz, tr_x_yyzz_xxyzz, tr_x_yyzz_xxzzz, tr_x_yyzz_xyyyz, tr_x_yyzz_xyyzz, tr_x_yyzz_xyzzz, tr_x_yyzz_xzzzz, tr_x_yyzz_yyyyz, tr_x_yyzz_yyyzz, tr_x_yyzz_yyzzz, tr_x_yyzz_yzzzz, tr_x_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyzz_xxxxx[i] = 3.0 * tr_x_yyzz_xxxxx[i] * fe_0 + tr_x_yyyzz_xxxxx[i] * pa_y[i];

        tr_x_yyyyzz_xxxxy[i] = tr_x_yyyy_xxxxy[i] * fe_0 + tr_x_yyyyz_xxxxy[i] * pa_z[i];

        tr_x_yyyyzz_xxxxz[i] = 3.0 * tr_x_yyzz_xxxxz[i] * fe_0 + tr_x_yyyzz_xxxxz[i] * pa_y[i];

        tr_x_yyyyzz_xxxyy[i] = tr_x_yyyy_xxxyy[i] * fe_0 + tr_x_yyyyz_xxxyy[i] * pa_z[i];

        tr_x_yyyyzz_xxxyz[i] = 3.0 * tr_x_yyzz_xxxyz[i] * fe_0 + tr_x_yyyzz_xxxz[i] * fe_0 + tr_x_yyyzz_xxxyz[i] * pa_y[i];

        tr_x_yyyyzz_xxxzz[i] = 3.0 * tr_x_yyzz_xxxzz[i] * fe_0 + tr_x_yyyzz_xxxzz[i] * pa_y[i];

        tr_x_yyyyzz_xxyyy[i] = tr_x_yyyy_xxyyy[i] * fe_0 + tr_x_yyyyz_xxyyy[i] * pa_z[i];

        tr_x_yyyyzz_xxyyz[i] = 3.0 * tr_x_yyzz_xxyyz[i] * fe_0 + 2.0 * tr_x_yyyzz_xxyz[i] * fe_0 + tr_x_yyyzz_xxyyz[i] * pa_y[i];

        tr_x_yyyyzz_xxyzz[i] = 3.0 * tr_x_yyzz_xxyzz[i] * fe_0 + tr_x_yyyzz_xxzz[i] * fe_0 + tr_x_yyyzz_xxyzz[i] * pa_y[i];

        tr_x_yyyyzz_xxzzz[i] = 3.0 * tr_x_yyzz_xxzzz[i] * fe_0 + tr_x_yyyzz_xxzzz[i] * pa_y[i];

        tr_x_yyyyzz_xyyyy[i] = tr_x_yyyy_xyyyy[i] * fe_0 + tr_x_yyyyz_xyyyy[i] * pa_z[i];

        tr_x_yyyyzz_xyyyz[i] = 3.0 * tr_x_yyzz_xyyyz[i] * fe_0 + 3.0 * tr_x_yyyzz_xyyz[i] * fe_0 + tr_x_yyyzz_xyyyz[i] * pa_y[i];

        tr_x_yyyyzz_xyyzz[i] = 3.0 * tr_x_yyzz_xyyzz[i] * fe_0 + 2.0 * tr_x_yyyzz_xyzz[i] * fe_0 + tr_x_yyyzz_xyyzz[i] * pa_y[i];

        tr_x_yyyyzz_xyzzz[i] = 3.0 * tr_x_yyzz_xyzzz[i] * fe_0 + tr_x_yyyzz_xzzz[i] * fe_0 + tr_x_yyyzz_xyzzz[i] * pa_y[i];

        tr_x_yyyyzz_xzzzz[i] = 3.0 * tr_x_yyzz_xzzzz[i] * fe_0 + tr_x_yyyzz_xzzzz[i] * pa_y[i];

        tr_x_yyyyzz_yyyyy[i] = tr_x_yyyy_yyyyy[i] * fe_0 + tr_x_yyyyz_yyyyy[i] * pa_z[i];

        tr_x_yyyyzz_yyyyz[i] = 3.0 * tr_x_yyzz_yyyyz[i] * fe_0 + 4.0 * tr_x_yyyzz_yyyz[i] * fe_0 + tr_x_yyyzz_yyyyz[i] * pa_y[i];

        tr_x_yyyyzz_yyyzz[i] = 3.0 * tr_x_yyzz_yyyzz[i] * fe_0 + 3.0 * tr_x_yyyzz_yyzz[i] * fe_0 + tr_x_yyyzz_yyyzz[i] * pa_y[i];

        tr_x_yyyyzz_yyzzz[i] = 3.0 * tr_x_yyzz_yyzzz[i] * fe_0 + 2.0 * tr_x_yyyzz_yzzz[i] * fe_0 + tr_x_yyyzz_yyzzz[i] * pa_y[i];

        tr_x_yyyyzz_yzzzz[i] = 3.0 * tr_x_yyzz_yzzzz[i] * fe_0 + tr_x_yyyzz_zzzz[i] * fe_0 + tr_x_yyyzz_yzzzz[i] * pa_y[i];

        tr_x_yyyyzz_zzzzz[i] = 3.0 * tr_x_yyzz_zzzzz[i] * fe_0 + tr_x_yyyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 504-525 components of targeted buffer : IH

    auto tr_x_yyyzzz_xxxxx = pbuffer.data(idx_dip_ih + 504);

    auto tr_x_yyyzzz_xxxxy = pbuffer.data(idx_dip_ih + 505);

    auto tr_x_yyyzzz_xxxxz = pbuffer.data(idx_dip_ih + 506);

    auto tr_x_yyyzzz_xxxyy = pbuffer.data(idx_dip_ih + 507);

    auto tr_x_yyyzzz_xxxyz = pbuffer.data(idx_dip_ih + 508);

    auto tr_x_yyyzzz_xxxzz = pbuffer.data(idx_dip_ih + 509);

    auto tr_x_yyyzzz_xxyyy = pbuffer.data(idx_dip_ih + 510);

    auto tr_x_yyyzzz_xxyyz = pbuffer.data(idx_dip_ih + 511);

    auto tr_x_yyyzzz_xxyzz = pbuffer.data(idx_dip_ih + 512);

    auto tr_x_yyyzzz_xxzzz = pbuffer.data(idx_dip_ih + 513);

    auto tr_x_yyyzzz_xyyyy = pbuffer.data(idx_dip_ih + 514);

    auto tr_x_yyyzzz_xyyyz = pbuffer.data(idx_dip_ih + 515);

    auto tr_x_yyyzzz_xyyzz = pbuffer.data(idx_dip_ih + 516);

    auto tr_x_yyyzzz_xyzzz = pbuffer.data(idx_dip_ih + 517);

    auto tr_x_yyyzzz_xzzzz = pbuffer.data(idx_dip_ih + 518);

    auto tr_x_yyyzzz_yyyyy = pbuffer.data(idx_dip_ih + 519);

    auto tr_x_yyyzzz_yyyyz = pbuffer.data(idx_dip_ih + 520);

    auto tr_x_yyyzzz_yyyzz = pbuffer.data(idx_dip_ih + 521);

    auto tr_x_yyyzzz_yyzzz = pbuffer.data(idx_dip_ih + 522);

    auto tr_x_yyyzzz_yzzzz = pbuffer.data(idx_dip_ih + 523);

    auto tr_x_yyyzzz_zzzzz = pbuffer.data(idx_dip_ih + 524);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyz_xxxxy, tr_x_yyyz_xxxyy, tr_x_yyyz_xxyyy, tr_x_yyyz_xyyyy, tr_x_yyyz_yyyyy, tr_x_yyyzz_xxxxy, tr_x_yyyzz_xxxyy, tr_x_yyyzz_xxyyy, tr_x_yyyzz_xyyyy, tr_x_yyyzz_yyyyy, tr_x_yyyzzz_xxxxx, tr_x_yyyzzz_xxxxy, tr_x_yyyzzz_xxxxz, tr_x_yyyzzz_xxxyy, tr_x_yyyzzz_xxxyz, tr_x_yyyzzz_xxxzz, tr_x_yyyzzz_xxyyy, tr_x_yyyzzz_xxyyz, tr_x_yyyzzz_xxyzz, tr_x_yyyzzz_xxzzz, tr_x_yyyzzz_xyyyy, tr_x_yyyzzz_xyyyz, tr_x_yyyzzz_xyyzz, tr_x_yyyzzz_xyzzz, tr_x_yyyzzz_xzzzz, tr_x_yyyzzz_yyyyy, tr_x_yyyzzz_yyyyz, tr_x_yyyzzz_yyyzz, tr_x_yyyzzz_yyzzz, tr_x_yyyzzz_yzzzz, tr_x_yyyzzz_zzzzz, tr_x_yyzzz_xxxxx, tr_x_yyzzz_xxxxz, tr_x_yyzzz_xxxyz, tr_x_yyzzz_xxxz, tr_x_yyzzz_xxxzz, tr_x_yyzzz_xxyyz, tr_x_yyzzz_xxyz, tr_x_yyzzz_xxyzz, tr_x_yyzzz_xxzz, tr_x_yyzzz_xxzzz, tr_x_yyzzz_xyyyz, tr_x_yyzzz_xyyz, tr_x_yyzzz_xyyzz, tr_x_yyzzz_xyzz, tr_x_yyzzz_xyzzz, tr_x_yyzzz_xzzz, tr_x_yyzzz_xzzzz, tr_x_yyzzz_yyyyz, tr_x_yyzzz_yyyz, tr_x_yyzzz_yyyzz, tr_x_yyzzz_yyzz, tr_x_yyzzz_yyzzz, tr_x_yyzzz_yzzz, tr_x_yyzzz_yzzzz, tr_x_yyzzz_zzzz, tr_x_yyzzz_zzzzz, tr_x_yzzz_xxxxx, tr_x_yzzz_xxxxz, tr_x_yzzz_xxxyz, tr_x_yzzz_xxxzz, tr_x_yzzz_xxyyz, tr_x_yzzz_xxyzz, tr_x_yzzz_xxzzz, tr_x_yzzz_xyyyz, tr_x_yzzz_xyyzz, tr_x_yzzz_xyzzz, tr_x_yzzz_xzzzz, tr_x_yzzz_yyyyz, tr_x_yzzz_yyyzz, tr_x_yzzz_yyzzz, tr_x_yzzz_yzzzz, tr_x_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzzz_xxxxx[i] = 2.0 * tr_x_yzzz_xxxxx[i] * fe_0 + tr_x_yyzzz_xxxxx[i] * pa_y[i];

        tr_x_yyyzzz_xxxxy[i] = 2.0 * tr_x_yyyz_xxxxy[i] * fe_0 + tr_x_yyyzz_xxxxy[i] * pa_z[i];

        tr_x_yyyzzz_xxxxz[i] = 2.0 * tr_x_yzzz_xxxxz[i] * fe_0 + tr_x_yyzzz_xxxxz[i] * pa_y[i];

        tr_x_yyyzzz_xxxyy[i] = 2.0 * tr_x_yyyz_xxxyy[i] * fe_0 + tr_x_yyyzz_xxxyy[i] * pa_z[i];

        tr_x_yyyzzz_xxxyz[i] = 2.0 * tr_x_yzzz_xxxyz[i] * fe_0 + tr_x_yyzzz_xxxz[i] * fe_0 + tr_x_yyzzz_xxxyz[i] * pa_y[i];

        tr_x_yyyzzz_xxxzz[i] = 2.0 * tr_x_yzzz_xxxzz[i] * fe_0 + tr_x_yyzzz_xxxzz[i] * pa_y[i];

        tr_x_yyyzzz_xxyyy[i] = 2.0 * tr_x_yyyz_xxyyy[i] * fe_0 + tr_x_yyyzz_xxyyy[i] * pa_z[i];

        tr_x_yyyzzz_xxyyz[i] = 2.0 * tr_x_yzzz_xxyyz[i] * fe_0 + 2.0 * tr_x_yyzzz_xxyz[i] * fe_0 + tr_x_yyzzz_xxyyz[i] * pa_y[i];

        tr_x_yyyzzz_xxyzz[i] = 2.0 * tr_x_yzzz_xxyzz[i] * fe_0 + tr_x_yyzzz_xxzz[i] * fe_0 + tr_x_yyzzz_xxyzz[i] * pa_y[i];

        tr_x_yyyzzz_xxzzz[i] = 2.0 * tr_x_yzzz_xxzzz[i] * fe_0 + tr_x_yyzzz_xxzzz[i] * pa_y[i];

        tr_x_yyyzzz_xyyyy[i] = 2.0 * tr_x_yyyz_xyyyy[i] * fe_0 + tr_x_yyyzz_xyyyy[i] * pa_z[i];

        tr_x_yyyzzz_xyyyz[i] = 2.0 * tr_x_yzzz_xyyyz[i] * fe_0 + 3.0 * tr_x_yyzzz_xyyz[i] * fe_0 + tr_x_yyzzz_xyyyz[i] * pa_y[i];

        tr_x_yyyzzz_xyyzz[i] = 2.0 * tr_x_yzzz_xyyzz[i] * fe_0 + 2.0 * tr_x_yyzzz_xyzz[i] * fe_0 + tr_x_yyzzz_xyyzz[i] * pa_y[i];

        tr_x_yyyzzz_xyzzz[i] = 2.0 * tr_x_yzzz_xyzzz[i] * fe_0 + tr_x_yyzzz_xzzz[i] * fe_0 + tr_x_yyzzz_xyzzz[i] * pa_y[i];

        tr_x_yyyzzz_xzzzz[i] = 2.0 * tr_x_yzzz_xzzzz[i] * fe_0 + tr_x_yyzzz_xzzzz[i] * pa_y[i];

        tr_x_yyyzzz_yyyyy[i] = 2.0 * tr_x_yyyz_yyyyy[i] * fe_0 + tr_x_yyyzz_yyyyy[i] * pa_z[i];

        tr_x_yyyzzz_yyyyz[i] = 2.0 * tr_x_yzzz_yyyyz[i] * fe_0 + 4.0 * tr_x_yyzzz_yyyz[i] * fe_0 + tr_x_yyzzz_yyyyz[i] * pa_y[i];

        tr_x_yyyzzz_yyyzz[i] = 2.0 * tr_x_yzzz_yyyzz[i] * fe_0 + 3.0 * tr_x_yyzzz_yyzz[i] * fe_0 + tr_x_yyzzz_yyyzz[i] * pa_y[i];

        tr_x_yyyzzz_yyzzz[i] = 2.0 * tr_x_yzzz_yyzzz[i] * fe_0 + 2.0 * tr_x_yyzzz_yzzz[i] * fe_0 + tr_x_yyzzz_yyzzz[i] * pa_y[i];

        tr_x_yyyzzz_yzzzz[i] = 2.0 * tr_x_yzzz_yzzzz[i] * fe_0 + tr_x_yyzzz_zzzz[i] * fe_0 + tr_x_yyzzz_yzzzz[i] * pa_y[i];

        tr_x_yyyzzz_zzzzz[i] = 2.0 * tr_x_yzzz_zzzzz[i] * fe_0 + tr_x_yyzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 525-546 components of targeted buffer : IH

    auto tr_x_yyzzzz_xxxxx = pbuffer.data(idx_dip_ih + 525);

    auto tr_x_yyzzzz_xxxxy = pbuffer.data(idx_dip_ih + 526);

    auto tr_x_yyzzzz_xxxxz = pbuffer.data(idx_dip_ih + 527);

    auto tr_x_yyzzzz_xxxyy = pbuffer.data(idx_dip_ih + 528);

    auto tr_x_yyzzzz_xxxyz = pbuffer.data(idx_dip_ih + 529);

    auto tr_x_yyzzzz_xxxzz = pbuffer.data(idx_dip_ih + 530);

    auto tr_x_yyzzzz_xxyyy = pbuffer.data(idx_dip_ih + 531);

    auto tr_x_yyzzzz_xxyyz = pbuffer.data(idx_dip_ih + 532);

    auto tr_x_yyzzzz_xxyzz = pbuffer.data(idx_dip_ih + 533);

    auto tr_x_yyzzzz_xxzzz = pbuffer.data(idx_dip_ih + 534);

    auto tr_x_yyzzzz_xyyyy = pbuffer.data(idx_dip_ih + 535);

    auto tr_x_yyzzzz_xyyyz = pbuffer.data(idx_dip_ih + 536);

    auto tr_x_yyzzzz_xyyzz = pbuffer.data(idx_dip_ih + 537);

    auto tr_x_yyzzzz_xyzzz = pbuffer.data(idx_dip_ih + 538);

    auto tr_x_yyzzzz_xzzzz = pbuffer.data(idx_dip_ih + 539);

    auto tr_x_yyzzzz_yyyyy = pbuffer.data(idx_dip_ih + 540);

    auto tr_x_yyzzzz_yyyyz = pbuffer.data(idx_dip_ih + 541);

    auto tr_x_yyzzzz_yyyzz = pbuffer.data(idx_dip_ih + 542);

    auto tr_x_yyzzzz_yyzzz = pbuffer.data(idx_dip_ih + 543);

    auto tr_x_yyzzzz_yzzzz = pbuffer.data(idx_dip_ih + 544);

    auto tr_x_yyzzzz_zzzzz = pbuffer.data(idx_dip_ih + 545);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyzz_xxxxy, tr_x_yyzz_xxxyy, tr_x_yyzz_xxyyy, tr_x_yyzz_xyyyy, tr_x_yyzz_yyyyy, tr_x_yyzzz_xxxxy, tr_x_yyzzz_xxxyy, tr_x_yyzzz_xxyyy, tr_x_yyzzz_xyyyy, tr_x_yyzzz_yyyyy, tr_x_yyzzzz_xxxxx, tr_x_yyzzzz_xxxxy, tr_x_yyzzzz_xxxxz, tr_x_yyzzzz_xxxyy, tr_x_yyzzzz_xxxyz, tr_x_yyzzzz_xxxzz, tr_x_yyzzzz_xxyyy, tr_x_yyzzzz_xxyyz, tr_x_yyzzzz_xxyzz, tr_x_yyzzzz_xxzzz, tr_x_yyzzzz_xyyyy, tr_x_yyzzzz_xyyyz, tr_x_yyzzzz_xyyzz, tr_x_yyzzzz_xyzzz, tr_x_yyzzzz_xzzzz, tr_x_yyzzzz_yyyyy, tr_x_yyzzzz_yyyyz, tr_x_yyzzzz_yyyzz, tr_x_yyzzzz_yyzzz, tr_x_yyzzzz_yzzzz, tr_x_yyzzzz_zzzzz, tr_x_yzzzz_xxxxx, tr_x_yzzzz_xxxxz, tr_x_yzzzz_xxxyz, tr_x_yzzzz_xxxz, tr_x_yzzzz_xxxzz, tr_x_yzzzz_xxyyz, tr_x_yzzzz_xxyz, tr_x_yzzzz_xxyzz, tr_x_yzzzz_xxzz, tr_x_yzzzz_xxzzz, tr_x_yzzzz_xyyyz, tr_x_yzzzz_xyyz, tr_x_yzzzz_xyyzz, tr_x_yzzzz_xyzz, tr_x_yzzzz_xyzzz, tr_x_yzzzz_xzzz, tr_x_yzzzz_xzzzz, tr_x_yzzzz_yyyyz, tr_x_yzzzz_yyyz, tr_x_yzzzz_yyyzz, tr_x_yzzzz_yyzz, tr_x_yzzzz_yyzzz, tr_x_yzzzz_yzzz, tr_x_yzzzz_yzzzz, tr_x_yzzzz_zzzz, tr_x_yzzzz_zzzzz, tr_x_zzzz_xxxxx, tr_x_zzzz_xxxxz, tr_x_zzzz_xxxyz, tr_x_zzzz_xxxzz, tr_x_zzzz_xxyyz, tr_x_zzzz_xxyzz, tr_x_zzzz_xxzzz, tr_x_zzzz_xyyyz, tr_x_zzzz_xyyzz, tr_x_zzzz_xyzzz, tr_x_zzzz_xzzzz, tr_x_zzzz_yyyyz, tr_x_zzzz_yyyzz, tr_x_zzzz_yyzzz, tr_x_zzzz_yzzzz, tr_x_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzzz_xxxxx[i] = tr_x_zzzz_xxxxx[i] * fe_0 + tr_x_yzzzz_xxxxx[i] * pa_y[i];

        tr_x_yyzzzz_xxxxy[i] = 3.0 * tr_x_yyzz_xxxxy[i] * fe_0 + tr_x_yyzzz_xxxxy[i] * pa_z[i];

        tr_x_yyzzzz_xxxxz[i] = tr_x_zzzz_xxxxz[i] * fe_0 + tr_x_yzzzz_xxxxz[i] * pa_y[i];

        tr_x_yyzzzz_xxxyy[i] = 3.0 * tr_x_yyzz_xxxyy[i] * fe_0 + tr_x_yyzzz_xxxyy[i] * pa_z[i];

        tr_x_yyzzzz_xxxyz[i] = tr_x_zzzz_xxxyz[i] * fe_0 + tr_x_yzzzz_xxxz[i] * fe_0 + tr_x_yzzzz_xxxyz[i] * pa_y[i];

        tr_x_yyzzzz_xxxzz[i] = tr_x_zzzz_xxxzz[i] * fe_0 + tr_x_yzzzz_xxxzz[i] * pa_y[i];

        tr_x_yyzzzz_xxyyy[i] = 3.0 * tr_x_yyzz_xxyyy[i] * fe_0 + tr_x_yyzzz_xxyyy[i] * pa_z[i];

        tr_x_yyzzzz_xxyyz[i] = tr_x_zzzz_xxyyz[i] * fe_0 + 2.0 * tr_x_yzzzz_xxyz[i] * fe_0 + tr_x_yzzzz_xxyyz[i] * pa_y[i];

        tr_x_yyzzzz_xxyzz[i] = tr_x_zzzz_xxyzz[i] * fe_0 + tr_x_yzzzz_xxzz[i] * fe_0 + tr_x_yzzzz_xxyzz[i] * pa_y[i];

        tr_x_yyzzzz_xxzzz[i] = tr_x_zzzz_xxzzz[i] * fe_0 + tr_x_yzzzz_xxzzz[i] * pa_y[i];

        tr_x_yyzzzz_xyyyy[i] = 3.0 * tr_x_yyzz_xyyyy[i] * fe_0 + tr_x_yyzzz_xyyyy[i] * pa_z[i];

        tr_x_yyzzzz_xyyyz[i] = tr_x_zzzz_xyyyz[i] * fe_0 + 3.0 * tr_x_yzzzz_xyyz[i] * fe_0 + tr_x_yzzzz_xyyyz[i] * pa_y[i];

        tr_x_yyzzzz_xyyzz[i] = tr_x_zzzz_xyyzz[i] * fe_0 + 2.0 * tr_x_yzzzz_xyzz[i] * fe_0 + tr_x_yzzzz_xyyzz[i] * pa_y[i];

        tr_x_yyzzzz_xyzzz[i] = tr_x_zzzz_xyzzz[i] * fe_0 + tr_x_yzzzz_xzzz[i] * fe_0 + tr_x_yzzzz_xyzzz[i] * pa_y[i];

        tr_x_yyzzzz_xzzzz[i] = tr_x_zzzz_xzzzz[i] * fe_0 + tr_x_yzzzz_xzzzz[i] * pa_y[i];

        tr_x_yyzzzz_yyyyy[i] = 3.0 * tr_x_yyzz_yyyyy[i] * fe_0 + tr_x_yyzzz_yyyyy[i] * pa_z[i];

        tr_x_yyzzzz_yyyyz[i] = tr_x_zzzz_yyyyz[i] * fe_0 + 4.0 * tr_x_yzzzz_yyyz[i] * fe_0 + tr_x_yzzzz_yyyyz[i] * pa_y[i];

        tr_x_yyzzzz_yyyzz[i] = tr_x_zzzz_yyyzz[i] * fe_0 + 3.0 * tr_x_yzzzz_yyzz[i] * fe_0 + tr_x_yzzzz_yyyzz[i] * pa_y[i];

        tr_x_yyzzzz_yyzzz[i] = tr_x_zzzz_yyzzz[i] * fe_0 + 2.0 * tr_x_yzzzz_yzzz[i] * fe_0 + tr_x_yzzzz_yyzzz[i] * pa_y[i];

        tr_x_yyzzzz_yzzzz[i] = tr_x_zzzz_yzzzz[i] * fe_0 + tr_x_yzzzz_zzzz[i] * fe_0 + tr_x_yzzzz_yzzzz[i] * pa_y[i];

        tr_x_yyzzzz_zzzzz[i] = tr_x_zzzz_zzzzz[i] * fe_0 + tr_x_yzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 546-567 components of targeted buffer : IH

    auto tr_x_yzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 546);

    auto tr_x_yzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 547);

    auto tr_x_yzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 548);

    auto tr_x_yzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 549);

    auto tr_x_yzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 550);

    auto tr_x_yzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 551);

    auto tr_x_yzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 552);

    auto tr_x_yzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 553);

    auto tr_x_yzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 554);

    auto tr_x_yzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 555);

    auto tr_x_yzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 556);

    auto tr_x_yzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 557);

    auto tr_x_yzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 558);

    auto tr_x_yzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 559);

    auto tr_x_yzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 560);

    auto tr_x_yzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 561);

    auto tr_x_yzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 562);

    auto tr_x_yzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 563);

    auto tr_x_yzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 564);

    auto tr_x_yzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 565);

    auto tr_x_yzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 566);

    #pragma omp simd aligned(pa_y, tr_x_yzzzzz_xxxxx, tr_x_yzzzzz_xxxxy, tr_x_yzzzzz_xxxxz, tr_x_yzzzzz_xxxyy, tr_x_yzzzzz_xxxyz, tr_x_yzzzzz_xxxzz, tr_x_yzzzzz_xxyyy, tr_x_yzzzzz_xxyyz, tr_x_yzzzzz_xxyzz, tr_x_yzzzzz_xxzzz, tr_x_yzzzzz_xyyyy, tr_x_yzzzzz_xyyyz, tr_x_yzzzzz_xyyzz, tr_x_yzzzzz_xyzzz, tr_x_yzzzzz_xzzzz, tr_x_yzzzzz_yyyyy, tr_x_yzzzzz_yyyyz, tr_x_yzzzzz_yyyzz, tr_x_yzzzzz_yyzzz, tr_x_yzzzzz_yzzzz, tr_x_yzzzzz_zzzzz, tr_x_zzzzz_xxxx, tr_x_zzzzz_xxxxx, tr_x_zzzzz_xxxxy, tr_x_zzzzz_xxxxz, tr_x_zzzzz_xxxy, tr_x_zzzzz_xxxyy, tr_x_zzzzz_xxxyz, tr_x_zzzzz_xxxz, tr_x_zzzzz_xxxzz, tr_x_zzzzz_xxyy, tr_x_zzzzz_xxyyy, tr_x_zzzzz_xxyyz, tr_x_zzzzz_xxyz, tr_x_zzzzz_xxyzz, tr_x_zzzzz_xxzz, tr_x_zzzzz_xxzzz, tr_x_zzzzz_xyyy, tr_x_zzzzz_xyyyy, tr_x_zzzzz_xyyyz, tr_x_zzzzz_xyyz, tr_x_zzzzz_xyyzz, tr_x_zzzzz_xyzz, tr_x_zzzzz_xyzzz, tr_x_zzzzz_xzzz, tr_x_zzzzz_xzzzz, tr_x_zzzzz_yyyy, tr_x_zzzzz_yyyyy, tr_x_zzzzz_yyyyz, tr_x_zzzzz_yyyz, tr_x_zzzzz_yyyzz, tr_x_zzzzz_yyzz, tr_x_zzzzz_yyzzz, tr_x_zzzzz_yzzz, tr_x_zzzzz_yzzzz, tr_x_zzzzz_zzzz, tr_x_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzzz_xxxxx[i] = tr_x_zzzzz_xxxxx[i] * pa_y[i];

        tr_x_yzzzzz_xxxxy[i] = tr_x_zzzzz_xxxx[i] * fe_0 + tr_x_zzzzz_xxxxy[i] * pa_y[i];

        tr_x_yzzzzz_xxxxz[i] = tr_x_zzzzz_xxxxz[i] * pa_y[i];

        tr_x_yzzzzz_xxxyy[i] = 2.0 * tr_x_zzzzz_xxxy[i] * fe_0 + tr_x_zzzzz_xxxyy[i] * pa_y[i];

        tr_x_yzzzzz_xxxyz[i] = tr_x_zzzzz_xxxz[i] * fe_0 + tr_x_zzzzz_xxxyz[i] * pa_y[i];

        tr_x_yzzzzz_xxxzz[i] = tr_x_zzzzz_xxxzz[i] * pa_y[i];

        tr_x_yzzzzz_xxyyy[i] = 3.0 * tr_x_zzzzz_xxyy[i] * fe_0 + tr_x_zzzzz_xxyyy[i] * pa_y[i];

        tr_x_yzzzzz_xxyyz[i] = 2.0 * tr_x_zzzzz_xxyz[i] * fe_0 + tr_x_zzzzz_xxyyz[i] * pa_y[i];

        tr_x_yzzzzz_xxyzz[i] = tr_x_zzzzz_xxzz[i] * fe_0 + tr_x_zzzzz_xxyzz[i] * pa_y[i];

        tr_x_yzzzzz_xxzzz[i] = tr_x_zzzzz_xxzzz[i] * pa_y[i];

        tr_x_yzzzzz_xyyyy[i] = 4.0 * tr_x_zzzzz_xyyy[i] * fe_0 + tr_x_zzzzz_xyyyy[i] * pa_y[i];

        tr_x_yzzzzz_xyyyz[i] = 3.0 * tr_x_zzzzz_xyyz[i] * fe_0 + tr_x_zzzzz_xyyyz[i] * pa_y[i];

        tr_x_yzzzzz_xyyzz[i] = 2.0 * tr_x_zzzzz_xyzz[i] * fe_0 + tr_x_zzzzz_xyyzz[i] * pa_y[i];

        tr_x_yzzzzz_xyzzz[i] = tr_x_zzzzz_xzzz[i] * fe_0 + tr_x_zzzzz_xyzzz[i] * pa_y[i];

        tr_x_yzzzzz_xzzzz[i] = tr_x_zzzzz_xzzzz[i] * pa_y[i];

        tr_x_yzzzzz_yyyyy[i] = 5.0 * tr_x_zzzzz_yyyy[i] * fe_0 + tr_x_zzzzz_yyyyy[i] * pa_y[i];

        tr_x_yzzzzz_yyyyz[i] = 4.0 * tr_x_zzzzz_yyyz[i] * fe_0 + tr_x_zzzzz_yyyyz[i] * pa_y[i];

        tr_x_yzzzzz_yyyzz[i] = 3.0 * tr_x_zzzzz_yyzz[i] * fe_0 + tr_x_zzzzz_yyyzz[i] * pa_y[i];

        tr_x_yzzzzz_yyzzz[i] = 2.0 * tr_x_zzzzz_yzzz[i] * fe_0 + tr_x_zzzzz_yyzzz[i] * pa_y[i];

        tr_x_yzzzzz_yzzzz[i] = tr_x_zzzzz_zzzz[i] * fe_0 + tr_x_zzzzz_yzzzz[i] * pa_y[i];

        tr_x_yzzzzz_zzzzz[i] = tr_x_zzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 567-588 components of targeted buffer : IH

    auto tr_x_zzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 567);

    auto tr_x_zzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 568);

    auto tr_x_zzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 569);

    auto tr_x_zzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 570);

    auto tr_x_zzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 571);

    auto tr_x_zzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 572);

    auto tr_x_zzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 573);

    auto tr_x_zzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 574);

    auto tr_x_zzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 575);

    auto tr_x_zzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 576);

    auto tr_x_zzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 577);

    auto tr_x_zzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 578);

    auto tr_x_zzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 579);

    auto tr_x_zzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 580);

    auto tr_x_zzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 581);

    auto tr_x_zzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 582);

    auto tr_x_zzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 583);

    auto tr_x_zzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 584);

    auto tr_x_zzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 585);

    auto tr_x_zzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 586);

    auto tr_x_zzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 587);

    #pragma omp simd aligned(pa_z, tr_x_zzzz_xxxxx, tr_x_zzzz_xxxxy, tr_x_zzzz_xxxxz, tr_x_zzzz_xxxyy, tr_x_zzzz_xxxyz, tr_x_zzzz_xxxzz, tr_x_zzzz_xxyyy, tr_x_zzzz_xxyyz, tr_x_zzzz_xxyzz, tr_x_zzzz_xxzzz, tr_x_zzzz_xyyyy, tr_x_zzzz_xyyyz, tr_x_zzzz_xyyzz, tr_x_zzzz_xyzzz, tr_x_zzzz_xzzzz, tr_x_zzzz_yyyyy, tr_x_zzzz_yyyyz, tr_x_zzzz_yyyzz, tr_x_zzzz_yyzzz, tr_x_zzzz_yzzzz, tr_x_zzzz_zzzzz, tr_x_zzzzz_xxxx, tr_x_zzzzz_xxxxx, tr_x_zzzzz_xxxxy, tr_x_zzzzz_xxxxz, tr_x_zzzzz_xxxy, tr_x_zzzzz_xxxyy, tr_x_zzzzz_xxxyz, tr_x_zzzzz_xxxz, tr_x_zzzzz_xxxzz, tr_x_zzzzz_xxyy, tr_x_zzzzz_xxyyy, tr_x_zzzzz_xxyyz, tr_x_zzzzz_xxyz, tr_x_zzzzz_xxyzz, tr_x_zzzzz_xxzz, tr_x_zzzzz_xxzzz, tr_x_zzzzz_xyyy, tr_x_zzzzz_xyyyy, tr_x_zzzzz_xyyyz, tr_x_zzzzz_xyyz, tr_x_zzzzz_xyyzz, tr_x_zzzzz_xyzz, tr_x_zzzzz_xyzzz, tr_x_zzzzz_xzzz, tr_x_zzzzz_xzzzz, tr_x_zzzzz_yyyy, tr_x_zzzzz_yyyyy, tr_x_zzzzz_yyyyz, tr_x_zzzzz_yyyz, tr_x_zzzzz_yyyzz, tr_x_zzzzz_yyzz, tr_x_zzzzz_yyzzz, tr_x_zzzzz_yzzz, tr_x_zzzzz_yzzzz, tr_x_zzzzz_zzzz, tr_x_zzzzz_zzzzz, tr_x_zzzzzz_xxxxx, tr_x_zzzzzz_xxxxy, tr_x_zzzzzz_xxxxz, tr_x_zzzzzz_xxxyy, tr_x_zzzzzz_xxxyz, tr_x_zzzzzz_xxxzz, tr_x_zzzzzz_xxyyy, tr_x_zzzzzz_xxyyz, tr_x_zzzzzz_xxyzz, tr_x_zzzzzz_xxzzz, tr_x_zzzzzz_xyyyy, tr_x_zzzzzz_xyyyz, tr_x_zzzzzz_xyyzz, tr_x_zzzzzz_xyzzz, tr_x_zzzzzz_xzzzz, tr_x_zzzzzz_yyyyy, tr_x_zzzzzz_yyyyz, tr_x_zzzzzz_yyyzz, tr_x_zzzzzz_yyzzz, tr_x_zzzzzz_yzzzz, tr_x_zzzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzzz_xxxxx[i] = 5.0 * tr_x_zzzz_xxxxx[i] * fe_0 + tr_x_zzzzz_xxxxx[i] * pa_z[i];

        tr_x_zzzzzz_xxxxy[i] = 5.0 * tr_x_zzzz_xxxxy[i] * fe_0 + tr_x_zzzzz_xxxxy[i] * pa_z[i];

        tr_x_zzzzzz_xxxxz[i] = 5.0 * tr_x_zzzz_xxxxz[i] * fe_0 + tr_x_zzzzz_xxxx[i] * fe_0 + tr_x_zzzzz_xxxxz[i] * pa_z[i];

        tr_x_zzzzzz_xxxyy[i] = 5.0 * tr_x_zzzz_xxxyy[i] * fe_0 + tr_x_zzzzz_xxxyy[i] * pa_z[i];

        tr_x_zzzzzz_xxxyz[i] = 5.0 * tr_x_zzzz_xxxyz[i] * fe_0 + tr_x_zzzzz_xxxy[i] * fe_0 + tr_x_zzzzz_xxxyz[i] * pa_z[i];

        tr_x_zzzzzz_xxxzz[i] = 5.0 * tr_x_zzzz_xxxzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xxxz[i] * fe_0 + tr_x_zzzzz_xxxzz[i] * pa_z[i];

        tr_x_zzzzzz_xxyyy[i] = 5.0 * tr_x_zzzz_xxyyy[i] * fe_0 + tr_x_zzzzz_xxyyy[i] * pa_z[i];

        tr_x_zzzzzz_xxyyz[i] = 5.0 * tr_x_zzzz_xxyyz[i] * fe_0 + tr_x_zzzzz_xxyy[i] * fe_0 + tr_x_zzzzz_xxyyz[i] * pa_z[i];

        tr_x_zzzzzz_xxyzz[i] = 5.0 * tr_x_zzzz_xxyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xxyz[i] * fe_0 + tr_x_zzzzz_xxyzz[i] * pa_z[i];

        tr_x_zzzzzz_xxzzz[i] = 5.0 * tr_x_zzzz_xxzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_xxzz[i] * fe_0 + tr_x_zzzzz_xxzzz[i] * pa_z[i];

        tr_x_zzzzzz_xyyyy[i] = 5.0 * tr_x_zzzz_xyyyy[i] * fe_0 + tr_x_zzzzz_xyyyy[i] * pa_z[i];

        tr_x_zzzzzz_xyyyz[i] = 5.0 * tr_x_zzzz_xyyyz[i] * fe_0 + tr_x_zzzzz_xyyy[i] * fe_0 + tr_x_zzzzz_xyyyz[i] * pa_z[i];

        tr_x_zzzzzz_xyyzz[i] = 5.0 * tr_x_zzzz_xyyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xyyz[i] * fe_0 + tr_x_zzzzz_xyyzz[i] * pa_z[i];

        tr_x_zzzzzz_xyzzz[i] = 5.0 * tr_x_zzzz_xyzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_xyzz[i] * fe_0 + tr_x_zzzzz_xyzzz[i] * pa_z[i];

        tr_x_zzzzzz_xzzzz[i] = 5.0 * tr_x_zzzz_xzzzz[i] * fe_0 + 4.0 * tr_x_zzzzz_xzzz[i] * fe_0 + tr_x_zzzzz_xzzzz[i] * pa_z[i];

        tr_x_zzzzzz_yyyyy[i] = 5.0 * tr_x_zzzz_yyyyy[i] * fe_0 + tr_x_zzzzz_yyyyy[i] * pa_z[i];

        tr_x_zzzzzz_yyyyz[i] = 5.0 * tr_x_zzzz_yyyyz[i] * fe_0 + tr_x_zzzzz_yyyy[i] * fe_0 + tr_x_zzzzz_yyyyz[i] * pa_z[i];

        tr_x_zzzzzz_yyyzz[i] = 5.0 * tr_x_zzzz_yyyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_yyyz[i] * fe_0 + tr_x_zzzzz_yyyzz[i] * pa_z[i];

        tr_x_zzzzzz_yyzzz[i] = 5.0 * tr_x_zzzz_yyzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_yyzz[i] * fe_0 + tr_x_zzzzz_yyzzz[i] * pa_z[i];

        tr_x_zzzzzz_yzzzz[i] = 5.0 * tr_x_zzzz_yzzzz[i] * fe_0 + 4.0 * tr_x_zzzzz_yzzz[i] * fe_0 + tr_x_zzzzz_yzzzz[i] * pa_z[i];

        tr_x_zzzzzz_zzzzz[i] = 5.0 * tr_x_zzzz_zzzzz[i] * fe_0 + 5.0 * tr_x_zzzzz_zzzz[i] * fe_0 + tr_x_zzzzz_zzzzz[i] * pa_z[i];
    }

    // Set up 588-609 components of targeted buffer : IH

    auto tr_y_xxxxxx_xxxxx = pbuffer.data(idx_dip_ih + 588);

    auto tr_y_xxxxxx_xxxxy = pbuffer.data(idx_dip_ih + 589);

    auto tr_y_xxxxxx_xxxxz = pbuffer.data(idx_dip_ih + 590);

    auto tr_y_xxxxxx_xxxyy = pbuffer.data(idx_dip_ih + 591);

    auto tr_y_xxxxxx_xxxyz = pbuffer.data(idx_dip_ih + 592);

    auto tr_y_xxxxxx_xxxzz = pbuffer.data(idx_dip_ih + 593);

    auto tr_y_xxxxxx_xxyyy = pbuffer.data(idx_dip_ih + 594);

    auto tr_y_xxxxxx_xxyyz = pbuffer.data(idx_dip_ih + 595);

    auto tr_y_xxxxxx_xxyzz = pbuffer.data(idx_dip_ih + 596);

    auto tr_y_xxxxxx_xxzzz = pbuffer.data(idx_dip_ih + 597);

    auto tr_y_xxxxxx_xyyyy = pbuffer.data(idx_dip_ih + 598);

    auto tr_y_xxxxxx_xyyyz = pbuffer.data(idx_dip_ih + 599);

    auto tr_y_xxxxxx_xyyzz = pbuffer.data(idx_dip_ih + 600);

    auto tr_y_xxxxxx_xyzzz = pbuffer.data(idx_dip_ih + 601);

    auto tr_y_xxxxxx_xzzzz = pbuffer.data(idx_dip_ih + 602);

    auto tr_y_xxxxxx_yyyyy = pbuffer.data(idx_dip_ih + 603);

    auto tr_y_xxxxxx_yyyyz = pbuffer.data(idx_dip_ih + 604);

    auto tr_y_xxxxxx_yyyzz = pbuffer.data(idx_dip_ih + 605);

    auto tr_y_xxxxxx_yyzzz = pbuffer.data(idx_dip_ih + 606);

    auto tr_y_xxxxxx_yzzzz = pbuffer.data(idx_dip_ih + 607);

    auto tr_y_xxxxxx_zzzzz = pbuffer.data(idx_dip_ih + 608);

    #pragma omp simd aligned(pa_x, tr_y_xxxx_xxxxx, tr_y_xxxx_xxxxy, tr_y_xxxx_xxxxz, tr_y_xxxx_xxxyy, tr_y_xxxx_xxxyz, tr_y_xxxx_xxxzz, tr_y_xxxx_xxyyy, tr_y_xxxx_xxyyz, tr_y_xxxx_xxyzz, tr_y_xxxx_xxzzz, tr_y_xxxx_xyyyy, tr_y_xxxx_xyyyz, tr_y_xxxx_xyyzz, tr_y_xxxx_xyzzz, tr_y_xxxx_xzzzz, tr_y_xxxx_yyyyy, tr_y_xxxx_yyyyz, tr_y_xxxx_yyyzz, tr_y_xxxx_yyzzz, tr_y_xxxx_yzzzz, tr_y_xxxx_zzzzz, tr_y_xxxxx_xxxx, tr_y_xxxxx_xxxxx, tr_y_xxxxx_xxxxy, tr_y_xxxxx_xxxxz, tr_y_xxxxx_xxxy, tr_y_xxxxx_xxxyy, tr_y_xxxxx_xxxyz, tr_y_xxxxx_xxxz, tr_y_xxxxx_xxxzz, tr_y_xxxxx_xxyy, tr_y_xxxxx_xxyyy, tr_y_xxxxx_xxyyz, tr_y_xxxxx_xxyz, tr_y_xxxxx_xxyzz, tr_y_xxxxx_xxzz, tr_y_xxxxx_xxzzz, tr_y_xxxxx_xyyy, tr_y_xxxxx_xyyyy, tr_y_xxxxx_xyyyz, tr_y_xxxxx_xyyz, tr_y_xxxxx_xyyzz, tr_y_xxxxx_xyzz, tr_y_xxxxx_xyzzz, tr_y_xxxxx_xzzz, tr_y_xxxxx_xzzzz, tr_y_xxxxx_yyyy, tr_y_xxxxx_yyyyy, tr_y_xxxxx_yyyyz, tr_y_xxxxx_yyyz, tr_y_xxxxx_yyyzz, tr_y_xxxxx_yyzz, tr_y_xxxxx_yyzzz, tr_y_xxxxx_yzzz, tr_y_xxxxx_yzzzz, tr_y_xxxxx_zzzz, tr_y_xxxxx_zzzzz, tr_y_xxxxxx_xxxxx, tr_y_xxxxxx_xxxxy, tr_y_xxxxxx_xxxxz, tr_y_xxxxxx_xxxyy, tr_y_xxxxxx_xxxyz, tr_y_xxxxxx_xxxzz, tr_y_xxxxxx_xxyyy, tr_y_xxxxxx_xxyyz, tr_y_xxxxxx_xxyzz, tr_y_xxxxxx_xxzzz, tr_y_xxxxxx_xyyyy, tr_y_xxxxxx_xyyyz, tr_y_xxxxxx_xyyzz, tr_y_xxxxxx_xyzzz, tr_y_xxxxxx_xzzzz, tr_y_xxxxxx_yyyyy, tr_y_xxxxxx_yyyyz, tr_y_xxxxxx_yyyzz, tr_y_xxxxxx_yyzzz, tr_y_xxxxxx_yzzzz, tr_y_xxxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxx_xxxxx[i] = 5.0 * tr_y_xxxx_xxxxx[i] * fe_0 + 5.0 * tr_y_xxxxx_xxxx[i] * fe_0 + tr_y_xxxxx_xxxxx[i] * pa_x[i];

        tr_y_xxxxxx_xxxxy[i] = 5.0 * tr_y_xxxx_xxxxy[i] * fe_0 + 4.0 * tr_y_xxxxx_xxxy[i] * fe_0 + tr_y_xxxxx_xxxxy[i] * pa_x[i];

        tr_y_xxxxxx_xxxxz[i] = 5.0 * tr_y_xxxx_xxxxz[i] * fe_0 + 4.0 * tr_y_xxxxx_xxxz[i] * fe_0 + tr_y_xxxxx_xxxxz[i] * pa_x[i];

        tr_y_xxxxxx_xxxyy[i] = 5.0 * tr_y_xxxx_xxxyy[i] * fe_0 + 3.0 * tr_y_xxxxx_xxyy[i] * fe_0 + tr_y_xxxxx_xxxyy[i] * pa_x[i];

        tr_y_xxxxxx_xxxyz[i] = 5.0 * tr_y_xxxx_xxxyz[i] * fe_0 + 3.0 * tr_y_xxxxx_xxyz[i] * fe_0 + tr_y_xxxxx_xxxyz[i] * pa_x[i];

        tr_y_xxxxxx_xxxzz[i] = 5.0 * tr_y_xxxx_xxxzz[i] * fe_0 + 3.0 * tr_y_xxxxx_xxzz[i] * fe_0 + tr_y_xxxxx_xxxzz[i] * pa_x[i];

        tr_y_xxxxxx_xxyyy[i] = 5.0 * tr_y_xxxx_xxyyy[i] * fe_0 + 2.0 * tr_y_xxxxx_xyyy[i] * fe_0 + tr_y_xxxxx_xxyyy[i] * pa_x[i];

        tr_y_xxxxxx_xxyyz[i] = 5.0 * tr_y_xxxx_xxyyz[i] * fe_0 + 2.0 * tr_y_xxxxx_xyyz[i] * fe_0 + tr_y_xxxxx_xxyyz[i] * pa_x[i];

        tr_y_xxxxxx_xxyzz[i] = 5.0 * tr_y_xxxx_xxyzz[i] * fe_0 + 2.0 * tr_y_xxxxx_xyzz[i] * fe_0 + tr_y_xxxxx_xxyzz[i] * pa_x[i];

        tr_y_xxxxxx_xxzzz[i] = 5.0 * tr_y_xxxx_xxzzz[i] * fe_0 + 2.0 * tr_y_xxxxx_xzzz[i] * fe_0 + tr_y_xxxxx_xxzzz[i] * pa_x[i];

        tr_y_xxxxxx_xyyyy[i] = 5.0 * tr_y_xxxx_xyyyy[i] * fe_0 + tr_y_xxxxx_yyyy[i] * fe_0 + tr_y_xxxxx_xyyyy[i] * pa_x[i];

        tr_y_xxxxxx_xyyyz[i] = 5.0 * tr_y_xxxx_xyyyz[i] * fe_0 + tr_y_xxxxx_yyyz[i] * fe_0 + tr_y_xxxxx_xyyyz[i] * pa_x[i];

        tr_y_xxxxxx_xyyzz[i] = 5.0 * tr_y_xxxx_xyyzz[i] * fe_0 + tr_y_xxxxx_yyzz[i] * fe_0 + tr_y_xxxxx_xyyzz[i] * pa_x[i];

        tr_y_xxxxxx_xyzzz[i] = 5.0 * tr_y_xxxx_xyzzz[i] * fe_0 + tr_y_xxxxx_yzzz[i] * fe_0 + tr_y_xxxxx_xyzzz[i] * pa_x[i];

        tr_y_xxxxxx_xzzzz[i] = 5.0 * tr_y_xxxx_xzzzz[i] * fe_0 + tr_y_xxxxx_zzzz[i] * fe_0 + tr_y_xxxxx_xzzzz[i] * pa_x[i];

        tr_y_xxxxxx_yyyyy[i] = 5.0 * tr_y_xxxx_yyyyy[i] * fe_0 + tr_y_xxxxx_yyyyy[i] * pa_x[i];

        tr_y_xxxxxx_yyyyz[i] = 5.0 * tr_y_xxxx_yyyyz[i] * fe_0 + tr_y_xxxxx_yyyyz[i] * pa_x[i];

        tr_y_xxxxxx_yyyzz[i] = 5.0 * tr_y_xxxx_yyyzz[i] * fe_0 + tr_y_xxxxx_yyyzz[i] * pa_x[i];

        tr_y_xxxxxx_yyzzz[i] = 5.0 * tr_y_xxxx_yyzzz[i] * fe_0 + tr_y_xxxxx_yyzzz[i] * pa_x[i];

        tr_y_xxxxxx_yzzzz[i] = 5.0 * tr_y_xxxx_yzzzz[i] * fe_0 + tr_y_xxxxx_yzzzz[i] * pa_x[i];

        tr_y_xxxxxx_zzzzz[i] = 5.0 * tr_y_xxxx_zzzzz[i] * fe_0 + tr_y_xxxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 609-630 components of targeted buffer : IH

    auto tr_y_xxxxxy_xxxxx = pbuffer.data(idx_dip_ih + 609);

    auto tr_y_xxxxxy_xxxxy = pbuffer.data(idx_dip_ih + 610);

    auto tr_y_xxxxxy_xxxxz = pbuffer.data(idx_dip_ih + 611);

    auto tr_y_xxxxxy_xxxyy = pbuffer.data(idx_dip_ih + 612);

    auto tr_y_xxxxxy_xxxyz = pbuffer.data(idx_dip_ih + 613);

    auto tr_y_xxxxxy_xxxzz = pbuffer.data(idx_dip_ih + 614);

    auto tr_y_xxxxxy_xxyyy = pbuffer.data(idx_dip_ih + 615);

    auto tr_y_xxxxxy_xxyyz = pbuffer.data(idx_dip_ih + 616);

    auto tr_y_xxxxxy_xxyzz = pbuffer.data(idx_dip_ih + 617);

    auto tr_y_xxxxxy_xxzzz = pbuffer.data(idx_dip_ih + 618);

    auto tr_y_xxxxxy_xyyyy = pbuffer.data(idx_dip_ih + 619);

    auto tr_y_xxxxxy_xyyyz = pbuffer.data(idx_dip_ih + 620);

    auto tr_y_xxxxxy_xyyzz = pbuffer.data(idx_dip_ih + 621);

    auto tr_y_xxxxxy_xyzzz = pbuffer.data(idx_dip_ih + 622);

    auto tr_y_xxxxxy_xzzzz = pbuffer.data(idx_dip_ih + 623);

    auto tr_y_xxxxxy_yyyyy = pbuffer.data(idx_dip_ih + 624);

    auto tr_y_xxxxxy_yyyyz = pbuffer.data(idx_dip_ih + 625);

    auto tr_y_xxxxxy_yyyzz = pbuffer.data(idx_dip_ih + 626);

    auto tr_y_xxxxxy_yyzzz = pbuffer.data(idx_dip_ih + 627);

    auto tr_y_xxxxxy_yzzzz = pbuffer.data(idx_dip_ih + 628);

    auto tr_y_xxxxxy_zzzzz = pbuffer.data(idx_dip_ih + 629);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xxxxx_xxxxx, tr_y_xxxxx_xxxxz, tr_y_xxxxx_xxxzz, tr_y_xxxxx_xxzzz, tr_y_xxxxx_xzzzz, tr_y_xxxxxy_xxxxx, tr_y_xxxxxy_xxxxy, tr_y_xxxxxy_xxxxz, tr_y_xxxxxy_xxxyy, tr_y_xxxxxy_xxxyz, tr_y_xxxxxy_xxxzz, tr_y_xxxxxy_xxyyy, tr_y_xxxxxy_xxyyz, tr_y_xxxxxy_xxyzz, tr_y_xxxxxy_xxzzz, tr_y_xxxxxy_xyyyy, tr_y_xxxxxy_xyyyz, tr_y_xxxxxy_xyyzz, tr_y_xxxxxy_xyzzz, tr_y_xxxxxy_xzzzz, tr_y_xxxxxy_yyyyy, tr_y_xxxxxy_yyyyz, tr_y_xxxxxy_yyyzz, tr_y_xxxxxy_yyzzz, tr_y_xxxxxy_yzzzz, tr_y_xxxxxy_zzzzz, tr_y_xxxxy_xxxxy, tr_y_xxxxy_xxxy, tr_y_xxxxy_xxxyy, tr_y_xxxxy_xxxyz, tr_y_xxxxy_xxyy, tr_y_xxxxy_xxyyy, tr_y_xxxxy_xxyyz, tr_y_xxxxy_xxyz, tr_y_xxxxy_xxyzz, tr_y_xxxxy_xyyy, tr_y_xxxxy_xyyyy, tr_y_xxxxy_xyyyz, tr_y_xxxxy_xyyz, tr_y_xxxxy_xyyzz, tr_y_xxxxy_xyzz, tr_y_xxxxy_xyzzz, tr_y_xxxxy_yyyy, tr_y_xxxxy_yyyyy, tr_y_xxxxy_yyyyz, tr_y_xxxxy_yyyz, tr_y_xxxxy_yyyzz, tr_y_xxxxy_yyzz, tr_y_xxxxy_yyzzz, tr_y_xxxxy_yzzz, tr_y_xxxxy_yzzzz, tr_y_xxxxy_zzzzz, tr_y_xxxy_xxxxy, tr_y_xxxy_xxxyy, tr_y_xxxy_xxxyz, tr_y_xxxy_xxyyy, tr_y_xxxy_xxyyz, tr_y_xxxy_xxyzz, tr_y_xxxy_xyyyy, tr_y_xxxy_xyyyz, tr_y_xxxy_xyyzz, tr_y_xxxy_xyzzz, tr_y_xxxy_yyyyy, tr_y_xxxy_yyyyz, tr_y_xxxy_yyyzz, tr_y_xxxy_yyzzz, tr_y_xxxy_yzzzz, tr_y_xxxy_zzzzz, ts_xxxxx_xxxxx, ts_xxxxx_xxxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxzzz, ts_xxxxx_xzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxy_xxxxx[i] = ts_xxxxx_xxxxx[i] * fe_0 + tr_y_xxxxx_xxxxx[i] * pa_y[i];

        tr_y_xxxxxy_xxxxy[i] = 4.0 * tr_y_xxxy_xxxxy[i] * fe_0 + 4.0 * tr_y_xxxxy_xxxy[i] * fe_0 + tr_y_xxxxy_xxxxy[i] * pa_x[i];

        tr_y_xxxxxy_xxxxz[i] = ts_xxxxx_xxxxz[i] * fe_0 + tr_y_xxxxx_xxxxz[i] * pa_y[i];

        tr_y_xxxxxy_xxxyy[i] = 4.0 * tr_y_xxxy_xxxyy[i] * fe_0 + 3.0 * tr_y_xxxxy_xxyy[i] * fe_0 + tr_y_xxxxy_xxxyy[i] * pa_x[i];

        tr_y_xxxxxy_xxxyz[i] = 4.0 * tr_y_xxxy_xxxyz[i] * fe_0 + 3.0 * tr_y_xxxxy_xxyz[i] * fe_0 + tr_y_xxxxy_xxxyz[i] * pa_x[i];

        tr_y_xxxxxy_xxxzz[i] = ts_xxxxx_xxxzz[i] * fe_0 + tr_y_xxxxx_xxxzz[i] * pa_y[i];

        tr_y_xxxxxy_xxyyy[i] = 4.0 * tr_y_xxxy_xxyyy[i] * fe_0 + 2.0 * tr_y_xxxxy_xyyy[i] * fe_0 + tr_y_xxxxy_xxyyy[i] * pa_x[i];

        tr_y_xxxxxy_xxyyz[i] = 4.0 * tr_y_xxxy_xxyyz[i] * fe_0 + 2.0 * tr_y_xxxxy_xyyz[i] * fe_0 + tr_y_xxxxy_xxyyz[i] * pa_x[i];

        tr_y_xxxxxy_xxyzz[i] = 4.0 * tr_y_xxxy_xxyzz[i] * fe_0 + 2.0 * tr_y_xxxxy_xyzz[i] * fe_0 + tr_y_xxxxy_xxyzz[i] * pa_x[i];

        tr_y_xxxxxy_xxzzz[i] = ts_xxxxx_xxzzz[i] * fe_0 + tr_y_xxxxx_xxzzz[i] * pa_y[i];

        tr_y_xxxxxy_xyyyy[i] = 4.0 * tr_y_xxxy_xyyyy[i] * fe_0 + tr_y_xxxxy_yyyy[i] * fe_0 + tr_y_xxxxy_xyyyy[i] * pa_x[i];

        tr_y_xxxxxy_xyyyz[i] = 4.0 * tr_y_xxxy_xyyyz[i] * fe_0 + tr_y_xxxxy_yyyz[i] * fe_0 + tr_y_xxxxy_xyyyz[i] * pa_x[i];

        tr_y_xxxxxy_xyyzz[i] = 4.0 * tr_y_xxxy_xyyzz[i] * fe_0 + tr_y_xxxxy_yyzz[i] * fe_0 + tr_y_xxxxy_xyyzz[i] * pa_x[i];

        tr_y_xxxxxy_xyzzz[i] = 4.0 * tr_y_xxxy_xyzzz[i] * fe_0 + tr_y_xxxxy_yzzz[i] * fe_0 + tr_y_xxxxy_xyzzz[i] * pa_x[i];

        tr_y_xxxxxy_xzzzz[i] = ts_xxxxx_xzzzz[i] * fe_0 + tr_y_xxxxx_xzzzz[i] * pa_y[i];

        tr_y_xxxxxy_yyyyy[i] = 4.0 * tr_y_xxxy_yyyyy[i] * fe_0 + tr_y_xxxxy_yyyyy[i] * pa_x[i];

        tr_y_xxxxxy_yyyyz[i] = 4.0 * tr_y_xxxy_yyyyz[i] * fe_0 + tr_y_xxxxy_yyyyz[i] * pa_x[i];

        tr_y_xxxxxy_yyyzz[i] = 4.0 * tr_y_xxxy_yyyzz[i] * fe_0 + tr_y_xxxxy_yyyzz[i] * pa_x[i];

        tr_y_xxxxxy_yyzzz[i] = 4.0 * tr_y_xxxy_yyzzz[i] * fe_0 + tr_y_xxxxy_yyzzz[i] * pa_x[i];

        tr_y_xxxxxy_yzzzz[i] = 4.0 * tr_y_xxxy_yzzzz[i] * fe_0 + tr_y_xxxxy_yzzzz[i] * pa_x[i];

        tr_y_xxxxxy_zzzzz[i] = 4.0 * tr_y_xxxy_zzzzz[i] * fe_0 + tr_y_xxxxy_zzzzz[i] * pa_x[i];
    }

    // Set up 630-651 components of targeted buffer : IH

    auto tr_y_xxxxxz_xxxxx = pbuffer.data(idx_dip_ih + 630);

    auto tr_y_xxxxxz_xxxxy = pbuffer.data(idx_dip_ih + 631);

    auto tr_y_xxxxxz_xxxxz = pbuffer.data(idx_dip_ih + 632);

    auto tr_y_xxxxxz_xxxyy = pbuffer.data(idx_dip_ih + 633);

    auto tr_y_xxxxxz_xxxyz = pbuffer.data(idx_dip_ih + 634);

    auto tr_y_xxxxxz_xxxzz = pbuffer.data(idx_dip_ih + 635);

    auto tr_y_xxxxxz_xxyyy = pbuffer.data(idx_dip_ih + 636);

    auto tr_y_xxxxxz_xxyyz = pbuffer.data(idx_dip_ih + 637);

    auto tr_y_xxxxxz_xxyzz = pbuffer.data(idx_dip_ih + 638);

    auto tr_y_xxxxxz_xxzzz = pbuffer.data(idx_dip_ih + 639);

    auto tr_y_xxxxxz_xyyyy = pbuffer.data(idx_dip_ih + 640);

    auto tr_y_xxxxxz_xyyyz = pbuffer.data(idx_dip_ih + 641);

    auto tr_y_xxxxxz_xyyzz = pbuffer.data(idx_dip_ih + 642);

    auto tr_y_xxxxxz_xyzzz = pbuffer.data(idx_dip_ih + 643);

    auto tr_y_xxxxxz_xzzzz = pbuffer.data(idx_dip_ih + 644);

    auto tr_y_xxxxxz_yyyyy = pbuffer.data(idx_dip_ih + 645);

    auto tr_y_xxxxxz_yyyyz = pbuffer.data(idx_dip_ih + 646);

    auto tr_y_xxxxxz_yyyzz = pbuffer.data(idx_dip_ih + 647);

    auto tr_y_xxxxxz_yyzzz = pbuffer.data(idx_dip_ih + 648);

    auto tr_y_xxxxxz_yzzzz = pbuffer.data(idx_dip_ih + 649);

    auto tr_y_xxxxxz_zzzzz = pbuffer.data(idx_dip_ih + 650);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxxx_xxxx, tr_y_xxxxx_xxxxx, tr_y_xxxxx_xxxxy, tr_y_xxxxx_xxxxz, tr_y_xxxxx_xxxy, tr_y_xxxxx_xxxyy, tr_y_xxxxx_xxxyz, tr_y_xxxxx_xxxz, tr_y_xxxxx_xxxzz, tr_y_xxxxx_xxyy, tr_y_xxxxx_xxyyy, tr_y_xxxxx_xxyyz, tr_y_xxxxx_xxyz, tr_y_xxxxx_xxyzz, tr_y_xxxxx_xxzz, tr_y_xxxxx_xxzzz, tr_y_xxxxx_xyyy, tr_y_xxxxx_xyyyy, tr_y_xxxxx_xyyyz, tr_y_xxxxx_xyyz, tr_y_xxxxx_xyyzz, tr_y_xxxxx_xyzz, tr_y_xxxxx_xyzzz, tr_y_xxxxx_xzzz, tr_y_xxxxx_xzzzz, tr_y_xxxxx_yyyyy, tr_y_xxxxxz_xxxxx, tr_y_xxxxxz_xxxxy, tr_y_xxxxxz_xxxxz, tr_y_xxxxxz_xxxyy, tr_y_xxxxxz_xxxyz, tr_y_xxxxxz_xxxzz, tr_y_xxxxxz_xxyyy, tr_y_xxxxxz_xxyyz, tr_y_xxxxxz_xxyzz, tr_y_xxxxxz_xxzzz, tr_y_xxxxxz_xyyyy, tr_y_xxxxxz_xyyyz, tr_y_xxxxxz_xyyzz, tr_y_xxxxxz_xyzzz, tr_y_xxxxxz_xzzzz, tr_y_xxxxxz_yyyyy, tr_y_xxxxxz_yyyyz, tr_y_xxxxxz_yyyzz, tr_y_xxxxxz_yyzzz, tr_y_xxxxxz_yzzzz, tr_y_xxxxxz_zzzzz, tr_y_xxxxz_yyyyz, tr_y_xxxxz_yyyzz, tr_y_xxxxz_yyzzz, tr_y_xxxxz_yzzzz, tr_y_xxxxz_zzzzz, tr_y_xxxz_yyyyz, tr_y_xxxz_yyyzz, tr_y_xxxz_yyzzz, tr_y_xxxz_yzzzz, tr_y_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxz_xxxxx[i] = tr_y_xxxxx_xxxxx[i] * pa_z[i];

        tr_y_xxxxxz_xxxxy[i] = tr_y_xxxxx_xxxxy[i] * pa_z[i];

        tr_y_xxxxxz_xxxxz[i] = tr_y_xxxxx_xxxx[i] * fe_0 + tr_y_xxxxx_xxxxz[i] * pa_z[i];

        tr_y_xxxxxz_xxxyy[i] = tr_y_xxxxx_xxxyy[i] * pa_z[i];

        tr_y_xxxxxz_xxxyz[i] = tr_y_xxxxx_xxxy[i] * fe_0 + tr_y_xxxxx_xxxyz[i] * pa_z[i];

        tr_y_xxxxxz_xxxzz[i] = 2.0 * tr_y_xxxxx_xxxz[i] * fe_0 + tr_y_xxxxx_xxxzz[i] * pa_z[i];

        tr_y_xxxxxz_xxyyy[i] = tr_y_xxxxx_xxyyy[i] * pa_z[i];

        tr_y_xxxxxz_xxyyz[i] = tr_y_xxxxx_xxyy[i] * fe_0 + tr_y_xxxxx_xxyyz[i] * pa_z[i];

        tr_y_xxxxxz_xxyzz[i] = 2.0 * tr_y_xxxxx_xxyz[i] * fe_0 + tr_y_xxxxx_xxyzz[i] * pa_z[i];

        tr_y_xxxxxz_xxzzz[i] = 3.0 * tr_y_xxxxx_xxzz[i] * fe_0 + tr_y_xxxxx_xxzzz[i] * pa_z[i];

        tr_y_xxxxxz_xyyyy[i] = tr_y_xxxxx_xyyyy[i] * pa_z[i];

        tr_y_xxxxxz_xyyyz[i] = tr_y_xxxxx_xyyy[i] * fe_0 + tr_y_xxxxx_xyyyz[i] * pa_z[i];

        tr_y_xxxxxz_xyyzz[i] = 2.0 * tr_y_xxxxx_xyyz[i] * fe_0 + tr_y_xxxxx_xyyzz[i] * pa_z[i];

        tr_y_xxxxxz_xyzzz[i] = 3.0 * tr_y_xxxxx_xyzz[i] * fe_0 + tr_y_xxxxx_xyzzz[i] * pa_z[i];

        tr_y_xxxxxz_xzzzz[i] = 4.0 * tr_y_xxxxx_xzzz[i] * fe_0 + tr_y_xxxxx_xzzzz[i] * pa_z[i];

        tr_y_xxxxxz_yyyyy[i] = tr_y_xxxxx_yyyyy[i] * pa_z[i];

        tr_y_xxxxxz_yyyyz[i] = 4.0 * tr_y_xxxz_yyyyz[i] * fe_0 + tr_y_xxxxz_yyyyz[i] * pa_x[i];

        tr_y_xxxxxz_yyyzz[i] = 4.0 * tr_y_xxxz_yyyzz[i] * fe_0 + tr_y_xxxxz_yyyzz[i] * pa_x[i];

        tr_y_xxxxxz_yyzzz[i] = 4.0 * tr_y_xxxz_yyzzz[i] * fe_0 + tr_y_xxxxz_yyzzz[i] * pa_x[i];

        tr_y_xxxxxz_yzzzz[i] = 4.0 * tr_y_xxxz_yzzzz[i] * fe_0 + tr_y_xxxxz_yzzzz[i] * pa_x[i];

        tr_y_xxxxxz_zzzzz[i] = 4.0 * tr_y_xxxz_zzzzz[i] * fe_0 + tr_y_xxxxz_zzzzz[i] * pa_x[i];
    }

    // Set up 651-672 components of targeted buffer : IH

    auto tr_y_xxxxyy_xxxxx = pbuffer.data(idx_dip_ih + 651);

    auto tr_y_xxxxyy_xxxxy = pbuffer.data(idx_dip_ih + 652);

    auto tr_y_xxxxyy_xxxxz = pbuffer.data(idx_dip_ih + 653);

    auto tr_y_xxxxyy_xxxyy = pbuffer.data(idx_dip_ih + 654);

    auto tr_y_xxxxyy_xxxyz = pbuffer.data(idx_dip_ih + 655);

    auto tr_y_xxxxyy_xxxzz = pbuffer.data(idx_dip_ih + 656);

    auto tr_y_xxxxyy_xxyyy = pbuffer.data(idx_dip_ih + 657);

    auto tr_y_xxxxyy_xxyyz = pbuffer.data(idx_dip_ih + 658);

    auto tr_y_xxxxyy_xxyzz = pbuffer.data(idx_dip_ih + 659);

    auto tr_y_xxxxyy_xxzzz = pbuffer.data(idx_dip_ih + 660);

    auto tr_y_xxxxyy_xyyyy = pbuffer.data(idx_dip_ih + 661);

    auto tr_y_xxxxyy_xyyyz = pbuffer.data(idx_dip_ih + 662);

    auto tr_y_xxxxyy_xyyzz = pbuffer.data(idx_dip_ih + 663);

    auto tr_y_xxxxyy_xyzzz = pbuffer.data(idx_dip_ih + 664);

    auto tr_y_xxxxyy_xzzzz = pbuffer.data(idx_dip_ih + 665);

    auto tr_y_xxxxyy_yyyyy = pbuffer.data(idx_dip_ih + 666);

    auto tr_y_xxxxyy_yyyyz = pbuffer.data(idx_dip_ih + 667);

    auto tr_y_xxxxyy_yyyzz = pbuffer.data(idx_dip_ih + 668);

    auto tr_y_xxxxyy_yyzzz = pbuffer.data(idx_dip_ih + 669);

    auto tr_y_xxxxyy_yzzzz = pbuffer.data(idx_dip_ih + 670);

    auto tr_y_xxxxyy_zzzzz = pbuffer.data(idx_dip_ih + 671);

    #pragma omp simd aligned(pa_x, tr_y_xxxxyy_xxxxx, tr_y_xxxxyy_xxxxy, tr_y_xxxxyy_xxxxz, tr_y_xxxxyy_xxxyy, tr_y_xxxxyy_xxxyz, tr_y_xxxxyy_xxxzz, tr_y_xxxxyy_xxyyy, tr_y_xxxxyy_xxyyz, tr_y_xxxxyy_xxyzz, tr_y_xxxxyy_xxzzz, tr_y_xxxxyy_xyyyy, tr_y_xxxxyy_xyyyz, tr_y_xxxxyy_xyyzz, tr_y_xxxxyy_xyzzz, tr_y_xxxxyy_xzzzz, tr_y_xxxxyy_yyyyy, tr_y_xxxxyy_yyyyz, tr_y_xxxxyy_yyyzz, tr_y_xxxxyy_yyzzz, tr_y_xxxxyy_yzzzz, tr_y_xxxxyy_zzzzz, tr_y_xxxyy_xxxx, tr_y_xxxyy_xxxxx, tr_y_xxxyy_xxxxy, tr_y_xxxyy_xxxxz, tr_y_xxxyy_xxxy, tr_y_xxxyy_xxxyy, tr_y_xxxyy_xxxyz, tr_y_xxxyy_xxxz, tr_y_xxxyy_xxxzz, tr_y_xxxyy_xxyy, tr_y_xxxyy_xxyyy, tr_y_xxxyy_xxyyz, tr_y_xxxyy_xxyz, tr_y_xxxyy_xxyzz, tr_y_xxxyy_xxzz, tr_y_xxxyy_xxzzz, tr_y_xxxyy_xyyy, tr_y_xxxyy_xyyyy, tr_y_xxxyy_xyyyz, tr_y_xxxyy_xyyz, tr_y_xxxyy_xyyzz, tr_y_xxxyy_xyzz, tr_y_xxxyy_xyzzz, tr_y_xxxyy_xzzz, tr_y_xxxyy_xzzzz, tr_y_xxxyy_yyyy, tr_y_xxxyy_yyyyy, tr_y_xxxyy_yyyyz, tr_y_xxxyy_yyyz, tr_y_xxxyy_yyyzz, tr_y_xxxyy_yyzz, tr_y_xxxyy_yyzzz, tr_y_xxxyy_yzzz, tr_y_xxxyy_yzzzz, tr_y_xxxyy_zzzz, tr_y_xxxyy_zzzzz, tr_y_xxyy_xxxxx, tr_y_xxyy_xxxxy, tr_y_xxyy_xxxxz, tr_y_xxyy_xxxyy, tr_y_xxyy_xxxyz, tr_y_xxyy_xxxzz, tr_y_xxyy_xxyyy, tr_y_xxyy_xxyyz, tr_y_xxyy_xxyzz, tr_y_xxyy_xxzzz, tr_y_xxyy_xyyyy, tr_y_xxyy_xyyyz, tr_y_xxyy_xyyzz, tr_y_xxyy_xyzzz, tr_y_xxyy_xzzzz, tr_y_xxyy_yyyyy, tr_y_xxyy_yyyyz, tr_y_xxyy_yyyzz, tr_y_xxyy_yyzzz, tr_y_xxyy_yzzzz, tr_y_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyy_xxxxx[i] = 3.0 * tr_y_xxyy_xxxxx[i] * fe_0 + 5.0 * tr_y_xxxyy_xxxx[i] * fe_0 + tr_y_xxxyy_xxxxx[i] * pa_x[i];

        tr_y_xxxxyy_xxxxy[i] = 3.0 * tr_y_xxyy_xxxxy[i] * fe_0 + 4.0 * tr_y_xxxyy_xxxy[i] * fe_0 + tr_y_xxxyy_xxxxy[i] * pa_x[i];

        tr_y_xxxxyy_xxxxz[i] = 3.0 * tr_y_xxyy_xxxxz[i] * fe_0 + 4.0 * tr_y_xxxyy_xxxz[i] * fe_0 + tr_y_xxxyy_xxxxz[i] * pa_x[i];

        tr_y_xxxxyy_xxxyy[i] = 3.0 * tr_y_xxyy_xxxyy[i] * fe_0 + 3.0 * tr_y_xxxyy_xxyy[i] * fe_0 + tr_y_xxxyy_xxxyy[i] * pa_x[i];

        tr_y_xxxxyy_xxxyz[i] = 3.0 * tr_y_xxyy_xxxyz[i] * fe_0 + 3.0 * tr_y_xxxyy_xxyz[i] * fe_0 + tr_y_xxxyy_xxxyz[i] * pa_x[i];

        tr_y_xxxxyy_xxxzz[i] = 3.0 * tr_y_xxyy_xxxzz[i] * fe_0 + 3.0 * tr_y_xxxyy_xxzz[i] * fe_0 + tr_y_xxxyy_xxxzz[i] * pa_x[i];

        tr_y_xxxxyy_xxyyy[i] = 3.0 * tr_y_xxyy_xxyyy[i] * fe_0 + 2.0 * tr_y_xxxyy_xyyy[i] * fe_0 + tr_y_xxxyy_xxyyy[i] * pa_x[i];

        tr_y_xxxxyy_xxyyz[i] = 3.0 * tr_y_xxyy_xxyyz[i] * fe_0 + 2.0 * tr_y_xxxyy_xyyz[i] * fe_0 + tr_y_xxxyy_xxyyz[i] * pa_x[i];

        tr_y_xxxxyy_xxyzz[i] = 3.0 * tr_y_xxyy_xxyzz[i] * fe_0 + 2.0 * tr_y_xxxyy_xyzz[i] * fe_0 + tr_y_xxxyy_xxyzz[i] * pa_x[i];

        tr_y_xxxxyy_xxzzz[i] = 3.0 * tr_y_xxyy_xxzzz[i] * fe_0 + 2.0 * tr_y_xxxyy_xzzz[i] * fe_0 + tr_y_xxxyy_xxzzz[i] * pa_x[i];

        tr_y_xxxxyy_xyyyy[i] = 3.0 * tr_y_xxyy_xyyyy[i] * fe_0 + tr_y_xxxyy_yyyy[i] * fe_0 + tr_y_xxxyy_xyyyy[i] * pa_x[i];

        tr_y_xxxxyy_xyyyz[i] = 3.0 * tr_y_xxyy_xyyyz[i] * fe_0 + tr_y_xxxyy_yyyz[i] * fe_0 + tr_y_xxxyy_xyyyz[i] * pa_x[i];

        tr_y_xxxxyy_xyyzz[i] = 3.0 * tr_y_xxyy_xyyzz[i] * fe_0 + tr_y_xxxyy_yyzz[i] * fe_0 + tr_y_xxxyy_xyyzz[i] * pa_x[i];

        tr_y_xxxxyy_xyzzz[i] = 3.0 * tr_y_xxyy_xyzzz[i] * fe_0 + tr_y_xxxyy_yzzz[i] * fe_0 + tr_y_xxxyy_xyzzz[i] * pa_x[i];

        tr_y_xxxxyy_xzzzz[i] = 3.0 * tr_y_xxyy_xzzzz[i] * fe_0 + tr_y_xxxyy_zzzz[i] * fe_0 + tr_y_xxxyy_xzzzz[i] * pa_x[i];

        tr_y_xxxxyy_yyyyy[i] = 3.0 * tr_y_xxyy_yyyyy[i] * fe_0 + tr_y_xxxyy_yyyyy[i] * pa_x[i];

        tr_y_xxxxyy_yyyyz[i] = 3.0 * tr_y_xxyy_yyyyz[i] * fe_0 + tr_y_xxxyy_yyyyz[i] * pa_x[i];

        tr_y_xxxxyy_yyyzz[i] = 3.0 * tr_y_xxyy_yyyzz[i] * fe_0 + tr_y_xxxyy_yyyzz[i] * pa_x[i];

        tr_y_xxxxyy_yyzzz[i] = 3.0 * tr_y_xxyy_yyzzz[i] * fe_0 + tr_y_xxxyy_yyzzz[i] * pa_x[i];

        tr_y_xxxxyy_yzzzz[i] = 3.0 * tr_y_xxyy_yzzzz[i] * fe_0 + tr_y_xxxyy_yzzzz[i] * pa_x[i];

        tr_y_xxxxyy_zzzzz[i] = 3.0 * tr_y_xxyy_zzzzz[i] * fe_0 + tr_y_xxxyy_zzzzz[i] * pa_x[i];
    }

    // Set up 672-693 components of targeted buffer : IH

    auto tr_y_xxxxyz_xxxxx = pbuffer.data(idx_dip_ih + 672);

    auto tr_y_xxxxyz_xxxxy = pbuffer.data(idx_dip_ih + 673);

    auto tr_y_xxxxyz_xxxxz = pbuffer.data(idx_dip_ih + 674);

    auto tr_y_xxxxyz_xxxyy = pbuffer.data(idx_dip_ih + 675);

    auto tr_y_xxxxyz_xxxyz = pbuffer.data(idx_dip_ih + 676);

    auto tr_y_xxxxyz_xxxzz = pbuffer.data(idx_dip_ih + 677);

    auto tr_y_xxxxyz_xxyyy = pbuffer.data(idx_dip_ih + 678);

    auto tr_y_xxxxyz_xxyyz = pbuffer.data(idx_dip_ih + 679);

    auto tr_y_xxxxyz_xxyzz = pbuffer.data(idx_dip_ih + 680);

    auto tr_y_xxxxyz_xxzzz = pbuffer.data(idx_dip_ih + 681);

    auto tr_y_xxxxyz_xyyyy = pbuffer.data(idx_dip_ih + 682);

    auto tr_y_xxxxyz_xyyyz = pbuffer.data(idx_dip_ih + 683);

    auto tr_y_xxxxyz_xyyzz = pbuffer.data(idx_dip_ih + 684);

    auto tr_y_xxxxyz_xyzzz = pbuffer.data(idx_dip_ih + 685);

    auto tr_y_xxxxyz_xzzzz = pbuffer.data(idx_dip_ih + 686);

    auto tr_y_xxxxyz_yyyyy = pbuffer.data(idx_dip_ih + 687);

    auto tr_y_xxxxyz_yyyyz = pbuffer.data(idx_dip_ih + 688);

    auto tr_y_xxxxyz_yyyzz = pbuffer.data(idx_dip_ih + 689);

    auto tr_y_xxxxyz_yyzzz = pbuffer.data(idx_dip_ih + 690);

    auto tr_y_xxxxyz_yzzzz = pbuffer.data(idx_dip_ih + 691);

    auto tr_y_xxxxyz_zzzzz = pbuffer.data(idx_dip_ih + 692);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxxy_xxxxx, tr_y_xxxxy_xxxxy, tr_y_xxxxy_xxxy, tr_y_xxxxy_xxxyy, tr_y_xxxxy_xxxyz, tr_y_xxxxy_xxyy, tr_y_xxxxy_xxyyy, tr_y_xxxxy_xxyyz, tr_y_xxxxy_xxyz, tr_y_xxxxy_xxyzz, tr_y_xxxxy_xyyy, tr_y_xxxxy_xyyyy, tr_y_xxxxy_xyyyz, tr_y_xxxxy_xyyz, tr_y_xxxxy_xyyzz, tr_y_xxxxy_xyzz, tr_y_xxxxy_xyzzz, tr_y_xxxxy_yyyyy, tr_y_xxxxyz_xxxxx, tr_y_xxxxyz_xxxxy, tr_y_xxxxyz_xxxxz, tr_y_xxxxyz_xxxyy, tr_y_xxxxyz_xxxyz, tr_y_xxxxyz_xxxzz, tr_y_xxxxyz_xxyyy, tr_y_xxxxyz_xxyyz, tr_y_xxxxyz_xxyzz, tr_y_xxxxyz_xxzzz, tr_y_xxxxyz_xyyyy, tr_y_xxxxyz_xyyyz, tr_y_xxxxyz_xyyzz, tr_y_xxxxyz_xyzzz, tr_y_xxxxyz_xzzzz, tr_y_xxxxyz_yyyyy, tr_y_xxxxyz_yyyyz, tr_y_xxxxyz_yyyzz, tr_y_xxxxyz_yyzzz, tr_y_xxxxyz_yzzzz, tr_y_xxxxyz_zzzzz, tr_y_xxxxz_xxxxz, tr_y_xxxxz_xxxzz, tr_y_xxxxz_xxzzz, tr_y_xxxxz_xzzzz, tr_y_xxxyz_yyyyz, tr_y_xxxyz_yyyzz, tr_y_xxxyz_yyzzz, tr_y_xxxyz_yzzzz, tr_y_xxxyz_zzzzz, tr_y_xxyz_yyyyz, tr_y_xxyz_yyyzz, tr_y_xxyz_yyzzz, tr_y_xxyz_yzzzz, tr_y_xxyz_zzzzz, ts_xxxxz_xxxxz, ts_xxxxz_xxxzz, ts_xxxxz_xxzzz, ts_xxxxz_xzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyz_xxxxx[i] = tr_y_xxxxy_xxxxx[i] * pa_z[i];

        tr_y_xxxxyz_xxxxy[i] = tr_y_xxxxy_xxxxy[i] * pa_z[i];

        tr_y_xxxxyz_xxxxz[i] = ts_xxxxz_xxxxz[i] * fe_0 + tr_y_xxxxz_xxxxz[i] * pa_y[i];

        tr_y_xxxxyz_xxxyy[i] = tr_y_xxxxy_xxxyy[i] * pa_z[i];

        tr_y_xxxxyz_xxxyz[i] = tr_y_xxxxy_xxxy[i] * fe_0 + tr_y_xxxxy_xxxyz[i] * pa_z[i];

        tr_y_xxxxyz_xxxzz[i] = ts_xxxxz_xxxzz[i] * fe_0 + tr_y_xxxxz_xxxzz[i] * pa_y[i];

        tr_y_xxxxyz_xxyyy[i] = tr_y_xxxxy_xxyyy[i] * pa_z[i];

        tr_y_xxxxyz_xxyyz[i] = tr_y_xxxxy_xxyy[i] * fe_0 + tr_y_xxxxy_xxyyz[i] * pa_z[i];

        tr_y_xxxxyz_xxyzz[i] = 2.0 * tr_y_xxxxy_xxyz[i] * fe_0 + tr_y_xxxxy_xxyzz[i] * pa_z[i];

        tr_y_xxxxyz_xxzzz[i] = ts_xxxxz_xxzzz[i] * fe_0 + tr_y_xxxxz_xxzzz[i] * pa_y[i];

        tr_y_xxxxyz_xyyyy[i] = tr_y_xxxxy_xyyyy[i] * pa_z[i];

        tr_y_xxxxyz_xyyyz[i] = tr_y_xxxxy_xyyy[i] * fe_0 + tr_y_xxxxy_xyyyz[i] * pa_z[i];

        tr_y_xxxxyz_xyyzz[i] = 2.0 * tr_y_xxxxy_xyyz[i] * fe_0 + tr_y_xxxxy_xyyzz[i] * pa_z[i];

        tr_y_xxxxyz_xyzzz[i] = 3.0 * tr_y_xxxxy_xyzz[i] * fe_0 + tr_y_xxxxy_xyzzz[i] * pa_z[i];

        tr_y_xxxxyz_xzzzz[i] = ts_xxxxz_xzzzz[i] * fe_0 + tr_y_xxxxz_xzzzz[i] * pa_y[i];

        tr_y_xxxxyz_yyyyy[i] = tr_y_xxxxy_yyyyy[i] * pa_z[i];

        tr_y_xxxxyz_yyyyz[i] = 3.0 * tr_y_xxyz_yyyyz[i] * fe_0 + tr_y_xxxyz_yyyyz[i] * pa_x[i];

        tr_y_xxxxyz_yyyzz[i] = 3.0 * tr_y_xxyz_yyyzz[i] * fe_0 + tr_y_xxxyz_yyyzz[i] * pa_x[i];

        tr_y_xxxxyz_yyzzz[i] = 3.0 * tr_y_xxyz_yyzzz[i] * fe_0 + tr_y_xxxyz_yyzzz[i] * pa_x[i];

        tr_y_xxxxyz_yzzzz[i] = 3.0 * tr_y_xxyz_yzzzz[i] * fe_0 + tr_y_xxxyz_yzzzz[i] * pa_x[i];

        tr_y_xxxxyz_zzzzz[i] = 3.0 * tr_y_xxyz_zzzzz[i] * fe_0 + tr_y_xxxyz_zzzzz[i] * pa_x[i];
    }

    // Set up 693-714 components of targeted buffer : IH

    auto tr_y_xxxxzz_xxxxx = pbuffer.data(idx_dip_ih + 693);

    auto tr_y_xxxxzz_xxxxy = pbuffer.data(idx_dip_ih + 694);

    auto tr_y_xxxxzz_xxxxz = pbuffer.data(idx_dip_ih + 695);

    auto tr_y_xxxxzz_xxxyy = pbuffer.data(idx_dip_ih + 696);

    auto tr_y_xxxxzz_xxxyz = pbuffer.data(idx_dip_ih + 697);

    auto tr_y_xxxxzz_xxxzz = pbuffer.data(idx_dip_ih + 698);

    auto tr_y_xxxxzz_xxyyy = pbuffer.data(idx_dip_ih + 699);

    auto tr_y_xxxxzz_xxyyz = pbuffer.data(idx_dip_ih + 700);

    auto tr_y_xxxxzz_xxyzz = pbuffer.data(idx_dip_ih + 701);

    auto tr_y_xxxxzz_xxzzz = pbuffer.data(idx_dip_ih + 702);

    auto tr_y_xxxxzz_xyyyy = pbuffer.data(idx_dip_ih + 703);

    auto tr_y_xxxxzz_xyyyz = pbuffer.data(idx_dip_ih + 704);

    auto tr_y_xxxxzz_xyyzz = pbuffer.data(idx_dip_ih + 705);

    auto tr_y_xxxxzz_xyzzz = pbuffer.data(idx_dip_ih + 706);

    auto tr_y_xxxxzz_xzzzz = pbuffer.data(idx_dip_ih + 707);

    auto tr_y_xxxxzz_yyyyy = pbuffer.data(idx_dip_ih + 708);

    auto tr_y_xxxxzz_yyyyz = pbuffer.data(idx_dip_ih + 709);

    auto tr_y_xxxxzz_yyyzz = pbuffer.data(idx_dip_ih + 710);

    auto tr_y_xxxxzz_yyzzz = pbuffer.data(idx_dip_ih + 711);

    auto tr_y_xxxxzz_yzzzz = pbuffer.data(idx_dip_ih + 712);

    auto tr_y_xxxxzz_zzzzz = pbuffer.data(idx_dip_ih + 713);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxx_xxxxx, tr_y_xxxx_xxxxy, tr_y_xxxx_xxxyy, tr_y_xxxx_xxyyy, tr_y_xxxx_xyyyy, tr_y_xxxxz_xxxxx, tr_y_xxxxz_xxxxy, tr_y_xxxxz_xxxyy, tr_y_xxxxz_xxyyy, tr_y_xxxxz_xyyyy, tr_y_xxxxzz_xxxxx, tr_y_xxxxzz_xxxxy, tr_y_xxxxzz_xxxxz, tr_y_xxxxzz_xxxyy, tr_y_xxxxzz_xxxyz, tr_y_xxxxzz_xxxzz, tr_y_xxxxzz_xxyyy, tr_y_xxxxzz_xxyyz, tr_y_xxxxzz_xxyzz, tr_y_xxxxzz_xxzzz, tr_y_xxxxzz_xyyyy, tr_y_xxxxzz_xyyyz, tr_y_xxxxzz_xyyzz, tr_y_xxxxzz_xyzzz, tr_y_xxxxzz_xzzzz, tr_y_xxxxzz_yyyyy, tr_y_xxxxzz_yyyyz, tr_y_xxxxzz_yyyzz, tr_y_xxxxzz_yyzzz, tr_y_xxxxzz_yzzzz, tr_y_xxxxzz_zzzzz, tr_y_xxxzz_xxxxz, tr_y_xxxzz_xxxyz, tr_y_xxxzz_xxxz, tr_y_xxxzz_xxxzz, tr_y_xxxzz_xxyyz, tr_y_xxxzz_xxyz, tr_y_xxxzz_xxyzz, tr_y_xxxzz_xxzz, tr_y_xxxzz_xxzzz, tr_y_xxxzz_xyyyz, tr_y_xxxzz_xyyz, tr_y_xxxzz_xyyzz, tr_y_xxxzz_xyzz, tr_y_xxxzz_xyzzz, tr_y_xxxzz_xzzz, tr_y_xxxzz_xzzzz, tr_y_xxxzz_yyyyy, tr_y_xxxzz_yyyyz, tr_y_xxxzz_yyyz, tr_y_xxxzz_yyyzz, tr_y_xxxzz_yyzz, tr_y_xxxzz_yyzzz, tr_y_xxxzz_yzzz, tr_y_xxxzz_yzzzz, tr_y_xxxzz_zzzz, tr_y_xxxzz_zzzzz, tr_y_xxzz_xxxxz, tr_y_xxzz_xxxyz, tr_y_xxzz_xxxzz, tr_y_xxzz_xxyyz, tr_y_xxzz_xxyzz, tr_y_xxzz_xxzzz, tr_y_xxzz_xyyyz, tr_y_xxzz_xyyzz, tr_y_xxzz_xyzzz, tr_y_xxzz_xzzzz, tr_y_xxzz_yyyyy, tr_y_xxzz_yyyyz, tr_y_xxzz_yyyzz, tr_y_xxzz_yyzzz, tr_y_xxzz_yzzzz, tr_y_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxzz_xxxxx[i] = tr_y_xxxx_xxxxx[i] * fe_0 + tr_y_xxxxz_xxxxx[i] * pa_z[i];

        tr_y_xxxxzz_xxxxy[i] = tr_y_xxxx_xxxxy[i] * fe_0 + tr_y_xxxxz_xxxxy[i] * pa_z[i];

        tr_y_xxxxzz_xxxxz[i] = 3.0 * tr_y_xxzz_xxxxz[i] * fe_0 + 4.0 * tr_y_xxxzz_xxxz[i] * fe_0 + tr_y_xxxzz_xxxxz[i] * pa_x[i];

        tr_y_xxxxzz_xxxyy[i] = tr_y_xxxx_xxxyy[i] * fe_0 + tr_y_xxxxz_xxxyy[i] * pa_z[i];

        tr_y_xxxxzz_xxxyz[i] = 3.0 * tr_y_xxzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xxxzz_xxyz[i] * fe_0 + tr_y_xxxzz_xxxyz[i] * pa_x[i];

        tr_y_xxxxzz_xxxzz[i] = 3.0 * tr_y_xxzz_xxxzz[i] * fe_0 + 3.0 * tr_y_xxxzz_xxzz[i] * fe_0 + tr_y_xxxzz_xxxzz[i] * pa_x[i];

        tr_y_xxxxzz_xxyyy[i] = tr_y_xxxx_xxyyy[i] * fe_0 + tr_y_xxxxz_xxyyy[i] * pa_z[i];

        tr_y_xxxxzz_xxyyz[i] = 3.0 * tr_y_xxzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xxxzz_xyyz[i] * fe_0 + tr_y_xxxzz_xxyyz[i] * pa_x[i];

        tr_y_xxxxzz_xxyzz[i] = 3.0 * tr_y_xxzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xxxzz_xyzz[i] * fe_0 + tr_y_xxxzz_xxyzz[i] * pa_x[i];

        tr_y_xxxxzz_xxzzz[i] = 3.0 * tr_y_xxzz_xxzzz[i] * fe_0 + 2.0 * tr_y_xxxzz_xzzz[i] * fe_0 + tr_y_xxxzz_xxzzz[i] * pa_x[i];

        tr_y_xxxxzz_xyyyy[i] = tr_y_xxxx_xyyyy[i] * fe_0 + tr_y_xxxxz_xyyyy[i] * pa_z[i];

        tr_y_xxxxzz_xyyyz[i] = 3.0 * tr_y_xxzz_xyyyz[i] * fe_0 + tr_y_xxxzz_yyyz[i] * fe_0 + tr_y_xxxzz_xyyyz[i] * pa_x[i];

        tr_y_xxxxzz_xyyzz[i] = 3.0 * tr_y_xxzz_xyyzz[i] * fe_0 + tr_y_xxxzz_yyzz[i] * fe_0 + tr_y_xxxzz_xyyzz[i] * pa_x[i];

        tr_y_xxxxzz_xyzzz[i] = 3.0 * tr_y_xxzz_xyzzz[i] * fe_0 + tr_y_xxxzz_yzzz[i] * fe_0 + tr_y_xxxzz_xyzzz[i] * pa_x[i];

        tr_y_xxxxzz_xzzzz[i] = 3.0 * tr_y_xxzz_xzzzz[i] * fe_0 + tr_y_xxxzz_zzzz[i] * fe_0 + tr_y_xxxzz_xzzzz[i] * pa_x[i];

        tr_y_xxxxzz_yyyyy[i] = 3.0 * tr_y_xxzz_yyyyy[i] * fe_0 + tr_y_xxxzz_yyyyy[i] * pa_x[i];

        tr_y_xxxxzz_yyyyz[i] = 3.0 * tr_y_xxzz_yyyyz[i] * fe_0 + tr_y_xxxzz_yyyyz[i] * pa_x[i];

        tr_y_xxxxzz_yyyzz[i] = 3.0 * tr_y_xxzz_yyyzz[i] * fe_0 + tr_y_xxxzz_yyyzz[i] * pa_x[i];

        tr_y_xxxxzz_yyzzz[i] = 3.0 * tr_y_xxzz_yyzzz[i] * fe_0 + tr_y_xxxzz_yyzzz[i] * pa_x[i];

        tr_y_xxxxzz_yzzzz[i] = 3.0 * tr_y_xxzz_yzzzz[i] * fe_0 + tr_y_xxxzz_yzzzz[i] * pa_x[i];

        tr_y_xxxxzz_zzzzz[i] = 3.0 * tr_y_xxzz_zzzzz[i] * fe_0 + tr_y_xxxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 714-735 components of targeted buffer : IH

    auto tr_y_xxxyyy_xxxxx = pbuffer.data(idx_dip_ih + 714);

    auto tr_y_xxxyyy_xxxxy = pbuffer.data(idx_dip_ih + 715);

    auto tr_y_xxxyyy_xxxxz = pbuffer.data(idx_dip_ih + 716);

    auto tr_y_xxxyyy_xxxyy = pbuffer.data(idx_dip_ih + 717);

    auto tr_y_xxxyyy_xxxyz = pbuffer.data(idx_dip_ih + 718);

    auto tr_y_xxxyyy_xxxzz = pbuffer.data(idx_dip_ih + 719);

    auto tr_y_xxxyyy_xxyyy = pbuffer.data(idx_dip_ih + 720);

    auto tr_y_xxxyyy_xxyyz = pbuffer.data(idx_dip_ih + 721);

    auto tr_y_xxxyyy_xxyzz = pbuffer.data(idx_dip_ih + 722);

    auto tr_y_xxxyyy_xxzzz = pbuffer.data(idx_dip_ih + 723);

    auto tr_y_xxxyyy_xyyyy = pbuffer.data(idx_dip_ih + 724);

    auto tr_y_xxxyyy_xyyyz = pbuffer.data(idx_dip_ih + 725);

    auto tr_y_xxxyyy_xyyzz = pbuffer.data(idx_dip_ih + 726);

    auto tr_y_xxxyyy_xyzzz = pbuffer.data(idx_dip_ih + 727);

    auto tr_y_xxxyyy_xzzzz = pbuffer.data(idx_dip_ih + 728);

    auto tr_y_xxxyyy_yyyyy = pbuffer.data(idx_dip_ih + 729);

    auto tr_y_xxxyyy_yyyyz = pbuffer.data(idx_dip_ih + 730);

    auto tr_y_xxxyyy_yyyzz = pbuffer.data(idx_dip_ih + 731);

    auto tr_y_xxxyyy_yyzzz = pbuffer.data(idx_dip_ih + 732);

    auto tr_y_xxxyyy_yzzzz = pbuffer.data(idx_dip_ih + 733);

    auto tr_y_xxxyyy_zzzzz = pbuffer.data(idx_dip_ih + 734);

    #pragma omp simd aligned(pa_x, tr_y_xxxyyy_xxxxx, tr_y_xxxyyy_xxxxy, tr_y_xxxyyy_xxxxz, tr_y_xxxyyy_xxxyy, tr_y_xxxyyy_xxxyz, tr_y_xxxyyy_xxxzz, tr_y_xxxyyy_xxyyy, tr_y_xxxyyy_xxyyz, tr_y_xxxyyy_xxyzz, tr_y_xxxyyy_xxzzz, tr_y_xxxyyy_xyyyy, tr_y_xxxyyy_xyyyz, tr_y_xxxyyy_xyyzz, tr_y_xxxyyy_xyzzz, tr_y_xxxyyy_xzzzz, tr_y_xxxyyy_yyyyy, tr_y_xxxyyy_yyyyz, tr_y_xxxyyy_yyyzz, tr_y_xxxyyy_yyzzz, tr_y_xxxyyy_yzzzz, tr_y_xxxyyy_zzzzz, tr_y_xxyyy_xxxx, tr_y_xxyyy_xxxxx, tr_y_xxyyy_xxxxy, tr_y_xxyyy_xxxxz, tr_y_xxyyy_xxxy, tr_y_xxyyy_xxxyy, tr_y_xxyyy_xxxyz, tr_y_xxyyy_xxxz, tr_y_xxyyy_xxxzz, tr_y_xxyyy_xxyy, tr_y_xxyyy_xxyyy, tr_y_xxyyy_xxyyz, tr_y_xxyyy_xxyz, tr_y_xxyyy_xxyzz, tr_y_xxyyy_xxzz, tr_y_xxyyy_xxzzz, tr_y_xxyyy_xyyy, tr_y_xxyyy_xyyyy, tr_y_xxyyy_xyyyz, tr_y_xxyyy_xyyz, tr_y_xxyyy_xyyzz, tr_y_xxyyy_xyzz, tr_y_xxyyy_xyzzz, tr_y_xxyyy_xzzz, tr_y_xxyyy_xzzzz, tr_y_xxyyy_yyyy, tr_y_xxyyy_yyyyy, tr_y_xxyyy_yyyyz, tr_y_xxyyy_yyyz, tr_y_xxyyy_yyyzz, tr_y_xxyyy_yyzz, tr_y_xxyyy_yyzzz, tr_y_xxyyy_yzzz, tr_y_xxyyy_yzzzz, tr_y_xxyyy_zzzz, tr_y_xxyyy_zzzzz, tr_y_xyyy_xxxxx, tr_y_xyyy_xxxxy, tr_y_xyyy_xxxxz, tr_y_xyyy_xxxyy, tr_y_xyyy_xxxyz, tr_y_xyyy_xxxzz, tr_y_xyyy_xxyyy, tr_y_xyyy_xxyyz, tr_y_xyyy_xxyzz, tr_y_xyyy_xxzzz, tr_y_xyyy_xyyyy, tr_y_xyyy_xyyyz, tr_y_xyyy_xyyzz, tr_y_xyyy_xyzzz, tr_y_xyyy_xzzzz, tr_y_xyyy_yyyyy, tr_y_xyyy_yyyyz, tr_y_xyyy_yyyzz, tr_y_xyyy_yyzzz, tr_y_xyyy_yzzzz, tr_y_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyy_xxxxx[i] = 2.0 * tr_y_xyyy_xxxxx[i] * fe_0 + 5.0 * tr_y_xxyyy_xxxx[i] * fe_0 + tr_y_xxyyy_xxxxx[i] * pa_x[i];

        tr_y_xxxyyy_xxxxy[i] = 2.0 * tr_y_xyyy_xxxxy[i] * fe_0 + 4.0 * tr_y_xxyyy_xxxy[i] * fe_0 + tr_y_xxyyy_xxxxy[i] * pa_x[i];

        tr_y_xxxyyy_xxxxz[i] = 2.0 * tr_y_xyyy_xxxxz[i] * fe_0 + 4.0 * tr_y_xxyyy_xxxz[i] * fe_0 + tr_y_xxyyy_xxxxz[i] * pa_x[i];

        tr_y_xxxyyy_xxxyy[i] = 2.0 * tr_y_xyyy_xxxyy[i] * fe_0 + 3.0 * tr_y_xxyyy_xxyy[i] * fe_0 + tr_y_xxyyy_xxxyy[i] * pa_x[i];

        tr_y_xxxyyy_xxxyz[i] = 2.0 * tr_y_xyyy_xxxyz[i] * fe_0 + 3.0 * tr_y_xxyyy_xxyz[i] * fe_0 + tr_y_xxyyy_xxxyz[i] * pa_x[i];

        tr_y_xxxyyy_xxxzz[i] = 2.0 * tr_y_xyyy_xxxzz[i] * fe_0 + 3.0 * tr_y_xxyyy_xxzz[i] * fe_0 + tr_y_xxyyy_xxxzz[i] * pa_x[i];

        tr_y_xxxyyy_xxyyy[i] = 2.0 * tr_y_xyyy_xxyyy[i] * fe_0 + 2.0 * tr_y_xxyyy_xyyy[i] * fe_0 + tr_y_xxyyy_xxyyy[i] * pa_x[i];

        tr_y_xxxyyy_xxyyz[i] = 2.0 * tr_y_xyyy_xxyyz[i] * fe_0 + 2.0 * tr_y_xxyyy_xyyz[i] * fe_0 + tr_y_xxyyy_xxyyz[i] * pa_x[i];

        tr_y_xxxyyy_xxyzz[i] = 2.0 * tr_y_xyyy_xxyzz[i] * fe_0 + 2.0 * tr_y_xxyyy_xyzz[i] * fe_0 + tr_y_xxyyy_xxyzz[i] * pa_x[i];

        tr_y_xxxyyy_xxzzz[i] = 2.0 * tr_y_xyyy_xxzzz[i] * fe_0 + 2.0 * tr_y_xxyyy_xzzz[i] * fe_0 + tr_y_xxyyy_xxzzz[i] * pa_x[i];

        tr_y_xxxyyy_xyyyy[i] = 2.0 * tr_y_xyyy_xyyyy[i] * fe_0 + tr_y_xxyyy_yyyy[i] * fe_0 + tr_y_xxyyy_xyyyy[i] * pa_x[i];

        tr_y_xxxyyy_xyyyz[i] = 2.0 * tr_y_xyyy_xyyyz[i] * fe_0 + tr_y_xxyyy_yyyz[i] * fe_0 + tr_y_xxyyy_xyyyz[i] * pa_x[i];

        tr_y_xxxyyy_xyyzz[i] = 2.0 * tr_y_xyyy_xyyzz[i] * fe_0 + tr_y_xxyyy_yyzz[i] * fe_0 + tr_y_xxyyy_xyyzz[i] * pa_x[i];

        tr_y_xxxyyy_xyzzz[i] = 2.0 * tr_y_xyyy_xyzzz[i] * fe_0 + tr_y_xxyyy_yzzz[i] * fe_0 + tr_y_xxyyy_xyzzz[i] * pa_x[i];

        tr_y_xxxyyy_xzzzz[i] = 2.0 * tr_y_xyyy_xzzzz[i] * fe_0 + tr_y_xxyyy_zzzz[i] * fe_0 + tr_y_xxyyy_xzzzz[i] * pa_x[i];

        tr_y_xxxyyy_yyyyy[i] = 2.0 * tr_y_xyyy_yyyyy[i] * fe_0 + tr_y_xxyyy_yyyyy[i] * pa_x[i];

        tr_y_xxxyyy_yyyyz[i] = 2.0 * tr_y_xyyy_yyyyz[i] * fe_0 + tr_y_xxyyy_yyyyz[i] * pa_x[i];

        tr_y_xxxyyy_yyyzz[i] = 2.0 * tr_y_xyyy_yyyzz[i] * fe_0 + tr_y_xxyyy_yyyzz[i] * pa_x[i];

        tr_y_xxxyyy_yyzzz[i] = 2.0 * tr_y_xyyy_yyzzz[i] * fe_0 + tr_y_xxyyy_yyzzz[i] * pa_x[i];

        tr_y_xxxyyy_yzzzz[i] = 2.0 * tr_y_xyyy_yzzzz[i] * fe_0 + tr_y_xxyyy_yzzzz[i] * pa_x[i];

        tr_y_xxxyyy_zzzzz[i] = 2.0 * tr_y_xyyy_zzzzz[i] * fe_0 + tr_y_xxyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 735-756 components of targeted buffer : IH

    auto tr_y_xxxyyz_xxxxx = pbuffer.data(idx_dip_ih + 735);

    auto tr_y_xxxyyz_xxxxy = pbuffer.data(idx_dip_ih + 736);

    auto tr_y_xxxyyz_xxxxz = pbuffer.data(idx_dip_ih + 737);

    auto tr_y_xxxyyz_xxxyy = pbuffer.data(idx_dip_ih + 738);

    auto tr_y_xxxyyz_xxxyz = pbuffer.data(idx_dip_ih + 739);

    auto tr_y_xxxyyz_xxxzz = pbuffer.data(idx_dip_ih + 740);

    auto tr_y_xxxyyz_xxyyy = pbuffer.data(idx_dip_ih + 741);

    auto tr_y_xxxyyz_xxyyz = pbuffer.data(idx_dip_ih + 742);

    auto tr_y_xxxyyz_xxyzz = pbuffer.data(idx_dip_ih + 743);

    auto tr_y_xxxyyz_xxzzz = pbuffer.data(idx_dip_ih + 744);

    auto tr_y_xxxyyz_xyyyy = pbuffer.data(idx_dip_ih + 745);

    auto tr_y_xxxyyz_xyyyz = pbuffer.data(idx_dip_ih + 746);

    auto tr_y_xxxyyz_xyyzz = pbuffer.data(idx_dip_ih + 747);

    auto tr_y_xxxyyz_xyzzz = pbuffer.data(idx_dip_ih + 748);

    auto tr_y_xxxyyz_xzzzz = pbuffer.data(idx_dip_ih + 749);

    auto tr_y_xxxyyz_yyyyy = pbuffer.data(idx_dip_ih + 750);

    auto tr_y_xxxyyz_yyyyz = pbuffer.data(idx_dip_ih + 751);

    auto tr_y_xxxyyz_yyyzz = pbuffer.data(idx_dip_ih + 752);

    auto tr_y_xxxyyz_yyzzz = pbuffer.data(idx_dip_ih + 753);

    auto tr_y_xxxyyz_yzzzz = pbuffer.data(idx_dip_ih + 754);

    auto tr_y_xxxyyz_zzzzz = pbuffer.data(idx_dip_ih + 755);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxyy_xxxx, tr_y_xxxyy_xxxxx, tr_y_xxxyy_xxxxy, tr_y_xxxyy_xxxxz, tr_y_xxxyy_xxxy, tr_y_xxxyy_xxxyy, tr_y_xxxyy_xxxyz, tr_y_xxxyy_xxxz, tr_y_xxxyy_xxxzz, tr_y_xxxyy_xxyy, tr_y_xxxyy_xxyyy, tr_y_xxxyy_xxyyz, tr_y_xxxyy_xxyz, tr_y_xxxyy_xxyzz, tr_y_xxxyy_xxzz, tr_y_xxxyy_xxzzz, tr_y_xxxyy_xyyy, tr_y_xxxyy_xyyyy, tr_y_xxxyy_xyyyz, tr_y_xxxyy_xyyz, tr_y_xxxyy_xyyzz, tr_y_xxxyy_xyzz, tr_y_xxxyy_xyzzz, tr_y_xxxyy_xzzz, tr_y_xxxyy_xzzzz, tr_y_xxxyy_yyyyy, tr_y_xxxyyz_xxxxx, tr_y_xxxyyz_xxxxy, tr_y_xxxyyz_xxxxz, tr_y_xxxyyz_xxxyy, tr_y_xxxyyz_xxxyz, tr_y_xxxyyz_xxxzz, tr_y_xxxyyz_xxyyy, tr_y_xxxyyz_xxyyz, tr_y_xxxyyz_xxyzz, tr_y_xxxyyz_xxzzz, tr_y_xxxyyz_xyyyy, tr_y_xxxyyz_xyyyz, tr_y_xxxyyz_xyyzz, tr_y_xxxyyz_xyzzz, tr_y_xxxyyz_xzzzz, tr_y_xxxyyz_yyyyy, tr_y_xxxyyz_yyyyz, tr_y_xxxyyz_yyyzz, tr_y_xxxyyz_yyzzz, tr_y_xxxyyz_yzzzz, tr_y_xxxyyz_zzzzz, tr_y_xxyyz_yyyyz, tr_y_xxyyz_yyyzz, tr_y_xxyyz_yyzzz, tr_y_xxyyz_yzzzz, tr_y_xxyyz_zzzzz, tr_y_xyyz_yyyyz, tr_y_xyyz_yyyzz, tr_y_xyyz_yyzzz, tr_y_xyyz_yzzzz, tr_y_xyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyz_xxxxx[i] = tr_y_xxxyy_xxxxx[i] * pa_z[i];

        tr_y_xxxyyz_xxxxy[i] = tr_y_xxxyy_xxxxy[i] * pa_z[i];

        tr_y_xxxyyz_xxxxz[i] = tr_y_xxxyy_xxxx[i] * fe_0 + tr_y_xxxyy_xxxxz[i] * pa_z[i];

        tr_y_xxxyyz_xxxyy[i] = tr_y_xxxyy_xxxyy[i] * pa_z[i];

        tr_y_xxxyyz_xxxyz[i] = tr_y_xxxyy_xxxy[i] * fe_0 + tr_y_xxxyy_xxxyz[i] * pa_z[i];

        tr_y_xxxyyz_xxxzz[i] = 2.0 * tr_y_xxxyy_xxxz[i] * fe_0 + tr_y_xxxyy_xxxzz[i] * pa_z[i];

        tr_y_xxxyyz_xxyyy[i] = tr_y_xxxyy_xxyyy[i] * pa_z[i];

        tr_y_xxxyyz_xxyyz[i] = tr_y_xxxyy_xxyy[i] * fe_0 + tr_y_xxxyy_xxyyz[i] * pa_z[i];

        tr_y_xxxyyz_xxyzz[i] = 2.0 * tr_y_xxxyy_xxyz[i] * fe_0 + tr_y_xxxyy_xxyzz[i] * pa_z[i];

        tr_y_xxxyyz_xxzzz[i] = 3.0 * tr_y_xxxyy_xxzz[i] * fe_0 + tr_y_xxxyy_xxzzz[i] * pa_z[i];

        tr_y_xxxyyz_xyyyy[i] = tr_y_xxxyy_xyyyy[i] * pa_z[i];

        tr_y_xxxyyz_xyyyz[i] = tr_y_xxxyy_xyyy[i] * fe_0 + tr_y_xxxyy_xyyyz[i] * pa_z[i];

        tr_y_xxxyyz_xyyzz[i] = 2.0 * tr_y_xxxyy_xyyz[i] * fe_0 + tr_y_xxxyy_xyyzz[i] * pa_z[i];

        tr_y_xxxyyz_xyzzz[i] = 3.0 * tr_y_xxxyy_xyzz[i] * fe_0 + tr_y_xxxyy_xyzzz[i] * pa_z[i];

        tr_y_xxxyyz_xzzzz[i] = 4.0 * tr_y_xxxyy_xzzz[i] * fe_0 + tr_y_xxxyy_xzzzz[i] * pa_z[i];

        tr_y_xxxyyz_yyyyy[i] = tr_y_xxxyy_yyyyy[i] * pa_z[i];

        tr_y_xxxyyz_yyyyz[i] = 2.0 * tr_y_xyyz_yyyyz[i] * fe_0 + tr_y_xxyyz_yyyyz[i] * pa_x[i];

        tr_y_xxxyyz_yyyzz[i] = 2.0 * tr_y_xyyz_yyyzz[i] * fe_0 + tr_y_xxyyz_yyyzz[i] * pa_x[i];

        tr_y_xxxyyz_yyzzz[i] = 2.0 * tr_y_xyyz_yyzzz[i] * fe_0 + tr_y_xxyyz_yyzzz[i] * pa_x[i];

        tr_y_xxxyyz_yzzzz[i] = 2.0 * tr_y_xyyz_yzzzz[i] * fe_0 + tr_y_xxyyz_yzzzz[i] * pa_x[i];

        tr_y_xxxyyz_zzzzz[i] = 2.0 * tr_y_xyyz_zzzzz[i] * fe_0 + tr_y_xxyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 756-777 components of targeted buffer : IH

    auto tr_y_xxxyzz_xxxxx = pbuffer.data(idx_dip_ih + 756);

    auto tr_y_xxxyzz_xxxxy = pbuffer.data(idx_dip_ih + 757);

    auto tr_y_xxxyzz_xxxxz = pbuffer.data(idx_dip_ih + 758);

    auto tr_y_xxxyzz_xxxyy = pbuffer.data(idx_dip_ih + 759);

    auto tr_y_xxxyzz_xxxyz = pbuffer.data(idx_dip_ih + 760);

    auto tr_y_xxxyzz_xxxzz = pbuffer.data(idx_dip_ih + 761);

    auto tr_y_xxxyzz_xxyyy = pbuffer.data(idx_dip_ih + 762);

    auto tr_y_xxxyzz_xxyyz = pbuffer.data(idx_dip_ih + 763);

    auto tr_y_xxxyzz_xxyzz = pbuffer.data(idx_dip_ih + 764);

    auto tr_y_xxxyzz_xxzzz = pbuffer.data(idx_dip_ih + 765);

    auto tr_y_xxxyzz_xyyyy = pbuffer.data(idx_dip_ih + 766);

    auto tr_y_xxxyzz_xyyyz = pbuffer.data(idx_dip_ih + 767);

    auto tr_y_xxxyzz_xyyzz = pbuffer.data(idx_dip_ih + 768);

    auto tr_y_xxxyzz_xyzzz = pbuffer.data(idx_dip_ih + 769);

    auto tr_y_xxxyzz_xzzzz = pbuffer.data(idx_dip_ih + 770);

    auto tr_y_xxxyzz_yyyyy = pbuffer.data(idx_dip_ih + 771);

    auto tr_y_xxxyzz_yyyyz = pbuffer.data(idx_dip_ih + 772);

    auto tr_y_xxxyzz_yyyzz = pbuffer.data(idx_dip_ih + 773);

    auto tr_y_xxxyzz_yyzzz = pbuffer.data(idx_dip_ih + 774);

    auto tr_y_xxxyzz_yzzzz = pbuffer.data(idx_dip_ih + 775);

    auto tr_y_xxxyzz_zzzzz = pbuffer.data(idx_dip_ih + 776);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxy_xxxxy, tr_y_xxxy_xxxyy, tr_y_xxxy_xxyyy, tr_y_xxxy_xyyyy, tr_y_xxxyz_xxxxy, tr_y_xxxyz_xxxyy, tr_y_xxxyz_xxyyy, tr_y_xxxyz_xyyyy, tr_y_xxxyzz_xxxxx, tr_y_xxxyzz_xxxxy, tr_y_xxxyzz_xxxxz, tr_y_xxxyzz_xxxyy, tr_y_xxxyzz_xxxyz, tr_y_xxxyzz_xxxzz, tr_y_xxxyzz_xxyyy, tr_y_xxxyzz_xxyyz, tr_y_xxxyzz_xxyzz, tr_y_xxxyzz_xxzzz, tr_y_xxxyzz_xyyyy, tr_y_xxxyzz_xyyyz, tr_y_xxxyzz_xyyzz, tr_y_xxxyzz_xyzzz, tr_y_xxxyzz_xzzzz, tr_y_xxxyzz_yyyyy, tr_y_xxxyzz_yyyyz, tr_y_xxxyzz_yyyzz, tr_y_xxxyzz_yyzzz, tr_y_xxxyzz_yzzzz, tr_y_xxxyzz_zzzzz, tr_y_xxxzz_xxxxx, tr_y_xxxzz_xxxxz, tr_y_xxxzz_xxxzz, tr_y_xxxzz_xxzzz, tr_y_xxxzz_xzzzz, tr_y_xxyzz_xxxyz, tr_y_xxyzz_xxyyz, tr_y_xxyzz_xxyz, tr_y_xxyzz_xxyzz, tr_y_xxyzz_xyyyz, tr_y_xxyzz_xyyz, tr_y_xxyzz_xyyzz, tr_y_xxyzz_xyzz, tr_y_xxyzz_xyzzz, tr_y_xxyzz_yyyyy, tr_y_xxyzz_yyyyz, tr_y_xxyzz_yyyz, tr_y_xxyzz_yyyzz, tr_y_xxyzz_yyzz, tr_y_xxyzz_yyzzz, tr_y_xxyzz_yzzz, tr_y_xxyzz_yzzzz, tr_y_xxyzz_zzzzz, tr_y_xyzz_xxxyz, tr_y_xyzz_xxyyz, tr_y_xyzz_xxyzz, tr_y_xyzz_xyyyz, tr_y_xyzz_xyyzz, tr_y_xyzz_xyzzz, tr_y_xyzz_yyyyy, tr_y_xyzz_yyyyz, tr_y_xyzz_yyyzz, tr_y_xyzz_yyzzz, tr_y_xyzz_yzzzz, tr_y_xyzz_zzzzz, ts_xxxzz_xxxxx, ts_xxxzz_xxxxz, ts_xxxzz_xxxzz, ts_xxxzz_xxzzz, ts_xxxzz_xzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyzz_xxxxx[i] = ts_xxxzz_xxxxx[i] * fe_0 + tr_y_xxxzz_xxxxx[i] * pa_y[i];

        tr_y_xxxyzz_xxxxy[i] = tr_y_xxxy_xxxxy[i] * fe_0 + tr_y_xxxyz_xxxxy[i] * pa_z[i];

        tr_y_xxxyzz_xxxxz[i] = ts_xxxzz_xxxxz[i] * fe_0 + tr_y_xxxzz_xxxxz[i] * pa_y[i];

        tr_y_xxxyzz_xxxyy[i] = tr_y_xxxy_xxxyy[i] * fe_0 + tr_y_xxxyz_xxxyy[i] * pa_z[i];

        tr_y_xxxyzz_xxxyz[i] = 2.0 * tr_y_xyzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xxyzz_xxyz[i] * fe_0 + tr_y_xxyzz_xxxyz[i] * pa_x[i];

        tr_y_xxxyzz_xxxzz[i] = ts_xxxzz_xxxzz[i] * fe_0 + tr_y_xxxzz_xxxzz[i] * pa_y[i];

        tr_y_xxxyzz_xxyyy[i] = tr_y_xxxy_xxyyy[i] * fe_0 + tr_y_xxxyz_xxyyy[i] * pa_z[i];

        tr_y_xxxyzz_xxyyz[i] = 2.0 * tr_y_xyzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xxyzz_xyyz[i] * fe_0 + tr_y_xxyzz_xxyyz[i] * pa_x[i];

        tr_y_xxxyzz_xxyzz[i] = 2.0 * tr_y_xyzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xxyzz_xyzz[i] * fe_0 + tr_y_xxyzz_xxyzz[i] * pa_x[i];

        tr_y_xxxyzz_xxzzz[i] = ts_xxxzz_xxzzz[i] * fe_0 + tr_y_xxxzz_xxzzz[i] * pa_y[i];

        tr_y_xxxyzz_xyyyy[i] = tr_y_xxxy_xyyyy[i] * fe_0 + tr_y_xxxyz_xyyyy[i] * pa_z[i];

        tr_y_xxxyzz_xyyyz[i] = 2.0 * tr_y_xyzz_xyyyz[i] * fe_0 + tr_y_xxyzz_yyyz[i] * fe_0 + tr_y_xxyzz_xyyyz[i] * pa_x[i];

        tr_y_xxxyzz_xyyzz[i] = 2.0 * tr_y_xyzz_xyyzz[i] * fe_0 + tr_y_xxyzz_yyzz[i] * fe_0 + tr_y_xxyzz_xyyzz[i] * pa_x[i];

        tr_y_xxxyzz_xyzzz[i] = 2.0 * tr_y_xyzz_xyzzz[i] * fe_0 + tr_y_xxyzz_yzzz[i] * fe_0 + tr_y_xxyzz_xyzzz[i] * pa_x[i];

        tr_y_xxxyzz_xzzzz[i] = ts_xxxzz_xzzzz[i] * fe_0 + tr_y_xxxzz_xzzzz[i] * pa_y[i];

        tr_y_xxxyzz_yyyyy[i] = 2.0 * tr_y_xyzz_yyyyy[i] * fe_0 + tr_y_xxyzz_yyyyy[i] * pa_x[i];

        tr_y_xxxyzz_yyyyz[i] = 2.0 * tr_y_xyzz_yyyyz[i] * fe_0 + tr_y_xxyzz_yyyyz[i] * pa_x[i];

        tr_y_xxxyzz_yyyzz[i] = 2.0 * tr_y_xyzz_yyyzz[i] * fe_0 + tr_y_xxyzz_yyyzz[i] * pa_x[i];

        tr_y_xxxyzz_yyzzz[i] = 2.0 * tr_y_xyzz_yyzzz[i] * fe_0 + tr_y_xxyzz_yyzzz[i] * pa_x[i];

        tr_y_xxxyzz_yzzzz[i] = 2.0 * tr_y_xyzz_yzzzz[i] * fe_0 + tr_y_xxyzz_yzzzz[i] * pa_x[i];

        tr_y_xxxyzz_zzzzz[i] = 2.0 * tr_y_xyzz_zzzzz[i] * fe_0 + tr_y_xxyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 777-798 components of targeted buffer : IH

    auto tr_y_xxxzzz_xxxxx = pbuffer.data(idx_dip_ih + 777);

    auto tr_y_xxxzzz_xxxxy = pbuffer.data(idx_dip_ih + 778);

    auto tr_y_xxxzzz_xxxxz = pbuffer.data(idx_dip_ih + 779);

    auto tr_y_xxxzzz_xxxyy = pbuffer.data(idx_dip_ih + 780);

    auto tr_y_xxxzzz_xxxyz = pbuffer.data(idx_dip_ih + 781);

    auto tr_y_xxxzzz_xxxzz = pbuffer.data(idx_dip_ih + 782);

    auto tr_y_xxxzzz_xxyyy = pbuffer.data(idx_dip_ih + 783);

    auto tr_y_xxxzzz_xxyyz = pbuffer.data(idx_dip_ih + 784);

    auto tr_y_xxxzzz_xxyzz = pbuffer.data(idx_dip_ih + 785);

    auto tr_y_xxxzzz_xxzzz = pbuffer.data(idx_dip_ih + 786);

    auto tr_y_xxxzzz_xyyyy = pbuffer.data(idx_dip_ih + 787);

    auto tr_y_xxxzzz_xyyyz = pbuffer.data(idx_dip_ih + 788);

    auto tr_y_xxxzzz_xyyzz = pbuffer.data(idx_dip_ih + 789);

    auto tr_y_xxxzzz_xyzzz = pbuffer.data(idx_dip_ih + 790);

    auto tr_y_xxxzzz_xzzzz = pbuffer.data(idx_dip_ih + 791);

    auto tr_y_xxxzzz_yyyyy = pbuffer.data(idx_dip_ih + 792);

    auto tr_y_xxxzzz_yyyyz = pbuffer.data(idx_dip_ih + 793);

    auto tr_y_xxxzzz_yyyzz = pbuffer.data(idx_dip_ih + 794);

    auto tr_y_xxxzzz_yyzzz = pbuffer.data(idx_dip_ih + 795);

    auto tr_y_xxxzzz_yzzzz = pbuffer.data(idx_dip_ih + 796);

    auto tr_y_xxxzzz_zzzzz = pbuffer.data(idx_dip_ih + 797);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxz_xxxxx, tr_y_xxxz_xxxxy, tr_y_xxxz_xxxyy, tr_y_xxxz_xxyyy, tr_y_xxxz_xyyyy, tr_y_xxxzz_xxxxx, tr_y_xxxzz_xxxxy, tr_y_xxxzz_xxxyy, tr_y_xxxzz_xxyyy, tr_y_xxxzz_xyyyy, tr_y_xxxzzz_xxxxx, tr_y_xxxzzz_xxxxy, tr_y_xxxzzz_xxxxz, tr_y_xxxzzz_xxxyy, tr_y_xxxzzz_xxxyz, tr_y_xxxzzz_xxxzz, tr_y_xxxzzz_xxyyy, tr_y_xxxzzz_xxyyz, tr_y_xxxzzz_xxyzz, tr_y_xxxzzz_xxzzz, tr_y_xxxzzz_xyyyy, tr_y_xxxzzz_xyyyz, tr_y_xxxzzz_xyyzz, tr_y_xxxzzz_xyzzz, tr_y_xxxzzz_xzzzz, tr_y_xxxzzz_yyyyy, tr_y_xxxzzz_yyyyz, tr_y_xxxzzz_yyyzz, tr_y_xxxzzz_yyzzz, tr_y_xxxzzz_yzzzz, tr_y_xxxzzz_zzzzz, tr_y_xxzzz_xxxxz, tr_y_xxzzz_xxxyz, tr_y_xxzzz_xxxz, tr_y_xxzzz_xxxzz, tr_y_xxzzz_xxyyz, tr_y_xxzzz_xxyz, tr_y_xxzzz_xxyzz, tr_y_xxzzz_xxzz, tr_y_xxzzz_xxzzz, tr_y_xxzzz_xyyyz, tr_y_xxzzz_xyyz, tr_y_xxzzz_xyyzz, tr_y_xxzzz_xyzz, tr_y_xxzzz_xyzzz, tr_y_xxzzz_xzzz, tr_y_xxzzz_xzzzz, tr_y_xxzzz_yyyyy, tr_y_xxzzz_yyyyz, tr_y_xxzzz_yyyz, tr_y_xxzzz_yyyzz, tr_y_xxzzz_yyzz, tr_y_xxzzz_yyzzz, tr_y_xxzzz_yzzz, tr_y_xxzzz_yzzzz, tr_y_xxzzz_zzzz, tr_y_xxzzz_zzzzz, tr_y_xzzz_xxxxz, tr_y_xzzz_xxxyz, tr_y_xzzz_xxxzz, tr_y_xzzz_xxyyz, tr_y_xzzz_xxyzz, tr_y_xzzz_xxzzz, tr_y_xzzz_xyyyz, tr_y_xzzz_xyyzz, tr_y_xzzz_xyzzz, tr_y_xzzz_xzzzz, tr_y_xzzz_yyyyy, tr_y_xzzz_yyyyz, tr_y_xzzz_yyyzz, tr_y_xzzz_yyzzz, tr_y_xzzz_yzzzz, tr_y_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzzz_xxxxx[i] = 2.0 * tr_y_xxxz_xxxxx[i] * fe_0 + tr_y_xxxzz_xxxxx[i] * pa_z[i];

        tr_y_xxxzzz_xxxxy[i] = 2.0 * tr_y_xxxz_xxxxy[i] * fe_0 + tr_y_xxxzz_xxxxy[i] * pa_z[i];

        tr_y_xxxzzz_xxxxz[i] = 2.0 * tr_y_xzzz_xxxxz[i] * fe_0 + 4.0 * tr_y_xxzzz_xxxz[i] * fe_0 + tr_y_xxzzz_xxxxz[i] * pa_x[i];

        tr_y_xxxzzz_xxxyy[i] = 2.0 * tr_y_xxxz_xxxyy[i] * fe_0 + tr_y_xxxzz_xxxyy[i] * pa_z[i];

        tr_y_xxxzzz_xxxyz[i] = 2.0 * tr_y_xzzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xxzzz_xxyz[i] * fe_0 + tr_y_xxzzz_xxxyz[i] * pa_x[i];

        tr_y_xxxzzz_xxxzz[i] = 2.0 * tr_y_xzzz_xxxzz[i] * fe_0 + 3.0 * tr_y_xxzzz_xxzz[i] * fe_0 + tr_y_xxzzz_xxxzz[i] * pa_x[i];

        tr_y_xxxzzz_xxyyy[i] = 2.0 * tr_y_xxxz_xxyyy[i] * fe_0 + tr_y_xxxzz_xxyyy[i] * pa_z[i];

        tr_y_xxxzzz_xxyyz[i] = 2.0 * tr_y_xzzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xxzzz_xyyz[i] * fe_0 + tr_y_xxzzz_xxyyz[i] * pa_x[i];

        tr_y_xxxzzz_xxyzz[i] = 2.0 * tr_y_xzzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xxzzz_xyzz[i] * fe_0 + tr_y_xxzzz_xxyzz[i] * pa_x[i];

        tr_y_xxxzzz_xxzzz[i] = 2.0 * tr_y_xzzz_xxzzz[i] * fe_0 + 2.0 * tr_y_xxzzz_xzzz[i] * fe_0 + tr_y_xxzzz_xxzzz[i] * pa_x[i];

        tr_y_xxxzzz_xyyyy[i] = 2.0 * tr_y_xxxz_xyyyy[i] * fe_0 + tr_y_xxxzz_xyyyy[i] * pa_z[i];

        tr_y_xxxzzz_xyyyz[i] = 2.0 * tr_y_xzzz_xyyyz[i] * fe_0 + tr_y_xxzzz_yyyz[i] * fe_0 + tr_y_xxzzz_xyyyz[i] * pa_x[i];

        tr_y_xxxzzz_xyyzz[i] = 2.0 * tr_y_xzzz_xyyzz[i] * fe_0 + tr_y_xxzzz_yyzz[i] * fe_0 + tr_y_xxzzz_xyyzz[i] * pa_x[i];

        tr_y_xxxzzz_xyzzz[i] = 2.0 * tr_y_xzzz_xyzzz[i] * fe_0 + tr_y_xxzzz_yzzz[i] * fe_0 + tr_y_xxzzz_xyzzz[i] * pa_x[i];

        tr_y_xxxzzz_xzzzz[i] = 2.0 * tr_y_xzzz_xzzzz[i] * fe_0 + tr_y_xxzzz_zzzz[i] * fe_0 + tr_y_xxzzz_xzzzz[i] * pa_x[i];

        tr_y_xxxzzz_yyyyy[i] = 2.0 * tr_y_xzzz_yyyyy[i] * fe_0 + tr_y_xxzzz_yyyyy[i] * pa_x[i];

        tr_y_xxxzzz_yyyyz[i] = 2.0 * tr_y_xzzz_yyyyz[i] * fe_0 + tr_y_xxzzz_yyyyz[i] * pa_x[i];

        tr_y_xxxzzz_yyyzz[i] = 2.0 * tr_y_xzzz_yyyzz[i] * fe_0 + tr_y_xxzzz_yyyzz[i] * pa_x[i];

        tr_y_xxxzzz_yyzzz[i] = 2.0 * tr_y_xzzz_yyzzz[i] * fe_0 + tr_y_xxzzz_yyzzz[i] * pa_x[i];

        tr_y_xxxzzz_yzzzz[i] = 2.0 * tr_y_xzzz_yzzzz[i] * fe_0 + tr_y_xxzzz_yzzzz[i] * pa_x[i];

        tr_y_xxxzzz_zzzzz[i] = 2.0 * tr_y_xzzz_zzzzz[i] * fe_0 + tr_y_xxzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 798-819 components of targeted buffer : IH

    auto tr_y_xxyyyy_xxxxx = pbuffer.data(idx_dip_ih + 798);

    auto tr_y_xxyyyy_xxxxy = pbuffer.data(idx_dip_ih + 799);

    auto tr_y_xxyyyy_xxxxz = pbuffer.data(idx_dip_ih + 800);

    auto tr_y_xxyyyy_xxxyy = pbuffer.data(idx_dip_ih + 801);

    auto tr_y_xxyyyy_xxxyz = pbuffer.data(idx_dip_ih + 802);

    auto tr_y_xxyyyy_xxxzz = pbuffer.data(idx_dip_ih + 803);

    auto tr_y_xxyyyy_xxyyy = pbuffer.data(idx_dip_ih + 804);

    auto tr_y_xxyyyy_xxyyz = pbuffer.data(idx_dip_ih + 805);

    auto tr_y_xxyyyy_xxyzz = pbuffer.data(idx_dip_ih + 806);

    auto tr_y_xxyyyy_xxzzz = pbuffer.data(idx_dip_ih + 807);

    auto tr_y_xxyyyy_xyyyy = pbuffer.data(idx_dip_ih + 808);

    auto tr_y_xxyyyy_xyyyz = pbuffer.data(idx_dip_ih + 809);

    auto tr_y_xxyyyy_xyyzz = pbuffer.data(idx_dip_ih + 810);

    auto tr_y_xxyyyy_xyzzz = pbuffer.data(idx_dip_ih + 811);

    auto tr_y_xxyyyy_xzzzz = pbuffer.data(idx_dip_ih + 812);

    auto tr_y_xxyyyy_yyyyy = pbuffer.data(idx_dip_ih + 813);

    auto tr_y_xxyyyy_yyyyz = pbuffer.data(idx_dip_ih + 814);

    auto tr_y_xxyyyy_yyyzz = pbuffer.data(idx_dip_ih + 815);

    auto tr_y_xxyyyy_yyzzz = pbuffer.data(idx_dip_ih + 816);

    auto tr_y_xxyyyy_yzzzz = pbuffer.data(idx_dip_ih + 817);

    auto tr_y_xxyyyy_zzzzz = pbuffer.data(idx_dip_ih + 818);

    #pragma omp simd aligned(pa_x, tr_y_xxyyyy_xxxxx, tr_y_xxyyyy_xxxxy, tr_y_xxyyyy_xxxxz, tr_y_xxyyyy_xxxyy, tr_y_xxyyyy_xxxyz, tr_y_xxyyyy_xxxzz, tr_y_xxyyyy_xxyyy, tr_y_xxyyyy_xxyyz, tr_y_xxyyyy_xxyzz, tr_y_xxyyyy_xxzzz, tr_y_xxyyyy_xyyyy, tr_y_xxyyyy_xyyyz, tr_y_xxyyyy_xyyzz, tr_y_xxyyyy_xyzzz, tr_y_xxyyyy_xzzzz, tr_y_xxyyyy_yyyyy, tr_y_xxyyyy_yyyyz, tr_y_xxyyyy_yyyzz, tr_y_xxyyyy_yyzzz, tr_y_xxyyyy_yzzzz, tr_y_xxyyyy_zzzzz, tr_y_xyyyy_xxxx, tr_y_xyyyy_xxxxx, tr_y_xyyyy_xxxxy, tr_y_xyyyy_xxxxz, tr_y_xyyyy_xxxy, tr_y_xyyyy_xxxyy, tr_y_xyyyy_xxxyz, tr_y_xyyyy_xxxz, tr_y_xyyyy_xxxzz, tr_y_xyyyy_xxyy, tr_y_xyyyy_xxyyy, tr_y_xyyyy_xxyyz, tr_y_xyyyy_xxyz, tr_y_xyyyy_xxyzz, tr_y_xyyyy_xxzz, tr_y_xyyyy_xxzzz, tr_y_xyyyy_xyyy, tr_y_xyyyy_xyyyy, tr_y_xyyyy_xyyyz, tr_y_xyyyy_xyyz, tr_y_xyyyy_xyyzz, tr_y_xyyyy_xyzz, tr_y_xyyyy_xyzzz, tr_y_xyyyy_xzzz, tr_y_xyyyy_xzzzz, tr_y_xyyyy_yyyy, tr_y_xyyyy_yyyyy, tr_y_xyyyy_yyyyz, tr_y_xyyyy_yyyz, tr_y_xyyyy_yyyzz, tr_y_xyyyy_yyzz, tr_y_xyyyy_yyzzz, tr_y_xyyyy_yzzz, tr_y_xyyyy_yzzzz, tr_y_xyyyy_zzzz, tr_y_xyyyy_zzzzz, tr_y_yyyy_xxxxx, tr_y_yyyy_xxxxy, tr_y_yyyy_xxxxz, tr_y_yyyy_xxxyy, tr_y_yyyy_xxxyz, tr_y_yyyy_xxxzz, tr_y_yyyy_xxyyy, tr_y_yyyy_xxyyz, tr_y_yyyy_xxyzz, tr_y_yyyy_xxzzz, tr_y_yyyy_xyyyy, tr_y_yyyy_xyyyz, tr_y_yyyy_xyyzz, tr_y_yyyy_xyzzz, tr_y_yyyy_xzzzz, tr_y_yyyy_yyyyy, tr_y_yyyy_yyyyz, tr_y_yyyy_yyyzz, tr_y_yyyy_yyzzz, tr_y_yyyy_yzzzz, tr_y_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyy_xxxxx[i] = tr_y_yyyy_xxxxx[i] * fe_0 + 5.0 * tr_y_xyyyy_xxxx[i] * fe_0 + tr_y_xyyyy_xxxxx[i] * pa_x[i];

        tr_y_xxyyyy_xxxxy[i] = tr_y_yyyy_xxxxy[i] * fe_0 + 4.0 * tr_y_xyyyy_xxxy[i] * fe_0 + tr_y_xyyyy_xxxxy[i] * pa_x[i];

        tr_y_xxyyyy_xxxxz[i] = tr_y_yyyy_xxxxz[i] * fe_0 + 4.0 * tr_y_xyyyy_xxxz[i] * fe_0 + tr_y_xyyyy_xxxxz[i] * pa_x[i];

        tr_y_xxyyyy_xxxyy[i] = tr_y_yyyy_xxxyy[i] * fe_0 + 3.0 * tr_y_xyyyy_xxyy[i] * fe_0 + tr_y_xyyyy_xxxyy[i] * pa_x[i];

        tr_y_xxyyyy_xxxyz[i] = tr_y_yyyy_xxxyz[i] * fe_0 + 3.0 * tr_y_xyyyy_xxyz[i] * fe_0 + tr_y_xyyyy_xxxyz[i] * pa_x[i];

        tr_y_xxyyyy_xxxzz[i] = tr_y_yyyy_xxxzz[i] * fe_0 + 3.0 * tr_y_xyyyy_xxzz[i] * fe_0 + tr_y_xyyyy_xxxzz[i] * pa_x[i];

        tr_y_xxyyyy_xxyyy[i] = tr_y_yyyy_xxyyy[i] * fe_0 + 2.0 * tr_y_xyyyy_xyyy[i] * fe_0 + tr_y_xyyyy_xxyyy[i] * pa_x[i];

        tr_y_xxyyyy_xxyyz[i] = tr_y_yyyy_xxyyz[i] * fe_0 + 2.0 * tr_y_xyyyy_xyyz[i] * fe_0 + tr_y_xyyyy_xxyyz[i] * pa_x[i];

        tr_y_xxyyyy_xxyzz[i] = tr_y_yyyy_xxyzz[i] * fe_0 + 2.0 * tr_y_xyyyy_xyzz[i] * fe_0 + tr_y_xyyyy_xxyzz[i] * pa_x[i];

        tr_y_xxyyyy_xxzzz[i] = tr_y_yyyy_xxzzz[i] * fe_0 + 2.0 * tr_y_xyyyy_xzzz[i] * fe_0 + tr_y_xyyyy_xxzzz[i] * pa_x[i];

        tr_y_xxyyyy_xyyyy[i] = tr_y_yyyy_xyyyy[i] * fe_0 + tr_y_xyyyy_yyyy[i] * fe_0 + tr_y_xyyyy_xyyyy[i] * pa_x[i];

        tr_y_xxyyyy_xyyyz[i] = tr_y_yyyy_xyyyz[i] * fe_0 + tr_y_xyyyy_yyyz[i] * fe_0 + tr_y_xyyyy_xyyyz[i] * pa_x[i];

        tr_y_xxyyyy_xyyzz[i] = tr_y_yyyy_xyyzz[i] * fe_0 + tr_y_xyyyy_yyzz[i] * fe_0 + tr_y_xyyyy_xyyzz[i] * pa_x[i];

        tr_y_xxyyyy_xyzzz[i] = tr_y_yyyy_xyzzz[i] * fe_0 + tr_y_xyyyy_yzzz[i] * fe_0 + tr_y_xyyyy_xyzzz[i] * pa_x[i];

        tr_y_xxyyyy_xzzzz[i] = tr_y_yyyy_xzzzz[i] * fe_0 + tr_y_xyyyy_zzzz[i] * fe_0 + tr_y_xyyyy_xzzzz[i] * pa_x[i];

        tr_y_xxyyyy_yyyyy[i] = tr_y_yyyy_yyyyy[i] * fe_0 + tr_y_xyyyy_yyyyy[i] * pa_x[i];

        tr_y_xxyyyy_yyyyz[i] = tr_y_yyyy_yyyyz[i] * fe_0 + tr_y_xyyyy_yyyyz[i] * pa_x[i];

        tr_y_xxyyyy_yyyzz[i] = tr_y_yyyy_yyyzz[i] * fe_0 + tr_y_xyyyy_yyyzz[i] * pa_x[i];

        tr_y_xxyyyy_yyzzz[i] = tr_y_yyyy_yyzzz[i] * fe_0 + tr_y_xyyyy_yyzzz[i] * pa_x[i];

        tr_y_xxyyyy_yzzzz[i] = tr_y_yyyy_yzzzz[i] * fe_0 + tr_y_xyyyy_yzzzz[i] * pa_x[i];

        tr_y_xxyyyy_zzzzz[i] = tr_y_yyyy_zzzzz[i] * fe_0 + tr_y_xyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 819-840 components of targeted buffer : IH

    auto tr_y_xxyyyz_xxxxx = pbuffer.data(idx_dip_ih + 819);

    auto tr_y_xxyyyz_xxxxy = pbuffer.data(idx_dip_ih + 820);

    auto tr_y_xxyyyz_xxxxz = pbuffer.data(idx_dip_ih + 821);

    auto tr_y_xxyyyz_xxxyy = pbuffer.data(idx_dip_ih + 822);

    auto tr_y_xxyyyz_xxxyz = pbuffer.data(idx_dip_ih + 823);

    auto tr_y_xxyyyz_xxxzz = pbuffer.data(idx_dip_ih + 824);

    auto tr_y_xxyyyz_xxyyy = pbuffer.data(idx_dip_ih + 825);

    auto tr_y_xxyyyz_xxyyz = pbuffer.data(idx_dip_ih + 826);

    auto tr_y_xxyyyz_xxyzz = pbuffer.data(idx_dip_ih + 827);

    auto tr_y_xxyyyz_xxzzz = pbuffer.data(idx_dip_ih + 828);

    auto tr_y_xxyyyz_xyyyy = pbuffer.data(idx_dip_ih + 829);

    auto tr_y_xxyyyz_xyyyz = pbuffer.data(idx_dip_ih + 830);

    auto tr_y_xxyyyz_xyyzz = pbuffer.data(idx_dip_ih + 831);

    auto tr_y_xxyyyz_xyzzz = pbuffer.data(idx_dip_ih + 832);

    auto tr_y_xxyyyz_xzzzz = pbuffer.data(idx_dip_ih + 833);

    auto tr_y_xxyyyz_yyyyy = pbuffer.data(idx_dip_ih + 834);

    auto tr_y_xxyyyz_yyyyz = pbuffer.data(idx_dip_ih + 835);

    auto tr_y_xxyyyz_yyyzz = pbuffer.data(idx_dip_ih + 836);

    auto tr_y_xxyyyz_yyzzz = pbuffer.data(idx_dip_ih + 837);

    auto tr_y_xxyyyz_yzzzz = pbuffer.data(idx_dip_ih + 838);

    auto tr_y_xxyyyz_zzzzz = pbuffer.data(idx_dip_ih + 839);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyyy_xxxx, tr_y_xxyyy_xxxxx, tr_y_xxyyy_xxxxy, tr_y_xxyyy_xxxxz, tr_y_xxyyy_xxxy, tr_y_xxyyy_xxxyy, tr_y_xxyyy_xxxyz, tr_y_xxyyy_xxxz, tr_y_xxyyy_xxxzz, tr_y_xxyyy_xxyy, tr_y_xxyyy_xxyyy, tr_y_xxyyy_xxyyz, tr_y_xxyyy_xxyz, tr_y_xxyyy_xxyzz, tr_y_xxyyy_xxzz, tr_y_xxyyy_xxzzz, tr_y_xxyyy_xyyy, tr_y_xxyyy_xyyyy, tr_y_xxyyy_xyyyz, tr_y_xxyyy_xyyz, tr_y_xxyyy_xyyzz, tr_y_xxyyy_xyzz, tr_y_xxyyy_xyzzz, tr_y_xxyyy_xzzz, tr_y_xxyyy_xzzzz, tr_y_xxyyy_yyyyy, tr_y_xxyyyz_xxxxx, tr_y_xxyyyz_xxxxy, tr_y_xxyyyz_xxxxz, tr_y_xxyyyz_xxxyy, tr_y_xxyyyz_xxxyz, tr_y_xxyyyz_xxxzz, tr_y_xxyyyz_xxyyy, tr_y_xxyyyz_xxyyz, tr_y_xxyyyz_xxyzz, tr_y_xxyyyz_xxzzz, tr_y_xxyyyz_xyyyy, tr_y_xxyyyz_xyyyz, tr_y_xxyyyz_xyyzz, tr_y_xxyyyz_xyzzz, tr_y_xxyyyz_xzzzz, tr_y_xxyyyz_yyyyy, tr_y_xxyyyz_yyyyz, tr_y_xxyyyz_yyyzz, tr_y_xxyyyz_yyzzz, tr_y_xxyyyz_yzzzz, tr_y_xxyyyz_zzzzz, tr_y_xyyyz_yyyyz, tr_y_xyyyz_yyyzz, tr_y_xyyyz_yyzzz, tr_y_xyyyz_yzzzz, tr_y_xyyyz_zzzzz, tr_y_yyyz_yyyyz, tr_y_yyyz_yyyzz, tr_y_yyyz_yyzzz, tr_y_yyyz_yzzzz, tr_y_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyz_xxxxx[i] = tr_y_xxyyy_xxxxx[i] * pa_z[i];

        tr_y_xxyyyz_xxxxy[i] = tr_y_xxyyy_xxxxy[i] * pa_z[i];

        tr_y_xxyyyz_xxxxz[i] = tr_y_xxyyy_xxxx[i] * fe_0 + tr_y_xxyyy_xxxxz[i] * pa_z[i];

        tr_y_xxyyyz_xxxyy[i] = tr_y_xxyyy_xxxyy[i] * pa_z[i];

        tr_y_xxyyyz_xxxyz[i] = tr_y_xxyyy_xxxy[i] * fe_0 + tr_y_xxyyy_xxxyz[i] * pa_z[i];

        tr_y_xxyyyz_xxxzz[i] = 2.0 * tr_y_xxyyy_xxxz[i] * fe_0 + tr_y_xxyyy_xxxzz[i] * pa_z[i];

        tr_y_xxyyyz_xxyyy[i] = tr_y_xxyyy_xxyyy[i] * pa_z[i];

        tr_y_xxyyyz_xxyyz[i] = tr_y_xxyyy_xxyy[i] * fe_0 + tr_y_xxyyy_xxyyz[i] * pa_z[i];

        tr_y_xxyyyz_xxyzz[i] = 2.0 * tr_y_xxyyy_xxyz[i] * fe_0 + tr_y_xxyyy_xxyzz[i] * pa_z[i];

        tr_y_xxyyyz_xxzzz[i] = 3.0 * tr_y_xxyyy_xxzz[i] * fe_0 + tr_y_xxyyy_xxzzz[i] * pa_z[i];

        tr_y_xxyyyz_xyyyy[i] = tr_y_xxyyy_xyyyy[i] * pa_z[i];

        tr_y_xxyyyz_xyyyz[i] = tr_y_xxyyy_xyyy[i] * fe_0 + tr_y_xxyyy_xyyyz[i] * pa_z[i];

        tr_y_xxyyyz_xyyzz[i] = 2.0 * tr_y_xxyyy_xyyz[i] * fe_0 + tr_y_xxyyy_xyyzz[i] * pa_z[i];

        tr_y_xxyyyz_xyzzz[i] = 3.0 * tr_y_xxyyy_xyzz[i] * fe_0 + tr_y_xxyyy_xyzzz[i] * pa_z[i];

        tr_y_xxyyyz_xzzzz[i] = 4.0 * tr_y_xxyyy_xzzz[i] * fe_0 + tr_y_xxyyy_xzzzz[i] * pa_z[i];

        tr_y_xxyyyz_yyyyy[i] = tr_y_xxyyy_yyyyy[i] * pa_z[i];

        tr_y_xxyyyz_yyyyz[i] = tr_y_yyyz_yyyyz[i] * fe_0 + tr_y_xyyyz_yyyyz[i] * pa_x[i];

        tr_y_xxyyyz_yyyzz[i] = tr_y_yyyz_yyyzz[i] * fe_0 + tr_y_xyyyz_yyyzz[i] * pa_x[i];

        tr_y_xxyyyz_yyzzz[i] = tr_y_yyyz_yyzzz[i] * fe_0 + tr_y_xyyyz_yyzzz[i] * pa_x[i];

        tr_y_xxyyyz_yzzzz[i] = tr_y_yyyz_yzzzz[i] * fe_0 + tr_y_xyyyz_yzzzz[i] * pa_x[i];

        tr_y_xxyyyz_zzzzz[i] = tr_y_yyyz_zzzzz[i] * fe_0 + tr_y_xyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 840-861 components of targeted buffer : IH

    auto tr_y_xxyyzz_xxxxx = pbuffer.data(idx_dip_ih + 840);

    auto tr_y_xxyyzz_xxxxy = pbuffer.data(idx_dip_ih + 841);

    auto tr_y_xxyyzz_xxxxz = pbuffer.data(idx_dip_ih + 842);

    auto tr_y_xxyyzz_xxxyy = pbuffer.data(idx_dip_ih + 843);

    auto tr_y_xxyyzz_xxxyz = pbuffer.data(idx_dip_ih + 844);

    auto tr_y_xxyyzz_xxxzz = pbuffer.data(idx_dip_ih + 845);

    auto tr_y_xxyyzz_xxyyy = pbuffer.data(idx_dip_ih + 846);

    auto tr_y_xxyyzz_xxyyz = pbuffer.data(idx_dip_ih + 847);

    auto tr_y_xxyyzz_xxyzz = pbuffer.data(idx_dip_ih + 848);

    auto tr_y_xxyyzz_xxzzz = pbuffer.data(idx_dip_ih + 849);

    auto tr_y_xxyyzz_xyyyy = pbuffer.data(idx_dip_ih + 850);

    auto tr_y_xxyyzz_xyyyz = pbuffer.data(idx_dip_ih + 851);

    auto tr_y_xxyyzz_xyyzz = pbuffer.data(idx_dip_ih + 852);

    auto tr_y_xxyyzz_xyzzz = pbuffer.data(idx_dip_ih + 853);

    auto tr_y_xxyyzz_xzzzz = pbuffer.data(idx_dip_ih + 854);

    auto tr_y_xxyyzz_yyyyy = pbuffer.data(idx_dip_ih + 855);

    auto tr_y_xxyyzz_yyyyz = pbuffer.data(idx_dip_ih + 856);

    auto tr_y_xxyyzz_yyyzz = pbuffer.data(idx_dip_ih + 857);

    auto tr_y_xxyyzz_yyzzz = pbuffer.data(idx_dip_ih + 858);

    auto tr_y_xxyyzz_yzzzz = pbuffer.data(idx_dip_ih + 859);

    auto tr_y_xxyyzz_zzzzz = pbuffer.data(idx_dip_ih + 860);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyy_xxxxx, tr_y_xxyy_xxxxy, tr_y_xxyy_xxxyy, tr_y_xxyy_xxyyy, tr_y_xxyy_xyyyy, tr_y_xxyyz_xxxxx, tr_y_xxyyz_xxxxy, tr_y_xxyyz_xxxyy, tr_y_xxyyz_xxyyy, tr_y_xxyyz_xyyyy, tr_y_xxyyzz_xxxxx, tr_y_xxyyzz_xxxxy, tr_y_xxyyzz_xxxxz, tr_y_xxyyzz_xxxyy, tr_y_xxyyzz_xxxyz, tr_y_xxyyzz_xxxzz, tr_y_xxyyzz_xxyyy, tr_y_xxyyzz_xxyyz, tr_y_xxyyzz_xxyzz, tr_y_xxyyzz_xxzzz, tr_y_xxyyzz_xyyyy, tr_y_xxyyzz_xyyyz, tr_y_xxyyzz_xyyzz, tr_y_xxyyzz_xyzzz, tr_y_xxyyzz_xzzzz, tr_y_xxyyzz_yyyyy, tr_y_xxyyzz_yyyyz, tr_y_xxyyzz_yyyzz, tr_y_xxyyzz_yyzzz, tr_y_xxyyzz_yzzzz, tr_y_xxyyzz_zzzzz, tr_y_xyyzz_xxxxz, tr_y_xyyzz_xxxyz, tr_y_xyyzz_xxxz, tr_y_xyyzz_xxxzz, tr_y_xyyzz_xxyyz, tr_y_xyyzz_xxyz, tr_y_xyyzz_xxyzz, tr_y_xyyzz_xxzz, tr_y_xyyzz_xxzzz, tr_y_xyyzz_xyyyz, tr_y_xyyzz_xyyz, tr_y_xyyzz_xyyzz, tr_y_xyyzz_xyzz, tr_y_xyyzz_xyzzz, tr_y_xyyzz_xzzz, tr_y_xyyzz_xzzzz, tr_y_xyyzz_yyyyy, tr_y_xyyzz_yyyyz, tr_y_xyyzz_yyyz, tr_y_xyyzz_yyyzz, tr_y_xyyzz_yyzz, tr_y_xyyzz_yyzzz, tr_y_xyyzz_yzzz, tr_y_xyyzz_yzzzz, tr_y_xyyzz_zzzz, tr_y_xyyzz_zzzzz, tr_y_yyzz_xxxxz, tr_y_yyzz_xxxyz, tr_y_yyzz_xxxzz, tr_y_yyzz_xxyyz, tr_y_yyzz_xxyzz, tr_y_yyzz_xxzzz, tr_y_yyzz_xyyyz, tr_y_yyzz_xyyzz, tr_y_yyzz_xyzzz, tr_y_yyzz_xzzzz, tr_y_yyzz_yyyyy, tr_y_yyzz_yyyyz, tr_y_yyzz_yyyzz, tr_y_yyzz_yyzzz, tr_y_yyzz_yzzzz, tr_y_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyzz_xxxxx[i] = tr_y_xxyy_xxxxx[i] * fe_0 + tr_y_xxyyz_xxxxx[i] * pa_z[i];

        tr_y_xxyyzz_xxxxy[i] = tr_y_xxyy_xxxxy[i] * fe_0 + tr_y_xxyyz_xxxxy[i] * pa_z[i];

        tr_y_xxyyzz_xxxxz[i] = tr_y_yyzz_xxxxz[i] * fe_0 + 4.0 * tr_y_xyyzz_xxxz[i] * fe_0 + tr_y_xyyzz_xxxxz[i] * pa_x[i];

        tr_y_xxyyzz_xxxyy[i] = tr_y_xxyy_xxxyy[i] * fe_0 + tr_y_xxyyz_xxxyy[i] * pa_z[i];

        tr_y_xxyyzz_xxxyz[i] = tr_y_yyzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xyyzz_xxyz[i] * fe_0 + tr_y_xyyzz_xxxyz[i] * pa_x[i];

        tr_y_xxyyzz_xxxzz[i] = tr_y_yyzz_xxxzz[i] * fe_0 + 3.0 * tr_y_xyyzz_xxzz[i] * fe_0 + tr_y_xyyzz_xxxzz[i] * pa_x[i];

        tr_y_xxyyzz_xxyyy[i] = tr_y_xxyy_xxyyy[i] * fe_0 + tr_y_xxyyz_xxyyy[i] * pa_z[i];

        tr_y_xxyyzz_xxyyz[i] = tr_y_yyzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xyyzz_xyyz[i] * fe_0 + tr_y_xyyzz_xxyyz[i] * pa_x[i];

        tr_y_xxyyzz_xxyzz[i] = tr_y_yyzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xyyzz_xyzz[i] * fe_0 + tr_y_xyyzz_xxyzz[i] * pa_x[i];

        tr_y_xxyyzz_xxzzz[i] = tr_y_yyzz_xxzzz[i] * fe_0 + 2.0 * tr_y_xyyzz_xzzz[i] * fe_0 + tr_y_xyyzz_xxzzz[i] * pa_x[i];

        tr_y_xxyyzz_xyyyy[i] = tr_y_xxyy_xyyyy[i] * fe_0 + tr_y_xxyyz_xyyyy[i] * pa_z[i];

        tr_y_xxyyzz_xyyyz[i] = tr_y_yyzz_xyyyz[i] * fe_0 + tr_y_xyyzz_yyyz[i] * fe_0 + tr_y_xyyzz_xyyyz[i] * pa_x[i];

        tr_y_xxyyzz_xyyzz[i] = tr_y_yyzz_xyyzz[i] * fe_0 + tr_y_xyyzz_yyzz[i] * fe_0 + tr_y_xyyzz_xyyzz[i] * pa_x[i];

        tr_y_xxyyzz_xyzzz[i] = tr_y_yyzz_xyzzz[i] * fe_0 + tr_y_xyyzz_yzzz[i] * fe_0 + tr_y_xyyzz_xyzzz[i] * pa_x[i];

        tr_y_xxyyzz_xzzzz[i] = tr_y_yyzz_xzzzz[i] * fe_0 + tr_y_xyyzz_zzzz[i] * fe_0 + tr_y_xyyzz_xzzzz[i] * pa_x[i];

        tr_y_xxyyzz_yyyyy[i] = tr_y_yyzz_yyyyy[i] * fe_0 + tr_y_xyyzz_yyyyy[i] * pa_x[i];

        tr_y_xxyyzz_yyyyz[i] = tr_y_yyzz_yyyyz[i] * fe_0 + tr_y_xyyzz_yyyyz[i] * pa_x[i];

        tr_y_xxyyzz_yyyzz[i] = tr_y_yyzz_yyyzz[i] * fe_0 + tr_y_xyyzz_yyyzz[i] * pa_x[i];

        tr_y_xxyyzz_yyzzz[i] = tr_y_yyzz_yyzzz[i] * fe_0 + tr_y_xyyzz_yyzzz[i] * pa_x[i];

        tr_y_xxyyzz_yzzzz[i] = tr_y_yyzz_yzzzz[i] * fe_0 + tr_y_xyyzz_yzzzz[i] * pa_x[i];

        tr_y_xxyyzz_zzzzz[i] = tr_y_yyzz_zzzzz[i] * fe_0 + tr_y_xyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 861-882 components of targeted buffer : IH

    auto tr_y_xxyzzz_xxxxx = pbuffer.data(idx_dip_ih + 861);

    auto tr_y_xxyzzz_xxxxy = pbuffer.data(idx_dip_ih + 862);

    auto tr_y_xxyzzz_xxxxz = pbuffer.data(idx_dip_ih + 863);

    auto tr_y_xxyzzz_xxxyy = pbuffer.data(idx_dip_ih + 864);

    auto tr_y_xxyzzz_xxxyz = pbuffer.data(idx_dip_ih + 865);

    auto tr_y_xxyzzz_xxxzz = pbuffer.data(idx_dip_ih + 866);

    auto tr_y_xxyzzz_xxyyy = pbuffer.data(idx_dip_ih + 867);

    auto tr_y_xxyzzz_xxyyz = pbuffer.data(idx_dip_ih + 868);

    auto tr_y_xxyzzz_xxyzz = pbuffer.data(idx_dip_ih + 869);

    auto tr_y_xxyzzz_xxzzz = pbuffer.data(idx_dip_ih + 870);

    auto tr_y_xxyzzz_xyyyy = pbuffer.data(idx_dip_ih + 871);

    auto tr_y_xxyzzz_xyyyz = pbuffer.data(idx_dip_ih + 872);

    auto tr_y_xxyzzz_xyyzz = pbuffer.data(idx_dip_ih + 873);

    auto tr_y_xxyzzz_xyzzz = pbuffer.data(idx_dip_ih + 874);

    auto tr_y_xxyzzz_xzzzz = pbuffer.data(idx_dip_ih + 875);

    auto tr_y_xxyzzz_yyyyy = pbuffer.data(idx_dip_ih + 876);

    auto tr_y_xxyzzz_yyyyz = pbuffer.data(idx_dip_ih + 877);

    auto tr_y_xxyzzz_yyyzz = pbuffer.data(idx_dip_ih + 878);

    auto tr_y_xxyzzz_yyzzz = pbuffer.data(idx_dip_ih + 879);

    auto tr_y_xxyzzz_yzzzz = pbuffer.data(idx_dip_ih + 880);

    auto tr_y_xxyzzz_zzzzz = pbuffer.data(idx_dip_ih + 881);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxyz_xxxxy, tr_y_xxyz_xxxyy, tr_y_xxyz_xxyyy, tr_y_xxyz_xyyyy, tr_y_xxyzz_xxxxy, tr_y_xxyzz_xxxyy, tr_y_xxyzz_xxyyy, tr_y_xxyzz_xyyyy, tr_y_xxyzzz_xxxxx, tr_y_xxyzzz_xxxxy, tr_y_xxyzzz_xxxxz, tr_y_xxyzzz_xxxyy, tr_y_xxyzzz_xxxyz, tr_y_xxyzzz_xxxzz, tr_y_xxyzzz_xxyyy, tr_y_xxyzzz_xxyyz, tr_y_xxyzzz_xxyzz, tr_y_xxyzzz_xxzzz, tr_y_xxyzzz_xyyyy, tr_y_xxyzzz_xyyyz, tr_y_xxyzzz_xyyzz, tr_y_xxyzzz_xyzzz, tr_y_xxyzzz_xzzzz, tr_y_xxyzzz_yyyyy, tr_y_xxyzzz_yyyyz, tr_y_xxyzzz_yyyzz, tr_y_xxyzzz_yyzzz, tr_y_xxyzzz_yzzzz, tr_y_xxyzzz_zzzzz, tr_y_xxzzz_xxxxx, tr_y_xxzzz_xxxxz, tr_y_xxzzz_xxxzz, tr_y_xxzzz_xxzzz, tr_y_xxzzz_xzzzz, tr_y_xyzzz_xxxyz, tr_y_xyzzz_xxyyz, tr_y_xyzzz_xxyz, tr_y_xyzzz_xxyzz, tr_y_xyzzz_xyyyz, tr_y_xyzzz_xyyz, tr_y_xyzzz_xyyzz, tr_y_xyzzz_xyzz, tr_y_xyzzz_xyzzz, tr_y_xyzzz_yyyyy, tr_y_xyzzz_yyyyz, tr_y_xyzzz_yyyz, tr_y_xyzzz_yyyzz, tr_y_xyzzz_yyzz, tr_y_xyzzz_yyzzz, tr_y_xyzzz_yzzz, tr_y_xyzzz_yzzzz, tr_y_xyzzz_zzzzz, tr_y_yzzz_xxxyz, tr_y_yzzz_xxyyz, tr_y_yzzz_xxyzz, tr_y_yzzz_xyyyz, tr_y_yzzz_xyyzz, tr_y_yzzz_xyzzz, tr_y_yzzz_yyyyy, tr_y_yzzz_yyyyz, tr_y_yzzz_yyyzz, tr_y_yzzz_yyzzz, tr_y_yzzz_yzzzz, tr_y_yzzz_zzzzz, ts_xxzzz_xxxxx, ts_xxzzz_xxxxz, ts_xxzzz_xxxzz, ts_xxzzz_xxzzz, ts_xxzzz_xzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzzz_xxxxx[i] = ts_xxzzz_xxxxx[i] * fe_0 + tr_y_xxzzz_xxxxx[i] * pa_y[i];

        tr_y_xxyzzz_xxxxy[i] = 2.0 * tr_y_xxyz_xxxxy[i] * fe_0 + tr_y_xxyzz_xxxxy[i] * pa_z[i];

        tr_y_xxyzzz_xxxxz[i] = ts_xxzzz_xxxxz[i] * fe_0 + tr_y_xxzzz_xxxxz[i] * pa_y[i];

        tr_y_xxyzzz_xxxyy[i] = 2.0 * tr_y_xxyz_xxxyy[i] * fe_0 + tr_y_xxyzz_xxxyy[i] * pa_z[i];

        tr_y_xxyzzz_xxxyz[i] = tr_y_yzzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xyzzz_xxyz[i] * fe_0 + tr_y_xyzzz_xxxyz[i] * pa_x[i];

        tr_y_xxyzzz_xxxzz[i] = ts_xxzzz_xxxzz[i] * fe_0 + tr_y_xxzzz_xxxzz[i] * pa_y[i];

        tr_y_xxyzzz_xxyyy[i] = 2.0 * tr_y_xxyz_xxyyy[i] * fe_0 + tr_y_xxyzz_xxyyy[i] * pa_z[i];

        tr_y_xxyzzz_xxyyz[i] = tr_y_yzzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xyzzz_xyyz[i] * fe_0 + tr_y_xyzzz_xxyyz[i] * pa_x[i];

        tr_y_xxyzzz_xxyzz[i] = tr_y_yzzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xyzzz_xyzz[i] * fe_0 + tr_y_xyzzz_xxyzz[i] * pa_x[i];

        tr_y_xxyzzz_xxzzz[i] = ts_xxzzz_xxzzz[i] * fe_0 + tr_y_xxzzz_xxzzz[i] * pa_y[i];

        tr_y_xxyzzz_xyyyy[i] = 2.0 * tr_y_xxyz_xyyyy[i] * fe_0 + tr_y_xxyzz_xyyyy[i] * pa_z[i];

        tr_y_xxyzzz_xyyyz[i] = tr_y_yzzz_xyyyz[i] * fe_0 + tr_y_xyzzz_yyyz[i] * fe_0 + tr_y_xyzzz_xyyyz[i] * pa_x[i];

        tr_y_xxyzzz_xyyzz[i] = tr_y_yzzz_xyyzz[i] * fe_0 + tr_y_xyzzz_yyzz[i] * fe_0 + tr_y_xyzzz_xyyzz[i] * pa_x[i];

        tr_y_xxyzzz_xyzzz[i] = tr_y_yzzz_xyzzz[i] * fe_0 + tr_y_xyzzz_yzzz[i] * fe_0 + tr_y_xyzzz_xyzzz[i] * pa_x[i];

        tr_y_xxyzzz_xzzzz[i] = ts_xxzzz_xzzzz[i] * fe_0 + tr_y_xxzzz_xzzzz[i] * pa_y[i];

        tr_y_xxyzzz_yyyyy[i] = tr_y_yzzz_yyyyy[i] * fe_0 + tr_y_xyzzz_yyyyy[i] * pa_x[i];

        tr_y_xxyzzz_yyyyz[i] = tr_y_yzzz_yyyyz[i] * fe_0 + tr_y_xyzzz_yyyyz[i] * pa_x[i];

        tr_y_xxyzzz_yyyzz[i] = tr_y_yzzz_yyyzz[i] * fe_0 + tr_y_xyzzz_yyyzz[i] * pa_x[i];

        tr_y_xxyzzz_yyzzz[i] = tr_y_yzzz_yyzzz[i] * fe_0 + tr_y_xyzzz_yyzzz[i] * pa_x[i];

        tr_y_xxyzzz_yzzzz[i] = tr_y_yzzz_yzzzz[i] * fe_0 + tr_y_xyzzz_yzzzz[i] * pa_x[i];

        tr_y_xxyzzz_zzzzz[i] = tr_y_yzzz_zzzzz[i] * fe_0 + tr_y_xyzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 882-903 components of targeted buffer : IH

    auto tr_y_xxzzzz_xxxxx = pbuffer.data(idx_dip_ih + 882);

    auto tr_y_xxzzzz_xxxxy = pbuffer.data(idx_dip_ih + 883);

    auto tr_y_xxzzzz_xxxxz = pbuffer.data(idx_dip_ih + 884);

    auto tr_y_xxzzzz_xxxyy = pbuffer.data(idx_dip_ih + 885);

    auto tr_y_xxzzzz_xxxyz = pbuffer.data(idx_dip_ih + 886);

    auto tr_y_xxzzzz_xxxzz = pbuffer.data(idx_dip_ih + 887);

    auto tr_y_xxzzzz_xxyyy = pbuffer.data(idx_dip_ih + 888);

    auto tr_y_xxzzzz_xxyyz = pbuffer.data(idx_dip_ih + 889);

    auto tr_y_xxzzzz_xxyzz = pbuffer.data(idx_dip_ih + 890);

    auto tr_y_xxzzzz_xxzzz = pbuffer.data(idx_dip_ih + 891);

    auto tr_y_xxzzzz_xyyyy = pbuffer.data(idx_dip_ih + 892);

    auto tr_y_xxzzzz_xyyyz = pbuffer.data(idx_dip_ih + 893);

    auto tr_y_xxzzzz_xyyzz = pbuffer.data(idx_dip_ih + 894);

    auto tr_y_xxzzzz_xyzzz = pbuffer.data(idx_dip_ih + 895);

    auto tr_y_xxzzzz_xzzzz = pbuffer.data(idx_dip_ih + 896);

    auto tr_y_xxzzzz_yyyyy = pbuffer.data(idx_dip_ih + 897);

    auto tr_y_xxzzzz_yyyyz = pbuffer.data(idx_dip_ih + 898);

    auto tr_y_xxzzzz_yyyzz = pbuffer.data(idx_dip_ih + 899);

    auto tr_y_xxzzzz_yyzzz = pbuffer.data(idx_dip_ih + 900);

    auto tr_y_xxzzzz_yzzzz = pbuffer.data(idx_dip_ih + 901);

    auto tr_y_xxzzzz_zzzzz = pbuffer.data(idx_dip_ih + 902);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxzz_xxxxx, tr_y_xxzz_xxxxy, tr_y_xxzz_xxxyy, tr_y_xxzz_xxyyy, tr_y_xxzz_xyyyy, tr_y_xxzzz_xxxxx, tr_y_xxzzz_xxxxy, tr_y_xxzzz_xxxyy, tr_y_xxzzz_xxyyy, tr_y_xxzzz_xyyyy, tr_y_xxzzzz_xxxxx, tr_y_xxzzzz_xxxxy, tr_y_xxzzzz_xxxxz, tr_y_xxzzzz_xxxyy, tr_y_xxzzzz_xxxyz, tr_y_xxzzzz_xxxzz, tr_y_xxzzzz_xxyyy, tr_y_xxzzzz_xxyyz, tr_y_xxzzzz_xxyzz, tr_y_xxzzzz_xxzzz, tr_y_xxzzzz_xyyyy, tr_y_xxzzzz_xyyyz, tr_y_xxzzzz_xyyzz, tr_y_xxzzzz_xyzzz, tr_y_xxzzzz_xzzzz, tr_y_xxzzzz_yyyyy, tr_y_xxzzzz_yyyyz, tr_y_xxzzzz_yyyzz, tr_y_xxzzzz_yyzzz, tr_y_xxzzzz_yzzzz, tr_y_xxzzzz_zzzzz, tr_y_xzzzz_xxxxz, tr_y_xzzzz_xxxyz, tr_y_xzzzz_xxxz, tr_y_xzzzz_xxxzz, tr_y_xzzzz_xxyyz, tr_y_xzzzz_xxyz, tr_y_xzzzz_xxyzz, tr_y_xzzzz_xxzz, tr_y_xzzzz_xxzzz, tr_y_xzzzz_xyyyz, tr_y_xzzzz_xyyz, tr_y_xzzzz_xyyzz, tr_y_xzzzz_xyzz, tr_y_xzzzz_xyzzz, tr_y_xzzzz_xzzz, tr_y_xzzzz_xzzzz, tr_y_xzzzz_yyyyy, tr_y_xzzzz_yyyyz, tr_y_xzzzz_yyyz, tr_y_xzzzz_yyyzz, tr_y_xzzzz_yyzz, tr_y_xzzzz_yyzzz, tr_y_xzzzz_yzzz, tr_y_xzzzz_yzzzz, tr_y_xzzzz_zzzz, tr_y_xzzzz_zzzzz, tr_y_zzzz_xxxxz, tr_y_zzzz_xxxyz, tr_y_zzzz_xxxzz, tr_y_zzzz_xxyyz, tr_y_zzzz_xxyzz, tr_y_zzzz_xxzzz, tr_y_zzzz_xyyyz, tr_y_zzzz_xyyzz, tr_y_zzzz_xyzzz, tr_y_zzzz_xzzzz, tr_y_zzzz_yyyyy, tr_y_zzzz_yyyyz, tr_y_zzzz_yyyzz, tr_y_zzzz_yyzzz, tr_y_zzzz_yzzzz, tr_y_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzzz_xxxxx[i] = 3.0 * tr_y_xxzz_xxxxx[i] * fe_0 + tr_y_xxzzz_xxxxx[i] * pa_z[i];

        tr_y_xxzzzz_xxxxy[i] = 3.0 * tr_y_xxzz_xxxxy[i] * fe_0 + tr_y_xxzzz_xxxxy[i] * pa_z[i];

        tr_y_xxzzzz_xxxxz[i] = tr_y_zzzz_xxxxz[i] * fe_0 + 4.0 * tr_y_xzzzz_xxxz[i] * fe_0 + tr_y_xzzzz_xxxxz[i] * pa_x[i];

        tr_y_xxzzzz_xxxyy[i] = 3.0 * tr_y_xxzz_xxxyy[i] * fe_0 + tr_y_xxzzz_xxxyy[i] * pa_z[i];

        tr_y_xxzzzz_xxxyz[i] = tr_y_zzzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xzzzz_xxyz[i] * fe_0 + tr_y_xzzzz_xxxyz[i] * pa_x[i];

        tr_y_xxzzzz_xxxzz[i] = tr_y_zzzz_xxxzz[i] * fe_0 + 3.0 * tr_y_xzzzz_xxzz[i] * fe_0 + tr_y_xzzzz_xxxzz[i] * pa_x[i];

        tr_y_xxzzzz_xxyyy[i] = 3.0 * tr_y_xxzz_xxyyy[i] * fe_0 + tr_y_xxzzz_xxyyy[i] * pa_z[i];

        tr_y_xxzzzz_xxyyz[i] = tr_y_zzzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xzzzz_xyyz[i] * fe_0 + tr_y_xzzzz_xxyyz[i] * pa_x[i];

        tr_y_xxzzzz_xxyzz[i] = tr_y_zzzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xzzzz_xyzz[i] * fe_0 + tr_y_xzzzz_xxyzz[i] * pa_x[i];

        tr_y_xxzzzz_xxzzz[i] = tr_y_zzzz_xxzzz[i] * fe_0 + 2.0 * tr_y_xzzzz_xzzz[i] * fe_0 + tr_y_xzzzz_xxzzz[i] * pa_x[i];

        tr_y_xxzzzz_xyyyy[i] = 3.0 * tr_y_xxzz_xyyyy[i] * fe_0 + tr_y_xxzzz_xyyyy[i] * pa_z[i];

        tr_y_xxzzzz_xyyyz[i] = tr_y_zzzz_xyyyz[i] * fe_0 + tr_y_xzzzz_yyyz[i] * fe_0 + tr_y_xzzzz_xyyyz[i] * pa_x[i];

        tr_y_xxzzzz_xyyzz[i] = tr_y_zzzz_xyyzz[i] * fe_0 + tr_y_xzzzz_yyzz[i] * fe_0 + tr_y_xzzzz_xyyzz[i] * pa_x[i];

        tr_y_xxzzzz_xyzzz[i] = tr_y_zzzz_xyzzz[i] * fe_0 + tr_y_xzzzz_yzzz[i] * fe_0 + tr_y_xzzzz_xyzzz[i] * pa_x[i];

        tr_y_xxzzzz_xzzzz[i] = tr_y_zzzz_xzzzz[i] * fe_0 + tr_y_xzzzz_zzzz[i] * fe_0 + tr_y_xzzzz_xzzzz[i] * pa_x[i];

        tr_y_xxzzzz_yyyyy[i] = tr_y_zzzz_yyyyy[i] * fe_0 + tr_y_xzzzz_yyyyy[i] * pa_x[i];

        tr_y_xxzzzz_yyyyz[i] = tr_y_zzzz_yyyyz[i] * fe_0 + tr_y_xzzzz_yyyyz[i] * pa_x[i];

        tr_y_xxzzzz_yyyzz[i] = tr_y_zzzz_yyyzz[i] * fe_0 + tr_y_xzzzz_yyyzz[i] * pa_x[i];

        tr_y_xxzzzz_yyzzz[i] = tr_y_zzzz_yyzzz[i] * fe_0 + tr_y_xzzzz_yyzzz[i] * pa_x[i];

        tr_y_xxzzzz_yzzzz[i] = tr_y_zzzz_yzzzz[i] * fe_0 + tr_y_xzzzz_yzzzz[i] * pa_x[i];

        tr_y_xxzzzz_zzzzz[i] = tr_y_zzzz_zzzzz[i] * fe_0 + tr_y_xzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 903-924 components of targeted buffer : IH

    auto tr_y_xyyyyy_xxxxx = pbuffer.data(idx_dip_ih + 903);

    auto tr_y_xyyyyy_xxxxy = pbuffer.data(idx_dip_ih + 904);

    auto tr_y_xyyyyy_xxxxz = pbuffer.data(idx_dip_ih + 905);

    auto tr_y_xyyyyy_xxxyy = pbuffer.data(idx_dip_ih + 906);

    auto tr_y_xyyyyy_xxxyz = pbuffer.data(idx_dip_ih + 907);

    auto tr_y_xyyyyy_xxxzz = pbuffer.data(idx_dip_ih + 908);

    auto tr_y_xyyyyy_xxyyy = pbuffer.data(idx_dip_ih + 909);

    auto tr_y_xyyyyy_xxyyz = pbuffer.data(idx_dip_ih + 910);

    auto tr_y_xyyyyy_xxyzz = pbuffer.data(idx_dip_ih + 911);

    auto tr_y_xyyyyy_xxzzz = pbuffer.data(idx_dip_ih + 912);

    auto tr_y_xyyyyy_xyyyy = pbuffer.data(idx_dip_ih + 913);

    auto tr_y_xyyyyy_xyyyz = pbuffer.data(idx_dip_ih + 914);

    auto tr_y_xyyyyy_xyyzz = pbuffer.data(idx_dip_ih + 915);

    auto tr_y_xyyyyy_xyzzz = pbuffer.data(idx_dip_ih + 916);

    auto tr_y_xyyyyy_xzzzz = pbuffer.data(idx_dip_ih + 917);

    auto tr_y_xyyyyy_yyyyy = pbuffer.data(idx_dip_ih + 918);

    auto tr_y_xyyyyy_yyyyz = pbuffer.data(idx_dip_ih + 919);

    auto tr_y_xyyyyy_yyyzz = pbuffer.data(idx_dip_ih + 920);

    auto tr_y_xyyyyy_yyzzz = pbuffer.data(idx_dip_ih + 921);

    auto tr_y_xyyyyy_yzzzz = pbuffer.data(idx_dip_ih + 922);

    auto tr_y_xyyyyy_zzzzz = pbuffer.data(idx_dip_ih + 923);

    #pragma omp simd aligned(pa_x, tr_y_xyyyyy_xxxxx, tr_y_xyyyyy_xxxxy, tr_y_xyyyyy_xxxxz, tr_y_xyyyyy_xxxyy, tr_y_xyyyyy_xxxyz, tr_y_xyyyyy_xxxzz, tr_y_xyyyyy_xxyyy, tr_y_xyyyyy_xxyyz, tr_y_xyyyyy_xxyzz, tr_y_xyyyyy_xxzzz, tr_y_xyyyyy_xyyyy, tr_y_xyyyyy_xyyyz, tr_y_xyyyyy_xyyzz, tr_y_xyyyyy_xyzzz, tr_y_xyyyyy_xzzzz, tr_y_xyyyyy_yyyyy, tr_y_xyyyyy_yyyyz, tr_y_xyyyyy_yyyzz, tr_y_xyyyyy_yyzzz, tr_y_xyyyyy_yzzzz, tr_y_xyyyyy_zzzzz, tr_y_yyyyy_xxxx, tr_y_yyyyy_xxxxx, tr_y_yyyyy_xxxxy, tr_y_yyyyy_xxxxz, tr_y_yyyyy_xxxy, tr_y_yyyyy_xxxyy, tr_y_yyyyy_xxxyz, tr_y_yyyyy_xxxz, tr_y_yyyyy_xxxzz, tr_y_yyyyy_xxyy, tr_y_yyyyy_xxyyy, tr_y_yyyyy_xxyyz, tr_y_yyyyy_xxyz, tr_y_yyyyy_xxyzz, tr_y_yyyyy_xxzz, tr_y_yyyyy_xxzzz, tr_y_yyyyy_xyyy, tr_y_yyyyy_xyyyy, tr_y_yyyyy_xyyyz, tr_y_yyyyy_xyyz, tr_y_yyyyy_xyyzz, tr_y_yyyyy_xyzz, tr_y_yyyyy_xyzzz, tr_y_yyyyy_xzzz, tr_y_yyyyy_xzzzz, tr_y_yyyyy_yyyy, tr_y_yyyyy_yyyyy, tr_y_yyyyy_yyyyz, tr_y_yyyyy_yyyz, tr_y_yyyyy_yyyzz, tr_y_yyyyy_yyzz, tr_y_yyyyy_yyzzz, tr_y_yyyyy_yzzz, tr_y_yyyyy_yzzzz, tr_y_yyyyy_zzzz, tr_y_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyy_xxxxx[i] = 5.0 * tr_y_yyyyy_xxxx[i] * fe_0 + tr_y_yyyyy_xxxxx[i] * pa_x[i];

        tr_y_xyyyyy_xxxxy[i] = 4.0 * tr_y_yyyyy_xxxy[i] * fe_0 + tr_y_yyyyy_xxxxy[i] * pa_x[i];

        tr_y_xyyyyy_xxxxz[i] = 4.0 * tr_y_yyyyy_xxxz[i] * fe_0 + tr_y_yyyyy_xxxxz[i] * pa_x[i];

        tr_y_xyyyyy_xxxyy[i] = 3.0 * tr_y_yyyyy_xxyy[i] * fe_0 + tr_y_yyyyy_xxxyy[i] * pa_x[i];

        tr_y_xyyyyy_xxxyz[i] = 3.0 * tr_y_yyyyy_xxyz[i] * fe_0 + tr_y_yyyyy_xxxyz[i] * pa_x[i];

        tr_y_xyyyyy_xxxzz[i] = 3.0 * tr_y_yyyyy_xxzz[i] * fe_0 + tr_y_yyyyy_xxxzz[i] * pa_x[i];

        tr_y_xyyyyy_xxyyy[i] = 2.0 * tr_y_yyyyy_xyyy[i] * fe_0 + tr_y_yyyyy_xxyyy[i] * pa_x[i];

        tr_y_xyyyyy_xxyyz[i] = 2.0 * tr_y_yyyyy_xyyz[i] * fe_0 + tr_y_yyyyy_xxyyz[i] * pa_x[i];

        tr_y_xyyyyy_xxyzz[i] = 2.0 * tr_y_yyyyy_xyzz[i] * fe_0 + tr_y_yyyyy_xxyzz[i] * pa_x[i];

        tr_y_xyyyyy_xxzzz[i] = 2.0 * tr_y_yyyyy_xzzz[i] * fe_0 + tr_y_yyyyy_xxzzz[i] * pa_x[i];

        tr_y_xyyyyy_xyyyy[i] = tr_y_yyyyy_yyyy[i] * fe_0 + tr_y_yyyyy_xyyyy[i] * pa_x[i];

        tr_y_xyyyyy_xyyyz[i] = tr_y_yyyyy_yyyz[i] * fe_0 + tr_y_yyyyy_xyyyz[i] * pa_x[i];

        tr_y_xyyyyy_xyyzz[i] = tr_y_yyyyy_yyzz[i] * fe_0 + tr_y_yyyyy_xyyzz[i] * pa_x[i];

        tr_y_xyyyyy_xyzzz[i] = tr_y_yyyyy_yzzz[i] * fe_0 + tr_y_yyyyy_xyzzz[i] * pa_x[i];

        tr_y_xyyyyy_xzzzz[i] = tr_y_yyyyy_zzzz[i] * fe_0 + tr_y_yyyyy_xzzzz[i] * pa_x[i];

        tr_y_xyyyyy_yyyyy[i] = tr_y_yyyyy_yyyyy[i] * pa_x[i];

        tr_y_xyyyyy_yyyyz[i] = tr_y_yyyyy_yyyyz[i] * pa_x[i];

        tr_y_xyyyyy_yyyzz[i] = tr_y_yyyyy_yyyzz[i] * pa_x[i];

        tr_y_xyyyyy_yyzzz[i] = tr_y_yyyyy_yyzzz[i] * pa_x[i];

        tr_y_xyyyyy_yzzzz[i] = tr_y_yyyyy_yzzzz[i] * pa_x[i];

        tr_y_xyyyyy_zzzzz[i] = tr_y_yyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 924-945 components of targeted buffer : IH

    auto tr_y_xyyyyz_xxxxx = pbuffer.data(idx_dip_ih + 924);

    auto tr_y_xyyyyz_xxxxy = pbuffer.data(idx_dip_ih + 925);

    auto tr_y_xyyyyz_xxxxz = pbuffer.data(idx_dip_ih + 926);

    auto tr_y_xyyyyz_xxxyy = pbuffer.data(idx_dip_ih + 927);

    auto tr_y_xyyyyz_xxxyz = pbuffer.data(idx_dip_ih + 928);

    auto tr_y_xyyyyz_xxxzz = pbuffer.data(idx_dip_ih + 929);

    auto tr_y_xyyyyz_xxyyy = pbuffer.data(idx_dip_ih + 930);

    auto tr_y_xyyyyz_xxyyz = pbuffer.data(idx_dip_ih + 931);

    auto tr_y_xyyyyz_xxyzz = pbuffer.data(idx_dip_ih + 932);

    auto tr_y_xyyyyz_xxzzz = pbuffer.data(idx_dip_ih + 933);

    auto tr_y_xyyyyz_xyyyy = pbuffer.data(idx_dip_ih + 934);

    auto tr_y_xyyyyz_xyyyz = pbuffer.data(idx_dip_ih + 935);

    auto tr_y_xyyyyz_xyyzz = pbuffer.data(idx_dip_ih + 936);

    auto tr_y_xyyyyz_xyzzz = pbuffer.data(idx_dip_ih + 937);

    auto tr_y_xyyyyz_xzzzz = pbuffer.data(idx_dip_ih + 938);

    auto tr_y_xyyyyz_yyyyy = pbuffer.data(idx_dip_ih + 939);

    auto tr_y_xyyyyz_yyyyz = pbuffer.data(idx_dip_ih + 940);

    auto tr_y_xyyyyz_yyyzz = pbuffer.data(idx_dip_ih + 941);

    auto tr_y_xyyyyz_yyzzz = pbuffer.data(idx_dip_ih + 942);

    auto tr_y_xyyyyz_yzzzz = pbuffer.data(idx_dip_ih + 943);

    auto tr_y_xyyyyz_zzzzz = pbuffer.data(idx_dip_ih + 944);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xyyyy_xxxxx, tr_y_xyyyy_xxxxy, tr_y_xyyyy_xxxyy, tr_y_xyyyy_xxyyy, tr_y_xyyyy_xyyyy, tr_y_xyyyyz_xxxxx, tr_y_xyyyyz_xxxxy, tr_y_xyyyyz_xxxxz, tr_y_xyyyyz_xxxyy, tr_y_xyyyyz_xxxyz, tr_y_xyyyyz_xxxzz, tr_y_xyyyyz_xxyyy, tr_y_xyyyyz_xxyyz, tr_y_xyyyyz_xxyzz, tr_y_xyyyyz_xxzzz, tr_y_xyyyyz_xyyyy, tr_y_xyyyyz_xyyyz, tr_y_xyyyyz_xyyzz, tr_y_xyyyyz_xyzzz, tr_y_xyyyyz_xzzzz, tr_y_xyyyyz_yyyyy, tr_y_xyyyyz_yyyyz, tr_y_xyyyyz_yyyzz, tr_y_xyyyyz_yyzzz, tr_y_xyyyyz_yzzzz, tr_y_xyyyyz_zzzzz, tr_y_yyyyz_xxxxz, tr_y_yyyyz_xxxyz, tr_y_yyyyz_xxxz, tr_y_yyyyz_xxxzz, tr_y_yyyyz_xxyyz, tr_y_yyyyz_xxyz, tr_y_yyyyz_xxyzz, tr_y_yyyyz_xxzz, tr_y_yyyyz_xxzzz, tr_y_yyyyz_xyyyz, tr_y_yyyyz_xyyz, tr_y_yyyyz_xyyzz, tr_y_yyyyz_xyzz, tr_y_yyyyz_xyzzz, tr_y_yyyyz_xzzz, tr_y_yyyyz_xzzzz, tr_y_yyyyz_yyyyy, tr_y_yyyyz_yyyyz, tr_y_yyyyz_yyyz, tr_y_yyyyz_yyyzz, tr_y_yyyyz_yyzz, tr_y_yyyyz_yyzzz, tr_y_yyyyz_yzzz, tr_y_yyyyz_yzzzz, tr_y_yyyyz_zzzz, tr_y_yyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyz_xxxxx[i] = tr_y_xyyyy_xxxxx[i] * pa_z[i];

        tr_y_xyyyyz_xxxxy[i] = tr_y_xyyyy_xxxxy[i] * pa_z[i];

        tr_y_xyyyyz_xxxxz[i] = 4.0 * tr_y_yyyyz_xxxz[i] * fe_0 + tr_y_yyyyz_xxxxz[i] * pa_x[i];

        tr_y_xyyyyz_xxxyy[i] = tr_y_xyyyy_xxxyy[i] * pa_z[i];

        tr_y_xyyyyz_xxxyz[i] = 3.0 * tr_y_yyyyz_xxyz[i] * fe_0 + tr_y_yyyyz_xxxyz[i] * pa_x[i];

        tr_y_xyyyyz_xxxzz[i] = 3.0 * tr_y_yyyyz_xxzz[i] * fe_0 + tr_y_yyyyz_xxxzz[i] * pa_x[i];

        tr_y_xyyyyz_xxyyy[i] = tr_y_xyyyy_xxyyy[i] * pa_z[i];

        tr_y_xyyyyz_xxyyz[i] = 2.0 * tr_y_yyyyz_xyyz[i] * fe_0 + tr_y_yyyyz_xxyyz[i] * pa_x[i];

        tr_y_xyyyyz_xxyzz[i] = 2.0 * tr_y_yyyyz_xyzz[i] * fe_0 + tr_y_yyyyz_xxyzz[i] * pa_x[i];

        tr_y_xyyyyz_xxzzz[i] = 2.0 * tr_y_yyyyz_xzzz[i] * fe_0 + tr_y_yyyyz_xxzzz[i] * pa_x[i];

        tr_y_xyyyyz_xyyyy[i] = tr_y_xyyyy_xyyyy[i] * pa_z[i];

        tr_y_xyyyyz_xyyyz[i] = tr_y_yyyyz_yyyz[i] * fe_0 + tr_y_yyyyz_xyyyz[i] * pa_x[i];

        tr_y_xyyyyz_xyyzz[i] = tr_y_yyyyz_yyzz[i] * fe_0 + tr_y_yyyyz_xyyzz[i] * pa_x[i];

        tr_y_xyyyyz_xyzzz[i] = tr_y_yyyyz_yzzz[i] * fe_0 + tr_y_yyyyz_xyzzz[i] * pa_x[i];

        tr_y_xyyyyz_xzzzz[i] = tr_y_yyyyz_zzzz[i] * fe_0 + tr_y_yyyyz_xzzzz[i] * pa_x[i];

        tr_y_xyyyyz_yyyyy[i] = tr_y_yyyyz_yyyyy[i] * pa_x[i];

        tr_y_xyyyyz_yyyyz[i] = tr_y_yyyyz_yyyyz[i] * pa_x[i];

        tr_y_xyyyyz_yyyzz[i] = tr_y_yyyyz_yyyzz[i] * pa_x[i];

        tr_y_xyyyyz_yyzzz[i] = tr_y_yyyyz_yyzzz[i] * pa_x[i];

        tr_y_xyyyyz_yzzzz[i] = tr_y_yyyyz_yzzzz[i] * pa_x[i];

        tr_y_xyyyyz_zzzzz[i] = tr_y_yyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 945-966 components of targeted buffer : IH

    auto tr_y_xyyyzz_xxxxx = pbuffer.data(idx_dip_ih + 945);

    auto tr_y_xyyyzz_xxxxy = pbuffer.data(idx_dip_ih + 946);

    auto tr_y_xyyyzz_xxxxz = pbuffer.data(idx_dip_ih + 947);

    auto tr_y_xyyyzz_xxxyy = pbuffer.data(idx_dip_ih + 948);

    auto tr_y_xyyyzz_xxxyz = pbuffer.data(idx_dip_ih + 949);

    auto tr_y_xyyyzz_xxxzz = pbuffer.data(idx_dip_ih + 950);

    auto tr_y_xyyyzz_xxyyy = pbuffer.data(idx_dip_ih + 951);

    auto tr_y_xyyyzz_xxyyz = pbuffer.data(idx_dip_ih + 952);

    auto tr_y_xyyyzz_xxyzz = pbuffer.data(idx_dip_ih + 953);

    auto tr_y_xyyyzz_xxzzz = pbuffer.data(idx_dip_ih + 954);

    auto tr_y_xyyyzz_xyyyy = pbuffer.data(idx_dip_ih + 955);

    auto tr_y_xyyyzz_xyyyz = pbuffer.data(idx_dip_ih + 956);

    auto tr_y_xyyyzz_xyyzz = pbuffer.data(idx_dip_ih + 957);

    auto tr_y_xyyyzz_xyzzz = pbuffer.data(idx_dip_ih + 958);

    auto tr_y_xyyyzz_xzzzz = pbuffer.data(idx_dip_ih + 959);

    auto tr_y_xyyyzz_yyyyy = pbuffer.data(idx_dip_ih + 960);

    auto tr_y_xyyyzz_yyyyz = pbuffer.data(idx_dip_ih + 961);

    auto tr_y_xyyyzz_yyyzz = pbuffer.data(idx_dip_ih + 962);

    auto tr_y_xyyyzz_yyzzz = pbuffer.data(idx_dip_ih + 963);

    auto tr_y_xyyyzz_yzzzz = pbuffer.data(idx_dip_ih + 964);

    auto tr_y_xyyyzz_zzzzz = pbuffer.data(idx_dip_ih + 965);

    #pragma omp simd aligned(pa_x, tr_y_xyyyzz_xxxxx, tr_y_xyyyzz_xxxxy, tr_y_xyyyzz_xxxxz, tr_y_xyyyzz_xxxyy, tr_y_xyyyzz_xxxyz, tr_y_xyyyzz_xxxzz, tr_y_xyyyzz_xxyyy, tr_y_xyyyzz_xxyyz, tr_y_xyyyzz_xxyzz, tr_y_xyyyzz_xxzzz, tr_y_xyyyzz_xyyyy, tr_y_xyyyzz_xyyyz, tr_y_xyyyzz_xyyzz, tr_y_xyyyzz_xyzzz, tr_y_xyyyzz_xzzzz, tr_y_xyyyzz_yyyyy, tr_y_xyyyzz_yyyyz, tr_y_xyyyzz_yyyzz, tr_y_xyyyzz_yyzzz, tr_y_xyyyzz_yzzzz, tr_y_xyyyzz_zzzzz, tr_y_yyyzz_xxxx, tr_y_yyyzz_xxxxx, tr_y_yyyzz_xxxxy, tr_y_yyyzz_xxxxz, tr_y_yyyzz_xxxy, tr_y_yyyzz_xxxyy, tr_y_yyyzz_xxxyz, tr_y_yyyzz_xxxz, tr_y_yyyzz_xxxzz, tr_y_yyyzz_xxyy, tr_y_yyyzz_xxyyy, tr_y_yyyzz_xxyyz, tr_y_yyyzz_xxyz, tr_y_yyyzz_xxyzz, tr_y_yyyzz_xxzz, tr_y_yyyzz_xxzzz, tr_y_yyyzz_xyyy, tr_y_yyyzz_xyyyy, tr_y_yyyzz_xyyyz, tr_y_yyyzz_xyyz, tr_y_yyyzz_xyyzz, tr_y_yyyzz_xyzz, tr_y_yyyzz_xyzzz, tr_y_yyyzz_xzzz, tr_y_yyyzz_xzzzz, tr_y_yyyzz_yyyy, tr_y_yyyzz_yyyyy, tr_y_yyyzz_yyyyz, tr_y_yyyzz_yyyz, tr_y_yyyzz_yyyzz, tr_y_yyyzz_yyzz, tr_y_yyyzz_yyzzz, tr_y_yyyzz_yzzz, tr_y_yyyzz_yzzzz, tr_y_yyyzz_zzzz, tr_y_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyzz_xxxxx[i] = 5.0 * tr_y_yyyzz_xxxx[i] * fe_0 + tr_y_yyyzz_xxxxx[i] * pa_x[i];

        tr_y_xyyyzz_xxxxy[i] = 4.0 * tr_y_yyyzz_xxxy[i] * fe_0 + tr_y_yyyzz_xxxxy[i] * pa_x[i];

        tr_y_xyyyzz_xxxxz[i] = 4.0 * tr_y_yyyzz_xxxz[i] * fe_0 + tr_y_yyyzz_xxxxz[i] * pa_x[i];

        tr_y_xyyyzz_xxxyy[i] = 3.0 * tr_y_yyyzz_xxyy[i] * fe_0 + tr_y_yyyzz_xxxyy[i] * pa_x[i];

        tr_y_xyyyzz_xxxyz[i] = 3.0 * tr_y_yyyzz_xxyz[i] * fe_0 + tr_y_yyyzz_xxxyz[i] * pa_x[i];

        tr_y_xyyyzz_xxxzz[i] = 3.0 * tr_y_yyyzz_xxzz[i] * fe_0 + tr_y_yyyzz_xxxzz[i] * pa_x[i];

        tr_y_xyyyzz_xxyyy[i] = 2.0 * tr_y_yyyzz_xyyy[i] * fe_0 + tr_y_yyyzz_xxyyy[i] * pa_x[i];

        tr_y_xyyyzz_xxyyz[i] = 2.0 * tr_y_yyyzz_xyyz[i] * fe_0 + tr_y_yyyzz_xxyyz[i] * pa_x[i];

        tr_y_xyyyzz_xxyzz[i] = 2.0 * tr_y_yyyzz_xyzz[i] * fe_0 + tr_y_yyyzz_xxyzz[i] * pa_x[i];

        tr_y_xyyyzz_xxzzz[i] = 2.0 * tr_y_yyyzz_xzzz[i] * fe_0 + tr_y_yyyzz_xxzzz[i] * pa_x[i];

        tr_y_xyyyzz_xyyyy[i] = tr_y_yyyzz_yyyy[i] * fe_0 + tr_y_yyyzz_xyyyy[i] * pa_x[i];

        tr_y_xyyyzz_xyyyz[i] = tr_y_yyyzz_yyyz[i] * fe_0 + tr_y_yyyzz_xyyyz[i] * pa_x[i];

        tr_y_xyyyzz_xyyzz[i] = tr_y_yyyzz_yyzz[i] * fe_0 + tr_y_yyyzz_xyyzz[i] * pa_x[i];

        tr_y_xyyyzz_xyzzz[i] = tr_y_yyyzz_yzzz[i] * fe_0 + tr_y_yyyzz_xyzzz[i] * pa_x[i];

        tr_y_xyyyzz_xzzzz[i] = tr_y_yyyzz_zzzz[i] * fe_0 + tr_y_yyyzz_xzzzz[i] * pa_x[i];

        tr_y_xyyyzz_yyyyy[i] = tr_y_yyyzz_yyyyy[i] * pa_x[i];

        tr_y_xyyyzz_yyyyz[i] = tr_y_yyyzz_yyyyz[i] * pa_x[i];

        tr_y_xyyyzz_yyyzz[i] = tr_y_yyyzz_yyyzz[i] * pa_x[i];

        tr_y_xyyyzz_yyzzz[i] = tr_y_yyyzz_yyzzz[i] * pa_x[i];

        tr_y_xyyyzz_yzzzz[i] = tr_y_yyyzz_yzzzz[i] * pa_x[i];

        tr_y_xyyyzz_zzzzz[i] = tr_y_yyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 966-987 components of targeted buffer : IH

    auto tr_y_xyyzzz_xxxxx = pbuffer.data(idx_dip_ih + 966);

    auto tr_y_xyyzzz_xxxxy = pbuffer.data(idx_dip_ih + 967);

    auto tr_y_xyyzzz_xxxxz = pbuffer.data(idx_dip_ih + 968);

    auto tr_y_xyyzzz_xxxyy = pbuffer.data(idx_dip_ih + 969);

    auto tr_y_xyyzzz_xxxyz = pbuffer.data(idx_dip_ih + 970);

    auto tr_y_xyyzzz_xxxzz = pbuffer.data(idx_dip_ih + 971);

    auto tr_y_xyyzzz_xxyyy = pbuffer.data(idx_dip_ih + 972);

    auto tr_y_xyyzzz_xxyyz = pbuffer.data(idx_dip_ih + 973);

    auto tr_y_xyyzzz_xxyzz = pbuffer.data(idx_dip_ih + 974);

    auto tr_y_xyyzzz_xxzzz = pbuffer.data(idx_dip_ih + 975);

    auto tr_y_xyyzzz_xyyyy = pbuffer.data(idx_dip_ih + 976);

    auto tr_y_xyyzzz_xyyyz = pbuffer.data(idx_dip_ih + 977);

    auto tr_y_xyyzzz_xyyzz = pbuffer.data(idx_dip_ih + 978);

    auto tr_y_xyyzzz_xyzzz = pbuffer.data(idx_dip_ih + 979);

    auto tr_y_xyyzzz_xzzzz = pbuffer.data(idx_dip_ih + 980);

    auto tr_y_xyyzzz_yyyyy = pbuffer.data(idx_dip_ih + 981);

    auto tr_y_xyyzzz_yyyyz = pbuffer.data(idx_dip_ih + 982);

    auto tr_y_xyyzzz_yyyzz = pbuffer.data(idx_dip_ih + 983);

    auto tr_y_xyyzzz_yyzzz = pbuffer.data(idx_dip_ih + 984);

    auto tr_y_xyyzzz_yzzzz = pbuffer.data(idx_dip_ih + 985);

    auto tr_y_xyyzzz_zzzzz = pbuffer.data(idx_dip_ih + 986);

    #pragma omp simd aligned(pa_x, tr_y_xyyzzz_xxxxx, tr_y_xyyzzz_xxxxy, tr_y_xyyzzz_xxxxz, tr_y_xyyzzz_xxxyy, tr_y_xyyzzz_xxxyz, tr_y_xyyzzz_xxxzz, tr_y_xyyzzz_xxyyy, tr_y_xyyzzz_xxyyz, tr_y_xyyzzz_xxyzz, tr_y_xyyzzz_xxzzz, tr_y_xyyzzz_xyyyy, tr_y_xyyzzz_xyyyz, tr_y_xyyzzz_xyyzz, tr_y_xyyzzz_xyzzz, tr_y_xyyzzz_xzzzz, tr_y_xyyzzz_yyyyy, tr_y_xyyzzz_yyyyz, tr_y_xyyzzz_yyyzz, tr_y_xyyzzz_yyzzz, tr_y_xyyzzz_yzzzz, tr_y_xyyzzz_zzzzz, tr_y_yyzzz_xxxx, tr_y_yyzzz_xxxxx, tr_y_yyzzz_xxxxy, tr_y_yyzzz_xxxxz, tr_y_yyzzz_xxxy, tr_y_yyzzz_xxxyy, tr_y_yyzzz_xxxyz, tr_y_yyzzz_xxxz, tr_y_yyzzz_xxxzz, tr_y_yyzzz_xxyy, tr_y_yyzzz_xxyyy, tr_y_yyzzz_xxyyz, tr_y_yyzzz_xxyz, tr_y_yyzzz_xxyzz, tr_y_yyzzz_xxzz, tr_y_yyzzz_xxzzz, tr_y_yyzzz_xyyy, tr_y_yyzzz_xyyyy, tr_y_yyzzz_xyyyz, tr_y_yyzzz_xyyz, tr_y_yyzzz_xyyzz, tr_y_yyzzz_xyzz, tr_y_yyzzz_xyzzz, tr_y_yyzzz_xzzz, tr_y_yyzzz_xzzzz, tr_y_yyzzz_yyyy, tr_y_yyzzz_yyyyy, tr_y_yyzzz_yyyyz, tr_y_yyzzz_yyyz, tr_y_yyzzz_yyyzz, tr_y_yyzzz_yyzz, tr_y_yyzzz_yyzzz, tr_y_yyzzz_yzzz, tr_y_yyzzz_yzzzz, tr_y_yyzzz_zzzz, tr_y_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzzz_xxxxx[i] = 5.0 * tr_y_yyzzz_xxxx[i] * fe_0 + tr_y_yyzzz_xxxxx[i] * pa_x[i];

        tr_y_xyyzzz_xxxxy[i] = 4.0 * tr_y_yyzzz_xxxy[i] * fe_0 + tr_y_yyzzz_xxxxy[i] * pa_x[i];

        tr_y_xyyzzz_xxxxz[i] = 4.0 * tr_y_yyzzz_xxxz[i] * fe_0 + tr_y_yyzzz_xxxxz[i] * pa_x[i];

        tr_y_xyyzzz_xxxyy[i] = 3.0 * tr_y_yyzzz_xxyy[i] * fe_0 + tr_y_yyzzz_xxxyy[i] * pa_x[i];

        tr_y_xyyzzz_xxxyz[i] = 3.0 * tr_y_yyzzz_xxyz[i] * fe_0 + tr_y_yyzzz_xxxyz[i] * pa_x[i];

        tr_y_xyyzzz_xxxzz[i] = 3.0 * tr_y_yyzzz_xxzz[i] * fe_0 + tr_y_yyzzz_xxxzz[i] * pa_x[i];

        tr_y_xyyzzz_xxyyy[i] = 2.0 * tr_y_yyzzz_xyyy[i] * fe_0 + tr_y_yyzzz_xxyyy[i] * pa_x[i];

        tr_y_xyyzzz_xxyyz[i] = 2.0 * tr_y_yyzzz_xyyz[i] * fe_0 + tr_y_yyzzz_xxyyz[i] * pa_x[i];

        tr_y_xyyzzz_xxyzz[i] = 2.0 * tr_y_yyzzz_xyzz[i] * fe_0 + tr_y_yyzzz_xxyzz[i] * pa_x[i];

        tr_y_xyyzzz_xxzzz[i] = 2.0 * tr_y_yyzzz_xzzz[i] * fe_0 + tr_y_yyzzz_xxzzz[i] * pa_x[i];

        tr_y_xyyzzz_xyyyy[i] = tr_y_yyzzz_yyyy[i] * fe_0 + tr_y_yyzzz_xyyyy[i] * pa_x[i];

        tr_y_xyyzzz_xyyyz[i] = tr_y_yyzzz_yyyz[i] * fe_0 + tr_y_yyzzz_xyyyz[i] * pa_x[i];

        tr_y_xyyzzz_xyyzz[i] = tr_y_yyzzz_yyzz[i] * fe_0 + tr_y_yyzzz_xyyzz[i] * pa_x[i];

        tr_y_xyyzzz_xyzzz[i] = tr_y_yyzzz_yzzz[i] * fe_0 + tr_y_yyzzz_xyzzz[i] * pa_x[i];

        tr_y_xyyzzz_xzzzz[i] = tr_y_yyzzz_zzzz[i] * fe_0 + tr_y_yyzzz_xzzzz[i] * pa_x[i];

        tr_y_xyyzzz_yyyyy[i] = tr_y_yyzzz_yyyyy[i] * pa_x[i];

        tr_y_xyyzzz_yyyyz[i] = tr_y_yyzzz_yyyyz[i] * pa_x[i];

        tr_y_xyyzzz_yyyzz[i] = tr_y_yyzzz_yyyzz[i] * pa_x[i];

        tr_y_xyyzzz_yyzzz[i] = tr_y_yyzzz_yyzzz[i] * pa_x[i];

        tr_y_xyyzzz_yzzzz[i] = tr_y_yyzzz_yzzzz[i] * pa_x[i];

        tr_y_xyyzzz_zzzzz[i] = tr_y_yyzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 987-1008 components of targeted buffer : IH

    auto tr_y_xyzzzz_xxxxx = pbuffer.data(idx_dip_ih + 987);

    auto tr_y_xyzzzz_xxxxy = pbuffer.data(idx_dip_ih + 988);

    auto tr_y_xyzzzz_xxxxz = pbuffer.data(idx_dip_ih + 989);

    auto tr_y_xyzzzz_xxxyy = pbuffer.data(idx_dip_ih + 990);

    auto tr_y_xyzzzz_xxxyz = pbuffer.data(idx_dip_ih + 991);

    auto tr_y_xyzzzz_xxxzz = pbuffer.data(idx_dip_ih + 992);

    auto tr_y_xyzzzz_xxyyy = pbuffer.data(idx_dip_ih + 993);

    auto tr_y_xyzzzz_xxyyz = pbuffer.data(idx_dip_ih + 994);

    auto tr_y_xyzzzz_xxyzz = pbuffer.data(idx_dip_ih + 995);

    auto tr_y_xyzzzz_xxzzz = pbuffer.data(idx_dip_ih + 996);

    auto tr_y_xyzzzz_xyyyy = pbuffer.data(idx_dip_ih + 997);

    auto tr_y_xyzzzz_xyyyz = pbuffer.data(idx_dip_ih + 998);

    auto tr_y_xyzzzz_xyyzz = pbuffer.data(idx_dip_ih + 999);

    auto tr_y_xyzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1000);

    auto tr_y_xyzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1001);

    auto tr_y_xyzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1002);

    auto tr_y_xyzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1003);

    auto tr_y_xyzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1004);

    auto tr_y_xyzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1005);

    auto tr_y_xyzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1006);

    auto tr_y_xyzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1007);

    #pragma omp simd aligned(pa_x, tr_y_xyzzzz_xxxxx, tr_y_xyzzzz_xxxxy, tr_y_xyzzzz_xxxxz, tr_y_xyzzzz_xxxyy, tr_y_xyzzzz_xxxyz, tr_y_xyzzzz_xxxzz, tr_y_xyzzzz_xxyyy, tr_y_xyzzzz_xxyyz, tr_y_xyzzzz_xxyzz, tr_y_xyzzzz_xxzzz, tr_y_xyzzzz_xyyyy, tr_y_xyzzzz_xyyyz, tr_y_xyzzzz_xyyzz, tr_y_xyzzzz_xyzzz, tr_y_xyzzzz_xzzzz, tr_y_xyzzzz_yyyyy, tr_y_xyzzzz_yyyyz, tr_y_xyzzzz_yyyzz, tr_y_xyzzzz_yyzzz, tr_y_xyzzzz_yzzzz, tr_y_xyzzzz_zzzzz, tr_y_yzzzz_xxxx, tr_y_yzzzz_xxxxx, tr_y_yzzzz_xxxxy, tr_y_yzzzz_xxxxz, tr_y_yzzzz_xxxy, tr_y_yzzzz_xxxyy, tr_y_yzzzz_xxxyz, tr_y_yzzzz_xxxz, tr_y_yzzzz_xxxzz, tr_y_yzzzz_xxyy, tr_y_yzzzz_xxyyy, tr_y_yzzzz_xxyyz, tr_y_yzzzz_xxyz, tr_y_yzzzz_xxyzz, tr_y_yzzzz_xxzz, tr_y_yzzzz_xxzzz, tr_y_yzzzz_xyyy, tr_y_yzzzz_xyyyy, tr_y_yzzzz_xyyyz, tr_y_yzzzz_xyyz, tr_y_yzzzz_xyyzz, tr_y_yzzzz_xyzz, tr_y_yzzzz_xyzzz, tr_y_yzzzz_xzzz, tr_y_yzzzz_xzzzz, tr_y_yzzzz_yyyy, tr_y_yzzzz_yyyyy, tr_y_yzzzz_yyyyz, tr_y_yzzzz_yyyz, tr_y_yzzzz_yyyzz, tr_y_yzzzz_yyzz, tr_y_yzzzz_yyzzz, tr_y_yzzzz_yzzz, tr_y_yzzzz_yzzzz, tr_y_yzzzz_zzzz, tr_y_yzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzzz_xxxxx[i] = 5.0 * tr_y_yzzzz_xxxx[i] * fe_0 + tr_y_yzzzz_xxxxx[i] * pa_x[i];

        tr_y_xyzzzz_xxxxy[i] = 4.0 * tr_y_yzzzz_xxxy[i] * fe_0 + tr_y_yzzzz_xxxxy[i] * pa_x[i];

        tr_y_xyzzzz_xxxxz[i] = 4.0 * tr_y_yzzzz_xxxz[i] * fe_0 + tr_y_yzzzz_xxxxz[i] * pa_x[i];

        tr_y_xyzzzz_xxxyy[i] = 3.0 * tr_y_yzzzz_xxyy[i] * fe_0 + tr_y_yzzzz_xxxyy[i] * pa_x[i];

        tr_y_xyzzzz_xxxyz[i] = 3.0 * tr_y_yzzzz_xxyz[i] * fe_0 + tr_y_yzzzz_xxxyz[i] * pa_x[i];

        tr_y_xyzzzz_xxxzz[i] = 3.0 * tr_y_yzzzz_xxzz[i] * fe_0 + tr_y_yzzzz_xxxzz[i] * pa_x[i];

        tr_y_xyzzzz_xxyyy[i] = 2.0 * tr_y_yzzzz_xyyy[i] * fe_0 + tr_y_yzzzz_xxyyy[i] * pa_x[i];

        tr_y_xyzzzz_xxyyz[i] = 2.0 * tr_y_yzzzz_xyyz[i] * fe_0 + tr_y_yzzzz_xxyyz[i] * pa_x[i];

        tr_y_xyzzzz_xxyzz[i] = 2.0 * tr_y_yzzzz_xyzz[i] * fe_0 + tr_y_yzzzz_xxyzz[i] * pa_x[i];

        tr_y_xyzzzz_xxzzz[i] = 2.0 * tr_y_yzzzz_xzzz[i] * fe_0 + tr_y_yzzzz_xxzzz[i] * pa_x[i];

        tr_y_xyzzzz_xyyyy[i] = tr_y_yzzzz_yyyy[i] * fe_0 + tr_y_yzzzz_xyyyy[i] * pa_x[i];

        tr_y_xyzzzz_xyyyz[i] = tr_y_yzzzz_yyyz[i] * fe_0 + tr_y_yzzzz_xyyyz[i] * pa_x[i];

        tr_y_xyzzzz_xyyzz[i] = tr_y_yzzzz_yyzz[i] * fe_0 + tr_y_yzzzz_xyyzz[i] * pa_x[i];

        tr_y_xyzzzz_xyzzz[i] = tr_y_yzzzz_yzzz[i] * fe_0 + tr_y_yzzzz_xyzzz[i] * pa_x[i];

        tr_y_xyzzzz_xzzzz[i] = tr_y_yzzzz_zzzz[i] * fe_0 + tr_y_yzzzz_xzzzz[i] * pa_x[i];

        tr_y_xyzzzz_yyyyy[i] = tr_y_yzzzz_yyyyy[i] * pa_x[i];

        tr_y_xyzzzz_yyyyz[i] = tr_y_yzzzz_yyyyz[i] * pa_x[i];

        tr_y_xyzzzz_yyyzz[i] = tr_y_yzzzz_yyyzz[i] * pa_x[i];

        tr_y_xyzzzz_yyzzz[i] = tr_y_yzzzz_yyzzz[i] * pa_x[i];

        tr_y_xyzzzz_yzzzz[i] = tr_y_yzzzz_yzzzz[i] * pa_x[i];

        tr_y_xyzzzz_zzzzz[i] = tr_y_yzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1008-1029 components of targeted buffer : IH

    auto tr_y_xzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1008);

    auto tr_y_xzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1009);

    auto tr_y_xzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1010);

    auto tr_y_xzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1011);

    auto tr_y_xzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1012);

    auto tr_y_xzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1013);

    auto tr_y_xzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1014);

    auto tr_y_xzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1015);

    auto tr_y_xzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1016);

    auto tr_y_xzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1017);

    auto tr_y_xzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1018);

    auto tr_y_xzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1019);

    auto tr_y_xzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1020);

    auto tr_y_xzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1021);

    auto tr_y_xzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1022);

    auto tr_y_xzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1023);

    auto tr_y_xzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1024);

    auto tr_y_xzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1025);

    auto tr_y_xzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1026);

    auto tr_y_xzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1027);

    auto tr_y_xzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1028);

    #pragma omp simd aligned(pa_x, tr_y_xzzzzz_xxxxx, tr_y_xzzzzz_xxxxy, tr_y_xzzzzz_xxxxz, tr_y_xzzzzz_xxxyy, tr_y_xzzzzz_xxxyz, tr_y_xzzzzz_xxxzz, tr_y_xzzzzz_xxyyy, tr_y_xzzzzz_xxyyz, tr_y_xzzzzz_xxyzz, tr_y_xzzzzz_xxzzz, tr_y_xzzzzz_xyyyy, tr_y_xzzzzz_xyyyz, tr_y_xzzzzz_xyyzz, tr_y_xzzzzz_xyzzz, tr_y_xzzzzz_xzzzz, tr_y_xzzzzz_yyyyy, tr_y_xzzzzz_yyyyz, tr_y_xzzzzz_yyyzz, tr_y_xzzzzz_yyzzz, tr_y_xzzzzz_yzzzz, tr_y_xzzzzz_zzzzz, tr_y_zzzzz_xxxx, tr_y_zzzzz_xxxxx, tr_y_zzzzz_xxxxy, tr_y_zzzzz_xxxxz, tr_y_zzzzz_xxxy, tr_y_zzzzz_xxxyy, tr_y_zzzzz_xxxyz, tr_y_zzzzz_xxxz, tr_y_zzzzz_xxxzz, tr_y_zzzzz_xxyy, tr_y_zzzzz_xxyyy, tr_y_zzzzz_xxyyz, tr_y_zzzzz_xxyz, tr_y_zzzzz_xxyzz, tr_y_zzzzz_xxzz, tr_y_zzzzz_xxzzz, tr_y_zzzzz_xyyy, tr_y_zzzzz_xyyyy, tr_y_zzzzz_xyyyz, tr_y_zzzzz_xyyz, tr_y_zzzzz_xyyzz, tr_y_zzzzz_xyzz, tr_y_zzzzz_xyzzz, tr_y_zzzzz_xzzz, tr_y_zzzzz_xzzzz, tr_y_zzzzz_yyyy, tr_y_zzzzz_yyyyy, tr_y_zzzzz_yyyyz, tr_y_zzzzz_yyyz, tr_y_zzzzz_yyyzz, tr_y_zzzzz_yyzz, tr_y_zzzzz_yyzzz, tr_y_zzzzz_yzzz, tr_y_zzzzz_yzzzz, tr_y_zzzzz_zzzz, tr_y_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzzz_xxxxx[i] = 5.0 * tr_y_zzzzz_xxxx[i] * fe_0 + tr_y_zzzzz_xxxxx[i] * pa_x[i];

        tr_y_xzzzzz_xxxxy[i] = 4.0 * tr_y_zzzzz_xxxy[i] * fe_0 + tr_y_zzzzz_xxxxy[i] * pa_x[i];

        tr_y_xzzzzz_xxxxz[i] = 4.0 * tr_y_zzzzz_xxxz[i] * fe_0 + tr_y_zzzzz_xxxxz[i] * pa_x[i];

        tr_y_xzzzzz_xxxyy[i] = 3.0 * tr_y_zzzzz_xxyy[i] * fe_0 + tr_y_zzzzz_xxxyy[i] * pa_x[i];

        tr_y_xzzzzz_xxxyz[i] = 3.0 * tr_y_zzzzz_xxyz[i] * fe_0 + tr_y_zzzzz_xxxyz[i] * pa_x[i];

        tr_y_xzzzzz_xxxzz[i] = 3.0 * tr_y_zzzzz_xxzz[i] * fe_0 + tr_y_zzzzz_xxxzz[i] * pa_x[i];

        tr_y_xzzzzz_xxyyy[i] = 2.0 * tr_y_zzzzz_xyyy[i] * fe_0 + tr_y_zzzzz_xxyyy[i] * pa_x[i];

        tr_y_xzzzzz_xxyyz[i] = 2.0 * tr_y_zzzzz_xyyz[i] * fe_0 + tr_y_zzzzz_xxyyz[i] * pa_x[i];

        tr_y_xzzzzz_xxyzz[i] = 2.0 * tr_y_zzzzz_xyzz[i] * fe_0 + tr_y_zzzzz_xxyzz[i] * pa_x[i];

        tr_y_xzzzzz_xxzzz[i] = 2.0 * tr_y_zzzzz_xzzz[i] * fe_0 + tr_y_zzzzz_xxzzz[i] * pa_x[i];

        tr_y_xzzzzz_xyyyy[i] = tr_y_zzzzz_yyyy[i] * fe_0 + tr_y_zzzzz_xyyyy[i] * pa_x[i];

        tr_y_xzzzzz_xyyyz[i] = tr_y_zzzzz_yyyz[i] * fe_0 + tr_y_zzzzz_xyyyz[i] * pa_x[i];

        tr_y_xzzzzz_xyyzz[i] = tr_y_zzzzz_yyzz[i] * fe_0 + tr_y_zzzzz_xyyzz[i] * pa_x[i];

        tr_y_xzzzzz_xyzzz[i] = tr_y_zzzzz_yzzz[i] * fe_0 + tr_y_zzzzz_xyzzz[i] * pa_x[i];

        tr_y_xzzzzz_xzzzz[i] = tr_y_zzzzz_zzzz[i] * fe_0 + tr_y_zzzzz_xzzzz[i] * pa_x[i];

        tr_y_xzzzzz_yyyyy[i] = tr_y_zzzzz_yyyyy[i] * pa_x[i];

        tr_y_xzzzzz_yyyyz[i] = tr_y_zzzzz_yyyyz[i] * pa_x[i];

        tr_y_xzzzzz_yyyzz[i] = tr_y_zzzzz_yyyzz[i] * pa_x[i];

        tr_y_xzzzzz_yyzzz[i] = tr_y_zzzzz_yyzzz[i] * pa_x[i];

        tr_y_xzzzzz_yzzzz[i] = tr_y_zzzzz_yzzzz[i] * pa_x[i];

        tr_y_xzzzzz_zzzzz[i] = tr_y_zzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1029-1050 components of targeted buffer : IH

    auto tr_y_yyyyyy_xxxxx = pbuffer.data(idx_dip_ih + 1029);

    auto tr_y_yyyyyy_xxxxy = pbuffer.data(idx_dip_ih + 1030);

    auto tr_y_yyyyyy_xxxxz = pbuffer.data(idx_dip_ih + 1031);

    auto tr_y_yyyyyy_xxxyy = pbuffer.data(idx_dip_ih + 1032);

    auto tr_y_yyyyyy_xxxyz = pbuffer.data(idx_dip_ih + 1033);

    auto tr_y_yyyyyy_xxxzz = pbuffer.data(idx_dip_ih + 1034);

    auto tr_y_yyyyyy_xxyyy = pbuffer.data(idx_dip_ih + 1035);

    auto tr_y_yyyyyy_xxyyz = pbuffer.data(idx_dip_ih + 1036);

    auto tr_y_yyyyyy_xxyzz = pbuffer.data(idx_dip_ih + 1037);

    auto tr_y_yyyyyy_xxzzz = pbuffer.data(idx_dip_ih + 1038);

    auto tr_y_yyyyyy_xyyyy = pbuffer.data(idx_dip_ih + 1039);

    auto tr_y_yyyyyy_xyyyz = pbuffer.data(idx_dip_ih + 1040);

    auto tr_y_yyyyyy_xyyzz = pbuffer.data(idx_dip_ih + 1041);

    auto tr_y_yyyyyy_xyzzz = pbuffer.data(idx_dip_ih + 1042);

    auto tr_y_yyyyyy_xzzzz = pbuffer.data(idx_dip_ih + 1043);

    auto tr_y_yyyyyy_yyyyy = pbuffer.data(idx_dip_ih + 1044);

    auto tr_y_yyyyyy_yyyyz = pbuffer.data(idx_dip_ih + 1045);

    auto tr_y_yyyyyy_yyyzz = pbuffer.data(idx_dip_ih + 1046);

    auto tr_y_yyyyyy_yyzzz = pbuffer.data(idx_dip_ih + 1047);

    auto tr_y_yyyyyy_yzzzz = pbuffer.data(idx_dip_ih + 1048);

    auto tr_y_yyyyyy_zzzzz = pbuffer.data(idx_dip_ih + 1049);

    #pragma omp simd aligned(pa_y, tr_y_yyyy_xxxxx, tr_y_yyyy_xxxxy, tr_y_yyyy_xxxxz, tr_y_yyyy_xxxyy, tr_y_yyyy_xxxyz, tr_y_yyyy_xxxzz, tr_y_yyyy_xxyyy, tr_y_yyyy_xxyyz, tr_y_yyyy_xxyzz, tr_y_yyyy_xxzzz, tr_y_yyyy_xyyyy, tr_y_yyyy_xyyyz, tr_y_yyyy_xyyzz, tr_y_yyyy_xyzzz, tr_y_yyyy_xzzzz, tr_y_yyyy_yyyyy, tr_y_yyyy_yyyyz, tr_y_yyyy_yyyzz, tr_y_yyyy_yyzzz, tr_y_yyyy_yzzzz, tr_y_yyyy_zzzzz, tr_y_yyyyy_xxxx, tr_y_yyyyy_xxxxx, tr_y_yyyyy_xxxxy, tr_y_yyyyy_xxxxz, tr_y_yyyyy_xxxy, tr_y_yyyyy_xxxyy, tr_y_yyyyy_xxxyz, tr_y_yyyyy_xxxz, tr_y_yyyyy_xxxzz, tr_y_yyyyy_xxyy, tr_y_yyyyy_xxyyy, tr_y_yyyyy_xxyyz, tr_y_yyyyy_xxyz, tr_y_yyyyy_xxyzz, tr_y_yyyyy_xxzz, tr_y_yyyyy_xxzzz, tr_y_yyyyy_xyyy, tr_y_yyyyy_xyyyy, tr_y_yyyyy_xyyyz, tr_y_yyyyy_xyyz, tr_y_yyyyy_xyyzz, tr_y_yyyyy_xyzz, tr_y_yyyyy_xyzzz, tr_y_yyyyy_xzzz, tr_y_yyyyy_xzzzz, tr_y_yyyyy_yyyy, tr_y_yyyyy_yyyyy, tr_y_yyyyy_yyyyz, tr_y_yyyyy_yyyz, tr_y_yyyyy_yyyzz, tr_y_yyyyy_yyzz, tr_y_yyyyy_yyzzz, tr_y_yyyyy_yzzz, tr_y_yyyyy_yzzzz, tr_y_yyyyy_zzzz, tr_y_yyyyy_zzzzz, tr_y_yyyyyy_xxxxx, tr_y_yyyyyy_xxxxy, tr_y_yyyyyy_xxxxz, tr_y_yyyyyy_xxxyy, tr_y_yyyyyy_xxxyz, tr_y_yyyyyy_xxxzz, tr_y_yyyyyy_xxyyy, tr_y_yyyyyy_xxyyz, tr_y_yyyyyy_xxyzz, tr_y_yyyyyy_xxzzz, tr_y_yyyyyy_xyyyy, tr_y_yyyyyy_xyyyz, tr_y_yyyyyy_xyyzz, tr_y_yyyyyy_xyzzz, tr_y_yyyyyy_xzzzz, tr_y_yyyyyy_yyyyy, tr_y_yyyyyy_yyyyz, tr_y_yyyyyy_yyyzz, tr_y_yyyyyy_yyzzz, tr_y_yyyyyy_yzzzz, tr_y_yyyyyy_zzzzz, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxzz, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyy_xxxxx[i] = 5.0 * tr_y_yyyy_xxxxx[i] * fe_0 + ts_yyyyy_xxxxx[i] * fe_0 + tr_y_yyyyy_xxxxx[i] * pa_y[i];

        tr_y_yyyyyy_xxxxy[i] = 5.0 * tr_y_yyyy_xxxxy[i] * fe_0 + tr_y_yyyyy_xxxx[i] * fe_0 + ts_yyyyy_xxxxy[i] * fe_0 + tr_y_yyyyy_xxxxy[i] * pa_y[i];

        tr_y_yyyyyy_xxxxz[i] = 5.0 * tr_y_yyyy_xxxxz[i] * fe_0 + ts_yyyyy_xxxxz[i] * fe_0 + tr_y_yyyyy_xxxxz[i] * pa_y[i];

        tr_y_yyyyyy_xxxyy[i] = 5.0 * tr_y_yyyy_xxxyy[i] * fe_0 + 2.0 * tr_y_yyyyy_xxxy[i] * fe_0 + ts_yyyyy_xxxyy[i] * fe_0 + tr_y_yyyyy_xxxyy[i] * pa_y[i];

        tr_y_yyyyyy_xxxyz[i] = 5.0 * tr_y_yyyy_xxxyz[i] * fe_0 + tr_y_yyyyy_xxxz[i] * fe_0 + ts_yyyyy_xxxyz[i] * fe_0 + tr_y_yyyyy_xxxyz[i] * pa_y[i];

        tr_y_yyyyyy_xxxzz[i] = 5.0 * tr_y_yyyy_xxxzz[i] * fe_0 + ts_yyyyy_xxxzz[i] * fe_0 + tr_y_yyyyy_xxxzz[i] * pa_y[i];

        tr_y_yyyyyy_xxyyy[i] = 5.0 * tr_y_yyyy_xxyyy[i] * fe_0 + 3.0 * tr_y_yyyyy_xxyy[i] * fe_0 + ts_yyyyy_xxyyy[i] * fe_0 + tr_y_yyyyy_xxyyy[i] * pa_y[i];

        tr_y_yyyyyy_xxyyz[i] = 5.0 * tr_y_yyyy_xxyyz[i] * fe_0 + 2.0 * tr_y_yyyyy_xxyz[i] * fe_0 + ts_yyyyy_xxyyz[i] * fe_0 + tr_y_yyyyy_xxyyz[i] * pa_y[i];

        tr_y_yyyyyy_xxyzz[i] = 5.0 * tr_y_yyyy_xxyzz[i] * fe_0 + tr_y_yyyyy_xxzz[i] * fe_0 + ts_yyyyy_xxyzz[i] * fe_0 + tr_y_yyyyy_xxyzz[i] * pa_y[i];

        tr_y_yyyyyy_xxzzz[i] = 5.0 * tr_y_yyyy_xxzzz[i] * fe_0 + ts_yyyyy_xxzzz[i] * fe_0 + tr_y_yyyyy_xxzzz[i] * pa_y[i];

        tr_y_yyyyyy_xyyyy[i] = 5.0 * tr_y_yyyy_xyyyy[i] * fe_0 + 4.0 * tr_y_yyyyy_xyyy[i] * fe_0 + ts_yyyyy_xyyyy[i] * fe_0 + tr_y_yyyyy_xyyyy[i] * pa_y[i];

        tr_y_yyyyyy_xyyyz[i] = 5.0 * tr_y_yyyy_xyyyz[i] * fe_0 + 3.0 * tr_y_yyyyy_xyyz[i] * fe_0 + ts_yyyyy_xyyyz[i] * fe_0 + tr_y_yyyyy_xyyyz[i] * pa_y[i];

        tr_y_yyyyyy_xyyzz[i] = 5.0 * tr_y_yyyy_xyyzz[i] * fe_0 + 2.0 * tr_y_yyyyy_xyzz[i] * fe_0 + ts_yyyyy_xyyzz[i] * fe_0 + tr_y_yyyyy_xyyzz[i] * pa_y[i];

        tr_y_yyyyyy_xyzzz[i] = 5.0 * tr_y_yyyy_xyzzz[i] * fe_0 + tr_y_yyyyy_xzzz[i] * fe_0 + ts_yyyyy_xyzzz[i] * fe_0 + tr_y_yyyyy_xyzzz[i] * pa_y[i];

        tr_y_yyyyyy_xzzzz[i] = 5.0 * tr_y_yyyy_xzzzz[i] * fe_0 + ts_yyyyy_xzzzz[i] * fe_0 + tr_y_yyyyy_xzzzz[i] * pa_y[i];

        tr_y_yyyyyy_yyyyy[i] = 5.0 * tr_y_yyyy_yyyyy[i] * fe_0 + 5.0 * tr_y_yyyyy_yyyy[i] * fe_0 + ts_yyyyy_yyyyy[i] * fe_0 + tr_y_yyyyy_yyyyy[i] * pa_y[i];

        tr_y_yyyyyy_yyyyz[i] = 5.0 * tr_y_yyyy_yyyyz[i] * fe_0 + 4.0 * tr_y_yyyyy_yyyz[i] * fe_0 + ts_yyyyy_yyyyz[i] * fe_0 + tr_y_yyyyy_yyyyz[i] * pa_y[i];

        tr_y_yyyyyy_yyyzz[i] = 5.0 * tr_y_yyyy_yyyzz[i] * fe_0 + 3.0 * tr_y_yyyyy_yyzz[i] * fe_0 + ts_yyyyy_yyyzz[i] * fe_0 + tr_y_yyyyy_yyyzz[i] * pa_y[i];

        tr_y_yyyyyy_yyzzz[i] = 5.0 * tr_y_yyyy_yyzzz[i] * fe_0 + 2.0 * tr_y_yyyyy_yzzz[i] * fe_0 + ts_yyyyy_yyzzz[i] * fe_0 + tr_y_yyyyy_yyzzz[i] * pa_y[i];

        tr_y_yyyyyy_yzzzz[i] = 5.0 * tr_y_yyyy_yzzzz[i] * fe_0 + tr_y_yyyyy_zzzz[i] * fe_0 + ts_yyyyy_yzzzz[i] * fe_0 + tr_y_yyyyy_yzzzz[i] * pa_y[i];

        tr_y_yyyyyy_zzzzz[i] = 5.0 * tr_y_yyyy_zzzzz[i] * fe_0 + ts_yyyyy_zzzzz[i] * fe_0 + tr_y_yyyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 1050-1071 components of targeted buffer : IH

    auto tr_y_yyyyyz_xxxxx = pbuffer.data(idx_dip_ih + 1050);

    auto tr_y_yyyyyz_xxxxy = pbuffer.data(idx_dip_ih + 1051);

    auto tr_y_yyyyyz_xxxxz = pbuffer.data(idx_dip_ih + 1052);

    auto tr_y_yyyyyz_xxxyy = pbuffer.data(idx_dip_ih + 1053);

    auto tr_y_yyyyyz_xxxyz = pbuffer.data(idx_dip_ih + 1054);

    auto tr_y_yyyyyz_xxxzz = pbuffer.data(idx_dip_ih + 1055);

    auto tr_y_yyyyyz_xxyyy = pbuffer.data(idx_dip_ih + 1056);

    auto tr_y_yyyyyz_xxyyz = pbuffer.data(idx_dip_ih + 1057);

    auto tr_y_yyyyyz_xxyzz = pbuffer.data(idx_dip_ih + 1058);

    auto tr_y_yyyyyz_xxzzz = pbuffer.data(idx_dip_ih + 1059);

    auto tr_y_yyyyyz_xyyyy = pbuffer.data(idx_dip_ih + 1060);

    auto tr_y_yyyyyz_xyyyz = pbuffer.data(idx_dip_ih + 1061);

    auto tr_y_yyyyyz_xyyzz = pbuffer.data(idx_dip_ih + 1062);

    auto tr_y_yyyyyz_xyzzz = pbuffer.data(idx_dip_ih + 1063);

    auto tr_y_yyyyyz_xzzzz = pbuffer.data(idx_dip_ih + 1064);

    auto tr_y_yyyyyz_yyyyy = pbuffer.data(idx_dip_ih + 1065);

    auto tr_y_yyyyyz_yyyyz = pbuffer.data(idx_dip_ih + 1066);

    auto tr_y_yyyyyz_yyyzz = pbuffer.data(idx_dip_ih + 1067);

    auto tr_y_yyyyyz_yyzzz = pbuffer.data(idx_dip_ih + 1068);

    auto tr_y_yyyyyz_yzzzz = pbuffer.data(idx_dip_ih + 1069);

    auto tr_y_yyyyyz_zzzzz = pbuffer.data(idx_dip_ih + 1070);

    #pragma omp simd aligned(pa_z, tr_y_yyyyy_xxxx, tr_y_yyyyy_xxxxx, tr_y_yyyyy_xxxxy, tr_y_yyyyy_xxxxz, tr_y_yyyyy_xxxy, tr_y_yyyyy_xxxyy, tr_y_yyyyy_xxxyz, tr_y_yyyyy_xxxz, tr_y_yyyyy_xxxzz, tr_y_yyyyy_xxyy, tr_y_yyyyy_xxyyy, tr_y_yyyyy_xxyyz, tr_y_yyyyy_xxyz, tr_y_yyyyy_xxyzz, tr_y_yyyyy_xxzz, tr_y_yyyyy_xxzzz, tr_y_yyyyy_xyyy, tr_y_yyyyy_xyyyy, tr_y_yyyyy_xyyyz, tr_y_yyyyy_xyyz, tr_y_yyyyy_xyyzz, tr_y_yyyyy_xyzz, tr_y_yyyyy_xyzzz, tr_y_yyyyy_xzzz, tr_y_yyyyy_xzzzz, tr_y_yyyyy_yyyy, tr_y_yyyyy_yyyyy, tr_y_yyyyy_yyyyz, tr_y_yyyyy_yyyz, tr_y_yyyyy_yyyzz, tr_y_yyyyy_yyzz, tr_y_yyyyy_yyzzz, tr_y_yyyyy_yzzz, tr_y_yyyyy_yzzzz, tr_y_yyyyy_zzzz, tr_y_yyyyy_zzzzz, tr_y_yyyyyz_xxxxx, tr_y_yyyyyz_xxxxy, tr_y_yyyyyz_xxxxz, tr_y_yyyyyz_xxxyy, tr_y_yyyyyz_xxxyz, tr_y_yyyyyz_xxxzz, tr_y_yyyyyz_xxyyy, tr_y_yyyyyz_xxyyz, tr_y_yyyyyz_xxyzz, tr_y_yyyyyz_xxzzz, tr_y_yyyyyz_xyyyy, tr_y_yyyyyz_xyyyz, tr_y_yyyyyz_xyyzz, tr_y_yyyyyz_xyzzz, tr_y_yyyyyz_xzzzz, tr_y_yyyyyz_yyyyy, tr_y_yyyyyz_yyyyz, tr_y_yyyyyz_yyyzz, tr_y_yyyyyz_yyzzz, tr_y_yyyyyz_yzzzz, tr_y_yyyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyz_xxxxx[i] = tr_y_yyyyy_xxxxx[i] * pa_z[i];

        tr_y_yyyyyz_xxxxy[i] = tr_y_yyyyy_xxxxy[i] * pa_z[i];

        tr_y_yyyyyz_xxxxz[i] = tr_y_yyyyy_xxxx[i] * fe_0 + tr_y_yyyyy_xxxxz[i] * pa_z[i];

        tr_y_yyyyyz_xxxyy[i] = tr_y_yyyyy_xxxyy[i] * pa_z[i];

        tr_y_yyyyyz_xxxyz[i] = tr_y_yyyyy_xxxy[i] * fe_0 + tr_y_yyyyy_xxxyz[i] * pa_z[i];

        tr_y_yyyyyz_xxxzz[i] = 2.0 * tr_y_yyyyy_xxxz[i] * fe_0 + tr_y_yyyyy_xxxzz[i] * pa_z[i];

        tr_y_yyyyyz_xxyyy[i] = tr_y_yyyyy_xxyyy[i] * pa_z[i];

        tr_y_yyyyyz_xxyyz[i] = tr_y_yyyyy_xxyy[i] * fe_0 + tr_y_yyyyy_xxyyz[i] * pa_z[i];

        tr_y_yyyyyz_xxyzz[i] = 2.0 * tr_y_yyyyy_xxyz[i] * fe_0 + tr_y_yyyyy_xxyzz[i] * pa_z[i];

        tr_y_yyyyyz_xxzzz[i] = 3.0 * tr_y_yyyyy_xxzz[i] * fe_0 + tr_y_yyyyy_xxzzz[i] * pa_z[i];

        tr_y_yyyyyz_xyyyy[i] = tr_y_yyyyy_xyyyy[i] * pa_z[i];

        tr_y_yyyyyz_xyyyz[i] = tr_y_yyyyy_xyyy[i] * fe_0 + tr_y_yyyyy_xyyyz[i] * pa_z[i];

        tr_y_yyyyyz_xyyzz[i] = 2.0 * tr_y_yyyyy_xyyz[i] * fe_0 + tr_y_yyyyy_xyyzz[i] * pa_z[i];

        tr_y_yyyyyz_xyzzz[i] = 3.0 * tr_y_yyyyy_xyzz[i] * fe_0 + tr_y_yyyyy_xyzzz[i] * pa_z[i];

        tr_y_yyyyyz_xzzzz[i] = 4.0 * tr_y_yyyyy_xzzz[i] * fe_0 + tr_y_yyyyy_xzzzz[i] * pa_z[i];

        tr_y_yyyyyz_yyyyy[i] = tr_y_yyyyy_yyyyy[i] * pa_z[i];

        tr_y_yyyyyz_yyyyz[i] = tr_y_yyyyy_yyyy[i] * fe_0 + tr_y_yyyyy_yyyyz[i] * pa_z[i];

        tr_y_yyyyyz_yyyzz[i] = 2.0 * tr_y_yyyyy_yyyz[i] * fe_0 + tr_y_yyyyy_yyyzz[i] * pa_z[i];

        tr_y_yyyyyz_yyzzz[i] = 3.0 * tr_y_yyyyy_yyzz[i] * fe_0 + tr_y_yyyyy_yyzzz[i] * pa_z[i];

        tr_y_yyyyyz_yzzzz[i] = 4.0 * tr_y_yyyyy_yzzz[i] * fe_0 + tr_y_yyyyy_yzzzz[i] * pa_z[i];

        tr_y_yyyyyz_zzzzz[i] = 5.0 * tr_y_yyyyy_zzzz[i] * fe_0 + tr_y_yyyyy_zzzzz[i] * pa_z[i];
    }

    // Set up 1071-1092 components of targeted buffer : IH

    auto tr_y_yyyyzz_xxxxx = pbuffer.data(idx_dip_ih + 1071);

    auto tr_y_yyyyzz_xxxxy = pbuffer.data(idx_dip_ih + 1072);

    auto tr_y_yyyyzz_xxxxz = pbuffer.data(idx_dip_ih + 1073);

    auto tr_y_yyyyzz_xxxyy = pbuffer.data(idx_dip_ih + 1074);

    auto tr_y_yyyyzz_xxxyz = pbuffer.data(idx_dip_ih + 1075);

    auto tr_y_yyyyzz_xxxzz = pbuffer.data(idx_dip_ih + 1076);

    auto tr_y_yyyyzz_xxyyy = pbuffer.data(idx_dip_ih + 1077);

    auto tr_y_yyyyzz_xxyyz = pbuffer.data(idx_dip_ih + 1078);

    auto tr_y_yyyyzz_xxyzz = pbuffer.data(idx_dip_ih + 1079);

    auto tr_y_yyyyzz_xxzzz = pbuffer.data(idx_dip_ih + 1080);

    auto tr_y_yyyyzz_xyyyy = pbuffer.data(idx_dip_ih + 1081);

    auto tr_y_yyyyzz_xyyyz = pbuffer.data(idx_dip_ih + 1082);

    auto tr_y_yyyyzz_xyyzz = pbuffer.data(idx_dip_ih + 1083);

    auto tr_y_yyyyzz_xyzzz = pbuffer.data(idx_dip_ih + 1084);

    auto tr_y_yyyyzz_xzzzz = pbuffer.data(idx_dip_ih + 1085);

    auto tr_y_yyyyzz_yyyyy = pbuffer.data(idx_dip_ih + 1086);

    auto tr_y_yyyyzz_yyyyz = pbuffer.data(idx_dip_ih + 1087);

    auto tr_y_yyyyzz_yyyzz = pbuffer.data(idx_dip_ih + 1088);

    auto tr_y_yyyyzz_yyzzz = pbuffer.data(idx_dip_ih + 1089);

    auto tr_y_yyyyzz_yzzzz = pbuffer.data(idx_dip_ih + 1090);

    auto tr_y_yyyyzz_zzzzz = pbuffer.data(idx_dip_ih + 1091);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyyy_xxxxx, tr_y_yyyy_xxxxy, tr_y_yyyy_xxxyy, tr_y_yyyy_xxxyz, tr_y_yyyy_xxyyy, tr_y_yyyy_xxyyz, tr_y_yyyy_xxyzz, tr_y_yyyy_xyyyy, tr_y_yyyy_xyyyz, tr_y_yyyy_xyyzz, tr_y_yyyy_xyzzz, tr_y_yyyy_yyyyy, tr_y_yyyy_yyyyz, tr_y_yyyy_yyyzz, tr_y_yyyy_yyzzz, tr_y_yyyy_yzzzz, tr_y_yyyyz_xxxxx, tr_y_yyyyz_xxxxy, tr_y_yyyyz_xxxy, tr_y_yyyyz_xxxyy, tr_y_yyyyz_xxxyz, tr_y_yyyyz_xxyy, tr_y_yyyyz_xxyyy, tr_y_yyyyz_xxyyz, tr_y_yyyyz_xxyz, tr_y_yyyyz_xxyzz, tr_y_yyyyz_xyyy, tr_y_yyyyz_xyyyy, tr_y_yyyyz_xyyyz, tr_y_yyyyz_xyyz, tr_y_yyyyz_xyyzz, tr_y_yyyyz_xyzz, tr_y_yyyyz_xyzzz, tr_y_yyyyz_yyyy, tr_y_yyyyz_yyyyy, tr_y_yyyyz_yyyyz, tr_y_yyyyz_yyyz, tr_y_yyyyz_yyyzz, tr_y_yyyyz_yyzz, tr_y_yyyyz_yyzzz, tr_y_yyyyz_yzzz, tr_y_yyyyz_yzzzz, tr_y_yyyyzz_xxxxx, tr_y_yyyyzz_xxxxy, tr_y_yyyyzz_xxxxz, tr_y_yyyyzz_xxxyy, tr_y_yyyyzz_xxxyz, tr_y_yyyyzz_xxxzz, tr_y_yyyyzz_xxyyy, tr_y_yyyyzz_xxyyz, tr_y_yyyyzz_xxyzz, tr_y_yyyyzz_xxzzz, tr_y_yyyyzz_xyyyy, tr_y_yyyyzz_xyyyz, tr_y_yyyyzz_xyyzz, tr_y_yyyyzz_xyzzz, tr_y_yyyyzz_xzzzz, tr_y_yyyyzz_yyyyy, tr_y_yyyyzz_yyyyz, tr_y_yyyyzz_yyyzz, tr_y_yyyyzz_yyzzz, tr_y_yyyyzz_yzzzz, tr_y_yyyyzz_zzzzz, tr_y_yyyzz_xxxxz, tr_y_yyyzz_xxxzz, tr_y_yyyzz_xxzzz, tr_y_yyyzz_xzzzz, tr_y_yyyzz_zzzzz, tr_y_yyzz_xxxxz, tr_y_yyzz_xxxzz, tr_y_yyzz_xxzzz, tr_y_yyzz_xzzzz, tr_y_yyzz_zzzzz, ts_yyyzz_xxxxz, ts_yyyzz_xxxzz, ts_yyyzz_xxzzz, ts_yyyzz_xzzzz, ts_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyzz_xxxxx[i] = tr_y_yyyy_xxxxx[i] * fe_0 + tr_y_yyyyz_xxxxx[i] * pa_z[i];

        tr_y_yyyyzz_xxxxy[i] = tr_y_yyyy_xxxxy[i] * fe_0 + tr_y_yyyyz_xxxxy[i] * pa_z[i];

        tr_y_yyyyzz_xxxxz[i] = 3.0 * tr_y_yyzz_xxxxz[i] * fe_0 + ts_yyyzz_xxxxz[i] * fe_0 + tr_y_yyyzz_xxxxz[i] * pa_y[i];

        tr_y_yyyyzz_xxxyy[i] = tr_y_yyyy_xxxyy[i] * fe_0 + tr_y_yyyyz_xxxyy[i] * pa_z[i];

        tr_y_yyyyzz_xxxyz[i] = tr_y_yyyy_xxxyz[i] * fe_0 + tr_y_yyyyz_xxxy[i] * fe_0 + tr_y_yyyyz_xxxyz[i] * pa_z[i];

        tr_y_yyyyzz_xxxzz[i] = 3.0 * tr_y_yyzz_xxxzz[i] * fe_0 + ts_yyyzz_xxxzz[i] * fe_0 + tr_y_yyyzz_xxxzz[i] * pa_y[i];

        tr_y_yyyyzz_xxyyy[i] = tr_y_yyyy_xxyyy[i] * fe_0 + tr_y_yyyyz_xxyyy[i] * pa_z[i];

        tr_y_yyyyzz_xxyyz[i] = tr_y_yyyy_xxyyz[i] * fe_0 + tr_y_yyyyz_xxyy[i] * fe_0 + tr_y_yyyyz_xxyyz[i] * pa_z[i];

        tr_y_yyyyzz_xxyzz[i] = tr_y_yyyy_xxyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_xxyz[i] * fe_0 + tr_y_yyyyz_xxyzz[i] * pa_z[i];

        tr_y_yyyyzz_xxzzz[i] = 3.0 * tr_y_yyzz_xxzzz[i] * fe_0 + ts_yyyzz_xxzzz[i] * fe_0 + tr_y_yyyzz_xxzzz[i] * pa_y[i];

        tr_y_yyyyzz_xyyyy[i] = tr_y_yyyy_xyyyy[i] * fe_0 + tr_y_yyyyz_xyyyy[i] * pa_z[i];

        tr_y_yyyyzz_xyyyz[i] = tr_y_yyyy_xyyyz[i] * fe_0 + tr_y_yyyyz_xyyy[i] * fe_0 + tr_y_yyyyz_xyyyz[i] * pa_z[i];

        tr_y_yyyyzz_xyyzz[i] = tr_y_yyyy_xyyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_xyyz[i] * fe_0 + tr_y_yyyyz_xyyzz[i] * pa_z[i];

        tr_y_yyyyzz_xyzzz[i] = tr_y_yyyy_xyzzz[i] * fe_0 + 3.0 * tr_y_yyyyz_xyzz[i] * fe_0 + tr_y_yyyyz_xyzzz[i] * pa_z[i];

        tr_y_yyyyzz_xzzzz[i] = 3.0 * tr_y_yyzz_xzzzz[i] * fe_0 + ts_yyyzz_xzzzz[i] * fe_0 + tr_y_yyyzz_xzzzz[i] * pa_y[i];

        tr_y_yyyyzz_yyyyy[i] = tr_y_yyyy_yyyyy[i] * fe_0 + tr_y_yyyyz_yyyyy[i] * pa_z[i];

        tr_y_yyyyzz_yyyyz[i] = tr_y_yyyy_yyyyz[i] * fe_0 + tr_y_yyyyz_yyyy[i] * fe_0 + tr_y_yyyyz_yyyyz[i] * pa_z[i];

        tr_y_yyyyzz_yyyzz[i] = tr_y_yyyy_yyyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_yyyz[i] * fe_0 + tr_y_yyyyz_yyyzz[i] * pa_z[i];

        tr_y_yyyyzz_yyzzz[i] = tr_y_yyyy_yyzzz[i] * fe_0 + 3.0 * tr_y_yyyyz_yyzz[i] * fe_0 + tr_y_yyyyz_yyzzz[i] * pa_z[i];

        tr_y_yyyyzz_yzzzz[i] = tr_y_yyyy_yzzzz[i] * fe_0 + 4.0 * tr_y_yyyyz_yzzz[i] * fe_0 + tr_y_yyyyz_yzzzz[i] * pa_z[i];

        tr_y_yyyyzz_zzzzz[i] = 3.0 * tr_y_yyzz_zzzzz[i] * fe_0 + ts_yyyzz_zzzzz[i] * fe_0 + tr_y_yyyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1092-1113 components of targeted buffer : IH

    auto tr_y_yyyzzz_xxxxx = pbuffer.data(idx_dip_ih + 1092);

    auto tr_y_yyyzzz_xxxxy = pbuffer.data(idx_dip_ih + 1093);

    auto tr_y_yyyzzz_xxxxz = pbuffer.data(idx_dip_ih + 1094);

    auto tr_y_yyyzzz_xxxyy = pbuffer.data(idx_dip_ih + 1095);

    auto tr_y_yyyzzz_xxxyz = pbuffer.data(idx_dip_ih + 1096);

    auto tr_y_yyyzzz_xxxzz = pbuffer.data(idx_dip_ih + 1097);

    auto tr_y_yyyzzz_xxyyy = pbuffer.data(idx_dip_ih + 1098);

    auto tr_y_yyyzzz_xxyyz = pbuffer.data(idx_dip_ih + 1099);

    auto tr_y_yyyzzz_xxyzz = pbuffer.data(idx_dip_ih + 1100);

    auto tr_y_yyyzzz_xxzzz = pbuffer.data(idx_dip_ih + 1101);

    auto tr_y_yyyzzz_xyyyy = pbuffer.data(idx_dip_ih + 1102);

    auto tr_y_yyyzzz_xyyyz = pbuffer.data(idx_dip_ih + 1103);

    auto tr_y_yyyzzz_xyyzz = pbuffer.data(idx_dip_ih + 1104);

    auto tr_y_yyyzzz_xyzzz = pbuffer.data(idx_dip_ih + 1105);

    auto tr_y_yyyzzz_xzzzz = pbuffer.data(idx_dip_ih + 1106);

    auto tr_y_yyyzzz_yyyyy = pbuffer.data(idx_dip_ih + 1107);

    auto tr_y_yyyzzz_yyyyz = pbuffer.data(idx_dip_ih + 1108);

    auto tr_y_yyyzzz_yyyzz = pbuffer.data(idx_dip_ih + 1109);

    auto tr_y_yyyzzz_yyzzz = pbuffer.data(idx_dip_ih + 1110);

    auto tr_y_yyyzzz_yzzzz = pbuffer.data(idx_dip_ih + 1111);

    auto tr_y_yyyzzz_zzzzz = pbuffer.data(idx_dip_ih + 1112);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyyz_xxxxx, tr_y_yyyz_xxxxy, tr_y_yyyz_xxxyy, tr_y_yyyz_xxxyz, tr_y_yyyz_xxyyy, tr_y_yyyz_xxyyz, tr_y_yyyz_xxyzz, tr_y_yyyz_xyyyy, tr_y_yyyz_xyyyz, tr_y_yyyz_xyyzz, tr_y_yyyz_xyzzz, tr_y_yyyz_yyyyy, tr_y_yyyz_yyyyz, tr_y_yyyz_yyyzz, tr_y_yyyz_yyzzz, tr_y_yyyz_yzzzz, tr_y_yyyzz_xxxxx, tr_y_yyyzz_xxxxy, tr_y_yyyzz_xxxy, tr_y_yyyzz_xxxyy, tr_y_yyyzz_xxxyz, tr_y_yyyzz_xxyy, tr_y_yyyzz_xxyyy, tr_y_yyyzz_xxyyz, tr_y_yyyzz_xxyz, tr_y_yyyzz_xxyzz, tr_y_yyyzz_xyyy, tr_y_yyyzz_xyyyy, tr_y_yyyzz_xyyyz, tr_y_yyyzz_xyyz, tr_y_yyyzz_xyyzz, tr_y_yyyzz_xyzz, tr_y_yyyzz_xyzzz, tr_y_yyyzz_yyyy, tr_y_yyyzz_yyyyy, tr_y_yyyzz_yyyyz, tr_y_yyyzz_yyyz, tr_y_yyyzz_yyyzz, tr_y_yyyzz_yyzz, tr_y_yyyzz_yyzzz, tr_y_yyyzz_yzzz, tr_y_yyyzz_yzzzz, tr_y_yyyzzz_xxxxx, tr_y_yyyzzz_xxxxy, tr_y_yyyzzz_xxxxz, tr_y_yyyzzz_xxxyy, tr_y_yyyzzz_xxxyz, tr_y_yyyzzz_xxxzz, tr_y_yyyzzz_xxyyy, tr_y_yyyzzz_xxyyz, tr_y_yyyzzz_xxyzz, tr_y_yyyzzz_xxzzz, tr_y_yyyzzz_xyyyy, tr_y_yyyzzz_xyyyz, tr_y_yyyzzz_xyyzz, tr_y_yyyzzz_xyzzz, tr_y_yyyzzz_xzzzz, tr_y_yyyzzz_yyyyy, tr_y_yyyzzz_yyyyz, tr_y_yyyzzz_yyyzz, tr_y_yyyzzz_yyzzz, tr_y_yyyzzz_yzzzz, tr_y_yyyzzz_zzzzz, tr_y_yyzzz_xxxxz, tr_y_yyzzz_xxxzz, tr_y_yyzzz_xxzzz, tr_y_yyzzz_xzzzz, tr_y_yyzzz_zzzzz, tr_y_yzzz_xxxxz, tr_y_yzzz_xxxzz, tr_y_yzzz_xxzzz, tr_y_yzzz_xzzzz, tr_y_yzzz_zzzzz, ts_yyzzz_xxxxz, ts_yyzzz_xxxzz, ts_yyzzz_xxzzz, ts_yyzzz_xzzzz, ts_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzzz_xxxxx[i] = 2.0 * tr_y_yyyz_xxxxx[i] * fe_0 + tr_y_yyyzz_xxxxx[i] * pa_z[i];

        tr_y_yyyzzz_xxxxy[i] = 2.0 * tr_y_yyyz_xxxxy[i] * fe_0 + tr_y_yyyzz_xxxxy[i] * pa_z[i];

        tr_y_yyyzzz_xxxxz[i] = 2.0 * tr_y_yzzz_xxxxz[i] * fe_0 + ts_yyzzz_xxxxz[i] * fe_0 + tr_y_yyzzz_xxxxz[i] * pa_y[i];

        tr_y_yyyzzz_xxxyy[i] = 2.0 * tr_y_yyyz_xxxyy[i] * fe_0 + tr_y_yyyzz_xxxyy[i] * pa_z[i];

        tr_y_yyyzzz_xxxyz[i] = 2.0 * tr_y_yyyz_xxxyz[i] * fe_0 + tr_y_yyyzz_xxxy[i] * fe_0 + tr_y_yyyzz_xxxyz[i] * pa_z[i];

        tr_y_yyyzzz_xxxzz[i] = 2.0 * tr_y_yzzz_xxxzz[i] * fe_0 + ts_yyzzz_xxxzz[i] * fe_0 + tr_y_yyzzz_xxxzz[i] * pa_y[i];

        tr_y_yyyzzz_xxyyy[i] = 2.0 * tr_y_yyyz_xxyyy[i] * fe_0 + tr_y_yyyzz_xxyyy[i] * pa_z[i];

        tr_y_yyyzzz_xxyyz[i] = 2.0 * tr_y_yyyz_xxyyz[i] * fe_0 + tr_y_yyyzz_xxyy[i] * fe_0 + tr_y_yyyzz_xxyyz[i] * pa_z[i];

        tr_y_yyyzzz_xxyzz[i] = 2.0 * tr_y_yyyz_xxyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_xxyz[i] * fe_0 + tr_y_yyyzz_xxyzz[i] * pa_z[i];

        tr_y_yyyzzz_xxzzz[i] = 2.0 * tr_y_yzzz_xxzzz[i] * fe_0 + ts_yyzzz_xxzzz[i] * fe_0 + tr_y_yyzzz_xxzzz[i] * pa_y[i];

        tr_y_yyyzzz_xyyyy[i] = 2.0 * tr_y_yyyz_xyyyy[i] * fe_0 + tr_y_yyyzz_xyyyy[i] * pa_z[i];

        tr_y_yyyzzz_xyyyz[i] = 2.0 * tr_y_yyyz_xyyyz[i] * fe_0 + tr_y_yyyzz_xyyy[i] * fe_0 + tr_y_yyyzz_xyyyz[i] * pa_z[i];

        tr_y_yyyzzz_xyyzz[i] = 2.0 * tr_y_yyyz_xyyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_xyyz[i] * fe_0 + tr_y_yyyzz_xyyzz[i] * pa_z[i];

        tr_y_yyyzzz_xyzzz[i] = 2.0 * tr_y_yyyz_xyzzz[i] * fe_0 + 3.0 * tr_y_yyyzz_xyzz[i] * fe_0 + tr_y_yyyzz_xyzzz[i] * pa_z[i];

        tr_y_yyyzzz_xzzzz[i] = 2.0 * tr_y_yzzz_xzzzz[i] * fe_0 + ts_yyzzz_xzzzz[i] * fe_0 + tr_y_yyzzz_xzzzz[i] * pa_y[i];

        tr_y_yyyzzz_yyyyy[i] = 2.0 * tr_y_yyyz_yyyyy[i] * fe_0 + tr_y_yyyzz_yyyyy[i] * pa_z[i];

        tr_y_yyyzzz_yyyyz[i] = 2.0 * tr_y_yyyz_yyyyz[i] * fe_0 + tr_y_yyyzz_yyyy[i] * fe_0 + tr_y_yyyzz_yyyyz[i] * pa_z[i];

        tr_y_yyyzzz_yyyzz[i] = 2.0 * tr_y_yyyz_yyyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_yyyz[i] * fe_0 + tr_y_yyyzz_yyyzz[i] * pa_z[i];

        tr_y_yyyzzz_yyzzz[i] = 2.0 * tr_y_yyyz_yyzzz[i] * fe_0 + 3.0 * tr_y_yyyzz_yyzz[i] * fe_0 + tr_y_yyyzz_yyzzz[i] * pa_z[i];

        tr_y_yyyzzz_yzzzz[i] = 2.0 * tr_y_yyyz_yzzzz[i] * fe_0 + 4.0 * tr_y_yyyzz_yzzz[i] * fe_0 + tr_y_yyyzz_yzzzz[i] * pa_z[i];

        tr_y_yyyzzz_zzzzz[i] = 2.0 * tr_y_yzzz_zzzzz[i] * fe_0 + ts_yyzzz_zzzzz[i] * fe_0 + tr_y_yyzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1113-1134 components of targeted buffer : IH

    auto tr_y_yyzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1113);

    auto tr_y_yyzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1114);

    auto tr_y_yyzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1115);

    auto tr_y_yyzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1116);

    auto tr_y_yyzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1117);

    auto tr_y_yyzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1118);

    auto tr_y_yyzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1119);

    auto tr_y_yyzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1120);

    auto tr_y_yyzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1121);

    auto tr_y_yyzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1122);

    auto tr_y_yyzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1123);

    auto tr_y_yyzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1124);

    auto tr_y_yyzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1125);

    auto tr_y_yyzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1126);

    auto tr_y_yyzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1127);

    auto tr_y_yyzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1128);

    auto tr_y_yyzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1129);

    auto tr_y_yyzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1130);

    auto tr_y_yyzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1131);

    auto tr_y_yyzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1132);

    auto tr_y_yyzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1133);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyzz_xxxxx, tr_y_yyzz_xxxxy, tr_y_yyzz_xxxyy, tr_y_yyzz_xxxyz, tr_y_yyzz_xxyyy, tr_y_yyzz_xxyyz, tr_y_yyzz_xxyzz, tr_y_yyzz_xyyyy, tr_y_yyzz_xyyyz, tr_y_yyzz_xyyzz, tr_y_yyzz_xyzzz, tr_y_yyzz_yyyyy, tr_y_yyzz_yyyyz, tr_y_yyzz_yyyzz, tr_y_yyzz_yyzzz, tr_y_yyzz_yzzzz, tr_y_yyzzz_xxxxx, tr_y_yyzzz_xxxxy, tr_y_yyzzz_xxxy, tr_y_yyzzz_xxxyy, tr_y_yyzzz_xxxyz, tr_y_yyzzz_xxyy, tr_y_yyzzz_xxyyy, tr_y_yyzzz_xxyyz, tr_y_yyzzz_xxyz, tr_y_yyzzz_xxyzz, tr_y_yyzzz_xyyy, tr_y_yyzzz_xyyyy, tr_y_yyzzz_xyyyz, tr_y_yyzzz_xyyz, tr_y_yyzzz_xyyzz, tr_y_yyzzz_xyzz, tr_y_yyzzz_xyzzz, tr_y_yyzzz_yyyy, tr_y_yyzzz_yyyyy, tr_y_yyzzz_yyyyz, tr_y_yyzzz_yyyz, tr_y_yyzzz_yyyzz, tr_y_yyzzz_yyzz, tr_y_yyzzz_yyzzz, tr_y_yyzzz_yzzz, tr_y_yyzzz_yzzzz, tr_y_yyzzzz_xxxxx, tr_y_yyzzzz_xxxxy, tr_y_yyzzzz_xxxxz, tr_y_yyzzzz_xxxyy, tr_y_yyzzzz_xxxyz, tr_y_yyzzzz_xxxzz, tr_y_yyzzzz_xxyyy, tr_y_yyzzzz_xxyyz, tr_y_yyzzzz_xxyzz, tr_y_yyzzzz_xxzzz, tr_y_yyzzzz_xyyyy, tr_y_yyzzzz_xyyyz, tr_y_yyzzzz_xyyzz, tr_y_yyzzzz_xyzzz, tr_y_yyzzzz_xzzzz, tr_y_yyzzzz_yyyyy, tr_y_yyzzzz_yyyyz, tr_y_yyzzzz_yyyzz, tr_y_yyzzzz_yyzzz, tr_y_yyzzzz_yzzzz, tr_y_yyzzzz_zzzzz, tr_y_yzzzz_xxxxz, tr_y_yzzzz_xxxzz, tr_y_yzzzz_xxzzz, tr_y_yzzzz_xzzzz, tr_y_yzzzz_zzzzz, tr_y_zzzz_xxxxz, tr_y_zzzz_xxxzz, tr_y_zzzz_xxzzz, tr_y_zzzz_xzzzz, tr_y_zzzz_zzzzz, ts_yzzzz_xxxxz, ts_yzzzz_xxxzz, ts_yzzzz_xxzzz, ts_yzzzz_xzzzz, ts_yzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzzz_xxxxx[i] = 3.0 * tr_y_yyzz_xxxxx[i] * fe_0 + tr_y_yyzzz_xxxxx[i] * pa_z[i];

        tr_y_yyzzzz_xxxxy[i] = 3.0 * tr_y_yyzz_xxxxy[i] * fe_0 + tr_y_yyzzz_xxxxy[i] * pa_z[i];

        tr_y_yyzzzz_xxxxz[i] = tr_y_zzzz_xxxxz[i] * fe_0 + ts_yzzzz_xxxxz[i] * fe_0 + tr_y_yzzzz_xxxxz[i] * pa_y[i];

        tr_y_yyzzzz_xxxyy[i] = 3.0 * tr_y_yyzz_xxxyy[i] * fe_0 + tr_y_yyzzz_xxxyy[i] * pa_z[i];

        tr_y_yyzzzz_xxxyz[i] = 3.0 * tr_y_yyzz_xxxyz[i] * fe_0 + tr_y_yyzzz_xxxy[i] * fe_0 + tr_y_yyzzz_xxxyz[i] * pa_z[i];

        tr_y_yyzzzz_xxxzz[i] = tr_y_zzzz_xxxzz[i] * fe_0 + ts_yzzzz_xxxzz[i] * fe_0 + tr_y_yzzzz_xxxzz[i] * pa_y[i];

        tr_y_yyzzzz_xxyyy[i] = 3.0 * tr_y_yyzz_xxyyy[i] * fe_0 + tr_y_yyzzz_xxyyy[i] * pa_z[i];

        tr_y_yyzzzz_xxyyz[i] = 3.0 * tr_y_yyzz_xxyyz[i] * fe_0 + tr_y_yyzzz_xxyy[i] * fe_0 + tr_y_yyzzz_xxyyz[i] * pa_z[i];

        tr_y_yyzzzz_xxyzz[i] = 3.0 * tr_y_yyzz_xxyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_xxyz[i] * fe_0 + tr_y_yyzzz_xxyzz[i] * pa_z[i];

        tr_y_yyzzzz_xxzzz[i] = tr_y_zzzz_xxzzz[i] * fe_0 + ts_yzzzz_xxzzz[i] * fe_0 + tr_y_yzzzz_xxzzz[i] * pa_y[i];

        tr_y_yyzzzz_xyyyy[i] = 3.0 * tr_y_yyzz_xyyyy[i] * fe_0 + tr_y_yyzzz_xyyyy[i] * pa_z[i];

        tr_y_yyzzzz_xyyyz[i] = 3.0 * tr_y_yyzz_xyyyz[i] * fe_0 + tr_y_yyzzz_xyyy[i] * fe_0 + tr_y_yyzzz_xyyyz[i] * pa_z[i];

        tr_y_yyzzzz_xyyzz[i] = 3.0 * tr_y_yyzz_xyyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_xyyz[i] * fe_0 + tr_y_yyzzz_xyyzz[i] * pa_z[i];

        tr_y_yyzzzz_xyzzz[i] = 3.0 * tr_y_yyzz_xyzzz[i] * fe_0 + 3.0 * tr_y_yyzzz_xyzz[i] * fe_0 + tr_y_yyzzz_xyzzz[i] * pa_z[i];

        tr_y_yyzzzz_xzzzz[i] = tr_y_zzzz_xzzzz[i] * fe_0 + ts_yzzzz_xzzzz[i] * fe_0 + tr_y_yzzzz_xzzzz[i] * pa_y[i];

        tr_y_yyzzzz_yyyyy[i] = 3.0 * tr_y_yyzz_yyyyy[i] * fe_0 + tr_y_yyzzz_yyyyy[i] * pa_z[i];

        tr_y_yyzzzz_yyyyz[i] = 3.0 * tr_y_yyzz_yyyyz[i] * fe_0 + tr_y_yyzzz_yyyy[i] * fe_0 + tr_y_yyzzz_yyyyz[i] * pa_z[i];

        tr_y_yyzzzz_yyyzz[i] = 3.0 * tr_y_yyzz_yyyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_yyyz[i] * fe_0 + tr_y_yyzzz_yyyzz[i] * pa_z[i];

        tr_y_yyzzzz_yyzzz[i] = 3.0 * tr_y_yyzz_yyzzz[i] * fe_0 + 3.0 * tr_y_yyzzz_yyzz[i] * fe_0 + tr_y_yyzzz_yyzzz[i] * pa_z[i];

        tr_y_yyzzzz_yzzzz[i] = 3.0 * tr_y_yyzz_yzzzz[i] * fe_0 + 4.0 * tr_y_yyzzz_yzzz[i] * fe_0 + tr_y_yyzzz_yzzzz[i] * pa_z[i];

        tr_y_yyzzzz_zzzzz[i] = tr_y_zzzz_zzzzz[i] * fe_0 + ts_yzzzz_zzzzz[i] * fe_0 + tr_y_yzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1134-1155 components of targeted buffer : IH

    auto tr_y_yzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1134);

    auto tr_y_yzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1135);

    auto tr_y_yzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1136);

    auto tr_y_yzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1137);

    auto tr_y_yzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1138);

    auto tr_y_yzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1139);

    auto tr_y_yzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1140);

    auto tr_y_yzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1141);

    auto tr_y_yzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1142);

    auto tr_y_yzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1143);

    auto tr_y_yzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1144);

    auto tr_y_yzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1145);

    auto tr_y_yzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1146);

    auto tr_y_yzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1147);

    auto tr_y_yzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1148);

    auto tr_y_yzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1149);

    auto tr_y_yzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1150);

    auto tr_y_yzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1151);

    auto tr_y_yzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1152);

    auto tr_y_yzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1153);

    auto tr_y_yzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1154);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yzzz_xxxxy, tr_y_yzzz_xxxyy, tr_y_yzzz_xxyyy, tr_y_yzzz_xyyyy, tr_y_yzzz_yyyyy, tr_y_yzzzz_xxxxy, tr_y_yzzzz_xxxyy, tr_y_yzzzz_xxyyy, tr_y_yzzzz_xyyyy, tr_y_yzzzz_yyyyy, tr_y_yzzzzz_xxxxx, tr_y_yzzzzz_xxxxy, tr_y_yzzzzz_xxxxz, tr_y_yzzzzz_xxxyy, tr_y_yzzzzz_xxxyz, tr_y_yzzzzz_xxxzz, tr_y_yzzzzz_xxyyy, tr_y_yzzzzz_xxyyz, tr_y_yzzzzz_xxyzz, tr_y_yzzzzz_xxzzz, tr_y_yzzzzz_xyyyy, tr_y_yzzzzz_xyyyz, tr_y_yzzzzz_xyyzz, tr_y_yzzzzz_xyzzz, tr_y_yzzzzz_xzzzz, tr_y_yzzzzz_yyyyy, tr_y_yzzzzz_yyyyz, tr_y_yzzzzz_yyyzz, tr_y_yzzzzz_yyzzz, tr_y_yzzzzz_yzzzz, tr_y_yzzzzz_zzzzz, tr_y_zzzzz_xxxxx, tr_y_zzzzz_xxxxz, tr_y_zzzzz_xxxyz, tr_y_zzzzz_xxxz, tr_y_zzzzz_xxxzz, tr_y_zzzzz_xxyyz, tr_y_zzzzz_xxyz, tr_y_zzzzz_xxyzz, tr_y_zzzzz_xxzz, tr_y_zzzzz_xxzzz, tr_y_zzzzz_xyyyz, tr_y_zzzzz_xyyz, tr_y_zzzzz_xyyzz, tr_y_zzzzz_xyzz, tr_y_zzzzz_xyzzz, tr_y_zzzzz_xzzz, tr_y_zzzzz_xzzzz, tr_y_zzzzz_yyyyz, tr_y_zzzzz_yyyz, tr_y_zzzzz_yyyzz, tr_y_zzzzz_yyzz, tr_y_zzzzz_yyzzz, tr_y_zzzzz_yzzz, tr_y_zzzzz_yzzzz, tr_y_zzzzz_zzzz, tr_y_zzzzz_zzzzz, ts_zzzzz_xxxxx, ts_zzzzz_xxxxz, ts_zzzzz_xxxyz, ts_zzzzz_xxxzz, ts_zzzzz_xxyyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzzz_xxxxx[i] = ts_zzzzz_xxxxx[i] * fe_0 + tr_y_zzzzz_xxxxx[i] * pa_y[i];

        tr_y_yzzzzz_xxxxy[i] = 4.0 * tr_y_yzzz_xxxxy[i] * fe_0 + tr_y_yzzzz_xxxxy[i] * pa_z[i];

        tr_y_yzzzzz_xxxxz[i] = ts_zzzzz_xxxxz[i] * fe_0 + tr_y_zzzzz_xxxxz[i] * pa_y[i];

        tr_y_yzzzzz_xxxyy[i] = 4.0 * tr_y_yzzz_xxxyy[i] * fe_0 + tr_y_yzzzz_xxxyy[i] * pa_z[i];

        tr_y_yzzzzz_xxxyz[i] = tr_y_zzzzz_xxxz[i] * fe_0 + ts_zzzzz_xxxyz[i] * fe_0 + tr_y_zzzzz_xxxyz[i] * pa_y[i];

        tr_y_yzzzzz_xxxzz[i] = ts_zzzzz_xxxzz[i] * fe_0 + tr_y_zzzzz_xxxzz[i] * pa_y[i];

        tr_y_yzzzzz_xxyyy[i] = 4.0 * tr_y_yzzz_xxyyy[i] * fe_0 + tr_y_yzzzz_xxyyy[i] * pa_z[i];

        tr_y_yzzzzz_xxyyz[i] = 2.0 * tr_y_zzzzz_xxyz[i] * fe_0 + ts_zzzzz_xxyyz[i] * fe_0 + tr_y_zzzzz_xxyyz[i] * pa_y[i];

        tr_y_yzzzzz_xxyzz[i] = tr_y_zzzzz_xxzz[i] * fe_0 + ts_zzzzz_xxyzz[i] * fe_0 + tr_y_zzzzz_xxyzz[i] * pa_y[i];

        tr_y_yzzzzz_xxzzz[i] = ts_zzzzz_xxzzz[i] * fe_0 + tr_y_zzzzz_xxzzz[i] * pa_y[i];

        tr_y_yzzzzz_xyyyy[i] = 4.0 * tr_y_yzzz_xyyyy[i] * fe_0 + tr_y_yzzzz_xyyyy[i] * pa_z[i];

        tr_y_yzzzzz_xyyyz[i] = 3.0 * tr_y_zzzzz_xyyz[i] * fe_0 + ts_zzzzz_xyyyz[i] * fe_0 + tr_y_zzzzz_xyyyz[i] * pa_y[i];

        tr_y_yzzzzz_xyyzz[i] = 2.0 * tr_y_zzzzz_xyzz[i] * fe_0 + ts_zzzzz_xyyzz[i] * fe_0 + tr_y_zzzzz_xyyzz[i] * pa_y[i];

        tr_y_yzzzzz_xyzzz[i] = tr_y_zzzzz_xzzz[i] * fe_0 + ts_zzzzz_xyzzz[i] * fe_0 + tr_y_zzzzz_xyzzz[i] * pa_y[i];

        tr_y_yzzzzz_xzzzz[i] = ts_zzzzz_xzzzz[i] * fe_0 + tr_y_zzzzz_xzzzz[i] * pa_y[i];

        tr_y_yzzzzz_yyyyy[i] = 4.0 * tr_y_yzzz_yyyyy[i] * fe_0 + tr_y_yzzzz_yyyyy[i] * pa_z[i];

        tr_y_yzzzzz_yyyyz[i] = 4.0 * tr_y_zzzzz_yyyz[i] * fe_0 + ts_zzzzz_yyyyz[i] * fe_0 + tr_y_zzzzz_yyyyz[i] * pa_y[i];

        tr_y_yzzzzz_yyyzz[i] = 3.0 * tr_y_zzzzz_yyzz[i] * fe_0 + ts_zzzzz_yyyzz[i] * fe_0 + tr_y_zzzzz_yyyzz[i] * pa_y[i];

        tr_y_yzzzzz_yyzzz[i] = 2.0 * tr_y_zzzzz_yzzz[i] * fe_0 + ts_zzzzz_yyzzz[i] * fe_0 + tr_y_zzzzz_yyzzz[i] * pa_y[i];

        tr_y_yzzzzz_yzzzz[i] = tr_y_zzzzz_zzzz[i] * fe_0 + ts_zzzzz_yzzzz[i] * fe_0 + tr_y_zzzzz_yzzzz[i] * pa_y[i];

        tr_y_yzzzzz_zzzzz[i] = ts_zzzzz_zzzzz[i] * fe_0 + tr_y_zzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1155-1176 components of targeted buffer : IH

    auto tr_y_zzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1155);

    auto tr_y_zzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1156);

    auto tr_y_zzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1157);

    auto tr_y_zzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1158);

    auto tr_y_zzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1159);

    auto tr_y_zzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1160);

    auto tr_y_zzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1161);

    auto tr_y_zzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1162);

    auto tr_y_zzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1163);

    auto tr_y_zzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1164);

    auto tr_y_zzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1165);

    auto tr_y_zzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1166);

    auto tr_y_zzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1167);

    auto tr_y_zzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1168);

    auto tr_y_zzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1169);

    auto tr_y_zzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1170);

    auto tr_y_zzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1171);

    auto tr_y_zzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1172);

    auto tr_y_zzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1173);

    auto tr_y_zzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1174);

    auto tr_y_zzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1175);

    #pragma omp simd aligned(pa_z, tr_y_zzzz_xxxxx, tr_y_zzzz_xxxxy, tr_y_zzzz_xxxxz, tr_y_zzzz_xxxyy, tr_y_zzzz_xxxyz, tr_y_zzzz_xxxzz, tr_y_zzzz_xxyyy, tr_y_zzzz_xxyyz, tr_y_zzzz_xxyzz, tr_y_zzzz_xxzzz, tr_y_zzzz_xyyyy, tr_y_zzzz_xyyyz, tr_y_zzzz_xyyzz, tr_y_zzzz_xyzzz, tr_y_zzzz_xzzzz, tr_y_zzzz_yyyyy, tr_y_zzzz_yyyyz, tr_y_zzzz_yyyzz, tr_y_zzzz_yyzzz, tr_y_zzzz_yzzzz, tr_y_zzzz_zzzzz, tr_y_zzzzz_xxxx, tr_y_zzzzz_xxxxx, tr_y_zzzzz_xxxxy, tr_y_zzzzz_xxxxz, tr_y_zzzzz_xxxy, tr_y_zzzzz_xxxyy, tr_y_zzzzz_xxxyz, tr_y_zzzzz_xxxz, tr_y_zzzzz_xxxzz, tr_y_zzzzz_xxyy, tr_y_zzzzz_xxyyy, tr_y_zzzzz_xxyyz, tr_y_zzzzz_xxyz, tr_y_zzzzz_xxyzz, tr_y_zzzzz_xxzz, tr_y_zzzzz_xxzzz, tr_y_zzzzz_xyyy, tr_y_zzzzz_xyyyy, tr_y_zzzzz_xyyyz, tr_y_zzzzz_xyyz, tr_y_zzzzz_xyyzz, tr_y_zzzzz_xyzz, tr_y_zzzzz_xyzzz, tr_y_zzzzz_xzzz, tr_y_zzzzz_xzzzz, tr_y_zzzzz_yyyy, tr_y_zzzzz_yyyyy, tr_y_zzzzz_yyyyz, tr_y_zzzzz_yyyz, tr_y_zzzzz_yyyzz, tr_y_zzzzz_yyzz, tr_y_zzzzz_yyzzz, tr_y_zzzzz_yzzz, tr_y_zzzzz_yzzzz, tr_y_zzzzz_zzzz, tr_y_zzzzz_zzzzz, tr_y_zzzzzz_xxxxx, tr_y_zzzzzz_xxxxy, tr_y_zzzzzz_xxxxz, tr_y_zzzzzz_xxxyy, tr_y_zzzzzz_xxxyz, tr_y_zzzzzz_xxxzz, tr_y_zzzzzz_xxyyy, tr_y_zzzzzz_xxyyz, tr_y_zzzzzz_xxyzz, tr_y_zzzzzz_xxzzz, tr_y_zzzzzz_xyyyy, tr_y_zzzzzz_xyyyz, tr_y_zzzzzz_xyyzz, tr_y_zzzzzz_xyzzz, tr_y_zzzzzz_xzzzz, tr_y_zzzzzz_yyyyy, tr_y_zzzzzz_yyyyz, tr_y_zzzzzz_yyyzz, tr_y_zzzzzz_yyzzz, tr_y_zzzzzz_yzzzz, tr_y_zzzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzzz_xxxxx[i] = 5.0 * tr_y_zzzz_xxxxx[i] * fe_0 + tr_y_zzzzz_xxxxx[i] * pa_z[i];

        tr_y_zzzzzz_xxxxy[i] = 5.0 * tr_y_zzzz_xxxxy[i] * fe_0 + tr_y_zzzzz_xxxxy[i] * pa_z[i];

        tr_y_zzzzzz_xxxxz[i] = 5.0 * tr_y_zzzz_xxxxz[i] * fe_0 + tr_y_zzzzz_xxxx[i] * fe_0 + tr_y_zzzzz_xxxxz[i] * pa_z[i];

        tr_y_zzzzzz_xxxyy[i] = 5.0 * tr_y_zzzz_xxxyy[i] * fe_0 + tr_y_zzzzz_xxxyy[i] * pa_z[i];

        tr_y_zzzzzz_xxxyz[i] = 5.0 * tr_y_zzzz_xxxyz[i] * fe_0 + tr_y_zzzzz_xxxy[i] * fe_0 + tr_y_zzzzz_xxxyz[i] * pa_z[i];

        tr_y_zzzzzz_xxxzz[i] = 5.0 * tr_y_zzzz_xxxzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xxxz[i] * fe_0 + tr_y_zzzzz_xxxzz[i] * pa_z[i];

        tr_y_zzzzzz_xxyyy[i] = 5.0 * tr_y_zzzz_xxyyy[i] * fe_0 + tr_y_zzzzz_xxyyy[i] * pa_z[i];

        tr_y_zzzzzz_xxyyz[i] = 5.0 * tr_y_zzzz_xxyyz[i] * fe_0 + tr_y_zzzzz_xxyy[i] * fe_0 + tr_y_zzzzz_xxyyz[i] * pa_z[i];

        tr_y_zzzzzz_xxyzz[i] = 5.0 * tr_y_zzzz_xxyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xxyz[i] * fe_0 + tr_y_zzzzz_xxyzz[i] * pa_z[i];

        tr_y_zzzzzz_xxzzz[i] = 5.0 * tr_y_zzzz_xxzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_xxzz[i] * fe_0 + tr_y_zzzzz_xxzzz[i] * pa_z[i];

        tr_y_zzzzzz_xyyyy[i] = 5.0 * tr_y_zzzz_xyyyy[i] * fe_0 + tr_y_zzzzz_xyyyy[i] * pa_z[i];

        tr_y_zzzzzz_xyyyz[i] = 5.0 * tr_y_zzzz_xyyyz[i] * fe_0 + tr_y_zzzzz_xyyy[i] * fe_0 + tr_y_zzzzz_xyyyz[i] * pa_z[i];

        tr_y_zzzzzz_xyyzz[i] = 5.0 * tr_y_zzzz_xyyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xyyz[i] * fe_0 + tr_y_zzzzz_xyyzz[i] * pa_z[i];

        tr_y_zzzzzz_xyzzz[i] = 5.0 * tr_y_zzzz_xyzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_xyzz[i] * fe_0 + tr_y_zzzzz_xyzzz[i] * pa_z[i];

        tr_y_zzzzzz_xzzzz[i] = 5.0 * tr_y_zzzz_xzzzz[i] * fe_0 + 4.0 * tr_y_zzzzz_xzzz[i] * fe_0 + tr_y_zzzzz_xzzzz[i] * pa_z[i];

        tr_y_zzzzzz_yyyyy[i] = 5.0 * tr_y_zzzz_yyyyy[i] * fe_0 + tr_y_zzzzz_yyyyy[i] * pa_z[i];

        tr_y_zzzzzz_yyyyz[i] = 5.0 * tr_y_zzzz_yyyyz[i] * fe_0 + tr_y_zzzzz_yyyy[i] * fe_0 + tr_y_zzzzz_yyyyz[i] * pa_z[i];

        tr_y_zzzzzz_yyyzz[i] = 5.0 * tr_y_zzzz_yyyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_yyyz[i] * fe_0 + tr_y_zzzzz_yyyzz[i] * pa_z[i];

        tr_y_zzzzzz_yyzzz[i] = 5.0 * tr_y_zzzz_yyzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_yyzz[i] * fe_0 + tr_y_zzzzz_yyzzz[i] * pa_z[i];

        tr_y_zzzzzz_yzzzz[i] = 5.0 * tr_y_zzzz_yzzzz[i] * fe_0 + 4.0 * tr_y_zzzzz_yzzz[i] * fe_0 + tr_y_zzzzz_yzzzz[i] * pa_z[i];

        tr_y_zzzzzz_zzzzz[i] = 5.0 * tr_y_zzzz_zzzzz[i] * fe_0 + 5.0 * tr_y_zzzzz_zzzz[i] * fe_0 + tr_y_zzzzz_zzzzz[i] * pa_z[i];
    }

    // Set up 1176-1197 components of targeted buffer : IH

    auto tr_z_xxxxxx_xxxxx = pbuffer.data(idx_dip_ih + 1176);

    auto tr_z_xxxxxx_xxxxy = pbuffer.data(idx_dip_ih + 1177);

    auto tr_z_xxxxxx_xxxxz = pbuffer.data(idx_dip_ih + 1178);

    auto tr_z_xxxxxx_xxxyy = pbuffer.data(idx_dip_ih + 1179);

    auto tr_z_xxxxxx_xxxyz = pbuffer.data(idx_dip_ih + 1180);

    auto tr_z_xxxxxx_xxxzz = pbuffer.data(idx_dip_ih + 1181);

    auto tr_z_xxxxxx_xxyyy = pbuffer.data(idx_dip_ih + 1182);

    auto tr_z_xxxxxx_xxyyz = pbuffer.data(idx_dip_ih + 1183);

    auto tr_z_xxxxxx_xxyzz = pbuffer.data(idx_dip_ih + 1184);

    auto tr_z_xxxxxx_xxzzz = pbuffer.data(idx_dip_ih + 1185);

    auto tr_z_xxxxxx_xyyyy = pbuffer.data(idx_dip_ih + 1186);

    auto tr_z_xxxxxx_xyyyz = pbuffer.data(idx_dip_ih + 1187);

    auto tr_z_xxxxxx_xyyzz = pbuffer.data(idx_dip_ih + 1188);

    auto tr_z_xxxxxx_xyzzz = pbuffer.data(idx_dip_ih + 1189);

    auto tr_z_xxxxxx_xzzzz = pbuffer.data(idx_dip_ih + 1190);

    auto tr_z_xxxxxx_yyyyy = pbuffer.data(idx_dip_ih + 1191);

    auto tr_z_xxxxxx_yyyyz = pbuffer.data(idx_dip_ih + 1192);

    auto tr_z_xxxxxx_yyyzz = pbuffer.data(idx_dip_ih + 1193);

    auto tr_z_xxxxxx_yyzzz = pbuffer.data(idx_dip_ih + 1194);

    auto tr_z_xxxxxx_yzzzz = pbuffer.data(idx_dip_ih + 1195);

    auto tr_z_xxxxxx_zzzzz = pbuffer.data(idx_dip_ih + 1196);

    #pragma omp simd aligned(pa_x, tr_z_xxxx_xxxxx, tr_z_xxxx_xxxxy, tr_z_xxxx_xxxxz, tr_z_xxxx_xxxyy, tr_z_xxxx_xxxyz, tr_z_xxxx_xxxzz, tr_z_xxxx_xxyyy, tr_z_xxxx_xxyyz, tr_z_xxxx_xxyzz, tr_z_xxxx_xxzzz, tr_z_xxxx_xyyyy, tr_z_xxxx_xyyyz, tr_z_xxxx_xyyzz, tr_z_xxxx_xyzzz, tr_z_xxxx_xzzzz, tr_z_xxxx_yyyyy, tr_z_xxxx_yyyyz, tr_z_xxxx_yyyzz, tr_z_xxxx_yyzzz, tr_z_xxxx_yzzzz, tr_z_xxxx_zzzzz, tr_z_xxxxx_xxxx, tr_z_xxxxx_xxxxx, tr_z_xxxxx_xxxxy, tr_z_xxxxx_xxxxz, tr_z_xxxxx_xxxy, tr_z_xxxxx_xxxyy, tr_z_xxxxx_xxxyz, tr_z_xxxxx_xxxz, tr_z_xxxxx_xxxzz, tr_z_xxxxx_xxyy, tr_z_xxxxx_xxyyy, tr_z_xxxxx_xxyyz, tr_z_xxxxx_xxyz, tr_z_xxxxx_xxyzz, tr_z_xxxxx_xxzz, tr_z_xxxxx_xxzzz, tr_z_xxxxx_xyyy, tr_z_xxxxx_xyyyy, tr_z_xxxxx_xyyyz, tr_z_xxxxx_xyyz, tr_z_xxxxx_xyyzz, tr_z_xxxxx_xyzz, tr_z_xxxxx_xyzzz, tr_z_xxxxx_xzzz, tr_z_xxxxx_xzzzz, tr_z_xxxxx_yyyy, tr_z_xxxxx_yyyyy, tr_z_xxxxx_yyyyz, tr_z_xxxxx_yyyz, tr_z_xxxxx_yyyzz, tr_z_xxxxx_yyzz, tr_z_xxxxx_yyzzz, tr_z_xxxxx_yzzz, tr_z_xxxxx_yzzzz, tr_z_xxxxx_zzzz, tr_z_xxxxx_zzzzz, tr_z_xxxxxx_xxxxx, tr_z_xxxxxx_xxxxy, tr_z_xxxxxx_xxxxz, tr_z_xxxxxx_xxxyy, tr_z_xxxxxx_xxxyz, tr_z_xxxxxx_xxxzz, tr_z_xxxxxx_xxyyy, tr_z_xxxxxx_xxyyz, tr_z_xxxxxx_xxyzz, tr_z_xxxxxx_xxzzz, tr_z_xxxxxx_xyyyy, tr_z_xxxxxx_xyyyz, tr_z_xxxxxx_xyyzz, tr_z_xxxxxx_xyzzz, tr_z_xxxxxx_xzzzz, tr_z_xxxxxx_yyyyy, tr_z_xxxxxx_yyyyz, tr_z_xxxxxx_yyyzz, tr_z_xxxxxx_yyzzz, tr_z_xxxxxx_yzzzz, tr_z_xxxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxx_xxxxx[i] = 5.0 * tr_z_xxxx_xxxxx[i] * fe_0 + 5.0 * tr_z_xxxxx_xxxx[i] * fe_0 + tr_z_xxxxx_xxxxx[i] * pa_x[i];

        tr_z_xxxxxx_xxxxy[i] = 5.0 * tr_z_xxxx_xxxxy[i] * fe_0 + 4.0 * tr_z_xxxxx_xxxy[i] * fe_0 + tr_z_xxxxx_xxxxy[i] * pa_x[i];

        tr_z_xxxxxx_xxxxz[i] = 5.0 * tr_z_xxxx_xxxxz[i] * fe_0 + 4.0 * tr_z_xxxxx_xxxz[i] * fe_0 + tr_z_xxxxx_xxxxz[i] * pa_x[i];

        tr_z_xxxxxx_xxxyy[i] = 5.0 * tr_z_xxxx_xxxyy[i] * fe_0 + 3.0 * tr_z_xxxxx_xxyy[i] * fe_0 + tr_z_xxxxx_xxxyy[i] * pa_x[i];

        tr_z_xxxxxx_xxxyz[i] = 5.0 * tr_z_xxxx_xxxyz[i] * fe_0 + 3.0 * tr_z_xxxxx_xxyz[i] * fe_0 + tr_z_xxxxx_xxxyz[i] * pa_x[i];

        tr_z_xxxxxx_xxxzz[i] = 5.0 * tr_z_xxxx_xxxzz[i] * fe_0 + 3.0 * tr_z_xxxxx_xxzz[i] * fe_0 + tr_z_xxxxx_xxxzz[i] * pa_x[i];

        tr_z_xxxxxx_xxyyy[i] = 5.0 * tr_z_xxxx_xxyyy[i] * fe_0 + 2.0 * tr_z_xxxxx_xyyy[i] * fe_0 + tr_z_xxxxx_xxyyy[i] * pa_x[i];

        tr_z_xxxxxx_xxyyz[i] = 5.0 * tr_z_xxxx_xxyyz[i] * fe_0 + 2.0 * tr_z_xxxxx_xyyz[i] * fe_0 + tr_z_xxxxx_xxyyz[i] * pa_x[i];

        tr_z_xxxxxx_xxyzz[i] = 5.0 * tr_z_xxxx_xxyzz[i] * fe_0 + 2.0 * tr_z_xxxxx_xyzz[i] * fe_0 + tr_z_xxxxx_xxyzz[i] * pa_x[i];

        tr_z_xxxxxx_xxzzz[i] = 5.0 * tr_z_xxxx_xxzzz[i] * fe_0 + 2.0 * tr_z_xxxxx_xzzz[i] * fe_0 + tr_z_xxxxx_xxzzz[i] * pa_x[i];

        tr_z_xxxxxx_xyyyy[i] = 5.0 * tr_z_xxxx_xyyyy[i] * fe_0 + tr_z_xxxxx_yyyy[i] * fe_0 + tr_z_xxxxx_xyyyy[i] * pa_x[i];

        tr_z_xxxxxx_xyyyz[i] = 5.0 * tr_z_xxxx_xyyyz[i] * fe_0 + tr_z_xxxxx_yyyz[i] * fe_0 + tr_z_xxxxx_xyyyz[i] * pa_x[i];

        tr_z_xxxxxx_xyyzz[i] = 5.0 * tr_z_xxxx_xyyzz[i] * fe_0 + tr_z_xxxxx_yyzz[i] * fe_0 + tr_z_xxxxx_xyyzz[i] * pa_x[i];

        tr_z_xxxxxx_xyzzz[i] = 5.0 * tr_z_xxxx_xyzzz[i] * fe_0 + tr_z_xxxxx_yzzz[i] * fe_0 + tr_z_xxxxx_xyzzz[i] * pa_x[i];

        tr_z_xxxxxx_xzzzz[i] = 5.0 * tr_z_xxxx_xzzzz[i] * fe_0 + tr_z_xxxxx_zzzz[i] * fe_0 + tr_z_xxxxx_xzzzz[i] * pa_x[i];

        tr_z_xxxxxx_yyyyy[i] = 5.0 * tr_z_xxxx_yyyyy[i] * fe_0 + tr_z_xxxxx_yyyyy[i] * pa_x[i];

        tr_z_xxxxxx_yyyyz[i] = 5.0 * tr_z_xxxx_yyyyz[i] * fe_0 + tr_z_xxxxx_yyyyz[i] * pa_x[i];

        tr_z_xxxxxx_yyyzz[i] = 5.0 * tr_z_xxxx_yyyzz[i] * fe_0 + tr_z_xxxxx_yyyzz[i] * pa_x[i];

        tr_z_xxxxxx_yyzzz[i] = 5.0 * tr_z_xxxx_yyzzz[i] * fe_0 + tr_z_xxxxx_yyzzz[i] * pa_x[i];

        tr_z_xxxxxx_yzzzz[i] = 5.0 * tr_z_xxxx_yzzzz[i] * fe_0 + tr_z_xxxxx_yzzzz[i] * pa_x[i];

        tr_z_xxxxxx_zzzzz[i] = 5.0 * tr_z_xxxx_zzzzz[i] * fe_0 + tr_z_xxxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 1197-1218 components of targeted buffer : IH

    auto tr_z_xxxxxy_xxxxx = pbuffer.data(idx_dip_ih + 1197);

    auto tr_z_xxxxxy_xxxxy = pbuffer.data(idx_dip_ih + 1198);

    auto tr_z_xxxxxy_xxxxz = pbuffer.data(idx_dip_ih + 1199);

    auto tr_z_xxxxxy_xxxyy = pbuffer.data(idx_dip_ih + 1200);

    auto tr_z_xxxxxy_xxxyz = pbuffer.data(idx_dip_ih + 1201);

    auto tr_z_xxxxxy_xxxzz = pbuffer.data(idx_dip_ih + 1202);

    auto tr_z_xxxxxy_xxyyy = pbuffer.data(idx_dip_ih + 1203);

    auto tr_z_xxxxxy_xxyyz = pbuffer.data(idx_dip_ih + 1204);

    auto tr_z_xxxxxy_xxyzz = pbuffer.data(idx_dip_ih + 1205);

    auto tr_z_xxxxxy_xxzzz = pbuffer.data(idx_dip_ih + 1206);

    auto tr_z_xxxxxy_xyyyy = pbuffer.data(idx_dip_ih + 1207);

    auto tr_z_xxxxxy_xyyyz = pbuffer.data(idx_dip_ih + 1208);

    auto tr_z_xxxxxy_xyyzz = pbuffer.data(idx_dip_ih + 1209);

    auto tr_z_xxxxxy_xyzzz = pbuffer.data(idx_dip_ih + 1210);

    auto tr_z_xxxxxy_xzzzz = pbuffer.data(idx_dip_ih + 1211);

    auto tr_z_xxxxxy_yyyyy = pbuffer.data(idx_dip_ih + 1212);

    auto tr_z_xxxxxy_yyyyz = pbuffer.data(idx_dip_ih + 1213);

    auto tr_z_xxxxxy_yyyzz = pbuffer.data(idx_dip_ih + 1214);

    auto tr_z_xxxxxy_yyzzz = pbuffer.data(idx_dip_ih + 1215);

    auto tr_z_xxxxxy_yzzzz = pbuffer.data(idx_dip_ih + 1216);

    auto tr_z_xxxxxy_zzzzz = pbuffer.data(idx_dip_ih + 1217);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxx_xxxx, tr_z_xxxxx_xxxxx, tr_z_xxxxx_xxxxy, tr_z_xxxxx_xxxxz, tr_z_xxxxx_xxxy, tr_z_xxxxx_xxxyy, tr_z_xxxxx_xxxyz, tr_z_xxxxx_xxxz, tr_z_xxxxx_xxxzz, tr_z_xxxxx_xxyy, tr_z_xxxxx_xxyyy, tr_z_xxxxx_xxyyz, tr_z_xxxxx_xxyz, tr_z_xxxxx_xxyzz, tr_z_xxxxx_xxzz, tr_z_xxxxx_xxzzz, tr_z_xxxxx_xyyy, tr_z_xxxxx_xyyyy, tr_z_xxxxx_xyyyz, tr_z_xxxxx_xyyz, tr_z_xxxxx_xyyzz, tr_z_xxxxx_xyzz, tr_z_xxxxx_xyzzz, tr_z_xxxxx_xzzz, tr_z_xxxxx_xzzzz, tr_z_xxxxx_zzzzz, tr_z_xxxxxy_xxxxx, tr_z_xxxxxy_xxxxy, tr_z_xxxxxy_xxxxz, tr_z_xxxxxy_xxxyy, tr_z_xxxxxy_xxxyz, tr_z_xxxxxy_xxxzz, tr_z_xxxxxy_xxyyy, tr_z_xxxxxy_xxyyz, tr_z_xxxxxy_xxyzz, tr_z_xxxxxy_xxzzz, tr_z_xxxxxy_xyyyy, tr_z_xxxxxy_xyyyz, tr_z_xxxxxy_xyyzz, tr_z_xxxxxy_xyzzz, tr_z_xxxxxy_xzzzz, tr_z_xxxxxy_yyyyy, tr_z_xxxxxy_yyyyz, tr_z_xxxxxy_yyyzz, tr_z_xxxxxy_yyzzz, tr_z_xxxxxy_yzzzz, tr_z_xxxxxy_zzzzz, tr_z_xxxxy_yyyyy, tr_z_xxxxy_yyyyz, tr_z_xxxxy_yyyzz, tr_z_xxxxy_yyzzz, tr_z_xxxxy_yzzzz, tr_z_xxxy_yyyyy, tr_z_xxxy_yyyyz, tr_z_xxxy_yyyzz, tr_z_xxxy_yyzzz, tr_z_xxxy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxy_xxxxx[i] = tr_z_xxxxx_xxxxx[i] * pa_y[i];

        tr_z_xxxxxy_xxxxy[i] = tr_z_xxxxx_xxxx[i] * fe_0 + tr_z_xxxxx_xxxxy[i] * pa_y[i];

        tr_z_xxxxxy_xxxxz[i] = tr_z_xxxxx_xxxxz[i] * pa_y[i];

        tr_z_xxxxxy_xxxyy[i] = 2.0 * tr_z_xxxxx_xxxy[i] * fe_0 + tr_z_xxxxx_xxxyy[i] * pa_y[i];

        tr_z_xxxxxy_xxxyz[i] = tr_z_xxxxx_xxxz[i] * fe_0 + tr_z_xxxxx_xxxyz[i] * pa_y[i];

        tr_z_xxxxxy_xxxzz[i] = tr_z_xxxxx_xxxzz[i] * pa_y[i];

        tr_z_xxxxxy_xxyyy[i] = 3.0 * tr_z_xxxxx_xxyy[i] * fe_0 + tr_z_xxxxx_xxyyy[i] * pa_y[i];

        tr_z_xxxxxy_xxyyz[i] = 2.0 * tr_z_xxxxx_xxyz[i] * fe_0 + tr_z_xxxxx_xxyyz[i] * pa_y[i];

        tr_z_xxxxxy_xxyzz[i] = tr_z_xxxxx_xxzz[i] * fe_0 + tr_z_xxxxx_xxyzz[i] * pa_y[i];

        tr_z_xxxxxy_xxzzz[i] = tr_z_xxxxx_xxzzz[i] * pa_y[i];

        tr_z_xxxxxy_xyyyy[i] = 4.0 * tr_z_xxxxx_xyyy[i] * fe_0 + tr_z_xxxxx_xyyyy[i] * pa_y[i];

        tr_z_xxxxxy_xyyyz[i] = 3.0 * tr_z_xxxxx_xyyz[i] * fe_0 + tr_z_xxxxx_xyyyz[i] * pa_y[i];

        tr_z_xxxxxy_xyyzz[i] = 2.0 * tr_z_xxxxx_xyzz[i] * fe_0 + tr_z_xxxxx_xyyzz[i] * pa_y[i];

        tr_z_xxxxxy_xyzzz[i] = tr_z_xxxxx_xzzz[i] * fe_0 + tr_z_xxxxx_xyzzz[i] * pa_y[i];

        tr_z_xxxxxy_xzzzz[i] = tr_z_xxxxx_xzzzz[i] * pa_y[i];

        tr_z_xxxxxy_yyyyy[i] = 4.0 * tr_z_xxxy_yyyyy[i] * fe_0 + tr_z_xxxxy_yyyyy[i] * pa_x[i];

        tr_z_xxxxxy_yyyyz[i] = 4.0 * tr_z_xxxy_yyyyz[i] * fe_0 + tr_z_xxxxy_yyyyz[i] * pa_x[i];

        tr_z_xxxxxy_yyyzz[i] = 4.0 * tr_z_xxxy_yyyzz[i] * fe_0 + tr_z_xxxxy_yyyzz[i] * pa_x[i];

        tr_z_xxxxxy_yyzzz[i] = 4.0 * tr_z_xxxy_yyzzz[i] * fe_0 + tr_z_xxxxy_yyzzz[i] * pa_x[i];

        tr_z_xxxxxy_yzzzz[i] = 4.0 * tr_z_xxxy_yzzzz[i] * fe_0 + tr_z_xxxxy_yzzzz[i] * pa_x[i];

        tr_z_xxxxxy_zzzzz[i] = tr_z_xxxxx_zzzzz[i] * pa_y[i];
    }

    // Set up 1218-1239 components of targeted buffer : IH

    auto tr_z_xxxxxz_xxxxx = pbuffer.data(idx_dip_ih + 1218);

    auto tr_z_xxxxxz_xxxxy = pbuffer.data(idx_dip_ih + 1219);

    auto tr_z_xxxxxz_xxxxz = pbuffer.data(idx_dip_ih + 1220);

    auto tr_z_xxxxxz_xxxyy = pbuffer.data(idx_dip_ih + 1221);

    auto tr_z_xxxxxz_xxxyz = pbuffer.data(idx_dip_ih + 1222);

    auto tr_z_xxxxxz_xxxzz = pbuffer.data(idx_dip_ih + 1223);

    auto tr_z_xxxxxz_xxyyy = pbuffer.data(idx_dip_ih + 1224);

    auto tr_z_xxxxxz_xxyyz = pbuffer.data(idx_dip_ih + 1225);

    auto tr_z_xxxxxz_xxyzz = pbuffer.data(idx_dip_ih + 1226);

    auto tr_z_xxxxxz_xxzzz = pbuffer.data(idx_dip_ih + 1227);

    auto tr_z_xxxxxz_xyyyy = pbuffer.data(idx_dip_ih + 1228);

    auto tr_z_xxxxxz_xyyyz = pbuffer.data(idx_dip_ih + 1229);

    auto tr_z_xxxxxz_xyyzz = pbuffer.data(idx_dip_ih + 1230);

    auto tr_z_xxxxxz_xyzzz = pbuffer.data(idx_dip_ih + 1231);

    auto tr_z_xxxxxz_xzzzz = pbuffer.data(idx_dip_ih + 1232);

    auto tr_z_xxxxxz_yyyyy = pbuffer.data(idx_dip_ih + 1233);

    auto tr_z_xxxxxz_yyyyz = pbuffer.data(idx_dip_ih + 1234);

    auto tr_z_xxxxxz_yyyzz = pbuffer.data(idx_dip_ih + 1235);

    auto tr_z_xxxxxz_yyzzz = pbuffer.data(idx_dip_ih + 1236);

    auto tr_z_xxxxxz_yzzzz = pbuffer.data(idx_dip_ih + 1237);

    auto tr_z_xxxxxz_zzzzz = pbuffer.data(idx_dip_ih + 1238);

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xxxxx_xxxxx, tr_z_xxxxx_xxxxy, tr_z_xxxxx_xxxyy, tr_z_xxxxx_xxyyy, tr_z_xxxxx_xyyyy, tr_z_xxxxxz_xxxxx, tr_z_xxxxxz_xxxxy, tr_z_xxxxxz_xxxxz, tr_z_xxxxxz_xxxyy, tr_z_xxxxxz_xxxyz, tr_z_xxxxxz_xxxzz, tr_z_xxxxxz_xxyyy, tr_z_xxxxxz_xxyyz, tr_z_xxxxxz_xxyzz, tr_z_xxxxxz_xxzzz, tr_z_xxxxxz_xyyyy, tr_z_xxxxxz_xyyyz, tr_z_xxxxxz_xyyzz, tr_z_xxxxxz_xyzzz, tr_z_xxxxxz_xzzzz, tr_z_xxxxxz_yyyyy, tr_z_xxxxxz_yyyyz, tr_z_xxxxxz_yyyzz, tr_z_xxxxxz_yyzzz, tr_z_xxxxxz_yzzzz, tr_z_xxxxxz_zzzzz, tr_z_xxxxz_xxxxz, tr_z_xxxxz_xxxyz, tr_z_xxxxz_xxxz, tr_z_xxxxz_xxxzz, tr_z_xxxxz_xxyyz, tr_z_xxxxz_xxyz, tr_z_xxxxz_xxyzz, tr_z_xxxxz_xxzz, tr_z_xxxxz_xxzzz, tr_z_xxxxz_xyyyz, tr_z_xxxxz_xyyz, tr_z_xxxxz_xyyzz, tr_z_xxxxz_xyzz, tr_z_xxxxz_xyzzz, tr_z_xxxxz_xzzz, tr_z_xxxxz_xzzzz, tr_z_xxxxz_yyyyy, tr_z_xxxxz_yyyyz, tr_z_xxxxz_yyyz, tr_z_xxxxz_yyyzz, tr_z_xxxxz_yyzz, tr_z_xxxxz_yyzzz, tr_z_xxxxz_yzzz, tr_z_xxxxz_yzzzz, tr_z_xxxxz_zzzz, tr_z_xxxxz_zzzzz, tr_z_xxxz_xxxxz, tr_z_xxxz_xxxyz, tr_z_xxxz_xxxzz, tr_z_xxxz_xxyyz, tr_z_xxxz_xxyzz, tr_z_xxxz_xxzzz, tr_z_xxxz_xyyyz, tr_z_xxxz_xyyzz, tr_z_xxxz_xyzzz, tr_z_xxxz_xzzzz, tr_z_xxxz_yyyyy, tr_z_xxxz_yyyyz, tr_z_xxxz_yyyzz, tr_z_xxxz_yyzzz, tr_z_xxxz_yzzzz, tr_z_xxxz_zzzzz, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxyyy, ts_xxxxx_xyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxz_xxxxx[i] = ts_xxxxx_xxxxx[i] * fe_0 + tr_z_xxxxx_xxxxx[i] * pa_z[i];

        tr_z_xxxxxz_xxxxy[i] = ts_xxxxx_xxxxy[i] * fe_0 + tr_z_xxxxx_xxxxy[i] * pa_z[i];

        tr_z_xxxxxz_xxxxz[i] = 4.0 * tr_z_xxxz_xxxxz[i] * fe_0 + 4.0 * tr_z_xxxxz_xxxz[i] * fe_0 + tr_z_xxxxz_xxxxz[i] * pa_x[i];

        tr_z_xxxxxz_xxxyy[i] = ts_xxxxx_xxxyy[i] * fe_0 + tr_z_xxxxx_xxxyy[i] * pa_z[i];

        tr_z_xxxxxz_xxxyz[i] = 4.0 * tr_z_xxxz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxxxz_xxyz[i] * fe_0 + tr_z_xxxxz_xxxyz[i] * pa_x[i];

        tr_z_xxxxxz_xxxzz[i] = 4.0 * tr_z_xxxz_xxxzz[i] * fe_0 + 3.0 * tr_z_xxxxz_xxzz[i] * fe_0 + tr_z_xxxxz_xxxzz[i] * pa_x[i];

        tr_z_xxxxxz_xxyyy[i] = ts_xxxxx_xxyyy[i] * fe_0 + tr_z_xxxxx_xxyyy[i] * pa_z[i];

        tr_z_xxxxxz_xxyyz[i] = 4.0 * tr_z_xxxz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxxxz_xyyz[i] * fe_0 + tr_z_xxxxz_xxyyz[i] * pa_x[i];

        tr_z_xxxxxz_xxyzz[i] = 4.0 * tr_z_xxxz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxxxz_xyzz[i] * fe_0 + tr_z_xxxxz_xxyzz[i] * pa_x[i];

        tr_z_xxxxxz_xxzzz[i] = 4.0 * tr_z_xxxz_xxzzz[i] * fe_0 + 2.0 * tr_z_xxxxz_xzzz[i] * fe_0 + tr_z_xxxxz_xxzzz[i] * pa_x[i];

        tr_z_xxxxxz_xyyyy[i] = ts_xxxxx_xyyyy[i] * fe_0 + tr_z_xxxxx_xyyyy[i] * pa_z[i];

        tr_z_xxxxxz_xyyyz[i] = 4.0 * tr_z_xxxz_xyyyz[i] * fe_0 + tr_z_xxxxz_yyyz[i] * fe_0 + tr_z_xxxxz_xyyyz[i] * pa_x[i];

        tr_z_xxxxxz_xyyzz[i] = 4.0 * tr_z_xxxz_xyyzz[i] * fe_0 + tr_z_xxxxz_yyzz[i] * fe_0 + tr_z_xxxxz_xyyzz[i] * pa_x[i];

        tr_z_xxxxxz_xyzzz[i] = 4.0 * tr_z_xxxz_xyzzz[i] * fe_0 + tr_z_xxxxz_yzzz[i] * fe_0 + tr_z_xxxxz_xyzzz[i] * pa_x[i];

        tr_z_xxxxxz_xzzzz[i] = 4.0 * tr_z_xxxz_xzzzz[i] * fe_0 + tr_z_xxxxz_zzzz[i] * fe_0 + tr_z_xxxxz_xzzzz[i] * pa_x[i];

        tr_z_xxxxxz_yyyyy[i] = 4.0 * tr_z_xxxz_yyyyy[i] * fe_0 + tr_z_xxxxz_yyyyy[i] * pa_x[i];

        tr_z_xxxxxz_yyyyz[i] = 4.0 * tr_z_xxxz_yyyyz[i] * fe_0 + tr_z_xxxxz_yyyyz[i] * pa_x[i];

        tr_z_xxxxxz_yyyzz[i] = 4.0 * tr_z_xxxz_yyyzz[i] * fe_0 + tr_z_xxxxz_yyyzz[i] * pa_x[i];

        tr_z_xxxxxz_yyzzz[i] = 4.0 * tr_z_xxxz_yyzzz[i] * fe_0 + tr_z_xxxxz_yyzzz[i] * pa_x[i];

        tr_z_xxxxxz_yzzzz[i] = 4.0 * tr_z_xxxz_yzzzz[i] * fe_0 + tr_z_xxxxz_yzzzz[i] * pa_x[i];

        tr_z_xxxxxz_zzzzz[i] = 4.0 * tr_z_xxxz_zzzzz[i] * fe_0 + tr_z_xxxxz_zzzzz[i] * pa_x[i];
    }

    // Set up 1239-1260 components of targeted buffer : IH

    auto tr_z_xxxxyy_xxxxx = pbuffer.data(idx_dip_ih + 1239);

    auto tr_z_xxxxyy_xxxxy = pbuffer.data(idx_dip_ih + 1240);

    auto tr_z_xxxxyy_xxxxz = pbuffer.data(idx_dip_ih + 1241);

    auto tr_z_xxxxyy_xxxyy = pbuffer.data(idx_dip_ih + 1242);

    auto tr_z_xxxxyy_xxxyz = pbuffer.data(idx_dip_ih + 1243);

    auto tr_z_xxxxyy_xxxzz = pbuffer.data(idx_dip_ih + 1244);

    auto tr_z_xxxxyy_xxyyy = pbuffer.data(idx_dip_ih + 1245);

    auto tr_z_xxxxyy_xxyyz = pbuffer.data(idx_dip_ih + 1246);

    auto tr_z_xxxxyy_xxyzz = pbuffer.data(idx_dip_ih + 1247);

    auto tr_z_xxxxyy_xxzzz = pbuffer.data(idx_dip_ih + 1248);

    auto tr_z_xxxxyy_xyyyy = pbuffer.data(idx_dip_ih + 1249);

    auto tr_z_xxxxyy_xyyyz = pbuffer.data(idx_dip_ih + 1250);

    auto tr_z_xxxxyy_xyyzz = pbuffer.data(idx_dip_ih + 1251);

    auto tr_z_xxxxyy_xyzzz = pbuffer.data(idx_dip_ih + 1252);

    auto tr_z_xxxxyy_xzzzz = pbuffer.data(idx_dip_ih + 1253);

    auto tr_z_xxxxyy_yyyyy = pbuffer.data(idx_dip_ih + 1254);

    auto tr_z_xxxxyy_yyyyz = pbuffer.data(idx_dip_ih + 1255);

    auto tr_z_xxxxyy_yyyzz = pbuffer.data(idx_dip_ih + 1256);

    auto tr_z_xxxxyy_yyzzz = pbuffer.data(idx_dip_ih + 1257);

    auto tr_z_xxxxyy_yzzzz = pbuffer.data(idx_dip_ih + 1258);

    auto tr_z_xxxxyy_zzzzz = pbuffer.data(idx_dip_ih + 1259);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxx_xxxxx, tr_z_xxxx_xxxxz, tr_z_xxxx_xxxzz, tr_z_xxxx_xxzzz, tr_z_xxxx_xzzzz, tr_z_xxxxy_xxxxx, tr_z_xxxxy_xxxxz, tr_z_xxxxy_xxxzz, tr_z_xxxxy_xxzzz, tr_z_xxxxy_xzzzz, tr_z_xxxxyy_xxxxx, tr_z_xxxxyy_xxxxy, tr_z_xxxxyy_xxxxz, tr_z_xxxxyy_xxxyy, tr_z_xxxxyy_xxxyz, tr_z_xxxxyy_xxxzz, tr_z_xxxxyy_xxyyy, tr_z_xxxxyy_xxyyz, tr_z_xxxxyy_xxyzz, tr_z_xxxxyy_xxzzz, tr_z_xxxxyy_xyyyy, tr_z_xxxxyy_xyyyz, tr_z_xxxxyy_xyyzz, tr_z_xxxxyy_xyzzz, tr_z_xxxxyy_xzzzz, tr_z_xxxxyy_yyyyy, tr_z_xxxxyy_yyyyz, tr_z_xxxxyy_yyyzz, tr_z_xxxxyy_yyzzz, tr_z_xxxxyy_yzzzz, tr_z_xxxxyy_zzzzz, tr_z_xxxyy_xxxxy, tr_z_xxxyy_xxxy, tr_z_xxxyy_xxxyy, tr_z_xxxyy_xxxyz, tr_z_xxxyy_xxyy, tr_z_xxxyy_xxyyy, tr_z_xxxyy_xxyyz, tr_z_xxxyy_xxyz, tr_z_xxxyy_xxyzz, tr_z_xxxyy_xyyy, tr_z_xxxyy_xyyyy, tr_z_xxxyy_xyyyz, tr_z_xxxyy_xyyz, tr_z_xxxyy_xyyzz, tr_z_xxxyy_xyzz, tr_z_xxxyy_xyzzz, tr_z_xxxyy_yyyy, tr_z_xxxyy_yyyyy, tr_z_xxxyy_yyyyz, tr_z_xxxyy_yyyz, tr_z_xxxyy_yyyzz, tr_z_xxxyy_yyzz, tr_z_xxxyy_yyzzz, tr_z_xxxyy_yzzz, tr_z_xxxyy_yzzzz, tr_z_xxxyy_zzzzz, tr_z_xxyy_xxxxy, tr_z_xxyy_xxxyy, tr_z_xxyy_xxxyz, tr_z_xxyy_xxyyy, tr_z_xxyy_xxyyz, tr_z_xxyy_xxyzz, tr_z_xxyy_xyyyy, tr_z_xxyy_xyyyz, tr_z_xxyy_xyyzz, tr_z_xxyy_xyzzz, tr_z_xxyy_yyyyy, tr_z_xxyy_yyyyz, tr_z_xxyy_yyyzz, tr_z_xxyy_yyzzz, tr_z_xxyy_yzzzz, tr_z_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyy_xxxxx[i] = tr_z_xxxx_xxxxx[i] * fe_0 + tr_z_xxxxy_xxxxx[i] * pa_y[i];

        tr_z_xxxxyy_xxxxy[i] = 3.0 * tr_z_xxyy_xxxxy[i] * fe_0 + 4.0 * tr_z_xxxyy_xxxy[i] * fe_0 + tr_z_xxxyy_xxxxy[i] * pa_x[i];

        tr_z_xxxxyy_xxxxz[i] = tr_z_xxxx_xxxxz[i] * fe_0 + tr_z_xxxxy_xxxxz[i] * pa_y[i];

        tr_z_xxxxyy_xxxyy[i] = 3.0 * tr_z_xxyy_xxxyy[i] * fe_0 + 3.0 * tr_z_xxxyy_xxyy[i] * fe_0 + tr_z_xxxyy_xxxyy[i] * pa_x[i];

        tr_z_xxxxyy_xxxyz[i] = 3.0 * tr_z_xxyy_xxxyz[i] * fe_0 + 3.0 * tr_z_xxxyy_xxyz[i] * fe_0 + tr_z_xxxyy_xxxyz[i] * pa_x[i];

        tr_z_xxxxyy_xxxzz[i] = tr_z_xxxx_xxxzz[i] * fe_0 + tr_z_xxxxy_xxxzz[i] * pa_y[i];

        tr_z_xxxxyy_xxyyy[i] = 3.0 * tr_z_xxyy_xxyyy[i] * fe_0 + 2.0 * tr_z_xxxyy_xyyy[i] * fe_0 + tr_z_xxxyy_xxyyy[i] * pa_x[i];

        tr_z_xxxxyy_xxyyz[i] = 3.0 * tr_z_xxyy_xxyyz[i] * fe_0 + 2.0 * tr_z_xxxyy_xyyz[i] * fe_0 + tr_z_xxxyy_xxyyz[i] * pa_x[i];

        tr_z_xxxxyy_xxyzz[i] = 3.0 * tr_z_xxyy_xxyzz[i] * fe_0 + 2.0 * tr_z_xxxyy_xyzz[i] * fe_0 + tr_z_xxxyy_xxyzz[i] * pa_x[i];

        tr_z_xxxxyy_xxzzz[i] = tr_z_xxxx_xxzzz[i] * fe_0 + tr_z_xxxxy_xxzzz[i] * pa_y[i];

        tr_z_xxxxyy_xyyyy[i] = 3.0 * tr_z_xxyy_xyyyy[i] * fe_0 + tr_z_xxxyy_yyyy[i] * fe_0 + tr_z_xxxyy_xyyyy[i] * pa_x[i];

        tr_z_xxxxyy_xyyyz[i] = 3.0 * tr_z_xxyy_xyyyz[i] * fe_0 + tr_z_xxxyy_yyyz[i] * fe_0 + tr_z_xxxyy_xyyyz[i] * pa_x[i];

        tr_z_xxxxyy_xyyzz[i] = 3.0 * tr_z_xxyy_xyyzz[i] * fe_0 + tr_z_xxxyy_yyzz[i] * fe_0 + tr_z_xxxyy_xyyzz[i] * pa_x[i];

        tr_z_xxxxyy_xyzzz[i] = 3.0 * tr_z_xxyy_xyzzz[i] * fe_0 + tr_z_xxxyy_yzzz[i] * fe_0 + tr_z_xxxyy_xyzzz[i] * pa_x[i];

        tr_z_xxxxyy_xzzzz[i] = tr_z_xxxx_xzzzz[i] * fe_0 + tr_z_xxxxy_xzzzz[i] * pa_y[i];

        tr_z_xxxxyy_yyyyy[i] = 3.0 * tr_z_xxyy_yyyyy[i] * fe_0 + tr_z_xxxyy_yyyyy[i] * pa_x[i];

        tr_z_xxxxyy_yyyyz[i] = 3.0 * tr_z_xxyy_yyyyz[i] * fe_0 + tr_z_xxxyy_yyyyz[i] * pa_x[i];

        tr_z_xxxxyy_yyyzz[i] = 3.0 * tr_z_xxyy_yyyzz[i] * fe_0 + tr_z_xxxyy_yyyzz[i] * pa_x[i];

        tr_z_xxxxyy_yyzzz[i] = 3.0 * tr_z_xxyy_yyzzz[i] * fe_0 + tr_z_xxxyy_yyzzz[i] * pa_x[i];

        tr_z_xxxxyy_yzzzz[i] = 3.0 * tr_z_xxyy_yzzzz[i] * fe_0 + tr_z_xxxyy_yzzzz[i] * pa_x[i];

        tr_z_xxxxyy_zzzzz[i] = 3.0 * tr_z_xxyy_zzzzz[i] * fe_0 + tr_z_xxxyy_zzzzz[i] * pa_x[i];
    }

    // Set up 1260-1281 components of targeted buffer : IH

    auto tr_z_xxxxyz_xxxxx = pbuffer.data(idx_dip_ih + 1260);

    auto tr_z_xxxxyz_xxxxy = pbuffer.data(idx_dip_ih + 1261);

    auto tr_z_xxxxyz_xxxxz = pbuffer.data(idx_dip_ih + 1262);

    auto tr_z_xxxxyz_xxxyy = pbuffer.data(idx_dip_ih + 1263);

    auto tr_z_xxxxyz_xxxyz = pbuffer.data(idx_dip_ih + 1264);

    auto tr_z_xxxxyz_xxxzz = pbuffer.data(idx_dip_ih + 1265);

    auto tr_z_xxxxyz_xxyyy = pbuffer.data(idx_dip_ih + 1266);

    auto tr_z_xxxxyz_xxyyz = pbuffer.data(idx_dip_ih + 1267);

    auto tr_z_xxxxyz_xxyzz = pbuffer.data(idx_dip_ih + 1268);

    auto tr_z_xxxxyz_xxzzz = pbuffer.data(idx_dip_ih + 1269);

    auto tr_z_xxxxyz_xyyyy = pbuffer.data(idx_dip_ih + 1270);

    auto tr_z_xxxxyz_xyyyz = pbuffer.data(idx_dip_ih + 1271);

    auto tr_z_xxxxyz_xyyzz = pbuffer.data(idx_dip_ih + 1272);

    auto tr_z_xxxxyz_xyzzz = pbuffer.data(idx_dip_ih + 1273);

    auto tr_z_xxxxyz_xzzzz = pbuffer.data(idx_dip_ih + 1274);

    auto tr_z_xxxxyz_yyyyy = pbuffer.data(idx_dip_ih + 1275);

    auto tr_z_xxxxyz_yyyyz = pbuffer.data(idx_dip_ih + 1276);

    auto tr_z_xxxxyz_yyyzz = pbuffer.data(idx_dip_ih + 1277);

    auto tr_z_xxxxyz_yyzzz = pbuffer.data(idx_dip_ih + 1278);

    auto tr_z_xxxxyz_yzzzz = pbuffer.data(idx_dip_ih + 1279);

    auto tr_z_xxxxyz_zzzzz = pbuffer.data(idx_dip_ih + 1280);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxyz_xxxxx, tr_z_xxxxyz_xxxxy, tr_z_xxxxyz_xxxxz, tr_z_xxxxyz_xxxyy, tr_z_xxxxyz_xxxyz, tr_z_xxxxyz_xxxzz, tr_z_xxxxyz_xxyyy, tr_z_xxxxyz_xxyyz, tr_z_xxxxyz_xxyzz, tr_z_xxxxyz_xxzzz, tr_z_xxxxyz_xyyyy, tr_z_xxxxyz_xyyyz, tr_z_xxxxyz_xyyzz, tr_z_xxxxyz_xyzzz, tr_z_xxxxyz_xzzzz, tr_z_xxxxyz_yyyyy, tr_z_xxxxyz_yyyyz, tr_z_xxxxyz_yyyzz, tr_z_xxxxyz_yyzzz, tr_z_xxxxyz_yzzzz, tr_z_xxxxyz_zzzzz, tr_z_xxxxz_xxxx, tr_z_xxxxz_xxxxx, tr_z_xxxxz_xxxxy, tr_z_xxxxz_xxxxz, tr_z_xxxxz_xxxy, tr_z_xxxxz_xxxyy, tr_z_xxxxz_xxxyz, tr_z_xxxxz_xxxz, tr_z_xxxxz_xxxzz, tr_z_xxxxz_xxyy, tr_z_xxxxz_xxyyy, tr_z_xxxxz_xxyyz, tr_z_xxxxz_xxyz, tr_z_xxxxz_xxyzz, tr_z_xxxxz_xxzz, tr_z_xxxxz_xxzzz, tr_z_xxxxz_xyyy, tr_z_xxxxz_xyyyy, tr_z_xxxxz_xyyyz, tr_z_xxxxz_xyyz, tr_z_xxxxz_xyyzz, tr_z_xxxxz_xyzz, tr_z_xxxxz_xyzzz, tr_z_xxxxz_xzzz, tr_z_xxxxz_xzzzz, tr_z_xxxxz_zzzzz, tr_z_xxxyz_yyyyy, tr_z_xxxyz_yyyyz, tr_z_xxxyz_yyyzz, tr_z_xxxyz_yyzzz, tr_z_xxxyz_yzzzz, tr_z_xxyz_yyyyy, tr_z_xxyz_yyyyz, tr_z_xxyz_yyyzz, tr_z_xxyz_yyzzz, tr_z_xxyz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyz_xxxxx[i] = tr_z_xxxxz_xxxxx[i] * pa_y[i];

        tr_z_xxxxyz_xxxxy[i] = tr_z_xxxxz_xxxx[i] * fe_0 + tr_z_xxxxz_xxxxy[i] * pa_y[i];

        tr_z_xxxxyz_xxxxz[i] = tr_z_xxxxz_xxxxz[i] * pa_y[i];

        tr_z_xxxxyz_xxxyy[i] = 2.0 * tr_z_xxxxz_xxxy[i] * fe_0 + tr_z_xxxxz_xxxyy[i] * pa_y[i];

        tr_z_xxxxyz_xxxyz[i] = tr_z_xxxxz_xxxz[i] * fe_0 + tr_z_xxxxz_xxxyz[i] * pa_y[i];

        tr_z_xxxxyz_xxxzz[i] = tr_z_xxxxz_xxxzz[i] * pa_y[i];

        tr_z_xxxxyz_xxyyy[i] = 3.0 * tr_z_xxxxz_xxyy[i] * fe_0 + tr_z_xxxxz_xxyyy[i] * pa_y[i];

        tr_z_xxxxyz_xxyyz[i] = 2.0 * tr_z_xxxxz_xxyz[i] * fe_0 + tr_z_xxxxz_xxyyz[i] * pa_y[i];

        tr_z_xxxxyz_xxyzz[i] = tr_z_xxxxz_xxzz[i] * fe_0 + tr_z_xxxxz_xxyzz[i] * pa_y[i];

        tr_z_xxxxyz_xxzzz[i] = tr_z_xxxxz_xxzzz[i] * pa_y[i];

        tr_z_xxxxyz_xyyyy[i] = 4.0 * tr_z_xxxxz_xyyy[i] * fe_0 + tr_z_xxxxz_xyyyy[i] * pa_y[i];

        tr_z_xxxxyz_xyyyz[i] = 3.0 * tr_z_xxxxz_xyyz[i] * fe_0 + tr_z_xxxxz_xyyyz[i] * pa_y[i];

        tr_z_xxxxyz_xyyzz[i] = 2.0 * tr_z_xxxxz_xyzz[i] * fe_0 + tr_z_xxxxz_xyyzz[i] * pa_y[i];

        tr_z_xxxxyz_xyzzz[i] = tr_z_xxxxz_xzzz[i] * fe_0 + tr_z_xxxxz_xyzzz[i] * pa_y[i];

        tr_z_xxxxyz_xzzzz[i] = tr_z_xxxxz_xzzzz[i] * pa_y[i];

        tr_z_xxxxyz_yyyyy[i] = 3.0 * tr_z_xxyz_yyyyy[i] * fe_0 + tr_z_xxxyz_yyyyy[i] * pa_x[i];

        tr_z_xxxxyz_yyyyz[i] = 3.0 * tr_z_xxyz_yyyyz[i] * fe_0 + tr_z_xxxyz_yyyyz[i] * pa_x[i];

        tr_z_xxxxyz_yyyzz[i] = 3.0 * tr_z_xxyz_yyyzz[i] * fe_0 + tr_z_xxxyz_yyyzz[i] * pa_x[i];

        tr_z_xxxxyz_yyzzz[i] = 3.0 * tr_z_xxyz_yyzzz[i] * fe_0 + tr_z_xxxyz_yyzzz[i] * pa_x[i];

        tr_z_xxxxyz_yzzzz[i] = 3.0 * tr_z_xxyz_yzzzz[i] * fe_0 + tr_z_xxxyz_yzzzz[i] * pa_x[i];

        tr_z_xxxxyz_zzzzz[i] = tr_z_xxxxz_zzzzz[i] * pa_y[i];
    }

    // Set up 1281-1302 components of targeted buffer : IH

    auto tr_z_xxxxzz_xxxxx = pbuffer.data(idx_dip_ih + 1281);

    auto tr_z_xxxxzz_xxxxy = pbuffer.data(idx_dip_ih + 1282);

    auto tr_z_xxxxzz_xxxxz = pbuffer.data(idx_dip_ih + 1283);

    auto tr_z_xxxxzz_xxxyy = pbuffer.data(idx_dip_ih + 1284);

    auto tr_z_xxxxzz_xxxyz = pbuffer.data(idx_dip_ih + 1285);

    auto tr_z_xxxxzz_xxxzz = pbuffer.data(idx_dip_ih + 1286);

    auto tr_z_xxxxzz_xxyyy = pbuffer.data(idx_dip_ih + 1287);

    auto tr_z_xxxxzz_xxyyz = pbuffer.data(idx_dip_ih + 1288);

    auto tr_z_xxxxzz_xxyzz = pbuffer.data(idx_dip_ih + 1289);

    auto tr_z_xxxxzz_xxzzz = pbuffer.data(idx_dip_ih + 1290);

    auto tr_z_xxxxzz_xyyyy = pbuffer.data(idx_dip_ih + 1291);

    auto tr_z_xxxxzz_xyyyz = pbuffer.data(idx_dip_ih + 1292);

    auto tr_z_xxxxzz_xyyzz = pbuffer.data(idx_dip_ih + 1293);

    auto tr_z_xxxxzz_xyzzz = pbuffer.data(idx_dip_ih + 1294);

    auto tr_z_xxxxzz_xzzzz = pbuffer.data(idx_dip_ih + 1295);

    auto tr_z_xxxxzz_yyyyy = pbuffer.data(idx_dip_ih + 1296);

    auto tr_z_xxxxzz_yyyyz = pbuffer.data(idx_dip_ih + 1297);

    auto tr_z_xxxxzz_yyyzz = pbuffer.data(idx_dip_ih + 1298);

    auto tr_z_xxxxzz_yyzzz = pbuffer.data(idx_dip_ih + 1299);

    auto tr_z_xxxxzz_yzzzz = pbuffer.data(idx_dip_ih + 1300);

    auto tr_z_xxxxzz_zzzzz = pbuffer.data(idx_dip_ih + 1301);

    #pragma omp simd aligned(pa_x, tr_z_xxxxzz_xxxxx, tr_z_xxxxzz_xxxxy, tr_z_xxxxzz_xxxxz, tr_z_xxxxzz_xxxyy, tr_z_xxxxzz_xxxyz, tr_z_xxxxzz_xxxzz, tr_z_xxxxzz_xxyyy, tr_z_xxxxzz_xxyyz, tr_z_xxxxzz_xxyzz, tr_z_xxxxzz_xxzzz, tr_z_xxxxzz_xyyyy, tr_z_xxxxzz_xyyyz, tr_z_xxxxzz_xyyzz, tr_z_xxxxzz_xyzzz, tr_z_xxxxzz_xzzzz, tr_z_xxxxzz_yyyyy, tr_z_xxxxzz_yyyyz, tr_z_xxxxzz_yyyzz, tr_z_xxxxzz_yyzzz, tr_z_xxxxzz_yzzzz, tr_z_xxxxzz_zzzzz, tr_z_xxxzz_xxxx, tr_z_xxxzz_xxxxx, tr_z_xxxzz_xxxxy, tr_z_xxxzz_xxxxz, tr_z_xxxzz_xxxy, tr_z_xxxzz_xxxyy, tr_z_xxxzz_xxxyz, tr_z_xxxzz_xxxz, tr_z_xxxzz_xxxzz, tr_z_xxxzz_xxyy, tr_z_xxxzz_xxyyy, tr_z_xxxzz_xxyyz, tr_z_xxxzz_xxyz, tr_z_xxxzz_xxyzz, tr_z_xxxzz_xxzz, tr_z_xxxzz_xxzzz, tr_z_xxxzz_xyyy, tr_z_xxxzz_xyyyy, tr_z_xxxzz_xyyyz, tr_z_xxxzz_xyyz, tr_z_xxxzz_xyyzz, tr_z_xxxzz_xyzz, tr_z_xxxzz_xyzzz, tr_z_xxxzz_xzzz, tr_z_xxxzz_xzzzz, tr_z_xxxzz_yyyy, tr_z_xxxzz_yyyyy, tr_z_xxxzz_yyyyz, tr_z_xxxzz_yyyz, tr_z_xxxzz_yyyzz, tr_z_xxxzz_yyzz, tr_z_xxxzz_yyzzz, tr_z_xxxzz_yzzz, tr_z_xxxzz_yzzzz, tr_z_xxxzz_zzzz, tr_z_xxxzz_zzzzz, tr_z_xxzz_xxxxx, tr_z_xxzz_xxxxy, tr_z_xxzz_xxxxz, tr_z_xxzz_xxxyy, tr_z_xxzz_xxxyz, tr_z_xxzz_xxxzz, tr_z_xxzz_xxyyy, tr_z_xxzz_xxyyz, tr_z_xxzz_xxyzz, tr_z_xxzz_xxzzz, tr_z_xxzz_xyyyy, tr_z_xxzz_xyyyz, tr_z_xxzz_xyyzz, tr_z_xxzz_xyzzz, tr_z_xxzz_xzzzz, tr_z_xxzz_yyyyy, tr_z_xxzz_yyyyz, tr_z_xxzz_yyyzz, tr_z_xxzz_yyzzz, tr_z_xxzz_yzzzz, tr_z_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxzz_xxxxx[i] = 3.0 * tr_z_xxzz_xxxxx[i] * fe_0 + 5.0 * tr_z_xxxzz_xxxx[i] * fe_0 + tr_z_xxxzz_xxxxx[i] * pa_x[i];

        tr_z_xxxxzz_xxxxy[i] = 3.0 * tr_z_xxzz_xxxxy[i] * fe_0 + 4.0 * tr_z_xxxzz_xxxy[i] * fe_0 + tr_z_xxxzz_xxxxy[i] * pa_x[i];

        tr_z_xxxxzz_xxxxz[i] = 3.0 * tr_z_xxzz_xxxxz[i] * fe_0 + 4.0 * tr_z_xxxzz_xxxz[i] * fe_0 + tr_z_xxxzz_xxxxz[i] * pa_x[i];

        tr_z_xxxxzz_xxxyy[i] = 3.0 * tr_z_xxzz_xxxyy[i] * fe_0 + 3.0 * tr_z_xxxzz_xxyy[i] * fe_0 + tr_z_xxxzz_xxxyy[i] * pa_x[i];

        tr_z_xxxxzz_xxxyz[i] = 3.0 * tr_z_xxzz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxxzz_xxyz[i] * fe_0 + tr_z_xxxzz_xxxyz[i] * pa_x[i];

        tr_z_xxxxzz_xxxzz[i] = 3.0 * tr_z_xxzz_xxxzz[i] * fe_0 + 3.0 * tr_z_xxxzz_xxzz[i] * fe_0 + tr_z_xxxzz_xxxzz[i] * pa_x[i];

        tr_z_xxxxzz_xxyyy[i] = 3.0 * tr_z_xxzz_xxyyy[i] * fe_0 + 2.0 * tr_z_xxxzz_xyyy[i] * fe_0 + tr_z_xxxzz_xxyyy[i] * pa_x[i];

        tr_z_xxxxzz_xxyyz[i] = 3.0 * tr_z_xxzz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxxzz_xyyz[i] * fe_0 + tr_z_xxxzz_xxyyz[i] * pa_x[i];

        tr_z_xxxxzz_xxyzz[i] = 3.0 * tr_z_xxzz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxxzz_xyzz[i] * fe_0 + tr_z_xxxzz_xxyzz[i] * pa_x[i];

        tr_z_xxxxzz_xxzzz[i] = 3.0 * tr_z_xxzz_xxzzz[i] * fe_0 + 2.0 * tr_z_xxxzz_xzzz[i] * fe_0 + tr_z_xxxzz_xxzzz[i] * pa_x[i];

        tr_z_xxxxzz_xyyyy[i] = 3.0 * tr_z_xxzz_xyyyy[i] * fe_0 + tr_z_xxxzz_yyyy[i] * fe_0 + tr_z_xxxzz_xyyyy[i] * pa_x[i];

        tr_z_xxxxzz_xyyyz[i] = 3.0 * tr_z_xxzz_xyyyz[i] * fe_0 + tr_z_xxxzz_yyyz[i] * fe_0 + tr_z_xxxzz_xyyyz[i] * pa_x[i];

        tr_z_xxxxzz_xyyzz[i] = 3.0 * tr_z_xxzz_xyyzz[i] * fe_0 + tr_z_xxxzz_yyzz[i] * fe_0 + tr_z_xxxzz_xyyzz[i] * pa_x[i];

        tr_z_xxxxzz_xyzzz[i] = 3.0 * tr_z_xxzz_xyzzz[i] * fe_0 + tr_z_xxxzz_yzzz[i] * fe_0 + tr_z_xxxzz_xyzzz[i] * pa_x[i];

        tr_z_xxxxzz_xzzzz[i] = 3.0 * tr_z_xxzz_xzzzz[i] * fe_0 + tr_z_xxxzz_zzzz[i] * fe_0 + tr_z_xxxzz_xzzzz[i] * pa_x[i];

        tr_z_xxxxzz_yyyyy[i] = 3.0 * tr_z_xxzz_yyyyy[i] * fe_0 + tr_z_xxxzz_yyyyy[i] * pa_x[i];

        tr_z_xxxxzz_yyyyz[i] = 3.0 * tr_z_xxzz_yyyyz[i] * fe_0 + tr_z_xxxzz_yyyyz[i] * pa_x[i];

        tr_z_xxxxzz_yyyzz[i] = 3.0 * tr_z_xxzz_yyyzz[i] * fe_0 + tr_z_xxxzz_yyyzz[i] * pa_x[i];

        tr_z_xxxxzz_yyzzz[i] = 3.0 * tr_z_xxzz_yyzzz[i] * fe_0 + tr_z_xxxzz_yyzzz[i] * pa_x[i];

        tr_z_xxxxzz_yzzzz[i] = 3.0 * tr_z_xxzz_yzzzz[i] * fe_0 + tr_z_xxxzz_yzzzz[i] * pa_x[i];

        tr_z_xxxxzz_zzzzz[i] = 3.0 * tr_z_xxzz_zzzzz[i] * fe_0 + tr_z_xxxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1302-1323 components of targeted buffer : IH

    auto tr_z_xxxyyy_xxxxx = pbuffer.data(idx_dip_ih + 1302);

    auto tr_z_xxxyyy_xxxxy = pbuffer.data(idx_dip_ih + 1303);

    auto tr_z_xxxyyy_xxxxz = pbuffer.data(idx_dip_ih + 1304);

    auto tr_z_xxxyyy_xxxyy = pbuffer.data(idx_dip_ih + 1305);

    auto tr_z_xxxyyy_xxxyz = pbuffer.data(idx_dip_ih + 1306);

    auto tr_z_xxxyyy_xxxzz = pbuffer.data(idx_dip_ih + 1307);

    auto tr_z_xxxyyy_xxyyy = pbuffer.data(idx_dip_ih + 1308);

    auto tr_z_xxxyyy_xxyyz = pbuffer.data(idx_dip_ih + 1309);

    auto tr_z_xxxyyy_xxyzz = pbuffer.data(idx_dip_ih + 1310);

    auto tr_z_xxxyyy_xxzzz = pbuffer.data(idx_dip_ih + 1311);

    auto tr_z_xxxyyy_xyyyy = pbuffer.data(idx_dip_ih + 1312);

    auto tr_z_xxxyyy_xyyyz = pbuffer.data(idx_dip_ih + 1313);

    auto tr_z_xxxyyy_xyyzz = pbuffer.data(idx_dip_ih + 1314);

    auto tr_z_xxxyyy_xyzzz = pbuffer.data(idx_dip_ih + 1315);

    auto tr_z_xxxyyy_xzzzz = pbuffer.data(idx_dip_ih + 1316);

    auto tr_z_xxxyyy_yyyyy = pbuffer.data(idx_dip_ih + 1317);

    auto tr_z_xxxyyy_yyyyz = pbuffer.data(idx_dip_ih + 1318);

    auto tr_z_xxxyyy_yyyzz = pbuffer.data(idx_dip_ih + 1319);

    auto tr_z_xxxyyy_yyzzz = pbuffer.data(idx_dip_ih + 1320);

    auto tr_z_xxxyyy_yzzzz = pbuffer.data(idx_dip_ih + 1321);

    auto tr_z_xxxyyy_zzzzz = pbuffer.data(idx_dip_ih + 1322);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxy_xxxxx, tr_z_xxxy_xxxxz, tr_z_xxxy_xxxzz, tr_z_xxxy_xxzzz, tr_z_xxxy_xzzzz, tr_z_xxxyy_xxxxx, tr_z_xxxyy_xxxxz, tr_z_xxxyy_xxxzz, tr_z_xxxyy_xxzzz, tr_z_xxxyy_xzzzz, tr_z_xxxyyy_xxxxx, tr_z_xxxyyy_xxxxy, tr_z_xxxyyy_xxxxz, tr_z_xxxyyy_xxxyy, tr_z_xxxyyy_xxxyz, tr_z_xxxyyy_xxxzz, tr_z_xxxyyy_xxyyy, tr_z_xxxyyy_xxyyz, tr_z_xxxyyy_xxyzz, tr_z_xxxyyy_xxzzz, tr_z_xxxyyy_xyyyy, tr_z_xxxyyy_xyyyz, tr_z_xxxyyy_xyyzz, tr_z_xxxyyy_xyzzz, tr_z_xxxyyy_xzzzz, tr_z_xxxyyy_yyyyy, tr_z_xxxyyy_yyyyz, tr_z_xxxyyy_yyyzz, tr_z_xxxyyy_yyzzz, tr_z_xxxyyy_yzzzz, tr_z_xxxyyy_zzzzz, tr_z_xxyyy_xxxxy, tr_z_xxyyy_xxxy, tr_z_xxyyy_xxxyy, tr_z_xxyyy_xxxyz, tr_z_xxyyy_xxyy, tr_z_xxyyy_xxyyy, tr_z_xxyyy_xxyyz, tr_z_xxyyy_xxyz, tr_z_xxyyy_xxyzz, tr_z_xxyyy_xyyy, tr_z_xxyyy_xyyyy, tr_z_xxyyy_xyyyz, tr_z_xxyyy_xyyz, tr_z_xxyyy_xyyzz, tr_z_xxyyy_xyzz, tr_z_xxyyy_xyzzz, tr_z_xxyyy_yyyy, tr_z_xxyyy_yyyyy, tr_z_xxyyy_yyyyz, tr_z_xxyyy_yyyz, tr_z_xxyyy_yyyzz, tr_z_xxyyy_yyzz, tr_z_xxyyy_yyzzz, tr_z_xxyyy_yzzz, tr_z_xxyyy_yzzzz, tr_z_xxyyy_zzzzz, tr_z_xyyy_xxxxy, tr_z_xyyy_xxxyy, tr_z_xyyy_xxxyz, tr_z_xyyy_xxyyy, tr_z_xyyy_xxyyz, tr_z_xyyy_xxyzz, tr_z_xyyy_xyyyy, tr_z_xyyy_xyyyz, tr_z_xyyy_xyyzz, tr_z_xyyy_xyzzz, tr_z_xyyy_yyyyy, tr_z_xyyy_yyyyz, tr_z_xyyy_yyyzz, tr_z_xyyy_yyzzz, tr_z_xyyy_yzzzz, tr_z_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyy_xxxxx[i] = 2.0 * tr_z_xxxy_xxxxx[i] * fe_0 + tr_z_xxxyy_xxxxx[i] * pa_y[i];

        tr_z_xxxyyy_xxxxy[i] = 2.0 * tr_z_xyyy_xxxxy[i] * fe_0 + 4.0 * tr_z_xxyyy_xxxy[i] * fe_0 + tr_z_xxyyy_xxxxy[i] * pa_x[i];

        tr_z_xxxyyy_xxxxz[i] = 2.0 * tr_z_xxxy_xxxxz[i] * fe_0 + tr_z_xxxyy_xxxxz[i] * pa_y[i];

        tr_z_xxxyyy_xxxyy[i] = 2.0 * tr_z_xyyy_xxxyy[i] * fe_0 + 3.0 * tr_z_xxyyy_xxyy[i] * fe_0 + tr_z_xxyyy_xxxyy[i] * pa_x[i];

        tr_z_xxxyyy_xxxyz[i] = 2.0 * tr_z_xyyy_xxxyz[i] * fe_0 + 3.0 * tr_z_xxyyy_xxyz[i] * fe_0 + tr_z_xxyyy_xxxyz[i] * pa_x[i];

        tr_z_xxxyyy_xxxzz[i] = 2.0 * tr_z_xxxy_xxxzz[i] * fe_0 + tr_z_xxxyy_xxxzz[i] * pa_y[i];

        tr_z_xxxyyy_xxyyy[i] = 2.0 * tr_z_xyyy_xxyyy[i] * fe_0 + 2.0 * tr_z_xxyyy_xyyy[i] * fe_0 + tr_z_xxyyy_xxyyy[i] * pa_x[i];

        tr_z_xxxyyy_xxyyz[i] = 2.0 * tr_z_xyyy_xxyyz[i] * fe_0 + 2.0 * tr_z_xxyyy_xyyz[i] * fe_0 + tr_z_xxyyy_xxyyz[i] * pa_x[i];

        tr_z_xxxyyy_xxyzz[i] = 2.0 * tr_z_xyyy_xxyzz[i] * fe_0 + 2.0 * tr_z_xxyyy_xyzz[i] * fe_0 + tr_z_xxyyy_xxyzz[i] * pa_x[i];

        tr_z_xxxyyy_xxzzz[i] = 2.0 * tr_z_xxxy_xxzzz[i] * fe_0 + tr_z_xxxyy_xxzzz[i] * pa_y[i];

        tr_z_xxxyyy_xyyyy[i] = 2.0 * tr_z_xyyy_xyyyy[i] * fe_0 + tr_z_xxyyy_yyyy[i] * fe_0 + tr_z_xxyyy_xyyyy[i] * pa_x[i];

        tr_z_xxxyyy_xyyyz[i] = 2.0 * tr_z_xyyy_xyyyz[i] * fe_0 + tr_z_xxyyy_yyyz[i] * fe_0 + tr_z_xxyyy_xyyyz[i] * pa_x[i];

        tr_z_xxxyyy_xyyzz[i] = 2.0 * tr_z_xyyy_xyyzz[i] * fe_0 + tr_z_xxyyy_yyzz[i] * fe_0 + tr_z_xxyyy_xyyzz[i] * pa_x[i];

        tr_z_xxxyyy_xyzzz[i] = 2.0 * tr_z_xyyy_xyzzz[i] * fe_0 + tr_z_xxyyy_yzzz[i] * fe_0 + tr_z_xxyyy_xyzzz[i] * pa_x[i];

        tr_z_xxxyyy_xzzzz[i] = 2.0 * tr_z_xxxy_xzzzz[i] * fe_0 + tr_z_xxxyy_xzzzz[i] * pa_y[i];

        tr_z_xxxyyy_yyyyy[i] = 2.0 * tr_z_xyyy_yyyyy[i] * fe_0 + tr_z_xxyyy_yyyyy[i] * pa_x[i];

        tr_z_xxxyyy_yyyyz[i] = 2.0 * tr_z_xyyy_yyyyz[i] * fe_0 + tr_z_xxyyy_yyyyz[i] * pa_x[i];

        tr_z_xxxyyy_yyyzz[i] = 2.0 * tr_z_xyyy_yyyzz[i] * fe_0 + tr_z_xxyyy_yyyzz[i] * pa_x[i];

        tr_z_xxxyyy_yyzzz[i] = 2.0 * tr_z_xyyy_yyzzz[i] * fe_0 + tr_z_xxyyy_yyzzz[i] * pa_x[i];

        tr_z_xxxyyy_yzzzz[i] = 2.0 * tr_z_xyyy_yzzzz[i] * fe_0 + tr_z_xxyyy_yzzzz[i] * pa_x[i];

        tr_z_xxxyyy_zzzzz[i] = 2.0 * tr_z_xyyy_zzzzz[i] * fe_0 + tr_z_xxyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 1323-1344 components of targeted buffer : IH

    auto tr_z_xxxyyz_xxxxx = pbuffer.data(idx_dip_ih + 1323);

    auto tr_z_xxxyyz_xxxxy = pbuffer.data(idx_dip_ih + 1324);

    auto tr_z_xxxyyz_xxxxz = pbuffer.data(idx_dip_ih + 1325);

    auto tr_z_xxxyyz_xxxyy = pbuffer.data(idx_dip_ih + 1326);

    auto tr_z_xxxyyz_xxxyz = pbuffer.data(idx_dip_ih + 1327);

    auto tr_z_xxxyyz_xxxzz = pbuffer.data(idx_dip_ih + 1328);

    auto tr_z_xxxyyz_xxyyy = pbuffer.data(idx_dip_ih + 1329);

    auto tr_z_xxxyyz_xxyyz = pbuffer.data(idx_dip_ih + 1330);

    auto tr_z_xxxyyz_xxyzz = pbuffer.data(idx_dip_ih + 1331);

    auto tr_z_xxxyyz_xxzzz = pbuffer.data(idx_dip_ih + 1332);

    auto tr_z_xxxyyz_xyyyy = pbuffer.data(idx_dip_ih + 1333);

    auto tr_z_xxxyyz_xyyyz = pbuffer.data(idx_dip_ih + 1334);

    auto tr_z_xxxyyz_xyyzz = pbuffer.data(idx_dip_ih + 1335);

    auto tr_z_xxxyyz_xyzzz = pbuffer.data(idx_dip_ih + 1336);

    auto tr_z_xxxyyz_xzzzz = pbuffer.data(idx_dip_ih + 1337);

    auto tr_z_xxxyyz_yyyyy = pbuffer.data(idx_dip_ih + 1338);

    auto tr_z_xxxyyz_yyyyz = pbuffer.data(idx_dip_ih + 1339);

    auto tr_z_xxxyyz_yyyzz = pbuffer.data(idx_dip_ih + 1340);

    auto tr_z_xxxyyz_yyzzz = pbuffer.data(idx_dip_ih + 1341);

    auto tr_z_xxxyyz_yzzzz = pbuffer.data(idx_dip_ih + 1342);

    auto tr_z_xxxyyz_zzzzz = pbuffer.data(idx_dip_ih + 1343);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxxyy_xxxxy, tr_z_xxxyy_xxxyy, tr_z_xxxyy_xxyyy, tr_z_xxxyy_xyyyy, tr_z_xxxyyz_xxxxx, tr_z_xxxyyz_xxxxy, tr_z_xxxyyz_xxxxz, tr_z_xxxyyz_xxxyy, tr_z_xxxyyz_xxxyz, tr_z_xxxyyz_xxxzz, tr_z_xxxyyz_xxyyy, tr_z_xxxyyz_xxyyz, tr_z_xxxyyz_xxyzz, tr_z_xxxyyz_xxzzz, tr_z_xxxyyz_xyyyy, tr_z_xxxyyz_xyyyz, tr_z_xxxyyz_xyyzz, tr_z_xxxyyz_xyzzz, tr_z_xxxyyz_xzzzz, tr_z_xxxyyz_yyyyy, tr_z_xxxyyz_yyyyz, tr_z_xxxyyz_yyyzz, tr_z_xxxyyz_yyzzz, tr_z_xxxyyz_yzzzz, tr_z_xxxyyz_zzzzz, tr_z_xxxyz_xxxxx, tr_z_xxxyz_xxxxz, tr_z_xxxyz_xxxzz, tr_z_xxxyz_xxzzz, tr_z_xxxyz_xzzzz, tr_z_xxxz_xxxxx, tr_z_xxxz_xxxxz, tr_z_xxxz_xxxzz, tr_z_xxxz_xxzzz, tr_z_xxxz_xzzzz, tr_z_xxyyz_xxxyz, tr_z_xxyyz_xxyyz, tr_z_xxyyz_xxyz, tr_z_xxyyz_xxyzz, tr_z_xxyyz_xyyyz, tr_z_xxyyz_xyyz, tr_z_xxyyz_xyyzz, tr_z_xxyyz_xyzz, tr_z_xxyyz_xyzzz, tr_z_xxyyz_yyyyy, tr_z_xxyyz_yyyyz, tr_z_xxyyz_yyyz, tr_z_xxyyz_yyyzz, tr_z_xxyyz_yyzz, tr_z_xxyyz_yyzzz, tr_z_xxyyz_yzzz, tr_z_xxyyz_yzzzz, tr_z_xxyyz_zzzzz, tr_z_xyyz_xxxyz, tr_z_xyyz_xxyyz, tr_z_xyyz_xxyzz, tr_z_xyyz_xyyyz, tr_z_xyyz_xyyzz, tr_z_xyyz_xyzzz, tr_z_xyyz_yyyyy, tr_z_xyyz_yyyyz, tr_z_xyyz_yyyzz, tr_z_xyyz_yyzzz, tr_z_xyyz_yzzzz, tr_z_xyyz_zzzzz, ts_xxxyy_xxxxy, ts_xxxyy_xxxyy, ts_xxxyy_xxyyy, ts_xxxyy_xyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyz_xxxxx[i] = tr_z_xxxz_xxxxx[i] * fe_0 + tr_z_xxxyz_xxxxx[i] * pa_y[i];

        tr_z_xxxyyz_xxxxy[i] = ts_xxxyy_xxxxy[i] * fe_0 + tr_z_xxxyy_xxxxy[i] * pa_z[i];

        tr_z_xxxyyz_xxxxz[i] = tr_z_xxxz_xxxxz[i] * fe_0 + tr_z_xxxyz_xxxxz[i] * pa_y[i];

        tr_z_xxxyyz_xxxyy[i] = ts_xxxyy_xxxyy[i] * fe_0 + tr_z_xxxyy_xxxyy[i] * pa_z[i];

        tr_z_xxxyyz_xxxyz[i] = 2.0 * tr_z_xyyz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxyyz_xxyz[i] * fe_0 + tr_z_xxyyz_xxxyz[i] * pa_x[i];

        tr_z_xxxyyz_xxxzz[i] = tr_z_xxxz_xxxzz[i] * fe_0 + tr_z_xxxyz_xxxzz[i] * pa_y[i];

        tr_z_xxxyyz_xxyyy[i] = ts_xxxyy_xxyyy[i] * fe_0 + tr_z_xxxyy_xxyyy[i] * pa_z[i];

        tr_z_xxxyyz_xxyyz[i] = 2.0 * tr_z_xyyz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxyyz_xyyz[i] * fe_0 + tr_z_xxyyz_xxyyz[i] * pa_x[i];

        tr_z_xxxyyz_xxyzz[i] = 2.0 * tr_z_xyyz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxyyz_xyzz[i] * fe_0 + tr_z_xxyyz_xxyzz[i] * pa_x[i];

        tr_z_xxxyyz_xxzzz[i] = tr_z_xxxz_xxzzz[i] * fe_0 + tr_z_xxxyz_xxzzz[i] * pa_y[i];

        tr_z_xxxyyz_xyyyy[i] = ts_xxxyy_xyyyy[i] * fe_0 + tr_z_xxxyy_xyyyy[i] * pa_z[i];

        tr_z_xxxyyz_xyyyz[i] = 2.0 * tr_z_xyyz_xyyyz[i] * fe_0 + tr_z_xxyyz_yyyz[i] * fe_0 + tr_z_xxyyz_xyyyz[i] * pa_x[i];

        tr_z_xxxyyz_xyyzz[i] = 2.0 * tr_z_xyyz_xyyzz[i] * fe_0 + tr_z_xxyyz_yyzz[i] * fe_0 + tr_z_xxyyz_xyyzz[i] * pa_x[i];

        tr_z_xxxyyz_xyzzz[i] = 2.0 * tr_z_xyyz_xyzzz[i] * fe_0 + tr_z_xxyyz_yzzz[i] * fe_0 + tr_z_xxyyz_xyzzz[i] * pa_x[i];

        tr_z_xxxyyz_xzzzz[i] = tr_z_xxxz_xzzzz[i] * fe_0 + tr_z_xxxyz_xzzzz[i] * pa_y[i];

        tr_z_xxxyyz_yyyyy[i] = 2.0 * tr_z_xyyz_yyyyy[i] * fe_0 + tr_z_xxyyz_yyyyy[i] * pa_x[i];

        tr_z_xxxyyz_yyyyz[i] = 2.0 * tr_z_xyyz_yyyyz[i] * fe_0 + tr_z_xxyyz_yyyyz[i] * pa_x[i];

        tr_z_xxxyyz_yyyzz[i] = 2.0 * tr_z_xyyz_yyyzz[i] * fe_0 + tr_z_xxyyz_yyyzz[i] * pa_x[i];

        tr_z_xxxyyz_yyzzz[i] = 2.0 * tr_z_xyyz_yyzzz[i] * fe_0 + tr_z_xxyyz_yyzzz[i] * pa_x[i];

        tr_z_xxxyyz_yzzzz[i] = 2.0 * tr_z_xyyz_yzzzz[i] * fe_0 + tr_z_xxyyz_yzzzz[i] * pa_x[i];

        tr_z_xxxyyz_zzzzz[i] = 2.0 * tr_z_xyyz_zzzzz[i] * fe_0 + tr_z_xxyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 1344-1365 components of targeted buffer : IH

    auto tr_z_xxxyzz_xxxxx = pbuffer.data(idx_dip_ih + 1344);

    auto tr_z_xxxyzz_xxxxy = pbuffer.data(idx_dip_ih + 1345);

    auto tr_z_xxxyzz_xxxxz = pbuffer.data(idx_dip_ih + 1346);

    auto tr_z_xxxyzz_xxxyy = pbuffer.data(idx_dip_ih + 1347);

    auto tr_z_xxxyzz_xxxyz = pbuffer.data(idx_dip_ih + 1348);

    auto tr_z_xxxyzz_xxxzz = pbuffer.data(idx_dip_ih + 1349);

    auto tr_z_xxxyzz_xxyyy = pbuffer.data(idx_dip_ih + 1350);

    auto tr_z_xxxyzz_xxyyz = pbuffer.data(idx_dip_ih + 1351);

    auto tr_z_xxxyzz_xxyzz = pbuffer.data(idx_dip_ih + 1352);

    auto tr_z_xxxyzz_xxzzz = pbuffer.data(idx_dip_ih + 1353);

    auto tr_z_xxxyzz_xyyyy = pbuffer.data(idx_dip_ih + 1354);

    auto tr_z_xxxyzz_xyyyz = pbuffer.data(idx_dip_ih + 1355);

    auto tr_z_xxxyzz_xyyzz = pbuffer.data(idx_dip_ih + 1356);

    auto tr_z_xxxyzz_xyzzz = pbuffer.data(idx_dip_ih + 1357);

    auto tr_z_xxxyzz_xzzzz = pbuffer.data(idx_dip_ih + 1358);

    auto tr_z_xxxyzz_yyyyy = pbuffer.data(idx_dip_ih + 1359);

    auto tr_z_xxxyzz_yyyyz = pbuffer.data(idx_dip_ih + 1360);

    auto tr_z_xxxyzz_yyyzz = pbuffer.data(idx_dip_ih + 1361);

    auto tr_z_xxxyzz_yyzzz = pbuffer.data(idx_dip_ih + 1362);

    auto tr_z_xxxyzz_yzzzz = pbuffer.data(idx_dip_ih + 1363);

    auto tr_z_xxxyzz_zzzzz = pbuffer.data(idx_dip_ih + 1364);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxyzz_xxxxx, tr_z_xxxyzz_xxxxy, tr_z_xxxyzz_xxxxz, tr_z_xxxyzz_xxxyy, tr_z_xxxyzz_xxxyz, tr_z_xxxyzz_xxxzz, tr_z_xxxyzz_xxyyy, tr_z_xxxyzz_xxyyz, tr_z_xxxyzz_xxyzz, tr_z_xxxyzz_xxzzz, tr_z_xxxyzz_xyyyy, tr_z_xxxyzz_xyyyz, tr_z_xxxyzz_xyyzz, tr_z_xxxyzz_xyzzz, tr_z_xxxyzz_xzzzz, tr_z_xxxyzz_yyyyy, tr_z_xxxyzz_yyyyz, tr_z_xxxyzz_yyyzz, tr_z_xxxyzz_yyzzz, tr_z_xxxyzz_yzzzz, tr_z_xxxyzz_zzzzz, tr_z_xxxzz_xxxx, tr_z_xxxzz_xxxxx, tr_z_xxxzz_xxxxy, tr_z_xxxzz_xxxxz, tr_z_xxxzz_xxxy, tr_z_xxxzz_xxxyy, tr_z_xxxzz_xxxyz, tr_z_xxxzz_xxxz, tr_z_xxxzz_xxxzz, tr_z_xxxzz_xxyy, tr_z_xxxzz_xxyyy, tr_z_xxxzz_xxyyz, tr_z_xxxzz_xxyz, tr_z_xxxzz_xxyzz, tr_z_xxxzz_xxzz, tr_z_xxxzz_xxzzz, tr_z_xxxzz_xyyy, tr_z_xxxzz_xyyyy, tr_z_xxxzz_xyyyz, tr_z_xxxzz_xyyz, tr_z_xxxzz_xyyzz, tr_z_xxxzz_xyzz, tr_z_xxxzz_xyzzz, tr_z_xxxzz_xzzz, tr_z_xxxzz_xzzzz, tr_z_xxxzz_zzzzz, tr_z_xxyzz_yyyyy, tr_z_xxyzz_yyyyz, tr_z_xxyzz_yyyzz, tr_z_xxyzz_yyzzz, tr_z_xxyzz_yzzzz, tr_z_xyzz_yyyyy, tr_z_xyzz_yyyyz, tr_z_xyzz_yyyzz, tr_z_xyzz_yyzzz, tr_z_xyzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyzz_xxxxx[i] = tr_z_xxxzz_xxxxx[i] * pa_y[i];

        tr_z_xxxyzz_xxxxy[i] = tr_z_xxxzz_xxxx[i] * fe_0 + tr_z_xxxzz_xxxxy[i] * pa_y[i];

        tr_z_xxxyzz_xxxxz[i] = tr_z_xxxzz_xxxxz[i] * pa_y[i];

        tr_z_xxxyzz_xxxyy[i] = 2.0 * tr_z_xxxzz_xxxy[i] * fe_0 + tr_z_xxxzz_xxxyy[i] * pa_y[i];

        tr_z_xxxyzz_xxxyz[i] = tr_z_xxxzz_xxxz[i] * fe_0 + tr_z_xxxzz_xxxyz[i] * pa_y[i];

        tr_z_xxxyzz_xxxzz[i] = tr_z_xxxzz_xxxzz[i] * pa_y[i];

        tr_z_xxxyzz_xxyyy[i] = 3.0 * tr_z_xxxzz_xxyy[i] * fe_0 + tr_z_xxxzz_xxyyy[i] * pa_y[i];

        tr_z_xxxyzz_xxyyz[i] = 2.0 * tr_z_xxxzz_xxyz[i] * fe_0 + tr_z_xxxzz_xxyyz[i] * pa_y[i];

        tr_z_xxxyzz_xxyzz[i] = tr_z_xxxzz_xxzz[i] * fe_0 + tr_z_xxxzz_xxyzz[i] * pa_y[i];

        tr_z_xxxyzz_xxzzz[i] = tr_z_xxxzz_xxzzz[i] * pa_y[i];

        tr_z_xxxyzz_xyyyy[i] = 4.0 * tr_z_xxxzz_xyyy[i] * fe_0 + tr_z_xxxzz_xyyyy[i] * pa_y[i];

        tr_z_xxxyzz_xyyyz[i] = 3.0 * tr_z_xxxzz_xyyz[i] * fe_0 + tr_z_xxxzz_xyyyz[i] * pa_y[i];

        tr_z_xxxyzz_xyyzz[i] = 2.0 * tr_z_xxxzz_xyzz[i] * fe_0 + tr_z_xxxzz_xyyzz[i] * pa_y[i];

        tr_z_xxxyzz_xyzzz[i] = tr_z_xxxzz_xzzz[i] * fe_0 + tr_z_xxxzz_xyzzz[i] * pa_y[i];

        tr_z_xxxyzz_xzzzz[i] = tr_z_xxxzz_xzzzz[i] * pa_y[i];

        tr_z_xxxyzz_yyyyy[i] = 2.0 * tr_z_xyzz_yyyyy[i] * fe_0 + tr_z_xxyzz_yyyyy[i] * pa_x[i];

        tr_z_xxxyzz_yyyyz[i] = 2.0 * tr_z_xyzz_yyyyz[i] * fe_0 + tr_z_xxyzz_yyyyz[i] * pa_x[i];

        tr_z_xxxyzz_yyyzz[i] = 2.0 * tr_z_xyzz_yyyzz[i] * fe_0 + tr_z_xxyzz_yyyzz[i] * pa_x[i];

        tr_z_xxxyzz_yyzzz[i] = 2.0 * tr_z_xyzz_yyzzz[i] * fe_0 + tr_z_xxyzz_yyzzz[i] * pa_x[i];

        tr_z_xxxyzz_yzzzz[i] = 2.0 * tr_z_xyzz_yzzzz[i] * fe_0 + tr_z_xxyzz_yzzzz[i] * pa_x[i];

        tr_z_xxxyzz_zzzzz[i] = tr_z_xxxzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1365-1386 components of targeted buffer : IH

    auto tr_z_xxxzzz_xxxxx = pbuffer.data(idx_dip_ih + 1365);

    auto tr_z_xxxzzz_xxxxy = pbuffer.data(idx_dip_ih + 1366);

    auto tr_z_xxxzzz_xxxxz = pbuffer.data(idx_dip_ih + 1367);

    auto tr_z_xxxzzz_xxxyy = pbuffer.data(idx_dip_ih + 1368);

    auto tr_z_xxxzzz_xxxyz = pbuffer.data(idx_dip_ih + 1369);

    auto tr_z_xxxzzz_xxxzz = pbuffer.data(idx_dip_ih + 1370);

    auto tr_z_xxxzzz_xxyyy = pbuffer.data(idx_dip_ih + 1371);

    auto tr_z_xxxzzz_xxyyz = pbuffer.data(idx_dip_ih + 1372);

    auto tr_z_xxxzzz_xxyzz = pbuffer.data(idx_dip_ih + 1373);

    auto tr_z_xxxzzz_xxzzz = pbuffer.data(idx_dip_ih + 1374);

    auto tr_z_xxxzzz_xyyyy = pbuffer.data(idx_dip_ih + 1375);

    auto tr_z_xxxzzz_xyyyz = pbuffer.data(idx_dip_ih + 1376);

    auto tr_z_xxxzzz_xyyzz = pbuffer.data(idx_dip_ih + 1377);

    auto tr_z_xxxzzz_xyzzz = pbuffer.data(idx_dip_ih + 1378);

    auto tr_z_xxxzzz_xzzzz = pbuffer.data(idx_dip_ih + 1379);

    auto tr_z_xxxzzz_yyyyy = pbuffer.data(idx_dip_ih + 1380);

    auto tr_z_xxxzzz_yyyyz = pbuffer.data(idx_dip_ih + 1381);

    auto tr_z_xxxzzz_yyyzz = pbuffer.data(idx_dip_ih + 1382);

    auto tr_z_xxxzzz_yyzzz = pbuffer.data(idx_dip_ih + 1383);

    auto tr_z_xxxzzz_yzzzz = pbuffer.data(idx_dip_ih + 1384);

    auto tr_z_xxxzzz_zzzzz = pbuffer.data(idx_dip_ih + 1385);

    #pragma omp simd aligned(pa_x, tr_z_xxxzzz_xxxxx, tr_z_xxxzzz_xxxxy, tr_z_xxxzzz_xxxxz, tr_z_xxxzzz_xxxyy, tr_z_xxxzzz_xxxyz, tr_z_xxxzzz_xxxzz, tr_z_xxxzzz_xxyyy, tr_z_xxxzzz_xxyyz, tr_z_xxxzzz_xxyzz, tr_z_xxxzzz_xxzzz, tr_z_xxxzzz_xyyyy, tr_z_xxxzzz_xyyyz, tr_z_xxxzzz_xyyzz, tr_z_xxxzzz_xyzzz, tr_z_xxxzzz_xzzzz, tr_z_xxxzzz_yyyyy, tr_z_xxxzzz_yyyyz, tr_z_xxxzzz_yyyzz, tr_z_xxxzzz_yyzzz, tr_z_xxxzzz_yzzzz, tr_z_xxxzzz_zzzzz, tr_z_xxzzz_xxxx, tr_z_xxzzz_xxxxx, tr_z_xxzzz_xxxxy, tr_z_xxzzz_xxxxz, tr_z_xxzzz_xxxy, tr_z_xxzzz_xxxyy, tr_z_xxzzz_xxxyz, tr_z_xxzzz_xxxz, tr_z_xxzzz_xxxzz, tr_z_xxzzz_xxyy, tr_z_xxzzz_xxyyy, tr_z_xxzzz_xxyyz, tr_z_xxzzz_xxyz, tr_z_xxzzz_xxyzz, tr_z_xxzzz_xxzz, tr_z_xxzzz_xxzzz, tr_z_xxzzz_xyyy, tr_z_xxzzz_xyyyy, tr_z_xxzzz_xyyyz, tr_z_xxzzz_xyyz, tr_z_xxzzz_xyyzz, tr_z_xxzzz_xyzz, tr_z_xxzzz_xyzzz, tr_z_xxzzz_xzzz, tr_z_xxzzz_xzzzz, tr_z_xxzzz_yyyy, tr_z_xxzzz_yyyyy, tr_z_xxzzz_yyyyz, tr_z_xxzzz_yyyz, tr_z_xxzzz_yyyzz, tr_z_xxzzz_yyzz, tr_z_xxzzz_yyzzz, tr_z_xxzzz_yzzz, tr_z_xxzzz_yzzzz, tr_z_xxzzz_zzzz, tr_z_xxzzz_zzzzz, tr_z_xzzz_xxxxx, tr_z_xzzz_xxxxy, tr_z_xzzz_xxxxz, tr_z_xzzz_xxxyy, tr_z_xzzz_xxxyz, tr_z_xzzz_xxxzz, tr_z_xzzz_xxyyy, tr_z_xzzz_xxyyz, tr_z_xzzz_xxyzz, tr_z_xzzz_xxzzz, tr_z_xzzz_xyyyy, tr_z_xzzz_xyyyz, tr_z_xzzz_xyyzz, tr_z_xzzz_xyzzz, tr_z_xzzz_xzzzz, tr_z_xzzz_yyyyy, tr_z_xzzz_yyyyz, tr_z_xzzz_yyyzz, tr_z_xzzz_yyzzz, tr_z_xzzz_yzzzz, tr_z_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzzz_xxxxx[i] = 2.0 * tr_z_xzzz_xxxxx[i] * fe_0 + 5.0 * tr_z_xxzzz_xxxx[i] * fe_0 + tr_z_xxzzz_xxxxx[i] * pa_x[i];

        tr_z_xxxzzz_xxxxy[i] = 2.0 * tr_z_xzzz_xxxxy[i] * fe_0 + 4.0 * tr_z_xxzzz_xxxy[i] * fe_0 + tr_z_xxzzz_xxxxy[i] * pa_x[i];

        tr_z_xxxzzz_xxxxz[i] = 2.0 * tr_z_xzzz_xxxxz[i] * fe_0 + 4.0 * tr_z_xxzzz_xxxz[i] * fe_0 + tr_z_xxzzz_xxxxz[i] * pa_x[i];

        tr_z_xxxzzz_xxxyy[i] = 2.0 * tr_z_xzzz_xxxyy[i] * fe_0 + 3.0 * tr_z_xxzzz_xxyy[i] * fe_0 + tr_z_xxzzz_xxxyy[i] * pa_x[i];

        tr_z_xxxzzz_xxxyz[i] = 2.0 * tr_z_xzzz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxzzz_xxyz[i] * fe_0 + tr_z_xxzzz_xxxyz[i] * pa_x[i];

        tr_z_xxxzzz_xxxzz[i] = 2.0 * tr_z_xzzz_xxxzz[i] * fe_0 + 3.0 * tr_z_xxzzz_xxzz[i] * fe_0 + tr_z_xxzzz_xxxzz[i] * pa_x[i];

        tr_z_xxxzzz_xxyyy[i] = 2.0 * tr_z_xzzz_xxyyy[i] * fe_0 + 2.0 * tr_z_xxzzz_xyyy[i] * fe_0 + tr_z_xxzzz_xxyyy[i] * pa_x[i];

        tr_z_xxxzzz_xxyyz[i] = 2.0 * tr_z_xzzz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxzzz_xyyz[i] * fe_0 + tr_z_xxzzz_xxyyz[i] * pa_x[i];

        tr_z_xxxzzz_xxyzz[i] = 2.0 * tr_z_xzzz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxzzz_xyzz[i] * fe_0 + tr_z_xxzzz_xxyzz[i] * pa_x[i];

        tr_z_xxxzzz_xxzzz[i] = 2.0 * tr_z_xzzz_xxzzz[i] * fe_0 + 2.0 * tr_z_xxzzz_xzzz[i] * fe_0 + tr_z_xxzzz_xxzzz[i] * pa_x[i];

        tr_z_xxxzzz_xyyyy[i] = 2.0 * tr_z_xzzz_xyyyy[i] * fe_0 + tr_z_xxzzz_yyyy[i] * fe_0 + tr_z_xxzzz_xyyyy[i] * pa_x[i];

        tr_z_xxxzzz_xyyyz[i] = 2.0 * tr_z_xzzz_xyyyz[i] * fe_0 + tr_z_xxzzz_yyyz[i] * fe_0 + tr_z_xxzzz_xyyyz[i] * pa_x[i];

        tr_z_xxxzzz_xyyzz[i] = 2.0 * tr_z_xzzz_xyyzz[i] * fe_0 + tr_z_xxzzz_yyzz[i] * fe_0 + tr_z_xxzzz_xyyzz[i] * pa_x[i];

        tr_z_xxxzzz_xyzzz[i] = 2.0 * tr_z_xzzz_xyzzz[i] * fe_0 + tr_z_xxzzz_yzzz[i] * fe_0 + tr_z_xxzzz_xyzzz[i] * pa_x[i];

        tr_z_xxxzzz_xzzzz[i] = 2.0 * tr_z_xzzz_xzzzz[i] * fe_0 + tr_z_xxzzz_zzzz[i] * fe_0 + tr_z_xxzzz_xzzzz[i] * pa_x[i];

        tr_z_xxxzzz_yyyyy[i] = 2.0 * tr_z_xzzz_yyyyy[i] * fe_0 + tr_z_xxzzz_yyyyy[i] * pa_x[i];

        tr_z_xxxzzz_yyyyz[i] = 2.0 * tr_z_xzzz_yyyyz[i] * fe_0 + tr_z_xxzzz_yyyyz[i] * pa_x[i];

        tr_z_xxxzzz_yyyzz[i] = 2.0 * tr_z_xzzz_yyyzz[i] * fe_0 + tr_z_xxzzz_yyyzz[i] * pa_x[i];

        tr_z_xxxzzz_yyzzz[i] = 2.0 * tr_z_xzzz_yyzzz[i] * fe_0 + tr_z_xxzzz_yyzzz[i] * pa_x[i];

        tr_z_xxxzzz_yzzzz[i] = 2.0 * tr_z_xzzz_yzzzz[i] * fe_0 + tr_z_xxzzz_yzzzz[i] * pa_x[i];

        tr_z_xxxzzz_zzzzz[i] = 2.0 * tr_z_xzzz_zzzzz[i] * fe_0 + tr_z_xxzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1386-1407 components of targeted buffer : IH

    auto tr_z_xxyyyy_xxxxx = pbuffer.data(idx_dip_ih + 1386);

    auto tr_z_xxyyyy_xxxxy = pbuffer.data(idx_dip_ih + 1387);

    auto tr_z_xxyyyy_xxxxz = pbuffer.data(idx_dip_ih + 1388);

    auto tr_z_xxyyyy_xxxyy = pbuffer.data(idx_dip_ih + 1389);

    auto tr_z_xxyyyy_xxxyz = pbuffer.data(idx_dip_ih + 1390);

    auto tr_z_xxyyyy_xxxzz = pbuffer.data(idx_dip_ih + 1391);

    auto tr_z_xxyyyy_xxyyy = pbuffer.data(idx_dip_ih + 1392);

    auto tr_z_xxyyyy_xxyyz = pbuffer.data(idx_dip_ih + 1393);

    auto tr_z_xxyyyy_xxyzz = pbuffer.data(idx_dip_ih + 1394);

    auto tr_z_xxyyyy_xxzzz = pbuffer.data(idx_dip_ih + 1395);

    auto tr_z_xxyyyy_xyyyy = pbuffer.data(idx_dip_ih + 1396);

    auto tr_z_xxyyyy_xyyyz = pbuffer.data(idx_dip_ih + 1397);

    auto tr_z_xxyyyy_xyyzz = pbuffer.data(idx_dip_ih + 1398);

    auto tr_z_xxyyyy_xyzzz = pbuffer.data(idx_dip_ih + 1399);

    auto tr_z_xxyyyy_xzzzz = pbuffer.data(idx_dip_ih + 1400);

    auto tr_z_xxyyyy_yyyyy = pbuffer.data(idx_dip_ih + 1401);

    auto tr_z_xxyyyy_yyyyz = pbuffer.data(idx_dip_ih + 1402);

    auto tr_z_xxyyyy_yyyzz = pbuffer.data(idx_dip_ih + 1403);

    auto tr_z_xxyyyy_yyzzz = pbuffer.data(idx_dip_ih + 1404);

    auto tr_z_xxyyyy_yzzzz = pbuffer.data(idx_dip_ih + 1405);

    auto tr_z_xxyyyy_zzzzz = pbuffer.data(idx_dip_ih + 1406);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyy_xxxxx, tr_z_xxyy_xxxxz, tr_z_xxyy_xxxzz, tr_z_xxyy_xxzzz, tr_z_xxyy_xzzzz, tr_z_xxyyy_xxxxx, tr_z_xxyyy_xxxxz, tr_z_xxyyy_xxxzz, tr_z_xxyyy_xxzzz, tr_z_xxyyy_xzzzz, tr_z_xxyyyy_xxxxx, tr_z_xxyyyy_xxxxy, tr_z_xxyyyy_xxxxz, tr_z_xxyyyy_xxxyy, tr_z_xxyyyy_xxxyz, tr_z_xxyyyy_xxxzz, tr_z_xxyyyy_xxyyy, tr_z_xxyyyy_xxyyz, tr_z_xxyyyy_xxyzz, tr_z_xxyyyy_xxzzz, tr_z_xxyyyy_xyyyy, tr_z_xxyyyy_xyyyz, tr_z_xxyyyy_xyyzz, tr_z_xxyyyy_xyzzz, tr_z_xxyyyy_xzzzz, tr_z_xxyyyy_yyyyy, tr_z_xxyyyy_yyyyz, tr_z_xxyyyy_yyyzz, tr_z_xxyyyy_yyzzz, tr_z_xxyyyy_yzzzz, tr_z_xxyyyy_zzzzz, tr_z_xyyyy_xxxxy, tr_z_xyyyy_xxxy, tr_z_xyyyy_xxxyy, tr_z_xyyyy_xxxyz, tr_z_xyyyy_xxyy, tr_z_xyyyy_xxyyy, tr_z_xyyyy_xxyyz, tr_z_xyyyy_xxyz, tr_z_xyyyy_xxyzz, tr_z_xyyyy_xyyy, tr_z_xyyyy_xyyyy, tr_z_xyyyy_xyyyz, tr_z_xyyyy_xyyz, tr_z_xyyyy_xyyzz, tr_z_xyyyy_xyzz, tr_z_xyyyy_xyzzz, tr_z_xyyyy_yyyy, tr_z_xyyyy_yyyyy, tr_z_xyyyy_yyyyz, tr_z_xyyyy_yyyz, tr_z_xyyyy_yyyzz, tr_z_xyyyy_yyzz, tr_z_xyyyy_yyzzz, tr_z_xyyyy_yzzz, tr_z_xyyyy_yzzzz, tr_z_xyyyy_zzzzz, tr_z_yyyy_xxxxy, tr_z_yyyy_xxxyy, tr_z_yyyy_xxxyz, tr_z_yyyy_xxyyy, tr_z_yyyy_xxyyz, tr_z_yyyy_xxyzz, tr_z_yyyy_xyyyy, tr_z_yyyy_xyyyz, tr_z_yyyy_xyyzz, tr_z_yyyy_xyzzz, tr_z_yyyy_yyyyy, tr_z_yyyy_yyyyz, tr_z_yyyy_yyyzz, tr_z_yyyy_yyzzz, tr_z_yyyy_yzzzz, tr_z_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyy_xxxxx[i] = 3.0 * tr_z_xxyy_xxxxx[i] * fe_0 + tr_z_xxyyy_xxxxx[i] * pa_y[i];

        tr_z_xxyyyy_xxxxy[i] = tr_z_yyyy_xxxxy[i] * fe_0 + 4.0 * tr_z_xyyyy_xxxy[i] * fe_0 + tr_z_xyyyy_xxxxy[i] * pa_x[i];

        tr_z_xxyyyy_xxxxz[i] = 3.0 * tr_z_xxyy_xxxxz[i] * fe_0 + tr_z_xxyyy_xxxxz[i] * pa_y[i];

        tr_z_xxyyyy_xxxyy[i] = tr_z_yyyy_xxxyy[i] * fe_0 + 3.0 * tr_z_xyyyy_xxyy[i] * fe_0 + tr_z_xyyyy_xxxyy[i] * pa_x[i];

        tr_z_xxyyyy_xxxyz[i] = tr_z_yyyy_xxxyz[i] * fe_0 + 3.0 * tr_z_xyyyy_xxyz[i] * fe_0 + tr_z_xyyyy_xxxyz[i] * pa_x[i];

        tr_z_xxyyyy_xxxzz[i] = 3.0 * tr_z_xxyy_xxxzz[i] * fe_0 + tr_z_xxyyy_xxxzz[i] * pa_y[i];

        tr_z_xxyyyy_xxyyy[i] = tr_z_yyyy_xxyyy[i] * fe_0 + 2.0 * tr_z_xyyyy_xyyy[i] * fe_0 + tr_z_xyyyy_xxyyy[i] * pa_x[i];

        tr_z_xxyyyy_xxyyz[i] = tr_z_yyyy_xxyyz[i] * fe_0 + 2.0 * tr_z_xyyyy_xyyz[i] * fe_0 + tr_z_xyyyy_xxyyz[i] * pa_x[i];

        tr_z_xxyyyy_xxyzz[i] = tr_z_yyyy_xxyzz[i] * fe_0 + 2.0 * tr_z_xyyyy_xyzz[i] * fe_0 + tr_z_xyyyy_xxyzz[i] * pa_x[i];

        tr_z_xxyyyy_xxzzz[i] = 3.0 * tr_z_xxyy_xxzzz[i] * fe_0 + tr_z_xxyyy_xxzzz[i] * pa_y[i];

        tr_z_xxyyyy_xyyyy[i] = tr_z_yyyy_xyyyy[i] * fe_0 + tr_z_xyyyy_yyyy[i] * fe_0 + tr_z_xyyyy_xyyyy[i] * pa_x[i];

        tr_z_xxyyyy_xyyyz[i] = tr_z_yyyy_xyyyz[i] * fe_0 + tr_z_xyyyy_yyyz[i] * fe_0 + tr_z_xyyyy_xyyyz[i] * pa_x[i];

        tr_z_xxyyyy_xyyzz[i] = tr_z_yyyy_xyyzz[i] * fe_0 + tr_z_xyyyy_yyzz[i] * fe_0 + tr_z_xyyyy_xyyzz[i] * pa_x[i];

        tr_z_xxyyyy_xyzzz[i] = tr_z_yyyy_xyzzz[i] * fe_0 + tr_z_xyyyy_yzzz[i] * fe_0 + tr_z_xyyyy_xyzzz[i] * pa_x[i];

        tr_z_xxyyyy_xzzzz[i] = 3.0 * tr_z_xxyy_xzzzz[i] * fe_0 + tr_z_xxyyy_xzzzz[i] * pa_y[i];

        tr_z_xxyyyy_yyyyy[i] = tr_z_yyyy_yyyyy[i] * fe_0 + tr_z_xyyyy_yyyyy[i] * pa_x[i];

        tr_z_xxyyyy_yyyyz[i] = tr_z_yyyy_yyyyz[i] * fe_0 + tr_z_xyyyy_yyyyz[i] * pa_x[i];

        tr_z_xxyyyy_yyyzz[i] = tr_z_yyyy_yyyzz[i] * fe_0 + tr_z_xyyyy_yyyzz[i] * pa_x[i];

        tr_z_xxyyyy_yyzzz[i] = tr_z_yyyy_yyzzz[i] * fe_0 + tr_z_xyyyy_yyzzz[i] * pa_x[i];

        tr_z_xxyyyy_yzzzz[i] = tr_z_yyyy_yzzzz[i] * fe_0 + tr_z_xyyyy_yzzzz[i] * pa_x[i];

        tr_z_xxyyyy_zzzzz[i] = tr_z_yyyy_zzzzz[i] * fe_0 + tr_z_xyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 1407-1428 components of targeted buffer : IH

    auto tr_z_xxyyyz_xxxxx = pbuffer.data(idx_dip_ih + 1407);

    auto tr_z_xxyyyz_xxxxy = pbuffer.data(idx_dip_ih + 1408);

    auto tr_z_xxyyyz_xxxxz = pbuffer.data(idx_dip_ih + 1409);

    auto tr_z_xxyyyz_xxxyy = pbuffer.data(idx_dip_ih + 1410);

    auto tr_z_xxyyyz_xxxyz = pbuffer.data(idx_dip_ih + 1411);

    auto tr_z_xxyyyz_xxxzz = pbuffer.data(idx_dip_ih + 1412);

    auto tr_z_xxyyyz_xxyyy = pbuffer.data(idx_dip_ih + 1413);

    auto tr_z_xxyyyz_xxyyz = pbuffer.data(idx_dip_ih + 1414);

    auto tr_z_xxyyyz_xxyzz = pbuffer.data(idx_dip_ih + 1415);

    auto tr_z_xxyyyz_xxzzz = pbuffer.data(idx_dip_ih + 1416);

    auto tr_z_xxyyyz_xyyyy = pbuffer.data(idx_dip_ih + 1417);

    auto tr_z_xxyyyz_xyyyz = pbuffer.data(idx_dip_ih + 1418);

    auto tr_z_xxyyyz_xyyzz = pbuffer.data(idx_dip_ih + 1419);

    auto tr_z_xxyyyz_xyzzz = pbuffer.data(idx_dip_ih + 1420);

    auto tr_z_xxyyyz_xzzzz = pbuffer.data(idx_dip_ih + 1421);

    auto tr_z_xxyyyz_yyyyy = pbuffer.data(idx_dip_ih + 1422);

    auto tr_z_xxyyyz_yyyyz = pbuffer.data(idx_dip_ih + 1423);

    auto tr_z_xxyyyz_yyyzz = pbuffer.data(idx_dip_ih + 1424);

    auto tr_z_xxyyyz_yyzzz = pbuffer.data(idx_dip_ih + 1425);

    auto tr_z_xxyyyz_yzzzz = pbuffer.data(idx_dip_ih + 1426);

    auto tr_z_xxyyyz_zzzzz = pbuffer.data(idx_dip_ih + 1427);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxyyy_xxxxy, tr_z_xxyyy_xxxyy, tr_z_xxyyy_xxyyy, tr_z_xxyyy_xyyyy, tr_z_xxyyyz_xxxxx, tr_z_xxyyyz_xxxxy, tr_z_xxyyyz_xxxxz, tr_z_xxyyyz_xxxyy, tr_z_xxyyyz_xxxyz, tr_z_xxyyyz_xxxzz, tr_z_xxyyyz_xxyyy, tr_z_xxyyyz_xxyyz, tr_z_xxyyyz_xxyzz, tr_z_xxyyyz_xxzzz, tr_z_xxyyyz_xyyyy, tr_z_xxyyyz_xyyyz, tr_z_xxyyyz_xyyzz, tr_z_xxyyyz_xyzzz, tr_z_xxyyyz_xzzzz, tr_z_xxyyyz_yyyyy, tr_z_xxyyyz_yyyyz, tr_z_xxyyyz_yyyzz, tr_z_xxyyyz_yyzzz, tr_z_xxyyyz_yzzzz, tr_z_xxyyyz_zzzzz, tr_z_xxyyz_xxxxx, tr_z_xxyyz_xxxxz, tr_z_xxyyz_xxxzz, tr_z_xxyyz_xxzzz, tr_z_xxyyz_xzzzz, tr_z_xxyz_xxxxx, tr_z_xxyz_xxxxz, tr_z_xxyz_xxxzz, tr_z_xxyz_xxzzz, tr_z_xxyz_xzzzz, tr_z_xyyyz_xxxyz, tr_z_xyyyz_xxyyz, tr_z_xyyyz_xxyz, tr_z_xyyyz_xxyzz, tr_z_xyyyz_xyyyz, tr_z_xyyyz_xyyz, tr_z_xyyyz_xyyzz, tr_z_xyyyz_xyzz, tr_z_xyyyz_xyzzz, tr_z_xyyyz_yyyyy, tr_z_xyyyz_yyyyz, tr_z_xyyyz_yyyz, tr_z_xyyyz_yyyzz, tr_z_xyyyz_yyzz, tr_z_xyyyz_yyzzz, tr_z_xyyyz_yzzz, tr_z_xyyyz_yzzzz, tr_z_xyyyz_zzzzz, tr_z_yyyz_xxxyz, tr_z_yyyz_xxyyz, tr_z_yyyz_xxyzz, tr_z_yyyz_xyyyz, tr_z_yyyz_xyyzz, tr_z_yyyz_xyzzz, tr_z_yyyz_yyyyy, tr_z_yyyz_yyyyz, tr_z_yyyz_yyyzz, tr_z_yyyz_yyzzz, tr_z_yyyz_yzzzz, tr_z_yyyz_zzzzz, ts_xxyyy_xxxxy, ts_xxyyy_xxxyy, ts_xxyyy_xxyyy, ts_xxyyy_xyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyz_xxxxx[i] = 2.0 * tr_z_xxyz_xxxxx[i] * fe_0 + tr_z_xxyyz_xxxxx[i] * pa_y[i];

        tr_z_xxyyyz_xxxxy[i] = ts_xxyyy_xxxxy[i] * fe_0 + tr_z_xxyyy_xxxxy[i] * pa_z[i];

        tr_z_xxyyyz_xxxxz[i] = 2.0 * tr_z_xxyz_xxxxz[i] * fe_0 + tr_z_xxyyz_xxxxz[i] * pa_y[i];

        tr_z_xxyyyz_xxxyy[i] = ts_xxyyy_xxxyy[i] * fe_0 + tr_z_xxyyy_xxxyy[i] * pa_z[i];

        tr_z_xxyyyz_xxxyz[i] = tr_z_yyyz_xxxyz[i] * fe_0 + 3.0 * tr_z_xyyyz_xxyz[i] * fe_0 + tr_z_xyyyz_xxxyz[i] * pa_x[i];

        tr_z_xxyyyz_xxxzz[i] = 2.0 * tr_z_xxyz_xxxzz[i] * fe_0 + tr_z_xxyyz_xxxzz[i] * pa_y[i];

        tr_z_xxyyyz_xxyyy[i] = ts_xxyyy_xxyyy[i] * fe_0 + tr_z_xxyyy_xxyyy[i] * pa_z[i];

        tr_z_xxyyyz_xxyyz[i] = tr_z_yyyz_xxyyz[i] * fe_0 + 2.0 * tr_z_xyyyz_xyyz[i] * fe_0 + tr_z_xyyyz_xxyyz[i] * pa_x[i];

        tr_z_xxyyyz_xxyzz[i] = tr_z_yyyz_xxyzz[i] * fe_0 + 2.0 * tr_z_xyyyz_xyzz[i] * fe_0 + tr_z_xyyyz_xxyzz[i] * pa_x[i];

        tr_z_xxyyyz_xxzzz[i] = 2.0 * tr_z_xxyz_xxzzz[i] * fe_0 + tr_z_xxyyz_xxzzz[i] * pa_y[i];

        tr_z_xxyyyz_xyyyy[i] = ts_xxyyy_xyyyy[i] * fe_0 + tr_z_xxyyy_xyyyy[i] * pa_z[i];

        tr_z_xxyyyz_xyyyz[i] = tr_z_yyyz_xyyyz[i] * fe_0 + tr_z_xyyyz_yyyz[i] * fe_0 + tr_z_xyyyz_xyyyz[i] * pa_x[i];

        tr_z_xxyyyz_xyyzz[i] = tr_z_yyyz_xyyzz[i] * fe_0 + tr_z_xyyyz_yyzz[i] * fe_0 + tr_z_xyyyz_xyyzz[i] * pa_x[i];

        tr_z_xxyyyz_xyzzz[i] = tr_z_yyyz_xyzzz[i] * fe_0 + tr_z_xyyyz_yzzz[i] * fe_0 + tr_z_xyyyz_xyzzz[i] * pa_x[i];

        tr_z_xxyyyz_xzzzz[i] = 2.0 * tr_z_xxyz_xzzzz[i] * fe_0 + tr_z_xxyyz_xzzzz[i] * pa_y[i];

        tr_z_xxyyyz_yyyyy[i] = tr_z_yyyz_yyyyy[i] * fe_0 + tr_z_xyyyz_yyyyy[i] * pa_x[i];

        tr_z_xxyyyz_yyyyz[i] = tr_z_yyyz_yyyyz[i] * fe_0 + tr_z_xyyyz_yyyyz[i] * pa_x[i];

        tr_z_xxyyyz_yyyzz[i] = tr_z_yyyz_yyyzz[i] * fe_0 + tr_z_xyyyz_yyyzz[i] * pa_x[i];

        tr_z_xxyyyz_yyzzz[i] = tr_z_yyyz_yyzzz[i] * fe_0 + tr_z_xyyyz_yyzzz[i] * pa_x[i];

        tr_z_xxyyyz_yzzzz[i] = tr_z_yyyz_yzzzz[i] * fe_0 + tr_z_xyyyz_yzzzz[i] * pa_x[i];

        tr_z_xxyyyz_zzzzz[i] = tr_z_yyyz_zzzzz[i] * fe_0 + tr_z_xyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 1428-1449 components of targeted buffer : IH

    auto tr_z_xxyyzz_xxxxx = pbuffer.data(idx_dip_ih + 1428);

    auto tr_z_xxyyzz_xxxxy = pbuffer.data(idx_dip_ih + 1429);

    auto tr_z_xxyyzz_xxxxz = pbuffer.data(idx_dip_ih + 1430);

    auto tr_z_xxyyzz_xxxyy = pbuffer.data(idx_dip_ih + 1431);

    auto tr_z_xxyyzz_xxxyz = pbuffer.data(idx_dip_ih + 1432);

    auto tr_z_xxyyzz_xxxzz = pbuffer.data(idx_dip_ih + 1433);

    auto tr_z_xxyyzz_xxyyy = pbuffer.data(idx_dip_ih + 1434);

    auto tr_z_xxyyzz_xxyyz = pbuffer.data(idx_dip_ih + 1435);

    auto tr_z_xxyyzz_xxyzz = pbuffer.data(idx_dip_ih + 1436);

    auto tr_z_xxyyzz_xxzzz = pbuffer.data(idx_dip_ih + 1437);

    auto tr_z_xxyyzz_xyyyy = pbuffer.data(idx_dip_ih + 1438);

    auto tr_z_xxyyzz_xyyyz = pbuffer.data(idx_dip_ih + 1439);

    auto tr_z_xxyyzz_xyyzz = pbuffer.data(idx_dip_ih + 1440);

    auto tr_z_xxyyzz_xyzzz = pbuffer.data(idx_dip_ih + 1441);

    auto tr_z_xxyyzz_xzzzz = pbuffer.data(idx_dip_ih + 1442);

    auto tr_z_xxyyzz_yyyyy = pbuffer.data(idx_dip_ih + 1443);

    auto tr_z_xxyyzz_yyyyz = pbuffer.data(idx_dip_ih + 1444);

    auto tr_z_xxyyzz_yyyzz = pbuffer.data(idx_dip_ih + 1445);

    auto tr_z_xxyyzz_yyzzz = pbuffer.data(idx_dip_ih + 1446);

    auto tr_z_xxyyzz_yzzzz = pbuffer.data(idx_dip_ih + 1447);

    auto tr_z_xxyyzz_zzzzz = pbuffer.data(idx_dip_ih + 1448);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyyzz_xxxxx, tr_z_xxyyzz_xxxxy, tr_z_xxyyzz_xxxxz, tr_z_xxyyzz_xxxyy, tr_z_xxyyzz_xxxyz, tr_z_xxyyzz_xxxzz, tr_z_xxyyzz_xxyyy, tr_z_xxyyzz_xxyyz, tr_z_xxyyzz_xxyzz, tr_z_xxyyzz_xxzzz, tr_z_xxyyzz_xyyyy, tr_z_xxyyzz_xyyyz, tr_z_xxyyzz_xyyzz, tr_z_xxyyzz_xyzzz, tr_z_xxyyzz_xzzzz, tr_z_xxyyzz_yyyyy, tr_z_xxyyzz_yyyyz, tr_z_xxyyzz_yyyzz, tr_z_xxyyzz_yyzzz, tr_z_xxyyzz_yzzzz, tr_z_xxyyzz_zzzzz, tr_z_xxyzz_xxxxx, tr_z_xxyzz_xxxxz, tr_z_xxyzz_xxxzz, tr_z_xxyzz_xxzzz, tr_z_xxyzz_xzzzz, tr_z_xxzz_xxxxx, tr_z_xxzz_xxxxz, tr_z_xxzz_xxxzz, tr_z_xxzz_xxzzz, tr_z_xxzz_xzzzz, tr_z_xyyzz_xxxxy, tr_z_xyyzz_xxxy, tr_z_xyyzz_xxxyy, tr_z_xyyzz_xxxyz, tr_z_xyyzz_xxyy, tr_z_xyyzz_xxyyy, tr_z_xyyzz_xxyyz, tr_z_xyyzz_xxyz, tr_z_xyyzz_xxyzz, tr_z_xyyzz_xyyy, tr_z_xyyzz_xyyyy, tr_z_xyyzz_xyyyz, tr_z_xyyzz_xyyz, tr_z_xyyzz_xyyzz, tr_z_xyyzz_xyzz, tr_z_xyyzz_xyzzz, tr_z_xyyzz_yyyy, tr_z_xyyzz_yyyyy, tr_z_xyyzz_yyyyz, tr_z_xyyzz_yyyz, tr_z_xyyzz_yyyzz, tr_z_xyyzz_yyzz, tr_z_xyyzz_yyzzz, tr_z_xyyzz_yzzz, tr_z_xyyzz_yzzzz, tr_z_xyyzz_zzzzz, tr_z_yyzz_xxxxy, tr_z_yyzz_xxxyy, tr_z_yyzz_xxxyz, tr_z_yyzz_xxyyy, tr_z_yyzz_xxyyz, tr_z_yyzz_xxyzz, tr_z_yyzz_xyyyy, tr_z_yyzz_xyyyz, tr_z_yyzz_xyyzz, tr_z_yyzz_xyzzz, tr_z_yyzz_yyyyy, tr_z_yyzz_yyyyz, tr_z_yyzz_yyyzz, tr_z_yyzz_yyzzz, tr_z_yyzz_yzzzz, tr_z_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyzz_xxxxx[i] = tr_z_xxzz_xxxxx[i] * fe_0 + tr_z_xxyzz_xxxxx[i] * pa_y[i];

        tr_z_xxyyzz_xxxxy[i] = tr_z_yyzz_xxxxy[i] * fe_0 + 4.0 * tr_z_xyyzz_xxxy[i] * fe_0 + tr_z_xyyzz_xxxxy[i] * pa_x[i];

        tr_z_xxyyzz_xxxxz[i] = tr_z_xxzz_xxxxz[i] * fe_0 + tr_z_xxyzz_xxxxz[i] * pa_y[i];

        tr_z_xxyyzz_xxxyy[i] = tr_z_yyzz_xxxyy[i] * fe_0 + 3.0 * tr_z_xyyzz_xxyy[i] * fe_0 + tr_z_xyyzz_xxxyy[i] * pa_x[i];

        tr_z_xxyyzz_xxxyz[i] = tr_z_yyzz_xxxyz[i] * fe_0 + 3.0 * tr_z_xyyzz_xxyz[i] * fe_0 + tr_z_xyyzz_xxxyz[i] * pa_x[i];

        tr_z_xxyyzz_xxxzz[i] = tr_z_xxzz_xxxzz[i] * fe_0 + tr_z_xxyzz_xxxzz[i] * pa_y[i];

        tr_z_xxyyzz_xxyyy[i] = tr_z_yyzz_xxyyy[i] * fe_0 + 2.0 * tr_z_xyyzz_xyyy[i] * fe_0 + tr_z_xyyzz_xxyyy[i] * pa_x[i];

        tr_z_xxyyzz_xxyyz[i] = tr_z_yyzz_xxyyz[i] * fe_0 + 2.0 * tr_z_xyyzz_xyyz[i] * fe_0 + tr_z_xyyzz_xxyyz[i] * pa_x[i];

        tr_z_xxyyzz_xxyzz[i] = tr_z_yyzz_xxyzz[i] * fe_0 + 2.0 * tr_z_xyyzz_xyzz[i] * fe_0 + tr_z_xyyzz_xxyzz[i] * pa_x[i];

        tr_z_xxyyzz_xxzzz[i] = tr_z_xxzz_xxzzz[i] * fe_0 + tr_z_xxyzz_xxzzz[i] * pa_y[i];

        tr_z_xxyyzz_xyyyy[i] = tr_z_yyzz_xyyyy[i] * fe_0 + tr_z_xyyzz_yyyy[i] * fe_0 + tr_z_xyyzz_xyyyy[i] * pa_x[i];

        tr_z_xxyyzz_xyyyz[i] = tr_z_yyzz_xyyyz[i] * fe_0 + tr_z_xyyzz_yyyz[i] * fe_0 + tr_z_xyyzz_xyyyz[i] * pa_x[i];

        tr_z_xxyyzz_xyyzz[i] = tr_z_yyzz_xyyzz[i] * fe_0 + tr_z_xyyzz_yyzz[i] * fe_0 + tr_z_xyyzz_xyyzz[i] * pa_x[i];

        tr_z_xxyyzz_xyzzz[i] = tr_z_yyzz_xyzzz[i] * fe_0 + tr_z_xyyzz_yzzz[i] * fe_0 + tr_z_xyyzz_xyzzz[i] * pa_x[i];

        tr_z_xxyyzz_xzzzz[i] = tr_z_xxzz_xzzzz[i] * fe_0 + tr_z_xxyzz_xzzzz[i] * pa_y[i];

        tr_z_xxyyzz_yyyyy[i] = tr_z_yyzz_yyyyy[i] * fe_0 + tr_z_xyyzz_yyyyy[i] * pa_x[i];

        tr_z_xxyyzz_yyyyz[i] = tr_z_yyzz_yyyyz[i] * fe_0 + tr_z_xyyzz_yyyyz[i] * pa_x[i];

        tr_z_xxyyzz_yyyzz[i] = tr_z_yyzz_yyyzz[i] * fe_0 + tr_z_xyyzz_yyyzz[i] * pa_x[i];

        tr_z_xxyyzz_yyzzz[i] = tr_z_yyzz_yyzzz[i] * fe_0 + tr_z_xyyzz_yyzzz[i] * pa_x[i];

        tr_z_xxyyzz_yzzzz[i] = tr_z_yyzz_yzzzz[i] * fe_0 + tr_z_xyyzz_yzzzz[i] * pa_x[i];

        tr_z_xxyyzz_zzzzz[i] = tr_z_yyzz_zzzzz[i] * fe_0 + tr_z_xyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1449-1470 components of targeted buffer : IH

    auto tr_z_xxyzzz_xxxxx = pbuffer.data(idx_dip_ih + 1449);

    auto tr_z_xxyzzz_xxxxy = pbuffer.data(idx_dip_ih + 1450);

    auto tr_z_xxyzzz_xxxxz = pbuffer.data(idx_dip_ih + 1451);

    auto tr_z_xxyzzz_xxxyy = pbuffer.data(idx_dip_ih + 1452);

    auto tr_z_xxyzzz_xxxyz = pbuffer.data(idx_dip_ih + 1453);

    auto tr_z_xxyzzz_xxxzz = pbuffer.data(idx_dip_ih + 1454);

    auto tr_z_xxyzzz_xxyyy = pbuffer.data(idx_dip_ih + 1455);

    auto tr_z_xxyzzz_xxyyz = pbuffer.data(idx_dip_ih + 1456);

    auto tr_z_xxyzzz_xxyzz = pbuffer.data(idx_dip_ih + 1457);

    auto tr_z_xxyzzz_xxzzz = pbuffer.data(idx_dip_ih + 1458);

    auto tr_z_xxyzzz_xyyyy = pbuffer.data(idx_dip_ih + 1459);

    auto tr_z_xxyzzz_xyyyz = pbuffer.data(idx_dip_ih + 1460);

    auto tr_z_xxyzzz_xyyzz = pbuffer.data(idx_dip_ih + 1461);

    auto tr_z_xxyzzz_xyzzz = pbuffer.data(idx_dip_ih + 1462);

    auto tr_z_xxyzzz_xzzzz = pbuffer.data(idx_dip_ih + 1463);

    auto tr_z_xxyzzz_yyyyy = pbuffer.data(idx_dip_ih + 1464);

    auto tr_z_xxyzzz_yyyyz = pbuffer.data(idx_dip_ih + 1465);

    auto tr_z_xxyzzz_yyyzz = pbuffer.data(idx_dip_ih + 1466);

    auto tr_z_xxyzzz_yyzzz = pbuffer.data(idx_dip_ih + 1467);

    auto tr_z_xxyzzz_yzzzz = pbuffer.data(idx_dip_ih + 1468);

    auto tr_z_xxyzzz_zzzzz = pbuffer.data(idx_dip_ih + 1469);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyzzz_xxxxx, tr_z_xxyzzz_xxxxy, tr_z_xxyzzz_xxxxz, tr_z_xxyzzz_xxxyy, tr_z_xxyzzz_xxxyz, tr_z_xxyzzz_xxxzz, tr_z_xxyzzz_xxyyy, tr_z_xxyzzz_xxyyz, tr_z_xxyzzz_xxyzz, tr_z_xxyzzz_xxzzz, tr_z_xxyzzz_xyyyy, tr_z_xxyzzz_xyyyz, tr_z_xxyzzz_xyyzz, tr_z_xxyzzz_xyzzz, tr_z_xxyzzz_xzzzz, tr_z_xxyzzz_yyyyy, tr_z_xxyzzz_yyyyz, tr_z_xxyzzz_yyyzz, tr_z_xxyzzz_yyzzz, tr_z_xxyzzz_yzzzz, tr_z_xxyzzz_zzzzz, tr_z_xxzzz_xxxx, tr_z_xxzzz_xxxxx, tr_z_xxzzz_xxxxy, tr_z_xxzzz_xxxxz, tr_z_xxzzz_xxxy, tr_z_xxzzz_xxxyy, tr_z_xxzzz_xxxyz, tr_z_xxzzz_xxxz, tr_z_xxzzz_xxxzz, tr_z_xxzzz_xxyy, tr_z_xxzzz_xxyyy, tr_z_xxzzz_xxyyz, tr_z_xxzzz_xxyz, tr_z_xxzzz_xxyzz, tr_z_xxzzz_xxzz, tr_z_xxzzz_xxzzz, tr_z_xxzzz_xyyy, tr_z_xxzzz_xyyyy, tr_z_xxzzz_xyyyz, tr_z_xxzzz_xyyz, tr_z_xxzzz_xyyzz, tr_z_xxzzz_xyzz, tr_z_xxzzz_xyzzz, tr_z_xxzzz_xzzz, tr_z_xxzzz_xzzzz, tr_z_xxzzz_zzzzz, tr_z_xyzzz_yyyyy, tr_z_xyzzz_yyyyz, tr_z_xyzzz_yyyzz, tr_z_xyzzz_yyzzz, tr_z_xyzzz_yzzzz, tr_z_yzzz_yyyyy, tr_z_yzzz_yyyyz, tr_z_yzzz_yyyzz, tr_z_yzzz_yyzzz, tr_z_yzzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzzz_xxxxx[i] = tr_z_xxzzz_xxxxx[i] * pa_y[i];

        tr_z_xxyzzz_xxxxy[i] = tr_z_xxzzz_xxxx[i] * fe_0 + tr_z_xxzzz_xxxxy[i] * pa_y[i];

        tr_z_xxyzzz_xxxxz[i] = tr_z_xxzzz_xxxxz[i] * pa_y[i];

        tr_z_xxyzzz_xxxyy[i] = 2.0 * tr_z_xxzzz_xxxy[i] * fe_0 + tr_z_xxzzz_xxxyy[i] * pa_y[i];

        tr_z_xxyzzz_xxxyz[i] = tr_z_xxzzz_xxxz[i] * fe_0 + tr_z_xxzzz_xxxyz[i] * pa_y[i];

        tr_z_xxyzzz_xxxzz[i] = tr_z_xxzzz_xxxzz[i] * pa_y[i];

        tr_z_xxyzzz_xxyyy[i] = 3.0 * tr_z_xxzzz_xxyy[i] * fe_0 + tr_z_xxzzz_xxyyy[i] * pa_y[i];

        tr_z_xxyzzz_xxyyz[i] = 2.0 * tr_z_xxzzz_xxyz[i] * fe_0 + tr_z_xxzzz_xxyyz[i] * pa_y[i];

        tr_z_xxyzzz_xxyzz[i] = tr_z_xxzzz_xxzz[i] * fe_0 + tr_z_xxzzz_xxyzz[i] * pa_y[i];

        tr_z_xxyzzz_xxzzz[i] = tr_z_xxzzz_xxzzz[i] * pa_y[i];

        tr_z_xxyzzz_xyyyy[i] = 4.0 * tr_z_xxzzz_xyyy[i] * fe_0 + tr_z_xxzzz_xyyyy[i] * pa_y[i];

        tr_z_xxyzzz_xyyyz[i] = 3.0 * tr_z_xxzzz_xyyz[i] * fe_0 + tr_z_xxzzz_xyyyz[i] * pa_y[i];

        tr_z_xxyzzz_xyyzz[i] = 2.0 * tr_z_xxzzz_xyzz[i] * fe_0 + tr_z_xxzzz_xyyzz[i] * pa_y[i];

        tr_z_xxyzzz_xyzzz[i] = tr_z_xxzzz_xzzz[i] * fe_0 + tr_z_xxzzz_xyzzz[i] * pa_y[i];

        tr_z_xxyzzz_xzzzz[i] = tr_z_xxzzz_xzzzz[i] * pa_y[i];

        tr_z_xxyzzz_yyyyy[i] = tr_z_yzzz_yyyyy[i] * fe_0 + tr_z_xyzzz_yyyyy[i] * pa_x[i];

        tr_z_xxyzzz_yyyyz[i] = tr_z_yzzz_yyyyz[i] * fe_0 + tr_z_xyzzz_yyyyz[i] * pa_x[i];

        tr_z_xxyzzz_yyyzz[i] = tr_z_yzzz_yyyzz[i] * fe_0 + tr_z_xyzzz_yyyzz[i] * pa_x[i];

        tr_z_xxyzzz_yyzzz[i] = tr_z_yzzz_yyzzz[i] * fe_0 + tr_z_xyzzz_yyzzz[i] * pa_x[i];

        tr_z_xxyzzz_yzzzz[i] = tr_z_yzzz_yzzzz[i] * fe_0 + tr_z_xyzzz_yzzzz[i] * pa_x[i];

        tr_z_xxyzzz_zzzzz[i] = tr_z_xxzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1470-1491 components of targeted buffer : IH

    auto tr_z_xxzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1470);

    auto tr_z_xxzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1471);

    auto tr_z_xxzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1472);

    auto tr_z_xxzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1473);

    auto tr_z_xxzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1474);

    auto tr_z_xxzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1475);

    auto tr_z_xxzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1476);

    auto tr_z_xxzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1477);

    auto tr_z_xxzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1478);

    auto tr_z_xxzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1479);

    auto tr_z_xxzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1480);

    auto tr_z_xxzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1481);

    auto tr_z_xxzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1482);

    auto tr_z_xxzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1483);

    auto tr_z_xxzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1484);

    auto tr_z_xxzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1485);

    auto tr_z_xxzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1486);

    auto tr_z_xxzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1487);

    auto tr_z_xxzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1488);

    auto tr_z_xxzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1489);

    auto tr_z_xxzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1490);

    #pragma omp simd aligned(pa_x, tr_z_xxzzzz_xxxxx, tr_z_xxzzzz_xxxxy, tr_z_xxzzzz_xxxxz, tr_z_xxzzzz_xxxyy, tr_z_xxzzzz_xxxyz, tr_z_xxzzzz_xxxzz, tr_z_xxzzzz_xxyyy, tr_z_xxzzzz_xxyyz, tr_z_xxzzzz_xxyzz, tr_z_xxzzzz_xxzzz, tr_z_xxzzzz_xyyyy, tr_z_xxzzzz_xyyyz, tr_z_xxzzzz_xyyzz, tr_z_xxzzzz_xyzzz, tr_z_xxzzzz_xzzzz, tr_z_xxzzzz_yyyyy, tr_z_xxzzzz_yyyyz, tr_z_xxzzzz_yyyzz, tr_z_xxzzzz_yyzzz, tr_z_xxzzzz_yzzzz, tr_z_xxzzzz_zzzzz, tr_z_xzzzz_xxxx, tr_z_xzzzz_xxxxx, tr_z_xzzzz_xxxxy, tr_z_xzzzz_xxxxz, tr_z_xzzzz_xxxy, tr_z_xzzzz_xxxyy, tr_z_xzzzz_xxxyz, tr_z_xzzzz_xxxz, tr_z_xzzzz_xxxzz, tr_z_xzzzz_xxyy, tr_z_xzzzz_xxyyy, tr_z_xzzzz_xxyyz, tr_z_xzzzz_xxyz, tr_z_xzzzz_xxyzz, tr_z_xzzzz_xxzz, tr_z_xzzzz_xxzzz, tr_z_xzzzz_xyyy, tr_z_xzzzz_xyyyy, tr_z_xzzzz_xyyyz, tr_z_xzzzz_xyyz, tr_z_xzzzz_xyyzz, tr_z_xzzzz_xyzz, tr_z_xzzzz_xyzzz, tr_z_xzzzz_xzzz, tr_z_xzzzz_xzzzz, tr_z_xzzzz_yyyy, tr_z_xzzzz_yyyyy, tr_z_xzzzz_yyyyz, tr_z_xzzzz_yyyz, tr_z_xzzzz_yyyzz, tr_z_xzzzz_yyzz, tr_z_xzzzz_yyzzz, tr_z_xzzzz_yzzz, tr_z_xzzzz_yzzzz, tr_z_xzzzz_zzzz, tr_z_xzzzz_zzzzz, tr_z_zzzz_xxxxx, tr_z_zzzz_xxxxy, tr_z_zzzz_xxxxz, tr_z_zzzz_xxxyy, tr_z_zzzz_xxxyz, tr_z_zzzz_xxxzz, tr_z_zzzz_xxyyy, tr_z_zzzz_xxyyz, tr_z_zzzz_xxyzz, tr_z_zzzz_xxzzz, tr_z_zzzz_xyyyy, tr_z_zzzz_xyyyz, tr_z_zzzz_xyyzz, tr_z_zzzz_xyzzz, tr_z_zzzz_xzzzz, tr_z_zzzz_yyyyy, tr_z_zzzz_yyyyz, tr_z_zzzz_yyyzz, tr_z_zzzz_yyzzz, tr_z_zzzz_yzzzz, tr_z_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzzz_xxxxx[i] = tr_z_zzzz_xxxxx[i] * fe_0 + 5.0 * tr_z_xzzzz_xxxx[i] * fe_0 + tr_z_xzzzz_xxxxx[i] * pa_x[i];

        tr_z_xxzzzz_xxxxy[i] = tr_z_zzzz_xxxxy[i] * fe_0 + 4.0 * tr_z_xzzzz_xxxy[i] * fe_0 + tr_z_xzzzz_xxxxy[i] * pa_x[i];

        tr_z_xxzzzz_xxxxz[i] = tr_z_zzzz_xxxxz[i] * fe_0 + 4.0 * tr_z_xzzzz_xxxz[i] * fe_0 + tr_z_xzzzz_xxxxz[i] * pa_x[i];

        tr_z_xxzzzz_xxxyy[i] = tr_z_zzzz_xxxyy[i] * fe_0 + 3.0 * tr_z_xzzzz_xxyy[i] * fe_0 + tr_z_xzzzz_xxxyy[i] * pa_x[i];

        tr_z_xxzzzz_xxxyz[i] = tr_z_zzzz_xxxyz[i] * fe_0 + 3.0 * tr_z_xzzzz_xxyz[i] * fe_0 + tr_z_xzzzz_xxxyz[i] * pa_x[i];

        tr_z_xxzzzz_xxxzz[i] = tr_z_zzzz_xxxzz[i] * fe_0 + 3.0 * tr_z_xzzzz_xxzz[i] * fe_0 + tr_z_xzzzz_xxxzz[i] * pa_x[i];

        tr_z_xxzzzz_xxyyy[i] = tr_z_zzzz_xxyyy[i] * fe_0 + 2.0 * tr_z_xzzzz_xyyy[i] * fe_0 + tr_z_xzzzz_xxyyy[i] * pa_x[i];

        tr_z_xxzzzz_xxyyz[i] = tr_z_zzzz_xxyyz[i] * fe_0 + 2.0 * tr_z_xzzzz_xyyz[i] * fe_0 + tr_z_xzzzz_xxyyz[i] * pa_x[i];

        tr_z_xxzzzz_xxyzz[i] = tr_z_zzzz_xxyzz[i] * fe_0 + 2.0 * tr_z_xzzzz_xyzz[i] * fe_0 + tr_z_xzzzz_xxyzz[i] * pa_x[i];

        tr_z_xxzzzz_xxzzz[i] = tr_z_zzzz_xxzzz[i] * fe_0 + 2.0 * tr_z_xzzzz_xzzz[i] * fe_0 + tr_z_xzzzz_xxzzz[i] * pa_x[i];

        tr_z_xxzzzz_xyyyy[i] = tr_z_zzzz_xyyyy[i] * fe_0 + tr_z_xzzzz_yyyy[i] * fe_0 + tr_z_xzzzz_xyyyy[i] * pa_x[i];

        tr_z_xxzzzz_xyyyz[i] = tr_z_zzzz_xyyyz[i] * fe_0 + tr_z_xzzzz_yyyz[i] * fe_0 + tr_z_xzzzz_xyyyz[i] * pa_x[i];

        tr_z_xxzzzz_xyyzz[i] = tr_z_zzzz_xyyzz[i] * fe_0 + tr_z_xzzzz_yyzz[i] * fe_0 + tr_z_xzzzz_xyyzz[i] * pa_x[i];

        tr_z_xxzzzz_xyzzz[i] = tr_z_zzzz_xyzzz[i] * fe_0 + tr_z_xzzzz_yzzz[i] * fe_0 + tr_z_xzzzz_xyzzz[i] * pa_x[i];

        tr_z_xxzzzz_xzzzz[i] = tr_z_zzzz_xzzzz[i] * fe_0 + tr_z_xzzzz_zzzz[i] * fe_0 + tr_z_xzzzz_xzzzz[i] * pa_x[i];

        tr_z_xxzzzz_yyyyy[i] = tr_z_zzzz_yyyyy[i] * fe_0 + tr_z_xzzzz_yyyyy[i] * pa_x[i];

        tr_z_xxzzzz_yyyyz[i] = tr_z_zzzz_yyyyz[i] * fe_0 + tr_z_xzzzz_yyyyz[i] * pa_x[i];

        tr_z_xxzzzz_yyyzz[i] = tr_z_zzzz_yyyzz[i] * fe_0 + tr_z_xzzzz_yyyzz[i] * pa_x[i];

        tr_z_xxzzzz_yyzzz[i] = tr_z_zzzz_yyzzz[i] * fe_0 + tr_z_xzzzz_yyzzz[i] * pa_x[i];

        tr_z_xxzzzz_yzzzz[i] = tr_z_zzzz_yzzzz[i] * fe_0 + tr_z_xzzzz_yzzzz[i] * pa_x[i];

        tr_z_xxzzzz_zzzzz[i] = tr_z_zzzz_zzzzz[i] * fe_0 + tr_z_xzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1491-1512 components of targeted buffer : IH

    auto tr_z_xyyyyy_xxxxx = pbuffer.data(idx_dip_ih + 1491);

    auto tr_z_xyyyyy_xxxxy = pbuffer.data(idx_dip_ih + 1492);

    auto tr_z_xyyyyy_xxxxz = pbuffer.data(idx_dip_ih + 1493);

    auto tr_z_xyyyyy_xxxyy = pbuffer.data(idx_dip_ih + 1494);

    auto tr_z_xyyyyy_xxxyz = pbuffer.data(idx_dip_ih + 1495);

    auto tr_z_xyyyyy_xxxzz = pbuffer.data(idx_dip_ih + 1496);

    auto tr_z_xyyyyy_xxyyy = pbuffer.data(idx_dip_ih + 1497);

    auto tr_z_xyyyyy_xxyyz = pbuffer.data(idx_dip_ih + 1498);

    auto tr_z_xyyyyy_xxyzz = pbuffer.data(idx_dip_ih + 1499);

    auto tr_z_xyyyyy_xxzzz = pbuffer.data(idx_dip_ih + 1500);

    auto tr_z_xyyyyy_xyyyy = pbuffer.data(idx_dip_ih + 1501);

    auto tr_z_xyyyyy_xyyyz = pbuffer.data(idx_dip_ih + 1502);

    auto tr_z_xyyyyy_xyyzz = pbuffer.data(idx_dip_ih + 1503);

    auto tr_z_xyyyyy_xyzzz = pbuffer.data(idx_dip_ih + 1504);

    auto tr_z_xyyyyy_xzzzz = pbuffer.data(idx_dip_ih + 1505);

    auto tr_z_xyyyyy_yyyyy = pbuffer.data(idx_dip_ih + 1506);

    auto tr_z_xyyyyy_yyyyz = pbuffer.data(idx_dip_ih + 1507);

    auto tr_z_xyyyyy_yyyzz = pbuffer.data(idx_dip_ih + 1508);

    auto tr_z_xyyyyy_yyzzz = pbuffer.data(idx_dip_ih + 1509);

    auto tr_z_xyyyyy_yzzzz = pbuffer.data(idx_dip_ih + 1510);

    auto tr_z_xyyyyy_zzzzz = pbuffer.data(idx_dip_ih + 1511);

    #pragma omp simd aligned(pa_x, tr_z_xyyyyy_xxxxx, tr_z_xyyyyy_xxxxy, tr_z_xyyyyy_xxxxz, tr_z_xyyyyy_xxxyy, tr_z_xyyyyy_xxxyz, tr_z_xyyyyy_xxxzz, tr_z_xyyyyy_xxyyy, tr_z_xyyyyy_xxyyz, tr_z_xyyyyy_xxyzz, tr_z_xyyyyy_xxzzz, tr_z_xyyyyy_xyyyy, tr_z_xyyyyy_xyyyz, tr_z_xyyyyy_xyyzz, tr_z_xyyyyy_xyzzz, tr_z_xyyyyy_xzzzz, tr_z_xyyyyy_yyyyy, tr_z_xyyyyy_yyyyz, tr_z_xyyyyy_yyyzz, tr_z_xyyyyy_yyzzz, tr_z_xyyyyy_yzzzz, tr_z_xyyyyy_zzzzz, tr_z_yyyyy_xxxx, tr_z_yyyyy_xxxxx, tr_z_yyyyy_xxxxy, tr_z_yyyyy_xxxxz, tr_z_yyyyy_xxxy, tr_z_yyyyy_xxxyy, tr_z_yyyyy_xxxyz, tr_z_yyyyy_xxxz, tr_z_yyyyy_xxxzz, tr_z_yyyyy_xxyy, tr_z_yyyyy_xxyyy, tr_z_yyyyy_xxyyz, tr_z_yyyyy_xxyz, tr_z_yyyyy_xxyzz, tr_z_yyyyy_xxzz, tr_z_yyyyy_xxzzz, tr_z_yyyyy_xyyy, tr_z_yyyyy_xyyyy, tr_z_yyyyy_xyyyz, tr_z_yyyyy_xyyz, tr_z_yyyyy_xyyzz, tr_z_yyyyy_xyzz, tr_z_yyyyy_xyzzz, tr_z_yyyyy_xzzz, tr_z_yyyyy_xzzzz, tr_z_yyyyy_yyyy, tr_z_yyyyy_yyyyy, tr_z_yyyyy_yyyyz, tr_z_yyyyy_yyyz, tr_z_yyyyy_yyyzz, tr_z_yyyyy_yyzz, tr_z_yyyyy_yyzzz, tr_z_yyyyy_yzzz, tr_z_yyyyy_yzzzz, tr_z_yyyyy_zzzz, tr_z_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyy_xxxxx[i] = 5.0 * tr_z_yyyyy_xxxx[i] * fe_0 + tr_z_yyyyy_xxxxx[i] * pa_x[i];

        tr_z_xyyyyy_xxxxy[i] = 4.0 * tr_z_yyyyy_xxxy[i] * fe_0 + tr_z_yyyyy_xxxxy[i] * pa_x[i];

        tr_z_xyyyyy_xxxxz[i] = 4.0 * tr_z_yyyyy_xxxz[i] * fe_0 + tr_z_yyyyy_xxxxz[i] * pa_x[i];

        tr_z_xyyyyy_xxxyy[i] = 3.0 * tr_z_yyyyy_xxyy[i] * fe_0 + tr_z_yyyyy_xxxyy[i] * pa_x[i];

        tr_z_xyyyyy_xxxyz[i] = 3.0 * tr_z_yyyyy_xxyz[i] * fe_0 + tr_z_yyyyy_xxxyz[i] * pa_x[i];

        tr_z_xyyyyy_xxxzz[i] = 3.0 * tr_z_yyyyy_xxzz[i] * fe_0 + tr_z_yyyyy_xxxzz[i] * pa_x[i];

        tr_z_xyyyyy_xxyyy[i] = 2.0 * tr_z_yyyyy_xyyy[i] * fe_0 + tr_z_yyyyy_xxyyy[i] * pa_x[i];

        tr_z_xyyyyy_xxyyz[i] = 2.0 * tr_z_yyyyy_xyyz[i] * fe_0 + tr_z_yyyyy_xxyyz[i] * pa_x[i];

        tr_z_xyyyyy_xxyzz[i] = 2.0 * tr_z_yyyyy_xyzz[i] * fe_0 + tr_z_yyyyy_xxyzz[i] * pa_x[i];

        tr_z_xyyyyy_xxzzz[i] = 2.0 * tr_z_yyyyy_xzzz[i] * fe_0 + tr_z_yyyyy_xxzzz[i] * pa_x[i];

        tr_z_xyyyyy_xyyyy[i] = tr_z_yyyyy_yyyy[i] * fe_0 + tr_z_yyyyy_xyyyy[i] * pa_x[i];

        tr_z_xyyyyy_xyyyz[i] = tr_z_yyyyy_yyyz[i] * fe_0 + tr_z_yyyyy_xyyyz[i] * pa_x[i];

        tr_z_xyyyyy_xyyzz[i] = tr_z_yyyyy_yyzz[i] * fe_0 + tr_z_yyyyy_xyyzz[i] * pa_x[i];

        tr_z_xyyyyy_xyzzz[i] = tr_z_yyyyy_yzzz[i] * fe_0 + tr_z_yyyyy_xyzzz[i] * pa_x[i];

        tr_z_xyyyyy_xzzzz[i] = tr_z_yyyyy_zzzz[i] * fe_0 + tr_z_yyyyy_xzzzz[i] * pa_x[i];

        tr_z_xyyyyy_yyyyy[i] = tr_z_yyyyy_yyyyy[i] * pa_x[i];

        tr_z_xyyyyy_yyyyz[i] = tr_z_yyyyy_yyyyz[i] * pa_x[i];

        tr_z_xyyyyy_yyyzz[i] = tr_z_yyyyy_yyyzz[i] * pa_x[i];

        tr_z_xyyyyy_yyzzz[i] = tr_z_yyyyy_yyzzz[i] * pa_x[i];

        tr_z_xyyyyy_yzzzz[i] = tr_z_yyyyy_yzzzz[i] * pa_x[i];

        tr_z_xyyyyy_zzzzz[i] = tr_z_yyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 1512-1533 components of targeted buffer : IH

    auto tr_z_xyyyyz_xxxxx = pbuffer.data(idx_dip_ih + 1512);

    auto tr_z_xyyyyz_xxxxy = pbuffer.data(idx_dip_ih + 1513);

    auto tr_z_xyyyyz_xxxxz = pbuffer.data(idx_dip_ih + 1514);

    auto tr_z_xyyyyz_xxxyy = pbuffer.data(idx_dip_ih + 1515);

    auto tr_z_xyyyyz_xxxyz = pbuffer.data(idx_dip_ih + 1516);

    auto tr_z_xyyyyz_xxxzz = pbuffer.data(idx_dip_ih + 1517);

    auto tr_z_xyyyyz_xxyyy = pbuffer.data(idx_dip_ih + 1518);

    auto tr_z_xyyyyz_xxyyz = pbuffer.data(idx_dip_ih + 1519);

    auto tr_z_xyyyyz_xxyzz = pbuffer.data(idx_dip_ih + 1520);

    auto tr_z_xyyyyz_xxzzz = pbuffer.data(idx_dip_ih + 1521);

    auto tr_z_xyyyyz_xyyyy = pbuffer.data(idx_dip_ih + 1522);

    auto tr_z_xyyyyz_xyyyz = pbuffer.data(idx_dip_ih + 1523);

    auto tr_z_xyyyyz_xyyzz = pbuffer.data(idx_dip_ih + 1524);

    auto tr_z_xyyyyz_xyzzz = pbuffer.data(idx_dip_ih + 1525);

    auto tr_z_xyyyyz_xzzzz = pbuffer.data(idx_dip_ih + 1526);

    auto tr_z_xyyyyz_yyyyy = pbuffer.data(idx_dip_ih + 1527);

    auto tr_z_xyyyyz_yyyyz = pbuffer.data(idx_dip_ih + 1528);

    auto tr_z_xyyyyz_yyyzz = pbuffer.data(idx_dip_ih + 1529);

    auto tr_z_xyyyyz_yyzzz = pbuffer.data(idx_dip_ih + 1530);

    auto tr_z_xyyyyz_yzzzz = pbuffer.data(idx_dip_ih + 1531);

    auto tr_z_xyyyyz_zzzzz = pbuffer.data(idx_dip_ih + 1532);

    #pragma omp simd aligned(pa_x, tr_z_xyyyyz_xxxxx, tr_z_xyyyyz_xxxxy, tr_z_xyyyyz_xxxxz, tr_z_xyyyyz_xxxyy, tr_z_xyyyyz_xxxyz, tr_z_xyyyyz_xxxzz, tr_z_xyyyyz_xxyyy, tr_z_xyyyyz_xxyyz, tr_z_xyyyyz_xxyzz, tr_z_xyyyyz_xxzzz, tr_z_xyyyyz_xyyyy, tr_z_xyyyyz_xyyyz, tr_z_xyyyyz_xyyzz, tr_z_xyyyyz_xyzzz, tr_z_xyyyyz_xzzzz, tr_z_xyyyyz_yyyyy, tr_z_xyyyyz_yyyyz, tr_z_xyyyyz_yyyzz, tr_z_xyyyyz_yyzzz, tr_z_xyyyyz_yzzzz, tr_z_xyyyyz_zzzzz, tr_z_yyyyz_xxxx, tr_z_yyyyz_xxxxx, tr_z_yyyyz_xxxxy, tr_z_yyyyz_xxxxz, tr_z_yyyyz_xxxy, tr_z_yyyyz_xxxyy, tr_z_yyyyz_xxxyz, tr_z_yyyyz_xxxz, tr_z_yyyyz_xxxzz, tr_z_yyyyz_xxyy, tr_z_yyyyz_xxyyy, tr_z_yyyyz_xxyyz, tr_z_yyyyz_xxyz, tr_z_yyyyz_xxyzz, tr_z_yyyyz_xxzz, tr_z_yyyyz_xxzzz, tr_z_yyyyz_xyyy, tr_z_yyyyz_xyyyy, tr_z_yyyyz_xyyyz, tr_z_yyyyz_xyyz, tr_z_yyyyz_xyyzz, tr_z_yyyyz_xyzz, tr_z_yyyyz_xyzzz, tr_z_yyyyz_xzzz, tr_z_yyyyz_xzzzz, tr_z_yyyyz_yyyy, tr_z_yyyyz_yyyyy, tr_z_yyyyz_yyyyz, tr_z_yyyyz_yyyz, tr_z_yyyyz_yyyzz, tr_z_yyyyz_yyzz, tr_z_yyyyz_yyzzz, tr_z_yyyyz_yzzz, tr_z_yyyyz_yzzzz, tr_z_yyyyz_zzzz, tr_z_yyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyz_xxxxx[i] = 5.0 * tr_z_yyyyz_xxxx[i] * fe_0 + tr_z_yyyyz_xxxxx[i] * pa_x[i];

        tr_z_xyyyyz_xxxxy[i] = 4.0 * tr_z_yyyyz_xxxy[i] * fe_0 + tr_z_yyyyz_xxxxy[i] * pa_x[i];

        tr_z_xyyyyz_xxxxz[i] = 4.0 * tr_z_yyyyz_xxxz[i] * fe_0 + tr_z_yyyyz_xxxxz[i] * pa_x[i];

        tr_z_xyyyyz_xxxyy[i] = 3.0 * tr_z_yyyyz_xxyy[i] * fe_0 + tr_z_yyyyz_xxxyy[i] * pa_x[i];

        tr_z_xyyyyz_xxxyz[i] = 3.0 * tr_z_yyyyz_xxyz[i] * fe_0 + tr_z_yyyyz_xxxyz[i] * pa_x[i];

        tr_z_xyyyyz_xxxzz[i] = 3.0 * tr_z_yyyyz_xxzz[i] * fe_0 + tr_z_yyyyz_xxxzz[i] * pa_x[i];

        tr_z_xyyyyz_xxyyy[i] = 2.0 * tr_z_yyyyz_xyyy[i] * fe_0 + tr_z_yyyyz_xxyyy[i] * pa_x[i];

        tr_z_xyyyyz_xxyyz[i] = 2.0 * tr_z_yyyyz_xyyz[i] * fe_0 + tr_z_yyyyz_xxyyz[i] * pa_x[i];

        tr_z_xyyyyz_xxyzz[i] = 2.0 * tr_z_yyyyz_xyzz[i] * fe_0 + tr_z_yyyyz_xxyzz[i] * pa_x[i];

        tr_z_xyyyyz_xxzzz[i] = 2.0 * tr_z_yyyyz_xzzz[i] * fe_0 + tr_z_yyyyz_xxzzz[i] * pa_x[i];

        tr_z_xyyyyz_xyyyy[i] = tr_z_yyyyz_yyyy[i] * fe_0 + tr_z_yyyyz_xyyyy[i] * pa_x[i];

        tr_z_xyyyyz_xyyyz[i] = tr_z_yyyyz_yyyz[i] * fe_0 + tr_z_yyyyz_xyyyz[i] * pa_x[i];

        tr_z_xyyyyz_xyyzz[i] = tr_z_yyyyz_yyzz[i] * fe_0 + tr_z_yyyyz_xyyzz[i] * pa_x[i];

        tr_z_xyyyyz_xyzzz[i] = tr_z_yyyyz_yzzz[i] * fe_0 + tr_z_yyyyz_xyzzz[i] * pa_x[i];

        tr_z_xyyyyz_xzzzz[i] = tr_z_yyyyz_zzzz[i] * fe_0 + tr_z_yyyyz_xzzzz[i] * pa_x[i];

        tr_z_xyyyyz_yyyyy[i] = tr_z_yyyyz_yyyyy[i] * pa_x[i];

        tr_z_xyyyyz_yyyyz[i] = tr_z_yyyyz_yyyyz[i] * pa_x[i];

        tr_z_xyyyyz_yyyzz[i] = tr_z_yyyyz_yyyzz[i] * pa_x[i];

        tr_z_xyyyyz_yyzzz[i] = tr_z_yyyyz_yyzzz[i] * pa_x[i];

        tr_z_xyyyyz_yzzzz[i] = tr_z_yyyyz_yzzzz[i] * pa_x[i];

        tr_z_xyyyyz_zzzzz[i] = tr_z_yyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 1533-1554 components of targeted buffer : IH

    auto tr_z_xyyyzz_xxxxx = pbuffer.data(idx_dip_ih + 1533);

    auto tr_z_xyyyzz_xxxxy = pbuffer.data(idx_dip_ih + 1534);

    auto tr_z_xyyyzz_xxxxz = pbuffer.data(idx_dip_ih + 1535);

    auto tr_z_xyyyzz_xxxyy = pbuffer.data(idx_dip_ih + 1536);

    auto tr_z_xyyyzz_xxxyz = pbuffer.data(idx_dip_ih + 1537);

    auto tr_z_xyyyzz_xxxzz = pbuffer.data(idx_dip_ih + 1538);

    auto tr_z_xyyyzz_xxyyy = pbuffer.data(idx_dip_ih + 1539);

    auto tr_z_xyyyzz_xxyyz = pbuffer.data(idx_dip_ih + 1540);

    auto tr_z_xyyyzz_xxyzz = pbuffer.data(idx_dip_ih + 1541);

    auto tr_z_xyyyzz_xxzzz = pbuffer.data(idx_dip_ih + 1542);

    auto tr_z_xyyyzz_xyyyy = pbuffer.data(idx_dip_ih + 1543);

    auto tr_z_xyyyzz_xyyyz = pbuffer.data(idx_dip_ih + 1544);

    auto tr_z_xyyyzz_xyyzz = pbuffer.data(idx_dip_ih + 1545);

    auto tr_z_xyyyzz_xyzzz = pbuffer.data(idx_dip_ih + 1546);

    auto tr_z_xyyyzz_xzzzz = pbuffer.data(idx_dip_ih + 1547);

    auto tr_z_xyyyzz_yyyyy = pbuffer.data(idx_dip_ih + 1548);

    auto tr_z_xyyyzz_yyyyz = pbuffer.data(idx_dip_ih + 1549);

    auto tr_z_xyyyzz_yyyzz = pbuffer.data(idx_dip_ih + 1550);

    auto tr_z_xyyyzz_yyzzz = pbuffer.data(idx_dip_ih + 1551);

    auto tr_z_xyyyzz_yzzzz = pbuffer.data(idx_dip_ih + 1552);

    auto tr_z_xyyyzz_zzzzz = pbuffer.data(idx_dip_ih + 1553);

    #pragma omp simd aligned(pa_x, tr_z_xyyyzz_xxxxx, tr_z_xyyyzz_xxxxy, tr_z_xyyyzz_xxxxz, tr_z_xyyyzz_xxxyy, tr_z_xyyyzz_xxxyz, tr_z_xyyyzz_xxxzz, tr_z_xyyyzz_xxyyy, tr_z_xyyyzz_xxyyz, tr_z_xyyyzz_xxyzz, tr_z_xyyyzz_xxzzz, tr_z_xyyyzz_xyyyy, tr_z_xyyyzz_xyyyz, tr_z_xyyyzz_xyyzz, tr_z_xyyyzz_xyzzz, tr_z_xyyyzz_xzzzz, tr_z_xyyyzz_yyyyy, tr_z_xyyyzz_yyyyz, tr_z_xyyyzz_yyyzz, tr_z_xyyyzz_yyzzz, tr_z_xyyyzz_yzzzz, tr_z_xyyyzz_zzzzz, tr_z_yyyzz_xxxx, tr_z_yyyzz_xxxxx, tr_z_yyyzz_xxxxy, tr_z_yyyzz_xxxxz, tr_z_yyyzz_xxxy, tr_z_yyyzz_xxxyy, tr_z_yyyzz_xxxyz, tr_z_yyyzz_xxxz, tr_z_yyyzz_xxxzz, tr_z_yyyzz_xxyy, tr_z_yyyzz_xxyyy, tr_z_yyyzz_xxyyz, tr_z_yyyzz_xxyz, tr_z_yyyzz_xxyzz, tr_z_yyyzz_xxzz, tr_z_yyyzz_xxzzz, tr_z_yyyzz_xyyy, tr_z_yyyzz_xyyyy, tr_z_yyyzz_xyyyz, tr_z_yyyzz_xyyz, tr_z_yyyzz_xyyzz, tr_z_yyyzz_xyzz, tr_z_yyyzz_xyzzz, tr_z_yyyzz_xzzz, tr_z_yyyzz_xzzzz, tr_z_yyyzz_yyyy, tr_z_yyyzz_yyyyy, tr_z_yyyzz_yyyyz, tr_z_yyyzz_yyyz, tr_z_yyyzz_yyyzz, tr_z_yyyzz_yyzz, tr_z_yyyzz_yyzzz, tr_z_yyyzz_yzzz, tr_z_yyyzz_yzzzz, tr_z_yyyzz_zzzz, tr_z_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyzz_xxxxx[i] = 5.0 * tr_z_yyyzz_xxxx[i] * fe_0 + tr_z_yyyzz_xxxxx[i] * pa_x[i];

        tr_z_xyyyzz_xxxxy[i] = 4.0 * tr_z_yyyzz_xxxy[i] * fe_0 + tr_z_yyyzz_xxxxy[i] * pa_x[i];

        tr_z_xyyyzz_xxxxz[i] = 4.0 * tr_z_yyyzz_xxxz[i] * fe_0 + tr_z_yyyzz_xxxxz[i] * pa_x[i];

        tr_z_xyyyzz_xxxyy[i] = 3.0 * tr_z_yyyzz_xxyy[i] * fe_0 + tr_z_yyyzz_xxxyy[i] * pa_x[i];

        tr_z_xyyyzz_xxxyz[i] = 3.0 * tr_z_yyyzz_xxyz[i] * fe_0 + tr_z_yyyzz_xxxyz[i] * pa_x[i];

        tr_z_xyyyzz_xxxzz[i] = 3.0 * tr_z_yyyzz_xxzz[i] * fe_0 + tr_z_yyyzz_xxxzz[i] * pa_x[i];

        tr_z_xyyyzz_xxyyy[i] = 2.0 * tr_z_yyyzz_xyyy[i] * fe_0 + tr_z_yyyzz_xxyyy[i] * pa_x[i];

        tr_z_xyyyzz_xxyyz[i] = 2.0 * tr_z_yyyzz_xyyz[i] * fe_0 + tr_z_yyyzz_xxyyz[i] * pa_x[i];

        tr_z_xyyyzz_xxyzz[i] = 2.0 * tr_z_yyyzz_xyzz[i] * fe_0 + tr_z_yyyzz_xxyzz[i] * pa_x[i];

        tr_z_xyyyzz_xxzzz[i] = 2.0 * tr_z_yyyzz_xzzz[i] * fe_0 + tr_z_yyyzz_xxzzz[i] * pa_x[i];

        tr_z_xyyyzz_xyyyy[i] = tr_z_yyyzz_yyyy[i] * fe_0 + tr_z_yyyzz_xyyyy[i] * pa_x[i];

        tr_z_xyyyzz_xyyyz[i] = tr_z_yyyzz_yyyz[i] * fe_0 + tr_z_yyyzz_xyyyz[i] * pa_x[i];

        tr_z_xyyyzz_xyyzz[i] = tr_z_yyyzz_yyzz[i] * fe_0 + tr_z_yyyzz_xyyzz[i] * pa_x[i];

        tr_z_xyyyzz_xyzzz[i] = tr_z_yyyzz_yzzz[i] * fe_0 + tr_z_yyyzz_xyzzz[i] * pa_x[i];

        tr_z_xyyyzz_xzzzz[i] = tr_z_yyyzz_zzzz[i] * fe_0 + tr_z_yyyzz_xzzzz[i] * pa_x[i];

        tr_z_xyyyzz_yyyyy[i] = tr_z_yyyzz_yyyyy[i] * pa_x[i];

        tr_z_xyyyzz_yyyyz[i] = tr_z_yyyzz_yyyyz[i] * pa_x[i];

        tr_z_xyyyzz_yyyzz[i] = tr_z_yyyzz_yyyzz[i] * pa_x[i];

        tr_z_xyyyzz_yyzzz[i] = tr_z_yyyzz_yyzzz[i] * pa_x[i];

        tr_z_xyyyzz_yzzzz[i] = tr_z_yyyzz_yzzzz[i] * pa_x[i];

        tr_z_xyyyzz_zzzzz[i] = tr_z_yyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1554-1575 components of targeted buffer : IH

    auto tr_z_xyyzzz_xxxxx = pbuffer.data(idx_dip_ih + 1554);

    auto tr_z_xyyzzz_xxxxy = pbuffer.data(idx_dip_ih + 1555);

    auto tr_z_xyyzzz_xxxxz = pbuffer.data(idx_dip_ih + 1556);

    auto tr_z_xyyzzz_xxxyy = pbuffer.data(idx_dip_ih + 1557);

    auto tr_z_xyyzzz_xxxyz = pbuffer.data(idx_dip_ih + 1558);

    auto tr_z_xyyzzz_xxxzz = pbuffer.data(idx_dip_ih + 1559);

    auto tr_z_xyyzzz_xxyyy = pbuffer.data(idx_dip_ih + 1560);

    auto tr_z_xyyzzz_xxyyz = pbuffer.data(idx_dip_ih + 1561);

    auto tr_z_xyyzzz_xxyzz = pbuffer.data(idx_dip_ih + 1562);

    auto tr_z_xyyzzz_xxzzz = pbuffer.data(idx_dip_ih + 1563);

    auto tr_z_xyyzzz_xyyyy = pbuffer.data(idx_dip_ih + 1564);

    auto tr_z_xyyzzz_xyyyz = pbuffer.data(idx_dip_ih + 1565);

    auto tr_z_xyyzzz_xyyzz = pbuffer.data(idx_dip_ih + 1566);

    auto tr_z_xyyzzz_xyzzz = pbuffer.data(idx_dip_ih + 1567);

    auto tr_z_xyyzzz_xzzzz = pbuffer.data(idx_dip_ih + 1568);

    auto tr_z_xyyzzz_yyyyy = pbuffer.data(idx_dip_ih + 1569);

    auto tr_z_xyyzzz_yyyyz = pbuffer.data(idx_dip_ih + 1570);

    auto tr_z_xyyzzz_yyyzz = pbuffer.data(idx_dip_ih + 1571);

    auto tr_z_xyyzzz_yyzzz = pbuffer.data(idx_dip_ih + 1572);

    auto tr_z_xyyzzz_yzzzz = pbuffer.data(idx_dip_ih + 1573);

    auto tr_z_xyyzzz_zzzzz = pbuffer.data(idx_dip_ih + 1574);

    #pragma omp simd aligned(pa_x, tr_z_xyyzzz_xxxxx, tr_z_xyyzzz_xxxxy, tr_z_xyyzzz_xxxxz, tr_z_xyyzzz_xxxyy, tr_z_xyyzzz_xxxyz, tr_z_xyyzzz_xxxzz, tr_z_xyyzzz_xxyyy, tr_z_xyyzzz_xxyyz, tr_z_xyyzzz_xxyzz, tr_z_xyyzzz_xxzzz, tr_z_xyyzzz_xyyyy, tr_z_xyyzzz_xyyyz, tr_z_xyyzzz_xyyzz, tr_z_xyyzzz_xyzzz, tr_z_xyyzzz_xzzzz, tr_z_xyyzzz_yyyyy, tr_z_xyyzzz_yyyyz, tr_z_xyyzzz_yyyzz, tr_z_xyyzzz_yyzzz, tr_z_xyyzzz_yzzzz, tr_z_xyyzzz_zzzzz, tr_z_yyzzz_xxxx, tr_z_yyzzz_xxxxx, tr_z_yyzzz_xxxxy, tr_z_yyzzz_xxxxz, tr_z_yyzzz_xxxy, tr_z_yyzzz_xxxyy, tr_z_yyzzz_xxxyz, tr_z_yyzzz_xxxz, tr_z_yyzzz_xxxzz, tr_z_yyzzz_xxyy, tr_z_yyzzz_xxyyy, tr_z_yyzzz_xxyyz, tr_z_yyzzz_xxyz, tr_z_yyzzz_xxyzz, tr_z_yyzzz_xxzz, tr_z_yyzzz_xxzzz, tr_z_yyzzz_xyyy, tr_z_yyzzz_xyyyy, tr_z_yyzzz_xyyyz, tr_z_yyzzz_xyyz, tr_z_yyzzz_xyyzz, tr_z_yyzzz_xyzz, tr_z_yyzzz_xyzzz, tr_z_yyzzz_xzzz, tr_z_yyzzz_xzzzz, tr_z_yyzzz_yyyy, tr_z_yyzzz_yyyyy, tr_z_yyzzz_yyyyz, tr_z_yyzzz_yyyz, tr_z_yyzzz_yyyzz, tr_z_yyzzz_yyzz, tr_z_yyzzz_yyzzz, tr_z_yyzzz_yzzz, tr_z_yyzzz_yzzzz, tr_z_yyzzz_zzzz, tr_z_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzzz_xxxxx[i] = 5.0 * tr_z_yyzzz_xxxx[i] * fe_0 + tr_z_yyzzz_xxxxx[i] * pa_x[i];

        tr_z_xyyzzz_xxxxy[i] = 4.0 * tr_z_yyzzz_xxxy[i] * fe_0 + tr_z_yyzzz_xxxxy[i] * pa_x[i];

        tr_z_xyyzzz_xxxxz[i] = 4.0 * tr_z_yyzzz_xxxz[i] * fe_0 + tr_z_yyzzz_xxxxz[i] * pa_x[i];

        tr_z_xyyzzz_xxxyy[i] = 3.0 * tr_z_yyzzz_xxyy[i] * fe_0 + tr_z_yyzzz_xxxyy[i] * pa_x[i];

        tr_z_xyyzzz_xxxyz[i] = 3.0 * tr_z_yyzzz_xxyz[i] * fe_0 + tr_z_yyzzz_xxxyz[i] * pa_x[i];

        tr_z_xyyzzz_xxxzz[i] = 3.0 * tr_z_yyzzz_xxzz[i] * fe_0 + tr_z_yyzzz_xxxzz[i] * pa_x[i];

        tr_z_xyyzzz_xxyyy[i] = 2.0 * tr_z_yyzzz_xyyy[i] * fe_0 + tr_z_yyzzz_xxyyy[i] * pa_x[i];

        tr_z_xyyzzz_xxyyz[i] = 2.0 * tr_z_yyzzz_xyyz[i] * fe_0 + tr_z_yyzzz_xxyyz[i] * pa_x[i];

        tr_z_xyyzzz_xxyzz[i] = 2.0 * tr_z_yyzzz_xyzz[i] * fe_0 + tr_z_yyzzz_xxyzz[i] * pa_x[i];

        tr_z_xyyzzz_xxzzz[i] = 2.0 * tr_z_yyzzz_xzzz[i] * fe_0 + tr_z_yyzzz_xxzzz[i] * pa_x[i];

        tr_z_xyyzzz_xyyyy[i] = tr_z_yyzzz_yyyy[i] * fe_0 + tr_z_yyzzz_xyyyy[i] * pa_x[i];

        tr_z_xyyzzz_xyyyz[i] = tr_z_yyzzz_yyyz[i] * fe_0 + tr_z_yyzzz_xyyyz[i] * pa_x[i];

        tr_z_xyyzzz_xyyzz[i] = tr_z_yyzzz_yyzz[i] * fe_0 + tr_z_yyzzz_xyyzz[i] * pa_x[i];

        tr_z_xyyzzz_xyzzz[i] = tr_z_yyzzz_yzzz[i] * fe_0 + tr_z_yyzzz_xyzzz[i] * pa_x[i];

        tr_z_xyyzzz_xzzzz[i] = tr_z_yyzzz_zzzz[i] * fe_0 + tr_z_yyzzz_xzzzz[i] * pa_x[i];

        tr_z_xyyzzz_yyyyy[i] = tr_z_yyzzz_yyyyy[i] * pa_x[i];

        tr_z_xyyzzz_yyyyz[i] = tr_z_yyzzz_yyyyz[i] * pa_x[i];

        tr_z_xyyzzz_yyyzz[i] = tr_z_yyzzz_yyyzz[i] * pa_x[i];

        tr_z_xyyzzz_yyzzz[i] = tr_z_yyzzz_yyzzz[i] * pa_x[i];

        tr_z_xyyzzz_yzzzz[i] = tr_z_yyzzz_yzzzz[i] * pa_x[i];

        tr_z_xyyzzz_zzzzz[i] = tr_z_yyzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1575-1596 components of targeted buffer : IH

    auto tr_z_xyzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1575);

    auto tr_z_xyzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1576);

    auto tr_z_xyzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1577);

    auto tr_z_xyzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1578);

    auto tr_z_xyzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1579);

    auto tr_z_xyzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1580);

    auto tr_z_xyzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1581);

    auto tr_z_xyzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1582);

    auto tr_z_xyzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1583);

    auto tr_z_xyzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1584);

    auto tr_z_xyzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1585);

    auto tr_z_xyzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1586);

    auto tr_z_xyzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1587);

    auto tr_z_xyzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1588);

    auto tr_z_xyzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1589);

    auto tr_z_xyzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1590);

    auto tr_z_xyzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1591);

    auto tr_z_xyzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1592);

    auto tr_z_xyzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1593);

    auto tr_z_xyzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1594);

    auto tr_z_xyzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1595);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyzzzz_xxxxx, tr_z_xyzzzz_xxxxy, tr_z_xyzzzz_xxxxz, tr_z_xyzzzz_xxxyy, tr_z_xyzzzz_xxxyz, tr_z_xyzzzz_xxxzz, tr_z_xyzzzz_xxyyy, tr_z_xyzzzz_xxyyz, tr_z_xyzzzz_xxyzz, tr_z_xyzzzz_xxzzz, tr_z_xyzzzz_xyyyy, tr_z_xyzzzz_xyyyz, tr_z_xyzzzz_xyyzz, tr_z_xyzzzz_xyzzz, tr_z_xyzzzz_xzzzz, tr_z_xyzzzz_yyyyy, tr_z_xyzzzz_yyyyz, tr_z_xyzzzz_yyyzz, tr_z_xyzzzz_yyzzz, tr_z_xyzzzz_yzzzz, tr_z_xyzzzz_zzzzz, tr_z_xzzzz_xxxxx, tr_z_xzzzz_xxxxz, tr_z_xzzzz_xxxzz, tr_z_xzzzz_xxzzz, tr_z_xzzzz_xzzzz, tr_z_yzzzz_xxxxy, tr_z_yzzzz_xxxy, tr_z_yzzzz_xxxyy, tr_z_yzzzz_xxxyz, tr_z_yzzzz_xxyy, tr_z_yzzzz_xxyyy, tr_z_yzzzz_xxyyz, tr_z_yzzzz_xxyz, tr_z_yzzzz_xxyzz, tr_z_yzzzz_xyyy, tr_z_yzzzz_xyyyy, tr_z_yzzzz_xyyyz, tr_z_yzzzz_xyyz, tr_z_yzzzz_xyyzz, tr_z_yzzzz_xyzz, tr_z_yzzzz_xyzzz, tr_z_yzzzz_yyyy, tr_z_yzzzz_yyyyy, tr_z_yzzzz_yyyyz, tr_z_yzzzz_yyyz, tr_z_yzzzz_yyyzz, tr_z_yzzzz_yyzz, tr_z_yzzzz_yyzzz, tr_z_yzzzz_yzzz, tr_z_yzzzz_yzzzz, tr_z_yzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzzz_xxxxx[i] = tr_z_xzzzz_xxxxx[i] * pa_y[i];

        tr_z_xyzzzz_xxxxy[i] = 4.0 * tr_z_yzzzz_xxxy[i] * fe_0 + tr_z_yzzzz_xxxxy[i] * pa_x[i];

        tr_z_xyzzzz_xxxxz[i] = tr_z_xzzzz_xxxxz[i] * pa_y[i];

        tr_z_xyzzzz_xxxyy[i] = 3.0 * tr_z_yzzzz_xxyy[i] * fe_0 + tr_z_yzzzz_xxxyy[i] * pa_x[i];

        tr_z_xyzzzz_xxxyz[i] = 3.0 * tr_z_yzzzz_xxyz[i] * fe_0 + tr_z_yzzzz_xxxyz[i] * pa_x[i];

        tr_z_xyzzzz_xxxzz[i] = tr_z_xzzzz_xxxzz[i] * pa_y[i];

        tr_z_xyzzzz_xxyyy[i] = 2.0 * tr_z_yzzzz_xyyy[i] * fe_0 + tr_z_yzzzz_xxyyy[i] * pa_x[i];

        tr_z_xyzzzz_xxyyz[i] = 2.0 * tr_z_yzzzz_xyyz[i] * fe_0 + tr_z_yzzzz_xxyyz[i] * pa_x[i];

        tr_z_xyzzzz_xxyzz[i] = 2.0 * tr_z_yzzzz_xyzz[i] * fe_0 + tr_z_yzzzz_xxyzz[i] * pa_x[i];

        tr_z_xyzzzz_xxzzz[i] = tr_z_xzzzz_xxzzz[i] * pa_y[i];

        tr_z_xyzzzz_xyyyy[i] = tr_z_yzzzz_yyyy[i] * fe_0 + tr_z_yzzzz_xyyyy[i] * pa_x[i];

        tr_z_xyzzzz_xyyyz[i] = tr_z_yzzzz_yyyz[i] * fe_0 + tr_z_yzzzz_xyyyz[i] * pa_x[i];

        tr_z_xyzzzz_xyyzz[i] = tr_z_yzzzz_yyzz[i] * fe_0 + tr_z_yzzzz_xyyzz[i] * pa_x[i];

        tr_z_xyzzzz_xyzzz[i] = tr_z_yzzzz_yzzz[i] * fe_0 + tr_z_yzzzz_xyzzz[i] * pa_x[i];

        tr_z_xyzzzz_xzzzz[i] = tr_z_xzzzz_xzzzz[i] * pa_y[i];

        tr_z_xyzzzz_yyyyy[i] = tr_z_yzzzz_yyyyy[i] * pa_x[i];

        tr_z_xyzzzz_yyyyz[i] = tr_z_yzzzz_yyyyz[i] * pa_x[i];

        tr_z_xyzzzz_yyyzz[i] = tr_z_yzzzz_yyyzz[i] * pa_x[i];

        tr_z_xyzzzz_yyzzz[i] = tr_z_yzzzz_yyzzz[i] * pa_x[i];

        tr_z_xyzzzz_yzzzz[i] = tr_z_yzzzz_yzzzz[i] * pa_x[i];

        tr_z_xyzzzz_zzzzz[i] = tr_z_yzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1596-1617 components of targeted buffer : IH

    auto tr_z_xzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1596);

    auto tr_z_xzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1597);

    auto tr_z_xzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1598);

    auto tr_z_xzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1599);

    auto tr_z_xzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1600);

    auto tr_z_xzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1601);

    auto tr_z_xzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1602);

    auto tr_z_xzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1603);

    auto tr_z_xzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1604);

    auto tr_z_xzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1605);

    auto tr_z_xzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1606);

    auto tr_z_xzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1607);

    auto tr_z_xzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1608);

    auto tr_z_xzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1609);

    auto tr_z_xzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1610);

    auto tr_z_xzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1611);

    auto tr_z_xzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1612);

    auto tr_z_xzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1613);

    auto tr_z_xzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1614);

    auto tr_z_xzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1615);

    auto tr_z_xzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1616);

    #pragma omp simd aligned(pa_x, tr_z_xzzzzz_xxxxx, tr_z_xzzzzz_xxxxy, tr_z_xzzzzz_xxxxz, tr_z_xzzzzz_xxxyy, tr_z_xzzzzz_xxxyz, tr_z_xzzzzz_xxxzz, tr_z_xzzzzz_xxyyy, tr_z_xzzzzz_xxyyz, tr_z_xzzzzz_xxyzz, tr_z_xzzzzz_xxzzz, tr_z_xzzzzz_xyyyy, tr_z_xzzzzz_xyyyz, tr_z_xzzzzz_xyyzz, tr_z_xzzzzz_xyzzz, tr_z_xzzzzz_xzzzz, tr_z_xzzzzz_yyyyy, tr_z_xzzzzz_yyyyz, tr_z_xzzzzz_yyyzz, tr_z_xzzzzz_yyzzz, tr_z_xzzzzz_yzzzz, tr_z_xzzzzz_zzzzz, tr_z_zzzzz_xxxx, tr_z_zzzzz_xxxxx, tr_z_zzzzz_xxxxy, tr_z_zzzzz_xxxxz, tr_z_zzzzz_xxxy, tr_z_zzzzz_xxxyy, tr_z_zzzzz_xxxyz, tr_z_zzzzz_xxxz, tr_z_zzzzz_xxxzz, tr_z_zzzzz_xxyy, tr_z_zzzzz_xxyyy, tr_z_zzzzz_xxyyz, tr_z_zzzzz_xxyz, tr_z_zzzzz_xxyzz, tr_z_zzzzz_xxzz, tr_z_zzzzz_xxzzz, tr_z_zzzzz_xyyy, tr_z_zzzzz_xyyyy, tr_z_zzzzz_xyyyz, tr_z_zzzzz_xyyz, tr_z_zzzzz_xyyzz, tr_z_zzzzz_xyzz, tr_z_zzzzz_xyzzz, tr_z_zzzzz_xzzz, tr_z_zzzzz_xzzzz, tr_z_zzzzz_yyyy, tr_z_zzzzz_yyyyy, tr_z_zzzzz_yyyyz, tr_z_zzzzz_yyyz, tr_z_zzzzz_yyyzz, tr_z_zzzzz_yyzz, tr_z_zzzzz_yyzzz, tr_z_zzzzz_yzzz, tr_z_zzzzz_yzzzz, tr_z_zzzzz_zzzz, tr_z_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzzz_xxxxx[i] = 5.0 * tr_z_zzzzz_xxxx[i] * fe_0 + tr_z_zzzzz_xxxxx[i] * pa_x[i];

        tr_z_xzzzzz_xxxxy[i] = 4.0 * tr_z_zzzzz_xxxy[i] * fe_0 + tr_z_zzzzz_xxxxy[i] * pa_x[i];

        tr_z_xzzzzz_xxxxz[i] = 4.0 * tr_z_zzzzz_xxxz[i] * fe_0 + tr_z_zzzzz_xxxxz[i] * pa_x[i];

        tr_z_xzzzzz_xxxyy[i] = 3.0 * tr_z_zzzzz_xxyy[i] * fe_0 + tr_z_zzzzz_xxxyy[i] * pa_x[i];

        tr_z_xzzzzz_xxxyz[i] = 3.0 * tr_z_zzzzz_xxyz[i] * fe_0 + tr_z_zzzzz_xxxyz[i] * pa_x[i];

        tr_z_xzzzzz_xxxzz[i] = 3.0 * tr_z_zzzzz_xxzz[i] * fe_0 + tr_z_zzzzz_xxxzz[i] * pa_x[i];

        tr_z_xzzzzz_xxyyy[i] = 2.0 * tr_z_zzzzz_xyyy[i] * fe_0 + tr_z_zzzzz_xxyyy[i] * pa_x[i];

        tr_z_xzzzzz_xxyyz[i] = 2.0 * tr_z_zzzzz_xyyz[i] * fe_0 + tr_z_zzzzz_xxyyz[i] * pa_x[i];

        tr_z_xzzzzz_xxyzz[i] = 2.0 * tr_z_zzzzz_xyzz[i] * fe_0 + tr_z_zzzzz_xxyzz[i] * pa_x[i];

        tr_z_xzzzzz_xxzzz[i] = 2.0 * tr_z_zzzzz_xzzz[i] * fe_0 + tr_z_zzzzz_xxzzz[i] * pa_x[i];

        tr_z_xzzzzz_xyyyy[i] = tr_z_zzzzz_yyyy[i] * fe_0 + tr_z_zzzzz_xyyyy[i] * pa_x[i];

        tr_z_xzzzzz_xyyyz[i] = tr_z_zzzzz_yyyz[i] * fe_0 + tr_z_zzzzz_xyyyz[i] * pa_x[i];

        tr_z_xzzzzz_xyyzz[i] = tr_z_zzzzz_yyzz[i] * fe_0 + tr_z_zzzzz_xyyzz[i] * pa_x[i];

        tr_z_xzzzzz_xyzzz[i] = tr_z_zzzzz_yzzz[i] * fe_0 + tr_z_zzzzz_xyzzz[i] * pa_x[i];

        tr_z_xzzzzz_xzzzz[i] = tr_z_zzzzz_zzzz[i] * fe_0 + tr_z_zzzzz_xzzzz[i] * pa_x[i];

        tr_z_xzzzzz_yyyyy[i] = tr_z_zzzzz_yyyyy[i] * pa_x[i];

        tr_z_xzzzzz_yyyyz[i] = tr_z_zzzzz_yyyyz[i] * pa_x[i];

        tr_z_xzzzzz_yyyzz[i] = tr_z_zzzzz_yyyzz[i] * pa_x[i];

        tr_z_xzzzzz_yyzzz[i] = tr_z_zzzzz_yyzzz[i] * pa_x[i];

        tr_z_xzzzzz_yzzzz[i] = tr_z_zzzzz_yzzzz[i] * pa_x[i];

        tr_z_xzzzzz_zzzzz[i] = tr_z_zzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1617-1638 components of targeted buffer : IH

    auto tr_z_yyyyyy_xxxxx = pbuffer.data(idx_dip_ih + 1617);

    auto tr_z_yyyyyy_xxxxy = pbuffer.data(idx_dip_ih + 1618);

    auto tr_z_yyyyyy_xxxxz = pbuffer.data(idx_dip_ih + 1619);

    auto tr_z_yyyyyy_xxxyy = pbuffer.data(idx_dip_ih + 1620);

    auto tr_z_yyyyyy_xxxyz = pbuffer.data(idx_dip_ih + 1621);

    auto tr_z_yyyyyy_xxxzz = pbuffer.data(idx_dip_ih + 1622);

    auto tr_z_yyyyyy_xxyyy = pbuffer.data(idx_dip_ih + 1623);

    auto tr_z_yyyyyy_xxyyz = pbuffer.data(idx_dip_ih + 1624);

    auto tr_z_yyyyyy_xxyzz = pbuffer.data(idx_dip_ih + 1625);

    auto tr_z_yyyyyy_xxzzz = pbuffer.data(idx_dip_ih + 1626);

    auto tr_z_yyyyyy_xyyyy = pbuffer.data(idx_dip_ih + 1627);

    auto tr_z_yyyyyy_xyyyz = pbuffer.data(idx_dip_ih + 1628);

    auto tr_z_yyyyyy_xyyzz = pbuffer.data(idx_dip_ih + 1629);

    auto tr_z_yyyyyy_xyzzz = pbuffer.data(idx_dip_ih + 1630);

    auto tr_z_yyyyyy_xzzzz = pbuffer.data(idx_dip_ih + 1631);

    auto tr_z_yyyyyy_yyyyy = pbuffer.data(idx_dip_ih + 1632);

    auto tr_z_yyyyyy_yyyyz = pbuffer.data(idx_dip_ih + 1633);

    auto tr_z_yyyyyy_yyyzz = pbuffer.data(idx_dip_ih + 1634);

    auto tr_z_yyyyyy_yyzzz = pbuffer.data(idx_dip_ih + 1635);

    auto tr_z_yyyyyy_yzzzz = pbuffer.data(idx_dip_ih + 1636);

    auto tr_z_yyyyyy_zzzzz = pbuffer.data(idx_dip_ih + 1637);

    #pragma omp simd aligned(pa_y, tr_z_yyyy_xxxxx, tr_z_yyyy_xxxxy, tr_z_yyyy_xxxxz, tr_z_yyyy_xxxyy, tr_z_yyyy_xxxyz, tr_z_yyyy_xxxzz, tr_z_yyyy_xxyyy, tr_z_yyyy_xxyyz, tr_z_yyyy_xxyzz, tr_z_yyyy_xxzzz, tr_z_yyyy_xyyyy, tr_z_yyyy_xyyyz, tr_z_yyyy_xyyzz, tr_z_yyyy_xyzzz, tr_z_yyyy_xzzzz, tr_z_yyyy_yyyyy, tr_z_yyyy_yyyyz, tr_z_yyyy_yyyzz, tr_z_yyyy_yyzzz, tr_z_yyyy_yzzzz, tr_z_yyyy_zzzzz, tr_z_yyyyy_xxxx, tr_z_yyyyy_xxxxx, tr_z_yyyyy_xxxxy, tr_z_yyyyy_xxxxz, tr_z_yyyyy_xxxy, tr_z_yyyyy_xxxyy, tr_z_yyyyy_xxxyz, tr_z_yyyyy_xxxz, tr_z_yyyyy_xxxzz, tr_z_yyyyy_xxyy, tr_z_yyyyy_xxyyy, tr_z_yyyyy_xxyyz, tr_z_yyyyy_xxyz, tr_z_yyyyy_xxyzz, tr_z_yyyyy_xxzz, tr_z_yyyyy_xxzzz, tr_z_yyyyy_xyyy, tr_z_yyyyy_xyyyy, tr_z_yyyyy_xyyyz, tr_z_yyyyy_xyyz, tr_z_yyyyy_xyyzz, tr_z_yyyyy_xyzz, tr_z_yyyyy_xyzzz, tr_z_yyyyy_xzzz, tr_z_yyyyy_xzzzz, tr_z_yyyyy_yyyy, tr_z_yyyyy_yyyyy, tr_z_yyyyy_yyyyz, tr_z_yyyyy_yyyz, tr_z_yyyyy_yyyzz, tr_z_yyyyy_yyzz, tr_z_yyyyy_yyzzz, tr_z_yyyyy_yzzz, tr_z_yyyyy_yzzzz, tr_z_yyyyy_zzzz, tr_z_yyyyy_zzzzz, tr_z_yyyyyy_xxxxx, tr_z_yyyyyy_xxxxy, tr_z_yyyyyy_xxxxz, tr_z_yyyyyy_xxxyy, tr_z_yyyyyy_xxxyz, tr_z_yyyyyy_xxxzz, tr_z_yyyyyy_xxyyy, tr_z_yyyyyy_xxyyz, tr_z_yyyyyy_xxyzz, tr_z_yyyyyy_xxzzz, tr_z_yyyyyy_xyyyy, tr_z_yyyyyy_xyyyz, tr_z_yyyyyy_xyyzz, tr_z_yyyyyy_xyzzz, tr_z_yyyyyy_xzzzz, tr_z_yyyyyy_yyyyy, tr_z_yyyyyy_yyyyz, tr_z_yyyyyy_yyyzz, tr_z_yyyyyy_yyzzz, tr_z_yyyyyy_yzzzz, tr_z_yyyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyy_xxxxx[i] = 5.0 * tr_z_yyyy_xxxxx[i] * fe_0 + tr_z_yyyyy_xxxxx[i] * pa_y[i];

        tr_z_yyyyyy_xxxxy[i] = 5.0 * tr_z_yyyy_xxxxy[i] * fe_0 + tr_z_yyyyy_xxxx[i] * fe_0 + tr_z_yyyyy_xxxxy[i] * pa_y[i];

        tr_z_yyyyyy_xxxxz[i] = 5.0 * tr_z_yyyy_xxxxz[i] * fe_0 + tr_z_yyyyy_xxxxz[i] * pa_y[i];

        tr_z_yyyyyy_xxxyy[i] = 5.0 * tr_z_yyyy_xxxyy[i] * fe_0 + 2.0 * tr_z_yyyyy_xxxy[i] * fe_0 + tr_z_yyyyy_xxxyy[i] * pa_y[i];

        tr_z_yyyyyy_xxxyz[i] = 5.0 * tr_z_yyyy_xxxyz[i] * fe_0 + tr_z_yyyyy_xxxz[i] * fe_0 + tr_z_yyyyy_xxxyz[i] * pa_y[i];

        tr_z_yyyyyy_xxxzz[i] = 5.0 * tr_z_yyyy_xxxzz[i] * fe_0 + tr_z_yyyyy_xxxzz[i] * pa_y[i];

        tr_z_yyyyyy_xxyyy[i] = 5.0 * tr_z_yyyy_xxyyy[i] * fe_0 + 3.0 * tr_z_yyyyy_xxyy[i] * fe_0 + tr_z_yyyyy_xxyyy[i] * pa_y[i];

        tr_z_yyyyyy_xxyyz[i] = 5.0 * tr_z_yyyy_xxyyz[i] * fe_0 + 2.0 * tr_z_yyyyy_xxyz[i] * fe_0 + tr_z_yyyyy_xxyyz[i] * pa_y[i];

        tr_z_yyyyyy_xxyzz[i] = 5.0 * tr_z_yyyy_xxyzz[i] * fe_0 + tr_z_yyyyy_xxzz[i] * fe_0 + tr_z_yyyyy_xxyzz[i] * pa_y[i];

        tr_z_yyyyyy_xxzzz[i] = 5.0 * tr_z_yyyy_xxzzz[i] * fe_0 + tr_z_yyyyy_xxzzz[i] * pa_y[i];

        tr_z_yyyyyy_xyyyy[i] = 5.0 * tr_z_yyyy_xyyyy[i] * fe_0 + 4.0 * tr_z_yyyyy_xyyy[i] * fe_0 + tr_z_yyyyy_xyyyy[i] * pa_y[i];

        tr_z_yyyyyy_xyyyz[i] = 5.0 * tr_z_yyyy_xyyyz[i] * fe_0 + 3.0 * tr_z_yyyyy_xyyz[i] * fe_0 + tr_z_yyyyy_xyyyz[i] * pa_y[i];

        tr_z_yyyyyy_xyyzz[i] = 5.0 * tr_z_yyyy_xyyzz[i] * fe_0 + 2.0 * tr_z_yyyyy_xyzz[i] * fe_0 + tr_z_yyyyy_xyyzz[i] * pa_y[i];

        tr_z_yyyyyy_xyzzz[i] = 5.0 * tr_z_yyyy_xyzzz[i] * fe_0 + tr_z_yyyyy_xzzz[i] * fe_0 + tr_z_yyyyy_xyzzz[i] * pa_y[i];

        tr_z_yyyyyy_xzzzz[i] = 5.0 * tr_z_yyyy_xzzzz[i] * fe_0 + tr_z_yyyyy_xzzzz[i] * pa_y[i];

        tr_z_yyyyyy_yyyyy[i] = 5.0 * tr_z_yyyy_yyyyy[i] * fe_0 + 5.0 * tr_z_yyyyy_yyyy[i] * fe_0 + tr_z_yyyyy_yyyyy[i] * pa_y[i];

        tr_z_yyyyyy_yyyyz[i] = 5.0 * tr_z_yyyy_yyyyz[i] * fe_0 + 4.0 * tr_z_yyyyy_yyyz[i] * fe_0 + tr_z_yyyyy_yyyyz[i] * pa_y[i];

        tr_z_yyyyyy_yyyzz[i] = 5.0 * tr_z_yyyy_yyyzz[i] * fe_0 + 3.0 * tr_z_yyyyy_yyzz[i] * fe_0 + tr_z_yyyyy_yyyzz[i] * pa_y[i];

        tr_z_yyyyyy_yyzzz[i] = 5.0 * tr_z_yyyy_yyzzz[i] * fe_0 + 2.0 * tr_z_yyyyy_yzzz[i] * fe_0 + tr_z_yyyyy_yyzzz[i] * pa_y[i];

        tr_z_yyyyyy_yzzzz[i] = 5.0 * tr_z_yyyy_yzzzz[i] * fe_0 + tr_z_yyyyy_zzzz[i] * fe_0 + tr_z_yyyyy_yzzzz[i] * pa_y[i];

        tr_z_yyyyyy_zzzzz[i] = 5.0 * tr_z_yyyy_zzzzz[i] * fe_0 + tr_z_yyyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 1638-1659 components of targeted buffer : IH

    auto tr_z_yyyyyz_xxxxx = pbuffer.data(idx_dip_ih + 1638);

    auto tr_z_yyyyyz_xxxxy = pbuffer.data(idx_dip_ih + 1639);

    auto tr_z_yyyyyz_xxxxz = pbuffer.data(idx_dip_ih + 1640);

    auto tr_z_yyyyyz_xxxyy = pbuffer.data(idx_dip_ih + 1641);

    auto tr_z_yyyyyz_xxxyz = pbuffer.data(idx_dip_ih + 1642);

    auto tr_z_yyyyyz_xxxzz = pbuffer.data(idx_dip_ih + 1643);

    auto tr_z_yyyyyz_xxyyy = pbuffer.data(idx_dip_ih + 1644);

    auto tr_z_yyyyyz_xxyyz = pbuffer.data(idx_dip_ih + 1645);

    auto tr_z_yyyyyz_xxyzz = pbuffer.data(idx_dip_ih + 1646);

    auto tr_z_yyyyyz_xxzzz = pbuffer.data(idx_dip_ih + 1647);

    auto tr_z_yyyyyz_xyyyy = pbuffer.data(idx_dip_ih + 1648);

    auto tr_z_yyyyyz_xyyyz = pbuffer.data(idx_dip_ih + 1649);

    auto tr_z_yyyyyz_xyyzz = pbuffer.data(idx_dip_ih + 1650);

    auto tr_z_yyyyyz_xyzzz = pbuffer.data(idx_dip_ih + 1651);

    auto tr_z_yyyyyz_xzzzz = pbuffer.data(idx_dip_ih + 1652);

    auto tr_z_yyyyyz_yyyyy = pbuffer.data(idx_dip_ih + 1653);

    auto tr_z_yyyyyz_yyyyz = pbuffer.data(idx_dip_ih + 1654);

    auto tr_z_yyyyyz_yyyzz = pbuffer.data(idx_dip_ih + 1655);

    auto tr_z_yyyyyz_yyzzz = pbuffer.data(idx_dip_ih + 1656);

    auto tr_z_yyyyyz_yzzzz = pbuffer.data(idx_dip_ih + 1657);

    auto tr_z_yyyyyz_zzzzz = pbuffer.data(idx_dip_ih + 1658);

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yyyyy_xxxxy, tr_z_yyyyy_xxxyy, tr_z_yyyyy_xxyyy, tr_z_yyyyy_xyyyy, tr_z_yyyyy_yyyyy, tr_z_yyyyyz_xxxxx, tr_z_yyyyyz_xxxxy, tr_z_yyyyyz_xxxxz, tr_z_yyyyyz_xxxyy, tr_z_yyyyyz_xxxyz, tr_z_yyyyyz_xxxzz, tr_z_yyyyyz_xxyyy, tr_z_yyyyyz_xxyyz, tr_z_yyyyyz_xxyzz, tr_z_yyyyyz_xxzzz, tr_z_yyyyyz_xyyyy, tr_z_yyyyyz_xyyyz, tr_z_yyyyyz_xyyzz, tr_z_yyyyyz_xyzzz, tr_z_yyyyyz_xzzzz, tr_z_yyyyyz_yyyyy, tr_z_yyyyyz_yyyyz, tr_z_yyyyyz_yyyzz, tr_z_yyyyyz_yyzzz, tr_z_yyyyyz_yzzzz, tr_z_yyyyyz_zzzzz, tr_z_yyyyz_xxxxx, tr_z_yyyyz_xxxxz, tr_z_yyyyz_xxxyz, tr_z_yyyyz_xxxz, tr_z_yyyyz_xxxzz, tr_z_yyyyz_xxyyz, tr_z_yyyyz_xxyz, tr_z_yyyyz_xxyzz, tr_z_yyyyz_xxzz, tr_z_yyyyz_xxzzz, tr_z_yyyyz_xyyyz, tr_z_yyyyz_xyyz, tr_z_yyyyz_xyyzz, tr_z_yyyyz_xyzz, tr_z_yyyyz_xyzzz, tr_z_yyyyz_xzzz, tr_z_yyyyz_xzzzz, tr_z_yyyyz_yyyyz, tr_z_yyyyz_yyyz, tr_z_yyyyz_yyyzz, tr_z_yyyyz_yyzz, tr_z_yyyyz_yyzzz, tr_z_yyyyz_yzzz, tr_z_yyyyz_yzzzz, tr_z_yyyyz_zzzz, tr_z_yyyyz_zzzzz, tr_z_yyyz_xxxxx, tr_z_yyyz_xxxxz, tr_z_yyyz_xxxyz, tr_z_yyyz_xxxzz, tr_z_yyyz_xxyyz, tr_z_yyyz_xxyzz, tr_z_yyyz_xxzzz, tr_z_yyyz_xyyyz, tr_z_yyyz_xyyzz, tr_z_yyyz_xyzzz, tr_z_yyyz_xzzzz, tr_z_yyyz_yyyyz, tr_z_yyyz_yyyzz, tr_z_yyyz_yyzzz, tr_z_yyyz_yzzzz, tr_z_yyyz_zzzzz, ts_yyyyy_xxxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxyyy, ts_yyyyy_xyyyy, ts_yyyyy_yyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyz_xxxxx[i] = 4.0 * tr_z_yyyz_xxxxx[i] * fe_0 + tr_z_yyyyz_xxxxx[i] * pa_y[i];

        tr_z_yyyyyz_xxxxy[i] = ts_yyyyy_xxxxy[i] * fe_0 + tr_z_yyyyy_xxxxy[i] * pa_z[i];

        tr_z_yyyyyz_xxxxz[i] = 4.0 * tr_z_yyyz_xxxxz[i] * fe_0 + tr_z_yyyyz_xxxxz[i] * pa_y[i];

        tr_z_yyyyyz_xxxyy[i] = ts_yyyyy_xxxyy[i] * fe_0 + tr_z_yyyyy_xxxyy[i] * pa_z[i];

        tr_z_yyyyyz_xxxyz[i] = 4.0 * tr_z_yyyz_xxxyz[i] * fe_0 + tr_z_yyyyz_xxxz[i] * fe_0 + tr_z_yyyyz_xxxyz[i] * pa_y[i];

        tr_z_yyyyyz_xxxzz[i] = 4.0 * tr_z_yyyz_xxxzz[i] * fe_0 + tr_z_yyyyz_xxxzz[i] * pa_y[i];

        tr_z_yyyyyz_xxyyy[i] = ts_yyyyy_xxyyy[i] * fe_0 + tr_z_yyyyy_xxyyy[i] * pa_z[i];

        tr_z_yyyyyz_xxyyz[i] = 4.0 * tr_z_yyyz_xxyyz[i] * fe_0 + 2.0 * tr_z_yyyyz_xxyz[i] * fe_0 + tr_z_yyyyz_xxyyz[i] * pa_y[i];

        tr_z_yyyyyz_xxyzz[i] = 4.0 * tr_z_yyyz_xxyzz[i] * fe_0 + tr_z_yyyyz_xxzz[i] * fe_0 + tr_z_yyyyz_xxyzz[i] * pa_y[i];

        tr_z_yyyyyz_xxzzz[i] = 4.0 * tr_z_yyyz_xxzzz[i] * fe_0 + tr_z_yyyyz_xxzzz[i] * pa_y[i];

        tr_z_yyyyyz_xyyyy[i] = ts_yyyyy_xyyyy[i] * fe_0 + tr_z_yyyyy_xyyyy[i] * pa_z[i];

        tr_z_yyyyyz_xyyyz[i] = 4.0 * tr_z_yyyz_xyyyz[i] * fe_0 + 3.0 * tr_z_yyyyz_xyyz[i] * fe_0 + tr_z_yyyyz_xyyyz[i] * pa_y[i];

        tr_z_yyyyyz_xyyzz[i] = 4.0 * tr_z_yyyz_xyyzz[i] * fe_0 + 2.0 * tr_z_yyyyz_xyzz[i] * fe_0 + tr_z_yyyyz_xyyzz[i] * pa_y[i];

        tr_z_yyyyyz_xyzzz[i] = 4.0 * tr_z_yyyz_xyzzz[i] * fe_0 + tr_z_yyyyz_xzzz[i] * fe_0 + tr_z_yyyyz_xyzzz[i] * pa_y[i];

        tr_z_yyyyyz_xzzzz[i] = 4.0 * tr_z_yyyz_xzzzz[i] * fe_0 + tr_z_yyyyz_xzzzz[i] * pa_y[i];

        tr_z_yyyyyz_yyyyy[i] = ts_yyyyy_yyyyy[i] * fe_0 + tr_z_yyyyy_yyyyy[i] * pa_z[i];

        tr_z_yyyyyz_yyyyz[i] = 4.0 * tr_z_yyyz_yyyyz[i] * fe_0 + 4.0 * tr_z_yyyyz_yyyz[i] * fe_0 + tr_z_yyyyz_yyyyz[i] * pa_y[i];

        tr_z_yyyyyz_yyyzz[i] = 4.0 * tr_z_yyyz_yyyzz[i] * fe_0 + 3.0 * tr_z_yyyyz_yyzz[i] * fe_0 + tr_z_yyyyz_yyyzz[i] * pa_y[i];

        tr_z_yyyyyz_yyzzz[i] = 4.0 * tr_z_yyyz_yyzzz[i] * fe_0 + 2.0 * tr_z_yyyyz_yzzz[i] * fe_0 + tr_z_yyyyz_yyzzz[i] * pa_y[i];

        tr_z_yyyyyz_yzzzz[i] = 4.0 * tr_z_yyyz_yzzzz[i] * fe_0 + tr_z_yyyyz_zzzz[i] * fe_0 + tr_z_yyyyz_yzzzz[i] * pa_y[i];

        tr_z_yyyyyz_zzzzz[i] = 4.0 * tr_z_yyyz_zzzzz[i] * fe_0 + tr_z_yyyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 1659-1680 components of targeted buffer : IH

    auto tr_z_yyyyzz_xxxxx = pbuffer.data(idx_dip_ih + 1659);

    auto tr_z_yyyyzz_xxxxy = pbuffer.data(idx_dip_ih + 1660);

    auto tr_z_yyyyzz_xxxxz = pbuffer.data(idx_dip_ih + 1661);

    auto tr_z_yyyyzz_xxxyy = pbuffer.data(idx_dip_ih + 1662);

    auto tr_z_yyyyzz_xxxyz = pbuffer.data(idx_dip_ih + 1663);

    auto tr_z_yyyyzz_xxxzz = pbuffer.data(idx_dip_ih + 1664);

    auto tr_z_yyyyzz_xxyyy = pbuffer.data(idx_dip_ih + 1665);

    auto tr_z_yyyyzz_xxyyz = pbuffer.data(idx_dip_ih + 1666);

    auto tr_z_yyyyzz_xxyzz = pbuffer.data(idx_dip_ih + 1667);

    auto tr_z_yyyyzz_xxzzz = pbuffer.data(idx_dip_ih + 1668);

    auto tr_z_yyyyzz_xyyyy = pbuffer.data(idx_dip_ih + 1669);

    auto tr_z_yyyyzz_xyyyz = pbuffer.data(idx_dip_ih + 1670);

    auto tr_z_yyyyzz_xyyzz = pbuffer.data(idx_dip_ih + 1671);

    auto tr_z_yyyyzz_xyzzz = pbuffer.data(idx_dip_ih + 1672);

    auto tr_z_yyyyzz_xzzzz = pbuffer.data(idx_dip_ih + 1673);

    auto tr_z_yyyyzz_yyyyy = pbuffer.data(idx_dip_ih + 1674);

    auto tr_z_yyyyzz_yyyyz = pbuffer.data(idx_dip_ih + 1675);

    auto tr_z_yyyyzz_yyyzz = pbuffer.data(idx_dip_ih + 1676);

    auto tr_z_yyyyzz_yyzzz = pbuffer.data(idx_dip_ih + 1677);

    auto tr_z_yyyyzz_yzzzz = pbuffer.data(idx_dip_ih + 1678);

    auto tr_z_yyyyzz_zzzzz = pbuffer.data(idx_dip_ih + 1679);

    #pragma omp simd aligned(pa_y, tr_z_yyyyzz_xxxxx, tr_z_yyyyzz_xxxxy, tr_z_yyyyzz_xxxxz, tr_z_yyyyzz_xxxyy, tr_z_yyyyzz_xxxyz, tr_z_yyyyzz_xxxzz, tr_z_yyyyzz_xxyyy, tr_z_yyyyzz_xxyyz, tr_z_yyyyzz_xxyzz, tr_z_yyyyzz_xxzzz, tr_z_yyyyzz_xyyyy, tr_z_yyyyzz_xyyyz, tr_z_yyyyzz_xyyzz, tr_z_yyyyzz_xyzzz, tr_z_yyyyzz_xzzzz, tr_z_yyyyzz_yyyyy, tr_z_yyyyzz_yyyyz, tr_z_yyyyzz_yyyzz, tr_z_yyyyzz_yyzzz, tr_z_yyyyzz_yzzzz, tr_z_yyyyzz_zzzzz, tr_z_yyyzz_xxxx, tr_z_yyyzz_xxxxx, tr_z_yyyzz_xxxxy, tr_z_yyyzz_xxxxz, tr_z_yyyzz_xxxy, tr_z_yyyzz_xxxyy, tr_z_yyyzz_xxxyz, tr_z_yyyzz_xxxz, tr_z_yyyzz_xxxzz, tr_z_yyyzz_xxyy, tr_z_yyyzz_xxyyy, tr_z_yyyzz_xxyyz, tr_z_yyyzz_xxyz, tr_z_yyyzz_xxyzz, tr_z_yyyzz_xxzz, tr_z_yyyzz_xxzzz, tr_z_yyyzz_xyyy, tr_z_yyyzz_xyyyy, tr_z_yyyzz_xyyyz, tr_z_yyyzz_xyyz, tr_z_yyyzz_xyyzz, tr_z_yyyzz_xyzz, tr_z_yyyzz_xyzzz, tr_z_yyyzz_xzzz, tr_z_yyyzz_xzzzz, tr_z_yyyzz_yyyy, tr_z_yyyzz_yyyyy, tr_z_yyyzz_yyyyz, tr_z_yyyzz_yyyz, tr_z_yyyzz_yyyzz, tr_z_yyyzz_yyzz, tr_z_yyyzz_yyzzz, tr_z_yyyzz_yzzz, tr_z_yyyzz_yzzzz, tr_z_yyyzz_zzzz, tr_z_yyyzz_zzzzz, tr_z_yyzz_xxxxx, tr_z_yyzz_xxxxy, tr_z_yyzz_xxxxz, tr_z_yyzz_xxxyy, tr_z_yyzz_xxxyz, tr_z_yyzz_xxxzz, tr_z_yyzz_xxyyy, tr_z_yyzz_xxyyz, tr_z_yyzz_xxyzz, tr_z_yyzz_xxzzz, tr_z_yyzz_xyyyy, tr_z_yyzz_xyyyz, tr_z_yyzz_xyyzz, tr_z_yyzz_xyzzz, tr_z_yyzz_xzzzz, tr_z_yyzz_yyyyy, tr_z_yyzz_yyyyz, tr_z_yyzz_yyyzz, tr_z_yyzz_yyzzz, tr_z_yyzz_yzzzz, tr_z_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyzz_xxxxx[i] = 3.0 * tr_z_yyzz_xxxxx[i] * fe_0 + tr_z_yyyzz_xxxxx[i] * pa_y[i];

        tr_z_yyyyzz_xxxxy[i] = 3.0 * tr_z_yyzz_xxxxy[i] * fe_0 + tr_z_yyyzz_xxxx[i] * fe_0 + tr_z_yyyzz_xxxxy[i] * pa_y[i];

        tr_z_yyyyzz_xxxxz[i] = 3.0 * tr_z_yyzz_xxxxz[i] * fe_0 + tr_z_yyyzz_xxxxz[i] * pa_y[i];

        tr_z_yyyyzz_xxxyy[i] = 3.0 * tr_z_yyzz_xxxyy[i] * fe_0 + 2.0 * tr_z_yyyzz_xxxy[i] * fe_0 + tr_z_yyyzz_xxxyy[i] * pa_y[i];

        tr_z_yyyyzz_xxxyz[i] = 3.0 * tr_z_yyzz_xxxyz[i] * fe_0 + tr_z_yyyzz_xxxz[i] * fe_0 + tr_z_yyyzz_xxxyz[i] * pa_y[i];

        tr_z_yyyyzz_xxxzz[i] = 3.0 * tr_z_yyzz_xxxzz[i] * fe_0 + tr_z_yyyzz_xxxzz[i] * pa_y[i];

        tr_z_yyyyzz_xxyyy[i] = 3.0 * tr_z_yyzz_xxyyy[i] * fe_0 + 3.0 * tr_z_yyyzz_xxyy[i] * fe_0 + tr_z_yyyzz_xxyyy[i] * pa_y[i];

        tr_z_yyyyzz_xxyyz[i] = 3.0 * tr_z_yyzz_xxyyz[i] * fe_0 + 2.0 * tr_z_yyyzz_xxyz[i] * fe_0 + tr_z_yyyzz_xxyyz[i] * pa_y[i];

        tr_z_yyyyzz_xxyzz[i] = 3.0 * tr_z_yyzz_xxyzz[i] * fe_0 + tr_z_yyyzz_xxzz[i] * fe_0 + tr_z_yyyzz_xxyzz[i] * pa_y[i];

        tr_z_yyyyzz_xxzzz[i] = 3.0 * tr_z_yyzz_xxzzz[i] * fe_0 + tr_z_yyyzz_xxzzz[i] * pa_y[i];

        tr_z_yyyyzz_xyyyy[i] = 3.0 * tr_z_yyzz_xyyyy[i] * fe_0 + 4.0 * tr_z_yyyzz_xyyy[i] * fe_0 + tr_z_yyyzz_xyyyy[i] * pa_y[i];

        tr_z_yyyyzz_xyyyz[i] = 3.0 * tr_z_yyzz_xyyyz[i] * fe_0 + 3.0 * tr_z_yyyzz_xyyz[i] * fe_0 + tr_z_yyyzz_xyyyz[i] * pa_y[i];

        tr_z_yyyyzz_xyyzz[i] = 3.0 * tr_z_yyzz_xyyzz[i] * fe_0 + 2.0 * tr_z_yyyzz_xyzz[i] * fe_0 + tr_z_yyyzz_xyyzz[i] * pa_y[i];

        tr_z_yyyyzz_xyzzz[i] = 3.0 * tr_z_yyzz_xyzzz[i] * fe_0 + tr_z_yyyzz_xzzz[i] * fe_0 + tr_z_yyyzz_xyzzz[i] * pa_y[i];

        tr_z_yyyyzz_xzzzz[i] = 3.0 * tr_z_yyzz_xzzzz[i] * fe_0 + tr_z_yyyzz_xzzzz[i] * pa_y[i];

        tr_z_yyyyzz_yyyyy[i] = 3.0 * tr_z_yyzz_yyyyy[i] * fe_0 + 5.0 * tr_z_yyyzz_yyyy[i] * fe_0 + tr_z_yyyzz_yyyyy[i] * pa_y[i];

        tr_z_yyyyzz_yyyyz[i] = 3.0 * tr_z_yyzz_yyyyz[i] * fe_0 + 4.0 * tr_z_yyyzz_yyyz[i] * fe_0 + tr_z_yyyzz_yyyyz[i] * pa_y[i];

        tr_z_yyyyzz_yyyzz[i] = 3.0 * tr_z_yyzz_yyyzz[i] * fe_0 + 3.0 * tr_z_yyyzz_yyzz[i] * fe_0 + tr_z_yyyzz_yyyzz[i] * pa_y[i];

        tr_z_yyyyzz_yyzzz[i] = 3.0 * tr_z_yyzz_yyzzz[i] * fe_0 + 2.0 * tr_z_yyyzz_yzzz[i] * fe_0 + tr_z_yyyzz_yyzzz[i] * pa_y[i];

        tr_z_yyyyzz_yzzzz[i] = 3.0 * tr_z_yyzz_yzzzz[i] * fe_0 + tr_z_yyyzz_zzzz[i] * fe_0 + tr_z_yyyzz_yzzzz[i] * pa_y[i];

        tr_z_yyyyzz_zzzzz[i] = 3.0 * tr_z_yyzz_zzzzz[i] * fe_0 + tr_z_yyyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1680-1701 components of targeted buffer : IH

    auto tr_z_yyyzzz_xxxxx = pbuffer.data(idx_dip_ih + 1680);

    auto tr_z_yyyzzz_xxxxy = pbuffer.data(idx_dip_ih + 1681);

    auto tr_z_yyyzzz_xxxxz = pbuffer.data(idx_dip_ih + 1682);

    auto tr_z_yyyzzz_xxxyy = pbuffer.data(idx_dip_ih + 1683);

    auto tr_z_yyyzzz_xxxyz = pbuffer.data(idx_dip_ih + 1684);

    auto tr_z_yyyzzz_xxxzz = pbuffer.data(idx_dip_ih + 1685);

    auto tr_z_yyyzzz_xxyyy = pbuffer.data(idx_dip_ih + 1686);

    auto tr_z_yyyzzz_xxyyz = pbuffer.data(idx_dip_ih + 1687);

    auto tr_z_yyyzzz_xxyzz = pbuffer.data(idx_dip_ih + 1688);

    auto tr_z_yyyzzz_xxzzz = pbuffer.data(idx_dip_ih + 1689);

    auto tr_z_yyyzzz_xyyyy = pbuffer.data(idx_dip_ih + 1690);

    auto tr_z_yyyzzz_xyyyz = pbuffer.data(idx_dip_ih + 1691);

    auto tr_z_yyyzzz_xyyzz = pbuffer.data(idx_dip_ih + 1692);

    auto tr_z_yyyzzz_xyzzz = pbuffer.data(idx_dip_ih + 1693);

    auto tr_z_yyyzzz_xzzzz = pbuffer.data(idx_dip_ih + 1694);

    auto tr_z_yyyzzz_yyyyy = pbuffer.data(idx_dip_ih + 1695);

    auto tr_z_yyyzzz_yyyyz = pbuffer.data(idx_dip_ih + 1696);

    auto tr_z_yyyzzz_yyyzz = pbuffer.data(idx_dip_ih + 1697);

    auto tr_z_yyyzzz_yyzzz = pbuffer.data(idx_dip_ih + 1698);

    auto tr_z_yyyzzz_yzzzz = pbuffer.data(idx_dip_ih + 1699);

    auto tr_z_yyyzzz_zzzzz = pbuffer.data(idx_dip_ih + 1700);

    #pragma omp simd aligned(pa_y, tr_z_yyyzzz_xxxxx, tr_z_yyyzzz_xxxxy, tr_z_yyyzzz_xxxxz, tr_z_yyyzzz_xxxyy, tr_z_yyyzzz_xxxyz, tr_z_yyyzzz_xxxzz, tr_z_yyyzzz_xxyyy, tr_z_yyyzzz_xxyyz, tr_z_yyyzzz_xxyzz, tr_z_yyyzzz_xxzzz, tr_z_yyyzzz_xyyyy, tr_z_yyyzzz_xyyyz, tr_z_yyyzzz_xyyzz, tr_z_yyyzzz_xyzzz, tr_z_yyyzzz_xzzzz, tr_z_yyyzzz_yyyyy, tr_z_yyyzzz_yyyyz, tr_z_yyyzzz_yyyzz, tr_z_yyyzzz_yyzzz, tr_z_yyyzzz_yzzzz, tr_z_yyyzzz_zzzzz, tr_z_yyzzz_xxxx, tr_z_yyzzz_xxxxx, tr_z_yyzzz_xxxxy, tr_z_yyzzz_xxxxz, tr_z_yyzzz_xxxy, tr_z_yyzzz_xxxyy, tr_z_yyzzz_xxxyz, tr_z_yyzzz_xxxz, tr_z_yyzzz_xxxzz, tr_z_yyzzz_xxyy, tr_z_yyzzz_xxyyy, tr_z_yyzzz_xxyyz, tr_z_yyzzz_xxyz, tr_z_yyzzz_xxyzz, tr_z_yyzzz_xxzz, tr_z_yyzzz_xxzzz, tr_z_yyzzz_xyyy, tr_z_yyzzz_xyyyy, tr_z_yyzzz_xyyyz, tr_z_yyzzz_xyyz, tr_z_yyzzz_xyyzz, tr_z_yyzzz_xyzz, tr_z_yyzzz_xyzzz, tr_z_yyzzz_xzzz, tr_z_yyzzz_xzzzz, tr_z_yyzzz_yyyy, tr_z_yyzzz_yyyyy, tr_z_yyzzz_yyyyz, tr_z_yyzzz_yyyz, tr_z_yyzzz_yyyzz, tr_z_yyzzz_yyzz, tr_z_yyzzz_yyzzz, tr_z_yyzzz_yzzz, tr_z_yyzzz_yzzzz, tr_z_yyzzz_zzzz, tr_z_yyzzz_zzzzz, tr_z_yzzz_xxxxx, tr_z_yzzz_xxxxy, tr_z_yzzz_xxxxz, tr_z_yzzz_xxxyy, tr_z_yzzz_xxxyz, tr_z_yzzz_xxxzz, tr_z_yzzz_xxyyy, tr_z_yzzz_xxyyz, tr_z_yzzz_xxyzz, tr_z_yzzz_xxzzz, tr_z_yzzz_xyyyy, tr_z_yzzz_xyyyz, tr_z_yzzz_xyyzz, tr_z_yzzz_xyzzz, tr_z_yzzz_xzzzz, tr_z_yzzz_yyyyy, tr_z_yzzz_yyyyz, tr_z_yzzz_yyyzz, tr_z_yzzz_yyzzz, tr_z_yzzz_yzzzz, tr_z_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzzz_xxxxx[i] = 2.0 * tr_z_yzzz_xxxxx[i] * fe_0 + tr_z_yyzzz_xxxxx[i] * pa_y[i];

        tr_z_yyyzzz_xxxxy[i] = 2.0 * tr_z_yzzz_xxxxy[i] * fe_0 + tr_z_yyzzz_xxxx[i] * fe_0 + tr_z_yyzzz_xxxxy[i] * pa_y[i];

        tr_z_yyyzzz_xxxxz[i] = 2.0 * tr_z_yzzz_xxxxz[i] * fe_0 + tr_z_yyzzz_xxxxz[i] * pa_y[i];

        tr_z_yyyzzz_xxxyy[i] = 2.0 * tr_z_yzzz_xxxyy[i] * fe_0 + 2.0 * tr_z_yyzzz_xxxy[i] * fe_0 + tr_z_yyzzz_xxxyy[i] * pa_y[i];

        tr_z_yyyzzz_xxxyz[i] = 2.0 * tr_z_yzzz_xxxyz[i] * fe_0 + tr_z_yyzzz_xxxz[i] * fe_0 + tr_z_yyzzz_xxxyz[i] * pa_y[i];

        tr_z_yyyzzz_xxxzz[i] = 2.0 * tr_z_yzzz_xxxzz[i] * fe_0 + tr_z_yyzzz_xxxzz[i] * pa_y[i];

        tr_z_yyyzzz_xxyyy[i] = 2.0 * tr_z_yzzz_xxyyy[i] * fe_0 + 3.0 * tr_z_yyzzz_xxyy[i] * fe_0 + tr_z_yyzzz_xxyyy[i] * pa_y[i];

        tr_z_yyyzzz_xxyyz[i] = 2.0 * tr_z_yzzz_xxyyz[i] * fe_0 + 2.0 * tr_z_yyzzz_xxyz[i] * fe_0 + tr_z_yyzzz_xxyyz[i] * pa_y[i];

        tr_z_yyyzzz_xxyzz[i] = 2.0 * tr_z_yzzz_xxyzz[i] * fe_0 + tr_z_yyzzz_xxzz[i] * fe_0 + tr_z_yyzzz_xxyzz[i] * pa_y[i];

        tr_z_yyyzzz_xxzzz[i] = 2.0 * tr_z_yzzz_xxzzz[i] * fe_0 + tr_z_yyzzz_xxzzz[i] * pa_y[i];

        tr_z_yyyzzz_xyyyy[i] = 2.0 * tr_z_yzzz_xyyyy[i] * fe_0 + 4.0 * tr_z_yyzzz_xyyy[i] * fe_0 + tr_z_yyzzz_xyyyy[i] * pa_y[i];

        tr_z_yyyzzz_xyyyz[i] = 2.0 * tr_z_yzzz_xyyyz[i] * fe_0 + 3.0 * tr_z_yyzzz_xyyz[i] * fe_0 + tr_z_yyzzz_xyyyz[i] * pa_y[i];

        tr_z_yyyzzz_xyyzz[i] = 2.0 * tr_z_yzzz_xyyzz[i] * fe_0 + 2.0 * tr_z_yyzzz_xyzz[i] * fe_0 + tr_z_yyzzz_xyyzz[i] * pa_y[i];

        tr_z_yyyzzz_xyzzz[i] = 2.0 * tr_z_yzzz_xyzzz[i] * fe_0 + tr_z_yyzzz_xzzz[i] * fe_0 + tr_z_yyzzz_xyzzz[i] * pa_y[i];

        tr_z_yyyzzz_xzzzz[i] = 2.0 * tr_z_yzzz_xzzzz[i] * fe_0 + tr_z_yyzzz_xzzzz[i] * pa_y[i];

        tr_z_yyyzzz_yyyyy[i] = 2.0 * tr_z_yzzz_yyyyy[i] * fe_0 + 5.0 * tr_z_yyzzz_yyyy[i] * fe_0 + tr_z_yyzzz_yyyyy[i] * pa_y[i];

        tr_z_yyyzzz_yyyyz[i] = 2.0 * tr_z_yzzz_yyyyz[i] * fe_0 + 4.0 * tr_z_yyzzz_yyyz[i] * fe_0 + tr_z_yyzzz_yyyyz[i] * pa_y[i];

        tr_z_yyyzzz_yyyzz[i] = 2.0 * tr_z_yzzz_yyyzz[i] * fe_0 + 3.0 * tr_z_yyzzz_yyzz[i] * fe_0 + tr_z_yyzzz_yyyzz[i] * pa_y[i];

        tr_z_yyyzzz_yyzzz[i] = 2.0 * tr_z_yzzz_yyzzz[i] * fe_0 + 2.0 * tr_z_yyzzz_yzzz[i] * fe_0 + tr_z_yyzzz_yyzzz[i] * pa_y[i];

        tr_z_yyyzzz_yzzzz[i] = 2.0 * tr_z_yzzz_yzzzz[i] * fe_0 + tr_z_yyzzz_zzzz[i] * fe_0 + tr_z_yyzzz_yzzzz[i] * pa_y[i];

        tr_z_yyyzzz_zzzzz[i] = 2.0 * tr_z_yzzz_zzzzz[i] * fe_0 + tr_z_yyzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1701-1722 components of targeted buffer : IH

    auto tr_z_yyzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1701);

    auto tr_z_yyzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1702);

    auto tr_z_yyzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1703);

    auto tr_z_yyzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1704);

    auto tr_z_yyzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1705);

    auto tr_z_yyzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1706);

    auto tr_z_yyzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1707);

    auto tr_z_yyzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1708);

    auto tr_z_yyzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1709);

    auto tr_z_yyzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1710);

    auto tr_z_yyzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1711);

    auto tr_z_yyzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1712);

    auto tr_z_yyzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1713);

    auto tr_z_yyzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1714);

    auto tr_z_yyzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1715);

    auto tr_z_yyzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1716);

    auto tr_z_yyzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1717);

    auto tr_z_yyzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1718);

    auto tr_z_yyzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1719);

    auto tr_z_yyzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1720);

    auto tr_z_yyzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1721);

    #pragma omp simd aligned(pa_y, tr_z_yyzzzz_xxxxx, tr_z_yyzzzz_xxxxy, tr_z_yyzzzz_xxxxz, tr_z_yyzzzz_xxxyy, tr_z_yyzzzz_xxxyz, tr_z_yyzzzz_xxxzz, tr_z_yyzzzz_xxyyy, tr_z_yyzzzz_xxyyz, tr_z_yyzzzz_xxyzz, tr_z_yyzzzz_xxzzz, tr_z_yyzzzz_xyyyy, tr_z_yyzzzz_xyyyz, tr_z_yyzzzz_xyyzz, tr_z_yyzzzz_xyzzz, tr_z_yyzzzz_xzzzz, tr_z_yyzzzz_yyyyy, tr_z_yyzzzz_yyyyz, tr_z_yyzzzz_yyyzz, tr_z_yyzzzz_yyzzz, tr_z_yyzzzz_yzzzz, tr_z_yyzzzz_zzzzz, tr_z_yzzzz_xxxx, tr_z_yzzzz_xxxxx, tr_z_yzzzz_xxxxy, tr_z_yzzzz_xxxxz, tr_z_yzzzz_xxxy, tr_z_yzzzz_xxxyy, tr_z_yzzzz_xxxyz, tr_z_yzzzz_xxxz, tr_z_yzzzz_xxxzz, tr_z_yzzzz_xxyy, tr_z_yzzzz_xxyyy, tr_z_yzzzz_xxyyz, tr_z_yzzzz_xxyz, tr_z_yzzzz_xxyzz, tr_z_yzzzz_xxzz, tr_z_yzzzz_xxzzz, tr_z_yzzzz_xyyy, tr_z_yzzzz_xyyyy, tr_z_yzzzz_xyyyz, tr_z_yzzzz_xyyz, tr_z_yzzzz_xyyzz, tr_z_yzzzz_xyzz, tr_z_yzzzz_xyzzz, tr_z_yzzzz_xzzz, tr_z_yzzzz_xzzzz, tr_z_yzzzz_yyyy, tr_z_yzzzz_yyyyy, tr_z_yzzzz_yyyyz, tr_z_yzzzz_yyyz, tr_z_yzzzz_yyyzz, tr_z_yzzzz_yyzz, tr_z_yzzzz_yyzzz, tr_z_yzzzz_yzzz, tr_z_yzzzz_yzzzz, tr_z_yzzzz_zzzz, tr_z_yzzzz_zzzzz, tr_z_zzzz_xxxxx, tr_z_zzzz_xxxxy, tr_z_zzzz_xxxxz, tr_z_zzzz_xxxyy, tr_z_zzzz_xxxyz, tr_z_zzzz_xxxzz, tr_z_zzzz_xxyyy, tr_z_zzzz_xxyyz, tr_z_zzzz_xxyzz, tr_z_zzzz_xxzzz, tr_z_zzzz_xyyyy, tr_z_zzzz_xyyyz, tr_z_zzzz_xyyzz, tr_z_zzzz_xyzzz, tr_z_zzzz_xzzzz, tr_z_zzzz_yyyyy, tr_z_zzzz_yyyyz, tr_z_zzzz_yyyzz, tr_z_zzzz_yyzzz, tr_z_zzzz_yzzzz, tr_z_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzzz_xxxxx[i] = tr_z_zzzz_xxxxx[i] * fe_0 + tr_z_yzzzz_xxxxx[i] * pa_y[i];

        tr_z_yyzzzz_xxxxy[i] = tr_z_zzzz_xxxxy[i] * fe_0 + tr_z_yzzzz_xxxx[i] * fe_0 + tr_z_yzzzz_xxxxy[i] * pa_y[i];

        tr_z_yyzzzz_xxxxz[i] = tr_z_zzzz_xxxxz[i] * fe_0 + tr_z_yzzzz_xxxxz[i] * pa_y[i];

        tr_z_yyzzzz_xxxyy[i] = tr_z_zzzz_xxxyy[i] * fe_0 + 2.0 * tr_z_yzzzz_xxxy[i] * fe_0 + tr_z_yzzzz_xxxyy[i] * pa_y[i];

        tr_z_yyzzzz_xxxyz[i] = tr_z_zzzz_xxxyz[i] * fe_0 + tr_z_yzzzz_xxxz[i] * fe_0 + tr_z_yzzzz_xxxyz[i] * pa_y[i];

        tr_z_yyzzzz_xxxzz[i] = tr_z_zzzz_xxxzz[i] * fe_0 + tr_z_yzzzz_xxxzz[i] * pa_y[i];

        tr_z_yyzzzz_xxyyy[i] = tr_z_zzzz_xxyyy[i] * fe_0 + 3.0 * tr_z_yzzzz_xxyy[i] * fe_0 + tr_z_yzzzz_xxyyy[i] * pa_y[i];

        tr_z_yyzzzz_xxyyz[i] = tr_z_zzzz_xxyyz[i] * fe_0 + 2.0 * tr_z_yzzzz_xxyz[i] * fe_0 + tr_z_yzzzz_xxyyz[i] * pa_y[i];

        tr_z_yyzzzz_xxyzz[i] = tr_z_zzzz_xxyzz[i] * fe_0 + tr_z_yzzzz_xxzz[i] * fe_0 + tr_z_yzzzz_xxyzz[i] * pa_y[i];

        tr_z_yyzzzz_xxzzz[i] = tr_z_zzzz_xxzzz[i] * fe_0 + tr_z_yzzzz_xxzzz[i] * pa_y[i];

        tr_z_yyzzzz_xyyyy[i] = tr_z_zzzz_xyyyy[i] * fe_0 + 4.0 * tr_z_yzzzz_xyyy[i] * fe_0 + tr_z_yzzzz_xyyyy[i] * pa_y[i];

        tr_z_yyzzzz_xyyyz[i] = tr_z_zzzz_xyyyz[i] * fe_0 + 3.0 * tr_z_yzzzz_xyyz[i] * fe_0 + tr_z_yzzzz_xyyyz[i] * pa_y[i];

        tr_z_yyzzzz_xyyzz[i] = tr_z_zzzz_xyyzz[i] * fe_0 + 2.0 * tr_z_yzzzz_xyzz[i] * fe_0 + tr_z_yzzzz_xyyzz[i] * pa_y[i];

        tr_z_yyzzzz_xyzzz[i] = tr_z_zzzz_xyzzz[i] * fe_0 + tr_z_yzzzz_xzzz[i] * fe_0 + tr_z_yzzzz_xyzzz[i] * pa_y[i];

        tr_z_yyzzzz_xzzzz[i] = tr_z_zzzz_xzzzz[i] * fe_0 + tr_z_yzzzz_xzzzz[i] * pa_y[i];

        tr_z_yyzzzz_yyyyy[i] = tr_z_zzzz_yyyyy[i] * fe_0 + 5.0 * tr_z_yzzzz_yyyy[i] * fe_0 + tr_z_yzzzz_yyyyy[i] * pa_y[i];

        tr_z_yyzzzz_yyyyz[i] = tr_z_zzzz_yyyyz[i] * fe_0 + 4.0 * tr_z_yzzzz_yyyz[i] * fe_0 + tr_z_yzzzz_yyyyz[i] * pa_y[i];

        tr_z_yyzzzz_yyyzz[i] = tr_z_zzzz_yyyzz[i] * fe_0 + 3.0 * tr_z_yzzzz_yyzz[i] * fe_0 + tr_z_yzzzz_yyyzz[i] * pa_y[i];

        tr_z_yyzzzz_yyzzz[i] = tr_z_zzzz_yyzzz[i] * fe_0 + 2.0 * tr_z_yzzzz_yzzz[i] * fe_0 + tr_z_yzzzz_yyzzz[i] * pa_y[i];

        tr_z_yyzzzz_yzzzz[i] = tr_z_zzzz_yzzzz[i] * fe_0 + tr_z_yzzzz_zzzz[i] * fe_0 + tr_z_yzzzz_yzzzz[i] * pa_y[i];

        tr_z_yyzzzz_zzzzz[i] = tr_z_zzzz_zzzzz[i] * fe_0 + tr_z_yzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1722-1743 components of targeted buffer : IH

    auto tr_z_yzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1722);

    auto tr_z_yzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1723);

    auto tr_z_yzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1724);

    auto tr_z_yzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1725);

    auto tr_z_yzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1726);

    auto tr_z_yzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1727);

    auto tr_z_yzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1728);

    auto tr_z_yzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1729);

    auto tr_z_yzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1730);

    auto tr_z_yzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1731);

    auto tr_z_yzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1732);

    auto tr_z_yzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1733);

    auto tr_z_yzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1734);

    auto tr_z_yzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1735);

    auto tr_z_yzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1736);

    auto tr_z_yzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1737);

    auto tr_z_yzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1738);

    auto tr_z_yzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1739);

    auto tr_z_yzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1740);

    auto tr_z_yzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1741);

    auto tr_z_yzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1742);

    #pragma omp simd aligned(pa_y, tr_z_yzzzzz_xxxxx, tr_z_yzzzzz_xxxxy, tr_z_yzzzzz_xxxxz, tr_z_yzzzzz_xxxyy, tr_z_yzzzzz_xxxyz, tr_z_yzzzzz_xxxzz, tr_z_yzzzzz_xxyyy, tr_z_yzzzzz_xxyyz, tr_z_yzzzzz_xxyzz, tr_z_yzzzzz_xxzzz, tr_z_yzzzzz_xyyyy, tr_z_yzzzzz_xyyyz, tr_z_yzzzzz_xyyzz, tr_z_yzzzzz_xyzzz, tr_z_yzzzzz_xzzzz, tr_z_yzzzzz_yyyyy, tr_z_yzzzzz_yyyyz, tr_z_yzzzzz_yyyzz, tr_z_yzzzzz_yyzzz, tr_z_yzzzzz_yzzzz, tr_z_yzzzzz_zzzzz, tr_z_zzzzz_xxxx, tr_z_zzzzz_xxxxx, tr_z_zzzzz_xxxxy, tr_z_zzzzz_xxxxz, tr_z_zzzzz_xxxy, tr_z_zzzzz_xxxyy, tr_z_zzzzz_xxxyz, tr_z_zzzzz_xxxz, tr_z_zzzzz_xxxzz, tr_z_zzzzz_xxyy, tr_z_zzzzz_xxyyy, tr_z_zzzzz_xxyyz, tr_z_zzzzz_xxyz, tr_z_zzzzz_xxyzz, tr_z_zzzzz_xxzz, tr_z_zzzzz_xxzzz, tr_z_zzzzz_xyyy, tr_z_zzzzz_xyyyy, tr_z_zzzzz_xyyyz, tr_z_zzzzz_xyyz, tr_z_zzzzz_xyyzz, tr_z_zzzzz_xyzz, tr_z_zzzzz_xyzzz, tr_z_zzzzz_xzzz, tr_z_zzzzz_xzzzz, tr_z_zzzzz_yyyy, tr_z_zzzzz_yyyyy, tr_z_zzzzz_yyyyz, tr_z_zzzzz_yyyz, tr_z_zzzzz_yyyzz, tr_z_zzzzz_yyzz, tr_z_zzzzz_yyzzz, tr_z_zzzzz_yzzz, tr_z_zzzzz_yzzzz, tr_z_zzzzz_zzzz, tr_z_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzzz_xxxxx[i] = tr_z_zzzzz_xxxxx[i] * pa_y[i];

        tr_z_yzzzzz_xxxxy[i] = tr_z_zzzzz_xxxx[i] * fe_0 + tr_z_zzzzz_xxxxy[i] * pa_y[i];

        tr_z_yzzzzz_xxxxz[i] = tr_z_zzzzz_xxxxz[i] * pa_y[i];

        tr_z_yzzzzz_xxxyy[i] = 2.0 * tr_z_zzzzz_xxxy[i] * fe_0 + tr_z_zzzzz_xxxyy[i] * pa_y[i];

        tr_z_yzzzzz_xxxyz[i] = tr_z_zzzzz_xxxz[i] * fe_0 + tr_z_zzzzz_xxxyz[i] * pa_y[i];

        tr_z_yzzzzz_xxxzz[i] = tr_z_zzzzz_xxxzz[i] * pa_y[i];

        tr_z_yzzzzz_xxyyy[i] = 3.0 * tr_z_zzzzz_xxyy[i] * fe_0 + tr_z_zzzzz_xxyyy[i] * pa_y[i];

        tr_z_yzzzzz_xxyyz[i] = 2.0 * tr_z_zzzzz_xxyz[i] * fe_0 + tr_z_zzzzz_xxyyz[i] * pa_y[i];

        tr_z_yzzzzz_xxyzz[i] = tr_z_zzzzz_xxzz[i] * fe_0 + tr_z_zzzzz_xxyzz[i] * pa_y[i];

        tr_z_yzzzzz_xxzzz[i] = tr_z_zzzzz_xxzzz[i] * pa_y[i];

        tr_z_yzzzzz_xyyyy[i] = 4.0 * tr_z_zzzzz_xyyy[i] * fe_0 + tr_z_zzzzz_xyyyy[i] * pa_y[i];

        tr_z_yzzzzz_xyyyz[i] = 3.0 * tr_z_zzzzz_xyyz[i] * fe_0 + tr_z_zzzzz_xyyyz[i] * pa_y[i];

        tr_z_yzzzzz_xyyzz[i] = 2.0 * tr_z_zzzzz_xyzz[i] * fe_0 + tr_z_zzzzz_xyyzz[i] * pa_y[i];

        tr_z_yzzzzz_xyzzz[i] = tr_z_zzzzz_xzzz[i] * fe_0 + tr_z_zzzzz_xyzzz[i] * pa_y[i];

        tr_z_yzzzzz_xzzzz[i] = tr_z_zzzzz_xzzzz[i] * pa_y[i];

        tr_z_yzzzzz_yyyyy[i] = 5.0 * tr_z_zzzzz_yyyy[i] * fe_0 + tr_z_zzzzz_yyyyy[i] * pa_y[i];

        tr_z_yzzzzz_yyyyz[i] = 4.0 * tr_z_zzzzz_yyyz[i] * fe_0 + tr_z_zzzzz_yyyyz[i] * pa_y[i];

        tr_z_yzzzzz_yyyzz[i] = 3.0 * tr_z_zzzzz_yyzz[i] * fe_0 + tr_z_zzzzz_yyyzz[i] * pa_y[i];

        tr_z_yzzzzz_yyzzz[i] = 2.0 * tr_z_zzzzz_yzzz[i] * fe_0 + tr_z_zzzzz_yyzzz[i] * pa_y[i];

        tr_z_yzzzzz_yzzzz[i] = tr_z_zzzzz_zzzz[i] * fe_0 + tr_z_zzzzz_yzzzz[i] * pa_y[i];

        tr_z_yzzzzz_zzzzz[i] = tr_z_zzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1743-1764 components of targeted buffer : IH

    auto tr_z_zzzzzz_xxxxx = pbuffer.data(idx_dip_ih + 1743);

    auto tr_z_zzzzzz_xxxxy = pbuffer.data(idx_dip_ih + 1744);

    auto tr_z_zzzzzz_xxxxz = pbuffer.data(idx_dip_ih + 1745);

    auto tr_z_zzzzzz_xxxyy = pbuffer.data(idx_dip_ih + 1746);

    auto tr_z_zzzzzz_xxxyz = pbuffer.data(idx_dip_ih + 1747);

    auto tr_z_zzzzzz_xxxzz = pbuffer.data(idx_dip_ih + 1748);

    auto tr_z_zzzzzz_xxyyy = pbuffer.data(idx_dip_ih + 1749);

    auto tr_z_zzzzzz_xxyyz = pbuffer.data(idx_dip_ih + 1750);

    auto tr_z_zzzzzz_xxyzz = pbuffer.data(idx_dip_ih + 1751);

    auto tr_z_zzzzzz_xxzzz = pbuffer.data(idx_dip_ih + 1752);

    auto tr_z_zzzzzz_xyyyy = pbuffer.data(idx_dip_ih + 1753);

    auto tr_z_zzzzzz_xyyyz = pbuffer.data(idx_dip_ih + 1754);

    auto tr_z_zzzzzz_xyyzz = pbuffer.data(idx_dip_ih + 1755);

    auto tr_z_zzzzzz_xyzzz = pbuffer.data(idx_dip_ih + 1756);

    auto tr_z_zzzzzz_xzzzz = pbuffer.data(idx_dip_ih + 1757);

    auto tr_z_zzzzzz_yyyyy = pbuffer.data(idx_dip_ih + 1758);

    auto tr_z_zzzzzz_yyyyz = pbuffer.data(idx_dip_ih + 1759);

    auto tr_z_zzzzzz_yyyzz = pbuffer.data(idx_dip_ih + 1760);

    auto tr_z_zzzzzz_yyzzz = pbuffer.data(idx_dip_ih + 1761);

    auto tr_z_zzzzzz_yzzzz = pbuffer.data(idx_dip_ih + 1762);

    auto tr_z_zzzzzz_zzzzz = pbuffer.data(idx_dip_ih + 1763);

    #pragma omp simd aligned(pa_z, tr_z_zzzz_xxxxx, tr_z_zzzz_xxxxy, tr_z_zzzz_xxxxz, tr_z_zzzz_xxxyy, tr_z_zzzz_xxxyz, tr_z_zzzz_xxxzz, tr_z_zzzz_xxyyy, tr_z_zzzz_xxyyz, tr_z_zzzz_xxyzz, tr_z_zzzz_xxzzz, tr_z_zzzz_xyyyy, tr_z_zzzz_xyyyz, tr_z_zzzz_xyyzz, tr_z_zzzz_xyzzz, tr_z_zzzz_xzzzz, tr_z_zzzz_yyyyy, tr_z_zzzz_yyyyz, tr_z_zzzz_yyyzz, tr_z_zzzz_yyzzz, tr_z_zzzz_yzzzz, tr_z_zzzz_zzzzz, tr_z_zzzzz_xxxx, tr_z_zzzzz_xxxxx, tr_z_zzzzz_xxxxy, tr_z_zzzzz_xxxxz, tr_z_zzzzz_xxxy, tr_z_zzzzz_xxxyy, tr_z_zzzzz_xxxyz, tr_z_zzzzz_xxxz, tr_z_zzzzz_xxxzz, tr_z_zzzzz_xxyy, tr_z_zzzzz_xxyyy, tr_z_zzzzz_xxyyz, tr_z_zzzzz_xxyz, tr_z_zzzzz_xxyzz, tr_z_zzzzz_xxzz, tr_z_zzzzz_xxzzz, tr_z_zzzzz_xyyy, tr_z_zzzzz_xyyyy, tr_z_zzzzz_xyyyz, tr_z_zzzzz_xyyz, tr_z_zzzzz_xyyzz, tr_z_zzzzz_xyzz, tr_z_zzzzz_xyzzz, tr_z_zzzzz_xzzz, tr_z_zzzzz_xzzzz, tr_z_zzzzz_yyyy, tr_z_zzzzz_yyyyy, tr_z_zzzzz_yyyyz, tr_z_zzzzz_yyyz, tr_z_zzzzz_yyyzz, tr_z_zzzzz_yyzz, tr_z_zzzzz_yyzzz, tr_z_zzzzz_yzzz, tr_z_zzzzz_yzzzz, tr_z_zzzzz_zzzz, tr_z_zzzzz_zzzzz, tr_z_zzzzzz_xxxxx, tr_z_zzzzzz_xxxxy, tr_z_zzzzzz_xxxxz, tr_z_zzzzzz_xxxyy, tr_z_zzzzzz_xxxyz, tr_z_zzzzzz_xxxzz, tr_z_zzzzzz_xxyyy, tr_z_zzzzzz_xxyyz, tr_z_zzzzzz_xxyzz, tr_z_zzzzzz_xxzzz, tr_z_zzzzzz_xyyyy, tr_z_zzzzzz_xyyyz, tr_z_zzzzzz_xyyzz, tr_z_zzzzzz_xyzzz, tr_z_zzzzzz_xzzzz, tr_z_zzzzzz_yyyyy, tr_z_zzzzzz_yyyyz, tr_z_zzzzzz_yyyzz, tr_z_zzzzzz_yyzzz, tr_z_zzzzzz_yzzzz, tr_z_zzzzzz_zzzzz, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxzz, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzzz_xxxxx[i] = 5.0 * tr_z_zzzz_xxxxx[i] * fe_0 + ts_zzzzz_xxxxx[i] * fe_0 + tr_z_zzzzz_xxxxx[i] * pa_z[i];

        tr_z_zzzzzz_xxxxy[i] = 5.0 * tr_z_zzzz_xxxxy[i] * fe_0 + ts_zzzzz_xxxxy[i] * fe_0 + tr_z_zzzzz_xxxxy[i] * pa_z[i];

        tr_z_zzzzzz_xxxxz[i] = 5.0 * tr_z_zzzz_xxxxz[i] * fe_0 + tr_z_zzzzz_xxxx[i] * fe_0 + ts_zzzzz_xxxxz[i] * fe_0 + tr_z_zzzzz_xxxxz[i] * pa_z[i];

        tr_z_zzzzzz_xxxyy[i] = 5.0 * tr_z_zzzz_xxxyy[i] * fe_0 + ts_zzzzz_xxxyy[i] * fe_0 + tr_z_zzzzz_xxxyy[i] * pa_z[i];

        tr_z_zzzzzz_xxxyz[i] = 5.0 * tr_z_zzzz_xxxyz[i] * fe_0 + tr_z_zzzzz_xxxy[i] * fe_0 + ts_zzzzz_xxxyz[i] * fe_0 + tr_z_zzzzz_xxxyz[i] * pa_z[i];

        tr_z_zzzzzz_xxxzz[i] = 5.0 * tr_z_zzzz_xxxzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xxxz[i] * fe_0 + ts_zzzzz_xxxzz[i] * fe_0 + tr_z_zzzzz_xxxzz[i] * pa_z[i];

        tr_z_zzzzzz_xxyyy[i] = 5.0 * tr_z_zzzz_xxyyy[i] * fe_0 + ts_zzzzz_xxyyy[i] * fe_0 + tr_z_zzzzz_xxyyy[i] * pa_z[i];

        tr_z_zzzzzz_xxyyz[i] = 5.0 * tr_z_zzzz_xxyyz[i] * fe_0 + tr_z_zzzzz_xxyy[i] * fe_0 + ts_zzzzz_xxyyz[i] * fe_0 + tr_z_zzzzz_xxyyz[i] * pa_z[i];

        tr_z_zzzzzz_xxyzz[i] = 5.0 * tr_z_zzzz_xxyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xxyz[i] * fe_0 + ts_zzzzz_xxyzz[i] * fe_0 + tr_z_zzzzz_xxyzz[i] * pa_z[i];

        tr_z_zzzzzz_xxzzz[i] = 5.0 * tr_z_zzzz_xxzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_xxzz[i] * fe_0 + ts_zzzzz_xxzzz[i] * fe_0 + tr_z_zzzzz_xxzzz[i] * pa_z[i];

        tr_z_zzzzzz_xyyyy[i] = 5.0 * tr_z_zzzz_xyyyy[i] * fe_0 + ts_zzzzz_xyyyy[i] * fe_0 + tr_z_zzzzz_xyyyy[i] * pa_z[i];

        tr_z_zzzzzz_xyyyz[i] = 5.0 * tr_z_zzzz_xyyyz[i] * fe_0 + tr_z_zzzzz_xyyy[i] * fe_0 + ts_zzzzz_xyyyz[i] * fe_0 + tr_z_zzzzz_xyyyz[i] * pa_z[i];

        tr_z_zzzzzz_xyyzz[i] = 5.0 * tr_z_zzzz_xyyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xyyz[i] * fe_0 + ts_zzzzz_xyyzz[i] * fe_0 + tr_z_zzzzz_xyyzz[i] * pa_z[i];

        tr_z_zzzzzz_xyzzz[i] = 5.0 * tr_z_zzzz_xyzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_xyzz[i] * fe_0 + ts_zzzzz_xyzzz[i] * fe_0 + tr_z_zzzzz_xyzzz[i] * pa_z[i];

        tr_z_zzzzzz_xzzzz[i] = 5.0 * tr_z_zzzz_xzzzz[i] * fe_0 + 4.0 * tr_z_zzzzz_xzzz[i] * fe_0 + ts_zzzzz_xzzzz[i] * fe_0 + tr_z_zzzzz_xzzzz[i] * pa_z[i];

        tr_z_zzzzzz_yyyyy[i] = 5.0 * tr_z_zzzz_yyyyy[i] * fe_0 + ts_zzzzz_yyyyy[i] * fe_0 + tr_z_zzzzz_yyyyy[i] * pa_z[i];

        tr_z_zzzzzz_yyyyz[i] = 5.0 * tr_z_zzzz_yyyyz[i] * fe_0 + tr_z_zzzzz_yyyy[i] * fe_0 + ts_zzzzz_yyyyz[i] * fe_0 + tr_z_zzzzz_yyyyz[i] * pa_z[i];

        tr_z_zzzzzz_yyyzz[i] = 5.0 * tr_z_zzzz_yyyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_yyyz[i] * fe_0 + ts_zzzzz_yyyzz[i] * fe_0 + tr_z_zzzzz_yyyzz[i] * pa_z[i];

        tr_z_zzzzzz_yyzzz[i] = 5.0 * tr_z_zzzz_yyzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_yyzz[i] * fe_0 + ts_zzzzz_yyzzz[i] * fe_0 + tr_z_zzzzz_yyzzz[i] * pa_z[i];

        tr_z_zzzzzz_yzzzz[i] = 5.0 * tr_z_zzzz_yzzzz[i] * fe_0 + 4.0 * tr_z_zzzzz_yzzz[i] * fe_0 + ts_zzzzz_yzzzz[i] * fe_0 + tr_z_zzzzz_yzzzz[i] * pa_z[i];

        tr_z_zzzzzz_zzzzz[i] = 5.0 * tr_z_zzzz_zzzzz[i] * fe_0 + 5.0 * tr_z_zzzzz_zzzz[i] * fe_0 + ts_zzzzz_zzzzz[i] * fe_0 + tr_z_zzzzz_zzzzz[i] * pa_z[i];
    }

}

} // diprec namespace

