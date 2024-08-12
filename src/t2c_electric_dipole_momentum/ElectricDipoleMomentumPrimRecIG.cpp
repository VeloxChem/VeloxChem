#include "ElectricDipoleMomentumPrimRecIG.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_ig(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_ig,
                                      const size_t              idx_dip_gg,
                                      const size_t              idx_dip_hf,
                                      const size_t              idx_ovl_hg,
                                      const size_t              idx_dip_hg,
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

    // Set up components of auxiliary buffer : GG

    auto tr_x_xxxx_xxxx = pbuffer.data(idx_dip_gg);

    auto tr_x_xxxx_xxxy = pbuffer.data(idx_dip_gg + 1);

    auto tr_x_xxxx_xxxz = pbuffer.data(idx_dip_gg + 2);

    auto tr_x_xxxx_xxyy = pbuffer.data(idx_dip_gg + 3);

    auto tr_x_xxxx_xxyz = pbuffer.data(idx_dip_gg + 4);

    auto tr_x_xxxx_xxzz = pbuffer.data(idx_dip_gg + 5);

    auto tr_x_xxxx_xyyy = pbuffer.data(idx_dip_gg + 6);

    auto tr_x_xxxx_xyyz = pbuffer.data(idx_dip_gg + 7);

    auto tr_x_xxxx_xyzz = pbuffer.data(idx_dip_gg + 8);

    auto tr_x_xxxx_xzzz = pbuffer.data(idx_dip_gg + 9);

    auto tr_x_xxxx_yyyy = pbuffer.data(idx_dip_gg + 10);

    auto tr_x_xxxx_yyyz = pbuffer.data(idx_dip_gg + 11);

    auto tr_x_xxxx_yyzz = pbuffer.data(idx_dip_gg + 12);

    auto tr_x_xxxx_yzzz = pbuffer.data(idx_dip_gg + 13);

    auto tr_x_xxxx_zzzz = pbuffer.data(idx_dip_gg + 14);

    auto tr_x_xxxy_xxxx = pbuffer.data(idx_dip_gg + 15);

    auto tr_x_xxxy_xxxy = pbuffer.data(idx_dip_gg + 16);

    auto tr_x_xxxy_xxxz = pbuffer.data(idx_dip_gg + 17);

    auto tr_x_xxxy_xxyy = pbuffer.data(idx_dip_gg + 18);

    auto tr_x_xxxy_xxyz = pbuffer.data(idx_dip_gg + 19);

    auto tr_x_xxxy_xxzz = pbuffer.data(idx_dip_gg + 20);

    auto tr_x_xxxy_xyyy = pbuffer.data(idx_dip_gg + 21);

    auto tr_x_xxxy_xyyz = pbuffer.data(idx_dip_gg + 22);

    auto tr_x_xxxy_xyzz = pbuffer.data(idx_dip_gg + 23);

    auto tr_x_xxxy_xzzz = pbuffer.data(idx_dip_gg + 24);

    auto tr_x_xxxy_zzzz = pbuffer.data(idx_dip_gg + 29);

    auto tr_x_xxxz_xxxx = pbuffer.data(idx_dip_gg + 30);

    auto tr_x_xxxz_xxxy = pbuffer.data(idx_dip_gg + 31);

    auto tr_x_xxxz_xxxz = pbuffer.data(idx_dip_gg + 32);

    auto tr_x_xxxz_xxyy = pbuffer.data(idx_dip_gg + 33);

    auto tr_x_xxxz_xxyz = pbuffer.data(idx_dip_gg + 34);

    auto tr_x_xxxz_xxzz = pbuffer.data(idx_dip_gg + 35);

    auto tr_x_xxxz_xyyy = pbuffer.data(idx_dip_gg + 36);

    auto tr_x_xxxz_xyyz = pbuffer.data(idx_dip_gg + 37);

    auto tr_x_xxxz_xyzz = pbuffer.data(idx_dip_gg + 38);

    auto tr_x_xxxz_xzzz = pbuffer.data(idx_dip_gg + 39);

    auto tr_x_xxxz_yyyy = pbuffer.data(idx_dip_gg + 40);

    auto tr_x_xxxz_zzzz = pbuffer.data(idx_dip_gg + 44);

    auto tr_x_xxyy_xxxx = pbuffer.data(idx_dip_gg + 45);

    auto tr_x_xxyy_xxxy = pbuffer.data(idx_dip_gg + 46);

    auto tr_x_xxyy_xxxz = pbuffer.data(idx_dip_gg + 47);

    auto tr_x_xxyy_xxyy = pbuffer.data(idx_dip_gg + 48);

    auto tr_x_xxyy_xxyz = pbuffer.data(idx_dip_gg + 49);

    auto tr_x_xxyy_xxzz = pbuffer.data(idx_dip_gg + 50);

    auto tr_x_xxyy_xyyy = pbuffer.data(idx_dip_gg + 51);

    auto tr_x_xxyy_xyyz = pbuffer.data(idx_dip_gg + 52);

    auto tr_x_xxyy_xyzz = pbuffer.data(idx_dip_gg + 53);

    auto tr_x_xxyy_xzzz = pbuffer.data(idx_dip_gg + 54);

    auto tr_x_xxyy_yyyy = pbuffer.data(idx_dip_gg + 55);

    auto tr_x_xxyy_yyyz = pbuffer.data(idx_dip_gg + 56);

    auto tr_x_xxyy_yyzz = pbuffer.data(idx_dip_gg + 57);

    auto tr_x_xxyy_yzzz = pbuffer.data(idx_dip_gg + 58);

    auto tr_x_xxyy_zzzz = pbuffer.data(idx_dip_gg + 59);

    auto tr_x_xxyz_xxxz = pbuffer.data(idx_dip_gg + 62);

    auto tr_x_xxyz_xxzz = pbuffer.data(idx_dip_gg + 65);

    auto tr_x_xxyz_xzzz = pbuffer.data(idx_dip_gg + 69);

    auto tr_x_xxyz_zzzz = pbuffer.data(idx_dip_gg + 74);

    auto tr_x_xxzz_xxxx = pbuffer.data(idx_dip_gg + 75);

    auto tr_x_xxzz_xxxy = pbuffer.data(idx_dip_gg + 76);

    auto tr_x_xxzz_xxxz = pbuffer.data(idx_dip_gg + 77);

    auto tr_x_xxzz_xxyy = pbuffer.data(idx_dip_gg + 78);

    auto tr_x_xxzz_xxyz = pbuffer.data(idx_dip_gg + 79);

    auto tr_x_xxzz_xxzz = pbuffer.data(idx_dip_gg + 80);

    auto tr_x_xxzz_xyyy = pbuffer.data(idx_dip_gg + 81);

    auto tr_x_xxzz_xyyz = pbuffer.data(idx_dip_gg + 82);

    auto tr_x_xxzz_xyzz = pbuffer.data(idx_dip_gg + 83);

    auto tr_x_xxzz_xzzz = pbuffer.data(idx_dip_gg + 84);

    auto tr_x_xxzz_yyyy = pbuffer.data(idx_dip_gg + 85);

    auto tr_x_xxzz_yyyz = pbuffer.data(idx_dip_gg + 86);

    auto tr_x_xxzz_yyzz = pbuffer.data(idx_dip_gg + 87);

    auto tr_x_xxzz_yzzz = pbuffer.data(idx_dip_gg + 88);

    auto tr_x_xxzz_zzzz = pbuffer.data(idx_dip_gg + 89);

    auto tr_x_xyyy_xxxx = pbuffer.data(idx_dip_gg + 90);

    auto tr_x_xyyy_xxxy = pbuffer.data(idx_dip_gg + 91);

    auto tr_x_xyyy_xxxz = pbuffer.data(idx_dip_gg + 92);

    auto tr_x_xyyy_xxyy = pbuffer.data(idx_dip_gg + 93);

    auto tr_x_xyyy_xxzz = pbuffer.data(idx_dip_gg + 95);

    auto tr_x_xyyy_xyyy = pbuffer.data(idx_dip_gg + 96);

    auto tr_x_xyyy_xzzz = pbuffer.data(idx_dip_gg + 99);

    auto tr_x_xyyy_yyyy = pbuffer.data(idx_dip_gg + 100);

    auto tr_x_xyyy_yyyz = pbuffer.data(idx_dip_gg + 101);

    auto tr_x_xyyy_yyzz = pbuffer.data(idx_dip_gg + 102);

    auto tr_x_xyyy_yzzz = pbuffer.data(idx_dip_gg + 103);

    auto tr_x_xyyz_xxxy = pbuffer.data(idx_dip_gg + 106);

    auto tr_x_xyyz_xxxz = pbuffer.data(idx_dip_gg + 107);

    auto tr_x_xyyz_xxyy = pbuffer.data(idx_dip_gg + 108);

    auto tr_x_xyyz_xxzz = pbuffer.data(idx_dip_gg + 110);

    auto tr_x_xyyz_xyyy = pbuffer.data(idx_dip_gg + 111);

    auto tr_x_xyyz_xzzz = pbuffer.data(idx_dip_gg + 114);

    auto tr_x_xyzz_xxxx = pbuffer.data(idx_dip_gg + 120);

    auto tr_x_xyzz_xxxz = pbuffer.data(idx_dip_gg + 122);

    auto tr_x_xyzz_xxzz = pbuffer.data(idx_dip_gg + 125);

    auto tr_x_xyzz_xzzz = pbuffer.data(idx_dip_gg + 129);

    auto tr_x_xzzz_xxxx = pbuffer.data(idx_dip_gg + 135);

    auto tr_x_xzzz_xxxy = pbuffer.data(idx_dip_gg + 136);

    auto tr_x_xzzz_xxxz = pbuffer.data(idx_dip_gg + 137);

    auto tr_x_xzzz_xxyy = pbuffer.data(idx_dip_gg + 138);

    auto tr_x_xzzz_xxzz = pbuffer.data(idx_dip_gg + 140);

    auto tr_x_xzzz_xyyy = pbuffer.data(idx_dip_gg + 141);

    auto tr_x_xzzz_xzzz = pbuffer.data(idx_dip_gg + 144);

    auto tr_x_xzzz_yyyz = pbuffer.data(idx_dip_gg + 146);

    auto tr_x_xzzz_yyzz = pbuffer.data(idx_dip_gg + 147);

    auto tr_x_xzzz_yzzz = pbuffer.data(idx_dip_gg + 148);

    auto tr_x_xzzz_zzzz = pbuffer.data(idx_dip_gg + 149);

    auto tr_x_yyyy_xxxx = pbuffer.data(idx_dip_gg + 150);

    auto tr_x_yyyy_xxxy = pbuffer.data(idx_dip_gg + 151);

    auto tr_x_yyyy_xxxz = pbuffer.data(idx_dip_gg + 152);

    auto tr_x_yyyy_xxyy = pbuffer.data(idx_dip_gg + 153);

    auto tr_x_yyyy_xxyz = pbuffer.data(idx_dip_gg + 154);

    auto tr_x_yyyy_xxzz = pbuffer.data(idx_dip_gg + 155);

    auto tr_x_yyyy_xyyy = pbuffer.data(idx_dip_gg + 156);

    auto tr_x_yyyy_xyyz = pbuffer.data(idx_dip_gg + 157);

    auto tr_x_yyyy_xyzz = pbuffer.data(idx_dip_gg + 158);

    auto tr_x_yyyy_xzzz = pbuffer.data(idx_dip_gg + 159);

    auto tr_x_yyyy_yyyy = pbuffer.data(idx_dip_gg + 160);

    auto tr_x_yyyy_yyyz = pbuffer.data(idx_dip_gg + 161);

    auto tr_x_yyyy_yyzz = pbuffer.data(idx_dip_gg + 162);

    auto tr_x_yyyy_yzzz = pbuffer.data(idx_dip_gg + 163);

    auto tr_x_yyyy_zzzz = pbuffer.data(idx_dip_gg + 164);

    auto tr_x_yyyz_xxxy = pbuffer.data(idx_dip_gg + 166);

    auto tr_x_yyyz_xxxz = pbuffer.data(idx_dip_gg + 167);

    auto tr_x_yyyz_xxyy = pbuffer.data(idx_dip_gg + 168);

    auto tr_x_yyyz_xxzz = pbuffer.data(idx_dip_gg + 170);

    auto tr_x_yyyz_xyyy = pbuffer.data(idx_dip_gg + 171);

    auto tr_x_yyyz_xzzz = pbuffer.data(idx_dip_gg + 174);

    auto tr_x_yyyz_yyyy = pbuffer.data(idx_dip_gg + 175);

    auto tr_x_yyyz_zzzz = pbuffer.data(idx_dip_gg + 179);

    auto tr_x_yyzz_xxxx = pbuffer.data(idx_dip_gg + 180);

    auto tr_x_yyzz_xxxy = pbuffer.data(idx_dip_gg + 181);

    auto tr_x_yyzz_xxxz = pbuffer.data(idx_dip_gg + 182);

    auto tr_x_yyzz_xxyy = pbuffer.data(idx_dip_gg + 183);

    auto tr_x_yyzz_xxyz = pbuffer.data(idx_dip_gg + 184);

    auto tr_x_yyzz_xxzz = pbuffer.data(idx_dip_gg + 185);

    auto tr_x_yyzz_xyyy = pbuffer.data(idx_dip_gg + 186);

    auto tr_x_yyzz_xyyz = pbuffer.data(idx_dip_gg + 187);

    auto tr_x_yyzz_xyzz = pbuffer.data(idx_dip_gg + 188);

    auto tr_x_yyzz_xzzz = pbuffer.data(idx_dip_gg + 189);

    auto tr_x_yyzz_yyyy = pbuffer.data(idx_dip_gg + 190);

    auto tr_x_yyzz_yyyz = pbuffer.data(idx_dip_gg + 191);

    auto tr_x_yyzz_yyzz = pbuffer.data(idx_dip_gg + 192);

    auto tr_x_yyzz_yzzz = pbuffer.data(idx_dip_gg + 193);

    auto tr_x_yyzz_zzzz = pbuffer.data(idx_dip_gg + 194);

    auto tr_x_yzzz_xxxx = pbuffer.data(idx_dip_gg + 195);

    auto tr_x_yzzz_xxxz = pbuffer.data(idx_dip_gg + 197);

    auto tr_x_yzzz_xxyz = pbuffer.data(idx_dip_gg + 199);

    auto tr_x_yzzz_xxzz = pbuffer.data(idx_dip_gg + 200);

    auto tr_x_yzzz_xyyz = pbuffer.data(idx_dip_gg + 202);

    auto tr_x_yzzz_xyzz = pbuffer.data(idx_dip_gg + 203);

    auto tr_x_yzzz_xzzz = pbuffer.data(idx_dip_gg + 204);

    auto tr_x_yzzz_yyyz = pbuffer.data(idx_dip_gg + 206);

    auto tr_x_yzzz_yyzz = pbuffer.data(idx_dip_gg + 207);

    auto tr_x_yzzz_yzzz = pbuffer.data(idx_dip_gg + 208);

    auto tr_x_yzzz_zzzz = pbuffer.data(idx_dip_gg + 209);

    auto tr_x_zzzz_xxxx = pbuffer.data(idx_dip_gg + 210);

    auto tr_x_zzzz_xxxy = pbuffer.data(idx_dip_gg + 211);

    auto tr_x_zzzz_xxxz = pbuffer.data(idx_dip_gg + 212);

    auto tr_x_zzzz_xxyy = pbuffer.data(idx_dip_gg + 213);

    auto tr_x_zzzz_xxyz = pbuffer.data(idx_dip_gg + 214);

    auto tr_x_zzzz_xxzz = pbuffer.data(idx_dip_gg + 215);

    auto tr_x_zzzz_xyyy = pbuffer.data(idx_dip_gg + 216);

    auto tr_x_zzzz_xyyz = pbuffer.data(idx_dip_gg + 217);

    auto tr_x_zzzz_xyzz = pbuffer.data(idx_dip_gg + 218);

    auto tr_x_zzzz_xzzz = pbuffer.data(idx_dip_gg + 219);

    auto tr_x_zzzz_yyyy = pbuffer.data(idx_dip_gg + 220);

    auto tr_x_zzzz_yyyz = pbuffer.data(idx_dip_gg + 221);

    auto tr_x_zzzz_yyzz = pbuffer.data(idx_dip_gg + 222);

    auto tr_x_zzzz_yzzz = pbuffer.data(idx_dip_gg + 223);

    auto tr_x_zzzz_zzzz = pbuffer.data(idx_dip_gg + 224);

    auto tr_y_xxxx_xxxx = pbuffer.data(idx_dip_gg + 225);

    auto tr_y_xxxx_xxxy = pbuffer.data(idx_dip_gg + 226);

    auto tr_y_xxxx_xxxz = pbuffer.data(idx_dip_gg + 227);

    auto tr_y_xxxx_xxyy = pbuffer.data(idx_dip_gg + 228);

    auto tr_y_xxxx_xxyz = pbuffer.data(idx_dip_gg + 229);

    auto tr_y_xxxx_xxzz = pbuffer.data(idx_dip_gg + 230);

    auto tr_y_xxxx_xyyy = pbuffer.data(idx_dip_gg + 231);

    auto tr_y_xxxx_xyyz = pbuffer.data(idx_dip_gg + 232);

    auto tr_y_xxxx_xyzz = pbuffer.data(idx_dip_gg + 233);

    auto tr_y_xxxx_xzzz = pbuffer.data(idx_dip_gg + 234);

    auto tr_y_xxxx_yyyy = pbuffer.data(idx_dip_gg + 235);

    auto tr_y_xxxx_yyyz = pbuffer.data(idx_dip_gg + 236);

    auto tr_y_xxxx_yyzz = pbuffer.data(idx_dip_gg + 237);

    auto tr_y_xxxx_yzzz = pbuffer.data(idx_dip_gg + 238);

    auto tr_y_xxxx_zzzz = pbuffer.data(idx_dip_gg + 239);

    auto tr_y_xxxy_xxxy = pbuffer.data(idx_dip_gg + 241);

    auto tr_y_xxxy_xxyy = pbuffer.data(idx_dip_gg + 243);

    auto tr_y_xxxy_xxyz = pbuffer.data(idx_dip_gg + 244);

    auto tr_y_xxxy_xyyy = pbuffer.data(idx_dip_gg + 246);

    auto tr_y_xxxy_xyyz = pbuffer.data(idx_dip_gg + 247);

    auto tr_y_xxxy_xyzz = pbuffer.data(idx_dip_gg + 248);

    auto tr_y_xxxy_yyyy = pbuffer.data(idx_dip_gg + 250);

    auto tr_y_xxxy_yyyz = pbuffer.data(idx_dip_gg + 251);

    auto tr_y_xxxy_yyzz = pbuffer.data(idx_dip_gg + 252);

    auto tr_y_xxxy_yzzz = pbuffer.data(idx_dip_gg + 253);

    auto tr_y_xxxy_zzzz = pbuffer.data(idx_dip_gg + 254);

    auto tr_y_xxxz_xxxx = pbuffer.data(idx_dip_gg + 255);

    auto tr_y_xxxz_xxxy = pbuffer.data(idx_dip_gg + 256);

    auto tr_y_xxxz_xxyy = pbuffer.data(idx_dip_gg + 258);

    auto tr_y_xxxz_xyyy = pbuffer.data(idx_dip_gg + 261);

    auto tr_y_xxxz_yyyz = pbuffer.data(idx_dip_gg + 266);

    auto tr_y_xxxz_yyzz = pbuffer.data(idx_dip_gg + 267);

    auto tr_y_xxxz_yzzz = pbuffer.data(idx_dip_gg + 268);

    auto tr_y_xxxz_zzzz = pbuffer.data(idx_dip_gg + 269);

    auto tr_y_xxyy_xxxx = pbuffer.data(idx_dip_gg + 270);

    auto tr_y_xxyy_xxxy = pbuffer.data(idx_dip_gg + 271);

    auto tr_y_xxyy_xxxz = pbuffer.data(idx_dip_gg + 272);

    auto tr_y_xxyy_xxyy = pbuffer.data(idx_dip_gg + 273);

    auto tr_y_xxyy_xxyz = pbuffer.data(idx_dip_gg + 274);

    auto tr_y_xxyy_xxzz = pbuffer.data(idx_dip_gg + 275);

    auto tr_y_xxyy_xyyy = pbuffer.data(idx_dip_gg + 276);

    auto tr_y_xxyy_xyyz = pbuffer.data(idx_dip_gg + 277);

    auto tr_y_xxyy_xyzz = pbuffer.data(idx_dip_gg + 278);

    auto tr_y_xxyy_xzzz = pbuffer.data(idx_dip_gg + 279);

    auto tr_y_xxyy_yyyy = pbuffer.data(idx_dip_gg + 280);

    auto tr_y_xxyy_yyyz = pbuffer.data(idx_dip_gg + 281);

    auto tr_y_xxyy_yyzz = pbuffer.data(idx_dip_gg + 282);

    auto tr_y_xxyy_yzzz = pbuffer.data(idx_dip_gg + 283);

    auto tr_y_xxyy_zzzz = pbuffer.data(idx_dip_gg + 284);

    auto tr_y_xxyz_xxxy = pbuffer.data(idx_dip_gg + 286);

    auto tr_y_xxyz_xxyy = pbuffer.data(idx_dip_gg + 288);

    auto tr_y_xxyz_xyyy = pbuffer.data(idx_dip_gg + 291);

    auto tr_y_xxyz_yyyz = pbuffer.data(idx_dip_gg + 296);

    auto tr_y_xxyz_yyzz = pbuffer.data(idx_dip_gg + 297);

    auto tr_y_xxyz_yzzz = pbuffer.data(idx_dip_gg + 298);

    auto tr_y_xxyz_zzzz = pbuffer.data(idx_dip_gg + 299);

    auto tr_y_xxzz_xxxx = pbuffer.data(idx_dip_gg + 300);

    auto tr_y_xxzz_xxxy = pbuffer.data(idx_dip_gg + 301);

    auto tr_y_xxzz_xxxz = pbuffer.data(idx_dip_gg + 302);

    auto tr_y_xxzz_xxyy = pbuffer.data(idx_dip_gg + 303);

    auto tr_y_xxzz_xxyz = pbuffer.data(idx_dip_gg + 304);

    auto tr_y_xxzz_xxzz = pbuffer.data(idx_dip_gg + 305);

    auto tr_y_xxzz_xyyy = pbuffer.data(idx_dip_gg + 306);

    auto tr_y_xxzz_xyyz = pbuffer.data(idx_dip_gg + 307);

    auto tr_y_xxzz_xyzz = pbuffer.data(idx_dip_gg + 308);

    auto tr_y_xxzz_xzzz = pbuffer.data(idx_dip_gg + 309);

    auto tr_y_xxzz_yyyy = pbuffer.data(idx_dip_gg + 310);

    auto tr_y_xxzz_yyyz = pbuffer.data(idx_dip_gg + 311);

    auto tr_y_xxzz_yyzz = pbuffer.data(idx_dip_gg + 312);

    auto tr_y_xxzz_yzzz = pbuffer.data(idx_dip_gg + 313);

    auto tr_y_xxzz_zzzz = pbuffer.data(idx_dip_gg + 314);

    auto tr_y_xyyy_xxxx = pbuffer.data(idx_dip_gg + 315);

    auto tr_y_xyyy_xxxy = pbuffer.data(idx_dip_gg + 316);

    auto tr_y_xyyy_xxxz = pbuffer.data(idx_dip_gg + 317);

    auto tr_y_xyyy_xxyy = pbuffer.data(idx_dip_gg + 318);

    auto tr_y_xyyy_xxyz = pbuffer.data(idx_dip_gg + 319);

    auto tr_y_xyyy_xxzz = pbuffer.data(idx_dip_gg + 320);

    auto tr_y_xyyy_xyyy = pbuffer.data(idx_dip_gg + 321);

    auto tr_y_xyyy_xyyz = pbuffer.data(idx_dip_gg + 322);

    auto tr_y_xyyy_xyzz = pbuffer.data(idx_dip_gg + 323);

    auto tr_y_xyyy_xzzz = pbuffer.data(idx_dip_gg + 324);

    auto tr_y_xyyy_yyyy = pbuffer.data(idx_dip_gg + 325);

    auto tr_y_xyyy_yyyz = pbuffer.data(idx_dip_gg + 326);

    auto tr_y_xyyy_yyzz = pbuffer.data(idx_dip_gg + 327);

    auto tr_y_xyyy_yzzz = pbuffer.data(idx_dip_gg + 328);

    auto tr_y_xyyy_zzzz = pbuffer.data(idx_dip_gg + 329);

    auto tr_y_xyyz_yyyz = pbuffer.data(idx_dip_gg + 341);

    auto tr_y_xyyz_yyzz = pbuffer.data(idx_dip_gg + 342);

    auto tr_y_xyyz_yzzz = pbuffer.data(idx_dip_gg + 343);

    auto tr_y_xyyz_zzzz = pbuffer.data(idx_dip_gg + 344);

    auto tr_y_xyzz_xxyz = pbuffer.data(idx_dip_gg + 349);

    auto tr_y_xyzz_xyyz = pbuffer.data(idx_dip_gg + 352);

    auto tr_y_xyzz_xyzz = pbuffer.data(idx_dip_gg + 353);

    auto tr_y_xyzz_yyyy = pbuffer.data(idx_dip_gg + 355);

    auto tr_y_xyzz_yyyz = pbuffer.data(idx_dip_gg + 356);

    auto tr_y_xyzz_yyzz = pbuffer.data(idx_dip_gg + 357);

    auto tr_y_xyzz_yzzz = pbuffer.data(idx_dip_gg + 358);

    auto tr_y_xyzz_zzzz = pbuffer.data(idx_dip_gg + 359);

    auto tr_y_xzzz_xxxz = pbuffer.data(idx_dip_gg + 362);

    auto tr_y_xzzz_xxyz = pbuffer.data(idx_dip_gg + 364);

    auto tr_y_xzzz_xxzz = pbuffer.data(idx_dip_gg + 365);

    auto tr_y_xzzz_xyyz = pbuffer.data(idx_dip_gg + 367);

    auto tr_y_xzzz_xyzz = pbuffer.data(idx_dip_gg + 368);

    auto tr_y_xzzz_xzzz = pbuffer.data(idx_dip_gg + 369);

    auto tr_y_xzzz_yyyy = pbuffer.data(idx_dip_gg + 370);

    auto tr_y_xzzz_yyyz = pbuffer.data(idx_dip_gg + 371);

    auto tr_y_xzzz_yyzz = pbuffer.data(idx_dip_gg + 372);

    auto tr_y_xzzz_yzzz = pbuffer.data(idx_dip_gg + 373);

    auto tr_y_xzzz_zzzz = pbuffer.data(idx_dip_gg + 374);

    auto tr_y_yyyy_xxxx = pbuffer.data(idx_dip_gg + 375);

    auto tr_y_yyyy_xxxy = pbuffer.data(idx_dip_gg + 376);

    auto tr_y_yyyy_xxxz = pbuffer.data(idx_dip_gg + 377);

    auto tr_y_yyyy_xxyy = pbuffer.data(idx_dip_gg + 378);

    auto tr_y_yyyy_xxyz = pbuffer.data(idx_dip_gg + 379);

    auto tr_y_yyyy_xxzz = pbuffer.data(idx_dip_gg + 380);

    auto tr_y_yyyy_xyyy = pbuffer.data(idx_dip_gg + 381);

    auto tr_y_yyyy_xyyz = pbuffer.data(idx_dip_gg + 382);

    auto tr_y_yyyy_xyzz = pbuffer.data(idx_dip_gg + 383);

    auto tr_y_yyyy_xzzz = pbuffer.data(idx_dip_gg + 384);

    auto tr_y_yyyy_yyyy = pbuffer.data(idx_dip_gg + 385);

    auto tr_y_yyyy_yyyz = pbuffer.data(idx_dip_gg + 386);

    auto tr_y_yyyy_yyzz = pbuffer.data(idx_dip_gg + 387);

    auto tr_y_yyyy_yzzz = pbuffer.data(idx_dip_gg + 388);

    auto tr_y_yyyy_zzzz = pbuffer.data(idx_dip_gg + 389);

    auto tr_y_yyyz_xxxx = pbuffer.data(idx_dip_gg + 390);

    auto tr_y_yyyz_xxxy = pbuffer.data(idx_dip_gg + 391);

    auto tr_y_yyyz_xxyy = pbuffer.data(idx_dip_gg + 393);

    auto tr_y_yyyz_xxyz = pbuffer.data(idx_dip_gg + 394);

    auto tr_y_yyyz_xyyy = pbuffer.data(idx_dip_gg + 396);

    auto tr_y_yyyz_xyyz = pbuffer.data(idx_dip_gg + 397);

    auto tr_y_yyyz_xyzz = pbuffer.data(idx_dip_gg + 398);

    auto tr_y_yyyz_yyyy = pbuffer.data(idx_dip_gg + 400);

    auto tr_y_yyyz_yyyz = pbuffer.data(idx_dip_gg + 401);

    auto tr_y_yyyz_yyzz = pbuffer.data(idx_dip_gg + 402);

    auto tr_y_yyyz_yzzz = pbuffer.data(idx_dip_gg + 403);

    auto tr_y_yyyz_zzzz = pbuffer.data(idx_dip_gg + 404);

    auto tr_y_yyzz_xxxx = pbuffer.data(idx_dip_gg + 405);

    auto tr_y_yyzz_xxxy = pbuffer.data(idx_dip_gg + 406);

    auto tr_y_yyzz_xxxz = pbuffer.data(idx_dip_gg + 407);

    auto tr_y_yyzz_xxyy = pbuffer.data(idx_dip_gg + 408);

    auto tr_y_yyzz_xxyz = pbuffer.data(idx_dip_gg + 409);

    auto tr_y_yyzz_xxzz = pbuffer.data(idx_dip_gg + 410);

    auto tr_y_yyzz_xyyy = pbuffer.data(idx_dip_gg + 411);

    auto tr_y_yyzz_xyyz = pbuffer.data(idx_dip_gg + 412);

    auto tr_y_yyzz_xyzz = pbuffer.data(idx_dip_gg + 413);

    auto tr_y_yyzz_xzzz = pbuffer.data(idx_dip_gg + 414);

    auto tr_y_yyzz_yyyy = pbuffer.data(idx_dip_gg + 415);

    auto tr_y_yyzz_yyyz = pbuffer.data(idx_dip_gg + 416);

    auto tr_y_yyzz_yyzz = pbuffer.data(idx_dip_gg + 417);

    auto tr_y_yyzz_yzzz = pbuffer.data(idx_dip_gg + 418);

    auto tr_y_yyzz_zzzz = pbuffer.data(idx_dip_gg + 419);

    auto tr_y_yzzz_xxxy = pbuffer.data(idx_dip_gg + 421);

    auto tr_y_yzzz_xxxz = pbuffer.data(idx_dip_gg + 422);

    auto tr_y_yzzz_xxyy = pbuffer.data(idx_dip_gg + 423);

    auto tr_y_yzzz_xxyz = pbuffer.data(idx_dip_gg + 424);

    auto tr_y_yzzz_xxzz = pbuffer.data(idx_dip_gg + 425);

    auto tr_y_yzzz_xyyy = pbuffer.data(idx_dip_gg + 426);

    auto tr_y_yzzz_xyyz = pbuffer.data(idx_dip_gg + 427);

    auto tr_y_yzzz_xyzz = pbuffer.data(idx_dip_gg + 428);

    auto tr_y_yzzz_xzzz = pbuffer.data(idx_dip_gg + 429);

    auto tr_y_yzzz_yyyy = pbuffer.data(idx_dip_gg + 430);

    auto tr_y_yzzz_yyyz = pbuffer.data(idx_dip_gg + 431);

    auto tr_y_yzzz_yyzz = pbuffer.data(idx_dip_gg + 432);

    auto tr_y_yzzz_yzzz = pbuffer.data(idx_dip_gg + 433);

    auto tr_y_yzzz_zzzz = pbuffer.data(idx_dip_gg + 434);

    auto tr_y_zzzz_xxxx = pbuffer.data(idx_dip_gg + 435);

    auto tr_y_zzzz_xxxy = pbuffer.data(idx_dip_gg + 436);

    auto tr_y_zzzz_xxxz = pbuffer.data(idx_dip_gg + 437);

    auto tr_y_zzzz_xxyy = pbuffer.data(idx_dip_gg + 438);

    auto tr_y_zzzz_xxyz = pbuffer.data(idx_dip_gg + 439);

    auto tr_y_zzzz_xxzz = pbuffer.data(idx_dip_gg + 440);

    auto tr_y_zzzz_xyyy = pbuffer.data(idx_dip_gg + 441);

    auto tr_y_zzzz_xyyz = pbuffer.data(idx_dip_gg + 442);

    auto tr_y_zzzz_xyzz = pbuffer.data(idx_dip_gg + 443);

    auto tr_y_zzzz_xzzz = pbuffer.data(idx_dip_gg + 444);

    auto tr_y_zzzz_yyyy = pbuffer.data(idx_dip_gg + 445);

    auto tr_y_zzzz_yyyz = pbuffer.data(idx_dip_gg + 446);

    auto tr_y_zzzz_yyzz = pbuffer.data(idx_dip_gg + 447);

    auto tr_y_zzzz_yzzz = pbuffer.data(idx_dip_gg + 448);

    auto tr_y_zzzz_zzzz = pbuffer.data(idx_dip_gg + 449);

    auto tr_z_xxxx_xxxx = pbuffer.data(idx_dip_gg + 450);

    auto tr_z_xxxx_xxxy = pbuffer.data(idx_dip_gg + 451);

    auto tr_z_xxxx_xxxz = pbuffer.data(idx_dip_gg + 452);

    auto tr_z_xxxx_xxyy = pbuffer.data(idx_dip_gg + 453);

    auto tr_z_xxxx_xxyz = pbuffer.data(idx_dip_gg + 454);

    auto tr_z_xxxx_xxzz = pbuffer.data(idx_dip_gg + 455);

    auto tr_z_xxxx_xyyy = pbuffer.data(idx_dip_gg + 456);

    auto tr_z_xxxx_xyyz = pbuffer.data(idx_dip_gg + 457);

    auto tr_z_xxxx_xyzz = pbuffer.data(idx_dip_gg + 458);

    auto tr_z_xxxx_xzzz = pbuffer.data(idx_dip_gg + 459);

    auto tr_z_xxxx_yyyy = pbuffer.data(idx_dip_gg + 460);

    auto tr_z_xxxx_yyyz = pbuffer.data(idx_dip_gg + 461);

    auto tr_z_xxxx_yyzz = pbuffer.data(idx_dip_gg + 462);

    auto tr_z_xxxx_yzzz = pbuffer.data(idx_dip_gg + 463);

    auto tr_z_xxxx_zzzz = pbuffer.data(idx_dip_gg + 464);

    auto tr_z_xxxy_xxxx = pbuffer.data(idx_dip_gg + 465);

    auto tr_z_xxxy_xxxz = pbuffer.data(idx_dip_gg + 467);

    auto tr_z_xxxy_xxzz = pbuffer.data(idx_dip_gg + 470);

    auto tr_z_xxxy_xzzz = pbuffer.data(idx_dip_gg + 474);

    auto tr_z_xxxy_yyyy = pbuffer.data(idx_dip_gg + 475);

    auto tr_z_xxxy_yyyz = pbuffer.data(idx_dip_gg + 476);

    auto tr_z_xxxy_yyzz = pbuffer.data(idx_dip_gg + 477);

    auto tr_z_xxxy_yzzz = pbuffer.data(idx_dip_gg + 478);

    auto tr_z_xxxz_xxxx = pbuffer.data(idx_dip_gg + 480);

    auto tr_z_xxxz_xxxz = pbuffer.data(idx_dip_gg + 482);

    auto tr_z_xxxz_xxyz = pbuffer.data(idx_dip_gg + 484);

    auto tr_z_xxxz_xxzz = pbuffer.data(idx_dip_gg + 485);

    auto tr_z_xxxz_xyyz = pbuffer.data(idx_dip_gg + 487);

    auto tr_z_xxxz_xyzz = pbuffer.data(idx_dip_gg + 488);

    auto tr_z_xxxz_xzzz = pbuffer.data(idx_dip_gg + 489);

    auto tr_z_xxxz_yyyy = pbuffer.data(idx_dip_gg + 490);

    auto tr_z_xxxz_yyyz = pbuffer.data(idx_dip_gg + 491);

    auto tr_z_xxxz_yyzz = pbuffer.data(idx_dip_gg + 492);

    auto tr_z_xxxz_yzzz = pbuffer.data(idx_dip_gg + 493);

    auto tr_z_xxxz_zzzz = pbuffer.data(idx_dip_gg + 494);

    auto tr_z_xxyy_xxxx = pbuffer.data(idx_dip_gg + 495);

    auto tr_z_xxyy_xxxy = pbuffer.data(idx_dip_gg + 496);

    auto tr_z_xxyy_xxxz = pbuffer.data(idx_dip_gg + 497);

    auto tr_z_xxyy_xxyy = pbuffer.data(idx_dip_gg + 498);

    auto tr_z_xxyy_xxyz = pbuffer.data(idx_dip_gg + 499);

    auto tr_z_xxyy_xxzz = pbuffer.data(idx_dip_gg + 500);

    auto tr_z_xxyy_xyyy = pbuffer.data(idx_dip_gg + 501);

    auto tr_z_xxyy_xyyz = pbuffer.data(idx_dip_gg + 502);

    auto tr_z_xxyy_xyzz = pbuffer.data(idx_dip_gg + 503);

    auto tr_z_xxyy_xzzz = pbuffer.data(idx_dip_gg + 504);

    auto tr_z_xxyy_yyyy = pbuffer.data(idx_dip_gg + 505);

    auto tr_z_xxyy_yyyz = pbuffer.data(idx_dip_gg + 506);

    auto tr_z_xxyy_yyzz = pbuffer.data(idx_dip_gg + 507);

    auto tr_z_xxyy_yzzz = pbuffer.data(idx_dip_gg + 508);

    auto tr_z_xxyy_zzzz = pbuffer.data(idx_dip_gg + 509);

    auto tr_z_xxyz_xxxx = pbuffer.data(idx_dip_gg + 510);

    auto tr_z_xxyz_xxxz = pbuffer.data(idx_dip_gg + 512);

    auto tr_z_xxyz_xxzz = pbuffer.data(idx_dip_gg + 515);

    auto tr_z_xxyz_xzzz = pbuffer.data(idx_dip_gg + 519);

    auto tr_z_xxyz_yyyy = pbuffer.data(idx_dip_gg + 520);

    auto tr_z_xxyz_yyyz = pbuffer.data(idx_dip_gg + 521);

    auto tr_z_xxyz_yyzz = pbuffer.data(idx_dip_gg + 522);

    auto tr_z_xxyz_yzzz = pbuffer.data(idx_dip_gg + 523);

    auto tr_z_xxzz_xxxx = pbuffer.data(idx_dip_gg + 525);

    auto tr_z_xxzz_xxxy = pbuffer.data(idx_dip_gg + 526);

    auto tr_z_xxzz_xxxz = pbuffer.data(idx_dip_gg + 527);

    auto tr_z_xxzz_xxyy = pbuffer.data(idx_dip_gg + 528);

    auto tr_z_xxzz_xxyz = pbuffer.data(idx_dip_gg + 529);

    auto tr_z_xxzz_xxzz = pbuffer.data(idx_dip_gg + 530);

    auto tr_z_xxzz_xyyy = pbuffer.data(idx_dip_gg + 531);

    auto tr_z_xxzz_xyyz = pbuffer.data(idx_dip_gg + 532);

    auto tr_z_xxzz_xyzz = pbuffer.data(idx_dip_gg + 533);

    auto tr_z_xxzz_xzzz = pbuffer.data(idx_dip_gg + 534);

    auto tr_z_xxzz_yyyy = pbuffer.data(idx_dip_gg + 535);

    auto tr_z_xxzz_yyyz = pbuffer.data(idx_dip_gg + 536);

    auto tr_z_xxzz_yyzz = pbuffer.data(idx_dip_gg + 537);

    auto tr_z_xxzz_yzzz = pbuffer.data(idx_dip_gg + 538);

    auto tr_z_xxzz_zzzz = pbuffer.data(idx_dip_gg + 539);

    auto tr_z_xyyy_xxxy = pbuffer.data(idx_dip_gg + 541);

    auto tr_z_xyyy_xxyy = pbuffer.data(idx_dip_gg + 543);

    auto tr_z_xyyy_xxyz = pbuffer.data(idx_dip_gg + 544);

    auto tr_z_xyyy_xyyy = pbuffer.data(idx_dip_gg + 546);

    auto tr_z_xyyy_xyyz = pbuffer.data(idx_dip_gg + 547);

    auto tr_z_xyyy_xyzz = pbuffer.data(idx_dip_gg + 548);

    auto tr_z_xyyy_yyyy = pbuffer.data(idx_dip_gg + 550);

    auto tr_z_xyyy_yyyz = pbuffer.data(idx_dip_gg + 551);

    auto tr_z_xyyy_yyzz = pbuffer.data(idx_dip_gg + 552);

    auto tr_z_xyyy_yzzz = pbuffer.data(idx_dip_gg + 553);

    auto tr_z_xyyy_zzzz = pbuffer.data(idx_dip_gg + 554);

    auto tr_z_xyyz_xxyz = pbuffer.data(idx_dip_gg + 559);

    auto tr_z_xyyz_xyyz = pbuffer.data(idx_dip_gg + 562);

    auto tr_z_xyyz_xyzz = pbuffer.data(idx_dip_gg + 563);

    auto tr_z_xyyz_yyyy = pbuffer.data(idx_dip_gg + 565);

    auto tr_z_xyyz_yyyz = pbuffer.data(idx_dip_gg + 566);

    auto tr_z_xyyz_yyzz = pbuffer.data(idx_dip_gg + 567);

    auto tr_z_xyyz_yzzz = pbuffer.data(idx_dip_gg + 568);

    auto tr_z_xyyz_zzzz = pbuffer.data(idx_dip_gg + 569);

    auto tr_z_xyzz_yyyy = pbuffer.data(idx_dip_gg + 580);

    auto tr_z_xyzz_yyyz = pbuffer.data(idx_dip_gg + 581);

    auto tr_z_xyzz_yyzz = pbuffer.data(idx_dip_gg + 582);

    auto tr_z_xyzz_yzzz = pbuffer.data(idx_dip_gg + 583);

    auto tr_z_xzzz_xxxx = pbuffer.data(idx_dip_gg + 585);

    auto tr_z_xzzz_xxxy = pbuffer.data(idx_dip_gg + 586);

    auto tr_z_xzzz_xxxz = pbuffer.data(idx_dip_gg + 587);

    auto tr_z_xzzz_xxyy = pbuffer.data(idx_dip_gg + 588);

    auto tr_z_xzzz_xxyz = pbuffer.data(idx_dip_gg + 589);

    auto tr_z_xzzz_xxzz = pbuffer.data(idx_dip_gg + 590);

    auto tr_z_xzzz_xyyy = pbuffer.data(idx_dip_gg + 591);

    auto tr_z_xzzz_xyyz = pbuffer.data(idx_dip_gg + 592);

    auto tr_z_xzzz_xyzz = pbuffer.data(idx_dip_gg + 593);

    auto tr_z_xzzz_xzzz = pbuffer.data(idx_dip_gg + 594);

    auto tr_z_xzzz_yyyy = pbuffer.data(idx_dip_gg + 595);

    auto tr_z_xzzz_yyyz = pbuffer.data(idx_dip_gg + 596);

    auto tr_z_xzzz_yyzz = pbuffer.data(idx_dip_gg + 597);

    auto tr_z_xzzz_yzzz = pbuffer.data(idx_dip_gg + 598);

    auto tr_z_xzzz_zzzz = pbuffer.data(idx_dip_gg + 599);

    auto tr_z_yyyy_xxxx = pbuffer.data(idx_dip_gg + 600);

    auto tr_z_yyyy_xxxy = pbuffer.data(idx_dip_gg + 601);

    auto tr_z_yyyy_xxxz = pbuffer.data(idx_dip_gg + 602);

    auto tr_z_yyyy_xxyy = pbuffer.data(idx_dip_gg + 603);

    auto tr_z_yyyy_xxyz = pbuffer.data(idx_dip_gg + 604);

    auto tr_z_yyyy_xxzz = pbuffer.data(idx_dip_gg + 605);

    auto tr_z_yyyy_xyyy = pbuffer.data(idx_dip_gg + 606);

    auto tr_z_yyyy_xyyz = pbuffer.data(idx_dip_gg + 607);

    auto tr_z_yyyy_xyzz = pbuffer.data(idx_dip_gg + 608);

    auto tr_z_yyyy_xzzz = pbuffer.data(idx_dip_gg + 609);

    auto tr_z_yyyy_yyyy = pbuffer.data(idx_dip_gg + 610);

    auto tr_z_yyyy_yyyz = pbuffer.data(idx_dip_gg + 611);

    auto tr_z_yyyy_yyzz = pbuffer.data(idx_dip_gg + 612);

    auto tr_z_yyyy_yzzz = pbuffer.data(idx_dip_gg + 613);

    auto tr_z_yyyy_zzzz = pbuffer.data(idx_dip_gg + 614);

    auto tr_z_yyyz_xxxx = pbuffer.data(idx_dip_gg + 615);

    auto tr_z_yyyz_xxxz = pbuffer.data(idx_dip_gg + 617);

    auto tr_z_yyyz_xxyz = pbuffer.data(idx_dip_gg + 619);

    auto tr_z_yyyz_xxzz = pbuffer.data(idx_dip_gg + 620);

    auto tr_z_yyyz_xyyz = pbuffer.data(idx_dip_gg + 622);

    auto tr_z_yyyz_xyzz = pbuffer.data(idx_dip_gg + 623);

    auto tr_z_yyyz_xzzz = pbuffer.data(idx_dip_gg + 624);

    auto tr_z_yyyz_yyyy = pbuffer.data(idx_dip_gg + 625);

    auto tr_z_yyyz_yyyz = pbuffer.data(idx_dip_gg + 626);

    auto tr_z_yyyz_yyzz = pbuffer.data(idx_dip_gg + 627);

    auto tr_z_yyyz_yzzz = pbuffer.data(idx_dip_gg + 628);

    auto tr_z_yyyz_zzzz = pbuffer.data(idx_dip_gg + 629);

    auto tr_z_yyzz_xxxx = pbuffer.data(idx_dip_gg + 630);

    auto tr_z_yyzz_xxxy = pbuffer.data(idx_dip_gg + 631);

    auto tr_z_yyzz_xxxz = pbuffer.data(idx_dip_gg + 632);

    auto tr_z_yyzz_xxyy = pbuffer.data(idx_dip_gg + 633);

    auto tr_z_yyzz_xxyz = pbuffer.data(idx_dip_gg + 634);

    auto tr_z_yyzz_xxzz = pbuffer.data(idx_dip_gg + 635);

    auto tr_z_yyzz_xyyy = pbuffer.data(idx_dip_gg + 636);

    auto tr_z_yyzz_xyyz = pbuffer.data(idx_dip_gg + 637);

    auto tr_z_yyzz_xyzz = pbuffer.data(idx_dip_gg + 638);

    auto tr_z_yyzz_xzzz = pbuffer.data(idx_dip_gg + 639);

    auto tr_z_yyzz_yyyy = pbuffer.data(idx_dip_gg + 640);

    auto tr_z_yyzz_yyyz = pbuffer.data(idx_dip_gg + 641);

    auto tr_z_yyzz_yyzz = pbuffer.data(idx_dip_gg + 642);

    auto tr_z_yyzz_yzzz = pbuffer.data(idx_dip_gg + 643);

    auto tr_z_yyzz_zzzz = pbuffer.data(idx_dip_gg + 644);

    auto tr_z_yzzz_xxxx = pbuffer.data(idx_dip_gg + 645);

    auto tr_z_yzzz_xxxy = pbuffer.data(idx_dip_gg + 646);

    auto tr_z_yzzz_xxxz = pbuffer.data(idx_dip_gg + 647);

    auto tr_z_yzzz_xxyy = pbuffer.data(idx_dip_gg + 648);

    auto tr_z_yzzz_xxyz = pbuffer.data(idx_dip_gg + 649);

    auto tr_z_yzzz_xxzz = pbuffer.data(idx_dip_gg + 650);

    auto tr_z_yzzz_xyyy = pbuffer.data(idx_dip_gg + 651);

    auto tr_z_yzzz_xyyz = pbuffer.data(idx_dip_gg + 652);

    auto tr_z_yzzz_xyzz = pbuffer.data(idx_dip_gg + 653);

    auto tr_z_yzzz_xzzz = pbuffer.data(idx_dip_gg + 654);

    auto tr_z_yzzz_yyyy = pbuffer.data(idx_dip_gg + 655);

    auto tr_z_yzzz_yyyz = pbuffer.data(idx_dip_gg + 656);

    auto tr_z_yzzz_yyzz = pbuffer.data(idx_dip_gg + 657);

    auto tr_z_yzzz_yzzz = pbuffer.data(idx_dip_gg + 658);

    auto tr_z_yzzz_zzzz = pbuffer.data(idx_dip_gg + 659);

    auto tr_z_zzzz_xxxx = pbuffer.data(idx_dip_gg + 660);

    auto tr_z_zzzz_xxxy = pbuffer.data(idx_dip_gg + 661);

    auto tr_z_zzzz_xxxz = pbuffer.data(idx_dip_gg + 662);

    auto tr_z_zzzz_xxyy = pbuffer.data(idx_dip_gg + 663);

    auto tr_z_zzzz_xxyz = pbuffer.data(idx_dip_gg + 664);

    auto tr_z_zzzz_xxzz = pbuffer.data(idx_dip_gg + 665);

    auto tr_z_zzzz_xyyy = pbuffer.data(idx_dip_gg + 666);

    auto tr_z_zzzz_xyyz = pbuffer.data(idx_dip_gg + 667);

    auto tr_z_zzzz_xyzz = pbuffer.data(idx_dip_gg + 668);

    auto tr_z_zzzz_xzzz = pbuffer.data(idx_dip_gg + 669);

    auto tr_z_zzzz_yyyy = pbuffer.data(idx_dip_gg + 670);

    auto tr_z_zzzz_yyyz = pbuffer.data(idx_dip_gg + 671);

    auto tr_z_zzzz_yyzz = pbuffer.data(idx_dip_gg + 672);

    auto tr_z_zzzz_yzzz = pbuffer.data(idx_dip_gg + 673);

    auto tr_z_zzzz_zzzz = pbuffer.data(idx_dip_gg + 674);

    // Set up components of auxiliary buffer : HF

    auto tr_x_xxxxx_xxx = pbuffer.data(idx_dip_hf);

    auto tr_x_xxxxx_xxy = pbuffer.data(idx_dip_hf + 1);

    auto tr_x_xxxxx_xxz = pbuffer.data(idx_dip_hf + 2);

    auto tr_x_xxxxx_xyy = pbuffer.data(idx_dip_hf + 3);

    auto tr_x_xxxxx_xyz = pbuffer.data(idx_dip_hf + 4);

    auto tr_x_xxxxx_xzz = pbuffer.data(idx_dip_hf + 5);

    auto tr_x_xxxxx_yyy = pbuffer.data(idx_dip_hf + 6);

    auto tr_x_xxxxx_yyz = pbuffer.data(idx_dip_hf + 7);

    auto tr_x_xxxxx_yzz = pbuffer.data(idx_dip_hf + 8);

    auto tr_x_xxxxx_zzz = pbuffer.data(idx_dip_hf + 9);

    auto tr_x_xxxxy_xxx = pbuffer.data(idx_dip_hf + 10);

    auto tr_x_xxxxy_xxy = pbuffer.data(idx_dip_hf + 11);

    auto tr_x_xxxxy_xxz = pbuffer.data(idx_dip_hf + 12);

    auto tr_x_xxxxy_xyy = pbuffer.data(idx_dip_hf + 13);

    auto tr_x_xxxxy_xyz = pbuffer.data(idx_dip_hf + 14);

    auto tr_x_xxxxy_xzz = pbuffer.data(idx_dip_hf + 15);

    auto tr_x_xxxxz_xxx = pbuffer.data(idx_dip_hf + 20);

    auto tr_x_xxxxz_xxy = pbuffer.data(idx_dip_hf + 21);

    auto tr_x_xxxxz_xxz = pbuffer.data(idx_dip_hf + 22);

    auto tr_x_xxxxz_xyy = pbuffer.data(idx_dip_hf + 23);

    auto tr_x_xxxxz_xyz = pbuffer.data(idx_dip_hf + 24);

    auto tr_x_xxxxz_xzz = pbuffer.data(idx_dip_hf + 25);

    auto tr_x_xxxxz_yyz = pbuffer.data(idx_dip_hf + 27);

    auto tr_x_xxxxz_yzz = pbuffer.data(idx_dip_hf + 28);

    auto tr_x_xxxxz_zzz = pbuffer.data(idx_dip_hf + 29);

    auto tr_x_xxxyy_xxx = pbuffer.data(idx_dip_hf + 30);

    auto tr_x_xxxyy_xxy = pbuffer.data(idx_dip_hf + 31);

    auto tr_x_xxxyy_xxz = pbuffer.data(idx_dip_hf + 32);

    auto tr_x_xxxyy_xyy = pbuffer.data(idx_dip_hf + 33);

    auto tr_x_xxxyy_xyz = pbuffer.data(idx_dip_hf + 34);

    auto tr_x_xxxyy_xzz = pbuffer.data(idx_dip_hf + 35);

    auto tr_x_xxxyy_yyy = pbuffer.data(idx_dip_hf + 36);

    auto tr_x_xxxyy_yyz = pbuffer.data(idx_dip_hf + 37);

    auto tr_x_xxxyy_yzz = pbuffer.data(idx_dip_hf + 38);

    auto tr_x_xxxzz_xxx = pbuffer.data(idx_dip_hf + 50);

    auto tr_x_xxxzz_xxy = pbuffer.data(idx_dip_hf + 51);

    auto tr_x_xxxzz_xxz = pbuffer.data(idx_dip_hf + 52);

    auto tr_x_xxxzz_xyy = pbuffer.data(idx_dip_hf + 53);

    auto tr_x_xxxzz_xyz = pbuffer.data(idx_dip_hf + 54);

    auto tr_x_xxxzz_xzz = pbuffer.data(idx_dip_hf + 55);

    auto tr_x_xxxzz_yyy = pbuffer.data(idx_dip_hf + 56);

    auto tr_x_xxxzz_yyz = pbuffer.data(idx_dip_hf + 57);

    auto tr_x_xxxzz_yzz = pbuffer.data(idx_dip_hf + 58);

    auto tr_x_xxxzz_zzz = pbuffer.data(idx_dip_hf + 59);

    auto tr_x_xxyyy_xxx = pbuffer.data(idx_dip_hf + 60);

    auto tr_x_xxyyy_xxy = pbuffer.data(idx_dip_hf + 61);

    auto tr_x_xxyyy_xxz = pbuffer.data(idx_dip_hf + 62);

    auto tr_x_xxyyy_xyy = pbuffer.data(idx_dip_hf + 63);

    auto tr_x_xxyyy_xyz = pbuffer.data(idx_dip_hf + 64);

    auto tr_x_xxyyy_xzz = pbuffer.data(idx_dip_hf + 65);

    auto tr_x_xxyyy_yyy = pbuffer.data(idx_dip_hf + 66);

    auto tr_x_xxyyy_yyz = pbuffer.data(idx_dip_hf + 67);

    auto tr_x_xxyyy_yzz = pbuffer.data(idx_dip_hf + 68);

    auto tr_x_xxyzz_xxz = pbuffer.data(idx_dip_hf + 82);

    auto tr_x_xxyzz_xyz = pbuffer.data(idx_dip_hf + 84);

    auto tr_x_xxyzz_xzz = pbuffer.data(idx_dip_hf + 85);

    auto tr_x_xxzzz_xxx = pbuffer.data(idx_dip_hf + 90);

    auto tr_x_xxzzz_xxy = pbuffer.data(idx_dip_hf + 91);

    auto tr_x_xxzzz_xxz = pbuffer.data(idx_dip_hf + 92);

    auto tr_x_xxzzz_xyy = pbuffer.data(idx_dip_hf + 93);

    auto tr_x_xxzzz_xyz = pbuffer.data(idx_dip_hf + 94);

    auto tr_x_xxzzz_xzz = pbuffer.data(idx_dip_hf + 95);

    auto tr_x_xxzzz_yyy = pbuffer.data(idx_dip_hf + 96);

    auto tr_x_xxzzz_yyz = pbuffer.data(idx_dip_hf + 97);

    auto tr_x_xxzzz_yzz = pbuffer.data(idx_dip_hf + 98);

    auto tr_x_xxzzz_zzz = pbuffer.data(idx_dip_hf + 99);

    auto tr_x_xyyyy_xxy = pbuffer.data(idx_dip_hf + 101);

    auto tr_x_xyyyy_xyy = pbuffer.data(idx_dip_hf + 103);

    auto tr_x_xyyyy_xyz = pbuffer.data(idx_dip_hf + 104);

    auto tr_x_xzzzz_xxx = pbuffer.data(idx_dip_hf + 140);

    auto tr_x_xzzzz_xxy = pbuffer.data(idx_dip_hf + 141);

    auto tr_x_xzzzz_xxz = pbuffer.data(idx_dip_hf + 142);

    auto tr_x_xzzzz_xyy = pbuffer.data(idx_dip_hf + 143);

    auto tr_x_xzzzz_xyz = pbuffer.data(idx_dip_hf + 144);

    auto tr_x_xzzzz_xzz = pbuffer.data(idx_dip_hf + 145);

    auto tr_x_yyyyy_xxx = pbuffer.data(idx_dip_hf + 150);

    auto tr_x_yyyyy_xxy = pbuffer.data(idx_dip_hf + 151);

    auto tr_x_yyyyy_xxz = pbuffer.data(idx_dip_hf + 152);

    auto tr_x_yyyyy_xyy = pbuffer.data(idx_dip_hf + 153);

    auto tr_x_yyyyy_xyz = pbuffer.data(idx_dip_hf + 154);

    auto tr_x_yyyyy_xzz = pbuffer.data(idx_dip_hf + 155);

    auto tr_x_yyyyy_yyy = pbuffer.data(idx_dip_hf + 156);

    auto tr_x_yyyyy_yyz = pbuffer.data(idx_dip_hf + 157);

    auto tr_x_yyyyy_yzz = pbuffer.data(idx_dip_hf + 158);

    auto tr_x_yyyyy_zzz = pbuffer.data(idx_dip_hf + 159);

    auto tr_x_yyyzz_xxz = pbuffer.data(idx_dip_hf + 172);

    auto tr_x_yyyzz_xyz = pbuffer.data(idx_dip_hf + 174);

    auto tr_x_yyyzz_xzz = pbuffer.data(idx_dip_hf + 175);

    auto tr_x_yyyzz_yyz = pbuffer.data(idx_dip_hf + 177);

    auto tr_x_yyyzz_yzz = pbuffer.data(idx_dip_hf + 178);

    auto tr_x_yyyzz_zzz = pbuffer.data(idx_dip_hf + 179);

    auto tr_x_yyzzz_xxz = pbuffer.data(idx_dip_hf + 182);

    auto tr_x_yyzzz_xyz = pbuffer.data(idx_dip_hf + 184);

    auto tr_x_yyzzz_xzz = pbuffer.data(idx_dip_hf + 185);

    auto tr_x_yyzzz_yyz = pbuffer.data(idx_dip_hf + 187);

    auto tr_x_yyzzz_yzz = pbuffer.data(idx_dip_hf + 188);

    auto tr_x_yyzzz_zzz = pbuffer.data(idx_dip_hf + 189);

    auto tr_x_yzzzz_xxz = pbuffer.data(idx_dip_hf + 192);

    auto tr_x_yzzzz_xyz = pbuffer.data(idx_dip_hf + 194);

    auto tr_x_yzzzz_xzz = pbuffer.data(idx_dip_hf + 195);

    auto tr_x_yzzzz_yyz = pbuffer.data(idx_dip_hf + 197);

    auto tr_x_yzzzz_yzz = pbuffer.data(idx_dip_hf + 198);

    auto tr_x_yzzzz_zzz = pbuffer.data(idx_dip_hf + 199);

    auto tr_x_zzzzz_xxx = pbuffer.data(idx_dip_hf + 200);

    auto tr_x_zzzzz_xxy = pbuffer.data(idx_dip_hf + 201);

    auto tr_x_zzzzz_xxz = pbuffer.data(idx_dip_hf + 202);

    auto tr_x_zzzzz_xyy = pbuffer.data(idx_dip_hf + 203);

    auto tr_x_zzzzz_xyz = pbuffer.data(idx_dip_hf + 204);

    auto tr_x_zzzzz_xzz = pbuffer.data(idx_dip_hf + 205);

    auto tr_x_zzzzz_yyy = pbuffer.data(idx_dip_hf + 206);

    auto tr_x_zzzzz_yyz = pbuffer.data(idx_dip_hf + 207);

    auto tr_x_zzzzz_yzz = pbuffer.data(idx_dip_hf + 208);

    auto tr_x_zzzzz_zzz = pbuffer.data(idx_dip_hf + 209);

    auto tr_y_xxxxx_xxx = pbuffer.data(idx_dip_hf + 210);

    auto tr_y_xxxxx_xxy = pbuffer.data(idx_dip_hf + 211);

    auto tr_y_xxxxx_xxz = pbuffer.data(idx_dip_hf + 212);

    auto tr_y_xxxxx_xyy = pbuffer.data(idx_dip_hf + 213);

    auto tr_y_xxxxx_xyz = pbuffer.data(idx_dip_hf + 214);

    auto tr_y_xxxxx_xzz = pbuffer.data(idx_dip_hf + 215);

    auto tr_y_xxxxx_yyy = pbuffer.data(idx_dip_hf + 216);

    auto tr_y_xxxxx_yyz = pbuffer.data(idx_dip_hf + 217);

    auto tr_y_xxxxx_yzz = pbuffer.data(idx_dip_hf + 218);

    auto tr_y_xxxxx_zzz = pbuffer.data(idx_dip_hf + 219);

    auto tr_y_xxxxy_xxy = pbuffer.data(idx_dip_hf + 221);

    auto tr_y_xxxxy_xyy = pbuffer.data(idx_dip_hf + 223);

    auto tr_y_xxxxy_xyz = pbuffer.data(idx_dip_hf + 224);

    auto tr_y_xxxxy_yyy = pbuffer.data(idx_dip_hf + 226);

    auto tr_y_xxxxy_yyz = pbuffer.data(idx_dip_hf + 227);

    auto tr_y_xxxxy_yzz = pbuffer.data(idx_dip_hf + 228);

    auto tr_y_xxxyy_xxx = pbuffer.data(idx_dip_hf + 240);

    auto tr_y_xxxyy_xxy = pbuffer.data(idx_dip_hf + 241);

    auto tr_y_xxxyy_xxz = pbuffer.data(idx_dip_hf + 242);

    auto tr_y_xxxyy_xyy = pbuffer.data(idx_dip_hf + 243);

    auto tr_y_xxxyy_xyz = pbuffer.data(idx_dip_hf + 244);

    auto tr_y_xxxyy_xzz = pbuffer.data(idx_dip_hf + 245);

    auto tr_y_xxxyy_yyy = pbuffer.data(idx_dip_hf + 246);

    auto tr_y_xxxyy_yyz = pbuffer.data(idx_dip_hf + 247);

    auto tr_y_xxxyy_yzz = pbuffer.data(idx_dip_hf + 248);

    auto tr_y_xxxyy_zzz = pbuffer.data(idx_dip_hf + 249);

    auto tr_y_xxxzz_xxz = pbuffer.data(idx_dip_hf + 262);

    auto tr_y_xxxzz_xyz = pbuffer.data(idx_dip_hf + 264);

    auto tr_y_xxxzz_xzz = pbuffer.data(idx_dip_hf + 265);

    auto tr_y_xxxzz_yyz = pbuffer.data(idx_dip_hf + 267);

    auto tr_y_xxxzz_yzz = pbuffer.data(idx_dip_hf + 268);

    auto tr_y_xxxzz_zzz = pbuffer.data(idx_dip_hf + 269);

    auto tr_y_xxyyy_xxx = pbuffer.data(idx_dip_hf + 270);

    auto tr_y_xxyyy_xxy = pbuffer.data(idx_dip_hf + 271);

    auto tr_y_xxyyy_xxz = pbuffer.data(idx_dip_hf + 272);

    auto tr_y_xxyyy_xyy = pbuffer.data(idx_dip_hf + 273);

    auto tr_y_xxyyy_xyz = pbuffer.data(idx_dip_hf + 274);

    auto tr_y_xxyyy_xzz = pbuffer.data(idx_dip_hf + 275);

    auto tr_y_xxyyy_yyy = pbuffer.data(idx_dip_hf + 276);

    auto tr_y_xxyyy_yyz = pbuffer.data(idx_dip_hf + 277);

    auto tr_y_xxyyy_yzz = pbuffer.data(idx_dip_hf + 278);

    auto tr_y_xxyyy_zzz = pbuffer.data(idx_dip_hf + 279);

    auto tr_y_xxyzz_xyz = pbuffer.data(idx_dip_hf + 294);

    auto tr_y_xxyzz_yyz = pbuffer.data(idx_dip_hf + 297);

    auto tr_y_xxyzz_yzz = pbuffer.data(idx_dip_hf + 298);

    auto tr_y_xxzzz_xxz = pbuffer.data(idx_dip_hf + 302);

    auto tr_y_xxzzz_xyz = pbuffer.data(idx_dip_hf + 304);

    auto tr_y_xxzzz_xzz = pbuffer.data(idx_dip_hf + 305);

    auto tr_y_xxzzz_yyz = pbuffer.data(idx_dip_hf + 307);

    auto tr_y_xxzzz_yzz = pbuffer.data(idx_dip_hf + 308);

    auto tr_y_xxzzz_zzz = pbuffer.data(idx_dip_hf + 309);

    auto tr_y_xyyyy_xxx = pbuffer.data(idx_dip_hf + 310);

    auto tr_y_xyyyy_xxy = pbuffer.data(idx_dip_hf + 311);

    auto tr_y_xyyyy_xxz = pbuffer.data(idx_dip_hf + 312);

    auto tr_y_xyyyy_xyy = pbuffer.data(idx_dip_hf + 313);

    auto tr_y_xyyyy_xyz = pbuffer.data(idx_dip_hf + 314);

    auto tr_y_xyyyy_xzz = pbuffer.data(idx_dip_hf + 315);

    auto tr_y_xyyyy_yyy = pbuffer.data(idx_dip_hf + 316);

    auto tr_y_xyyyy_yyz = pbuffer.data(idx_dip_hf + 317);

    auto tr_y_xyyyy_yzz = pbuffer.data(idx_dip_hf + 318);

    auto tr_y_xyyyy_zzz = pbuffer.data(idx_dip_hf + 319);

    auto tr_y_xyyzz_xxz = pbuffer.data(idx_dip_hf + 332);

    auto tr_y_xyyzz_xyz = pbuffer.data(idx_dip_hf + 334);

    auto tr_y_xyyzz_xzz = pbuffer.data(idx_dip_hf + 335);

    auto tr_y_xyyzz_yyz = pbuffer.data(idx_dip_hf + 337);

    auto tr_y_xyyzz_yzz = pbuffer.data(idx_dip_hf + 338);

    auto tr_y_xyyzz_zzz = pbuffer.data(idx_dip_hf + 339);

    auto tr_y_xyzzz_xyz = pbuffer.data(idx_dip_hf + 344);

    auto tr_y_xyzzz_yyz = pbuffer.data(idx_dip_hf + 347);

    auto tr_y_xyzzz_yzz = pbuffer.data(idx_dip_hf + 348);

    auto tr_y_xzzzz_xxz = pbuffer.data(idx_dip_hf + 352);

    auto tr_y_xzzzz_xyz = pbuffer.data(idx_dip_hf + 354);

    auto tr_y_xzzzz_xzz = pbuffer.data(idx_dip_hf + 355);

    auto tr_y_xzzzz_yyz = pbuffer.data(idx_dip_hf + 357);

    auto tr_y_xzzzz_yzz = pbuffer.data(idx_dip_hf + 358);

    auto tr_y_xzzzz_zzz = pbuffer.data(idx_dip_hf + 359);

    auto tr_y_yyyyy_xxx = pbuffer.data(idx_dip_hf + 360);

    auto tr_y_yyyyy_xxy = pbuffer.data(idx_dip_hf + 361);

    auto tr_y_yyyyy_xxz = pbuffer.data(idx_dip_hf + 362);

    auto tr_y_yyyyy_xyy = pbuffer.data(idx_dip_hf + 363);

    auto tr_y_yyyyy_xyz = pbuffer.data(idx_dip_hf + 364);

    auto tr_y_yyyyy_xzz = pbuffer.data(idx_dip_hf + 365);

    auto tr_y_yyyyy_yyy = pbuffer.data(idx_dip_hf + 366);

    auto tr_y_yyyyy_yyz = pbuffer.data(idx_dip_hf + 367);

    auto tr_y_yyyyy_yzz = pbuffer.data(idx_dip_hf + 368);

    auto tr_y_yyyyy_zzz = pbuffer.data(idx_dip_hf + 369);

    auto tr_y_yyyyz_xxy = pbuffer.data(idx_dip_hf + 371);

    auto tr_y_yyyyz_xxz = pbuffer.data(idx_dip_hf + 372);

    auto tr_y_yyyyz_xyy = pbuffer.data(idx_dip_hf + 373);

    auto tr_y_yyyyz_xyz = pbuffer.data(idx_dip_hf + 374);

    auto tr_y_yyyyz_xzz = pbuffer.data(idx_dip_hf + 375);

    auto tr_y_yyyyz_yyy = pbuffer.data(idx_dip_hf + 376);

    auto tr_y_yyyyz_yyz = pbuffer.data(idx_dip_hf + 377);

    auto tr_y_yyyyz_yzz = pbuffer.data(idx_dip_hf + 378);

    auto tr_y_yyyyz_zzz = pbuffer.data(idx_dip_hf + 379);

    auto tr_y_yyyzz_xxx = pbuffer.data(idx_dip_hf + 380);

    auto tr_y_yyyzz_xxy = pbuffer.data(idx_dip_hf + 381);

    auto tr_y_yyyzz_xxz = pbuffer.data(idx_dip_hf + 382);

    auto tr_y_yyyzz_xyy = pbuffer.data(idx_dip_hf + 383);

    auto tr_y_yyyzz_xyz = pbuffer.data(idx_dip_hf + 384);

    auto tr_y_yyyzz_xzz = pbuffer.data(idx_dip_hf + 385);

    auto tr_y_yyyzz_yyy = pbuffer.data(idx_dip_hf + 386);

    auto tr_y_yyyzz_yyz = pbuffer.data(idx_dip_hf + 387);

    auto tr_y_yyyzz_yzz = pbuffer.data(idx_dip_hf + 388);

    auto tr_y_yyyzz_zzz = pbuffer.data(idx_dip_hf + 389);

    auto tr_y_yyzzz_xxx = pbuffer.data(idx_dip_hf + 390);

    auto tr_y_yyzzz_xxy = pbuffer.data(idx_dip_hf + 391);

    auto tr_y_yyzzz_xxz = pbuffer.data(idx_dip_hf + 392);

    auto tr_y_yyzzz_xyy = pbuffer.data(idx_dip_hf + 393);

    auto tr_y_yyzzz_xyz = pbuffer.data(idx_dip_hf + 394);

    auto tr_y_yyzzz_xzz = pbuffer.data(idx_dip_hf + 395);

    auto tr_y_yyzzz_yyy = pbuffer.data(idx_dip_hf + 396);

    auto tr_y_yyzzz_yyz = pbuffer.data(idx_dip_hf + 397);

    auto tr_y_yyzzz_yzz = pbuffer.data(idx_dip_hf + 398);

    auto tr_y_yyzzz_zzz = pbuffer.data(idx_dip_hf + 399);

    auto tr_y_yzzzz_xxx = pbuffer.data(idx_dip_hf + 400);

    auto tr_y_yzzzz_xxy = pbuffer.data(idx_dip_hf + 401);

    auto tr_y_yzzzz_xxz = pbuffer.data(idx_dip_hf + 402);

    auto tr_y_yzzzz_xyy = pbuffer.data(idx_dip_hf + 403);

    auto tr_y_yzzzz_xyz = pbuffer.data(idx_dip_hf + 404);

    auto tr_y_yzzzz_xzz = pbuffer.data(idx_dip_hf + 405);

    auto tr_y_yzzzz_yyy = pbuffer.data(idx_dip_hf + 406);

    auto tr_y_yzzzz_yyz = pbuffer.data(idx_dip_hf + 407);

    auto tr_y_yzzzz_yzz = pbuffer.data(idx_dip_hf + 408);

    auto tr_y_yzzzz_zzz = pbuffer.data(idx_dip_hf + 409);

    auto tr_y_zzzzz_xxx = pbuffer.data(idx_dip_hf + 410);

    auto tr_y_zzzzz_xxy = pbuffer.data(idx_dip_hf + 411);

    auto tr_y_zzzzz_xxz = pbuffer.data(idx_dip_hf + 412);

    auto tr_y_zzzzz_xyy = pbuffer.data(idx_dip_hf + 413);

    auto tr_y_zzzzz_xyz = pbuffer.data(idx_dip_hf + 414);

    auto tr_y_zzzzz_xzz = pbuffer.data(idx_dip_hf + 415);

    auto tr_y_zzzzz_yyy = pbuffer.data(idx_dip_hf + 416);

    auto tr_y_zzzzz_yyz = pbuffer.data(idx_dip_hf + 417);

    auto tr_y_zzzzz_yzz = pbuffer.data(idx_dip_hf + 418);

    auto tr_y_zzzzz_zzz = pbuffer.data(idx_dip_hf + 419);

    auto tr_z_xxxxx_xxx = pbuffer.data(idx_dip_hf + 420);

    auto tr_z_xxxxx_xxy = pbuffer.data(idx_dip_hf + 421);

    auto tr_z_xxxxx_xxz = pbuffer.data(idx_dip_hf + 422);

    auto tr_z_xxxxx_xyy = pbuffer.data(idx_dip_hf + 423);

    auto tr_z_xxxxx_xyz = pbuffer.data(idx_dip_hf + 424);

    auto tr_z_xxxxx_xzz = pbuffer.data(idx_dip_hf + 425);

    auto tr_z_xxxxx_yyy = pbuffer.data(idx_dip_hf + 426);

    auto tr_z_xxxxx_yyz = pbuffer.data(idx_dip_hf + 427);

    auto tr_z_xxxxx_yzz = pbuffer.data(idx_dip_hf + 428);

    auto tr_z_xxxxx_zzz = pbuffer.data(idx_dip_hf + 429);

    auto tr_z_xxxxz_xxx = pbuffer.data(idx_dip_hf + 440);

    auto tr_z_xxxxz_xxy = pbuffer.data(idx_dip_hf + 441);

    auto tr_z_xxxxz_xxz = pbuffer.data(idx_dip_hf + 442);

    auto tr_z_xxxxz_xyy = pbuffer.data(idx_dip_hf + 443);

    auto tr_z_xxxxz_xyz = pbuffer.data(idx_dip_hf + 444);

    auto tr_z_xxxxz_xzz = pbuffer.data(idx_dip_hf + 445);

    auto tr_z_xxxxz_yyz = pbuffer.data(idx_dip_hf + 447);

    auto tr_z_xxxxz_yzz = pbuffer.data(idx_dip_hf + 448);

    auto tr_z_xxxxz_zzz = pbuffer.data(idx_dip_hf + 449);

    auto tr_z_xxxyy_xxy = pbuffer.data(idx_dip_hf + 451);

    auto tr_z_xxxyy_xyy = pbuffer.data(idx_dip_hf + 453);

    auto tr_z_xxxyy_xyz = pbuffer.data(idx_dip_hf + 454);

    auto tr_z_xxxyy_yyy = pbuffer.data(idx_dip_hf + 456);

    auto tr_z_xxxyy_yyz = pbuffer.data(idx_dip_hf + 457);

    auto tr_z_xxxyy_yzz = pbuffer.data(idx_dip_hf + 458);

    auto tr_z_xxxzz_xxx = pbuffer.data(idx_dip_hf + 470);

    auto tr_z_xxxzz_xxy = pbuffer.data(idx_dip_hf + 471);

    auto tr_z_xxxzz_xxz = pbuffer.data(idx_dip_hf + 472);

    auto tr_z_xxxzz_xyy = pbuffer.data(idx_dip_hf + 473);

    auto tr_z_xxxzz_xyz = pbuffer.data(idx_dip_hf + 474);

    auto tr_z_xxxzz_xzz = pbuffer.data(idx_dip_hf + 475);

    auto tr_z_xxxzz_yyy = pbuffer.data(idx_dip_hf + 476);

    auto tr_z_xxxzz_yyz = pbuffer.data(idx_dip_hf + 477);

    auto tr_z_xxxzz_yzz = pbuffer.data(idx_dip_hf + 478);

    auto tr_z_xxxzz_zzz = pbuffer.data(idx_dip_hf + 479);

    auto tr_z_xxyyy_xxy = pbuffer.data(idx_dip_hf + 481);

    auto tr_z_xxyyy_xyy = pbuffer.data(idx_dip_hf + 483);

    auto tr_z_xxyyy_xyz = pbuffer.data(idx_dip_hf + 484);

    auto tr_z_xxyyy_yyy = pbuffer.data(idx_dip_hf + 486);

    auto tr_z_xxyyy_yyz = pbuffer.data(idx_dip_hf + 487);

    auto tr_z_xxyyy_yzz = pbuffer.data(idx_dip_hf + 488);

    auto tr_z_xxyyz_xyz = pbuffer.data(idx_dip_hf + 494);

    auto tr_z_xxyyz_yyz = pbuffer.data(idx_dip_hf + 497);

    auto tr_z_xxyyz_yzz = pbuffer.data(idx_dip_hf + 498);

    auto tr_z_xxzzz_xxx = pbuffer.data(idx_dip_hf + 510);

    auto tr_z_xxzzz_xxy = pbuffer.data(idx_dip_hf + 511);

    auto tr_z_xxzzz_xxz = pbuffer.data(idx_dip_hf + 512);

    auto tr_z_xxzzz_xyy = pbuffer.data(idx_dip_hf + 513);

    auto tr_z_xxzzz_xyz = pbuffer.data(idx_dip_hf + 514);

    auto tr_z_xxzzz_xzz = pbuffer.data(idx_dip_hf + 515);

    auto tr_z_xxzzz_yyy = pbuffer.data(idx_dip_hf + 516);

    auto tr_z_xxzzz_yyz = pbuffer.data(idx_dip_hf + 517);

    auto tr_z_xxzzz_yzz = pbuffer.data(idx_dip_hf + 518);

    auto tr_z_xxzzz_zzz = pbuffer.data(idx_dip_hf + 519);

    auto tr_z_xyyyy_xxy = pbuffer.data(idx_dip_hf + 521);

    auto tr_z_xyyyy_xyy = pbuffer.data(idx_dip_hf + 523);

    auto tr_z_xyyyy_xyz = pbuffer.data(idx_dip_hf + 524);

    auto tr_z_xyyyy_yyy = pbuffer.data(idx_dip_hf + 526);

    auto tr_z_xyyyy_yyz = pbuffer.data(idx_dip_hf + 527);

    auto tr_z_xyyyy_yzz = pbuffer.data(idx_dip_hf + 528);

    auto tr_z_xyyyz_xyz = pbuffer.data(idx_dip_hf + 534);

    auto tr_z_xyyyz_yyz = pbuffer.data(idx_dip_hf + 537);

    auto tr_z_xyyyz_yzz = pbuffer.data(idx_dip_hf + 538);

    auto tr_z_xyyzz_xxy = pbuffer.data(idx_dip_hf + 541);

    auto tr_z_xyyzz_xyy = pbuffer.data(idx_dip_hf + 543);

    auto tr_z_xyyzz_xyz = pbuffer.data(idx_dip_hf + 544);

    auto tr_z_xyyzz_yyy = pbuffer.data(idx_dip_hf + 546);

    auto tr_z_xyyzz_yyz = pbuffer.data(idx_dip_hf + 547);

    auto tr_z_xyyzz_yzz = pbuffer.data(idx_dip_hf + 548);

    auto tr_z_xzzzz_xxx = pbuffer.data(idx_dip_hf + 560);

    auto tr_z_xzzzz_xxy = pbuffer.data(idx_dip_hf + 561);

    auto tr_z_xzzzz_xxz = pbuffer.data(idx_dip_hf + 562);

    auto tr_z_xzzzz_xyy = pbuffer.data(idx_dip_hf + 563);

    auto tr_z_xzzzz_xyz = pbuffer.data(idx_dip_hf + 564);

    auto tr_z_xzzzz_xzz = pbuffer.data(idx_dip_hf + 565);

    auto tr_z_xzzzz_yyy = pbuffer.data(idx_dip_hf + 566);

    auto tr_z_xzzzz_yyz = pbuffer.data(idx_dip_hf + 567);

    auto tr_z_xzzzz_yzz = pbuffer.data(idx_dip_hf + 568);

    auto tr_z_xzzzz_zzz = pbuffer.data(idx_dip_hf + 569);

    auto tr_z_yyyyy_xxx = pbuffer.data(idx_dip_hf + 570);

    auto tr_z_yyyyy_xxy = pbuffer.data(idx_dip_hf + 571);

    auto tr_z_yyyyy_xxz = pbuffer.data(idx_dip_hf + 572);

    auto tr_z_yyyyy_xyy = pbuffer.data(idx_dip_hf + 573);

    auto tr_z_yyyyy_xyz = pbuffer.data(idx_dip_hf + 574);

    auto tr_z_yyyyy_xzz = pbuffer.data(idx_dip_hf + 575);

    auto tr_z_yyyyy_yyy = pbuffer.data(idx_dip_hf + 576);

    auto tr_z_yyyyy_yyz = pbuffer.data(idx_dip_hf + 577);

    auto tr_z_yyyyy_yzz = pbuffer.data(idx_dip_hf + 578);

    auto tr_z_yyyyy_zzz = pbuffer.data(idx_dip_hf + 579);

    auto tr_z_yyyyz_xxx = pbuffer.data(idx_dip_hf + 580);

    auto tr_z_yyyyz_xxy = pbuffer.data(idx_dip_hf + 581);

    auto tr_z_yyyyz_xxz = pbuffer.data(idx_dip_hf + 582);

    auto tr_z_yyyyz_xyy = pbuffer.data(idx_dip_hf + 583);

    auto tr_z_yyyyz_xyz = pbuffer.data(idx_dip_hf + 584);

    auto tr_z_yyyyz_xzz = pbuffer.data(idx_dip_hf + 585);

    auto tr_z_yyyyz_yyy = pbuffer.data(idx_dip_hf + 586);

    auto tr_z_yyyyz_yyz = pbuffer.data(idx_dip_hf + 587);

    auto tr_z_yyyyz_yzz = pbuffer.data(idx_dip_hf + 588);

    auto tr_z_yyyyz_zzz = pbuffer.data(idx_dip_hf + 589);

    auto tr_z_yyyzz_xxx = pbuffer.data(idx_dip_hf + 590);

    auto tr_z_yyyzz_xxy = pbuffer.data(idx_dip_hf + 591);

    auto tr_z_yyyzz_xxz = pbuffer.data(idx_dip_hf + 592);

    auto tr_z_yyyzz_xyy = pbuffer.data(idx_dip_hf + 593);

    auto tr_z_yyyzz_xyz = pbuffer.data(idx_dip_hf + 594);

    auto tr_z_yyyzz_xzz = pbuffer.data(idx_dip_hf + 595);

    auto tr_z_yyyzz_yyy = pbuffer.data(idx_dip_hf + 596);

    auto tr_z_yyyzz_yyz = pbuffer.data(idx_dip_hf + 597);

    auto tr_z_yyyzz_yzz = pbuffer.data(idx_dip_hf + 598);

    auto tr_z_yyyzz_zzz = pbuffer.data(idx_dip_hf + 599);

    auto tr_z_yyzzz_xxx = pbuffer.data(idx_dip_hf + 600);

    auto tr_z_yyzzz_xxy = pbuffer.data(idx_dip_hf + 601);

    auto tr_z_yyzzz_xxz = pbuffer.data(idx_dip_hf + 602);

    auto tr_z_yyzzz_xyy = pbuffer.data(idx_dip_hf + 603);

    auto tr_z_yyzzz_xyz = pbuffer.data(idx_dip_hf + 604);

    auto tr_z_yyzzz_xzz = pbuffer.data(idx_dip_hf + 605);

    auto tr_z_yyzzz_yyy = pbuffer.data(idx_dip_hf + 606);

    auto tr_z_yyzzz_yyz = pbuffer.data(idx_dip_hf + 607);

    auto tr_z_yyzzz_yzz = pbuffer.data(idx_dip_hf + 608);

    auto tr_z_yyzzz_zzz = pbuffer.data(idx_dip_hf + 609);

    auto tr_z_yzzzz_xxx = pbuffer.data(idx_dip_hf + 610);

    auto tr_z_yzzzz_xxy = pbuffer.data(idx_dip_hf + 611);

    auto tr_z_yzzzz_xxz = pbuffer.data(idx_dip_hf + 612);

    auto tr_z_yzzzz_xyy = pbuffer.data(idx_dip_hf + 613);

    auto tr_z_yzzzz_xyz = pbuffer.data(idx_dip_hf + 614);

    auto tr_z_yzzzz_xzz = pbuffer.data(idx_dip_hf + 615);

    auto tr_z_yzzzz_yyy = pbuffer.data(idx_dip_hf + 616);

    auto tr_z_yzzzz_yyz = pbuffer.data(idx_dip_hf + 617);

    auto tr_z_yzzzz_yzz = pbuffer.data(idx_dip_hf + 618);

    auto tr_z_yzzzz_zzz = pbuffer.data(idx_dip_hf + 619);

    auto tr_z_zzzzz_xxx = pbuffer.data(idx_dip_hf + 620);

    auto tr_z_zzzzz_xxy = pbuffer.data(idx_dip_hf + 621);

    auto tr_z_zzzzz_xxz = pbuffer.data(idx_dip_hf + 622);

    auto tr_z_zzzzz_xyy = pbuffer.data(idx_dip_hf + 623);

    auto tr_z_zzzzz_xyz = pbuffer.data(idx_dip_hf + 624);

    auto tr_z_zzzzz_xzz = pbuffer.data(idx_dip_hf + 625);

    auto tr_z_zzzzz_yyy = pbuffer.data(idx_dip_hf + 626);

    auto tr_z_zzzzz_yyz = pbuffer.data(idx_dip_hf + 627);

    auto tr_z_zzzzz_yzz = pbuffer.data(idx_dip_hf + 628);

    auto tr_z_zzzzz_zzz = pbuffer.data(idx_dip_hf + 629);

    // Set up components of auxiliary buffer : HG

    auto ts_xxxxx_xxxx = pbuffer.data(idx_ovl_hg);

    auto ts_xxxxx_xxxy = pbuffer.data(idx_ovl_hg + 1);

    auto ts_xxxxx_xxxz = pbuffer.data(idx_ovl_hg + 2);

    auto ts_xxxxx_xxyy = pbuffer.data(idx_ovl_hg + 3);

    auto ts_xxxxx_xxyz = pbuffer.data(idx_ovl_hg + 4);

    auto ts_xxxxx_xxzz = pbuffer.data(idx_ovl_hg + 5);

    auto ts_xxxxx_xyyy = pbuffer.data(idx_ovl_hg + 6);

    auto ts_xxxxx_xyyz = pbuffer.data(idx_ovl_hg + 7);

    auto ts_xxxxx_xyzz = pbuffer.data(idx_ovl_hg + 8);

    auto ts_xxxxx_xzzz = pbuffer.data(idx_ovl_hg + 9);

    auto ts_xxxxx_yyyy = pbuffer.data(idx_ovl_hg + 10);

    auto ts_xxxxx_yyyz = pbuffer.data(idx_ovl_hg + 11);

    auto ts_xxxxx_yyzz = pbuffer.data(idx_ovl_hg + 12);

    auto ts_xxxxx_yzzz = pbuffer.data(idx_ovl_hg + 13);

    auto ts_xxxxx_zzzz = pbuffer.data(idx_ovl_hg + 14);

    auto ts_xxxxz_xxxz = pbuffer.data(idx_ovl_hg + 32);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_ovl_hg + 35);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_ovl_hg + 39);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_ovl_hg + 46);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_ovl_hg + 48);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_ovl_hg + 51);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_ovl_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_ovl_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_ovl_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_ovl_hg + 58);

    auto ts_xxxzz_xxxx = pbuffer.data(idx_ovl_hg + 75);

    auto ts_xxxzz_xxxz = pbuffer.data(idx_ovl_hg + 77);

    auto ts_xxxzz_xxzz = pbuffer.data(idx_ovl_hg + 80);

    auto ts_xxxzz_xzzz = pbuffer.data(idx_ovl_hg + 84);

    auto ts_xxxzz_yyyz = pbuffer.data(idx_ovl_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_ovl_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_ovl_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_ovl_hg + 89);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_ovl_hg + 91);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_ovl_hg + 93);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_ovl_hg + 96);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_ovl_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_ovl_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_ovl_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_ovl_hg + 103);

    auto ts_xxzzz_xxxx = pbuffer.data(idx_ovl_hg + 135);

    auto ts_xxzzz_xxxz = pbuffer.data(idx_ovl_hg + 137);

    auto ts_xxzzz_xxzz = pbuffer.data(idx_ovl_hg + 140);

    auto ts_xxzzz_xzzz = pbuffer.data(idx_ovl_hg + 144);

    auto ts_xxzzz_yyyz = pbuffer.data(idx_ovl_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_ovl_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_ovl_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_ovl_hg + 149);

    auto ts_xyyyy_yyyy = pbuffer.data(idx_ovl_hg + 160);

    auto ts_xyyyy_yyyz = pbuffer.data(idx_ovl_hg + 161);

    auto ts_xyyyy_yyzz = pbuffer.data(idx_ovl_hg + 162);

    auto ts_xyyyy_yzzz = pbuffer.data(idx_ovl_hg + 163);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_ovl_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_ovl_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_ovl_hg + 193);

    auto ts_xzzzz_yyyz = pbuffer.data(idx_ovl_hg + 221);

    auto ts_xzzzz_yyzz = pbuffer.data(idx_ovl_hg + 222);

    auto ts_xzzzz_yzzz = pbuffer.data(idx_ovl_hg + 223);

    auto ts_xzzzz_zzzz = pbuffer.data(idx_ovl_hg + 224);

    auto ts_yyyyy_xxxx = pbuffer.data(idx_ovl_hg + 225);

    auto ts_yyyyy_xxxy = pbuffer.data(idx_ovl_hg + 226);

    auto ts_yyyyy_xxxz = pbuffer.data(idx_ovl_hg + 227);

    auto ts_yyyyy_xxyy = pbuffer.data(idx_ovl_hg + 228);

    auto ts_yyyyy_xxyz = pbuffer.data(idx_ovl_hg + 229);

    auto ts_yyyyy_xxzz = pbuffer.data(idx_ovl_hg + 230);

    auto ts_yyyyy_xyyy = pbuffer.data(idx_ovl_hg + 231);

    auto ts_yyyyy_xyyz = pbuffer.data(idx_ovl_hg + 232);

    auto ts_yyyyy_xyzz = pbuffer.data(idx_ovl_hg + 233);

    auto ts_yyyyy_xzzz = pbuffer.data(idx_ovl_hg + 234);

    auto ts_yyyyy_yyyy = pbuffer.data(idx_ovl_hg + 235);

    auto ts_yyyyy_yyyz = pbuffer.data(idx_ovl_hg + 236);

    auto ts_yyyyy_yyzz = pbuffer.data(idx_ovl_hg + 237);

    auto ts_yyyyy_yzzz = pbuffer.data(idx_ovl_hg + 238);

    auto ts_yyyyy_zzzz = pbuffer.data(idx_ovl_hg + 239);

    auto ts_yyyyz_yyyz = pbuffer.data(idx_ovl_hg + 251);

    auto ts_yyyyz_yyzz = pbuffer.data(idx_ovl_hg + 252);

    auto ts_yyyyz_yzzz = pbuffer.data(idx_ovl_hg + 253);

    auto ts_yyyyz_zzzz = pbuffer.data(idx_ovl_hg + 254);

    auto ts_yyyzz_xxxz = pbuffer.data(idx_ovl_hg + 257);

    auto ts_yyyzz_xxyz = pbuffer.data(idx_ovl_hg + 259);

    auto ts_yyyzz_xxzz = pbuffer.data(idx_ovl_hg + 260);

    auto ts_yyyzz_xyyz = pbuffer.data(idx_ovl_hg + 262);

    auto ts_yyyzz_xyzz = pbuffer.data(idx_ovl_hg + 263);

    auto ts_yyyzz_xzzz = pbuffer.data(idx_ovl_hg + 264);

    auto ts_yyyzz_yyyy = pbuffer.data(idx_ovl_hg + 265);

    auto ts_yyyzz_yyyz = pbuffer.data(idx_ovl_hg + 266);

    auto ts_yyyzz_yyzz = pbuffer.data(idx_ovl_hg + 267);

    auto ts_yyyzz_yzzz = pbuffer.data(idx_ovl_hg + 268);

    auto ts_yyyzz_zzzz = pbuffer.data(idx_ovl_hg + 269);

    auto ts_yyzzz_xxxz = pbuffer.data(idx_ovl_hg + 272);

    auto ts_yyzzz_xxyz = pbuffer.data(idx_ovl_hg + 274);

    auto ts_yyzzz_xxzz = pbuffer.data(idx_ovl_hg + 275);

    auto ts_yyzzz_xyyz = pbuffer.data(idx_ovl_hg + 277);

    auto ts_yyzzz_xyzz = pbuffer.data(idx_ovl_hg + 278);

    auto ts_yyzzz_xzzz = pbuffer.data(idx_ovl_hg + 279);

    auto ts_yyzzz_yyyy = pbuffer.data(idx_ovl_hg + 280);

    auto ts_yyzzz_yyyz = pbuffer.data(idx_ovl_hg + 281);

    auto ts_yyzzz_yyzz = pbuffer.data(idx_ovl_hg + 282);

    auto ts_yyzzz_yzzz = pbuffer.data(idx_ovl_hg + 283);

    auto ts_yyzzz_zzzz = pbuffer.data(idx_ovl_hg + 284);

    auto ts_yzzzz_xxxz = pbuffer.data(idx_ovl_hg + 287);

    auto ts_yzzzz_xxzz = pbuffer.data(idx_ovl_hg + 290);

    auto ts_yzzzz_xzzz = pbuffer.data(idx_ovl_hg + 294);

    auto ts_yzzzz_yyyy = pbuffer.data(idx_ovl_hg + 295);

    auto ts_yzzzz_yyyz = pbuffer.data(idx_ovl_hg + 296);

    auto ts_yzzzz_yyzz = pbuffer.data(idx_ovl_hg + 297);

    auto ts_yzzzz_yzzz = pbuffer.data(idx_ovl_hg + 298);

    auto ts_yzzzz_zzzz = pbuffer.data(idx_ovl_hg + 299);

    auto ts_zzzzz_xxxx = pbuffer.data(idx_ovl_hg + 300);

    auto ts_zzzzz_xxxy = pbuffer.data(idx_ovl_hg + 301);

    auto ts_zzzzz_xxxz = pbuffer.data(idx_ovl_hg + 302);

    auto ts_zzzzz_xxyy = pbuffer.data(idx_ovl_hg + 303);

    auto ts_zzzzz_xxyz = pbuffer.data(idx_ovl_hg + 304);

    auto ts_zzzzz_xxzz = pbuffer.data(idx_ovl_hg + 305);

    auto ts_zzzzz_xyyy = pbuffer.data(idx_ovl_hg + 306);

    auto ts_zzzzz_xyyz = pbuffer.data(idx_ovl_hg + 307);

    auto ts_zzzzz_xyzz = pbuffer.data(idx_ovl_hg + 308);

    auto ts_zzzzz_xzzz = pbuffer.data(idx_ovl_hg + 309);

    auto ts_zzzzz_yyyy = pbuffer.data(idx_ovl_hg + 310);

    auto ts_zzzzz_yyyz = pbuffer.data(idx_ovl_hg + 311);

    auto ts_zzzzz_yyzz = pbuffer.data(idx_ovl_hg + 312);

    auto ts_zzzzz_yzzz = pbuffer.data(idx_ovl_hg + 313);

    auto ts_zzzzz_zzzz = pbuffer.data(idx_ovl_hg + 314);

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

    auto tr_x_xxxxy_yyyy = pbuffer.data(idx_dip_hg + 25);

    auto tr_x_xxxxy_zzzz = pbuffer.data(idx_dip_hg + 29);

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

    auto tr_x_xxxxz_yyyy = pbuffer.data(idx_dip_hg + 40);

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

    auto tr_x_xxxyy_zzzz = pbuffer.data(idx_dip_hg + 59);

    auto tr_x_xxxyz_xxxz = pbuffer.data(idx_dip_hg + 62);

    auto tr_x_xxxyz_xxzz = pbuffer.data(idx_dip_hg + 65);

    auto tr_x_xxxyz_xzzz = pbuffer.data(idx_dip_hg + 69);

    auto tr_x_xxxyz_zzzz = pbuffer.data(idx_dip_hg + 74);

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

    auto tr_x_xxyyy_zzzz = pbuffer.data(idx_dip_hg + 104);

    auto tr_x_xxyyz_xxxy = pbuffer.data(idx_dip_hg + 106);

    auto tr_x_xxyyz_xxxz = pbuffer.data(idx_dip_hg + 107);

    auto tr_x_xxyyz_xxyy = pbuffer.data(idx_dip_hg + 108);

    auto tr_x_xxyyz_xxzz = pbuffer.data(idx_dip_hg + 110);

    auto tr_x_xxyyz_xyyy = pbuffer.data(idx_dip_hg + 111);

    auto tr_x_xxyyz_xzzz = pbuffer.data(idx_dip_hg + 114);

    auto tr_x_xxyyz_yyyy = pbuffer.data(idx_dip_hg + 115);

    auto tr_x_xxyyz_zzzz = pbuffer.data(idx_dip_hg + 119);

    auto tr_x_xxyzz_xxxx = pbuffer.data(idx_dip_hg + 120);

    auto tr_x_xxyzz_xxxz = pbuffer.data(idx_dip_hg + 122);

    auto tr_x_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 124);

    auto tr_x_xxyzz_xxzz = pbuffer.data(idx_dip_hg + 125);

    auto tr_x_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 127);

    auto tr_x_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 128);

    auto tr_x_xxyzz_xzzz = pbuffer.data(idx_dip_hg + 129);

    auto tr_x_xxyzz_zzzz = pbuffer.data(idx_dip_hg + 134);

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

    auto tr_x_xyyyy_xxxx = pbuffer.data(idx_dip_hg + 150);

    auto tr_x_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 151);

    auto tr_x_xyyyy_xxxz = pbuffer.data(idx_dip_hg + 152);

    auto tr_x_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 153);

    auto tr_x_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 154);

    auto tr_x_xyyyy_xxzz = pbuffer.data(idx_dip_hg + 155);

    auto tr_x_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 156);

    auto tr_x_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 157);

    auto tr_x_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 158);

    auto tr_x_xyyyy_xzzz = pbuffer.data(idx_dip_hg + 159);

    auto tr_x_xyyyy_yyyy = pbuffer.data(idx_dip_hg + 160);

    auto tr_x_xyyyy_yyyz = pbuffer.data(idx_dip_hg + 161);

    auto tr_x_xyyyy_yyzz = pbuffer.data(idx_dip_hg + 162);

    auto tr_x_xyyyy_yzzz = pbuffer.data(idx_dip_hg + 163);

    auto tr_x_xyyyz_xxxy = pbuffer.data(idx_dip_hg + 166);

    auto tr_x_xyyyz_xxxz = pbuffer.data(idx_dip_hg + 167);

    auto tr_x_xyyyz_xxyy = pbuffer.data(idx_dip_hg + 168);

    auto tr_x_xyyyz_xxzz = pbuffer.data(idx_dip_hg + 170);

    auto tr_x_xyyyz_xyyy = pbuffer.data(idx_dip_hg + 171);

    auto tr_x_xyyyz_xzzz = pbuffer.data(idx_dip_hg + 174);

    auto tr_x_xyyzz_xxxx = pbuffer.data(idx_dip_hg + 180);

    auto tr_x_xyyzz_xxxy = pbuffer.data(idx_dip_hg + 181);

    auto tr_x_xyyzz_xxxz = pbuffer.data(idx_dip_hg + 182);

    auto tr_x_xyyzz_xxyy = pbuffer.data(idx_dip_hg + 183);

    auto tr_x_xyyzz_xxzz = pbuffer.data(idx_dip_hg + 185);

    auto tr_x_xyyzz_xyyy = pbuffer.data(idx_dip_hg + 186);

    auto tr_x_xyyzz_xzzz = pbuffer.data(idx_dip_hg + 189);

    auto tr_x_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 191);

    auto tr_x_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 192);

    auto tr_x_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 193);

    auto tr_x_xyzzz_xxxx = pbuffer.data(idx_dip_hg + 195);

    auto tr_x_xyzzz_xxxz = pbuffer.data(idx_dip_hg + 197);

    auto tr_x_xyzzz_xxzz = pbuffer.data(idx_dip_hg + 200);

    auto tr_x_xyzzz_xzzz = pbuffer.data(idx_dip_hg + 204);

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

    auto tr_x_xzzzz_yyyz = pbuffer.data(idx_dip_hg + 221);

    auto tr_x_xzzzz_yyzz = pbuffer.data(idx_dip_hg + 222);

    auto tr_x_xzzzz_yzzz = pbuffer.data(idx_dip_hg + 223);

    auto tr_x_xzzzz_zzzz = pbuffer.data(idx_dip_hg + 224);

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

    auto tr_x_yyyyz_xxxy = pbuffer.data(idx_dip_hg + 241);

    auto tr_x_yyyyz_xxxz = pbuffer.data(idx_dip_hg + 242);

    auto tr_x_yyyyz_xxyy = pbuffer.data(idx_dip_hg + 243);

    auto tr_x_yyyyz_xxzz = pbuffer.data(idx_dip_hg + 245);

    auto tr_x_yyyyz_xyyy = pbuffer.data(idx_dip_hg + 246);

    auto tr_x_yyyyz_xzzz = pbuffer.data(idx_dip_hg + 249);

    auto tr_x_yyyyz_yyyy = pbuffer.data(idx_dip_hg + 250);

    auto tr_x_yyyyz_yyyz = pbuffer.data(idx_dip_hg + 251);

    auto tr_x_yyyyz_yyzz = pbuffer.data(idx_dip_hg + 252);

    auto tr_x_yyyyz_yzzz = pbuffer.data(idx_dip_hg + 253);

    auto tr_x_yyyyz_zzzz = pbuffer.data(idx_dip_hg + 254);

    auto tr_x_yyyzz_xxxx = pbuffer.data(idx_dip_hg + 255);

    auto tr_x_yyyzz_xxxy = pbuffer.data(idx_dip_hg + 256);

    auto tr_x_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 257);

    auto tr_x_yyyzz_xxyy = pbuffer.data(idx_dip_hg + 258);

    auto tr_x_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 259);

    auto tr_x_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 260);

    auto tr_x_yyyzz_xyyy = pbuffer.data(idx_dip_hg + 261);

    auto tr_x_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 262);

    auto tr_x_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 263);

    auto tr_x_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 264);

    auto tr_x_yyyzz_yyyy = pbuffer.data(idx_dip_hg + 265);

    auto tr_x_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 266);

    auto tr_x_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 267);

    auto tr_x_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 268);

    auto tr_x_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 269);

    auto tr_x_yyzzz_xxxx = pbuffer.data(idx_dip_hg + 270);

    auto tr_x_yyzzz_xxxy = pbuffer.data(idx_dip_hg + 271);

    auto tr_x_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 272);

    auto tr_x_yyzzz_xxyy = pbuffer.data(idx_dip_hg + 273);

    auto tr_x_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 274);

    auto tr_x_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 275);

    auto tr_x_yyzzz_xyyy = pbuffer.data(idx_dip_hg + 276);

    auto tr_x_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 277);

    auto tr_x_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 278);

    auto tr_x_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 279);

    auto tr_x_yyzzz_yyyy = pbuffer.data(idx_dip_hg + 280);

    auto tr_x_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 281);

    auto tr_x_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 282);

    auto tr_x_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 283);

    auto tr_x_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 284);

    auto tr_x_yzzzz_xxxx = pbuffer.data(idx_dip_hg + 285);

    auto tr_x_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 287);

    auto tr_x_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 289);

    auto tr_x_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 290);

    auto tr_x_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 292);

    auto tr_x_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 293);

    auto tr_x_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 294);

    auto tr_x_yzzzz_yyyy = pbuffer.data(idx_dip_hg + 295);

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

    auto tr_y_xxxxy_xxxx = pbuffer.data(idx_dip_hg + 330);

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

    auto tr_y_xxxxy_zzzz = pbuffer.data(idx_dip_hg + 344);

    auto tr_y_xxxxz_xxxx = pbuffer.data(idx_dip_hg + 345);

    auto tr_y_xxxxz_xxxy = pbuffer.data(idx_dip_hg + 346);

    auto tr_y_xxxxz_xxxz = pbuffer.data(idx_dip_hg + 347);

    auto tr_y_xxxxz_xxyy = pbuffer.data(idx_dip_hg + 348);

    auto tr_y_xxxxz_xxzz = pbuffer.data(idx_dip_hg + 350);

    auto tr_y_xxxxz_xyyy = pbuffer.data(idx_dip_hg + 351);

    auto tr_y_xxxxz_xzzz = pbuffer.data(idx_dip_hg + 354);

    auto tr_y_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 356);

    auto tr_y_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 357);

    auto tr_y_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 358);

    auto tr_y_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 359);

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

    auto tr_y_xxxyz_xxxy = pbuffer.data(idx_dip_hg + 376);

    auto tr_y_xxxyz_xxyy = pbuffer.data(idx_dip_hg + 378);

    auto tr_y_xxxyz_xyyy = pbuffer.data(idx_dip_hg + 381);

    auto tr_y_xxxyz_yyyz = pbuffer.data(idx_dip_hg + 386);

    auto tr_y_xxxyz_yyzz = pbuffer.data(idx_dip_hg + 387);

    auto tr_y_xxxyz_yzzz = pbuffer.data(idx_dip_hg + 388);

    auto tr_y_xxxyz_zzzz = pbuffer.data(idx_dip_hg + 389);

    auto tr_y_xxxzz_xxxx = pbuffer.data(idx_dip_hg + 390);

    auto tr_y_xxxzz_xxxy = pbuffer.data(idx_dip_hg + 391);

    auto tr_y_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 392);

    auto tr_y_xxxzz_xxyy = pbuffer.data(idx_dip_hg + 393);

    auto tr_y_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 394);

    auto tr_y_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 395);

    auto tr_y_xxxzz_xyyy = pbuffer.data(idx_dip_hg + 396);

    auto tr_y_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 397);

    auto tr_y_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 398);

    auto tr_y_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 399);

    auto tr_y_xxxzz_yyyy = pbuffer.data(idx_dip_hg + 400);

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

    auto tr_y_xxyyz_xxxx = pbuffer.data(idx_dip_hg + 420);

    auto tr_y_xxyyz_xxxy = pbuffer.data(idx_dip_hg + 421);

    auto tr_y_xxyyz_xxyy = pbuffer.data(idx_dip_hg + 423);

    auto tr_y_xxyyz_xyyy = pbuffer.data(idx_dip_hg + 426);

    auto tr_y_xxyyz_yyyz = pbuffer.data(idx_dip_hg + 431);

    auto tr_y_xxyyz_yyzz = pbuffer.data(idx_dip_hg + 432);

    auto tr_y_xxyyz_yzzz = pbuffer.data(idx_dip_hg + 433);

    auto tr_y_xxyyz_zzzz = pbuffer.data(idx_dip_hg + 434);

    auto tr_y_xxyzz_xxxy = pbuffer.data(idx_dip_hg + 436);

    auto tr_y_xxyzz_xxyy = pbuffer.data(idx_dip_hg + 438);

    auto tr_y_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 439);

    auto tr_y_xxyzz_xyyy = pbuffer.data(idx_dip_hg + 441);

    auto tr_y_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 442);

    auto tr_y_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 443);

    auto tr_y_xxyzz_yyyy = pbuffer.data(idx_dip_hg + 445);

    auto tr_y_xxyzz_yyyz = pbuffer.data(idx_dip_hg + 446);

    auto tr_y_xxyzz_yyzz = pbuffer.data(idx_dip_hg + 447);

    auto tr_y_xxyzz_yzzz = pbuffer.data(idx_dip_hg + 448);

    auto tr_y_xxyzz_zzzz = pbuffer.data(idx_dip_hg + 449);

    auto tr_y_xxzzz_xxxx = pbuffer.data(idx_dip_hg + 450);

    auto tr_y_xxzzz_xxxy = pbuffer.data(idx_dip_hg + 451);

    auto tr_y_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 452);

    auto tr_y_xxzzz_xxyy = pbuffer.data(idx_dip_hg + 453);

    auto tr_y_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 454);

    auto tr_y_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 455);

    auto tr_y_xxzzz_xyyy = pbuffer.data(idx_dip_hg + 456);

    auto tr_y_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 457);

    auto tr_y_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 458);

    auto tr_y_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 459);

    auto tr_y_xxzzz_yyyy = pbuffer.data(idx_dip_hg + 460);

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

    auto tr_y_xyyyz_yyyz = pbuffer.data(idx_dip_hg + 491);

    auto tr_y_xyyyz_yyzz = pbuffer.data(idx_dip_hg + 492);

    auto tr_y_xyyyz_yzzz = pbuffer.data(idx_dip_hg + 493);

    auto tr_y_xyyyz_zzzz = pbuffer.data(idx_dip_hg + 494);

    auto tr_y_xyyzz_xxxz = pbuffer.data(idx_dip_hg + 497);

    auto tr_y_xyyzz_xxyz = pbuffer.data(idx_dip_hg + 499);

    auto tr_y_xyyzz_xxzz = pbuffer.data(idx_dip_hg + 500);

    auto tr_y_xyyzz_xyyz = pbuffer.data(idx_dip_hg + 502);

    auto tr_y_xyyzz_xyzz = pbuffer.data(idx_dip_hg + 503);

    auto tr_y_xyyzz_xzzz = pbuffer.data(idx_dip_hg + 504);

    auto tr_y_xyyzz_yyyy = pbuffer.data(idx_dip_hg + 505);

    auto tr_y_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 506);

    auto tr_y_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 507);

    auto tr_y_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 508);

    auto tr_y_xyyzz_zzzz = pbuffer.data(idx_dip_hg + 509);

    auto tr_y_xyzzz_xxyz = pbuffer.data(idx_dip_hg + 514);

    auto tr_y_xyzzz_xyyz = pbuffer.data(idx_dip_hg + 517);

    auto tr_y_xyzzz_xyzz = pbuffer.data(idx_dip_hg + 518);

    auto tr_y_xyzzz_yyyy = pbuffer.data(idx_dip_hg + 520);

    auto tr_y_xyzzz_yyyz = pbuffer.data(idx_dip_hg + 521);

    auto tr_y_xyzzz_yyzz = pbuffer.data(idx_dip_hg + 522);

    auto tr_y_xyzzz_yzzz = pbuffer.data(idx_dip_hg + 523);

    auto tr_y_xyzzz_zzzz = pbuffer.data(idx_dip_hg + 524);

    auto tr_y_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 527);

    auto tr_y_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 529);

    auto tr_y_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 530);

    auto tr_y_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 532);

    auto tr_y_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 533);

    auto tr_y_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 534);

    auto tr_y_xzzzz_yyyy = pbuffer.data(idx_dip_hg + 535);

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

    auto tr_y_yyyyz_xxxx = pbuffer.data(idx_dip_hg + 555);

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

    auto tr_z_xxxxy_xxxx = pbuffer.data(idx_dip_hg + 645);

    auto tr_z_xxxxy_xxxz = pbuffer.data(idx_dip_hg + 647);

    auto tr_z_xxxxy_xxzz = pbuffer.data(idx_dip_hg + 650);

    auto tr_z_xxxxy_xzzz = pbuffer.data(idx_dip_hg + 654);

    auto tr_z_xxxxy_yyyy = pbuffer.data(idx_dip_hg + 655);

    auto tr_z_xxxxy_yyyz = pbuffer.data(idx_dip_hg + 656);

    auto tr_z_xxxxy_yyzz = pbuffer.data(idx_dip_hg + 657);

    auto tr_z_xxxxy_yzzz = pbuffer.data(idx_dip_hg + 658);

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

    auto tr_z_xxxxz_yyyy = pbuffer.data(idx_dip_hg + 670);

    auto tr_z_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 671);

    auto tr_z_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 672);

    auto tr_z_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 673);

    auto tr_z_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 674);

    auto tr_z_xxxyy_xxxx = pbuffer.data(idx_dip_hg + 675);

    auto tr_z_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 676);

    auto tr_z_xxxyy_xxxz = pbuffer.data(idx_dip_hg + 677);

    auto tr_z_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 678);

    auto tr_z_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 679);

    auto tr_z_xxxyy_xxzz = pbuffer.data(idx_dip_hg + 680);

    auto tr_z_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 681);

    auto tr_z_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 682);

    auto tr_z_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 683);

    auto tr_z_xxxyy_xzzz = pbuffer.data(idx_dip_hg + 684);

    auto tr_z_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 685);

    auto tr_z_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 686);

    auto tr_z_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 687);

    auto tr_z_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 688);

    auto tr_z_xxxyy_zzzz = pbuffer.data(idx_dip_hg + 689);

    auto tr_z_xxxyz_xxxx = pbuffer.data(idx_dip_hg + 690);

    auto tr_z_xxxyz_xxxz = pbuffer.data(idx_dip_hg + 692);

    auto tr_z_xxxyz_xxzz = pbuffer.data(idx_dip_hg + 695);

    auto tr_z_xxxyz_xzzz = pbuffer.data(idx_dip_hg + 699);

    auto tr_z_xxxyz_yyyy = pbuffer.data(idx_dip_hg + 700);

    auto tr_z_xxxyz_yyyz = pbuffer.data(idx_dip_hg + 701);

    auto tr_z_xxxyz_yyzz = pbuffer.data(idx_dip_hg + 702);

    auto tr_z_xxxyz_yzzz = pbuffer.data(idx_dip_hg + 703);

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

    auto tr_z_xxyyy_xxxx = pbuffer.data(idx_dip_hg + 720);

    auto tr_z_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 721);

    auto tr_z_xxyyy_xxxz = pbuffer.data(idx_dip_hg + 722);

    auto tr_z_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 723);

    auto tr_z_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 724);

    auto tr_z_xxyyy_xxzz = pbuffer.data(idx_dip_hg + 725);

    auto tr_z_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 726);

    auto tr_z_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 727);

    auto tr_z_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 728);

    auto tr_z_xxyyy_xzzz = pbuffer.data(idx_dip_hg + 729);

    auto tr_z_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 730);

    auto tr_z_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 731);

    auto tr_z_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 732);

    auto tr_z_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 733);

    auto tr_z_xxyyy_zzzz = pbuffer.data(idx_dip_hg + 734);

    auto tr_z_xxyyz_xxxx = pbuffer.data(idx_dip_hg + 735);

    auto tr_z_xxyyz_xxxz = pbuffer.data(idx_dip_hg + 737);

    auto tr_z_xxyyz_xxyz = pbuffer.data(idx_dip_hg + 739);

    auto tr_z_xxyyz_xxzz = pbuffer.data(idx_dip_hg + 740);

    auto tr_z_xxyyz_xyyz = pbuffer.data(idx_dip_hg + 742);

    auto tr_z_xxyyz_xyzz = pbuffer.data(idx_dip_hg + 743);

    auto tr_z_xxyyz_xzzz = pbuffer.data(idx_dip_hg + 744);

    auto tr_z_xxyyz_yyyy = pbuffer.data(idx_dip_hg + 745);

    auto tr_z_xxyyz_yyyz = pbuffer.data(idx_dip_hg + 746);

    auto tr_z_xxyyz_yyzz = pbuffer.data(idx_dip_hg + 747);

    auto tr_z_xxyyz_yzzz = pbuffer.data(idx_dip_hg + 748);

    auto tr_z_xxyyz_zzzz = pbuffer.data(idx_dip_hg + 749);

    auto tr_z_xxyzz_xxxx = pbuffer.data(idx_dip_hg + 750);

    auto tr_z_xxyzz_xxxz = pbuffer.data(idx_dip_hg + 752);

    auto tr_z_xxyzz_xxzz = pbuffer.data(idx_dip_hg + 755);

    auto tr_z_xxyzz_xzzz = pbuffer.data(idx_dip_hg + 759);

    auto tr_z_xxyzz_yyyy = pbuffer.data(idx_dip_hg + 760);

    auto tr_z_xxyzz_yyyz = pbuffer.data(idx_dip_hg + 761);

    auto tr_z_xxyzz_yyzz = pbuffer.data(idx_dip_hg + 762);

    auto tr_z_xxyzz_yzzz = pbuffer.data(idx_dip_hg + 763);

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

    auto tr_z_xyyyy_zzzz = pbuffer.data(idx_dip_hg + 794);

    auto tr_z_xyyyz_xxyz = pbuffer.data(idx_dip_hg + 799);

    auto tr_z_xyyyz_xyyz = pbuffer.data(idx_dip_hg + 802);

    auto tr_z_xyyyz_xyzz = pbuffer.data(idx_dip_hg + 803);

    auto tr_z_xyyyz_yyyy = pbuffer.data(idx_dip_hg + 805);

    auto tr_z_xyyyz_yyyz = pbuffer.data(idx_dip_hg + 806);

    auto tr_z_xyyyz_yyzz = pbuffer.data(idx_dip_hg + 807);

    auto tr_z_xyyyz_yzzz = pbuffer.data(idx_dip_hg + 808);

    auto tr_z_xyyyz_zzzz = pbuffer.data(idx_dip_hg + 809);

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

    auto tr_z_xyyzz_zzzz = pbuffer.data(idx_dip_hg + 824);

    auto tr_z_xyzzz_yyyy = pbuffer.data(idx_dip_hg + 835);

    auto tr_z_xyzzz_yyyz = pbuffer.data(idx_dip_hg + 836);

    auto tr_z_xyzzz_yyzz = pbuffer.data(idx_dip_hg + 837);

    auto tr_z_xyzzz_yzzz = pbuffer.data(idx_dip_hg + 838);

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

    // Set up 0-15 components of targeted buffer : IG

    auto tr_x_xxxxxx_xxxx = pbuffer.data(idx_dip_ig);

    auto tr_x_xxxxxx_xxxy = pbuffer.data(idx_dip_ig + 1);

    auto tr_x_xxxxxx_xxxz = pbuffer.data(idx_dip_ig + 2);

    auto tr_x_xxxxxx_xxyy = pbuffer.data(idx_dip_ig + 3);

    auto tr_x_xxxxxx_xxyz = pbuffer.data(idx_dip_ig + 4);

    auto tr_x_xxxxxx_xxzz = pbuffer.data(idx_dip_ig + 5);

    auto tr_x_xxxxxx_xyyy = pbuffer.data(idx_dip_ig + 6);

    auto tr_x_xxxxxx_xyyz = pbuffer.data(idx_dip_ig + 7);

    auto tr_x_xxxxxx_xyzz = pbuffer.data(idx_dip_ig + 8);

    auto tr_x_xxxxxx_xzzz = pbuffer.data(idx_dip_ig + 9);

    auto tr_x_xxxxxx_yyyy = pbuffer.data(idx_dip_ig + 10);

    auto tr_x_xxxxxx_yyyz = pbuffer.data(idx_dip_ig + 11);

    auto tr_x_xxxxxx_yyzz = pbuffer.data(idx_dip_ig + 12);

    auto tr_x_xxxxxx_yzzz = pbuffer.data(idx_dip_ig + 13);

    auto tr_x_xxxxxx_zzzz = pbuffer.data(idx_dip_ig + 14);

#pragma omp simd aligned(pa_x,                 \
                             tr_x_xxxx_xxxx,   \
                             tr_x_xxxx_xxxy,   \
                             tr_x_xxxx_xxxz,   \
                             tr_x_xxxx_xxyy,   \
                             tr_x_xxxx_xxyz,   \
                             tr_x_xxxx_xxzz,   \
                             tr_x_xxxx_xyyy,   \
                             tr_x_xxxx_xyyz,   \
                             tr_x_xxxx_xyzz,   \
                             tr_x_xxxx_xzzz,   \
                             tr_x_xxxx_yyyy,   \
                             tr_x_xxxx_yyyz,   \
                             tr_x_xxxx_yyzz,   \
                             tr_x_xxxx_yzzz,   \
                             tr_x_xxxx_zzzz,   \
                             tr_x_xxxxx_xxx,   \
                             tr_x_xxxxx_xxxx,  \
                             tr_x_xxxxx_xxxy,  \
                             tr_x_xxxxx_xxxz,  \
                             tr_x_xxxxx_xxy,   \
                             tr_x_xxxxx_xxyy,  \
                             tr_x_xxxxx_xxyz,  \
                             tr_x_xxxxx_xxz,   \
                             tr_x_xxxxx_xxzz,  \
                             tr_x_xxxxx_xyy,   \
                             tr_x_xxxxx_xyyy,  \
                             tr_x_xxxxx_xyyz,  \
                             tr_x_xxxxx_xyz,   \
                             tr_x_xxxxx_xyzz,  \
                             tr_x_xxxxx_xzz,   \
                             tr_x_xxxxx_xzzz,  \
                             tr_x_xxxxx_yyy,   \
                             tr_x_xxxxx_yyyy,  \
                             tr_x_xxxxx_yyyz,  \
                             tr_x_xxxxx_yyz,   \
                             tr_x_xxxxx_yyzz,  \
                             tr_x_xxxxx_yzz,   \
                             tr_x_xxxxx_yzzz,  \
                             tr_x_xxxxx_zzz,   \
                             tr_x_xxxxx_zzzz,  \
                             tr_x_xxxxxx_xxxx, \
                             tr_x_xxxxxx_xxxy, \
                             tr_x_xxxxxx_xxxz, \
                             tr_x_xxxxxx_xxyy, \
                             tr_x_xxxxxx_xxyz, \
                             tr_x_xxxxxx_xxzz, \
                             tr_x_xxxxxx_xyyy, \
                             tr_x_xxxxxx_xyyz, \
                             tr_x_xxxxxx_xyzz, \
                             tr_x_xxxxxx_xzzz, \
                             tr_x_xxxxxx_yyyy, \
                             tr_x_xxxxxx_yyyz, \
                             tr_x_xxxxxx_yyzz, \
                             tr_x_xxxxxx_yzzz, \
                             tr_x_xxxxxx_zzzz, \
                             ts_xxxxx_xxxx,    \
                             ts_xxxxx_xxxy,    \
                             ts_xxxxx_xxxz,    \
                             ts_xxxxx_xxyy,    \
                             ts_xxxxx_xxyz,    \
                             ts_xxxxx_xxzz,    \
                             ts_xxxxx_xyyy,    \
                             ts_xxxxx_xyyz,    \
                             ts_xxxxx_xyzz,    \
                             ts_xxxxx_xzzz,    \
                             ts_xxxxx_yyyy,    \
                             ts_xxxxx_yyyz,    \
                             ts_xxxxx_yyzz,    \
                             ts_xxxxx_yzzz,    \
                             ts_xxxxx_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_xxxx[i] =
            5.0 * tr_x_xxxx_xxxx[i] * fe_0 + 4.0 * tr_x_xxxxx_xxx[i] * fe_0 + ts_xxxxx_xxxx[i] * fe_0 + tr_x_xxxxx_xxxx[i] * pa_x[i];

        tr_x_xxxxxx_xxxy[i] =
            5.0 * tr_x_xxxx_xxxy[i] * fe_0 + 3.0 * tr_x_xxxxx_xxy[i] * fe_0 + ts_xxxxx_xxxy[i] * fe_0 + tr_x_xxxxx_xxxy[i] * pa_x[i];

        tr_x_xxxxxx_xxxz[i] =
            5.0 * tr_x_xxxx_xxxz[i] * fe_0 + 3.0 * tr_x_xxxxx_xxz[i] * fe_0 + ts_xxxxx_xxxz[i] * fe_0 + tr_x_xxxxx_xxxz[i] * pa_x[i];

        tr_x_xxxxxx_xxyy[i] =
            5.0 * tr_x_xxxx_xxyy[i] * fe_0 + 2.0 * tr_x_xxxxx_xyy[i] * fe_0 + ts_xxxxx_xxyy[i] * fe_0 + tr_x_xxxxx_xxyy[i] * pa_x[i];

        tr_x_xxxxxx_xxyz[i] =
            5.0 * tr_x_xxxx_xxyz[i] * fe_0 + 2.0 * tr_x_xxxxx_xyz[i] * fe_0 + ts_xxxxx_xxyz[i] * fe_0 + tr_x_xxxxx_xxyz[i] * pa_x[i];

        tr_x_xxxxxx_xxzz[i] =
            5.0 * tr_x_xxxx_xxzz[i] * fe_0 + 2.0 * tr_x_xxxxx_xzz[i] * fe_0 + ts_xxxxx_xxzz[i] * fe_0 + tr_x_xxxxx_xxzz[i] * pa_x[i];

        tr_x_xxxxxx_xyyy[i] = 5.0 * tr_x_xxxx_xyyy[i] * fe_0 + tr_x_xxxxx_yyy[i] * fe_0 + ts_xxxxx_xyyy[i] * fe_0 + tr_x_xxxxx_xyyy[i] * pa_x[i];

        tr_x_xxxxxx_xyyz[i] = 5.0 * tr_x_xxxx_xyyz[i] * fe_0 + tr_x_xxxxx_yyz[i] * fe_0 + ts_xxxxx_xyyz[i] * fe_0 + tr_x_xxxxx_xyyz[i] * pa_x[i];

        tr_x_xxxxxx_xyzz[i] = 5.0 * tr_x_xxxx_xyzz[i] * fe_0 + tr_x_xxxxx_yzz[i] * fe_0 + ts_xxxxx_xyzz[i] * fe_0 + tr_x_xxxxx_xyzz[i] * pa_x[i];

        tr_x_xxxxxx_xzzz[i] = 5.0 * tr_x_xxxx_xzzz[i] * fe_0 + tr_x_xxxxx_zzz[i] * fe_0 + ts_xxxxx_xzzz[i] * fe_0 + tr_x_xxxxx_xzzz[i] * pa_x[i];

        tr_x_xxxxxx_yyyy[i] = 5.0 * tr_x_xxxx_yyyy[i] * fe_0 + ts_xxxxx_yyyy[i] * fe_0 + tr_x_xxxxx_yyyy[i] * pa_x[i];

        tr_x_xxxxxx_yyyz[i] = 5.0 * tr_x_xxxx_yyyz[i] * fe_0 + ts_xxxxx_yyyz[i] * fe_0 + tr_x_xxxxx_yyyz[i] * pa_x[i];

        tr_x_xxxxxx_yyzz[i] = 5.0 * tr_x_xxxx_yyzz[i] * fe_0 + ts_xxxxx_yyzz[i] * fe_0 + tr_x_xxxxx_yyzz[i] * pa_x[i];

        tr_x_xxxxxx_yzzz[i] = 5.0 * tr_x_xxxx_yzzz[i] * fe_0 + ts_xxxxx_yzzz[i] * fe_0 + tr_x_xxxxx_yzzz[i] * pa_x[i];

        tr_x_xxxxxx_zzzz[i] = 5.0 * tr_x_xxxx_zzzz[i] * fe_0 + ts_xxxxx_zzzz[i] * fe_0 + tr_x_xxxxx_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : IG

    auto tr_x_xxxxxy_xxxx = pbuffer.data(idx_dip_ig + 15);

    auto tr_x_xxxxxy_xxxy = pbuffer.data(idx_dip_ig + 16);

    auto tr_x_xxxxxy_xxxz = pbuffer.data(idx_dip_ig + 17);

    auto tr_x_xxxxxy_xxyy = pbuffer.data(idx_dip_ig + 18);

    auto tr_x_xxxxxy_xxyz = pbuffer.data(idx_dip_ig + 19);

    auto tr_x_xxxxxy_xxzz = pbuffer.data(idx_dip_ig + 20);

    auto tr_x_xxxxxy_xyyy = pbuffer.data(idx_dip_ig + 21);

    auto tr_x_xxxxxy_xyyz = pbuffer.data(idx_dip_ig + 22);

    auto tr_x_xxxxxy_xyzz = pbuffer.data(idx_dip_ig + 23);

    auto tr_x_xxxxxy_xzzz = pbuffer.data(idx_dip_ig + 24);

    auto tr_x_xxxxxy_yyyy = pbuffer.data(idx_dip_ig + 25);

    auto tr_x_xxxxxy_yyyz = pbuffer.data(idx_dip_ig + 26);

    auto tr_x_xxxxxy_yyzz = pbuffer.data(idx_dip_ig + 27);

    auto tr_x_xxxxxy_yzzz = pbuffer.data(idx_dip_ig + 28);

    auto tr_x_xxxxxy_zzzz = pbuffer.data(idx_dip_ig + 29);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_xxxxx_xxx,   \
                             tr_x_xxxxx_xxxx,  \
                             tr_x_xxxxx_xxxy,  \
                             tr_x_xxxxx_xxxz,  \
                             tr_x_xxxxx_xxy,   \
                             tr_x_xxxxx_xxyy,  \
                             tr_x_xxxxx_xxyz,  \
                             tr_x_xxxxx_xxz,   \
                             tr_x_xxxxx_xxzz,  \
                             tr_x_xxxxx_xyy,   \
                             tr_x_xxxxx_xyyy,  \
                             tr_x_xxxxx_xyyz,  \
                             tr_x_xxxxx_xyz,   \
                             tr_x_xxxxx_xyzz,  \
                             tr_x_xxxxx_xzz,   \
                             tr_x_xxxxx_xzzz,  \
                             tr_x_xxxxx_yyy,   \
                             tr_x_xxxxx_yyyy,  \
                             tr_x_xxxxx_yyyz,  \
                             tr_x_xxxxx_yyz,   \
                             tr_x_xxxxx_yyzz,  \
                             tr_x_xxxxx_yzz,   \
                             tr_x_xxxxx_yzzz,  \
                             tr_x_xxxxx_zzz,   \
                             tr_x_xxxxx_zzzz,  \
                             tr_x_xxxxxy_xxxx, \
                             tr_x_xxxxxy_xxxy, \
                             tr_x_xxxxxy_xxxz, \
                             tr_x_xxxxxy_xxyy, \
                             tr_x_xxxxxy_xxyz, \
                             tr_x_xxxxxy_xxzz, \
                             tr_x_xxxxxy_xyyy, \
                             tr_x_xxxxxy_xyyz, \
                             tr_x_xxxxxy_xyzz, \
                             tr_x_xxxxxy_xzzz, \
                             tr_x_xxxxxy_yyyy, \
                             tr_x_xxxxxy_yyyz, \
                             tr_x_xxxxxy_yyzz, \
                             tr_x_xxxxxy_yzzz, \
                             tr_x_xxxxxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxy_xxxx[i] = tr_x_xxxxx_xxxx[i] * pa_y[i];

        tr_x_xxxxxy_xxxy[i] = tr_x_xxxxx_xxx[i] * fe_0 + tr_x_xxxxx_xxxy[i] * pa_y[i];

        tr_x_xxxxxy_xxxz[i] = tr_x_xxxxx_xxxz[i] * pa_y[i];

        tr_x_xxxxxy_xxyy[i] = 2.0 * tr_x_xxxxx_xxy[i] * fe_0 + tr_x_xxxxx_xxyy[i] * pa_y[i];

        tr_x_xxxxxy_xxyz[i] = tr_x_xxxxx_xxz[i] * fe_0 + tr_x_xxxxx_xxyz[i] * pa_y[i];

        tr_x_xxxxxy_xxzz[i] = tr_x_xxxxx_xxzz[i] * pa_y[i];

        tr_x_xxxxxy_xyyy[i] = 3.0 * tr_x_xxxxx_xyy[i] * fe_0 + tr_x_xxxxx_xyyy[i] * pa_y[i];

        tr_x_xxxxxy_xyyz[i] = 2.0 * tr_x_xxxxx_xyz[i] * fe_0 + tr_x_xxxxx_xyyz[i] * pa_y[i];

        tr_x_xxxxxy_xyzz[i] = tr_x_xxxxx_xzz[i] * fe_0 + tr_x_xxxxx_xyzz[i] * pa_y[i];

        tr_x_xxxxxy_xzzz[i] = tr_x_xxxxx_xzzz[i] * pa_y[i];

        tr_x_xxxxxy_yyyy[i] = 4.0 * tr_x_xxxxx_yyy[i] * fe_0 + tr_x_xxxxx_yyyy[i] * pa_y[i];

        tr_x_xxxxxy_yyyz[i] = 3.0 * tr_x_xxxxx_yyz[i] * fe_0 + tr_x_xxxxx_yyyz[i] * pa_y[i];

        tr_x_xxxxxy_yyzz[i] = 2.0 * tr_x_xxxxx_yzz[i] * fe_0 + tr_x_xxxxx_yyzz[i] * pa_y[i];

        tr_x_xxxxxy_yzzz[i] = tr_x_xxxxx_zzz[i] * fe_0 + tr_x_xxxxx_yzzz[i] * pa_y[i];

        tr_x_xxxxxy_zzzz[i] = tr_x_xxxxx_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : IG

    auto tr_x_xxxxxz_xxxx = pbuffer.data(idx_dip_ig + 30);

    auto tr_x_xxxxxz_xxxy = pbuffer.data(idx_dip_ig + 31);

    auto tr_x_xxxxxz_xxxz = pbuffer.data(idx_dip_ig + 32);

    auto tr_x_xxxxxz_xxyy = pbuffer.data(idx_dip_ig + 33);

    auto tr_x_xxxxxz_xxyz = pbuffer.data(idx_dip_ig + 34);

    auto tr_x_xxxxxz_xxzz = pbuffer.data(idx_dip_ig + 35);

    auto tr_x_xxxxxz_xyyy = pbuffer.data(idx_dip_ig + 36);

    auto tr_x_xxxxxz_xyyz = pbuffer.data(idx_dip_ig + 37);

    auto tr_x_xxxxxz_xyzz = pbuffer.data(idx_dip_ig + 38);

    auto tr_x_xxxxxz_xzzz = pbuffer.data(idx_dip_ig + 39);

    auto tr_x_xxxxxz_yyyy = pbuffer.data(idx_dip_ig + 40);

    auto tr_x_xxxxxz_yyyz = pbuffer.data(idx_dip_ig + 41);

    auto tr_x_xxxxxz_yyzz = pbuffer.data(idx_dip_ig + 42);

    auto tr_x_xxxxxz_yzzz = pbuffer.data(idx_dip_ig + 43);

    auto tr_x_xxxxxz_zzzz = pbuffer.data(idx_dip_ig + 44);

#pragma omp simd aligned(pa_z,                 \
                             tr_x_xxxxx_xxx,   \
                             tr_x_xxxxx_xxxx,  \
                             tr_x_xxxxx_xxxy,  \
                             tr_x_xxxxx_xxxz,  \
                             tr_x_xxxxx_xxy,   \
                             tr_x_xxxxx_xxyy,  \
                             tr_x_xxxxx_xxyz,  \
                             tr_x_xxxxx_xxz,   \
                             tr_x_xxxxx_xxzz,  \
                             tr_x_xxxxx_xyy,   \
                             tr_x_xxxxx_xyyy,  \
                             tr_x_xxxxx_xyyz,  \
                             tr_x_xxxxx_xyz,   \
                             tr_x_xxxxx_xyzz,  \
                             tr_x_xxxxx_xzz,   \
                             tr_x_xxxxx_xzzz,  \
                             tr_x_xxxxx_yyy,   \
                             tr_x_xxxxx_yyyy,  \
                             tr_x_xxxxx_yyyz,  \
                             tr_x_xxxxx_yyz,   \
                             tr_x_xxxxx_yyzz,  \
                             tr_x_xxxxx_yzz,   \
                             tr_x_xxxxx_yzzz,  \
                             tr_x_xxxxx_zzz,   \
                             tr_x_xxxxx_zzzz,  \
                             tr_x_xxxxxz_xxxx, \
                             tr_x_xxxxxz_xxxy, \
                             tr_x_xxxxxz_xxxz, \
                             tr_x_xxxxxz_xxyy, \
                             tr_x_xxxxxz_xxyz, \
                             tr_x_xxxxxz_xxzz, \
                             tr_x_xxxxxz_xyyy, \
                             tr_x_xxxxxz_xyyz, \
                             tr_x_xxxxxz_xyzz, \
                             tr_x_xxxxxz_xzzz, \
                             tr_x_xxxxxz_yyyy, \
                             tr_x_xxxxxz_yyyz, \
                             tr_x_xxxxxz_yyzz, \
                             tr_x_xxxxxz_yzzz, \
                             tr_x_xxxxxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxz_xxxx[i] = tr_x_xxxxx_xxxx[i] * pa_z[i];

        tr_x_xxxxxz_xxxy[i] = tr_x_xxxxx_xxxy[i] * pa_z[i];

        tr_x_xxxxxz_xxxz[i] = tr_x_xxxxx_xxx[i] * fe_0 + tr_x_xxxxx_xxxz[i] * pa_z[i];

        tr_x_xxxxxz_xxyy[i] = tr_x_xxxxx_xxyy[i] * pa_z[i];

        tr_x_xxxxxz_xxyz[i] = tr_x_xxxxx_xxy[i] * fe_0 + tr_x_xxxxx_xxyz[i] * pa_z[i];

        tr_x_xxxxxz_xxzz[i] = 2.0 * tr_x_xxxxx_xxz[i] * fe_0 + tr_x_xxxxx_xxzz[i] * pa_z[i];

        tr_x_xxxxxz_xyyy[i] = tr_x_xxxxx_xyyy[i] * pa_z[i];

        tr_x_xxxxxz_xyyz[i] = tr_x_xxxxx_xyy[i] * fe_0 + tr_x_xxxxx_xyyz[i] * pa_z[i];

        tr_x_xxxxxz_xyzz[i] = 2.0 * tr_x_xxxxx_xyz[i] * fe_0 + tr_x_xxxxx_xyzz[i] * pa_z[i];

        tr_x_xxxxxz_xzzz[i] = 3.0 * tr_x_xxxxx_xzz[i] * fe_0 + tr_x_xxxxx_xzzz[i] * pa_z[i];

        tr_x_xxxxxz_yyyy[i] = tr_x_xxxxx_yyyy[i] * pa_z[i];

        tr_x_xxxxxz_yyyz[i] = tr_x_xxxxx_yyy[i] * fe_0 + tr_x_xxxxx_yyyz[i] * pa_z[i];

        tr_x_xxxxxz_yyzz[i] = 2.0 * tr_x_xxxxx_yyz[i] * fe_0 + tr_x_xxxxx_yyzz[i] * pa_z[i];

        tr_x_xxxxxz_yzzz[i] = 3.0 * tr_x_xxxxx_yzz[i] * fe_0 + tr_x_xxxxx_yzzz[i] * pa_z[i];

        tr_x_xxxxxz_zzzz[i] = 4.0 * tr_x_xxxxx_zzz[i] * fe_0 + tr_x_xxxxx_zzzz[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : IG

    auto tr_x_xxxxyy_xxxx = pbuffer.data(idx_dip_ig + 45);

    auto tr_x_xxxxyy_xxxy = pbuffer.data(idx_dip_ig + 46);

    auto tr_x_xxxxyy_xxxz = pbuffer.data(idx_dip_ig + 47);

    auto tr_x_xxxxyy_xxyy = pbuffer.data(idx_dip_ig + 48);

    auto tr_x_xxxxyy_xxyz = pbuffer.data(idx_dip_ig + 49);

    auto tr_x_xxxxyy_xxzz = pbuffer.data(idx_dip_ig + 50);

    auto tr_x_xxxxyy_xyyy = pbuffer.data(idx_dip_ig + 51);

    auto tr_x_xxxxyy_xyyz = pbuffer.data(idx_dip_ig + 52);

    auto tr_x_xxxxyy_xyzz = pbuffer.data(idx_dip_ig + 53);

    auto tr_x_xxxxyy_xzzz = pbuffer.data(idx_dip_ig + 54);

    auto tr_x_xxxxyy_yyyy = pbuffer.data(idx_dip_ig + 55);

    auto tr_x_xxxxyy_yyyz = pbuffer.data(idx_dip_ig + 56);

    auto tr_x_xxxxyy_yyzz = pbuffer.data(idx_dip_ig + 57);

    auto tr_x_xxxxyy_yzzz = pbuffer.data(idx_dip_ig + 58);

    auto tr_x_xxxxyy_zzzz = pbuffer.data(idx_dip_ig + 59);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xxxx_xxxx,   \
                             tr_x_xxxx_xxxy,   \
                             tr_x_xxxx_xxxz,   \
                             tr_x_xxxx_xxyy,   \
                             tr_x_xxxx_xxyz,   \
                             tr_x_xxxx_xxzz,   \
                             tr_x_xxxx_xyyy,   \
                             tr_x_xxxx_xyyz,   \
                             tr_x_xxxx_xyzz,   \
                             tr_x_xxxx_xzzz,   \
                             tr_x_xxxx_zzzz,   \
                             tr_x_xxxxy_xxx,   \
                             tr_x_xxxxy_xxxx,  \
                             tr_x_xxxxy_xxxy,  \
                             tr_x_xxxxy_xxxz,  \
                             tr_x_xxxxy_xxy,   \
                             tr_x_xxxxy_xxyy,  \
                             tr_x_xxxxy_xxyz,  \
                             tr_x_xxxxy_xxz,   \
                             tr_x_xxxxy_xxzz,  \
                             tr_x_xxxxy_xyy,   \
                             tr_x_xxxxy_xyyy,  \
                             tr_x_xxxxy_xyyz,  \
                             tr_x_xxxxy_xyz,   \
                             tr_x_xxxxy_xyzz,  \
                             tr_x_xxxxy_xzz,   \
                             tr_x_xxxxy_xzzz,  \
                             tr_x_xxxxy_zzzz,  \
                             tr_x_xxxxyy_xxxx, \
                             tr_x_xxxxyy_xxxy, \
                             tr_x_xxxxyy_xxxz, \
                             tr_x_xxxxyy_xxyy, \
                             tr_x_xxxxyy_xxyz, \
                             tr_x_xxxxyy_xxzz, \
                             tr_x_xxxxyy_xyyy, \
                             tr_x_xxxxyy_xyyz, \
                             tr_x_xxxxyy_xyzz, \
                             tr_x_xxxxyy_xzzz, \
                             tr_x_xxxxyy_yyyy, \
                             tr_x_xxxxyy_yyyz, \
                             tr_x_xxxxyy_yyzz, \
                             tr_x_xxxxyy_yzzz, \
                             tr_x_xxxxyy_zzzz, \
                             tr_x_xxxyy_yyyy,  \
                             tr_x_xxxyy_yyyz,  \
                             tr_x_xxxyy_yyzz,  \
                             tr_x_xxxyy_yzzz,  \
                             tr_x_xxyy_yyyy,   \
                             tr_x_xxyy_yyyz,   \
                             tr_x_xxyy_yyzz,   \
                             tr_x_xxyy_yzzz,   \
                             ts_xxxyy_yyyy,    \
                             ts_xxxyy_yyyz,    \
                             ts_xxxyy_yyzz,    \
                             ts_xxxyy_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyy_xxxx[i] = tr_x_xxxx_xxxx[i] * fe_0 + tr_x_xxxxy_xxxx[i] * pa_y[i];

        tr_x_xxxxyy_xxxy[i] = tr_x_xxxx_xxxy[i] * fe_0 + tr_x_xxxxy_xxx[i] * fe_0 + tr_x_xxxxy_xxxy[i] * pa_y[i];

        tr_x_xxxxyy_xxxz[i] = tr_x_xxxx_xxxz[i] * fe_0 + tr_x_xxxxy_xxxz[i] * pa_y[i];

        tr_x_xxxxyy_xxyy[i] = tr_x_xxxx_xxyy[i] * fe_0 + 2.0 * tr_x_xxxxy_xxy[i] * fe_0 + tr_x_xxxxy_xxyy[i] * pa_y[i];

        tr_x_xxxxyy_xxyz[i] = tr_x_xxxx_xxyz[i] * fe_0 + tr_x_xxxxy_xxz[i] * fe_0 + tr_x_xxxxy_xxyz[i] * pa_y[i];

        tr_x_xxxxyy_xxzz[i] = tr_x_xxxx_xxzz[i] * fe_0 + tr_x_xxxxy_xxzz[i] * pa_y[i];

        tr_x_xxxxyy_xyyy[i] = tr_x_xxxx_xyyy[i] * fe_0 + 3.0 * tr_x_xxxxy_xyy[i] * fe_0 + tr_x_xxxxy_xyyy[i] * pa_y[i];

        tr_x_xxxxyy_xyyz[i] = tr_x_xxxx_xyyz[i] * fe_0 + 2.0 * tr_x_xxxxy_xyz[i] * fe_0 + tr_x_xxxxy_xyyz[i] * pa_y[i];

        tr_x_xxxxyy_xyzz[i] = tr_x_xxxx_xyzz[i] * fe_0 + tr_x_xxxxy_xzz[i] * fe_0 + tr_x_xxxxy_xyzz[i] * pa_y[i];

        tr_x_xxxxyy_xzzz[i] = tr_x_xxxx_xzzz[i] * fe_0 + tr_x_xxxxy_xzzz[i] * pa_y[i];

        tr_x_xxxxyy_yyyy[i] = 3.0 * tr_x_xxyy_yyyy[i] * fe_0 + ts_xxxyy_yyyy[i] * fe_0 + tr_x_xxxyy_yyyy[i] * pa_x[i];

        tr_x_xxxxyy_yyyz[i] = 3.0 * tr_x_xxyy_yyyz[i] * fe_0 + ts_xxxyy_yyyz[i] * fe_0 + tr_x_xxxyy_yyyz[i] * pa_x[i];

        tr_x_xxxxyy_yyzz[i] = 3.0 * tr_x_xxyy_yyzz[i] * fe_0 + ts_xxxyy_yyzz[i] * fe_0 + tr_x_xxxyy_yyzz[i] * pa_x[i];

        tr_x_xxxxyy_yzzz[i] = 3.0 * tr_x_xxyy_yzzz[i] * fe_0 + ts_xxxyy_yzzz[i] * fe_0 + tr_x_xxxyy_yzzz[i] * pa_x[i];

        tr_x_xxxxyy_zzzz[i] = tr_x_xxxx_zzzz[i] * fe_0 + tr_x_xxxxy_zzzz[i] * pa_y[i];
    }

    // Set up 60-75 components of targeted buffer : IG

    auto tr_x_xxxxyz_xxxx = pbuffer.data(idx_dip_ig + 60);

    auto tr_x_xxxxyz_xxxy = pbuffer.data(idx_dip_ig + 61);

    auto tr_x_xxxxyz_xxxz = pbuffer.data(idx_dip_ig + 62);

    auto tr_x_xxxxyz_xxyy = pbuffer.data(idx_dip_ig + 63);

    auto tr_x_xxxxyz_xxyz = pbuffer.data(idx_dip_ig + 64);

    auto tr_x_xxxxyz_xxzz = pbuffer.data(idx_dip_ig + 65);

    auto tr_x_xxxxyz_xyyy = pbuffer.data(idx_dip_ig + 66);

    auto tr_x_xxxxyz_xyyz = pbuffer.data(idx_dip_ig + 67);

    auto tr_x_xxxxyz_xyzz = pbuffer.data(idx_dip_ig + 68);

    auto tr_x_xxxxyz_xzzz = pbuffer.data(idx_dip_ig + 69);

    auto tr_x_xxxxyz_yyyy = pbuffer.data(idx_dip_ig + 70);

    auto tr_x_xxxxyz_yyyz = pbuffer.data(idx_dip_ig + 71);

    auto tr_x_xxxxyz_yyzz = pbuffer.data(idx_dip_ig + 72);

    auto tr_x_xxxxyz_yzzz = pbuffer.data(idx_dip_ig + 73);

    auto tr_x_xxxxyz_zzzz = pbuffer.data(idx_dip_ig + 74);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_xxxxy_xxxy,  \
                             tr_x_xxxxy_xxyy,  \
                             tr_x_xxxxy_xyyy,  \
                             tr_x_xxxxy_yyyy,  \
                             tr_x_xxxxyz_xxxx, \
                             tr_x_xxxxyz_xxxy, \
                             tr_x_xxxxyz_xxxz, \
                             tr_x_xxxxyz_xxyy, \
                             tr_x_xxxxyz_xxyz, \
                             tr_x_xxxxyz_xxzz, \
                             tr_x_xxxxyz_xyyy, \
                             tr_x_xxxxyz_xyyz, \
                             tr_x_xxxxyz_xyzz, \
                             tr_x_xxxxyz_xzzz, \
                             tr_x_xxxxyz_yyyy, \
                             tr_x_xxxxyz_yyyz, \
                             tr_x_xxxxyz_yyzz, \
                             tr_x_xxxxyz_yzzz, \
                             tr_x_xxxxyz_zzzz, \
                             tr_x_xxxxz_xxxx,  \
                             tr_x_xxxxz_xxxz,  \
                             tr_x_xxxxz_xxyz,  \
                             tr_x_xxxxz_xxz,   \
                             tr_x_xxxxz_xxzz,  \
                             tr_x_xxxxz_xyyz,  \
                             tr_x_xxxxz_xyz,   \
                             tr_x_xxxxz_xyzz,  \
                             tr_x_xxxxz_xzz,   \
                             tr_x_xxxxz_xzzz,  \
                             tr_x_xxxxz_yyyz,  \
                             tr_x_xxxxz_yyz,   \
                             tr_x_xxxxz_yyzz,  \
                             tr_x_xxxxz_yzz,   \
                             tr_x_xxxxz_yzzz,  \
                             tr_x_xxxxz_zzz,   \
                             tr_x_xxxxz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyz_xxxx[i] = tr_x_xxxxz_xxxx[i] * pa_y[i];

        tr_x_xxxxyz_xxxy[i] = tr_x_xxxxy_xxxy[i] * pa_z[i];

        tr_x_xxxxyz_xxxz[i] = tr_x_xxxxz_xxxz[i] * pa_y[i];

        tr_x_xxxxyz_xxyy[i] = tr_x_xxxxy_xxyy[i] * pa_z[i];

        tr_x_xxxxyz_xxyz[i] = tr_x_xxxxz_xxz[i] * fe_0 + tr_x_xxxxz_xxyz[i] * pa_y[i];

        tr_x_xxxxyz_xxzz[i] = tr_x_xxxxz_xxzz[i] * pa_y[i];

        tr_x_xxxxyz_xyyy[i] = tr_x_xxxxy_xyyy[i] * pa_z[i];

        tr_x_xxxxyz_xyyz[i] = 2.0 * tr_x_xxxxz_xyz[i] * fe_0 + tr_x_xxxxz_xyyz[i] * pa_y[i];

        tr_x_xxxxyz_xyzz[i] = tr_x_xxxxz_xzz[i] * fe_0 + tr_x_xxxxz_xyzz[i] * pa_y[i];

        tr_x_xxxxyz_xzzz[i] = tr_x_xxxxz_xzzz[i] * pa_y[i];

        tr_x_xxxxyz_yyyy[i] = tr_x_xxxxy_yyyy[i] * pa_z[i];

        tr_x_xxxxyz_yyyz[i] = 3.0 * tr_x_xxxxz_yyz[i] * fe_0 + tr_x_xxxxz_yyyz[i] * pa_y[i];

        tr_x_xxxxyz_yyzz[i] = 2.0 * tr_x_xxxxz_yzz[i] * fe_0 + tr_x_xxxxz_yyzz[i] * pa_y[i];

        tr_x_xxxxyz_yzzz[i] = tr_x_xxxxz_zzz[i] * fe_0 + tr_x_xxxxz_yzzz[i] * pa_y[i];

        tr_x_xxxxyz_zzzz[i] = tr_x_xxxxz_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : IG

    auto tr_x_xxxxzz_xxxx = pbuffer.data(idx_dip_ig + 75);

    auto tr_x_xxxxzz_xxxy = pbuffer.data(idx_dip_ig + 76);

    auto tr_x_xxxxzz_xxxz = pbuffer.data(idx_dip_ig + 77);

    auto tr_x_xxxxzz_xxyy = pbuffer.data(idx_dip_ig + 78);

    auto tr_x_xxxxzz_xxyz = pbuffer.data(idx_dip_ig + 79);

    auto tr_x_xxxxzz_xxzz = pbuffer.data(idx_dip_ig + 80);

    auto tr_x_xxxxzz_xyyy = pbuffer.data(idx_dip_ig + 81);

    auto tr_x_xxxxzz_xyyz = pbuffer.data(idx_dip_ig + 82);

    auto tr_x_xxxxzz_xyzz = pbuffer.data(idx_dip_ig + 83);

    auto tr_x_xxxxzz_xzzz = pbuffer.data(idx_dip_ig + 84);

    auto tr_x_xxxxzz_yyyy = pbuffer.data(idx_dip_ig + 85);

    auto tr_x_xxxxzz_yyyz = pbuffer.data(idx_dip_ig + 86);

    auto tr_x_xxxxzz_yyzz = pbuffer.data(idx_dip_ig + 87);

    auto tr_x_xxxxzz_yzzz = pbuffer.data(idx_dip_ig + 88);

    auto tr_x_xxxxzz_zzzz = pbuffer.data(idx_dip_ig + 89);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xxxx_xxxx,   \
                             tr_x_xxxx_xxxy,   \
                             tr_x_xxxx_xxxz,   \
                             tr_x_xxxx_xxyy,   \
                             tr_x_xxxx_xxyz,   \
                             tr_x_xxxx_xxzz,   \
                             tr_x_xxxx_xyyy,   \
                             tr_x_xxxx_xyyz,   \
                             tr_x_xxxx_xyzz,   \
                             tr_x_xxxx_xzzz,   \
                             tr_x_xxxx_yyyy,   \
                             tr_x_xxxxz_xxx,   \
                             tr_x_xxxxz_xxxx,  \
                             tr_x_xxxxz_xxxy,  \
                             tr_x_xxxxz_xxxz,  \
                             tr_x_xxxxz_xxy,   \
                             tr_x_xxxxz_xxyy,  \
                             tr_x_xxxxz_xxyz,  \
                             tr_x_xxxxz_xxz,   \
                             tr_x_xxxxz_xxzz,  \
                             tr_x_xxxxz_xyy,   \
                             tr_x_xxxxz_xyyy,  \
                             tr_x_xxxxz_xyyz,  \
                             tr_x_xxxxz_xyz,   \
                             tr_x_xxxxz_xyzz,  \
                             tr_x_xxxxz_xzz,   \
                             tr_x_xxxxz_xzzz,  \
                             tr_x_xxxxz_yyyy,  \
                             tr_x_xxxxzz_xxxx, \
                             tr_x_xxxxzz_xxxy, \
                             tr_x_xxxxzz_xxxz, \
                             tr_x_xxxxzz_xxyy, \
                             tr_x_xxxxzz_xxyz, \
                             tr_x_xxxxzz_xxzz, \
                             tr_x_xxxxzz_xyyy, \
                             tr_x_xxxxzz_xyyz, \
                             tr_x_xxxxzz_xyzz, \
                             tr_x_xxxxzz_xzzz, \
                             tr_x_xxxxzz_yyyy, \
                             tr_x_xxxxzz_yyyz, \
                             tr_x_xxxxzz_yyzz, \
                             tr_x_xxxxzz_yzzz, \
                             tr_x_xxxxzz_zzzz, \
                             tr_x_xxxzz_yyyz,  \
                             tr_x_xxxzz_yyzz,  \
                             tr_x_xxxzz_yzzz,  \
                             tr_x_xxxzz_zzzz,  \
                             tr_x_xxzz_yyyz,   \
                             tr_x_xxzz_yyzz,   \
                             tr_x_xxzz_yzzz,   \
                             tr_x_xxzz_zzzz,   \
                             ts_xxxzz_yyyz,    \
                             ts_xxxzz_yyzz,    \
                             ts_xxxzz_yzzz,    \
                             ts_xxxzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxzz_xxxx[i] = tr_x_xxxx_xxxx[i] * fe_0 + tr_x_xxxxz_xxxx[i] * pa_z[i];

        tr_x_xxxxzz_xxxy[i] = tr_x_xxxx_xxxy[i] * fe_0 + tr_x_xxxxz_xxxy[i] * pa_z[i];

        tr_x_xxxxzz_xxxz[i] = tr_x_xxxx_xxxz[i] * fe_0 + tr_x_xxxxz_xxx[i] * fe_0 + tr_x_xxxxz_xxxz[i] * pa_z[i];

        tr_x_xxxxzz_xxyy[i] = tr_x_xxxx_xxyy[i] * fe_0 + tr_x_xxxxz_xxyy[i] * pa_z[i];

        tr_x_xxxxzz_xxyz[i] = tr_x_xxxx_xxyz[i] * fe_0 + tr_x_xxxxz_xxy[i] * fe_0 + tr_x_xxxxz_xxyz[i] * pa_z[i];

        tr_x_xxxxzz_xxzz[i] = tr_x_xxxx_xxzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xxz[i] * fe_0 + tr_x_xxxxz_xxzz[i] * pa_z[i];

        tr_x_xxxxzz_xyyy[i] = tr_x_xxxx_xyyy[i] * fe_0 + tr_x_xxxxz_xyyy[i] * pa_z[i];

        tr_x_xxxxzz_xyyz[i] = tr_x_xxxx_xyyz[i] * fe_0 + tr_x_xxxxz_xyy[i] * fe_0 + tr_x_xxxxz_xyyz[i] * pa_z[i];

        tr_x_xxxxzz_xyzz[i] = tr_x_xxxx_xyzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xyz[i] * fe_0 + tr_x_xxxxz_xyzz[i] * pa_z[i];

        tr_x_xxxxzz_xzzz[i] = tr_x_xxxx_xzzz[i] * fe_0 + 3.0 * tr_x_xxxxz_xzz[i] * fe_0 + tr_x_xxxxz_xzzz[i] * pa_z[i];

        tr_x_xxxxzz_yyyy[i] = tr_x_xxxx_yyyy[i] * fe_0 + tr_x_xxxxz_yyyy[i] * pa_z[i];

        tr_x_xxxxzz_yyyz[i] = 3.0 * tr_x_xxzz_yyyz[i] * fe_0 + ts_xxxzz_yyyz[i] * fe_0 + tr_x_xxxzz_yyyz[i] * pa_x[i];

        tr_x_xxxxzz_yyzz[i] = 3.0 * tr_x_xxzz_yyzz[i] * fe_0 + ts_xxxzz_yyzz[i] * fe_0 + tr_x_xxxzz_yyzz[i] * pa_x[i];

        tr_x_xxxxzz_yzzz[i] = 3.0 * tr_x_xxzz_yzzz[i] * fe_0 + ts_xxxzz_yzzz[i] * fe_0 + tr_x_xxxzz_yzzz[i] * pa_x[i];

        tr_x_xxxxzz_zzzz[i] = 3.0 * tr_x_xxzz_zzzz[i] * fe_0 + ts_xxxzz_zzzz[i] * fe_0 + tr_x_xxxzz_zzzz[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : IG

    auto tr_x_xxxyyy_xxxx = pbuffer.data(idx_dip_ig + 90);

    auto tr_x_xxxyyy_xxxy = pbuffer.data(idx_dip_ig + 91);

    auto tr_x_xxxyyy_xxxz = pbuffer.data(idx_dip_ig + 92);

    auto tr_x_xxxyyy_xxyy = pbuffer.data(idx_dip_ig + 93);

    auto tr_x_xxxyyy_xxyz = pbuffer.data(idx_dip_ig + 94);

    auto tr_x_xxxyyy_xxzz = pbuffer.data(idx_dip_ig + 95);

    auto tr_x_xxxyyy_xyyy = pbuffer.data(idx_dip_ig + 96);

    auto tr_x_xxxyyy_xyyz = pbuffer.data(idx_dip_ig + 97);

    auto tr_x_xxxyyy_xyzz = pbuffer.data(idx_dip_ig + 98);

    auto tr_x_xxxyyy_xzzz = pbuffer.data(idx_dip_ig + 99);

    auto tr_x_xxxyyy_yyyy = pbuffer.data(idx_dip_ig + 100);

    auto tr_x_xxxyyy_yyyz = pbuffer.data(idx_dip_ig + 101);

    auto tr_x_xxxyyy_yyzz = pbuffer.data(idx_dip_ig + 102);

    auto tr_x_xxxyyy_yzzz = pbuffer.data(idx_dip_ig + 103);

    auto tr_x_xxxyyy_zzzz = pbuffer.data(idx_dip_ig + 104);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xxxy_xxxx,   \
                             tr_x_xxxy_xxxy,   \
                             tr_x_xxxy_xxxz,   \
                             tr_x_xxxy_xxyy,   \
                             tr_x_xxxy_xxyz,   \
                             tr_x_xxxy_xxzz,   \
                             tr_x_xxxy_xyyy,   \
                             tr_x_xxxy_xyyz,   \
                             tr_x_xxxy_xyzz,   \
                             tr_x_xxxy_xzzz,   \
                             tr_x_xxxy_zzzz,   \
                             tr_x_xxxyy_xxx,   \
                             tr_x_xxxyy_xxxx,  \
                             tr_x_xxxyy_xxxy,  \
                             tr_x_xxxyy_xxxz,  \
                             tr_x_xxxyy_xxy,   \
                             tr_x_xxxyy_xxyy,  \
                             tr_x_xxxyy_xxyz,  \
                             tr_x_xxxyy_xxz,   \
                             tr_x_xxxyy_xxzz,  \
                             tr_x_xxxyy_xyy,   \
                             tr_x_xxxyy_xyyy,  \
                             tr_x_xxxyy_xyyz,  \
                             tr_x_xxxyy_xyz,   \
                             tr_x_xxxyy_xyzz,  \
                             tr_x_xxxyy_xzz,   \
                             tr_x_xxxyy_xzzz,  \
                             tr_x_xxxyy_zzzz,  \
                             tr_x_xxxyyy_xxxx, \
                             tr_x_xxxyyy_xxxy, \
                             tr_x_xxxyyy_xxxz, \
                             tr_x_xxxyyy_xxyy, \
                             tr_x_xxxyyy_xxyz, \
                             tr_x_xxxyyy_xxzz, \
                             tr_x_xxxyyy_xyyy, \
                             tr_x_xxxyyy_xyyz, \
                             tr_x_xxxyyy_xyzz, \
                             tr_x_xxxyyy_xzzz, \
                             tr_x_xxxyyy_yyyy, \
                             tr_x_xxxyyy_yyyz, \
                             tr_x_xxxyyy_yyzz, \
                             tr_x_xxxyyy_yzzz, \
                             tr_x_xxxyyy_zzzz, \
                             tr_x_xxyyy_yyyy,  \
                             tr_x_xxyyy_yyyz,  \
                             tr_x_xxyyy_yyzz,  \
                             tr_x_xxyyy_yzzz,  \
                             tr_x_xyyy_yyyy,   \
                             tr_x_xyyy_yyyz,   \
                             tr_x_xyyy_yyzz,   \
                             tr_x_xyyy_yzzz,   \
                             ts_xxyyy_yyyy,    \
                             ts_xxyyy_yyyz,    \
                             ts_xxyyy_yyzz,    \
                             ts_xxyyy_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyy_xxxx[i] = 2.0 * tr_x_xxxy_xxxx[i] * fe_0 + tr_x_xxxyy_xxxx[i] * pa_y[i];

        tr_x_xxxyyy_xxxy[i] = 2.0 * tr_x_xxxy_xxxy[i] * fe_0 + tr_x_xxxyy_xxx[i] * fe_0 + tr_x_xxxyy_xxxy[i] * pa_y[i];

        tr_x_xxxyyy_xxxz[i] = 2.0 * tr_x_xxxy_xxxz[i] * fe_0 + tr_x_xxxyy_xxxz[i] * pa_y[i];

        tr_x_xxxyyy_xxyy[i] = 2.0 * tr_x_xxxy_xxyy[i] * fe_0 + 2.0 * tr_x_xxxyy_xxy[i] * fe_0 + tr_x_xxxyy_xxyy[i] * pa_y[i];

        tr_x_xxxyyy_xxyz[i] = 2.0 * tr_x_xxxy_xxyz[i] * fe_0 + tr_x_xxxyy_xxz[i] * fe_0 + tr_x_xxxyy_xxyz[i] * pa_y[i];

        tr_x_xxxyyy_xxzz[i] = 2.0 * tr_x_xxxy_xxzz[i] * fe_0 + tr_x_xxxyy_xxzz[i] * pa_y[i];

        tr_x_xxxyyy_xyyy[i] = 2.0 * tr_x_xxxy_xyyy[i] * fe_0 + 3.0 * tr_x_xxxyy_xyy[i] * fe_0 + tr_x_xxxyy_xyyy[i] * pa_y[i];

        tr_x_xxxyyy_xyyz[i] = 2.0 * tr_x_xxxy_xyyz[i] * fe_0 + 2.0 * tr_x_xxxyy_xyz[i] * fe_0 + tr_x_xxxyy_xyyz[i] * pa_y[i];

        tr_x_xxxyyy_xyzz[i] = 2.0 * tr_x_xxxy_xyzz[i] * fe_0 + tr_x_xxxyy_xzz[i] * fe_0 + tr_x_xxxyy_xyzz[i] * pa_y[i];

        tr_x_xxxyyy_xzzz[i] = 2.0 * tr_x_xxxy_xzzz[i] * fe_0 + tr_x_xxxyy_xzzz[i] * pa_y[i];

        tr_x_xxxyyy_yyyy[i] = 2.0 * tr_x_xyyy_yyyy[i] * fe_0 + ts_xxyyy_yyyy[i] * fe_0 + tr_x_xxyyy_yyyy[i] * pa_x[i];

        tr_x_xxxyyy_yyyz[i] = 2.0 * tr_x_xyyy_yyyz[i] * fe_0 + ts_xxyyy_yyyz[i] * fe_0 + tr_x_xxyyy_yyyz[i] * pa_x[i];

        tr_x_xxxyyy_yyzz[i] = 2.0 * tr_x_xyyy_yyzz[i] * fe_0 + ts_xxyyy_yyzz[i] * fe_0 + tr_x_xxyyy_yyzz[i] * pa_x[i];

        tr_x_xxxyyy_yzzz[i] = 2.0 * tr_x_xyyy_yzzz[i] * fe_0 + ts_xxyyy_yzzz[i] * fe_0 + tr_x_xxyyy_yzzz[i] * pa_x[i];

        tr_x_xxxyyy_zzzz[i] = 2.0 * tr_x_xxxy_zzzz[i] * fe_0 + tr_x_xxxyy_zzzz[i] * pa_y[i];
    }

    // Set up 105-120 components of targeted buffer : IG

    auto tr_x_xxxyyz_xxxx = pbuffer.data(idx_dip_ig + 105);

    auto tr_x_xxxyyz_xxxy = pbuffer.data(idx_dip_ig + 106);

    auto tr_x_xxxyyz_xxxz = pbuffer.data(idx_dip_ig + 107);

    auto tr_x_xxxyyz_xxyy = pbuffer.data(idx_dip_ig + 108);

    auto tr_x_xxxyyz_xxyz = pbuffer.data(idx_dip_ig + 109);

    auto tr_x_xxxyyz_xxzz = pbuffer.data(idx_dip_ig + 110);

    auto tr_x_xxxyyz_xyyy = pbuffer.data(idx_dip_ig + 111);

    auto tr_x_xxxyyz_xyyz = pbuffer.data(idx_dip_ig + 112);

    auto tr_x_xxxyyz_xyzz = pbuffer.data(idx_dip_ig + 113);

    auto tr_x_xxxyyz_xzzz = pbuffer.data(idx_dip_ig + 114);

    auto tr_x_xxxyyz_yyyy = pbuffer.data(idx_dip_ig + 115);

    auto tr_x_xxxyyz_yyyz = pbuffer.data(idx_dip_ig + 116);

    auto tr_x_xxxyyz_yyzz = pbuffer.data(idx_dip_ig + 117);

    auto tr_x_xxxyyz_yzzz = pbuffer.data(idx_dip_ig + 118);

    auto tr_x_xxxyyz_zzzz = pbuffer.data(idx_dip_ig + 119);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_xxxyy_xxxx,  \
                             tr_x_xxxyy_xxxy,  \
                             tr_x_xxxyy_xxy,   \
                             tr_x_xxxyy_xxyy,  \
                             tr_x_xxxyy_xxyz,  \
                             tr_x_xxxyy_xyy,   \
                             tr_x_xxxyy_xyyy,  \
                             tr_x_xxxyy_xyyz,  \
                             tr_x_xxxyy_xyz,   \
                             tr_x_xxxyy_xyzz,  \
                             tr_x_xxxyy_yyy,   \
                             tr_x_xxxyy_yyyy,  \
                             tr_x_xxxyy_yyyz,  \
                             tr_x_xxxyy_yyz,   \
                             tr_x_xxxyy_yyzz,  \
                             tr_x_xxxyy_yzz,   \
                             tr_x_xxxyy_yzzz,  \
                             tr_x_xxxyyz_xxxx, \
                             tr_x_xxxyyz_xxxy, \
                             tr_x_xxxyyz_xxxz, \
                             tr_x_xxxyyz_xxyy, \
                             tr_x_xxxyyz_xxyz, \
                             tr_x_xxxyyz_xxzz, \
                             tr_x_xxxyyz_xyyy, \
                             tr_x_xxxyyz_xyyz, \
                             tr_x_xxxyyz_xyzz, \
                             tr_x_xxxyyz_xzzz, \
                             tr_x_xxxyyz_yyyy, \
                             tr_x_xxxyyz_yyyz, \
                             tr_x_xxxyyz_yyzz, \
                             tr_x_xxxyyz_yzzz, \
                             tr_x_xxxyyz_zzzz, \
                             tr_x_xxxyz_xxxz,  \
                             tr_x_xxxyz_xxzz,  \
                             tr_x_xxxyz_xzzz,  \
                             tr_x_xxxyz_zzzz,  \
                             tr_x_xxxz_xxxz,   \
                             tr_x_xxxz_xxzz,   \
                             tr_x_xxxz_xzzz,   \
                             tr_x_xxxz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyz_xxxx[i] = tr_x_xxxyy_xxxx[i] * pa_z[i];

        tr_x_xxxyyz_xxxy[i] = tr_x_xxxyy_xxxy[i] * pa_z[i];

        tr_x_xxxyyz_xxxz[i] = tr_x_xxxz_xxxz[i] * fe_0 + tr_x_xxxyz_xxxz[i] * pa_y[i];

        tr_x_xxxyyz_xxyy[i] = tr_x_xxxyy_xxyy[i] * pa_z[i];

        tr_x_xxxyyz_xxyz[i] = tr_x_xxxyy_xxy[i] * fe_0 + tr_x_xxxyy_xxyz[i] * pa_z[i];

        tr_x_xxxyyz_xxzz[i] = tr_x_xxxz_xxzz[i] * fe_0 + tr_x_xxxyz_xxzz[i] * pa_y[i];

        tr_x_xxxyyz_xyyy[i] = tr_x_xxxyy_xyyy[i] * pa_z[i];

        tr_x_xxxyyz_xyyz[i] = tr_x_xxxyy_xyy[i] * fe_0 + tr_x_xxxyy_xyyz[i] * pa_z[i];

        tr_x_xxxyyz_xyzz[i] = 2.0 * tr_x_xxxyy_xyz[i] * fe_0 + tr_x_xxxyy_xyzz[i] * pa_z[i];

        tr_x_xxxyyz_xzzz[i] = tr_x_xxxz_xzzz[i] * fe_0 + tr_x_xxxyz_xzzz[i] * pa_y[i];

        tr_x_xxxyyz_yyyy[i] = tr_x_xxxyy_yyyy[i] * pa_z[i];

        tr_x_xxxyyz_yyyz[i] = tr_x_xxxyy_yyy[i] * fe_0 + tr_x_xxxyy_yyyz[i] * pa_z[i];

        tr_x_xxxyyz_yyzz[i] = 2.0 * tr_x_xxxyy_yyz[i] * fe_0 + tr_x_xxxyy_yyzz[i] * pa_z[i];

        tr_x_xxxyyz_yzzz[i] = 3.0 * tr_x_xxxyy_yzz[i] * fe_0 + tr_x_xxxyy_yzzz[i] * pa_z[i];

        tr_x_xxxyyz_zzzz[i] = tr_x_xxxz_zzzz[i] * fe_0 + tr_x_xxxyz_zzzz[i] * pa_y[i];
    }

    // Set up 120-135 components of targeted buffer : IG

    auto tr_x_xxxyzz_xxxx = pbuffer.data(idx_dip_ig + 120);

    auto tr_x_xxxyzz_xxxy = pbuffer.data(idx_dip_ig + 121);

    auto tr_x_xxxyzz_xxxz = pbuffer.data(idx_dip_ig + 122);

    auto tr_x_xxxyzz_xxyy = pbuffer.data(idx_dip_ig + 123);

    auto tr_x_xxxyzz_xxyz = pbuffer.data(idx_dip_ig + 124);

    auto tr_x_xxxyzz_xxzz = pbuffer.data(idx_dip_ig + 125);

    auto tr_x_xxxyzz_xyyy = pbuffer.data(idx_dip_ig + 126);

    auto tr_x_xxxyzz_xyyz = pbuffer.data(idx_dip_ig + 127);

    auto tr_x_xxxyzz_xyzz = pbuffer.data(idx_dip_ig + 128);

    auto tr_x_xxxyzz_xzzz = pbuffer.data(idx_dip_ig + 129);

    auto tr_x_xxxyzz_yyyy = pbuffer.data(idx_dip_ig + 130);

    auto tr_x_xxxyzz_yyyz = pbuffer.data(idx_dip_ig + 131);

    auto tr_x_xxxyzz_yyzz = pbuffer.data(idx_dip_ig + 132);

    auto tr_x_xxxyzz_yzzz = pbuffer.data(idx_dip_ig + 133);

    auto tr_x_xxxyzz_zzzz = pbuffer.data(idx_dip_ig + 134);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_xxxyzz_xxxx, \
                             tr_x_xxxyzz_xxxy, \
                             tr_x_xxxyzz_xxxz, \
                             tr_x_xxxyzz_xxyy, \
                             tr_x_xxxyzz_xxyz, \
                             tr_x_xxxyzz_xxzz, \
                             tr_x_xxxyzz_xyyy, \
                             tr_x_xxxyzz_xyyz, \
                             tr_x_xxxyzz_xyzz, \
                             tr_x_xxxyzz_xzzz, \
                             tr_x_xxxyzz_yyyy, \
                             tr_x_xxxyzz_yyyz, \
                             tr_x_xxxyzz_yyzz, \
                             tr_x_xxxyzz_yzzz, \
                             tr_x_xxxyzz_zzzz, \
                             tr_x_xxxzz_xxx,   \
                             tr_x_xxxzz_xxxx,  \
                             tr_x_xxxzz_xxxy,  \
                             tr_x_xxxzz_xxxz,  \
                             tr_x_xxxzz_xxy,   \
                             tr_x_xxxzz_xxyy,  \
                             tr_x_xxxzz_xxyz,  \
                             tr_x_xxxzz_xxz,   \
                             tr_x_xxxzz_xxzz,  \
                             tr_x_xxxzz_xyy,   \
                             tr_x_xxxzz_xyyy,  \
                             tr_x_xxxzz_xyyz,  \
                             tr_x_xxxzz_xyz,   \
                             tr_x_xxxzz_xyzz,  \
                             tr_x_xxxzz_xzz,   \
                             tr_x_xxxzz_xzzz,  \
                             tr_x_xxxzz_yyy,   \
                             tr_x_xxxzz_yyyy,  \
                             tr_x_xxxzz_yyyz,  \
                             tr_x_xxxzz_yyz,   \
                             tr_x_xxxzz_yyzz,  \
                             tr_x_xxxzz_yzz,   \
                             tr_x_xxxzz_yzzz,  \
                             tr_x_xxxzz_zzz,   \
                             tr_x_xxxzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyzz_xxxx[i] = tr_x_xxxzz_xxxx[i] * pa_y[i];

        tr_x_xxxyzz_xxxy[i] = tr_x_xxxzz_xxx[i] * fe_0 + tr_x_xxxzz_xxxy[i] * pa_y[i];

        tr_x_xxxyzz_xxxz[i] = tr_x_xxxzz_xxxz[i] * pa_y[i];

        tr_x_xxxyzz_xxyy[i] = 2.0 * tr_x_xxxzz_xxy[i] * fe_0 + tr_x_xxxzz_xxyy[i] * pa_y[i];

        tr_x_xxxyzz_xxyz[i] = tr_x_xxxzz_xxz[i] * fe_0 + tr_x_xxxzz_xxyz[i] * pa_y[i];

        tr_x_xxxyzz_xxzz[i] = tr_x_xxxzz_xxzz[i] * pa_y[i];

        tr_x_xxxyzz_xyyy[i] = 3.0 * tr_x_xxxzz_xyy[i] * fe_0 + tr_x_xxxzz_xyyy[i] * pa_y[i];

        tr_x_xxxyzz_xyyz[i] = 2.0 * tr_x_xxxzz_xyz[i] * fe_0 + tr_x_xxxzz_xyyz[i] * pa_y[i];

        tr_x_xxxyzz_xyzz[i] = tr_x_xxxzz_xzz[i] * fe_0 + tr_x_xxxzz_xyzz[i] * pa_y[i];

        tr_x_xxxyzz_xzzz[i] = tr_x_xxxzz_xzzz[i] * pa_y[i];

        tr_x_xxxyzz_yyyy[i] = 4.0 * tr_x_xxxzz_yyy[i] * fe_0 + tr_x_xxxzz_yyyy[i] * pa_y[i];

        tr_x_xxxyzz_yyyz[i] = 3.0 * tr_x_xxxzz_yyz[i] * fe_0 + tr_x_xxxzz_yyyz[i] * pa_y[i];

        tr_x_xxxyzz_yyzz[i] = 2.0 * tr_x_xxxzz_yzz[i] * fe_0 + tr_x_xxxzz_yyzz[i] * pa_y[i];

        tr_x_xxxyzz_yzzz[i] = tr_x_xxxzz_zzz[i] * fe_0 + tr_x_xxxzz_yzzz[i] * pa_y[i];

        tr_x_xxxyzz_zzzz[i] = tr_x_xxxzz_zzzz[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : IG

    auto tr_x_xxxzzz_xxxx = pbuffer.data(idx_dip_ig + 135);

    auto tr_x_xxxzzz_xxxy = pbuffer.data(idx_dip_ig + 136);

    auto tr_x_xxxzzz_xxxz = pbuffer.data(idx_dip_ig + 137);

    auto tr_x_xxxzzz_xxyy = pbuffer.data(idx_dip_ig + 138);

    auto tr_x_xxxzzz_xxyz = pbuffer.data(idx_dip_ig + 139);

    auto tr_x_xxxzzz_xxzz = pbuffer.data(idx_dip_ig + 140);

    auto tr_x_xxxzzz_xyyy = pbuffer.data(idx_dip_ig + 141);

    auto tr_x_xxxzzz_xyyz = pbuffer.data(idx_dip_ig + 142);

    auto tr_x_xxxzzz_xyzz = pbuffer.data(idx_dip_ig + 143);

    auto tr_x_xxxzzz_xzzz = pbuffer.data(idx_dip_ig + 144);

    auto tr_x_xxxzzz_yyyy = pbuffer.data(idx_dip_ig + 145);

    auto tr_x_xxxzzz_yyyz = pbuffer.data(idx_dip_ig + 146);

    auto tr_x_xxxzzz_yyzz = pbuffer.data(idx_dip_ig + 147);

    auto tr_x_xxxzzz_yzzz = pbuffer.data(idx_dip_ig + 148);

    auto tr_x_xxxzzz_zzzz = pbuffer.data(idx_dip_ig + 149);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xxxz_xxxx,   \
                             tr_x_xxxz_xxxy,   \
                             tr_x_xxxz_xxxz,   \
                             tr_x_xxxz_xxyy,   \
                             tr_x_xxxz_xxyz,   \
                             tr_x_xxxz_xxzz,   \
                             tr_x_xxxz_xyyy,   \
                             tr_x_xxxz_xyyz,   \
                             tr_x_xxxz_xyzz,   \
                             tr_x_xxxz_xzzz,   \
                             tr_x_xxxz_yyyy,   \
                             tr_x_xxxzz_xxx,   \
                             tr_x_xxxzz_xxxx,  \
                             tr_x_xxxzz_xxxy,  \
                             tr_x_xxxzz_xxxz,  \
                             tr_x_xxxzz_xxy,   \
                             tr_x_xxxzz_xxyy,  \
                             tr_x_xxxzz_xxyz,  \
                             tr_x_xxxzz_xxz,   \
                             tr_x_xxxzz_xxzz,  \
                             tr_x_xxxzz_xyy,   \
                             tr_x_xxxzz_xyyy,  \
                             tr_x_xxxzz_xyyz,  \
                             tr_x_xxxzz_xyz,   \
                             tr_x_xxxzz_xyzz,  \
                             tr_x_xxxzz_xzz,   \
                             tr_x_xxxzz_xzzz,  \
                             tr_x_xxxzz_yyyy,  \
                             tr_x_xxxzzz_xxxx, \
                             tr_x_xxxzzz_xxxy, \
                             tr_x_xxxzzz_xxxz, \
                             tr_x_xxxzzz_xxyy, \
                             tr_x_xxxzzz_xxyz, \
                             tr_x_xxxzzz_xxzz, \
                             tr_x_xxxzzz_xyyy, \
                             tr_x_xxxzzz_xyyz, \
                             tr_x_xxxzzz_xyzz, \
                             tr_x_xxxzzz_xzzz, \
                             tr_x_xxxzzz_yyyy, \
                             tr_x_xxxzzz_yyyz, \
                             tr_x_xxxzzz_yyzz, \
                             tr_x_xxxzzz_yzzz, \
                             tr_x_xxxzzz_zzzz, \
                             tr_x_xxzzz_yyyz,  \
                             tr_x_xxzzz_yyzz,  \
                             tr_x_xxzzz_yzzz,  \
                             tr_x_xxzzz_zzzz,  \
                             tr_x_xzzz_yyyz,   \
                             tr_x_xzzz_yyzz,   \
                             tr_x_xzzz_yzzz,   \
                             tr_x_xzzz_zzzz,   \
                             ts_xxzzz_yyyz,    \
                             ts_xxzzz_yyzz,    \
                             ts_xxzzz_yzzz,    \
                             ts_xxzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzzz_xxxx[i] = 2.0 * tr_x_xxxz_xxxx[i] * fe_0 + tr_x_xxxzz_xxxx[i] * pa_z[i];

        tr_x_xxxzzz_xxxy[i] = 2.0 * tr_x_xxxz_xxxy[i] * fe_0 + tr_x_xxxzz_xxxy[i] * pa_z[i];

        tr_x_xxxzzz_xxxz[i] = 2.0 * tr_x_xxxz_xxxz[i] * fe_0 + tr_x_xxxzz_xxx[i] * fe_0 + tr_x_xxxzz_xxxz[i] * pa_z[i];

        tr_x_xxxzzz_xxyy[i] = 2.0 * tr_x_xxxz_xxyy[i] * fe_0 + tr_x_xxxzz_xxyy[i] * pa_z[i];

        tr_x_xxxzzz_xxyz[i] = 2.0 * tr_x_xxxz_xxyz[i] * fe_0 + tr_x_xxxzz_xxy[i] * fe_0 + tr_x_xxxzz_xxyz[i] * pa_z[i];

        tr_x_xxxzzz_xxzz[i] = 2.0 * tr_x_xxxz_xxzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xxz[i] * fe_0 + tr_x_xxxzz_xxzz[i] * pa_z[i];

        tr_x_xxxzzz_xyyy[i] = 2.0 * tr_x_xxxz_xyyy[i] * fe_0 + tr_x_xxxzz_xyyy[i] * pa_z[i];

        tr_x_xxxzzz_xyyz[i] = 2.0 * tr_x_xxxz_xyyz[i] * fe_0 + tr_x_xxxzz_xyy[i] * fe_0 + tr_x_xxxzz_xyyz[i] * pa_z[i];

        tr_x_xxxzzz_xyzz[i] = 2.0 * tr_x_xxxz_xyzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xyz[i] * fe_0 + tr_x_xxxzz_xyzz[i] * pa_z[i];

        tr_x_xxxzzz_xzzz[i] = 2.0 * tr_x_xxxz_xzzz[i] * fe_0 + 3.0 * tr_x_xxxzz_xzz[i] * fe_0 + tr_x_xxxzz_xzzz[i] * pa_z[i];

        tr_x_xxxzzz_yyyy[i] = 2.0 * tr_x_xxxz_yyyy[i] * fe_0 + tr_x_xxxzz_yyyy[i] * pa_z[i];

        tr_x_xxxzzz_yyyz[i] = 2.0 * tr_x_xzzz_yyyz[i] * fe_0 + ts_xxzzz_yyyz[i] * fe_0 + tr_x_xxzzz_yyyz[i] * pa_x[i];

        tr_x_xxxzzz_yyzz[i] = 2.0 * tr_x_xzzz_yyzz[i] * fe_0 + ts_xxzzz_yyzz[i] * fe_0 + tr_x_xxzzz_yyzz[i] * pa_x[i];

        tr_x_xxxzzz_yzzz[i] = 2.0 * tr_x_xzzz_yzzz[i] * fe_0 + ts_xxzzz_yzzz[i] * fe_0 + tr_x_xxzzz_yzzz[i] * pa_x[i];

        tr_x_xxxzzz_zzzz[i] = 2.0 * tr_x_xzzz_zzzz[i] * fe_0 + ts_xxzzz_zzzz[i] * fe_0 + tr_x_xxzzz_zzzz[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : IG

    auto tr_x_xxyyyy_xxxx = pbuffer.data(idx_dip_ig + 150);

    auto tr_x_xxyyyy_xxxy = pbuffer.data(idx_dip_ig + 151);

    auto tr_x_xxyyyy_xxxz = pbuffer.data(idx_dip_ig + 152);

    auto tr_x_xxyyyy_xxyy = pbuffer.data(idx_dip_ig + 153);

    auto tr_x_xxyyyy_xxyz = pbuffer.data(idx_dip_ig + 154);

    auto tr_x_xxyyyy_xxzz = pbuffer.data(idx_dip_ig + 155);

    auto tr_x_xxyyyy_xyyy = pbuffer.data(idx_dip_ig + 156);

    auto tr_x_xxyyyy_xyyz = pbuffer.data(idx_dip_ig + 157);

    auto tr_x_xxyyyy_xyzz = pbuffer.data(idx_dip_ig + 158);

    auto tr_x_xxyyyy_xzzz = pbuffer.data(idx_dip_ig + 159);

    auto tr_x_xxyyyy_yyyy = pbuffer.data(idx_dip_ig + 160);

    auto tr_x_xxyyyy_yyyz = pbuffer.data(idx_dip_ig + 161);

    auto tr_x_xxyyyy_yyzz = pbuffer.data(idx_dip_ig + 162);

    auto tr_x_xxyyyy_yzzz = pbuffer.data(idx_dip_ig + 163);

    auto tr_x_xxyyyy_zzzz = pbuffer.data(idx_dip_ig + 164);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xxyy_xxxx,   \
                             tr_x_xxyy_xxxy,   \
                             tr_x_xxyy_xxxz,   \
                             tr_x_xxyy_xxyy,   \
                             tr_x_xxyy_xxyz,   \
                             tr_x_xxyy_xxzz,   \
                             tr_x_xxyy_xyyy,   \
                             tr_x_xxyy_xyyz,   \
                             tr_x_xxyy_xyzz,   \
                             tr_x_xxyy_xzzz,   \
                             tr_x_xxyy_zzzz,   \
                             tr_x_xxyyy_xxx,   \
                             tr_x_xxyyy_xxxx,  \
                             tr_x_xxyyy_xxxy,  \
                             tr_x_xxyyy_xxxz,  \
                             tr_x_xxyyy_xxy,   \
                             tr_x_xxyyy_xxyy,  \
                             tr_x_xxyyy_xxyz,  \
                             tr_x_xxyyy_xxz,   \
                             tr_x_xxyyy_xxzz,  \
                             tr_x_xxyyy_xyy,   \
                             tr_x_xxyyy_xyyy,  \
                             tr_x_xxyyy_xyyz,  \
                             tr_x_xxyyy_xyz,   \
                             tr_x_xxyyy_xyzz,  \
                             tr_x_xxyyy_xzz,   \
                             tr_x_xxyyy_xzzz,  \
                             tr_x_xxyyy_zzzz,  \
                             tr_x_xxyyyy_xxxx, \
                             tr_x_xxyyyy_xxxy, \
                             tr_x_xxyyyy_xxxz, \
                             tr_x_xxyyyy_xxyy, \
                             tr_x_xxyyyy_xxyz, \
                             tr_x_xxyyyy_xxzz, \
                             tr_x_xxyyyy_xyyy, \
                             tr_x_xxyyyy_xyyz, \
                             tr_x_xxyyyy_xyzz, \
                             tr_x_xxyyyy_xzzz, \
                             tr_x_xxyyyy_yyyy, \
                             tr_x_xxyyyy_yyyz, \
                             tr_x_xxyyyy_yyzz, \
                             tr_x_xxyyyy_yzzz, \
                             tr_x_xxyyyy_zzzz, \
                             tr_x_xyyyy_yyyy,  \
                             tr_x_xyyyy_yyyz,  \
                             tr_x_xyyyy_yyzz,  \
                             tr_x_xyyyy_yzzz,  \
                             tr_x_yyyy_yyyy,   \
                             tr_x_yyyy_yyyz,   \
                             tr_x_yyyy_yyzz,   \
                             tr_x_yyyy_yzzz,   \
                             ts_xyyyy_yyyy,    \
                             ts_xyyyy_yyyz,    \
                             ts_xyyyy_yyzz,    \
                             ts_xyyyy_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyy_xxxx[i] = 3.0 * tr_x_xxyy_xxxx[i] * fe_0 + tr_x_xxyyy_xxxx[i] * pa_y[i];

        tr_x_xxyyyy_xxxy[i] = 3.0 * tr_x_xxyy_xxxy[i] * fe_0 + tr_x_xxyyy_xxx[i] * fe_0 + tr_x_xxyyy_xxxy[i] * pa_y[i];

        tr_x_xxyyyy_xxxz[i] = 3.0 * tr_x_xxyy_xxxz[i] * fe_0 + tr_x_xxyyy_xxxz[i] * pa_y[i];

        tr_x_xxyyyy_xxyy[i] = 3.0 * tr_x_xxyy_xxyy[i] * fe_0 + 2.0 * tr_x_xxyyy_xxy[i] * fe_0 + tr_x_xxyyy_xxyy[i] * pa_y[i];

        tr_x_xxyyyy_xxyz[i] = 3.0 * tr_x_xxyy_xxyz[i] * fe_0 + tr_x_xxyyy_xxz[i] * fe_0 + tr_x_xxyyy_xxyz[i] * pa_y[i];

        tr_x_xxyyyy_xxzz[i] = 3.0 * tr_x_xxyy_xxzz[i] * fe_0 + tr_x_xxyyy_xxzz[i] * pa_y[i];

        tr_x_xxyyyy_xyyy[i] = 3.0 * tr_x_xxyy_xyyy[i] * fe_0 + 3.0 * tr_x_xxyyy_xyy[i] * fe_0 + tr_x_xxyyy_xyyy[i] * pa_y[i];

        tr_x_xxyyyy_xyyz[i] = 3.0 * tr_x_xxyy_xyyz[i] * fe_0 + 2.0 * tr_x_xxyyy_xyz[i] * fe_0 + tr_x_xxyyy_xyyz[i] * pa_y[i];

        tr_x_xxyyyy_xyzz[i] = 3.0 * tr_x_xxyy_xyzz[i] * fe_0 + tr_x_xxyyy_xzz[i] * fe_0 + tr_x_xxyyy_xyzz[i] * pa_y[i];

        tr_x_xxyyyy_xzzz[i] = 3.0 * tr_x_xxyy_xzzz[i] * fe_0 + tr_x_xxyyy_xzzz[i] * pa_y[i];

        tr_x_xxyyyy_yyyy[i] = tr_x_yyyy_yyyy[i] * fe_0 + ts_xyyyy_yyyy[i] * fe_0 + tr_x_xyyyy_yyyy[i] * pa_x[i];

        tr_x_xxyyyy_yyyz[i] = tr_x_yyyy_yyyz[i] * fe_0 + ts_xyyyy_yyyz[i] * fe_0 + tr_x_xyyyy_yyyz[i] * pa_x[i];

        tr_x_xxyyyy_yyzz[i] = tr_x_yyyy_yyzz[i] * fe_0 + ts_xyyyy_yyzz[i] * fe_0 + tr_x_xyyyy_yyzz[i] * pa_x[i];

        tr_x_xxyyyy_yzzz[i] = tr_x_yyyy_yzzz[i] * fe_0 + ts_xyyyy_yzzz[i] * fe_0 + tr_x_xyyyy_yzzz[i] * pa_x[i];

        tr_x_xxyyyy_zzzz[i] = 3.0 * tr_x_xxyy_zzzz[i] * fe_0 + tr_x_xxyyy_zzzz[i] * pa_y[i];
    }

    // Set up 165-180 components of targeted buffer : IG

    auto tr_x_xxyyyz_xxxx = pbuffer.data(idx_dip_ig + 165);

    auto tr_x_xxyyyz_xxxy = pbuffer.data(idx_dip_ig + 166);

    auto tr_x_xxyyyz_xxxz = pbuffer.data(idx_dip_ig + 167);

    auto tr_x_xxyyyz_xxyy = pbuffer.data(idx_dip_ig + 168);

    auto tr_x_xxyyyz_xxyz = pbuffer.data(idx_dip_ig + 169);

    auto tr_x_xxyyyz_xxzz = pbuffer.data(idx_dip_ig + 170);

    auto tr_x_xxyyyz_xyyy = pbuffer.data(idx_dip_ig + 171);

    auto tr_x_xxyyyz_xyyz = pbuffer.data(idx_dip_ig + 172);

    auto tr_x_xxyyyz_xyzz = pbuffer.data(idx_dip_ig + 173);

    auto tr_x_xxyyyz_xzzz = pbuffer.data(idx_dip_ig + 174);

    auto tr_x_xxyyyz_yyyy = pbuffer.data(idx_dip_ig + 175);

    auto tr_x_xxyyyz_yyyz = pbuffer.data(idx_dip_ig + 176);

    auto tr_x_xxyyyz_yyzz = pbuffer.data(idx_dip_ig + 177);

    auto tr_x_xxyyyz_yzzz = pbuffer.data(idx_dip_ig + 178);

    auto tr_x_xxyyyz_zzzz = pbuffer.data(idx_dip_ig + 179);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_xxyyy_xxxx,  \
                             tr_x_xxyyy_xxxy,  \
                             tr_x_xxyyy_xxy,   \
                             tr_x_xxyyy_xxyy,  \
                             tr_x_xxyyy_xxyz,  \
                             tr_x_xxyyy_xyy,   \
                             tr_x_xxyyy_xyyy,  \
                             tr_x_xxyyy_xyyz,  \
                             tr_x_xxyyy_xyz,   \
                             tr_x_xxyyy_xyzz,  \
                             tr_x_xxyyy_yyy,   \
                             tr_x_xxyyy_yyyy,  \
                             tr_x_xxyyy_yyyz,  \
                             tr_x_xxyyy_yyz,   \
                             tr_x_xxyyy_yyzz,  \
                             tr_x_xxyyy_yzz,   \
                             tr_x_xxyyy_yzzz,  \
                             tr_x_xxyyyz_xxxx, \
                             tr_x_xxyyyz_xxxy, \
                             tr_x_xxyyyz_xxxz, \
                             tr_x_xxyyyz_xxyy, \
                             tr_x_xxyyyz_xxyz, \
                             tr_x_xxyyyz_xxzz, \
                             tr_x_xxyyyz_xyyy, \
                             tr_x_xxyyyz_xyyz, \
                             tr_x_xxyyyz_xyzz, \
                             tr_x_xxyyyz_xzzz, \
                             tr_x_xxyyyz_yyyy, \
                             tr_x_xxyyyz_yyyz, \
                             tr_x_xxyyyz_yyzz, \
                             tr_x_xxyyyz_yzzz, \
                             tr_x_xxyyyz_zzzz, \
                             tr_x_xxyyz_xxxz,  \
                             tr_x_xxyyz_xxzz,  \
                             tr_x_xxyyz_xzzz,  \
                             tr_x_xxyyz_zzzz,  \
                             tr_x_xxyz_xxxz,   \
                             tr_x_xxyz_xxzz,   \
                             tr_x_xxyz_xzzz,   \
                             tr_x_xxyz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyz_xxxx[i] = tr_x_xxyyy_xxxx[i] * pa_z[i];

        tr_x_xxyyyz_xxxy[i] = tr_x_xxyyy_xxxy[i] * pa_z[i];

        tr_x_xxyyyz_xxxz[i] = 2.0 * tr_x_xxyz_xxxz[i] * fe_0 + tr_x_xxyyz_xxxz[i] * pa_y[i];

        tr_x_xxyyyz_xxyy[i] = tr_x_xxyyy_xxyy[i] * pa_z[i];

        tr_x_xxyyyz_xxyz[i] = tr_x_xxyyy_xxy[i] * fe_0 + tr_x_xxyyy_xxyz[i] * pa_z[i];

        tr_x_xxyyyz_xxzz[i] = 2.0 * tr_x_xxyz_xxzz[i] * fe_0 + tr_x_xxyyz_xxzz[i] * pa_y[i];

        tr_x_xxyyyz_xyyy[i] = tr_x_xxyyy_xyyy[i] * pa_z[i];

        tr_x_xxyyyz_xyyz[i] = tr_x_xxyyy_xyy[i] * fe_0 + tr_x_xxyyy_xyyz[i] * pa_z[i];

        tr_x_xxyyyz_xyzz[i] = 2.0 * tr_x_xxyyy_xyz[i] * fe_0 + tr_x_xxyyy_xyzz[i] * pa_z[i];

        tr_x_xxyyyz_xzzz[i] = 2.0 * tr_x_xxyz_xzzz[i] * fe_0 + tr_x_xxyyz_xzzz[i] * pa_y[i];

        tr_x_xxyyyz_yyyy[i] = tr_x_xxyyy_yyyy[i] * pa_z[i];

        tr_x_xxyyyz_yyyz[i] = tr_x_xxyyy_yyy[i] * fe_0 + tr_x_xxyyy_yyyz[i] * pa_z[i];

        tr_x_xxyyyz_yyzz[i] = 2.0 * tr_x_xxyyy_yyz[i] * fe_0 + tr_x_xxyyy_yyzz[i] * pa_z[i];

        tr_x_xxyyyz_yzzz[i] = 3.0 * tr_x_xxyyy_yzz[i] * fe_0 + tr_x_xxyyy_yzzz[i] * pa_z[i];

        tr_x_xxyyyz_zzzz[i] = 2.0 * tr_x_xxyz_zzzz[i] * fe_0 + tr_x_xxyyz_zzzz[i] * pa_y[i];
    }

    // Set up 180-195 components of targeted buffer : IG

    auto tr_x_xxyyzz_xxxx = pbuffer.data(idx_dip_ig + 180);

    auto tr_x_xxyyzz_xxxy = pbuffer.data(idx_dip_ig + 181);

    auto tr_x_xxyyzz_xxxz = pbuffer.data(idx_dip_ig + 182);

    auto tr_x_xxyyzz_xxyy = pbuffer.data(idx_dip_ig + 183);

    auto tr_x_xxyyzz_xxyz = pbuffer.data(idx_dip_ig + 184);

    auto tr_x_xxyyzz_xxzz = pbuffer.data(idx_dip_ig + 185);

    auto tr_x_xxyyzz_xyyy = pbuffer.data(idx_dip_ig + 186);

    auto tr_x_xxyyzz_xyyz = pbuffer.data(idx_dip_ig + 187);

    auto tr_x_xxyyzz_xyzz = pbuffer.data(idx_dip_ig + 188);

    auto tr_x_xxyyzz_xzzz = pbuffer.data(idx_dip_ig + 189);

    auto tr_x_xxyyzz_yyyy = pbuffer.data(idx_dip_ig + 190);

    auto tr_x_xxyyzz_yyyz = pbuffer.data(idx_dip_ig + 191);

    auto tr_x_xxyyzz_yyzz = pbuffer.data(idx_dip_ig + 192);

    auto tr_x_xxyyzz_yzzz = pbuffer.data(idx_dip_ig + 193);

    auto tr_x_xxyyzz_zzzz = pbuffer.data(idx_dip_ig + 194);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_x_xxyy_xxxy,   \
                             tr_x_xxyy_xxyy,   \
                             tr_x_xxyy_xyyy,   \
                             tr_x_xxyy_yyyy,   \
                             tr_x_xxyyz_xxxy,  \
                             tr_x_xxyyz_xxyy,  \
                             tr_x_xxyyz_xyyy,  \
                             tr_x_xxyyz_yyyy,  \
                             tr_x_xxyyzz_xxxx, \
                             tr_x_xxyyzz_xxxy, \
                             tr_x_xxyyzz_xxxz, \
                             tr_x_xxyyzz_xxyy, \
                             tr_x_xxyyzz_xxyz, \
                             tr_x_xxyyzz_xxzz, \
                             tr_x_xxyyzz_xyyy, \
                             tr_x_xxyyzz_xyyz, \
                             tr_x_xxyyzz_xyzz, \
                             tr_x_xxyyzz_xzzz, \
                             tr_x_xxyyzz_yyyy, \
                             tr_x_xxyyzz_yyyz, \
                             tr_x_xxyyzz_yyzz, \
                             tr_x_xxyyzz_yzzz, \
                             tr_x_xxyyzz_zzzz, \
                             tr_x_xxyzz_xxxx,  \
                             tr_x_xxyzz_xxxz,  \
                             tr_x_xxyzz_xxyz,  \
                             tr_x_xxyzz_xxz,   \
                             tr_x_xxyzz_xxzz,  \
                             tr_x_xxyzz_xyyz,  \
                             tr_x_xxyzz_xyz,   \
                             tr_x_xxyzz_xyzz,  \
                             tr_x_xxyzz_xzz,   \
                             tr_x_xxyzz_xzzz,  \
                             tr_x_xxyzz_zzzz,  \
                             tr_x_xxzz_xxxx,   \
                             tr_x_xxzz_xxxz,   \
                             tr_x_xxzz_xxyz,   \
                             tr_x_xxzz_xxzz,   \
                             tr_x_xxzz_xyyz,   \
                             tr_x_xxzz_xyzz,   \
                             tr_x_xxzz_xzzz,   \
                             tr_x_xxzz_zzzz,   \
                             tr_x_xyyzz_yyyz,  \
                             tr_x_xyyzz_yyzz,  \
                             tr_x_xyyzz_yzzz,  \
                             tr_x_yyzz_yyyz,   \
                             tr_x_yyzz_yyzz,   \
                             tr_x_yyzz_yzzz,   \
                             ts_xyyzz_yyyz,    \
                             ts_xyyzz_yyzz,    \
                             ts_xyyzz_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyzz_xxxx[i] = tr_x_xxzz_xxxx[i] * fe_0 + tr_x_xxyzz_xxxx[i] * pa_y[i];

        tr_x_xxyyzz_xxxy[i] = tr_x_xxyy_xxxy[i] * fe_0 + tr_x_xxyyz_xxxy[i] * pa_z[i];

        tr_x_xxyyzz_xxxz[i] = tr_x_xxzz_xxxz[i] * fe_0 + tr_x_xxyzz_xxxz[i] * pa_y[i];

        tr_x_xxyyzz_xxyy[i] = tr_x_xxyy_xxyy[i] * fe_0 + tr_x_xxyyz_xxyy[i] * pa_z[i];

        tr_x_xxyyzz_xxyz[i] = tr_x_xxzz_xxyz[i] * fe_0 + tr_x_xxyzz_xxz[i] * fe_0 + tr_x_xxyzz_xxyz[i] * pa_y[i];

        tr_x_xxyyzz_xxzz[i] = tr_x_xxzz_xxzz[i] * fe_0 + tr_x_xxyzz_xxzz[i] * pa_y[i];

        tr_x_xxyyzz_xyyy[i] = tr_x_xxyy_xyyy[i] * fe_0 + tr_x_xxyyz_xyyy[i] * pa_z[i];

        tr_x_xxyyzz_xyyz[i] = tr_x_xxzz_xyyz[i] * fe_0 + 2.0 * tr_x_xxyzz_xyz[i] * fe_0 + tr_x_xxyzz_xyyz[i] * pa_y[i];

        tr_x_xxyyzz_xyzz[i] = tr_x_xxzz_xyzz[i] * fe_0 + tr_x_xxyzz_xzz[i] * fe_0 + tr_x_xxyzz_xyzz[i] * pa_y[i];

        tr_x_xxyyzz_xzzz[i] = tr_x_xxzz_xzzz[i] * fe_0 + tr_x_xxyzz_xzzz[i] * pa_y[i];

        tr_x_xxyyzz_yyyy[i] = tr_x_xxyy_yyyy[i] * fe_0 + tr_x_xxyyz_yyyy[i] * pa_z[i];

        tr_x_xxyyzz_yyyz[i] = tr_x_yyzz_yyyz[i] * fe_0 + ts_xyyzz_yyyz[i] * fe_0 + tr_x_xyyzz_yyyz[i] * pa_x[i];

        tr_x_xxyyzz_yyzz[i] = tr_x_yyzz_yyzz[i] * fe_0 + ts_xyyzz_yyzz[i] * fe_0 + tr_x_xyyzz_yyzz[i] * pa_x[i];

        tr_x_xxyyzz_yzzz[i] = tr_x_yyzz_yzzz[i] * fe_0 + ts_xyyzz_yzzz[i] * fe_0 + tr_x_xyyzz_yzzz[i] * pa_x[i];

        tr_x_xxyyzz_zzzz[i] = tr_x_xxzz_zzzz[i] * fe_0 + tr_x_xxyzz_zzzz[i] * pa_y[i];
    }

    // Set up 195-210 components of targeted buffer : IG

    auto tr_x_xxyzzz_xxxx = pbuffer.data(idx_dip_ig + 195);

    auto tr_x_xxyzzz_xxxy = pbuffer.data(idx_dip_ig + 196);

    auto tr_x_xxyzzz_xxxz = pbuffer.data(idx_dip_ig + 197);

    auto tr_x_xxyzzz_xxyy = pbuffer.data(idx_dip_ig + 198);

    auto tr_x_xxyzzz_xxyz = pbuffer.data(idx_dip_ig + 199);

    auto tr_x_xxyzzz_xxzz = pbuffer.data(idx_dip_ig + 200);

    auto tr_x_xxyzzz_xyyy = pbuffer.data(idx_dip_ig + 201);

    auto tr_x_xxyzzz_xyyz = pbuffer.data(idx_dip_ig + 202);

    auto tr_x_xxyzzz_xyzz = pbuffer.data(idx_dip_ig + 203);

    auto tr_x_xxyzzz_xzzz = pbuffer.data(idx_dip_ig + 204);

    auto tr_x_xxyzzz_yyyy = pbuffer.data(idx_dip_ig + 205);

    auto tr_x_xxyzzz_yyyz = pbuffer.data(idx_dip_ig + 206);

    auto tr_x_xxyzzz_yyzz = pbuffer.data(idx_dip_ig + 207);

    auto tr_x_xxyzzz_yzzz = pbuffer.data(idx_dip_ig + 208);

    auto tr_x_xxyzzz_zzzz = pbuffer.data(idx_dip_ig + 209);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_xxyzzz_xxxx, \
                             tr_x_xxyzzz_xxxy, \
                             tr_x_xxyzzz_xxxz, \
                             tr_x_xxyzzz_xxyy, \
                             tr_x_xxyzzz_xxyz, \
                             tr_x_xxyzzz_xxzz, \
                             tr_x_xxyzzz_xyyy, \
                             tr_x_xxyzzz_xyyz, \
                             tr_x_xxyzzz_xyzz, \
                             tr_x_xxyzzz_xzzz, \
                             tr_x_xxyzzz_yyyy, \
                             tr_x_xxyzzz_yyyz, \
                             tr_x_xxyzzz_yyzz, \
                             tr_x_xxyzzz_yzzz, \
                             tr_x_xxyzzz_zzzz, \
                             tr_x_xxzzz_xxx,   \
                             tr_x_xxzzz_xxxx,  \
                             tr_x_xxzzz_xxxy,  \
                             tr_x_xxzzz_xxxz,  \
                             tr_x_xxzzz_xxy,   \
                             tr_x_xxzzz_xxyy,  \
                             tr_x_xxzzz_xxyz,  \
                             tr_x_xxzzz_xxz,   \
                             tr_x_xxzzz_xxzz,  \
                             tr_x_xxzzz_xyy,   \
                             tr_x_xxzzz_xyyy,  \
                             tr_x_xxzzz_xyyz,  \
                             tr_x_xxzzz_xyz,   \
                             tr_x_xxzzz_xyzz,  \
                             tr_x_xxzzz_xzz,   \
                             tr_x_xxzzz_xzzz,  \
                             tr_x_xxzzz_yyy,   \
                             tr_x_xxzzz_yyyy,  \
                             tr_x_xxzzz_yyyz,  \
                             tr_x_xxzzz_yyz,   \
                             tr_x_xxzzz_yyzz,  \
                             tr_x_xxzzz_yzz,   \
                             tr_x_xxzzz_yzzz,  \
                             tr_x_xxzzz_zzz,   \
                             tr_x_xxzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzzz_xxxx[i] = tr_x_xxzzz_xxxx[i] * pa_y[i];

        tr_x_xxyzzz_xxxy[i] = tr_x_xxzzz_xxx[i] * fe_0 + tr_x_xxzzz_xxxy[i] * pa_y[i];

        tr_x_xxyzzz_xxxz[i] = tr_x_xxzzz_xxxz[i] * pa_y[i];

        tr_x_xxyzzz_xxyy[i] = 2.0 * tr_x_xxzzz_xxy[i] * fe_0 + tr_x_xxzzz_xxyy[i] * pa_y[i];

        tr_x_xxyzzz_xxyz[i] = tr_x_xxzzz_xxz[i] * fe_0 + tr_x_xxzzz_xxyz[i] * pa_y[i];

        tr_x_xxyzzz_xxzz[i] = tr_x_xxzzz_xxzz[i] * pa_y[i];

        tr_x_xxyzzz_xyyy[i] = 3.0 * tr_x_xxzzz_xyy[i] * fe_0 + tr_x_xxzzz_xyyy[i] * pa_y[i];

        tr_x_xxyzzz_xyyz[i] = 2.0 * tr_x_xxzzz_xyz[i] * fe_0 + tr_x_xxzzz_xyyz[i] * pa_y[i];

        tr_x_xxyzzz_xyzz[i] = tr_x_xxzzz_xzz[i] * fe_0 + tr_x_xxzzz_xyzz[i] * pa_y[i];

        tr_x_xxyzzz_xzzz[i] = tr_x_xxzzz_xzzz[i] * pa_y[i];

        tr_x_xxyzzz_yyyy[i] = 4.0 * tr_x_xxzzz_yyy[i] * fe_0 + tr_x_xxzzz_yyyy[i] * pa_y[i];

        tr_x_xxyzzz_yyyz[i] = 3.0 * tr_x_xxzzz_yyz[i] * fe_0 + tr_x_xxzzz_yyyz[i] * pa_y[i];

        tr_x_xxyzzz_yyzz[i] = 2.0 * tr_x_xxzzz_yzz[i] * fe_0 + tr_x_xxzzz_yyzz[i] * pa_y[i];

        tr_x_xxyzzz_yzzz[i] = tr_x_xxzzz_zzz[i] * fe_0 + tr_x_xxzzz_yzzz[i] * pa_y[i];

        tr_x_xxyzzz_zzzz[i] = tr_x_xxzzz_zzzz[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : IG

    auto tr_x_xxzzzz_xxxx = pbuffer.data(idx_dip_ig + 210);

    auto tr_x_xxzzzz_xxxy = pbuffer.data(idx_dip_ig + 211);

    auto tr_x_xxzzzz_xxxz = pbuffer.data(idx_dip_ig + 212);

    auto tr_x_xxzzzz_xxyy = pbuffer.data(idx_dip_ig + 213);

    auto tr_x_xxzzzz_xxyz = pbuffer.data(idx_dip_ig + 214);

    auto tr_x_xxzzzz_xxzz = pbuffer.data(idx_dip_ig + 215);

    auto tr_x_xxzzzz_xyyy = pbuffer.data(idx_dip_ig + 216);

    auto tr_x_xxzzzz_xyyz = pbuffer.data(idx_dip_ig + 217);

    auto tr_x_xxzzzz_xyzz = pbuffer.data(idx_dip_ig + 218);

    auto tr_x_xxzzzz_xzzz = pbuffer.data(idx_dip_ig + 219);

    auto tr_x_xxzzzz_yyyy = pbuffer.data(idx_dip_ig + 220);

    auto tr_x_xxzzzz_yyyz = pbuffer.data(idx_dip_ig + 221);

    auto tr_x_xxzzzz_yyzz = pbuffer.data(idx_dip_ig + 222);

    auto tr_x_xxzzzz_yzzz = pbuffer.data(idx_dip_ig + 223);

    auto tr_x_xxzzzz_zzzz = pbuffer.data(idx_dip_ig + 224);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xxzz_xxxx,   \
                             tr_x_xxzz_xxxy,   \
                             tr_x_xxzz_xxxz,   \
                             tr_x_xxzz_xxyy,   \
                             tr_x_xxzz_xxyz,   \
                             tr_x_xxzz_xxzz,   \
                             tr_x_xxzz_xyyy,   \
                             tr_x_xxzz_xyyz,   \
                             tr_x_xxzz_xyzz,   \
                             tr_x_xxzz_xzzz,   \
                             tr_x_xxzz_yyyy,   \
                             tr_x_xxzzz_xxx,   \
                             tr_x_xxzzz_xxxx,  \
                             tr_x_xxzzz_xxxy,  \
                             tr_x_xxzzz_xxxz,  \
                             tr_x_xxzzz_xxy,   \
                             tr_x_xxzzz_xxyy,  \
                             tr_x_xxzzz_xxyz,  \
                             tr_x_xxzzz_xxz,   \
                             tr_x_xxzzz_xxzz,  \
                             tr_x_xxzzz_xyy,   \
                             tr_x_xxzzz_xyyy,  \
                             tr_x_xxzzz_xyyz,  \
                             tr_x_xxzzz_xyz,   \
                             tr_x_xxzzz_xyzz,  \
                             tr_x_xxzzz_xzz,   \
                             tr_x_xxzzz_xzzz,  \
                             tr_x_xxzzz_yyyy,  \
                             tr_x_xxzzzz_xxxx, \
                             tr_x_xxzzzz_xxxy, \
                             tr_x_xxzzzz_xxxz, \
                             tr_x_xxzzzz_xxyy, \
                             tr_x_xxzzzz_xxyz, \
                             tr_x_xxzzzz_xxzz, \
                             tr_x_xxzzzz_xyyy, \
                             tr_x_xxzzzz_xyyz, \
                             tr_x_xxzzzz_xyzz, \
                             tr_x_xxzzzz_xzzz, \
                             tr_x_xxzzzz_yyyy, \
                             tr_x_xxzzzz_yyyz, \
                             tr_x_xxzzzz_yyzz, \
                             tr_x_xxzzzz_yzzz, \
                             tr_x_xxzzzz_zzzz, \
                             tr_x_xzzzz_yyyz,  \
                             tr_x_xzzzz_yyzz,  \
                             tr_x_xzzzz_yzzz,  \
                             tr_x_xzzzz_zzzz,  \
                             tr_x_zzzz_yyyz,   \
                             tr_x_zzzz_yyzz,   \
                             tr_x_zzzz_yzzz,   \
                             tr_x_zzzz_zzzz,   \
                             ts_xzzzz_yyyz,    \
                             ts_xzzzz_yyzz,    \
                             ts_xzzzz_yzzz,    \
                             ts_xzzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzzz_xxxx[i] = 3.0 * tr_x_xxzz_xxxx[i] * fe_0 + tr_x_xxzzz_xxxx[i] * pa_z[i];

        tr_x_xxzzzz_xxxy[i] = 3.0 * tr_x_xxzz_xxxy[i] * fe_0 + tr_x_xxzzz_xxxy[i] * pa_z[i];

        tr_x_xxzzzz_xxxz[i] = 3.0 * tr_x_xxzz_xxxz[i] * fe_0 + tr_x_xxzzz_xxx[i] * fe_0 + tr_x_xxzzz_xxxz[i] * pa_z[i];

        tr_x_xxzzzz_xxyy[i] = 3.0 * tr_x_xxzz_xxyy[i] * fe_0 + tr_x_xxzzz_xxyy[i] * pa_z[i];

        tr_x_xxzzzz_xxyz[i] = 3.0 * tr_x_xxzz_xxyz[i] * fe_0 + tr_x_xxzzz_xxy[i] * fe_0 + tr_x_xxzzz_xxyz[i] * pa_z[i];

        tr_x_xxzzzz_xxzz[i] = 3.0 * tr_x_xxzz_xxzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xxz[i] * fe_0 + tr_x_xxzzz_xxzz[i] * pa_z[i];

        tr_x_xxzzzz_xyyy[i] = 3.0 * tr_x_xxzz_xyyy[i] * fe_0 + tr_x_xxzzz_xyyy[i] * pa_z[i];

        tr_x_xxzzzz_xyyz[i] = 3.0 * tr_x_xxzz_xyyz[i] * fe_0 + tr_x_xxzzz_xyy[i] * fe_0 + tr_x_xxzzz_xyyz[i] * pa_z[i];

        tr_x_xxzzzz_xyzz[i] = 3.0 * tr_x_xxzz_xyzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xyz[i] * fe_0 + tr_x_xxzzz_xyzz[i] * pa_z[i];

        tr_x_xxzzzz_xzzz[i] = 3.0 * tr_x_xxzz_xzzz[i] * fe_0 + 3.0 * tr_x_xxzzz_xzz[i] * fe_0 + tr_x_xxzzz_xzzz[i] * pa_z[i];

        tr_x_xxzzzz_yyyy[i] = 3.0 * tr_x_xxzz_yyyy[i] * fe_0 + tr_x_xxzzz_yyyy[i] * pa_z[i];

        tr_x_xxzzzz_yyyz[i] = tr_x_zzzz_yyyz[i] * fe_0 + ts_xzzzz_yyyz[i] * fe_0 + tr_x_xzzzz_yyyz[i] * pa_x[i];

        tr_x_xxzzzz_yyzz[i] = tr_x_zzzz_yyzz[i] * fe_0 + ts_xzzzz_yyzz[i] * fe_0 + tr_x_xzzzz_yyzz[i] * pa_x[i];

        tr_x_xxzzzz_yzzz[i] = tr_x_zzzz_yzzz[i] * fe_0 + ts_xzzzz_yzzz[i] * fe_0 + tr_x_xzzzz_yzzz[i] * pa_x[i];

        tr_x_xxzzzz_zzzz[i] = tr_x_zzzz_zzzz[i] * fe_0 + ts_xzzzz_zzzz[i] * fe_0 + tr_x_xzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : IG

    auto tr_x_xyyyyy_xxxx = pbuffer.data(idx_dip_ig + 225);

    auto tr_x_xyyyyy_xxxy = pbuffer.data(idx_dip_ig + 226);

    auto tr_x_xyyyyy_xxxz = pbuffer.data(idx_dip_ig + 227);

    auto tr_x_xyyyyy_xxyy = pbuffer.data(idx_dip_ig + 228);

    auto tr_x_xyyyyy_xxyz = pbuffer.data(idx_dip_ig + 229);

    auto tr_x_xyyyyy_xxzz = pbuffer.data(idx_dip_ig + 230);

    auto tr_x_xyyyyy_xyyy = pbuffer.data(idx_dip_ig + 231);

    auto tr_x_xyyyyy_xyyz = pbuffer.data(idx_dip_ig + 232);

    auto tr_x_xyyyyy_xyzz = pbuffer.data(idx_dip_ig + 233);

    auto tr_x_xyyyyy_xzzz = pbuffer.data(idx_dip_ig + 234);

    auto tr_x_xyyyyy_yyyy = pbuffer.data(idx_dip_ig + 235);

    auto tr_x_xyyyyy_yyyz = pbuffer.data(idx_dip_ig + 236);

    auto tr_x_xyyyyy_yyzz = pbuffer.data(idx_dip_ig + 237);

    auto tr_x_xyyyyy_yzzz = pbuffer.data(idx_dip_ig + 238);

    auto tr_x_xyyyyy_zzzz = pbuffer.data(idx_dip_ig + 239);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xyyy_xxxx,   \
                             tr_x_xyyy_xxxz,   \
                             tr_x_xyyy_xxzz,   \
                             tr_x_xyyy_xzzz,   \
                             tr_x_xyyyy_xxxx,  \
                             tr_x_xyyyy_xxxz,  \
                             tr_x_xyyyy_xxzz,  \
                             tr_x_xyyyy_xzzz,  \
                             tr_x_xyyyyy_xxxx, \
                             tr_x_xyyyyy_xxxy, \
                             tr_x_xyyyyy_xxxz, \
                             tr_x_xyyyyy_xxyy, \
                             tr_x_xyyyyy_xxyz, \
                             tr_x_xyyyyy_xxzz, \
                             tr_x_xyyyyy_xyyy, \
                             tr_x_xyyyyy_xyyz, \
                             tr_x_xyyyyy_xyzz, \
                             tr_x_xyyyyy_xzzz, \
                             tr_x_xyyyyy_yyyy, \
                             tr_x_xyyyyy_yyyz, \
                             tr_x_xyyyyy_yyzz, \
                             tr_x_xyyyyy_yzzz, \
                             tr_x_xyyyyy_zzzz, \
                             tr_x_yyyyy_xxxy,  \
                             tr_x_yyyyy_xxy,   \
                             tr_x_yyyyy_xxyy,  \
                             tr_x_yyyyy_xxyz,  \
                             tr_x_yyyyy_xyy,   \
                             tr_x_yyyyy_xyyy,  \
                             tr_x_yyyyy_xyyz,  \
                             tr_x_yyyyy_xyz,   \
                             tr_x_yyyyy_xyzz,  \
                             tr_x_yyyyy_yyy,   \
                             tr_x_yyyyy_yyyy,  \
                             tr_x_yyyyy_yyyz,  \
                             tr_x_yyyyy_yyz,   \
                             tr_x_yyyyy_yyzz,  \
                             tr_x_yyyyy_yzz,   \
                             tr_x_yyyyy_yzzz,  \
                             tr_x_yyyyy_zzzz,  \
                             ts_yyyyy_xxxy,    \
                             ts_yyyyy_xxyy,    \
                             ts_yyyyy_xxyz,    \
                             ts_yyyyy_xyyy,    \
                             ts_yyyyy_xyyz,    \
                             ts_yyyyy_xyzz,    \
                             ts_yyyyy_yyyy,    \
                             ts_yyyyy_yyyz,    \
                             ts_yyyyy_yyzz,    \
                             ts_yyyyy_yzzz,    \
                             ts_yyyyy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyy_xxxx[i] = 4.0 * tr_x_xyyy_xxxx[i] * fe_0 + tr_x_xyyyy_xxxx[i] * pa_y[i];

        tr_x_xyyyyy_xxxy[i] = 3.0 * tr_x_yyyyy_xxy[i] * fe_0 + ts_yyyyy_xxxy[i] * fe_0 + tr_x_yyyyy_xxxy[i] * pa_x[i];

        tr_x_xyyyyy_xxxz[i] = 4.0 * tr_x_xyyy_xxxz[i] * fe_0 + tr_x_xyyyy_xxxz[i] * pa_y[i];

        tr_x_xyyyyy_xxyy[i] = 2.0 * tr_x_yyyyy_xyy[i] * fe_0 + ts_yyyyy_xxyy[i] * fe_0 + tr_x_yyyyy_xxyy[i] * pa_x[i];

        tr_x_xyyyyy_xxyz[i] = 2.0 * tr_x_yyyyy_xyz[i] * fe_0 + ts_yyyyy_xxyz[i] * fe_0 + tr_x_yyyyy_xxyz[i] * pa_x[i];

        tr_x_xyyyyy_xxzz[i] = 4.0 * tr_x_xyyy_xxzz[i] * fe_0 + tr_x_xyyyy_xxzz[i] * pa_y[i];

        tr_x_xyyyyy_xyyy[i] = tr_x_yyyyy_yyy[i] * fe_0 + ts_yyyyy_xyyy[i] * fe_0 + tr_x_yyyyy_xyyy[i] * pa_x[i];

        tr_x_xyyyyy_xyyz[i] = tr_x_yyyyy_yyz[i] * fe_0 + ts_yyyyy_xyyz[i] * fe_0 + tr_x_yyyyy_xyyz[i] * pa_x[i];

        tr_x_xyyyyy_xyzz[i] = tr_x_yyyyy_yzz[i] * fe_0 + ts_yyyyy_xyzz[i] * fe_0 + tr_x_yyyyy_xyzz[i] * pa_x[i];

        tr_x_xyyyyy_xzzz[i] = 4.0 * tr_x_xyyy_xzzz[i] * fe_0 + tr_x_xyyyy_xzzz[i] * pa_y[i];

        tr_x_xyyyyy_yyyy[i] = ts_yyyyy_yyyy[i] * fe_0 + tr_x_yyyyy_yyyy[i] * pa_x[i];

        tr_x_xyyyyy_yyyz[i] = ts_yyyyy_yyyz[i] * fe_0 + tr_x_yyyyy_yyyz[i] * pa_x[i];

        tr_x_xyyyyy_yyzz[i] = ts_yyyyy_yyzz[i] * fe_0 + tr_x_yyyyy_yyzz[i] * pa_x[i];

        tr_x_xyyyyy_yzzz[i] = ts_yyyyy_yzzz[i] * fe_0 + tr_x_yyyyy_yzzz[i] * pa_x[i];

        tr_x_xyyyyy_zzzz[i] = ts_yyyyy_zzzz[i] * fe_0 + tr_x_yyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 240-255 components of targeted buffer : IG

    auto tr_x_xyyyyz_xxxx = pbuffer.data(idx_dip_ig + 240);

    auto tr_x_xyyyyz_xxxy = pbuffer.data(idx_dip_ig + 241);

    auto tr_x_xyyyyz_xxxz = pbuffer.data(idx_dip_ig + 242);

    auto tr_x_xyyyyz_xxyy = pbuffer.data(idx_dip_ig + 243);

    auto tr_x_xyyyyz_xxyz = pbuffer.data(idx_dip_ig + 244);

    auto tr_x_xyyyyz_xxzz = pbuffer.data(idx_dip_ig + 245);

    auto tr_x_xyyyyz_xyyy = pbuffer.data(idx_dip_ig + 246);

    auto tr_x_xyyyyz_xyyz = pbuffer.data(idx_dip_ig + 247);

    auto tr_x_xyyyyz_xyzz = pbuffer.data(idx_dip_ig + 248);

    auto tr_x_xyyyyz_xzzz = pbuffer.data(idx_dip_ig + 249);

    auto tr_x_xyyyyz_yyyy = pbuffer.data(idx_dip_ig + 250);

    auto tr_x_xyyyyz_yyyz = pbuffer.data(idx_dip_ig + 251);

    auto tr_x_xyyyyz_yyzz = pbuffer.data(idx_dip_ig + 252);

    auto tr_x_xyyyyz_yzzz = pbuffer.data(idx_dip_ig + 253);

    auto tr_x_xyyyyz_zzzz = pbuffer.data(idx_dip_ig + 254);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_x_xyyyy_xxxx,  \
                             tr_x_xyyyy_xxxy,  \
                             tr_x_xyyyy_xxy,   \
                             tr_x_xyyyy_xxyy,  \
                             tr_x_xyyyy_xxyz,  \
                             tr_x_xyyyy_xyy,   \
                             tr_x_xyyyy_xyyy,  \
                             tr_x_xyyyy_xyyz,  \
                             tr_x_xyyyy_xyz,   \
                             tr_x_xyyyy_xyzz,  \
                             tr_x_xyyyy_yyyy,  \
                             tr_x_xyyyyz_xxxx, \
                             tr_x_xyyyyz_xxxy, \
                             tr_x_xyyyyz_xxxz, \
                             tr_x_xyyyyz_xxyy, \
                             tr_x_xyyyyz_xxyz, \
                             tr_x_xyyyyz_xxzz, \
                             tr_x_xyyyyz_xyyy, \
                             tr_x_xyyyyz_xyyz, \
                             tr_x_xyyyyz_xyzz, \
                             tr_x_xyyyyz_xzzz, \
                             tr_x_xyyyyz_yyyy, \
                             tr_x_xyyyyz_yyyz, \
                             tr_x_xyyyyz_yyzz, \
                             tr_x_xyyyyz_yzzz, \
                             tr_x_xyyyyz_zzzz, \
                             tr_x_xyyyz_xxxz,  \
                             tr_x_xyyyz_xxzz,  \
                             tr_x_xyyyz_xzzz,  \
                             tr_x_xyyz_xxxz,   \
                             tr_x_xyyz_xxzz,   \
                             tr_x_xyyz_xzzz,   \
                             tr_x_yyyyz_yyyz,  \
                             tr_x_yyyyz_yyzz,  \
                             tr_x_yyyyz_yzzz,  \
                             tr_x_yyyyz_zzzz,  \
                             ts_yyyyz_yyyz,    \
                             ts_yyyyz_yyzz,    \
                             ts_yyyyz_yzzz,    \
                             ts_yyyyz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyz_xxxx[i] = tr_x_xyyyy_xxxx[i] * pa_z[i];

        tr_x_xyyyyz_xxxy[i] = tr_x_xyyyy_xxxy[i] * pa_z[i];

        tr_x_xyyyyz_xxxz[i] = 3.0 * tr_x_xyyz_xxxz[i] * fe_0 + tr_x_xyyyz_xxxz[i] * pa_y[i];

        tr_x_xyyyyz_xxyy[i] = tr_x_xyyyy_xxyy[i] * pa_z[i];

        tr_x_xyyyyz_xxyz[i] = tr_x_xyyyy_xxy[i] * fe_0 + tr_x_xyyyy_xxyz[i] * pa_z[i];

        tr_x_xyyyyz_xxzz[i] = 3.0 * tr_x_xyyz_xxzz[i] * fe_0 + tr_x_xyyyz_xxzz[i] * pa_y[i];

        tr_x_xyyyyz_xyyy[i] = tr_x_xyyyy_xyyy[i] * pa_z[i];

        tr_x_xyyyyz_xyyz[i] = tr_x_xyyyy_xyy[i] * fe_0 + tr_x_xyyyy_xyyz[i] * pa_z[i];

        tr_x_xyyyyz_xyzz[i] = 2.0 * tr_x_xyyyy_xyz[i] * fe_0 + tr_x_xyyyy_xyzz[i] * pa_z[i];

        tr_x_xyyyyz_xzzz[i] = 3.0 * tr_x_xyyz_xzzz[i] * fe_0 + tr_x_xyyyz_xzzz[i] * pa_y[i];

        tr_x_xyyyyz_yyyy[i] = tr_x_xyyyy_yyyy[i] * pa_z[i];

        tr_x_xyyyyz_yyyz[i] = ts_yyyyz_yyyz[i] * fe_0 + tr_x_yyyyz_yyyz[i] * pa_x[i];

        tr_x_xyyyyz_yyzz[i] = ts_yyyyz_yyzz[i] * fe_0 + tr_x_yyyyz_yyzz[i] * pa_x[i];

        tr_x_xyyyyz_yzzz[i] = ts_yyyyz_yzzz[i] * fe_0 + tr_x_yyyyz_yzzz[i] * pa_x[i];

        tr_x_xyyyyz_zzzz[i] = ts_yyyyz_zzzz[i] * fe_0 + tr_x_yyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 255-270 components of targeted buffer : IG

    auto tr_x_xyyyzz_xxxx = pbuffer.data(idx_dip_ig + 255);

    auto tr_x_xyyyzz_xxxy = pbuffer.data(idx_dip_ig + 256);

    auto tr_x_xyyyzz_xxxz = pbuffer.data(idx_dip_ig + 257);

    auto tr_x_xyyyzz_xxyy = pbuffer.data(idx_dip_ig + 258);

    auto tr_x_xyyyzz_xxyz = pbuffer.data(idx_dip_ig + 259);

    auto tr_x_xyyyzz_xxzz = pbuffer.data(idx_dip_ig + 260);

    auto tr_x_xyyyzz_xyyy = pbuffer.data(idx_dip_ig + 261);

    auto tr_x_xyyyzz_xyyz = pbuffer.data(idx_dip_ig + 262);

    auto tr_x_xyyyzz_xyzz = pbuffer.data(idx_dip_ig + 263);

    auto tr_x_xyyyzz_xzzz = pbuffer.data(idx_dip_ig + 264);

    auto tr_x_xyyyzz_yyyy = pbuffer.data(idx_dip_ig + 265);

    auto tr_x_xyyyzz_yyyz = pbuffer.data(idx_dip_ig + 266);

    auto tr_x_xyyyzz_yyzz = pbuffer.data(idx_dip_ig + 267);

    auto tr_x_xyyyzz_yzzz = pbuffer.data(idx_dip_ig + 268);

    auto tr_x_xyyyzz_zzzz = pbuffer.data(idx_dip_ig + 269);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_x_xyyy_xxxy,   \
                             tr_x_xyyy_xxyy,   \
                             tr_x_xyyy_xyyy,   \
                             tr_x_xyyyz_xxxy,  \
                             tr_x_xyyyz_xxyy,  \
                             tr_x_xyyyz_xyyy,  \
                             tr_x_xyyyzz_xxxx, \
                             tr_x_xyyyzz_xxxy, \
                             tr_x_xyyyzz_xxxz, \
                             tr_x_xyyyzz_xxyy, \
                             tr_x_xyyyzz_xxyz, \
                             tr_x_xyyyzz_xxzz, \
                             tr_x_xyyyzz_xyyy, \
                             tr_x_xyyyzz_xyyz, \
                             tr_x_xyyyzz_xyzz, \
                             tr_x_xyyyzz_xzzz, \
                             tr_x_xyyyzz_yyyy, \
                             tr_x_xyyyzz_yyyz, \
                             tr_x_xyyyzz_yyzz, \
                             tr_x_xyyyzz_yzzz, \
                             tr_x_xyyyzz_zzzz, \
                             tr_x_xyyzz_xxxx,  \
                             tr_x_xyyzz_xxxz,  \
                             tr_x_xyyzz_xxzz,  \
                             tr_x_xyyzz_xzzz,  \
                             tr_x_xyzz_xxxx,   \
                             tr_x_xyzz_xxxz,   \
                             tr_x_xyzz_xxzz,   \
                             tr_x_xyzz_xzzz,   \
                             tr_x_yyyzz_xxyz,  \
                             tr_x_yyyzz_xyyz,  \
                             tr_x_yyyzz_xyz,   \
                             tr_x_yyyzz_xyzz,  \
                             tr_x_yyyzz_yyyy,  \
                             tr_x_yyyzz_yyyz,  \
                             tr_x_yyyzz_yyz,   \
                             tr_x_yyyzz_yyzz,  \
                             tr_x_yyyzz_yzz,   \
                             tr_x_yyyzz_yzzz,  \
                             tr_x_yyyzz_zzzz,  \
                             ts_yyyzz_xxyz,    \
                             ts_yyyzz_xyyz,    \
                             ts_yyyzz_xyzz,    \
                             ts_yyyzz_yyyy,    \
                             ts_yyyzz_yyyz,    \
                             ts_yyyzz_yyzz,    \
                             ts_yyyzz_yzzz,    \
                             ts_yyyzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyzz_xxxx[i] = 2.0 * tr_x_xyzz_xxxx[i] * fe_0 + tr_x_xyyzz_xxxx[i] * pa_y[i];

        tr_x_xyyyzz_xxxy[i] = tr_x_xyyy_xxxy[i] * fe_0 + tr_x_xyyyz_xxxy[i] * pa_z[i];

        tr_x_xyyyzz_xxxz[i] = 2.0 * tr_x_xyzz_xxxz[i] * fe_0 + tr_x_xyyzz_xxxz[i] * pa_y[i];

        tr_x_xyyyzz_xxyy[i] = tr_x_xyyy_xxyy[i] * fe_0 + tr_x_xyyyz_xxyy[i] * pa_z[i];

        tr_x_xyyyzz_xxyz[i] = 2.0 * tr_x_yyyzz_xyz[i] * fe_0 + ts_yyyzz_xxyz[i] * fe_0 + tr_x_yyyzz_xxyz[i] * pa_x[i];

        tr_x_xyyyzz_xxzz[i] = 2.0 * tr_x_xyzz_xxzz[i] * fe_0 + tr_x_xyyzz_xxzz[i] * pa_y[i];

        tr_x_xyyyzz_xyyy[i] = tr_x_xyyy_xyyy[i] * fe_0 + tr_x_xyyyz_xyyy[i] * pa_z[i];

        tr_x_xyyyzz_xyyz[i] = tr_x_yyyzz_yyz[i] * fe_0 + ts_yyyzz_xyyz[i] * fe_0 + tr_x_yyyzz_xyyz[i] * pa_x[i];

        tr_x_xyyyzz_xyzz[i] = tr_x_yyyzz_yzz[i] * fe_0 + ts_yyyzz_xyzz[i] * fe_0 + tr_x_yyyzz_xyzz[i] * pa_x[i];

        tr_x_xyyyzz_xzzz[i] = 2.0 * tr_x_xyzz_xzzz[i] * fe_0 + tr_x_xyyzz_xzzz[i] * pa_y[i];

        tr_x_xyyyzz_yyyy[i] = ts_yyyzz_yyyy[i] * fe_0 + tr_x_yyyzz_yyyy[i] * pa_x[i];

        tr_x_xyyyzz_yyyz[i] = ts_yyyzz_yyyz[i] * fe_0 + tr_x_yyyzz_yyyz[i] * pa_x[i];

        tr_x_xyyyzz_yyzz[i] = ts_yyyzz_yyzz[i] * fe_0 + tr_x_yyyzz_yyzz[i] * pa_x[i];

        tr_x_xyyyzz_yzzz[i] = ts_yyyzz_yzzz[i] * fe_0 + tr_x_yyyzz_yzzz[i] * pa_x[i];

        tr_x_xyyyzz_zzzz[i] = ts_yyyzz_zzzz[i] * fe_0 + tr_x_yyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 270-285 components of targeted buffer : IG

    auto tr_x_xyyzzz_xxxx = pbuffer.data(idx_dip_ig + 270);

    auto tr_x_xyyzzz_xxxy = pbuffer.data(idx_dip_ig + 271);

    auto tr_x_xyyzzz_xxxz = pbuffer.data(idx_dip_ig + 272);

    auto tr_x_xyyzzz_xxyy = pbuffer.data(idx_dip_ig + 273);

    auto tr_x_xyyzzz_xxyz = pbuffer.data(idx_dip_ig + 274);

    auto tr_x_xyyzzz_xxzz = pbuffer.data(idx_dip_ig + 275);

    auto tr_x_xyyzzz_xyyy = pbuffer.data(idx_dip_ig + 276);

    auto tr_x_xyyzzz_xyyz = pbuffer.data(idx_dip_ig + 277);

    auto tr_x_xyyzzz_xyzz = pbuffer.data(idx_dip_ig + 278);

    auto tr_x_xyyzzz_xzzz = pbuffer.data(idx_dip_ig + 279);

    auto tr_x_xyyzzz_yyyy = pbuffer.data(idx_dip_ig + 280);

    auto tr_x_xyyzzz_yyyz = pbuffer.data(idx_dip_ig + 281);

    auto tr_x_xyyzzz_yyzz = pbuffer.data(idx_dip_ig + 282);

    auto tr_x_xyyzzz_yzzz = pbuffer.data(idx_dip_ig + 283);

    auto tr_x_xyyzzz_zzzz = pbuffer.data(idx_dip_ig + 284);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_x_xyyz_xxxy,   \
                             tr_x_xyyz_xxyy,   \
                             tr_x_xyyz_xyyy,   \
                             tr_x_xyyzz_xxxy,  \
                             tr_x_xyyzz_xxyy,  \
                             tr_x_xyyzz_xyyy,  \
                             tr_x_xyyzzz_xxxx, \
                             tr_x_xyyzzz_xxxy, \
                             tr_x_xyyzzz_xxxz, \
                             tr_x_xyyzzz_xxyy, \
                             tr_x_xyyzzz_xxyz, \
                             tr_x_xyyzzz_xxzz, \
                             tr_x_xyyzzz_xyyy, \
                             tr_x_xyyzzz_xyyz, \
                             tr_x_xyyzzz_xyzz, \
                             tr_x_xyyzzz_xzzz, \
                             tr_x_xyyzzz_yyyy, \
                             tr_x_xyyzzz_yyyz, \
                             tr_x_xyyzzz_yyzz, \
                             tr_x_xyyzzz_yzzz, \
                             tr_x_xyyzzz_zzzz, \
                             tr_x_xyzzz_xxxx,  \
                             tr_x_xyzzz_xxxz,  \
                             tr_x_xyzzz_xxzz,  \
                             tr_x_xyzzz_xzzz,  \
                             tr_x_xzzz_xxxx,   \
                             tr_x_xzzz_xxxz,   \
                             tr_x_xzzz_xxzz,   \
                             tr_x_xzzz_xzzz,   \
                             tr_x_yyzzz_xxyz,  \
                             tr_x_yyzzz_xyyz,  \
                             tr_x_yyzzz_xyz,   \
                             tr_x_yyzzz_xyzz,  \
                             tr_x_yyzzz_yyyy,  \
                             tr_x_yyzzz_yyyz,  \
                             tr_x_yyzzz_yyz,   \
                             tr_x_yyzzz_yyzz,  \
                             tr_x_yyzzz_yzz,   \
                             tr_x_yyzzz_yzzz,  \
                             tr_x_yyzzz_zzzz,  \
                             ts_yyzzz_xxyz,    \
                             ts_yyzzz_xyyz,    \
                             ts_yyzzz_xyzz,    \
                             ts_yyzzz_yyyy,    \
                             ts_yyzzz_yyyz,    \
                             ts_yyzzz_yyzz,    \
                             ts_yyzzz_yzzz,    \
                             ts_yyzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzzz_xxxx[i] = tr_x_xzzz_xxxx[i] * fe_0 + tr_x_xyzzz_xxxx[i] * pa_y[i];

        tr_x_xyyzzz_xxxy[i] = 2.0 * tr_x_xyyz_xxxy[i] * fe_0 + tr_x_xyyzz_xxxy[i] * pa_z[i];

        tr_x_xyyzzz_xxxz[i] = tr_x_xzzz_xxxz[i] * fe_0 + tr_x_xyzzz_xxxz[i] * pa_y[i];

        tr_x_xyyzzz_xxyy[i] = 2.0 * tr_x_xyyz_xxyy[i] * fe_0 + tr_x_xyyzz_xxyy[i] * pa_z[i];

        tr_x_xyyzzz_xxyz[i] = 2.0 * tr_x_yyzzz_xyz[i] * fe_0 + ts_yyzzz_xxyz[i] * fe_0 + tr_x_yyzzz_xxyz[i] * pa_x[i];

        tr_x_xyyzzz_xxzz[i] = tr_x_xzzz_xxzz[i] * fe_0 + tr_x_xyzzz_xxzz[i] * pa_y[i];

        tr_x_xyyzzz_xyyy[i] = 2.0 * tr_x_xyyz_xyyy[i] * fe_0 + tr_x_xyyzz_xyyy[i] * pa_z[i];

        tr_x_xyyzzz_xyyz[i] = tr_x_yyzzz_yyz[i] * fe_0 + ts_yyzzz_xyyz[i] * fe_0 + tr_x_yyzzz_xyyz[i] * pa_x[i];

        tr_x_xyyzzz_xyzz[i] = tr_x_yyzzz_yzz[i] * fe_0 + ts_yyzzz_xyzz[i] * fe_0 + tr_x_yyzzz_xyzz[i] * pa_x[i];

        tr_x_xyyzzz_xzzz[i] = tr_x_xzzz_xzzz[i] * fe_0 + tr_x_xyzzz_xzzz[i] * pa_y[i];

        tr_x_xyyzzz_yyyy[i] = ts_yyzzz_yyyy[i] * fe_0 + tr_x_yyzzz_yyyy[i] * pa_x[i];

        tr_x_xyyzzz_yyyz[i] = ts_yyzzz_yyyz[i] * fe_0 + tr_x_yyzzz_yyyz[i] * pa_x[i];

        tr_x_xyyzzz_yyzz[i] = ts_yyzzz_yyzz[i] * fe_0 + tr_x_yyzzz_yyzz[i] * pa_x[i];

        tr_x_xyyzzz_yzzz[i] = ts_yyzzz_yzzz[i] * fe_0 + tr_x_yyzzz_yzzz[i] * pa_x[i];

        tr_x_xyyzzz_zzzz[i] = ts_yyzzz_zzzz[i] * fe_0 + tr_x_yyzzz_zzzz[i] * pa_x[i];
    }

    // Set up 285-300 components of targeted buffer : IG

    auto tr_x_xyzzzz_xxxx = pbuffer.data(idx_dip_ig + 285);

    auto tr_x_xyzzzz_xxxy = pbuffer.data(idx_dip_ig + 286);

    auto tr_x_xyzzzz_xxxz = pbuffer.data(idx_dip_ig + 287);

    auto tr_x_xyzzzz_xxyy = pbuffer.data(idx_dip_ig + 288);

    auto tr_x_xyzzzz_xxyz = pbuffer.data(idx_dip_ig + 289);

    auto tr_x_xyzzzz_xxzz = pbuffer.data(idx_dip_ig + 290);

    auto tr_x_xyzzzz_xyyy = pbuffer.data(idx_dip_ig + 291);

    auto tr_x_xyzzzz_xyyz = pbuffer.data(idx_dip_ig + 292);

    auto tr_x_xyzzzz_xyzz = pbuffer.data(idx_dip_ig + 293);

    auto tr_x_xyzzzz_xzzz = pbuffer.data(idx_dip_ig + 294);

    auto tr_x_xyzzzz_yyyy = pbuffer.data(idx_dip_ig + 295);

    auto tr_x_xyzzzz_yyyz = pbuffer.data(idx_dip_ig + 296);

    auto tr_x_xyzzzz_yyzz = pbuffer.data(idx_dip_ig + 297);

    auto tr_x_xyzzzz_yzzz = pbuffer.data(idx_dip_ig + 298);

    auto tr_x_xyzzzz_zzzz = pbuffer.data(idx_dip_ig + 299);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xyzzzz_xxxx, \
                             tr_x_xyzzzz_xxxy, \
                             tr_x_xyzzzz_xxxz, \
                             tr_x_xyzzzz_xxyy, \
                             tr_x_xyzzzz_xxyz, \
                             tr_x_xyzzzz_xxzz, \
                             tr_x_xyzzzz_xyyy, \
                             tr_x_xyzzzz_xyyz, \
                             tr_x_xyzzzz_xyzz, \
                             tr_x_xyzzzz_xzzz, \
                             tr_x_xyzzzz_yyyy, \
                             tr_x_xyzzzz_yyyz, \
                             tr_x_xyzzzz_yyzz, \
                             tr_x_xyzzzz_yzzz, \
                             tr_x_xyzzzz_zzzz, \
                             tr_x_xzzzz_xxx,   \
                             tr_x_xzzzz_xxxx,  \
                             tr_x_xzzzz_xxxy,  \
                             tr_x_xzzzz_xxxz,  \
                             tr_x_xzzzz_xxy,   \
                             tr_x_xzzzz_xxyy,  \
                             tr_x_xzzzz_xxyz,  \
                             tr_x_xzzzz_xxz,   \
                             tr_x_xzzzz_xxzz,  \
                             tr_x_xzzzz_xyy,   \
                             tr_x_xzzzz_xyyy,  \
                             tr_x_xzzzz_xyyz,  \
                             tr_x_xzzzz_xyz,   \
                             tr_x_xzzzz_xyzz,  \
                             tr_x_xzzzz_xzz,   \
                             tr_x_xzzzz_xzzz,  \
                             tr_x_xzzzz_zzzz,  \
                             tr_x_yzzzz_yyyy,  \
                             tr_x_yzzzz_yyyz,  \
                             tr_x_yzzzz_yyzz,  \
                             tr_x_yzzzz_yzzz,  \
                             ts_yzzzz_yyyy,    \
                             ts_yzzzz_yyyz,    \
                             ts_yzzzz_yyzz,    \
                             ts_yzzzz_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzzz_xxxx[i] = tr_x_xzzzz_xxxx[i] * pa_y[i];

        tr_x_xyzzzz_xxxy[i] = tr_x_xzzzz_xxx[i] * fe_0 + tr_x_xzzzz_xxxy[i] * pa_y[i];

        tr_x_xyzzzz_xxxz[i] = tr_x_xzzzz_xxxz[i] * pa_y[i];

        tr_x_xyzzzz_xxyy[i] = 2.0 * tr_x_xzzzz_xxy[i] * fe_0 + tr_x_xzzzz_xxyy[i] * pa_y[i];

        tr_x_xyzzzz_xxyz[i] = tr_x_xzzzz_xxz[i] * fe_0 + tr_x_xzzzz_xxyz[i] * pa_y[i];

        tr_x_xyzzzz_xxzz[i] = tr_x_xzzzz_xxzz[i] * pa_y[i];

        tr_x_xyzzzz_xyyy[i] = 3.0 * tr_x_xzzzz_xyy[i] * fe_0 + tr_x_xzzzz_xyyy[i] * pa_y[i];

        tr_x_xyzzzz_xyyz[i] = 2.0 * tr_x_xzzzz_xyz[i] * fe_0 + tr_x_xzzzz_xyyz[i] * pa_y[i];

        tr_x_xyzzzz_xyzz[i] = tr_x_xzzzz_xzz[i] * fe_0 + tr_x_xzzzz_xyzz[i] * pa_y[i];

        tr_x_xyzzzz_xzzz[i] = tr_x_xzzzz_xzzz[i] * pa_y[i];

        tr_x_xyzzzz_yyyy[i] = ts_yzzzz_yyyy[i] * fe_0 + tr_x_yzzzz_yyyy[i] * pa_x[i];

        tr_x_xyzzzz_yyyz[i] = ts_yzzzz_yyyz[i] * fe_0 + tr_x_yzzzz_yyyz[i] * pa_x[i];

        tr_x_xyzzzz_yyzz[i] = ts_yzzzz_yyzz[i] * fe_0 + tr_x_yzzzz_yyzz[i] * pa_x[i];

        tr_x_xyzzzz_yzzz[i] = ts_yzzzz_yzzz[i] * fe_0 + tr_x_yzzzz_yzzz[i] * pa_x[i];

        tr_x_xyzzzz_zzzz[i] = tr_x_xzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 300-315 components of targeted buffer : IG

    auto tr_x_xzzzzz_xxxx = pbuffer.data(idx_dip_ig + 300);

    auto tr_x_xzzzzz_xxxy = pbuffer.data(idx_dip_ig + 301);

    auto tr_x_xzzzzz_xxxz = pbuffer.data(idx_dip_ig + 302);

    auto tr_x_xzzzzz_xxyy = pbuffer.data(idx_dip_ig + 303);

    auto tr_x_xzzzzz_xxyz = pbuffer.data(idx_dip_ig + 304);

    auto tr_x_xzzzzz_xxzz = pbuffer.data(idx_dip_ig + 305);

    auto tr_x_xzzzzz_xyyy = pbuffer.data(idx_dip_ig + 306);

    auto tr_x_xzzzzz_xyyz = pbuffer.data(idx_dip_ig + 307);

    auto tr_x_xzzzzz_xyzz = pbuffer.data(idx_dip_ig + 308);

    auto tr_x_xzzzzz_xzzz = pbuffer.data(idx_dip_ig + 309);

    auto tr_x_xzzzzz_yyyy = pbuffer.data(idx_dip_ig + 310);

    auto tr_x_xzzzzz_yyyz = pbuffer.data(idx_dip_ig + 311);

    auto tr_x_xzzzzz_yyzz = pbuffer.data(idx_dip_ig + 312);

    auto tr_x_xzzzzz_yzzz = pbuffer.data(idx_dip_ig + 313);

    auto tr_x_xzzzzz_zzzz = pbuffer.data(idx_dip_ig + 314);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xzzz_xxxx,   \
                             tr_x_xzzz_xxxy,   \
                             tr_x_xzzz_xxyy,   \
                             tr_x_xzzz_xyyy,   \
                             tr_x_xzzzz_xxxx,  \
                             tr_x_xzzzz_xxxy,  \
                             tr_x_xzzzz_xxyy,  \
                             tr_x_xzzzz_xyyy,  \
                             tr_x_xzzzzz_xxxx, \
                             tr_x_xzzzzz_xxxy, \
                             tr_x_xzzzzz_xxxz, \
                             tr_x_xzzzzz_xxyy, \
                             tr_x_xzzzzz_xxyz, \
                             tr_x_xzzzzz_xxzz, \
                             tr_x_xzzzzz_xyyy, \
                             tr_x_xzzzzz_xyyz, \
                             tr_x_xzzzzz_xyzz, \
                             tr_x_xzzzzz_xzzz, \
                             tr_x_xzzzzz_yyyy, \
                             tr_x_xzzzzz_yyyz, \
                             tr_x_xzzzzz_yyzz, \
                             tr_x_xzzzzz_yzzz, \
                             tr_x_xzzzzz_zzzz, \
                             tr_x_zzzzz_xxxz,  \
                             tr_x_zzzzz_xxyz,  \
                             tr_x_zzzzz_xxz,   \
                             tr_x_zzzzz_xxzz,  \
                             tr_x_zzzzz_xyyz,  \
                             tr_x_zzzzz_xyz,   \
                             tr_x_zzzzz_xyzz,  \
                             tr_x_zzzzz_xzz,   \
                             tr_x_zzzzz_xzzz,  \
                             tr_x_zzzzz_yyyy,  \
                             tr_x_zzzzz_yyyz,  \
                             tr_x_zzzzz_yyz,   \
                             tr_x_zzzzz_yyzz,  \
                             tr_x_zzzzz_yzz,   \
                             tr_x_zzzzz_yzzz,  \
                             tr_x_zzzzz_zzz,   \
                             tr_x_zzzzz_zzzz,  \
                             ts_zzzzz_xxxz,    \
                             ts_zzzzz_xxyz,    \
                             ts_zzzzz_xxzz,    \
                             ts_zzzzz_xyyz,    \
                             ts_zzzzz_xyzz,    \
                             ts_zzzzz_xzzz,    \
                             ts_zzzzz_yyyy,    \
                             ts_zzzzz_yyyz,    \
                             ts_zzzzz_yyzz,    \
                             ts_zzzzz_yzzz,    \
                             ts_zzzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzzz_xxxx[i] = 4.0 * tr_x_xzzz_xxxx[i] * fe_0 + tr_x_xzzzz_xxxx[i] * pa_z[i];

        tr_x_xzzzzz_xxxy[i] = 4.0 * tr_x_xzzz_xxxy[i] * fe_0 + tr_x_xzzzz_xxxy[i] * pa_z[i];

        tr_x_xzzzzz_xxxz[i] = 3.0 * tr_x_zzzzz_xxz[i] * fe_0 + ts_zzzzz_xxxz[i] * fe_0 + tr_x_zzzzz_xxxz[i] * pa_x[i];

        tr_x_xzzzzz_xxyy[i] = 4.0 * tr_x_xzzz_xxyy[i] * fe_0 + tr_x_xzzzz_xxyy[i] * pa_z[i];

        tr_x_xzzzzz_xxyz[i] = 2.0 * tr_x_zzzzz_xyz[i] * fe_0 + ts_zzzzz_xxyz[i] * fe_0 + tr_x_zzzzz_xxyz[i] * pa_x[i];

        tr_x_xzzzzz_xxzz[i] = 2.0 * tr_x_zzzzz_xzz[i] * fe_0 + ts_zzzzz_xxzz[i] * fe_0 + tr_x_zzzzz_xxzz[i] * pa_x[i];

        tr_x_xzzzzz_xyyy[i] = 4.0 * tr_x_xzzz_xyyy[i] * fe_0 + tr_x_xzzzz_xyyy[i] * pa_z[i];

        tr_x_xzzzzz_xyyz[i] = tr_x_zzzzz_yyz[i] * fe_0 + ts_zzzzz_xyyz[i] * fe_0 + tr_x_zzzzz_xyyz[i] * pa_x[i];

        tr_x_xzzzzz_xyzz[i] = tr_x_zzzzz_yzz[i] * fe_0 + ts_zzzzz_xyzz[i] * fe_0 + tr_x_zzzzz_xyzz[i] * pa_x[i];

        tr_x_xzzzzz_xzzz[i] = tr_x_zzzzz_zzz[i] * fe_0 + ts_zzzzz_xzzz[i] * fe_0 + tr_x_zzzzz_xzzz[i] * pa_x[i];

        tr_x_xzzzzz_yyyy[i] = ts_zzzzz_yyyy[i] * fe_0 + tr_x_zzzzz_yyyy[i] * pa_x[i];

        tr_x_xzzzzz_yyyz[i] = ts_zzzzz_yyyz[i] * fe_0 + tr_x_zzzzz_yyyz[i] * pa_x[i];

        tr_x_xzzzzz_yyzz[i] = ts_zzzzz_yyzz[i] * fe_0 + tr_x_zzzzz_yyzz[i] * pa_x[i];

        tr_x_xzzzzz_yzzz[i] = ts_zzzzz_yzzz[i] * fe_0 + tr_x_zzzzz_yzzz[i] * pa_x[i];

        tr_x_xzzzzz_zzzz[i] = ts_zzzzz_zzzz[i] * fe_0 + tr_x_zzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 315-330 components of targeted buffer : IG

    auto tr_x_yyyyyy_xxxx = pbuffer.data(idx_dip_ig + 315);

    auto tr_x_yyyyyy_xxxy = pbuffer.data(idx_dip_ig + 316);

    auto tr_x_yyyyyy_xxxz = pbuffer.data(idx_dip_ig + 317);

    auto tr_x_yyyyyy_xxyy = pbuffer.data(idx_dip_ig + 318);

    auto tr_x_yyyyyy_xxyz = pbuffer.data(idx_dip_ig + 319);

    auto tr_x_yyyyyy_xxzz = pbuffer.data(idx_dip_ig + 320);

    auto tr_x_yyyyyy_xyyy = pbuffer.data(idx_dip_ig + 321);

    auto tr_x_yyyyyy_xyyz = pbuffer.data(idx_dip_ig + 322);

    auto tr_x_yyyyyy_xyzz = pbuffer.data(idx_dip_ig + 323);

    auto tr_x_yyyyyy_xzzz = pbuffer.data(idx_dip_ig + 324);

    auto tr_x_yyyyyy_yyyy = pbuffer.data(idx_dip_ig + 325);

    auto tr_x_yyyyyy_yyyz = pbuffer.data(idx_dip_ig + 326);

    auto tr_x_yyyyyy_yyzz = pbuffer.data(idx_dip_ig + 327);

    auto tr_x_yyyyyy_yzzz = pbuffer.data(idx_dip_ig + 328);

    auto tr_x_yyyyyy_zzzz = pbuffer.data(idx_dip_ig + 329);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_yyyy_xxxx,   \
                             tr_x_yyyy_xxxy,   \
                             tr_x_yyyy_xxxz,   \
                             tr_x_yyyy_xxyy,   \
                             tr_x_yyyy_xxyz,   \
                             tr_x_yyyy_xxzz,   \
                             tr_x_yyyy_xyyy,   \
                             tr_x_yyyy_xyyz,   \
                             tr_x_yyyy_xyzz,   \
                             tr_x_yyyy_xzzz,   \
                             tr_x_yyyy_yyyy,   \
                             tr_x_yyyy_yyyz,   \
                             tr_x_yyyy_yyzz,   \
                             tr_x_yyyy_yzzz,   \
                             tr_x_yyyy_zzzz,   \
                             tr_x_yyyyy_xxx,   \
                             tr_x_yyyyy_xxxx,  \
                             tr_x_yyyyy_xxxy,  \
                             tr_x_yyyyy_xxxz,  \
                             tr_x_yyyyy_xxy,   \
                             tr_x_yyyyy_xxyy,  \
                             tr_x_yyyyy_xxyz,  \
                             tr_x_yyyyy_xxz,   \
                             tr_x_yyyyy_xxzz,  \
                             tr_x_yyyyy_xyy,   \
                             tr_x_yyyyy_xyyy,  \
                             tr_x_yyyyy_xyyz,  \
                             tr_x_yyyyy_xyz,   \
                             tr_x_yyyyy_xyzz,  \
                             tr_x_yyyyy_xzz,   \
                             tr_x_yyyyy_xzzz,  \
                             tr_x_yyyyy_yyy,   \
                             tr_x_yyyyy_yyyy,  \
                             tr_x_yyyyy_yyyz,  \
                             tr_x_yyyyy_yyz,   \
                             tr_x_yyyyy_yyzz,  \
                             tr_x_yyyyy_yzz,   \
                             tr_x_yyyyy_yzzz,  \
                             tr_x_yyyyy_zzz,   \
                             tr_x_yyyyy_zzzz,  \
                             tr_x_yyyyyy_xxxx, \
                             tr_x_yyyyyy_xxxy, \
                             tr_x_yyyyyy_xxxz, \
                             tr_x_yyyyyy_xxyy, \
                             tr_x_yyyyyy_xxyz, \
                             tr_x_yyyyyy_xxzz, \
                             tr_x_yyyyyy_xyyy, \
                             tr_x_yyyyyy_xyyz, \
                             tr_x_yyyyyy_xyzz, \
                             tr_x_yyyyyy_xzzz, \
                             tr_x_yyyyyy_yyyy, \
                             tr_x_yyyyyy_yyyz, \
                             tr_x_yyyyyy_yyzz, \
                             tr_x_yyyyyy_yzzz, \
                             tr_x_yyyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyy_xxxx[i] = 5.0 * tr_x_yyyy_xxxx[i] * fe_0 + tr_x_yyyyy_xxxx[i] * pa_y[i];

        tr_x_yyyyyy_xxxy[i] = 5.0 * tr_x_yyyy_xxxy[i] * fe_0 + tr_x_yyyyy_xxx[i] * fe_0 + tr_x_yyyyy_xxxy[i] * pa_y[i];

        tr_x_yyyyyy_xxxz[i] = 5.0 * tr_x_yyyy_xxxz[i] * fe_0 + tr_x_yyyyy_xxxz[i] * pa_y[i];

        tr_x_yyyyyy_xxyy[i] = 5.0 * tr_x_yyyy_xxyy[i] * fe_0 + 2.0 * tr_x_yyyyy_xxy[i] * fe_0 + tr_x_yyyyy_xxyy[i] * pa_y[i];

        tr_x_yyyyyy_xxyz[i] = 5.0 * tr_x_yyyy_xxyz[i] * fe_0 + tr_x_yyyyy_xxz[i] * fe_0 + tr_x_yyyyy_xxyz[i] * pa_y[i];

        tr_x_yyyyyy_xxzz[i] = 5.0 * tr_x_yyyy_xxzz[i] * fe_0 + tr_x_yyyyy_xxzz[i] * pa_y[i];

        tr_x_yyyyyy_xyyy[i] = 5.0 * tr_x_yyyy_xyyy[i] * fe_0 + 3.0 * tr_x_yyyyy_xyy[i] * fe_0 + tr_x_yyyyy_xyyy[i] * pa_y[i];

        tr_x_yyyyyy_xyyz[i] = 5.0 * tr_x_yyyy_xyyz[i] * fe_0 + 2.0 * tr_x_yyyyy_xyz[i] * fe_0 + tr_x_yyyyy_xyyz[i] * pa_y[i];

        tr_x_yyyyyy_xyzz[i] = 5.0 * tr_x_yyyy_xyzz[i] * fe_0 + tr_x_yyyyy_xzz[i] * fe_0 + tr_x_yyyyy_xyzz[i] * pa_y[i];

        tr_x_yyyyyy_xzzz[i] = 5.0 * tr_x_yyyy_xzzz[i] * fe_0 + tr_x_yyyyy_xzzz[i] * pa_y[i];

        tr_x_yyyyyy_yyyy[i] = 5.0 * tr_x_yyyy_yyyy[i] * fe_0 + 4.0 * tr_x_yyyyy_yyy[i] * fe_0 + tr_x_yyyyy_yyyy[i] * pa_y[i];

        tr_x_yyyyyy_yyyz[i] = 5.0 * tr_x_yyyy_yyyz[i] * fe_0 + 3.0 * tr_x_yyyyy_yyz[i] * fe_0 + tr_x_yyyyy_yyyz[i] * pa_y[i];

        tr_x_yyyyyy_yyzz[i] = 5.0 * tr_x_yyyy_yyzz[i] * fe_0 + 2.0 * tr_x_yyyyy_yzz[i] * fe_0 + tr_x_yyyyy_yyzz[i] * pa_y[i];

        tr_x_yyyyyy_yzzz[i] = 5.0 * tr_x_yyyy_yzzz[i] * fe_0 + tr_x_yyyyy_zzz[i] * fe_0 + tr_x_yyyyy_yzzz[i] * pa_y[i];

        tr_x_yyyyyy_zzzz[i] = 5.0 * tr_x_yyyy_zzzz[i] * fe_0 + tr_x_yyyyy_zzzz[i] * pa_y[i];
    }

    // Set up 330-345 components of targeted buffer : IG

    auto tr_x_yyyyyz_xxxx = pbuffer.data(idx_dip_ig + 330);

    auto tr_x_yyyyyz_xxxy = pbuffer.data(idx_dip_ig + 331);

    auto tr_x_yyyyyz_xxxz = pbuffer.data(idx_dip_ig + 332);

    auto tr_x_yyyyyz_xxyy = pbuffer.data(idx_dip_ig + 333);

    auto tr_x_yyyyyz_xxyz = pbuffer.data(idx_dip_ig + 334);

    auto tr_x_yyyyyz_xxzz = pbuffer.data(idx_dip_ig + 335);

    auto tr_x_yyyyyz_xyyy = pbuffer.data(idx_dip_ig + 336);

    auto tr_x_yyyyyz_xyyz = pbuffer.data(idx_dip_ig + 337);

    auto tr_x_yyyyyz_xyzz = pbuffer.data(idx_dip_ig + 338);

    auto tr_x_yyyyyz_xzzz = pbuffer.data(idx_dip_ig + 339);

    auto tr_x_yyyyyz_yyyy = pbuffer.data(idx_dip_ig + 340);

    auto tr_x_yyyyyz_yyyz = pbuffer.data(idx_dip_ig + 341);

    auto tr_x_yyyyyz_yyzz = pbuffer.data(idx_dip_ig + 342);

    auto tr_x_yyyyyz_yzzz = pbuffer.data(idx_dip_ig + 343);

    auto tr_x_yyyyyz_zzzz = pbuffer.data(idx_dip_ig + 344);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyyyy_xxxx,  \
                             tr_x_yyyyy_xxxy,  \
                             tr_x_yyyyy_xxy,   \
                             tr_x_yyyyy_xxyy,  \
                             tr_x_yyyyy_xxyz,  \
                             tr_x_yyyyy_xyy,   \
                             tr_x_yyyyy_xyyy,  \
                             tr_x_yyyyy_xyyz,  \
                             tr_x_yyyyy_xyz,   \
                             tr_x_yyyyy_xyzz,  \
                             tr_x_yyyyy_yyy,   \
                             tr_x_yyyyy_yyyy,  \
                             tr_x_yyyyy_yyyz,  \
                             tr_x_yyyyy_yyz,   \
                             tr_x_yyyyy_yyzz,  \
                             tr_x_yyyyy_yzz,   \
                             tr_x_yyyyy_yzzz,  \
                             tr_x_yyyyyz_xxxx, \
                             tr_x_yyyyyz_xxxy, \
                             tr_x_yyyyyz_xxxz, \
                             tr_x_yyyyyz_xxyy, \
                             tr_x_yyyyyz_xxyz, \
                             tr_x_yyyyyz_xxzz, \
                             tr_x_yyyyyz_xyyy, \
                             tr_x_yyyyyz_xyyz, \
                             tr_x_yyyyyz_xyzz, \
                             tr_x_yyyyyz_xzzz, \
                             tr_x_yyyyyz_yyyy, \
                             tr_x_yyyyyz_yyyz, \
                             tr_x_yyyyyz_yyzz, \
                             tr_x_yyyyyz_yzzz, \
                             tr_x_yyyyyz_zzzz, \
                             tr_x_yyyyz_xxxz,  \
                             tr_x_yyyyz_xxzz,  \
                             tr_x_yyyyz_xzzz,  \
                             tr_x_yyyyz_zzzz,  \
                             tr_x_yyyz_xxxz,   \
                             tr_x_yyyz_xxzz,   \
                             tr_x_yyyz_xzzz,   \
                             tr_x_yyyz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyz_xxxx[i] = tr_x_yyyyy_xxxx[i] * pa_z[i];

        tr_x_yyyyyz_xxxy[i] = tr_x_yyyyy_xxxy[i] * pa_z[i];

        tr_x_yyyyyz_xxxz[i] = 4.0 * tr_x_yyyz_xxxz[i] * fe_0 + tr_x_yyyyz_xxxz[i] * pa_y[i];

        tr_x_yyyyyz_xxyy[i] = tr_x_yyyyy_xxyy[i] * pa_z[i];

        tr_x_yyyyyz_xxyz[i] = tr_x_yyyyy_xxy[i] * fe_0 + tr_x_yyyyy_xxyz[i] * pa_z[i];

        tr_x_yyyyyz_xxzz[i] = 4.0 * tr_x_yyyz_xxzz[i] * fe_0 + tr_x_yyyyz_xxzz[i] * pa_y[i];

        tr_x_yyyyyz_xyyy[i] = tr_x_yyyyy_xyyy[i] * pa_z[i];

        tr_x_yyyyyz_xyyz[i] = tr_x_yyyyy_xyy[i] * fe_0 + tr_x_yyyyy_xyyz[i] * pa_z[i];

        tr_x_yyyyyz_xyzz[i] = 2.0 * tr_x_yyyyy_xyz[i] * fe_0 + tr_x_yyyyy_xyzz[i] * pa_z[i];

        tr_x_yyyyyz_xzzz[i] = 4.0 * tr_x_yyyz_xzzz[i] * fe_0 + tr_x_yyyyz_xzzz[i] * pa_y[i];

        tr_x_yyyyyz_yyyy[i] = tr_x_yyyyy_yyyy[i] * pa_z[i];

        tr_x_yyyyyz_yyyz[i] = tr_x_yyyyy_yyy[i] * fe_0 + tr_x_yyyyy_yyyz[i] * pa_z[i];

        tr_x_yyyyyz_yyzz[i] = 2.0 * tr_x_yyyyy_yyz[i] * fe_0 + tr_x_yyyyy_yyzz[i] * pa_z[i];

        tr_x_yyyyyz_yzzz[i] = 3.0 * tr_x_yyyyy_yzz[i] * fe_0 + tr_x_yyyyy_yzzz[i] * pa_z[i];

        tr_x_yyyyyz_zzzz[i] = 4.0 * tr_x_yyyz_zzzz[i] * fe_0 + tr_x_yyyyz_zzzz[i] * pa_y[i];
    }

    // Set up 345-360 components of targeted buffer : IG

    auto tr_x_yyyyzz_xxxx = pbuffer.data(idx_dip_ig + 345);

    auto tr_x_yyyyzz_xxxy = pbuffer.data(idx_dip_ig + 346);

    auto tr_x_yyyyzz_xxxz = pbuffer.data(idx_dip_ig + 347);

    auto tr_x_yyyyzz_xxyy = pbuffer.data(idx_dip_ig + 348);

    auto tr_x_yyyyzz_xxyz = pbuffer.data(idx_dip_ig + 349);

    auto tr_x_yyyyzz_xxzz = pbuffer.data(idx_dip_ig + 350);

    auto tr_x_yyyyzz_xyyy = pbuffer.data(idx_dip_ig + 351);

    auto tr_x_yyyyzz_xyyz = pbuffer.data(idx_dip_ig + 352);

    auto tr_x_yyyyzz_xyzz = pbuffer.data(idx_dip_ig + 353);

    auto tr_x_yyyyzz_xzzz = pbuffer.data(idx_dip_ig + 354);

    auto tr_x_yyyyzz_yyyy = pbuffer.data(idx_dip_ig + 355);

    auto tr_x_yyyyzz_yyyz = pbuffer.data(idx_dip_ig + 356);

    auto tr_x_yyyyzz_yyzz = pbuffer.data(idx_dip_ig + 357);

    auto tr_x_yyyyzz_yzzz = pbuffer.data(idx_dip_ig + 358);

    auto tr_x_yyyyzz_zzzz = pbuffer.data(idx_dip_ig + 359);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyyy_xxxy,   \
                             tr_x_yyyy_xxyy,   \
                             tr_x_yyyy_xyyy,   \
                             tr_x_yyyy_yyyy,   \
                             tr_x_yyyyz_xxxy,  \
                             tr_x_yyyyz_xxyy,  \
                             tr_x_yyyyz_xyyy,  \
                             tr_x_yyyyz_yyyy,  \
                             tr_x_yyyyzz_xxxx, \
                             tr_x_yyyyzz_xxxy, \
                             tr_x_yyyyzz_xxxz, \
                             tr_x_yyyyzz_xxyy, \
                             tr_x_yyyyzz_xxyz, \
                             tr_x_yyyyzz_xxzz, \
                             tr_x_yyyyzz_xyyy, \
                             tr_x_yyyyzz_xyyz, \
                             tr_x_yyyyzz_xyzz, \
                             tr_x_yyyyzz_xzzz, \
                             tr_x_yyyyzz_yyyy, \
                             tr_x_yyyyzz_yyyz, \
                             tr_x_yyyyzz_yyzz, \
                             tr_x_yyyyzz_yzzz, \
                             tr_x_yyyyzz_zzzz, \
                             tr_x_yyyzz_xxxx,  \
                             tr_x_yyyzz_xxxz,  \
                             tr_x_yyyzz_xxyz,  \
                             tr_x_yyyzz_xxz,   \
                             tr_x_yyyzz_xxzz,  \
                             tr_x_yyyzz_xyyz,  \
                             tr_x_yyyzz_xyz,   \
                             tr_x_yyyzz_xyzz,  \
                             tr_x_yyyzz_xzz,   \
                             tr_x_yyyzz_xzzz,  \
                             tr_x_yyyzz_yyyz,  \
                             tr_x_yyyzz_yyz,   \
                             tr_x_yyyzz_yyzz,  \
                             tr_x_yyyzz_yzz,   \
                             tr_x_yyyzz_yzzz,  \
                             tr_x_yyyzz_zzz,   \
                             tr_x_yyyzz_zzzz,  \
                             tr_x_yyzz_xxxx,   \
                             tr_x_yyzz_xxxz,   \
                             tr_x_yyzz_xxyz,   \
                             tr_x_yyzz_xxzz,   \
                             tr_x_yyzz_xyyz,   \
                             tr_x_yyzz_xyzz,   \
                             tr_x_yyzz_xzzz,   \
                             tr_x_yyzz_yyyz,   \
                             tr_x_yyzz_yyzz,   \
                             tr_x_yyzz_yzzz,   \
                             tr_x_yyzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyzz_xxxx[i] = 3.0 * tr_x_yyzz_xxxx[i] * fe_0 + tr_x_yyyzz_xxxx[i] * pa_y[i];

        tr_x_yyyyzz_xxxy[i] = tr_x_yyyy_xxxy[i] * fe_0 + tr_x_yyyyz_xxxy[i] * pa_z[i];

        tr_x_yyyyzz_xxxz[i] = 3.0 * tr_x_yyzz_xxxz[i] * fe_0 + tr_x_yyyzz_xxxz[i] * pa_y[i];

        tr_x_yyyyzz_xxyy[i] = tr_x_yyyy_xxyy[i] * fe_0 + tr_x_yyyyz_xxyy[i] * pa_z[i];

        tr_x_yyyyzz_xxyz[i] = 3.0 * tr_x_yyzz_xxyz[i] * fe_0 + tr_x_yyyzz_xxz[i] * fe_0 + tr_x_yyyzz_xxyz[i] * pa_y[i];

        tr_x_yyyyzz_xxzz[i] = 3.0 * tr_x_yyzz_xxzz[i] * fe_0 + tr_x_yyyzz_xxzz[i] * pa_y[i];

        tr_x_yyyyzz_xyyy[i] = tr_x_yyyy_xyyy[i] * fe_0 + tr_x_yyyyz_xyyy[i] * pa_z[i];

        tr_x_yyyyzz_xyyz[i] = 3.0 * tr_x_yyzz_xyyz[i] * fe_0 + 2.0 * tr_x_yyyzz_xyz[i] * fe_0 + tr_x_yyyzz_xyyz[i] * pa_y[i];

        tr_x_yyyyzz_xyzz[i] = 3.0 * tr_x_yyzz_xyzz[i] * fe_0 + tr_x_yyyzz_xzz[i] * fe_0 + tr_x_yyyzz_xyzz[i] * pa_y[i];

        tr_x_yyyyzz_xzzz[i] = 3.0 * tr_x_yyzz_xzzz[i] * fe_0 + tr_x_yyyzz_xzzz[i] * pa_y[i];

        tr_x_yyyyzz_yyyy[i] = tr_x_yyyy_yyyy[i] * fe_0 + tr_x_yyyyz_yyyy[i] * pa_z[i];

        tr_x_yyyyzz_yyyz[i] = 3.0 * tr_x_yyzz_yyyz[i] * fe_0 + 3.0 * tr_x_yyyzz_yyz[i] * fe_0 + tr_x_yyyzz_yyyz[i] * pa_y[i];

        tr_x_yyyyzz_yyzz[i] = 3.0 * tr_x_yyzz_yyzz[i] * fe_0 + 2.0 * tr_x_yyyzz_yzz[i] * fe_0 + tr_x_yyyzz_yyzz[i] * pa_y[i];

        tr_x_yyyyzz_yzzz[i] = 3.0 * tr_x_yyzz_yzzz[i] * fe_0 + tr_x_yyyzz_zzz[i] * fe_0 + tr_x_yyyzz_yzzz[i] * pa_y[i];

        tr_x_yyyyzz_zzzz[i] = 3.0 * tr_x_yyzz_zzzz[i] * fe_0 + tr_x_yyyzz_zzzz[i] * pa_y[i];
    }

    // Set up 360-375 components of targeted buffer : IG

    auto tr_x_yyyzzz_xxxx = pbuffer.data(idx_dip_ig + 360);

    auto tr_x_yyyzzz_xxxy = pbuffer.data(idx_dip_ig + 361);

    auto tr_x_yyyzzz_xxxz = pbuffer.data(idx_dip_ig + 362);

    auto tr_x_yyyzzz_xxyy = pbuffer.data(idx_dip_ig + 363);

    auto tr_x_yyyzzz_xxyz = pbuffer.data(idx_dip_ig + 364);

    auto tr_x_yyyzzz_xxzz = pbuffer.data(idx_dip_ig + 365);

    auto tr_x_yyyzzz_xyyy = pbuffer.data(idx_dip_ig + 366);

    auto tr_x_yyyzzz_xyyz = pbuffer.data(idx_dip_ig + 367);

    auto tr_x_yyyzzz_xyzz = pbuffer.data(idx_dip_ig + 368);

    auto tr_x_yyyzzz_xzzz = pbuffer.data(idx_dip_ig + 369);

    auto tr_x_yyyzzz_yyyy = pbuffer.data(idx_dip_ig + 370);

    auto tr_x_yyyzzz_yyyz = pbuffer.data(idx_dip_ig + 371);

    auto tr_x_yyyzzz_yyzz = pbuffer.data(idx_dip_ig + 372);

    auto tr_x_yyyzzz_yzzz = pbuffer.data(idx_dip_ig + 373);

    auto tr_x_yyyzzz_zzzz = pbuffer.data(idx_dip_ig + 374);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyyz_xxxy,   \
                             tr_x_yyyz_xxyy,   \
                             tr_x_yyyz_xyyy,   \
                             tr_x_yyyz_yyyy,   \
                             tr_x_yyyzz_xxxy,  \
                             tr_x_yyyzz_xxyy,  \
                             tr_x_yyyzz_xyyy,  \
                             tr_x_yyyzz_yyyy,  \
                             tr_x_yyyzzz_xxxx, \
                             tr_x_yyyzzz_xxxy, \
                             tr_x_yyyzzz_xxxz, \
                             tr_x_yyyzzz_xxyy, \
                             tr_x_yyyzzz_xxyz, \
                             tr_x_yyyzzz_xxzz, \
                             tr_x_yyyzzz_xyyy, \
                             tr_x_yyyzzz_xyyz, \
                             tr_x_yyyzzz_xyzz, \
                             tr_x_yyyzzz_xzzz, \
                             tr_x_yyyzzz_yyyy, \
                             tr_x_yyyzzz_yyyz, \
                             tr_x_yyyzzz_yyzz, \
                             tr_x_yyyzzz_yzzz, \
                             tr_x_yyyzzz_zzzz, \
                             tr_x_yyzzz_xxxx,  \
                             tr_x_yyzzz_xxxz,  \
                             tr_x_yyzzz_xxyz,  \
                             tr_x_yyzzz_xxz,   \
                             tr_x_yyzzz_xxzz,  \
                             tr_x_yyzzz_xyyz,  \
                             tr_x_yyzzz_xyz,   \
                             tr_x_yyzzz_xyzz,  \
                             tr_x_yyzzz_xzz,   \
                             tr_x_yyzzz_xzzz,  \
                             tr_x_yyzzz_yyyz,  \
                             tr_x_yyzzz_yyz,   \
                             tr_x_yyzzz_yyzz,  \
                             tr_x_yyzzz_yzz,   \
                             tr_x_yyzzz_yzzz,  \
                             tr_x_yyzzz_zzz,   \
                             tr_x_yyzzz_zzzz,  \
                             tr_x_yzzz_xxxx,   \
                             tr_x_yzzz_xxxz,   \
                             tr_x_yzzz_xxyz,   \
                             tr_x_yzzz_xxzz,   \
                             tr_x_yzzz_xyyz,   \
                             tr_x_yzzz_xyzz,   \
                             tr_x_yzzz_xzzz,   \
                             tr_x_yzzz_yyyz,   \
                             tr_x_yzzz_yyzz,   \
                             tr_x_yzzz_yzzz,   \
                             tr_x_yzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzzz_xxxx[i] = 2.0 * tr_x_yzzz_xxxx[i] * fe_0 + tr_x_yyzzz_xxxx[i] * pa_y[i];

        tr_x_yyyzzz_xxxy[i] = 2.0 * tr_x_yyyz_xxxy[i] * fe_0 + tr_x_yyyzz_xxxy[i] * pa_z[i];

        tr_x_yyyzzz_xxxz[i] = 2.0 * tr_x_yzzz_xxxz[i] * fe_0 + tr_x_yyzzz_xxxz[i] * pa_y[i];

        tr_x_yyyzzz_xxyy[i] = 2.0 * tr_x_yyyz_xxyy[i] * fe_0 + tr_x_yyyzz_xxyy[i] * pa_z[i];

        tr_x_yyyzzz_xxyz[i] = 2.0 * tr_x_yzzz_xxyz[i] * fe_0 + tr_x_yyzzz_xxz[i] * fe_0 + tr_x_yyzzz_xxyz[i] * pa_y[i];

        tr_x_yyyzzz_xxzz[i] = 2.0 * tr_x_yzzz_xxzz[i] * fe_0 + tr_x_yyzzz_xxzz[i] * pa_y[i];

        tr_x_yyyzzz_xyyy[i] = 2.0 * tr_x_yyyz_xyyy[i] * fe_0 + tr_x_yyyzz_xyyy[i] * pa_z[i];

        tr_x_yyyzzz_xyyz[i] = 2.0 * tr_x_yzzz_xyyz[i] * fe_0 + 2.0 * tr_x_yyzzz_xyz[i] * fe_0 + tr_x_yyzzz_xyyz[i] * pa_y[i];

        tr_x_yyyzzz_xyzz[i] = 2.0 * tr_x_yzzz_xyzz[i] * fe_0 + tr_x_yyzzz_xzz[i] * fe_0 + tr_x_yyzzz_xyzz[i] * pa_y[i];

        tr_x_yyyzzz_xzzz[i] = 2.0 * tr_x_yzzz_xzzz[i] * fe_0 + tr_x_yyzzz_xzzz[i] * pa_y[i];

        tr_x_yyyzzz_yyyy[i] = 2.0 * tr_x_yyyz_yyyy[i] * fe_0 + tr_x_yyyzz_yyyy[i] * pa_z[i];

        tr_x_yyyzzz_yyyz[i] = 2.0 * tr_x_yzzz_yyyz[i] * fe_0 + 3.0 * tr_x_yyzzz_yyz[i] * fe_0 + tr_x_yyzzz_yyyz[i] * pa_y[i];

        tr_x_yyyzzz_yyzz[i] = 2.0 * tr_x_yzzz_yyzz[i] * fe_0 + 2.0 * tr_x_yyzzz_yzz[i] * fe_0 + tr_x_yyzzz_yyzz[i] * pa_y[i];

        tr_x_yyyzzz_yzzz[i] = 2.0 * tr_x_yzzz_yzzz[i] * fe_0 + tr_x_yyzzz_zzz[i] * fe_0 + tr_x_yyzzz_yzzz[i] * pa_y[i];

        tr_x_yyyzzz_zzzz[i] = 2.0 * tr_x_yzzz_zzzz[i] * fe_0 + tr_x_yyzzz_zzzz[i] * pa_y[i];
    }

    // Set up 375-390 components of targeted buffer : IG

    auto tr_x_yyzzzz_xxxx = pbuffer.data(idx_dip_ig + 375);

    auto tr_x_yyzzzz_xxxy = pbuffer.data(idx_dip_ig + 376);

    auto tr_x_yyzzzz_xxxz = pbuffer.data(idx_dip_ig + 377);

    auto tr_x_yyzzzz_xxyy = pbuffer.data(idx_dip_ig + 378);

    auto tr_x_yyzzzz_xxyz = pbuffer.data(idx_dip_ig + 379);

    auto tr_x_yyzzzz_xxzz = pbuffer.data(idx_dip_ig + 380);

    auto tr_x_yyzzzz_xyyy = pbuffer.data(idx_dip_ig + 381);

    auto tr_x_yyzzzz_xyyz = pbuffer.data(idx_dip_ig + 382);

    auto tr_x_yyzzzz_xyzz = pbuffer.data(idx_dip_ig + 383);

    auto tr_x_yyzzzz_xzzz = pbuffer.data(idx_dip_ig + 384);

    auto tr_x_yyzzzz_yyyy = pbuffer.data(idx_dip_ig + 385);

    auto tr_x_yyzzzz_yyyz = pbuffer.data(idx_dip_ig + 386);

    auto tr_x_yyzzzz_yyzz = pbuffer.data(idx_dip_ig + 387);

    auto tr_x_yyzzzz_yzzz = pbuffer.data(idx_dip_ig + 388);

    auto tr_x_yyzzzz_zzzz = pbuffer.data(idx_dip_ig + 389);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyzz_xxxy,   \
                             tr_x_yyzz_xxyy,   \
                             tr_x_yyzz_xyyy,   \
                             tr_x_yyzz_yyyy,   \
                             tr_x_yyzzz_xxxy,  \
                             tr_x_yyzzz_xxyy,  \
                             tr_x_yyzzz_xyyy,  \
                             tr_x_yyzzz_yyyy,  \
                             tr_x_yyzzzz_xxxx, \
                             tr_x_yyzzzz_xxxy, \
                             tr_x_yyzzzz_xxxz, \
                             tr_x_yyzzzz_xxyy, \
                             tr_x_yyzzzz_xxyz, \
                             tr_x_yyzzzz_xxzz, \
                             tr_x_yyzzzz_xyyy, \
                             tr_x_yyzzzz_xyyz, \
                             tr_x_yyzzzz_xyzz, \
                             tr_x_yyzzzz_xzzz, \
                             tr_x_yyzzzz_yyyy, \
                             tr_x_yyzzzz_yyyz, \
                             tr_x_yyzzzz_yyzz, \
                             tr_x_yyzzzz_yzzz, \
                             tr_x_yyzzzz_zzzz, \
                             tr_x_yzzzz_xxxx,  \
                             tr_x_yzzzz_xxxz,  \
                             tr_x_yzzzz_xxyz,  \
                             tr_x_yzzzz_xxz,   \
                             tr_x_yzzzz_xxzz,  \
                             tr_x_yzzzz_xyyz,  \
                             tr_x_yzzzz_xyz,   \
                             tr_x_yzzzz_xyzz,  \
                             tr_x_yzzzz_xzz,   \
                             tr_x_yzzzz_xzzz,  \
                             tr_x_yzzzz_yyyz,  \
                             tr_x_yzzzz_yyz,   \
                             tr_x_yzzzz_yyzz,  \
                             tr_x_yzzzz_yzz,   \
                             tr_x_yzzzz_yzzz,  \
                             tr_x_yzzzz_zzz,   \
                             tr_x_yzzzz_zzzz,  \
                             tr_x_zzzz_xxxx,   \
                             tr_x_zzzz_xxxz,   \
                             tr_x_zzzz_xxyz,   \
                             tr_x_zzzz_xxzz,   \
                             tr_x_zzzz_xyyz,   \
                             tr_x_zzzz_xyzz,   \
                             tr_x_zzzz_xzzz,   \
                             tr_x_zzzz_yyyz,   \
                             tr_x_zzzz_yyzz,   \
                             tr_x_zzzz_yzzz,   \
                             tr_x_zzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzzz_xxxx[i] = tr_x_zzzz_xxxx[i] * fe_0 + tr_x_yzzzz_xxxx[i] * pa_y[i];

        tr_x_yyzzzz_xxxy[i] = 3.0 * tr_x_yyzz_xxxy[i] * fe_0 + tr_x_yyzzz_xxxy[i] * pa_z[i];

        tr_x_yyzzzz_xxxz[i] = tr_x_zzzz_xxxz[i] * fe_0 + tr_x_yzzzz_xxxz[i] * pa_y[i];

        tr_x_yyzzzz_xxyy[i] = 3.0 * tr_x_yyzz_xxyy[i] * fe_0 + tr_x_yyzzz_xxyy[i] * pa_z[i];

        tr_x_yyzzzz_xxyz[i] = tr_x_zzzz_xxyz[i] * fe_0 + tr_x_yzzzz_xxz[i] * fe_0 + tr_x_yzzzz_xxyz[i] * pa_y[i];

        tr_x_yyzzzz_xxzz[i] = tr_x_zzzz_xxzz[i] * fe_0 + tr_x_yzzzz_xxzz[i] * pa_y[i];

        tr_x_yyzzzz_xyyy[i] = 3.0 * tr_x_yyzz_xyyy[i] * fe_0 + tr_x_yyzzz_xyyy[i] * pa_z[i];

        tr_x_yyzzzz_xyyz[i] = tr_x_zzzz_xyyz[i] * fe_0 + 2.0 * tr_x_yzzzz_xyz[i] * fe_0 + tr_x_yzzzz_xyyz[i] * pa_y[i];

        tr_x_yyzzzz_xyzz[i] = tr_x_zzzz_xyzz[i] * fe_0 + tr_x_yzzzz_xzz[i] * fe_0 + tr_x_yzzzz_xyzz[i] * pa_y[i];

        tr_x_yyzzzz_xzzz[i] = tr_x_zzzz_xzzz[i] * fe_0 + tr_x_yzzzz_xzzz[i] * pa_y[i];

        tr_x_yyzzzz_yyyy[i] = 3.0 * tr_x_yyzz_yyyy[i] * fe_0 + tr_x_yyzzz_yyyy[i] * pa_z[i];

        tr_x_yyzzzz_yyyz[i] = tr_x_zzzz_yyyz[i] * fe_0 + 3.0 * tr_x_yzzzz_yyz[i] * fe_0 + tr_x_yzzzz_yyyz[i] * pa_y[i];

        tr_x_yyzzzz_yyzz[i] = tr_x_zzzz_yyzz[i] * fe_0 + 2.0 * tr_x_yzzzz_yzz[i] * fe_0 + tr_x_yzzzz_yyzz[i] * pa_y[i];

        tr_x_yyzzzz_yzzz[i] = tr_x_zzzz_yzzz[i] * fe_0 + tr_x_yzzzz_zzz[i] * fe_0 + tr_x_yzzzz_yzzz[i] * pa_y[i];

        tr_x_yyzzzz_zzzz[i] = tr_x_zzzz_zzzz[i] * fe_0 + tr_x_yzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 390-405 components of targeted buffer : IG

    auto tr_x_yzzzzz_xxxx = pbuffer.data(idx_dip_ig + 390);

    auto tr_x_yzzzzz_xxxy = pbuffer.data(idx_dip_ig + 391);

    auto tr_x_yzzzzz_xxxz = pbuffer.data(idx_dip_ig + 392);

    auto tr_x_yzzzzz_xxyy = pbuffer.data(idx_dip_ig + 393);

    auto tr_x_yzzzzz_xxyz = pbuffer.data(idx_dip_ig + 394);

    auto tr_x_yzzzzz_xxzz = pbuffer.data(idx_dip_ig + 395);

    auto tr_x_yzzzzz_xyyy = pbuffer.data(idx_dip_ig + 396);

    auto tr_x_yzzzzz_xyyz = pbuffer.data(idx_dip_ig + 397);

    auto tr_x_yzzzzz_xyzz = pbuffer.data(idx_dip_ig + 398);

    auto tr_x_yzzzzz_xzzz = pbuffer.data(idx_dip_ig + 399);

    auto tr_x_yzzzzz_yyyy = pbuffer.data(idx_dip_ig + 400);

    auto tr_x_yzzzzz_yyyz = pbuffer.data(idx_dip_ig + 401);

    auto tr_x_yzzzzz_yyzz = pbuffer.data(idx_dip_ig + 402);

    auto tr_x_yzzzzz_yzzz = pbuffer.data(idx_dip_ig + 403);

    auto tr_x_yzzzzz_zzzz = pbuffer.data(idx_dip_ig + 404);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_yzzzzz_xxxx, \
                             tr_x_yzzzzz_xxxy, \
                             tr_x_yzzzzz_xxxz, \
                             tr_x_yzzzzz_xxyy, \
                             tr_x_yzzzzz_xxyz, \
                             tr_x_yzzzzz_xxzz, \
                             tr_x_yzzzzz_xyyy, \
                             tr_x_yzzzzz_xyyz, \
                             tr_x_yzzzzz_xyzz, \
                             tr_x_yzzzzz_xzzz, \
                             tr_x_yzzzzz_yyyy, \
                             tr_x_yzzzzz_yyyz, \
                             tr_x_yzzzzz_yyzz, \
                             tr_x_yzzzzz_yzzz, \
                             tr_x_yzzzzz_zzzz, \
                             tr_x_zzzzz_xxx,   \
                             tr_x_zzzzz_xxxx,  \
                             tr_x_zzzzz_xxxy,  \
                             tr_x_zzzzz_xxxz,  \
                             tr_x_zzzzz_xxy,   \
                             tr_x_zzzzz_xxyy,  \
                             tr_x_zzzzz_xxyz,  \
                             tr_x_zzzzz_xxz,   \
                             tr_x_zzzzz_xxzz,  \
                             tr_x_zzzzz_xyy,   \
                             tr_x_zzzzz_xyyy,  \
                             tr_x_zzzzz_xyyz,  \
                             tr_x_zzzzz_xyz,   \
                             tr_x_zzzzz_xyzz,  \
                             tr_x_zzzzz_xzz,   \
                             tr_x_zzzzz_xzzz,  \
                             tr_x_zzzzz_yyy,   \
                             tr_x_zzzzz_yyyy,  \
                             tr_x_zzzzz_yyyz,  \
                             tr_x_zzzzz_yyz,   \
                             tr_x_zzzzz_yyzz,  \
                             tr_x_zzzzz_yzz,   \
                             tr_x_zzzzz_yzzz,  \
                             tr_x_zzzzz_zzz,   \
                             tr_x_zzzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzzz_xxxx[i] = tr_x_zzzzz_xxxx[i] * pa_y[i];

        tr_x_yzzzzz_xxxy[i] = tr_x_zzzzz_xxx[i] * fe_0 + tr_x_zzzzz_xxxy[i] * pa_y[i];

        tr_x_yzzzzz_xxxz[i] = tr_x_zzzzz_xxxz[i] * pa_y[i];

        tr_x_yzzzzz_xxyy[i] = 2.0 * tr_x_zzzzz_xxy[i] * fe_0 + tr_x_zzzzz_xxyy[i] * pa_y[i];

        tr_x_yzzzzz_xxyz[i] = tr_x_zzzzz_xxz[i] * fe_0 + tr_x_zzzzz_xxyz[i] * pa_y[i];

        tr_x_yzzzzz_xxzz[i] = tr_x_zzzzz_xxzz[i] * pa_y[i];

        tr_x_yzzzzz_xyyy[i] = 3.0 * tr_x_zzzzz_xyy[i] * fe_0 + tr_x_zzzzz_xyyy[i] * pa_y[i];

        tr_x_yzzzzz_xyyz[i] = 2.0 * tr_x_zzzzz_xyz[i] * fe_0 + tr_x_zzzzz_xyyz[i] * pa_y[i];

        tr_x_yzzzzz_xyzz[i] = tr_x_zzzzz_xzz[i] * fe_0 + tr_x_zzzzz_xyzz[i] * pa_y[i];

        tr_x_yzzzzz_xzzz[i] = tr_x_zzzzz_xzzz[i] * pa_y[i];

        tr_x_yzzzzz_yyyy[i] = 4.0 * tr_x_zzzzz_yyy[i] * fe_0 + tr_x_zzzzz_yyyy[i] * pa_y[i];

        tr_x_yzzzzz_yyyz[i] = 3.0 * tr_x_zzzzz_yyz[i] * fe_0 + tr_x_zzzzz_yyyz[i] * pa_y[i];

        tr_x_yzzzzz_yyzz[i] = 2.0 * tr_x_zzzzz_yzz[i] * fe_0 + tr_x_zzzzz_yyzz[i] * pa_y[i];

        tr_x_yzzzzz_yzzz[i] = tr_x_zzzzz_zzz[i] * fe_0 + tr_x_zzzzz_yzzz[i] * pa_y[i];

        tr_x_yzzzzz_zzzz[i] = tr_x_zzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 405-420 components of targeted buffer : IG

    auto tr_x_zzzzzz_xxxx = pbuffer.data(idx_dip_ig + 405);

    auto tr_x_zzzzzz_xxxy = pbuffer.data(idx_dip_ig + 406);

    auto tr_x_zzzzzz_xxxz = pbuffer.data(idx_dip_ig + 407);

    auto tr_x_zzzzzz_xxyy = pbuffer.data(idx_dip_ig + 408);

    auto tr_x_zzzzzz_xxyz = pbuffer.data(idx_dip_ig + 409);

    auto tr_x_zzzzzz_xxzz = pbuffer.data(idx_dip_ig + 410);

    auto tr_x_zzzzzz_xyyy = pbuffer.data(idx_dip_ig + 411);

    auto tr_x_zzzzzz_xyyz = pbuffer.data(idx_dip_ig + 412);

    auto tr_x_zzzzzz_xyzz = pbuffer.data(idx_dip_ig + 413);

    auto tr_x_zzzzzz_xzzz = pbuffer.data(idx_dip_ig + 414);

    auto tr_x_zzzzzz_yyyy = pbuffer.data(idx_dip_ig + 415);

    auto tr_x_zzzzzz_yyyz = pbuffer.data(idx_dip_ig + 416);

    auto tr_x_zzzzzz_yyzz = pbuffer.data(idx_dip_ig + 417);

    auto tr_x_zzzzzz_yzzz = pbuffer.data(idx_dip_ig + 418);

    auto tr_x_zzzzzz_zzzz = pbuffer.data(idx_dip_ig + 419);

#pragma omp simd aligned(pa_z,                 \
                             tr_x_zzzz_xxxx,   \
                             tr_x_zzzz_xxxy,   \
                             tr_x_zzzz_xxxz,   \
                             tr_x_zzzz_xxyy,   \
                             tr_x_zzzz_xxyz,   \
                             tr_x_zzzz_xxzz,   \
                             tr_x_zzzz_xyyy,   \
                             tr_x_zzzz_xyyz,   \
                             tr_x_zzzz_xyzz,   \
                             tr_x_zzzz_xzzz,   \
                             tr_x_zzzz_yyyy,   \
                             tr_x_zzzz_yyyz,   \
                             tr_x_zzzz_yyzz,   \
                             tr_x_zzzz_yzzz,   \
                             tr_x_zzzz_zzzz,   \
                             tr_x_zzzzz_xxx,   \
                             tr_x_zzzzz_xxxx,  \
                             tr_x_zzzzz_xxxy,  \
                             tr_x_zzzzz_xxxz,  \
                             tr_x_zzzzz_xxy,   \
                             tr_x_zzzzz_xxyy,  \
                             tr_x_zzzzz_xxyz,  \
                             tr_x_zzzzz_xxz,   \
                             tr_x_zzzzz_xxzz,  \
                             tr_x_zzzzz_xyy,   \
                             tr_x_zzzzz_xyyy,  \
                             tr_x_zzzzz_xyyz,  \
                             tr_x_zzzzz_xyz,   \
                             tr_x_zzzzz_xyzz,  \
                             tr_x_zzzzz_xzz,   \
                             tr_x_zzzzz_xzzz,  \
                             tr_x_zzzzz_yyy,   \
                             tr_x_zzzzz_yyyy,  \
                             tr_x_zzzzz_yyyz,  \
                             tr_x_zzzzz_yyz,   \
                             tr_x_zzzzz_yyzz,  \
                             tr_x_zzzzz_yzz,   \
                             tr_x_zzzzz_yzzz,  \
                             tr_x_zzzzz_zzz,   \
                             tr_x_zzzzz_zzzz,  \
                             tr_x_zzzzzz_xxxx, \
                             tr_x_zzzzzz_xxxy, \
                             tr_x_zzzzzz_xxxz, \
                             tr_x_zzzzzz_xxyy, \
                             tr_x_zzzzzz_xxyz, \
                             tr_x_zzzzzz_xxzz, \
                             tr_x_zzzzzz_xyyy, \
                             tr_x_zzzzzz_xyyz, \
                             tr_x_zzzzzz_xyzz, \
                             tr_x_zzzzzz_xzzz, \
                             tr_x_zzzzzz_yyyy, \
                             tr_x_zzzzzz_yyyz, \
                             tr_x_zzzzzz_yyzz, \
                             tr_x_zzzzzz_yzzz, \
                             tr_x_zzzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzzz_xxxx[i] = 5.0 * tr_x_zzzz_xxxx[i] * fe_0 + tr_x_zzzzz_xxxx[i] * pa_z[i];

        tr_x_zzzzzz_xxxy[i] = 5.0 * tr_x_zzzz_xxxy[i] * fe_0 + tr_x_zzzzz_xxxy[i] * pa_z[i];

        tr_x_zzzzzz_xxxz[i] = 5.0 * tr_x_zzzz_xxxz[i] * fe_0 + tr_x_zzzzz_xxx[i] * fe_0 + tr_x_zzzzz_xxxz[i] * pa_z[i];

        tr_x_zzzzzz_xxyy[i] = 5.0 * tr_x_zzzz_xxyy[i] * fe_0 + tr_x_zzzzz_xxyy[i] * pa_z[i];

        tr_x_zzzzzz_xxyz[i] = 5.0 * tr_x_zzzz_xxyz[i] * fe_0 + tr_x_zzzzz_xxy[i] * fe_0 + tr_x_zzzzz_xxyz[i] * pa_z[i];

        tr_x_zzzzzz_xxzz[i] = 5.0 * tr_x_zzzz_xxzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xxz[i] * fe_0 + tr_x_zzzzz_xxzz[i] * pa_z[i];

        tr_x_zzzzzz_xyyy[i] = 5.0 * tr_x_zzzz_xyyy[i] * fe_0 + tr_x_zzzzz_xyyy[i] * pa_z[i];

        tr_x_zzzzzz_xyyz[i] = 5.0 * tr_x_zzzz_xyyz[i] * fe_0 + tr_x_zzzzz_xyy[i] * fe_0 + tr_x_zzzzz_xyyz[i] * pa_z[i];

        tr_x_zzzzzz_xyzz[i] = 5.0 * tr_x_zzzz_xyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xyz[i] * fe_0 + tr_x_zzzzz_xyzz[i] * pa_z[i];

        tr_x_zzzzzz_xzzz[i] = 5.0 * tr_x_zzzz_xzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_xzz[i] * fe_0 + tr_x_zzzzz_xzzz[i] * pa_z[i];

        tr_x_zzzzzz_yyyy[i] = 5.0 * tr_x_zzzz_yyyy[i] * fe_0 + tr_x_zzzzz_yyyy[i] * pa_z[i];

        tr_x_zzzzzz_yyyz[i] = 5.0 * tr_x_zzzz_yyyz[i] * fe_0 + tr_x_zzzzz_yyy[i] * fe_0 + tr_x_zzzzz_yyyz[i] * pa_z[i];

        tr_x_zzzzzz_yyzz[i] = 5.0 * tr_x_zzzz_yyzz[i] * fe_0 + 2.0 * tr_x_zzzzz_yyz[i] * fe_0 + tr_x_zzzzz_yyzz[i] * pa_z[i];

        tr_x_zzzzzz_yzzz[i] = 5.0 * tr_x_zzzz_yzzz[i] * fe_0 + 3.0 * tr_x_zzzzz_yzz[i] * fe_0 + tr_x_zzzzz_yzzz[i] * pa_z[i];

        tr_x_zzzzzz_zzzz[i] = 5.0 * tr_x_zzzz_zzzz[i] * fe_0 + 4.0 * tr_x_zzzzz_zzz[i] * fe_0 + tr_x_zzzzz_zzzz[i] * pa_z[i];
    }

    // Set up 420-435 components of targeted buffer : IG

    auto tr_y_xxxxxx_xxxx = pbuffer.data(idx_dip_ig + 420);

    auto tr_y_xxxxxx_xxxy = pbuffer.data(idx_dip_ig + 421);

    auto tr_y_xxxxxx_xxxz = pbuffer.data(idx_dip_ig + 422);

    auto tr_y_xxxxxx_xxyy = pbuffer.data(idx_dip_ig + 423);

    auto tr_y_xxxxxx_xxyz = pbuffer.data(idx_dip_ig + 424);

    auto tr_y_xxxxxx_xxzz = pbuffer.data(idx_dip_ig + 425);

    auto tr_y_xxxxxx_xyyy = pbuffer.data(idx_dip_ig + 426);

    auto tr_y_xxxxxx_xyyz = pbuffer.data(idx_dip_ig + 427);

    auto tr_y_xxxxxx_xyzz = pbuffer.data(idx_dip_ig + 428);

    auto tr_y_xxxxxx_xzzz = pbuffer.data(idx_dip_ig + 429);

    auto tr_y_xxxxxx_yyyy = pbuffer.data(idx_dip_ig + 430);

    auto tr_y_xxxxxx_yyyz = pbuffer.data(idx_dip_ig + 431);

    auto tr_y_xxxxxx_yyzz = pbuffer.data(idx_dip_ig + 432);

    auto tr_y_xxxxxx_yzzz = pbuffer.data(idx_dip_ig + 433);

    auto tr_y_xxxxxx_zzzz = pbuffer.data(idx_dip_ig + 434);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxxx_xxxx,   \
                             tr_y_xxxx_xxxy,   \
                             tr_y_xxxx_xxxz,   \
                             tr_y_xxxx_xxyy,   \
                             tr_y_xxxx_xxyz,   \
                             tr_y_xxxx_xxzz,   \
                             tr_y_xxxx_xyyy,   \
                             tr_y_xxxx_xyyz,   \
                             tr_y_xxxx_xyzz,   \
                             tr_y_xxxx_xzzz,   \
                             tr_y_xxxx_yyyy,   \
                             tr_y_xxxx_yyyz,   \
                             tr_y_xxxx_yyzz,   \
                             tr_y_xxxx_yzzz,   \
                             tr_y_xxxx_zzzz,   \
                             tr_y_xxxxx_xxx,   \
                             tr_y_xxxxx_xxxx,  \
                             tr_y_xxxxx_xxxy,  \
                             tr_y_xxxxx_xxxz,  \
                             tr_y_xxxxx_xxy,   \
                             tr_y_xxxxx_xxyy,  \
                             tr_y_xxxxx_xxyz,  \
                             tr_y_xxxxx_xxz,   \
                             tr_y_xxxxx_xxzz,  \
                             tr_y_xxxxx_xyy,   \
                             tr_y_xxxxx_xyyy,  \
                             tr_y_xxxxx_xyyz,  \
                             tr_y_xxxxx_xyz,   \
                             tr_y_xxxxx_xyzz,  \
                             tr_y_xxxxx_xzz,   \
                             tr_y_xxxxx_xzzz,  \
                             tr_y_xxxxx_yyy,   \
                             tr_y_xxxxx_yyyy,  \
                             tr_y_xxxxx_yyyz,  \
                             tr_y_xxxxx_yyz,   \
                             tr_y_xxxxx_yyzz,  \
                             tr_y_xxxxx_yzz,   \
                             tr_y_xxxxx_yzzz,  \
                             tr_y_xxxxx_zzz,   \
                             tr_y_xxxxx_zzzz,  \
                             tr_y_xxxxxx_xxxx, \
                             tr_y_xxxxxx_xxxy, \
                             tr_y_xxxxxx_xxxz, \
                             tr_y_xxxxxx_xxyy, \
                             tr_y_xxxxxx_xxyz, \
                             tr_y_xxxxxx_xxzz, \
                             tr_y_xxxxxx_xyyy, \
                             tr_y_xxxxxx_xyyz, \
                             tr_y_xxxxxx_xyzz, \
                             tr_y_xxxxxx_xzzz, \
                             tr_y_xxxxxx_yyyy, \
                             tr_y_xxxxxx_yyyz, \
                             tr_y_xxxxxx_yyzz, \
                             tr_y_xxxxxx_yzzz, \
                             tr_y_xxxxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxx_xxxx[i] = 5.0 * tr_y_xxxx_xxxx[i] * fe_0 + 4.0 * tr_y_xxxxx_xxx[i] * fe_0 + tr_y_xxxxx_xxxx[i] * pa_x[i];

        tr_y_xxxxxx_xxxy[i] = 5.0 * tr_y_xxxx_xxxy[i] * fe_0 + 3.0 * tr_y_xxxxx_xxy[i] * fe_0 + tr_y_xxxxx_xxxy[i] * pa_x[i];

        tr_y_xxxxxx_xxxz[i] = 5.0 * tr_y_xxxx_xxxz[i] * fe_0 + 3.0 * tr_y_xxxxx_xxz[i] * fe_0 + tr_y_xxxxx_xxxz[i] * pa_x[i];

        tr_y_xxxxxx_xxyy[i] = 5.0 * tr_y_xxxx_xxyy[i] * fe_0 + 2.0 * tr_y_xxxxx_xyy[i] * fe_0 + tr_y_xxxxx_xxyy[i] * pa_x[i];

        tr_y_xxxxxx_xxyz[i] = 5.0 * tr_y_xxxx_xxyz[i] * fe_0 + 2.0 * tr_y_xxxxx_xyz[i] * fe_0 + tr_y_xxxxx_xxyz[i] * pa_x[i];

        tr_y_xxxxxx_xxzz[i] = 5.0 * tr_y_xxxx_xxzz[i] * fe_0 + 2.0 * tr_y_xxxxx_xzz[i] * fe_0 + tr_y_xxxxx_xxzz[i] * pa_x[i];

        tr_y_xxxxxx_xyyy[i] = 5.0 * tr_y_xxxx_xyyy[i] * fe_0 + tr_y_xxxxx_yyy[i] * fe_0 + tr_y_xxxxx_xyyy[i] * pa_x[i];

        tr_y_xxxxxx_xyyz[i] = 5.0 * tr_y_xxxx_xyyz[i] * fe_0 + tr_y_xxxxx_yyz[i] * fe_0 + tr_y_xxxxx_xyyz[i] * pa_x[i];

        tr_y_xxxxxx_xyzz[i] = 5.0 * tr_y_xxxx_xyzz[i] * fe_0 + tr_y_xxxxx_yzz[i] * fe_0 + tr_y_xxxxx_xyzz[i] * pa_x[i];

        tr_y_xxxxxx_xzzz[i] = 5.0 * tr_y_xxxx_xzzz[i] * fe_0 + tr_y_xxxxx_zzz[i] * fe_0 + tr_y_xxxxx_xzzz[i] * pa_x[i];

        tr_y_xxxxxx_yyyy[i] = 5.0 * tr_y_xxxx_yyyy[i] * fe_0 + tr_y_xxxxx_yyyy[i] * pa_x[i];

        tr_y_xxxxxx_yyyz[i] = 5.0 * tr_y_xxxx_yyyz[i] * fe_0 + tr_y_xxxxx_yyyz[i] * pa_x[i];

        tr_y_xxxxxx_yyzz[i] = 5.0 * tr_y_xxxx_yyzz[i] * fe_0 + tr_y_xxxxx_yyzz[i] * pa_x[i];

        tr_y_xxxxxx_yzzz[i] = 5.0 * tr_y_xxxx_yzzz[i] * fe_0 + tr_y_xxxxx_yzzz[i] * pa_x[i];

        tr_y_xxxxxx_zzzz[i] = 5.0 * tr_y_xxxx_zzzz[i] * fe_0 + tr_y_xxxxx_zzzz[i] * pa_x[i];
    }

    // Set up 435-450 components of targeted buffer : IG

    auto tr_y_xxxxxy_xxxx = pbuffer.data(idx_dip_ig + 435);

    auto tr_y_xxxxxy_xxxy = pbuffer.data(idx_dip_ig + 436);

    auto tr_y_xxxxxy_xxxz = pbuffer.data(idx_dip_ig + 437);

    auto tr_y_xxxxxy_xxyy = pbuffer.data(idx_dip_ig + 438);

    auto tr_y_xxxxxy_xxyz = pbuffer.data(idx_dip_ig + 439);

    auto tr_y_xxxxxy_xxzz = pbuffer.data(idx_dip_ig + 440);

    auto tr_y_xxxxxy_xyyy = pbuffer.data(idx_dip_ig + 441);

    auto tr_y_xxxxxy_xyyz = pbuffer.data(idx_dip_ig + 442);

    auto tr_y_xxxxxy_xyzz = pbuffer.data(idx_dip_ig + 443);

    auto tr_y_xxxxxy_xzzz = pbuffer.data(idx_dip_ig + 444);

    auto tr_y_xxxxxy_yyyy = pbuffer.data(idx_dip_ig + 445);

    auto tr_y_xxxxxy_yyyz = pbuffer.data(idx_dip_ig + 446);

    auto tr_y_xxxxxy_yyzz = pbuffer.data(idx_dip_ig + 447);

    auto tr_y_xxxxxy_yzzz = pbuffer.data(idx_dip_ig + 448);

    auto tr_y_xxxxxy_zzzz = pbuffer.data(idx_dip_ig + 449);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_y_xxxxx_xxxx,  \
                             tr_y_xxxxx_xxxz,  \
                             tr_y_xxxxx_xxzz,  \
                             tr_y_xxxxx_xzzz,  \
                             tr_y_xxxxxy_xxxx, \
                             tr_y_xxxxxy_xxxy, \
                             tr_y_xxxxxy_xxxz, \
                             tr_y_xxxxxy_xxyy, \
                             tr_y_xxxxxy_xxyz, \
                             tr_y_xxxxxy_xxzz, \
                             tr_y_xxxxxy_xyyy, \
                             tr_y_xxxxxy_xyyz, \
                             tr_y_xxxxxy_xyzz, \
                             tr_y_xxxxxy_xzzz, \
                             tr_y_xxxxxy_yyyy, \
                             tr_y_xxxxxy_yyyz, \
                             tr_y_xxxxxy_yyzz, \
                             tr_y_xxxxxy_yzzz, \
                             tr_y_xxxxxy_zzzz, \
                             tr_y_xxxxy_xxxy,  \
                             tr_y_xxxxy_xxy,   \
                             tr_y_xxxxy_xxyy,  \
                             tr_y_xxxxy_xxyz,  \
                             tr_y_xxxxy_xyy,   \
                             tr_y_xxxxy_xyyy,  \
                             tr_y_xxxxy_xyyz,  \
                             tr_y_xxxxy_xyz,   \
                             tr_y_xxxxy_xyzz,  \
                             tr_y_xxxxy_yyy,   \
                             tr_y_xxxxy_yyyy,  \
                             tr_y_xxxxy_yyyz,  \
                             tr_y_xxxxy_yyz,   \
                             tr_y_xxxxy_yyzz,  \
                             tr_y_xxxxy_yzz,   \
                             tr_y_xxxxy_yzzz,  \
                             tr_y_xxxxy_zzzz,  \
                             tr_y_xxxy_xxxy,   \
                             tr_y_xxxy_xxyy,   \
                             tr_y_xxxy_xxyz,   \
                             tr_y_xxxy_xyyy,   \
                             tr_y_xxxy_xyyz,   \
                             tr_y_xxxy_xyzz,   \
                             tr_y_xxxy_yyyy,   \
                             tr_y_xxxy_yyyz,   \
                             tr_y_xxxy_yyzz,   \
                             tr_y_xxxy_yzzz,   \
                             tr_y_xxxy_zzzz,   \
                             ts_xxxxx_xxxx,    \
                             ts_xxxxx_xxxz,    \
                             ts_xxxxx_xxzz,    \
                             ts_xxxxx_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxy_xxxx[i] = ts_xxxxx_xxxx[i] * fe_0 + tr_y_xxxxx_xxxx[i] * pa_y[i];

        tr_y_xxxxxy_xxxy[i] = 4.0 * tr_y_xxxy_xxxy[i] * fe_0 + 3.0 * tr_y_xxxxy_xxy[i] * fe_0 + tr_y_xxxxy_xxxy[i] * pa_x[i];

        tr_y_xxxxxy_xxxz[i] = ts_xxxxx_xxxz[i] * fe_0 + tr_y_xxxxx_xxxz[i] * pa_y[i];

        tr_y_xxxxxy_xxyy[i] = 4.0 * tr_y_xxxy_xxyy[i] * fe_0 + 2.0 * tr_y_xxxxy_xyy[i] * fe_0 + tr_y_xxxxy_xxyy[i] * pa_x[i];

        tr_y_xxxxxy_xxyz[i] = 4.0 * tr_y_xxxy_xxyz[i] * fe_0 + 2.0 * tr_y_xxxxy_xyz[i] * fe_0 + tr_y_xxxxy_xxyz[i] * pa_x[i];

        tr_y_xxxxxy_xxzz[i] = ts_xxxxx_xxzz[i] * fe_0 + tr_y_xxxxx_xxzz[i] * pa_y[i];

        tr_y_xxxxxy_xyyy[i] = 4.0 * tr_y_xxxy_xyyy[i] * fe_0 + tr_y_xxxxy_yyy[i] * fe_0 + tr_y_xxxxy_xyyy[i] * pa_x[i];

        tr_y_xxxxxy_xyyz[i] = 4.0 * tr_y_xxxy_xyyz[i] * fe_0 + tr_y_xxxxy_yyz[i] * fe_0 + tr_y_xxxxy_xyyz[i] * pa_x[i];

        tr_y_xxxxxy_xyzz[i] = 4.0 * tr_y_xxxy_xyzz[i] * fe_0 + tr_y_xxxxy_yzz[i] * fe_0 + tr_y_xxxxy_xyzz[i] * pa_x[i];

        tr_y_xxxxxy_xzzz[i] = ts_xxxxx_xzzz[i] * fe_0 + tr_y_xxxxx_xzzz[i] * pa_y[i];

        tr_y_xxxxxy_yyyy[i] = 4.0 * tr_y_xxxy_yyyy[i] * fe_0 + tr_y_xxxxy_yyyy[i] * pa_x[i];

        tr_y_xxxxxy_yyyz[i] = 4.0 * tr_y_xxxy_yyyz[i] * fe_0 + tr_y_xxxxy_yyyz[i] * pa_x[i];

        tr_y_xxxxxy_yyzz[i] = 4.0 * tr_y_xxxy_yyzz[i] * fe_0 + tr_y_xxxxy_yyzz[i] * pa_x[i];

        tr_y_xxxxxy_yzzz[i] = 4.0 * tr_y_xxxy_yzzz[i] * fe_0 + tr_y_xxxxy_yzzz[i] * pa_x[i];

        tr_y_xxxxxy_zzzz[i] = 4.0 * tr_y_xxxy_zzzz[i] * fe_0 + tr_y_xxxxy_zzzz[i] * pa_x[i];
    }

    // Set up 450-465 components of targeted buffer : IG

    auto tr_y_xxxxxz_xxxx = pbuffer.data(idx_dip_ig + 450);

    auto tr_y_xxxxxz_xxxy = pbuffer.data(idx_dip_ig + 451);

    auto tr_y_xxxxxz_xxxz = pbuffer.data(idx_dip_ig + 452);

    auto tr_y_xxxxxz_xxyy = pbuffer.data(idx_dip_ig + 453);

    auto tr_y_xxxxxz_xxyz = pbuffer.data(idx_dip_ig + 454);

    auto tr_y_xxxxxz_xxzz = pbuffer.data(idx_dip_ig + 455);

    auto tr_y_xxxxxz_xyyy = pbuffer.data(idx_dip_ig + 456);

    auto tr_y_xxxxxz_xyyz = pbuffer.data(idx_dip_ig + 457);

    auto tr_y_xxxxxz_xyzz = pbuffer.data(idx_dip_ig + 458);

    auto tr_y_xxxxxz_xzzz = pbuffer.data(idx_dip_ig + 459);

    auto tr_y_xxxxxz_yyyy = pbuffer.data(idx_dip_ig + 460);

    auto tr_y_xxxxxz_yyyz = pbuffer.data(idx_dip_ig + 461);

    auto tr_y_xxxxxz_yyzz = pbuffer.data(idx_dip_ig + 462);

    auto tr_y_xxxxxz_yzzz = pbuffer.data(idx_dip_ig + 463);

    auto tr_y_xxxxxz_zzzz = pbuffer.data(idx_dip_ig + 464);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxxxx_xxx,   \
                             tr_y_xxxxx_xxxx,  \
                             tr_y_xxxxx_xxxy,  \
                             tr_y_xxxxx_xxxz,  \
                             tr_y_xxxxx_xxy,   \
                             tr_y_xxxxx_xxyy,  \
                             tr_y_xxxxx_xxyz,  \
                             tr_y_xxxxx_xxz,   \
                             tr_y_xxxxx_xxzz,  \
                             tr_y_xxxxx_xyy,   \
                             tr_y_xxxxx_xyyy,  \
                             tr_y_xxxxx_xyyz,  \
                             tr_y_xxxxx_xyz,   \
                             tr_y_xxxxx_xyzz,  \
                             tr_y_xxxxx_xzz,   \
                             tr_y_xxxxx_xzzz,  \
                             tr_y_xxxxx_yyyy,  \
                             tr_y_xxxxxz_xxxx, \
                             tr_y_xxxxxz_xxxy, \
                             tr_y_xxxxxz_xxxz, \
                             tr_y_xxxxxz_xxyy, \
                             tr_y_xxxxxz_xxyz, \
                             tr_y_xxxxxz_xxzz, \
                             tr_y_xxxxxz_xyyy, \
                             tr_y_xxxxxz_xyyz, \
                             tr_y_xxxxxz_xyzz, \
                             tr_y_xxxxxz_xzzz, \
                             tr_y_xxxxxz_yyyy, \
                             tr_y_xxxxxz_yyyz, \
                             tr_y_xxxxxz_yyzz, \
                             tr_y_xxxxxz_yzzz, \
                             tr_y_xxxxxz_zzzz, \
                             tr_y_xxxxz_yyyz,  \
                             tr_y_xxxxz_yyzz,  \
                             tr_y_xxxxz_yzzz,  \
                             tr_y_xxxxz_zzzz,  \
                             tr_y_xxxz_yyyz,   \
                             tr_y_xxxz_yyzz,   \
                             tr_y_xxxz_yzzz,   \
                             tr_y_xxxz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxz_xxxx[i] = tr_y_xxxxx_xxxx[i] * pa_z[i];

        tr_y_xxxxxz_xxxy[i] = tr_y_xxxxx_xxxy[i] * pa_z[i];

        tr_y_xxxxxz_xxxz[i] = tr_y_xxxxx_xxx[i] * fe_0 + tr_y_xxxxx_xxxz[i] * pa_z[i];

        tr_y_xxxxxz_xxyy[i] = tr_y_xxxxx_xxyy[i] * pa_z[i];

        tr_y_xxxxxz_xxyz[i] = tr_y_xxxxx_xxy[i] * fe_0 + tr_y_xxxxx_xxyz[i] * pa_z[i];

        tr_y_xxxxxz_xxzz[i] = 2.0 * tr_y_xxxxx_xxz[i] * fe_0 + tr_y_xxxxx_xxzz[i] * pa_z[i];

        tr_y_xxxxxz_xyyy[i] = tr_y_xxxxx_xyyy[i] * pa_z[i];

        tr_y_xxxxxz_xyyz[i] = tr_y_xxxxx_xyy[i] * fe_0 + tr_y_xxxxx_xyyz[i] * pa_z[i];

        tr_y_xxxxxz_xyzz[i] = 2.0 * tr_y_xxxxx_xyz[i] * fe_0 + tr_y_xxxxx_xyzz[i] * pa_z[i];

        tr_y_xxxxxz_xzzz[i] = 3.0 * tr_y_xxxxx_xzz[i] * fe_0 + tr_y_xxxxx_xzzz[i] * pa_z[i];

        tr_y_xxxxxz_yyyy[i] = tr_y_xxxxx_yyyy[i] * pa_z[i];

        tr_y_xxxxxz_yyyz[i] = 4.0 * tr_y_xxxz_yyyz[i] * fe_0 + tr_y_xxxxz_yyyz[i] * pa_x[i];

        tr_y_xxxxxz_yyzz[i] = 4.0 * tr_y_xxxz_yyzz[i] * fe_0 + tr_y_xxxxz_yyzz[i] * pa_x[i];

        tr_y_xxxxxz_yzzz[i] = 4.0 * tr_y_xxxz_yzzz[i] * fe_0 + tr_y_xxxxz_yzzz[i] * pa_x[i];

        tr_y_xxxxxz_zzzz[i] = 4.0 * tr_y_xxxz_zzzz[i] * fe_0 + tr_y_xxxxz_zzzz[i] * pa_x[i];
    }

    // Set up 465-480 components of targeted buffer : IG

    auto tr_y_xxxxyy_xxxx = pbuffer.data(idx_dip_ig + 465);

    auto tr_y_xxxxyy_xxxy = pbuffer.data(idx_dip_ig + 466);

    auto tr_y_xxxxyy_xxxz = pbuffer.data(idx_dip_ig + 467);

    auto tr_y_xxxxyy_xxyy = pbuffer.data(idx_dip_ig + 468);

    auto tr_y_xxxxyy_xxyz = pbuffer.data(idx_dip_ig + 469);

    auto tr_y_xxxxyy_xxzz = pbuffer.data(idx_dip_ig + 470);

    auto tr_y_xxxxyy_xyyy = pbuffer.data(idx_dip_ig + 471);

    auto tr_y_xxxxyy_xyyz = pbuffer.data(idx_dip_ig + 472);

    auto tr_y_xxxxyy_xyzz = pbuffer.data(idx_dip_ig + 473);

    auto tr_y_xxxxyy_xzzz = pbuffer.data(idx_dip_ig + 474);

    auto tr_y_xxxxyy_yyyy = pbuffer.data(idx_dip_ig + 475);

    auto tr_y_xxxxyy_yyyz = pbuffer.data(idx_dip_ig + 476);

    auto tr_y_xxxxyy_yyzz = pbuffer.data(idx_dip_ig + 477);

    auto tr_y_xxxxyy_yzzz = pbuffer.data(idx_dip_ig + 478);

    auto tr_y_xxxxyy_zzzz = pbuffer.data(idx_dip_ig + 479);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxxxyy_xxxx, \
                             tr_y_xxxxyy_xxxy, \
                             tr_y_xxxxyy_xxxz, \
                             tr_y_xxxxyy_xxyy, \
                             tr_y_xxxxyy_xxyz, \
                             tr_y_xxxxyy_xxzz, \
                             tr_y_xxxxyy_xyyy, \
                             tr_y_xxxxyy_xyyz, \
                             tr_y_xxxxyy_xyzz, \
                             tr_y_xxxxyy_xzzz, \
                             tr_y_xxxxyy_yyyy, \
                             tr_y_xxxxyy_yyyz, \
                             tr_y_xxxxyy_yyzz, \
                             tr_y_xxxxyy_yzzz, \
                             tr_y_xxxxyy_zzzz, \
                             tr_y_xxxyy_xxx,   \
                             tr_y_xxxyy_xxxx,  \
                             tr_y_xxxyy_xxxy,  \
                             tr_y_xxxyy_xxxz,  \
                             tr_y_xxxyy_xxy,   \
                             tr_y_xxxyy_xxyy,  \
                             tr_y_xxxyy_xxyz,  \
                             tr_y_xxxyy_xxz,   \
                             tr_y_xxxyy_xxzz,  \
                             tr_y_xxxyy_xyy,   \
                             tr_y_xxxyy_xyyy,  \
                             tr_y_xxxyy_xyyz,  \
                             tr_y_xxxyy_xyz,   \
                             tr_y_xxxyy_xyzz,  \
                             tr_y_xxxyy_xzz,   \
                             tr_y_xxxyy_xzzz,  \
                             tr_y_xxxyy_yyy,   \
                             tr_y_xxxyy_yyyy,  \
                             tr_y_xxxyy_yyyz,  \
                             tr_y_xxxyy_yyz,   \
                             tr_y_xxxyy_yyzz,  \
                             tr_y_xxxyy_yzz,   \
                             tr_y_xxxyy_yzzz,  \
                             tr_y_xxxyy_zzz,   \
                             tr_y_xxxyy_zzzz,  \
                             tr_y_xxyy_xxxx,   \
                             tr_y_xxyy_xxxy,   \
                             tr_y_xxyy_xxxz,   \
                             tr_y_xxyy_xxyy,   \
                             tr_y_xxyy_xxyz,   \
                             tr_y_xxyy_xxzz,   \
                             tr_y_xxyy_xyyy,   \
                             tr_y_xxyy_xyyz,   \
                             tr_y_xxyy_xyzz,   \
                             tr_y_xxyy_xzzz,   \
                             tr_y_xxyy_yyyy,   \
                             tr_y_xxyy_yyyz,   \
                             tr_y_xxyy_yyzz,   \
                             tr_y_xxyy_yzzz,   \
                             tr_y_xxyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyy_xxxx[i] = 3.0 * tr_y_xxyy_xxxx[i] * fe_0 + 4.0 * tr_y_xxxyy_xxx[i] * fe_0 + tr_y_xxxyy_xxxx[i] * pa_x[i];

        tr_y_xxxxyy_xxxy[i] = 3.0 * tr_y_xxyy_xxxy[i] * fe_0 + 3.0 * tr_y_xxxyy_xxy[i] * fe_0 + tr_y_xxxyy_xxxy[i] * pa_x[i];

        tr_y_xxxxyy_xxxz[i] = 3.0 * tr_y_xxyy_xxxz[i] * fe_0 + 3.0 * tr_y_xxxyy_xxz[i] * fe_0 + tr_y_xxxyy_xxxz[i] * pa_x[i];

        tr_y_xxxxyy_xxyy[i] = 3.0 * tr_y_xxyy_xxyy[i] * fe_0 + 2.0 * tr_y_xxxyy_xyy[i] * fe_0 + tr_y_xxxyy_xxyy[i] * pa_x[i];

        tr_y_xxxxyy_xxyz[i] = 3.0 * tr_y_xxyy_xxyz[i] * fe_0 + 2.0 * tr_y_xxxyy_xyz[i] * fe_0 + tr_y_xxxyy_xxyz[i] * pa_x[i];

        tr_y_xxxxyy_xxzz[i] = 3.0 * tr_y_xxyy_xxzz[i] * fe_0 + 2.0 * tr_y_xxxyy_xzz[i] * fe_0 + tr_y_xxxyy_xxzz[i] * pa_x[i];

        tr_y_xxxxyy_xyyy[i] = 3.0 * tr_y_xxyy_xyyy[i] * fe_0 + tr_y_xxxyy_yyy[i] * fe_0 + tr_y_xxxyy_xyyy[i] * pa_x[i];

        tr_y_xxxxyy_xyyz[i] = 3.0 * tr_y_xxyy_xyyz[i] * fe_0 + tr_y_xxxyy_yyz[i] * fe_0 + tr_y_xxxyy_xyyz[i] * pa_x[i];

        tr_y_xxxxyy_xyzz[i] = 3.0 * tr_y_xxyy_xyzz[i] * fe_0 + tr_y_xxxyy_yzz[i] * fe_0 + tr_y_xxxyy_xyzz[i] * pa_x[i];

        tr_y_xxxxyy_xzzz[i] = 3.0 * tr_y_xxyy_xzzz[i] * fe_0 + tr_y_xxxyy_zzz[i] * fe_0 + tr_y_xxxyy_xzzz[i] * pa_x[i];

        tr_y_xxxxyy_yyyy[i] = 3.0 * tr_y_xxyy_yyyy[i] * fe_0 + tr_y_xxxyy_yyyy[i] * pa_x[i];

        tr_y_xxxxyy_yyyz[i] = 3.0 * tr_y_xxyy_yyyz[i] * fe_0 + tr_y_xxxyy_yyyz[i] * pa_x[i];

        tr_y_xxxxyy_yyzz[i] = 3.0 * tr_y_xxyy_yyzz[i] * fe_0 + tr_y_xxxyy_yyzz[i] * pa_x[i];

        tr_y_xxxxyy_yzzz[i] = 3.0 * tr_y_xxyy_yzzz[i] * fe_0 + tr_y_xxxyy_yzzz[i] * pa_x[i];

        tr_y_xxxxyy_zzzz[i] = 3.0 * tr_y_xxyy_zzzz[i] * fe_0 + tr_y_xxxyy_zzzz[i] * pa_x[i];
    }

    // Set up 480-495 components of targeted buffer : IG

    auto tr_y_xxxxyz_xxxx = pbuffer.data(idx_dip_ig + 480);

    auto tr_y_xxxxyz_xxxy = pbuffer.data(idx_dip_ig + 481);

    auto tr_y_xxxxyz_xxxz = pbuffer.data(idx_dip_ig + 482);

    auto tr_y_xxxxyz_xxyy = pbuffer.data(idx_dip_ig + 483);

    auto tr_y_xxxxyz_xxyz = pbuffer.data(idx_dip_ig + 484);

    auto tr_y_xxxxyz_xxzz = pbuffer.data(idx_dip_ig + 485);

    auto tr_y_xxxxyz_xyyy = pbuffer.data(idx_dip_ig + 486);

    auto tr_y_xxxxyz_xyyz = pbuffer.data(idx_dip_ig + 487);

    auto tr_y_xxxxyz_xyzz = pbuffer.data(idx_dip_ig + 488);

    auto tr_y_xxxxyz_xzzz = pbuffer.data(idx_dip_ig + 489);

    auto tr_y_xxxxyz_yyyy = pbuffer.data(idx_dip_ig + 490);

    auto tr_y_xxxxyz_yyyz = pbuffer.data(idx_dip_ig + 491);

    auto tr_y_xxxxyz_yyzz = pbuffer.data(idx_dip_ig + 492);

    auto tr_y_xxxxyz_yzzz = pbuffer.data(idx_dip_ig + 493);

    auto tr_y_xxxxyz_zzzz = pbuffer.data(idx_dip_ig + 494);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_y_xxxxy_xxxx,  \
                             tr_y_xxxxy_xxxy,  \
                             tr_y_xxxxy_xxy,   \
                             tr_y_xxxxy_xxyy,  \
                             tr_y_xxxxy_xxyz,  \
                             tr_y_xxxxy_xyy,   \
                             tr_y_xxxxy_xyyy,  \
                             tr_y_xxxxy_xyyz,  \
                             tr_y_xxxxy_xyz,   \
                             tr_y_xxxxy_xyzz,  \
                             tr_y_xxxxy_yyyy,  \
                             tr_y_xxxxyz_xxxx, \
                             tr_y_xxxxyz_xxxy, \
                             tr_y_xxxxyz_xxxz, \
                             tr_y_xxxxyz_xxyy, \
                             tr_y_xxxxyz_xxyz, \
                             tr_y_xxxxyz_xxzz, \
                             tr_y_xxxxyz_xyyy, \
                             tr_y_xxxxyz_xyyz, \
                             tr_y_xxxxyz_xyzz, \
                             tr_y_xxxxyz_xzzz, \
                             tr_y_xxxxyz_yyyy, \
                             tr_y_xxxxyz_yyyz, \
                             tr_y_xxxxyz_yyzz, \
                             tr_y_xxxxyz_yzzz, \
                             tr_y_xxxxyz_zzzz, \
                             tr_y_xxxxz_xxxz,  \
                             tr_y_xxxxz_xxzz,  \
                             tr_y_xxxxz_xzzz,  \
                             tr_y_xxxyz_yyyz,  \
                             tr_y_xxxyz_yyzz,  \
                             tr_y_xxxyz_yzzz,  \
                             tr_y_xxxyz_zzzz,  \
                             tr_y_xxyz_yyyz,   \
                             tr_y_xxyz_yyzz,   \
                             tr_y_xxyz_yzzz,   \
                             tr_y_xxyz_zzzz,   \
                             ts_xxxxz_xxxz,    \
                             ts_xxxxz_xxzz,    \
                             ts_xxxxz_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyz_xxxx[i] = tr_y_xxxxy_xxxx[i] * pa_z[i];

        tr_y_xxxxyz_xxxy[i] = tr_y_xxxxy_xxxy[i] * pa_z[i];

        tr_y_xxxxyz_xxxz[i] = ts_xxxxz_xxxz[i] * fe_0 + tr_y_xxxxz_xxxz[i] * pa_y[i];

        tr_y_xxxxyz_xxyy[i] = tr_y_xxxxy_xxyy[i] * pa_z[i];

        tr_y_xxxxyz_xxyz[i] = tr_y_xxxxy_xxy[i] * fe_0 + tr_y_xxxxy_xxyz[i] * pa_z[i];

        tr_y_xxxxyz_xxzz[i] = ts_xxxxz_xxzz[i] * fe_0 + tr_y_xxxxz_xxzz[i] * pa_y[i];

        tr_y_xxxxyz_xyyy[i] = tr_y_xxxxy_xyyy[i] * pa_z[i];

        tr_y_xxxxyz_xyyz[i] = tr_y_xxxxy_xyy[i] * fe_0 + tr_y_xxxxy_xyyz[i] * pa_z[i];

        tr_y_xxxxyz_xyzz[i] = 2.0 * tr_y_xxxxy_xyz[i] * fe_0 + tr_y_xxxxy_xyzz[i] * pa_z[i];

        tr_y_xxxxyz_xzzz[i] = ts_xxxxz_xzzz[i] * fe_0 + tr_y_xxxxz_xzzz[i] * pa_y[i];

        tr_y_xxxxyz_yyyy[i] = tr_y_xxxxy_yyyy[i] * pa_z[i];

        tr_y_xxxxyz_yyyz[i] = 3.0 * tr_y_xxyz_yyyz[i] * fe_0 + tr_y_xxxyz_yyyz[i] * pa_x[i];

        tr_y_xxxxyz_yyzz[i] = 3.0 * tr_y_xxyz_yyzz[i] * fe_0 + tr_y_xxxyz_yyzz[i] * pa_x[i];

        tr_y_xxxxyz_yzzz[i] = 3.0 * tr_y_xxyz_yzzz[i] * fe_0 + tr_y_xxxyz_yzzz[i] * pa_x[i];

        tr_y_xxxxyz_zzzz[i] = 3.0 * tr_y_xxyz_zzzz[i] * fe_0 + tr_y_xxxyz_zzzz[i] * pa_x[i];
    }

    // Set up 495-510 components of targeted buffer : IG

    auto tr_y_xxxxzz_xxxx = pbuffer.data(idx_dip_ig + 495);

    auto tr_y_xxxxzz_xxxy = pbuffer.data(idx_dip_ig + 496);

    auto tr_y_xxxxzz_xxxz = pbuffer.data(idx_dip_ig + 497);

    auto tr_y_xxxxzz_xxyy = pbuffer.data(idx_dip_ig + 498);

    auto tr_y_xxxxzz_xxyz = pbuffer.data(idx_dip_ig + 499);

    auto tr_y_xxxxzz_xxzz = pbuffer.data(idx_dip_ig + 500);

    auto tr_y_xxxxzz_xyyy = pbuffer.data(idx_dip_ig + 501);

    auto tr_y_xxxxzz_xyyz = pbuffer.data(idx_dip_ig + 502);

    auto tr_y_xxxxzz_xyzz = pbuffer.data(idx_dip_ig + 503);

    auto tr_y_xxxxzz_xzzz = pbuffer.data(idx_dip_ig + 504);

    auto tr_y_xxxxzz_yyyy = pbuffer.data(idx_dip_ig + 505);

    auto tr_y_xxxxzz_yyyz = pbuffer.data(idx_dip_ig + 506);

    auto tr_y_xxxxzz_yyzz = pbuffer.data(idx_dip_ig + 507);

    auto tr_y_xxxxzz_yzzz = pbuffer.data(idx_dip_ig + 508);

    auto tr_y_xxxxzz_zzzz = pbuffer.data(idx_dip_ig + 509);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxxx_xxxx,   \
                             tr_y_xxxx_xxxy,   \
                             tr_y_xxxx_xxyy,   \
                             tr_y_xxxx_xyyy,   \
                             tr_y_xxxxz_xxxx,  \
                             tr_y_xxxxz_xxxy,  \
                             tr_y_xxxxz_xxyy,  \
                             tr_y_xxxxz_xyyy,  \
                             tr_y_xxxxzz_xxxx, \
                             tr_y_xxxxzz_xxxy, \
                             tr_y_xxxxzz_xxxz, \
                             tr_y_xxxxzz_xxyy, \
                             tr_y_xxxxzz_xxyz, \
                             tr_y_xxxxzz_xxzz, \
                             tr_y_xxxxzz_xyyy, \
                             tr_y_xxxxzz_xyyz, \
                             tr_y_xxxxzz_xyzz, \
                             tr_y_xxxxzz_xzzz, \
                             tr_y_xxxxzz_yyyy, \
                             tr_y_xxxxzz_yyyz, \
                             tr_y_xxxxzz_yyzz, \
                             tr_y_xxxxzz_yzzz, \
                             tr_y_xxxxzz_zzzz, \
                             tr_y_xxxzz_xxxz,  \
                             tr_y_xxxzz_xxyz,  \
                             tr_y_xxxzz_xxz,   \
                             tr_y_xxxzz_xxzz,  \
                             tr_y_xxxzz_xyyz,  \
                             tr_y_xxxzz_xyz,   \
                             tr_y_xxxzz_xyzz,  \
                             tr_y_xxxzz_xzz,   \
                             tr_y_xxxzz_xzzz,  \
                             tr_y_xxxzz_yyyy,  \
                             tr_y_xxxzz_yyyz,  \
                             tr_y_xxxzz_yyz,   \
                             tr_y_xxxzz_yyzz,  \
                             tr_y_xxxzz_yzz,   \
                             tr_y_xxxzz_yzzz,  \
                             tr_y_xxxzz_zzz,   \
                             tr_y_xxxzz_zzzz,  \
                             tr_y_xxzz_xxxz,   \
                             tr_y_xxzz_xxyz,   \
                             tr_y_xxzz_xxzz,   \
                             tr_y_xxzz_xyyz,   \
                             tr_y_xxzz_xyzz,   \
                             tr_y_xxzz_xzzz,   \
                             tr_y_xxzz_yyyy,   \
                             tr_y_xxzz_yyyz,   \
                             tr_y_xxzz_yyzz,   \
                             tr_y_xxzz_yzzz,   \
                             tr_y_xxzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxzz_xxxx[i] = tr_y_xxxx_xxxx[i] * fe_0 + tr_y_xxxxz_xxxx[i] * pa_z[i];

        tr_y_xxxxzz_xxxy[i] = tr_y_xxxx_xxxy[i] * fe_0 + tr_y_xxxxz_xxxy[i] * pa_z[i];

        tr_y_xxxxzz_xxxz[i] = 3.0 * tr_y_xxzz_xxxz[i] * fe_0 + 3.0 * tr_y_xxxzz_xxz[i] * fe_0 + tr_y_xxxzz_xxxz[i] * pa_x[i];

        tr_y_xxxxzz_xxyy[i] = tr_y_xxxx_xxyy[i] * fe_0 + tr_y_xxxxz_xxyy[i] * pa_z[i];

        tr_y_xxxxzz_xxyz[i] = 3.0 * tr_y_xxzz_xxyz[i] * fe_0 + 2.0 * tr_y_xxxzz_xyz[i] * fe_0 + tr_y_xxxzz_xxyz[i] * pa_x[i];

        tr_y_xxxxzz_xxzz[i] = 3.0 * tr_y_xxzz_xxzz[i] * fe_0 + 2.0 * tr_y_xxxzz_xzz[i] * fe_0 + tr_y_xxxzz_xxzz[i] * pa_x[i];

        tr_y_xxxxzz_xyyy[i] = tr_y_xxxx_xyyy[i] * fe_0 + tr_y_xxxxz_xyyy[i] * pa_z[i];

        tr_y_xxxxzz_xyyz[i] = 3.0 * tr_y_xxzz_xyyz[i] * fe_0 + tr_y_xxxzz_yyz[i] * fe_0 + tr_y_xxxzz_xyyz[i] * pa_x[i];

        tr_y_xxxxzz_xyzz[i] = 3.0 * tr_y_xxzz_xyzz[i] * fe_0 + tr_y_xxxzz_yzz[i] * fe_0 + tr_y_xxxzz_xyzz[i] * pa_x[i];

        tr_y_xxxxzz_xzzz[i] = 3.0 * tr_y_xxzz_xzzz[i] * fe_0 + tr_y_xxxzz_zzz[i] * fe_0 + tr_y_xxxzz_xzzz[i] * pa_x[i];

        tr_y_xxxxzz_yyyy[i] = 3.0 * tr_y_xxzz_yyyy[i] * fe_0 + tr_y_xxxzz_yyyy[i] * pa_x[i];

        tr_y_xxxxzz_yyyz[i] = 3.0 * tr_y_xxzz_yyyz[i] * fe_0 + tr_y_xxxzz_yyyz[i] * pa_x[i];

        tr_y_xxxxzz_yyzz[i] = 3.0 * tr_y_xxzz_yyzz[i] * fe_0 + tr_y_xxxzz_yyzz[i] * pa_x[i];

        tr_y_xxxxzz_yzzz[i] = 3.0 * tr_y_xxzz_yzzz[i] * fe_0 + tr_y_xxxzz_yzzz[i] * pa_x[i];

        tr_y_xxxxzz_zzzz[i] = 3.0 * tr_y_xxzz_zzzz[i] * fe_0 + tr_y_xxxzz_zzzz[i] * pa_x[i];
    }

    // Set up 510-525 components of targeted buffer : IG

    auto tr_y_xxxyyy_xxxx = pbuffer.data(idx_dip_ig + 510);

    auto tr_y_xxxyyy_xxxy = pbuffer.data(idx_dip_ig + 511);

    auto tr_y_xxxyyy_xxxz = pbuffer.data(idx_dip_ig + 512);

    auto tr_y_xxxyyy_xxyy = pbuffer.data(idx_dip_ig + 513);

    auto tr_y_xxxyyy_xxyz = pbuffer.data(idx_dip_ig + 514);

    auto tr_y_xxxyyy_xxzz = pbuffer.data(idx_dip_ig + 515);

    auto tr_y_xxxyyy_xyyy = pbuffer.data(idx_dip_ig + 516);

    auto tr_y_xxxyyy_xyyz = pbuffer.data(idx_dip_ig + 517);

    auto tr_y_xxxyyy_xyzz = pbuffer.data(idx_dip_ig + 518);

    auto tr_y_xxxyyy_xzzz = pbuffer.data(idx_dip_ig + 519);

    auto tr_y_xxxyyy_yyyy = pbuffer.data(idx_dip_ig + 520);

    auto tr_y_xxxyyy_yyyz = pbuffer.data(idx_dip_ig + 521);

    auto tr_y_xxxyyy_yyzz = pbuffer.data(idx_dip_ig + 522);

    auto tr_y_xxxyyy_yzzz = pbuffer.data(idx_dip_ig + 523);

    auto tr_y_xxxyyy_zzzz = pbuffer.data(idx_dip_ig + 524);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxxyyy_xxxx, \
                             tr_y_xxxyyy_xxxy, \
                             tr_y_xxxyyy_xxxz, \
                             tr_y_xxxyyy_xxyy, \
                             tr_y_xxxyyy_xxyz, \
                             tr_y_xxxyyy_xxzz, \
                             tr_y_xxxyyy_xyyy, \
                             tr_y_xxxyyy_xyyz, \
                             tr_y_xxxyyy_xyzz, \
                             tr_y_xxxyyy_xzzz, \
                             tr_y_xxxyyy_yyyy, \
                             tr_y_xxxyyy_yyyz, \
                             tr_y_xxxyyy_yyzz, \
                             tr_y_xxxyyy_yzzz, \
                             tr_y_xxxyyy_zzzz, \
                             tr_y_xxyyy_xxx,   \
                             tr_y_xxyyy_xxxx,  \
                             tr_y_xxyyy_xxxy,  \
                             tr_y_xxyyy_xxxz,  \
                             tr_y_xxyyy_xxy,   \
                             tr_y_xxyyy_xxyy,  \
                             tr_y_xxyyy_xxyz,  \
                             tr_y_xxyyy_xxz,   \
                             tr_y_xxyyy_xxzz,  \
                             tr_y_xxyyy_xyy,   \
                             tr_y_xxyyy_xyyy,  \
                             tr_y_xxyyy_xyyz,  \
                             tr_y_xxyyy_xyz,   \
                             tr_y_xxyyy_xyzz,  \
                             tr_y_xxyyy_xzz,   \
                             tr_y_xxyyy_xzzz,  \
                             tr_y_xxyyy_yyy,   \
                             tr_y_xxyyy_yyyy,  \
                             tr_y_xxyyy_yyyz,  \
                             tr_y_xxyyy_yyz,   \
                             tr_y_xxyyy_yyzz,  \
                             tr_y_xxyyy_yzz,   \
                             tr_y_xxyyy_yzzz,  \
                             tr_y_xxyyy_zzz,   \
                             tr_y_xxyyy_zzzz,  \
                             tr_y_xyyy_xxxx,   \
                             tr_y_xyyy_xxxy,   \
                             tr_y_xyyy_xxxz,   \
                             tr_y_xyyy_xxyy,   \
                             tr_y_xyyy_xxyz,   \
                             tr_y_xyyy_xxzz,   \
                             tr_y_xyyy_xyyy,   \
                             tr_y_xyyy_xyyz,   \
                             tr_y_xyyy_xyzz,   \
                             tr_y_xyyy_xzzz,   \
                             tr_y_xyyy_yyyy,   \
                             tr_y_xyyy_yyyz,   \
                             tr_y_xyyy_yyzz,   \
                             tr_y_xyyy_yzzz,   \
                             tr_y_xyyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyy_xxxx[i] = 2.0 * tr_y_xyyy_xxxx[i] * fe_0 + 4.0 * tr_y_xxyyy_xxx[i] * fe_0 + tr_y_xxyyy_xxxx[i] * pa_x[i];

        tr_y_xxxyyy_xxxy[i] = 2.0 * tr_y_xyyy_xxxy[i] * fe_0 + 3.0 * tr_y_xxyyy_xxy[i] * fe_0 + tr_y_xxyyy_xxxy[i] * pa_x[i];

        tr_y_xxxyyy_xxxz[i] = 2.0 * tr_y_xyyy_xxxz[i] * fe_0 + 3.0 * tr_y_xxyyy_xxz[i] * fe_0 + tr_y_xxyyy_xxxz[i] * pa_x[i];

        tr_y_xxxyyy_xxyy[i] = 2.0 * tr_y_xyyy_xxyy[i] * fe_0 + 2.0 * tr_y_xxyyy_xyy[i] * fe_0 + tr_y_xxyyy_xxyy[i] * pa_x[i];

        tr_y_xxxyyy_xxyz[i] = 2.0 * tr_y_xyyy_xxyz[i] * fe_0 + 2.0 * tr_y_xxyyy_xyz[i] * fe_0 + tr_y_xxyyy_xxyz[i] * pa_x[i];

        tr_y_xxxyyy_xxzz[i] = 2.0 * tr_y_xyyy_xxzz[i] * fe_0 + 2.0 * tr_y_xxyyy_xzz[i] * fe_0 + tr_y_xxyyy_xxzz[i] * pa_x[i];

        tr_y_xxxyyy_xyyy[i] = 2.0 * tr_y_xyyy_xyyy[i] * fe_0 + tr_y_xxyyy_yyy[i] * fe_0 + tr_y_xxyyy_xyyy[i] * pa_x[i];

        tr_y_xxxyyy_xyyz[i] = 2.0 * tr_y_xyyy_xyyz[i] * fe_0 + tr_y_xxyyy_yyz[i] * fe_0 + tr_y_xxyyy_xyyz[i] * pa_x[i];

        tr_y_xxxyyy_xyzz[i] = 2.0 * tr_y_xyyy_xyzz[i] * fe_0 + tr_y_xxyyy_yzz[i] * fe_0 + tr_y_xxyyy_xyzz[i] * pa_x[i];

        tr_y_xxxyyy_xzzz[i] = 2.0 * tr_y_xyyy_xzzz[i] * fe_0 + tr_y_xxyyy_zzz[i] * fe_0 + tr_y_xxyyy_xzzz[i] * pa_x[i];

        tr_y_xxxyyy_yyyy[i] = 2.0 * tr_y_xyyy_yyyy[i] * fe_0 + tr_y_xxyyy_yyyy[i] * pa_x[i];

        tr_y_xxxyyy_yyyz[i] = 2.0 * tr_y_xyyy_yyyz[i] * fe_0 + tr_y_xxyyy_yyyz[i] * pa_x[i];

        tr_y_xxxyyy_yyzz[i] = 2.0 * tr_y_xyyy_yyzz[i] * fe_0 + tr_y_xxyyy_yyzz[i] * pa_x[i];

        tr_y_xxxyyy_yzzz[i] = 2.0 * tr_y_xyyy_yzzz[i] * fe_0 + tr_y_xxyyy_yzzz[i] * pa_x[i];

        tr_y_xxxyyy_zzzz[i] = 2.0 * tr_y_xyyy_zzzz[i] * fe_0 + tr_y_xxyyy_zzzz[i] * pa_x[i];
    }

    // Set up 525-540 components of targeted buffer : IG

    auto tr_y_xxxyyz_xxxx = pbuffer.data(idx_dip_ig + 525);

    auto tr_y_xxxyyz_xxxy = pbuffer.data(idx_dip_ig + 526);

    auto tr_y_xxxyyz_xxxz = pbuffer.data(idx_dip_ig + 527);

    auto tr_y_xxxyyz_xxyy = pbuffer.data(idx_dip_ig + 528);

    auto tr_y_xxxyyz_xxyz = pbuffer.data(idx_dip_ig + 529);

    auto tr_y_xxxyyz_xxzz = pbuffer.data(idx_dip_ig + 530);

    auto tr_y_xxxyyz_xyyy = pbuffer.data(idx_dip_ig + 531);

    auto tr_y_xxxyyz_xyyz = pbuffer.data(idx_dip_ig + 532);

    auto tr_y_xxxyyz_xyzz = pbuffer.data(idx_dip_ig + 533);

    auto tr_y_xxxyyz_xzzz = pbuffer.data(idx_dip_ig + 534);

    auto tr_y_xxxyyz_yyyy = pbuffer.data(idx_dip_ig + 535);

    auto tr_y_xxxyyz_yyyz = pbuffer.data(idx_dip_ig + 536);

    auto tr_y_xxxyyz_yyzz = pbuffer.data(idx_dip_ig + 537);

    auto tr_y_xxxyyz_yzzz = pbuffer.data(idx_dip_ig + 538);

    auto tr_y_xxxyyz_zzzz = pbuffer.data(idx_dip_ig + 539);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxxyy_xxx,   \
                             tr_y_xxxyy_xxxx,  \
                             tr_y_xxxyy_xxxy,  \
                             tr_y_xxxyy_xxxz,  \
                             tr_y_xxxyy_xxy,   \
                             tr_y_xxxyy_xxyy,  \
                             tr_y_xxxyy_xxyz,  \
                             tr_y_xxxyy_xxz,   \
                             tr_y_xxxyy_xxzz,  \
                             tr_y_xxxyy_xyy,   \
                             tr_y_xxxyy_xyyy,  \
                             tr_y_xxxyy_xyyz,  \
                             tr_y_xxxyy_xyz,   \
                             tr_y_xxxyy_xyzz,  \
                             tr_y_xxxyy_xzz,   \
                             tr_y_xxxyy_xzzz,  \
                             tr_y_xxxyy_yyyy,  \
                             tr_y_xxxyyz_xxxx, \
                             tr_y_xxxyyz_xxxy, \
                             tr_y_xxxyyz_xxxz, \
                             tr_y_xxxyyz_xxyy, \
                             tr_y_xxxyyz_xxyz, \
                             tr_y_xxxyyz_xxzz, \
                             tr_y_xxxyyz_xyyy, \
                             tr_y_xxxyyz_xyyz, \
                             tr_y_xxxyyz_xyzz, \
                             tr_y_xxxyyz_xzzz, \
                             tr_y_xxxyyz_yyyy, \
                             tr_y_xxxyyz_yyyz, \
                             tr_y_xxxyyz_yyzz, \
                             tr_y_xxxyyz_yzzz, \
                             tr_y_xxxyyz_zzzz, \
                             tr_y_xxyyz_yyyz,  \
                             tr_y_xxyyz_yyzz,  \
                             tr_y_xxyyz_yzzz,  \
                             tr_y_xxyyz_zzzz,  \
                             tr_y_xyyz_yyyz,   \
                             tr_y_xyyz_yyzz,   \
                             tr_y_xyyz_yzzz,   \
                             tr_y_xyyz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyz_xxxx[i] = tr_y_xxxyy_xxxx[i] * pa_z[i];

        tr_y_xxxyyz_xxxy[i] = tr_y_xxxyy_xxxy[i] * pa_z[i];

        tr_y_xxxyyz_xxxz[i] = tr_y_xxxyy_xxx[i] * fe_0 + tr_y_xxxyy_xxxz[i] * pa_z[i];

        tr_y_xxxyyz_xxyy[i] = tr_y_xxxyy_xxyy[i] * pa_z[i];

        tr_y_xxxyyz_xxyz[i] = tr_y_xxxyy_xxy[i] * fe_0 + tr_y_xxxyy_xxyz[i] * pa_z[i];

        tr_y_xxxyyz_xxzz[i] = 2.0 * tr_y_xxxyy_xxz[i] * fe_0 + tr_y_xxxyy_xxzz[i] * pa_z[i];

        tr_y_xxxyyz_xyyy[i] = tr_y_xxxyy_xyyy[i] * pa_z[i];

        tr_y_xxxyyz_xyyz[i] = tr_y_xxxyy_xyy[i] * fe_0 + tr_y_xxxyy_xyyz[i] * pa_z[i];

        tr_y_xxxyyz_xyzz[i] = 2.0 * tr_y_xxxyy_xyz[i] * fe_0 + tr_y_xxxyy_xyzz[i] * pa_z[i];

        tr_y_xxxyyz_xzzz[i] = 3.0 * tr_y_xxxyy_xzz[i] * fe_0 + tr_y_xxxyy_xzzz[i] * pa_z[i];

        tr_y_xxxyyz_yyyy[i] = tr_y_xxxyy_yyyy[i] * pa_z[i];

        tr_y_xxxyyz_yyyz[i] = 2.0 * tr_y_xyyz_yyyz[i] * fe_0 + tr_y_xxyyz_yyyz[i] * pa_x[i];

        tr_y_xxxyyz_yyzz[i] = 2.0 * tr_y_xyyz_yyzz[i] * fe_0 + tr_y_xxyyz_yyzz[i] * pa_x[i];

        tr_y_xxxyyz_yzzz[i] = 2.0 * tr_y_xyyz_yzzz[i] * fe_0 + tr_y_xxyyz_yzzz[i] * pa_x[i];

        tr_y_xxxyyz_zzzz[i] = 2.0 * tr_y_xyyz_zzzz[i] * fe_0 + tr_y_xxyyz_zzzz[i] * pa_x[i];
    }

    // Set up 540-555 components of targeted buffer : IG

    auto tr_y_xxxyzz_xxxx = pbuffer.data(idx_dip_ig + 540);

    auto tr_y_xxxyzz_xxxy = pbuffer.data(idx_dip_ig + 541);

    auto tr_y_xxxyzz_xxxz = pbuffer.data(idx_dip_ig + 542);

    auto tr_y_xxxyzz_xxyy = pbuffer.data(idx_dip_ig + 543);

    auto tr_y_xxxyzz_xxyz = pbuffer.data(idx_dip_ig + 544);

    auto tr_y_xxxyzz_xxzz = pbuffer.data(idx_dip_ig + 545);

    auto tr_y_xxxyzz_xyyy = pbuffer.data(idx_dip_ig + 546);

    auto tr_y_xxxyzz_xyyz = pbuffer.data(idx_dip_ig + 547);

    auto tr_y_xxxyzz_xyzz = pbuffer.data(idx_dip_ig + 548);

    auto tr_y_xxxyzz_xzzz = pbuffer.data(idx_dip_ig + 549);

    auto tr_y_xxxyzz_yyyy = pbuffer.data(idx_dip_ig + 550);

    auto tr_y_xxxyzz_yyyz = pbuffer.data(idx_dip_ig + 551);

    auto tr_y_xxxyzz_yyzz = pbuffer.data(idx_dip_ig + 552);

    auto tr_y_xxxyzz_yzzz = pbuffer.data(idx_dip_ig + 553);

    auto tr_y_xxxyzz_zzzz = pbuffer.data(idx_dip_ig + 554);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_y_xxxy_xxxy,   \
                             tr_y_xxxy_xxyy,   \
                             tr_y_xxxy_xyyy,   \
                             tr_y_xxxyz_xxxy,  \
                             tr_y_xxxyz_xxyy,  \
                             tr_y_xxxyz_xyyy,  \
                             tr_y_xxxyzz_xxxx, \
                             tr_y_xxxyzz_xxxy, \
                             tr_y_xxxyzz_xxxz, \
                             tr_y_xxxyzz_xxyy, \
                             tr_y_xxxyzz_xxyz, \
                             tr_y_xxxyzz_xxzz, \
                             tr_y_xxxyzz_xyyy, \
                             tr_y_xxxyzz_xyyz, \
                             tr_y_xxxyzz_xyzz, \
                             tr_y_xxxyzz_xzzz, \
                             tr_y_xxxyzz_yyyy, \
                             tr_y_xxxyzz_yyyz, \
                             tr_y_xxxyzz_yyzz, \
                             tr_y_xxxyzz_yzzz, \
                             tr_y_xxxyzz_zzzz, \
                             tr_y_xxxzz_xxxx,  \
                             tr_y_xxxzz_xxxz,  \
                             tr_y_xxxzz_xxzz,  \
                             tr_y_xxxzz_xzzz,  \
                             tr_y_xxyzz_xxyz,  \
                             tr_y_xxyzz_xyyz,  \
                             tr_y_xxyzz_xyz,   \
                             tr_y_xxyzz_xyzz,  \
                             tr_y_xxyzz_yyyy,  \
                             tr_y_xxyzz_yyyz,  \
                             tr_y_xxyzz_yyz,   \
                             tr_y_xxyzz_yyzz,  \
                             tr_y_xxyzz_yzz,   \
                             tr_y_xxyzz_yzzz,  \
                             tr_y_xxyzz_zzzz,  \
                             tr_y_xyzz_xxyz,   \
                             tr_y_xyzz_xyyz,   \
                             tr_y_xyzz_xyzz,   \
                             tr_y_xyzz_yyyy,   \
                             tr_y_xyzz_yyyz,   \
                             tr_y_xyzz_yyzz,   \
                             tr_y_xyzz_yzzz,   \
                             tr_y_xyzz_zzzz,   \
                             ts_xxxzz_xxxx,    \
                             ts_xxxzz_xxxz,    \
                             ts_xxxzz_xxzz,    \
                             ts_xxxzz_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyzz_xxxx[i] = ts_xxxzz_xxxx[i] * fe_0 + tr_y_xxxzz_xxxx[i] * pa_y[i];

        tr_y_xxxyzz_xxxy[i] = tr_y_xxxy_xxxy[i] * fe_0 + tr_y_xxxyz_xxxy[i] * pa_z[i];

        tr_y_xxxyzz_xxxz[i] = ts_xxxzz_xxxz[i] * fe_0 + tr_y_xxxzz_xxxz[i] * pa_y[i];

        tr_y_xxxyzz_xxyy[i] = tr_y_xxxy_xxyy[i] * fe_0 + tr_y_xxxyz_xxyy[i] * pa_z[i];

        tr_y_xxxyzz_xxyz[i] = 2.0 * tr_y_xyzz_xxyz[i] * fe_0 + 2.0 * tr_y_xxyzz_xyz[i] * fe_0 + tr_y_xxyzz_xxyz[i] * pa_x[i];

        tr_y_xxxyzz_xxzz[i] = ts_xxxzz_xxzz[i] * fe_0 + tr_y_xxxzz_xxzz[i] * pa_y[i];

        tr_y_xxxyzz_xyyy[i] = tr_y_xxxy_xyyy[i] * fe_0 + tr_y_xxxyz_xyyy[i] * pa_z[i];

        tr_y_xxxyzz_xyyz[i] = 2.0 * tr_y_xyzz_xyyz[i] * fe_0 + tr_y_xxyzz_yyz[i] * fe_0 + tr_y_xxyzz_xyyz[i] * pa_x[i];

        tr_y_xxxyzz_xyzz[i] = 2.0 * tr_y_xyzz_xyzz[i] * fe_0 + tr_y_xxyzz_yzz[i] * fe_0 + tr_y_xxyzz_xyzz[i] * pa_x[i];

        tr_y_xxxyzz_xzzz[i] = ts_xxxzz_xzzz[i] * fe_0 + tr_y_xxxzz_xzzz[i] * pa_y[i];

        tr_y_xxxyzz_yyyy[i] = 2.0 * tr_y_xyzz_yyyy[i] * fe_0 + tr_y_xxyzz_yyyy[i] * pa_x[i];

        tr_y_xxxyzz_yyyz[i] = 2.0 * tr_y_xyzz_yyyz[i] * fe_0 + tr_y_xxyzz_yyyz[i] * pa_x[i];

        tr_y_xxxyzz_yyzz[i] = 2.0 * tr_y_xyzz_yyzz[i] * fe_0 + tr_y_xxyzz_yyzz[i] * pa_x[i];

        tr_y_xxxyzz_yzzz[i] = 2.0 * tr_y_xyzz_yzzz[i] * fe_0 + tr_y_xxyzz_yzzz[i] * pa_x[i];

        tr_y_xxxyzz_zzzz[i] = 2.0 * tr_y_xyzz_zzzz[i] * fe_0 + tr_y_xxyzz_zzzz[i] * pa_x[i];
    }

    // Set up 555-570 components of targeted buffer : IG

    auto tr_y_xxxzzz_xxxx = pbuffer.data(idx_dip_ig + 555);

    auto tr_y_xxxzzz_xxxy = pbuffer.data(idx_dip_ig + 556);

    auto tr_y_xxxzzz_xxxz = pbuffer.data(idx_dip_ig + 557);

    auto tr_y_xxxzzz_xxyy = pbuffer.data(idx_dip_ig + 558);

    auto tr_y_xxxzzz_xxyz = pbuffer.data(idx_dip_ig + 559);

    auto tr_y_xxxzzz_xxzz = pbuffer.data(idx_dip_ig + 560);

    auto tr_y_xxxzzz_xyyy = pbuffer.data(idx_dip_ig + 561);

    auto tr_y_xxxzzz_xyyz = pbuffer.data(idx_dip_ig + 562);

    auto tr_y_xxxzzz_xyzz = pbuffer.data(idx_dip_ig + 563);

    auto tr_y_xxxzzz_xzzz = pbuffer.data(idx_dip_ig + 564);

    auto tr_y_xxxzzz_yyyy = pbuffer.data(idx_dip_ig + 565);

    auto tr_y_xxxzzz_yyyz = pbuffer.data(idx_dip_ig + 566);

    auto tr_y_xxxzzz_yyzz = pbuffer.data(idx_dip_ig + 567);

    auto tr_y_xxxzzz_yzzz = pbuffer.data(idx_dip_ig + 568);

    auto tr_y_xxxzzz_zzzz = pbuffer.data(idx_dip_ig + 569);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxxz_xxxx,   \
                             tr_y_xxxz_xxxy,   \
                             tr_y_xxxz_xxyy,   \
                             tr_y_xxxz_xyyy,   \
                             tr_y_xxxzz_xxxx,  \
                             tr_y_xxxzz_xxxy,  \
                             tr_y_xxxzz_xxyy,  \
                             tr_y_xxxzz_xyyy,  \
                             tr_y_xxxzzz_xxxx, \
                             tr_y_xxxzzz_xxxy, \
                             tr_y_xxxzzz_xxxz, \
                             tr_y_xxxzzz_xxyy, \
                             tr_y_xxxzzz_xxyz, \
                             tr_y_xxxzzz_xxzz, \
                             tr_y_xxxzzz_xyyy, \
                             tr_y_xxxzzz_xyyz, \
                             tr_y_xxxzzz_xyzz, \
                             tr_y_xxxzzz_xzzz, \
                             tr_y_xxxzzz_yyyy, \
                             tr_y_xxxzzz_yyyz, \
                             tr_y_xxxzzz_yyzz, \
                             tr_y_xxxzzz_yzzz, \
                             tr_y_xxxzzz_zzzz, \
                             tr_y_xxzzz_xxxz,  \
                             tr_y_xxzzz_xxyz,  \
                             tr_y_xxzzz_xxz,   \
                             tr_y_xxzzz_xxzz,  \
                             tr_y_xxzzz_xyyz,  \
                             tr_y_xxzzz_xyz,   \
                             tr_y_xxzzz_xyzz,  \
                             tr_y_xxzzz_xzz,   \
                             tr_y_xxzzz_xzzz,  \
                             tr_y_xxzzz_yyyy,  \
                             tr_y_xxzzz_yyyz,  \
                             tr_y_xxzzz_yyz,   \
                             tr_y_xxzzz_yyzz,  \
                             tr_y_xxzzz_yzz,   \
                             tr_y_xxzzz_yzzz,  \
                             tr_y_xxzzz_zzz,   \
                             tr_y_xxzzz_zzzz,  \
                             tr_y_xzzz_xxxz,   \
                             tr_y_xzzz_xxyz,   \
                             tr_y_xzzz_xxzz,   \
                             tr_y_xzzz_xyyz,   \
                             tr_y_xzzz_xyzz,   \
                             tr_y_xzzz_xzzz,   \
                             tr_y_xzzz_yyyy,   \
                             tr_y_xzzz_yyyz,   \
                             tr_y_xzzz_yyzz,   \
                             tr_y_xzzz_yzzz,   \
                             tr_y_xzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzzz_xxxx[i] = 2.0 * tr_y_xxxz_xxxx[i] * fe_0 + tr_y_xxxzz_xxxx[i] * pa_z[i];

        tr_y_xxxzzz_xxxy[i] = 2.0 * tr_y_xxxz_xxxy[i] * fe_0 + tr_y_xxxzz_xxxy[i] * pa_z[i];

        tr_y_xxxzzz_xxxz[i] = 2.0 * tr_y_xzzz_xxxz[i] * fe_0 + 3.0 * tr_y_xxzzz_xxz[i] * fe_0 + tr_y_xxzzz_xxxz[i] * pa_x[i];

        tr_y_xxxzzz_xxyy[i] = 2.0 * tr_y_xxxz_xxyy[i] * fe_0 + tr_y_xxxzz_xxyy[i] * pa_z[i];

        tr_y_xxxzzz_xxyz[i] = 2.0 * tr_y_xzzz_xxyz[i] * fe_0 + 2.0 * tr_y_xxzzz_xyz[i] * fe_0 + tr_y_xxzzz_xxyz[i] * pa_x[i];

        tr_y_xxxzzz_xxzz[i] = 2.0 * tr_y_xzzz_xxzz[i] * fe_0 + 2.0 * tr_y_xxzzz_xzz[i] * fe_0 + tr_y_xxzzz_xxzz[i] * pa_x[i];

        tr_y_xxxzzz_xyyy[i] = 2.0 * tr_y_xxxz_xyyy[i] * fe_0 + tr_y_xxxzz_xyyy[i] * pa_z[i];

        tr_y_xxxzzz_xyyz[i] = 2.0 * tr_y_xzzz_xyyz[i] * fe_0 + tr_y_xxzzz_yyz[i] * fe_0 + tr_y_xxzzz_xyyz[i] * pa_x[i];

        tr_y_xxxzzz_xyzz[i] = 2.0 * tr_y_xzzz_xyzz[i] * fe_0 + tr_y_xxzzz_yzz[i] * fe_0 + tr_y_xxzzz_xyzz[i] * pa_x[i];

        tr_y_xxxzzz_xzzz[i] = 2.0 * tr_y_xzzz_xzzz[i] * fe_0 + tr_y_xxzzz_zzz[i] * fe_0 + tr_y_xxzzz_xzzz[i] * pa_x[i];

        tr_y_xxxzzz_yyyy[i] = 2.0 * tr_y_xzzz_yyyy[i] * fe_0 + tr_y_xxzzz_yyyy[i] * pa_x[i];

        tr_y_xxxzzz_yyyz[i] = 2.0 * tr_y_xzzz_yyyz[i] * fe_0 + tr_y_xxzzz_yyyz[i] * pa_x[i];

        tr_y_xxxzzz_yyzz[i] = 2.0 * tr_y_xzzz_yyzz[i] * fe_0 + tr_y_xxzzz_yyzz[i] * pa_x[i];

        tr_y_xxxzzz_yzzz[i] = 2.0 * tr_y_xzzz_yzzz[i] * fe_0 + tr_y_xxzzz_yzzz[i] * pa_x[i];

        tr_y_xxxzzz_zzzz[i] = 2.0 * tr_y_xzzz_zzzz[i] * fe_0 + tr_y_xxzzz_zzzz[i] * pa_x[i];
    }

    // Set up 570-585 components of targeted buffer : IG

    auto tr_y_xxyyyy_xxxx = pbuffer.data(idx_dip_ig + 570);

    auto tr_y_xxyyyy_xxxy = pbuffer.data(idx_dip_ig + 571);

    auto tr_y_xxyyyy_xxxz = pbuffer.data(idx_dip_ig + 572);

    auto tr_y_xxyyyy_xxyy = pbuffer.data(idx_dip_ig + 573);

    auto tr_y_xxyyyy_xxyz = pbuffer.data(idx_dip_ig + 574);

    auto tr_y_xxyyyy_xxzz = pbuffer.data(idx_dip_ig + 575);

    auto tr_y_xxyyyy_xyyy = pbuffer.data(idx_dip_ig + 576);

    auto tr_y_xxyyyy_xyyz = pbuffer.data(idx_dip_ig + 577);

    auto tr_y_xxyyyy_xyzz = pbuffer.data(idx_dip_ig + 578);

    auto tr_y_xxyyyy_xzzz = pbuffer.data(idx_dip_ig + 579);

    auto tr_y_xxyyyy_yyyy = pbuffer.data(idx_dip_ig + 580);

    auto tr_y_xxyyyy_yyyz = pbuffer.data(idx_dip_ig + 581);

    auto tr_y_xxyyyy_yyzz = pbuffer.data(idx_dip_ig + 582);

    auto tr_y_xxyyyy_yzzz = pbuffer.data(idx_dip_ig + 583);

    auto tr_y_xxyyyy_zzzz = pbuffer.data(idx_dip_ig + 584);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxyyyy_xxxx, \
                             tr_y_xxyyyy_xxxy, \
                             tr_y_xxyyyy_xxxz, \
                             tr_y_xxyyyy_xxyy, \
                             tr_y_xxyyyy_xxyz, \
                             tr_y_xxyyyy_xxzz, \
                             tr_y_xxyyyy_xyyy, \
                             tr_y_xxyyyy_xyyz, \
                             tr_y_xxyyyy_xyzz, \
                             tr_y_xxyyyy_xzzz, \
                             tr_y_xxyyyy_yyyy, \
                             tr_y_xxyyyy_yyyz, \
                             tr_y_xxyyyy_yyzz, \
                             tr_y_xxyyyy_yzzz, \
                             tr_y_xxyyyy_zzzz, \
                             tr_y_xyyyy_xxx,   \
                             tr_y_xyyyy_xxxx,  \
                             tr_y_xyyyy_xxxy,  \
                             tr_y_xyyyy_xxxz,  \
                             tr_y_xyyyy_xxy,   \
                             tr_y_xyyyy_xxyy,  \
                             tr_y_xyyyy_xxyz,  \
                             tr_y_xyyyy_xxz,   \
                             tr_y_xyyyy_xxzz,  \
                             tr_y_xyyyy_xyy,   \
                             tr_y_xyyyy_xyyy,  \
                             tr_y_xyyyy_xyyz,  \
                             tr_y_xyyyy_xyz,   \
                             tr_y_xyyyy_xyzz,  \
                             tr_y_xyyyy_xzz,   \
                             tr_y_xyyyy_xzzz,  \
                             tr_y_xyyyy_yyy,   \
                             tr_y_xyyyy_yyyy,  \
                             tr_y_xyyyy_yyyz,  \
                             tr_y_xyyyy_yyz,   \
                             tr_y_xyyyy_yyzz,  \
                             tr_y_xyyyy_yzz,   \
                             tr_y_xyyyy_yzzz,  \
                             tr_y_xyyyy_zzz,   \
                             tr_y_xyyyy_zzzz,  \
                             tr_y_yyyy_xxxx,   \
                             tr_y_yyyy_xxxy,   \
                             tr_y_yyyy_xxxz,   \
                             tr_y_yyyy_xxyy,   \
                             tr_y_yyyy_xxyz,   \
                             tr_y_yyyy_xxzz,   \
                             tr_y_yyyy_xyyy,   \
                             tr_y_yyyy_xyyz,   \
                             tr_y_yyyy_xyzz,   \
                             tr_y_yyyy_xzzz,   \
                             tr_y_yyyy_yyyy,   \
                             tr_y_yyyy_yyyz,   \
                             tr_y_yyyy_yyzz,   \
                             tr_y_yyyy_yzzz,   \
                             tr_y_yyyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyy_xxxx[i] = tr_y_yyyy_xxxx[i] * fe_0 + 4.0 * tr_y_xyyyy_xxx[i] * fe_0 + tr_y_xyyyy_xxxx[i] * pa_x[i];

        tr_y_xxyyyy_xxxy[i] = tr_y_yyyy_xxxy[i] * fe_0 + 3.0 * tr_y_xyyyy_xxy[i] * fe_0 + tr_y_xyyyy_xxxy[i] * pa_x[i];

        tr_y_xxyyyy_xxxz[i] = tr_y_yyyy_xxxz[i] * fe_0 + 3.0 * tr_y_xyyyy_xxz[i] * fe_0 + tr_y_xyyyy_xxxz[i] * pa_x[i];

        tr_y_xxyyyy_xxyy[i] = tr_y_yyyy_xxyy[i] * fe_0 + 2.0 * tr_y_xyyyy_xyy[i] * fe_0 + tr_y_xyyyy_xxyy[i] * pa_x[i];

        tr_y_xxyyyy_xxyz[i] = tr_y_yyyy_xxyz[i] * fe_0 + 2.0 * tr_y_xyyyy_xyz[i] * fe_0 + tr_y_xyyyy_xxyz[i] * pa_x[i];

        tr_y_xxyyyy_xxzz[i] = tr_y_yyyy_xxzz[i] * fe_0 + 2.0 * tr_y_xyyyy_xzz[i] * fe_0 + tr_y_xyyyy_xxzz[i] * pa_x[i];

        tr_y_xxyyyy_xyyy[i] = tr_y_yyyy_xyyy[i] * fe_0 + tr_y_xyyyy_yyy[i] * fe_0 + tr_y_xyyyy_xyyy[i] * pa_x[i];

        tr_y_xxyyyy_xyyz[i] = tr_y_yyyy_xyyz[i] * fe_0 + tr_y_xyyyy_yyz[i] * fe_0 + tr_y_xyyyy_xyyz[i] * pa_x[i];

        tr_y_xxyyyy_xyzz[i] = tr_y_yyyy_xyzz[i] * fe_0 + tr_y_xyyyy_yzz[i] * fe_0 + tr_y_xyyyy_xyzz[i] * pa_x[i];

        tr_y_xxyyyy_xzzz[i] = tr_y_yyyy_xzzz[i] * fe_0 + tr_y_xyyyy_zzz[i] * fe_0 + tr_y_xyyyy_xzzz[i] * pa_x[i];

        tr_y_xxyyyy_yyyy[i] = tr_y_yyyy_yyyy[i] * fe_0 + tr_y_xyyyy_yyyy[i] * pa_x[i];

        tr_y_xxyyyy_yyyz[i] = tr_y_yyyy_yyyz[i] * fe_0 + tr_y_xyyyy_yyyz[i] * pa_x[i];

        tr_y_xxyyyy_yyzz[i] = tr_y_yyyy_yyzz[i] * fe_0 + tr_y_xyyyy_yyzz[i] * pa_x[i];

        tr_y_xxyyyy_yzzz[i] = tr_y_yyyy_yzzz[i] * fe_0 + tr_y_xyyyy_yzzz[i] * pa_x[i];

        tr_y_xxyyyy_zzzz[i] = tr_y_yyyy_zzzz[i] * fe_0 + tr_y_xyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 585-600 components of targeted buffer : IG

    auto tr_y_xxyyyz_xxxx = pbuffer.data(idx_dip_ig + 585);

    auto tr_y_xxyyyz_xxxy = pbuffer.data(idx_dip_ig + 586);

    auto tr_y_xxyyyz_xxxz = pbuffer.data(idx_dip_ig + 587);

    auto tr_y_xxyyyz_xxyy = pbuffer.data(idx_dip_ig + 588);

    auto tr_y_xxyyyz_xxyz = pbuffer.data(idx_dip_ig + 589);

    auto tr_y_xxyyyz_xxzz = pbuffer.data(idx_dip_ig + 590);

    auto tr_y_xxyyyz_xyyy = pbuffer.data(idx_dip_ig + 591);

    auto tr_y_xxyyyz_xyyz = pbuffer.data(idx_dip_ig + 592);

    auto tr_y_xxyyyz_xyzz = pbuffer.data(idx_dip_ig + 593);

    auto tr_y_xxyyyz_xzzz = pbuffer.data(idx_dip_ig + 594);

    auto tr_y_xxyyyz_yyyy = pbuffer.data(idx_dip_ig + 595);

    auto tr_y_xxyyyz_yyyz = pbuffer.data(idx_dip_ig + 596);

    auto tr_y_xxyyyz_yyzz = pbuffer.data(idx_dip_ig + 597);

    auto tr_y_xxyyyz_yzzz = pbuffer.data(idx_dip_ig + 598);

    auto tr_y_xxyyyz_zzzz = pbuffer.data(idx_dip_ig + 599);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxyyy_xxx,   \
                             tr_y_xxyyy_xxxx,  \
                             tr_y_xxyyy_xxxy,  \
                             tr_y_xxyyy_xxxz,  \
                             tr_y_xxyyy_xxy,   \
                             tr_y_xxyyy_xxyy,  \
                             tr_y_xxyyy_xxyz,  \
                             tr_y_xxyyy_xxz,   \
                             tr_y_xxyyy_xxzz,  \
                             tr_y_xxyyy_xyy,   \
                             tr_y_xxyyy_xyyy,  \
                             tr_y_xxyyy_xyyz,  \
                             tr_y_xxyyy_xyz,   \
                             tr_y_xxyyy_xyzz,  \
                             tr_y_xxyyy_xzz,   \
                             tr_y_xxyyy_xzzz,  \
                             tr_y_xxyyy_yyyy,  \
                             tr_y_xxyyyz_xxxx, \
                             tr_y_xxyyyz_xxxy, \
                             tr_y_xxyyyz_xxxz, \
                             tr_y_xxyyyz_xxyy, \
                             tr_y_xxyyyz_xxyz, \
                             tr_y_xxyyyz_xxzz, \
                             tr_y_xxyyyz_xyyy, \
                             tr_y_xxyyyz_xyyz, \
                             tr_y_xxyyyz_xyzz, \
                             tr_y_xxyyyz_xzzz, \
                             tr_y_xxyyyz_yyyy, \
                             tr_y_xxyyyz_yyyz, \
                             tr_y_xxyyyz_yyzz, \
                             tr_y_xxyyyz_yzzz, \
                             tr_y_xxyyyz_zzzz, \
                             tr_y_xyyyz_yyyz,  \
                             tr_y_xyyyz_yyzz,  \
                             tr_y_xyyyz_yzzz,  \
                             tr_y_xyyyz_zzzz,  \
                             tr_y_yyyz_yyyz,   \
                             tr_y_yyyz_yyzz,   \
                             tr_y_yyyz_yzzz,   \
                             tr_y_yyyz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyz_xxxx[i] = tr_y_xxyyy_xxxx[i] * pa_z[i];

        tr_y_xxyyyz_xxxy[i] = tr_y_xxyyy_xxxy[i] * pa_z[i];

        tr_y_xxyyyz_xxxz[i] = tr_y_xxyyy_xxx[i] * fe_0 + tr_y_xxyyy_xxxz[i] * pa_z[i];

        tr_y_xxyyyz_xxyy[i] = tr_y_xxyyy_xxyy[i] * pa_z[i];

        tr_y_xxyyyz_xxyz[i] = tr_y_xxyyy_xxy[i] * fe_0 + tr_y_xxyyy_xxyz[i] * pa_z[i];

        tr_y_xxyyyz_xxzz[i] = 2.0 * tr_y_xxyyy_xxz[i] * fe_0 + tr_y_xxyyy_xxzz[i] * pa_z[i];

        tr_y_xxyyyz_xyyy[i] = tr_y_xxyyy_xyyy[i] * pa_z[i];

        tr_y_xxyyyz_xyyz[i] = tr_y_xxyyy_xyy[i] * fe_0 + tr_y_xxyyy_xyyz[i] * pa_z[i];

        tr_y_xxyyyz_xyzz[i] = 2.0 * tr_y_xxyyy_xyz[i] * fe_0 + tr_y_xxyyy_xyzz[i] * pa_z[i];

        tr_y_xxyyyz_xzzz[i] = 3.0 * tr_y_xxyyy_xzz[i] * fe_0 + tr_y_xxyyy_xzzz[i] * pa_z[i];

        tr_y_xxyyyz_yyyy[i] = tr_y_xxyyy_yyyy[i] * pa_z[i];

        tr_y_xxyyyz_yyyz[i] = tr_y_yyyz_yyyz[i] * fe_0 + tr_y_xyyyz_yyyz[i] * pa_x[i];

        tr_y_xxyyyz_yyzz[i] = tr_y_yyyz_yyzz[i] * fe_0 + tr_y_xyyyz_yyzz[i] * pa_x[i];

        tr_y_xxyyyz_yzzz[i] = tr_y_yyyz_yzzz[i] * fe_0 + tr_y_xyyyz_yzzz[i] * pa_x[i];

        tr_y_xxyyyz_zzzz[i] = tr_y_yyyz_zzzz[i] * fe_0 + tr_y_xyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 600-615 components of targeted buffer : IG

    auto tr_y_xxyyzz_xxxx = pbuffer.data(idx_dip_ig + 600);

    auto tr_y_xxyyzz_xxxy = pbuffer.data(idx_dip_ig + 601);

    auto tr_y_xxyyzz_xxxz = pbuffer.data(idx_dip_ig + 602);

    auto tr_y_xxyyzz_xxyy = pbuffer.data(idx_dip_ig + 603);

    auto tr_y_xxyyzz_xxyz = pbuffer.data(idx_dip_ig + 604);

    auto tr_y_xxyyzz_xxzz = pbuffer.data(idx_dip_ig + 605);

    auto tr_y_xxyyzz_xyyy = pbuffer.data(idx_dip_ig + 606);

    auto tr_y_xxyyzz_xyyz = pbuffer.data(idx_dip_ig + 607);

    auto tr_y_xxyyzz_xyzz = pbuffer.data(idx_dip_ig + 608);

    auto tr_y_xxyyzz_xzzz = pbuffer.data(idx_dip_ig + 609);

    auto tr_y_xxyyzz_yyyy = pbuffer.data(idx_dip_ig + 610);

    auto tr_y_xxyyzz_yyyz = pbuffer.data(idx_dip_ig + 611);

    auto tr_y_xxyyzz_yyzz = pbuffer.data(idx_dip_ig + 612);

    auto tr_y_xxyyzz_yzzz = pbuffer.data(idx_dip_ig + 613);

    auto tr_y_xxyyzz_zzzz = pbuffer.data(idx_dip_ig + 614);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxyy_xxxx,   \
                             tr_y_xxyy_xxxy,   \
                             tr_y_xxyy_xxyy,   \
                             tr_y_xxyy_xyyy,   \
                             tr_y_xxyyz_xxxx,  \
                             tr_y_xxyyz_xxxy,  \
                             tr_y_xxyyz_xxyy,  \
                             tr_y_xxyyz_xyyy,  \
                             tr_y_xxyyzz_xxxx, \
                             tr_y_xxyyzz_xxxy, \
                             tr_y_xxyyzz_xxxz, \
                             tr_y_xxyyzz_xxyy, \
                             tr_y_xxyyzz_xxyz, \
                             tr_y_xxyyzz_xxzz, \
                             tr_y_xxyyzz_xyyy, \
                             tr_y_xxyyzz_xyyz, \
                             tr_y_xxyyzz_xyzz, \
                             tr_y_xxyyzz_xzzz, \
                             tr_y_xxyyzz_yyyy, \
                             tr_y_xxyyzz_yyyz, \
                             tr_y_xxyyzz_yyzz, \
                             tr_y_xxyyzz_yzzz, \
                             tr_y_xxyyzz_zzzz, \
                             tr_y_xyyzz_xxxz,  \
                             tr_y_xyyzz_xxyz,  \
                             tr_y_xyyzz_xxz,   \
                             tr_y_xyyzz_xxzz,  \
                             tr_y_xyyzz_xyyz,  \
                             tr_y_xyyzz_xyz,   \
                             tr_y_xyyzz_xyzz,  \
                             tr_y_xyyzz_xzz,   \
                             tr_y_xyyzz_xzzz,  \
                             tr_y_xyyzz_yyyy,  \
                             tr_y_xyyzz_yyyz,  \
                             tr_y_xyyzz_yyz,   \
                             tr_y_xyyzz_yyzz,  \
                             tr_y_xyyzz_yzz,   \
                             tr_y_xyyzz_yzzz,  \
                             tr_y_xyyzz_zzz,   \
                             tr_y_xyyzz_zzzz,  \
                             tr_y_yyzz_xxxz,   \
                             tr_y_yyzz_xxyz,   \
                             tr_y_yyzz_xxzz,   \
                             tr_y_yyzz_xyyz,   \
                             tr_y_yyzz_xyzz,   \
                             tr_y_yyzz_xzzz,   \
                             tr_y_yyzz_yyyy,   \
                             tr_y_yyzz_yyyz,   \
                             tr_y_yyzz_yyzz,   \
                             tr_y_yyzz_yzzz,   \
                             tr_y_yyzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyzz_xxxx[i] = tr_y_xxyy_xxxx[i] * fe_0 + tr_y_xxyyz_xxxx[i] * pa_z[i];

        tr_y_xxyyzz_xxxy[i] = tr_y_xxyy_xxxy[i] * fe_0 + tr_y_xxyyz_xxxy[i] * pa_z[i];

        tr_y_xxyyzz_xxxz[i] = tr_y_yyzz_xxxz[i] * fe_0 + 3.0 * tr_y_xyyzz_xxz[i] * fe_0 + tr_y_xyyzz_xxxz[i] * pa_x[i];

        tr_y_xxyyzz_xxyy[i] = tr_y_xxyy_xxyy[i] * fe_0 + tr_y_xxyyz_xxyy[i] * pa_z[i];

        tr_y_xxyyzz_xxyz[i] = tr_y_yyzz_xxyz[i] * fe_0 + 2.0 * tr_y_xyyzz_xyz[i] * fe_0 + tr_y_xyyzz_xxyz[i] * pa_x[i];

        tr_y_xxyyzz_xxzz[i] = tr_y_yyzz_xxzz[i] * fe_0 + 2.0 * tr_y_xyyzz_xzz[i] * fe_0 + tr_y_xyyzz_xxzz[i] * pa_x[i];

        tr_y_xxyyzz_xyyy[i] = tr_y_xxyy_xyyy[i] * fe_0 + tr_y_xxyyz_xyyy[i] * pa_z[i];

        tr_y_xxyyzz_xyyz[i] = tr_y_yyzz_xyyz[i] * fe_0 + tr_y_xyyzz_yyz[i] * fe_0 + tr_y_xyyzz_xyyz[i] * pa_x[i];

        tr_y_xxyyzz_xyzz[i] = tr_y_yyzz_xyzz[i] * fe_0 + tr_y_xyyzz_yzz[i] * fe_0 + tr_y_xyyzz_xyzz[i] * pa_x[i];

        tr_y_xxyyzz_xzzz[i] = tr_y_yyzz_xzzz[i] * fe_0 + tr_y_xyyzz_zzz[i] * fe_0 + tr_y_xyyzz_xzzz[i] * pa_x[i];

        tr_y_xxyyzz_yyyy[i] = tr_y_yyzz_yyyy[i] * fe_0 + tr_y_xyyzz_yyyy[i] * pa_x[i];

        tr_y_xxyyzz_yyyz[i] = tr_y_yyzz_yyyz[i] * fe_0 + tr_y_xyyzz_yyyz[i] * pa_x[i];

        tr_y_xxyyzz_yyzz[i] = tr_y_yyzz_yyzz[i] * fe_0 + tr_y_xyyzz_yyzz[i] * pa_x[i];

        tr_y_xxyyzz_yzzz[i] = tr_y_yyzz_yzzz[i] * fe_0 + tr_y_xyyzz_yzzz[i] * pa_x[i];

        tr_y_xxyyzz_zzzz[i] = tr_y_yyzz_zzzz[i] * fe_0 + tr_y_xyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 615-630 components of targeted buffer : IG

    auto tr_y_xxyzzz_xxxx = pbuffer.data(idx_dip_ig + 615);

    auto tr_y_xxyzzz_xxxy = pbuffer.data(idx_dip_ig + 616);

    auto tr_y_xxyzzz_xxxz = pbuffer.data(idx_dip_ig + 617);

    auto tr_y_xxyzzz_xxyy = pbuffer.data(idx_dip_ig + 618);

    auto tr_y_xxyzzz_xxyz = pbuffer.data(idx_dip_ig + 619);

    auto tr_y_xxyzzz_xxzz = pbuffer.data(idx_dip_ig + 620);

    auto tr_y_xxyzzz_xyyy = pbuffer.data(idx_dip_ig + 621);

    auto tr_y_xxyzzz_xyyz = pbuffer.data(idx_dip_ig + 622);

    auto tr_y_xxyzzz_xyzz = pbuffer.data(idx_dip_ig + 623);

    auto tr_y_xxyzzz_xzzz = pbuffer.data(idx_dip_ig + 624);

    auto tr_y_xxyzzz_yyyy = pbuffer.data(idx_dip_ig + 625);

    auto tr_y_xxyzzz_yyyz = pbuffer.data(idx_dip_ig + 626);

    auto tr_y_xxyzzz_yyzz = pbuffer.data(idx_dip_ig + 627);

    auto tr_y_xxyzzz_yzzz = pbuffer.data(idx_dip_ig + 628);

    auto tr_y_xxyzzz_zzzz = pbuffer.data(idx_dip_ig + 629);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_y_xxyz_xxxy,   \
                             tr_y_xxyz_xxyy,   \
                             tr_y_xxyz_xyyy,   \
                             tr_y_xxyzz_xxxy,  \
                             tr_y_xxyzz_xxyy,  \
                             tr_y_xxyzz_xyyy,  \
                             tr_y_xxyzzz_xxxx, \
                             tr_y_xxyzzz_xxxy, \
                             tr_y_xxyzzz_xxxz, \
                             tr_y_xxyzzz_xxyy, \
                             tr_y_xxyzzz_xxyz, \
                             tr_y_xxyzzz_xxzz, \
                             tr_y_xxyzzz_xyyy, \
                             tr_y_xxyzzz_xyyz, \
                             tr_y_xxyzzz_xyzz, \
                             tr_y_xxyzzz_xzzz, \
                             tr_y_xxyzzz_yyyy, \
                             tr_y_xxyzzz_yyyz, \
                             tr_y_xxyzzz_yyzz, \
                             tr_y_xxyzzz_yzzz, \
                             tr_y_xxyzzz_zzzz, \
                             tr_y_xxzzz_xxxx,  \
                             tr_y_xxzzz_xxxz,  \
                             tr_y_xxzzz_xxzz,  \
                             tr_y_xxzzz_xzzz,  \
                             tr_y_xyzzz_xxyz,  \
                             tr_y_xyzzz_xyyz,  \
                             tr_y_xyzzz_xyz,   \
                             tr_y_xyzzz_xyzz,  \
                             tr_y_xyzzz_yyyy,  \
                             tr_y_xyzzz_yyyz,  \
                             tr_y_xyzzz_yyz,   \
                             tr_y_xyzzz_yyzz,  \
                             tr_y_xyzzz_yzz,   \
                             tr_y_xyzzz_yzzz,  \
                             tr_y_xyzzz_zzzz,  \
                             tr_y_yzzz_xxyz,   \
                             tr_y_yzzz_xyyz,   \
                             tr_y_yzzz_xyzz,   \
                             tr_y_yzzz_yyyy,   \
                             tr_y_yzzz_yyyz,   \
                             tr_y_yzzz_yyzz,   \
                             tr_y_yzzz_yzzz,   \
                             tr_y_yzzz_zzzz,   \
                             ts_xxzzz_xxxx,    \
                             ts_xxzzz_xxxz,    \
                             ts_xxzzz_xxzz,    \
                             ts_xxzzz_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzzz_xxxx[i] = ts_xxzzz_xxxx[i] * fe_0 + tr_y_xxzzz_xxxx[i] * pa_y[i];

        tr_y_xxyzzz_xxxy[i] = 2.0 * tr_y_xxyz_xxxy[i] * fe_0 + tr_y_xxyzz_xxxy[i] * pa_z[i];

        tr_y_xxyzzz_xxxz[i] = ts_xxzzz_xxxz[i] * fe_0 + tr_y_xxzzz_xxxz[i] * pa_y[i];

        tr_y_xxyzzz_xxyy[i] = 2.0 * tr_y_xxyz_xxyy[i] * fe_0 + tr_y_xxyzz_xxyy[i] * pa_z[i];

        tr_y_xxyzzz_xxyz[i] = tr_y_yzzz_xxyz[i] * fe_0 + 2.0 * tr_y_xyzzz_xyz[i] * fe_0 + tr_y_xyzzz_xxyz[i] * pa_x[i];

        tr_y_xxyzzz_xxzz[i] = ts_xxzzz_xxzz[i] * fe_0 + tr_y_xxzzz_xxzz[i] * pa_y[i];

        tr_y_xxyzzz_xyyy[i] = 2.0 * tr_y_xxyz_xyyy[i] * fe_0 + tr_y_xxyzz_xyyy[i] * pa_z[i];

        tr_y_xxyzzz_xyyz[i] = tr_y_yzzz_xyyz[i] * fe_0 + tr_y_xyzzz_yyz[i] * fe_0 + tr_y_xyzzz_xyyz[i] * pa_x[i];

        tr_y_xxyzzz_xyzz[i] = tr_y_yzzz_xyzz[i] * fe_0 + tr_y_xyzzz_yzz[i] * fe_0 + tr_y_xyzzz_xyzz[i] * pa_x[i];

        tr_y_xxyzzz_xzzz[i] = ts_xxzzz_xzzz[i] * fe_0 + tr_y_xxzzz_xzzz[i] * pa_y[i];

        tr_y_xxyzzz_yyyy[i] = tr_y_yzzz_yyyy[i] * fe_0 + tr_y_xyzzz_yyyy[i] * pa_x[i];

        tr_y_xxyzzz_yyyz[i] = tr_y_yzzz_yyyz[i] * fe_0 + tr_y_xyzzz_yyyz[i] * pa_x[i];

        tr_y_xxyzzz_yyzz[i] = tr_y_yzzz_yyzz[i] * fe_0 + tr_y_xyzzz_yyzz[i] * pa_x[i];

        tr_y_xxyzzz_yzzz[i] = tr_y_yzzz_yzzz[i] * fe_0 + tr_y_xyzzz_yzzz[i] * pa_x[i];

        tr_y_xxyzzz_zzzz[i] = tr_y_yzzz_zzzz[i] * fe_0 + tr_y_xyzzz_zzzz[i] * pa_x[i];
    }

    // Set up 630-645 components of targeted buffer : IG

    auto tr_y_xxzzzz_xxxx = pbuffer.data(idx_dip_ig + 630);

    auto tr_y_xxzzzz_xxxy = pbuffer.data(idx_dip_ig + 631);

    auto tr_y_xxzzzz_xxxz = pbuffer.data(idx_dip_ig + 632);

    auto tr_y_xxzzzz_xxyy = pbuffer.data(idx_dip_ig + 633);

    auto tr_y_xxzzzz_xxyz = pbuffer.data(idx_dip_ig + 634);

    auto tr_y_xxzzzz_xxzz = pbuffer.data(idx_dip_ig + 635);

    auto tr_y_xxzzzz_xyyy = pbuffer.data(idx_dip_ig + 636);

    auto tr_y_xxzzzz_xyyz = pbuffer.data(idx_dip_ig + 637);

    auto tr_y_xxzzzz_xyzz = pbuffer.data(idx_dip_ig + 638);

    auto tr_y_xxzzzz_xzzz = pbuffer.data(idx_dip_ig + 639);

    auto tr_y_xxzzzz_yyyy = pbuffer.data(idx_dip_ig + 640);

    auto tr_y_xxzzzz_yyyz = pbuffer.data(idx_dip_ig + 641);

    auto tr_y_xxzzzz_yyzz = pbuffer.data(idx_dip_ig + 642);

    auto tr_y_xxzzzz_yzzz = pbuffer.data(idx_dip_ig + 643);

    auto tr_y_xxzzzz_zzzz = pbuffer.data(idx_dip_ig + 644);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxzz_xxxx,   \
                             tr_y_xxzz_xxxy,   \
                             tr_y_xxzz_xxyy,   \
                             tr_y_xxzz_xyyy,   \
                             tr_y_xxzzz_xxxx,  \
                             tr_y_xxzzz_xxxy,  \
                             tr_y_xxzzz_xxyy,  \
                             tr_y_xxzzz_xyyy,  \
                             tr_y_xxzzzz_xxxx, \
                             tr_y_xxzzzz_xxxy, \
                             tr_y_xxzzzz_xxxz, \
                             tr_y_xxzzzz_xxyy, \
                             tr_y_xxzzzz_xxyz, \
                             tr_y_xxzzzz_xxzz, \
                             tr_y_xxzzzz_xyyy, \
                             tr_y_xxzzzz_xyyz, \
                             tr_y_xxzzzz_xyzz, \
                             tr_y_xxzzzz_xzzz, \
                             tr_y_xxzzzz_yyyy, \
                             tr_y_xxzzzz_yyyz, \
                             tr_y_xxzzzz_yyzz, \
                             tr_y_xxzzzz_yzzz, \
                             tr_y_xxzzzz_zzzz, \
                             tr_y_xzzzz_xxxz,  \
                             tr_y_xzzzz_xxyz,  \
                             tr_y_xzzzz_xxz,   \
                             tr_y_xzzzz_xxzz,  \
                             tr_y_xzzzz_xyyz,  \
                             tr_y_xzzzz_xyz,   \
                             tr_y_xzzzz_xyzz,  \
                             tr_y_xzzzz_xzz,   \
                             tr_y_xzzzz_xzzz,  \
                             tr_y_xzzzz_yyyy,  \
                             tr_y_xzzzz_yyyz,  \
                             tr_y_xzzzz_yyz,   \
                             tr_y_xzzzz_yyzz,  \
                             tr_y_xzzzz_yzz,   \
                             tr_y_xzzzz_yzzz,  \
                             tr_y_xzzzz_zzz,   \
                             tr_y_xzzzz_zzzz,  \
                             tr_y_zzzz_xxxz,   \
                             tr_y_zzzz_xxyz,   \
                             tr_y_zzzz_xxzz,   \
                             tr_y_zzzz_xyyz,   \
                             tr_y_zzzz_xyzz,   \
                             tr_y_zzzz_xzzz,   \
                             tr_y_zzzz_yyyy,   \
                             tr_y_zzzz_yyyz,   \
                             tr_y_zzzz_yyzz,   \
                             tr_y_zzzz_yzzz,   \
                             tr_y_zzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzzz_xxxx[i] = 3.0 * tr_y_xxzz_xxxx[i] * fe_0 + tr_y_xxzzz_xxxx[i] * pa_z[i];

        tr_y_xxzzzz_xxxy[i] = 3.0 * tr_y_xxzz_xxxy[i] * fe_0 + tr_y_xxzzz_xxxy[i] * pa_z[i];

        tr_y_xxzzzz_xxxz[i] = tr_y_zzzz_xxxz[i] * fe_0 + 3.0 * tr_y_xzzzz_xxz[i] * fe_0 + tr_y_xzzzz_xxxz[i] * pa_x[i];

        tr_y_xxzzzz_xxyy[i] = 3.0 * tr_y_xxzz_xxyy[i] * fe_0 + tr_y_xxzzz_xxyy[i] * pa_z[i];

        tr_y_xxzzzz_xxyz[i] = tr_y_zzzz_xxyz[i] * fe_0 + 2.0 * tr_y_xzzzz_xyz[i] * fe_0 + tr_y_xzzzz_xxyz[i] * pa_x[i];

        tr_y_xxzzzz_xxzz[i] = tr_y_zzzz_xxzz[i] * fe_0 + 2.0 * tr_y_xzzzz_xzz[i] * fe_0 + tr_y_xzzzz_xxzz[i] * pa_x[i];

        tr_y_xxzzzz_xyyy[i] = 3.0 * tr_y_xxzz_xyyy[i] * fe_0 + tr_y_xxzzz_xyyy[i] * pa_z[i];

        tr_y_xxzzzz_xyyz[i] = tr_y_zzzz_xyyz[i] * fe_0 + tr_y_xzzzz_yyz[i] * fe_0 + tr_y_xzzzz_xyyz[i] * pa_x[i];

        tr_y_xxzzzz_xyzz[i] = tr_y_zzzz_xyzz[i] * fe_0 + tr_y_xzzzz_yzz[i] * fe_0 + tr_y_xzzzz_xyzz[i] * pa_x[i];

        tr_y_xxzzzz_xzzz[i] = tr_y_zzzz_xzzz[i] * fe_0 + tr_y_xzzzz_zzz[i] * fe_0 + tr_y_xzzzz_xzzz[i] * pa_x[i];

        tr_y_xxzzzz_yyyy[i] = tr_y_zzzz_yyyy[i] * fe_0 + tr_y_xzzzz_yyyy[i] * pa_x[i];

        tr_y_xxzzzz_yyyz[i] = tr_y_zzzz_yyyz[i] * fe_0 + tr_y_xzzzz_yyyz[i] * pa_x[i];

        tr_y_xxzzzz_yyzz[i] = tr_y_zzzz_yyzz[i] * fe_0 + tr_y_xzzzz_yyzz[i] * pa_x[i];

        tr_y_xxzzzz_yzzz[i] = tr_y_zzzz_yzzz[i] * fe_0 + tr_y_xzzzz_yzzz[i] * pa_x[i];

        tr_y_xxzzzz_zzzz[i] = tr_y_zzzz_zzzz[i] * fe_0 + tr_y_xzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 645-660 components of targeted buffer : IG

    auto tr_y_xyyyyy_xxxx = pbuffer.data(idx_dip_ig + 645);

    auto tr_y_xyyyyy_xxxy = pbuffer.data(idx_dip_ig + 646);

    auto tr_y_xyyyyy_xxxz = pbuffer.data(idx_dip_ig + 647);

    auto tr_y_xyyyyy_xxyy = pbuffer.data(idx_dip_ig + 648);

    auto tr_y_xyyyyy_xxyz = pbuffer.data(idx_dip_ig + 649);

    auto tr_y_xyyyyy_xxzz = pbuffer.data(idx_dip_ig + 650);

    auto tr_y_xyyyyy_xyyy = pbuffer.data(idx_dip_ig + 651);

    auto tr_y_xyyyyy_xyyz = pbuffer.data(idx_dip_ig + 652);

    auto tr_y_xyyyyy_xyzz = pbuffer.data(idx_dip_ig + 653);

    auto tr_y_xyyyyy_xzzz = pbuffer.data(idx_dip_ig + 654);

    auto tr_y_xyyyyy_yyyy = pbuffer.data(idx_dip_ig + 655);

    auto tr_y_xyyyyy_yyyz = pbuffer.data(idx_dip_ig + 656);

    auto tr_y_xyyyyy_yyzz = pbuffer.data(idx_dip_ig + 657);

    auto tr_y_xyyyyy_yzzz = pbuffer.data(idx_dip_ig + 658);

    auto tr_y_xyyyyy_zzzz = pbuffer.data(idx_dip_ig + 659);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyyyyy_xxxx, \
                             tr_y_xyyyyy_xxxy, \
                             tr_y_xyyyyy_xxxz, \
                             tr_y_xyyyyy_xxyy, \
                             tr_y_xyyyyy_xxyz, \
                             tr_y_xyyyyy_xxzz, \
                             tr_y_xyyyyy_xyyy, \
                             tr_y_xyyyyy_xyyz, \
                             tr_y_xyyyyy_xyzz, \
                             tr_y_xyyyyy_xzzz, \
                             tr_y_xyyyyy_yyyy, \
                             tr_y_xyyyyy_yyyz, \
                             tr_y_xyyyyy_yyzz, \
                             tr_y_xyyyyy_yzzz, \
                             tr_y_xyyyyy_zzzz, \
                             tr_y_yyyyy_xxx,   \
                             tr_y_yyyyy_xxxx,  \
                             tr_y_yyyyy_xxxy,  \
                             tr_y_yyyyy_xxxz,  \
                             tr_y_yyyyy_xxy,   \
                             tr_y_yyyyy_xxyy,  \
                             tr_y_yyyyy_xxyz,  \
                             tr_y_yyyyy_xxz,   \
                             tr_y_yyyyy_xxzz,  \
                             tr_y_yyyyy_xyy,   \
                             tr_y_yyyyy_xyyy,  \
                             tr_y_yyyyy_xyyz,  \
                             tr_y_yyyyy_xyz,   \
                             tr_y_yyyyy_xyzz,  \
                             tr_y_yyyyy_xzz,   \
                             tr_y_yyyyy_xzzz,  \
                             tr_y_yyyyy_yyy,   \
                             tr_y_yyyyy_yyyy,  \
                             tr_y_yyyyy_yyyz,  \
                             tr_y_yyyyy_yyz,   \
                             tr_y_yyyyy_yyzz,  \
                             tr_y_yyyyy_yzz,   \
                             tr_y_yyyyy_yzzz,  \
                             tr_y_yyyyy_zzz,   \
                             tr_y_yyyyy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyy_xxxx[i] = 4.0 * tr_y_yyyyy_xxx[i] * fe_0 + tr_y_yyyyy_xxxx[i] * pa_x[i];

        tr_y_xyyyyy_xxxy[i] = 3.0 * tr_y_yyyyy_xxy[i] * fe_0 + tr_y_yyyyy_xxxy[i] * pa_x[i];

        tr_y_xyyyyy_xxxz[i] = 3.0 * tr_y_yyyyy_xxz[i] * fe_0 + tr_y_yyyyy_xxxz[i] * pa_x[i];

        tr_y_xyyyyy_xxyy[i] = 2.0 * tr_y_yyyyy_xyy[i] * fe_0 + tr_y_yyyyy_xxyy[i] * pa_x[i];

        tr_y_xyyyyy_xxyz[i] = 2.0 * tr_y_yyyyy_xyz[i] * fe_0 + tr_y_yyyyy_xxyz[i] * pa_x[i];

        tr_y_xyyyyy_xxzz[i] = 2.0 * tr_y_yyyyy_xzz[i] * fe_0 + tr_y_yyyyy_xxzz[i] * pa_x[i];

        tr_y_xyyyyy_xyyy[i] = tr_y_yyyyy_yyy[i] * fe_0 + tr_y_yyyyy_xyyy[i] * pa_x[i];

        tr_y_xyyyyy_xyyz[i] = tr_y_yyyyy_yyz[i] * fe_0 + tr_y_yyyyy_xyyz[i] * pa_x[i];

        tr_y_xyyyyy_xyzz[i] = tr_y_yyyyy_yzz[i] * fe_0 + tr_y_yyyyy_xyzz[i] * pa_x[i];

        tr_y_xyyyyy_xzzz[i] = tr_y_yyyyy_zzz[i] * fe_0 + tr_y_yyyyy_xzzz[i] * pa_x[i];

        tr_y_xyyyyy_yyyy[i] = tr_y_yyyyy_yyyy[i] * pa_x[i];

        tr_y_xyyyyy_yyyz[i] = tr_y_yyyyy_yyyz[i] * pa_x[i];

        tr_y_xyyyyy_yyzz[i] = tr_y_yyyyy_yyzz[i] * pa_x[i];

        tr_y_xyyyyy_yzzz[i] = tr_y_yyyyy_yzzz[i] * pa_x[i];

        tr_y_xyyyyy_zzzz[i] = tr_y_yyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 660-675 components of targeted buffer : IG

    auto tr_y_xyyyyz_xxxx = pbuffer.data(idx_dip_ig + 660);

    auto tr_y_xyyyyz_xxxy = pbuffer.data(idx_dip_ig + 661);

    auto tr_y_xyyyyz_xxxz = pbuffer.data(idx_dip_ig + 662);

    auto tr_y_xyyyyz_xxyy = pbuffer.data(idx_dip_ig + 663);

    auto tr_y_xyyyyz_xxyz = pbuffer.data(idx_dip_ig + 664);

    auto tr_y_xyyyyz_xxzz = pbuffer.data(idx_dip_ig + 665);

    auto tr_y_xyyyyz_xyyy = pbuffer.data(idx_dip_ig + 666);

    auto tr_y_xyyyyz_xyyz = pbuffer.data(idx_dip_ig + 667);

    auto tr_y_xyyyyz_xyzz = pbuffer.data(idx_dip_ig + 668);

    auto tr_y_xyyyyz_xzzz = pbuffer.data(idx_dip_ig + 669);

    auto tr_y_xyyyyz_yyyy = pbuffer.data(idx_dip_ig + 670);

    auto tr_y_xyyyyz_yyyz = pbuffer.data(idx_dip_ig + 671);

    auto tr_y_xyyyyz_yyzz = pbuffer.data(idx_dip_ig + 672);

    auto tr_y_xyyyyz_yzzz = pbuffer.data(idx_dip_ig + 673);

    auto tr_y_xyyyyz_zzzz = pbuffer.data(idx_dip_ig + 674);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xyyyy_xxxx,  \
                             tr_y_xyyyy_xxxy,  \
                             tr_y_xyyyy_xxyy,  \
                             tr_y_xyyyy_xyyy,  \
                             tr_y_xyyyyz_xxxx, \
                             tr_y_xyyyyz_xxxy, \
                             tr_y_xyyyyz_xxxz, \
                             tr_y_xyyyyz_xxyy, \
                             tr_y_xyyyyz_xxyz, \
                             tr_y_xyyyyz_xxzz, \
                             tr_y_xyyyyz_xyyy, \
                             tr_y_xyyyyz_xyyz, \
                             tr_y_xyyyyz_xyzz, \
                             tr_y_xyyyyz_xzzz, \
                             tr_y_xyyyyz_yyyy, \
                             tr_y_xyyyyz_yyyz, \
                             tr_y_xyyyyz_yyzz, \
                             tr_y_xyyyyz_yzzz, \
                             tr_y_xyyyyz_zzzz, \
                             tr_y_yyyyz_xxxz,  \
                             tr_y_yyyyz_xxyz,  \
                             tr_y_yyyyz_xxz,   \
                             tr_y_yyyyz_xxzz,  \
                             tr_y_yyyyz_xyyz,  \
                             tr_y_yyyyz_xyz,   \
                             tr_y_yyyyz_xyzz,  \
                             tr_y_yyyyz_xzz,   \
                             tr_y_yyyyz_xzzz,  \
                             tr_y_yyyyz_yyyy,  \
                             tr_y_yyyyz_yyyz,  \
                             tr_y_yyyyz_yyz,   \
                             tr_y_yyyyz_yyzz,  \
                             tr_y_yyyyz_yzz,   \
                             tr_y_yyyyz_yzzz,  \
                             tr_y_yyyyz_zzz,   \
                             tr_y_yyyyz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyz_xxxx[i] = tr_y_xyyyy_xxxx[i] * pa_z[i];

        tr_y_xyyyyz_xxxy[i] = tr_y_xyyyy_xxxy[i] * pa_z[i];

        tr_y_xyyyyz_xxxz[i] = 3.0 * tr_y_yyyyz_xxz[i] * fe_0 + tr_y_yyyyz_xxxz[i] * pa_x[i];

        tr_y_xyyyyz_xxyy[i] = tr_y_xyyyy_xxyy[i] * pa_z[i];

        tr_y_xyyyyz_xxyz[i] = 2.0 * tr_y_yyyyz_xyz[i] * fe_0 + tr_y_yyyyz_xxyz[i] * pa_x[i];

        tr_y_xyyyyz_xxzz[i] = 2.0 * tr_y_yyyyz_xzz[i] * fe_0 + tr_y_yyyyz_xxzz[i] * pa_x[i];

        tr_y_xyyyyz_xyyy[i] = tr_y_xyyyy_xyyy[i] * pa_z[i];

        tr_y_xyyyyz_xyyz[i] = tr_y_yyyyz_yyz[i] * fe_0 + tr_y_yyyyz_xyyz[i] * pa_x[i];

        tr_y_xyyyyz_xyzz[i] = tr_y_yyyyz_yzz[i] * fe_0 + tr_y_yyyyz_xyzz[i] * pa_x[i];

        tr_y_xyyyyz_xzzz[i] = tr_y_yyyyz_zzz[i] * fe_0 + tr_y_yyyyz_xzzz[i] * pa_x[i];

        tr_y_xyyyyz_yyyy[i] = tr_y_yyyyz_yyyy[i] * pa_x[i];

        tr_y_xyyyyz_yyyz[i] = tr_y_yyyyz_yyyz[i] * pa_x[i];

        tr_y_xyyyyz_yyzz[i] = tr_y_yyyyz_yyzz[i] * pa_x[i];

        tr_y_xyyyyz_yzzz[i] = tr_y_yyyyz_yzzz[i] * pa_x[i];

        tr_y_xyyyyz_zzzz[i] = tr_y_yyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 675-690 components of targeted buffer : IG

    auto tr_y_xyyyzz_xxxx = pbuffer.data(idx_dip_ig + 675);

    auto tr_y_xyyyzz_xxxy = pbuffer.data(idx_dip_ig + 676);

    auto tr_y_xyyyzz_xxxz = pbuffer.data(idx_dip_ig + 677);

    auto tr_y_xyyyzz_xxyy = pbuffer.data(idx_dip_ig + 678);

    auto tr_y_xyyyzz_xxyz = pbuffer.data(idx_dip_ig + 679);

    auto tr_y_xyyyzz_xxzz = pbuffer.data(idx_dip_ig + 680);

    auto tr_y_xyyyzz_xyyy = pbuffer.data(idx_dip_ig + 681);

    auto tr_y_xyyyzz_xyyz = pbuffer.data(idx_dip_ig + 682);

    auto tr_y_xyyyzz_xyzz = pbuffer.data(idx_dip_ig + 683);

    auto tr_y_xyyyzz_xzzz = pbuffer.data(idx_dip_ig + 684);

    auto tr_y_xyyyzz_yyyy = pbuffer.data(idx_dip_ig + 685);

    auto tr_y_xyyyzz_yyyz = pbuffer.data(idx_dip_ig + 686);

    auto tr_y_xyyyzz_yyzz = pbuffer.data(idx_dip_ig + 687);

    auto tr_y_xyyyzz_yzzz = pbuffer.data(idx_dip_ig + 688);

    auto tr_y_xyyyzz_zzzz = pbuffer.data(idx_dip_ig + 689);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyyyzz_xxxx, \
                             tr_y_xyyyzz_xxxy, \
                             tr_y_xyyyzz_xxxz, \
                             tr_y_xyyyzz_xxyy, \
                             tr_y_xyyyzz_xxyz, \
                             tr_y_xyyyzz_xxzz, \
                             tr_y_xyyyzz_xyyy, \
                             tr_y_xyyyzz_xyyz, \
                             tr_y_xyyyzz_xyzz, \
                             tr_y_xyyyzz_xzzz, \
                             tr_y_xyyyzz_yyyy, \
                             tr_y_xyyyzz_yyyz, \
                             tr_y_xyyyzz_yyzz, \
                             tr_y_xyyyzz_yzzz, \
                             tr_y_xyyyzz_zzzz, \
                             tr_y_yyyzz_xxx,   \
                             tr_y_yyyzz_xxxx,  \
                             tr_y_yyyzz_xxxy,  \
                             tr_y_yyyzz_xxxz,  \
                             tr_y_yyyzz_xxy,   \
                             tr_y_yyyzz_xxyy,  \
                             tr_y_yyyzz_xxyz,  \
                             tr_y_yyyzz_xxz,   \
                             tr_y_yyyzz_xxzz,  \
                             tr_y_yyyzz_xyy,   \
                             tr_y_yyyzz_xyyy,  \
                             tr_y_yyyzz_xyyz,  \
                             tr_y_yyyzz_xyz,   \
                             tr_y_yyyzz_xyzz,  \
                             tr_y_yyyzz_xzz,   \
                             tr_y_yyyzz_xzzz,  \
                             tr_y_yyyzz_yyy,   \
                             tr_y_yyyzz_yyyy,  \
                             tr_y_yyyzz_yyyz,  \
                             tr_y_yyyzz_yyz,   \
                             tr_y_yyyzz_yyzz,  \
                             tr_y_yyyzz_yzz,   \
                             tr_y_yyyzz_yzzz,  \
                             tr_y_yyyzz_zzz,   \
                             tr_y_yyyzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyzz_xxxx[i] = 4.0 * tr_y_yyyzz_xxx[i] * fe_0 + tr_y_yyyzz_xxxx[i] * pa_x[i];

        tr_y_xyyyzz_xxxy[i] = 3.0 * tr_y_yyyzz_xxy[i] * fe_0 + tr_y_yyyzz_xxxy[i] * pa_x[i];

        tr_y_xyyyzz_xxxz[i] = 3.0 * tr_y_yyyzz_xxz[i] * fe_0 + tr_y_yyyzz_xxxz[i] * pa_x[i];

        tr_y_xyyyzz_xxyy[i] = 2.0 * tr_y_yyyzz_xyy[i] * fe_0 + tr_y_yyyzz_xxyy[i] * pa_x[i];

        tr_y_xyyyzz_xxyz[i] = 2.0 * tr_y_yyyzz_xyz[i] * fe_0 + tr_y_yyyzz_xxyz[i] * pa_x[i];

        tr_y_xyyyzz_xxzz[i] = 2.0 * tr_y_yyyzz_xzz[i] * fe_0 + tr_y_yyyzz_xxzz[i] * pa_x[i];

        tr_y_xyyyzz_xyyy[i] = tr_y_yyyzz_yyy[i] * fe_0 + tr_y_yyyzz_xyyy[i] * pa_x[i];

        tr_y_xyyyzz_xyyz[i] = tr_y_yyyzz_yyz[i] * fe_0 + tr_y_yyyzz_xyyz[i] * pa_x[i];

        tr_y_xyyyzz_xyzz[i] = tr_y_yyyzz_yzz[i] * fe_0 + tr_y_yyyzz_xyzz[i] * pa_x[i];

        tr_y_xyyyzz_xzzz[i] = tr_y_yyyzz_zzz[i] * fe_0 + tr_y_yyyzz_xzzz[i] * pa_x[i];

        tr_y_xyyyzz_yyyy[i] = tr_y_yyyzz_yyyy[i] * pa_x[i];

        tr_y_xyyyzz_yyyz[i] = tr_y_yyyzz_yyyz[i] * pa_x[i];

        tr_y_xyyyzz_yyzz[i] = tr_y_yyyzz_yyzz[i] * pa_x[i];

        tr_y_xyyyzz_yzzz[i] = tr_y_yyyzz_yzzz[i] * pa_x[i];

        tr_y_xyyyzz_zzzz[i] = tr_y_yyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 690-705 components of targeted buffer : IG

    auto tr_y_xyyzzz_xxxx = pbuffer.data(idx_dip_ig + 690);

    auto tr_y_xyyzzz_xxxy = pbuffer.data(idx_dip_ig + 691);

    auto tr_y_xyyzzz_xxxz = pbuffer.data(idx_dip_ig + 692);

    auto tr_y_xyyzzz_xxyy = pbuffer.data(idx_dip_ig + 693);

    auto tr_y_xyyzzz_xxyz = pbuffer.data(idx_dip_ig + 694);

    auto tr_y_xyyzzz_xxzz = pbuffer.data(idx_dip_ig + 695);

    auto tr_y_xyyzzz_xyyy = pbuffer.data(idx_dip_ig + 696);

    auto tr_y_xyyzzz_xyyz = pbuffer.data(idx_dip_ig + 697);

    auto tr_y_xyyzzz_xyzz = pbuffer.data(idx_dip_ig + 698);

    auto tr_y_xyyzzz_xzzz = pbuffer.data(idx_dip_ig + 699);

    auto tr_y_xyyzzz_yyyy = pbuffer.data(idx_dip_ig + 700);

    auto tr_y_xyyzzz_yyyz = pbuffer.data(idx_dip_ig + 701);

    auto tr_y_xyyzzz_yyzz = pbuffer.data(idx_dip_ig + 702);

    auto tr_y_xyyzzz_yzzz = pbuffer.data(idx_dip_ig + 703);

    auto tr_y_xyyzzz_zzzz = pbuffer.data(idx_dip_ig + 704);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyyzzz_xxxx, \
                             tr_y_xyyzzz_xxxy, \
                             tr_y_xyyzzz_xxxz, \
                             tr_y_xyyzzz_xxyy, \
                             tr_y_xyyzzz_xxyz, \
                             tr_y_xyyzzz_xxzz, \
                             tr_y_xyyzzz_xyyy, \
                             tr_y_xyyzzz_xyyz, \
                             tr_y_xyyzzz_xyzz, \
                             tr_y_xyyzzz_xzzz, \
                             tr_y_xyyzzz_yyyy, \
                             tr_y_xyyzzz_yyyz, \
                             tr_y_xyyzzz_yyzz, \
                             tr_y_xyyzzz_yzzz, \
                             tr_y_xyyzzz_zzzz, \
                             tr_y_yyzzz_xxx,   \
                             tr_y_yyzzz_xxxx,  \
                             tr_y_yyzzz_xxxy,  \
                             tr_y_yyzzz_xxxz,  \
                             tr_y_yyzzz_xxy,   \
                             tr_y_yyzzz_xxyy,  \
                             tr_y_yyzzz_xxyz,  \
                             tr_y_yyzzz_xxz,   \
                             tr_y_yyzzz_xxzz,  \
                             tr_y_yyzzz_xyy,   \
                             tr_y_yyzzz_xyyy,  \
                             tr_y_yyzzz_xyyz,  \
                             tr_y_yyzzz_xyz,   \
                             tr_y_yyzzz_xyzz,  \
                             tr_y_yyzzz_xzz,   \
                             tr_y_yyzzz_xzzz,  \
                             tr_y_yyzzz_yyy,   \
                             tr_y_yyzzz_yyyy,  \
                             tr_y_yyzzz_yyyz,  \
                             tr_y_yyzzz_yyz,   \
                             tr_y_yyzzz_yyzz,  \
                             tr_y_yyzzz_yzz,   \
                             tr_y_yyzzz_yzzz,  \
                             tr_y_yyzzz_zzz,   \
                             tr_y_yyzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzzz_xxxx[i] = 4.0 * tr_y_yyzzz_xxx[i] * fe_0 + tr_y_yyzzz_xxxx[i] * pa_x[i];

        tr_y_xyyzzz_xxxy[i] = 3.0 * tr_y_yyzzz_xxy[i] * fe_0 + tr_y_yyzzz_xxxy[i] * pa_x[i];

        tr_y_xyyzzz_xxxz[i] = 3.0 * tr_y_yyzzz_xxz[i] * fe_0 + tr_y_yyzzz_xxxz[i] * pa_x[i];

        tr_y_xyyzzz_xxyy[i] = 2.0 * tr_y_yyzzz_xyy[i] * fe_0 + tr_y_yyzzz_xxyy[i] * pa_x[i];

        tr_y_xyyzzz_xxyz[i] = 2.0 * tr_y_yyzzz_xyz[i] * fe_0 + tr_y_yyzzz_xxyz[i] * pa_x[i];

        tr_y_xyyzzz_xxzz[i] = 2.0 * tr_y_yyzzz_xzz[i] * fe_0 + tr_y_yyzzz_xxzz[i] * pa_x[i];

        tr_y_xyyzzz_xyyy[i] = tr_y_yyzzz_yyy[i] * fe_0 + tr_y_yyzzz_xyyy[i] * pa_x[i];

        tr_y_xyyzzz_xyyz[i] = tr_y_yyzzz_yyz[i] * fe_0 + tr_y_yyzzz_xyyz[i] * pa_x[i];

        tr_y_xyyzzz_xyzz[i] = tr_y_yyzzz_yzz[i] * fe_0 + tr_y_yyzzz_xyzz[i] * pa_x[i];

        tr_y_xyyzzz_xzzz[i] = tr_y_yyzzz_zzz[i] * fe_0 + tr_y_yyzzz_xzzz[i] * pa_x[i];

        tr_y_xyyzzz_yyyy[i] = tr_y_yyzzz_yyyy[i] * pa_x[i];

        tr_y_xyyzzz_yyyz[i] = tr_y_yyzzz_yyyz[i] * pa_x[i];

        tr_y_xyyzzz_yyzz[i] = tr_y_yyzzz_yyzz[i] * pa_x[i];

        tr_y_xyyzzz_yzzz[i] = tr_y_yyzzz_yzzz[i] * pa_x[i];

        tr_y_xyyzzz_zzzz[i] = tr_y_yyzzz_zzzz[i] * pa_x[i];
    }

    // Set up 705-720 components of targeted buffer : IG

    auto tr_y_xyzzzz_xxxx = pbuffer.data(idx_dip_ig + 705);

    auto tr_y_xyzzzz_xxxy = pbuffer.data(idx_dip_ig + 706);

    auto tr_y_xyzzzz_xxxz = pbuffer.data(idx_dip_ig + 707);

    auto tr_y_xyzzzz_xxyy = pbuffer.data(idx_dip_ig + 708);

    auto tr_y_xyzzzz_xxyz = pbuffer.data(idx_dip_ig + 709);

    auto tr_y_xyzzzz_xxzz = pbuffer.data(idx_dip_ig + 710);

    auto tr_y_xyzzzz_xyyy = pbuffer.data(idx_dip_ig + 711);

    auto tr_y_xyzzzz_xyyz = pbuffer.data(idx_dip_ig + 712);

    auto tr_y_xyzzzz_xyzz = pbuffer.data(idx_dip_ig + 713);

    auto tr_y_xyzzzz_xzzz = pbuffer.data(idx_dip_ig + 714);

    auto tr_y_xyzzzz_yyyy = pbuffer.data(idx_dip_ig + 715);

    auto tr_y_xyzzzz_yyyz = pbuffer.data(idx_dip_ig + 716);

    auto tr_y_xyzzzz_yyzz = pbuffer.data(idx_dip_ig + 717);

    auto tr_y_xyzzzz_yzzz = pbuffer.data(idx_dip_ig + 718);

    auto tr_y_xyzzzz_zzzz = pbuffer.data(idx_dip_ig + 719);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyzzzz_xxxx, \
                             tr_y_xyzzzz_xxxy, \
                             tr_y_xyzzzz_xxxz, \
                             tr_y_xyzzzz_xxyy, \
                             tr_y_xyzzzz_xxyz, \
                             tr_y_xyzzzz_xxzz, \
                             tr_y_xyzzzz_xyyy, \
                             tr_y_xyzzzz_xyyz, \
                             tr_y_xyzzzz_xyzz, \
                             tr_y_xyzzzz_xzzz, \
                             tr_y_xyzzzz_yyyy, \
                             tr_y_xyzzzz_yyyz, \
                             tr_y_xyzzzz_yyzz, \
                             tr_y_xyzzzz_yzzz, \
                             tr_y_xyzzzz_zzzz, \
                             tr_y_yzzzz_xxx,   \
                             tr_y_yzzzz_xxxx,  \
                             tr_y_yzzzz_xxxy,  \
                             tr_y_yzzzz_xxxz,  \
                             tr_y_yzzzz_xxy,   \
                             tr_y_yzzzz_xxyy,  \
                             tr_y_yzzzz_xxyz,  \
                             tr_y_yzzzz_xxz,   \
                             tr_y_yzzzz_xxzz,  \
                             tr_y_yzzzz_xyy,   \
                             tr_y_yzzzz_xyyy,  \
                             tr_y_yzzzz_xyyz,  \
                             tr_y_yzzzz_xyz,   \
                             tr_y_yzzzz_xyzz,  \
                             tr_y_yzzzz_xzz,   \
                             tr_y_yzzzz_xzzz,  \
                             tr_y_yzzzz_yyy,   \
                             tr_y_yzzzz_yyyy,  \
                             tr_y_yzzzz_yyyz,  \
                             tr_y_yzzzz_yyz,   \
                             tr_y_yzzzz_yyzz,  \
                             tr_y_yzzzz_yzz,   \
                             tr_y_yzzzz_yzzz,  \
                             tr_y_yzzzz_zzz,   \
                             tr_y_yzzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzzz_xxxx[i] = 4.0 * tr_y_yzzzz_xxx[i] * fe_0 + tr_y_yzzzz_xxxx[i] * pa_x[i];

        tr_y_xyzzzz_xxxy[i] = 3.0 * tr_y_yzzzz_xxy[i] * fe_0 + tr_y_yzzzz_xxxy[i] * pa_x[i];

        tr_y_xyzzzz_xxxz[i] = 3.0 * tr_y_yzzzz_xxz[i] * fe_0 + tr_y_yzzzz_xxxz[i] * pa_x[i];

        tr_y_xyzzzz_xxyy[i] = 2.0 * tr_y_yzzzz_xyy[i] * fe_0 + tr_y_yzzzz_xxyy[i] * pa_x[i];

        tr_y_xyzzzz_xxyz[i] = 2.0 * tr_y_yzzzz_xyz[i] * fe_0 + tr_y_yzzzz_xxyz[i] * pa_x[i];

        tr_y_xyzzzz_xxzz[i] = 2.0 * tr_y_yzzzz_xzz[i] * fe_0 + tr_y_yzzzz_xxzz[i] * pa_x[i];

        tr_y_xyzzzz_xyyy[i] = tr_y_yzzzz_yyy[i] * fe_0 + tr_y_yzzzz_xyyy[i] * pa_x[i];

        tr_y_xyzzzz_xyyz[i] = tr_y_yzzzz_yyz[i] * fe_0 + tr_y_yzzzz_xyyz[i] * pa_x[i];

        tr_y_xyzzzz_xyzz[i] = tr_y_yzzzz_yzz[i] * fe_0 + tr_y_yzzzz_xyzz[i] * pa_x[i];

        tr_y_xyzzzz_xzzz[i] = tr_y_yzzzz_zzz[i] * fe_0 + tr_y_yzzzz_xzzz[i] * pa_x[i];

        tr_y_xyzzzz_yyyy[i] = tr_y_yzzzz_yyyy[i] * pa_x[i];

        tr_y_xyzzzz_yyyz[i] = tr_y_yzzzz_yyyz[i] * pa_x[i];

        tr_y_xyzzzz_yyzz[i] = tr_y_yzzzz_yyzz[i] * pa_x[i];

        tr_y_xyzzzz_yzzz[i] = tr_y_yzzzz_yzzz[i] * pa_x[i];

        tr_y_xyzzzz_zzzz[i] = tr_y_yzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 720-735 components of targeted buffer : IG

    auto tr_y_xzzzzz_xxxx = pbuffer.data(idx_dip_ig + 720);

    auto tr_y_xzzzzz_xxxy = pbuffer.data(idx_dip_ig + 721);

    auto tr_y_xzzzzz_xxxz = pbuffer.data(idx_dip_ig + 722);

    auto tr_y_xzzzzz_xxyy = pbuffer.data(idx_dip_ig + 723);

    auto tr_y_xzzzzz_xxyz = pbuffer.data(idx_dip_ig + 724);

    auto tr_y_xzzzzz_xxzz = pbuffer.data(idx_dip_ig + 725);

    auto tr_y_xzzzzz_xyyy = pbuffer.data(idx_dip_ig + 726);

    auto tr_y_xzzzzz_xyyz = pbuffer.data(idx_dip_ig + 727);

    auto tr_y_xzzzzz_xyzz = pbuffer.data(idx_dip_ig + 728);

    auto tr_y_xzzzzz_xzzz = pbuffer.data(idx_dip_ig + 729);

    auto tr_y_xzzzzz_yyyy = pbuffer.data(idx_dip_ig + 730);

    auto tr_y_xzzzzz_yyyz = pbuffer.data(idx_dip_ig + 731);

    auto tr_y_xzzzzz_yyzz = pbuffer.data(idx_dip_ig + 732);

    auto tr_y_xzzzzz_yzzz = pbuffer.data(idx_dip_ig + 733);

    auto tr_y_xzzzzz_zzzz = pbuffer.data(idx_dip_ig + 734);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xzzzzz_xxxx, \
                             tr_y_xzzzzz_xxxy, \
                             tr_y_xzzzzz_xxxz, \
                             tr_y_xzzzzz_xxyy, \
                             tr_y_xzzzzz_xxyz, \
                             tr_y_xzzzzz_xxzz, \
                             tr_y_xzzzzz_xyyy, \
                             tr_y_xzzzzz_xyyz, \
                             tr_y_xzzzzz_xyzz, \
                             tr_y_xzzzzz_xzzz, \
                             tr_y_xzzzzz_yyyy, \
                             tr_y_xzzzzz_yyyz, \
                             tr_y_xzzzzz_yyzz, \
                             tr_y_xzzzzz_yzzz, \
                             tr_y_xzzzzz_zzzz, \
                             tr_y_zzzzz_xxx,   \
                             tr_y_zzzzz_xxxx,  \
                             tr_y_zzzzz_xxxy,  \
                             tr_y_zzzzz_xxxz,  \
                             tr_y_zzzzz_xxy,   \
                             tr_y_zzzzz_xxyy,  \
                             tr_y_zzzzz_xxyz,  \
                             tr_y_zzzzz_xxz,   \
                             tr_y_zzzzz_xxzz,  \
                             tr_y_zzzzz_xyy,   \
                             tr_y_zzzzz_xyyy,  \
                             tr_y_zzzzz_xyyz,  \
                             tr_y_zzzzz_xyz,   \
                             tr_y_zzzzz_xyzz,  \
                             tr_y_zzzzz_xzz,   \
                             tr_y_zzzzz_xzzz,  \
                             tr_y_zzzzz_yyy,   \
                             tr_y_zzzzz_yyyy,  \
                             tr_y_zzzzz_yyyz,  \
                             tr_y_zzzzz_yyz,   \
                             tr_y_zzzzz_yyzz,  \
                             tr_y_zzzzz_yzz,   \
                             tr_y_zzzzz_yzzz,  \
                             tr_y_zzzzz_zzz,   \
                             tr_y_zzzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzzz_xxxx[i] = 4.0 * tr_y_zzzzz_xxx[i] * fe_0 + tr_y_zzzzz_xxxx[i] * pa_x[i];

        tr_y_xzzzzz_xxxy[i] = 3.0 * tr_y_zzzzz_xxy[i] * fe_0 + tr_y_zzzzz_xxxy[i] * pa_x[i];

        tr_y_xzzzzz_xxxz[i] = 3.0 * tr_y_zzzzz_xxz[i] * fe_0 + tr_y_zzzzz_xxxz[i] * pa_x[i];

        tr_y_xzzzzz_xxyy[i] = 2.0 * tr_y_zzzzz_xyy[i] * fe_0 + tr_y_zzzzz_xxyy[i] * pa_x[i];

        tr_y_xzzzzz_xxyz[i] = 2.0 * tr_y_zzzzz_xyz[i] * fe_0 + tr_y_zzzzz_xxyz[i] * pa_x[i];

        tr_y_xzzzzz_xxzz[i] = 2.0 * tr_y_zzzzz_xzz[i] * fe_0 + tr_y_zzzzz_xxzz[i] * pa_x[i];

        tr_y_xzzzzz_xyyy[i] = tr_y_zzzzz_yyy[i] * fe_0 + tr_y_zzzzz_xyyy[i] * pa_x[i];

        tr_y_xzzzzz_xyyz[i] = tr_y_zzzzz_yyz[i] * fe_0 + tr_y_zzzzz_xyyz[i] * pa_x[i];

        tr_y_xzzzzz_xyzz[i] = tr_y_zzzzz_yzz[i] * fe_0 + tr_y_zzzzz_xyzz[i] * pa_x[i];

        tr_y_xzzzzz_xzzz[i] = tr_y_zzzzz_zzz[i] * fe_0 + tr_y_zzzzz_xzzz[i] * pa_x[i];

        tr_y_xzzzzz_yyyy[i] = tr_y_zzzzz_yyyy[i] * pa_x[i];

        tr_y_xzzzzz_yyyz[i] = tr_y_zzzzz_yyyz[i] * pa_x[i];

        tr_y_xzzzzz_yyzz[i] = tr_y_zzzzz_yyzz[i] * pa_x[i];

        tr_y_xzzzzz_yzzz[i] = tr_y_zzzzz_yzzz[i] * pa_x[i];

        tr_y_xzzzzz_zzzz[i] = tr_y_zzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 735-750 components of targeted buffer : IG

    auto tr_y_yyyyyy_xxxx = pbuffer.data(idx_dip_ig + 735);

    auto tr_y_yyyyyy_xxxy = pbuffer.data(idx_dip_ig + 736);

    auto tr_y_yyyyyy_xxxz = pbuffer.data(idx_dip_ig + 737);

    auto tr_y_yyyyyy_xxyy = pbuffer.data(idx_dip_ig + 738);

    auto tr_y_yyyyyy_xxyz = pbuffer.data(idx_dip_ig + 739);

    auto tr_y_yyyyyy_xxzz = pbuffer.data(idx_dip_ig + 740);

    auto tr_y_yyyyyy_xyyy = pbuffer.data(idx_dip_ig + 741);

    auto tr_y_yyyyyy_xyyz = pbuffer.data(idx_dip_ig + 742);

    auto tr_y_yyyyyy_xyzz = pbuffer.data(idx_dip_ig + 743);

    auto tr_y_yyyyyy_xzzz = pbuffer.data(idx_dip_ig + 744);

    auto tr_y_yyyyyy_yyyy = pbuffer.data(idx_dip_ig + 745);

    auto tr_y_yyyyyy_yyyz = pbuffer.data(idx_dip_ig + 746);

    auto tr_y_yyyyyy_yyzz = pbuffer.data(idx_dip_ig + 747);

    auto tr_y_yyyyyy_yzzz = pbuffer.data(idx_dip_ig + 748);

    auto tr_y_yyyyyy_zzzz = pbuffer.data(idx_dip_ig + 749);

#pragma omp simd aligned(pa_y,                 \
                             tr_y_yyyy_xxxx,   \
                             tr_y_yyyy_xxxy,   \
                             tr_y_yyyy_xxxz,   \
                             tr_y_yyyy_xxyy,   \
                             tr_y_yyyy_xxyz,   \
                             tr_y_yyyy_xxzz,   \
                             tr_y_yyyy_xyyy,   \
                             tr_y_yyyy_xyyz,   \
                             tr_y_yyyy_xyzz,   \
                             tr_y_yyyy_xzzz,   \
                             tr_y_yyyy_yyyy,   \
                             tr_y_yyyy_yyyz,   \
                             tr_y_yyyy_yyzz,   \
                             tr_y_yyyy_yzzz,   \
                             tr_y_yyyy_zzzz,   \
                             tr_y_yyyyy_xxx,   \
                             tr_y_yyyyy_xxxx,  \
                             tr_y_yyyyy_xxxy,  \
                             tr_y_yyyyy_xxxz,  \
                             tr_y_yyyyy_xxy,   \
                             tr_y_yyyyy_xxyy,  \
                             tr_y_yyyyy_xxyz,  \
                             tr_y_yyyyy_xxz,   \
                             tr_y_yyyyy_xxzz,  \
                             tr_y_yyyyy_xyy,   \
                             tr_y_yyyyy_xyyy,  \
                             tr_y_yyyyy_xyyz,  \
                             tr_y_yyyyy_xyz,   \
                             tr_y_yyyyy_xyzz,  \
                             tr_y_yyyyy_xzz,   \
                             tr_y_yyyyy_xzzz,  \
                             tr_y_yyyyy_yyy,   \
                             tr_y_yyyyy_yyyy,  \
                             tr_y_yyyyy_yyyz,  \
                             tr_y_yyyyy_yyz,   \
                             tr_y_yyyyy_yyzz,  \
                             tr_y_yyyyy_yzz,   \
                             tr_y_yyyyy_yzzz,  \
                             tr_y_yyyyy_zzz,   \
                             tr_y_yyyyy_zzzz,  \
                             tr_y_yyyyyy_xxxx, \
                             tr_y_yyyyyy_xxxy, \
                             tr_y_yyyyyy_xxxz, \
                             tr_y_yyyyyy_xxyy, \
                             tr_y_yyyyyy_xxyz, \
                             tr_y_yyyyyy_xxzz, \
                             tr_y_yyyyyy_xyyy, \
                             tr_y_yyyyyy_xyyz, \
                             tr_y_yyyyyy_xyzz, \
                             tr_y_yyyyyy_xzzz, \
                             tr_y_yyyyyy_yyyy, \
                             tr_y_yyyyyy_yyyz, \
                             tr_y_yyyyyy_yyzz, \
                             tr_y_yyyyyy_yzzz, \
                             tr_y_yyyyyy_zzzz, \
                             ts_yyyyy_xxxx,    \
                             ts_yyyyy_xxxy,    \
                             ts_yyyyy_xxxz,    \
                             ts_yyyyy_xxyy,    \
                             ts_yyyyy_xxyz,    \
                             ts_yyyyy_xxzz,    \
                             ts_yyyyy_xyyy,    \
                             ts_yyyyy_xyyz,    \
                             ts_yyyyy_xyzz,    \
                             ts_yyyyy_xzzz,    \
                             ts_yyyyy_yyyy,    \
                             ts_yyyyy_yyyz,    \
                             ts_yyyyy_yyzz,    \
                             ts_yyyyy_yzzz,    \
                             ts_yyyyy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyy_xxxx[i] = 5.0 * tr_y_yyyy_xxxx[i] * fe_0 + ts_yyyyy_xxxx[i] * fe_0 + tr_y_yyyyy_xxxx[i] * pa_y[i];

        tr_y_yyyyyy_xxxy[i] = 5.0 * tr_y_yyyy_xxxy[i] * fe_0 + tr_y_yyyyy_xxx[i] * fe_0 + ts_yyyyy_xxxy[i] * fe_0 + tr_y_yyyyy_xxxy[i] * pa_y[i];

        tr_y_yyyyyy_xxxz[i] = 5.0 * tr_y_yyyy_xxxz[i] * fe_0 + ts_yyyyy_xxxz[i] * fe_0 + tr_y_yyyyy_xxxz[i] * pa_y[i];

        tr_y_yyyyyy_xxyy[i] =
            5.0 * tr_y_yyyy_xxyy[i] * fe_0 + 2.0 * tr_y_yyyyy_xxy[i] * fe_0 + ts_yyyyy_xxyy[i] * fe_0 + tr_y_yyyyy_xxyy[i] * pa_y[i];

        tr_y_yyyyyy_xxyz[i] = 5.0 * tr_y_yyyy_xxyz[i] * fe_0 + tr_y_yyyyy_xxz[i] * fe_0 + ts_yyyyy_xxyz[i] * fe_0 + tr_y_yyyyy_xxyz[i] * pa_y[i];

        tr_y_yyyyyy_xxzz[i] = 5.0 * tr_y_yyyy_xxzz[i] * fe_0 + ts_yyyyy_xxzz[i] * fe_0 + tr_y_yyyyy_xxzz[i] * pa_y[i];

        tr_y_yyyyyy_xyyy[i] =
            5.0 * tr_y_yyyy_xyyy[i] * fe_0 + 3.0 * tr_y_yyyyy_xyy[i] * fe_0 + ts_yyyyy_xyyy[i] * fe_0 + tr_y_yyyyy_xyyy[i] * pa_y[i];

        tr_y_yyyyyy_xyyz[i] =
            5.0 * tr_y_yyyy_xyyz[i] * fe_0 + 2.0 * tr_y_yyyyy_xyz[i] * fe_0 + ts_yyyyy_xyyz[i] * fe_0 + tr_y_yyyyy_xyyz[i] * pa_y[i];

        tr_y_yyyyyy_xyzz[i] = 5.0 * tr_y_yyyy_xyzz[i] * fe_0 + tr_y_yyyyy_xzz[i] * fe_0 + ts_yyyyy_xyzz[i] * fe_0 + tr_y_yyyyy_xyzz[i] * pa_y[i];

        tr_y_yyyyyy_xzzz[i] = 5.0 * tr_y_yyyy_xzzz[i] * fe_0 + ts_yyyyy_xzzz[i] * fe_0 + tr_y_yyyyy_xzzz[i] * pa_y[i];

        tr_y_yyyyyy_yyyy[i] =
            5.0 * tr_y_yyyy_yyyy[i] * fe_0 + 4.0 * tr_y_yyyyy_yyy[i] * fe_0 + ts_yyyyy_yyyy[i] * fe_0 + tr_y_yyyyy_yyyy[i] * pa_y[i];

        tr_y_yyyyyy_yyyz[i] =
            5.0 * tr_y_yyyy_yyyz[i] * fe_0 + 3.0 * tr_y_yyyyy_yyz[i] * fe_0 + ts_yyyyy_yyyz[i] * fe_0 + tr_y_yyyyy_yyyz[i] * pa_y[i];

        tr_y_yyyyyy_yyzz[i] =
            5.0 * tr_y_yyyy_yyzz[i] * fe_0 + 2.0 * tr_y_yyyyy_yzz[i] * fe_0 + ts_yyyyy_yyzz[i] * fe_0 + tr_y_yyyyy_yyzz[i] * pa_y[i];

        tr_y_yyyyyy_yzzz[i] = 5.0 * tr_y_yyyy_yzzz[i] * fe_0 + tr_y_yyyyy_zzz[i] * fe_0 + ts_yyyyy_yzzz[i] * fe_0 + tr_y_yyyyy_yzzz[i] * pa_y[i];

        tr_y_yyyyyy_zzzz[i] = 5.0 * tr_y_yyyy_zzzz[i] * fe_0 + ts_yyyyy_zzzz[i] * fe_0 + tr_y_yyyyy_zzzz[i] * pa_y[i];
    }

    // Set up 750-765 components of targeted buffer : IG

    auto tr_y_yyyyyz_xxxx = pbuffer.data(idx_dip_ig + 750);

    auto tr_y_yyyyyz_xxxy = pbuffer.data(idx_dip_ig + 751);

    auto tr_y_yyyyyz_xxxz = pbuffer.data(idx_dip_ig + 752);

    auto tr_y_yyyyyz_xxyy = pbuffer.data(idx_dip_ig + 753);

    auto tr_y_yyyyyz_xxyz = pbuffer.data(idx_dip_ig + 754);

    auto tr_y_yyyyyz_xxzz = pbuffer.data(idx_dip_ig + 755);

    auto tr_y_yyyyyz_xyyy = pbuffer.data(idx_dip_ig + 756);

    auto tr_y_yyyyyz_xyyz = pbuffer.data(idx_dip_ig + 757);

    auto tr_y_yyyyyz_xyzz = pbuffer.data(idx_dip_ig + 758);

    auto tr_y_yyyyyz_xzzz = pbuffer.data(idx_dip_ig + 759);

    auto tr_y_yyyyyz_yyyy = pbuffer.data(idx_dip_ig + 760);

    auto tr_y_yyyyyz_yyyz = pbuffer.data(idx_dip_ig + 761);

    auto tr_y_yyyyyz_yyzz = pbuffer.data(idx_dip_ig + 762);

    auto tr_y_yyyyyz_yzzz = pbuffer.data(idx_dip_ig + 763);

    auto tr_y_yyyyyz_zzzz = pbuffer.data(idx_dip_ig + 764);

#pragma omp simd aligned(pa_z,                 \
                             tr_y_yyyyy_xxx,   \
                             tr_y_yyyyy_xxxx,  \
                             tr_y_yyyyy_xxxy,  \
                             tr_y_yyyyy_xxxz,  \
                             tr_y_yyyyy_xxy,   \
                             tr_y_yyyyy_xxyy,  \
                             tr_y_yyyyy_xxyz,  \
                             tr_y_yyyyy_xxz,   \
                             tr_y_yyyyy_xxzz,  \
                             tr_y_yyyyy_xyy,   \
                             tr_y_yyyyy_xyyy,  \
                             tr_y_yyyyy_xyyz,  \
                             tr_y_yyyyy_xyz,   \
                             tr_y_yyyyy_xyzz,  \
                             tr_y_yyyyy_xzz,   \
                             tr_y_yyyyy_xzzz,  \
                             tr_y_yyyyy_yyy,   \
                             tr_y_yyyyy_yyyy,  \
                             tr_y_yyyyy_yyyz,  \
                             tr_y_yyyyy_yyz,   \
                             tr_y_yyyyy_yyzz,  \
                             tr_y_yyyyy_yzz,   \
                             tr_y_yyyyy_yzzz,  \
                             tr_y_yyyyy_zzz,   \
                             tr_y_yyyyy_zzzz,  \
                             tr_y_yyyyyz_xxxx, \
                             tr_y_yyyyyz_xxxy, \
                             tr_y_yyyyyz_xxxz, \
                             tr_y_yyyyyz_xxyy, \
                             tr_y_yyyyyz_xxyz, \
                             tr_y_yyyyyz_xxzz, \
                             tr_y_yyyyyz_xyyy, \
                             tr_y_yyyyyz_xyyz, \
                             tr_y_yyyyyz_xyzz, \
                             tr_y_yyyyyz_xzzz, \
                             tr_y_yyyyyz_yyyy, \
                             tr_y_yyyyyz_yyyz, \
                             tr_y_yyyyyz_yyzz, \
                             tr_y_yyyyyz_yzzz, \
                             tr_y_yyyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyz_xxxx[i] = tr_y_yyyyy_xxxx[i] * pa_z[i];

        tr_y_yyyyyz_xxxy[i] = tr_y_yyyyy_xxxy[i] * pa_z[i];

        tr_y_yyyyyz_xxxz[i] = tr_y_yyyyy_xxx[i] * fe_0 + tr_y_yyyyy_xxxz[i] * pa_z[i];

        tr_y_yyyyyz_xxyy[i] = tr_y_yyyyy_xxyy[i] * pa_z[i];

        tr_y_yyyyyz_xxyz[i] = tr_y_yyyyy_xxy[i] * fe_0 + tr_y_yyyyy_xxyz[i] * pa_z[i];

        tr_y_yyyyyz_xxzz[i] = 2.0 * tr_y_yyyyy_xxz[i] * fe_0 + tr_y_yyyyy_xxzz[i] * pa_z[i];

        tr_y_yyyyyz_xyyy[i] = tr_y_yyyyy_xyyy[i] * pa_z[i];

        tr_y_yyyyyz_xyyz[i] = tr_y_yyyyy_xyy[i] * fe_0 + tr_y_yyyyy_xyyz[i] * pa_z[i];

        tr_y_yyyyyz_xyzz[i] = 2.0 * tr_y_yyyyy_xyz[i] * fe_0 + tr_y_yyyyy_xyzz[i] * pa_z[i];

        tr_y_yyyyyz_xzzz[i] = 3.0 * tr_y_yyyyy_xzz[i] * fe_0 + tr_y_yyyyy_xzzz[i] * pa_z[i];

        tr_y_yyyyyz_yyyy[i] = tr_y_yyyyy_yyyy[i] * pa_z[i];

        tr_y_yyyyyz_yyyz[i] = tr_y_yyyyy_yyy[i] * fe_0 + tr_y_yyyyy_yyyz[i] * pa_z[i];

        tr_y_yyyyyz_yyzz[i] = 2.0 * tr_y_yyyyy_yyz[i] * fe_0 + tr_y_yyyyy_yyzz[i] * pa_z[i];

        tr_y_yyyyyz_yzzz[i] = 3.0 * tr_y_yyyyy_yzz[i] * fe_0 + tr_y_yyyyy_yzzz[i] * pa_z[i];

        tr_y_yyyyyz_zzzz[i] = 4.0 * tr_y_yyyyy_zzz[i] * fe_0 + tr_y_yyyyy_zzzz[i] * pa_z[i];
    }

    // Set up 765-780 components of targeted buffer : IG

    auto tr_y_yyyyzz_xxxx = pbuffer.data(idx_dip_ig + 765);

    auto tr_y_yyyyzz_xxxy = pbuffer.data(idx_dip_ig + 766);

    auto tr_y_yyyyzz_xxxz = pbuffer.data(idx_dip_ig + 767);

    auto tr_y_yyyyzz_xxyy = pbuffer.data(idx_dip_ig + 768);

    auto tr_y_yyyyzz_xxyz = pbuffer.data(idx_dip_ig + 769);

    auto tr_y_yyyyzz_xxzz = pbuffer.data(idx_dip_ig + 770);

    auto tr_y_yyyyzz_xyyy = pbuffer.data(idx_dip_ig + 771);

    auto tr_y_yyyyzz_xyyz = pbuffer.data(idx_dip_ig + 772);

    auto tr_y_yyyyzz_xyzz = pbuffer.data(idx_dip_ig + 773);

    auto tr_y_yyyyzz_xzzz = pbuffer.data(idx_dip_ig + 774);

    auto tr_y_yyyyzz_yyyy = pbuffer.data(idx_dip_ig + 775);

    auto tr_y_yyyyzz_yyyz = pbuffer.data(idx_dip_ig + 776);

    auto tr_y_yyyyzz_yyzz = pbuffer.data(idx_dip_ig + 777);

    auto tr_y_yyyyzz_yzzz = pbuffer.data(idx_dip_ig + 778);

    auto tr_y_yyyyzz_zzzz = pbuffer.data(idx_dip_ig + 779);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yyyy_xxxx,   \
                             tr_y_yyyy_xxxy,   \
                             tr_y_yyyy_xxyy,   \
                             tr_y_yyyy_xxyz,   \
                             tr_y_yyyy_xyyy,   \
                             tr_y_yyyy_xyyz,   \
                             tr_y_yyyy_xyzz,   \
                             tr_y_yyyy_yyyy,   \
                             tr_y_yyyy_yyyz,   \
                             tr_y_yyyy_yyzz,   \
                             tr_y_yyyy_yzzz,   \
                             tr_y_yyyyz_xxxx,  \
                             tr_y_yyyyz_xxxy,  \
                             tr_y_yyyyz_xxy,   \
                             tr_y_yyyyz_xxyy,  \
                             tr_y_yyyyz_xxyz,  \
                             tr_y_yyyyz_xyy,   \
                             tr_y_yyyyz_xyyy,  \
                             tr_y_yyyyz_xyyz,  \
                             tr_y_yyyyz_xyz,   \
                             tr_y_yyyyz_xyzz,  \
                             tr_y_yyyyz_yyy,   \
                             tr_y_yyyyz_yyyy,  \
                             tr_y_yyyyz_yyyz,  \
                             tr_y_yyyyz_yyz,   \
                             tr_y_yyyyz_yyzz,  \
                             tr_y_yyyyz_yzz,   \
                             tr_y_yyyyz_yzzz,  \
                             tr_y_yyyyzz_xxxx, \
                             tr_y_yyyyzz_xxxy, \
                             tr_y_yyyyzz_xxxz, \
                             tr_y_yyyyzz_xxyy, \
                             tr_y_yyyyzz_xxyz, \
                             tr_y_yyyyzz_xxzz, \
                             tr_y_yyyyzz_xyyy, \
                             tr_y_yyyyzz_xyyz, \
                             tr_y_yyyyzz_xyzz, \
                             tr_y_yyyyzz_xzzz, \
                             tr_y_yyyyzz_yyyy, \
                             tr_y_yyyyzz_yyyz, \
                             tr_y_yyyyzz_yyzz, \
                             tr_y_yyyyzz_yzzz, \
                             tr_y_yyyyzz_zzzz, \
                             tr_y_yyyzz_xxxz,  \
                             tr_y_yyyzz_xxzz,  \
                             tr_y_yyyzz_xzzz,  \
                             tr_y_yyyzz_zzzz,  \
                             tr_y_yyzz_xxxz,   \
                             tr_y_yyzz_xxzz,   \
                             tr_y_yyzz_xzzz,   \
                             tr_y_yyzz_zzzz,   \
                             ts_yyyzz_xxxz,    \
                             ts_yyyzz_xxzz,    \
                             ts_yyyzz_xzzz,    \
                             ts_yyyzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyzz_xxxx[i] = tr_y_yyyy_xxxx[i] * fe_0 + tr_y_yyyyz_xxxx[i] * pa_z[i];

        tr_y_yyyyzz_xxxy[i] = tr_y_yyyy_xxxy[i] * fe_0 + tr_y_yyyyz_xxxy[i] * pa_z[i];

        tr_y_yyyyzz_xxxz[i] = 3.0 * tr_y_yyzz_xxxz[i] * fe_0 + ts_yyyzz_xxxz[i] * fe_0 + tr_y_yyyzz_xxxz[i] * pa_y[i];

        tr_y_yyyyzz_xxyy[i] = tr_y_yyyy_xxyy[i] * fe_0 + tr_y_yyyyz_xxyy[i] * pa_z[i];

        tr_y_yyyyzz_xxyz[i] = tr_y_yyyy_xxyz[i] * fe_0 + tr_y_yyyyz_xxy[i] * fe_0 + tr_y_yyyyz_xxyz[i] * pa_z[i];

        tr_y_yyyyzz_xxzz[i] = 3.0 * tr_y_yyzz_xxzz[i] * fe_0 + ts_yyyzz_xxzz[i] * fe_0 + tr_y_yyyzz_xxzz[i] * pa_y[i];

        tr_y_yyyyzz_xyyy[i] = tr_y_yyyy_xyyy[i] * fe_0 + tr_y_yyyyz_xyyy[i] * pa_z[i];

        tr_y_yyyyzz_xyyz[i] = tr_y_yyyy_xyyz[i] * fe_0 + tr_y_yyyyz_xyy[i] * fe_0 + tr_y_yyyyz_xyyz[i] * pa_z[i];

        tr_y_yyyyzz_xyzz[i] = tr_y_yyyy_xyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_xyz[i] * fe_0 + tr_y_yyyyz_xyzz[i] * pa_z[i];

        tr_y_yyyyzz_xzzz[i] = 3.0 * tr_y_yyzz_xzzz[i] * fe_0 + ts_yyyzz_xzzz[i] * fe_0 + tr_y_yyyzz_xzzz[i] * pa_y[i];

        tr_y_yyyyzz_yyyy[i] = tr_y_yyyy_yyyy[i] * fe_0 + tr_y_yyyyz_yyyy[i] * pa_z[i];

        tr_y_yyyyzz_yyyz[i] = tr_y_yyyy_yyyz[i] * fe_0 + tr_y_yyyyz_yyy[i] * fe_0 + tr_y_yyyyz_yyyz[i] * pa_z[i];

        tr_y_yyyyzz_yyzz[i] = tr_y_yyyy_yyzz[i] * fe_0 + 2.0 * tr_y_yyyyz_yyz[i] * fe_0 + tr_y_yyyyz_yyzz[i] * pa_z[i];

        tr_y_yyyyzz_yzzz[i] = tr_y_yyyy_yzzz[i] * fe_0 + 3.0 * tr_y_yyyyz_yzz[i] * fe_0 + tr_y_yyyyz_yzzz[i] * pa_z[i];

        tr_y_yyyyzz_zzzz[i] = 3.0 * tr_y_yyzz_zzzz[i] * fe_0 + ts_yyyzz_zzzz[i] * fe_0 + tr_y_yyyzz_zzzz[i] * pa_y[i];
    }

    // Set up 780-795 components of targeted buffer : IG

    auto tr_y_yyyzzz_xxxx = pbuffer.data(idx_dip_ig + 780);

    auto tr_y_yyyzzz_xxxy = pbuffer.data(idx_dip_ig + 781);

    auto tr_y_yyyzzz_xxxz = pbuffer.data(idx_dip_ig + 782);

    auto tr_y_yyyzzz_xxyy = pbuffer.data(idx_dip_ig + 783);

    auto tr_y_yyyzzz_xxyz = pbuffer.data(idx_dip_ig + 784);

    auto tr_y_yyyzzz_xxzz = pbuffer.data(idx_dip_ig + 785);

    auto tr_y_yyyzzz_xyyy = pbuffer.data(idx_dip_ig + 786);

    auto tr_y_yyyzzz_xyyz = pbuffer.data(idx_dip_ig + 787);

    auto tr_y_yyyzzz_xyzz = pbuffer.data(idx_dip_ig + 788);

    auto tr_y_yyyzzz_xzzz = pbuffer.data(idx_dip_ig + 789);

    auto tr_y_yyyzzz_yyyy = pbuffer.data(idx_dip_ig + 790);

    auto tr_y_yyyzzz_yyyz = pbuffer.data(idx_dip_ig + 791);

    auto tr_y_yyyzzz_yyzz = pbuffer.data(idx_dip_ig + 792);

    auto tr_y_yyyzzz_yzzz = pbuffer.data(idx_dip_ig + 793);

    auto tr_y_yyyzzz_zzzz = pbuffer.data(idx_dip_ig + 794);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yyyz_xxxx,   \
                             tr_y_yyyz_xxxy,   \
                             tr_y_yyyz_xxyy,   \
                             tr_y_yyyz_xxyz,   \
                             tr_y_yyyz_xyyy,   \
                             tr_y_yyyz_xyyz,   \
                             tr_y_yyyz_xyzz,   \
                             tr_y_yyyz_yyyy,   \
                             tr_y_yyyz_yyyz,   \
                             tr_y_yyyz_yyzz,   \
                             tr_y_yyyz_yzzz,   \
                             tr_y_yyyzz_xxxx,  \
                             tr_y_yyyzz_xxxy,  \
                             tr_y_yyyzz_xxy,   \
                             tr_y_yyyzz_xxyy,  \
                             tr_y_yyyzz_xxyz,  \
                             tr_y_yyyzz_xyy,   \
                             tr_y_yyyzz_xyyy,  \
                             tr_y_yyyzz_xyyz,  \
                             tr_y_yyyzz_xyz,   \
                             tr_y_yyyzz_xyzz,  \
                             tr_y_yyyzz_yyy,   \
                             tr_y_yyyzz_yyyy,  \
                             tr_y_yyyzz_yyyz,  \
                             tr_y_yyyzz_yyz,   \
                             tr_y_yyyzz_yyzz,  \
                             tr_y_yyyzz_yzz,   \
                             tr_y_yyyzz_yzzz,  \
                             tr_y_yyyzzz_xxxx, \
                             tr_y_yyyzzz_xxxy, \
                             tr_y_yyyzzz_xxxz, \
                             tr_y_yyyzzz_xxyy, \
                             tr_y_yyyzzz_xxyz, \
                             tr_y_yyyzzz_xxzz, \
                             tr_y_yyyzzz_xyyy, \
                             tr_y_yyyzzz_xyyz, \
                             tr_y_yyyzzz_xyzz, \
                             tr_y_yyyzzz_xzzz, \
                             tr_y_yyyzzz_yyyy, \
                             tr_y_yyyzzz_yyyz, \
                             tr_y_yyyzzz_yyzz, \
                             tr_y_yyyzzz_yzzz, \
                             tr_y_yyyzzz_zzzz, \
                             tr_y_yyzzz_xxxz,  \
                             tr_y_yyzzz_xxzz,  \
                             tr_y_yyzzz_xzzz,  \
                             tr_y_yyzzz_zzzz,  \
                             tr_y_yzzz_xxxz,   \
                             tr_y_yzzz_xxzz,   \
                             tr_y_yzzz_xzzz,   \
                             tr_y_yzzz_zzzz,   \
                             ts_yyzzz_xxxz,    \
                             ts_yyzzz_xxzz,    \
                             ts_yyzzz_xzzz,    \
                             ts_yyzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzzz_xxxx[i] = 2.0 * tr_y_yyyz_xxxx[i] * fe_0 + tr_y_yyyzz_xxxx[i] * pa_z[i];

        tr_y_yyyzzz_xxxy[i] = 2.0 * tr_y_yyyz_xxxy[i] * fe_0 + tr_y_yyyzz_xxxy[i] * pa_z[i];

        tr_y_yyyzzz_xxxz[i] = 2.0 * tr_y_yzzz_xxxz[i] * fe_0 + ts_yyzzz_xxxz[i] * fe_0 + tr_y_yyzzz_xxxz[i] * pa_y[i];

        tr_y_yyyzzz_xxyy[i] = 2.0 * tr_y_yyyz_xxyy[i] * fe_0 + tr_y_yyyzz_xxyy[i] * pa_z[i];

        tr_y_yyyzzz_xxyz[i] = 2.0 * tr_y_yyyz_xxyz[i] * fe_0 + tr_y_yyyzz_xxy[i] * fe_0 + tr_y_yyyzz_xxyz[i] * pa_z[i];

        tr_y_yyyzzz_xxzz[i] = 2.0 * tr_y_yzzz_xxzz[i] * fe_0 + ts_yyzzz_xxzz[i] * fe_0 + tr_y_yyzzz_xxzz[i] * pa_y[i];

        tr_y_yyyzzz_xyyy[i] = 2.0 * tr_y_yyyz_xyyy[i] * fe_0 + tr_y_yyyzz_xyyy[i] * pa_z[i];

        tr_y_yyyzzz_xyyz[i] = 2.0 * tr_y_yyyz_xyyz[i] * fe_0 + tr_y_yyyzz_xyy[i] * fe_0 + tr_y_yyyzz_xyyz[i] * pa_z[i];

        tr_y_yyyzzz_xyzz[i] = 2.0 * tr_y_yyyz_xyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_xyz[i] * fe_0 + tr_y_yyyzz_xyzz[i] * pa_z[i];

        tr_y_yyyzzz_xzzz[i] = 2.0 * tr_y_yzzz_xzzz[i] * fe_0 + ts_yyzzz_xzzz[i] * fe_0 + tr_y_yyzzz_xzzz[i] * pa_y[i];

        tr_y_yyyzzz_yyyy[i] = 2.0 * tr_y_yyyz_yyyy[i] * fe_0 + tr_y_yyyzz_yyyy[i] * pa_z[i];

        tr_y_yyyzzz_yyyz[i] = 2.0 * tr_y_yyyz_yyyz[i] * fe_0 + tr_y_yyyzz_yyy[i] * fe_0 + tr_y_yyyzz_yyyz[i] * pa_z[i];

        tr_y_yyyzzz_yyzz[i] = 2.0 * tr_y_yyyz_yyzz[i] * fe_0 + 2.0 * tr_y_yyyzz_yyz[i] * fe_0 + tr_y_yyyzz_yyzz[i] * pa_z[i];

        tr_y_yyyzzz_yzzz[i] = 2.0 * tr_y_yyyz_yzzz[i] * fe_0 + 3.0 * tr_y_yyyzz_yzz[i] * fe_0 + tr_y_yyyzz_yzzz[i] * pa_z[i];

        tr_y_yyyzzz_zzzz[i] = 2.0 * tr_y_yzzz_zzzz[i] * fe_0 + ts_yyzzz_zzzz[i] * fe_0 + tr_y_yyzzz_zzzz[i] * pa_y[i];
    }

    // Set up 795-810 components of targeted buffer : IG

    auto tr_y_yyzzzz_xxxx = pbuffer.data(idx_dip_ig + 795);

    auto tr_y_yyzzzz_xxxy = pbuffer.data(idx_dip_ig + 796);

    auto tr_y_yyzzzz_xxxz = pbuffer.data(idx_dip_ig + 797);

    auto tr_y_yyzzzz_xxyy = pbuffer.data(idx_dip_ig + 798);

    auto tr_y_yyzzzz_xxyz = pbuffer.data(idx_dip_ig + 799);

    auto tr_y_yyzzzz_xxzz = pbuffer.data(idx_dip_ig + 800);

    auto tr_y_yyzzzz_xyyy = pbuffer.data(idx_dip_ig + 801);

    auto tr_y_yyzzzz_xyyz = pbuffer.data(idx_dip_ig + 802);

    auto tr_y_yyzzzz_xyzz = pbuffer.data(idx_dip_ig + 803);

    auto tr_y_yyzzzz_xzzz = pbuffer.data(idx_dip_ig + 804);

    auto tr_y_yyzzzz_yyyy = pbuffer.data(idx_dip_ig + 805);

    auto tr_y_yyzzzz_yyyz = pbuffer.data(idx_dip_ig + 806);

    auto tr_y_yyzzzz_yyzz = pbuffer.data(idx_dip_ig + 807);

    auto tr_y_yyzzzz_yzzz = pbuffer.data(idx_dip_ig + 808);

    auto tr_y_yyzzzz_zzzz = pbuffer.data(idx_dip_ig + 809);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yyzz_xxxx,   \
                             tr_y_yyzz_xxxy,   \
                             tr_y_yyzz_xxyy,   \
                             tr_y_yyzz_xxyz,   \
                             tr_y_yyzz_xyyy,   \
                             tr_y_yyzz_xyyz,   \
                             tr_y_yyzz_xyzz,   \
                             tr_y_yyzz_yyyy,   \
                             tr_y_yyzz_yyyz,   \
                             tr_y_yyzz_yyzz,   \
                             tr_y_yyzz_yzzz,   \
                             tr_y_yyzzz_xxxx,  \
                             tr_y_yyzzz_xxxy,  \
                             tr_y_yyzzz_xxy,   \
                             tr_y_yyzzz_xxyy,  \
                             tr_y_yyzzz_xxyz,  \
                             tr_y_yyzzz_xyy,   \
                             tr_y_yyzzz_xyyy,  \
                             tr_y_yyzzz_xyyz,  \
                             tr_y_yyzzz_xyz,   \
                             tr_y_yyzzz_xyzz,  \
                             tr_y_yyzzz_yyy,   \
                             tr_y_yyzzz_yyyy,  \
                             tr_y_yyzzz_yyyz,  \
                             tr_y_yyzzz_yyz,   \
                             tr_y_yyzzz_yyzz,  \
                             tr_y_yyzzz_yzz,   \
                             tr_y_yyzzz_yzzz,  \
                             tr_y_yyzzzz_xxxx, \
                             tr_y_yyzzzz_xxxy, \
                             tr_y_yyzzzz_xxxz, \
                             tr_y_yyzzzz_xxyy, \
                             tr_y_yyzzzz_xxyz, \
                             tr_y_yyzzzz_xxzz, \
                             tr_y_yyzzzz_xyyy, \
                             tr_y_yyzzzz_xyyz, \
                             tr_y_yyzzzz_xyzz, \
                             tr_y_yyzzzz_xzzz, \
                             tr_y_yyzzzz_yyyy, \
                             tr_y_yyzzzz_yyyz, \
                             tr_y_yyzzzz_yyzz, \
                             tr_y_yyzzzz_yzzz, \
                             tr_y_yyzzzz_zzzz, \
                             tr_y_yzzzz_xxxz,  \
                             tr_y_yzzzz_xxzz,  \
                             tr_y_yzzzz_xzzz,  \
                             tr_y_yzzzz_zzzz,  \
                             tr_y_zzzz_xxxz,   \
                             tr_y_zzzz_xxzz,   \
                             tr_y_zzzz_xzzz,   \
                             tr_y_zzzz_zzzz,   \
                             ts_yzzzz_xxxz,    \
                             ts_yzzzz_xxzz,    \
                             ts_yzzzz_xzzz,    \
                             ts_yzzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzzz_xxxx[i] = 3.0 * tr_y_yyzz_xxxx[i] * fe_0 + tr_y_yyzzz_xxxx[i] * pa_z[i];

        tr_y_yyzzzz_xxxy[i] = 3.0 * tr_y_yyzz_xxxy[i] * fe_0 + tr_y_yyzzz_xxxy[i] * pa_z[i];

        tr_y_yyzzzz_xxxz[i] = tr_y_zzzz_xxxz[i] * fe_0 + ts_yzzzz_xxxz[i] * fe_0 + tr_y_yzzzz_xxxz[i] * pa_y[i];

        tr_y_yyzzzz_xxyy[i] = 3.0 * tr_y_yyzz_xxyy[i] * fe_0 + tr_y_yyzzz_xxyy[i] * pa_z[i];

        tr_y_yyzzzz_xxyz[i] = 3.0 * tr_y_yyzz_xxyz[i] * fe_0 + tr_y_yyzzz_xxy[i] * fe_0 + tr_y_yyzzz_xxyz[i] * pa_z[i];

        tr_y_yyzzzz_xxzz[i] = tr_y_zzzz_xxzz[i] * fe_0 + ts_yzzzz_xxzz[i] * fe_0 + tr_y_yzzzz_xxzz[i] * pa_y[i];

        tr_y_yyzzzz_xyyy[i] = 3.0 * tr_y_yyzz_xyyy[i] * fe_0 + tr_y_yyzzz_xyyy[i] * pa_z[i];

        tr_y_yyzzzz_xyyz[i] = 3.0 * tr_y_yyzz_xyyz[i] * fe_0 + tr_y_yyzzz_xyy[i] * fe_0 + tr_y_yyzzz_xyyz[i] * pa_z[i];

        tr_y_yyzzzz_xyzz[i] = 3.0 * tr_y_yyzz_xyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_xyz[i] * fe_0 + tr_y_yyzzz_xyzz[i] * pa_z[i];

        tr_y_yyzzzz_xzzz[i] = tr_y_zzzz_xzzz[i] * fe_0 + ts_yzzzz_xzzz[i] * fe_0 + tr_y_yzzzz_xzzz[i] * pa_y[i];

        tr_y_yyzzzz_yyyy[i] = 3.0 * tr_y_yyzz_yyyy[i] * fe_0 + tr_y_yyzzz_yyyy[i] * pa_z[i];

        tr_y_yyzzzz_yyyz[i] = 3.0 * tr_y_yyzz_yyyz[i] * fe_0 + tr_y_yyzzz_yyy[i] * fe_0 + tr_y_yyzzz_yyyz[i] * pa_z[i];

        tr_y_yyzzzz_yyzz[i] = 3.0 * tr_y_yyzz_yyzz[i] * fe_0 + 2.0 * tr_y_yyzzz_yyz[i] * fe_0 + tr_y_yyzzz_yyzz[i] * pa_z[i];

        tr_y_yyzzzz_yzzz[i] = 3.0 * tr_y_yyzz_yzzz[i] * fe_0 + 3.0 * tr_y_yyzzz_yzz[i] * fe_0 + tr_y_yyzzz_yzzz[i] * pa_z[i];

        tr_y_yyzzzz_zzzz[i] = tr_y_zzzz_zzzz[i] * fe_0 + ts_yzzzz_zzzz[i] * fe_0 + tr_y_yzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 810-825 components of targeted buffer : IG

    auto tr_y_yzzzzz_xxxx = pbuffer.data(idx_dip_ig + 810);

    auto tr_y_yzzzzz_xxxy = pbuffer.data(idx_dip_ig + 811);

    auto tr_y_yzzzzz_xxxz = pbuffer.data(idx_dip_ig + 812);

    auto tr_y_yzzzzz_xxyy = pbuffer.data(idx_dip_ig + 813);

    auto tr_y_yzzzzz_xxyz = pbuffer.data(idx_dip_ig + 814);

    auto tr_y_yzzzzz_xxzz = pbuffer.data(idx_dip_ig + 815);

    auto tr_y_yzzzzz_xyyy = pbuffer.data(idx_dip_ig + 816);

    auto tr_y_yzzzzz_xyyz = pbuffer.data(idx_dip_ig + 817);

    auto tr_y_yzzzzz_xyzz = pbuffer.data(idx_dip_ig + 818);

    auto tr_y_yzzzzz_xzzz = pbuffer.data(idx_dip_ig + 819);

    auto tr_y_yzzzzz_yyyy = pbuffer.data(idx_dip_ig + 820);

    auto tr_y_yzzzzz_yyyz = pbuffer.data(idx_dip_ig + 821);

    auto tr_y_yzzzzz_yyzz = pbuffer.data(idx_dip_ig + 822);

    auto tr_y_yzzzzz_yzzz = pbuffer.data(idx_dip_ig + 823);

    auto tr_y_yzzzzz_zzzz = pbuffer.data(idx_dip_ig + 824);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yzzz_xxxy,   \
                             tr_y_yzzz_xxyy,   \
                             tr_y_yzzz_xyyy,   \
                             tr_y_yzzz_yyyy,   \
                             tr_y_yzzzz_xxxy,  \
                             tr_y_yzzzz_xxyy,  \
                             tr_y_yzzzz_xyyy,  \
                             tr_y_yzzzz_yyyy,  \
                             tr_y_yzzzzz_xxxx, \
                             tr_y_yzzzzz_xxxy, \
                             tr_y_yzzzzz_xxxz, \
                             tr_y_yzzzzz_xxyy, \
                             tr_y_yzzzzz_xxyz, \
                             tr_y_yzzzzz_xxzz, \
                             tr_y_yzzzzz_xyyy, \
                             tr_y_yzzzzz_xyyz, \
                             tr_y_yzzzzz_xyzz, \
                             tr_y_yzzzzz_xzzz, \
                             tr_y_yzzzzz_yyyy, \
                             tr_y_yzzzzz_yyyz, \
                             tr_y_yzzzzz_yyzz, \
                             tr_y_yzzzzz_yzzz, \
                             tr_y_yzzzzz_zzzz, \
                             tr_y_zzzzz_xxxx,  \
                             tr_y_zzzzz_xxxz,  \
                             tr_y_zzzzz_xxyz,  \
                             tr_y_zzzzz_xxz,   \
                             tr_y_zzzzz_xxzz,  \
                             tr_y_zzzzz_xyyz,  \
                             tr_y_zzzzz_xyz,   \
                             tr_y_zzzzz_xyzz,  \
                             tr_y_zzzzz_xzz,   \
                             tr_y_zzzzz_xzzz,  \
                             tr_y_zzzzz_yyyz,  \
                             tr_y_zzzzz_yyz,   \
                             tr_y_zzzzz_yyzz,  \
                             tr_y_zzzzz_yzz,   \
                             tr_y_zzzzz_yzzz,  \
                             tr_y_zzzzz_zzz,   \
                             tr_y_zzzzz_zzzz,  \
                             ts_zzzzz_xxxx,    \
                             ts_zzzzz_xxxz,    \
                             ts_zzzzz_xxyz,    \
                             ts_zzzzz_xxzz,    \
                             ts_zzzzz_xyyz,    \
                             ts_zzzzz_xyzz,    \
                             ts_zzzzz_xzzz,    \
                             ts_zzzzz_yyyz,    \
                             ts_zzzzz_yyzz,    \
                             ts_zzzzz_yzzz,    \
                             ts_zzzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzzz_xxxx[i] = ts_zzzzz_xxxx[i] * fe_0 + tr_y_zzzzz_xxxx[i] * pa_y[i];

        tr_y_yzzzzz_xxxy[i] = 4.0 * tr_y_yzzz_xxxy[i] * fe_0 + tr_y_yzzzz_xxxy[i] * pa_z[i];

        tr_y_yzzzzz_xxxz[i] = ts_zzzzz_xxxz[i] * fe_0 + tr_y_zzzzz_xxxz[i] * pa_y[i];

        tr_y_yzzzzz_xxyy[i] = 4.0 * tr_y_yzzz_xxyy[i] * fe_0 + tr_y_yzzzz_xxyy[i] * pa_z[i];

        tr_y_yzzzzz_xxyz[i] = tr_y_zzzzz_xxz[i] * fe_0 + ts_zzzzz_xxyz[i] * fe_0 + tr_y_zzzzz_xxyz[i] * pa_y[i];

        tr_y_yzzzzz_xxzz[i] = ts_zzzzz_xxzz[i] * fe_0 + tr_y_zzzzz_xxzz[i] * pa_y[i];

        tr_y_yzzzzz_xyyy[i] = 4.0 * tr_y_yzzz_xyyy[i] * fe_0 + tr_y_yzzzz_xyyy[i] * pa_z[i];

        tr_y_yzzzzz_xyyz[i] = 2.0 * tr_y_zzzzz_xyz[i] * fe_0 + ts_zzzzz_xyyz[i] * fe_0 + tr_y_zzzzz_xyyz[i] * pa_y[i];

        tr_y_yzzzzz_xyzz[i] = tr_y_zzzzz_xzz[i] * fe_0 + ts_zzzzz_xyzz[i] * fe_0 + tr_y_zzzzz_xyzz[i] * pa_y[i];

        tr_y_yzzzzz_xzzz[i] = ts_zzzzz_xzzz[i] * fe_0 + tr_y_zzzzz_xzzz[i] * pa_y[i];

        tr_y_yzzzzz_yyyy[i] = 4.0 * tr_y_yzzz_yyyy[i] * fe_0 + tr_y_yzzzz_yyyy[i] * pa_z[i];

        tr_y_yzzzzz_yyyz[i] = 3.0 * tr_y_zzzzz_yyz[i] * fe_0 + ts_zzzzz_yyyz[i] * fe_0 + tr_y_zzzzz_yyyz[i] * pa_y[i];

        tr_y_yzzzzz_yyzz[i] = 2.0 * tr_y_zzzzz_yzz[i] * fe_0 + ts_zzzzz_yyzz[i] * fe_0 + tr_y_zzzzz_yyzz[i] * pa_y[i];

        tr_y_yzzzzz_yzzz[i] = tr_y_zzzzz_zzz[i] * fe_0 + ts_zzzzz_yzzz[i] * fe_0 + tr_y_zzzzz_yzzz[i] * pa_y[i];

        tr_y_yzzzzz_zzzz[i] = ts_zzzzz_zzzz[i] * fe_0 + tr_y_zzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 825-840 components of targeted buffer : IG

    auto tr_y_zzzzzz_xxxx = pbuffer.data(idx_dip_ig + 825);

    auto tr_y_zzzzzz_xxxy = pbuffer.data(idx_dip_ig + 826);

    auto tr_y_zzzzzz_xxxz = pbuffer.data(idx_dip_ig + 827);

    auto tr_y_zzzzzz_xxyy = pbuffer.data(idx_dip_ig + 828);

    auto tr_y_zzzzzz_xxyz = pbuffer.data(idx_dip_ig + 829);

    auto tr_y_zzzzzz_xxzz = pbuffer.data(idx_dip_ig + 830);

    auto tr_y_zzzzzz_xyyy = pbuffer.data(idx_dip_ig + 831);

    auto tr_y_zzzzzz_xyyz = pbuffer.data(idx_dip_ig + 832);

    auto tr_y_zzzzzz_xyzz = pbuffer.data(idx_dip_ig + 833);

    auto tr_y_zzzzzz_xzzz = pbuffer.data(idx_dip_ig + 834);

    auto tr_y_zzzzzz_yyyy = pbuffer.data(idx_dip_ig + 835);

    auto tr_y_zzzzzz_yyyz = pbuffer.data(idx_dip_ig + 836);

    auto tr_y_zzzzzz_yyzz = pbuffer.data(idx_dip_ig + 837);

    auto tr_y_zzzzzz_yzzz = pbuffer.data(idx_dip_ig + 838);

    auto tr_y_zzzzzz_zzzz = pbuffer.data(idx_dip_ig + 839);

#pragma omp simd aligned(pa_z,                 \
                             tr_y_zzzz_xxxx,   \
                             tr_y_zzzz_xxxy,   \
                             tr_y_zzzz_xxxz,   \
                             tr_y_zzzz_xxyy,   \
                             tr_y_zzzz_xxyz,   \
                             tr_y_zzzz_xxzz,   \
                             tr_y_zzzz_xyyy,   \
                             tr_y_zzzz_xyyz,   \
                             tr_y_zzzz_xyzz,   \
                             tr_y_zzzz_xzzz,   \
                             tr_y_zzzz_yyyy,   \
                             tr_y_zzzz_yyyz,   \
                             tr_y_zzzz_yyzz,   \
                             tr_y_zzzz_yzzz,   \
                             tr_y_zzzz_zzzz,   \
                             tr_y_zzzzz_xxx,   \
                             tr_y_zzzzz_xxxx,  \
                             tr_y_zzzzz_xxxy,  \
                             tr_y_zzzzz_xxxz,  \
                             tr_y_zzzzz_xxy,   \
                             tr_y_zzzzz_xxyy,  \
                             tr_y_zzzzz_xxyz,  \
                             tr_y_zzzzz_xxz,   \
                             tr_y_zzzzz_xxzz,  \
                             tr_y_zzzzz_xyy,   \
                             tr_y_zzzzz_xyyy,  \
                             tr_y_zzzzz_xyyz,  \
                             tr_y_zzzzz_xyz,   \
                             tr_y_zzzzz_xyzz,  \
                             tr_y_zzzzz_xzz,   \
                             tr_y_zzzzz_xzzz,  \
                             tr_y_zzzzz_yyy,   \
                             tr_y_zzzzz_yyyy,  \
                             tr_y_zzzzz_yyyz,  \
                             tr_y_zzzzz_yyz,   \
                             tr_y_zzzzz_yyzz,  \
                             tr_y_zzzzz_yzz,   \
                             tr_y_zzzzz_yzzz,  \
                             tr_y_zzzzz_zzz,   \
                             tr_y_zzzzz_zzzz,  \
                             tr_y_zzzzzz_xxxx, \
                             tr_y_zzzzzz_xxxy, \
                             tr_y_zzzzzz_xxxz, \
                             tr_y_zzzzzz_xxyy, \
                             tr_y_zzzzzz_xxyz, \
                             tr_y_zzzzzz_xxzz, \
                             tr_y_zzzzzz_xyyy, \
                             tr_y_zzzzzz_xyyz, \
                             tr_y_zzzzzz_xyzz, \
                             tr_y_zzzzzz_xzzz, \
                             tr_y_zzzzzz_yyyy, \
                             tr_y_zzzzzz_yyyz, \
                             tr_y_zzzzzz_yyzz, \
                             tr_y_zzzzzz_yzzz, \
                             tr_y_zzzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzzz_xxxx[i] = 5.0 * tr_y_zzzz_xxxx[i] * fe_0 + tr_y_zzzzz_xxxx[i] * pa_z[i];

        tr_y_zzzzzz_xxxy[i] = 5.0 * tr_y_zzzz_xxxy[i] * fe_0 + tr_y_zzzzz_xxxy[i] * pa_z[i];

        tr_y_zzzzzz_xxxz[i] = 5.0 * tr_y_zzzz_xxxz[i] * fe_0 + tr_y_zzzzz_xxx[i] * fe_0 + tr_y_zzzzz_xxxz[i] * pa_z[i];

        tr_y_zzzzzz_xxyy[i] = 5.0 * tr_y_zzzz_xxyy[i] * fe_0 + tr_y_zzzzz_xxyy[i] * pa_z[i];

        tr_y_zzzzzz_xxyz[i] = 5.0 * tr_y_zzzz_xxyz[i] * fe_0 + tr_y_zzzzz_xxy[i] * fe_0 + tr_y_zzzzz_xxyz[i] * pa_z[i];

        tr_y_zzzzzz_xxzz[i] = 5.0 * tr_y_zzzz_xxzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xxz[i] * fe_0 + tr_y_zzzzz_xxzz[i] * pa_z[i];

        tr_y_zzzzzz_xyyy[i] = 5.0 * tr_y_zzzz_xyyy[i] * fe_0 + tr_y_zzzzz_xyyy[i] * pa_z[i];

        tr_y_zzzzzz_xyyz[i] = 5.0 * tr_y_zzzz_xyyz[i] * fe_0 + tr_y_zzzzz_xyy[i] * fe_0 + tr_y_zzzzz_xyyz[i] * pa_z[i];

        tr_y_zzzzzz_xyzz[i] = 5.0 * tr_y_zzzz_xyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xyz[i] * fe_0 + tr_y_zzzzz_xyzz[i] * pa_z[i];

        tr_y_zzzzzz_xzzz[i] = 5.0 * tr_y_zzzz_xzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_xzz[i] * fe_0 + tr_y_zzzzz_xzzz[i] * pa_z[i];

        tr_y_zzzzzz_yyyy[i] = 5.0 * tr_y_zzzz_yyyy[i] * fe_0 + tr_y_zzzzz_yyyy[i] * pa_z[i];

        tr_y_zzzzzz_yyyz[i] = 5.0 * tr_y_zzzz_yyyz[i] * fe_0 + tr_y_zzzzz_yyy[i] * fe_0 + tr_y_zzzzz_yyyz[i] * pa_z[i];

        tr_y_zzzzzz_yyzz[i] = 5.0 * tr_y_zzzz_yyzz[i] * fe_0 + 2.0 * tr_y_zzzzz_yyz[i] * fe_0 + tr_y_zzzzz_yyzz[i] * pa_z[i];

        tr_y_zzzzzz_yzzz[i] = 5.0 * tr_y_zzzz_yzzz[i] * fe_0 + 3.0 * tr_y_zzzzz_yzz[i] * fe_0 + tr_y_zzzzz_yzzz[i] * pa_z[i];

        tr_y_zzzzzz_zzzz[i] = 5.0 * tr_y_zzzz_zzzz[i] * fe_0 + 4.0 * tr_y_zzzzz_zzz[i] * fe_0 + tr_y_zzzzz_zzzz[i] * pa_z[i];
    }

    // Set up 840-855 components of targeted buffer : IG

    auto tr_z_xxxxxx_xxxx = pbuffer.data(idx_dip_ig + 840);

    auto tr_z_xxxxxx_xxxy = pbuffer.data(idx_dip_ig + 841);

    auto tr_z_xxxxxx_xxxz = pbuffer.data(idx_dip_ig + 842);

    auto tr_z_xxxxxx_xxyy = pbuffer.data(idx_dip_ig + 843);

    auto tr_z_xxxxxx_xxyz = pbuffer.data(idx_dip_ig + 844);

    auto tr_z_xxxxxx_xxzz = pbuffer.data(idx_dip_ig + 845);

    auto tr_z_xxxxxx_xyyy = pbuffer.data(idx_dip_ig + 846);

    auto tr_z_xxxxxx_xyyz = pbuffer.data(idx_dip_ig + 847);

    auto tr_z_xxxxxx_xyzz = pbuffer.data(idx_dip_ig + 848);

    auto tr_z_xxxxxx_xzzz = pbuffer.data(idx_dip_ig + 849);

    auto tr_z_xxxxxx_yyyy = pbuffer.data(idx_dip_ig + 850);

    auto tr_z_xxxxxx_yyyz = pbuffer.data(idx_dip_ig + 851);

    auto tr_z_xxxxxx_yyzz = pbuffer.data(idx_dip_ig + 852);

    auto tr_z_xxxxxx_yzzz = pbuffer.data(idx_dip_ig + 853);

    auto tr_z_xxxxxx_zzzz = pbuffer.data(idx_dip_ig + 854);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxxx_xxxx,   \
                             tr_z_xxxx_xxxy,   \
                             tr_z_xxxx_xxxz,   \
                             tr_z_xxxx_xxyy,   \
                             tr_z_xxxx_xxyz,   \
                             tr_z_xxxx_xxzz,   \
                             tr_z_xxxx_xyyy,   \
                             tr_z_xxxx_xyyz,   \
                             tr_z_xxxx_xyzz,   \
                             tr_z_xxxx_xzzz,   \
                             tr_z_xxxx_yyyy,   \
                             tr_z_xxxx_yyyz,   \
                             tr_z_xxxx_yyzz,   \
                             tr_z_xxxx_yzzz,   \
                             tr_z_xxxx_zzzz,   \
                             tr_z_xxxxx_xxx,   \
                             tr_z_xxxxx_xxxx,  \
                             tr_z_xxxxx_xxxy,  \
                             tr_z_xxxxx_xxxz,  \
                             tr_z_xxxxx_xxy,   \
                             tr_z_xxxxx_xxyy,  \
                             tr_z_xxxxx_xxyz,  \
                             tr_z_xxxxx_xxz,   \
                             tr_z_xxxxx_xxzz,  \
                             tr_z_xxxxx_xyy,   \
                             tr_z_xxxxx_xyyy,  \
                             tr_z_xxxxx_xyyz,  \
                             tr_z_xxxxx_xyz,   \
                             tr_z_xxxxx_xyzz,  \
                             tr_z_xxxxx_xzz,   \
                             tr_z_xxxxx_xzzz,  \
                             tr_z_xxxxx_yyy,   \
                             tr_z_xxxxx_yyyy,  \
                             tr_z_xxxxx_yyyz,  \
                             tr_z_xxxxx_yyz,   \
                             tr_z_xxxxx_yyzz,  \
                             tr_z_xxxxx_yzz,   \
                             tr_z_xxxxx_yzzz,  \
                             tr_z_xxxxx_zzz,   \
                             tr_z_xxxxx_zzzz,  \
                             tr_z_xxxxxx_xxxx, \
                             tr_z_xxxxxx_xxxy, \
                             tr_z_xxxxxx_xxxz, \
                             tr_z_xxxxxx_xxyy, \
                             tr_z_xxxxxx_xxyz, \
                             tr_z_xxxxxx_xxzz, \
                             tr_z_xxxxxx_xyyy, \
                             tr_z_xxxxxx_xyyz, \
                             tr_z_xxxxxx_xyzz, \
                             tr_z_xxxxxx_xzzz, \
                             tr_z_xxxxxx_yyyy, \
                             tr_z_xxxxxx_yyyz, \
                             tr_z_xxxxxx_yyzz, \
                             tr_z_xxxxxx_yzzz, \
                             tr_z_xxxxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxx_xxxx[i] = 5.0 * tr_z_xxxx_xxxx[i] * fe_0 + 4.0 * tr_z_xxxxx_xxx[i] * fe_0 + tr_z_xxxxx_xxxx[i] * pa_x[i];

        tr_z_xxxxxx_xxxy[i] = 5.0 * tr_z_xxxx_xxxy[i] * fe_0 + 3.0 * tr_z_xxxxx_xxy[i] * fe_0 + tr_z_xxxxx_xxxy[i] * pa_x[i];

        tr_z_xxxxxx_xxxz[i] = 5.0 * tr_z_xxxx_xxxz[i] * fe_0 + 3.0 * tr_z_xxxxx_xxz[i] * fe_0 + tr_z_xxxxx_xxxz[i] * pa_x[i];

        tr_z_xxxxxx_xxyy[i] = 5.0 * tr_z_xxxx_xxyy[i] * fe_0 + 2.0 * tr_z_xxxxx_xyy[i] * fe_0 + tr_z_xxxxx_xxyy[i] * pa_x[i];

        tr_z_xxxxxx_xxyz[i] = 5.0 * tr_z_xxxx_xxyz[i] * fe_0 + 2.0 * tr_z_xxxxx_xyz[i] * fe_0 + tr_z_xxxxx_xxyz[i] * pa_x[i];

        tr_z_xxxxxx_xxzz[i] = 5.0 * tr_z_xxxx_xxzz[i] * fe_0 + 2.0 * tr_z_xxxxx_xzz[i] * fe_0 + tr_z_xxxxx_xxzz[i] * pa_x[i];

        tr_z_xxxxxx_xyyy[i] = 5.0 * tr_z_xxxx_xyyy[i] * fe_0 + tr_z_xxxxx_yyy[i] * fe_0 + tr_z_xxxxx_xyyy[i] * pa_x[i];

        tr_z_xxxxxx_xyyz[i] = 5.0 * tr_z_xxxx_xyyz[i] * fe_0 + tr_z_xxxxx_yyz[i] * fe_0 + tr_z_xxxxx_xyyz[i] * pa_x[i];

        tr_z_xxxxxx_xyzz[i] = 5.0 * tr_z_xxxx_xyzz[i] * fe_0 + tr_z_xxxxx_yzz[i] * fe_0 + tr_z_xxxxx_xyzz[i] * pa_x[i];

        tr_z_xxxxxx_xzzz[i] = 5.0 * tr_z_xxxx_xzzz[i] * fe_0 + tr_z_xxxxx_zzz[i] * fe_0 + tr_z_xxxxx_xzzz[i] * pa_x[i];

        tr_z_xxxxxx_yyyy[i] = 5.0 * tr_z_xxxx_yyyy[i] * fe_0 + tr_z_xxxxx_yyyy[i] * pa_x[i];

        tr_z_xxxxxx_yyyz[i] = 5.0 * tr_z_xxxx_yyyz[i] * fe_0 + tr_z_xxxxx_yyyz[i] * pa_x[i];

        tr_z_xxxxxx_yyzz[i] = 5.0 * tr_z_xxxx_yyzz[i] * fe_0 + tr_z_xxxxx_yyzz[i] * pa_x[i];

        tr_z_xxxxxx_yzzz[i] = 5.0 * tr_z_xxxx_yzzz[i] * fe_0 + tr_z_xxxxx_yzzz[i] * pa_x[i];

        tr_z_xxxxxx_zzzz[i] = 5.0 * tr_z_xxxx_zzzz[i] * fe_0 + tr_z_xxxxx_zzzz[i] * pa_x[i];
    }

    // Set up 855-870 components of targeted buffer : IG

    auto tr_z_xxxxxy_xxxx = pbuffer.data(idx_dip_ig + 855);

    auto tr_z_xxxxxy_xxxy = pbuffer.data(idx_dip_ig + 856);

    auto tr_z_xxxxxy_xxxz = pbuffer.data(idx_dip_ig + 857);

    auto tr_z_xxxxxy_xxyy = pbuffer.data(idx_dip_ig + 858);

    auto tr_z_xxxxxy_xxyz = pbuffer.data(idx_dip_ig + 859);

    auto tr_z_xxxxxy_xxzz = pbuffer.data(idx_dip_ig + 860);

    auto tr_z_xxxxxy_xyyy = pbuffer.data(idx_dip_ig + 861);

    auto tr_z_xxxxxy_xyyz = pbuffer.data(idx_dip_ig + 862);

    auto tr_z_xxxxxy_xyzz = pbuffer.data(idx_dip_ig + 863);

    auto tr_z_xxxxxy_xzzz = pbuffer.data(idx_dip_ig + 864);

    auto tr_z_xxxxxy_yyyy = pbuffer.data(idx_dip_ig + 865);

    auto tr_z_xxxxxy_yyyz = pbuffer.data(idx_dip_ig + 866);

    auto tr_z_xxxxxy_yyzz = pbuffer.data(idx_dip_ig + 867);

    auto tr_z_xxxxxy_yzzz = pbuffer.data(idx_dip_ig + 868);

    auto tr_z_xxxxxy_zzzz = pbuffer.data(idx_dip_ig + 869);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxxx_xxx,   \
                             tr_z_xxxxx_xxxx,  \
                             tr_z_xxxxx_xxxy,  \
                             tr_z_xxxxx_xxxz,  \
                             tr_z_xxxxx_xxy,   \
                             tr_z_xxxxx_xxyy,  \
                             tr_z_xxxxx_xxyz,  \
                             tr_z_xxxxx_xxz,   \
                             tr_z_xxxxx_xxzz,  \
                             tr_z_xxxxx_xyy,   \
                             tr_z_xxxxx_xyyy,  \
                             tr_z_xxxxx_xyyz,  \
                             tr_z_xxxxx_xyz,   \
                             tr_z_xxxxx_xyzz,  \
                             tr_z_xxxxx_xzz,   \
                             tr_z_xxxxx_xzzz,  \
                             tr_z_xxxxx_zzzz,  \
                             tr_z_xxxxxy_xxxx, \
                             tr_z_xxxxxy_xxxy, \
                             tr_z_xxxxxy_xxxz, \
                             tr_z_xxxxxy_xxyy, \
                             tr_z_xxxxxy_xxyz, \
                             tr_z_xxxxxy_xxzz, \
                             tr_z_xxxxxy_xyyy, \
                             tr_z_xxxxxy_xyyz, \
                             tr_z_xxxxxy_xyzz, \
                             tr_z_xxxxxy_xzzz, \
                             tr_z_xxxxxy_yyyy, \
                             tr_z_xxxxxy_yyyz, \
                             tr_z_xxxxxy_yyzz, \
                             tr_z_xxxxxy_yzzz, \
                             tr_z_xxxxxy_zzzz, \
                             tr_z_xxxxy_yyyy,  \
                             tr_z_xxxxy_yyyz,  \
                             tr_z_xxxxy_yyzz,  \
                             tr_z_xxxxy_yzzz,  \
                             tr_z_xxxy_yyyy,   \
                             tr_z_xxxy_yyyz,   \
                             tr_z_xxxy_yyzz,   \
                             tr_z_xxxy_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxy_xxxx[i] = tr_z_xxxxx_xxxx[i] * pa_y[i];

        tr_z_xxxxxy_xxxy[i] = tr_z_xxxxx_xxx[i] * fe_0 + tr_z_xxxxx_xxxy[i] * pa_y[i];

        tr_z_xxxxxy_xxxz[i] = tr_z_xxxxx_xxxz[i] * pa_y[i];

        tr_z_xxxxxy_xxyy[i] = 2.0 * tr_z_xxxxx_xxy[i] * fe_0 + tr_z_xxxxx_xxyy[i] * pa_y[i];

        tr_z_xxxxxy_xxyz[i] = tr_z_xxxxx_xxz[i] * fe_0 + tr_z_xxxxx_xxyz[i] * pa_y[i];

        tr_z_xxxxxy_xxzz[i] = tr_z_xxxxx_xxzz[i] * pa_y[i];

        tr_z_xxxxxy_xyyy[i] = 3.0 * tr_z_xxxxx_xyy[i] * fe_0 + tr_z_xxxxx_xyyy[i] * pa_y[i];

        tr_z_xxxxxy_xyyz[i] = 2.0 * tr_z_xxxxx_xyz[i] * fe_0 + tr_z_xxxxx_xyyz[i] * pa_y[i];

        tr_z_xxxxxy_xyzz[i] = tr_z_xxxxx_xzz[i] * fe_0 + tr_z_xxxxx_xyzz[i] * pa_y[i];

        tr_z_xxxxxy_xzzz[i] = tr_z_xxxxx_xzzz[i] * pa_y[i];

        tr_z_xxxxxy_yyyy[i] = 4.0 * tr_z_xxxy_yyyy[i] * fe_0 + tr_z_xxxxy_yyyy[i] * pa_x[i];

        tr_z_xxxxxy_yyyz[i] = 4.0 * tr_z_xxxy_yyyz[i] * fe_0 + tr_z_xxxxy_yyyz[i] * pa_x[i];

        tr_z_xxxxxy_yyzz[i] = 4.0 * tr_z_xxxy_yyzz[i] * fe_0 + tr_z_xxxxy_yyzz[i] * pa_x[i];

        tr_z_xxxxxy_yzzz[i] = 4.0 * tr_z_xxxy_yzzz[i] * fe_0 + tr_z_xxxxy_yzzz[i] * pa_x[i];

        tr_z_xxxxxy_zzzz[i] = tr_z_xxxxx_zzzz[i] * pa_y[i];
    }

    // Set up 870-885 components of targeted buffer : IG

    auto tr_z_xxxxxz_xxxx = pbuffer.data(idx_dip_ig + 870);

    auto tr_z_xxxxxz_xxxy = pbuffer.data(idx_dip_ig + 871);

    auto tr_z_xxxxxz_xxxz = pbuffer.data(idx_dip_ig + 872);

    auto tr_z_xxxxxz_xxyy = pbuffer.data(idx_dip_ig + 873);

    auto tr_z_xxxxxz_xxyz = pbuffer.data(idx_dip_ig + 874);

    auto tr_z_xxxxxz_xxzz = pbuffer.data(idx_dip_ig + 875);

    auto tr_z_xxxxxz_xyyy = pbuffer.data(idx_dip_ig + 876);

    auto tr_z_xxxxxz_xyyz = pbuffer.data(idx_dip_ig + 877);

    auto tr_z_xxxxxz_xyzz = pbuffer.data(idx_dip_ig + 878);

    auto tr_z_xxxxxz_xzzz = pbuffer.data(idx_dip_ig + 879);

    auto tr_z_xxxxxz_yyyy = pbuffer.data(idx_dip_ig + 880);

    auto tr_z_xxxxxz_yyyz = pbuffer.data(idx_dip_ig + 881);

    auto tr_z_xxxxxz_yyzz = pbuffer.data(idx_dip_ig + 882);

    auto tr_z_xxxxxz_yzzz = pbuffer.data(idx_dip_ig + 883);

    auto tr_z_xxxxxz_zzzz = pbuffer.data(idx_dip_ig + 884);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_z_xxxxx_xxxx,  \
                             tr_z_xxxxx_xxxy,  \
                             tr_z_xxxxx_xxyy,  \
                             tr_z_xxxxx_xyyy,  \
                             tr_z_xxxxxz_xxxx, \
                             tr_z_xxxxxz_xxxy, \
                             tr_z_xxxxxz_xxxz, \
                             tr_z_xxxxxz_xxyy, \
                             tr_z_xxxxxz_xxyz, \
                             tr_z_xxxxxz_xxzz, \
                             tr_z_xxxxxz_xyyy, \
                             tr_z_xxxxxz_xyyz, \
                             tr_z_xxxxxz_xyzz, \
                             tr_z_xxxxxz_xzzz, \
                             tr_z_xxxxxz_yyyy, \
                             tr_z_xxxxxz_yyyz, \
                             tr_z_xxxxxz_yyzz, \
                             tr_z_xxxxxz_yzzz, \
                             tr_z_xxxxxz_zzzz, \
                             tr_z_xxxxz_xxxz,  \
                             tr_z_xxxxz_xxyz,  \
                             tr_z_xxxxz_xxz,   \
                             tr_z_xxxxz_xxzz,  \
                             tr_z_xxxxz_xyyz,  \
                             tr_z_xxxxz_xyz,   \
                             tr_z_xxxxz_xyzz,  \
                             tr_z_xxxxz_xzz,   \
                             tr_z_xxxxz_xzzz,  \
                             tr_z_xxxxz_yyyy,  \
                             tr_z_xxxxz_yyyz,  \
                             tr_z_xxxxz_yyz,   \
                             tr_z_xxxxz_yyzz,  \
                             tr_z_xxxxz_yzz,   \
                             tr_z_xxxxz_yzzz,  \
                             tr_z_xxxxz_zzz,   \
                             tr_z_xxxxz_zzzz,  \
                             tr_z_xxxz_xxxz,   \
                             tr_z_xxxz_xxyz,   \
                             tr_z_xxxz_xxzz,   \
                             tr_z_xxxz_xyyz,   \
                             tr_z_xxxz_xyzz,   \
                             tr_z_xxxz_xzzz,   \
                             tr_z_xxxz_yyyy,   \
                             tr_z_xxxz_yyyz,   \
                             tr_z_xxxz_yyzz,   \
                             tr_z_xxxz_yzzz,   \
                             tr_z_xxxz_zzzz,   \
                             ts_xxxxx_xxxx,    \
                             ts_xxxxx_xxxy,    \
                             ts_xxxxx_xxyy,    \
                             ts_xxxxx_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxz_xxxx[i] = ts_xxxxx_xxxx[i] * fe_0 + tr_z_xxxxx_xxxx[i] * pa_z[i];

        tr_z_xxxxxz_xxxy[i] = ts_xxxxx_xxxy[i] * fe_0 + tr_z_xxxxx_xxxy[i] * pa_z[i];

        tr_z_xxxxxz_xxxz[i] = 4.0 * tr_z_xxxz_xxxz[i] * fe_0 + 3.0 * tr_z_xxxxz_xxz[i] * fe_0 + tr_z_xxxxz_xxxz[i] * pa_x[i];

        tr_z_xxxxxz_xxyy[i] = ts_xxxxx_xxyy[i] * fe_0 + tr_z_xxxxx_xxyy[i] * pa_z[i];

        tr_z_xxxxxz_xxyz[i] = 4.0 * tr_z_xxxz_xxyz[i] * fe_0 + 2.0 * tr_z_xxxxz_xyz[i] * fe_0 + tr_z_xxxxz_xxyz[i] * pa_x[i];

        tr_z_xxxxxz_xxzz[i] = 4.0 * tr_z_xxxz_xxzz[i] * fe_0 + 2.0 * tr_z_xxxxz_xzz[i] * fe_0 + tr_z_xxxxz_xxzz[i] * pa_x[i];

        tr_z_xxxxxz_xyyy[i] = ts_xxxxx_xyyy[i] * fe_0 + tr_z_xxxxx_xyyy[i] * pa_z[i];

        tr_z_xxxxxz_xyyz[i] = 4.0 * tr_z_xxxz_xyyz[i] * fe_0 + tr_z_xxxxz_yyz[i] * fe_0 + tr_z_xxxxz_xyyz[i] * pa_x[i];

        tr_z_xxxxxz_xyzz[i] = 4.0 * tr_z_xxxz_xyzz[i] * fe_0 + tr_z_xxxxz_yzz[i] * fe_0 + tr_z_xxxxz_xyzz[i] * pa_x[i];

        tr_z_xxxxxz_xzzz[i] = 4.0 * tr_z_xxxz_xzzz[i] * fe_0 + tr_z_xxxxz_zzz[i] * fe_0 + tr_z_xxxxz_xzzz[i] * pa_x[i];

        tr_z_xxxxxz_yyyy[i] = 4.0 * tr_z_xxxz_yyyy[i] * fe_0 + tr_z_xxxxz_yyyy[i] * pa_x[i];

        tr_z_xxxxxz_yyyz[i] = 4.0 * tr_z_xxxz_yyyz[i] * fe_0 + tr_z_xxxxz_yyyz[i] * pa_x[i];

        tr_z_xxxxxz_yyzz[i] = 4.0 * tr_z_xxxz_yyzz[i] * fe_0 + tr_z_xxxxz_yyzz[i] * pa_x[i];

        tr_z_xxxxxz_yzzz[i] = 4.0 * tr_z_xxxz_yzzz[i] * fe_0 + tr_z_xxxxz_yzzz[i] * pa_x[i];

        tr_z_xxxxxz_zzzz[i] = 4.0 * tr_z_xxxz_zzzz[i] * fe_0 + tr_z_xxxxz_zzzz[i] * pa_x[i];
    }

    // Set up 885-900 components of targeted buffer : IG

    auto tr_z_xxxxyy_xxxx = pbuffer.data(idx_dip_ig + 885);

    auto tr_z_xxxxyy_xxxy = pbuffer.data(idx_dip_ig + 886);

    auto tr_z_xxxxyy_xxxz = pbuffer.data(idx_dip_ig + 887);

    auto tr_z_xxxxyy_xxyy = pbuffer.data(idx_dip_ig + 888);

    auto tr_z_xxxxyy_xxyz = pbuffer.data(idx_dip_ig + 889);

    auto tr_z_xxxxyy_xxzz = pbuffer.data(idx_dip_ig + 890);

    auto tr_z_xxxxyy_xyyy = pbuffer.data(idx_dip_ig + 891);

    auto tr_z_xxxxyy_xyyz = pbuffer.data(idx_dip_ig + 892);

    auto tr_z_xxxxyy_xyzz = pbuffer.data(idx_dip_ig + 893);

    auto tr_z_xxxxyy_xzzz = pbuffer.data(idx_dip_ig + 894);

    auto tr_z_xxxxyy_yyyy = pbuffer.data(idx_dip_ig + 895);

    auto tr_z_xxxxyy_yyyz = pbuffer.data(idx_dip_ig + 896);

    auto tr_z_xxxxyy_yyzz = pbuffer.data(idx_dip_ig + 897);

    auto tr_z_xxxxyy_yzzz = pbuffer.data(idx_dip_ig + 898);

    auto tr_z_xxxxyy_zzzz = pbuffer.data(idx_dip_ig + 899);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxx_xxxx,   \
                             tr_z_xxxx_xxxz,   \
                             tr_z_xxxx_xxzz,   \
                             tr_z_xxxx_xzzz,   \
                             tr_z_xxxxy_xxxx,  \
                             tr_z_xxxxy_xxxz,  \
                             tr_z_xxxxy_xxzz,  \
                             tr_z_xxxxy_xzzz,  \
                             tr_z_xxxxyy_xxxx, \
                             tr_z_xxxxyy_xxxy, \
                             tr_z_xxxxyy_xxxz, \
                             tr_z_xxxxyy_xxyy, \
                             tr_z_xxxxyy_xxyz, \
                             tr_z_xxxxyy_xxzz, \
                             tr_z_xxxxyy_xyyy, \
                             tr_z_xxxxyy_xyyz, \
                             tr_z_xxxxyy_xyzz, \
                             tr_z_xxxxyy_xzzz, \
                             tr_z_xxxxyy_yyyy, \
                             tr_z_xxxxyy_yyyz, \
                             tr_z_xxxxyy_yyzz, \
                             tr_z_xxxxyy_yzzz, \
                             tr_z_xxxxyy_zzzz, \
                             tr_z_xxxyy_xxxy,  \
                             tr_z_xxxyy_xxy,   \
                             tr_z_xxxyy_xxyy,  \
                             tr_z_xxxyy_xxyz,  \
                             tr_z_xxxyy_xyy,   \
                             tr_z_xxxyy_xyyy,  \
                             tr_z_xxxyy_xyyz,  \
                             tr_z_xxxyy_xyz,   \
                             tr_z_xxxyy_xyzz,  \
                             tr_z_xxxyy_yyy,   \
                             tr_z_xxxyy_yyyy,  \
                             tr_z_xxxyy_yyyz,  \
                             tr_z_xxxyy_yyz,   \
                             tr_z_xxxyy_yyzz,  \
                             tr_z_xxxyy_yzz,   \
                             tr_z_xxxyy_yzzz,  \
                             tr_z_xxxyy_zzzz,  \
                             tr_z_xxyy_xxxy,   \
                             tr_z_xxyy_xxyy,   \
                             tr_z_xxyy_xxyz,   \
                             tr_z_xxyy_xyyy,   \
                             tr_z_xxyy_xyyz,   \
                             tr_z_xxyy_xyzz,   \
                             tr_z_xxyy_yyyy,   \
                             tr_z_xxyy_yyyz,   \
                             tr_z_xxyy_yyzz,   \
                             tr_z_xxyy_yzzz,   \
                             tr_z_xxyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyy_xxxx[i] = tr_z_xxxx_xxxx[i] * fe_0 + tr_z_xxxxy_xxxx[i] * pa_y[i];

        tr_z_xxxxyy_xxxy[i] = 3.0 * tr_z_xxyy_xxxy[i] * fe_0 + 3.0 * tr_z_xxxyy_xxy[i] * fe_0 + tr_z_xxxyy_xxxy[i] * pa_x[i];

        tr_z_xxxxyy_xxxz[i] = tr_z_xxxx_xxxz[i] * fe_0 + tr_z_xxxxy_xxxz[i] * pa_y[i];

        tr_z_xxxxyy_xxyy[i] = 3.0 * tr_z_xxyy_xxyy[i] * fe_0 + 2.0 * tr_z_xxxyy_xyy[i] * fe_0 + tr_z_xxxyy_xxyy[i] * pa_x[i];

        tr_z_xxxxyy_xxyz[i] = 3.0 * tr_z_xxyy_xxyz[i] * fe_0 + 2.0 * tr_z_xxxyy_xyz[i] * fe_0 + tr_z_xxxyy_xxyz[i] * pa_x[i];

        tr_z_xxxxyy_xxzz[i] = tr_z_xxxx_xxzz[i] * fe_0 + tr_z_xxxxy_xxzz[i] * pa_y[i];

        tr_z_xxxxyy_xyyy[i] = 3.0 * tr_z_xxyy_xyyy[i] * fe_0 + tr_z_xxxyy_yyy[i] * fe_0 + tr_z_xxxyy_xyyy[i] * pa_x[i];

        tr_z_xxxxyy_xyyz[i] = 3.0 * tr_z_xxyy_xyyz[i] * fe_0 + tr_z_xxxyy_yyz[i] * fe_0 + tr_z_xxxyy_xyyz[i] * pa_x[i];

        tr_z_xxxxyy_xyzz[i] = 3.0 * tr_z_xxyy_xyzz[i] * fe_0 + tr_z_xxxyy_yzz[i] * fe_0 + tr_z_xxxyy_xyzz[i] * pa_x[i];

        tr_z_xxxxyy_xzzz[i] = tr_z_xxxx_xzzz[i] * fe_0 + tr_z_xxxxy_xzzz[i] * pa_y[i];

        tr_z_xxxxyy_yyyy[i] = 3.0 * tr_z_xxyy_yyyy[i] * fe_0 + tr_z_xxxyy_yyyy[i] * pa_x[i];

        tr_z_xxxxyy_yyyz[i] = 3.0 * tr_z_xxyy_yyyz[i] * fe_0 + tr_z_xxxyy_yyyz[i] * pa_x[i];

        tr_z_xxxxyy_yyzz[i] = 3.0 * tr_z_xxyy_yyzz[i] * fe_0 + tr_z_xxxyy_yyzz[i] * pa_x[i];

        tr_z_xxxxyy_yzzz[i] = 3.0 * tr_z_xxyy_yzzz[i] * fe_0 + tr_z_xxxyy_yzzz[i] * pa_x[i];

        tr_z_xxxxyy_zzzz[i] = 3.0 * tr_z_xxyy_zzzz[i] * fe_0 + tr_z_xxxyy_zzzz[i] * pa_x[i];
    }

    // Set up 900-915 components of targeted buffer : IG

    auto tr_z_xxxxyz_xxxx = pbuffer.data(idx_dip_ig + 900);

    auto tr_z_xxxxyz_xxxy = pbuffer.data(idx_dip_ig + 901);

    auto tr_z_xxxxyz_xxxz = pbuffer.data(idx_dip_ig + 902);

    auto tr_z_xxxxyz_xxyy = pbuffer.data(idx_dip_ig + 903);

    auto tr_z_xxxxyz_xxyz = pbuffer.data(idx_dip_ig + 904);

    auto tr_z_xxxxyz_xxzz = pbuffer.data(idx_dip_ig + 905);

    auto tr_z_xxxxyz_xyyy = pbuffer.data(idx_dip_ig + 906);

    auto tr_z_xxxxyz_xyyz = pbuffer.data(idx_dip_ig + 907);

    auto tr_z_xxxxyz_xyzz = pbuffer.data(idx_dip_ig + 908);

    auto tr_z_xxxxyz_xzzz = pbuffer.data(idx_dip_ig + 909);

    auto tr_z_xxxxyz_yyyy = pbuffer.data(idx_dip_ig + 910);

    auto tr_z_xxxxyz_yyyz = pbuffer.data(idx_dip_ig + 911);

    auto tr_z_xxxxyz_yyzz = pbuffer.data(idx_dip_ig + 912);

    auto tr_z_xxxxyz_yzzz = pbuffer.data(idx_dip_ig + 913);

    auto tr_z_xxxxyz_zzzz = pbuffer.data(idx_dip_ig + 914);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxxyz_xxxx, \
                             tr_z_xxxxyz_xxxy, \
                             tr_z_xxxxyz_xxxz, \
                             tr_z_xxxxyz_xxyy, \
                             tr_z_xxxxyz_xxyz, \
                             tr_z_xxxxyz_xxzz, \
                             tr_z_xxxxyz_xyyy, \
                             tr_z_xxxxyz_xyyz, \
                             tr_z_xxxxyz_xyzz, \
                             tr_z_xxxxyz_xzzz, \
                             tr_z_xxxxyz_yyyy, \
                             tr_z_xxxxyz_yyyz, \
                             tr_z_xxxxyz_yyzz, \
                             tr_z_xxxxyz_yzzz, \
                             tr_z_xxxxyz_zzzz, \
                             tr_z_xxxxz_xxx,   \
                             tr_z_xxxxz_xxxx,  \
                             tr_z_xxxxz_xxxy,  \
                             tr_z_xxxxz_xxxz,  \
                             tr_z_xxxxz_xxy,   \
                             tr_z_xxxxz_xxyy,  \
                             tr_z_xxxxz_xxyz,  \
                             tr_z_xxxxz_xxz,   \
                             tr_z_xxxxz_xxzz,  \
                             tr_z_xxxxz_xyy,   \
                             tr_z_xxxxz_xyyy,  \
                             tr_z_xxxxz_xyyz,  \
                             tr_z_xxxxz_xyz,   \
                             tr_z_xxxxz_xyzz,  \
                             tr_z_xxxxz_xzz,   \
                             tr_z_xxxxz_xzzz,  \
                             tr_z_xxxxz_zzzz,  \
                             tr_z_xxxyz_yyyy,  \
                             tr_z_xxxyz_yyyz,  \
                             tr_z_xxxyz_yyzz,  \
                             tr_z_xxxyz_yzzz,  \
                             tr_z_xxyz_yyyy,   \
                             tr_z_xxyz_yyyz,   \
                             tr_z_xxyz_yyzz,   \
                             tr_z_xxyz_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyz_xxxx[i] = tr_z_xxxxz_xxxx[i] * pa_y[i];

        tr_z_xxxxyz_xxxy[i] = tr_z_xxxxz_xxx[i] * fe_0 + tr_z_xxxxz_xxxy[i] * pa_y[i];

        tr_z_xxxxyz_xxxz[i] = tr_z_xxxxz_xxxz[i] * pa_y[i];

        tr_z_xxxxyz_xxyy[i] = 2.0 * tr_z_xxxxz_xxy[i] * fe_0 + tr_z_xxxxz_xxyy[i] * pa_y[i];

        tr_z_xxxxyz_xxyz[i] = tr_z_xxxxz_xxz[i] * fe_0 + tr_z_xxxxz_xxyz[i] * pa_y[i];

        tr_z_xxxxyz_xxzz[i] = tr_z_xxxxz_xxzz[i] * pa_y[i];

        tr_z_xxxxyz_xyyy[i] = 3.0 * tr_z_xxxxz_xyy[i] * fe_0 + tr_z_xxxxz_xyyy[i] * pa_y[i];

        tr_z_xxxxyz_xyyz[i] = 2.0 * tr_z_xxxxz_xyz[i] * fe_0 + tr_z_xxxxz_xyyz[i] * pa_y[i];

        tr_z_xxxxyz_xyzz[i] = tr_z_xxxxz_xzz[i] * fe_0 + tr_z_xxxxz_xyzz[i] * pa_y[i];

        tr_z_xxxxyz_xzzz[i] = tr_z_xxxxz_xzzz[i] * pa_y[i];

        tr_z_xxxxyz_yyyy[i] = 3.0 * tr_z_xxyz_yyyy[i] * fe_0 + tr_z_xxxyz_yyyy[i] * pa_x[i];

        tr_z_xxxxyz_yyyz[i] = 3.0 * tr_z_xxyz_yyyz[i] * fe_0 + tr_z_xxxyz_yyyz[i] * pa_x[i];

        tr_z_xxxxyz_yyzz[i] = 3.0 * tr_z_xxyz_yyzz[i] * fe_0 + tr_z_xxxyz_yyzz[i] * pa_x[i];

        tr_z_xxxxyz_yzzz[i] = 3.0 * tr_z_xxyz_yzzz[i] * fe_0 + tr_z_xxxyz_yzzz[i] * pa_x[i];

        tr_z_xxxxyz_zzzz[i] = tr_z_xxxxz_zzzz[i] * pa_y[i];
    }

    // Set up 915-930 components of targeted buffer : IG

    auto tr_z_xxxxzz_xxxx = pbuffer.data(idx_dip_ig + 915);

    auto tr_z_xxxxzz_xxxy = pbuffer.data(idx_dip_ig + 916);

    auto tr_z_xxxxzz_xxxz = pbuffer.data(idx_dip_ig + 917);

    auto tr_z_xxxxzz_xxyy = pbuffer.data(idx_dip_ig + 918);

    auto tr_z_xxxxzz_xxyz = pbuffer.data(idx_dip_ig + 919);

    auto tr_z_xxxxzz_xxzz = pbuffer.data(idx_dip_ig + 920);

    auto tr_z_xxxxzz_xyyy = pbuffer.data(idx_dip_ig + 921);

    auto tr_z_xxxxzz_xyyz = pbuffer.data(idx_dip_ig + 922);

    auto tr_z_xxxxzz_xyzz = pbuffer.data(idx_dip_ig + 923);

    auto tr_z_xxxxzz_xzzz = pbuffer.data(idx_dip_ig + 924);

    auto tr_z_xxxxzz_yyyy = pbuffer.data(idx_dip_ig + 925);

    auto tr_z_xxxxzz_yyyz = pbuffer.data(idx_dip_ig + 926);

    auto tr_z_xxxxzz_yyzz = pbuffer.data(idx_dip_ig + 927);

    auto tr_z_xxxxzz_yzzz = pbuffer.data(idx_dip_ig + 928);

    auto tr_z_xxxxzz_zzzz = pbuffer.data(idx_dip_ig + 929);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxxxzz_xxxx, \
                             tr_z_xxxxzz_xxxy, \
                             tr_z_xxxxzz_xxxz, \
                             tr_z_xxxxzz_xxyy, \
                             tr_z_xxxxzz_xxyz, \
                             tr_z_xxxxzz_xxzz, \
                             tr_z_xxxxzz_xyyy, \
                             tr_z_xxxxzz_xyyz, \
                             tr_z_xxxxzz_xyzz, \
                             tr_z_xxxxzz_xzzz, \
                             tr_z_xxxxzz_yyyy, \
                             tr_z_xxxxzz_yyyz, \
                             tr_z_xxxxzz_yyzz, \
                             tr_z_xxxxzz_yzzz, \
                             tr_z_xxxxzz_zzzz, \
                             tr_z_xxxzz_xxx,   \
                             tr_z_xxxzz_xxxx,  \
                             tr_z_xxxzz_xxxy,  \
                             tr_z_xxxzz_xxxz,  \
                             tr_z_xxxzz_xxy,   \
                             tr_z_xxxzz_xxyy,  \
                             tr_z_xxxzz_xxyz,  \
                             tr_z_xxxzz_xxz,   \
                             tr_z_xxxzz_xxzz,  \
                             tr_z_xxxzz_xyy,   \
                             tr_z_xxxzz_xyyy,  \
                             tr_z_xxxzz_xyyz,  \
                             tr_z_xxxzz_xyz,   \
                             tr_z_xxxzz_xyzz,  \
                             tr_z_xxxzz_xzz,   \
                             tr_z_xxxzz_xzzz,  \
                             tr_z_xxxzz_yyy,   \
                             tr_z_xxxzz_yyyy,  \
                             tr_z_xxxzz_yyyz,  \
                             tr_z_xxxzz_yyz,   \
                             tr_z_xxxzz_yyzz,  \
                             tr_z_xxxzz_yzz,   \
                             tr_z_xxxzz_yzzz,  \
                             tr_z_xxxzz_zzz,   \
                             tr_z_xxxzz_zzzz,  \
                             tr_z_xxzz_xxxx,   \
                             tr_z_xxzz_xxxy,   \
                             tr_z_xxzz_xxxz,   \
                             tr_z_xxzz_xxyy,   \
                             tr_z_xxzz_xxyz,   \
                             tr_z_xxzz_xxzz,   \
                             tr_z_xxzz_xyyy,   \
                             tr_z_xxzz_xyyz,   \
                             tr_z_xxzz_xyzz,   \
                             tr_z_xxzz_xzzz,   \
                             tr_z_xxzz_yyyy,   \
                             tr_z_xxzz_yyyz,   \
                             tr_z_xxzz_yyzz,   \
                             tr_z_xxzz_yzzz,   \
                             tr_z_xxzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxzz_xxxx[i] = 3.0 * tr_z_xxzz_xxxx[i] * fe_0 + 4.0 * tr_z_xxxzz_xxx[i] * fe_0 + tr_z_xxxzz_xxxx[i] * pa_x[i];

        tr_z_xxxxzz_xxxy[i] = 3.0 * tr_z_xxzz_xxxy[i] * fe_0 + 3.0 * tr_z_xxxzz_xxy[i] * fe_0 + tr_z_xxxzz_xxxy[i] * pa_x[i];

        tr_z_xxxxzz_xxxz[i] = 3.0 * tr_z_xxzz_xxxz[i] * fe_0 + 3.0 * tr_z_xxxzz_xxz[i] * fe_0 + tr_z_xxxzz_xxxz[i] * pa_x[i];

        tr_z_xxxxzz_xxyy[i] = 3.0 * tr_z_xxzz_xxyy[i] * fe_0 + 2.0 * tr_z_xxxzz_xyy[i] * fe_0 + tr_z_xxxzz_xxyy[i] * pa_x[i];

        tr_z_xxxxzz_xxyz[i] = 3.0 * tr_z_xxzz_xxyz[i] * fe_0 + 2.0 * tr_z_xxxzz_xyz[i] * fe_0 + tr_z_xxxzz_xxyz[i] * pa_x[i];

        tr_z_xxxxzz_xxzz[i] = 3.0 * tr_z_xxzz_xxzz[i] * fe_0 + 2.0 * tr_z_xxxzz_xzz[i] * fe_0 + tr_z_xxxzz_xxzz[i] * pa_x[i];

        tr_z_xxxxzz_xyyy[i] = 3.0 * tr_z_xxzz_xyyy[i] * fe_0 + tr_z_xxxzz_yyy[i] * fe_0 + tr_z_xxxzz_xyyy[i] * pa_x[i];

        tr_z_xxxxzz_xyyz[i] = 3.0 * tr_z_xxzz_xyyz[i] * fe_0 + tr_z_xxxzz_yyz[i] * fe_0 + tr_z_xxxzz_xyyz[i] * pa_x[i];

        tr_z_xxxxzz_xyzz[i] = 3.0 * tr_z_xxzz_xyzz[i] * fe_0 + tr_z_xxxzz_yzz[i] * fe_0 + tr_z_xxxzz_xyzz[i] * pa_x[i];

        tr_z_xxxxzz_xzzz[i] = 3.0 * tr_z_xxzz_xzzz[i] * fe_0 + tr_z_xxxzz_zzz[i] * fe_0 + tr_z_xxxzz_xzzz[i] * pa_x[i];

        tr_z_xxxxzz_yyyy[i] = 3.0 * tr_z_xxzz_yyyy[i] * fe_0 + tr_z_xxxzz_yyyy[i] * pa_x[i];

        tr_z_xxxxzz_yyyz[i] = 3.0 * tr_z_xxzz_yyyz[i] * fe_0 + tr_z_xxxzz_yyyz[i] * pa_x[i];

        tr_z_xxxxzz_yyzz[i] = 3.0 * tr_z_xxzz_yyzz[i] * fe_0 + tr_z_xxxzz_yyzz[i] * pa_x[i];

        tr_z_xxxxzz_yzzz[i] = 3.0 * tr_z_xxzz_yzzz[i] * fe_0 + tr_z_xxxzz_yzzz[i] * pa_x[i];

        tr_z_xxxxzz_zzzz[i] = 3.0 * tr_z_xxzz_zzzz[i] * fe_0 + tr_z_xxxzz_zzzz[i] * pa_x[i];
    }

    // Set up 930-945 components of targeted buffer : IG

    auto tr_z_xxxyyy_xxxx = pbuffer.data(idx_dip_ig + 930);

    auto tr_z_xxxyyy_xxxy = pbuffer.data(idx_dip_ig + 931);

    auto tr_z_xxxyyy_xxxz = pbuffer.data(idx_dip_ig + 932);

    auto tr_z_xxxyyy_xxyy = pbuffer.data(idx_dip_ig + 933);

    auto tr_z_xxxyyy_xxyz = pbuffer.data(idx_dip_ig + 934);

    auto tr_z_xxxyyy_xxzz = pbuffer.data(idx_dip_ig + 935);

    auto tr_z_xxxyyy_xyyy = pbuffer.data(idx_dip_ig + 936);

    auto tr_z_xxxyyy_xyyz = pbuffer.data(idx_dip_ig + 937);

    auto tr_z_xxxyyy_xyzz = pbuffer.data(idx_dip_ig + 938);

    auto tr_z_xxxyyy_xzzz = pbuffer.data(idx_dip_ig + 939);

    auto tr_z_xxxyyy_yyyy = pbuffer.data(idx_dip_ig + 940);

    auto tr_z_xxxyyy_yyyz = pbuffer.data(idx_dip_ig + 941);

    auto tr_z_xxxyyy_yyzz = pbuffer.data(idx_dip_ig + 942);

    auto tr_z_xxxyyy_yzzz = pbuffer.data(idx_dip_ig + 943);

    auto tr_z_xxxyyy_zzzz = pbuffer.data(idx_dip_ig + 944);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxy_xxxx,   \
                             tr_z_xxxy_xxxz,   \
                             tr_z_xxxy_xxzz,   \
                             tr_z_xxxy_xzzz,   \
                             tr_z_xxxyy_xxxx,  \
                             tr_z_xxxyy_xxxz,  \
                             tr_z_xxxyy_xxzz,  \
                             tr_z_xxxyy_xzzz,  \
                             tr_z_xxxyyy_xxxx, \
                             tr_z_xxxyyy_xxxy, \
                             tr_z_xxxyyy_xxxz, \
                             tr_z_xxxyyy_xxyy, \
                             tr_z_xxxyyy_xxyz, \
                             tr_z_xxxyyy_xxzz, \
                             tr_z_xxxyyy_xyyy, \
                             tr_z_xxxyyy_xyyz, \
                             tr_z_xxxyyy_xyzz, \
                             tr_z_xxxyyy_xzzz, \
                             tr_z_xxxyyy_yyyy, \
                             tr_z_xxxyyy_yyyz, \
                             tr_z_xxxyyy_yyzz, \
                             tr_z_xxxyyy_yzzz, \
                             tr_z_xxxyyy_zzzz, \
                             tr_z_xxyyy_xxxy,  \
                             tr_z_xxyyy_xxy,   \
                             tr_z_xxyyy_xxyy,  \
                             tr_z_xxyyy_xxyz,  \
                             tr_z_xxyyy_xyy,   \
                             tr_z_xxyyy_xyyy,  \
                             tr_z_xxyyy_xyyz,  \
                             tr_z_xxyyy_xyz,   \
                             tr_z_xxyyy_xyzz,  \
                             tr_z_xxyyy_yyy,   \
                             tr_z_xxyyy_yyyy,  \
                             tr_z_xxyyy_yyyz,  \
                             tr_z_xxyyy_yyz,   \
                             tr_z_xxyyy_yyzz,  \
                             tr_z_xxyyy_yzz,   \
                             tr_z_xxyyy_yzzz,  \
                             tr_z_xxyyy_zzzz,  \
                             tr_z_xyyy_xxxy,   \
                             tr_z_xyyy_xxyy,   \
                             tr_z_xyyy_xxyz,   \
                             tr_z_xyyy_xyyy,   \
                             tr_z_xyyy_xyyz,   \
                             tr_z_xyyy_xyzz,   \
                             tr_z_xyyy_yyyy,   \
                             tr_z_xyyy_yyyz,   \
                             tr_z_xyyy_yyzz,   \
                             tr_z_xyyy_yzzz,   \
                             tr_z_xyyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyy_xxxx[i] = 2.0 * tr_z_xxxy_xxxx[i] * fe_0 + tr_z_xxxyy_xxxx[i] * pa_y[i];

        tr_z_xxxyyy_xxxy[i] = 2.0 * tr_z_xyyy_xxxy[i] * fe_0 + 3.0 * tr_z_xxyyy_xxy[i] * fe_0 + tr_z_xxyyy_xxxy[i] * pa_x[i];

        tr_z_xxxyyy_xxxz[i] = 2.0 * tr_z_xxxy_xxxz[i] * fe_0 + tr_z_xxxyy_xxxz[i] * pa_y[i];

        tr_z_xxxyyy_xxyy[i] = 2.0 * tr_z_xyyy_xxyy[i] * fe_0 + 2.0 * tr_z_xxyyy_xyy[i] * fe_0 + tr_z_xxyyy_xxyy[i] * pa_x[i];

        tr_z_xxxyyy_xxyz[i] = 2.0 * tr_z_xyyy_xxyz[i] * fe_0 + 2.0 * tr_z_xxyyy_xyz[i] * fe_0 + tr_z_xxyyy_xxyz[i] * pa_x[i];

        tr_z_xxxyyy_xxzz[i] = 2.0 * tr_z_xxxy_xxzz[i] * fe_0 + tr_z_xxxyy_xxzz[i] * pa_y[i];

        tr_z_xxxyyy_xyyy[i] = 2.0 * tr_z_xyyy_xyyy[i] * fe_0 + tr_z_xxyyy_yyy[i] * fe_0 + tr_z_xxyyy_xyyy[i] * pa_x[i];

        tr_z_xxxyyy_xyyz[i] = 2.0 * tr_z_xyyy_xyyz[i] * fe_0 + tr_z_xxyyy_yyz[i] * fe_0 + tr_z_xxyyy_xyyz[i] * pa_x[i];

        tr_z_xxxyyy_xyzz[i] = 2.0 * tr_z_xyyy_xyzz[i] * fe_0 + tr_z_xxyyy_yzz[i] * fe_0 + tr_z_xxyyy_xyzz[i] * pa_x[i];

        tr_z_xxxyyy_xzzz[i] = 2.0 * tr_z_xxxy_xzzz[i] * fe_0 + tr_z_xxxyy_xzzz[i] * pa_y[i];

        tr_z_xxxyyy_yyyy[i] = 2.0 * tr_z_xyyy_yyyy[i] * fe_0 + tr_z_xxyyy_yyyy[i] * pa_x[i];

        tr_z_xxxyyy_yyyz[i] = 2.0 * tr_z_xyyy_yyyz[i] * fe_0 + tr_z_xxyyy_yyyz[i] * pa_x[i];

        tr_z_xxxyyy_yyzz[i] = 2.0 * tr_z_xyyy_yyzz[i] * fe_0 + tr_z_xxyyy_yyzz[i] * pa_x[i];

        tr_z_xxxyyy_yzzz[i] = 2.0 * tr_z_xyyy_yzzz[i] * fe_0 + tr_z_xxyyy_yzzz[i] * pa_x[i];

        tr_z_xxxyyy_zzzz[i] = 2.0 * tr_z_xyyy_zzzz[i] * fe_0 + tr_z_xxyyy_zzzz[i] * pa_x[i];
    }

    // Set up 945-960 components of targeted buffer : IG

    auto tr_z_xxxyyz_xxxx = pbuffer.data(idx_dip_ig + 945);

    auto tr_z_xxxyyz_xxxy = pbuffer.data(idx_dip_ig + 946);

    auto tr_z_xxxyyz_xxxz = pbuffer.data(idx_dip_ig + 947);

    auto tr_z_xxxyyz_xxyy = pbuffer.data(idx_dip_ig + 948);

    auto tr_z_xxxyyz_xxyz = pbuffer.data(idx_dip_ig + 949);

    auto tr_z_xxxyyz_xxzz = pbuffer.data(idx_dip_ig + 950);

    auto tr_z_xxxyyz_xyyy = pbuffer.data(idx_dip_ig + 951);

    auto tr_z_xxxyyz_xyyz = pbuffer.data(idx_dip_ig + 952);

    auto tr_z_xxxyyz_xyzz = pbuffer.data(idx_dip_ig + 953);

    auto tr_z_xxxyyz_xzzz = pbuffer.data(idx_dip_ig + 954);

    auto tr_z_xxxyyz_yyyy = pbuffer.data(idx_dip_ig + 955);

    auto tr_z_xxxyyz_yyyz = pbuffer.data(idx_dip_ig + 956);

    auto tr_z_xxxyyz_yyzz = pbuffer.data(idx_dip_ig + 957);

    auto tr_z_xxxyyz_yzzz = pbuffer.data(idx_dip_ig + 958);

    auto tr_z_xxxyyz_zzzz = pbuffer.data(idx_dip_ig + 959);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_z_xxxyy_xxxy,  \
                             tr_z_xxxyy_xxyy,  \
                             tr_z_xxxyy_xyyy,  \
                             tr_z_xxxyyz_xxxx, \
                             tr_z_xxxyyz_xxxy, \
                             tr_z_xxxyyz_xxxz, \
                             tr_z_xxxyyz_xxyy, \
                             tr_z_xxxyyz_xxyz, \
                             tr_z_xxxyyz_xxzz, \
                             tr_z_xxxyyz_xyyy, \
                             tr_z_xxxyyz_xyyz, \
                             tr_z_xxxyyz_xyzz, \
                             tr_z_xxxyyz_xzzz, \
                             tr_z_xxxyyz_yyyy, \
                             tr_z_xxxyyz_yyyz, \
                             tr_z_xxxyyz_yyzz, \
                             tr_z_xxxyyz_yzzz, \
                             tr_z_xxxyyz_zzzz, \
                             tr_z_xxxyz_xxxx,  \
                             tr_z_xxxyz_xxxz,  \
                             tr_z_xxxyz_xxzz,  \
                             tr_z_xxxyz_xzzz,  \
                             tr_z_xxxz_xxxx,   \
                             tr_z_xxxz_xxxz,   \
                             tr_z_xxxz_xxzz,   \
                             tr_z_xxxz_xzzz,   \
                             tr_z_xxyyz_xxyz,  \
                             tr_z_xxyyz_xyyz,  \
                             tr_z_xxyyz_xyz,   \
                             tr_z_xxyyz_xyzz,  \
                             tr_z_xxyyz_yyyy,  \
                             tr_z_xxyyz_yyyz,  \
                             tr_z_xxyyz_yyz,   \
                             tr_z_xxyyz_yyzz,  \
                             tr_z_xxyyz_yzz,   \
                             tr_z_xxyyz_yzzz,  \
                             tr_z_xxyyz_zzzz,  \
                             tr_z_xyyz_xxyz,   \
                             tr_z_xyyz_xyyz,   \
                             tr_z_xyyz_xyzz,   \
                             tr_z_xyyz_yyyy,   \
                             tr_z_xyyz_yyyz,   \
                             tr_z_xyyz_yyzz,   \
                             tr_z_xyyz_yzzz,   \
                             tr_z_xyyz_zzzz,   \
                             ts_xxxyy_xxxy,    \
                             ts_xxxyy_xxyy,    \
                             ts_xxxyy_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyz_xxxx[i] = tr_z_xxxz_xxxx[i] * fe_0 + tr_z_xxxyz_xxxx[i] * pa_y[i];

        tr_z_xxxyyz_xxxy[i] = ts_xxxyy_xxxy[i] * fe_0 + tr_z_xxxyy_xxxy[i] * pa_z[i];

        tr_z_xxxyyz_xxxz[i] = tr_z_xxxz_xxxz[i] * fe_0 + tr_z_xxxyz_xxxz[i] * pa_y[i];

        tr_z_xxxyyz_xxyy[i] = ts_xxxyy_xxyy[i] * fe_0 + tr_z_xxxyy_xxyy[i] * pa_z[i];

        tr_z_xxxyyz_xxyz[i] = 2.0 * tr_z_xyyz_xxyz[i] * fe_0 + 2.0 * tr_z_xxyyz_xyz[i] * fe_0 + tr_z_xxyyz_xxyz[i] * pa_x[i];

        tr_z_xxxyyz_xxzz[i] = tr_z_xxxz_xxzz[i] * fe_0 + tr_z_xxxyz_xxzz[i] * pa_y[i];

        tr_z_xxxyyz_xyyy[i] = ts_xxxyy_xyyy[i] * fe_0 + tr_z_xxxyy_xyyy[i] * pa_z[i];

        tr_z_xxxyyz_xyyz[i] = 2.0 * tr_z_xyyz_xyyz[i] * fe_0 + tr_z_xxyyz_yyz[i] * fe_0 + tr_z_xxyyz_xyyz[i] * pa_x[i];

        tr_z_xxxyyz_xyzz[i] = 2.0 * tr_z_xyyz_xyzz[i] * fe_0 + tr_z_xxyyz_yzz[i] * fe_0 + tr_z_xxyyz_xyzz[i] * pa_x[i];

        tr_z_xxxyyz_xzzz[i] = tr_z_xxxz_xzzz[i] * fe_0 + tr_z_xxxyz_xzzz[i] * pa_y[i];

        tr_z_xxxyyz_yyyy[i] = 2.0 * tr_z_xyyz_yyyy[i] * fe_0 + tr_z_xxyyz_yyyy[i] * pa_x[i];

        tr_z_xxxyyz_yyyz[i] = 2.0 * tr_z_xyyz_yyyz[i] * fe_0 + tr_z_xxyyz_yyyz[i] * pa_x[i];

        tr_z_xxxyyz_yyzz[i] = 2.0 * tr_z_xyyz_yyzz[i] * fe_0 + tr_z_xxyyz_yyzz[i] * pa_x[i];

        tr_z_xxxyyz_yzzz[i] = 2.0 * tr_z_xyyz_yzzz[i] * fe_0 + tr_z_xxyyz_yzzz[i] * pa_x[i];

        tr_z_xxxyyz_zzzz[i] = 2.0 * tr_z_xyyz_zzzz[i] * fe_0 + tr_z_xxyyz_zzzz[i] * pa_x[i];
    }

    // Set up 960-975 components of targeted buffer : IG

    auto tr_z_xxxyzz_xxxx = pbuffer.data(idx_dip_ig + 960);

    auto tr_z_xxxyzz_xxxy = pbuffer.data(idx_dip_ig + 961);

    auto tr_z_xxxyzz_xxxz = pbuffer.data(idx_dip_ig + 962);

    auto tr_z_xxxyzz_xxyy = pbuffer.data(idx_dip_ig + 963);

    auto tr_z_xxxyzz_xxyz = pbuffer.data(idx_dip_ig + 964);

    auto tr_z_xxxyzz_xxzz = pbuffer.data(idx_dip_ig + 965);

    auto tr_z_xxxyzz_xyyy = pbuffer.data(idx_dip_ig + 966);

    auto tr_z_xxxyzz_xyyz = pbuffer.data(idx_dip_ig + 967);

    auto tr_z_xxxyzz_xyzz = pbuffer.data(idx_dip_ig + 968);

    auto tr_z_xxxyzz_xzzz = pbuffer.data(idx_dip_ig + 969);

    auto tr_z_xxxyzz_yyyy = pbuffer.data(idx_dip_ig + 970);

    auto tr_z_xxxyzz_yyyz = pbuffer.data(idx_dip_ig + 971);

    auto tr_z_xxxyzz_yyzz = pbuffer.data(idx_dip_ig + 972);

    auto tr_z_xxxyzz_yzzz = pbuffer.data(idx_dip_ig + 973);

    auto tr_z_xxxyzz_zzzz = pbuffer.data(idx_dip_ig + 974);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxyzz_xxxx, \
                             tr_z_xxxyzz_xxxy, \
                             tr_z_xxxyzz_xxxz, \
                             tr_z_xxxyzz_xxyy, \
                             tr_z_xxxyzz_xxyz, \
                             tr_z_xxxyzz_xxzz, \
                             tr_z_xxxyzz_xyyy, \
                             tr_z_xxxyzz_xyyz, \
                             tr_z_xxxyzz_xyzz, \
                             tr_z_xxxyzz_xzzz, \
                             tr_z_xxxyzz_yyyy, \
                             tr_z_xxxyzz_yyyz, \
                             tr_z_xxxyzz_yyzz, \
                             tr_z_xxxyzz_yzzz, \
                             tr_z_xxxyzz_zzzz, \
                             tr_z_xxxzz_xxx,   \
                             tr_z_xxxzz_xxxx,  \
                             tr_z_xxxzz_xxxy,  \
                             tr_z_xxxzz_xxxz,  \
                             tr_z_xxxzz_xxy,   \
                             tr_z_xxxzz_xxyy,  \
                             tr_z_xxxzz_xxyz,  \
                             tr_z_xxxzz_xxz,   \
                             tr_z_xxxzz_xxzz,  \
                             tr_z_xxxzz_xyy,   \
                             tr_z_xxxzz_xyyy,  \
                             tr_z_xxxzz_xyyz,  \
                             tr_z_xxxzz_xyz,   \
                             tr_z_xxxzz_xyzz,  \
                             tr_z_xxxzz_xzz,   \
                             tr_z_xxxzz_xzzz,  \
                             tr_z_xxxzz_zzzz,  \
                             tr_z_xxyzz_yyyy,  \
                             tr_z_xxyzz_yyyz,  \
                             tr_z_xxyzz_yyzz,  \
                             tr_z_xxyzz_yzzz,  \
                             tr_z_xyzz_yyyy,   \
                             tr_z_xyzz_yyyz,   \
                             tr_z_xyzz_yyzz,   \
                             tr_z_xyzz_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyzz_xxxx[i] = tr_z_xxxzz_xxxx[i] * pa_y[i];

        tr_z_xxxyzz_xxxy[i] = tr_z_xxxzz_xxx[i] * fe_0 + tr_z_xxxzz_xxxy[i] * pa_y[i];

        tr_z_xxxyzz_xxxz[i] = tr_z_xxxzz_xxxz[i] * pa_y[i];

        tr_z_xxxyzz_xxyy[i] = 2.0 * tr_z_xxxzz_xxy[i] * fe_0 + tr_z_xxxzz_xxyy[i] * pa_y[i];

        tr_z_xxxyzz_xxyz[i] = tr_z_xxxzz_xxz[i] * fe_0 + tr_z_xxxzz_xxyz[i] * pa_y[i];

        tr_z_xxxyzz_xxzz[i] = tr_z_xxxzz_xxzz[i] * pa_y[i];

        tr_z_xxxyzz_xyyy[i] = 3.0 * tr_z_xxxzz_xyy[i] * fe_0 + tr_z_xxxzz_xyyy[i] * pa_y[i];

        tr_z_xxxyzz_xyyz[i] = 2.0 * tr_z_xxxzz_xyz[i] * fe_0 + tr_z_xxxzz_xyyz[i] * pa_y[i];

        tr_z_xxxyzz_xyzz[i] = tr_z_xxxzz_xzz[i] * fe_0 + tr_z_xxxzz_xyzz[i] * pa_y[i];

        tr_z_xxxyzz_xzzz[i] = tr_z_xxxzz_xzzz[i] * pa_y[i];

        tr_z_xxxyzz_yyyy[i] = 2.0 * tr_z_xyzz_yyyy[i] * fe_0 + tr_z_xxyzz_yyyy[i] * pa_x[i];

        tr_z_xxxyzz_yyyz[i] = 2.0 * tr_z_xyzz_yyyz[i] * fe_0 + tr_z_xxyzz_yyyz[i] * pa_x[i];

        tr_z_xxxyzz_yyzz[i] = 2.0 * tr_z_xyzz_yyzz[i] * fe_0 + tr_z_xxyzz_yyzz[i] * pa_x[i];

        tr_z_xxxyzz_yzzz[i] = 2.0 * tr_z_xyzz_yzzz[i] * fe_0 + tr_z_xxyzz_yzzz[i] * pa_x[i];

        tr_z_xxxyzz_zzzz[i] = tr_z_xxxzz_zzzz[i] * pa_y[i];
    }

    // Set up 975-990 components of targeted buffer : IG

    auto tr_z_xxxzzz_xxxx = pbuffer.data(idx_dip_ig + 975);

    auto tr_z_xxxzzz_xxxy = pbuffer.data(idx_dip_ig + 976);

    auto tr_z_xxxzzz_xxxz = pbuffer.data(idx_dip_ig + 977);

    auto tr_z_xxxzzz_xxyy = pbuffer.data(idx_dip_ig + 978);

    auto tr_z_xxxzzz_xxyz = pbuffer.data(idx_dip_ig + 979);

    auto tr_z_xxxzzz_xxzz = pbuffer.data(idx_dip_ig + 980);

    auto tr_z_xxxzzz_xyyy = pbuffer.data(idx_dip_ig + 981);

    auto tr_z_xxxzzz_xyyz = pbuffer.data(idx_dip_ig + 982);

    auto tr_z_xxxzzz_xyzz = pbuffer.data(idx_dip_ig + 983);

    auto tr_z_xxxzzz_xzzz = pbuffer.data(idx_dip_ig + 984);

    auto tr_z_xxxzzz_yyyy = pbuffer.data(idx_dip_ig + 985);

    auto tr_z_xxxzzz_yyyz = pbuffer.data(idx_dip_ig + 986);

    auto tr_z_xxxzzz_yyzz = pbuffer.data(idx_dip_ig + 987);

    auto tr_z_xxxzzz_yzzz = pbuffer.data(idx_dip_ig + 988);

    auto tr_z_xxxzzz_zzzz = pbuffer.data(idx_dip_ig + 989);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxxzzz_xxxx, \
                             tr_z_xxxzzz_xxxy, \
                             tr_z_xxxzzz_xxxz, \
                             tr_z_xxxzzz_xxyy, \
                             tr_z_xxxzzz_xxyz, \
                             tr_z_xxxzzz_xxzz, \
                             tr_z_xxxzzz_xyyy, \
                             tr_z_xxxzzz_xyyz, \
                             tr_z_xxxzzz_xyzz, \
                             tr_z_xxxzzz_xzzz, \
                             tr_z_xxxzzz_yyyy, \
                             tr_z_xxxzzz_yyyz, \
                             tr_z_xxxzzz_yyzz, \
                             tr_z_xxxzzz_yzzz, \
                             tr_z_xxxzzz_zzzz, \
                             tr_z_xxzzz_xxx,   \
                             tr_z_xxzzz_xxxx,  \
                             tr_z_xxzzz_xxxy,  \
                             tr_z_xxzzz_xxxz,  \
                             tr_z_xxzzz_xxy,   \
                             tr_z_xxzzz_xxyy,  \
                             tr_z_xxzzz_xxyz,  \
                             tr_z_xxzzz_xxz,   \
                             tr_z_xxzzz_xxzz,  \
                             tr_z_xxzzz_xyy,   \
                             tr_z_xxzzz_xyyy,  \
                             tr_z_xxzzz_xyyz,  \
                             tr_z_xxzzz_xyz,   \
                             tr_z_xxzzz_xyzz,  \
                             tr_z_xxzzz_xzz,   \
                             tr_z_xxzzz_xzzz,  \
                             tr_z_xxzzz_yyy,   \
                             tr_z_xxzzz_yyyy,  \
                             tr_z_xxzzz_yyyz,  \
                             tr_z_xxzzz_yyz,   \
                             tr_z_xxzzz_yyzz,  \
                             tr_z_xxzzz_yzz,   \
                             tr_z_xxzzz_yzzz,  \
                             tr_z_xxzzz_zzz,   \
                             tr_z_xxzzz_zzzz,  \
                             tr_z_xzzz_xxxx,   \
                             tr_z_xzzz_xxxy,   \
                             tr_z_xzzz_xxxz,   \
                             tr_z_xzzz_xxyy,   \
                             tr_z_xzzz_xxyz,   \
                             tr_z_xzzz_xxzz,   \
                             tr_z_xzzz_xyyy,   \
                             tr_z_xzzz_xyyz,   \
                             tr_z_xzzz_xyzz,   \
                             tr_z_xzzz_xzzz,   \
                             tr_z_xzzz_yyyy,   \
                             tr_z_xzzz_yyyz,   \
                             tr_z_xzzz_yyzz,   \
                             tr_z_xzzz_yzzz,   \
                             tr_z_xzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzzz_xxxx[i] = 2.0 * tr_z_xzzz_xxxx[i] * fe_0 + 4.0 * tr_z_xxzzz_xxx[i] * fe_0 + tr_z_xxzzz_xxxx[i] * pa_x[i];

        tr_z_xxxzzz_xxxy[i] = 2.0 * tr_z_xzzz_xxxy[i] * fe_0 + 3.0 * tr_z_xxzzz_xxy[i] * fe_0 + tr_z_xxzzz_xxxy[i] * pa_x[i];

        tr_z_xxxzzz_xxxz[i] = 2.0 * tr_z_xzzz_xxxz[i] * fe_0 + 3.0 * tr_z_xxzzz_xxz[i] * fe_0 + tr_z_xxzzz_xxxz[i] * pa_x[i];

        tr_z_xxxzzz_xxyy[i] = 2.0 * tr_z_xzzz_xxyy[i] * fe_0 + 2.0 * tr_z_xxzzz_xyy[i] * fe_0 + tr_z_xxzzz_xxyy[i] * pa_x[i];

        tr_z_xxxzzz_xxyz[i] = 2.0 * tr_z_xzzz_xxyz[i] * fe_0 + 2.0 * tr_z_xxzzz_xyz[i] * fe_0 + tr_z_xxzzz_xxyz[i] * pa_x[i];

        tr_z_xxxzzz_xxzz[i] = 2.0 * tr_z_xzzz_xxzz[i] * fe_0 + 2.0 * tr_z_xxzzz_xzz[i] * fe_0 + tr_z_xxzzz_xxzz[i] * pa_x[i];

        tr_z_xxxzzz_xyyy[i] = 2.0 * tr_z_xzzz_xyyy[i] * fe_0 + tr_z_xxzzz_yyy[i] * fe_0 + tr_z_xxzzz_xyyy[i] * pa_x[i];

        tr_z_xxxzzz_xyyz[i] = 2.0 * tr_z_xzzz_xyyz[i] * fe_0 + tr_z_xxzzz_yyz[i] * fe_0 + tr_z_xxzzz_xyyz[i] * pa_x[i];

        tr_z_xxxzzz_xyzz[i] = 2.0 * tr_z_xzzz_xyzz[i] * fe_0 + tr_z_xxzzz_yzz[i] * fe_0 + tr_z_xxzzz_xyzz[i] * pa_x[i];

        tr_z_xxxzzz_xzzz[i] = 2.0 * tr_z_xzzz_xzzz[i] * fe_0 + tr_z_xxzzz_zzz[i] * fe_0 + tr_z_xxzzz_xzzz[i] * pa_x[i];

        tr_z_xxxzzz_yyyy[i] = 2.0 * tr_z_xzzz_yyyy[i] * fe_0 + tr_z_xxzzz_yyyy[i] * pa_x[i];

        tr_z_xxxzzz_yyyz[i] = 2.0 * tr_z_xzzz_yyyz[i] * fe_0 + tr_z_xxzzz_yyyz[i] * pa_x[i];

        tr_z_xxxzzz_yyzz[i] = 2.0 * tr_z_xzzz_yyzz[i] * fe_0 + tr_z_xxzzz_yyzz[i] * pa_x[i];

        tr_z_xxxzzz_yzzz[i] = 2.0 * tr_z_xzzz_yzzz[i] * fe_0 + tr_z_xxzzz_yzzz[i] * pa_x[i];

        tr_z_xxxzzz_zzzz[i] = 2.0 * tr_z_xzzz_zzzz[i] * fe_0 + tr_z_xxzzz_zzzz[i] * pa_x[i];
    }

    // Set up 990-1005 components of targeted buffer : IG

    auto tr_z_xxyyyy_xxxx = pbuffer.data(idx_dip_ig + 990);

    auto tr_z_xxyyyy_xxxy = pbuffer.data(idx_dip_ig + 991);

    auto tr_z_xxyyyy_xxxz = pbuffer.data(idx_dip_ig + 992);

    auto tr_z_xxyyyy_xxyy = pbuffer.data(idx_dip_ig + 993);

    auto tr_z_xxyyyy_xxyz = pbuffer.data(idx_dip_ig + 994);

    auto tr_z_xxyyyy_xxzz = pbuffer.data(idx_dip_ig + 995);

    auto tr_z_xxyyyy_xyyy = pbuffer.data(idx_dip_ig + 996);

    auto tr_z_xxyyyy_xyyz = pbuffer.data(idx_dip_ig + 997);

    auto tr_z_xxyyyy_xyzz = pbuffer.data(idx_dip_ig + 998);

    auto tr_z_xxyyyy_xzzz = pbuffer.data(idx_dip_ig + 999);

    auto tr_z_xxyyyy_yyyy = pbuffer.data(idx_dip_ig + 1000);

    auto tr_z_xxyyyy_yyyz = pbuffer.data(idx_dip_ig + 1001);

    auto tr_z_xxyyyy_yyzz = pbuffer.data(idx_dip_ig + 1002);

    auto tr_z_xxyyyy_yzzz = pbuffer.data(idx_dip_ig + 1003);

    auto tr_z_xxyyyy_zzzz = pbuffer.data(idx_dip_ig + 1004);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxyy_xxxx,   \
                             tr_z_xxyy_xxxz,   \
                             tr_z_xxyy_xxzz,   \
                             tr_z_xxyy_xzzz,   \
                             tr_z_xxyyy_xxxx,  \
                             tr_z_xxyyy_xxxz,  \
                             tr_z_xxyyy_xxzz,  \
                             tr_z_xxyyy_xzzz,  \
                             tr_z_xxyyyy_xxxx, \
                             tr_z_xxyyyy_xxxy, \
                             tr_z_xxyyyy_xxxz, \
                             tr_z_xxyyyy_xxyy, \
                             tr_z_xxyyyy_xxyz, \
                             tr_z_xxyyyy_xxzz, \
                             tr_z_xxyyyy_xyyy, \
                             tr_z_xxyyyy_xyyz, \
                             tr_z_xxyyyy_xyzz, \
                             tr_z_xxyyyy_xzzz, \
                             tr_z_xxyyyy_yyyy, \
                             tr_z_xxyyyy_yyyz, \
                             tr_z_xxyyyy_yyzz, \
                             tr_z_xxyyyy_yzzz, \
                             tr_z_xxyyyy_zzzz, \
                             tr_z_xyyyy_xxxy,  \
                             tr_z_xyyyy_xxy,   \
                             tr_z_xyyyy_xxyy,  \
                             tr_z_xyyyy_xxyz,  \
                             tr_z_xyyyy_xyy,   \
                             tr_z_xyyyy_xyyy,  \
                             tr_z_xyyyy_xyyz,  \
                             tr_z_xyyyy_xyz,   \
                             tr_z_xyyyy_xyzz,  \
                             tr_z_xyyyy_yyy,   \
                             tr_z_xyyyy_yyyy,  \
                             tr_z_xyyyy_yyyz,  \
                             tr_z_xyyyy_yyz,   \
                             tr_z_xyyyy_yyzz,  \
                             tr_z_xyyyy_yzz,   \
                             tr_z_xyyyy_yzzz,  \
                             tr_z_xyyyy_zzzz,  \
                             tr_z_yyyy_xxxy,   \
                             tr_z_yyyy_xxyy,   \
                             tr_z_yyyy_xxyz,   \
                             tr_z_yyyy_xyyy,   \
                             tr_z_yyyy_xyyz,   \
                             tr_z_yyyy_xyzz,   \
                             tr_z_yyyy_yyyy,   \
                             tr_z_yyyy_yyyz,   \
                             tr_z_yyyy_yyzz,   \
                             tr_z_yyyy_yzzz,   \
                             tr_z_yyyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyy_xxxx[i] = 3.0 * tr_z_xxyy_xxxx[i] * fe_0 + tr_z_xxyyy_xxxx[i] * pa_y[i];

        tr_z_xxyyyy_xxxy[i] = tr_z_yyyy_xxxy[i] * fe_0 + 3.0 * tr_z_xyyyy_xxy[i] * fe_0 + tr_z_xyyyy_xxxy[i] * pa_x[i];

        tr_z_xxyyyy_xxxz[i] = 3.0 * tr_z_xxyy_xxxz[i] * fe_0 + tr_z_xxyyy_xxxz[i] * pa_y[i];

        tr_z_xxyyyy_xxyy[i] = tr_z_yyyy_xxyy[i] * fe_0 + 2.0 * tr_z_xyyyy_xyy[i] * fe_0 + tr_z_xyyyy_xxyy[i] * pa_x[i];

        tr_z_xxyyyy_xxyz[i] = tr_z_yyyy_xxyz[i] * fe_0 + 2.0 * tr_z_xyyyy_xyz[i] * fe_0 + tr_z_xyyyy_xxyz[i] * pa_x[i];

        tr_z_xxyyyy_xxzz[i] = 3.0 * tr_z_xxyy_xxzz[i] * fe_0 + tr_z_xxyyy_xxzz[i] * pa_y[i];

        tr_z_xxyyyy_xyyy[i] = tr_z_yyyy_xyyy[i] * fe_0 + tr_z_xyyyy_yyy[i] * fe_0 + tr_z_xyyyy_xyyy[i] * pa_x[i];

        tr_z_xxyyyy_xyyz[i] = tr_z_yyyy_xyyz[i] * fe_0 + tr_z_xyyyy_yyz[i] * fe_0 + tr_z_xyyyy_xyyz[i] * pa_x[i];

        tr_z_xxyyyy_xyzz[i] = tr_z_yyyy_xyzz[i] * fe_0 + tr_z_xyyyy_yzz[i] * fe_0 + tr_z_xyyyy_xyzz[i] * pa_x[i];

        tr_z_xxyyyy_xzzz[i] = 3.0 * tr_z_xxyy_xzzz[i] * fe_0 + tr_z_xxyyy_xzzz[i] * pa_y[i];

        tr_z_xxyyyy_yyyy[i] = tr_z_yyyy_yyyy[i] * fe_0 + tr_z_xyyyy_yyyy[i] * pa_x[i];

        tr_z_xxyyyy_yyyz[i] = tr_z_yyyy_yyyz[i] * fe_0 + tr_z_xyyyy_yyyz[i] * pa_x[i];

        tr_z_xxyyyy_yyzz[i] = tr_z_yyyy_yyzz[i] * fe_0 + tr_z_xyyyy_yyzz[i] * pa_x[i];

        tr_z_xxyyyy_yzzz[i] = tr_z_yyyy_yzzz[i] * fe_0 + tr_z_xyyyy_yzzz[i] * pa_x[i];

        tr_z_xxyyyy_zzzz[i] = tr_z_yyyy_zzzz[i] * fe_0 + tr_z_xyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 1005-1020 components of targeted buffer : IG

    auto tr_z_xxyyyz_xxxx = pbuffer.data(idx_dip_ig + 1005);

    auto tr_z_xxyyyz_xxxy = pbuffer.data(idx_dip_ig + 1006);

    auto tr_z_xxyyyz_xxxz = pbuffer.data(idx_dip_ig + 1007);

    auto tr_z_xxyyyz_xxyy = pbuffer.data(idx_dip_ig + 1008);

    auto tr_z_xxyyyz_xxyz = pbuffer.data(idx_dip_ig + 1009);

    auto tr_z_xxyyyz_xxzz = pbuffer.data(idx_dip_ig + 1010);

    auto tr_z_xxyyyz_xyyy = pbuffer.data(idx_dip_ig + 1011);

    auto tr_z_xxyyyz_xyyz = pbuffer.data(idx_dip_ig + 1012);

    auto tr_z_xxyyyz_xyzz = pbuffer.data(idx_dip_ig + 1013);

    auto tr_z_xxyyyz_xzzz = pbuffer.data(idx_dip_ig + 1014);

    auto tr_z_xxyyyz_yyyy = pbuffer.data(idx_dip_ig + 1015);

    auto tr_z_xxyyyz_yyyz = pbuffer.data(idx_dip_ig + 1016);

    auto tr_z_xxyyyz_yyzz = pbuffer.data(idx_dip_ig + 1017);

    auto tr_z_xxyyyz_yzzz = pbuffer.data(idx_dip_ig + 1018);

    auto tr_z_xxyyyz_zzzz = pbuffer.data(idx_dip_ig + 1019);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_z_xxyyy_xxxy,  \
                             tr_z_xxyyy_xxyy,  \
                             tr_z_xxyyy_xyyy,  \
                             tr_z_xxyyyz_xxxx, \
                             tr_z_xxyyyz_xxxy, \
                             tr_z_xxyyyz_xxxz, \
                             tr_z_xxyyyz_xxyy, \
                             tr_z_xxyyyz_xxyz, \
                             tr_z_xxyyyz_xxzz, \
                             tr_z_xxyyyz_xyyy, \
                             tr_z_xxyyyz_xyyz, \
                             tr_z_xxyyyz_xyzz, \
                             tr_z_xxyyyz_xzzz, \
                             tr_z_xxyyyz_yyyy, \
                             tr_z_xxyyyz_yyyz, \
                             tr_z_xxyyyz_yyzz, \
                             tr_z_xxyyyz_yzzz, \
                             tr_z_xxyyyz_zzzz, \
                             tr_z_xxyyz_xxxx,  \
                             tr_z_xxyyz_xxxz,  \
                             tr_z_xxyyz_xxzz,  \
                             tr_z_xxyyz_xzzz,  \
                             tr_z_xxyz_xxxx,   \
                             tr_z_xxyz_xxxz,   \
                             tr_z_xxyz_xxzz,   \
                             tr_z_xxyz_xzzz,   \
                             tr_z_xyyyz_xxyz,  \
                             tr_z_xyyyz_xyyz,  \
                             tr_z_xyyyz_xyz,   \
                             tr_z_xyyyz_xyzz,  \
                             tr_z_xyyyz_yyyy,  \
                             tr_z_xyyyz_yyyz,  \
                             tr_z_xyyyz_yyz,   \
                             tr_z_xyyyz_yyzz,  \
                             tr_z_xyyyz_yzz,   \
                             tr_z_xyyyz_yzzz,  \
                             tr_z_xyyyz_zzzz,  \
                             tr_z_yyyz_xxyz,   \
                             tr_z_yyyz_xyyz,   \
                             tr_z_yyyz_xyzz,   \
                             tr_z_yyyz_yyyy,   \
                             tr_z_yyyz_yyyz,   \
                             tr_z_yyyz_yyzz,   \
                             tr_z_yyyz_yzzz,   \
                             tr_z_yyyz_zzzz,   \
                             ts_xxyyy_xxxy,    \
                             ts_xxyyy_xxyy,    \
                             ts_xxyyy_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyz_xxxx[i] = 2.0 * tr_z_xxyz_xxxx[i] * fe_0 + tr_z_xxyyz_xxxx[i] * pa_y[i];

        tr_z_xxyyyz_xxxy[i] = ts_xxyyy_xxxy[i] * fe_0 + tr_z_xxyyy_xxxy[i] * pa_z[i];

        tr_z_xxyyyz_xxxz[i] = 2.0 * tr_z_xxyz_xxxz[i] * fe_0 + tr_z_xxyyz_xxxz[i] * pa_y[i];

        tr_z_xxyyyz_xxyy[i] = ts_xxyyy_xxyy[i] * fe_0 + tr_z_xxyyy_xxyy[i] * pa_z[i];

        tr_z_xxyyyz_xxyz[i] = tr_z_yyyz_xxyz[i] * fe_0 + 2.0 * tr_z_xyyyz_xyz[i] * fe_0 + tr_z_xyyyz_xxyz[i] * pa_x[i];

        tr_z_xxyyyz_xxzz[i] = 2.0 * tr_z_xxyz_xxzz[i] * fe_0 + tr_z_xxyyz_xxzz[i] * pa_y[i];

        tr_z_xxyyyz_xyyy[i] = ts_xxyyy_xyyy[i] * fe_0 + tr_z_xxyyy_xyyy[i] * pa_z[i];

        tr_z_xxyyyz_xyyz[i] = tr_z_yyyz_xyyz[i] * fe_0 + tr_z_xyyyz_yyz[i] * fe_0 + tr_z_xyyyz_xyyz[i] * pa_x[i];

        tr_z_xxyyyz_xyzz[i] = tr_z_yyyz_xyzz[i] * fe_0 + tr_z_xyyyz_yzz[i] * fe_0 + tr_z_xyyyz_xyzz[i] * pa_x[i];

        tr_z_xxyyyz_xzzz[i] = 2.0 * tr_z_xxyz_xzzz[i] * fe_0 + tr_z_xxyyz_xzzz[i] * pa_y[i];

        tr_z_xxyyyz_yyyy[i] = tr_z_yyyz_yyyy[i] * fe_0 + tr_z_xyyyz_yyyy[i] * pa_x[i];

        tr_z_xxyyyz_yyyz[i] = tr_z_yyyz_yyyz[i] * fe_0 + tr_z_xyyyz_yyyz[i] * pa_x[i];

        tr_z_xxyyyz_yyzz[i] = tr_z_yyyz_yyzz[i] * fe_0 + tr_z_xyyyz_yyzz[i] * pa_x[i];

        tr_z_xxyyyz_yzzz[i] = tr_z_yyyz_yzzz[i] * fe_0 + tr_z_xyyyz_yzzz[i] * pa_x[i];

        tr_z_xxyyyz_zzzz[i] = tr_z_yyyz_zzzz[i] * fe_0 + tr_z_xyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 1020-1035 components of targeted buffer : IG

    auto tr_z_xxyyzz_xxxx = pbuffer.data(idx_dip_ig + 1020);

    auto tr_z_xxyyzz_xxxy = pbuffer.data(idx_dip_ig + 1021);

    auto tr_z_xxyyzz_xxxz = pbuffer.data(idx_dip_ig + 1022);

    auto tr_z_xxyyzz_xxyy = pbuffer.data(idx_dip_ig + 1023);

    auto tr_z_xxyyzz_xxyz = pbuffer.data(idx_dip_ig + 1024);

    auto tr_z_xxyyzz_xxzz = pbuffer.data(idx_dip_ig + 1025);

    auto tr_z_xxyyzz_xyyy = pbuffer.data(idx_dip_ig + 1026);

    auto tr_z_xxyyzz_xyyz = pbuffer.data(idx_dip_ig + 1027);

    auto tr_z_xxyyzz_xyzz = pbuffer.data(idx_dip_ig + 1028);

    auto tr_z_xxyyzz_xzzz = pbuffer.data(idx_dip_ig + 1029);

    auto tr_z_xxyyzz_yyyy = pbuffer.data(idx_dip_ig + 1030);

    auto tr_z_xxyyzz_yyyz = pbuffer.data(idx_dip_ig + 1031);

    auto tr_z_xxyyzz_yyzz = pbuffer.data(idx_dip_ig + 1032);

    auto tr_z_xxyyzz_yzzz = pbuffer.data(idx_dip_ig + 1033);

    auto tr_z_xxyyzz_zzzz = pbuffer.data(idx_dip_ig + 1034);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxyyzz_xxxx, \
                             tr_z_xxyyzz_xxxy, \
                             tr_z_xxyyzz_xxxz, \
                             tr_z_xxyyzz_xxyy, \
                             tr_z_xxyyzz_xxyz, \
                             tr_z_xxyyzz_xxzz, \
                             tr_z_xxyyzz_xyyy, \
                             tr_z_xxyyzz_xyyz, \
                             tr_z_xxyyzz_xyzz, \
                             tr_z_xxyyzz_xzzz, \
                             tr_z_xxyyzz_yyyy, \
                             tr_z_xxyyzz_yyyz, \
                             tr_z_xxyyzz_yyzz, \
                             tr_z_xxyyzz_yzzz, \
                             tr_z_xxyyzz_zzzz, \
                             tr_z_xxyzz_xxxx,  \
                             tr_z_xxyzz_xxxz,  \
                             tr_z_xxyzz_xxzz,  \
                             tr_z_xxyzz_xzzz,  \
                             tr_z_xxzz_xxxx,   \
                             tr_z_xxzz_xxxz,   \
                             tr_z_xxzz_xxzz,   \
                             tr_z_xxzz_xzzz,   \
                             tr_z_xyyzz_xxxy,  \
                             tr_z_xyyzz_xxy,   \
                             tr_z_xyyzz_xxyy,  \
                             tr_z_xyyzz_xxyz,  \
                             tr_z_xyyzz_xyy,   \
                             tr_z_xyyzz_xyyy,  \
                             tr_z_xyyzz_xyyz,  \
                             tr_z_xyyzz_xyz,   \
                             tr_z_xyyzz_xyzz,  \
                             tr_z_xyyzz_yyy,   \
                             tr_z_xyyzz_yyyy,  \
                             tr_z_xyyzz_yyyz,  \
                             tr_z_xyyzz_yyz,   \
                             tr_z_xyyzz_yyzz,  \
                             tr_z_xyyzz_yzz,   \
                             tr_z_xyyzz_yzzz,  \
                             tr_z_xyyzz_zzzz,  \
                             tr_z_yyzz_xxxy,   \
                             tr_z_yyzz_xxyy,   \
                             tr_z_yyzz_xxyz,   \
                             tr_z_yyzz_xyyy,   \
                             tr_z_yyzz_xyyz,   \
                             tr_z_yyzz_xyzz,   \
                             tr_z_yyzz_yyyy,   \
                             tr_z_yyzz_yyyz,   \
                             tr_z_yyzz_yyzz,   \
                             tr_z_yyzz_yzzz,   \
                             tr_z_yyzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyzz_xxxx[i] = tr_z_xxzz_xxxx[i] * fe_0 + tr_z_xxyzz_xxxx[i] * pa_y[i];

        tr_z_xxyyzz_xxxy[i] = tr_z_yyzz_xxxy[i] * fe_0 + 3.0 * tr_z_xyyzz_xxy[i] * fe_0 + tr_z_xyyzz_xxxy[i] * pa_x[i];

        tr_z_xxyyzz_xxxz[i] = tr_z_xxzz_xxxz[i] * fe_0 + tr_z_xxyzz_xxxz[i] * pa_y[i];

        tr_z_xxyyzz_xxyy[i] = tr_z_yyzz_xxyy[i] * fe_0 + 2.0 * tr_z_xyyzz_xyy[i] * fe_0 + tr_z_xyyzz_xxyy[i] * pa_x[i];

        tr_z_xxyyzz_xxyz[i] = tr_z_yyzz_xxyz[i] * fe_0 + 2.0 * tr_z_xyyzz_xyz[i] * fe_0 + tr_z_xyyzz_xxyz[i] * pa_x[i];

        tr_z_xxyyzz_xxzz[i] = tr_z_xxzz_xxzz[i] * fe_0 + tr_z_xxyzz_xxzz[i] * pa_y[i];

        tr_z_xxyyzz_xyyy[i] = tr_z_yyzz_xyyy[i] * fe_0 + tr_z_xyyzz_yyy[i] * fe_0 + tr_z_xyyzz_xyyy[i] * pa_x[i];

        tr_z_xxyyzz_xyyz[i] = tr_z_yyzz_xyyz[i] * fe_0 + tr_z_xyyzz_yyz[i] * fe_0 + tr_z_xyyzz_xyyz[i] * pa_x[i];

        tr_z_xxyyzz_xyzz[i] = tr_z_yyzz_xyzz[i] * fe_0 + tr_z_xyyzz_yzz[i] * fe_0 + tr_z_xyyzz_xyzz[i] * pa_x[i];

        tr_z_xxyyzz_xzzz[i] = tr_z_xxzz_xzzz[i] * fe_0 + tr_z_xxyzz_xzzz[i] * pa_y[i];

        tr_z_xxyyzz_yyyy[i] = tr_z_yyzz_yyyy[i] * fe_0 + tr_z_xyyzz_yyyy[i] * pa_x[i];

        tr_z_xxyyzz_yyyz[i] = tr_z_yyzz_yyyz[i] * fe_0 + tr_z_xyyzz_yyyz[i] * pa_x[i];

        tr_z_xxyyzz_yyzz[i] = tr_z_yyzz_yyzz[i] * fe_0 + tr_z_xyyzz_yyzz[i] * pa_x[i];

        tr_z_xxyyzz_yzzz[i] = tr_z_yyzz_yzzz[i] * fe_0 + tr_z_xyyzz_yzzz[i] * pa_x[i];

        tr_z_xxyyzz_zzzz[i] = tr_z_yyzz_zzzz[i] * fe_0 + tr_z_xyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 1035-1050 components of targeted buffer : IG

    auto tr_z_xxyzzz_xxxx = pbuffer.data(idx_dip_ig + 1035);

    auto tr_z_xxyzzz_xxxy = pbuffer.data(idx_dip_ig + 1036);

    auto tr_z_xxyzzz_xxxz = pbuffer.data(idx_dip_ig + 1037);

    auto tr_z_xxyzzz_xxyy = pbuffer.data(idx_dip_ig + 1038);

    auto tr_z_xxyzzz_xxyz = pbuffer.data(idx_dip_ig + 1039);

    auto tr_z_xxyzzz_xxzz = pbuffer.data(idx_dip_ig + 1040);

    auto tr_z_xxyzzz_xyyy = pbuffer.data(idx_dip_ig + 1041);

    auto tr_z_xxyzzz_xyyz = pbuffer.data(idx_dip_ig + 1042);

    auto tr_z_xxyzzz_xyzz = pbuffer.data(idx_dip_ig + 1043);

    auto tr_z_xxyzzz_xzzz = pbuffer.data(idx_dip_ig + 1044);

    auto tr_z_xxyzzz_yyyy = pbuffer.data(idx_dip_ig + 1045);

    auto tr_z_xxyzzz_yyyz = pbuffer.data(idx_dip_ig + 1046);

    auto tr_z_xxyzzz_yyzz = pbuffer.data(idx_dip_ig + 1047);

    auto tr_z_xxyzzz_yzzz = pbuffer.data(idx_dip_ig + 1048);

    auto tr_z_xxyzzz_zzzz = pbuffer.data(idx_dip_ig + 1049);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxyzzz_xxxx, \
                             tr_z_xxyzzz_xxxy, \
                             tr_z_xxyzzz_xxxz, \
                             tr_z_xxyzzz_xxyy, \
                             tr_z_xxyzzz_xxyz, \
                             tr_z_xxyzzz_xxzz, \
                             tr_z_xxyzzz_xyyy, \
                             tr_z_xxyzzz_xyyz, \
                             tr_z_xxyzzz_xyzz, \
                             tr_z_xxyzzz_xzzz, \
                             tr_z_xxyzzz_yyyy, \
                             tr_z_xxyzzz_yyyz, \
                             tr_z_xxyzzz_yyzz, \
                             tr_z_xxyzzz_yzzz, \
                             tr_z_xxyzzz_zzzz, \
                             tr_z_xxzzz_xxx,   \
                             tr_z_xxzzz_xxxx,  \
                             tr_z_xxzzz_xxxy,  \
                             tr_z_xxzzz_xxxz,  \
                             tr_z_xxzzz_xxy,   \
                             tr_z_xxzzz_xxyy,  \
                             tr_z_xxzzz_xxyz,  \
                             tr_z_xxzzz_xxz,   \
                             tr_z_xxzzz_xxzz,  \
                             tr_z_xxzzz_xyy,   \
                             tr_z_xxzzz_xyyy,  \
                             tr_z_xxzzz_xyyz,  \
                             tr_z_xxzzz_xyz,   \
                             tr_z_xxzzz_xyzz,  \
                             tr_z_xxzzz_xzz,   \
                             tr_z_xxzzz_xzzz,  \
                             tr_z_xxzzz_zzzz,  \
                             tr_z_xyzzz_yyyy,  \
                             tr_z_xyzzz_yyyz,  \
                             tr_z_xyzzz_yyzz,  \
                             tr_z_xyzzz_yzzz,  \
                             tr_z_yzzz_yyyy,   \
                             tr_z_yzzz_yyyz,   \
                             tr_z_yzzz_yyzz,   \
                             tr_z_yzzz_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzzz_xxxx[i] = tr_z_xxzzz_xxxx[i] * pa_y[i];

        tr_z_xxyzzz_xxxy[i] = tr_z_xxzzz_xxx[i] * fe_0 + tr_z_xxzzz_xxxy[i] * pa_y[i];

        tr_z_xxyzzz_xxxz[i] = tr_z_xxzzz_xxxz[i] * pa_y[i];

        tr_z_xxyzzz_xxyy[i] = 2.0 * tr_z_xxzzz_xxy[i] * fe_0 + tr_z_xxzzz_xxyy[i] * pa_y[i];

        tr_z_xxyzzz_xxyz[i] = tr_z_xxzzz_xxz[i] * fe_0 + tr_z_xxzzz_xxyz[i] * pa_y[i];

        tr_z_xxyzzz_xxzz[i] = tr_z_xxzzz_xxzz[i] * pa_y[i];

        tr_z_xxyzzz_xyyy[i] = 3.0 * tr_z_xxzzz_xyy[i] * fe_0 + tr_z_xxzzz_xyyy[i] * pa_y[i];

        tr_z_xxyzzz_xyyz[i] = 2.0 * tr_z_xxzzz_xyz[i] * fe_0 + tr_z_xxzzz_xyyz[i] * pa_y[i];

        tr_z_xxyzzz_xyzz[i] = tr_z_xxzzz_xzz[i] * fe_0 + tr_z_xxzzz_xyzz[i] * pa_y[i];

        tr_z_xxyzzz_xzzz[i] = tr_z_xxzzz_xzzz[i] * pa_y[i];

        tr_z_xxyzzz_yyyy[i] = tr_z_yzzz_yyyy[i] * fe_0 + tr_z_xyzzz_yyyy[i] * pa_x[i];

        tr_z_xxyzzz_yyyz[i] = tr_z_yzzz_yyyz[i] * fe_0 + tr_z_xyzzz_yyyz[i] * pa_x[i];

        tr_z_xxyzzz_yyzz[i] = tr_z_yzzz_yyzz[i] * fe_0 + tr_z_xyzzz_yyzz[i] * pa_x[i];

        tr_z_xxyzzz_yzzz[i] = tr_z_yzzz_yzzz[i] * fe_0 + tr_z_xyzzz_yzzz[i] * pa_x[i];

        tr_z_xxyzzz_zzzz[i] = tr_z_xxzzz_zzzz[i] * pa_y[i];
    }

    // Set up 1050-1065 components of targeted buffer : IG

    auto tr_z_xxzzzz_xxxx = pbuffer.data(idx_dip_ig + 1050);

    auto tr_z_xxzzzz_xxxy = pbuffer.data(idx_dip_ig + 1051);

    auto tr_z_xxzzzz_xxxz = pbuffer.data(idx_dip_ig + 1052);

    auto tr_z_xxzzzz_xxyy = pbuffer.data(idx_dip_ig + 1053);

    auto tr_z_xxzzzz_xxyz = pbuffer.data(idx_dip_ig + 1054);

    auto tr_z_xxzzzz_xxzz = pbuffer.data(idx_dip_ig + 1055);

    auto tr_z_xxzzzz_xyyy = pbuffer.data(idx_dip_ig + 1056);

    auto tr_z_xxzzzz_xyyz = pbuffer.data(idx_dip_ig + 1057);

    auto tr_z_xxzzzz_xyzz = pbuffer.data(idx_dip_ig + 1058);

    auto tr_z_xxzzzz_xzzz = pbuffer.data(idx_dip_ig + 1059);

    auto tr_z_xxzzzz_yyyy = pbuffer.data(idx_dip_ig + 1060);

    auto tr_z_xxzzzz_yyyz = pbuffer.data(idx_dip_ig + 1061);

    auto tr_z_xxzzzz_yyzz = pbuffer.data(idx_dip_ig + 1062);

    auto tr_z_xxzzzz_yzzz = pbuffer.data(idx_dip_ig + 1063);

    auto tr_z_xxzzzz_zzzz = pbuffer.data(idx_dip_ig + 1064);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxzzzz_xxxx, \
                             tr_z_xxzzzz_xxxy, \
                             tr_z_xxzzzz_xxxz, \
                             tr_z_xxzzzz_xxyy, \
                             tr_z_xxzzzz_xxyz, \
                             tr_z_xxzzzz_xxzz, \
                             tr_z_xxzzzz_xyyy, \
                             tr_z_xxzzzz_xyyz, \
                             tr_z_xxzzzz_xyzz, \
                             tr_z_xxzzzz_xzzz, \
                             tr_z_xxzzzz_yyyy, \
                             tr_z_xxzzzz_yyyz, \
                             tr_z_xxzzzz_yyzz, \
                             tr_z_xxzzzz_yzzz, \
                             tr_z_xxzzzz_zzzz, \
                             tr_z_xzzzz_xxx,   \
                             tr_z_xzzzz_xxxx,  \
                             tr_z_xzzzz_xxxy,  \
                             tr_z_xzzzz_xxxz,  \
                             tr_z_xzzzz_xxy,   \
                             tr_z_xzzzz_xxyy,  \
                             tr_z_xzzzz_xxyz,  \
                             tr_z_xzzzz_xxz,   \
                             tr_z_xzzzz_xxzz,  \
                             tr_z_xzzzz_xyy,   \
                             tr_z_xzzzz_xyyy,  \
                             tr_z_xzzzz_xyyz,  \
                             tr_z_xzzzz_xyz,   \
                             tr_z_xzzzz_xyzz,  \
                             tr_z_xzzzz_xzz,   \
                             tr_z_xzzzz_xzzz,  \
                             tr_z_xzzzz_yyy,   \
                             tr_z_xzzzz_yyyy,  \
                             tr_z_xzzzz_yyyz,  \
                             tr_z_xzzzz_yyz,   \
                             tr_z_xzzzz_yyzz,  \
                             tr_z_xzzzz_yzz,   \
                             tr_z_xzzzz_yzzz,  \
                             tr_z_xzzzz_zzz,   \
                             tr_z_xzzzz_zzzz,  \
                             tr_z_zzzz_xxxx,   \
                             tr_z_zzzz_xxxy,   \
                             tr_z_zzzz_xxxz,   \
                             tr_z_zzzz_xxyy,   \
                             tr_z_zzzz_xxyz,   \
                             tr_z_zzzz_xxzz,   \
                             tr_z_zzzz_xyyy,   \
                             tr_z_zzzz_xyyz,   \
                             tr_z_zzzz_xyzz,   \
                             tr_z_zzzz_xzzz,   \
                             tr_z_zzzz_yyyy,   \
                             tr_z_zzzz_yyyz,   \
                             tr_z_zzzz_yyzz,   \
                             tr_z_zzzz_yzzz,   \
                             tr_z_zzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzzz_xxxx[i] = tr_z_zzzz_xxxx[i] * fe_0 + 4.0 * tr_z_xzzzz_xxx[i] * fe_0 + tr_z_xzzzz_xxxx[i] * pa_x[i];

        tr_z_xxzzzz_xxxy[i] = tr_z_zzzz_xxxy[i] * fe_0 + 3.0 * tr_z_xzzzz_xxy[i] * fe_0 + tr_z_xzzzz_xxxy[i] * pa_x[i];

        tr_z_xxzzzz_xxxz[i] = tr_z_zzzz_xxxz[i] * fe_0 + 3.0 * tr_z_xzzzz_xxz[i] * fe_0 + tr_z_xzzzz_xxxz[i] * pa_x[i];

        tr_z_xxzzzz_xxyy[i] = tr_z_zzzz_xxyy[i] * fe_0 + 2.0 * tr_z_xzzzz_xyy[i] * fe_0 + tr_z_xzzzz_xxyy[i] * pa_x[i];

        tr_z_xxzzzz_xxyz[i] = tr_z_zzzz_xxyz[i] * fe_0 + 2.0 * tr_z_xzzzz_xyz[i] * fe_0 + tr_z_xzzzz_xxyz[i] * pa_x[i];

        tr_z_xxzzzz_xxzz[i] = tr_z_zzzz_xxzz[i] * fe_0 + 2.0 * tr_z_xzzzz_xzz[i] * fe_0 + tr_z_xzzzz_xxzz[i] * pa_x[i];

        tr_z_xxzzzz_xyyy[i] = tr_z_zzzz_xyyy[i] * fe_0 + tr_z_xzzzz_yyy[i] * fe_0 + tr_z_xzzzz_xyyy[i] * pa_x[i];

        tr_z_xxzzzz_xyyz[i] = tr_z_zzzz_xyyz[i] * fe_0 + tr_z_xzzzz_yyz[i] * fe_0 + tr_z_xzzzz_xyyz[i] * pa_x[i];

        tr_z_xxzzzz_xyzz[i] = tr_z_zzzz_xyzz[i] * fe_0 + tr_z_xzzzz_yzz[i] * fe_0 + tr_z_xzzzz_xyzz[i] * pa_x[i];

        tr_z_xxzzzz_xzzz[i] = tr_z_zzzz_xzzz[i] * fe_0 + tr_z_xzzzz_zzz[i] * fe_0 + tr_z_xzzzz_xzzz[i] * pa_x[i];

        tr_z_xxzzzz_yyyy[i] = tr_z_zzzz_yyyy[i] * fe_0 + tr_z_xzzzz_yyyy[i] * pa_x[i];

        tr_z_xxzzzz_yyyz[i] = tr_z_zzzz_yyyz[i] * fe_0 + tr_z_xzzzz_yyyz[i] * pa_x[i];

        tr_z_xxzzzz_yyzz[i] = tr_z_zzzz_yyzz[i] * fe_0 + tr_z_xzzzz_yyzz[i] * pa_x[i];

        tr_z_xxzzzz_yzzz[i] = tr_z_zzzz_yzzz[i] * fe_0 + tr_z_xzzzz_yzzz[i] * pa_x[i];

        tr_z_xxzzzz_zzzz[i] = tr_z_zzzz_zzzz[i] * fe_0 + tr_z_xzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 1065-1080 components of targeted buffer : IG

    auto tr_z_xyyyyy_xxxx = pbuffer.data(idx_dip_ig + 1065);

    auto tr_z_xyyyyy_xxxy = pbuffer.data(idx_dip_ig + 1066);

    auto tr_z_xyyyyy_xxxz = pbuffer.data(idx_dip_ig + 1067);

    auto tr_z_xyyyyy_xxyy = pbuffer.data(idx_dip_ig + 1068);

    auto tr_z_xyyyyy_xxyz = pbuffer.data(idx_dip_ig + 1069);

    auto tr_z_xyyyyy_xxzz = pbuffer.data(idx_dip_ig + 1070);

    auto tr_z_xyyyyy_xyyy = pbuffer.data(idx_dip_ig + 1071);

    auto tr_z_xyyyyy_xyyz = pbuffer.data(idx_dip_ig + 1072);

    auto tr_z_xyyyyy_xyzz = pbuffer.data(idx_dip_ig + 1073);

    auto tr_z_xyyyyy_xzzz = pbuffer.data(idx_dip_ig + 1074);

    auto tr_z_xyyyyy_yyyy = pbuffer.data(idx_dip_ig + 1075);

    auto tr_z_xyyyyy_yyyz = pbuffer.data(idx_dip_ig + 1076);

    auto tr_z_xyyyyy_yyzz = pbuffer.data(idx_dip_ig + 1077);

    auto tr_z_xyyyyy_yzzz = pbuffer.data(idx_dip_ig + 1078);

    auto tr_z_xyyyyy_zzzz = pbuffer.data(idx_dip_ig + 1079);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyyyy_xxxx, \
                             tr_z_xyyyyy_xxxy, \
                             tr_z_xyyyyy_xxxz, \
                             tr_z_xyyyyy_xxyy, \
                             tr_z_xyyyyy_xxyz, \
                             tr_z_xyyyyy_xxzz, \
                             tr_z_xyyyyy_xyyy, \
                             tr_z_xyyyyy_xyyz, \
                             tr_z_xyyyyy_xyzz, \
                             tr_z_xyyyyy_xzzz, \
                             tr_z_xyyyyy_yyyy, \
                             tr_z_xyyyyy_yyyz, \
                             tr_z_xyyyyy_yyzz, \
                             tr_z_xyyyyy_yzzz, \
                             tr_z_xyyyyy_zzzz, \
                             tr_z_yyyyy_xxx,   \
                             tr_z_yyyyy_xxxx,  \
                             tr_z_yyyyy_xxxy,  \
                             tr_z_yyyyy_xxxz,  \
                             tr_z_yyyyy_xxy,   \
                             tr_z_yyyyy_xxyy,  \
                             tr_z_yyyyy_xxyz,  \
                             tr_z_yyyyy_xxz,   \
                             tr_z_yyyyy_xxzz,  \
                             tr_z_yyyyy_xyy,   \
                             tr_z_yyyyy_xyyy,  \
                             tr_z_yyyyy_xyyz,  \
                             tr_z_yyyyy_xyz,   \
                             tr_z_yyyyy_xyzz,  \
                             tr_z_yyyyy_xzz,   \
                             tr_z_yyyyy_xzzz,  \
                             tr_z_yyyyy_yyy,   \
                             tr_z_yyyyy_yyyy,  \
                             tr_z_yyyyy_yyyz,  \
                             tr_z_yyyyy_yyz,   \
                             tr_z_yyyyy_yyzz,  \
                             tr_z_yyyyy_yzz,   \
                             tr_z_yyyyy_yzzz,  \
                             tr_z_yyyyy_zzz,   \
                             tr_z_yyyyy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyy_xxxx[i] = 4.0 * tr_z_yyyyy_xxx[i] * fe_0 + tr_z_yyyyy_xxxx[i] * pa_x[i];

        tr_z_xyyyyy_xxxy[i] = 3.0 * tr_z_yyyyy_xxy[i] * fe_0 + tr_z_yyyyy_xxxy[i] * pa_x[i];

        tr_z_xyyyyy_xxxz[i] = 3.0 * tr_z_yyyyy_xxz[i] * fe_0 + tr_z_yyyyy_xxxz[i] * pa_x[i];

        tr_z_xyyyyy_xxyy[i] = 2.0 * tr_z_yyyyy_xyy[i] * fe_0 + tr_z_yyyyy_xxyy[i] * pa_x[i];

        tr_z_xyyyyy_xxyz[i] = 2.0 * tr_z_yyyyy_xyz[i] * fe_0 + tr_z_yyyyy_xxyz[i] * pa_x[i];

        tr_z_xyyyyy_xxzz[i] = 2.0 * tr_z_yyyyy_xzz[i] * fe_0 + tr_z_yyyyy_xxzz[i] * pa_x[i];

        tr_z_xyyyyy_xyyy[i] = tr_z_yyyyy_yyy[i] * fe_0 + tr_z_yyyyy_xyyy[i] * pa_x[i];

        tr_z_xyyyyy_xyyz[i] = tr_z_yyyyy_yyz[i] * fe_0 + tr_z_yyyyy_xyyz[i] * pa_x[i];

        tr_z_xyyyyy_xyzz[i] = tr_z_yyyyy_yzz[i] * fe_0 + tr_z_yyyyy_xyzz[i] * pa_x[i];

        tr_z_xyyyyy_xzzz[i] = tr_z_yyyyy_zzz[i] * fe_0 + tr_z_yyyyy_xzzz[i] * pa_x[i];

        tr_z_xyyyyy_yyyy[i] = tr_z_yyyyy_yyyy[i] * pa_x[i];

        tr_z_xyyyyy_yyyz[i] = tr_z_yyyyy_yyyz[i] * pa_x[i];

        tr_z_xyyyyy_yyzz[i] = tr_z_yyyyy_yyzz[i] * pa_x[i];

        tr_z_xyyyyy_yzzz[i] = tr_z_yyyyy_yzzz[i] * pa_x[i];

        tr_z_xyyyyy_zzzz[i] = tr_z_yyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 1080-1095 components of targeted buffer : IG

    auto tr_z_xyyyyz_xxxx = pbuffer.data(idx_dip_ig + 1080);

    auto tr_z_xyyyyz_xxxy = pbuffer.data(idx_dip_ig + 1081);

    auto tr_z_xyyyyz_xxxz = pbuffer.data(idx_dip_ig + 1082);

    auto tr_z_xyyyyz_xxyy = pbuffer.data(idx_dip_ig + 1083);

    auto tr_z_xyyyyz_xxyz = pbuffer.data(idx_dip_ig + 1084);

    auto tr_z_xyyyyz_xxzz = pbuffer.data(idx_dip_ig + 1085);

    auto tr_z_xyyyyz_xyyy = pbuffer.data(idx_dip_ig + 1086);

    auto tr_z_xyyyyz_xyyz = pbuffer.data(idx_dip_ig + 1087);

    auto tr_z_xyyyyz_xyzz = pbuffer.data(idx_dip_ig + 1088);

    auto tr_z_xyyyyz_xzzz = pbuffer.data(idx_dip_ig + 1089);

    auto tr_z_xyyyyz_yyyy = pbuffer.data(idx_dip_ig + 1090);

    auto tr_z_xyyyyz_yyyz = pbuffer.data(idx_dip_ig + 1091);

    auto tr_z_xyyyyz_yyzz = pbuffer.data(idx_dip_ig + 1092);

    auto tr_z_xyyyyz_yzzz = pbuffer.data(idx_dip_ig + 1093);

    auto tr_z_xyyyyz_zzzz = pbuffer.data(idx_dip_ig + 1094);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyyyz_xxxx, \
                             tr_z_xyyyyz_xxxy, \
                             tr_z_xyyyyz_xxxz, \
                             tr_z_xyyyyz_xxyy, \
                             tr_z_xyyyyz_xxyz, \
                             tr_z_xyyyyz_xxzz, \
                             tr_z_xyyyyz_xyyy, \
                             tr_z_xyyyyz_xyyz, \
                             tr_z_xyyyyz_xyzz, \
                             tr_z_xyyyyz_xzzz, \
                             tr_z_xyyyyz_yyyy, \
                             tr_z_xyyyyz_yyyz, \
                             tr_z_xyyyyz_yyzz, \
                             tr_z_xyyyyz_yzzz, \
                             tr_z_xyyyyz_zzzz, \
                             tr_z_yyyyz_xxx,   \
                             tr_z_yyyyz_xxxx,  \
                             tr_z_yyyyz_xxxy,  \
                             tr_z_yyyyz_xxxz,  \
                             tr_z_yyyyz_xxy,   \
                             tr_z_yyyyz_xxyy,  \
                             tr_z_yyyyz_xxyz,  \
                             tr_z_yyyyz_xxz,   \
                             tr_z_yyyyz_xxzz,  \
                             tr_z_yyyyz_xyy,   \
                             tr_z_yyyyz_xyyy,  \
                             tr_z_yyyyz_xyyz,  \
                             tr_z_yyyyz_xyz,   \
                             tr_z_yyyyz_xyzz,  \
                             tr_z_yyyyz_xzz,   \
                             tr_z_yyyyz_xzzz,  \
                             tr_z_yyyyz_yyy,   \
                             tr_z_yyyyz_yyyy,  \
                             tr_z_yyyyz_yyyz,  \
                             tr_z_yyyyz_yyz,   \
                             tr_z_yyyyz_yyzz,  \
                             tr_z_yyyyz_yzz,   \
                             tr_z_yyyyz_yzzz,  \
                             tr_z_yyyyz_zzz,   \
                             tr_z_yyyyz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyz_xxxx[i] = 4.0 * tr_z_yyyyz_xxx[i] * fe_0 + tr_z_yyyyz_xxxx[i] * pa_x[i];

        tr_z_xyyyyz_xxxy[i] = 3.0 * tr_z_yyyyz_xxy[i] * fe_0 + tr_z_yyyyz_xxxy[i] * pa_x[i];

        tr_z_xyyyyz_xxxz[i] = 3.0 * tr_z_yyyyz_xxz[i] * fe_0 + tr_z_yyyyz_xxxz[i] * pa_x[i];

        tr_z_xyyyyz_xxyy[i] = 2.0 * tr_z_yyyyz_xyy[i] * fe_0 + tr_z_yyyyz_xxyy[i] * pa_x[i];

        tr_z_xyyyyz_xxyz[i] = 2.0 * tr_z_yyyyz_xyz[i] * fe_0 + tr_z_yyyyz_xxyz[i] * pa_x[i];

        tr_z_xyyyyz_xxzz[i] = 2.0 * tr_z_yyyyz_xzz[i] * fe_0 + tr_z_yyyyz_xxzz[i] * pa_x[i];

        tr_z_xyyyyz_xyyy[i] = tr_z_yyyyz_yyy[i] * fe_0 + tr_z_yyyyz_xyyy[i] * pa_x[i];

        tr_z_xyyyyz_xyyz[i] = tr_z_yyyyz_yyz[i] * fe_0 + tr_z_yyyyz_xyyz[i] * pa_x[i];

        tr_z_xyyyyz_xyzz[i] = tr_z_yyyyz_yzz[i] * fe_0 + tr_z_yyyyz_xyzz[i] * pa_x[i];

        tr_z_xyyyyz_xzzz[i] = tr_z_yyyyz_zzz[i] * fe_0 + tr_z_yyyyz_xzzz[i] * pa_x[i];

        tr_z_xyyyyz_yyyy[i] = tr_z_yyyyz_yyyy[i] * pa_x[i];

        tr_z_xyyyyz_yyyz[i] = tr_z_yyyyz_yyyz[i] * pa_x[i];

        tr_z_xyyyyz_yyzz[i] = tr_z_yyyyz_yyzz[i] * pa_x[i];

        tr_z_xyyyyz_yzzz[i] = tr_z_yyyyz_yzzz[i] * pa_x[i];

        tr_z_xyyyyz_zzzz[i] = tr_z_yyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 1095-1110 components of targeted buffer : IG

    auto tr_z_xyyyzz_xxxx = pbuffer.data(idx_dip_ig + 1095);

    auto tr_z_xyyyzz_xxxy = pbuffer.data(idx_dip_ig + 1096);

    auto tr_z_xyyyzz_xxxz = pbuffer.data(idx_dip_ig + 1097);

    auto tr_z_xyyyzz_xxyy = pbuffer.data(idx_dip_ig + 1098);

    auto tr_z_xyyyzz_xxyz = pbuffer.data(idx_dip_ig + 1099);

    auto tr_z_xyyyzz_xxzz = pbuffer.data(idx_dip_ig + 1100);

    auto tr_z_xyyyzz_xyyy = pbuffer.data(idx_dip_ig + 1101);

    auto tr_z_xyyyzz_xyyz = pbuffer.data(idx_dip_ig + 1102);

    auto tr_z_xyyyzz_xyzz = pbuffer.data(idx_dip_ig + 1103);

    auto tr_z_xyyyzz_xzzz = pbuffer.data(idx_dip_ig + 1104);

    auto tr_z_xyyyzz_yyyy = pbuffer.data(idx_dip_ig + 1105);

    auto tr_z_xyyyzz_yyyz = pbuffer.data(idx_dip_ig + 1106);

    auto tr_z_xyyyzz_yyzz = pbuffer.data(idx_dip_ig + 1107);

    auto tr_z_xyyyzz_yzzz = pbuffer.data(idx_dip_ig + 1108);

    auto tr_z_xyyyzz_zzzz = pbuffer.data(idx_dip_ig + 1109);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyyzz_xxxx, \
                             tr_z_xyyyzz_xxxy, \
                             tr_z_xyyyzz_xxxz, \
                             tr_z_xyyyzz_xxyy, \
                             tr_z_xyyyzz_xxyz, \
                             tr_z_xyyyzz_xxzz, \
                             tr_z_xyyyzz_xyyy, \
                             tr_z_xyyyzz_xyyz, \
                             tr_z_xyyyzz_xyzz, \
                             tr_z_xyyyzz_xzzz, \
                             tr_z_xyyyzz_yyyy, \
                             tr_z_xyyyzz_yyyz, \
                             tr_z_xyyyzz_yyzz, \
                             tr_z_xyyyzz_yzzz, \
                             tr_z_xyyyzz_zzzz, \
                             tr_z_yyyzz_xxx,   \
                             tr_z_yyyzz_xxxx,  \
                             tr_z_yyyzz_xxxy,  \
                             tr_z_yyyzz_xxxz,  \
                             tr_z_yyyzz_xxy,   \
                             tr_z_yyyzz_xxyy,  \
                             tr_z_yyyzz_xxyz,  \
                             tr_z_yyyzz_xxz,   \
                             tr_z_yyyzz_xxzz,  \
                             tr_z_yyyzz_xyy,   \
                             tr_z_yyyzz_xyyy,  \
                             tr_z_yyyzz_xyyz,  \
                             tr_z_yyyzz_xyz,   \
                             tr_z_yyyzz_xyzz,  \
                             tr_z_yyyzz_xzz,   \
                             tr_z_yyyzz_xzzz,  \
                             tr_z_yyyzz_yyy,   \
                             tr_z_yyyzz_yyyy,  \
                             tr_z_yyyzz_yyyz,  \
                             tr_z_yyyzz_yyz,   \
                             tr_z_yyyzz_yyzz,  \
                             tr_z_yyyzz_yzz,   \
                             tr_z_yyyzz_yzzz,  \
                             tr_z_yyyzz_zzz,   \
                             tr_z_yyyzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyzz_xxxx[i] = 4.0 * tr_z_yyyzz_xxx[i] * fe_0 + tr_z_yyyzz_xxxx[i] * pa_x[i];

        tr_z_xyyyzz_xxxy[i] = 3.0 * tr_z_yyyzz_xxy[i] * fe_0 + tr_z_yyyzz_xxxy[i] * pa_x[i];

        tr_z_xyyyzz_xxxz[i] = 3.0 * tr_z_yyyzz_xxz[i] * fe_0 + tr_z_yyyzz_xxxz[i] * pa_x[i];

        tr_z_xyyyzz_xxyy[i] = 2.0 * tr_z_yyyzz_xyy[i] * fe_0 + tr_z_yyyzz_xxyy[i] * pa_x[i];

        tr_z_xyyyzz_xxyz[i] = 2.0 * tr_z_yyyzz_xyz[i] * fe_0 + tr_z_yyyzz_xxyz[i] * pa_x[i];

        tr_z_xyyyzz_xxzz[i] = 2.0 * tr_z_yyyzz_xzz[i] * fe_0 + tr_z_yyyzz_xxzz[i] * pa_x[i];

        tr_z_xyyyzz_xyyy[i] = tr_z_yyyzz_yyy[i] * fe_0 + tr_z_yyyzz_xyyy[i] * pa_x[i];

        tr_z_xyyyzz_xyyz[i] = tr_z_yyyzz_yyz[i] * fe_0 + tr_z_yyyzz_xyyz[i] * pa_x[i];

        tr_z_xyyyzz_xyzz[i] = tr_z_yyyzz_yzz[i] * fe_0 + tr_z_yyyzz_xyzz[i] * pa_x[i];

        tr_z_xyyyzz_xzzz[i] = tr_z_yyyzz_zzz[i] * fe_0 + tr_z_yyyzz_xzzz[i] * pa_x[i];

        tr_z_xyyyzz_yyyy[i] = tr_z_yyyzz_yyyy[i] * pa_x[i];

        tr_z_xyyyzz_yyyz[i] = tr_z_yyyzz_yyyz[i] * pa_x[i];

        tr_z_xyyyzz_yyzz[i] = tr_z_yyyzz_yyzz[i] * pa_x[i];

        tr_z_xyyyzz_yzzz[i] = tr_z_yyyzz_yzzz[i] * pa_x[i];

        tr_z_xyyyzz_zzzz[i] = tr_z_yyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 1110-1125 components of targeted buffer : IG

    auto tr_z_xyyzzz_xxxx = pbuffer.data(idx_dip_ig + 1110);

    auto tr_z_xyyzzz_xxxy = pbuffer.data(idx_dip_ig + 1111);

    auto tr_z_xyyzzz_xxxz = pbuffer.data(idx_dip_ig + 1112);

    auto tr_z_xyyzzz_xxyy = pbuffer.data(idx_dip_ig + 1113);

    auto tr_z_xyyzzz_xxyz = pbuffer.data(idx_dip_ig + 1114);

    auto tr_z_xyyzzz_xxzz = pbuffer.data(idx_dip_ig + 1115);

    auto tr_z_xyyzzz_xyyy = pbuffer.data(idx_dip_ig + 1116);

    auto tr_z_xyyzzz_xyyz = pbuffer.data(idx_dip_ig + 1117);

    auto tr_z_xyyzzz_xyzz = pbuffer.data(idx_dip_ig + 1118);

    auto tr_z_xyyzzz_xzzz = pbuffer.data(idx_dip_ig + 1119);

    auto tr_z_xyyzzz_yyyy = pbuffer.data(idx_dip_ig + 1120);

    auto tr_z_xyyzzz_yyyz = pbuffer.data(idx_dip_ig + 1121);

    auto tr_z_xyyzzz_yyzz = pbuffer.data(idx_dip_ig + 1122);

    auto tr_z_xyyzzz_yzzz = pbuffer.data(idx_dip_ig + 1123);

    auto tr_z_xyyzzz_zzzz = pbuffer.data(idx_dip_ig + 1124);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyzzz_xxxx, \
                             tr_z_xyyzzz_xxxy, \
                             tr_z_xyyzzz_xxxz, \
                             tr_z_xyyzzz_xxyy, \
                             tr_z_xyyzzz_xxyz, \
                             tr_z_xyyzzz_xxzz, \
                             tr_z_xyyzzz_xyyy, \
                             tr_z_xyyzzz_xyyz, \
                             tr_z_xyyzzz_xyzz, \
                             tr_z_xyyzzz_xzzz, \
                             tr_z_xyyzzz_yyyy, \
                             tr_z_xyyzzz_yyyz, \
                             tr_z_xyyzzz_yyzz, \
                             tr_z_xyyzzz_yzzz, \
                             tr_z_xyyzzz_zzzz, \
                             tr_z_yyzzz_xxx,   \
                             tr_z_yyzzz_xxxx,  \
                             tr_z_yyzzz_xxxy,  \
                             tr_z_yyzzz_xxxz,  \
                             tr_z_yyzzz_xxy,   \
                             tr_z_yyzzz_xxyy,  \
                             tr_z_yyzzz_xxyz,  \
                             tr_z_yyzzz_xxz,   \
                             tr_z_yyzzz_xxzz,  \
                             tr_z_yyzzz_xyy,   \
                             tr_z_yyzzz_xyyy,  \
                             tr_z_yyzzz_xyyz,  \
                             tr_z_yyzzz_xyz,   \
                             tr_z_yyzzz_xyzz,  \
                             tr_z_yyzzz_xzz,   \
                             tr_z_yyzzz_xzzz,  \
                             tr_z_yyzzz_yyy,   \
                             tr_z_yyzzz_yyyy,  \
                             tr_z_yyzzz_yyyz,  \
                             tr_z_yyzzz_yyz,   \
                             tr_z_yyzzz_yyzz,  \
                             tr_z_yyzzz_yzz,   \
                             tr_z_yyzzz_yzzz,  \
                             tr_z_yyzzz_zzz,   \
                             tr_z_yyzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzzz_xxxx[i] = 4.0 * tr_z_yyzzz_xxx[i] * fe_0 + tr_z_yyzzz_xxxx[i] * pa_x[i];

        tr_z_xyyzzz_xxxy[i] = 3.0 * tr_z_yyzzz_xxy[i] * fe_0 + tr_z_yyzzz_xxxy[i] * pa_x[i];

        tr_z_xyyzzz_xxxz[i] = 3.0 * tr_z_yyzzz_xxz[i] * fe_0 + tr_z_yyzzz_xxxz[i] * pa_x[i];

        tr_z_xyyzzz_xxyy[i] = 2.0 * tr_z_yyzzz_xyy[i] * fe_0 + tr_z_yyzzz_xxyy[i] * pa_x[i];

        tr_z_xyyzzz_xxyz[i] = 2.0 * tr_z_yyzzz_xyz[i] * fe_0 + tr_z_yyzzz_xxyz[i] * pa_x[i];

        tr_z_xyyzzz_xxzz[i] = 2.0 * tr_z_yyzzz_xzz[i] * fe_0 + tr_z_yyzzz_xxzz[i] * pa_x[i];

        tr_z_xyyzzz_xyyy[i] = tr_z_yyzzz_yyy[i] * fe_0 + tr_z_yyzzz_xyyy[i] * pa_x[i];

        tr_z_xyyzzz_xyyz[i] = tr_z_yyzzz_yyz[i] * fe_0 + tr_z_yyzzz_xyyz[i] * pa_x[i];

        tr_z_xyyzzz_xyzz[i] = tr_z_yyzzz_yzz[i] * fe_0 + tr_z_yyzzz_xyzz[i] * pa_x[i];

        tr_z_xyyzzz_xzzz[i] = tr_z_yyzzz_zzz[i] * fe_0 + tr_z_yyzzz_xzzz[i] * pa_x[i];

        tr_z_xyyzzz_yyyy[i] = tr_z_yyzzz_yyyy[i] * pa_x[i];

        tr_z_xyyzzz_yyyz[i] = tr_z_yyzzz_yyyz[i] * pa_x[i];

        tr_z_xyyzzz_yyzz[i] = tr_z_yyzzz_yyzz[i] * pa_x[i];

        tr_z_xyyzzz_yzzz[i] = tr_z_yyzzz_yzzz[i] * pa_x[i];

        tr_z_xyyzzz_zzzz[i] = tr_z_yyzzz_zzzz[i] * pa_x[i];
    }

    // Set up 1125-1140 components of targeted buffer : IG

    auto tr_z_xyzzzz_xxxx = pbuffer.data(idx_dip_ig + 1125);

    auto tr_z_xyzzzz_xxxy = pbuffer.data(idx_dip_ig + 1126);

    auto tr_z_xyzzzz_xxxz = pbuffer.data(idx_dip_ig + 1127);

    auto tr_z_xyzzzz_xxyy = pbuffer.data(idx_dip_ig + 1128);

    auto tr_z_xyzzzz_xxyz = pbuffer.data(idx_dip_ig + 1129);

    auto tr_z_xyzzzz_xxzz = pbuffer.data(idx_dip_ig + 1130);

    auto tr_z_xyzzzz_xyyy = pbuffer.data(idx_dip_ig + 1131);

    auto tr_z_xyzzzz_xyyz = pbuffer.data(idx_dip_ig + 1132);

    auto tr_z_xyzzzz_xyzz = pbuffer.data(idx_dip_ig + 1133);

    auto tr_z_xyzzzz_xzzz = pbuffer.data(idx_dip_ig + 1134);

    auto tr_z_xyzzzz_yyyy = pbuffer.data(idx_dip_ig + 1135);

    auto tr_z_xyzzzz_yyyz = pbuffer.data(idx_dip_ig + 1136);

    auto tr_z_xyzzzz_yyzz = pbuffer.data(idx_dip_ig + 1137);

    auto tr_z_xyzzzz_yzzz = pbuffer.data(idx_dip_ig + 1138);

    auto tr_z_xyzzzz_zzzz = pbuffer.data(idx_dip_ig + 1139);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xyzzzz_xxxx, \
                             tr_z_xyzzzz_xxxy, \
                             tr_z_xyzzzz_xxxz, \
                             tr_z_xyzzzz_xxyy, \
                             tr_z_xyzzzz_xxyz, \
                             tr_z_xyzzzz_xxzz, \
                             tr_z_xyzzzz_xyyy, \
                             tr_z_xyzzzz_xyyz, \
                             tr_z_xyzzzz_xyzz, \
                             tr_z_xyzzzz_xzzz, \
                             tr_z_xyzzzz_yyyy, \
                             tr_z_xyzzzz_yyyz, \
                             tr_z_xyzzzz_yyzz, \
                             tr_z_xyzzzz_yzzz, \
                             tr_z_xyzzzz_zzzz, \
                             tr_z_xzzzz_xxxx,  \
                             tr_z_xzzzz_xxxz,  \
                             tr_z_xzzzz_xxzz,  \
                             tr_z_xzzzz_xzzz,  \
                             tr_z_yzzzz_xxxy,  \
                             tr_z_yzzzz_xxy,   \
                             tr_z_yzzzz_xxyy,  \
                             tr_z_yzzzz_xxyz,  \
                             tr_z_yzzzz_xyy,   \
                             tr_z_yzzzz_xyyy,  \
                             tr_z_yzzzz_xyyz,  \
                             tr_z_yzzzz_xyz,   \
                             tr_z_yzzzz_xyzz,  \
                             tr_z_yzzzz_yyy,   \
                             tr_z_yzzzz_yyyy,  \
                             tr_z_yzzzz_yyyz,  \
                             tr_z_yzzzz_yyz,   \
                             tr_z_yzzzz_yyzz,  \
                             tr_z_yzzzz_yzz,   \
                             tr_z_yzzzz_yzzz,  \
                             tr_z_yzzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzzz_xxxx[i] = tr_z_xzzzz_xxxx[i] * pa_y[i];

        tr_z_xyzzzz_xxxy[i] = 3.0 * tr_z_yzzzz_xxy[i] * fe_0 + tr_z_yzzzz_xxxy[i] * pa_x[i];

        tr_z_xyzzzz_xxxz[i] = tr_z_xzzzz_xxxz[i] * pa_y[i];

        tr_z_xyzzzz_xxyy[i] = 2.0 * tr_z_yzzzz_xyy[i] * fe_0 + tr_z_yzzzz_xxyy[i] * pa_x[i];

        tr_z_xyzzzz_xxyz[i] = 2.0 * tr_z_yzzzz_xyz[i] * fe_0 + tr_z_yzzzz_xxyz[i] * pa_x[i];

        tr_z_xyzzzz_xxzz[i] = tr_z_xzzzz_xxzz[i] * pa_y[i];

        tr_z_xyzzzz_xyyy[i] = tr_z_yzzzz_yyy[i] * fe_0 + tr_z_yzzzz_xyyy[i] * pa_x[i];

        tr_z_xyzzzz_xyyz[i] = tr_z_yzzzz_yyz[i] * fe_0 + tr_z_yzzzz_xyyz[i] * pa_x[i];

        tr_z_xyzzzz_xyzz[i] = tr_z_yzzzz_yzz[i] * fe_0 + tr_z_yzzzz_xyzz[i] * pa_x[i];

        tr_z_xyzzzz_xzzz[i] = tr_z_xzzzz_xzzz[i] * pa_y[i];

        tr_z_xyzzzz_yyyy[i] = tr_z_yzzzz_yyyy[i] * pa_x[i];

        tr_z_xyzzzz_yyyz[i] = tr_z_yzzzz_yyyz[i] * pa_x[i];

        tr_z_xyzzzz_yyzz[i] = tr_z_yzzzz_yyzz[i] * pa_x[i];

        tr_z_xyzzzz_yzzz[i] = tr_z_yzzzz_yzzz[i] * pa_x[i];

        tr_z_xyzzzz_zzzz[i] = tr_z_yzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 1140-1155 components of targeted buffer : IG

    auto tr_z_xzzzzz_xxxx = pbuffer.data(idx_dip_ig + 1140);

    auto tr_z_xzzzzz_xxxy = pbuffer.data(idx_dip_ig + 1141);

    auto tr_z_xzzzzz_xxxz = pbuffer.data(idx_dip_ig + 1142);

    auto tr_z_xzzzzz_xxyy = pbuffer.data(idx_dip_ig + 1143);

    auto tr_z_xzzzzz_xxyz = pbuffer.data(idx_dip_ig + 1144);

    auto tr_z_xzzzzz_xxzz = pbuffer.data(idx_dip_ig + 1145);

    auto tr_z_xzzzzz_xyyy = pbuffer.data(idx_dip_ig + 1146);

    auto tr_z_xzzzzz_xyyz = pbuffer.data(idx_dip_ig + 1147);

    auto tr_z_xzzzzz_xyzz = pbuffer.data(idx_dip_ig + 1148);

    auto tr_z_xzzzzz_xzzz = pbuffer.data(idx_dip_ig + 1149);

    auto tr_z_xzzzzz_yyyy = pbuffer.data(idx_dip_ig + 1150);

    auto tr_z_xzzzzz_yyyz = pbuffer.data(idx_dip_ig + 1151);

    auto tr_z_xzzzzz_yyzz = pbuffer.data(idx_dip_ig + 1152);

    auto tr_z_xzzzzz_yzzz = pbuffer.data(idx_dip_ig + 1153);

    auto tr_z_xzzzzz_zzzz = pbuffer.data(idx_dip_ig + 1154);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xzzzzz_xxxx, \
                             tr_z_xzzzzz_xxxy, \
                             tr_z_xzzzzz_xxxz, \
                             tr_z_xzzzzz_xxyy, \
                             tr_z_xzzzzz_xxyz, \
                             tr_z_xzzzzz_xxzz, \
                             tr_z_xzzzzz_xyyy, \
                             tr_z_xzzzzz_xyyz, \
                             tr_z_xzzzzz_xyzz, \
                             tr_z_xzzzzz_xzzz, \
                             tr_z_xzzzzz_yyyy, \
                             tr_z_xzzzzz_yyyz, \
                             tr_z_xzzzzz_yyzz, \
                             tr_z_xzzzzz_yzzz, \
                             tr_z_xzzzzz_zzzz, \
                             tr_z_zzzzz_xxx,   \
                             tr_z_zzzzz_xxxx,  \
                             tr_z_zzzzz_xxxy,  \
                             tr_z_zzzzz_xxxz,  \
                             tr_z_zzzzz_xxy,   \
                             tr_z_zzzzz_xxyy,  \
                             tr_z_zzzzz_xxyz,  \
                             tr_z_zzzzz_xxz,   \
                             tr_z_zzzzz_xxzz,  \
                             tr_z_zzzzz_xyy,   \
                             tr_z_zzzzz_xyyy,  \
                             tr_z_zzzzz_xyyz,  \
                             tr_z_zzzzz_xyz,   \
                             tr_z_zzzzz_xyzz,  \
                             tr_z_zzzzz_xzz,   \
                             tr_z_zzzzz_xzzz,  \
                             tr_z_zzzzz_yyy,   \
                             tr_z_zzzzz_yyyy,  \
                             tr_z_zzzzz_yyyz,  \
                             tr_z_zzzzz_yyz,   \
                             tr_z_zzzzz_yyzz,  \
                             tr_z_zzzzz_yzz,   \
                             tr_z_zzzzz_yzzz,  \
                             tr_z_zzzzz_zzz,   \
                             tr_z_zzzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzzz_xxxx[i] = 4.0 * tr_z_zzzzz_xxx[i] * fe_0 + tr_z_zzzzz_xxxx[i] * pa_x[i];

        tr_z_xzzzzz_xxxy[i] = 3.0 * tr_z_zzzzz_xxy[i] * fe_0 + tr_z_zzzzz_xxxy[i] * pa_x[i];

        tr_z_xzzzzz_xxxz[i] = 3.0 * tr_z_zzzzz_xxz[i] * fe_0 + tr_z_zzzzz_xxxz[i] * pa_x[i];

        tr_z_xzzzzz_xxyy[i] = 2.0 * tr_z_zzzzz_xyy[i] * fe_0 + tr_z_zzzzz_xxyy[i] * pa_x[i];

        tr_z_xzzzzz_xxyz[i] = 2.0 * tr_z_zzzzz_xyz[i] * fe_0 + tr_z_zzzzz_xxyz[i] * pa_x[i];

        tr_z_xzzzzz_xxzz[i] = 2.0 * tr_z_zzzzz_xzz[i] * fe_0 + tr_z_zzzzz_xxzz[i] * pa_x[i];

        tr_z_xzzzzz_xyyy[i] = tr_z_zzzzz_yyy[i] * fe_0 + tr_z_zzzzz_xyyy[i] * pa_x[i];

        tr_z_xzzzzz_xyyz[i] = tr_z_zzzzz_yyz[i] * fe_0 + tr_z_zzzzz_xyyz[i] * pa_x[i];

        tr_z_xzzzzz_xyzz[i] = tr_z_zzzzz_yzz[i] * fe_0 + tr_z_zzzzz_xyzz[i] * pa_x[i];

        tr_z_xzzzzz_xzzz[i] = tr_z_zzzzz_zzz[i] * fe_0 + tr_z_zzzzz_xzzz[i] * pa_x[i];

        tr_z_xzzzzz_yyyy[i] = tr_z_zzzzz_yyyy[i] * pa_x[i];

        tr_z_xzzzzz_yyyz[i] = tr_z_zzzzz_yyyz[i] * pa_x[i];

        tr_z_xzzzzz_yyzz[i] = tr_z_zzzzz_yyzz[i] * pa_x[i];

        tr_z_xzzzzz_yzzz[i] = tr_z_zzzzz_yzzz[i] * pa_x[i];

        tr_z_xzzzzz_zzzz[i] = tr_z_zzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 1155-1170 components of targeted buffer : IG

    auto tr_z_yyyyyy_xxxx = pbuffer.data(idx_dip_ig + 1155);

    auto tr_z_yyyyyy_xxxy = pbuffer.data(idx_dip_ig + 1156);

    auto tr_z_yyyyyy_xxxz = pbuffer.data(idx_dip_ig + 1157);

    auto tr_z_yyyyyy_xxyy = pbuffer.data(idx_dip_ig + 1158);

    auto tr_z_yyyyyy_xxyz = pbuffer.data(idx_dip_ig + 1159);

    auto tr_z_yyyyyy_xxzz = pbuffer.data(idx_dip_ig + 1160);

    auto tr_z_yyyyyy_xyyy = pbuffer.data(idx_dip_ig + 1161);

    auto tr_z_yyyyyy_xyyz = pbuffer.data(idx_dip_ig + 1162);

    auto tr_z_yyyyyy_xyzz = pbuffer.data(idx_dip_ig + 1163);

    auto tr_z_yyyyyy_xzzz = pbuffer.data(idx_dip_ig + 1164);

    auto tr_z_yyyyyy_yyyy = pbuffer.data(idx_dip_ig + 1165);

    auto tr_z_yyyyyy_yyyz = pbuffer.data(idx_dip_ig + 1166);

    auto tr_z_yyyyyy_yyzz = pbuffer.data(idx_dip_ig + 1167);

    auto tr_z_yyyyyy_yzzz = pbuffer.data(idx_dip_ig + 1168);

    auto tr_z_yyyyyy_zzzz = pbuffer.data(idx_dip_ig + 1169);

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyyy_xxxx,   \
                             tr_z_yyyy_xxxy,   \
                             tr_z_yyyy_xxxz,   \
                             tr_z_yyyy_xxyy,   \
                             tr_z_yyyy_xxyz,   \
                             tr_z_yyyy_xxzz,   \
                             tr_z_yyyy_xyyy,   \
                             tr_z_yyyy_xyyz,   \
                             tr_z_yyyy_xyzz,   \
                             tr_z_yyyy_xzzz,   \
                             tr_z_yyyy_yyyy,   \
                             tr_z_yyyy_yyyz,   \
                             tr_z_yyyy_yyzz,   \
                             tr_z_yyyy_yzzz,   \
                             tr_z_yyyy_zzzz,   \
                             tr_z_yyyyy_xxx,   \
                             tr_z_yyyyy_xxxx,  \
                             tr_z_yyyyy_xxxy,  \
                             tr_z_yyyyy_xxxz,  \
                             tr_z_yyyyy_xxy,   \
                             tr_z_yyyyy_xxyy,  \
                             tr_z_yyyyy_xxyz,  \
                             tr_z_yyyyy_xxz,   \
                             tr_z_yyyyy_xxzz,  \
                             tr_z_yyyyy_xyy,   \
                             tr_z_yyyyy_xyyy,  \
                             tr_z_yyyyy_xyyz,  \
                             tr_z_yyyyy_xyz,   \
                             tr_z_yyyyy_xyzz,  \
                             tr_z_yyyyy_xzz,   \
                             tr_z_yyyyy_xzzz,  \
                             tr_z_yyyyy_yyy,   \
                             tr_z_yyyyy_yyyy,  \
                             tr_z_yyyyy_yyyz,  \
                             tr_z_yyyyy_yyz,   \
                             tr_z_yyyyy_yyzz,  \
                             tr_z_yyyyy_yzz,   \
                             tr_z_yyyyy_yzzz,  \
                             tr_z_yyyyy_zzz,   \
                             tr_z_yyyyy_zzzz,  \
                             tr_z_yyyyyy_xxxx, \
                             tr_z_yyyyyy_xxxy, \
                             tr_z_yyyyyy_xxxz, \
                             tr_z_yyyyyy_xxyy, \
                             tr_z_yyyyyy_xxyz, \
                             tr_z_yyyyyy_xxzz, \
                             tr_z_yyyyyy_xyyy, \
                             tr_z_yyyyyy_xyyz, \
                             tr_z_yyyyyy_xyzz, \
                             tr_z_yyyyyy_xzzz, \
                             tr_z_yyyyyy_yyyy, \
                             tr_z_yyyyyy_yyyz, \
                             tr_z_yyyyyy_yyzz, \
                             tr_z_yyyyyy_yzzz, \
                             tr_z_yyyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyy_xxxx[i] = 5.0 * tr_z_yyyy_xxxx[i] * fe_0 + tr_z_yyyyy_xxxx[i] * pa_y[i];

        tr_z_yyyyyy_xxxy[i] = 5.0 * tr_z_yyyy_xxxy[i] * fe_0 + tr_z_yyyyy_xxx[i] * fe_0 + tr_z_yyyyy_xxxy[i] * pa_y[i];

        tr_z_yyyyyy_xxxz[i] = 5.0 * tr_z_yyyy_xxxz[i] * fe_0 + tr_z_yyyyy_xxxz[i] * pa_y[i];

        tr_z_yyyyyy_xxyy[i] = 5.0 * tr_z_yyyy_xxyy[i] * fe_0 + 2.0 * tr_z_yyyyy_xxy[i] * fe_0 + tr_z_yyyyy_xxyy[i] * pa_y[i];

        tr_z_yyyyyy_xxyz[i] = 5.0 * tr_z_yyyy_xxyz[i] * fe_0 + tr_z_yyyyy_xxz[i] * fe_0 + tr_z_yyyyy_xxyz[i] * pa_y[i];

        tr_z_yyyyyy_xxzz[i] = 5.0 * tr_z_yyyy_xxzz[i] * fe_0 + tr_z_yyyyy_xxzz[i] * pa_y[i];

        tr_z_yyyyyy_xyyy[i] = 5.0 * tr_z_yyyy_xyyy[i] * fe_0 + 3.0 * tr_z_yyyyy_xyy[i] * fe_0 + tr_z_yyyyy_xyyy[i] * pa_y[i];

        tr_z_yyyyyy_xyyz[i] = 5.0 * tr_z_yyyy_xyyz[i] * fe_0 + 2.0 * tr_z_yyyyy_xyz[i] * fe_0 + tr_z_yyyyy_xyyz[i] * pa_y[i];

        tr_z_yyyyyy_xyzz[i] = 5.0 * tr_z_yyyy_xyzz[i] * fe_0 + tr_z_yyyyy_xzz[i] * fe_0 + tr_z_yyyyy_xyzz[i] * pa_y[i];

        tr_z_yyyyyy_xzzz[i] = 5.0 * tr_z_yyyy_xzzz[i] * fe_0 + tr_z_yyyyy_xzzz[i] * pa_y[i];

        tr_z_yyyyyy_yyyy[i] = 5.0 * tr_z_yyyy_yyyy[i] * fe_0 + 4.0 * tr_z_yyyyy_yyy[i] * fe_0 + tr_z_yyyyy_yyyy[i] * pa_y[i];

        tr_z_yyyyyy_yyyz[i] = 5.0 * tr_z_yyyy_yyyz[i] * fe_0 + 3.0 * tr_z_yyyyy_yyz[i] * fe_0 + tr_z_yyyyy_yyyz[i] * pa_y[i];

        tr_z_yyyyyy_yyzz[i] = 5.0 * tr_z_yyyy_yyzz[i] * fe_0 + 2.0 * tr_z_yyyyy_yzz[i] * fe_0 + tr_z_yyyyy_yyzz[i] * pa_y[i];

        tr_z_yyyyyy_yzzz[i] = 5.0 * tr_z_yyyy_yzzz[i] * fe_0 + tr_z_yyyyy_zzz[i] * fe_0 + tr_z_yyyyy_yzzz[i] * pa_y[i];

        tr_z_yyyyyy_zzzz[i] = 5.0 * tr_z_yyyy_zzzz[i] * fe_0 + tr_z_yyyyy_zzzz[i] * pa_y[i];
    }

    // Set up 1170-1185 components of targeted buffer : IG

    auto tr_z_yyyyyz_xxxx = pbuffer.data(idx_dip_ig + 1170);

    auto tr_z_yyyyyz_xxxy = pbuffer.data(idx_dip_ig + 1171);

    auto tr_z_yyyyyz_xxxz = pbuffer.data(idx_dip_ig + 1172);

    auto tr_z_yyyyyz_xxyy = pbuffer.data(idx_dip_ig + 1173);

    auto tr_z_yyyyyz_xxyz = pbuffer.data(idx_dip_ig + 1174);

    auto tr_z_yyyyyz_xxzz = pbuffer.data(idx_dip_ig + 1175);

    auto tr_z_yyyyyz_xyyy = pbuffer.data(idx_dip_ig + 1176);

    auto tr_z_yyyyyz_xyyz = pbuffer.data(idx_dip_ig + 1177);

    auto tr_z_yyyyyz_xyzz = pbuffer.data(idx_dip_ig + 1178);

    auto tr_z_yyyyyz_xzzz = pbuffer.data(idx_dip_ig + 1179);

    auto tr_z_yyyyyz_yyyy = pbuffer.data(idx_dip_ig + 1180);

    auto tr_z_yyyyyz_yyyz = pbuffer.data(idx_dip_ig + 1181);

    auto tr_z_yyyyyz_yyzz = pbuffer.data(idx_dip_ig + 1182);

    auto tr_z_yyyyyz_yzzz = pbuffer.data(idx_dip_ig + 1183);

    auto tr_z_yyyyyz_zzzz = pbuffer.data(idx_dip_ig + 1184);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_z_yyyyy_xxxy,  \
                             tr_z_yyyyy_xxyy,  \
                             tr_z_yyyyy_xyyy,  \
                             tr_z_yyyyy_yyyy,  \
                             tr_z_yyyyyz_xxxx, \
                             tr_z_yyyyyz_xxxy, \
                             tr_z_yyyyyz_xxxz, \
                             tr_z_yyyyyz_xxyy, \
                             tr_z_yyyyyz_xxyz, \
                             tr_z_yyyyyz_xxzz, \
                             tr_z_yyyyyz_xyyy, \
                             tr_z_yyyyyz_xyyz, \
                             tr_z_yyyyyz_xyzz, \
                             tr_z_yyyyyz_xzzz, \
                             tr_z_yyyyyz_yyyy, \
                             tr_z_yyyyyz_yyyz, \
                             tr_z_yyyyyz_yyzz, \
                             tr_z_yyyyyz_yzzz, \
                             tr_z_yyyyyz_zzzz, \
                             tr_z_yyyyz_xxxx,  \
                             tr_z_yyyyz_xxxz,  \
                             tr_z_yyyyz_xxyz,  \
                             tr_z_yyyyz_xxz,   \
                             tr_z_yyyyz_xxzz,  \
                             tr_z_yyyyz_xyyz,  \
                             tr_z_yyyyz_xyz,   \
                             tr_z_yyyyz_xyzz,  \
                             tr_z_yyyyz_xzz,   \
                             tr_z_yyyyz_xzzz,  \
                             tr_z_yyyyz_yyyz,  \
                             tr_z_yyyyz_yyz,   \
                             tr_z_yyyyz_yyzz,  \
                             tr_z_yyyyz_yzz,   \
                             tr_z_yyyyz_yzzz,  \
                             tr_z_yyyyz_zzz,   \
                             tr_z_yyyyz_zzzz,  \
                             tr_z_yyyz_xxxx,   \
                             tr_z_yyyz_xxxz,   \
                             tr_z_yyyz_xxyz,   \
                             tr_z_yyyz_xxzz,   \
                             tr_z_yyyz_xyyz,   \
                             tr_z_yyyz_xyzz,   \
                             tr_z_yyyz_xzzz,   \
                             tr_z_yyyz_yyyz,   \
                             tr_z_yyyz_yyzz,   \
                             tr_z_yyyz_yzzz,   \
                             tr_z_yyyz_zzzz,   \
                             ts_yyyyy_xxxy,    \
                             ts_yyyyy_xxyy,    \
                             ts_yyyyy_xyyy,    \
                             ts_yyyyy_yyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyz_xxxx[i] = 4.0 * tr_z_yyyz_xxxx[i] * fe_0 + tr_z_yyyyz_xxxx[i] * pa_y[i];

        tr_z_yyyyyz_xxxy[i] = ts_yyyyy_xxxy[i] * fe_0 + tr_z_yyyyy_xxxy[i] * pa_z[i];

        tr_z_yyyyyz_xxxz[i] = 4.0 * tr_z_yyyz_xxxz[i] * fe_0 + tr_z_yyyyz_xxxz[i] * pa_y[i];

        tr_z_yyyyyz_xxyy[i] = ts_yyyyy_xxyy[i] * fe_0 + tr_z_yyyyy_xxyy[i] * pa_z[i];

        tr_z_yyyyyz_xxyz[i] = 4.0 * tr_z_yyyz_xxyz[i] * fe_0 + tr_z_yyyyz_xxz[i] * fe_0 + tr_z_yyyyz_xxyz[i] * pa_y[i];

        tr_z_yyyyyz_xxzz[i] = 4.0 * tr_z_yyyz_xxzz[i] * fe_0 + tr_z_yyyyz_xxzz[i] * pa_y[i];

        tr_z_yyyyyz_xyyy[i] = ts_yyyyy_xyyy[i] * fe_0 + tr_z_yyyyy_xyyy[i] * pa_z[i];

        tr_z_yyyyyz_xyyz[i] = 4.0 * tr_z_yyyz_xyyz[i] * fe_0 + 2.0 * tr_z_yyyyz_xyz[i] * fe_0 + tr_z_yyyyz_xyyz[i] * pa_y[i];

        tr_z_yyyyyz_xyzz[i] = 4.0 * tr_z_yyyz_xyzz[i] * fe_0 + tr_z_yyyyz_xzz[i] * fe_0 + tr_z_yyyyz_xyzz[i] * pa_y[i];

        tr_z_yyyyyz_xzzz[i] = 4.0 * tr_z_yyyz_xzzz[i] * fe_0 + tr_z_yyyyz_xzzz[i] * pa_y[i];

        tr_z_yyyyyz_yyyy[i] = ts_yyyyy_yyyy[i] * fe_0 + tr_z_yyyyy_yyyy[i] * pa_z[i];

        tr_z_yyyyyz_yyyz[i] = 4.0 * tr_z_yyyz_yyyz[i] * fe_0 + 3.0 * tr_z_yyyyz_yyz[i] * fe_0 + tr_z_yyyyz_yyyz[i] * pa_y[i];

        tr_z_yyyyyz_yyzz[i] = 4.0 * tr_z_yyyz_yyzz[i] * fe_0 + 2.0 * tr_z_yyyyz_yzz[i] * fe_0 + tr_z_yyyyz_yyzz[i] * pa_y[i];

        tr_z_yyyyyz_yzzz[i] = 4.0 * tr_z_yyyz_yzzz[i] * fe_0 + tr_z_yyyyz_zzz[i] * fe_0 + tr_z_yyyyz_yzzz[i] * pa_y[i];

        tr_z_yyyyyz_zzzz[i] = 4.0 * tr_z_yyyz_zzzz[i] * fe_0 + tr_z_yyyyz_zzzz[i] * pa_y[i];
    }

    // Set up 1185-1200 components of targeted buffer : IG

    auto tr_z_yyyyzz_xxxx = pbuffer.data(idx_dip_ig + 1185);

    auto tr_z_yyyyzz_xxxy = pbuffer.data(idx_dip_ig + 1186);

    auto tr_z_yyyyzz_xxxz = pbuffer.data(idx_dip_ig + 1187);

    auto tr_z_yyyyzz_xxyy = pbuffer.data(idx_dip_ig + 1188);

    auto tr_z_yyyyzz_xxyz = pbuffer.data(idx_dip_ig + 1189);

    auto tr_z_yyyyzz_xxzz = pbuffer.data(idx_dip_ig + 1190);

    auto tr_z_yyyyzz_xyyy = pbuffer.data(idx_dip_ig + 1191);

    auto tr_z_yyyyzz_xyyz = pbuffer.data(idx_dip_ig + 1192);

    auto tr_z_yyyyzz_xyzz = pbuffer.data(idx_dip_ig + 1193);

    auto tr_z_yyyyzz_xzzz = pbuffer.data(idx_dip_ig + 1194);

    auto tr_z_yyyyzz_yyyy = pbuffer.data(idx_dip_ig + 1195);

    auto tr_z_yyyyzz_yyyz = pbuffer.data(idx_dip_ig + 1196);

    auto tr_z_yyyyzz_yyzz = pbuffer.data(idx_dip_ig + 1197);

    auto tr_z_yyyyzz_yzzz = pbuffer.data(idx_dip_ig + 1198);

    auto tr_z_yyyyzz_zzzz = pbuffer.data(idx_dip_ig + 1199);

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyyyzz_xxxx, \
                             tr_z_yyyyzz_xxxy, \
                             tr_z_yyyyzz_xxxz, \
                             tr_z_yyyyzz_xxyy, \
                             tr_z_yyyyzz_xxyz, \
                             tr_z_yyyyzz_xxzz, \
                             tr_z_yyyyzz_xyyy, \
                             tr_z_yyyyzz_xyyz, \
                             tr_z_yyyyzz_xyzz, \
                             tr_z_yyyyzz_xzzz, \
                             tr_z_yyyyzz_yyyy, \
                             tr_z_yyyyzz_yyyz, \
                             tr_z_yyyyzz_yyzz, \
                             tr_z_yyyyzz_yzzz, \
                             tr_z_yyyyzz_zzzz, \
                             tr_z_yyyzz_xxx,   \
                             tr_z_yyyzz_xxxx,  \
                             tr_z_yyyzz_xxxy,  \
                             tr_z_yyyzz_xxxz,  \
                             tr_z_yyyzz_xxy,   \
                             tr_z_yyyzz_xxyy,  \
                             tr_z_yyyzz_xxyz,  \
                             tr_z_yyyzz_xxz,   \
                             tr_z_yyyzz_xxzz,  \
                             tr_z_yyyzz_xyy,   \
                             tr_z_yyyzz_xyyy,  \
                             tr_z_yyyzz_xyyz,  \
                             tr_z_yyyzz_xyz,   \
                             tr_z_yyyzz_xyzz,  \
                             tr_z_yyyzz_xzz,   \
                             tr_z_yyyzz_xzzz,  \
                             tr_z_yyyzz_yyy,   \
                             tr_z_yyyzz_yyyy,  \
                             tr_z_yyyzz_yyyz,  \
                             tr_z_yyyzz_yyz,   \
                             tr_z_yyyzz_yyzz,  \
                             tr_z_yyyzz_yzz,   \
                             tr_z_yyyzz_yzzz,  \
                             tr_z_yyyzz_zzz,   \
                             tr_z_yyyzz_zzzz,  \
                             tr_z_yyzz_xxxx,   \
                             tr_z_yyzz_xxxy,   \
                             tr_z_yyzz_xxxz,   \
                             tr_z_yyzz_xxyy,   \
                             tr_z_yyzz_xxyz,   \
                             tr_z_yyzz_xxzz,   \
                             tr_z_yyzz_xyyy,   \
                             tr_z_yyzz_xyyz,   \
                             tr_z_yyzz_xyzz,   \
                             tr_z_yyzz_xzzz,   \
                             tr_z_yyzz_yyyy,   \
                             tr_z_yyzz_yyyz,   \
                             tr_z_yyzz_yyzz,   \
                             tr_z_yyzz_yzzz,   \
                             tr_z_yyzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyzz_xxxx[i] = 3.0 * tr_z_yyzz_xxxx[i] * fe_0 + tr_z_yyyzz_xxxx[i] * pa_y[i];

        tr_z_yyyyzz_xxxy[i] = 3.0 * tr_z_yyzz_xxxy[i] * fe_0 + tr_z_yyyzz_xxx[i] * fe_0 + tr_z_yyyzz_xxxy[i] * pa_y[i];

        tr_z_yyyyzz_xxxz[i] = 3.0 * tr_z_yyzz_xxxz[i] * fe_0 + tr_z_yyyzz_xxxz[i] * pa_y[i];

        tr_z_yyyyzz_xxyy[i] = 3.0 * tr_z_yyzz_xxyy[i] * fe_0 + 2.0 * tr_z_yyyzz_xxy[i] * fe_0 + tr_z_yyyzz_xxyy[i] * pa_y[i];

        tr_z_yyyyzz_xxyz[i] = 3.0 * tr_z_yyzz_xxyz[i] * fe_0 + tr_z_yyyzz_xxz[i] * fe_0 + tr_z_yyyzz_xxyz[i] * pa_y[i];

        tr_z_yyyyzz_xxzz[i] = 3.0 * tr_z_yyzz_xxzz[i] * fe_0 + tr_z_yyyzz_xxzz[i] * pa_y[i];

        tr_z_yyyyzz_xyyy[i] = 3.0 * tr_z_yyzz_xyyy[i] * fe_0 + 3.0 * tr_z_yyyzz_xyy[i] * fe_0 + tr_z_yyyzz_xyyy[i] * pa_y[i];

        tr_z_yyyyzz_xyyz[i] = 3.0 * tr_z_yyzz_xyyz[i] * fe_0 + 2.0 * tr_z_yyyzz_xyz[i] * fe_0 + tr_z_yyyzz_xyyz[i] * pa_y[i];

        tr_z_yyyyzz_xyzz[i] = 3.0 * tr_z_yyzz_xyzz[i] * fe_0 + tr_z_yyyzz_xzz[i] * fe_0 + tr_z_yyyzz_xyzz[i] * pa_y[i];

        tr_z_yyyyzz_xzzz[i] = 3.0 * tr_z_yyzz_xzzz[i] * fe_0 + tr_z_yyyzz_xzzz[i] * pa_y[i];

        tr_z_yyyyzz_yyyy[i] = 3.0 * tr_z_yyzz_yyyy[i] * fe_0 + 4.0 * tr_z_yyyzz_yyy[i] * fe_0 + tr_z_yyyzz_yyyy[i] * pa_y[i];

        tr_z_yyyyzz_yyyz[i] = 3.0 * tr_z_yyzz_yyyz[i] * fe_0 + 3.0 * tr_z_yyyzz_yyz[i] * fe_0 + tr_z_yyyzz_yyyz[i] * pa_y[i];

        tr_z_yyyyzz_yyzz[i] = 3.0 * tr_z_yyzz_yyzz[i] * fe_0 + 2.0 * tr_z_yyyzz_yzz[i] * fe_0 + tr_z_yyyzz_yyzz[i] * pa_y[i];

        tr_z_yyyyzz_yzzz[i] = 3.0 * tr_z_yyzz_yzzz[i] * fe_0 + tr_z_yyyzz_zzz[i] * fe_0 + tr_z_yyyzz_yzzz[i] * pa_y[i];

        tr_z_yyyyzz_zzzz[i] = 3.0 * tr_z_yyzz_zzzz[i] * fe_0 + tr_z_yyyzz_zzzz[i] * pa_y[i];
    }

    // Set up 1200-1215 components of targeted buffer : IG

    auto tr_z_yyyzzz_xxxx = pbuffer.data(idx_dip_ig + 1200);

    auto tr_z_yyyzzz_xxxy = pbuffer.data(idx_dip_ig + 1201);

    auto tr_z_yyyzzz_xxxz = pbuffer.data(idx_dip_ig + 1202);

    auto tr_z_yyyzzz_xxyy = pbuffer.data(idx_dip_ig + 1203);

    auto tr_z_yyyzzz_xxyz = pbuffer.data(idx_dip_ig + 1204);

    auto tr_z_yyyzzz_xxzz = pbuffer.data(idx_dip_ig + 1205);

    auto tr_z_yyyzzz_xyyy = pbuffer.data(idx_dip_ig + 1206);

    auto tr_z_yyyzzz_xyyz = pbuffer.data(idx_dip_ig + 1207);

    auto tr_z_yyyzzz_xyzz = pbuffer.data(idx_dip_ig + 1208);

    auto tr_z_yyyzzz_xzzz = pbuffer.data(idx_dip_ig + 1209);

    auto tr_z_yyyzzz_yyyy = pbuffer.data(idx_dip_ig + 1210);

    auto tr_z_yyyzzz_yyyz = pbuffer.data(idx_dip_ig + 1211);

    auto tr_z_yyyzzz_yyzz = pbuffer.data(idx_dip_ig + 1212);

    auto tr_z_yyyzzz_yzzz = pbuffer.data(idx_dip_ig + 1213);

    auto tr_z_yyyzzz_zzzz = pbuffer.data(idx_dip_ig + 1214);

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyyzzz_xxxx, \
                             tr_z_yyyzzz_xxxy, \
                             tr_z_yyyzzz_xxxz, \
                             tr_z_yyyzzz_xxyy, \
                             tr_z_yyyzzz_xxyz, \
                             tr_z_yyyzzz_xxzz, \
                             tr_z_yyyzzz_xyyy, \
                             tr_z_yyyzzz_xyyz, \
                             tr_z_yyyzzz_xyzz, \
                             tr_z_yyyzzz_xzzz, \
                             tr_z_yyyzzz_yyyy, \
                             tr_z_yyyzzz_yyyz, \
                             tr_z_yyyzzz_yyzz, \
                             tr_z_yyyzzz_yzzz, \
                             tr_z_yyyzzz_zzzz, \
                             tr_z_yyzzz_xxx,   \
                             tr_z_yyzzz_xxxx,  \
                             tr_z_yyzzz_xxxy,  \
                             tr_z_yyzzz_xxxz,  \
                             tr_z_yyzzz_xxy,   \
                             tr_z_yyzzz_xxyy,  \
                             tr_z_yyzzz_xxyz,  \
                             tr_z_yyzzz_xxz,   \
                             tr_z_yyzzz_xxzz,  \
                             tr_z_yyzzz_xyy,   \
                             tr_z_yyzzz_xyyy,  \
                             tr_z_yyzzz_xyyz,  \
                             tr_z_yyzzz_xyz,   \
                             tr_z_yyzzz_xyzz,  \
                             tr_z_yyzzz_xzz,   \
                             tr_z_yyzzz_xzzz,  \
                             tr_z_yyzzz_yyy,   \
                             tr_z_yyzzz_yyyy,  \
                             tr_z_yyzzz_yyyz,  \
                             tr_z_yyzzz_yyz,   \
                             tr_z_yyzzz_yyzz,  \
                             tr_z_yyzzz_yzz,   \
                             tr_z_yyzzz_yzzz,  \
                             tr_z_yyzzz_zzz,   \
                             tr_z_yyzzz_zzzz,  \
                             tr_z_yzzz_xxxx,   \
                             tr_z_yzzz_xxxy,   \
                             tr_z_yzzz_xxxz,   \
                             tr_z_yzzz_xxyy,   \
                             tr_z_yzzz_xxyz,   \
                             tr_z_yzzz_xxzz,   \
                             tr_z_yzzz_xyyy,   \
                             tr_z_yzzz_xyyz,   \
                             tr_z_yzzz_xyzz,   \
                             tr_z_yzzz_xzzz,   \
                             tr_z_yzzz_yyyy,   \
                             tr_z_yzzz_yyyz,   \
                             tr_z_yzzz_yyzz,   \
                             tr_z_yzzz_yzzz,   \
                             tr_z_yzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzzz_xxxx[i] = 2.0 * tr_z_yzzz_xxxx[i] * fe_0 + tr_z_yyzzz_xxxx[i] * pa_y[i];

        tr_z_yyyzzz_xxxy[i] = 2.0 * tr_z_yzzz_xxxy[i] * fe_0 + tr_z_yyzzz_xxx[i] * fe_0 + tr_z_yyzzz_xxxy[i] * pa_y[i];

        tr_z_yyyzzz_xxxz[i] = 2.0 * tr_z_yzzz_xxxz[i] * fe_0 + tr_z_yyzzz_xxxz[i] * pa_y[i];

        tr_z_yyyzzz_xxyy[i] = 2.0 * tr_z_yzzz_xxyy[i] * fe_0 + 2.0 * tr_z_yyzzz_xxy[i] * fe_0 + tr_z_yyzzz_xxyy[i] * pa_y[i];

        tr_z_yyyzzz_xxyz[i] = 2.0 * tr_z_yzzz_xxyz[i] * fe_0 + tr_z_yyzzz_xxz[i] * fe_0 + tr_z_yyzzz_xxyz[i] * pa_y[i];

        tr_z_yyyzzz_xxzz[i] = 2.0 * tr_z_yzzz_xxzz[i] * fe_0 + tr_z_yyzzz_xxzz[i] * pa_y[i];

        tr_z_yyyzzz_xyyy[i] = 2.0 * tr_z_yzzz_xyyy[i] * fe_0 + 3.0 * tr_z_yyzzz_xyy[i] * fe_0 + tr_z_yyzzz_xyyy[i] * pa_y[i];

        tr_z_yyyzzz_xyyz[i] = 2.0 * tr_z_yzzz_xyyz[i] * fe_0 + 2.0 * tr_z_yyzzz_xyz[i] * fe_0 + tr_z_yyzzz_xyyz[i] * pa_y[i];

        tr_z_yyyzzz_xyzz[i] = 2.0 * tr_z_yzzz_xyzz[i] * fe_0 + tr_z_yyzzz_xzz[i] * fe_0 + tr_z_yyzzz_xyzz[i] * pa_y[i];

        tr_z_yyyzzz_xzzz[i] = 2.0 * tr_z_yzzz_xzzz[i] * fe_0 + tr_z_yyzzz_xzzz[i] * pa_y[i];

        tr_z_yyyzzz_yyyy[i] = 2.0 * tr_z_yzzz_yyyy[i] * fe_0 + 4.0 * tr_z_yyzzz_yyy[i] * fe_0 + tr_z_yyzzz_yyyy[i] * pa_y[i];

        tr_z_yyyzzz_yyyz[i] = 2.0 * tr_z_yzzz_yyyz[i] * fe_0 + 3.0 * tr_z_yyzzz_yyz[i] * fe_0 + tr_z_yyzzz_yyyz[i] * pa_y[i];

        tr_z_yyyzzz_yyzz[i] = 2.0 * tr_z_yzzz_yyzz[i] * fe_0 + 2.0 * tr_z_yyzzz_yzz[i] * fe_0 + tr_z_yyzzz_yyzz[i] * pa_y[i];

        tr_z_yyyzzz_yzzz[i] = 2.0 * tr_z_yzzz_yzzz[i] * fe_0 + tr_z_yyzzz_zzz[i] * fe_0 + tr_z_yyzzz_yzzz[i] * pa_y[i];

        tr_z_yyyzzz_zzzz[i] = 2.0 * tr_z_yzzz_zzzz[i] * fe_0 + tr_z_yyzzz_zzzz[i] * pa_y[i];
    }

    // Set up 1215-1230 components of targeted buffer : IG

    auto tr_z_yyzzzz_xxxx = pbuffer.data(idx_dip_ig + 1215);

    auto tr_z_yyzzzz_xxxy = pbuffer.data(idx_dip_ig + 1216);

    auto tr_z_yyzzzz_xxxz = pbuffer.data(idx_dip_ig + 1217);

    auto tr_z_yyzzzz_xxyy = pbuffer.data(idx_dip_ig + 1218);

    auto tr_z_yyzzzz_xxyz = pbuffer.data(idx_dip_ig + 1219);

    auto tr_z_yyzzzz_xxzz = pbuffer.data(idx_dip_ig + 1220);

    auto tr_z_yyzzzz_xyyy = pbuffer.data(idx_dip_ig + 1221);

    auto tr_z_yyzzzz_xyyz = pbuffer.data(idx_dip_ig + 1222);

    auto tr_z_yyzzzz_xyzz = pbuffer.data(idx_dip_ig + 1223);

    auto tr_z_yyzzzz_xzzz = pbuffer.data(idx_dip_ig + 1224);

    auto tr_z_yyzzzz_yyyy = pbuffer.data(idx_dip_ig + 1225);

    auto tr_z_yyzzzz_yyyz = pbuffer.data(idx_dip_ig + 1226);

    auto tr_z_yyzzzz_yyzz = pbuffer.data(idx_dip_ig + 1227);

    auto tr_z_yyzzzz_yzzz = pbuffer.data(idx_dip_ig + 1228);

    auto tr_z_yyzzzz_zzzz = pbuffer.data(idx_dip_ig + 1229);

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyzzzz_xxxx, \
                             tr_z_yyzzzz_xxxy, \
                             tr_z_yyzzzz_xxxz, \
                             tr_z_yyzzzz_xxyy, \
                             tr_z_yyzzzz_xxyz, \
                             tr_z_yyzzzz_xxzz, \
                             tr_z_yyzzzz_xyyy, \
                             tr_z_yyzzzz_xyyz, \
                             tr_z_yyzzzz_xyzz, \
                             tr_z_yyzzzz_xzzz, \
                             tr_z_yyzzzz_yyyy, \
                             tr_z_yyzzzz_yyyz, \
                             tr_z_yyzzzz_yyzz, \
                             tr_z_yyzzzz_yzzz, \
                             tr_z_yyzzzz_zzzz, \
                             tr_z_yzzzz_xxx,   \
                             tr_z_yzzzz_xxxx,  \
                             tr_z_yzzzz_xxxy,  \
                             tr_z_yzzzz_xxxz,  \
                             tr_z_yzzzz_xxy,   \
                             tr_z_yzzzz_xxyy,  \
                             tr_z_yzzzz_xxyz,  \
                             tr_z_yzzzz_xxz,   \
                             tr_z_yzzzz_xxzz,  \
                             tr_z_yzzzz_xyy,   \
                             tr_z_yzzzz_xyyy,  \
                             tr_z_yzzzz_xyyz,  \
                             tr_z_yzzzz_xyz,   \
                             tr_z_yzzzz_xyzz,  \
                             tr_z_yzzzz_xzz,   \
                             tr_z_yzzzz_xzzz,  \
                             tr_z_yzzzz_yyy,   \
                             tr_z_yzzzz_yyyy,  \
                             tr_z_yzzzz_yyyz,  \
                             tr_z_yzzzz_yyz,   \
                             tr_z_yzzzz_yyzz,  \
                             tr_z_yzzzz_yzz,   \
                             tr_z_yzzzz_yzzz,  \
                             tr_z_yzzzz_zzz,   \
                             tr_z_yzzzz_zzzz,  \
                             tr_z_zzzz_xxxx,   \
                             tr_z_zzzz_xxxy,   \
                             tr_z_zzzz_xxxz,   \
                             tr_z_zzzz_xxyy,   \
                             tr_z_zzzz_xxyz,   \
                             tr_z_zzzz_xxzz,   \
                             tr_z_zzzz_xyyy,   \
                             tr_z_zzzz_xyyz,   \
                             tr_z_zzzz_xyzz,   \
                             tr_z_zzzz_xzzz,   \
                             tr_z_zzzz_yyyy,   \
                             tr_z_zzzz_yyyz,   \
                             tr_z_zzzz_yyzz,   \
                             tr_z_zzzz_yzzz,   \
                             tr_z_zzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzzz_xxxx[i] = tr_z_zzzz_xxxx[i] * fe_0 + tr_z_yzzzz_xxxx[i] * pa_y[i];

        tr_z_yyzzzz_xxxy[i] = tr_z_zzzz_xxxy[i] * fe_0 + tr_z_yzzzz_xxx[i] * fe_0 + tr_z_yzzzz_xxxy[i] * pa_y[i];

        tr_z_yyzzzz_xxxz[i] = tr_z_zzzz_xxxz[i] * fe_0 + tr_z_yzzzz_xxxz[i] * pa_y[i];

        tr_z_yyzzzz_xxyy[i] = tr_z_zzzz_xxyy[i] * fe_0 + 2.0 * tr_z_yzzzz_xxy[i] * fe_0 + tr_z_yzzzz_xxyy[i] * pa_y[i];

        tr_z_yyzzzz_xxyz[i] = tr_z_zzzz_xxyz[i] * fe_0 + tr_z_yzzzz_xxz[i] * fe_0 + tr_z_yzzzz_xxyz[i] * pa_y[i];

        tr_z_yyzzzz_xxzz[i] = tr_z_zzzz_xxzz[i] * fe_0 + tr_z_yzzzz_xxzz[i] * pa_y[i];

        tr_z_yyzzzz_xyyy[i] = tr_z_zzzz_xyyy[i] * fe_0 + 3.0 * tr_z_yzzzz_xyy[i] * fe_0 + tr_z_yzzzz_xyyy[i] * pa_y[i];

        tr_z_yyzzzz_xyyz[i] = tr_z_zzzz_xyyz[i] * fe_0 + 2.0 * tr_z_yzzzz_xyz[i] * fe_0 + tr_z_yzzzz_xyyz[i] * pa_y[i];

        tr_z_yyzzzz_xyzz[i] = tr_z_zzzz_xyzz[i] * fe_0 + tr_z_yzzzz_xzz[i] * fe_0 + tr_z_yzzzz_xyzz[i] * pa_y[i];

        tr_z_yyzzzz_xzzz[i] = tr_z_zzzz_xzzz[i] * fe_0 + tr_z_yzzzz_xzzz[i] * pa_y[i];

        tr_z_yyzzzz_yyyy[i] = tr_z_zzzz_yyyy[i] * fe_0 + 4.0 * tr_z_yzzzz_yyy[i] * fe_0 + tr_z_yzzzz_yyyy[i] * pa_y[i];

        tr_z_yyzzzz_yyyz[i] = tr_z_zzzz_yyyz[i] * fe_0 + 3.0 * tr_z_yzzzz_yyz[i] * fe_0 + tr_z_yzzzz_yyyz[i] * pa_y[i];

        tr_z_yyzzzz_yyzz[i] = tr_z_zzzz_yyzz[i] * fe_0 + 2.0 * tr_z_yzzzz_yzz[i] * fe_0 + tr_z_yzzzz_yyzz[i] * pa_y[i];

        tr_z_yyzzzz_yzzz[i] = tr_z_zzzz_yzzz[i] * fe_0 + tr_z_yzzzz_zzz[i] * fe_0 + tr_z_yzzzz_yzzz[i] * pa_y[i];

        tr_z_yyzzzz_zzzz[i] = tr_z_zzzz_zzzz[i] * fe_0 + tr_z_yzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 1230-1245 components of targeted buffer : IG

    auto tr_z_yzzzzz_xxxx = pbuffer.data(idx_dip_ig + 1230);

    auto tr_z_yzzzzz_xxxy = pbuffer.data(idx_dip_ig + 1231);

    auto tr_z_yzzzzz_xxxz = pbuffer.data(idx_dip_ig + 1232);

    auto tr_z_yzzzzz_xxyy = pbuffer.data(idx_dip_ig + 1233);

    auto tr_z_yzzzzz_xxyz = pbuffer.data(idx_dip_ig + 1234);

    auto tr_z_yzzzzz_xxzz = pbuffer.data(idx_dip_ig + 1235);

    auto tr_z_yzzzzz_xyyy = pbuffer.data(idx_dip_ig + 1236);

    auto tr_z_yzzzzz_xyyz = pbuffer.data(idx_dip_ig + 1237);

    auto tr_z_yzzzzz_xyzz = pbuffer.data(idx_dip_ig + 1238);

    auto tr_z_yzzzzz_xzzz = pbuffer.data(idx_dip_ig + 1239);

    auto tr_z_yzzzzz_yyyy = pbuffer.data(idx_dip_ig + 1240);

    auto tr_z_yzzzzz_yyyz = pbuffer.data(idx_dip_ig + 1241);

    auto tr_z_yzzzzz_yyzz = pbuffer.data(idx_dip_ig + 1242);

    auto tr_z_yzzzzz_yzzz = pbuffer.data(idx_dip_ig + 1243);

    auto tr_z_yzzzzz_zzzz = pbuffer.data(idx_dip_ig + 1244);

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yzzzzz_xxxx, \
                             tr_z_yzzzzz_xxxy, \
                             tr_z_yzzzzz_xxxz, \
                             tr_z_yzzzzz_xxyy, \
                             tr_z_yzzzzz_xxyz, \
                             tr_z_yzzzzz_xxzz, \
                             tr_z_yzzzzz_xyyy, \
                             tr_z_yzzzzz_xyyz, \
                             tr_z_yzzzzz_xyzz, \
                             tr_z_yzzzzz_xzzz, \
                             tr_z_yzzzzz_yyyy, \
                             tr_z_yzzzzz_yyyz, \
                             tr_z_yzzzzz_yyzz, \
                             tr_z_yzzzzz_yzzz, \
                             tr_z_yzzzzz_zzzz, \
                             tr_z_zzzzz_xxx,   \
                             tr_z_zzzzz_xxxx,  \
                             tr_z_zzzzz_xxxy,  \
                             tr_z_zzzzz_xxxz,  \
                             tr_z_zzzzz_xxy,   \
                             tr_z_zzzzz_xxyy,  \
                             tr_z_zzzzz_xxyz,  \
                             tr_z_zzzzz_xxz,   \
                             tr_z_zzzzz_xxzz,  \
                             tr_z_zzzzz_xyy,   \
                             tr_z_zzzzz_xyyy,  \
                             tr_z_zzzzz_xyyz,  \
                             tr_z_zzzzz_xyz,   \
                             tr_z_zzzzz_xyzz,  \
                             tr_z_zzzzz_xzz,   \
                             tr_z_zzzzz_xzzz,  \
                             tr_z_zzzzz_yyy,   \
                             tr_z_zzzzz_yyyy,  \
                             tr_z_zzzzz_yyyz,  \
                             tr_z_zzzzz_yyz,   \
                             tr_z_zzzzz_yyzz,  \
                             tr_z_zzzzz_yzz,   \
                             tr_z_zzzzz_yzzz,  \
                             tr_z_zzzzz_zzz,   \
                             tr_z_zzzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzzz_xxxx[i] = tr_z_zzzzz_xxxx[i] * pa_y[i];

        tr_z_yzzzzz_xxxy[i] = tr_z_zzzzz_xxx[i] * fe_0 + tr_z_zzzzz_xxxy[i] * pa_y[i];

        tr_z_yzzzzz_xxxz[i] = tr_z_zzzzz_xxxz[i] * pa_y[i];

        tr_z_yzzzzz_xxyy[i] = 2.0 * tr_z_zzzzz_xxy[i] * fe_0 + tr_z_zzzzz_xxyy[i] * pa_y[i];

        tr_z_yzzzzz_xxyz[i] = tr_z_zzzzz_xxz[i] * fe_0 + tr_z_zzzzz_xxyz[i] * pa_y[i];

        tr_z_yzzzzz_xxzz[i] = tr_z_zzzzz_xxzz[i] * pa_y[i];

        tr_z_yzzzzz_xyyy[i] = 3.0 * tr_z_zzzzz_xyy[i] * fe_0 + tr_z_zzzzz_xyyy[i] * pa_y[i];

        tr_z_yzzzzz_xyyz[i] = 2.0 * tr_z_zzzzz_xyz[i] * fe_0 + tr_z_zzzzz_xyyz[i] * pa_y[i];

        tr_z_yzzzzz_xyzz[i] = tr_z_zzzzz_xzz[i] * fe_0 + tr_z_zzzzz_xyzz[i] * pa_y[i];

        tr_z_yzzzzz_xzzz[i] = tr_z_zzzzz_xzzz[i] * pa_y[i];

        tr_z_yzzzzz_yyyy[i] = 4.0 * tr_z_zzzzz_yyy[i] * fe_0 + tr_z_zzzzz_yyyy[i] * pa_y[i];

        tr_z_yzzzzz_yyyz[i] = 3.0 * tr_z_zzzzz_yyz[i] * fe_0 + tr_z_zzzzz_yyyz[i] * pa_y[i];

        tr_z_yzzzzz_yyzz[i] = 2.0 * tr_z_zzzzz_yzz[i] * fe_0 + tr_z_zzzzz_yyzz[i] * pa_y[i];

        tr_z_yzzzzz_yzzz[i] = tr_z_zzzzz_zzz[i] * fe_0 + tr_z_zzzzz_yzzz[i] * pa_y[i];

        tr_z_yzzzzz_zzzz[i] = tr_z_zzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 1245-1260 components of targeted buffer : IG

    auto tr_z_zzzzzz_xxxx = pbuffer.data(idx_dip_ig + 1245);

    auto tr_z_zzzzzz_xxxy = pbuffer.data(idx_dip_ig + 1246);

    auto tr_z_zzzzzz_xxxz = pbuffer.data(idx_dip_ig + 1247);

    auto tr_z_zzzzzz_xxyy = pbuffer.data(idx_dip_ig + 1248);

    auto tr_z_zzzzzz_xxyz = pbuffer.data(idx_dip_ig + 1249);

    auto tr_z_zzzzzz_xxzz = pbuffer.data(idx_dip_ig + 1250);

    auto tr_z_zzzzzz_xyyy = pbuffer.data(idx_dip_ig + 1251);

    auto tr_z_zzzzzz_xyyz = pbuffer.data(idx_dip_ig + 1252);

    auto tr_z_zzzzzz_xyzz = pbuffer.data(idx_dip_ig + 1253);

    auto tr_z_zzzzzz_xzzz = pbuffer.data(idx_dip_ig + 1254);

    auto tr_z_zzzzzz_yyyy = pbuffer.data(idx_dip_ig + 1255);

    auto tr_z_zzzzzz_yyyz = pbuffer.data(idx_dip_ig + 1256);

    auto tr_z_zzzzzz_yyzz = pbuffer.data(idx_dip_ig + 1257);

    auto tr_z_zzzzzz_yzzz = pbuffer.data(idx_dip_ig + 1258);

    auto tr_z_zzzzzz_zzzz = pbuffer.data(idx_dip_ig + 1259);

#pragma omp simd aligned(pa_z,                 \
                             tr_z_zzzz_xxxx,   \
                             tr_z_zzzz_xxxy,   \
                             tr_z_zzzz_xxxz,   \
                             tr_z_zzzz_xxyy,   \
                             tr_z_zzzz_xxyz,   \
                             tr_z_zzzz_xxzz,   \
                             tr_z_zzzz_xyyy,   \
                             tr_z_zzzz_xyyz,   \
                             tr_z_zzzz_xyzz,   \
                             tr_z_zzzz_xzzz,   \
                             tr_z_zzzz_yyyy,   \
                             tr_z_zzzz_yyyz,   \
                             tr_z_zzzz_yyzz,   \
                             tr_z_zzzz_yzzz,   \
                             tr_z_zzzz_zzzz,   \
                             tr_z_zzzzz_xxx,   \
                             tr_z_zzzzz_xxxx,  \
                             tr_z_zzzzz_xxxy,  \
                             tr_z_zzzzz_xxxz,  \
                             tr_z_zzzzz_xxy,   \
                             tr_z_zzzzz_xxyy,  \
                             tr_z_zzzzz_xxyz,  \
                             tr_z_zzzzz_xxz,   \
                             tr_z_zzzzz_xxzz,  \
                             tr_z_zzzzz_xyy,   \
                             tr_z_zzzzz_xyyy,  \
                             tr_z_zzzzz_xyyz,  \
                             tr_z_zzzzz_xyz,   \
                             tr_z_zzzzz_xyzz,  \
                             tr_z_zzzzz_xzz,   \
                             tr_z_zzzzz_xzzz,  \
                             tr_z_zzzzz_yyy,   \
                             tr_z_zzzzz_yyyy,  \
                             tr_z_zzzzz_yyyz,  \
                             tr_z_zzzzz_yyz,   \
                             tr_z_zzzzz_yyzz,  \
                             tr_z_zzzzz_yzz,   \
                             tr_z_zzzzz_yzzz,  \
                             tr_z_zzzzz_zzz,   \
                             tr_z_zzzzz_zzzz,  \
                             tr_z_zzzzzz_xxxx, \
                             tr_z_zzzzzz_xxxy, \
                             tr_z_zzzzzz_xxxz, \
                             tr_z_zzzzzz_xxyy, \
                             tr_z_zzzzzz_xxyz, \
                             tr_z_zzzzzz_xxzz, \
                             tr_z_zzzzzz_xyyy, \
                             tr_z_zzzzzz_xyyz, \
                             tr_z_zzzzzz_xyzz, \
                             tr_z_zzzzzz_xzzz, \
                             tr_z_zzzzzz_yyyy, \
                             tr_z_zzzzzz_yyyz, \
                             tr_z_zzzzzz_yyzz, \
                             tr_z_zzzzzz_yzzz, \
                             tr_z_zzzzzz_zzzz, \
                             ts_zzzzz_xxxx,    \
                             ts_zzzzz_xxxy,    \
                             ts_zzzzz_xxxz,    \
                             ts_zzzzz_xxyy,    \
                             ts_zzzzz_xxyz,    \
                             ts_zzzzz_xxzz,    \
                             ts_zzzzz_xyyy,    \
                             ts_zzzzz_xyyz,    \
                             ts_zzzzz_xyzz,    \
                             ts_zzzzz_xzzz,    \
                             ts_zzzzz_yyyy,    \
                             ts_zzzzz_yyyz,    \
                             ts_zzzzz_yyzz,    \
                             ts_zzzzz_yzzz,    \
                             ts_zzzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzzz_xxxx[i] = 5.0 * tr_z_zzzz_xxxx[i] * fe_0 + ts_zzzzz_xxxx[i] * fe_0 + tr_z_zzzzz_xxxx[i] * pa_z[i];

        tr_z_zzzzzz_xxxy[i] = 5.0 * tr_z_zzzz_xxxy[i] * fe_0 + ts_zzzzz_xxxy[i] * fe_0 + tr_z_zzzzz_xxxy[i] * pa_z[i];

        tr_z_zzzzzz_xxxz[i] = 5.0 * tr_z_zzzz_xxxz[i] * fe_0 + tr_z_zzzzz_xxx[i] * fe_0 + ts_zzzzz_xxxz[i] * fe_0 + tr_z_zzzzz_xxxz[i] * pa_z[i];

        tr_z_zzzzzz_xxyy[i] = 5.0 * tr_z_zzzz_xxyy[i] * fe_0 + ts_zzzzz_xxyy[i] * fe_0 + tr_z_zzzzz_xxyy[i] * pa_z[i];

        tr_z_zzzzzz_xxyz[i] = 5.0 * tr_z_zzzz_xxyz[i] * fe_0 + tr_z_zzzzz_xxy[i] * fe_0 + ts_zzzzz_xxyz[i] * fe_0 + tr_z_zzzzz_xxyz[i] * pa_z[i];

        tr_z_zzzzzz_xxzz[i] =
            5.0 * tr_z_zzzz_xxzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xxz[i] * fe_0 + ts_zzzzz_xxzz[i] * fe_0 + tr_z_zzzzz_xxzz[i] * pa_z[i];

        tr_z_zzzzzz_xyyy[i] = 5.0 * tr_z_zzzz_xyyy[i] * fe_0 + ts_zzzzz_xyyy[i] * fe_0 + tr_z_zzzzz_xyyy[i] * pa_z[i];

        tr_z_zzzzzz_xyyz[i] = 5.0 * tr_z_zzzz_xyyz[i] * fe_0 + tr_z_zzzzz_xyy[i] * fe_0 + ts_zzzzz_xyyz[i] * fe_0 + tr_z_zzzzz_xyyz[i] * pa_z[i];

        tr_z_zzzzzz_xyzz[i] =
            5.0 * tr_z_zzzz_xyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xyz[i] * fe_0 + ts_zzzzz_xyzz[i] * fe_0 + tr_z_zzzzz_xyzz[i] * pa_z[i];

        tr_z_zzzzzz_xzzz[i] =
            5.0 * tr_z_zzzz_xzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_xzz[i] * fe_0 + ts_zzzzz_xzzz[i] * fe_0 + tr_z_zzzzz_xzzz[i] * pa_z[i];

        tr_z_zzzzzz_yyyy[i] = 5.0 * tr_z_zzzz_yyyy[i] * fe_0 + ts_zzzzz_yyyy[i] * fe_0 + tr_z_zzzzz_yyyy[i] * pa_z[i];

        tr_z_zzzzzz_yyyz[i] = 5.0 * tr_z_zzzz_yyyz[i] * fe_0 + tr_z_zzzzz_yyy[i] * fe_0 + ts_zzzzz_yyyz[i] * fe_0 + tr_z_zzzzz_yyyz[i] * pa_z[i];

        tr_z_zzzzzz_yyzz[i] =
            5.0 * tr_z_zzzz_yyzz[i] * fe_0 + 2.0 * tr_z_zzzzz_yyz[i] * fe_0 + ts_zzzzz_yyzz[i] * fe_0 + tr_z_zzzzz_yyzz[i] * pa_z[i];

        tr_z_zzzzzz_yzzz[i] =
            5.0 * tr_z_zzzz_yzzz[i] * fe_0 + 3.0 * tr_z_zzzzz_yzz[i] * fe_0 + ts_zzzzz_yzzz[i] * fe_0 + tr_z_zzzzz_yzzz[i] * pa_z[i];

        tr_z_zzzzzz_zzzz[i] =
            5.0 * tr_z_zzzz_zzzz[i] * fe_0 + 4.0 * tr_z_zzzzz_zzz[i] * fe_0 + ts_zzzzz_zzzz[i] * fe_0 + tr_z_zzzzz_zzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
