#include "ElectricDipoleMomentumPrimRecGH.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_gh(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_gh,
                                      const size_t idx_dip_dh,
                                      const size_t idx_dip_fg,
                                      const size_t idx_ovl_fh,
                                      const size_t idx_dip_fh,
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

    // Set up components of auxiliary buffer : DH

    auto tr_x_xx_xxxxx = pbuffer.data(idx_dip_dh);

    auto tr_x_xx_xxxxy = pbuffer.data(idx_dip_dh + 1);

    auto tr_x_xx_xxxxz = pbuffer.data(idx_dip_dh + 2);

    auto tr_x_xx_xxxyy = pbuffer.data(idx_dip_dh + 3);

    auto tr_x_xx_xxxyz = pbuffer.data(idx_dip_dh + 4);

    auto tr_x_xx_xxxzz = pbuffer.data(idx_dip_dh + 5);

    auto tr_x_xx_xxyyy = pbuffer.data(idx_dip_dh + 6);

    auto tr_x_xx_xxyyz = pbuffer.data(idx_dip_dh + 7);

    auto tr_x_xx_xxyzz = pbuffer.data(idx_dip_dh + 8);

    auto tr_x_xx_xxzzz = pbuffer.data(idx_dip_dh + 9);

    auto tr_x_xx_xyyyy = pbuffer.data(idx_dip_dh + 10);

    auto tr_x_xx_xyyyz = pbuffer.data(idx_dip_dh + 11);

    auto tr_x_xx_xyyzz = pbuffer.data(idx_dip_dh + 12);

    auto tr_x_xx_xyzzz = pbuffer.data(idx_dip_dh + 13);

    auto tr_x_xx_xzzzz = pbuffer.data(idx_dip_dh + 14);

    auto tr_x_xx_yyyyy = pbuffer.data(idx_dip_dh + 15);

    auto tr_x_xx_yyyyz = pbuffer.data(idx_dip_dh + 16);

    auto tr_x_xx_yyyzz = pbuffer.data(idx_dip_dh + 17);

    auto tr_x_xx_yyzzz = pbuffer.data(idx_dip_dh + 18);

    auto tr_x_xx_yzzzz = pbuffer.data(idx_dip_dh + 19);

    auto tr_x_xx_zzzzz = pbuffer.data(idx_dip_dh + 20);

    auto tr_x_xy_xxxxx = pbuffer.data(idx_dip_dh + 21);

    auto tr_x_xy_xxxxz = pbuffer.data(idx_dip_dh + 23);

    auto tr_x_xy_xxxzz = pbuffer.data(idx_dip_dh + 26);

    auto tr_x_xy_xxzzz = pbuffer.data(idx_dip_dh + 30);

    auto tr_x_xy_xzzzz = pbuffer.data(idx_dip_dh + 35);

    auto tr_x_xz_xxxxx = pbuffer.data(idx_dip_dh + 42);

    auto tr_x_xz_xxxxy = pbuffer.data(idx_dip_dh + 43);

    auto tr_x_xz_xxxxz = pbuffer.data(idx_dip_dh + 44);

    auto tr_x_xz_xxxyy = pbuffer.data(idx_dip_dh + 45);

    auto tr_x_xz_xxxzz = pbuffer.data(idx_dip_dh + 47);

    auto tr_x_xz_xxyyy = pbuffer.data(idx_dip_dh + 48);

    auto tr_x_xz_xxzzz = pbuffer.data(idx_dip_dh + 51);

    auto tr_x_xz_xyyyy = pbuffer.data(idx_dip_dh + 52);

    auto tr_x_xz_xzzzz = pbuffer.data(idx_dip_dh + 56);

    auto tr_x_yy_xxxxx = pbuffer.data(idx_dip_dh + 63);

    auto tr_x_yy_xxxxy = pbuffer.data(idx_dip_dh + 64);

    auto tr_x_yy_xxxxz = pbuffer.data(idx_dip_dh + 65);

    auto tr_x_yy_xxxyy = pbuffer.data(idx_dip_dh + 66);

    auto tr_x_yy_xxxyz = pbuffer.data(idx_dip_dh + 67);

    auto tr_x_yy_xxxzz = pbuffer.data(idx_dip_dh + 68);

    auto tr_x_yy_xxyyy = pbuffer.data(idx_dip_dh + 69);

    auto tr_x_yy_xxyyz = pbuffer.data(idx_dip_dh + 70);

    auto tr_x_yy_xxyzz = pbuffer.data(idx_dip_dh + 71);

    auto tr_x_yy_xxzzz = pbuffer.data(idx_dip_dh + 72);

    auto tr_x_yy_xyyyy = pbuffer.data(idx_dip_dh + 73);

    auto tr_x_yy_xyyyz = pbuffer.data(idx_dip_dh + 74);

    auto tr_x_yy_xyyzz = pbuffer.data(idx_dip_dh + 75);

    auto tr_x_yy_xyzzz = pbuffer.data(idx_dip_dh + 76);

    auto tr_x_yy_xzzzz = pbuffer.data(idx_dip_dh + 77);

    auto tr_x_yy_yyyyy = pbuffer.data(idx_dip_dh + 78);

    auto tr_x_yy_yyyyz = pbuffer.data(idx_dip_dh + 79);

    auto tr_x_yy_yyyzz = pbuffer.data(idx_dip_dh + 80);

    auto tr_x_yy_yyzzz = pbuffer.data(idx_dip_dh + 81);

    auto tr_x_yy_yzzzz = pbuffer.data(idx_dip_dh + 82);

    auto tr_x_yy_zzzzz = pbuffer.data(idx_dip_dh + 83);

    auto tr_x_yz_xxxxz = pbuffer.data(idx_dip_dh + 86);

    auto tr_x_yz_xxxzz = pbuffer.data(idx_dip_dh + 89);

    auto tr_x_yz_xxzzz = pbuffer.data(idx_dip_dh + 93);

    auto tr_x_yz_xzzzz = pbuffer.data(idx_dip_dh + 98);

    auto tr_x_yz_zzzzz = pbuffer.data(idx_dip_dh + 104);

    auto tr_x_zz_xxxxx = pbuffer.data(idx_dip_dh + 105);

    auto tr_x_zz_xxxxy = pbuffer.data(idx_dip_dh + 106);

    auto tr_x_zz_xxxxz = pbuffer.data(idx_dip_dh + 107);

    auto tr_x_zz_xxxyy = pbuffer.data(idx_dip_dh + 108);

    auto tr_x_zz_xxxyz = pbuffer.data(idx_dip_dh + 109);

    auto tr_x_zz_xxxzz = pbuffer.data(idx_dip_dh + 110);

    auto tr_x_zz_xxyyy = pbuffer.data(idx_dip_dh + 111);

    auto tr_x_zz_xxyyz = pbuffer.data(idx_dip_dh + 112);

    auto tr_x_zz_xxyzz = pbuffer.data(idx_dip_dh + 113);

    auto tr_x_zz_xxzzz = pbuffer.data(idx_dip_dh + 114);

    auto tr_x_zz_xyyyy = pbuffer.data(idx_dip_dh + 115);

    auto tr_x_zz_xyyyz = pbuffer.data(idx_dip_dh + 116);

    auto tr_x_zz_xyyzz = pbuffer.data(idx_dip_dh + 117);

    auto tr_x_zz_xyzzz = pbuffer.data(idx_dip_dh + 118);

    auto tr_x_zz_xzzzz = pbuffer.data(idx_dip_dh + 119);

    auto tr_x_zz_yyyyy = pbuffer.data(idx_dip_dh + 120);

    auto tr_x_zz_yyyyz = pbuffer.data(idx_dip_dh + 121);

    auto tr_x_zz_yyyzz = pbuffer.data(idx_dip_dh + 122);

    auto tr_x_zz_yyzzz = pbuffer.data(idx_dip_dh + 123);

    auto tr_x_zz_yzzzz = pbuffer.data(idx_dip_dh + 124);

    auto tr_x_zz_zzzzz = pbuffer.data(idx_dip_dh + 125);

    auto tr_y_xx_xxxxx = pbuffer.data(idx_dip_dh + 126);

    auto tr_y_xx_xxxxy = pbuffer.data(idx_dip_dh + 127);

    auto tr_y_xx_xxxxz = pbuffer.data(idx_dip_dh + 128);

    auto tr_y_xx_xxxyy = pbuffer.data(idx_dip_dh + 129);

    auto tr_y_xx_xxxyz = pbuffer.data(idx_dip_dh + 130);

    auto tr_y_xx_xxxzz = pbuffer.data(idx_dip_dh + 131);

    auto tr_y_xx_xxyyy = pbuffer.data(idx_dip_dh + 132);

    auto tr_y_xx_xxyyz = pbuffer.data(idx_dip_dh + 133);

    auto tr_y_xx_xxyzz = pbuffer.data(idx_dip_dh + 134);

    auto tr_y_xx_xxzzz = pbuffer.data(idx_dip_dh + 135);

    auto tr_y_xx_xyyyy = pbuffer.data(idx_dip_dh + 136);

    auto tr_y_xx_xyyyz = pbuffer.data(idx_dip_dh + 137);

    auto tr_y_xx_xyyzz = pbuffer.data(idx_dip_dh + 138);

    auto tr_y_xx_xyzzz = pbuffer.data(idx_dip_dh + 139);

    auto tr_y_xx_xzzzz = pbuffer.data(idx_dip_dh + 140);

    auto tr_y_xx_yyyyy = pbuffer.data(idx_dip_dh + 141);

    auto tr_y_xx_yyyyz = pbuffer.data(idx_dip_dh + 142);

    auto tr_y_xx_yyyzz = pbuffer.data(idx_dip_dh + 143);

    auto tr_y_xx_yyzzz = pbuffer.data(idx_dip_dh + 144);

    auto tr_y_xx_yzzzz = pbuffer.data(idx_dip_dh + 145);

    auto tr_y_xx_zzzzz = pbuffer.data(idx_dip_dh + 146);

    auto tr_y_xy_xxxxy = pbuffer.data(idx_dip_dh + 148);

    auto tr_y_xy_xxxyy = pbuffer.data(idx_dip_dh + 150);

    auto tr_y_xy_xxxyz = pbuffer.data(idx_dip_dh + 151);

    auto tr_y_xy_xxyyy = pbuffer.data(idx_dip_dh + 153);

    auto tr_y_xy_xxyyz = pbuffer.data(idx_dip_dh + 154);

    auto tr_y_xy_xxyzz = pbuffer.data(idx_dip_dh + 155);

    auto tr_y_xy_xyyyy = pbuffer.data(idx_dip_dh + 157);

    auto tr_y_xy_xyyyz = pbuffer.data(idx_dip_dh + 158);

    auto tr_y_xy_xyyzz = pbuffer.data(idx_dip_dh + 159);

    auto tr_y_xy_xyzzz = pbuffer.data(idx_dip_dh + 160);

    auto tr_y_xy_yyyyy = pbuffer.data(idx_dip_dh + 162);

    auto tr_y_xy_yyyyz = pbuffer.data(idx_dip_dh + 163);

    auto tr_y_xy_yyyzz = pbuffer.data(idx_dip_dh + 164);

    auto tr_y_xy_yyzzz = pbuffer.data(idx_dip_dh + 165);

    auto tr_y_xy_yzzzz = pbuffer.data(idx_dip_dh + 166);

    auto tr_y_xy_zzzzz = pbuffer.data(idx_dip_dh + 167);

    auto tr_y_xz_yyyyz = pbuffer.data(idx_dip_dh + 184);

    auto tr_y_xz_yyyzz = pbuffer.data(idx_dip_dh + 185);

    auto tr_y_xz_yyzzz = pbuffer.data(idx_dip_dh + 186);

    auto tr_y_xz_yzzzz = pbuffer.data(idx_dip_dh + 187);

    auto tr_y_xz_zzzzz = pbuffer.data(idx_dip_dh + 188);

    auto tr_y_yy_xxxxx = pbuffer.data(idx_dip_dh + 189);

    auto tr_y_yy_xxxxy = pbuffer.data(idx_dip_dh + 190);

    auto tr_y_yy_xxxxz = pbuffer.data(idx_dip_dh + 191);

    auto tr_y_yy_xxxyy = pbuffer.data(idx_dip_dh + 192);

    auto tr_y_yy_xxxyz = pbuffer.data(idx_dip_dh + 193);

    auto tr_y_yy_xxxzz = pbuffer.data(idx_dip_dh + 194);

    auto tr_y_yy_xxyyy = pbuffer.data(idx_dip_dh + 195);

    auto tr_y_yy_xxyyz = pbuffer.data(idx_dip_dh + 196);

    auto tr_y_yy_xxyzz = pbuffer.data(idx_dip_dh + 197);

    auto tr_y_yy_xxzzz = pbuffer.data(idx_dip_dh + 198);

    auto tr_y_yy_xyyyy = pbuffer.data(idx_dip_dh + 199);

    auto tr_y_yy_xyyyz = pbuffer.data(idx_dip_dh + 200);

    auto tr_y_yy_xyyzz = pbuffer.data(idx_dip_dh + 201);

    auto tr_y_yy_xyzzz = pbuffer.data(idx_dip_dh + 202);

    auto tr_y_yy_xzzzz = pbuffer.data(idx_dip_dh + 203);

    auto tr_y_yy_yyyyy = pbuffer.data(idx_dip_dh + 204);

    auto tr_y_yy_yyyyz = pbuffer.data(idx_dip_dh + 205);

    auto tr_y_yy_yyyzz = pbuffer.data(idx_dip_dh + 206);

    auto tr_y_yy_yyzzz = pbuffer.data(idx_dip_dh + 207);

    auto tr_y_yy_yzzzz = pbuffer.data(idx_dip_dh + 208);

    auto tr_y_yy_zzzzz = pbuffer.data(idx_dip_dh + 209);

    auto tr_y_yz_xxxxy = pbuffer.data(idx_dip_dh + 211);

    auto tr_y_yz_xxxyy = pbuffer.data(idx_dip_dh + 213);

    auto tr_y_yz_xxyyy = pbuffer.data(idx_dip_dh + 216);

    auto tr_y_yz_xyyyy = pbuffer.data(idx_dip_dh + 220);

    auto tr_y_yz_yyyyy = pbuffer.data(idx_dip_dh + 225);

    auto tr_y_yz_yyyyz = pbuffer.data(idx_dip_dh + 226);

    auto tr_y_yz_yyyzz = pbuffer.data(idx_dip_dh + 227);

    auto tr_y_yz_yyzzz = pbuffer.data(idx_dip_dh + 228);

    auto tr_y_yz_yzzzz = pbuffer.data(idx_dip_dh + 229);

    auto tr_y_yz_zzzzz = pbuffer.data(idx_dip_dh + 230);

    auto tr_y_zz_xxxxx = pbuffer.data(idx_dip_dh + 231);

    auto tr_y_zz_xxxxy = pbuffer.data(idx_dip_dh + 232);

    auto tr_y_zz_xxxxz = pbuffer.data(idx_dip_dh + 233);

    auto tr_y_zz_xxxyy = pbuffer.data(idx_dip_dh + 234);

    auto tr_y_zz_xxxyz = pbuffer.data(idx_dip_dh + 235);

    auto tr_y_zz_xxxzz = pbuffer.data(idx_dip_dh + 236);

    auto tr_y_zz_xxyyy = pbuffer.data(idx_dip_dh + 237);

    auto tr_y_zz_xxyyz = pbuffer.data(idx_dip_dh + 238);

    auto tr_y_zz_xxyzz = pbuffer.data(idx_dip_dh + 239);

    auto tr_y_zz_xxzzz = pbuffer.data(idx_dip_dh + 240);

    auto tr_y_zz_xyyyy = pbuffer.data(idx_dip_dh + 241);

    auto tr_y_zz_xyyyz = pbuffer.data(idx_dip_dh + 242);

    auto tr_y_zz_xyyzz = pbuffer.data(idx_dip_dh + 243);

    auto tr_y_zz_xyzzz = pbuffer.data(idx_dip_dh + 244);

    auto tr_y_zz_xzzzz = pbuffer.data(idx_dip_dh + 245);

    auto tr_y_zz_yyyyy = pbuffer.data(idx_dip_dh + 246);

    auto tr_y_zz_yyyyz = pbuffer.data(idx_dip_dh + 247);

    auto tr_y_zz_yyyzz = pbuffer.data(idx_dip_dh + 248);

    auto tr_y_zz_yyzzz = pbuffer.data(idx_dip_dh + 249);

    auto tr_y_zz_yzzzz = pbuffer.data(idx_dip_dh + 250);

    auto tr_y_zz_zzzzz = pbuffer.data(idx_dip_dh + 251);

    auto tr_z_xx_xxxxx = pbuffer.data(idx_dip_dh + 252);

    auto tr_z_xx_xxxxy = pbuffer.data(idx_dip_dh + 253);

    auto tr_z_xx_xxxxz = pbuffer.data(idx_dip_dh + 254);

    auto tr_z_xx_xxxyy = pbuffer.data(idx_dip_dh + 255);

    auto tr_z_xx_xxxyz = pbuffer.data(idx_dip_dh + 256);

    auto tr_z_xx_xxxzz = pbuffer.data(idx_dip_dh + 257);

    auto tr_z_xx_xxyyy = pbuffer.data(idx_dip_dh + 258);

    auto tr_z_xx_xxyyz = pbuffer.data(idx_dip_dh + 259);

    auto tr_z_xx_xxyzz = pbuffer.data(idx_dip_dh + 260);

    auto tr_z_xx_xxzzz = pbuffer.data(idx_dip_dh + 261);

    auto tr_z_xx_xyyyy = pbuffer.data(idx_dip_dh + 262);

    auto tr_z_xx_xyyyz = pbuffer.data(idx_dip_dh + 263);

    auto tr_z_xx_xyyzz = pbuffer.data(idx_dip_dh + 264);

    auto tr_z_xx_xyzzz = pbuffer.data(idx_dip_dh + 265);

    auto tr_z_xx_xzzzz = pbuffer.data(idx_dip_dh + 266);

    auto tr_z_xx_yyyyy = pbuffer.data(idx_dip_dh + 267);

    auto tr_z_xx_yyyyz = pbuffer.data(idx_dip_dh + 268);

    auto tr_z_xx_yyyzz = pbuffer.data(idx_dip_dh + 269);

    auto tr_z_xx_yyzzz = pbuffer.data(idx_dip_dh + 270);

    auto tr_z_xx_yzzzz = pbuffer.data(idx_dip_dh + 271);

    auto tr_z_xx_zzzzz = pbuffer.data(idx_dip_dh + 272);

    auto tr_z_xy_yyyyy = pbuffer.data(idx_dip_dh + 288);

    auto tr_z_xy_yyyyz = pbuffer.data(idx_dip_dh + 289);

    auto tr_z_xy_yyyzz = pbuffer.data(idx_dip_dh + 290);

    auto tr_z_xy_yyzzz = pbuffer.data(idx_dip_dh + 291);

    auto tr_z_xy_yzzzz = pbuffer.data(idx_dip_dh + 292);

    auto tr_z_xz_xxxxz = pbuffer.data(idx_dip_dh + 296);

    auto tr_z_xz_xxxyz = pbuffer.data(idx_dip_dh + 298);

    auto tr_z_xz_xxxzz = pbuffer.data(idx_dip_dh + 299);

    auto tr_z_xz_xxyyz = pbuffer.data(idx_dip_dh + 301);

    auto tr_z_xz_xxyzz = pbuffer.data(idx_dip_dh + 302);

    auto tr_z_xz_xxzzz = pbuffer.data(idx_dip_dh + 303);

    auto tr_z_xz_xyyyz = pbuffer.data(idx_dip_dh + 305);

    auto tr_z_xz_xyyzz = pbuffer.data(idx_dip_dh + 306);

    auto tr_z_xz_xyzzz = pbuffer.data(idx_dip_dh + 307);

    auto tr_z_xz_xzzzz = pbuffer.data(idx_dip_dh + 308);

    auto tr_z_xz_yyyyy = pbuffer.data(idx_dip_dh + 309);

    auto tr_z_xz_yyyyz = pbuffer.data(idx_dip_dh + 310);

    auto tr_z_xz_yyyzz = pbuffer.data(idx_dip_dh + 311);

    auto tr_z_xz_yyzzz = pbuffer.data(idx_dip_dh + 312);

    auto tr_z_xz_yzzzz = pbuffer.data(idx_dip_dh + 313);

    auto tr_z_xz_zzzzz = pbuffer.data(idx_dip_dh + 314);

    auto tr_z_yy_xxxxx = pbuffer.data(idx_dip_dh + 315);

    auto tr_z_yy_xxxxy = pbuffer.data(idx_dip_dh + 316);

    auto tr_z_yy_xxxxz = pbuffer.data(idx_dip_dh + 317);

    auto tr_z_yy_xxxyy = pbuffer.data(idx_dip_dh + 318);

    auto tr_z_yy_xxxyz = pbuffer.data(idx_dip_dh + 319);

    auto tr_z_yy_xxxzz = pbuffer.data(idx_dip_dh + 320);

    auto tr_z_yy_xxyyy = pbuffer.data(idx_dip_dh + 321);

    auto tr_z_yy_xxyyz = pbuffer.data(idx_dip_dh + 322);

    auto tr_z_yy_xxyzz = pbuffer.data(idx_dip_dh + 323);

    auto tr_z_yy_xxzzz = pbuffer.data(idx_dip_dh + 324);

    auto tr_z_yy_xyyyy = pbuffer.data(idx_dip_dh + 325);

    auto tr_z_yy_xyyyz = pbuffer.data(idx_dip_dh + 326);

    auto tr_z_yy_xyyzz = pbuffer.data(idx_dip_dh + 327);

    auto tr_z_yy_xyzzz = pbuffer.data(idx_dip_dh + 328);

    auto tr_z_yy_xzzzz = pbuffer.data(idx_dip_dh + 329);

    auto tr_z_yy_yyyyy = pbuffer.data(idx_dip_dh + 330);

    auto tr_z_yy_yyyyz = pbuffer.data(idx_dip_dh + 331);

    auto tr_z_yy_yyyzz = pbuffer.data(idx_dip_dh + 332);

    auto tr_z_yy_yyzzz = pbuffer.data(idx_dip_dh + 333);

    auto tr_z_yy_yzzzz = pbuffer.data(idx_dip_dh + 334);

    auto tr_z_yy_zzzzz = pbuffer.data(idx_dip_dh + 335);

    auto tr_z_yz_xxxxx = pbuffer.data(idx_dip_dh + 336);

    auto tr_z_yz_xxxxz = pbuffer.data(idx_dip_dh + 338);

    auto tr_z_yz_xxxyz = pbuffer.data(idx_dip_dh + 340);

    auto tr_z_yz_xxxzz = pbuffer.data(idx_dip_dh + 341);

    auto tr_z_yz_xxyyz = pbuffer.data(idx_dip_dh + 343);

    auto tr_z_yz_xxyzz = pbuffer.data(idx_dip_dh + 344);

    auto tr_z_yz_xxzzz = pbuffer.data(idx_dip_dh + 345);

    auto tr_z_yz_xyyyz = pbuffer.data(idx_dip_dh + 347);

    auto tr_z_yz_xyyzz = pbuffer.data(idx_dip_dh + 348);

    auto tr_z_yz_xyzzz = pbuffer.data(idx_dip_dh + 349);

    auto tr_z_yz_xzzzz = pbuffer.data(idx_dip_dh + 350);

    auto tr_z_yz_yyyyy = pbuffer.data(idx_dip_dh + 351);

    auto tr_z_yz_yyyyz = pbuffer.data(idx_dip_dh + 352);

    auto tr_z_yz_yyyzz = pbuffer.data(idx_dip_dh + 353);

    auto tr_z_yz_yyzzz = pbuffer.data(idx_dip_dh + 354);

    auto tr_z_yz_yzzzz = pbuffer.data(idx_dip_dh + 355);

    auto tr_z_yz_zzzzz = pbuffer.data(idx_dip_dh + 356);

    auto tr_z_zz_xxxxx = pbuffer.data(idx_dip_dh + 357);

    auto tr_z_zz_xxxxy = pbuffer.data(idx_dip_dh + 358);

    auto tr_z_zz_xxxxz = pbuffer.data(idx_dip_dh + 359);

    auto tr_z_zz_xxxyy = pbuffer.data(idx_dip_dh + 360);

    auto tr_z_zz_xxxyz = pbuffer.data(idx_dip_dh + 361);

    auto tr_z_zz_xxxzz = pbuffer.data(idx_dip_dh + 362);

    auto tr_z_zz_xxyyy = pbuffer.data(idx_dip_dh + 363);

    auto tr_z_zz_xxyyz = pbuffer.data(idx_dip_dh + 364);

    auto tr_z_zz_xxyzz = pbuffer.data(idx_dip_dh + 365);

    auto tr_z_zz_xxzzz = pbuffer.data(idx_dip_dh + 366);

    auto tr_z_zz_xyyyy = pbuffer.data(idx_dip_dh + 367);

    auto tr_z_zz_xyyyz = pbuffer.data(idx_dip_dh + 368);

    auto tr_z_zz_xyyzz = pbuffer.data(idx_dip_dh + 369);

    auto tr_z_zz_xyzzz = pbuffer.data(idx_dip_dh + 370);

    auto tr_z_zz_xzzzz = pbuffer.data(idx_dip_dh + 371);

    auto tr_z_zz_yyyyy = pbuffer.data(idx_dip_dh + 372);

    auto tr_z_zz_yyyyz = pbuffer.data(idx_dip_dh + 373);

    auto tr_z_zz_yyyzz = pbuffer.data(idx_dip_dh + 374);

    auto tr_z_zz_yyzzz = pbuffer.data(idx_dip_dh + 375);

    auto tr_z_zz_yzzzz = pbuffer.data(idx_dip_dh + 376);

    auto tr_z_zz_zzzzz = pbuffer.data(idx_dip_dh + 377);

    // Set up components of auxiliary buffer : FG

    auto tr_x_xxx_xxxx = pbuffer.data(idx_dip_fg);

    auto tr_x_xxx_xxxy = pbuffer.data(idx_dip_fg + 1);

    auto tr_x_xxx_xxxz = pbuffer.data(idx_dip_fg + 2);

    auto tr_x_xxx_xxyy = pbuffer.data(idx_dip_fg + 3);

    auto tr_x_xxx_xxyz = pbuffer.data(idx_dip_fg + 4);

    auto tr_x_xxx_xxzz = pbuffer.data(idx_dip_fg + 5);

    auto tr_x_xxx_xyyy = pbuffer.data(idx_dip_fg + 6);

    auto tr_x_xxx_xyyz = pbuffer.data(idx_dip_fg + 7);

    auto tr_x_xxx_xyzz = pbuffer.data(idx_dip_fg + 8);

    auto tr_x_xxx_xzzz = pbuffer.data(idx_dip_fg + 9);

    auto tr_x_xxx_yyyy = pbuffer.data(idx_dip_fg + 10);

    auto tr_x_xxx_yyyz = pbuffer.data(idx_dip_fg + 11);

    auto tr_x_xxx_yyzz = pbuffer.data(idx_dip_fg + 12);

    auto tr_x_xxx_yzzz = pbuffer.data(idx_dip_fg + 13);

    auto tr_x_xxx_zzzz = pbuffer.data(idx_dip_fg + 14);

    auto tr_x_xxy_xxxx = pbuffer.data(idx_dip_fg + 15);

    auto tr_x_xxy_xxxy = pbuffer.data(idx_dip_fg + 16);

    auto tr_x_xxy_xxxz = pbuffer.data(idx_dip_fg + 17);

    auto tr_x_xxy_xxyy = pbuffer.data(idx_dip_fg + 18);

    auto tr_x_xxy_xxyz = pbuffer.data(idx_dip_fg + 19);

    auto tr_x_xxy_xxzz = pbuffer.data(idx_dip_fg + 20);

    auto tr_x_xxy_xyyy = pbuffer.data(idx_dip_fg + 21);

    auto tr_x_xxy_xyyz = pbuffer.data(idx_dip_fg + 22);

    auto tr_x_xxy_xyzz = pbuffer.data(idx_dip_fg + 23);

    auto tr_x_xxy_xzzz = pbuffer.data(idx_dip_fg + 24);

    auto tr_x_xxz_xxxx = pbuffer.data(idx_dip_fg + 30);

    auto tr_x_xxz_xxxy = pbuffer.data(idx_dip_fg + 31);

    auto tr_x_xxz_xxxz = pbuffer.data(idx_dip_fg + 32);

    auto tr_x_xxz_xxyy = pbuffer.data(idx_dip_fg + 33);

    auto tr_x_xxz_xxyz = pbuffer.data(idx_dip_fg + 34);

    auto tr_x_xxz_xxzz = pbuffer.data(idx_dip_fg + 35);

    auto tr_x_xxz_xyyy = pbuffer.data(idx_dip_fg + 36);

    auto tr_x_xxz_xyyz = pbuffer.data(idx_dip_fg + 37);

    auto tr_x_xxz_xyzz = pbuffer.data(idx_dip_fg + 38);

    auto tr_x_xxz_xzzz = pbuffer.data(idx_dip_fg + 39);

    auto tr_x_xxz_yyyz = pbuffer.data(idx_dip_fg + 41);

    auto tr_x_xxz_yyzz = pbuffer.data(idx_dip_fg + 42);

    auto tr_x_xxz_yzzz = pbuffer.data(idx_dip_fg + 43);

    auto tr_x_xxz_zzzz = pbuffer.data(idx_dip_fg + 44);

    auto tr_x_xyy_xxxy = pbuffer.data(idx_dip_fg + 46);

    auto tr_x_xyy_xxyy = pbuffer.data(idx_dip_fg + 48);

    auto tr_x_xyy_xxyz = pbuffer.data(idx_dip_fg + 49);

    auto tr_x_xyy_xyyy = pbuffer.data(idx_dip_fg + 51);

    auto tr_x_xyy_xyyz = pbuffer.data(idx_dip_fg + 52);

    auto tr_x_xyy_xyzz = pbuffer.data(idx_dip_fg + 53);

    auto tr_x_xzz_xxxx = pbuffer.data(idx_dip_fg + 75);

    auto tr_x_xzz_xxxy = pbuffer.data(idx_dip_fg + 76);

    auto tr_x_xzz_xxxz = pbuffer.data(idx_dip_fg + 77);

    auto tr_x_xzz_xxyy = pbuffer.data(idx_dip_fg + 78);

    auto tr_x_xzz_xxyz = pbuffer.data(idx_dip_fg + 79);

    auto tr_x_xzz_xxzz = pbuffer.data(idx_dip_fg + 80);

    auto tr_x_xzz_xyyy = pbuffer.data(idx_dip_fg + 81);

    auto tr_x_xzz_xyyz = pbuffer.data(idx_dip_fg + 82);

    auto tr_x_xzz_xyzz = pbuffer.data(idx_dip_fg + 83);

    auto tr_x_xzz_xzzz = pbuffer.data(idx_dip_fg + 84);

    auto tr_x_yyy_xxxx = pbuffer.data(idx_dip_fg + 90);

    auto tr_x_yyy_xxxy = pbuffer.data(idx_dip_fg + 91);

    auto tr_x_yyy_xxxz = pbuffer.data(idx_dip_fg + 92);

    auto tr_x_yyy_xxyy = pbuffer.data(idx_dip_fg + 93);

    auto tr_x_yyy_xxyz = pbuffer.data(idx_dip_fg + 94);

    auto tr_x_yyy_xxzz = pbuffer.data(idx_dip_fg + 95);

    auto tr_x_yyy_xyyy = pbuffer.data(idx_dip_fg + 96);

    auto tr_x_yyy_xyyz = pbuffer.data(idx_dip_fg + 97);

    auto tr_x_yyy_xyzz = pbuffer.data(idx_dip_fg + 98);

    auto tr_x_yyy_xzzz = pbuffer.data(idx_dip_fg + 99);

    auto tr_x_yyy_yyyy = pbuffer.data(idx_dip_fg + 100);

    auto tr_x_yyy_yyyz = pbuffer.data(idx_dip_fg + 101);

    auto tr_x_yyy_yyzz = pbuffer.data(idx_dip_fg + 102);

    auto tr_x_yyy_yzzz = pbuffer.data(idx_dip_fg + 103);

    auto tr_x_yyy_zzzz = pbuffer.data(idx_dip_fg + 104);

    auto tr_x_yzz_xxxz = pbuffer.data(idx_dip_fg + 122);

    auto tr_x_yzz_xxyz = pbuffer.data(idx_dip_fg + 124);

    auto tr_x_yzz_xxzz = pbuffer.data(idx_dip_fg + 125);

    auto tr_x_yzz_xyyz = pbuffer.data(idx_dip_fg + 127);

    auto tr_x_yzz_xyzz = pbuffer.data(idx_dip_fg + 128);

    auto tr_x_yzz_xzzz = pbuffer.data(idx_dip_fg + 129);

    auto tr_x_yzz_yyyz = pbuffer.data(idx_dip_fg + 131);

    auto tr_x_yzz_yyzz = pbuffer.data(idx_dip_fg + 132);

    auto tr_x_yzz_yzzz = pbuffer.data(idx_dip_fg + 133);

    auto tr_x_yzz_zzzz = pbuffer.data(idx_dip_fg + 134);

    auto tr_x_zzz_xxxx = pbuffer.data(idx_dip_fg + 135);

    auto tr_x_zzz_xxxy = pbuffer.data(idx_dip_fg + 136);

    auto tr_x_zzz_xxxz = pbuffer.data(idx_dip_fg + 137);

    auto tr_x_zzz_xxyy = pbuffer.data(idx_dip_fg + 138);

    auto tr_x_zzz_xxyz = pbuffer.data(idx_dip_fg + 139);

    auto tr_x_zzz_xxzz = pbuffer.data(idx_dip_fg + 140);

    auto tr_x_zzz_xyyy = pbuffer.data(idx_dip_fg + 141);

    auto tr_x_zzz_xyyz = pbuffer.data(idx_dip_fg + 142);

    auto tr_x_zzz_xyzz = pbuffer.data(idx_dip_fg + 143);

    auto tr_x_zzz_xzzz = pbuffer.data(idx_dip_fg + 144);

    auto tr_x_zzz_yyyy = pbuffer.data(idx_dip_fg + 145);

    auto tr_x_zzz_yyyz = pbuffer.data(idx_dip_fg + 146);

    auto tr_x_zzz_yyzz = pbuffer.data(idx_dip_fg + 147);

    auto tr_x_zzz_yzzz = pbuffer.data(idx_dip_fg + 148);

    auto tr_x_zzz_zzzz = pbuffer.data(idx_dip_fg + 149);

    auto tr_y_xxx_xxxx = pbuffer.data(idx_dip_fg + 150);

    auto tr_y_xxx_xxxy = pbuffer.data(idx_dip_fg + 151);

    auto tr_y_xxx_xxxz = pbuffer.data(idx_dip_fg + 152);

    auto tr_y_xxx_xxyy = pbuffer.data(idx_dip_fg + 153);

    auto tr_y_xxx_xxyz = pbuffer.data(idx_dip_fg + 154);

    auto tr_y_xxx_xxzz = pbuffer.data(idx_dip_fg + 155);

    auto tr_y_xxx_xyyy = pbuffer.data(idx_dip_fg + 156);

    auto tr_y_xxx_xyyz = pbuffer.data(idx_dip_fg + 157);

    auto tr_y_xxx_xyzz = pbuffer.data(idx_dip_fg + 158);

    auto tr_y_xxx_xzzz = pbuffer.data(idx_dip_fg + 159);

    auto tr_y_xxx_yyyy = pbuffer.data(idx_dip_fg + 160);

    auto tr_y_xxx_yyyz = pbuffer.data(idx_dip_fg + 161);

    auto tr_y_xxx_yyzz = pbuffer.data(idx_dip_fg + 162);

    auto tr_y_xxx_yzzz = pbuffer.data(idx_dip_fg + 163);

    auto tr_y_xxx_zzzz = pbuffer.data(idx_dip_fg + 164);

    auto tr_y_xxy_xxxy = pbuffer.data(idx_dip_fg + 166);

    auto tr_y_xxy_xxyy = pbuffer.data(idx_dip_fg + 168);

    auto tr_y_xxy_xxyz = pbuffer.data(idx_dip_fg + 169);

    auto tr_y_xxy_xyyy = pbuffer.data(idx_dip_fg + 171);

    auto tr_y_xxy_xyyz = pbuffer.data(idx_dip_fg + 172);

    auto tr_y_xxy_xyzz = pbuffer.data(idx_dip_fg + 173);

    auto tr_y_xxy_yyyy = pbuffer.data(idx_dip_fg + 175);

    auto tr_y_xxy_yyyz = pbuffer.data(idx_dip_fg + 176);

    auto tr_y_xxy_yyzz = pbuffer.data(idx_dip_fg + 177);

    auto tr_y_xxy_yzzz = pbuffer.data(idx_dip_fg + 178);

    auto tr_y_xyy_xxxx = pbuffer.data(idx_dip_fg + 195);

    auto tr_y_xyy_xxxy = pbuffer.data(idx_dip_fg + 196);

    auto tr_y_xyy_xxxz = pbuffer.data(idx_dip_fg + 197);

    auto tr_y_xyy_xxyy = pbuffer.data(idx_dip_fg + 198);

    auto tr_y_xyy_xxyz = pbuffer.data(idx_dip_fg + 199);

    auto tr_y_xyy_xxzz = pbuffer.data(idx_dip_fg + 200);

    auto tr_y_xyy_xyyy = pbuffer.data(idx_dip_fg + 201);

    auto tr_y_xyy_xyyz = pbuffer.data(idx_dip_fg + 202);

    auto tr_y_xyy_xyzz = pbuffer.data(idx_dip_fg + 203);

    auto tr_y_xyy_xzzz = pbuffer.data(idx_dip_fg + 204);

    auto tr_y_xyy_yyyy = pbuffer.data(idx_dip_fg + 205);

    auto tr_y_xyy_yyyz = pbuffer.data(idx_dip_fg + 206);

    auto tr_y_xyy_yyzz = pbuffer.data(idx_dip_fg + 207);

    auto tr_y_xyy_yzzz = pbuffer.data(idx_dip_fg + 208);

    auto tr_y_xyy_zzzz = pbuffer.data(idx_dip_fg + 209);

    auto tr_y_xzz_xxxz = pbuffer.data(idx_dip_fg + 227);

    auto tr_y_xzz_xxyz = pbuffer.data(idx_dip_fg + 229);

    auto tr_y_xzz_xxzz = pbuffer.data(idx_dip_fg + 230);

    auto tr_y_xzz_xyyz = pbuffer.data(idx_dip_fg + 232);

    auto tr_y_xzz_xyzz = pbuffer.data(idx_dip_fg + 233);

    auto tr_y_xzz_xzzz = pbuffer.data(idx_dip_fg + 234);

    auto tr_y_xzz_yyyz = pbuffer.data(idx_dip_fg + 236);

    auto tr_y_xzz_yyzz = pbuffer.data(idx_dip_fg + 237);

    auto tr_y_xzz_yzzz = pbuffer.data(idx_dip_fg + 238);

    auto tr_y_xzz_zzzz = pbuffer.data(idx_dip_fg + 239);

    auto tr_y_yyy_xxxx = pbuffer.data(idx_dip_fg + 240);

    auto tr_y_yyy_xxxy = pbuffer.data(idx_dip_fg + 241);

    auto tr_y_yyy_xxxz = pbuffer.data(idx_dip_fg + 242);

    auto tr_y_yyy_xxyy = pbuffer.data(idx_dip_fg + 243);

    auto tr_y_yyy_xxyz = pbuffer.data(idx_dip_fg + 244);

    auto tr_y_yyy_xxzz = pbuffer.data(idx_dip_fg + 245);

    auto tr_y_yyy_xyyy = pbuffer.data(idx_dip_fg + 246);

    auto tr_y_yyy_xyyz = pbuffer.data(idx_dip_fg + 247);

    auto tr_y_yyy_xyzz = pbuffer.data(idx_dip_fg + 248);

    auto tr_y_yyy_xzzz = pbuffer.data(idx_dip_fg + 249);

    auto tr_y_yyy_yyyy = pbuffer.data(idx_dip_fg + 250);

    auto tr_y_yyy_yyyz = pbuffer.data(idx_dip_fg + 251);

    auto tr_y_yyy_yyzz = pbuffer.data(idx_dip_fg + 252);

    auto tr_y_yyy_yzzz = pbuffer.data(idx_dip_fg + 253);

    auto tr_y_yyy_zzzz = pbuffer.data(idx_dip_fg + 254);

    auto tr_y_yyz_xxxy = pbuffer.data(idx_dip_fg + 256);

    auto tr_y_yyz_xxxz = pbuffer.data(idx_dip_fg + 257);

    auto tr_y_yyz_xxyy = pbuffer.data(idx_dip_fg + 258);

    auto tr_y_yyz_xxyz = pbuffer.data(idx_dip_fg + 259);

    auto tr_y_yyz_xxzz = pbuffer.data(idx_dip_fg + 260);

    auto tr_y_yyz_xyyy = pbuffer.data(idx_dip_fg + 261);

    auto tr_y_yyz_xyyz = pbuffer.data(idx_dip_fg + 262);

    auto tr_y_yyz_xyzz = pbuffer.data(idx_dip_fg + 263);

    auto tr_y_yyz_xzzz = pbuffer.data(idx_dip_fg + 264);

    auto tr_y_yyz_yyyy = pbuffer.data(idx_dip_fg + 265);

    auto tr_y_yyz_yyyz = pbuffer.data(idx_dip_fg + 266);

    auto tr_y_yyz_yyzz = pbuffer.data(idx_dip_fg + 267);

    auto tr_y_yyz_yzzz = pbuffer.data(idx_dip_fg + 268);

    auto tr_y_yyz_zzzz = pbuffer.data(idx_dip_fg + 269);

    auto tr_y_yzz_xxxx = pbuffer.data(idx_dip_fg + 270);

    auto tr_y_yzz_xxxy = pbuffer.data(idx_dip_fg + 271);

    auto tr_y_yzz_xxxz = pbuffer.data(idx_dip_fg + 272);

    auto tr_y_yzz_xxyy = pbuffer.data(idx_dip_fg + 273);

    auto tr_y_yzz_xxyz = pbuffer.data(idx_dip_fg + 274);

    auto tr_y_yzz_xxzz = pbuffer.data(idx_dip_fg + 275);

    auto tr_y_yzz_xyyy = pbuffer.data(idx_dip_fg + 276);

    auto tr_y_yzz_xyyz = pbuffer.data(idx_dip_fg + 277);

    auto tr_y_yzz_xyzz = pbuffer.data(idx_dip_fg + 278);

    auto tr_y_yzz_xzzz = pbuffer.data(idx_dip_fg + 279);

    auto tr_y_yzz_yyyy = pbuffer.data(idx_dip_fg + 280);

    auto tr_y_yzz_yyyz = pbuffer.data(idx_dip_fg + 281);

    auto tr_y_yzz_yyzz = pbuffer.data(idx_dip_fg + 282);

    auto tr_y_yzz_yzzz = pbuffer.data(idx_dip_fg + 283);

    auto tr_y_yzz_zzzz = pbuffer.data(idx_dip_fg + 284);

    auto tr_y_zzz_xxxx = pbuffer.data(idx_dip_fg + 285);

    auto tr_y_zzz_xxxy = pbuffer.data(idx_dip_fg + 286);

    auto tr_y_zzz_xxxz = pbuffer.data(idx_dip_fg + 287);

    auto tr_y_zzz_xxyy = pbuffer.data(idx_dip_fg + 288);

    auto tr_y_zzz_xxyz = pbuffer.data(idx_dip_fg + 289);

    auto tr_y_zzz_xxzz = pbuffer.data(idx_dip_fg + 290);

    auto tr_y_zzz_xyyy = pbuffer.data(idx_dip_fg + 291);

    auto tr_y_zzz_xyyz = pbuffer.data(idx_dip_fg + 292);

    auto tr_y_zzz_xyzz = pbuffer.data(idx_dip_fg + 293);

    auto tr_y_zzz_xzzz = pbuffer.data(idx_dip_fg + 294);

    auto tr_y_zzz_yyyy = pbuffer.data(idx_dip_fg + 295);

    auto tr_y_zzz_yyyz = pbuffer.data(idx_dip_fg + 296);

    auto tr_y_zzz_yyzz = pbuffer.data(idx_dip_fg + 297);

    auto tr_y_zzz_yzzz = pbuffer.data(idx_dip_fg + 298);

    auto tr_y_zzz_zzzz = pbuffer.data(idx_dip_fg + 299);

    auto tr_z_xxx_xxxx = pbuffer.data(idx_dip_fg + 300);

    auto tr_z_xxx_xxxy = pbuffer.data(idx_dip_fg + 301);

    auto tr_z_xxx_xxxz = pbuffer.data(idx_dip_fg + 302);

    auto tr_z_xxx_xxyy = pbuffer.data(idx_dip_fg + 303);

    auto tr_z_xxx_xxyz = pbuffer.data(idx_dip_fg + 304);

    auto tr_z_xxx_xxzz = pbuffer.data(idx_dip_fg + 305);

    auto tr_z_xxx_xyyy = pbuffer.data(idx_dip_fg + 306);

    auto tr_z_xxx_xyyz = pbuffer.data(idx_dip_fg + 307);

    auto tr_z_xxx_xyzz = pbuffer.data(idx_dip_fg + 308);

    auto tr_z_xxx_xzzz = pbuffer.data(idx_dip_fg + 309);

    auto tr_z_xxx_yyyy = pbuffer.data(idx_dip_fg + 310);

    auto tr_z_xxx_yyyz = pbuffer.data(idx_dip_fg + 311);

    auto tr_z_xxx_yyzz = pbuffer.data(idx_dip_fg + 312);

    auto tr_z_xxx_yzzz = pbuffer.data(idx_dip_fg + 313);

    auto tr_z_xxx_zzzz = pbuffer.data(idx_dip_fg + 314);

    auto tr_z_xxz_xxxx = pbuffer.data(idx_dip_fg + 330);

    auto tr_z_xxz_xxxy = pbuffer.data(idx_dip_fg + 331);

    auto tr_z_xxz_xxxz = pbuffer.data(idx_dip_fg + 332);

    auto tr_z_xxz_xxyy = pbuffer.data(idx_dip_fg + 333);

    auto tr_z_xxz_xxyz = pbuffer.data(idx_dip_fg + 334);

    auto tr_z_xxz_xxzz = pbuffer.data(idx_dip_fg + 335);

    auto tr_z_xxz_xyyy = pbuffer.data(idx_dip_fg + 336);

    auto tr_z_xxz_xyyz = pbuffer.data(idx_dip_fg + 337);

    auto tr_z_xxz_xyzz = pbuffer.data(idx_dip_fg + 338);

    auto tr_z_xxz_xzzz = pbuffer.data(idx_dip_fg + 339);

    auto tr_z_xxz_yyyz = pbuffer.data(idx_dip_fg + 341);

    auto tr_z_xxz_yyzz = pbuffer.data(idx_dip_fg + 342);

    auto tr_z_xxz_yzzz = pbuffer.data(idx_dip_fg + 343);

    auto tr_z_xxz_zzzz = pbuffer.data(idx_dip_fg + 344);

    auto tr_z_xyy_xxxy = pbuffer.data(idx_dip_fg + 346);

    auto tr_z_xyy_xxyy = pbuffer.data(idx_dip_fg + 348);

    auto tr_z_xyy_xxyz = pbuffer.data(idx_dip_fg + 349);

    auto tr_z_xyy_xyyy = pbuffer.data(idx_dip_fg + 351);

    auto tr_z_xyy_xyyz = pbuffer.data(idx_dip_fg + 352);

    auto tr_z_xyy_xyzz = pbuffer.data(idx_dip_fg + 353);

    auto tr_z_xyy_yyyy = pbuffer.data(idx_dip_fg + 355);

    auto tr_z_xyy_yyyz = pbuffer.data(idx_dip_fg + 356);

    auto tr_z_xyy_yyzz = pbuffer.data(idx_dip_fg + 357);

    auto tr_z_xyy_yzzz = pbuffer.data(idx_dip_fg + 358);

    auto tr_z_xzz_xxxx = pbuffer.data(idx_dip_fg + 375);

    auto tr_z_xzz_xxxy = pbuffer.data(idx_dip_fg + 376);

    auto tr_z_xzz_xxxz = pbuffer.data(idx_dip_fg + 377);

    auto tr_z_xzz_xxyy = pbuffer.data(idx_dip_fg + 378);

    auto tr_z_xzz_xxyz = pbuffer.data(idx_dip_fg + 379);

    auto tr_z_xzz_xxzz = pbuffer.data(idx_dip_fg + 380);

    auto tr_z_xzz_xyyy = pbuffer.data(idx_dip_fg + 381);

    auto tr_z_xzz_xyyz = pbuffer.data(idx_dip_fg + 382);

    auto tr_z_xzz_xyzz = pbuffer.data(idx_dip_fg + 383);

    auto tr_z_xzz_xzzz = pbuffer.data(idx_dip_fg + 384);

    auto tr_z_xzz_yyyy = pbuffer.data(idx_dip_fg + 385);

    auto tr_z_xzz_yyyz = pbuffer.data(idx_dip_fg + 386);

    auto tr_z_xzz_yyzz = pbuffer.data(idx_dip_fg + 387);

    auto tr_z_xzz_yzzz = pbuffer.data(idx_dip_fg + 388);

    auto tr_z_xzz_zzzz = pbuffer.data(idx_dip_fg + 389);

    auto tr_z_yyy_xxxx = pbuffer.data(idx_dip_fg + 390);

    auto tr_z_yyy_xxxy = pbuffer.data(idx_dip_fg + 391);

    auto tr_z_yyy_xxxz = pbuffer.data(idx_dip_fg + 392);

    auto tr_z_yyy_xxyy = pbuffer.data(idx_dip_fg + 393);

    auto tr_z_yyy_xxyz = pbuffer.data(idx_dip_fg + 394);

    auto tr_z_yyy_xxzz = pbuffer.data(idx_dip_fg + 395);

    auto tr_z_yyy_xyyy = pbuffer.data(idx_dip_fg + 396);

    auto tr_z_yyy_xyyz = pbuffer.data(idx_dip_fg + 397);

    auto tr_z_yyy_xyzz = pbuffer.data(idx_dip_fg + 398);

    auto tr_z_yyy_xzzz = pbuffer.data(idx_dip_fg + 399);

    auto tr_z_yyy_yyyy = pbuffer.data(idx_dip_fg + 400);

    auto tr_z_yyy_yyyz = pbuffer.data(idx_dip_fg + 401);

    auto tr_z_yyy_yyzz = pbuffer.data(idx_dip_fg + 402);

    auto tr_z_yyy_yzzz = pbuffer.data(idx_dip_fg + 403);

    auto tr_z_yyy_zzzz = pbuffer.data(idx_dip_fg + 404);

    auto tr_z_yyz_xxxx = pbuffer.data(idx_dip_fg + 405);

    auto tr_z_yyz_xxxy = pbuffer.data(idx_dip_fg + 406);

    auto tr_z_yyz_xxxz = pbuffer.data(idx_dip_fg + 407);

    auto tr_z_yyz_xxyy = pbuffer.data(idx_dip_fg + 408);

    auto tr_z_yyz_xxyz = pbuffer.data(idx_dip_fg + 409);

    auto tr_z_yyz_xxzz = pbuffer.data(idx_dip_fg + 410);

    auto tr_z_yyz_xyyy = pbuffer.data(idx_dip_fg + 411);

    auto tr_z_yyz_xyyz = pbuffer.data(idx_dip_fg + 412);

    auto tr_z_yyz_xyzz = pbuffer.data(idx_dip_fg + 413);

    auto tr_z_yyz_xzzz = pbuffer.data(idx_dip_fg + 414);

    auto tr_z_yyz_yyyy = pbuffer.data(idx_dip_fg + 415);

    auto tr_z_yyz_yyyz = pbuffer.data(idx_dip_fg + 416);

    auto tr_z_yyz_yyzz = pbuffer.data(idx_dip_fg + 417);

    auto tr_z_yyz_yzzz = pbuffer.data(idx_dip_fg + 418);

    auto tr_z_yyz_zzzz = pbuffer.data(idx_dip_fg + 419);

    auto tr_z_yzz_xxxx = pbuffer.data(idx_dip_fg + 420);

    auto tr_z_yzz_xxxy = pbuffer.data(idx_dip_fg + 421);

    auto tr_z_yzz_xxxz = pbuffer.data(idx_dip_fg + 422);

    auto tr_z_yzz_xxyy = pbuffer.data(idx_dip_fg + 423);

    auto tr_z_yzz_xxyz = pbuffer.data(idx_dip_fg + 424);

    auto tr_z_yzz_xxzz = pbuffer.data(idx_dip_fg + 425);

    auto tr_z_yzz_xyyy = pbuffer.data(idx_dip_fg + 426);

    auto tr_z_yzz_xyyz = pbuffer.data(idx_dip_fg + 427);

    auto tr_z_yzz_xyzz = pbuffer.data(idx_dip_fg + 428);

    auto tr_z_yzz_xzzz = pbuffer.data(idx_dip_fg + 429);

    auto tr_z_yzz_yyyy = pbuffer.data(idx_dip_fg + 430);

    auto tr_z_yzz_yyyz = pbuffer.data(idx_dip_fg + 431);

    auto tr_z_yzz_yyzz = pbuffer.data(idx_dip_fg + 432);

    auto tr_z_yzz_yzzz = pbuffer.data(idx_dip_fg + 433);

    auto tr_z_yzz_zzzz = pbuffer.data(idx_dip_fg + 434);

    auto tr_z_zzz_xxxx = pbuffer.data(idx_dip_fg + 435);

    auto tr_z_zzz_xxxy = pbuffer.data(idx_dip_fg + 436);

    auto tr_z_zzz_xxxz = pbuffer.data(idx_dip_fg + 437);

    auto tr_z_zzz_xxyy = pbuffer.data(idx_dip_fg + 438);

    auto tr_z_zzz_xxyz = pbuffer.data(idx_dip_fg + 439);

    auto tr_z_zzz_xxzz = pbuffer.data(idx_dip_fg + 440);

    auto tr_z_zzz_xyyy = pbuffer.data(idx_dip_fg + 441);

    auto tr_z_zzz_xyyz = pbuffer.data(idx_dip_fg + 442);

    auto tr_z_zzz_xyzz = pbuffer.data(idx_dip_fg + 443);

    auto tr_z_zzz_xzzz = pbuffer.data(idx_dip_fg + 444);

    auto tr_z_zzz_yyyy = pbuffer.data(idx_dip_fg + 445);

    auto tr_z_zzz_yyyz = pbuffer.data(idx_dip_fg + 446);

    auto tr_z_zzz_yyzz = pbuffer.data(idx_dip_fg + 447);

    auto tr_z_zzz_yzzz = pbuffer.data(idx_dip_fg + 448);

    auto tr_z_zzz_zzzz = pbuffer.data(idx_dip_fg + 449);

    // Set up components of auxiliary buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_ovl_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_ovl_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_ovl_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_ovl_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_ovl_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_ovl_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_ovl_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_ovl_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_ovl_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_ovl_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_ovl_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_ovl_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_ovl_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_ovl_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_ovl_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_ovl_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_ovl_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_ovl_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_ovl_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_ovl_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_ovl_fh + 20);

    auto ts_xxz_xxxxz = pbuffer.data(idx_ovl_fh + 44);

    auto ts_xxz_xxxzz = pbuffer.data(idx_ovl_fh + 47);

    auto ts_xxz_xxzzz = pbuffer.data(idx_ovl_fh + 51);

    auto ts_xxz_xzzzz = pbuffer.data(idx_ovl_fh + 56);

    auto ts_xyy_yyyyy = pbuffer.data(idx_ovl_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_ovl_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_ovl_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_ovl_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_ovl_fh + 82);

    auto ts_xzz_yyyyz = pbuffer.data(idx_ovl_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_ovl_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_ovl_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_ovl_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_ovl_fh + 125);

    auto ts_yyy_xxxxx = pbuffer.data(idx_ovl_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_ovl_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_ovl_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_ovl_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_ovl_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_ovl_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_ovl_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_ovl_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_ovl_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_ovl_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_ovl_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_ovl_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_ovl_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_ovl_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_ovl_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_ovl_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_ovl_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_ovl_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_ovl_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_ovl_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_ovl_fh + 146);

    auto ts_yyz_yyyyz = pbuffer.data(idx_ovl_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_ovl_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_ovl_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_ovl_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_ovl_fh + 167);

    auto ts_yzz_xxxxz = pbuffer.data(idx_ovl_fh + 170);

    auto ts_yzz_xxxzz = pbuffer.data(idx_ovl_fh + 173);

    auto ts_yzz_xxzzz = pbuffer.data(idx_ovl_fh + 177);

    auto ts_yzz_xzzzz = pbuffer.data(idx_ovl_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_ovl_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_ovl_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_ovl_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_ovl_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_ovl_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_ovl_fh + 188);

    auto ts_zzz_xxxxx = pbuffer.data(idx_ovl_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_ovl_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_ovl_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_ovl_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_ovl_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_ovl_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_ovl_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_ovl_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_ovl_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_ovl_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_ovl_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_ovl_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_ovl_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_ovl_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_ovl_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_ovl_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_ovl_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_ovl_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_ovl_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_ovl_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_ovl_fh + 209);

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

    auto tr_x_xxy_yyyyy = pbuffer.data(idx_dip_fh + 36);

    auto tr_x_xxy_zzzzz = pbuffer.data(idx_dip_fh + 41);

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

    auto tr_x_xxz_yyyyy = pbuffer.data(idx_dip_fh + 57);

    auto tr_x_xxz_yyyyz = pbuffer.data(idx_dip_fh + 58);

    auto tr_x_xxz_yyyzz = pbuffer.data(idx_dip_fh + 59);

    auto tr_x_xxz_yyzzz = pbuffer.data(idx_dip_fh + 60);

    auto tr_x_xxz_yzzzz = pbuffer.data(idx_dip_fh + 61);

    auto tr_x_xxz_zzzzz = pbuffer.data(idx_dip_fh + 62);

    auto tr_x_xyy_xxxxx = pbuffer.data(idx_dip_fh + 63);

    auto tr_x_xyy_xxxxy = pbuffer.data(idx_dip_fh + 64);

    auto tr_x_xyy_xxxxz = pbuffer.data(idx_dip_fh + 65);

    auto tr_x_xyy_xxxyy = pbuffer.data(idx_dip_fh + 66);

    auto tr_x_xyy_xxxyz = pbuffer.data(idx_dip_fh + 67);

    auto tr_x_xyy_xxxzz = pbuffer.data(idx_dip_fh + 68);

    auto tr_x_xyy_xxyyy = pbuffer.data(idx_dip_fh + 69);

    auto tr_x_xyy_xxyyz = pbuffer.data(idx_dip_fh + 70);

    auto tr_x_xyy_xxyzz = pbuffer.data(idx_dip_fh + 71);

    auto tr_x_xyy_xxzzz = pbuffer.data(idx_dip_fh + 72);

    auto tr_x_xyy_xyyyy = pbuffer.data(idx_dip_fh + 73);

    auto tr_x_xyy_xyyyz = pbuffer.data(idx_dip_fh + 74);

    auto tr_x_xyy_xyyzz = pbuffer.data(idx_dip_fh + 75);

    auto tr_x_xyy_xyzzz = pbuffer.data(idx_dip_fh + 76);

    auto tr_x_xyy_xzzzz = pbuffer.data(idx_dip_fh + 77);

    auto tr_x_xyy_yyyyy = pbuffer.data(idx_dip_fh + 78);

    auto tr_x_xyy_yyyyz = pbuffer.data(idx_dip_fh + 79);

    auto tr_x_xyy_yyyzz = pbuffer.data(idx_dip_fh + 80);

    auto tr_x_xyy_yyzzz = pbuffer.data(idx_dip_fh + 81);

    auto tr_x_xyy_yzzzz = pbuffer.data(idx_dip_fh + 82);

    auto tr_x_xyz_xxxxz = pbuffer.data(idx_dip_fh + 86);

    auto tr_x_xyz_xxxzz = pbuffer.data(idx_dip_fh + 89);

    auto tr_x_xyz_xxzzz = pbuffer.data(idx_dip_fh + 93);

    auto tr_x_xyz_xzzzz = pbuffer.data(idx_dip_fh + 98);

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

    auto tr_x_xzz_yyyyz = pbuffer.data(idx_dip_fh + 121);

    auto tr_x_xzz_yyyzz = pbuffer.data(idx_dip_fh + 122);

    auto tr_x_xzz_yyzzz = pbuffer.data(idx_dip_fh + 123);

    auto tr_x_xzz_yzzzz = pbuffer.data(idx_dip_fh + 124);

    auto tr_x_xzz_zzzzz = pbuffer.data(idx_dip_fh + 125);

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

    auto tr_x_yyz_xxxxy = pbuffer.data(idx_dip_fh + 148);

    auto tr_x_yyz_xxxxz = pbuffer.data(idx_dip_fh + 149);

    auto tr_x_yyz_xxxyy = pbuffer.data(idx_dip_fh + 150);

    auto tr_x_yyz_xxxzz = pbuffer.data(idx_dip_fh + 152);

    auto tr_x_yyz_xxyyy = pbuffer.data(idx_dip_fh + 153);

    auto tr_x_yyz_xxzzz = pbuffer.data(idx_dip_fh + 156);

    auto tr_x_yyz_xyyyy = pbuffer.data(idx_dip_fh + 157);

    auto tr_x_yyz_xzzzz = pbuffer.data(idx_dip_fh + 161);

    auto tr_x_yyz_yyyyy = pbuffer.data(idx_dip_fh + 162);

    auto tr_x_yyz_yyyyz = pbuffer.data(idx_dip_fh + 163);

    auto tr_x_yyz_yyyzz = pbuffer.data(idx_dip_fh + 164);

    auto tr_x_yyz_yyzzz = pbuffer.data(idx_dip_fh + 165);

    auto tr_x_yyz_yzzzz = pbuffer.data(idx_dip_fh + 166);

    auto tr_x_yyz_zzzzz = pbuffer.data(idx_dip_fh + 167);

    auto tr_x_yzz_xxxxx = pbuffer.data(idx_dip_fh + 168);

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

    auto tr_x_yzz_yyyyy = pbuffer.data(idx_dip_fh + 183);

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

    auto tr_y_xxy_xxxxx = pbuffer.data(idx_dip_fh + 231);

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

    auto tr_y_xxy_zzzzz = pbuffer.data(idx_dip_fh + 251);

    auto tr_y_xxz_xxxxx = pbuffer.data(idx_dip_fh + 252);

    auto tr_y_xxz_xxxxy = pbuffer.data(idx_dip_fh + 253);

    auto tr_y_xxz_xxxxz = pbuffer.data(idx_dip_fh + 254);

    auto tr_y_xxz_xxxyy = pbuffer.data(idx_dip_fh + 255);

    auto tr_y_xxz_xxxzz = pbuffer.data(idx_dip_fh + 257);

    auto tr_y_xxz_xxyyy = pbuffer.data(idx_dip_fh + 258);

    auto tr_y_xxz_xxzzz = pbuffer.data(idx_dip_fh + 261);

    auto tr_y_xxz_xyyyy = pbuffer.data(idx_dip_fh + 262);

    auto tr_y_xxz_xzzzz = pbuffer.data(idx_dip_fh + 266);

    auto tr_y_xxz_yyyyz = pbuffer.data(idx_dip_fh + 268);

    auto tr_y_xxz_yyyzz = pbuffer.data(idx_dip_fh + 269);

    auto tr_y_xxz_yyzzz = pbuffer.data(idx_dip_fh + 270);

    auto tr_y_xxz_yzzzz = pbuffer.data(idx_dip_fh + 271);

    auto tr_y_xxz_zzzzz = pbuffer.data(idx_dip_fh + 272);

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

    auto tr_y_xyz_yyyyz = pbuffer.data(idx_dip_fh + 310);

    auto tr_y_xyz_yyyzz = pbuffer.data(idx_dip_fh + 311);

    auto tr_y_xyz_yyzzz = pbuffer.data(idx_dip_fh + 312);

    auto tr_y_xyz_yzzzz = pbuffer.data(idx_dip_fh + 313);

    auto tr_y_xyz_zzzzz = pbuffer.data(idx_dip_fh + 314);

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

    auto tr_y_xzz_yyyyy = pbuffer.data(idx_dip_fh + 330);

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

    auto tr_y_yyz_xxxxx = pbuffer.data(idx_dip_fh + 357);

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

    auto tr_z_xxy_xxxxx = pbuffer.data(idx_dip_fh + 441);

    auto tr_z_xxy_xxxxz = pbuffer.data(idx_dip_fh + 443);

    auto tr_z_xxy_xxxzz = pbuffer.data(idx_dip_fh + 446);

    auto tr_z_xxy_xxzzz = pbuffer.data(idx_dip_fh + 450);

    auto tr_z_xxy_xzzzz = pbuffer.data(idx_dip_fh + 455);

    auto tr_z_xxy_yyyyy = pbuffer.data(idx_dip_fh + 456);

    auto tr_z_xxy_yyyyz = pbuffer.data(idx_dip_fh + 457);

    auto tr_z_xxy_yyyzz = pbuffer.data(idx_dip_fh + 458);

    auto tr_z_xxy_yyzzz = pbuffer.data(idx_dip_fh + 459);

    auto tr_z_xxy_yzzzz = pbuffer.data(idx_dip_fh + 460);

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

    auto tr_z_xxz_yyyyy = pbuffer.data(idx_dip_fh + 477);

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

    auto tr_z_xyy_zzzzz = pbuffer.data(idx_dip_fh + 503);

    auto tr_z_xyz_yyyyy = pbuffer.data(idx_dip_fh + 519);

    auto tr_z_xyz_yyyyz = pbuffer.data(idx_dip_fh + 520);

    auto tr_z_xyz_yyyzz = pbuffer.data(idx_dip_fh + 521);

    auto tr_z_xyz_yyzzz = pbuffer.data(idx_dip_fh + 522);

    auto tr_z_xyz_yzzzz = pbuffer.data(idx_dip_fh + 523);

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

    // Set up 0-21 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_x_xx_xxxxx, tr_x_xx_xxxxy, tr_x_xx_xxxxz, tr_x_xx_xxxyy, tr_x_xx_xxxyz, tr_x_xx_xxxzz, tr_x_xx_xxyyy, tr_x_xx_xxyyz, tr_x_xx_xxyzz, tr_x_xx_xxzzz, tr_x_xx_xyyyy, tr_x_xx_xyyyz, tr_x_xx_xyyzz, tr_x_xx_xyzzz, tr_x_xx_xzzzz, tr_x_xx_yyyyy, tr_x_xx_yyyyz, tr_x_xx_yyyzz, tr_x_xx_yyzzz, tr_x_xx_yzzzz, tr_x_xx_zzzzz, tr_x_xxx_xxxx, tr_x_xxx_xxxxx, tr_x_xxx_xxxxy, tr_x_xxx_xxxxz, tr_x_xxx_xxxy, tr_x_xxx_xxxyy, tr_x_xxx_xxxyz, tr_x_xxx_xxxz, tr_x_xxx_xxxzz, tr_x_xxx_xxyy, tr_x_xxx_xxyyy, tr_x_xxx_xxyyz, tr_x_xxx_xxyz, tr_x_xxx_xxyzz, tr_x_xxx_xxzz, tr_x_xxx_xxzzz, tr_x_xxx_xyyy, tr_x_xxx_xyyyy, tr_x_xxx_xyyyz, tr_x_xxx_xyyz, tr_x_xxx_xyyzz, tr_x_xxx_xyzz, tr_x_xxx_xyzzz, tr_x_xxx_xzzz, tr_x_xxx_xzzzz, tr_x_xxx_yyyy, tr_x_xxx_yyyyy, tr_x_xxx_yyyyz, tr_x_xxx_yyyz, tr_x_xxx_yyyzz, tr_x_xxx_yyzz, tr_x_xxx_yyzzz, tr_x_xxx_yzzz, tr_x_xxx_yzzzz, tr_x_xxx_zzzz, tr_x_xxx_zzzzz, tr_x_xxxx_xxxxx, tr_x_xxxx_xxxxy, tr_x_xxxx_xxxxz, tr_x_xxxx_xxxyy, tr_x_xxxx_xxxyz, tr_x_xxxx_xxxzz, tr_x_xxxx_xxyyy, tr_x_xxxx_xxyyz, tr_x_xxxx_xxyzz, tr_x_xxxx_xxzzz, tr_x_xxxx_xyyyy, tr_x_xxxx_xyyyz, tr_x_xxxx_xyyzz, tr_x_xxxx_xyzzz, tr_x_xxxx_xzzzz, tr_x_xxxx_yyyyy, tr_x_xxxx_yyyyz, tr_x_xxxx_yyyzz, tr_x_xxxx_yyzzz, tr_x_xxxx_yzzzz, tr_x_xxxx_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxzz, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyzz, ts_xxx_xxzzz, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyzz, ts_xxx_xyzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyzz, ts_xxx_yyzzz, ts_xxx_yzzzz, ts_xxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_xxxxx[i] = 3.0 * tr_x_xx_xxxxx[i] * fe_0 + 5.0 * tr_x_xxx_xxxx[i] * fe_0 + ts_xxx_xxxxx[i] * fe_0 + tr_x_xxx_xxxxx[i] * pa_x[i];

        tr_x_xxxx_xxxxy[i] = 3.0 * tr_x_xx_xxxxy[i] * fe_0 + 4.0 * tr_x_xxx_xxxy[i] * fe_0 + ts_xxx_xxxxy[i] * fe_0 + tr_x_xxx_xxxxy[i] * pa_x[i];

        tr_x_xxxx_xxxxz[i] = 3.0 * tr_x_xx_xxxxz[i] * fe_0 + 4.0 * tr_x_xxx_xxxz[i] * fe_0 + ts_xxx_xxxxz[i] * fe_0 + tr_x_xxx_xxxxz[i] * pa_x[i];

        tr_x_xxxx_xxxyy[i] = 3.0 * tr_x_xx_xxxyy[i] * fe_0 + 3.0 * tr_x_xxx_xxyy[i] * fe_0 + ts_xxx_xxxyy[i] * fe_0 + tr_x_xxx_xxxyy[i] * pa_x[i];

        tr_x_xxxx_xxxyz[i] = 3.0 * tr_x_xx_xxxyz[i] * fe_0 + 3.0 * tr_x_xxx_xxyz[i] * fe_0 + ts_xxx_xxxyz[i] * fe_0 + tr_x_xxx_xxxyz[i] * pa_x[i];

        tr_x_xxxx_xxxzz[i] = 3.0 * tr_x_xx_xxxzz[i] * fe_0 + 3.0 * tr_x_xxx_xxzz[i] * fe_0 + ts_xxx_xxxzz[i] * fe_0 + tr_x_xxx_xxxzz[i] * pa_x[i];

        tr_x_xxxx_xxyyy[i] = 3.0 * tr_x_xx_xxyyy[i] * fe_0 + 2.0 * tr_x_xxx_xyyy[i] * fe_0 + ts_xxx_xxyyy[i] * fe_0 + tr_x_xxx_xxyyy[i] * pa_x[i];

        tr_x_xxxx_xxyyz[i] = 3.0 * tr_x_xx_xxyyz[i] * fe_0 + 2.0 * tr_x_xxx_xyyz[i] * fe_0 + ts_xxx_xxyyz[i] * fe_0 + tr_x_xxx_xxyyz[i] * pa_x[i];

        tr_x_xxxx_xxyzz[i] = 3.0 * tr_x_xx_xxyzz[i] * fe_0 + 2.0 * tr_x_xxx_xyzz[i] * fe_0 + ts_xxx_xxyzz[i] * fe_0 + tr_x_xxx_xxyzz[i] * pa_x[i];

        tr_x_xxxx_xxzzz[i] = 3.0 * tr_x_xx_xxzzz[i] * fe_0 + 2.0 * tr_x_xxx_xzzz[i] * fe_0 + ts_xxx_xxzzz[i] * fe_0 + tr_x_xxx_xxzzz[i] * pa_x[i];

        tr_x_xxxx_xyyyy[i] = 3.0 * tr_x_xx_xyyyy[i] * fe_0 + tr_x_xxx_yyyy[i] * fe_0 + ts_xxx_xyyyy[i] * fe_0 + tr_x_xxx_xyyyy[i] * pa_x[i];

        tr_x_xxxx_xyyyz[i] = 3.0 * tr_x_xx_xyyyz[i] * fe_0 + tr_x_xxx_yyyz[i] * fe_0 + ts_xxx_xyyyz[i] * fe_0 + tr_x_xxx_xyyyz[i] * pa_x[i];

        tr_x_xxxx_xyyzz[i] = 3.0 * tr_x_xx_xyyzz[i] * fe_0 + tr_x_xxx_yyzz[i] * fe_0 + ts_xxx_xyyzz[i] * fe_0 + tr_x_xxx_xyyzz[i] * pa_x[i];

        tr_x_xxxx_xyzzz[i] = 3.0 * tr_x_xx_xyzzz[i] * fe_0 + tr_x_xxx_yzzz[i] * fe_0 + ts_xxx_xyzzz[i] * fe_0 + tr_x_xxx_xyzzz[i] * pa_x[i];

        tr_x_xxxx_xzzzz[i] = 3.0 * tr_x_xx_xzzzz[i] * fe_0 + tr_x_xxx_zzzz[i] * fe_0 + ts_xxx_xzzzz[i] * fe_0 + tr_x_xxx_xzzzz[i] * pa_x[i];

        tr_x_xxxx_yyyyy[i] = 3.0 * tr_x_xx_yyyyy[i] * fe_0 + ts_xxx_yyyyy[i] * fe_0 + tr_x_xxx_yyyyy[i] * pa_x[i];

        tr_x_xxxx_yyyyz[i] = 3.0 * tr_x_xx_yyyyz[i] * fe_0 + ts_xxx_yyyyz[i] * fe_0 + tr_x_xxx_yyyyz[i] * pa_x[i];

        tr_x_xxxx_yyyzz[i] = 3.0 * tr_x_xx_yyyzz[i] * fe_0 + ts_xxx_yyyzz[i] * fe_0 + tr_x_xxx_yyyzz[i] * pa_x[i];

        tr_x_xxxx_yyzzz[i] = 3.0 * tr_x_xx_yyzzz[i] * fe_0 + ts_xxx_yyzzz[i] * fe_0 + tr_x_xxx_yyzzz[i] * pa_x[i];

        tr_x_xxxx_yzzzz[i] = 3.0 * tr_x_xx_yzzzz[i] * fe_0 + ts_xxx_yzzzz[i] * fe_0 + tr_x_xxx_yzzzz[i] * pa_x[i];

        tr_x_xxxx_zzzzz[i] = 3.0 * tr_x_xx_zzzzz[i] * fe_0 + ts_xxx_zzzzz[i] * fe_0 + tr_x_xxx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : GH

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

    auto tr_x_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 36);

    auto tr_x_xxxy_yyyyz = pbuffer.data(idx_dip_gh + 37);

    auto tr_x_xxxy_yyyzz = pbuffer.data(idx_dip_gh + 38);

    auto tr_x_xxxy_yyzzz = pbuffer.data(idx_dip_gh + 39);

    auto tr_x_xxxy_yzzzz = pbuffer.data(idx_dip_gh + 40);

    auto tr_x_xxxy_zzzzz = pbuffer.data(idx_dip_gh + 41);

    #pragma omp simd aligned(pa_y, tr_x_xxx_xxxx, tr_x_xxx_xxxxx, tr_x_xxx_xxxxy, tr_x_xxx_xxxxz, tr_x_xxx_xxxy, tr_x_xxx_xxxyy, tr_x_xxx_xxxyz, tr_x_xxx_xxxz, tr_x_xxx_xxxzz, tr_x_xxx_xxyy, tr_x_xxx_xxyyy, tr_x_xxx_xxyyz, tr_x_xxx_xxyz, tr_x_xxx_xxyzz, tr_x_xxx_xxzz, tr_x_xxx_xxzzz, tr_x_xxx_xyyy, tr_x_xxx_xyyyy, tr_x_xxx_xyyyz, tr_x_xxx_xyyz, tr_x_xxx_xyyzz, tr_x_xxx_xyzz, tr_x_xxx_xyzzz, tr_x_xxx_xzzz, tr_x_xxx_xzzzz, tr_x_xxx_yyyy, tr_x_xxx_yyyyy, tr_x_xxx_yyyyz, tr_x_xxx_yyyz, tr_x_xxx_yyyzz, tr_x_xxx_yyzz, tr_x_xxx_yyzzz, tr_x_xxx_yzzz, tr_x_xxx_yzzzz, tr_x_xxx_zzzz, tr_x_xxx_zzzzz, tr_x_xxxy_xxxxx, tr_x_xxxy_xxxxy, tr_x_xxxy_xxxxz, tr_x_xxxy_xxxyy, tr_x_xxxy_xxxyz, tr_x_xxxy_xxxzz, tr_x_xxxy_xxyyy, tr_x_xxxy_xxyyz, tr_x_xxxy_xxyzz, tr_x_xxxy_xxzzz, tr_x_xxxy_xyyyy, tr_x_xxxy_xyyyz, tr_x_xxxy_xyyzz, tr_x_xxxy_xyzzz, tr_x_xxxy_xzzzz, tr_x_xxxy_yyyyy, tr_x_xxxy_yyyyz, tr_x_xxxy_yyyzz, tr_x_xxxy_yyzzz, tr_x_xxxy_yzzzz, tr_x_xxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxy_xxxxx[i] = tr_x_xxx_xxxxx[i] * pa_y[i];

        tr_x_xxxy_xxxxy[i] = tr_x_xxx_xxxx[i] * fe_0 + tr_x_xxx_xxxxy[i] * pa_y[i];

        tr_x_xxxy_xxxxz[i] = tr_x_xxx_xxxxz[i] * pa_y[i];

        tr_x_xxxy_xxxyy[i] = 2.0 * tr_x_xxx_xxxy[i] * fe_0 + tr_x_xxx_xxxyy[i] * pa_y[i];

        tr_x_xxxy_xxxyz[i] = tr_x_xxx_xxxz[i] * fe_0 + tr_x_xxx_xxxyz[i] * pa_y[i];

        tr_x_xxxy_xxxzz[i] = tr_x_xxx_xxxzz[i] * pa_y[i];

        tr_x_xxxy_xxyyy[i] = 3.0 * tr_x_xxx_xxyy[i] * fe_0 + tr_x_xxx_xxyyy[i] * pa_y[i];

        tr_x_xxxy_xxyyz[i] = 2.0 * tr_x_xxx_xxyz[i] * fe_0 + tr_x_xxx_xxyyz[i] * pa_y[i];

        tr_x_xxxy_xxyzz[i] = tr_x_xxx_xxzz[i] * fe_0 + tr_x_xxx_xxyzz[i] * pa_y[i];

        tr_x_xxxy_xxzzz[i] = tr_x_xxx_xxzzz[i] * pa_y[i];

        tr_x_xxxy_xyyyy[i] = 4.0 * tr_x_xxx_xyyy[i] * fe_0 + tr_x_xxx_xyyyy[i] * pa_y[i];

        tr_x_xxxy_xyyyz[i] = 3.0 * tr_x_xxx_xyyz[i] * fe_0 + tr_x_xxx_xyyyz[i] * pa_y[i];

        tr_x_xxxy_xyyzz[i] = 2.0 * tr_x_xxx_xyzz[i] * fe_0 + tr_x_xxx_xyyzz[i] * pa_y[i];

        tr_x_xxxy_xyzzz[i] = tr_x_xxx_xzzz[i] * fe_0 + tr_x_xxx_xyzzz[i] * pa_y[i];

        tr_x_xxxy_xzzzz[i] = tr_x_xxx_xzzzz[i] * pa_y[i];

        tr_x_xxxy_yyyyy[i] = 5.0 * tr_x_xxx_yyyy[i] * fe_0 + tr_x_xxx_yyyyy[i] * pa_y[i];

        tr_x_xxxy_yyyyz[i] = 4.0 * tr_x_xxx_yyyz[i] * fe_0 + tr_x_xxx_yyyyz[i] * pa_y[i];

        tr_x_xxxy_yyyzz[i] = 3.0 * tr_x_xxx_yyzz[i] * fe_0 + tr_x_xxx_yyyzz[i] * pa_y[i];

        tr_x_xxxy_yyzzz[i] = 2.0 * tr_x_xxx_yzzz[i] * fe_0 + tr_x_xxx_yyzzz[i] * pa_y[i];

        tr_x_xxxy_yzzzz[i] = tr_x_xxx_zzzz[i] * fe_0 + tr_x_xxx_yzzzz[i] * pa_y[i];

        tr_x_xxxy_zzzzz[i] = tr_x_xxx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : GH

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

    auto tr_x_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 58);

    auto tr_x_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 59);

    auto tr_x_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 60);

    auto tr_x_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 61);

    auto tr_x_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 62);

    #pragma omp simd aligned(pa_z, tr_x_xxx_xxxx, tr_x_xxx_xxxxx, tr_x_xxx_xxxxy, tr_x_xxx_xxxxz, tr_x_xxx_xxxy, tr_x_xxx_xxxyy, tr_x_xxx_xxxyz, tr_x_xxx_xxxz, tr_x_xxx_xxxzz, tr_x_xxx_xxyy, tr_x_xxx_xxyyy, tr_x_xxx_xxyyz, tr_x_xxx_xxyz, tr_x_xxx_xxyzz, tr_x_xxx_xxzz, tr_x_xxx_xxzzz, tr_x_xxx_xyyy, tr_x_xxx_xyyyy, tr_x_xxx_xyyyz, tr_x_xxx_xyyz, tr_x_xxx_xyyzz, tr_x_xxx_xyzz, tr_x_xxx_xyzzz, tr_x_xxx_xzzz, tr_x_xxx_xzzzz, tr_x_xxx_yyyy, tr_x_xxx_yyyyy, tr_x_xxx_yyyyz, tr_x_xxx_yyyz, tr_x_xxx_yyyzz, tr_x_xxx_yyzz, tr_x_xxx_yyzzz, tr_x_xxx_yzzz, tr_x_xxx_yzzzz, tr_x_xxx_zzzz, tr_x_xxx_zzzzz, tr_x_xxxz_xxxxx, tr_x_xxxz_xxxxy, tr_x_xxxz_xxxxz, tr_x_xxxz_xxxyy, tr_x_xxxz_xxxyz, tr_x_xxxz_xxxzz, tr_x_xxxz_xxyyy, tr_x_xxxz_xxyyz, tr_x_xxxz_xxyzz, tr_x_xxxz_xxzzz, tr_x_xxxz_xyyyy, tr_x_xxxz_xyyyz, tr_x_xxxz_xyyzz, tr_x_xxxz_xyzzz, tr_x_xxxz_xzzzz, tr_x_xxxz_yyyyy, tr_x_xxxz_yyyyz, tr_x_xxxz_yyyzz, tr_x_xxxz_yyzzz, tr_x_xxxz_yzzzz, tr_x_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxz_xxxxx[i] = tr_x_xxx_xxxxx[i] * pa_z[i];

        tr_x_xxxz_xxxxy[i] = tr_x_xxx_xxxxy[i] * pa_z[i];

        tr_x_xxxz_xxxxz[i] = tr_x_xxx_xxxx[i] * fe_0 + tr_x_xxx_xxxxz[i] * pa_z[i];

        tr_x_xxxz_xxxyy[i] = tr_x_xxx_xxxyy[i] * pa_z[i];

        tr_x_xxxz_xxxyz[i] = tr_x_xxx_xxxy[i] * fe_0 + tr_x_xxx_xxxyz[i] * pa_z[i];

        tr_x_xxxz_xxxzz[i] = 2.0 * tr_x_xxx_xxxz[i] * fe_0 + tr_x_xxx_xxxzz[i] * pa_z[i];

        tr_x_xxxz_xxyyy[i] = tr_x_xxx_xxyyy[i] * pa_z[i];

        tr_x_xxxz_xxyyz[i] = tr_x_xxx_xxyy[i] * fe_0 + tr_x_xxx_xxyyz[i] * pa_z[i];

        tr_x_xxxz_xxyzz[i] = 2.0 * tr_x_xxx_xxyz[i] * fe_0 + tr_x_xxx_xxyzz[i] * pa_z[i];

        tr_x_xxxz_xxzzz[i] = 3.0 * tr_x_xxx_xxzz[i] * fe_0 + tr_x_xxx_xxzzz[i] * pa_z[i];

        tr_x_xxxz_xyyyy[i] = tr_x_xxx_xyyyy[i] * pa_z[i];

        tr_x_xxxz_xyyyz[i] = tr_x_xxx_xyyy[i] * fe_0 + tr_x_xxx_xyyyz[i] * pa_z[i];

        tr_x_xxxz_xyyzz[i] = 2.0 * tr_x_xxx_xyyz[i] * fe_0 + tr_x_xxx_xyyzz[i] * pa_z[i];

        tr_x_xxxz_xyzzz[i] = 3.0 * tr_x_xxx_xyzz[i] * fe_0 + tr_x_xxx_xyzzz[i] * pa_z[i];

        tr_x_xxxz_xzzzz[i] = 4.0 * tr_x_xxx_xzzz[i] * fe_0 + tr_x_xxx_xzzzz[i] * pa_z[i];

        tr_x_xxxz_yyyyy[i] = tr_x_xxx_yyyyy[i] * pa_z[i];

        tr_x_xxxz_yyyyz[i] = tr_x_xxx_yyyy[i] * fe_0 + tr_x_xxx_yyyyz[i] * pa_z[i];

        tr_x_xxxz_yyyzz[i] = 2.0 * tr_x_xxx_yyyz[i] * fe_0 + tr_x_xxx_yyyzz[i] * pa_z[i];

        tr_x_xxxz_yyzzz[i] = 3.0 * tr_x_xxx_yyzz[i] * fe_0 + tr_x_xxx_yyzzz[i] * pa_z[i];

        tr_x_xxxz_yzzzz[i] = 4.0 * tr_x_xxx_yzzz[i] * fe_0 + tr_x_xxx_yzzzz[i] * pa_z[i];

        tr_x_xxxz_zzzzz[i] = 5.0 * tr_x_xxx_zzzz[i] * fe_0 + tr_x_xxx_zzzzz[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xx_xxxxx, tr_x_xx_xxxxy, tr_x_xx_xxxxz, tr_x_xx_xxxyy, tr_x_xx_xxxyz, tr_x_xx_xxxzz, tr_x_xx_xxyyy, tr_x_xx_xxyyz, tr_x_xx_xxyzz, tr_x_xx_xxzzz, tr_x_xx_xyyyy, tr_x_xx_xyyyz, tr_x_xx_xyyzz, tr_x_xx_xyzzz, tr_x_xx_xzzzz, tr_x_xx_zzzzz, tr_x_xxy_xxxx, tr_x_xxy_xxxxx, tr_x_xxy_xxxxy, tr_x_xxy_xxxxz, tr_x_xxy_xxxy, tr_x_xxy_xxxyy, tr_x_xxy_xxxyz, tr_x_xxy_xxxz, tr_x_xxy_xxxzz, tr_x_xxy_xxyy, tr_x_xxy_xxyyy, tr_x_xxy_xxyyz, tr_x_xxy_xxyz, tr_x_xxy_xxyzz, tr_x_xxy_xxzz, tr_x_xxy_xxzzz, tr_x_xxy_xyyy, tr_x_xxy_xyyyy, tr_x_xxy_xyyyz, tr_x_xxy_xyyz, tr_x_xxy_xyyzz, tr_x_xxy_xyzz, tr_x_xxy_xyzzz, tr_x_xxy_xzzz, tr_x_xxy_xzzzz, tr_x_xxy_zzzzz, tr_x_xxyy_xxxxx, tr_x_xxyy_xxxxy, tr_x_xxyy_xxxxz, tr_x_xxyy_xxxyy, tr_x_xxyy_xxxyz, tr_x_xxyy_xxxzz, tr_x_xxyy_xxyyy, tr_x_xxyy_xxyyz, tr_x_xxyy_xxyzz, tr_x_xxyy_xxzzz, tr_x_xxyy_xyyyy, tr_x_xxyy_xyyyz, tr_x_xxyy_xyyzz, tr_x_xxyy_xyzzz, tr_x_xxyy_xzzzz, tr_x_xxyy_yyyyy, tr_x_xxyy_yyyyz, tr_x_xxyy_yyyzz, tr_x_xxyy_yyzzz, tr_x_xxyy_yzzzz, tr_x_xxyy_zzzzz, tr_x_xyy_yyyyy, tr_x_xyy_yyyyz, tr_x_xyy_yyyzz, tr_x_xyy_yyzzz, tr_x_xyy_yzzzz, tr_x_yy_yyyyy, tr_x_yy_yyyyz, tr_x_yy_yyyzz, tr_x_yy_yyzzz, tr_x_yy_yzzzz, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyzz, ts_xyy_yyzzz, ts_xyy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyy_xxxxx[i] = tr_x_xx_xxxxx[i] * fe_0 + tr_x_xxy_xxxxx[i] * pa_y[i];

        tr_x_xxyy_xxxxy[i] = tr_x_xx_xxxxy[i] * fe_0 + tr_x_xxy_xxxx[i] * fe_0 + tr_x_xxy_xxxxy[i] * pa_y[i];

        tr_x_xxyy_xxxxz[i] = tr_x_xx_xxxxz[i] * fe_0 + tr_x_xxy_xxxxz[i] * pa_y[i];

        tr_x_xxyy_xxxyy[i] = tr_x_xx_xxxyy[i] * fe_0 + 2.0 * tr_x_xxy_xxxy[i] * fe_0 + tr_x_xxy_xxxyy[i] * pa_y[i];

        tr_x_xxyy_xxxyz[i] = tr_x_xx_xxxyz[i] * fe_0 + tr_x_xxy_xxxz[i] * fe_0 + tr_x_xxy_xxxyz[i] * pa_y[i];

        tr_x_xxyy_xxxzz[i] = tr_x_xx_xxxzz[i] * fe_0 + tr_x_xxy_xxxzz[i] * pa_y[i];

        tr_x_xxyy_xxyyy[i] = tr_x_xx_xxyyy[i] * fe_0 + 3.0 * tr_x_xxy_xxyy[i] * fe_0 + tr_x_xxy_xxyyy[i] * pa_y[i];

        tr_x_xxyy_xxyyz[i] = tr_x_xx_xxyyz[i] * fe_0 + 2.0 * tr_x_xxy_xxyz[i] * fe_0 + tr_x_xxy_xxyyz[i] * pa_y[i];

        tr_x_xxyy_xxyzz[i] = tr_x_xx_xxyzz[i] * fe_0 + tr_x_xxy_xxzz[i] * fe_0 + tr_x_xxy_xxyzz[i] * pa_y[i];

        tr_x_xxyy_xxzzz[i] = tr_x_xx_xxzzz[i] * fe_0 + tr_x_xxy_xxzzz[i] * pa_y[i];

        tr_x_xxyy_xyyyy[i] = tr_x_xx_xyyyy[i] * fe_0 + 4.0 * tr_x_xxy_xyyy[i] * fe_0 + tr_x_xxy_xyyyy[i] * pa_y[i];

        tr_x_xxyy_xyyyz[i] = tr_x_xx_xyyyz[i] * fe_0 + 3.0 * tr_x_xxy_xyyz[i] * fe_0 + tr_x_xxy_xyyyz[i] * pa_y[i];

        tr_x_xxyy_xyyzz[i] = tr_x_xx_xyyzz[i] * fe_0 + 2.0 * tr_x_xxy_xyzz[i] * fe_0 + tr_x_xxy_xyyzz[i] * pa_y[i];

        tr_x_xxyy_xyzzz[i] = tr_x_xx_xyzzz[i] * fe_0 + tr_x_xxy_xzzz[i] * fe_0 + tr_x_xxy_xyzzz[i] * pa_y[i];

        tr_x_xxyy_xzzzz[i] = tr_x_xx_xzzzz[i] * fe_0 + tr_x_xxy_xzzzz[i] * pa_y[i];

        tr_x_xxyy_yyyyy[i] = tr_x_yy_yyyyy[i] * fe_0 + ts_xyy_yyyyy[i] * fe_0 + tr_x_xyy_yyyyy[i] * pa_x[i];

        tr_x_xxyy_yyyyz[i] = tr_x_yy_yyyyz[i] * fe_0 + ts_xyy_yyyyz[i] * fe_0 + tr_x_xyy_yyyyz[i] * pa_x[i];

        tr_x_xxyy_yyyzz[i] = tr_x_yy_yyyzz[i] * fe_0 + ts_xyy_yyyzz[i] * fe_0 + tr_x_xyy_yyyzz[i] * pa_x[i];

        tr_x_xxyy_yyzzz[i] = tr_x_yy_yyzzz[i] * fe_0 + ts_xyy_yyzzz[i] * fe_0 + tr_x_xyy_yyzzz[i] * pa_x[i];

        tr_x_xxyy_yzzzz[i] = tr_x_yy_yzzzz[i] * fe_0 + ts_xyy_yzzzz[i] * fe_0 + tr_x_xyy_yzzzz[i] * pa_x[i];

        tr_x_xxyy_zzzzz[i] = tr_x_xx_zzzzz[i] * fe_0 + tr_x_xxy_zzzzz[i] * pa_y[i];
    }

    // Set up 84-105 components of targeted buffer : GH

    auto tr_x_xxyz_xxxxx = pbuffer.data(idx_dip_gh + 84);

    auto tr_x_xxyz_xxxxy = pbuffer.data(idx_dip_gh + 85);

    auto tr_x_xxyz_xxxxz = pbuffer.data(idx_dip_gh + 86);

    auto tr_x_xxyz_xxxyy = pbuffer.data(idx_dip_gh + 87);

    auto tr_x_xxyz_xxxyz = pbuffer.data(idx_dip_gh + 88);

    auto tr_x_xxyz_xxxzz = pbuffer.data(idx_dip_gh + 89);

    auto tr_x_xxyz_xxyyy = pbuffer.data(idx_dip_gh + 90);

    auto tr_x_xxyz_xxyyz = pbuffer.data(idx_dip_gh + 91);

    auto tr_x_xxyz_xxyzz = pbuffer.data(idx_dip_gh + 92);

    auto tr_x_xxyz_xxzzz = pbuffer.data(idx_dip_gh + 93);

    auto tr_x_xxyz_xyyyy = pbuffer.data(idx_dip_gh + 94);

    auto tr_x_xxyz_xyyyz = pbuffer.data(idx_dip_gh + 95);

    auto tr_x_xxyz_xyyzz = pbuffer.data(idx_dip_gh + 96);

    auto tr_x_xxyz_xyzzz = pbuffer.data(idx_dip_gh + 97);

    auto tr_x_xxyz_xzzzz = pbuffer.data(idx_dip_gh + 98);

    auto tr_x_xxyz_yyyyy = pbuffer.data(idx_dip_gh + 99);

    auto tr_x_xxyz_yyyyz = pbuffer.data(idx_dip_gh + 100);

    auto tr_x_xxyz_yyyzz = pbuffer.data(idx_dip_gh + 101);

    auto tr_x_xxyz_yyzzz = pbuffer.data(idx_dip_gh + 102);

    auto tr_x_xxyz_yzzzz = pbuffer.data(idx_dip_gh + 103);

    auto tr_x_xxyz_zzzzz = pbuffer.data(idx_dip_gh + 104);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxy_xxxxy, tr_x_xxy_xxxyy, tr_x_xxy_xxyyy, tr_x_xxy_xyyyy, tr_x_xxy_yyyyy, tr_x_xxyz_xxxxx, tr_x_xxyz_xxxxy, tr_x_xxyz_xxxxz, tr_x_xxyz_xxxyy, tr_x_xxyz_xxxyz, tr_x_xxyz_xxxzz, tr_x_xxyz_xxyyy, tr_x_xxyz_xxyyz, tr_x_xxyz_xxyzz, tr_x_xxyz_xxzzz, tr_x_xxyz_xyyyy, tr_x_xxyz_xyyyz, tr_x_xxyz_xyyzz, tr_x_xxyz_xyzzz, tr_x_xxyz_xzzzz, tr_x_xxyz_yyyyy, tr_x_xxyz_yyyyz, tr_x_xxyz_yyyzz, tr_x_xxyz_yyzzz, tr_x_xxyz_yzzzz, tr_x_xxyz_zzzzz, tr_x_xxz_xxxxx, tr_x_xxz_xxxxz, tr_x_xxz_xxxyz, tr_x_xxz_xxxz, tr_x_xxz_xxxzz, tr_x_xxz_xxyyz, tr_x_xxz_xxyz, tr_x_xxz_xxyzz, tr_x_xxz_xxzz, tr_x_xxz_xxzzz, tr_x_xxz_xyyyz, tr_x_xxz_xyyz, tr_x_xxz_xyyzz, tr_x_xxz_xyzz, tr_x_xxz_xyzzz, tr_x_xxz_xzzz, tr_x_xxz_xzzzz, tr_x_xxz_yyyyz, tr_x_xxz_yyyz, tr_x_xxz_yyyzz, tr_x_xxz_yyzz, tr_x_xxz_yyzzz, tr_x_xxz_yzzz, tr_x_xxz_yzzzz, tr_x_xxz_zzzz, tr_x_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyz_xxxxx[i] = tr_x_xxz_xxxxx[i] * pa_y[i];

        tr_x_xxyz_xxxxy[i] = tr_x_xxy_xxxxy[i] * pa_z[i];

        tr_x_xxyz_xxxxz[i] = tr_x_xxz_xxxxz[i] * pa_y[i];

        tr_x_xxyz_xxxyy[i] = tr_x_xxy_xxxyy[i] * pa_z[i];

        tr_x_xxyz_xxxyz[i] = tr_x_xxz_xxxz[i] * fe_0 + tr_x_xxz_xxxyz[i] * pa_y[i];

        tr_x_xxyz_xxxzz[i] = tr_x_xxz_xxxzz[i] * pa_y[i];

        tr_x_xxyz_xxyyy[i] = tr_x_xxy_xxyyy[i] * pa_z[i];

        tr_x_xxyz_xxyyz[i] = 2.0 * tr_x_xxz_xxyz[i] * fe_0 + tr_x_xxz_xxyyz[i] * pa_y[i];

        tr_x_xxyz_xxyzz[i] = tr_x_xxz_xxzz[i] * fe_0 + tr_x_xxz_xxyzz[i] * pa_y[i];

        tr_x_xxyz_xxzzz[i] = tr_x_xxz_xxzzz[i] * pa_y[i];

        tr_x_xxyz_xyyyy[i] = tr_x_xxy_xyyyy[i] * pa_z[i];

        tr_x_xxyz_xyyyz[i] = 3.0 * tr_x_xxz_xyyz[i] * fe_0 + tr_x_xxz_xyyyz[i] * pa_y[i];

        tr_x_xxyz_xyyzz[i] = 2.0 * tr_x_xxz_xyzz[i] * fe_0 + tr_x_xxz_xyyzz[i] * pa_y[i];

        tr_x_xxyz_xyzzz[i] = tr_x_xxz_xzzz[i] * fe_0 + tr_x_xxz_xyzzz[i] * pa_y[i];

        tr_x_xxyz_xzzzz[i] = tr_x_xxz_xzzzz[i] * pa_y[i];

        tr_x_xxyz_yyyyy[i] = tr_x_xxy_yyyyy[i] * pa_z[i];

        tr_x_xxyz_yyyyz[i] = 4.0 * tr_x_xxz_yyyz[i] * fe_0 + tr_x_xxz_yyyyz[i] * pa_y[i];

        tr_x_xxyz_yyyzz[i] = 3.0 * tr_x_xxz_yyzz[i] * fe_0 + tr_x_xxz_yyyzz[i] * pa_y[i];

        tr_x_xxyz_yyzzz[i] = 2.0 * tr_x_xxz_yzzz[i] * fe_0 + tr_x_xxz_yyzzz[i] * pa_y[i];

        tr_x_xxyz_yzzzz[i] = tr_x_xxz_zzzz[i] * fe_0 + tr_x_xxz_yzzzz[i] * pa_y[i];

        tr_x_xxyz_zzzzz[i] = tr_x_xxz_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xx_xxxxx, tr_x_xx_xxxxy, tr_x_xx_xxxxz, tr_x_xx_xxxyy, tr_x_xx_xxxyz, tr_x_xx_xxxzz, tr_x_xx_xxyyy, tr_x_xx_xxyyz, tr_x_xx_xxyzz, tr_x_xx_xxzzz, tr_x_xx_xyyyy, tr_x_xx_xyyyz, tr_x_xx_xyyzz, tr_x_xx_xyzzz, tr_x_xx_xzzzz, tr_x_xx_yyyyy, tr_x_xxz_xxxx, tr_x_xxz_xxxxx, tr_x_xxz_xxxxy, tr_x_xxz_xxxxz, tr_x_xxz_xxxy, tr_x_xxz_xxxyy, tr_x_xxz_xxxyz, tr_x_xxz_xxxz, tr_x_xxz_xxxzz, tr_x_xxz_xxyy, tr_x_xxz_xxyyy, tr_x_xxz_xxyyz, tr_x_xxz_xxyz, tr_x_xxz_xxyzz, tr_x_xxz_xxzz, tr_x_xxz_xxzzz, tr_x_xxz_xyyy, tr_x_xxz_xyyyy, tr_x_xxz_xyyyz, tr_x_xxz_xyyz, tr_x_xxz_xyyzz, tr_x_xxz_xyzz, tr_x_xxz_xyzzz, tr_x_xxz_xzzz, tr_x_xxz_xzzzz, tr_x_xxz_yyyyy, tr_x_xxzz_xxxxx, tr_x_xxzz_xxxxy, tr_x_xxzz_xxxxz, tr_x_xxzz_xxxyy, tr_x_xxzz_xxxyz, tr_x_xxzz_xxxzz, tr_x_xxzz_xxyyy, tr_x_xxzz_xxyyz, tr_x_xxzz_xxyzz, tr_x_xxzz_xxzzz, tr_x_xxzz_xyyyy, tr_x_xxzz_xyyyz, tr_x_xxzz_xyyzz, tr_x_xxzz_xyzzz, tr_x_xxzz_xzzzz, tr_x_xxzz_yyyyy, tr_x_xxzz_yyyyz, tr_x_xxzz_yyyzz, tr_x_xxzz_yyzzz, tr_x_xxzz_yzzzz, tr_x_xxzz_zzzzz, tr_x_xzz_yyyyz, tr_x_xzz_yyyzz, tr_x_xzz_yyzzz, tr_x_xzz_yzzzz, tr_x_xzz_zzzzz, tr_x_zz_yyyyz, tr_x_zz_yyyzz, tr_x_zz_yyzzz, tr_x_zz_yzzzz, tr_x_zz_zzzzz, ts_xzz_yyyyz, ts_xzz_yyyzz, ts_xzz_yyzzz, ts_xzz_yzzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzz_xxxxx[i] = tr_x_xx_xxxxx[i] * fe_0 + tr_x_xxz_xxxxx[i] * pa_z[i];

        tr_x_xxzz_xxxxy[i] = tr_x_xx_xxxxy[i] * fe_0 + tr_x_xxz_xxxxy[i] * pa_z[i];

        tr_x_xxzz_xxxxz[i] = tr_x_xx_xxxxz[i] * fe_0 + tr_x_xxz_xxxx[i] * fe_0 + tr_x_xxz_xxxxz[i] * pa_z[i];

        tr_x_xxzz_xxxyy[i] = tr_x_xx_xxxyy[i] * fe_0 + tr_x_xxz_xxxyy[i] * pa_z[i];

        tr_x_xxzz_xxxyz[i] = tr_x_xx_xxxyz[i] * fe_0 + tr_x_xxz_xxxy[i] * fe_0 + tr_x_xxz_xxxyz[i] * pa_z[i];

        tr_x_xxzz_xxxzz[i] = tr_x_xx_xxxzz[i] * fe_0 + 2.0 * tr_x_xxz_xxxz[i] * fe_0 + tr_x_xxz_xxxzz[i] * pa_z[i];

        tr_x_xxzz_xxyyy[i] = tr_x_xx_xxyyy[i] * fe_0 + tr_x_xxz_xxyyy[i] * pa_z[i];

        tr_x_xxzz_xxyyz[i] = tr_x_xx_xxyyz[i] * fe_0 + tr_x_xxz_xxyy[i] * fe_0 + tr_x_xxz_xxyyz[i] * pa_z[i];

        tr_x_xxzz_xxyzz[i] = tr_x_xx_xxyzz[i] * fe_0 + 2.0 * tr_x_xxz_xxyz[i] * fe_0 + tr_x_xxz_xxyzz[i] * pa_z[i];

        tr_x_xxzz_xxzzz[i] = tr_x_xx_xxzzz[i] * fe_0 + 3.0 * tr_x_xxz_xxzz[i] * fe_0 + tr_x_xxz_xxzzz[i] * pa_z[i];

        tr_x_xxzz_xyyyy[i] = tr_x_xx_xyyyy[i] * fe_0 + tr_x_xxz_xyyyy[i] * pa_z[i];

        tr_x_xxzz_xyyyz[i] = tr_x_xx_xyyyz[i] * fe_0 + tr_x_xxz_xyyy[i] * fe_0 + tr_x_xxz_xyyyz[i] * pa_z[i];

        tr_x_xxzz_xyyzz[i] = tr_x_xx_xyyzz[i] * fe_0 + 2.0 * tr_x_xxz_xyyz[i] * fe_0 + tr_x_xxz_xyyzz[i] * pa_z[i];

        tr_x_xxzz_xyzzz[i] = tr_x_xx_xyzzz[i] * fe_0 + 3.0 * tr_x_xxz_xyzz[i] * fe_0 + tr_x_xxz_xyzzz[i] * pa_z[i];

        tr_x_xxzz_xzzzz[i] = tr_x_xx_xzzzz[i] * fe_0 + 4.0 * tr_x_xxz_xzzz[i] * fe_0 + tr_x_xxz_xzzzz[i] * pa_z[i];

        tr_x_xxzz_yyyyy[i] = tr_x_xx_yyyyy[i] * fe_0 + tr_x_xxz_yyyyy[i] * pa_z[i];

        tr_x_xxzz_yyyyz[i] = tr_x_zz_yyyyz[i] * fe_0 + ts_xzz_yyyyz[i] * fe_0 + tr_x_xzz_yyyyz[i] * pa_x[i];

        tr_x_xxzz_yyyzz[i] = tr_x_zz_yyyzz[i] * fe_0 + ts_xzz_yyyzz[i] * fe_0 + tr_x_xzz_yyyzz[i] * pa_x[i];

        tr_x_xxzz_yyzzz[i] = tr_x_zz_yyzzz[i] * fe_0 + ts_xzz_yyzzz[i] * fe_0 + tr_x_xzz_yyzzz[i] * pa_x[i];

        tr_x_xxzz_yzzzz[i] = tr_x_zz_yzzzz[i] * fe_0 + ts_xzz_yzzzz[i] * fe_0 + tr_x_xzz_yzzzz[i] * pa_x[i];

        tr_x_xxzz_zzzzz[i] = tr_x_zz_zzzzz[i] * fe_0 + ts_xzz_zzzzz[i] * fe_0 + tr_x_xzz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : GH

    auto tr_x_xyyy_xxxxx = pbuffer.data(idx_dip_gh + 126);

    auto tr_x_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 127);

    auto tr_x_xyyy_xxxxz = pbuffer.data(idx_dip_gh + 128);

    auto tr_x_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 129);

    auto tr_x_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 130);

    auto tr_x_xyyy_xxxzz = pbuffer.data(idx_dip_gh + 131);

    auto tr_x_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 132);

    auto tr_x_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 133);

    auto tr_x_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 134);

    auto tr_x_xyyy_xxzzz = pbuffer.data(idx_dip_gh + 135);

    auto tr_x_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 136);

    auto tr_x_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 137);

    auto tr_x_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 138);

    auto tr_x_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 139);

    auto tr_x_xyyy_xzzzz = pbuffer.data(idx_dip_gh + 140);

    auto tr_x_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 141);

    auto tr_x_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 142);

    auto tr_x_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 143);

    auto tr_x_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 144);

    auto tr_x_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 145);

    auto tr_x_xyyy_zzzzz = pbuffer.data(idx_dip_gh + 146);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xy_xxxxx, tr_x_xy_xxxxz, tr_x_xy_xxxzz, tr_x_xy_xxzzz, tr_x_xy_xzzzz, tr_x_xyy_xxxxx, tr_x_xyy_xxxxz, tr_x_xyy_xxxzz, tr_x_xyy_xxzzz, tr_x_xyy_xzzzz, tr_x_xyyy_xxxxx, tr_x_xyyy_xxxxy, tr_x_xyyy_xxxxz, tr_x_xyyy_xxxyy, tr_x_xyyy_xxxyz, tr_x_xyyy_xxxzz, tr_x_xyyy_xxyyy, tr_x_xyyy_xxyyz, tr_x_xyyy_xxyzz, tr_x_xyyy_xxzzz, tr_x_xyyy_xyyyy, tr_x_xyyy_xyyyz, tr_x_xyyy_xyyzz, tr_x_xyyy_xyzzz, tr_x_xyyy_xzzzz, tr_x_xyyy_yyyyy, tr_x_xyyy_yyyyz, tr_x_xyyy_yyyzz, tr_x_xyyy_yyzzz, tr_x_xyyy_yzzzz, tr_x_xyyy_zzzzz, tr_x_yyy_xxxxy, tr_x_yyy_xxxy, tr_x_yyy_xxxyy, tr_x_yyy_xxxyz, tr_x_yyy_xxyy, tr_x_yyy_xxyyy, tr_x_yyy_xxyyz, tr_x_yyy_xxyz, tr_x_yyy_xxyzz, tr_x_yyy_xyyy, tr_x_yyy_xyyyy, tr_x_yyy_xyyyz, tr_x_yyy_xyyz, tr_x_yyy_xyyzz, tr_x_yyy_xyzz, tr_x_yyy_xyzzz, tr_x_yyy_yyyy, tr_x_yyy_yyyyy, tr_x_yyy_yyyyz, tr_x_yyy_yyyz, tr_x_yyy_yyyzz, tr_x_yyy_yyzz, tr_x_yyy_yyzzz, tr_x_yyy_yzzz, tr_x_yyy_yzzzz, tr_x_yyy_zzzzz, ts_yyy_xxxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyy_xxxxx[i] = 2.0 * tr_x_xy_xxxxx[i] * fe_0 + tr_x_xyy_xxxxx[i] * pa_y[i];

        tr_x_xyyy_xxxxy[i] = 4.0 * tr_x_yyy_xxxy[i] * fe_0 + ts_yyy_xxxxy[i] * fe_0 + tr_x_yyy_xxxxy[i] * pa_x[i];

        tr_x_xyyy_xxxxz[i] = 2.0 * tr_x_xy_xxxxz[i] * fe_0 + tr_x_xyy_xxxxz[i] * pa_y[i];

        tr_x_xyyy_xxxyy[i] = 3.0 * tr_x_yyy_xxyy[i] * fe_0 + ts_yyy_xxxyy[i] * fe_0 + tr_x_yyy_xxxyy[i] * pa_x[i];

        tr_x_xyyy_xxxyz[i] = 3.0 * tr_x_yyy_xxyz[i] * fe_0 + ts_yyy_xxxyz[i] * fe_0 + tr_x_yyy_xxxyz[i] * pa_x[i];

        tr_x_xyyy_xxxzz[i] = 2.0 * tr_x_xy_xxxzz[i] * fe_0 + tr_x_xyy_xxxzz[i] * pa_y[i];

        tr_x_xyyy_xxyyy[i] = 2.0 * tr_x_yyy_xyyy[i] * fe_0 + ts_yyy_xxyyy[i] * fe_0 + tr_x_yyy_xxyyy[i] * pa_x[i];

        tr_x_xyyy_xxyyz[i] = 2.0 * tr_x_yyy_xyyz[i] * fe_0 + ts_yyy_xxyyz[i] * fe_0 + tr_x_yyy_xxyyz[i] * pa_x[i];

        tr_x_xyyy_xxyzz[i] = 2.0 * tr_x_yyy_xyzz[i] * fe_0 + ts_yyy_xxyzz[i] * fe_0 + tr_x_yyy_xxyzz[i] * pa_x[i];

        tr_x_xyyy_xxzzz[i] = 2.0 * tr_x_xy_xxzzz[i] * fe_0 + tr_x_xyy_xxzzz[i] * pa_y[i];

        tr_x_xyyy_xyyyy[i] = tr_x_yyy_yyyy[i] * fe_0 + ts_yyy_xyyyy[i] * fe_0 + tr_x_yyy_xyyyy[i] * pa_x[i];

        tr_x_xyyy_xyyyz[i] = tr_x_yyy_yyyz[i] * fe_0 + ts_yyy_xyyyz[i] * fe_0 + tr_x_yyy_xyyyz[i] * pa_x[i];

        tr_x_xyyy_xyyzz[i] = tr_x_yyy_yyzz[i] * fe_0 + ts_yyy_xyyzz[i] * fe_0 + tr_x_yyy_xyyzz[i] * pa_x[i];

        tr_x_xyyy_xyzzz[i] = tr_x_yyy_yzzz[i] * fe_0 + ts_yyy_xyzzz[i] * fe_0 + tr_x_yyy_xyzzz[i] * pa_x[i];

        tr_x_xyyy_xzzzz[i] = 2.0 * tr_x_xy_xzzzz[i] * fe_0 + tr_x_xyy_xzzzz[i] * pa_y[i];

        tr_x_xyyy_yyyyy[i] = ts_yyy_yyyyy[i] * fe_0 + tr_x_yyy_yyyyy[i] * pa_x[i];

        tr_x_xyyy_yyyyz[i] = ts_yyy_yyyyz[i] * fe_0 + tr_x_yyy_yyyyz[i] * pa_x[i];

        tr_x_xyyy_yyyzz[i] = ts_yyy_yyyzz[i] * fe_0 + tr_x_yyy_yyyzz[i] * pa_x[i];

        tr_x_xyyy_yyzzz[i] = ts_yyy_yyzzz[i] * fe_0 + tr_x_yyy_yyzzz[i] * pa_x[i];

        tr_x_xyyy_yzzzz[i] = ts_yyy_yzzzz[i] * fe_0 + tr_x_yyy_yzzzz[i] * pa_x[i];

        tr_x_xyyy_zzzzz[i] = ts_yyy_zzzzz[i] * fe_0 + tr_x_yyy_zzzzz[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : GH

    auto tr_x_xyyz_xxxxx = pbuffer.data(idx_dip_gh + 147);

    auto tr_x_xyyz_xxxxy = pbuffer.data(idx_dip_gh + 148);

    auto tr_x_xyyz_xxxxz = pbuffer.data(idx_dip_gh + 149);

    auto tr_x_xyyz_xxxyy = pbuffer.data(idx_dip_gh + 150);

    auto tr_x_xyyz_xxxyz = pbuffer.data(idx_dip_gh + 151);

    auto tr_x_xyyz_xxxzz = pbuffer.data(idx_dip_gh + 152);

    auto tr_x_xyyz_xxyyy = pbuffer.data(idx_dip_gh + 153);

    auto tr_x_xyyz_xxyyz = pbuffer.data(idx_dip_gh + 154);

    auto tr_x_xyyz_xxyzz = pbuffer.data(idx_dip_gh + 155);

    auto tr_x_xyyz_xxzzz = pbuffer.data(idx_dip_gh + 156);

    auto tr_x_xyyz_xyyyy = pbuffer.data(idx_dip_gh + 157);

    auto tr_x_xyyz_xyyyz = pbuffer.data(idx_dip_gh + 158);

    auto tr_x_xyyz_xyyzz = pbuffer.data(idx_dip_gh + 159);

    auto tr_x_xyyz_xyzzz = pbuffer.data(idx_dip_gh + 160);

    auto tr_x_xyyz_xzzzz = pbuffer.data(idx_dip_gh + 161);

    auto tr_x_xyyz_yyyyy = pbuffer.data(idx_dip_gh + 162);

    auto tr_x_xyyz_yyyyz = pbuffer.data(idx_dip_gh + 163);

    auto tr_x_xyyz_yyyzz = pbuffer.data(idx_dip_gh + 164);

    auto tr_x_xyyz_yyzzz = pbuffer.data(idx_dip_gh + 165);

    auto tr_x_xyyz_yzzzz = pbuffer.data(idx_dip_gh + 166);

    auto tr_x_xyyz_zzzzz = pbuffer.data(idx_dip_gh + 167);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyy_xxxxx, tr_x_xyy_xxxxy, tr_x_xyy_xxxy, tr_x_xyy_xxxyy, tr_x_xyy_xxxyz, tr_x_xyy_xxyy, tr_x_xyy_xxyyy, tr_x_xyy_xxyyz, tr_x_xyy_xxyz, tr_x_xyy_xxyzz, tr_x_xyy_xyyy, tr_x_xyy_xyyyy, tr_x_xyy_xyyyz, tr_x_xyy_xyyz, tr_x_xyy_xyyzz, tr_x_xyy_xyzz, tr_x_xyy_xyzzz, tr_x_xyy_yyyyy, tr_x_xyyz_xxxxx, tr_x_xyyz_xxxxy, tr_x_xyyz_xxxxz, tr_x_xyyz_xxxyy, tr_x_xyyz_xxxyz, tr_x_xyyz_xxxzz, tr_x_xyyz_xxyyy, tr_x_xyyz_xxyyz, tr_x_xyyz_xxyzz, tr_x_xyyz_xxzzz, tr_x_xyyz_xyyyy, tr_x_xyyz_xyyyz, tr_x_xyyz_xyyzz, tr_x_xyyz_xyzzz, tr_x_xyyz_xzzzz, tr_x_xyyz_yyyyy, tr_x_xyyz_yyyyz, tr_x_xyyz_yyyzz, tr_x_xyyz_yyzzz, tr_x_xyyz_yzzzz, tr_x_xyyz_zzzzz, tr_x_xyz_xxxxz, tr_x_xyz_xxxzz, tr_x_xyz_xxzzz, tr_x_xyz_xzzzz, tr_x_xz_xxxxz, tr_x_xz_xxxzz, tr_x_xz_xxzzz, tr_x_xz_xzzzz, tr_x_yyz_yyyyz, tr_x_yyz_yyyzz, tr_x_yyz_yyzzz, tr_x_yyz_yzzzz, tr_x_yyz_zzzzz, ts_yyz_yyyyz, ts_yyz_yyyzz, ts_yyz_yyzzz, ts_yyz_yzzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyz_xxxxx[i] = tr_x_xyy_xxxxx[i] * pa_z[i];

        tr_x_xyyz_xxxxy[i] = tr_x_xyy_xxxxy[i] * pa_z[i];

        tr_x_xyyz_xxxxz[i] = tr_x_xz_xxxxz[i] * fe_0 + tr_x_xyz_xxxxz[i] * pa_y[i];

        tr_x_xyyz_xxxyy[i] = tr_x_xyy_xxxyy[i] * pa_z[i];

        tr_x_xyyz_xxxyz[i] = tr_x_xyy_xxxy[i] * fe_0 + tr_x_xyy_xxxyz[i] * pa_z[i];

        tr_x_xyyz_xxxzz[i] = tr_x_xz_xxxzz[i] * fe_0 + tr_x_xyz_xxxzz[i] * pa_y[i];

        tr_x_xyyz_xxyyy[i] = tr_x_xyy_xxyyy[i] * pa_z[i];

        tr_x_xyyz_xxyyz[i] = tr_x_xyy_xxyy[i] * fe_0 + tr_x_xyy_xxyyz[i] * pa_z[i];

        tr_x_xyyz_xxyzz[i] = 2.0 * tr_x_xyy_xxyz[i] * fe_0 + tr_x_xyy_xxyzz[i] * pa_z[i];

        tr_x_xyyz_xxzzz[i] = tr_x_xz_xxzzz[i] * fe_0 + tr_x_xyz_xxzzz[i] * pa_y[i];

        tr_x_xyyz_xyyyy[i] = tr_x_xyy_xyyyy[i] * pa_z[i];

        tr_x_xyyz_xyyyz[i] = tr_x_xyy_xyyy[i] * fe_0 + tr_x_xyy_xyyyz[i] * pa_z[i];

        tr_x_xyyz_xyyzz[i] = 2.0 * tr_x_xyy_xyyz[i] * fe_0 + tr_x_xyy_xyyzz[i] * pa_z[i];

        tr_x_xyyz_xyzzz[i] = 3.0 * tr_x_xyy_xyzz[i] * fe_0 + tr_x_xyy_xyzzz[i] * pa_z[i];

        tr_x_xyyz_xzzzz[i] = tr_x_xz_xzzzz[i] * fe_0 + tr_x_xyz_xzzzz[i] * pa_y[i];

        tr_x_xyyz_yyyyy[i] = tr_x_xyy_yyyyy[i] * pa_z[i];

        tr_x_xyyz_yyyyz[i] = ts_yyz_yyyyz[i] * fe_0 + tr_x_yyz_yyyyz[i] * pa_x[i];

        tr_x_xyyz_yyyzz[i] = ts_yyz_yyyzz[i] * fe_0 + tr_x_yyz_yyyzz[i] * pa_x[i];

        tr_x_xyyz_yyzzz[i] = ts_yyz_yyzzz[i] * fe_0 + tr_x_yyz_yyzzz[i] * pa_x[i];

        tr_x_xyyz_yzzzz[i] = ts_yyz_yzzzz[i] * fe_0 + tr_x_yyz_yzzzz[i] * pa_x[i];

        tr_x_xyyz_zzzzz[i] = ts_yyz_zzzzz[i] * fe_0 + tr_x_yyz_zzzzz[i] * pa_x[i];
    }

    // Set up 168-189 components of targeted buffer : GH

    auto tr_x_xyzz_xxxxx = pbuffer.data(idx_dip_gh + 168);

    auto tr_x_xyzz_xxxxy = pbuffer.data(idx_dip_gh + 169);

    auto tr_x_xyzz_xxxxz = pbuffer.data(idx_dip_gh + 170);

    auto tr_x_xyzz_xxxyy = pbuffer.data(idx_dip_gh + 171);

    auto tr_x_xyzz_xxxyz = pbuffer.data(idx_dip_gh + 172);

    auto tr_x_xyzz_xxxzz = pbuffer.data(idx_dip_gh + 173);

    auto tr_x_xyzz_xxyyy = pbuffer.data(idx_dip_gh + 174);

    auto tr_x_xyzz_xxyyz = pbuffer.data(idx_dip_gh + 175);

    auto tr_x_xyzz_xxyzz = pbuffer.data(idx_dip_gh + 176);

    auto tr_x_xyzz_xxzzz = pbuffer.data(idx_dip_gh + 177);

    auto tr_x_xyzz_xyyyy = pbuffer.data(idx_dip_gh + 178);

    auto tr_x_xyzz_xyyyz = pbuffer.data(idx_dip_gh + 179);

    auto tr_x_xyzz_xyyzz = pbuffer.data(idx_dip_gh + 180);

    auto tr_x_xyzz_xyzzz = pbuffer.data(idx_dip_gh + 181);

    auto tr_x_xyzz_xzzzz = pbuffer.data(idx_dip_gh + 182);

    auto tr_x_xyzz_yyyyy = pbuffer.data(idx_dip_gh + 183);

    auto tr_x_xyzz_yyyyz = pbuffer.data(idx_dip_gh + 184);

    auto tr_x_xyzz_yyyzz = pbuffer.data(idx_dip_gh + 185);

    auto tr_x_xyzz_yyzzz = pbuffer.data(idx_dip_gh + 186);

    auto tr_x_xyzz_yzzzz = pbuffer.data(idx_dip_gh + 187);

    auto tr_x_xyzz_zzzzz = pbuffer.data(idx_dip_gh + 188);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyzz_xxxxx, tr_x_xyzz_xxxxy, tr_x_xyzz_xxxxz, tr_x_xyzz_xxxyy, tr_x_xyzz_xxxyz, tr_x_xyzz_xxxzz, tr_x_xyzz_xxyyy, tr_x_xyzz_xxyyz, tr_x_xyzz_xxyzz, tr_x_xyzz_xxzzz, tr_x_xyzz_xyyyy, tr_x_xyzz_xyyyz, tr_x_xyzz_xyyzz, tr_x_xyzz_xyzzz, tr_x_xyzz_xzzzz, tr_x_xyzz_yyyyy, tr_x_xyzz_yyyyz, tr_x_xyzz_yyyzz, tr_x_xyzz_yyzzz, tr_x_xyzz_yzzzz, tr_x_xyzz_zzzzz, tr_x_xzz_xxxx, tr_x_xzz_xxxxx, tr_x_xzz_xxxxy, tr_x_xzz_xxxxz, tr_x_xzz_xxxy, tr_x_xzz_xxxyy, tr_x_xzz_xxxyz, tr_x_xzz_xxxz, tr_x_xzz_xxxzz, tr_x_xzz_xxyy, tr_x_xzz_xxyyy, tr_x_xzz_xxyyz, tr_x_xzz_xxyz, tr_x_xzz_xxyzz, tr_x_xzz_xxzz, tr_x_xzz_xxzzz, tr_x_xzz_xyyy, tr_x_xzz_xyyyy, tr_x_xzz_xyyyz, tr_x_xzz_xyyz, tr_x_xzz_xyyzz, tr_x_xzz_xyzz, tr_x_xzz_xyzzz, tr_x_xzz_xzzz, tr_x_xzz_xzzzz, tr_x_xzz_zzzzz, tr_x_yzz_yyyyy, tr_x_yzz_yyyyz, tr_x_yzz_yyyzz, tr_x_yzz_yyzzz, tr_x_yzz_yzzzz, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzz_xxxxx[i] = tr_x_xzz_xxxxx[i] * pa_y[i];

        tr_x_xyzz_xxxxy[i] = tr_x_xzz_xxxx[i] * fe_0 + tr_x_xzz_xxxxy[i] * pa_y[i];

        tr_x_xyzz_xxxxz[i] = tr_x_xzz_xxxxz[i] * pa_y[i];

        tr_x_xyzz_xxxyy[i] = 2.0 * tr_x_xzz_xxxy[i] * fe_0 + tr_x_xzz_xxxyy[i] * pa_y[i];

        tr_x_xyzz_xxxyz[i] = tr_x_xzz_xxxz[i] * fe_0 + tr_x_xzz_xxxyz[i] * pa_y[i];

        tr_x_xyzz_xxxzz[i] = tr_x_xzz_xxxzz[i] * pa_y[i];

        tr_x_xyzz_xxyyy[i] = 3.0 * tr_x_xzz_xxyy[i] * fe_0 + tr_x_xzz_xxyyy[i] * pa_y[i];

        tr_x_xyzz_xxyyz[i] = 2.0 * tr_x_xzz_xxyz[i] * fe_0 + tr_x_xzz_xxyyz[i] * pa_y[i];

        tr_x_xyzz_xxyzz[i] = tr_x_xzz_xxzz[i] * fe_0 + tr_x_xzz_xxyzz[i] * pa_y[i];

        tr_x_xyzz_xxzzz[i] = tr_x_xzz_xxzzz[i] * pa_y[i];

        tr_x_xyzz_xyyyy[i] = 4.0 * tr_x_xzz_xyyy[i] * fe_0 + tr_x_xzz_xyyyy[i] * pa_y[i];

        tr_x_xyzz_xyyyz[i] = 3.0 * tr_x_xzz_xyyz[i] * fe_0 + tr_x_xzz_xyyyz[i] * pa_y[i];

        tr_x_xyzz_xyyzz[i] = 2.0 * tr_x_xzz_xyzz[i] * fe_0 + tr_x_xzz_xyyzz[i] * pa_y[i];

        tr_x_xyzz_xyzzz[i] = tr_x_xzz_xzzz[i] * fe_0 + tr_x_xzz_xyzzz[i] * pa_y[i];

        tr_x_xyzz_xzzzz[i] = tr_x_xzz_xzzzz[i] * pa_y[i];

        tr_x_xyzz_yyyyy[i] = ts_yzz_yyyyy[i] * fe_0 + tr_x_yzz_yyyyy[i] * pa_x[i];

        tr_x_xyzz_yyyyz[i] = ts_yzz_yyyyz[i] * fe_0 + tr_x_yzz_yyyyz[i] * pa_x[i];

        tr_x_xyzz_yyyzz[i] = ts_yzz_yyyzz[i] * fe_0 + tr_x_yzz_yyyzz[i] * pa_x[i];

        tr_x_xyzz_yyzzz[i] = ts_yzz_yyzzz[i] * fe_0 + tr_x_yzz_yyzzz[i] * pa_x[i];

        tr_x_xyzz_yzzzz[i] = ts_yzz_yzzzz[i] * fe_0 + tr_x_yzz_yzzzz[i] * pa_x[i];

        tr_x_xyzz_zzzzz[i] = tr_x_xzz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : GH

    auto tr_x_xzzz_xxxxx = pbuffer.data(idx_dip_gh + 189);

    auto tr_x_xzzz_xxxxy = pbuffer.data(idx_dip_gh + 190);

    auto tr_x_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 191);

    auto tr_x_xzzz_xxxyy = pbuffer.data(idx_dip_gh + 192);

    auto tr_x_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 193);

    auto tr_x_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 194);

    auto tr_x_xzzz_xxyyy = pbuffer.data(idx_dip_gh + 195);

    auto tr_x_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 196);

    auto tr_x_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 197);

    auto tr_x_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 198);

    auto tr_x_xzzz_xyyyy = pbuffer.data(idx_dip_gh + 199);

    auto tr_x_xzzz_xyyyz = pbuffer.data(idx_dip_gh + 200);

    auto tr_x_xzzz_xyyzz = pbuffer.data(idx_dip_gh + 201);

    auto tr_x_xzzz_xyzzz = pbuffer.data(idx_dip_gh + 202);

    auto tr_x_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 203);

    auto tr_x_xzzz_yyyyy = pbuffer.data(idx_dip_gh + 204);

    auto tr_x_xzzz_yyyyz = pbuffer.data(idx_dip_gh + 205);

    auto tr_x_xzzz_yyyzz = pbuffer.data(idx_dip_gh + 206);

    auto tr_x_xzzz_yyzzz = pbuffer.data(idx_dip_gh + 207);

    auto tr_x_xzzz_yzzzz = pbuffer.data(idx_dip_gh + 208);

    auto tr_x_xzzz_zzzzz = pbuffer.data(idx_dip_gh + 209);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xz_xxxxx, tr_x_xz_xxxxy, tr_x_xz_xxxyy, tr_x_xz_xxyyy, tr_x_xz_xyyyy, tr_x_xzz_xxxxx, tr_x_xzz_xxxxy, tr_x_xzz_xxxyy, tr_x_xzz_xxyyy, tr_x_xzz_xyyyy, tr_x_xzzz_xxxxx, tr_x_xzzz_xxxxy, tr_x_xzzz_xxxxz, tr_x_xzzz_xxxyy, tr_x_xzzz_xxxyz, tr_x_xzzz_xxxzz, tr_x_xzzz_xxyyy, tr_x_xzzz_xxyyz, tr_x_xzzz_xxyzz, tr_x_xzzz_xxzzz, tr_x_xzzz_xyyyy, tr_x_xzzz_xyyyz, tr_x_xzzz_xyyzz, tr_x_xzzz_xyzzz, tr_x_xzzz_xzzzz, tr_x_xzzz_yyyyy, tr_x_xzzz_yyyyz, tr_x_xzzz_yyyzz, tr_x_xzzz_yyzzz, tr_x_xzzz_yzzzz, tr_x_xzzz_zzzzz, tr_x_zzz_xxxxz, tr_x_zzz_xxxyz, tr_x_zzz_xxxz, tr_x_zzz_xxxzz, tr_x_zzz_xxyyz, tr_x_zzz_xxyz, tr_x_zzz_xxyzz, tr_x_zzz_xxzz, tr_x_zzz_xxzzz, tr_x_zzz_xyyyz, tr_x_zzz_xyyz, tr_x_zzz_xyyzz, tr_x_zzz_xyzz, tr_x_zzz_xyzzz, tr_x_zzz_xzzz, tr_x_zzz_xzzzz, tr_x_zzz_yyyyy, tr_x_zzz_yyyyz, tr_x_zzz_yyyz, tr_x_zzz_yyyzz, tr_x_zzz_yyzz, tr_x_zzz_yyzzz, tr_x_zzz_yzzz, tr_x_zzz_yzzzz, tr_x_zzz_zzzz, tr_x_zzz_zzzzz, ts_zzz_xxxxz, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzz_xxxxx[i] = 2.0 * tr_x_xz_xxxxx[i] * fe_0 + tr_x_xzz_xxxxx[i] * pa_z[i];

        tr_x_xzzz_xxxxy[i] = 2.0 * tr_x_xz_xxxxy[i] * fe_0 + tr_x_xzz_xxxxy[i] * pa_z[i];

        tr_x_xzzz_xxxxz[i] = 4.0 * tr_x_zzz_xxxz[i] * fe_0 + ts_zzz_xxxxz[i] * fe_0 + tr_x_zzz_xxxxz[i] * pa_x[i];

        tr_x_xzzz_xxxyy[i] = 2.0 * tr_x_xz_xxxyy[i] * fe_0 + tr_x_xzz_xxxyy[i] * pa_z[i];

        tr_x_xzzz_xxxyz[i] = 3.0 * tr_x_zzz_xxyz[i] * fe_0 + ts_zzz_xxxyz[i] * fe_0 + tr_x_zzz_xxxyz[i] * pa_x[i];

        tr_x_xzzz_xxxzz[i] = 3.0 * tr_x_zzz_xxzz[i] * fe_0 + ts_zzz_xxxzz[i] * fe_0 + tr_x_zzz_xxxzz[i] * pa_x[i];

        tr_x_xzzz_xxyyy[i] = 2.0 * tr_x_xz_xxyyy[i] * fe_0 + tr_x_xzz_xxyyy[i] * pa_z[i];

        tr_x_xzzz_xxyyz[i] = 2.0 * tr_x_zzz_xyyz[i] * fe_0 + ts_zzz_xxyyz[i] * fe_0 + tr_x_zzz_xxyyz[i] * pa_x[i];

        tr_x_xzzz_xxyzz[i] = 2.0 * tr_x_zzz_xyzz[i] * fe_0 + ts_zzz_xxyzz[i] * fe_0 + tr_x_zzz_xxyzz[i] * pa_x[i];

        tr_x_xzzz_xxzzz[i] = 2.0 * tr_x_zzz_xzzz[i] * fe_0 + ts_zzz_xxzzz[i] * fe_0 + tr_x_zzz_xxzzz[i] * pa_x[i];

        tr_x_xzzz_xyyyy[i] = 2.0 * tr_x_xz_xyyyy[i] * fe_0 + tr_x_xzz_xyyyy[i] * pa_z[i];

        tr_x_xzzz_xyyyz[i] = tr_x_zzz_yyyz[i] * fe_0 + ts_zzz_xyyyz[i] * fe_0 + tr_x_zzz_xyyyz[i] * pa_x[i];

        tr_x_xzzz_xyyzz[i] = tr_x_zzz_yyzz[i] * fe_0 + ts_zzz_xyyzz[i] * fe_0 + tr_x_zzz_xyyzz[i] * pa_x[i];

        tr_x_xzzz_xyzzz[i] = tr_x_zzz_yzzz[i] * fe_0 + ts_zzz_xyzzz[i] * fe_0 + tr_x_zzz_xyzzz[i] * pa_x[i];

        tr_x_xzzz_xzzzz[i] = tr_x_zzz_zzzz[i] * fe_0 + ts_zzz_xzzzz[i] * fe_0 + tr_x_zzz_xzzzz[i] * pa_x[i];

        tr_x_xzzz_yyyyy[i] = ts_zzz_yyyyy[i] * fe_0 + tr_x_zzz_yyyyy[i] * pa_x[i];

        tr_x_xzzz_yyyyz[i] = ts_zzz_yyyyz[i] * fe_0 + tr_x_zzz_yyyyz[i] * pa_x[i];

        tr_x_xzzz_yyyzz[i] = ts_zzz_yyyzz[i] * fe_0 + tr_x_zzz_yyyzz[i] * pa_x[i];

        tr_x_xzzz_yyzzz[i] = ts_zzz_yyzzz[i] * fe_0 + tr_x_zzz_yyzzz[i] * pa_x[i];

        tr_x_xzzz_yzzzz[i] = ts_zzz_yzzzz[i] * fe_0 + tr_x_zzz_yzzzz[i] * pa_x[i];

        tr_x_xzzz_zzzzz[i] = ts_zzz_zzzzz[i] * fe_0 + tr_x_zzz_zzzzz[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, tr_x_yy_xxxxx, tr_x_yy_xxxxy, tr_x_yy_xxxxz, tr_x_yy_xxxyy, tr_x_yy_xxxyz, tr_x_yy_xxxzz, tr_x_yy_xxyyy, tr_x_yy_xxyyz, tr_x_yy_xxyzz, tr_x_yy_xxzzz, tr_x_yy_xyyyy, tr_x_yy_xyyyz, tr_x_yy_xyyzz, tr_x_yy_xyzzz, tr_x_yy_xzzzz, tr_x_yy_yyyyy, tr_x_yy_yyyyz, tr_x_yy_yyyzz, tr_x_yy_yyzzz, tr_x_yy_yzzzz, tr_x_yy_zzzzz, tr_x_yyy_xxxx, tr_x_yyy_xxxxx, tr_x_yyy_xxxxy, tr_x_yyy_xxxxz, tr_x_yyy_xxxy, tr_x_yyy_xxxyy, tr_x_yyy_xxxyz, tr_x_yyy_xxxz, tr_x_yyy_xxxzz, tr_x_yyy_xxyy, tr_x_yyy_xxyyy, tr_x_yyy_xxyyz, tr_x_yyy_xxyz, tr_x_yyy_xxyzz, tr_x_yyy_xxzz, tr_x_yyy_xxzzz, tr_x_yyy_xyyy, tr_x_yyy_xyyyy, tr_x_yyy_xyyyz, tr_x_yyy_xyyz, tr_x_yyy_xyyzz, tr_x_yyy_xyzz, tr_x_yyy_xyzzz, tr_x_yyy_xzzz, tr_x_yyy_xzzzz, tr_x_yyy_yyyy, tr_x_yyy_yyyyy, tr_x_yyy_yyyyz, tr_x_yyy_yyyz, tr_x_yyy_yyyzz, tr_x_yyy_yyzz, tr_x_yyy_yyzzz, tr_x_yyy_yzzz, tr_x_yyy_yzzzz, tr_x_yyy_zzzz, tr_x_yyy_zzzzz, tr_x_yyyy_xxxxx, tr_x_yyyy_xxxxy, tr_x_yyyy_xxxxz, tr_x_yyyy_xxxyy, tr_x_yyyy_xxxyz, tr_x_yyyy_xxxzz, tr_x_yyyy_xxyyy, tr_x_yyyy_xxyyz, tr_x_yyyy_xxyzz, tr_x_yyyy_xxzzz, tr_x_yyyy_xyyyy, tr_x_yyyy_xyyyz, tr_x_yyyy_xyyzz, tr_x_yyyy_xyzzz, tr_x_yyyy_xzzzz, tr_x_yyyy_yyyyy, tr_x_yyyy_yyyyz, tr_x_yyyy_yyyzz, tr_x_yyyy_yyzzz, tr_x_yyyy_yzzzz, tr_x_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyy_xxxxx[i] = 3.0 * tr_x_yy_xxxxx[i] * fe_0 + tr_x_yyy_xxxxx[i] * pa_y[i];

        tr_x_yyyy_xxxxy[i] = 3.0 * tr_x_yy_xxxxy[i] * fe_0 + tr_x_yyy_xxxx[i] * fe_0 + tr_x_yyy_xxxxy[i] * pa_y[i];

        tr_x_yyyy_xxxxz[i] = 3.0 * tr_x_yy_xxxxz[i] * fe_0 + tr_x_yyy_xxxxz[i] * pa_y[i];

        tr_x_yyyy_xxxyy[i] = 3.0 * tr_x_yy_xxxyy[i] * fe_0 + 2.0 * tr_x_yyy_xxxy[i] * fe_0 + tr_x_yyy_xxxyy[i] * pa_y[i];

        tr_x_yyyy_xxxyz[i] = 3.0 * tr_x_yy_xxxyz[i] * fe_0 + tr_x_yyy_xxxz[i] * fe_0 + tr_x_yyy_xxxyz[i] * pa_y[i];

        tr_x_yyyy_xxxzz[i] = 3.0 * tr_x_yy_xxxzz[i] * fe_0 + tr_x_yyy_xxxzz[i] * pa_y[i];

        tr_x_yyyy_xxyyy[i] = 3.0 * tr_x_yy_xxyyy[i] * fe_0 + 3.0 * tr_x_yyy_xxyy[i] * fe_0 + tr_x_yyy_xxyyy[i] * pa_y[i];

        tr_x_yyyy_xxyyz[i] = 3.0 * tr_x_yy_xxyyz[i] * fe_0 + 2.0 * tr_x_yyy_xxyz[i] * fe_0 + tr_x_yyy_xxyyz[i] * pa_y[i];

        tr_x_yyyy_xxyzz[i] = 3.0 * tr_x_yy_xxyzz[i] * fe_0 + tr_x_yyy_xxzz[i] * fe_0 + tr_x_yyy_xxyzz[i] * pa_y[i];

        tr_x_yyyy_xxzzz[i] = 3.0 * tr_x_yy_xxzzz[i] * fe_0 + tr_x_yyy_xxzzz[i] * pa_y[i];

        tr_x_yyyy_xyyyy[i] = 3.0 * tr_x_yy_xyyyy[i] * fe_0 + 4.0 * tr_x_yyy_xyyy[i] * fe_0 + tr_x_yyy_xyyyy[i] * pa_y[i];

        tr_x_yyyy_xyyyz[i] = 3.0 * tr_x_yy_xyyyz[i] * fe_0 + 3.0 * tr_x_yyy_xyyz[i] * fe_0 + tr_x_yyy_xyyyz[i] * pa_y[i];

        tr_x_yyyy_xyyzz[i] = 3.0 * tr_x_yy_xyyzz[i] * fe_0 + 2.0 * tr_x_yyy_xyzz[i] * fe_0 + tr_x_yyy_xyyzz[i] * pa_y[i];

        tr_x_yyyy_xyzzz[i] = 3.0 * tr_x_yy_xyzzz[i] * fe_0 + tr_x_yyy_xzzz[i] * fe_0 + tr_x_yyy_xyzzz[i] * pa_y[i];

        tr_x_yyyy_xzzzz[i] = 3.0 * tr_x_yy_xzzzz[i] * fe_0 + tr_x_yyy_xzzzz[i] * pa_y[i];

        tr_x_yyyy_yyyyy[i] = 3.0 * tr_x_yy_yyyyy[i] * fe_0 + 5.0 * tr_x_yyy_yyyy[i] * fe_0 + tr_x_yyy_yyyyy[i] * pa_y[i];

        tr_x_yyyy_yyyyz[i] = 3.0 * tr_x_yy_yyyyz[i] * fe_0 + 4.0 * tr_x_yyy_yyyz[i] * fe_0 + tr_x_yyy_yyyyz[i] * pa_y[i];

        tr_x_yyyy_yyyzz[i] = 3.0 * tr_x_yy_yyyzz[i] * fe_0 + 3.0 * tr_x_yyy_yyzz[i] * fe_0 + tr_x_yyy_yyyzz[i] * pa_y[i];

        tr_x_yyyy_yyzzz[i] = 3.0 * tr_x_yy_yyzzz[i] * fe_0 + 2.0 * tr_x_yyy_yzzz[i] * fe_0 + tr_x_yyy_yyzzz[i] * pa_y[i];

        tr_x_yyyy_yzzzz[i] = 3.0 * tr_x_yy_yzzzz[i] * fe_0 + tr_x_yyy_zzzz[i] * fe_0 + tr_x_yyy_yzzzz[i] * pa_y[i];

        tr_x_yyyy_zzzzz[i] = 3.0 * tr_x_yy_zzzzz[i] * fe_0 + tr_x_yyy_zzzzz[i] * pa_y[i];
    }

    // Set up 231-252 components of targeted buffer : GH

    auto tr_x_yyyz_xxxxx = pbuffer.data(idx_dip_gh + 231);

    auto tr_x_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 232);

    auto tr_x_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 233);

    auto tr_x_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 234);

    auto tr_x_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 235);

    auto tr_x_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 236);

    auto tr_x_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 237);

    auto tr_x_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 238);

    auto tr_x_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 239);

    auto tr_x_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 240);

    auto tr_x_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 241);

    auto tr_x_yyyz_xyyyz = pbuffer.data(idx_dip_gh + 242);

    auto tr_x_yyyz_xyyzz = pbuffer.data(idx_dip_gh + 243);

    auto tr_x_yyyz_xyzzz = pbuffer.data(idx_dip_gh + 244);

    auto tr_x_yyyz_xzzzz = pbuffer.data(idx_dip_gh + 245);

    auto tr_x_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 246);

    auto tr_x_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 247);

    auto tr_x_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 248);

    auto tr_x_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 249);

    auto tr_x_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 250);

    auto tr_x_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 251);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyy_xxxxx, tr_x_yyy_xxxxy, tr_x_yyy_xxxy, tr_x_yyy_xxxyy, tr_x_yyy_xxxyz, tr_x_yyy_xxyy, tr_x_yyy_xxyyy, tr_x_yyy_xxyyz, tr_x_yyy_xxyz, tr_x_yyy_xxyzz, tr_x_yyy_xyyy, tr_x_yyy_xyyyy, tr_x_yyy_xyyyz, tr_x_yyy_xyyz, tr_x_yyy_xyyzz, tr_x_yyy_xyzz, tr_x_yyy_xyzzz, tr_x_yyy_yyyy, tr_x_yyy_yyyyy, tr_x_yyy_yyyyz, tr_x_yyy_yyyz, tr_x_yyy_yyyzz, tr_x_yyy_yyzz, tr_x_yyy_yyzzz, tr_x_yyy_yzzz, tr_x_yyy_yzzzz, tr_x_yyyz_xxxxx, tr_x_yyyz_xxxxy, tr_x_yyyz_xxxxz, tr_x_yyyz_xxxyy, tr_x_yyyz_xxxyz, tr_x_yyyz_xxxzz, tr_x_yyyz_xxyyy, tr_x_yyyz_xxyyz, tr_x_yyyz_xxyzz, tr_x_yyyz_xxzzz, tr_x_yyyz_xyyyy, tr_x_yyyz_xyyyz, tr_x_yyyz_xyyzz, tr_x_yyyz_xyzzz, tr_x_yyyz_xzzzz, tr_x_yyyz_yyyyy, tr_x_yyyz_yyyyz, tr_x_yyyz_yyyzz, tr_x_yyyz_yyzzz, tr_x_yyyz_yzzzz, tr_x_yyyz_zzzzz, tr_x_yyz_xxxxz, tr_x_yyz_xxxzz, tr_x_yyz_xxzzz, tr_x_yyz_xzzzz, tr_x_yyz_zzzzz, tr_x_yz_xxxxz, tr_x_yz_xxxzz, tr_x_yz_xxzzz, tr_x_yz_xzzzz, tr_x_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyz_xxxxx[i] = tr_x_yyy_xxxxx[i] * pa_z[i];

        tr_x_yyyz_xxxxy[i] = tr_x_yyy_xxxxy[i] * pa_z[i];

        tr_x_yyyz_xxxxz[i] = 2.0 * tr_x_yz_xxxxz[i] * fe_0 + tr_x_yyz_xxxxz[i] * pa_y[i];

        tr_x_yyyz_xxxyy[i] = tr_x_yyy_xxxyy[i] * pa_z[i];

        tr_x_yyyz_xxxyz[i] = tr_x_yyy_xxxy[i] * fe_0 + tr_x_yyy_xxxyz[i] * pa_z[i];

        tr_x_yyyz_xxxzz[i] = 2.0 * tr_x_yz_xxxzz[i] * fe_0 + tr_x_yyz_xxxzz[i] * pa_y[i];

        tr_x_yyyz_xxyyy[i] = tr_x_yyy_xxyyy[i] * pa_z[i];

        tr_x_yyyz_xxyyz[i] = tr_x_yyy_xxyy[i] * fe_0 + tr_x_yyy_xxyyz[i] * pa_z[i];

        tr_x_yyyz_xxyzz[i] = 2.0 * tr_x_yyy_xxyz[i] * fe_0 + tr_x_yyy_xxyzz[i] * pa_z[i];

        tr_x_yyyz_xxzzz[i] = 2.0 * tr_x_yz_xxzzz[i] * fe_0 + tr_x_yyz_xxzzz[i] * pa_y[i];

        tr_x_yyyz_xyyyy[i] = tr_x_yyy_xyyyy[i] * pa_z[i];

        tr_x_yyyz_xyyyz[i] = tr_x_yyy_xyyy[i] * fe_0 + tr_x_yyy_xyyyz[i] * pa_z[i];

        tr_x_yyyz_xyyzz[i] = 2.0 * tr_x_yyy_xyyz[i] * fe_0 + tr_x_yyy_xyyzz[i] * pa_z[i];

        tr_x_yyyz_xyzzz[i] = 3.0 * tr_x_yyy_xyzz[i] * fe_0 + tr_x_yyy_xyzzz[i] * pa_z[i];

        tr_x_yyyz_xzzzz[i] = 2.0 * tr_x_yz_xzzzz[i] * fe_0 + tr_x_yyz_xzzzz[i] * pa_y[i];

        tr_x_yyyz_yyyyy[i] = tr_x_yyy_yyyyy[i] * pa_z[i];

        tr_x_yyyz_yyyyz[i] = tr_x_yyy_yyyy[i] * fe_0 + tr_x_yyy_yyyyz[i] * pa_z[i];

        tr_x_yyyz_yyyzz[i] = 2.0 * tr_x_yyy_yyyz[i] * fe_0 + tr_x_yyy_yyyzz[i] * pa_z[i];

        tr_x_yyyz_yyzzz[i] = 3.0 * tr_x_yyy_yyzz[i] * fe_0 + tr_x_yyy_yyzzz[i] * pa_z[i];

        tr_x_yyyz_yzzzz[i] = 4.0 * tr_x_yyy_yzzz[i] * fe_0 + tr_x_yyy_yzzzz[i] * pa_z[i];

        tr_x_yyyz_zzzzz[i] = 2.0 * tr_x_yz_zzzzz[i] * fe_0 + tr_x_yyz_zzzzz[i] * pa_y[i];
    }

    // Set up 252-273 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yy_xxxxy, tr_x_yy_xxxyy, tr_x_yy_xxyyy, tr_x_yy_xyyyy, tr_x_yy_yyyyy, tr_x_yyz_xxxxy, tr_x_yyz_xxxyy, tr_x_yyz_xxyyy, tr_x_yyz_xyyyy, tr_x_yyz_yyyyy, tr_x_yyzz_xxxxx, tr_x_yyzz_xxxxy, tr_x_yyzz_xxxxz, tr_x_yyzz_xxxyy, tr_x_yyzz_xxxyz, tr_x_yyzz_xxxzz, tr_x_yyzz_xxyyy, tr_x_yyzz_xxyyz, tr_x_yyzz_xxyzz, tr_x_yyzz_xxzzz, tr_x_yyzz_xyyyy, tr_x_yyzz_xyyyz, tr_x_yyzz_xyyzz, tr_x_yyzz_xyzzz, tr_x_yyzz_xzzzz, tr_x_yyzz_yyyyy, tr_x_yyzz_yyyyz, tr_x_yyzz_yyyzz, tr_x_yyzz_yyzzz, tr_x_yyzz_yzzzz, tr_x_yyzz_zzzzz, tr_x_yzz_xxxxx, tr_x_yzz_xxxxz, tr_x_yzz_xxxyz, tr_x_yzz_xxxz, tr_x_yzz_xxxzz, tr_x_yzz_xxyyz, tr_x_yzz_xxyz, tr_x_yzz_xxyzz, tr_x_yzz_xxzz, tr_x_yzz_xxzzz, tr_x_yzz_xyyyz, tr_x_yzz_xyyz, tr_x_yzz_xyyzz, tr_x_yzz_xyzz, tr_x_yzz_xyzzz, tr_x_yzz_xzzz, tr_x_yzz_xzzzz, tr_x_yzz_yyyyz, tr_x_yzz_yyyz, tr_x_yzz_yyyzz, tr_x_yzz_yyzz, tr_x_yzz_yyzzz, tr_x_yzz_yzzz, tr_x_yzz_yzzzz, tr_x_yzz_zzzz, tr_x_yzz_zzzzz, tr_x_zz_xxxxx, tr_x_zz_xxxxz, tr_x_zz_xxxyz, tr_x_zz_xxxzz, tr_x_zz_xxyyz, tr_x_zz_xxyzz, tr_x_zz_xxzzz, tr_x_zz_xyyyz, tr_x_zz_xyyzz, tr_x_zz_xyzzz, tr_x_zz_xzzzz, tr_x_zz_yyyyz, tr_x_zz_yyyzz, tr_x_zz_yyzzz, tr_x_zz_yzzzz, tr_x_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzz_xxxxx[i] = tr_x_zz_xxxxx[i] * fe_0 + tr_x_yzz_xxxxx[i] * pa_y[i];

        tr_x_yyzz_xxxxy[i] = tr_x_yy_xxxxy[i] * fe_0 + tr_x_yyz_xxxxy[i] * pa_z[i];

        tr_x_yyzz_xxxxz[i] = tr_x_zz_xxxxz[i] * fe_0 + tr_x_yzz_xxxxz[i] * pa_y[i];

        tr_x_yyzz_xxxyy[i] = tr_x_yy_xxxyy[i] * fe_0 + tr_x_yyz_xxxyy[i] * pa_z[i];

        tr_x_yyzz_xxxyz[i] = tr_x_zz_xxxyz[i] * fe_0 + tr_x_yzz_xxxz[i] * fe_0 + tr_x_yzz_xxxyz[i] * pa_y[i];

        tr_x_yyzz_xxxzz[i] = tr_x_zz_xxxzz[i] * fe_0 + tr_x_yzz_xxxzz[i] * pa_y[i];

        tr_x_yyzz_xxyyy[i] = tr_x_yy_xxyyy[i] * fe_0 + tr_x_yyz_xxyyy[i] * pa_z[i];

        tr_x_yyzz_xxyyz[i] = tr_x_zz_xxyyz[i] * fe_0 + 2.0 * tr_x_yzz_xxyz[i] * fe_0 + tr_x_yzz_xxyyz[i] * pa_y[i];

        tr_x_yyzz_xxyzz[i] = tr_x_zz_xxyzz[i] * fe_0 + tr_x_yzz_xxzz[i] * fe_0 + tr_x_yzz_xxyzz[i] * pa_y[i];

        tr_x_yyzz_xxzzz[i] = tr_x_zz_xxzzz[i] * fe_0 + tr_x_yzz_xxzzz[i] * pa_y[i];

        tr_x_yyzz_xyyyy[i] = tr_x_yy_xyyyy[i] * fe_0 + tr_x_yyz_xyyyy[i] * pa_z[i];

        tr_x_yyzz_xyyyz[i] = tr_x_zz_xyyyz[i] * fe_0 + 3.0 * tr_x_yzz_xyyz[i] * fe_0 + tr_x_yzz_xyyyz[i] * pa_y[i];

        tr_x_yyzz_xyyzz[i] = tr_x_zz_xyyzz[i] * fe_0 + 2.0 * tr_x_yzz_xyzz[i] * fe_0 + tr_x_yzz_xyyzz[i] * pa_y[i];

        tr_x_yyzz_xyzzz[i] = tr_x_zz_xyzzz[i] * fe_0 + tr_x_yzz_xzzz[i] * fe_0 + tr_x_yzz_xyzzz[i] * pa_y[i];

        tr_x_yyzz_xzzzz[i] = tr_x_zz_xzzzz[i] * fe_0 + tr_x_yzz_xzzzz[i] * pa_y[i];

        tr_x_yyzz_yyyyy[i] = tr_x_yy_yyyyy[i] * fe_0 + tr_x_yyz_yyyyy[i] * pa_z[i];

        tr_x_yyzz_yyyyz[i] = tr_x_zz_yyyyz[i] * fe_0 + 4.0 * tr_x_yzz_yyyz[i] * fe_0 + tr_x_yzz_yyyyz[i] * pa_y[i];

        tr_x_yyzz_yyyzz[i] = tr_x_zz_yyyzz[i] * fe_0 + 3.0 * tr_x_yzz_yyzz[i] * fe_0 + tr_x_yzz_yyyzz[i] * pa_y[i];

        tr_x_yyzz_yyzzz[i] = tr_x_zz_yyzzz[i] * fe_0 + 2.0 * tr_x_yzz_yzzz[i] * fe_0 + tr_x_yzz_yyzzz[i] * pa_y[i];

        tr_x_yyzz_yzzzz[i] = tr_x_zz_yzzzz[i] * fe_0 + tr_x_yzz_zzzz[i] * fe_0 + tr_x_yzz_yzzzz[i] * pa_y[i];

        tr_x_yyzz_zzzzz[i] = tr_x_zz_zzzzz[i] * fe_0 + tr_x_yzz_zzzzz[i] * pa_y[i];
    }

    // Set up 273-294 components of targeted buffer : GH

    auto tr_x_yzzz_xxxxx = pbuffer.data(idx_dip_gh + 273);

    auto tr_x_yzzz_xxxxy = pbuffer.data(idx_dip_gh + 274);

    auto tr_x_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 275);

    auto tr_x_yzzz_xxxyy = pbuffer.data(idx_dip_gh + 276);

    auto tr_x_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 277);

    auto tr_x_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 278);

    auto tr_x_yzzz_xxyyy = pbuffer.data(idx_dip_gh + 279);

    auto tr_x_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 280);

    auto tr_x_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 281);

    auto tr_x_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 282);

    auto tr_x_yzzz_xyyyy = pbuffer.data(idx_dip_gh + 283);

    auto tr_x_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 284);

    auto tr_x_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 285);

    auto tr_x_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 286);

    auto tr_x_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 287);

    auto tr_x_yzzz_yyyyy = pbuffer.data(idx_dip_gh + 288);

    auto tr_x_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 289);

    auto tr_x_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 290);

    auto tr_x_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 291);

    auto tr_x_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 292);

    auto tr_x_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 293);

    #pragma omp simd aligned(pa_y, tr_x_yzzz_xxxxx, tr_x_yzzz_xxxxy, tr_x_yzzz_xxxxz, tr_x_yzzz_xxxyy, tr_x_yzzz_xxxyz, tr_x_yzzz_xxxzz, tr_x_yzzz_xxyyy, tr_x_yzzz_xxyyz, tr_x_yzzz_xxyzz, tr_x_yzzz_xxzzz, tr_x_yzzz_xyyyy, tr_x_yzzz_xyyyz, tr_x_yzzz_xyyzz, tr_x_yzzz_xyzzz, tr_x_yzzz_xzzzz, tr_x_yzzz_yyyyy, tr_x_yzzz_yyyyz, tr_x_yzzz_yyyzz, tr_x_yzzz_yyzzz, tr_x_yzzz_yzzzz, tr_x_yzzz_zzzzz, tr_x_zzz_xxxx, tr_x_zzz_xxxxx, tr_x_zzz_xxxxy, tr_x_zzz_xxxxz, tr_x_zzz_xxxy, tr_x_zzz_xxxyy, tr_x_zzz_xxxyz, tr_x_zzz_xxxz, tr_x_zzz_xxxzz, tr_x_zzz_xxyy, tr_x_zzz_xxyyy, tr_x_zzz_xxyyz, tr_x_zzz_xxyz, tr_x_zzz_xxyzz, tr_x_zzz_xxzz, tr_x_zzz_xxzzz, tr_x_zzz_xyyy, tr_x_zzz_xyyyy, tr_x_zzz_xyyyz, tr_x_zzz_xyyz, tr_x_zzz_xyyzz, tr_x_zzz_xyzz, tr_x_zzz_xyzzz, tr_x_zzz_xzzz, tr_x_zzz_xzzzz, tr_x_zzz_yyyy, tr_x_zzz_yyyyy, tr_x_zzz_yyyyz, tr_x_zzz_yyyz, tr_x_zzz_yyyzz, tr_x_zzz_yyzz, tr_x_zzz_yyzzz, tr_x_zzz_yzzz, tr_x_zzz_yzzzz, tr_x_zzz_zzzz, tr_x_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzz_xxxxx[i] = tr_x_zzz_xxxxx[i] * pa_y[i];

        tr_x_yzzz_xxxxy[i] = tr_x_zzz_xxxx[i] * fe_0 + tr_x_zzz_xxxxy[i] * pa_y[i];

        tr_x_yzzz_xxxxz[i] = tr_x_zzz_xxxxz[i] * pa_y[i];

        tr_x_yzzz_xxxyy[i] = 2.0 * tr_x_zzz_xxxy[i] * fe_0 + tr_x_zzz_xxxyy[i] * pa_y[i];

        tr_x_yzzz_xxxyz[i] = tr_x_zzz_xxxz[i] * fe_0 + tr_x_zzz_xxxyz[i] * pa_y[i];

        tr_x_yzzz_xxxzz[i] = tr_x_zzz_xxxzz[i] * pa_y[i];

        tr_x_yzzz_xxyyy[i] = 3.0 * tr_x_zzz_xxyy[i] * fe_0 + tr_x_zzz_xxyyy[i] * pa_y[i];

        tr_x_yzzz_xxyyz[i] = 2.0 * tr_x_zzz_xxyz[i] * fe_0 + tr_x_zzz_xxyyz[i] * pa_y[i];

        tr_x_yzzz_xxyzz[i] = tr_x_zzz_xxzz[i] * fe_0 + tr_x_zzz_xxyzz[i] * pa_y[i];

        tr_x_yzzz_xxzzz[i] = tr_x_zzz_xxzzz[i] * pa_y[i];

        tr_x_yzzz_xyyyy[i] = 4.0 * tr_x_zzz_xyyy[i] * fe_0 + tr_x_zzz_xyyyy[i] * pa_y[i];

        tr_x_yzzz_xyyyz[i] = 3.0 * tr_x_zzz_xyyz[i] * fe_0 + tr_x_zzz_xyyyz[i] * pa_y[i];

        tr_x_yzzz_xyyzz[i] = 2.0 * tr_x_zzz_xyzz[i] * fe_0 + tr_x_zzz_xyyzz[i] * pa_y[i];

        tr_x_yzzz_xyzzz[i] = tr_x_zzz_xzzz[i] * fe_0 + tr_x_zzz_xyzzz[i] * pa_y[i];

        tr_x_yzzz_xzzzz[i] = tr_x_zzz_xzzzz[i] * pa_y[i];

        tr_x_yzzz_yyyyy[i] = 5.0 * tr_x_zzz_yyyy[i] * fe_0 + tr_x_zzz_yyyyy[i] * pa_y[i];

        tr_x_yzzz_yyyyz[i] = 4.0 * tr_x_zzz_yyyz[i] * fe_0 + tr_x_zzz_yyyyz[i] * pa_y[i];

        tr_x_yzzz_yyyzz[i] = 3.0 * tr_x_zzz_yyzz[i] * fe_0 + tr_x_zzz_yyyzz[i] * pa_y[i];

        tr_x_yzzz_yyzzz[i] = 2.0 * tr_x_zzz_yzzz[i] * fe_0 + tr_x_zzz_yyzzz[i] * pa_y[i];

        tr_x_yzzz_yzzzz[i] = tr_x_zzz_zzzz[i] * fe_0 + tr_x_zzz_yzzzz[i] * pa_y[i];

        tr_x_yzzz_zzzzz[i] = tr_x_zzz_zzzzz[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_z, tr_x_zz_xxxxx, tr_x_zz_xxxxy, tr_x_zz_xxxxz, tr_x_zz_xxxyy, tr_x_zz_xxxyz, tr_x_zz_xxxzz, tr_x_zz_xxyyy, tr_x_zz_xxyyz, tr_x_zz_xxyzz, tr_x_zz_xxzzz, tr_x_zz_xyyyy, tr_x_zz_xyyyz, tr_x_zz_xyyzz, tr_x_zz_xyzzz, tr_x_zz_xzzzz, tr_x_zz_yyyyy, tr_x_zz_yyyyz, tr_x_zz_yyyzz, tr_x_zz_yyzzz, tr_x_zz_yzzzz, tr_x_zz_zzzzz, tr_x_zzz_xxxx, tr_x_zzz_xxxxx, tr_x_zzz_xxxxy, tr_x_zzz_xxxxz, tr_x_zzz_xxxy, tr_x_zzz_xxxyy, tr_x_zzz_xxxyz, tr_x_zzz_xxxz, tr_x_zzz_xxxzz, tr_x_zzz_xxyy, tr_x_zzz_xxyyy, tr_x_zzz_xxyyz, tr_x_zzz_xxyz, tr_x_zzz_xxyzz, tr_x_zzz_xxzz, tr_x_zzz_xxzzz, tr_x_zzz_xyyy, tr_x_zzz_xyyyy, tr_x_zzz_xyyyz, tr_x_zzz_xyyz, tr_x_zzz_xyyzz, tr_x_zzz_xyzz, tr_x_zzz_xyzzz, tr_x_zzz_xzzz, tr_x_zzz_xzzzz, tr_x_zzz_yyyy, tr_x_zzz_yyyyy, tr_x_zzz_yyyyz, tr_x_zzz_yyyz, tr_x_zzz_yyyzz, tr_x_zzz_yyzz, tr_x_zzz_yyzzz, tr_x_zzz_yzzz, tr_x_zzz_yzzzz, tr_x_zzz_zzzz, tr_x_zzz_zzzzz, tr_x_zzzz_xxxxx, tr_x_zzzz_xxxxy, tr_x_zzzz_xxxxz, tr_x_zzzz_xxxyy, tr_x_zzzz_xxxyz, tr_x_zzzz_xxxzz, tr_x_zzzz_xxyyy, tr_x_zzzz_xxyyz, tr_x_zzzz_xxyzz, tr_x_zzzz_xxzzz, tr_x_zzzz_xyyyy, tr_x_zzzz_xyyyz, tr_x_zzzz_xyyzz, tr_x_zzzz_xyzzz, tr_x_zzzz_xzzzz, tr_x_zzzz_yyyyy, tr_x_zzzz_yyyyz, tr_x_zzzz_yyyzz, tr_x_zzzz_yyzzz, tr_x_zzzz_yzzzz, tr_x_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzz_xxxxx[i] = 3.0 * tr_x_zz_xxxxx[i] * fe_0 + tr_x_zzz_xxxxx[i] * pa_z[i];

        tr_x_zzzz_xxxxy[i] = 3.0 * tr_x_zz_xxxxy[i] * fe_0 + tr_x_zzz_xxxxy[i] * pa_z[i];

        tr_x_zzzz_xxxxz[i] = 3.0 * tr_x_zz_xxxxz[i] * fe_0 + tr_x_zzz_xxxx[i] * fe_0 + tr_x_zzz_xxxxz[i] * pa_z[i];

        tr_x_zzzz_xxxyy[i] = 3.0 * tr_x_zz_xxxyy[i] * fe_0 + tr_x_zzz_xxxyy[i] * pa_z[i];

        tr_x_zzzz_xxxyz[i] = 3.0 * tr_x_zz_xxxyz[i] * fe_0 + tr_x_zzz_xxxy[i] * fe_0 + tr_x_zzz_xxxyz[i] * pa_z[i];

        tr_x_zzzz_xxxzz[i] = 3.0 * tr_x_zz_xxxzz[i] * fe_0 + 2.0 * tr_x_zzz_xxxz[i] * fe_0 + tr_x_zzz_xxxzz[i] * pa_z[i];

        tr_x_zzzz_xxyyy[i] = 3.0 * tr_x_zz_xxyyy[i] * fe_0 + tr_x_zzz_xxyyy[i] * pa_z[i];

        tr_x_zzzz_xxyyz[i] = 3.0 * tr_x_zz_xxyyz[i] * fe_0 + tr_x_zzz_xxyy[i] * fe_0 + tr_x_zzz_xxyyz[i] * pa_z[i];

        tr_x_zzzz_xxyzz[i] = 3.0 * tr_x_zz_xxyzz[i] * fe_0 + 2.0 * tr_x_zzz_xxyz[i] * fe_0 + tr_x_zzz_xxyzz[i] * pa_z[i];

        tr_x_zzzz_xxzzz[i] = 3.0 * tr_x_zz_xxzzz[i] * fe_0 + 3.0 * tr_x_zzz_xxzz[i] * fe_0 + tr_x_zzz_xxzzz[i] * pa_z[i];

        tr_x_zzzz_xyyyy[i] = 3.0 * tr_x_zz_xyyyy[i] * fe_0 + tr_x_zzz_xyyyy[i] * pa_z[i];

        tr_x_zzzz_xyyyz[i] = 3.0 * tr_x_zz_xyyyz[i] * fe_0 + tr_x_zzz_xyyy[i] * fe_0 + tr_x_zzz_xyyyz[i] * pa_z[i];

        tr_x_zzzz_xyyzz[i] = 3.0 * tr_x_zz_xyyzz[i] * fe_0 + 2.0 * tr_x_zzz_xyyz[i] * fe_0 + tr_x_zzz_xyyzz[i] * pa_z[i];

        tr_x_zzzz_xyzzz[i] = 3.0 * tr_x_zz_xyzzz[i] * fe_0 + 3.0 * tr_x_zzz_xyzz[i] * fe_0 + tr_x_zzz_xyzzz[i] * pa_z[i];

        tr_x_zzzz_xzzzz[i] = 3.0 * tr_x_zz_xzzzz[i] * fe_0 + 4.0 * tr_x_zzz_xzzz[i] * fe_0 + tr_x_zzz_xzzzz[i] * pa_z[i];

        tr_x_zzzz_yyyyy[i] = 3.0 * tr_x_zz_yyyyy[i] * fe_0 + tr_x_zzz_yyyyy[i] * pa_z[i];

        tr_x_zzzz_yyyyz[i] = 3.0 * tr_x_zz_yyyyz[i] * fe_0 + tr_x_zzz_yyyy[i] * fe_0 + tr_x_zzz_yyyyz[i] * pa_z[i];

        tr_x_zzzz_yyyzz[i] = 3.0 * tr_x_zz_yyyzz[i] * fe_0 + 2.0 * tr_x_zzz_yyyz[i] * fe_0 + tr_x_zzz_yyyzz[i] * pa_z[i];

        tr_x_zzzz_yyzzz[i] = 3.0 * tr_x_zz_yyzzz[i] * fe_0 + 3.0 * tr_x_zzz_yyzz[i] * fe_0 + tr_x_zzz_yyzzz[i] * pa_z[i];

        tr_x_zzzz_yzzzz[i] = 3.0 * tr_x_zz_yzzzz[i] * fe_0 + 4.0 * tr_x_zzz_yzzz[i] * fe_0 + tr_x_zzz_yzzzz[i] * pa_z[i];

        tr_x_zzzz_zzzzz[i] = 3.0 * tr_x_zz_zzzzz[i] * fe_0 + 5.0 * tr_x_zzz_zzzz[i] * fe_0 + tr_x_zzz_zzzzz[i] * pa_z[i];
    }

    // Set up 315-336 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_y_xx_xxxxx, tr_y_xx_xxxxy, tr_y_xx_xxxxz, tr_y_xx_xxxyy, tr_y_xx_xxxyz, tr_y_xx_xxxzz, tr_y_xx_xxyyy, tr_y_xx_xxyyz, tr_y_xx_xxyzz, tr_y_xx_xxzzz, tr_y_xx_xyyyy, tr_y_xx_xyyyz, tr_y_xx_xyyzz, tr_y_xx_xyzzz, tr_y_xx_xzzzz, tr_y_xx_yyyyy, tr_y_xx_yyyyz, tr_y_xx_yyyzz, tr_y_xx_yyzzz, tr_y_xx_yzzzz, tr_y_xx_zzzzz, tr_y_xxx_xxxx, tr_y_xxx_xxxxx, tr_y_xxx_xxxxy, tr_y_xxx_xxxxz, tr_y_xxx_xxxy, tr_y_xxx_xxxyy, tr_y_xxx_xxxyz, tr_y_xxx_xxxz, tr_y_xxx_xxxzz, tr_y_xxx_xxyy, tr_y_xxx_xxyyy, tr_y_xxx_xxyyz, tr_y_xxx_xxyz, tr_y_xxx_xxyzz, tr_y_xxx_xxzz, tr_y_xxx_xxzzz, tr_y_xxx_xyyy, tr_y_xxx_xyyyy, tr_y_xxx_xyyyz, tr_y_xxx_xyyz, tr_y_xxx_xyyzz, tr_y_xxx_xyzz, tr_y_xxx_xyzzz, tr_y_xxx_xzzz, tr_y_xxx_xzzzz, tr_y_xxx_yyyy, tr_y_xxx_yyyyy, tr_y_xxx_yyyyz, tr_y_xxx_yyyz, tr_y_xxx_yyyzz, tr_y_xxx_yyzz, tr_y_xxx_yyzzz, tr_y_xxx_yzzz, tr_y_xxx_yzzzz, tr_y_xxx_zzzz, tr_y_xxx_zzzzz, tr_y_xxxx_xxxxx, tr_y_xxxx_xxxxy, tr_y_xxxx_xxxxz, tr_y_xxxx_xxxyy, tr_y_xxxx_xxxyz, tr_y_xxxx_xxxzz, tr_y_xxxx_xxyyy, tr_y_xxxx_xxyyz, tr_y_xxxx_xxyzz, tr_y_xxxx_xxzzz, tr_y_xxxx_xyyyy, tr_y_xxxx_xyyyz, tr_y_xxxx_xyyzz, tr_y_xxxx_xyzzz, tr_y_xxxx_xzzzz, tr_y_xxxx_yyyyy, tr_y_xxxx_yyyyz, tr_y_xxxx_yyyzz, tr_y_xxxx_yyzzz, tr_y_xxxx_yzzzz, tr_y_xxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxx_xxxxx[i] = 3.0 * tr_y_xx_xxxxx[i] * fe_0 + 5.0 * tr_y_xxx_xxxx[i] * fe_0 + tr_y_xxx_xxxxx[i] * pa_x[i];

        tr_y_xxxx_xxxxy[i] = 3.0 * tr_y_xx_xxxxy[i] * fe_0 + 4.0 * tr_y_xxx_xxxy[i] * fe_0 + tr_y_xxx_xxxxy[i] * pa_x[i];

        tr_y_xxxx_xxxxz[i] = 3.0 * tr_y_xx_xxxxz[i] * fe_0 + 4.0 * tr_y_xxx_xxxz[i] * fe_0 + tr_y_xxx_xxxxz[i] * pa_x[i];

        tr_y_xxxx_xxxyy[i] = 3.0 * tr_y_xx_xxxyy[i] * fe_0 + 3.0 * tr_y_xxx_xxyy[i] * fe_0 + tr_y_xxx_xxxyy[i] * pa_x[i];

        tr_y_xxxx_xxxyz[i] = 3.0 * tr_y_xx_xxxyz[i] * fe_0 + 3.0 * tr_y_xxx_xxyz[i] * fe_0 + tr_y_xxx_xxxyz[i] * pa_x[i];

        tr_y_xxxx_xxxzz[i] = 3.0 * tr_y_xx_xxxzz[i] * fe_0 + 3.0 * tr_y_xxx_xxzz[i] * fe_0 + tr_y_xxx_xxxzz[i] * pa_x[i];

        tr_y_xxxx_xxyyy[i] = 3.0 * tr_y_xx_xxyyy[i] * fe_0 + 2.0 * tr_y_xxx_xyyy[i] * fe_0 + tr_y_xxx_xxyyy[i] * pa_x[i];

        tr_y_xxxx_xxyyz[i] = 3.0 * tr_y_xx_xxyyz[i] * fe_0 + 2.0 * tr_y_xxx_xyyz[i] * fe_0 + tr_y_xxx_xxyyz[i] * pa_x[i];

        tr_y_xxxx_xxyzz[i] = 3.0 * tr_y_xx_xxyzz[i] * fe_0 + 2.0 * tr_y_xxx_xyzz[i] * fe_0 + tr_y_xxx_xxyzz[i] * pa_x[i];

        tr_y_xxxx_xxzzz[i] = 3.0 * tr_y_xx_xxzzz[i] * fe_0 + 2.0 * tr_y_xxx_xzzz[i] * fe_0 + tr_y_xxx_xxzzz[i] * pa_x[i];

        tr_y_xxxx_xyyyy[i] = 3.0 * tr_y_xx_xyyyy[i] * fe_0 + tr_y_xxx_yyyy[i] * fe_0 + tr_y_xxx_xyyyy[i] * pa_x[i];

        tr_y_xxxx_xyyyz[i] = 3.0 * tr_y_xx_xyyyz[i] * fe_0 + tr_y_xxx_yyyz[i] * fe_0 + tr_y_xxx_xyyyz[i] * pa_x[i];

        tr_y_xxxx_xyyzz[i] = 3.0 * tr_y_xx_xyyzz[i] * fe_0 + tr_y_xxx_yyzz[i] * fe_0 + tr_y_xxx_xyyzz[i] * pa_x[i];

        tr_y_xxxx_xyzzz[i] = 3.0 * tr_y_xx_xyzzz[i] * fe_0 + tr_y_xxx_yzzz[i] * fe_0 + tr_y_xxx_xyzzz[i] * pa_x[i];

        tr_y_xxxx_xzzzz[i] = 3.0 * tr_y_xx_xzzzz[i] * fe_0 + tr_y_xxx_zzzz[i] * fe_0 + tr_y_xxx_xzzzz[i] * pa_x[i];

        tr_y_xxxx_yyyyy[i] = 3.0 * tr_y_xx_yyyyy[i] * fe_0 + tr_y_xxx_yyyyy[i] * pa_x[i];

        tr_y_xxxx_yyyyz[i] = 3.0 * tr_y_xx_yyyyz[i] * fe_0 + tr_y_xxx_yyyyz[i] * pa_x[i];

        tr_y_xxxx_yyyzz[i] = 3.0 * tr_y_xx_yyyzz[i] * fe_0 + tr_y_xxx_yyyzz[i] * pa_x[i];

        tr_y_xxxx_yyzzz[i] = 3.0 * tr_y_xx_yyzzz[i] * fe_0 + tr_y_xxx_yyzzz[i] * pa_x[i];

        tr_y_xxxx_yzzzz[i] = 3.0 * tr_y_xx_yzzzz[i] * fe_0 + tr_y_xxx_yzzzz[i] * pa_x[i];

        tr_y_xxxx_zzzzz[i] = 3.0 * tr_y_xx_zzzzz[i] * fe_0 + tr_y_xxx_zzzzz[i] * pa_x[i];
    }

    // Set up 336-357 components of targeted buffer : GH

    auto tr_y_xxxy_xxxxx = pbuffer.data(idx_dip_gh + 336);

    auto tr_y_xxxy_xxxxy = pbuffer.data(idx_dip_gh + 337);

    auto tr_y_xxxy_xxxxz = pbuffer.data(idx_dip_gh + 338);

    auto tr_y_xxxy_xxxyy = pbuffer.data(idx_dip_gh + 339);

    auto tr_y_xxxy_xxxyz = pbuffer.data(idx_dip_gh + 340);

    auto tr_y_xxxy_xxxzz = pbuffer.data(idx_dip_gh + 341);

    auto tr_y_xxxy_xxyyy = pbuffer.data(idx_dip_gh + 342);

    auto tr_y_xxxy_xxyyz = pbuffer.data(idx_dip_gh + 343);

    auto tr_y_xxxy_xxyzz = pbuffer.data(idx_dip_gh + 344);

    auto tr_y_xxxy_xxzzz = pbuffer.data(idx_dip_gh + 345);

    auto tr_y_xxxy_xyyyy = pbuffer.data(idx_dip_gh + 346);

    auto tr_y_xxxy_xyyyz = pbuffer.data(idx_dip_gh + 347);

    auto tr_y_xxxy_xyyzz = pbuffer.data(idx_dip_gh + 348);

    auto tr_y_xxxy_xyzzz = pbuffer.data(idx_dip_gh + 349);

    auto tr_y_xxxy_xzzzz = pbuffer.data(idx_dip_gh + 350);

    auto tr_y_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 351);

    auto tr_y_xxxy_yyyyz = pbuffer.data(idx_dip_gh + 352);

    auto tr_y_xxxy_yyyzz = pbuffer.data(idx_dip_gh + 353);

    auto tr_y_xxxy_yyzzz = pbuffer.data(idx_dip_gh + 354);

    auto tr_y_xxxy_yzzzz = pbuffer.data(idx_dip_gh + 355);

    auto tr_y_xxxy_zzzzz = pbuffer.data(idx_dip_gh + 356);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xxx_xxxxx, tr_y_xxx_xxxxz, tr_y_xxx_xxxzz, tr_y_xxx_xxzzz, tr_y_xxx_xzzzz, tr_y_xxxy_xxxxx, tr_y_xxxy_xxxxy, tr_y_xxxy_xxxxz, tr_y_xxxy_xxxyy, tr_y_xxxy_xxxyz, tr_y_xxxy_xxxzz, tr_y_xxxy_xxyyy, tr_y_xxxy_xxyyz, tr_y_xxxy_xxyzz, tr_y_xxxy_xxzzz, tr_y_xxxy_xyyyy, tr_y_xxxy_xyyyz, tr_y_xxxy_xyyzz, tr_y_xxxy_xyzzz, tr_y_xxxy_xzzzz, tr_y_xxxy_yyyyy, tr_y_xxxy_yyyyz, tr_y_xxxy_yyyzz, tr_y_xxxy_yyzzz, tr_y_xxxy_yzzzz, tr_y_xxxy_zzzzz, tr_y_xxy_xxxxy, tr_y_xxy_xxxy, tr_y_xxy_xxxyy, tr_y_xxy_xxxyz, tr_y_xxy_xxyy, tr_y_xxy_xxyyy, tr_y_xxy_xxyyz, tr_y_xxy_xxyz, tr_y_xxy_xxyzz, tr_y_xxy_xyyy, tr_y_xxy_xyyyy, tr_y_xxy_xyyyz, tr_y_xxy_xyyz, tr_y_xxy_xyyzz, tr_y_xxy_xyzz, tr_y_xxy_xyzzz, tr_y_xxy_yyyy, tr_y_xxy_yyyyy, tr_y_xxy_yyyyz, tr_y_xxy_yyyz, tr_y_xxy_yyyzz, tr_y_xxy_yyzz, tr_y_xxy_yyzzz, tr_y_xxy_yzzz, tr_y_xxy_yzzzz, tr_y_xxy_zzzzz, tr_y_xy_xxxxy, tr_y_xy_xxxyy, tr_y_xy_xxxyz, tr_y_xy_xxyyy, tr_y_xy_xxyyz, tr_y_xy_xxyzz, tr_y_xy_xyyyy, tr_y_xy_xyyyz, tr_y_xy_xyyzz, tr_y_xy_xyzzz, tr_y_xy_yyyyy, tr_y_xy_yyyyz, tr_y_xy_yyyzz, tr_y_xy_yyzzz, tr_y_xy_yzzzz, tr_y_xy_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxz, ts_xxx_xxxzz, ts_xxx_xxzzz, ts_xxx_xzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxy_xxxxx[i] = ts_xxx_xxxxx[i] * fe_0 + tr_y_xxx_xxxxx[i] * pa_y[i];

        tr_y_xxxy_xxxxy[i] = 2.0 * tr_y_xy_xxxxy[i] * fe_0 + 4.0 * tr_y_xxy_xxxy[i] * fe_0 + tr_y_xxy_xxxxy[i] * pa_x[i];

        tr_y_xxxy_xxxxz[i] = ts_xxx_xxxxz[i] * fe_0 + tr_y_xxx_xxxxz[i] * pa_y[i];

        tr_y_xxxy_xxxyy[i] = 2.0 * tr_y_xy_xxxyy[i] * fe_0 + 3.0 * tr_y_xxy_xxyy[i] * fe_0 + tr_y_xxy_xxxyy[i] * pa_x[i];

        tr_y_xxxy_xxxyz[i] = 2.0 * tr_y_xy_xxxyz[i] * fe_0 + 3.0 * tr_y_xxy_xxyz[i] * fe_0 + tr_y_xxy_xxxyz[i] * pa_x[i];

        tr_y_xxxy_xxxzz[i] = ts_xxx_xxxzz[i] * fe_0 + tr_y_xxx_xxxzz[i] * pa_y[i];

        tr_y_xxxy_xxyyy[i] = 2.0 * tr_y_xy_xxyyy[i] * fe_0 + 2.0 * tr_y_xxy_xyyy[i] * fe_0 + tr_y_xxy_xxyyy[i] * pa_x[i];

        tr_y_xxxy_xxyyz[i] = 2.0 * tr_y_xy_xxyyz[i] * fe_0 + 2.0 * tr_y_xxy_xyyz[i] * fe_0 + tr_y_xxy_xxyyz[i] * pa_x[i];

        tr_y_xxxy_xxyzz[i] = 2.0 * tr_y_xy_xxyzz[i] * fe_0 + 2.0 * tr_y_xxy_xyzz[i] * fe_0 + tr_y_xxy_xxyzz[i] * pa_x[i];

        tr_y_xxxy_xxzzz[i] = ts_xxx_xxzzz[i] * fe_0 + tr_y_xxx_xxzzz[i] * pa_y[i];

        tr_y_xxxy_xyyyy[i] = 2.0 * tr_y_xy_xyyyy[i] * fe_0 + tr_y_xxy_yyyy[i] * fe_0 + tr_y_xxy_xyyyy[i] * pa_x[i];

        tr_y_xxxy_xyyyz[i] = 2.0 * tr_y_xy_xyyyz[i] * fe_0 + tr_y_xxy_yyyz[i] * fe_0 + tr_y_xxy_xyyyz[i] * pa_x[i];

        tr_y_xxxy_xyyzz[i] = 2.0 * tr_y_xy_xyyzz[i] * fe_0 + tr_y_xxy_yyzz[i] * fe_0 + tr_y_xxy_xyyzz[i] * pa_x[i];

        tr_y_xxxy_xyzzz[i] = 2.0 * tr_y_xy_xyzzz[i] * fe_0 + tr_y_xxy_yzzz[i] * fe_0 + tr_y_xxy_xyzzz[i] * pa_x[i];

        tr_y_xxxy_xzzzz[i] = ts_xxx_xzzzz[i] * fe_0 + tr_y_xxx_xzzzz[i] * pa_y[i];

        tr_y_xxxy_yyyyy[i] = 2.0 * tr_y_xy_yyyyy[i] * fe_0 + tr_y_xxy_yyyyy[i] * pa_x[i];

        tr_y_xxxy_yyyyz[i] = 2.0 * tr_y_xy_yyyyz[i] * fe_0 + tr_y_xxy_yyyyz[i] * pa_x[i];

        tr_y_xxxy_yyyzz[i] = 2.0 * tr_y_xy_yyyzz[i] * fe_0 + tr_y_xxy_yyyzz[i] * pa_x[i];

        tr_y_xxxy_yyzzz[i] = 2.0 * tr_y_xy_yyzzz[i] * fe_0 + tr_y_xxy_yyzzz[i] * pa_x[i];

        tr_y_xxxy_yzzzz[i] = 2.0 * tr_y_xy_yzzzz[i] * fe_0 + tr_y_xxy_yzzzz[i] * pa_x[i];

        tr_y_xxxy_zzzzz[i] = 2.0 * tr_y_xy_zzzzz[i] * fe_0 + tr_y_xxy_zzzzz[i] * pa_x[i];
    }

    // Set up 357-378 components of targeted buffer : GH

    auto tr_y_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 357);

    auto tr_y_xxxz_xxxxy = pbuffer.data(idx_dip_gh + 358);

    auto tr_y_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 359);

    auto tr_y_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 360);

    auto tr_y_xxxz_xxxyz = pbuffer.data(idx_dip_gh + 361);

    auto tr_y_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 362);

    auto tr_y_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 363);

    auto tr_y_xxxz_xxyyz = pbuffer.data(idx_dip_gh + 364);

    auto tr_y_xxxz_xxyzz = pbuffer.data(idx_dip_gh + 365);

    auto tr_y_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 366);

    auto tr_y_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 367);

    auto tr_y_xxxz_xyyyz = pbuffer.data(idx_dip_gh + 368);

    auto tr_y_xxxz_xyyzz = pbuffer.data(idx_dip_gh + 369);

    auto tr_y_xxxz_xyzzz = pbuffer.data(idx_dip_gh + 370);

    auto tr_y_xxxz_xzzzz = pbuffer.data(idx_dip_gh + 371);

    auto tr_y_xxxz_yyyyy = pbuffer.data(idx_dip_gh + 372);

    auto tr_y_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 373);

    auto tr_y_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 374);

    auto tr_y_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 375);

    auto tr_y_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 376);

    auto tr_y_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 377);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxx_xxxx, tr_y_xxx_xxxxx, tr_y_xxx_xxxxy, tr_y_xxx_xxxxz, tr_y_xxx_xxxy, tr_y_xxx_xxxyy, tr_y_xxx_xxxyz, tr_y_xxx_xxxz, tr_y_xxx_xxxzz, tr_y_xxx_xxyy, tr_y_xxx_xxyyy, tr_y_xxx_xxyyz, tr_y_xxx_xxyz, tr_y_xxx_xxyzz, tr_y_xxx_xxzz, tr_y_xxx_xxzzz, tr_y_xxx_xyyy, tr_y_xxx_xyyyy, tr_y_xxx_xyyyz, tr_y_xxx_xyyz, tr_y_xxx_xyyzz, tr_y_xxx_xyzz, tr_y_xxx_xyzzz, tr_y_xxx_xzzz, tr_y_xxx_xzzzz, tr_y_xxx_yyyyy, tr_y_xxxz_xxxxx, tr_y_xxxz_xxxxy, tr_y_xxxz_xxxxz, tr_y_xxxz_xxxyy, tr_y_xxxz_xxxyz, tr_y_xxxz_xxxzz, tr_y_xxxz_xxyyy, tr_y_xxxz_xxyyz, tr_y_xxxz_xxyzz, tr_y_xxxz_xxzzz, tr_y_xxxz_xyyyy, tr_y_xxxz_xyyyz, tr_y_xxxz_xyyzz, tr_y_xxxz_xyzzz, tr_y_xxxz_xzzzz, tr_y_xxxz_yyyyy, tr_y_xxxz_yyyyz, tr_y_xxxz_yyyzz, tr_y_xxxz_yyzzz, tr_y_xxxz_yzzzz, tr_y_xxxz_zzzzz, tr_y_xxz_yyyyz, tr_y_xxz_yyyzz, tr_y_xxz_yyzzz, tr_y_xxz_yzzzz, tr_y_xxz_zzzzz, tr_y_xz_yyyyz, tr_y_xz_yyyzz, tr_y_xz_yyzzz, tr_y_xz_yzzzz, tr_y_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxz_xxxxx[i] = tr_y_xxx_xxxxx[i] * pa_z[i];

        tr_y_xxxz_xxxxy[i] = tr_y_xxx_xxxxy[i] * pa_z[i];

        tr_y_xxxz_xxxxz[i] = tr_y_xxx_xxxx[i] * fe_0 + tr_y_xxx_xxxxz[i] * pa_z[i];

        tr_y_xxxz_xxxyy[i] = tr_y_xxx_xxxyy[i] * pa_z[i];

        tr_y_xxxz_xxxyz[i] = tr_y_xxx_xxxy[i] * fe_0 + tr_y_xxx_xxxyz[i] * pa_z[i];

        tr_y_xxxz_xxxzz[i] = 2.0 * tr_y_xxx_xxxz[i] * fe_0 + tr_y_xxx_xxxzz[i] * pa_z[i];

        tr_y_xxxz_xxyyy[i] = tr_y_xxx_xxyyy[i] * pa_z[i];

        tr_y_xxxz_xxyyz[i] = tr_y_xxx_xxyy[i] * fe_0 + tr_y_xxx_xxyyz[i] * pa_z[i];

        tr_y_xxxz_xxyzz[i] = 2.0 * tr_y_xxx_xxyz[i] * fe_0 + tr_y_xxx_xxyzz[i] * pa_z[i];

        tr_y_xxxz_xxzzz[i] = 3.0 * tr_y_xxx_xxzz[i] * fe_0 + tr_y_xxx_xxzzz[i] * pa_z[i];

        tr_y_xxxz_xyyyy[i] = tr_y_xxx_xyyyy[i] * pa_z[i];

        tr_y_xxxz_xyyyz[i] = tr_y_xxx_xyyy[i] * fe_0 + tr_y_xxx_xyyyz[i] * pa_z[i];

        tr_y_xxxz_xyyzz[i] = 2.0 * tr_y_xxx_xyyz[i] * fe_0 + tr_y_xxx_xyyzz[i] * pa_z[i];

        tr_y_xxxz_xyzzz[i] = 3.0 * tr_y_xxx_xyzz[i] * fe_0 + tr_y_xxx_xyzzz[i] * pa_z[i];

        tr_y_xxxz_xzzzz[i] = 4.0 * tr_y_xxx_xzzz[i] * fe_0 + tr_y_xxx_xzzzz[i] * pa_z[i];

        tr_y_xxxz_yyyyy[i] = tr_y_xxx_yyyyy[i] * pa_z[i];

        tr_y_xxxz_yyyyz[i] = 2.0 * tr_y_xz_yyyyz[i] * fe_0 + tr_y_xxz_yyyyz[i] * pa_x[i];

        tr_y_xxxz_yyyzz[i] = 2.0 * tr_y_xz_yyyzz[i] * fe_0 + tr_y_xxz_yyyzz[i] * pa_x[i];

        tr_y_xxxz_yyzzz[i] = 2.0 * tr_y_xz_yyzzz[i] * fe_0 + tr_y_xxz_yyzzz[i] * pa_x[i];

        tr_y_xxxz_yzzzz[i] = 2.0 * tr_y_xz_yzzzz[i] * fe_0 + tr_y_xxz_yzzzz[i] * pa_x[i];

        tr_y_xxxz_zzzzz[i] = 2.0 * tr_y_xz_zzzzz[i] * fe_0 + tr_y_xxz_zzzzz[i] * pa_x[i];
    }

    // Set up 378-399 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_y_xxyy_xxxxx, tr_y_xxyy_xxxxy, tr_y_xxyy_xxxxz, tr_y_xxyy_xxxyy, tr_y_xxyy_xxxyz, tr_y_xxyy_xxxzz, tr_y_xxyy_xxyyy, tr_y_xxyy_xxyyz, tr_y_xxyy_xxyzz, tr_y_xxyy_xxzzz, tr_y_xxyy_xyyyy, tr_y_xxyy_xyyyz, tr_y_xxyy_xyyzz, tr_y_xxyy_xyzzz, tr_y_xxyy_xzzzz, tr_y_xxyy_yyyyy, tr_y_xxyy_yyyyz, tr_y_xxyy_yyyzz, tr_y_xxyy_yyzzz, tr_y_xxyy_yzzzz, tr_y_xxyy_zzzzz, tr_y_xyy_xxxx, tr_y_xyy_xxxxx, tr_y_xyy_xxxxy, tr_y_xyy_xxxxz, tr_y_xyy_xxxy, tr_y_xyy_xxxyy, tr_y_xyy_xxxyz, tr_y_xyy_xxxz, tr_y_xyy_xxxzz, tr_y_xyy_xxyy, tr_y_xyy_xxyyy, tr_y_xyy_xxyyz, tr_y_xyy_xxyz, tr_y_xyy_xxyzz, tr_y_xyy_xxzz, tr_y_xyy_xxzzz, tr_y_xyy_xyyy, tr_y_xyy_xyyyy, tr_y_xyy_xyyyz, tr_y_xyy_xyyz, tr_y_xyy_xyyzz, tr_y_xyy_xyzz, tr_y_xyy_xyzzz, tr_y_xyy_xzzz, tr_y_xyy_xzzzz, tr_y_xyy_yyyy, tr_y_xyy_yyyyy, tr_y_xyy_yyyyz, tr_y_xyy_yyyz, tr_y_xyy_yyyzz, tr_y_xyy_yyzz, tr_y_xyy_yyzzz, tr_y_xyy_yzzz, tr_y_xyy_yzzzz, tr_y_xyy_zzzz, tr_y_xyy_zzzzz, tr_y_yy_xxxxx, tr_y_yy_xxxxy, tr_y_yy_xxxxz, tr_y_yy_xxxyy, tr_y_yy_xxxyz, tr_y_yy_xxxzz, tr_y_yy_xxyyy, tr_y_yy_xxyyz, tr_y_yy_xxyzz, tr_y_yy_xxzzz, tr_y_yy_xyyyy, tr_y_yy_xyyyz, tr_y_yy_xyyzz, tr_y_yy_xyzzz, tr_y_yy_xzzzz, tr_y_yy_yyyyy, tr_y_yy_yyyyz, tr_y_yy_yyyzz, tr_y_yy_yyzzz, tr_y_yy_yzzzz, tr_y_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyy_xxxxx[i] = tr_y_yy_xxxxx[i] * fe_0 + 5.0 * tr_y_xyy_xxxx[i] * fe_0 + tr_y_xyy_xxxxx[i] * pa_x[i];

        tr_y_xxyy_xxxxy[i] = tr_y_yy_xxxxy[i] * fe_0 + 4.0 * tr_y_xyy_xxxy[i] * fe_0 + tr_y_xyy_xxxxy[i] * pa_x[i];

        tr_y_xxyy_xxxxz[i] = tr_y_yy_xxxxz[i] * fe_0 + 4.0 * tr_y_xyy_xxxz[i] * fe_0 + tr_y_xyy_xxxxz[i] * pa_x[i];

        tr_y_xxyy_xxxyy[i] = tr_y_yy_xxxyy[i] * fe_0 + 3.0 * tr_y_xyy_xxyy[i] * fe_0 + tr_y_xyy_xxxyy[i] * pa_x[i];

        tr_y_xxyy_xxxyz[i] = tr_y_yy_xxxyz[i] * fe_0 + 3.0 * tr_y_xyy_xxyz[i] * fe_0 + tr_y_xyy_xxxyz[i] * pa_x[i];

        tr_y_xxyy_xxxzz[i] = tr_y_yy_xxxzz[i] * fe_0 + 3.0 * tr_y_xyy_xxzz[i] * fe_0 + tr_y_xyy_xxxzz[i] * pa_x[i];

        tr_y_xxyy_xxyyy[i] = tr_y_yy_xxyyy[i] * fe_0 + 2.0 * tr_y_xyy_xyyy[i] * fe_0 + tr_y_xyy_xxyyy[i] * pa_x[i];

        tr_y_xxyy_xxyyz[i] = tr_y_yy_xxyyz[i] * fe_0 + 2.0 * tr_y_xyy_xyyz[i] * fe_0 + tr_y_xyy_xxyyz[i] * pa_x[i];

        tr_y_xxyy_xxyzz[i] = tr_y_yy_xxyzz[i] * fe_0 + 2.0 * tr_y_xyy_xyzz[i] * fe_0 + tr_y_xyy_xxyzz[i] * pa_x[i];

        tr_y_xxyy_xxzzz[i] = tr_y_yy_xxzzz[i] * fe_0 + 2.0 * tr_y_xyy_xzzz[i] * fe_0 + tr_y_xyy_xxzzz[i] * pa_x[i];

        tr_y_xxyy_xyyyy[i] = tr_y_yy_xyyyy[i] * fe_0 + tr_y_xyy_yyyy[i] * fe_0 + tr_y_xyy_xyyyy[i] * pa_x[i];

        tr_y_xxyy_xyyyz[i] = tr_y_yy_xyyyz[i] * fe_0 + tr_y_xyy_yyyz[i] * fe_0 + tr_y_xyy_xyyyz[i] * pa_x[i];

        tr_y_xxyy_xyyzz[i] = tr_y_yy_xyyzz[i] * fe_0 + tr_y_xyy_yyzz[i] * fe_0 + tr_y_xyy_xyyzz[i] * pa_x[i];

        tr_y_xxyy_xyzzz[i] = tr_y_yy_xyzzz[i] * fe_0 + tr_y_xyy_yzzz[i] * fe_0 + tr_y_xyy_xyzzz[i] * pa_x[i];

        tr_y_xxyy_xzzzz[i] = tr_y_yy_xzzzz[i] * fe_0 + tr_y_xyy_zzzz[i] * fe_0 + tr_y_xyy_xzzzz[i] * pa_x[i];

        tr_y_xxyy_yyyyy[i] = tr_y_yy_yyyyy[i] * fe_0 + tr_y_xyy_yyyyy[i] * pa_x[i];

        tr_y_xxyy_yyyyz[i] = tr_y_yy_yyyyz[i] * fe_0 + tr_y_xyy_yyyyz[i] * pa_x[i];

        tr_y_xxyy_yyyzz[i] = tr_y_yy_yyyzz[i] * fe_0 + tr_y_xyy_yyyzz[i] * pa_x[i];

        tr_y_xxyy_yyzzz[i] = tr_y_yy_yyzzz[i] * fe_0 + tr_y_xyy_yyzzz[i] * pa_x[i];

        tr_y_xxyy_yzzzz[i] = tr_y_yy_yzzzz[i] * fe_0 + tr_y_xyy_yzzzz[i] * pa_x[i];

        tr_y_xxyy_zzzzz[i] = tr_y_yy_zzzzz[i] * fe_0 + tr_y_xyy_zzzzz[i] * pa_x[i];
    }

    // Set up 399-420 components of targeted buffer : GH

    auto tr_y_xxyz_xxxxx = pbuffer.data(idx_dip_gh + 399);

    auto tr_y_xxyz_xxxxy = pbuffer.data(idx_dip_gh + 400);

    auto tr_y_xxyz_xxxxz = pbuffer.data(idx_dip_gh + 401);

    auto tr_y_xxyz_xxxyy = pbuffer.data(idx_dip_gh + 402);

    auto tr_y_xxyz_xxxyz = pbuffer.data(idx_dip_gh + 403);

    auto tr_y_xxyz_xxxzz = pbuffer.data(idx_dip_gh + 404);

    auto tr_y_xxyz_xxyyy = pbuffer.data(idx_dip_gh + 405);

    auto tr_y_xxyz_xxyyz = pbuffer.data(idx_dip_gh + 406);

    auto tr_y_xxyz_xxyzz = pbuffer.data(idx_dip_gh + 407);

    auto tr_y_xxyz_xxzzz = pbuffer.data(idx_dip_gh + 408);

    auto tr_y_xxyz_xyyyy = pbuffer.data(idx_dip_gh + 409);

    auto tr_y_xxyz_xyyyz = pbuffer.data(idx_dip_gh + 410);

    auto tr_y_xxyz_xyyzz = pbuffer.data(idx_dip_gh + 411);

    auto tr_y_xxyz_xyzzz = pbuffer.data(idx_dip_gh + 412);

    auto tr_y_xxyz_xzzzz = pbuffer.data(idx_dip_gh + 413);

    auto tr_y_xxyz_yyyyy = pbuffer.data(idx_dip_gh + 414);

    auto tr_y_xxyz_yyyyz = pbuffer.data(idx_dip_gh + 415);

    auto tr_y_xxyz_yyyzz = pbuffer.data(idx_dip_gh + 416);

    auto tr_y_xxyz_yyzzz = pbuffer.data(idx_dip_gh + 417);

    auto tr_y_xxyz_yzzzz = pbuffer.data(idx_dip_gh + 418);

    auto tr_y_xxyz_zzzzz = pbuffer.data(idx_dip_gh + 419);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxy_xxxxx, tr_y_xxy_xxxxy, tr_y_xxy_xxxy, tr_y_xxy_xxxyy, tr_y_xxy_xxxyz, tr_y_xxy_xxyy, tr_y_xxy_xxyyy, tr_y_xxy_xxyyz, tr_y_xxy_xxyz, tr_y_xxy_xxyzz, tr_y_xxy_xyyy, tr_y_xxy_xyyyy, tr_y_xxy_xyyyz, tr_y_xxy_xyyz, tr_y_xxy_xyyzz, tr_y_xxy_xyzz, tr_y_xxy_xyzzz, tr_y_xxy_yyyyy, tr_y_xxyz_xxxxx, tr_y_xxyz_xxxxy, tr_y_xxyz_xxxxz, tr_y_xxyz_xxxyy, tr_y_xxyz_xxxyz, tr_y_xxyz_xxxzz, tr_y_xxyz_xxyyy, tr_y_xxyz_xxyyz, tr_y_xxyz_xxyzz, tr_y_xxyz_xxzzz, tr_y_xxyz_xyyyy, tr_y_xxyz_xyyyz, tr_y_xxyz_xyyzz, tr_y_xxyz_xyzzz, tr_y_xxyz_xzzzz, tr_y_xxyz_yyyyy, tr_y_xxyz_yyyyz, tr_y_xxyz_yyyzz, tr_y_xxyz_yyzzz, tr_y_xxyz_yzzzz, tr_y_xxyz_zzzzz, tr_y_xxz_xxxxz, tr_y_xxz_xxxzz, tr_y_xxz_xxzzz, tr_y_xxz_xzzzz, tr_y_xyz_yyyyz, tr_y_xyz_yyyzz, tr_y_xyz_yyzzz, tr_y_xyz_yzzzz, tr_y_xyz_zzzzz, tr_y_yz_yyyyz, tr_y_yz_yyyzz, tr_y_yz_yyzzz, tr_y_yz_yzzzz, tr_y_yz_zzzzz, ts_xxz_xxxxz, ts_xxz_xxxzz, ts_xxz_xxzzz, ts_xxz_xzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyz_xxxxx[i] = tr_y_xxy_xxxxx[i] * pa_z[i];

        tr_y_xxyz_xxxxy[i] = tr_y_xxy_xxxxy[i] * pa_z[i];

        tr_y_xxyz_xxxxz[i] = ts_xxz_xxxxz[i] * fe_0 + tr_y_xxz_xxxxz[i] * pa_y[i];

        tr_y_xxyz_xxxyy[i] = tr_y_xxy_xxxyy[i] * pa_z[i];

        tr_y_xxyz_xxxyz[i] = tr_y_xxy_xxxy[i] * fe_0 + tr_y_xxy_xxxyz[i] * pa_z[i];

        tr_y_xxyz_xxxzz[i] = ts_xxz_xxxzz[i] * fe_0 + tr_y_xxz_xxxzz[i] * pa_y[i];

        tr_y_xxyz_xxyyy[i] = tr_y_xxy_xxyyy[i] * pa_z[i];

        tr_y_xxyz_xxyyz[i] = tr_y_xxy_xxyy[i] * fe_0 + tr_y_xxy_xxyyz[i] * pa_z[i];

        tr_y_xxyz_xxyzz[i] = 2.0 * tr_y_xxy_xxyz[i] * fe_0 + tr_y_xxy_xxyzz[i] * pa_z[i];

        tr_y_xxyz_xxzzz[i] = ts_xxz_xxzzz[i] * fe_0 + tr_y_xxz_xxzzz[i] * pa_y[i];

        tr_y_xxyz_xyyyy[i] = tr_y_xxy_xyyyy[i] * pa_z[i];

        tr_y_xxyz_xyyyz[i] = tr_y_xxy_xyyy[i] * fe_0 + tr_y_xxy_xyyyz[i] * pa_z[i];

        tr_y_xxyz_xyyzz[i] = 2.0 * tr_y_xxy_xyyz[i] * fe_0 + tr_y_xxy_xyyzz[i] * pa_z[i];

        tr_y_xxyz_xyzzz[i] = 3.0 * tr_y_xxy_xyzz[i] * fe_0 + tr_y_xxy_xyzzz[i] * pa_z[i];

        tr_y_xxyz_xzzzz[i] = ts_xxz_xzzzz[i] * fe_0 + tr_y_xxz_xzzzz[i] * pa_y[i];

        tr_y_xxyz_yyyyy[i] = tr_y_xxy_yyyyy[i] * pa_z[i];

        tr_y_xxyz_yyyyz[i] = tr_y_yz_yyyyz[i] * fe_0 + tr_y_xyz_yyyyz[i] * pa_x[i];

        tr_y_xxyz_yyyzz[i] = tr_y_yz_yyyzz[i] * fe_0 + tr_y_xyz_yyyzz[i] * pa_x[i];

        tr_y_xxyz_yyzzz[i] = tr_y_yz_yyzzz[i] * fe_0 + tr_y_xyz_yyzzz[i] * pa_x[i];

        tr_y_xxyz_yzzzz[i] = tr_y_yz_yzzzz[i] * fe_0 + tr_y_xyz_yzzzz[i] * pa_x[i];

        tr_y_xxyz_zzzzz[i] = tr_y_yz_zzzzz[i] * fe_0 + tr_y_xyz_zzzzz[i] * pa_x[i];
    }

    // Set up 420-441 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xx_xxxxx, tr_y_xx_xxxxy, tr_y_xx_xxxyy, tr_y_xx_xxyyy, tr_y_xx_xyyyy, tr_y_xxz_xxxxx, tr_y_xxz_xxxxy, tr_y_xxz_xxxyy, tr_y_xxz_xxyyy, tr_y_xxz_xyyyy, tr_y_xxzz_xxxxx, tr_y_xxzz_xxxxy, tr_y_xxzz_xxxxz, tr_y_xxzz_xxxyy, tr_y_xxzz_xxxyz, tr_y_xxzz_xxxzz, tr_y_xxzz_xxyyy, tr_y_xxzz_xxyyz, tr_y_xxzz_xxyzz, tr_y_xxzz_xxzzz, tr_y_xxzz_xyyyy, tr_y_xxzz_xyyyz, tr_y_xxzz_xyyzz, tr_y_xxzz_xyzzz, tr_y_xxzz_xzzzz, tr_y_xxzz_yyyyy, tr_y_xxzz_yyyyz, tr_y_xxzz_yyyzz, tr_y_xxzz_yyzzz, tr_y_xxzz_yzzzz, tr_y_xxzz_zzzzz, tr_y_xzz_xxxxz, tr_y_xzz_xxxyz, tr_y_xzz_xxxz, tr_y_xzz_xxxzz, tr_y_xzz_xxyyz, tr_y_xzz_xxyz, tr_y_xzz_xxyzz, tr_y_xzz_xxzz, tr_y_xzz_xxzzz, tr_y_xzz_xyyyz, tr_y_xzz_xyyz, tr_y_xzz_xyyzz, tr_y_xzz_xyzz, tr_y_xzz_xyzzz, tr_y_xzz_xzzz, tr_y_xzz_xzzzz, tr_y_xzz_yyyyy, tr_y_xzz_yyyyz, tr_y_xzz_yyyz, tr_y_xzz_yyyzz, tr_y_xzz_yyzz, tr_y_xzz_yyzzz, tr_y_xzz_yzzz, tr_y_xzz_yzzzz, tr_y_xzz_zzzz, tr_y_xzz_zzzzz, tr_y_zz_xxxxz, tr_y_zz_xxxyz, tr_y_zz_xxxzz, tr_y_zz_xxyyz, tr_y_zz_xxyzz, tr_y_zz_xxzzz, tr_y_zz_xyyyz, tr_y_zz_xyyzz, tr_y_zz_xyzzz, tr_y_zz_xzzzz, tr_y_zz_yyyyy, tr_y_zz_yyyyz, tr_y_zz_yyyzz, tr_y_zz_yyzzz, tr_y_zz_yzzzz, tr_y_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzz_xxxxx[i] = tr_y_xx_xxxxx[i] * fe_0 + tr_y_xxz_xxxxx[i] * pa_z[i];

        tr_y_xxzz_xxxxy[i] = tr_y_xx_xxxxy[i] * fe_0 + tr_y_xxz_xxxxy[i] * pa_z[i];

        tr_y_xxzz_xxxxz[i] = tr_y_zz_xxxxz[i] * fe_0 + 4.0 * tr_y_xzz_xxxz[i] * fe_0 + tr_y_xzz_xxxxz[i] * pa_x[i];

        tr_y_xxzz_xxxyy[i] = tr_y_xx_xxxyy[i] * fe_0 + tr_y_xxz_xxxyy[i] * pa_z[i];

        tr_y_xxzz_xxxyz[i] = tr_y_zz_xxxyz[i] * fe_0 + 3.0 * tr_y_xzz_xxyz[i] * fe_0 + tr_y_xzz_xxxyz[i] * pa_x[i];

        tr_y_xxzz_xxxzz[i] = tr_y_zz_xxxzz[i] * fe_0 + 3.0 * tr_y_xzz_xxzz[i] * fe_0 + tr_y_xzz_xxxzz[i] * pa_x[i];

        tr_y_xxzz_xxyyy[i] = tr_y_xx_xxyyy[i] * fe_0 + tr_y_xxz_xxyyy[i] * pa_z[i];

        tr_y_xxzz_xxyyz[i] = tr_y_zz_xxyyz[i] * fe_0 + 2.0 * tr_y_xzz_xyyz[i] * fe_0 + tr_y_xzz_xxyyz[i] * pa_x[i];

        tr_y_xxzz_xxyzz[i] = tr_y_zz_xxyzz[i] * fe_0 + 2.0 * tr_y_xzz_xyzz[i] * fe_0 + tr_y_xzz_xxyzz[i] * pa_x[i];

        tr_y_xxzz_xxzzz[i] = tr_y_zz_xxzzz[i] * fe_0 + 2.0 * tr_y_xzz_xzzz[i] * fe_0 + tr_y_xzz_xxzzz[i] * pa_x[i];

        tr_y_xxzz_xyyyy[i] = tr_y_xx_xyyyy[i] * fe_0 + tr_y_xxz_xyyyy[i] * pa_z[i];

        tr_y_xxzz_xyyyz[i] = tr_y_zz_xyyyz[i] * fe_0 + tr_y_xzz_yyyz[i] * fe_0 + tr_y_xzz_xyyyz[i] * pa_x[i];

        tr_y_xxzz_xyyzz[i] = tr_y_zz_xyyzz[i] * fe_0 + tr_y_xzz_yyzz[i] * fe_0 + tr_y_xzz_xyyzz[i] * pa_x[i];

        tr_y_xxzz_xyzzz[i] = tr_y_zz_xyzzz[i] * fe_0 + tr_y_xzz_yzzz[i] * fe_0 + tr_y_xzz_xyzzz[i] * pa_x[i];

        tr_y_xxzz_xzzzz[i] = tr_y_zz_xzzzz[i] * fe_0 + tr_y_xzz_zzzz[i] * fe_0 + tr_y_xzz_xzzzz[i] * pa_x[i];

        tr_y_xxzz_yyyyy[i] = tr_y_zz_yyyyy[i] * fe_0 + tr_y_xzz_yyyyy[i] * pa_x[i];

        tr_y_xxzz_yyyyz[i] = tr_y_zz_yyyyz[i] * fe_0 + tr_y_xzz_yyyyz[i] * pa_x[i];

        tr_y_xxzz_yyyzz[i] = tr_y_zz_yyyzz[i] * fe_0 + tr_y_xzz_yyyzz[i] * pa_x[i];

        tr_y_xxzz_yyzzz[i] = tr_y_zz_yyzzz[i] * fe_0 + tr_y_xzz_yyzzz[i] * pa_x[i];

        tr_y_xxzz_yzzzz[i] = tr_y_zz_yzzzz[i] * fe_0 + tr_y_xzz_yzzzz[i] * pa_x[i];

        tr_y_xxzz_zzzzz[i] = tr_y_zz_zzzzz[i] * fe_0 + tr_y_xzz_zzzzz[i] * pa_x[i];
    }

    // Set up 441-462 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_y_xyyy_xxxxx, tr_y_xyyy_xxxxy, tr_y_xyyy_xxxxz, tr_y_xyyy_xxxyy, tr_y_xyyy_xxxyz, tr_y_xyyy_xxxzz, tr_y_xyyy_xxyyy, tr_y_xyyy_xxyyz, tr_y_xyyy_xxyzz, tr_y_xyyy_xxzzz, tr_y_xyyy_xyyyy, tr_y_xyyy_xyyyz, tr_y_xyyy_xyyzz, tr_y_xyyy_xyzzz, tr_y_xyyy_xzzzz, tr_y_xyyy_yyyyy, tr_y_xyyy_yyyyz, tr_y_xyyy_yyyzz, tr_y_xyyy_yyzzz, tr_y_xyyy_yzzzz, tr_y_xyyy_zzzzz, tr_y_yyy_xxxx, tr_y_yyy_xxxxx, tr_y_yyy_xxxxy, tr_y_yyy_xxxxz, tr_y_yyy_xxxy, tr_y_yyy_xxxyy, tr_y_yyy_xxxyz, tr_y_yyy_xxxz, tr_y_yyy_xxxzz, tr_y_yyy_xxyy, tr_y_yyy_xxyyy, tr_y_yyy_xxyyz, tr_y_yyy_xxyz, tr_y_yyy_xxyzz, tr_y_yyy_xxzz, tr_y_yyy_xxzzz, tr_y_yyy_xyyy, tr_y_yyy_xyyyy, tr_y_yyy_xyyyz, tr_y_yyy_xyyz, tr_y_yyy_xyyzz, tr_y_yyy_xyzz, tr_y_yyy_xyzzz, tr_y_yyy_xzzz, tr_y_yyy_xzzzz, tr_y_yyy_yyyy, tr_y_yyy_yyyyy, tr_y_yyy_yyyyz, tr_y_yyy_yyyz, tr_y_yyy_yyyzz, tr_y_yyy_yyzz, tr_y_yyy_yyzzz, tr_y_yyy_yzzz, tr_y_yyy_yzzzz, tr_y_yyy_zzzz, tr_y_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyy_xxxxx[i] = 5.0 * tr_y_yyy_xxxx[i] * fe_0 + tr_y_yyy_xxxxx[i] * pa_x[i];

        tr_y_xyyy_xxxxy[i] = 4.0 * tr_y_yyy_xxxy[i] * fe_0 + tr_y_yyy_xxxxy[i] * pa_x[i];

        tr_y_xyyy_xxxxz[i] = 4.0 * tr_y_yyy_xxxz[i] * fe_0 + tr_y_yyy_xxxxz[i] * pa_x[i];

        tr_y_xyyy_xxxyy[i] = 3.0 * tr_y_yyy_xxyy[i] * fe_0 + tr_y_yyy_xxxyy[i] * pa_x[i];

        tr_y_xyyy_xxxyz[i] = 3.0 * tr_y_yyy_xxyz[i] * fe_0 + tr_y_yyy_xxxyz[i] * pa_x[i];

        tr_y_xyyy_xxxzz[i] = 3.0 * tr_y_yyy_xxzz[i] * fe_0 + tr_y_yyy_xxxzz[i] * pa_x[i];

        tr_y_xyyy_xxyyy[i] = 2.0 * tr_y_yyy_xyyy[i] * fe_0 + tr_y_yyy_xxyyy[i] * pa_x[i];

        tr_y_xyyy_xxyyz[i] = 2.0 * tr_y_yyy_xyyz[i] * fe_0 + tr_y_yyy_xxyyz[i] * pa_x[i];

        tr_y_xyyy_xxyzz[i] = 2.0 * tr_y_yyy_xyzz[i] * fe_0 + tr_y_yyy_xxyzz[i] * pa_x[i];

        tr_y_xyyy_xxzzz[i] = 2.0 * tr_y_yyy_xzzz[i] * fe_0 + tr_y_yyy_xxzzz[i] * pa_x[i];

        tr_y_xyyy_xyyyy[i] = tr_y_yyy_yyyy[i] * fe_0 + tr_y_yyy_xyyyy[i] * pa_x[i];

        tr_y_xyyy_xyyyz[i] = tr_y_yyy_yyyz[i] * fe_0 + tr_y_yyy_xyyyz[i] * pa_x[i];

        tr_y_xyyy_xyyzz[i] = tr_y_yyy_yyzz[i] * fe_0 + tr_y_yyy_xyyzz[i] * pa_x[i];

        tr_y_xyyy_xyzzz[i] = tr_y_yyy_yzzz[i] * fe_0 + tr_y_yyy_xyzzz[i] * pa_x[i];

        tr_y_xyyy_xzzzz[i] = tr_y_yyy_zzzz[i] * fe_0 + tr_y_yyy_xzzzz[i] * pa_x[i];

        tr_y_xyyy_yyyyy[i] = tr_y_yyy_yyyyy[i] * pa_x[i];

        tr_y_xyyy_yyyyz[i] = tr_y_yyy_yyyyz[i] * pa_x[i];

        tr_y_xyyy_yyyzz[i] = tr_y_yyy_yyyzz[i] * pa_x[i];

        tr_y_xyyy_yyzzz[i] = tr_y_yyy_yyzzz[i] * pa_x[i];

        tr_y_xyyy_yzzzz[i] = tr_y_yyy_yzzzz[i] * pa_x[i];

        tr_y_xyyy_zzzzz[i] = tr_y_yyy_zzzzz[i] * pa_x[i];
    }

    // Set up 462-483 components of targeted buffer : GH

    auto tr_y_xyyz_xxxxx = pbuffer.data(idx_dip_gh + 462);

    auto tr_y_xyyz_xxxxy = pbuffer.data(idx_dip_gh + 463);

    auto tr_y_xyyz_xxxxz = pbuffer.data(idx_dip_gh + 464);

    auto tr_y_xyyz_xxxyy = pbuffer.data(idx_dip_gh + 465);

    auto tr_y_xyyz_xxxyz = pbuffer.data(idx_dip_gh + 466);

    auto tr_y_xyyz_xxxzz = pbuffer.data(idx_dip_gh + 467);

    auto tr_y_xyyz_xxyyy = pbuffer.data(idx_dip_gh + 468);

    auto tr_y_xyyz_xxyyz = pbuffer.data(idx_dip_gh + 469);

    auto tr_y_xyyz_xxyzz = pbuffer.data(idx_dip_gh + 470);

    auto tr_y_xyyz_xxzzz = pbuffer.data(idx_dip_gh + 471);

    auto tr_y_xyyz_xyyyy = pbuffer.data(idx_dip_gh + 472);

    auto tr_y_xyyz_xyyyz = pbuffer.data(idx_dip_gh + 473);

    auto tr_y_xyyz_xyyzz = pbuffer.data(idx_dip_gh + 474);

    auto tr_y_xyyz_xyzzz = pbuffer.data(idx_dip_gh + 475);

    auto tr_y_xyyz_xzzzz = pbuffer.data(idx_dip_gh + 476);

    auto tr_y_xyyz_yyyyy = pbuffer.data(idx_dip_gh + 477);

    auto tr_y_xyyz_yyyyz = pbuffer.data(idx_dip_gh + 478);

    auto tr_y_xyyz_yyyzz = pbuffer.data(idx_dip_gh + 479);

    auto tr_y_xyyz_yyzzz = pbuffer.data(idx_dip_gh + 480);

    auto tr_y_xyyz_yzzzz = pbuffer.data(idx_dip_gh + 481);

    auto tr_y_xyyz_zzzzz = pbuffer.data(idx_dip_gh + 482);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xyy_xxxxx, tr_y_xyy_xxxxy, tr_y_xyy_xxxyy, tr_y_xyy_xxyyy, tr_y_xyy_xyyyy, tr_y_xyyz_xxxxx, tr_y_xyyz_xxxxy, tr_y_xyyz_xxxxz, tr_y_xyyz_xxxyy, tr_y_xyyz_xxxyz, tr_y_xyyz_xxxzz, tr_y_xyyz_xxyyy, tr_y_xyyz_xxyyz, tr_y_xyyz_xxyzz, tr_y_xyyz_xxzzz, tr_y_xyyz_xyyyy, tr_y_xyyz_xyyyz, tr_y_xyyz_xyyzz, tr_y_xyyz_xyzzz, tr_y_xyyz_xzzzz, tr_y_xyyz_yyyyy, tr_y_xyyz_yyyyz, tr_y_xyyz_yyyzz, tr_y_xyyz_yyzzz, tr_y_xyyz_yzzzz, tr_y_xyyz_zzzzz, tr_y_yyz_xxxxz, tr_y_yyz_xxxyz, tr_y_yyz_xxxz, tr_y_yyz_xxxzz, tr_y_yyz_xxyyz, tr_y_yyz_xxyz, tr_y_yyz_xxyzz, tr_y_yyz_xxzz, tr_y_yyz_xxzzz, tr_y_yyz_xyyyz, tr_y_yyz_xyyz, tr_y_yyz_xyyzz, tr_y_yyz_xyzz, tr_y_yyz_xyzzz, tr_y_yyz_xzzz, tr_y_yyz_xzzzz, tr_y_yyz_yyyyy, tr_y_yyz_yyyyz, tr_y_yyz_yyyz, tr_y_yyz_yyyzz, tr_y_yyz_yyzz, tr_y_yyz_yyzzz, tr_y_yyz_yzzz, tr_y_yyz_yzzzz, tr_y_yyz_zzzz, tr_y_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyz_xxxxx[i] = tr_y_xyy_xxxxx[i] * pa_z[i];

        tr_y_xyyz_xxxxy[i] = tr_y_xyy_xxxxy[i] * pa_z[i];

        tr_y_xyyz_xxxxz[i] = 4.0 * tr_y_yyz_xxxz[i] * fe_0 + tr_y_yyz_xxxxz[i] * pa_x[i];

        tr_y_xyyz_xxxyy[i] = tr_y_xyy_xxxyy[i] * pa_z[i];

        tr_y_xyyz_xxxyz[i] = 3.0 * tr_y_yyz_xxyz[i] * fe_0 + tr_y_yyz_xxxyz[i] * pa_x[i];

        tr_y_xyyz_xxxzz[i] = 3.0 * tr_y_yyz_xxzz[i] * fe_0 + tr_y_yyz_xxxzz[i] * pa_x[i];

        tr_y_xyyz_xxyyy[i] = tr_y_xyy_xxyyy[i] * pa_z[i];

        tr_y_xyyz_xxyyz[i] = 2.0 * tr_y_yyz_xyyz[i] * fe_0 + tr_y_yyz_xxyyz[i] * pa_x[i];

        tr_y_xyyz_xxyzz[i] = 2.0 * tr_y_yyz_xyzz[i] * fe_0 + tr_y_yyz_xxyzz[i] * pa_x[i];

        tr_y_xyyz_xxzzz[i] = 2.0 * tr_y_yyz_xzzz[i] * fe_0 + tr_y_yyz_xxzzz[i] * pa_x[i];

        tr_y_xyyz_xyyyy[i] = tr_y_xyy_xyyyy[i] * pa_z[i];

        tr_y_xyyz_xyyyz[i] = tr_y_yyz_yyyz[i] * fe_0 + tr_y_yyz_xyyyz[i] * pa_x[i];

        tr_y_xyyz_xyyzz[i] = tr_y_yyz_yyzz[i] * fe_0 + tr_y_yyz_xyyzz[i] * pa_x[i];

        tr_y_xyyz_xyzzz[i] = tr_y_yyz_yzzz[i] * fe_0 + tr_y_yyz_xyzzz[i] * pa_x[i];

        tr_y_xyyz_xzzzz[i] = tr_y_yyz_zzzz[i] * fe_0 + tr_y_yyz_xzzzz[i] * pa_x[i];

        tr_y_xyyz_yyyyy[i] = tr_y_yyz_yyyyy[i] * pa_x[i];

        tr_y_xyyz_yyyyz[i] = tr_y_yyz_yyyyz[i] * pa_x[i];

        tr_y_xyyz_yyyzz[i] = tr_y_yyz_yyyzz[i] * pa_x[i];

        tr_y_xyyz_yyzzz[i] = tr_y_yyz_yyzzz[i] * pa_x[i];

        tr_y_xyyz_yzzzz[i] = tr_y_yyz_yzzzz[i] * pa_x[i];

        tr_y_xyyz_zzzzz[i] = tr_y_yyz_zzzzz[i] * pa_x[i];
    }

    // Set up 483-504 components of targeted buffer : GH

    auto tr_y_xyzz_xxxxx = pbuffer.data(idx_dip_gh + 483);

    auto tr_y_xyzz_xxxxy = pbuffer.data(idx_dip_gh + 484);

    auto tr_y_xyzz_xxxxz = pbuffer.data(idx_dip_gh + 485);

    auto tr_y_xyzz_xxxyy = pbuffer.data(idx_dip_gh + 486);

    auto tr_y_xyzz_xxxyz = pbuffer.data(idx_dip_gh + 487);

    auto tr_y_xyzz_xxxzz = pbuffer.data(idx_dip_gh + 488);

    auto tr_y_xyzz_xxyyy = pbuffer.data(idx_dip_gh + 489);

    auto tr_y_xyzz_xxyyz = pbuffer.data(idx_dip_gh + 490);

    auto tr_y_xyzz_xxyzz = pbuffer.data(idx_dip_gh + 491);

    auto tr_y_xyzz_xxzzz = pbuffer.data(idx_dip_gh + 492);

    auto tr_y_xyzz_xyyyy = pbuffer.data(idx_dip_gh + 493);

    auto tr_y_xyzz_xyyyz = pbuffer.data(idx_dip_gh + 494);

    auto tr_y_xyzz_xyyzz = pbuffer.data(idx_dip_gh + 495);

    auto tr_y_xyzz_xyzzz = pbuffer.data(idx_dip_gh + 496);

    auto tr_y_xyzz_xzzzz = pbuffer.data(idx_dip_gh + 497);

    auto tr_y_xyzz_yyyyy = pbuffer.data(idx_dip_gh + 498);

    auto tr_y_xyzz_yyyyz = pbuffer.data(idx_dip_gh + 499);

    auto tr_y_xyzz_yyyzz = pbuffer.data(idx_dip_gh + 500);

    auto tr_y_xyzz_yyzzz = pbuffer.data(idx_dip_gh + 501);

    auto tr_y_xyzz_yzzzz = pbuffer.data(idx_dip_gh + 502);

    auto tr_y_xyzz_zzzzz = pbuffer.data(idx_dip_gh + 503);

    #pragma omp simd aligned(pa_x, tr_y_xyzz_xxxxx, tr_y_xyzz_xxxxy, tr_y_xyzz_xxxxz, tr_y_xyzz_xxxyy, tr_y_xyzz_xxxyz, tr_y_xyzz_xxxzz, tr_y_xyzz_xxyyy, tr_y_xyzz_xxyyz, tr_y_xyzz_xxyzz, tr_y_xyzz_xxzzz, tr_y_xyzz_xyyyy, tr_y_xyzz_xyyyz, tr_y_xyzz_xyyzz, tr_y_xyzz_xyzzz, tr_y_xyzz_xzzzz, tr_y_xyzz_yyyyy, tr_y_xyzz_yyyyz, tr_y_xyzz_yyyzz, tr_y_xyzz_yyzzz, tr_y_xyzz_yzzzz, tr_y_xyzz_zzzzz, tr_y_yzz_xxxx, tr_y_yzz_xxxxx, tr_y_yzz_xxxxy, tr_y_yzz_xxxxz, tr_y_yzz_xxxy, tr_y_yzz_xxxyy, tr_y_yzz_xxxyz, tr_y_yzz_xxxz, tr_y_yzz_xxxzz, tr_y_yzz_xxyy, tr_y_yzz_xxyyy, tr_y_yzz_xxyyz, tr_y_yzz_xxyz, tr_y_yzz_xxyzz, tr_y_yzz_xxzz, tr_y_yzz_xxzzz, tr_y_yzz_xyyy, tr_y_yzz_xyyyy, tr_y_yzz_xyyyz, tr_y_yzz_xyyz, tr_y_yzz_xyyzz, tr_y_yzz_xyzz, tr_y_yzz_xyzzz, tr_y_yzz_xzzz, tr_y_yzz_xzzzz, tr_y_yzz_yyyy, tr_y_yzz_yyyyy, tr_y_yzz_yyyyz, tr_y_yzz_yyyz, tr_y_yzz_yyyzz, tr_y_yzz_yyzz, tr_y_yzz_yyzzz, tr_y_yzz_yzzz, tr_y_yzz_yzzzz, tr_y_yzz_zzzz, tr_y_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzz_xxxxx[i] = 5.0 * tr_y_yzz_xxxx[i] * fe_0 + tr_y_yzz_xxxxx[i] * pa_x[i];

        tr_y_xyzz_xxxxy[i] = 4.0 * tr_y_yzz_xxxy[i] * fe_0 + tr_y_yzz_xxxxy[i] * pa_x[i];

        tr_y_xyzz_xxxxz[i] = 4.0 * tr_y_yzz_xxxz[i] * fe_0 + tr_y_yzz_xxxxz[i] * pa_x[i];

        tr_y_xyzz_xxxyy[i] = 3.0 * tr_y_yzz_xxyy[i] * fe_0 + tr_y_yzz_xxxyy[i] * pa_x[i];

        tr_y_xyzz_xxxyz[i] = 3.0 * tr_y_yzz_xxyz[i] * fe_0 + tr_y_yzz_xxxyz[i] * pa_x[i];

        tr_y_xyzz_xxxzz[i] = 3.0 * tr_y_yzz_xxzz[i] * fe_0 + tr_y_yzz_xxxzz[i] * pa_x[i];

        tr_y_xyzz_xxyyy[i] = 2.0 * tr_y_yzz_xyyy[i] * fe_0 + tr_y_yzz_xxyyy[i] * pa_x[i];

        tr_y_xyzz_xxyyz[i] = 2.0 * tr_y_yzz_xyyz[i] * fe_0 + tr_y_yzz_xxyyz[i] * pa_x[i];

        tr_y_xyzz_xxyzz[i] = 2.0 * tr_y_yzz_xyzz[i] * fe_0 + tr_y_yzz_xxyzz[i] * pa_x[i];

        tr_y_xyzz_xxzzz[i] = 2.0 * tr_y_yzz_xzzz[i] * fe_0 + tr_y_yzz_xxzzz[i] * pa_x[i];

        tr_y_xyzz_xyyyy[i] = tr_y_yzz_yyyy[i] * fe_0 + tr_y_yzz_xyyyy[i] * pa_x[i];

        tr_y_xyzz_xyyyz[i] = tr_y_yzz_yyyz[i] * fe_0 + tr_y_yzz_xyyyz[i] * pa_x[i];

        tr_y_xyzz_xyyzz[i] = tr_y_yzz_yyzz[i] * fe_0 + tr_y_yzz_xyyzz[i] * pa_x[i];

        tr_y_xyzz_xyzzz[i] = tr_y_yzz_yzzz[i] * fe_0 + tr_y_yzz_xyzzz[i] * pa_x[i];

        tr_y_xyzz_xzzzz[i] = tr_y_yzz_zzzz[i] * fe_0 + tr_y_yzz_xzzzz[i] * pa_x[i];

        tr_y_xyzz_yyyyy[i] = tr_y_yzz_yyyyy[i] * pa_x[i];

        tr_y_xyzz_yyyyz[i] = tr_y_yzz_yyyyz[i] * pa_x[i];

        tr_y_xyzz_yyyzz[i] = tr_y_yzz_yyyzz[i] * pa_x[i];

        tr_y_xyzz_yyzzz[i] = tr_y_yzz_yyzzz[i] * pa_x[i];

        tr_y_xyzz_yzzzz[i] = tr_y_yzz_yzzzz[i] * pa_x[i];

        tr_y_xyzz_zzzzz[i] = tr_y_yzz_zzzzz[i] * pa_x[i];
    }

    // Set up 504-525 components of targeted buffer : GH

    auto tr_y_xzzz_xxxxx = pbuffer.data(idx_dip_gh + 504);

    auto tr_y_xzzz_xxxxy = pbuffer.data(idx_dip_gh + 505);

    auto tr_y_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 506);

    auto tr_y_xzzz_xxxyy = pbuffer.data(idx_dip_gh + 507);

    auto tr_y_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 508);

    auto tr_y_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 509);

    auto tr_y_xzzz_xxyyy = pbuffer.data(idx_dip_gh + 510);

    auto tr_y_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 511);

    auto tr_y_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 512);

    auto tr_y_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 513);

    auto tr_y_xzzz_xyyyy = pbuffer.data(idx_dip_gh + 514);

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

    #pragma omp simd aligned(pa_x, tr_y_xzzz_xxxxx, tr_y_xzzz_xxxxy, tr_y_xzzz_xxxxz, tr_y_xzzz_xxxyy, tr_y_xzzz_xxxyz, tr_y_xzzz_xxxzz, tr_y_xzzz_xxyyy, tr_y_xzzz_xxyyz, tr_y_xzzz_xxyzz, tr_y_xzzz_xxzzz, tr_y_xzzz_xyyyy, tr_y_xzzz_xyyyz, tr_y_xzzz_xyyzz, tr_y_xzzz_xyzzz, tr_y_xzzz_xzzzz, tr_y_xzzz_yyyyy, tr_y_xzzz_yyyyz, tr_y_xzzz_yyyzz, tr_y_xzzz_yyzzz, tr_y_xzzz_yzzzz, tr_y_xzzz_zzzzz, tr_y_zzz_xxxx, tr_y_zzz_xxxxx, tr_y_zzz_xxxxy, tr_y_zzz_xxxxz, tr_y_zzz_xxxy, tr_y_zzz_xxxyy, tr_y_zzz_xxxyz, tr_y_zzz_xxxz, tr_y_zzz_xxxzz, tr_y_zzz_xxyy, tr_y_zzz_xxyyy, tr_y_zzz_xxyyz, tr_y_zzz_xxyz, tr_y_zzz_xxyzz, tr_y_zzz_xxzz, tr_y_zzz_xxzzz, tr_y_zzz_xyyy, tr_y_zzz_xyyyy, tr_y_zzz_xyyyz, tr_y_zzz_xyyz, tr_y_zzz_xyyzz, tr_y_zzz_xyzz, tr_y_zzz_xyzzz, tr_y_zzz_xzzz, tr_y_zzz_xzzzz, tr_y_zzz_yyyy, tr_y_zzz_yyyyy, tr_y_zzz_yyyyz, tr_y_zzz_yyyz, tr_y_zzz_yyyzz, tr_y_zzz_yyzz, tr_y_zzz_yyzzz, tr_y_zzz_yzzz, tr_y_zzz_yzzzz, tr_y_zzz_zzzz, tr_y_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzz_xxxxx[i] = 5.0 * tr_y_zzz_xxxx[i] * fe_0 + tr_y_zzz_xxxxx[i] * pa_x[i];

        tr_y_xzzz_xxxxy[i] = 4.0 * tr_y_zzz_xxxy[i] * fe_0 + tr_y_zzz_xxxxy[i] * pa_x[i];

        tr_y_xzzz_xxxxz[i] = 4.0 * tr_y_zzz_xxxz[i] * fe_0 + tr_y_zzz_xxxxz[i] * pa_x[i];

        tr_y_xzzz_xxxyy[i] = 3.0 * tr_y_zzz_xxyy[i] * fe_0 + tr_y_zzz_xxxyy[i] * pa_x[i];

        tr_y_xzzz_xxxyz[i] = 3.0 * tr_y_zzz_xxyz[i] * fe_0 + tr_y_zzz_xxxyz[i] * pa_x[i];

        tr_y_xzzz_xxxzz[i] = 3.0 * tr_y_zzz_xxzz[i] * fe_0 + tr_y_zzz_xxxzz[i] * pa_x[i];

        tr_y_xzzz_xxyyy[i] = 2.0 * tr_y_zzz_xyyy[i] * fe_0 + tr_y_zzz_xxyyy[i] * pa_x[i];

        tr_y_xzzz_xxyyz[i] = 2.0 * tr_y_zzz_xyyz[i] * fe_0 + tr_y_zzz_xxyyz[i] * pa_x[i];

        tr_y_xzzz_xxyzz[i] = 2.0 * tr_y_zzz_xyzz[i] * fe_0 + tr_y_zzz_xxyzz[i] * pa_x[i];

        tr_y_xzzz_xxzzz[i] = 2.0 * tr_y_zzz_xzzz[i] * fe_0 + tr_y_zzz_xxzzz[i] * pa_x[i];

        tr_y_xzzz_xyyyy[i] = tr_y_zzz_yyyy[i] * fe_0 + tr_y_zzz_xyyyy[i] * pa_x[i];

        tr_y_xzzz_xyyyz[i] = tr_y_zzz_yyyz[i] * fe_0 + tr_y_zzz_xyyyz[i] * pa_x[i];

        tr_y_xzzz_xyyzz[i] = tr_y_zzz_yyzz[i] * fe_0 + tr_y_zzz_xyyzz[i] * pa_x[i];

        tr_y_xzzz_xyzzz[i] = tr_y_zzz_yzzz[i] * fe_0 + tr_y_zzz_xyzzz[i] * pa_x[i];

        tr_y_xzzz_xzzzz[i] = tr_y_zzz_zzzz[i] * fe_0 + tr_y_zzz_xzzzz[i] * pa_x[i];

        tr_y_xzzz_yyyyy[i] = tr_y_zzz_yyyyy[i] * pa_x[i];

        tr_y_xzzz_yyyyz[i] = tr_y_zzz_yyyyz[i] * pa_x[i];

        tr_y_xzzz_yyyzz[i] = tr_y_zzz_yyyzz[i] * pa_x[i];

        tr_y_xzzz_yyzzz[i] = tr_y_zzz_yyzzz[i] * pa_x[i];

        tr_y_xzzz_yzzzz[i] = tr_y_zzz_yzzzz[i] * pa_x[i];

        tr_y_xzzz_zzzzz[i] = tr_y_zzz_zzzzz[i] * pa_x[i];
    }

    // Set up 525-546 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, tr_y_yy_xxxxx, tr_y_yy_xxxxy, tr_y_yy_xxxxz, tr_y_yy_xxxyy, tr_y_yy_xxxyz, tr_y_yy_xxxzz, tr_y_yy_xxyyy, tr_y_yy_xxyyz, tr_y_yy_xxyzz, tr_y_yy_xxzzz, tr_y_yy_xyyyy, tr_y_yy_xyyyz, tr_y_yy_xyyzz, tr_y_yy_xyzzz, tr_y_yy_xzzzz, tr_y_yy_yyyyy, tr_y_yy_yyyyz, tr_y_yy_yyyzz, tr_y_yy_yyzzz, tr_y_yy_yzzzz, tr_y_yy_zzzzz, tr_y_yyy_xxxx, tr_y_yyy_xxxxx, tr_y_yyy_xxxxy, tr_y_yyy_xxxxz, tr_y_yyy_xxxy, tr_y_yyy_xxxyy, tr_y_yyy_xxxyz, tr_y_yyy_xxxz, tr_y_yyy_xxxzz, tr_y_yyy_xxyy, tr_y_yyy_xxyyy, tr_y_yyy_xxyyz, tr_y_yyy_xxyz, tr_y_yyy_xxyzz, tr_y_yyy_xxzz, tr_y_yyy_xxzzz, tr_y_yyy_xyyy, tr_y_yyy_xyyyy, tr_y_yyy_xyyyz, tr_y_yyy_xyyz, tr_y_yyy_xyyzz, tr_y_yyy_xyzz, tr_y_yyy_xyzzz, tr_y_yyy_xzzz, tr_y_yyy_xzzzz, tr_y_yyy_yyyy, tr_y_yyy_yyyyy, tr_y_yyy_yyyyz, tr_y_yyy_yyyz, tr_y_yyy_yyyzz, tr_y_yyy_yyzz, tr_y_yyy_yyzzz, tr_y_yyy_yzzz, tr_y_yyy_yzzzz, tr_y_yyy_zzzz, tr_y_yyy_zzzzz, tr_y_yyyy_xxxxx, tr_y_yyyy_xxxxy, tr_y_yyyy_xxxxz, tr_y_yyyy_xxxyy, tr_y_yyyy_xxxyz, tr_y_yyyy_xxxzz, tr_y_yyyy_xxyyy, tr_y_yyyy_xxyyz, tr_y_yyyy_xxyzz, tr_y_yyyy_xxzzz, tr_y_yyyy_xyyyy, tr_y_yyyy_xyyyz, tr_y_yyyy_xyyzz, tr_y_yyyy_xyzzz, tr_y_yyyy_xzzzz, tr_y_yyyy_yyyyy, tr_y_yyyy_yyyyz, tr_y_yyyy_yyyzz, tr_y_yyyy_yyzzz, tr_y_yyyy_yzzzz, tr_y_yyyy_zzzzz, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxzz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xxzzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_xzzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyy_xxxxx[i] = 3.0 * tr_y_yy_xxxxx[i] * fe_0 + ts_yyy_xxxxx[i] * fe_0 + tr_y_yyy_xxxxx[i] * pa_y[i];

        tr_y_yyyy_xxxxy[i] = 3.0 * tr_y_yy_xxxxy[i] * fe_0 + tr_y_yyy_xxxx[i] * fe_0 + ts_yyy_xxxxy[i] * fe_0 + tr_y_yyy_xxxxy[i] * pa_y[i];

        tr_y_yyyy_xxxxz[i] = 3.0 * tr_y_yy_xxxxz[i] * fe_0 + ts_yyy_xxxxz[i] * fe_0 + tr_y_yyy_xxxxz[i] * pa_y[i];

        tr_y_yyyy_xxxyy[i] = 3.0 * tr_y_yy_xxxyy[i] * fe_0 + 2.0 * tr_y_yyy_xxxy[i] * fe_0 + ts_yyy_xxxyy[i] * fe_0 + tr_y_yyy_xxxyy[i] * pa_y[i];

        tr_y_yyyy_xxxyz[i] = 3.0 * tr_y_yy_xxxyz[i] * fe_0 + tr_y_yyy_xxxz[i] * fe_0 + ts_yyy_xxxyz[i] * fe_0 + tr_y_yyy_xxxyz[i] * pa_y[i];

        tr_y_yyyy_xxxzz[i] = 3.0 * tr_y_yy_xxxzz[i] * fe_0 + ts_yyy_xxxzz[i] * fe_0 + tr_y_yyy_xxxzz[i] * pa_y[i];

        tr_y_yyyy_xxyyy[i] = 3.0 * tr_y_yy_xxyyy[i] * fe_0 + 3.0 * tr_y_yyy_xxyy[i] * fe_0 + ts_yyy_xxyyy[i] * fe_0 + tr_y_yyy_xxyyy[i] * pa_y[i];

        tr_y_yyyy_xxyyz[i] = 3.0 * tr_y_yy_xxyyz[i] * fe_0 + 2.0 * tr_y_yyy_xxyz[i] * fe_0 + ts_yyy_xxyyz[i] * fe_0 + tr_y_yyy_xxyyz[i] * pa_y[i];

        tr_y_yyyy_xxyzz[i] = 3.0 * tr_y_yy_xxyzz[i] * fe_0 + tr_y_yyy_xxzz[i] * fe_0 + ts_yyy_xxyzz[i] * fe_0 + tr_y_yyy_xxyzz[i] * pa_y[i];

        tr_y_yyyy_xxzzz[i] = 3.0 * tr_y_yy_xxzzz[i] * fe_0 + ts_yyy_xxzzz[i] * fe_0 + tr_y_yyy_xxzzz[i] * pa_y[i];

        tr_y_yyyy_xyyyy[i] = 3.0 * tr_y_yy_xyyyy[i] * fe_0 + 4.0 * tr_y_yyy_xyyy[i] * fe_0 + ts_yyy_xyyyy[i] * fe_0 + tr_y_yyy_xyyyy[i] * pa_y[i];

        tr_y_yyyy_xyyyz[i] = 3.0 * tr_y_yy_xyyyz[i] * fe_0 + 3.0 * tr_y_yyy_xyyz[i] * fe_0 + ts_yyy_xyyyz[i] * fe_0 + tr_y_yyy_xyyyz[i] * pa_y[i];

        tr_y_yyyy_xyyzz[i] = 3.0 * tr_y_yy_xyyzz[i] * fe_0 + 2.0 * tr_y_yyy_xyzz[i] * fe_0 + ts_yyy_xyyzz[i] * fe_0 + tr_y_yyy_xyyzz[i] * pa_y[i];

        tr_y_yyyy_xyzzz[i] = 3.0 * tr_y_yy_xyzzz[i] * fe_0 + tr_y_yyy_xzzz[i] * fe_0 + ts_yyy_xyzzz[i] * fe_0 + tr_y_yyy_xyzzz[i] * pa_y[i];

        tr_y_yyyy_xzzzz[i] = 3.0 * tr_y_yy_xzzzz[i] * fe_0 + ts_yyy_xzzzz[i] * fe_0 + tr_y_yyy_xzzzz[i] * pa_y[i];

        tr_y_yyyy_yyyyy[i] = 3.0 * tr_y_yy_yyyyy[i] * fe_0 + 5.0 * tr_y_yyy_yyyy[i] * fe_0 + ts_yyy_yyyyy[i] * fe_0 + tr_y_yyy_yyyyy[i] * pa_y[i];

        tr_y_yyyy_yyyyz[i] = 3.0 * tr_y_yy_yyyyz[i] * fe_0 + 4.0 * tr_y_yyy_yyyz[i] * fe_0 + ts_yyy_yyyyz[i] * fe_0 + tr_y_yyy_yyyyz[i] * pa_y[i];

        tr_y_yyyy_yyyzz[i] = 3.0 * tr_y_yy_yyyzz[i] * fe_0 + 3.0 * tr_y_yyy_yyzz[i] * fe_0 + ts_yyy_yyyzz[i] * fe_0 + tr_y_yyy_yyyzz[i] * pa_y[i];

        tr_y_yyyy_yyzzz[i] = 3.0 * tr_y_yy_yyzzz[i] * fe_0 + 2.0 * tr_y_yyy_yzzz[i] * fe_0 + ts_yyy_yyzzz[i] * fe_0 + tr_y_yyy_yyzzz[i] * pa_y[i];

        tr_y_yyyy_yzzzz[i] = 3.0 * tr_y_yy_yzzzz[i] * fe_0 + tr_y_yyy_zzzz[i] * fe_0 + ts_yyy_yzzzz[i] * fe_0 + tr_y_yyy_yzzzz[i] * pa_y[i];

        tr_y_yyyy_zzzzz[i] = 3.0 * tr_y_yy_zzzzz[i] * fe_0 + ts_yyy_zzzzz[i] * fe_0 + tr_y_yyy_zzzzz[i] * pa_y[i];
    }

    // Set up 546-567 components of targeted buffer : GH

    auto tr_y_yyyz_xxxxx = pbuffer.data(idx_dip_gh + 546);

    auto tr_y_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 547);

    auto tr_y_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 548);

    auto tr_y_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 549);

    auto tr_y_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 550);

    auto tr_y_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 551);

    auto tr_y_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 552);

    auto tr_y_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 553);

    auto tr_y_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 554);

    auto tr_y_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 555);

    auto tr_y_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 556);

    auto tr_y_yyyz_xyyyz = pbuffer.data(idx_dip_gh + 557);

    auto tr_y_yyyz_xyyzz = pbuffer.data(idx_dip_gh + 558);

    auto tr_y_yyyz_xyzzz = pbuffer.data(idx_dip_gh + 559);

    auto tr_y_yyyz_xzzzz = pbuffer.data(idx_dip_gh + 560);

    auto tr_y_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 561);

    auto tr_y_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 562);

    auto tr_y_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 563);

    auto tr_y_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 564);

    auto tr_y_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 565);

    auto tr_y_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 566);

    #pragma omp simd aligned(pa_z, tr_y_yyy_xxxx, tr_y_yyy_xxxxx, tr_y_yyy_xxxxy, tr_y_yyy_xxxxz, tr_y_yyy_xxxy, tr_y_yyy_xxxyy, tr_y_yyy_xxxyz, tr_y_yyy_xxxz, tr_y_yyy_xxxzz, tr_y_yyy_xxyy, tr_y_yyy_xxyyy, tr_y_yyy_xxyyz, tr_y_yyy_xxyz, tr_y_yyy_xxyzz, tr_y_yyy_xxzz, tr_y_yyy_xxzzz, tr_y_yyy_xyyy, tr_y_yyy_xyyyy, tr_y_yyy_xyyyz, tr_y_yyy_xyyz, tr_y_yyy_xyyzz, tr_y_yyy_xyzz, tr_y_yyy_xyzzz, tr_y_yyy_xzzz, tr_y_yyy_xzzzz, tr_y_yyy_yyyy, tr_y_yyy_yyyyy, tr_y_yyy_yyyyz, tr_y_yyy_yyyz, tr_y_yyy_yyyzz, tr_y_yyy_yyzz, tr_y_yyy_yyzzz, tr_y_yyy_yzzz, tr_y_yyy_yzzzz, tr_y_yyy_zzzz, tr_y_yyy_zzzzz, tr_y_yyyz_xxxxx, tr_y_yyyz_xxxxy, tr_y_yyyz_xxxxz, tr_y_yyyz_xxxyy, tr_y_yyyz_xxxyz, tr_y_yyyz_xxxzz, tr_y_yyyz_xxyyy, tr_y_yyyz_xxyyz, tr_y_yyyz_xxyzz, tr_y_yyyz_xxzzz, tr_y_yyyz_xyyyy, tr_y_yyyz_xyyyz, tr_y_yyyz_xyyzz, tr_y_yyyz_xyzzz, tr_y_yyyz_xzzzz, tr_y_yyyz_yyyyy, tr_y_yyyz_yyyyz, tr_y_yyyz_yyyzz, tr_y_yyyz_yyzzz, tr_y_yyyz_yzzzz, tr_y_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyz_xxxxx[i] = tr_y_yyy_xxxxx[i] * pa_z[i];

        tr_y_yyyz_xxxxy[i] = tr_y_yyy_xxxxy[i] * pa_z[i];

        tr_y_yyyz_xxxxz[i] = tr_y_yyy_xxxx[i] * fe_0 + tr_y_yyy_xxxxz[i] * pa_z[i];

        tr_y_yyyz_xxxyy[i] = tr_y_yyy_xxxyy[i] * pa_z[i];

        tr_y_yyyz_xxxyz[i] = tr_y_yyy_xxxy[i] * fe_0 + tr_y_yyy_xxxyz[i] * pa_z[i];

        tr_y_yyyz_xxxzz[i] = 2.0 * tr_y_yyy_xxxz[i] * fe_0 + tr_y_yyy_xxxzz[i] * pa_z[i];

        tr_y_yyyz_xxyyy[i] = tr_y_yyy_xxyyy[i] * pa_z[i];

        tr_y_yyyz_xxyyz[i] = tr_y_yyy_xxyy[i] * fe_0 + tr_y_yyy_xxyyz[i] * pa_z[i];

        tr_y_yyyz_xxyzz[i] = 2.0 * tr_y_yyy_xxyz[i] * fe_0 + tr_y_yyy_xxyzz[i] * pa_z[i];

        tr_y_yyyz_xxzzz[i] = 3.0 * tr_y_yyy_xxzz[i] * fe_0 + tr_y_yyy_xxzzz[i] * pa_z[i];

        tr_y_yyyz_xyyyy[i] = tr_y_yyy_xyyyy[i] * pa_z[i];

        tr_y_yyyz_xyyyz[i] = tr_y_yyy_xyyy[i] * fe_0 + tr_y_yyy_xyyyz[i] * pa_z[i];

        tr_y_yyyz_xyyzz[i] = 2.0 * tr_y_yyy_xyyz[i] * fe_0 + tr_y_yyy_xyyzz[i] * pa_z[i];

        tr_y_yyyz_xyzzz[i] = 3.0 * tr_y_yyy_xyzz[i] * fe_0 + tr_y_yyy_xyzzz[i] * pa_z[i];

        tr_y_yyyz_xzzzz[i] = 4.0 * tr_y_yyy_xzzz[i] * fe_0 + tr_y_yyy_xzzzz[i] * pa_z[i];

        tr_y_yyyz_yyyyy[i] = tr_y_yyy_yyyyy[i] * pa_z[i];

        tr_y_yyyz_yyyyz[i] = tr_y_yyy_yyyy[i] * fe_0 + tr_y_yyy_yyyyz[i] * pa_z[i];

        tr_y_yyyz_yyyzz[i] = 2.0 * tr_y_yyy_yyyz[i] * fe_0 + tr_y_yyy_yyyzz[i] * pa_z[i];

        tr_y_yyyz_yyzzz[i] = 3.0 * tr_y_yyy_yyzz[i] * fe_0 + tr_y_yyy_yyzzz[i] * pa_z[i];

        tr_y_yyyz_yzzzz[i] = 4.0 * tr_y_yyy_yzzz[i] * fe_0 + tr_y_yyy_yzzzz[i] * pa_z[i];

        tr_y_yyyz_zzzzz[i] = 5.0 * tr_y_yyy_zzzz[i] * fe_0 + tr_y_yyy_zzzzz[i] * pa_z[i];
    }

    // Set up 567-588 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yy_xxxxx, tr_y_yy_xxxxy, tr_y_yy_xxxyy, tr_y_yy_xxxyz, tr_y_yy_xxyyy, tr_y_yy_xxyyz, tr_y_yy_xxyzz, tr_y_yy_xyyyy, tr_y_yy_xyyyz, tr_y_yy_xyyzz, tr_y_yy_xyzzz, tr_y_yy_yyyyy, tr_y_yy_yyyyz, tr_y_yy_yyyzz, tr_y_yy_yyzzz, tr_y_yy_yzzzz, tr_y_yyz_xxxxx, tr_y_yyz_xxxxy, tr_y_yyz_xxxy, tr_y_yyz_xxxyy, tr_y_yyz_xxxyz, tr_y_yyz_xxyy, tr_y_yyz_xxyyy, tr_y_yyz_xxyyz, tr_y_yyz_xxyz, tr_y_yyz_xxyzz, tr_y_yyz_xyyy, tr_y_yyz_xyyyy, tr_y_yyz_xyyyz, tr_y_yyz_xyyz, tr_y_yyz_xyyzz, tr_y_yyz_xyzz, tr_y_yyz_xyzzz, tr_y_yyz_yyyy, tr_y_yyz_yyyyy, tr_y_yyz_yyyyz, tr_y_yyz_yyyz, tr_y_yyz_yyyzz, tr_y_yyz_yyzz, tr_y_yyz_yyzzz, tr_y_yyz_yzzz, tr_y_yyz_yzzzz, tr_y_yyzz_xxxxx, tr_y_yyzz_xxxxy, tr_y_yyzz_xxxxz, tr_y_yyzz_xxxyy, tr_y_yyzz_xxxyz, tr_y_yyzz_xxxzz, tr_y_yyzz_xxyyy, tr_y_yyzz_xxyyz, tr_y_yyzz_xxyzz, tr_y_yyzz_xxzzz, tr_y_yyzz_xyyyy, tr_y_yyzz_xyyyz, tr_y_yyzz_xyyzz, tr_y_yyzz_xyzzz, tr_y_yyzz_xzzzz, tr_y_yyzz_yyyyy, tr_y_yyzz_yyyyz, tr_y_yyzz_yyyzz, tr_y_yyzz_yyzzz, tr_y_yyzz_yzzzz, tr_y_yyzz_zzzzz, tr_y_yzz_xxxxz, tr_y_yzz_xxxzz, tr_y_yzz_xxzzz, tr_y_yzz_xzzzz, tr_y_yzz_zzzzz, tr_y_zz_xxxxz, tr_y_zz_xxxzz, tr_y_zz_xxzzz, tr_y_zz_xzzzz, tr_y_zz_zzzzz, ts_yzz_xxxxz, ts_yzz_xxxzz, ts_yzz_xxzzz, ts_yzz_xzzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzz_xxxxx[i] = tr_y_yy_xxxxx[i] * fe_0 + tr_y_yyz_xxxxx[i] * pa_z[i];

        tr_y_yyzz_xxxxy[i] = tr_y_yy_xxxxy[i] * fe_0 + tr_y_yyz_xxxxy[i] * pa_z[i];

        tr_y_yyzz_xxxxz[i] = tr_y_zz_xxxxz[i] * fe_0 + ts_yzz_xxxxz[i] * fe_0 + tr_y_yzz_xxxxz[i] * pa_y[i];

        tr_y_yyzz_xxxyy[i] = tr_y_yy_xxxyy[i] * fe_0 + tr_y_yyz_xxxyy[i] * pa_z[i];

        tr_y_yyzz_xxxyz[i] = tr_y_yy_xxxyz[i] * fe_0 + tr_y_yyz_xxxy[i] * fe_0 + tr_y_yyz_xxxyz[i] * pa_z[i];

        tr_y_yyzz_xxxzz[i] = tr_y_zz_xxxzz[i] * fe_0 + ts_yzz_xxxzz[i] * fe_0 + tr_y_yzz_xxxzz[i] * pa_y[i];

        tr_y_yyzz_xxyyy[i] = tr_y_yy_xxyyy[i] * fe_0 + tr_y_yyz_xxyyy[i] * pa_z[i];

        tr_y_yyzz_xxyyz[i] = tr_y_yy_xxyyz[i] * fe_0 + tr_y_yyz_xxyy[i] * fe_0 + tr_y_yyz_xxyyz[i] * pa_z[i];

        tr_y_yyzz_xxyzz[i] = tr_y_yy_xxyzz[i] * fe_0 + 2.0 * tr_y_yyz_xxyz[i] * fe_0 + tr_y_yyz_xxyzz[i] * pa_z[i];

        tr_y_yyzz_xxzzz[i] = tr_y_zz_xxzzz[i] * fe_0 + ts_yzz_xxzzz[i] * fe_0 + tr_y_yzz_xxzzz[i] * pa_y[i];

        tr_y_yyzz_xyyyy[i] = tr_y_yy_xyyyy[i] * fe_0 + tr_y_yyz_xyyyy[i] * pa_z[i];

        tr_y_yyzz_xyyyz[i] = tr_y_yy_xyyyz[i] * fe_0 + tr_y_yyz_xyyy[i] * fe_0 + tr_y_yyz_xyyyz[i] * pa_z[i];

        tr_y_yyzz_xyyzz[i] = tr_y_yy_xyyzz[i] * fe_0 + 2.0 * tr_y_yyz_xyyz[i] * fe_0 + tr_y_yyz_xyyzz[i] * pa_z[i];

        tr_y_yyzz_xyzzz[i] = tr_y_yy_xyzzz[i] * fe_0 + 3.0 * tr_y_yyz_xyzz[i] * fe_0 + tr_y_yyz_xyzzz[i] * pa_z[i];

        tr_y_yyzz_xzzzz[i] = tr_y_zz_xzzzz[i] * fe_0 + ts_yzz_xzzzz[i] * fe_0 + tr_y_yzz_xzzzz[i] * pa_y[i];

        tr_y_yyzz_yyyyy[i] = tr_y_yy_yyyyy[i] * fe_0 + tr_y_yyz_yyyyy[i] * pa_z[i];

        tr_y_yyzz_yyyyz[i] = tr_y_yy_yyyyz[i] * fe_0 + tr_y_yyz_yyyy[i] * fe_0 + tr_y_yyz_yyyyz[i] * pa_z[i];

        tr_y_yyzz_yyyzz[i] = tr_y_yy_yyyzz[i] * fe_0 + 2.0 * tr_y_yyz_yyyz[i] * fe_0 + tr_y_yyz_yyyzz[i] * pa_z[i];

        tr_y_yyzz_yyzzz[i] = tr_y_yy_yyzzz[i] * fe_0 + 3.0 * tr_y_yyz_yyzz[i] * fe_0 + tr_y_yyz_yyzzz[i] * pa_z[i];

        tr_y_yyzz_yzzzz[i] = tr_y_yy_yzzzz[i] * fe_0 + 4.0 * tr_y_yyz_yzzz[i] * fe_0 + tr_y_yyz_yzzzz[i] * pa_z[i];

        tr_y_yyzz_zzzzz[i] = tr_y_zz_zzzzz[i] * fe_0 + ts_yzz_zzzzz[i] * fe_0 + tr_y_yzz_zzzzz[i] * pa_y[i];
    }

    // Set up 588-609 components of targeted buffer : GH

    auto tr_y_yzzz_xxxxx = pbuffer.data(idx_dip_gh + 588);

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

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yz_xxxxy, tr_y_yz_xxxyy, tr_y_yz_xxyyy, tr_y_yz_xyyyy, tr_y_yz_yyyyy, tr_y_yzz_xxxxy, tr_y_yzz_xxxyy, tr_y_yzz_xxyyy, tr_y_yzz_xyyyy, tr_y_yzz_yyyyy, tr_y_yzzz_xxxxx, tr_y_yzzz_xxxxy, tr_y_yzzz_xxxxz, tr_y_yzzz_xxxyy, tr_y_yzzz_xxxyz, tr_y_yzzz_xxxzz, tr_y_yzzz_xxyyy, tr_y_yzzz_xxyyz, tr_y_yzzz_xxyzz, tr_y_yzzz_xxzzz, tr_y_yzzz_xyyyy, tr_y_yzzz_xyyyz, tr_y_yzzz_xyyzz, tr_y_yzzz_xyzzz, tr_y_yzzz_xzzzz, tr_y_yzzz_yyyyy, tr_y_yzzz_yyyyz, tr_y_yzzz_yyyzz, tr_y_yzzz_yyzzz, tr_y_yzzz_yzzzz, tr_y_yzzz_zzzzz, tr_y_zzz_xxxxx, tr_y_zzz_xxxxz, tr_y_zzz_xxxyz, tr_y_zzz_xxxz, tr_y_zzz_xxxzz, tr_y_zzz_xxyyz, tr_y_zzz_xxyz, tr_y_zzz_xxyzz, tr_y_zzz_xxzz, tr_y_zzz_xxzzz, tr_y_zzz_xyyyz, tr_y_zzz_xyyz, tr_y_zzz_xyyzz, tr_y_zzz_xyzz, tr_y_zzz_xyzzz, tr_y_zzz_xzzz, tr_y_zzz_xzzzz, tr_y_zzz_yyyyz, tr_y_zzz_yyyz, tr_y_zzz_yyyzz, tr_y_zzz_yyzz, tr_y_zzz_yyzzz, tr_y_zzz_yzzz, tr_y_zzz_yzzzz, tr_y_zzz_zzzz, tr_y_zzz_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxz, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzz_xxxxx[i] = ts_zzz_xxxxx[i] * fe_0 + tr_y_zzz_xxxxx[i] * pa_y[i];

        tr_y_yzzz_xxxxy[i] = 2.0 * tr_y_yz_xxxxy[i] * fe_0 + tr_y_yzz_xxxxy[i] * pa_z[i];

        tr_y_yzzz_xxxxz[i] = ts_zzz_xxxxz[i] * fe_0 + tr_y_zzz_xxxxz[i] * pa_y[i];

        tr_y_yzzz_xxxyy[i] = 2.0 * tr_y_yz_xxxyy[i] * fe_0 + tr_y_yzz_xxxyy[i] * pa_z[i];

        tr_y_yzzz_xxxyz[i] = tr_y_zzz_xxxz[i] * fe_0 + ts_zzz_xxxyz[i] * fe_0 + tr_y_zzz_xxxyz[i] * pa_y[i];

        tr_y_yzzz_xxxzz[i] = ts_zzz_xxxzz[i] * fe_0 + tr_y_zzz_xxxzz[i] * pa_y[i];

        tr_y_yzzz_xxyyy[i] = 2.0 * tr_y_yz_xxyyy[i] * fe_0 + tr_y_yzz_xxyyy[i] * pa_z[i];

        tr_y_yzzz_xxyyz[i] = 2.0 * tr_y_zzz_xxyz[i] * fe_0 + ts_zzz_xxyyz[i] * fe_0 + tr_y_zzz_xxyyz[i] * pa_y[i];

        tr_y_yzzz_xxyzz[i] = tr_y_zzz_xxzz[i] * fe_0 + ts_zzz_xxyzz[i] * fe_0 + tr_y_zzz_xxyzz[i] * pa_y[i];

        tr_y_yzzz_xxzzz[i] = ts_zzz_xxzzz[i] * fe_0 + tr_y_zzz_xxzzz[i] * pa_y[i];

        tr_y_yzzz_xyyyy[i] = 2.0 * tr_y_yz_xyyyy[i] * fe_0 + tr_y_yzz_xyyyy[i] * pa_z[i];

        tr_y_yzzz_xyyyz[i] = 3.0 * tr_y_zzz_xyyz[i] * fe_0 + ts_zzz_xyyyz[i] * fe_0 + tr_y_zzz_xyyyz[i] * pa_y[i];

        tr_y_yzzz_xyyzz[i] = 2.0 * tr_y_zzz_xyzz[i] * fe_0 + ts_zzz_xyyzz[i] * fe_0 + tr_y_zzz_xyyzz[i] * pa_y[i];

        tr_y_yzzz_xyzzz[i] = tr_y_zzz_xzzz[i] * fe_0 + ts_zzz_xyzzz[i] * fe_0 + tr_y_zzz_xyzzz[i] * pa_y[i];

        tr_y_yzzz_xzzzz[i] = ts_zzz_xzzzz[i] * fe_0 + tr_y_zzz_xzzzz[i] * pa_y[i];

        tr_y_yzzz_yyyyy[i] = 2.0 * tr_y_yz_yyyyy[i] * fe_0 + tr_y_yzz_yyyyy[i] * pa_z[i];

        tr_y_yzzz_yyyyz[i] = 4.0 * tr_y_zzz_yyyz[i] * fe_0 + ts_zzz_yyyyz[i] * fe_0 + tr_y_zzz_yyyyz[i] * pa_y[i];

        tr_y_yzzz_yyyzz[i] = 3.0 * tr_y_zzz_yyzz[i] * fe_0 + ts_zzz_yyyzz[i] * fe_0 + tr_y_zzz_yyyzz[i] * pa_y[i];

        tr_y_yzzz_yyzzz[i] = 2.0 * tr_y_zzz_yzzz[i] * fe_0 + ts_zzz_yyzzz[i] * fe_0 + tr_y_zzz_yyzzz[i] * pa_y[i];

        tr_y_yzzz_yzzzz[i] = tr_y_zzz_zzzz[i] * fe_0 + ts_zzz_yzzzz[i] * fe_0 + tr_y_zzz_yzzzz[i] * pa_y[i];

        tr_y_yzzz_zzzzz[i] = ts_zzz_zzzzz[i] * fe_0 + tr_y_zzz_zzzzz[i] * pa_y[i];
    }

    // Set up 609-630 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_z, tr_y_zz_xxxxx, tr_y_zz_xxxxy, tr_y_zz_xxxxz, tr_y_zz_xxxyy, tr_y_zz_xxxyz, tr_y_zz_xxxzz, tr_y_zz_xxyyy, tr_y_zz_xxyyz, tr_y_zz_xxyzz, tr_y_zz_xxzzz, tr_y_zz_xyyyy, tr_y_zz_xyyyz, tr_y_zz_xyyzz, tr_y_zz_xyzzz, tr_y_zz_xzzzz, tr_y_zz_yyyyy, tr_y_zz_yyyyz, tr_y_zz_yyyzz, tr_y_zz_yyzzz, tr_y_zz_yzzzz, tr_y_zz_zzzzz, tr_y_zzz_xxxx, tr_y_zzz_xxxxx, tr_y_zzz_xxxxy, tr_y_zzz_xxxxz, tr_y_zzz_xxxy, tr_y_zzz_xxxyy, tr_y_zzz_xxxyz, tr_y_zzz_xxxz, tr_y_zzz_xxxzz, tr_y_zzz_xxyy, tr_y_zzz_xxyyy, tr_y_zzz_xxyyz, tr_y_zzz_xxyz, tr_y_zzz_xxyzz, tr_y_zzz_xxzz, tr_y_zzz_xxzzz, tr_y_zzz_xyyy, tr_y_zzz_xyyyy, tr_y_zzz_xyyyz, tr_y_zzz_xyyz, tr_y_zzz_xyyzz, tr_y_zzz_xyzz, tr_y_zzz_xyzzz, tr_y_zzz_xzzz, tr_y_zzz_xzzzz, tr_y_zzz_yyyy, tr_y_zzz_yyyyy, tr_y_zzz_yyyyz, tr_y_zzz_yyyz, tr_y_zzz_yyyzz, tr_y_zzz_yyzz, tr_y_zzz_yyzzz, tr_y_zzz_yzzz, tr_y_zzz_yzzzz, tr_y_zzz_zzzz, tr_y_zzz_zzzzz, tr_y_zzzz_xxxxx, tr_y_zzzz_xxxxy, tr_y_zzzz_xxxxz, tr_y_zzzz_xxxyy, tr_y_zzzz_xxxyz, tr_y_zzzz_xxxzz, tr_y_zzzz_xxyyy, tr_y_zzzz_xxyyz, tr_y_zzzz_xxyzz, tr_y_zzzz_xxzzz, tr_y_zzzz_xyyyy, tr_y_zzzz_xyyyz, tr_y_zzzz_xyyzz, tr_y_zzzz_xyzzz, tr_y_zzzz_xzzzz, tr_y_zzzz_yyyyy, tr_y_zzzz_yyyyz, tr_y_zzzz_yyyzz, tr_y_zzzz_yyzzz, tr_y_zzzz_yzzzz, tr_y_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzz_xxxxx[i] = 3.0 * tr_y_zz_xxxxx[i] * fe_0 + tr_y_zzz_xxxxx[i] * pa_z[i];

        tr_y_zzzz_xxxxy[i] = 3.0 * tr_y_zz_xxxxy[i] * fe_0 + tr_y_zzz_xxxxy[i] * pa_z[i];

        tr_y_zzzz_xxxxz[i] = 3.0 * tr_y_zz_xxxxz[i] * fe_0 + tr_y_zzz_xxxx[i] * fe_0 + tr_y_zzz_xxxxz[i] * pa_z[i];

        tr_y_zzzz_xxxyy[i] = 3.0 * tr_y_zz_xxxyy[i] * fe_0 + tr_y_zzz_xxxyy[i] * pa_z[i];

        tr_y_zzzz_xxxyz[i] = 3.0 * tr_y_zz_xxxyz[i] * fe_0 + tr_y_zzz_xxxy[i] * fe_0 + tr_y_zzz_xxxyz[i] * pa_z[i];

        tr_y_zzzz_xxxzz[i] = 3.0 * tr_y_zz_xxxzz[i] * fe_0 + 2.0 * tr_y_zzz_xxxz[i] * fe_0 + tr_y_zzz_xxxzz[i] * pa_z[i];

        tr_y_zzzz_xxyyy[i] = 3.0 * tr_y_zz_xxyyy[i] * fe_0 + tr_y_zzz_xxyyy[i] * pa_z[i];

        tr_y_zzzz_xxyyz[i] = 3.0 * tr_y_zz_xxyyz[i] * fe_0 + tr_y_zzz_xxyy[i] * fe_0 + tr_y_zzz_xxyyz[i] * pa_z[i];

        tr_y_zzzz_xxyzz[i] = 3.0 * tr_y_zz_xxyzz[i] * fe_0 + 2.0 * tr_y_zzz_xxyz[i] * fe_0 + tr_y_zzz_xxyzz[i] * pa_z[i];

        tr_y_zzzz_xxzzz[i] = 3.0 * tr_y_zz_xxzzz[i] * fe_0 + 3.0 * tr_y_zzz_xxzz[i] * fe_0 + tr_y_zzz_xxzzz[i] * pa_z[i];

        tr_y_zzzz_xyyyy[i] = 3.0 * tr_y_zz_xyyyy[i] * fe_0 + tr_y_zzz_xyyyy[i] * pa_z[i];

        tr_y_zzzz_xyyyz[i] = 3.0 * tr_y_zz_xyyyz[i] * fe_0 + tr_y_zzz_xyyy[i] * fe_0 + tr_y_zzz_xyyyz[i] * pa_z[i];

        tr_y_zzzz_xyyzz[i] = 3.0 * tr_y_zz_xyyzz[i] * fe_0 + 2.0 * tr_y_zzz_xyyz[i] * fe_0 + tr_y_zzz_xyyzz[i] * pa_z[i];

        tr_y_zzzz_xyzzz[i] = 3.0 * tr_y_zz_xyzzz[i] * fe_0 + 3.0 * tr_y_zzz_xyzz[i] * fe_0 + tr_y_zzz_xyzzz[i] * pa_z[i];

        tr_y_zzzz_xzzzz[i] = 3.0 * tr_y_zz_xzzzz[i] * fe_0 + 4.0 * tr_y_zzz_xzzz[i] * fe_0 + tr_y_zzz_xzzzz[i] * pa_z[i];

        tr_y_zzzz_yyyyy[i] = 3.0 * tr_y_zz_yyyyy[i] * fe_0 + tr_y_zzz_yyyyy[i] * pa_z[i];

        tr_y_zzzz_yyyyz[i] = 3.0 * tr_y_zz_yyyyz[i] * fe_0 + tr_y_zzz_yyyy[i] * fe_0 + tr_y_zzz_yyyyz[i] * pa_z[i];

        tr_y_zzzz_yyyzz[i] = 3.0 * tr_y_zz_yyyzz[i] * fe_0 + 2.0 * tr_y_zzz_yyyz[i] * fe_0 + tr_y_zzz_yyyzz[i] * pa_z[i];

        tr_y_zzzz_yyzzz[i] = 3.0 * tr_y_zz_yyzzz[i] * fe_0 + 3.0 * tr_y_zzz_yyzz[i] * fe_0 + tr_y_zzz_yyzzz[i] * pa_z[i];

        tr_y_zzzz_yzzzz[i] = 3.0 * tr_y_zz_yzzzz[i] * fe_0 + 4.0 * tr_y_zzz_yzzz[i] * fe_0 + tr_y_zzz_yzzzz[i] * pa_z[i];

        tr_y_zzzz_zzzzz[i] = 3.0 * tr_y_zz_zzzzz[i] * fe_0 + 5.0 * tr_y_zzz_zzzz[i] * fe_0 + tr_y_zzz_zzzzz[i] * pa_z[i];
    }

    // Set up 630-651 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_z_xx_xxxxx, tr_z_xx_xxxxy, tr_z_xx_xxxxz, tr_z_xx_xxxyy, tr_z_xx_xxxyz, tr_z_xx_xxxzz, tr_z_xx_xxyyy, tr_z_xx_xxyyz, tr_z_xx_xxyzz, tr_z_xx_xxzzz, tr_z_xx_xyyyy, tr_z_xx_xyyyz, tr_z_xx_xyyzz, tr_z_xx_xyzzz, tr_z_xx_xzzzz, tr_z_xx_yyyyy, tr_z_xx_yyyyz, tr_z_xx_yyyzz, tr_z_xx_yyzzz, tr_z_xx_yzzzz, tr_z_xx_zzzzz, tr_z_xxx_xxxx, tr_z_xxx_xxxxx, tr_z_xxx_xxxxy, tr_z_xxx_xxxxz, tr_z_xxx_xxxy, tr_z_xxx_xxxyy, tr_z_xxx_xxxyz, tr_z_xxx_xxxz, tr_z_xxx_xxxzz, tr_z_xxx_xxyy, tr_z_xxx_xxyyy, tr_z_xxx_xxyyz, tr_z_xxx_xxyz, tr_z_xxx_xxyzz, tr_z_xxx_xxzz, tr_z_xxx_xxzzz, tr_z_xxx_xyyy, tr_z_xxx_xyyyy, tr_z_xxx_xyyyz, tr_z_xxx_xyyz, tr_z_xxx_xyyzz, tr_z_xxx_xyzz, tr_z_xxx_xyzzz, tr_z_xxx_xzzz, tr_z_xxx_xzzzz, tr_z_xxx_yyyy, tr_z_xxx_yyyyy, tr_z_xxx_yyyyz, tr_z_xxx_yyyz, tr_z_xxx_yyyzz, tr_z_xxx_yyzz, tr_z_xxx_yyzzz, tr_z_xxx_yzzz, tr_z_xxx_yzzzz, tr_z_xxx_zzzz, tr_z_xxx_zzzzz, tr_z_xxxx_xxxxx, tr_z_xxxx_xxxxy, tr_z_xxxx_xxxxz, tr_z_xxxx_xxxyy, tr_z_xxxx_xxxyz, tr_z_xxxx_xxxzz, tr_z_xxxx_xxyyy, tr_z_xxxx_xxyyz, tr_z_xxxx_xxyzz, tr_z_xxxx_xxzzz, tr_z_xxxx_xyyyy, tr_z_xxxx_xyyyz, tr_z_xxxx_xyyzz, tr_z_xxxx_xyzzz, tr_z_xxxx_xzzzz, tr_z_xxxx_yyyyy, tr_z_xxxx_yyyyz, tr_z_xxxx_yyyzz, tr_z_xxxx_yyzzz, tr_z_xxxx_yzzzz, tr_z_xxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxx_xxxxx[i] = 3.0 * tr_z_xx_xxxxx[i] * fe_0 + 5.0 * tr_z_xxx_xxxx[i] * fe_0 + tr_z_xxx_xxxxx[i] * pa_x[i];

        tr_z_xxxx_xxxxy[i] = 3.0 * tr_z_xx_xxxxy[i] * fe_0 + 4.0 * tr_z_xxx_xxxy[i] * fe_0 + tr_z_xxx_xxxxy[i] * pa_x[i];

        tr_z_xxxx_xxxxz[i] = 3.0 * tr_z_xx_xxxxz[i] * fe_0 + 4.0 * tr_z_xxx_xxxz[i] * fe_0 + tr_z_xxx_xxxxz[i] * pa_x[i];

        tr_z_xxxx_xxxyy[i] = 3.0 * tr_z_xx_xxxyy[i] * fe_0 + 3.0 * tr_z_xxx_xxyy[i] * fe_0 + tr_z_xxx_xxxyy[i] * pa_x[i];

        tr_z_xxxx_xxxyz[i] = 3.0 * tr_z_xx_xxxyz[i] * fe_0 + 3.0 * tr_z_xxx_xxyz[i] * fe_0 + tr_z_xxx_xxxyz[i] * pa_x[i];

        tr_z_xxxx_xxxzz[i] = 3.0 * tr_z_xx_xxxzz[i] * fe_0 + 3.0 * tr_z_xxx_xxzz[i] * fe_0 + tr_z_xxx_xxxzz[i] * pa_x[i];

        tr_z_xxxx_xxyyy[i] = 3.0 * tr_z_xx_xxyyy[i] * fe_0 + 2.0 * tr_z_xxx_xyyy[i] * fe_0 + tr_z_xxx_xxyyy[i] * pa_x[i];

        tr_z_xxxx_xxyyz[i] = 3.0 * tr_z_xx_xxyyz[i] * fe_0 + 2.0 * tr_z_xxx_xyyz[i] * fe_0 + tr_z_xxx_xxyyz[i] * pa_x[i];

        tr_z_xxxx_xxyzz[i] = 3.0 * tr_z_xx_xxyzz[i] * fe_0 + 2.0 * tr_z_xxx_xyzz[i] * fe_0 + tr_z_xxx_xxyzz[i] * pa_x[i];

        tr_z_xxxx_xxzzz[i] = 3.0 * tr_z_xx_xxzzz[i] * fe_0 + 2.0 * tr_z_xxx_xzzz[i] * fe_0 + tr_z_xxx_xxzzz[i] * pa_x[i];

        tr_z_xxxx_xyyyy[i] = 3.0 * tr_z_xx_xyyyy[i] * fe_0 + tr_z_xxx_yyyy[i] * fe_0 + tr_z_xxx_xyyyy[i] * pa_x[i];

        tr_z_xxxx_xyyyz[i] = 3.0 * tr_z_xx_xyyyz[i] * fe_0 + tr_z_xxx_yyyz[i] * fe_0 + tr_z_xxx_xyyyz[i] * pa_x[i];

        tr_z_xxxx_xyyzz[i] = 3.0 * tr_z_xx_xyyzz[i] * fe_0 + tr_z_xxx_yyzz[i] * fe_0 + tr_z_xxx_xyyzz[i] * pa_x[i];

        tr_z_xxxx_xyzzz[i] = 3.0 * tr_z_xx_xyzzz[i] * fe_0 + tr_z_xxx_yzzz[i] * fe_0 + tr_z_xxx_xyzzz[i] * pa_x[i];

        tr_z_xxxx_xzzzz[i] = 3.0 * tr_z_xx_xzzzz[i] * fe_0 + tr_z_xxx_zzzz[i] * fe_0 + tr_z_xxx_xzzzz[i] * pa_x[i];

        tr_z_xxxx_yyyyy[i] = 3.0 * tr_z_xx_yyyyy[i] * fe_0 + tr_z_xxx_yyyyy[i] * pa_x[i];

        tr_z_xxxx_yyyyz[i] = 3.0 * tr_z_xx_yyyyz[i] * fe_0 + tr_z_xxx_yyyyz[i] * pa_x[i];

        tr_z_xxxx_yyyzz[i] = 3.0 * tr_z_xx_yyyzz[i] * fe_0 + tr_z_xxx_yyyzz[i] * pa_x[i];

        tr_z_xxxx_yyzzz[i] = 3.0 * tr_z_xx_yyzzz[i] * fe_0 + tr_z_xxx_yyzzz[i] * pa_x[i];

        tr_z_xxxx_yzzzz[i] = 3.0 * tr_z_xx_yzzzz[i] * fe_0 + tr_z_xxx_yzzzz[i] * pa_x[i];

        tr_z_xxxx_zzzzz[i] = 3.0 * tr_z_xx_zzzzz[i] * fe_0 + tr_z_xxx_zzzzz[i] * pa_x[i];
    }

    // Set up 651-672 components of targeted buffer : GH

    auto tr_z_xxxy_xxxxx = pbuffer.data(idx_dip_gh + 651);

    auto tr_z_xxxy_xxxxy = pbuffer.data(idx_dip_gh + 652);

    auto tr_z_xxxy_xxxxz = pbuffer.data(idx_dip_gh + 653);

    auto tr_z_xxxy_xxxyy = pbuffer.data(idx_dip_gh + 654);

    auto tr_z_xxxy_xxxyz = pbuffer.data(idx_dip_gh + 655);

    auto tr_z_xxxy_xxxzz = pbuffer.data(idx_dip_gh + 656);

    auto tr_z_xxxy_xxyyy = pbuffer.data(idx_dip_gh + 657);

    auto tr_z_xxxy_xxyyz = pbuffer.data(idx_dip_gh + 658);

    auto tr_z_xxxy_xxyzz = pbuffer.data(idx_dip_gh + 659);

    auto tr_z_xxxy_xxzzz = pbuffer.data(idx_dip_gh + 660);

    auto tr_z_xxxy_xyyyy = pbuffer.data(idx_dip_gh + 661);

    auto tr_z_xxxy_xyyyz = pbuffer.data(idx_dip_gh + 662);

    auto tr_z_xxxy_xyyzz = pbuffer.data(idx_dip_gh + 663);

    auto tr_z_xxxy_xyzzz = pbuffer.data(idx_dip_gh + 664);

    auto tr_z_xxxy_xzzzz = pbuffer.data(idx_dip_gh + 665);

    auto tr_z_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 666);

    auto tr_z_xxxy_yyyyz = pbuffer.data(idx_dip_gh + 667);

    auto tr_z_xxxy_yyyzz = pbuffer.data(idx_dip_gh + 668);

    auto tr_z_xxxy_yyzzz = pbuffer.data(idx_dip_gh + 669);

    auto tr_z_xxxy_yzzzz = pbuffer.data(idx_dip_gh + 670);

    auto tr_z_xxxy_zzzzz = pbuffer.data(idx_dip_gh + 671);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxx_xxxx, tr_z_xxx_xxxxx, tr_z_xxx_xxxxy, tr_z_xxx_xxxxz, tr_z_xxx_xxxy, tr_z_xxx_xxxyy, tr_z_xxx_xxxyz, tr_z_xxx_xxxz, tr_z_xxx_xxxzz, tr_z_xxx_xxyy, tr_z_xxx_xxyyy, tr_z_xxx_xxyyz, tr_z_xxx_xxyz, tr_z_xxx_xxyzz, tr_z_xxx_xxzz, tr_z_xxx_xxzzz, tr_z_xxx_xyyy, tr_z_xxx_xyyyy, tr_z_xxx_xyyyz, tr_z_xxx_xyyz, tr_z_xxx_xyyzz, tr_z_xxx_xyzz, tr_z_xxx_xyzzz, tr_z_xxx_xzzz, tr_z_xxx_xzzzz, tr_z_xxx_zzzzz, tr_z_xxxy_xxxxx, tr_z_xxxy_xxxxy, tr_z_xxxy_xxxxz, tr_z_xxxy_xxxyy, tr_z_xxxy_xxxyz, tr_z_xxxy_xxxzz, tr_z_xxxy_xxyyy, tr_z_xxxy_xxyyz, tr_z_xxxy_xxyzz, tr_z_xxxy_xxzzz, tr_z_xxxy_xyyyy, tr_z_xxxy_xyyyz, tr_z_xxxy_xyyzz, tr_z_xxxy_xyzzz, tr_z_xxxy_xzzzz, tr_z_xxxy_yyyyy, tr_z_xxxy_yyyyz, tr_z_xxxy_yyyzz, tr_z_xxxy_yyzzz, tr_z_xxxy_yzzzz, tr_z_xxxy_zzzzz, tr_z_xxy_yyyyy, tr_z_xxy_yyyyz, tr_z_xxy_yyyzz, tr_z_xxy_yyzzz, tr_z_xxy_yzzzz, tr_z_xy_yyyyy, tr_z_xy_yyyyz, tr_z_xy_yyyzz, tr_z_xy_yyzzz, tr_z_xy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxy_xxxxx[i] = tr_z_xxx_xxxxx[i] * pa_y[i];

        tr_z_xxxy_xxxxy[i] = tr_z_xxx_xxxx[i] * fe_0 + tr_z_xxx_xxxxy[i] * pa_y[i];

        tr_z_xxxy_xxxxz[i] = tr_z_xxx_xxxxz[i] * pa_y[i];

        tr_z_xxxy_xxxyy[i] = 2.0 * tr_z_xxx_xxxy[i] * fe_0 + tr_z_xxx_xxxyy[i] * pa_y[i];

        tr_z_xxxy_xxxyz[i] = tr_z_xxx_xxxz[i] * fe_0 + tr_z_xxx_xxxyz[i] * pa_y[i];

        tr_z_xxxy_xxxzz[i] = tr_z_xxx_xxxzz[i] * pa_y[i];

        tr_z_xxxy_xxyyy[i] = 3.0 * tr_z_xxx_xxyy[i] * fe_0 + tr_z_xxx_xxyyy[i] * pa_y[i];

        tr_z_xxxy_xxyyz[i] = 2.0 * tr_z_xxx_xxyz[i] * fe_0 + tr_z_xxx_xxyyz[i] * pa_y[i];

        tr_z_xxxy_xxyzz[i] = tr_z_xxx_xxzz[i] * fe_0 + tr_z_xxx_xxyzz[i] * pa_y[i];

        tr_z_xxxy_xxzzz[i] = tr_z_xxx_xxzzz[i] * pa_y[i];

        tr_z_xxxy_xyyyy[i] = 4.0 * tr_z_xxx_xyyy[i] * fe_0 + tr_z_xxx_xyyyy[i] * pa_y[i];

        tr_z_xxxy_xyyyz[i] = 3.0 * tr_z_xxx_xyyz[i] * fe_0 + tr_z_xxx_xyyyz[i] * pa_y[i];

        tr_z_xxxy_xyyzz[i] = 2.0 * tr_z_xxx_xyzz[i] * fe_0 + tr_z_xxx_xyyzz[i] * pa_y[i];

        tr_z_xxxy_xyzzz[i] = tr_z_xxx_xzzz[i] * fe_0 + tr_z_xxx_xyzzz[i] * pa_y[i];

        tr_z_xxxy_xzzzz[i] = tr_z_xxx_xzzzz[i] * pa_y[i];

        tr_z_xxxy_yyyyy[i] = 2.0 * tr_z_xy_yyyyy[i] * fe_0 + tr_z_xxy_yyyyy[i] * pa_x[i];

        tr_z_xxxy_yyyyz[i] = 2.0 * tr_z_xy_yyyyz[i] * fe_0 + tr_z_xxy_yyyyz[i] * pa_x[i];

        tr_z_xxxy_yyyzz[i] = 2.0 * tr_z_xy_yyyzz[i] * fe_0 + tr_z_xxy_yyyzz[i] * pa_x[i];

        tr_z_xxxy_yyzzz[i] = 2.0 * tr_z_xy_yyzzz[i] * fe_0 + tr_z_xxy_yyzzz[i] * pa_x[i];

        tr_z_xxxy_yzzzz[i] = 2.0 * tr_z_xy_yzzzz[i] * fe_0 + tr_z_xxy_yzzzz[i] * pa_x[i];

        tr_z_xxxy_zzzzz[i] = tr_z_xxx_zzzzz[i] * pa_y[i];
    }

    // Set up 672-693 components of targeted buffer : GH

    auto tr_z_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 672);

    auto tr_z_xxxz_xxxxy = pbuffer.data(idx_dip_gh + 673);

    auto tr_z_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 674);

    auto tr_z_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 675);

    auto tr_z_xxxz_xxxyz = pbuffer.data(idx_dip_gh + 676);

    auto tr_z_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 677);

    auto tr_z_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 678);

    auto tr_z_xxxz_xxyyz = pbuffer.data(idx_dip_gh + 679);

    auto tr_z_xxxz_xxyzz = pbuffer.data(idx_dip_gh + 680);

    auto tr_z_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 681);

    auto tr_z_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 682);

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

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xxx_xxxxx, tr_z_xxx_xxxxy, tr_z_xxx_xxxyy, tr_z_xxx_xxyyy, tr_z_xxx_xyyyy, tr_z_xxxz_xxxxx, tr_z_xxxz_xxxxy, tr_z_xxxz_xxxxz, tr_z_xxxz_xxxyy, tr_z_xxxz_xxxyz, tr_z_xxxz_xxxzz, tr_z_xxxz_xxyyy, tr_z_xxxz_xxyyz, tr_z_xxxz_xxyzz, tr_z_xxxz_xxzzz, tr_z_xxxz_xyyyy, tr_z_xxxz_xyyyz, tr_z_xxxz_xyyzz, tr_z_xxxz_xyzzz, tr_z_xxxz_xzzzz, tr_z_xxxz_yyyyy, tr_z_xxxz_yyyyz, tr_z_xxxz_yyyzz, tr_z_xxxz_yyzzz, tr_z_xxxz_yzzzz, tr_z_xxxz_zzzzz, tr_z_xxz_xxxxz, tr_z_xxz_xxxyz, tr_z_xxz_xxxz, tr_z_xxz_xxxzz, tr_z_xxz_xxyyz, tr_z_xxz_xxyz, tr_z_xxz_xxyzz, tr_z_xxz_xxzz, tr_z_xxz_xxzzz, tr_z_xxz_xyyyz, tr_z_xxz_xyyz, tr_z_xxz_xyyzz, tr_z_xxz_xyzz, tr_z_xxz_xyzzz, tr_z_xxz_xzzz, tr_z_xxz_xzzzz, tr_z_xxz_yyyyy, tr_z_xxz_yyyyz, tr_z_xxz_yyyz, tr_z_xxz_yyyzz, tr_z_xxz_yyzz, tr_z_xxz_yyzzz, tr_z_xxz_yzzz, tr_z_xxz_yzzzz, tr_z_xxz_zzzz, tr_z_xxz_zzzzz, tr_z_xz_xxxxz, tr_z_xz_xxxyz, tr_z_xz_xxxzz, tr_z_xz_xxyyz, tr_z_xz_xxyzz, tr_z_xz_xxzzz, tr_z_xz_xyyyz, tr_z_xz_xyyzz, tr_z_xz_xyzzz, tr_z_xz_xzzzz, tr_z_xz_yyyyy, tr_z_xz_yyyyz, tr_z_xz_yyyzz, tr_z_xz_yyzzz, tr_z_xz_yzzzz, tr_z_xz_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxyy, ts_xxx_xxyyy, ts_xxx_xyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxz_xxxxx[i] = ts_xxx_xxxxx[i] * fe_0 + tr_z_xxx_xxxxx[i] * pa_z[i];

        tr_z_xxxz_xxxxy[i] = ts_xxx_xxxxy[i] * fe_0 + tr_z_xxx_xxxxy[i] * pa_z[i];

        tr_z_xxxz_xxxxz[i] = 2.0 * tr_z_xz_xxxxz[i] * fe_0 + 4.0 * tr_z_xxz_xxxz[i] * fe_0 + tr_z_xxz_xxxxz[i] * pa_x[i];

        tr_z_xxxz_xxxyy[i] = ts_xxx_xxxyy[i] * fe_0 + tr_z_xxx_xxxyy[i] * pa_z[i];

        tr_z_xxxz_xxxyz[i] = 2.0 * tr_z_xz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxz_xxyz[i] * fe_0 + tr_z_xxz_xxxyz[i] * pa_x[i];

        tr_z_xxxz_xxxzz[i] = 2.0 * tr_z_xz_xxxzz[i] * fe_0 + 3.0 * tr_z_xxz_xxzz[i] * fe_0 + tr_z_xxz_xxxzz[i] * pa_x[i];

        tr_z_xxxz_xxyyy[i] = ts_xxx_xxyyy[i] * fe_0 + tr_z_xxx_xxyyy[i] * pa_z[i];

        tr_z_xxxz_xxyyz[i] = 2.0 * tr_z_xz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxz_xyyz[i] * fe_0 + tr_z_xxz_xxyyz[i] * pa_x[i];

        tr_z_xxxz_xxyzz[i] = 2.0 * tr_z_xz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxz_xyzz[i] * fe_0 + tr_z_xxz_xxyzz[i] * pa_x[i];

        tr_z_xxxz_xxzzz[i] = 2.0 * tr_z_xz_xxzzz[i] * fe_0 + 2.0 * tr_z_xxz_xzzz[i] * fe_0 + tr_z_xxz_xxzzz[i] * pa_x[i];

        tr_z_xxxz_xyyyy[i] = ts_xxx_xyyyy[i] * fe_0 + tr_z_xxx_xyyyy[i] * pa_z[i];

        tr_z_xxxz_xyyyz[i] = 2.0 * tr_z_xz_xyyyz[i] * fe_0 + tr_z_xxz_yyyz[i] * fe_0 + tr_z_xxz_xyyyz[i] * pa_x[i];

        tr_z_xxxz_xyyzz[i] = 2.0 * tr_z_xz_xyyzz[i] * fe_0 + tr_z_xxz_yyzz[i] * fe_0 + tr_z_xxz_xyyzz[i] * pa_x[i];

        tr_z_xxxz_xyzzz[i] = 2.0 * tr_z_xz_xyzzz[i] * fe_0 + tr_z_xxz_yzzz[i] * fe_0 + tr_z_xxz_xyzzz[i] * pa_x[i];

        tr_z_xxxz_xzzzz[i] = 2.0 * tr_z_xz_xzzzz[i] * fe_0 + tr_z_xxz_zzzz[i] * fe_0 + tr_z_xxz_xzzzz[i] * pa_x[i];

        tr_z_xxxz_yyyyy[i] = 2.0 * tr_z_xz_yyyyy[i] * fe_0 + tr_z_xxz_yyyyy[i] * pa_x[i];

        tr_z_xxxz_yyyyz[i] = 2.0 * tr_z_xz_yyyyz[i] * fe_0 + tr_z_xxz_yyyyz[i] * pa_x[i];

        tr_z_xxxz_yyyzz[i] = 2.0 * tr_z_xz_yyyzz[i] * fe_0 + tr_z_xxz_yyyzz[i] * pa_x[i];

        tr_z_xxxz_yyzzz[i] = 2.0 * tr_z_xz_yyzzz[i] * fe_0 + tr_z_xxz_yyzzz[i] * pa_x[i];

        tr_z_xxxz_yzzzz[i] = 2.0 * tr_z_xz_yzzzz[i] * fe_0 + tr_z_xxz_yzzzz[i] * pa_x[i];

        tr_z_xxxz_zzzzz[i] = 2.0 * tr_z_xz_zzzzz[i] * fe_0 + tr_z_xxz_zzzzz[i] * pa_x[i];
    }

    // Set up 693-714 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xx_xxxxx, tr_z_xx_xxxxz, tr_z_xx_xxxzz, tr_z_xx_xxzzz, tr_z_xx_xzzzz, tr_z_xxy_xxxxx, tr_z_xxy_xxxxz, tr_z_xxy_xxxzz, tr_z_xxy_xxzzz, tr_z_xxy_xzzzz, tr_z_xxyy_xxxxx, tr_z_xxyy_xxxxy, tr_z_xxyy_xxxxz, tr_z_xxyy_xxxyy, tr_z_xxyy_xxxyz, tr_z_xxyy_xxxzz, tr_z_xxyy_xxyyy, tr_z_xxyy_xxyyz, tr_z_xxyy_xxyzz, tr_z_xxyy_xxzzz, tr_z_xxyy_xyyyy, tr_z_xxyy_xyyyz, tr_z_xxyy_xyyzz, tr_z_xxyy_xyzzz, tr_z_xxyy_xzzzz, tr_z_xxyy_yyyyy, tr_z_xxyy_yyyyz, tr_z_xxyy_yyyzz, tr_z_xxyy_yyzzz, tr_z_xxyy_yzzzz, tr_z_xxyy_zzzzz, tr_z_xyy_xxxxy, tr_z_xyy_xxxy, tr_z_xyy_xxxyy, tr_z_xyy_xxxyz, tr_z_xyy_xxyy, tr_z_xyy_xxyyy, tr_z_xyy_xxyyz, tr_z_xyy_xxyz, tr_z_xyy_xxyzz, tr_z_xyy_xyyy, tr_z_xyy_xyyyy, tr_z_xyy_xyyyz, tr_z_xyy_xyyz, tr_z_xyy_xyyzz, tr_z_xyy_xyzz, tr_z_xyy_xyzzz, tr_z_xyy_yyyy, tr_z_xyy_yyyyy, tr_z_xyy_yyyyz, tr_z_xyy_yyyz, tr_z_xyy_yyyzz, tr_z_xyy_yyzz, tr_z_xyy_yyzzz, tr_z_xyy_yzzz, tr_z_xyy_yzzzz, tr_z_xyy_zzzzz, tr_z_yy_xxxxy, tr_z_yy_xxxyy, tr_z_yy_xxxyz, tr_z_yy_xxyyy, tr_z_yy_xxyyz, tr_z_yy_xxyzz, tr_z_yy_xyyyy, tr_z_yy_xyyyz, tr_z_yy_xyyzz, tr_z_yy_xyzzz, tr_z_yy_yyyyy, tr_z_yy_yyyyz, tr_z_yy_yyyzz, tr_z_yy_yyzzz, tr_z_yy_yzzzz, tr_z_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyy_xxxxx[i] = tr_z_xx_xxxxx[i] * fe_0 + tr_z_xxy_xxxxx[i] * pa_y[i];

        tr_z_xxyy_xxxxy[i] = tr_z_yy_xxxxy[i] * fe_0 + 4.0 * tr_z_xyy_xxxy[i] * fe_0 + tr_z_xyy_xxxxy[i] * pa_x[i];

        tr_z_xxyy_xxxxz[i] = tr_z_xx_xxxxz[i] * fe_0 + tr_z_xxy_xxxxz[i] * pa_y[i];

        tr_z_xxyy_xxxyy[i] = tr_z_yy_xxxyy[i] * fe_0 + 3.0 * tr_z_xyy_xxyy[i] * fe_0 + tr_z_xyy_xxxyy[i] * pa_x[i];

        tr_z_xxyy_xxxyz[i] = tr_z_yy_xxxyz[i] * fe_0 + 3.0 * tr_z_xyy_xxyz[i] * fe_0 + tr_z_xyy_xxxyz[i] * pa_x[i];

        tr_z_xxyy_xxxzz[i] = tr_z_xx_xxxzz[i] * fe_0 + tr_z_xxy_xxxzz[i] * pa_y[i];

        tr_z_xxyy_xxyyy[i] = tr_z_yy_xxyyy[i] * fe_0 + 2.0 * tr_z_xyy_xyyy[i] * fe_0 + tr_z_xyy_xxyyy[i] * pa_x[i];

        tr_z_xxyy_xxyyz[i] = tr_z_yy_xxyyz[i] * fe_0 + 2.0 * tr_z_xyy_xyyz[i] * fe_0 + tr_z_xyy_xxyyz[i] * pa_x[i];

        tr_z_xxyy_xxyzz[i] = tr_z_yy_xxyzz[i] * fe_0 + 2.0 * tr_z_xyy_xyzz[i] * fe_0 + tr_z_xyy_xxyzz[i] * pa_x[i];

        tr_z_xxyy_xxzzz[i] = tr_z_xx_xxzzz[i] * fe_0 + tr_z_xxy_xxzzz[i] * pa_y[i];

        tr_z_xxyy_xyyyy[i] = tr_z_yy_xyyyy[i] * fe_0 + tr_z_xyy_yyyy[i] * fe_0 + tr_z_xyy_xyyyy[i] * pa_x[i];

        tr_z_xxyy_xyyyz[i] = tr_z_yy_xyyyz[i] * fe_0 + tr_z_xyy_yyyz[i] * fe_0 + tr_z_xyy_xyyyz[i] * pa_x[i];

        tr_z_xxyy_xyyzz[i] = tr_z_yy_xyyzz[i] * fe_0 + tr_z_xyy_yyzz[i] * fe_0 + tr_z_xyy_xyyzz[i] * pa_x[i];

        tr_z_xxyy_xyzzz[i] = tr_z_yy_xyzzz[i] * fe_0 + tr_z_xyy_yzzz[i] * fe_0 + tr_z_xyy_xyzzz[i] * pa_x[i];

        tr_z_xxyy_xzzzz[i] = tr_z_xx_xzzzz[i] * fe_0 + tr_z_xxy_xzzzz[i] * pa_y[i];

        tr_z_xxyy_yyyyy[i] = tr_z_yy_yyyyy[i] * fe_0 + tr_z_xyy_yyyyy[i] * pa_x[i];

        tr_z_xxyy_yyyyz[i] = tr_z_yy_yyyyz[i] * fe_0 + tr_z_xyy_yyyyz[i] * pa_x[i];

        tr_z_xxyy_yyyzz[i] = tr_z_yy_yyyzz[i] * fe_0 + tr_z_xyy_yyyzz[i] * pa_x[i];

        tr_z_xxyy_yyzzz[i] = tr_z_yy_yyzzz[i] * fe_0 + tr_z_xyy_yyzzz[i] * pa_x[i];

        tr_z_xxyy_yzzzz[i] = tr_z_yy_yzzzz[i] * fe_0 + tr_z_xyy_yzzzz[i] * pa_x[i];

        tr_z_xxyy_zzzzz[i] = tr_z_yy_zzzzz[i] * fe_0 + tr_z_xyy_zzzzz[i] * pa_x[i];
    }

    // Set up 714-735 components of targeted buffer : GH

    auto tr_z_xxyz_xxxxx = pbuffer.data(idx_dip_gh + 714);

    auto tr_z_xxyz_xxxxy = pbuffer.data(idx_dip_gh + 715);

    auto tr_z_xxyz_xxxxz = pbuffer.data(idx_dip_gh + 716);

    auto tr_z_xxyz_xxxyy = pbuffer.data(idx_dip_gh + 717);

    auto tr_z_xxyz_xxxyz = pbuffer.data(idx_dip_gh + 718);

    auto tr_z_xxyz_xxxzz = pbuffer.data(idx_dip_gh + 719);

    auto tr_z_xxyz_xxyyy = pbuffer.data(idx_dip_gh + 720);

    auto tr_z_xxyz_xxyyz = pbuffer.data(idx_dip_gh + 721);

    auto tr_z_xxyz_xxyzz = pbuffer.data(idx_dip_gh + 722);

    auto tr_z_xxyz_xxzzz = pbuffer.data(idx_dip_gh + 723);

    auto tr_z_xxyz_xyyyy = pbuffer.data(idx_dip_gh + 724);

    auto tr_z_xxyz_xyyyz = pbuffer.data(idx_dip_gh + 725);

    auto tr_z_xxyz_xyyzz = pbuffer.data(idx_dip_gh + 726);

    auto tr_z_xxyz_xyzzz = pbuffer.data(idx_dip_gh + 727);

    auto tr_z_xxyz_xzzzz = pbuffer.data(idx_dip_gh + 728);

    auto tr_z_xxyz_yyyyy = pbuffer.data(idx_dip_gh + 729);

    auto tr_z_xxyz_yyyyz = pbuffer.data(idx_dip_gh + 730);

    auto tr_z_xxyz_yyyzz = pbuffer.data(idx_dip_gh + 731);

    auto tr_z_xxyz_yyzzz = pbuffer.data(idx_dip_gh + 732);

    auto tr_z_xxyz_yzzzz = pbuffer.data(idx_dip_gh + 733);

    auto tr_z_xxyz_zzzzz = pbuffer.data(idx_dip_gh + 734);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyz_xxxxx, tr_z_xxyz_xxxxy, tr_z_xxyz_xxxxz, tr_z_xxyz_xxxyy, tr_z_xxyz_xxxyz, tr_z_xxyz_xxxzz, tr_z_xxyz_xxyyy, tr_z_xxyz_xxyyz, tr_z_xxyz_xxyzz, tr_z_xxyz_xxzzz, tr_z_xxyz_xyyyy, tr_z_xxyz_xyyyz, tr_z_xxyz_xyyzz, tr_z_xxyz_xyzzz, tr_z_xxyz_xzzzz, tr_z_xxyz_yyyyy, tr_z_xxyz_yyyyz, tr_z_xxyz_yyyzz, tr_z_xxyz_yyzzz, tr_z_xxyz_yzzzz, tr_z_xxyz_zzzzz, tr_z_xxz_xxxx, tr_z_xxz_xxxxx, tr_z_xxz_xxxxy, tr_z_xxz_xxxxz, tr_z_xxz_xxxy, tr_z_xxz_xxxyy, tr_z_xxz_xxxyz, tr_z_xxz_xxxz, tr_z_xxz_xxxzz, tr_z_xxz_xxyy, tr_z_xxz_xxyyy, tr_z_xxz_xxyyz, tr_z_xxz_xxyz, tr_z_xxz_xxyzz, tr_z_xxz_xxzz, tr_z_xxz_xxzzz, tr_z_xxz_xyyy, tr_z_xxz_xyyyy, tr_z_xxz_xyyyz, tr_z_xxz_xyyz, tr_z_xxz_xyyzz, tr_z_xxz_xyzz, tr_z_xxz_xyzzz, tr_z_xxz_xzzz, tr_z_xxz_xzzzz, tr_z_xxz_zzzzz, tr_z_xyz_yyyyy, tr_z_xyz_yyyyz, tr_z_xyz_yyyzz, tr_z_xyz_yyzzz, tr_z_xyz_yzzzz, tr_z_yz_yyyyy, tr_z_yz_yyyyz, tr_z_yz_yyyzz, tr_z_yz_yyzzz, tr_z_yz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyz_xxxxx[i] = tr_z_xxz_xxxxx[i] * pa_y[i];

        tr_z_xxyz_xxxxy[i] = tr_z_xxz_xxxx[i] * fe_0 + tr_z_xxz_xxxxy[i] * pa_y[i];

        tr_z_xxyz_xxxxz[i] = tr_z_xxz_xxxxz[i] * pa_y[i];

        tr_z_xxyz_xxxyy[i] = 2.0 * tr_z_xxz_xxxy[i] * fe_0 + tr_z_xxz_xxxyy[i] * pa_y[i];

        tr_z_xxyz_xxxyz[i] = tr_z_xxz_xxxz[i] * fe_0 + tr_z_xxz_xxxyz[i] * pa_y[i];

        tr_z_xxyz_xxxzz[i] = tr_z_xxz_xxxzz[i] * pa_y[i];

        tr_z_xxyz_xxyyy[i] = 3.0 * tr_z_xxz_xxyy[i] * fe_0 + tr_z_xxz_xxyyy[i] * pa_y[i];

        tr_z_xxyz_xxyyz[i] = 2.0 * tr_z_xxz_xxyz[i] * fe_0 + tr_z_xxz_xxyyz[i] * pa_y[i];

        tr_z_xxyz_xxyzz[i] = tr_z_xxz_xxzz[i] * fe_0 + tr_z_xxz_xxyzz[i] * pa_y[i];

        tr_z_xxyz_xxzzz[i] = tr_z_xxz_xxzzz[i] * pa_y[i];

        tr_z_xxyz_xyyyy[i] = 4.0 * tr_z_xxz_xyyy[i] * fe_0 + tr_z_xxz_xyyyy[i] * pa_y[i];

        tr_z_xxyz_xyyyz[i] = 3.0 * tr_z_xxz_xyyz[i] * fe_0 + tr_z_xxz_xyyyz[i] * pa_y[i];

        tr_z_xxyz_xyyzz[i] = 2.0 * tr_z_xxz_xyzz[i] * fe_0 + tr_z_xxz_xyyzz[i] * pa_y[i];

        tr_z_xxyz_xyzzz[i] = tr_z_xxz_xzzz[i] * fe_0 + tr_z_xxz_xyzzz[i] * pa_y[i];

        tr_z_xxyz_xzzzz[i] = tr_z_xxz_xzzzz[i] * pa_y[i];

        tr_z_xxyz_yyyyy[i] = tr_z_yz_yyyyy[i] * fe_0 + tr_z_xyz_yyyyy[i] * pa_x[i];

        tr_z_xxyz_yyyyz[i] = tr_z_yz_yyyyz[i] * fe_0 + tr_z_xyz_yyyyz[i] * pa_x[i];

        tr_z_xxyz_yyyzz[i] = tr_z_yz_yyyzz[i] * fe_0 + tr_z_xyz_yyyzz[i] * pa_x[i];

        tr_z_xxyz_yyzzz[i] = tr_z_yz_yyzzz[i] * fe_0 + tr_z_xyz_yyzzz[i] * pa_x[i];

        tr_z_xxyz_yzzzz[i] = tr_z_yz_yzzzz[i] * fe_0 + tr_z_xyz_yzzzz[i] * pa_x[i];

        tr_z_xxyz_zzzzz[i] = tr_z_xxz_zzzzz[i] * pa_y[i];
    }

    // Set up 735-756 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_z_xxzz_xxxxx, tr_z_xxzz_xxxxy, tr_z_xxzz_xxxxz, tr_z_xxzz_xxxyy, tr_z_xxzz_xxxyz, tr_z_xxzz_xxxzz, tr_z_xxzz_xxyyy, tr_z_xxzz_xxyyz, tr_z_xxzz_xxyzz, tr_z_xxzz_xxzzz, tr_z_xxzz_xyyyy, tr_z_xxzz_xyyyz, tr_z_xxzz_xyyzz, tr_z_xxzz_xyzzz, tr_z_xxzz_xzzzz, tr_z_xxzz_yyyyy, tr_z_xxzz_yyyyz, tr_z_xxzz_yyyzz, tr_z_xxzz_yyzzz, tr_z_xxzz_yzzzz, tr_z_xxzz_zzzzz, tr_z_xzz_xxxx, tr_z_xzz_xxxxx, tr_z_xzz_xxxxy, tr_z_xzz_xxxxz, tr_z_xzz_xxxy, tr_z_xzz_xxxyy, tr_z_xzz_xxxyz, tr_z_xzz_xxxz, tr_z_xzz_xxxzz, tr_z_xzz_xxyy, tr_z_xzz_xxyyy, tr_z_xzz_xxyyz, tr_z_xzz_xxyz, tr_z_xzz_xxyzz, tr_z_xzz_xxzz, tr_z_xzz_xxzzz, tr_z_xzz_xyyy, tr_z_xzz_xyyyy, tr_z_xzz_xyyyz, tr_z_xzz_xyyz, tr_z_xzz_xyyzz, tr_z_xzz_xyzz, tr_z_xzz_xyzzz, tr_z_xzz_xzzz, tr_z_xzz_xzzzz, tr_z_xzz_yyyy, tr_z_xzz_yyyyy, tr_z_xzz_yyyyz, tr_z_xzz_yyyz, tr_z_xzz_yyyzz, tr_z_xzz_yyzz, tr_z_xzz_yyzzz, tr_z_xzz_yzzz, tr_z_xzz_yzzzz, tr_z_xzz_zzzz, tr_z_xzz_zzzzz, tr_z_zz_xxxxx, tr_z_zz_xxxxy, tr_z_zz_xxxxz, tr_z_zz_xxxyy, tr_z_zz_xxxyz, tr_z_zz_xxxzz, tr_z_zz_xxyyy, tr_z_zz_xxyyz, tr_z_zz_xxyzz, tr_z_zz_xxzzz, tr_z_zz_xyyyy, tr_z_zz_xyyyz, tr_z_zz_xyyzz, tr_z_zz_xyzzz, tr_z_zz_xzzzz, tr_z_zz_yyyyy, tr_z_zz_yyyyz, tr_z_zz_yyyzz, tr_z_zz_yyzzz, tr_z_zz_yzzzz, tr_z_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzz_xxxxx[i] = tr_z_zz_xxxxx[i] * fe_0 + 5.0 * tr_z_xzz_xxxx[i] * fe_0 + tr_z_xzz_xxxxx[i] * pa_x[i];

        tr_z_xxzz_xxxxy[i] = tr_z_zz_xxxxy[i] * fe_0 + 4.0 * tr_z_xzz_xxxy[i] * fe_0 + tr_z_xzz_xxxxy[i] * pa_x[i];

        tr_z_xxzz_xxxxz[i] = tr_z_zz_xxxxz[i] * fe_0 + 4.0 * tr_z_xzz_xxxz[i] * fe_0 + tr_z_xzz_xxxxz[i] * pa_x[i];

        tr_z_xxzz_xxxyy[i] = tr_z_zz_xxxyy[i] * fe_0 + 3.0 * tr_z_xzz_xxyy[i] * fe_0 + tr_z_xzz_xxxyy[i] * pa_x[i];

        tr_z_xxzz_xxxyz[i] = tr_z_zz_xxxyz[i] * fe_0 + 3.0 * tr_z_xzz_xxyz[i] * fe_0 + tr_z_xzz_xxxyz[i] * pa_x[i];

        tr_z_xxzz_xxxzz[i] = tr_z_zz_xxxzz[i] * fe_0 + 3.0 * tr_z_xzz_xxzz[i] * fe_0 + tr_z_xzz_xxxzz[i] * pa_x[i];

        tr_z_xxzz_xxyyy[i] = tr_z_zz_xxyyy[i] * fe_0 + 2.0 * tr_z_xzz_xyyy[i] * fe_0 + tr_z_xzz_xxyyy[i] * pa_x[i];

        tr_z_xxzz_xxyyz[i] = tr_z_zz_xxyyz[i] * fe_0 + 2.0 * tr_z_xzz_xyyz[i] * fe_0 + tr_z_xzz_xxyyz[i] * pa_x[i];

        tr_z_xxzz_xxyzz[i] = tr_z_zz_xxyzz[i] * fe_0 + 2.0 * tr_z_xzz_xyzz[i] * fe_0 + tr_z_xzz_xxyzz[i] * pa_x[i];

        tr_z_xxzz_xxzzz[i] = tr_z_zz_xxzzz[i] * fe_0 + 2.0 * tr_z_xzz_xzzz[i] * fe_0 + tr_z_xzz_xxzzz[i] * pa_x[i];

        tr_z_xxzz_xyyyy[i] = tr_z_zz_xyyyy[i] * fe_0 + tr_z_xzz_yyyy[i] * fe_0 + tr_z_xzz_xyyyy[i] * pa_x[i];

        tr_z_xxzz_xyyyz[i] = tr_z_zz_xyyyz[i] * fe_0 + tr_z_xzz_yyyz[i] * fe_0 + tr_z_xzz_xyyyz[i] * pa_x[i];

        tr_z_xxzz_xyyzz[i] = tr_z_zz_xyyzz[i] * fe_0 + tr_z_xzz_yyzz[i] * fe_0 + tr_z_xzz_xyyzz[i] * pa_x[i];

        tr_z_xxzz_xyzzz[i] = tr_z_zz_xyzzz[i] * fe_0 + tr_z_xzz_yzzz[i] * fe_0 + tr_z_xzz_xyzzz[i] * pa_x[i];

        tr_z_xxzz_xzzzz[i] = tr_z_zz_xzzzz[i] * fe_0 + tr_z_xzz_zzzz[i] * fe_0 + tr_z_xzz_xzzzz[i] * pa_x[i];

        tr_z_xxzz_yyyyy[i] = tr_z_zz_yyyyy[i] * fe_0 + tr_z_xzz_yyyyy[i] * pa_x[i];

        tr_z_xxzz_yyyyz[i] = tr_z_zz_yyyyz[i] * fe_0 + tr_z_xzz_yyyyz[i] * pa_x[i];

        tr_z_xxzz_yyyzz[i] = tr_z_zz_yyyzz[i] * fe_0 + tr_z_xzz_yyyzz[i] * pa_x[i];

        tr_z_xxzz_yyzzz[i] = tr_z_zz_yyzzz[i] * fe_0 + tr_z_xzz_yyzzz[i] * pa_x[i];

        tr_z_xxzz_yzzzz[i] = tr_z_zz_yzzzz[i] * fe_0 + tr_z_xzz_yzzzz[i] * pa_x[i];

        tr_z_xxzz_zzzzz[i] = tr_z_zz_zzzzz[i] * fe_0 + tr_z_xzz_zzzzz[i] * pa_x[i];
    }

    // Set up 756-777 components of targeted buffer : GH

    auto tr_z_xyyy_xxxxx = pbuffer.data(idx_dip_gh + 756);

    auto tr_z_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 757);

    auto tr_z_xyyy_xxxxz = pbuffer.data(idx_dip_gh + 758);

    auto tr_z_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 759);

    auto tr_z_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 760);

    auto tr_z_xyyy_xxxzz = pbuffer.data(idx_dip_gh + 761);

    auto tr_z_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 762);

    auto tr_z_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 763);

    auto tr_z_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 764);

    auto tr_z_xyyy_xxzzz = pbuffer.data(idx_dip_gh + 765);

    auto tr_z_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 766);

    auto tr_z_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 767);

    auto tr_z_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 768);

    auto tr_z_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 769);

    auto tr_z_xyyy_xzzzz = pbuffer.data(idx_dip_gh + 770);

    auto tr_z_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 771);

    auto tr_z_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 772);

    auto tr_z_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 773);

    auto tr_z_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 774);

    auto tr_z_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 775);

    auto tr_z_xyyy_zzzzz = pbuffer.data(idx_dip_gh + 776);

    #pragma omp simd aligned(pa_x, tr_z_xyyy_xxxxx, tr_z_xyyy_xxxxy, tr_z_xyyy_xxxxz, tr_z_xyyy_xxxyy, tr_z_xyyy_xxxyz, tr_z_xyyy_xxxzz, tr_z_xyyy_xxyyy, tr_z_xyyy_xxyyz, tr_z_xyyy_xxyzz, tr_z_xyyy_xxzzz, tr_z_xyyy_xyyyy, tr_z_xyyy_xyyyz, tr_z_xyyy_xyyzz, tr_z_xyyy_xyzzz, tr_z_xyyy_xzzzz, tr_z_xyyy_yyyyy, tr_z_xyyy_yyyyz, tr_z_xyyy_yyyzz, tr_z_xyyy_yyzzz, tr_z_xyyy_yzzzz, tr_z_xyyy_zzzzz, tr_z_yyy_xxxx, tr_z_yyy_xxxxx, tr_z_yyy_xxxxy, tr_z_yyy_xxxxz, tr_z_yyy_xxxy, tr_z_yyy_xxxyy, tr_z_yyy_xxxyz, tr_z_yyy_xxxz, tr_z_yyy_xxxzz, tr_z_yyy_xxyy, tr_z_yyy_xxyyy, tr_z_yyy_xxyyz, tr_z_yyy_xxyz, tr_z_yyy_xxyzz, tr_z_yyy_xxzz, tr_z_yyy_xxzzz, tr_z_yyy_xyyy, tr_z_yyy_xyyyy, tr_z_yyy_xyyyz, tr_z_yyy_xyyz, tr_z_yyy_xyyzz, tr_z_yyy_xyzz, tr_z_yyy_xyzzz, tr_z_yyy_xzzz, tr_z_yyy_xzzzz, tr_z_yyy_yyyy, tr_z_yyy_yyyyy, tr_z_yyy_yyyyz, tr_z_yyy_yyyz, tr_z_yyy_yyyzz, tr_z_yyy_yyzz, tr_z_yyy_yyzzz, tr_z_yyy_yzzz, tr_z_yyy_yzzzz, tr_z_yyy_zzzz, tr_z_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyy_xxxxx[i] = 5.0 * tr_z_yyy_xxxx[i] * fe_0 + tr_z_yyy_xxxxx[i] * pa_x[i];

        tr_z_xyyy_xxxxy[i] = 4.0 * tr_z_yyy_xxxy[i] * fe_0 + tr_z_yyy_xxxxy[i] * pa_x[i];

        tr_z_xyyy_xxxxz[i] = 4.0 * tr_z_yyy_xxxz[i] * fe_0 + tr_z_yyy_xxxxz[i] * pa_x[i];

        tr_z_xyyy_xxxyy[i] = 3.0 * tr_z_yyy_xxyy[i] * fe_0 + tr_z_yyy_xxxyy[i] * pa_x[i];

        tr_z_xyyy_xxxyz[i] = 3.0 * tr_z_yyy_xxyz[i] * fe_0 + tr_z_yyy_xxxyz[i] * pa_x[i];

        tr_z_xyyy_xxxzz[i] = 3.0 * tr_z_yyy_xxzz[i] * fe_0 + tr_z_yyy_xxxzz[i] * pa_x[i];

        tr_z_xyyy_xxyyy[i] = 2.0 * tr_z_yyy_xyyy[i] * fe_0 + tr_z_yyy_xxyyy[i] * pa_x[i];

        tr_z_xyyy_xxyyz[i] = 2.0 * tr_z_yyy_xyyz[i] * fe_0 + tr_z_yyy_xxyyz[i] * pa_x[i];

        tr_z_xyyy_xxyzz[i] = 2.0 * tr_z_yyy_xyzz[i] * fe_0 + tr_z_yyy_xxyzz[i] * pa_x[i];

        tr_z_xyyy_xxzzz[i] = 2.0 * tr_z_yyy_xzzz[i] * fe_0 + tr_z_yyy_xxzzz[i] * pa_x[i];

        tr_z_xyyy_xyyyy[i] = tr_z_yyy_yyyy[i] * fe_0 + tr_z_yyy_xyyyy[i] * pa_x[i];

        tr_z_xyyy_xyyyz[i] = tr_z_yyy_yyyz[i] * fe_0 + tr_z_yyy_xyyyz[i] * pa_x[i];

        tr_z_xyyy_xyyzz[i] = tr_z_yyy_yyzz[i] * fe_0 + tr_z_yyy_xyyzz[i] * pa_x[i];

        tr_z_xyyy_xyzzz[i] = tr_z_yyy_yzzz[i] * fe_0 + tr_z_yyy_xyzzz[i] * pa_x[i];

        tr_z_xyyy_xzzzz[i] = tr_z_yyy_zzzz[i] * fe_0 + tr_z_yyy_xzzzz[i] * pa_x[i];

        tr_z_xyyy_yyyyy[i] = tr_z_yyy_yyyyy[i] * pa_x[i];

        tr_z_xyyy_yyyyz[i] = tr_z_yyy_yyyyz[i] * pa_x[i];

        tr_z_xyyy_yyyzz[i] = tr_z_yyy_yyyzz[i] * pa_x[i];

        tr_z_xyyy_yyzzz[i] = tr_z_yyy_yyzzz[i] * pa_x[i];

        tr_z_xyyy_yzzzz[i] = tr_z_yyy_yzzzz[i] * pa_x[i];

        tr_z_xyyy_zzzzz[i] = tr_z_yyy_zzzzz[i] * pa_x[i];
    }

    // Set up 777-798 components of targeted buffer : GH

    auto tr_z_xyyz_xxxxx = pbuffer.data(idx_dip_gh + 777);

    auto tr_z_xyyz_xxxxy = pbuffer.data(idx_dip_gh + 778);

    auto tr_z_xyyz_xxxxz = pbuffer.data(idx_dip_gh + 779);

    auto tr_z_xyyz_xxxyy = pbuffer.data(idx_dip_gh + 780);

    auto tr_z_xyyz_xxxyz = pbuffer.data(idx_dip_gh + 781);

    auto tr_z_xyyz_xxxzz = pbuffer.data(idx_dip_gh + 782);

    auto tr_z_xyyz_xxyyy = pbuffer.data(idx_dip_gh + 783);

    auto tr_z_xyyz_xxyyz = pbuffer.data(idx_dip_gh + 784);

    auto tr_z_xyyz_xxyzz = pbuffer.data(idx_dip_gh + 785);

    auto tr_z_xyyz_xxzzz = pbuffer.data(idx_dip_gh + 786);

    auto tr_z_xyyz_xyyyy = pbuffer.data(idx_dip_gh + 787);

    auto tr_z_xyyz_xyyyz = pbuffer.data(idx_dip_gh + 788);

    auto tr_z_xyyz_xyyzz = pbuffer.data(idx_dip_gh + 789);

    auto tr_z_xyyz_xyzzz = pbuffer.data(idx_dip_gh + 790);

    auto tr_z_xyyz_xzzzz = pbuffer.data(idx_dip_gh + 791);

    auto tr_z_xyyz_yyyyy = pbuffer.data(idx_dip_gh + 792);

    auto tr_z_xyyz_yyyyz = pbuffer.data(idx_dip_gh + 793);

    auto tr_z_xyyz_yyyzz = pbuffer.data(idx_dip_gh + 794);

    auto tr_z_xyyz_yyzzz = pbuffer.data(idx_dip_gh + 795);

    auto tr_z_xyyz_yzzzz = pbuffer.data(idx_dip_gh + 796);

    auto tr_z_xyyz_zzzzz = pbuffer.data(idx_dip_gh + 797);

    #pragma omp simd aligned(pa_x, tr_z_xyyz_xxxxx, tr_z_xyyz_xxxxy, tr_z_xyyz_xxxxz, tr_z_xyyz_xxxyy, tr_z_xyyz_xxxyz, tr_z_xyyz_xxxzz, tr_z_xyyz_xxyyy, tr_z_xyyz_xxyyz, tr_z_xyyz_xxyzz, tr_z_xyyz_xxzzz, tr_z_xyyz_xyyyy, tr_z_xyyz_xyyyz, tr_z_xyyz_xyyzz, tr_z_xyyz_xyzzz, tr_z_xyyz_xzzzz, tr_z_xyyz_yyyyy, tr_z_xyyz_yyyyz, tr_z_xyyz_yyyzz, tr_z_xyyz_yyzzz, tr_z_xyyz_yzzzz, tr_z_xyyz_zzzzz, tr_z_yyz_xxxx, tr_z_yyz_xxxxx, tr_z_yyz_xxxxy, tr_z_yyz_xxxxz, tr_z_yyz_xxxy, tr_z_yyz_xxxyy, tr_z_yyz_xxxyz, tr_z_yyz_xxxz, tr_z_yyz_xxxzz, tr_z_yyz_xxyy, tr_z_yyz_xxyyy, tr_z_yyz_xxyyz, tr_z_yyz_xxyz, tr_z_yyz_xxyzz, tr_z_yyz_xxzz, tr_z_yyz_xxzzz, tr_z_yyz_xyyy, tr_z_yyz_xyyyy, tr_z_yyz_xyyyz, tr_z_yyz_xyyz, tr_z_yyz_xyyzz, tr_z_yyz_xyzz, tr_z_yyz_xyzzz, tr_z_yyz_xzzz, tr_z_yyz_xzzzz, tr_z_yyz_yyyy, tr_z_yyz_yyyyy, tr_z_yyz_yyyyz, tr_z_yyz_yyyz, tr_z_yyz_yyyzz, tr_z_yyz_yyzz, tr_z_yyz_yyzzz, tr_z_yyz_yzzz, tr_z_yyz_yzzzz, tr_z_yyz_zzzz, tr_z_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyz_xxxxx[i] = 5.0 * tr_z_yyz_xxxx[i] * fe_0 + tr_z_yyz_xxxxx[i] * pa_x[i];

        tr_z_xyyz_xxxxy[i] = 4.0 * tr_z_yyz_xxxy[i] * fe_0 + tr_z_yyz_xxxxy[i] * pa_x[i];

        tr_z_xyyz_xxxxz[i] = 4.0 * tr_z_yyz_xxxz[i] * fe_0 + tr_z_yyz_xxxxz[i] * pa_x[i];

        tr_z_xyyz_xxxyy[i] = 3.0 * tr_z_yyz_xxyy[i] * fe_0 + tr_z_yyz_xxxyy[i] * pa_x[i];

        tr_z_xyyz_xxxyz[i] = 3.0 * tr_z_yyz_xxyz[i] * fe_0 + tr_z_yyz_xxxyz[i] * pa_x[i];

        tr_z_xyyz_xxxzz[i] = 3.0 * tr_z_yyz_xxzz[i] * fe_0 + tr_z_yyz_xxxzz[i] * pa_x[i];

        tr_z_xyyz_xxyyy[i] = 2.0 * tr_z_yyz_xyyy[i] * fe_0 + tr_z_yyz_xxyyy[i] * pa_x[i];

        tr_z_xyyz_xxyyz[i] = 2.0 * tr_z_yyz_xyyz[i] * fe_0 + tr_z_yyz_xxyyz[i] * pa_x[i];

        tr_z_xyyz_xxyzz[i] = 2.0 * tr_z_yyz_xyzz[i] * fe_0 + tr_z_yyz_xxyzz[i] * pa_x[i];

        tr_z_xyyz_xxzzz[i] = 2.0 * tr_z_yyz_xzzz[i] * fe_0 + tr_z_yyz_xxzzz[i] * pa_x[i];

        tr_z_xyyz_xyyyy[i] = tr_z_yyz_yyyy[i] * fe_0 + tr_z_yyz_xyyyy[i] * pa_x[i];

        tr_z_xyyz_xyyyz[i] = tr_z_yyz_yyyz[i] * fe_0 + tr_z_yyz_xyyyz[i] * pa_x[i];

        tr_z_xyyz_xyyzz[i] = tr_z_yyz_yyzz[i] * fe_0 + tr_z_yyz_xyyzz[i] * pa_x[i];

        tr_z_xyyz_xyzzz[i] = tr_z_yyz_yzzz[i] * fe_0 + tr_z_yyz_xyzzz[i] * pa_x[i];

        tr_z_xyyz_xzzzz[i] = tr_z_yyz_zzzz[i] * fe_0 + tr_z_yyz_xzzzz[i] * pa_x[i];

        tr_z_xyyz_yyyyy[i] = tr_z_yyz_yyyyy[i] * pa_x[i];

        tr_z_xyyz_yyyyz[i] = tr_z_yyz_yyyyz[i] * pa_x[i];

        tr_z_xyyz_yyyzz[i] = tr_z_yyz_yyyzz[i] * pa_x[i];

        tr_z_xyyz_yyzzz[i] = tr_z_yyz_yyzzz[i] * pa_x[i];

        tr_z_xyyz_yzzzz[i] = tr_z_yyz_yzzzz[i] * pa_x[i];

        tr_z_xyyz_zzzzz[i] = tr_z_yyz_zzzzz[i] * pa_x[i];
    }

    // Set up 798-819 components of targeted buffer : GH

    auto tr_z_xyzz_xxxxx = pbuffer.data(idx_dip_gh + 798);

    auto tr_z_xyzz_xxxxy = pbuffer.data(idx_dip_gh + 799);

    auto tr_z_xyzz_xxxxz = pbuffer.data(idx_dip_gh + 800);

    auto tr_z_xyzz_xxxyy = pbuffer.data(idx_dip_gh + 801);

    auto tr_z_xyzz_xxxyz = pbuffer.data(idx_dip_gh + 802);

    auto tr_z_xyzz_xxxzz = pbuffer.data(idx_dip_gh + 803);

    auto tr_z_xyzz_xxyyy = pbuffer.data(idx_dip_gh + 804);

    auto tr_z_xyzz_xxyyz = pbuffer.data(idx_dip_gh + 805);

    auto tr_z_xyzz_xxyzz = pbuffer.data(idx_dip_gh + 806);

    auto tr_z_xyzz_xxzzz = pbuffer.data(idx_dip_gh + 807);

    auto tr_z_xyzz_xyyyy = pbuffer.data(idx_dip_gh + 808);

    auto tr_z_xyzz_xyyyz = pbuffer.data(idx_dip_gh + 809);

    auto tr_z_xyzz_xyyzz = pbuffer.data(idx_dip_gh + 810);

    auto tr_z_xyzz_xyzzz = pbuffer.data(idx_dip_gh + 811);

    auto tr_z_xyzz_xzzzz = pbuffer.data(idx_dip_gh + 812);

    auto tr_z_xyzz_yyyyy = pbuffer.data(idx_dip_gh + 813);

    auto tr_z_xyzz_yyyyz = pbuffer.data(idx_dip_gh + 814);

    auto tr_z_xyzz_yyyzz = pbuffer.data(idx_dip_gh + 815);

    auto tr_z_xyzz_yyzzz = pbuffer.data(idx_dip_gh + 816);

    auto tr_z_xyzz_yzzzz = pbuffer.data(idx_dip_gh + 817);

    auto tr_z_xyzz_zzzzz = pbuffer.data(idx_dip_gh + 818);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyzz_xxxxx, tr_z_xyzz_xxxxy, tr_z_xyzz_xxxxz, tr_z_xyzz_xxxyy, tr_z_xyzz_xxxyz, tr_z_xyzz_xxxzz, tr_z_xyzz_xxyyy, tr_z_xyzz_xxyyz, tr_z_xyzz_xxyzz, tr_z_xyzz_xxzzz, tr_z_xyzz_xyyyy, tr_z_xyzz_xyyyz, tr_z_xyzz_xyyzz, tr_z_xyzz_xyzzz, tr_z_xyzz_xzzzz, tr_z_xyzz_yyyyy, tr_z_xyzz_yyyyz, tr_z_xyzz_yyyzz, tr_z_xyzz_yyzzz, tr_z_xyzz_yzzzz, tr_z_xyzz_zzzzz, tr_z_xzz_xxxxx, tr_z_xzz_xxxxz, tr_z_xzz_xxxzz, tr_z_xzz_xxzzz, tr_z_xzz_xzzzz, tr_z_yzz_xxxxy, tr_z_yzz_xxxy, tr_z_yzz_xxxyy, tr_z_yzz_xxxyz, tr_z_yzz_xxyy, tr_z_yzz_xxyyy, tr_z_yzz_xxyyz, tr_z_yzz_xxyz, tr_z_yzz_xxyzz, tr_z_yzz_xyyy, tr_z_yzz_xyyyy, tr_z_yzz_xyyyz, tr_z_yzz_xyyz, tr_z_yzz_xyyzz, tr_z_yzz_xyzz, tr_z_yzz_xyzzz, tr_z_yzz_yyyy, tr_z_yzz_yyyyy, tr_z_yzz_yyyyz, tr_z_yzz_yyyz, tr_z_yzz_yyyzz, tr_z_yzz_yyzz, tr_z_yzz_yyzzz, tr_z_yzz_yzzz, tr_z_yzz_yzzzz, tr_z_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzz_xxxxx[i] = tr_z_xzz_xxxxx[i] * pa_y[i];

        tr_z_xyzz_xxxxy[i] = 4.0 * tr_z_yzz_xxxy[i] * fe_0 + tr_z_yzz_xxxxy[i] * pa_x[i];

        tr_z_xyzz_xxxxz[i] = tr_z_xzz_xxxxz[i] * pa_y[i];

        tr_z_xyzz_xxxyy[i] = 3.0 * tr_z_yzz_xxyy[i] * fe_0 + tr_z_yzz_xxxyy[i] * pa_x[i];

        tr_z_xyzz_xxxyz[i] = 3.0 * tr_z_yzz_xxyz[i] * fe_0 + tr_z_yzz_xxxyz[i] * pa_x[i];

        tr_z_xyzz_xxxzz[i] = tr_z_xzz_xxxzz[i] * pa_y[i];

        tr_z_xyzz_xxyyy[i] = 2.0 * tr_z_yzz_xyyy[i] * fe_0 + tr_z_yzz_xxyyy[i] * pa_x[i];

        tr_z_xyzz_xxyyz[i] = 2.0 * tr_z_yzz_xyyz[i] * fe_0 + tr_z_yzz_xxyyz[i] * pa_x[i];

        tr_z_xyzz_xxyzz[i] = 2.0 * tr_z_yzz_xyzz[i] * fe_0 + tr_z_yzz_xxyzz[i] * pa_x[i];

        tr_z_xyzz_xxzzz[i] = tr_z_xzz_xxzzz[i] * pa_y[i];

        tr_z_xyzz_xyyyy[i] = tr_z_yzz_yyyy[i] * fe_0 + tr_z_yzz_xyyyy[i] * pa_x[i];

        tr_z_xyzz_xyyyz[i] = tr_z_yzz_yyyz[i] * fe_0 + tr_z_yzz_xyyyz[i] * pa_x[i];

        tr_z_xyzz_xyyzz[i] = tr_z_yzz_yyzz[i] * fe_0 + tr_z_yzz_xyyzz[i] * pa_x[i];

        tr_z_xyzz_xyzzz[i] = tr_z_yzz_yzzz[i] * fe_0 + tr_z_yzz_xyzzz[i] * pa_x[i];

        tr_z_xyzz_xzzzz[i] = tr_z_xzz_xzzzz[i] * pa_y[i];

        tr_z_xyzz_yyyyy[i] = tr_z_yzz_yyyyy[i] * pa_x[i];

        tr_z_xyzz_yyyyz[i] = tr_z_yzz_yyyyz[i] * pa_x[i];

        tr_z_xyzz_yyyzz[i] = tr_z_yzz_yyyzz[i] * pa_x[i];

        tr_z_xyzz_yyzzz[i] = tr_z_yzz_yyzzz[i] * pa_x[i];

        tr_z_xyzz_yzzzz[i] = tr_z_yzz_yzzzz[i] * pa_x[i];

        tr_z_xyzz_zzzzz[i] = tr_z_yzz_zzzzz[i] * pa_x[i];
    }

    // Set up 819-840 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, tr_z_xzzz_xxxxx, tr_z_xzzz_xxxxy, tr_z_xzzz_xxxxz, tr_z_xzzz_xxxyy, tr_z_xzzz_xxxyz, tr_z_xzzz_xxxzz, tr_z_xzzz_xxyyy, tr_z_xzzz_xxyyz, tr_z_xzzz_xxyzz, tr_z_xzzz_xxzzz, tr_z_xzzz_xyyyy, tr_z_xzzz_xyyyz, tr_z_xzzz_xyyzz, tr_z_xzzz_xyzzz, tr_z_xzzz_xzzzz, tr_z_xzzz_yyyyy, tr_z_xzzz_yyyyz, tr_z_xzzz_yyyzz, tr_z_xzzz_yyzzz, tr_z_xzzz_yzzzz, tr_z_xzzz_zzzzz, tr_z_zzz_xxxx, tr_z_zzz_xxxxx, tr_z_zzz_xxxxy, tr_z_zzz_xxxxz, tr_z_zzz_xxxy, tr_z_zzz_xxxyy, tr_z_zzz_xxxyz, tr_z_zzz_xxxz, tr_z_zzz_xxxzz, tr_z_zzz_xxyy, tr_z_zzz_xxyyy, tr_z_zzz_xxyyz, tr_z_zzz_xxyz, tr_z_zzz_xxyzz, tr_z_zzz_xxzz, tr_z_zzz_xxzzz, tr_z_zzz_xyyy, tr_z_zzz_xyyyy, tr_z_zzz_xyyyz, tr_z_zzz_xyyz, tr_z_zzz_xyyzz, tr_z_zzz_xyzz, tr_z_zzz_xyzzz, tr_z_zzz_xzzz, tr_z_zzz_xzzzz, tr_z_zzz_yyyy, tr_z_zzz_yyyyy, tr_z_zzz_yyyyz, tr_z_zzz_yyyz, tr_z_zzz_yyyzz, tr_z_zzz_yyzz, tr_z_zzz_yyzzz, tr_z_zzz_yzzz, tr_z_zzz_yzzzz, tr_z_zzz_zzzz, tr_z_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzz_xxxxx[i] = 5.0 * tr_z_zzz_xxxx[i] * fe_0 + tr_z_zzz_xxxxx[i] * pa_x[i];

        tr_z_xzzz_xxxxy[i] = 4.0 * tr_z_zzz_xxxy[i] * fe_0 + tr_z_zzz_xxxxy[i] * pa_x[i];

        tr_z_xzzz_xxxxz[i] = 4.0 * tr_z_zzz_xxxz[i] * fe_0 + tr_z_zzz_xxxxz[i] * pa_x[i];

        tr_z_xzzz_xxxyy[i] = 3.0 * tr_z_zzz_xxyy[i] * fe_0 + tr_z_zzz_xxxyy[i] * pa_x[i];

        tr_z_xzzz_xxxyz[i] = 3.0 * tr_z_zzz_xxyz[i] * fe_0 + tr_z_zzz_xxxyz[i] * pa_x[i];

        tr_z_xzzz_xxxzz[i] = 3.0 * tr_z_zzz_xxzz[i] * fe_0 + tr_z_zzz_xxxzz[i] * pa_x[i];

        tr_z_xzzz_xxyyy[i] = 2.0 * tr_z_zzz_xyyy[i] * fe_0 + tr_z_zzz_xxyyy[i] * pa_x[i];

        tr_z_xzzz_xxyyz[i] = 2.0 * tr_z_zzz_xyyz[i] * fe_0 + tr_z_zzz_xxyyz[i] * pa_x[i];

        tr_z_xzzz_xxyzz[i] = 2.0 * tr_z_zzz_xyzz[i] * fe_0 + tr_z_zzz_xxyzz[i] * pa_x[i];

        tr_z_xzzz_xxzzz[i] = 2.0 * tr_z_zzz_xzzz[i] * fe_0 + tr_z_zzz_xxzzz[i] * pa_x[i];

        tr_z_xzzz_xyyyy[i] = tr_z_zzz_yyyy[i] * fe_0 + tr_z_zzz_xyyyy[i] * pa_x[i];

        tr_z_xzzz_xyyyz[i] = tr_z_zzz_yyyz[i] * fe_0 + tr_z_zzz_xyyyz[i] * pa_x[i];

        tr_z_xzzz_xyyzz[i] = tr_z_zzz_yyzz[i] * fe_0 + tr_z_zzz_xyyzz[i] * pa_x[i];

        tr_z_xzzz_xyzzz[i] = tr_z_zzz_yzzz[i] * fe_0 + tr_z_zzz_xyzzz[i] * pa_x[i];

        tr_z_xzzz_xzzzz[i] = tr_z_zzz_zzzz[i] * fe_0 + tr_z_zzz_xzzzz[i] * pa_x[i];

        tr_z_xzzz_yyyyy[i] = tr_z_zzz_yyyyy[i] * pa_x[i];

        tr_z_xzzz_yyyyz[i] = tr_z_zzz_yyyyz[i] * pa_x[i];

        tr_z_xzzz_yyyzz[i] = tr_z_zzz_yyyzz[i] * pa_x[i];

        tr_z_xzzz_yyzzz[i] = tr_z_zzz_yyzzz[i] * pa_x[i];

        tr_z_xzzz_yzzzz[i] = tr_z_zzz_yzzzz[i] * pa_x[i];

        tr_z_xzzz_zzzzz[i] = tr_z_zzz_zzzzz[i] * pa_x[i];
    }

    // Set up 840-861 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, tr_z_yy_xxxxx, tr_z_yy_xxxxy, tr_z_yy_xxxxz, tr_z_yy_xxxyy, tr_z_yy_xxxyz, tr_z_yy_xxxzz, tr_z_yy_xxyyy, tr_z_yy_xxyyz, tr_z_yy_xxyzz, tr_z_yy_xxzzz, tr_z_yy_xyyyy, tr_z_yy_xyyyz, tr_z_yy_xyyzz, tr_z_yy_xyzzz, tr_z_yy_xzzzz, tr_z_yy_yyyyy, tr_z_yy_yyyyz, tr_z_yy_yyyzz, tr_z_yy_yyzzz, tr_z_yy_yzzzz, tr_z_yy_zzzzz, tr_z_yyy_xxxx, tr_z_yyy_xxxxx, tr_z_yyy_xxxxy, tr_z_yyy_xxxxz, tr_z_yyy_xxxy, tr_z_yyy_xxxyy, tr_z_yyy_xxxyz, tr_z_yyy_xxxz, tr_z_yyy_xxxzz, tr_z_yyy_xxyy, tr_z_yyy_xxyyy, tr_z_yyy_xxyyz, tr_z_yyy_xxyz, tr_z_yyy_xxyzz, tr_z_yyy_xxzz, tr_z_yyy_xxzzz, tr_z_yyy_xyyy, tr_z_yyy_xyyyy, tr_z_yyy_xyyyz, tr_z_yyy_xyyz, tr_z_yyy_xyyzz, tr_z_yyy_xyzz, tr_z_yyy_xyzzz, tr_z_yyy_xzzz, tr_z_yyy_xzzzz, tr_z_yyy_yyyy, tr_z_yyy_yyyyy, tr_z_yyy_yyyyz, tr_z_yyy_yyyz, tr_z_yyy_yyyzz, tr_z_yyy_yyzz, tr_z_yyy_yyzzz, tr_z_yyy_yzzz, tr_z_yyy_yzzzz, tr_z_yyy_zzzz, tr_z_yyy_zzzzz, tr_z_yyyy_xxxxx, tr_z_yyyy_xxxxy, tr_z_yyyy_xxxxz, tr_z_yyyy_xxxyy, tr_z_yyyy_xxxyz, tr_z_yyyy_xxxzz, tr_z_yyyy_xxyyy, tr_z_yyyy_xxyyz, tr_z_yyyy_xxyzz, tr_z_yyyy_xxzzz, tr_z_yyyy_xyyyy, tr_z_yyyy_xyyyz, tr_z_yyyy_xyyzz, tr_z_yyyy_xyzzz, tr_z_yyyy_xzzzz, tr_z_yyyy_yyyyy, tr_z_yyyy_yyyyz, tr_z_yyyy_yyyzz, tr_z_yyyy_yyzzz, tr_z_yyyy_yzzzz, tr_z_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyy_xxxxx[i] = 3.0 * tr_z_yy_xxxxx[i] * fe_0 + tr_z_yyy_xxxxx[i] * pa_y[i];

        tr_z_yyyy_xxxxy[i] = 3.0 * tr_z_yy_xxxxy[i] * fe_0 + tr_z_yyy_xxxx[i] * fe_0 + tr_z_yyy_xxxxy[i] * pa_y[i];

        tr_z_yyyy_xxxxz[i] = 3.0 * tr_z_yy_xxxxz[i] * fe_0 + tr_z_yyy_xxxxz[i] * pa_y[i];

        tr_z_yyyy_xxxyy[i] = 3.0 * tr_z_yy_xxxyy[i] * fe_0 + 2.0 * tr_z_yyy_xxxy[i] * fe_0 + tr_z_yyy_xxxyy[i] * pa_y[i];

        tr_z_yyyy_xxxyz[i] = 3.0 * tr_z_yy_xxxyz[i] * fe_0 + tr_z_yyy_xxxz[i] * fe_0 + tr_z_yyy_xxxyz[i] * pa_y[i];

        tr_z_yyyy_xxxzz[i] = 3.0 * tr_z_yy_xxxzz[i] * fe_0 + tr_z_yyy_xxxzz[i] * pa_y[i];

        tr_z_yyyy_xxyyy[i] = 3.0 * tr_z_yy_xxyyy[i] * fe_0 + 3.0 * tr_z_yyy_xxyy[i] * fe_0 + tr_z_yyy_xxyyy[i] * pa_y[i];

        tr_z_yyyy_xxyyz[i] = 3.0 * tr_z_yy_xxyyz[i] * fe_0 + 2.0 * tr_z_yyy_xxyz[i] * fe_0 + tr_z_yyy_xxyyz[i] * pa_y[i];

        tr_z_yyyy_xxyzz[i] = 3.0 * tr_z_yy_xxyzz[i] * fe_0 + tr_z_yyy_xxzz[i] * fe_0 + tr_z_yyy_xxyzz[i] * pa_y[i];

        tr_z_yyyy_xxzzz[i] = 3.0 * tr_z_yy_xxzzz[i] * fe_0 + tr_z_yyy_xxzzz[i] * pa_y[i];

        tr_z_yyyy_xyyyy[i] = 3.0 * tr_z_yy_xyyyy[i] * fe_0 + 4.0 * tr_z_yyy_xyyy[i] * fe_0 + tr_z_yyy_xyyyy[i] * pa_y[i];

        tr_z_yyyy_xyyyz[i] = 3.0 * tr_z_yy_xyyyz[i] * fe_0 + 3.0 * tr_z_yyy_xyyz[i] * fe_0 + tr_z_yyy_xyyyz[i] * pa_y[i];

        tr_z_yyyy_xyyzz[i] = 3.0 * tr_z_yy_xyyzz[i] * fe_0 + 2.0 * tr_z_yyy_xyzz[i] * fe_0 + tr_z_yyy_xyyzz[i] * pa_y[i];

        tr_z_yyyy_xyzzz[i] = 3.0 * tr_z_yy_xyzzz[i] * fe_0 + tr_z_yyy_xzzz[i] * fe_0 + tr_z_yyy_xyzzz[i] * pa_y[i];

        tr_z_yyyy_xzzzz[i] = 3.0 * tr_z_yy_xzzzz[i] * fe_0 + tr_z_yyy_xzzzz[i] * pa_y[i];

        tr_z_yyyy_yyyyy[i] = 3.0 * tr_z_yy_yyyyy[i] * fe_0 + 5.0 * tr_z_yyy_yyyy[i] * fe_0 + tr_z_yyy_yyyyy[i] * pa_y[i];

        tr_z_yyyy_yyyyz[i] = 3.0 * tr_z_yy_yyyyz[i] * fe_0 + 4.0 * tr_z_yyy_yyyz[i] * fe_0 + tr_z_yyy_yyyyz[i] * pa_y[i];

        tr_z_yyyy_yyyzz[i] = 3.0 * tr_z_yy_yyyzz[i] * fe_0 + 3.0 * tr_z_yyy_yyzz[i] * fe_0 + tr_z_yyy_yyyzz[i] * pa_y[i];

        tr_z_yyyy_yyzzz[i] = 3.0 * tr_z_yy_yyzzz[i] * fe_0 + 2.0 * tr_z_yyy_yzzz[i] * fe_0 + tr_z_yyy_yyzzz[i] * pa_y[i];

        tr_z_yyyy_yzzzz[i] = 3.0 * tr_z_yy_yzzzz[i] * fe_0 + tr_z_yyy_zzzz[i] * fe_0 + tr_z_yyy_yzzzz[i] * pa_y[i];

        tr_z_yyyy_zzzzz[i] = 3.0 * tr_z_yy_zzzzz[i] * fe_0 + tr_z_yyy_zzzzz[i] * pa_y[i];
    }

    // Set up 861-882 components of targeted buffer : GH

    auto tr_z_yyyz_xxxxx = pbuffer.data(idx_dip_gh + 861);

    auto tr_z_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 862);

    auto tr_z_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 863);

    auto tr_z_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 864);

    auto tr_z_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 865);

    auto tr_z_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 866);

    auto tr_z_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 867);

    auto tr_z_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 868);

    auto tr_z_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 869);

    auto tr_z_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 870);

    auto tr_z_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 871);

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

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yyy_xxxxy, tr_z_yyy_xxxyy, tr_z_yyy_xxyyy, tr_z_yyy_xyyyy, tr_z_yyy_yyyyy, tr_z_yyyz_xxxxx, tr_z_yyyz_xxxxy, tr_z_yyyz_xxxxz, tr_z_yyyz_xxxyy, tr_z_yyyz_xxxyz, tr_z_yyyz_xxxzz, tr_z_yyyz_xxyyy, tr_z_yyyz_xxyyz, tr_z_yyyz_xxyzz, tr_z_yyyz_xxzzz, tr_z_yyyz_xyyyy, tr_z_yyyz_xyyyz, tr_z_yyyz_xyyzz, tr_z_yyyz_xyzzz, tr_z_yyyz_xzzzz, tr_z_yyyz_yyyyy, tr_z_yyyz_yyyyz, tr_z_yyyz_yyyzz, tr_z_yyyz_yyzzz, tr_z_yyyz_yzzzz, tr_z_yyyz_zzzzz, tr_z_yyz_xxxxx, tr_z_yyz_xxxxz, tr_z_yyz_xxxyz, tr_z_yyz_xxxz, tr_z_yyz_xxxzz, tr_z_yyz_xxyyz, tr_z_yyz_xxyz, tr_z_yyz_xxyzz, tr_z_yyz_xxzz, tr_z_yyz_xxzzz, tr_z_yyz_xyyyz, tr_z_yyz_xyyz, tr_z_yyz_xyyzz, tr_z_yyz_xyzz, tr_z_yyz_xyzzz, tr_z_yyz_xzzz, tr_z_yyz_xzzzz, tr_z_yyz_yyyyz, tr_z_yyz_yyyz, tr_z_yyz_yyyzz, tr_z_yyz_yyzz, tr_z_yyz_yyzzz, tr_z_yyz_yzzz, tr_z_yyz_yzzzz, tr_z_yyz_zzzz, tr_z_yyz_zzzzz, tr_z_yz_xxxxx, tr_z_yz_xxxxz, tr_z_yz_xxxyz, tr_z_yz_xxxzz, tr_z_yz_xxyyz, tr_z_yz_xxyzz, tr_z_yz_xxzzz, tr_z_yz_xyyyz, tr_z_yz_xyyzz, tr_z_yz_xyzzz, tr_z_yz_xzzzz, tr_z_yz_yyyyz, tr_z_yz_yyyzz, tr_z_yz_yyzzz, tr_z_yz_yzzzz, tr_z_yz_zzzzz, ts_yyy_xxxxy, ts_yyy_xxxyy, ts_yyy_xxyyy, ts_yyy_xyyyy, ts_yyy_yyyyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyz_xxxxx[i] = 2.0 * tr_z_yz_xxxxx[i] * fe_0 + tr_z_yyz_xxxxx[i] * pa_y[i];

        tr_z_yyyz_xxxxy[i] = ts_yyy_xxxxy[i] * fe_0 + tr_z_yyy_xxxxy[i] * pa_z[i];

        tr_z_yyyz_xxxxz[i] = 2.0 * tr_z_yz_xxxxz[i] * fe_0 + tr_z_yyz_xxxxz[i] * pa_y[i];

        tr_z_yyyz_xxxyy[i] = ts_yyy_xxxyy[i] * fe_0 + tr_z_yyy_xxxyy[i] * pa_z[i];

        tr_z_yyyz_xxxyz[i] = 2.0 * tr_z_yz_xxxyz[i] * fe_0 + tr_z_yyz_xxxz[i] * fe_0 + tr_z_yyz_xxxyz[i] * pa_y[i];

        tr_z_yyyz_xxxzz[i] = 2.0 * tr_z_yz_xxxzz[i] * fe_0 + tr_z_yyz_xxxzz[i] * pa_y[i];

        tr_z_yyyz_xxyyy[i] = ts_yyy_xxyyy[i] * fe_0 + tr_z_yyy_xxyyy[i] * pa_z[i];

        tr_z_yyyz_xxyyz[i] = 2.0 * tr_z_yz_xxyyz[i] * fe_0 + 2.0 * tr_z_yyz_xxyz[i] * fe_0 + tr_z_yyz_xxyyz[i] * pa_y[i];

        tr_z_yyyz_xxyzz[i] = 2.0 * tr_z_yz_xxyzz[i] * fe_0 + tr_z_yyz_xxzz[i] * fe_0 + tr_z_yyz_xxyzz[i] * pa_y[i];

        tr_z_yyyz_xxzzz[i] = 2.0 * tr_z_yz_xxzzz[i] * fe_0 + tr_z_yyz_xxzzz[i] * pa_y[i];

        tr_z_yyyz_xyyyy[i] = ts_yyy_xyyyy[i] * fe_0 + tr_z_yyy_xyyyy[i] * pa_z[i];

        tr_z_yyyz_xyyyz[i] = 2.0 * tr_z_yz_xyyyz[i] * fe_0 + 3.0 * tr_z_yyz_xyyz[i] * fe_0 + tr_z_yyz_xyyyz[i] * pa_y[i];

        tr_z_yyyz_xyyzz[i] = 2.0 * tr_z_yz_xyyzz[i] * fe_0 + 2.0 * tr_z_yyz_xyzz[i] * fe_0 + tr_z_yyz_xyyzz[i] * pa_y[i];

        tr_z_yyyz_xyzzz[i] = 2.0 * tr_z_yz_xyzzz[i] * fe_0 + tr_z_yyz_xzzz[i] * fe_0 + tr_z_yyz_xyzzz[i] * pa_y[i];

        tr_z_yyyz_xzzzz[i] = 2.0 * tr_z_yz_xzzzz[i] * fe_0 + tr_z_yyz_xzzzz[i] * pa_y[i];

        tr_z_yyyz_yyyyy[i] = ts_yyy_yyyyy[i] * fe_0 + tr_z_yyy_yyyyy[i] * pa_z[i];

        tr_z_yyyz_yyyyz[i] = 2.0 * tr_z_yz_yyyyz[i] * fe_0 + 4.0 * tr_z_yyz_yyyz[i] * fe_0 + tr_z_yyz_yyyyz[i] * pa_y[i];

        tr_z_yyyz_yyyzz[i] = 2.0 * tr_z_yz_yyyzz[i] * fe_0 + 3.0 * tr_z_yyz_yyzz[i] * fe_0 + tr_z_yyz_yyyzz[i] * pa_y[i];

        tr_z_yyyz_yyzzz[i] = 2.0 * tr_z_yz_yyzzz[i] * fe_0 + 2.0 * tr_z_yyz_yzzz[i] * fe_0 + tr_z_yyz_yyzzz[i] * pa_y[i];

        tr_z_yyyz_yzzzz[i] = 2.0 * tr_z_yz_yzzzz[i] * fe_0 + tr_z_yyz_zzzz[i] * fe_0 + tr_z_yyz_yzzzz[i] * pa_y[i];

        tr_z_yyyz_zzzzz[i] = 2.0 * tr_z_yz_zzzzz[i] * fe_0 + tr_z_yyz_zzzzz[i] * pa_y[i];
    }

    // Set up 882-903 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, tr_z_yyzz_xxxxx, tr_z_yyzz_xxxxy, tr_z_yyzz_xxxxz, tr_z_yyzz_xxxyy, tr_z_yyzz_xxxyz, tr_z_yyzz_xxxzz, tr_z_yyzz_xxyyy, tr_z_yyzz_xxyyz, tr_z_yyzz_xxyzz, tr_z_yyzz_xxzzz, tr_z_yyzz_xyyyy, tr_z_yyzz_xyyyz, tr_z_yyzz_xyyzz, tr_z_yyzz_xyzzz, tr_z_yyzz_xzzzz, tr_z_yyzz_yyyyy, tr_z_yyzz_yyyyz, tr_z_yyzz_yyyzz, tr_z_yyzz_yyzzz, tr_z_yyzz_yzzzz, tr_z_yyzz_zzzzz, tr_z_yzz_xxxx, tr_z_yzz_xxxxx, tr_z_yzz_xxxxy, tr_z_yzz_xxxxz, tr_z_yzz_xxxy, tr_z_yzz_xxxyy, tr_z_yzz_xxxyz, tr_z_yzz_xxxz, tr_z_yzz_xxxzz, tr_z_yzz_xxyy, tr_z_yzz_xxyyy, tr_z_yzz_xxyyz, tr_z_yzz_xxyz, tr_z_yzz_xxyzz, tr_z_yzz_xxzz, tr_z_yzz_xxzzz, tr_z_yzz_xyyy, tr_z_yzz_xyyyy, tr_z_yzz_xyyyz, tr_z_yzz_xyyz, tr_z_yzz_xyyzz, tr_z_yzz_xyzz, tr_z_yzz_xyzzz, tr_z_yzz_xzzz, tr_z_yzz_xzzzz, tr_z_yzz_yyyy, tr_z_yzz_yyyyy, tr_z_yzz_yyyyz, tr_z_yzz_yyyz, tr_z_yzz_yyyzz, tr_z_yzz_yyzz, tr_z_yzz_yyzzz, tr_z_yzz_yzzz, tr_z_yzz_yzzzz, tr_z_yzz_zzzz, tr_z_yzz_zzzzz, tr_z_zz_xxxxx, tr_z_zz_xxxxy, tr_z_zz_xxxxz, tr_z_zz_xxxyy, tr_z_zz_xxxyz, tr_z_zz_xxxzz, tr_z_zz_xxyyy, tr_z_zz_xxyyz, tr_z_zz_xxyzz, tr_z_zz_xxzzz, tr_z_zz_xyyyy, tr_z_zz_xyyyz, tr_z_zz_xyyzz, tr_z_zz_xyzzz, tr_z_zz_xzzzz, tr_z_zz_yyyyy, tr_z_zz_yyyyz, tr_z_zz_yyyzz, tr_z_zz_yyzzz, tr_z_zz_yzzzz, tr_z_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzz_xxxxx[i] = tr_z_zz_xxxxx[i] * fe_0 + tr_z_yzz_xxxxx[i] * pa_y[i];

        tr_z_yyzz_xxxxy[i] = tr_z_zz_xxxxy[i] * fe_0 + tr_z_yzz_xxxx[i] * fe_0 + tr_z_yzz_xxxxy[i] * pa_y[i];

        tr_z_yyzz_xxxxz[i] = tr_z_zz_xxxxz[i] * fe_0 + tr_z_yzz_xxxxz[i] * pa_y[i];

        tr_z_yyzz_xxxyy[i] = tr_z_zz_xxxyy[i] * fe_0 + 2.0 * tr_z_yzz_xxxy[i] * fe_0 + tr_z_yzz_xxxyy[i] * pa_y[i];

        tr_z_yyzz_xxxyz[i] = tr_z_zz_xxxyz[i] * fe_0 + tr_z_yzz_xxxz[i] * fe_0 + tr_z_yzz_xxxyz[i] * pa_y[i];

        tr_z_yyzz_xxxzz[i] = tr_z_zz_xxxzz[i] * fe_0 + tr_z_yzz_xxxzz[i] * pa_y[i];

        tr_z_yyzz_xxyyy[i] = tr_z_zz_xxyyy[i] * fe_0 + 3.0 * tr_z_yzz_xxyy[i] * fe_0 + tr_z_yzz_xxyyy[i] * pa_y[i];

        tr_z_yyzz_xxyyz[i] = tr_z_zz_xxyyz[i] * fe_0 + 2.0 * tr_z_yzz_xxyz[i] * fe_0 + tr_z_yzz_xxyyz[i] * pa_y[i];

        tr_z_yyzz_xxyzz[i] = tr_z_zz_xxyzz[i] * fe_0 + tr_z_yzz_xxzz[i] * fe_0 + tr_z_yzz_xxyzz[i] * pa_y[i];

        tr_z_yyzz_xxzzz[i] = tr_z_zz_xxzzz[i] * fe_0 + tr_z_yzz_xxzzz[i] * pa_y[i];

        tr_z_yyzz_xyyyy[i] = tr_z_zz_xyyyy[i] * fe_0 + 4.0 * tr_z_yzz_xyyy[i] * fe_0 + tr_z_yzz_xyyyy[i] * pa_y[i];

        tr_z_yyzz_xyyyz[i] = tr_z_zz_xyyyz[i] * fe_0 + 3.0 * tr_z_yzz_xyyz[i] * fe_0 + tr_z_yzz_xyyyz[i] * pa_y[i];

        tr_z_yyzz_xyyzz[i] = tr_z_zz_xyyzz[i] * fe_0 + 2.0 * tr_z_yzz_xyzz[i] * fe_0 + tr_z_yzz_xyyzz[i] * pa_y[i];

        tr_z_yyzz_xyzzz[i] = tr_z_zz_xyzzz[i] * fe_0 + tr_z_yzz_xzzz[i] * fe_0 + tr_z_yzz_xyzzz[i] * pa_y[i];

        tr_z_yyzz_xzzzz[i] = tr_z_zz_xzzzz[i] * fe_0 + tr_z_yzz_xzzzz[i] * pa_y[i];

        tr_z_yyzz_yyyyy[i] = tr_z_zz_yyyyy[i] * fe_0 + 5.0 * tr_z_yzz_yyyy[i] * fe_0 + tr_z_yzz_yyyyy[i] * pa_y[i];

        tr_z_yyzz_yyyyz[i] = tr_z_zz_yyyyz[i] * fe_0 + 4.0 * tr_z_yzz_yyyz[i] * fe_0 + tr_z_yzz_yyyyz[i] * pa_y[i];

        tr_z_yyzz_yyyzz[i] = tr_z_zz_yyyzz[i] * fe_0 + 3.0 * tr_z_yzz_yyzz[i] * fe_0 + tr_z_yzz_yyyzz[i] * pa_y[i];

        tr_z_yyzz_yyzzz[i] = tr_z_zz_yyzzz[i] * fe_0 + 2.0 * tr_z_yzz_yzzz[i] * fe_0 + tr_z_yzz_yyzzz[i] * pa_y[i];

        tr_z_yyzz_yzzzz[i] = tr_z_zz_yzzzz[i] * fe_0 + tr_z_yzz_zzzz[i] * fe_0 + tr_z_yzz_yzzzz[i] * pa_y[i];

        tr_z_yyzz_zzzzz[i] = tr_z_zz_zzzzz[i] * fe_0 + tr_z_yzz_zzzzz[i] * pa_y[i];
    }

    // Set up 903-924 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, tr_z_yzzz_xxxxx, tr_z_yzzz_xxxxy, tr_z_yzzz_xxxxz, tr_z_yzzz_xxxyy, tr_z_yzzz_xxxyz, tr_z_yzzz_xxxzz, tr_z_yzzz_xxyyy, tr_z_yzzz_xxyyz, tr_z_yzzz_xxyzz, tr_z_yzzz_xxzzz, tr_z_yzzz_xyyyy, tr_z_yzzz_xyyyz, tr_z_yzzz_xyyzz, tr_z_yzzz_xyzzz, tr_z_yzzz_xzzzz, tr_z_yzzz_yyyyy, tr_z_yzzz_yyyyz, tr_z_yzzz_yyyzz, tr_z_yzzz_yyzzz, tr_z_yzzz_yzzzz, tr_z_yzzz_zzzzz, tr_z_zzz_xxxx, tr_z_zzz_xxxxx, tr_z_zzz_xxxxy, tr_z_zzz_xxxxz, tr_z_zzz_xxxy, tr_z_zzz_xxxyy, tr_z_zzz_xxxyz, tr_z_zzz_xxxz, tr_z_zzz_xxxzz, tr_z_zzz_xxyy, tr_z_zzz_xxyyy, tr_z_zzz_xxyyz, tr_z_zzz_xxyz, tr_z_zzz_xxyzz, tr_z_zzz_xxzz, tr_z_zzz_xxzzz, tr_z_zzz_xyyy, tr_z_zzz_xyyyy, tr_z_zzz_xyyyz, tr_z_zzz_xyyz, tr_z_zzz_xyyzz, tr_z_zzz_xyzz, tr_z_zzz_xyzzz, tr_z_zzz_xzzz, tr_z_zzz_xzzzz, tr_z_zzz_yyyy, tr_z_zzz_yyyyy, tr_z_zzz_yyyyz, tr_z_zzz_yyyz, tr_z_zzz_yyyzz, tr_z_zzz_yyzz, tr_z_zzz_yyzzz, tr_z_zzz_yzzz, tr_z_zzz_yzzzz, tr_z_zzz_zzzz, tr_z_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzz_xxxxx[i] = tr_z_zzz_xxxxx[i] * pa_y[i];

        tr_z_yzzz_xxxxy[i] = tr_z_zzz_xxxx[i] * fe_0 + tr_z_zzz_xxxxy[i] * pa_y[i];

        tr_z_yzzz_xxxxz[i] = tr_z_zzz_xxxxz[i] * pa_y[i];

        tr_z_yzzz_xxxyy[i] = 2.0 * tr_z_zzz_xxxy[i] * fe_0 + tr_z_zzz_xxxyy[i] * pa_y[i];

        tr_z_yzzz_xxxyz[i] = tr_z_zzz_xxxz[i] * fe_0 + tr_z_zzz_xxxyz[i] * pa_y[i];

        tr_z_yzzz_xxxzz[i] = tr_z_zzz_xxxzz[i] * pa_y[i];

        tr_z_yzzz_xxyyy[i] = 3.0 * tr_z_zzz_xxyy[i] * fe_0 + tr_z_zzz_xxyyy[i] * pa_y[i];

        tr_z_yzzz_xxyyz[i] = 2.0 * tr_z_zzz_xxyz[i] * fe_0 + tr_z_zzz_xxyyz[i] * pa_y[i];

        tr_z_yzzz_xxyzz[i] = tr_z_zzz_xxzz[i] * fe_0 + tr_z_zzz_xxyzz[i] * pa_y[i];

        tr_z_yzzz_xxzzz[i] = tr_z_zzz_xxzzz[i] * pa_y[i];

        tr_z_yzzz_xyyyy[i] = 4.0 * tr_z_zzz_xyyy[i] * fe_0 + tr_z_zzz_xyyyy[i] * pa_y[i];

        tr_z_yzzz_xyyyz[i] = 3.0 * tr_z_zzz_xyyz[i] * fe_0 + tr_z_zzz_xyyyz[i] * pa_y[i];

        tr_z_yzzz_xyyzz[i] = 2.0 * tr_z_zzz_xyzz[i] * fe_0 + tr_z_zzz_xyyzz[i] * pa_y[i];

        tr_z_yzzz_xyzzz[i] = tr_z_zzz_xzzz[i] * fe_0 + tr_z_zzz_xyzzz[i] * pa_y[i];

        tr_z_yzzz_xzzzz[i] = tr_z_zzz_xzzzz[i] * pa_y[i];

        tr_z_yzzz_yyyyy[i] = 5.0 * tr_z_zzz_yyyy[i] * fe_0 + tr_z_zzz_yyyyy[i] * pa_y[i];

        tr_z_yzzz_yyyyz[i] = 4.0 * tr_z_zzz_yyyz[i] * fe_0 + tr_z_zzz_yyyyz[i] * pa_y[i];

        tr_z_yzzz_yyyzz[i] = 3.0 * tr_z_zzz_yyzz[i] * fe_0 + tr_z_zzz_yyyzz[i] * pa_y[i];

        tr_z_yzzz_yyzzz[i] = 2.0 * tr_z_zzz_yzzz[i] * fe_0 + tr_z_zzz_yyzzz[i] * pa_y[i];

        tr_z_yzzz_yzzzz[i] = tr_z_zzz_zzzz[i] * fe_0 + tr_z_zzz_yzzzz[i] * pa_y[i];

        tr_z_yzzz_zzzzz[i] = tr_z_zzz_zzzzz[i] * pa_y[i];
    }

    // Set up 924-945 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_z, tr_z_zz_xxxxx, tr_z_zz_xxxxy, tr_z_zz_xxxxz, tr_z_zz_xxxyy, tr_z_zz_xxxyz, tr_z_zz_xxxzz, tr_z_zz_xxyyy, tr_z_zz_xxyyz, tr_z_zz_xxyzz, tr_z_zz_xxzzz, tr_z_zz_xyyyy, tr_z_zz_xyyyz, tr_z_zz_xyyzz, tr_z_zz_xyzzz, tr_z_zz_xzzzz, tr_z_zz_yyyyy, tr_z_zz_yyyyz, tr_z_zz_yyyzz, tr_z_zz_yyzzz, tr_z_zz_yzzzz, tr_z_zz_zzzzz, tr_z_zzz_xxxx, tr_z_zzz_xxxxx, tr_z_zzz_xxxxy, tr_z_zzz_xxxxz, tr_z_zzz_xxxy, tr_z_zzz_xxxyy, tr_z_zzz_xxxyz, tr_z_zzz_xxxz, tr_z_zzz_xxxzz, tr_z_zzz_xxyy, tr_z_zzz_xxyyy, tr_z_zzz_xxyyz, tr_z_zzz_xxyz, tr_z_zzz_xxyzz, tr_z_zzz_xxzz, tr_z_zzz_xxzzz, tr_z_zzz_xyyy, tr_z_zzz_xyyyy, tr_z_zzz_xyyyz, tr_z_zzz_xyyz, tr_z_zzz_xyyzz, tr_z_zzz_xyzz, tr_z_zzz_xyzzz, tr_z_zzz_xzzz, tr_z_zzz_xzzzz, tr_z_zzz_yyyy, tr_z_zzz_yyyyy, tr_z_zzz_yyyyz, tr_z_zzz_yyyz, tr_z_zzz_yyyzz, tr_z_zzz_yyzz, tr_z_zzz_yyzzz, tr_z_zzz_yzzz, tr_z_zzz_yzzzz, tr_z_zzz_zzzz, tr_z_zzz_zzzzz, tr_z_zzzz_xxxxx, tr_z_zzzz_xxxxy, tr_z_zzzz_xxxxz, tr_z_zzzz_xxxyy, tr_z_zzzz_xxxyz, tr_z_zzzz_xxxzz, tr_z_zzzz_xxyyy, tr_z_zzzz_xxyyz, tr_z_zzzz_xxyzz, tr_z_zzzz_xxzzz, tr_z_zzzz_xyyyy, tr_z_zzzz_xyyyz, tr_z_zzzz_xyyzz, tr_z_zzzz_xyzzz, tr_z_zzzz_xzzzz, tr_z_zzzz_yyyyy, tr_z_zzzz_yyyyz, tr_z_zzzz_yyyzz, tr_z_zzzz_yyzzz, tr_z_zzzz_yzzzz, tr_z_zzzz_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzz_xxxxx[i] = 3.0 * tr_z_zz_xxxxx[i] * fe_0 + ts_zzz_xxxxx[i] * fe_0 + tr_z_zzz_xxxxx[i] * pa_z[i];

        tr_z_zzzz_xxxxy[i] = 3.0 * tr_z_zz_xxxxy[i] * fe_0 + ts_zzz_xxxxy[i] * fe_0 + tr_z_zzz_xxxxy[i] * pa_z[i];

        tr_z_zzzz_xxxxz[i] = 3.0 * tr_z_zz_xxxxz[i] * fe_0 + tr_z_zzz_xxxx[i] * fe_0 + ts_zzz_xxxxz[i] * fe_0 + tr_z_zzz_xxxxz[i] * pa_z[i];

        tr_z_zzzz_xxxyy[i] = 3.0 * tr_z_zz_xxxyy[i] * fe_0 + ts_zzz_xxxyy[i] * fe_0 + tr_z_zzz_xxxyy[i] * pa_z[i];

        tr_z_zzzz_xxxyz[i] = 3.0 * tr_z_zz_xxxyz[i] * fe_0 + tr_z_zzz_xxxy[i] * fe_0 + ts_zzz_xxxyz[i] * fe_0 + tr_z_zzz_xxxyz[i] * pa_z[i];

        tr_z_zzzz_xxxzz[i] = 3.0 * tr_z_zz_xxxzz[i] * fe_0 + 2.0 * tr_z_zzz_xxxz[i] * fe_0 + ts_zzz_xxxzz[i] * fe_0 + tr_z_zzz_xxxzz[i] * pa_z[i];

        tr_z_zzzz_xxyyy[i] = 3.0 * tr_z_zz_xxyyy[i] * fe_0 + ts_zzz_xxyyy[i] * fe_0 + tr_z_zzz_xxyyy[i] * pa_z[i];

        tr_z_zzzz_xxyyz[i] = 3.0 * tr_z_zz_xxyyz[i] * fe_0 + tr_z_zzz_xxyy[i] * fe_0 + ts_zzz_xxyyz[i] * fe_0 + tr_z_zzz_xxyyz[i] * pa_z[i];

        tr_z_zzzz_xxyzz[i] = 3.0 * tr_z_zz_xxyzz[i] * fe_0 + 2.0 * tr_z_zzz_xxyz[i] * fe_0 + ts_zzz_xxyzz[i] * fe_0 + tr_z_zzz_xxyzz[i] * pa_z[i];

        tr_z_zzzz_xxzzz[i] = 3.0 * tr_z_zz_xxzzz[i] * fe_0 + 3.0 * tr_z_zzz_xxzz[i] * fe_0 + ts_zzz_xxzzz[i] * fe_0 + tr_z_zzz_xxzzz[i] * pa_z[i];

        tr_z_zzzz_xyyyy[i] = 3.0 * tr_z_zz_xyyyy[i] * fe_0 + ts_zzz_xyyyy[i] * fe_0 + tr_z_zzz_xyyyy[i] * pa_z[i];

        tr_z_zzzz_xyyyz[i] = 3.0 * tr_z_zz_xyyyz[i] * fe_0 + tr_z_zzz_xyyy[i] * fe_0 + ts_zzz_xyyyz[i] * fe_0 + tr_z_zzz_xyyyz[i] * pa_z[i];

        tr_z_zzzz_xyyzz[i] = 3.0 * tr_z_zz_xyyzz[i] * fe_0 + 2.0 * tr_z_zzz_xyyz[i] * fe_0 + ts_zzz_xyyzz[i] * fe_0 + tr_z_zzz_xyyzz[i] * pa_z[i];

        tr_z_zzzz_xyzzz[i] = 3.0 * tr_z_zz_xyzzz[i] * fe_0 + 3.0 * tr_z_zzz_xyzz[i] * fe_0 + ts_zzz_xyzzz[i] * fe_0 + tr_z_zzz_xyzzz[i] * pa_z[i];

        tr_z_zzzz_xzzzz[i] = 3.0 * tr_z_zz_xzzzz[i] * fe_0 + 4.0 * tr_z_zzz_xzzz[i] * fe_0 + ts_zzz_xzzzz[i] * fe_0 + tr_z_zzz_xzzzz[i] * pa_z[i];

        tr_z_zzzz_yyyyy[i] = 3.0 * tr_z_zz_yyyyy[i] * fe_0 + ts_zzz_yyyyy[i] * fe_0 + tr_z_zzz_yyyyy[i] * pa_z[i];

        tr_z_zzzz_yyyyz[i] = 3.0 * tr_z_zz_yyyyz[i] * fe_0 + tr_z_zzz_yyyy[i] * fe_0 + ts_zzz_yyyyz[i] * fe_0 + tr_z_zzz_yyyyz[i] * pa_z[i];

        tr_z_zzzz_yyyzz[i] = 3.0 * tr_z_zz_yyyzz[i] * fe_0 + 2.0 * tr_z_zzz_yyyz[i] * fe_0 + ts_zzz_yyyzz[i] * fe_0 + tr_z_zzz_yyyzz[i] * pa_z[i];

        tr_z_zzzz_yyzzz[i] = 3.0 * tr_z_zz_yyzzz[i] * fe_0 + 3.0 * tr_z_zzz_yyzz[i] * fe_0 + ts_zzz_yyzzz[i] * fe_0 + tr_z_zzz_yyzzz[i] * pa_z[i];

        tr_z_zzzz_yzzzz[i] = 3.0 * tr_z_zz_yzzzz[i] * fe_0 + 4.0 * tr_z_zzz_yzzz[i] * fe_0 + ts_zzz_yzzzz[i] * fe_0 + tr_z_zzz_yzzzz[i] * pa_z[i];

        tr_z_zzzz_zzzzz[i] = 3.0 * tr_z_zz_zzzzz[i] * fe_0 + 5.0 * tr_z_zzz_zzzz[i] * fe_0 + ts_zzz_zzzzz[i] * fe_0 + tr_z_zzz_zzzzz[i] * pa_z[i];
    }

}

} // diprec namespace

