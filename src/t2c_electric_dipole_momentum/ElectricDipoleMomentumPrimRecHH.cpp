#include "ElectricDipoleMomentumPrimRecHH.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_hh(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_hh,
                                      const size_t              idx_dip_fh,
                                      const size_t              idx_dip_gg,
                                      const size_t              idx_ovl_gh,
                                      const size_t              idx_dip_gh,
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

    auto tr_x_xxz_zzzzz = pbuffer.data(idx_dip_fh + 62);

    auto tr_x_xyy_xxxxx = pbuffer.data(idx_dip_fh + 63);

    auto tr_x_xyy_xxxxy = pbuffer.data(idx_dip_fh + 64);

    auto tr_x_xyy_xxxxz = pbuffer.data(idx_dip_fh + 65);

    auto tr_x_xyy_xxxyy = pbuffer.data(idx_dip_fh + 66);

    auto tr_x_xyy_xxxzz = pbuffer.data(idx_dip_fh + 68);

    auto tr_x_xyy_xxyyy = pbuffer.data(idx_dip_fh + 69);

    auto tr_x_xyy_xxzzz = pbuffer.data(idx_dip_fh + 72);

    auto tr_x_xyy_xyyyy = pbuffer.data(idx_dip_fh + 73);

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

    auto tr_x_xzz_xxxzz = pbuffer.data(idx_dip_fh + 110);

    auto tr_x_xzz_xxyyy = pbuffer.data(idx_dip_fh + 111);

    auto tr_x_xzz_xxzzz = pbuffer.data(idx_dip_fh + 114);

    auto tr_x_xzz_xyyyy = pbuffer.data(idx_dip_fh + 115);

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

    auto tr_y_xxy_zzzzz = pbuffer.data(idx_dip_fh + 251);

    auto tr_y_xxz_xxxxx = pbuffer.data(idx_dip_fh + 252);

    auto tr_y_xxz_xxxxy = pbuffer.data(idx_dip_fh + 253);

    auto tr_y_xxz_xxxyy = pbuffer.data(idx_dip_fh + 255);

    auto tr_y_xxz_xxyyy = pbuffer.data(idx_dip_fh + 258);

    auto tr_y_xxz_xyyyy = pbuffer.data(idx_dip_fh + 262);

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

    auto tr_y_yyz_xxxyy = pbuffer.data(idx_dip_fh + 360);

    auto tr_y_yyz_xxxyz = pbuffer.data(idx_dip_fh + 361);

    auto tr_y_yyz_xxyyy = pbuffer.data(idx_dip_fh + 363);

    auto tr_y_yyz_xxyyz = pbuffer.data(idx_dip_fh + 364);

    auto tr_y_yyz_xxyzz = pbuffer.data(idx_dip_fh + 365);

    auto tr_y_yyz_xyyyy = pbuffer.data(idx_dip_fh + 367);

    auto tr_y_yyz_xyyyz = pbuffer.data(idx_dip_fh + 368);

    auto tr_y_yyz_xyyzz = pbuffer.data(idx_dip_fh + 369);

    auto tr_y_yyz_xyzzz = pbuffer.data(idx_dip_fh + 370);

    auto tr_y_yyz_yyyyy = pbuffer.data(idx_dip_fh + 372);

    auto tr_y_yyz_yyyyz = pbuffer.data(idx_dip_fh + 373);

    auto tr_y_yyz_yyyzz = pbuffer.data(idx_dip_fh + 374);

    auto tr_y_yyz_yyzzz = pbuffer.data(idx_dip_fh + 375);

    auto tr_y_yyz_yzzzz = pbuffer.data(idx_dip_fh + 376);

    auto tr_y_yyz_zzzzz = pbuffer.data(idx_dip_fh + 377);

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

    auto tr_z_xxz_xxxxz = pbuffer.data(idx_dip_fh + 464);

    auto tr_z_xxz_xxxyz = pbuffer.data(idx_dip_fh + 466);

    auto tr_z_xxz_xxxzz = pbuffer.data(idx_dip_fh + 467);

    auto tr_z_xxz_xxyyz = pbuffer.data(idx_dip_fh + 469);

    auto tr_z_xxz_xxyzz = pbuffer.data(idx_dip_fh + 470);

    auto tr_z_xxz_xxzzz = pbuffer.data(idx_dip_fh + 471);

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

    auto tr_z_yyz_xxxxz = pbuffer.data(idx_dip_fh + 569);

    auto tr_z_yyz_xxxyz = pbuffer.data(idx_dip_fh + 571);

    auto tr_z_yyz_xxxzz = pbuffer.data(idx_dip_fh + 572);

    auto tr_z_yyz_xxyyz = pbuffer.data(idx_dip_fh + 574);

    auto tr_z_yyz_xxyzz = pbuffer.data(idx_dip_fh + 575);

    auto tr_z_yyz_xxzzz = pbuffer.data(idx_dip_fh + 576);

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

    auto tr_x_xxxz_yyyz = pbuffer.data(idx_dip_gg + 41);

    auto tr_x_xxxz_yyzz = pbuffer.data(idx_dip_gg + 42);

    auto tr_x_xxxz_yzzz = pbuffer.data(idx_dip_gg + 43);

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

    auto tr_x_xyyy_xxxy = pbuffer.data(idx_dip_gg + 91);

    auto tr_x_xyyy_xxyy = pbuffer.data(idx_dip_gg + 93);

    auto tr_x_xyyy_xxyz = pbuffer.data(idx_dip_gg + 94);

    auto tr_x_xyyy_xyyy = pbuffer.data(idx_dip_gg + 96);

    auto tr_x_xyyy_xyyz = pbuffer.data(idx_dip_gg + 97);

    auto tr_x_xyyy_xyzz = pbuffer.data(idx_dip_gg + 98);

    auto tr_x_xzzz_xxxx = pbuffer.data(idx_dip_gg + 135);

    auto tr_x_xzzz_xxxy = pbuffer.data(idx_dip_gg + 136);

    auto tr_x_xzzz_xxxz = pbuffer.data(idx_dip_gg + 137);

    auto tr_x_xzzz_xxyy = pbuffer.data(idx_dip_gg + 138);

    auto tr_x_xzzz_xxyz = pbuffer.data(idx_dip_gg + 139);

    auto tr_x_xzzz_xxzz = pbuffer.data(idx_dip_gg + 140);

    auto tr_x_xzzz_xyyy = pbuffer.data(idx_dip_gg + 141);

    auto tr_x_xzzz_xyyz = pbuffer.data(idx_dip_gg + 142);

    auto tr_x_xzzz_xyzz = pbuffer.data(idx_dip_gg + 143);

    auto tr_x_xzzz_xzzz = pbuffer.data(idx_dip_gg + 144);

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

    auto tr_x_yyzz_xxxz = pbuffer.data(idx_dip_gg + 182);

    auto tr_x_yyzz_xxyz = pbuffer.data(idx_dip_gg + 184);

    auto tr_x_yyzz_xxzz = pbuffer.data(idx_dip_gg + 185);

    auto tr_x_yyzz_xyyz = pbuffer.data(idx_dip_gg + 187);

    auto tr_x_yyzz_xyzz = pbuffer.data(idx_dip_gg + 188);

    auto tr_x_yyzz_xzzz = pbuffer.data(idx_dip_gg + 189);

    auto tr_x_yyzz_yyyz = pbuffer.data(idx_dip_gg + 191);

    auto tr_x_yyzz_yyzz = pbuffer.data(idx_dip_gg + 192);

    auto tr_x_yyzz_yzzz = pbuffer.data(idx_dip_gg + 193);

    auto tr_x_yyzz_zzzz = pbuffer.data(idx_dip_gg + 194);

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

    auto tr_y_xxzz_xxxz = pbuffer.data(idx_dip_gg + 302);

    auto tr_y_xxzz_xxyz = pbuffer.data(idx_dip_gg + 304);

    auto tr_y_xxzz_xxzz = pbuffer.data(idx_dip_gg + 305);

    auto tr_y_xxzz_xyyz = pbuffer.data(idx_dip_gg + 307);

    auto tr_y_xxzz_xyzz = pbuffer.data(idx_dip_gg + 308);

    auto tr_y_xxzz_xzzz = pbuffer.data(idx_dip_gg + 309);

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

    auto tr_y_xyzz_xxyz = pbuffer.data(idx_dip_gg + 349);

    auto tr_y_xyzz_xyyz = pbuffer.data(idx_dip_gg + 352);

    auto tr_y_xyzz_xyzz = pbuffer.data(idx_dip_gg + 353);

    auto tr_y_xyzz_yyyz = pbuffer.data(idx_dip_gg + 356);

    auto tr_y_xyzz_yyzz = pbuffer.data(idx_dip_gg + 357);

    auto tr_y_xyzz_yzzz = pbuffer.data(idx_dip_gg + 358);

    auto tr_y_xzzz_xxxz = pbuffer.data(idx_dip_gg + 362);

    auto tr_y_xzzz_xxyz = pbuffer.data(idx_dip_gg + 364);

    auto tr_y_xzzz_xxzz = pbuffer.data(idx_dip_gg + 365);

    auto tr_y_xzzz_xyyz = pbuffer.data(idx_dip_gg + 367);

    auto tr_y_xzzz_xyzz = pbuffer.data(idx_dip_gg + 368);

    auto tr_y_xzzz_xzzz = pbuffer.data(idx_dip_gg + 369);

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

    auto tr_y_yyyz_xxxy = pbuffer.data(idx_dip_gg + 391);

    auto tr_y_yyyz_xxxz = pbuffer.data(idx_dip_gg + 392);

    auto tr_y_yyyz_xxyy = pbuffer.data(idx_dip_gg + 393);

    auto tr_y_yyyz_xxyz = pbuffer.data(idx_dip_gg + 394);

    auto tr_y_yyyz_xxzz = pbuffer.data(idx_dip_gg + 395);

    auto tr_y_yyyz_xyyy = pbuffer.data(idx_dip_gg + 396);

    auto tr_y_yyyz_xyyz = pbuffer.data(idx_dip_gg + 397);

    auto tr_y_yyyz_xyzz = pbuffer.data(idx_dip_gg + 398);

    auto tr_y_yyyz_xzzz = pbuffer.data(idx_dip_gg + 399);

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

    auto tr_y_yzzz_xxxx = pbuffer.data(idx_dip_gg + 420);

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

    auto tr_z_xxxz_xxxx = pbuffer.data(idx_dip_gg + 480);

    auto tr_z_xxxz_xxxy = pbuffer.data(idx_dip_gg + 481);

    auto tr_z_xxxz_xxxz = pbuffer.data(idx_dip_gg + 482);

    auto tr_z_xxxz_xxyy = pbuffer.data(idx_dip_gg + 483);

    auto tr_z_xxxz_xxyz = pbuffer.data(idx_dip_gg + 484);

    auto tr_z_xxxz_xxzz = pbuffer.data(idx_dip_gg + 485);

    auto tr_z_xxxz_xyyy = pbuffer.data(idx_dip_gg + 486);

    auto tr_z_xxxz_xyyz = pbuffer.data(idx_dip_gg + 487);

    auto tr_z_xxxz_xyzz = pbuffer.data(idx_dip_gg + 488);

    auto tr_z_xxxz_xzzz = pbuffer.data(idx_dip_gg + 489);

    auto tr_z_xxxz_yyyz = pbuffer.data(idx_dip_gg + 491);

    auto tr_z_xxxz_yyzz = pbuffer.data(idx_dip_gg + 492);

    auto tr_z_xxxz_yzzz = pbuffer.data(idx_dip_gg + 493);

    auto tr_z_xxxz_zzzz = pbuffer.data(idx_dip_gg + 494);

    auto tr_z_xxyy_xxxy = pbuffer.data(idx_dip_gg + 496);

    auto tr_z_xxyy_xxyy = pbuffer.data(idx_dip_gg + 498);

    auto tr_z_xxyy_xxyz = pbuffer.data(idx_dip_gg + 499);

    auto tr_z_xxyy_xyyy = pbuffer.data(idx_dip_gg + 501);

    auto tr_z_xxyy_xyyz = pbuffer.data(idx_dip_gg + 502);

    auto tr_z_xxyy_xyzz = pbuffer.data(idx_dip_gg + 503);

    auto tr_z_xxyy_yyyy = pbuffer.data(idx_dip_gg + 505);

    auto tr_z_xxyy_yyyz = pbuffer.data(idx_dip_gg + 506);

    auto tr_z_xxyy_yyzz = pbuffer.data(idx_dip_gg + 507);

    auto tr_z_xxyy_yzzz = pbuffer.data(idx_dip_gg + 508);

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

    auto tr_z_xyyz_xxyz = pbuffer.data(idx_dip_gg + 559);

    auto tr_z_xyyz_xyyz = pbuffer.data(idx_dip_gg + 562);

    auto tr_z_xyyz_xyzz = pbuffer.data(idx_dip_gg + 563);

    auto tr_z_xyyz_yyyz = pbuffer.data(idx_dip_gg + 566);

    auto tr_z_xyyz_yyzz = pbuffer.data(idx_dip_gg + 567);

    auto tr_z_xyyz_yzzz = pbuffer.data(idx_dip_gg + 568);

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

    auto tr_z_yyyz_xxxy = pbuffer.data(idx_dip_gg + 616);

    auto tr_z_yyyz_xxxz = pbuffer.data(idx_dip_gg + 617);

    auto tr_z_yyyz_xxyy = pbuffer.data(idx_dip_gg + 618);

    auto tr_z_yyyz_xxyz = pbuffer.data(idx_dip_gg + 619);

    auto tr_z_yyyz_xxzz = pbuffer.data(idx_dip_gg + 620);

    auto tr_z_yyyz_xyyy = pbuffer.data(idx_dip_gg + 621);

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

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_ovl_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_ovl_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_ovl_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_ovl_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_ovl_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_ovl_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_ovl_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_ovl_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_ovl_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_ovl_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_ovl_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_ovl_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_ovl_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_ovl_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_ovl_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_ovl_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_ovl_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_ovl_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_ovl_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_ovl_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_ovl_gh + 20);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_ovl_gh + 44);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_ovl_gh + 47);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_ovl_gh + 51);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_ovl_gh + 56);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_ovl_gh + 64);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_ovl_gh + 66);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_ovl_gh + 69);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_ovl_gh + 73);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_ovl_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_ovl_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_ovl_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_ovl_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_ovl_gh + 82);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_ovl_gh + 105);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_ovl_gh + 107);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_ovl_gh + 110);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_ovl_gh + 114);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_ovl_gh + 119);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_ovl_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_ovl_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_ovl_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_ovl_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_ovl_gh + 125);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_ovl_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_ovl_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_ovl_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_ovl_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_ovl_gh + 145);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_ovl_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_ovl_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_ovl_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_ovl_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_ovl_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_ovl_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_ovl_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_ovl_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_ovl_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_ovl_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_ovl_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_ovl_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_ovl_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_ovl_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_ovl_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_ovl_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_ovl_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_ovl_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_ovl_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_ovl_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_ovl_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_ovl_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_ovl_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_ovl_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_ovl_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_ovl_gh + 230);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_ovl_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_ovl_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_ovl_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_ovl_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_ovl_gh + 251);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_ovl_gh + 254);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_ovl_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_ovl_gh + 257);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_ovl_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_ovl_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_ovl_gh + 261);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_ovl_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_ovl_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_ovl_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_ovl_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_ovl_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_ovl_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_ovl_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_ovl_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_ovl_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_ovl_gh + 272);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_ovl_gh + 275);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_ovl_gh + 278);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_ovl_gh + 282);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_ovl_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_ovl_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_ovl_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_ovl_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_ovl_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_ovl_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_ovl_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_ovl_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_ovl_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_ovl_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_ovl_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_ovl_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_ovl_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_ovl_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_ovl_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_ovl_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_ovl_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_ovl_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_ovl_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_ovl_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_ovl_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_ovl_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_ovl_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_ovl_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_ovl_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_ovl_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_ovl_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_ovl_gh + 314);

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

    auto tr_x_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 36);

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

    auto tr_x_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 58);

    auto tr_x_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 59);

    auto tr_x_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 60);

    auto tr_x_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 61);

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

    auto tr_x_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 247);

    auto tr_x_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 248);

    auto tr_x_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 249);

    auto tr_x_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 250);

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

    auto tr_x_yzzz_yyyyy = pbuffer.data(idx_dip_gh + 288);

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

    auto tr_y_xxxy_xxxxx = pbuffer.data(idx_dip_gh + 336);

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

    auto tr_y_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 359);

    auto tr_y_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 360);

    auto tr_y_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 362);

    auto tr_y_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 363);

    auto tr_y_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 366);

    auto tr_y_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 367);

    auto tr_y_xxxz_xzzzz = pbuffer.data(idx_dip_gh + 371);

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

    // Set up 0-21 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_x_xxx_xxxxx,   \
                             tr_x_xxx_xxxxy,   \
                             tr_x_xxx_xxxxz,   \
                             tr_x_xxx_xxxyy,   \
                             tr_x_xxx_xxxyz,   \
                             tr_x_xxx_xxxzz,   \
                             tr_x_xxx_xxyyy,   \
                             tr_x_xxx_xxyyz,   \
                             tr_x_xxx_xxyzz,   \
                             tr_x_xxx_xxzzz,   \
                             tr_x_xxx_xyyyy,   \
                             tr_x_xxx_xyyyz,   \
                             tr_x_xxx_xyyzz,   \
                             tr_x_xxx_xyzzz,   \
                             tr_x_xxx_xzzzz,   \
                             tr_x_xxx_yyyyy,   \
                             tr_x_xxx_yyyyz,   \
                             tr_x_xxx_yyyzz,   \
                             tr_x_xxx_yyzzz,   \
                             tr_x_xxx_yzzzz,   \
                             tr_x_xxx_zzzzz,   \
                             tr_x_xxxx_xxxx,   \
                             tr_x_xxxx_xxxxx,  \
                             tr_x_xxxx_xxxxy,  \
                             tr_x_xxxx_xxxxz,  \
                             tr_x_xxxx_xxxy,   \
                             tr_x_xxxx_xxxyy,  \
                             tr_x_xxxx_xxxyz,  \
                             tr_x_xxxx_xxxz,   \
                             tr_x_xxxx_xxxzz,  \
                             tr_x_xxxx_xxyy,   \
                             tr_x_xxxx_xxyyy,  \
                             tr_x_xxxx_xxyyz,  \
                             tr_x_xxxx_xxyz,   \
                             tr_x_xxxx_xxyzz,  \
                             tr_x_xxxx_xxzz,   \
                             tr_x_xxxx_xxzzz,  \
                             tr_x_xxxx_xyyy,   \
                             tr_x_xxxx_xyyyy,  \
                             tr_x_xxxx_xyyyz,  \
                             tr_x_xxxx_xyyz,   \
                             tr_x_xxxx_xyyzz,  \
                             tr_x_xxxx_xyzz,   \
                             tr_x_xxxx_xyzzz,  \
                             tr_x_xxxx_xzzz,   \
                             tr_x_xxxx_xzzzz,  \
                             tr_x_xxxx_yyyy,   \
                             tr_x_xxxx_yyyyy,  \
                             tr_x_xxxx_yyyyz,  \
                             tr_x_xxxx_yyyz,   \
                             tr_x_xxxx_yyyzz,  \
                             tr_x_xxxx_yyzz,   \
                             tr_x_xxxx_yyzzz,  \
                             tr_x_xxxx_yzzz,   \
                             tr_x_xxxx_yzzzz,  \
                             tr_x_xxxx_zzzz,   \
                             tr_x_xxxx_zzzzz,  \
                             tr_x_xxxxx_xxxxx, \
                             tr_x_xxxxx_xxxxy, \
                             tr_x_xxxxx_xxxxz, \
                             tr_x_xxxxx_xxxyy, \
                             tr_x_xxxxx_xxxyz, \
                             tr_x_xxxxx_xxxzz, \
                             tr_x_xxxxx_xxyyy, \
                             tr_x_xxxxx_xxyyz, \
                             tr_x_xxxxx_xxyzz, \
                             tr_x_xxxxx_xxzzz, \
                             tr_x_xxxxx_xyyyy, \
                             tr_x_xxxxx_xyyyz, \
                             tr_x_xxxxx_xyyzz, \
                             tr_x_xxxxx_xyzzz, \
                             tr_x_xxxxx_xzzzz, \
                             tr_x_xxxxx_yyyyy, \
                             tr_x_xxxxx_yyyyz, \
                             tr_x_xxxxx_yyyzz, \
                             tr_x_xxxxx_yyzzz, \
                             tr_x_xxxxx_yzzzz, \
                             tr_x_xxxxx_zzzzz, \
                             ts_xxxx_xxxxx,    \
                             ts_xxxx_xxxxy,    \
                             ts_xxxx_xxxxz,    \
                             ts_xxxx_xxxyy,    \
                             ts_xxxx_xxxyz,    \
                             ts_xxxx_xxxzz,    \
                             ts_xxxx_xxyyy,    \
                             ts_xxxx_xxyyz,    \
                             ts_xxxx_xxyzz,    \
                             ts_xxxx_xxzzz,    \
                             ts_xxxx_xyyyy,    \
                             ts_xxxx_xyyyz,    \
                             ts_xxxx_xyyzz,    \
                             ts_xxxx_xyzzz,    \
                             ts_xxxx_xzzzz,    \
                             ts_xxxx_yyyyy,    \
                             ts_xxxx_yyyyz,    \
                             ts_xxxx_yyyzz,    \
                             ts_xxxx_yyzzz,    \
                             ts_xxxx_yzzzz,    \
                             ts_xxxx_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_xxxxx[i] =
            4.0 * tr_x_xxx_xxxxx[i] * fe_0 + 5.0 * tr_x_xxxx_xxxx[i] * fe_0 + ts_xxxx_xxxxx[i] * fe_0 + tr_x_xxxx_xxxxx[i] * pa_x[i];

        tr_x_xxxxx_xxxxy[i] =
            4.0 * tr_x_xxx_xxxxy[i] * fe_0 + 4.0 * tr_x_xxxx_xxxy[i] * fe_0 + ts_xxxx_xxxxy[i] * fe_0 + tr_x_xxxx_xxxxy[i] * pa_x[i];

        tr_x_xxxxx_xxxxz[i] =
            4.0 * tr_x_xxx_xxxxz[i] * fe_0 + 4.0 * tr_x_xxxx_xxxz[i] * fe_0 + ts_xxxx_xxxxz[i] * fe_0 + tr_x_xxxx_xxxxz[i] * pa_x[i];

        tr_x_xxxxx_xxxyy[i] =
            4.0 * tr_x_xxx_xxxyy[i] * fe_0 + 3.0 * tr_x_xxxx_xxyy[i] * fe_0 + ts_xxxx_xxxyy[i] * fe_0 + tr_x_xxxx_xxxyy[i] * pa_x[i];

        tr_x_xxxxx_xxxyz[i] =
            4.0 * tr_x_xxx_xxxyz[i] * fe_0 + 3.0 * tr_x_xxxx_xxyz[i] * fe_0 + ts_xxxx_xxxyz[i] * fe_0 + tr_x_xxxx_xxxyz[i] * pa_x[i];

        tr_x_xxxxx_xxxzz[i] =
            4.0 * tr_x_xxx_xxxzz[i] * fe_0 + 3.0 * tr_x_xxxx_xxzz[i] * fe_0 + ts_xxxx_xxxzz[i] * fe_0 + tr_x_xxxx_xxxzz[i] * pa_x[i];

        tr_x_xxxxx_xxyyy[i] =
            4.0 * tr_x_xxx_xxyyy[i] * fe_0 + 2.0 * tr_x_xxxx_xyyy[i] * fe_0 + ts_xxxx_xxyyy[i] * fe_0 + tr_x_xxxx_xxyyy[i] * pa_x[i];

        tr_x_xxxxx_xxyyz[i] =
            4.0 * tr_x_xxx_xxyyz[i] * fe_0 + 2.0 * tr_x_xxxx_xyyz[i] * fe_0 + ts_xxxx_xxyyz[i] * fe_0 + tr_x_xxxx_xxyyz[i] * pa_x[i];

        tr_x_xxxxx_xxyzz[i] =
            4.0 * tr_x_xxx_xxyzz[i] * fe_0 + 2.0 * tr_x_xxxx_xyzz[i] * fe_0 + ts_xxxx_xxyzz[i] * fe_0 + tr_x_xxxx_xxyzz[i] * pa_x[i];

        tr_x_xxxxx_xxzzz[i] =
            4.0 * tr_x_xxx_xxzzz[i] * fe_0 + 2.0 * tr_x_xxxx_xzzz[i] * fe_0 + ts_xxxx_xxzzz[i] * fe_0 + tr_x_xxxx_xxzzz[i] * pa_x[i];

        tr_x_xxxxx_xyyyy[i] = 4.0 * tr_x_xxx_xyyyy[i] * fe_0 + tr_x_xxxx_yyyy[i] * fe_0 + ts_xxxx_xyyyy[i] * fe_0 + tr_x_xxxx_xyyyy[i] * pa_x[i];

        tr_x_xxxxx_xyyyz[i] = 4.0 * tr_x_xxx_xyyyz[i] * fe_0 + tr_x_xxxx_yyyz[i] * fe_0 + ts_xxxx_xyyyz[i] * fe_0 + tr_x_xxxx_xyyyz[i] * pa_x[i];

        tr_x_xxxxx_xyyzz[i] = 4.0 * tr_x_xxx_xyyzz[i] * fe_0 + tr_x_xxxx_yyzz[i] * fe_0 + ts_xxxx_xyyzz[i] * fe_0 + tr_x_xxxx_xyyzz[i] * pa_x[i];

        tr_x_xxxxx_xyzzz[i] = 4.0 * tr_x_xxx_xyzzz[i] * fe_0 + tr_x_xxxx_yzzz[i] * fe_0 + ts_xxxx_xyzzz[i] * fe_0 + tr_x_xxxx_xyzzz[i] * pa_x[i];

        tr_x_xxxxx_xzzzz[i] = 4.0 * tr_x_xxx_xzzzz[i] * fe_0 + tr_x_xxxx_zzzz[i] * fe_0 + ts_xxxx_xzzzz[i] * fe_0 + tr_x_xxxx_xzzzz[i] * pa_x[i];

        tr_x_xxxxx_yyyyy[i] = 4.0 * tr_x_xxx_yyyyy[i] * fe_0 + ts_xxxx_yyyyy[i] * fe_0 + tr_x_xxxx_yyyyy[i] * pa_x[i];

        tr_x_xxxxx_yyyyz[i] = 4.0 * tr_x_xxx_yyyyz[i] * fe_0 + ts_xxxx_yyyyz[i] * fe_0 + tr_x_xxxx_yyyyz[i] * pa_x[i];

        tr_x_xxxxx_yyyzz[i] = 4.0 * tr_x_xxx_yyyzz[i] * fe_0 + ts_xxxx_yyyzz[i] * fe_0 + tr_x_xxxx_yyyzz[i] * pa_x[i];

        tr_x_xxxxx_yyzzz[i] = 4.0 * tr_x_xxx_yyzzz[i] * fe_0 + ts_xxxx_yyzzz[i] * fe_0 + tr_x_xxxx_yyzzz[i] * pa_x[i];

        tr_x_xxxxx_yzzzz[i] = 4.0 * tr_x_xxx_yzzzz[i] * fe_0 + ts_xxxx_yzzzz[i] * fe_0 + tr_x_xxxx_yzzzz[i] * pa_x[i];

        tr_x_xxxxx_zzzzz[i] = 4.0 * tr_x_xxx_zzzzz[i] * fe_0 + ts_xxxx_zzzzz[i] * fe_0 + tr_x_xxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : HH

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

    auto tr_x_xxxxy_yyyyz = pbuffer.data(idx_dip_hh + 37);

    auto tr_x_xxxxy_yyyzz = pbuffer.data(idx_dip_hh + 38);

    auto tr_x_xxxxy_yyzzz = pbuffer.data(idx_dip_hh + 39);

    auto tr_x_xxxxy_yzzzz = pbuffer.data(idx_dip_hh + 40);

    auto tr_x_xxxxy_zzzzz = pbuffer.data(idx_dip_hh + 41);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_xxxx_xxxx,   \
                             tr_x_xxxx_xxxxx,  \
                             tr_x_xxxx_xxxxy,  \
                             tr_x_xxxx_xxxxz,  \
                             tr_x_xxxx_xxxy,   \
                             tr_x_xxxx_xxxyy,  \
                             tr_x_xxxx_xxxyz,  \
                             tr_x_xxxx_xxxz,   \
                             tr_x_xxxx_xxxzz,  \
                             tr_x_xxxx_xxyy,   \
                             tr_x_xxxx_xxyyy,  \
                             tr_x_xxxx_xxyyz,  \
                             tr_x_xxxx_xxyz,   \
                             tr_x_xxxx_xxyzz,  \
                             tr_x_xxxx_xxzz,   \
                             tr_x_xxxx_xxzzz,  \
                             tr_x_xxxx_xyyy,   \
                             tr_x_xxxx_xyyyy,  \
                             tr_x_xxxx_xyyyz,  \
                             tr_x_xxxx_xyyz,   \
                             tr_x_xxxx_xyyzz,  \
                             tr_x_xxxx_xyzz,   \
                             tr_x_xxxx_xyzzz,  \
                             tr_x_xxxx_xzzz,   \
                             tr_x_xxxx_xzzzz,  \
                             tr_x_xxxx_yyyy,   \
                             tr_x_xxxx_yyyyy,  \
                             tr_x_xxxx_yyyyz,  \
                             tr_x_xxxx_yyyz,   \
                             tr_x_xxxx_yyyzz,  \
                             tr_x_xxxx_yyzz,   \
                             tr_x_xxxx_yyzzz,  \
                             tr_x_xxxx_yzzz,   \
                             tr_x_xxxx_yzzzz,  \
                             tr_x_xxxx_zzzz,   \
                             tr_x_xxxx_zzzzz,  \
                             tr_x_xxxxy_xxxxx, \
                             tr_x_xxxxy_xxxxy, \
                             tr_x_xxxxy_xxxxz, \
                             tr_x_xxxxy_xxxyy, \
                             tr_x_xxxxy_xxxyz, \
                             tr_x_xxxxy_xxxzz, \
                             tr_x_xxxxy_xxyyy, \
                             tr_x_xxxxy_xxyyz, \
                             tr_x_xxxxy_xxyzz, \
                             tr_x_xxxxy_xxzzz, \
                             tr_x_xxxxy_xyyyy, \
                             tr_x_xxxxy_xyyyz, \
                             tr_x_xxxxy_xyyzz, \
                             tr_x_xxxxy_xyzzz, \
                             tr_x_xxxxy_xzzzz, \
                             tr_x_xxxxy_yyyyy, \
                             tr_x_xxxxy_yyyyz, \
                             tr_x_xxxxy_yyyzz, \
                             tr_x_xxxxy_yyzzz, \
                             tr_x_xxxxy_yzzzz, \
                             tr_x_xxxxy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxy_xxxxx[i] = tr_x_xxxx_xxxxx[i] * pa_y[i];

        tr_x_xxxxy_xxxxy[i] = tr_x_xxxx_xxxx[i] * fe_0 + tr_x_xxxx_xxxxy[i] * pa_y[i];

        tr_x_xxxxy_xxxxz[i] = tr_x_xxxx_xxxxz[i] * pa_y[i];

        tr_x_xxxxy_xxxyy[i] = 2.0 * tr_x_xxxx_xxxy[i] * fe_0 + tr_x_xxxx_xxxyy[i] * pa_y[i];

        tr_x_xxxxy_xxxyz[i] = tr_x_xxxx_xxxz[i] * fe_0 + tr_x_xxxx_xxxyz[i] * pa_y[i];

        tr_x_xxxxy_xxxzz[i] = tr_x_xxxx_xxxzz[i] * pa_y[i];

        tr_x_xxxxy_xxyyy[i] = 3.0 * tr_x_xxxx_xxyy[i] * fe_0 + tr_x_xxxx_xxyyy[i] * pa_y[i];

        tr_x_xxxxy_xxyyz[i] = 2.0 * tr_x_xxxx_xxyz[i] * fe_0 + tr_x_xxxx_xxyyz[i] * pa_y[i];

        tr_x_xxxxy_xxyzz[i] = tr_x_xxxx_xxzz[i] * fe_0 + tr_x_xxxx_xxyzz[i] * pa_y[i];

        tr_x_xxxxy_xxzzz[i] = tr_x_xxxx_xxzzz[i] * pa_y[i];

        tr_x_xxxxy_xyyyy[i] = 4.0 * tr_x_xxxx_xyyy[i] * fe_0 + tr_x_xxxx_xyyyy[i] * pa_y[i];

        tr_x_xxxxy_xyyyz[i] = 3.0 * tr_x_xxxx_xyyz[i] * fe_0 + tr_x_xxxx_xyyyz[i] * pa_y[i];

        tr_x_xxxxy_xyyzz[i] = 2.0 * tr_x_xxxx_xyzz[i] * fe_0 + tr_x_xxxx_xyyzz[i] * pa_y[i];

        tr_x_xxxxy_xyzzz[i] = tr_x_xxxx_xzzz[i] * fe_0 + tr_x_xxxx_xyzzz[i] * pa_y[i];

        tr_x_xxxxy_xzzzz[i] = tr_x_xxxx_xzzzz[i] * pa_y[i];

        tr_x_xxxxy_yyyyy[i] = 5.0 * tr_x_xxxx_yyyy[i] * fe_0 + tr_x_xxxx_yyyyy[i] * pa_y[i];

        tr_x_xxxxy_yyyyz[i] = 4.0 * tr_x_xxxx_yyyz[i] * fe_0 + tr_x_xxxx_yyyyz[i] * pa_y[i];

        tr_x_xxxxy_yyyzz[i] = 3.0 * tr_x_xxxx_yyzz[i] * fe_0 + tr_x_xxxx_yyyzz[i] * pa_y[i];

        tr_x_xxxxy_yyzzz[i] = 2.0 * tr_x_xxxx_yzzz[i] * fe_0 + tr_x_xxxx_yyzzz[i] * pa_y[i];

        tr_x_xxxxy_yzzzz[i] = tr_x_xxxx_zzzz[i] * fe_0 + tr_x_xxxx_yzzzz[i] * pa_y[i];

        tr_x_xxxxy_zzzzz[i] = tr_x_xxxx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_z,                 \
                             tr_x_xxxx_xxxx,   \
                             tr_x_xxxx_xxxxx,  \
                             tr_x_xxxx_xxxxy,  \
                             tr_x_xxxx_xxxxz,  \
                             tr_x_xxxx_xxxy,   \
                             tr_x_xxxx_xxxyy,  \
                             tr_x_xxxx_xxxyz,  \
                             tr_x_xxxx_xxxz,   \
                             tr_x_xxxx_xxxzz,  \
                             tr_x_xxxx_xxyy,   \
                             tr_x_xxxx_xxyyy,  \
                             tr_x_xxxx_xxyyz,  \
                             tr_x_xxxx_xxyz,   \
                             tr_x_xxxx_xxyzz,  \
                             tr_x_xxxx_xxzz,   \
                             tr_x_xxxx_xxzzz,  \
                             tr_x_xxxx_xyyy,   \
                             tr_x_xxxx_xyyyy,  \
                             tr_x_xxxx_xyyyz,  \
                             tr_x_xxxx_xyyz,   \
                             tr_x_xxxx_xyyzz,  \
                             tr_x_xxxx_xyzz,   \
                             tr_x_xxxx_xyzzz,  \
                             tr_x_xxxx_xzzz,   \
                             tr_x_xxxx_xzzzz,  \
                             tr_x_xxxx_yyyy,   \
                             tr_x_xxxx_yyyyy,  \
                             tr_x_xxxx_yyyyz,  \
                             tr_x_xxxx_yyyz,   \
                             tr_x_xxxx_yyyzz,  \
                             tr_x_xxxx_yyzz,   \
                             tr_x_xxxx_yyzzz,  \
                             tr_x_xxxx_yzzz,   \
                             tr_x_xxxx_yzzzz,  \
                             tr_x_xxxx_zzzz,   \
                             tr_x_xxxx_zzzzz,  \
                             tr_x_xxxxz_xxxxx, \
                             tr_x_xxxxz_xxxxy, \
                             tr_x_xxxxz_xxxxz, \
                             tr_x_xxxxz_xxxyy, \
                             tr_x_xxxxz_xxxyz, \
                             tr_x_xxxxz_xxxzz, \
                             tr_x_xxxxz_xxyyy, \
                             tr_x_xxxxz_xxyyz, \
                             tr_x_xxxxz_xxyzz, \
                             tr_x_xxxxz_xxzzz, \
                             tr_x_xxxxz_xyyyy, \
                             tr_x_xxxxz_xyyyz, \
                             tr_x_xxxxz_xyyzz, \
                             tr_x_xxxxz_xyzzz, \
                             tr_x_xxxxz_xzzzz, \
                             tr_x_xxxxz_yyyyy, \
                             tr_x_xxxxz_yyyyz, \
                             tr_x_xxxxz_yyyzz, \
                             tr_x_xxxxz_yyzzz, \
                             tr_x_xxxxz_yzzzz, \
                             tr_x_xxxxz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxz_xxxxx[i] = tr_x_xxxx_xxxxx[i] * pa_z[i];

        tr_x_xxxxz_xxxxy[i] = tr_x_xxxx_xxxxy[i] * pa_z[i];

        tr_x_xxxxz_xxxxz[i] = tr_x_xxxx_xxxx[i] * fe_0 + tr_x_xxxx_xxxxz[i] * pa_z[i];

        tr_x_xxxxz_xxxyy[i] = tr_x_xxxx_xxxyy[i] * pa_z[i];

        tr_x_xxxxz_xxxyz[i] = tr_x_xxxx_xxxy[i] * fe_0 + tr_x_xxxx_xxxyz[i] * pa_z[i];

        tr_x_xxxxz_xxxzz[i] = 2.0 * tr_x_xxxx_xxxz[i] * fe_0 + tr_x_xxxx_xxxzz[i] * pa_z[i];

        tr_x_xxxxz_xxyyy[i] = tr_x_xxxx_xxyyy[i] * pa_z[i];

        tr_x_xxxxz_xxyyz[i] = tr_x_xxxx_xxyy[i] * fe_0 + tr_x_xxxx_xxyyz[i] * pa_z[i];

        tr_x_xxxxz_xxyzz[i] = 2.0 * tr_x_xxxx_xxyz[i] * fe_0 + tr_x_xxxx_xxyzz[i] * pa_z[i];

        tr_x_xxxxz_xxzzz[i] = 3.0 * tr_x_xxxx_xxzz[i] * fe_0 + tr_x_xxxx_xxzzz[i] * pa_z[i];

        tr_x_xxxxz_xyyyy[i] = tr_x_xxxx_xyyyy[i] * pa_z[i];

        tr_x_xxxxz_xyyyz[i] = tr_x_xxxx_xyyy[i] * fe_0 + tr_x_xxxx_xyyyz[i] * pa_z[i];

        tr_x_xxxxz_xyyzz[i] = 2.0 * tr_x_xxxx_xyyz[i] * fe_0 + tr_x_xxxx_xyyzz[i] * pa_z[i];

        tr_x_xxxxz_xyzzz[i] = 3.0 * tr_x_xxxx_xyzz[i] * fe_0 + tr_x_xxxx_xyzzz[i] * pa_z[i];

        tr_x_xxxxz_xzzzz[i] = 4.0 * tr_x_xxxx_xzzz[i] * fe_0 + tr_x_xxxx_xzzzz[i] * pa_z[i];

        tr_x_xxxxz_yyyyy[i] = tr_x_xxxx_yyyyy[i] * pa_z[i];

        tr_x_xxxxz_yyyyz[i] = tr_x_xxxx_yyyy[i] * fe_0 + tr_x_xxxx_yyyyz[i] * pa_z[i];

        tr_x_xxxxz_yyyzz[i] = 2.0 * tr_x_xxxx_yyyz[i] * fe_0 + tr_x_xxxx_yyyzz[i] * pa_z[i];

        tr_x_xxxxz_yyzzz[i] = 3.0 * tr_x_xxxx_yyzz[i] * fe_0 + tr_x_xxxx_yyzzz[i] * pa_z[i];

        tr_x_xxxxz_yzzzz[i] = 4.0 * tr_x_xxxx_yzzz[i] * fe_0 + tr_x_xxxx_yzzzz[i] * pa_z[i];

        tr_x_xxxxz_zzzzz[i] = 5.0 * tr_x_xxxx_zzzz[i] * fe_0 + tr_x_xxxx_zzzzz[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xxx_xxxxx,   \
                             tr_x_xxx_xxxxy,   \
                             tr_x_xxx_xxxxz,   \
                             tr_x_xxx_xxxyy,   \
                             tr_x_xxx_xxxyz,   \
                             tr_x_xxx_xxxzz,   \
                             tr_x_xxx_xxyyy,   \
                             tr_x_xxx_xxyyz,   \
                             tr_x_xxx_xxyzz,   \
                             tr_x_xxx_xxzzz,   \
                             tr_x_xxx_xyyyy,   \
                             tr_x_xxx_xyyyz,   \
                             tr_x_xxx_xyyzz,   \
                             tr_x_xxx_xyzzz,   \
                             tr_x_xxx_xzzzz,   \
                             tr_x_xxx_zzzzz,   \
                             tr_x_xxxy_xxxx,   \
                             tr_x_xxxy_xxxxx,  \
                             tr_x_xxxy_xxxxy,  \
                             tr_x_xxxy_xxxxz,  \
                             tr_x_xxxy_xxxy,   \
                             tr_x_xxxy_xxxyy,  \
                             tr_x_xxxy_xxxyz,  \
                             tr_x_xxxy_xxxz,   \
                             tr_x_xxxy_xxxzz,  \
                             tr_x_xxxy_xxyy,   \
                             tr_x_xxxy_xxyyy,  \
                             tr_x_xxxy_xxyyz,  \
                             tr_x_xxxy_xxyz,   \
                             tr_x_xxxy_xxyzz,  \
                             tr_x_xxxy_xxzz,   \
                             tr_x_xxxy_xxzzz,  \
                             tr_x_xxxy_xyyy,   \
                             tr_x_xxxy_xyyyy,  \
                             tr_x_xxxy_xyyyz,  \
                             tr_x_xxxy_xyyz,   \
                             tr_x_xxxy_xyyzz,  \
                             tr_x_xxxy_xyzz,   \
                             tr_x_xxxy_xyzzz,  \
                             tr_x_xxxy_xzzz,   \
                             tr_x_xxxy_xzzzz,  \
                             tr_x_xxxy_zzzzz,  \
                             tr_x_xxxyy_xxxxx, \
                             tr_x_xxxyy_xxxxy, \
                             tr_x_xxxyy_xxxxz, \
                             tr_x_xxxyy_xxxyy, \
                             tr_x_xxxyy_xxxyz, \
                             tr_x_xxxyy_xxxzz, \
                             tr_x_xxxyy_xxyyy, \
                             tr_x_xxxyy_xxyyz, \
                             tr_x_xxxyy_xxyzz, \
                             tr_x_xxxyy_xxzzz, \
                             tr_x_xxxyy_xyyyy, \
                             tr_x_xxxyy_xyyyz, \
                             tr_x_xxxyy_xyyzz, \
                             tr_x_xxxyy_xyzzz, \
                             tr_x_xxxyy_xzzzz, \
                             tr_x_xxxyy_yyyyy, \
                             tr_x_xxxyy_yyyyz, \
                             tr_x_xxxyy_yyyzz, \
                             tr_x_xxxyy_yyzzz, \
                             tr_x_xxxyy_yzzzz, \
                             tr_x_xxxyy_zzzzz, \
                             tr_x_xxyy_yyyyy,  \
                             tr_x_xxyy_yyyyz,  \
                             tr_x_xxyy_yyyzz,  \
                             tr_x_xxyy_yyzzz,  \
                             tr_x_xxyy_yzzzz,  \
                             tr_x_xyy_yyyyy,   \
                             tr_x_xyy_yyyyz,   \
                             tr_x_xyy_yyyzz,   \
                             tr_x_xyy_yyzzz,   \
                             tr_x_xyy_yzzzz,   \
                             ts_xxyy_yyyyy,    \
                             ts_xxyy_yyyyz,    \
                             ts_xxyy_yyyzz,    \
                             ts_xxyy_yyzzz,    \
                             ts_xxyy_yzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyy_xxxxx[i] = tr_x_xxx_xxxxx[i] * fe_0 + tr_x_xxxy_xxxxx[i] * pa_y[i];

        tr_x_xxxyy_xxxxy[i] = tr_x_xxx_xxxxy[i] * fe_0 + tr_x_xxxy_xxxx[i] * fe_0 + tr_x_xxxy_xxxxy[i] * pa_y[i];

        tr_x_xxxyy_xxxxz[i] = tr_x_xxx_xxxxz[i] * fe_0 + tr_x_xxxy_xxxxz[i] * pa_y[i];

        tr_x_xxxyy_xxxyy[i] = tr_x_xxx_xxxyy[i] * fe_0 + 2.0 * tr_x_xxxy_xxxy[i] * fe_0 + tr_x_xxxy_xxxyy[i] * pa_y[i];

        tr_x_xxxyy_xxxyz[i] = tr_x_xxx_xxxyz[i] * fe_0 + tr_x_xxxy_xxxz[i] * fe_0 + tr_x_xxxy_xxxyz[i] * pa_y[i];

        tr_x_xxxyy_xxxzz[i] = tr_x_xxx_xxxzz[i] * fe_0 + tr_x_xxxy_xxxzz[i] * pa_y[i];

        tr_x_xxxyy_xxyyy[i] = tr_x_xxx_xxyyy[i] * fe_0 + 3.0 * tr_x_xxxy_xxyy[i] * fe_0 + tr_x_xxxy_xxyyy[i] * pa_y[i];

        tr_x_xxxyy_xxyyz[i] = tr_x_xxx_xxyyz[i] * fe_0 + 2.0 * tr_x_xxxy_xxyz[i] * fe_0 + tr_x_xxxy_xxyyz[i] * pa_y[i];

        tr_x_xxxyy_xxyzz[i] = tr_x_xxx_xxyzz[i] * fe_0 + tr_x_xxxy_xxzz[i] * fe_0 + tr_x_xxxy_xxyzz[i] * pa_y[i];

        tr_x_xxxyy_xxzzz[i] = tr_x_xxx_xxzzz[i] * fe_0 + tr_x_xxxy_xxzzz[i] * pa_y[i];

        tr_x_xxxyy_xyyyy[i] = tr_x_xxx_xyyyy[i] * fe_0 + 4.0 * tr_x_xxxy_xyyy[i] * fe_0 + tr_x_xxxy_xyyyy[i] * pa_y[i];

        tr_x_xxxyy_xyyyz[i] = tr_x_xxx_xyyyz[i] * fe_0 + 3.0 * tr_x_xxxy_xyyz[i] * fe_0 + tr_x_xxxy_xyyyz[i] * pa_y[i];

        tr_x_xxxyy_xyyzz[i] = tr_x_xxx_xyyzz[i] * fe_0 + 2.0 * tr_x_xxxy_xyzz[i] * fe_0 + tr_x_xxxy_xyyzz[i] * pa_y[i];

        tr_x_xxxyy_xyzzz[i] = tr_x_xxx_xyzzz[i] * fe_0 + tr_x_xxxy_xzzz[i] * fe_0 + tr_x_xxxy_xyzzz[i] * pa_y[i];

        tr_x_xxxyy_xzzzz[i] = tr_x_xxx_xzzzz[i] * fe_0 + tr_x_xxxy_xzzzz[i] * pa_y[i];

        tr_x_xxxyy_yyyyy[i] = 2.0 * tr_x_xyy_yyyyy[i] * fe_0 + ts_xxyy_yyyyy[i] * fe_0 + tr_x_xxyy_yyyyy[i] * pa_x[i];

        tr_x_xxxyy_yyyyz[i] = 2.0 * tr_x_xyy_yyyyz[i] * fe_0 + ts_xxyy_yyyyz[i] * fe_0 + tr_x_xxyy_yyyyz[i] * pa_x[i];

        tr_x_xxxyy_yyyzz[i] = 2.0 * tr_x_xyy_yyyzz[i] * fe_0 + ts_xxyy_yyyzz[i] * fe_0 + tr_x_xxyy_yyyzz[i] * pa_x[i];

        tr_x_xxxyy_yyzzz[i] = 2.0 * tr_x_xyy_yyzzz[i] * fe_0 + ts_xxyy_yyzzz[i] * fe_0 + tr_x_xxyy_yyzzz[i] * pa_x[i];

        tr_x_xxxyy_yzzzz[i] = 2.0 * tr_x_xyy_yzzzz[i] * fe_0 + ts_xxyy_yzzzz[i] * fe_0 + tr_x_xxyy_yzzzz[i] * pa_x[i];

        tr_x_xxxyy_zzzzz[i] = tr_x_xxx_zzzzz[i] * fe_0 + tr_x_xxxy_zzzzz[i] * pa_y[i];
    }

    // Set up 84-105 components of targeted buffer : HH

    auto tr_x_xxxyz_xxxxx = pbuffer.data(idx_dip_hh + 84);

    auto tr_x_xxxyz_xxxxy = pbuffer.data(idx_dip_hh + 85);

    auto tr_x_xxxyz_xxxxz = pbuffer.data(idx_dip_hh + 86);

    auto tr_x_xxxyz_xxxyy = pbuffer.data(idx_dip_hh + 87);

    auto tr_x_xxxyz_xxxyz = pbuffer.data(idx_dip_hh + 88);

    auto tr_x_xxxyz_xxxzz = pbuffer.data(idx_dip_hh + 89);

    auto tr_x_xxxyz_xxyyy = pbuffer.data(idx_dip_hh + 90);

    auto tr_x_xxxyz_xxyyz = pbuffer.data(idx_dip_hh + 91);

    auto tr_x_xxxyz_xxyzz = pbuffer.data(idx_dip_hh + 92);

    auto tr_x_xxxyz_xxzzz = pbuffer.data(idx_dip_hh + 93);

    auto tr_x_xxxyz_xyyyy = pbuffer.data(idx_dip_hh + 94);

    auto tr_x_xxxyz_xyyyz = pbuffer.data(idx_dip_hh + 95);

    auto tr_x_xxxyz_xyyzz = pbuffer.data(idx_dip_hh + 96);

    auto tr_x_xxxyz_xyzzz = pbuffer.data(idx_dip_hh + 97);

    auto tr_x_xxxyz_xzzzz = pbuffer.data(idx_dip_hh + 98);

    auto tr_x_xxxyz_yyyyy = pbuffer.data(idx_dip_hh + 99);

    auto tr_x_xxxyz_yyyyz = pbuffer.data(idx_dip_hh + 100);

    auto tr_x_xxxyz_yyyzz = pbuffer.data(idx_dip_hh + 101);

    auto tr_x_xxxyz_yyzzz = pbuffer.data(idx_dip_hh + 102);

    auto tr_x_xxxyz_yzzzz = pbuffer.data(idx_dip_hh + 103);

    auto tr_x_xxxyz_zzzzz = pbuffer.data(idx_dip_hh + 104);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_xxxy_xxxxy,  \
                             tr_x_xxxy_xxxyy,  \
                             tr_x_xxxy_xxyyy,  \
                             tr_x_xxxy_xyyyy,  \
                             tr_x_xxxy_yyyyy,  \
                             tr_x_xxxyz_xxxxx, \
                             tr_x_xxxyz_xxxxy, \
                             tr_x_xxxyz_xxxxz, \
                             tr_x_xxxyz_xxxyy, \
                             tr_x_xxxyz_xxxyz, \
                             tr_x_xxxyz_xxxzz, \
                             tr_x_xxxyz_xxyyy, \
                             tr_x_xxxyz_xxyyz, \
                             tr_x_xxxyz_xxyzz, \
                             tr_x_xxxyz_xxzzz, \
                             tr_x_xxxyz_xyyyy, \
                             tr_x_xxxyz_xyyyz, \
                             tr_x_xxxyz_xyyzz, \
                             tr_x_xxxyz_xyzzz, \
                             tr_x_xxxyz_xzzzz, \
                             tr_x_xxxyz_yyyyy, \
                             tr_x_xxxyz_yyyyz, \
                             tr_x_xxxyz_yyyzz, \
                             tr_x_xxxyz_yyzzz, \
                             tr_x_xxxyz_yzzzz, \
                             tr_x_xxxyz_zzzzz, \
                             tr_x_xxxz_xxxxx,  \
                             tr_x_xxxz_xxxxz,  \
                             tr_x_xxxz_xxxyz,  \
                             tr_x_xxxz_xxxz,   \
                             tr_x_xxxz_xxxzz,  \
                             tr_x_xxxz_xxyyz,  \
                             tr_x_xxxz_xxyz,   \
                             tr_x_xxxz_xxyzz,  \
                             tr_x_xxxz_xxzz,   \
                             tr_x_xxxz_xxzzz,  \
                             tr_x_xxxz_xyyyz,  \
                             tr_x_xxxz_xyyz,   \
                             tr_x_xxxz_xyyzz,  \
                             tr_x_xxxz_xyzz,   \
                             tr_x_xxxz_xyzzz,  \
                             tr_x_xxxz_xzzz,   \
                             tr_x_xxxz_xzzzz,  \
                             tr_x_xxxz_yyyyz,  \
                             tr_x_xxxz_yyyz,   \
                             tr_x_xxxz_yyyzz,  \
                             tr_x_xxxz_yyzz,   \
                             tr_x_xxxz_yyzzz,  \
                             tr_x_xxxz_yzzz,   \
                             tr_x_xxxz_yzzzz,  \
                             tr_x_xxxz_zzzz,   \
                             tr_x_xxxz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyz_xxxxx[i] = tr_x_xxxz_xxxxx[i] * pa_y[i];

        tr_x_xxxyz_xxxxy[i] = tr_x_xxxy_xxxxy[i] * pa_z[i];

        tr_x_xxxyz_xxxxz[i] = tr_x_xxxz_xxxxz[i] * pa_y[i];

        tr_x_xxxyz_xxxyy[i] = tr_x_xxxy_xxxyy[i] * pa_z[i];

        tr_x_xxxyz_xxxyz[i] = tr_x_xxxz_xxxz[i] * fe_0 + tr_x_xxxz_xxxyz[i] * pa_y[i];

        tr_x_xxxyz_xxxzz[i] = tr_x_xxxz_xxxzz[i] * pa_y[i];

        tr_x_xxxyz_xxyyy[i] = tr_x_xxxy_xxyyy[i] * pa_z[i];

        tr_x_xxxyz_xxyyz[i] = 2.0 * tr_x_xxxz_xxyz[i] * fe_0 + tr_x_xxxz_xxyyz[i] * pa_y[i];

        tr_x_xxxyz_xxyzz[i] = tr_x_xxxz_xxzz[i] * fe_0 + tr_x_xxxz_xxyzz[i] * pa_y[i];

        tr_x_xxxyz_xxzzz[i] = tr_x_xxxz_xxzzz[i] * pa_y[i];

        tr_x_xxxyz_xyyyy[i] = tr_x_xxxy_xyyyy[i] * pa_z[i];

        tr_x_xxxyz_xyyyz[i] = 3.0 * tr_x_xxxz_xyyz[i] * fe_0 + tr_x_xxxz_xyyyz[i] * pa_y[i];

        tr_x_xxxyz_xyyzz[i] = 2.0 * tr_x_xxxz_xyzz[i] * fe_0 + tr_x_xxxz_xyyzz[i] * pa_y[i];

        tr_x_xxxyz_xyzzz[i] = tr_x_xxxz_xzzz[i] * fe_0 + tr_x_xxxz_xyzzz[i] * pa_y[i];

        tr_x_xxxyz_xzzzz[i] = tr_x_xxxz_xzzzz[i] * pa_y[i];

        tr_x_xxxyz_yyyyy[i] = tr_x_xxxy_yyyyy[i] * pa_z[i];

        tr_x_xxxyz_yyyyz[i] = 4.0 * tr_x_xxxz_yyyz[i] * fe_0 + tr_x_xxxz_yyyyz[i] * pa_y[i];

        tr_x_xxxyz_yyyzz[i] = 3.0 * tr_x_xxxz_yyzz[i] * fe_0 + tr_x_xxxz_yyyzz[i] * pa_y[i];

        tr_x_xxxyz_yyzzz[i] = 2.0 * tr_x_xxxz_yzzz[i] * fe_0 + tr_x_xxxz_yyzzz[i] * pa_y[i];

        tr_x_xxxyz_yzzzz[i] = tr_x_xxxz_zzzz[i] * fe_0 + tr_x_xxxz_yzzzz[i] * pa_y[i];

        tr_x_xxxyz_zzzzz[i] = tr_x_xxxz_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xxx_xxxxx,   \
                             tr_x_xxx_xxxxy,   \
                             tr_x_xxx_xxxxz,   \
                             tr_x_xxx_xxxyy,   \
                             tr_x_xxx_xxxyz,   \
                             tr_x_xxx_xxxzz,   \
                             tr_x_xxx_xxyyy,   \
                             tr_x_xxx_xxyyz,   \
                             tr_x_xxx_xxyzz,   \
                             tr_x_xxx_xxzzz,   \
                             tr_x_xxx_xyyyy,   \
                             tr_x_xxx_xyyyz,   \
                             tr_x_xxx_xyyzz,   \
                             tr_x_xxx_xyzzz,   \
                             tr_x_xxx_xzzzz,   \
                             tr_x_xxx_yyyyy,   \
                             tr_x_xxxz_xxxx,   \
                             tr_x_xxxz_xxxxx,  \
                             tr_x_xxxz_xxxxy,  \
                             tr_x_xxxz_xxxxz,  \
                             tr_x_xxxz_xxxy,   \
                             tr_x_xxxz_xxxyy,  \
                             tr_x_xxxz_xxxyz,  \
                             tr_x_xxxz_xxxz,   \
                             tr_x_xxxz_xxxzz,  \
                             tr_x_xxxz_xxyy,   \
                             tr_x_xxxz_xxyyy,  \
                             tr_x_xxxz_xxyyz,  \
                             tr_x_xxxz_xxyz,   \
                             tr_x_xxxz_xxyzz,  \
                             tr_x_xxxz_xxzz,   \
                             tr_x_xxxz_xxzzz,  \
                             tr_x_xxxz_xyyy,   \
                             tr_x_xxxz_xyyyy,  \
                             tr_x_xxxz_xyyyz,  \
                             tr_x_xxxz_xyyz,   \
                             tr_x_xxxz_xyyzz,  \
                             tr_x_xxxz_xyzz,   \
                             tr_x_xxxz_xyzzz,  \
                             tr_x_xxxz_xzzz,   \
                             tr_x_xxxz_xzzzz,  \
                             tr_x_xxxz_yyyyy,  \
                             tr_x_xxxzz_xxxxx, \
                             tr_x_xxxzz_xxxxy, \
                             tr_x_xxxzz_xxxxz, \
                             tr_x_xxxzz_xxxyy, \
                             tr_x_xxxzz_xxxyz, \
                             tr_x_xxxzz_xxxzz, \
                             tr_x_xxxzz_xxyyy, \
                             tr_x_xxxzz_xxyyz, \
                             tr_x_xxxzz_xxyzz, \
                             tr_x_xxxzz_xxzzz, \
                             tr_x_xxxzz_xyyyy, \
                             tr_x_xxxzz_xyyyz, \
                             tr_x_xxxzz_xyyzz, \
                             tr_x_xxxzz_xyzzz, \
                             tr_x_xxxzz_xzzzz, \
                             tr_x_xxxzz_yyyyy, \
                             tr_x_xxxzz_yyyyz, \
                             tr_x_xxxzz_yyyzz, \
                             tr_x_xxxzz_yyzzz, \
                             tr_x_xxxzz_yzzzz, \
                             tr_x_xxxzz_zzzzz, \
                             tr_x_xxzz_yyyyz,  \
                             tr_x_xxzz_yyyzz,  \
                             tr_x_xxzz_yyzzz,  \
                             tr_x_xxzz_yzzzz,  \
                             tr_x_xxzz_zzzzz,  \
                             tr_x_xzz_yyyyz,   \
                             tr_x_xzz_yyyzz,   \
                             tr_x_xzz_yyzzz,   \
                             tr_x_xzz_yzzzz,   \
                             tr_x_xzz_zzzzz,   \
                             ts_xxzz_yyyyz,    \
                             ts_xxzz_yyyzz,    \
                             ts_xxzz_yyzzz,    \
                             ts_xxzz_yzzzz,    \
                             ts_xxzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzz_xxxxx[i] = tr_x_xxx_xxxxx[i] * fe_0 + tr_x_xxxz_xxxxx[i] * pa_z[i];

        tr_x_xxxzz_xxxxy[i] = tr_x_xxx_xxxxy[i] * fe_0 + tr_x_xxxz_xxxxy[i] * pa_z[i];

        tr_x_xxxzz_xxxxz[i] = tr_x_xxx_xxxxz[i] * fe_0 + tr_x_xxxz_xxxx[i] * fe_0 + tr_x_xxxz_xxxxz[i] * pa_z[i];

        tr_x_xxxzz_xxxyy[i] = tr_x_xxx_xxxyy[i] * fe_0 + tr_x_xxxz_xxxyy[i] * pa_z[i];

        tr_x_xxxzz_xxxyz[i] = tr_x_xxx_xxxyz[i] * fe_0 + tr_x_xxxz_xxxy[i] * fe_0 + tr_x_xxxz_xxxyz[i] * pa_z[i];

        tr_x_xxxzz_xxxzz[i] = tr_x_xxx_xxxzz[i] * fe_0 + 2.0 * tr_x_xxxz_xxxz[i] * fe_0 + tr_x_xxxz_xxxzz[i] * pa_z[i];

        tr_x_xxxzz_xxyyy[i] = tr_x_xxx_xxyyy[i] * fe_0 + tr_x_xxxz_xxyyy[i] * pa_z[i];

        tr_x_xxxzz_xxyyz[i] = tr_x_xxx_xxyyz[i] * fe_0 + tr_x_xxxz_xxyy[i] * fe_0 + tr_x_xxxz_xxyyz[i] * pa_z[i];

        tr_x_xxxzz_xxyzz[i] = tr_x_xxx_xxyzz[i] * fe_0 + 2.0 * tr_x_xxxz_xxyz[i] * fe_0 + tr_x_xxxz_xxyzz[i] * pa_z[i];

        tr_x_xxxzz_xxzzz[i] = tr_x_xxx_xxzzz[i] * fe_0 + 3.0 * tr_x_xxxz_xxzz[i] * fe_0 + tr_x_xxxz_xxzzz[i] * pa_z[i];

        tr_x_xxxzz_xyyyy[i] = tr_x_xxx_xyyyy[i] * fe_0 + tr_x_xxxz_xyyyy[i] * pa_z[i];

        tr_x_xxxzz_xyyyz[i] = tr_x_xxx_xyyyz[i] * fe_0 + tr_x_xxxz_xyyy[i] * fe_0 + tr_x_xxxz_xyyyz[i] * pa_z[i];

        tr_x_xxxzz_xyyzz[i] = tr_x_xxx_xyyzz[i] * fe_0 + 2.0 * tr_x_xxxz_xyyz[i] * fe_0 + tr_x_xxxz_xyyzz[i] * pa_z[i];

        tr_x_xxxzz_xyzzz[i] = tr_x_xxx_xyzzz[i] * fe_0 + 3.0 * tr_x_xxxz_xyzz[i] * fe_0 + tr_x_xxxz_xyzzz[i] * pa_z[i];

        tr_x_xxxzz_xzzzz[i] = tr_x_xxx_xzzzz[i] * fe_0 + 4.0 * tr_x_xxxz_xzzz[i] * fe_0 + tr_x_xxxz_xzzzz[i] * pa_z[i];

        tr_x_xxxzz_yyyyy[i] = tr_x_xxx_yyyyy[i] * fe_0 + tr_x_xxxz_yyyyy[i] * pa_z[i];

        tr_x_xxxzz_yyyyz[i] = 2.0 * tr_x_xzz_yyyyz[i] * fe_0 + ts_xxzz_yyyyz[i] * fe_0 + tr_x_xxzz_yyyyz[i] * pa_x[i];

        tr_x_xxxzz_yyyzz[i] = 2.0 * tr_x_xzz_yyyzz[i] * fe_0 + ts_xxzz_yyyzz[i] * fe_0 + tr_x_xxzz_yyyzz[i] * pa_x[i];

        tr_x_xxxzz_yyzzz[i] = 2.0 * tr_x_xzz_yyzzz[i] * fe_0 + ts_xxzz_yyzzz[i] * fe_0 + tr_x_xxzz_yyzzz[i] * pa_x[i];

        tr_x_xxxzz_yzzzz[i] = 2.0 * tr_x_xzz_yzzzz[i] * fe_0 + ts_xxzz_yzzzz[i] * fe_0 + tr_x_xxzz_yzzzz[i] * pa_x[i];

        tr_x_xxxzz_zzzzz[i] = 2.0 * tr_x_xzz_zzzzz[i] * fe_0 + ts_xxzz_zzzzz[i] * fe_0 + tr_x_xxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xxy_xxxxx,   \
                             tr_x_xxy_xxxxy,   \
                             tr_x_xxy_xxxxz,   \
                             tr_x_xxy_xxxyy,   \
                             tr_x_xxy_xxxyz,   \
                             tr_x_xxy_xxxzz,   \
                             tr_x_xxy_xxyyy,   \
                             tr_x_xxy_xxyyz,   \
                             tr_x_xxy_xxyzz,   \
                             tr_x_xxy_xxzzz,   \
                             tr_x_xxy_xyyyy,   \
                             tr_x_xxy_xyyyz,   \
                             tr_x_xxy_xyyzz,   \
                             tr_x_xxy_xyzzz,   \
                             tr_x_xxy_xzzzz,   \
                             tr_x_xxy_zzzzz,   \
                             tr_x_xxyy_xxxx,   \
                             tr_x_xxyy_xxxxx,  \
                             tr_x_xxyy_xxxxy,  \
                             tr_x_xxyy_xxxxz,  \
                             tr_x_xxyy_xxxy,   \
                             tr_x_xxyy_xxxyy,  \
                             tr_x_xxyy_xxxyz,  \
                             tr_x_xxyy_xxxz,   \
                             tr_x_xxyy_xxxzz,  \
                             tr_x_xxyy_xxyy,   \
                             tr_x_xxyy_xxyyy,  \
                             tr_x_xxyy_xxyyz,  \
                             tr_x_xxyy_xxyz,   \
                             tr_x_xxyy_xxyzz,  \
                             tr_x_xxyy_xxzz,   \
                             tr_x_xxyy_xxzzz,  \
                             tr_x_xxyy_xyyy,   \
                             tr_x_xxyy_xyyyy,  \
                             tr_x_xxyy_xyyyz,  \
                             tr_x_xxyy_xyyz,   \
                             tr_x_xxyy_xyyzz,  \
                             tr_x_xxyy_xyzz,   \
                             tr_x_xxyy_xyzzz,  \
                             tr_x_xxyy_xzzz,   \
                             tr_x_xxyy_xzzzz,  \
                             tr_x_xxyy_zzzzz,  \
                             tr_x_xxyyy_xxxxx, \
                             tr_x_xxyyy_xxxxy, \
                             tr_x_xxyyy_xxxxz, \
                             tr_x_xxyyy_xxxyy, \
                             tr_x_xxyyy_xxxyz, \
                             tr_x_xxyyy_xxxzz, \
                             tr_x_xxyyy_xxyyy, \
                             tr_x_xxyyy_xxyyz, \
                             tr_x_xxyyy_xxyzz, \
                             tr_x_xxyyy_xxzzz, \
                             tr_x_xxyyy_xyyyy, \
                             tr_x_xxyyy_xyyyz, \
                             tr_x_xxyyy_xyyzz, \
                             tr_x_xxyyy_xyzzz, \
                             tr_x_xxyyy_xzzzz, \
                             tr_x_xxyyy_yyyyy, \
                             tr_x_xxyyy_yyyyz, \
                             tr_x_xxyyy_yyyzz, \
                             tr_x_xxyyy_yyzzz, \
                             tr_x_xxyyy_yzzzz, \
                             tr_x_xxyyy_zzzzz, \
                             tr_x_xyyy_yyyyy,  \
                             tr_x_xyyy_yyyyz,  \
                             tr_x_xyyy_yyyzz,  \
                             tr_x_xyyy_yyzzz,  \
                             tr_x_xyyy_yzzzz,  \
                             tr_x_yyy_yyyyy,   \
                             tr_x_yyy_yyyyz,   \
                             tr_x_yyy_yyyzz,   \
                             tr_x_yyy_yyzzz,   \
                             tr_x_yyy_yzzzz,   \
                             ts_xyyy_yyyyy,    \
                             ts_xyyy_yyyyz,    \
                             ts_xyyy_yyyzz,    \
                             ts_xyyy_yyzzz,    \
                             ts_xyyy_yzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyy_xxxxx[i] = 2.0 * tr_x_xxy_xxxxx[i] * fe_0 + tr_x_xxyy_xxxxx[i] * pa_y[i];

        tr_x_xxyyy_xxxxy[i] = 2.0 * tr_x_xxy_xxxxy[i] * fe_0 + tr_x_xxyy_xxxx[i] * fe_0 + tr_x_xxyy_xxxxy[i] * pa_y[i];

        tr_x_xxyyy_xxxxz[i] = 2.0 * tr_x_xxy_xxxxz[i] * fe_0 + tr_x_xxyy_xxxxz[i] * pa_y[i];

        tr_x_xxyyy_xxxyy[i] = 2.0 * tr_x_xxy_xxxyy[i] * fe_0 + 2.0 * tr_x_xxyy_xxxy[i] * fe_0 + tr_x_xxyy_xxxyy[i] * pa_y[i];

        tr_x_xxyyy_xxxyz[i] = 2.0 * tr_x_xxy_xxxyz[i] * fe_0 + tr_x_xxyy_xxxz[i] * fe_0 + tr_x_xxyy_xxxyz[i] * pa_y[i];

        tr_x_xxyyy_xxxzz[i] = 2.0 * tr_x_xxy_xxxzz[i] * fe_0 + tr_x_xxyy_xxxzz[i] * pa_y[i];

        tr_x_xxyyy_xxyyy[i] = 2.0 * tr_x_xxy_xxyyy[i] * fe_0 + 3.0 * tr_x_xxyy_xxyy[i] * fe_0 + tr_x_xxyy_xxyyy[i] * pa_y[i];

        tr_x_xxyyy_xxyyz[i] = 2.0 * tr_x_xxy_xxyyz[i] * fe_0 + 2.0 * tr_x_xxyy_xxyz[i] * fe_0 + tr_x_xxyy_xxyyz[i] * pa_y[i];

        tr_x_xxyyy_xxyzz[i] = 2.0 * tr_x_xxy_xxyzz[i] * fe_0 + tr_x_xxyy_xxzz[i] * fe_0 + tr_x_xxyy_xxyzz[i] * pa_y[i];

        tr_x_xxyyy_xxzzz[i] = 2.0 * tr_x_xxy_xxzzz[i] * fe_0 + tr_x_xxyy_xxzzz[i] * pa_y[i];

        tr_x_xxyyy_xyyyy[i] = 2.0 * tr_x_xxy_xyyyy[i] * fe_0 + 4.0 * tr_x_xxyy_xyyy[i] * fe_0 + tr_x_xxyy_xyyyy[i] * pa_y[i];

        tr_x_xxyyy_xyyyz[i] = 2.0 * tr_x_xxy_xyyyz[i] * fe_0 + 3.0 * tr_x_xxyy_xyyz[i] * fe_0 + tr_x_xxyy_xyyyz[i] * pa_y[i];

        tr_x_xxyyy_xyyzz[i] = 2.0 * tr_x_xxy_xyyzz[i] * fe_0 + 2.0 * tr_x_xxyy_xyzz[i] * fe_0 + tr_x_xxyy_xyyzz[i] * pa_y[i];

        tr_x_xxyyy_xyzzz[i] = 2.0 * tr_x_xxy_xyzzz[i] * fe_0 + tr_x_xxyy_xzzz[i] * fe_0 + tr_x_xxyy_xyzzz[i] * pa_y[i];

        tr_x_xxyyy_xzzzz[i] = 2.0 * tr_x_xxy_xzzzz[i] * fe_0 + tr_x_xxyy_xzzzz[i] * pa_y[i];

        tr_x_xxyyy_yyyyy[i] = tr_x_yyy_yyyyy[i] * fe_0 + ts_xyyy_yyyyy[i] * fe_0 + tr_x_xyyy_yyyyy[i] * pa_x[i];

        tr_x_xxyyy_yyyyz[i] = tr_x_yyy_yyyyz[i] * fe_0 + ts_xyyy_yyyyz[i] * fe_0 + tr_x_xyyy_yyyyz[i] * pa_x[i];

        tr_x_xxyyy_yyyzz[i] = tr_x_yyy_yyyzz[i] * fe_0 + ts_xyyy_yyyzz[i] * fe_0 + tr_x_xyyy_yyyzz[i] * pa_x[i];

        tr_x_xxyyy_yyzzz[i] = tr_x_yyy_yyzzz[i] * fe_0 + ts_xyyy_yyzzz[i] * fe_0 + tr_x_xyyy_yyzzz[i] * pa_x[i];

        tr_x_xxyyy_yzzzz[i] = tr_x_yyy_yzzzz[i] * fe_0 + ts_xyyy_yzzzz[i] * fe_0 + tr_x_xyyy_yzzzz[i] * pa_x[i];

        tr_x_xxyyy_zzzzz[i] = 2.0 * tr_x_xxy_zzzzz[i] * fe_0 + tr_x_xxyy_zzzzz[i] * pa_y[i];
    }

    // Set up 147-168 components of targeted buffer : HH

    auto tr_x_xxyyz_xxxxx = pbuffer.data(idx_dip_hh + 147);

    auto tr_x_xxyyz_xxxxy = pbuffer.data(idx_dip_hh + 148);

    auto tr_x_xxyyz_xxxxz = pbuffer.data(idx_dip_hh + 149);

    auto tr_x_xxyyz_xxxyy = pbuffer.data(idx_dip_hh + 150);

    auto tr_x_xxyyz_xxxyz = pbuffer.data(idx_dip_hh + 151);

    auto tr_x_xxyyz_xxxzz = pbuffer.data(idx_dip_hh + 152);

    auto tr_x_xxyyz_xxyyy = pbuffer.data(idx_dip_hh + 153);

    auto tr_x_xxyyz_xxyyz = pbuffer.data(idx_dip_hh + 154);

    auto tr_x_xxyyz_xxyzz = pbuffer.data(idx_dip_hh + 155);

    auto tr_x_xxyyz_xxzzz = pbuffer.data(idx_dip_hh + 156);

    auto tr_x_xxyyz_xyyyy = pbuffer.data(idx_dip_hh + 157);

    auto tr_x_xxyyz_xyyyz = pbuffer.data(idx_dip_hh + 158);

    auto tr_x_xxyyz_xyyzz = pbuffer.data(idx_dip_hh + 159);

    auto tr_x_xxyyz_xyzzz = pbuffer.data(idx_dip_hh + 160);

    auto tr_x_xxyyz_xzzzz = pbuffer.data(idx_dip_hh + 161);

    auto tr_x_xxyyz_yyyyy = pbuffer.data(idx_dip_hh + 162);

    auto tr_x_xxyyz_yyyyz = pbuffer.data(idx_dip_hh + 163);

    auto tr_x_xxyyz_yyyzz = pbuffer.data(idx_dip_hh + 164);

    auto tr_x_xxyyz_yyzzz = pbuffer.data(idx_dip_hh + 165);

    auto tr_x_xxyyz_yzzzz = pbuffer.data(idx_dip_hh + 166);

    auto tr_x_xxyyz_zzzzz = pbuffer.data(idx_dip_hh + 167);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_xxyy_xxxxx,  \
                             tr_x_xxyy_xxxxy,  \
                             tr_x_xxyy_xxxy,   \
                             tr_x_xxyy_xxxyy,  \
                             tr_x_xxyy_xxxyz,  \
                             tr_x_xxyy_xxyy,   \
                             tr_x_xxyy_xxyyy,  \
                             tr_x_xxyy_xxyyz,  \
                             tr_x_xxyy_xxyz,   \
                             tr_x_xxyy_xxyzz,  \
                             tr_x_xxyy_xyyy,   \
                             tr_x_xxyy_xyyyy,  \
                             tr_x_xxyy_xyyyz,  \
                             tr_x_xxyy_xyyz,   \
                             tr_x_xxyy_xyyzz,  \
                             tr_x_xxyy_xyzz,   \
                             tr_x_xxyy_xyzzz,  \
                             tr_x_xxyy_yyyy,   \
                             tr_x_xxyy_yyyyy,  \
                             tr_x_xxyy_yyyyz,  \
                             tr_x_xxyy_yyyz,   \
                             tr_x_xxyy_yyyzz,  \
                             tr_x_xxyy_yyzz,   \
                             tr_x_xxyy_yyzzz,  \
                             tr_x_xxyy_yzzz,   \
                             tr_x_xxyy_yzzzz,  \
                             tr_x_xxyyz_xxxxx, \
                             tr_x_xxyyz_xxxxy, \
                             tr_x_xxyyz_xxxxz, \
                             tr_x_xxyyz_xxxyy, \
                             tr_x_xxyyz_xxxyz, \
                             tr_x_xxyyz_xxxzz, \
                             tr_x_xxyyz_xxyyy, \
                             tr_x_xxyyz_xxyyz, \
                             tr_x_xxyyz_xxyzz, \
                             tr_x_xxyyz_xxzzz, \
                             tr_x_xxyyz_xyyyy, \
                             tr_x_xxyyz_xyyyz, \
                             tr_x_xxyyz_xyyzz, \
                             tr_x_xxyyz_xyzzz, \
                             tr_x_xxyyz_xzzzz, \
                             tr_x_xxyyz_yyyyy, \
                             tr_x_xxyyz_yyyyz, \
                             tr_x_xxyyz_yyyzz, \
                             tr_x_xxyyz_yyzzz, \
                             tr_x_xxyyz_yzzzz, \
                             tr_x_xxyyz_zzzzz, \
                             tr_x_xxyz_xxxxz,  \
                             tr_x_xxyz_xxxzz,  \
                             tr_x_xxyz_xxzzz,  \
                             tr_x_xxyz_xzzzz,  \
                             tr_x_xxyz_zzzzz,  \
                             tr_x_xxz_xxxxz,   \
                             tr_x_xxz_xxxzz,   \
                             tr_x_xxz_xxzzz,   \
                             tr_x_xxz_xzzzz,   \
                             tr_x_xxz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyz_xxxxx[i] = tr_x_xxyy_xxxxx[i] * pa_z[i];

        tr_x_xxyyz_xxxxy[i] = tr_x_xxyy_xxxxy[i] * pa_z[i];

        tr_x_xxyyz_xxxxz[i] = tr_x_xxz_xxxxz[i] * fe_0 + tr_x_xxyz_xxxxz[i] * pa_y[i];

        tr_x_xxyyz_xxxyy[i] = tr_x_xxyy_xxxyy[i] * pa_z[i];

        tr_x_xxyyz_xxxyz[i] = tr_x_xxyy_xxxy[i] * fe_0 + tr_x_xxyy_xxxyz[i] * pa_z[i];

        tr_x_xxyyz_xxxzz[i] = tr_x_xxz_xxxzz[i] * fe_0 + tr_x_xxyz_xxxzz[i] * pa_y[i];

        tr_x_xxyyz_xxyyy[i] = tr_x_xxyy_xxyyy[i] * pa_z[i];

        tr_x_xxyyz_xxyyz[i] = tr_x_xxyy_xxyy[i] * fe_0 + tr_x_xxyy_xxyyz[i] * pa_z[i];

        tr_x_xxyyz_xxyzz[i] = 2.0 * tr_x_xxyy_xxyz[i] * fe_0 + tr_x_xxyy_xxyzz[i] * pa_z[i];

        tr_x_xxyyz_xxzzz[i] = tr_x_xxz_xxzzz[i] * fe_0 + tr_x_xxyz_xxzzz[i] * pa_y[i];

        tr_x_xxyyz_xyyyy[i] = tr_x_xxyy_xyyyy[i] * pa_z[i];

        tr_x_xxyyz_xyyyz[i] = tr_x_xxyy_xyyy[i] * fe_0 + tr_x_xxyy_xyyyz[i] * pa_z[i];

        tr_x_xxyyz_xyyzz[i] = 2.0 * tr_x_xxyy_xyyz[i] * fe_0 + tr_x_xxyy_xyyzz[i] * pa_z[i];

        tr_x_xxyyz_xyzzz[i] = 3.0 * tr_x_xxyy_xyzz[i] * fe_0 + tr_x_xxyy_xyzzz[i] * pa_z[i];

        tr_x_xxyyz_xzzzz[i] = tr_x_xxz_xzzzz[i] * fe_0 + tr_x_xxyz_xzzzz[i] * pa_y[i];

        tr_x_xxyyz_yyyyy[i] = tr_x_xxyy_yyyyy[i] * pa_z[i];

        tr_x_xxyyz_yyyyz[i] = tr_x_xxyy_yyyy[i] * fe_0 + tr_x_xxyy_yyyyz[i] * pa_z[i];

        tr_x_xxyyz_yyyzz[i] = 2.0 * tr_x_xxyy_yyyz[i] * fe_0 + tr_x_xxyy_yyyzz[i] * pa_z[i];

        tr_x_xxyyz_yyzzz[i] = 3.0 * tr_x_xxyy_yyzz[i] * fe_0 + tr_x_xxyy_yyzzz[i] * pa_z[i];

        tr_x_xxyyz_yzzzz[i] = 4.0 * tr_x_xxyy_yzzz[i] * fe_0 + tr_x_xxyy_yzzzz[i] * pa_z[i];

        tr_x_xxyyz_zzzzz[i] = tr_x_xxz_zzzzz[i] * fe_0 + tr_x_xxyz_zzzzz[i] * pa_y[i];
    }

    // Set up 168-189 components of targeted buffer : HH

    auto tr_x_xxyzz_xxxxx = pbuffer.data(idx_dip_hh + 168);

    auto tr_x_xxyzz_xxxxy = pbuffer.data(idx_dip_hh + 169);

    auto tr_x_xxyzz_xxxxz = pbuffer.data(idx_dip_hh + 170);

    auto tr_x_xxyzz_xxxyy = pbuffer.data(idx_dip_hh + 171);

    auto tr_x_xxyzz_xxxyz = pbuffer.data(idx_dip_hh + 172);

    auto tr_x_xxyzz_xxxzz = pbuffer.data(idx_dip_hh + 173);

    auto tr_x_xxyzz_xxyyy = pbuffer.data(idx_dip_hh + 174);

    auto tr_x_xxyzz_xxyyz = pbuffer.data(idx_dip_hh + 175);

    auto tr_x_xxyzz_xxyzz = pbuffer.data(idx_dip_hh + 176);

    auto tr_x_xxyzz_xxzzz = pbuffer.data(idx_dip_hh + 177);

    auto tr_x_xxyzz_xyyyy = pbuffer.data(idx_dip_hh + 178);

    auto tr_x_xxyzz_xyyyz = pbuffer.data(idx_dip_hh + 179);

    auto tr_x_xxyzz_xyyzz = pbuffer.data(idx_dip_hh + 180);

    auto tr_x_xxyzz_xyzzz = pbuffer.data(idx_dip_hh + 181);

    auto tr_x_xxyzz_xzzzz = pbuffer.data(idx_dip_hh + 182);

    auto tr_x_xxyzz_yyyyy = pbuffer.data(idx_dip_hh + 183);

    auto tr_x_xxyzz_yyyyz = pbuffer.data(idx_dip_hh + 184);

    auto tr_x_xxyzz_yyyzz = pbuffer.data(idx_dip_hh + 185);

    auto tr_x_xxyzz_yyzzz = pbuffer.data(idx_dip_hh + 186);

    auto tr_x_xxyzz_yzzzz = pbuffer.data(idx_dip_hh + 187);

    auto tr_x_xxyzz_zzzzz = pbuffer.data(idx_dip_hh + 188);

#pragma omp simd aligned(pa_y,                 \
                             tr_x_xxyzz_xxxxx, \
                             tr_x_xxyzz_xxxxy, \
                             tr_x_xxyzz_xxxxz, \
                             tr_x_xxyzz_xxxyy, \
                             tr_x_xxyzz_xxxyz, \
                             tr_x_xxyzz_xxxzz, \
                             tr_x_xxyzz_xxyyy, \
                             tr_x_xxyzz_xxyyz, \
                             tr_x_xxyzz_xxyzz, \
                             tr_x_xxyzz_xxzzz, \
                             tr_x_xxyzz_xyyyy, \
                             tr_x_xxyzz_xyyyz, \
                             tr_x_xxyzz_xyyzz, \
                             tr_x_xxyzz_xyzzz, \
                             tr_x_xxyzz_xzzzz, \
                             tr_x_xxyzz_yyyyy, \
                             tr_x_xxyzz_yyyyz, \
                             tr_x_xxyzz_yyyzz, \
                             tr_x_xxyzz_yyzzz, \
                             tr_x_xxyzz_yzzzz, \
                             tr_x_xxyzz_zzzzz, \
                             tr_x_xxzz_xxxx,   \
                             tr_x_xxzz_xxxxx,  \
                             tr_x_xxzz_xxxxy,  \
                             tr_x_xxzz_xxxxz,  \
                             tr_x_xxzz_xxxy,   \
                             tr_x_xxzz_xxxyy,  \
                             tr_x_xxzz_xxxyz,  \
                             tr_x_xxzz_xxxz,   \
                             tr_x_xxzz_xxxzz,  \
                             tr_x_xxzz_xxyy,   \
                             tr_x_xxzz_xxyyy,  \
                             tr_x_xxzz_xxyyz,  \
                             tr_x_xxzz_xxyz,   \
                             tr_x_xxzz_xxyzz,  \
                             tr_x_xxzz_xxzz,   \
                             tr_x_xxzz_xxzzz,  \
                             tr_x_xxzz_xyyy,   \
                             tr_x_xxzz_xyyyy,  \
                             tr_x_xxzz_xyyyz,  \
                             tr_x_xxzz_xyyz,   \
                             tr_x_xxzz_xyyzz,  \
                             tr_x_xxzz_xyzz,   \
                             tr_x_xxzz_xyzzz,  \
                             tr_x_xxzz_xzzz,   \
                             tr_x_xxzz_xzzzz,  \
                             tr_x_xxzz_yyyy,   \
                             tr_x_xxzz_yyyyy,  \
                             tr_x_xxzz_yyyyz,  \
                             tr_x_xxzz_yyyz,   \
                             tr_x_xxzz_yyyzz,  \
                             tr_x_xxzz_yyzz,   \
                             tr_x_xxzz_yyzzz,  \
                             tr_x_xxzz_yzzz,   \
                             tr_x_xxzz_yzzzz,  \
                             tr_x_xxzz_zzzz,   \
                             tr_x_xxzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzz_xxxxx[i] = tr_x_xxzz_xxxxx[i] * pa_y[i];

        tr_x_xxyzz_xxxxy[i] = tr_x_xxzz_xxxx[i] * fe_0 + tr_x_xxzz_xxxxy[i] * pa_y[i];

        tr_x_xxyzz_xxxxz[i] = tr_x_xxzz_xxxxz[i] * pa_y[i];

        tr_x_xxyzz_xxxyy[i] = 2.0 * tr_x_xxzz_xxxy[i] * fe_0 + tr_x_xxzz_xxxyy[i] * pa_y[i];

        tr_x_xxyzz_xxxyz[i] = tr_x_xxzz_xxxz[i] * fe_0 + tr_x_xxzz_xxxyz[i] * pa_y[i];

        tr_x_xxyzz_xxxzz[i] = tr_x_xxzz_xxxzz[i] * pa_y[i];

        tr_x_xxyzz_xxyyy[i] = 3.0 * tr_x_xxzz_xxyy[i] * fe_0 + tr_x_xxzz_xxyyy[i] * pa_y[i];

        tr_x_xxyzz_xxyyz[i] = 2.0 * tr_x_xxzz_xxyz[i] * fe_0 + tr_x_xxzz_xxyyz[i] * pa_y[i];

        tr_x_xxyzz_xxyzz[i] = tr_x_xxzz_xxzz[i] * fe_0 + tr_x_xxzz_xxyzz[i] * pa_y[i];

        tr_x_xxyzz_xxzzz[i] = tr_x_xxzz_xxzzz[i] * pa_y[i];

        tr_x_xxyzz_xyyyy[i] = 4.0 * tr_x_xxzz_xyyy[i] * fe_0 + tr_x_xxzz_xyyyy[i] * pa_y[i];

        tr_x_xxyzz_xyyyz[i] = 3.0 * tr_x_xxzz_xyyz[i] * fe_0 + tr_x_xxzz_xyyyz[i] * pa_y[i];

        tr_x_xxyzz_xyyzz[i] = 2.0 * tr_x_xxzz_xyzz[i] * fe_0 + tr_x_xxzz_xyyzz[i] * pa_y[i];

        tr_x_xxyzz_xyzzz[i] = tr_x_xxzz_xzzz[i] * fe_0 + tr_x_xxzz_xyzzz[i] * pa_y[i];

        tr_x_xxyzz_xzzzz[i] = tr_x_xxzz_xzzzz[i] * pa_y[i];

        tr_x_xxyzz_yyyyy[i] = 5.0 * tr_x_xxzz_yyyy[i] * fe_0 + tr_x_xxzz_yyyyy[i] * pa_y[i];

        tr_x_xxyzz_yyyyz[i] = 4.0 * tr_x_xxzz_yyyz[i] * fe_0 + tr_x_xxzz_yyyyz[i] * pa_y[i];

        tr_x_xxyzz_yyyzz[i] = 3.0 * tr_x_xxzz_yyzz[i] * fe_0 + tr_x_xxzz_yyyzz[i] * pa_y[i];

        tr_x_xxyzz_yyzzz[i] = 2.0 * tr_x_xxzz_yzzz[i] * fe_0 + tr_x_xxzz_yyzzz[i] * pa_y[i];

        tr_x_xxyzz_yzzzz[i] = tr_x_xxzz_zzzz[i] * fe_0 + tr_x_xxzz_yzzzz[i] * pa_y[i];

        tr_x_xxyzz_zzzzz[i] = tr_x_xxzz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xxz_xxxxx,   \
                             tr_x_xxz_xxxxy,   \
                             tr_x_xxz_xxxxz,   \
                             tr_x_xxz_xxxyy,   \
                             tr_x_xxz_xxxyz,   \
                             tr_x_xxz_xxxzz,   \
                             tr_x_xxz_xxyyy,   \
                             tr_x_xxz_xxyyz,   \
                             tr_x_xxz_xxyzz,   \
                             tr_x_xxz_xxzzz,   \
                             tr_x_xxz_xyyyy,   \
                             tr_x_xxz_xyyyz,   \
                             tr_x_xxz_xyyzz,   \
                             tr_x_xxz_xyzzz,   \
                             tr_x_xxz_xzzzz,   \
                             tr_x_xxz_yyyyy,   \
                             tr_x_xxzz_xxxx,   \
                             tr_x_xxzz_xxxxx,  \
                             tr_x_xxzz_xxxxy,  \
                             tr_x_xxzz_xxxxz,  \
                             tr_x_xxzz_xxxy,   \
                             tr_x_xxzz_xxxyy,  \
                             tr_x_xxzz_xxxyz,  \
                             tr_x_xxzz_xxxz,   \
                             tr_x_xxzz_xxxzz,  \
                             tr_x_xxzz_xxyy,   \
                             tr_x_xxzz_xxyyy,  \
                             tr_x_xxzz_xxyyz,  \
                             tr_x_xxzz_xxyz,   \
                             tr_x_xxzz_xxyzz,  \
                             tr_x_xxzz_xxzz,   \
                             tr_x_xxzz_xxzzz,  \
                             tr_x_xxzz_xyyy,   \
                             tr_x_xxzz_xyyyy,  \
                             tr_x_xxzz_xyyyz,  \
                             tr_x_xxzz_xyyz,   \
                             tr_x_xxzz_xyyzz,  \
                             tr_x_xxzz_xyzz,   \
                             tr_x_xxzz_xyzzz,  \
                             tr_x_xxzz_xzzz,   \
                             tr_x_xxzz_xzzzz,  \
                             tr_x_xxzz_yyyyy,  \
                             tr_x_xxzzz_xxxxx, \
                             tr_x_xxzzz_xxxxy, \
                             tr_x_xxzzz_xxxxz, \
                             tr_x_xxzzz_xxxyy, \
                             tr_x_xxzzz_xxxyz, \
                             tr_x_xxzzz_xxxzz, \
                             tr_x_xxzzz_xxyyy, \
                             tr_x_xxzzz_xxyyz, \
                             tr_x_xxzzz_xxyzz, \
                             tr_x_xxzzz_xxzzz, \
                             tr_x_xxzzz_xyyyy, \
                             tr_x_xxzzz_xyyyz, \
                             tr_x_xxzzz_xyyzz, \
                             tr_x_xxzzz_xyzzz, \
                             tr_x_xxzzz_xzzzz, \
                             tr_x_xxzzz_yyyyy, \
                             tr_x_xxzzz_yyyyz, \
                             tr_x_xxzzz_yyyzz, \
                             tr_x_xxzzz_yyzzz, \
                             tr_x_xxzzz_yzzzz, \
                             tr_x_xxzzz_zzzzz, \
                             tr_x_xzzz_yyyyz,  \
                             tr_x_xzzz_yyyzz,  \
                             tr_x_xzzz_yyzzz,  \
                             tr_x_xzzz_yzzzz,  \
                             tr_x_xzzz_zzzzz,  \
                             tr_x_zzz_yyyyz,   \
                             tr_x_zzz_yyyzz,   \
                             tr_x_zzz_yyzzz,   \
                             tr_x_zzz_yzzzz,   \
                             tr_x_zzz_zzzzz,   \
                             ts_xzzz_yyyyz,    \
                             ts_xzzz_yyyzz,    \
                             ts_xzzz_yyzzz,    \
                             ts_xzzz_yzzzz,    \
                             ts_xzzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzz_xxxxx[i] = 2.0 * tr_x_xxz_xxxxx[i] * fe_0 + tr_x_xxzz_xxxxx[i] * pa_z[i];

        tr_x_xxzzz_xxxxy[i] = 2.0 * tr_x_xxz_xxxxy[i] * fe_0 + tr_x_xxzz_xxxxy[i] * pa_z[i];

        tr_x_xxzzz_xxxxz[i] = 2.0 * tr_x_xxz_xxxxz[i] * fe_0 + tr_x_xxzz_xxxx[i] * fe_0 + tr_x_xxzz_xxxxz[i] * pa_z[i];

        tr_x_xxzzz_xxxyy[i] = 2.0 * tr_x_xxz_xxxyy[i] * fe_0 + tr_x_xxzz_xxxyy[i] * pa_z[i];

        tr_x_xxzzz_xxxyz[i] = 2.0 * tr_x_xxz_xxxyz[i] * fe_0 + tr_x_xxzz_xxxy[i] * fe_0 + tr_x_xxzz_xxxyz[i] * pa_z[i];

        tr_x_xxzzz_xxxzz[i] = 2.0 * tr_x_xxz_xxxzz[i] * fe_0 + 2.0 * tr_x_xxzz_xxxz[i] * fe_0 + tr_x_xxzz_xxxzz[i] * pa_z[i];

        tr_x_xxzzz_xxyyy[i] = 2.0 * tr_x_xxz_xxyyy[i] * fe_0 + tr_x_xxzz_xxyyy[i] * pa_z[i];

        tr_x_xxzzz_xxyyz[i] = 2.0 * tr_x_xxz_xxyyz[i] * fe_0 + tr_x_xxzz_xxyy[i] * fe_0 + tr_x_xxzz_xxyyz[i] * pa_z[i];

        tr_x_xxzzz_xxyzz[i] = 2.0 * tr_x_xxz_xxyzz[i] * fe_0 + 2.0 * tr_x_xxzz_xxyz[i] * fe_0 + tr_x_xxzz_xxyzz[i] * pa_z[i];

        tr_x_xxzzz_xxzzz[i] = 2.0 * tr_x_xxz_xxzzz[i] * fe_0 + 3.0 * tr_x_xxzz_xxzz[i] * fe_0 + tr_x_xxzz_xxzzz[i] * pa_z[i];

        tr_x_xxzzz_xyyyy[i] = 2.0 * tr_x_xxz_xyyyy[i] * fe_0 + tr_x_xxzz_xyyyy[i] * pa_z[i];

        tr_x_xxzzz_xyyyz[i] = 2.0 * tr_x_xxz_xyyyz[i] * fe_0 + tr_x_xxzz_xyyy[i] * fe_0 + tr_x_xxzz_xyyyz[i] * pa_z[i];

        tr_x_xxzzz_xyyzz[i] = 2.0 * tr_x_xxz_xyyzz[i] * fe_0 + 2.0 * tr_x_xxzz_xyyz[i] * fe_0 + tr_x_xxzz_xyyzz[i] * pa_z[i];

        tr_x_xxzzz_xyzzz[i] = 2.0 * tr_x_xxz_xyzzz[i] * fe_0 + 3.0 * tr_x_xxzz_xyzz[i] * fe_0 + tr_x_xxzz_xyzzz[i] * pa_z[i];

        tr_x_xxzzz_xzzzz[i] = 2.0 * tr_x_xxz_xzzzz[i] * fe_0 + 4.0 * tr_x_xxzz_xzzz[i] * fe_0 + tr_x_xxzz_xzzzz[i] * pa_z[i];

        tr_x_xxzzz_yyyyy[i] = 2.0 * tr_x_xxz_yyyyy[i] * fe_0 + tr_x_xxzz_yyyyy[i] * pa_z[i];

        tr_x_xxzzz_yyyyz[i] = tr_x_zzz_yyyyz[i] * fe_0 + ts_xzzz_yyyyz[i] * fe_0 + tr_x_xzzz_yyyyz[i] * pa_x[i];

        tr_x_xxzzz_yyyzz[i] = tr_x_zzz_yyyzz[i] * fe_0 + ts_xzzz_yyyzz[i] * fe_0 + tr_x_xzzz_yyyzz[i] * pa_x[i];

        tr_x_xxzzz_yyzzz[i] = tr_x_zzz_yyzzz[i] * fe_0 + ts_xzzz_yyzzz[i] * fe_0 + tr_x_xzzz_yyzzz[i] * pa_x[i];

        tr_x_xxzzz_yzzzz[i] = tr_x_zzz_yzzzz[i] * fe_0 + ts_xzzz_yzzzz[i] * fe_0 + tr_x_xzzz_yzzzz[i] * pa_x[i];

        tr_x_xxzzz_zzzzz[i] = tr_x_zzz_zzzzz[i] * fe_0 + ts_xzzz_zzzzz[i] * fe_0 + tr_x_xzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : HH

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

    auto tr_x_xyyyy_zzzzz = pbuffer.data(idx_dip_hh + 230);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xyy_xxxxx,   \
                             tr_x_xyy_xxxxz,   \
                             tr_x_xyy_xxxzz,   \
                             tr_x_xyy_xxzzz,   \
                             tr_x_xyy_xzzzz,   \
                             tr_x_xyyy_xxxxx,  \
                             tr_x_xyyy_xxxxz,  \
                             tr_x_xyyy_xxxzz,  \
                             tr_x_xyyy_xxzzz,  \
                             tr_x_xyyy_xzzzz,  \
                             tr_x_xyyyy_xxxxx, \
                             tr_x_xyyyy_xxxxy, \
                             tr_x_xyyyy_xxxxz, \
                             tr_x_xyyyy_xxxyy, \
                             tr_x_xyyyy_xxxyz, \
                             tr_x_xyyyy_xxxzz, \
                             tr_x_xyyyy_xxyyy, \
                             tr_x_xyyyy_xxyyz, \
                             tr_x_xyyyy_xxyzz, \
                             tr_x_xyyyy_xxzzz, \
                             tr_x_xyyyy_xyyyy, \
                             tr_x_xyyyy_xyyyz, \
                             tr_x_xyyyy_xyyzz, \
                             tr_x_xyyyy_xyzzz, \
                             tr_x_xyyyy_xzzzz, \
                             tr_x_xyyyy_yyyyy, \
                             tr_x_xyyyy_yyyyz, \
                             tr_x_xyyyy_yyyzz, \
                             tr_x_xyyyy_yyzzz, \
                             tr_x_xyyyy_yzzzz, \
                             tr_x_xyyyy_zzzzz, \
                             tr_x_yyyy_xxxxy,  \
                             tr_x_yyyy_xxxy,   \
                             tr_x_yyyy_xxxyy,  \
                             tr_x_yyyy_xxxyz,  \
                             tr_x_yyyy_xxyy,   \
                             tr_x_yyyy_xxyyy,  \
                             tr_x_yyyy_xxyyz,  \
                             tr_x_yyyy_xxyz,   \
                             tr_x_yyyy_xxyzz,  \
                             tr_x_yyyy_xyyy,   \
                             tr_x_yyyy_xyyyy,  \
                             tr_x_yyyy_xyyyz,  \
                             tr_x_yyyy_xyyz,   \
                             tr_x_yyyy_xyyzz,  \
                             tr_x_yyyy_xyzz,   \
                             tr_x_yyyy_xyzzz,  \
                             tr_x_yyyy_yyyy,   \
                             tr_x_yyyy_yyyyy,  \
                             tr_x_yyyy_yyyyz,  \
                             tr_x_yyyy_yyyz,   \
                             tr_x_yyyy_yyyzz,  \
                             tr_x_yyyy_yyzz,   \
                             tr_x_yyyy_yyzzz,  \
                             tr_x_yyyy_yzzz,   \
                             tr_x_yyyy_yzzzz,  \
                             tr_x_yyyy_zzzzz,  \
                             ts_yyyy_xxxxy,    \
                             ts_yyyy_xxxyy,    \
                             ts_yyyy_xxxyz,    \
                             ts_yyyy_xxyyy,    \
                             ts_yyyy_xxyyz,    \
                             ts_yyyy_xxyzz,    \
                             ts_yyyy_xyyyy,    \
                             ts_yyyy_xyyyz,    \
                             ts_yyyy_xyyzz,    \
                             ts_yyyy_xyzzz,    \
                             ts_yyyy_yyyyy,    \
                             ts_yyyy_yyyyz,    \
                             ts_yyyy_yyyzz,    \
                             ts_yyyy_yyzzz,    \
                             ts_yyyy_yzzzz,    \
                             ts_yyyy_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyy_xxxxx[i] = 3.0 * tr_x_xyy_xxxxx[i] * fe_0 + tr_x_xyyy_xxxxx[i] * pa_y[i];

        tr_x_xyyyy_xxxxy[i] = 4.0 * tr_x_yyyy_xxxy[i] * fe_0 + ts_yyyy_xxxxy[i] * fe_0 + tr_x_yyyy_xxxxy[i] * pa_x[i];

        tr_x_xyyyy_xxxxz[i] = 3.0 * tr_x_xyy_xxxxz[i] * fe_0 + tr_x_xyyy_xxxxz[i] * pa_y[i];

        tr_x_xyyyy_xxxyy[i] = 3.0 * tr_x_yyyy_xxyy[i] * fe_0 + ts_yyyy_xxxyy[i] * fe_0 + tr_x_yyyy_xxxyy[i] * pa_x[i];

        tr_x_xyyyy_xxxyz[i] = 3.0 * tr_x_yyyy_xxyz[i] * fe_0 + ts_yyyy_xxxyz[i] * fe_0 + tr_x_yyyy_xxxyz[i] * pa_x[i];

        tr_x_xyyyy_xxxzz[i] = 3.0 * tr_x_xyy_xxxzz[i] * fe_0 + tr_x_xyyy_xxxzz[i] * pa_y[i];

        tr_x_xyyyy_xxyyy[i] = 2.0 * tr_x_yyyy_xyyy[i] * fe_0 + ts_yyyy_xxyyy[i] * fe_0 + tr_x_yyyy_xxyyy[i] * pa_x[i];

        tr_x_xyyyy_xxyyz[i] = 2.0 * tr_x_yyyy_xyyz[i] * fe_0 + ts_yyyy_xxyyz[i] * fe_0 + tr_x_yyyy_xxyyz[i] * pa_x[i];

        tr_x_xyyyy_xxyzz[i] = 2.0 * tr_x_yyyy_xyzz[i] * fe_0 + ts_yyyy_xxyzz[i] * fe_0 + tr_x_yyyy_xxyzz[i] * pa_x[i];

        tr_x_xyyyy_xxzzz[i] = 3.0 * tr_x_xyy_xxzzz[i] * fe_0 + tr_x_xyyy_xxzzz[i] * pa_y[i];

        tr_x_xyyyy_xyyyy[i] = tr_x_yyyy_yyyy[i] * fe_0 + ts_yyyy_xyyyy[i] * fe_0 + tr_x_yyyy_xyyyy[i] * pa_x[i];

        tr_x_xyyyy_xyyyz[i] = tr_x_yyyy_yyyz[i] * fe_0 + ts_yyyy_xyyyz[i] * fe_0 + tr_x_yyyy_xyyyz[i] * pa_x[i];

        tr_x_xyyyy_xyyzz[i] = tr_x_yyyy_yyzz[i] * fe_0 + ts_yyyy_xyyzz[i] * fe_0 + tr_x_yyyy_xyyzz[i] * pa_x[i];

        tr_x_xyyyy_xyzzz[i] = tr_x_yyyy_yzzz[i] * fe_0 + ts_yyyy_xyzzz[i] * fe_0 + tr_x_yyyy_xyzzz[i] * pa_x[i];

        tr_x_xyyyy_xzzzz[i] = 3.0 * tr_x_xyy_xzzzz[i] * fe_0 + tr_x_xyyy_xzzzz[i] * pa_y[i];

        tr_x_xyyyy_yyyyy[i] = ts_yyyy_yyyyy[i] * fe_0 + tr_x_yyyy_yyyyy[i] * pa_x[i];

        tr_x_xyyyy_yyyyz[i] = ts_yyyy_yyyyz[i] * fe_0 + tr_x_yyyy_yyyyz[i] * pa_x[i];

        tr_x_xyyyy_yyyzz[i] = ts_yyyy_yyyzz[i] * fe_0 + tr_x_yyyy_yyyzz[i] * pa_x[i];

        tr_x_xyyyy_yyzzz[i] = ts_yyyy_yyzzz[i] * fe_0 + tr_x_yyyy_yyzzz[i] * pa_x[i];

        tr_x_xyyyy_yzzzz[i] = ts_yyyy_yzzzz[i] * fe_0 + tr_x_yyyy_yzzzz[i] * pa_x[i];

        tr_x_xyyyy_zzzzz[i] = ts_yyyy_zzzzz[i] * fe_0 + tr_x_yyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : HH

    auto tr_x_xyyyz_xxxxx = pbuffer.data(idx_dip_hh + 231);

    auto tr_x_xyyyz_xxxxy = pbuffer.data(idx_dip_hh + 232);

    auto tr_x_xyyyz_xxxxz = pbuffer.data(idx_dip_hh + 233);

    auto tr_x_xyyyz_xxxyy = pbuffer.data(idx_dip_hh + 234);

    auto tr_x_xyyyz_xxxyz = pbuffer.data(idx_dip_hh + 235);

    auto tr_x_xyyyz_xxxzz = pbuffer.data(idx_dip_hh + 236);

    auto tr_x_xyyyz_xxyyy = pbuffer.data(idx_dip_hh + 237);

    auto tr_x_xyyyz_xxyyz = pbuffer.data(idx_dip_hh + 238);

    auto tr_x_xyyyz_xxyzz = pbuffer.data(idx_dip_hh + 239);

    auto tr_x_xyyyz_xxzzz = pbuffer.data(idx_dip_hh + 240);

    auto tr_x_xyyyz_xyyyy = pbuffer.data(idx_dip_hh + 241);

    auto tr_x_xyyyz_xyyyz = pbuffer.data(idx_dip_hh + 242);

    auto tr_x_xyyyz_xyyzz = pbuffer.data(idx_dip_hh + 243);

    auto tr_x_xyyyz_xyzzz = pbuffer.data(idx_dip_hh + 244);

    auto tr_x_xyyyz_xzzzz = pbuffer.data(idx_dip_hh + 245);

    auto tr_x_xyyyz_yyyyy = pbuffer.data(idx_dip_hh + 246);

    auto tr_x_xyyyz_yyyyz = pbuffer.data(idx_dip_hh + 247);

    auto tr_x_xyyyz_yyyzz = pbuffer.data(idx_dip_hh + 248);

    auto tr_x_xyyyz_yyzzz = pbuffer.data(idx_dip_hh + 249);

    auto tr_x_xyyyz_yzzzz = pbuffer.data(idx_dip_hh + 250);

    auto tr_x_xyyyz_zzzzz = pbuffer.data(idx_dip_hh + 251);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_x_xyyy_xxxxx,  \
                             tr_x_xyyy_xxxxy,  \
                             tr_x_xyyy_xxxy,   \
                             tr_x_xyyy_xxxyy,  \
                             tr_x_xyyy_xxxyz,  \
                             tr_x_xyyy_xxyy,   \
                             tr_x_xyyy_xxyyy,  \
                             tr_x_xyyy_xxyyz,  \
                             tr_x_xyyy_xxyz,   \
                             tr_x_xyyy_xxyzz,  \
                             tr_x_xyyy_xyyy,   \
                             tr_x_xyyy_xyyyy,  \
                             tr_x_xyyy_xyyyz,  \
                             tr_x_xyyy_xyyz,   \
                             tr_x_xyyy_xyyzz,  \
                             tr_x_xyyy_xyzz,   \
                             tr_x_xyyy_xyzzz,  \
                             tr_x_xyyy_yyyyy,  \
                             tr_x_xyyyz_xxxxx, \
                             tr_x_xyyyz_xxxxy, \
                             tr_x_xyyyz_xxxxz, \
                             tr_x_xyyyz_xxxyy, \
                             tr_x_xyyyz_xxxyz, \
                             tr_x_xyyyz_xxxzz, \
                             tr_x_xyyyz_xxyyy, \
                             tr_x_xyyyz_xxyyz, \
                             tr_x_xyyyz_xxyzz, \
                             tr_x_xyyyz_xxzzz, \
                             tr_x_xyyyz_xyyyy, \
                             tr_x_xyyyz_xyyyz, \
                             tr_x_xyyyz_xyyzz, \
                             tr_x_xyyyz_xyzzz, \
                             tr_x_xyyyz_xzzzz, \
                             tr_x_xyyyz_yyyyy, \
                             tr_x_xyyyz_yyyyz, \
                             tr_x_xyyyz_yyyzz, \
                             tr_x_xyyyz_yyzzz, \
                             tr_x_xyyyz_yzzzz, \
                             tr_x_xyyyz_zzzzz, \
                             tr_x_xyyz_xxxxz,  \
                             tr_x_xyyz_xxxzz,  \
                             tr_x_xyyz_xxzzz,  \
                             tr_x_xyyz_xzzzz,  \
                             tr_x_xyz_xxxxz,   \
                             tr_x_xyz_xxxzz,   \
                             tr_x_xyz_xxzzz,   \
                             tr_x_xyz_xzzzz,   \
                             tr_x_yyyz_yyyyz,  \
                             tr_x_yyyz_yyyzz,  \
                             tr_x_yyyz_yyzzz,  \
                             tr_x_yyyz_yzzzz,  \
                             tr_x_yyyz_zzzzz,  \
                             ts_yyyz_yyyyz,    \
                             ts_yyyz_yyyzz,    \
                             ts_yyyz_yyzzz,    \
                             ts_yyyz_yzzzz,    \
                             ts_yyyz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyz_xxxxx[i] = tr_x_xyyy_xxxxx[i] * pa_z[i];

        tr_x_xyyyz_xxxxy[i] = tr_x_xyyy_xxxxy[i] * pa_z[i];

        tr_x_xyyyz_xxxxz[i] = 2.0 * tr_x_xyz_xxxxz[i] * fe_0 + tr_x_xyyz_xxxxz[i] * pa_y[i];

        tr_x_xyyyz_xxxyy[i] = tr_x_xyyy_xxxyy[i] * pa_z[i];

        tr_x_xyyyz_xxxyz[i] = tr_x_xyyy_xxxy[i] * fe_0 + tr_x_xyyy_xxxyz[i] * pa_z[i];

        tr_x_xyyyz_xxxzz[i] = 2.0 * tr_x_xyz_xxxzz[i] * fe_0 + tr_x_xyyz_xxxzz[i] * pa_y[i];

        tr_x_xyyyz_xxyyy[i] = tr_x_xyyy_xxyyy[i] * pa_z[i];

        tr_x_xyyyz_xxyyz[i] = tr_x_xyyy_xxyy[i] * fe_0 + tr_x_xyyy_xxyyz[i] * pa_z[i];

        tr_x_xyyyz_xxyzz[i] = 2.0 * tr_x_xyyy_xxyz[i] * fe_0 + tr_x_xyyy_xxyzz[i] * pa_z[i];

        tr_x_xyyyz_xxzzz[i] = 2.0 * tr_x_xyz_xxzzz[i] * fe_0 + tr_x_xyyz_xxzzz[i] * pa_y[i];

        tr_x_xyyyz_xyyyy[i] = tr_x_xyyy_xyyyy[i] * pa_z[i];

        tr_x_xyyyz_xyyyz[i] = tr_x_xyyy_xyyy[i] * fe_0 + tr_x_xyyy_xyyyz[i] * pa_z[i];

        tr_x_xyyyz_xyyzz[i] = 2.0 * tr_x_xyyy_xyyz[i] * fe_0 + tr_x_xyyy_xyyzz[i] * pa_z[i];

        tr_x_xyyyz_xyzzz[i] = 3.0 * tr_x_xyyy_xyzz[i] * fe_0 + tr_x_xyyy_xyzzz[i] * pa_z[i];

        tr_x_xyyyz_xzzzz[i] = 2.0 * tr_x_xyz_xzzzz[i] * fe_0 + tr_x_xyyz_xzzzz[i] * pa_y[i];

        tr_x_xyyyz_yyyyy[i] = tr_x_xyyy_yyyyy[i] * pa_z[i];

        tr_x_xyyyz_yyyyz[i] = ts_yyyz_yyyyz[i] * fe_0 + tr_x_yyyz_yyyyz[i] * pa_x[i];

        tr_x_xyyyz_yyyzz[i] = ts_yyyz_yyyzz[i] * fe_0 + tr_x_yyyz_yyyzz[i] * pa_x[i];

        tr_x_xyyyz_yyzzz[i] = ts_yyyz_yyzzz[i] * fe_0 + tr_x_yyyz_yyzzz[i] * pa_x[i];

        tr_x_xyyyz_yzzzz[i] = ts_yyyz_yzzzz[i] * fe_0 + tr_x_yyyz_yzzzz[i] * pa_x[i];

        tr_x_xyyyz_zzzzz[i] = ts_yyyz_zzzzz[i] * fe_0 + tr_x_yyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 252-273 components of targeted buffer : HH

    auto tr_x_xyyzz_xxxxx = pbuffer.data(idx_dip_hh + 252);

    auto tr_x_xyyzz_xxxxy = pbuffer.data(idx_dip_hh + 253);

    auto tr_x_xyyzz_xxxxz = pbuffer.data(idx_dip_hh + 254);

    auto tr_x_xyyzz_xxxyy = pbuffer.data(idx_dip_hh + 255);

    auto tr_x_xyyzz_xxxyz = pbuffer.data(idx_dip_hh + 256);

    auto tr_x_xyyzz_xxxzz = pbuffer.data(idx_dip_hh + 257);

    auto tr_x_xyyzz_xxyyy = pbuffer.data(idx_dip_hh + 258);

    auto tr_x_xyyzz_xxyyz = pbuffer.data(idx_dip_hh + 259);

    auto tr_x_xyyzz_xxyzz = pbuffer.data(idx_dip_hh + 260);

    auto tr_x_xyyzz_xxzzz = pbuffer.data(idx_dip_hh + 261);

    auto tr_x_xyyzz_xyyyy = pbuffer.data(idx_dip_hh + 262);

    auto tr_x_xyyzz_xyyyz = pbuffer.data(idx_dip_hh + 263);

    auto tr_x_xyyzz_xyyzz = pbuffer.data(idx_dip_hh + 264);

    auto tr_x_xyyzz_xyzzz = pbuffer.data(idx_dip_hh + 265);

    auto tr_x_xyyzz_xzzzz = pbuffer.data(idx_dip_hh + 266);

    auto tr_x_xyyzz_yyyyy = pbuffer.data(idx_dip_hh + 267);

    auto tr_x_xyyzz_yyyyz = pbuffer.data(idx_dip_hh + 268);

    auto tr_x_xyyzz_yyyzz = pbuffer.data(idx_dip_hh + 269);

    auto tr_x_xyyzz_yyzzz = pbuffer.data(idx_dip_hh + 270);

    auto tr_x_xyyzz_yzzzz = pbuffer.data(idx_dip_hh + 271);

    auto tr_x_xyyzz_zzzzz = pbuffer.data(idx_dip_hh + 272);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_x_xyy_xxxxy,   \
                             tr_x_xyy_xxxyy,   \
                             tr_x_xyy_xxyyy,   \
                             tr_x_xyy_xyyyy,   \
                             tr_x_xyyz_xxxxy,  \
                             tr_x_xyyz_xxxyy,  \
                             tr_x_xyyz_xxyyy,  \
                             tr_x_xyyz_xyyyy,  \
                             tr_x_xyyzz_xxxxx, \
                             tr_x_xyyzz_xxxxy, \
                             tr_x_xyyzz_xxxxz, \
                             tr_x_xyyzz_xxxyy, \
                             tr_x_xyyzz_xxxyz, \
                             tr_x_xyyzz_xxxzz, \
                             tr_x_xyyzz_xxyyy, \
                             tr_x_xyyzz_xxyyz, \
                             tr_x_xyyzz_xxyzz, \
                             tr_x_xyyzz_xxzzz, \
                             tr_x_xyyzz_xyyyy, \
                             tr_x_xyyzz_xyyyz, \
                             tr_x_xyyzz_xyyzz, \
                             tr_x_xyyzz_xyzzz, \
                             tr_x_xyyzz_xzzzz, \
                             tr_x_xyyzz_yyyyy, \
                             tr_x_xyyzz_yyyyz, \
                             tr_x_xyyzz_yyyzz, \
                             tr_x_xyyzz_yyzzz, \
                             tr_x_xyyzz_yzzzz, \
                             tr_x_xyyzz_zzzzz, \
                             tr_x_xyzz_xxxxx,  \
                             tr_x_xyzz_xxxxz,  \
                             tr_x_xyzz_xxxzz,  \
                             tr_x_xyzz_xxzzz,  \
                             tr_x_xyzz_xzzzz,  \
                             tr_x_xzz_xxxxx,   \
                             tr_x_xzz_xxxxz,   \
                             tr_x_xzz_xxxzz,   \
                             tr_x_xzz_xxzzz,   \
                             tr_x_xzz_xzzzz,   \
                             tr_x_yyzz_xxxyz,  \
                             tr_x_yyzz_xxyyz,  \
                             tr_x_yyzz_xxyz,   \
                             tr_x_yyzz_xxyzz,  \
                             tr_x_yyzz_xyyyz,  \
                             tr_x_yyzz_xyyz,   \
                             tr_x_yyzz_xyyzz,  \
                             tr_x_yyzz_xyzz,   \
                             tr_x_yyzz_xyzzz,  \
                             tr_x_yyzz_yyyyy,  \
                             tr_x_yyzz_yyyyz,  \
                             tr_x_yyzz_yyyz,   \
                             tr_x_yyzz_yyyzz,  \
                             tr_x_yyzz_yyzz,   \
                             tr_x_yyzz_yyzzz,  \
                             tr_x_yyzz_yzzz,   \
                             tr_x_yyzz_yzzzz,  \
                             tr_x_yyzz_zzzzz,  \
                             ts_yyzz_xxxyz,    \
                             ts_yyzz_xxyyz,    \
                             ts_yyzz_xxyzz,    \
                             ts_yyzz_xyyyz,    \
                             ts_yyzz_xyyzz,    \
                             ts_yyzz_xyzzz,    \
                             ts_yyzz_yyyyy,    \
                             ts_yyzz_yyyyz,    \
                             ts_yyzz_yyyzz,    \
                             ts_yyzz_yyzzz,    \
                             ts_yyzz_yzzzz,    \
                             ts_yyzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzz_xxxxx[i] = tr_x_xzz_xxxxx[i] * fe_0 + tr_x_xyzz_xxxxx[i] * pa_y[i];

        tr_x_xyyzz_xxxxy[i] = tr_x_xyy_xxxxy[i] * fe_0 + tr_x_xyyz_xxxxy[i] * pa_z[i];

        tr_x_xyyzz_xxxxz[i] = tr_x_xzz_xxxxz[i] * fe_0 + tr_x_xyzz_xxxxz[i] * pa_y[i];

        tr_x_xyyzz_xxxyy[i] = tr_x_xyy_xxxyy[i] * fe_0 + tr_x_xyyz_xxxyy[i] * pa_z[i];

        tr_x_xyyzz_xxxyz[i] = 3.0 * tr_x_yyzz_xxyz[i] * fe_0 + ts_yyzz_xxxyz[i] * fe_0 + tr_x_yyzz_xxxyz[i] * pa_x[i];

        tr_x_xyyzz_xxxzz[i] = tr_x_xzz_xxxzz[i] * fe_0 + tr_x_xyzz_xxxzz[i] * pa_y[i];

        tr_x_xyyzz_xxyyy[i] = tr_x_xyy_xxyyy[i] * fe_0 + tr_x_xyyz_xxyyy[i] * pa_z[i];

        tr_x_xyyzz_xxyyz[i] = 2.0 * tr_x_yyzz_xyyz[i] * fe_0 + ts_yyzz_xxyyz[i] * fe_0 + tr_x_yyzz_xxyyz[i] * pa_x[i];

        tr_x_xyyzz_xxyzz[i] = 2.0 * tr_x_yyzz_xyzz[i] * fe_0 + ts_yyzz_xxyzz[i] * fe_0 + tr_x_yyzz_xxyzz[i] * pa_x[i];

        tr_x_xyyzz_xxzzz[i] = tr_x_xzz_xxzzz[i] * fe_0 + tr_x_xyzz_xxzzz[i] * pa_y[i];

        tr_x_xyyzz_xyyyy[i] = tr_x_xyy_xyyyy[i] * fe_0 + tr_x_xyyz_xyyyy[i] * pa_z[i];

        tr_x_xyyzz_xyyyz[i] = tr_x_yyzz_yyyz[i] * fe_0 + ts_yyzz_xyyyz[i] * fe_0 + tr_x_yyzz_xyyyz[i] * pa_x[i];

        tr_x_xyyzz_xyyzz[i] = tr_x_yyzz_yyzz[i] * fe_0 + ts_yyzz_xyyzz[i] * fe_0 + tr_x_yyzz_xyyzz[i] * pa_x[i];

        tr_x_xyyzz_xyzzz[i] = tr_x_yyzz_yzzz[i] * fe_0 + ts_yyzz_xyzzz[i] * fe_0 + tr_x_yyzz_xyzzz[i] * pa_x[i];

        tr_x_xyyzz_xzzzz[i] = tr_x_xzz_xzzzz[i] * fe_0 + tr_x_xyzz_xzzzz[i] * pa_y[i];

        tr_x_xyyzz_yyyyy[i] = ts_yyzz_yyyyy[i] * fe_0 + tr_x_yyzz_yyyyy[i] * pa_x[i];

        tr_x_xyyzz_yyyyz[i] = ts_yyzz_yyyyz[i] * fe_0 + tr_x_yyzz_yyyyz[i] * pa_x[i];

        tr_x_xyyzz_yyyzz[i] = ts_yyzz_yyyzz[i] * fe_0 + tr_x_yyzz_yyyzz[i] * pa_x[i];

        tr_x_xyyzz_yyzzz[i] = ts_yyzz_yyzzz[i] * fe_0 + tr_x_yyzz_yyzzz[i] * pa_x[i];

        tr_x_xyyzz_yzzzz[i] = ts_yyzz_yzzzz[i] * fe_0 + tr_x_yyzz_yzzzz[i] * pa_x[i];

        tr_x_xyyzz_zzzzz[i] = ts_yyzz_zzzzz[i] * fe_0 + tr_x_yyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : HH

    auto tr_x_xyzzz_xxxxx = pbuffer.data(idx_dip_hh + 273);

    auto tr_x_xyzzz_xxxxy = pbuffer.data(idx_dip_hh + 274);

    auto tr_x_xyzzz_xxxxz = pbuffer.data(idx_dip_hh + 275);

    auto tr_x_xyzzz_xxxyy = pbuffer.data(idx_dip_hh + 276);

    auto tr_x_xyzzz_xxxyz = pbuffer.data(idx_dip_hh + 277);

    auto tr_x_xyzzz_xxxzz = pbuffer.data(idx_dip_hh + 278);

    auto tr_x_xyzzz_xxyyy = pbuffer.data(idx_dip_hh + 279);

    auto tr_x_xyzzz_xxyyz = pbuffer.data(idx_dip_hh + 280);

    auto tr_x_xyzzz_xxyzz = pbuffer.data(idx_dip_hh + 281);

    auto tr_x_xyzzz_xxzzz = pbuffer.data(idx_dip_hh + 282);

    auto tr_x_xyzzz_xyyyy = pbuffer.data(idx_dip_hh + 283);

    auto tr_x_xyzzz_xyyyz = pbuffer.data(idx_dip_hh + 284);

    auto tr_x_xyzzz_xyyzz = pbuffer.data(idx_dip_hh + 285);

    auto tr_x_xyzzz_xyzzz = pbuffer.data(idx_dip_hh + 286);

    auto tr_x_xyzzz_xzzzz = pbuffer.data(idx_dip_hh + 287);

    auto tr_x_xyzzz_yyyyy = pbuffer.data(idx_dip_hh + 288);

    auto tr_x_xyzzz_yyyyz = pbuffer.data(idx_dip_hh + 289);

    auto tr_x_xyzzz_yyyzz = pbuffer.data(idx_dip_hh + 290);

    auto tr_x_xyzzz_yyzzz = pbuffer.data(idx_dip_hh + 291);

    auto tr_x_xyzzz_yzzzz = pbuffer.data(idx_dip_hh + 292);

    auto tr_x_xyzzz_zzzzz = pbuffer.data(idx_dip_hh + 293);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_x_xyzzz_xxxxx, \
                             tr_x_xyzzz_xxxxy, \
                             tr_x_xyzzz_xxxxz, \
                             tr_x_xyzzz_xxxyy, \
                             tr_x_xyzzz_xxxyz, \
                             tr_x_xyzzz_xxxzz, \
                             tr_x_xyzzz_xxyyy, \
                             tr_x_xyzzz_xxyyz, \
                             tr_x_xyzzz_xxyzz, \
                             tr_x_xyzzz_xxzzz, \
                             tr_x_xyzzz_xyyyy, \
                             tr_x_xyzzz_xyyyz, \
                             tr_x_xyzzz_xyyzz, \
                             tr_x_xyzzz_xyzzz, \
                             tr_x_xyzzz_xzzzz, \
                             tr_x_xyzzz_yyyyy, \
                             tr_x_xyzzz_yyyyz, \
                             tr_x_xyzzz_yyyzz, \
                             tr_x_xyzzz_yyzzz, \
                             tr_x_xyzzz_yzzzz, \
                             tr_x_xyzzz_zzzzz, \
                             tr_x_xzzz_xxxx,   \
                             tr_x_xzzz_xxxxx,  \
                             tr_x_xzzz_xxxxy,  \
                             tr_x_xzzz_xxxxz,  \
                             tr_x_xzzz_xxxy,   \
                             tr_x_xzzz_xxxyy,  \
                             tr_x_xzzz_xxxyz,  \
                             tr_x_xzzz_xxxz,   \
                             tr_x_xzzz_xxxzz,  \
                             tr_x_xzzz_xxyy,   \
                             tr_x_xzzz_xxyyy,  \
                             tr_x_xzzz_xxyyz,  \
                             tr_x_xzzz_xxyz,   \
                             tr_x_xzzz_xxyzz,  \
                             tr_x_xzzz_xxzz,   \
                             tr_x_xzzz_xxzzz,  \
                             tr_x_xzzz_xyyy,   \
                             tr_x_xzzz_xyyyy,  \
                             tr_x_xzzz_xyyyz,  \
                             tr_x_xzzz_xyyz,   \
                             tr_x_xzzz_xyyzz,  \
                             tr_x_xzzz_xyzz,   \
                             tr_x_xzzz_xyzzz,  \
                             tr_x_xzzz_xzzz,   \
                             tr_x_xzzz_xzzzz,  \
                             tr_x_xzzz_zzzzz,  \
                             tr_x_yzzz_yyyyy,  \
                             tr_x_yzzz_yyyyz,  \
                             tr_x_yzzz_yyyzz,  \
                             tr_x_yzzz_yyzzz,  \
                             tr_x_yzzz_yzzzz,  \
                             ts_yzzz_yyyyy,    \
                             ts_yzzz_yyyyz,    \
                             ts_yzzz_yyyzz,    \
                             ts_yzzz_yyzzz,    \
                             ts_yzzz_yzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzz_xxxxx[i] = tr_x_xzzz_xxxxx[i] * pa_y[i];

        tr_x_xyzzz_xxxxy[i] = tr_x_xzzz_xxxx[i] * fe_0 + tr_x_xzzz_xxxxy[i] * pa_y[i];

        tr_x_xyzzz_xxxxz[i] = tr_x_xzzz_xxxxz[i] * pa_y[i];

        tr_x_xyzzz_xxxyy[i] = 2.0 * tr_x_xzzz_xxxy[i] * fe_0 + tr_x_xzzz_xxxyy[i] * pa_y[i];

        tr_x_xyzzz_xxxyz[i] = tr_x_xzzz_xxxz[i] * fe_0 + tr_x_xzzz_xxxyz[i] * pa_y[i];

        tr_x_xyzzz_xxxzz[i] = tr_x_xzzz_xxxzz[i] * pa_y[i];

        tr_x_xyzzz_xxyyy[i] = 3.0 * tr_x_xzzz_xxyy[i] * fe_0 + tr_x_xzzz_xxyyy[i] * pa_y[i];

        tr_x_xyzzz_xxyyz[i] = 2.0 * tr_x_xzzz_xxyz[i] * fe_0 + tr_x_xzzz_xxyyz[i] * pa_y[i];

        tr_x_xyzzz_xxyzz[i] = tr_x_xzzz_xxzz[i] * fe_0 + tr_x_xzzz_xxyzz[i] * pa_y[i];

        tr_x_xyzzz_xxzzz[i] = tr_x_xzzz_xxzzz[i] * pa_y[i];

        tr_x_xyzzz_xyyyy[i] = 4.0 * tr_x_xzzz_xyyy[i] * fe_0 + tr_x_xzzz_xyyyy[i] * pa_y[i];

        tr_x_xyzzz_xyyyz[i] = 3.0 * tr_x_xzzz_xyyz[i] * fe_0 + tr_x_xzzz_xyyyz[i] * pa_y[i];

        tr_x_xyzzz_xyyzz[i] = 2.0 * tr_x_xzzz_xyzz[i] * fe_0 + tr_x_xzzz_xyyzz[i] * pa_y[i];

        tr_x_xyzzz_xyzzz[i] = tr_x_xzzz_xzzz[i] * fe_0 + tr_x_xzzz_xyzzz[i] * pa_y[i];

        tr_x_xyzzz_xzzzz[i] = tr_x_xzzz_xzzzz[i] * pa_y[i];

        tr_x_xyzzz_yyyyy[i] = ts_yzzz_yyyyy[i] * fe_0 + tr_x_yzzz_yyyyy[i] * pa_x[i];

        tr_x_xyzzz_yyyyz[i] = ts_yzzz_yyyyz[i] * fe_0 + tr_x_yzzz_yyyyz[i] * pa_x[i];

        tr_x_xyzzz_yyyzz[i] = ts_yzzz_yyyzz[i] * fe_0 + tr_x_yzzz_yyyzz[i] * pa_x[i];

        tr_x_xyzzz_yyzzz[i] = ts_yzzz_yyzzz[i] * fe_0 + tr_x_yzzz_yyzzz[i] * pa_x[i];

        tr_x_xyzzz_yzzzz[i] = ts_yzzz_yzzzz[i] * fe_0 + tr_x_yzzz_yzzzz[i] * pa_x[i];

        tr_x_xyzzz_zzzzz[i] = tr_x_xzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : HH

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

    auto tr_x_xzzzz_yyyyy = pbuffer.data(idx_dip_hh + 309);

    auto tr_x_xzzzz_yyyyz = pbuffer.data(idx_dip_hh + 310);

    auto tr_x_xzzzz_yyyzz = pbuffer.data(idx_dip_hh + 311);

    auto tr_x_xzzzz_yyzzz = pbuffer.data(idx_dip_hh + 312);

    auto tr_x_xzzzz_yzzzz = pbuffer.data(idx_dip_hh + 313);

    auto tr_x_xzzzz_zzzzz = pbuffer.data(idx_dip_hh + 314);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_x_xzz_xxxxx,   \
                             tr_x_xzz_xxxxy,   \
                             tr_x_xzz_xxxyy,   \
                             tr_x_xzz_xxyyy,   \
                             tr_x_xzz_xyyyy,   \
                             tr_x_xzzz_xxxxx,  \
                             tr_x_xzzz_xxxxy,  \
                             tr_x_xzzz_xxxyy,  \
                             tr_x_xzzz_xxyyy,  \
                             tr_x_xzzz_xyyyy,  \
                             tr_x_xzzzz_xxxxx, \
                             tr_x_xzzzz_xxxxy, \
                             tr_x_xzzzz_xxxxz, \
                             tr_x_xzzzz_xxxyy, \
                             tr_x_xzzzz_xxxyz, \
                             tr_x_xzzzz_xxxzz, \
                             tr_x_xzzzz_xxyyy, \
                             tr_x_xzzzz_xxyyz, \
                             tr_x_xzzzz_xxyzz, \
                             tr_x_xzzzz_xxzzz, \
                             tr_x_xzzzz_xyyyy, \
                             tr_x_xzzzz_xyyyz, \
                             tr_x_xzzzz_xyyzz, \
                             tr_x_xzzzz_xyzzz, \
                             tr_x_xzzzz_xzzzz, \
                             tr_x_xzzzz_yyyyy, \
                             tr_x_xzzzz_yyyyz, \
                             tr_x_xzzzz_yyyzz, \
                             tr_x_xzzzz_yyzzz, \
                             tr_x_xzzzz_yzzzz, \
                             tr_x_xzzzz_zzzzz, \
                             tr_x_zzzz_xxxxz,  \
                             tr_x_zzzz_xxxyz,  \
                             tr_x_zzzz_xxxz,   \
                             tr_x_zzzz_xxxzz,  \
                             tr_x_zzzz_xxyyz,  \
                             tr_x_zzzz_xxyz,   \
                             tr_x_zzzz_xxyzz,  \
                             tr_x_zzzz_xxzz,   \
                             tr_x_zzzz_xxzzz,  \
                             tr_x_zzzz_xyyyz,  \
                             tr_x_zzzz_xyyz,   \
                             tr_x_zzzz_xyyzz,  \
                             tr_x_zzzz_xyzz,   \
                             tr_x_zzzz_xyzzz,  \
                             tr_x_zzzz_xzzz,   \
                             tr_x_zzzz_xzzzz,  \
                             tr_x_zzzz_yyyyy,  \
                             tr_x_zzzz_yyyyz,  \
                             tr_x_zzzz_yyyz,   \
                             tr_x_zzzz_yyyzz,  \
                             tr_x_zzzz_yyzz,   \
                             tr_x_zzzz_yyzzz,  \
                             tr_x_zzzz_yzzz,   \
                             tr_x_zzzz_yzzzz,  \
                             tr_x_zzzz_zzzz,   \
                             tr_x_zzzz_zzzzz,  \
                             ts_zzzz_xxxxz,    \
                             ts_zzzz_xxxyz,    \
                             ts_zzzz_xxxzz,    \
                             ts_zzzz_xxyyz,    \
                             ts_zzzz_xxyzz,    \
                             ts_zzzz_xxzzz,    \
                             ts_zzzz_xyyyz,    \
                             ts_zzzz_xyyzz,    \
                             ts_zzzz_xyzzz,    \
                             ts_zzzz_xzzzz,    \
                             ts_zzzz_yyyyy,    \
                             ts_zzzz_yyyyz,    \
                             ts_zzzz_yyyzz,    \
                             ts_zzzz_yyzzz,    \
                             ts_zzzz_yzzzz,    \
                             ts_zzzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzz_xxxxx[i] = 3.0 * tr_x_xzz_xxxxx[i] * fe_0 + tr_x_xzzz_xxxxx[i] * pa_z[i];

        tr_x_xzzzz_xxxxy[i] = 3.0 * tr_x_xzz_xxxxy[i] * fe_0 + tr_x_xzzz_xxxxy[i] * pa_z[i];

        tr_x_xzzzz_xxxxz[i] = 4.0 * tr_x_zzzz_xxxz[i] * fe_0 + ts_zzzz_xxxxz[i] * fe_0 + tr_x_zzzz_xxxxz[i] * pa_x[i];

        tr_x_xzzzz_xxxyy[i] = 3.0 * tr_x_xzz_xxxyy[i] * fe_0 + tr_x_xzzz_xxxyy[i] * pa_z[i];

        tr_x_xzzzz_xxxyz[i] = 3.0 * tr_x_zzzz_xxyz[i] * fe_0 + ts_zzzz_xxxyz[i] * fe_0 + tr_x_zzzz_xxxyz[i] * pa_x[i];

        tr_x_xzzzz_xxxzz[i] = 3.0 * tr_x_zzzz_xxzz[i] * fe_0 + ts_zzzz_xxxzz[i] * fe_0 + tr_x_zzzz_xxxzz[i] * pa_x[i];

        tr_x_xzzzz_xxyyy[i] = 3.0 * tr_x_xzz_xxyyy[i] * fe_0 + tr_x_xzzz_xxyyy[i] * pa_z[i];

        tr_x_xzzzz_xxyyz[i] = 2.0 * tr_x_zzzz_xyyz[i] * fe_0 + ts_zzzz_xxyyz[i] * fe_0 + tr_x_zzzz_xxyyz[i] * pa_x[i];

        tr_x_xzzzz_xxyzz[i] = 2.0 * tr_x_zzzz_xyzz[i] * fe_0 + ts_zzzz_xxyzz[i] * fe_0 + tr_x_zzzz_xxyzz[i] * pa_x[i];

        tr_x_xzzzz_xxzzz[i] = 2.0 * tr_x_zzzz_xzzz[i] * fe_0 + ts_zzzz_xxzzz[i] * fe_0 + tr_x_zzzz_xxzzz[i] * pa_x[i];

        tr_x_xzzzz_xyyyy[i] = 3.0 * tr_x_xzz_xyyyy[i] * fe_0 + tr_x_xzzz_xyyyy[i] * pa_z[i];

        tr_x_xzzzz_xyyyz[i] = tr_x_zzzz_yyyz[i] * fe_0 + ts_zzzz_xyyyz[i] * fe_0 + tr_x_zzzz_xyyyz[i] * pa_x[i];

        tr_x_xzzzz_xyyzz[i] = tr_x_zzzz_yyzz[i] * fe_0 + ts_zzzz_xyyzz[i] * fe_0 + tr_x_zzzz_xyyzz[i] * pa_x[i];

        tr_x_xzzzz_xyzzz[i] = tr_x_zzzz_yzzz[i] * fe_0 + ts_zzzz_xyzzz[i] * fe_0 + tr_x_zzzz_xyzzz[i] * pa_x[i];

        tr_x_xzzzz_xzzzz[i] = tr_x_zzzz_zzzz[i] * fe_0 + ts_zzzz_xzzzz[i] * fe_0 + tr_x_zzzz_xzzzz[i] * pa_x[i];

        tr_x_xzzzz_yyyyy[i] = ts_zzzz_yyyyy[i] * fe_0 + tr_x_zzzz_yyyyy[i] * pa_x[i];

        tr_x_xzzzz_yyyyz[i] = ts_zzzz_yyyyz[i] * fe_0 + tr_x_zzzz_yyyyz[i] * pa_x[i];

        tr_x_xzzzz_yyyzz[i] = ts_zzzz_yyyzz[i] * fe_0 + tr_x_zzzz_yyyzz[i] * pa_x[i];

        tr_x_xzzzz_yyzzz[i] = ts_zzzz_yyzzz[i] * fe_0 + tr_x_zzzz_yyzzz[i] * pa_x[i];

        tr_x_xzzzz_yzzzz[i] = ts_zzzz_yzzzz[i] * fe_0 + tr_x_zzzz_yzzzz[i] * pa_x[i];

        tr_x_xzzzz_zzzzz[i] = ts_zzzz_zzzzz[i] * fe_0 + tr_x_zzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             tr_x_yyy_xxxxx,   \
                             tr_x_yyy_xxxxy,   \
                             tr_x_yyy_xxxxz,   \
                             tr_x_yyy_xxxyy,   \
                             tr_x_yyy_xxxyz,   \
                             tr_x_yyy_xxxzz,   \
                             tr_x_yyy_xxyyy,   \
                             tr_x_yyy_xxyyz,   \
                             tr_x_yyy_xxyzz,   \
                             tr_x_yyy_xxzzz,   \
                             tr_x_yyy_xyyyy,   \
                             tr_x_yyy_xyyyz,   \
                             tr_x_yyy_xyyzz,   \
                             tr_x_yyy_xyzzz,   \
                             tr_x_yyy_xzzzz,   \
                             tr_x_yyy_yyyyy,   \
                             tr_x_yyy_yyyyz,   \
                             tr_x_yyy_yyyzz,   \
                             tr_x_yyy_yyzzz,   \
                             tr_x_yyy_yzzzz,   \
                             tr_x_yyy_zzzzz,   \
                             tr_x_yyyy_xxxx,   \
                             tr_x_yyyy_xxxxx,  \
                             tr_x_yyyy_xxxxy,  \
                             tr_x_yyyy_xxxxz,  \
                             tr_x_yyyy_xxxy,   \
                             tr_x_yyyy_xxxyy,  \
                             tr_x_yyyy_xxxyz,  \
                             tr_x_yyyy_xxxz,   \
                             tr_x_yyyy_xxxzz,  \
                             tr_x_yyyy_xxyy,   \
                             tr_x_yyyy_xxyyy,  \
                             tr_x_yyyy_xxyyz,  \
                             tr_x_yyyy_xxyz,   \
                             tr_x_yyyy_xxyzz,  \
                             tr_x_yyyy_xxzz,   \
                             tr_x_yyyy_xxzzz,  \
                             tr_x_yyyy_xyyy,   \
                             tr_x_yyyy_xyyyy,  \
                             tr_x_yyyy_xyyyz,  \
                             tr_x_yyyy_xyyz,   \
                             tr_x_yyyy_xyyzz,  \
                             tr_x_yyyy_xyzz,   \
                             tr_x_yyyy_xyzzz,  \
                             tr_x_yyyy_xzzz,   \
                             tr_x_yyyy_xzzzz,  \
                             tr_x_yyyy_yyyy,   \
                             tr_x_yyyy_yyyyy,  \
                             tr_x_yyyy_yyyyz,  \
                             tr_x_yyyy_yyyz,   \
                             tr_x_yyyy_yyyzz,  \
                             tr_x_yyyy_yyzz,   \
                             tr_x_yyyy_yyzzz,  \
                             tr_x_yyyy_yzzz,   \
                             tr_x_yyyy_yzzzz,  \
                             tr_x_yyyy_zzzz,   \
                             tr_x_yyyy_zzzzz,  \
                             tr_x_yyyyy_xxxxx, \
                             tr_x_yyyyy_xxxxy, \
                             tr_x_yyyyy_xxxxz, \
                             tr_x_yyyyy_xxxyy, \
                             tr_x_yyyyy_xxxyz, \
                             tr_x_yyyyy_xxxzz, \
                             tr_x_yyyyy_xxyyy, \
                             tr_x_yyyyy_xxyyz, \
                             tr_x_yyyyy_xxyzz, \
                             tr_x_yyyyy_xxzzz, \
                             tr_x_yyyyy_xyyyy, \
                             tr_x_yyyyy_xyyyz, \
                             tr_x_yyyyy_xyyzz, \
                             tr_x_yyyyy_xyzzz, \
                             tr_x_yyyyy_xzzzz, \
                             tr_x_yyyyy_yyyyy, \
                             tr_x_yyyyy_yyyyz, \
                             tr_x_yyyyy_yyyzz, \
                             tr_x_yyyyy_yyzzz, \
                             tr_x_yyyyy_yzzzz, \
                             tr_x_yyyyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyy_xxxxx[i] = 4.0 * tr_x_yyy_xxxxx[i] * fe_0 + tr_x_yyyy_xxxxx[i] * pa_y[i];

        tr_x_yyyyy_xxxxy[i] = 4.0 * tr_x_yyy_xxxxy[i] * fe_0 + tr_x_yyyy_xxxx[i] * fe_0 + tr_x_yyyy_xxxxy[i] * pa_y[i];

        tr_x_yyyyy_xxxxz[i] = 4.0 * tr_x_yyy_xxxxz[i] * fe_0 + tr_x_yyyy_xxxxz[i] * pa_y[i];

        tr_x_yyyyy_xxxyy[i] = 4.0 * tr_x_yyy_xxxyy[i] * fe_0 + 2.0 * tr_x_yyyy_xxxy[i] * fe_0 + tr_x_yyyy_xxxyy[i] * pa_y[i];

        tr_x_yyyyy_xxxyz[i] = 4.0 * tr_x_yyy_xxxyz[i] * fe_0 + tr_x_yyyy_xxxz[i] * fe_0 + tr_x_yyyy_xxxyz[i] * pa_y[i];

        tr_x_yyyyy_xxxzz[i] = 4.0 * tr_x_yyy_xxxzz[i] * fe_0 + tr_x_yyyy_xxxzz[i] * pa_y[i];

        tr_x_yyyyy_xxyyy[i] = 4.0 * tr_x_yyy_xxyyy[i] * fe_0 + 3.0 * tr_x_yyyy_xxyy[i] * fe_0 + tr_x_yyyy_xxyyy[i] * pa_y[i];

        tr_x_yyyyy_xxyyz[i] = 4.0 * tr_x_yyy_xxyyz[i] * fe_0 + 2.0 * tr_x_yyyy_xxyz[i] * fe_0 + tr_x_yyyy_xxyyz[i] * pa_y[i];

        tr_x_yyyyy_xxyzz[i] = 4.0 * tr_x_yyy_xxyzz[i] * fe_0 + tr_x_yyyy_xxzz[i] * fe_0 + tr_x_yyyy_xxyzz[i] * pa_y[i];

        tr_x_yyyyy_xxzzz[i] = 4.0 * tr_x_yyy_xxzzz[i] * fe_0 + tr_x_yyyy_xxzzz[i] * pa_y[i];

        tr_x_yyyyy_xyyyy[i] = 4.0 * tr_x_yyy_xyyyy[i] * fe_0 + 4.0 * tr_x_yyyy_xyyy[i] * fe_0 + tr_x_yyyy_xyyyy[i] * pa_y[i];

        tr_x_yyyyy_xyyyz[i] = 4.0 * tr_x_yyy_xyyyz[i] * fe_0 + 3.0 * tr_x_yyyy_xyyz[i] * fe_0 + tr_x_yyyy_xyyyz[i] * pa_y[i];

        tr_x_yyyyy_xyyzz[i] = 4.0 * tr_x_yyy_xyyzz[i] * fe_0 + 2.0 * tr_x_yyyy_xyzz[i] * fe_0 + tr_x_yyyy_xyyzz[i] * pa_y[i];

        tr_x_yyyyy_xyzzz[i] = 4.0 * tr_x_yyy_xyzzz[i] * fe_0 + tr_x_yyyy_xzzz[i] * fe_0 + tr_x_yyyy_xyzzz[i] * pa_y[i];

        tr_x_yyyyy_xzzzz[i] = 4.0 * tr_x_yyy_xzzzz[i] * fe_0 + tr_x_yyyy_xzzzz[i] * pa_y[i];

        tr_x_yyyyy_yyyyy[i] = 4.0 * tr_x_yyy_yyyyy[i] * fe_0 + 5.0 * tr_x_yyyy_yyyy[i] * fe_0 + tr_x_yyyy_yyyyy[i] * pa_y[i];

        tr_x_yyyyy_yyyyz[i] = 4.0 * tr_x_yyy_yyyyz[i] * fe_0 + 4.0 * tr_x_yyyy_yyyz[i] * fe_0 + tr_x_yyyy_yyyyz[i] * pa_y[i];

        tr_x_yyyyy_yyyzz[i] = 4.0 * tr_x_yyy_yyyzz[i] * fe_0 + 3.0 * tr_x_yyyy_yyzz[i] * fe_0 + tr_x_yyyy_yyyzz[i] * pa_y[i];

        tr_x_yyyyy_yyzzz[i] = 4.0 * tr_x_yyy_yyzzz[i] * fe_0 + 2.0 * tr_x_yyyy_yzzz[i] * fe_0 + tr_x_yyyy_yyzzz[i] * pa_y[i];

        tr_x_yyyyy_yzzzz[i] = 4.0 * tr_x_yyy_yzzzz[i] * fe_0 + tr_x_yyyy_zzzz[i] * fe_0 + tr_x_yyyy_yzzzz[i] * pa_y[i];

        tr_x_yyyyy_zzzzz[i] = 4.0 * tr_x_yyy_zzzzz[i] * fe_0 + tr_x_yyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 336-357 components of targeted buffer : HH

    auto tr_x_yyyyz_xxxxx = pbuffer.data(idx_dip_hh + 336);

    auto tr_x_yyyyz_xxxxy = pbuffer.data(idx_dip_hh + 337);

    auto tr_x_yyyyz_xxxxz = pbuffer.data(idx_dip_hh + 338);

    auto tr_x_yyyyz_xxxyy = pbuffer.data(idx_dip_hh + 339);

    auto tr_x_yyyyz_xxxyz = pbuffer.data(idx_dip_hh + 340);

    auto tr_x_yyyyz_xxxzz = pbuffer.data(idx_dip_hh + 341);

    auto tr_x_yyyyz_xxyyy = pbuffer.data(idx_dip_hh + 342);

    auto tr_x_yyyyz_xxyyz = pbuffer.data(idx_dip_hh + 343);

    auto tr_x_yyyyz_xxyzz = pbuffer.data(idx_dip_hh + 344);

    auto tr_x_yyyyz_xxzzz = pbuffer.data(idx_dip_hh + 345);

    auto tr_x_yyyyz_xyyyy = pbuffer.data(idx_dip_hh + 346);

    auto tr_x_yyyyz_xyyyz = pbuffer.data(idx_dip_hh + 347);

    auto tr_x_yyyyz_xyyzz = pbuffer.data(idx_dip_hh + 348);

    auto tr_x_yyyyz_xyzzz = pbuffer.data(idx_dip_hh + 349);

    auto tr_x_yyyyz_xzzzz = pbuffer.data(idx_dip_hh + 350);

    auto tr_x_yyyyz_yyyyy = pbuffer.data(idx_dip_hh + 351);

    auto tr_x_yyyyz_yyyyz = pbuffer.data(idx_dip_hh + 352);

    auto tr_x_yyyyz_yyyzz = pbuffer.data(idx_dip_hh + 353);

    auto tr_x_yyyyz_yyzzz = pbuffer.data(idx_dip_hh + 354);

    auto tr_x_yyyyz_yzzzz = pbuffer.data(idx_dip_hh + 355);

    auto tr_x_yyyyz_zzzzz = pbuffer.data(idx_dip_hh + 356);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyyy_xxxxx,  \
                             tr_x_yyyy_xxxxy,  \
                             tr_x_yyyy_xxxy,   \
                             tr_x_yyyy_xxxyy,  \
                             tr_x_yyyy_xxxyz,  \
                             tr_x_yyyy_xxyy,   \
                             tr_x_yyyy_xxyyy,  \
                             tr_x_yyyy_xxyyz,  \
                             tr_x_yyyy_xxyz,   \
                             tr_x_yyyy_xxyzz,  \
                             tr_x_yyyy_xyyy,   \
                             tr_x_yyyy_xyyyy,  \
                             tr_x_yyyy_xyyyz,  \
                             tr_x_yyyy_xyyz,   \
                             tr_x_yyyy_xyyzz,  \
                             tr_x_yyyy_xyzz,   \
                             tr_x_yyyy_xyzzz,  \
                             tr_x_yyyy_yyyy,   \
                             tr_x_yyyy_yyyyy,  \
                             tr_x_yyyy_yyyyz,  \
                             tr_x_yyyy_yyyz,   \
                             tr_x_yyyy_yyyzz,  \
                             tr_x_yyyy_yyzz,   \
                             tr_x_yyyy_yyzzz,  \
                             tr_x_yyyy_yzzz,   \
                             tr_x_yyyy_yzzzz,  \
                             tr_x_yyyyz_xxxxx, \
                             tr_x_yyyyz_xxxxy, \
                             tr_x_yyyyz_xxxxz, \
                             tr_x_yyyyz_xxxyy, \
                             tr_x_yyyyz_xxxyz, \
                             tr_x_yyyyz_xxxzz, \
                             tr_x_yyyyz_xxyyy, \
                             tr_x_yyyyz_xxyyz, \
                             tr_x_yyyyz_xxyzz, \
                             tr_x_yyyyz_xxzzz, \
                             tr_x_yyyyz_xyyyy, \
                             tr_x_yyyyz_xyyyz, \
                             tr_x_yyyyz_xyyzz, \
                             tr_x_yyyyz_xyzzz, \
                             tr_x_yyyyz_xzzzz, \
                             tr_x_yyyyz_yyyyy, \
                             tr_x_yyyyz_yyyyz, \
                             tr_x_yyyyz_yyyzz, \
                             tr_x_yyyyz_yyzzz, \
                             tr_x_yyyyz_yzzzz, \
                             tr_x_yyyyz_zzzzz, \
                             tr_x_yyyz_xxxxz,  \
                             tr_x_yyyz_xxxzz,  \
                             tr_x_yyyz_xxzzz,  \
                             tr_x_yyyz_xzzzz,  \
                             tr_x_yyyz_zzzzz,  \
                             tr_x_yyz_xxxxz,   \
                             tr_x_yyz_xxxzz,   \
                             tr_x_yyz_xxzzz,   \
                             tr_x_yyz_xzzzz,   \
                             tr_x_yyz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyz_xxxxx[i] = tr_x_yyyy_xxxxx[i] * pa_z[i];

        tr_x_yyyyz_xxxxy[i] = tr_x_yyyy_xxxxy[i] * pa_z[i];

        tr_x_yyyyz_xxxxz[i] = 3.0 * tr_x_yyz_xxxxz[i] * fe_0 + tr_x_yyyz_xxxxz[i] * pa_y[i];

        tr_x_yyyyz_xxxyy[i] = tr_x_yyyy_xxxyy[i] * pa_z[i];

        tr_x_yyyyz_xxxyz[i] = tr_x_yyyy_xxxy[i] * fe_0 + tr_x_yyyy_xxxyz[i] * pa_z[i];

        tr_x_yyyyz_xxxzz[i] = 3.0 * tr_x_yyz_xxxzz[i] * fe_0 + tr_x_yyyz_xxxzz[i] * pa_y[i];

        tr_x_yyyyz_xxyyy[i] = tr_x_yyyy_xxyyy[i] * pa_z[i];

        tr_x_yyyyz_xxyyz[i] = tr_x_yyyy_xxyy[i] * fe_0 + tr_x_yyyy_xxyyz[i] * pa_z[i];

        tr_x_yyyyz_xxyzz[i] = 2.0 * tr_x_yyyy_xxyz[i] * fe_0 + tr_x_yyyy_xxyzz[i] * pa_z[i];

        tr_x_yyyyz_xxzzz[i] = 3.0 * tr_x_yyz_xxzzz[i] * fe_0 + tr_x_yyyz_xxzzz[i] * pa_y[i];

        tr_x_yyyyz_xyyyy[i] = tr_x_yyyy_xyyyy[i] * pa_z[i];

        tr_x_yyyyz_xyyyz[i] = tr_x_yyyy_xyyy[i] * fe_0 + tr_x_yyyy_xyyyz[i] * pa_z[i];

        tr_x_yyyyz_xyyzz[i] = 2.0 * tr_x_yyyy_xyyz[i] * fe_0 + tr_x_yyyy_xyyzz[i] * pa_z[i];

        tr_x_yyyyz_xyzzz[i] = 3.0 * tr_x_yyyy_xyzz[i] * fe_0 + tr_x_yyyy_xyzzz[i] * pa_z[i];

        tr_x_yyyyz_xzzzz[i] = 3.0 * tr_x_yyz_xzzzz[i] * fe_0 + tr_x_yyyz_xzzzz[i] * pa_y[i];

        tr_x_yyyyz_yyyyy[i] = tr_x_yyyy_yyyyy[i] * pa_z[i];

        tr_x_yyyyz_yyyyz[i] = tr_x_yyyy_yyyy[i] * fe_0 + tr_x_yyyy_yyyyz[i] * pa_z[i];

        tr_x_yyyyz_yyyzz[i] = 2.0 * tr_x_yyyy_yyyz[i] * fe_0 + tr_x_yyyy_yyyzz[i] * pa_z[i];

        tr_x_yyyyz_yyzzz[i] = 3.0 * tr_x_yyyy_yyzz[i] * fe_0 + tr_x_yyyy_yyzzz[i] * pa_z[i];

        tr_x_yyyyz_yzzzz[i] = 4.0 * tr_x_yyyy_yzzz[i] * fe_0 + tr_x_yyyy_yzzzz[i] * pa_z[i];

        tr_x_yyyyz_zzzzz[i] = 3.0 * tr_x_yyz_zzzzz[i] * fe_0 + tr_x_yyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 357-378 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyy_xxxxy,   \
                             tr_x_yyy_xxxyy,   \
                             tr_x_yyy_xxyyy,   \
                             tr_x_yyy_xyyyy,   \
                             tr_x_yyy_yyyyy,   \
                             tr_x_yyyz_xxxxy,  \
                             tr_x_yyyz_xxxyy,  \
                             tr_x_yyyz_xxyyy,  \
                             tr_x_yyyz_xyyyy,  \
                             tr_x_yyyz_yyyyy,  \
                             tr_x_yyyzz_xxxxx, \
                             tr_x_yyyzz_xxxxy, \
                             tr_x_yyyzz_xxxxz, \
                             tr_x_yyyzz_xxxyy, \
                             tr_x_yyyzz_xxxyz, \
                             tr_x_yyyzz_xxxzz, \
                             tr_x_yyyzz_xxyyy, \
                             tr_x_yyyzz_xxyyz, \
                             tr_x_yyyzz_xxyzz, \
                             tr_x_yyyzz_xxzzz, \
                             tr_x_yyyzz_xyyyy, \
                             tr_x_yyyzz_xyyyz, \
                             tr_x_yyyzz_xyyzz, \
                             tr_x_yyyzz_xyzzz, \
                             tr_x_yyyzz_xzzzz, \
                             tr_x_yyyzz_yyyyy, \
                             tr_x_yyyzz_yyyyz, \
                             tr_x_yyyzz_yyyzz, \
                             tr_x_yyyzz_yyzzz, \
                             tr_x_yyyzz_yzzzz, \
                             tr_x_yyyzz_zzzzz, \
                             tr_x_yyzz_xxxxx,  \
                             tr_x_yyzz_xxxxz,  \
                             tr_x_yyzz_xxxyz,  \
                             tr_x_yyzz_xxxz,   \
                             tr_x_yyzz_xxxzz,  \
                             tr_x_yyzz_xxyyz,  \
                             tr_x_yyzz_xxyz,   \
                             tr_x_yyzz_xxyzz,  \
                             tr_x_yyzz_xxzz,   \
                             tr_x_yyzz_xxzzz,  \
                             tr_x_yyzz_xyyyz,  \
                             tr_x_yyzz_xyyz,   \
                             tr_x_yyzz_xyyzz,  \
                             tr_x_yyzz_xyzz,   \
                             tr_x_yyzz_xyzzz,  \
                             tr_x_yyzz_xzzz,   \
                             tr_x_yyzz_xzzzz,  \
                             tr_x_yyzz_yyyyz,  \
                             tr_x_yyzz_yyyz,   \
                             tr_x_yyzz_yyyzz,  \
                             tr_x_yyzz_yyzz,   \
                             tr_x_yyzz_yyzzz,  \
                             tr_x_yyzz_yzzz,   \
                             tr_x_yyzz_yzzzz,  \
                             tr_x_yyzz_zzzz,   \
                             tr_x_yyzz_zzzzz,  \
                             tr_x_yzz_xxxxx,   \
                             tr_x_yzz_xxxxz,   \
                             tr_x_yzz_xxxyz,   \
                             tr_x_yzz_xxxzz,   \
                             tr_x_yzz_xxyyz,   \
                             tr_x_yzz_xxyzz,   \
                             tr_x_yzz_xxzzz,   \
                             tr_x_yzz_xyyyz,   \
                             tr_x_yzz_xyyzz,   \
                             tr_x_yzz_xyzzz,   \
                             tr_x_yzz_xzzzz,   \
                             tr_x_yzz_yyyyz,   \
                             tr_x_yzz_yyyzz,   \
                             tr_x_yzz_yyzzz,   \
                             tr_x_yzz_yzzzz,   \
                             tr_x_yzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzz_xxxxx[i] = 2.0 * tr_x_yzz_xxxxx[i] * fe_0 + tr_x_yyzz_xxxxx[i] * pa_y[i];

        tr_x_yyyzz_xxxxy[i] = tr_x_yyy_xxxxy[i] * fe_0 + tr_x_yyyz_xxxxy[i] * pa_z[i];

        tr_x_yyyzz_xxxxz[i] = 2.0 * tr_x_yzz_xxxxz[i] * fe_0 + tr_x_yyzz_xxxxz[i] * pa_y[i];

        tr_x_yyyzz_xxxyy[i] = tr_x_yyy_xxxyy[i] * fe_0 + tr_x_yyyz_xxxyy[i] * pa_z[i];

        tr_x_yyyzz_xxxyz[i] = 2.0 * tr_x_yzz_xxxyz[i] * fe_0 + tr_x_yyzz_xxxz[i] * fe_0 + tr_x_yyzz_xxxyz[i] * pa_y[i];

        tr_x_yyyzz_xxxzz[i] = 2.0 * tr_x_yzz_xxxzz[i] * fe_0 + tr_x_yyzz_xxxzz[i] * pa_y[i];

        tr_x_yyyzz_xxyyy[i] = tr_x_yyy_xxyyy[i] * fe_0 + tr_x_yyyz_xxyyy[i] * pa_z[i];

        tr_x_yyyzz_xxyyz[i] = 2.0 * tr_x_yzz_xxyyz[i] * fe_0 + 2.0 * tr_x_yyzz_xxyz[i] * fe_0 + tr_x_yyzz_xxyyz[i] * pa_y[i];

        tr_x_yyyzz_xxyzz[i] = 2.0 * tr_x_yzz_xxyzz[i] * fe_0 + tr_x_yyzz_xxzz[i] * fe_0 + tr_x_yyzz_xxyzz[i] * pa_y[i];

        tr_x_yyyzz_xxzzz[i] = 2.0 * tr_x_yzz_xxzzz[i] * fe_0 + tr_x_yyzz_xxzzz[i] * pa_y[i];

        tr_x_yyyzz_xyyyy[i] = tr_x_yyy_xyyyy[i] * fe_0 + tr_x_yyyz_xyyyy[i] * pa_z[i];

        tr_x_yyyzz_xyyyz[i] = 2.0 * tr_x_yzz_xyyyz[i] * fe_0 + 3.0 * tr_x_yyzz_xyyz[i] * fe_0 + tr_x_yyzz_xyyyz[i] * pa_y[i];

        tr_x_yyyzz_xyyzz[i] = 2.0 * tr_x_yzz_xyyzz[i] * fe_0 + 2.0 * tr_x_yyzz_xyzz[i] * fe_0 + tr_x_yyzz_xyyzz[i] * pa_y[i];

        tr_x_yyyzz_xyzzz[i] = 2.0 * tr_x_yzz_xyzzz[i] * fe_0 + tr_x_yyzz_xzzz[i] * fe_0 + tr_x_yyzz_xyzzz[i] * pa_y[i];

        tr_x_yyyzz_xzzzz[i] = 2.0 * tr_x_yzz_xzzzz[i] * fe_0 + tr_x_yyzz_xzzzz[i] * pa_y[i];

        tr_x_yyyzz_yyyyy[i] = tr_x_yyy_yyyyy[i] * fe_0 + tr_x_yyyz_yyyyy[i] * pa_z[i];

        tr_x_yyyzz_yyyyz[i] = 2.0 * tr_x_yzz_yyyyz[i] * fe_0 + 4.0 * tr_x_yyzz_yyyz[i] * fe_0 + tr_x_yyzz_yyyyz[i] * pa_y[i];

        tr_x_yyyzz_yyyzz[i] = 2.0 * tr_x_yzz_yyyzz[i] * fe_0 + 3.0 * tr_x_yyzz_yyzz[i] * fe_0 + tr_x_yyzz_yyyzz[i] * pa_y[i];

        tr_x_yyyzz_yyzzz[i] = 2.0 * tr_x_yzz_yyzzz[i] * fe_0 + 2.0 * tr_x_yyzz_yzzz[i] * fe_0 + tr_x_yyzz_yyzzz[i] * pa_y[i];

        tr_x_yyyzz_yzzzz[i] = 2.0 * tr_x_yzz_yzzzz[i] * fe_0 + tr_x_yyzz_zzzz[i] * fe_0 + tr_x_yyzz_yzzzz[i] * pa_y[i];

        tr_x_yyyzz_zzzzz[i] = 2.0 * tr_x_yzz_zzzzz[i] * fe_0 + tr_x_yyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 378-399 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_x_yyz_xxxxy,   \
                             tr_x_yyz_xxxyy,   \
                             tr_x_yyz_xxyyy,   \
                             tr_x_yyz_xyyyy,   \
                             tr_x_yyz_yyyyy,   \
                             tr_x_yyzz_xxxxy,  \
                             tr_x_yyzz_xxxyy,  \
                             tr_x_yyzz_xxyyy,  \
                             tr_x_yyzz_xyyyy,  \
                             tr_x_yyzz_yyyyy,  \
                             tr_x_yyzzz_xxxxx, \
                             tr_x_yyzzz_xxxxy, \
                             tr_x_yyzzz_xxxxz, \
                             tr_x_yyzzz_xxxyy, \
                             tr_x_yyzzz_xxxyz, \
                             tr_x_yyzzz_xxxzz, \
                             tr_x_yyzzz_xxyyy, \
                             tr_x_yyzzz_xxyyz, \
                             tr_x_yyzzz_xxyzz, \
                             tr_x_yyzzz_xxzzz, \
                             tr_x_yyzzz_xyyyy, \
                             tr_x_yyzzz_xyyyz, \
                             tr_x_yyzzz_xyyzz, \
                             tr_x_yyzzz_xyzzz, \
                             tr_x_yyzzz_xzzzz, \
                             tr_x_yyzzz_yyyyy, \
                             tr_x_yyzzz_yyyyz, \
                             tr_x_yyzzz_yyyzz, \
                             tr_x_yyzzz_yyzzz, \
                             tr_x_yyzzz_yzzzz, \
                             tr_x_yyzzz_zzzzz, \
                             tr_x_yzzz_xxxxx,  \
                             tr_x_yzzz_xxxxz,  \
                             tr_x_yzzz_xxxyz,  \
                             tr_x_yzzz_xxxz,   \
                             tr_x_yzzz_xxxzz,  \
                             tr_x_yzzz_xxyyz,  \
                             tr_x_yzzz_xxyz,   \
                             tr_x_yzzz_xxyzz,  \
                             tr_x_yzzz_xxzz,   \
                             tr_x_yzzz_xxzzz,  \
                             tr_x_yzzz_xyyyz,  \
                             tr_x_yzzz_xyyz,   \
                             tr_x_yzzz_xyyzz,  \
                             tr_x_yzzz_xyzz,   \
                             tr_x_yzzz_xyzzz,  \
                             tr_x_yzzz_xzzz,   \
                             tr_x_yzzz_xzzzz,  \
                             tr_x_yzzz_yyyyz,  \
                             tr_x_yzzz_yyyz,   \
                             tr_x_yzzz_yyyzz,  \
                             tr_x_yzzz_yyzz,   \
                             tr_x_yzzz_yyzzz,  \
                             tr_x_yzzz_yzzz,   \
                             tr_x_yzzz_yzzzz,  \
                             tr_x_yzzz_zzzz,   \
                             tr_x_yzzz_zzzzz,  \
                             tr_x_zzz_xxxxx,   \
                             tr_x_zzz_xxxxz,   \
                             tr_x_zzz_xxxyz,   \
                             tr_x_zzz_xxxzz,   \
                             tr_x_zzz_xxyyz,   \
                             tr_x_zzz_xxyzz,   \
                             tr_x_zzz_xxzzz,   \
                             tr_x_zzz_xyyyz,   \
                             tr_x_zzz_xyyzz,   \
                             tr_x_zzz_xyzzz,   \
                             tr_x_zzz_xzzzz,   \
                             tr_x_zzz_yyyyz,   \
                             tr_x_zzz_yyyzz,   \
                             tr_x_zzz_yyzzz,   \
                             tr_x_zzz_yzzzz,   \
                             tr_x_zzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzz_xxxxx[i] = tr_x_zzz_xxxxx[i] * fe_0 + tr_x_yzzz_xxxxx[i] * pa_y[i];

        tr_x_yyzzz_xxxxy[i] = 2.0 * tr_x_yyz_xxxxy[i] * fe_0 + tr_x_yyzz_xxxxy[i] * pa_z[i];

        tr_x_yyzzz_xxxxz[i] = tr_x_zzz_xxxxz[i] * fe_0 + tr_x_yzzz_xxxxz[i] * pa_y[i];

        tr_x_yyzzz_xxxyy[i] = 2.0 * tr_x_yyz_xxxyy[i] * fe_0 + tr_x_yyzz_xxxyy[i] * pa_z[i];

        tr_x_yyzzz_xxxyz[i] = tr_x_zzz_xxxyz[i] * fe_0 + tr_x_yzzz_xxxz[i] * fe_0 + tr_x_yzzz_xxxyz[i] * pa_y[i];

        tr_x_yyzzz_xxxzz[i] = tr_x_zzz_xxxzz[i] * fe_0 + tr_x_yzzz_xxxzz[i] * pa_y[i];

        tr_x_yyzzz_xxyyy[i] = 2.0 * tr_x_yyz_xxyyy[i] * fe_0 + tr_x_yyzz_xxyyy[i] * pa_z[i];

        tr_x_yyzzz_xxyyz[i] = tr_x_zzz_xxyyz[i] * fe_0 + 2.0 * tr_x_yzzz_xxyz[i] * fe_0 + tr_x_yzzz_xxyyz[i] * pa_y[i];

        tr_x_yyzzz_xxyzz[i] = tr_x_zzz_xxyzz[i] * fe_0 + tr_x_yzzz_xxzz[i] * fe_0 + tr_x_yzzz_xxyzz[i] * pa_y[i];

        tr_x_yyzzz_xxzzz[i] = tr_x_zzz_xxzzz[i] * fe_0 + tr_x_yzzz_xxzzz[i] * pa_y[i];

        tr_x_yyzzz_xyyyy[i] = 2.0 * tr_x_yyz_xyyyy[i] * fe_0 + tr_x_yyzz_xyyyy[i] * pa_z[i];

        tr_x_yyzzz_xyyyz[i] = tr_x_zzz_xyyyz[i] * fe_0 + 3.0 * tr_x_yzzz_xyyz[i] * fe_0 + tr_x_yzzz_xyyyz[i] * pa_y[i];

        tr_x_yyzzz_xyyzz[i] = tr_x_zzz_xyyzz[i] * fe_0 + 2.0 * tr_x_yzzz_xyzz[i] * fe_0 + tr_x_yzzz_xyyzz[i] * pa_y[i];

        tr_x_yyzzz_xyzzz[i] = tr_x_zzz_xyzzz[i] * fe_0 + tr_x_yzzz_xzzz[i] * fe_0 + tr_x_yzzz_xyzzz[i] * pa_y[i];

        tr_x_yyzzz_xzzzz[i] = tr_x_zzz_xzzzz[i] * fe_0 + tr_x_yzzz_xzzzz[i] * pa_y[i];

        tr_x_yyzzz_yyyyy[i] = 2.0 * tr_x_yyz_yyyyy[i] * fe_0 + tr_x_yyzz_yyyyy[i] * pa_z[i];

        tr_x_yyzzz_yyyyz[i] = tr_x_zzz_yyyyz[i] * fe_0 + 4.0 * tr_x_yzzz_yyyz[i] * fe_0 + tr_x_yzzz_yyyyz[i] * pa_y[i];

        tr_x_yyzzz_yyyzz[i] = tr_x_zzz_yyyzz[i] * fe_0 + 3.0 * tr_x_yzzz_yyzz[i] * fe_0 + tr_x_yzzz_yyyzz[i] * pa_y[i];

        tr_x_yyzzz_yyzzz[i] = tr_x_zzz_yyzzz[i] * fe_0 + 2.0 * tr_x_yzzz_yzzz[i] * fe_0 + tr_x_yzzz_yyzzz[i] * pa_y[i];

        tr_x_yyzzz_yzzzz[i] = tr_x_zzz_yzzzz[i] * fe_0 + tr_x_yzzz_zzzz[i] * fe_0 + tr_x_yzzz_yzzzz[i] * pa_y[i];

        tr_x_yyzzz_zzzzz[i] = tr_x_zzz_zzzzz[i] * fe_0 + tr_x_yzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 399-420 components of targeted buffer : HH

    auto tr_x_yzzzz_xxxxx = pbuffer.data(idx_dip_hh + 399);

    auto tr_x_yzzzz_xxxxy = pbuffer.data(idx_dip_hh + 400);

    auto tr_x_yzzzz_xxxxz = pbuffer.data(idx_dip_hh + 401);

    auto tr_x_yzzzz_xxxyy = pbuffer.data(idx_dip_hh + 402);

    auto tr_x_yzzzz_xxxyz = pbuffer.data(idx_dip_hh + 403);

    auto tr_x_yzzzz_xxxzz = pbuffer.data(idx_dip_hh + 404);

    auto tr_x_yzzzz_xxyyy = pbuffer.data(idx_dip_hh + 405);

    auto tr_x_yzzzz_xxyyz = pbuffer.data(idx_dip_hh + 406);

    auto tr_x_yzzzz_xxyzz = pbuffer.data(idx_dip_hh + 407);

    auto tr_x_yzzzz_xxzzz = pbuffer.data(idx_dip_hh + 408);

    auto tr_x_yzzzz_xyyyy = pbuffer.data(idx_dip_hh + 409);

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

#pragma omp simd aligned(pa_y,                 \
                             tr_x_yzzzz_xxxxx, \
                             tr_x_yzzzz_xxxxy, \
                             tr_x_yzzzz_xxxxz, \
                             tr_x_yzzzz_xxxyy, \
                             tr_x_yzzzz_xxxyz, \
                             tr_x_yzzzz_xxxzz, \
                             tr_x_yzzzz_xxyyy, \
                             tr_x_yzzzz_xxyyz, \
                             tr_x_yzzzz_xxyzz, \
                             tr_x_yzzzz_xxzzz, \
                             tr_x_yzzzz_xyyyy, \
                             tr_x_yzzzz_xyyyz, \
                             tr_x_yzzzz_xyyzz, \
                             tr_x_yzzzz_xyzzz, \
                             tr_x_yzzzz_xzzzz, \
                             tr_x_yzzzz_yyyyy, \
                             tr_x_yzzzz_yyyyz, \
                             tr_x_yzzzz_yyyzz, \
                             tr_x_yzzzz_yyzzz, \
                             tr_x_yzzzz_yzzzz, \
                             tr_x_yzzzz_zzzzz, \
                             tr_x_zzzz_xxxx,   \
                             tr_x_zzzz_xxxxx,  \
                             tr_x_zzzz_xxxxy,  \
                             tr_x_zzzz_xxxxz,  \
                             tr_x_zzzz_xxxy,   \
                             tr_x_zzzz_xxxyy,  \
                             tr_x_zzzz_xxxyz,  \
                             tr_x_zzzz_xxxz,   \
                             tr_x_zzzz_xxxzz,  \
                             tr_x_zzzz_xxyy,   \
                             tr_x_zzzz_xxyyy,  \
                             tr_x_zzzz_xxyyz,  \
                             tr_x_zzzz_xxyz,   \
                             tr_x_zzzz_xxyzz,  \
                             tr_x_zzzz_xxzz,   \
                             tr_x_zzzz_xxzzz,  \
                             tr_x_zzzz_xyyy,   \
                             tr_x_zzzz_xyyyy,  \
                             tr_x_zzzz_xyyyz,  \
                             tr_x_zzzz_xyyz,   \
                             tr_x_zzzz_xyyzz,  \
                             tr_x_zzzz_xyzz,   \
                             tr_x_zzzz_xyzzz,  \
                             tr_x_zzzz_xzzz,   \
                             tr_x_zzzz_xzzzz,  \
                             tr_x_zzzz_yyyy,   \
                             tr_x_zzzz_yyyyy,  \
                             tr_x_zzzz_yyyyz,  \
                             tr_x_zzzz_yyyz,   \
                             tr_x_zzzz_yyyzz,  \
                             tr_x_zzzz_yyzz,   \
                             tr_x_zzzz_yyzzz,  \
                             tr_x_zzzz_yzzz,   \
                             tr_x_zzzz_yzzzz,  \
                             tr_x_zzzz_zzzz,   \
                             tr_x_zzzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzz_xxxxx[i] = tr_x_zzzz_xxxxx[i] * pa_y[i];

        tr_x_yzzzz_xxxxy[i] = tr_x_zzzz_xxxx[i] * fe_0 + tr_x_zzzz_xxxxy[i] * pa_y[i];

        tr_x_yzzzz_xxxxz[i] = tr_x_zzzz_xxxxz[i] * pa_y[i];

        tr_x_yzzzz_xxxyy[i] = 2.0 * tr_x_zzzz_xxxy[i] * fe_0 + tr_x_zzzz_xxxyy[i] * pa_y[i];

        tr_x_yzzzz_xxxyz[i] = tr_x_zzzz_xxxz[i] * fe_0 + tr_x_zzzz_xxxyz[i] * pa_y[i];

        tr_x_yzzzz_xxxzz[i] = tr_x_zzzz_xxxzz[i] * pa_y[i];

        tr_x_yzzzz_xxyyy[i] = 3.0 * tr_x_zzzz_xxyy[i] * fe_0 + tr_x_zzzz_xxyyy[i] * pa_y[i];

        tr_x_yzzzz_xxyyz[i] = 2.0 * tr_x_zzzz_xxyz[i] * fe_0 + tr_x_zzzz_xxyyz[i] * pa_y[i];

        tr_x_yzzzz_xxyzz[i] = tr_x_zzzz_xxzz[i] * fe_0 + tr_x_zzzz_xxyzz[i] * pa_y[i];

        tr_x_yzzzz_xxzzz[i] = tr_x_zzzz_xxzzz[i] * pa_y[i];

        tr_x_yzzzz_xyyyy[i] = 4.0 * tr_x_zzzz_xyyy[i] * fe_0 + tr_x_zzzz_xyyyy[i] * pa_y[i];

        tr_x_yzzzz_xyyyz[i] = 3.0 * tr_x_zzzz_xyyz[i] * fe_0 + tr_x_zzzz_xyyyz[i] * pa_y[i];

        tr_x_yzzzz_xyyzz[i] = 2.0 * tr_x_zzzz_xyzz[i] * fe_0 + tr_x_zzzz_xyyzz[i] * pa_y[i];

        tr_x_yzzzz_xyzzz[i] = tr_x_zzzz_xzzz[i] * fe_0 + tr_x_zzzz_xyzzz[i] * pa_y[i];

        tr_x_yzzzz_xzzzz[i] = tr_x_zzzz_xzzzz[i] * pa_y[i];

        tr_x_yzzzz_yyyyy[i] = 5.0 * tr_x_zzzz_yyyy[i] * fe_0 + tr_x_zzzz_yyyyy[i] * pa_y[i];

        tr_x_yzzzz_yyyyz[i] = 4.0 * tr_x_zzzz_yyyz[i] * fe_0 + tr_x_zzzz_yyyyz[i] * pa_y[i];

        tr_x_yzzzz_yyyzz[i] = 3.0 * tr_x_zzzz_yyzz[i] * fe_0 + tr_x_zzzz_yyyzz[i] * pa_y[i];

        tr_x_yzzzz_yyzzz[i] = 2.0 * tr_x_zzzz_yzzz[i] * fe_0 + tr_x_zzzz_yyzzz[i] * pa_y[i];

        tr_x_yzzzz_yzzzz[i] = tr_x_zzzz_zzzz[i] * fe_0 + tr_x_zzzz_yzzzz[i] * pa_y[i];

        tr_x_yzzzz_zzzzz[i] = tr_x_zzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 420-441 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_z,                 \
                             tr_x_zzz_xxxxx,   \
                             tr_x_zzz_xxxxy,   \
                             tr_x_zzz_xxxxz,   \
                             tr_x_zzz_xxxyy,   \
                             tr_x_zzz_xxxyz,   \
                             tr_x_zzz_xxxzz,   \
                             tr_x_zzz_xxyyy,   \
                             tr_x_zzz_xxyyz,   \
                             tr_x_zzz_xxyzz,   \
                             tr_x_zzz_xxzzz,   \
                             tr_x_zzz_xyyyy,   \
                             tr_x_zzz_xyyyz,   \
                             tr_x_zzz_xyyzz,   \
                             tr_x_zzz_xyzzz,   \
                             tr_x_zzz_xzzzz,   \
                             tr_x_zzz_yyyyy,   \
                             tr_x_zzz_yyyyz,   \
                             tr_x_zzz_yyyzz,   \
                             tr_x_zzz_yyzzz,   \
                             tr_x_zzz_yzzzz,   \
                             tr_x_zzz_zzzzz,   \
                             tr_x_zzzz_xxxx,   \
                             tr_x_zzzz_xxxxx,  \
                             tr_x_zzzz_xxxxy,  \
                             tr_x_zzzz_xxxxz,  \
                             tr_x_zzzz_xxxy,   \
                             tr_x_zzzz_xxxyy,  \
                             tr_x_zzzz_xxxyz,  \
                             tr_x_zzzz_xxxz,   \
                             tr_x_zzzz_xxxzz,  \
                             tr_x_zzzz_xxyy,   \
                             tr_x_zzzz_xxyyy,  \
                             tr_x_zzzz_xxyyz,  \
                             tr_x_zzzz_xxyz,   \
                             tr_x_zzzz_xxyzz,  \
                             tr_x_zzzz_xxzz,   \
                             tr_x_zzzz_xxzzz,  \
                             tr_x_zzzz_xyyy,   \
                             tr_x_zzzz_xyyyy,  \
                             tr_x_zzzz_xyyyz,  \
                             tr_x_zzzz_xyyz,   \
                             tr_x_zzzz_xyyzz,  \
                             tr_x_zzzz_xyzz,   \
                             tr_x_zzzz_xyzzz,  \
                             tr_x_zzzz_xzzz,   \
                             tr_x_zzzz_xzzzz,  \
                             tr_x_zzzz_yyyy,   \
                             tr_x_zzzz_yyyyy,  \
                             tr_x_zzzz_yyyyz,  \
                             tr_x_zzzz_yyyz,   \
                             tr_x_zzzz_yyyzz,  \
                             tr_x_zzzz_yyzz,   \
                             tr_x_zzzz_yyzzz,  \
                             tr_x_zzzz_yzzz,   \
                             tr_x_zzzz_yzzzz,  \
                             tr_x_zzzz_zzzz,   \
                             tr_x_zzzz_zzzzz,  \
                             tr_x_zzzzz_xxxxx, \
                             tr_x_zzzzz_xxxxy, \
                             tr_x_zzzzz_xxxxz, \
                             tr_x_zzzzz_xxxyy, \
                             tr_x_zzzzz_xxxyz, \
                             tr_x_zzzzz_xxxzz, \
                             tr_x_zzzzz_xxyyy, \
                             tr_x_zzzzz_xxyyz, \
                             tr_x_zzzzz_xxyzz, \
                             tr_x_zzzzz_xxzzz, \
                             tr_x_zzzzz_xyyyy, \
                             tr_x_zzzzz_xyyyz, \
                             tr_x_zzzzz_xyyzz, \
                             tr_x_zzzzz_xyzzz, \
                             tr_x_zzzzz_xzzzz, \
                             tr_x_zzzzz_yyyyy, \
                             tr_x_zzzzz_yyyyz, \
                             tr_x_zzzzz_yyyzz, \
                             tr_x_zzzzz_yyzzz, \
                             tr_x_zzzzz_yzzzz, \
                             tr_x_zzzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzz_xxxxx[i] = 4.0 * tr_x_zzz_xxxxx[i] * fe_0 + tr_x_zzzz_xxxxx[i] * pa_z[i];

        tr_x_zzzzz_xxxxy[i] = 4.0 * tr_x_zzz_xxxxy[i] * fe_0 + tr_x_zzzz_xxxxy[i] * pa_z[i];

        tr_x_zzzzz_xxxxz[i] = 4.0 * tr_x_zzz_xxxxz[i] * fe_0 + tr_x_zzzz_xxxx[i] * fe_0 + tr_x_zzzz_xxxxz[i] * pa_z[i];

        tr_x_zzzzz_xxxyy[i] = 4.0 * tr_x_zzz_xxxyy[i] * fe_0 + tr_x_zzzz_xxxyy[i] * pa_z[i];

        tr_x_zzzzz_xxxyz[i] = 4.0 * tr_x_zzz_xxxyz[i] * fe_0 + tr_x_zzzz_xxxy[i] * fe_0 + tr_x_zzzz_xxxyz[i] * pa_z[i];

        tr_x_zzzzz_xxxzz[i] = 4.0 * tr_x_zzz_xxxzz[i] * fe_0 + 2.0 * tr_x_zzzz_xxxz[i] * fe_0 + tr_x_zzzz_xxxzz[i] * pa_z[i];

        tr_x_zzzzz_xxyyy[i] = 4.0 * tr_x_zzz_xxyyy[i] * fe_0 + tr_x_zzzz_xxyyy[i] * pa_z[i];

        tr_x_zzzzz_xxyyz[i] = 4.0 * tr_x_zzz_xxyyz[i] * fe_0 + tr_x_zzzz_xxyy[i] * fe_0 + tr_x_zzzz_xxyyz[i] * pa_z[i];

        tr_x_zzzzz_xxyzz[i] = 4.0 * tr_x_zzz_xxyzz[i] * fe_0 + 2.0 * tr_x_zzzz_xxyz[i] * fe_0 + tr_x_zzzz_xxyzz[i] * pa_z[i];

        tr_x_zzzzz_xxzzz[i] = 4.0 * tr_x_zzz_xxzzz[i] * fe_0 + 3.0 * tr_x_zzzz_xxzz[i] * fe_0 + tr_x_zzzz_xxzzz[i] * pa_z[i];

        tr_x_zzzzz_xyyyy[i] = 4.0 * tr_x_zzz_xyyyy[i] * fe_0 + tr_x_zzzz_xyyyy[i] * pa_z[i];

        tr_x_zzzzz_xyyyz[i] = 4.0 * tr_x_zzz_xyyyz[i] * fe_0 + tr_x_zzzz_xyyy[i] * fe_0 + tr_x_zzzz_xyyyz[i] * pa_z[i];

        tr_x_zzzzz_xyyzz[i] = 4.0 * tr_x_zzz_xyyzz[i] * fe_0 + 2.0 * tr_x_zzzz_xyyz[i] * fe_0 + tr_x_zzzz_xyyzz[i] * pa_z[i];

        tr_x_zzzzz_xyzzz[i] = 4.0 * tr_x_zzz_xyzzz[i] * fe_0 + 3.0 * tr_x_zzzz_xyzz[i] * fe_0 + tr_x_zzzz_xyzzz[i] * pa_z[i];

        tr_x_zzzzz_xzzzz[i] = 4.0 * tr_x_zzz_xzzzz[i] * fe_0 + 4.0 * tr_x_zzzz_xzzz[i] * fe_0 + tr_x_zzzz_xzzzz[i] * pa_z[i];

        tr_x_zzzzz_yyyyy[i] = 4.0 * tr_x_zzz_yyyyy[i] * fe_0 + tr_x_zzzz_yyyyy[i] * pa_z[i];

        tr_x_zzzzz_yyyyz[i] = 4.0 * tr_x_zzz_yyyyz[i] * fe_0 + tr_x_zzzz_yyyy[i] * fe_0 + tr_x_zzzz_yyyyz[i] * pa_z[i];

        tr_x_zzzzz_yyyzz[i] = 4.0 * tr_x_zzz_yyyzz[i] * fe_0 + 2.0 * tr_x_zzzz_yyyz[i] * fe_0 + tr_x_zzzz_yyyzz[i] * pa_z[i];

        tr_x_zzzzz_yyzzz[i] = 4.0 * tr_x_zzz_yyzzz[i] * fe_0 + 3.0 * tr_x_zzzz_yyzz[i] * fe_0 + tr_x_zzzz_yyzzz[i] * pa_z[i];

        tr_x_zzzzz_yzzzz[i] = 4.0 * tr_x_zzz_yzzzz[i] * fe_0 + 4.0 * tr_x_zzzz_yzzz[i] * fe_0 + tr_x_zzzz_yzzzz[i] * pa_z[i];

        tr_x_zzzzz_zzzzz[i] = 4.0 * tr_x_zzz_zzzzz[i] * fe_0 + 5.0 * tr_x_zzzz_zzzz[i] * fe_0 + tr_x_zzzz_zzzzz[i] * pa_z[i];
    }

    // Set up 441-462 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxx_xxxxx,   \
                             tr_y_xxx_xxxxy,   \
                             tr_y_xxx_xxxxz,   \
                             tr_y_xxx_xxxyy,   \
                             tr_y_xxx_xxxyz,   \
                             tr_y_xxx_xxxzz,   \
                             tr_y_xxx_xxyyy,   \
                             tr_y_xxx_xxyyz,   \
                             tr_y_xxx_xxyzz,   \
                             tr_y_xxx_xxzzz,   \
                             tr_y_xxx_xyyyy,   \
                             tr_y_xxx_xyyyz,   \
                             tr_y_xxx_xyyzz,   \
                             tr_y_xxx_xyzzz,   \
                             tr_y_xxx_xzzzz,   \
                             tr_y_xxx_yyyyy,   \
                             tr_y_xxx_yyyyz,   \
                             tr_y_xxx_yyyzz,   \
                             tr_y_xxx_yyzzz,   \
                             tr_y_xxx_yzzzz,   \
                             tr_y_xxx_zzzzz,   \
                             tr_y_xxxx_xxxx,   \
                             tr_y_xxxx_xxxxx,  \
                             tr_y_xxxx_xxxxy,  \
                             tr_y_xxxx_xxxxz,  \
                             tr_y_xxxx_xxxy,   \
                             tr_y_xxxx_xxxyy,  \
                             tr_y_xxxx_xxxyz,  \
                             tr_y_xxxx_xxxz,   \
                             tr_y_xxxx_xxxzz,  \
                             tr_y_xxxx_xxyy,   \
                             tr_y_xxxx_xxyyy,  \
                             tr_y_xxxx_xxyyz,  \
                             tr_y_xxxx_xxyz,   \
                             tr_y_xxxx_xxyzz,  \
                             tr_y_xxxx_xxzz,   \
                             tr_y_xxxx_xxzzz,  \
                             tr_y_xxxx_xyyy,   \
                             tr_y_xxxx_xyyyy,  \
                             tr_y_xxxx_xyyyz,  \
                             tr_y_xxxx_xyyz,   \
                             tr_y_xxxx_xyyzz,  \
                             tr_y_xxxx_xyzz,   \
                             tr_y_xxxx_xyzzz,  \
                             tr_y_xxxx_xzzz,   \
                             tr_y_xxxx_xzzzz,  \
                             tr_y_xxxx_yyyy,   \
                             tr_y_xxxx_yyyyy,  \
                             tr_y_xxxx_yyyyz,  \
                             tr_y_xxxx_yyyz,   \
                             tr_y_xxxx_yyyzz,  \
                             tr_y_xxxx_yyzz,   \
                             tr_y_xxxx_yyzzz,  \
                             tr_y_xxxx_yzzz,   \
                             tr_y_xxxx_yzzzz,  \
                             tr_y_xxxx_zzzz,   \
                             tr_y_xxxx_zzzzz,  \
                             tr_y_xxxxx_xxxxx, \
                             tr_y_xxxxx_xxxxy, \
                             tr_y_xxxxx_xxxxz, \
                             tr_y_xxxxx_xxxyy, \
                             tr_y_xxxxx_xxxyz, \
                             tr_y_xxxxx_xxxzz, \
                             tr_y_xxxxx_xxyyy, \
                             tr_y_xxxxx_xxyyz, \
                             tr_y_xxxxx_xxyzz, \
                             tr_y_xxxxx_xxzzz, \
                             tr_y_xxxxx_xyyyy, \
                             tr_y_xxxxx_xyyyz, \
                             tr_y_xxxxx_xyyzz, \
                             tr_y_xxxxx_xyzzz, \
                             tr_y_xxxxx_xzzzz, \
                             tr_y_xxxxx_yyyyy, \
                             tr_y_xxxxx_yyyyz, \
                             tr_y_xxxxx_yyyzz, \
                             tr_y_xxxxx_yyzzz, \
                             tr_y_xxxxx_yzzzz, \
                             tr_y_xxxxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxx_xxxxx[i] = 4.0 * tr_y_xxx_xxxxx[i] * fe_0 + 5.0 * tr_y_xxxx_xxxx[i] * fe_0 + tr_y_xxxx_xxxxx[i] * pa_x[i];

        tr_y_xxxxx_xxxxy[i] = 4.0 * tr_y_xxx_xxxxy[i] * fe_0 + 4.0 * tr_y_xxxx_xxxy[i] * fe_0 + tr_y_xxxx_xxxxy[i] * pa_x[i];

        tr_y_xxxxx_xxxxz[i] = 4.0 * tr_y_xxx_xxxxz[i] * fe_0 + 4.0 * tr_y_xxxx_xxxz[i] * fe_0 + tr_y_xxxx_xxxxz[i] * pa_x[i];

        tr_y_xxxxx_xxxyy[i] = 4.0 * tr_y_xxx_xxxyy[i] * fe_0 + 3.0 * tr_y_xxxx_xxyy[i] * fe_0 + tr_y_xxxx_xxxyy[i] * pa_x[i];

        tr_y_xxxxx_xxxyz[i] = 4.0 * tr_y_xxx_xxxyz[i] * fe_0 + 3.0 * tr_y_xxxx_xxyz[i] * fe_0 + tr_y_xxxx_xxxyz[i] * pa_x[i];

        tr_y_xxxxx_xxxzz[i] = 4.0 * tr_y_xxx_xxxzz[i] * fe_0 + 3.0 * tr_y_xxxx_xxzz[i] * fe_0 + tr_y_xxxx_xxxzz[i] * pa_x[i];

        tr_y_xxxxx_xxyyy[i] = 4.0 * tr_y_xxx_xxyyy[i] * fe_0 + 2.0 * tr_y_xxxx_xyyy[i] * fe_0 + tr_y_xxxx_xxyyy[i] * pa_x[i];

        tr_y_xxxxx_xxyyz[i] = 4.0 * tr_y_xxx_xxyyz[i] * fe_0 + 2.0 * tr_y_xxxx_xyyz[i] * fe_0 + tr_y_xxxx_xxyyz[i] * pa_x[i];

        tr_y_xxxxx_xxyzz[i] = 4.0 * tr_y_xxx_xxyzz[i] * fe_0 + 2.0 * tr_y_xxxx_xyzz[i] * fe_0 + tr_y_xxxx_xxyzz[i] * pa_x[i];

        tr_y_xxxxx_xxzzz[i] = 4.0 * tr_y_xxx_xxzzz[i] * fe_0 + 2.0 * tr_y_xxxx_xzzz[i] * fe_0 + tr_y_xxxx_xxzzz[i] * pa_x[i];

        tr_y_xxxxx_xyyyy[i] = 4.0 * tr_y_xxx_xyyyy[i] * fe_0 + tr_y_xxxx_yyyy[i] * fe_0 + tr_y_xxxx_xyyyy[i] * pa_x[i];

        tr_y_xxxxx_xyyyz[i] = 4.0 * tr_y_xxx_xyyyz[i] * fe_0 + tr_y_xxxx_yyyz[i] * fe_0 + tr_y_xxxx_xyyyz[i] * pa_x[i];

        tr_y_xxxxx_xyyzz[i] = 4.0 * tr_y_xxx_xyyzz[i] * fe_0 + tr_y_xxxx_yyzz[i] * fe_0 + tr_y_xxxx_xyyzz[i] * pa_x[i];

        tr_y_xxxxx_xyzzz[i] = 4.0 * tr_y_xxx_xyzzz[i] * fe_0 + tr_y_xxxx_yzzz[i] * fe_0 + tr_y_xxxx_xyzzz[i] * pa_x[i];

        tr_y_xxxxx_xzzzz[i] = 4.0 * tr_y_xxx_xzzzz[i] * fe_0 + tr_y_xxxx_zzzz[i] * fe_0 + tr_y_xxxx_xzzzz[i] * pa_x[i];

        tr_y_xxxxx_yyyyy[i] = 4.0 * tr_y_xxx_yyyyy[i] * fe_0 + tr_y_xxxx_yyyyy[i] * pa_x[i];

        tr_y_xxxxx_yyyyz[i] = 4.0 * tr_y_xxx_yyyyz[i] * fe_0 + tr_y_xxxx_yyyyz[i] * pa_x[i];

        tr_y_xxxxx_yyyzz[i] = 4.0 * tr_y_xxx_yyyzz[i] * fe_0 + tr_y_xxxx_yyyzz[i] * pa_x[i];

        tr_y_xxxxx_yyzzz[i] = 4.0 * tr_y_xxx_yyzzz[i] * fe_0 + tr_y_xxxx_yyzzz[i] * pa_x[i];

        tr_y_xxxxx_yzzzz[i] = 4.0 * tr_y_xxx_yzzzz[i] * fe_0 + tr_y_xxxx_yzzzz[i] * pa_x[i];

        tr_y_xxxxx_zzzzz[i] = 4.0 * tr_y_xxx_zzzzz[i] * fe_0 + tr_y_xxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 462-483 components of targeted buffer : HH

    auto tr_y_xxxxy_xxxxx = pbuffer.data(idx_dip_hh + 462);

    auto tr_y_xxxxy_xxxxy = pbuffer.data(idx_dip_hh + 463);

    auto tr_y_xxxxy_xxxxz = pbuffer.data(idx_dip_hh + 464);

    auto tr_y_xxxxy_xxxyy = pbuffer.data(idx_dip_hh + 465);

    auto tr_y_xxxxy_xxxyz = pbuffer.data(idx_dip_hh + 466);

    auto tr_y_xxxxy_xxxzz = pbuffer.data(idx_dip_hh + 467);

    auto tr_y_xxxxy_xxyyy = pbuffer.data(idx_dip_hh + 468);

    auto tr_y_xxxxy_xxyyz = pbuffer.data(idx_dip_hh + 469);

    auto tr_y_xxxxy_xxyzz = pbuffer.data(idx_dip_hh + 470);

    auto tr_y_xxxxy_xxzzz = pbuffer.data(idx_dip_hh + 471);

    auto tr_y_xxxxy_xyyyy = pbuffer.data(idx_dip_hh + 472);

    auto tr_y_xxxxy_xyyyz = pbuffer.data(idx_dip_hh + 473);

    auto tr_y_xxxxy_xyyzz = pbuffer.data(idx_dip_hh + 474);

    auto tr_y_xxxxy_xyzzz = pbuffer.data(idx_dip_hh + 475);

    auto tr_y_xxxxy_xzzzz = pbuffer.data(idx_dip_hh + 476);

    auto tr_y_xxxxy_yyyyy = pbuffer.data(idx_dip_hh + 477);

    auto tr_y_xxxxy_yyyyz = pbuffer.data(idx_dip_hh + 478);

    auto tr_y_xxxxy_yyyzz = pbuffer.data(idx_dip_hh + 479);

    auto tr_y_xxxxy_yyzzz = pbuffer.data(idx_dip_hh + 480);

    auto tr_y_xxxxy_yzzzz = pbuffer.data(idx_dip_hh + 481);

    auto tr_y_xxxxy_zzzzz = pbuffer.data(idx_dip_hh + 482);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_y_xxxx_xxxxx,  \
                             tr_y_xxxx_xxxxz,  \
                             tr_y_xxxx_xxxzz,  \
                             tr_y_xxxx_xxzzz,  \
                             tr_y_xxxx_xzzzz,  \
                             tr_y_xxxxy_xxxxx, \
                             tr_y_xxxxy_xxxxy, \
                             tr_y_xxxxy_xxxxz, \
                             tr_y_xxxxy_xxxyy, \
                             tr_y_xxxxy_xxxyz, \
                             tr_y_xxxxy_xxxzz, \
                             tr_y_xxxxy_xxyyy, \
                             tr_y_xxxxy_xxyyz, \
                             tr_y_xxxxy_xxyzz, \
                             tr_y_xxxxy_xxzzz, \
                             tr_y_xxxxy_xyyyy, \
                             tr_y_xxxxy_xyyyz, \
                             tr_y_xxxxy_xyyzz, \
                             tr_y_xxxxy_xyzzz, \
                             tr_y_xxxxy_xzzzz, \
                             tr_y_xxxxy_yyyyy, \
                             tr_y_xxxxy_yyyyz, \
                             tr_y_xxxxy_yyyzz, \
                             tr_y_xxxxy_yyzzz, \
                             tr_y_xxxxy_yzzzz, \
                             tr_y_xxxxy_zzzzz, \
                             tr_y_xxxy_xxxxy,  \
                             tr_y_xxxy_xxxy,   \
                             tr_y_xxxy_xxxyy,  \
                             tr_y_xxxy_xxxyz,  \
                             tr_y_xxxy_xxyy,   \
                             tr_y_xxxy_xxyyy,  \
                             tr_y_xxxy_xxyyz,  \
                             tr_y_xxxy_xxyz,   \
                             tr_y_xxxy_xxyzz,  \
                             tr_y_xxxy_xyyy,   \
                             tr_y_xxxy_xyyyy,  \
                             tr_y_xxxy_xyyyz,  \
                             tr_y_xxxy_xyyz,   \
                             tr_y_xxxy_xyyzz,  \
                             tr_y_xxxy_xyzz,   \
                             tr_y_xxxy_xyzzz,  \
                             tr_y_xxxy_yyyy,   \
                             tr_y_xxxy_yyyyy,  \
                             tr_y_xxxy_yyyyz,  \
                             tr_y_xxxy_yyyz,   \
                             tr_y_xxxy_yyyzz,  \
                             tr_y_xxxy_yyzz,   \
                             tr_y_xxxy_yyzzz,  \
                             tr_y_xxxy_yzzz,   \
                             tr_y_xxxy_yzzzz,  \
                             tr_y_xxxy_zzzzz,  \
                             tr_y_xxy_xxxxy,   \
                             tr_y_xxy_xxxyy,   \
                             tr_y_xxy_xxxyz,   \
                             tr_y_xxy_xxyyy,   \
                             tr_y_xxy_xxyyz,   \
                             tr_y_xxy_xxyzz,   \
                             tr_y_xxy_xyyyy,   \
                             tr_y_xxy_xyyyz,   \
                             tr_y_xxy_xyyzz,   \
                             tr_y_xxy_xyzzz,   \
                             tr_y_xxy_yyyyy,   \
                             tr_y_xxy_yyyyz,   \
                             tr_y_xxy_yyyzz,   \
                             tr_y_xxy_yyzzz,   \
                             tr_y_xxy_yzzzz,   \
                             tr_y_xxy_zzzzz,   \
                             ts_xxxx_xxxxx,    \
                             ts_xxxx_xxxxz,    \
                             ts_xxxx_xxxzz,    \
                             ts_xxxx_xxzzz,    \
                             ts_xxxx_xzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxy_xxxxx[i] = ts_xxxx_xxxxx[i] * fe_0 + tr_y_xxxx_xxxxx[i] * pa_y[i];

        tr_y_xxxxy_xxxxy[i] = 3.0 * tr_y_xxy_xxxxy[i] * fe_0 + 4.0 * tr_y_xxxy_xxxy[i] * fe_0 + tr_y_xxxy_xxxxy[i] * pa_x[i];

        tr_y_xxxxy_xxxxz[i] = ts_xxxx_xxxxz[i] * fe_0 + tr_y_xxxx_xxxxz[i] * pa_y[i];

        tr_y_xxxxy_xxxyy[i] = 3.0 * tr_y_xxy_xxxyy[i] * fe_0 + 3.0 * tr_y_xxxy_xxyy[i] * fe_0 + tr_y_xxxy_xxxyy[i] * pa_x[i];

        tr_y_xxxxy_xxxyz[i] = 3.0 * tr_y_xxy_xxxyz[i] * fe_0 + 3.0 * tr_y_xxxy_xxyz[i] * fe_0 + tr_y_xxxy_xxxyz[i] * pa_x[i];

        tr_y_xxxxy_xxxzz[i] = ts_xxxx_xxxzz[i] * fe_0 + tr_y_xxxx_xxxzz[i] * pa_y[i];

        tr_y_xxxxy_xxyyy[i] = 3.0 * tr_y_xxy_xxyyy[i] * fe_0 + 2.0 * tr_y_xxxy_xyyy[i] * fe_0 + tr_y_xxxy_xxyyy[i] * pa_x[i];

        tr_y_xxxxy_xxyyz[i] = 3.0 * tr_y_xxy_xxyyz[i] * fe_0 + 2.0 * tr_y_xxxy_xyyz[i] * fe_0 + tr_y_xxxy_xxyyz[i] * pa_x[i];

        tr_y_xxxxy_xxyzz[i] = 3.0 * tr_y_xxy_xxyzz[i] * fe_0 + 2.0 * tr_y_xxxy_xyzz[i] * fe_0 + tr_y_xxxy_xxyzz[i] * pa_x[i];

        tr_y_xxxxy_xxzzz[i] = ts_xxxx_xxzzz[i] * fe_0 + tr_y_xxxx_xxzzz[i] * pa_y[i];

        tr_y_xxxxy_xyyyy[i] = 3.0 * tr_y_xxy_xyyyy[i] * fe_0 + tr_y_xxxy_yyyy[i] * fe_0 + tr_y_xxxy_xyyyy[i] * pa_x[i];

        tr_y_xxxxy_xyyyz[i] = 3.0 * tr_y_xxy_xyyyz[i] * fe_0 + tr_y_xxxy_yyyz[i] * fe_0 + tr_y_xxxy_xyyyz[i] * pa_x[i];

        tr_y_xxxxy_xyyzz[i] = 3.0 * tr_y_xxy_xyyzz[i] * fe_0 + tr_y_xxxy_yyzz[i] * fe_0 + tr_y_xxxy_xyyzz[i] * pa_x[i];

        tr_y_xxxxy_xyzzz[i] = 3.0 * tr_y_xxy_xyzzz[i] * fe_0 + tr_y_xxxy_yzzz[i] * fe_0 + tr_y_xxxy_xyzzz[i] * pa_x[i];

        tr_y_xxxxy_xzzzz[i] = ts_xxxx_xzzzz[i] * fe_0 + tr_y_xxxx_xzzzz[i] * pa_y[i];

        tr_y_xxxxy_yyyyy[i] = 3.0 * tr_y_xxy_yyyyy[i] * fe_0 + tr_y_xxxy_yyyyy[i] * pa_x[i];

        tr_y_xxxxy_yyyyz[i] = 3.0 * tr_y_xxy_yyyyz[i] * fe_0 + tr_y_xxxy_yyyyz[i] * pa_x[i];

        tr_y_xxxxy_yyyzz[i] = 3.0 * tr_y_xxy_yyyzz[i] * fe_0 + tr_y_xxxy_yyyzz[i] * pa_x[i];

        tr_y_xxxxy_yyzzz[i] = 3.0 * tr_y_xxy_yyzzz[i] * fe_0 + tr_y_xxxy_yyzzz[i] * pa_x[i];

        tr_y_xxxxy_yzzzz[i] = 3.0 * tr_y_xxy_yzzzz[i] * fe_0 + tr_y_xxxy_yzzzz[i] * pa_x[i];

        tr_y_xxxxy_zzzzz[i] = 3.0 * tr_y_xxy_zzzzz[i] * fe_0 + tr_y_xxxy_zzzzz[i] * pa_x[i];
    }

    // Set up 483-504 components of targeted buffer : HH

    auto tr_y_xxxxz_xxxxx = pbuffer.data(idx_dip_hh + 483);

    auto tr_y_xxxxz_xxxxy = pbuffer.data(idx_dip_hh + 484);

    auto tr_y_xxxxz_xxxxz = pbuffer.data(idx_dip_hh + 485);

    auto tr_y_xxxxz_xxxyy = pbuffer.data(idx_dip_hh + 486);

    auto tr_y_xxxxz_xxxyz = pbuffer.data(idx_dip_hh + 487);

    auto tr_y_xxxxz_xxxzz = pbuffer.data(idx_dip_hh + 488);

    auto tr_y_xxxxz_xxyyy = pbuffer.data(idx_dip_hh + 489);

    auto tr_y_xxxxz_xxyyz = pbuffer.data(idx_dip_hh + 490);

    auto tr_y_xxxxz_xxyzz = pbuffer.data(idx_dip_hh + 491);

    auto tr_y_xxxxz_xxzzz = pbuffer.data(idx_dip_hh + 492);

    auto tr_y_xxxxz_xyyyy = pbuffer.data(idx_dip_hh + 493);

    auto tr_y_xxxxz_xyyyz = pbuffer.data(idx_dip_hh + 494);

    auto tr_y_xxxxz_xyyzz = pbuffer.data(idx_dip_hh + 495);

    auto tr_y_xxxxz_xyzzz = pbuffer.data(idx_dip_hh + 496);

    auto tr_y_xxxxz_xzzzz = pbuffer.data(idx_dip_hh + 497);

    auto tr_y_xxxxz_yyyyy = pbuffer.data(idx_dip_hh + 498);

    auto tr_y_xxxxz_yyyyz = pbuffer.data(idx_dip_hh + 499);

    auto tr_y_xxxxz_yyyzz = pbuffer.data(idx_dip_hh + 500);

    auto tr_y_xxxxz_yyzzz = pbuffer.data(idx_dip_hh + 501);

    auto tr_y_xxxxz_yzzzz = pbuffer.data(idx_dip_hh + 502);

    auto tr_y_xxxxz_zzzzz = pbuffer.data(idx_dip_hh + 503);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxxx_xxxx,   \
                             tr_y_xxxx_xxxxx,  \
                             tr_y_xxxx_xxxxy,  \
                             tr_y_xxxx_xxxxz,  \
                             tr_y_xxxx_xxxy,   \
                             tr_y_xxxx_xxxyy,  \
                             tr_y_xxxx_xxxyz,  \
                             tr_y_xxxx_xxxz,   \
                             tr_y_xxxx_xxxzz,  \
                             tr_y_xxxx_xxyy,   \
                             tr_y_xxxx_xxyyy,  \
                             tr_y_xxxx_xxyyz,  \
                             tr_y_xxxx_xxyz,   \
                             tr_y_xxxx_xxyzz,  \
                             tr_y_xxxx_xxzz,   \
                             tr_y_xxxx_xxzzz,  \
                             tr_y_xxxx_xyyy,   \
                             tr_y_xxxx_xyyyy,  \
                             tr_y_xxxx_xyyyz,  \
                             tr_y_xxxx_xyyz,   \
                             tr_y_xxxx_xyyzz,  \
                             tr_y_xxxx_xyzz,   \
                             tr_y_xxxx_xyzzz,  \
                             tr_y_xxxx_xzzz,   \
                             tr_y_xxxx_xzzzz,  \
                             tr_y_xxxx_yyyyy,  \
                             tr_y_xxxxz_xxxxx, \
                             tr_y_xxxxz_xxxxy, \
                             tr_y_xxxxz_xxxxz, \
                             tr_y_xxxxz_xxxyy, \
                             tr_y_xxxxz_xxxyz, \
                             tr_y_xxxxz_xxxzz, \
                             tr_y_xxxxz_xxyyy, \
                             tr_y_xxxxz_xxyyz, \
                             tr_y_xxxxz_xxyzz, \
                             tr_y_xxxxz_xxzzz, \
                             tr_y_xxxxz_xyyyy, \
                             tr_y_xxxxz_xyyyz, \
                             tr_y_xxxxz_xyyzz, \
                             tr_y_xxxxz_xyzzz, \
                             tr_y_xxxxz_xzzzz, \
                             tr_y_xxxxz_yyyyy, \
                             tr_y_xxxxz_yyyyz, \
                             tr_y_xxxxz_yyyzz, \
                             tr_y_xxxxz_yyzzz, \
                             tr_y_xxxxz_yzzzz, \
                             tr_y_xxxxz_zzzzz, \
                             tr_y_xxxz_yyyyz,  \
                             tr_y_xxxz_yyyzz,  \
                             tr_y_xxxz_yyzzz,  \
                             tr_y_xxxz_yzzzz,  \
                             tr_y_xxxz_zzzzz,  \
                             tr_y_xxz_yyyyz,   \
                             tr_y_xxz_yyyzz,   \
                             tr_y_xxz_yyzzz,   \
                             tr_y_xxz_yzzzz,   \
                             tr_y_xxz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxz_xxxxx[i] = tr_y_xxxx_xxxxx[i] * pa_z[i];

        tr_y_xxxxz_xxxxy[i] = tr_y_xxxx_xxxxy[i] * pa_z[i];

        tr_y_xxxxz_xxxxz[i] = tr_y_xxxx_xxxx[i] * fe_0 + tr_y_xxxx_xxxxz[i] * pa_z[i];

        tr_y_xxxxz_xxxyy[i] = tr_y_xxxx_xxxyy[i] * pa_z[i];

        tr_y_xxxxz_xxxyz[i] = tr_y_xxxx_xxxy[i] * fe_0 + tr_y_xxxx_xxxyz[i] * pa_z[i];

        tr_y_xxxxz_xxxzz[i] = 2.0 * tr_y_xxxx_xxxz[i] * fe_0 + tr_y_xxxx_xxxzz[i] * pa_z[i];

        tr_y_xxxxz_xxyyy[i] = tr_y_xxxx_xxyyy[i] * pa_z[i];

        tr_y_xxxxz_xxyyz[i] = tr_y_xxxx_xxyy[i] * fe_0 + tr_y_xxxx_xxyyz[i] * pa_z[i];

        tr_y_xxxxz_xxyzz[i] = 2.0 * tr_y_xxxx_xxyz[i] * fe_0 + tr_y_xxxx_xxyzz[i] * pa_z[i];

        tr_y_xxxxz_xxzzz[i] = 3.0 * tr_y_xxxx_xxzz[i] * fe_0 + tr_y_xxxx_xxzzz[i] * pa_z[i];

        tr_y_xxxxz_xyyyy[i] = tr_y_xxxx_xyyyy[i] * pa_z[i];

        tr_y_xxxxz_xyyyz[i] = tr_y_xxxx_xyyy[i] * fe_0 + tr_y_xxxx_xyyyz[i] * pa_z[i];

        tr_y_xxxxz_xyyzz[i] = 2.0 * tr_y_xxxx_xyyz[i] * fe_0 + tr_y_xxxx_xyyzz[i] * pa_z[i];

        tr_y_xxxxz_xyzzz[i] = 3.0 * tr_y_xxxx_xyzz[i] * fe_0 + tr_y_xxxx_xyzzz[i] * pa_z[i];

        tr_y_xxxxz_xzzzz[i] = 4.0 * tr_y_xxxx_xzzz[i] * fe_0 + tr_y_xxxx_xzzzz[i] * pa_z[i];

        tr_y_xxxxz_yyyyy[i] = tr_y_xxxx_yyyyy[i] * pa_z[i];

        tr_y_xxxxz_yyyyz[i] = 3.0 * tr_y_xxz_yyyyz[i] * fe_0 + tr_y_xxxz_yyyyz[i] * pa_x[i];

        tr_y_xxxxz_yyyzz[i] = 3.0 * tr_y_xxz_yyyzz[i] * fe_0 + tr_y_xxxz_yyyzz[i] * pa_x[i];

        tr_y_xxxxz_yyzzz[i] = 3.0 * tr_y_xxz_yyzzz[i] * fe_0 + tr_y_xxxz_yyzzz[i] * pa_x[i];

        tr_y_xxxxz_yzzzz[i] = 3.0 * tr_y_xxz_yzzzz[i] * fe_0 + tr_y_xxxz_yzzzz[i] * pa_x[i];

        tr_y_xxxxz_zzzzz[i] = 3.0 * tr_y_xxz_zzzzz[i] * fe_0 + tr_y_xxxz_zzzzz[i] * pa_x[i];
    }

    // Set up 504-525 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxxyy_xxxxx, \
                             tr_y_xxxyy_xxxxy, \
                             tr_y_xxxyy_xxxxz, \
                             tr_y_xxxyy_xxxyy, \
                             tr_y_xxxyy_xxxyz, \
                             tr_y_xxxyy_xxxzz, \
                             tr_y_xxxyy_xxyyy, \
                             tr_y_xxxyy_xxyyz, \
                             tr_y_xxxyy_xxyzz, \
                             tr_y_xxxyy_xxzzz, \
                             tr_y_xxxyy_xyyyy, \
                             tr_y_xxxyy_xyyyz, \
                             tr_y_xxxyy_xyyzz, \
                             tr_y_xxxyy_xyzzz, \
                             tr_y_xxxyy_xzzzz, \
                             tr_y_xxxyy_yyyyy, \
                             tr_y_xxxyy_yyyyz, \
                             tr_y_xxxyy_yyyzz, \
                             tr_y_xxxyy_yyzzz, \
                             tr_y_xxxyy_yzzzz, \
                             tr_y_xxxyy_zzzzz, \
                             tr_y_xxyy_xxxx,   \
                             tr_y_xxyy_xxxxx,  \
                             tr_y_xxyy_xxxxy,  \
                             tr_y_xxyy_xxxxz,  \
                             tr_y_xxyy_xxxy,   \
                             tr_y_xxyy_xxxyy,  \
                             tr_y_xxyy_xxxyz,  \
                             tr_y_xxyy_xxxz,   \
                             tr_y_xxyy_xxxzz,  \
                             tr_y_xxyy_xxyy,   \
                             tr_y_xxyy_xxyyy,  \
                             tr_y_xxyy_xxyyz,  \
                             tr_y_xxyy_xxyz,   \
                             tr_y_xxyy_xxyzz,  \
                             tr_y_xxyy_xxzz,   \
                             tr_y_xxyy_xxzzz,  \
                             tr_y_xxyy_xyyy,   \
                             tr_y_xxyy_xyyyy,  \
                             tr_y_xxyy_xyyyz,  \
                             tr_y_xxyy_xyyz,   \
                             tr_y_xxyy_xyyzz,  \
                             tr_y_xxyy_xyzz,   \
                             tr_y_xxyy_xyzzz,  \
                             tr_y_xxyy_xzzz,   \
                             tr_y_xxyy_xzzzz,  \
                             tr_y_xxyy_yyyy,   \
                             tr_y_xxyy_yyyyy,  \
                             tr_y_xxyy_yyyyz,  \
                             tr_y_xxyy_yyyz,   \
                             tr_y_xxyy_yyyzz,  \
                             tr_y_xxyy_yyzz,   \
                             tr_y_xxyy_yyzzz,  \
                             tr_y_xxyy_yzzz,   \
                             tr_y_xxyy_yzzzz,  \
                             tr_y_xxyy_zzzz,   \
                             tr_y_xxyy_zzzzz,  \
                             tr_y_xyy_xxxxx,   \
                             tr_y_xyy_xxxxy,   \
                             tr_y_xyy_xxxxz,   \
                             tr_y_xyy_xxxyy,   \
                             tr_y_xyy_xxxyz,   \
                             tr_y_xyy_xxxzz,   \
                             tr_y_xyy_xxyyy,   \
                             tr_y_xyy_xxyyz,   \
                             tr_y_xyy_xxyzz,   \
                             tr_y_xyy_xxzzz,   \
                             tr_y_xyy_xyyyy,   \
                             tr_y_xyy_xyyyz,   \
                             tr_y_xyy_xyyzz,   \
                             tr_y_xyy_xyzzz,   \
                             tr_y_xyy_xzzzz,   \
                             tr_y_xyy_yyyyy,   \
                             tr_y_xyy_yyyyz,   \
                             tr_y_xyy_yyyzz,   \
                             tr_y_xyy_yyzzz,   \
                             tr_y_xyy_yzzzz,   \
                             tr_y_xyy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyy_xxxxx[i] = 2.0 * tr_y_xyy_xxxxx[i] * fe_0 + 5.0 * tr_y_xxyy_xxxx[i] * fe_0 + tr_y_xxyy_xxxxx[i] * pa_x[i];

        tr_y_xxxyy_xxxxy[i] = 2.0 * tr_y_xyy_xxxxy[i] * fe_0 + 4.0 * tr_y_xxyy_xxxy[i] * fe_0 + tr_y_xxyy_xxxxy[i] * pa_x[i];

        tr_y_xxxyy_xxxxz[i] = 2.0 * tr_y_xyy_xxxxz[i] * fe_0 + 4.0 * tr_y_xxyy_xxxz[i] * fe_0 + tr_y_xxyy_xxxxz[i] * pa_x[i];

        tr_y_xxxyy_xxxyy[i] = 2.0 * tr_y_xyy_xxxyy[i] * fe_0 + 3.0 * tr_y_xxyy_xxyy[i] * fe_0 + tr_y_xxyy_xxxyy[i] * pa_x[i];

        tr_y_xxxyy_xxxyz[i] = 2.0 * tr_y_xyy_xxxyz[i] * fe_0 + 3.0 * tr_y_xxyy_xxyz[i] * fe_0 + tr_y_xxyy_xxxyz[i] * pa_x[i];

        tr_y_xxxyy_xxxzz[i] = 2.0 * tr_y_xyy_xxxzz[i] * fe_0 + 3.0 * tr_y_xxyy_xxzz[i] * fe_0 + tr_y_xxyy_xxxzz[i] * pa_x[i];

        tr_y_xxxyy_xxyyy[i] = 2.0 * tr_y_xyy_xxyyy[i] * fe_0 + 2.0 * tr_y_xxyy_xyyy[i] * fe_0 + tr_y_xxyy_xxyyy[i] * pa_x[i];

        tr_y_xxxyy_xxyyz[i] = 2.0 * tr_y_xyy_xxyyz[i] * fe_0 + 2.0 * tr_y_xxyy_xyyz[i] * fe_0 + tr_y_xxyy_xxyyz[i] * pa_x[i];

        tr_y_xxxyy_xxyzz[i] = 2.0 * tr_y_xyy_xxyzz[i] * fe_0 + 2.0 * tr_y_xxyy_xyzz[i] * fe_0 + tr_y_xxyy_xxyzz[i] * pa_x[i];

        tr_y_xxxyy_xxzzz[i] = 2.0 * tr_y_xyy_xxzzz[i] * fe_0 + 2.0 * tr_y_xxyy_xzzz[i] * fe_0 + tr_y_xxyy_xxzzz[i] * pa_x[i];

        tr_y_xxxyy_xyyyy[i] = 2.0 * tr_y_xyy_xyyyy[i] * fe_0 + tr_y_xxyy_yyyy[i] * fe_0 + tr_y_xxyy_xyyyy[i] * pa_x[i];

        tr_y_xxxyy_xyyyz[i] = 2.0 * tr_y_xyy_xyyyz[i] * fe_0 + tr_y_xxyy_yyyz[i] * fe_0 + tr_y_xxyy_xyyyz[i] * pa_x[i];

        tr_y_xxxyy_xyyzz[i] = 2.0 * tr_y_xyy_xyyzz[i] * fe_0 + tr_y_xxyy_yyzz[i] * fe_0 + tr_y_xxyy_xyyzz[i] * pa_x[i];

        tr_y_xxxyy_xyzzz[i] = 2.0 * tr_y_xyy_xyzzz[i] * fe_0 + tr_y_xxyy_yzzz[i] * fe_0 + tr_y_xxyy_xyzzz[i] * pa_x[i];

        tr_y_xxxyy_xzzzz[i] = 2.0 * tr_y_xyy_xzzzz[i] * fe_0 + tr_y_xxyy_zzzz[i] * fe_0 + tr_y_xxyy_xzzzz[i] * pa_x[i];

        tr_y_xxxyy_yyyyy[i] = 2.0 * tr_y_xyy_yyyyy[i] * fe_0 + tr_y_xxyy_yyyyy[i] * pa_x[i];

        tr_y_xxxyy_yyyyz[i] = 2.0 * tr_y_xyy_yyyyz[i] * fe_0 + tr_y_xxyy_yyyyz[i] * pa_x[i];

        tr_y_xxxyy_yyyzz[i] = 2.0 * tr_y_xyy_yyyzz[i] * fe_0 + tr_y_xxyy_yyyzz[i] * pa_x[i];

        tr_y_xxxyy_yyzzz[i] = 2.0 * tr_y_xyy_yyzzz[i] * fe_0 + tr_y_xxyy_yyzzz[i] * pa_x[i];

        tr_y_xxxyy_yzzzz[i] = 2.0 * tr_y_xyy_yzzzz[i] * fe_0 + tr_y_xxyy_yzzzz[i] * pa_x[i];

        tr_y_xxxyy_zzzzz[i] = 2.0 * tr_y_xyy_zzzzz[i] * fe_0 + tr_y_xxyy_zzzzz[i] * pa_x[i];
    }

    // Set up 525-546 components of targeted buffer : HH

    auto tr_y_xxxyz_xxxxx = pbuffer.data(idx_dip_hh + 525);

    auto tr_y_xxxyz_xxxxy = pbuffer.data(idx_dip_hh + 526);

    auto tr_y_xxxyz_xxxxz = pbuffer.data(idx_dip_hh + 527);

    auto tr_y_xxxyz_xxxyy = pbuffer.data(idx_dip_hh + 528);

    auto tr_y_xxxyz_xxxyz = pbuffer.data(idx_dip_hh + 529);

    auto tr_y_xxxyz_xxxzz = pbuffer.data(idx_dip_hh + 530);

    auto tr_y_xxxyz_xxyyy = pbuffer.data(idx_dip_hh + 531);

    auto tr_y_xxxyz_xxyyz = pbuffer.data(idx_dip_hh + 532);

    auto tr_y_xxxyz_xxyzz = pbuffer.data(idx_dip_hh + 533);

    auto tr_y_xxxyz_xxzzz = pbuffer.data(idx_dip_hh + 534);

    auto tr_y_xxxyz_xyyyy = pbuffer.data(idx_dip_hh + 535);

    auto tr_y_xxxyz_xyyyz = pbuffer.data(idx_dip_hh + 536);

    auto tr_y_xxxyz_xyyzz = pbuffer.data(idx_dip_hh + 537);

    auto tr_y_xxxyz_xyzzz = pbuffer.data(idx_dip_hh + 538);

    auto tr_y_xxxyz_xzzzz = pbuffer.data(idx_dip_hh + 539);

    auto tr_y_xxxyz_yyyyy = pbuffer.data(idx_dip_hh + 540);

    auto tr_y_xxxyz_yyyyz = pbuffer.data(idx_dip_hh + 541);

    auto tr_y_xxxyz_yyyzz = pbuffer.data(idx_dip_hh + 542);

    auto tr_y_xxxyz_yyzzz = pbuffer.data(idx_dip_hh + 543);

    auto tr_y_xxxyz_yzzzz = pbuffer.data(idx_dip_hh + 544);

    auto tr_y_xxxyz_zzzzz = pbuffer.data(idx_dip_hh + 545);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_y_xxxy_xxxxx,  \
                             tr_y_xxxy_xxxxy,  \
                             tr_y_xxxy_xxxy,   \
                             tr_y_xxxy_xxxyy,  \
                             tr_y_xxxy_xxxyz,  \
                             tr_y_xxxy_xxyy,   \
                             tr_y_xxxy_xxyyy,  \
                             tr_y_xxxy_xxyyz,  \
                             tr_y_xxxy_xxyz,   \
                             tr_y_xxxy_xxyzz,  \
                             tr_y_xxxy_xyyy,   \
                             tr_y_xxxy_xyyyy,  \
                             tr_y_xxxy_xyyyz,  \
                             tr_y_xxxy_xyyz,   \
                             tr_y_xxxy_xyyzz,  \
                             tr_y_xxxy_xyzz,   \
                             tr_y_xxxy_xyzzz,  \
                             tr_y_xxxy_yyyyy,  \
                             tr_y_xxxyz_xxxxx, \
                             tr_y_xxxyz_xxxxy, \
                             tr_y_xxxyz_xxxxz, \
                             tr_y_xxxyz_xxxyy, \
                             tr_y_xxxyz_xxxyz, \
                             tr_y_xxxyz_xxxzz, \
                             tr_y_xxxyz_xxyyy, \
                             tr_y_xxxyz_xxyyz, \
                             tr_y_xxxyz_xxyzz, \
                             tr_y_xxxyz_xxzzz, \
                             tr_y_xxxyz_xyyyy, \
                             tr_y_xxxyz_xyyyz, \
                             tr_y_xxxyz_xyyzz, \
                             tr_y_xxxyz_xyzzz, \
                             tr_y_xxxyz_xzzzz, \
                             tr_y_xxxyz_yyyyy, \
                             tr_y_xxxyz_yyyyz, \
                             tr_y_xxxyz_yyyzz, \
                             tr_y_xxxyz_yyzzz, \
                             tr_y_xxxyz_yzzzz, \
                             tr_y_xxxyz_zzzzz, \
                             tr_y_xxxz_xxxxz,  \
                             tr_y_xxxz_xxxzz,  \
                             tr_y_xxxz_xxzzz,  \
                             tr_y_xxxz_xzzzz,  \
                             tr_y_xxyz_yyyyz,  \
                             tr_y_xxyz_yyyzz,  \
                             tr_y_xxyz_yyzzz,  \
                             tr_y_xxyz_yzzzz,  \
                             tr_y_xxyz_zzzzz,  \
                             tr_y_xyz_yyyyz,   \
                             tr_y_xyz_yyyzz,   \
                             tr_y_xyz_yyzzz,   \
                             tr_y_xyz_yzzzz,   \
                             tr_y_xyz_zzzzz,   \
                             ts_xxxz_xxxxz,    \
                             ts_xxxz_xxxzz,    \
                             ts_xxxz_xxzzz,    \
                             ts_xxxz_xzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyz_xxxxx[i] = tr_y_xxxy_xxxxx[i] * pa_z[i];

        tr_y_xxxyz_xxxxy[i] = tr_y_xxxy_xxxxy[i] * pa_z[i];

        tr_y_xxxyz_xxxxz[i] = ts_xxxz_xxxxz[i] * fe_0 + tr_y_xxxz_xxxxz[i] * pa_y[i];

        tr_y_xxxyz_xxxyy[i] = tr_y_xxxy_xxxyy[i] * pa_z[i];

        tr_y_xxxyz_xxxyz[i] = tr_y_xxxy_xxxy[i] * fe_0 + tr_y_xxxy_xxxyz[i] * pa_z[i];

        tr_y_xxxyz_xxxzz[i] = ts_xxxz_xxxzz[i] * fe_0 + tr_y_xxxz_xxxzz[i] * pa_y[i];

        tr_y_xxxyz_xxyyy[i] = tr_y_xxxy_xxyyy[i] * pa_z[i];

        tr_y_xxxyz_xxyyz[i] = tr_y_xxxy_xxyy[i] * fe_0 + tr_y_xxxy_xxyyz[i] * pa_z[i];

        tr_y_xxxyz_xxyzz[i] = 2.0 * tr_y_xxxy_xxyz[i] * fe_0 + tr_y_xxxy_xxyzz[i] * pa_z[i];

        tr_y_xxxyz_xxzzz[i] = ts_xxxz_xxzzz[i] * fe_0 + tr_y_xxxz_xxzzz[i] * pa_y[i];

        tr_y_xxxyz_xyyyy[i] = tr_y_xxxy_xyyyy[i] * pa_z[i];

        tr_y_xxxyz_xyyyz[i] = tr_y_xxxy_xyyy[i] * fe_0 + tr_y_xxxy_xyyyz[i] * pa_z[i];

        tr_y_xxxyz_xyyzz[i] = 2.0 * tr_y_xxxy_xyyz[i] * fe_0 + tr_y_xxxy_xyyzz[i] * pa_z[i];

        tr_y_xxxyz_xyzzz[i] = 3.0 * tr_y_xxxy_xyzz[i] * fe_0 + tr_y_xxxy_xyzzz[i] * pa_z[i];

        tr_y_xxxyz_xzzzz[i] = ts_xxxz_xzzzz[i] * fe_0 + tr_y_xxxz_xzzzz[i] * pa_y[i];

        tr_y_xxxyz_yyyyy[i] = tr_y_xxxy_yyyyy[i] * pa_z[i];

        tr_y_xxxyz_yyyyz[i] = 2.0 * tr_y_xyz_yyyyz[i] * fe_0 + tr_y_xxyz_yyyyz[i] * pa_x[i];

        tr_y_xxxyz_yyyzz[i] = 2.0 * tr_y_xyz_yyyzz[i] * fe_0 + tr_y_xxyz_yyyzz[i] * pa_x[i];

        tr_y_xxxyz_yyzzz[i] = 2.0 * tr_y_xyz_yyzzz[i] * fe_0 + tr_y_xxyz_yyzzz[i] * pa_x[i];

        tr_y_xxxyz_yzzzz[i] = 2.0 * tr_y_xyz_yzzzz[i] * fe_0 + tr_y_xxyz_yzzzz[i] * pa_x[i];

        tr_y_xxxyz_zzzzz[i] = 2.0 * tr_y_xyz_zzzzz[i] * fe_0 + tr_y_xxyz_zzzzz[i] * pa_x[i];
    }

    // Set up 546-567 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxx_xxxxx,   \
                             tr_y_xxx_xxxxy,   \
                             tr_y_xxx_xxxyy,   \
                             tr_y_xxx_xxyyy,   \
                             tr_y_xxx_xyyyy,   \
                             tr_y_xxxz_xxxxx,  \
                             tr_y_xxxz_xxxxy,  \
                             tr_y_xxxz_xxxyy,  \
                             tr_y_xxxz_xxyyy,  \
                             tr_y_xxxz_xyyyy,  \
                             tr_y_xxxzz_xxxxx, \
                             tr_y_xxxzz_xxxxy, \
                             tr_y_xxxzz_xxxxz, \
                             tr_y_xxxzz_xxxyy, \
                             tr_y_xxxzz_xxxyz, \
                             tr_y_xxxzz_xxxzz, \
                             tr_y_xxxzz_xxyyy, \
                             tr_y_xxxzz_xxyyz, \
                             tr_y_xxxzz_xxyzz, \
                             tr_y_xxxzz_xxzzz, \
                             tr_y_xxxzz_xyyyy, \
                             tr_y_xxxzz_xyyyz, \
                             tr_y_xxxzz_xyyzz, \
                             tr_y_xxxzz_xyzzz, \
                             tr_y_xxxzz_xzzzz, \
                             tr_y_xxxzz_yyyyy, \
                             tr_y_xxxzz_yyyyz, \
                             tr_y_xxxzz_yyyzz, \
                             tr_y_xxxzz_yyzzz, \
                             tr_y_xxxzz_yzzzz, \
                             tr_y_xxxzz_zzzzz, \
                             tr_y_xxzz_xxxxz,  \
                             tr_y_xxzz_xxxyz,  \
                             tr_y_xxzz_xxxz,   \
                             tr_y_xxzz_xxxzz,  \
                             tr_y_xxzz_xxyyz,  \
                             tr_y_xxzz_xxyz,   \
                             tr_y_xxzz_xxyzz,  \
                             tr_y_xxzz_xxzz,   \
                             tr_y_xxzz_xxzzz,  \
                             tr_y_xxzz_xyyyz,  \
                             tr_y_xxzz_xyyz,   \
                             tr_y_xxzz_xyyzz,  \
                             tr_y_xxzz_xyzz,   \
                             tr_y_xxzz_xyzzz,  \
                             tr_y_xxzz_xzzz,   \
                             tr_y_xxzz_xzzzz,  \
                             tr_y_xxzz_yyyyy,  \
                             tr_y_xxzz_yyyyz,  \
                             tr_y_xxzz_yyyz,   \
                             tr_y_xxzz_yyyzz,  \
                             tr_y_xxzz_yyzz,   \
                             tr_y_xxzz_yyzzz,  \
                             tr_y_xxzz_yzzz,   \
                             tr_y_xxzz_yzzzz,  \
                             tr_y_xxzz_zzzz,   \
                             tr_y_xxzz_zzzzz,  \
                             tr_y_xzz_xxxxz,   \
                             tr_y_xzz_xxxyz,   \
                             tr_y_xzz_xxxzz,   \
                             tr_y_xzz_xxyyz,   \
                             tr_y_xzz_xxyzz,   \
                             tr_y_xzz_xxzzz,   \
                             tr_y_xzz_xyyyz,   \
                             tr_y_xzz_xyyzz,   \
                             tr_y_xzz_xyzzz,   \
                             tr_y_xzz_xzzzz,   \
                             tr_y_xzz_yyyyy,   \
                             tr_y_xzz_yyyyz,   \
                             tr_y_xzz_yyyzz,   \
                             tr_y_xzz_yyzzz,   \
                             tr_y_xzz_yzzzz,   \
                             tr_y_xzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzz_xxxxx[i] = tr_y_xxx_xxxxx[i] * fe_0 + tr_y_xxxz_xxxxx[i] * pa_z[i];

        tr_y_xxxzz_xxxxy[i] = tr_y_xxx_xxxxy[i] * fe_0 + tr_y_xxxz_xxxxy[i] * pa_z[i];

        tr_y_xxxzz_xxxxz[i] = 2.0 * tr_y_xzz_xxxxz[i] * fe_0 + 4.0 * tr_y_xxzz_xxxz[i] * fe_0 + tr_y_xxzz_xxxxz[i] * pa_x[i];

        tr_y_xxxzz_xxxyy[i] = tr_y_xxx_xxxyy[i] * fe_0 + tr_y_xxxz_xxxyy[i] * pa_z[i];

        tr_y_xxxzz_xxxyz[i] = 2.0 * tr_y_xzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xxzz_xxyz[i] * fe_0 + tr_y_xxzz_xxxyz[i] * pa_x[i];

        tr_y_xxxzz_xxxzz[i] = 2.0 * tr_y_xzz_xxxzz[i] * fe_0 + 3.0 * tr_y_xxzz_xxzz[i] * fe_0 + tr_y_xxzz_xxxzz[i] * pa_x[i];

        tr_y_xxxzz_xxyyy[i] = tr_y_xxx_xxyyy[i] * fe_0 + tr_y_xxxz_xxyyy[i] * pa_z[i];

        tr_y_xxxzz_xxyyz[i] = 2.0 * tr_y_xzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xxzz_xyyz[i] * fe_0 + tr_y_xxzz_xxyyz[i] * pa_x[i];

        tr_y_xxxzz_xxyzz[i] = 2.0 * tr_y_xzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xxzz_xyzz[i] * fe_0 + tr_y_xxzz_xxyzz[i] * pa_x[i];

        tr_y_xxxzz_xxzzz[i] = 2.0 * tr_y_xzz_xxzzz[i] * fe_0 + 2.0 * tr_y_xxzz_xzzz[i] * fe_0 + tr_y_xxzz_xxzzz[i] * pa_x[i];

        tr_y_xxxzz_xyyyy[i] = tr_y_xxx_xyyyy[i] * fe_0 + tr_y_xxxz_xyyyy[i] * pa_z[i];

        tr_y_xxxzz_xyyyz[i] = 2.0 * tr_y_xzz_xyyyz[i] * fe_0 + tr_y_xxzz_yyyz[i] * fe_0 + tr_y_xxzz_xyyyz[i] * pa_x[i];

        tr_y_xxxzz_xyyzz[i] = 2.0 * tr_y_xzz_xyyzz[i] * fe_0 + tr_y_xxzz_yyzz[i] * fe_0 + tr_y_xxzz_xyyzz[i] * pa_x[i];

        tr_y_xxxzz_xyzzz[i] = 2.0 * tr_y_xzz_xyzzz[i] * fe_0 + tr_y_xxzz_yzzz[i] * fe_0 + tr_y_xxzz_xyzzz[i] * pa_x[i];

        tr_y_xxxzz_xzzzz[i] = 2.0 * tr_y_xzz_xzzzz[i] * fe_0 + tr_y_xxzz_zzzz[i] * fe_0 + tr_y_xxzz_xzzzz[i] * pa_x[i];

        tr_y_xxxzz_yyyyy[i] = 2.0 * tr_y_xzz_yyyyy[i] * fe_0 + tr_y_xxzz_yyyyy[i] * pa_x[i];

        tr_y_xxxzz_yyyyz[i] = 2.0 * tr_y_xzz_yyyyz[i] * fe_0 + tr_y_xxzz_yyyyz[i] * pa_x[i];

        tr_y_xxxzz_yyyzz[i] = 2.0 * tr_y_xzz_yyyzz[i] * fe_0 + tr_y_xxzz_yyyzz[i] * pa_x[i];

        tr_y_xxxzz_yyzzz[i] = 2.0 * tr_y_xzz_yyzzz[i] * fe_0 + tr_y_xxzz_yyzzz[i] * pa_x[i];

        tr_y_xxxzz_yzzzz[i] = 2.0 * tr_y_xzz_yzzzz[i] * fe_0 + tr_y_xxzz_yzzzz[i] * pa_x[i];

        tr_y_xxxzz_zzzzz[i] = 2.0 * tr_y_xzz_zzzzz[i] * fe_0 + tr_y_xxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 567-588 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xxyyy_xxxxx, \
                             tr_y_xxyyy_xxxxy, \
                             tr_y_xxyyy_xxxxz, \
                             tr_y_xxyyy_xxxyy, \
                             tr_y_xxyyy_xxxyz, \
                             tr_y_xxyyy_xxxzz, \
                             tr_y_xxyyy_xxyyy, \
                             tr_y_xxyyy_xxyyz, \
                             tr_y_xxyyy_xxyzz, \
                             tr_y_xxyyy_xxzzz, \
                             tr_y_xxyyy_xyyyy, \
                             tr_y_xxyyy_xyyyz, \
                             tr_y_xxyyy_xyyzz, \
                             tr_y_xxyyy_xyzzz, \
                             tr_y_xxyyy_xzzzz, \
                             tr_y_xxyyy_yyyyy, \
                             tr_y_xxyyy_yyyyz, \
                             tr_y_xxyyy_yyyzz, \
                             tr_y_xxyyy_yyzzz, \
                             tr_y_xxyyy_yzzzz, \
                             tr_y_xxyyy_zzzzz, \
                             tr_y_xyyy_xxxx,   \
                             tr_y_xyyy_xxxxx,  \
                             tr_y_xyyy_xxxxy,  \
                             tr_y_xyyy_xxxxz,  \
                             tr_y_xyyy_xxxy,   \
                             tr_y_xyyy_xxxyy,  \
                             tr_y_xyyy_xxxyz,  \
                             tr_y_xyyy_xxxz,   \
                             tr_y_xyyy_xxxzz,  \
                             tr_y_xyyy_xxyy,   \
                             tr_y_xyyy_xxyyy,  \
                             tr_y_xyyy_xxyyz,  \
                             tr_y_xyyy_xxyz,   \
                             tr_y_xyyy_xxyzz,  \
                             tr_y_xyyy_xxzz,   \
                             tr_y_xyyy_xxzzz,  \
                             tr_y_xyyy_xyyy,   \
                             tr_y_xyyy_xyyyy,  \
                             tr_y_xyyy_xyyyz,  \
                             tr_y_xyyy_xyyz,   \
                             tr_y_xyyy_xyyzz,  \
                             tr_y_xyyy_xyzz,   \
                             tr_y_xyyy_xyzzz,  \
                             tr_y_xyyy_xzzz,   \
                             tr_y_xyyy_xzzzz,  \
                             tr_y_xyyy_yyyy,   \
                             tr_y_xyyy_yyyyy,  \
                             tr_y_xyyy_yyyyz,  \
                             tr_y_xyyy_yyyz,   \
                             tr_y_xyyy_yyyzz,  \
                             tr_y_xyyy_yyzz,   \
                             tr_y_xyyy_yyzzz,  \
                             tr_y_xyyy_yzzz,   \
                             tr_y_xyyy_yzzzz,  \
                             tr_y_xyyy_zzzz,   \
                             tr_y_xyyy_zzzzz,  \
                             tr_y_yyy_xxxxx,   \
                             tr_y_yyy_xxxxy,   \
                             tr_y_yyy_xxxxz,   \
                             tr_y_yyy_xxxyy,   \
                             tr_y_yyy_xxxyz,   \
                             tr_y_yyy_xxxzz,   \
                             tr_y_yyy_xxyyy,   \
                             tr_y_yyy_xxyyz,   \
                             tr_y_yyy_xxyzz,   \
                             tr_y_yyy_xxzzz,   \
                             tr_y_yyy_xyyyy,   \
                             tr_y_yyy_xyyyz,   \
                             tr_y_yyy_xyyzz,   \
                             tr_y_yyy_xyzzz,   \
                             tr_y_yyy_xzzzz,   \
                             tr_y_yyy_yyyyy,   \
                             tr_y_yyy_yyyyz,   \
                             tr_y_yyy_yyyzz,   \
                             tr_y_yyy_yyzzz,   \
                             tr_y_yyy_yzzzz,   \
                             tr_y_yyy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyy_xxxxx[i] = tr_y_yyy_xxxxx[i] * fe_0 + 5.0 * tr_y_xyyy_xxxx[i] * fe_0 + tr_y_xyyy_xxxxx[i] * pa_x[i];

        tr_y_xxyyy_xxxxy[i] = tr_y_yyy_xxxxy[i] * fe_0 + 4.0 * tr_y_xyyy_xxxy[i] * fe_0 + tr_y_xyyy_xxxxy[i] * pa_x[i];

        tr_y_xxyyy_xxxxz[i] = tr_y_yyy_xxxxz[i] * fe_0 + 4.0 * tr_y_xyyy_xxxz[i] * fe_0 + tr_y_xyyy_xxxxz[i] * pa_x[i];

        tr_y_xxyyy_xxxyy[i] = tr_y_yyy_xxxyy[i] * fe_0 + 3.0 * tr_y_xyyy_xxyy[i] * fe_0 + tr_y_xyyy_xxxyy[i] * pa_x[i];

        tr_y_xxyyy_xxxyz[i] = tr_y_yyy_xxxyz[i] * fe_0 + 3.0 * tr_y_xyyy_xxyz[i] * fe_0 + tr_y_xyyy_xxxyz[i] * pa_x[i];

        tr_y_xxyyy_xxxzz[i] = tr_y_yyy_xxxzz[i] * fe_0 + 3.0 * tr_y_xyyy_xxzz[i] * fe_0 + tr_y_xyyy_xxxzz[i] * pa_x[i];

        tr_y_xxyyy_xxyyy[i] = tr_y_yyy_xxyyy[i] * fe_0 + 2.0 * tr_y_xyyy_xyyy[i] * fe_0 + tr_y_xyyy_xxyyy[i] * pa_x[i];

        tr_y_xxyyy_xxyyz[i] = tr_y_yyy_xxyyz[i] * fe_0 + 2.0 * tr_y_xyyy_xyyz[i] * fe_0 + tr_y_xyyy_xxyyz[i] * pa_x[i];

        tr_y_xxyyy_xxyzz[i] = tr_y_yyy_xxyzz[i] * fe_0 + 2.0 * tr_y_xyyy_xyzz[i] * fe_0 + tr_y_xyyy_xxyzz[i] * pa_x[i];

        tr_y_xxyyy_xxzzz[i] = tr_y_yyy_xxzzz[i] * fe_0 + 2.0 * tr_y_xyyy_xzzz[i] * fe_0 + tr_y_xyyy_xxzzz[i] * pa_x[i];

        tr_y_xxyyy_xyyyy[i] = tr_y_yyy_xyyyy[i] * fe_0 + tr_y_xyyy_yyyy[i] * fe_0 + tr_y_xyyy_xyyyy[i] * pa_x[i];

        tr_y_xxyyy_xyyyz[i] = tr_y_yyy_xyyyz[i] * fe_0 + tr_y_xyyy_yyyz[i] * fe_0 + tr_y_xyyy_xyyyz[i] * pa_x[i];

        tr_y_xxyyy_xyyzz[i] = tr_y_yyy_xyyzz[i] * fe_0 + tr_y_xyyy_yyzz[i] * fe_0 + tr_y_xyyy_xyyzz[i] * pa_x[i];

        tr_y_xxyyy_xyzzz[i] = tr_y_yyy_xyzzz[i] * fe_0 + tr_y_xyyy_yzzz[i] * fe_0 + tr_y_xyyy_xyzzz[i] * pa_x[i];

        tr_y_xxyyy_xzzzz[i] = tr_y_yyy_xzzzz[i] * fe_0 + tr_y_xyyy_zzzz[i] * fe_0 + tr_y_xyyy_xzzzz[i] * pa_x[i];

        tr_y_xxyyy_yyyyy[i] = tr_y_yyy_yyyyy[i] * fe_0 + tr_y_xyyy_yyyyy[i] * pa_x[i];

        tr_y_xxyyy_yyyyz[i] = tr_y_yyy_yyyyz[i] * fe_0 + tr_y_xyyy_yyyyz[i] * pa_x[i];

        tr_y_xxyyy_yyyzz[i] = tr_y_yyy_yyyzz[i] * fe_0 + tr_y_xyyy_yyyzz[i] * pa_x[i];

        tr_y_xxyyy_yyzzz[i] = tr_y_yyy_yyzzz[i] * fe_0 + tr_y_xyyy_yyzzz[i] * pa_x[i];

        tr_y_xxyyy_yzzzz[i] = tr_y_yyy_yzzzz[i] * fe_0 + tr_y_xyyy_yzzzz[i] * pa_x[i];

        tr_y_xxyyy_zzzzz[i] = tr_y_yyy_zzzzz[i] * fe_0 + tr_y_xyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 588-609 components of targeted buffer : HH

    auto tr_y_xxyyz_xxxxx = pbuffer.data(idx_dip_hh + 588);

    auto tr_y_xxyyz_xxxxy = pbuffer.data(idx_dip_hh + 589);

    auto tr_y_xxyyz_xxxxz = pbuffer.data(idx_dip_hh + 590);

    auto tr_y_xxyyz_xxxyy = pbuffer.data(idx_dip_hh + 591);

    auto tr_y_xxyyz_xxxyz = pbuffer.data(idx_dip_hh + 592);

    auto tr_y_xxyyz_xxxzz = pbuffer.data(idx_dip_hh + 593);

    auto tr_y_xxyyz_xxyyy = pbuffer.data(idx_dip_hh + 594);

    auto tr_y_xxyyz_xxyyz = pbuffer.data(idx_dip_hh + 595);

    auto tr_y_xxyyz_xxyzz = pbuffer.data(idx_dip_hh + 596);

    auto tr_y_xxyyz_xxzzz = pbuffer.data(idx_dip_hh + 597);

    auto tr_y_xxyyz_xyyyy = pbuffer.data(idx_dip_hh + 598);

    auto tr_y_xxyyz_xyyyz = pbuffer.data(idx_dip_hh + 599);

    auto tr_y_xxyyz_xyyzz = pbuffer.data(idx_dip_hh + 600);

    auto tr_y_xxyyz_xyzzz = pbuffer.data(idx_dip_hh + 601);

    auto tr_y_xxyyz_xzzzz = pbuffer.data(idx_dip_hh + 602);

    auto tr_y_xxyyz_yyyyy = pbuffer.data(idx_dip_hh + 603);

    auto tr_y_xxyyz_yyyyz = pbuffer.data(idx_dip_hh + 604);

    auto tr_y_xxyyz_yyyzz = pbuffer.data(idx_dip_hh + 605);

    auto tr_y_xxyyz_yyzzz = pbuffer.data(idx_dip_hh + 606);

    auto tr_y_xxyyz_yzzzz = pbuffer.data(idx_dip_hh + 607);

    auto tr_y_xxyyz_zzzzz = pbuffer.data(idx_dip_hh + 608);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxyy_xxxx,   \
                             tr_y_xxyy_xxxxx,  \
                             tr_y_xxyy_xxxxy,  \
                             tr_y_xxyy_xxxxz,  \
                             tr_y_xxyy_xxxy,   \
                             tr_y_xxyy_xxxyy,  \
                             tr_y_xxyy_xxxyz,  \
                             tr_y_xxyy_xxxz,   \
                             tr_y_xxyy_xxxzz,  \
                             tr_y_xxyy_xxyy,   \
                             tr_y_xxyy_xxyyy,  \
                             tr_y_xxyy_xxyyz,  \
                             tr_y_xxyy_xxyz,   \
                             tr_y_xxyy_xxyzz,  \
                             tr_y_xxyy_xxzz,   \
                             tr_y_xxyy_xxzzz,  \
                             tr_y_xxyy_xyyy,   \
                             tr_y_xxyy_xyyyy,  \
                             tr_y_xxyy_xyyyz,  \
                             tr_y_xxyy_xyyz,   \
                             tr_y_xxyy_xyyzz,  \
                             tr_y_xxyy_xyzz,   \
                             tr_y_xxyy_xyzzz,  \
                             tr_y_xxyy_xzzz,   \
                             tr_y_xxyy_xzzzz,  \
                             tr_y_xxyy_yyyyy,  \
                             tr_y_xxyyz_xxxxx, \
                             tr_y_xxyyz_xxxxy, \
                             tr_y_xxyyz_xxxxz, \
                             tr_y_xxyyz_xxxyy, \
                             tr_y_xxyyz_xxxyz, \
                             tr_y_xxyyz_xxxzz, \
                             tr_y_xxyyz_xxyyy, \
                             tr_y_xxyyz_xxyyz, \
                             tr_y_xxyyz_xxyzz, \
                             tr_y_xxyyz_xxzzz, \
                             tr_y_xxyyz_xyyyy, \
                             tr_y_xxyyz_xyyyz, \
                             tr_y_xxyyz_xyyzz, \
                             tr_y_xxyyz_xyzzz, \
                             tr_y_xxyyz_xzzzz, \
                             tr_y_xxyyz_yyyyy, \
                             tr_y_xxyyz_yyyyz, \
                             tr_y_xxyyz_yyyzz, \
                             tr_y_xxyyz_yyzzz, \
                             tr_y_xxyyz_yzzzz, \
                             tr_y_xxyyz_zzzzz, \
                             tr_y_xyyz_yyyyz,  \
                             tr_y_xyyz_yyyzz,  \
                             tr_y_xyyz_yyzzz,  \
                             tr_y_xyyz_yzzzz,  \
                             tr_y_xyyz_zzzzz,  \
                             tr_y_yyz_yyyyz,   \
                             tr_y_yyz_yyyzz,   \
                             tr_y_yyz_yyzzz,   \
                             tr_y_yyz_yzzzz,   \
                             tr_y_yyz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyz_xxxxx[i] = tr_y_xxyy_xxxxx[i] * pa_z[i];

        tr_y_xxyyz_xxxxy[i] = tr_y_xxyy_xxxxy[i] * pa_z[i];

        tr_y_xxyyz_xxxxz[i] = tr_y_xxyy_xxxx[i] * fe_0 + tr_y_xxyy_xxxxz[i] * pa_z[i];

        tr_y_xxyyz_xxxyy[i] = tr_y_xxyy_xxxyy[i] * pa_z[i];

        tr_y_xxyyz_xxxyz[i] = tr_y_xxyy_xxxy[i] * fe_0 + tr_y_xxyy_xxxyz[i] * pa_z[i];

        tr_y_xxyyz_xxxzz[i] = 2.0 * tr_y_xxyy_xxxz[i] * fe_0 + tr_y_xxyy_xxxzz[i] * pa_z[i];

        tr_y_xxyyz_xxyyy[i] = tr_y_xxyy_xxyyy[i] * pa_z[i];

        tr_y_xxyyz_xxyyz[i] = tr_y_xxyy_xxyy[i] * fe_0 + tr_y_xxyy_xxyyz[i] * pa_z[i];

        tr_y_xxyyz_xxyzz[i] = 2.0 * tr_y_xxyy_xxyz[i] * fe_0 + tr_y_xxyy_xxyzz[i] * pa_z[i];

        tr_y_xxyyz_xxzzz[i] = 3.0 * tr_y_xxyy_xxzz[i] * fe_0 + tr_y_xxyy_xxzzz[i] * pa_z[i];

        tr_y_xxyyz_xyyyy[i] = tr_y_xxyy_xyyyy[i] * pa_z[i];

        tr_y_xxyyz_xyyyz[i] = tr_y_xxyy_xyyy[i] * fe_0 + tr_y_xxyy_xyyyz[i] * pa_z[i];

        tr_y_xxyyz_xyyzz[i] = 2.0 * tr_y_xxyy_xyyz[i] * fe_0 + tr_y_xxyy_xyyzz[i] * pa_z[i];

        tr_y_xxyyz_xyzzz[i] = 3.0 * tr_y_xxyy_xyzz[i] * fe_0 + tr_y_xxyy_xyzzz[i] * pa_z[i];

        tr_y_xxyyz_xzzzz[i] = 4.0 * tr_y_xxyy_xzzz[i] * fe_0 + tr_y_xxyy_xzzzz[i] * pa_z[i];

        tr_y_xxyyz_yyyyy[i] = tr_y_xxyy_yyyyy[i] * pa_z[i];

        tr_y_xxyyz_yyyyz[i] = tr_y_yyz_yyyyz[i] * fe_0 + tr_y_xyyz_yyyyz[i] * pa_x[i];

        tr_y_xxyyz_yyyzz[i] = tr_y_yyz_yyyzz[i] * fe_0 + tr_y_xyyz_yyyzz[i] * pa_x[i];

        tr_y_xxyyz_yyzzz[i] = tr_y_yyz_yyzzz[i] * fe_0 + tr_y_xyyz_yyzzz[i] * pa_x[i];

        tr_y_xxyyz_yzzzz[i] = tr_y_yyz_yzzzz[i] * fe_0 + tr_y_xyyz_yzzzz[i] * pa_x[i];

        tr_y_xxyyz_zzzzz[i] = tr_y_yyz_zzzzz[i] * fe_0 + tr_y_xyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 609-630 components of targeted buffer : HH

    auto tr_y_xxyzz_xxxxx = pbuffer.data(idx_dip_hh + 609);

    auto tr_y_xxyzz_xxxxy = pbuffer.data(idx_dip_hh + 610);

    auto tr_y_xxyzz_xxxxz = pbuffer.data(idx_dip_hh + 611);

    auto tr_y_xxyzz_xxxyy = pbuffer.data(idx_dip_hh + 612);

    auto tr_y_xxyzz_xxxyz = pbuffer.data(idx_dip_hh + 613);

    auto tr_y_xxyzz_xxxzz = pbuffer.data(idx_dip_hh + 614);

    auto tr_y_xxyzz_xxyyy = pbuffer.data(idx_dip_hh + 615);

    auto tr_y_xxyzz_xxyyz = pbuffer.data(idx_dip_hh + 616);

    auto tr_y_xxyzz_xxyzz = pbuffer.data(idx_dip_hh + 617);

    auto tr_y_xxyzz_xxzzz = pbuffer.data(idx_dip_hh + 618);

    auto tr_y_xxyzz_xyyyy = pbuffer.data(idx_dip_hh + 619);

    auto tr_y_xxyzz_xyyyz = pbuffer.data(idx_dip_hh + 620);

    auto tr_y_xxyzz_xyyzz = pbuffer.data(idx_dip_hh + 621);

    auto tr_y_xxyzz_xyzzz = pbuffer.data(idx_dip_hh + 622);

    auto tr_y_xxyzz_xzzzz = pbuffer.data(idx_dip_hh + 623);

    auto tr_y_xxyzz_yyyyy = pbuffer.data(idx_dip_hh + 624);

    auto tr_y_xxyzz_yyyyz = pbuffer.data(idx_dip_hh + 625);

    auto tr_y_xxyzz_yyyzz = pbuffer.data(idx_dip_hh + 626);

    auto tr_y_xxyzz_yyzzz = pbuffer.data(idx_dip_hh + 627);

    auto tr_y_xxyzz_yzzzz = pbuffer.data(idx_dip_hh + 628);

    auto tr_y_xxyzz_zzzzz = pbuffer.data(idx_dip_hh + 629);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_y_xxy_xxxxy,   \
                             tr_y_xxy_xxxyy,   \
                             tr_y_xxy_xxyyy,   \
                             tr_y_xxy_xyyyy,   \
                             tr_y_xxyz_xxxxy,  \
                             tr_y_xxyz_xxxyy,  \
                             tr_y_xxyz_xxyyy,  \
                             tr_y_xxyz_xyyyy,  \
                             tr_y_xxyzz_xxxxx, \
                             tr_y_xxyzz_xxxxy, \
                             tr_y_xxyzz_xxxxz, \
                             tr_y_xxyzz_xxxyy, \
                             tr_y_xxyzz_xxxyz, \
                             tr_y_xxyzz_xxxzz, \
                             tr_y_xxyzz_xxyyy, \
                             tr_y_xxyzz_xxyyz, \
                             tr_y_xxyzz_xxyzz, \
                             tr_y_xxyzz_xxzzz, \
                             tr_y_xxyzz_xyyyy, \
                             tr_y_xxyzz_xyyyz, \
                             tr_y_xxyzz_xyyzz, \
                             tr_y_xxyzz_xyzzz, \
                             tr_y_xxyzz_xzzzz, \
                             tr_y_xxyzz_yyyyy, \
                             tr_y_xxyzz_yyyyz, \
                             tr_y_xxyzz_yyyzz, \
                             tr_y_xxyzz_yyzzz, \
                             tr_y_xxyzz_yzzzz, \
                             tr_y_xxyzz_zzzzz, \
                             tr_y_xxzz_xxxxx,  \
                             tr_y_xxzz_xxxxz,  \
                             tr_y_xxzz_xxxzz,  \
                             tr_y_xxzz_xxzzz,  \
                             tr_y_xxzz_xzzzz,  \
                             tr_y_xyzz_xxxyz,  \
                             tr_y_xyzz_xxyyz,  \
                             tr_y_xyzz_xxyz,   \
                             tr_y_xyzz_xxyzz,  \
                             tr_y_xyzz_xyyyz,  \
                             tr_y_xyzz_xyyz,   \
                             tr_y_xyzz_xyyzz,  \
                             tr_y_xyzz_xyzz,   \
                             tr_y_xyzz_xyzzz,  \
                             tr_y_xyzz_yyyyy,  \
                             tr_y_xyzz_yyyyz,  \
                             tr_y_xyzz_yyyz,   \
                             tr_y_xyzz_yyyzz,  \
                             tr_y_xyzz_yyzz,   \
                             tr_y_xyzz_yyzzz,  \
                             tr_y_xyzz_yzzz,   \
                             tr_y_xyzz_yzzzz,  \
                             tr_y_xyzz_zzzzz,  \
                             tr_y_yzz_xxxyz,   \
                             tr_y_yzz_xxyyz,   \
                             tr_y_yzz_xxyzz,   \
                             tr_y_yzz_xyyyz,   \
                             tr_y_yzz_xyyzz,   \
                             tr_y_yzz_xyzzz,   \
                             tr_y_yzz_yyyyy,   \
                             tr_y_yzz_yyyyz,   \
                             tr_y_yzz_yyyzz,   \
                             tr_y_yzz_yyzzz,   \
                             tr_y_yzz_yzzzz,   \
                             tr_y_yzz_zzzzz,   \
                             ts_xxzz_xxxxx,    \
                             ts_xxzz_xxxxz,    \
                             ts_xxzz_xxxzz,    \
                             ts_xxzz_xxzzz,    \
                             ts_xxzz_xzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzz_xxxxx[i] = ts_xxzz_xxxxx[i] * fe_0 + tr_y_xxzz_xxxxx[i] * pa_y[i];

        tr_y_xxyzz_xxxxy[i] = tr_y_xxy_xxxxy[i] * fe_0 + tr_y_xxyz_xxxxy[i] * pa_z[i];

        tr_y_xxyzz_xxxxz[i] = ts_xxzz_xxxxz[i] * fe_0 + tr_y_xxzz_xxxxz[i] * pa_y[i];

        tr_y_xxyzz_xxxyy[i] = tr_y_xxy_xxxyy[i] * fe_0 + tr_y_xxyz_xxxyy[i] * pa_z[i];

        tr_y_xxyzz_xxxyz[i] = tr_y_yzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xyzz_xxyz[i] * fe_0 + tr_y_xyzz_xxxyz[i] * pa_x[i];

        tr_y_xxyzz_xxxzz[i] = ts_xxzz_xxxzz[i] * fe_0 + tr_y_xxzz_xxxzz[i] * pa_y[i];

        tr_y_xxyzz_xxyyy[i] = tr_y_xxy_xxyyy[i] * fe_0 + tr_y_xxyz_xxyyy[i] * pa_z[i];

        tr_y_xxyzz_xxyyz[i] = tr_y_yzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xyzz_xyyz[i] * fe_0 + tr_y_xyzz_xxyyz[i] * pa_x[i];

        tr_y_xxyzz_xxyzz[i] = tr_y_yzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xyzz_xyzz[i] * fe_0 + tr_y_xyzz_xxyzz[i] * pa_x[i];

        tr_y_xxyzz_xxzzz[i] = ts_xxzz_xxzzz[i] * fe_0 + tr_y_xxzz_xxzzz[i] * pa_y[i];

        tr_y_xxyzz_xyyyy[i] = tr_y_xxy_xyyyy[i] * fe_0 + tr_y_xxyz_xyyyy[i] * pa_z[i];

        tr_y_xxyzz_xyyyz[i] = tr_y_yzz_xyyyz[i] * fe_0 + tr_y_xyzz_yyyz[i] * fe_0 + tr_y_xyzz_xyyyz[i] * pa_x[i];

        tr_y_xxyzz_xyyzz[i] = tr_y_yzz_xyyzz[i] * fe_0 + tr_y_xyzz_yyzz[i] * fe_0 + tr_y_xyzz_xyyzz[i] * pa_x[i];

        tr_y_xxyzz_xyzzz[i] = tr_y_yzz_xyzzz[i] * fe_0 + tr_y_xyzz_yzzz[i] * fe_0 + tr_y_xyzz_xyzzz[i] * pa_x[i];

        tr_y_xxyzz_xzzzz[i] = ts_xxzz_xzzzz[i] * fe_0 + tr_y_xxzz_xzzzz[i] * pa_y[i];

        tr_y_xxyzz_yyyyy[i] = tr_y_yzz_yyyyy[i] * fe_0 + tr_y_xyzz_yyyyy[i] * pa_x[i];

        tr_y_xxyzz_yyyyz[i] = tr_y_yzz_yyyyz[i] * fe_0 + tr_y_xyzz_yyyyz[i] * pa_x[i];

        tr_y_xxyzz_yyyzz[i] = tr_y_yzz_yyyzz[i] * fe_0 + tr_y_xyzz_yyyzz[i] * pa_x[i];

        tr_y_xxyzz_yyzzz[i] = tr_y_yzz_yyzzz[i] * fe_0 + tr_y_xyzz_yyzzz[i] * pa_x[i];

        tr_y_xxyzz_yzzzz[i] = tr_y_yzz_yzzzz[i] * fe_0 + tr_y_xyzz_yzzzz[i] * pa_x[i];

        tr_y_xxyzz_zzzzz[i] = tr_y_yzz_zzzzz[i] * fe_0 + tr_y_xyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 630-651 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xxz_xxxxx,   \
                             tr_y_xxz_xxxxy,   \
                             tr_y_xxz_xxxyy,   \
                             tr_y_xxz_xxyyy,   \
                             tr_y_xxz_xyyyy,   \
                             tr_y_xxzz_xxxxx,  \
                             tr_y_xxzz_xxxxy,  \
                             tr_y_xxzz_xxxyy,  \
                             tr_y_xxzz_xxyyy,  \
                             tr_y_xxzz_xyyyy,  \
                             tr_y_xxzzz_xxxxx, \
                             tr_y_xxzzz_xxxxy, \
                             tr_y_xxzzz_xxxxz, \
                             tr_y_xxzzz_xxxyy, \
                             tr_y_xxzzz_xxxyz, \
                             tr_y_xxzzz_xxxzz, \
                             tr_y_xxzzz_xxyyy, \
                             tr_y_xxzzz_xxyyz, \
                             tr_y_xxzzz_xxyzz, \
                             tr_y_xxzzz_xxzzz, \
                             tr_y_xxzzz_xyyyy, \
                             tr_y_xxzzz_xyyyz, \
                             tr_y_xxzzz_xyyzz, \
                             tr_y_xxzzz_xyzzz, \
                             tr_y_xxzzz_xzzzz, \
                             tr_y_xxzzz_yyyyy, \
                             tr_y_xxzzz_yyyyz, \
                             tr_y_xxzzz_yyyzz, \
                             tr_y_xxzzz_yyzzz, \
                             tr_y_xxzzz_yzzzz, \
                             tr_y_xxzzz_zzzzz, \
                             tr_y_xzzz_xxxxz,  \
                             tr_y_xzzz_xxxyz,  \
                             tr_y_xzzz_xxxz,   \
                             tr_y_xzzz_xxxzz,  \
                             tr_y_xzzz_xxyyz,  \
                             tr_y_xzzz_xxyz,   \
                             tr_y_xzzz_xxyzz,  \
                             tr_y_xzzz_xxzz,   \
                             tr_y_xzzz_xxzzz,  \
                             tr_y_xzzz_xyyyz,  \
                             tr_y_xzzz_xyyz,   \
                             tr_y_xzzz_xyyzz,  \
                             tr_y_xzzz_xyzz,   \
                             tr_y_xzzz_xyzzz,  \
                             tr_y_xzzz_xzzz,   \
                             tr_y_xzzz_xzzzz,  \
                             tr_y_xzzz_yyyyy,  \
                             tr_y_xzzz_yyyyz,  \
                             tr_y_xzzz_yyyz,   \
                             tr_y_xzzz_yyyzz,  \
                             tr_y_xzzz_yyzz,   \
                             tr_y_xzzz_yyzzz,  \
                             tr_y_xzzz_yzzz,   \
                             tr_y_xzzz_yzzzz,  \
                             tr_y_xzzz_zzzz,   \
                             tr_y_xzzz_zzzzz,  \
                             tr_y_zzz_xxxxz,   \
                             tr_y_zzz_xxxyz,   \
                             tr_y_zzz_xxxzz,   \
                             tr_y_zzz_xxyyz,   \
                             tr_y_zzz_xxyzz,   \
                             tr_y_zzz_xxzzz,   \
                             tr_y_zzz_xyyyz,   \
                             tr_y_zzz_xyyzz,   \
                             tr_y_zzz_xyzzz,   \
                             tr_y_zzz_xzzzz,   \
                             tr_y_zzz_yyyyy,   \
                             tr_y_zzz_yyyyz,   \
                             tr_y_zzz_yyyzz,   \
                             tr_y_zzz_yyzzz,   \
                             tr_y_zzz_yzzzz,   \
                             tr_y_zzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzz_xxxxx[i] = 2.0 * tr_y_xxz_xxxxx[i] * fe_0 + tr_y_xxzz_xxxxx[i] * pa_z[i];

        tr_y_xxzzz_xxxxy[i] = 2.0 * tr_y_xxz_xxxxy[i] * fe_0 + tr_y_xxzz_xxxxy[i] * pa_z[i];

        tr_y_xxzzz_xxxxz[i] = tr_y_zzz_xxxxz[i] * fe_0 + 4.0 * tr_y_xzzz_xxxz[i] * fe_0 + tr_y_xzzz_xxxxz[i] * pa_x[i];

        tr_y_xxzzz_xxxyy[i] = 2.0 * tr_y_xxz_xxxyy[i] * fe_0 + tr_y_xxzz_xxxyy[i] * pa_z[i];

        tr_y_xxzzz_xxxyz[i] = tr_y_zzz_xxxyz[i] * fe_0 + 3.0 * tr_y_xzzz_xxyz[i] * fe_0 + tr_y_xzzz_xxxyz[i] * pa_x[i];

        tr_y_xxzzz_xxxzz[i] = tr_y_zzz_xxxzz[i] * fe_0 + 3.0 * tr_y_xzzz_xxzz[i] * fe_0 + tr_y_xzzz_xxxzz[i] * pa_x[i];

        tr_y_xxzzz_xxyyy[i] = 2.0 * tr_y_xxz_xxyyy[i] * fe_0 + tr_y_xxzz_xxyyy[i] * pa_z[i];

        tr_y_xxzzz_xxyyz[i] = tr_y_zzz_xxyyz[i] * fe_0 + 2.0 * tr_y_xzzz_xyyz[i] * fe_0 + tr_y_xzzz_xxyyz[i] * pa_x[i];

        tr_y_xxzzz_xxyzz[i] = tr_y_zzz_xxyzz[i] * fe_0 + 2.0 * tr_y_xzzz_xyzz[i] * fe_0 + tr_y_xzzz_xxyzz[i] * pa_x[i];

        tr_y_xxzzz_xxzzz[i] = tr_y_zzz_xxzzz[i] * fe_0 + 2.0 * tr_y_xzzz_xzzz[i] * fe_0 + tr_y_xzzz_xxzzz[i] * pa_x[i];

        tr_y_xxzzz_xyyyy[i] = 2.0 * tr_y_xxz_xyyyy[i] * fe_0 + tr_y_xxzz_xyyyy[i] * pa_z[i];

        tr_y_xxzzz_xyyyz[i] = tr_y_zzz_xyyyz[i] * fe_0 + tr_y_xzzz_yyyz[i] * fe_0 + tr_y_xzzz_xyyyz[i] * pa_x[i];

        tr_y_xxzzz_xyyzz[i] = tr_y_zzz_xyyzz[i] * fe_0 + tr_y_xzzz_yyzz[i] * fe_0 + tr_y_xzzz_xyyzz[i] * pa_x[i];

        tr_y_xxzzz_xyzzz[i] = tr_y_zzz_xyzzz[i] * fe_0 + tr_y_xzzz_yzzz[i] * fe_0 + tr_y_xzzz_xyzzz[i] * pa_x[i];

        tr_y_xxzzz_xzzzz[i] = tr_y_zzz_xzzzz[i] * fe_0 + tr_y_xzzz_zzzz[i] * fe_0 + tr_y_xzzz_xzzzz[i] * pa_x[i];

        tr_y_xxzzz_yyyyy[i] = tr_y_zzz_yyyyy[i] * fe_0 + tr_y_xzzz_yyyyy[i] * pa_x[i];

        tr_y_xxzzz_yyyyz[i] = tr_y_zzz_yyyyz[i] * fe_0 + tr_y_xzzz_yyyyz[i] * pa_x[i];

        tr_y_xxzzz_yyyzz[i] = tr_y_zzz_yyyzz[i] * fe_0 + tr_y_xzzz_yyyzz[i] * pa_x[i];

        tr_y_xxzzz_yyzzz[i] = tr_y_zzz_yyzzz[i] * fe_0 + tr_y_xzzz_yyzzz[i] * pa_x[i];

        tr_y_xxzzz_yzzzz[i] = tr_y_zzz_yzzzz[i] * fe_0 + tr_y_xzzz_yzzzz[i] * pa_x[i];

        tr_y_xxzzz_zzzzz[i] = tr_y_zzz_zzzzz[i] * fe_0 + tr_y_xzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 651-672 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyyyy_xxxxx, \
                             tr_y_xyyyy_xxxxy, \
                             tr_y_xyyyy_xxxxz, \
                             tr_y_xyyyy_xxxyy, \
                             tr_y_xyyyy_xxxyz, \
                             tr_y_xyyyy_xxxzz, \
                             tr_y_xyyyy_xxyyy, \
                             tr_y_xyyyy_xxyyz, \
                             tr_y_xyyyy_xxyzz, \
                             tr_y_xyyyy_xxzzz, \
                             tr_y_xyyyy_xyyyy, \
                             tr_y_xyyyy_xyyyz, \
                             tr_y_xyyyy_xyyzz, \
                             tr_y_xyyyy_xyzzz, \
                             tr_y_xyyyy_xzzzz, \
                             tr_y_xyyyy_yyyyy, \
                             tr_y_xyyyy_yyyyz, \
                             tr_y_xyyyy_yyyzz, \
                             tr_y_xyyyy_yyzzz, \
                             tr_y_xyyyy_yzzzz, \
                             tr_y_xyyyy_zzzzz, \
                             tr_y_yyyy_xxxx,   \
                             tr_y_yyyy_xxxxx,  \
                             tr_y_yyyy_xxxxy,  \
                             tr_y_yyyy_xxxxz,  \
                             tr_y_yyyy_xxxy,   \
                             tr_y_yyyy_xxxyy,  \
                             tr_y_yyyy_xxxyz,  \
                             tr_y_yyyy_xxxz,   \
                             tr_y_yyyy_xxxzz,  \
                             tr_y_yyyy_xxyy,   \
                             tr_y_yyyy_xxyyy,  \
                             tr_y_yyyy_xxyyz,  \
                             tr_y_yyyy_xxyz,   \
                             tr_y_yyyy_xxyzz,  \
                             tr_y_yyyy_xxzz,   \
                             tr_y_yyyy_xxzzz,  \
                             tr_y_yyyy_xyyy,   \
                             tr_y_yyyy_xyyyy,  \
                             tr_y_yyyy_xyyyz,  \
                             tr_y_yyyy_xyyz,   \
                             tr_y_yyyy_xyyzz,  \
                             tr_y_yyyy_xyzz,   \
                             tr_y_yyyy_xyzzz,  \
                             tr_y_yyyy_xzzz,   \
                             tr_y_yyyy_xzzzz,  \
                             tr_y_yyyy_yyyy,   \
                             tr_y_yyyy_yyyyy,  \
                             tr_y_yyyy_yyyyz,  \
                             tr_y_yyyy_yyyz,   \
                             tr_y_yyyy_yyyzz,  \
                             tr_y_yyyy_yyzz,   \
                             tr_y_yyyy_yyzzz,  \
                             tr_y_yyyy_yzzz,   \
                             tr_y_yyyy_yzzzz,  \
                             tr_y_yyyy_zzzz,   \
                             tr_y_yyyy_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyy_xxxxx[i] = 5.0 * tr_y_yyyy_xxxx[i] * fe_0 + tr_y_yyyy_xxxxx[i] * pa_x[i];

        tr_y_xyyyy_xxxxy[i] = 4.0 * tr_y_yyyy_xxxy[i] * fe_0 + tr_y_yyyy_xxxxy[i] * pa_x[i];

        tr_y_xyyyy_xxxxz[i] = 4.0 * tr_y_yyyy_xxxz[i] * fe_0 + tr_y_yyyy_xxxxz[i] * pa_x[i];

        tr_y_xyyyy_xxxyy[i] = 3.0 * tr_y_yyyy_xxyy[i] * fe_0 + tr_y_yyyy_xxxyy[i] * pa_x[i];

        tr_y_xyyyy_xxxyz[i] = 3.0 * tr_y_yyyy_xxyz[i] * fe_0 + tr_y_yyyy_xxxyz[i] * pa_x[i];

        tr_y_xyyyy_xxxzz[i] = 3.0 * tr_y_yyyy_xxzz[i] * fe_0 + tr_y_yyyy_xxxzz[i] * pa_x[i];

        tr_y_xyyyy_xxyyy[i] = 2.0 * tr_y_yyyy_xyyy[i] * fe_0 + tr_y_yyyy_xxyyy[i] * pa_x[i];

        tr_y_xyyyy_xxyyz[i] = 2.0 * tr_y_yyyy_xyyz[i] * fe_0 + tr_y_yyyy_xxyyz[i] * pa_x[i];

        tr_y_xyyyy_xxyzz[i] = 2.0 * tr_y_yyyy_xyzz[i] * fe_0 + tr_y_yyyy_xxyzz[i] * pa_x[i];

        tr_y_xyyyy_xxzzz[i] = 2.0 * tr_y_yyyy_xzzz[i] * fe_0 + tr_y_yyyy_xxzzz[i] * pa_x[i];

        tr_y_xyyyy_xyyyy[i] = tr_y_yyyy_yyyy[i] * fe_0 + tr_y_yyyy_xyyyy[i] * pa_x[i];

        tr_y_xyyyy_xyyyz[i] = tr_y_yyyy_yyyz[i] * fe_0 + tr_y_yyyy_xyyyz[i] * pa_x[i];

        tr_y_xyyyy_xyyzz[i] = tr_y_yyyy_yyzz[i] * fe_0 + tr_y_yyyy_xyyzz[i] * pa_x[i];

        tr_y_xyyyy_xyzzz[i] = tr_y_yyyy_yzzz[i] * fe_0 + tr_y_yyyy_xyzzz[i] * pa_x[i];

        tr_y_xyyyy_xzzzz[i] = tr_y_yyyy_zzzz[i] * fe_0 + tr_y_yyyy_xzzzz[i] * pa_x[i];

        tr_y_xyyyy_yyyyy[i] = tr_y_yyyy_yyyyy[i] * pa_x[i];

        tr_y_xyyyy_yyyyz[i] = tr_y_yyyy_yyyyz[i] * pa_x[i];

        tr_y_xyyyy_yyyzz[i] = tr_y_yyyy_yyyzz[i] * pa_x[i];

        tr_y_xyyyy_yyzzz[i] = tr_y_yyyy_yyzzz[i] * pa_x[i];

        tr_y_xyyyy_yzzzz[i] = tr_y_yyyy_yzzzz[i] * pa_x[i];

        tr_y_xyyyy_zzzzz[i] = tr_y_yyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 672-693 components of targeted buffer : HH

    auto tr_y_xyyyz_xxxxx = pbuffer.data(idx_dip_hh + 672);

    auto tr_y_xyyyz_xxxxy = pbuffer.data(idx_dip_hh + 673);

    auto tr_y_xyyyz_xxxxz = pbuffer.data(idx_dip_hh + 674);

    auto tr_y_xyyyz_xxxyy = pbuffer.data(idx_dip_hh + 675);

    auto tr_y_xyyyz_xxxyz = pbuffer.data(idx_dip_hh + 676);

    auto tr_y_xyyyz_xxxzz = pbuffer.data(idx_dip_hh + 677);

    auto tr_y_xyyyz_xxyyy = pbuffer.data(idx_dip_hh + 678);

    auto tr_y_xyyyz_xxyyz = pbuffer.data(idx_dip_hh + 679);

    auto tr_y_xyyyz_xxyzz = pbuffer.data(idx_dip_hh + 680);

    auto tr_y_xyyyz_xxzzz = pbuffer.data(idx_dip_hh + 681);

    auto tr_y_xyyyz_xyyyy = pbuffer.data(idx_dip_hh + 682);

    auto tr_y_xyyyz_xyyyz = pbuffer.data(idx_dip_hh + 683);

    auto tr_y_xyyyz_xyyzz = pbuffer.data(idx_dip_hh + 684);

    auto tr_y_xyyyz_xyzzz = pbuffer.data(idx_dip_hh + 685);

    auto tr_y_xyyyz_xzzzz = pbuffer.data(idx_dip_hh + 686);

    auto tr_y_xyyyz_yyyyy = pbuffer.data(idx_dip_hh + 687);

    auto tr_y_xyyyz_yyyyz = pbuffer.data(idx_dip_hh + 688);

    auto tr_y_xyyyz_yyyzz = pbuffer.data(idx_dip_hh + 689);

    auto tr_y_xyyyz_yyzzz = pbuffer.data(idx_dip_hh + 690);

    auto tr_y_xyyyz_yzzzz = pbuffer.data(idx_dip_hh + 691);

    auto tr_y_xyyyz_zzzzz = pbuffer.data(idx_dip_hh + 692);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_y_xyyy_xxxxx,  \
                             tr_y_xyyy_xxxxy,  \
                             tr_y_xyyy_xxxyy,  \
                             tr_y_xyyy_xxyyy,  \
                             tr_y_xyyy_xyyyy,  \
                             tr_y_xyyyz_xxxxx, \
                             tr_y_xyyyz_xxxxy, \
                             tr_y_xyyyz_xxxxz, \
                             tr_y_xyyyz_xxxyy, \
                             tr_y_xyyyz_xxxyz, \
                             tr_y_xyyyz_xxxzz, \
                             tr_y_xyyyz_xxyyy, \
                             tr_y_xyyyz_xxyyz, \
                             tr_y_xyyyz_xxyzz, \
                             tr_y_xyyyz_xxzzz, \
                             tr_y_xyyyz_xyyyy, \
                             tr_y_xyyyz_xyyyz, \
                             tr_y_xyyyz_xyyzz, \
                             tr_y_xyyyz_xyzzz, \
                             tr_y_xyyyz_xzzzz, \
                             tr_y_xyyyz_yyyyy, \
                             tr_y_xyyyz_yyyyz, \
                             tr_y_xyyyz_yyyzz, \
                             tr_y_xyyyz_yyzzz, \
                             tr_y_xyyyz_yzzzz, \
                             tr_y_xyyyz_zzzzz, \
                             tr_y_yyyz_xxxxz,  \
                             tr_y_yyyz_xxxyz,  \
                             tr_y_yyyz_xxxz,   \
                             tr_y_yyyz_xxxzz,  \
                             tr_y_yyyz_xxyyz,  \
                             tr_y_yyyz_xxyz,   \
                             tr_y_yyyz_xxyzz,  \
                             tr_y_yyyz_xxzz,   \
                             tr_y_yyyz_xxzzz,  \
                             tr_y_yyyz_xyyyz,  \
                             tr_y_yyyz_xyyz,   \
                             tr_y_yyyz_xyyzz,  \
                             tr_y_yyyz_xyzz,   \
                             tr_y_yyyz_xyzzz,  \
                             tr_y_yyyz_xzzz,   \
                             tr_y_yyyz_xzzzz,  \
                             tr_y_yyyz_yyyyy,  \
                             tr_y_yyyz_yyyyz,  \
                             tr_y_yyyz_yyyz,   \
                             tr_y_yyyz_yyyzz,  \
                             tr_y_yyyz_yyzz,   \
                             tr_y_yyyz_yyzzz,  \
                             tr_y_yyyz_yzzz,   \
                             tr_y_yyyz_yzzzz,  \
                             tr_y_yyyz_zzzz,   \
                             tr_y_yyyz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyz_xxxxx[i] = tr_y_xyyy_xxxxx[i] * pa_z[i];

        tr_y_xyyyz_xxxxy[i] = tr_y_xyyy_xxxxy[i] * pa_z[i];

        tr_y_xyyyz_xxxxz[i] = 4.0 * tr_y_yyyz_xxxz[i] * fe_0 + tr_y_yyyz_xxxxz[i] * pa_x[i];

        tr_y_xyyyz_xxxyy[i] = tr_y_xyyy_xxxyy[i] * pa_z[i];

        tr_y_xyyyz_xxxyz[i] = 3.0 * tr_y_yyyz_xxyz[i] * fe_0 + tr_y_yyyz_xxxyz[i] * pa_x[i];

        tr_y_xyyyz_xxxzz[i] = 3.0 * tr_y_yyyz_xxzz[i] * fe_0 + tr_y_yyyz_xxxzz[i] * pa_x[i];

        tr_y_xyyyz_xxyyy[i] = tr_y_xyyy_xxyyy[i] * pa_z[i];

        tr_y_xyyyz_xxyyz[i] = 2.0 * tr_y_yyyz_xyyz[i] * fe_0 + tr_y_yyyz_xxyyz[i] * pa_x[i];

        tr_y_xyyyz_xxyzz[i] = 2.0 * tr_y_yyyz_xyzz[i] * fe_0 + tr_y_yyyz_xxyzz[i] * pa_x[i];

        tr_y_xyyyz_xxzzz[i] = 2.0 * tr_y_yyyz_xzzz[i] * fe_0 + tr_y_yyyz_xxzzz[i] * pa_x[i];

        tr_y_xyyyz_xyyyy[i] = tr_y_xyyy_xyyyy[i] * pa_z[i];

        tr_y_xyyyz_xyyyz[i] = tr_y_yyyz_yyyz[i] * fe_0 + tr_y_yyyz_xyyyz[i] * pa_x[i];

        tr_y_xyyyz_xyyzz[i] = tr_y_yyyz_yyzz[i] * fe_0 + tr_y_yyyz_xyyzz[i] * pa_x[i];

        tr_y_xyyyz_xyzzz[i] = tr_y_yyyz_yzzz[i] * fe_0 + tr_y_yyyz_xyzzz[i] * pa_x[i];

        tr_y_xyyyz_xzzzz[i] = tr_y_yyyz_zzzz[i] * fe_0 + tr_y_yyyz_xzzzz[i] * pa_x[i];

        tr_y_xyyyz_yyyyy[i] = tr_y_yyyz_yyyyy[i] * pa_x[i];

        tr_y_xyyyz_yyyyz[i] = tr_y_yyyz_yyyyz[i] * pa_x[i];

        tr_y_xyyyz_yyyzz[i] = tr_y_yyyz_yyyzz[i] * pa_x[i];

        tr_y_xyyyz_yyzzz[i] = tr_y_yyyz_yyzzz[i] * pa_x[i];

        tr_y_xyyyz_yzzzz[i] = tr_y_yyyz_yzzzz[i] * pa_x[i];

        tr_y_xyyyz_zzzzz[i] = tr_y_yyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 693-714 components of targeted buffer : HH

    auto tr_y_xyyzz_xxxxx = pbuffer.data(idx_dip_hh + 693);

    auto tr_y_xyyzz_xxxxy = pbuffer.data(idx_dip_hh + 694);

    auto tr_y_xyyzz_xxxxz = pbuffer.data(idx_dip_hh + 695);

    auto tr_y_xyyzz_xxxyy = pbuffer.data(idx_dip_hh + 696);

    auto tr_y_xyyzz_xxxyz = pbuffer.data(idx_dip_hh + 697);

    auto tr_y_xyyzz_xxxzz = pbuffer.data(idx_dip_hh + 698);

    auto tr_y_xyyzz_xxyyy = pbuffer.data(idx_dip_hh + 699);

    auto tr_y_xyyzz_xxyyz = pbuffer.data(idx_dip_hh + 700);

    auto tr_y_xyyzz_xxyzz = pbuffer.data(idx_dip_hh + 701);

    auto tr_y_xyyzz_xxzzz = pbuffer.data(idx_dip_hh + 702);

    auto tr_y_xyyzz_xyyyy = pbuffer.data(idx_dip_hh + 703);

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

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyyzz_xxxxx, \
                             tr_y_xyyzz_xxxxy, \
                             tr_y_xyyzz_xxxxz, \
                             tr_y_xyyzz_xxxyy, \
                             tr_y_xyyzz_xxxyz, \
                             tr_y_xyyzz_xxxzz, \
                             tr_y_xyyzz_xxyyy, \
                             tr_y_xyyzz_xxyyz, \
                             tr_y_xyyzz_xxyzz, \
                             tr_y_xyyzz_xxzzz, \
                             tr_y_xyyzz_xyyyy, \
                             tr_y_xyyzz_xyyyz, \
                             tr_y_xyyzz_xyyzz, \
                             tr_y_xyyzz_xyzzz, \
                             tr_y_xyyzz_xzzzz, \
                             tr_y_xyyzz_yyyyy, \
                             tr_y_xyyzz_yyyyz, \
                             tr_y_xyyzz_yyyzz, \
                             tr_y_xyyzz_yyzzz, \
                             tr_y_xyyzz_yzzzz, \
                             tr_y_xyyzz_zzzzz, \
                             tr_y_yyzz_xxxx,   \
                             tr_y_yyzz_xxxxx,  \
                             tr_y_yyzz_xxxxy,  \
                             tr_y_yyzz_xxxxz,  \
                             tr_y_yyzz_xxxy,   \
                             tr_y_yyzz_xxxyy,  \
                             tr_y_yyzz_xxxyz,  \
                             tr_y_yyzz_xxxz,   \
                             tr_y_yyzz_xxxzz,  \
                             tr_y_yyzz_xxyy,   \
                             tr_y_yyzz_xxyyy,  \
                             tr_y_yyzz_xxyyz,  \
                             tr_y_yyzz_xxyz,   \
                             tr_y_yyzz_xxyzz,  \
                             tr_y_yyzz_xxzz,   \
                             tr_y_yyzz_xxzzz,  \
                             tr_y_yyzz_xyyy,   \
                             tr_y_yyzz_xyyyy,  \
                             tr_y_yyzz_xyyyz,  \
                             tr_y_yyzz_xyyz,   \
                             tr_y_yyzz_xyyzz,  \
                             tr_y_yyzz_xyzz,   \
                             tr_y_yyzz_xyzzz,  \
                             tr_y_yyzz_xzzz,   \
                             tr_y_yyzz_xzzzz,  \
                             tr_y_yyzz_yyyy,   \
                             tr_y_yyzz_yyyyy,  \
                             tr_y_yyzz_yyyyz,  \
                             tr_y_yyzz_yyyz,   \
                             tr_y_yyzz_yyyzz,  \
                             tr_y_yyzz_yyzz,   \
                             tr_y_yyzz_yyzzz,  \
                             tr_y_yyzz_yzzz,   \
                             tr_y_yyzz_yzzzz,  \
                             tr_y_yyzz_zzzz,   \
                             tr_y_yyzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzz_xxxxx[i] = 5.0 * tr_y_yyzz_xxxx[i] * fe_0 + tr_y_yyzz_xxxxx[i] * pa_x[i];

        tr_y_xyyzz_xxxxy[i] = 4.0 * tr_y_yyzz_xxxy[i] * fe_0 + tr_y_yyzz_xxxxy[i] * pa_x[i];

        tr_y_xyyzz_xxxxz[i] = 4.0 * tr_y_yyzz_xxxz[i] * fe_0 + tr_y_yyzz_xxxxz[i] * pa_x[i];

        tr_y_xyyzz_xxxyy[i] = 3.0 * tr_y_yyzz_xxyy[i] * fe_0 + tr_y_yyzz_xxxyy[i] * pa_x[i];

        tr_y_xyyzz_xxxyz[i] = 3.0 * tr_y_yyzz_xxyz[i] * fe_0 + tr_y_yyzz_xxxyz[i] * pa_x[i];

        tr_y_xyyzz_xxxzz[i] = 3.0 * tr_y_yyzz_xxzz[i] * fe_0 + tr_y_yyzz_xxxzz[i] * pa_x[i];

        tr_y_xyyzz_xxyyy[i] = 2.0 * tr_y_yyzz_xyyy[i] * fe_0 + tr_y_yyzz_xxyyy[i] * pa_x[i];

        tr_y_xyyzz_xxyyz[i] = 2.0 * tr_y_yyzz_xyyz[i] * fe_0 + tr_y_yyzz_xxyyz[i] * pa_x[i];

        tr_y_xyyzz_xxyzz[i] = 2.0 * tr_y_yyzz_xyzz[i] * fe_0 + tr_y_yyzz_xxyzz[i] * pa_x[i];

        tr_y_xyyzz_xxzzz[i] = 2.0 * tr_y_yyzz_xzzz[i] * fe_0 + tr_y_yyzz_xxzzz[i] * pa_x[i];

        tr_y_xyyzz_xyyyy[i] = tr_y_yyzz_yyyy[i] * fe_0 + tr_y_yyzz_xyyyy[i] * pa_x[i];

        tr_y_xyyzz_xyyyz[i] = tr_y_yyzz_yyyz[i] * fe_0 + tr_y_yyzz_xyyyz[i] * pa_x[i];

        tr_y_xyyzz_xyyzz[i] = tr_y_yyzz_yyzz[i] * fe_0 + tr_y_yyzz_xyyzz[i] * pa_x[i];

        tr_y_xyyzz_xyzzz[i] = tr_y_yyzz_yzzz[i] * fe_0 + tr_y_yyzz_xyzzz[i] * pa_x[i];

        tr_y_xyyzz_xzzzz[i] = tr_y_yyzz_zzzz[i] * fe_0 + tr_y_yyzz_xzzzz[i] * pa_x[i];

        tr_y_xyyzz_yyyyy[i] = tr_y_yyzz_yyyyy[i] * pa_x[i];

        tr_y_xyyzz_yyyyz[i] = tr_y_yyzz_yyyyz[i] * pa_x[i];

        tr_y_xyyzz_yyyzz[i] = tr_y_yyzz_yyyzz[i] * pa_x[i];

        tr_y_xyyzz_yyzzz[i] = tr_y_yyzz_yyzzz[i] * pa_x[i];

        tr_y_xyyzz_yzzzz[i] = tr_y_yyzz_yzzzz[i] * pa_x[i];

        tr_y_xyyzz_zzzzz[i] = tr_y_yyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 714-735 components of targeted buffer : HH

    auto tr_y_xyzzz_xxxxx = pbuffer.data(idx_dip_hh + 714);

    auto tr_y_xyzzz_xxxxy = pbuffer.data(idx_dip_hh + 715);

    auto tr_y_xyzzz_xxxxz = pbuffer.data(idx_dip_hh + 716);

    auto tr_y_xyzzz_xxxyy = pbuffer.data(idx_dip_hh + 717);

    auto tr_y_xyzzz_xxxyz = pbuffer.data(idx_dip_hh + 718);

    auto tr_y_xyzzz_xxxzz = pbuffer.data(idx_dip_hh + 719);

    auto tr_y_xyzzz_xxyyy = pbuffer.data(idx_dip_hh + 720);

    auto tr_y_xyzzz_xxyyz = pbuffer.data(idx_dip_hh + 721);

    auto tr_y_xyzzz_xxyzz = pbuffer.data(idx_dip_hh + 722);

    auto tr_y_xyzzz_xxzzz = pbuffer.data(idx_dip_hh + 723);

    auto tr_y_xyzzz_xyyyy = pbuffer.data(idx_dip_hh + 724);

    auto tr_y_xyzzz_xyyyz = pbuffer.data(idx_dip_hh + 725);

    auto tr_y_xyzzz_xyyzz = pbuffer.data(idx_dip_hh + 726);

    auto tr_y_xyzzz_xyzzz = pbuffer.data(idx_dip_hh + 727);

    auto tr_y_xyzzz_xzzzz = pbuffer.data(idx_dip_hh + 728);

    auto tr_y_xyzzz_yyyyy = pbuffer.data(idx_dip_hh + 729);

    auto tr_y_xyzzz_yyyyz = pbuffer.data(idx_dip_hh + 730);

    auto tr_y_xyzzz_yyyzz = pbuffer.data(idx_dip_hh + 731);

    auto tr_y_xyzzz_yyzzz = pbuffer.data(idx_dip_hh + 732);

    auto tr_y_xyzzz_yzzzz = pbuffer.data(idx_dip_hh + 733);

    auto tr_y_xyzzz_zzzzz = pbuffer.data(idx_dip_hh + 734);

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xyzzz_xxxxx, \
                             tr_y_xyzzz_xxxxy, \
                             tr_y_xyzzz_xxxxz, \
                             tr_y_xyzzz_xxxyy, \
                             tr_y_xyzzz_xxxyz, \
                             tr_y_xyzzz_xxxzz, \
                             tr_y_xyzzz_xxyyy, \
                             tr_y_xyzzz_xxyyz, \
                             tr_y_xyzzz_xxyzz, \
                             tr_y_xyzzz_xxzzz, \
                             tr_y_xyzzz_xyyyy, \
                             tr_y_xyzzz_xyyyz, \
                             tr_y_xyzzz_xyyzz, \
                             tr_y_xyzzz_xyzzz, \
                             tr_y_xyzzz_xzzzz, \
                             tr_y_xyzzz_yyyyy, \
                             tr_y_xyzzz_yyyyz, \
                             tr_y_xyzzz_yyyzz, \
                             tr_y_xyzzz_yyzzz, \
                             tr_y_xyzzz_yzzzz, \
                             tr_y_xyzzz_zzzzz, \
                             tr_y_yzzz_xxxx,   \
                             tr_y_yzzz_xxxxx,  \
                             tr_y_yzzz_xxxxy,  \
                             tr_y_yzzz_xxxxz,  \
                             tr_y_yzzz_xxxy,   \
                             tr_y_yzzz_xxxyy,  \
                             tr_y_yzzz_xxxyz,  \
                             tr_y_yzzz_xxxz,   \
                             tr_y_yzzz_xxxzz,  \
                             tr_y_yzzz_xxyy,   \
                             tr_y_yzzz_xxyyy,  \
                             tr_y_yzzz_xxyyz,  \
                             tr_y_yzzz_xxyz,   \
                             tr_y_yzzz_xxyzz,  \
                             tr_y_yzzz_xxzz,   \
                             tr_y_yzzz_xxzzz,  \
                             tr_y_yzzz_xyyy,   \
                             tr_y_yzzz_xyyyy,  \
                             tr_y_yzzz_xyyyz,  \
                             tr_y_yzzz_xyyz,   \
                             tr_y_yzzz_xyyzz,  \
                             tr_y_yzzz_xyzz,   \
                             tr_y_yzzz_xyzzz,  \
                             tr_y_yzzz_xzzz,   \
                             tr_y_yzzz_xzzzz,  \
                             tr_y_yzzz_yyyy,   \
                             tr_y_yzzz_yyyyy,  \
                             tr_y_yzzz_yyyyz,  \
                             tr_y_yzzz_yyyz,   \
                             tr_y_yzzz_yyyzz,  \
                             tr_y_yzzz_yyzz,   \
                             tr_y_yzzz_yyzzz,  \
                             tr_y_yzzz_yzzz,   \
                             tr_y_yzzz_yzzzz,  \
                             tr_y_yzzz_zzzz,   \
                             tr_y_yzzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzz_xxxxx[i] = 5.0 * tr_y_yzzz_xxxx[i] * fe_0 + tr_y_yzzz_xxxxx[i] * pa_x[i];

        tr_y_xyzzz_xxxxy[i] = 4.0 * tr_y_yzzz_xxxy[i] * fe_0 + tr_y_yzzz_xxxxy[i] * pa_x[i];

        tr_y_xyzzz_xxxxz[i] = 4.0 * tr_y_yzzz_xxxz[i] * fe_0 + tr_y_yzzz_xxxxz[i] * pa_x[i];

        tr_y_xyzzz_xxxyy[i] = 3.0 * tr_y_yzzz_xxyy[i] * fe_0 + tr_y_yzzz_xxxyy[i] * pa_x[i];

        tr_y_xyzzz_xxxyz[i] = 3.0 * tr_y_yzzz_xxyz[i] * fe_0 + tr_y_yzzz_xxxyz[i] * pa_x[i];

        tr_y_xyzzz_xxxzz[i] = 3.0 * tr_y_yzzz_xxzz[i] * fe_0 + tr_y_yzzz_xxxzz[i] * pa_x[i];

        tr_y_xyzzz_xxyyy[i] = 2.0 * tr_y_yzzz_xyyy[i] * fe_0 + tr_y_yzzz_xxyyy[i] * pa_x[i];

        tr_y_xyzzz_xxyyz[i] = 2.0 * tr_y_yzzz_xyyz[i] * fe_0 + tr_y_yzzz_xxyyz[i] * pa_x[i];

        tr_y_xyzzz_xxyzz[i] = 2.0 * tr_y_yzzz_xyzz[i] * fe_0 + tr_y_yzzz_xxyzz[i] * pa_x[i];

        tr_y_xyzzz_xxzzz[i] = 2.0 * tr_y_yzzz_xzzz[i] * fe_0 + tr_y_yzzz_xxzzz[i] * pa_x[i];

        tr_y_xyzzz_xyyyy[i] = tr_y_yzzz_yyyy[i] * fe_0 + tr_y_yzzz_xyyyy[i] * pa_x[i];

        tr_y_xyzzz_xyyyz[i] = tr_y_yzzz_yyyz[i] * fe_0 + tr_y_yzzz_xyyyz[i] * pa_x[i];

        tr_y_xyzzz_xyyzz[i] = tr_y_yzzz_yyzz[i] * fe_0 + tr_y_yzzz_xyyzz[i] * pa_x[i];

        tr_y_xyzzz_xyzzz[i] = tr_y_yzzz_yzzz[i] * fe_0 + tr_y_yzzz_xyzzz[i] * pa_x[i];

        tr_y_xyzzz_xzzzz[i] = tr_y_yzzz_zzzz[i] * fe_0 + tr_y_yzzz_xzzzz[i] * pa_x[i];

        tr_y_xyzzz_yyyyy[i] = tr_y_yzzz_yyyyy[i] * pa_x[i];

        tr_y_xyzzz_yyyyz[i] = tr_y_yzzz_yyyyz[i] * pa_x[i];

        tr_y_xyzzz_yyyzz[i] = tr_y_yzzz_yyyzz[i] * pa_x[i];

        tr_y_xyzzz_yyzzz[i] = tr_y_yzzz_yyzzz[i] * pa_x[i];

        tr_y_xyzzz_yzzzz[i] = tr_y_yzzz_yzzzz[i] * pa_x[i];

        tr_y_xyzzz_zzzzz[i] = tr_y_yzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 735-756 components of targeted buffer : HH

    auto tr_y_xzzzz_xxxxx = pbuffer.data(idx_dip_hh + 735);

    auto tr_y_xzzzz_xxxxy = pbuffer.data(idx_dip_hh + 736);

    auto tr_y_xzzzz_xxxxz = pbuffer.data(idx_dip_hh + 737);

    auto tr_y_xzzzz_xxxyy = pbuffer.data(idx_dip_hh + 738);

    auto tr_y_xzzzz_xxxyz = pbuffer.data(idx_dip_hh + 739);

    auto tr_y_xzzzz_xxxzz = pbuffer.data(idx_dip_hh + 740);

    auto tr_y_xzzzz_xxyyy = pbuffer.data(idx_dip_hh + 741);

    auto tr_y_xzzzz_xxyyz = pbuffer.data(idx_dip_hh + 742);

    auto tr_y_xzzzz_xxyzz = pbuffer.data(idx_dip_hh + 743);

    auto tr_y_xzzzz_xxzzz = pbuffer.data(idx_dip_hh + 744);

    auto tr_y_xzzzz_xyyyy = pbuffer.data(idx_dip_hh + 745);

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

#pragma omp simd aligned(pa_x,                 \
                             tr_y_xzzzz_xxxxx, \
                             tr_y_xzzzz_xxxxy, \
                             tr_y_xzzzz_xxxxz, \
                             tr_y_xzzzz_xxxyy, \
                             tr_y_xzzzz_xxxyz, \
                             tr_y_xzzzz_xxxzz, \
                             tr_y_xzzzz_xxyyy, \
                             tr_y_xzzzz_xxyyz, \
                             tr_y_xzzzz_xxyzz, \
                             tr_y_xzzzz_xxzzz, \
                             tr_y_xzzzz_xyyyy, \
                             tr_y_xzzzz_xyyyz, \
                             tr_y_xzzzz_xyyzz, \
                             tr_y_xzzzz_xyzzz, \
                             tr_y_xzzzz_xzzzz, \
                             tr_y_xzzzz_yyyyy, \
                             tr_y_xzzzz_yyyyz, \
                             tr_y_xzzzz_yyyzz, \
                             tr_y_xzzzz_yyzzz, \
                             tr_y_xzzzz_yzzzz, \
                             tr_y_xzzzz_zzzzz, \
                             tr_y_zzzz_xxxx,   \
                             tr_y_zzzz_xxxxx,  \
                             tr_y_zzzz_xxxxy,  \
                             tr_y_zzzz_xxxxz,  \
                             tr_y_zzzz_xxxy,   \
                             tr_y_zzzz_xxxyy,  \
                             tr_y_zzzz_xxxyz,  \
                             tr_y_zzzz_xxxz,   \
                             tr_y_zzzz_xxxzz,  \
                             tr_y_zzzz_xxyy,   \
                             tr_y_zzzz_xxyyy,  \
                             tr_y_zzzz_xxyyz,  \
                             tr_y_zzzz_xxyz,   \
                             tr_y_zzzz_xxyzz,  \
                             tr_y_zzzz_xxzz,   \
                             tr_y_zzzz_xxzzz,  \
                             tr_y_zzzz_xyyy,   \
                             tr_y_zzzz_xyyyy,  \
                             tr_y_zzzz_xyyyz,  \
                             tr_y_zzzz_xyyz,   \
                             tr_y_zzzz_xyyzz,  \
                             tr_y_zzzz_xyzz,   \
                             tr_y_zzzz_xyzzz,  \
                             tr_y_zzzz_xzzz,   \
                             tr_y_zzzz_xzzzz,  \
                             tr_y_zzzz_yyyy,   \
                             tr_y_zzzz_yyyyy,  \
                             tr_y_zzzz_yyyyz,  \
                             tr_y_zzzz_yyyz,   \
                             tr_y_zzzz_yyyzz,  \
                             tr_y_zzzz_yyzz,   \
                             tr_y_zzzz_yyzzz,  \
                             tr_y_zzzz_yzzz,   \
                             tr_y_zzzz_yzzzz,  \
                             tr_y_zzzz_zzzz,   \
                             tr_y_zzzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzz_xxxxx[i] = 5.0 * tr_y_zzzz_xxxx[i] * fe_0 + tr_y_zzzz_xxxxx[i] * pa_x[i];

        tr_y_xzzzz_xxxxy[i] = 4.0 * tr_y_zzzz_xxxy[i] * fe_0 + tr_y_zzzz_xxxxy[i] * pa_x[i];

        tr_y_xzzzz_xxxxz[i] = 4.0 * tr_y_zzzz_xxxz[i] * fe_0 + tr_y_zzzz_xxxxz[i] * pa_x[i];

        tr_y_xzzzz_xxxyy[i] = 3.0 * tr_y_zzzz_xxyy[i] * fe_0 + tr_y_zzzz_xxxyy[i] * pa_x[i];

        tr_y_xzzzz_xxxyz[i] = 3.0 * tr_y_zzzz_xxyz[i] * fe_0 + tr_y_zzzz_xxxyz[i] * pa_x[i];

        tr_y_xzzzz_xxxzz[i] = 3.0 * tr_y_zzzz_xxzz[i] * fe_0 + tr_y_zzzz_xxxzz[i] * pa_x[i];

        tr_y_xzzzz_xxyyy[i] = 2.0 * tr_y_zzzz_xyyy[i] * fe_0 + tr_y_zzzz_xxyyy[i] * pa_x[i];

        tr_y_xzzzz_xxyyz[i] = 2.0 * tr_y_zzzz_xyyz[i] * fe_0 + tr_y_zzzz_xxyyz[i] * pa_x[i];

        tr_y_xzzzz_xxyzz[i] = 2.0 * tr_y_zzzz_xyzz[i] * fe_0 + tr_y_zzzz_xxyzz[i] * pa_x[i];

        tr_y_xzzzz_xxzzz[i] = 2.0 * tr_y_zzzz_xzzz[i] * fe_0 + tr_y_zzzz_xxzzz[i] * pa_x[i];

        tr_y_xzzzz_xyyyy[i] = tr_y_zzzz_yyyy[i] * fe_0 + tr_y_zzzz_xyyyy[i] * pa_x[i];

        tr_y_xzzzz_xyyyz[i] = tr_y_zzzz_yyyz[i] * fe_0 + tr_y_zzzz_xyyyz[i] * pa_x[i];

        tr_y_xzzzz_xyyzz[i] = tr_y_zzzz_yyzz[i] * fe_0 + tr_y_zzzz_xyyzz[i] * pa_x[i];

        tr_y_xzzzz_xyzzz[i] = tr_y_zzzz_yzzz[i] * fe_0 + tr_y_zzzz_xyzzz[i] * pa_x[i];

        tr_y_xzzzz_xzzzz[i] = tr_y_zzzz_zzzz[i] * fe_0 + tr_y_zzzz_xzzzz[i] * pa_x[i];

        tr_y_xzzzz_yyyyy[i] = tr_y_zzzz_yyyyy[i] * pa_x[i];

        tr_y_xzzzz_yyyyz[i] = tr_y_zzzz_yyyyz[i] * pa_x[i];

        tr_y_xzzzz_yyyzz[i] = tr_y_zzzz_yyyzz[i] * pa_x[i];

        tr_y_xzzzz_yyzzz[i] = tr_y_zzzz_yyzzz[i] * pa_x[i];

        tr_y_xzzzz_yzzzz[i] = tr_y_zzzz_yzzzz[i] * pa_x[i];

        tr_y_xzzzz_zzzzz[i] = tr_y_zzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 756-777 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             tr_y_yyy_xxxxx,   \
                             tr_y_yyy_xxxxy,   \
                             tr_y_yyy_xxxxz,   \
                             tr_y_yyy_xxxyy,   \
                             tr_y_yyy_xxxyz,   \
                             tr_y_yyy_xxxzz,   \
                             tr_y_yyy_xxyyy,   \
                             tr_y_yyy_xxyyz,   \
                             tr_y_yyy_xxyzz,   \
                             tr_y_yyy_xxzzz,   \
                             tr_y_yyy_xyyyy,   \
                             tr_y_yyy_xyyyz,   \
                             tr_y_yyy_xyyzz,   \
                             tr_y_yyy_xyzzz,   \
                             tr_y_yyy_xzzzz,   \
                             tr_y_yyy_yyyyy,   \
                             tr_y_yyy_yyyyz,   \
                             tr_y_yyy_yyyzz,   \
                             tr_y_yyy_yyzzz,   \
                             tr_y_yyy_yzzzz,   \
                             tr_y_yyy_zzzzz,   \
                             tr_y_yyyy_xxxx,   \
                             tr_y_yyyy_xxxxx,  \
                             tr_y_yyyy_xxxxy,  \
                             tr_y_yyyy_xxxxz,  \
                             tr_y_yyyy_xxxy,   \
                             tr_y_yyyy_xxxyy,  \
                             tr_y_yyyy_xxxyz,  \
                             tr_y_yyyy_xxxz,   \
                             tr_y_yyyy_xxxzz,  \
                             tr_y_yyyy_xxyy,   \
                             tr_y_yyyy_xxyyy,  \
                             tr_y_yyyy_xxyyz,  \
                             tr_y_yyyy_xxyz,   \
                             tr_y_yyyy_xxyzz,  \
                             tr_y_yyyy_xxzz,   \
                             tr_y_yyyy_xxzzz,  \
                             tr_y_yyyy_xyyy,   \
                             tr_y_yyyy_xyyyy,  \
                             tr_y_yyyy_xyyyz,  \
                             tr_y_yyyy_xyyz,   \
                             tr_y_yyyy_xyyzz,  \
                             tr_y_yyyy_xyzz,   \
                             tr_y_yyyy_xyzzz,  \
                             tr_y_yyyy_xzzz,   \
                             tr_y_yyyy_xzzzz,  \
                             tr_y_yyyy_yyyy,   \
                             tr_y_yyyy_yyyyy,  \
                             tr_y_yyyy_yyyyz,  \
                             tr_y_yyyy_yyyz,   \
                             tr_y_yyyy_yyyzz,  \
                             tr_y_yyyy_yyzz,   \
                             tr_y_yyyy_yyzzz,  \
                             tr_y_yyyy_yzzz,   \
                             tr_y_yyyy_yzzzz,  \
                             tr_y_yyyy_zzzz,   \
                             tr_y_yyyy_zzzzz,  \
                             tr_y_yyyyy_xxxxx, \
                             tr_y_yyyyy_xxxxy, \
                             tr_y_yyyyy_xxxxz, \
                             tr_y_yyyyy_xxxyy, \
                             tr_y_yyyyy_xxxyz, \
                             tr_y_yyyyy_xxxzz, \
                             tr_y_yyyyy_xxyyy, \
                             tr_y_yyyyy_xxyyz, \
                             tr_y_yyyyy_xxyzz, \
                             tr_y_yyyyy_xxzzz, \
                             tr_y_yyyyy_xyyyy, \
                             tr_y_yyyyy_xyyyz, \
                             tr_y_yyyyy_xyyzz, \
                             tr_y_yyyyy_xyzzz, \
                             tr_y_yyyyy_xzzzz, \
                             tr_y_yyyyy_yyyyy, \
                             tr_y_yyyyy_yyyyz, \
                             tr_y_yyyyy_yyyzz, \
                             tr_y_yyyyy_yyzzz, \
                             tr_y_yyyyy_yzzzz, \
                             tr_y_yyyyy_zzzzz, \
                             ts_yyyy_xxxxx,    \
                             ts_yyyy_xxxxy,    \
                             ts_yyyy_xxxxz,    \
                             ts_yyyy_xxxyy,    \
                             ts_yyyy_xxxyz,    \
                             ts_yyyy_xxxzz,    \
                             ts_yyyy_xxyyy,    \
                             ts_yyyy_xxyyz,    \
                             ts_yyyy_xxyzz,    \
                             ts_yyyy_xxzzz,    \
                             ts_yyyy_xyyyy,    \
                             ts_yyyy_xyyyz,    \
                             ts_yyyy_xyyzz,    \
                             ts_yyyy_xyzzz,    \
                             ts_yyyy_xzzzz,    \
                             ts_yyyy_yyyyy,    \
                             ts_yyyy_yyyyz,    \
                             ts_yyyy_yyyzz,    \
                             ts_yyyy_yyzzz,    \
                             ts_yyyy_yzzzz,    \
                             ts_yyyy_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyy_xxxxx[i] = 4.0 * tr_y_yyy_xxxxx[i] * fe_0 + ts_yyyy_xxxxx[i] * fe_0 + tr_y_yyyy_xxxxx[i] * pa_y[i];

        tr_y_yyyyy_xxxxy[i] = 4.0 * tr_y_yyy_xxxxy[i] * fe_0 + tr_y_yyyy_xxxx[i] * fe_0 + ts_yyyy_xxxxy[i] * fe_0 + tr_y_yyyy_xxxxy[i] * pa_y[i];

        tr_y_yyyyy_xxxxz[i] = 4.0 * tr_y_yyy_xxxxz[i] * fe_0 + ts_yyyy_xxxxz[i] * fe_0 + tr_y_yyyy_xxxxz[i] * pa_y[i];

        tr_y_yyyyy_xxxyy[i] =
            4.0 * tr_y_yyy_xxxyy[i] * fe_0 + 2.0 * tr_y_yyyy_xxxy[i] * fe_0 + ts_yyyy_xxxyy[i] * fe_0 + tr_y_yyyy_xxxyy[i] * pa_y[i];

        tr_y_yyyyy_xxxyz[i] = 4.0 * tr_y_yyy_xxxyz[i] * fe_0 + tr_y_yyyy_xxxz[i] * fe_0 + ts_yyyy_xxxyz[i] * fe_0 + tr_y_yyyy_xxxyz[i] * pa_y[i];

        tr_y_yyyyy_xxxzz[i] = 4.0 * tr_y_yyy_xxxzz[i] * fe_0 + ts_yyyy_xxxzz[i] * fe_0 + tr_y_yyyy_xxxzz[i] * pa_y[i];

        tr_y_yyyyy_xxyyy[i] =
            4.0 * tr_y_yyy_xxyyy[i] * fe_0 + 3.0 * tr_y_yyyy_xxyy[i] * fe_0 + ts_yyyy_xxyyy[i] * fe_0 + tr_y_yyyy_xxyyy[i] * pa_y[i];

        tr_y_yyyyy_xxyyz[i] =
            4.0 * tr_y_yyy_xxyyz[i] * fe_0 + 2.0 * tr_y_yyyy_xxyz[i] * fe_0 + ts_yyyy_xxyyz[i] * fe_0 + tr_y_yyyy_xxyyz[i] * pa_y[i];

        tr_y_yyyyy_xxyzz[i] = 4.0 * tr_y_yyy_xxyzz[i] * fe_0 + tr_y_yyyy_xxzz[i] * fe_0 + ts_yyyy_xxyzz[i] * fe_0 + tr_y_yyyy_xxyzz[i] * pa_y[i];

        tr_y_yyyyy_xxzzz[i] = 4.0 * tr_y_yyy_xxzzz[i] * fe_0 + ts_yyyy_xxzzz[i] * fe_0 + tr_y_yyyy_xxzzz[i] * pa_y[i];

        tr_y_yyyyy_xyyyy[i] =
            4.0 * tr_y_yyy_xyyyy[i] * fe_0 + 4.0 * tr_y_yyyy_xyyy[i] * fe_0 + ts_yyyy_xyyyy[i] * fe_0 + tr_y_yyyy_xyyyy[i] * pa_y[i];

        tr_y_yyyyy_xyyyz[i] =
            4.0 * tr_y_yyy_xyyyz[i] * fe_0 + 3.0 * tr_y_yyyy_xyyz[i] * fe_0 + ts_yyyy_xyyyz[i] * fe_0 + tr_y_yyyy_xyyyz[i] * pa_y[i];

        tr_y_yyyyy_xyyzz[i] =
            4.0 * tr_y_yyy_xyyzz[i] * fe_0 + 2.0 * tr_y_yyyy_xyzz[i] * fe_0 + ts_yyyy_xyyzz[i] * fe_0 + tr_y_yyyy_xyyzz[i] * pa_y[i];

        tr_y_yyyyy_xyzzz[i] = 4.0 * tr_y_yyy_xyzzz[i] * fe_0 + tr_y_yyyy_xzzz[i] * fe_0 + ts_yyyy_xyzzz[i] * fe_0 + tr_y_yyyy_xyzzz[i] * pa_y[i];

        tr_y_yyyyy_xzzzz[i] = 4.0 * tr_y_yyy_xzzzz[i] * fe_0 + ts_yyyy_xzzzz[i] * fe_0 + tr_y_yyyy_xzzzz[i] * pa_y[i];

        tr_y_yyyyy_yyyyy[i] =
            4.0 * tr_y_yyy_yyyyy[i] * fe_0 + 5.0 * tr_y_yyyy_yyyy[i] * fe_0 + ts_yyyy_yyyyy[i] * fe_0 + tr_y_yyyy_yyyyy[i] * pa_y[i];

        tr_y_yyyyy_yyyyz[i] =
            4.0 * tr_y_yyy_yyyyz[i] * fe_0 + 4.0 * tr_y_yyyy_yyyz[i] * fe_0 + ts_yyyy_yyyyz[i] * fe_0 + tr_y_yyyy_yyyyz[i] * pa_y[i];

        tr_y_yyyyy_yyyzz[i] =
            4.0 * tr_y_yyy_yyyzz[i] * fe_0 + 3.0 * tr_y_yyyy_yyzz[i] * fe_0 + ts_yyyy_yyyzz[i] * fe_0 + tr_y_yyyy_yyyzz[i] * pa_y[i];

        tr_y_yyyyy_yyzzz[i] =
            4.0 * tr_y_yyy_yyzzz[i] * fe_0 + 2.0 * tr_y_yyyy_yzzz[i] * fe_0 + ts_yyyy_yyzzz[i] * fe_0 + tr_y_yyyy_yyzzz[i] * pa_y[i];

        tr_y_yyyyy_yzzzz[i] = 4.0 * tr_y_yyy_yzzzz[i] * fe_0 + tr_y_yyyy_zzzz[i] * fe_0 + ts_yyyy_yzzzz[i] * fe_0 + tr_y_yyyy_yzzzz[i] * pa_y[i];

        tr_y_yyyyy_zzzzz[i] = 4.0 * tr_y_yyy_zzzzz[i] * fe_0 + ts_yyyy_zzzzz[i] * fe_0 + tr_y_yyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 777-798 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_z,                 \
                             tr_y_yyyy_xxxx,   \
                             tr_y_yyyy_xxxxx,  \
                             tr_y_yyyy_xxxxy,  \
                             tr_y_yyyy_xxxxz,  \
                             tr_y_yyyy_xxxy,   \
                             tr_y_yyyy_xxxyy,  \
                             tr_y_yyyy_xxxyz,  \
                             tr_y_yyyy_xxxz,   \
                             tr_y_yyyy_xxxzz,  \
                             tr_y_yyyy_xxyy,   \
                             tr_y_yyyy_xxyyy,  \
                             tr_y_yyyy_xxyyz,  \
                             tr_y_yyyy_xxyz,   \
                             tr_y_yyyy_xxyzz,  \
                             tr_y_yyyy_xxzz,   \
                             tr_y_yyyy_xxzzz,  \
                             tr_y_yyyy_xyyy,   \
                             tr_y_yyyy_xyyyy,  \
                             tr_y_yyyy_xyyyz,  \
                             tr_y_yyyy_xyyz,   \
                             tr_y_yyyy_xyyzz,  \
                             tr_y_yyyy_xyzz,   \
                             tr_y_yyyy_xyzzz,  \
                             tr_y_yyyy_xzzz,   \
                             tr_y_yyyy_xzzzz,  \
                             tr_y_yyyy_yyyy,   \
                             tr_y_yyyy_yyyyy,  \
                             tr_y_yyyy_yyyyz,  \
                             tr_y_yyyy_yyyz,   \
                             tr_y_yyyy_yyyzz,  \
                             tr_y_yyyy_yyzz,   \
                             tr_y_yyyy_yyzzz,  \
                             tr_y_yyyy_yzzz,   \
                             tr_y_yyyy_yzzzz,  \
                             tr_y_yyyy_zzzz,   \
                             tr_y_yyyy_zzzzz,  \
                             tr_y_yyyyz_xxxxx, \
                             tr_y_yyyyz_xxxxy, \
                             tr_y_yyyyz_xxxxz, \
                             tr_y_yyyyz_xxxyy, \
                             tr_y_yyyyz_xxxyz, \
                             tr_y_yyyyz_xxxzz, \
                             tr_y_yyyyz_xxyyy, \
                             tr_y_yyyyz_xxyyz, \
                             tr_y_yyyyz_xxyzz, \
                             tr_y_yyyyz_xxzzz, \
                             tr_y_yyyyz_xyyyy, \
                             tr_y_yyyyz_xyyyz, \
                             tr_y_yyyyz_xyyzz, \
                             tr_y_yyyyz_xyzzz, \
                             tr_y_yyyyz_xzzzz, \
                             tr_y_yyyyz_yyyyy, \
                             tr_y_yyyyz_yyyyz, \
                             tr_y_yyyyz_yyyzz, \
                             tr_y_yyyyz_yyzzz, \
                             tr_y_yyyyz_yzzzz, \
                             tr_y_yyyyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyz_xxxxx[i] = tr_y_yyyy_xxxxx[i] * pa_z[i];

        tr_y_yyyyz_xxxxy[i] = tr_y_yyyy_xxxxy[i] * pa_z[i];

        tr_y_yyyyz_xxxxz[i] = tr_y_yyyy_xxxx[i] * fe_0 + tr_y_yyyy_xxxxz[i] * pa_z[i];

        tr_y_yyyyz_xxxyy[i] = tr_y_yyyy_xxxyy[i] * pa_z[i];

        tr_y_yyyyz_xxxyz[i] = tr_y_yyyy_xxxy[i] * fe_0 + tr_y_yyyy_xxxyz[i] * pa_z[i];

        tr_y_yyyyz_xxxzz[i] = 2.0 * tr_y_yyyy_xxxz[i] * fe_0 + tr_y_yyyy_xxxzz[i] * pa_z[i];

        tr_y_yyyyz_xxyyy[i] = tr_y_yyyy_xxyyy[i] * pa_z[i];

        tr_y_yyyyz_xxyyz[i] = tr_y_yyyy_xxyy[i] * fe_0 + tr_y_yyyy_xxyyz[i] * pa_z[i];

        tr_y_yyyyz_xxyzz[i] = 2.0 * tr_y_yyyy_xxyz[i] * fe_0 + tr_y_yyyy_xxyzz[i] * pa_z[i];

        tr_y_yyyyz_xxzzz[i] = 3.0 * tr_y_yyyy_xxzz[i] * fe_0 + tr_y_yyyy_xxzzz[i] * pa_z[i];

        tr_y_yyyyz_xyyyy[i] = tr_y_yyyy_xyyyy[i] * pa_z[i];

        tr_y_yyyyz_xyyyz[i] = tr_y_yyyy_xyyy[i] * fe_0 + tr_y_yyyy_xyyyz[i] * pa_z[i];

        tr_y_yyyyz_xyyzz[i] = 2.0 * tr_y_yyyy_xyyz[i] * fe_0 + tr_y_yyyy_xyyzz[i] * pa_z[i];

        tr_y_yyyyz_xyzzz[i] = 3.0 * tr_y_yyyy_xyzz[i] * fe_0 + tr_y_yyyy_xyzzz[i] * pa_z[i];

        tr_y_yyyyz_xzzzz[i] = 4.0 * tr_y_yyyy_xzzz[i] * fe_0 + tr_y_yyyy_xzzzz[i] * pa_z[i];

        tr_y_yyyyz_yyyyy[i] = tr_y_yyyy_yyyyy[i] * pa_z[i];

        tr_y_yyyyz_yyyyz[i] = tr_y_yyyy_yyyy[i] * fe_0 + tr_y_yyyy_yyyyz[i] * pa_z[i];

        tr_y_yyyyz_yyyzz[i] = 2.0 * tr_y_yyyy_yyyz[i] * fe_0 + tr_y_yyyy_yyyzz[i] * pa_z[i];

        tr_y_yyyyz_yyzzz[i] = 3.0 * tr_y_yyyy_yyzz[i] * fe_0 + tr_y_yyyy_yyzzz[i] * pa_z[i];

        tr_y_yyyyz_yzzzz[i] = 4.0 * tr_y_yyyy_yzzz[i] * fe_0 + tr_y_yyyy_yzzzz[i] * pa_z[i];

        tr_y_yyyyz_zzzzz[i] = 5.0 * tr_y_yyyy_zzzz[i] * fe_0 + tr_y_yyyy_zzzzz[i] * pa_z[i];
    }

    // Set up 798-819 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yyy_xxxxx,   \
                             tr_y_yyy_xxxxy,   \
                             tr_y_yyy_xxxyy,   \
                             tr_y_yyy_xxxyz,   \
                             tr_y_yyy_xxyyy,   \
                             tr_y_yyy_xxyyz,   \
                             tr_y_yyy_xxyzz,   \
                             tr_y_yyy_xyyyy,   \
                             tr_y_yyy_xyyyz,   \
                             tr_y_yyy_xyyzz,   \
                             tr_y_yyy_xyzzz,   \
                             tr_y_yyy_yyyyy,   \
                             tr_y_yyy_yyyyz,   \
                             tr_y_yyy_yyyzz,   \
                             tr_y_yyy_yyzzz,   \
                             tr_y_yyy_yzzzz,   \
                             tr_y_yyyz_xxxxx,  \
                             tr_y_yyyz_xxxxy,  \
                             tr_y_yyyz_xxxy,   \
                             tr_y_yyyz_xxxyy,  \
                             tr_y_yyyz_xxxyz,  \
                             tr_y_yyyz_xxyy,   \
                             tr_y_yyyz_xxyyy,  \
                             tr_y_yyyz_xxyyz,  \
                             tr_y_yyyz_xxyz,   \
                             tr_y_yyyz_xxyzz,  \
                             tr_y_yyyz_xyyy,   \
                             tr_y_yyyz_xyyyy,  \
                             tr_y_yyyz_xyyyz,  \
                             tr_y_yyyz_xyyz,   \
                             tr_y_yyyz_xyyzz,  \
                             tr_y_yyyz_xyzz,   \
                             tr_y_yyyz_xyzzz,  \
                             tr_y_yyyz_yyyy,   \
                             tr_y_yyyz_yyyyy,  \
                             tr_y_yyyz_yyyyz,  \
                             tr_y_yyyz_yyyz,   \
                             tr_y_yyyz_yyyzz,  \
                             tr_y_yyyz_yyzz,   \
                             tr_y_yyyz_yyzzz,  \
                             tr_y_yyyz_yzzz,   \
                             tr_y_yyyz_yzzzz,  \
                             tr_y_yyyzz_xxxxx, \
                             tr_y_yyyzz_xxxxy, \
                             tr_y_yyyzz_xxxxz, \
                             tr_y_yyyzz_xxxyy, \
                             tr_y_yyyzz_xxxyz, \
                             tr_y_yyyzz_xxxzz, \
                             tr_y_yyyzz_xxyyy, \
                             tr_y_yyyzz_xxyyz, \
                             tr_y_yyyzz_xxyzz, \
                             tr_y_yyyzz_xxzzz, \
                             tr_y_yyyzz_xyyyy, \
                             tr_y_yyyzz_xyyyz, \
                             tr_y_yyyzz_xyyzz, \
                             tr_y_yyyzz_xyzzz, \
                             tr_y_yyyzz_xzzzz, \
                             tr_y_yyyzz_yyyyy, \
                             tr_y_yyyzz_yyyyz, \
                             tr_y_yyyzz_yyyzz, \
                             tr_y_yyyzz_yyzzz, \
                             tr_y_yyyzz_yzzzz, \
                             tr_y_yyyzz_zzzzz, \
                             tr_y_yyzz_xxxxz,  \
                             tr_y_yyzz_xxxzz,  \
                             tr_y_yyzz_xxzzz,  \
                             tr_y_yyzz_xzzzz,  \
                             tr_y_yyzz_zzzzz,  \
                             tr_y_yzz_xxxxz,   \
                             tr_y_yzz_xxxzz,   \
                             tr_y_yzz_xxzzz,   \
                             tr_y_yzz_xzzzz,   \
                             tr_y_yzz_zzzzz,   \
                             ts_yyzz_xxxxz,    \
                             ts_yyzz_xxxzz,    \
                             ts_yyzz_xxzzz,    \
                             ts_yyzz_xzzzz,    \
                             ts_yyzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzz_xxxxx[i] = tr_y_yyy_xxxxx[i] * fe_0 + tr_y_yyyz_xxxxx[i] * pa_z[i];

        tr_y_yyyzz_xxxxy[i] = tr_y_yyy_xxxxy[i] * fe_0 + tr_y_yyyz_xxxxy[i] * pa_z[i];

        tr_y_yyyzz_xxxxz[i] = 2.0 * tr_y_yzz_xxxxz[i] * fe_0 + ts_yyzz_xxxxz[i] * fe_0 + tr_y_yyzz_xxxxz[i] * pa_y[i];

        tr_y_yyyzz_xxxyy[i] = tr_y_yyy_xxxyy[i] * fe_0 + tr_y_yyyz_xxxyy[i] * pa_z[i];

        tr_y_yyyzz_xxxyz[i] = tr_y_yyy_xxxyz[i] * fe_0 + tr_y_yyyz_xxxy[i] * fe_0 + tr_y_yyyz_xxxyz[i] * pa_z[i];

        tr_y_yyyzz_xxxzz[i] = 2.0 * tr_y_yzz_xxxzz[i] * fe_0 + ts_yyzz_xxxzz[i] * fe_0 + tr_y_yyzz_xxxzz[i] * pa_y[i];

        tr_y_yyyzz_xxyyy[i] = tr_y_yyy_xxyyy[i] * fe_0 + tr_y_yyyz_xxyyy[i] * pa_z[i];

        tr_y_yyyzz_xxyyz[i] = tr_y_yyy_xxyyz[i] * fe_0 + tr_y_yyyz_xxyy[i] * fe_0 + tr_y_yyyz_xxyyz[i] * pa_z[i];

        tr_y_yyyzz_xxyzz[i] = tr_y_yyy_xxyzz[i] * fe_0 + 2.0 * tr_y_yyyz_xxyz[i] * fe_0 + tr_y_yyyz_xxyzz[i] * pa_z[i];

        tr_y_yyyzz_xxzzz[i] = 2.0 * tr_y_yzz_xxzzz[i] * fe_0 + ts_yyzz_xxzzz[i] * fe_0 + tr_y_yyzz_xxzzz[i] * pa_y[i];

        tr_y_yyyzz_xyyyy[i] = tr_y_yyy_xyyyy[i] * fe_0 + tr_y_yyyz_xyyyy[i] * pa_z[i];

        tr_y_yyyzz_xyyyz[i] = tr_y_yyy_xyyyz[i] * fe_0 + tr_y_yyyz_xyyy[i] * fe_0 + tr_y_yyyz_xyyyz[i] * pa_z[i];

        tr_y_yyyzz_xyyzz[i] = tr_y_yyy_xyyzz[i] * fe_0 + 2.0 * tr_y_yyyz_xyyz[i] * fe_0 + tr_y_yyyz_xyyzz[i] * pa_z[i];

        tr_y_yyyzz_xyzzz[i] = tr_y_yyy_xyzzz[i] * fe_0 + 3.0 * tr_y_yyyz_xyzz[i] * fe_0 + tr_y_yyyz_xyzzz[i] * pa_z[i];

        tr_y_yyyzz_xzzzz[i] = 2.0 * tr_y_yzz_xzzzz[i] * fe_0 + ts_yyzz_xzzzz[i] * fe_0 + tr_y_yyzz_xzzzz[i] * pa_y[i];

        tr_y_yyyzz_yyyyy[i] = tr_y_yyy_yyyyy[i] * fe_0 + tr_y_yyyz_yyyyy[i] * pa_z[i];

        tr_y_yyyzz_yyyyz[i] = tr_y_yyy_yyyyz[i] * fe_0 + tr_y_yyyz_yyyy[i] * fe_0 + tr_y_yyyz_yyyyz[i] * pa_z[i];

        tr_y_yyyzz_yyyzz[i] = tr_y_yyy_yyyzz[i] * fe_0 + 2.0 * tr_y_yyyz_yyyz[i] * fe_0 + tr_y_yyyz_yyyzz[i] * pa_z[i];

        tr_y_yyyzz_yyzzz[i] = tr_y_yyy_yyzzz[i] * fe_0 + 3.0 * tr_y_yyyz_yyzz[i] * fe_0 + tr_y_yyyz_yyzzz[i] * pa_z[i];

        tr_y_yyyzz_yzzzz[i] = tr_y_yyy_yzzzz[i] * fe_0 + 4.0 * tr_y_yyyz_yzzz[i] * fe_0 + tr_y_yyyz_yzzzz[i] * pa_z[i];

        tr_y_yyyzz_zzzzz[i] = 2.0 * tr_y_yzz_zzzzz[i] * fe_0 + ts_yyzz_zzzzz[i] * fe_0 + tr_y_yyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 819-840 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yyz_xxxxx,   \
                             tr_y_yyz_xxxxy,   \
                             tr_y_yyz_xxxyy,   \
                             tr_y_yyz_xxxyz,   \
                             tr_y_yyz_xxyyy,   \
                             tr_y_yyz_xxyyz,   \
                             tr_y_yyz_xxyzz,   \
                             tr_y_yyz_xyyyy,   \
                             tr_y_yyz_xyyyz,   \
                             tr_y_yyz_xyyzz,   \
                             tr_y_yyz_xyzzz,   \
                             tr_y_yyz_yyyyy,   \
                             tr_y_yyz_yyyyz,   \
                             tr_y_yyz_yyyzz,   \
                             tr_y_yyz_yyzzz,   \
                             tr_y_yyz_yzzzz,   \
                             tr_y_yyzz_xxxxx,  \
                             tr_y_yyzz_xxxxy,  \
                             tr_y_yyzz_xxxy,   \
                             tr_y_yyzz_xxxyy,  \
                             tr_y_yyzz_xxxyz,  \
                             tr_y_yyzz_xxyy,   \
                             tr_y_yyzz_xxyyy,  \
                             tr_y_yyzz_xxyyz,  \
                             tr_y_yyzz_xxyz,   \
                             tr_y_yyzz_xxyzz,  \
                             tr_y_yyzz_xyyy,   \
                             tr_y_yyzz_xyyyy,  \
                             tr_y_yyzz_xyyyz,  \
                             tr_y_yyzz_xyyz,   \
                             tr_y_yyzz_xyyzz,  \
                             tr_y_yyzz_xyzz,   \
                             tr_y_yyzz_xyzzz,  \
                             tr_y_yyzz_yyyy,   \
                             tr_y_yyzz_yyyyy,  \
                             tr_y_yyzz_yyyyz,  \
                             tr_y_yyzz_yyyz,   \
                             tr_y_yyzz_yyyzz,  \
                             tr_y_yyzz_yyzz,   \
                             tr_y_yyzz_yyzzz,  \
                             tr_y_yyzz_yzzz,   \
                             tr_y_yyzz_yzzzz,  \
                             tr_y_yyzzz_xxxxx, \
                             tr_y_yyzzz_xxxxy, \
                             tr_y_yyzzz_xxxxz, \
                             tr_y_yyzzz_xxxyy, \
                             tr_y_yyzzz_xxxyz, \
                             tr_y_yyzzz_xxxzz, \
                             tr_y_yyzzz_xxyyy, \
                             tr_y_yyzzz_xxyyz, \
                             tr_y_yyzzz_xxyzz, \
                             tr_y_yyzzz_xxzzz, \
                             tr_y_yyzzz_xyyyy, \
                             tr_y_yyzzz_xyyyz, \
                             tr_y_yyzzz_xyyzz, \
                             tr_y_yyzzz_xyzzz, \
                             tr_y_yyzzz_xzzzz, \
                             tr_y_yyzzz_yyyyy, \
                             tr_y_yyzzz_yyyyz, \
                             tr_y_yyzzz_yyyzz, \
                             tr_y_yyzzz_yyzzz, \
                             tr_y_yyzzz_yzzzz, \
                             tr_y_yyzzz_zzzzz, \
                             tr_y_yzzz_xxxxz,  \
                             tr_y_yzzz_xxxzz,  \
                             tr_y_yzzz_xxzzz,  \
                             tr_y_yzzz_xzzzz,  \
                             tr_y_yzzz_zzzzz,  \
                             tr_y_zzz_xxxxz,   \
                             tr_y_zzz_xxxzz,   \
                             tr_y_zzz_xxzzz,   \
                             tr_y_zzz_xzzzz,   \
                             tr_y_zzz_zzzzz,   \
                             ts_yzzz_xxxxz,    \
                             ts_yzzz_xxxzz,    \
                             ts_yzzz_xxzzz,    \
                             ts_yzzz_xzzzz,    \
                             ts_yzzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzz_xxxxx[i] = 2.0 * tr_y_yyz_xxxxx[i] * fe_0 + tr_y_yyzz_xxxxx[i] * pa_z[i];

        tr_y_yyzzz_xxxxy[i] = 2.0 * tr_y_yyz_xxxxy[i] * fe_0 + tr_y_yyzz_xxxxy[i] * pa_z[i];

        tr_y_yyzzz_xxxxz[i] = tr_y_zzz_xxxxz[i] * fe_0 + ts_yzzz_xxxxz[i] * fe_0 + tr_y_yzzz_xxxxz[i] * pa_y[i];

        tr_y_yyzzz_xxxyy[i] = 2.0 * tr_y_yyz_xxxyy[i] * fe_0 + tr_y_yyzz_xxxyy[i] * pa_z[i];

        tr_y_yyzzz_xxxyz[i] = 2.0 * tr_y_yyz_xxxyz[i] * fe_0 + tr_y_yyzz_xxxy[i] * fe_0 + tr_y_yyzz_xxxyz[i] * pa_z[i];

        tr_y_yyzzz_xxxzz[i] = tr_y_zzz_xxxzz[i] * fe_0 + ts_yzzz_xxxzz[i] * fe_0 + tr_y_yzzz_xxxzz[i] * pa_y[i];

        tr_y_yyzzz_xxyyy[i] = 2.0 * tr_y_yyz_xxyyy[i] * fe_0 + tr_y_yyzz_xxyyy[i] * pa_z[i];

        tr_y_yyzzz_xxyyz[i] = 2.0 * tr_y_yyz_xxyyz[i] * fe_0 + tr_y_yyzz_xxyy[i] * fe_0 + tr_y_yyzz_xxyyz[i] * pa_z[i];

        tr_y_yyzzz_xxyzz[i] = 2.0 * tr_y_yyz_xxyzz[i] * fe_0 + 2.0 * tr_y_yyzz_xxyz[i] * fe_0 + tr_y_yyzz_xxyzz[i] * pa_z[i];

        tr_y_yyzzz_xxzzz[i] = tr_y_zzz_xxzzz[i] * fe_0 + ts_yzzz_xxzzz[i] * fe_0 + tr_y_yzzz_xxzzz[i] * pa_y[i];

        tr_y_yyzzz_xyyyy[i] = 2.0 * tr_y_yyz_xyyyy[i] * fe_0 + tr_y_yyzz_xyyyy[i] * pa_z[i];

        tr_y_yyzzz_xyyyz[i] = 2.0 * tr_y_yyz_xyyyz[i] * fe_0 + tr_y_yyzz_xyyy[i] * fe_0 + tr_y_yyzz_xyyyz[i] * pa_z[i];

        tr_y_yyzzz_xyyzz[i] = 2.0 * tr_y_yyz_xyyzz[i] * fe_0 + 2.0 * tr_y_yyzz_xyyz[i] * fe_0 + tr_y_yyzz_xyyzz[i] * pa_z[i];

        tr_y_yyzzz_xyzzz[i] = 2.0 * tr_y_yyz_xyzzz[i] * fe_0 + 3.0 * tr_y_yyzz_xyzz[i] * fe_0 + tr_y_yyzz_xyzzz[i] * pa_z[i];

        tr_y_yyzzz_xzzzz[i] = tr_y_zzz_xzzzz[i] * fe_0 + ts_yzzz_xzzzz[i] * fe_0 + tr_y_yzzz_xzzzz[i] * pa_y[i];

        tr_y_yyzzz_yyyyy[i] = 2.0 * tr_y_yyz_yyyyy[i] * fe_0 + tr_y_yyzz_yyyyy[i] * pa_z[i];

        tr_y_yyzzz_yyyyz[i] = 2.0 * tr_y_yyz_yyyyz[i] * fe_0 + tr_y_yyzz_yyyy[i] * fe_0 + tr_y_yyzz_yyyyz[i] * pa_z[i];

        tr_y_yyzzz_yyyzz[i] = 2.0 * tr_y_yyz_yyyzz[i] * fe_0 + 2.0 * tr_y_yyzz_yyyz[i] * fe_0 + tr_y_yyzz_yyyzz[i] * pa_z[i];

        tr_y_yyzzz_yyzzz[i] = 2.0 * tr_y_yyz_yyzzz[i] * fe_0 + 3.0 * tr_y_yyzz_yyzz[i] * fe_0 + tr_y_yyzz_yyzzz[i] * pa_z[i];

        tr_y_yyzzz_yzzzz[i] = 2.0 * tr_y_yyz_yzzzz[i] * fe_0 + 4.0 * tr_y_yyzz_yzzz[i] * fe_0 + tr_y_yyzz_yzzzz[i] * pa_z[i];

        tr_y_yyzzz_zzzzz[i] = tr_y_zzz_zzzzz[i] * fe_0 + ts_yzzz_zzzzz[i] * fe_0 + tr_y_yzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 840-861 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_y_yzz_xxxxy,   \
                             tr_y_yzz_xxxyy,   \
                             tr_y_yzz_xxyyy,   \
                             tr_y_yzz_xyyyy,   \
                             tr_y_yzz_yyyyy,   \
                             tr_y_yzzz_xxxxy,  \
                             tr_y_yzzz_xxxyy,  \
                             tr_y_yzzz_xxyyy,  \
                             tr_y_yzzz_xyyyy,  \
                             tr_y_yzzz_yyyyy,  \
                             tr_y_yzzzz_xxxxx, \
                             tr_y_yzzzz_xxxxy, \
                             tr_y_yzzzz_xxxxz, \
                             tr_y_yzzzz_xxxyy, \
                             tr_y_yzzzz_xxxyz, \
                             tr_y_yzzzz_xxxzz, \
                             tr_y_yzzzz_xxyyy, \
                             tr_y_yzzzz_xxyyz, \
                             tr_y_yzzzz_xxyzz, \
                             tr_y_yzzzz_xxzzz, \
                             tr_y_yzzzz_xyyyy, \
                             tr_y_yzzzz_xyyyz, \
                             tr_y_yzzzz_xyyzz, \
                             tr_y_yzzzz_xyzzz, \
                             tr_y_yzzzz_xzzzz, \
                             tr_y_yzzzz_yyyyy, \
                             tr_y_yzzzz_yyyyz, \
                             tr_y_yzzzz_yyyzz, \
                             tr_y_yzzzz_yyzzz, \
                             tr_y_yzzzz_yzzzz, \
                             tr_y_yzzzz_zzzzz, \
                             tr_y_zzzz_xxxxx,  \
                             tr_y_zzzz_xxxxz,  \
                             tr_y_zzzz_xxxyz,  \
                             tr_y_zzzz_xxxz,   \
                             tr_y_zzzz_xxxzz,  \
                             tr_y_zzzz_xxyyz,  \
                             tr_y_zzzz_xxyz,   \
                             tr_y_zzzz_xxyzz,  \
                             tr_y_zzzz_xxzz,   \
                             tr_y_zzzz_xxzzz,  \
                             tr_y_zzzz_xyyyz,  \
                             tr_y_zzzz_xyyz,   \
                             tr_y_zzzz_xyyzz,  \
                             tr_y_zzzz_xyzz,   \
                             tr_y_zzzz_xyzzz,  \
                             tr_y_zzzz_xzzz,   \
                             tr_y_zzzz_xzzzz,  \
                             tr_y_zzzz_yyyyz,  \
                             tr_y_zzzz_yyyz,   \
                             tr_y_zzzz_yyyzz,  \
                             tr_y_zzzz_yyzz,   \
                             tr_y_zzzz_yyzzz,  \
                             tr_y_zzzz_yzzz,   \
                             tr_y_zzzz_yzzzz,  \
                             tr_y_zzzz_zzzz,   \
                             tr_y_zzzz_zzzzz,  \
                             ts_zzzz_xxxxx,    \
                             ts_zzzz_xxxxz,    \
                             ts_zzzz_xxxyz,    \
                             ts_zzzz_xxxzz,    \
                             ts_zzzz_xxyyz,    \
                             ts_zzzz_xxyzz,    \
                             ts_zzzz_xxzzz,    \
                             ts_zzzz_xyyyz,    \
                             ts_zzzz_xyyzz,    \
                             ts_zzzz_xyzzz,    \
                             ts_zzzz_xzzzz,    \
                             ts_zzzz_yyyyz,    \
                             ts_zzzz_yyyzz,    \
                             ts_zzzz_yyzzz,    \
                             ts_zzzz_yzzzz,    \
                             ts_zzzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzz_xxxxx[i] = ts_zzzz_xxxxx[i] * fe_0 + tr_y_zzzz_xxxxx[i] * pa_y[i];

        tr_y_yzzzz_xxxxy[i] = 3.0 * tr_y_yzz_xxxxy[i] * fe_0 + tr_y_yzzz_xxxxy[i] * pa_z[i];

        tr_y_yzzzz_xxxxz[i] = ts_zzzz_xxxxz[i] * fe_0 + tr_y_zzzz_xxxxz[i] * pa_y[i];

        tr_y_yzzzz_xxxyy[i] = 3.0 * tr_y_yzz_xxxyy[i] * fe_0 + tr_y_yzzz_xxxyy[i] * pa_z[i];

        tr_y_yzzzz_xxxyz[i] = tr_y_zzzz_xxxz[i] * fe_0 + ts_zzzz_xxxyz[i] * fe_0 + tr_y_zzzz_xxxyz[i] * pa_y[i];

        tr_y_yzzzz_xxxzz[i] = ts_zzzz_xxxzz[i] * fe_0 + tr_y_zzzz_xxxzz[i] * pa_y[i];

        tr_y_yzzzz_xxyyy[i] = 3.0 * tr_y_yzz_xxyyy[i] * fe_0 + tr_y_yzzz_xxyyy[i] * pa_z[i];

        tr_y_yzzzz_xxyyz[i] = 2.0 * tr_y_zzzz_xxyz[i] * fe_0 + ts_zzzz_xxyyz[i] * fe_0 + tr_y_zzzz_xxyyz[i] * pa_y[i];

        tr_y_yzzzz_xxyzz[i] = tr_y_zzzz_xxzz[i] * fe_0 + ts_zzzz_xxyzz[i] * fe_0 + tr_y_zzzz_xxyzz[i] * pa_y[i];

        tr_y_yzzzz_xxzzz[i] = ts_zzzz_xxzzz[i] * fe_0 + tr_y_zzzz_xxzzz[i] * pa_y[i];

        tr_y_yzzzz_xyyyy[i] = 3.0 * tr_y_yzz_xyyyy[i] * fe_0 + tr_y_yzzz_xyyyy[i] * pa_z[i];

        tr_y_yzzzz_xyyyz[i] = 3.0 * tr_y_zzzz_xyyz[i] * fe_0 + ts_zzzz_xyyyz[i] * fe_0 + tr_y_zzzz_xyyyz[i] * pa_y[i];

        tr_y_yzzzz_xyyzz[i] = 2.0 * tr_y_zzzz_xyzz[i] * fe_0 + ts_zzzz_xyyzz[i] * fe_0 + tr_y_zzzz_xyyzz[i] * pa_y[i];

        tr_y_yzzzz_xyzzz[i] = tr_y_zzzz_xzzz[i] * fe_0 + ts_zzzz_xyzzz[i] * fe_0 + tr_y_zzzz_xyzzz[i] * pa_y[i];

        tr_y_yzzzz_xzzzz[i] = ts_zzzz_xzzzz[i] * fe_0 + tr_y_zzzz_xzzzz[i] * pa_y[i];

        tr_y_yzzzz_yyyyy[i] = 3.0 * tr_y_yzz_yyyyy[i] * fe_0 + tr_y_yzzz_yyyyy[i] * pa_z[i];

        tr_y_yzzzz_yyyyz[i] = 4.0 * tr_y_zzzz_yyyz[i] * fe_0 + ts_zzzz_yyyyz[i] * fe_0 + tr_y_zzzz_yyyyz[i] * pa_y[i];

        tr_y_yzzzz_yyyzz[i] = 3.0 * tr_y_zzzz_yyzz[i] * fe_0 + ts_zzzz_yyyzz[i] * fe_0 + tr_y_zzzz_yyyzz[i] * pa_y[i];

        tr_y_yzzzz_yyzzz[i] = 2.0 * tr_y_zzzz_yzzz[i] * fe_0 + ts_zzzz_yyzzz[i] * fe_0 + tr_y_zzzz_yyzzz[i] * pa_y[i];

        tr_y_yzzzz_yzzzz[i] = tr_y_zzzz_zzzz[i] * fe_0 + ts_zzzz_yzzzz[i] * fe_0 + tr_y_zzzz_yzzzz[i] * pa_y[i];

        tr_y_yzzzz_zzzzz[i] = ts_zzzz_zzzzz[i] * fe_0 + tr_y_zzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 861-882 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_z,                 \
                             tr_y_zzz_xxxxx,   \
                             tr_y_zzz_xxxxy,   \
                             tr_y_zzz_xxxxz,   \
                             tr_y_zzz_xxxyy,   \
                             tr_y_zzz_xxxyz,   \
                             tr_y_zzz_xxxzz,   \
                             tr_y_zzz_xxyyy,   \
                             tr_y_zzz_xxyyz,   \
                             tr_y_zzz_xxyzz,   \
                             tr_y_zzz_xxzzz,   \
                             tr_y_zzz_xyyyy,   \
                             tr_y_zzz_xyyyz,   \
                             tr_y_zzz_xyyzz,   \
                             tr_y_zzz_xyzzz,   \
                             tr_y_zzz_xzzzz,   \
                             tr_y_zzz_yyyyy,   \
                             tr_y_zzz_yyyyz,   \
                             tr_y_zzz_yyyzz,   \
                             tr_y_zzz_yyzzz,   \
                             tr_y_zzz_yzzzz,   \
                             tr_y_zzz_zzzzz,   \
                             tr_y_zzzz_xxxx,   \
                             tr_y_zzzz_xxxxx,  \
                             tr_y_zzzz_xxxxy,  \
                             tr_y_zzzz_xxxxz,  \
                             tr_y_zzzz_xxxy,   \
                             tr_y_zzzz_xxxyy,  \
                             tr_y_zzzz_xxxyz,  \
                             tr_y_zzzz_xxxz,   \
                             tr_y_zzzz_xxxzz,  \
                             tr_y_zzzz_xxyy,   \
                             tr_y_zzzz_xxyyy,  \
                             tr_y_zzzz_xxyyz,  \
                             tr_y_zzzz_xxyz,   \
                             tr_y_zzzz_xxyzz,  \
                             tr_y_zzzz_xxzz,   \
                             tr_y_zzzz_xxzzz,  \
                             tr_y_zzzz_xyyy,   \
                             tr_y_zzzz_xyyyy,  \
                             tr_y_zzzz_xyyyz,  \
                             tr_y_zzzz_xyyz,   \
                             tr_y_zzzz_xyyzz,  \
                             tr_y_zzzz_xyzz,   \
                             tr_y_zzzz_xyzzz,  \
                             tr_y_zzzz_xzzz,   \
                             tr_y_zzzz_xzzzz,  \
                             tr_y_zzzz_yyyy,   \
                             tr_y_zzzz_yyyyy,  \
                             tr_y_zzzz_yyyyz,  \
                             tr_y_zzzz_yyyz,   \
                             tr_y_zzzz_yyyzz,  \
                             tr_y_zzzz_yyzz,   \
                             tr_y_zzzz_yyzzz,  \
                             tr_y_zzzz_yzzz,   \
                             tr_y_zzzz_yzzzz,  \
                             tr_y_zzzz_zzzz,   \
                             tr_y_zzzz_zzzzz,  \
                             tr_y_zzzzz_xxxxx, \
                             tr_y_zzzzz_xxxxy, \
                             tr_y_zzzzz_xxxxz, \
                             tr_y_zzzzz_xxxyy, \
                             tr_y_zzzzz_xxxyz, \
                             tr_y_zzzzz_xxxzz, \
                             tr_y_zzzzz_xxyyy, \
                             tr_y_zzzzz_xxyyz, \
                             tr_y_zzzzz_xxyzz, \
                             tr_y_zzzzz_xxzzz, \
                             tr_y_zzzzz_xyyyy, \
                             tr_y_zzzzz_xyyyz, \
                             tr_y_zzzzz_xyyzz, \
                             tr_y_zzzzz_xyzzz, \
                             tr_y_zzzzz_xzzzz, \
                             tr_y_zzzzz_yyyyy, \
                             tr_y_zzzzz_yyyyz, \
                             tr_y_zzzzz_yyyzz, \
                             tr_y_zzzzz_yyzzz, \
                             tr_y_zzzzz_yzzzz, \
                             tr_y_zzzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzz_xxxxx[i] = 4.0 * tr_y_zzz_xxxxx[i] * fe_0 + tr_y_zzzz_xxxxx[i] * pa_z[i];

        tr_y_zzzzz_xxxxy[i] = 4.0 * tr_y_zzz_xxxxy[i] * fe_0 + tr_y_zzzz_xxxxy[i] * pa_z[i];

        tr_y_zzzzz_xxxxz[i] = 4.0 * tr_y_zzz_xxxxz[i] * fe_0 + tr_y_zzzz_xxxx[i] * fe_0 + tr_y_zzzz_xxxxz[i] * pa_z[i];

        tr_y_zzzzz_xxxyy[i] = 4.0 * tr_y_zzz_xxxyy[i] * fe_0 + tr_y_zzzz_xxxyy[i] * pa_z[i];

        tr_y_zzzzz_xxxyz[i] = 4.0 * tr_y_zzz_xxxyz[i] * fe_0 + tr_y_zzzz_xxxy[i] * fe_0 + tr_y_zzzz_xxxyz[i] * pa_z[i];

        tr_y_zzzzz_xxxzz[i] = 4.0 * tr_y_zzz_xxxzz[i] * fe_0 + 2.0 * tr_y_zzzz_xxxz[i] * fe_0 + tr_y_zzzz_xxxzz[i] * pa_z[i];

        tr_y_zzzzz_xxyyy[i] = 4.0 * tr_y_zzz_xxyyy[i] * fe_0 + tr_y_zzzz_xxyyy[i] * pa_z[i];

        tr_y_zzzzz_xxyyz[i] = 4.0 * tr_y_zzz_xxyyz[i] * fe_0 + tr_y_zzzz_xxyy[i] * fe_0 + tr_y_zzzz_xxyyz[i] * pa_z[i];

        tr_y_zzzzz_xxyzz[i] = 4.0 * tr_y_zzz_xxyzz[i] * fe_0 + 2.0 * tr_y_zzzz_xxyz[i] * fe_0 + tr_y_zzzz_xxyzz[i] * pa_z[i];

        tr_y_zzzzz_xxzzz[i] = 4.0 * tr_y_zzz_xxzzz[i] * fe_0 + 3.0 * tr_y_zzzz_xxzz[i] * fe_0 + tr_y_zzzz_xxzzz[i] * pa_z[i];

        tr_y_zzzzz_xyyyy[i] = 4.0 * tr_y_zzz_xyyyy[i] * fe_0 + tr_y_zzzz_xyyyy[i] * pa_z[i];

        tr_y_zzzzz_xyyyz[i] = 4.0 * tr_y_zzz_xyyyz[i] * fe_0 + tr_y_zzzz_xyyy[i] * fe_0 + tr_y_zzzz_xyyyz[i] * pa_z[i];

        tr_y_zzzzz_xyyzz[i] = 4.0 * tr_y_zzz_xyyzz[i] * fe_0 + 2.0 * tr_y_zzzz_xyyz[i] * fe_0 + tr_y_zzzz_xyyzz[i] * pa_z[i];

        tr_y_zzzzz_xyzzz[i] = 4.0 * tr_y_zzz_xyzzz[i] * fe_0 + 3.0 * tr_y_zzzz_xyzz[i] * fe_0 + tr_y_zzzz_xyzzz[i] * pa_z[i];

        tr_y_zzzzz_xzzzz[i] = 4.0 * tr_y_zzz_xzzzz[i] * fe_0 + 4.0 * tr_y_zzzz_xzzz[i] * fe_0 + tr_y_zzzz_xzzzz[i] * pa_z[i];

        tr_y_zzzzz_yyyyy[i] = 4.0 * tr_y_zzz_yyyyy[i] * fe_0 + tr_y_zzzz_yyyyy[i] * pa_z[i];

        tr_y_zzzzz_yyyyz[i] = 4.0 * tr_y_zzz_yyyyz[i] * fe_0 + tr_y_zzzz_yyyy[i] * fe_0 + tr_y_zzzz_yyyyz[i] * pa_z[i];

        tr_y_zzzzz_yyyzz[i] = 4.0 * tr_y_zzz_yyyzz[i] * fe_0 + 2.0 * tr_y_zzzz_yyyz[i] * fe_0 + tr_y_zzzz_yyyzz[i] * pa_z[i];

        tr_y_zzzzz_yyzzz[i] = 4.0 * tr_y_zzz_yyzzz[i] * fe_0 + 3.0 * tr_y_zzzz_yyzz[i] * fe_0 + tr_y_zzzz_yyzzz[i] * pa_z[i];

        tr_y_zzzzz_yzzzz[i] = 4.0 * tr_y_zzz_yzzzz[i] * fe_0 + 4.0 * tr_y_zzzz_yzzz[i] * fe_0 + tr_y_zzzz_yzzzz[i] * pa_z[i];

        tr_y_zzzzz_zzzzz[i] = 4.0 * tr_y_zzz_zzzzz[i] * fe_0 + 5.0 * tr_y_zzzz_zzzz[i] * fe_0 + tr_y_zzzz_zzzzz[i] * pa_z[i];
    }

    // Set up 882-903 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxx_xxxxx,   \
                             tr_z_xxx_xxxxy,   \
                             tr_z_xxx_xxxxz,   \
                             tr_z_xxx_xxxyy,   \
                             tr_z_xxx_xxxyz,   \
                             tr_z_xxx_xxxzz,   \
                             tr_z_xxx_xxyyy,   \
                             tr_z_xxx_xxyyz,   \
                             tr_z_xxx_xxyzz,   \
                             tr_z_xxx_xxzzz,   \
                             tr_z_xxx_xyyyy,   \
                             tr_z_xxx_xyyyz,   \
                             tr_z_xxx_xyyzz,   \
                             tr_z_xxx_xyzzz,   \
                             tr_z_xxx_xzzzz,   \
                             tr_z_xxx_yyyyy,   \
                             tr_z_xxx_yyyyz,   \
                             tr_z_xxx_yyyzz,   \
                             tr_z_xxx_yyzzz,   \
                             tr_z_xxx_yzzzz,   \
                             tr_z_xxx_zzzzz,   \
                             tr_z_xxxx_xxxx,   \
                             tr_z_xxxx_xxxxx,  \
                             tr_z_xxxx_xxxxy,  \
                             tr_z_xxxx_xxxxz,  \
                             tr_z_xxxx_xxxy,   \
                             tr_z_xxxx_xxxyy,  \
                             tr_z_xxxx_xxxyz,  \
                             tr_z_xxxx_xxxz,   \
                             tr_z_xxxx_xxxzz,  \
                             tr_z_xxxx_xxyy,   \
                             tr_z_xxxx_xxyyy,  \
                             tr_z_xxxx_xxyyz,  \
                             tr_z_xxxx_xxyz,   \
                             tr_z_xxxx_xxyzz,  \
                             tr_z_xxxx_xxzz,   \
                             tr_z_xxxx_xxzzz,  \
                             tr_z_xxxx_xyyy,   \
                             tr_z_xxxx_xyyyy,  \
                             tr_z_xxxx_xyyyz,  \
                             tr_z_xxxx_xyyz,   \
                             tr_z_xxxx_xyyzz,  \
                             tr_z_xxxx_xyzz,   \
                             tr_z_xxxx_xyzzz,  \
                             tr_z_xxxx_xzzz,   \
                             tr_z_xxxx_xzzzz,  \
                             tr_z_xxxx_yyyy,   \
                             tr_z_xxxx_yyyyy,  \
                             tr_z_xxxx_yyyyz,  \
                             tr_z_xxxx_yyyz,   \
                             tr_z_xxxx_yyyzz,  \
                             tr_z_xxxx_yyzz,   \
                             tr_z_xxxx_yyzzz,  \
                             tr_z_xxxx_yzzz,   \
                             tr_z_xxxx_yzzzz,  \
                             tr_z_xxxx_zzzz,   \
                             tr_z_xxxx_zzzzz,  \
                             tr_z_xxxxx_xxxxx, \
                             tr_z_xxxxx_xxxxy, \
                             tr_z_xxxxx_xxxxz, \
                             tr_z_xxxxx_xxxyy, \
                             tr_z_xxxxx_xxxyz, \
                             tr_z_xxxxx_xxxzz, \
                             tr_z_xxxxx_xxyyy, \
                             tr_z_xxxxx_xxyyz, \
                             tr_z_xxxxx_xxyzz, \
                             tr_z_xxxxx_xxzzz, \
                             tr_z_xxxxx_xyyyy, \
                             tr_z_xxxxx_xyyyz, \
                             tr_z_xxxxx_xyyzz, \
                             tr_z_xxxxx_xyzzz, \
                             tr_z_xxxxx_xzzzz, \
                             tr_z_xxxxx_yyyyy, \
                             tr_z_xxxxx_yyyyz, \
                             tr_z_xxxxx_yyyzz, \
                             tr_z_xxxxx_yyzzz, \
                             tr_z_xxxxx_yzzzz, \
                             tr_z_xxxxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxx_xxxxx[i] = 4.0 * tr_z_xxx_xxxxx[i] * fe_0 + 5.0 * tr_z_xxxx_xxxx[i] * fe_0 + tr_z_xxxx_xxxxx[i] * pa_x[i];

        tr_z_xxxxx_xxxxy[i] = 4.0 * tr_z_xxx_xxxxy[i] * fe_0 + 4.0 * tr_z_xxxx_xxxy[i] * fe_0 + tr_z_xxxx_xxxxy[i] * pa_x[i];

        tr_z_xxxxx_xxxxz[i] = 4.0 * tr_z_xxx_xxxxz[i] * fe_0 + 4.0 * tr_z_xxxx_xxxz[i] * fe_0 + tr_z_xxxx_xxxxz[i] * pa_x[i];

        tr_z_xxxxx_xxxyy[i] = 4.0 * tr_z_xxx_xxxyy[i] * fe_0 + 3.0 * tr_z_xxxx_xxyy[i] * fe_0 + tr_z_xxxx_xxxyy[i] * pa_x[i];

        tr_z_xxxxx_xxxyz[i] = 4.0 * tr_z_xxx_xxxyz[i] * fe_0 + 3.0 * tr_z_xxxx_xxyz[i] * fe_0 + tr_z_xxxx_xxxyz[i] * pa_x[i];

        tr_z_xxxxx_xxxzz[i] = 4.0 * tr_z_xxx_xxxzz[i] * fe_0 + 3.0 * tr_z_xxxx_xxzz[i] * fe_0 + tr_z_xxxx_xxxzz[i] * pa_x[i];

        tr_z_xxxxx_xxyyy[i] = 4.0 * tr_z_xxx_xxyyy[i] * fe_0 + 2.0 * tr_z_xxxx_xyyy[i] * fe_0 + tr_z_xxxx_xxyyy[i] * pa_x[i];

        tr_z_xxxxx_xxyyz[i] = 4.0 * tr_z_xxx_xxyyz[i] * fe_0 + 2.0 * tr_z_xxxx_xyyz[i] * fe_0 + tr_z_xxxx_xxyyz[i] * pa_x[i];

        tr_z_xxxxx_xxyzz[i] = 4.0 * tr_z_xxx_xxyzz[i] * fe_0 + 2.0 * tr_z_xxxx_xyzz[i] * fe_0 + tr_z_xxxx_xxyzz[i] * pa_x[i];

        tr_z_xxxxx_xxzzz[i] = 4.0 * tr_z_xxx_xxzzz[i] * fe_0 + 2.0 * tr_z_xxxx_xzzz[i] * fe_0 + tr_z_xxxx_xxzzz[i] * pa_x[i];

        tr_z_xxxxx_xyyyy[i] = 4.0 * tr_z_xxx_xyyyy[i] * fe_0 + tr_z_xxxx_yyyy[i] * fe_0 + tr_z_xxxx_xyyyy[i] * pa_x[i];

        tr_z_xxxxx_xyyyz[i] = 4.0 * tr_z_xxx_xyyyz[i] * fe_0 + tr_z_xxxx_yyyz[i] * fe_0 + tr_z_xxxx_xyyyz[i] * pa_x[i];

        tr_z_xxxxx_xyyzz[i] = 4.0 * tr_z_xxx_xyyzz[i] * fe_0 + tr_z_xxxx_yyzz[i] * fe_0 + tr_z_xxxx_xyyzz[i] * pa_x[i];

        tr_z_xxxxx_xyzzz[i] = 4.0 * tr_z_xxx_xyzzz[i] * fe_0 + tr_z_xxxx_yzzz[i] * fe_0 + tr_z_xxxx_xyzzz[i] * pa_x[i];

        tr_z_xxxxx_xzzzz[i] = 4.0 * tr_z_xxx_xzzzz[i] * fe_0 + tr_z_xxxx_zzzz[i] * fe_0 + tr_z_xxxx_xzzzz[i] * pa_x[i];

        tr_z_xxxxx_yyyyy[i] = 4.0 * tr_z_xxx_yyyyy[i] * fe_0 + tr_z_xxxx_yyyyy[i] * pa_x[i];

        tr_z_xxxxx_yyyyz[i] = 4.0 * tr_z_xxx_yyyyz[i] * fe_0 + tr_z_xxxx_yyyyz[i] * pa_x[i];

        tr_z_xxxxx_yyyzz[i] = 4.0 * tr_z_xxx_yyyzz[i] * fe_0 + tr_z_xxxx_yyyzz[i] * pa_x[i];

        tr_z_xxxxx_yyzzz[i] = 4.0 * tr_z_xxx_yyzzz[i] * fe_0 + tr_z_xxxx_yyzzz[i] * pa_x[i];

        tr_z_xxxxx_yzzzz[i] = 4.0 * tr_z_xxx_yzzzz[i] * fe_0 + tr_z_xxxx_yzzzz[i] * pa_x[i];

        tr_z_xxxxx_zzzzz[i] = 4.0 * tr_z_xxx_zzzzz[i] * fe_0 + tr_z_xxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 903-924 components of targeted buffer : HH

    auto tr_z_xxxxy_xxxxx = pbuffer.data(idx_dip_hh + 903);

    auto tr_z_xxxxy_xxxxy = pbuffer.data(idx_dip_hh + 904);

    auto tr_z_xxxxy_xxxxz = pbuffer.data(idx_dip_hh + 905);

    auto tr_z_xxxxy_xxxyy = pbuffer.data(idx_dip_hh + 906);

    auto tr_z_xxxxy_xxxyz = pbuffer.data(idx_dip_hh + 907);

    auto tr_z_xxxxy_xxxzz = pbuffer.data(idx_dip_hh + 908);

    auto tr_z_xxxxy_xxyyy = pbuffer.data(idx_dip_hh + 909);

    auto tr_z_xxxxy_xxyyz = pbuffer.data(idx_dip_hh + 910);

    auto tr_z_xxxxy_xxyzz = pbuffer.data(idx_dip_hh + 911);

    auto tr_z_xxxxy_xxzzz = pbuffer.data(idx_dip_hh + 912);

    auto tr_z_xxxxy_xyyyy = pbuffer.data(idx_dip_hh + 913);

    auto tr_z_xxxxy_xyyyz = pbuffer.data(idx_dip_hh + 914);

    auto tr_z_xxxxy_xyyzz = pbuffer.data(idx_dip_hh + 915);

    auto tr_z_xxxxy_xyzzz = pbuffer.data(idx_dip_hh + 916);

    auto tr_z_xxxxy_xzzzz = pbuffer.data(idx_dip_hh + 917);

    auto tr_z_xxxxy_yyyyy = pbuffer.data(idx_dip_hh + 918);

    auto tr_z_xxxxy_yyyyz = pbuffer.data(idx_dip_hh + 919);

    auto tr_z_xxxxy_yyyzz = pbuffer.data(idx_dip_hh + 920);

    auto tr_z_xxxxy_yyzzz = pbuffer.data(idx_dip_hh + 921);

    auto tr_z_xxxxy_yzzzz = pbuffer.data(idx_dip_hh + 922);

    auto tr_z_xxxxy_zzzzz = pbuffer.data(idx_dip_hh + 923);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxx_xxxx,   \
                             tr_z_xxxx_xxxxx,  \
                             tr_z_xxxx_xxxxy,  \
                             tr_z_xxxx_xxxxz,  \
                             tr_z_xxxx_xxxy,   \
                             tr_z_xxxx_xxxyy,  \
                             tr_z_xxxx_xxxyz,  \
                             tr_z_xxxx_xxxz,   \
                             tr_z_xxxx_xxxzz,  \
                             tr_z_xxxx_xxyy,   \
                             tr_z_xxxx_xxyyy,  \
                             tr_z_xxxx_xxyyz,  \
                             tr_z_xxxx_xxyz,   \
                             tr_z_xxxx_xxyzz,  \
                             tr_z_xxxx_xxzz,   \
                             tr_z_xxxx_xxzzz,  \
                             tr_z_xxxx_xyyy,   \
                             tr_z_xxxx_xyyyy,  \
                             tr_z_xxxx_xyyyz,  \
                             tr_z_xxxx_xyyz,   \
                             tr_z_xxxx_xyyzz,  \
                             tr_z_xxxx_xyzz,   \
                             tr_z_xxxx_xyzzz,  \
                             tr_z_xxxx_xzzz,   \
                             tr_z_xxxx_xzzzz,  \
                             tr_z_xxxx_zzzzz,  \
                             tr_z_xxxxy_xxxxx, \
                             tr_z_xxxxy_xxxxy, \
                             tr_z_xxxxy_xxxxz, \
                             tr_z_xxxxy_xxxyy, \
                             tr_z_xxxxy_xxxyz, \
                             tr_z_xxxxy_xxxzz, \
                             tr_z_xxxxy_xxyyy, \
                             tr_z_xxxxy_xxyyz, \
                             tr_z_xxxxy_xxyzz, \
                             tr_z_xxxxy_xxzzz, \
                             tr_z_xxxxy_xyyyy, \
                             tr_z_xxxxy_xyyyz, \
                             tr_z_xxxxy_xyyzz, \
                             tr_z_xxxxy_xyzzz, \
                             tr_z_xxxxy_xzzzz, \
                             tr_z_xxxxy_yyyyy, \
                             tr_z_xxxxy_yyyyz, \
                             tr_z_xxxxy_yyyzz, \
                             tr_z_xxxxy_yyzzz, \
                             tr_z_xxxxy_yzzzz, \
                             tr_z_xxxxy_zzzzz, \
                             tr_z_xxxy_yyyyy,  \
                             tr_z_xxxy_yyyyz,  \
                             tr_z_xxxy_yyyzz,  \
                             tr_z_xxxy_yyzzz,  \
                             tr_z_xxxy_yzzzz,  \
                             tr_z_xxy_yyyyy,   \
                             tr_z_xxy_yyyyz,   \
                             tr_z_xxy_yyyzz,   \
                             tr_z_xxy_yyzzz,   \
                             tr_z_xxy_yzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxy_xxxxx[i] = tr_z_xxxx_xxxxx[i] * pa_y[i];

        tr_z_xxxxy_xxxxy[i] = tr_z_xxxx_xxxx[i] * fe_0 + tr_z_xxxx_xxxxy[i] * pa_y[i];

        tr_z_xxxxy_xxxxz[i] = tr_z_xxxx_xxxxz[i] * pa_y[i];

        tr_z_xxxxy_xxxyy[i] = 2.0 * tr_z_xxxx_xxxy[i] * fe_0 + tr_z_xxxx_xxxyy[i] * pa_y[i];

        tr_z_xxxxy_xxxyz[i] = tr_z_xxxx_xxxz[i] * fe_0 + tr_z_xxxx_xxxyz[i] * pa_y[i];

        tr_z_xxxxy_xxxzz[i] = tr_z_xxxx_xxxzz[i] * pa_y[i];

        tr_z_xxxxy_xxyyy[i] = 3.0 * tr_z_xxxx_xxyy[i] * fe_0 + tr_z_xxxx_xxyyy[i] * pa_y[i];

        tr_z_xxxxy_xxyyz[i] = 2.0 * tr_z_xxxx_xxyz[i] * fe_0 + tr_z_xxxx_xxyyz[i] * pa_y[i];

        tr_z_xxxxy_xxyzz[i] = tr_z_xxxx_xxzz[i] * fe_0 + tr_z_xxxx_xxyzz[i] * pa_y[i];

        tr_z_xxxxy_xxzzz[i] = tr_z_xxxx_xxzzz[i] * pa_y[i];

        tr_z_xxxxy_xyyyy[i] = 4.0 * tr_z_xxxx_xyyy[i] * fe_0 + tr_z_xxxx_xyyyy[i] * pa_y[i];

        tr_z_xxxxy_xyyyz[i] = 3.0 * tr_z_xxxx_xyyz[i] * fe_0 + tr_z_xxxx_xyyyz[i] * pa_y[i];

        tr_z_xxxxy_xyyzz[i] = 2.0 * tr_z_xxxx_xyzz[i] * fe_0 + tr_z_xxxx_xyyzz[i] * pa_y[i];

        tr_z_xxxxy_xyzzz[i] = tr_z_xxxx_xzzz[i] * fe_0 + tr_z_xxxx_xyzzz[i] * pa_y[i];

        tr_z_xxxxy_xzzzz[i] = tr_z_xxxx_xzzzz[i] * pa_y[i];

        tr_z_xxxxy_yyyyy[i] = 3.0 * tr_z_xxy_yyyyy[i] * fe_0 + tr_z_xxxy_yyyyy[i] * pa_x[i];

        tr_z_xxxxy_yyyyz[i] = 3.0 * tr_z_xxy_yyyyz[i] * fe_0 + tr_z_xxxy_yyyyz[i] * pa_x[i];

        tr_z_xxxxy_yyyzz[i] = 3.0 * tr_z_xxy_yyyzz[i] * fe_0 + tr_z_xxxy_yyyzz[i] * pa_x[i];

        tr_z_xxxxy_yyzzz[i] = 3.0 * tr_z_xxy_yyzzz[i] * fe_0 + tr_z_xxxy_yyzzz[i] * pa_x[i];

        tr_z_xxxxy_yzzzz[i] = 3.0 * tr_z_xxy_yzzzz[i] * fe_0 + tr_z_xxxy_yzzzz[i] * pa_x[i];

        tr_z_xxxxy_zzzzz[i] = tr_z_xxxx_zzzzz[i] * pa_y[i];
    }

    // Set up 924-945 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             tr_z_xxxx_xxxxx,  \
                             tr_z_xxxx_xxxxy,  \
                             tr_z_xxxx_xxxyy,  \
                             tr_z_xxxx_xxyyy,  \
                             tr_z_xxxx_xyyyy,  \
                             tr_z_xxxxz_xxxxx, \
                             tr_z_xxxxz_xxxxy, \
                             tr_z_xxxxz_xxxxz, \
                             tr_z_xxxxz_xxxyy, \
                             tr_z_xxxxz_xxxyz, \
                             tr_z_xxxxz_xxxzz, \
                             tr_z_xxxxz_xxyyy, \
                             tr_z_xxxxz_xxyyz, \
                             tr_z_xxxxz_xxyzz, \
                             tr_z_xxxxz_xxzzz, \
                             tr_z_xxxxz_xyyyy, \
                             tr_z_xxxxz_xyyyz, \
                             tr_z_xxxxz_xyyzz, \
                             tr_z_xxxxz_xyzzz, \
                             tr_z_xxxxz_xzzzz, \
                             tr_z_xxxxz_yyyyy, \
                             tr_z_xxxxz_yyyyz, \
                             tr_z_xxxxz_yyyzz, \
                             tr_z_xxxxz_yyzzz, \
                             tr_z_xxxxz_yzzzz, \
                             tr_z_xxxxz_zzzzz, \
                             tr_z_xxxz_xxxxz,  \
                             tr_z_xxxz_xxxyz,  \
                             tr_z_xxxz_xxxz,   \
                             tr_z_xxxz_xxxzz,  \
                             tr_z_xxxz_xxyyz,  \
                             tr_z_xxxz_xxyz,   \
                             tr_z_xxxz_xxyzz,  \
                             tr_z_xxxz_xxzz,   \
                             tr_z_xxxz_xxzzz,  \
                             tr_z_xxxz_xyyyz,  \
                             tr_z_xxxz_xyyz,   \
                             tr_z_xxxz_xyyzz,  \
                             tr_z_xxxz_xyzz,   \
                             tr_z_xxxz_xyzzz,  \
                             tr_z_xxxz_xzzz,   \
                             tr_z_xxxz_xzzzz,  \
                             tr_z_xxxz_yyyyy,  \
                             tr_z_xxxz_yyyyz,  \
                             tr_z_xxxz_yyyz,   \
                             tr_z_xxxz_yyyzz,  \
                             tr_z_xxxz_yyzz,   \
                             tr_z_xxxz_yyzzz,  \
                             tr_z_xxxz_yzzz,   \
                             tr_z_xxxz_yzzzz,  \
                             tr_z_xxxz_zzzz,   \
                             tr_z_xxxz_zzzzz,  \
                             tr_z_xxz_xxxxz,   \
                             tr_z_xxz_xxxyz,   \
                             tr_z_xxz_xxxzz,   \
                             tr_z_xxz_xxyyz,   \
                             tr_z_xxz_xxyzz,   \
                             tr_z_xxz_xxzzz,   \
                             tr_z_xxz_xyyyz,   \
                             tr_z_xxz_xyyzz,   \
                             tr_z_xxz_xyzzz,   \
                             tr_z_xxz_xzzzz,   \
                             tr_z_xxz_yyyyy,   \
                             tr_z_xxz_yyyyz,   \
                             tr_z_xxz_yyyzz,   \
                             tr_z_xxz_yyzzz,   \
                             tr_z_xxz_yzzzz,   \
                             tr_z_xxz_zzzzz,   \
                             ts_xxxx_xxxxx,    \
                             ts_xxxx_xxxxy,    \
                             ts_xxxx_xxxyy,    \
                             ts_xxxx_xxyyy,    \
                             ts_xxxx_xyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxz_xxxxx[i] = ts_xxxx_xxxxx[i] * fe_0 + tr_z_xxxx_xxxxx[i] * pa_z[i];

        tr_z_xxxxz_xxxxy[i] = ts_xxxx_xxxxy[i] * fe_0 + tr_z_xxxx_xxxxy[i] * pa_z[i];

        tr_z_xxxxz_xxxxz[i] = 3.0 * tr_z_xxz_xxxxz[i] * fe_0 + 4.0 * tr_z_xxxz_xxxz[i] * fe_0 + tr_z_xxxz_xxxxz[i] * pa_x[i];

        tr_z_xxxxz_xxxyy[i] = ts_xxxx_xxxyy[i] * fe_0 + tr_z_xxxx_xxxyy[i] * pa_z[i];

        tr_z_xxxxz_xxxyz[i] = 3.0 * tr_z_xxz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxxz_xxyz[i] * fe_0 + tr_z_xxxz_xxxyz[i] * pa_x[i];

        tr_z_xxxxz_xxxzz[i] = 3.0 * tr_z_xxz_xxxzz[i] * fe_0 + 3.0 * tr_z_xxxz_xxzz[i] * fe_0 + tr_z_xxxz_xxxzz[i] * pa_x[i];

        tr_z_xxxxz_xxyyy[i] = ts_xxxx_xxyyy[i] * fe_0 + tr_z_xxxx_xxyyy[i] * pa_z[i];

        tr_z_xxxxz_xxyyz[i] = 3.0 * tr_z_xxz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxxz_xyyz[i] * fe_0 + tr_z_xxxz_xxyyz[i] * pa_x[i];

        tr_z_xxxxz_xxyzz[i] = 3.0 * tr_z_xxz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxxz_xyzz[i] * fe_0 + tr_z_xxxz_xxyzz[i] * pa_x[i];

        tr_z_xxxxz_xxzzz[i] = 3.0 * tr_z_xxz_xxzzz[i] * fe_0 + 2.0 * tr_z_xxxz_xzzz[i] * fe_0 + tr_z_xxxz_xxzzz[i] * pa_x[i];

        tr_z_xxxxz_xyyyy[i] = ts_xxxx_xyyyy[i] * fe_0 + tr_z_xxxx_xyyyy[i] * pa_z[i];

        tr_z_xxxxz_xyyyz[i] = 3.0 * tr_z_xxz_xyyyz[i] * fe_0 + tr_z_xxxz_yyyz[i] * fe_0 + tr_z_xxxz_xyyyz[i] * pa_x[i];

        tr_z_xxxxz_xyyzz[i] = 3.0 * tr_z_xxz_xyyzz[i] * fe_0 + tr_z_xxxz_yyzz[i] * fe_0 + tr_z_xxxz_xyyzz[i] * pa_x[i];

        tr_z_xxxxz_xyzzz[i] = 3.0 * tr_z_xxz_xyzzz[i] * fe_0 + tr_z_xxxz_yzzz[i] * fe_0 + tr_z_xxxz_xyzzz[i] * pa_x[i];

        tr_z_xxxxz_xzzzz[i] = 3.0 * tr_z_xxz_xzzzz[i] * fe_0 + tr_z_xxxz_zzzz[i] * fe_0 + tr_z_xxxz_xzzzz[i] * pa_x[i];

        tr_z_xxxxz_yyyyy[i] = 3.0 * tr_z_xxz_yyyyy[i] * fe_0 + tr_z_xxxz_yyyyy[i] * pa_x[i];

        tr_z_xxxxz_yyyyz[i] = 3.0 * tr_z_xxz_yyyyz[i] * fe_0 + tr_z_xxxz_yyyyz[i] * pa_x[i];

        tr_z_xxxxz_yyyzz[i] = 3.0 * tr_z_xxz_yyyzz[i] * fe_0 + tr_z_xxxz_yyyzz[i] * pa_x[i];

        tr_z_xxxxz_yyzzz[i] = 3.0 * tr_z_xxz_yyzzz[i] * fe_0 + tr_z_xxxz_yyzzz[i] * pa_x[i];

        tr_z_xxxxz_yzzzz[i] = 3.0 * tr_z_xxz_yzzzz[i] * fe_0 + tr_z_xxxz_yzzzz[i] * pa_x[i];

        tr_z_xxxxz_zzzzz[i] = 3.0 * tr_z_xxz_zzzzz[i] * fe_0 + tr_z_xxxz_zzzzz[i] * pa_x[i];
    }

    // Set up 945-966 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxx_xxxxx,   \
                             tr_z_xxx_xxxxz,   \
                             tr_z_xxx_xxxzz,   \
                             tr_z_xxx_xxzzz,   \
                             tr_z_xxx_xzzzz,   \
                             tr_z_xxxy_xxxxx,  \
                             tr_z_xxxy_xxxxz,  \
                             tr_z_xxxy_xxxzz,  \
                             tr_z_xxxy_xxzzz,  \
                             tr_z_xxxy_xzzzz,  \
                             tr_z_xxxyy_xxxxx, \
                             tr_z_xxxyy_xxxxy, \
                             tr_z_xxxyy_xxxxz, \
                             tr_z_xxxyy_xxxyy, \
                             tr_z_xxxyy_xxxyz, \
                             tr_z_xxxyy_xxxzz, \
                             tr_z_xxxyy_xxyyy, \
                             tr_z_xxxyy_xxyyz, \
                             tr_z_xxxyy_xxyzz, \
                             tr_z_xxxyy_xxzzz, \
                             tr_z_xxxyy_xyyyy, \
                             tr_z_xxxyy_xyyyz, \
                             tr_z_xxxyy_xyyzz, \
                             tr_z_xxxyy_xyzzz, \
                             tr_z_xxxyy_xzzzz, \
                             tr_z_xxxyy_yyyyy, \
                             tr_z_xxxyy_yyyyz, \
                             tr_z_xxxyy_yyyzz, \
                             tr_z_xxxyy_yyzzz, \
                             tr_z_xxxyy_yzzzz, \
                             tr_z_xxxyy_zzzzz, \
                             tr_z_xxyy_xxxxy,  \
                             tr_z_xxyy_xxxy,   \
                             tr_z_xxyy_xxxyy,  \
                             tr_z_xxyy_xxxyz,  \
                             tr_z_xxyy_xxyy,   \
                             tr_z_xxyy_xxyyy,  \
                             tr_z_xxyy_xxyyz,  \
                             tr_z_xxyy_xxyz,   \
                             tr_z_xxyy_xxyzz,  \
                             tr_z_xxyy_xyyy,   \
                             tr_z_xxyy_xyyyy,  \
                             tr_z_xxyy_xyyyz,  \
                             tr_z_xxyy_xyyz,   \
                             tr_z_xxyy_xyyzz,  \
                             tr_z_xxyy_xyzz,   \
                             tr_z_xxyy_xyzzz,  \
                             tr_z_xxyy_yyyy,   \
                             tr_z_xxyy_yyyyy,  \
                             tr_z_xxyy_yyyyz,  \
                             tr_z_xxyy_yyyz,   \
                             tr_z_xxyy_yyyzz,  \
                             tr_z_xxyy_yyzz,   \
                             tr_z_xxyy_yyzzz,  \
                             tr_z_xxyy_yzzz,   \
                             tr_z_xxyy_yzzzz,  \
                             tr_z_xxyy_zzzzz,  \
                             tr_z_xyy_xxxxy,   \
                             tr_z_xyy_xxxyy,   \
                             tr_z_xyy_xxxyz,   \
                             tr_z_xyy_xxyyy,   \
                             tr_z_xyy_xxyyz,   \
                             tr_z_xyy_xxyzz,   \
                             tr_z_xyy_xyyyy,   \
                             tr_z_xyy_xyyyz,   \
                             tr_z_xyy_xyyzz,   \
                             tr_z_xyy_xyzzz,   \
                             tr_z_xyy_yyyyy,   \
                             tr_z_xyy_yyyyz,   \
                             tr_z_xyy_yyyzz,   \
                             tr_z_xyy_yyzzz,   \
                             tr_z_xyy_yzzzz,   \
                             tr_z_xyy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyy_xxxxx[i] = tr_z_xxx_xxxxx[i] * fe_0 + tr_z_xxxy_xxxxx[i] * pa_y[i];

        tr_z_xxxyy_xxxxy[i] = 2.0 * tr_z_xyy_xxxxy[i] * fe_0 + 4.0 * tr_z_xxyy_xxxy[i] * fe_0 + tr_z_xxyy_xxxxy[i] * pa_x[i];

        tr_z_xxxyy_xxxxz[i] = tr_z_xxx_xxxxz[i] * fe_0 + tr_z_xxxy_xxxxz[i] * pa_y[i];

        tr_z_xxxyy_xxxyy[i] = 2.0 * tr_z_xyy_xxxyy[i] * fe_0 + 3.0 * tr_z_xxyy_xxyy[i] * fe_0 + tr_z_xxyy_xxxyy[i] * pa_x[i];

        tr_z_xxxyy_xxxyz[i] = 2.0 * tr_z_xyy_xxxyz[i] * fe_0 + 3.0 * tr_z_xxyy_xxyz[i] * fe_0 + tr_z_xxyy_xxxyz[i] * pa_x[i];

        tr_z_xxxyy_xxxzz[i] = tr_z_xxx_xxxzz[i] * fe_0 + tr_z_xxxy_xxxzz[i] * pa_y[i];

        tr_z_xxxyy_xxyyy[i] = 2.0 * tr_z_xyy_xxyyy[i] * fe_0 + 2.0 * tr_z_xxyy_xyyy[i] * fe_0 + tr_z_xxyy_xxyyy[i] * pa_x[i];

        tr_z_xxxyy_xxyyz[i] = 2.0 * tr_z_xyy_xxyyz[i] * fe_0 + 2.0 * tr_z_xxyy_xyyz[i] * fe_0 + tr_z_xxyy_xxyyz[i] * pa_x[i];

        tr_z_xxxyy_xxyzz[i] = 2.0 * tr_z_xyy_xxyzz[i] * fe_0 + 2.0 * tr_z_xxyy_xyzz[i] * fe_0 + tr_z_xxyy_xxyzz[i] * pa_x[i];

        tr_z_xxxyy_xxzzz[i] = tr_z_xxx_xxzzz[i] * fe_0 + tr_z_xxxy_xxzzz[i] * pa_y[i];

        tr_z_xxxyy_xyyyy[i] = 2.0 * tr_z_xyy_xyyyy[i] * fe_0 + tr_z_xxyy_yyyy[i] * fe_0 + tr_z_xxyy_xyyyy[i] * pa_x[i];

        tr_z_xxxyy_xyyyz[i] = 2.0 * tr_z_xyy_xyyyz[i] * fe_0 + tr_z_xxyy_yyyz[i] * fe_0 + tr_z_xxyy_xyyyz[i] * pa_x[i];

        tr_z_xxxyy_xyyzz[i] = 2.0 * tr_z_xyy_xyyzz[i] * fe_0 + tr_z_xxyy_yyzz[i] * fe_0 + tr_z_xxyy_xyyzz[i] * pa_x[i];

        tr_z_xxxyy_xyzzz[i] = 2.0 * tr_z_xyy_xyzzz[i] * fe_0 + tr_z_xxyy_yzzz[i] * fe_0 + tr_z_xxyy_xyzzz[i] * pa_x[i];

        tr_z_xxxyy_xzzzz[i] = tr_z_xxx_xzzzz[i] * fe_0 + tr_z_xxxy_xzzzz[i] * pa_y[i];

        tr_z_xxxyy_yyyyy[i] = 2.0 * tr_z_xyy_yyyyy[i] * fe_0 + tr_z_xxyy_yyyyy[i] * pa_x[i];

        tr_z_xxxyy_yyyyz[i] = 2.0 * tr_z_xyy_yyyyz[i] * fe_0 + tr_z_xxyy_yyyyz[i] * pa_x[i];

        tr_z_xxxyy_yyyzz[i] = 2.0 * tr_z_xyy_yyyzz[i] * fe_0 + tr_z_xxyy_yyyzz[i] * pa_x[i];

        tr_z_xxxyy_yyzzz[i] = 2.0 * tr_z_xyy_yyzzz[i] * fe_0 + tr_z_xxyy_yyzzz[i] * pa_x[i];

        tr_z_xxxyy_yzzzz[i] = 2.0 * tr_z_xyy_yzzzz[i] * fe_0 + tr_z_xxyy_yzzzz[i] * pa_x[i];

        tr_z_xxxyy_zzzzz[i] = 2.0 * tr_z_xyy_zzzzz[i] * fe_0 + tr_z_xxyy_zzzzz[i] * pa_x[i];
    }

    // Set up 966-987 components of targeted buffer : HH

    auto tr_z_xxxyz_xxxxx = pbuffer.data(idx_dip_hh + 966);

    auto tr_z_xxxyz_xxxxy = pbuffer.data(idx_dip_hh + 967);

    auto tr_z_xxxyz_xxxxz = pbuffer.data(idx_dip_hh + 968);

    auto tr_z_xxxyz_xxxyy = pbuffer.data(idx_dip_hh + 969);

    auto tr_z_xxxyz_xxxyz = pbuffer.data(idx_dip_hh + 970);

    auto tr_z_xxxyz_xxxzz = pbuffer.data(idx_dip_hh + 971);

    auto tr_z_xxxyz_xxyyy = pbuffer.data(idx_dip_hh + 972);

    auto tr_z_xxxyz_xxyyz = pbuffer.data(idx_dip_hh + 973);

    auto tr_z_xxxyz_xxyzz = pbuffer.data(idx_dip_hh + 974);

    auto tr_z_xxxyz_xxzzz = pbuffer.data(idx_dip_hh + 975);

    auto tr_z_xxxyz_xyyyy = pbuffer.data(idx_dip_hh + 976);

    auto tr_z_xxxyz_xyyyz = pbuffer.data(idx_dip_hh + 977);

    auto tr_z_xxxyz_xyyzz = pbuffer.data(idx_dip_hh + 978);

    auto tr_z_xxxyz_xyzzz = pbuffer.data(idx_dip_hh + 979);

    auto tr_z_xxxyz_xzzzz = pbuffer.data(idx_dip_hh + 980);

    auto tr_z_xxxyz_yyyyy = pbuffer.data(idx_dip_hh + 981);

    auto tr_z_xxxyz_yyyyz = pbuffer.data(idx_dip_hh + 982);

    auto tr_z_xxxyz_yyyzz = pbuffer.data(idx_dip_hh + 983);

    auto tr_z_xxxyz_yyzzz = pbuffer.data(idx_dip_hh + 984);

    auto tr_z_xxxyz_yzzzz = pbuffer.data(idx_dip_hh + 985);

    auto tr_z_xxxyz_zzzzz = pbuffer.data(idx_dip_hh + 986);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxxyz_xxxxx, \
                             tr_z_xxxyz_xxxxy, \
                             tr_z_xxxyz_xxxxz, \
                             tr_z_xxxyz_xxxyy, \
                             tr_z_xxxyz_xxxyz, \
                             tr_z_xxxyz_xxxzz, \
                             tr_z_xxxyz_xxyyy, \
                             tr_z_xxxyz_xxyyz, \
                             tr_z_xxxyz_xxyzz, \
                             tr_z_xxxyz_xxzzz, \
                             tr_z_xxxyz_xyyyy, \
                             tr_z_xxxyz_xyyyz, \
                             tr_z_xxxyz_xyyzz, \
                             tr_z_xxxyz_xyzzz, \
                             tr_z_xxxyz_xzzzz, \
                             tr_z_xxxyz_yyyyy, \
                             tr_z_xxxyz_yyyyz, \
                             tr_z_xxxyz_yyyzz, \
                             tr_z_xxxyz_yyzzz, \
                             tr_z_xxxyz_yzzzz, \
                             tr_z_xxxyz_zzzzz, \
                             tr_z_xxxz_xxxx,   \
                             tr_z_xxxz_xxxxx,  \
                             tr_z_xxxz_xxxxy,  \
                             tr_z_xxxz_xxxxz,  \
                             tr_z_xxxz_xxxy,   \
                             tr_z_xxxz_xxxyy,  \
                             tr_z_xxxz_xxxyz,  \
                             tr_z_xxxz_xxxz,   \
                             tr_z_xxxz_xxxzz,  \
                             tr_z_xxxz_xxyy,   \
                             tr_z_xxxz_xxyyy,  \
                             tr_z_xxxz_xxyyz,  \
                             tr_z_xxxz_xxyz,   \
                             tr_z_xxxz_xxyzz,  \
                             tr_z_xxxz_xxzz,   \
                             tr_z_xxxz_xxzzz,  \
                             tr_z_xxxz_xyyy,   \
                             tr_z_xxxz_xyyyy,  \
                             tr_z_xxxz_xyyyz,  \
                             tr_z_xxxz_xyyz,   \
                             tr_z_xxxz_xyyzz,  \
                             tr_z_xxxz_xyzz,   \
                             tr_z_xxxz_xyzzz,  \
                             tr_z_xxxz_xzzz,   \
                             tr_z_xxxz_xzzzz,  \
                             tr_z_xxxz_zzzzz,  \
                             tr_z_xxyz_yyyyy,  \
                             tr_z_xxyz_yyyyz,  \
                             tr_z_xxyz_yyyzz,  \
                             tr_z_xxyz_yyzzz,  \
                             tr_z_xxyz_yzzzz,  \
                             tr_z_xyz_yyyyy,   \
                             tr_z_xyz_yyyyz,   \
                             tr_z_xyz_yyyzz,   \
                             tr_z_xyz_yyzzz,   \
                             tr_z_xyz_yzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyz_xxxxx[i] = tr_z_xxxz_xxxxx[i] * pa_y[i];

        tr_z_xxxyz_xxxxy[i] = tr_z_xxxz_xxxx[i] * fe_0 + tr_z_xxxz_xxxxy[i] * pa_y[i];

        tr_z_xxxyz_xxxxz[i] = tr_z_xxxz_xxxxz[i] * pa_y[i];

        tr_z_xxxyz_xxxyy[i] = 2.0 * tr_z_xxxz_xxxy[i] * fe_0 + tr_z_xxxz_xxxyy[i] * pa_y[i];

        tr_z_xxxyz_xxxyz[i] = tr_z_xxxz_xxxz[i] * fe_0 + tr_z_xxxz_xxxyz[i] * pa_y[i];

        tr_z_xxxyz_xxxzz[i] = tr_z_xxxz_xxxzz[i] * pa_y[i];

        tr_z_xxxyz_xxyyy[i] = 3.0 * tr_z_xxxz_xxyy[i] * fe_0 + tr_z_xxxz_xxyyy[i] * pa_y[i];

        tr_z_xxxyz_xxyyz[i] = 2.0 * tr_z_xxxz_xxyz[i] * fe_0 + tr_z_xxxz_xxyyz[i] * pa_y[i];

        tr_z_xxxyz_xxyzz[i] = tr_z_xxxz_xxzz[i] * fe_0 + tr_z_xxxz_xxyzz[i] * pa_y[i];

        tr_z_xxxyz_xxzzz[i] = tr_z_xxxz_xxzzz[i] * pa_y[i];

        tr_z_xxxyz_xyyyy[i] = 4.0 * tr_z_xxxz_xyyy[i] * fe_0 + tr_z_xxxz_xyyyy[i] * pa_y[i];

        tr_z_xxxyz_xyyyz[i] = 3.0 * tr_z_xxxz_xyyz[i] * fe_0 + tr_z_xxxz_xyyyz[i] * pa_y[i];

        tr_z_xxxyz_xyyzz[i] = 2.0 * tr_z_xxxz_xyzz[i] * fe_0 + tr_z_xxxz_xyyzz[i] * pa_y[i];

        tr_z_xxxyz_xyzzz[i] = tr_z_xxxz_xzzz[i] * fe_0 + tr_z_xxxz_xyzzz[i] * pa_y[i];

        tr_z_xxxyz_xzzzz[i] = tr_z_xxxz_xzzzz[i] * pa_y[i];

        tr_z_xxxyz_yyyyy[i] = 2.0 * tr_z_xyz_yyyyy[i] * fe_0 + tr_z_xxyz_yyyyy[i] * pa_x[i];

        tr_z_xxxyz_yyyyz[i] = 2.0 * tr_z_xyz_yyyyz[i] * fe_0 + tr_z_xxyz_yyyyz[i] * pa_x[i];

        tr_z_xxxyz_yyyzz[i] = 2.0 * tr_z_xyz_yyyzz[i] * fe_0 + tr_z_xxyz_yyyzz[i] * pa_x[i];

        tr_z_xxxyz_yyzzz[i] = 2.0 * tr_z_xyz_yyzzz[i] * fe_0 + tr_z_xxyz_yyzzz[i] * pa_x[i];

        tr_z_xxxyz_yzzzz[i] = 2.0 * tr_z_xyz_yzzzz[i] * fe_0 + tr_z_xxyz_yzzzz[i] * pa_x[i];

        tr_z_xxxyz_zzzzz[i] = tr_z_xxxz_zzzzz[i] * pa_y[i];
    }

    // Set up 987-1008 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxxzz_xxxxx, \
                             tr_z_xxxzz_xxxxy, \
                             tr_z_xxxzz_xxxxz, \
                             tr_z_xxxzz_xxxyy, \
                             tr_z_xxxzz_xxxyz, \
                             tr_z_xxxzz_xxxzz, \
                             tr_z_xxxzz_xxyyy, \
                             tr_z_xxxzz_xxyyz, \
                             tr_z_xxxzz_xxyzz, \
                             tr_z_xxxzz_xxzzz, \
                             tr_z_xxxzz_xyyyy, \
                             tr_z_xxxzz_xyyyz, \
                             tr_z_xxxzz_xyyzz, \
                             tr_z_xxxzz_xyzzz, \
                             tr_z_xxxzz_xzzzz, \
                             tr_z_xxxzz_yyyyy, \
                             tr_z_xxxzz_yyyyz, \
                             tr_z_xxxzz_yyyzz, \
                             tr_z_xxxzz_yyzzz, \
                             tr_z_xxxzz_yzzzz, \
                             tr_z_xxxzz_zzzzz, \
                             tr_z_xxzz_xxxx,   \
                             tr_z_xxzz_xxxxx,  \
                             tr_z_xxzz_xxxxy,  \
                             tr_z_xxzz_xxxxz,  \
                             tr_z_xxzz_xxxy,   \
                             tr_z_xxzz_xxxyy,  \
                             tr_z_xxzz_xxxyz,  \
                             tr_z_xxzz_xxxz,   \
                             tr_z_xxzz_xxxzz,  \
                             tr_z_xxzz_xxyy,   \
                             tr_z_xxzz_xxyyy,  \
                             tr_z_xxzz_xxyyz,  \
                             tr_z_xxzz_xxyz,   \
                             tr_z_xxzz_xxyzz,  \
                             tr_z_xxzz_xxzz,   \
                             tr_z_xxzz_xxzzz,  \
                             tr_z_xxzz_xyyy,   \
                             tr_z_xxzz_xyyyy,  \
                             tr_z_xxzz_xyyyz,  \
                             tr_z_xxzz_xyyz,   \
                             tr_z_xxzz_xyyzz,  \
                             tr_z_xxzz_xyzz,   \
                             tr_z_xxzz_xyzzz,  \
                             tr_z_xxzz_xzzz,   \
                             tr_z_xxzz_xzzzz,  \
                             tr_z_xxzz_yyyy,   \
                             tr_z_xxzz_yyyyy,  \
                             tr_z_xxzz_yyyyz,  \
                             tr_z_xxzz_yyyz,   \
                             tr_z_xxzz_yyyzz,  \
                             tr_z_xxzz_yyzz,   \
                             tr_z_xxzz_yyzzz,  \
                             tr_z_xxzz_yzzz,   \
                             tr_z_xxzz_yzzzz,  \
                             tr_z_xxzz_zzzz,   \
                             tr_z_xxzz_zzzzz,  \
                             tr_z_xzz_xxxxx,   \
                             tr_z_xzz_xxxxy,   \
                             tr_z_xzz_xxxxz,   \
                             tr_z_xzz_xxxyy,   \
                             tr_z_xzz_xxxyz,   \
                             tr_z_xzz_xxxzz,   \
                             tr_z_xzz_xxyyy,   \
                             tr_z_xzz_xxyyz,   \
                             tr_z_xzz_xxyzz,   \
                             tr_z_xzz_xxzzz,   \
                             tr_z_xzz_xyyyy,   \
                             tr_z_xzz_xyyyz,   \
                             tr_z_xzz_xyyzz,   \
                             tr_z_xzz_xyzzz,   \
                             tr_z_xzz_xzzzz,   \
                             tr_z_xzz_yyyyy,   \
                             tr_z_xzz_yyyyz,   \
                             tr_z_xzz_yyyzz,   \
                             tr_z_xzz_yyzzz,   \
                             tr_z_xzz_yzzzz,   \
                             tr_z_xzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzz_xxxxx[i] = 2.0 * tr_z_xzz_xxxxx[i] * fe_0 + 5.0 * tr_z_xxzz_xxxx[i] * fe_0 + tr_z_xxzz_xxxxx[i] * pa_x[i];

        tr_z_xxxzz_xxxxy[i] = 2.0 * tr_z_xzz_xxxxy[i] * fe_0 + 4.0 * tr_z_xxzz_xxxy[i] * fe_0 + tr_z_xxzz_xxxxy[i] * pa_x[i];

        tr_z_xxxzz_xxxxz[i] = 2.0 * tr_z_xzz_xxxxz[i] * fe_0 + 4.0 * tr_z_xxzz_xxxz[i] * fe_0 + tr_z_xxzz_xxxxz[i] * pa_x[i];

        tr_z_xxxzz_xxxyy[i] = 2.0 * tr_z_xzz_xxxyy[i] * fe_0 + 3.0 * tr_z_xxzz_xxyy[i] * fe_0 + tr_z_xxzz_xxxyy[i] * pa_x[i];

        tr_z_xxxzz_xxxyz[i] = 2.0 * tr_z_xzz_xxxyz[i] * fe_0 + 3.0 * tr_z_xxzz_xxyz[i] * fe_0 + tr_z_xxzz_xxxyz[i] * pa_x[i];

        tr_z_xxxzz_xxxzz[i] = 2.0 * tr_z_xzz_xxxzz[i] * fe_0 + 3.0 * tr_z_xxzz_xxzz[i] * fe_0 + tr_z_xxzz_xxxzz[i] * pa_x[i];

        tr_z_xxxzz_xxyyy[i] = 2.0 * tr_z_xzz_xxyyy[i] * fe_0 + 2.0 * tr_z_xxzz_xyyy[i] * fe_0 + tr_z_xxzz_xxyyy[i] * pa_x[i];

        tr_z_xxxzz_xxyyz[i] = 2.0 * tr_z_xzz_xxyyz[i] * fe_0 + 2.0 * tr_z_xxzz_xyyz[i] * fe_0 + tr_z_xxzz_xxyyz[i] * pa_x[i];

        tr_z_xxxzz_xxyzz[i] = 2.0 * tr_z_xzz_xxyzz[i] * fe_0 + 2.0 * tr_z_xxzz_xyzz[i] * fe_0 + tr_z_xxzz_xxyzz[i] * pa_x[i];

        tr_z_xxxzz_xxzzz[i] = 2.0 * tr_z_xzz_xxzzz[i] * fe_0 + 2.0 * tr_z_xxzz_xzzz[i] * fe_0 + tr_z_xxzz_xxzzz[i] * pa_x[i];

        tr_z_xxxzz_xyyyy[i] = 2.0 * tr_z_xzz_xyyyy[i] * fe_0 + tr_z_xxzz_yyyy[i] * fe_0 + tr_z_xxzz_xyyyy[i] * pa_x[i];

        tr_z_xxxzz_xyyyz[i] = 2.0 * tr_z_xzz_xyyyz[i] * fe_0 + tr_z_xxzz_yyyz[i] * fe_0 + tr_z_xxzz_xyyyz[i] * pa_x[i];

        tr_z_xxxzz_xyyzz[i] = 2.0 * tr_z_xzz_xyyzz[i] * fe_0 + tr_z_xxzz_yyzz[i] * fe_0 + tr_z_xxzz_xyyzz[i] * pa_x[i];

        tr_z_xxxzz_xyzzz[i] = 2.0 * tr_z_xzz_xyzzz[i] * fe_0 + tr_z_xxzz_yzzz[i] * fe_0 + tr_z_xxzz_xyzzz[i] * pa_x[i];

        tr_z_xxxzz_xzzzz[i] = 2.0 * tr_z_xzz_xzzzz[i] * fe_0 + tr_z_xxzz_zzzz[i] * fe_0 + tr_z_xxzz_xzzzz[i] * pa_x[i];

        tr_z_xxxzz_yyyyy[i] = 2.0 * tr_z_xzz_yyyyy[i] * fe_0 + tr_z_xxzz_yyyyy[i] * pa_x[i];

        tr_z_xxxzz_yyyyz[i] = 2.0 * tr_z_xzz_yyyyz[i] * fe_0 + tr_z_xxzz_yyyyz[i] * pa_x[i];

        tr_z_xxxzz_yyyzz[i] = 2.0 * tr_z_xzz_yyyzz[i] * fe_0 + tr_z_xxzz_yyyzz[i] * pa_x[i];

        tr_z_xxxzz_yyzzz[i] = 2.0 * tr_z_xzz_yyzzz[i] * fe_0 + tr_z_xxzz_yyzzz[i] * pa_x[i];

        tr_z_xxxzz_yzzzz[i] = 2.0 * tr_z_xzz_yzzzz[i] * fe_0 + tr_z_xxzz_yzzzz[i] * pa_x[i];

        tr_z_xxxzz_zzzzz[i] = 2.0 * tr_z_xzz_zzzzz[i] * fe_0 + tr_z_xxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1008-1029 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxy_xxxxx,   \
                             tr_z_xxy_xxxxz,   \
                             tr_z_xxy_xxxzz,   \
                             tr_z_xxy_xxzzz,   \
                             tr_z_xxy_xzzzz,   \
                             tr_z_xxyy_xxxxx,  \
                             tr_z_xxyy_xxxxz,  \
                             tr_z_xxyy_xxxzz,  \
                             tr_z_xxyy_xxzzz,  \
                             tr_z_xxyy_xzzzz,  \
                             tr_z_xxyyy_xxxxx, \
                             tr_z_xxyyy_xxxxy, \
                             tr_z_xxyyy_xxxxz, \
                             tr_z_xxyyy_xxxyy, \
                             tr_z_xxyyy_xxxyz, \
                             tr_z_xxyyy_xxxzz, \
                             tr_z_xxyyy_xxyyy, \
                             tr_z_xxyyy_xxyyz, \
                             tr_z_xxyyy_xxyzz, \
                             tr_z_xxyyy_xxzzz, \
                             tr_z_xxyyy_xyyyy, \
                             tr_z_xxyyy_xyyyz, \
                             tr_z_xxyyy_xyyzz, \
                             tr_z_xxyyy_xyzzz, \
                             tr_z_xxyyy_xzzzz, \
                             tr_z_xxyyy_yyyyy, \
                             tr_z_xxyyy_yyyyz, \
                             tr_z_xxyyy_yyyzz, \
                             tr_z_xxyyy_yyzzz, \
                             tr_z_xxyyy_yzzzz, \
                             tr_z_xxyyy_zzzzz, \
                             tr_z_xyyy_xxxxy,  \
                             tr_z_xyyy_xxxy,   \
                             tr_z_xyyy_xxxyy,  \
                             tr_z_xyyy_xxxyz,  \
                             tr_z_xyyy_xxyy,   \
                             tr_z_xyyy_xxyyy,  \
                             tr_z_xyyy_xxyyz,  \
                             tr_z_xyyy_xxyz,   \
                             tr_z_xyyy_xxyzz,  \
                             tr_z_xyyy_xyyy,   \
                             tr_z_xyyy_xyyyy,  \
                             tr_z_xyyy_xyyyz,  \
                             tr_z_xyyy_xyyz,   \
                             tr_z_xyyy_xyyzz,  \
                             tr_z_xyyy_xyzz,   \
                             tr_z_xyyy_xyzzz,  \
                             tr_z_xyyy_yyyy,   \
                             tr_z_xyyy_yyyyy,  \
                             tr_z_xyyy_yyyyz,  \
                             tr_z_xyyy_yyyz,   \
                             tr_z_xyyy_yyyzz,  \
                             tr_z_xyyy_yyzz,   \
                             tr_z_xyyy_yyzzz,  \
                             tr_z_xyyy_yzzz,   \
                             tr_z_xyyy_yzzzz,  \
                             tr_z_xyyy_zzzzz,  \
                             tr_z_yyy_xxxxy,   \
                             tr_z_yyy_xxxyy,   \
                             tr_z_yyy_xxxyz,   \
                             tr_z_yyy_xxyyy,   \
                             tr_z_yyy_xxyyz,   \
                             tr_z_yyy_xxyzz,   \
                             tr_z_yyy_xyyyy,   \
                             tr_z_yyy_xyyyz,   \
                             tr_z_yyy_xyyzz,   \
                             tr_z_yyy_xyzzz,   \
                             tr_z_yyy_yyyyy,   \
                             tr_z_yyy_yyyyz,   \
                             tr_z_yyy_yyyzz,   \
                             tr_z_yyy_yyzzz,   \
                             tr_z_yyy_yzzzz,   \
                             tr_z_yyy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyy_xxxxx[i] = 2.0 * tr_z_xxy_xxxxx[i] * fe_0 + tr_z_xxyy_xxxxx[i] * pa_y[i];

        tr_z_xxyyy_xxxxy[i] = tr_z_yyy_xxxxy[i] * fe_0 + 4.0 * tr_z_xyyy_xxxy[i] * fe_0 + tr_z_xyyy_xxxxy[i] * pa_x[i];

        tr_z_xxyyy_xxxxz[i] = 2.0 * tr_z_xxy_xxxxz[i] * fe_0 + tr_z_xxyy_xxxxz[i] * pa_y[i];

        tr_z_xxyyy_xxxyy[i] = tr_z_yyy_xxxyy[i] * fe_0 + 3.0 * tr_z_xyyy_xxyy[i] * fe_0 + tr_z_xyyy_xxxyy[i] * pa_x[i];

        tr_z_xxyyy_xxxyz[i] = tr_z_yyy_xxxyz[i] * fe_0 + 3.0 * tr_z_xyyy_xxyz[i] * fe_0 + tr_z_xyyy_xxxyz[i] * pa_x[i];

        tr_z_xxyyy_xxxzz[i] = 2.0 * tr_z_xxy_xxxzz[i] * fe_0 + tr_z_xxyy_xxxzz[i] * pa_y[i];

        tr_z_xxyyy_xxyyy[i] = tr_z_yyy_xxyyy[i] * fe_0 + 2.0 * tr_z_xyyy_xyyy[i] * fe_0 + tr_z_xyyy_xxyyy[i] * pa_x[i];

        tr_z_xxyyy_xxyyz[i] = tr_z_yyy_xxyyz[i] * fe_0 + 2.0 * tr_z_xyyy_xyyz[i] * fe_0 + tr_z_xyyy_xxyyz[i] * pa_x[i];

        tr_z_xxyyy_xxyzz[i] = tr_z_yyy_xxyzz[i] * fe_0 + 2.0 * tr_z_xyyy_xyzz[i] * fe_0 + tr_z_xyyy_xxyzz[i] * pa_x[i];

        tr_z_xxyyy_xxzzz[i] = 2.0 * tr_z_xxy_xxzzz[i] * fe_0 + tr_z_xxyy_xxzzz[i] * pa_y[i];

        tr_z_xxyyy_xyyyy[i] = tr_z_yyy_xyyyy[i] * fe_0 + tr_z_xyyy_yyyy[i] * fe_0 + tr_z_xyyy_xyyyy[i] * pa_x[i];

        tr_z_xxyyy_xyyyz[i] = tr_z_yyy_xyyyz[i] * fe_0 + tr_z_xyyy_yyyz[i] * fe_0 + tr_z_xyyy_xyyyz[i] * pa_x[i];

        tr_z_xxyyy_xyyzz[i] = tr_z_yyy_xyyzz[i] * fe_0 + tr_z_xyyy_yyzz[i] * fe_0 + tr_z_xyyy_xyyzz[i] * pa_x[i];

        tr_z_xxyyy_xyzzz[i] = tr_z_yyy_xyzzz[i] * fe_0 + tr_z_xyyy_yzzz[i] * fe_0 + tr_z_xyyy_xyzzz[i] * pa_x[i];

        tr_z_xxyyy_xzzzz[i] = 2.0 * tr_z_xxy_xzzzz[i] * fe_0 + tr_z_xxyy_xzzzz[i] * pa_y[i];

        tr_z_xxyyy_yyyyy[i] = tr_z_yyy_yyyyy[i] * fe_0 + tr_z_xyyy_yyyyy[i] * pa_x[i];

        tr_z_xxyyy_yyyyz[i] = tr_z_yyy_yyyyz[i] * fe_0 + tr_z_xyyy_yyyyz[i] * pa_x[i];

        tr_z_xxyyy_yyyzz[i] = tr_z_yyy_yyyzz[i] * fe_0 + tr_z_xyyy_yyyzz[i] * pa_x[i];

        tr_z_xxyyy_yyzzz[i] = tr_z_yyy_yyzzz[i] * fe_0 + tr_z_xyyy_yyzzz[i] * pa_x[i];

        tr_z_xxyyy_yzzzz[i] = tr_z_yyy_yzzzz[i] * fe_0 + tr_z_xyyy_yzzzz[i] * pa_x[i];

        tr_z_xxyyy_zzzzz[i] = tr_z_yyy_zzzzz[i] * fe_0 + tr_z_xyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 1029-1050 components of targeted buffer : HH

    auto tr_z_xxyyz_xxxxx = pbuffer.data(idx_dip_hh + 1029);

    auto tr_z_xxyyz_xxxxy = pbuffer.data(idx_dip_hh + 1030);

    auto tr_z_xxyyz_xxxxz = pbuffer.data(idx_dip_hh + 1031);

    auto tr_z_xxyyz_xxxyy = pbuffer.data(idx_dip_hh + 1032);

    auto tr_z_xxyyz_xxxyz = pbuffer.data(idx_dip_hh + 1033);

    auto tr_z_xxyyz_xxxzz = pbuffer.data(idx_dip_hh + 1034);

    auto tr_z_xxyyz_xxyyy = pbuffer.data(idx_dip_hh + 1035);

    auto tr_z_xxyyz_xxyyz = pbuffer.data(idx_dip_hh + 1036);

    auto tr_z_xxyyz_xxyzz = pbuffer.data(idx_dip_hh + 1037);

    auto tr_z_xxyyz_xxzzz = pbuffer.data(idx_dip_hh + 1038);

    auto tr_z_xxyyz_xyyyy = pbuffer.data(idx_dip_hh + 1039);

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

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             tr_z_xxyy_xxxxy,  \
                             tr_z_xxyy_xxxyy,  \
                             tr_z_xxyy_xxyyy,  \
                             tr_z_xxyy_xyyyy,  \
                             tr_z_xxyyz_xxxxx, \
                             tr_z_xxyyz_xxxxy, \
                             tr_z_xxyyz_xxxxz, \
                             tr_z_xxyyz_xxxyy, \
                             tr_z_xxyyz_xxxyz, \
                             tr_z_xxyyz_xxxzz, \
                             tr_z_xxyyz_xxyyy, \
                             tr_z_xxyyz_xxyyz, \
                             tr_z_xxyyz_xxyzz, \
                             tr_z_xxyyz_xxzzz, \
                             tr_z_xxyyz_xyyyy, \
                             tr_z_xxyyz_xyyyz, \
                             tr_z_xxyyz_xyyzz, \
                             tr_z_xxyyz_xyzzz, \
                             tr_z_xxyyz_xzzzz, \
                             tr_z_xxyyz_yyyyy, \
                             tr_z_xxyyz_yyyyz, \
                             tr_z_xxyyz_yyyzz, \
                             tr_z_xxyyz_yyzzz, \
                             tr_z_xxyyz_yzzzz, \
                             tr_z_xxyyz_zzzzz, \
                             tr_z_xxyz_xxxxx,  \
                             tr_z_xxyz_xxxxz,  \
                             tr_z_xxyz_xxxzz,  \
                             tr_z_xxyz_xxzzz,  \
                             tr_z_xxyz_xzzzz,  \
                             tr_z_xxz_xxxxx,   \
                             tr_z_xxz_xxxxz,   \
                             tr_z_xxz_xxxzz,   \
                             tr_z_xxz_xxzzz,   \
                             tr_z_xxz_xzzzz,   \
                             tr_z_xyyz_xxxyz,  \
                             tr_z_xyyz_xxyyz,  \
                             tr_z_xyyz_xxyz,   \
                             tr_z_xyyz_xxyzz,  \
                             tr_z_xyyz_xyyyz,  \
                             tr_z_xyyz_xyyz,   \
                             tr_z_xyyz_xyyzz,  \
                             tr_z_xyyz_xyzz,   \
                             tr_z_xyyz_xyzzz,  \
                             tr_z_xyyz_yyyyy,  \
                             tr_z_xyyz_yyyyz,  \
                             tr_z_xyyz_yyyz,   \
                             tr_z_xyyz_yyyzz,  \
                             tr_z_xyyz_yyzz,   \
                             tr_z_xyyz_yyzzz,  \
                             tr_z_xyyz_yzzz,   \
                             tr_z_xyyz_yzzzz,  \
                             tr_z_xyyz_zzzzz,  \
                             tr_z_yyz_xxxyz,   \
                             tr_z_yyz_xxyyz,   \
                             tr_z_yyz_xxyzz,   \
                             tr_z_yyz_xyyyz,   \
                             tr_z_yyz_xyyzz,   \
                             tr_z_yyz_xyzzz,   \
                             tr_z_yyz_yyyyy,   \
                             tr_z_yyz_yyyyz,   \
                             tr_z_yyz_yyyzz,   \
                             tr_z_yyz_yyzzz,   \
                             tr_z_yyz_yzzzz,   \
                             tr_z_yyz_zzzzz,   \
                             ts_xxyy_xxxxy,    \
                             ts_xxyy_xxxyy,    \
                             ts_xxyy_xxyyy,    \
                             ts_xxyy_xyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyz_xxxxx[i] = tr_z_xxz_xxxxx[i] * fe_0 + tr_z_xxyz_xxxxx[i] * pa_y[i];

        tr_z_xxyyz_xxxxy[i] = ts_xxyy_xxxxy[i] * fe_0 + tr_z_xxyy_xxxxy[i] * pa_z[i];

        tr_z_xxyyz_xxxxz[i] = tr_z_xxz_xxxxz[i] * fe_0 + tr_z_xxyz_xxxxz[i] * pa_y[i];

        tr_z_xxyyz_xxxyy[i] = ts_xxyy_xxxyy[i] * fe_0 + tr_z_xxyy_xxxyy[i] * pa_z[i];

        tr_z_xxyyz_xxxyz[i] = tr_z_yyz_xxxyz[i] * fe_0 + 3.0 * tr_z_xyyz_xxyz[i] * fe_0 + tr_z_xyyz_xxxyz[i] * pa_x[i];

        tr_z_xxyyz_xxxzz[i] = tr_z_xxz_xxxzz[i] * fe_0 + tr_z_xxyz_xxxzz[i] * pa_y[i];

        tr_z_xxyyz_xxyyy[i] = ts_xxyy_xxyyy[i] * fe_0 + tr_z_xxyy_xxyyy[i] * pa_z[i];

        tr_z_xxyyz_xxyyz[i] = tr_z_yyz_xxyyz[i] * fe_0 + 2.0 * tr_z_xyyz_xyyz[i] * fe_0 + tr_z_xyyz_xxyyz[i] * pa_x[i];

        tr_z_xxyyz_xxyzz[i] = tr_z_yyz_xxyzz[i] * fe_0 + 2.0 * tr_z_xyyz_xyzz[i] * fe_0 + tr_z_xyyz_xxyzz[i] * pa_x[i];

        tr_z_xxyyz_xxzzz[i] = tr_z_xxz_xxzzz[i] * fe_0 + tr_z_xxyz_xxzzz[i] * pa_y[i];

        tr_z_xxyyz_xyyyy[i] = ts_xxyy_xyyyy[i] * fe_0 + tr_z_xxyy_xyyyy[i] * pa_z[i];

        tr_z_xxyyz_xyyyz[i] = tr_z_yyz_xyyyz[i] * fe_0 + tr_z_xyyz_yyyz[i] * fe_0 + tr_z_xyyz_xyyyz[i] * pa_x[i];

        tr_z_xxyyz_xyyzz[i] = tr_z_yyz_xyyzz[i] * fe_0 + tr_z_xyyz_yyzz[i] * fe_0 + tr_z_xyyz_xyyzz[i] * pa_x[i];

        tr_z_xxyyz_xyzzz[i] = tr_z_yyz_xyzzz[i] * fe_0 + tr_z_xyyz_yzzz[i] * fe_0 + tr_z_xyyz_xyzzz[i] * pa_x[i];

        tr_z_xxyyz_xzzzz[i] = tr_z_xxz_xzzzz[i] * fe_0 + tr_z_xxyz_xzzzz[i] * pa_y[i];

        tr_z_xxyyz_yyyyy[i] = tr_z_yyz_yyyyy[i] * fe_0 + tr_z_xyyz_yyyyy[i] * pa_x[i];

        tr_z_xxyyz_yyyyz[i] = tr_z_yyz_yyyyz[i] * fe_0 + tr_z_xyyz_yyyyz[i] * pa_x[i];

        tr_z_xxyyz_yyyzz[i] = tr_z_yyz_yyyzz[i] * fe_0 + tr_z_xyyz_yyyzz[i] * pa_x[i];

        tr_z_xxyyz_yyzzz[i] = tr_z_yyz_yyzzz[i] * fe_0 + tr_z_xyyz_yyzzz[i] * pa_x[i];

        tr_z_xxyyz_yzzzz[i] = tr_z_yyz_yzzzz[i] * fe_0 + tr_z_xyyz_yzzzz[i] * pa_x[i];

        tr_z_xxyyz_zzzzz[i] = tr_z_yyz_zzzzz[i] * fe_0 + tr_z_xyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 1050-1071 components of targeted buffer : HH

    auto tr_z_xxyzz_xxxxx = pbuffer.data(idx_dip_hh + 1050);

    auto tr_z_xxyzz_xxxxy = pbuffer.data(idx_dip_hh + 1051);

    auto tr_z_xxyzz_xxxxz = pbuffer.data(idx_dip_hh + 1052);

    auto tr_z_xxyzz_xxxyy = pbuffer.data(idx_dip_hh + 1053);

    auto tr_z_xxyzz_xxxyz = pbuffer.data(idx_dip_hh + 1054);

    auto tr_z_xxyzz_xxxzz = pbuffer.data(idx_dip_hh + 1055);

    auto tr_z_xxyzz_xxyyy = pbuffer.data(idx_dip_hh + 1056);

    auto tr_z_xxyzz_xxyyz = pbuffer.data(idx_dip_hh + 1057);

    auto tr_z_xxyzz_xxyzz = pbuffer.data(idx_dip_hh + 1058);

    auto tr_z_xxyzz_xxzzz = pbuffer.data(idx_dip_hh + 1059);

    auto tr_z_xxyzz_xyyyy = pbuffer.data(idx_dip_hh + 1060);

    auto tr_z_xxyzz_xyyyz = pbuffer.data(idx_dip_hh + 1061);

    auto tr_z_xxyzz_xyyzz = pbuffer.data(idx_dip_hh + 1062);

    auto tr_z_xxyzz_xyzzz = pbuffer.data(idx_dip_hh + 1063);

    auto tr_z_xxyzz_xzzzz = pbuffer.data(idx_dip_hh + 1064);

    auto tr_z_xxyzz_yyyyy = pbuffer.data(idx_dip_hh + 1065);

    auto tr_z_xxyzz_yyyyz = pbuffer.data(idx_dip_hh + 1066);

    auto tr_z_xxyzz_yyyzz = pbuffer.data(idx_dip_hh + 1067);

    auto tr_z_xxyzz_yyzzz = pbuffer.data(idx_dip_hh + 1068);

    auto tr_z_xxyzz_yzzzz = pbuffer.data(idx_dip_hh + 1069);

    auto tr_z_xxyzz_zzzzz = pbuffer.data(idx_dip_hh + 1070);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xxyzz_xxxxx, \
                             tr_z_xxyzz_xxxxy, \
                             tr_z_xxyzz_xxxxz, \
                             tr_z_xxyzz_xxxyy, \
                             tr_z_xxyzz_xxxyz, \
                             tr_z_xxyzz_xxxzz, \
                             tr_z_xxyzz_xxyyy, \
                             tr_z_xxyzz_xxyyz, \
                             tr_z_xxyzz_xxyzz, \
                             tr_z_xxyzz_xxzzz, \
                             tr_z_xxyzz_xyyyy, \
                             tr_z_xxyzz_xyyyz, \
                             tr_z_xxyzz_xyyzz, \
                             tr_z_xxyzz_xyzzz, \
                             tr_z_xxyzz_xzzzz, \
                             tr_z_xxyzz_yyyyy, \
                             tr_z_xxyzz_yyyyz, \
                             tr_z_xxyzz_yyyzz, \
                             tr_z_xxyzz_yyzzz, \
                             tr_z_xxyzz_yzzzz, \
                             tr_z_xxyzz_zzzzz, \
                             tr_z_xxzz_xxxx,   \
                             tr_z_xxzz_xxxxx,  \
                             tr_z_xxzz_xxxxy,  \
                             tr_z_xxzz_xxxxz,  \
                             tr_z_xxzz_xxxy,   \
                             tr_z_xxzz_xxxyy,  \
                             tr_z_xxzz_xxxyz,  \
                             tr_z_xxzz_xxxz,   \
                             tr_z_xxzz_xxxzz,  \
                             tr_z_xxzz_xxyy,   \
                             tr_z_xxzz_xxyyy,  \
                             tr_z_xxzz_xxyyz,  \
                             tr_z_xxzz_xxyz,   \
                             tr_z_xxzz_xxyzz,  \
                             tr_z_xxzz_xxzz,   \
                             tr_z_xxzz_xxzzz,  \
                             tr_z_xxzz_xyyy,   \
                             tr_z_xxzz_xyyyy,  \
                             tr_z_xxzz_xyyyz,  \
                             tr_z_xxzz_xyyz,   \
                             tr_z_xxzz_xyyzz,  \
                             tr_z_xxzz_xyzz,   \
                             tr_z_xxzz_xyzzz,  \
                             tr_z_xxzz_xzzz,   \
                             tr_z_xxzz_xzzzz,  \
                             tr_z_xxzz_zzzzz,  \
                             tr_z_xyzz_yyyyy,  \
                             tr_z_xyzz_yyyyz,  \
                             tr_z_xyzz_yyyzz,  \
                             tr_z_xyzz_yyzzz,  \
                             tr_z_xyzz_yzzzz,  \
                             tr_z_yzz_yyyyy,   \
                             tr_z_yzz_yyyyz,   \
                             tr_z_yzz_yyyzz,   \
                             tr_z_yzz_yyzzz,   \
                             tr_z_yzz_yzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzz_xxxxx[i] = tr_z_xxzz_xxxxx[i] * pa_y[i];

        tr_z_xxyzz_xxxxy[i] = tr_z_xxzz_xxxx[i] * fe_0 + tr_z_xxzz_xxxxy[i] * pa_y[i];

        tr_z_xxyzz_xxxxz[i] = tr_z_xxzz_xxxxz[i] * pa_y[i];

        tr_z_xxyzz_xxxyy[i] = 2.0 * tr_z_xxzz_xxxy[i] * fe_0 + tr_z_xxzz_xxxyy[i] * pa_y[i];

        tr_z_xxyzz_xxxyz[i] = tr_z_xxzz_xxxz[i] * fe_0 + tr_z_xxzz_xxxyz[i] * pa_y[i];

        tr_z_xxyzz_xxxzz[i] = tr_z_xxzz_xxxzz[i] * pa_y[i];

        tr_z_xxyzz_xxyyy[i] = 3.0 * tr_z_xxzz_xxyy[i] * fe_0 + tr_z_xxzz_xxyyy[i] * pa_y[i];

        tr_z_xxyzz_xxyyz[i] = 2.0 * tr_z_xxzz_xxyz[i] * fe_0 + tr_z_xxzz_xxyyz[i] * pa_y[i];

        tr_z_xxyzz_xxyzz[i] = tr_z_xxzz_xxzz[i] * fe_0 + tr_z_xxzz_xxyzz[i] * pa_y[i];

        tr_z_xxyzz_xxzzz[i] = tr_z_xxzz_xxzzz[i] * pa_y[i];

        tr_z_xxyzz_xyyyy[i] = 4.0 * tr_z_xxzz_xyyy[i] * fe_0 + tr_z_xxzz_xyyyy[i] * pa_y[i];

        tr_z_xxyzz_xyyyz[i] = 3.0 * tr_z_xxzz_xyyz[i] * fe_0 + tr_z_xxzz_xyyyz[i] * pa_y[i];

        tr_z_xxyzz_xyyzz[i] = 2.0 * tr_z_xxzz_xyzz[i] * fe_0 + tr_z_xxzz_xyyzz[i] * pa_y[i];

        tr_z_xxyzz_xyzzz[i] = tr_z_xxzz_xzzz[i] * fe_0 + tr_z_xxzz_xyzzz[i] * pa_y[i];

        tr_z_xxyzz_xzzzz[i] = tr_z_xxzz_xzzzz[i] * pa_y[i];

        tr_z_xxyzz_yyyyy[i] = tr_z_yzz_yyyyy[i] * fe_0 + tr_z_xyzz_yyyyy[i] * pa_x[i];

        tr_z_xxyzz_yyyyz[i] = tr_z_yzz_yyyyz[i] * fe_0 + tr_z_xyzz_yyyyz[i] * pa_x[i];

        tr_z_xxyzz_yyyzz[i] = tr_z_yzz_yyyzz[i] * fe_0 + tr_z_xyzz_yyyzz[i] * pa_x[i];

        tr_z_xxyzz_yyzzz[i] = tr_z_yzz_yyzzz[i] * fe_0 + tr_z_xyzz_yyzzz[i] * pa_x[i];

        tr_z_xxyzz_yzzzz[i] = tr_z_yzz_yzzzz[i] * fe_0 + tr_z_xyzz_yzzzz[i] * pa_x[i];

        tr_z_xxyzz_zzzzz[i] = tr_z_xxzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1071-1092 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xxzzz_xxxxx, \
                             tr_z_xxzzz_xxxxy, \
                             tr_z_xxzzz_xxxxz, \
                             tr_z_xxzzz_xxxyy, \
                             tr_z_xxzzz_xxxyz, \
                             tr_z_xxzzz_xxxzz, \
                             tr_z_xxzzz_xxyyy, \
                             tr_z_xxzzz_xxyyz, \
                             tr_z_xxzzz_xxyzz, \
                             tr_z_xxzzz_xxzzz, \
                             tr_z_xxzzz_xyyyy, \
                             tr_z_xxzzz_xyyyz, \
                             tr_z_xxzzz_xyyzz, \
                             tr_z_xxzzz_xyzzz, \
                             tr_z_xxzzz_xzzzz, \
                             tr_z_xxzzz_yyyyy, \
                             tr_z_xxzzz_yyyyz, \
                             tr_z_xxzzz_yyyzz, \
                             tr_z_xxzzz_yyzzz, \
                             tr_z_xxzzz_yzzzz, \
                             tr_z_xxzzz_zzzzz, \
                             tr_z_xzzz_xxxx,   \
                             tr_z_xzzz_xxxxx,  \
                             tr_z_xzzz_xxxxy,  \
                             tr_z_xzzz_xxxxz,  \
                             tr_z_xzzz_xxxy,   \
                             tr_z_xzzz_xxxyy,  \
                             tr_z_xzzz_xxxyz,  \
                             tr_z_xzzz_xxxz,   \
                             tr_z_xzzz_xxxzz,  \
                             tr_z_xzzz_xxyy,   \
                             tr_z_xzzz_xxyyy,  \
                             tr_z_xzzz_xxyyz,  \
                             tr_z_xzzz_xxyz,   \
                             tr_z_xzzz_xxyzz,  \
                             tr_z_xzzz_xxzz,   \
                             tr_z_xzzz_xxzzz,  \
                             tr_z_xzzz_xyyy,   \
                             tr_z_xzzz_xyyyy,  \
                             tr_z_xzzz_xyyyz,  \
                             tr_z_xzzz_xyyz,   \
                             tr_z_xzzz_xyyzz,  \
                             tr_z_xzzz_xyzz,   \
                             tr_z_xzzz_xyzzz,  \
                             tr_z_xzzz_xzzz,   \
                             tr_z_xzzz_xzzzz,  \
                             tr_z_xzzz_yyyy,   \
                             tr_z_xzzz_yyyyy,  \
                             tr_z_xzzz_yyyyz,  \
                             tr_z_xzzz_yyyz,   \
                             tr_z_xzzz_yyyzz,  \
                             tr_z_xzzz_yyzz,   \
                             tr_z_xzzz_yyzzz,  \
                             tr_z_xzzz_yzzz,   \
                             tr_z_xzzz_yzzzz,  \
                             tr_z_xzzz_zzzz,   \
                             tr_z_xzzz_zzzzz,  \
                             tr_z_zzz_xxxxx,   \
                             tr_z_zzz_xxxxy,   \
                             tr_z_zzz_xxxxz,   \
                             tr_z_zzz_xxxyy,   \
                             tr_z_zzz_xxxyz,   \
                             tr_z_zzz_xxxzz,   \
                             tr_z_zzz_xxyyy,   \
                             tr_z_zzz_xxyyz,   \
                             tr_z_zzz_xxyzz,   \
                             tr_z_zzz_xxzzz,   \
                             tr_z_zzz_xyyyy,   \
                             tr_z_zzz_xyyyz,   \
                             tr_z_zzz_xyyzz,   \
                             tr_z_zzz_xyzzz,   \
                             tr_z_zzz_xzzzz,   \
                             tr_z_zzz_yyyyy,   \
                             tr_z_zzz_yyyyz,   \
                             tr_z_zzz_yyyzz,   \
                             tr_z_zzz_yyzzz,   \
                             tr_z_zzz_yzzzz,   \
                             tr_z_zzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzz_xxxxx[i] = tr_z_zzz_xxxxx[i] * fe_0 + 5.0 * tr_z_xzzz_xxxx[i] * fe_0 + tr_z_xzzz_xxxxx[i] * pa_x[i];

        tr_z_xxzzz_xxxxy[i] = tr_z_zzz_xxxxy[i] * fe_0 + 4.0 * tr_z_xzzz_xxxy[i] * fe_0 + tr_z_xzzz_xxxxy[i] * pa_x[i];

        tr_z_xxzzz_xxxxz[i] = tr_z_zzz_xxxxz[i] * fe_0 + 4.0 * tr_z_xzzz_xxxz[i] * fe_0 + tr_z_xzzz_xxxxz[i] * pa_x[i];

        tr_z_xxzzz_xxxyy[i] = tr_z_zzz_xxxyy[i] * fe_0 + 3.0 * tr_z_xzzz_xxyy[i] * fe_0 + tr_z_xzzz_xxxyy[i] * pa_x[i];

        tr_z_xxzzz_xxxyz[i] = tr_z_zzz_xxxyz[i] * fe_0 + 3.0 * tr_z_xzzz_xxyz[i] * fe_0 + tr_z_xzzz_xxxyz[i] * pa_x[i];

        tr_z_xxzzz_xxxzz[i] = tr_z_zzz_xxxzz[i] * fe_0 + 3.0 * tr_z_xzzz_xxzz[i] * fe_0 + tr_z_xzzz_xxxzz[i] * pa_x[i];

        tr_z_xxzzz_xxyyy[i] = tr_z_zzz_xxyyy[i] * fe_0 + 2.0 * tr_z_xzzz_xyyy[i] * fe_0 + tr_z_xzzz_xxyyy[i] * pa_x[i];

        tr_z_xxzzz_xxyyz[i] = tr_z_zzz_xxyyz[i] * fe_0 + 2.0 * tr_z_xzzz_xyyz[i] * fe_0 + tr_z_xzzz_xxyyz[i] * pa_x[i];

        tr_z_xxzzz_xxyzz[i] = tr_z_zzz_xxyzz[i] * fe_0 + 2.0 * tr_z_xzzz_xyzz[i] * fe_0 + tr_z_xzzz_xxyzz[i] * pa_x[i];

        tr_z_xxzzz_xxzzz[i] = tr_z_zzz_xxzzz[i] * fe_0 + 2.0 * tr_z_xzzz_xzzz[i] * fe_0 + tr_z_xzzz_xxzzz[i] * pa_x[i];

        tr_z_xxzzz_xyyyy[i] = tr_z_zzz_xyyyy[i] * fe_0 + tr_z_xzzz_yyyy[i] * fe_0 + tr_z_xzzz_xyyyy[i] * pa_x[i];

        tr_z_xxzzz_xyyyz[i] = tr_z_zzz_xyyyz[i] * fe_0 + tr_z_xzzz_yyyz[i] * fe_0 + tr_z_xzzz_xyyyz[i] * pa_x[i];

        tr_z_xxzzz_xyyzz[i] = tr_z_zzz_xyyzz[i] * fe_0 + tr_z_xzzz_yyzz[i] * fe_0 + tr_z_xzzz_xyyzz[i] * pa_x[i];

        tr_z_xxzzz_xyzzz[i] = tr_z_zzz_xyzzz[i] * fe_0 + tr_z_xzzz_yzzz[i] * fe_0 + tr_z_xzzz_xyzzz[i] * pa_x[i];

        tr_z_xxzzz_xzzzz[i] = tr_z_zzz_xzzzz[i] * fe_0 + tr_z_xzzz_zzzz[i] * fe_0 + tr_z_xzzz_xzzzz[i] * pa_x[i];

        tr_z_xxzzz_yyyyy[i] = tr_z_zzz_yyyyy[i] * fe_0 + tr_z_xzzz_yyyyy[i] * pa_x[i];

        tr_z_xxzzz_yyyyz[i] = tr_z_zzz_yyyyz[i] * fe_0 + tr_z_xzzz_yyyyz[i] * pa_x[i];

        tr_z_xxzzz_yyyzz[i] = tr_z_zzz_yyyzz[i] * fe_0 + tr_z_xzzz_yyyzz[i] * pa_x[i];

        tr_z_xxzzz_yyzzz[i] = tr_z_zzz_yyzzz[i] * fe_0 + tr_z_xzzz_yyzzz[i] * pa_x[i];

        tr_z_xxzzz_yzzzz[i] = tr_z_zzz_yzzzz[i] * fe_0 + tr_z_xzzz_yzzzz[i] * pa_x[i];

        tr_z_xxzzz_zzzzz[i] = tr_z_zzz_zzzzz[i] * fe_0 + tr_z_xzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1092-1113 components of targeted buffer : HH

    auto tr_z_xyyyy_xxxxx = pbuffer.data(idx_dip_hh + 1092);

    auto tr_z_xyyyy_xxxxy = pbuffer.data(idx_dip_hh + 1093);

    auto tr_z_xyyyy_xxxxz = pbuffer.data(idx_dip_hh + 1094);

    auto tr_z_xyyyy_xxxyy = pbuffer.data(idx_dip_hh + 1095);

    auto tr_z_xyyyy_xxxyz = pbuffer.data(idx_dip_hh + 1096);

    auto tr_z_xyyyy_xxxzz = pbuffer.data(idx_dip_hh + 1097);

    auto tr_z_xyyyy_xxyyy = pbuffer.data(idx_dip_hh + 1098);

    auto tr_z_xyyyy_xxyyz = pbuffer.data(idx_dip_hh + 1099);

    auto tr_z_xyyyy_xxyzz = pbuffer.data(idx_dip_hh + 1100);

    auto tr_z_xyyyy_xxzzz = pbuffer.data(idx_dip_hh + 1101);

    auto tr_z_xyyyy_xyyyy = pbuffer.data(idx_dip_hh + 1102);

    auto tr_z_xyyyy_xyyyz = pbuffer.data(idx_dip_hh + 1103);

    auto tr_z_xyyyy_xyyzz = pbuffer.data(idx_dip_hh + 1104);

    auto tr_z_xyyyy_xyzzz = pbuffer.data(idx_dip_hh + 1105);

    auto tr_z_xyyyy_xzzzz = pbuffer.data(idx_dip_hh + 1106);

    auto tr_z_xyyyy_yyyyy = pbuffer.data(idx_dip_hh + 1107);

    auto tr_z_xyyyy_yyyyz = pbuffer.data(idx_dip_hh + 1108);

    auto tr_z_xyyyy_yyyzz = pbuffer.data(idx_dip_hh + 1109);

    auto tr_z_xyyyy_yyzzz = pbuffer.data(idx_dip_hh + 1110);

    auto tr_z_xyyyy_yzzzz = pbuffer.data(idx_dip_hh + 1111);

    auto tr_z_xyyyy_zzzzz = pbuffer.data(idx_dip_hh + 1112);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyyy_xxxxx, \
                             tr_z_xyyyy_xxxxy, \
                             tr_z_xyyyy_xxxxz, \
                             tr_z_xyyyy_xxxyy, \
                             tr_z_xyyyy_xxxyz, \
                             tr_z_xyyyy_xxxzz, \
                             tr_z_xyyyy_xxyyy, \
                             tr_z_xyyyy_xxyyz, \
                             tr_z_xyyyy_xxyzz, \
                             tr_z_xyyyy_xxzzz, \
                             tr_z_xyyyy_xyyyy, \
                             tr_z_xyyyy_xyyyz, \
                             tr_z_xyyyy_xyyzz, \
                             tr_z_xyyyy_xyzzz, \
                             tr_z_xyyyy_xzzzz, \
                             tr_z_xyyyy_yyyyy, \
                             tr_z_xyyyy_yyyyz, \
                             tr_z_xyyyy_yyyzz, \
                             tr_z_xyyyy_yyzzz, \
                             tr_z_xyyyy_yzzzz, \
                             tr_z_xyyyy_zzzzz, \
                             tr_z_yyyy_xxxx,   \
                             tr_z_yyyy_xxxxx,  \
                             tr_z_yyyy_xxxxy,  \
                             tr_z_yyyy_xxxxz,  \
                             tr_z_yyyy_xxxy,   \
                             tr_z_yyyy_xxxyy,  \
                             tr_z_yyyy_xxxyz,  \
                             tr_z_yyyy_xxxz,   \
                             tr_z_yyyy_xxxzz,  \
                             tr_z_yyyy_xxyy,   \
                             tr_z_yyyy_xxyyy,  \
                             tr_z_yyyy_xxyyz,  \
                             tr_z_yyyy_xxyz,   \
                             tr_z_yyyy_xxyzz,  \
                             tr_z_yyyy_xxzz,   \
                             tr_z_yyyy_xxzzz,  \
                             tr_z_yyyy_xyyy,   \
                             tr_z_yyyy_xyyyy,  \
                             tr_z_yyyy_xyyyz,  \
                             tr_z_yyyy_xyyz,   \
                             tr_z_yyyy_xyyzz,  \
                             tr_z_yyyy_xyzz,   \
                             tr_z_yyyy_xyzzz,  \
                             tr_z_yyyy_xzzz,   \
                             tr_z_yyyy_xzzzz,  \
                             tr_z_yyyy_yyyy,   \
                             tr_z_yyyy_yyyyy,  \
                             tr_z_yyyy_yyyyz,  \
                             tr_z_yyyy_yyyz,   \
                             tr_z_yyyy_yyyzz,  \
                             tr_z_yyyy_yyzz,   \
                             tr_z_yyyy_yyzzz,  \
                             tr_z_yyyy_yzzz,   \
                             tr_z_yyyy_yzzzz,  \
                             tr_z_yyyy_zzzz,   \
                             tr_z_yyyy_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyy_xxxxx[i] = 5.0 * tr_z_yyyy_xxxx[i] * fe_0 + tr_z_yyyy_xxxxx[i] * pa_x[i];

        tr_z_xyyyy_xxxxy[i] = 4.0 * tr_z_yyyy_xxxy[i] * fe_0 + tr_z_yyyy_xxxxy[i] * pa_x[i];

        tr_z_xyyyy_xxxxz[i] = 4.0 * tr_z_yyyy_xxxz[i] * fe_0 + tr_z_yyyy_xxxxz[i] * pa_x[i];

        tr_z_xyyyy_xxxyy[i] = 3.0 * tr_z_yyyy_xxyy[i] * fe_0 + tr_z_yyyy_xxxyy[i] * pa_x[i];

        tr_z_xyyyy_xxxyz[i] = 3.0 * tr_z_yyyy_xxyz[i] * fe_0 + tr_z_yyyy_xxxyz[i] * pa_x[i];

        tr_z_xyyyy_xxxzz[i] = 3.0 * tr_z_yyyy_xxzz[i] * fe_0 + tr_z_yyyy_xxxzz[i] * pa_x[i];

        tr_z_xyyyy_xxyyy[i] = 2.0 * tr_z_yyyy_xyyy[i] * fe_0 + tr_z_yyyy_xxyyy[i] * pa_x[i];

        tr_z_xyyyy_xxyyz[i] = 2.0 * tr_z_yyyy_xyyz[i] * fe_0 + tr_z_yyyy_xxyyz[i] * pa_x[i];

        tr_z_xyyyy_xxyzz[i] = 2.0 * tr_z_yyyy_xyzz[i] * fe_0 + tr_z_yyyy_xxyzz[i] * pa_x[i];

        tr_z_xyyyy_xxzzz[i] = 2.0 * tr_z_yyyy_xzzz[i] * fe_0 + tr_z_yyyy_xxzzz[i] * pa_x[i];

        tr_z_xyyyy_xyyyy[i] = tr_z_yyyy_yyyy[i] * fe_0 + tr_z_yyyy_xyyyy[i] * pa_x[i];

        tr_z_xyyyy_xyyyz[i] = tr_z_yyyy_yyyz[i] * fe_0 + tr_z_yyyy_xyyyz[i] * pa_x[i];

        tr_z_xyyyy_xyyzz[i] = tr_z_yyyy_yyzz[i] * fe_0 + tr_z_yyyy_xyyzz[i] * pa_x[i];

        tr_z_xyyyy_xyzzz[i] = tr_z_yyyy_yzzz[i] * fe_0 + tr_z_yyyy_xyzzz[i] * pa_x[i];

        tr_z_xyyyy_xzzzz[i] = tr_z_yyyy_zzzz[i] * fe_0 + tr_z_yyyy_xzzzz[i] * pa_x[i];

        tr_z_xyyyy_yyyyy[i] = tr_z_yyyy_yyyyy[i] * pa_x[i];

        tr_z_xyyyy_yyyyz[i] = tr_z_yyyy_yyyyz[i] * pa_x[i];

        tr_z_xyyyy_yyyzz[i] = tr_z_yyyy_yyyzz[i] * pa_x[i];

        tr_z_xyyyy_yyzzz[i] = tr_z_yyyy_yyzzz[i] * pa_x[i];

        tr_z_xyyyy_yzzzz[i] = tr_z_yyyy_yzzzz[i] * pa_x[i];

        tr_z_xyyyy_zzzzz[i] = tr_z_yyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 1113-1134 components of targeted buffer : HH

    auto tr_z_xyyyz_xxxxx = pbuffer.data(idx_dip_hh + 1113);

    auto tr_z_xyyyz_xxxxy = pbuffer.data(idx_dip_hh + 1114);

    auto tr_z_xyyyz_xxxxz = pbuffer.data(idx_dip_hh + 1115);

    auto tr_z_xyyyz_xxxyy = pbuffer.data(idx_dip_hh + 1116);

    auto tr_z_xyyyz_xxxyz = pbuffer.data(idx_dip_hh + 1117);

    auto tr_z_xyyyz_xxxzz = pbuffer.data(idx_dip_hh + 1118);

    auto tr_z_xyyyz_xxyyy = pbuffer.data(idx_dip_hh + 1119);

    auto tr_z_xyyyz_xxyyz = pbuffer.data(idx_dip_hh + 1120);

    auto tr_z_xyyyz_xxyzz = pbuffer.data(idx_dip_hh + 1121);

    auto tr_z_xyyyz_xxzzz = pbuffer.data(idx_dip_hh + 1122);

    auto tr_z_xyyyz_xyyyy = pbuffer.data(idx_dip_hh + 1123);

    auto tr_z_xyyyz_xyyyz = pbuffer.data(idx_dip_hh + 1124);

    auto tr_z_xyyyz_xyyzz = pbuffer.data(idx_dip_hh + 1125);

    auto tr_z_xyyyz_xyzzz = pbuffer.data(idx_dip_hh + 1126);

    auto tr_z_xyyyz_xzzzz = pbuffer.data(idx_dip_hh + 1127);

    auto tr_z_xyyyz_yyyyy = pbuffer.data(idx_dip_hh + 1128);

    auto tr_z_xyyyz_yyyyz = pbuffer.data(idx_dip_hh + 1129);

    auto tr_z_xyyyz_yyyzz = pbuffer.data(idx_dip_hh + 1130);

    auto tr_z_xyyyz_yyzzz = pbuffer.data(idx_dip_hh + 1131);

    auto tr_z_xyyyz_yzzzz = pbuffer.data(idx_dip_hh + 1132);

    auto tr_z_xyyyz_zzzzz = pbuffer.data(idx_dip_hh + 1133);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyyz_xxxxx, \
                             tr_z_xyyyz_xxxxy, \
                             tr_z_xyyyz_xxxxz, \
                             tr_z_xyyyz_xxxyy, \
                             tr_z_xyyyz_xxxyz, \
                             tr_z_xyyyz_xxxzz, \
                             tr_z_xyyyz_xxyyy, \
                             tr_z_xyyyz_xxyyz, \
                             tr_z_xyyyz_xxyzz, \
                             tr_z_xyyyz_xxzzz, \
                             tr_z_xyyyz_xyyyy, \
                             tr_z_xyyyz_xyyyz, \
                             tr_z_xyyyz_xyyzz, \
                             tr_z_xyyyz_xyzzz, \
                             tr_z_xyyyz_xzzzz, \
                             tr_z_xyyyz_yyyyy, \
                             tr_z_xyyyz_yyyyz, \
                             tr_z_xyyyz_yyyzz, \
                             tr_z_xyyyz_yyzzz, \
                             tr_z_xyyyz_yzzzz, \
                             tr_z_xyyyz_zzzzz, \
                             tr_z_yyyz_xxxx,   \
                             tr_z_yyyz_xxxxx,  \
                             tr_z_yyyz_xxxxy,  \
                             tr_z_yyyz_xxxxz,  \
                             tr_z_yyyz_xxxy,   \
                             tr_z_yyyz_xxxyy,  \
                             tr_z_yyyz_xxxyz,  \
                             tr_z_yyyz_xxxz,   \
                             tr_z_yyyz_xxxzz,  \
                             tr_z_yyyz_xxyy,   \
                             tr_z_yyyz_xxyyy,  \
                             tr_z_yyyz_xxyyz,  \
                             tr_z_yyyz_xxyz,   \
                             tr_z_yyyz_xxyzz,  \
                             tr_z_yyyz_xxzz,   \
                             tr_z_yyyz_xxzzz,  \
                             tr_z_yyyz_xyyy,   \
                             tr_z_yyyz_xyyyy,  \
                             tr_z_yyyz_xyyyz,  \
                             tr_z_yyyz_xyyz,   \
                             tr_z_yyyz_xyyzz,  \
                             tr_z_yyyz_xyzz,   \
                             tr_z_yyyz_xyzzz,  \
                             tr_z_yyyz_xzzz,   \
                             tr_z_yyyz_xzzzz,  \
                             tr_z_yyyz_yyyy,   \
                             tr_z_yyyz_yyyyy,  \
                             tr_z_yyyz_yyyyz,  \
                             tr_z_yyyz_yyyz,   \
                             tr_z_yyyz_yyyzz,  \
                             tr_z_yyyz_yyzz,   \
                             tr_z_yyyz_yyzzz,  \
                             tr_z_yyyz_yzzz,   \
                             tr_z_yyyz_yzzzz,  \
                             tr_z_yyyz_zzzz,   \
                             tr_z_yyyz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyz_xxxxx[i] = 5.0 * tr_z_yyyz_xxxx[i] * fe_0 + tr_z_yyyz_xxxxx[i] * pa_x[i];

        tr_z_xyyyz_xxxxy[i] = 4.0 * tr_z_yyyz_xxxy[i] * fe_0 + tr_z_yyyz_xxxxy[i] * pa_x[i];

        tr_z_xyyyz_xxxxz[i] = 4.0 * tr_z_yyyz_xxxz[i] * fe_0 + tr_z_yyyz_xxxxz[i] * pa_x[i];

        tr_z_xyyyz_xxxyy[i] = 3.0 * tr_z_yyyz_xxyy[i] * fe_0 + tr_z_yyyz_xxxyy[i] * pa_x[i];

        tr_z_xyyyz_xxxyz[i] = 3.0 * tr_z_yyyz_xxyz[i] * fe_0 + tr_z_yyyz_xxxyz[i] * pa_x[i];

        tr_z_xyyyz_xxxzz[i] = 3.0 * tr_z_yyyz_xxzz[i] * fe_0 + tr_z_yyyz_xxxzz[i] * pa_x[i];

        tr_z_xyyyz_xxyyy[i] = 2.0 * tr_z_yyyz_xyyy[i] * fe_0 + tr_z_yyyz_xxyyy[i] * pa_x[i];

        tr_z_xyyyz_xxyyz[i] = 2.0 * tr_z_yyyz_xyyz[i] * fe_0 + tr_z_yyyz_xxyyz[i] * pa_x[i];

        tr_z_xyyyz_xxyzz[i] = 2.0 * tr_z_yyyz_xyzz[i] * fe_0 + tr_z_yyyz_xxyzz[i] * pa_x[i];

        tr_z_xyyyz_xxzzz[i] = 2.0 * tr_z_yyyz_xzzz[i] * fe_0 + tr_z_yyyz_xxzzz[i] * pa_x[i];

        tr_z_xyyyz_xyyyy[i] = tr_z_yyyz_yyyy[i] * fe_0 + tr_z_yyyz_xyyyy[i] * pa_x[i];

        tr_z_xyyyz_xyyyz[i] = tr_z_yyyz_yyyz[i] * fe_0 + tr_z_yyyz_xyyyz[i] * pa_x[i];

        tr_z_xyyyz_xyyzz[i] = tr_z_yyyz_yyzz[i] * fe_0 + tr_z_yyyz_xyyzz[i] * pa_x[i];

        tr_z_xyyyz_xyzzz[i] = tr_z_yyyz_yzzz[i] * fe_0 + tr_z_yyyz_xyzzz[i] * pa_x[i];

        tr_z_xyyyz_xzzzz[i] = tr_z_yyyz_zzzz[i] * fe_0 + tr_z_yyyz_xzzzz[i] * pa_x[i];

        tr_z_xyyyz_yyyyy[i] = tr_z_yyyz_yyyyy[i] * pa_x[i];

        tr_z_xyyyz_yyyyz[i] = tr_z_yyyz_yyyyz[i] * pa_x[i];

        tr_z_xyyyz_yyyzz[i] = tr_z_yyyz_yyyzz[i] * pa_x[i];

        tr_z_xyyyz_yyzzz[i] = tr_z_yyyz_yyzzz[i] * pa_x[i];

        tr_z_xyyyz_yzzzz[i] = tr_z_yyyz_yzzzz[i] * pa_x[i];

        tr_z_xyyyz_zzzzz[i] = tr_z_yyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 1134-1155 components of targeted buffer : HH

    auto tr_z_xyyzz_xxxxx = pbuffer.data(idx_dip_hh + 1134);

    auto tr_z_xyyzz_xxxxy = pbuffer.data(idx_dip_hh + 1135);

    auto tr_z_xyyzz_xxxxz = pbuffer.data(idx_dip_hh + 1136);

    auto tr_z_xyyzz_xxxyy = pbuffer.data(idx_dip_hh + 1137);

    auto tr_z_xyyzz_xxxyz = pbuffer.data(idx_dip_hh + 1138);

    auto tr_z_xyyzz_xxxzz = pbuffer.data(idx_dip_hh + 1139);

    auto tr_z_xyyzz_xxyyy = pbuffer.data(idx_dip_hh + 1140);

    auto tr_z_xyyzz_xxyyz = pbuffer.data(idx_dip_hh + 1141);

    auto tr_z_xyyzz_xxyzz = pbuffer.data(idx_dip_hh + 1142);

    auto tr_z_xyyzz_xxzzz = pbuffer.data(idx_dip_hh + 1143);

    auto tr_z_xyyzz_xyyyy = pbuffer.data(idx_dip_hh + 1144);

    auto tr_z_xyyzz_xyyyz = pbuffer.data(idx_dip_hh + 1145);

    auto tr_z_xyyzz_xyyzz = pbuffer.data(idx_dip_hh + 1146);

    auto tr_z_xyyzz_xyzzz = pbuffer.data(idx_dip_hh + 1147);

    auto tr_z_xyyzz_xzzzz = pbuffer.data(idx_dip_hh + 1148);

    auto tr_z_xyyzz_yyyyy = pbuffer.data(idx_dip_hh + 1149);

    auto tr_z_xyyzz_yyyyz = pbuffer.data(idx_dip_hh + 1150);

    auto tr_z_xyyzz_yyyzz = pbuffer.data(idx_dip_hh + 1151);

    auto tr_z_xyyzz_yyzzz = pbuffer.data(idx_dip_hh + 1152);

    auto tr_z_xyyzz_yzzzz = pbuffer.data(idx_dip_hh + 1153);

    auto tr_z_xyyzz_zzzzz = pbuffer.data(idx_dip_hh + 1154);

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xyyzz_xxxxx, \
                             tr_z_xyyzz_xxxxy, \
                             tr_z_xyyzz_xxxxz, \
                             tr_z_xyyzz_xxxyy, \
                             tr_z_xyyzz_xxxyz, \
                             tr_z_xyyzz_xxxzz, \
                             tr_z_xyyzz_xxyyy, \
                             tr_z_xyyzz_xxyyz, \
                             tr_z_xyyzz_xxyzz, \
                             tr_z_xyyzz_xxzzz, \
                             tr_z_xyyzz_xyyyy, \
                             tr_z_xyyzz_xyyyz, \
                             tr_z_xyyzz_xyyzz, \
                             tr_z_xyyzz_xyzzz, \
                             tr_z_xyyzz_xzzzz, \
                             tr_z_xyyzz_yyyyy, \
                             tr_z_xyyzz_yyyyz, \
                             tr_z_xyyzz_yyyzz, \
                             tr_z_xyyzz_yyzzz, \
                             tr_z_xyyzz_yzzzz, \
                             tr_z_xyyzz_zzzzz, \
                             tr_z_yyzz_xxxx,   \
                             tr_z_yyzz_xxxxx,  \
                             tr_z_yyzz_xxxxy,  \
                             tr_z_yyzz_xxxxz,  \
                             tr_z_yyzz_xxxy,   \
                             tr_z_yyzz_xxxyy,  \
                             tr_z_yyzz_xxxyz,  \
                             tr_z_yyzz_xxxz,   \
                             tr_z_yyzz_xxxzz,  \
                             tr_z_yyzz_xxyy,   \
                             tr_z_yyzz_xxyyy,  \
                             tr_z_yyzz_xxyyz,  \
                             tr_z_yyzz_xxyz,   \
                             tr_z_yyzz_xxyzz,  \
                             tr_z_yyzz_xxzz,   \
                             tr_z_yyzz_xxzzz,  \
                             tr_z_yyzz_xyyy,   \
                             tr_z_yyzz_xyyyy,  \
                             tr_z_yyzz_xyyyz,  \
                             tr_z_yyzz_xyyz,   \
                             tr_z_yyzz_xyyzz,  \
                             tr_z_yyzz_xyzz,   \
                             tr_z_yyzz_xyzzz,  \
                             tr_z_yyzz_xzzz,   \
                             tr_z_yyzz_xzzzz,  \
                             tr_z_yyzz_yyyy,   \
                             tr_z_yyzz_yyyyy,  \
                             tr_z_yyzz_yyyyz,  \
                             tr_z_yyzz_yyyz,   \
                             tr_z_yyzz_yyyzz,  \
                             tr_z_yyzz_yyzz,   \
                             tr_z_yyzz_yyzzz,  \
                             tr_z_yyzz_yzzz,   \
                             tr_z_yyzz_yzzzz,  \
                             tr_z_yyzz_zzzz,   \
                             tr_z_yyzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzz_xxxxx[i] = 5.0 * tr_z_yyzz_xxxx[i] * fe_0 + tr_z_yyzz_xxxxx[i] * pa_x[i];

        tr_z_xyyzz_xxxxy[i] = 4.0 * tr_z_yyzz_xxxy[i] * fe_0 + tr_z_yyzz_xxxxy[i] * pa_x[i];

        tr_z_xyyzz_xxxxz[i] = 4.0 * tr_z_yyzz_xxxz[i] * fe_0 + tr_z_yyzz_xxxxz[i] * pa_x[i];

        tr_z_xyyzz_xxxyy[i] = 3.0 * tr_z_yyzz_xxyy[i] * fe_0 + tr_z_yyzz_xxxyy[i] * pa_x[i];

        tr_z_xyyzz_xxxyz[i] = 3.0 * tr_z_yyzz_xxyz[i] * fe_0 + tr_z_yyzz_xxxyz[i] * pa_x[i];

        tr_z_xyyzz_xxxzz[i] = 3.0 * tr_z_yyzz_xxzz[i] * fe_0 + tr_z_yyzz_xxxzz[i] * pa_x[i];

        tr_z_xyyzz_xxyyy[i] = 2.0 * tr_z_yyzz_xyyy[i] * fe_0 + tr_z_yyzz_xxyyy[i] * pa_x[i];

        tr_z_xyyzz_xxyyz[i] = 2.0 * tr_z_yyzz_xyyz[i] * fe_0 + tr_z_yyzz_xxyyz[i] * pa_x[i];

        tr_z_xyyzz_xxyzz[i] = 2.0 * tr_z_yyzz_xyzz[i] * fe_0 + tr_z_yyzz_xxyzz[i] * pa_x[i];

        tr_z_xyyzz_xxzzz[i] = 2.0 * tr_z_yyzz_xzzz[i] * fe_0 + tr_z_yyzz_xxzzz[i] * pa_x[i];

        tr_z_xyyzz_xyyyy[i] = tr_z_yyzz_yyyy[i] * fe_0 + tr_z_yyzz_xyyyy[i] * pa_x[i];

        tr_z_xyyzz_xyyyz[i] = tr_z_yyzz_yyyz[i] * fe_0 + tr_z_yyzz_xyyyz[i] * pa_x[i];

        tr_z_xyyzz_xyyzz[i] = tr_z_yyzz_yyzz[i] * fe_0 + tr_z_yyzz_xyyzz[i] * pa_x[i];

        tr_z_xyyzz_xyzzz[i] = tr_z_yyzz_yzzz[i] * fe_0 + tr_z_yyzz_xyzzz[i] * pa_x[i];

        tr_z_xyyzz_xzzzz[i] = tr_z_yyzz_zzzz[i] * fe_0 + tr_z_yyzz_xzzzz[i] * pa_x[i];

        tr_z_xyyzz_yyyyy[i] = tr_z_yyzz_yyyyy[i] * pa_x[i];

        tr_z_xyyzz_yyyyz[i] = tr_z_yyzz_yyyyz[i] * pa_x[i];

        tr_z_xyyzz_yyyzz[i] = tr_z_yyzz_yyyzz[i] * pa_x[i];

        tr_z_xyyzz_yyzzz[i] = tr_z_yyzz_yyzzz[i] * pa_x[i];

        tr_z_xyyzz_yzzzz[i] = tr_z_yyzz_yzzzz[i] * pa_x[i];

        tr_z_xyyzz_zzzzz[i] = tr_z_yyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1155-1176 components of targeted buffer : HH

    auto tr_z_xyzzz_xxxxx = pbuffer.data(idx_dip_hh + 1155);

    auto tr_z_xyzzz_xxxxy = pbuffer.data(idx_dip_hh + 1156);

    auto tr_z_xyzzz_xxxxz = pbuffer.data(idx_dip_hh + 1157);

    auto tr_z_xyzzz_xxxyy = pbuffer.data(idx_dip_hh + 1158);

    auto tr_z_xyzzz_xxxyz = pbuffer.data(idx_dip_hh + 1159);

    auto tr_z_xyzzz_xxxzz = pbuffer.data(idx_dip_hh + 1160);

    auto tr_z_xyzzz_xxyyy = pbuffer.data(idx_dip_hh + 1161);

    auto tr_z_xyzzz_xxyyz = pbuffer.data(idx_dip_hh + 1162);

    auto tr_z_xyzzz_xxyzz = pbuffer.data(idx_dip_hh + 1163);

    auto tr_z_xyzzz_xxzzz = pbuffer.data(idx_dip_hh + 1164);

    auto tr_z_xyzzz_xyyyy = pbuffer.data(idx_dip_hh + 1165);

    auto tr_z_xyzzz_xyyyz = pbuffer.data(idx_dip_hh + 1166);

    auto tr_z_xyzzz_xyyzz = pbuffer.data(idx_dip_hh + 1167);

    auto tr_z_xyzzz_xyzzz = pbuffer.data(idx_dip_hh + 1168);

    auto tr_z_xyzzz_xzzzz = pbuffer.data(idx_dip_hh + 1169);

    auto tr_z_xyzzz_yyyyy = pbuffer.data(idx_dip_hh + 1170);

    auto tr_z_xyzzz_yyyyz = pbuffer.data(idx_dip_hh + 1171);

    auto tr_z_xyzzz_yyyzz = pbuffer.data(idx_dip_hh + 1172);

    auto tr_z_xyzzz_yyzzz = pbuffer.data(idx_dip_hh + 1173);

    auto tr_z_xyzzz_yzzzz = pbuffer.data(idx_dip_hh + 1174);

    auto tr_z_xyzzz_zzzzz = pbuffer.data(idx_dip_hh + 1175);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             tr_z_xyzzz_xxxxx, \
                             tr_z_xyzzz_xxxxy, \
                             tr_z_xyzzz_xxxxz, \
                             tr_z_xyzzz_xxxyy, \
                             tr_z_xyzzz_xxxyz, \
                             tr_z_xyzzz_xxxzz, \
                             tr_z_xyzzz_xxyyy, \
                             tr_z_xyzzz_xxyyz, \
                             tr_z_xyzzz_xxyzz, \
                             tr_z_xyzzz_xxzzz, \
                             tr_z_xyzzz_xyyyy, \
                             tr_z_xyzzz_xyyyz, \
                             tr_z_xyzzz_xyyzz, \
                             tr_z_xyzzz_xyzzz, \
                             tr_z_xyzzz_xzzzz, \
                             tr_z_xyzzz_yyyyy, \
                             tr_z_xyzzz_yyyyz, \
                             tr_z_xyzzz_yyyzz, \
                             tr_z_xyzzz_yyzzz, \
                             tr_z_xyzzz_yzzzz, \
                             tr_z_xyzzz_zzzzz, \
                             tr_z_xzzz_xxxxx,  \
                             tr_z_xzzz_xxxxz,  \
                             tr_z_xzzz_xxxzz,  \
                             tr_z_xzzz_xxzzz,  \
                             tr_z_xzzz_xzzzz,  \
                             tr_z_yzzz_xxxxy,  \
                             tr_z_yzzz_xxxy,   \
                             tr_z_yzzz_xxxyy,  \
                             tr_z_yzzz_xxxyz,  \
                             tr_z_yzzz_xxyy,   \
                             tr_z_yzzz_xxyyy,  \
                             tr_z_yzzz_xxyyz,  \
                             tr_z_yzzz_xxyz,   \
                             tr_z_yzzz_xxyzz,  \
                             tr_z_yzzz_xyyy,   \
                             tr_z_yzzz_xyyyy,  \
                             tr_z_yzzz_xyyyz,  \
                             tr_z_yzzz_xyyz,   \
                             tr_z_yzzz_xyyzz,  \
                             tr_z_yzzz_xyzz,   \
                             tr_z_yzzz_xyzzz,  \
                             tr_z_yzzz_yyyy,   \
                             tr_z_yzzz_yyyyy,  \
                             tr_z_yzzz_yyyyz,  \
                             tr_z_yzzz_yyyz,   \
                             tr_z_yzzz_yyyzz,  \
                             tr_z_yzzz_yyzz,   \
                             tr_z_yzzz_yyzzz,  \
                             tr_z_yzzz_yzzz,   \
                             tr_z_yzzz_yzzzz,  \
                             tr_z_yzzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzz_xxxxx[i] = tr_z_xzzz_xxxxx[i] * pa_y[i];

        tr_z_xyzzz_xxxxy[i] = 4.0 * tr_z_yzzz_xxxy[i] * fe_0 + tr_z_yzzz_xxxxy[i] * pa_x[i];

        tr_z_xyzzz_xxxxz[i] = tr_z_xzzz_xxxxz[i] * pa_y[i];

        tr_z_xyzzz_xxxyy[i] = 3.0 * tr_z_yzzz_xxyy[i] * fe_0 + tr_z_yzzz_xxxyy[i] * pa_x[i];

        tr_z_xyzzz_xxxyz[i] = 3.0 * tr_z_yzzz_xxyz[i] * fe_0 + tr_z_yzzz_xxxyz[i] * pa_x[i];

        tr_z_xyzzz_xxxzz[i] = tr_z_xzzz_xxxzz[i] * pa_y[i];

        tr_z_xyzzz_xxyyy[i] = 2.0 * tr_z_yzzz_xyyy[i] * fe_0 + tr_z_yzzz_xxyyy[i] * pa_x[i];

        tr_z_xyzzz_xxyyz[i] = 2.0 * tr_z_yzzz_xyyz[i] * fe_0 + tr_z_yzzz_xxyyz[i] * pa_x[i];

        tr_z_xyzzz_xxyzz[i] = 2.0 * tr_z_yzzz_xyzz[i] * fe_0 + tr_z_yzzz_xxyzz[i] * pa_x[i];

        tr_z_xyzzz_xxzzz[i] = tr_z_xzzz_xxzzz[i] * pa_y[i];

        tr_z_xyzzz_xyyyy[i] = tr_z_yzzz_yyyy[i] * fe_0 + tr_z_yzzz_xyyyy[i] * pa_x[i];

        tr_z_xyzzz_xyyyz[i] = tr_z_yzzz_yyyz[i] * fe_0 + tr_z_yzzz_xyyyz[i] * pa_x[i];

        tr_z_xyzzz_xyyzz[i] = tr_z_yzzz_yyzz[i] * fe_0 + tr_z_yzzz_xyyzz[i] * pa_x[i];

        tr_z_xyzzz_xyzzz[i] = tr_z_yzzz_yzzz[i] * fe_0 + tr_z_yzzz_xyzzz[i] * pa_x[i];

        tr_z_xyzzz_xzzzz[i] = tr_z_xzzz_xzzzz[i] * pa_y[i];

        tr_z_xyzzz_yyyyy[i] = tr_z_yzzz_yyyyy[i] * pa_x[i];

        tr_z_xyzzz_yyyyz[i] = tr_z_yzzz_yyyyz[i] * pa_x[i];

        tr_z_xyzzz_yyyzz[i] = tr_z_yzzz_yyyzz[i] * pa_x[i];

        tr_z_xyzzz_yyzzz[i] = tr_z_yzzz_yyzzz[i] * pa_x[i];

        tr_z_xyzzz_yzzzz[i] = tr_z_yzzz_yzzzz[i] * pa_x[i];

        tr_z_xyzzz_zzzzz[i] = tr_z_yzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1176-1197 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_x,                 \
                             tr_z_xzzzz_xxxxx, \
                             tr_z_xzzzz_xxxxy, \
                             tr_z_xzzzz_xxxxz, \
                             tr_z_xzzzz_xxxyy, \
                             tr_z_xzzzz_xxxyz, \
                             tr_z_xzzzz_xxxzz, \
                             tr_z_xzzzz_xxyyy, \
                             tr_z_xzzzz_xxyyz, \
                             tr_z_xzzzz_xxyzz, \
                             tr_z_xzzzz_xxzzz, \
                             tr_z_xzzzz_xyyyy, \
                             tr_z_xzzzz_xyyyz, \
                             tr_z_xzzzz_xyyzz, \
                             tr_z_xzzzz_xyzzz, \
                             tr_z_xzzzz_xzzzz, \
                             tr_z_xzzzz_yyyyy, \
                             tr_z_xzzzz_yyyyz, \
                             tr_z_xzzzz_yyyzz, \
                             tr_z_xzzzz_yyzzz, \
                             tr_z_xzzzz_yzzzz, \
                             tr_z_xzzzz_zzzzz, \
                             tr_z_zzzz_xxxx,   \
                             tr_z_zzzz_xxxxx,  \
                             tr_z_zzzz_xxxxy,  \
                             tr_z_zzzz_xxxxz,  \
                             tr_z_zzzz_xxxy,   \
                             tr_z_zzzz_xxxyy,  \
                             tr_z_zzzz_xxxyz,  \
                             tr_z_zzzz_xxxz,   \
                             tr_z_zzzz_xxxzz,  \
                             tr_z_zzzz_xxyy,   \
                             tr_z_zzzz_xxyyy,  \
                             tr_z_zzzz_xxyyz,  \
                             tr_z_zzzz_xxyz,   \
                             tr_z_zzzz_xxyzz,  \
                             tr_z_zzzz_xxzz,   \
                             tr_z_zzzz_xxzzz,  \
                             tr_z_zzzz_xyyy,   \
                             tr_z_zzzz_xyyyy,  \
                             tr_z_zzzz_xyyyz,  \
                             tr_z_zzzz_xyyz,   \
                             tr_z_zzzz_xyyzz,  \
                             tr_z_zzzz_xyzz,   \
                             tr_z_zzzz_xyzzz,  \
                             tr_z_zzzz_xzzz,   \
                             tr_z_zzzz_xzzzz,  \
                             tr_z_zzzz_yyyy,   \
                             tr_z_zzzz_yyyyy,  \
                             tr_z_zzzz_yyyyz,  \
                             tr_z_zzzz_yyyz,   \
                             tr_z_zzzz_yyyzz,  \
                             tr_z_zzzz_yyzz,   \
                             tr_z_zzzz_yyzzz,  \
                             tr_z_zzzz_yzzz,   \
                             tr_z_zzzz_yzzzz,  \
                             tr_z_zzzz_zzzz,   \
                             tr_z_zzzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzz_xxxxx[i] = 5.0 * tr_z_zzzz_xxxx[i] * fe_0 + tr_z_zzzz_xxxxx[i] * pa_x[i];

        tr_z_xzzzz_xxxxy[i] = 4.0 * tr_z_zzzz_xxxy[i] * fe_0 + tr_z_zzzz_xxxxy[i] * pa_x[i];

        tr_z_xzzzz_xxxxz[i] = 4.0 * tr_z_zzzz_xxxz[i] * fe_0 + tr_z_zzzz_xxxxz[i] * pa_x[i];

        tr_z_xzzzz_xxxyy[i] = 3.0 * tr_z_zzzz_xxyy[i] * fe_0 + tr_z_zzzz_xxxyy[i] * pa_x[i];

        tr_z_xzzzz_xxxyz[i] = 3.0 * tr_z_zzzz_xxyz[i] * fe_0 + tr_z_zzzz_xxxyz[i] * pa_x[i];

        tr_z_xzzzz_xxxzz[i] = 3.0 * tr_z_zzzz_xxzz[i] * fe_0 + tr_z_zzzz_xxxzz[i] * pa_x[i];

        tr_z_xzzzz_xxyyy[i] = 2.0 * tr_z_zzzz_xyyy[i] * fe_0 + tr_z_zzzz_xxyyy[i] * pa_x[i];

        tr_z_xzzzz_xxyyz[i] = 2.0 * tr_z_zzzz_xyyz[i] * fe_0 + tr_z_zzzz_xxyyz[i] * pa_x[i];

        tr_z_xzzzz_xxyzz[i] = 2.0 * tr_z_zzzz_xyzz[i] * fe_0 + tr_z_zzzz_xxyzz[i] * pa_x[i];

        tr_z_xzzzz_xxzzz[i] = 2.0 * tr_z_zzzz_xzzz[i] * fe_0 + tr_z_zzzz_xxzzz[i] * pa_x[i];

        tr_z_xzzzz_xyyyy[i] = tr_z_zzzz_yyyy[i] * fe_0 + tr_z_zzzz_xyyyy[i] * pa_x[i];

        tr_z_xzzzz_xyyyz[i] = tr_z_zzzz_yyyz[i] * fe_0 + tr_z_zzzz_xyyyz[i] * pa_x[i];

        tr_z_xzzzz_xyyzz[i] = tr_z_zzzz_yyzz[i] * fe_0 + tr_z_zzzz_xyyzz[i] * pa_x[i];

        tr_z_xzzzz_xyzzz[i] = tr_z_zzzz_yzzz[i] * fe_0 + tr_z_zzzz_xyzzz[i] * pa_x[i];

        tr_z_xzzzz_xzzzz[i] = tr_z_zzzz_zzzz[i] * fe_0 + tr_z_zzzz_xzzzz[i] * pa_x[i];

        tr_z_xzzzz_yyyyy[i] = tr_z_zzzz_yyyyy[i] * pa_x[i];

        tr_z_xzzzz_yyyyz[i] = tr_z_zzzz_yyyyz[i] * pa_x[i];

        tr_z_xzzzz_yyyzz[i] = tr_z_zzzz_yyyzz[i] * pa_x[i];

        tr_z_xzzzz_yyzzz[i] = tr_z_zzzz_yyzzz[i] * pa_x[i];

        tr_z_xzzzz_yzzzz[i] = tr_z_zzzz_yzzzz[i] * pa_x[i];

        tr_z_xzzzz_zzzzz[i] = tr_z_zzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 1197-1218 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyy_xxxxx,   \
                             tr_z_yyy_xxxxy,   \
                             tr_z_yyy_xxxxz,   \
                             tr_z_yyy_xxxyy,   \
                             tr_z_yyy_xxxyz,   \
                             tr_z_yyy_xxxzz,   \
                             tr_z_yyy_xxyyy,   \
                             tr_z_yyy_xxyyz,   \
                             tr_z_yyy_xxyzz,   \
                             tr_z_yyy_xxzzz,   \
                             tr_z_yyy_xyyyy,   \
                             tr_z_yyy_xyyyz,   \
                             tr_z_yyy_xyyzz,   \
                             tr_z_yyy_xyzzz,   \
                             tr_z_yyy_xzzzz,   \
                             tr_z_yyy_yyyyy,   \
                             tr_z_yyy_yyyyz,   \
                             tr_z_yyy_yyyzz,   \
                             tr_z_yyy_yyzzz,   \
                             tr_z_yyy_yzzzz,   \
                             tr_z_yyy_zzzzz,   \
                             tr_z_yyyy_xxxx,   \
                             tr_z_yyyy_xxxxx,  \
                             tr_z_yyyy_xxxxy,  \
                             tr_z_yyyy_xxxxz,  \
                             tr_z_yyyy_xxxy,   \
                             tr_z_yyyy_xxxyy,  \
                             tr_z_yyyy_xxxyz,  \
                             tr_z_yyyy_xxxz,   \
                             tr_z_yyyy_xxxzz,  \
                             tr_z_yyyy_xxyy,   \
                             tr_z_yyyy_xxyyy,  \
                             tr_z_yyyy_xxyyz,  \
                             tr_z_yyyy_xxyz,   \
                             tr_z_yyyy_xxyzz,  \
                             tr_z_yyyy_xxzz,   \
                             tr_z_yyyy_xxzzz,  \
                             tr_z_yyyy_xyyy,   \
                             tr_z_yyyy_xyyyy,  \
                             tr_z_yyyy_xyyyz,  \
                             tr_z_yyyy_xyyz,   \
                             tr_z_yyyy_xyyzz,  \
                             tr_z_yyyy_xyzz,   \
                             tr_z_yyyy_xyzzz,  \
                             tr_z_yyyy_xzzz,   \
                             tr_z_yyyy_xzzzz,  \
                             tr_z_yyyy_yyyy,   \
                             tr_z_yyyy_yyyyy,  \
                             tr_z_yyyy_yyyyz,  \
                             tr_z_yyyy_yyyz,   \
                             tr_z_yyyy_yyyzz,  \
                             tr_z_yyyy_yyzz,   \
                             tr_z_yyyy_yyzzz,  \
                             tr_z_yyyy_yzzz,   \
                             tr_z_yyyy_yzzzz,  \
                             tr_z_yyyy_zzzz,   \
                             tr_z_yyyy_zzzzz,  \
                             tr_z_yyyyy_xxxxx, \
                             tr_z_yyyyy_xxxxy, \
                             tr_z_yyyyy_xxxxz, \
                             tr_z_yyyyy_xxxyy, \
                             tr_z_yyyyy_xxxyz, \
                             tr_z_yyyyy_xxxzz, \
                             tr_z_yyyyy_xxyyy, \
                             tr_z_yyyyy_xxyyz, \
                             tr_z_yyyyy_xxyzz, \
                             tr_z_yyyyy_xxzzz, \
                             tr_z_yyyyy_xyyyy, \
                             tr_z_yyyyy_xyyyz, \
                             tr_z_yyyyy_xyyzz, \
                             tr_z_yyyyy_xyzzz, \
                             tr_z_yyyyy_xzzzz, \
                             tr_z_yyyyy_yyyyy, \
                             tr_z_yyyyy_yyyyz, \
                             tr_z_yyyyy_yyyzz, \
                             tr_z_yyyyy_yyzzz, \
                             tr_z_yyyyy_yzzzz, \
                             tr_z_yyyyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyy_xxxxx[i] = 4.0 * tr_z_yyy_xxxxx[i] * fe_0 + tr_z_yyyy_xxxxx[i] * pa_y[i];

        tr_z_yyyyy_xxxxy[i] = 4.0 * tr_z_yyy_xxxxy[i] * fe_0 + tr_z_yyyy_xxxx[i] * fe_0 + tr_z_yyyy_xxxxy[i] * pa_y[i];

        tr_z_yyyyy_xxxxz[i] = 4.0 * tr_z_yyy_xxxxz[i] * fe_0 + tr_z_yyyy_xxxxz[i] * pa_y[i];

        tr_z_yyyyy_xxxyy[i] = 4.0 * tr_z_yyy_xxxyy[i] * fe_0 + 2.0 * tr_z_yyyy_xxxy[i] * fe_0 + tr_z_yyyy_xxxyy[i] * pa_y[i];

        tr_z_yyyyy_xxxyz[i] = 4.0 * tr_z_yyy_xxxyz[i] * fe_0 + tr_z_yyyy_xxxz[i] * fe_0 + tr_z_yyyy_xxxyz[i] * pa_y[i];

        tr_z_yyyyy_xxxzz[i] = 4.0 * tr_z_yyy_xxxzz[i] * fe_0 + tr_z_yyyy_xxxzz[i] * pa_y[i];

        tr_z_yyyyy_xxyyy[i] = 4.0 * tr_z_yyy_xxyyy[i] * fe_0 + 3.0 * tr_z_yyyy_xxyy[i] * fe_0 + tr_z_yyyy_xxyyy[i] * pa_y[i];

        tr_z_yyyyy_xxyyz[i] = 4.0 * tr_z_yyy_xxyyz[i] * fe_0 + 2.0 * tr_z_yyyy_xxyz[i] * fe_0 + tr_z_yyyy_xxyyz[i] * pa_y[i];

        tr_z_yyyyy_xxyzz[i] = 4.0 * tr_z_yyy_xxyzz[i] * fe_0 + tr_z_yyyy_xxzz[i] * fe_0 + tr_z_yyyy_xxyzz[i] * pa_y[i];

        tr_z_yyyyy_xxzzz[i] = 4.0 * tr_z_yyy_xxzzz[i] * fe_0 + tr_z_yyyy_xxzzz[i] * pa_y[i];

        tr_z_yyyyy_xyyyy[i] = 4.0 * tr_z_yyy_xyyyy[i] * fe_0 + 4.0 * tr_z_yyyy_xyyy[i] * fe_0 + tr_z_yyyy_xyyyy[i] * pa_y[i];

        tr_z_yyyyy_xyyyz[i] = 4.0 * tr_z_yyy_xyyyz[i] * fe_0 + 3.0 * tr_z_yyyy_xyyz[i] * fe_0 + tr_z_yyyy_xyyyz[i] * pa_y[i];

        tr_z_yyyyy_xyyzz[i] = 4.0 * tr_z_yyy_xyyzz[i] * fe_0 + 2.0 * tr_z_yyyy_xyzz[i] * fe_0 + tr_z_yyyy_xyyzz[i] * pa_y[i];

        tr_z_yyyyy_xyzzz[i] = 4.0 * tr_z_yyy_xyzzz[i] * fe_0 + tr_z_yyyy_xzzz[i] * fe_0 + tr_z_yyyy_xyzzz[i] * pa_y[i];

        tr_z_yyyyy_xzzzz[i] = 4.0 * tr_z_yyy_xzzzz[i] * fe_0 + tr_z_yyyy_xzzzz[i] * pa_y[i];

        tr_z_yyyyy_yyyyy[i] = 4.0 * tr_z_yyy_yyyyy[i] * fe_0 + 5.0 * tr_z_yyyy_yyyy[i] * fe_0 + tr_z_yyyy_yyyyy[i] * pa_y[i];

        tr_z_yyyyy_yyyyz[i] = 4.0 * tr_z_yyy_yyyyz[i] * fe_0 + 4.0 * tr_z_yyyy_yyyz[i] * fe_0 + tr_z_yyyy_yyyyz[i] * pa_y[i];

        tr_z_yyyyy_yyyzz[i] = 4.0 * tr_z_yyy_yyyzz[i] * fe_0 + 3.0 * tr_z_yyyy_yyzz[i] * fe_0 + tr_z_yyyy_yyyzz[i] * pa_y[i];

        tr_z_yyyyy_yyzzz[i] = 4.0 * tr_z_yyy_yyzzz[i] * fe_0 + 2.0 * tr_z_yyyy_yzzz[i] * fe_0 + tr_z_yyyy_yyzzz[i] * pa_y[i];

        tr_z_yyyyy_yzzzz[i] = 4.0 * tr_z_yyy_yzzzz[i] * fe_0 + tr_z_yyyy_zzzz[i] * fe_0 + tr_z_yyyy_yzzzz[i] * pa_y[i];

        tr_z_yyyyy_zzzzz[i] = 4.0 * tr_z_yyy_zzzzz[i] * fe_0 + tr_z_yyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 1218-1239 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             tr_z_yyyy_xxxxy,  \
                             tr_z_yyyy_xxxyy,  \
                             tr_z_yyyy_xxyyy,  \
                             tr_z_yyyy_xyyyy,  \
                             tr_z_yyyy_yyyyy,  \
                             tr_z_yyyyz_xxxxx, \
                             tr_z_yyyyz_xxxxy, \
                             tr_z_yyyyz_xxxxz, \
                             tr_z_yyyyz_xxxyy, \
                             tr_z_yyyyz_xxxyz, \
                             tr_z_yyyyz_xxxzz, \
                             tr_z_yyyyz_xxyyy, \
                             tr_z_yyyyz_xxyyz, \
                             tr_z_yyyyz_xxyzz, \
                             tr_z_yyyyz_xxzzz, \
                             tr_z_yyyyz_xyyyy, \
                             tr_z_yyyyz_xyyyz, \
                             tr_z_yyyyz_xyyzz, \
                             tr_z_yyyyz_xyzzz, \
                             tr_z_yyyyz_xzzzz, \
                             tr_z_yyyyz_yyyyy, \
                             tr_z_yyyyz_yyyyz, \
                             tr_z_yyyyz_yyyzz, \
                             tr_z_yyyyz_yyzzz, \
                             tr_z_yyyyz_yzzzz, \
                             tr_z_yyyyz_zzzzz, \
                             tr_z_yyyz_xxxxx,  \
                             tr_z_yyyz_xxxxz,  \
                             tr_z_yyyz_xxxyz,  \
                             tr_z_yyyz_xxxz,   \
                             tr_z_yyyz_xxxzz,  \
                             tr_z_yyyz_xxyyz,  \
                             tr_z_yyyz_xxyz,   \
                             tr_z_yyyz_xxyzz,  \
                             tr_z_yyyz_xxzz,   \
                             tr_z_yyyz_xxzzz,  \
                             tr_z_yyyz_xyyyz,  \
                             tr_z_yyyz_xyyz,   \
                             tr_z_yyyz_xyyzz,  \
                             tr_z_yyyz_xyzz,   \
                             tr_z_yyyz_xyzzz,  \
                             tr_z_yyyz_xzzz,   \
                             tr_z_yyyz_xzzzz,  \
                             tr_z_yyyz_yyyyz,  \
                             tr_z_yyyz_yyyz,   \
                             tr_z_yyyz_yyyzz,  \
                             tr_z_yyyz_yyzz,   \
                             tr_z_yyyz_yyzzz,  \
                             tr_z_yyyz_yzzz,   \
                             tr_z_yyyz_yzzzz,  \
                             tr_z_yyyz_zzzz,   \
                             tr_z_yyyz_zzzzz,  \
                             tr_z_yyz_xxxxx,   \
                             tr_z_yyz_xxxxz,   \
                             tr_z_yyz_xxxyz,   \
                             tr_z_yyz_xxxzz,   \
                             tr_z_yyz_xxyyz,   \
                             tr_z_yyz_xxyzz,   \
                             tr_z_yyz_xxzzz,   \
                             tr_z_yyz_xyyyz,   \
                             tr_z_yyz_xyyzz,   \
                             tr_z_yyz_xyzzz,   \
                             tr_z_yyz_xzzzz,   \
                             tr_z_yyz_yyyyz,   \
                             tr_z_yyz_yyyzz,   \
                             tr_z_yyz_yyzzz,   \
                             tr_z_yyz_yzzzz,   \
                             tr_z_yyz_zzzzz,   \
                             ts_yyyy_xxxxy,    \
                             ts_yyyy_xxxyy,    \
                             ts_yyyy_xxyyy,    \
                             ts_yyyy_xyyyy,    \
                             ts_yyyy_yyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyz_xxxxx[i] = 3.0 * tr_z_yyz_xxxxx[i] * fe_0 + tr_z_yyyz_xxxxx[i] * pa_y[i];

        tr_z_yyyyz_xxxxy[i] = ts_yyyy_xxxxy[i] * fe_0 + tr_z_yyyy_xxxxy[i] * pa_z[i];

        tr_z_yyyyz_xxxxz[i] = 3.0 * tr_z_yyz_xxxxz[i] * fe_0 + tr_z_yyyz_xxxxz[i] * pa_y[i];

        tr_z_yyyyz_xxxyy[i] = ts_yyyy_xxxyy[i] * fe_0 + tr_z_yyyy_xxxyy[i] * pa_z[i];

        tr_z_yyyyz_xxxyz[i] = 3.0 * tr_z_yyz_xxxyz[i] * fe_0 + tr_z_yyyz_xxxz[i] * fe_0 + tr_z_yyyz_xxxyz[i] * pa_y[i];

        tr_z_yyyyz_xxxzz[i] = 3.0 * tr_z_yyz_xxxzz[i] * fe_0 + tr_z_yyyz_xxxzz[i] * pa_y[i];

        tr_z_yyyyz_xxyyy[i] = ts_yyyy_xxyyy[i] * fe_0 + tr_z_yyyy_xxyyy[i] * pa_z[i];

        tr_z_yyyyz_xxyyz[i] = 3.0 * tr_z_yyz_xxyyz[i] * fe_0 + 2.0 * tr_z_yyyz_xxyz[i] * fe_0 + tr_z_yyyz_xxyyz[i] * pa_y[i];

        tr_z_yyyyz_xxyzz[i] = 3.0 * tr_z_yyz_xxyzz[i] * fe_0 + tr_z_yyyz_xxzz[i] * fe_0 + tr_z_yyyz_xxyzz[i] * pa_y[i];

        tr_z_yyyyz_xxzzz[i] = 3.0 * tr_z_yyz_xxzzz[i] * fe_0 + tr_z_yyyz_xxzzz[i] * pa_y[i];

        tr_z_yyyyz_xyyyy[i] = ts_yyyy_xyyyy[i] * fe_0 + tr_z_yyyy_xyyyy[i] * pa_z[i];

        tr_z_yyyyz_xyyyz[i] = 3.0 * tr_z_yyz_xyyyz[i] * fe_0 + 3.0 * tr_z_yyyz_xyyz[i] * fe_0 + tr_z_yyyz_xyyyz[i] * pa_y[i];

        tr_z_yyyyz_xyyzz[i] = 3.0 * tr_z_yyz_xyyzz[i] * fe_0 + 2.0 * tr_z_yyyz_xyzz[i] * fe_0 + tr_z_yyyz_xyyzz[i] * pa_y[i];

        tr_z_yyyyz_xyzzz[i] = 3.0 * tr_z_yyz_xyzzz[i] * fe_0 + tr_z_yyyz_xzzz[i] * fe_0 + tr_z_yyyz_xyzzz[i] * pa_y[i];

        tr_z_yyyyz_xzzzz[i] = 3.0 * tr_z_yyz_xzzzz[i] * fe_0 + tr_z_yyyz_xzzzz[i] * pa_y[i];

        tr_z_yyyyz_yyyyy[i] = ts_yyyy_yyyyy[i] * fe_0 + tr_z_yyyy_yyyyy[i] * pa_z[i];

        tr_z_yyyyz_yyyyz[i] = 3.0 * tr_z_yyz_yyyyz[i] * fe_0 + 4.0 * tr_z_yyyz_yyyz[i] * fe_0 + tr_z_yyyz_yyyyz[i] * pa_y[i];

        tr_z_yyyyz_yyyzz[i] = 3.0 * tr_z_yyz_yyyzz[i] * fe_0 + 3.0 * tr_z_yyyz_yyzz[i] * fe_0 + tr_z_yyyz_yyyzz[i] * pa_y[i];

        tr_z_yyyyz_yyzzz[i] = 3.0 * tr_z_yyz_yyzzz[i] * fe_0 + 2.0 * tr_z_yyyz_yzzz[i] * fe_0 + tr_z_yyyz_yyzzz[i] * pa_y[i];

        tr_z_yyyyz_yzzzz[i] = 3.0 * tr_z_yyz_yzzzz[i] * fe_0 + tr_z_yyyz_zzzz[i] * fe_0 + tr_z_yyyz_yzzzz[i] * pa_y[i];

        tr_z_yyyyz_zzzzz[i] = 3.0 * tr_z_yyz_zzzzz[i] * fe_0 + tr_z_yyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 1239-1260 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyyzz_xxxxx, \
                             tr_z_yyyzz_xxxxy, \
                             tr_z_yyyzz_xxxxz, \
                             tr_z_yyyzz_xxxyy, \
                             tr_z_yyyzz_xxxyz, \
                             tr_z_yyyzz_xxxzz, \
                             tr_z_yyyzz_xxyyy, \
                             tr_z_yyyzz_xxyyz, \
                             tr_z_yyyzz_xxyzz, \
                             tr_z_yyyzz_xxzzz, \
                             tr_z_yyyzz_xyyyy, \
                             tr_z_yyyzz_xyyyz, \
                             tr_z_yyyzz_xyyzz, \
                             tr_z_yyyzz_xyzzz, \
                             tr_z_yyyzz_xzzzz, \
                             tr_z_yyyzz_yyyyy, \
                             tr_z_yyyzz_yyyyz, \
                             tr_z_yyyzz_yyyzz, \
                             tr_z_yyyzz_yyzzz, \
                             tr_z_yyyzz_yzzzz, \
                             tr_z_yyyzz_zzzzz, \
                             tr_z_yyzz_xxxx,   \
                             tr_z_yyzz_xxxxx,  \
                             tr_z_yyzz_xxxxy,  \
                             tr_z_yyzz_xxxxz,  \
                             tr_z_yyzz_xxxy,   \
                             tr_z_yyzz_xxxyy,  \
                             tr_z_yyzz_xxxyz,  \
                             tr_z_yyzz_xxxz,   \
                             tr_z_yyzz_xxxzz,  \
                             tr_z_yyzz_xxyy,   \
                             tr_z_yyzz_xxyyy,  \
                             tr_z_yyzz_xxyyz,  \
                             tr_z_yyzz_xxyz,   \
                             tr_z_yyzz_xxyzz,  \
                             tr_z_yyzz_xxzz,   \
                             tr_z_yyzz_xxzzz,  \
                             tr_z_yyzz_xyyy,   \
                             tr_z_yyzz_xyyyy,  \
                             tr_z_yyzz_xyyyz,  \
                             tr_z_yyzz_xyyz,   \
                             tr_z_yyzz_xyyzz,  \
                             tr_z_yyzz_xyzz,   \
                             tr_z_yyzz_xyzzz,  \
                             tr_z_yyzz_xzzz,   \
                             tr_z_yyzz_xzzzz,  \
                             tr_z_yyzz_yyyy,   \
                             tr_z_yyzz_yyyyy,  \
                             tr_z_yyzz_yyyyz,  \
                             tr_z_yyzz_yyyz,   \
                             tr_z_yyzz_yyyzz,  \
                             tr_z_yyzz_yyzz,   \
                             tr_z_yyzz_yyzzz,  \
                             tr_z_yyzz_yzzz,   \
                             tr_z_yyzz_yzzzz,  \
                             tr_z_yyzz_zzzz,   \
                             tr_z_yyzz_zzzzz,  \
                             tr_z_yzz_xxxxx,   \
                             tr_z_yzz_xxxxy,   \
                             tr_z_yzz_xxxxz,   \
                             tr_z_yzz_xxxyy,   \
                             tr_z_yzz_xxxyz,   \
                             tr_z_yzz_xxxzz,   \
                             tr_z_yzz_xxyyy,   \
                             tr_z_yzz_xxyyz,   \
                             tr_z_yzz_xxyzz,   \
                             tr_z_yzz_xxzzz,   \
                             tr_z_yzz_xyyyy,   \
                             tr_z_yzz_xyyyz,   \
                             tr_z_yzz_xyyzz,   \
                             tr_z_yzz_xyzzz,   \
                             tr_z_yzz_xzzzz,   \
                             tr_z_yzz_yyyyy,   \
                             tr_z_yzz_yyyyz,   \
                             tr_z_yzz_yyyzz,   \
                             tr_z_yzz_yyzzz,   \
                             tr_z_yzz_yzzzz,   \
                             tr_z_yzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzz_xxxxx[i] = 2.0 * tr_z_yzz_xxxxx[i] * fe_0 + tr_z_yyzz_xxxxx[i] * pa_y[i];

        tr_z_yyyzz_xxxxy[i] = 2.0 * tr_z_yzz_xxxxy[i] * fe_0 + tr_z_yyzz_xxxx[i] * fe_0 + tr_z_yyzz_xxxxy[i] * pa_y[i];

        tr_z_yyyzz_xxxxz[i] = 2.0 * tr_z_yzz_xxxxz[i] * fe_0 + tr_z_yyzz_xxxxz[i] * pa_y[i];

        tr_z_yyyzz_xxxyy[i] = 2.0 * tr_z_yzz_xxxyy[i] * fe_0 + 2.0 * tr_z_yyzz_xxxy[i] * fe_0 + tr_z_yyzz_xxxyy[i] * pa_y[i];

        tr_z_yyyzz_xxxyz[i] = 2.0 * tr_z_yzz_xxxyz[i] * fe_0 + tr_z_yyzz_xxxz[i] * fe_0 + tr_z_yyzz_xxxyz[i] * pa_y[i];

        tr_z_yyyzz_xxxzz[i] = 2.0 * tr_z_yzz_xxxzz[i] * fe_0 + tr_z_yyzz_xxxzz[i] * pa_y[i];

        tr_z_yyyzz_xxyyy[i] = 2.0 * tr_z_yzz_xxyyy[i] * fe_0 + 3.0 * tr_z_yyzz_xxyy[i] * fe_0 + tr_z_yyzz_xxyyy[i] * pa_y[i];

        tr_z_yyyzz_xxyyz[i] = 2.0 * tr_z_yzz_xxyyz[i] * fe_0 + 2.0 * tr_z_yyzz_xxyz[i] * fe_0 + tr_z_yyzz_xxyyz[i] * pa_y[i];

        tr_z_yyyzz_xxyzz[i] = 2.0 * tr_z_yzz_xxyzz[i] * fe_0 + tr_z_yyzz_xxzz[i] * fe_0 + tr_z_yyzz_xxyzz[i] * pa_y[i];

        tr_z_yyyzz_xxzzz[i] = 2.0 * tr_z_yzz_xxzzz[i] * fe_0 + tr_z_yyzz_xxzzz[i] * pa_y[i];

        tr_z_yyyzz_xyyyy[i] = 2.0 * tr_z_yzz_xyyyy[i] * fe_0 + 4.0 * tr_z_yyzz_xyyy[i] * fe_0 + tr_z_yyzz_xyyyy[i] * pa_y[i];

        tr_z_yyyzz_xyyyz[i] = 2.0 * tr_z_yzz_xyyyz[i] * fe_0 + 3.0 * tr_z_yyzz_xyyz[i] * fe_0 + tr_z_yyzz_xyyyz[i] * pa_y[i];

        tr_z_yyyzz_xyyzz[i] = 2.0 * tr_z_yzz_xyyzz[i] * fe_0 + 2.0 * tr_z_yyzz_xyzz[i] * fe_0 + tr_z_yyzz_xyyzz[i] * pa_y[i];

        tr_z_yyyzz_xyzzz[i] = 2.0 * tr_z_yzz_xyzzz[i] * fe_0 + tr_z_yyzz_xzzz[i] * fe_0 + tr_z_yyzz_xyzzz[i] * pa_y[i];

        tr_z_yyyzz_xzzzz[i] = 2.0 * tr_z_yzz_xzzzz[i] * fe_0 + tr_z_yyzz_xzzzz[i] * pa_y[i];

        tr_z_yyyzz_yyyyy[i] = 2.0 * tr_z_yzz_yyyyy[i] * fe_0 + 5.0 * tr_z_yyzz_yyyy[i] * fe_0 + tr_z_yyzz_yyyyy[i] * pa_y[i];

        tr_z_yyyzz_yyyyz[i] = 2.0 * tr_z_yzz_yyyyz[i] * fe_0 + 4.0 * tr_z_yyzz_yyyz[i] * fe_0 + tr_z_yyzz_yyyyz[i] * pa_y[i];

        tr_z_yyyzz_yyyzz[i] = 2.0 * tr_z_yzz_yyyzz[i] * fe_0 + 3.0 * tr_z_yyzz_yyzz[i] * fe_0 + tr_z_yyzz_yyyzz[i] * pa_y[i];

        tr_z_yyyzz_yyzzz[i] = 2.0 * tr_z_yzz_yyzzz[i] * fe_0 + 2.0 * tr_z_yyzz_yzzz[i] * fe_0 + tr_z_yyzz_yyzzz[i] * pa_y[i];

        tr_z_yyyzz_yzzzz[i] = 2.0 * tr_z_yzz_yzzzz[i] * fe_0 + tr_z_yyzz_zzzz[i] * fe_0 + tr_z_yyzz_yzzzz[i] * pa_y[i];

        tr_z_yyyzz_zzzzz[i] = 2.0 * tr_z_yzz_zzzzz[i] * fe_0 + tr_z_yyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1260-1281 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yyzzz_xxxxx, \
                             tr_z_yyzzz_xxxxy, \
                             tr_z_yyzzz_xxxxz, \
                             tr_z_yyzzz_xxxyy, \
                             tr_z_yyzzz_xxxyz, \
                             tr_z_yyzzz_xxxzz, \
                             tr_z_yyzzz_xxyyy, \
                             tr_z_yyzzz_xxyyz, \
                             tr_z_yyzzz_xxyzz, \
                             tr_z_yyzzz_xxzzz, \
                             tr_z_yyzzz_xyyyy, \
                             tr_z_yyzzz_xyyyz, \
                             tr_z_yyzzz_xyyzz, \
                             tr_z_yyzzz_xyzzz, \
                             tr_z_yyzzz_xzzzz, \
                             tr_z_yyzzz_yyyyy, \
                             tr_z_yyzzz_yyyyz, \
                             tr_z_yyzzz_yyyzz, \
                             tr_z_yyzzz_yyzzz, \
                             tr_z_yyzzz_yzzzz, \
                             tr_z_yyzzz_zzzzz, \
                             tr_z_yzzz_xxxx,   \
                             tr_z_yzzz_xxxxx,  \
                             tr_z_yzzz_xxxxy,  \
                             tr_z_yzzz_xxxxz,  \
                             tr_z_yzzz_xxxy,   \
                             tr_z_yzzz_xxxyy,  \
                             tr_z_yzzz_xxxyz,  \
                             tr_z_yzzz_xxxz,   \
                             tr_z_yzzz_xxxzz,  \
                             tr_z_yzzz_xxyy,   \
                             tr_z_yzzz_xxyyy,  \
                             tr_z_yzzz_xxyyz,  \
                             tr_z_yzzz_xxyz,   \
                             tr_z_yzzz_xxyzz,  \
                             tr_z_yzzz_xxzz,   \
                             tr_z_yzzz_xxzzz,  \
                             tr_z_yzzz_xyyy,   \
                             tr_z_yzzz_xyyyy,  \
                             tr_z_yzzz_xyyyz,  \
                             tr_z_yzzz_xyyz,   \
                             tr_z_yzzz_xyyzz,  \
                             tr_z_yzzz_xyzz,   \
                             tr_z_yzzz_xyzzz,  \
                             tr_z_yzzz_xzzz,   \
                             tr_z_yzzz_xzzzz,  \
                             tr_z_yzzz_yyyy,   \
                             tr_z_yzzz_yyyyy,  \
                             tr_z_yzzz_yyyyz,  \
                             tr_z_yzzz_yyyz,   \
                             tr_z_yzzz_yyyzz,  \
                             tr_z_yzzz_yyzz,   \
                             tr_z_yzzz_yyzzz,  \
                             tr_z_yzzz_yzzz,   \
                             tr_z_yzzz_yzzzz,  \
                             tr_z_yzzz_zzzz,   \
                             tr_z_yzzz_zzzzz,  \
                             tr_z_zzz_xxxxx,   \
                             tr_z_zzz_xxxxy,   \
                             tr_z_zzz_xxxxz,   \
                             tr_z_zzz_xxxyy,   \
                             tr_z_zzz_xxxyz,   \
                             tr_z_zzz_xxxzz,   \
                             tr_z_zzz_xxyyy,   \
                             tr_z_zzz_xxyyz,   \
                             tr_z_zzz_xxyzz,   \
                             tr_z_zzz_xxzzz,   \
                             tr_z_zzz_xyyyy,   \
                             tr_z_zzz_xyyyz,   \
                             tr_z_zzz_xyyzz,   \
                             tr_z_zzz_xyzzz,   \
                             tr_z_zzz_xzzzz,   \
                             tr_z_zzz_yyyyy,   \
                             tr_z_zzz_yyyyz,   \
                             tr_z_zzz_yyyzz,   \
                             tr_z_zzz_yyzzz,   \
                             tr_z_zzz_yzzzz,   \
                             tr_z_zzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzz_xxxxx[i] = tr_z_zzz_xxxxx[i] * fe_0 + tr_z_yzzz_xxxxx[i] * pa_y[i];

        tr_z_yyzzz_xxxxy[i] = tr_z_zzz_xxxxy[i] * fe_0 + tr_z_yzzz_xxxx[i] * fe_0 + tr_z_yzzz_xxxxy[i] * pa_y[i];

        tr_z_yyzzz_xxxxz[i] = tr_z_zzz_xxxxz[i] * fe_0 + tr_z_yzzz_xxxxz[i] * pa_y[i];

        tr_z_yyzzz_xxxyy[i] = tr_z_zzz_xxxyy[i] * fe_0 + 2.0 * tr_z_yzzz_xxxy[i] * fe_0 + tr_z_yzzz_xxxyy[i] * pa_y[i];

        tr_z_yyzzz_xxxyz[i] = tr_z_zzz_xxxyz[i] * fe_0 + tr_z_yzzz_xxxz[i] * fe_0 + tr_z_yzzz_xxxyz[i] * pa_y[i];

        tr_z_yyzzz_xxxzz[i] = tr_z_zzz_xxxzz[i] * fe_0 + tr_z_yzzz_xxxzz[i] * pa_y[i];

        tr_z_yyzzz_xxyyy[i] = tr_z_zzz_xxyyy[i] * fe_0 + 3.0 * tr_z_yzzz_xxyy[i] * fe_0 + tr_z_yzzz_xxyyy[i] * pa_y[i];

        tr_z_yyzzz_xxyyz[i] = tr_z_zzz_xxyyz[i] * fe_0 + 2.0 * tr_z_yzzz_xxyz[i] * fe_0 + tr_z_yzzz_xxyyz[i] * pa_y[i];

        tr_z_yyzzz_xxyzz[i] = tr_z_zzz_xxyzz[i] * fe_0 + tr_z_yzzz_xxzz[i] * fe_0 + tr_z_yzzz_xxyzz[i] * pa_y[i];

        tr_z_yyzzz_xxzzz[i] = tr_z_zzz_xxzzz[i] * fe_0 + tr_z_yzzz_xxzzz[i] * pa_y[i];

        tr_z_yyzzz_xyyyy[i] = tr_z_zzz_xyyyy[i] * fe_0 + 4.0 * tr_z_yzzz_xyyy[i] * fe_0 + tr_z_yzzz_xyyyy[i] * pa_y[i];

        tr_z_yyzzz_xyyyz[i] = tr_z_zzz_xyyyz[i] * fe_0 + 3.0 * tr_z_yzzz_xyyz[i] * fe_0 + tr_z_yzzz_xyyyz[i] * pa_y[i];

        tr_z_yyzzz_xyyzz[i] = tr_z_zzz_xyyzz[i] * fe_0 + 2.0 * tr_z_yzzz_xyzz[i] * fe_0 + tr_z_yzzz_xyyzz[i] * pa_y[i];

        tr_z_yyzzz_xyzzz[i] = tr_z_zzz_xyzzz[i] * fe_0 + tr_z_yzzz_xzzz[i] * fe_0 + tr_z_yzzz_xyzzz[i] * pa_y[i];

        tr_z_yyzzz_xzzzz[i] = tr_z_zzz_xzzzz[i] * fe_0 + tr_z_yzzz_xzzzz[i] * pa_y[i];

        tr_z_yyzzz_yyyyy[i] = tr_z_zzz_yyyyy[i] * fe_0 + 5.0 * tr_z_yzzz_yyyy[i] * fe_0 + tr_z_yzzz_yyyyy[i] * pa_y[i];

        tr_z_yyzzz_yyyyz[i] = tr_z_zzz_yyyyz[i] * fe_0 + 4.0 * tr_z_yzzz_yyyz[i] * fe_0 + tr_z_yzzz_yyyyz[i] * pa_y[i];

        tr_z_yyzzz_yyyzz[i] = tr_z_zzz_yyyzz[i] * fe_0 + 3.0 * tr_z_yzzz_yyzz[i] * fe_0 + tr_z_yzzz_yyyzz[i] * pa_y[i];

        tr_z_yyzzz_yyzzz[i] = tr_z_zzz_yyzzz[i] * fe_0 + 2.0 * tr_z_yzzz_yzzz[i] * fe_0 + tr_z_yzzz_yyzzz[i] * pa_y[i];

        tr_z_yyzzz_yzzzz[i] = tr_z_zzz_yzzzz[i] * fe_0 + tr_z_yzzz_zzzz[i] * fe_0 + tr_z_yzzz_yzzzz[i] * pa_y[i];

        tr_z_yyzzz_zzzzz[i] = tr_z_zzz_zzzzz[i] * fe_0 + tr_z_yzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1281-1302 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_y,                 \
                             tr_z_yzzzz_xxxxx, \
                             tr_z_yzzzz_xxxxy, \
                             tr_z_yzzzz_xxxxz, \
                             tr_z_yzzzz_xxxyy, \
                             tr_z_yzzzz_xxxyz, \
                             tr_z_yzzzz_xxxzz, \
                             tr_z_yzzzz_xxyyy, \
                             tr_z_yzzzz_xxyyz, \
                             tr_z_yzzzz_xxyzz, \
                             tr_z_yzzzz_xxzzz, \
                             tr_z_yzzzz_xyyyy, \
                             tr_z_yzzzz_xyyyz, \
                             tr_z_yzzzz_xyyzz, \
                             tr_z_yzzzz_xyzzz, \
                             tr_z_yzzzz_xzzzz, \
                             tr_z_yzzzz_yyyyy, \
                             tr_z_yzzzz_yyyyz, \
                             tr_z_yzzzz_yyyzz, \
                             tr_z_yzzzz_yyzzz, \
                             tr_z_yzzzz_yzzzz, \
                             tr_z_yzzzz_zzzzz, \
                             tr_z_zzzz_xxxx,   \
                             tr_z_zzzz_xxxxx,  \
                             tr_z_zzzz_xxxxy,  \
                             tr_z_zzzz_xxxxz,  \
                             tr_z_zzzz_xxxy,   \
                             tr_z_zzzz_xxxyy,  \
                             tr_z_zzzz_xxxyz,  \
                             tr_z_zzzz_xxxz,   \
                             tr_z_zzzz_xxxzz,  \
                             tr_z_zzzz_xxyy,   \
                             tr_z_zzzz_xxyyy,  \
                             tr_z_zzzz_xxyyz,  \
                             tr_z_zzzz_xxyz,   \
                             tr_z_zzzz_xxyzz,  \
                             tr_z_zzzz_xxzz,   \
                             tr_z_zzzz_xxzzz,  \
                             tr_z_zzzz_xyyy,   \
                             tr_z_zzzz_xyyyy,  \
                             tr_z_zzzz_xyyyz,  \
                             tr_z_zzzz_xyyz,   \
                             tr_z_zzzz_xyyzz,  \
                             tr_z_zzzz_xyzz,   \
                             tr_z_zzzz_xyzzz,  \
                             tr_z_zzzz_xzzz,   \
                             tr_z_zzzz_xzzzz,  \
                             tr_z_zzzz_yyyy,   \
                             tr_z_zzzz_yyyyy,  \
                             tr_z_zzzz_yyyyz,  \
                             tr_z_zzzz_yyyz,   \
                             tr_z_zzzz_yyyzz,  \
                             tr_z_zzzz_yyzz,   \
                             tr_z_zzzz_yyzzz,  \
                             tr_z_zzzz_yzzz,   \
                             tr_z_zzzz_yzzzz,  \
                             tr_z_zzzz_zzzz,   \
                             tr_z_zzzz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzz_xxxxx[i] = tr_z_zzzz_xxxxx[i] * pa_y[i];

        tr_z_yzzzz_xxxxy[i] = tr_z_zzzz_xxxx[i] * fe_0 + tr_z_zzzz_xxxxy[i] * pa_y[i];

        tr_z_yzzzz_xxxxz[i] = tr_z_zzzz_xxxxz[i] * pa_y[i];

        tr_z_yzzzz_xxxyy[i] = 2.0 * tr_z_zzzz_xxxy[i] * fe_0 + tr_z_zzzz_xxxyy[i] * pa_y[i];

        tr_z_yzzzz_xxxyz[i] = tr_z_zzzz_xxxz[i] * fe_0 + tr_z_zzzz_xxxyz[i] * pa_y[i];

        tr_z_yzzzz_xxxzz[i] = tr_z_zzzz_xxxzz[i] * pa_y[i];

        tr_z_yzzzz_xxyyy[i] = 3.0 * tr_z_zzzz_xxyy[i] * fe_0 + tr_z_zzzz_xxyyy[i] * pa_y[i];

        tr_z_yzzzz_xxyyz[i] = 2.0 * tr_z_zzzz_xxyz[i] * fe_0 + tr_z_zzzz_xxyyz[i] * pa_y[i];

        tr_z_yzzzz_xxyzz[i] = tr_z_zzzz_xxzz[i] * fe_0 + tr_z_zzzz_xxyzz[i] * pa_y[i];

        tr_z_yzzzz_xxzzz[i] = tr_z_zzzz_xxzzz[i] * pa_y[i];

        tr_z_yzzzz_xyyyy[i] = 4.0 * tr_z_zzzz_xyyy[i] * fe_0 + tr_z_zzzz_xyyyy[i] * pa_y[i];

        tr_z_yzzzz_xyyyz[i] = 3.0 * tr_z_zzzz_xyyz[i] * fe_0 + tr_z_zzzz_xyyyz[i] * pa_y[i];

        tr_z_yzzzz_xyyzz[i] = 2.0 * tr_z_zzzz_xyzz[i] * fe_0 + tr_z_zzzz_xyyzz[i] * pa_y[i];

        tr_z_yzzzz_xyzzz[i] = tr_z_zzzz_xzzz[i] * fe_0 + tr_z_zzzz_xyzzz[i] * pa_y[i];

        tr_z_yzzzz_xzzzz[i] = tr_z_zzzz_xzzzz[i] * pa_y[i];

        tr_z_yzzzz_yyyyy[i] = 5.0 * tr_z_zzzz_yyyy[i] * fe_0 + tr_z_zzzz_yyyyy[i] * pa_y[i];

        tr_z_yzzzz_yyyyz[i] = 4.0 * tr_z_zzzz_yyyz[i] * fe_0 + tr_z_zzzz_yyyyz[i] * pa_y[i];

        tr_z_yzzzz_yyyzz[i] = 3.0 * tr_z_zzzz_yyzz[i] * fe_0 + tr_z_zzzz_yyyzz[i] * pa_y[i];

        tr_z_yzzzz_yyzzz[i] = 2.0 * tr_z_zzzz_yzzz[i] * fe_0 + tr_z_zzzz_yyzzz[i] * pa_y[i];

        tr_z_yzzzz_yzzzz[i] = tr_z_zzzz_zzzz[i] * fe_0 + tr_z_zzzz_yzzzz[i] * pa_y[i];

        tr_z_yzzzz_zzzzz[i] = tr_z_zzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 1302-1323 components of targeted buffer : HH

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

#pragma omp simd aligned(pa_z,                 \
                             tr_z_zzz_xxxxx,   \
                             tr_z_zzz_xxxxy,   \
                             tr_z_zzz_xxxxz,   \
                             tr_z_zzz_xxxyy,   \
                             tr_z_zzz_xxxyz,   \
                             tr_z_zzz_xxxzz,   \
                             tr_z_zzz_xxyyy,   \
                             tr_z_zzz_xxyyz,   \
                             tr_z_zzz_xxyzz,   \
                             tr_z_zzz_xxzzz,   \
                             tr_z_zzz_xyyyy,   \
                             tr_z_zzz_xyyyz,   \
                             tr_z_zzz_xyyzz,   \
                             tr_z_zzz_xyzzz,   \
                             tr_z_zzz_xzzzz,   \
                             tr_z_zzz_yyyyy,   \
                             tr_z_zzz_yyyyz,   \
                             tr_z_zzz_yyyzz,   \
                             tr_z_zzz_yyzzz,   \
                             tr_z_zzz_yzzzz,   \
                             tr_z_zzz_zzzzz,   \
                             tr_z_zzzz_xxxx,   \
                             tr_z_zzzz_xxxxx,  \
                             tr_z_zzzz_xxxxy,  \
                             tr_z_zzzz_xxxxz,  \
                             tr_z_zzzz_xxxy,   \
                             tr_z_zzzz_xxxyy,  \
                             tr_z_zzzz_xxxyz,  \
                             tr_z_zzzz_xxxz,   \
                             tr_z_zzzz_xxxzz,  \
                             tr_z_zzzz_xxyy,   \
                             tr_z_zzzz_xxyyy,  \
                             tr_z_zzzz_xxyyz,  \
                             tr_z_zzzz_xxyz,   \
                             tr_z_zzzz_xxyzz,  \
                             tr_z_zzzz_xxzz,   \
                             tr_z_zzzz_xxzzz,  \
                             tr_z_zzzz_xyyy,   \
                             tr_z_zzzz_xyyyy,  \
                             tr_z_zzzz_xyyyz,  \
                             tr_z_zzzz_xyyz,   \
                             tr_z_zzzz_xyyzz,  \
                             tr_z_zzzz_xyzz,   \
                             tr_z_zzzz_xyzzz,  \
                             tr_z_zzzz_xzzz,   \
                             tr_z_zzzz_xzzzz,  \
                             tr_z_zzzz_yyyy,   \
                             tr_z_zzzz_yyyyy,  \
                             tr_z_zzzz_yyyyz,  \
                             tr_z_zzzz_yyyz,   \
                             tr_z_zzzz_yyyzz,  \
                             tr_z_zzzz_yyzz,   \
                             tr_z_zzzz_yyzzz,  \
                             tr_z_zzzz_yzzz,   \
                             tr_z_zzzz_yzzzz,  \
                             tr_z_zzzz_zzzz,   \
                             tr_z_zzzz_zzzzz,  \
                             tr_z_zzzzz_xxxxx, \
                             tr_z_zzzzz_xxxxy, \
                             tr_z_zzzzz_xxxxz, \
                             tr_z_zzzzz_xxxyy, \
                             tr_z_zzzzz_xxxyz, \
                             tr_z_zzzzz_xxxzz, \
                             tr_z_zzzzz_xxyyy, \
                             tr_z_zzzzz_xxyyz, \
                             tr_z_zzzzz_xxyzz, \
                             tr_z_zzzzz_xxzzz, \
                             tr_z_zzzzz_xyyyy, \
                             tr_z_zzzzz_xyyyz, \
                             tr_z_zzzzz_xyyzz, \
                             tr_z_zzzzz_xyzzz, \
                             tr_z_zzzzz_xzzzz, \
                             tr_z_zzzzz_yyyyy, \
                             tr_z_zzzzz_yyyyz, \
                             tr_z_zzzzz_yyyzz, \
                             tr_z_zzzzz_yyzzz, \
                             tr_z_zzzzz_yzzzz, \
                             tr_z_zzzzz_zzzzz, \
                             ts_zzzz_xxxxx,    \
                             ts_zzzz_xxxxy,    \
                             ts_zzzz_xxxxz,    \
                             ts_zzzz_xxxyy,    \
                             ts_zzzz_xxxyz,    \
                             ts_zzzz_xxxzz,    \
                             ts_zzzz_xxyyy,    \
                             ts_zzzz_xxyyz,    \
                             ts_zzzz_xxyzz,    \
                             ts_zzzz_xxzzz,    \
                             ts_zzzz_xyyyy,    \
                             ts_zzzz_xyyyz,    \
                             ts_zzzz_xyyzz,    \
                             ts_zzzz_xyzzz,    \
                             ts_zzzz_xzzzz,    \
                             ts_zzzz_yyyyy,    \
                             ts_zzzz_yyyyz,    \
                             ts_zzzz_yyyzz,    \
                             ts_zzzz_yyzzz,    \
                             ts_zzzz_yzzzz,    \
                             ts_zzzz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzz_xxxxx[i] = 4.0 * tr_z_zzz_xxxxx[i] * fe_0 + ts_zzzz_xxxxx[i] * fe_0 + tr_z_zzzz_xxxxx[i] * pa_z[i];

        tr_z_zzzzz_xxxxy[i] = 4.0 * tr_z_zzz_xxxxy[i] * fe_0 + ts_zzzz_xxxxy[i] * fe_0 + tr_z_zzzz_xxxxy[i] * pa_z[i];

        tr_z_zzzzz_xxxxz[i] = 4.0 * tr_z_zzz_xxxxz[i] * fe_0 + tr_z_zzzz_xxxx[i] * fe_0 + ts_zzzz_xxxxz[i] * fe_0 + tr_z_zzzz_xxxxz[i] * pa_z[i];

        tr_z_zzzzz_xxxyy[i] = 4.0 * tr_z_zzz_xxxyy[i] * fe_0 + ts_zzzz_xxxyy[i] * fe_0 + tr_z_zzzz_xxxyy[i] * pa_z[i];

        tr_z_zzzzz_xxxyz[i] = 4.0 * tr_z_zzz_xxxyz[i] * fe_0 + tr_z_zzzz_xxxy[i] * fe_0 + ts_zzzz_xxxyz[i] * fe_0 + tr_z_zzzz_xxxyz[i] * pa_z[i];

        tr_z_zzzzz_xxxzz[i] =
            4.0 * tr_z_zzz_xxxzz[i] * fe_0 + 2.0 * tr_z_zzzz_xxxz[i] * fe_0 + ts_zzzz_xxxzz[i] * fe_0 + tr_z_zzzz_xxxzz[i] * pa_z[i];

        tr_z_zzzzz_xxyyy[i] = 4.0 * tr_z_zzz_xxyyy[i] * fe_0 + ts_zzzz_xxyyy[i] * fe_0 + tr_z_zzzz_xxyyy[i] * pa_z[i];

        tr_z_zzzzz_xxyyz[i] = 4.0 * tr_z_zzz_xxyyz[i] * fe_0 + tr_z_zzzz_xxyy[i] * fe_0 + ts_zzzz_xxyyz[i] * fe_0 + tr_z_zzzz_xxyyz[i] * pa_z[i];

        tr_z_zzzzz_xxyzz[i] =
            4.0 * tr_z_zzz_xxyzz[i] * fe_0 + 2.0 * tr_z_zzzz_xxyz[i] * fe_0 + ts_zzzz_xxyzz[i] * fe_0 + tr_z_zzzz_xxyzz[i] * pa_z[i];

        tr_z_zzzzz_xxzzz[i] =
            4.0 * tr_z_zzz_xxzzz[i] * fe_0 + 3.0 * tr_z_zzzz_xxzz[i] * fe_0 + ts_zzzz_xxzzz[i] * fe_0 + tr_z_zzzz_xxzzz[i] * pa_z[i];

        tr_z_zzzzz_xyyyy[i] = 4.0 * tr_z_zzz_xyyyy[i] * fe_0 + ts_zzzz_xyyyy[i] * fe_0 + tr_z_zzzz_xyyyy[i] * pa_z[i];

        tr_z_zzzzz_xyyyz[i] = 4.0 * tr_z_zzz_xyyyz[i] * fe_0 + tr_z_zzzz_xyyy[i] * fe_0 + ts_zzzz_xyyyz[i] * fe_0 + tr_z_zzzz_xyyyz[i] * pa_z[i];

        tr_z_zzzzz_xyyzz[i] =
            4.0 * tr_z_zzz_xyyzz[i] * fe_0 + 2.0 * tr_z_zzzz_xyyz[i] * fe_0 + ts_zzzz_xyyzz[i] * fe_0 + tr_z_zzzz_xyyzz[i] * pa_z[i];

        tr_z_zzzzz_xyzzz[i] =
            4.0 * tr_z_zzz_xyzzz[i] * fe_0 + 3.0 * tr_z_zzzz_xyzz[i] * fe_0 + ts_zzzz_xyzzz[i] * fe_0 + tr_z_zzzz_xyzzz[i] * pa_z[i];

        tr_z_zzzzz_xzzzz[i] =
            4.0 * tr_z_zzz_xzzzz[i] * fe_0 + 4.0 * tr_z_zzzz_xzzz[i] * fe_0 + ts_zzzz_xzzzz[i] * fe_0 + tr_z_zzzz_xzzzz[i] * pa_z[i];

        tr_z_zzzzz_yyyyy[i] = 4.0 * tr_z_zzz_yyyyy[i] * fe_0 + ts_zzzz_yyyyy[i] * fe_0 + tr_z_zzzz_yyyyy[i] * pa_z[i];

        tr_z_zzzzz_yyyyz[i] = 4.0 * tr_z_zzz_yyyyz[i] * fe_0 + tr_z_zzzz_yyyy[i] * fe_0 + ts_zzzz_yyyyz[i] * fe_0 + tr_z_zzzz_yyyyz[i] * pa_z[i];

        tr_z_zzzzz_yyyzz[i] =
            4.0 * tr_z_zzz_yyyzz[i] * fe_0 + 2.0 * tr_z_zzzz_yyyz[i] * fe_0 + ts_zzzz_yyyzz[i] * fe_0 + tr_z_zzzz_yyyzz[i] * pa_z[i];

        tr_z_zzzzz_yyzzz[i] =
            4.0 * tr_z_zzz_yyzzz[i] * fe_0 + 3.0 * tr_z_zzzz_yyzz[i] * fe_0 + ts_zzzz_yyzzz[i] * fe_0 + tr_z_zzzz_yyzzz[i] * pa_z[i];

        tr_z_zzzzz_yzzzz[i] =
            4.0 * tr_z_zzz_yzzzz[i] * fe_0 + 4.0 * tr_z_zzzz_yzzz[i] * fe_0 + ts_zzzz_yzzzz[i] * fe_0 + tr_z_zzzz_yzzzz[i] * pa_z[i];

        tr_z_zzzzz_zzzzz[i] =
            4.0 * tr_z_zzz_zzzzz[i] * fe_0 + 5.0 * tr_z_zzzz_zzzz[i] * fe_0 + ts_zzzz_zzzzz[i] * fe_0 + tr_z_zzzz_zzzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
