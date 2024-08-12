#include "ElectricDipoleMomentumPrimRecIF.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_if(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_if,
                                      const size_t idx_dip_gf,
                                      const size_t idx_dip_hd,
                                      const size_t idx_ovl_hf,
                                      const size_t idx_dip_hf,
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

    // Set up components of auxiliary buffer : GF

    auto tr_x_xxxx_xxx = pbuffer.data(idx_dip_gf);

    auto tr_x_xxxx_xxy = pbuffer.data(idx_dip_gf + 1);

    auto tr_x_xxxx_xxz = pbuffer.data(idx_dip_gf + 2);

    auto tr_x_xxxx_xyy = pbuffer.data(idx_dip_gf + 3);

    auto tr_x_xxxx_xyz = pbuffer.data(idx_dip_gf + 4);

    auto tr_x_xxxx_xzz = pbuffer.data(idx_dip_gf + 5);

    auto tr_x_xxxx_yyy = pbuffer.data(idx_dip_gf + 6);

    auto tr_x_xxxx_yyz = pbuffer.data(idx_dip_gf + 7);

    auto tr_x_xxxx_yzz = pbuffer.data(idx_dip_gf + 8);

    auto tr_x_xxxx_zzz = pbuffer.data(idx_dip_gf + 9);

    auto tr_x_xxxy_xxx = pbuffer.data(idx_dip_gf + 10);

    auto tr_x_xxxy_xxy = pbuffer.data(idx_dip_gf + 11);

    auto tr_x_xxxy_xxz = pbuffer.data(idx_dip_gf + 12);

    auto tr_x_xxxy_xyy = pbuffer.data(idx_dip_gf + 13);

    auto tr_x_xxxy_xyz = pbuffer.data(idx_dip_gf + 14);

    auto tr_x_xxxy_xzz = pbuffer.data(idx_dip_gf + 15);

    auto tr_x_xxxy_zzz = pbuffer.data(idx_dip_gf + 19);

    auto tr_x_xxxz_xxx = pbuffer.data(idx_dip_gf + 20);

    auto tr_x_xxxz_xxy = pbuffer.data(idx_dip_gf + 21);

    auto tr_x_xxxz_xxz = pbuffer.data(idx_dip_gf + 22);

    auto tr_x_xxxz_xyy = pbuffer.data(idx_dip_gf + 23);

    auto tr_x_xxxz_xyz = pbuffer.data(idx_dip_gf + 24);

    auto tr_x_xxxz_xzz = pbuffer.data(idx_dip_gf + 25);

    auto tr_x_xxxz_yyy = pbuffer.data(idx_dip_gf + 26);

    auto tr_x_xxxz_zzz = pbuffer.data(idx_dip_gf + 29);

    auto tr_x_xxyy_xxx = pbuffer.data(idx_dip_gf + 30);

    auto tr_x_xxyy_xxy = pbuffer.data(idx_dip_gf + 31);

    auto tr_x_xxyy_xxz = pbuffer.data(idx_dip_gf + 32);

    auto tr_x_xxyy_xyy = pbuffer.data(idx_dip_gf + 33);

    auto tr_x_xxyy_xyz = pbuffer.data(idx_dip_gf + 34);

    auto tr_x_xxyy_xzz = pbuffer.data(idx_dip_gf + 35);

    auto tr_x_xxyy_yyy = pbuffer.data(idx_dip_gf + 36);

    auto tr_x_xxyy_yyz = pbuffer.data(idx_dip_gf + 37);

    auto tr_x_xxyy_yzz = pbuffer.data(idx_dip_gf + 38);

    auto tr_x_xxyy_zzz = pbuffer.data(idx_dip_gf + 39);

    auto tr_x_xxyz_xxz = pbuffer.data(idx_dip_gf + 42);

    auto tr_x_xxyz_xzz = pbuffer.data(idx_dip_gf + 45);

    auto tr_x_xxyz_zzz = pbuffer.data(idx_dip_gf + 49);

    auto tr_x_xxzz_xxx = pbuffer.data(idx_dip_gf + 50);

    auto tr_x_xxzz_xxy = pbuffer.data(idx_dip_gf + 51);

    auto tr_x_xxzz_xxz = pbuffer.data(idx_dip_gf + 52);

    auto tr_x_xxzz_xyy = pbuffer.data(idx_dip_gf + 53);

    auto tr_x_xxzz_xyz = pbuffer.data(idx_dip_gf + 54);

    auto tr_x_xxzz_xzz = pbuffer.data(idx_dip_gf + 55);

    auto tr_x_xxzz_yyy = pbuffer.data(idx_dip_gf + 56);

    auto tr_x_xxzz_yyz = pbuffer.data(idx_dip_gf + 57);

    auto tr_x_xxzz_yzz = pbuffer.data(idx_dip_gf + 58);

    auto tr_x_xxzz_zzz = pbuffer.data(idx_dip_gf + 59);

    auto tr_x_xyyy_xxx = pbuffer.data(idx_dip_gf + 60);

    auto tr_x_xyyy_xxy = pbuffer.data(idx_dip_gf + 61);

    auto tr_x_xyyy_xxz = pbuffer.data(idx_dip_gf + 62);

    auto tr_x_xyyy_xyy = pbuffer.data(idx_dip_gf + 63);

    auto tr_x_xyyy_xzz = pbuffer.data(idx_dip_gf + 65);

    auto tr_x_xyyy_yyy = pbuffer.data(idx_dip_gf + 66);

    auto tr_x_xyyy_yyz = pbuffer.data(idx_dip_gf + 67);

    auto tr_x_xyyy_yzz = pbuffer.data(idx_dip_gf + 68);

    auto tr_x_xyyz_xxy = pbuffer.data(idx_dip_gf + 71);

    auto tr_x_xyyz_xxz = pbuffer.data(idx_dip_gf + 72);

    auto tr_x_xyyz_xyy = pbuffer.data(idx_dip_gf + 73);

    auto tr_x_xyyz_xzz = pbuffer.data(idx_dip_gf + 75);

    auto tr_x_xyzz_xxx = pbuffer.data(idx_dip_gf + 80);

    auto tr_x_xyzz_xxz = pbuffer.data(idx_dip_gf + 82);

    auto tr_x_xyzz_xzz = pbuffer.data(idx_dip_gf + 85);

    auto tr_x_xzzz_xxx = pbuffer.data(idx_dip_gf + 90);

    auto tr_x_xzzz_xxy = pbuffer.data(idx_dip_gf + 91);

    auto tr_x_xzzz_xxz = pbuffer.data(idx_dip_gf + 92);

    auto tr_x_xzzz_xyy = pbuffer.data(idx_dip_gf + 93);

    auto tr_x_xzzz_xzz = pbuffer.data(idx_dip_gf + 95);

    auto tr_x_xzzz_yyz = pbuffer.data(idx_dip_gf + 97);

    auto tr_x_xzzz_yzz = pbuffer.data(idx_dip_gf + 98);

    auto tr_x_xzzz_zzz = pbuffer.data(idx_dip_gf + 99);

    auto tr_x_yyyy_xxx = pbuffer.data(idx_dip_gf + 100);

    auto tr_x_yyyy_xxy = pbuffer.data(idx_dip_gf + 101);

    auto tr_x_yyyy_xxz = pbuffer.data(idx_dip_gf + 102);

    auto tr_x_yyyy_xyy = pbuffer.data(idx_dip_gf + 103);

    auto tr_x_yyyy_xyz = pbuffer.data(idx_dip_gf + 104);

    auto tr_x_yyyy_xzz = pbuffer.data(idx_dip_gf + 105);

    auto tr_x_yyyy_yyy = pbuffer.data(idx_dip_gf + 106);

    auto tr_x_yyyy_yyz = pbuffer.data(idx_dip_gf + 107);

    auto tr_x_yyyy_yzz = pbuffer.data(idx_dip_gf + 108);

    auto tr_x_yyyy_zzz = pbuffer.data(idx_dip_gf + 109);

    auto tr_x_yyyz_xxy = pbuffer.data(idx_dip_gf + 111);

    auto tr_x_yyyz_xxz = pbuffer.data(idx_dip_gf + 112);

    auto tr_x_yyyz_xyy = pbuffer.data(idx_dip_gf + 113);

    auto tr_x_yyyz_xzz = pbuffer.data(idx_dip_gf + 115);

    auto tr_x_yyyz_yyy = pbuffer.data(idx_dip_gf + 116);

    auto tr_x_yyyz_zzz = pbuffer.data(idx_dip_gf + 119);

    auto tr_x_yyzz_xxx = pbuffer.data(idx_dip_gf + 120);

    auto tr_x_yyzz_xxy = pbuffer.data(idx_dip_gf + 121);

    auto tr_x_yyzz_xxz = pbuffer.data(idx_dip_gf + 122);

    auto tr_x_yyzz_xyy = pbuffer.data(idx_dip_gf + 123);

    auto tr_x_yyzz_xyz = pbuffer.data(idx_dip_gf + 124);

    auto tr_x_yyzz_xzz = pbuffer.data(idx_dip_gf + 125);

    auto tr_x_yyzz_yyy = pbuffer.data(idx_dip_gf + 126);

    auto tr_x_yyzz_yyz = pbuffer.data(idx_dip_gf + 127);

    auto tr_x_yyzz_yzz = pbuffer.data(idx_dip_gf + 128);

    auto tr_x_yyzz_zzz = pbuffer.data(idx_dip_gf + 129);

    auto tr_x_yzzz_xxx = pbuffer.data(idx_dip_gf + 130);

    auto tr_x_yzzz_xxz = pbuffer.data(idx_dip_gf + 132);

    auto tr_x_yzzz_xyz = pbuffer.data(idx_dip_gf + 134);

    auto tr_x_yzzz_xzz = pbuffer.data(idx_dip_gf + 135);

    auto tr_x_yzzz_yyz = pbuffer.data(idx_dip_gf + 137);

    auto tr_x_yzzz_yzz = pbuffer.data(idx_dip_gf + 138);

    auto tr_x_yzzz_zzz = pbuffer.data(idx_dip_gf + 139);

    auto tr_x_zzzz_xxx = pbuffer.data(idx_dip_gf + 140);

    auto tr_x_zzzz_xxy = pbuffer.data(idx_dip_gf + 141);

    auto tr_x_zzzz_xxz = pbuffer.data(idx_dip_gf + 142);

    auto tr_x_zzzz_xyy = pbuffer.data(idx_dip_gf + 143);

    auto tr_x_zzzz_xyz = pbuffer.data(idx_dip_gf + 144);

    auto tr_x_zzzz_xzz = pbuffer.data(idx_dip_gf + 145);

    auto tr_x_zzzz_yyy = pbuffer.data(idx_dip_gf + 146);

    auto tr_x_zzzz_yyz = pbuffer.data(idx_dip_gf + 147);

    auto tr_x_zzzz_yzz = pbuffer.data(idx_dip_gf + 148);

    auto tr_x_zzzz_zzz = pbuffer.data(idx_dip_gf + 149);

    auto tr_y_xxxx_xxx = pbuffer.data(idx_dip_gf + 150);

    auto tr_y_xxxx_xxy = pbuffer.data(idx_dip_gf + 151);

    auto tr_y_xxxx_xxz = pbuffer.data(idx_dip_gf + 152);

    auto tr_y_xxxx_xyy = pbuffer.data(idx_dip_gf + 153);

    auto tr_y_xxxx_xyz = pbuffer.data(idx_dip_gf + 154);

    auto tr_y_xxxx_xzz = pbuffer.data(idx_dip_gf + 155);

    auto tr_y_xxxx_yyy = pbuffer.data(idx_dip_gf + 156);

    auto tr_y_xxxx_yyz = pbuffer.data(idx_dip_gf + 157);

    auto tr_y_xxxx_yzz = pbuffer.data(idx_dip_gf + 158);

    auto tr_y_xxxx_zzz = pbuffer.data(idx_dip_gf + 159);

    auto tr_y_xxxy_xxy = pbuffer.data(idx_dip_gf + 161);

    auto tr_y_xxxy_xyy = pbuffer.data(idx_dip_gf + 163);

    auto tr_y_xxxy_xyz = pbuffer.data(idx_dip_gf + 164);

    auto tr_y_xxxy_yyy = pbuffer.data(idx_dip_gf + 166);

    auto tr_y_xxxy_yyz = pbuffer.data(idx_dip_gf + 167);

    auto tr_y_xxxy_yzz = pbuffer.data(idx_dip_gf + 168);

    auto tr_y_xxxy_zzz = pbuffer.data(idx_dip_gf + 169);

    auto tr_y_xxxz_xxx = pbuffer.data(idx_dip_gf + 170);

    auto tr_y_xxxz_xxy = pbuffer.data(idx_dip_gf + 171);

    auto tr_y_xxxz_xyy = pbuffer.data(idx_dip_gf + 173);

    auto tr_y_xxxz_yyz = pbuffer.data(idx_dip_gf + 177);

    auto tr_y_xxxz_yzz = pbuffer.data(idx_dip_gf + 178);

    auto tr_y_xxxz_zzz = pbuffer.data(idx_dip_gf + 179);

    auto tr_y_xxyy_xxx = pbuffer.data(idx_dip_gf + 180);

    auto tr_y_xxyy_xxy = pbuffer.data(idx_dip_gf + 181);

    auto tr_y_xxyy_xxz = pbuffer.data(idx_dip_gf + 182);

    auto tr_y_xxyy_xyy = pbuffer.data(idx_dip_gf + 183);

    auto tr_y_xxyy_xyz = pbuffer.data(idx_dip_gf + 184);

    auto tr_y_xxyy_xzz = pbuffer.data(idx_dip_gf + 185);

    auto tr_y_xxyy_yyy = pbuffer.data(idx_dip_gf + 186);

    auto tr_y_xxyy_yyz = pbuffer.data(idx_dip_gf + 187);

    auto tr_y_xxyy_yzz = pbuffer.data(idx_dip_gf + 188);

    auto tr_y_xxyy_zzz = pbuffer.data(idx_dip_gf + 189);

    auto tr_y_xxyz_xxy = pbuffer.data(idx_dip_gf + 191);

    auto tr_y_xxyz_xyy = pbuffer.data(idx_dip_gf + 193);

    auto tr_y_xxyz_yyz = pbuffer.data(idx_dip_gf + 197);

    auto tr_y_xxyz_yzz = pbuffer.data(idx_dip_gf + 198);

    auto tr_y_xxyz_zzz = pbuffer.data(idx_dip_gf + 199);

    auto tr_y_xxzz_xxx = pbuffer.data(idx_dip_gf + 200);

    auto tr_y_xxzz_xxy = pbuffer.data(idx_dip_gf + 201);

    auto tr_y_xxzz_xxz = pbuffer.data(idx_dip_gf + 202);

    auto tr_y_xxzz_xyy = pbuffer.data(idx_dip_gf + 203);

    auto tr_y_xxzz_xyz = pbuffer.data(idx_dip_gf + 204);

    auto tr_y_xxzz_xzz = pbuffer.data(idx_dip_gf + 205);

    auto tr_y_xxzz_yyy = pbuffer.data(idx_dip_gf + 206);

    auto tr_y_xxzz_yyz = pbuffer.data(idx_dip_gf + 207);

    auto tr_y_xxzz_yzz = pbuffer.data(idx_dip_gf + 208);

    auto tr_y_xxzz_zzz = pbuffer.data(idx_dip_gf + 209);

    auto tr_y_xyyy_xxx = pbuffer.data(idx_dip_gf + 210);

    auto tr_y_xyyy_xxy = pbuffer.data(idx_dip_gf + 211);

    auto tr_y_xyyy_xxz = pbuffer.data(idx_dip_gf + 212);

    auto tr_y_xyyy_xyy = pbuffer.data(idx_dip_gf + 213);

    auto tr_y_xyyy_xyz = pbuffer.data(idx_dip_gf + 214);

    auto tr_y_xyyy_xzz = pbuffer.data(idx_dip_gf + 215);

    auto tr_y_xyyy_yyy = pbuffer.data(idx_dip_gf + 216);

    auto tr_y_xyyy_yyz = pbuffer.data(idx_dip_gf + 217);

    auto tr_y_xyyy_yzz = pbuffer.data(idx_dip_gf + 218);

    auto tr_y_xyyy_zzz = pbuffer.data(idx_dip_gf + 219);

    auto tr_y_xyyz_yyz = pbuffer.data(idx_dip_gf + 227);

    auto tr_y_xyyz_yzz = pbuffer.data(idx_dip_gf + 228);

    auto tr_y_xyyz_zzz = pbuffer.data(idx_dip_gf + 229);

    auto tr_y_xyzz_xyz = pbuffer.data(idx_dip_gf + 234);

    auto tr_y_xyzz_yyy = pbuffer.data(idx_dip_gf + 236);

    auto tr_y_xyzz_yyz = pbuffer.data(idx_dip_gf + 237);

    auto tr_y_xyzz_yzz = pbuffer.data(idx_dip_gf + 238);

    auto tr_y_xyzz_zzz = pbuffer.data(idx_dip_gf + 239);

    auto tr_y_xzzz_xxz = pbuffer.data(idx_dip_gf + 242);

    auto tr_y_xzzz_xyz = pbuffer.data(idx_dip_gf + 244);

    auto tr_y_xzzz_xzz = pbuffer.data(idx_dip_gf + 245);

    auto tr_y_xzzz_yyy = pbuffer.data(idx_dip_gf + 246);

    auto tr_y_xzzz_yyz = pbuffer.data(idx_dip_gf + 247);

    auto tr_y_xzzz_yzz = pbuffer.data(idx_dip_gf + 248);

    auto tr_y_xzzz_zzz = pbuffer.data(idx_dip_gf + 249);

    auto tr_y_yyyy_xxx = pbuffer.data(idx_dip_gf + 250);

    auto tr_y_yyyy_xxy = pbuffer.data(idx_dip_gf + 251);

    auto tr_y_yyyy_xxz = pbuffer.data(idx_dip_gf + 252);

    auto tr_y_yyyy_xyy = pbuffer.data(idx_dip_gf + 253);

    auto tr_y_yyyy_xyz = pbuffer.data(idx_dip_gf + 254);

    auto tr_y_yyyy_xzz = pbuffer.data(idx_dip_gf + 255);

    auto tr_y_yyyy_yyy = pbuffer.data(idx_dip_gf + 256);

    auto tr_y_yyyy_yyz = pbuffer.data(idx_dip_gf + 257);

    auto tr_y_yyyy_yzz = pbuffer.data(idx_dip_gf + 258);

    auto tr_y_yyyy_zzz = pbuffer.data(idx_dip_gf + 259);

    auto tr_y_yyyz_xxx = pbuffer.data(idx_dip_gf + 260);

    auto tr_y_yyyz_xxy = pbuffer.data(idx_dip_gf + 261);

    auto tr_y_yyyz_xyy = pbuffer.data(idx_dip_gf + 263);

    auto tr_y_yyyz_xyz = pbuffer.data(idx_dip_gf + 264);

    auto tr_y_yyyz_yyy = pbuffer.data(idx_dip_gf + 266);

    auto tr_y_yyyz_yyz = pbuffer.data(idx_dip_gf + 267);

    auto tr_y_yyyz_yzz = pbuffer.data(idx_dip_gf + 268);

    auto tr_y_yyyz_zzz = pbuffer.data(idx_dip_gf + 269);

    auto tr_y_yyzz_xxx = pbuffer.data(idx_dip_gf + 270);

    auto tr_y_yyzz_xxy = pbuffer.data(idx_dip_gf + 271);

    auto tr_y_yyzz_xxz = pbuffer.data(idx_dip_gf + 272);

    auto tr_y_yyzz_xyy = pbuffer.data(idx_dip_gf + 273);

    auto tr_y_yyzz_xyz = pbuffer.data(idx_dip_gf + 274);

    auto tr_y_yyzz_xzz = pbuffer.data(idx_dip_gf + 275);

    auto tr_y_yyzz_yyy = pbuffer.data(idx_dip_gf + 276);

    auto tr_y_yyzz_yyz = pbuffer.data(idx_dip_gf + 277);

    auto tr_y_yyzz_yzz = pbuffer.data(idx_dip_gf + 278);

    auto tr_y_yyzz_zzz = pbuffer.data(idx_dip_gf + 279);

    auto tr_y_yzzz_xxy = pbuffer.data(idx_dip_gf + 281);

    auto tr_y_yzzz_xxz = pbuffer.data(idx_dip_gf + 282);

    auto tr_y_yzzz_xyy = pbuffer.data(idx_dip_gf + 283);

    auto tr_y_yzzz_xyz = pbuffer.data(idx_dip_gf + 284);

    auto tr_y_yzzz_xzz = pbuffer.data(idx_dip_gf + 285);

    auto tr_y_yzzz_yyy = pbuffer.data(idx_dip_gf + 286);

    auto tr_y_yzzz_yyz = pbuffer.data(idx_dip_gf + 287);

    auto tr_y_yzzz_yzz = pbuffer.data(idx_dip_gf + 288);

    auto tr_y_yzzz_zzz = pbuffer.data(idx_dip_gf + 289);

    auto tr_y_zzzz_xxx = pbuffer.data(idx_dip_gf + 290);

    auto tr_y_zzzz_xxy = pbuffer.data(idx_dip_gf + 291);

    auto tr_y_zzzz_xxz = pbuffer.data(idx_dip_gf + 292);

    auto tr_y_zzzz_xyy = pbuffer.data(idx_dip_gf + 293);

    auto tr_y_zzzz_xyz = pbuffer.data(idx_dip_gf + 294);

    auto tr_y_zzzz_xzz = pbuffer.data(idx_dip_gf + 295);

    auto tr_y_zzzz_yyy = pbuffer.data(idx_dip_gf + 296);

    auto tr_y_zzzz_yyz = pbuffer.data(idx_dip_gf + 297);

    auto tr_y_zzzz_yzz = pbuffer.data(idx_dip_gf + 298);

    auto tr_y_zzzz_zzz = pbuffer.data(idx_dip_gf + 299);

    auto tr_z_xxxx_xxx = pbuffer.data(idx_dip_gf + 300);

    auto tr_z_xxxx_xxy = pbuffer.data(idx_dip_gf + 301);

    auto tr_z_xxxx_xxz = pbuffer.data(idx_dip_gf + 302);

    auto tr_z_xxxx_xyy = pbuffer.data(idx_dip_gf + 303);

    auto tr_z_xxxx_xyz = pbuffer.data(idx_dip_gf + 304);

    auto tr_z_xxxx_xzz = pbuffer.data(idx_dip_gf + 305);

    auto tr_z_xxxx_yyy = pbuffer.data(idx_dip_gf + 306);

    auto tr_z_xxxx_yyz = pbuffer.data(idx_dip_gf + 307);

    auto tr_z_xxxx_yzz = pbuffer.data(idx_dip_gf + 308);

    auto tr_z_xxxx_zzz = pbuffer.data(idx_dip_gf + 309);

    auto tr_z_xxxy_xxx = pbuffer.data(idx_dip_gf + 310);

    auto tr_z_xxxy_xxz = pbuffer.data(idx_dip_gf + 312);

    auto tr_z_xxxy_xzz = pbuffer.data(idx_dip_gf + 315);

    auto tr_z_xxxy_yyy = pbuffer.data(idx_dip_gf + 316);

    auto tr_z_xxxy_yyz = pbuffer.data(idx_dip_gf + 317);

    auto tr_z_xxxy_yzz = pbuffer.data(idx_dip_gf + 318);

    auto tr_z_xxxz_xxx = pbuffer.data(idx_dip_gf + 320);

    auto tr_z_xxxz_xxz = pbuffer.data(idx_dip_gf + 322);

    auto tr_z_xxxz_xyz = pbuffer.data(idx_dip_gf + 324);

    auto tr_z_xxxz_xzz = pbuffer.data(idx_dip_gf + 325);

    auto tr_z_xxxz_yyy = pbuffer.data(idx_dip_gf + 326);

    auto tr_z_xxxz_yyz = pbuffer.data(idx_dip_gf + 327);

    auto tr_z_xxxz_yzz = pbuffer.data(idx_dip_gf + 328);

    auto tr_z_xxxz_zzz = pbuffer.data(idx_dip_gf + 329);

    auto tr_z_xxyy_xxx = pbuffer.data(idx_dip_gf + 330);

    auto tr_z_xxyy_xxy = pbuffer.data(idx_dip_gf + 331);

    auto tr_z_xxyy_xxz = pbuffer.data(idx_dip_gf + 332);

    auto tr_z_xxyy_xyy = pbuffer.data(idx_dip_gf + 333);

    auto tr_z_xxyy_xyz = pbuffer.data(idx_dip_gf + 334);

    auto tr_z_xxyy_xzz = pbuffer.data(idx_dip_gf + 335);

    auto tr_z_xxyy_yyy = pbuffer.data(idx_dip_gf + 336);

    auto tr_z_xxyy_yyz = pbuffer.data(idx_dip_gf + 337);

    auto tr_z_xxyy_yzz = pbuffer.data(idx_dip_gf + 338);

    auto tr_z_xxyy_zzz = pbuffer.data(idx_dip_gf + 339);

    auto tr_z_xxyz_xxx = pbuffer.data(idx_dip_gf + 340);

    auto tr_z_xxyz_xxz = pbuffer.data(idx_dip_gf + 342);

    auto tr_z_xxyz_xzz = pbuffer.data(idx_dip_gf + 345);

    auto tr_z_xxyz_yyy = pbuffer.data(idx_dip_gf + 346);

    auto tr_z_xxyz_yyz = pbuffer.data(idx_dip_gf + 347);

    auto tr_z_xxyz_yzz = pbuffer.data(idx_dip_gf + 348);

    auto tr_z_xxzz_xxx = pbuffer.data(idx_dip_gf + 350);

    auto tr_z_xxzz_xxy = pbuffer.data(idx_dip_gf + 351);

    auto tr_z_xxzz_xxz = pbuffer.data(idx_dip_gf + 352);

    auto tr_z_xxzz_xyy = pbuffer.data(idx_dip_gf + 353);

    auto tr_z_xxzz_xyz = pbuffer.data(idx_dip_gf + 354);

    auto tr_z_xxzz_xzz = pbuffer.data(idx_dip_gf + 355);

    auto tr_z_xxzz_yyy = pbuffer.data(idx_dip_gf + 356);

    auto tr_z_xxzz_yyz = pbuffer.data(idx_dip_gf + 357);

    auto tr_z_xxzz_yzz = pbuffer.data(idx_dip_gf + 358);

    auto tr_z_xxzz_zzz = pbuffer.data(idx_dip_gf + 359);

    auto tr_z_xyyy_xxy = pbuffer.data(idx_dip_gf + 361);

    auto tr_z_xyyy_xyy = pbuffer.data(idx_dip_gf + 363);

    auto tr_z_xyyy_xyz = pbuffer.data(idx_dip_gf + 364);

    auto tr_z_xyyy_yyy = pbuffer.data(idx_dip_gf + 366);

    auto tr_z_xyyy_yyz = pbuffer.data(idx_dip_gf + 367);

    auto tr_z_xyyy_yzz = pbuffer.data(idx_dip_gf + 368);

    auto tr_z_xyyy_zzz = pbuffer.data(idx_dip_gf + 369);

    auto tr_z_xyyz_xyz = pbuffer.data(idx_dip_gf + 374);

    auto tr_z_xyyz_yyy = pbuffer.data(idx_dip_gf + 376);

    auto tr_z_xyyz_yyz = pbuffer.data(idx_dip_gf + 377);

    auto tr_z_xyyz_yzz = pbuffer.data(idx_dip_gf + 378);

    auto tr_z_xyyz_zzz = pbuffer.data(idx_dip_gf + 379);

    auto tr_z_xyzz_yyy = pbuffer.data(idx_dip_gf + 386);

    auto tr_z_xyzz_yyz = pbuffer.data(idx_dip_gf + 387);

    auto tr_z_xyzz_yzz = pbuffer.data(idx_dip_gf + 388);

    auto tr_z_xzzz_xxx = pbuffer.data(idx_dip_gf + 390);

    auto tr_z_xzzz_xxy = pbuffer.data(idx_dip_gf + 391);

    auto tr_z_xzzz_xxz = pbuffer.data(idx_dip_gf + 392);

    auto tr_z_xzzz_xyy = pbuffer.data(idx_dip_gf + 393);

    auto tr_z_xzzz_xyz = pbuffer.data(idx_dip_gf + 394);

    auto tr_z_xzzz_xzz = pbuffer.data(idx_dip_gf + 395);

    auto tr_z_xzzz_yyy = pbuffer.data(idx_dip_gf + 396);

    auto tr_z_xzzz_yyz = pbuffer.data(idx_dip_gf + 397);

    auto tr_z_xzzz_yzz = pbuffer.data(idx_dip_gf + 398);

    auto tr_z_xzzz_zzz = pbuffer.data(idx_dip_gf + 399);

    auto tr_z_yyyy_xxx = pbuffer.data(idx_dip_gf + 400);

    auto tr_z_yyyy_xxy = pbuffer.data(idx_dip_gf + 401);

    auto tr_z_yyyy_xxz = pbuffer.data(idx_dip_gf + 402);

    auto tr_z_yyyy_xyy = pbuffer.data(idx_dip_gf + 403);

    auto tr_z_yyyy_xyz = pbuffer.data(idx_dip_gf + 404);

    auto tr_z_yyyy_xzz = pbuffer.data(idx_dip_gf + 405);

    auto tr_z_yyyy_yyy = pbuffer.data(idx_dip_gf + 406);

    auto tr_z_yyyy_yyz = pbuffer.data(idx_dip_gf + 407);

    auto tr_z_yyyy_yzz = pbuffer.data(idx_dip_gf + 408);

    auto tr_z_yyyy_zzz = pbuffer.data(idx_dip_gf + 409);

    auto tr_z_yyyz_xxx = pbuffer.data(idx_dip_gf + 410);

    auto tr_z_yyyz_xxz = pbuffer.data(idx_dip_gf + 412);

    auto tr_z_yyyz_xyz = pbuffer.data(idx_dip_gf + 414);

    auto tr_z_yyyz_xzz = pbuffer.data(idx_dip_gf + 415);

    auto tr_z_yyyz_yyy = pbuffer.data(idx_dip_gf + 416);

    auto tr_z_yyyz_yyz = pbuffer.data(idx_dip_gf + 417);

    auto tr_z_yyyz_yzz = pbuffer.data(idx_dip_gf + 418);

    auto tr_z_yyyz_zzz = pbuffer.data(idx_dip_gf + 419);

    auto tr_z_yyzz_xxx = pbuffer.data(idx_dip_gf + 420);

    auto tr_z_yyzz_xxy = pbuffer.data(idx_dip_gf + 421);

    auto tr_z_yyzz_xxz = pbuffer.data(idx_dip_gf + 422);

    auto tr_z_yyzz_xyy = pbuffer.data(idx_dip_gf + 423);

    auto tr_z_yyzz_xyz = pbuffer.data(idx_dip_gf + 424);

    auto tr_z_yyzz_xzz = pbuffer.data(idx_dip_gf + 425);

    auto tr_z_yyzz_yyy = pbuffer.data(idx_dip_gf + 426);

    auto tr_z_yyzz_yyz = pbuffer.data(idx_dip_gf + 427);

    auto tr_z_yyzz_yzz = pbuffer.data(idx_dip_gf + 428);

    auto tr_z_yyzz_zzz = pbuffer.data(idx_dip_gf + 429);

    auto tr_z_yzzz_xxx = pbuffer.data(idx_dip_gf + 430);

    auto tr_z_yzzz_xxy = pbuffer.data(idx_dip_gf + 431);

    auto tr_z_yzzz_xxz = pbuffer.data(idx_dip_gf + 432);

    auto tr_z_yzzz_xyy = pbuffer.data(idx_dip_gf + 433);

    auto tr_z_yzzz_xyz = pbuffer.data(idx_dip_gf + 434);

    auto tr_z_yzzz_xzz = pbuffer.data(idx_dip_gf + 435);

    auto tr_z_yzzz_yyy = pbuffer.data(idx_dip_gf + 436);

    auto tr_z_yzzz_yyz = pbuffer.data(idx_dip_gf + 437);

    auto tr_z_yzzz_yzz = pbuffer.data(idx_dip_gf + 438);

    auto tr_z_yzzz_zzz = pbuffer.data(idx_dip_gf + 439);

    auto tr_z_zzzz_xxx = pbuffer.data(idx_dip_gf + 440);

    auto tr_z_zzzz_xxy = pbuffer.data(idx_dip_gf + 441);

    auto tr_z_zzzz_xxz = pbuffer.data(idx_dip_gf + 442);

    auto tr_z_zzzz_xyy = pbuffer.data(idx_dip_gf + 443);

    auto tr_z_zzzz_xyz = pbuffer.data(idx_dip_gf + 444);

    auto tr_z_zzzz_xzz = pbuffer.data(idx_dip_gf + 445);

    auto tr_z_zzzz_yyy = pbuffer.data(idx_dip_gf + 446);

    auto tr_z_zzzz_yyz = pbuffer.data(idx_dip_gf + 447);

    auto tr_z_zzzz_yzz = pbuffer.data(idx_dip_gf + 448);

    auto tr_z_zzzz_zzz = pbuffer.data(idx_dip_gf + 449);

    // Set up components of auxiliary buffer : HD

    auto tr_x_xxxxx_xx = pbuffer.data(idx_dip_hd);

    auto tr_x_xxxxx_xy = pbuffer.data(idx_dip_hd + 1);

    auto tr_x_xxxxx_xz = pbuffer.data(idx_dip_hd + 2);

    auto tr_x_xxxxx_yy = pbuffer.data(idx_dip_hd + 3);

    auto tr_x_xxxxx_yz = pbuffer.data(idx_dip_hd + 4);

    auto tr_x_xxxxx_zz = pbuffer.data(idx_dip_hd + 5);

    auto tr_x_xxxxy_xx = pbuffer.data(idx_dip_hd + 6);

    auto tr_x_xxxxy_xy = pbuffer.data(idx_dip_hd + 7);

    auto tr_x_xxxxy_xz = pbuffer.data(idx_dip_hd + 8);

    auto tr_x_xxxxz_xx = pbuffer.data(idx_dip_hd + 12);

    auto tr_x_xxxxz_xy = pbuffer.data(idx_dip_hd + 13);

    auto tr_x_xxxxz_xz = pbuffer.data(idx_dip_hd + 14);

    auto tr_x_xxxxz_yz = pbuffer.data(idx_dip_hd + 16);

    auto tr_x_xxxxz_zz = pbuffer.data(idx_dip_hd + 17);

    auto tr_x_xxxyy_xx = pbuffer.data(idx_dip_hd + 18);

    auto tr_x_xxxyy_xy = pbuffer.data(idx_dip_hd + 19);

    auto tr_x_xxxyy_xz = pbuffer.data(idx_dip_hd + 20);

    auto tr_x_xxxyy_yy = pbuffer.data(idx_dip_hd + 21);

    auto tr_x_xxxyy_yz = pbuffer.data(idx_dip_hd + 22);

    auto tr_x_xxxzz_xx = pbuffer.data(idx_dip_hd + 30);

    auto tr_x_xxxzz_xy = pbuffer.data(idx_dip_hd + 31);

    auto tr_x_xxxzz_xz = pbuffer.data(idx_dip_hd + 32);

    auto tr_x_xxxzz_yy = pbuffer.data(idx_dip_hd + 33);

    auto tr_x_xxxzz_yz = pbuffer.data(idx_dip_hd + 34);

    auto tr_x_xxxzz_zz = pbuffer.data(idx_dip_hd + 35);

    auto tr_x_xxyyy_xx = pbuffer.data(idx_dip_hd + 36);

    auto tr_x_xxyyy_xy = pbuffer.data(idx_dip_hd + 37);

    auto tr_x_xxyyy_xz = pbuffer.data(idx_dip_hd + 38);

    auto tr_x_xxyyy_yy = pbuffer.data(idx_dip_hd + 39);

    auto tr_x_xxyyy_yz = pbuffer.data(idx_dip_hd + 40);

    auto tr_x_xxyzz_xz = pbuffer.data(idx_dip_hd + 50);

    auto tr_x_xxzzz_xx = pbuffer.data(idx_dip_hd + 54);

    auto tr_x_xxzzz_xy = pbuffer.data(idx_dip_hd + 55);

    auto tr_x_xxzzz_xz = pbuffer.data(idx_dip_hd + 56);

    auto tr_x_xxzzz_yy = pbuffer.data(idx_dip_hd + 57);

    auto tr_x_xxzzz_yz = pbuffer.data(idx_dip_hd + 58);

    auto tr_x_xxzzz_zz = pbuffer.data(idx_dip_hd + 59);

    auto tr_x_xyyyy_xy = pbuffer.data(idx_dip_hd + 61);

    auto tr_x_xzzzz_xx = pbuffer.data(idx_dip_hd + 84);

    auto tr_x_xzzzz_xy = pbuffer.data(idx_dip_hd + 85);

    auto tr_x_xzzzz_xz = pbuffer.data(idx_dip_hd + 86);

    auto tr_x_yyyyy_xx = pbuffer.data(idx_dip_hd + 90);

    auto tr_x_yyyyy_xy = pbuffer.data(idx_dip_hd + 91);

    auto tr_x_yyyyy_xz = pbuffer.data(idx_dip_hd + 92);

    auto tr_x_yyyyy_yy = pbuffer.data(idx_dip_hd + 93);

    auto tr_x_yyyyy_yz = pbuffer.data(idx_dip_hd + 94);

    auto tr_x_yyyyy_zz = pbuffer.data(idx_dip_hd + 95);

    auto tr_x_yyyzz_xz = pbuffer.data(idx_dip_hd + 104);

    auto tr_x_yyyzz_yz = pbuffer.data(idx_dip_hd + 106);

    auto tr_x_yyyzz_zz = pbuffer.data(idx_dip_hd + 107);

    auto tr_x_yyzzz_xz = pbuffer.data(idx_dip_hd + 110);

    auto tr_x_yyzzz_yz = pbuffer.data(idx_dip_hd + 112);

    auto tr_x_yyzzz_zz = pbuffer.data(idx_dip_hd + 113);

    auto tr_x_yzzzz_xz = pbuffer.data(idx_dip_hd + 116);

    auto tr_x_yzzzz_yz = pbuffer.data(idx_dip_hd + 118);

    auto tr_x_yzzzz_zz = pbuffer.data(idx_dip_hd + 119);

    auto tr_x_zzzzz_xx = pbuffer.data(idx_dip_hd + 120);

    auto tr_x_zzzzz_xy = pbuffer.data(idx_dip_hd + 121);

    auto tr_x_zzzzz_xz = pbuffer.data(idx_dip_hd + 122);

    auto tr_x_zzzzz_yy = pbuffer.data(idx_dip_hd + 123);

    auto tr_x_zzzzz_yz = pbuffer.data(idx_dip_hd + 124);

    auto tr_x_zzzzz_zz = pbuffer.data(idx_dip_hd + 125);

    auto tr_y_xxxxx_xx = pbuffer.data(idx_dip_hd + 126);

    auto tr_y_xxxxx_xy = pbuffer.data(idx_dip_hd + 127);

    auto tr_y_xxxxx_xz = pbuffer.data(idx_dip_hd + 128);

    auto tr_y_xxxxx_yy = pbuffer.data(idx_dip_hd + 129);

    auto tr_y_xxxxx_yz = pbuffer.data(idx_dip_hd + 130);

    auto tr_y_xxxxx_zz = pbuffer.data(idx_dip_hd + 131);

    auto tr_y_xxxxy_xy = pbuffer.data(idx_dip_hd + 133);

    auto tr_y_xxxxy_yy = pbuffer.data(idx_dip_hd + 135);

    auto tr_y_xxxxy_yz = pbuffer.data(idx_dip_hd + 136);

    auto tr_y_xxxyy_xx = pbuffer.data(idx_dip_hd + 144);

    auto tr_y_xxxyy_xy = pbuffer.data(idx_dip_hd + 145);

    auto tr_y_xxxyy_xz = pbuffer.data(idx_dip_hd + 146);

    auto tr_y_xxxyy_yy = pbuffer.data(idx_dip_hd + 147);

    auto tr_y_xxxyy_yz = pbuffer.data(idx_dip_hd + 148);

    auto tr_y_xxxyy_zz = pbuffer.data(idx_dip_hd + 149);

    auto tr_y_xxxzz_xz = pbuffer.data(idx_dip_hd + 158);

    auto tr_y_xxxzz_yz = pbuffer.data(idx_dip_hd + 160);

    auto tr_y_xxxzz_zz = pbuffer.data(idx_dip_hd + 161);

    auto tr_y_xxyyy_xx = pbuffer.data(idx_dip_hd + 162);

    auto tr_y_xxyyy_xy = pbuffer.data(idx_dip_hd + 163);

    auto tr_y_xxyyy_xz = pbuffer.data(idx_dip_hd + 164);

    auto tr_y_xxyyy_yy = pbuffer.data(idx_dip_hd + 165);

    auto tr_y_xxyyy_yz = pbuffer.data(idx_dip_hd + 166);

    auto tr_y_xxyyy_zz = pbuffer.data(idx_dip_hd + 167);

    auto tr_y_xxyzz_yz = pbuffer.data(idx_dip_hd + 178);

    auto tr_y_xxzzz_xz = pbuffer.data(idx_dip_hd + 182);

    auto tr_y_xxzzz_yz = pbuffer.data(idx_dip_hd + 184);

    auto tr_y_xxzzz_zz = pbuffer.data(idx_dip_hd + 185);

    auto tr_y_xyyyy_xx = pbuffer.data(idx_dip_hd + 186);

    auto tr_y_xyyyy_xy = pbuffer.data(idx_dip_hd + 187);

    auto tr_y_xyyyy_xz = pbuffer.data(idx_dip_hd + 188);

    auto tr_y_xyyyy_yy = pbuffer.data(idx_dip_hd + 189);

    auto tr_y_xyyyy_yz = pbuffer.data(idx_dip_hd + 190);

    auto tr_y_xyyyy_zz = pbuffer.data(idx_dip_hd + 191);

    auto tr_y_xyyzz_xz = pbuffer.data(idx_dip_hd + 200);

    auto tr_y_xyyzz_yz = pbuffer.data(idx_dip_hd + 202);

    auto tr_y_xyyzz_zz = pbuffer.data(idx_dip_hd + 203);

    auto tr_y_xyzzz_yz = pbuffer.data(idx_dip_hd + 208);

    auto tr_y_xzzzz_xz = pbuffer.data(idx_dip_hd + 212);

    auto tr_y_xzzzz_yz = pbuffer.data(idx_dip_hd + 214);

    auto tr_y_xzzzz_zz = pbuffer.data(idx_dip_hd + 215);

    auto tr_y_yyyyy_xx = pbuffer.data(idx_dip_hd + 216);

    auto tr_y_yyyyy_xy = pbuffer.data(idx_dip_hd + 217);

    auto tr_y_yyyyy_xz = pbuffer.data(idx_dip_hd + 218);

    auto tr_y_yyyyy_yy = pbuffer.data(idx_dip_hd + 219);

    auto tr_y_yyyyy_yz = pbuffer.data(idx_dip_hd + 220);

    auto tr_y_yyyyy_zz = pbuffer.data(idx_dip_hd + 221);

    auto tr_y_yyyyz_xy = pbuffer.data(idx_dip_hd + 223);

    auto tr_y_yyyyz_xz = pbuffer.data(idx_dip_hd + 224);

    auto tr_y_yyyyz_yy = pbuffer.data(idx_dip_hd + 225);

    auto tr_y_yyyyz_yz = pbuffer.data(idx_dip_hd + 226);

    auto tr_y_yyyyz_zz = pbuffer.data(idx_dip_hd + 227);

    auto tr_y_yyyzz_xx = pbuffer.data(idx_dip_hd + 228);

    auto tr_y_yyyzz_xy = pbuffer.data(idx_dip_hd + 229);

    auto tr_y_yyyzz_xz = pbuffer.data(idx_dip_hd + 230);

    auto tr_y_yyyzz_yy = pbuffer.data(idx_dip_hd + 231);

    auto tr_y_yyyzz_yz = pbuffer.data(idx_dip_hd + 232);

    auto tr_y_yyyzz_zz = pbuffer.data(idx_dip_hd + 233);

    auto tr_y_yyzzz_xx = pbuffer.data(idx_dip_hd + 234);

    auto tr_y_yyzzz_xy = pbuffer.data(idx_dip_hd + 235);

    auto tr_y_yyzzz_xz = pbuffer.data(idx_dip_hd + 236);

    auto tr_y_yyzzz_yy = pbuffer.data(idx_dip_hd + 237);

    auto tr_y_yyzzz_yz = pbuffer.data(idx_dip_hd + 238);

    auto tr_y_yyzzz_zz = pbuffer.data(idx_dip_hd + 239);

    auto tr_y_yzzzz_xx = pbuffer.data(idx_dip_hd + 240);

    auto tr_y_yzzzz_xy = pbuffer.data(idx_dip_hd + 241);

    auto tr_y_yzzzz_xz = pbuffer.data(idx_dip_hd + 242);

    auto tr_y_yzzzz_yy = pbuffer.data(idx_dip_hd + 243);

    auto tr_y_yzzzz_yz = pbuffer.data(idx_dip_hd + 244);

    auto tr_y_yzzzz_zz = pbuffer.data(idx_dip_hd + 245);

    auto tr_y_zzzzz_xx = pbuffer.data(idx_dip_hd + 246);

    auto tr_y_zzzzz_xy = pbuffer.data(idx_dip_hd + 247);

    auto tr_y_zzzzz_xz = pbuffer.data(idx_dip_hd + 248);

    auto tr_y_zzzzz_yy = pbuffer.data(idx_dip_hd + 249);

    auto tr_y_zzzzz_yz = pbuffer.data(idx_dip_hd + 250);

    auto tr_y_zzzzz_zz = pbuffer.data(idx_dip_hd + 251);

    auto tr_z_xxxxx_xx = pbuffer.data(idx_dip_hd + 252);

    auto tr_z_xxxxx_xy = pbuffer.data(idx_dip_hd + 253);

    auto tr_z_xxxxx_xz = pbuffer.data(idx_dip_hd + 254);

    auto tr_z_xxxxx_yy = pbuffer.data(idx_dip_hd + 255);

    auto tr_z_xxxxx_yz = pbuffer.data(idx_dip_hd + 256);

    auto tr_z_xxxxx_zz = pbuffer.data(idx_dip_hd + 257);

    auto tr_z_xxxxz_xx = pbuffer.data(idx_dip_hd + 264);

    auto tr_z_xxxxz_xy = pbuffer.data(idx_dip_hd + 265);

    auto tr_z_xxxxz_xz = pbuffer.data(idx_dip_hd + 266);

    auto tr_z_xxxxz_yz = pbuffer.data(idx_dip_hd + 268);

    auto tr_z_xxxxz_zz = pbuffer.data(idx_dip_hd + 269);

    auto tr_z_xxxyy_xy = pbuffer.data(idx_dip_hd + 271);

    auto tr_z_xxxyy_yy = pbuffer.data(idx_dip_hd + 273);

    auto tr_z_xxxyy_yz = pbuffer.data(idx_dip_hd + 274);

    auto tr_z_xxxzz_xx = pbuffer.data(idx_dip_hd + 282);

    auto tr_z_xxxzz_xy = pbuffer.data(idx_dip_hd + 283);

    auto tr_z_xxxzz_xz = pbuffer.data(idx_dip_hd + 284);

    auto tr_z_xxxzz_yy = pbuffer.data(idx_dip_hd + 285);

    auto tr_z_xxxzz_yz = pbuffer.data(idx_dip_hd + 286);

    auto tr_z_xxxzz_zz = pbuffer.data(idx_dip_hd + 287);

    auto tr_z_xxyyy_xy = pbuffer.data(idx_dip_hd + 289);

    auto tr_z_xxyyy_yy = pbuffer.data(idx_dip_hd + 291);

    auto tr_z_xxyyy_yz = pbuffer.data(idx_dip_hd + 292);

    auto tr_z_xxyyz_yz = pbuffer.data(idx_dip_hd + 298);

    auto tr_z_xxzzz_xx = pbuffer.data(idx_dip_hd + 306);

    auto tr_z_xxzzz_xy = pbuffer.data(idx_dip_hd + 307);

    auto tr_z_xxzzz_xz = pbuffer.data(idx_dip_hd + 308);

    auto tr_z_xxzzz_yy = pbuffer.data(idx_dip_hd + 309);

    auto tr_z_xxzzz_yz = pbuffer.data(idx_dip_hd + 310);

    auto tr_z_xxzzz_zz = pbuffer.data(idx_dip_hd + 311);

    auto tr_z_xyyyy_xy = pbuffer.data(idx_dip_hd + 313);

    auto tr_z_xyyyy_yy = pbuffer.data(idx_dip_hd + 315);

    auto tr_z_xyyyy_yz = pbuffer.data(idx_dip_hd + 316);

    auto tr_z_xyyyz_yz = pbuffer.data(idx_dip_hd + 322);

    auto tr_z_xyyzz_xy = pbuffer.data(idx_dip_hd + 325);

    auto tr_z_xyyzz_yy = pbuffer.data(idx_dip_hd + 327);

    auto tr_z_xyyzz_yz = pbuffer.data(idx_dip_hd + 328);

    auto tr_z_xzzzz_xx = pbuffer.data(idx_dip_hd + 336);

    auto tr_z_xzzzz_xy = pbuffer.data(idx_dip_hd + 337);

    auto tr_z_xzzzz_xz = pbuffer.data(idx_dip_hd + 338);

    auto tr_z_xzzzz_yy = pbuffer.data(idx_dip_hd + 339);

    auto tr_z_xzzzz_yz = pbuffer.data(idx_dip_hd + 340);

    auto tr_z_xzzzz_zz = pbuffer.data(idx_dip_hd + 341);

    auto tr_z_yyyyy_xx = pbuffer.data(idx_dip_hd + 342);

    auto tr_z_yyyyy_xy = pbuffer.data(idx_dip_hd + 343);

    auto tr_z_yyyyy_xz = pbuffer.data(idx_dip_hd + 344);

    auto tr_z_yyyyy_yy = pbuffer.data(idx_dip_hd + 345);

    auto tr_z_yyyyy_yz = pbuffer.data(idx_dip_hd + 346);

    auto tr_z_yyyyy_zz = pbuffer.data(idx_dip_hd + 347);

    auto tr_z_yyyyz_xx = pbuffer.data(idx_dip_hd + 348);

    auto tr_z_yyyyz_xy = pbuffer.data(idx_dip_hd + 349);

    auto tr_z_yyyyz_xz = pbuffer.data(idx_dip_hd + 350);

    auto tr_z_yyyyz_yy = pbuffer.data(idx_dip_hd + 351);

    auto tr_z_yyyyz_yz = pbuffer.data(idx_dip_hd + 352);

    auto tr_z_yyyyz_zz = pbuffer.data(idx_dip_hd + 353);

    auto tr_z_yyyzz_xx = pbuffer.data(idx_dip_hd + 354);

    auto tr_z_yyyzz_xy = pbuffer.data(idx_dip_hd + 355);

    auto tr_z_yyyzz_xz = pbuffer.data(idx_dip_hd + 356);

    auto tr_z_yyyzz_yy = pbuffer.data(idx_dip_hd + 357);

    auto tr_z_yyyzz_yz = pbuffer.data(idx_dip_hd + 358);

    auto tr_z_yyyzz_zz = pbuffer.data(idx_dip_hd + 359);

    auto tr_z_yyzzz_xx = pbuffer.data(idx_dip_hd + 360);

    auto tr_z_yyzzz_xy = pbuffer.data(idx_dip_hd + 361);

    auto tr_z_yyzzz_xz = pbuffer.data(idx_dip_hd + 362);

    auto tr_z_yyzzz_yy = pbuffer.data(idx_dip_hd + 363);

    auto tr_z_yyzzz_yz = pbuffer.data(idx_dip_hd + 364);

    auto tr_z_yyzzz_zz = pbuffer.data(idx_dip_hd + 365);

    auto tr_z_yzzzz_xx = pbuffer.data(idx_dip_hd + 366);

    auto tr_z_yzzzz_xy = pbuffer.data(idx_dip_hd + 367);

    auto tr_z_yzzzz_xz = pbuffer.data(idx_dip_hd + 368);

    auto tr_z_yzzzz_yy = pbuffer.data(idx_dip_hd + 369);

    auto tr_z_yzzzz_yz = pbuffer.data(idx_dip_hd + 370);

    auto tr_z_yzzzz_zz = pbuffer.data(idx_dip_hd + 371);

    auto tr_z_zzzzz_xx = pbuffer.data(idx_dip_hd + 372);

    auto tr_z_zzzzz_xy = pbuffer.data(idx_dip_hd + 373);

    auto tr_z_zzzzz_xz = pbuffer.data(idx_dip_hd + 374);

    auto tr_z_zzzzz_yy = pbuffer.data(idx_dip_hd + 375);

    auto tr_z_zzzzz_yz = pbuffer.data(idx_dip_hd + 376);

    auto tr_z_zzzzz_zz = pbuffer.data(idx_dip_hd + 377);

    // Set up components of auxiliary buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_ovl_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_ovl_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_ovl_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_ovl_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_ovl_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_ovl_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_ovl_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_ovl_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_ovl_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_ovl_hf + 9);

    auto ts_xxxxz_xxz = pbuffer.data(idx_ovl_hf + 22);

    auto ts_xxxxz_xzz = pbuffer.data(idx_ovl_hf + 25);

    auto ts_xxxyy_xxy = pbuffer.data(idx_ovl_hf + 31);

    auto ts_xxxyy_xyy = pbuffer.data(idx_ovl_hf + 33);

    auto ts_xxxyy_yyy = pbuffer.data(idx_ovl_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_ovl_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_ovl_hf + 38);

    auto ts_xxxzz_xxx = pbuffer.data(idx_ovl_hf + 50);

    auto ts_xxxzz_xxz = pbuffer.data(idx_ovl_hf + 52);

    auto ts_xxxzz_xzz = pbuffer.data(idx_ovl_hf + 55);

    auto ts_xxxzz_yyz = pbuffer.data(idx_ovl_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_ovl_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_ovl_hf + 59);

    auto ts_xxyyy_xxy = pbuffer.data(idx_ovl_hf + 61);

    auto ts_xxyyy_xyy = pbuffer.data(idx_ovl_hf + 63);

    auto ts_xxyyy_yyy = pbuffer.data(idx_ovl_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_ovl_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_ovl_hf + 68);

    auto ts_xxzzz_xxx = pbuffer.data(idx_ovl_hf + 90);

    auto ts_xxzzz_xxz = pbuffer.data(idx_ovl_hf + 92);

    auto ts_xxzzz_xzz = pbuffer.data(idx_ovl_hf + 95);

    auto ts_xxzzz_yyz = pbuffer.data(idx_ovl_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_ovl_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_ovl_hf + 99);

    auto ts_xyyyy_yyy = pbuffer.data(idx_ovl_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_ovl_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_ovl_hf + 108);

    auto ts_xyyzz_yyz = pbuffer.data(idx_ovl_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_ovl_hf + 128);

    auto ts_xzzzz_yyz = pbuffer.data(idx_ovl_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_ovl_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_ovl_hf + 149);

    auto ts_yyyyy_xxx = pbuffer.data(idx_ovl_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_ovl_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_ovl_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_ovl_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_ovl_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_ovl_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_ovl_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_ovl_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_ovl_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_ovl_hf + 159);

    auto ts_yyyyz_yyz = pbuffer.data(idx_ovl_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_ovl_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_ovl_hf + 169);

    auto ts_yyyzz_xxz = pbuffer.data(idx_ovl_hf + 172);

    auto ts_yyyzz_xyz = pbuffer.data(idx_ovl_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_ovl_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_ovl_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_ovl_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_ovl_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_ovl_hf + 179);

    auto ts_yyzzz_xxz = pbuffer.data(idx_ovl_hf + 182);

    auto ts_yyzzz_xyz = pbuffer.data(idx_ovl_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_ovl_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_ovl_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_ovl_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_ovl_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_ovl_hf + 189);

    auto ts_yzzzz_xxz = pbuffer.data(idx_ovl_hf + 192);

    auto ts_yzzzz_xzz = pbuffer.data(idx_ovl_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_ovl_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_ovl_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_ovl_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_ovl_hf + 199);

    auto ts_zzzzz_xxx = pbuffer.data(idx_ovl_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_ovl_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_ovl_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_ovl_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_ovl_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_ovl_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_ovl_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_ovl_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_ovl_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_ovl_hf + 209);

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

    auto tr_x_xxxxy_yyy = pbuffer.data(idx_dip_hf + 16);

    auto tr_x_xxxxy_zzz = pbuffer.data(idx_dip_hf + 19);

    auto tr_x_xxxxz_xxx = pbuffer.data(idx_dip_hf + 20);

    auto tr_x_xxxxz_xxy = pbuffer.data(idx_dip_hf + 21);

    auto tr_x_xxxxz_xxz = pbuffer.data(idx_dip_hf + 22);

    auto tr_x_xxxxz_xyy = pbuffer.data(idx_dip_hf + 23);

    auto tr_x_xxxxz_xyz = pbuffer.data(idx_dip_hf + 24);

    auto tr_x_xxxxz_xzz = pbuffer.data(idx_dip_hf + 25);

    auto tr_x_xxxxz_yyy = pbuffer.data(idx_dip_hf + 26);

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

    auto tr_x_xxxyy_zzz = pbuffer.data(idx_dip_hf + 39);

    auto tr_x_xxxyz_xxz = pbuffer.data(idx_dip_hf + 42);

    auto tr_x_xxxyz_xzz = pbuffer.data(idx_dip_hf + 45);

    auto tr_x_xxxyz_zzz = pbuffer.data(idx_dip_hf + 49);

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

    auto tr_x_xxyyy_zzz = pbuffer.data(idx_dip_hf + 69);

    auto tr_x_xxyyz_xxy = pbuffer.data(idx_dip_hf + 71);

    auto tr_x_xxyyz_xxz = pbuffer.data(idx_dip_hf + 72);

    auto tr_x_xxyyz_xyy = pbuffer.data(idx_dip_hf + 73);

    auto tr_x_xxyyz_xzz = pbuffer.data(idx_dip_hf + 75);

    auto tr_x_xxyyz_yyy = pbuffer.data(idx_dip_hf + 76);

    auto tr_x_xxyyz_zzz = pbuffer.data(idx_dip_hf + 79);

    auto tr_x_xxyzz_xxx = pbuffer.data(idx_dip_hf + 80);

    auto tr_x_xxyzz_xxz = pbuffer.data(idx_dip_hf + 82);

    auto tr_x_xxyzz_xyz = pbuffer.data(idx_dip_hf + 84);

    auto tr_x_xxyzz_xzz = pbuffer.data(idx_dip_hf + 85);

    auto tr_x_xxyzz_zzz = pbuffer.data(idx_dip_hf + 89);

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

    auto tr_x_xyyyy_xxx = pbuffer.data(idx_dip_hf + 100);

    auto tr_x_xyyyy_xxy = pbuffer.data(idx_dip_hf + 101);

    auto tr_x_xyyyy_xxz = pbuffer.data(idx_dip_hf + 102);

    auto tr_x_xyyyy_xyy = pbuffer.data(idx_dip_hf + 103);

    auto tr_x_xyyyy_xyz = pbuffer.data(idx_dip_hf + 104);

    auto tr_x_xyyyy_xzz = pbuffer.data(idx_dip_hf + 105);

    auto tr_x_xyyyy_yyy = pbuffer.data(idx_dip_hf + 106);

    auto tr_x_xyyyy_yyz = pbuffer.data(idx_dip_hf + 107);

    auto tr_x_xyyyy_yzz = pbuffer.data(idx_dip_hf + 108);

    auto tr_x_xyyyz_xxy = pbuffer.data(idx_dip_hf + 111);

    auto tr_x_xyyyz_xxz = pbuffer.data(idx_dip_hf + 112);

    auto tr_x_xyyyz_xyy = pbuffer.data(idx_dip_hf + 113);

    auto tr_x_xyyyz_xzz = pbuffer.data(idx_dip_hf + 115);

    auto tr_x_xyyzz_xxx = pbuffer.data(idx_dip_hf + 120);

    auto tr_x_xyyzz_xxy = pbuffer.data(idx_dip_hf + 121);

    auto tr_x_xyyzz_xxz = pbuffer.data(idx_dip_hf + 122);

    auto tr_x_xyyzz_xyy = pbuffer.data(idx_dip_hf + 123);

    auto tr_x_xyyzz_xzz = pbuffer.data(idx_dip_hf + 125);

    auto tr_x_xyyzz_yyz = pbuffer.data(idx_dip_hf + 127);

    auto tr_x_xyyzz_yzz = pbuffer.data(idx_dip_hf + 128);

    auto tr_x_xyzzz_xxx = pbuffer.data(idx_dip_hf + 130);

    auto tr_x_xyzzz_xxz = pbuffer.data(idx_dip_hf + 132);

    auto tr_x_xyzzz_xzz = pbuffer.data(idx_dip_hf + 135);

    auto tr_x_xzzzz_xxx = pbuffer.data(idx_dip_hf + 140);

    auto tr_x_xzzzz_xxy = pbuffer.data(idx_dip_hf + 141);

    auto tr_x_xzzzz_xxz = pbuffer.data(idx_dip_hf + 142);

    auto tr_x_xzzzz_xyy = pbuffer.data(idx_dip_hf + 143);

    auto tr_x_xzzzz_xyz = pbuffer.data(idx_dip_hf + 144);

    auto tr_x_xzzzz_xzz = pbuffer.data(idx_dip_hf + 145);

    auto tr_x_xzzzz_yyz = pbuffer.data(idx_dip_hf + 147);

    auto tr_x_xzzzz_yzz = pbuffer.data(idx_dip_hf + 148);

    auto tr_x_xzzzz_zzz = pbuffer.data(idx_dip_hf + 149);

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

    auto tr_x_yyyyz_xxy = pbuffer.data(idx_dip_hf + 161);

    auto tr_x_yyyyz_xxz = pbuffer.data(idx_dip_hf + 162);

    auto tr_x_yyyyz_xyy = pbuffer.data(idx_dip_hf + 163);

    auto tr_x_yyyyz_xzz = pbuffer.data(idx_dip_hf + 165);

    auto tr_x_yyyyz_yyy = pbuffer.data(idx_dip_hf + 166);

    auto tr_x_yyyyz_yyz = pbuffer.data(idx_dip_hf + 167);

    auto tr_x_yyyyz_yzz = pbuffer.data(idx_dip_hf + 168);

    auto tr_x_yyyyz_zzz = pbuffer.data(idx_dip_hf + 169);

    auto tr_x_yyyzz_xxx = pbuffer.data(idx_dip_hf + 170);

    auto tr_x_yyyzz_xxy = pbuffer.data(idx_dip_hf + 171);

    auto tr_x_yyyzz_xxz = pbuffer.data(idx_dip_hf + 172);

    auto tr_x_yyyzz_xyy = pbuffer.data(idx_dip_hf + 173);

    auto tr_x_yyyzz_xyz = pbuffer.data(idx_dip_hf + 174);

    auto tr_x_yyyzz_xzz = pbuffer.data(idx_dip_hf + 175);

    auto tr_x_yyyzz_yyy = pbuffer.data(idx_dip_hf + 176);

    auto tr_x_yyyzz_yyz = pbuffer.data(idx_dip_hf + 177);

    auto tr_x_yyyzz_yzz = pbuffer.data(idx_dip_hf + 178);

    auto tr_x_yyyzz_zzz = pbuffer.data(idx_dip_hf + 179);

    auto tr_x_yyzzz_xxx = pbuffer.data(idx_dip_hf + 180);

    auto tr_x_yyzzz_xxy = pbuffer.data(idx_dip_hf + 181);

    auto tr_x_yyzzz_xxz = pbuffer.data(idx_dip_hf + 182);

    auto tr_x_yyzzz_xyy = pbuffer.data(idx_dip_hf + 183);

    auto tr_x_yyzzz_xyz = pbuffer.data(idx_dip_hf + 184);

    auto tr_x_yyzzz_xzz = pbuffer.data(idx_dip_hf + 185);

    auto tr_x_yyzzz_yyy = pbuffer.data(idx_dip_hf + 186);

    auto tr_x_yyzzz_yyz = pbuffer.data(idx_dip_hf + 187);

    auto tr_x_yyzzz_yzz = pbuffer.data(idx_dip_hf + 188);

    auto tr_x_yyzzz_zzz = pbuffer.data(idx_dip_hf + 189);

    auto tr_x_yzzzz_xxx = pbuffer.data(idx_dip_hf + 190);

    auto tr_x_yzzzz_xxz = pbuffer.data(idx_dip_hf + 192);

    auto tr_x_yzzzz_xyz = pbuffer.data(idx_dip_hf + 194);

    auto tr_x_yzzzz_xzz = pbuffer.data(idx_dip_hf + 195);

    auto tr_x_yzzzz_yyy = pbuffer.data(idx_dip_hf + 196);

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

    auto tr_y_xxxxy_xxx = pbuffer.data(idx_dip_hf + 220);

    auto tr_y_xxxxy_xxy = pbuffer.data(idx_dip_hf + 221);

    auto tr_y_xxxxy_xyy = pbuffer.data(idx_dip_hf + 223);

    auto tr_y_xxxxy_xyz = pbuffer.data(idx_dip_hf + 224);

    auto tr_y_xxxxy_yyy = pbuffer.data(idx_dip_hf + 226);

    auto tr_y_xxxxy_yyz = pbuffer.data(idx_dip_hf + 227);

    auto tr_y_xxxxy_yzz = pbuffer.data(idx_dip_hf + 228);

    auto tr_y_xxxxy_zzz = pbuffer.data(idx_dip_hf + 229);

    auto tr_y_xxxxz_xxx = pbuffer.data(idx_dip_hf + 230);

    auto tr_y_xxxxz_xxy = pbuffer.data(idx_dip_hf + 231);

    auto tr_y_xxxxz_xxz = pbuffer.data(idx_dip_hf + 232);

    auto tr_y_xxxxz_xyy = pbuffer.data(idx_dip_hf + 233);

    auto tr_y_xxxxz_xzz = pbuffer.data(idx_dip_hf + 235);

    auto tr_y_xxxxz_yyz = pbuffer.data(idx_dip_hf + 237);

    auto tr_y_xxxxz_yzz = pbuffer.data(idx_dip_hf + 238);

    auto tr_y_xxxxz_zzz = pbuffer.data(idx_dip_hf + 239);

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

    auto tr_y_xxxyz_xxy = pbuffer.data(idx_dip_hf + 251);

    auto tr_y_xxxyz_xyy = pbuffer.data(idx_dip_hf + 253);

    auto tr_y_xxxyz_yyz = pbuffer.data(idx_dip_hf + 257);

    auto tr_y_xxxyz_yzz = pbuffer.data(idx_dip_hf + 258);

    auto tr_y_xxxyz_zzz = pbuffer.data(idx_dip_hf + 259);

    auto tr_y_xxxzz_xxx = pbuffer.data(idx_dip_hf + 260);

    auto tr_y_xxxzz_xxy = pbuffer.data(idx_dip_hf + 261);

    auto tr_y_xxxzz_xxz = pbuffer.data(idx_dip_hf + 262);

    auto tr_y_xxxzz_xyy = pbuffer.data(idx_dip_hf + 263);

    auto tr_y_xxxzz_xyz = pbuffer.data(idx_dip_hf + 264);

    auto tr_y_xxxzz_xzz = pbuffer.data(idx_dip_hf + 265);

    auto tr_y_xxxzz_yyy = pbuffer.data(idx_dip_hf + 266);

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

    auto tr_y_xxyyz_xxx = pbuffer.data(idx_dip_hf + 280);

    auto tr_y_xxyyz_xxy = pbuffer.data(idx_dip_hf + 281);

    auto tr_y_xxyyz_xyy = pbuffer.data(idx_dip_hf + 283);

    auto tr_y_xxyyz_yyz = pbuffer.data(idx_dip_hf + 287);

    auto tr_y_xxyyz_yzz = pbuffer.data(idx_dip_hf + 288);

    auto tr_y_xxyyz_zzz = pbuffer.data(idx_dip_hf + 289);

    auto tr_y_xxyzz_xxy = pbuffer.data(idx_dip_hf + 291);

    auto tr_y_xxyzz_xyy = pbuffer.data(idx_dip_hf + 293);

    auto tr_y_xxyzz_xyz = pbuffer.data(idx_dip_hf + 294);

    auto tr_y_xxyzz_yyy = pbuffer.data(idx_dip_hf + 296);

    auto tr_y_xxyzz_yyz = pbuffer.data(idx_dip_hf + 297);

    auto tr_y_xxyzz_yzz = pbuffer.data(idx_dip_hf + 298);

    auto tr_y_xxyzz_zzz = pbuffer.data(idx_dip_hf + 299);

    auto tr_y_xxzzz_xxx = pbuffer.data(idx_dip_hf + 300);

    auto tr_y_xxzzz_xxy = pbuffer.data(idx_dip_hf + 301);

    auto tr_y_xxzzz_xxz = pbuffer.data(idx_dip_hf + 302);

    auto tr_y_xxzzz_xyy = pbuffer.data(idx_dip_hf + 303);

    auto tr_y_xxzzz_xyz = pbuffer.data(idx_dip_hf + 304);

    auto tr_y_xxzzz_xzz = pbuffer.data(idx_dip_hf + 305);

    auto tr_y_xxzzz_yyy = pbuffer.data(idx_dip_hf + 306);

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

    auto tr_y_xyyyz_yyz = pbuffer.data(idx_dip_hf + 327);

    auto tr_y_xyyyz_yzz = pbuffer.data(idx_dip_hf + 328);

    auto tr_y_xyyyz_zzz = pbuffer.data(idx_dip_hf + 329);

    auto tr_y_xyyzz_xxz = pbuffer.data(idx_dip_hf + 332);

    auto tr_y_xyyzz_xyz = pbuffer.data(idx_dip_hf + 334);

    auto tr_y_xyyzz_xzz = pbuffer.data(idx_dip_hf + 335);

    auto tr_y_xyyzz_yyy = pbuffer.data(idx_dip_hf + 336);

    auto tr_y_xyyzz_yyz = pbuffer.data(idx_dip_hf + 337);

    auto tr_y_xyyzz_yzz = pbuffer.data(idx_dip_hf + 338);

    auto tr_y_xyyzz_zzz = pbuffer.data(idx_dip_hf + 339);

    auto tr_y_xyzzz_xyz = pbuffer.data(idx_dip_hf + 344);

    auto tr_y_xyzzz_yyy = pbuffer.data(idx_dip_hf + 346);

    auto tr_y_xyzzz_yyz = pbuffer.data(idx_dip_hf + 347);

    auto tr_y_xyzzz_yzz = pbuffer.data(idx_dip_hf + 348);

    auto tr_y_xyzzz_zzz = pbuffer.data(idx_dip_hf + 349);

    auto tr_y_xzzzz_xxz = pbuffer.data(idx_dip_hf + 352);

    auto tr_y_xzzzz_xyz = pbuffer.data(idx_dip_hf + 354);

    auto tr_y_xzzzz_xzz = pbuffer.data(idx_dip_hf + 355);

    auto tr_y_xzzzz_yyy = pbuffer.data(idx_dip_hf + 356);

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

    auto tr_y_yyyyz_xxx = pbuffer.data(idx_dip_hf + 370);

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

    auto tr_z_xxxxy_xxx = pbuffer.data(idx_dip_hf + 430);

    auto tr_z_xxxxy_xxz = pbuffer.data(idx_dip_hf + 432);

    auto tr_z_xxxxy_xzz = pbuffer.data(idx_dip_hf + 435);

    auto tr_z_xxxxy_yyy = pbuffer.data(idx_dip_hf + 436);

    auto tr_z_xxxxy_yyz = pbuffer.data(idx_dip_hf + 437);

    auto tr_z_xxxxy_yzz = pbuffer.data(idx_dip_hf + 438);

    auto tr_z_xxxxz_xxx = pbuffer.data(idx_dip_hf + 440);

    auto tr_z_xxxxz_xxy = pbuffer.data(idx_dip_hf + 441);

    auto tr_z_xxxxz_xxz = pbuffer.data(idx_dip_hf + 442);

    auto tr_z_xxxxz_xyy = pbuffer.data(idx_dip_hf + 443);

    auto tr_z_xxxxz_xyz = pbuffer.data(idx_dip_hf + 444);

    auto tr_z_xxxxz_xzz = pbuffer.data(idx_dip_hf + 445);

    auto tr_z_xxxxz_yyy = pbuffer.data(idx_dip_hf + 446);

    auto tr_z_xxxxz_yyz = pbuffer.data(idx_dip_hf + 447);

    auto tr_z_xxxxz_yzz = pbuffer.data(idx_dip_hf + 448);

    auto tr_z_xxxxz_zzz = pbuffer.data(idx_dip_hf + 449);

    auto tr_z_xxxyy_xxx = pbuffer.data(idx_dip_hf + 450);

    auto tr_z_xxxyy_xxy = pbuffer.data(idx_dip_hf + 451);

    auto tr_z_xxxyy_xxz = pbuffer.data(idx_dip_hf + 452);

    auto tr_z_xxxyy_xyy = pbuffer.data(idx_dip_hf + 453);

    auto tr_z_xxxyy_xyz = pbuffer.data(idx_dip_hf + 454);

    auto tr_z_xxxyy_xzz = pbuffer.data(idx_dip_hf + 455);

    auto tr_z_xxxyy_yyy = pbuffer.data(idx_dip_hf + 456);

    auto tr_z_xxxyy_yyz = pbuffer.data(idx_dip_hf + 457);

    auto tr_z_xxxyy_yzz = pbuffer.data(idx_dip_hf + 458);

    auto tr_z_xxxyy_zzz = pbuffer.data(idx_dip_hf + 459);

    auto tr_z_xxxyz_xxx = pbuffer.data(idx_dip_hf + 460);

    auto tr_z_xxxyz_xxz = pbuffer.data(idx_dip_hf + 462);

    auto tr_z_xxxyz_xzz = pbuffer.data(idx_dip_hf + 465);

    auto tr_z_xxxyz_yyy = pbuffer.data(idx_dip_hf + 466);

    auto tr_z_xxxyz_yyz = pbuffer.data(idx_dip_hf + 467);

    auto tr_z_xxxyz_yzz = pbuffer.data(idx_dip_hf + 468);

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

    auto tr_z_xxyyy_xxx = pbuffer.data(idx_dip_hf + 480);

    auto tr_z_xxyyy_xxy = pbuffer.data(idx_dip_hf + 481);

    auto tr_z_xxyyy_xxz = pbuffer.data(idx_dip_hf + 482);

    auto tr_z_xxyyy_xyy = pbuffer.data(idx_dip_hf + 483);

    auto tr_z_xxyyy_xyz = pbuffer.data(idx_dip_hf + 484);

    auto tr_z_xxyyy_xzz = pbuffer.data(idx_dip_hf + 485);

    auto tr_z_xxyyy_yyy = pbuffer.data(idx_dip_hf + 486);

    auto tr_z_xxyyy_yyz = pbuffer.data(idx_dip_hf + 487);

    auto tr_z_xxyyy_yzz = pbuffer.data(idx_dip_hf + 488);

    auto tr_z_xxyyy_zzz = pbuffer.data(idx_dip_hf + 489);

    auto tr_z_xxyyz_xxx = pbuffer.data(idx_dip_hf + 490);

    auto tr_z_xxyyz_xxz = pbuffer.data(idx_dip_hf + 492);

    auto tr_z_xxyyz_xyz = pbuffer.data(idx_dip_hf + 494);

    auto tr_z_xxyyz_xzz = pbuffer.data(idx_dip_hf + 495);

    auto tr_z_xxyyz_yyy = pbuffer.data(idx_dip_hf + 496);

    auto tr_z_xxyyz_yyz = pbuffer.data(idx_dip_hf + 497);

    auto tr_z_xxyyz_yzz = pbuffer.data(idx_dip_hf + 498);

    auto tr_z_xxyyz_zzz = pbuffer.data(idx_dip_hf + 499);

    auto tr_z_xxyzz_xxx = pbuffer.data(idx_dip_hf + 500);

    auto tr_z_xxyzz_xxz = pbuffer.data(idx_dip_hf + 502);

    auto tr_z_xxyzz_xzz = pbuffer.data(idx_dip_hf + 505);

    auto tr_z_xxyzz_yyy = pbuffer.data(idx_dip_hf + 506);

    auto tr_z_xxyzz_yyz = pbuffer.data(idx_dip_hf + 507);

    auto tr_z_xxyzz_yzz = pbuffer.data(idx_dip_hf + 508);

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

    auto tr_z_xyyyy_zzz = pbuffer.data(idx_dip_hf + 529);

    auto tr_z_xyyyz_xyz = pbuffer.data(idx_dip_hf + 534);

    auto tr_z_xyyyz_yyy = pbuffer.data(idx_dip_hf + 536);

    auto tr_z_xyyyz_yyz = pbuffer.data(idx_dip_hf + 537);

    auto tr_z_xyyyz_yzz = pbuffer.data(idx_dip_hf + 538);

    auto tr_z_xyyyz_zzz = pbuffer.data(idx_dip_hf + 539);

    auto tr_z_xyyzz_xxy = pbuffer.data(idx_dip_hf + 541);

    auto tr_z_xyyzz_xyy = pbuffer.data(idx_dip_hf + 543);

    auto tr_z_xyyzz_xyz = pbuffer.data(idx_dip_hf + 544);

    auto tr_z_xyyzz_yyy = pbuffer.data(idx_dip_hf + 546);

    auto tr_z_xyyzz_yyz = pbuffer.data(idx_dip_hf + 547);

    auto tr_z_xyyzz_yzz = pbuffer.data(idx_dip_hf + 548);

    auto tr_z_xyyzz_zzz = pbuffer.data(idx_dip_hf + 549);

    auto tr_z_xyzzz_yyy = pbuffer.data(idx_dip_hf + 556);

    auto tr_z_xyzzz_yyz = pbuffer.data(idx_dip_hf + 557);

    auto tr_z_xyzzz_yzz = pbuffer.data(idx_dip_hf + 558);

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

    // Set up 0-10 components of targeted buffer : IF

    auto tr_x_xxxxxx_xxx = pbuffer.data(idx_dip_if);

    auto tr_x_xxxxxx_xxy = pbuffer.data(idx_dip_if + 1);

    auto tr_x_xxxxxx_xxz = pbuffer.data(idx_dip_if + 2);

    auto tr_x_xxxxxx_xyy = pbuffer.data(idx_dip_if + 3);

    auto tr_x_xxxxxx_xyz = pbuffer.data(idx_dip_if + 4);

    auto tr_x_xxxxxx_xzz = pbuffer.data(idx_dip_if + 5);

    auto tr_x_xxxxxx_yyy = pbuffer.data(idx_dip_if + 6);

    auto tr_x_xxxxxx_yyz = pbuffer.data(idx_dip_if + 7);

    auto tr_x_xxxxxx_yzz = pbuffer.data(idx_dip_if + 8);

    auto tr_x_xxxxxx_zzz = pbuffer.data(idx_dip_if + 9);

    #pragma omp simd aligned(pa_x, tr_x_xxxx_xxx, tr_x_xxxx_xxy, tr_x_xxxx_xxz, tr_x_xxxx_xyy, tr_x_xxxx_xyz, tr_x_xxxx_xzz, tr_x_xxxx_yyy, tr_x_xxxx_yyz, tr_x_xxxx_yzz, tr_x_xxxx_zzz, tr_x_xxxxx_xx, tr_x_xxxxx_xxx, tr_x_xxxxx_xxy, tr_x_xxxxx_xxz, tr_x_xxxxx_xy, tr_x_xxxxx_xyy, tr_x_xxxxx_xyz, tr_x_xxxxx_xz, tr_x_xxxxx_xzz, tr_x_xxxxx_yy, tr_x_xxxxx_yyy, tr_x_xxxxx_yyz, tr_x_xxxxx_yz, tr_x_xxxxx_yzz, tr_x_xxxxx_zz, tr_x_xxxxx_zzz, tr_x_xxxxxx_xxx, tr_x_xxxxxx_xxy, tr_x_xxxxxx_xxz, tr_x_xxxxxx_xyy, tr_x_xxxxxx_xyz, tr_x_xxxxxx_xzz, tr_x_xxxxxx_yyy, tr_x_xxxxxx_yyz, tr_x_xxxxxx_yzz, tr_x_xxxxxx_zzz, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xzz, ts_xxxxx_yyy, ts_xxxxx_yyz, ts_xxxxx_yzz, ts_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_xxx[i] = 5.0 * tr_x_xxxx_xxx[i] * fe_0 + 3.0 * tr_x_xxxxx_xx[i] * fe_0 + ts_xxxxx_xxx[i] * fe_0 + tr_x_xxxxx_xxx[i] * pa_x[i];

        tr_x_xxxxxx_xxy[i] = 5.0 * tr_x_xxxx_xxy[i] * fe_0 + 2.0 * tr_x_xxxxx_xy[i] * fe_0 + ts_xxxxx_xxy[i] * fe_0 + tr_x_xxxxx_xxy[i] * pa_x[i];

        tr_x_xxxxxx_xxz[i] = 5.0 * tr_x_xxxx_xxz[i] * fe_0 + 2.0 * tr_x_xxxxx_xz[i] * fe_0 + ts_xxxxx_xxz[i] * fe_0 + tr_x_xxxxx_xxz[i] * pa_x[i];

        tr_x_xxxxxx_xyy[i] = 5.0 * tr_x_xxxx_xyy[i] * fe_0 + tr_x_xxxxx_yy[i] * fe_0 + ts_xxxxx_xyy[i] * fe_0 + tr_x_xxxxx_xyy[i] * pa_x[i];

        tr_x_xxxxxx_xyz[i] = 5.0 * tr_x_xxxx_xyz[i] * fe_0 + tr_x_xxxxx_yz[i] * fe_0 + ts_xxxxx_xyz[i] * fe_0 + tr_x_xxxxx_xyz[i] * pa_x[i];

        tr_x_xxxxxx_xzz[i] = 5.0 * tr_x_xxxx_xzz[i] * fe_0 + tr_x_xxxxx_zz[i] * fe_0 + ts_xxxxx_xzz[i] * fe_0 + tr_x_xxxxx_xzz[i] * pa_x[i];

        tr_x_xxxxxx_yyy[i] = 5.0 * tr_x_xxxx_yyy[i] * fe_0 + ts_xxxxx_yyy[i] * fe_0 + tr_x_xxxxx_yyy[i] * pa_x[i];

        tr_x_xxxxxx_yyz[i] = 5.0 * tr_x_xxxx_yyz[i] * fe_0 + ts_xxxxx_yyz[i] * fe_0 + tr_x_xxxxx_yyz[i] * pa_x[i];

        tr_x_xxxxxx_yzz[i] = 5.0 * tr_x_xxxx_yzz[i] * fe_0 + ts_xxxxx_yzz[i] * fe_0 + tr_x_xxxxx_yzz[i] * pa_x[i];

        tr_x_xxxxxx_zzz[i] = 5.0 * tr_x_xxxx_zzz[i] * fe_0 + ts_xxxxx_zzz[i] * fe_0 + tr_x_xxxxx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : IF

    auto tr_x_xxxxxy_xxx = pbuffer.data(idx_dip_if + 10);

    auto tr_x_xxxxxy_xxy = pbuffer.data(idx_dip_if + 11);

    auto tr_x_xxxxxy_xxz = pbuffer.data(idx_dip_if + 12);

    auto tr_x_xxxxxy_xyy = pbuffer.data(idx_dip_if + 13);

    auto tr_x_xxxxxy_xyz = pbuffer.data(idx_dip_if + 14);

    auto tr_x_xxxxxy_xzz = pbuffer.data(idx_dip_if + 15);

    auto tr_x_xxxxxy_yyy = pbuffer.data(idx_dip_if + 16);

    auto tr_x_xxxxxy_yyz = pbuffer.data(idx_dip_if + 17);

    auto tr_x_xxxxxy_yzz = pbuffer.data(idx_dip_if + 18);

    auto tr_x_xxxxxy_zzz = pbuffer.data(idx_dip_if + 19);

    #pragma omp simd aligned(pa_y, tr_x_xxxxx_xx, tr_x_xxxxx_xxx, tr_x_xxxxx_xxy, tr_x_xxxxx_xxz, tr_x_xxxxx_xy, tr_x_xxxxx_xyy, tr_x_xxxxx_xyz, tr_x_xxxxx_xz, tr_x_xxxxx_xzz, tr_x_xxxxx_yy, tr_x_xxxxx_yyy, tr_x_xxxxx_yyz, tr_x_xxxxx_yz, tr_x_xxxxx_yzz, tr_x_xxxxx_zz, tr_x_xxxxx_zzz, tr_x_xxxxxy_xxx, tr_x_xxxxxy_xxy, tr_x_xxxxxy_xxz, tr_x_xxxxxy_xyy, tr_x_xxxxxy_xyz, tr_x_xxxxxy_xzz, tr_x_xxxxxy_yyy, tr_x_xxxxxy_yyz, tr_x_xxxxxy_yzz, tr_x_xxxxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxy_xxx[i] = tr_x_xxxxx_xxx[i] * pa_y[i];

        tr_x_xxxxxy_xxy[i] = tr_x_xxxxx_xx[i] * fe_0 + tr_x_xxxxx_xxy[i] * pa_y[i];

        tr_x_xxxxxy_xxz[i] = tr_x_xxxxx_xxz[i] * pa_y[i];

        tr_x_xxxxxy_xyy[i] = 2.0 * tr_x_xxxxx_xy[i] * fe_0 + tr_x_xxxxx_xyy[i] * pa_y[i];

        tr_x_xxxxxy_xyz[i] = tr_x_xxxxx_xz[i] * fe_0 + tr_x_xxxxx_xyz[i] * pa_y[i];

        tr_x_xxxxxy_xzz[i] = tr_x_xxxxx_xzz[i] * pa_y[i];

        tr_x_xxxxxy_yyy[i] = 3.0 * tr_x_xxxxx_yy[i] * fe_0 + tr_x_xxxxx_yyy[i] * pa_y[i];

        tr_x_xxxxxy_yyz[i] = 2.0 * tr_x_xxxxx_yz[i] * fe_0 + tr_x_xxxxx_yyz[i] * pa_y[i];

        tr_x_xxxxxy_yzz[i] = tr_x_xxxxx_zz[i] * fe_0 + tr_x_xxxxx_yzz[i] * pa_y[i];

        tr_x_xxxxxy_zzz[i] = tr_x_xxxxx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : IF

    auto tr_x_xxxxxz_xxx = pbuffer.data(idx_dip_if + 20);

    auto tr_x_xxxxxz_xxy = pbuffer.data(idx_dip_if + 21);

    auto tr_x_xxxxxz_xxz = pbuffer.data(idx_dip_if + 22);

    auto tr_x_xxxxxz_xyy = pbuffer.data(idx_dip_if + 23);

    auto tr_x_xxxxxz_xyz = pbuffer.data(idx_dip_if + 24);

    auto tr_x_xxxxxz_xzz = pbuffer.data(idx_dip_if + 25);

    auto tr_x_xxxxxz_yyy = pbuffer.data(idx_dip_if + 26);

    auto tr_x_xxxxxz_yyz = pbuffer.data(idx_dip_if + 27);

    auto tr_x_xxxxxz_yzz = pbuffer.data(idx_dip_if + 28);

    auto tr_x_xxxxxz_zzz = pbuffer.data(idx_dip_if + 29);

    #pragma omp simd aligned(pa_z, tr_x_xxxxx_xx, tr_x_xxxxx_xxx, tr_x_xxxxx_xxy, tr_x_xxxxx_xxz, tr_x_xxxxx_xy, tr_x_xxxxx_xyy, tr_x_xxxxx_xyz, tr_x_xxxxx_xz, tr_x_xxxxx_xzz, tr_x_xxxxx_yy, tr_x_xxxxx_yyy, tr_x_xxxxx_yyz, tr_x_xxxxx_yz, tr_x_xxxxx_yzz, tr_x_xxxxx_zz, tr_x_xxxxx_zzz, tr_x_xxxxxz_xxx, tr_x_xxxxxz_xxy, tr_x_xxxxxz_xxz, tr_x_xxxxxz_xyy, tr_x_xxxxxz_xyz, tr_x_xxxxxz_xzz, tr_x_xxxxxz_yyy, tr_x_xxxxxz_yyz, tr_x_xxxxxz_yzz, tr_x_xxxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxz_xxx[i] = tr_x_xxxxx_xxx[i] * pa_z[i];

        tr_x_xxxxxz_xxy[i] = tr_x_xxxxx_xxy[i] * pa_z[i];

        tr_x_xxxxxz_xxz[i] = tr_x_xxxxx_xx[i] * fe_0 + tr_x_xxxxx_xxz[i] * pa_z[i];

        tr_x_xxxxxz_xyy[i] = tr_x_xxxxx_xyy[i] * pa_z[i];

        tr_x_xxxxxz_xyz[i] = tr_x_xxxxx_xy[i] * fe_0 + tr_x_xxxxx_xyz[i] * pa_z[i];

        tr_x_xxxxxz_xzz[i] = 2.0 * tr_x_xxxxx_xz[i] * fe_0 + tr_x_xxxxx_xzz[i] * pa_z[i];

        tr_x_xxxxxz_yyy[i] = tr_x_xxxxx_yyy[i] * pa_z[i];

        tr_x_xxxxxz_yyz[i] = tr_x_xxxxx_yy[i] * fe_0 + tr_x_xxxxx_yyz[i] * pa_z[i];

        tr_x_xxxxxz_yzz[i] = 2.0 * tr_x_xxxxx_yz[i] * fe_0 + tr_x_xxxxx_yzz[i] * pa_z[i];

        tr_x_xxxxxz_zzz[i] = 3.0 * tr_x_xxxxx_zz[i] * fe_0 + tr_x_xxxxx_zzz[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : IF

    auto tr_x_xxxxyy_xxx = pbuffer.data(idx_dip_if + 30);

    auto tr_x_xxxxyy_xxy = pbuffer.data(idx_dip_if + 31);

    auto tr_x_xxxxyy_xxz = pbuffer.data(idx_dip_if + 32);

    auto tr_x_xxxxyy_xyy = pbuffer.data(idx_dip_if + 33);

    auto tr_x_xxxxyy_xyz = pbuffer.data(idx_dip_if + 34);

    auto tr_x_xxxxyy_xzz = pbuffer.data(idx_dip_if + 35);

    auto tr_x_xxxxyy_yyy = pbuffer.data(idx_dip_if + 36);

    auto tr_x_xxxxyy_yyz = pbuffer.data(idx_dip_if + 37);

    auto tr_x_xxxxyy_yzz = pbuffer.data(idx_dip_if + 38);

    auto tr_x_xxxxyy_zzz = pbuffer.data(idx_dip_if + 39);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxxx_xxx, tr_x_xxxx_xxy, tr_x_xxxx_xxz, tr_x_xxxx_xyy, tr_x_xxxx_xyz, tr_x_xxxx_xzz, tr_x_xxxx_zzz, tr_x_xxxxy_xx, tr_x_xxxxy_xxx, tr_x_xxxxy_xxy, tr_x_xxxxy_xxz, tr_x_xxxxy_xy, tr_x_xxxxy_xyy, tr_x_xxxxy_xyz, tr_x_xxxxy_xz, tr_x_xxxxy_xzz, tr_x_xxxxy_zzz, tr_x_xxxxyy_xxx, tr_x_xxxxyy_xxy, tr_x_xxxxyy_xxz, tr_x_xxxxyy_xyy, tr_x_xxxxyy_xyz, tr_x_xxxxyy_xzz, tr_x_xxxxyy_yyy, tr_x_xxxxyy_yyz, tr_x_xxxxyy_yzz, tr_x_xxxxyy_zzz, tr_x_xxxyy_yyy, tr_x_xxxyy_yyz, tr_x_xxxyy_yzz, tr_x_xxyy_yyy, tr_x_xxyy_yyz, tr_x_xxyy_yzz, ts_xxxyy_yyy, ts_xxxyy_yyz, ts_xxxyy_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyy_xxx[i] = tr_x_xxxx_xxx[i] * fe_0 + tr_x_xxxxy_xxx[i] * pa_y[i];

        tr_x_xxxxyy_xxy[i] = tr_x_xxxx_xxy[i] * fe_0 + tr_x_xxxxy_xx[i] * fe_0 + tr_x_xxxxy_xxy[i] * pa_y[i];

        tr_x_xxxxyy_xxz[i] = tr_x_xxxx_xxz[i] * fe_0 + tr_x_xxxxy_xxz[i] * pa_y[i];

        tr_x_xxxxyy_xyy[i] = tr_x_xxxx_xyy[i] * fe_0 + 2.0 * tr_x_xxxxy_xy[i] * fe_0 + tr_x_xxxxy_xyy[i] * pa_y[i];

        tr_x_xxxxyy_xyz[i] = tr_x_xxxx_xyz[i] * fe_0 + tr_x_xxxxy_xz[i] * fe_0 + tr_x_xxxxy_xyz[i] * pa_y[i];

        tr_x_xxxxyy_xzz[i] = tr_x_xxxx_xzz[i] * fe_0 + tr_x_xxxxy_xzz[i] * pa_y[i];

        tr_x_xxxxyy_yyy[i] = 3.0 * tr_x_xxyy_yyy[i] * fe_0 + ts_xxxyy_yyy[i] * fe_0 + tr_x_xxxyy_yyy[i] * pa_x[i];

        tr_x_xxxxyy_yyz[i] = 3.0 * tr_x_xxyy_yyz[i] * fe_0 + ts_xxxyy_yyz[i] * fe_0 + tr_x_xxxyy_yyz[i] * pa_x[i];

        tr_x_xxxxyy_yzz[i] = 3.0 * tr_x_xxyy_yzz[i] * fe_0 + ts_xxxyy_yzz[i] * fe_0 + tr_x_xxxyy_yzz[i] * pa_x[i];

        tr_x_xxxxyy_zzz[i] = tr_x_xxxx_zzz[i] * fe_0 + tr_x_xxxxy_zzz[i] * pa_y[i];
    }

    // Set up 40-50 components of targeted buffer : IF

    auto tr_x_xxxxyz_xxx = pbuffer.data(idx_dip_if + 40);

    auto tr_x_xxxxyz_xxy = pbuffer.data(idx_dip_if + 41);

    auto tr_x_xxxxyz_xxz = pbuffer.data(idx_dip_if + 42);

    auto tr_x_xxxxyz_xyy = pbuffer.data(idx_dip_if + 43);

    auto tr_x_xxxxyz_xyz = pbuffer.data(idx_dip_if + 44);

    auto tr_x_xxxxyz_xzz = pbuffer.data(idx_dip_if + 45);

    auto tr_x_xxxxyz_yyy = pbuffer.data(idx_dip_if + 46);

    auto tr_x_xxxxyz_yyz = pbuffer.data(idx_dip_if + 47);

    auto tr_x_xxxxyz_yzz = pbuffer.data(idx_dip_if + 48);

    auto tr_x_xxxxyz_zzz = pbuffer.data(idx_dip_if + 49);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxxy_xxy, tr_x_xxxxy_xyy, tr_x_xxxxy_yyy, tr_x_xxxxyz_xxx, tr_x_xxxxyz_xxy, tr_x_xxxxyz_xxz, tr_x_xxxxyz_xyy, tr_x_xxxxyz_xyz, tr_x_xxxxyz_xzz, tr_x_xxxxyz_yyy, tr_x_xxxxyz_yyz, tr_x_xxxxyz_yzz, tr_x_xxxxyz_zzz, tr_x_xxxxz_xxx, tr_x_xxxxz_xxz, tr_x_xxxxz_xyz, tr_x_xxxxz_xz, tr_x_xxxxz_xzz, tr_x_xxxxz_yyz, tr_x_xxxxz_yz, tr_x_xxxxz_yzz, tr_x_xxxxz_zz, tr_x_xxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyz_xxx[i] = tr_x_xxxxz_xxx[i] * pa_y[i];

        tr_x_xxxxyz_xxy[i] = tr_x_xxxxy_xxy[i] * pa_z[i];

        tr_x_xxxxyz_xxz[i] = tr_x_xxxxz_xxz[i] * pa_y[i];

        tr_x_xxxxyz_xyy[i] = tr_x_xxxxy_xyy[i] * pa_z[i];

        tr_x_xxxxyz_xyz[i] = tr_x_xxxxz_xz[i] * fe_0 + tr_x_xxxxz_xyz[i] * pa_y[i];

        tr_x_xxxxyz_xzz[i] = tr_x_xxxxz_xzz[i] * pa_y[i];

        tr_x_xxxxyz_yyy[i] = tr_x_xxxxy_yyy[i] * pa_z[i];

        tr_x_xxxxyz_yyz[i] = 2.0 * tr_x_xxxxz_yz[i] * fe_0 + tr_x_xxxxz_yyz[i] * pa_y[i];

        tr_x_xxxxyz_yzz[i] = tr_x_xxxxz_zz[i] * fe_0 + tr_x_xxxxz_yzz[i] * pa_y[i];

        tr_x_xxxxyz_zzz[i] = tr_x_xxxxz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : IF

    auto tr_x_xxxxzz_xxx = pbuffer.data(idx_dip_if + 50);

    auto tr_x_xxxxzz_xxy = pbuffer.data(idx_dip_if + 51);

    auto tr_x_xxxxzz_xxz = pbuffer.data(idx_dip_if + 52);

    auto tr_x_xxxxzz_xyy = pbuffer.data(idx_dip_if + 53);

    auto tr_x_xxxxzz_xyz = pbuffer.data(idx_dip_if + 54);

    auto tr_x_xxxxzz_xzz = pbuffer.data(idx_dip_if + 55);

    auto tr_x_xxxxzz_yyy = pbuffer.data(idx_dip_if + 56);

    auto tr_x_xxxxzz_yyz = pbuffer.data(idx_dip_if + 57);

    auto tr_x_xxxxzz_yzz = pbuffer.data(idx_dip_if + 58);

    auto tr_x_xxxxzz_zzz = pbuffer.data(idx_dip_if + 59);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxxx_xxx, tr_x_xxxx_xxy, tr_x_xxxx_xxz, tr_x_xxxx_xyy, tr_x_xxxx_xyz, tr_x_xxxx_xzz, tr_x_xxxx_yyy, tr_x_xxxxz_xx, tr_x_xxxxz_xxx, tr_x_xxxxz_xxy, tr_x_xxxxz_xxz, tr_x_xxxxz_xy, tr_x_xxxxz_xyy, tr_x_xxxxz_xyz, tr_x_xxxxz_xz, tr_x_xxxxz_xzz, tr_x_xxxxz_yyy, tr_x_xxxxzz_xxx, tr_x_xxxxzz_xxy, tr_x_xxxxzz_xxz, tr_x_xxxxzz_xyy, tr_x_xxxxzz_xyz, tr_x_xxxxzz_xzz, tr_x_xxxxzz_yyy, tr_x_xxxxzz_yyz, tr_x_xxxxzz_yzz, tr_x_xxxxzz_zzz, tr_x_xxxzz_yyz, tr_x_xxxzz_yzz, tr_x_xxxzz_zzz, tr_x_xxzz_yyz, tr_x_xxzz_yzz, tr_x_xxzz_zzz, ts_xxxzz_yyz, ts_xxxzz_yzz, ts_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxzz_xxx[i] = tr_x_xxxx_xxx[i] * fe_0 + tr_x_xxxxz_xxx[i] * pa_z[i];

        tr_x_xxxxzz_xxy[i] = tr_x_xxxx_xxy[i] * fe_0 + tr_x_xxxxz_xxy[i] * pa_z[i];

        tr_x_xxxxzz_xxz[i] = tr_x_xxxx_xxz[i] * fe_0 + tr_x_xxxxz_xx[i] * fe_0 + tr_x_xxxxz_xxz[i] * pa_z[i];

        tr_x_xxxxzz_xyy[i] = tr_x_xxxx_xyy[i] * fe_0 + tr_x_xxxxz_xyy[i] * pa_z[i];

        tr_x_xxxxzz_xyz[i] = tr_x_xxxx_xyz[i] * fe_0 + tr_x_xxxxz_xy[i] * fe_0 + tr_x_xxxxz_xyz[i] * pa_z[i];

        tr_x_xxxxzz_xzz[i] = tr_x_xxxx_xzz[i] * fe_0 + 2.0 * tr_x_xxxxz_xz[i] * fe_0 + tr_x_xxxxz_xzz[i] * pa_z[i];

        tr_x_xxxxzz_yyy[i] = tr_x_xxxx_yyy[i] * fe_0 + tr_x_xxxxz_yyy[i] * pa_z[i];

        tr_x_xxxxzz_yyz[i] = 3.0 * tr_x_xxzz_yyz[i] * fe_0 + ts_xxxzz_yyz[i] * fe_0 + tr_x_xxxzz_yyz[i] * pa_x[i];

        tr_x_xxxxzz_yzz[i] = 3.0 * tr_x_xxzz_yzz[i] * fe_0 + ts_xxxzz_yzz[i] * fe_0 + tr_x_xxxzz_yzz[i] * pa_x[i];

        tr_x_xxxxzz_zzz[i] = 3.0 * tr_x_xxzz_zzz[i] * fe_0 + ts_xxxzz_zzz[i] * fe_0 + tr_x_xxxzz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : IF

    auto tr_x_xxxyyy_xxx = pbuffer.data(idx_dip_if + 60);

    auto tr_x_xxxyyy_xxy = pbuffer.data(idx_dip_if + 61);

    auto tr_x_xxxyyy_xxz = pbuffer.data(idx_dip_if + 62);

    auto tr_x_xxxyyy_xyy = pbuffer.data(idx_dip_if + 63);

    auto tr_x_xxxyyy_xyz = pbuffer.data(idx_dip_if + 64);

    auto tr_x_xxxyyy_xzz = pbuffer.data(idx_dip_if + 65);

    auto tr_x_xxxyyy_yyy = pbuffer.data(idx_dip_if + 66);

    auto tr_x_xxxyyy_yyz = pbuffer.data(idx_dip_if + 67);

    auto tr_x_xxxyyy_yzz = pbuffer.data(idx_dip_if + 68);

    auto tr_x_xxxyyy_zzz = pbuffer.data(idx_dip_if + 69);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxxy_xxx, tr_x_xxxy_xxy, tr_x_xxxy_xxz, tr_x_xxxy_xyy, tr_x_xxxy_xyz, tr_x_xxxy_xzz, tr_x_xxxy_zzz, tr_x_xxxyy_xx, tr_x_xxxyy_xxx, tr_x_xxxyy_xxy, tr_x_xxxyy_xxz, tr_x_xxxyy_xy, tr_x_xxxyy_xyy, tr_x_xxxyy_xyz, tr_x_xxxyy_xz, tr_x_xxxyy_xzz, tr_x_xxxyy_zzz, tr_x_xxxyyy_xxx, tr_x_xxxyyy_xxy, tr_x_xxxyyy_xxz, tr_x_xxxyyy_xyy, tr_x_xxxyyy_xyz, tr_x_xxxyyy_xzz, tr_x_xxxyyy_yyy, tr_x_xxxyyy_yyz, tr_x_xxxyyy_yzz, tr_x_xxxyyy_zzz, tr_x_xxyyy_yyy, tr_x_xxyyy_yyz, tr_x_xxyyy_yzz, tr_x_xyyy_yyy, tr_x_xyyy_yyz, tr_x_xyyy_yzz, ts_xxyyy_yyy, ts_xxyyy_yyz, ts_xxyyy_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyy_xxx[i] = 2.0 * tr_x_xxxy_xxx[i] * fe_0 + tr_x_xxxyy_xxx[i] * pa_y[i];

        tr_x_xxxyyy_xxy[i] = 2.0 * tr_x_xxxy_xxy[i] * fe_0 + tr_x_xxxyy_xx[i] * fe_0 + tr_x_xxxyy_xxy[i] * pa_y[i];

        tr_x_xxxyyy_xxz[i] = 2.0 * tr_x_xxxy_xxz[i] * fe_0 + tr_x_xxxyy_xxz[i] * pa_y[i];

        tr_x_xxxyyy_xyy[i] = 2.0 * tr_x_xxxy_xyy[i] * fe_0 + 2.0 * tr_x_xxxyy_xy[i] * fe_0 + tr_x_xxxyy_xyy[i] * pa_y[i];

        tr_x_xxxyyy_xyz[i] = 2.0 * tr_x_xxxy_xyz[i] * fe_0 + tr_x_xxxyy_xz[i] * fe_0 + tr_x_xxxyy_xyz[i] * pa_y[i];

        tr_x_xxxyyy_xzz[i] = 2.0 * tr_x_xxxy_xzz[i] * fe_0 + tr_x_xxxyy_xzz[i] * pa_y[i];

        tr_x_xxxyyy_yyy[i] = 2.0 * tr_x_xyyy_yyy[i] * fe_0 + ts_xxyyy_yyy[i] * fe_0 + tr_x_xxyyy_yyy[i] * pa_x[i];

        tr_x_xxxyyy_yyz[i] = 2.0 * tr_x_xyyy_yyz[i] * fe_0 + ts_xxyyy_yyz[i] * fe_0 + tr_x_xxyyy_yyz[i] * pa_x[i];

        tr_x_xxxyyy_yzz[i] = 2.0 * tr_x_xyyy_yzz[i] * fe_0 + ts_xxyyy_yzz[i] * fe_0 + tr_x_xxyyy_yzz[i] * pa_x[i];

        tr_x_xxxyyy_zzz[i] = 2.0 * tr_x_xxxy_zzz[i] * fe_0 + tr_x_xxxyy_zzz[i] * pa_y[i];
    }

    // Set up 70-80 components of targeted buffer : IF

    auto tr_x_xxxyyz_xxx = pbuffer.data(idx_dip_if + 70);

    auto tr_x_xxxyyz_xxy = pbuffer.data(idx_dip_if + 71);

    auto tr_x_xxxyyz_xxz = pbuffer.data(idx_dip_if + 72);

    auto tr_x_xxxyyz_xyy = pbuffer.data(idx_dip_if + 73);

    auto tr_x_xxxyyz_xyz = pbuffer.data(idx_dip_if + 74);

    auto tr_x_xxxyyz_xzz = pbuffer.data(idx_dip_if + 75);

    auto tr_x_xxxyyz_yyy = pbuffer.data(idx_dip_if + 76);

    auto tr_x_xxxyyz_yyz = pbuffer.data(idx_dip_if + 77);

    auto tr_x_xxxyyz_yzz = pbuffer.data(idx_dip_if + 78);

    auto tr_x_xxxyyz_zzz = pbuffer.data(idx_dip_if + 79);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxyy_xxx, tr_x_xxxyy_xxy, tr_x_xxxyy_xy, tr_x_xxxyy_xyy, tr_x_xxxyy_xyz, tr_x_xxxyy_yy, tr_x_xxxyy_yyy, tr_x_xxxyy_yyz, tr_x_xxxyy_yz, tr_x_xxxyy_yzz, tr_x_xxxyyz_xxx, tr_x_xxxyyz_xxy, tr_x_xxxyyz_xxz, tr_x_xxxyyz_xyy, tr_x_xxxyyz_xyz, tr_x_xxxyyz_xzz, tr_x_xxxyyz_yyy, tr_x_xxxyyz_yyz, tr_x_xxxyyz_yzz, tr_x_xxxyyz_zzz, tr_x_xxxyz_xxz, tr_x_xxxyz_xzz, tr_x_xxxyz_zzz, tr_x_xxxz_xxz, tr_x_xxxz_xzz, tr_x_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyz_xxx[i] = tr_x_xxxyy_xxx[i] * pa_z[i];

        tr_x_xxxyyz_xxy[i] = tr_x_xxxyy_xxy[i] * pa_z[i];

        tr_x_xxxyyz_xxz[i] = tr_x_xxxz_xxz[i] * fe_0 + tr_x_xxxyz_xxz[i] * pa_y[i];

        tr_x_xxxyyz_xyy[i] = tr_x_xxxyy_xyy[i] * pa_z[i];

        tr_x_xxxyyz_xyz[i] = tr_x_xxxyy_xy[i] * fe_0 + tr_x_xxxyy_xyz[i] * pa_z[i];

        tr_x_xxxyyz_xzz[i] = tr_x_xxxz_xzz[i] * fe_0 + tr_x_xxxyz_xzz[i] * pa_y[i];

        tr_x_xxxyyz_yyy[i] = tr_x_xxxyy_yyy[i] * pa_z[i];

        tr_x_xxxyyz_yyz[i] = tr_x_xxxyy_yy[i] * fe_0 + tr_x_xxxyy_yyz[i] * pa_z[i];

        tr_x_xxxyyz_yzz[i] = 2.0 * tr_x_xxxyy_yz[i] * fe_0 + tr_x_xxxyy_yzz[i] * pa_z[i];

        tr_x_xxxyyz_zzz[i] = tr_x_xxxz_zzz[i] * fe_0 + tr_x_xxxyz_zzz[i] * pa_y[i];
    }

    // Set up 80-90 components of targeted buffer : IF

    auto tr_x_xxxyzz_xxx = pbuffer.data(idx_dip_if + 80);

    auto tr_x_xxxyzz_xxy = pbuffer.data(idx_dip_if + 81);

    auto tr_x_xxxyzz_xxz = pbuffer.data(idx_dip_if + 82);

    auto tr_x_xxxyzz_xyy = pbuffer.data(idx_dip_if + 83);

    auto tr_x_xxxyzz_xyz = pbuffer.data(idx_dip_if + 84);

    auto tr_x_xxxyzz_xzz = pbuffer.data(idx_dip_if + 85);

    auto tr_x_xxxyzz_yyy = pbuffer.data(idx_dip_if + 86);

    auto tr_x_xxxyzz_yyz = pbuffer.data(idx_dip_if + 87);

    auto tr_x_xxxyzz_yzz = pbuffer.data(idx_dip_if + 88);

    auto tr_x_xxxyzz_zzz = pbuffer.data(idx_dip_if + 89);

    #pragma omp simd aligned(pa_y, tr_x_xxxyzz_xxx, tr_x_xxxyzz_xxy, tr_x_xxxyzz_xxz, tr_x_xxxyzz_xyy, tr_x_xxxyzz_xyz, tr_x_xxxyzz_xzz, tr_x_xxxyzz_yyy, tr_x_xxxyzz_yyz, tr_x_xxxyzz_yzz, tr_x_xxxyzz_zzz, tr_x_xxxzz_xx, tr_x_xxxzz_xxx, tr_x_xxxzz_xxy, tr_x_xxxzz_xxz, tr_x_xxxzz_xy, tr_x_xxxzz_xyy, tr_x_xxxzz_xyz, tr_x_xxxzz_xz, tr_x_xxxzz_xzz, tr_x_xxxzz_yy, tr_x_xxxzz_yyy, tr_x_xxxzz_yyz, tr_x_xxxzz_yz, tr_x_xxxzz_yzz, tr_x_xxxzz_zz, tr_x_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyzz_xxx[i] = tr_x_xxxzz_xxx[i] * pa_y[i];

        tr_x_xxxyzz_xxy[i] = tr_x_xxxzz_xx[i] * fe_0 + tr_x_xxxzz_xxy[i] * pa_y[i];

        tr_x_xxxyzz_xxz[i] = tr_x_xxxzz_xxz[i] * pa_y[i];

        tr_x_xxxyzz_xyy[i] = 2.0 * tr_x_xxxzz_xy[i] * fe_0 + tr_x_xxxzz_xyy[i] * pa_y[i];

        tr_x_xxxyzz_xyz[i] = tr_x_xxxzz_xz[i] * fe_0 + tr_x_xxxzz_xyz[i] * pa_y[i];

        tr_x_xxxyzz_xzz[i] = tr_x_xxxzz_xzz[i] * pa_y[i];

        tr_x_xxxyzz_yyy[i] = 3.0 * tr_x_xxxzz_yy[i] * fe_0 + tr_x_xxxzz_yyy[i] * pa_y[i];

        tr_x_xxxyzz_yyz[i] = 2.0 * tr_x_xxxzz_yz[i] * fe_0 + tr_x_xxxzz_yyz[i] * pa_y[i];

        tr_x_xxxyzz_yzz[i] = tr_x_xxxzz_zz[i] * fe_0 + tr_x_xxxzz_yzz[i] * pa_y[i];

        tr_x_xxxyzz_zzz[i] = tr_x_xxxzz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : IF

    auto tr_x_xxxzzz_xxx = pbuffer.data(idx_dip_if + 90);

    auto tr_x_xxxzzz_xxy = pbuffer.data(idx_dip_if + 91);

    auto tr_x_xxxzzz_xxz = pbuffer.data(idx_dip_if + 92);

    auto tr_x_xxxzzz_xyy = pbuffer.data(idx_dip_if + 93);

    auto tr_x_xxxzzz_xyz = pbuffer.data(idx_dip_if + 94);

    auto tr_x_xxxzzz_xzz = pbuffer.data(idx_dip_if + 95);

    auto tr_x_xxxzzz_yyy = pbuffer.data(idx_dip_if + 96);

    auto tr_x_xxxzzz_yyz = pbuffer.data(idx_dip_if + 97);

    auto tr_x_xxxzzz_yzz = pbuffer.data(idx_dip_if + 98);

    auto tr_x_xxxzzz_zzz = pbuffer.data(idx_dip_if + 99);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxxz_xxx, tr_x_xxxz_xxy, tr_x_xxxz_xxz, tr_x_xxxz_xyy, tr_x_xxxz_xyz, tr_x_xxxz_xzz, tr_x_xxxz_yyy, tr_x_xxxzz_xx, tr_x_xxxzz_xxx, tr_x_xxxzz_xxy, tr_x_xxxzz_xxz, tr_x_xxxzz_xy, tr_x_xxxzz_xyy, tr_x_xxxzz_xyz, tr_x_xxxzz_xz, tr_x_xxxzz_xzz, tr_x_xxxzz_yyy, tr_x_xxxzzz_xxx, tr_x_xxxzzz_xxy, tr_x_xxxzzz_xxz, tr_x_xxxzzz_xyy, tr_x_xxxzzz_xyz, tr_x_xxxzzz_xzz, tr_x_xxxzzz_yyy, tr_x_xxxzzz_yyz, tr_x_xxxzzz_yzz, tr_x_xxxzzz_zzz, tr_x_xxzzz_yyz, tr_x_xxzzz_yzz, tr_x_xxzzz_zzz, tr_x_xzzz_yyz, tr_x_xzzz_yzz, tr_x_xzzz_zzz, ts_xxzzz_yyz, ts_xxzzz_yzz, ts_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzzz_xxx[i] = 2.0 * tr_x_xxxz_xxx[i] * fe_0 + tr_x_xxxzz_xxx[i] * pa_z[i];

        tr_x_xxxzzz_xxy[i] = 2.0 * tr_x_xxxz_xxy[i] * fe_0 + tr_x_xxxzz_xxy[i] * pa_z[i];

        tr_x_xxxzzz_xxz[i] = 2.0 * tr_x_xxxz_xxz[i] * fe_0 + tr_x_xxxzz_xx[i] * fe_0 + tr_x_xxxzz_xxz[i] * pa_z[i];

        tr_x_xxxzzz_xyy[i] = 2.0 * tr_x_xxxz_xyy[i] * fe_0 + tr_x_xxxzz_xyy[i] * pa_z[i];

        tr_x_xxxzzz_xyz[i] = 2.0 * tr_x_xxxz_xyz[i] * fe_0 + tr_x_xxxzz_xy[i] * fe_0 + tr_x_xxxzz_xyz[i] * pa_z[i];

        tr_x_xxxzzz_xzz[i] = 2.0 * tr_x_xxxz_xzz[i] * fe_0 + 2.0 * tr_x_xxxzz_xz[i] * fe_0 + tr_x_xxxzz_xzz[i] * pa_z[i];

        tr_x_xxxzzz_yyy[i] = 2.0 * tr_x_xxxz_yyy[i] * fe_0 + tr_x_xxxzz_yyy[i] * pa_z[i];

        tr_x_xxxzzz_yyz[i] = 2.0 * tr_x_xzzz_yyz[i] * fe_0 + ts_xxzzz_yyz[i] * fe_0 + tr_x_xxzzz_yyz[i] * pa_x[i];

        tr_x_xxxzzz_yzz[i] = 2.0 * tr_x_xzzz_yzz[i] * fe_0 + ts_xxzzz_yzz[i] * fe_0 + tr_x_xxzzz_yzz[i] * pa_x[i];

        tr_x_xxxzzz_zzz[i] = 2.0 * tr_x_xzzz_zzz[i] * fe_0 + ts_xxzzz_zzz[i] * fe_0 + tr_x_xxzzz_zzz[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : IF

    auto tr_x_xxyyyy_xxx = pbuffer.data(idx_dip_if + 100);

    auto tr_x_xxyyyy_xxy = pbuffer.data(idx_dip_if + 101);

    auto tr_x_xxyyyy_xxz = pbuffer.data(idx_dip_if + 102);

    auto tr_x_xxyyyy_xyy = pbuffer.data(idx_dip_if + 103);

    auto tr_x_xxyyyy_xyz = pbuffer.data(idx_dip_if + 104);

    auto tr_x_xxyyyy_xzz = pbuffer.data(idx_dip_if + 105);

    auto tr_x_xxyyyy_yyy = pbuffer.data(idx_dip_if + 106);

    auto tr_x_xxyyyy_yyz = pbuffer.data(idx_dip_if + 107);

    auto tr_x_xxyyyy_yzz = pbuffer.data(idx_dip_if + 108);

    auto tr_x_xxyyyy_zzz = pbuffer.data(idx_dip_if + 109);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxyy_xxx, tr_x_xxyy_xxy, tr_x_xxyy_xxz, tr_x_xxyy_xyy, tr_x_xxyy_xyz, tr_x_xxyy_xzz, tr_x_xxyy_zzz, tr_x_xxyyy_xx, tr_x_xxyyy_xxx, tr_x_xxyyy_xxy, tr_x_xxyyy_xxz, tr_x_xxyyy_xy, tr_x_xxyyy_xyy, tr_x_xxyyy_xyz, tr_x_xxyyy_xz, tr_x_xxyyy_xzz, tr_x_xxyyy_zzz, tr_x_xxyyyy_xxx, tr_x_xxyyyy_xxy, tr_x_xxyyyy_xxz, tr_x_xxyyyy_xyy, tr_x_xxyyyy_xyz, tr_x_xxyyyy_xzz, tr_x_xxyyyy_yyy, tr_x_xxyyyy_yyz, tr_x_xxyyyy_yzz, tr_x_xxyyyy_zzz, tr_x_xyyyy_yyy, tr_x_xyyyy_yyz, tr_x_xyyyy_yzz, tr_x_yyyy_yyy, tr_x_yyyy_yyz, tr_x_yyyy_yzz, ts_xyyyy_yyy, ts_xyyyy_yyz, ts_xyyyy_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyy_xxx[i] = 3.0 * tr_x_xxyy_xxx[i] * fe_0 + tr_x_xxyyy_xxx[i] * pa_y[i];

        tr_x_xxyyyy_xxy[i] = 3.0 * tr_x_xxyy_xxy[i] * fe_0 + tr_x_xxyyy_xx[i] * fe_0 + tr_x_xxyyy_xxy[i] * pa_y[i];

        tr_x_xxyyyy_xxz[i] = 3.0 * tr_x_xxyy_xxz[i] * fe_0 + tr_x_xxyyy_xxz[i] * pa_y[i];

        tr_x_xxyyyy_xyy[i] = 3.0 * tr_x_xxyy_xyy[i] * fe_0 + 2.0 * tr_x_xxyyy_xy[i] * fe_0 + tr_x_xxyyy_xyy[i] * pa_y[i];

        tr_x_xxyyyy_xyz[i] = 3.0 * tr_x_xxyy_xyz[i] * fe_0 + tr_x_xxyyy_xz[i] * fe_0 + tr_x_xxyyy_xyz[i] * pa_y[i];

        tr_x_xxyyyy_xzz[i] = 3.0 * tr_x_xxyy_xzz[i] * fe_0 + tr_x_xxyyy_xzz[i] * pa_y[i];

        tr_x_xxyyyy_yyy[i] = tr_x_yyyy_yyy[i] * fe_0 + ts_xyyyy_yyy[i] * fe_0 + tr_x_xyyyy_yyy[i] * pa_x[i];

        tr_x_xxyyyy_yyz[i] = tr_x_yyyy_yyz[i] * fe_0 + ts_xyyyy_yyz[i] * fe_0 + tr_x_xyyyy_yyz[i] * pa_x[i];

        tr_x_xxyyyy_yzz[i] = tr_x_yyyy_yzz[i] * fe_0 + ts_xyyyy_yzz[i] * fe_0 + tr_x_xyyyy_yzz[i] * pa_x[i];

        tr_x_xxyyyy_zzz[i] = 3.0 * tr_x_xxyy_zzz[i] * fe_0 + tr_x_xxyyy_zzz[i] * pa_y[i];
    }

    // Set up 110-120 components of targeted buffer : IF

    auto tr_x_xxyyyz_xxx = pbuffer.data(idx_dip_if + 110);

    auto tr_x_xxyyyz_xxy = pbuffer.data(idx_dip_if + 111);

    auto tr_x_xxyyyz_xxz = pbuffer.data(idx_dip_if + 112);

    auto tr_x_xxyyyz_xyy = pbuffer.data(idx_dip_if + 113);

    auto tr_x_xxyyyz_xyz = pbuffer.data(idx_dip_if + 114);

    auto tr_x_xxyyyz_xzz = pbuffer.data(idx_dip_if + 115);

    auto tr_x_xxyyyz_yyy = pbuffer.data(idx_dip_if + 116);

    auto tr_x_xxyyyz_yyz = pbuffer.data(idx_dip_if + 117);

    auto tr_x_xxyyyz_yzz = pbuffer.data(idx_dip_if + 118);

    auto tr_x_xxyyyz_zzz = pbuffer.data(idx_dip_if + 119);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxyyy_xxx, tr_x_xxyyy_xxy, tr_x_xxyyy_xy, tr_x_xxyyy_xyy, tr_x_xxyyy_xyz, tr_x_xxyyy_yy, tr_x_xxyyy_yyy, tr_x_xxyyy_yyz, tr_x_xxyyy_yz, tr_x_xxyyy_yzz, tr_x_xxyyyz_xxx, tr_x_xxyyyz_xxy, tr_x_xxyyyz_xxz, tr_x_xxyyyz_xyy, tr_x_xxyyyz_xyz, tr_x_xxyyyz_xzz, tr_x_xxyyyz_yyy, tr_x_xxyyyz_yyz, tr_x_xxyyyz_yzz, tr_x_xxyyyz_zzz, tr_x_xxyyz_xxz, tr_x_xxyyz_xzz, tr_x_xxyyz_zzz, tr_x_xxyz_xxz, tr_x_xxyz_xzz, tr_x_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyz_xxx[i] = tr_x_xxyyy_xxx[i] * pa_z[i];

        tr_x_xxyyyz_xxy[i] = tr_x_xxyyy_xxy[i] * pa_z[i];

        tr_x_xxyyyz_xxz[i] = 2.0 * tr_x_xxyz_xxz[i] * fe_0 + tr_x_xxyyz_xxz[i] * pa_y[i];

        tr_x_xxyyyz_xyy[i] = tr_x_xxyyy_xyy[i] * pa_z[i];

        tr_x_xxyyyz_xyz[i] = tr_x_xxyyy_xy[i] * fe_0 + tr_x_xxyyy_xyz[i] * pa_z[i];

        tr_x_xxyyyz_xzz[i] = 2.0 * tr_x_xxyz_xzz[i] * fe_0 + tr_x_xxyyz_xzz[i] * pa_y[i];

        tr_x_xxyyyz_yyy[i] = tr_x_xxyyy_yyy[i] * pa_z[i];

        tr_x_xxyyyz_yyz[i] = tr_x_xxyyy_yy[i] * fe_0 + tr_x_xxyyy_yyz[i] * pa_z[i];

        tr_x_xxyyyz_yzz[i] = 2.0 * tr_x_xxyyy_yz[i] * fe_0 + tr_x_xxyyy_yzz[i] * pa_z[i];

        tr_x_xxyyyz_zzz[i] = 2.0 * tr_x_xxyz_zzz[i] * fe_0 + tr_x_xxyyz_zzz[i] * pa_y[i];
    }

    // Set up 120-130 components of targeted buffer : IF

    auto tr_x_xxyyzz_xxx = pbuffer.data(idx_dip_if + 120);

    auto tr_x_xxyyzz_xxy = pbuffer.data(idx_dip_if + 121);

    auto tr_x_xxyyzz_xxz = pbuffer.data(idx_dip_if + 122);

    auto tr_x_xxyyzz_xyy = pbuffer.data(idx_dip_if + 123);

    auto tr_x_xxyyzz_xyz = pbuffer.data(idx_dip_if + 124);

    auto tr_x_xxyyzz_xzz = pbuffer.data(idx_dip_if + 125);

    auto tr_x_xxyyzz_yyy = pbuffer.data(idx_dip_if + 126);

    auto tr_x_xxyyzz_yyz = pbuffer.data(idx_dip_if + 127);

    auto tr_x_xxyyzz_yzz = pbuffer.data(idx_dip_if + 128);

    auto tr_x_xxyyzz_zzz = pbuffer.data(idx_dip_if + 129);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xxyy_xxy, tr_x_xxyy_xyy, tr_x_xxyy_yyy, tr_x_xxyyz_xxy, tr_x_xxyyz_xyy, tr_x_xxyyz_yyy, tr_x_xxyyzz_xxx, tr_x_xxyyzz_xxy, tr_x_xxyyzz_xxz, tr_x_xxyyzz_xyy, tr_x_xxyyzz_xyz, tr_x_xxyyzz_xzz, tr_x_xxyyzz_yyy, tr_x_xxyyzz_yyz, tr_x_xxyyzz_yzz, tr_x_xxyyzz_zzz, tr_x_xxyzz_xxx, tr_x_xxyzz_xxz, tr_x_xxyzz_xyz, tr_x_xxyzz_xz, tr_x_xxyzz_xzz, tr_x_xxyzz_zzz, tr_x_xxzz_xxx, tr_x_xxzz_xxz, tr_x_xxzz_xyz, tr_x_xxzz_xzz, tr_x_xxzz_zzz, tr_x_xyyzz_yyz, tr_x_xyyzz_yzz, tr_x_yyzz_yyz, tr_x_yyzz_yzz, ts_xyyzz_yyz, ts_xyyzz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyzz_xxx[i] = tr_x_xxzz_xxx[i] * fe_0 + tr_x_xxyzz_xxx[i] * pa_y[i];

        tr_x_xxyyzz_xxy[i] = tr_x_xxyy_xxy[i] * fe_0 + tr_x_xxyyz_xxy[i] * pa_z[i];

        tr_x_xxyyzz_xxz[i] = tr_x_xxzz_xxz[i] * fe_0 + tr_x_xxyzz_xxz[i] * pa_y[i];

        tr_x_xxyyzz_xyy[i] = tr_x_xxyy_xyy[i] * fe_0 + tr_x_xxyyz_xyy[i] * pa_z[i];

        tr_x_xxyyzz_xyz[i] = tr_x_xxzz_xyz[i] * fe_0 + tr_x_xxyzz_xz[i] * fe_0 + tr_x_xxyzz_xyz[i] * pa_y[i];

        tr_x_xxyyzz_xzz[i] = tr_x_xxzz_xzz[i] * fe_0 + tr_x_xxyzz_xzz[i] * pa_y[i];

        tr_x_xxyyzz_yyy[i] = tr_x_xxyy_yyy[i] * fe_0 + tr_x_xxyyz_yyy[i] * pa_z[i];

        tr_x_xxyyzz_yyz[i] = tr_x_yyzz_yyz[i] * fe_0 + ts_xyyzz_yyz[i] * fe_0 + tr_x_xyyzz_yyz[i] * pa_x[i];

        tr_x_xxyyzz_yzz[i] = tr_x_yyzz_yzz[i] * fe_0 + ts_xyyzz_yzz[i] * fe_0 + tr_x_xyyzz_yzz[i] * pa_x[i];

        tr_x_xxyyzz_zzz[i] = tr_x_xxzz_zzz[i] * fe_0 + tr_x_xxyzz_zzz[i] * pa_y[i];
    }

    // Set up 130-140 components of targeted buffer : IF

    auto tr_x_xxyzzz_xxx = pbuffer.data(idx_dip_if + 130);

    auto tr_x_xxyzzz_xxy = pbuffer.data(idx_dip_if + 131);

    auto tr_x_xxyzzz_xxz = pbuffer.data(idx_dip_if + 132);

    auto tr_x_xxyzzz_xyy = pbuffer.data(idx_dip_if + 133);

    auto tr_x_xxyzzz_xyz = pbuffer.data(idx_dip_if + 134);

    auto tr_x_xxyzzz_xzz = pbuffer.data(idx_dip_if + 135);

    auto tr_x_xxyzzz_yyy = pbuffer.data(idx_dip_if + 136);

    auto tr_x_xxyzzz_yyz = pbuffer.data(idx_dip_if + 137);

    auto tr_x_xxyzzz_yzz = pbuffer.data(idx_dip_if + 138);

    auto tr_x_xxyzzz_zzz = pbuffer.data(idx_dip_if + 139);

    #pragma omp simd aligned(pa_y, tr_x_xxyzzz_xxx, tr_x_xxyzzz_xxy, tr_x_xxyzzz_xxz, tr_x_xxyzzz_xyy, tr_x_xxyzzz_xyz, tr_x_xxyzzz_xzz, tr_x_xxyzzz_yyy, tr_x_xxyzzz_yyz, tr_x_xxyzzz_yzz, tr_x_xxyzzz_zzz, tr_x_xxzzz_xx, tr_x_xxzzz_xxx, tr_x_xxzzz_xxy, tr_x_xxzzz_xxz, tr_x_xxzzz_xy, tr_x_xxzzz_xyy, tr_x_xxzzz_xyz, tr_x_xxzzz_xz, tr_x_xxzzz_xzz, tr_x_xxzzz_yy, tr_x_xxzzz_yyy, tr_x_xxzzz_yyz, tr_x_xxzzz_yz, tr_x_xxzzz_yzz, tr_x_xxzzz_zz, tr_x_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzzz_xxx[i] = tr_x_xxzzz_xxx[i] * pa_y[i];

        tr_x_xxyzzz_xxy[i] = tr_x_xxzzz_xx[i] * fe_0 + tr_x_xxzzz_xxy[i] * pa_y[i];

        tr_x_xxyzzz_xxz[i] = tr_x_xxzzz_xxz[i] * pa_y[i];

        tr_x_xxyzzz_xyy[i] = 2.0 * tr_x_xxzzz_xy[i] * fe_0 + tr_x_xxzzz_xyy[i] * pa_y[i];

        tr_x_xxyzzz_xyz[i] = tr_x_xxzzz_xz[i] * fe_0 + tr_x_xxzzz_xyz[i] * pa_y[i];

        tr_x_xxyzzz_xzz[i] = tr_x_xxzzz_xzz[i] * pa_y[i];

        tr_x_xxyzzz_yyy[i] = 3.0 * tr_x_xxzzz_yy[i] * fe_0 + tr_x_xxzzz_yyy[i] * pa_y[i];

        tr_x_xxyzzz_yyz[i] = 2.0 * tr_x_xxzzz_yz[i] * fe_0 + tr_x_xxzzz_yyz[i] * pa_y[i];

        tr_x_xxyzzz_yzz[i] = tr_x_xxzzz_zz[i] * fe_0 + tr_x_xxzzz_yzz[i] * pa_y[i];

        tr_x_xxyzzz_zzz[i] = tr_x_xxzzz_zzz[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : IF

    auto tr_x_xxzzzz_xxx = pbuffer.data(idx_dip_if + 140);

    auto tr_x_xxzzzz_xxy = pbuffer.data(idx_dip_if + 141);

    auto tr_x_xxzzzz_xxz = pbuffer.data(idx_dip_if + 142);

    auto tr_x_xxzzzz_xyy = pbuffer.data(idx_dip_if + 143);

    auto tr_x_xxzzzz_xyz = pbuffer.data(idx_dip_if + 144);

    auto tr_x_xxzzzz_xzz = pbuffer.data(idx_dip_if + 145);

    auto tr_x_xxzzzz_yyy = pbuffer.data(idx_dip_if + 146);

    auto tr_x_xxzzzz_yyz = pbuffer.data(idx_dip_if + 147);

    auto tr_x_xxzzzz_yzz = pbuffer.data(idx_dip_if + 148);

    auto tr_x_xxzzzz_zzz = pbuffer.data(idx_dip_if + 149);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxzz_xxx, tr_x_xxzz_xxy, tr_x_xxzz_xxz, tr_x_xxzz_xyy, tr_x_xxzz_xyz, tr_x_xxzz_xzz, tr_x_xxzz_yyy, tr_x_xxzzz_xx, tr_x_xxzzz_xxx, tr_x_xxzzz_xxy, tr_x_xxzzz_xxz, tr_x_xxzzz_xy, tr_x_xxzzz_xyy, tr_x_xxzzz_xyz, tr_x_xxzzz_xz, tr_x_xxzzz_xzz, tr_x_xxzzz_yyy, tr_x_xxzzzz_xxx, tr_x_xxzzzz_xxy, tr_x_xxzzzz_xxz, tr_x_xxzzzz_xyy, tr_x_xxzzzz_xyz, tr_x_xxzzzz_xzz, tr_x_xxzzzz_yyy, tr_x_xxzzzz_yyz, tr_x_xxzzzz_yzz, tr_x_xxzzzz_zzz, tr_x_xzzzz_yyz, tr_x_xzzzz_yzz, tr_x_xzzzz_zzz, tr_x_zzzz_yyz, tr_x_zzzz_yzz, tr_x_zzzz_zzz, ts_xzzzz_yyz, ts_xzzzz_yzz, ts_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzzz_xxx[i] = 3.0 * tr_x_xxzz_xxx[i] * fe_0 + tr_x_xxzzz_xxx[i] * pa_z[i];

        tr_x_xxzzzz_xxy[i] = 3.0 * tr_x_xxzz_xxy[i] * fe_0 + tr_x_xxzzz_xxy[i] * pa_z[i];

        tr_x_xxzzzz_xxz[i] = 3.0 * tr_x_xxzz_xxz[i] * fe_0 + tr_x_xxzzz_xx[i] * fe_0 + tr_x_xxzzz_xxz[i] * pa_z[i];

        tr_x_xxzzzz_xyy[i] = 3.0 * tr_x_xxzz_xyy[i] * fe_0 + tr_x_xxzzz_xyy[i] * pa_z[i];

        tr_x_xxzzzz_xyz[i] = 3.0 * tr_x_xxzz_xyz[i] * fe_0 + tr_x_xxzzz_xy[i] * fe_0 + tr_x_xxzzz_xyz[i] * pa_z[i];

        tr_x_xxzzzz_xzz[i] = 3.0 * tr_x_xxzz_xzz[i] * fe_0 + 2.0 * tr_x_xxzzz_xz[i] * fe_0 + tr_x_xxzzz_xzz[i] * pa_z[i];

        tr_x_xxzzzz_yyy[i] = 3.0 * tr_x_xxzz_yyy[i] * fe_0 + tr_x_xxzzz_yyy[i] * pa_z[i];

        tr_x_xxzzzz_yyz[i] = tr_x_zzzz_yyz[i] * fe_0 + ts_xzzzz_yyz[i] * fe_0 + tr_x_xzzzz_yyz[i] * pa_x[i];

        tr_x_xxzzzz_yzz[i] = tr_x_zzzz_yzz[i] * fe_0 + ts_xzzzz_yzz[i] * fe_0 + tr_x_xzzzz_yzz[i] * pa_x[i];

        tr_x_xxzzzz_zzz[i] = tr_x_zzzz_zzz[i] * fe_0 + ts_xzzzz_zzz[i] * fe_0 + tr_x_xzzzz_zzz[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : IF

    auto tr_x_xyyyyy_xxx = pbuffer.data(idx_dip_if + 150);

    auto tr_x_xyyyyy_xxy = pbuffer.data(idx_dip_if + 151);

    auto tr_x_xyyyyy_xxz = pbuffer.data(idx_dip_if + 152);

    auto tr_x_xyyyyy_xyy = pbuffer.data(idx_dip_if + 153);

    auto tr_x_xyyyyy_xyz = pbuffer.data(idx_dip_if + 154);

    auto tr_x_xyyyyy_xzz = pbuffer.data(idx_dip_if + 155);

    auto tr_x_xyyyyy_yyy = pbuffer.data(idx_dip_if + 156);

    auto tr_x_xyyyyy_yyz = pbuffer.data(idx_dip_if + 157);

    auto tr_x_xyyyyy_yzz = pbuffer.data(idx_dip_if + 158);

    auto tr_x_xyyyyy_zzz = pbuffer.data(idx_dip_if + 159);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyyy_xxx, tr_x_xyyy_xxz, tr_x_xyyy_xzz, tr_x_xyyyy_xxx, tr_x_xyyyy_xxz, tr_x_xyyyy_xzz, tr_x_xyyyyy_xxx, tr_x_xyyyyy_xxy, tr_x_xyyyyy_xxz, tr_x_xyyyyy_xyy, tr_x_xyyyyy_xyz, tr_x_xyyyyy_xzz, tr_x_xyyyyy_yyy, tr_x_xyyyyy_yyz, tr_x_xyyyyy_yzz, tr_x_xyyyyy_zzz, tr_x_yyyyy_xxy, tr_x_yyyyy_xy, tr_x_yyyyy_xyy, tr_x_yyyyy_xyz, tr_x_yyyyy_yy, tr_x_yyyyy_yyy, tr_x_yyyyy_yyz, tr_x_yyyyy_yz, tr_x_yyyyy_yzz, tr_x_yyyyy_zzz, ts_yyyyy_xxy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yzz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyy_xxx[i] = 4.0 * tr_x_xyyy_xxx[i] * fe_0 + tr_x_xyyyy_xxx[i] * pa_y[i];

        tr_x_xyyyyy_xxy[i] = 2.0 * tr_x_yyyyy_xy[i] * fe_0 + ts_yyyyy_xxy[i] * fe_0 + tr_x_yyyyy_xxy[i] * pa_x[i];

        tr_x_xyyyyy_xxz[i] = 4.0 * tr_x_xyyy_xxz[i] * fe_0 + tr_x_xyyyy_xxz[i] * pa_y[i];

        tr_x_xyyyyy_xyy[i] = tr_x_yyyyy_yy[i] * fe_0 + ts_yyyyy_xyy[i] * fe_0 + tr_x_yyyyy_xyy[i] * pa_x[i];

        tr_x_xyyyyy_xyz[i] = tr_x_yyyyy_yz[i] * fe_0 + ts_yyyyy_xyz[i] * fe_0 + tr_x_yyyyy_xyz[i] * pa_x[i];

        tr_x_xyyyyy_xzz[i] = 4.0 * tr_x_xyyy_xzz[i] * fe_0 + tr_x_xyyyy_xzz[i] * pa_y[i];

        tr_x_xyyyyy_yyy[i] = ts_yyyyy_yyy[i] * fe_0 + tr_x_yyyyy_yyy[i] * pa_x[i];

        tr_x_xyyyyy_yyz[i] = ts_yyyyy_yyz[i] * fe_0 + tr_x_yyyyy_yyz[i] * pa_x[i];

        tr_x_xyyyyy_yzz[i] = ts_yyyyy_yzz[i] * fe_0 + tr_x_yyyyy_yzz[i] * pa_x[i];

        tr_x_xyyyyy_zzz[i] = ts_yyyyy_zzz[i] * fe_0 + tr_x_yyyyy_zzz[i] * pa_x[i];
    }

    // Set up 160-170 components of targeted buffer : IF

    auto tr_x_xyyyyz_xxx = pbuffer.data(idx_dip_if + 160);

    auto tr_x_xyyyyz_xxy = pbuffer.data(idx_dip_if + 161);

    auto tr_x_xyyyyz_xxz = pbuffer.data(idx_dip_if + 162);

    auto tr_x_xyyyyz_xyy = pbuffer.data(idx_dip_if + 163);

    auto tr_x_xyyyyz_xyz = pbuffer.data(idx_dip_if + 164);

    auto tr_x_xyyyyz_xzz = pbuffer.data(idx_dip_if + 165);

    auto tr_x_xyyyyz_yyy = pbuffer.data(idx_dip_if + 166);

    auto tr_x_xyyyyz_yyz = pbuffer.data(idx_dip_if + 167);

    auto tr_x_xyyyyz_yzz = pbuffer.data(idx_dip_if + 168);

    auto tr_x_xyyyyz_zzz = pbuffer.data(idx_dip_if + 169);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyyy_xxx, tr_x_xyyyy_xxy, tr_x_xyyyy_xy, tr_x_xyyyy_xyy, tr_x_xyyyy_xyz, tr_x_xyyyy_yyy, tr_x_xyyyyz_xxx, tr_x_xyyyyz_xxy, tr_x_xyyyyz_xxz, tr_x_xyyyyz_xyy, tr_x_xyyyyz_xyz, tr_x_xyyyyz_xzz, tr_x_xyyyyz_yyy, tr_x_xyyyyz_yyz, tr_x_xyyyyz_yzz, tr_x_xyyyyz_zzz, tr_x_xyyyz_xxz, tr_x_xyyyz_xzz, tr_x_xyyz_xxz, tr_x_xyyz_xzz, tr_x_yyyyz_yyz, tr_x_yyyyz_yzz, tr_x_yyyyz_zzz, ts_yyyyz_yyz, ts_yyyyz_yzz, ts_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyz_xxx[i] = tr_x_xyyyy_xxx[i] * pa_z[i];

        tr_x_xyyyyz_xxy[i] = tr_x_xyyyy_xxy[i] * pa_z[i];

        tr_x_xyyyyz_xxz[i] = 3.0 * tr_x_xyyz_xxz[i] * fe_0 + tr_x_xyyyz_xxz[i] * pa_y[i];

        tr_x_xyyyyz_xyy[i] = tr_x_xyyyy_xyy[i] * pa_z[i];

        tr_x_xyyyyz_xyz[i] = tr_x_xyyyy_xy[i] * fe_0 + tr_x_xyyyy_xyz[i] * pa_z[i];

        tr_x_xyyyyz_xzz[i] = 3.0 * tr_x_xyyz_xzz[i] * fe_0 + tr_x_xyyyz_xzz[i] * pa_y[i];

        tr_x_xyyyyz_yyy[i] = tr_x_xyyyy_yyy[i] * pa_z[i];

        tr_x_xyyyyz_yyz[i] = ts_yyyyz_yyz[i] * fe_0 + tr_x_yyyyz_yyz[i] * pa_x[i];

        tr_x_xyyyyz_yzz[i] = ts_yyyyz_yzz[i] * fe_0 + tr_x_yyyyz_yzz[i] * pa_x[i];

        tr_x_xyyyyz_zzz[i] = ts_yyyyz_zzz[i] * fe_0 + tr_x_yyyyz_zzz[i] * pa_x[i];
    }

    // Set up 170-180 components of targeted buffer : IF

    auto tr_x_xyyyzz_xxx = pbuffer.data(idx_dip_if + 170);

    auto tr_x_xyyyzz_xxy = pbuffer.data(idx_dip_if + 171);

    auto tr_x_xyyyzz_xxz = pbuffer.data(idx_dip_if + 172);

    auto tr_x_xyyyzz_xyy = pbuffer.data(idx_dip_if + 173);

    auto tr_x_xyyyzz_xyz = pbuffer.data(idx_dip_if + 174);

    auto tr_x_xyyyzz_xzz = pbuffer.data(idx_dip_if + 175);

    auto tr_x_xyyyzz_yyy = pbuffer.data(idx_dip_if + 176);

    auto tr_x_xyyyzz_yyz = pbuffer.data(idx_dip_if + 177);

    auto tr_x_xyyyzz_yzz = pbuffer.data(idx_dip_if + 178);

    auto tr_x_xyyyzz_zzz = pbuffer.data(idx_dip_if + 179);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyy_xxy, tr_x_xyyy_xyy, tr_x_xyyyz_xxy, tr_x_xyyyz_xyy, tr_x_xyyyzz_xxx, tr_x_xyyyzz_xxy, tr_x_xyyyzz_xxz, tr_x_xyyyzz_xyy, tr_x_xyyyzz_xyz, tr_x_xyyyzz_xzz, tr_x_xyyyzz_yyy, tr_x_xyyyzz_yyz, tr_x_xyyyzz_yzz, tr_x_xyyyzz_zzz, tr_x_xyyzz_xxx, tr_x_xyyzz_xxz, tr_x_xyyzz_xzz, tr_x_xyzz_xxx, tr_x_xyzz_xxz, tr_x_xyzz_xzz, tr_x_yyyzz_xyz, tr_x_yyyzz_yyy, tr_x_yyyzz_yyz, tr_x_yyyzz_yz, tr_x_yyyzz_yzz, tr_x_yyyzz_zzz, ts_yyyzz_xyz, ts_yyyzz_yyy, ts_yyyzz_yyz, ts_yyyzz_yzz, ts_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyzz_xxx[i] = 2.0 * tr_x_xyzz_xxx[i] * fe_0 + tr_x_xyyzz_xxx[i] * pa_y[i];

        tr_x_xyyyzz_xxy[i] = tr_x_xyyy_xxy[i] * fe_0 + tr_x_xyyyz_xxy[i] * pa_z[i];

        tr_x_xyyyzz_xxz[i] = 2.0 * tr_x_xyzz_xxz[i] * fe_0 + tr_x_xyyzz_xxz[i] * pa_y[i];

        tr_x_xyyyzz_xyy[i] = tr_x_xyyy_xyy[i] * fe_0 + tr_x_xyyyz_xyy[i] * pa_z[i];

        tr_x_xyyyzz_xyz[i] = tr_x_yyyzz_yz[i] * fe_0 + ts_yyyzz_xyz[i] * fe_0 + tr_x_yyyzz_xyz[i] * pa_x[i];

        tr_x_xyyyzz_xzz[i] = 2.0 * tr_x_xyzz_xzz[i] * fe_0 + tr_x_xyyzz_xzz[i] * pa_y[i];

        tr_x_xyyyzz_yyy[i] = ts_yyyzz_yyy[i] * fe_0 + tr_x_yyyzz_yyy[i] * pa_x[i];

        tr_x_xyyyzz_yyz[i] = ts_yyyzz_yyz[i] * fe_0 + tr_x_yyyzz_yyz[i] * pa_x[i];

        tr_x_xyyyzz_yzz[i] = ts_yyyzz_yzz[i] * fe_0 + tr_x_yyyzz_yzz[i] * pa_x[i];

        tr_x_xyyyzz_zzz[i] = ts_yyyzz_zzz[i] * fe_0 + tr_x_yyyzz_zzz[i] * pa_x[i];
    }

    // Set up 180-190 components of targeted buffer : IF

    auto tr_x_xyyzzz_xxx = pbuffer.data(idx_dip_if + 180);

    auto tr_x_xyyzzz_xxy = pbuffer.data(idx_dip_if + 181);

    auto tr_x_xyyzzz_xxz = pbuffer.data(idx_dip_if + 182);

    auto tr_x_xyyzzz_xyy = pbuffer.data(idx_dip_if + 183);

    auto tr_x_xyyzzz_xyz = pbuffer.data(idx_dip_if + 184);

    auto tr_x_xyyzzz_xzz = pbuffer.data(idx_dip_if + 185);

    auto tr_x_xyyzzz_yyy = pbuffer.data(idx_dip_if + 186);

    auto tr_x_xyyzzz_yyz = pbuffer.data(idx_dip_if + 187);

    auto tr_x_xyyzzz_yzz = pbuffer.data(idx_dip_if + 188);

    auto tr_x_xyyzzz_zzz = pbuffer.data(idx_dip_if + 189);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyz_xxy, tr_x_xyyz_xyy, tr_x_xyyzz_xxy, tr_x_xyyzz_xyy, tr_x_xyyzzz_xxx, tr_x_xyyzzz_xxy, tr_x_xyyzzz_xxz, tr_x_xyyzzz_xyy, tr_x_xyyzzz_xyz, tr_x_xyyzzz_xzz, tr_x_xyyzzz_yyy, tr_x_xyyzzz_yyz, tr_x_xyyzzz_yzz, tr_x_xyyzzz_zzz, tr_x_xyzzz_xxx, tr_x_xyzzz_xxz, tr_x_xyzzz_xzz, tr_x_xzzz_xxx, tr_x_xzzz_xxz, tr_x_xzzz_xzz, tr_x_yyzzz_xyz, tr_x_yyzzz_yyy, tr_x_yyzzz_yyz, tr_x_yyzzz_yz, tr_x_yyzzz_yzz, tr_x_yyzzz_zzz, ts_yyzzz_xyz, ts_yyzzz_yyy, ts_yyzzz_yyz, ts_yyzzz_yzz, ts_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzzz_xxx[i] = tr_x_xzzz_xxx[i] * fe_0 + tr_x_xyzzz_xxx[i] * pa_y[i];

        tr_x_xyyzzz_xxy[i] = 2.0 * tr_x_xyyz_xxy[i] * fe_0 + tr_x_xyyzz_xxy[i] * pa_z[i];

        tr_x_xyyzzz_xxz[i] = tr_x_xzzz_xxz[i] * fe_0 + tr_x_xyzzz_xxz[i] * pa_y[i];

        tr_x_xyyzzz_xyy[i] = 2.0 * tr_x_xyyz_xyy[i] * fe_0 + tr_x_xyyzz_xyy[i] * pa_z[i];

        tr_x_xyyzzz_xyz[i] = tr_x_yyzzz_yz[i] * fe_0 + ts_yyzzz_xyz[i] * fe_0 + tr_x_yyzzz_xyz[i] * pa_x[i];

        tr_x_xyyzzz_xzz[i] = tr_x_xzzz_xzz[i] * fe_0 + tr_x_xyzzz_xzz[i] * pa_y[i];

        tr_x_xyyzzz_yyy[i] = ts_yyzzz_yyy[i] * fe_0 + tr_x_yyzzz_yyy[i] * pa_x[i];

        tr_x_xyyzzz_yyz[i] = ts_yyzzz_yyz[i] * fe_0 + tr_x_yyzzz_yyz[i] * pa_x[i];

        tr_x_xyyzzz_yzz[i] = ts_yyzzz_yzz[i] * fe_0 + tr_x_yyzzz_yzz[i] * pa_x[i];

        tr_x_xyyzzz_zzz[i] = ts_yyzzz_zzz[i] * fe_0 + tr_x_yyzzz_zzz[i] * pa_x[i];
    }

    // Set up 190-200 components of targeted buffer : IF

    auto tr_x_xyzzzz_xxx = pbuffer.data(idx_dip_if + 190);

    auto tr_x_xyzzzz_xxy = pbuffer.data(idx_dip_if + 191);

    auto tr_x_xyzzzz_xxz = pbuffer.data(idx_dip_if + 192);

    auto tr_x_xyzzzz_xyy = pbuffer.data(idx_dip_if + 193);

    auto tr_x_xyzzzz_xyz = pbuffer.data(idx_dip_if + 194);

    auto tr_x_xyzzzz_xzz = pbuffer.data(idx_dip_if + 195);

    auto tr_x_xyzzzz_yyy = pbuffer.data(idx_dip_if + 196);

    auto tr_x_xyzzzz_yyz = pbuffer.data(idx_dip_if + 197);

    auto tr_x_xyzzzz_yzz = pbuffer.data(idx_dip_if + 198);

    auto tr_x_xyzzzz_zzz = pbuffer.data(idx_dip_if + 199);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyzzzz_xxx, tr_x_xyzzzz_xxy, tr_x_xyzzzz_xxz, tr_x_xyzzzz_xyy, tr_x_xyzzzz_xyz, tr_x_xyzzzz_xzz, tr_x_xyzzzz_yyy, tr_x_xyzzzz_yyz, tr_x_xyzzzz_yzz, tr_x_xyzzzz_zzz, tr_x_xzzzz_xx, tr_x_xzzzz_xxx, tr_x_xzzzz_xxy, tr_x_xzzzz_xxz, tr_x_xzzzz_xy, tr_x_xzzzz_xyy, tr_x_xzzzz_xyz, tr_x_xzzzz_xz, tr_x_xzzzz_xzz, tr_x_xzzzz_zzz, tr_x_yzzzz_yyy, tr_x_yzzzz_yyz, tr_x_yzzzz_yzz, ts_yzzzz_yyy, ts_yzzzz_yyz, ts_yzzzz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzzz_xxx[i] = tr_x_xzzzz_xxx[i] * pa_y[i];

        tr_x_xyzzzz_xxy[i] = tr_x_xzzzz_xx[i] * fe_0 + tr_x_xzzzz_xxy[i] * pa_y[i];

        tr_x_xyzzzz_xxz[i] = tr_x_xzzzz_xxz[i] * pa_y[i];

        tr_x_xyzzzz_xyy[i] = 2.0 * tr_x_xzzzz_xy[i] * fe_0 + tr_x_xzzzz_xyy[i] * pa_y[i];

        tr_x_xyzzzz_xyz[i] = tr_x_xzzzz_xz[i] * fe_0 + tr_x_xzzzz_xyz[i] * pa_y[i];

        tr_x_xyzzzz_xzz[i] = tr_x_xzzzz_xzz[i] * pa_y[i];

        tr_x_xyzzzz_yyy[i] = ts_yzzzz_yyy[i] * fe_0 + tr_x_yzzzz_yyy[i] * pa_x[i];

        tr_x_xyzzzz_yyz[i] = ts_yzzzz_yyz[i] * fe_0 + tr_x_yzzzz_yyz[i] * pa_x[i];

        tr_x_xyzzzz_yzz[i] = ts_yzzzz_yzz[i] * fe_0 + tr_x_yzzzz_yzz[i] * pa_x[i];

        tr_x_xyzzzz_zzz[i] = tr_x_xzzzz_zzz[i] * pa_y[i];
    }

    // Set up 200-210 components of targeted buffer : IF

    auto tr_x_xzzzzz_xxx = pbuffer.data(idx_dip_if + 200);

    auto tr_x_xzzzzz_xxy = pbuffer.data(idx_dip_if + 201);

    auto tr_x_xzzzzz_xxz = pbuffer.data(idx_dip_if + 202);

    auto tr_x_xzzzzz_xyy = pbuffer.data(idx_dip_if + 203);

    auto tr_x_xzzzzz_xyz = pbuffer.data(idx_dip_if + 204);

    auto tr_x_xzzzzz_xzz = pbuffer.data(idx_dip_if + 205);

    auto tr_x_xzzzzz_yyy = pbuffer.data(idx_dip_if + 206);

    auto tr_x_xzzzzz_yyz = pbuffer.data(idx_dip_if + 207);

    auto tr_x_xzzzzz_yzz = pbuffer.data(idx_dip_if + 208);

    auto tr_x_xzzzzz_zzz = pbuffer.data(idx_dip_if + 209);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xzzz_xxx, tr_x_xzzz_xxy, tr_x_xzzz_xyy, tr_x_xzzzz_xxx, tr_x_xzzzz_xxy, tr_x_xzzzz_xyy, tr_x_xzzzzz_xxx, tr_x_xzzzzz_xxy, tr_x_xzzzzz_xxz, tr_x_xzzzzz_xyy, tr_x_xzzzzz_xyz, tr_x_xzzzzz_xzz, tr_x_xzzzzz_yyy, tr_x_xzzzzz_yyz, tr_x_xzzzzz_yzz, tr_x_xzzzzz_zzz, tr_x_zzzzz_xxz, tr_x_zzzzz_xyz, tr_x_zzzzz_xz, tr_x_zzzzz_xzz, tr_x_zzzzz_yyy, tr_x_zzzzz_yyz, tr_x_zzzzz_yz, tr_x_zzzzz_yzz, tr_x_zzzzz_zz, tr_x_zzzzz_zzz, ts_zzzzz_xxz, ts_zzzzz_xyz, ts_zzzzz_xzz, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yzz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzzz_xxx[i] = 4.0 * tr_x_xzzz_xxx[i] * fe_0 + tr_x_xzzzz_xxx[i] * pa_z[i];

        tr_x_xzzzzz_xxy[i] = 4.0 * tr_x_xzzz_xxy[i] * fe_0 + tr_x_xzzzz_xxy[i] * pa_z[i];

        tr_x_xzzzzz_xxz[i] = 2.0 * tr_x_zzzzz_xz[i] * fe_0 + ts_zzzzz_xxz[i] * fe_0 + tr_x_zzzzz_xxz[i] * pa_x[i];

        tr_x_xzzzzz_xyy[i] = 4.0 * tr_x_xzzz_xyy[i] * fe_0 + tr_x_xzzzz_xyy[i] * pa_z[i];

        tr_x_xzzzzz_xyz[i] = tr_x_zzzzz_yz[i] * fe_0 + ts_zzzzz_xyz[i] * fe_0 + tr_x_zzzzz_xyz[i] * pa_x[i];

        tr_x_xzzzzz_xzz[i] = tr_x_zzzzz_zz[i] * fe_0 + ts_zzzzz_xzz[i] * fe_0 + tr_x_zzzzz_xzz[i] * pa_x[i];

        tr_x_xzzzzz_yyy[i] = ts_zzzzz_yyy[i] * fe_0 + tr_x_zzzzz_yyy[i] * pa_x[i];

        tr_x_xzzzzz_yyz[i] = ts_zzzzz_yyz[i] * fe_0 + tr_x_zzzzz_yyz[i] * pa_x[i];

        tr_x_xzzzzz_yzz[i] = ts_zzzzz_yzz[i] * fe_0 + tr_x_zzzzz_yzz[i] * pa_x[i];

        tr_x_xzzzzz_zzz[i] = ts_zzzzz_zzz[i] * fe_0 + tr_x_zzzzz_zzz[i] * pa_x[i];
    }

    // Set up 210-220 components of targeted buffer : IF

    auto tr_x_yyyyyy_xxx = pbuffer.data(idx_dip_if + 210);

    auto tr_x_yyyyyy_xxy = pbuffer.data(idx_dip_if + 211);

    auto tr_x_yyyyyy_xxz = pbuffer.data(idx_dip_if + 212);

    auto tr_x_yyyyyy_xyy = pbuffer.data(idx_dip_if + 213);

    auto tr_x_yyyyyy_xyz = pbuffer.data(idx_dip_if + 214);

    auto tr_x_yyyyyy_xzz = pbuffer.data(idx_dip_if + 215);

    auto tr_x_yyyyyy_yyy = pbuffer.data(idx_dip_if + 216);

    auto tr_x_yyyyyy_yyz = pbuffer.data(idx_dip_if + 217);

    auto tr_x_yyyyyy_yzz = pbuffer.data(idx_dip_if + 218);

    auto tr_x_yyyyyy_zzz = pbuffer.data(idx_dip_if + 219);

    #pragma omp simd aligned(pa_y, tr_x_yyyy_xxx, tr_x_yyyy_xxy, tr_x_yyyy_xxz, tr_x_yyyy_xyy, tr_x_yyyy_xyz, tr_x_yyyy_xzz, tr_x_yyyy_yyy, tr_x_yyyy_yyz, tr_x_yyyy_yzz, tr_x_yyyy_zzz, tr_x_yyyyy_xx, tr_x_yyyyy_xxx, tr_x_yyyyy_xxy, tr_x_yyyyy_xxz, tr_x_yyyyy_xy, tr_x_yyyyy_xyy, tr_x_yyyyy_xyz, tr_x_yyyyy_xz, tr_x_yyyyy_xzz, tr_x_yyyyy_yy, tr_x_yyyyy_yyy, tr_x_yyyyy_yyz, tr_x_yyyyy_yz, tr_x_yyyyy_yzz, tr_x_yyyyy_zz, tr_x_yyyyy_zzz, tr_x_yyyyyy_xxx, tr_x_yyyyyy_xxy, tr_x_yyyyyy_xxz, tr_x_yyyyyy_xyy, tr_x_yyyyyy_xyz, tr_x_yyyyyy_xzz, tr_x_yyyyyy_yyy, tr_x_yyyyyy_yyz, tr_x_yyyyyy_yzz, tr_x_yyyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyy_xxx[i] = 5.0 * tr_x_yyyy_xxx[i] * fe_0 + tr_x_yyyyy_xxx[i] * pa_y[i];

        tr_x_yyyyyy_xxy[i] = 5.0 * tr_x_yyyy_xxy[i] * fe_0 + tr_x_yyyyy_xx[i] * fe_0 + tr_x_yyyyy_xxy[i] * pa_y[i];

        tr_x_yyyyyy_xxz[i] = 5.0 * tr_x_yyyy_xxz[i] * fe_0 + tr_x_yyyyy_xxz[i] * pa_y[i];

        tr_x_yyyyyy_xyy[i] = 5.0 * tr_x_yyyy_xyy[i] * fe_0 + 2.0 * tr_x_yyyyy_xy[i] * fe_0 + tr_x_yyyyy_xyy[i] * pa_y[i];

        tr_x_yyyyyy_xyz[i] = 5.0 * tr_x_yyyy_xyz[i] * fe_0 + tr_x_yyyyy_xz[i] * fe_0 + tr_x_yyyyy_xyz[i] * pa_y[i];

        tr_x_yyyyyy_xzz[i] = 5.0 * tr_x_yyyy_xzz[i] * fe_0 + tr_x_yyyyy_xzz[i] * pa_y[i];

        tr_x_yyyyyy_yyy[i] = 5.0 * tr_x_yyyy_yyy[i] * fe_0 + 3.0 * tr_x_yyyyy_yy[i] * fe_0 + tr_x_yyyyy_yyy[i] * pa_y[i];

        tr_x_yyyyyy_yyz[i] = 5.0 * tr_x_yyyy_yyz[i] * fe_0 + 2.0 * tr_x_yyyyy_yz[i] * fe_0 + tr_x_yyyyy_yyz[i] * pa_y[i];

        tr_x_yyyyyy_yzz[i] = 5.0 * tr_x_yyyy_yzz[i] * fe_0 + tr_x_yyyyy_zz[i] * fe_0 + tr_x_yyyyy_yzz[i] * pa_y[i];

        tr_x_yyyyyy_zzz[i] = 5.0 * tr_x_yyyy_zzz[i] * fe_0 + tr_x_yyyyy_zzz[i] * pa_y[i];
    }

    // Set up 220-230 components of targeted buffer : IF

    auto tr_x_yyyyyz_xxx = pbuffer.data(idx_dip_if + 220);

    auto tr_x_yyyyyz_xxy = pbuffer.data(idx_dip_if + 221);

    auto tr_x_yyyyyz_xxz = pbuffer.data(idx_dip_if + 222);

    auto tr_x_yyyyyz_xyy = pbuffer.data(idx_dip_if + 223);

    auto tr_x_yyyyyz_xyz = pbuffer.data(idx_dip_if + 224);

    auto tr_x_yyyyyz_xzz = pbuffer.data(idx_dip_if + 225);

    auto tr_x_yyyyyz_yyy = pbuffer.data(idx_dip_if + 226);

    auto tr_x_yyyyyz_yyz = pbuffer.data(idx_dip_if + 227);

    auto tr_x_yyyyyz_yzz = pbuffer.data(idx_dip_if + 228);

    auto tr_x_yyyyyz_zzz = pbuffer.data(idx_dip_if + 229);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyyy_xxx, tr_x_yyyyy_xxy, tr_x_yyyyy_xy, tr_x_yyyyy_xyy, tr_x_yyyyy_xyz, tr_x_yyyyy_yy, tr_x_yyyyy_yyy, tr_x_yyyyy_yyz, tr_x_yyyyy_yz, tr_x_yyyyy_yzz, tr_x_yyyyyz_xxx, tr_x_yyyyyz_xxy, tr_x_yyyyyz_xxz, tr_x_yyyyyz_xyy, tr_x_yyyyyz_xyz, tr_x_yyyyyz_xzz, tr_x_yyyyyz_yyy, tr_x_yyyyyz_yyz, tr_x_yyyyyz_yzz, tr_x_yyyyyz_zzz, tr_x_yyyyz_xxz, tr_x_yyyyz_xzz, tr_x_yyyyz_zzz, tr_x_yyyz_xxz, tr_x_yyyz_xzz, tr_x_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyz_xxx[i] = tr_x_yyyyy_xxx[i] * pa_z[i];

        tr_x_yyyyyz_xxy[i] = tr_x_yyyyy_xxy[i] * pa_z[i];

        tr_x_yyyyyz_xxz[i] = 4.0 * tr_x_yyyz_xxz[i] * fe_0 + tr_x_yyyyz_xxz[i] * pa_y[i];

        tr_x_yyyyyz_xyy[i] = tr_x_yyyyy_xyy[i] * pa_z[i];

        tr_x_yyyyyz_xyz[i] = tr_x_yyyyy_xy[i] * fe_0 + tr_x_yyyyy_xyz[i] * pa_z[i];

        tr_x_yyyyyz_xzz[i] = 4.0 * tr_x_yyyz_xzz[i] * fe_0 + tr_x_yyyyz_xzz[i] * pa_y[i];

        tr_x_yyyyyz_yyy[i] = tr_x_yyyyy_yyy[i] * pa_z[i];

        tr_x_yyyyyz_yyz[i] = tr_x_yyyyy_yy[i] * fe_0 + tr_x_yyyyy_yyz[i] * pa_z[i];

        tr_x_yyyyyz_yzz[i] = 2.0 * tr_x_yyyyy_yz[i] * fe_0 + tr_x_yyyyy_yzz[i] * pa_z[i];

        tr_x_yyyyyz_zzz[i] = 4.0 * tr_x_yyyz_zzz[i] * fe_0 + tr_x_yyyyz_zzz[i] * pa_y[i];
    }

    // Set up 230-240 components of targeted buffer : IF

    auto tr_x_yyyyzz_xxx = pbuffer.data(idx_dip_if + 230);

    auto tr_x_yyyyzz_xxy = pbuffer.data(idx_dip_if + 231);

    auto tr_x_yyyyzz_xxz = pbuffer.data(idx_dip_if + 232);

    auto tr_x_yyyyzz_xyy = pbuffer.data(idx_dip_if + 233);

    auto tr_x_yyyyzz_xyz = pbuffer.data(idx_dip_if + 234);

    auto tr_x_yyyyzz_xzz = pbuffer.data(idx_dip_if + 235);

    auto tr_x_yyyyzz_yyy = pbuffer.data(idx_dip_if + 236);

    auto tr_x_yyyyzz_yyz = pbuffer.data(idx_dip_if + 237);

    auto tr_x_yyyyzz_yzz = pbuffer.data(idx_dip_if + 238);

    auto tr_x_yyyyzz_zzz = pbuffer.data(idx_dip_if + 239);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyy_xxy, tr_x_yyyy_xyy, tr_x_yyyy_yyy, tr_x_yyyyz_xxy, tr_x_yyyyz_xyy, tr_x_yyyyz_yyy, tr_x_yyyyzz_xxx, tr_x_yyyyzz_xxy, tr_x_yyyyzz_xxz, tr_x_yyyyzz_xyy, tr_x_yyyyzz_xyz, tr_x_yyyyzz_xzz, tr_x_yyyyzz_yyy, tr_x_yyyyzz_yyz, tr_x_yyyyzz_yzz, tr_x_yyyyzz_zzz, tr_x_yyyzz_xxx, tr_x_yyyzz_xxz, tr_x_yyyzz_xyz, tr_x_yyyzz_xz, tr_x_yyyzz_xzz, tr_x_yyyzz_yyz, tr_x_yyyzz_yz, tr_x_yyyzz_yzz, tr_x_yyyzz_zz, tr_x_yyyzz_zzz, tr_x_yyzz_xxx, tr_x_yyzz_xxz, tr_x_yyzz_xyz, tr_x_yyzz_xzz, tr_x_yyzz_yyz, tr_x_yyzz_yzz, tr_x_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyzz_xxx[i] = 3.0 * tr_x_yyzz_xxx[i] * fe_0 + tr_x_yyyzz_xxx[i] * pa_y[i];

        tr_x_yyyyzz_xxy[i] = tr_x_yyyy_xxy[i] * fe_0 + tr_x_yyyyz_xxy[i] * pa_z[i];

        tr_x_yyyyzz_xxz[i] = 3.0 * tr_x_yyzz_xxz[i] * fe_0 + tr_x_yyyzz_xxz[i] * pa_y[i];

        tr_x_yyyyzz_xyy[i] = tr_x_yyyy_xyy[i] * fe_0 + tr_x_yyyyz_xyy[i] * pa_z[i];

        tr_x_yyyyzz_xyz[i] = 3.0 * tr_x_yyzz_xyz[i] * fe_0 + tr_x_yyyzz_xz[i] * fe_0 + tr_x_yyyzz_xyz[i] * pa_y[i];

        tr_x_yyyyzz_xzz[i] = 3.0 * tr_x_yyzz_xzz[i] * fe_0 + tr_x_yyyzz_xzz[i] * pa_y[i];

        tr_x_yyyyzz_yyy[i] = tr_x_yyyy_yyy[i] * fe_0 + tr_x_yyyyz_yyy[i] * pa_z[i];

        tr_x_yyyyzz_yyz[i] = 3.0 * tr_x_yyzz_yyz[i] * fe_0 + 2.0 * tr_x_yyyzz_yz[i] * fe_0 + tr_x_yyyzz_yyz[i] * pa_y[i];

        tr_x_yyyyzz_yzz[i] = 3.0 * tr_x_yyzz_yzz[i] * fe_0 + tr_x_yyyzz_zz[i] * fe_0 + tr_x_yyyzz_yzz[i] * pa_y[i];

        tr_x_yyyyzz_zzz[i] = 3.0 * tr_x_yyzz_zzz[i] * fe_0 + tr_x_yyyzz_zzz[i] * pa_y[i];
    }

    // Set up 240-250 components of targeted buffer : IF

    auto tr_x_yyyzzz_xxx = pbuffer.data(idx_dip_if + 240);

    auto tr_x_yyyzzz_xxy = pbuffer.data(idx_dip_if + 241);

    auto tr_x_yyyzzz_xxz = pbuffer.data(idx_dip_if + 242);

    auto tr_x_yyyzzz_xyy = pbuffer.data(idx_dip_if + 243);

    auto tr_x_yyyzzz_xyz = pbuffer.data(idx_dip_if + 244);

    auto tr_x_yyyzzz_xzz = pbuffer.data(idx_dip_if + 245);

    auto tr_x_yyyzzz_yyy = pbuffer.data(idx_dip_if + 246);

    auto tr_x_yyyzzz_yyz = pbuffer.data(idx_dip_if + 247);

    auto tr_x_yyyzzz_yzz = pbuffer.data(idx_dip_if + 248);

    auto tr_x_yyyzzz_zzz = pbuffer.data(idx_dip_if + 249);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyz_xxy, tr_x_yyyz_xyy, tr_x_yyyz_yyy, tr_x_yyyzz_xxy, tr_x_yyyzz_xyy, tr_x_yyyzz_yyy, tr_x_yyyzzz_xxx, tr_x_yyyzzz_xxy, tr_x_yyyzzz_xxz, tr_x_yyyzzz_xyy, tr_x_yyyzzz_xyz, tr_x_yyyzzz_xzz, tr_x_yyyzzz_yyy, tr_x_yyyzzz_yyz, tr_x_yyyzzz_yzz, tr_x_yyyzzz_zzz, tr_x_yyzzz_xxx, tr_x_yyzzz_xxz, tr_x_yyzzz_xyz, tr_x_yyzzz_xz, tr_x_yyzzz_xzz, tr_x_yyzzz_yyz, tr_x_yyzzz_yz, tr_x_yyzzz_yzz, tr_x_yyzzz_zz, tr_x_yyzzz_zzz, tr_x_yzzz_xxx, tr_x_yzzz_xxz, tr_x_yzzz_xyz, tr_x_yzzz_xzz, tr_x_yzzz_yyz, tr_x_yzzz_yzz, tr_x_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzzz_xxx[i] = 2.0 * tr_x_yzzz_xxx[i] * fe_0 + tr_x_yyzzz_xxx[i] * pa_y[i];

        tr_x_yyyzzz_xxy[i] = 2.0 * tr_x_yyyz_xxy[i] * fe_0 + tr_x_yyyzz_xxy[i] * pa_z[i];

        tr_x_yyyzzz_xxz[i] = 2.0 * tr_x_yzzz_xxz[i] * fe_0 + tr_x_yyzzz_xxz[i] * pa_y[i];

        tr_x_yyyzzz_xyy[i] = 2.0 * tr_x_yyyz_xyy[i] * fe_0 + tr_x_yyyzz_xyy[i] * pa_z[i];

        tr_x_yyyzzz_xyz[i] = 2.0 * tr_x_yzzz_xyz[i] * fe_0 + tr_x_yyzzz_xz[i] * fe_0 + tr_x_yyzzz_xyz[i] * pa_y[i];

        tr_x_yyyzzz_xzz[i] = 2.0 * tr_x_yzzz_xzz[i] * fe_0 + tr_x_yyzzz_xzz[i] * pa_y[i];

        tr_x_yyyzzz_yyy[i] = 2.0 * tr_x_yyyz_yyy[i] * fe_0 + tr_x_yyyzz_yyy[i] * pa_z[i];

        tr_x_yyyzzz_yyz[i] = 2.0 * tr_x_yzzz_yyz[i] * fe_0 + 2.0 * tr_x_yyzzz_yz[i] * fe_0 + tr_x_yyzzz_yyz[i] * pa_y[i];

        tr_x_yyyzzz_yzz[i] = 2.0 * tr_x_yzzz_yzz[i] * fe_0 + tr_x_yyzzz_zz[i] * fe_0 + tr_x_yyzzz_yzz[i] * pa_y[i];

        tr_x_yyyzzz_zzz[i] = 2.0 * tr_x_yzzz_zzz[i] * fe_0 + tr_x_yyzzz_zzz[i] * pa_y[i];
    }

    // Set up 250-260 components of targeted buffer : IF

    auto tr_x_yyzzzz_xxx = pbuffer.data(idx_dip_if + 250);

    auto tr_x_yyzzzz_xxy = pbuffer.data(idx_dip_if + 251);

    auto tr_x_yyzzzz_xxz = pbuffer.data(idx_dip_if + 252);

    auto tr_x_yyzzzz_xyy = pbuffer.data(idx_dip_if + 253);

    auto tr_x_yyzzzz_xyz = pbuffer.data(idx_dip_if + 254);

    auto tr_x_yyzzzz_xzz = pbuffer.data(idx_dip_if + 255);

    auto tr_x_yyzzzz_yyy = pbuffer.data(idx_dip_if + 256);

    auto tr_x_yyzzzz_yyz = pbuffer.data(idx_dip_if + 257);

    auto tr_x_yyzzzz_yzz = pbuffer.data(idx_dip_if + 258);

    auto tr_x_yyzzzz_zzz = pbuffer.data(idx_dip_if + 259);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyzz_xxy, tr_x_yyzz_xyy, tr_x_yyzz_yyy, tr_x_yyzzz_xxy, tr_x_yyzzz_xyy, tr_x_yyzzz_yyy, tr_x_yyzzzz_xxx, tr_x_yyzzzz_xxy, tr_x_yyzzzz_xxz, tr_x_yyzzzz_xyy, tr_x_yyzzzz_xyz, tr_x_yyzzzz_xzz, tr_x_yyzzzz_yyy, tr_x_yyzzzz_yyz, tr_x_yyzzzz_yzz, tr_x_yyzzzz_zzz, tr_x_yzzzz_xxx, tr_x_yzzzz_xxz, tr_x_yzzzz_xyz, tr_x_yzzzz_xz, tr_x_yzzzz_xzz, tr_x_yzzzz_yyz, tr_x_yzzzz_yz, tr_x_yzzzz_yzz, tr_x_yzzzz_zz, tr_x_yzzzz_zzz, tr_x_zzzz_xxx, tr_x_zzzz_xxz, tr_x_zzzz_xyz, tr_x_zzzz_xzz, tr_x_zzzz_yyz, tr_x_zzzz_yzz, tr_x_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzzz_xxx[i] = tr_x_zzzz_xxx[i] * fe_0 + tr_x_yzzzz_xxx[i] * pa_y[i];

        tr_x_yyzzzz_xxy[i] = 3.0 * tr_x_yyzz_xxy[i] * fe_0 + tr_x_yyzzz_xxy[i] * pa_z[i];

        tr_x_yyzzzz_xxz[i] = tr_x_zzzz_xxz[i] * fe_0 + tr_x_yzzzz_xxz[i] * pa_y[i];

        tr_x_yyzzzz_xyy[i] = 3.0 * tr_x_yyzz_xyy[i] * fe_0 + tr_x_yyzzz_xyy[i] * pa_z[i];

        tr_x_yyzzzz_xyz[i] = tr_x_zzzz_xyz[i] * fe_0 + tr_x_yzzzz_xz[i] * fe_0 + tr_x_yzzzz_xyz[i] * pa_y[i];

        tr_x_yyzzzz_xzz[i] = tr_x_zzzz_xzz[i] * fe_0 + tr_x_yzzzz_xzz[i] * pa_y[i];

        tr_x_yyzzzz_yyy[i] = 3.0 * tr_x_yyzz_yyy[i] * fe_0 + tr_x_yyzzz_yyy[i] * pa_z[i];

        tr_x_yyzzzz_yyz[i] = tr_x_zzzz_yyz[i] * fe_0 + 2.0 * tr_x_yzzzz_yz[i] * fe_0 + tr_x_yzzzz_yyz[i] * pa_y[i];

        tr_x_yyzzzz_yzz[i] = tr_x_zzzz_yzz[i] * fe_0 + tr_x_yzzzz_zz[i] * fe_0 + tr_x_yzzzz_yzz[i] * pa_y[i];

        tr_x_yyzzzz_zzz[i] = tr_x_zzzz_zzz[i] * fe_0 + tr_x_yzzzz_zzz[i] * pa_y[i];
    }

    // Set up 260-270 components of targeted buffer : IF

    auto tr_x_yzzzzz_xxx = pbuffer.data(idx_dip_if + 260);

    auto tr_x_yzzzzz_xxy = pbuffer.data(idx_dip_if + 261);

    auto tr_x_yzzzzz_xxz = pbuffer.data(idx_dip_if + 262);

    auto tr_x_yzzzzz_xyy = pbuffer.data(idx_dip_if + 263);

    auto tr_x_yzzzzz_xyz = pbuffer.data(idx_dip_if + 264);

    auto tr_x_yzzzzz_xzz = pbuffer.data(idx_dip_if + 265);

    auto tr_x_yzzzzz_yyy = pbuffer.data(idx_dip_if + 266);

    auto tr_x_yzzzzz_yyz = pbuffer.data(idx_dip_if + 267);

    auto tr_x_yzzzzz_yzz = pbuffer.data(idx_dip_if + 268);

    auto tr_x_yzzzzz_zzz = pbuffer.data(idx_dip_if + 269);

    #pragma omp simd aligned(pa_y, tr_x_yzzzzz_xxx, tr_x_yzzzzz_xxy, tr_x_yzzzzz_xxz, tr_x_yzzzzz_xyy, tr_x_yzzzzz_xyz, tr_x_yzzzzz_xzz, tr_x_yzzzzz_yyy, tr_x_yzzzzz_yyz, tr_x_yzzzzz_yzz, tr_x_yzzzzz_zzz, tr_x_zzzzz_xx, tr_x_zzzzz_xxx, tr_x_zzzzz_xxy, tr_x_zzzzz_xxz, tr_x_zzzzz_xy, tr_x_zzzzz_xyy, tr_x_zzzzz_xyz, tr_x_zzzzz_xz, tr_x_zzzzz_xzz, tr_x_zzzzz_yy, tr_x_zzzzz_yyy, tr_x_zzzzz_yyz, tr_x_zzzzz_yz, tr_x_zzzzz_yzz, tr_x_zzzzz_zz, tr_x_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzzz_xxx[i] = tr_x_zzzzz_xxx[i] * pa_y[i];

        tr_x_yzzzzz_xxy[i] = tr_x_zzzzz_xx[i] * fe_0 + tr_x_zzzzz_xxy[i] * pa_y[i];

        tr_x_yzzzzz_xxz[i] = tr_x_zzzzz_xxz[i] * pa_y[i];

        tr_x_yzzzzz_xyy[i] = 2.0 * tr_x_zzzzz_xy[i] * fe_0 + tr_x_zzzzz_xyy[i] * pa_y[i];

        tr_x_yzzzzz_xyz[i] = tr_x_zzzzz_xz[i] * fe_0 + tr_x_zzzzz_xyz[i] * pa_y[i];

        tr_x_yzzzzz_xzz[i] = tr_x_zzzzz_xzz[i] * pa_y[i];

        tr_x_yzzzzz_yyy[i] = 3.0 * tr_x_zzzzz_yy[i] * fe_0 + tr_x_zzzzz_yyy[i] * pa_y[i];

        tr_x_yzzzzz_yyz[i] = 2.0 * tr_x_zzzzz_yz[i] * fe_0 + tr_x_zzzzz_yyz[i] * pa_y[i];

        tr_x_yzzzzz_yzz[i] = tr_x_zzzzz_zz[i] * fe_0 + tr_x_zzzzz_yzz[i] * pa_y[i];

        tr_x_yzzzzz_zzz[i] = tr_x_zzzzz_zzz[i] * pa_y[i];
    }

    // Set up 270-280 components of targeted buffer : IF

    auto tr_x_zzzzzz_xxx = pbuffer.data(idx_dip_if + 270);

    auto tr_x_zzzzzz_xxy = pbuffer.data(idx_dip_if + 271);

    auto tr_x_zzzzzz_xxz = pbuffer.data(idx_dip_if + 272);

    auto tr_x_zzzzzz_xyy = pbuffer.data(idx_dip_if + 273);

    auto tr_x_zzzzzz_xyz = pbuffer.data(idx_dip_if + 274);

    auto tr_x_zzzzzz_xzz = pbuffer.data(idx_dip_if + 275);

    auto tr_x_zzzzzz_yyy = pbuffer.data(idx_dip_if + 276);

    auto tr_x_zzzzzz_yyz = pbuffer.data(idx_dip_if + 277);

    auto tr_x_zzzzzz_yzz = pbuffer.data(idx_dip_if + 278);

    auto tr_x_zzzzzz_zzz = pbuffer.data(idx_dip_if + 279);

    #pragma omp simd aligned(pa_z, tr_x_zzzz_xxx, tr_x_zzzz_xxy, tr_x_zzzz_xxz, tr_x_zzzz_xyy, tr_x_zzzz_xyz, tr_x_zzzz_xzz, tr_x_zzzz_yyy, tr_x_zzzz_yyz, tr_x_zzzz_yzz, tr_x_zzzz_zzz, tr_x_zzzzz_xx, tr_x_zzzzz_xxx, tr_x_zzzzz_xxy, tr_x_zzzzz_xxz, tr_x_zzzzz_xy, tr_x_zzzzz_xyy, tr_x_zzzzz_xyz, tr_x_zzzzz_xz, tr_x_zzzzz_xzz, tr_x_zzzzz_yy, tr_x_zzzzz_yyy, tr_x_zzzzz_yyz, tr_x_zzzzz_yz, tr_x_zzzzz_yzz, tr_x_zzzzz_zz, tr_x_zzzzz_zzz, tr_x_zzzzzz_xxx, tr_x_zzzzzz_xxy, tr_x_zzzzzz_xxz, tr_x_zzzzzz_xyy, tr_x_zzzzzz_xyz, tr_x_zzzzzz_xzz, tr_x_zzzzzz_yyy, tr_x_zzzzzz_yyz, tr_x_zzzzzz_yzz, tr_x_zzzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzzz_xxx[i] = 5.0 * tr_x_zzzz_xxx[i] * fe_0 + tr_x_zzzzz_xxx[i] * pa_z[i];

        tr_x_zzzzzz_xxy[i] = 5.0 * tr_x_zzzz_xxy[i] * fe_0 + tr_x_zzzzz_xxy[i] * pa_z[i];

        tr_x_zzzzzz_xxz[i] = 5.0 * tr_x_zzzz_xxz[i] * fe_0 + tr_x_zzzzz_xx[i] * fe_0 + tr_x_zzzzz_xxz[i] * pa_z[i];

        tr_x_zzzzzz_xyy[i] = 5.0 * tr_x_zzzz_xyy[i] * fe_0 + tr_x_zzzzz_xyy[i] * pa_z[i];

        tr_x_zzzzzz_xyz[i] = 5.0 * tr_x_zzzz_xyz[i] * fe_0 + tr_x_zzzzz_xy[i] * fe_0 + tr_x_zzzzz_xyz[i] * pa_z[i];

        tr_x_zzzzzz_xzz[i] = 5.0 * tr_x_zzzz_xzz[i] * fe_0 + 2.0 * tr_x_zzzzz_xz[i] * fe_0 + tr_x_zzzzz_xzz[i] * pa_z[i];

        tr_x_zzzzzz_yyy[i] = 5.0 * tr_x_zzzz_yyy[i] * fe_0 + tr_x_zzzzz_yyy[i] * pa_z[i];

        tr_x_zzzzzz_yyz[i] = 5.0 * tr_x_zzzz_yyz[i] * fe_0 + tr_x_zzzzz_yy[i] * fe_0 + tr_x_zzzzz_yyz[i] * pa_z[i];

        tr_x_zzzzzz_yzz[i] = 5.0 * tr_x_zzzz_yzz[i] * fe_0 + 2.0 * tr_x_zzzzz_yz[i] * fe_0 + tr_x_zzzzz_yzz[i] * pa_z[i];

        tr_x_zzzzzz_zzz[i] = 5.0 * tr_x_zzzz_zzz[i] * fe_0 + 3.0 * tr_x_zzzzz_zz[i] * fe_0 + tr_x_zzzzz_zzz[i] * pa_z[i];
    }

    // Set up 280-290 components of targeted buffer : IF

    auto tr_y_xxxxxx_xxx = pbuffer.data(idx_dip_if + 280);

    auto tr_y_xxxxxx_xxy = pbuffer.data(idx_dip_if + 281);

    auto tr_y_xxxxxx_xxz = pbuffer.data(idx_dip_if + 282);

    auto tr_y_xxxxxx_xyy = pbuffer.data(idx_dip_if + 283);

    auto tr_y_xxxxxx_xyz = pbuffer.data(idx_dip_if + 284);

    auto tr_y_xxxxxx_xzz = pbuffer.data(idx_dip_if + 285);

    auto tr_y_xxxxxx_yyy = pbuffer.data(idx_dip_if + 286);

    auto tr_y_xxxxxx_yyz = pbuffer.data(idx_dip_if + 287);

    auto tr_y_xxxxxx_yzz = pbuffer.data(idx_dip_if + 288);

    auto tr_y_xxxxxx_zzz = pbuffer.data(idx_dip_if + 289);

    #pragma omp simd aligned(pa_x, tr_y_xxxx_xxx, tr_y_xxxx_xxy, tr_y_xxxx_xxz, tr_y_xxxx_xyy, tr_y_xxxx_xyz, tr_y_xxxx_xzz, tr_y_xxxx_yyy, tr_y_xxxx_yyz, tr_y_xxxx_yzz, tr_y_xxxx_zzz, tr_y_xxxxx_xx, tr_y_xxxxx_xxx, tr_y_xxxxx_xxy, tr_y_xxxxx_xxz, tr_y_xxxxx_xy, tr_y_xxxxx_xyy, tr_y_xxxxx_xyz, tr_y_xxxxx_xz, tr_y_xxxxx_xzz, tr_y_xxxxx_yy, tr_y_xxxxx_yyy, tr_y_xxxxx_yyz, tr_y_xxxxx_yz, tr_y_xxxxx_yzz, tr_y_xxxxx_zz, tr_y_xxxxx_zzz, tr_y_xxxxxx_xxx, tr_y_xxxxxx_xxy, tr_y_xxxxxx_xxz, tr_y_xxxxxx_xyy, tr_y_xxxxxx_xyz, tr_y_xxxxxx_xzz, tr_y_xxxxxx_yyy, tr_y_xxxxxx_yyz, tr_y_xxxxxx_yzz, tr_y_xxxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxx_xxx[i] = 5.0 * tr_y_xxxx_xxx[i] * fe_0 + 3.0 * tr_y_xxxxx_xx[i] * fe_0 + tr_y_xxxxx_xxx[i] * pa_x[i];

        tr_y_xxxxxx_xxy[i] = 5.0 * tr_y_xxxx_xxy[i] * fe_0 + 2.0 * tr_y_xxxxx_xy[i] * fe_0 + tr_y_xxxxx_xxy[i] * pa_x[i];

        tr_y_xxxxxx_xxz[i] = 5.0 * tr_y_xxxx_xxz[i] * fe_0 + 2.0 * tr_y_xxxxx_xz[i] * fe_0 + tr_y_xxxxx_xxz[i] * pa_x[i];

        tr_y_xxxxxx_xyy[i] = 5.0 * tr_y_xxxx_xyy[i] * fe_0 + tr_y_xxxxx_yy[i] * fe_0 + tr_y_xxxxx_xyy[i] * pa_x[i];

        tr_y_xxxxxx_xyz[i] = 5.0 * tr_y_xxxx_xyz[i] * fe_0 + tr_y_xxxxx_yz[i] * fe_0 + tr_y_xxxxx_xyz[i] * pa_x[i];

        tr_y_xxxxxx_xzz[i] = 5.0 * tr_y_xxxx_xzz[i] * fe_0 + tr_y_xxxxx_zz[i] * fe_0 + tr_y_xxxxx_xzz[i] * pa_x[i];

        tr_y_xxxxxx_yyy[i] = 5.0 * tr_y_xxxx_yyy[i] * fe_0 + tr_y_xxxxx_yyy[i] * pa_x[i];

        tr_y_xxxxxx_yyz[i] = 5.0 * tr_y_xxxx_yyz[i] * fe_0 + tr_y_xxxxx_yyz[i] * pa_x[i];

        tr_y_xxxxxx_yzz[i] = 5.0 * tr_y_xxxx_yzz[i] * fe_0 + tr_y_xxxxx_yzz[i] * pa_x[i];

        tr_y_xxxxxx_zzz[i] = 5.0 * tr_y_xxxx_zzz[i] * fe_0 + tr_y_xxxxx_zzz[i] * pa_x[i];
    }

    // Set up 290-300 components of targeted buffer : IF

    auto tr_y_xxxxxy_xxx = pbuffer.data(idx_dip_if + 290);

    auto tr_y_xxxxxy_xxy = pbuffer.data(idx_dip_if + 291);

    auto tr_y_xxxxxy_xxz = pbuffer.data(idx_dip_if + 292);

    auto tr_y_xxxxxy_xyy = pbuffer.data(idx_dip_if + 293);

    auto tr_y_xxxxxy_xyz = pbuffer.data(idx_dip_if + 294);

    auto tr_y_xxxxxy_xzz = pbuffer.data(idx_dip_if + 295);

    auto tr_y_xxxxxy_yyy = pbuffer.data(idx_dip_if + 296);

    auto tr_y_xxxxxy_yyz = pbuffer.data(idx_dip_if + 297);

    auto tr_y_xxxxxy_yzz = pbuffer.data(idx_dip_if + 298);

    auto tr_y_xxxxxy_zzz = pbuffer.data(idx_dip_if + 299);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xxxxx_xxx, tr_y_xxxxx_xxz, tr_y_xxxxx_xzz, tr_y_xxxxxy_xxx, tr_y_xxxxxy_xxy, tr_y_xxxxxy_xxz, tr_y_xxxxxy_xyy, tr_y_xxxxxy_xyz, tr_y_xxxxxy_xzz, tr_y_xxxxxy_yyy, tr_y_xxxxxy_yyz, tr_y_xxxxxy_yzz, tr_y_xxxxxy_zzz, tr_y_xxxxy_xxy, tr_y_xxxxy_xy, tr_y_xxxxy_xyy, tr_y_xxxxy_xyz, tr_y_xxxxy_yy, tr_y_xxxxy_yyy, tr_y_xxxxy_yyz, tr_y_xxxxy_yz, tr_y_xxxxy_yzz, tr_y_xxxxy_zzz, tr_y_xxxy_xxy, tr_y_xxxy_xyy, tr_y_xxxy_xyz, tr_y_xxxy_yyy, tr_y_xxxy_yyz, tr_y_xxxy_yzz, tr_y_xxxy_zzz, ts_xxxxx_xxx, ts_xxxxx_xxz, ts_xxxxx_xzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxy_xxx[i] = ts_xxxxx_xxx[i] * fe_0 + tr_y_xxxxx_xxx[i] * pa_y[i];

        tr_y_xxxxxy_xxy[i] = 4.0 * tr_y_xxxy_xxy[i] * fe_0 + 2.0 * tr_y_xxxxy_xy[i] * fe_0 + tr_y_xxxxy_xxy[i] * pa_x[i];

        tr_y_xxxxxy_xxz[i] = ts_xxxxx_xxz[i] * fe_0 + tr_y_xxxxx_xxz[i] * pa_y[i];

        tr_y_xxxxxy_xyy[i] = 4.0 * tr_y_xxxy_xyy[i] * fe_0 + tr_y_xxxxy_yy[i] * fe_0 + tr_y_xxxxy_xyy[i] * pa_x[i];

        tr_y_xxxxxy_xyz[i] = 4.0 * tr_y_xxxy_xyz[i] * fe_0 + tr_y_xxxxy_yz[i] * fe_0 + tr_y_xxxxy_xyz[i] * pa_x[i];

        tr_y_xxxxxy_xzz[i] = ts_xxxxx_xzz[i] * fe_0 + tr_y_xxxxx_xzz[i] * pa_y[i];

        tr_y_xxxxxy_yyy[i] = 4.0 * tr_y_xxxy_yyy[i] * fe_0 + tr_y_xxxxy_yyy[i] * pa_x[i];

        tr_y_xxxxxy_yyz[i] = 4.0 * tr_y_xxxy_yyz[i] * fe_0 + tr_y_xxxxy_yyz[i] * pa_x[i];

        tr_y_xxxxxy_yzz[i] = 4.0 * tr_y_xxxy_yzz[i] * fe_0 + tr_y_xxxxy_yzz[i] * pa_x[i];

        tr_y_xxxxxy_zzz[i] = 4.0 * tr_y_xxxy_zzz[i] * fe_0 + tr_y_xxxxy_zzz[i] * pa_x[i];
    }

    // Set up 300-310 components of targeted buffer : IF

    auto tr_y_xxxxxz_xxx = pbuffer.data(idx_dip_if + 300);

    auto tr_y_xxxxxz_xxy = pbuffer.data(idx_dip_if + 301);

    auto tr_y_xxxxxz_xxz = pbuffer.data(idx_dip_if + 302);

    auto tr_y_xxxxxz_xyy = pbuffer.data(idx_dip_if + 303);

    auto tr_y_xxxxxz_xyz = pbuffer.data(idx_dip_if + 304);

    auto tr_y_xxxxxz_xzz = pbuffer.data(idx_dip_if + 305);

    auto tr_y_xxxxxz_yyy = pbuffer.data(idx_dip_if + 306);

    auto tr_y_xxxxxz_yyz = pbuffer.data(idx_dip_if + 307);

    auto tr_y_xxxxxz_yzz = pbuffer.data(idx_dip_if + 308);

    auto tr_y_xxxxxz_zzz = pbuffer.data(idx_dip_if + 309);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxxx_xx, tr_y_xxxxx_xxx, tr_y_xxxxx_xxy, tr_y_xxxxx_xxz, tr_y_xxxxx_xy, tr_y_xxxxx_xyy, tr_y_xxxxx_xyz, tr_y_xxxxx_xz, tr_y_xxxxx_xzz, tr_y_xxxxx_yyy, tr_y_xxxxxz_xxx, tr_y_xxxxxz_xxy, tr_y_xxxxxz_xxz, tr_y_xxxxxz_xyy, tr_y_xxxxxz_xyz, tr_y_xxxxxz_xzz, tr_y_xxxxxz_yyy, tr_y_xxxxxz_yyz, tr_y_xxxxxz_yzz, tr_y_xxxxxz_zzz, tr_y_xxxxz_yyz, tr_y_xxxxz_yzz, tr_y_xxxxz_zzz, tr_y_xxxz_yyz, tr_y_xxxz_yzz, tr_y_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxz_xxx[i] = tr_y_xxxxx_xxx[i] * pa_z[i];

        tr_y_xxxxxz_xxy[i] = tr_y_xxxxx_xxy[i] * pa_z[i];

        tr_y_xxxxxz_xxz[i] = tr_y_xxxxx_xx[i] * fe_0 + tr_y_xxxxx_xxz[i] * pa_z[i];

        tr_y_xxxxxz_xyy[i] = tr_y_xxxxx_xyy[i] * pa_z[i];

        tr_y_xxxxxz_xyz[i] = tr_y_xxxxx_xy[i] * fe_0 + tr_y_xxxxx_xyz[i] * pa_z[i];

        tr_y_xxxxxz_xzz[i] = 2.0 * tr_y_xxxxx_xz[i] * fe_0 + tr_y_xxxxx_xzz[i] * pa_z[i];

        tr_y_xxxxxz_yyy[i] = tr_y_xxxxx_yyy[i] * pa_z[i];

        tr_y_xxxxxz_yyz[i] = 4.0 * tr_y_xxxz_yyz[i] * fe_0 + tr_y_xxxxz_yyz[i] * pa_x[i];

        tr_y_xxxxxz_yzz[i] = 4.0 * tr_y_xxxz_yzz[i] * fe_0 + tr_y_xxxxz_yzz[i] * pa_x[i];

        tr_y_xxxxxz_zzz[i] = 4.0 * tr_y_xxxz_zzz[i] * fe_0 + tr_y_xxxxz_zzz[i] * pa_x[i];
    }

    // Set up 310-320 components of targeted buffer : IF

    auto tr_y_xxxxyy_xxx = pbuffer.data(idx_dip_if + 310);

    auto tr_y_xxxxyy_xxy = pbuffer.data(idx_dip_if + 311);

    auto tr_y_xxxxyy_xxz = pbuffer.data(idx_dip_if + 312);

    auto tr_y_xxxxyy_xyy = pbuffer.data(idx_dip_if + 313);

    auto tr_y_xxxxyy_xyz = pbuffer.data(idx_dip_if + 314);

    auto tr_y_xxxxyy_xzz = pbuffer.data(idx_dip_if + 315);

    auto tr_y_xxxxyy_yyy = pbuffer.data(idx_dip_if + 316);

    auto tr_y_xxxxyy_yyz = pbuffer.data(idx_dip_if + 317);

    auto tr_y_xxxxyy_yzz = pbuffer.data(idx_dip_if + 318);

    auto tr_y_xxxxyy_zzz = pbuffer.data(idx_dip_if + 319);

    #pragma omp simd aligned(pa_x, tr_y_xxxxyy_xxx, tr_y_xxxxyy_xxy, tr_y_xxxxyy_xxz, tr_y_xxxxyy_xyy, tr_y_xxxxyy_xyz, tr_y_xxxxyy_xzz, tr_y_xxxxyy_yyy, tr_y_xxxxyy_yyz, tr_y_xxxxyy_yzz, tr_y_xxxxyy_zzz, tr_y_xxxyy_xx, tr_y_xxxyy_xxx, tr_y_xxxyy_xxy, tr_y_xxxyy_xxz, tr_y_xxxyy_xy, tr_y_xxxyy_xyy, tr_y_xxxyy_xyz, tr_y_xxxyy_xz, tr_y_xxxyy_xzz, tr_y_xxxyy_yy, tr_y_xxxyy_yyy, tr_y_xxxyy_yyz, tr_y_xxxyy_yz, tr_y_xxxyy_yzz, tr_y_xxxyy_zz, tr_y_xxxyy_zzz, tr_y_xxyy_xxx, tr_y_xxyy_xxy, tr_y_xxyy_xxz, tr_y_xxyy_xyy, tr_y_xxyy_xyz, tr_y_xxyy_xzz, tr_y_xxyy_yyy, tr_y_xxyy_yyz, tr_y_xxyy_yzz, tr_y_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyy_xxx[i] = 3.0 * tr_y_xxyy_xxx[i] * fe_0 + 3.0 * tr_y_xxxyy_xx[i] * fe_0 + tr_y_xxxyy_xxx[i] * pa_x[i];

        tr_y_xxxxyy_xxy[i] = 3.0 * tr_y_xxyy_xxy[i] * fe_0 + 2.0 * tr_y_xxxyy_xy[i] * fe_0 + tr_y_xxxyy_xxy[i] * pa_x[i];

        tr_y_xxxxyy_xxz[i] = 3.0 * tr_y_xxyy_xxz[i] * fe_0 + 2.0 * tr_y_xxxyy_xz[i] * fe_0 + tr_y_xxxyy_xxz[i] * pa_x[i];

        tr_y_xxxxyy_xyy[i] = 3.0 * tr_y_xxyy_xyy[i] * fe_0 + tr_y_xxxyy_yy[i] * fe_0 + tr_y_xxxyy_xyy[i] * pa_x[i];

        tr_y_xxxxyy_xyz[i] = 3.0 * tr_y_xxyy_xyz[i] * fe_0 + tr_y_xxxyy_yz[i] * fe_0 + tr_y_xxxyy_xyz[i] * pa_x[i];

        tr_y_xxxxyy_xzz[i] = 3.0 * tr_y_xxyy_xzz[i] * fe_0 + tr_y_xxxyy_zz[i] * fe_0 + tr_y_xxxyy_xzz[i] * pa_x[i];

        tr_y_xxxxyy_yyy[i] = 3.0 * tr_y_xxyy_yyy[i] * fe_0 + tr_y_xxxyy_yyy[i] * pa_x[i];

        tr_y_xxxxyy_yyz[i] = 3.0 * tr_y_xxyy_yyz[i] * fe_0 + tr_y_xxxyy_yyz[i] * pa_x[i];

        tr_y_xxxxyy_yzz[i] = 3.0 * tr_y_xxyy_yzz[i] * fe_0 + tr_y_xxxyy_yzz[i] * pa_x[i];

        tr_y_xxxxyy_zzz[i] = 3.0 * tr_y_xxyy_zzz[i] * fe_0 + tr_y_xxxyy_zzz[i] * pa_x[i];
    }

    // Set up 320-330 components of targeted buffer : IF

    auto tr_y_xxxxyz_xxx = pbuffer.data(idx_dip_if + 320);

    auto tr_y_xxxxyz_xxy = pbuffer.data(idx_dip_if + 321);

    auto tr_y_xxxxyz_xxz = pbuffer.data(idx_dip_if + 322);

    auto tr_y_xxxxyz_xyy = pbuffer.data(idx_dip_if + 323);

    auto tr_y_xxxxyz_xyz = pbuffer.data(idx_dip_if + 324);

    auto tr_y_xxxxyz_xzz = pbuffer.data(idx_dip_if + 325);

    auto tr_y_xxxxyz_yyy = pbuffer.data(idx_dip_if + 326);

    auto tr_y_xxxxyz_yyz = pbuffer.data(idx_dip_if + 327);

    auto tr_y_xxxxyz_yzz = pbuffer.data(idx_dip_if + 328);

    auto tr_y_xxxxyz_zzz = pbuffer.data(idx_dip_if + 329);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxxy_xxx, tr_y_xxxxy_xxy, tr_y_xxxxy_xy, tr_y_xxxxy_xyy, tr_y_xxxxy_xyz, tr_y_xxxxy_yyy, tr_y_xxxxyz_xxx, tr_y_xxxxyz_xxy, tr_y_xxxxyz_xxz, tr_y_xxxxyz_xyy, tr_y_xxxxyz_xyz, tr_y_xxxxyz_xzz, tr_y_xxxxyz_yyy, tr_y_xxxxyz_yyz, tr_y_xxxxyz_yzz, tr_y_xxxxyz_zzz, tr_y_xxxxz_xxz, tr_y_xxxxz_xzz, tr_y_xxxyz_yyz, tr_y_xxxyz_yzz, tr_y_xxxyz_zzz, tr_y_xxyz_yyz, tr_y_xxyz_yzz, tr_y_xxyz_zzz, ts_xxxxz_xxz, ts_xxxxz_xzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyz_xxx[i] = tr_y_xxxxy_xxx[i] * pa_z[i];

        tr_y_xxxxyz_xxy[i] = tr_y_xxxxy_xxy[i] * pa_z[i];

        tr_y_xxxxyz_xxz[i] = ts_xxxxz_xxz[i] * fe_0 + tr_y_xxxxz_xxz[i] * pa_y[i];

        tr_y_xxxxyz_xyy[i] = tr_y_xxxxy_xyy[i] * pa_z[i];

        tr_y_xxxxyz_xyz[i] = tr_y_xxxxy_xy[i] * fe_0 + tr_y_xxxxy_xyz[i] * pa_z[i];

        tr_y_xxxxyz_xzz[i] = ts_xxxxz_xzz[i] * fe_0 + tr_y_xxxxz_xzz[i] * pa_y[i];

        tr_y_xxxxyz_yyy[i] = tr_y_xxxxy_yyy[i] * pa_z[i];

        tr_y_xxxxyz_yyz[i] = 3.0 * tr_y_xxyz_yyz[i] * fe_0 + tr_y_xxxyz_yyz[i] * pa_x[i];

        tr_y_xxxxyz_yzz[i] = 3.0 * tr_y_xxyz_yzz[i] * fe_0 + tr_y_xxxyz_yzz[i] * pa_x[i];

        tr_y_xxxxyz_zzz[i] = 3.0 * tr_y_xxyz_zzz[i] * fe_0 + tr_y_xxxyz_zzz[i] * pa_x[i];
    }

    // Set up 330-340 components of targeted buffer : IF

    auto tr_y_xxxxzz_xxx = pbuffer.data(idx_dip_if + 330);

    auto tr_y_xxxxzz_xxy = pbuffer.data(idx_dip_if + 331);

    auto tr_y_xxxxzz_xxz = pbuffer.data(idx_dip_if + 332);

    auto tr_y_xxxxzz_xyy = pbuffer.data(idx_dip_if + 333);

    auto tr_y_xxxxzz_xyz = pbuffer.data(idx_dip_if + 334);

    auto tr_y_xxxxzz_xzz = pbuffer.data(idx_dip_if + 335);

    auto tr_y_xxxxzz_yyy = pbuffer.data(idx_dip_if + 336);

    auto tr_y_xxxxzz_yyz = pbuffer.data(idx_dip_if + 337);

    auto tr_y_xxxxzz_yzz = pbuffer.data(idx_dip_if + 338);

    auto tr_y_xxxxzz_zzz = pbuffer.data(idx_dip_if + 339);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxx_xxx, tr_y_xxxx_xxy, tr_y_xxxx_xyy, tr_y_xxxxz_xxx, tr_y_xxxxz_xxy, tr_y_xxxxz_xyy, tr_y_xxxxzz_xxx, tr_y_xxxxzz_xxy, tr_y_xxxxzz_xxz, tr_y_xxxxzz_xyy, tr_y_xxxxzz_xyz, tr_y_xxxxzz_xzz, tr_y_xxxxzz_yyy, tr_y_xxxxzz_yyz, tr_y_xxxxzz_yzz, tr_y_xxxxzz_zzz, tr_y_xxxzz_xxz, tr_y_xxxzz_xyz, tr_y_xxxzz_xz, tr_y_xxxzz_xzz, tr_y_xxxzz_yyy, tr_y_xxxzz_yyz, tr_y_xxxzz_yz, tr_y_xxxzz_yzz, tr_y_xxxzz_zz, tr_y_xxxzz_zzz, tr_y_xxzz_xxz, tr_y_xxzz_xyz, tr_y_xxzz_xzz, tr_y_xxzz_yyy, tr_y_xxzz_yyz, tr_y_xxzz_yzz, tr_y_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxzz_xxx[i] = tr_y_xxxx_xxx[i] * fe_0 + tr_y_xxxxz_xxx[i] * pa_z[i];

        tr_y_xxxxzz_xxy[i] = tr_y_xxxx_xxy[i] * fe_0 + tr_y_xxxxz_xxy[i] * pa_z[i];

        tr_y_xxxxzz_xxz[i] = 3.0 * tr_y_xxzz_xxz[i] * fe_0 + 2.0 * tr_y_xxxzz_xz[i] * fe_0 + tr_y_xxxzz_xxz[i] * pa_x[i];

        tr_y_xxxxzz_xyy[i] = tr_y_xxxx_xyy[i] * fe_0 + tr_y_xxxxz_xyy[i] * pa_z[i];

        tr_y_xxxxzz_xyz[i] = 3.0 * tr_y_xxzz_xyz[i] * fe_0 + tr_y_xxxzz_yz[i] * fe_0 + tr_y_xxxzz_xyz[i] * pa_x[i];

        tr_y_xxxxzz_xzz[i] = 3.0 * tr_y_xxzz_xzz[i] * fe_0 + tr_y_xxxzz_zz[i] * fe_0 + tr_y_xxxzz_xzz[i] * pa_x[i];

        tr_y_xxxxzz_yyy[i] = 3.0 * tr_y_xxzz_yyy[i] * fe_0 + tr_y_xxxzz_yyy[i] * pa_x[i];

        tr_y_xxxxzz_yyz[i] = 3.0 * tr_y_xxzz_yyz[i] * fe_0 + tr_y_xxxzz_yyz[i] * pa_x[i];

        tr_y_xxxxzz_yzz[i] = 3.0 * tr_y_xxzz_yzz[i] * fe_0 + tr_y_xxxzz_yzz[i] * pa_x[i];

        tr_y_xxxxzz_zzz[i] = 3.0 * tr_y_xxzz_zzz[i] * fe_0 + tr_y_xxxzz_zzz[i] * pa_x[i];
    }

    // Set up 340-350 components of targeted buffer : IF

    auto tr_y_xxxyyy_xxx = pbuffer.data(idx_dip_if + 340);

    auto tr_y_xxxyyy_xxy = pbuffer.data(idx_dip_if + 341);

    auto tr_y_xxxyyy_xxz = pbuffer.data(idx_dip_if + 342);

    auto tr_y_xxxyyy_xyy = pbuffer.data(idx_dip_if + 343);

    auto tr_y_xxxyyy_xyz = pbuffer.data(idx_dip_if + 344);

    auto tr_y_xxxyyy_xzz = pbuffer.data(idx_dip_if + 345);

    auto tr_y_xxxyyy_yyy = pbuffer.data(idx_dip_if + 346);

    auto tr_y_xxxyyy_yyz = pbuffer.data(idx_dip_if + 347);

    auto tr_y_xxxyyy_yzz = pbuffer.data(idx_dip_if + 348);

    auto tr_y_xxxyyy_zzz = pbuffer.data(idx_dip_if + 349);

    #pragma omp simd aligned(pa_x, tr_y_xxxyyy_xxx, tr_y_xxxyyy_xxy, tr_y_xxxyyy_xxz, tr_y_xxxyyy_xyy, tr_y_xxxyyy_xyz, tr_y_xxxyyy_xzz, tr_y_xxxyyy_yyy, tr_y_xxxyyy_yyz, tr_y_xxxyyy_yzz, tr_y_xxxyyy_zzz, tr_y_xxyyy_xx, tr_y_xxyyy_xxx, tr_y_xxyyy_xxy, tr_y_xxyyy_xxz, tr_y_xxyyy_xy, tr_y_xxyyy_xyy, tr_y_xxyyy_xyz, tr_y_xxyyy_xz, tr_y_xxyyy_xzz, tr_y_xxyyy_yy, tr_y_xxyyy_yyy, tr_y_xxyyy_yyz, tr_y_xxyyy_yz, tr_y_xxyyy_yzz, tr_y_xxyyy_zz, tr_y_xxyyy_zzz, tr_y_xyyy_xxx, tr_y_xyyy_xxy, tr_y_xyyy_xxz, tr_y_xyyy_xyy, tr_y_xyyy_xyz, tr_y_xyyy_xzz, tr_y_xyyy_yyy, tr_y_xyyy_yyz, tr_y_xyyy_yzz, tr_y_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyy_xxx[i] = 2.0 * tr_y_xyyy_xxx[i] * fe_0 + 3.0 * tr_y_xxyyy_xx[i] * fe_0 + tr_y_xxyyy_xxx[i] * pa_x[i];

        tr_y_xxxyyy_xxy[i] = 2.0 * tr_y_xyyy_xxy[i] * fe_0 + 2.0 * tr_y_xxyyy_xy[i] * fe_0 + tr_y_xxyyy_xxy[i] * pa_x[i];

        tr_y_xxxyyy_xxz[i] = 2.0 * tr_y_xyyy_xxz[i] * fe_0 + 2.0 * tr_y_xxyyy_xz[i] * fe_0 + tr_y_xxyyy_xxz[i] * pa_x[i];

        tr_y_xxxyyy_xyy[i] = 2.0 * tr_y_xyyy_xyy[i] * fe_0 + tr_y_xxyyy_yy[i] * fe_0 + tr_y_xxyyy_xyy[i] * pa_x[i];

        tr_y_xxxyyy_xyz[i] = 2.0 * tr_y_xyyy_xyz[i] * fe_0 + tr_y_xxyyy_yz[i] * fe_0 + tr_y_xxyyy_xyz[i] * pa_x[i];

        tr_y_xxxyyy_xzz[i] = 2.0 * tr_y_xyyy_xzz[i] * fe_0 + tr_y_xxyyy_zz[i] * fe_0 + tr_y_xxyyy_xzz[i] * pa_x[i];

        tr_y_xxxyyy_yyy[i] = 2.0 * tr_y_xyyy_yyy[i] * fe_0 + tr_y_xxyyy_yyy[i] * pa_x[i];

        tr_y_xxxyyy_yyz[i] = 2.0 * tr_y_xyyy_yyz[i] * fe_0 + tr_y_xxyyy_yyz[i] * pa_x[i];

        tr_y_xxxyyy_yzz[i] = 2.0 * tr_y_xyyy_yzz[i] * fe_0 + tr_y_xxyyy_yzz[i] * pa_x[i];

        tr_y_xxxyyy_zzz[i] = 2.0 * tr_y_xyyy_zzz[i] * fe_0 + tr_y_xxyyy_zzz[i] * pa_x[i];
    }

    // Set up 350-360 components of targeted buffer : IF

    auto tr_y_xxxyyz_xxx = pbuffer.data(idx_dip_if + 350);

    auto tr_y_xxxyyz_xxy = pbuffer.data(idx_dip_if + 351);

    auto tr_y_xxxyyz_xxz = pbuffer.data(idx_dip_if + 352);

    auto tr_y_xxxyyz_xyy = pbuffer.data(idx_dip_if + 353);

    auto tr_y_xxxyyz_xyz = pbuffer.data(idx_dip_if + 354);

    auto tr_y_xxxyyz_xzz = pbuffer.data(idx_dip_if + 355);

    auto tr_y_xxxyyz_yyy = pbuffer.data(idx_dip_if + 356);

    auto tr_y_xxxyyz_yyz = pbuffer.data(idx_dip_if + 357);

    auto tr_y_xxxyyz_yzz = pbuffer.data(idx_dip_if + 358);

    auto tr_y_xxxyyz_zzz = pbuffer.data(idx_dip_if + 359);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxyy_xx, tr_y_xxxyy_xxx, tr_y_xxxyy_xxy, tr_y_xxxyy_xxz, tr_y_xxxyy_xy, tr_y_xxxyy_xyy, tr_y_xxxyy_xyz, tr_y_xxxyy_xz, tr_y_xxxyy_xzz, tr_y_xxxyy_yyy, tr_y_xxxyyz_xxx, tr_y_xxxyyz_xxy, tr_y_xxxyyz_xxz, tr_y_xxxyyz_xyy, tr_y_xxxyyz_xyz, tr_y_xxxyyz_xzz, tr_y_xxxyyz_yyy, tr_y_xxxyyz_yyz, tr_y_xxxyyz_yzz, tr_y_xxxyyz_zzz, tr_y_xxyyz_yyz, tr_y_xxyyz_yzz, tr_y_xxyyz_zzz, tr_y_xyyz_yyz, tr_y_xyyz_yzz, tr_y_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyz_xxx[i] = tr_y_xxxyy_xxx[i] * pa_z[i];

        tr_y_xxxyyz_xxy[i] = tr_y_xxxyy_xxy[i] * pa_z[i];

        tr_y_xxxyyz_xxz[i] = tr_y_xxxyy_xx[i] * fe_0 + tr_y_xxxyy_xxz[i] * pa_z[i];

        tr_y_xxxyyz_xyy[i] = tr_y_xxxyy_xyy[i] * pa_z[i];

        tr_y_xxxyyz_xyz[i] = tr_y_xxxyy_xy[i] * fe_0 + tr_y_xxxyy_xyz[i] * pa_z[i];

        tr_y_xxxyyz_xzz[i] = 2.0 * tr_y_xxxyy_xz[i] * fe_0 + tr_y_xxxyy_xzz[i] * pa_z[i];

        tr_y_xxxyyz_yyy[i] = tr_y_xxxyy_yyy[i] * pa_z[i];

        tr_y_xxxyyz_yyz[i] = 2.0 * tr_y_xyyz_yyz[i] * fe_0 + tr_y_xxyyz_yyz[i] * pa_x[i];

        tr_y_xxxyyz_yzz[i] = 2.0 * tr_y_xyyz_yzz[i] * fe_0 + tr_y_xxyyz_yzz[i] * pa_x[i];

        tr_y_xxxyyz_zzz[i] = 2.0 * tr_y_xyyz_zzz[i] * fe_0 + tr_y_xxyyz_zzz[i] * pa_x[i];
    }

    // Set up 360-370 components of targeted buffer : IF

    auto tr_y_xxxyzz_xxx = pbuffer.data(idx_dip_if + 360);

    auto tr_y_xxxyzz_xxy = pbuffer.data(idx_dip_if + 361);

    auto tr_y_xxxyzz_xxz = pbuffer.data(idx_dip_if + 362);

    auto tr_y_xxxyzz_xyy = pbuffer.data(idx_dip_if + 363);

    auto tr_y_xxxyzz_xyz = pbuffer.data(idx_dip_if + 364);

    auto tr_y_xxxyzz_xzz = pbuffer.data(idx_dip_if + 365);

    auto tr_y_xxxyzz_yyy = pbuffer.data(idx_dip_if + 366);

    auto tr_y_xxxyzz_yyz = pbuffer.data(idx_dip_if + 367);

    auto tr_y_xxxyzz_yzz = pbuffer.data(idx_dip_if + 368);

    auto tr_y_xxxyzz_zzz = pbuffer.data(idx_dip_if + 369);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxy_xxy, tr_y_xxxy_xyy, tr_y_xxxyz_xxy, tr_y_xxxyz_xyy, tr_y_xxxyzz_xxx, tr_y_xxxyzz_xxy, tr_y_xxxyzz_xxz, tr_y_xxxyzz_xyy, tr_y_xxxyzz_xyz, tr_y_xxxyzz_xzz, tr_y_xxxyzz_yyy, tr_y_xxxyzz_yyz, tr_y_xxxyzz_yzz, tr_y_xxxyzz_zzz, tr_y_xxxzz_xxx, tr_y_xxxzz_xxz, tr_y_xxxzz_xzz, tr_y_xxyzz_xyz, tr_y_xxyzz_yyy, tr_y_xxyzz_yyz, tr_y_xxyzz_yz, tr_y_xxyzz_yzz, tr_y_xxyzz_zzz, tr_y_xyzz_xyz, tr_y_xyzz_yyy, tr_y_xyzz_yyz, tr_y_xyzz_yzz, tr_y_xyzz_zzz, ts_xxxzz_xxx, ts_xxxzz_xxz, ts_xxxzz_xzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyzz_xxx[i] = ts_xxxzz_xxx[i] * fe_0 + tr_y_xxxzz_xxx[i] * pa_y[i];

        tr_y_xxxyzz_xxy[i] = tr_y_xxxy_xxy[i] * fe_0 + tr_y_xxxyz_xxy[i] * pa_z[i];

        tr_y_xxxyzz_xxz[i] = ts_xxxzz_xxz[i] * fe_0 + tr_y_xxxzz_xxz[i] * pa_y[i];

        tr_y_xxxyzz_xyy[i] = tr_y_xxxy_xyy[i] * fe_0 + tr_y_xxxyz_xyy[i] * pa_z[i];

        tr_y_xxxyzz_xyz[i] = 2.0 * tr_y_xyzz_xyz[i] * fe_0 + tr_y_xxyzz_yz[i] * fe_0 + tr_y_xxyzz_xyz[i] * pa_x[i];

        tr_y_xxxyzz_xzz[i] = ts_xxxzz_xzz[i] * fe_0 + tr_y_xxxzz_xzz[i] * pa_y[i];

        tr_y_xxxyzz_yyy[i] = 2.0 * tr_y_xyzz_yyy[i] * fe_0 + tr_y_xxyzz_yyy[i] * pa_x[i];

        tr_y_xxxyzz_yyz[i] = 2.0 * tr_y_xyzz_yyz[i] * fe_0 + tr_y_xxyzz_yyz[i] * pa_x[i];

        tr_y_xxxyzz_yzz[i] = 2.0 * tr_y_xyzz_yzz[i] * fe_0 + tr_y_xxyzz_yzz[i] * pa_x[i];

        tr_y_xxxyzz_zzz[i] = 2.0 * tr_y_xyzz_zzz[i] * fe_0 + tr_y_xxyzz_zzz[i] * pa_x[i];
    }

    // Set up 370-380 components of targeted buffer : IF

    auto tr_y_xxxzzz_xxx = pbuffer.data(idx_dip_if + 370);

    auto tr_y_xxxzzz_xxy = pbuffer.data(idx_dip_if + 371);

    auto tr_y_xxxzzz_xxz = pbuffer.data(idx_dip_if + 372);

    auto tr_y_xxxzzz_xyy = pbuffer.data(idx_dip_if + 373);

    auto tr_y_xxxzzz_xyz = pbuffer.data(idx_dip_if + 374);

    auto tr_y_xxxzzz_xzz = pbuffer.data(idx_dip_if + 375);

    auto tr_y_xxxzzz_yyy = pbuffer.data(idx_dip_if + 376);

    auto tr_y_xxxzzz_yyz = pbuffer.data(idx_dip_if + 377);

    auto tr_y_xxxzzz_yzz = pbuffer.data(idx_dip_if + 378);

    auto tr_y_xxxzzz_zzz = pbuffer.data(idx_dip_if + 379);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxz_xxx, tr_y_xxxz_xxy, tr_y_xxxz_xyy, tr_y_xxxzz_xxx, tr_y_xxxzz_xxy, tr_y_xxxzz_xyy, tr_y_xxxzzz_xxx, tr_y_xxxzzz_xxy, tr_y_xxxzzz_xxz, tr_y_xxxzzz_xyy, tr_y_xxxzzz_xyz, tr_y_xxxzzz_xzz, tr_y_xxxzzz_yyy, tr_y_xxxzzz_yyz, tr_y_xxxzzz_yzz, tr_y_xxxzzz_zzz, tr_y_xxzzz_xxz, tr_y_xxzzz_xyz, tr_y_xxzzz_xz, tr_y_xxzzz_xzz, tr_y_xxzzz_yyy, tr_y_xxzzz_yyz, tr_y_xxzzz_yz, tr_y_xxzzz_yzz, tr_y_xxzzz_zz, tr_y_xxzzz_zzz, tr_y_xzzz_xxz, tr_y_xzzz_xyz, tr_y_xzzz_xzz, tr_y_xzzz_yyy, tr_y_xzzz_yyz, tr_y_xzzz_yzz, tr_y_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzzz_xxx[i] = 2.0 * tr_y_xxxz_xxx[i] * fe_0 + tr_y_xxxzz_xxx[i] * pa_z[i];

        tr_y_xxxzzz_xxy[i] = 2.0 * tr_y_xxxz_xxy[i] * fe_0 + tr_y_xxxzz_xxy[i] * pa_z[i];

        tr_y_xxxzzz_xxz[i] = 2.0 * tr_y_xzzz_xxz[i] * fe_0 + 2.0 * tr_y_xxzzz_xz[i] * fe_0 + tr_y_xxzzz_xxz[i] * pa_x[i];

        tr_y_xxxzzz_xyy[i] = 2.0 * tr_y_xxxz_xyy[i] * fe_0 + tr_y_xxxzz_xyy[i] * pa_z[i];

        tr_y_xxxzzz_xyz[i] = 2.0 * tr_y_xzzz_xyz[i] * fe_0 + tr_y_xxzzz_yz[i] * fe_0 + tr_y_xxzzz_xyz[i] * pa_x[i];

        tr_y_xxxzzz_xzz[i] = 2.0 * tr_y_xzzz_xzz[i] * fe_0 + tr_y_xxzzz_zz[i] * fe_0 + tr_y_xxzzz_xzz[i] * pa_x[i];

        tr_y_xxxzzz_yyy[i] = 2.0 * tr_y_xzzz_yyy[i] * fe_0 + tr_y_xxzzz_yyy[i] * pa_x[i];

        tr_y_xxxzzz_yyz[i] = 2.0 * tr_y_xzzz_yyz[i] * fe_0 + tr_y_xxzzz_yyz[i] * pa_x[i];

        tr_y_xxxzzz_yzz[i] = 2.0 * tr_y_xzzz_yzz[i] * fe_0 + tr_y_xxzzz_yzz[i] * pa_x[i];

        tr_y_xxxzzz_zzz[i] = 2.0 * tr_y_xzzz_zzz[i] * fe_0 + tr_y_xxzzz_zzz[i] * pa_x[i];
    }

    // Set up 380-390 components of targeted buffer : IF

    auto tr_y_xxyyyy_xxx = pbuffer.data(idx_dip_if + 380);

    auto tr_y_xxyyyy_xxy = pbuffer.data(idx_dip_if + 381);

    auto tr_y_xxyyyy_xxz = pbuffer.data(idx_dip_if + 382);

    auto tr_y_xxyyyy_xyy = pbuffer.data(idx_dip_if + 383);

    auto tr_y_xxyyyy_xyz = pbuffer.data(idx_dip_if + 384);

    auto tr_y_xxyyyy_xzz = pbuffer.data(idx_dip_if + 385);

    auto tr_y_xxyyyy_yyy = pbuffer.data(idx_dip_if + 386);

    auto tr_y_xxyyyy_yyz = pbuffer.data(idx_dip_if + 387);

    auto tr_y_xxyyyy_yzz = pbuffer.data(idx_dip_if + 388);

    auto tr_y_xxyyyy_zzz = pbuffer.data(idx_dip_if + 389);

    #pragma omp simd aligned(pa_x, tr_y_xxyyyy_xxx, tr_y_xxyyyy_xxy, tr_y_xxyyyy_xxz, tr_y_xxyyyy_xyy, tr_y_xxyyyy_xyz, tr_y_xxyyyy_xzz, tr_y_xxyyyy_yyy, tr_y_xxyyyy_yyz, tr_y_xxyyyy_yzz, tr_y_xxyyyy_zzz, tr_y_xyyyy_xx, tr_y_xyyyy_xxx, tr_y_xyyyy_xxy, tr_y_xyyyy_xxz, tr_y_xyyyy_xy, tr_y_xyyyy_xyy, tr_y_xyyyy_xyz, tr_y_xyyyy_xz, tr_y_xyyyy_xzz, tr_y_xyyyy_yy, tr_y_xyyyy_yyy, tr_y_xyyyy_yyz, tr_y_xyyyy_yz, tr_y_xyyyy_yzz, tr_y_xyyyy_zz, tr_y_xyyyy_zzz, tr_y_yyyy_xxx, tr_y_yyyy_xxy, tr_y_yyyy_xxz, tr_y_yyyy_xyy, tr_y_yyyy_xyz, tr_y_yyyy_xzz, tr_y_yyyy_yyy, tr_y_yyyy_yyz, tr_y_yyyy_yzz, tr_y_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyy_xxx[i] = tr_y_yyyy_xxx[i] * fe_0 + 3.0 * tr_y_xyyyy_xx[i] * fe_0 + tr_y_xyyyy_xxx[i] * pa_x[i];

        tr_y_xxyyyy_xxy[i] = tr_y_yyyy_xxy[i] * fe_0 + 2.0 * tr_y_xyyyy_xy[i] * fe_0 + tr_y_xyyyy_xxy[i] * pa_x[i];

        tr_y_xxyyyy_xxz[i] = tr_y_yyyy_xxz[i] * fe_0 + 2.0 * tr_y_xyyyy_xz[i] * fe_0 + tr_y_xyyyy_xxz[i] * pa_x[i];

        tr_y_xxyyyy_xyy[i] = tr_y_yyyy_xyy[i] * fe_0 + tr_y_xyyyy_yy[i] * fe_0 + tr_y_xyyyy_xyy[i] * pa_x[i];

        tr_y_xxyyyy_xyz[i] = tr_y_yyyy_xyz[i] * fe_0 + tr_y_xyyyy_yz[i] * fe_0 + tr_y_xyyyy_xyz[i] * pa_x[i];

        tr_y_xxyyyy_xzz[i] = tr_y_yyyy_xzz[i] * fe_0 + tr_y_xyyyy_zz[i] * fe_0 + tr_y_xyyyy_xzz[i] * pa_x[i];

        tr_y_xxyyyy_yyy[i] = tr_y_yyyy_yyy[i] * fe_0 + tr_y_xyyyy_yyy[i] * pa_x[i];

        tr_y_xxyyyy_yyz[i] = tr_y_yyyy_yyz[i] * fe_0 + tr_y_xyyyy_yyz[i] * pa_x[i];

        tr_y_xxyyyy_yzz[i] = tr_y_yyyy_yzz[i] * fe_0 + tr_y_xyyyy_yzz[i] * pa_x[i];

        tr_y_xxyyyy_zzz[i] = tr_y_yyyy_zzz[i] * fe_0 + tr_y_xyyyy_zzz[i] * pa_x[i];
    }

    // Set up 390-400 components of targeted buffer : IF

    auto tr_y_xxyyyz_xxx = pbuffer.data(idx_dip_if + 390);

    auto tr_y_xxyyyz_xxy = pbuffer.data(idx_dip_if + 391);

    auto tr_y_xxyyyz_xxz = pbuffer.data(idx_dip_if + 392);

    auto tr_y_xxyyyz_xyy = pbuffer.data(idx_dip_if + 393);

    auto tr_y_xxyyyz_xyz = pbuffer.data(idx_dip_if + 394);

    auto tr_y_xxyyyz_xzz = pbuffer.data(idx_dip_if + 395);

    auto tr_y_xxyyyz_yyy = pbuffer.data(idx_dip_if + 396);

    auto tr_y_xxyyyz_yyz = pbuffer.data(idx_dip_if + 397);

    auto tr_y_xxyyyz_yzz = pbuffer.data(idx_dip_if + 398);

    auto tr_y_xxyyyz_zzz = pbuffer.data(idx_dip_if + 399);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyyy_xx, tr_y_xxyyy_xxx, tr_y_xxyyy_xxy, tr_y_xxyyy_xxz, tr_y_xxyyy_xy, tr_y_xxyyy_xyy, tr_y_xxyyy_xyz, tr_y_xxyyy_xz, tr_y_xxyyy_xzz, tr_y_xxyyy_yyy, tr_y_xxyyyz_xxx, tr_y_xxyyyz_xxy, tr_y_xxyyyz_xxz, tr_y_xxyyyz_xyy, tr_y_xxyyyz_xyz, tr_y_xxyyyz_xzz, tr_y_xxyyyz_yyy, tr_y_xxyyyz_yyz, tr_y_xxyyyz_yzz, tr_y_xxyyyz_zzz, tr_y_xyyyz_yyz, tr_y_xyyyz_yzz, tr_y_xyyyz_zzz, tr_y_yyyz_yyz, tr_y_yyyz_yzz, tr_y_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyz_xxx[i] = tr_y_xxyyy_xxx[i] * pa_z[i];

        tr_y_xxyyyz_xxy[i] = tr_y_xxyyy_xxy[i] * pa_z[i];

        tr_y_xxyyyz_xxz[i] = tr_y_xxyyy_xx[i] * fe_0 + tr_y_xxyyy_xxz[i] * pa_z[i];

        tr_y_xxyyyz_xyy[i] = tr_y_xxyyy_xyy[i] * pa_z[i];

        tr_y_xxyyyz_xyz[i] = tr_y_xxyyy_xy[i] * fe_0 + tr_y_xxyyy_xyz[i] * pa_z[i];

        tr_y_xxyyyz_xzz[i] = 2.0 * tr_y_xxyyy_xz[i] * fe_0 + tr_y_xxyyy_xzz[i] * pa_z[i];

        tr_y_xxyyyz_yyy[i] = tr_y_xxyyy_yyy[i] * pa_z[i];

        tr_y_xxyyyz_yyz[i] = tr_y_yyyz_yyz[i] * fe_0 + tr_y_xyyyz_yyz[i] * pa_x[i];

        tr_y_xxyyyz_yzz[i] = tr_y_yyyz_yzz[i] * fe_0 + tr_y_xyyyz_yzz[i] * pa_x[i];

        tr_y_xxyyyz_zzz[i] = tr_y_yyyz_zzz[i] * fe_0 + tr_y_xyyyz_zzz[i] * pa_x[i];
    }

    // Set up 400-410 components of targeted buffer : IF

    auto tr_y_xxyyzz_xxx = pbuffer.data(idx_dip_if + 400);

    auto tr_y_xxyyzz_xxy = pbuffer.data(idx_dip_if + 401);

    auto tr_y_xxyyzz_xxz = pbuffer.data(idx_dip_if + 402);

    auto tr_y_xxyyzz_xyy = pbuffer.data(idx_dip_if + 403);

    auto tr_y_xxyyzz_xyz = pbuffer.data(idx_dip_if + 404);

    auto tr_y_xxyyzz_xzz = pbuffer.data(idx_dip_if + 405);

    auto tr_y_xxyyzz_yyy = pbuffer.data(idx_dip_if + 406);

    auto tr_y_xxyyzz_yyz = pbuffer.data(idx_dip_if + 407);

    auto tr_y_xxyyzz_yzz = pbuffer.data(idx_dip_if + 408);

    auto tr_y_xxyyzz_zzz = pbuffer.data(idx_dip_if + 409);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyy_xxx, tr_y_xxyy_xxy, tr_y_xxyy_xyy, tr_y_xxyyz_xxx, tr_y_xxyyz_xxy, tr_y_xxyyz_xyy, tr_y_xxyyzz_xxx, tr_y_xxyyzz_xxy, tr_y_xxyyzz_xxz, tr_y_xxyyzz_xyy, tr_y_xxyyzz_xyz, tr_y_xxyyzz_xzz, tr_y_xxyyzz_yyy, tr_y_xxyyzz_yyz, tr_y_xxyyzz_yzz, tr_y_xxyyzz_zzz, tr_y_xyyzz_xxz, tr_y_xyyzz_xyz, tr_y_xyyzz_xz, tr_y_xyyzz_xzz, tr_y_xyyzz_yyy, tr_y_xyyzz_yyz, tr_y_xyyzz_yz, tr_y_xyyzz_yzz, tr_y_xyyzz_zz, tr_y_xyyzz_zzz, tr_y_yyzz_xxz, tr_y_yyzz_xyz, tr_y_yyzz_xzz, tr_y_yyzz_yyy, tr_y_yyzz_yyz, tr_y_yyzz_yzz, tr_y_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyzz_xxx[i] = tr_y_xxyy_xxx[i] * fe_0 + tr_y_xxyyz_xxx[i] * pa_z[i];

        tr_y_xxyyzz_xxy[i] = tr_y_xxyy_xxy[i] * fe_0 + tr_y_xxyyz_xxy[i] * pa_z[i];

        tr_y_xxyyzz_xxz[i] = tr_y_yyzz_xxz[i] * fe_0 + 2.0 * tr_y_xyyzz_xz[i] * fe_0 + tr_y_xyyzz_xxz[i] * pa_x[i];

        tr_y_xxyyzz_xyy[i] = tr_y_xxyy_xyy[i] * fe_0 + tr_y_xxyyz_xyy[i] * pa_z[i];

        tr_y_xxyyzz_xyz[i] = tr_y_yyzz_xyz[i] * fe_0 + tr_y_xyyzz_yz[i] * fe_0 + tr_y_xyyzz_xyz[i] * pa_x[i];

        tr_y_xxyyzz_xzz[i] = tr_y_yyzz_xzz[i] * fe_0 + tr_y_xyyzz_zz[i] * fe_0 + tr_y_xyyzz_xzz[i] * pa_x[i];

        tr_y_xxyyzz_yyy[i] = tr_y_yyzz_yyy[i] * fe_0 + tr_y_xyyzz_yyy[i] * pa_x[i];

        tr_y_xxyyzz_yyz[i] = tr_y_yyzz_yyz[i] * fe_0 + tr_y_xyyzz_yyz[i] * pa_x[i];

        tr_y_xxyyzz_yzz[i] = tr_y_yyzz_yzz[i] * fe_0 + tr_y_xyyzz_yzz[i] * pa_x[i];

        tr_y_xxyyzz_zzz[i] = tr_y_yyzz_zzz[i] * fe_0 + tr_y_xyyzz_zzz[i] * pa_x[i];
    }

    // Set up 410-420 components of targeted buffer : IF

    auto tr_y_xxyzzz_xxx = pbuffer.data(idx_dip_if + 410);

    auto tr_y_xxyzzz_xxy = pbuffer.data(idx_dip_if + 411);

    auto tr_y_xxyzzz_xxz = pbuffer.data(idx_dip_if + 412);

    auto tr_y_xxyzzz_xyy = pbuffer.data(idx_dip_if + 413);

    auto tr_y_xxyzzz_xyz = pbuffer.data(idx_dip_if + 414);

    auto tr_y_xxyzzz_xzz = pbuffer.data(idx_dip_if + 415);

    auto tr_y_xxyzzz_yyy = pbuffer.data(idx_dip_if + 416);

    auto tr_y_xxyzzz_yyz = pbuffer.data(idx_dip_if + 417);

    auto tr_y_xxyzzz_yzz = pbuffer.data(idx_dip_if + 418);

    auto tr_y_xxyzzz_zzz = pbuffer.data(idx_dip_if + 419);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxyz_xxy, tr_y_xxyz_xyy, tr_y_xxyzz_xxy, tr_y_xxyzz_xyy, tr_y_xxyzzz_xxx, tr_y_xxyzzz_xxy, tr_y_xxyzzz_xxz, tr_y_xxyzzz_xyy, tr_y_xxyzzz_xyz, tr_y_xxyzzz_xzz, tr_y_xxyzzz_yyy, tr_y_xxyzzz_yyz, tr_y_xxyzzz_yzz, tr_y_xxyzzz_zzz, tr_y_xxzzz_xxx, tr_y_xxzzz_xxz, tr_y_xxzzz_xzz, tr_y_xyzzz_xyz, tr_y_xyzzz_yyy, tr_y_xyzzz_yyz, tr_y_xyzzz_yz, tr_y_xyzzz_yzz, tr_y_xyzzz_zzz, tr_y_yzzz_xyz, tr_y_yzzz_yyy, tr_y_yzzz_yyz, tr_y_yzzz_yzz, tr_y_yzzz_zzz, ts_xxzzz_xxx, ts_xxzzz_xxz, ts_xxzzz_xzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzzz_xxx[i] = ts_xxzzz_xxx[i] * fe_0 + tr_y_xxzzz_xxx[i] * pa_y[i];

        tr_y_xxyzzz_xxy[i] = 2.0 * tr_y_xxyz_xxy[i] * fe_0 + tr_y_xxyzz_xxy[i] * pa_z[i];

        tr_y_xxyzzz_xxz[i] = ts_xxzzz_xxz[i] * fe_0 + tr_y_xxzzz_xxz[i] * pa_y[i];

        tr_y_xxyzzz_xyy[i] = 2.0 * tr_y_xxyz_xyy[i] * fe_0 + tr_y_xxyzz_xyy[i] * pa_z[i];

        tr_y_xxyzzz_xyz[i] = tr_y_yzzz_xyz[i] * fe_0 + tr_y_xyzzz_yz[i] * fe_0 + tr_y_xyzzz_xyz[i] * pa_x[i];

        tr_y_xxyzzz_xzz[i] = ts_xxzzz_xzz[i] * fe_0 + tr_y_xxzzz_xzz[i] * pa_y[i];

        tr_y_xxyzzz_yyy[i] = tr_y_yzzz_yyy[i] * fe_0 + tr_y_xyzzz_yyy[i] * pa_x[i];

        tr_y_xxyzzz_yyz[i] = tr_y_yzzz_yyz[i] * fe_0 + tr_y_xyzzz_yyz[i] * pa_x[i];

        tr_y_xxyzzz_yzz[i] = tr_y_yzzz_yzz[i] * fe_0 + tr_y_xyzzz_yzz[i] * pa_x[i];

        tr_y_xxyzzz_zzz[i] = tr_y_yzzz_zzz[i] * fe_0 + tr_y_xyzzz_zzz[i] * pa_x[i];
    }

    // Set up 420-430 components of targeted buffer : IF

    auto tr_y_xxzzzz_xxx = pbuffer.data(idx_dip_if + 420);

    auto tr_y_xxzzzz_xxy = pbuffer.data(idx_dip_if + 421);

    auto tr_y_xxzzzz_xxz = pbuffer.data(idx_dip_if + 422);

    auto tr_y_xxzzzz_xyy = pbuffer.data(idx_dip_if + 423);

    auto tr_y_xxzzzz_xyz = pbuffer.data(idx_dip_if + 424);

    auto tr_y_xxzzzz_xzz = pbuffer.data(idx_dip_if + 425);

    auto tr_y_xxzzzz_yyy = pbuffer.data(idx_dip_if + 426);

    auto tr_y_xxzzzz_yyz = pbuffer.data(idx_dip_if + 427);

    auto tr_y_xxzzzz_yzz = pbuffer.data(idx_dip_if + 428);

    auto tr_y_xxzzzz_zzz = pbuffer.data(idx_dip_if + 429);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxzz_xxx, tr_y_xxzz_xxy, tr_y_xxzz_xyy, tr_y_xxzzz_xxx, tr_y_xxzzz_xxy, tr_y_xxzzz_xyy, tr_y_xxzzzz_xxx, tr_y_xxzzzz_xxy, tr_y_xxzzzz_xxz, tr_y_xxzzzz_xyy, tr_y_xxzzzz_xyz, tr_y_xxzzzz_xzz, tr_y_xxzzzz_yyy, tr_y_xxzzzz_yyz, tr_y_xxzzzz_yzz, tr_y_xxzzzz_zzz, tr_y_xzzzz_xxz, tr_y_xzzzz_xyz, tr_y_xzzzz_xz, tr_y_xzzzz_xzz, tr_y_xzzzz_yyy, tr_y_xzzzz_yyz, tr_y_xzzzz_yz, tr_y_xzzzz_yzz, tr_y_xzzzz_zz, tr_y_xzzzz_zzz, tr_y_zzzz_xxz, tr_y_zzzz_xyz, tr_y_zzzz_xzz, tr_y_zzzz_yyy, tr_y_zzzz_yyz, tr_y_zzzz_yzz, tr_y_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzzz_xxx[i] = 3.0 * tr_y_xxzz_xxx[i] * fe_0 + tr_y_xxzzz_xxx[i] * pa_z[i];

        tr_y_xxzzzz_xxy[i] = 3.0 * tr_y_xxzz_xxy[i] * fe_0 + tr_y_xxzzz_xxy[i] * pa_z[i];

        tr_y_xxzzzz_xxz[i] = tr_y_zzzz_xxz[i] * fe_0 + 2.0 * tr_y_xzzzz_xz[i] * fe_0 + tr_y_xzzzz_xxz[i] * pa_x[i];

        tr_y_xxzzzz_xyy[i] = 3.0 * tr_y_xxzz_xyy[i] * fe_0 + tr_y_xxzzz_xyy[i] * pa_z[i];

        tr_y_xxzzzz_xyz[i] = tr_y_zzzz_xyz[i] * fe_0 + tr_y_xzzzz_yz[i] * fe_0 + tr_y_xzzzz_xyz[i] * pa_x[i];

        tr_y_xxzzzz_xzz[i] = tr_y_zzzz_xzz[i] * fe_0 + tr_y_xzzzz_zz[i] * fe_0 + tr_y_xzzzz_xzz[i] * pa_x[i];

        tr_y_xxzzzz_yyy[i] = tr_y_zzzz_yyy[i] * fe_0 + tr_y_xzzzz_yyy[i] * pa_x[i];

        tr_y_xxzzzz_yyz[i] = tr_y_zzzz_yyz[i] * fe_0 + tr_y_xzzzz_yyz[i] * pa_x[i];

        tr_y_xxzzzz_yzz[i] = tr_y_zzzz_yzz[i] * fe_0 + tr_y_xzzzz_yzz[i] * pa_x[i];

        tr_y_xxzzzz_zzz[i] = tr_y_zzzz_zzz[i] * fe_0 + tr_y_xzzzz_zzz[i] * pa_x[i];
    }

    // Set up 430-440 components of targeted buffer : IF

    auto tr_y_xyyyyy_xxx = pbuffer.data(idx_dip_if + 430);

    auto tr_y_xyyyyy_xxy = pbuffer.data(idx_dip_if + 431);

    auto tr_y_xyyyyy_xxz = pbuffer.data(idx_dip_if + 432);

    auto tr_y_xyyyyy_xyy = pbuffer.data(idx_dip_if + 433);

    auto tr_y_xyyyyy_xyz = pbuffer.data(idx_dip_if + 434);

    auto tr_y_xyyyyy_xzz = pbuffer.data(idx_dip_if + 435);

    auto tr_y_xyyyyy_yyy = pbuffer.data(idx_dip_if + 436);

    auto tr_y_xyyyyy_yyz = pbuffer.data(idx_dip_if + 437);

    auto tr_y_xyyyyy_yzz = pbuffer.data(idx_dip_if + 438);

    auto tr_y_xyyyyy_zzz = pbuffer.data(idx_dip_if + 439);

    #pragma omp simd aligned(pa_x, tr_y_xyyyyy_xxx, tr_y_xyyyyy_xxy, tr_y_xyyyyy_xxz, tr_y_xyyyyy_xyy, tr_y_xyyyyy_xyz, tr_y_xyyyyy_xzz, tr_y_xyyyyy_yyy, tr_y_xyyyyy_yyz, tr_y_xyyyyy_yzz, tr_y_xyyyyy_zzz, tr_y_yyyyy_xx, tr_y_yyyyy_xxx, tr_y_yyyyy_xxy, tr_y_yyyyy_xxz, tr_y_yyyyy_xy, tr_y_yyyyy_xyy, tr_y_yyyyy_xyz, tr_y_yyyyy_xz, tr_y_yyyyy_xzz, tr_y_yyyyy_yy, tr_y_yyyyy_yyy, tr_y_yyyyy_yyz, tr_y_yyyyy_yz, tr_y_yyyyy_yzz, tr_y_yyyyy_zz, tr_y_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyy_xxx[i] = 3.0 * tr_y_yyyyy_xx[i] * fe_0 + tr_y_yyyyy_xxx[i] * pa_x[i];

        tr_y_xyyyyy_xxy[i] = 2.0 * tr_y_yyyyy_xy[i] * fe_0 + tr_y_yyyyy_xxy[i] * pa_x[i];

        tr_y_xyyyyy_xxz[i] = 2.0 * tr_y_yyyyy_xz[i] * fe_0 + tr_y_yyyyy_xxz[i] * pa_x[i];

        tr_y_xyyyyy_xyy[i] = tr_y_yyyyy_yy[i] * fe_0 + tr_y_yyyyy_xyy[i] * pa_x[i];

        tr_y_xyyyyy_xyz[i] = tr_y_yyyyy_yz[i] * fe_0 + tr_y_yyyyy_xyz[i] * pa_x[i];

        tr_y_xyyyyy_xzz[i] = tr_y_yyyyy_zz[i] * fe_0 + tr_y_yyyyy_xzz[i] * pa_x[i];

        tr_y_xyyyyy_yyy[i] = tr_y_yyyyy_yyy[i] * pa_x[i];

        tr_y_xyyyyy_yyz[i] = tr_y_yyyyy_yyz[i] * pa_x[i];

        tr_y_xyyyyy_yzz[i] = tr_y_yyyyy_yzz[i] * pa_x[i];

        tr_y_xyyyyy_zzz[i] = tr_y_yyyyy_zzz[i] * pa_x[i];
    }

    // Set up 440-450 components of targeted buffer : IF

    auto tr_y_xyyyyz_xxx = pbuffer.data(idx_dip_if + 440);

    auto tr_y_xyyyyz_xxy = pbuffer.data(idx_dip_if + 441);

    auto tr_y_xyyyyz_xxz = pbuffer.data(idx_dip_if + 442);

    auto tr_y_xyyyyz_xyy = pbuffer.data(idx_dip_if + 443);

    auto tr_y_xyyyyz_xyz = pbuffer.data(idx_dip_if + 444);

    auto tr_y_xyyyyz_xzz = pbuffer.data(idx_dip_if + 445);

    auto tr_y_xyyyyz_yyy = pbuffer.data(idx_dip_if + 446);

    auto tr_y_xyyyyz_yyz = pbuffer.data(idx_dip_if + 447);

    auto tr_y_xyyyyz_yzz = pbuffer.data(idx_dip_if + 448);

    auto tr_y_xyyyyz_zzz = pbuffer.data(idx_dip_if + 449);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xyyyy_xxx, tr_y_xyyyy_xxy, tr_y_xyyyy_xyy, tr_y_xyyyyz_xxx, tr_y_xyyyyz_xxy, tr_y_xyyyyz_xxz, tr_y_xyyyyz_xyy, tr_y_xyyyyz_xyz, tr_y_xyyyyz_xzz, tr_y_xyyyyz_yyy, tr_y_xyyyyz_yyz, tr_y_xyyyyz_yzz, tr_y_xyyyyz_zzz, tr_y_yyyyz_xxz, tr_y_yyyyz_xyz, tr_y_yyyyz_xz, tr_y_yyyyz_xzz, tr_y_yyyyz_yyy, tr_y_yyyyz_yyz, tr_y_yyyyz_yz, tr_y_yyyyz_yzz, tr_y_yyyyz_zz, tr_y_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyz_xxx[i] = tr_y_xyyyy_xxx[i] * pa_z[i];

        tr_y_xyyyyz_xxy[i] = tr_y_xyyyy_xxy[i] * pa_z[i];

        tr_y_xyyyyz_xxz[i] = 2.0 * tr_y_yyyyz_xz[i] * fe_0 + tr_y_yyyyz_xxz[i] * pa_x[i];

        tr_y_xyyyyz_xyy[i] = tr_y_xyyyy_xyy[i] * pa_z[i];

        tr_y_xyyyyz_xyz[i] = tr_y_yyyyz_yz[i] * fe_0 + tr_y_yyyyz_xyz[i] * pa_x[i];

        tr_y_xyyyyz_xzz[i] = tr_y_yyyyz_zz[i] * fe_0 + tr_y_yyyyz_xzz[i] * pa_x[i];

        tr_y_xyyyyz_yyy[i] = tr_y_yyyyz_yyy[i] * pa_x[i];

        tr_y_xyyyyz_yyz[i] = tr_y_yyyyz_yyz[i] * pa_x[i];

        tr_y_xyyyyz_yzz[i] = tr_y_yyyyz_yzz[i] * pa_x[i];

        tr_y_xyyyyz_zzz[i] = tr_y_yyyyz_zzz[i] * pa_x[i];
    }

    // Set up 450-460 components of targeted buffer : IF

    auto tr_y_xyyyzz_xxx = pbuffer.data(idx_dip_if + 450);

    auto tr_y_xyyyzz_xxy = pbuffer.data(idx_dip_if + 451);

    auto tr_y_xyyyzz_xxz = pbuffer.data(idx_dip_if + 452);

    auto tr_y_xyyyzz_xyy = pbuffer.data(idx_dip_if + 453);

    auto tr_y_xyyyzz_xyz = pbuffer.data(idx_dip_if + 454);

    auto tr_y_xyyyzz_xzz = pbuffer.data(idx_dip_if + 455);

    auto tr_y_xyyyzz_yyy = pbuffer.data(idx_dip_if + 456);

    auto tr_y_xyyyzz_yyz = pbuffer.data(idx_dip_if + 457);

    auto tr_y_xyyyzz_yzz = pbuffer.data(idx_dip_if + 458);

    auto tr_y_xyyyzz_zzz = pbuffer.data(idx_dip_if + 459);

    #pragma omp simd aligned(pa_x, tr_y_xyyyzz_xxx, tr_y_xyyyzz_xxy, tr_y_xyyyzz_xxz, tr_y_xyyyzz_xyy, tr_y_xyyyzz_xyz, tr_y_xyyyzz_xzz, tr_y_xyyyzz_yyy, tr_y_xyyyzz_yyz, tr_y_xyyyzz_yzz, tr_y_xyyyzz_zzz, tr_y_yyyzz_xx, tr_y_yyyzz_xxx, tr_y_yyyzz_xxy, tr_y_yyyzz_xxz, tr_y_yyyzz_xy, tr_y_yyyzz_xyy, tr_y_yyyzz_xyz, tr_y_yyyzz_xz, tr_y_yyyzz_xzz, tr_y_yyyzz_yy, tr_y_yyyzz_yyy, tr_y_yyyzz_yyz, tr_y_yyyzz_yz, tr_y_yyyzz_yzz, tr_y_yyyzz_zz, tr_y_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyzz_xxx[i] = 3.0 * tr_y_yyyzz_xx[i] * fe_0 + tr_y_yyyzz_xxx[i] * pa_x[i];

        tr_y_xyyyzz_xxy[i] = 2.0 * tr_y_yyyzz_xy[i] * fe_0 + tr_y_yyyzz_xxy[i] * pa_x[i];

        tr_y_xyyyzz_xxz[i] = 2.0 * tr_y_yyyzz_xz[i] * fe_0 + tr_y_yyyzz_xxz[i] * pa_x[i];

        tr_y_xyyyzz_xyy[i] = tr_y_yyyzz_yy[i] * fe_0 + tr_y_yyyzz_xyy[i] * pa_x[i];

        tr_y_xyyyzz_xyz[i] = tr_y_yyyzz_yz[i] * fe_0 + tr_y_yyyzz_xyz[i] * pa_x[i];

        tr_y_xyyyzz_xzz[i] = tr_y_yyyzz_zz[i] * fe_0 + tr_y_yyyzz_xzz[i] * pa_x[i];

        tr_y_xyyyzz_yyy[i] = tr_y_yyyzz_yyy[i] * pa_x[i];

        tr_y_xyyyzz_yyz[i] = tr_y_yyyzz_yyz[i] * pa_x[i];

        tr_y_xyyyzz_yzz[i] = tr_y_yyyzz_yzz[i] * pa_x[i];

        tr_y_xyyyzz_zzz[i] = tr_y_yyyzz_zzz[i] * pa_x[i];
    }

    // Set up 460-470 components of targeted buffer : IF

    auto tr_y_xyyzzz_xxx = pbuffer.data(idx_dip_if + 460);

    auto tr_y_xyyzzz_xxy = pbuffer.data(idx_dip_if + 461);

    auto tr_y_xyyzzz_xxz = pbuffer.data(idx_dip_if + 462);

    auto tr_y_xyyzzz_xyy = pbuffer.data(idx_dip_if + 463);

    auto tr_y_xyyzzz_xyz = pbuffer.data(idx_dip_if + 464);

    auto tr_y_xyyzzz_xzz = pbuffer.data(idx_dip_if + 465);

    auto tr_y_xyyzzz_yyy = pbuffer.data(idx_dip_if + 466);

    auto tr_y_xyyzzz_yyz = pbuffer.data(idx_dip_if + 467);

    auto tr_y_xyyzzz_yzz = pbuffer.data(idx_dip_if + 468);

    auto tr_y_xyyzzz_zzz = pbuffer.data(idx_dip_if + 469);

    #pragma omp simd aligned(pa_x, tr_y_xyyzzz_xxx, tr_y_xyyzzz_xxy, tr_y_xyyzzz_xxz, tr_y_xyyzzz_xyy, tr_y_xyyzzz_xyz, tr_y_xyyzzz_xzz, tr_y_xyyzzz_yyy, tr_y_xyyzzz_yyz, tr_y_xyyzzz_yzz, tr_y_xyyzzz_zzz, tr_y_yyzzz_xx, tr_y_yyzzz_xxx, tr_y_yyzzz_xxy, tr_y_yyzzz_xxz, tr_y_yyzzz_xy, tr_y_yyzzz_xyy, tr_y_yyzzz_xyz, tr_y_yyzzz_xz, tr_y_yyzzz_xzz, tr_y_yyzzz_yy, tr_y_yyzzz_yyy, tr_y_yyzzz_yyz, tr_y_yyzzz_yz, tr_y_yyzzz_yzz, tr_y_yyzzz_zz, tr_y_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzzz_xxx[i] = 3.0 * tr_y_yyzzz_xx[i] * fe_0 + tr_y_yyzzz_xxx[i] * pa_x[i];

        tr_y_xyyzzz_xxy[i] = 2.0 * tr_y_yyzzz_xy[i] * fe_0 + tr_y_yyzzz_xxy[i] * pa_x[i];

        tr_y_xyyzzz_xxz[i] = 2.0 * tr_y_yyzzz_xz[i] * fe_0 + tr_y_yyzzz_xxz[i] * pa_x[i];

        tr_y_xyyzzz_xyy[i] = tr_y_yyzzz_yy[i] * fe_0 + tr_y_yyzzz_xyy[i] * pa_x[i];

        tr_y_xyyzzz_xyz[i] = tr_y_yyzzz_yz[i] * fe_0 + tr_y_yyzzz_xyz[i] * pa_x[i];

        tr_y_xyyzzz_xzz[i] = tr_y_yyzzz_zz[i] * fe_0 + tr_y_yyzzz_xzz[i] * pa_x[i];

        tr_y_xyyzzz_yyy[i] = tr_y_yyzzz_yyy[i] * pa_x[i];

        tr_y_xyyzzz_yyz[i] = tr_y_yyzzz_yyz[i] * pa_x[i];

        tr_y_xyyzzz_yzz[i] = tr_y_yyzzz_yzz[i] * pa_x[i];

        tr_y_xyyzzz_zzz[i] = tr_y_yyzzz_zzz[i] * pa_x[i];
    }

    // Set up 470-480 components of targeted buffer : IF

    auto tr_y_xyzzzz_xxx = pbuffer.data(idx_dip_if + 470);

    auto tr_y_xyzzzz_xxy = pbuffer.data(idx_dip_if + 471);

    auto tr_y_xyzzzz_xxz = pbuffer.data(idx_dip_if + 472);

    auto tr_y_xyzzzz_xyy = pbuffer.data(idx_dip_if + 473);

    auto tr_y_xyzzzz_xyz = pbuffer.data(idx_dip_if + 474);

    auto tr_y_xyzzzz_xzz = pbuffer.data(idx_dip_if + 475);

    auto tr_y_xyzzzz_yyy = pbuffer.data(idx_dip_if + 476);

    auto tr_y_xyzzzz_yyz = pbuffer.data(idx_dip_if + 477);

    auto tr_y_xyzzzz_yzz = pbuffer.data(idx_dip_if + 478);

    auto tr_y_xyzzzz_zzz = pbuffer.data(idx_dip_if + 479);

    #pragma omp simd aligned(pa_x, tr_y_xyzzzz_xxx, tr_y_xyzzzz_xxy, tr_y_xyzzzz_xxz, tr_y_xyzzzz_xyy, tr_y_xyzzzz_xyz, tr_y_xyzzzz_xzz, tr_y_xyzzzz_yyy, tr_y_xyzzzz_yyz, tr_y_xyzzzz_yzz, tr_y_xyzzzz_zzz, tr_y_yzzzz_xx, tr_y_yzzzz_xxx, tr_y_yzzzz_xxy, tr_y_yzzzz_xxz, tr_y_yzzzz_xy, tr_y_yzzzz_xyy, tr_y_yzzzz_xyz, tr_y_yzzzz_xz, tr_y_yzzzz_xzz, tr_y_yzzzz_yy, tr_y_yzzzz_yyy, tr_y_yzzzz_yyz, tr_y_yzzzz_yz, tr_y_yzzzz_yzz, tr_y_yzzzz_zz, tr_y_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzzz_xxx[i] = 3.0 * tr_y_yzzzz_xx[i] * fe_0 + tr_y_yzzzz_xxx[i] * pa_x[i];

        tr_y_xyzzzz_xxy[i] = 2.0 * tr_y_yzzzz_xy[i] * fe_0 + tr_y_yzzzz_xxy[i] * pa_x[i];

        tr_y_xyzzzz_xxz[i] = 2.0 * tr_y_yzzzz_xz[i] * fe_0 + tr_y_yzzzz_xxz[i] * pa_x[i];

        tr_y_xyzzzz_xyy[i] = tr_y_yzzzz_yy[i] * fe_0 + tr_y_yzzzz_xyy[i] * pa_x[i];

        tr_y_xyzzzz_xyz[i] = tr_y_yzzzz_yz[i] * fe_0 + tr_y_yzzzz_xyz[i] * pa_x[i];

        tr_y_xyzzzz_xzz[i] = tr_y_yzzzz_zz[i] * fe_0 + tr_y_yzzzz_xzz[i] * pa_x[i];

        tr_y_xyzzzz_yyy[i] = tr_y_yzzzz_yyy[i] * pa_x[i];

        tr_y_xyzzzz_yyz[i] = tr_y_yzzzz_yyz[i] * pa_x[i];

        tr_y_xyzzzz_yzz[i] = tr_y_yzzzz_yzz[i] * pa_x[i];

        tr_y_xyzzzz_zzz[i] = tr_y_yzzzz_zzz[i] * pa_x[i];
    }

    // Set up 480-490 components of targeted buffer : IF

    auto tr_y_xzzzzz_xxx = pbuffer.data(idx_dip_if + 480);

    auto tr_y_xzzzzz_xxy = pbuffer.data(idx_dip_if + 481);

    auto tr_y_xzzzzz_xxz = pbuffer.data(idx_dip_if + 482);

    auto tr_y_xzzzzz_xyy = pbuffer.data(idx_dip_if + 483);

    auto tr_y_xzzzzz_xyz = pbuffer.data(idx_dip_if + 484);

    auto tr_y_xzzzzz_xzz = pbuffer.data(idx_dip_if + 485);

    auto tr_y_xzzzzz_yyy = pbuffer.data(idx_dip_if + 486);

    auto tr_y_xzzzzz_yyz = pbuffer.data(idx_dip_if + 487);

    auto tr_y_xzzzzz_yzz = pbuffer.data(idx_dip_if + 488);

    auto tr_y_xzzzzz_zzz = pbuffer.data(idx_dip_if + 489);

    #pragma omp simd aligned(pa_x, tr_y_xzzzzz_xxx, tr_y_xzzzzz_xxy, tr_y_xzzzzz_xxz, tr_y_xzzzzz_xyy, tr_y_xzzzzz_xyz, tr_y_xzzzzz_xzz, tr_y_xzzzzz_yyy, tr_y_xzzzzz_yyz, tr_y_xzzzzz_yzz, tr_y_xzzzzz_zzz, tr_y_zzzzz_xx, tr_y_zzzzz_xxx, tr_y_zzzzz_xxy, tr_y_zzzzz_xxz, tr_y_zzzzz_xy, tr_y_zzzzz_xyy, tr_y_zzzzz_xyz, tr_y_zzzzz_xz, tr_y_zzzzz_xzz, tr_y_zzzzz_yy, tr_y_zzzzz_yyy, tr_y_zzzzz_yyz, tr_y_zzzzz_yz, tr_y_zzzzz_yzz, tr_y_zzzzz_zz, tr_y_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzzz_xxx[i] = 3.0 * tr_y_zzzzz_xx[i] * fe_0 + tr_y_zzzzz_xxx[i] * pa_x[i];

        tr_y_xzzzzz_xxy[i] = 2.0 * tr_y_zzzzz_xy[i] * fe_0 + tr_y_zzzzz_xxy[i] * pa_x[i];

        tr_y_xzzzzz_xxz[i] = 2.0 * tr_y_zzzzz_xz[i] * fe_0 + tr_y_zzzzz_xxz[i] * pa_x[i];

        tr_y_xzzzzz_xyy[i] = tr_y_zzzzz_yy[i] * fe_0 + tr_y_zzzzz_xyy[i] * pa_x[i];

        tr_y_xzzzzz_xyz[i] = tr_y_zzzzz_yz[i] * fe_0 + tr_y_zzzzz_xyz[i] * pa_x[i];

        tr_y_xzzzzz_xzz[i] = tr_y_zzzzz_zz[i] * fe_0 + tr_y_zzzzz_xzz[i] * pa_x[i];

        tr_y_xzzzzz_yyy[i] = tr_y_zzzzz_yyy[i] * pa_x[i];

        tr_y_xzzzzz_yyz[i] = tr_y_zzzzz_yyz[i] * pa_x[i];

        tr_y_xzzzzz_yzz[i] = tr_y_zzzzz_yzz[i] * pa_x[i];

        tr_y_xzzzzz_zzz[i] = tr_y_zzzzz_zzz[i] * pa_x[i];
    }

    // Set up 490-500 components of targeted buffer : IF

    auto tr_y_yyyyyy_xxx = pbuffer.data(idx_dip_if + 490);

    auto tr_y_yyyyyy_xxy = pbuffer.data(idx_dip_if + 491);

    auto tr_y_yyyyyy_xxz = pbuffer.data(idx_dip_if + 492);

    auto tr_y_yyyyyy_xyy = pbuffer.data(idx_dip_if + 493);

    auto tr_y_yyyyyy_xyz = pbuffer.data(idx_dip_if + 494);

    auto tr_y_yyyyyy_xzz = pbuffer.data(idx_dip_if + 495);

    auto tr_y_yyyyyy_yyy = pbuffer.data(idx_dip_if + 496);

    auto tr_y_yyyyyy_yyz = pbuffer.data(idx_dip_if + 497);

    auto tr_y_yyyyyy_yzz = pbuffer.data(idx_dip_if + 498);

    auto tr_y_yyyyyy_zzz = pbuffer.data(idx_dip_if + 499);

    #pragma omp simd aligned(pa_y, tr_y_yyyy_xxx, tr_y_yyyy_xxy, tr_y_yyyy_xxz, tr_y_yyyy_xyy, tr_y_yyyy_xyz, tr_y_yyyy_xzz, tr_y_yyyy_yyy, tr_y_yyyy_yyz, tr_y_yyyy_yzz, tr_y_yyyy_zzz, tr_y_yyyyy_xx, tr_y_yyyyy_xxx, tr_y_yyyyy_xxy, tr_y_yyyyy_xxz, tr_y_yyyyy_xy, tr_y_yyyyy_xyy, tr_y_yyyyy_xyz, tr_y_yyyyy_xz, tr_y_yyyyy_xzz, tr_y_yyyyy_yy, tr_y_yyyyy_yyy, tr_y_yyyyy_yyz, tr_y_yyyyy_yz, tr_y_yyyyy_yzz, tr_y_yyyyy_zz, tr_y_yyyyy_zzz, tr_y_yyyyyy_xxx, tr_y_yyyyyy_xxy, tr_y_yyyyyy_xxz, tr_y_yyyyyy_xyy, tr_y_yyyyyy_xyz, tr_y_yyyyyy_xzz, tr_y_yyyyyy_yyy, tr_y_yyyyyy_yyz, tr_y_yyyyyy_yzz, tr_y_yyyyyy_zzz, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xzz, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yzz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyy_xxx[i] = 5.0 * tr_y_yyyy_xxx[i] * fe_0 + ts_yyyyy_xxx[i] * fe_0 + tr_y_yyyyy_xxx[i] * pa_y[i];

        tr_y_yyyyyy_xxy[i] = 5.0 * tr_y_yyyy_xxy[i] * fe_0 + tr_y_yyyyy_xx[i] * fe_0 + ts_yyyyy_xxy[i] * fe_0 + tr_y_yyyyy_xxy[i] * pa_y[i];

        tr_y_yyyyyy_xxz[i] = 5.0 * tr_y_yyyy_xxz[i] * fe_0 + ts_yyyyy_xxz[i] * fe_0 + tr_y_yyyyy_xxz[i] * pa_y[i];

        tr_y_yyyyyy_xyy[i] = 5.0 * tr_y_yyyy_xyy[i] * fe_0 + 2.0 * tr_y_yyyyy_xy[i] * fe_0 + ts_yyyyy_xyy[i] * fe_0 + tr_y_yyyyy_xyy[i] * pa_y[i];

        tr_y_yyyyyy_xyz[i] = 5.0 * tr_y_yyyy_xyz[i] * fe_0 + tr_y_yyyyy_xz[i] * fe_0 + ts_yyyyy_xyz[i] * fe_0 + tr_y_yyyyy_xyz[i] * pa_y[i];

        tr_y_yyyyyy_xzz[i] = 5.0 * tr_y_yyyy_xzz[i] * fe_0 + ts_yyyyy_xzz[i] * fe_0 + tr_y_yyyyy_xzz[i] * pa_y[i];

        tr_y_yyyyyy_yyy[i] = 5.0 * tr_y_yyyy_yyy[i] * fe_0 + 3.0 * tr_y_yyyyy_yy[i] * fe_0 + ts_yyyyy_yyy[i] * fe_0 + tr_y_yyyyy_yyy[i] * pa_y[i];

        tr_y_yyyyyy_yyz[i] = 5.0 * tr_y_yyyy_yyz[i] * fe_0 + 2.0 * tr_y_yyyyy_yz[i] * fe_0 + ts_yyyyy_yyz[i] * fe_0 + tr_y_yyyyy_yyz[i] * pa_y[i];

        tr_y_yyyyyy_yzz[i] = 5.0 * tr_y_yyyy_yzz[i] * fe_0 + tr_y_yyyyy_zz[i] * fe_0 + ts_yyyyy_yzz[i] * fe_0 + tr_y_yyyyy_yzz[i] * pa_y[i];

        tr_y_yyyyyy_zzz[i] = 5.0 * tr_y_yyyy_zzz[i] * fe_0 + ts_yyyyy_zzz[i] * fe_0 + tr_y_yyyyy_zzz[i] * pa_y[i];
    }

    // Set up 500-510 components of targeted buffer : IF

    auto tr_y_yyyyyz_xxx = pbuffer.data(idx_dip_if + 500);

    auto tr_y_yyyyyz_xxy = pbuffer.data(idx_dip_if + 501);

    auto tr_y_yyyyyz_xxz = pbuffer.data(idx_dip_if + 502);

    auto tr_y_yyyyyz_xyy = pbuffer.data(idx_dip_if + 503);

    auto tr_y_yyyyyz_xyz = pbuffer.data(idx_dip_if + 504);

    auto tr_y_yyyyyz_xzz = pbuffer.data(idx_dip_if + 505);

    auto tr_y_yyyyyz_yyy = pbuffer.data(idx_dip_if + 506);

    auto tr_y_yyyyyz_yyz = pbuffer.data(idx_dip_if + 507);

    auto tr_y_yyyyyz_yzz = pbuffer.data(idx_dip_if + 508);

    auto tr_y_yyyyyz_zzz = pbuffer.data(idx_dip_if + 509);

    #pragma omp simd aligned(pa_z, tr_y_yyyyy_xx, tr_y_yyyyy_xxx, tr_y_yyyyy_xxy, tr_y_yyyyy_xxz, tr_y_yyyyy_xy, tr_y_yyyyy_xyy, tr_y_yyyyy_xyz, tr_y_yyyyy_xz, tr_y_yyyyy_xzz, tr_y_yyyyy_yy, tr_y_yyyyy_yyy, tr_y_yyyyy_yyz, tr_y_yyyyy_yz, tr_y_yyyyy_yzz, tr_y_yyyyy_zz, tr_y_yyyyy_zzz, tr_y_yyyyyz_xxx, tr_y_yyyyyz_xxy, tr_y_yyyyyz_xxz, tr_y_yyyyyz_xyy, tr_y_yyyyyz_xyz, tr_y_yyyyyz_xzz, tr_y_yyyyyz_yyy, tr_y_yyyyyz_yyz, tr_y_yyyyyz_yzz, tr_y_yyyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyz_xxx[i] = tr_y_yyyyy_xxx[i] * pa_z[i];

        tr_y_yyyyyz_xxy[i] = tr_y_yyyyy_xxy[i] * pa_z[i];

        tr_y_yyyyyz_xxz[i] = tr_y_yyyyy_xx[i] * fe_0 + tr_y_yyyyy_xxz[i] * pa_z[i];

        tr_y_yyyyyz_xyy[i] = tr_y_yyyyy_xyy[i] * pa_z[i];

        tr_y_yyyyyz_xyz[i] = tr_y_yyyyy_xy[i] * fe_0 + tr_y_yyyyy_xyz[i] * pa_z[i];

        tr_y_yyyyyz_xzz[i] = 2.0 * tr_y_yyyyy_xz[i] * fe_0 + tr_y_yyyyy_xzz[i] * pa_z[i];

        tr_y_yyyyyz_yyy[i] = tr_y_yyyyy_yyy[i] * pa_z[i];

        tr_y_yyyyyz_yyz[i] = tr_y_yyyyy_yy[i] * fe_0 + tr_y_yyyyy_yyz[i] * pa_z[i];

        tr_y_yyyyyz_yzz[i] = 2.0 * tr_y_yyyyy_yz[i] * fe_0 + tr_y_yyyyy_yzz[i] * pa_z[i];

        tr_y_yyyyyz_zzz[i] = 3.0 * tr_y_yyyyy_zz[i] * fe_0 + tr_y_yyyyy_zzz[i] * pa_z[i];
    }

    // Set up 510-520 components of targeted buffer : IF

    auto tr_y_yyyyzz_xxx = pbuffer.data(idx_dip_if + 510);

    auto tr_y_yyyyzz_xxy = pbuffer.data(idx_dip_if + 511);

    auto tr_y_yyyyzz_xxz = pbuffer.data(idx_dip_if + 512);

    auto tr_y_yyyyzz_xyy = pbuffer.data(idx_dip_if + 513);

    auto tr_y_yyyyzz_xyz = pbuffer.data(idx_dip_if + 514);

    auto tr_y_yyyyzz_xzz = pbuffer.data(idx_dip_if + 515);

    auto tr_y_yyyyzz_yyy = pbuffer.data(idx_dip_if + 516);

    auto tr_y_yyyyzz_yyz = pbuffer.data(idx_dip_if + 517);

    auto tr_y_yyyyzz_yzz = pbuffer.data(idx_dip_if + 518);

    auto tr_y_yyyyzz_zzz = pbuffer.data(idx_dip_if + 519);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyyy_xxx, tr_y_yyyy_xxy, tr_y_yyyy_xyy, tr_y_yyyy_xyz, tr_y_yyyy_yyy, tr_y_yyyy_yyz, tr_y_yyyy_yzz, tr_y_yyyyz_xxx, tr_y_yyyyz_xxy, tr_y_yyyyz_xy, tr_y_yyyyz_xyy, tr_y_yyyyz_xyz, tr_y_yyyyz_yy, tr_y_yyyyz_yyy, tr_y_yyyyz_yyz, tr_y_yyyyz_yz, tr_y_yyyyz_yzz, tr_y_yyyyzz_xxx, tr_y_yyyyzz_xxy, tr_y_yyyyzz_xxz, tr_y_yyyyzz_xyy, tr_y_yyyyzz_xyz, tr_y_yyyyzz_xzz, tr_y_yyyyzz_yyy, tr_y_yyyyzz_yyz, tr_y_yyyyzz_yzz, tr_y_yyyyzz_zzz, tr_y_yyyzz_xxz, tr_y_yyyzz_xzz, tr_y_yyyzz_zzz, tr_y_yyzz_xxz, tr_y_yyzz_xzz, tr_y_yyzz_zzz, ts_yyyzz_xxz, ts_yyyzz_xzz, ts_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyzz_xxx[i] = tr_y_yyyy_xxx[i] * fe_0 + tr_y_yyyyz_xxx[i] * pa_z[i];

        tr_y_yyyyzz_xxy[i] = tr_y_yyyy_xxy[i] * fe_0 + tr_y_yyyyz_xxy[i] * pa_z[i];

        tr_y_yyyyzz_xxz[i] = 3.0 * tr_y_yyzz_xxz[i] * fe_0 + ts_yyyzz_xxz[i] * fe_0 + tr_y_yyyzz_xxz[i] * pa_y[i];

        tr_y_yyyyzz_xyy[i] = tr_y_yyyy_xyy[i] * fe_0 + tr_y_yyyyz_xyy[i] * pa_z[i];

        tr_y_yyyyzz_xyz[i] = tr_y_yyyy_xyz[i] * fe_0 + tr_y_yyyyz_xy[i] * fe_0 + tr_y_yyyyz_xyz[i] * pa_z[i];

        tr_y_yyyyzz_xzz[i] = 3.0 * tr_y_yyzz_xzz[i] * fe_0 + ts_yyyzz_xzz[i] * fe_0 + tr_y_yyyzz_xzz[i] * pa_y[i];

        tr_y_yyyyzz_yyy[i] = tr_y_yyyy_yyy[i] * fe_0 + tr_y_yyyyz_yyy[i] * pa_z[i];

        tr_y_yyyyzz_yyz[i] = tr_y_yyyy_yyz[i] * fe_0 + tr_y_yyyyz_yy[i] * fe_0 + tr_y_yyyyz_yyz[i] * pa_z[i];

        tr_y_yyyyzz_yzz[i] = tr_y_yyyy_yzz[i] * fe_0 + 2.0 * tr_y_yyyyz_yz[i] * fe_0 + tr_y_yyyyz_yzz[i] * pa_z[i];

        tr_y_yyyyzz_zzz[i] = 3.0 * tr_y_yyzz_zzz[i] * fe_0 + ts_yyyzz_zzz[i] * fe_0 + tr_y_yyyzz_zzz[i] * pa_y[i];
    }

    // Set up 520-530 components of targeted buffer : IF

    auto tr_y_yyyzzz_xxx = pbuffer.data(idx_dip_if + 520);

    auto tr_y_yyyzzz_xxy = pbuffer.data(idx_dip_if + 521);

    auto tr_y_yyyzzz_xxz = pbuffer.data(idx_dip_if + 522);

    auto tr_y_yyyzzz_xyy = pbuffer.data(idx_dip_if + 523);

    auto tr_y_yyyzzz_xyz = pbuffer.data(idx_dip_if + 524);

    auto tr_y_yyyzzz_xzz = pbuffer.data(idx_dip_if + 525);

    auto tr_y_yyyzzz_yyy = pbuffer.data(idx_dip_if + 526);

    auto tr_y_yyyzzz_yyz = pbuffer.data(idx_dip_if + 527);

    auto tr_y_yyyzzz_yzz = pbuffer.data(idx_dip_if + 528);

    auto tr_y_yyyzzz_zzz = pbuffer.data(idx_dip_if + 529);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyyz_xxx, tr_y_yyyz_xxy, tr_y_yyyz_xyy, tr_y_yyyz_xyz, tr_y_yyyz_yyy, tr_y_yyyz_yyz, tr_y_yyyz_yzz, tr_y_yyyzz_xxx, tr_y_yyyzz_xxy, tr_y_yyyzz_xy, tr_y_yyyzz_xyy, tr_y_yyyzz_xyz, tr_y_yyyzz_yy, tr_y_yyyzz_yyy, tr_y_yyyzz_yyz, tr_y_yyyzz_yz, tr_y_yyyzz_yzz, tr_y_yyyzzz_xxx, tr_y_yyyzzz_xxy, tr_y_yyyzzz_xxz, tr_y_yyyzzz_xyy, tr_y_yyyzzz_xyz, tr_y_yyyzzz_xzz, tr_y_yyyzzz_yyy, tr_y_yyyzzz_yyz, tr_y_yyyzzz_yzz, tr_y_yyyzzz_zzz, tr_y_yyzzz_xxz, tr_y_yyzzz_xzz, tr_y_yyzzz_zzz, tr_y_yzzz_xxz, tr_y_yzzz_xzz, tr_y_yzzz_zzz, ts_yyzzz_xxz, ts_yyzzz_xzz, ts_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzzz_xxx[i] = 2.0 * tr_y_yyyz_xxx[i] * fe_0 + tr_y_yyyzz_xxx[i] * pa_z[i];

        tr_y_yyyzzz_xxy[i] = 2.0 * tr_y_yyyz_xxy[i] * fe_0 + tr_y_yyyzz_xxy[i] * pa_z[i];

        tr_y_yyyzzz_xxz[i] = 2.0 * tr_y_yzzz_xxz[i] * fe_0 + ts_yyzzz_xxz[i] * fe_0 + tr_y_yyzzz_xxz[i] * pa_y[i];

        tr_y_yyyzzz_xyy[i] = 2.0 * tr_y_yyyz_xyy[i] * fe_0 + tr_y_yyyzz_xyy[i] * pa_z[i];

        tr_y_yyyzzz_xyz[i] = 2.0 * tr_y_yyyz_xyz[i] * fe_0 + tr_y_yyyzz_xy[i] * fe_0 + tr_y_yyyzz_xyz[i] * pa_z[i];

        tr_y_yyyzzz_xzz[i] = 2.0 * tr_y_yzzz_xzz[i] * fe_0 + ts_yyzzz_xzz[i] * fe_0 + tr_y_yyzzz_xzz[i] * pa_y[i];

        tr_y_yyyzzz_yyy[i] = 2.0 * tr_y_yyyz_yyy[i] * fe_0 + tr_y_yyyzz_yyy[i] * pa_z[i];

        tr_y_yyyzzz_yyz[i] = 2.0 * tr_y_yyyz_yyz[i] * fe_0 + tr_y_yyyzz_yy[i] * fe_0 + tr_y_yyyzz_yyz[i] * pa_z[i];

        tr_y_yyyzzz_yzz[i] = 2.0 * tr_y_yyyz_yzz[i] * fe_0 + 2.0 * tr_y_yyyzz_yz[i] * fe_0 + tr_y_yyyzz_yzz[i] * pa_z[i];

        tr_y_yyyzzz_zzz[i] = 2.0 * tr_y_yzzz_zzz[i] * fe_0 + ts_yyzzz_zzz[i] * fe_0 + tr_y_yyzzz_zzz[i] * pa_y[i];
    }

    // Set up 530-540 components of targeted buffer : IF

    auto tr_y_yyzzzz_xxx = pbuffer.data(idx_dip_if + 530);

    auto tr_y_yyzzzz_xxy = pbuffer.data(idx_dip_if + 531);

    auto tr_y_yyzzzz_xxz = pbuffer.data(idx_dip_if + 532);

    auto tr_y_yyzzzz_xyy = pbuffer.data(idx_dip_if + 533);

    auto tr_y_yyzzzz_xyz = pbuffer.data(idx_dip_if + 534);

    auto tr_y_yyzzzz_xzz = pbuffer.data(idx_dip_if + 535);

    auto tr_y_yyzzzz_yyy = pbuffer.data(idx_dip_if + 536);

    auto tr_y_yyzzzz_yyz = pbuffer.data(idx_dip_if + 537);

    auto tr_y_yyzzzz_yzz = pbuffer.data(idx_dip_if + 538);

    auto tr_y_yyzzzz_zzz = pbuffer.data(idx_dip_if + 539);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyzz_xxx, tr_y_yyzz_xxy, tr_y_yyzz_xyy, tr_y_yyzz_xyz, tr_y_yyzz_yyy, tr_y_yyzz_yyz, tr_y_yyzz_yzz, tr_y_yyzzz_xxx, tr_y_yyzzz_xxy, tr_y_yyzzz_xy, tr_y_yyzzz_xyy, tr_y_yyzzz_xyz, tr_y_yyzzz_yy, tr_y_yyzzz_yyy, tr_y_yyzzz_yyz, tr_y_yyzzz_yz, tr_y_yyzzz_yzz, tr_y_yyzzzz_xxx, tr_y_yyzzzz_xxy, tr_y_yyzzzz_xxz, tr_y_yyzzzz_xyy, tr_y_yyzzzz_xyz, tr_y_yyzzzz_xzz, tr_y_yyzzzz_yyy, tr_y_yyzzzz_yyz, tr_y_yyzzzz_yzz, tr_y_yyzzzz_zzz, tr_y_yzzzz_xxz, tr_y_yzzzz_xzz, tr_y_yzzzz_zzz, tr_y_zzzz_xxz, tr_y_zzzz_xzz, tr_y_zzzz_zzz, ts_yzzzz_xxz, ts_yzzzz_xzz, ts_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzzz_xxx[i] = 3.0 * tr_y_yyzz_xxx[i] * fe_0 + tr_y_yyzzz_xxx[i] * pa_z[i];

        tr_y_yyzzzz_xxy[i] = 3.0 * tr_y_yyzz_xxy[i] * fe_0 + tr_y_yyzzz_xxy[i] * pa_z[i];

        tr_y_yyzzzz_xxz[i] = tr_y_zzzz_xxz[i] * fe_0 + ts_yzzzz_xxz[i] * fe_0 + tr_y_yzzzz_xxz[i] * pa_y[i];

        tr_y_yyzzzz_xyy[i] = 3.0 * tr_y_yyzz_xyy[i] * fe_0 + tr_y_yyzzz_xyy[i] * pa_z[i];

        tr_y_yyzzzz_xyz[i] = 3.0 * tr_y_yyzz_xyz[i] * fe_0 + tr_y_yyzzz_xy[i] * fe_0 + tr_y_yyzzz_xyz[i] * pa_z[i];

        tr_y_yyzzzz_xzz[i] = tr_y_zzzz_xzz[i] * fe_0 + ts_yzzzz_xzz[i] * fe_0 + tr_y_yzzzz_xzz[i] * pa_y[i];

        tr_y_yyzzzz_yyy[i] = 3.0 * tr_y_yyzz_yyy[i] * fe_0 + tr_y_yyzzz_yyy[i] * pa_z[i];

        tr_y_yyzzzz_yyz[i] = 3.0 * tr_y_yyzz_yyz[i] * fe_0 + tr_y_yyzzz_yy[i] * fe_0 + tr_y_yyzzz_yyz[i] * pa_z[i];

        tr_y_yyzzzz_yzz[i] = 3.0 * tr_y_yyzz_yzz[i] * fe_0 + 2.0 * tr_y_yyzzz_yz[i] * fe_0 + tr_y_yyzzz_yzz[i] * pa_z[i];

        tr_y_yyzzzz_zzz[i] = tr_y_zzzz_zzz[i] * fe_0 + ts_yzzzz_zzz[i] * fe_0 + tr_y_yzzzz_zzz[i] * pa_y[i];
    }

    // Set up 540-550 components of targeted buffer : IF

    auto tr_y_yzzzzz_xxx = pbuffer.data(idx_dip_if + 540);

    auto tr_y_yzzzzz_xxy = pbuffer.data(idx_dip_if + 541);

    auto tr_y_yzzzzz_xxz = pbuffer.data(idx_dip_if + 542);

    auto tr_y_yzzzzz_xyy = pbuffer.data(idx_dip_if + 543);

    auto tr_y_yzzzzz_xyz = pbuffer.data(idx_dip_if + 544);

    auto tr_y_yzzzzz_xzz = pbuffer.data(idx_dip_if + 545);

    auto tr_y_yzzzzz_yyy = pbuffer.data(idx_dip_if + 546);

    auto tr_y_yzzzzz_yyz = pbuffer.data(idx_dip_if + 547);

    auto tr_y_yzzzzz_yzz = pbuffer.data(idx_dip_if + 548);

    auto tr_y_yzzzzz_zzz = pbuffer.data(idx_dip_if + 549);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yzzz_xxy, tr_y_yzzz_xyy, tr_y_yzzz_yyy, tr_y_yzzzz_xxy, tr_y_yzzzz_xyy, tr_y_yzzzz_yyy, tr_y_yzzzzz_xxx, tr_y_yzzzzz_xxy, tr_y_yzzzzz_xxz, tr_y_yzzzzz_xyy, tr_y_yzzzzz_xyz, tr_y_yzzzzz_xzz, tr_y_yzzzzz_yyy, tr_y_yzzzzz_yyz, tr_y_yzzzzz_yzz, tr_y_yzzzzz_zzz, tr_y_zzzzz_xxx, tr_y_zzzzz_xxz, tr_y_zzzzz_xyz, tr_y_zzzzz_xz, tr_y_zzzzz_xzz, tr_y_zzzzz_yyz, tr_y_zzzzz_yz, tr_y_zzzzz_yzz, tr_y_zzzzz_zz, tr_y_zzzzz_zzz, ts_zzzzz_xxx, ts_zzzzz_xxz, ts_zzzzz_xyz, ts_zzzzz_xzz, ts_zzzzz_yyz, ts_zzzzz_yzz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzzz_xxx[i] = ts_zzzzz_xxx[i] * fe_0 + tr_y_zzzzz_xxx[i] * pa_y[i];

        tr_y_yzzzzz_xxy[i] = 4.0 * tr_y_yzzz_xxy[i] * fe_0 + tr_y_yzzzz_xxy[i] * pa_z[i];

        tr_y_yzzzzz_xxz[i] = ts_zzzzz_xxz[i] * fe_0 + tr_y_zzzzz_xxz[i] * pa_y[i];

        tr_y_yzzzzz_xyy[i] = 4.0 * tr_y_yzzz_xyy[i] * fe_0 + tr_y_yzzzz_xyy[i] * pa_z[i];

        tr_y_yzzzzz_xyz[i] = tr_y_zzzzz_xz[i] * fe_0 + ts_zzzzz_xyz[i] * fe_0 + tr_y_zzzzz_xyz[i] * pa_y[i];

        tr_y_yzzzzz_xzz[i] = ts_zzzzz_xzz[i] * fe_0 + tr_y_zzzzz_xzz[i] * pa_y[i];

        tr_y_yzzzzz_yyy[i] = 4.0 * tr_y_yzzz_yyy[i] * fe_0 + tr_y_yzzzz_yyy[i] * pa_z[i];

        tr_y_yzzzzz_yyz[i] = 2.0 * tr_y_zzzzz_yz[i] * fe_0 + ts_zzzzz_yyz[i] * fe_0 + tr_y_zzzzz_yyz[i] * pa_y[i];

        tr_y_yzzzzz_yzz[i] = tr_y_zzzzz_zz[i] * fe_0 + ts_zzzzz_yzz[i] * fe_0 + tr_y_zzzzz_yzz[i] * pa_y[i];

        tr_y_yzzzzz_zzz[i] = ts_zzzzz_zzz[i] * fe_0 + tr_y_zzzzz_zzz[i] * pa_y[i];
    }

    // Set up 550-560 components of targeted buffer : IF

    auto tr_y_zzzzzz_xxx = pbuffer.data(idx_dip_if + 550);

    auto tr_y_zzzzzz_xxy = pbuffer.data(idx_dip_if + 551);

    auto tr_y_zzzzzz_xxz = pbuffer.data(idx_dip_if + 552);

    auto tr_y_zzzzzz_xyy = pbuffer.data(idx_dip_if + 553);

    auto tr_y_zzzzzz_xyz = pbuffer.data(idx_dip_if + 554);

    auto tr_y_zzzzzz_xzz = pbuffer.data(idx_dip_if + 555);

    auto tr_y_zzzzzz_yyy = pbuffer.data(idx_dip_if + 556);

    auto tr_y_zzzzzz_yyz = pbuffer.data(idx_dip_if + 557);

    auto tr_y_zzzzzz_yzz = pbuffer.data(idx_dip_if + 558);

    auto tr_y_zzzzzz_zzz = pbuffer.data(idx_dip_if + 559);

    #pragma omp simd aligned(pa_z, tr_y_zzzz_xxx, tr_y_zzzz_xxy, tr_y_zzzz_xxz, tr_y_zzzz_xyy, tr_y_zzzz_xyz, tr_y_zzzz_xzz, tr_y_zzzz_yyy, tr_y_zzzz_yyz, tr_y_zzzz_yzz, tr_y_zzzz_zzz, tr_y_zzzzz_xx, tr_y_zzzzz_xxx, tr_y_zzzzz_xxy, tr_y_zzzzz_xxz, tr_y_zzzzz_xy, tr_y_zzzzz_xyy, tr_y_zzzzz_xyz, tr_y_zzzzz_xz, tr_y_zzzzz_xzz, tr_y_zzzzz_yy, tr_y_zzzzz_yyy, tr_y_zzzzz_yyz, tr_y_zzzzz_yz, tr_y_zzzzz_yzz, tr_y_zzzzz_zz, tr_y_zzzzz_zzz, tr_y_zzzzzz_xxx, tr_y_zzzzzz_xxy, tr_y_zzzzzz_xxz, tr_y_zzzzzz_xyy, tr_y_zzzzzz_xyz, tr_y_zzzzzz_xzz, tr_y_zzzzzz_yyy, tr_y_zzzzzz_yyz, tr_y_zzzzzz_yzz, tr_y_zzzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzzz_xxx[i] = 5.0 * tr_y_zzzz_xxx[i] * fe_0 + tr_y_zzzzz_xxx[i] * pa_z[i];

        tr_y_zzzzzz_xxy[i] = 5.0 * tr_y_zzzz_xxy[i] * fe_0 + tr_y_zzzzz_xxy[i] * pa_z[i];

        tr_y_zzzzzz_xxz[i] = 5.0 * tr_y_zzzz_xxz[i] * fe_0 + tr_y_zzzzz_xx[i] * fe_0 + tr_y_zzzzz_xxz[i] * pa_z[i];

        tr_y_zzzzzz_xyy[i] = 5.0 * tr_y_zzzz_xyy[i] * fe_0 + tr_y_zzzzz_xyy[i] * pa_z[i];

        tr_y_zzzzzz_xyz[i] = 5.0 * tr_y_zzzz_xyz[i] * fe_0 + tr_y_zzzzz_xy[i] * fe_0 + tr_y_zzzzz_xyz[i] * pa_z[i];

        tr_y_zzzzzz_xzz[i] = 5.0 * tr_y_zzzz_xzz[i] * fe_0 + 2.0 * tr_y_zzzzz_xz[i] * fe_0 + tr_y_zzzzz_xzz[i] * pa_z[i];

        tr_y_zzzzzz_yyy[i] = 5.0 * tr_y_zzzz_yyy[i] * fe_0 + tr_y_zzzzz_yyy[i] * pa_z[i];

        tr_y_zzzzzz_yyz[i] = 5.0 * tr_y_zzzz_yyz[i] * fe_0 + tr_y_zzzzz_yy[i] * fe_0 + tr_y_zzzzz_yyz[i] * pa_z[i];

        tr_y_zzzzzz_yzz[i] = 5.0 * tr_y_zzzz_yzz[i] * fe_0 + 2.0 * tr_y_zzzzz_yz[i] * fe_0 + tr_y_zzzzz_yzz[i] * pa_z[i];

        tr_y_zzzzzz_zzz[i] = 5.0 * tr_y_zzzz_zzz[i] * fe_0 + 3.0 * tr_y_zzzzz_zz[i] * fe_0 + tr_y_zzzzz_zzz[i] * pa_z[i];
    }

    // Set up 560-570 components of targeted buffer : IF

    auto tr_z_xxxxxx_xxx = pbuffer.data(idx_dip_if + 560);

    auto tr_z_xxxxxx_xxy = pbuffer.data(idx_dip_if + 561);

    auto tr_z_xxxxxx_xxz = pbuffer.data(idx_dip_if + 562);

    auto tr_z_xxxxxx_xyy = pbuffer.data(idx_dip_if + 563);

    auto tr_z_xxxxxx_xyz = pbuffer.data(idx_dip_if + 564);

    auto tr_z_xxxxxx_xzz = pbuffer.data(idx_dip_if + 565);

    auto tr_z_xxxxxx_yyy = pbuffer.data(idx_dip_if + 566);

    auto tr_z_xxxxxx_yyz = pbuffer.data(idx_dip_if + 567);

    auto tr_z_xxxxxx_yzz = pbuffer.data(idx_dip_if + 568);

    auto tr_z_xxxxxx_zzz = pbuffer.data(idx_dip_if + 569);

    #pragma omp simd aligned(pa_x, tr_z_xxxx_xxx, tr_z_xxxx_xxy, tr_z_xxxx_xxz, tr_z_xxxx_xyy, tr_z_xxxx_xyz, tr_z_xxxx_xzz, tr_z_xxxx_yyy, tr_z_xxxx_yyz, tr_z_xxxx_yzz, tr_z_xxxx_zzz, tr_z_xxxxx_xx, tr_z_xxxxx_xxx, tr_z_xxxxx_xxy, tr_z_xxxxx_xxz, tr_z_xxxxx_xy, tr_z_xxxxx_xyy, tr_z_xxxxx_xyz, tr_z_xxxxx_xz, tr_z_xxxxx_xzz, tr_z_xxxxx_yy, tr_z_xxxxx_yyy, tr_z_xxxxx_yyz, tr_z_xxxxx_yz, tr_z_xxxxx_yzz, tr_z_xxxxx_zz, tr_z_xxxxx_zzz, tr_z_xxxxxx_xxx, tr_z_xxxxxx_xxy, tr_z_xxxxxx_xxz, tr_z_xxxxxx_xyy, tr_z_xxxxxx_xyz, tr_z_xxxxxx_xzz, tr_z_xxxxxx_yyy, tr_z_xxxxxx_yyz, tr_z_xxxxxx_yzz, tr_z_xxxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxx_xxx[i] = 5.0 * tr_z_xxxx_xxx[i] * fe_0 + 3.0 * tr_z_xxxxx_xx[i] * fe_0 + tr_z_xxxxx_xxx[i] * pa_x[i];

        tr_z_xxxxxx_xxy[i] = 5.0 * tr_z_xxxx_xxy[i] * fe_0 + 2.0 * tr_z_xxxxx_xy[i] * fe_0 + tr_z_xxxxx_xxy[i] * pa_x[i];

        tr_z_xxxxxx_xxz[i] = 5.0 * tr_z_xxxx_xxz[i] * fe_0 + 2.0 * tr_z_xxxxx_xz[i] * fe_0 + tr_z_xxxxx_xxz[i] * pa_x[i];

        tr_z_xxxxxx_xyy[i] = 5.0 * tr_z_xxxx_xyy[i] * fe_0 + tr_z_xxxxx_yy[i] * fe_0 + tr_z_xxxxx_xyy[i] * pa_x[i];

        tr_z_xxxxxx_xyz[i] = 5.0 * tr_z_xxxx_xyz[i] * fe_0 + tr_z_xxxxx_yz[i] * fe_0 + tr_z_xxxxx_xyz[i] * pa_x[i];

        tr_z_xxxxxx_xzz[i] = 5.0 * tr_z_xxxx_xzz[i] * fe_0 + tr_z_xxxxx_zz[i] * fe_0 + tr_z_xxxxx_xzz[i] * pa_x[i];

        tr_z_xxxxxx_yyy[i] = 5.0 * tr_z_xxxx_yyy[i] * fe_0 + tr_z_xxxxx_yyy[i] * pa_x[i];

        tr_z_xxxxxx_yyz[i] = 5.0 * tr_z_xxxx_yyz[i] * fe_0 + tr_z_xxxxx_yyz[i] * pa_x[i];

        tr_z_xxxxxx_yzz[i] = 5.0 * tr_z_xxxx_yzz[i] * fe_0 + tr_z_xxxxx_yzz[i] * pa_x[i];

        tr_z_xxxxxx_zzz[i] = 5.0 * tr_z_xxxx_zzz[i] * fe_0 + tr_z_xxxxx_zzz[i] * pa_x[i];
    }

    // Set up 570-580 components of targeted buffer : IF

    auto tr_z_xxxxxy_xxx = pbuffer.data(idx_dip_if + 570);

    auto tr_z_xxxxxy_xxy = pbuffer.data(idx_dip_if + 571);

    auto tr_z_xxxxxy_xxz = pbuffer.data(idx_dip_if + 572);

    auto tr_z_xxxxxy_xyy = pbuffer.data(idx_dip_if + 573);

    auto tr_z_xxxxxy_xyz = pbuffer.data(idx_dip_if + 574);

    auto tr_z_xxxxxy_xzz = pbuffer.data(idx_dip_if + 575);

    auto tr_z_xxxxxy_yyy = pbuffer.data(idx_dip_if + 576);

    auto tr_z_xxxxxy_yyz = pbuffer.data(idx_dip_if + 577);

    auto tr_z_xxxxxy_yzz = pbuffer.data(idx_dip_if + 578);

    auto tr_z_xxxxxy_zzz = pbuffer.data(idx_dip_if + 579);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxx_xx, tr_z_xxxxx_xxx, tr_z_xxxxx_xxy, tr_z_xxxxx_xxz, tr_z_xxxxx_xy, tr_z_xxxxx_xyy, tr_z_xxxxx_xyz, tr_z_xxxxx_xz, tr_z_xxxxx_xzz, tr_z_xxxxx_zzz, tr_z_xxxxxy_xxx, tr_z_xxxxxy_xxy, tr_z_xxxxxy_xxz, tr_z_xxxxxy_xyy, tr_z_xxxxxy_xyz, tr_z_xxxxxy_xzz, tr_z_xxxxxy_yyy, tr_z_xxxxxy_yyz, tr_z_xxxxxy_yzz, tr_z_xxxxxy_zzz, tr_z_xxxxy_yyy, tr_z_xxxxy_yyz, tr_z_xxxxy_yzz, tr_z_xxxy_yyy, tr_z_xxxy_yyz, tr_z_xxxy_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxy_xxx[i] = tr_z_xxxxx_xxx[i] * pa_y[i];

        tr_z_xxxxxy_xxy[i] = tr_z_xxxxx_xx[i] * fe_0 + tr_z_xxxxx_xxy[i] * pa_y[i];

        tr_z_xxxxxy_xxz[i] = tr_z_xxxxx_xxz[i] * pa_y[i];

        tr_z_xxxxxy_xyy[i] = 2.0 * tr_z_xxxxx_xy[i] * fe_0 + tr_z_xxxxx_xyy[i] * pa_y[i];

        tr_z_xxxxxy_xyz[i] = tr_z_xxxxx_xz[i] * fe_0 + tr_z_xxxxx_xyz[i] * pa_y[i];

        tr_z_xxxxxy_xzz[i] = tr_z_xxxxx_xzz[i] * pa_y[i];

        tr_z_xxxxxy_yyy[i] = 4.0 * tr_z_xxxy_yyy[i] * fe_0 + tr_z_xxxxy_yyy[i] * pa_x[i];

        tr_z_xxxxxy_yyz[i] = 4.0 * tr_z_xxxy_yyz[i] * fe_0 + tr_z_xxxxy_yyz[i] * pa_x[i];

        tr_z_xxxxxy_yzz[i] = 4.0 * tr_z_xxxy_yzz[i] * fe_0 + tr_z_xxxxy_yzz[i] * pa_x[i];

        tr_z_xxxxxy_zzz[i] = tr_z_xxxxx_zzz[i] * pa_y[i];
    }

    // Set up 580-590 components of targeted buffer : IF

    auto tr_z_xxxxxz_xxx = pbuffer.data(idx_dip_if + 580);

    auto tr_z_xxxxxz_xxy = pbuffer.data(idx_dip_if + 581);

    auto tr_z_xxxxxz_xxz = pbuffer.data(idx_dip_if + 582);

    auto tr_z_xxxxxz_xyy = pbuffer.data(idx_dip_if + 583);

    auto tr_z_xxxxxz_xyz = pbuffer.data(idx_dip_if + 584);

    auto tr_z_xxxxxz_xzz = pbuffer.data(idx_dip_if + 585);

    auto tr_z_xxxxxz_yyy = pbuffer.data(idx_dip_if + 586);

    auto tr_z_xxxxxz_yyz = pbuffer.data(idx_dip_if + 587);

    auto tr_z_xxxxxz_yzz = pbuffer.data(idx_dip_if + 588);

    auto tr_z_xxxxxz_zzz = pbuffer.data(idx_dip_if + 589);

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xxxxx_xxx, tr_z_xxxxx_xxy, tr_z_xxxxx_xyy, tr_z_xxxxxz_xxx, tr_z_xxxxxz_xxy, tr_z_xxxxxz_xxz, tr_z_xxxxxz_xyy, tr_z_xxxxxz_xyz, tr_z_xxxxxz_xzz, tr_z_xxxxxz_yyy, tr_z_xxxxxz_yyz, tr_z_xxxxxz_yzz, tr_z_xxxxxz_zzz, tr_z_xxxxz_xxz, tr_z_xxxxz_xyz, tr_z_xxxxz_xz, tr_z_xxxxz_xzz, tr_z_xxxxz_yyy, tr_z_xxxxz_yyz, tr_z_xxxxz_yz, tr_z_xxxxz_yzz, tr_z_xxxxz_zz, tr_z_xxxxz_zzz, tr_z_xxxz_xxz, tr_z_xxxz_xyz, tr_z_xxxz_xzz, tr_z_xxxz_yyy, tr_z_xxxz_yyz, tr_z_xxxz_yzz, tr_z_xxxz_zzz, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxz_xxx[i] = ts_xxxxx_xxx[i] * fe_0 + tr_z_xxxxx_xxx[i] * pa_z[i];

        tr_z_xxxxxz_xxy[i] = ts_xxxxx_xxy[i] * fe_0 + tr_z_xxxxx_xxy[i] * pa_z[i];

        tr_z_xxxxxz_xxz[i] = 4.0 * tr_z_xxxz_xxz[i] * fe_0 + 2.0 * tr_z_xxxxz_xz[i] * fe_0 + tr_z_xxxxz_xxz[i] * pa_x[i];

        tr_z_xxxxxz_xyy[i] = ts_xxxxx_xyy[i] * fe_0 + tr_z_xxxxx_xyy[i] * pa_z[i];

        tr_z_xxxxxz_xyz[i] = 4.0 * tr_z_xxxz_xyz[i] * fe_0 + tr_z_xxxxz_yz[i] * fe_0 + tr_z_xxxxz_xyz[i] * pa_x[i];

        tr_z_xxxxxz_xzz[i] = 4.0 * tr_z_xxxz_xzz[i] * fe_0 + tr_z_xxxxz_zz[i] * fe_0 + tr_z_xxxxz_xzz[i] * pa_x[i];

        tr_z_xxxxxz_yyy[i] = 4.0 * tr_z_xxxz_yyy[i] * fe_0 + tr_z_xxxxz_yyy[i] * pa_x[i];

        tr_z_xxxxxz_yyz[i] = 4.0 * tr_z_xxxz_yyz[i] * fe_0 + tr_z_xxxxz_yyz[i] * pa_x[i];

        tr_z_xxxxxz_yzz[i] = 4.0 * tr_z_xxxz_yzz[i] * fe_0 + tr_z_xxxxz_yzz[i] * pa_x[i];

        tr_z_xxxxxz_zzz[i] = 4.0 * tr_z_xxxz_zzz[i] * fe_0 + tr_z_xxxxz_zzz[i] * pa_x[i];
    }

    // Set up 590-600 components of targeted buffer : IF

    auto tr_z_xxxxyy_xxx = pbuffer.data(idx_dip_if + 590);

    auto tr_z_xxxxyy_xxy = pbuffer.data(idx_dip_if + 591);

    auto tr_z_xxxxyy_xxz = pbuffer.data(idx_dip_if + 592);

    auto tr_z_xxxxyy_xyy = pbuffer.data(idx_dip_if + 593);

    auto tr_z_xxxxyy_xyz = pbuffer.data(idx_dip_if + 594);

    auto tr_z_xxxxyy_xzz = pbuffer.data(idx_dip_if + 595);

    auto tr_z_xxxxyy_yyy = pbuffer.data(idx_dip_if + 596);

    auto tr_z_xxxxyy_yyz = pbuffer.data(idx_dip_if + 597);

    auto tr_z_xxxxyy_yzz = pbuffer.data(idx_dip_if + 598);

    auto tr_z_xxxxyy_zzz = pbuffer.data(idx_dip_if + 599);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxx_xxx, tr_z_xxxx_xxz, tr_z_xxxx_xzz, tr_z_xxxxy_xxx, tr_z_xxxxy_xxz, tr_z_xxxxy_xzz, tr_z_xxxxyy_xxx, tr_z_xxxxyy_xxy, tr_z_xxxxyy_xxz, tr_z_xxxxyy_xyy, tr_z_xxxxyy_xyz, tr_z_xxxxyy_xzz, tr_z_xxxxyy_yyy, tr_z_xxxxyy_yyz, tr_z_xxxxyy_yzz, tr_z_xxxxyy_zzz, tr_z_xxxyy_xxy, tr_z_xxxyy_xy, tr_z_xxxyy_xyy, tr_z_xxxyy_xyz, tr_z_xxxyy_yy, tr_z_xxxyy_yyy, tr_z_xxxyy_yyz, tr_z_xxxyy_yz, tr_z_xxxyy_yzz, tr_z_xxxyy_zzz, tr_z_xxyy_xxy, tr_z_xxyy_xyy, tr_z_xxyy_xyz, tr_z_xxyy_yyy, tr_z_xxyy_yyz, tr_z_xxyy_yzz, tr_z_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyy_xxx[i] = tr_z_xxxx_xxx[i] * fe_0 + tr_z_xxxxy_xxx[i] * pa_y[i];

        tr_z_xxxxyy_xxy[i] = 3.0 * tr_z_xxyy_xxy[i] * fe_0 + 2.0 * tr_z_xxxyy_xy[i] * fe_0 + tr_z_xxxyy_xxy[i] * pa_x[i];

        tr_z_xxxxyy_xxz[i] = tr_z_xxxx_xxz[i] * fe_0 + tr_z_xxxxy_xxz[i] * pa_y[i];

        tr_z_xxxxyy_xyy[i] = 3.0 * tr_z_xxyy_xyy[i] * fe_0 + tr_z_xxxyy_yy[i] * fe_0 + tr_z_xxxyy_xyy[i] * pa_x[i];

        tr_z_xxxxyy_xyz[i] = 3.0 * tr_z_xxyy_xyz[i] * fe_0 + tr_z_xxxyy_yz[i] * fe_0 + tr_z_xxxyy_xyz[i] * pa_x[i];

        tr_z_xxxxyy_xzz[i] = tr_z_xxxx_xzz[i] * fe_0 + tr_z_xxxxy_xzz[i] * pa_y[i];

        tr_z_xxxxyy_yyy[i] = 3.0 * tr_z_xxyy_yyy[i] * fe_0 + tr_z_xxxyy_yyy[i] * pa_x[i];

        tr_z_xxxxyy_yyz[i] = 3.0 * tr_z_xxyy_yyz[i] * fe_0 + tr_z_xxxyy_yyz[i] * pa_x[i];

        tr_z_xxxxyy_yzz[i] = 3.0 * tr_z_xxyy_yzz[i] * fe_0 + tr_z_xxxyy_yzz[i] * pa_x[i];

        tr_z_xxxxyy_zzz[i] = 3.0 * tr_z_xxyy_zzz[i] * fe_0 + tr_z_xxxyy_zzz[i] * pa_x[i];
    }

    // Set up 600-610 components of targeted buffer : IF

    auto tr_z_xxxxyz_xxx = pbuffer.data(idx_dip_if + 600);

    auto tr_z_xxxxyz_xxy = pbuffer.data(idx_dip_if + 601);

    auto tr_z_xxxxyz_xxz = pbuffer.data(idx_dip_if + 602);

    auto tr_z_xxxxyz_xyy = pbuffer.data(idx_dip_if + 603);

    auto tr_z_xxxxyz_xyz = pbuffer.data(idx_dip_if + 604);

    auto tr_z_xxxxyz_xzz = pbuffer.data(idx_dip_if + 605);

    auto tr_z_xxxxyz_yyy = pbuffer.data(idx_dip_if + 606);

    auto tr_z_xxxxyz_yyz = pbuffer.data(idx_dip_if + 607);

    auto tr_z_xxxxyz_yzz = pbuffer.data(idx_dip_if + 608);

    auto tr_z_xxxxyz_zzz = pbuffer.data(idx_dip_if + 609);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxyz_xxx, tr_z_xxxxyz_xxy, tr_z_xxxxyz_xxz, tr_z_xxxxyz_xyy, tr_z_xxxxyz_xyz, tr_z_xxxxyz_xzz, tr_z_xxxxyz_yyy, tr_z_xxxxyz_yyz, tr_z_xxxxyz_yzz, tr_z_xxxxyz_zzz, tr_z_xxxxz_xx, tr_z_xxxxz_xxx, tr_z_xxxxz_xxy, tr_z_xxxxz_xxz, tr_z_xxxxz_xy, tr_z_xxxxz_xyy, tr_z_xxxxz_xyz, tr_z_xxxxz_xz, tr_z_xxxxz_xzz, tr_z_xxxxz_zzz, tr_z_xxxyz_yyy, tr_z_xxxyz_yyz, tr_z_xxxyz_yzz, tr_z_xxyz_yyy, tr_z_xxyz_yyz, tr_z_xxyz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyz_xxx[i] = tr_z_xxxxz_xxx[i] * pa_y[i];

        tr_z_xxxxyz_xxy[i] = tr_z_xxxxz_xx[i] * fe_0 + tr_z_xxxxz_xxy[i] * pa_y[i];

        tr_z_xxxxyz_xxz[i] = tr_z_xxxxz_xxz[i] * pa_y[i];

        tr_z_xxxxyz_xyy[i] = 2.0 * tr_z_xxxxz_xy[i] * fe_0 + tr_z_xxxxz_xyy[i] * pa_y[i];

        tr_z_xxxxyz_xyz[i] = tr_z_xxxxz_xz[i] * fe_0 + tr_z_xxxxz_xyz[i] * pa_y[i];

        tr_z_xxxxyz_xzz[i] = tr_z_xxxxz_xzz[i] * pa_y[i];

        tr_z_xxxxyz_yyy[i] = 3.0 * tr_z_xxyz_yyy[i] * fe_0 + tr_z_xxxyz_yyy[i] * pa_x[i];

        tr_z_xxxxyz_yyz[i] = 3.0 * tr_z_xxyz_yyz[i] * fe_0 + tr_z_xxxyz_yyz[i] * pa_x[i];

        tr_z_xxxxyz_yzz[i] = 3.0 * tr_z_xxyz_yzz[i] * fe_0 + tr_z_xxxyz_yzz[i] * pa_x[i];

        tr_z_xxxxyz_zzz[i] = tr_z_xxxxz_zzz[i] * pa_y[i];
    }

    // Set up 610-620 components of targeted buffer : IF

    auto tr_z_xxxxzz_xxx = pbuffer.data(idx_dip_if + 610);

    auto tr_z_xxxxzz_xxy = pbuffer.data(idx_dip_if + 611);

    auto tr_z_xxxxzz_xxz = pbuffer.data(idx_dip_if + 612);

    auto tr_z_xxxxzz_xyy = pbuffer.data(idx_dip_if + 613);

    auto tr_z_xxxxzz_xyz = pbuffer.data(idx_dip_if + 614);

    auto tr_z_xxxxzz_xzz = pbuffer.data(idx_dip_if + 615);

    auto tr_z_xxxxzz_yyy = pbuffer.data(idx_dip_if + 616);

    auto tr_z_xxxxzz_yyz = pbuffer.data(idx_dip_if + 617);

    auto tr_z_xxxxzz_yzz = pbuffer.data(idx_dip_if + 618);

    auto tr_z_xxxxzz_zzz = pbuffer.data(idx_dip_if + 619);

    #pragma omp simd aligned(pa_x, tr_z_xxxxzz_xxx, tr_z_xxxxzz_xxy, tr_z_xxxxzz_xxz, tr_z_xxxxzz_xyy, tr_z_xxxxzz_xyz, tr_z_xxxxzz_xzz, tr_z_xxxxzz_yyy, tr_z_xxxxzz_yyz, tr_z_xxxxzz_yzz, tr_z_xxxxzz_zzz, tr_z_xxxzz_xx, tr_z_xxxzz_xxx, tr_z_xxxzz_xxy, tr_z_xxxzz_xxz, tr_z_xxxzz_xy, tr_z_xxxzz_xyy, tr_z_xxxzz_xyz, tr_z_xxxzz_xz, tr_z_xxxzz_xzz, tr_z_xxxzz_yy, tr_z_xxxzz_yyy, tr_z_xxxzz_yyz, tr_z_xxxzz_yz, tr_z_xxxzz_yzz, tr_z_xxxzz_zz, tr_z_xxxzz_zzz, tr_z_xxzz_xxx, tr_z_xxzz_xxy, tr_z_xxzz_xxz, tr_z_xxzz_xyy, tr_z_xxzz_xyz, tr_z_xxzz_xzz, tr_z_xxzz_yyy, tr_z_xxzz_yyz, tr_z_xxzz_yzz, tr_z_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxzz_xxx[i] = 3.0 * tr_z_xxzz_xxx[i] * fe_0 + 3.0 * tr_z_xxxzz_xx[i] * fe_0 + tr_z_xxxzz_xxx[i] * pa_x[i];

        tr_z_xxxxzz_xxy[i] = 3.0 * tr_z_xxzz_xxy[i] * fe_0 + 2.0 * tr_z_xxxzz_xy[i] * fe_0 + tr_z_xxxzz_xxy[i] * pa_x[i];

        tr_z_xxxxzz_xxz[i] = 3.0 * tr_z_xxzz_xxz[i] * fe_0 + 2.0 * tr_z_xxxzz_xz[i] * fe_0 + tr_z_xxxzz_xxz[i] * pa_x[i];

        tr_z_xxxxzz_xyy[i] = 3.0 * tr_z_xxzz_xyy[i] * fe_0 + tr_z_xxxzz_yy[i] * fe_0 + tr_z_xxxzz_xyy[i] * pa_x[i];

        tr_z_xxxxzz_xyz[i] = 3.0 * tr_z_xxzz_xyz[i] * fe_0 + tr_z_xxxzz_yz[i] * fe_0 + tr_z_xxxzz_xyz[i] * pa_x[i];

        tr_z_xxxxzz_xzz[i] = 3.0 * tr_z_xxzz_xzz[i] * fe_0 + tr_z_xxxzz_zz[i] * fe_0 + tr_z_xxxzz_xzz[i] * pa_x[i];

        tr_z_xxxxzz_yyy[i] = 3.0 * tr_z_xxzz_yyy[i] * fe_0 + tr_z_xxxzz_yyy[i] * pa_x[i];

        tr_z_xxxxzz_yyz[i] = 3.0 * tr_z_xxzz_yyz[i] * fe_0 + tr_z_xxxzz_yyz[i] * pa_x[i];

        tr_z_xxxxzz_yzz[i] = 3.0 * tr_z_xxzz_yzz[i] * fe_0 + tr_z_xxxzz_yzz[i] * pa_x[i];

        tr_z_xxxxzz_zzz[i] = 3.0 * tr_z_xxzz_zzz[i] * fe_0 + tr_z_xxxzz_zzz[i] * pa_x[i];
    }

    // Set up 620-630 components of targeted buffer : IF

    auto tr_z_xxxyyy_xxx = pbuffer.data(idx_dip_if + 620);

    auto tr_z_xxxyyy_xxy = pbuffer.data(idx_dip_if + 621);

    auto tr_z_xxxyyy_xxz = pbuffer.data(idx_dip_if + 622);

    auto tr_z_xxxyyy_xyy = pbuffer.data(idx_dip_if + 623);

    auto tr_z_xxxyyy_xyz = pbuffer.data(idx_dip_if + 624);

    auto tr_z_xxxyyy_xzz = pbuffer.data(idx_dip_if + 625);

    auto tr_z_xxxyyy_yyy = pbuffer.data(idx_dip_if + 626);

    auto tr_z_xxxyyy_yyz = pbuffer.data(idx_dip_if + 627);

    auto tr_z_xxxyyy_yzz = pbuffer.data(idx_dip_if + 628);

    auto tr_z_xxxyyy_zzz = pbuffer.data(idx_dip_if + 629);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxy_xxx, tr_z_xxxy_xxz, tr_z_xxxy_xzz, tr_z_xxxyy_xxx, tr_z_xxxyy_xxz, tr_z_xxxyy_xzz, tr_z_xxxyyy_xxx, tr_z_xxxyyy_xxy, tr_z_xxxyyy_xxz, tr_z_xxxyyy_xyy, tr_z_xxxyyy_xyz, tr_z_xxxyyy_xzz, tr_z_xxxyyy_yyy, tr_z_xxxyyy_yyz, tr_z_xxxyyy_yzz, tr_z_xxxyyy_zzz, tr_z_xxyyy_xxy, tr_z_xxyyy_xy, tr_z_xxyyy_xyy, tr_z_xxyyy_xyz, tr_z_xxyyy_yy, tr_z_xxyyy_yyy, tr_z_xxyyy_yyz, tr_z_xxyyy_yz, tr_z_xxyyy_yzz, tr_z_xxyyy_zzz, tr_z_xyyy_xxy, tr_z_xyyy_xyy, tr_z_xyyy_xyz, tr_z_xyyy_yyy, tr_z_xyyy_yyz, tr_z_xyyy_yzz, tr_z_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyy_xxx[i] = 2.0 * tr_z_xxxy_xxx[i] * fe_0 + tr_z_xxxyy_xxx[i] * pa_y[i];

        tr_z_xxxyyy_xxy[i] = 2.0 * tr_z_xyyy_xxy[i] * fe_0 + 2.0 * tr_z_xxyyy_xy[i] * fe_0 + tr_z_xxyyy_xxy[i] * pa_x[i];

        tr_z_xxxyyy_xxz[i] = 2.0 * tr_z_xxxy_xxz[i] * fe_0 + tr_z_xxxyy_xxz[i] * pa_y[i];

        tr_z_xxxyyy_xyy[i] = 2.0 * tr_z_xyyy_xyy[i] * fe_0 + tr_z_xxyyy_yy[i] * fe_0 + tr_z_xxyyy_xyy[i] * pa_x[i];

        tr_z_xxxyyy_xyz[i] = 2.0 * tr_z_xyyy_xyz[i] * fe_0 + tr_z_xxyyy_yz[i] * fe_0 + tr_z_xxyyy_xyz[i] * pa_x[i];

        tr_z_xxxyyy_xzz[i] = 2.0 * tr_z_xxxy_xzz[i] * fe_0 + tr_z_xxxyy_xzz[i] * pa_y[i];

        tr_z_xxxyyy_yyy[i] = 2.0 * tr_z_xyyy_yyy[i] * fe_0 + tr_z_xxyyy_yyy[i] * pa_x[i];

        tr_z_xxxyyy_yyz[i] = 2.0 * tr_z_xyyy_yyz[i] * fe_0 + tr_z_xxyyy_yyz[i] * pa_x[i];

        tr_z_xxxyyy_yzz[i] = 2.0 * tr_z_xyyy_yzz[i] * fe_0 + tr_z_xxyyy_yzz[i] * pa_x[i];

        tr_z_xxxyyy_zzz[i] = 2.0 * tr_z_xyyy_zzz[i] * fe_0 + tr_z_xxyyy_zzz[i] * pa_x[i];
    }

    // Set up 630-640 components of targeted buffer : IF

    auto tr_z_xxxyyz_xxx = pbuffer.data(idx_dip_if + 630);

    auto tr_z_xxxyyz_xxy = pbuffer.data(idx_dip_if + 631);

    auto tr_z_xxxyyz_xxz = pbuffer.data(idx_dip_if + 632);

    auto tr_z_xxxyyz_xyy = pbuffer.data(idx_dip_if + 633);

    auto tr_z_xxxyyz_xyz = pbuffer.data(idx_dip_if + 634);

    auto tr_z_xxxyyz_xzz = pbuffer.data(idx_dip_if + 635);

    auto tr_z_xxxyyz_yyy = pbuffer.data(idx_dip_if + 636);

    auto tr_z_xxxyyz_yyz = pbuffer.data(idx_dip_if + 637);

    auto tr_z_xxxyyz_yzz = pbuffer.data(idx_dip_if + 638);

    auto tr_z_xxxyyz_zzz = pbuffer.data(idx_dip_if + 639);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxxyy_xxy, tr_z_xxxyy_xyy, tr_z_xxxyyz_xxx, tr_z_xxxyyz_xxy, tr_z_xxxyyz_xxz, tr_z_xxxyyz_xyy, tr_z_xxxyyz_xyz, tr_z_xxxyyz_xzz, tr_z_xxxyyz_yyy, tr_z_xxxyyz_yyz, tr_z_xxxyyz_yzz, tr_z_xxxyyz_zzz, tr_z_xxxyz_xxx, tr_z_xxxyz_xxz, tr_z_xxxyz_xzz, tr_z_xxxz_xxx, tr_z_xxxz_xxz, tr_z_xxxz_xzz, tr_z_xxyyz_xyz, tr_z_xxyyz_yyy, tr_z_xxyyz_yyz, tr_z_xxyyz_yz, tr_z_xxyyz_yzz, tr_z_xxyyz_zzz, tr_z_xyyz_xyz, tr_z_xyyz_yyy, tr_z_xyyz_yyz, tr_z_xyyz_yzz, tr_z_xyyz_zzz, ts_xxxyy_xxy, ts_xxxyy_xyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyz_xxx[i] = tr_z_xxxz_xxx[i] * fe_0 + tr_z_xxxyz_xxx[i] * pa_y[i];

        tr_z_xxxyyz_xxy[i] = ts_xxxyy_xxy[i] * fe_0 + tr_z_xxxyy_xxy[i] * pa_z[i];

        tr_z_xxxyyz_xxz[i] = tr_z_xxxz_xxz[i] * fe_0 + tr_z_xxxyz_xxz[i] * pa_y[i];

        tr_z_xxxyyz_xyy[i] = ts_xxxyy_xyy[i] * fe_0 + tr_z_xxxyy_xyy[i] * pa_z[i];

        tr_z_xxxyyz_xyz[i] = 2.0 * tr_z_xyyz_xyz[i] * fe_0 + tr_z_xxyyz_yz[i] * fe_0 + tr_z_xxyyz_xyz[i] * pa_x[i];

        tr_z_xxxyyz_xzz[i] = tr_z_xxxz_xzz[i] * fe_0 + tr_z_xxxyz_xzz[i] * pa_y[i];

        tr_z_xxxyyz_yyy[i] = 2.0 * tr_z_xyyz_yyy[i] * fe_0 + tr_z_xxyyz_yyy[i] * pa_x[i];

        tr_z_xxxyyz_yyz[i] = 2.0 * tr_z_xyyz_yyz[i] * fe_0 + tr_z_xxyyz_yyz[i] * pa_x[i];

        tr_z_xxxyyz_yzz[i] = 2.0 * tr_z_xyyz_yzz[i] * fe_0 + tr_z_xxyyz_yzz[i] * pa_x[i];

        tr_z_xxxyyz_zzz[i] = 2.0 * tr_z_xyyz_zzz[i] * fe_0 + tr_z_xxyyz_zzz[i] * pa_x[i];
    }

    // Set up 640-650 components of targeted buffer : IF

    auto tr_z_xxxyzz_xxx = pbuffer.data(idx_dip_if + 640);

    auto tr_z_xxxyzz_xxy = pbuffer.data(idx_dip_if + 641);

    auto tr_z_xxxyzz_xxz = pbuffer.data(idx_dip_if + 642);

    auto tr_z_xxxyzz_xyy = pbuffer.data(idx_dip_if + 643);

    auto tr_z_xxxyzz_xyz = pbuffer.data(idx_dip_if + 644);

    auto tr_z_xxxyzz_xzz = pbuffer.data(idx_dip_if + 645);

    auto tr_z_xxxyzz_yyy = pbuffer.data(idx_dip_if + 646);

    auto tr_z_xxxyzz_yyz = pbuffer.data(idx_dip_if + 647);

    auto tr_z_xxxyzz_yzz = pbuffer.data(idx_dip_if + 648);

    auto tr_z_xxxyzz_zzz = pbuffer.data(idx_dip_if + 649);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxyzz_xxx, tr_z_xxxyzz_xxy, tr_z_xxxyzz_xxz, tr_z_xxxyzz_xyy, tr_z_xxxyzz_xyz, tr_z_xxxyzz_xzz, tr_z_xxxyzz_yyy, tr_z_xxxyzz_yyz, tr_z_xxxyzz_yzz, tr_z_xxxyzz_zzz, tr_z_xxxzz_xx, tr_z_xxxzz_xxx, tr_z_xxxzz_xxy, tr_z_xxxzz_xxz, tr_z_xxxzz_xy, tr_z_xxxzz_xyy, tr_z_xxxzz_xyz, tr_z_xxxzz_xz, tr_z_xxxzz_xzz, tr_z_xxxzz_zzz, tr_z_xxyzz_yyy, tr_z_xxyzz_yyz, tr_z_xxyzz_yzz, tr_z_xyzz_yyy, tr_z_xyzz_yyz, tr_z_xyzz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyzz_xxx[i] = tr_z_xxxzz_xxx[i] * pa_y[i];

        tr_z_xxxyzz_xxy[i] = tr_z_xxxzz_xx[i] * fe_0 + tr_z_xxxzz_xxy[i] * pa_y[i];

        tr_z_xxxyzz_xxz[i] = tr_z_xxxzz_xxz[i] * pa_y[i];

        tr_z_xxxyzz_xyy[i] = 2.0 * tr_z_xxxzz_xy[i] * fe_0 + tr_z_xxxzz_xyy[i] * pa_y[i];

        tr_z_xxxyzz_xyz[i] = tr_z_xxxzz_xz[i] * fe_0 + tr_z_xxxzz_xyz[i] * pa_y[i];

        tr_z_xxxyzz_xzz[i] = tr_z_xxxzz_xzz[i] * pa_y[i];

        tr_z_xxxyzz_yyy[i] = 2.0 * tr_z_xyzz_yyy[i] * fe_0 + tr_z_xxyzz_yyy[i] * pa_x[i];

        tr_z_xxxyzz_yyz[i] = 2.0 * tr_z_xyzz_yyz[i] * fe_0 + tr_z_xxyzz_yyz[i] * pa_x[i];

        tr_z_xxxyzz_yzz[i] = 2.0 * tr_z_xyzz_yzz[i] * fe_0 + tr_z_xxyzz_yzz[i] * pa_x[i];

        tr_z_xxxyzz_zzz[i] = tr_z_xxxzz_zzz[i] * pa_y[i];
    }

    // Set up 650-660 components of targeted buffer : IF

    auto tr_z_xxxzzz_xxx = pbuffer.data(idx_dip_if + 650);

    auto tr_z_xxxzzz_xxy = pbuffer.data(idx_dip_if + 651);

    auto tr_z_xxxzzz_xxz = pbuffer.data(idx_dip_if + 652);

    auto tr_z_xxxzzz_xyy = pbuffer.data(idx_dip_if + 653);

    auto tr_z_xxxzzz_xyz = pbuffer.data(idx_dip_if + 654);

    auto tr_z_xxxzzz_xzz = pbuffer.data(idx_dip_if + 655);

    auto tr_z_xxxzzz_yyy = pbuffer.data(idx_dip_if + 656);

    auto tr_z_xxxzzz_yyz = pbuffer.data(idx_dip_if + 657);

    auto tr_z_xxxzzz_yzz = pbuffer.data(idx_dip_if + 658);

    auto tr_z_xxxzzz_zzz = pbuffer.data(idx_dip_if + 659);

    #pragma omp simd aligned(pa_x, tr_z_xxxzzz_xxx, tr_z_xxxzzz_xxy, tr_z_xxxzzz_xxz, tr_z_xxxzzz_xyy, tr_z_xxxzzz_xyz, tr_z_xxxzzz_xzz, tr_z_xxxzzz_yyy, tr_z_xxxzzz_yyz, tr_z_xxxzzz_yzz, tr_z_xxxzzz_zzz, tr_z_xxzzz_xx, tr_z_xxzzz_xxx, tr_z_xxzzz_xxy, tr_z_xxzzz_xxz, tr_z_xxzzz_xy, tr_z_xxzzz_xyy, tr_z_xxzzz_xyz, tr_z_xxzzz_xz, tr_z_xxzzz_xzz, tr_z_xxzzz_yy, tr_z_xxzzz_yyy, tr_z_xxzzz_yyz, tr_z_xxzzz_yz, tr_z_xxzzz_yzz, tr_z_xxzzz_zz, tr_z_xxzzz_zzz, tr_z_xzzz_xxx, tr_z_xzzz_xxy, tr_z_xzzz_xxz, tr_z_xzzz_xyy, tr_z_xzzz_xyz, tr_z_xzzz_xzz, tr_z_xzzz_yyy, tr_z_xzzz_yyz, tr_z_xzzz_yzz, tr_z_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzzz_xxx[i] = 2.0 * tr_z_xzzz_xxx[i] * fe_0 + 3.0 * tr_z_xxzzz_xx[i] * fe_0 + tr_z_xxzzz_xxx[i] * pa_x[i];

        tr_z_xxxzzz_xxy[i] = 2.0 * tr_z_xzzz_xxy[i] * fe_0 + 2.0 * tr_z_xxzzz_xy[i] * fe_0 + tr_z_xxzzz_xxy[i] * pa_x[i];

        tr_z_xxxzzz_xxz[i] = 2.0 * tr_z_xzzz_xxz[i] * fe_0 + 2.0 * tr_z_xxzzz_xz[i] * fe_0 + tr_z_xxzzz_xxz[i] * pa_x[i];

        tr_z_xxxzzz_xyy[i] = 2.0 * tr_z_xzzz_xyy[i] * fe_0 + tr_z_xxzzz_yy[i] * fe_0 + tr_z_xxzzz_xyy[i] * pa_x[i];

        tr_z_xxxzzz_xyz[i] = 2.0 * tr_z_xzzz_xyz[i] * fe_0 + tr_z_xxzzz_yz[i] * fe_0 + tr_z_xxzzz_xyz[i] * pa_x[i];

        tr_z_xxxzzz_xzz[i] = 2.0 * tr_z_xzzz_xzz[i] * fe_0 + tr_z_xxzzz_zz[i] * fe_0 + tr_z_xxzzz_xzz[i] * pa_x[i];

        tr_z_xxxzzz_yyy[i] = 2.0 * tr_z_xzzz_yyy[i] * fe_0 + tr_z_xxzzz_yyy[i] * pa_x[i];

        tr_z_xxxzzz_yyz[i] = 2.0 * tr_z_xzzz_yyz[i] * fe_0 + tr_z_xxzzz_yyz[i] * pa_x[i];

        tr_z_xxxzzz_yzz[i] = 2.0 * tr_z_xzzz_yzz[i] * fe_0 + tr_z_xxzzz_yzz[i] * pa_x[i];

        tr_z_xxxzzz_zzz[i] = 2.0 * tr_z_xzzz_zzz[i] * fe_0 + tr_z_xxzzz_zzz[i] * pa_x[i];
    }

    // Set up 660-670 components of targeted buffer : IF

    auto tr_z_xxyyyy_xxx = pbuffer.data(idx_dip_if + 660);

    auto tr_z_xxyyyy_xxy = pbuffer.data(idx_dip_if + 661);

    auto tr_z_xxyyyy_xxz = pbuffer.data(idx_dip_if + 662);

    auto tr_z_xxyyyy_xyy = pbuffer.data(idx_dip_if + 663);

    auto tr_z_xxyyyy_xyz = pbuffer.data(idx_dip_if + 664);

    auto tr_z_xxyyyy_xzz = pbuffer.data(idx_dip_if + 665);

    auto tr_z_xxyyyy_yyy = pbuffer.data(idx_dip_if + 666);

    auto tr_z_xxyyyy_yyz = pbuffer.data(idx_dip_if + 667);

    auto tr_z_xxyyyy_yzz = pbuffer.data(idx_dip_if + 668);

    auto tr_z_xxyyyy_zzz = pbuffer.data(idx_dip_if + 669);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyy_xxx, tr_z_xxyy_xxz, tr_z_xxyy_xzz, tr_z_xxyyy_xxx, tr_z_xxyyy_xxz, tr_z_xxyyy_xzz, tr_z_xxyyyy_xxx, tr_z_xxyyyy_xxy, tr_z_xxyyyy_xxz, tr_z_xxyyyy_xyy, tr_z_xxyyyy_xyz, tr_z_xxyyyy_xzz, tr_z_xxyyyy_yyy, tr_z_xxyyyy_yyz, tr_z_xxyyyy_yzz, tr_z_xxyyyy_zzz, tr_z_xyyyy_xxy, tr_z_xyyyy_xy, tr_z_xyyyy_xyy, tr_z_xyyyy_xyz, tr_z_xyyyy_yy, tr_z_xyyyy_yyy, tr_z_xyyyy_yyz, tr_z_xyyyy_yz, tr_z_xyyyy_yzz, tr_z_xyyyy_zzz, tr_z_yyyy_xxy, tr_z_yyyy_xyy, tr_z_yyyy_xyz, tr_z_yyyy_yyy, tr_z_yyyy_yyz, tr_z_yyyy_yzz, tr_z_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyy_xxx[i] = 3.0 * tr_z_xxyy_xxx[i] * fe_0 + tr_z_xxyyy_xxx[i] * pa_y[i];

        tr_z_xxyyyy_xxy[i] = tr_z_yyyy_xxy[i] * fe_0 + 2.0 * tr_z_xyyyy_xy[i] * fe_0 + tr_z_xyyyy_xxy[i] * pa_x[i];

        tr_z_xxyyyy_xxz[i] = 3.0 * tr_z_xxyy_xxz[i] * fe_0 + tr_z_xxyyy_xxz[i] * pa_y[i];

        tr_z_xxyyyy_xyy[i] = tr_z_yyyy_xyy[i] * fe_0 + tr_z_xyyyy_yy[i] * fe_0 + tr_z_xyyyy_xyy[i] * pa_x[i];

        tr_z_xxyyyy_xyz[i] = tr_z_yyyy_xyz[i] * fe_0 + tr_z_xyyyy_yz[i] * fe_0 + tr_z_xyyyy_xyz[i] * pa_x[i];

        tr_z_xxyyyy_xzz[i] = 3.0 * tr_z_xxyy_xzz[i] * fe_0 + tr_z_xxyyy_xzz[i] * pa_y[i];

        tr_z_xxyyyy_yyy[i] = tr_z_yyyy_yyy[i] * fe_0 + tr_z_xyyyy_yyy[i] * pa_x[i];

        tr_z_xxyyyy_yyz[i] = tr_z_yyyy_yyz[i] * fe_0 + tr_z_xyyyy_yyz[i] * pa_x[i];

        tr_z_xxyyyy_yzz[i] = tr_z_yyyy_yzz[i] * fe_0 + tr_z_xyyyy_yzz[i] * pa_x[i];

        tr_z_xxyyyy_zzz[i] = tr_z_yyyy_zzz[i] * fe_0 + tr_z_xyyyy_zzz[i] * pa_x[i];
    }

    // Set up 670-680 components of targeted buffer : IF

    auto tr_z_xxyyyz_xxx = pbuffer.data(idx_dip_if + 670);

    auto tr_z_xxyyyz_xxy = pbuffer.data(idx_dip_if + 671);

    auto tr_z_xxyyyz_xxz = pbuffer.data(idx_dip_if + 672);

    auto tr_z_xxyyyz_xyy = pbuffer.data(idx_dip_if + 673);

    auto tr_z_xxyyyz_xyz = pbuffer.data(idx_dip_if + 674);

    auto tr_z_xxyyyz_xzz = pbuffer.data(idx_dip_if + 675);

    auto tr_z_xxyyyz_yyy = pbuffer.data(idx_dip_if + 676);

    auto tr_z_xxyyyz_yyz = pbuffer.data(idx_dip_if + 677);

    auto tr_z_xxyyyz_yzz = pbuffer.data(idx_dip_if + 678);

    auto tr_z_xxyyyz_zzz = pbuffer.data(idx_dip_if + 679);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxyyy_xxy, tr_z_xxyyy_xyy, tr_z_xxyyyz_xxx, tr_z_xxyyyz_xxy, tr_z_xxyyyz_xxz, tr_z_xxyyyz_xyy, tr_z_xxyyyz_xyz, tr_z_xxyyyz_xzz, tr_z_xxyyyz_yyy, tr_z_xxyyyz_yyz, tr_z_xxyyyz_yzz, tr_z_xxyyyz_zzz, tr_z_xxyyz_xxx, tr_z_xxyyz_xxz, tr_z_xxyyz_xzz, tr_z_xxyz_xxx, tr_z_xxyz_xxz, tr_z_xxyz_xzz, tr_z_xyyyz_xyz, tr_z_xyyyz_yyy, tr_z_xyyyz_yyz, tr_z_xyyyz_yz, tr_z_xyyyz_yzz, tr_z_xyyyz_zzz, tr_z_yyyz_xyz, tr_z_yyyz_yyy, tr_z_yyyz_yyz, tr_z_yyyz_yzz, tr_z_yyyz_zzz, ts_xxyyy_xxy, ts_xxyyy_xyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyz_xxx[i] = 2.0 * tr_z_xxyz_xxx[i] * fe_0 + tr_z_xxyyz_xxx[i] * pa_y[i];

        tr_z_xxyyyz_xxy[i] = ts_xxyyy_xxy[i] * fe_0 + tr_z_xxyyy_xxy[i] * pa_z[i];

        tr_z_xxyyyz_xxz[i] = 2.0 * tr_z_xxyz_xxz[i] * fe_0 + tr_z_xxyyz_xxz[i] * pa_y[i];

        tr_z_xxyyyz_xyy[i] = ts_xxyyy_xyy[i] * fe_0 + tr_z_xxyyy_xyy[i] * pa_z[i];

        tr_z_xxyyyz_xyz[i] = tr_z_yyyz_xyz[i] * fe_0 + tr_z_xyyyz_yz[i] * fe_0 + tr_z_xyyyz_xyz[i] * pa_x[i];

        tr_z_xxyyyz_xzz[i] = 2.0 * tr_z_xxyz_xzz[i] * fe_0 + tr_z_xxyyz_xzz[i] * pa_y[i];

        tr_z_xxyyyz_yyy[i] = tr_z_yyyz_yyy[i] * fe_0 + tr_z_xyyyz_yyy[i] * pa_x[i];

        tr_z_xxyyyz_yyz[i] = tr_z_yyyz_yyz[i] * fe_0 + tr_z_xyyyz_yyz[i] * pa_x[i];

        tr_z_xxyyyz_yzz[i] = tr_z_yyyz_yzz[i] * fe_0 + tr_z_xyyyz_yzz[i] * pa_x[i];

        tr_z_xxyyyz_zzz[i] = tr_z_yyyz_zzz[i] * fe_0 + tr_z_xyyyz_zzz[i] * pa_x[i];
    }

    // Set up 680-690 components of targeted buffer : IF

    auto tr_z_xxyyzz_xxx = pbuffer.data(idx_dip_if + 680);

    auto tr_z_xxyyzz_xxy = pbuffer.data(idx_dip_if + 681);

    auto tr_z_xxyyzz_xxz = pbuffer.data(idx_dip_if + 682);

    auto tr_z_xxyyzz_xyy = pbuffer.data(idx_dip_if + 683);

    auto tr_z_xxyyzz_xyz = pbuffer.data(idx_dip_if + 684);

    auto tr_z_xxyyzz_xzz = pbuffer.data(idx_dip_if + 685);

    auto tr_z_xxyyzz_yyy = pbuffer.data(idx_dip_if + 686);

    auto tr_z_xxyyzz_yyz = pbuffer.data(idx_dip_if + 687);

    auto tr_z_xxyyzz_yzz = pbuffer.data(idx_dip_if + 688);

    auto tr_z_xxyyzz_zzz = pbuffer.data(idx_dip_if + 689);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyyzz_xxx, tr_z_xxyyzz_xxy, tr_z_xxyyzz_xxz, tr_z_xxyyzz_xyy, tr_z_xxyyzz_xyz, tr_z_xxyyzz_xzz, tr_z_xxyyzz_yyy, tr_z_xxyyzz_yyz, tr_z_xxyyzz_yzz, tr_z_xxyyzz_zzz, tr_z_xxyzz_xxx, tr_z_xxyzz_xxz, tr_z_xxyzz_xzz, tr_z_xxzz_xxx, tr_z_xxzz_xxz, tr_z_xxzz_xzz, tr_z_xyyzz_xxy, tr_z_xyyzz_xy, tr_z_xyyzz_xyy, tr_z_xyyzz_xyz, tr_z_xyyzz_yy, tr_z_xyyzz_yyy, tr_z_xyyzz_yyz, tr_z_xyyzz_yz, tr_z_xyyzz_yzz, tr_z_xyyzz_zzz, tr_z_yyzz_xxy, tr_z_yyzz_xyy, tr_z_yyzz_xyz, tr_z_yyzz_yyy, tr_z_yyzz_yyz, tr_z_yyzz_yzz, tr_z_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyzz_xxx[i] = tr_z_xxzz_xxx[i] * fe_0 + tr_z_xxyzz_xxx[i] * pa_y[i];

        tr_z_xxyyzz_xxy[i] = tr_z_yyzz_xxy[i] * fe_0 + 2.0 * tr_z_xyyzz_xy[i] * fe_0 + tr_z_xyyzz_xxy[i] * pa_x[i];

        tr_z_xxyyzz_xxz[i] = tr_z_xxzz_xxz[i] * fe_0 + tr_z_xxyzz_xxz[i] * pa_y[i];

        tr_z_xxyyzz_xyy[i] = tr_z_yyzz_xyy[i] * fe_0 + tr_z_xyyzz_yy[i] * fe_0 + tr_z_xyyzz_xyy[i] * pa_x[i];

        tr_z_xxyyzz_xyz[i] = tr_z_yyzz_xyz[i] * fe_0 + tr_z_xyyzz_yz[i] * fe_0 + tr_z_xyyzz_xyz[i] * pa_x[i];

        tr_z_xxyyzz_xzz[i] = tr_z_xxzz_xzz[i] * fe_0 + tr_z_xxyzz_xzz[i] * pa_y[i];

        tr_z_xxyyzz_yyy[i] = tr_z_yyzz_yyy[i] * fe_0 + tr_z_xyyzz_yyy[i] * pa_x[i];

        tr_z_xxyyzz_yyz[i] = tr_z_yyzz_yyz[i] * fe_0 + tr_z_xyyzz_yyz[i] * pa_x[i];

        tr_z_xxyyzz_yzz[i] = tr_z_yyzz_yzz[i] * fe_0 + tr_z_xyyzz_yzz[i] * pa_x[i];

        tr_z_xxyyzz_zzz[i] = tr_z_yyzz_zzz[i] * fe_0 + tr_z_xyyzz_zzz[i] * pa_x[i];
    }

    // Set up 690-700 components of targeted buffer : IF

    auto tr_z_xxyzzz_xxx = pbuffer.data(idx_dip_if + 690);

    auto tr_z_xxyzzz_xxy = pbuffer.data(idx_dip_if + 691);

    auto tr_z_xxyzzz_xxz = pbuffer.data(idx_dip_if + 692);

    auto tr_z_xxyzzz_xyy = pbuffer.data(idx_dip_if + 693);

    auto tr_z_xxyzzz_xyz = pbuffer.data(idx_dip_if + 694);

    auto tr_z_xxyzzz_xzz = pbuffer.data(idx_dip_if + 695);

    auto tr_z_xxyzzz_yyy = pbuffer.data(idx_dip_if + 696);

    auto tr_z_xxyzzz_yyz = pbuffer.data(idx_dip_if + 697);

    auto tr_z_xxyzzz_yzz = pbuffer.data(idx_dip_if + 698);

    auto tr_z_xxyzzz_zzz = pbuffer.data(idx_dip_if + 699);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyzzz_xxx, tr_z_xxyzzz_xxy, tr_z_xxyzzz_xxz, tr_z_xxyzzz_xyy, tr_z_xxyzzz_xyz, tr_z_xxyzzz_xzz, tr_z_xxyzzz_yyy, tr_z_xxyzzz_yyz, tr_z_xxyzzz_yzz, tr_z_xxyzzz_zzz, tr_z_xxzzz_xx, tr_z_xxzzz_xxx, tr_z_xxzzz_xxy, tr_z_xxzzz_xxz, tr_z_xxzzz_xy, tr_z_xxzzz_xyy, tr_z_xxzzz_xyz, tr_z_xxzzz_xz, tr_z_xxzzz_xzz, tr_z_xxzzz_zzz, tr_z_xyzzz_yyy, tr_z_xyzzz_yyz, tr_z_xyzzz_yzz, tr_z_yzzz_yyy, tr_z_yzzz_yyz, tr_z_yzzz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzzz_xxx[i] = tr_z_xxzzz_xxx[i] * pa_y[i];

        tr_z_xxyzzz_xxy[i] = tr_z_xxzzz_xx[i] * fe_0 + tr_z_xxzzz_xxy[i] * pa_y[i];

        tr_z_xxyzzz_xxz[i] = tr_z_xxzzz_xxz[i] * pa_y[i];

        tr_z_xxyzzz_xyy[i] = 2.0 * tr_z_xxzzz_xy[i] * fe_0 + tr_z_xxzzz_xyy[i] * pa_y[i];

        tr_z_xxyzzz_xyz[i] = tr_z_xxzzz_xz[i] * fe_0 + tr_z_xxzzz_xyz[i] * pa_y[i];

        tr_z_xxyzzz_xzz[i] = tr_z_xxzzz_xzz[i] * pa_y[i];

        tr_z_xxyzzz_yyy[i] = tr_z_yzzz_yyy[i] * fe_0 + tr_z_xyzzz_yyy[i] * pa_x[i];

        tr_z_xxyzzz_yyz[i] = tr_z_yzzz_yyz[i] * fe_0 + tr_z_xyzzz_yyz[i] * pa_x[i];

        tr_z_xxyzzz_yzz[i] = tr_z_yzzz_yzz[i] * fe_0 + tr_z_xyzzz_yzz[i] * pa_x[i];

        tr_z_xxyzzz_zzz[i] = tr_z_xxzzz_zzz[i] * pa_y[i];
    }

    // Set up 700-710 components of targeted buffer : IF

    auto tr_z_xxzzzz_xxx = pbuffer.data(idx_dip_if + 700);

    auto tr_z_xxzzzz_xxy = pbuffer.data(idx_dip_if + 701);

    auto tr_z_xxzzzz_xxz = pbuffer.data(idx_dip_if + 702);

    auto tr_z_xxzzzz_xyy = pbuffer.data(idx_dip_if + 703);

    auto tr_z_xxzzzz_xyz = pbuffer.data(idx_dip_if + 704);

    auto tr_z_xxzzzz_xzz = pbuffer.data(idx_dip_if + 705);

    auto tr_z_xxzzzz_yyy = pbuffer.data(idx_dip_if + 706);

    auto tr_z_xxzzzz_yyz = pbuffer.data(idx_dip_if + 707);

    auto tr_z_xxzzzz_yzz = pbuffer.data(idx_dip_if + 708);

    auto tr_z_xxzzzz_zzz = pbuffer.data(idx_dip_if + 709);

    #pragma omp simd aligned(pa_x, tr_z_xxzzzz_xxx, tr_z_xxzzzz_xxy, tr_z_xxzzzz_xxz, tr_z_xxzzzz_xyy, tr_z_xxzzzz_xyz, tr_z_xxzzzz_xzz, tr_z_xxzzzz_yyy, tr_z_xxzzzz_yyz, tr_z_xxzzzz_yzz, tr_z_xxzzzz_zzz, tr_z_xzzzz_xx, tr_z_xzzzz_xxx, tr_z_xzzzz_xxy, tr_z_xzzzz_xxz, tr_z_xzzzz_xy, tr_z_xzzzz_xyy, tr_z_xzzzz_xyz, tr_z_xzzzz_xz, tr_z_xzzzz_xzz, tr_z_xzzzz_yy, tr_z_xzzzz_yyy, tr_z_xzzzz_yyz, tr_z_xzzzz_yz, tr_z_xzzzz_yzz, tr_z_xzzzz_zz, tr_z_xzzzz_zzz, tr_z_zzzz_xxx, tr_z_zzzz_xxy, tr_z_zzzz_xxz, tr_z_zzzz_xyy, tr_z_zzzz_xyz, tr_z_zzzz_xzz, tr_z_zzzz_yyy, tr_z_zzzz_yyz, tr_z_zzzz_yzz, tr_z_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzzz_xxx[i] = tr_z_zzzz_xxx[i] * fe_0 + 3.0 * tr_z_xzzzz_xx[i] * fe_0 + tr_z_xzzzz_xxx[i] * pa_x[i];

        tr_z_xxzzzz_xxy[i] = tr_z_zzzz_xxy[i] * fe_0 + 2.0 * tr_z_xzzzz_xy[i] * fe_0 + tr_z_xzzzz_xxy[i] * pa_x[i];

        tr_z_xxzzzz_xxz[i] = tr_z_zzzz_xxz[i] * fe_0 + 2.0 * tr_z_xzzzz_xz[i] * fe_0 + tr_z_xzzzz_xxz[i] * pa_x[i];

        tr_z_xxzzzz_xyy[i] = tr_z_zzzz_xyy[i] * fe_0 + tr_z_xzzzz_yy[i] * fe_0 + tr_z_xzzzz_xyy[i] * pa_x[i];

        tr_z_xxzzzz_xyz[i] = tr_z_zzzz_xyz[i] * fe_0 + tr_z_xzzzz_yz[i] * fe_0 + tr_z_xzzzz_xyz[i] * pa_x[i];

        tr_z_xxzzzz_xzz[i] = tr_z_zzzz_xzz[i] * fe_0 + tr_z_xzzzz_zz[i] * fe_0 + tr_z_xzzzz_xzz[i] * pa_x[i];

        tr_z_xxzzzz_yyy[i] = tr_z_zzzz_yyy[i] * fe_0 + tr_z_xzzzz_yyy[i] * pa_x[i];

        tr_z_xxzzzz_yyz[i] = tr_z_zzzz_yyz[i] * fe_0 + tr_z_xzzzz_yyz[i] * pa_x[i];

        tr_z_xxzzzz_yzz[i] = tr_z_zzzz_yzz[i] * fe_0 + tr_z_xzzzz_yzz[i] * pa_x[i];

        tr_z_xxzzzz_zzz[i] = tr_z_zzzz_zzz[i] * fe_0 + tr_z_xzzzz_zzz[i] * pa_x[i];
    }

    // Set up 710-720 components of targeted buffer : IF

    auto tr_z_xyyyyy_xxx = pbuffer.data(idx_dip_if + 710);

    auto tr_z_xyyyyy_xxy = pbuffer.data(idx_dip_if + 711);

    auto tr_z_xyyyyy_xxz = pbuffer.data(idx_dip_if + 712);

    auto tr_z_xyyyyy_xyy = pbuffer.data(idx_dip_if + 713);

    auto tr_z_xyyyyy_xyz = pbuffer.data(idx_dip_if + 714);

    auto tr_z_xyyyyy_xzz = pbuffer.data(idx_dip_if + 715);

    auto tr_z_xyyyyy_yyy = pbuffer.data(idx_dip_if + 716);

    auto tr_z_xyyyyy_yyz = pbuffer.data(idx_dip_if + 717);

    auto tr_z_xyyyyy_yzz = pbuffer.data(idx_dip_if + 718);

    auto tr_z_xyyyyy_zzz = pbuffer.data(idx_dip_if + 719);

    #pragma omp simd aligned(pa_x, tr_z_xyyyyy_xxx, tr_z_xyyyyy_xxy, tr_z_xyyyyy_xxz, tr_z_xyyyyy_xyy, tr_z_xyyyyy_xyz, tr_z_xyyyyy_xzz, tr_z_xyyyyy_yyy, tr_z_xyyyyy_yyz, tr_z_xyyyyy_yzz, tr_z_xyyyyy_zzz, tr_z_yyyyy_xx, tr_z_yyyyy_xxx, tr_z_yyyyy_xxy, tr_z_yyyyy_xxz, tr_z_yyyyy_xy, tr_z_yyyyy_xyy, tr_z_yyyyy_xyz, tr_z_yyyyy_xz, tr_z_yyyyy_xzz, tr_z_yyyyy_yy, tr_z_yyyyy_yyy, tr_z_yyyyy_yyz, tr_z_yyyyy_yz, tr_z_yyyyy_yzz, tr_z_yyyyy_zz, tr_z_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyy_xxx[i] = 3.0 * tr_z_yyyyy_xx[i] * fe_0 + tr_z_yyyyy_xxx[i] * pa_x[i];

        tr_z_xyyyyy_xxy[i] = 2.0 * tr_z_yyyyy_xy[i] * fe_0 + tr_z_yyyyy_xxy[i] * pa_x[i];

        tr_z_xyyyyy_xxz[i] = 2.0 * tr_z_yyyyy_xz[i] * fe_0 + tr_z_yyyyy_xxz[i] * pa_x[i];

        tr_z_xyyyyy_xyy[i] = tr_z_yyyyy_yy[i] * fe_0 + tr_z_yyyyy_xyy[i] * pa_x[i];

        tr_z_xyyyyy_xyz[i] = tr_z_yyyyy_yz[i] * fe_0 + tr_z_yyyyy_xyz[i] * pa_x[i];

        tr_z_xyyyyy_xzz[i] = tr_z_yyyyy_zz[i] * fe_0 + tr_z_yyyyy_xzz[i] * pa_x[i];

        tr_z_xyyyyy_yyy[i] = tr_z_yyyyy_yyy[i] * pa_x[i];

        tr_z_xyyyyy_yyz[i] = tr_z_yyyyy_yyz[i] * pa_x[i];

        tr_z_xyyyyy_yzz[i] = tr_z_yyyyy_yzz[i] * pa_x[i];

        tr_z_xyyyyy_zzz[i] = tr_z_yyyyy_zzz[i] * pa_x[i];
    }

    // Set up 720-730 components of targeted buffer : IF

    auto tr_z_xyyyyz_xxx = pbuffer.data(idx_dip_if + 720);

    auto tr_z_xyyyyz_xxy = pbuffer.data(idx_dip_if + 721);

    auto tr_z_xyyyyz_xxz = pbuffer.data(idx_dip_if + 722);

    auto tr_z_xyyyyz_xyy = pbuffer.data(idx_dip_if + 723);

    auto tr_z_xyyyyz_xyz = pbuffer.data(idx_dip_if + 724);

    auto tr_z_xyyyyz_xzz = pbuffer.data(idx_dip_if + 725);

    auto tr_z_xyyyyz_yyy = pbuffer.data(idx_dip_if + 726);

    auto tr_z_xyyyyz_yyz = pbuffer.data(idx_dip_if + 727);

    auto tr_z_xyyyyz_yzz = pbuffer.data(idx_dip_if + 728);

    auto tr_z_xyyyyz_zzz = pbuffer.data(idx_dip_if + 729);

    #pragma omp simd aligned(pa_x, tr_z_xyyyyz_xxx, tr_z_xyyyyz_xxy, tr_z_xyyyyz_xxz, tr_z_xyyyyz_xyy, tr_z_xyyyyz_xyz, tr_z_xyyyyz_xzz, tr_z_xyyyyz_yyy, tr_z_xyyyyz_yyz, tr_z_xyyyyz_yzz, tr_z_xyyyyz_zzz, tr_z_yyyyz_xx, tr_z_yyyyz_xxx, tr_z_yyyyz_xxy, tr_z_yyyyz_xxz, tr_z_yyyyz_xy, tr_z_yyyyz_xyy, tr_z_yyyyz_xyz, tr_z_yyyyz_xz, tr_z_yyyyz_xzz, tr_z_yyyyz_yy, tr_z_yyyyz_yyy, tr_z_yyyyz_yyz, tr_z_yyyyz_yz, tr_z_yyyyz_yzz, tr_z_yyyyz_zz, tr_z_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyz_xxx[i] = 3.0 * tr_z_yyyyz_xx[i] * fe_0 + tr_z_yyyyz_xxx[i] * pa_x[i];

        tr_z_xyyyyz_xxy[i] = 2.0 * tr_z_yyyyz_xy[i] * fe_0 + tr_z_yyyyz_xxy[i] * pa_x[i];

        tr_z_xyyyyz_xxz[i] = 2.0 * tr_z_yyyyz_xz[i] * fe_0 + tr_z_yyyyz_xxz[i] * pa_x[i];

        tr_z_xyyyyz_xyy[i] = tr_z_yyyyz_yy[i] * fe_0 + tr_z_yyyyz_xyy[i] * pa_x[i];

        tr_z_xyyyyz_xyz[i] = tr_z_yyyyz_yz[i] * fe_0 + tr_z_yyyyz_xyz[i] * pa_x[i];

        tr_z_xyyyyz_xzz[i] = tr_z_yyyyz_zz[i] * fe_0 + tr_z_yyyyz_xzz[i] * pa_x[i];

        tr_z_xyyyyz_yyy[i] = tr_z_yyyyz_yyy[i] * pa_x[i];

        tr_z_xyyyyz_yyz[i] = tr_z_yyyyz_yyz[i] * pa_x[i];

        tr_z_xyyyyz_yzz[i] = tr_z_yyyyz_yzz[i] * pa_x[i];

        tr_z_xyyyyz_zzz[i] = tr_z_yyyyz_zzz[i] * pa_x[i];
    }

    // Set up 730-740 components of targeted buffer : IF

    auto tr_z_xyyyzz_xxx = pbuffer.data(idx_dip_if + 730);

    auto tr_z_xyyyzz_xxy = pbuffer.data(idx_dip_if + 731);

    auto tr_z_xyyyzz_xxz = pbuffer.data(idx_dip_if + 732);

    auto tr_z_xyyyzz_xyy = pbuffer.data(idx_dip_if + 733);

    auto tr_z_xyyyzz_xyz = pbuffer.data(idx_dip_if + 734);

    auto tr_z_xyyyzz_xzz = pbuffer.data(idx_dip_if + 735);

    auto tr_z_xyyyzz_yyy = pbuffer.data(idx_dip_if + 736);

    auto tr_z_xyyyzz_yyz = pbuffer.data(idx_dip_if + 737);

    auto tr_z_xyyyzz_yzz = pbuffer.data(idx_dip_if + 738);

    auto tr_z_xyyyzz_zzz = pbuffer.data(idx_dip_if + 739);

    #pragma omp simd aligned(pa_x, tr_z_xyyyzz_xxx, tr_z_xyyyzz_xxy, tr_z_xyyyzz_xxz, tr_z_xyyyzz_xyy, tr_z_xyyyzz_xyz, tr_z_xyyyzz_xzz, tr_z_xyyyzz_yyy, tr_z_xyyyzz_yyz, tr_z_xyyyzz_yzz, tr_z_xyyyzz_zzz, tr_z_yyyzz_xx, tr_z_yyyzz_xxx, tr_z_yyyzz_xxy, tr_z_yyyzz_xxz, tr_z_yyyzz_xy, tr_z_yyyzz_xyy, tr_z_yyyzz_xyz, tr_z_yyyzz_xz, tr_z_yyyzz_xzz, tr_z_yyyzz_yy, tr_z_yyyzz_yyy, tr_z_yyyzz_yyz, tr_z_yyyzz_yz, tr_z_yyyzz_yzz, tr_z_yyyzz_zz, tr_z_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyzz_xxx[i] = 3.0 * tr_z_yyyzz_xx[i] * fe_0 + tr_z_yyyzz_xxx[i] * pa_x[i];

        tr_z_xyyyzz_xxy[i] = 2.0 * tr_z_yyyzz_xy[i] * fe_0 + tr_z_yyyzz_xxy[i] * pa_x[i];

        tr_z_xyyyzz_xxz[i] = 2.0 * tr_z_yyyzz_xz[i] * fe_0 + tr_z_yyyzz_xxz[i] * pa_x[i];

        tr_z_xyyyzz_xyy[i] = tr_z_yyyzz_yy[i] * fe_0 + tr_z_yyyzz_xyy[i] * pa_x[i];

        tr_z_xyyyzz_xyz[i] = tr_z_yyyzz_yz[i] * fe_0 + tr_z_yyyzz_xyz[i] * pa_x[i];

        tr_z_xyyyzz_xzz[i] = tr_z_yyyzz_zz[i] * fe_0 + tr_z_yyyzz_xzz[i] * pa_x[i];

        tr_z_xyyyzz_yyy[i] = tr_z_yyyzz_yyy[i] * pa_x[i];

        tr_z_xyyyzz_yyz[i] = tr_z_yyyzz_yyz[i] * pa_x[i];

        tr_z_xyyyzz_yzz[i] = tr_z_yyyzz_yzz[i] * pa_x[i];

        tr_z_xyyyzz_zzz[i] = tr_z_yyyzz_zzz[i] * pa_x[i];
    }

    // Set up 740-750 components of targeted buffer : IF

    auto tr_z_xyyzzz_xxx = pbuffer.data(idx_dip_if + 740);

    auto tr_z_xyyzzz_xxy = pbuffer.data(idx_dip_if + 741);

    auto tr_z_xyyzzz_xxz = pbuffer.data(idx_dip_if + 742);

    auto tr_z_xyyzzz_xyy = pbuffer.data(idx_dip_if + 743);

    auto tr_z_xyyzzz_xyz = pbuffer.data(idx_dip_if + 744);

    auto tr_z_xyyzzz_xzz = pbuffer.data(idx_dip_if + 745);

    auto tr_z_xyyzzz_yyy = pbuffer.data(idx_dip_if + 746);

    auto tr_z_xyyzzz_yyz = pbuffer.data(idx_dip_if + 747);

    auto tr_z_xyyzzz_yzz = pbuffer.data(idx_dip_if + 748);

    auto tr_z_xyyzzz_zzz = pbuffer.data(idx_dip_if + 749);

    #pragma omp simd aligned(pa_x, tr_z_xyyzzz_xxx, tr_z_xyyzzz_xxy, tr_z_xyyzzz_xxz, tr_z_xyyzzz_xyy, tr_z_xyyzzz_xyz, tr_z_xyyzzz_xzz, tr_z_xyyzzz_yyy, tr_z_xyyzzz_yyz, tr_z_xyyzzz_yzz, tr_z_xyyzzz_zzz, tr_z_yyzzz_xx, tr_z_yyzzz_xxx, tr_z_yyzzz_xxy, tr_z_yyzzz_xxz, tr_z_yyzzz_xy, tr_z_yyzzz_xyy, tr_z_yyzzz_xyz, tr_z_yyzzz_xz, tr_z_yyzzz_xzz, tr_z_yyzzz_yy, tr_z_yyzzz_yyy, tr_z_yyzzz_yyz, tr_z_yyzzz_yz, tr_z_yyzzz_yzz, tr_z_yyzzz_zz, tr_z_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzzz_xxx[i] = 3.0 * tr_z_yyzzz_xx[i] * fe_0 + tr_z_yyzzz_xxx[i] * pa_x[i];

        tr_z_xyyzzz_xxy[i] = 2.0 * tr_z_yyzzz_xy[i] * fe_0 + tr_z_yyzzz_xxy[i] * pa_x[i];

        tr_z_xyyzzz_xxz[i] = 2.0 * tr_z_yyzzz_xz[i] * fe_0 + tr_z_yyzzz_xxz[i] * pa_x[i];

        tr_z_xyyzzz_xyy[i] = tr_z_yyzzz_yy[i] * fe_0 + tr_z_yyzzz_xyy[i] * pa_x[i];

        tr_z_xyyzzz_xyz[i] = tr_z_yyzzz_yz[i] * fe_0 + tr_z_yyzzz_xyz[i] * pa_x[i];

        tr_z_xyyzzz_xzz[i] = tr_z_yyzzz_zz[i] * fe_0 + tr_z_yyzzz_xzz[i] * pa_x[i];

        tr_z_xyyzzz_yyy[i] = tr_z_yyzzz_yyy[i] * pa_x[i];

        tr_z_xyyzzz_yyz[i] = tr_z_yyzzz_yyz[i] * pa_x[i];

        tr_z_xyyzzz_yzz[i] = tr_z_yyzzz_yzz[i] * pa_x[i];

        tr_z_xyyzzz_zzz[i] = tr_z_yyzzz_zzz[i] * pa_x[i];
    }

    // Set up 750-760 components of targeted buffer : IF

    auto tr_z_xyzzzz_xxx = pbuffer.data(idx_dip_if + 750);

    auto tr_z_xyzzzz_xxy = pbuffer.data(idx_dip_if + 751);

    auto tr_z_xyzzzz_xxz = pbuffer.data(idx_dip_if + 752);

    auto tr_z_xyzzzz_xyy = pbuffer.data(idx_dip_if + 753);

    auto tr_z_xyzzzz_xyz = pbuffer.data(idx_dip_if + 754);

    auto tr_z_xyzzzz_xzz = pbuffer.data(idx_dip_if + 755);

    auto tr_z_xyzzzz_yyy = pbuffer.data(idx_dip_if + 756);

    auto tr_z_xyzzzz_yyz = pbuffer.data(idx_dip_if + 757);

    auto tr_z_xyzzzz_yzz = pbuffer.data(idx_dip_if + 758);

    auto tr_z_xyzzzz_zzz = pbuffer.data(idx_dip_if + 759);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyzzzz_xxx, tr_z_xyzzzz_xxy, tr_z_xyzzzz_xxz, tr_z_xyzzzz_xyy, tr_z_xyzzzz_xyz, tr_z_xyzzzz_xzz, tr_z_xyzzzz_yyy, tr_z_xyzzzz_yyz, tr_z_xyzzzz_yzz, tr_z_xyzzzz_zzz, tr_z_xzzzz_xxx, tr_z_xzzzz_xxz, tr_z_xzzzz_xzz, tr_z_yzzzz_xxy, tr_z_yzzzz_xy, tr_z_yzzzz_xyy, tr_z_yzzzz_xyz, tr_z_yzzzz_yy, tr_z_yzzzz_yyy, tr_z_yzzzz_yyz, tr_z_yzzzz_yz, tr_z_yzzzz_yzz, tr_z_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzzz_xxx[i] = tr_z_xzzzz_xxx[i] * pa_y[i];

        tr_z_xyzzzz_xxy[i] = 2.0 * tr_z_yzzzz_xy[i] * fe_0 + tr_z_yzzzz_xxy[i] * pa_x[i];

        tr_z_xyzzzz_xxz[i] = tr_z_xzzzz_xxz[i] * pa_y[i];

        tr_z_xyzzzz_xyy[i] = tr_z_yzzzz_yy[i] * fe_0 + tr_z_yzzzz_xyy[i] * pa_x[i];

        tr_z_xyzzzz_xyz[i] = tr_z_yzzzz_yz[i] * fe_0 + tr_z_yzzzz_xyz[i] * pa_x[i];

        tr_z_xyzzzz_xzz[i] = tr_z_xzzzz_xzz[i] * pa_y[i];

        tr_z_xyzzzz_yyy[i] = tr_z_yzzzz_yyy[i] * pa_x[i];

        tr_z_xyzzzz_yyz[i] = tr_z_yzzzz_yyz[i] * pa_x[i];

        tr_z_xyzzzz_yzz[i] = tr_z_yzzzz_yzz[i] * pa_x[i];

        tr_z_xyzzzz_zzz[i] = tr_z_yzzzz_zzz[i] * pa_x[i];
    }

    // Set up 760-770 components of targeted buffer : IF

    auto tr_z_xzzzzz_xxx = pbuffer.data(idx_dip_if + 760);

    auto tr_z_xzzzzz_xxy = pbuffer.data(idx_dip_if + 761);

    auto tr_z_xzzzzz_xxz = pbuffer.data(idx_dip_if + 762);

    auto tr_z_xzzzzz_xyy = pbuffer.data(idx_dip_if + 763);

    auto tr_z_xzzzzz_xyz = pbuffer.data(idx_dip_if + 764);

    auto tr_z_xzzzzz_xzz = pbuffer.data(idx_dip_if + 765);

    auto tr_z_xzzzzz_yyy = pbuffer.data(idx_dip_if + 766);

    auto tr_z_xzzzzz_yyz = pbuffer.data(idx_dip_if + 767);

    auto tr_z_xzzzzz_yzz = pbuffer.data(idx_dip_if + 768);

    auto tr_z_xzzzzz_zzz = pbuffer.data(idx_dip_if + 769);

    #pragma omp simd aligned(pa_x, tr_z_xzzzzz_xxx, tr_z_xzzzzz_xxy, tr_z_xzzzzz_xxz, tr_z_xzzzzz_xyy, tr_z_xzzzzz_xyz, tr_z_xzzzzz_xzz, tr_z_xzzzzz_yyy, tr_z_xzzzzz_yyz, tr_z_xzzzzz_yzz, tr_z_xzzzzz_zzz, tr_z_zzzzz_xx, tr_z_zzzzz_xxx, tr_z_zzzzz_xxy, tr_z_zzzzz_xxz, tr_z_zzzzz_xy, tr_z_zzzzz_xyy, tr_z_zzzzz_xyz, tr_z_zzzzz_xz, tr_z_zzzzz_xzz, tr_z_zzzzz_yy, tr_z_zzzzz_yyy, tr_z_zzzzz_yyz, tr_z_zzzzz_yz, tr_z_zzzzz_yzz, tr_z_zzzzz_zz, tr_z_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzzz_xxx[i] = 3.0 * tr_z_zzzzz_xx[i] * fe_0 + tr_z_zzzzz_xxx[i] * pa_x[i];

        tr_z_xzzzzz_xxy[i] = 2.0 * tr_z_zzzzz_xy[i] * fe_0 + tr_z_zzzzz_xxy[i] * pa_x[i];

        tr_z_xzzzzz_xxz[i] = 2.0 * tr_z_zzzzz_xz[i] * fe_0 + tr_z_zzzzz_xxz[i] * pa_x[i];

        tr_z_xzzzzz_xyy[i] = tr_z_zzzzz_yy[i] * fe_0 + tr_z_zzzzz_xyy[i] * pa_x[i];

        tr_z_xzzzzz_xyz[i] = tr_z_zzzzz_yz[i] * fe_0 + tr_z_zzzzz_xyz[i] * pa_x[i];

        tr_z_xzzzzz_xzz[i] = tr_z_zzzzz_zz[i] * fe_0 + tr_z_zzzzz_xzz[i] * pa_x[i];

        tr_z_xzzzzz_yyy[i] = tr_z_zzzzz_yyy[i] * pa_x[i];

        tr_z_xzzzzz_yyz[i] = tr_z_zzzzz_yyz[i] * pa_x[i];

        tr_z_xzzzzz_yzz[i] = tr_z_zzzzz_yzz[i] * pa_x[i];

        tr_z_xzzzzz_zzz[i] = tr_z_zzzzz_zzz[i] * pa_x[i];
    }

    // Set up 770-780 components of targeted buffer : IF

    auto tr_z_yyyyyy_xxx = pbuffer.data(idx_dip_if + 770);

    auto tr_z_yyyyyy_xxy = pbuffer.data(idx_dip_if + 771);

    auto tr_z_yyyyyy_xxz = pbuffer.data(idx_dip_if + 772);

    auto tr_z_yyyyyy_xyy = pbuffer.data(idx_dip_if + 773);

    auto tr_z_yyyyyy_xyz = pbuffer.data(idx_dip_if + 774);

    auto tr_z_yyyyyy_xzz = pbuffer.data(idx_dip_if + 775);

    auto tr_z_yyyyyy_yyy = pbuffer.data(idx_dip_if + 776);

    auto tr_z_yyyyyy_yyz = pbuffer.data(idx_dip_if + 777);

    auto tr_z_yyyyyy_yzz = pbuffer.data(idx_dip_if + 778);

    auto tr_z_yyyyyy_zzz = pbuffer.data(idx_dip_if + 779);

    #pragma omp simd aligned(pa_y, tr_z_yyyy_xxx, tr_z_yyyy_xxy, tr_z_yyyy_xxz, tr_z_yyyy_xyy, tr_z_yyyy_xyz, tr_z_yyyy_xzz, tr_z_yyyy_yyy, tr_z_yyyy_yyz, tr_z_yyyy_yzz, tr_z_yyyy_zzz, tr_z_yyyyy_xx, tr_z_yyyyy_xxx, tr_z_yyyyy_xxy, tr_z_yyyyy_xxz, tr_z_yyyyy_xy, tr_z_yyyyy_xyy, tr_z_yyyyy_xyz, tr_z_yyyyy_xz, tr_z_yyyyy_xzz, tr_z_yyyyy_yy, tr_z_yyyyy_yyy, tr_z_yyyyy_yyz, tr_z_yyyyy_yz, tr_z_yyyyy_yzz, tr_z_yyyyy_zz, tr_z_yyyyy_zzz, tr_z_yyyyyy_xxx, tr_z_yyyyyy_xxy, tr_z_yyyyyy_xxz, tr_z_yyyyyy_xyy, tr_z_yyyyyy_xyz, tr_z_yyyyyy_xzz, tr_z_yyyyyy_yyy, tr_z_yyyyyy_yyz, tr_z_yyyyyy_yzz, tr_z_yyyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyy_xxx[i] = 5.0 * tr_z_yyyy_xxx[i] * fe_0 + tr_z_yyyyy_xxx[i] * pa_y[i];

        tr_z_yyyyyy_xxy[i] = 5.0 * tr_z_yyyy_xxy[i] * fe_0 + tr_z_yyyyy_xx[i] * fe_0 + tr_z_yyyyy_xxy[i] * pa_y[i];

        tr_z_yyyyyy_xxz[i] = 5.0 * tr_z_yyyy_xxz[i] * fe_0 + tr_z_yyyyy_xxz[i] * pa_y[i];

        tr_z_yyyyyy_xyy[i] = 5.0 * tr_z_yyyy_xyy[i] * fe_0 + 2.0 * tr_z_yyyyy_xy[i] * fe_0 + tr_z_yyyyy_xyy[i] * pa_y[i];

        tr_z_yyyyyy_xyz[i] = 5.0 * tr_z_yyyy_xyz[i] * fe_0 + tr_z_yyyyy_xz[i] * fe_0 + tr_z_yyyyy_xyz[i] * pa_y[i];

        tr_z_yyyyyy_xzz[i] = 5.0 * tr_z_yyyy_xzz[i] * fe_0 + tr_z_yyyyy_xzz[i] * pa_y[i];

        tr_z_yyyyyy_yyy[i] = 5.0 * tr_z_yyyy_yyy[i] * fe_0 + 3.0 * tr_z_yyyyy_yy[i] * fe_0 + tr_z_yyyyy_yyy[i] * pa_y[i];

        tr_z_yyyyyy_yyz[i] = 5.0 * tr_z_yyyy_yyz[i] * fe_0 + 2.0 * tr_z_yyyyy_yz[i] * fe_0 + tr_z_yyyyy_yyz[i] * pa_y[i];

        tr_z_yyyyyy_yzz[i] = 5.0 * tr_z_yyyy_yzz[i] * fe_0 + tr_z_yyyyy_zz[i] * fe_0 + tr_z_yyyyy_yzz[i] * pa_y[i];

        tr_z_yyyyyy_zzz[i] = 5.0 * tr_z_yyyy_zzz[i] * fe_0 + tr_z_yyyyy_zzz[i] * pa_y[i];
    }

    // Set up 780-790 components of targeted buffer : IF

    auto tr_z_yyyyyz_xxx = pbuffer.data(idx_dip_if + 780);

    auto tr_z_yyyyyz_xxy = pbuffer.data(idx_dip_if + 781);

    auto tr_z_yyyyyz_xxz = pbuffer.data(idx_dip_if + 782);

    auto tr_z_yyyyyz_xyy = pbuffer.data(idx_dip_if + 783);

    auto tr_z_yyyyyz_xyz = pbuffer.data(idx_dip_if + 784);

    auto tr_z_yyyyyz_xzz = pbuffer.data(idx_dip_if + 785);

    auto tr_z_yyyyyz_yyy = pbuffer.data(idx_dip_if + 786);

    auto tr_z_yyyyyz_yyz = pbuffer.data(idx_dip_if + 787);

    auto tr_z_yyyyyz_yzz = pbuffer.data(idx_dip_if + 788);

    auto tr_z_yyyyyz_zzz = pbuffer.data(idx_dip_if + 789);

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yyyyy_xxy, tr_z_yyyyy_xyy, tr_z_yyyyy_yyy, tr_z_yyyyyz_xxx, tr_z_yyyyyz_xxy, tr_z_yyyyyz_xxz, tr_z_yyyyyz_xyy, tr_z_yyyyyz_xyz, tr_z_yyyyyz_xzz, tr_z_yyyyyz_yyy, tr_z_yyyyyz_yyz, tr_z_yyyyyz_yzz, tr_z_yyyyyz_zzz, tr_z_yyyyz_xxx, tr_z_yyyyz_xxz, tr_z_yyyyz_xyz, tr_z_yyyyz_xz, tr_z_yyyyz_xzz, tr_z_yyyyz_yyz, tr_z_yyyyz_yz, tr_z_yyyyz_yzz, tr_z_yyyyz_zz, tr_z_yyyyz_zzz, tr_z_yyyz_xxx, tr_z_yyyz_xxz, tr_z_yyyz_xyz, tr_z_yyyz_xzz, tr_z_yyyz_yyz, tr_z_yyyz_yzz, tr_z_yyyz_zzz, ts_yyyyy_xxy, ts_yyyyy_xyy, ts_yyyyy_yyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyz_xxx[i] = 4.0 * tr_z_yyyz_xxx[i] * fe_0 + tr_z_yyyyz_xxx[i] * pa_y[i];

        tr_z_yyyyyz_xxy[i] = ts_yyyyy_xxy[i] * fe_0 + tr_z_yyyyy_xxy[i] * pa_z[i];

        tr_z_yyyyyz_xxz[i] = 4.0 * tr_z_yyyz_xxz[i] * fe_0 + tr_z_yyyyz_xxz[i] * pa_y[i];

        tr_z_yyyyyz_xyy[i] = ts_yyyyy_xyy[i] * fe_0 + tr_z_yyyyy_xyy[i] * pa_z[i];

        tr_z_yyyyyz_xyz[i] = 4.0 * tr_z_yyyz_xyz[i] * fe_0 + tr_z_yyyyz_xz[i] * fe_0 + tr_z_yyyyz_xyz[i] * pa_y[i];

        tr_z_yyyyyz_xzz[i] = 4.0 * tr_z_yyyz_xzz[i] * fe_0 + tr_z_yyyyz_xzz[i] * pa_y[i];

        tr_z_yyyyyz_yyy[i] = ts_yyyyy_yyy[i] * fe_0 + tr_z_yyyyy_yyy[i] * pa_z[i];

        tr_z_yyyyyz_yyz[i] = 4.0 * tr_z_yyyz_yyz[i] * fe_0 + 2.0 * tr_z_yyyyz_yz[i] * fe_0 + tr_z_yyyyz_yyz[i] * pa_y[i];

        tr_z_yyyyyz_yzz[i] = 4.0 * tr_z_yyyz_yzz[i] * fe_0 + tr_z_yyyyz_zz[i] * fe_0 + tr_z_yyyyz_yzz[i] * pa_y[i];

        tr_z_yyyyyz_zzz[i] = 4.0 * tr_z_yyyz_zzz[i] * fe_0 + tr_z_yyyyz_zzz[i] * pa_y[i];
    }

    // Set up 790-800 components of targeted buffer : IF

    auto tr_z_yyyyzz_xxx = pbuffer.data(idx_dip_if + 790);

    auto tr_z_yyyyzz_xxy = pbuffer.data(idx_dip_if + 791);

    auto tr_z_yyyyzz_xxz = pbuffer.data(idx_dip_if + 792);

    auto tr_z_yyyyzz_xyy = pbuffer.data(idx_dip_if + 793);

    auto tr_z_yyyyzz_xyz = pbuffer.data(idx_dip_if + 794);

    auto tr_z_yyyyzz_xzz = pbuffer.data(idx_dip_if + 795);

    auto tr_z_yyyyzz_yyy = pbuffer.data(idx_dip_if + 796);

    auto tr_z_yyyyzz_yyz = pbuffer.data(idx_dip_if + 797);

    auto tr_z_yyyyzz_yzz = pbuffer.data(idx_dip_if + 798);

    auto tr_z_yyyyzz_zzz = pbuffer.data(idx_dip_if + 799);

    #pragma omp simd aligned(pa_y, tr_z_yyyyzz_xxx, tr_z_yyyyzz_xxy, tr_z_yyyyzz_xxz, tr_z_yyyyzz_xyy, tr_z_yyyyzz_xyz, tr_z_yyyyzz_xzz, tr_z_yyyyzz_yyy, tr_z_yyyyzz_yyz, tr_z_yyyyzz_yzz, tr_z_yyyyzz_zzz, tr_z_yyyzz_xx, tr_z_yyyzz_xxx, tr_z_yyyzz_xxy, tr_z_yyyzz_xxz, tr_z_yyyzz_xy, tr_z_yyyzz_xyy, tr_z_yyyzz_xyz, tr_z_yyyzz_xz, tr_z_yyyzz_xzz, tr_z_yyyzz_yy, tr_z_yyyzz_yyy, tr_z_yyyzz_yyz, tr_z_yyyzz_yz, tr_z_yyyzz_yzz, tr_z_yyyzz_zz, tr_z_yyyzz_zzz, tr_z_yyzz_xxx, tr_z_yyzz_xxy, tr_z_yyzz_xxz, tr_z_yyzz_xyy, tr_z_yyzz_xyz, tr_z_yyzz_xzz, tr_z_yyzz_yyy, tr_z_yyzz_yyz, tr_z_yyzz_yzz, tr_z_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyzz_xxx[i] = 3.0 * tr_z_yyzz_xxx[i] * fe_0 + tr_z_yyyzz_xxx[i] * pa_y[i];

        tr_z_yyyyzz_xxy[i] = 3.0 * tr_z_yyzz_xxy[i] * fe_0 + tr_z_yyyzz_xx[i] * fe_0 + tr_z_yyyzz_xxy[i] * pa_y[i];

        tr_z_yyyyzz_xxz[i] = 3.0 * tr_z_yyzz_xxz[i] * fe_0 + tr_z_yyyzz_xxz[i] * pa_y[i];

        tr_z_yyyyzz_xyy[i] = 3.0 * tr_z_yyzz_xyy[i] * fe_0 + 2.0 * tr_z_yyyzz_xy[i] * fe_0 + tr_z_yyyzz_xyy[i] * pa_y[i];

        tr_z_yyyyzz_xyz[i] = 3.0 * tr_z_yyzz_xyz[i] * fe_0 + tr_z_yyyzz_xz[i] * fe_0 + tr_z_yyyzz_xyz[i] * pa_y[i];

        tr_z_yyyyzz_xzz[i] = 3.0 * tr_z_yyzz_xzz[i] * fe_0 + tr_z_yyyzz_xzz[i] * pa_y[i];

        tr_z_yyyyzz_yyy[i] = 3.0 * tr_z_yyzz_yyy[i] * fe_0 + 3.0 * tr_z_yyyzz_yy[i] * fe_0 + tr_z_yyyzz_yyy[i] * pa_y[i];

        tr_z_yyyyzz_yyz[i] = 3.0 * tr_z_yyzz_yyz[i] * fe_0 + 2.0 * tr_z_yyyzz_yz[i] * fe_0 + tr_z_yyyzz_yyz[i] * pa_y[i];

        tr_z_yyyyzz_yzz[i] = 3.0 * tr_z_yyzz_yzz[i] * fe_0 + tr_z_yyyzz_zz[i] * fe_0 + tr_z_yyyzz_yzz[i] * pa_y[i];

        tr_z_yyyyzz_zzz[i] = 3.0 * tr_z_yyzz_zzz[i] * fe_0 + tr_z_yyyzz_zzz[i] * pa_y[i];
    }

    // Set up 800-810 components of targeted buffer : IF

    auto tr_z_yyyzzz_xxx = pbuffer.data(idx_dip_if + 800);

    auto tr_z_yyyzzz_xxy = pbuffer.data(idx_dip_if + 801);

    auto tr_z_yyyzzz_xxz = pbuffer.data(idx_dip_if + 802);

    auto tr_z_yyyzzz_xyy = pbuffer.data(idx_dip_if + 803);

    auto tr_z_yyyzzz_xyz = pbuffer.data(idx_dip_if + 804);

    auto tr_z_yyyzzz_xzz = pbuffer.data(idx_dip_if + 805);

    auto tr_z_yyyzzz_yyy = pbuffer.data(idx_dip_if + 806);

    auto tr_z_yyyzzz_yyz = pbuffer.data(idx_dip_if + 807);

    auto tr_z_yyyzzz_yzz = pbuffer.data(idx_dip_if + 808);

    auto tr_z_yyyzzz_zzz = pbuffer.data(idx_dip_if + 809);

    #pragma omp simd aligned(pa_y, tr_z_yyyzzz_xxx, tr_z_yyyzzz_xxy, tr_z_yyyzzz_xxz, tr_z_yyyzzz_xyy, tr_z_yyyzzz_xyz, tr_z_yyyzzz_xzz, tr_z_yyyzzz_yyy, tr_z_yyyzzz_yyz, tr_z_yyyzzz_yzz, tr_z_yyyzzz_zzz, tr_z_yyzzz_xx, tr_z_yyzzz_xxx, tr_z_yyzzz_xxy, tr_z_yyzzz_xxz, tr_z_yyzzz_xy, tr_z_yyzzz_xyy, tr_z_yyzzz_xyz, tr_z_yyzzz_xz, tr_z_yyzzz_xzz, tr_z_yyzzz_yy, tr_z_yyzzz_yyy, tr_z_yyzzz_yyz, tr_z_yyzzz_yz, tr_z_yyzzz_yzz, tr_z_yyzzz_zz, tr_z_yyzzz_zzz, tr_z_yzzz_xxx, tr_z_yzzz_xxy, tr_z_yzzz_xxz, tr_z_yzzz_xyy, tr_z_yzzz_xyz, tr_z_yzzz_xzz, tr_z_yzzz_yyy, tr_z_yzzz_yyz, tr_z_yzzz_yzz, tr_z_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzzz_xxx[i] = 2.0 * tr_z_yzzz_xxx[i] * fe_0 + tr_z_yyzzz_xxx[i] * pa_y[i];

        tr_z_yyyzzz_xxy[i] = 2.0 * tr_z_yzzz_xxy[i] * fe_0 + tr_z_yyzzz_xx[i] * fe_0 + tr_z_yyzzz_xxy[i] * pa_y[i];

        tr_z_yyyzzz_xxz[i] = 2.0 * tr_z_yzzz_xxz[i] * fe_0 + tr_z_yyzzz_xxz[i] * pa_y[i];

        tr_z_yyyzzz_xyy[i] = 2.0 * tr_z_yzzz_xyy[i] * fe_0 + 2.0 * tr_z_yyzzz_xy[i] * fe_0 + tr_z_yyzzz_xyy[i] * pa_y[i];

        tr_z_yyyzzz_xyz[i] = 2.0 * tr_z_yzzz_xyz[i] * fe_0 + tr_z_yyzzz_xz[i] * fe_0 + tr_z_yyzzz_xyz[i] * pa_y[i];

        tr_z_yyyzzz_xzz[i] = 2.0 * tr_z_yzzz_xzz[i] * fe_0 + tr_z_yyzzz_xzz[i] * pa_y[i];

        tr_z_yyyzzz_yyy[i] = 2.0 * tr_z_yzzz_yyy[i] * fe_0 + 3.0 * tr_z_yyzzz_yy[i] * fe_0 + tr_z_yyzzz_yyy[i] * pa_y[i];

        tr_z_yyyzzz_yyz[i] = 2.0 * tr_z_yzzz_yyz[i] * fe_0 + 2.0 * tr_z_yyzzz_yz[i] * fe_0 + tr_z_yyzzz_yyz[i] * pa_y[i];

        tr_z_yyyzzz_yzz[i] = 2.0 * tr_z_yzzz_yzz[i] * fe_0 + tr_z_yyzzz_zz[i] * fe_0 + tr_z_yyzzz_yzz[i] * pa_y[i];

        tr_z_yyyzzz_zzz[i] = 2.0 * tr_z_yzzz_zzz[i] * fe_0 + tr_z_yyzzz_zzz[i] * pa_y[i];
    }

    // Set up 810-820 components of targeted buffer : IF

    auto tr_z_yyzzzz_xxx = pbuffer.data(idx_dip_if + 810);

    auto tr_z_yyzzzz_xxy = pbuffer.data(idx_dip_if + 811);

    auto tr_z_yyzzzz_xxz = pbuffer.data(idx_dip_if + 812);

    auto tr_z_yyzzzz_xyy = pbuffer.data(idx_dip_if + 813);

    auto tr_z_yyzzzz_xyz = pbuffer.data(idx_dip_if + 814);

    auto tr_z_yyzzzz_xzz = pbuffer.data(idx_dip_if + 815);

    auto tr_z_yyzzzz_yyy = pbuffer.data(idx_dip_if + 816);

    auto tr_z_yyzzzz_yyz = pbuffer.data(idx_dip_if + 817);

    auto tr_z_yyzzzz_yzz = pbuffer.data(idx_dip_if + 818);

    auto tr_z_yyzzzz_zzz = pbuffer.data(idx_dip_if + 819);

    #pragma omp simd aligned(pa_y, tr_z_yyzzzz_xxx, tr_z_yyzzzz_xxy, tr_z_yyzzzz_xxz, tr_z_yyzzzz_xyy, tr_z_yyzzzz_xyz, tr_z_yyzzzz_xzz, tr_z_yyzzzz_yyy, tr_z_yyzzzz_yyz, tr_z_yyzzzz_yzz, tr_z_yyzzzz_zzz, tr_z_yzzzz_xx, tr_z_yzzzz_xxx, tr_z_yzzzz_xxy, tr_z_yzzzz_xxz, tr_z_yzzzz_xy, tr_z_yzzzz_xyy, tr_z_yzzzz_xyz, tr_z_yzzzz_xz, tr_z_yzzzz_xzz, tr_z_yzzzz_yy, tr_z_yzzzz_yyy, tr_z_yzzzz_yyz, tr_z_yzzzz_yz, tr_z_yzzzz_yzz, tr_z_yzzzz_zz, tr_z_yzzzz_zzz, tr_z_zzzz_xxx, tr_z_zzzz_xxy, tr_z_zzzz_xxz, tr_z_zzzz_xyy, tr_z_zzzz_xyz, tr_z_zzzz_xzz, tr_z_zzzz_yyy, tr_z_zzzz_yyz, tr_z_zzzz_yzz, tr_z_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzzz_xxx[i] = tr_z_zzzz_xxx[i] * fe_0 + tr_z_yzzzz_xxx[i] * pa_y[i];

        tr_z_yyzzzz_xxy[i] = tr_z_zzzz_xxy[i] * fe_0 + tr_z_yzzzz_xx[i] * fe_0 + tr_z_yzzzz_xxy[i] * pa_y[i];

        tr_z_yyzzzz_xxz[i] = tr_z_zzzz_xxz[i] * fe_0 + tr_z_yzzzz_xxz[i] * pa_y[i];

        tr_z_yyzzzz_xyy[i] = tr_z_zzzz_xyy[i] * fe_0 + 2.0 * tr_z_yzzzz_xy[i] * fe_0 + tr_z_yzzzz_xyy[i] * pa_y[i];

        tr_z_yyzzzz_xyz[i] = tr_z_zzzz_xyz[i] * fe_0 + tr_z_yzzzz_xz[i] * fe_0 + tr_z_yzzzz_xyz[i] * pa_y[i];

        tr_z_yyzzzz_xzz[i] = tr_z_zzzz_xzz[i] * fe_0 + tr_z_yzzzz_xzz[i] * pa_y[i];

        tr_z_yyzzzz_yyy[i] = tr_z_zzzz_yyy[i] * fe_0 + 3.0 * tr_z_yzzzz_yy[i] * fe_0 + tr_z_yzzzz_yyy[i] * pa_y[i];

        tr_z_yyzzzz_yyz[i] = tr_z_zzzz_yyz[i] * fe_0 + 2.0 * tr_z_yzzzz_yz[i] * fe_0 + tr_z_yzzzz_yyz[i] * pa_y[i];

        tr_z_yyzzzz_yzz[i] = tr_z_zzzz_yzz[i] * fe_0 + tr_z_yzzzz_zz[i] * fe_0 + tr_z_yzzzz_yzz[i] * pa_y[i];

        tr_z_yyzzzz_zzz[i] = tr_z_zzzz_zzz[i] * fe_0 + tr_z_yzzzz_zzz[i] * pa_y[i];
    }

    // Set up 820-830 components of targeted buffer : IF

    auto tr_z_yzzzzz_xxx = pbuffer.data(idx_dip_if + 820);

    auto tr_z_yzzzzz_xxy = pbuffer.data(idx_dip_if + 821);

    auto tr_z_yzzzzz_xxz = pbuffer.data(idx_dip_if + 822);

    auto tr_z_yzzzzz_xyy = pbuffer.data(idx_dip_if + 823);

    auto tr_z_yzzzzz_xyz = pbuffer.data(idx_dip_if + 824);

    auto tr_z_yzzzzz_xzz = pbuffer.data(idx_dip_if + 825);

    auto tr_z_yzzzzz_yyy = pbuffer.data(idx_dip_if + 826);

    auto tr_z_yzzzzz_yyz = pbuffer.data(idx_dip_if + 827);

    auto tr_z_yzzzzz_yzz = pbuffer.data(idx_dip_if + 828);

    auto tr_z_yzzzzz_zzz = pbuffer.data(idx_dip_if + 829);

    #pragma omp simd aligned(pa_y, tr_z_yzzzzz_xxx, tr_z_yzzzzz_xxy, tr_z_yzzzzz_xxz, tr_z_yzzzzz_xyy, tr_z_yzzzzz_xyz, tr_z_yzzzzz_xzz, tr_z_yzzzzz_yyy, tr_z_yzzzzz_yyz, tr_z_yzzzzz_yzz, tr_z_yzzzzz_zzz, tr_z_zzzzz_xx, tr_z_zzzzz_xxx, tr_z_zzzzz_xxy, tr_z_zzzzz_xxz, tr_z_zzzzz_xy, tr_z_zzzzz_xyy, tr_z_zzzzz_xyz, tr_z_zzzzz_xz, tr_z_zzzzz_xzz, tr_z_zzzzz_yy, tr_z_zzzzz_yyy, tr_z_zzzzz_yyz, tr_z_zzzzz_yz, tr_z_zzzzz_yzz, tr_z_zzzzz_zz, tr_z_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzzz_xxx[i] = tr_z_zzzzz_xxx[i] * pa_y[i];

        tr_z_yzzzzz_xxy[i] = tr_z_zzzzz_xx[i] * fe_0 + tr_z_zzzzz_xxy[i] * pa_y[i];

        tr_z_yzzzzz_xxz[i] = tr_z_zzzzz_xxz[i] * pa_y[i];

        tr_z_yzzzzz_xyy[i] = 2.0 * tr_z_zzzzz_xy[i] * fe_0 + tr_z_zzzzz_xyy[i] * pa_y[i];

        tr_z_yzzzzz_xyz[i] = tr_z_zzzzz_xz[i] * fe_0 + tr_z_zzzzz_xyz[i] * pa_y[i];

        tr_z_yzzzzz_xzz[i] = tr_z_zzzzz_xzz[i] * pa_y[i];

        tr_z_yzzzzz_yyy[i] = 3.0 * tr_z_zzzzz_yy[i] * fe_0 + tr_z_zzzzz_yyy[i] * pa_y[i];

        tr_z_yzzzzz_yyz[i] = 2.0 * tr_z_zzzzz_yz[i] * fe_0 + tr_z_zzzzz_yyz[i] * pa_y[i];

        tr_z_yzzzzz_yzz[i] = tr_z_zzzzz_zz[i] * fe_0 + tr_z_zzzzz_yzz[i] * pa_y[i];

        tr_z_yzzzzz_zzz[i] = tr_z_zzzzz_zzz[i] * pa_y[i];
    }

    // Set up 830-840 components of targeted buffer : IF

    auto tr_z_zzzzzz_xxx = pbuffer.data(idx_dip_if + 830);

    auto tr_z_zzzzzz_xxy = pbuffer.data(idx_dip_if + 831);

    auto tr_z_zzzzzz_xxz = pbuffer.data(idx_dip_if + 832);

    auto tr_z_zzzzzz_xyy = pbuffer.data(idx_dip_if + 833);

    auto tr_z_zzzzzz_xyz = pbuffer.data(idx_dip_if + 834);

    auto tr_z_zzzzzz_xzz = pbuffer.data(idx_dip_if + 835);

    auto tr_z_zzzzzz_yyy = pbuffer.data(idx_dip_if + 836);

    auto tr_z_zzzzzz_yyz = pbuffer.data(idx_dip_if + 837);

    auto tr_z_zzzzzz_yzz = pbuffer.data(idx_dip_if + 838);

    auto tr_z_zzzzzz_zzz = pbuffer.data(idx_dip_if + 839);

    #pragma omp simd aligned(pa_z, tr_z_zzzz_xxx, tr_z_zzzz_xxy, tr_z_zzzz_xxz, tr_z_zzzz_xyy, tr_z_zzzz_xyz, tr_z_zzzz_xzz, tr_z_zzzz_yyy, tr_z_zzzz_yyz, tr_z_zzzz_yzz, tr_z_zzzz_zzz, tr_z_zzzzz_xx, tr_z_zzzzz_xxx, tr_z_zzzzz_xxy, tr_z_zzzzz_xxz, tr_z_zzzzz_xy, tr_z_zzzzz_xyy, tr_z_zzzzz_xyz, tr_z_zzzzz_xz, tr_z_zzzzz_xzz, tr_z_zzzzz_yy, tr_z_zzzzz_yyy, tr_z_zzzzz_yyz, tr_z_zzzzz_yz, tr_z_zzzzz_yzz, tr_z_zzzzz_zz, tr_z_zzzzz_zzz, tr_z_zzzzzz_xxx, tr_z_zzzzzz_xxy, tr_z_zzzzzz_xxz, tr_z_zzzzzz_xyy, tr_z_zzzzzz_xyz, tr_z_zzzzzz_xzz, tr_z_zzzzzz_yyy, tr_z_zzzzzz_yyz, tr_z_zzzzzz_yzz, tr_z_zzzzzz_zzz, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xzz, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yzz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzzz_xxx[i] = 5.0 * tr_z_zzzz_xxx[i] * fe_0 + ts_zzzzz_xxx[i] * fe_0 + tr_z_zzzzz_xxx[i] * pa_z[i];

        tr_z_zzzzzz_xxy[i] = 5.0 * tr_z_zzzz_xxy[i] * fe_0 + ts_zzzzz_xxy[i] * fe_0 + tr_z_zzzzz_xxy[i] * pa_z[i];

        tr_z_zzzzzz_xxz[i] = 5.0 * tr_z_zzzz_xxz[i] * fe_0 + tr_z_zzzzz_xx[i] * fe_0 + ts_zzzzz_xxz[i] * fe_0 + tr_z_zzzzz_xxz[i] * pa_z[i];

        tr_z_zzzzzz_xyy[i] = 5.0 * tr_z_zzzz_xyy[i] * fe_0 + ts_zzzzz_xyy[i] * fe_0 + tr_z_zzzzz_xyy[i] * pa_z[i];

        tr_z_zzzzzz_xyz[i] = 5.0 * tr_z_zzzz_xyz[i] * fe_0 + tr_z_zzzzz_xy[i] * fe_0 + ts_zzzzz_xyz[i] * fe_0 + tr_z_zzzzz_xyz[i] * pa_z[i];

        tr_z_zzzzzz_xzz[i] = 5.0 * tr_z_zzzz_xzz[i] * fe_0 + 2.0 * tr_z_zzzzz_xz[i] * fe_0 + ts_zzzzz_xzz[i] * fe_0 + tr_z_zzzzz_xzz[i] * pa_z[i];

        tr_z_zzzzzz_yyy[i] = 5.0 * tr_z_zzzz_yyy[i] * fe_0 + ts_zzzzz_yyy[i] * fe_0 + tr_z_zzzzz_yyy[i] * pa_z[i];

        tr_z_zzzzzz_yyz[i] = 5.0 * tr_z_zzzz_yyz[i] * fe_0 + tr_z_zzzzz_yy[i] * fe_0 + ts_zzzzz_yyz[i] * fe_0 + tr_z_zzzzz_yyz[i] * pa_z[i];

        tr_z_zzzzzz_yzz[i] = 5.0 * tr_z_zzzz_yzz[i] * fe_0 + 2.0 * tr_z_zzzzz_yz[i] * fe_0 + ts_zzzzz_yzz[i] * fe_0 + tr_z_zzzzz_yzz[i] * pa_z[i];

        tr_z_zzzzzz_zzz[i] = 5.0 * tr_z_zzzz_zzz[i] * fe_0 + 3.0 * tr_z_zzzzz_zz[i] * fe_0 + ts_zzzzz_zzz[i] * fe_0 + tr_z_zzzzz_zzz[i] * pa_z[i];
    }

}

} // diprec namespace

