#include "ElectricDipoleMomentumPrimRecID.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_id(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_id,
                                      const size_t              idx_dip_gd,
                                      const size_t              idx_dip_hp,
                                      const size_t              idx_ovl_hd,
                                      const size_t              idx_dip_hd,
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

    // Set up components of auxiliary buffer : GD

    auto tr_x_xxxx_xx = pbuffer.data(idx_dip_gd);

    auto tr_x_xxxx_xy = pbuffer.data(idx_dip_gd + 1);

    auto tr_x_xxxx_xz = pbuffer.data(idx_dip_gd + 2);

    auto tr_x_xxxx_yy = pbuffer.data(idx_dip_gd + 3);

    auto tr_x_xxxx_yz = pbuffer.data(idx_dip_gd + 4);

    auto tr_x_xxxx_zz = pbuffer.data(idx_dip_gd + 5);

    auto tr_x_xxxy_xx = pbuffer.data(idx_dip_gd + 6);

    auto tr_x_xxxy_xy = pbuffer.data(idx_dip_gd + 7);

    auto tr_x_xxxy_xz = pbuffer.data(idx_dip_gd + 8);

    auto tr_x_xxxy_zz = pbuffer.data(idx_dip_gd + 11);

    auto tr_x_xxxz_xx = pbuffer.data(idx_dip_gd + 12);

    auto tr_x_xxxz_xy = pbuffer.data(idx_dip_gd + 13);

    auto tr_x_xxxz_xz = pbuffer.data(idx_dip_gd + 14);

    auto tr_x_xxxz_yy = pbuffer.data(idx_dip_gd + 15);

    auto tr_x_xxxz_zz = pbuffer.data(idx_dip_gd + 17);

    auto tr_x_xxyy_xx = pbuffer.data(idx_dip_gd + 18);

    auto tr_x_xxyy_xy = pbuffer.data(idx_dip_gd + 19);

    auto tr_x_xxyy_xz = pbuffer.data(idx_dip_gd + 20);

    auto tr_x_xxyy_yy = pbuffer.data(idx_dip_gd + 21);

    auto tr_x_xxyy_yz = pbuffer.data(idx_dip_gd + 22);

    auto tr_x_xxyy_zz = pbuffer.data(idx_dip_gd + 23);

    auto tr_x_xxyz_xz = pbuffer.data(idx_dip_gd + 26);

    auto tr_x_xxyz_zz = pbuffer.data(idx_dip_gd + 29);

    auto tr_x_xxzz_xx = pbuffer.data(idx_dip_gd + 30);

    auto tr_x_xxzz_xy = pbuffer.data(idx_dip_gd + 31);

    auto tr_x_xxzz_xz = pbuffer.data(idx_dip_gd + 32);

    auto tr_x_xxzz_yy = pbuffer.data(idx_dip_gd + 33);

    auto tr_x_xxzz_yz = pbuffer.data(idx_dip_gd + 34);

    auto tr_x_xxzz_zz = pbuffer.data(idx_dip_gd + 35);

    auto tr_x_xyyy_xx = pbuffer.data(idx_dip_gd + 36);

    auto tr_x_xyyy_xy = pbuffer.data(idx_dip_gd + 37);

    auto tr_x_xyyy_xz = pbuffer.data(idx_dip_gd + 38);

    auto tr_x_xyyy_yy = pbuffer.data(idx_dip_gd + 39);

    auto tr_x_xyyy_yz = pbuffer.data(idx_dip_gd + 40);

    auto tr_x_xyyz_xy = pbuffer.data(idx_dip_gd + 43);

    auto tr_x_xyyz_xz = pbuffer.data(idx_dip_gd + 44);

    auto tr_x_xyzz_xx = pbuffer.data(idx_dip_gd + 48);

    auto tr_x_xyzz_xz = pbuffer.data(idx_dip_gd + 50);

    auto tr_x_xzzz_xx = pbuffer.data(idx_dip_gd + 54);

    auto tr_x_xzzz_xy = pbuffer.data(idx_dip_gd + 55);

    auto tr_x_xzzz_xz = pbuffer.data(idx_dip_gd + 56);

    auto tr_x_xzzz_yz = pbuffer.data(idx_dip_gd + 58);

    auto tr_x_xzzz_zz = pbuffer.data(idx_dip_gd + 59);

    auto tr_x_yyyy_xx = pbuffer.data(idx_dip_gd + 60);

    auto tr_x_yyyy_xy = pbuffer.data(idx_dip_gd + 61);

    auto tr_x_yyyy_xz = pbuffer.data(idx_dip_gd + 62);

    auto tr_x_yyyy_yy = pbuffer.data(idx_dip_gd + 63);

    auto tr_x_yyyy_yz = pbuffer.data(idx_dip_gd + 64);

    auto tr_x_yyyy_zz = pbuffer.data(idx_dip_gd + 65);

    auto tr_x_yyyz_xy = pbuffer.data(idx_dip_gd + 67);

    auto tr_x_yyyz_xz = pbuffer.data(idx_dip_gd + 68);

    auto tr_x_yyyz_yy = pbuffer.data(idx_dip_gd + 69);

    auto tr_x_yyyz_zz = pbuffer.data(idx_dip_gd + 71);

    auto tr_x_yyzz_xx = pbuffer.data(idx_dip_gd + 72);

    auto tr_x_yyzz_xy = pbuffer.data(idx_dip_gd + 73);

    auto tr_x_yyzz_xz = pbuffer.data(idx_dip_gd + 74);

    auto tr_x_yyzz_yy = pbuffer.data(idx_dip_gd + 75);

    auto tr_x_yyzz_yz = pbuffer.data(idx_dip_gd + 76);

    auto tr_x_yyzz_zz = pbuffer.data(idx_dip_gd + 77);

    auto tr_x_yzzz_xx = pbuffer.data(idx_dip_gd + 78);

    auto tr_x_yzzz_xz = pbuffer.data(idx_dip_gd + 80);

    auto tr_x_yzzz_yz = pbuffer.data(idx_dip_gd + 82);

    auto tr_x_yzzz_zz = pbuffer.data(idx_dip_gd + 83);

    auto tr_x_zzzz_xx = pbuffer.data(idx_dip_gd + 84);

    auto tr_x_zzzz_xy = pbuffer.data(idx_dip_gd + 85);

    auto tr_x_zzzz_xz = pbuffer.data(idx_dip_gd + 86);

    auto tr_x_zzzz_yy = pbuffer.data(idx_dip_gd + 87);

    auto tr_x_zzzz_yz = pbuffer.data(idx_dip_gd + 88);

    auto tr_x_zzzz_zz = pbuffer.data(idx_dip_gd + 89);

    auto tr_y_xxxx_xx = pbuffer.data(idx_dip_gd + 90);

    auto tr_y_xxxx_xy = pbuffer.data(idx_dip_gd + 91);

    auto tr_y_xxxx_xz = pbuffer.data(idx_dip_gd + 92);

    auto tr_y_xxxx_yy = pbuffer.data(idx_dip_gd + 93);

    auto tr_y_xxxx_yz = pbuffer.data(idx_dip_gd + 94);

    auto tr_y_xxxx_zz = pbuffer.data(idx_dip_gd + 95);

    auto tr_y_xxxy_xy = pbuffer.data(idx_dip_gd + 97);

    auto tr_y_xxxy_yy = pbuffer.data(idx_dip_gd + 99);

    auto tr_y_xxxy_yz = pbuffer.data(idx_dip_gd + 100);

    auto tr_y_xxxy_zz = pbuffer.data(idx_dip_gd + 101);

    auto tr_y_xxxz_xx = pbuffer.data(idx_dip_gd + 102);

    auto tr_y_xxxz_xy = pbuffer.data(idx_dip_gd + 103);

    auto tr_y_xxxz_yz = pbuffer.data(idx_dip_gd + 106);

    auto tr_y_xxxz_zz = pbuffer.data(idx_dip_gd + 107);

    auto tr_y_xxyy_xx = pbuffer.data(idx_dip_gd + 108);

    auto tr_y_xxyy_xy = pbuffer.data(idx_dip_gd + 109);

    auto tr_y_xxyy_xz = pbuffer.data(idx_dip_gd + 110);

    auto tr_y_xxyy_yy = pbuffer.data(idx_dip_gd + 111);

    auto tr_y_xxyy_yz = pbuffer.data(idx_dip_gd + 112);

    auto tr_y_xxyy_zz = pbuffer.data(idx_dip_gd + 113);

    auto tr_y_xxyz_xy = pbuffer.data(idx_dip_gd + 115);

    auto tr_y_xxyz_yz = pbuffer.data(idx_dip_gd + 118);

    auto tr_y_xxyz_zz = pbuffer.data(idx_dip_gd + 119);

    auto tr_y_xxzz_xx = pbuffer.data(idx_dip_gd + 120);

    auto tr_y_xxzz_xy = pbuffer.data(idx_dip_gd + 121);

    auto tr_y_xxzz_xz = pbuffer.data(idx_dip_gd + 122);

    auto tr_y_xxzz_yy = pbuffer.data(idx_dip_gd + 123);

    auto tr_y_xxzz_yz = pbuffer.data(idx_dip_gd + 124);

    auto tr_y_xxzz_zz = pbuffer.data(idx_dip_gd + 125);

    auto tr_y_xyyy_xx = pbuffer.data(idx_dip_gd + 126);

    auto tr_y_xyyy_xy = pbuffer.data(idx_dip_gd + 127);

    auto tr_y_xyyy_xz = pbuffer.data(idx_dip_gd + 128);

    auto tr_y_xyyy_yy = pbuffer.data(idx_dip_gd + 129);

    auto tr_y_xyyy_yz = pbuffer.data(idx_dip_gd + 130);

    auto tr_y_xyyy_zz = pbuffer.data(idx_dip_gd + 131);

    auto tr_y_xyyz_yz = pbuffer.data(idx_dip_gd + 136);

    auto tr_y_xyyz_zz = pbuffer.data(idx_dip_gd + 137);

    auto tr_y_xyzz_yy = pbuffer.data(idx_dip_gd + 141);

    auto tr_y_xyzz_yz = pbuffer.data(idx_dip_gd + 142);

    auto tr_y_xyzz_zz = pbuffer.data(idx_dip_gd + 143);

    auto tr_y_xzzz_xz = pbuffer.data(idx_dip_gd + 146);

    auto tr_y_xzzz_yy = pbuffer.data(idx_dip_gd + 147);

    auto tr_y_xzzz_yz = pbuffer.data(idx_dip_gd + 148);

    auto tr_y_xzzz_zz = pbuffer.data(idx_dip_gd + 149);

    auto tr_y_yyyy_xx = pbuffer.data(idx_dip_gd + 150);

    auto tr_y_yyyy_xy = pbuffer.data(idx_dip_gd + 151);

    auto tr_y_yyyy_xz = pbuffer.data(idx_dip_gd + 152);

    auto tr_y_yyyy_yy = pbuffer.data(idx_dip_gd + 153);

    auto tr_y_yyyy_yz = pbuffer.data(idx_dip_gd + 154);

    auto tr_y_yyyy_zz = pbuffer.data(idx_dip_gd + 155);

    auto tr_y_yyyz_xx = pbuffer.data(idx_dip_gd + 156);

    auto tr_y_yyyz_xy = pbuffer.data(idx_dip_gd + 157);

    auto tr_y_yyyz_yy = pbuffer.data(idx_dip_gd + 159);

    auto tr_y_yyyz_yz = pbuffer.data(idx_dip_gd + 160);

    auto tr_y_yyyz_zz = pbuffer.data(idx_dip_gd + 161);

    auto tr_y_yyzz_xx = pbuffer.data(idx_dip_gd + 162);

    auto tr_y_yyzz_xy = pbuffer.data(idx_dip_gd + 163);

    auto tr_y_yyzz_xz = pbuffer.data(idx_dip_gd + 164);

    auto tr_y_yyzz_yy = pbuffer.data(idx_dip_gd + 165);

    auto tr_y_yyzz_yz = pbuffer.data(idx_dip_gd + 166);

    auto tr_y_yyzz_zz = pbuffer.data(idx_dip_gd + 167);

    auto tr_y_yzzz_xy = pbuffer.data(idx_dip_gd + 169);

    auto tr_y_yzzz_xz = pbuffer.data(idx_dip_gd + 170);

    auto tr_y_yzzz_yy = pbuffer.data(idx_dip_gd + 171);

    auto tr_y_yzzz_yz = pbuffer.data(idx_dip_gd + 172);

    auto tr_y_yzzz_zz = pbuffer.data(idx_dip_gd + 173);

    auto tr_y_zzzz_xx = pbuffer.data(idx_dip_gd + 174);

    auto tr_y_zzzz_xy = pbuffer.data(idx_dip_gd + 175);

    auto tr_y_zzzz_xz = pbuffer.data(idx_dip_gd + 176);

    auto tr_y_zzzz_yy = pbuffer.data(idx_dip_gd + 177);

    auto tr_y_zzzz_yz = pbuffer.data(idx_dip_gd + 178);

    auto tr_y_zzzz_zz = pbuffer.data(idx_dip_gd + 179);

    auto tr_z_xxxx_xx = pbuffer.data(idx_dip_gd + 180);

    auto tr_z_xxxx_xy = pbuffer.data(idx_dip_gd + 181);

    auto tr_z_xxxx_xz = pbuffer.data(idx_dip_gd + 182);

    auto tr_z_xxxx_yy = pbuffer.data(idx_dip_gd + 183);

    auto tr_z_xxxx_yz = pbuffer.data(idx_dip_gd + 184);

    auto tr_z_xxxx_zz = pbuffer.data(idx_dip_gd + 185);

    auto tr_z_xxxy_xx = pbuffer.data(idx_dip_gd + 186);

    auto tr_z_xxxy_xz = pbuffer.data(idx_dip_gd + 188);

    auto tr_z_xxxy_yy = pbuffer.data(idx_dip_gd + 189);

    auto tr_z_xxxy_yz = pbuffer.data(idx_dip_gd + 190);

    auto tr_z_xxxz_xx = pbuffer.data(idx_dip_gd + 192);

    auto tr_z_xxxz_xz = pbuffer.data(idx_dip_gd + 194);

    auto tr_z_xxxz_yy = pbuffer.data(idx_dip_gd + 195);

    auto tr_z_xxxz_yz = pbuffer.data(idx_dip_gd + 196);

    auto tr_z_xxxz_zz = pbuffer.data(idx_dip_gd + 197);

    auto tr_z_xxyy_xx = pbuffer.data(idx_dip_gd + 198);

    auto tr_z_xxyy_xy = pbuffer.data(idx_dip_gd + 199);

    auto tr_z_xxyy_xz = pbuffer.data(idx_dip_gd + 200);

    auto tr_z_xxyy_yy = pbuffer.data(idx_dip_gd + 201);

    auto tr_z_xxyy_yz = pbuffer.data(idx_dip_gd + 202);

    auto tr_z_xxyy_zz = pbuffer.data(idx_dip_gd + 203);

    auto tr_z_xxyz_xx = pbuffer.data(idx_dip_gd + 204);

    auto tr_z_xxyz_xz = pbuffer.data(idx_dip_gd + 206);

    auto tr_z_xxyz_yy = pbuffer.data(idx_dip_gd + 207);

    auto tr_z_xxyz_yz = pbuffer.data(idx_dip_gd + 208);

    auto tr_z_xxzz_xx = pbuffer.data(idx_dip_gd + 210);

    auto tr_z_xxzz_xy = pbuffer.data(idx_dip_gd + 211);

    auto tr_z_xxzz_xz = pbuffer.data(idx_dip_gd + 212);

    auto tr_z_xxzz_yy = pbuffer.data(idx_dip_gd + 213);

    auto tr_z_xxzz_yz = pbuffer.data(idx_dip_gd + 214);

    auto tr_z_xxzz_zz = pbuffer.data(idx_dip_gd + 215);

    auto tr_z_xyyy_xy = pbuffer.data(idx_dip_gd + 217);

    auto tr_z_xyyy_yy = pbuffer.data(idx_dip_gd + 219);

    auto tr_z_xyyy_yz = pbuffer.data(idx_dip_gd + 220);

    auto tr_z_xyyy_zz = pbuffer.data(idx_dip_gd + 221);

    auto tr_z_xyyz_yy = pbuffer.data(idx_dip_gd + 225);

    auto tr_z_xyyz_yz = pbuffer.data(idx_dip_gd + 226);

    auto tr_z_xyyz_zz = pbuffer.data(idx_dip_gd + 227);

    auto tr_z_xyzz_yy = pbuffer.data(idx_dip_gd + 231);

    auto tr_z_xyzz_yz = pbuffer.data(idx_dip_gd + 232);

    auto tr_z_xzzz_xx = pbuffer.data(idx_dip_gd + 234);

    auto tr_z_xzzz_xy = pbuffer.data(idx_dip_gd + 235);

    auto tr_z_xzzz_xz = pbuffer.data(idx_dip_gd + 236);

    auto tr_z_xzzz_yy = pbuffer.data(idx_dip_gd + 237);

    auto tr_z_xzzz_yz = pbuffer.data(idx_dip_gd + 238);

    auto tr_z_xzzz_zz = pbuffer.data(idx_dip_gd + 239);

    auto tr_z_yyyy_xx = pbuffer.data(idx_dip_gd + 240);

    auto tr_z_yyyy_xy = pbuffer.data(idx_dip_gd + 241);

    auto tr_z_yyyy_xz = pbuffer.data(idx_dip_gd + 242);

    auto tr_z_yyyy_yy = pbuffer.data(idx_dip_gd + 243);

    auto tr_z_yyyy_yz = pbuffer.data(idx_dip_gd + 244);

    auto tr_z_yyyy_zz = pbuffer.data(idx_dip_gd + 245);

    auto tr_z_yyyz_xx = pbuffer.data(idx_dip_gd + 246);

    auto tr_z_yyyz_xz = pbuffer.data(idx_dip_gd + 248);

    auto tr_z_yyyz_yy = pbuffer.data(idx_dip_gd + 249);

    auto tr_z_yyyz_yz = pbuffer.data(idx_dip_gd + 250);

    auto tr_z_yyyz_zz = pbuffer.data(idx_dip_gd + 251);

    auto tr_z_yyzz_xx = pbuffer.data(idx_dip_gd + 252);

    auto tr_z_yyzz_xy = pbuffer.data(idx_dip_gd + 253);

    auto tr_z_yyzz_xz = pbuffer.data(idx_dip_gd + 254);

    auto tr_z_yyzz_yy = pbuffer.data(idx_dip_gd + 255);

    auto tr_z_yyzz_yz = pbuffer.data(idx_dip_gd + 256);

    auto tr_z_yyzz_zz = pbuffer.data(idx_dip_gd + 257);

    auto tr_z_yzzz_xx = pbuffer.data(idx_dip_gd + 258);

    auto tr_z_yzzz_xy = pbuffer.data(idx_dip_gd + 259);

    auto tr_z_yzzz_xz = pbuffer.data(idx_dip_gd + 260);

    auto tr_z_yzzz_yy = pbuffer.data(idx_dip_gd + 261);

    auto tr_z_yzzz_yz = pbuffer.data(idx_dip_gd + 262);

    auto tr_z_yzzz_zz = pbuffer.data(idx_dip_gd + 263);

    auto tr_z_zzzz_xx = pbuffer.data(idx_dip_gd + 264);

    auto tr_z_zzzz_xy = pbuffer.data(idx_dip_gd + 265);

    auto tr_z_zzzz_xz = pbuffer.data(idx_dip_gd + 266);

    auto tr_z_zzzz_yy = pbuffer.data(idx_dip_gd + 267);

    auto tr_z_zzzz_yz = pbuffer.data(idx_dip_gd + 268);

    auto tr_z_zzzz_zz = pbuffer.data(idx_dip_gd + 269);

    // Set up components of auxiliary buffer : HP

    auto tr_x_xxxxx_x = pbuffer.data(idx_dip_hp);

    auto tr_x_xxxxx_y = pbuffer.data(idx_dip_hp + 1);

    auto tr_x_xxxxx_z = pbuffer.data(idx_dip_hp + 2);

    auto tr_x_xxxxy_x = pbuffer.data(idx_dip_hp + 3);

    auto tr_x_xxxxz_x = pbuffer.data(idx_dip_hp + 6);

    auto tr_x_xxxxz_z = pbuffer.data(idx_dip_hp + 8);

    auto tr_x_xxxyy_x = pbuffer.data(idx_dip_hp + 9);

    auto tr_x_xxxyy_y = pbuffer.data(idx_dip_hp + 10);

    auto tr_x_xxxzz_x = pbuffer.data(idx_dip_hp + 15);

    auto tr_x_xxxzz_y = pbuffer.data(idx_dip_hp + 16);

    auto tr_x_xxxzz_z = pbuffer.data(idx_dip_hp + 17);

    auto tr_x_xxyyy_x = pbuffer.data(idx_dip_hp + 18);

    auto tr_x_xxyyy_y = pbuffer.data(idx_dip_hp + 19);

    auto tr_x_xxzzz_x = pbuffer.data(idx_dip_hp + 27);

    auto tr_x_xxzzz_y = pbuffer.data(idx_dip_hp + 28);

    auto tr_x_xxzzz_z = pbuffer.data(idx_dip_hp + 29);

    auto tr_x_xzzzz_x = pbuffer.data(idx_dip_hp + 42);

    auto tr_x_yyyyy_x = pbuffer.data(idx_dip_hp + 45);

    auto tr_x_yyyyy_y = pbuffer.data(idx_dip_hp + 46);

    auto tr_x_yyyyy_z = pbuffer.data(idx_dip_hp + 47);

    auto tr_x_yyyzz_z = pbuffer.data(idx_dip_hp + 53);

    auto tr_x_yyzzz_z = pbuffer.data(idx_dip_hp + 56);

    auto tr_x_yzzzz_z = pbuffer.data(idx_dip_hp + 59);

    auto tr_x_zzzzz_x = pbuffer.data(idx_dip_hp + 60);

    auto tr_x_zzzzz_y = pbuffer.data(idx_dip_hp + 61);

    auto tr_x_zzzzz_z = pbuffer.data(idx_dip_hp + 62);

    auto tr_y_xxxxx_x = pbuffer.data(idx_dip_hp + 63);

    auto tr_y_xxxxx_y = pbuffer.data(idx_dip_hp + 64);

    auto tr_y_xxxxx_z = pbuffer.data(idx_dip_hp + 65);

    auto tr_y_xxxxy_y = pbuffer.data(idx_dip_hp + 67);

    auto tr_y_xxxyy_x = pbuffer.data(idx_dip_hp + 72);

    auto tr_y_xxxyy_y = pbuffer.data(idx_dip_hp + 73);

    auto tr_y_xxxyy_z = pbuffer.data(idx_dip_hp + 74);

    auto tr_y_xxxzz_z = pbuffer.data(idx_dip_hp + 80);

    auto tr_y_xxyyy_x = pbuffer.data(idx_dip_hp + 81);

    auto tr_y_xxyyy_y = pbuffer.data(idx_dip_hp + 82);

    auto tr_y_xxyyy_z = pbuffer.data(idx_dip_hp + 83);

    auto tr_y_xxzzz_z = pbuffer.data(idx_dip_hp + 92);

    auto tr_y_xyyyy_x = pbuffer.data(idx_dip_hp + 93);

    auto tr_y_xyyyy_y = pbuffer.data(idx_dip_hp + 94);

    auto tr_y_xyyyy_z = pbuffer.data(idx_dip_hp + 95);

    auto tr_y_xyyzz_z = pbuffer.data(idx_dip_hp + 101);

    auto tr_y_xzzzz_z = pbuffer.data(idx_dip_hp + 107);

    auto tr_y_yyyyy_x = pbuffer.data(idx_dip_hp + 108);

    auto tr_y_yyyyy_y = pbuffer.data(idx_dip_hp + 109);

    auto tr_y_yyyyy_z = pbuffer.data(idx_dip_hp + 110);

    auto tr_y_yyyyz_y = pbuffer.data(idx_dip_hp + 112);

    auto tr_y_yyyyz_z = pbuffer.data(idx_dip_hp + 113);

    auto tr_y_yyyzz_x = pbuffer.data(idx_dip_hp + 114);

    auto tr_y_yyyzz_y = pbuffer.data(idx_dip_hp + 115);

    auto tr_y_yyyzz_z = pbuffer.data(idx_dip_hp + 116);

    auto tr_y_yyzzz_x = pbuffer.data(idx_dip_hp + 117);

    auto tr_y_yyzzz_y = pbuffer.data(idx_dip_hp + 118);

    auto tr_y_yyzzz_z = pbuffer.data(idx_dip_hp + 119);

    auto tr_y_yzzzz_x = pbuffer.data(idx_dip_hp + 120);

    auto tr_y_yzzzz_y = pbuffer.data(idx_dip_hp + 121);

    auto tr_y_yzzzz_z = pbuffer.data(idx_dip_hp + 122);

    auto tr_y_zzzzz_x = pbuffer.data(idx_dip_hp + 123);

    auto tr_y_zzzzz_y = pbuffer.data(idx_dip_hp + 124);

    auto tr_y_zzzzz_z = pbuffer.data(idx_dip_hp + 125);

    auto tr_z_xxxxx_x = pbuffer.data(idx_dip_hp + 126);

    auto tr_z_xxxxx_y = pbuffer.data(idx_dip_hp + 127);

    auto tr_z_xxxxx_z = pbuffer.data(idx_dip_hp + 128);

    auto tr_z_xxxxz_x = pbuffer.data(idx_dip_hp + 132);

    auto tr_z_xxxxz_z = pbuffer.data(idx_dip_hp + 134);

    auto tr_z_xxxyy_y = pbuffer.data(idx_dip_hp + 136);

    auto tr_z_xxxzz_x = pbuffer.data(idx_dip_hp + 141);

    auto tr_z_xxxzz_y = pbuffer.data(idx_dip_hp + 142);

    auto tr_z_xxxzz_z = pbuffer.data(idx_dip_hp + 143);

    auto tr_z_xxyyy_y = pbuffer.data(idx_dip_hp + 145);

    auto tr_z_xxzzz_x = pbuffer.data(idx_dip_hp + 153);

    auto tr_z_xxzzz_y = pbuffer.data(idx_dip_hp + 154);

    auto tr_z_xxzzz_z = pbuffer.data(idx_dip_hp + 155);

    auto tr_z_xyyyy_y = pbuffer.data(idx_dip_hp + 157);

    auto tr_z_xyyzz_y = pbuffer.data(idx_dip_hp + 163);

    auto tr_z_xzzzz_x = pbuffer.data(idx_dip_hp + 168);

    auto tr_z_xzzzz_y = pbuffer.data(idx_dip_hp + 169);

    auto tr_z_xzzzz_z = pbuffer.data(idx_dip_hp + 170);

    auto tr_z_yyyyy_x = pbuffer.data(idx_dip_hp + 171);

    auto tr_z_yyyyy_y = pbuffer.data(idx_dip_hp + 172);

    auto tr_z_yyyyy_z = pbuffer.data(idx_dip_hp + 173);

    auto tr_z_yyyyz_x = pbuffer.data(idx_dip_hp + 174);

    auto tr_z_yyyyz_y = pbuffer.data(idx_dip_hp + 175);

    auto tr_z_yyyyz_z = pbuffer.data(idx_dip_hp + 176);

    auto tr_z_yyyzz_x = pbuffer.data(idx_dip_hp + 177);

    auto tr_z_yyyzz_y = pbuffer.data(idx_dip_hp + 178);

    auto tr_z_yyyzz_z = pbuffer.data(idx_dip_hp + 179);

    auto tr_z_yyzzz_x = pbuffer.data(idx_dip_hp + 180);

    auto tr_z_yyzzz_y = pbuffer.data(idx_dip_hp + 181);

    auto tr_z_yyzzz_z = pbuffer.data(idx_dip_hp + 182);

    auto tr_z_yzzzz_x = pbuffer.data(idx_dip_hp + 183);

    auto tr_z_yzzzz_y = pbuffer.data(idx_dip_hp + 184);

    auto tr_z_yzzzz_z = pbuffer.data(idx_dip_hp + 185);

    auto tr_z_zzzzz_x = pbuffer.data(idx_dip_hp + 186);

    auto tr_z_zzzzz_y = pbuffer.data(idx_dip_hp + 187);

    auto tr_z_zzzzz_z = pbuffer.data(idx_dip_hp + 188);

    // Set up components of auxiliary buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_ovl_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_ovl_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_ovl_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_ovl_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_ovl_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_ovl_hd + 5);

    auto ts_xxxxz_xz = pbuffer.data(idx_ovl_hd + 14);

    auto ts_xxxyy_xy = pbuffer.data(idx_ovl_hd + 19);

    auto ts_xxxyy_yy = pbuffer.data(idx_ovl_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_ovl_hd + 22);

    auto ts_xxxzz_xx = pbuffer.data(idx_ovl_hd + 30);

    auto ts_xxxzz_xz = pbuffer.data(idx_ovl_hd + 32);

    auto ts_xxxzz_yz = pbuffer.data(idx_ovl_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_ovl_hd + 35);

    auto ts_xxyyy_xy = pbuffer.data(idx_ovl_hd + 37);

    auto ts_xxyyy_yy = pbuffer.data(idx_ovl_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_ovl_hd + 40);

    auto ts_xxzzz_xx = pbuffer.data(idx_ovl_hd + 54);

    auto ts_xxzzz_xz = pbuffer.data(idx_ovl_hd + 56);

    auto ts_xxzzz_yz = pbuffer.data(idx_ovl_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_ovl_hd + 59);

    auto ts_xyyyy_yy = pbuffer.data(idx_ovl_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_ovl_hd + 64);

    auto ts_xyyzz_yz = pbuffer.data(idx_ovl_hd + 76);

    auto ts_xzzzz_yz = pbuffer.data(idx_ovl_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_ovl_hd + 89);

    auto ts_yyyyy_xx = pbuffer.data(idx_ovl_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_ovl_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_ovl_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_ovl_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_ovl_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_ovl_hd + 95);

    auto ts_yyyyz_yz = pbuffer.data(idx_ovl_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_ovl_hd + 101);

    auto ts_yyyzz_xz = pbuffer.data(idx_ovl_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_ovl_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_ovl_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_ovl_hd + 107);

    auto ts_yyzzz_xz = pbuffer.data(idx_ovl_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_ovl_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_ovl_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_ovl_hd + 113);

    auto ts_yzzzz_xz = pbuffer.data(idx_ovl_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_ovl_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_ovl_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_ovl_hd + 119);

    auto ts_zzzzz_xx = pbuffer.data(idx_ovl_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_ovl_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_ovl_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_ovl_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_ovl_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_ovl_hd + 125);

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

    auto tr_x_xxxxy_yy = pbuffer.data(idx_dip_hd + 9);

    auto tr_x_xxxxy_zz = pbuffer.data(idx_dip_hd + 11);

    auto tr_x_xxxxz_xx = pbuffer.data(idx_dip_hd + 12);

    auto tr_x_xxxxz_xy = pbuffer.data(idx_dip_hd + 13);

    auto tr_x_xxxxz_xz = pbuffer.data(idx_dip_hd + 14);

    auto tr_x_xxxxz_yy = pbuffer.data(idx_dip_hd + 15);

    auto tr_x_xxxxz_yz = pbuffer.data(idx_dip_hd + 16);

    auto tr_x_xxxxz_zz = pbuffer.data(idx_dip_hd + 17);

    auto tr_x_xxxyy_xx = pbuffer.data(idx_dip_hd + 18);

    auto tr_x_xxxyy_xy = pbuffer.data(idx_dip_hd + 19);

    auto tr_x_xxxyy_xz = pbuffer.data(idx_dip_hd + 20);

    auto tr_x_xxxyy_yy = pbuffer.data(idx_dip_hd + 21);

    auto tr_x_xxxyy_yz = pbuffer.data(idx_dip_hd + 22);

    auto tr_x_xxxyy_zz = pbuffer.data(idx_dip_hd + 23);

    auto tr_x_xxxyz_xz = pbuffer.data(idx_dip_hd + 26);

    auto tr_x_xxxyz_zz = pbuffer.data(idx_dip_hd + 29);

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

    auto tr_x_xxyyy_zz = pbuffer.data(idx_dip_hd + 41);

    auto tr_x_xxyyz_xy = pbuffer.data(idx_dip_hd + 43);

    auto tr_x_xxyyz_xz = pbuffer.data(idx_dip_hd + 44);

    auto tr_x_xxyyz_yy = pbuffer.data(idx_dip_hd + 45);

    auto tr_x_xxyyz_zz = pbuffer.data(idx_dip_hd + 47);

    auto tr_x_xxyzz_xx = pbuffer.data(idx_dip_hd + 48);

    auto tr_x_xxyzz_xz = pbuffer.data(idx_dip_hd + 50);

    auto tr_x_xxyzz_zz = pbuffer.data(idx_dip_hd + 53);

    auto tr_x_xxzzz_xx = pbuffer.data(idx_dip_hd + 54);

    auto tr_x_xxzzz_xy = pbuffer.data(idx_dip_hd + 55);

    auto tr_x_xxzzz_xz = pbuffer.data(idx_dip_hd + 56);

    auto tr_x_xxzzz_yy = pbuffer.data(idx_dip_hd + 57);

    auto tr_x_xxzzz_yz = pbuffer.data(idx_dip_hd + 58);

    auto tr_x_xxzzz_zz = pbuffer.data(idx_dip_hd + 59);

    auto tr_x_xyyyy_xx = pbuffer.data(idx_dip_hd + 60);

    auto tr_x_xyyyy_xy = pbuffer.data(idx_dip_hd + 61);

    auto tr_x_xyyyy_xz = pbuffer.data(idx_dip_hd + 62);

    auto tr_x_xyyyy_yy = pbuffer.data(idx_dip_hd + 63);

    auto tr_x_xyyyy_yz = pbuffer.data(idx_dip_hd + 64);

    auto tr_x_xyyyz_xy = pbuffer.data(idx_dip_hd + 67);

    auto tr_x_xyyyz_xz = pbuffer.data(idx_dip_hd + 68);

    auto tr_x_xyyzz_xx = pbuffer.data(idx_dip_hd + 72);

    auto tr_x_xyyzz_xy = pbuffer.data(idx_dip_hd + 73);

    auto tr_x_xyyzz_xz = pbuffer.data(idx_dip_hd + 74);

    auto tr_x_xyyzz_yz = pbuffer.data(idx_dip_hd + 76);

    auto tr_x_xyzzz_xx = pbuffer.data(idx_dip_hd + 78);

    auto tr_x_xyzzz_xz = pbuffer.data(idx_dip_hd + 80);

    auto tr_x_xzzzz_xx = pbuffer.data(idx_dip_hd + 84);

    auto tr_x_xzzzz_xy = pbuffer.data(idx_dip_hd + 85);

    auto tr_x_xzzzz_xz = pbuffer.data(idx_dip_hd + 86);

    auto tr_x_xzzzz_yz = pbuffer.data(idx_dip_hd + 88);

    auto tr_x_xzzzz_zz = pbuffer.data(idx_dip_hd + 89);

    auto tr_x_yyyyy_xx = pbuffer.data(idx_dip_hd + 90);

    auto tr_x_yyyyy_xy = pbuffer.data(idx_dip_hd + 91);

    auto tr_x_yyyyy_xz = pbuffer.data(idx_dip_hd + 92);

    auto tr_x_yyyyy_yy = pbuffer.data(idx_dip_hd + 93);

    auto tr_x_yyyyy_yz = pbuffer.data(idx_dip_hd + 94);

    auto tr_x_yyyyy_zz = pbuffer.data(idx_dip_hd + 95);

    auto tr_x_yyyyz_xy = pbuffer.data(idx_dip_hd + 97);

    auto tr_x_yyyyz_xz = pbuffer.data(idx_dip_hd + 98);

    auto tr_x_yyyyz_yy = pbuffer.data(idx_dip_hd + 99);

    auto tr_x_yyyyz_yz = pbuffer.data(idx_dip_hd + 100);

    auto tr_x_yyyyz_zz = pbuffer.data(idx_dip_hd + 101);

    auto tr_x_yyyzz_xx = pbuffer.data(idx_dip_hd + 102);

    auto tr_x_yyyzz_xy = pbuffer.data(idx_dip_hd + 103);

    auto tr_x_yyyzz_xz = pbuffer.data(idx_dip_hd + 104);

    auto tr_x_yyyzz_yy = pbuffer.data(idx_dip_hd + 105);

    auto tr_x_yyyzz_yz = pbuffer.data(idx_dip_hd + 106);

    auto tr_x_yyyzz_zz = pbuffer.data(idx_dip_hd + 107);

    auto tr_x_yyzzz_xx = pbuffer.data(idx_dip_hd + 108);

    auto tr_x_yyzzz_xy = pbuffer.data(idx_dip_hd + 109);

    auto tr_x_yyzzz_xz = pbuffer.data(idx_dip_hd + 110);

    auto tr_x_yyzzz_yy = pbuffer.data(idx_dip_hd + 111);

    auto tr_x_yyzzz_yz = pbuffer.data(idx_dip_hd + 112);

    auto tr_x_yyzzz_zz = pbuffer.data(idx_dip_hd + 113);

    auto tr_x_yzzzz_xx = pbuffer.data(idx_dip_hd + 114);

    auto tr_x_yzzzz_xz = pbuffer.data(idx_dip_hd + 116);

    auto tr_x_yzzzz_yy = pbuffer.data(idx_dip_hd + 117);

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

    auto tr_y_xxxxy_xx = pbuffer.data(idx_dip_hd + 132);

    auto tr_y_xxxxy_xy = pbuffer.data(idx_dip_hd + 133);

    auto tr_y_xxxxy_yy = pbuffer.data(idx_dip_hd + 135);

    auto tr_y_xxxxy_yz = pbuffer.data(idx_dip_hd + 136);

    auto tr_y_xxxxy_zz = pbuffer.data(idx_dip_hd + 137);

    auto tr_y_xxxxz_xx = pbuffer.data(idx_dip_hd + 138);

    auto tr_y_xxxxz_xy = pbuffer.data(idx_dip_hd + 139);

    auto tr_y_xxxxz_xz = pbuffer.data(idx_dip_hd + 140);

    auto tr_y_xxxxz_yz = pbuffer.data(idx_dip_hd + 142);

    auto tr_y_xxxxz_zz = pbuffer.data(idx_dip_hd + 143);

    auto tr_y_xxxyy_xx = pbuffer.data(idx_dip_hd + 144);

    auto tr_y_xxxyy_xy = pbuffer.data(idx_dip_hd + 145);

    auto tr_y_xxxyy_xz = pbuffer.data(idx_dip_hd + 146);

    auto tr_y_xxxyy_yy = pbuffer.data(idx_dip_hd + 147);

    auto tr_y_xxxyy_yz = pbuffer.data(idx_dip_hd + 148);

    auto tr_y_xxxyy_zz = pbuffer.data(idx_dip_hd + 149);

    auto tr_y_xxxyz_xy = pbuffer.data(idx_dip_hd + 151);

    auto tr_y_xxxyz_yz = pbuffer.data(idx_dip_hd + 154);

    auto tr_y_xxxyz_zz = pbuffer.data(idx_dip_hd + 155);

    auto tr_y_xxxzz_xx = pbuffer.data(idx_dip_hd + 156);

    auto tr_y_xxxzz_xy = pbuffer.data(idx_dip_hd + 157);

    auto tr_y_xxxzz_xz = pbuffer.data(idx_dip_hd + 158);

    auto tr_y_xxxzz_yy = pbuffer.data(idx_dip_hd + 159);

    auto tr_y_xxxzz_yz = pbuffer.data(idx_dip_hd + 160);

    auto tr_y_xxxzz_zz = pbuffer.data(idx_dip_hd + 161);

    auto tr_y_xxyyy_xx = pbuffer.data(idx_dip_hd + 162);

    auto tr_y_xxyyy_xy = pbuffer.data(idx_dip_hd + 163);

    auto tr_y_xxyyy_xz = pbuffer.data(idx_dip_hd + 164);

    auto tr_y_xxyyy_yy = pbuffer.data(idx_dip_hd + 165);

    auto tr_y_xxyyy_yz = pbuffer.data(idx_dip_hd + 166);

    auto tr_y_xxyyy_zz = pbuffer.data(idx_dip_hd + 167);

    auto tr_y_xxyyz_xx = pbuffer.data(idx_dip_hd + 168);

    auto tr_y_xxyyz_xy = pbuffer.data(idx_dip_hd + 169);

    auto tr_y_xxyyz_yz = pbuffer.data(idx_dip_hd + 172);

    auto tr_y_xxyyz_zz = pbuffer.data(idx_dip_hd + 173);

    auto tr_y_xxyzz_xy = pbuffer.data(idx_dip_hd + 175);

    auto tr_y_xxyzz_yy = pbuffer.data(idx_dip_hd + 177);

    auto tr_y_xxyzz_yz = pbuffer.data(idx_dip_hd + 178);

    auto tr_y_xxyzz_zz = pbuffer.data(idx_dip_hd + 179);

    auto tr_y_xxzzz_xx = pbuffer.data(idx_dip_hd + 180);

    auto tr_y_xxzzz_xy = pbuffer.data(idx_dip_hd + 181);

    auto tr_y_xxzzz_xz = pbuffer.data(idx_dip_hd + 182);

    auto tr_y_xxzzz_yy = pbuffer.data(idx_dip_hd + 183);

    auto tr_y_xxzzz_yz = pbuffer.data(idx_dip_hd + 184);

    auto tr_y_xxzzz_zz = pbuffer.data(idx_dip_hd + 185);

    auto tr_y_xyyyy_xx = pbuffer.data(idx_dip_hd + 186);

    auto tr_y_xyyyy_xy = pbuffer.data(idx_dip_hd + 187);

    auto tr_y_xyyyy_xz = pbuffer.data(idx_dip_hd + 188);

    auto tr_y_xyyyy_yy = pbuffer.data(idx_dip_hd + 189);

    auto tr_y_xyyyy_yz = pbuffer.data(idx_dip_hd + 190);

    auto tr_y_xyyyy_zz = pbuffer.data(idx_dip_hd + 191);

    auto tr_y_xyyyz_yz = pbuffer.data(idx_dip_hd + 196);

    auto tr_y_xyyyz_zz = pbuffer.data(idx_dip_hd + 197);

    auto tr_y_xyyzz_xz = pbuffer.data(idx_dip_hd + 200);

    auto tr_y_xyyzz_yy = pbuffer.data(idx_dip_hd + 201);

    auto tr_y_xyyzz_yz = pbuffer.data(idx_dip_hd + 202);

    auto tr_y_xyyzz_zz = pbuffer.data(idx_dip_hd + 203);

    auto tr_y_xyzzz_yy = pbuffer.data(idx_dip_hd + 207);

    auto tr_y_xyzzz_yz = pbuffer.data(idx_dip_hd + 208);

    auto tr_y_xyzzz_zz = pbuffer.data(idx_dip_hd + 209);

    auto tr_y_xzzzz_xz = pbuffer.data(idx_dip_hd + 212);

    auto tr_y_xzzzz_yy = pbuffer.data(idx_dip_hd + 213);

    auto tr_y_xzzzz_yz = pbuffer.data(idx_dip_hd + 214);

    auto tr_y_xzzzz_zz = pbuffer.data(idx_dip_hd + 215);

    auto tr_y_yyyyy_xx = pbuffer.data(idx_dip_hd + 216);

    auto tr_y_yyyyy_xy = pbuffer.data(idx_dip_hd + 217);

    auto tr_y_yyyyy_xz = pbuffer.data(idx_dip_hd + 218);

    auto tr_y_yyyyy_yy = pbuffer.data(idx_dip_hd + 219);

    auto tr_y_yyyyy_yz = pbuffer.data(idx_dip_hd + 220);

    auto tr_y_yyyyy_zz = pbuffer.data(idx_dip_hd + 221);

    auto tr_y_yyyyz_xx = pbuffer.data(idx_dip_hd + 222);

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

    auto tr_z_xxxxy_xx = pbuffer.data(idx_dip_hd + 258);

    auto tr_z_xxxxy_xz = pbuffer.data(idx_dip_hd + 260);

    auto tr_z_xxxxy_yy = pbuffer.data(idx_dip_hd + 261);

    auto tr_z_xxxxy_yz = pbuffer.data(idx_dip_hd + 262);

    auto tr_z_xxxxz_xx = pbuffer.data(idx_dip_hd + 264);

    auto tr_z_xxxxz_xy = pbuffer.data(idx_dip_hd + 265);

    auto tr_z_xxxxz_xz = pbuffer.data(idx_dip_hd + 266);

    auto tr_z_xxxxz_yy = pbuffer.data(idx_dip_hd + 267);

    auto tr_z_xxxxz_yz = pbuffer.data(idx_dip_hd + 268);

    auto tr_z_xxxxz_zz = pbuffer.data(idx_dip_hd + 269);

    auto tr_z_xxxyy_xx = pbuffer.data(idx_dip_hd + 270);

    auto tr_z_xxxyy_xy = pbuffer.data(idx_dip_hd + 271);

    auto tr_z_xxxyy_xz = pbuffer.data(idx_dip_hd + 272);

    auto tr_z_xxxyy_yy = pbuffer.data(idx_dip_hd + 273);

    auto tr_z_xxxyy_yz = pbuffer.data(idx_dip_hd + 274);

    auto tr_z_xxxyy_zz = pbuffer.data(idx_dip_hd + 275);

    auto tr_z_xxxyz_xx = pbuffer.data(idx_dip_hd + 276);

    auto tr_z_xxxyz_xz = pbuffer.data(idx_dip_hd + 278);

    auto tr_z_xxxyz_yy = pbuffer.data(idx_dip_hd + 279);

    auto tr_z_xxxyz_yz = pbuffer.data(idx_dip_hd + 280);

    auto tr_z_xxxzz_xx = pbuffer.data(idx_dip_hd + 282);

    auto tr_z_xxxzz_xy = pbuffer.data(idx_dip_hd + 283);

    auto tr_z_xxxzz_xz = pbuffer.data(idx_dip_hd + 284);

    auto tr_z_xxxzz_yy = pbuffer.data(idx_dip_hd + 285);

    auto tr_z_xxxzz_yz = pbuffer.data(idx_dip_hd + 286);

    auto tr_z_xxxzz_zz = pbuffer.data(idx_dip_hd + 287);

    auto tr_z_xxyyy_xx = pbuffer.data(idx_dip_hd + 288);

    auto tr_z_xxyyy_xy = pbuffer.data(idx_dip_hd + 289);

    auto tr_z_xxyyy_xz = pbuffer.data(idx_dip_hd + 290);

    auto tr_z_xxyyy_yy = pbuffer.data(idx_dip_hd + 291);

    auto tr_z_xxyyy_yz = pbuffer.data(idx_dip_hd + 292);

    auto tr_z_xxyyy_zz = pbuffer.data(idx_dip_hd + 293);

    auto tr_z_xxyyz_xx = pbuffer.data(idx_dip_hd + 294);

    auto tr_z_xxyyz_xz = pbuffer.data(idx_dip_hd + 296);

    auto tr_z_xxyyz_yy = pbuffer.data(idx_dip_hd + 297);

    auto tr_z_xxyyz_yz = pbuffer.data(idx_dip_hd + 298);

    auto tr_z_xxyyz_zz = pbuffer.data(idx_dip_hd + 299);

    auto tr_z_xxyzz_xx = pbuffer.data(idx_dip_hd + 300);

    auto tr_z_xxyzz_xz = pbuffer.data(idx_dip_hd + 302);

    auto tr_z_xxyzz_yy = pbuffer.data(idx_dip_hd + 303);

    auto tr_z_xxyzz_yz = pbuffer.data(idx_dip_hd + 304);

    auto tr_z_xxzzz_xx = pbuffer.data(idx_dip_hd + 306);

    auto tr_z_xxzzz_xy = pbuffer.data(idx_dip_hd + 307);

    auto tr_z_xxzzz_xz = pbuffer.data(idx_dip_hd + 308);

    auto tr_z_xxzzz_yy = pbuffer.data(idx_dip_hd + 309);

    auto tr_z_xxzzz_yz = pbuffer.data(idx_dip_hd + 310);

    auto tr_z_xxzzz_zz = pbuffer.data(idx_dip_hd + 311);

    auto tr_z_xyyyy_xy = pbuffer.data(idx_dip_hd + 313);

    auto tr_z_xyyyy_yy = pbuffer.data(idx_dip_hd + 315);

    auto tr_z_xyyyy_yz = pbuffer.data(idx_dip_hd + 316);

    auto tr_z_xyyyy_zz = pbuffer.data(idx_dip_hd + 317);

    auto tr_z_xyyyz_yy = pbuffer.data(idx_dip_hd + 321);

    auto tr_z_xyyyz_yz = pbuffer.data(idx_dip_hd + 322);

    auto tr_z_xyyyz_zz = pbuffer.data(idx_dip_hd + 323);

    auto tr_z_xyyzz_xy = pbuffer.data(idx_dip_hd + 325);

    auto tr_z_xyyzz_yy = pbuffer.data(idx_dip_hd + 327);

    auto tr_z_xyyzz_yz = pbuffer.data(idx_dip_hd + 328);

    auto tr_z_xyyzz_zz = pbuffer.data(idx_dip_hd + 329);

    auto tr_z_xyzzz_yy = pbuffer.data(idx_dip_hd + 333);

    auto tr_z_xyzzz_yz = pbuffer.data(idx_dip_hd + 334);

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

    // Set up 0-6 components of targeted buffer : ID

    auto tr_x_xxxxxx_xx = pbuffer.data(idx_dip_id);

    auto tr_x_xxxxxx_xy = pbuffer.data(idx_dip_id + 1);

    auto tr_x_xxxxxx_xz = pbuffer.data(idx_dip_id + 2);

    auto tr_x_xxxxxx_yy = pbuffer.data(idx_dip_id + 3);

    auto tr_x_xxxxxx_yz = pbuffer.data(idx_dip_id + 4);

    auto tr_x_xxxxxx_zz = pbuffer.data(idx_dip_id + 5);

#pragma omp simd aligned(pa_x,               \
                             tr_x_xxxx_xx,   \
                             tr_x_xxxx_xy,   \
                             tr_x_xxxx_xz,   \
                             tr_x_xxxx_yy,   \
                             tr_x_xxxx_yz,   \
                             tr_x_xxxx_zz,   \
                             tr_x_xxxxx_x,   \
                             tr_x_xxxxx_xx,  \
                             tr_x_xxxxx_xy,  \
                             tr_x_xxxxx_xz,  \
                             tr_x_xxxxx_y,   \
                             tr_x_xxxxx_yy,  \
                             tr_x_xxxxx_yz,  \
                             tr_x_xxxxx_z,   \
                             tr_x_xxxxx_zz,  \
                             tr_x_xxxxxx_xx, \
                             tr_x_xxxxxx_xy, \
                             tr_x_xxxxxx_xz, \
                             tr_x_xxxxxx_yy, \
                             tr_x_xxxxxx_yz, \
                             tr_x_xxxxxx_zz, \
                             ts_xxxxx_xx,    \
                             ts_xxxxx_xy,    \
                             ts_xxxxx_xz,    \
                             ts_xxxxx_yy,    \
                             ts_xxxxx_yz,    \
                             ts_xxxxx_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_xx[i] = 5.0 * tr_x_xxxx_xx[i] * fe_0 + 2.0 * tr_x_xxxxx_x[i] * fe_0 + ts_xxxxx_xx[i] * fe_0 + tr_x_xxxxx_xx[i] * pa_x[i];

        tr_x_xxxxxx_xy[i] = 5.0 * tr_x_xxxx_xy[i] * fe_0 + tr_x_xxxxx_y[i] * fe_0 + ts_xxxxx_xy[i] * fe_0 + tr_x_xxxxx_xy[i] * pa_x[i];

        tr_x_xxxxxx_xz[i] = 5.0 * tr_x_xxxx_xz[i] * fe_0 + tr_x_xxxxx_z[i] * fe_0 + ts_xxxxx_xz[i] * fe_0 + tr_x_xxxxx_xz[i] * pa_x[i];

        tr_x_xxxxxx_yy[i] = 5.0 * tr_x_xxxx_yy[i] * fe_0 + ts_xxxxx_yy[i] * fe_0 + tr_x_xxxxx_yy[i] * pa_x[i];

        tr_x_xxxxxx_yz[i] = 5.0 * tr_x_xxxx_yz[i] * fe_0 + ts_xxxxx_yz[i] * fe_0 + tr_x_xxxxx_yz[i] * pa_x[i];

        tr_x_xxxxxx_zz[i] = 5.0 * tr_x_xxxx_zz[i] * fe_0 + ts_xxxxx_zz[i] * fe_0 + tr_x_xxxxx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : ID

    auto tr_x_xxxxxy_xx = pbuffer.data(idx_dip_id + 6);

    auto tr_x_xxxxxy_xy = pbuffer.data(idx_dip_id + 7);

    auto tr_x_xxxxxy_xz = pbuffer.data(idx_dip_id + 8);

    auto tr_x_xxxxxy_yy = pbuffer.data(idx_dip_id + 9);

    auto tr_x_xxxxxy_yz = pbuffer.data(idx_dip_id + 10);

    auto tr_x_xxxxxy_zz = pbuffer.data(idx_dip_id + 11);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xxxxx_x,   \
                             tr_x_xxxxx_xx,  \
                             tr_x_xxxxx_xy,  \
                             tr_x_xxxxx_xz,  \
                             tr_x_xxxxx_y,   \
                             tr_x_xxxxx_yy,  \
                             tr_x_xxxxx_yz,  \
                             tr_x_xxxxx_z,   \
                             tr_x_xxxxx_zz,  \
                             tr_x_xxxxxy_xx, \
                             tr_x_xxxxxy_xy, \
                             tr_x_xxxxxy_xz, \
                             tr_x_xxxxxy_yy, \
                             tr_x_xxxxxy_yz, \
                             tr_x_xxxxxy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxy_xx[i] = tr_x_xxxxx_xx[i] * pa_y[i];

        tr_x_xxxxxy_xy[i] = tr_x_xxxxx_x[i] * fe_0 + tr_x_xxxxx_xy[i] * pa_y[i];

        tr_x_xxxxxy_xz[i] = tr_x_xxxxx_xz[i] * pa_y[i];

        tr_x_xxxxxy_yy[i] = 2.0 * tr_x_xxxxx_y[i] * fe_0 + tr_x_xxxxx_yy[i] * pa_y[i];

        tr_x_xxxxxy_yz[i] = tr_x_xxxxx_z[i] * fe_0 + tr_x_xxxxx_yz[i] * pa_y[i];

        tr_x_xxxxxy_zz[i] = tr_x_xxxxx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : ID

    auto tr_x_xxxxxz_xx = pbuffer.data(idx_dip_id + 12);

    auto tr_x_xxxxxz_xy = pbuffer.data(idx_dip_id + 13);

    auto tr_x_xxxxxz_xz = pbuffer.data(idx_dip_id + 14);

    auto tr_x_xxxxxz_yy = pbuffer.data(idx_dip_id + 15);

    auto tr_x_xxxxxz_yz = pbuffer.data(idx_dip_id + 16);

    auto tr_x_xxxxxz_zz = pbuffer.data(idx_dip_id + 17);

#pragma omp simd aligned(pa_z,               \
                             tr_x_xxxxx_x,   \
                             tr_x_xxxxx_xx,  \
                             tr_x_xxxxx_xy,  \
                             tr_x_xxxxx_xz,  \
                             tr_x_xxxxx_y,   \
                             tr_x_xxxxx_yy,  \
                             tr_x_xxxxx_yz,  \
                             tr_x_xxxxx_z,   \
                             tr_x_xxxxx_zz,  \
                             tr_x_xxxxxz_xx, \
                             tr_x_xxxxxz_xy, \
                             tr_x_xxxxxz_xz, \
                             tr_x_xxxxxz_yy, \
                             tr_x_xxxxxz_yz, \
                             tr_x_xxxxxz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxz_xx[i] = tr_x_xxxxx_xx[i] * pa_z[i];

        tr_x_xxxxxz_xy[i] = tr_x_xxxxx_xy[i] * pa_z[i];

        tr_x_xxxxxz_xz[i] = tr_x_xxxxx_x[i] * fe_0 + tr_x_xxxxx_xz[i] * pa_z[i];

        tr_x_xxxxxz_yy[i] = tr_x_xxxxx_yy[i] * pa_z[i];

        tr_x_xxxxxz_yz[i] = tr_x_xxxxx_y[i] * fe_0 + tr_x_xxxxx_yz[i] * pa_z[i];

        tr_x_xxxxxz_zz[i] = 2.0 * tr_x_xxxxx_z[i] * fe_0 + tr_x_xxxxx_zz[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : ID

    auto tr_x_xxxxyy_xx = pbuffer.data(idx_dip_id + 18);

    auto tr_x_xxxxyy_xy = pbuffer.data(idx_dip_id + 19);

    auto tr_x_xxxxyy_xz = pbuffer.data(idx_dip_id + 20);

    auto tr_x_xxxxyy_yy = pbuffer.data(idx_dip_id + 21);

    auto tr_x_xxxxyy_yz = pbuffer.data(idx_dip_id + 22);

    auto tr_x_xxxxyy_zz = pbuffer.data(idx_dip_id + 23);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xxxx_xx,   \
                             tr_x_xxxx_xy,   \
                             tr_x_xxxx_xz,   \
                             tr_x_xxxx_zz,   \
                             tr_x_xxxxy_x,   \
                             tr_x_xxxxy_xx,  \
                             tr_x_xxxxy_xy,  \
                             tr_x_xxxxy_xz,  \
                             tr_x_xxxxy_zz,  \
                             tr_x_xxxxyy_xx, \
                             tr_x_xxxxyy_xy, \
                             tr_x_xxxxyy_xz, \
                             tr_x_xxxxyy_yy, \
                             tr_x_xxxxyy_yz, \
                             tr_x_xxxxyy_zz, \
                             tr_x_xxxyy_yy,  \
                             tr_x_xxxyy_yz,  \
                             tr_x_xxyy_yy,   \
                             tr_x_xxyy_yz,   \
                             ts_xxxyy_yy,    \
                             ts_xxxyy_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyy_xx[i] = tr_x_xxxx_xx[i] * fe_0 + tr_x_xxxxy_xx[i] * pa_y[i];

        tr_x_xxxxyy_xy[i] = tr_x_xxxx_xy[i] * fe_0 + tr_x_xxxxy_x[i] * fe_0 + tr_x_xxxxy_xy[i] * pa_y[i];

        tr_x_xxxxyy_xz[i] = tr_x_xxxx_xz[i] * fe_0 + tr_x_xxxxy_xz[i] * pa_y[i];

        tr_x_xxxxyy_yy[i] = 3.0 * tr_x_xxyy_yy[i] * fe_0 + ts_xxxyy_yy[i] * fe_0 + tr_x_xxxyy_yy[i] * pa_x[i];

        tr_x_xxxxyy_yz[i] = 3.0 * tr_x_xxyy_yz[i] * fe_0 + ts_xxxyy_yz[i] * fe_0 + tr_x_xxxyy_yz[i] * pa_x[i];

        tr_x_xxxxyy_zz[i] = tr_x_xxxx_zz[i] * fe_0 + tr_x_xxxxy_zz[i] * pa_y[i];
    }

    // Set up 24-30 components of targeted buffer : ID

    auto tr_x_xxxxyz_xx = pbuffer.data(idx_dip_id + 24);

    auto tr_x_xxxxyz_xy = pbuffer.data(idx_dip_id + 25);

    auto tr_x_xxxxyz_xz = pbuffer.data(idx_dip_id + 26);

    auto tr_x_xxxxyz_yy = pbuffer.data(idx_dip_id + 27);

    auto tr_x_xxxxyz_yz = pbuffer.data(idx_dip_id + 28);

    auto tr_x_xxxxyz_zz = pbuffer.data(idx_dip_id + 29);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_xxxxy_xy,  \
                             tr_x_xxxxy_yy,  \
                             tr_x_xxxxyz_xx, \
                             tr_x_xxxxyz_xy, \
                             tr_x_xxxxyz_xz, \
                             tr_x_xxxxyz_yy, \
                             tr_x_xxxxyz_yz, \
                             tr_x_xxxxyz_zz, \
                             tr_x_xxxxz_xx,  \
                             tr_x_xxxxz_xz,  \
                             tr_x_xxxxz_yz,  \
                             tr_x_xxxxz_z,   \
                             tr_x_xxxxz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyz_xx[i] = tr_x_xxxxz_xx[i] * pa_y[i];

        tr_x_xxxxyz_xy[i] = tr_x_xxxxy_xy[i] * pa_z[i];

        tr_x_xxxxyz_xz[i] = tr_x_xxxxz_xz[i] * pa_y[i];

        tr_x_xxxxyz_yy[i] = tr_x_xxxxy_yy[i] * pa_z[i];

        tr_x_xxxxyz_yz[i] = tr_x_xxxxz_z[i] * fe_0 + tr_x_xxxxz_yz[i] * pa_y[i];

        tr_x_xxxxyz_zz[i] = tr_x_xxxxz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : ID

    auto tr_x_xxxxzz_xx = pbuffer.data(idx_dip_id + 30);

    auto tr_x_xxxxzz_xy = pbuffer.data(idx_dip_id + 31);

    auto tr_x_xxxxzz_xz = pbuffer.data(idx_dip_id + 32);

    auto tr_x_xxxxzz_yy = pbuffer.data(idx_dip_id + 33);

    auto tr_x_xxxxzz_yz = pbuffer.data(idx_dip_id + 34);

    auto tr_x_xxxxzz_zz = pbuffer.data(idx_dip_id + 35);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xxxx_xx,   \
                             tr_x_xxxx_xy,   \
                             tr_x_xxxx_xz,   \
                             tr_x_xxxx_yy,   \
                             tr_x_xxxxz_x,   \
                             tr_x_xxxxz_xx,  \
                             tr_x_xxxxz_xy,  \
                             tr_x_xxxxz_xz,  \
                             tr_x_xxxxz_yy,  \
                             tr_x_xxxxzz_xx, \
                             tr_x_xxxxzz_xy, \
                             tr_x_xxxxzz_xz, \
                             tr_x_xxxxzz_yy, \
                             tr_x_xxxxzz_yz, \
                             tr_x_xxxxzz_zz, \
                             tr_x_xxxzz_yz,  \
                             tr_x_xxxzz_zz,  \
                             tr_x_xxzz_yz,   \
                             tr_x_xxzz_zz,   \
                             ts_xxxzz_yz,    \
                             ts_xxxzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxzz_xx[i] = tr_x_xxxx_xx[i] * fe_0 + tr_x_xxxxz_xx[i] * pa_z[i];

        tr_x_xxxxzz_xy[i] = tr_x_xxxx_xy[i] * fe_0 + tr_x_xxxxz_xy[i] * pa_z[i];

        tr_x_xxxxzz_xz[i] = tr_x_xxxx_xz[i] * fe_0 + tr_x_xxxxz_x[i] * fe_0 + tr_x_xxxxz_xz[i] * pa_z[i];

        tr_x_xxxxzz_yy[i] = tr_x_xxxx_yy[i] * fe_0 + tr_x_xxxxz_yy[i] * pa_z[i];

        tr_x_xxxxzz_yz[i] = 3.0 * tr_x_xxzz_yz[i] * fe_0 + ts_xxxzz_yz[i] * fe_0 + tr_x_xxxzz_yz[i] * pa_x[i];

        tr_x_xxxxzz_zz[i] = 3.0 * tr_x_xxzz_zz[i] * fe_0 + ts_xxxzz_zz[i] * fe_0 + tr_x_xxxzz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : ID

    auto tr_x_xxxyyy_xx = pbuffer.data(idx_dip_id + 36);

    auto tr_x_xxxyyy_xy = pbuffer.data(idx_dip_id + 37);

    auto tr_x_xxxyyy_xz = pbuffer.data(idx_dip_id + 38);

    auto tr_x_xxxyyy_yy = pbuffer.data(idx_dip_id + 39);

    auto tr_x_xxxyyy_yz = pbuffer.data(idx_dip_id + 40);

    auto tr_x_xxxyyy_zz = pbuffer.data(idx_dip_id + 41);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xxxy_xx,   \
                             tr_x_xxxy_xy,   \
                             tr_x_xxxy_xz,   \
                             tr_x_xxxy_zz,   \
                             tr_x_xxxyy_x,   \
                             tr_x_xxxyy_xx,  \
                             tr_x_xxxyy_xy,  \
                             tr_x_xxxyy_xz,  \
                             tr_x_xxxyy_zz,  \
                             tr_x_xxxyyy_xx, \
                             tr_x_xxxyyy_xy, \
                             tr_x_xxxyyy_xz, \
                             tr_x_xxxyyy_yy, \
                             tr_x_xxxyyy_yz, \
                             tr_x_xxxyyy_zz, \
                             tr_x_xxyyy_yy,  \
                             tr_x_xxyyy_yz,  \
                             tr_x_xyyy_yy,   \
                             tr_x_xyyy_yz,   \
                             ts_xxyyy_yy,    \
                             ts_xxyyy_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyy_xx[i] = 2.0 * tr_x_xxxy_xx[i] * fe_0 + tr_x_xxxyy_xx[i] * pa_y[i];

        tr_x_xxxyyy_xy[i] = 2.0 * tr_x_xxxy_xy[i] * fe_0 + tr_x_xxxyy_x[i] * fe_0 + tr_x_xxxyy_xy[i] * pa_y[i];

        tr_x_xxxyyy_xz[i] = 2.0 * tr_x_xxxy_xz[i] * fe_0 + tr_x_xxxyy_xz[i] * pa_y[i];

        tr_x_xxxyyy_yy[i] = 2.0 * tr_x_xyyy_yy[i] * fe_0 + ts_xxyyy_yy[i] * fe_0 + tr_x_xxyyy_yy[i] * pa_x[i];

        tr_x_xxxyyy_yz[i] = 2.0 * tr_x_xyyy_yz[i] * fe_0 + ts_xxyyy_yz[i] * fe_0 + tr_x_xxyyy_yz[i] * pa_x[i];

        tr_x_xxxyyy_zz[i] = 2.0 * tr_x_xxxy_zz[i] * fe_0 + tr_x_xxxyy_zz[i] * pa_y[i];
    }

    // Set up 42-48 components of targeted buffer : ID

    auto tr_x_xxxyyz_xx = pbuffer.data(idx_dip_id + 42);

    auto tr_x_xxxyyz_xy = pbuffer.data(idx_dip_id + 43);

    auto tr_x_xxxyyz_xz = pbuffer.data(idx_dip_id + 44);

    auto tr_x_xxxyyz_yy = pbuffer.data(idx_dip_id + 45);

    auto tr_x_xxxyyz_yz = pbuffer.data(idx_dip_id + 46);

    auto tr_x_xxxyyz_zz = pbuffer.data(idx_dip_id + 47);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_xxxyy_xx,  \
                             tr_x_xxxyy_xy,  \
                             tr_x_xxxyy_y,   \
                             tr_x_xxxyy_yy,  \
                             tr_x_xxxyy_yz,  \
                             tr_x_xxxyyz_xx, \
                             tr_x_xxxyyz_xy, \
                             tr_x_xxxyyz_xz, \
                             tr_x_xxxyyz_yy, \
                             tr_x_xxxyyz_yz, \
                             tr_x_xxxyyz_zz, \
                             tr_x_xxxyz_xz,  \
                             tr_x_xxxyz_zz,  \
                             tr_x_xxxz_xz,   \
                             tr_x_xxxz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyz_xx[i] = tr_x_xxxyy_xx[i] * pa_z[i];

        tr_x_xxxyyz_xy[i] = tr_x_xxxyy_xy[i] * pa_z[i];

        tr_x_xxxyyz_xz[i] = tr_x_xxxz_xz[i] * fe_0 + tr_x_xxxyz_xz[i] * pa_y[i];

        tr_x_xxxyyz_yy[i] = tr_x_xxxyy_yy[i] * pa_z[i];

        tr_x_xxxyyz_yz[i] = tr_x_xxxyy_y[i] * fe_0 + tr_x_xxxyy_yz[i] * pa_z[i];

        tr_x_xxxyyz_zz[i] = tr_x_xxxz_zz[i] * fe_0 + tr_x_xxxyz_zz[i] * pa_y[i];
    }

    // Set up 48-54 components of targeted buffer : ID

    auto tr_x_xxxyzz_xx = pbuffer.data(idx_dip_id + 48);

    auto tr_x_xxxyzz_xy = pbuffer.data(idx_dip_id + 49);

    auto tr_x_xxxyzz_xz = pbuffer.data(idx_dip_id + 50);

    auto tr_x_xxxyzz_yy = pbuffer.data(idx_dip_id + 51);

    auto tr_x_xxxyzz_yz = pbuffer.data(idx_dip_id + 52);

    auto tr_x_xxxyzz_zz = pbuffer.data(idx_dip_id + 53);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xxxyzz_xx, \
                             tr_x_xxxyzz_xy, \
                             tr_x_xxxyzz_xz, \
                             tr_x_xxxyzz_yy, \
                             tr_x_xxxyzz_yz, \
                             tr_x_xxxyzz_zz, \
                             tr_x_xxxzz_x,   \
                             tr_x_xxxzz_xx,  \
                             tr_x_xxxzz_xy,  \
                             tr_x_xxxzz_xz,  \
                             tr_x_xxxzz_y,   \
                             tr_x_xxxzz_yy,  \
                             tr_x_xxxzz_yz,  \
                             tr_x_xxxzz_z,   \
                             tr_x_xxxzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyzz_xx[i] = tr_x_xxxzz_xx[i] * pa_y[i];

        tr_x_xxxyzz_xy[i] = tr_x_xxxzz_x[i] * fe_0 + tr_x_xxxzz_xy[i] * pa_y[i];

        tr_x_xxxyzz_xz[i] = tr_x_xxxzz_xz[i] * pa_y[i];

        tr_x_xxxyzz_yy[i] = 2.0 * tr_x_xxxzz_y[i] * fe_0 + tr_x_xxxzz_yy[i] * pa_y[i];

        tr_x_xxxyzz_yz[i] = tr_x_xxxzz_z[i] * fe_0 + tr_x_xxxzz_yz[i] * pa_y[i];

        tr_x_xxxyzz_zz[i] = tr_x_xxxzz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : ID

    auto tr_x_xxxzzz_xx = pbuffer.data(idx_dip_id + 54);

    auto tr_x_xxxzzz_xy = pbuffer.data(idx_dip_id + 55);

    auto tr_x_xxxzzz_xz = pbuffer.data(idx_dip_id + 56);

    auto tr_x_xxxzzz_yy = pbuffer.data(idx_dip_id + 57);

    auto tr_x_xxxzzz_yz = pbuffer.data(idx_dip_id + 58);

    auto tr_x_xxxzzz_zz = pbuffer.data(idx_dip_id + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xxxz_xx,   \
                             tr_x_xxxz_xy,   \
                             tr_x_xxxz_xz,   \
                             tr_x_xxxz_yy,   \
                             tr_x_xxxzz_x,   \
                             tr_x_xxxzz_xx,  \
                             tr_x_xxxzz_xy,  \
                             tr_x_xxxzz_xz,  \
                             tr_x_xxxzz_yy,  \
                             tr_x_xxxzzz_xx, \
                             tr_x_xxxzzz_xy, \
                             tr_x_xxxzzz_xz, \
                             tr_x_xxxzzz_yy, \
                             tr_x_xxxzzz_yz, \
                             tr_x_xxxzzz_zz, \
                             tr_x_xxzzz_yz,  \
                             tr_x_xxzzz_zz,  \
                             tr_x_xzzz_yz,   \
                             tr_x_xzzz_zz,   \
                             ts_xxzzz_yz,    \
                             ts_xxzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzzz_xx[i] = 2.0 * tr_x_xxxz_xx[i] * fe_0 + tr_x_xxxzz_xx[i] * pa_z[i];

        tr_x_xxxzzz_xy[i] = 2.0 * tr_x_xxxz_xy[i] * fe_0 + tr_x_xxxzz_xy[i] * pa_z[i];

        tr_x_xxxzzz_xz[i] = 2.0 * tr_x_xxxz_xz[i] * fe_0 + tr_x_xxxzz_x[i] * fe_0 + tr_x_xxxzz_xz[i] * pa_z[i];

        tr_x_xxxzzz_yy[i] = 2.0 * tr_x_xxxz_yy[i] * fe_0 + tr_x_xxxzz_yy[i] * pa_z[i];

        tr_x_xxxzzz_yz[i] = 2.0 * tr_x_xzzz_yz[i] * fe_0 + ts_xxzzz_yz[i] * fe_0 + tr_x_xxzzz_yz[i] * pa_x[i];

        tr_x_xxxzzz_zz[i] = 2.0 * tr_x_xzzz_zz[i] * fe_0 + ts_xxzzz_zz[i] * fe_0 + tr_x_xxzzz_zz[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : ID

    auto tr_x_xxyyyy_xx = pbuffer.data(idx_dip_id + 60);

    auto tr_x_xxyyyy_xy = pbuffer.data(idx_dip_id + 61);

    auto tr_x_xxyyyy_xz = pbuffer.data(idx_dip_id + 62);

    auto tr_x_xxyyyy_yy = pbuffer.data(idx_dip_id + 63);

    auto tr_x_xxyyyy_yz = pbuffer.data(idx_dip_id + 64);

    auto tr_x_xxyyyy_zz = pbuffer.data(idx_dip_id + 65);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xxyy_xx,   \
                             tr_x_xxyy_xy,   \
                             tr_x_xxyy_xz,   \
                             tr_x_xxyy_zz,   \
                             tr_x_xxyyy_x,   \
                             tr_x_xxyyy_xx,  \
                             tr_x_xxyyy_xy,  \
                             tr_x_xxyyy_xz,  \
                             tr_x_xxyyy_zz,  \
                             tr_x_xxyyyy_xx, \
                             tr_x_xxyyyy_xy, \
                             tr_x_xxyyyy_xz, \
                             tr_x_xxyyyy_yy, \
                             tr_x_xxyyyy_yz, \
                             tr_x_xxyyyy_zz, \
                             tr_x_xyyyy_yy,  \
                             tr_x_xyyyy_yz,  \
                             tr_x_yyyy_yy,   \
                             tr_x_yyyy_yz,   \
                             ts_xyyyy_yy,    \
                             ts_xyyyy_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyy_xx[i] = 3.0 * tr_x_xxyy_xx[i] * fe_0 + tr_x_xxyyy_xx[i] * pa_y[i];

        tr_x_xxyyyy_xy[i] = 3.0 * tr_x_xxyy_xy[i] * fe_0 + tr_x_xxyyy_x[i] * fe_0 + tr_x_xxyyy_xy[i] * pa_y[i];

        tr_x_xxyyyy_xz[i] = 3.0 * tr_x_xxyy_xz[i] * fe_0 + tr_x_xxyyy_xz[i] * pa_y[i];

        tr_x_xxyyyy_yy[i] = tr_x_yyyy_yy[i] * fe_0 + ts_xyyyy_yy[i] * fe_0 + tr_x_xyyyy_yy[i] * pa_x[i];

        tr_x_xxyyyy_yz[i] = tr_x_yyyy_yz[i] * fe_0 + ts_xyyyy_yz[i] * fe_0 + tr_x_xyyyy_yz[i] * pa_x[i];

        tr_x_xxyyyy_zz[i] = 3.0 * tr_x_xxyy_zz[i] * fe_0 + tr_x_xxyyy_zz[i] * pa_y[i];
    }

    // Set up 66-72 components of targeted buffer : ID

    auto tr_x_xxyyyz_xx = pbuffer.data(idx_dip_id + 66);

    auto tr_x_xxyyyz_xy = pbuffer.data(idx_dip_id + 67);

    auto tr_x_xxyyyz_xz = pbuffer.data(idx_dip_id + 68);

    auto tr_x_xxyyyz_yy = pbuffer.data(idx_dip_id + 69);

    auto tr_x_xxyyyz_yz = pbuffer.data(idx_dip_id + 70);

    auto tr_x_xxyyyz_zz = pbuffer.data(idx_dip_id + 71);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_xxyyy_xx,  \
                             tr_x_xxyyy_xy,  \
                             tr_x_xxyyy_y,   \
                             tr_x_xxyyy_yy,  \
                             tr_x_xxyyy_yz,  \
                             tr_x_xxyyyz_xx, \
                             tr_x_xxyyyz_xy, \
                             tr_x_xxyyyz_xz, \
                             tr_x_xxyyyz_yy, \
                             tr_x_xxyyyz_yz, \
                             tr_x_xxyyyz_zz, \
                             tr_x_xxyyz_xz,  \
                             tr_x_xxyyz_zz,  \
                             tr_x_xxyz_xz,   \
                             tr_x_xxyz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyz_xx[i] = tr_x_xxyyy_xx[i] * pa_z[i];

        tr_x_xxyyyz_xy[i] = tr_x_xxyyy_xy[i] * pa_z[i];

        tr_x_xxyyyz_xz[i] = 2.0 * tr_x_xxyz_xz[i] * fe_0 + tr_x_xxyyz_xz[i] * pa_y[i];

        tr_x_xxyyyz_yy[i] = tr_x_xxyyy_yy[i] * pa_z[i];

        tr_x_xxyyyz_yz[i] = tr_x_xxyyy_y[i] * fe_0 + tr_x_xxyyy_yz[i] * pa_z[i];

        tr_x_xxyyyz_zz[i] = 2.0 * tr_x_xxyz_zz[i] * fe_0 + tr_x_xxyyz_zz[i] * pa_y[i];
    }

    // Set up 72-78 components of targeted buffer : ID

    auto tr_x_xxyyzz_xx = pbuffer.data(idx_dip_id + 72);

    auto tr_x_xxyyzz_xy = pbuffer.data(idx_dip_id + 73);

    auto tr_x_xxyyzz_xz = pbuffer.data(idx_dip_id + 74);

    auto tr_x_xxyyzz_yy = pbuffer.data(idx_dip_id + 75);

    auto tr_x_xxyyzz_yz = pbuffer.data(idx_dip_id + 76);

    auto tr_x_xxyyzz_zz = pbuffer.data(idx_dip_id + 77);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xxyy_xy,   \
                             tr_x_xxyy_yy,   \
                             tr_x_xxyyz_xy,  \
                             tr_x_xxyyz_yy,  \
                             tr_x_xxyyzz_xx, \
                             tr_x_xxyyzz_xy, \
                             tr_x_xxyyzz_xz, \
                             tr_x_xxyyzz_yy, \
                             tr_x_xxyyzz_yz, \
                             tr_x_xxyyzz_zz, \
                             tr_x_xxyzz_xx,  \
                             tr_x_xxyzz_xz,  \
                             tr_x_xxyzz_zz,  \
                             tr_x_xxzz_xx,   \
                             tr_x_xxzz_xz,   \
                             tr_x_xxzz_zz,   \
                             tr_x_xyyzz_yz,  \
                             tr_x_yyzz_yz,   \
                             ts_xyyzz_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyzz_xx[i] = tr_x_xxzz_xx[i] * fe_0 + tr_x_xxyzz_xx[i] * pa_y[i];

        tr_x_xxyyzz_xy[i] = tr_x_xxyy_xy[i] * fe_0 + tr_x_xxyyz_xy[i] * pa_z[i];

        tr_x_xxyyzz_xz[i] = tr_x_xxzz_xz[i] * fe_0 + tr_x_xxyzz_xz[i] * pa_y[i];

        tr_x_xxyyzz_yy[i] = tr_x_xxyy_yy[i] * fe_0 + tr_x_xxyyz_yy[i] * pa_z[i];

        tr_x_xxyyzz_yz[i] = tr_x_yyzz_yz[i] * fe_0 + ts_xyyzz_yz[i] * fe_0 + tr_x_xyyzz_yz[i] * pa_x[i];

        tr_x_xxyyzz_zz[i] = tr_x_xxzz_zz[i] * fe_0 + tr_x_xxyzz_zz[i] * pa_y[i];
    }

    // Set up 78-84 components of targeted buffer : ID

    auto tr_x_xxyzzz_xx = pbuffer.data(idx_dip_id + 78);

    auto tr_x_xxyzzz_xy = pbuffer.data(idx_dip_id + 79);

    auto tr_x_xxyzzz_xz = pbuffer.data(idx_dip_id + 80);

    auto tr_x_xxyzzz_yy = pbuffer.data(idx_dip_id + 81);

    auto tr_x_xxyzzz_yz = pbuffer.data(idx_dip_id + 82);

    auto tr_x_xxyzzz_zz = pbuffer.data(idx_dip_id + 83);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xxyzzz_xx, \
                             tr_x_xxyzzz_xy, \
                             tr_x_xxyzzz_xz, \
                             tr_x_xxyzzz_yy, \
                             tr_x_xxyzzz_yz, \
                             tr_x_xxyzzz_zz, \
                             tr_x_xxzzz_x,   \
                             tr_x_xxzzz_xx,  \
                             tr_x_xxzzz_xy,  \
                             tr_x_xxzzz_xz,  \
                             tr_x_xxzzz_y,   \
                             tr_x_xxzzz_yy,  \
                             tr_x_xxzzz_yz,  \
                             tr_x_xxzzz_z,   \
                             tr_x_xxzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzzz_xx[i] = tr_x_xxzzz_xx[i] * pa_y[i];

        tr_x_xxyzzz_xy[i] = tr_x_xxzzz_x[i] * fe_0 + tr_x_xxzzz_xy[i] * pa_y[i];

        tr_x_xxyzzz_xz[i] = tr_x_xxzzz_xz[i] * pa_y[i];

        tr_x_xxyzzz_yy[i] = 2.0 * tr_x_xxzzz_y[i] * fe_0 + tr_x_xxzzz_yy[i] * pa_y[i];

        tr_x_xxyzzz_yz[i] = tr_x_xxzzz_z[i] * fe_0 + tr_x_xxzzz_yz[i] * pa_y[i];

        tr_x_xxyzzz_zz[i] = tr_x_xxzzz_zz[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : ID

    auto tr_x_xxzzzz_xx = pbuffer.data(idx_dip_id + 84);

    auto tr_x_xxzzzz_xy = pbuffer.data(idx_dip_id + 85);

    auto tr_x_xxzzzz_xz = pbuffer.data(idx_dip_id + 86);

    auto tr_x_xxzzzz_yy = pbuffer.data(idx_dip_id + 87);

    auto tr_x_xxzzzz_yz = pbuffer.data(idx_dip_id + 88);

    auto tr_x_xxzzzz_zz = pbuffer.data(idx_dip_id + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xxzz_xx,   \
                             tr_x_xxzz_xy,   \
                             tr_x_xxzz_xz,   \
                             tr_x_xxzz_yy,   \
                             tr_x_xxzzz_x,   \
                             tr_x_xxzzz_xx,  \
                             tr_x_xxzzz_xy,  \
                             tr_x_xxzzz_xz,  \
                             tr_x_xxzzz_yy,  \
                             tr_x_xxzzzz_xx, \
                             tr_x_xxzzzz_xy, \
                             tr_x_xxzzzz_xz, \
                             tr_x_xxzzzz_yy, \
                             tr_x_xxzzzz_yz, \
                             tr_x_xxzzzz_zz, \
                             tr_x_xzzzz_yz,  \
                             tr_x_xzzzz_zz,  \
                             tr_x_zzzz_yz,   \
                             tr_x_zzzz_zz,   \
                             ts_xzzzz_yz,    \
                             ts_xzzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzzz_xx[i] = 3.0 * tr_x_xxzz_xx[i] * fe_0 + tr_x_xxzzz_xx[i] * pa_z[i];

        tr_x_xxzzzz_xy[i] = 3.0 * tr_x_xxzz_xy[i] * fe_0 + tr_x_xxzzz_xy[i] * pa_z[i];

        tr_x_xxzzzz_xz[i] = 3.0 * tr_x_xxzz_xz[i] * fe_0 + tr_x_xxzzz_x[i] * fe_0 + tr_x_xxzzz_xz[i] * pa_z[i];

        tr_x_xxzzzz_yy[i] = 3.0 * tr_x_xxzz_yy[i] * fe_0 + tr_x_xxzzz_yy[i] * pa_z[i];

        tr_x_xxzzzz_yz[i] = tr_x_zzzz_yz[i] * fe_0 + ts_xzzzz_yz[i] * fe_0 + tr_x_xzzzz_yz[i] * pa_x[i];

        tr_x_xxzzzz_zz[i] = tr_x_zzzz_zz[i] * fe_0 + ts_xzzzz_zz[i] * fe_0 + tr_x_xzzzz_zz[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : ID

    auto tr_x_xyyyyy_xx = pbuffer.data(idx_dip_id + 90);

    auto tr_x_xyyyyy_xy = pbuffer.data(idx_dip_id + 91);

    auto tr_x_xyyyyy_xz = pbuffer.data(idx_dip_id + 92);

    auto tr_x_xyyyyy_yy = pbuffer.data(idx_dip_id + 93);

    auto tr_x_xyyyyy_yz = pbuffer.data(idx_dip_id + 94);

    auto tr_x_xyyyyy_zz = pbuffer.data(idx_dip_id + 95);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xyyy_xx,   \
                             tr_x_xyyy_xz,   \
                             tr_x_xyyyy_xx,  \
                             tr_x_xyyyy_xz,  \
                             tr_x_xyyyyy_xx, \
                             tr_x_xyyyyy_xy, \
                             tr_x_xyyyyy_xz, \
                             tr_x_xyyyyy_yy, \
                             tr_x_xyyyyy_yz, \
                             tr_x_xyyyyy_zz, \
                             tr_x_yyyyy_xy,  \
                             tr_x_yyyyy_y,   \
                             tr_x_yyyyy_yy,  \
                             tr_x_yyyyy_yz,  \
                             tr_x_yyyyy_zz,  \
                             ts_yyyyy_xy,    \
                             ts_yyyyy_yy,    \
                             ts_yyyyy_yz,    \
                             ts_yyyyy_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyy_xx[i] = 4.0 * tr_x_xyyy_xx[i] * fe_0 + tr_x_xyyyy_xx[i] * pa_y[i];

        tr_x_xyyyyy_xy[i] = tr_x_yyyyy_y[i] * fe_0 + ts_yyyyy_xy[i] * fe_0 + tr_x_yyyyy_xy[i] * pa_x[i];

        tr_x_xyyyyy_xz[i] = 4.0 * tr_x_xyyy_xz[i] * fe_0 + tr_x_xyyyy_xz[i] * pa_y[i];

        tr_x_xyyyyy_yy[i] = ts_yyyyy_yy[i] * fe_0 + tr_x_yyyyy_yy[i] * pa_x[i];

        tr_x_xyyyyy_yz[i] = ts_yyyyy_yz[i] * fe_0 + tr_x_yyyyy_yz[i] * pa_x[i];

        tr_x_xyyyyy_zz[i] = ts_yyyyy_zz[i] * fe_0 + tr_x_yyyyy_zz[i] * pa_x[i];
    }

    // Set up 96-102 components of targeted buffer : ID

    auto tr_x_xyyyyz_xx = pbuffer.data(idx_dip_id + 96);

    auto tr_x_xyyyyz_xy = pbuffer.data(idx_dip_id + 97);

    auto tr_x_xyyyyz_xz = pbuffer.data(idx_dip_id + 98);

    auto tr_x_xyyyyz_yy = pbuffer.data(idx_dip_id + 99);

    auto tr_x_xyyyyz_yz = pbuffer.data(idx_dip_id + 100);

    auto tr_x_xyyyyz_zz = pbuffer.data(idx_dip_id + 101);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xyyyy_xx,  \
                             tr_x_xyyyy_xy,  \
                             tr_x_xyyyy_yy,  \
                             tr_x_xyyyyz_xx, \
                             tr_x_xyyyyz_xy, \
                             tr_x_xyyyyz_xz, \
                             tr_x_xyyyyz_yy, \
                             tr_x_xyyyyz_yz, \
                             tr_x_xyyyyz_zz, \
                             tr_x_xyyyz_xz,  \
                             tr_x_xyyz_xz,   \
                             tr_x_yyyyz_yz,  \
                             tr_x_yyyyz_zz,  \
                             ts_yyyyz_yz,    \
                             ts_yyyyz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyz_xx[i] = tr_x_xyyyy_xx[i] * pa_z[i];

        tr_x_xyyyyz_xy[i] = tr_x_xyyyy_xy[i] * pa_z[i];

        tr_x_xyyyyz_xz[i] = 3.0 * tr_x_xyyz_xz[i] * fe_0 + tr_x_xyyyz_xz[i] * pa_y[i];

        tr_x_xyyyyz_yy[i] = tr_x_xyyyy_yy[i] * pa_z[i];

        tr_x_xyyyyz_yz[i] = ts_yyyyz_yz[i] * fe_0 + tr_x_yyyyz_yz[i] * pa_x[i];

        tr_x_xyyyyz_zz[i] = ts_yyyyz_zz[i] * fe_0 + tr_x_yyyyz_zz[i] * pa_x[i];
    }

    // Set up 102-108 components of targeted buffer : ID

    auto tr_x_xyyyzz_xx = pbuffer.data(idx_dip_id + 102);

    auto tr_x_xyyyzz_xy = pbuffer.data(idx_dip_id + 103);

    auto tr_x_xyyyzz_xz = pbuffer.data(idx_dip_id + 104);

    auto tr_x_xyyyzz_yy = pbuffer.data(idx_dip_id + 105);

    auto tr_x_xyyyzz_yz = pbuffer.data(idx_dip_id + 106);

    auto tr_x_xyyyzz_zz = pbuffer.data(idx_dip_id + 107);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xyyy_xy,   \
                             tr_x_xyyyz_xy,  \
                             tr_x_xyyyzz_xx, \
                             tr_x_xyyyzz_xy, \
                             tr_x_xyyyzz_xz, \
                             tr_x_xyyyzz_yy, \
                             tr_x_xyyyzz_yz, \
                             tr_x_xyyyzz_zz, \
                             tr_x_xyyzz_xx,  \
                             tr_x_xyyzz_xz,  \
                             tr_x_xyzz_xx,   \
                             tr_x_xyzz_xz,   \
                             tr_x_yyyzz_yy,  \
                             tr_x_yyyzz_yz,  \
                             tr_x_yyyzz_zz,  \
                             ts_yyyzz_yy,    \
                             ts_yyyzz_yz,    \
                             ts_yyyzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyzz_xx[i] = 2.0 * tr_x_xyzz_xx[i] * fe_0 + tr_x_xyyzz_xx[i] * pa_y[i];

        tr_x_xyyyzz_xy[i] = tr_x_xyyy_xy[i] * fe_0 + tr_x_xyyyz_xy[i] * pa_z[i];

        tr_x_xyyyzz_xz[i] = 2.0 * tr_x_xyzz_xz[i] * fe_0 + tr_x_xyyzz_xz[i] * pa_y[i];

        tr_x_xyyyzz_yy[i] = ts_yyyzz_yy[i] * fe_0 + tr_x_yyyzz_yy[i] * pa_x[i];

        tr_x_xyyyzz_yz[i] = ts_yyyzz_yz[i] * fe_0 + tr_x_yyyzz_yz[i] * pa_x[i];

        tr_x_xyyyzz_zz[i] = ts_yyyzz_zz[i] * fe_0 + tr_x_yyyzz_zz[i] * pa_x[i];
    }

    // Set up 108-114 components of targeted buffer : ID

    auto tr_x_xyyzzz_xx = pbuffer.data(idx_dip_id + 108);

    auto tr_x_xyyzzz_xy = pbuffer.data(idx_dip_id + 109);

    auto tr_x_xyyzzz_xz = pbuffer.data(idx_dip_id + 110);

    auto tr_x_xyyzzz_yy = pbuffer.data(idx_dip_id + 111);

    auto tr_x_xyyzzz_yz = pbuffer.data(idx_dip_id + 112);

    auto tr_x_xyyzzz_zz = pbuffer.data(idx_dip_id + 113);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xyyz_xy,   \
                             tr_x_xyyzz_xy,  \
                             tr_x_xyyzzz_xx, \
                             tr_x_xyyzzz_xy, \
                             tr_x_xyyzzz_xz, \
                             tr_x_xyyzzz_yy, \
                             tr_x_xyyzzz_yz, \
                             tr_x_xyyzzz_zz, \
                             tr_x_xyzzz_xx,  \
                             tr_x_xyzzz_xz,  \
                             tr_x_xzzz_xx,   \
                             tr_x_xzzz_xz,   \
                             tr_x_yyzzz_yy,  \
                             tr_x_yyzzz_yz,  \
                             tr_x_yyzzz_zz,  \
                             ts_yyzzz_yy,    \
                             ts_yyzzz_yz,    \
                             ts_yyzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzzz_xx[i] = tr_x_xzzz_xx[i] * fe_0 + tr_x_xyzzz_xx[i] * pa_y[i];

        tr_x_xyyzzz_xy[i] = 2.0 * tr_x_xyyz_xy[i] * fe_0 + tr_x_xyyzz_xy[i] * pa_z[i];

        tr_x_xyyzzz_xz[i] = tr_x_xzzz_xz[i] * fe_0 + tr_x_xyzzz_xz[i] * pa_y[i];

        tr_x_xyyzzz_yy[i] = ts_yyzzz_yy[i] * fe_0 + tr_x_yyzzz_yy[i] * pa_x[i];

        tr_x_xyyzzz_yz[i] = ts_yyzzz_yz[i] * fe_0 + tr_x_yyzzz_yz[i] * pa_x[i];

        tr_x_xyyzzz_zz[i] = ts_yyzzz_zz[i] * fe_0 + tr_x_yyzzz_zz[i] * pa_x[i];
    }

    // Set up 114-120 components of targeted buffer : ID

    auto tr_x_xyzzzz_xx = pbuffer.data(idx_dip_id + 114);

    auto tr_x_xyzzzz_xy = pbuffer.data(idx_dip_id + 115);

    auto tr_x_xyzzzz_xz = pbuffer.data(idx_dip_id + 116);

    auto tr_x_xyzzzz_yy = pbuffer.data(idx_dip_id + 117);

    auto tr_x_xyzzzz_yz = pbuffer.data(idx_dip_id + 118);

    auto tr_x_xyzzzz_zz = pbuffer.data(idx_dip_id + 119);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xyzzzz_xx, \
                             tr_x_xyzzzz_xy, \
                             tr_x_xyzzzz_xz, \
                             tr_x_xyzzzz_yy, \
                             tr_x_xyzzzz_yz, \
                             tr_x_xyzzzz_zz, \
                             tr_x_xzzzz_x,   \
                             tr_x_xzzzz_xx,  \
                             tr_x_xzzzz_xy,  \
                             tr_x_xzzzz_xz,  \
                             tr_x_xzzzz_zz,  \
                             tr_x_yzzzz_yy,  \
                             tr_x_yzzzz_yz,  \
                             ts_yzzzz_yy,    \
                             ts_yzzzz_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzzz_xx[i] = tr_x_xzzzz_xx[i] * pa_y[i];

        tr_x_xyzzzz_xy[i] = tr_x_xzzzz_x[i] * fe_0 + tr_x_xzzzz_xy[i] * pa_y[i];

        tr_x_xyzzzz_xz[i] = tr_x_xzzzz_xz[i] * pa_y[i];

        tr_x_xyzzzz_yy[i] = ts_yzzzz_yy[i] * fe_0 + tr_x_yzzzz_yy[i] * pa_x[i];

        tr_x_xyzzzz_yz[i] = ts_yzzzz_yz[i] * fe_0 + tr_x_yzzzz_yz[i] * pa_x[i];

        tr_x_xyzzzz_zz[i] = tr_x_xzzzz_zz[i] * pa_y[i];
    }

    // Set up 120-126 components of targeted buffer : ID

    auto tr_x_xzzzzz_xx = pbuffer.data(idx_dip_id + 120);

    auto tr_x_xzzzzz_xy = pbuffer.data(idx_dip_id + 121);

    auto tr_x_xzzzzz_xz = pbuffer.data(idx_dip_id + 122);

    auto tr_x_xzzzzz_yy = pbuffer.data(idx_dip_id + 123);

    auto tr_x_xzzzzz_yz = pbuffer.data(idx_dip_id + 124);

    auto tr_x_xzzzzz_zz = pbuffer.data(idx_dip_id + 125);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xzzz_xx,   \
                             tr_x_xzzz_xy,   \
                             tr_x_xzzzz_xx,  \
                             tr_x_xzzzz_xy,  \
                             tr_x_xzzzzz_xx, \
                             tr_x_xzzzzz_xy, \
                             tr_x_xzzzzz_xz, \
                             tr_x_xzzzzz_yy, \
                             tr_x_xzzzzz_yz, \
                             tr_x_xzzzzz_zz, \
                             tr_x_zzzzz_xz,  \
                             tr_x_zzzzz_yy,  \
                             tr_x_zzzzz_yz,  \
                             tr_x_zzzzz_z,   \
                             tr_x_zzzzz_zz,  \
                             ts_zzzzz_xz,    \
                             ts_zzzzz_yy,    \
                             ts_zzzzz_yz,    \
                             ts_zzzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzzz_xx[i] = 4.0 * tr_x_xzzz_xx[i] * fe_0 + tr_x_xzzzz_xx[i] * pa_z[i];

        tr_x_xzzzzz_xy[i] = 4.0 * tr_x_xzzz_xy[i] * fe_0 + tr_x_xzzzz_xy[i] * pa_z[i];

        tr_x_xzzzzz_xz[i] = tr_x_zzzzz_z[i] * fe_0 + ts_zzzzz_xz[i] * fe_0 + tr_x_zzzzz_xz[i] * pa_x[i];

        tr_x_xzzzzz_yy[i] = ts_zzzzz_yy[i] * fe_0 + tr_x_zzzzz_yy[i] * pa_x[i];

        tr_x_xzzzzz_yz[i] = ts_zzzzz_yz[i] * fe_0 + tr_x_zzzzz_yz[i] * pa_x[i];

        tr_x_xzzzzz_zz[i] = ts_zzzzz_zz[i] * fe_0 + tr_x_zzzzz_zz[i] * pa_x[i];
    }

    // Set up 126-132 components of targeted buffer : ID

    auto tr_x_yyyyyy_xx = pbuffer.data(idx_dip_id + 126);

    auto tr_x_yyyyyy_xy = pbuffer.data(idx_dip_id + 127);

    auto tr_x_yyyyyy_xz = pbuffer.data(idx_dip_id + 128);

    auto tr_x_yyyyyy_yy = pbuffer.data(idx_dip_id + 129);

    auto tr_x_yyyyyy_yz = pbuffer.data(idx_dip_id + 130);

    auto tr_x_yyyyyy_zz = pbuffer.data(idx_dip_id + 131);

#pragma omp simd aligned(pa_y,               \
                             tr_x_yyyy_xx,   \
                             tr_x_yyyy_xy,   \
                             tr_x_yyyy_xz,   \
                             tr_x_yyyy_yy,   \
                             tr_x_yyyy_yz,   \
                             tr_x_yyyy_zz,   \
                             tr_x_yyyyy_x,   \
                             tr_x_yyyyy_xx,  \
                             tr_x_yyyyy_xy,  \
                             tr_x_yyyyy_xz,  \
                             tr_x_yyyyy_y,   \
                             tr_x_yyyyy_yy,  \
                             tr_x_yyyyy_yz,  \
                             tr_x_yyyyy_z,   \
                             tr_x_yyyyy_zz,  \
                             tr_x_yyyyyy_xx, \
                             tr_x_yyyyyy_xy, \
                             tr_x_yyyyyy_xz, \
                             tr_x_yyyyyy_yy, \
                             tr_x_yyyyyy_yz, \
                             tr_x_yyyyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyy_xx[i] = 5.0 * tr_x_yyyy_xx[i] * fe_0 + tr_x_yyyyy_xx[i] * pa_y[i];

        tr_x_yyyyyy_xy[i] = 5.0 * tr_x_yyyy_xy[i] * fe_0 + tr_x_yyyyy_x[i] * fe_0 + tr_x_yyyyy_xy[i] * pa_y[i];

        tr_x_yyyyyy_xz[i] = 5.0 * tr_x_yyyy_xz[i] * fe_0 + tr_x_yyyyy_xz[i] * pa_y[i];

        tr_x_yyyyyy_yy[i] = 5.0 * tr_x_yyyy_yy[i] * fe_0 + 2.0 * tr_x_yyyyy_y[i] * fe_0 + tr_x_yyyyy_yy[i] * pa_y[i];

        tr_x_yyyyyy_yz[i] = 5.0 * tr_x_yyyy_yz[i] * fe_0 + tr_x_yyyyy_z[i] * fe_0 + tr_x_yyyyy_yz[i] * pa_y[i];

        tr_x_yyyyyy_zz[i] = 5.0 * tr_x_yyyy_zz[i] * fe_0 + tr_x_yyyyy_zz[i] * pa_y[i];
    }

    // Set up 132-138 components of targeted buffer : ID

    auto tr_x_yyyyyz_xx = pbuffer.data(idx_dip_id + 132);

    auto tr_x_yyyyyz_xy = pbuffer.data(idx_dip_id + 133);

    auto tr_x_yyyyyz_xz = pbuffer.data(idx_dip_id + 134);

    auto tr_x_yyyyyz_yy = pbuffer.data(idx_dip_id + 135);

    auto tr_x_yyyyyz_yz = pbuffer.data(idx_dip_id + 136);

    auto tr_x_yyyyyz_zz = pbuffer.data(idx_dip_id + 137);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyyyy_xx,  \
                             tr_x_yyyyy_xy,  \
                             tr_x_yyyyy_y,   \
                             tr_x_yyyyy_yy,  \
                             tr_x_yyyyy_yz,  \
                             tr_x_yyyyyz_xx, \
                             tr_x_yyyyyz_xy, \
                             tr_x_yyyyyz_xz, \
                             tr_x_yyyyyz_yy, \
                             tr_x_yyyyyz_yz, \
                             tr_x_yyyyyz_zz, \
                             tr_x_yyyyz_xz,  \
                             tr_x_yyyyz_zz,  \
                             tr_x_yyyz_xz,   \
                             tr_x_yyyz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyz_xx[i] = tr_x_yyyyy_xx[i] * pa_z[i];

        tr_x_yyyyyz_xy[i] = tr_x_yyyyy_xy[i] * pa_z[i];

        tr_x_yyyyyz_xz[i] = 4.0 * tr_x_yyyz_xz[i] * fe_0 + tr_x_yyyyz_xz[i] * pa_y[i];

        tr_x_yyyyyz_yy[i] = tr_x_yyyyy_yy[i] * pa_z[i];

        tr_x_yyyyyz_yz[i] = tr_x_yyyyy_y[i] * fe_0 + tr_x_yyyyy_yz[i] * pa_z[i];

        tr_x_yyyyyz_zz[i] = 4.0 * tr_x_yyyz_zz[i] * fe_0 + tr_x_yyyyz_zz[i] * pa_y[i];
    }

    // Set up 138-144 components of targeted buffer : ID

    auto tr_x_yyyyzz_xx = pbuffer.data(idx_dip_id + 138);

    auto tr_x_yyyyzz_xy = pbuffer.data(idx_dip_id + 139);

    auto tr_x_yyyyzz_xz = pbuffer.data(idx_dip_id + 140);

    auto tr_x_yyyyzz_yy = pbuffer.data(idx_dip_id + 141);

    auto tr_x_yyyyzz_yz = pbuffer.data(idx_dip_id + 142);

    auto tr_x_yyyyzz_zz = pbuffer.data(idx_dip_id + 143);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyyy_xy,   \
                             tr_x_yyyy_yy,   \
                             tr_x_yyyyz_xy,  \
                             tr_x_yyyyz_yy,  \
                             tr_x_yyyyzz_xx, \
                             tr_x_yyyyzz_xy, \
                             tr_x_yyyyzz_xz, \
                             tr_x_yyyyzz_yy, \
                             tr_x_yyyyzz_yz, \
                             tr_x_yyyyzz_zz, \
                             tr_x_yyyzz_xx,  \
                             tr_x_yyyzz_xz,  \
                             tr_x_yyyzz_yz,  \
                             tr_x_yyyzz_z,   \
                             tr_x_yyyzz_zz,  \
                             tr_x_yyzz_xx,   \
                             tr_x_yyzz_xz,   \
                             tr_x_yyzz_yz,   \
                             tr_x_yyzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyzz_xx[i] = 3.0 * tr_x_yyzz_xx[i] * fe_0 + tr_x_yyyzz_xx[i] * pa_y[i];

        tr_x_yyyyzz_xy[i] = tr_x_yyyy_xy[i] * fe_0 + tr_x_yyyyz_xy[i] * pa_z[i];

        tr_x_yyyyzz_xz[i] = 3.0 * tr_x_yyzz_xz[i] * fe_0 + tr_x_yyyzz_xz[i] * pa_y[i];

        tr_x_yyyyzz_yy[i] = tr_x_yyyy_yy[i] * fe_0 + tr_x_yyyyz_yy[i] * pa_z[i];

        tr_x_yyyyzz_yz[i] = 3.0 * tr_x_yyzz_yz[i] * fe_0 + tr_x_yyyzz_z[i] * fe_0 + tr_x_yyyzz_yz[i] * pa_y[i];

        tr_x_yyyyzz_zz[i] = 3.0 * tr_x_yyzz_zz[i] * fe_0 + tr_x_yyyzz_zz[i] * pa_y[i];
    }

    // Set up 144-150 components of targeted buffer : ID

    auto tr_x_yyyzzz_xx = pbuffer.data(idx_dip_id + 144);

    auto tr_x_yyyzzz_xy = pbuffer.data(idx_dip_id + 145);

    auto tr_x_yyyzzz_xz = pbuffer.data(idx_dip_id + 146);

    auto tr_x_yyyzzz_yy = pbuffer.data(idx_dip_id + 147);

    auto tr_x_yyyzzz_yz = pbuffer.data(idx_dip_id + 148);

    auto tr_x_yyyzzz_zz = pbuffer.data(idx_dip_id + 149);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyyz_xy,   \
                             tr_x_yyyz_yy,   \
                             tr_x_yyyzz_xy,  \
                             tr_x_yyyzz_yy,  \
                             tr_x_yyyzzz_xx, \
                             tr_x_yyyzzz_xy, \
                             tr_x_yyyzzz_xz, \
                             tr_x_yyyzzz_yy, \
                             tr_x_yyyzzz_yz, \
                             tr_x_yyyzzz_zz, \
                             tr_x_yyzzz_xx,  \
                             tr_x_yyzzz_xz,  \
                             tr_x_yyzzz_yz,  \
                             tr_x_yyzzz_z,   \
                             tr_x_yyzzz_zz,  \
                             tr_x_yzzz_xx,   \
                             tr_x_yzzz_xz,   \
                             tr_x_yzzz_yz,   \
                             tr_x_yzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzzz_xx[i] = 2.0 * tr_x_yzzz_xx[i] * fe_0 + tr_x_yyzzz_xx[i] * pa_y[i];

        tr_x_yyyzzz_xy[i] = 2.0 * tr_x_yyyz_xy[i] * fe_0 + tr_x_yyyzz_xy[i] * pa_z[i];

        tr_x_yyyzzz_xz[i] = 2.0 * tr_x_yzzz_xz[i] * fe_0 + tr_x_yyzzz_xz[i] * pa_y[i];

        tr_x_yyyzzz_yy[i] = 2.0 * tr_x_yyyz_yy[i] * fe_0 + tr_x_yyyzz_yy[i] * pa_z[i];

        tr_x_yyyzzz_yz[i] = 2.0 * tr_x_yzzz_yz[i] * fe_0 + tr_x_yyzzz_z[i] * fe_0 + tr_x_yyzzz_yz[i] * pa_y[i];

        tr_x_yyyzzz_zz[i] = 2.0 * tr_x_yzzz_zz[i] * fe_0 + tr_x_yyzzz_zz[i] * pa_y[i];
    }

    // Set up 150-156 components of targeted buffer : ID

    auto tr_x_yyzzzz_xx = pbuffer.data(idx_dip_id + 150);

    auto tr_x_yyzzzz_xy = pbuffer.data(idx_dip_id + 151);

    auto tr_x_yyzzzz_xz = pbuffer.data(idx_dip_id + 152);

    auto tr_x_yyzzzz_yy = pbuffer.data(idx_dip_id + 153);

    auto tr_x_yyzzzz_yz = pbuffer.data(idx_dip_id + 154);

    auto tr_x_yyzzzz_zz = pbuffer.data(idx_dip_id + 155);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyzz_xy,   \
                             tr_x_yyzz_yy,   \
                             tr_x_yyzzz_xy,  \
                             tr_x_yyzzz_yy,  \
                             tr_x_yyzzzz_xx, \
                             tr_x_yyzzzz_xy, \
                             tr_x_yyzzzz_xz, \
                             tr_x_yyzzzz_yy, \
                             tr_x_yyzzzz_yz, \
                             tr_x_yyzzzz_zz, \
                             tr_x_yzzzz_xx,  \
                             tr_x_yzzzz_xz,  \
                             tr_x_yzzzz_yz,  \
                             tr_x_yzzzz_z,   \
                             tr_x_yzzzz_zz,  \
                             tr_x_zzzz_xx,   \
                             tr_x_zzzz_xz,   \
                             tr_x_zzzz_yz,   \
                             tr_x_zzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzzz_xx[i] = tr_x_zzzz_xx[i] * fe_0 + tr_x_yzzzz_xx[i] * pa_y[i];

        tr_x_yyzzzz_xy[i] = 3.0 * tr_x_yyzz_xy[i] * fe_0 + tr_x_yyzzz_xy[i] * pa_z[i];

        tr_x_yyzzzz_xz[i] = tr_x_zzzz_xz[i] * fe_0 + tr_x_yzzzz_xz[i] * pa_y[i];

        tr_x_yyzzzz_yy[i] = 3.0 * tr_x_yyzz_yy[i] * fe_0 + tr_x_yyzzz_yy[i] * pa_z[i];

        tr_x_yyzzzz_yz[i] = tr_x_zzzz_yz[i] * fe_0 + tr_x_yzzzz_z[i] * fe_0 + tr_x_yzzzz_yz[i] * pa_y[i];

        tr_x_yyzzzz_zz[i] = tr_x_zzzz_zz[i] * fe_0 + tr_x_yzzzz_zz[i] * pa_y[i];
    }

    // Set up 156-162 components of targeted buffer : ID

    auto tr_x_yzzzzz_xx = pbuffer.data(idx_dip_id + 156);

    auto tr_x_yzzzzz_xy = pbuffer.data(idx_dip_id + 157);

    auto tr_x_yzzzzz_xz = pbuffer.data(idx_dip_id + 158);

    auto tr_x_yzzzzz_yy = pbuffer.data(idx_dip_id + 159);

    auto tr_x_yzzzzz_yz = pbuffer.data(idx_dip_id + 160);

    auto tr_x_yzzzzz_zz = pbuffer.data(idx_dip_id + 161);

#pragma omp simd aligned(pa_y,               \
                             tr_x_yzzzzz_xx, \
                             tr_x_yzzzzz_xy, \
                             tr_x_yzzzzz_xz, \
                             tr_x_yzzzzz_yy, \
                             tr_x_yzzzzz_yz, \
                             tr_x_yzzzzz_zz, \
                             tr_x_zzzzz_x,   \
                             tr_x_zzzzz_xx,  \
                             tr_x_zzzzz_xy,  \
                             tr_x_zzzzz_xz,  \
                             tr_x_zzzzz_y,   \
                             tr_x_zzzzz_yy,  \
                             tr_x_zzzzz_yz,  \
                             tr_x_zzzzz_z,   \
                             tr_x_zzzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzzz_xx[i] = tr_x_zzzzz_xx[i] * pa_y[i];

        tr_x_yzzzzz_xy[i] = tr_x_zzzzz_x[i] * fe_0 + tr_x_zzzzz_xy[i] * pa_y[i];

        tr_x_yzzzzz_xz[i] = tr_x_zzzzz_xz[i] * pa_y[i];

        tr_x_yzzzzz_yy[i] = 2.0 * tr_x_zzzzz_y[i] * fe_0 + tr_x_zzzzz_yy[i] * pa_y[i];

        tr_x_yzzzzz_yz[i] = tr_x_zzzzz_z[i] * fe_0 + tr_x_zzzzz_yz[i] * pa_y[i];

        tr_x_yzzzzz_zz[i] = tr_x_zzzzz_zz[i] * pa_y[i];
    }

    // Set up 162-168 components of targeted buffer : ID

    auto tr_x_zzzzzz_xx = pbuffer.data(idx_dip_id + 162);

    auto tr_x_zzzzzz_xy = pbuffer.data(idx_dip_id + 163);

    auto tr_x_zzzzzz_xz = pbuffer.data(idx_dip_id + 164);

    auto tr_x_zzzzzz_yy = pbuffer.data(idx_dip_id + 165);

    auto tr_x_zzzzzz_yz = pbuffer.data(idx_dip_id + 166);

    auto tr_x_zzzzzz_zz = pbuffer.data(idx_dip_id + 167);

#pragma omp simd aligned(pa_z,               \
                             tr_x_zzzz_xx,   \
                             tr_x_zzzz_xy,   \
                             tr_x_zzzz_xz,   \
                             tr_x_zzzz_yy,   \
                             tr_x_zzzz_yz,   \
                             tr_x_zzzz_zz,   \
                             tr_x_zzzzz_x,   \
                             tr_x_zzzzz_xx,  \
                             tr_x_zzzzz_xy,  \
                             tr_x_zzzzz_xz,  \
                             tr_x_zzzzz_y,   \
                             tr_x_zzzzz_yy,  \
                             tr_x_zzzzz_yz,  \
                             tr_x_zzzzz_z,   \
                             tr_x_zzzzz_zz,  \
                             tr_x_zzzzzz_xx, \
                             tr_x_zzzzzz_xy, \
                             tr_x_zzzzzz_xz, \
                             tr_x_zzzzzz_yy, \
                             tr_x_zzzzzz_yz, \
                             tr_x_zzzzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzzz_xx[i] = 5.0 * tr_x_zzzz_xx[i] * fe_0 + tr_x_zzzzz_xx[i] * pa_z[i];

        tr_x_zzzzzz_xy[i] = 5.0 * tr_x_zzzz_xy[i] * fe_0 + tr_x_zzzzz_xy[i] * pa_z[i];

        tr_x_zzzzzz_xz[i] = 5.0 * tr_x_zzzz_xz[i] * fe_0 + tr_x_zzzzz_x[i] * fe_0 + tr_x_zzzzz_xz[i] * pa_z[i];

        tr_x_zzzzzz_yy[i] = 5.0 * tr_x_zzzz_yy[i] * fe_0 + tr_x_zzzzz_yy[i] * pa_z[i];

        tr_x_zzzzzz_yz[i] = 5.0 * tr_x_zzzz_yz[i] * fe_0 + tr_x_zzzzz_y[i] * fe_0 + tr_x_zzzzz_yz[i] * pa_z[i];

        tr_x_zzzzzz_zz[i] = 5.0 * tr_x_zzzz_zz[i] * fe_0 + 2.0 * tr_x_zzzzz_z[i] * fe_0 + tr_x_zzzzz_zz[i] * pa_z[i];
    }

    // Set up 168-174 components of targeted buffer : ID

    auto tr_y_xxxxxx_xx = pbuffer.data(idx_dip_id + 168);

    auto tr_y_xxxxxx_xy = pbuffer.data(idx_dip_id + 169);

    auto tr_y_xxxxxx_xz = pbuffer.data(idx_dip_id + 170);

    auto tr_y_xxxxxx_yy = pbuffer.data(idx_dip_id + 171);

    auto tr_y_xxxxxx_yz = pbuffer.data(idx_dip_id + 172);

    auto tr_y_xxxxxx_zz = pbuffer.data(idx_dip_id + 173);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxxx_xx,   \
                             tr_y_xxxx_xy,   \
                             tr_y_xxxx_xz,   \
                             tr_y_xxxx_yy,   \
                             tr_y_xxxx_yz,   \
                             tr_y_xxxx_zz,   \
                             tr_y_xxxxx_x,   \
                             tr_y_xxxxx_xx,  \
                             tr_y_xxxxx_xy,  \
                             tr_y_xxxxx_xz,  \
                             tr_y_xxxxx_y,   \
                             tr_y_xxxxx_yy,  \
                             tr_y_xxxxx_yz,  \
                             tr_y_xxxxx_z,   \
                             tr_y_xxxxx_zz,  \
                             tr_y_xxxxxx_xx, \
                             tr_y_xxxxxx_xy, \
                             tr_y_xxxxxx_xz, \
                             tr_y_xxxxxx_yy, \
                             tr_y_xxxxxx_yz, \
                             tr_y_xxxxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxx_xx[i] = 5.0 * tr_y_xxxx_xx[i] * fe_0 + 2.0 * tr_y_xxxxx_x[i] * fe_0 + tr_y_xxxxx_xx[i] * pa_x[i];

        tr_y_xxxxxx_xy[i] = 5.0 * tr_y_xxxx_xy[i] * fe_0 + tr_y_xxxxx_y[i] * fe_0 + tr_y_xxxxx_xy[i] * pa_x[i];

        tr_y_xxxxxx_xz[i] = 5.0 * tr_y_xxxx_xz[i] * fe_0 + tr_y_xxxxx_z[i] * fe_0 + tr_y_xxxxx_xz[i] * pa_x[i];

        tr_y_xxxxxx_yy[i] = 5.0 * tr_y_xxxx_yy[i] * fe_0 + tr_y_xxxxx_yy[i] * pa_x[i];

        tr_y_xxxxxx_yz[i] = 5.0 * tr_y_xxxx_yz[i] * fe_0 + tr_y_xxxxx_yz[i] * pa_x[i];

        tr_y_xxxxxx_zz[i] = 5.0 * tr_y_xxxx_zz[i] * fe_0 + tr_y_xxxxx_zz[i] * pa_x[i];
    }

    // Set up 174-180 components of targeted buffer : ID

    auto tr_y_xxxxxy_xx = pbuffer.data(idx_dip_id + 174);

    auto tr_y_xxxxxy_xy = pbuffer.data(idx_dip_id + 175);

    auto tr_y_xxxxxy_xz = pbuffer.data(idx_dip_id + 176);

    auto tr_y_xxxxxy_yy = pbuffer.data(idx_dip_id + 177);

    auto tr_y_xxxxxy_yz = pbuffer.data(idx_dip_id + 178);

    auto tr_y_xxxxxy_zz = pbuffer.data(idx_dip_id + 179);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_y_xxxxx_xx,  \
                             tr_y_xxxxx_xz,  \
                             tr_y_xxxxxy_xx, \
                             tr_y_xxxxxy_xy, \
                             tr_y_xxxxxy_xz, \
                             tr_y_xxxxxy_yy, \
                             tr_y_xxxxxy_yz, \
                             tr_y_xxxxxy_zz, \
                             tr_y_xxxxy_xy,  \
                             tr_y_xxxxy_y,   \
                             tr_y_xxxxy_yy,  \
                             tr_y_xxxxy_yz,  \
                             tr_y_xxxxy_zz,  \
                             tr_y_xxxy_xy,   \
                             tr_y_xxxy_yy,   \
                             tr_y_xxxy_yz,   \
                             tr_y_xxxy_zz,   \
                             ts_xxxxx_xx,    \
                             ts_xxxxx_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxy_xx[i] = ts_xxxxx_xx[i] * fe_0 + tr_y_xxxxx_xx[i] * pa_y[i];

        tr_y_xxxxxy_xy[i] = 4.0 * tr_y_xxxy_xy[i] * fe_0 + tr_y_xxxxy_y[i] * fe_0 + tr_y_xxxxy_xy[i] * pa_x[i];

        tr_y_xxxxxy_xz[i] = ts_xxxxx_xz[i] * fe_0 + tr_y_xxxxx_xz[i] * pa_y[i];

        tr_y_xxxxxy_yy[i] = 4.0 * tr_y_xxxy_yy[i] * fe_0 + tr_y_xxxxy_yy[i] * pa_x[i];

        tr_y_xxxxxy_yz[i] = 4.0 * tr_y_xxxy_yz[i] * fe_0 + tr_y_xxxxy_yz[i] * pa_x[i];

        tr_y_xxxxxy_zz[i] = 4.0 * tr_y_xxxy_zz[i] * fe_0 + tr_y_xxxxy_zz[i] * pa_x[i];
    }

    // Set up 180-186 components of targeted buffer : ID

    auto tr_y_xxxxxz_xx = pbuffer.data(idx_dip_id + 180);

    auto tr_y_xxxxxz_xy = pbuffer.data(idx_dip_id + 181);

    auto tr_y_xxxxxz_xz = pbuffer.data(idx_dip_id + 182);

    auto tr_y_xxxxxz_yy = pbuffer.data(idx_dip_id + 183);

    auto tr_y_xxxxxz_yz = pbuffer.data(idx_dip_id + 184);

    auto tr_y_xxxxxz_zz = pbuffer.data(idx_dip_id + 185);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxxxx_x,   \
                             tr_y_xxxxx_xx,  \
                             tr_y_xxxxx_xy,  \
                             tr_y_xxxxx_xz,  \
                             tr_y_xxxxx_yy,  \
                             tr_y_xxxxxz_xx, \
                             tr_y_xxxxxz_xy, \
                             tr_y_xxxxxz_xz, \
                             tr_y_xxxxxz_yy, \
                             tr_y_xxxxxz_yz, \
                             tr_y_xxxxxz_zz, \
                             tr_y_xxxxz_yz,  \
                             tr_y_xxxxz_zz,  \
                             tr_y_xxxz_yz,   \
                             tr_y_xxxz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxz_xx[i] = tr_y_xxxxx_xx[i] * pa_z[i];

        tr_y_xxxxxz_xy[i] = tr_y_xxxxx_xy[i] * pa_z[i];

        tr_y_xxxxxz_xz[i] = tr_y_xxxxx_x[i] * fe_0 + tr_y_xxxxx_xz[i] * pa_z[i];

        tr_y_xxxxxz_yy[i] = tr_y_xxxxx_yy[i] * pa_z[i];

        tr_y_xxxxxz_yz[i] = 4.0 * tr_y_xxxz_yz[i] * fe_0 + tr_y_xxxxz_yz[i] * pa_x[i];

        tr_y_xxxxxz_zz[i] = 4.0 * tr_y_xxxz_zz[i] * fe_0 + tr_y_xxxxz_zz[i] * pa_x[i];
    }

    // Set up 186-192 components of targeted buffer : ID

    auto tr_y_xxxxyy_xx = pbuffer.data(idx_dip_id + 186);

    auto tr_y_xxxxyy_xy = pbuffer.data(idx_dip_id + 187);

    auto tr_y_xxxxyy_xz = pbuffer.data(idx_dip_id + 188);

    auto tr_y_xxxxyy_yy = pbuffer.data(idx_dip_id + 189);

    auto tr_y_xxxxyy_yz = pbuffer.data(idx_dip_id + 190);

    auto tr_y_xxxxyy_zz = pbuffer.data(idx_dip_id + 191);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxxxyy_xx, \
                             tr_y_xxxxyy_xy, \
                             tr_y_xxxxyy_xz, \
                             tr_y_xxxxyy_yy, \
                             tr_y_xxxxyy_yz, \
                             tr_y_xxxxyy_zz, \
                             tr_y_xxxyy_x,   \
                             tr_y_xxxyy_xx,  \
                             tr_y_xxxyy_xy,  \
                             tr_y_xxxyy_xz,  \
                             tr_y_xxxyy_y,   \
                             tr_y_xxxyy_yy,  \
                             tr_y_xxxyy_yz,  \
                             tr_y_xxxyy_z,   \
                             tr_y_xxxyy_zz,  \
                             tr_y_xxyy_xx,   \
                             tr_y_xxyy_xy,   \
                             tr_y_xxyy_xz,   \
                             tr_y_xxyy_yy,   \
                             tr_y_xxyy_yz,   \
                             tr_y_xxyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyy_xx[i] = 3.0 * tr_y_xxyy_xx[i] * fe_0 + 2.0 * tr_y_xxxyy_x[i] * fe_0 + tr_y_xxxyy_xx[i] * pa_x[i];

        tr_y_xxxxyy_xy[i] = 3.0 * tr_y_xxyy_xy[i] * fe_0 + tr_y_xxxyy_y[i] * fe_0 + tr_y_xxxyy_xy[i] * pa_x[i];

        tr_y_xxxxyy_xz[i] = 3.0 * tr_y_xxyy_xz[i] * fe_0 + tr_y_xxxyy_z[i] * fe_0 + tr_y_xxxyy_xz[i] * pa_x[i];

        tr_y_xxxxyy_yy[i] = 3.0 * tr_y_xxyy_yy[i] * fe_0 + tr_y_xxxyy_yy[i] * pa_x[i];

        tr_y_xxxxyy_yz[i] = 3.0 * tr_y_xxyy_yz[i] * fe_0 + tr_y_xxxyy_yz[i] * pa_x[i];

        tr_y_xxxxyy_zz[i] = 3.0 * tr_y_xxyy_zz[i] * fe_0 + tr_y_xxxyy_zz[i] * pa_x[i];
    }

    // Set up 192-198 components of targeted buffer : ID

    auto tr_y_xxxxyz_xx = pbuffer.data(idx_dip_id + 192);

    auto tr_y_xxxxyz_xy = pbuffer.data(idx_dip_id + 193);

    auto tr_y_xxxxyz_xz = pbuffer.data(idx_dip_id + 194);

    auto tr_y_xxxxyz_yy = pbuffer.data(idx_dip_id + 195);

    auto tr_y_xxxxyz_yz = pbuffer.data(idx_dip_id + 196);

    auto tr_y_xxxxyz_zz = pbuffer.data(idx_dip_id + 197);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_y_xxxxy_xx,  \
                             tr_y_xxxxy_xy,  \
                             tr_y_xxxxy_yy,  \
                             tr_y_xxxxyz_xx, \
                             tr_y_xxxxyz_xy, \
                             tr_y_xxxxyz_xz, \
                             tr_y_xxxxyz_yy, \
                             tr_y_xxxxyz_yz, \
                             tr_y_xxxxyz_zz, \
                             tr_y_xxxxz_xz,  \
                             tr_y_xxxyz_yz,  \
                             tr_y_xxxyz_zz,  \
                             tr_y_xxyz_yz,   \
                             tr_y_xxyz_zz,   \
                             ts_xxxxz_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyz_xx[i] = tr_y_xxxxy_xx[i] * pa_z[i];

        tr_y_xxxxyz_xy[i] = tr_y_xxxxy_xy[i] * pa_z[i];

        tr_y_xxxxyz_xz[i] = ts_xxxxz_xz[i] * fe_0 + tr_y_xxxxz_xz[i] * pa_y[i];

        tr_y_xxxxyz_yy[i] = tr_y_xxxxy_yy[i] * pa_z[i];

        tr_y_xxxxyz_yz[i] = 3.0 * tr_y_xxyz_yz[i] * fe_0 + tr_y_xxxyz_yz[i] * pa_x[i];

        tr_y_xxxxyz_zz[i] = 3.0 * tr_y_xxyz_zz[i] * fe_0 + tr_y_xxxyz_zz[i] * pa_x[i];
    }

    // Set up 198-204 components of targeted buffer : ID

    auto tr_y_xxxxzz_xx = pbuffer.data(idx_dip_id + 198);

    auto tr_y_xxxxzz_xy = pbuffer.data(idx_dip_id + 199);

    auto tr_y_xxxxzz_xz = pbuffer.data(idx_dip_id + 200);

    auto tr_y_xxxxzz_yy = pbuffer.data(idx_dip_id + 201);

    auto tr_y_xxxxzz_yz = pbuffer.data(idx_dip_id + 202);

    auto tr_y_xxxxzz_zz = pbuffer.data(idx_dip_id + 203);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxxx_xx,   \
                             tr_y_xxxx_xy,   \
                             tr_y_xxxxz_xx,  \
                             tr_y_xxxxz_xy,  \
                             tr_y_xxxxzz_xx, \
                             tr_y_xxxxzz_xy, \
                             tr_y_xxxxzz_xz, \
                             tr_y_xxxxzz_yy, \
                             tr_y_xxxxzz_yz, \
                             tr_y_xxxxzz_zz, \
                             tr_y_xxxzz_xz,  \
                             tr_y_xxxzz_yy,  \
                             tr_y_xxxzz_yz,  \
                             tr_y_xxxzz_z,   \
                             tr_y_xxxzz_zz,  \
                             tr_y_xxzz_xz,   \
                             tr_y_xxzz_yy,   \
                             tr_y_xxzz_yz,   \
                             tr_y_xxzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxzz_xx[i] = tr_y_xxxx_xx[i] * fe_0 + tr_y_xxxxz_xx[i] * pa_z[i];

        tr_y_xxxxzz_xy[i] = tr_y_xxxx_xy[i] * fe_0 + tr_y_xxxxz_xy[i] * pa_z[i];

        tr_y_xxxxzz_xz[i] = 3.0 * tr_y_xxzz_xz[i] * fe_0 + tr_y_xxxzz_z[i] * fe_0 + tr_y_xxxzz_xz[i] * pa_x[i];

        tr_y_xxxxzz_yy[i] = 3.0 * tr_y_xxzz_yy[i] * fe_0 + tr_y_xxxzz_yy[i] * pa_x[i];

        tr_y_xxxxzz_yz[i] = 3.0 * tr_y_xxzz_yz[i] * fe_0 + tr_y_xxxzz_yz[i] * pa_x[i];

        tr_y_xxxxzz_zz[i] = 3.0 * tr_y_xxzz_zz[i] * fe_0 + tr_y_xxxzz_zz[i] * pa_x[i];
    }

    // Set up 204-210 components of targeted buffer : ID

    auto tr_y_xxxyyy_xx = pbuffer.data(idx_dip_id + 204);

    auto tr_y_xxxyyy_xy = pbuffer.data(idx_dip_id + 205);

    auto tr_y_xxxyyy_xz = pbuffer.data(idx_dip_id + 206);

    auto tr_y_xxxyyy_yy = pbuffer.data(idx_dip_id + 207);

    auto tr_y_xxxyyy_yz = pbuffer.data(idx_dip_id + 208);

    auto tr_y_xxxyyy_zz = pbuffer.data(idx_dip_id + 209);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxxyyy_xx, \
                             tr_y_xxxyyy_xy, \
                             tr_y_xxxyyy_xz, \
                             tr_y_xxxyyy_yy, \
                             tr_y_xxxyyy_yz, \
                             tr_y_xxxyyy_zz, \
                             tr_y_xxyyy_x,   \
                             tr_y_xxyyy_xx,  \
                             tr_y_xxyyy_xy,  \
                             tr_y_xxyyy_xz,  \
                             tr_y_xxyyy_y,   \
                             tr_y_xxyyy_yy,  \
                             tr_y_xxyyy_yz,  \
                             tr_y_xxyyy_z,   \
                             tr_y_xxyyy_zz,  \
                             tr_y_xyyy_xx,   \
                             tr_y_xyyy_xy,   \
                             tr_y_xyyy_xz,   \
                             tr_y_xyyy_yy,   \
                             tr_y_xyyy_yz,   \
                             tr_y_xyyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyy_xx[i] = 2.0 * tr_y_xyyy_xx[i] * fe_0 + 2.0 * tr_y_xxyyy_x[i] * fe_0 + tr_y_xxyyy_xx[i] * pa_x[i];

        tr_y_xxxyyy_xy[i] = 2.0 * tr_y_xyyy_xy[i] * fe_0 + tr_y_xxyyy_y[i] * fe_0 + tr_y_xxyyy_xy[i] * pa_x[i];

        tr_y_xxxyyy_xz[i] = 2.0 * tr_y_xyyy_xz[i] * fe_0 + tr_y_xxyyy_z[i] * fe_0 + tr_y_xxyyy_xz[i] * pa_x[i];

        tr_y_xxxyyy_yy[i] = 2.0 * tr_y_xyyy_yy[i] * fe_0 + tr_y_xxyyy_yy[i] * pa_x[i];

        tr_y_xxxyyy_yz[i] = 2.0 * tr_y_xyyy_yz[i] * fe_0 + tr_y_xxyyy_yz[i] * pa_x[i];

        tr_y_xxxyyy_zz[i] = 2.0 * tr_y_xyyy_zz[i] * fe_0 + tr_y_xxyyy_zz[i] * pa_x[i];
    }

    // Set up 210-216 components of targeted buffer : ID

    auto tr_y_xxxyyz_xx = pbuffer.data(idx_dip_id + 210);

    auto tr_y_xxxyyz_xy = pbuffer.data(idx_dip_id + 211);

    auto tr_y_xxxyyz_xz = pbuffer.data(idx_dip_id + 212);

    auto tr_y_xxxyyz_yy = pbuffer.data(idx_dip_id + 213);

    auto tr_y_xxxyyz_yz = pbuffer.data(idx_dip_id + 214);

    auto tr_y_xxxyyz_zz = pbuffer.data(idx_dip_id + 215);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxxyy_x,   \
                             tr_y_xxxyy_xx,  \
                             tr_y_xxxyy_xy,  \
                             tr_y_xxxyy_xz,  \
                             tr_y_xxxyy_yy,  \
                             tr_y_xxxyyz_xx, \
                             tr_y_xxxyyz_xy, \
                             tr_y_xxxyyz_xz, \
                             tr_y_xxxyyz_yy, \
                             tr_y_xxxyyz_yz, \
                             tr_y_xxxyyz_zz, \
                             tr_y_xxyyz_yz,  \
                             tr_y_xxyyz_zz,  \
                             tr_y_xyyz_yz,   \
                             tr_y_xyyz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyz_xx[i] = tr_y_xxxyy_xx[i] * pa_z[i];

        tr_y_xxxyyz_xy[i] = tr_y_xxxyy_xy[i] * pa_z[i];

        tr_y_xxxyyz_xz[i] = tr_y_xxxyy_x[i] * fe_0 + tr_y_xxxyy_xz[i] * pa_z[i];

        tr_y_xxxyyz_yy[i] = tr_y_xxxyy_yy[i] * pa_z[i];

        tr_y_xxxyyz_yz[i] = 2.0 * tr_y_xyyz_yz[i] * fe_0 + tr_y_xxyyz_yz[i] * pa_x[i];

        tr_y_xxxyyz_zz[i] = 2.0 * tr_y_xyyz_zz[i] * fe_0 + tr_y_xxyyz_zz[i] * pa_x[i];
    }

    // Set up 216-222 components of targeted buffer : ID

    auto tr_y_xxxyzz_xx = pbuffer.data(idx_dip_id + 216);

    auto tr_y_xxxyzz_xy = pbuffer.data(idx_dip_id + 217);

    auto tr_y_xxxyzz_xz = pbuffer.data(idx_dip_id + 218);

    auto tr_y_xxxyzz_yy = pbuffer.data(idx_dip_id + 219);

    auto tr_y_xxxyzz_yz = pbuffer.data(idx_dip_id + 220);

    auto tr_y_xxxyzz_zz = pbuffer.data(idx_dip_id + 221);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_y_xxxy_xy,   \
                             tr_y_xxxyz_xy,  \
                             tr_y_xxxyzz_xx, \
                             tr_y_xxxyzz_xy, \
                             tr_y_xxxyzz_xz, \
                             tr_y_xxxyzz_yy, \
                             tr_y_xxxyzz_yz, \
                             tr_y_xxxyzz_zz, \
                             tr_y_xxxzz_xx,  \
                             tr_y_xxxzz_xz,  \
                             tr_y_xxyzz_yy,  \
                             tr_y_xxyzz_yz,  \
                             tr_y_xxyzz_zz,  \
                             tr_y_xyzz_yy,   \
                             tr_y_xyzz_yz,   \
                             tr_y_xyzz_zz,   \
                             ts_xxxzz_xx,    \
                             ts_xxxzz_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyzz_xx[i] = ts_xxxzz_xx[i] * fe_0 + tr_y_xxxzz_xx[i] * pa_y[i];

        tr_y_xxxyzz_xy[i] = tr_y_xxxy_xy[i] * fe_0 + tr_y_xxxyz_xy[i] * pa_z[i];

        tr_y_xxxyzz_xz[i] = ts_xxxzz_xz[i] * fe_0 + tr_y_xxxzz_xz[i] * pa_y[i];

        tr_y_xxxyzz_yy[i] = 2.0 * tr_y_xyzz_yy[i] * fe_0 + tr_y_xxyzz_yy[i] * pa_x[i];

        tr_y_xxxyzz_yz[i] = 2.0 * tr_y_xyzz_yz[i] * fe_0 + tr_y_xxyzz_yz[i] * pa_x[i];

        tr_y_xxxyzz_zz[i] = 2.0 * tr_y_xyzz_zz[i] * fe_0 + tr_y_xxyzz_zz[i] * pa_x[i];
    }

    // Set up 222-228 components of targeted buffer : ID

    auto tr_y_xxxzzz_xx = pbuffer.data(idx_dip_id + 222);

    auto tr_y_xxxzzz_xy = pbuffer.data(idx_dip_id + 223);

    auto tr_y_xxxzzz_xz = pbuffer.data(idx_dip_id + 224);

    auto tr_y_xxxzzz_yy = pbuffer.data(idx_dip_id + 225);

    auto tr_y_xxxzzz_yz = pbuffer.data(idx_dip_id + 226);

    auto tr_y_xxxzzz_zz = pbuffer.data(idx_dip_id + 227);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxxz_xx,   \
                             tr_y_xxxz_xy,   \
                             tr_y_xxxzz_xx,  \
                             tr_y_xxxzz_xy,  \
                             tr_y_xxxzzz_xx, \
                             tr_y_xxxzzz_xy, \
                             tr_y_xxxzzz_xz, \
                             tr_y_xxxzzz_yy, \
                             tr_y_xxxzzz_yz, \
                             tr_y_xxxzzz_zz, \
                             tr_y_xxzzz_xz,  \
                             tr_y_xxzzz_yy,  \
                             tr_y_xxzzz_yz,  \
                             tr_y_xxzzz_z,   \
                             tr_y_xxzzz_zz,  \
                             tr_y_xzzz_xz,   \
                             tr_y_xzzz_yy,   \
                             tr_y_xzzz_yz,   \
                             tr_y_xzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzzz_xx[i] = 2.0 * tr_y_xxxz_xx[i] * fe_0 + tr_y_xxxzz_xx[i] * pa_z[i];

        tr_y_xxxzzz_xy[i] = 2.0 * tr_y_xxxz_xy[i] * fe_0 + tr_y_xxxzz_xy[i] * pa_z[i];

        tr_y_xxxzzz_xz[i] = 2.0 * tr_y_xzzz_xz[i] * fe_0 + tr_y_xxzzz_z[i] * fe_0 + tr_y_xxzzz_xz[i] * pa_x[i];

        tr_y_xxxzzz_yy[i] = 2.0 * tr_y_xzzz_yy[i] * fe_0 + tr_y_xxzzz_yy[i] * pa_x[i];

        tr_y_xxxzzz_yz[i] = 2.0 * tr_y_xzzz_yz[i] * fe_0 + tr_y_xxzzz_yz[i] * pa_x[i];

        tr_y_xxxzzz_zz[i] = 2.0 * tr_y_xzzz_zz[i] * fe_0 + tr_y_xxzzz_zz[i] * pa_x[i];
    }

    // Set up 228-234 components of targeted buffer : ID

    auto tr_y_xxyyyy_xx = pbuffer.data(idx_dip_id + 228);

    auto tr_y_xxyyyy_xy = pbuffer.data(idx_dip_id + 229);

    auto tr_y_xxyyyy_xz = pbuffer.data(idx_dip_id + 230);

    auto tr_y_xxyyyy_yy = pbuffer.data(idx_dip_id + 231);

    auto tr_y_xxyyyy_yz = pbuffer.data(idx_dip_id + 232);

    auto tr_y_xxyyyy_zz = pbuffer.data(idx_dip_id + 233);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxyyyy_xx, \
                             tr_y_xxyyyy_xy, \
                             tr_y_xxyyyy_xz, \
                             tr_y_xxyyyy_yy, \
                             tr_y_xxyyyy_yz, \
                             tr_y_xxyyyy_zz, \
                             tr_y_xyyyy_x,   \
                             tr_y_xyyyy_xx,  \
                             tr_y_xyyyy_xy,  \
                             tr_y_xyyyy_xz,  \
                             tr_y_xyyyy_y,   \
                             tr_y_xyyyy_yy,  \
                             tr_y_xyyyy_yz,  \
                             tr_y_xyyyy_z,   \
                             tr_y_xyyyy_zz,  \
                             tr_y_yyyy_xx,   \
                             tr_y_yyyy_xy,   \
                             tr_y_yyyy_xz,   \
                             tr_y_yyyy_yy,   \
                             tr_y_yyyy_yz,   \
                             tr_y_yyyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyy_xx[i] = tr_y_yyyy_xx[i] * fe_0 + 2.0 * tr_y_xyyyy_x[i] * fe_0 + tr_y_xyyyy_xx[i] * pa_x[i];

        tr_y_xxyyyy_xy[i] = tr_y_yyyy_xy[i] * fe_0 + tr_y_xyyyy_y[i] * fe_0 + tr_y_xyyyy_xy[i] * pa_x[i];

        tr_y_xxyyyy_xz[i] = tr_y_yyyy_xz[i] * fe_0 + tr_y_xyyyy_z[i] * fe_0 + tr_y_xyyyy_xz[i] * pa_x[i];

        tr_y_xxyyyy_yy[i] = tr_y_yyyy_yy[i] * fe_0 + tr_y_xyyyy_yy[i] * pa_x[i];

        tr_y_xxyyyy_yz[i] = tr_y_yyyy_yz[i] * fe_0 + tr_y_xyyyy_yz[i] * pa_x[i];

        tr_y_xxyyyy_zz[i] = tr_y_yyyy_zz[i] * fe_0 + tr_y_xyyyy_zz[i] * pa_x[i];
    }

    // Set up 234-240 components of targeted buffer : ID

    auto tr_y_xxyyyz_xx = pbuffer.data(idx_dip_id + 234);

    auto tr_y_xxyyyz_xy = pbuffer.data(idx_dip_id + 235);

    auto tr_y_xxyyyz_xz = pbuffer.data(idx_dip_id + 236);

    auto tr_y_xxyyyz_yy = pbuffer.data(idx_dip_id + 237);

    auto tr_y_xxyyyz_yz = pbuffer.data(idx_dip_id + 238);

    auto tr_y_xxyyyz_zz = pbuffer.data(idx_dip_id + 239);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxyyy_x,   \
                             tr_y_xxyyy_xx,  \
                             tr_y_xxyyy_xy,  \
                             tr_y_xxyyy_xz,  \
                             tr_y_xxyyy_yy,  \
                             tr_y_xxyyyz_xx, \
                             tr_y_xxyyyz_xy, \
                             tr_y_xxyyyz_xz, \
                             tr_y_xxyyyz_yy, \
                             tr_y_xxyyyz_yz, \
                             tr_y_xxyyyz_zz, \
                             tr_y_xyyyz_yz,  \
                             tr_y_xyyyz_zz,  \
                             tr_y_yyyz_yz,   \
                             tr_y_yyyz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyz_xx[i] = tr_y_xxyyy_xx[i] * pa_z[i];

        tr_y_xxyyyz_xy[i] = tr_y_xxyyy_xy[i] * pa_z[i];

        tr_y_xxyyyz_xz[i] = tr_y_xxyyy_x[i] * fe_0 + tr_y_xxyyy_xz[i] * pa_z[i];

        tr_y_xxyyyz_yy[i] = tr_y_xxyyy_yy[i] * pa_z[i];

        tr_y_xxyyyz_yz[i] = tr_y_yyyz_yz[i] * fe_0 + tr_y_xyyyz_yz[i] * pa_x[i];

        tr_y_xxyyyz_zz[i] = tr_y_yyyz_zz[i] * fe_0 + tr_y_xyyyz_zz[i] * pa_x[i];
    }

    // Set up 240-246 components of targeted buffer : ID

    auto tr_y_xxyyzz_xx = pbuffer.data(idx_dip_id + 240);

    auto tr_y_xxyyzz_xy = pbuffer.data(idx_dip_id + 241);

    auto tr_y_xxyyzz_xz = pbuffer.data(idx_dip_id + 242);

    auto tr_y_xxyyzz_yy = pbuffer.data(idx_dip_id + 243);

    auto tr_y_xxyyzz_yz = pbuffer.data(idx_dip_id + 244);

    auto tr_y_xxyyzz_zz = pbuffer.data(idx_dip_id + 245);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxyy_xx,   \
                             tr_y_xxyy_xy,   \
                             tr_y_xxyyz_xx,  \
                             tr_y_xxyyz_xy,  \
                             tr_y_xxyyzz_xx, \
                             tr_y_xxyyzz_xy, \
                             tr_y_xxyyzz_xz, \
                             tr_y_xxyyzz_yy, \
                             tr_y_xxyyzz_yz, \
                             tr_y_xxyyzz_zz, \
                             tr_y_xyyzz_xz,  \
                             tr_y_xyyzz_yy,  \
                             tr_y_xyyzz_yz,  \
                             tr_y_xyyzz_z,   \
                             tr_y_xyyzz_zz,  \
                             tr_y_yyzz_xz,   \
                             tr_y_yyzz_yy,   \
                             tr_y_yyzz_yz,   \
                             tr_y_yyzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyzz_xx[i] = tr_y_xxyy_xx[i] * fe_0 + tr_y_xxyyz_xx[i] * pa_z[i];

        tr_y_xxyyzz_xy[i] = tr_y_xxyy_xy[i] * fe_0 + tr_y_xxyyz_xy[i] * pa_z[i];

        tr_y_xxyyzz_xz[i] = tr_y_yyzz_xz[i] * fe_0 + tr_y_xyyzz_z[i] * fe_0 + tr_y_xyyzz_xz[i] * pa_x[i];

        tr_y_xxyyzz_yy[i] = tr_y_yyzz_yy[i] * fe_0 + tr_y_xyyzz_yy[i] * pa_x[i];

        tr_y_xxyyzz_yz[i] = tr_y_yyzz_yz[i] * fe_0 + tr_y_xyyzz_yz[i] * pa_x[i];

        tr_y_xxyyzz_zz[i] = tr_y_yyzz_zz[i] * fe_0 + tr_y_xyyzz_zz[i] * pa_x[i];
    }

    // Set up 246-252 components of targeted buffer : ID

    auto tr_y_xxyzzz_xx = pbuffer.data(idx_dip_id + 246);

    auto tr_y_xxyzzz_xy = pbuffer.data(idx_dip_id + 247);

    auto tr_y_xxyzzz_xz = pbuffer.data(idx_dip_id + 248);

    auto tr_y_xxyzzz_yy = pbuffer.data(idx_dip_id + 249);

    auto tr_y_xxyzzz_yz = pbuffer.data(idx_dip_id + 250);

    auto tr_y_xxyzzz_zz = pbuffer.data(idx_dip_id + 251);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_y_xxyz_xy,   \
                             tr_y_xxyzz_xy,  \
                             tr_y_xxyzzz_xx, \
                             tr_y_xxyzzz_xy, \
                             tr_y_xxyzzz_xz, \
                             tr_y_xxyzzz_yy, \
                             tr_y_xxyzzz_yz, \
                             tr_y_xxyzzz_zz, \
                             tr_y_xxzzz_xx,  \
                             tr_y_xxzzz_xz,  \
                             tr_y_xyzzz_yy,  \
                             tr_y_xyzzz_yz,  \
                             tr_y_xyzzz_zz,  \
                             tr_y_yzzz_yy,   \
                             tr_y_yzzz_yz,   \
                             tr_y_yzzz_zz,   \
                             ts_xxzzz_xx,    \
                             ts_xxzzz_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzzz_xx[i] = ts_xxzzz_xx[i] * fe_0 + tr_y_xxzzz_xx[i] * pa_y[i];

        tr_y_xxyzzz_xy[i] = 2.0 * tr_y_xxyz_xy[i] * fe_0 + tr_y_xxyzz_xy[i] * pa_z[i];

        tr_y_xxyzzz_xz[i] = ts_xxzzz_xz[i] * fe_0 + tr_y_xxzzz_xz[i] * pa_y[i];

        tr_y_xxyzzz_yy[i] = tr_y_yzzz_yy[i] * fe_0 + tr_y_xyzzz_yy[i] * pa_x[i];

        tr_y_xxyzzz_yz[i] = tr_y_yzzz_yz[i] * fe_0 + tr_y_xyzzz_yz[i] * pa_x[i];

        tr_y_xxyzzz_zz[i] = tr_y_yzzz_zz[i] * fe_0 + tr_y_xyzzz_zz[i] * pa_x[i];
    }

    // Set up 252-258 components of targeted buffer : ID

    auto tr_y_xxzzzz_xx = pbuffer.data(idx_dip_id + 252);

    auto tr_y_xxzzzz_xy = pbuffer.data(idx_dip_id + 253);

    auto tr_y_xxzzzz_xz = pbuffer.data(idx_dip_id + 254);

    auto tr_y_xxzzzz_yy = pbuffer.data(idx_dip_id + 255);

    auto tr_y_xxzzzz_yz = pbuffer.data(idx_dip_id + 256);

    auto tr_y_xxzzzz_zz = pbuffer.data(idx_dip_id + 257);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxzz_xx,   \
                             tr_y_xxzz_xy,   \
                             tr_y_xxzzz_xx,  \
                             tr_y_xxzzz_xy,  \
                             tr_y_xxzzzz_xx, \
                             tr_y_xxzzzz_xy, \
                             tr_y_xxzzzz_xz, \
                             tr_y_xxzzzz_yy, \
                             tr_y_xxzzzz_yz, \
                             tr_y_xxzzzz_zz, \
                             tr_y_xzzzz_xz,  \
                             tr_y_xzzzz_yy,  \
                             tr_y_xzzzz_yz,  \
                             tr_y_xzzzz_z,   \
                             tr_y_xzzzz_zz,  \
                             tr_y_zzzz_xz,   \
                             tr_y_zzzz_yy,   \
                             tr_y_zzzz_yz,   \
                             tr_y_zzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzzz_xx[i] = 3.0 * tr_y_xxzz_xx[i] * fe_0 + tr_y_xxzzz_xx[i] * pa_z[i];

        tr_y_xxzzzz_xy[i] = 3.0 * tr_y_xxzz_xy[i] * fe_0 + tr_y_xxzzz_xy[i] * pa_z[i];

        tr_y_xxzzzz_xz[i] = tr_y_zzzz_xz[i] * fe_0 + tr_y_xzzzz_z[i] * fe_0 + tr_y_xzzzz_xz[i] * pa_x[i];

        tr_y_xxzzzz_yy[i] = tr_y_zzzz_yy[i] * fe_0 + tr_y_xzzzz_yy[i] * pa_x[i];

        tr_y_xxzzzz_yz[i] = tr_y_zzzz_yz[i] * fe_0 + tr_y_xzzzz_yz[i] * pa_x[i];

        tr_y_xxzzzz_zz[i] = tr_y_zzzz_zz[i] * fe_0 + tr_y_xzzzz_zz[i] * pa_x[i];
    }

    // Set up 258-264 components of targeted buffer : ID

    auto tr_y_xyyyyy_xx = pbuffer.data(idx_dip_id + 258);

    auto tr_y_xyyyyy_xy = pbuffer.data(idx_dip_id + 259);

    auto tr_y_xyyyyy_xz = pbuffer.data(idx_dip_id + 260);

    auto tr_y_xyyyyy_yy = pbuffer.data(idx_dip_id + 261);

    auto tr_y_xyyyyy_yz = pbuffer.data(idx_dip_id + 262);

    auto tr_y_xyyyyy_zz = pbuffer.data(idx_dip_id + 263);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyyyyy_xx, \
                             tr_y_xyyyyy_xy, \
                             tr_y_xyyyyy_xz, \
                             tr_y_xyyyyy_yy, \
                             tr_y_xyyyyy_yz, \
                             tr_y_xyyyyy_zz, \
                             tr_y_yyyyy_x,   \
                             tr_y_yyyyy_xx,  \
                             tr_y_yyyyy_xy,  \
                             tr_y_yyyyy_xz,  \
                             tr_y_yyyyy_y,   \
                             tr_y_yyyyy_yy,  \
                             tr_y_yyyyy_yz,  \
                             tr_y_yyyyy_z,   \
                             tr_y_yyyyy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyy_xx[i] = 2.0 * tr_y_yyyyy_x[i] * fe_0 + tr_y_yyyyy_xx[i] * pa_x[i];

        tr_y_xyyyyy_xy[i] = tr_y_yyyyy_y[i] * fe_0 + tr_y_yyyyy_xy[i] * pa_x[i];

        tr_y_xyyyyy_xz[i] = tr_y_yyyyy_z[i] * fe_0 + tr_y_yyyyy_xz[i] * pa_x[i];

        tr_y_xyyyyy_yy[i] = tr_y_yyyyy_yy[i] * pa_x[i];

        tr_y_xyyyyy_yz[i] = tr_y_yyyyy_yz[i] * pa_x[i];

        tr_y_xyyyyy_zz[i] = tr_y_yyyyy_zz[i] * pa_x[i];
    }

    // Set up 264-270 components of targeted buffer : ID

    auto tr_y_xyyyyz_xx = pbuffer.data(idx_dip_id + 264);

    auto tr_y_xyyyyz_xy = pbuffer.data(idx_dip_id + 265);

    auto tr_y_xyyyyz_xz = pbuffer.data(idx_dip_id + 266);

    auto tr_y_xyyyyz_yy = pbuffer.data(idx_dip_id + 267);

    auto tr_y_xyyyyz_yz = pbuffer.data(idx_dip_id + 268);

    auto tr_y_xyyyyz_zz = pbuffer.data(idx_dip_id + 269);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xyyyy_xx,  \
                             tr_y_xyyyy_xy,  \
                             tr_y_xyyyyz_xx, \
                             tr_y_xyyyyz_xy, \
                             tr_y_xyyyyz_xz, \
                             tr_y_xyyyyz_yy, \
                             tr_y_xyyyyz_yz, \
                             tr_y_xyyyyz_zz, \
                             tr_y_yyyyz_xz,  \
                             tr_y_yyyyz_yy,  \
                             tr_y_yyyyz_yz,  \
                             tr_y_yyyyz_z,   \
                             tr_y_yyyyz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyz_xx[i] = tr_y_xyyyy_xx[i] * pa_z[i];

        tr_y_xyyyyz_xy[i] = tr_y_xyyyy_xy[i] * pa_z[i];

        tr_y_xyyyyz_xz[i] = tr_y_yyyyz_z[i] * fe_0 + tr_y_yyyyz_xz[i] * pa_x[i];

        tr_y_xyyyyz_yy[i] = tr_y_yyyyz_yy[i] * pa_x[i];

        tr_y_xyyyyz_yz[i] = tr_y_yyyyz_yz[i] * pa_x[i];

        tr_y_xyyyyz_zz[i] = tr_y_yyyyz_zz[i] * pa_x[i];
    }

    // Set up 270-276 components of targeted buffer : ID

    auto tr_y_xyyyzz_xx = pbuffer.data(idx_dip_id + 270);

    auto tr_y_xyyyzz_xy = pbuffer.data(idx_dip_id + 271);

    auto tr_y_xyyyzz_xz = pbuffer.data(idx_dip_id + 272);

    auto tr_y_xyyyzz_yy = pbuffer.data(idx_dip_id + 273);

    auto tr_y_xyyyzz_yz = pbuffer.data(idx_dip_id + 274);

    auto tr_y_xyyyzz_zz = pbuffer.data(idx_dip_id + 275);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyyyzz_xx, \
                             tr_y_xyyyzz_xy, \
                             tr_y_xyyyzz_xz, \
                             tr_y_xyyyzz_yy, \
                             tr_y_xyyyzz_yz, \
                             tr_y_xyyyzz_zz, \
                             tr_y_yyyzz_x,   \
                             tr_y_yyyzz_xx,  \
                             tr_y_yyyzz_xy,  \
                             tr_y_yyyzz_xz,  \
                             tr_y_yyyzz_y,   \
                             tr_y_yyyzz_yy,  \
                             tr_y_yyyzz_yz,  \
                             tr_y_yyyzz_z,   \
                             tr_y_yyyzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyzz_xx[i] = 2.0 * tr_y_yyyzz_x[i] * fe_0 + tr_y_yyyzz_xx[i] * pa_x[i];

        tr_y_xyyyzz_xy[i] = tr_y_yyyzz_y[i] * fe_0 + tr_y_yyyzz_xy[i] * pa_x[i];

        tr_y_xyyyzz_xz[i] = tr_y_yyyzz_z[i] * fe_0 + tr_y_yyyzz_xz[i] * pa_x[i];

        tr_y_xyyyzz_yy[i] = tr_y_yyyzz_yy[i] * pa_x[i];

        tr_y_xyyyzz_yz[i] = tr_y_yyyzz_yz[i] * pa_x[i];

        tr_y_xyyyzz_zz[i] = tr_y_yyyzz_zz[i] * pa_x[i];
    }

    // Set up 276-282 components of targeted buffer : ID

    auto tr_y_xyyzzz_xx = pbuffer.data(idx_dip_id + 276);

    auto tr_y_xyyzzz_xy = pbuffer.data(idx_dip_id + 277);

    auto tr_y_xyyzzz_xz = pbuffer.data(idx_dip_id + 278);

    auto tr_y_xyyzzz_yy = pbuffer.data(idx_dip_id + 279);

    auto tr_y_xyyzzz_yz = pbuffer.data(idx_dip_id + 280);

    auto tr_y_xyyzzz_zz = pbuffer.data(idx_dip_id + 281);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyyzzz_xx, \
                             tr_y_xyyzzz_xy, \
                             tr_y_xyyzzz_xz, \
                             tr_y_xyyzzz_yy, \
                             tr_y_xyyzzz_yz, \
                             tr_y_xyyzzz_zz, \
                             tr_y_yyzzz_x,   \
                             tr_y_yyzzz_xx,  \
                             tr_y_yyzzz_xy,  \
                             tr_y_yyzzz_xz,  \
                             tr_y_yyzzz_y,   \
                             tr_y_yyzzz_yy,  \
                             tr_y_yyzzz_yz,  \
                             tr_y_yyzzz_z,   \
                             tr_y_yyzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzzz_xx[i] = 2.0 * tr_y_yyzzz_x[i] * fe_0 + tr_y_yyzzz_xx[i] * pa_x[i];

        tr_y_xyyzzz_xy[i] = tr_y_yyzzz_y[i] * fe_0 + tr_y_yyzzz_xy[i] * pa_x[i];

        tr_y_xyyzzz_xz[i] = tr_y_yyzzz_z[i] * fe_0 + tr_y_yyzzz_xz[i] * pa_x[i];

        tr_y_xyyzzz_yy[i] = tr_y_yyzzz_yy[i] * pa_x[i];

        tr_y_xyyzzz_yz[i] = tr_y_yyzzz_yz[i] * pa_x[i];

        tr_y_xyyzzz_zz[i] = tr_y_yyzzz_zz[i] * pa_x[i];
    }

    // Set up 282-288 components of targeted buffer : ID

    auto tr_y_xyzzzz_xx = pbuffer.data(idx_dip_id + 282);

    auto tr_y_xyzzzz_xy = pbuffer.data(idx_dip_id + 283);

    auto tr_y_xyzzzz_xz = pbuffer.data(idx_dip_id + 284);

    auto tr_y_xyzzzz_yy = pbuffer.data(idx_dip_id + 285);

    auto tr_y_xyzzzz_yz = pbuffer.data(idx_dip_id + 286);

    auto tr_y_xyzzzz_zz = pbuffer.data(idx_dip_id + 287);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyzzzz_xx, \
                             tr_y_xyzzzz_xy, \
                             tr_y_xyzzzz_xz, \
                             tr_y_xyzzzz_yy, \
                             tr_y_xyzzzz_yz, \
                             tr_y_xyzzzz_zz, \
                             tr_y_yzzzz_x,   \
                             tr_y_yzzzz_xx,  \
                             tr_y_yzzzz_xy,  \
                             tr_y_yzzzz_xz,  \
                             tr_y_yzzzz_y,   \
                             tr_y_yzzzz_yy,  \
                             tr_y_yzzzz_yz,  \
                             tr_y_yzzzz_z,   \
                             tr_y_yzzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzzz_xx[i] = 2.0 * tr_y_yzzzz_x[i] * fe_0 + tr_y_yzzzz_xx[i] * pa_x[i];

        tr_y_xyzzzz_xy[i] = tr_y_yzzzz_y[i] * fe_0 + tr_y_yzzzz_xy[i] * pa_x[i];

        tr_y_xyzzzz_xz[i] = tr_y_yzzzz_z[i] * fe_0 + tr_y_yzzzz_xz[i] * pa_x[i];

        tr_y_xyzzzz_yy[i] = tr_y_yzzzz_yy[i] * pa_x[i];

        tr_y_xyzzzz_yz[i] = tr_y_yzzzz_yz[i] * pa_x[i];

        tr_y_xyzzzz_zz[i] = tr_y_yzzzz_zz[i] * pa_x[i];
    }

    // Set up 288-294 components of targeted buffer : ID

    auto tr_y_xzzzzz_xx = pbuffer.data(idx_dip_id + 288);

    auto tr_y_xzzzzz_xy = pbuffer.data(idx_dip_id + 289);

    auto tr_y_xzzzzz_xz = pbuffer.data(idx_dip_id + 290);

    auto tr_y_xzzzzz_yy = pbuffer.data(idx_dip_id + 291);

    auto tr_y_xzzzzz_yz = pbuffer.data(idx_dip_id + 292);

    auto tr_y_xzzzzz_zz = pbuffer.data(idx_dip_id + 293);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xzzzzz_xx, \
                             tr_y_xzzzzz_xy, \
                             tr_y_xzzzzz_xz, \
                             tr_y_xzzzzz_yy, \
                             tr_y_xzzzzz_yz, \
                             tr_y_xzzzzz_zz, \
                             tr_y_zzzzz_x,   \
                             tr_y_zzzzz_xx,  \
                             tr_y_zzzzz_xy,  \
                             tr_y_zzzzz_xz,  \
                             tr_y_zzzzz_y,   \
                             tr_y_zzzzz_yy,  \
                             tr_y_zzzzz_yz,  \
                             tr_y_zzzzz_z,   \
                             tr_y_zzzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzzz_xx[i] = 2.0 * tr_y_zzzzz_x[i] * fe_0 + tr_y_zzzzz_xx[i] * pa_x[i];

        tr_y_xzzzzz_xy[i] = tr_y_zzzzz_y[i] * fe_0 + tr_y_zzzzz_xy[i] * pa_x[i];

        tr_y_xzzzzz_xz[i] = tr_y_zzzzz_z[i] * fe_0 + tr_y_zzzzz_xz[i] * pa_x[i];

        tr_y_xzzzzz_yy[i] = tr_y_zzzzz_yy[i] * pa_x[i];

        tr_y_xzzzzz_yz[i] = tr_y_zzzzz_yz[i] * pa_x[i];

        tr_y_xzzzzz_zz[i] = tr_y_zzzzz_zz[i] * pa_x[i];
    }

    // Set up 294-300 components of targeted buffer : ID

    auto tr_y_yyyyyy_xx = pbuffer.data(idx_dip_id + 294);

    auto tr_y_yyyyyy_xy = pbuffer.data(idx_dip_id + 295);

    auto tr_y_yyyyyy_xz = pbuffer.data(idx_dip_id + 296);

    auto tr_y_yyyyyy_yy = pbuffer.data(idx_dip_id + 297);

    auto tr_y_yyyyyy_yz = pbuffer.data(idx_dip_id + 298);

    auto tr_y_yyyyyy_zz = pbuffer.data(idx_dip_id + 299);

#pragma omp simd aligned(pa_y,               \
                             tr_y_yyyy_xx,   \
                             tr_y_yyyy_xy,   \
                             tr_y_yyyy_xz,   \
                             tr_y_yyyy_yy,   \
                             tr_y_yyyy_yz,   \
                             tr_y_yyyy_zz,   \
                             tr_y_yyyyy_x,   \
                             tr_y_yyyyy_xx,  \
                             tr_y_yyyyy_xy,  \
                             tr_y_yyyyy_xz,  \
                             tr_y_yyyyy_y,   \
                             tr_y_yyyyy_yy,  \
                             tr_y_yyyyy_yz,  \
                             tr_y_yyyyy_z,   \
                             tr_y_yyyyy_zz,  \
                             tr_y_yyyyyy_xx, \
                             tr_y_yyyyyy_xy, \
                             tr_y_yyyyyy_xz, \
                             tr_y_yyyyyy_yy, \
                             tr_y_yyyyyy_yz, \
                             tr_y_yyyyyy_zz, \
                             ts_yyyyy_xx,    \
                             ts_yyyyy_xy,    \
                             ts_yyyyy_xz,    \
                             ts_yyyyy_yy,    \
                             ts_yyyyy_yz,    \
                             ts_yyyyy_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyy_xx[i] = 5.0 * tr_y_yyyy_xx[i] * fe_0 + ts_yyyyy_xx[i] * fe_0 + tr_y_yyyyy_xx[i] * pa_y[i];

        tr_y_yyyyyy_xy[i] = 5.0 * tr_y_yyyy_xy[i] * fe_0 + tr_y_yyyyy_x[i] * fe_0 + ts_yyyyy_xy[i] * fe_0 + tr_y_yyyyy_xy[i] * pa_y[i];

        tr_y_yyyyyy_xz[i] = 5.0 * tr_y_yyyy_xz[i] * fe_0 + ts_yyyyy_xz[i] * fe_0 + tr_y_yyyyy_xz[i] * pa_y[i];

        tr_y_yyyyyy_yy[i] = 5.0 * tr_y_yyyy_yy[i] * fe_0 + 2.0 * tr_y_yyyyy_y[i] * fe_0 + ts_yyyyy_yy[i] * fe_0 + tr_y_yyyyy_yy[i] * pa_y[i];

        tr_y_yyyyyy_yz[i] = 5.0 * tr_y_yyyy_yz[i] * fe_0 + tr_y_yyyyy_z[i] * fe_0 + ts_yyyyy_yz[i] * fe_0 + tr_y_yyyyy_yz[i] * pa_y[i];

        tr_y_yyyyyy_zz[i] = 5.0 * tr_y_yyyy_zz[i] * fe_0 + ts_yyyyy_zz[i] * fe_0 + tr_y_yyyyy_zz[i] * pa_y[i];
    }

    // Set up 300-306 components of targeted buffer : ID

    auto tr_y_yyyyyz_xx = pbuffer.data(idx_dip_id + 300);

    auto tr_y_yyyyyz_xy = pbuffer.data(idx_dip_id + 301);

    auto tr_y_yyyyyz_xz = pbuffer.data(idx_dip_id + 302);

    auto tr_y_yyyyyz_yy = pbuffer.data(idx_dip_id + 303);

    auto tr_y_yyyyyz_yz = pbuffer.data(idx_dip_id + 304);

    auto tr_y_yyyyyz_zz = pbuffer.data(idx_dip_id + 305);

#pragma omp simd aligned(pa_z,               \
                             tr_y_yyyyy_x,   \
                             tr_y_yyyyy_xx,  \
                             tr_y_yyyyy_xy,  \
                             tr_y_yyyyy_xz,  \
                             tr_y_yyyyy_y,   \
                             tr_y_yyyyy_yy,  \
                             tr_y_yyyyy_yz,  \
                             tr_y_yyyyy_z,   \
                             tr_y_yyyyy_zz,  \
                             tr_y_yyyyyz_xx, \
                             tr_y_yyyyyz_xy, \
                             tr_y_yyyyyz_xz, \
                             tr_y_yyyyyz_yy, \
                             tr_y_yyyyyz_yz, \
                             tr_y_yyyyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyz_xx[i] = tr_y_yyyyy_xx[i] * pa_z[i];

        tr_y_yyyyyz_xy[i] = tr_y_yyyyy_xy[i] * pa_z[i];

        tr_y_yyyyyz_xz[i] = tr_y_yyyyy_x[i] * fe_0 + tr_y_yyyyy_xz[i] * pa_z[i];

        tr_y_yyyyyz_yy[i] = tr_y_yyyyy_yy[i] * pa_z[i];

        tr_y_yyyyyz_yz[i] = tr_y_yyyyy_y[i] * fe_0 + tr_y_yyyyy_yz[i] * pa_z[i];

        tr_y_yyyyyz_zz[i] = 2.0 * tr_y_yyyyy_z[i] * fe_0 + tr_y_yyyyy_zz[i] * pa_z[i];
    }

    // Set up 306-312 components of targeted buffer : ID

    auto tr_y_yyyyzz_xx = pbuffer.data(idx_dip_id + 306);

    auto tr_y_yyyyzz_xy = pbuffer.data(idx_dip_id + 307);

    auto tr_y_yyyyzz_xz = pbuffer.data(idx_dip_id + 308);

    auto tr_y_yyyyzz_yy = pbuffer.data(idx_dip_id + 309);

    auto tr_y_yyyyzz_yz = pbuffer.data(idx_dip_id + 310);

    auto tr_y_yyyyzz_zz = pbuffer.data(idx_dip_id + 311);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yyyy_xx,   \
                             tr_y_yyyy_xy,   \
                             tr_y_yyyy_yy,   \
                             tr_y_yyyy_yz,   \
                             tr_y_yyyyz_xx,  \
                             tr_y_yyyyz_xy,  \
                             tr_y_yyyyz_y,   \
                             tr_y_yyyyz_yy,  \
                             tr_y_yyyyz_yz,  \
                             tr_y_yyyyzz_xx, \
                             tr_y_yyyyzz_xy, \
                             tr_y_yyyyzz_xz, \
                             tr_y_yyyyzz_yy, \
                             tr_y_yyyyzz_yz, \
                             tr_y_yyyyzz_zz, \
                             tr_y_yyyzz_xz,  \
                             tr_y_yyyzz_zz,  \
                             tr_y_yyzz_xz,   \
                             tr_y_yyzz_zz,   \
                             ts_yyyzz_xz,    \
                             ts_yyyzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyzz_xx[i] = tr_y_yyyy_xx[i] * fe_0 + tr_y_yyyyz_xx[i] * pa_z[i];

        tr_y_yyyyzz_xy[i] = tr_y_yyyy_xy[i] * fe_0 + tr_y_yyyyz_xy[i] * pa_z[i];

        tr_y_yyyyzz_xz[i] = 3.0 * tr_y_yyzz_xz[i] * fe_0 + ts_yyyzz_xz[i] * fe_0 + tr_y_yyyzz_xz[i] * pa_y[i];

        tr_y_yyyyzz_yy[i] = tr_y_yyyy_yy[i] * fe_0 + tr_y_yyyyz_yy[i] * pa_z[i];

        tr_y_yyyyzz_yz[i] = tr_y_yyyy_yz[i] * fe_0 + tr_y_yyyyz_y[i] * fe_0 + tr_y_yyyyz_yz[i] * pa_z[i];

        tr_y_yyyyzz_zz[i] = 3.0 * tr_y_yyzz_zz[i] * fe_0 + ts_yyyzz_zz[i] * fe_0 + tr_y_yyyzz_zz[i] * pa_y[i];
    }

    // Set up 312-318 components of targeted buffer : ID

    auto tr_y_yyyzzz_xx = pbuffer.data(idx_dip_id + 312);

    auto tr_y_yyyzzz_xy = pbuffer.data(idx_dip_id + 313);

    auto tr_y_yyyzzz_xz = pbuffer.data(idx_dip_id + 314);

    auto tr_y_yyyzzz_yy = pbuffer.data(idx_dip_id + 315);

    auto tr_y_yyyzzz_yz = pbuffer.data(idx_dip_id + 316);

    auto tr_y_yyyzzz_zz = pbuffer.data(idx_dip_id + 317);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yyyz_xx,   \
                             tr_y_yyyz_xy,   \
                             tr_y_yyyz_yy,   \
                             tr_y_yyyz_yz,   \
                             tr_y_yyyzz_xx,  \
                             tr_y_yyyzz_xy,  \
                             tr_y_yyyzz_y,   \
                             tr_y_yyyzz_yy,  \
                             tr_y_yyyzz_yz,  \
                             tr_y_yyyzzz_xx, \
                             tr_y_yyyzzz_xy, \
                             tr_y_yyyzzz_xz, \
                             tr_y_yyyzzz_yy, \
                             tr_y_yyyzzz_yz, \
                             tr_y_yyyzzz_zz, \
                             tr_y_yyzzz_xz,  \
                             tr_y_yyzzz_zz,  \
                             tr_y_yzzz_xz,   \
                             tr_y_yzzz_zz,   \
                             ts_yyzzz_xz,    \
                             ts_yyzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzzz_xx[i] = 2.0 * tr_y_yyyz_xx[i] * fe_0 + tr_y_yyyzz_xx[i] * pa_z[i];

        tr_y_yyyzzz_xy[i] = 2.0 * tr_y_yyyz_xy[i] * fe_0 + tr_y_yyyzz_xy[i] * pa_z[i];

        tr_y_yyyzzz_xz[i] = 2.0 * tr_y_yzzz_xz[i] * fe_0 + ts_yyzzz_xz[i] * fe_0 + tr_y_yyzzz_xz[i] * pa_y[i];

        tr_y_yyyzzz_yy[i] = 2.0 * tr_y_yyyz_yy[i] * fe_0 + tr_y_yyyzz_yy[i] * pa_z[i];

        tr_y_yyyzzz_yz[i] = 2.0 * tr_y_yyyz_yz[i] * fe_0 + tr_y_yyyzz_y[i] * fe_0 + tr_y_yyyzz_yz[i] * pa_z[i];

        tr_y_yyyzzz_zz[i] = 2.0 * tr_y_yzzz_zz[i] * fe_0 + ts_yyzzz_zz[i] * fe_0 + tr_y_yyzzz_zz[i] * pa_y[i];
    }

    // Set up 318-324 components of targeted buffer : ID

    auto tr_y_yyzzzz_xx = pbuffer.data(idx_dip_id + 318);

    auto tr_y_yyzzzz_xy = pbuffer.data(idx_dip_id + 319);

    auto tr_y_yyzzzz_xz = pbuffer.data(idx_dip_id + 320);

    auto tr_y_yyzzzz_yy = pbuffer.data(idx_dip_id + 321);

    auto tr_y_yyzzzz_yz = pbuffer.data(idx_dip_id + 322);

    auto tr_y_yyzzzz_zz = pbuffer.data(idx_dip_id + 323);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yyzz_xx,   \
                             tr_y_yyzz_xy,   \
                             tr_y_yyzz_yy,   \
                             tr_y_yyzz_yz,   \
                             tr_y_yyzzz_xx,  \
                             tr_y_yyzzz_xy,  \
                             tr_y_yyzzz_y,   \
                             tr_y_yyzzz_yy,  \
                             tr_y_yyzzz_yz,  \
                             tr_y_yyzzzz_xx, \
                             tr_y_yyzzzz_xy, \
                             tr_y_yyzzzz_xz, \
                             tr_y_yyzzzz_yy, \
                             tr_y_yyzzzz_yz, \
                             tr_y_yyzzzz_zz, \
                             tr_y_yzzzz_xz,  \
                             tr_y_yzzzz_zz,  \
                             tr_y_zzzz_xz,   \
                             tr_y_zzzz_zz,   \
                             ts_yzzzz_xz,    \
                             ts_yzzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzzz_xx[i] = 3.0 * tr_y_yyzz_xx[i] * fe_0 + tr_y_yyzzz_xx[i] * pa_z[i];

        tr_y_yyzzzz_xy[i] = 3.0 * tr_y_yyzz_xy[i] * fe_0 + tr_y_yyzzz_xy[i] * pa_z[i];

        tr_y_yyzzzz_xz[i] = tr_y_zzzz_xz[i] * fe_0 + ts_yzzzz_xz[i] * fe_0 + tr_y_yzzzz_xz[i] * pa_y[i];

        tr_y_yyzzzz_yy[i] = 3.0 * tr_y_yyzz_yy[i] * fe_0 + tr_y_yyzzz_yy[i] * pa_z[i];

        tr_y_yyzzzz_yz[i] = 3.0 * tr_y_yyzz_yz[i] * fe_0 + tr_y_yyzzz_y[i] * fe_0 + tr_y_yyzzz_yz[i] * pa_z[i];

        tr_y_yyzzzz_zz[i] = tr_y_zzzz_zz[i] * fe_0 + ts_yzzzz_zz[i] * fe_0 + tr_y_yzzzz_zz[i] * pa_y[i];
    }

    // Set up 324-330 components of targeted buffer : ID

    auto tr_y_yzzzzz_xx = pbuffer.data(idx_dip_id + 324);

    auto tr_y_yzzzzz_xy = pbuffer.data(idx_dip_id + 325);

    auto tr_y_yzzzzz_xz = pbuffer.data(idx_dip_id + 326);

    auto tr_y_yzzzzz_yy = pbuffer.data(idx_dip_id + 327);

    auto tr_y_yzzzzz_yz = pbuffer.data(idx_dip_id + 328);

    auto tr_y_yzzzzz_zz = pbuffer.data(idx_dip_id + 329);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yzzz_xy,   \
                             tr_y_yzzz_yy,   \
                             tr_y_yzzzz_xy,  \
                             tr_y_yzzzz_yy,  \
                             tr_y_yzzzzz_xx, \
                             tr_y_yzzzzz_xy, \
                             tr_y_yzzzzz_xz, \
                             tr_y_yzzzzz_yy, \
                             tr_y_yzzzzz_yz, \
                             tr_y_yzzzzz_zz, \
                             tr_y_zzzzz_xx,  \
                             tr_y_zzzzz_xz,  \
                             tr_y_zzzzz_yz,  \
                             tr_y_zzzzz_z,   \
                             tr_y_zzzzz_zz,  \
                             ts_zzzzz_xx,    \
                             ts_zzzzz_xz,    \
                             ts_zzzzz_yz,    \
                             ts_zzzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzzz_xx[i] = ts_zzzzz_xx[i] * fe_0 + tr_y_zzzzz_xx[i] * pa_y[i];

        tr_y_yzzzzz_xy[i] = 4.0 * tr_y_yzzz_xy[i] * fe_0 + tr_y_yzzzz_xy[i] * pa_z[i];

        tr_y_yzzzzz_xz[i] = ts_zzzzz_xz[i] * fe_0 + tr_y_zzzzz_xz[i] * pa_y[i];

        tr_y_yzzzzz_yy[i] = 4.0 * tr_y_yzzz_yy[i] * fe_0 + tr_y_yzzzz_yy[i] * pa_z[i];

        tr_y_yzzzzz_yz[i] = tr_y_zzzzz_z[i] * fe_0 + ts_zzzzz_yz[i] * fe_0 + tr_y_zzzzz_yz[i] * pa_y[i];

        tr_y_yzzzzz_zz[i] = ts_zzzzz_zz[i] * fe_0 + tr_y_zzzzz_zz[i] * pa_y[i];
    }

    // Set up 330-336 components of targeted buffer : ID

    auto tr_y_zzzzzz_xx = pbuffer.data(idx_dip_id + 330);

    auto tr_y_zzzzzz_xy = pbuffer.data(idx_dip_id + 331);

    auto tr_y_zzzzzz_xz = pbuffer.data(idx_dip_id + 332);

    auto tr_y_zzzzzz_yy = pbuffer.data(idx_dip_id + 333);

    auto tr_y_zzzzzz_yz = pbuffer.data(idx_dip_id + 334);

    auto tr_y_zzzzzz_zz = pbuffer.data(idx_dip_id + 335);

#pragma omp simd aligned(pa_z,               \
                             tr_y_zzzz_xx,   \
                             tr_y_zzzz_xy,   \
                             tr_y_zzzz_xz,   \
                             tr_y_zzzz_yy,   \
                             tr_y_zzzz_yz,   \
                             tr_y_zzzz_zz,   \
                             tr_y_zzzzz_x,   \
                             tr_y_zzzzz_xx,  \
                             tr_y_zzzzz_xy,  \
                             tr_y_zzzzz_xz,  \
                             tr_y_zzzzz_y,   \
                             tr_y_zzzzz_yy,  \
                             tr_y_zzzzz_yz,  \
                             tr_y_zzzzz_z,   \
                             tr_y_zzzzz_zz,  \
                             tr_y_zzzzzz_xx, \
                             tr_y_zzzzzz_xy, \
                             tr_y_zzzzzz_xz, \
                             tr_y_zzzzzz_yy, \
                             tr_y_zzzzzz_yz, \
                             tr_y_zzzzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzzz_xx[i] = 5.0 * tr_y_zzzz_xx[i] * fe_0 + tr_y_zzzzz_xx[i] * pa_z[i];

        tr_y_zzzzzz_xy[i] = 5.0 * tr_y_zzzz_xy[i] * fe_0 + tr_y_zzzzz_xy[i] * pa_z[i];

        tr_y_zzzzzz_xz[i] = 5.0 * tr_y_zzzz_xz[i] * fe_0 + tr_y_zzzzz_x[i] * fe_0 + tr_y_zzzzz_xz[i] * pa_z[i];

        tr_y_zzzzzz_yy[i] = 5.0 * tr_y_zzzz_yy[i] * fe_0 + tr_y_zzzzz_yy[i] * pa_z[i];

        tr_y_zzzzzz_yz[i] = 5.0 * tr_y_zzzz_yz[i] * fe_0 + tr_y_zzzzz_y[i] * fe_0 + tr_y_zzzzz_yz[i] * pa_z[i];

        tr_y_zzzzzz_zz[i] = 5.0 * tr_y_zzzz_zz[i] * fe_0 + 2.0 * tr_y_zzzzz_z[i] * fe_0 + tr_y_zzzzz_zz[i] * pa_z[i];
    }

    // Set up 336-342 components of targeted buffer : ID

    auto tr_z_xxxxxx_xx = pbuffer.data(idx_dip_id + 336);

    auto tr_z_xxxxxx_xy = pbuffer.data(idx_dip_id + 337);

    auto tr_z_xxxxxx_xz = pbuffer.data(idx_dip_id + 338);

    auto tr_z_xxxxxx_yy = pbuffer.data(idx_dip_id + 339);

    auto tr_z_xxxxxx_yz = pbuffer.data(idx_dip_id + 340);

    auto tr_z_xxxxxx_zz = pbuffer.data(idx_dip_id + 341);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxxx_xx,   \
                             tr_z_xxxx_xy,   \
                             tr_z_xxxx_xz,   \
                             tr_z_xxxx_yy,   \
                             tr_z_xxxx_yz,   \
                             tr_z_xxxx_zz,   \
                             tr_z_xxxxx_x,   \
                             tr_z_xxxxx_xx,  \
                             tr_z_xxxxx_xy,  \
                             tr_z_xxxxx_xz,  \
                             tr_z_xxxxx_y,   \
                             tr_z_xxxxx_yy,  \
                             tr_z_xxxxx_yz,  \
                             tr_z_xxxxx_z,   \
                             tr_z_xxxxx_zz,  \
                             tr_z_xxxxxx_xx, \
                             tr_z_xxxxxx_xy, \
                             tr_z_xxxxxx_xz, \
                             tr_z_xxxxxx_yy, \
                             tr_z_xxxxxx_yz, \
                             tr_z_xxxxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxx_xx[i] = 5.0 * tr_z_xxxx_xx[i] * fe_0 + 2.0 * tr_z_xxxxx_x[i] * fe_0 + tr_z_xxxxx_xx[i] * pa_x[i];

        tr_z_xxxxxx_xy[i] = 5.0 * tr_z_xxxx_xy[i] * fe_0 + tr_z_xxxxx_y[i] * fe_0 + tr_z_xxxxx_xy[i] * pa_x[i];

        tr_z_xxxxxx_xz[i] = 5.0 * tr_z_xxxx_xz[i] * fe_0 + tr_z_xxxxx_z[i] * fe_0 + tr_z_xxxxx_xz[i] * pa_x[i];

        tr_z_xxxxxx_yy[i] = 5.0 * tr_z_xxxx_yy[i] * fe_0 + tr_z_xxxxx_yy[i] * pa_x[i];

        tr_z_xxxxxx_yz[i] = 5.0 * tr_z_xxxx_yz[i] * fe_0 + tr_z_xxxxx_yz[i] * pa_x[i];

        tr_z_xxxxxx_zz[i] = 5.0 * tr_z_xxxx_zz[i] * fe_0 + tr_z_xxxxx_zz[i] * pa_x[i];
    }

    // Set up 342-348 components of targeted buffer : ID

    auto tr_z_xxxxxy_xx = pbuffer.data(idx_dip_id + 342);

    auto tr_z_xxxxxy_xy = pbuffer.data(idx_dip_id + 343);

    auto tr_z_xxxxxy_xz = pbuffer.data(idx_dip_id + 344);

    auto tr_z_xxxxxy_yy = pbuffer.data(idx_dip_id + 345);

    auto tr_z_xxxxxy_yz = pbuffer.data(idx_dip_id + 346);

    auto tr_z_xxxxxy_zz = pbuffer.data(idx_dip_id + 347);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxxx_x,   \
                             tr_z_xxxxx_xx,  \
                             tr_z_xxxxx_xy,  \
                             tr_z_xxxxx_xz,  \
                             tr_z_xxxxx_zz,  \
                             tr_z_xxxxxy_xx, \
                             tr_z_xxxxxy_xy, \
                             tr_z_xxxxxy_xz, \
                             tr_z_xxxxxy_yy, \
                             tr_z_xxxxxy_yz, \
                             tr_z_xxxxxy_zz, \
                             tr_z_xxxxy_yy,  \
                             tr_z_xxxxy_yz,  \
                             tr_z_xxxy_yy,   \
                             tr_z_xxxy_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxy_xx[i] = tr_z_xxxxx_xx[i] * pa_y[i];

        tr_z_xxxxxy_xy[i] = tr_z_xxxxx_x[i] * fe_0 + tr_z_xxxxx_xy[i] * pa_y[i];

        tr_z_xxxxxy_xz[i] = tr_z_xxxxx_xz[i] * pa_y[i];

        tr_z_xxxxxy_yy[i] = 4.0 * tr_z_xxxy_yy[i] * fe_0 + tr_z_xxxxy_yy[i] * pa_x[i];

        tr_z_xxxxxy_yz[i] = 4.0 * tr_z_xxxy_yz[i] * fe_0 + tr_z_xxxxy_yz[i] * pa_x[i];

        tr_z_xxxxxy_zz[i] = tr_z_xxxxx_zz[i] * pa_y[i];
    }

    // Set up 348-354 components of targeted buffer : ID

    auto tr_z_xxxxxz_xx = pbuffer.data(idx_dip_id + 348);

    auto tr_z_xxxxxz_xy = pbuffer.data(idx_dip_id + 349);

    auto tr_z_xxxxxz_xz = pbuffer.data(idx_dip_id + 350);

    auto tr_z_xxxxxz_yy = pbuffer.data(idx_dip_id + 351);

    auto tr_z_xxxxxz_yz = pbuffer.data(idx_dip_id + 352);

    auto tr_z_xxxxxz_zz = pbuffer.data(idx_dip_id + 353);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_z_xxxxx_xx,  \
                             tr_z_xxxxx_xy,  \
                             tr_z_xxxxxz_xx, \
                             tr_z_xxxxxz_xy, \
                             tr_z_xxxxxz_xz, \
                             tr_z_xxxxxz_yy, \
                             tr_z_xxxxxz_yz, \
                             tr_z_xxxxxz_zz, \
                             tr_z_xxxxz_xz,  \
                             tr_z_xxxxz_yy,  \
                             tr_z_xxxxz_yz,  \
                             tr_z_xxxxz_z,   \
                             tr_z_xxxxz_zz,  \
                             tr_z_xxxz_xz,   \
                             tr_z_xxxz_yy,   \
                             tr_z_xxxz_yz,   \
                             tr_z_xxxz_zz,   \
                             ts_xxxxx_xx,    \
                             ts_xxxxx_xy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxz_xx[i] = ts_xxxxx_xx[i] * fe_0 + tr_z_xxxxx_xx[i] * pa_z[i];

        tr_z_xxxxxz_xy[i] = ts_xxxxx_xy[i] * fe_0 + tr_z_xxxxx_xy[i] * pa_z[i];

        tr_z_xxxxxz_xz[i] = 4.0 * tr_z_xxxz_xz[i] * fe_0 + tr_z_xxxxz_z[i] * fe_0 + tr_z_xxxxz_xz[i] * pa_x[i];

        tr_z_xxxxxz_yy[i] = 4.0 * tr_z_xxxz_yy[i] * fe_0 + tr_z_xxxxz_yy[i] * pa_x[i];

        tr_z_xxxxxz_yz[i] = 4.0 * tr_z_xxxz_yz[i] * fe_0 + tr_z_xxxxz_yz[i] * pa_x[i];

        tr_z_xxxxxz_zz[i] = 4.0 * tr_z_xxxz_zz[i] * fe_0 + tr_z_xxxxz_zz[i] * pa_x[i];
    }

    // Set up 354-360 components of targeted buffer : ID

    auto tr_z_xxxxyy_xx = pbuffer.data(idx_dip_id + 354);

    auto tr_z_xxxxyy_xy = pbuffer.data(idx_dip_id + 355);

    auto tr_z_xxxxyy_xz = pbuffer.data(idx_dip_id + 356);

    auto tr_z_xxxxyy_yy = pbuffer.data(idx_dip_id + 357);

    auto tr_z_xxxxyy_yz = pbuffer.data(idx_dip_id + 358);

    auto tr_z_xxxxyy_zz = pbuffer.data(idx_dip_id + 359);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxx_xx,   \
                             tr_z_xxxx_xz,   \
                             tr_z_xxxxy_xx,  \
                             tr_z_xxxxy_xz,  \
                             tr_z_xxxxyy_xx, \
                             tr_z_xxxxyy_xy, \
                             tr_z_xxxxyy_xz, \
                             tr_z_xxxxyy_yy, \
                             tr_z_xxxxyy_yz, \
                             tr_z_xxxxyy_zz, \
                             tr_z_xxxyy_xy,  \
                             tr_z_xxxyy_y,   \
                             tr_z_xxxyy_yy,  \
                             tr_z_xxxyy_yz,  \
                             tr_z_xxxyy_zz,  \
                             tr_z_xxyy_xy,   \
                             tr_z_xxyy_yy,   \
                             tr_z_xxyy_yz,   \
                             tr_z_xxyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyy_xx[i] = tr_z_xxxx_xx[i] * fe_0 + tr_z_xxxxy_xx[i] * pa_y[i];

        tr_z_xxxxyy_xy[i] = 3.0 * tr_z_xxyy_xy[i] * fe_0 + tr_z_xxxyy_y[i] * fe_0 + tr_z_xxxyy_xy[i] * pa_x[i];

        tr_z_xxxxyy_xz[i] = tr_z_xxxx_xz[i] * fe_0 + tr_z_xxxxy_xz[i] * pa_y[i];

        tr_z_xxxxyy_yy[i] = 3.0 * tr_z_xxyy_yy[i] * fe_0 + tr_z_xxxyy_yy[i] * pa_x[i];

        tr_z_xxxxyy_yz[i] = 3.0 * tr_z_xxyy_yz[i] * fe_0 + tr_z_xxxyy_yz[i] * pa_x[i];

        tr_z_xxxxyy_zz[i] = 3.0 * tr_z_xxyy_zz[i] * fe_0 + tr_z_xxxyy_zz[i] * pa_x[i];
    }

    // Set up 360-366 components of targeted buffer : ID

    auto tr_z_xxxxyz_xx = pbuffer.data(idx_dip_id + 360);

    auto tr_z_xxxxyz_xy = pbuffer.data(idx_dip_id + 361);

    auto tr_z_xxxxyz_xz = pbuffer.data(idx_dip_id + 362);

    auto tr_z_xxxxyz_yy = pbuffer.data(idx_dip_id + 363);

    auto tr_z_xxxxyz_yz = pbuffer.data(idx_dip_id + 364);

    auto tr_z_xxxxyz_zz = pbuffer.data(idx_dip_id + 365);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxxyz_xx, \
                             tr_z_xxxxyz_xy, \
                             tr_z_xxxxyz_xz, \
                             tr_z_xxxxyz_yy, \
                             tr_z_xxxxyz_yz, \
                             tr_z_xxxxyz_zz, \
                             tr_z_xxxxz_x,   \
                             tr_z_xxxxz_xx,  \
                             tr_z_xxxxz_xy,  \
                             tr_z_xxxxz_xz,  \
                             tr_z_xxxxz_zz,  \
                             tr_z_xxxyz_yy,  \
                             tr_z_xxxyz_yz,  \
                             tr_z_xxyz_yy,   \
                             tr_z_xxyz_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyz_xx[i] = tr_z_xxxxz_xx[i] * pa_y[i];

        tr_z_xxxxyz_xy[i] = tr_z_xxxxz_x[i] * fe_0 + tr_z_xxxxz_xy[i] * pa_y[i];

        tr_z_xxxxyz_xz[i] = tr_z_xxxxz_xz[i] * pa_y[i];

        tr_z_xxxxyz_yy[i] = 3.0 * tr_z_xxyz_yy[i] * fe_0 + tr_z_xxxyz_yy[i] * pa_x[i];

        tr_z_xxxxyz_yz[i] = 3.0 * tr_z_xxyz_yz[i] * fe_0 + tr_z_xxxyz_yz[i] * pa_x[i];

        tr_z_xxxxyz_zz[i] = tr_z_xxxxz_zz[i] * pa_y[i];
    }

    // Set up 366-372 components of targeted buffer : ID

    auto tr_z_xxxxzz_xx = pbuffer.data(idx_dip_id + 366);

    auto tr_z_xxxxzz_xy = pbuffer.data(idx_dip_id + 367);

    auto tr_z_xxxxzz_xz = pbuffer.data(idx_dip_id + 368);

    auto tr_z_xxxxzz_yy = pbuffer.data(idx_dip_id + 369);

    auto tr_z_xxxxzz_yz = pbuffer.data(idx_dip_id + 370);

    auto tr_z_xxxxzz_zz = pbuffer.data(idx_dip_id + 371);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxxxzz_xx, \
                             tr_z_xxxxzz_xy, \
                             tr_z_xxxxzz_xz, \
                             tr_z_xxxxzz_yy, \
                             tr_z_xxxxzz_yz, \
                             tr_z_xxxxzz_zz, \
                             tr_z_xxxzz_x,   \
                             tr_z_xxxzz_xx,  \
                             tr_z_xxxzz_xy,  \
                             tr_z_xxxzz_xz,  \
                             tr_z_xxxzz_y,   \
                             tr_z_xxxzz_yy,  \
                             tr_z_xxxzz_yz,  \
                             tr_z_xxxzz_z,   \
                             tr_z_xxxzz_zz,  \
                             tr_z_xxzz_xx,   \
                             tr_z_xxzz_xy,   \
                             tr_z_xxzz_xz,   \
                             tr_z_xxzz_yy,   \
                             tr_z_xxzz_yz,   \
                             tr_z_xxzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxzz_xx[i] = 3.0 * tr_z_xxzz_xx[i] * fe_0 + 2.0 * tr_z_xxxzz_x[i] * fe_0 + tr_z_xxxzz_xx[i] * pa_x[i];

        tr_z_xxxxzz_xy[i] = 3.0 * tr_z_xxzz_xy[i] * fe_0 + tr_z_xxxzz_y[i] * fe_0 + tr_z_xxxzz_xy[i] * pa_x[i];

        tr_z_xxxxzz_xz[i] = 3.0 * tr_z_xxzz_xz[i] * fe_0 + tr_z_xxxzz_z[i] * fe_0 + tr_z_xxxzz_xz[i] * pa_x[i];

        tr_z_xxxxzz_yy[i] = 3.0 * tr_z_xxzz_yy[i] * fe_0 + tr_z_xxxzz_yy[i] * pa_x[i];

        tr_z_xxxxzz_yz[i] = 3.0 * tr_z_xxzz_yz[i] * fe_0 + tr_z_xxxzz_yz[i] * pa_x[i];

        tr_z_xxxxzz_zz[i] = 3.0 * tr_z_xxzz_zz[i] * fe_0 + tr_z_xxxzz_zz[i] * pa_x[i];
    }

    // Set up 372-378 components of targeted buffer : ID

    auto tr_z_xxxyyy_xx = pbuffer.data(idx_dip_id + 372);

    auto tr_z_xxxyyy_xy = pbuffer.data(idx_dip_id + 373);

    auto tr_z_xxxyyy_xz = pbuffer.data(idx_dip_id + 374);

    auto tr_z_xxxyyy_yy = pbuffer.data(idx_dip_id + 375);

    auto tr_z_xxxyyy_yz = pbuffer.data(idx_dip_id + 376);

    auto tr_z_xxxyyy_zz = pbuffer.data(idx_dip_id + 377);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxy_xx,   \
                             tr_z_xxxy_xz,   \
                             tr_z_xxxyy_xx,  \
                             tr_z_xxxyy_xz,  \
                             tr_z_xxxyyy_xx, \
                             tr_z_xxxyyy_xy, \
                             tr_z_xxxyyy_xz, \
                             tr_z_xxxyyy_yy, \
                             tr_z_xxxyyy_yz, \
                             tr_z_xxxyyy_zz, \
                             tr_z_xxyyy_xy,  \
                             tr_z_xxyyy_y,   \
                             tr_z_xxyyy_yy,  \
                             tr_z_xxyyy_yz,  \
                             tr_z_xxyyy_zz,  \
                             tr_z_xyyy_xy,   \
                             tr_z_xyyy_yy,   \
                             tr_z_xyyy_yz,   \
                             tr_z_xyyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyy_xx[i] = 2.0 * tr_z_xxxy_xx[i] * fe_0 + tr_z_xxxyy_xx[i] * pa_y[i];

        tr_z_xxxyyy_xy[i] = 2.0 * tr_z_xyyy_xy[i] * fe_0 + tr_z_xxyyy_y[i] * fe_0 + tr_z_xxyyy_xy[i] * pa_x[i];

        tr_z_xxxyyy_xz[i] = 2.0 * tr_z_xxxy_xz[i] * fe_0 + tr_z_xxxyy_xz[i] * pa_y[i];

        tr_z_xxxyyy_yy[i] = 2.0 * tr_z_xyyy_yy[i] * fe_0 + tr_z_xxyyy_yy[i] * pa_x[i];

        tr_z_xxxyyy_yz[i] = 2.0 * tr_z_xyyy_yz[i] * fe_0 + tr_z_xxyyy_yz[i] * pa_x[i];

        tr_z_xxxyyy_zz[i] = 2.0 * tr_z_xyyy_zz[i] * fe_0 + tr_z_xxyyy_zz[i] * pa_x[i];
    }

    // Set up 378-384 components of targeted buffer : ID

    auto tr_z_xxxyyz_xx = pbuffer.data(idx_dip_id + 378);

    auto tr_z_xxxyyz_xy = pbuffer.data(idx_dip_id + 379);

    auto tr_z_xxxyyz_xz = pbuffer.data(idx_dip_id + 380);

    auto tr_z_xxxyyz_yy = pbuffer.data(idx_dip_id + 381);

    auto tr_z_xxxyyz_yz = pbuffer.data(idx_dip_id + 382);

    auto tr_z_xxxyyz_zz = pbuffer.data(idx_dip_id + 383);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_z_xxxyy_xy,  \
                             tr_z_xxxyyz_xx, \
                             tr_z_xxxyyz_xy, \
                             tr_z_xxxyyz_xz, \
                             tr_z_xxxyyz_yy, \
                             tr_z_xxxyyz_yz, \
                             tr_z_xxxyyz_zz, \
                             tr_z_xxxyz_xx,  \
                             tr_z_xxxyz_xz,  \
                             tr_z_xxxz_xx,   \
                             tr_z_xxxz_xz,   \
                             tr_z_xxyyz_yy,  \
                             tr_z_xxyyz_yz,  \
                             tr_z_xxyyz_zz,  \
                             tr_z_xyyz_yy,   \
                             tr_z_xyyz_yz,   \
                             tr_z_xyyz_zz,   \
                             ts_xxxyy_xy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyz_xx[i] = tr_z_xxxz_xx[i] * fe_0 + tr_z_xxxyz_xx[i] * pa_y[i];

        tr_z_xxxyyz_xy[i] = ts_xxxyy_xy[i] * fe_0 + tr_z_xxxyy_xy[i] * pa_z[i];

        tr_z_xxxyyz_xz[i] = tr_z_xxxz_xz[i] * fe_0 + tr_z_xxxyz_xz[i] * pa_y[i];

        tr_z_xxxyyz_yy[i] = 2.0 * tr_z_xyyz_yy[i] * fe_0 + tr_z_xxyyz_yy[i] * pa_x[i];

        tr_z_xxxyyz_yz[i] = 2.0 * tr_z_xyyz_yz[i] * fe_0 + tr_z_xxyyz_yz[i] * pa_x[i];

        tr_z_xxxyyz_zz[i] = 2.0 * tr_z_xyyz_zz[i] * fe_0 + tr_z_xxyyz_zz[i] * pa_x[i];
    }

    // Set up 384-390 components of targeted buffer : ID

    auto tr_z_xxxyzz_xx = pbuffer.data(idx_dip_id + 384);

    auto tr_z_xxxyzz_xy = pbuffer.data(idx_dip_id + 385);

    auto tr_z_xxxyzz_xz = pbuffer.data(idx_dip_id + 386);

    auto tr_z_xxxyzz_yy = pbuffer.data(idx_dip_id + 387);

    auto tr_z_xxxyzz_yz = pbuffer.data(idx_dip_id + 388);

    auto tr_z_xxxyzz_zz = pbuffer.data(idx_dip_id + 389);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxyzz_xx, \
                             tr_z_xxxyzz_xy, \
                             tr_z_xxxyzz_xz, \
                             tr_z_xxxyzz_yy, \
                             tr_z_xxxyzz_yz, \
                             tr_z_xxxyzz_zz, \
                             tr_z_xxxzz_x,   \
                             tr_z_xxxzz_xx,  \
                             tr_z_xxxzz_xy,  \
                             tr_z_xxxzz_xz,  \
                             tr_z_xxxzz_zz,  \
                             tr_z_xxyzz_yy,  \
                             tr_z_xxyzz_yz,  \
                             tr_z_xyzz_yy,   \
                             tr_z_xyzz_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyzz_xx[i] = tr_z_xxxzz_xx[i] * pa_y[i];

        tr_z_xxxyzz_xy[i] = tr_z_xxxzz_x[i] * fe_0 + tr_z_xxxzz_xy[i] * pa_y[i];

        tr_z_xxxyzz_xz[i] = tr_z_xxxzz_xz[i] * pa_y[i];

        tr_z_xxxyzz_yy[i] = 2.0 * tr_z_xyzz_yy[i] * fe_0 + tr_z_xxyzz_yy[i] * pa_x[i];

        tr_z_xxxyzz_yz[i] = 2.0 * tr_z_xyzz_yz[i] * fe_0 + tr_z_xxyzz_yz[i] * pa_x[i];

        tr_z_xxxyzz_zz[i] = tr_z_xxxzz_zz[i] * pa_y[i];
    }

    // Set up 390-396 components of targeted buffer : ID

    auto tr_z_xxxzzz_xx = pbuffer.data(idx_dip_id + 390);

    auto tr_z_xxxzzz_xy = pbuffer.data(idx_dip_id + 391);

    auto tr_z_xxxzzz_xz = pbuffer.data(idx_dip_id + 392);

    auto tr_z_xxxzzz_yy = pbuffer.data(idx_dip_id + 393);

    auto tr_z_xxxzzz_yz = pbuffer.data(idx_dip_id + 394);

    auto tr_z_xxxzzz_zz = pbuffer.data(idx_dip_id + 395);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxxzzz_xx, \
                             tr_z_xxxzzz_xy, \
                             tr_z_xxxzzz_xz, \
                             tr_z_xxxzzz_yy, \
                             tr_z_xxxzzz_yz, \
                             tr_z_xxxzzz_zz, \
                             tr_z_xxzzz_x,   \
                             tr_z_xxzzz_xx,  \
                             tr_z_xxzzz_xy,  \
                             tr_z_xxzzz_xz,  \
                             tr_z_xxzzz_y,   \
                             tr_z_xxzzz_yy,  \
                             tr_z_xxzzz_yz,  \
                             tr_z_xxzzz_z,   \
                             tr_z_xxzzz_zz,  \
                             tr_z_xzzz_xx,   \
                             tr_z_xzzz_xy,   \
                             tr_z_xzzz_xz,   \
                             tr_z_xzzz_yy,   \
                             tr_z_xzzz_yz,   \
                             tr_z_xzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzzz_xx[i] = 2.0 * tr_z_xzzz_xx[i] * fe_0 + 2.0 * tr_z_xxzzz_x[i] * fe_0 + tr_z_xxzzz_xx[i] * pa_x[i];

        tr_z_xxxzzz_xy[i] = 2.0 * tr_z_xzzz_xy[i] * fe_0 + tr_z_xxzzz_y[i] * fe_0 + tr_z_xxzzz_xy[i] * pa_x[i];

        tr_z_xxxzzz_xz[i] = 2.0 * tr_z_xzzz_xz[i] * fe_0 + tr_z_xxzzz_z[i] * fe_0 + tr_z_xxzzz_xz[i] * pa_x[i];

        tr_z_xxxzzz_yy[i] = 2.0 * tr_z_xzzz_yy[i] * fe_0 + tr_z_xxzzz_yy[i] * pa_x[i];

        tr_z_xxxzzz_yz[i] = 2.0 * tr_z_xzzz_yz[i] * fe_0 + tr_z_xxzzz_yz[i] * pa_x[i];

        tr_z_xxxzzz_zz[i] = 2.0 * tr_z_xzzz_zz[i] * fe_0 + tr_z_xxzzz_zz[i] * pa_x[i];
    }

    // Set up 396-402 components of targeted buffer : ID

    auto tr_z_xxyyyy_xx = pbuffer.data(idx_dip_id + 396);

    auto tr_z_xxyyyy_xy = pbuffer.data(idx_dip_id + 397);

    auto tr_z_xxyyyy_xz = pbuffer.data(idx_dip_id + 398);

    auto tr_z_xxyyyy_yy = pbuffer.data(idx_dip_id + 399);

    auto tr_z_xxyyyy_yz = pbuffer.data(idx_dip_id + 400);

    auto tr_z_xxyyyy_zz = pbuffer.data(idx_dip_id + 401);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxyy_xx,   \
                             tr_z_xxyy_xz,   \
                             tr_z_xxyyy_xx,  \
                             tr_z_xxyyy_xz,  \
                             tr_z_xxyyyy_xx, \
                             tr_z_xxyyyy_xy, \
                             tr_z_xxyyyy_xz, \
                             tr_z_xxyyyy_yy, \
                             tr_z_xxyyyy_yz, \
                             tr_z_xxyyyy_zz, \
                             tr_z_xyyyy_xy,  \
                             tr_z_xyyyy_y,   \
                             tr_z_xyyyy_yy,  \
                             tr_z_xyyyy_yz,  \
                             tr_z_xyyyy_zz,  \
                             tr_z_yyyy_xy,   \
                             tr_z_yyyy_yy,   \
                             tr_z_yyyy_yz,   \
                             tr_z_yyyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyy_xx[i] = 3.0 * tr_z_xxyy_xx[i] * fe_0 + tr_z_xxyyy_xx[i] * pa_y[i];

        tr_z_xxyyyy_xy[i] = tr_z_yyyy_xy[i] * fe_0 + tr_z_xyyyy_y[i] * fe_0 + tr_z_xyyyy_xy[i] * pa_x[i];

        tr_z_xxyyyy_xz[i] = 3.0 * tr_z_xxyy_xz[i] * fe_0 + tr_z_xxyyy_xz[i] * pa_y[i];

        tr_z_xxyyyy_yy[i] = tr_z_yyyy_yy[i] * fe_0 + tr_z_xyyyy_yy[i] * pa_x[i];

        tr_z_xxyyyy_yz[i] = tr_z_yyyy_yz[i] * fe_0 + tr_z_xyyyy_yz[i] * pa_x[i];

        tr_z_xxyyyy_zz[i] = tr_z_yyyy_zz[i] * fe_0 + tr_z_xyyyy_zz[i] * pa_x[i];
    }

    // Set up 402-408 components of targeted buffer : ID

    auto tr_z_xxyyyz_xx = pbuffer.data(idx_dip_id + 402);

    auto tr_z_xxyyyz_xy = pbuffer.data(idx_dip_id + 403);

    auto tr_z_xxyyyz_xz = pbuffer.data(idx_dip_id + 404);

    auto tr_z_xxyyyz_yy = pbuffer.data(idx_dip_id + 405);

    auto tr_z_xxyyyz_yz = pbuffer.data(idx_dip_id + 406);

    auto tr_z_xxyyyz_zz = pbuffer.data(idx_dip_id + 407);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_z_xxyyy_xy,  \
                             tr_z_xxyyyz_xx, \
                             tr_z_xxyyyz_xy, \
                             tr_z_xxyyyz_xz, \
                             tr_z_xxyyyz_yy, \
                             tr_z_xxyyyz_yz, \
                             tr_z_xxyyyz_zz, \
                             tr_z_xxyyz_xx,  \
                             tr_z_xxyyz_xz,  \
                             tr_z_xxyz_xx,   \
                             tr_z_xxyz_xz,   \
                             tr_z_xyyyz_yy,  \
                             tr_z_xyyyz_yz,  \
                             tr_z_xyyyz_zz,  \
                             tr_z_yyyz_yy,   \
                             tr_z_yyyz_yz,   \
                             tr_z_yyyz_zz,   \
                             ts_xxyyy_xy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyz_xx[i] = 2.0 * tr_z_xxyz_xx[i] * fe_0 + tr_z_xxyyz_xx[i] * pa_y[i];

        tr_z_xxyyyz_xy[i] = ts_xxyyy_xy[i] * fe_0 + tr_z_xxyyy_xy[i] * pa_z[i];

        tr_z_xxyyyz_xz[i] = 2.0 * tr_z_xxyz_xz[i] * fe_0 + tr_z_xxyyz_xz[i] * pa_y[i];

        tr_z_xxyyyz_yy[i] = tr_z_yyyz_yy[i] * fe_0 + tr_z_xyyyz_yy[i] * pa_x[i];

        tr_z_xxyyyz_yz[i] = tr_z_yyyz_yz[i] * fe_0 + tr_z_xyyyz_yz[i] * pa_x[i];

        tr_z_xxyyyz_zz[i] = tr_z_yyyz_zz[i] * fe_0 + tr_z_xyyyz_zz[i] * pa_x[i];
    }

    // Set up 408-414 components of targeted buffer : ID

    auto tr_z_xxyyzz_xx = pbuffer.data(idx_dip_id + 408);

    auto tr_z_xxyyzz_xy = pbuffer.data(idx_dip_id + 409);

    auto tr_z_xxyyzz_xz = pbuffer.data(idx_dip_id + 410);

    auto tr_z_xxyyzz_yy = pbuffer.data(idx_dip_id + 411);

    auto tr_z_xxyyzz_yz = pbuffer.data(idx_dip_id + 412);

    auto tr_z_xxyyzz_zz = pbuffer.data(idx_dip_id + 413);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxyyzz_xx, \
                             tr_z_xxyyzz_xy, \
                             tr_z_xxyyzz_xz, \
                             tr_z_xxyyzz_yy, \
                             tr_z_xxyyzz_yz, \
                             tr_z_xxyyzz_zz, \
                             tr_z_xxyzz_xx,  \
                             tr_z_xxyzz_xz,  \
                             tr_z_xxzz_xx,   \
                             tr_z_xxzz_xz,   \
                             tr_z_xyyzz_xy,  \
                             tr_z_xyyzz_y,   \
                             tr_z_xyyzz_yy,  \
                             tr_z_xyyzz_yz,  \
                             tr_z_xyyzz_zz,  \
                             tr_z_yyzz_xy,   \
                             tr_z_yyzz_yy,   \
                             tr_z_yyzz_yz,   \
                             tr_z_yyzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyzz_xx[i] = tr_z_xxzz_xx[i] * fe_0 + tr_z_xxyzz_xx[i] * pa_y[i];

        tr_z_xxyyzz_xy[i] = tr_z_yyzz_xy[i] * fe_0 + tr_z_xyyzz_y[i] * fe_0 + tr_z_xyyzz_xy[i] * pa_x[i];

        tr_z_xxyyzz_xz[i] = tr_z_xxzz_xz[i] * fe_0 + tr_z_xxyzz_xz[i] * pa_y[i];

        tr_z_xxyyzz_yy[i] = tr_z_yyzz_yy[i] * fe_0 + tr_z_xyyzz_yy[i] * pa_x[i];

        tr_z_xxyyzz_yz[i] = tr_z_yyzz_yz[i] * fe_0 + tr_z_xyyzz_yz[i] * pa_x[i];

        tr_z_xxyyzz_zz[i] = tr_z_yyzz_zz[i] * fe_0 + tr_z_xyyzz_zz[i] * pa_x[i];
    }

    // Set up 414-420 components of targeted buffer : ID

    auto tr_z_xxyzzz_xx = pbuffer.data(idx_dip_id + 414);

    auto tr_z_xxyzzz_xy = pbuffer.data(idx_dip_id + 415);

    auto tr_z_xxyzzz_xz = pbuffer.data(idx_dip_id + 416);

    auto tr_z_xxyzzz_yy = pbuffer.data(idx_dip_id + 417);

    auto tr_z_xxyzzz_yz = pbuffer.data(idx_dip_id + 418);

    auto tr_z_xxyzzz_zz = pbuffer.data(idx_dip_id + 419);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxyzzz_xx, \
                             tr_z_xxyzzz_xy, \
                             tr_z_xxyzzz_xz, \
                             tr_z_xxyzzz_yy, \
                             tr_z_xxyzzz_yz, \
                             tr_z_xxyzzz_zz, \
                             tr_z_xxzzz_x,   \
                             tr_z_xxzzz_xx,  \
                             tr_z_xxzzz_xy,  \
                             tr_z_xxzzz_xz,  \
                             tr_z_xxzzz_zz,  \
                             tr_z_xyzzz_yy,  \
                             tr_z_xyzzz_yz,  \
                             tr_z_yzzz_yy,   \
                             tr_z_yzzz_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzzz_xx[i] = tr_z_xxzzz_xx[i] * pa_y[i];

        tr_z_xxyzzz_xy[i] = tr_z_xxzzz_x[i] * fe_0 + tr_z_xxzzz_xy[i] * pa_y[i];

        tr_z_xxyzzz_xz[i] = tr_z_xxzzz_xz[i] * pa_y[i];

        tr_z_xxyzzz_yy[i] = tr_z_yzzz_yy[i] * fe_0 + tr_z_xyzzz_yy[i] * pa_x[i];

        tr_z_xxyzzz_yz[i] = tr_z_yzzz_yz[i] * fe_0 + tr_z_xyzzz_yz[i] * pa_x[i];

        tr_z_xxyzzz_zz[i] = tr_z_xxzzz_zz[i] * pa_y[i];
    }

    // Set up 420-426 components of targeted buffer : ID

    auto tr_z_xxzzzz_xx = pbuffer.data(idx_dip_id + 420);

    auto tr_z_xxzzzz_xy = pbuffer.data(idx_dip_id + 421);

    auto tr_z_xxzzzz_xz = pbuffer.data(idx_dip_id + 422);

    auto tr_z_xxzzzz_yy = pbuffer.data(idx_dip_id + 423);

    auto tr_z_xxzzzz_yz = pbuffer.data(idx_dip_id + 424);

    auto tr_z_xxzzzz_zz = pbuffer.data(idx_dip_id + 425);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxzzzz_xx, \
                             tr_z_xxzzzz_xy, \
                             tr_z_xxzzzz_xz, \
                             tr_z_xxzzzz_yy, \
                             tr_z_xxzzzz_yz, \
                             tr_z_xxzzzz_zz, \
                             tr_z_xzzzz_x,   \
                             tr_z_xzzzz_xx,  \
                             tr_z_xzzzz_xy,  \
                             tr_z_xzzzz_xz,  \
                             tr_z_xzzzz_y,   \
                             tr_z_xzzzz_yy,  \
                             tr_z_xzzzz_yz,  \
                             tr_z_xzzzz_z,   \
                             tr_z_xzzzz_zz,  \
                             tr_z_zzzz_xx,   \
                             tr_z_zzzz_xy,   \
                             tr_z_zzzz_xz,   \
                             tr_z_zzzz_yy,   \
                             tr_z_zzzz_yz,   \
                             tr_z_zzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzzz_xx[i] = tr_z_zzzz_xx[i] * fe_0 + 2.0 * tr_z_xzzzz_x[i] * fe_0 + tr_z_xzzzz_xx[i] * pa_x[i];

        tr_z_xxzzzz_xy[i] = tr_z_zzzz_xy[i] * fe_0 + tr_z_xzzzz_y[i] * fe_0 + tr_z_xzzzz_xy[i] * pa_x[i];

        tr_z_xxzzzz_xz[i] = tr_z_zzzz_xz[i] * fe_0 + tr_z_xzzzz_z[i] * fe_0 + tr_z_xzzzz_xz[i] * pa_x[i];

        tr_z_xxzzzz_yy[i] = tr_z_zzzz_yy[i] * fe_0 + tr_z_xzzzz_yy[i] * pa_x[i];

        tr_z_xxzzzz_yz[i] = tr_z_zzzz_yz[i] * fe_0 + tr_z_xzzzz_yz[i] * pa_x[i];

        tr_z_xxzzzz_zz[i] = tr_z_zzzz_zz[i] * fe_0 + tr_z_xzzzz_zz[i] * pa_x[i];
    }

    // Set up 426-432 components of targeted buffer : ID

    auto tr_z_xyyyyy_xx = pbuffer.data(idx_dip_id + 426);

    auto tr_z_xyyyyy_xy = pbuffer.data(idx_dip_id + 427);

    auto tr_z_xyyyyy_xz = pbuffer.data(idx_dip_id + 428);

    auto tr_z_xyyyyy_yy = pbuffer.data(idx_dip_id + 429);

    auto tr_z_xyyyyy_yz = pbuffer.data(idx_dip_id + 430);

    auto tr_z_xyyyyy_zz = pbuffer.data(idx_dip_id + 431);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyyyy_xx, \
                             tr_z_xyyyyy_xy, \
                             tr_z_xyyyyy_xz, \
                             tr_z_xyyyyy_yy, \
                             tr_z_xyyyyy_yz, \
                             tr_z_xyyyyy_zz, \
                             tr_z_yyyyy_x,   \
                             tr_z_yyyyy_xx,  \
                             tr_z_yyyyy_xy,  \
                             tr_z_yyyyy_xz,  \
                             tr_z_yyyyy_y,   \
                             tr_z_yyyyy_yy,  \
                             tr_z_yyyyy_yz,  \
                             tr_z_yyyyy_z,   \
                             tr_z_yyyyy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyy_xx[i] = 2.0 * tr_z_yyyyy_x[i] * fe_0 + tr_z_yyyyy_xx[i] * pa_x[i];

        tr_z_xyyyyy_xy[i] = tr_z_yyyyy_y[i] * fe_0 + tr_z_yyyyy_xy[i] * pa_x[i];

        tr_z_xyyyyy_xz[i] = tr_z_yyyyy_z[i] * fe_0 + tr_z_yyyyy_xz[i] * pa_x[i];

        tr_z_xyyyyy_yy[i] = tr_z_yyyyy_yy[i] * pa_x[i];

        tr_z_xyyyyy_yz[i] = tr_z_yyyyy_yz[i] * pa_x[i];

        tr_z_xyyyyy_zz[i] = tr_z_yyyyy_zz[i] * pa_x[i];
    }

    // Set up 432-438 components of targeted buffer : ID

    auto tr_z_xyyyyz_xx = pbuffer.data(idx_dip_id + 432);

    auto tr_z_xyyyyz_xy = pbuffer.data(idx_dip_id + 433);

    auto tr_z_xyyyyz_xz = pbuffer.data(idx_dip_id + 434);

    auto tr_z_xyyyyz_yy = pbuffer.data(idx_dip_id + 435);

    auto tr_z_xyyyyz_yz = pbuffer.data(idx_dip_id + 436);

    auto tr_z_xyyyyz_zz = pbuffer.data(idx_dip_id + 437);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyyyz_xx, \
                             tr_z_xyyyyz_xy, \
                             tr_z_xyyyyz_xz, \
                             tr_z_xyyyyz_yy, \
                             tr_z_xyyyyz_yz, \
                             tr_z_xyyyyz_zz, \
                             tr_z_yyyyz_x,   \
                             tr_z_yyyyz_xx,  \
                             tr_z_yyyyz_xy,  \
                             tr_z_yyyyz_xz,  \
                             tr_z_yyyyz_y,   \
                             tr_z_yyyyz_yy,  \
                             tr_z_yyyyz_yz,  \
                             tr_z_yyyyz_z,   \
                             tr_z_yyyyz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyz_xx[i] = 2.0 * tr_z_yyyyz_x[i] * fe_0 + tr_z_yyyyz_xx[i] * pa_x[i];

        tr_z_xyyyyz_xy[i] = tr_z_yyyyz_y[i] * fe_0 + tr_z_yyyyz_xy[i] * pa_x[i];

        tr_z_xyyyyz_xz[i] = tr_z_yyyyz_z[i] * fe_0 + tr_z_yyyyz_xz[i] * pa_x[i];

        tr_z_xyyyyz_yy[i] = tr_z_yyyyz_yy[i] * pa_x[i];

        tr_z_xyyyyz_yz[i] = tr_z_yyyyz_yz[i] * pa_x[i];

        tr_z_xyyyyz_zz[i] = tr_z_yyyyz_zz[i] * pa_x[i];
    }

    // Set up 438-444 components of targeted buffer : ID

    auto tr_z_xyyyzz_xx = pbuffer.data(idx_dip_id + 438);

    auto tr_z_xyyyzz_xy = pbuffer.data(idx_dip_id + 439);

    auto tr_z_xyyyzz_xz = pbuffer.data(idx_dip_id + 440);

    auto tr_z_xyyyzz_yy = pbuffer.data(idx_dip_id + 441);

    auto tr_z_xyyyzz_yz = pbuffer.data(idx_dip_id + 442);

    auto tr_z_xyyyzz_zz = pbuffer.data(idx_dip_id + 443);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyyzz_xx, \
                             tr_z_xyyyzz_xy, \
                             tr_z_xyyyzz_xz, \
                             tr_z_xyyyzz_yy, \
                             tr_z_xyyyzz_yz, \
                             tr_z_xyyyzz_zz, \
                             tr_z_yyyzz_x,   \
                             tr_z_yyyzz_xx,  \
                             tr_z_yyyzz_xy,  \
                             tr_z_yyyzz_xz,  \
                             tr_z_yyyzz_y,   \
                             tr_z_yyyzz_yy,  \
                             tr_z_yyyzz_yz,  \
                             tr_z_yyyzz_z,   \
                             tr_z_yyyzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyzz_xx[i] = 2.0 * tr_z_yyyzz_x[i] * fe_0 + tr_z_yyyzz_xx[i] * pa_x[i];

        tr_z_xyyyzz_xy[i] = tr_z_yyyzz_y[i] * fe_0 + tr_z_yyyzz_xy[i] * pa_x[i];

        tr_z_xyyyzz_xz[i] = tr_z_yyyzz_z[i] * fe_0 + tr_z_yyyzz_xz[i] * pa_x[i];

        tr_z_xyyyzz_yy[i] = tr_z_yyyzz_yy[i] * pa_x[i];

        tr_z_xyyyzz_yz[i] = tr_z_yyyzz_yz[i] * pa_x[i];

        tr_z_xyyyzz_zz[i] = tr_z_yyyzz_zz[i] * pa_x[i];
    }

    // Set up 444-450 components of targeted buffer : ID

    auto tr_z_xyyzzz_xx = pbuffer.data(idx_dip_id + 444);

    auto tr_z_xyyzzz_xy = pbuffer.data(idx_dip_id + 445);

    auto tr_z_xyyzzz_xz = pbuffer.data(idx_dip_id + 446);

    auto tr_z_xyyzzz_yy = pbuffer.data(idx_dip_id + 447);

    auto tr_z_xyyzzz_yz = pbuffer.data(idx_dip_id + 448);

    auto tr_z_xyyzzz_zz = pbuffer.data(idx_dip_id + 449);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyzzz_xx, \
                             tr_z_xyyzzz_xy, \
                             tr_z_xyyzzz_xz, \
                             tr_z_xyyzzz_yy, \
                             tr_z_xyyzzz_yz, \
                             tr_z_xyyzzz_zz, \
                             tr_z_yyzzz_x,   \
                             tr_z_yyzzz_xx,  \
                             tr_z_yyzzz_xy,  \
                             tr_z_yyzzz_xz,  \
                             tr_z_yyzzz_y,   \
                             tr_z_yyzzz_yy,  \
                             tr_z_yyzzz_yz,  \
                             tr_z_yyzzz_z,   \
                             tr_z_yyzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzzz_xx[i] = 2.0 * tr_z_yyzzz_x[i] * fe_0 + tr_z_yyzzz_xx[i] * pa_x[i];

        tr_z_xyyzzz_xy[i] = tr_z_yyzzz_y[i] * fe_0 + tr_z_yyzzz_xy[i] * pa_x[i];

        tr_z_xyyzzz_xz[i] = tr_z_yyzzz_z[i] * fe_0 + tr_z_yyzzz_xz[i] * pa_x[i];

        tr_z_xyyzzz_yy[i] = tr_z_yyzzz_yy[i] * pa_x[i];

        tr_z_xyyzzz_yz[i] = tr_z_yyzzz_yz[i] * pa_x[i];

        tr_z_xyyzzz_zz[i] = tr_z_yyzzz_zz[i] * pa_x[i];
    }

    // Set up 450-456 components of targeted buffer : ID

    auto tr_z_xyzzzz_xx = pbuffer.data(idx_dip_id + 450);

    auto tr_z_xyzzzz_xy = pbuffer.data(idx_dip_id + 451);

    auto tr_z_xyzzzz_xz = pbuffer.data(idx_dip_id + 452);

    auto tr_z_xyzzzz_yy = pbuffer.data(idx_dip_id + 453);

    auto tr_z_xyzzzz_yz = pbuffer.data(idx_dip_id + 454);

    auto tr_z_xyzzzz_zz = pbuffer.data(idx_dip_id + 455);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xyzzzz_xx, \
                             tr_z_xyzzzz_xy, \
                             tr_z_xyzzzz_xz, \
                             tr_z_xyzzzz_yy, \
                             tr_z_xyzzzz_yz, \
                             tr_z_xyzzzz_zz, \
                             tr_z_xzzzz_xx,  \
                             tr_z_xzzzz_xz,  \
                             tr_z_yzzzz_xy,  \
                             tr_z_yzzzz_y,   \
                             tr_z_yzzzz_yy,  \
                             tr_z_yzzzz_yz,  \
                             tr_z_yzzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzzz_xx[i] = tr_z_xzzzz_xx[i] * pa_y[i];

        tr_z_xyzzzz_xy[i] = tr_z_yzzzz_y[i] * fe_0 + tr_z_yzzzz_xy[i] * pa_x[i];

        tr_z_xyzzzz_xz[i] = tr_z_xzzzz_xz[i] * pa_y[i];

        tr_z_xyzzzz_yy[i] = tr_z_yzzzz_yy[i] * pa_x[i];

        tr_z_xyzzzz_yz[i] = tr_z_yzzzz_yz[i] * pa_x[i];

        tr_z_xyzzzz_zz[i] = tr_z_yzzzz_zz[i] * pa_x[i];
    }

    // Set up 456-462 components of targeted buffer : ID

    auto tr_z_xzzzzz_xx = pbuffer.data(idx_dip_id + 456);

    auto tr_z_xzzzzz_xy = pbuffer.data(idx_dip_id + 457);

    auto tr_z_xzzzzz_xz = pbuffer.data(idx_dip_id + 458);

    auto tr_z_xzzzzz_yy = pbuffer.data(idx_dip_id + 459);

    auto tr_z_xzzzzz_yz = pbuffer.data(idx_dip_id + 460);

    auto tr_z_xzzzzz_zz = pbuffer.data(idx_dip_id + 461);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xzzzzz_xx, \
                             tr_z_xzzzzz_xy, \
                             tr_z_xzzzzz_xz, \
                             tr_z_xzzzzz_yy, \
                             tr_z_xzzzzz_yz, \
                             tr_z_xzzzzz_zz, \
                             tr_z_zzzzz_x,   \
                             tr_z_zzzzz_xx,  \
                             tr_z_zzzzz_xy,  \
                             tr_z_zzzzz_xz,  \
                             tr_z_zzzzz_y,   \
                             tr_z_zzzzz_yy,  \
                             tr_z_zzzzz_yz,  \
                             tr_z_zzzzz_z,   \
                             tr_z_zzzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzzz_xx[i] = 2.0 * tr_z_zzzzz_x[i] * fe_0 + tr_z_zzzzz_xx[i] * pa_x[i];

        tr_z_xzzzzz_xy[i] = tr_z_zzzzz_y[i] * fe_0 + tr_z_zzzzz_xy[i] * pa_x[i];

        tr_z_xzzzzz_xz[i] = tr_z_zzzzz_z[i] * fe_0 + tr_z_zzzzz_xz[i] * pa_x[i];

        tr_z_xzzzzz_yy[i] = tr_z_zzzzz_yy[i] * pa_x[i];

        tr_z_xzzzzz_yz[i] = tr_z_zzzzz_yz[i] * pa_x[i];

        tr_z_xzzzzz_zz[i] = tr_z_zzzzz_zz[i] * pa_x[i];
    }

    // Set up 462-468 components of targeted buffer : ID

    auto tr_z_yyyyyy_xx = pbuffer.data(idx_dip_id + 462);

    auto tr_z_yyyyyy_xy = pbuffer.data(idx_dip_id + 463);

    auto tr_z_yyyyyy_xz = pbuffer.data(idx_dip_id + 464);

    auto tr_z_yyyyyy_yy = pbuffer.data(idx_dip_id + 465);

    auto tr_z_yyyyyy_yz = pbuffer.data(idx_dip_id + 466);

    auto tr_z_yyyyyy_zz = pbuffer.data(idx_dip_id + 467);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyyy_xx,   \
                             tr_z_yyyy_xy,   \
                             tr_z_yyyy_xz,   \
                             tr_z_yyyy_yy,   \
                             tr_z_yyyy_yz,   \
                             tr_z_yyyy_zz,   \
                             tr_z_yyyyy_x,   \
                             tr_z_yyyyy_xx,  \
                             tr_z_yyyyy_xy,  \
                             tr_z_yyyyy_xz,  \
                             tr_z_yyyyy_y,   \
                             tr_z_yyyyy_yy,  \
                             tr_z_yyyyy_yz,  \
                             tr_z_yyyyy_z,   \
                             tr_z_yyyyy_zz,  \
                             tr_z_yyyyyy_xx, \
                             tr_z_yyyyyy_xy, \
                             tr_z_yyyyyy_xz, \
                             tr_z_yyyyyy_yy, \
                             tr_z_yyyyyy_yz, \
                             tr_z_yyyyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyy_xx[i] = 5.0 * tr_z_yyyy_xx[i] * fe_0 + tr_z_yyyyy_xx[i] * pa_y[i];

        tr_z_yyyyyy_xy[i] = 5.0 * tr_z_yyyy_xy[i] * fe_0 + tr_z_yyyyy_x[i] * fe_0 + tr_z_yyyyy_xy[i] * pa_y[i];

        tr_z_yyyyyy_xz[i] = 5.0 * tr_z_yyyy_xz[i] * fe_0 + tr_z_yyyyy_xz[i] * pa_y[i];

        tr_z_yyyyyy_yy[i] = 5.0 * tr_z_yyyy_yy[i] * fe_0 + 2.0 * tr_z_yyyyy_y[i] * fe_0 + tr_z_yyyyy_yy[i] * pa_y[i];

        tr_z_yyyyyy_yz[i] = 5.0 * tr_z_yyyy_yz[i] * fe_0 + tr_z_yyyyy_z[i] * fe_0 + tr_z_yyyyy_yz[i] * pa_y[i];

        tr_z_yyyyyy_zz[i] = 5.0 * tr_z_yyyy_zz[i] * fe_0 + tr_z_yyyyy_zz[i] * pa_y[i];
    }

    // Set up 468-474 components of targeted buffer : ID

    auto tr_z_yyyyyz_xx = pbuffer.data(idx_dip_id + 468);

    auto tr_z_yyyyyz_xy = pbuffer.data(idx_dip_id + 469);

    auto tr_z_yyyyyz_xz = pbuffer.data(idx_dip_id + 470);

    auto tr_z_yyyyyz_yy = pbuffer.data(idx_dip_id + 471);

    auto tr_z_yyyyyz_yz = pbuffer.data(idx_dip_id + 472);

    auto tr_z_yyyyyz_zz = pbuffer.data(idx_dip_id + 473);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_z_yyyyy_xy,  \
                             tr_z_yyyyy_yy,  \
                             tr_z_yyyyyz_xx, \
                             tr_z_yyyyyz_xy, \
                             tr_z_yyyyyz_xz, \
                             tr_z_yyyyyz_yy, \
                             tr_z_yyyyyz_yz, \
                             tr_z_yyyyyz_zz, \
                             tr_z_yyyyz_xx,  \
                             tr_z_yyyyz_xz,  \
                             tr_z_yyyyz_yz,  \
                             tr_z_yyyyz_z,   \
                             tr_z_yyyyz_zz,  \
                             tr_z_yyyz_xx,   \
                             tr_z_yyyz_xz,   \
                             tr_z_yyyz_yz,   \
                             tr_z_yyyz_zz,   \
                             ts_yyyyy_xy,    \
                             ts_yyyyy_yy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyz_xx[i] = 4.0 * tr_z_yyyz_xx[i] * fe_0 + tr_z_yyyyz_xx[i] * pa_y[i];

        tr_z_yyyyyz_xy[i] = ts_yyyyy_xy[i] * fe_0 + tr_z_yyyyy_xy[i] * pa_z[i];

        tr_z_yyyyyz_xz[i] = 4.0 * tr_z_yyyz_xz[i] * fe_0 + tr_z_yyyyz_xz[i] * pa_y[i];

        tr_z_yyyyyz_yy[i] = ts_yyyyy_yy[i] * fe_0 + tr_z_yyyyy_yy[i] * pa_z[i];

        tr_z_yyyyyz_yz[i] = 4.0 * tr_z_yyyz_yz[i] * fe_0 + tr_z_yyyyz_z[i] * fe_0 + tr_z_yyyyz_yz[i] * pa_y[i];

        tr_z_yyyyyz_zz[i] = 4.0 * tr_z_yyyz_zz[i] * fe_0 + tr_z_yyyyz_zz[i] * pa_y[i];
    }

    // Set up 474-480 components of targeted buffer : ID

    auto tr_z_yyyyzz_xx = pbuffer.data(idx_dip_id + 474);

    auto tr_z_yyyyzz_xy = pbuffer.data(idx_dip_id + 475);

    auto tr_z_yyyyzz_xz = pbuffer.data(idx_dip_id + 476);

    auto tr_z_yyyyzz_yy = pbuffer.data(idx_dip_id + 477);

    auto tr_z_yyyyzz_yz = pbuffer.data(idx_dip_id + 478);

    auto tr_z_yyyyzz_zz = pbuffer.data(idx_dip_id + 479);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyyyzz_xx, \
                             tr_z_yyyyzz_xy, \
                             tr_z_yyyyzz_xz, \
                             tr_z_yyyyzz_yy, \
                             tr_z_yyyyzz_yz, \
                             tr_z_yyyyzz_zz, \
                             tr_z_yyyzz_x,   \
                             tr_z_yyyzz_xx,  \
                             tr_z_yyyzz_xy,  \
                             tr_z_yyyzz_xz,  \
                             tr_z_yyyzz_y,   \
                             tr_z_yyyzz_yy,  \
                             tr_z_yyyzz_yz,  \
                             tr_z_yyyzz_z,   \
                             tr_z_yyyzz_zz,  \
                             tr_z_yyzz_xx,   \
                             tr_z_yyzz_xy,   \
                             tr_z_yyzz_xz,   \
                             tr_z_yyzz_yy,   \
                             tr_z_yyzz_yz,   \
                             tr_z_yyzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyzz_xx[i] = 3.0 * tr_z_yyzz_xx[i] * fe_0 + tr_z_yyyzz_xx[i] * pa_y[i];

        tr_z_yyyyzz_xy[i] = 3.0 * tr_z_yyzz_xy[i] * fe_0 + tr_z_yyyzz_x[i] * fe_0 + tr_z_yyyzz_xy[i] * pa_y[i];

        tr_z_yyyyzz_xz[i] = 3.0 * tr_z_yyzz_xz[i] * fe_0 + tr_z_yyyzz_xz[i] * pa_y[i];

        tr_z_yyyyzz_yy[i] = 3.0 * tr_z_yyzz_yy[i] * fe_0 + 2.0 * tr_z_yyyzz_y[i] * fe_0 + tr_z_yyyzz_yy[i] * pa_y[i];

        tr_z_yyyyzz_yz[i] = 3.0 * tr_z_yyzz_yz[i] * fe_0 + tr_z_yyyzz_z[i] * fe_0 + tr_z_yyyzz_yz[i] * pa_y[i];

        tr_z_yyyyzz_zz[i] = 3.0 * tr_z_yyzz_zz[i] * fe_0 + tr_z_yyyzz_zz[i] * pa_y[i];
    }

    // Set up 480-486 components of targeted buffer : ID

    auto tr_z_yyyzzz_xx = pbuffer.data(idx_dip_id + 480);

    auto tr_z_yyyzzz_xy = pbuffer.data(idx_dip_id + 481);

    auto tr_z_yyyzzz_xz = pbuffer.data(idx_dip_id + 482);

    auto tr_z_yyyzzz_yy = pbuffer.data(idx_dip_id + 483);

    auto tr_z_yyyzzz_yz = pbuffer.data(idx_dip_id + 484);

    auto tr_z_yyyzzz_zz = pbuffer.data(idx_dip_id + 485);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyyzzz_xx, \
                             tr_z_yyyzzz_xy, \
                             tr_z_yyyzzz_xz, \
                             tr_z_yyyzzz_yy, \
                             tr_z_yyyzzz_yz, \
                             tr_z_yyyzzz_zz, \
                             tr_z_yyzzz_x,   \
                             tr_z_yyzzz_xx,  \
                             tr_z_yyzzz_xy,  \
                             tr_z_yyzzz_xz,  \
                             tr_z_yyzzz_y,   \
                             tr_z_yyzzz_yy,  \
                             tr_z_yyzzz_yz,  \
                             tr_z_yyzzz_z,   \
                             tr_z_yyzzz_zz,  \
                             tr_z_yzzz_xx,   \
                             tr_z_yzzz_xy,   \
                             tr_z_yzzz_xz,   \
                             tr_z_yzzz_yy,   \
                             tr_z_yzzz_yz,   \
                             tr_z_yzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzzz_xx[i] = 2.0 * tr_z_yzzz_xx[i] * fe_0 + tr_z_yyzzz_xx[i] * pa_y[i];

        tr_z_yyyzzz_xy[i] = 2.0 * tr_z_yzzz_xy[i] * fe_0 + tr_z_yyzzz_x[i] * fe_0 + tr_z_yyzzz_xy[i] * pa_y[i];

        tr_z_yyyzzz_xz[i] = 2.0 * tr_z_yzzz_xz[i] * fe_0 + tr_z_yyzzz_xz[i] * pa_y[i];

        tr_z_yyyzzz_yy[i] = 2.0 * tr_z_yzzz_yy[i] * fe_0 + 2.0 * tr_z_yyzzz_y[i] * fe_0 + tr_z_yyzzz_yy[i] * pa_y[i];

        tr_z_yyyzzz_yz[i] = 2.0 * tr_z_yzzz_yz[i] * fe_0 + tr_z_yyzzz_z[i] * fe_0 + tr_z_yyzzz_yz[i] * pa_y[i];

        tr_z_yyyzzz_zz[i] = 2.0 * tr_z_yzzz_zz[i] * fe_0 + tr_z_yyzzz_zz[i] * pa_y[i];
    }

    // Set up 486-492 components of targeted buffer : ID

    auto tr_z_yyzzzz_xx = pbuffer.data(idx_dip_id + 486);

    auto tr_z_yyzzzz_xy = pbuffer.data(idx_dip_id + 487);

    auto tr_z_yyzzzz_xz = pbuffer.data(idx_dip_id + 488);

    auto tr_z_yyzzzz_yy = pbuffer.data(idx_dip_id + 489);

    auto tr_z_yyzzzz_yz = pbuffer.data(idx_dip_id + 490);

    auto tr_z_yyzzzz_zz = pbuffer.data(idx_dip_id + 491);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyzzzz_xx, \
                             tr_z_yyzzzz_xy, \
                             tr_z_yyzzzz_xz, \
                             tr_z_yyzzzz_yy, \
                             tr_z_yyzzzz_yz, \
                             tr_z_yyzzzz_zz, \
                             tr_z_yzzzz_x,   \
                             tr_z_yzzzz_xx,  \
                             tr_z_yzzzz_xy,  \
                             tr_z_yzzzz_xz,  \
                             tr_z_yzzzz_y,   \
                             tr_z_yzzzz_yy,  \
                             tr_z_yzzzz_yz,  \
                             tr_z_yzzzz_z,   \
                             tr_z_yzzzz_zz,  \
                             tr_z_zzzz_xx,   \
                             tr_z_zzzz_xy,   \
                             tr_z_zzzz_xz,   \
                             tr_z_zzzz_yy,   \
                             tr_z_zzzz_yz,   \
                             tr_z_zzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzzz_xx[i] = tr_z_zzzz_xx[i] * fe_0 + tr_z_yzzzz_xx[i] * pa_y[i];

        tr_z_yyzzzz_xy[i] = tr_z_zzzz_xy[i] * fe_0 + tr_z_yzzzz_x[i] * fe_0 + tr_z_yzzzz_xy[i] * pa_y[i];

        tr_z_yyzzzz_xz[i] = tr_z_zzzz_xz[i] * fe_0 + tr_z_yzzzz_xz[i] * pa_y[i];

        tr_z_yyzzzz_yy[i] = tr_z_zzzz_yy[i] * fe_0 + 2.0 * tr_z_yzzzz_y[i] * fe_0 + tr_z_yzzzz_yy[i] * pa_y[i];

        tr_z_yyzzzz_yz[i] = tr_z_zzzz_yz[i] * fe_0 + tr_z_yzzzz_z[i] * fe_0 + tr_z_yzzzz_yz[i] * pa_y[i];

        tr_z_yyzzzz_zz[i] = tr_z_zzzz_zz[i] * fe_0 + tr_z_yzzzz_zz[i] * pa_y[i];
    }

    // Set up 492-498 components of targeted buffer : ID

    auto tr_z_yzzzzz_xx = pbuffer.data(idx_dip_id + 492);

    auto tr_z_yzzzzz_xy = pbuffer.data(idx_dip_id + 493);

    auto tr_z_yzzzzz_xz = pbuffer.data(idx_dip_id + 494);

    auto tr_z_yzzzzz_yy = pbuffer.data(idx_dip_id + 495);

    auto tr_z_yzzzzz_yz = pbuffer.data(idx_dip_id + 496);

    auto tr_z_yzzzzz_zz = pbuffer.data(idx_dip_id + 497);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yzzzzz_xx, \
                             tr_z_yzzzzz_xy, \
                             tr_z_yzzzzz_xz, \
                             tr_z_yzzzzz_yy, \
                             tr_z_yzzzzz_yz, \
                             tr_z_yzzzzz_zz, \
                             tr_z_zzzzz_x,   \
                             tr_z_zzzzz_xx,  \
                             tr_z_zzzzz_xy,  \
                             tr_z_zzzzz_xz,  \
                             tr_z_zzzzz_y,   \
                             tr_z_zzzzz_yy,  \
                             tr_z_zzzzz_yz,  \
                             tr_z_zzzzz_z,   \
                             tr_z_zzzzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzzz_xx[i] = tr_z_zzzzz_xx[i] * pa_y[i];

        tr_z_yzzzzz_xy[i] = tr_z_zzzzz_x[i] * fe_0 + tr_z_zzzzz_xy[i] * pa_y[i];

        tr_z_yzzzzz_xz[i] = tr_z_zzzzz_xz[i] * pa_y[i];

        tr_z_yzzzzz_yy[i] = 2.0 * tr_z_zzzzz_y[i] * fe_0 + tr_z_zzzzz_yy[i] * pa_y[i];

        tr_z_yzzzzz_yz[i] = tr_z_zzzzz_z[i] * fe_0 + tr_z_zzzzz_yz[i] * pa_y[i];

        tr_z_yzzzzz_zz[i] = tr_z_zzzzz_zz[i] * pa_y[i];
    }

    // Set up 498-504 components of targeted buffer : ID

    auto tr_z_zzzzzz_xx = pbuffer.data(idx_dip_id + 498);

    auto tr_z_zzzzzz_xy = pbuffer.data(idx_dip_id + 499);

    auto tr_z_zzzzzz_xz = pbuffer.data(idx_dip_id + 500);

    auto tr_z_zzzzzz_yy = pbuffer.data(idx_dip_id + 501);

    auto tr_z_zzzzzz_yz = pbuffer.data(idx_dip_id + 502);

    auto tr_z_zzzzzz_zz = pbuffer.data(idx_dip_id + 503);

#pragma omp simd aligned(pa_z,               \
                             tr_z_zzzz_xx,   \
                             tr_z_zzzz_xy,   \
                             tr_z_zzzz_xz,   \
                             tr_z_zzzz_yy,   \
                             tr_z_zzzz_yz,   \
                             tr_z_zzzz_zz,   \
                             tr_z_zzzzz_x,   \
                             tr_z_zzzzz_xx,  \
                             tr_z_zzzzz_xy,  \
                             tr_z_zzzzz_xz,  \
                             tr_z_zzzzz_y,   \
                             tr_z_zzzzz_yy,  \
                             tr_z_zzzzz_yz,  \
                             tr_z_zzzzz_z,   \
                             tr_z_zzzzz_zz,  \
                             tr_z_zzzzzz_xx, \
                             tr_z_zzzzzz_xy, \
                             tr_z_zzzzzz_xz, \
                             tr_z_zzzzzz_yy, \
                             tr_z_zzzzzz_yz, \
                             tr_z_zzzzzz_zz, \
                             ts_zzzzz_xx,    \
                             ts_zzzzz_xy,    \
                             ts_zzzzz_xz,    \
                             ts_zzzzz_yy,    \
                             ts_zzzzz_yz,    \
                             ts_zzzzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzzz_xx[i] = 5.0 * tr_z_zzzz_xx[i] * fe_0 + ts_zzzzz_xx[i] * fe_0 + tr_z_zzzzz_xx[i] * pa_z[i];

        tr_z_zzzzzz_xy[i] = 5.0 * tr_z_zzzz_xy[i] * fe_0 + ts_zzzzz_xy[i] * fe_0 + tr_z_zzzzz_xy[i] * pa_z[i];

        tr_z_zzzzzz_xz[i] = 5.0 * tr_z_zzzz_xz[i] * fe_0 + tr_z_zzzzz_x[i] * fe_0 + ts_zzzzz_xz[i] * fe_0 + tr_z_zzzzz_xz[i] * pa_z[i];

        tr_z_zzzzzz_yy[i] = 5.0 * tr_z_zzzz_yy[i] * fe_0 + ts_zzzzz_yy[i] * fe_0 + tr_z_zzzzz_yy[i] * pa_z[i];

        tr_z_zzzzzz_yz[i] = 5.0 * tr_z_zzzz_yz[i] * fe_0 + tr_z_zzzzz_y[i] * fe_0 + ts_zzzzz_yz[i] * fe_0 + tr_z_zzzzz_yz[i] * pa_z[i];

        tr_z_zzzzzz_zz[i] = 5.0 * tr_z_zzzz_zz[i] * fe_0 + 2.0 * tr_z_zzzzz_z[i] * fe_0 + ts_zzzzz_zz[i] * fe_0 + tr_z_zzzzz_zz[i] * pa_z[i];
    }
}

}  // namespace diprec
