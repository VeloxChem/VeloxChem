#include "NuclearPotentialGeom010PrimRecHG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_hg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_hg,
                                        const size_t              idx_npot_geom_010_0_fg,
                                        const size_t              idx_npot_geom_010_1_fg,
                                        const size_t              idx_npot_geom_010_0_gf,
                                        const size_t              idx_npot_geom_010_1_gf,
                                        const size_t              idx_npot_1_gg,
                                        const size_t              idx_npot_geom_010_0_gg,
                                        const size_t              idx_npot_geom_010_1_gg,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
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

    // Set up components of auxiliary buffer : FG

    auto ta1_x_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg);

    auto ta1_x_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 1);

    auto ta1_x_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 2);

    auto ta1_x_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 3);

    auto ta1_x_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 4);

    auto ta1_x_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 5);

    auto ta1_x_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 6);

    auto ta1_x_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 7);

    auto ta1_x_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 8);

    auto ta1_x_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 9);

    auto ta1_x_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 10);

    auto ta1_x_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 11);

    auto ta1_x_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 12);

    auto ta1_x_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 13);

    auto ta1_x_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 14);

    auto ta1_x_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 15);

    auto ta1_x_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 16);

    auto ta1_x_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 17);

    auto ta1_x_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 18);

    auto ta1_x_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 19);

    auto ta1_x_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 20);

    auto ta1_x_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 21);

    auto ta1_x_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 22);

    auto ta1_x_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 23);

    auto ta1_x_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 24);

    auto ta1_x_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 29);

    auto ta1_x_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 30);

    auto ta1_x_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 31);

    auto ta1_x_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 32);

    auto ta1_x_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 33);

    auto ta1_x_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 34);

    auto ta1_x_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 35);

    auto ta1_x_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 36);

    auto ta1_x_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 37);

    auto ta1_x_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 38);

    auto ta1_x_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 39);

    auto ta1_x_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 40);

    auto ta1_x_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 44);

    auto ta1_x_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 45);

    auto ta1_x_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 46);

    auto ta1_x_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 47);

    auto ta1_x_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 48);

    auto ta1_x_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 50);

    auto ta1_x_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 51);

    auto ta1_x_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 54);

    auto ta1_x_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 55);

    auto ta1_x_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 56);

    auto ta1_x_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 57);

    auto ta1_x_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 58);

    auto ta1_x_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 62);

    auto ta1_x_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 65);

    auto ta1_x_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 69);

    auto ta1_x_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 75);

    auto ta1_x_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 76);

    auto ta1_x_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 77);

    auto ta1_x_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 78);

    auto ta1_x_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 80);

    auto ta1_x_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 81);

    auto ta1_x_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 84);

    auto ta1_x_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 86);

    auto ta1_x_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 87);

    auto ta1_x_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 88);

    auto ta1_x_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 89);

    auto ta1_x_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 90);

    auto ta1_x_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 91);

    auto ta1_x_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 92);

    auto ta1_x_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 93);

    auto ta1_x_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 94);

    auto ta1_x_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 95);

    auto ta1_x_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 96);

    auto ta1_x_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 97);

    auto ta1_x_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 98);

    auto ta1_x_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 99);

    auto ta1_x_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 100);

    auto ta1_x_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 101);

    auto ta1_x_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 102);

    auto ta1_x_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 103);

    auto ta1_x_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 104);

    auto ta1_x_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 106);

    auto ta1_x_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 107);

    auto ta1_x_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 108);

    auto ta1_x_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 110);

    auto ta1_x_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 111);

    auto ta1_x_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 114);

    auto ta1_x_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 115);

    auto ta1_x_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 119);

    auto ta1_x_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 120);

    auto ta1_x_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 122);

    auto ta1_x_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 124);

    auto ta1_x_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 125);

    auto ta1_x_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 127);

    auto ta1_x_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 128);

    auto ta1_x_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 129);

    auto ta1_x_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 131);

    auto ta1_x_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 132);

    auto ta1_x_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 133);

    auto ta1_x_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 134);

    auto ta1_x_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 135);

    auto ta1_x_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 136);

    auto ta1_x_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 137);

    auto ta1_x_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 138);

    auto ta1_x_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 139);

    auto ta1_x_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 140);

    auto ta1_x_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 141);

    auto ta1_x_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 142);

    auto ta1_x_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 143);

    auto ta1_x_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 144);

    auto ta1_x_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 145);

    auto ta1_x_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 146);

    auto ta1_x_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 147);

    auto ta1_x_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 148);

    auto ta1_x_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 149);

    auto ta1_y_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 150);

    auto ta1_y_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 151);

    auto ta1_y_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 152);

    auto ta1_y_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 153);

    auto ta1_y_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 154);

    auto ta1_y_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 155);

    auto ta1_y_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 156);

    auto ta1_y_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 157);

    auto ta1_y_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 158);

    auto ta1_y_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 159);

    auto ta1_y_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 160);

    auto ta1_y_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 161);

    auto ta1_y_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 162);

    auto ta1_y_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 163);

    auto ta1_y_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 164);

    auto ta1_y_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 165);

    auto ta1_y_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 166);

    auto ta1_y_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 167);

    auto ta1_y_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 168);

    auto ta1_y_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 170);

    auto ta1_y_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 171);

    auto ta1_y_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 174);

    auto ta1_y_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 175);

    auto ta1_y_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 176);

    auto ta1_y_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 177);

    auto ta1_y_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 178);

    auto ta1_y_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 180);

    auto ta1_y_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 181);

    auto ta1_y_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 183);

    auto ta1_y_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 186);

    auto ta1_y_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 191);

    auto ta1_y_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 192);

    auto ta1_y_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 193);

    auto ta1_y_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 194);

    auto ta1_y_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 196);

    auto ta1_y_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 198);

    auto ta1_y_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 199);

    auto ta1_y_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 201);

    auto ta1_y_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 202);

    auto ta1_y_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 203);

    auto ta1_y_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 205);

    auto ta1_y_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 206);

    auto ta1_y_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 207);

    auto ta1_y_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 208);

    auto ta1_y_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 209);

    auto ta1_y_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 221);

    auto ta1_y_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 222);

    auto ta1_y_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 223);

    auto ta1_y_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 227);

    auto ta1_y_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 229);

    auto ta1_y_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 230);

    auto ta1_y_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 232);

    auto ta1_y_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 233);

    auto ta1_y_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 234);

    auto ta1_y_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 235);

    auto ta1_y_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 236);

    auto ta1_y_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 237);

    auto ta1_y_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 238);

    auto ta1_y_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 239);

    auto ta1_y_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 240);

    auto ta1_y_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 241);

    auto ta1_y_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 242);

    auto ta1_y_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 243);

    auto ta1_y_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 244);

    auto ta1_y_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 245);

    auto ta1_y_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 246);

    auto ta1_y_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 247);

    auto ta1_y_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 248);

    auto ta1_y_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 249);

    auto ta1_y_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 250);

    auto ta1_y_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 251);

    auto ta1_y_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 252);

    auto ta1_y_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 253);

    auto ta1_y_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 254);

    auto ta1_y_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 255);

    auto ta1_y_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 256);

    auto ta1_y_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 258);

    auto ta1_y_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 259);

    auto ta1_y_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 261);

    auto ta1_y_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 262);

    auto ta1_y_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 263);

    auto ta1_y_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 265);

    auto ta1_y_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 266);

    auto ta1_y_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 267);

    auto ta1_y_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 268);

    auto ta1_y_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 269);

    auto ta1_y_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 271);

    auto ta1_y_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 272);

    auto ta1_y_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 273);

    auto ta1_y_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 275);

    auto ta1_y_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 276);

    auto ta1_y_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 279);

    auto ta1_y_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 280);

    auto ta1_y_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 281);

    auto ta1_y_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 282);

    auto ta1_y_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 283);

    auto ta1_y_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 284);

    auto ta1_y_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 285);

    auto ta1_y_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 286);

    auto ta1_y_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 287);

    auto ta1_y_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 288);

    auto ta1_y_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 289);

    auto ta1_y_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 290);

    auto ta1_y_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 291);

    auto ta1_y_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 292);

    auto ta1_y_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 293);

    auto ta1_y_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 294);

    auto ta1_y_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 295);

    auto ta1_y_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 296);

    auto ta1_y_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 297);

    auto ta1_y_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 298);

    auto ta1_y_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 299);

    auto ta1_z_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 300);

    auto ta1_z_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 301);

    auto ta1_z_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 302);

    auto ta1_z_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 303);

    auto ta1_z_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 304);

    auto ta1_z_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 305);

    auto ta1_z_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 306);

    auto ta1_z_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 307);

    auto ta1_z_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 308);

    auto ta1_z_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 309);

    auto ta1_z_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 310);

    auto ta1_z_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 311);

    auto ta1_z_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 312);

    auto ta1_z_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 313);

    auto ta1_z_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 314);

    auto ta1_z_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 315);

    auto ta1_z_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 317);

    auto ta1_z_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 320);

    auto ta1_z_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 324);

    auto ta1_z_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 325);

    auto ta1_z_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 326);

    auto ta1_z_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 327);

    auto ta1_z_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 328);

    auto ta1_z_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 330);

    auto ta1_z_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 331);

    auto ta1_z_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 332);

    auto ta1_z_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 333);

    auto ta1_z_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 335);

    auto ta1_z_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 336);

    auto ta1_z_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 339);

    auto ta1_z_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 341);

    auto ta1_z_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 342);

    auto ta1_z_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 343);

    auto ta1_z_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 344);

    auto ta1_z_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 346);

    auto ta1_z_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 348);

    auto ta1_z_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 349);

    auto ta1_z_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 351);

    auto ta1_z_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 352);

    auto ta1_z_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 353);

    auto ta1_z_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 355);

    auto ta1_z_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 356);

    auto ta1_z_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 357);

    auto ta1_z_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 358);

    auto ta1_z_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 359);

    auto ta1_z_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 371);

    auto ta1_z_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 372);

    auto ta1_z_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 373);

    auto ta1_z_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 377);

    auto ta1_z_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 379);

    auto ta1_z_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 380);

    auto ta1_z_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 382);

    auto ta1_z_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 383);

    auto ta1_z_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 384);

    auto ta1_z_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 385);

    auto ta1_z_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 386);

    auto ta1_z_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 387);

    auto ta1_z_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 388);

    auto ta1_z_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 389);

    auto ta1_z_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 390);

    auto ta1_z_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 391);

    auto ta1_z_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 392);

    auto ta1_z_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 393);

    auto ta1_z_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 394);

    auto ta1_z_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 395);

    auto ta1_z_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 396);

    auto ta1_z_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 397);

    auto ta1_z_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 398);

    auto ta1_z_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 399);

    auto ta1_z_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 400);

    auto ta1_z_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 401);

    auto ta1_z_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 402);

    auto ta1_z_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 403);

    auto ta1_z_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 404);

    auto ta1_z_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 406);

    auto ta1_z_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 407);

    auto ta1_z_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 408);

    auto ta1_z_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 410);

    auto ta1_z_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 411);

    auto ta1_z_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 414);

    auto ta1_z_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 415);

    auto ta1_z_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 416);

    auto ta1_z_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 417);

    auto ta1_z_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 418);

    auto ta1_z_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 419);

    auto ta1_z_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 420);

    auto ta1_z_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 422);

    auto ta1_z_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 424);

    auto ta1_z_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 425);

    auto ta1_z_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 427);

    auto ta1_z_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 428);

    auto ta1_z_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 429);

    auto ta1_z_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 430);

    auto ta1_z_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 431);

    auto ta1_z_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 432);

    auto ta1_z_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 433);

    auto ta1_z_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 434);

    auto ta1_z_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 435);

    auto ta1_z_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 436);

    auto ta1_z_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 437);

    auto ta1_z_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 438);

    auto ta1_z_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 439);

    auto ta1_z_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 440);

    auto ta1_z_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 441);

    auto ta1_z_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 442);

    auto ta1_z_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 443);

    auto ta1_z_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 444);

    auto ta1_z_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 445);

    auto ta1_z_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 446);

    auto ta1_z_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 447);

    auto ta1_z_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 448);

    auto ta1_z_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 449);

    // Set up components of auxiliary buffer : FG

    auto ta1_x_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg);

    auto ta1_x_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 1);

    auto ta1_x_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 2);

    auto ta1_x_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 3);

    auto ta1_x_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 4);

    auto ta1_x_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 5);

    auto ta1_x_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 6);

    auto ta1_x_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 7);

    auto ta1_x_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 8);

    auto ta1_x_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 9);

    auto ta1_x_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 10);

    auto ta1_x_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 11);

    auto ta1_x_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 12);

    auto ta1_x_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 13);

    auto ta1_x_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 14);

    auto ta1_x_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 15);

    auto ta1_x_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 16);

    auto ta1_x_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 17);

    auto ta1_x_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 18);

    auto ta1_x_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 19);

    auto ta1_x_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 20);

    auto ta1_x_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 21);

    auto ta1_x_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 22);

    auto ta1_x_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 23);

    auto ta1_x_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 24);

    auto ta1_x_xxy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 29);

    auto ta1_x_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 30);

    auto ta1_x_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 31);

    auto ta1_x_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 32);

    auto ta1_x_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 33);

    auto ta1_x_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 34);

    auto ta1_x_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 35);

    auto ta1_x_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 36);

    auto ta1_x_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 37);

    auto ta1_x_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 38);

    auto ta1_x_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 39);

    auto ta1_x_xxz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 40);

    auto ta1_x_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 44);

    auto ta1_x_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 45);

    auto ta1_x_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 46);

    auto ta1_x_xyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 47);

    auto ta1_x_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 48);

    auto ta1_x_xyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 50);

    auto ta1_x_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 51);

    auto ta1_x_xyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 54);

    auto ta1_x_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 55);

    auto ta1_x_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 56);

    auto ta1_x_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 57);

    auto ta1_x_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 58);

    auto ta1_x_xyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 62);

    auto ta1_x_xyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 65);

    auto ta1_x_xyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 69);

    auto ta1_x_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 75);

    auto ta1_x_xzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 76);

    auto ta1_x_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 77);

    auto ta1_x_xzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 78);

    auto ta1_x_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 80);

    auto ta1_x_xzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 81);

    auto ta1_x_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 84);

    auto ta1_x_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 86);

    auto ta1_x_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 87);

    auto ta1_x_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 88);

    auto ta1_x_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 89);

    auto ta1_x_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 90);

    auto ta1_x_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 91);

    auto ta1_x_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 92);

    auto ta1_x_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 93);

    auto ta1_x_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 94);

    auto ta1_x_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 95);

    auto ta1_x_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 96);

    auto ta1_x_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 97);

    auto ta1_x_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 98);

    auto ta1_x_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 99);

    auto ta1_x_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 100);

    auto ta1_x_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 101);

    auto ta1_x_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 102);

    auto ta1_x_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 103);

    auto ta1_x_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 104);

    auto ta1_x_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 106);

    auto ta1_x_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 107);

    auto ta1_x_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 108);

    auto ta1_x_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 110);

    auto ta1_x_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 111);

    auto ta1_x_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 114);

    auto ta1_x_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 115);

    auto ta1_x_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 119);

    auto ta1_x_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 120);

    auto ta1_x_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 122);

    auto ta1_x_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 124);

    auto ta1_x_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 125);

    auto ta1_x_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 127);

    auto ta1_x_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 128);

    auto ta1_x_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 129);

    auto ta1_x_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 131);

    auto ta1_x_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 132);

    auto ta1_x_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 133);

    auto ta1_x_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 134);

    auto ta1_x_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 135);

    auto ta1_x_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 136);

    auto ta1_x_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 137);

    auto ta1_x_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 138);

    auto ta1_x_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 139);

    auto ta1_x_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 140);

    auto ta1_x_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 141);

    auto ta1_x_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 142);

    auto ta1_x_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 143);

    auto ta1_x_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 144);

    auto ta1_x_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 145);

    auto ta1_x_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 146);

    auto ta1_x_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 147);

    auto ta1_x_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 148);

    auto ta1_x_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 149);

    auto ta1_y_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 150);

    auto ta1_y_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 151);

    auto ta1_y_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 152);

    auto ta1_y_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 153);

    auto ta1_y_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 154);

    auto ta1_y_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 155);

    auto ta1_y_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 156);

    auto ta1_y_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 157);

    auto ta1_y_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 158);

    auto ta1_y_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 159);

    auto ta1_y_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 160);

    auto ta1_y_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 161);

    auto ta1_y_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 162);

    auto ta1_y_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 163);

    auto ta1_y_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 164);

    auto ta1_y_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 165);

    auto ta1_y_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 166);

    auto ta1_y_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 167);

    auto ta1_y_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 168);

    auto ta1_y_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 170);

    auto ta1_y_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 171);

    auto ta1_y_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 174);

    auto ta1_y_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 175);

    auto ta1_y_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 176);

    auto ta1_y_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 177);

    auto ta1_y_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 178);

    auto ta1_y_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 180);

    auto ta1_y_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 181);

    auto ta1_y_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 183);

    auto ta1_y_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 186);

    auto ta1_y_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 191);

    auto ta1_y_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 192);

    auto ta1_y_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 193);

    auto ta1_y_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 194);

    auto ta1_y_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 196);

    auto ta1_y_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 198);

    auto ta1_y_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 199);

    auto ta1_y_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 201);

    auto ta1_y_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 202);

    auto ta1_y_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 203);

    auto ta1_y_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 205);

    auto ta1_y_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 206);

    auto ta1_y_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 207);

    auto ta1_y_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 208);

    auto ta1_y_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 209);

    auto ta1_y_xyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 221);

    auto ta1_y_xyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 222);

    auto ta1_y_xyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 223);

    auto ta1_y_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 227);

    auto ta1_y_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 229);

    auto ta1_y_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 230);

    auto ta1_y_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 232);

    auto ta1_y_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 233);

    auto ta1_y_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 234);

    auto ta1_y_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 235);

    auto ta1_y_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 236);

    auto ta1_y_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 237);

    auto ta1_y_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 238);

    auto ta1_y_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 239);

    auto ta1_y_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 240);

    auto ta1_y_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 241);

    auto ta1_y_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 242);

    auto ta1_y_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 243);

    auto ta1_y_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 244);

    auto ta1_y_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 245);

    auto ta1_y_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 246);

    auto ta1_y_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 247);

    auto ta1_y_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 248);

    auto ta1_y_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 249);

    auto ta1_y_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 250);

    auto ta1_y_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 251);

    auto ta1_y_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 252);

    auto ta1_y_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 253);

    auto ta1_y_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 254);

    auto ta1_y_yyz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 255);

    auto ta1_y_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 256);

    auto ta1_y_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 258);

    auto ta1_y_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 259);

    auto ta1_y_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 261);

    auto ta1_y_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 262);

    auto ta1_y_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 263);

    auto ta1_y_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 265);

    auto ta1_y_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 266);

    auto ta1_y_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 267);

    auto ta1_y_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 268);

    auto ta1_y_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 269);

    auto ta1_y_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 271);

    auto ta1_y_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 272);

    auto ta1_y_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 273);

    auto ta1_y_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 275);

    auto ta1_y_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 276);

    auto ta1_y_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 279);

    auto ta1_y_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 280);

    auto ta1_y_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 281);

    auto ta1_y_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 282);

    auto ta1_y_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 283);

    auto ta1_y_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 284);

    auto ta1_y_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 285);

    auto ta1_y_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 286);

    auto ta1_y_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 287);

    auto ta1_y_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 288);

    auto ta1_y_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 289);

    auto ta1_y_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 290);

    auto ta1_y_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 291);

    auto ta1_y_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 292);

    auto ta1_y_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 293);

    auto ta1_y_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 294);

    auto ta1_y_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 295);

    auto ta1_y_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 296);

    auto ta1_y_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 297);

    auto ta1_y_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 298);

    auto ta1_y_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 299);

    auto ta1_z_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 300);

    auto ta1_z_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 301);

    auto ta1_z_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 302);

    auto ta1_z_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 303);

    auto ta1_z_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 304);

    auto ta1_z_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 305);

    auto ta1_z_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 306);

    auto ta1_z_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 307);

    auto ta1_z_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 308);

    auto ta1_z_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 309);

    auto ta1_z_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 310);

    auto ta1_z_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 311);

    auto ta1_z_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 312);

    auto ta1_z_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 313);

    auto ta1_z_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 314);

    auto ta1_z_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 315);

    auto ta1_z_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 317);

    auto ta1_z_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 320);

    auto ta1_z_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 324);

    auto ta1_z_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 325);

    auto ta1_z_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 326);

    auto ta1_z_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 327);

    auto ta1_z_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 328);

    auto ta1_z_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 330);

    auto ta1_z_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 331);

    auto ta1_z_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 332);

    auto ta1_z_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 333);

    auto ta1_z_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 335);

    auto ta1_z_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 336);

    auto ta1_z_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 339);

    auto ta1_z_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 341);

    auto ta1_z_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 342);

    auto ta1_z_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 343);

    auto ta1_z_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 344);

    auto ta1_z_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 346);

    auto ta1_z_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 348);

    auto ta1_z_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 349);

    auto ta1_z_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 351);

    auto ta1_z_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 352);

    auto ta1_z_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 353);

    auto ta1_z_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 355);

    auto ta1_z_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 356);

    auto ta1_z_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 357);

    auto ta1_z_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 358);

    auto ta1_z_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 359);

    auto ta1_z_xyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 371);

    auto ta1_z_xyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 372);

    auto ta1_z_xyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 373);

    auto ta1_z_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 377);

    auto ta1_z_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 379);

    auto ta1_z_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 380);

    auto ta1_z_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 382);

    auto ta1_z_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 383);

    auto ta1_z_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 384);

    auto ta1_z_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 385);

    auto ta1_z_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 386);

    auto ta1_z_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 387);

    auto ta1_z_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 388);

    auto ta1_z_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 389);

    auto ta1_z_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 390);

    auto ta1_z_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 391);

    auto ta1_z_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 392);

    auto ta1_z_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 393);

    auto ta1_z_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 394);

    auto ta1_z_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 395);

    auto ta1_z_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 396);

    auto ta1_z_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 397);

    auto ta1_z_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 398);

    auto ta1_z_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 399);

    auto ta1_z_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 400);

    auto ta1_z_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 401);

    auto ta1_z_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 402);

    auto ta1_z_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 403);

    auto ta1_z_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 404);

    auto ta1_z_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 406);

    auto ta1_z_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 407);

    auto ta1_z_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 408);

    auto ta1_z_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 410);

    auto ta1_z_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 411);

    auto ta1_z_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 414);

    auto ta1_z_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 415);

    auto ta1_z_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 416);

    auto ta1_z_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 417);

    auto ta1_z_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 418);

    auto ta1_z_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 419);

    auto ta1_z_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 420);

    auto ta1_z_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 422);

    auto ta1_z_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 424);

    auto ta1_z_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 425);

    auto ta1_z_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 427);

    auto ta1_z_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 428);

    auto ta1_z_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 429);

    auto ta1_z_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 430);

    auto ta1_z_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 431);

    auto ta1_z_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 432);

    auto ta1_z_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 433);

    auto ta1_z_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 434);

    auto ta1_z_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 435);

    auto ta1_z_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 436);

    auto ta1_z_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 437);

    auto ta1_z_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 438);

    auto ta1_z_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 439);

    auto ta1_z_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 440);

    auto ta1_z_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 441);

    auto ta1_z_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 442);

    auto ta1_z_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 443);

    auto ta1_z_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 444);

    auto ta1_z_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 445);

    auto ta1_z_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 446);

    auto ta1_z_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 447);

    auto ta1_z_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 448);

    auto ta1_z_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 449);

    // Set up components of auxiliary buffer : GF

    auto ta1_x_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf);

    auto ta1_x_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 1);

    auto ta1_x_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 2);

    auto ta1_x_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 3);

    auto ta1_x_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 4);

    auto ta1_x_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 5);

    auto ta1_x_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 6);

    auto ta1_x_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 7);

    auto ta1_x_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 8);

    auto ta1_x_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 9);

    auto ta1_x_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 10);

    auto ta1_x_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 11);

    auto ta1_x_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 12);

    auto ta1_x_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 13);

    auto ta1_x_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 14);

    auto ta1_x_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 15);

    auto ta1_x_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 20);

    auto ta1_x_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 21);

    auto ta1_x_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 22);

    auto ta1_x_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 23);

    auto ta1_x_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 24);

    auto ta1_x_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 25);

    auto ta1_x_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 27);

    auto ta1_x_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 28);

    auto ta1_x_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 29);

    auto ta1_x_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 30);

    auto ta1_x_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 31);

    auto ta1_x_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 32);

    auto ta1_x_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 33);

    auto ta1_x_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 34);

    auto ta1_x_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 35);

    auto ta1_x_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 36);

    auto ta1_x_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 37);

    auto ta1_x_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 38);

    auto ta1_x_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 50);

    auto ta1_x_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 51);

    auto ta1_x_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 52);

    auto ta1_x_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 53);

    auto ta1_x_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 54);

    auto ta1_x_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 55);

    auto ta1_x_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 56);

    auto ta1_x_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 57);

    auto ta1_x_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 58);

    auto ta1_x_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 59);

    auto ta1_x_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 61);

    auto ta1_x_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 63);

    auto ta1_x_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 64);

    auto ta1_x_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 90);

    auto ta1_x_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 91);

    auto ta1_x_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 92);

    auto ta1_x_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 93);

    auto ta1_x_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 94);

    auto ta1_x_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 95);

    auto ta1_x_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 100);

    auto ta1_x_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 101);

    auto ta1_x_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 102);

    auto ta1_x_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 103);

    auto ta1_x_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 104);

    auto ta1_x_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 105);

    auto ta1_x_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 106);

    auto ta1_x_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 107);

    auto ta1_x_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 108);

    auto ta1_x_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 109);

    auto ta1_x_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 122);

    auto ta1_x_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 124);

    auto ta1_x_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 125);

    auto ta1_x_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 127);

    auto ta1_x_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 128);

    auto ta1_x_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 129);

    auto ta1_x_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 132);

    auto ta1_x_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 134);

    auto ta1_x_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 135);

    auto ta1_x_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 137);

    auto ta1_x_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 138);

    auto ta1_x_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 139);

    auto ta1_x_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 140);

    auto ta1_x_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 141);

    auto ta1_x_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 142);

    auto ta1_x_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 143);

    auto ta1_x_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 144);

    auto ta1_x_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 145);

    auto ta1_x_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 146);

    auto ta1_x_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 147);

    auto ta1_x_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 148);

    auto ta1_x_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 149);

    auto ta1_y_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 150);

    auto ta1_y_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 151);

    auto ta1_y_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 152);

    auto ta1_y_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 153);

    auto ta1_y_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 154);

    auto ta1_y_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 155);

    auto ta1_y_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 156);

    auto ta1_y_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 157);

    auto ta1_y_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 158);

    auto ta1_y_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 159);

    auto ta1_y_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 161);

    auto ta1_y_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 163);

    auto ta1_y_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 164);

    auto ta1_y_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 180);

    auto ta1_y_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 181);

    auto ta1_y_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 182);

    auto ta1_y_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 183);

    auto ta1_y_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 184);

    auto ta1_y_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 185);

    auto ta1_y_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 186);

    auto ta1_y_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 187);

    auto ta1_y_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 188);

    auto ta1_y_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 202);

    auto ta1_y_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 204);

    auto ta1_y_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 205);

    auto ta1_y_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 207);

    auto ta1_y_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 208);

    auto ta1_y_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 209);

    auto ta1_y_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 211);

    auto ta1_y_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 213);

    auto ta1_y_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 214);

    auto ta1_y_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 216);

    auto ta1_y_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 217);

    auto ta1_y_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 218);

    auto ta1_y_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 242);

    auto ta1_y_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 244);

    auto ta1_y_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 245);

    auto ta1_y_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 247);

    auto ta1_y_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 248);

    auto ta1_y_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 249);

    auto ta1_y_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 250);

    auto ta1_y_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 251);

    auto ta1_y_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 252);

    auto ta1_y_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 253);

    auto ta1_y_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 254);

    auto ta1_y_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 255);

    auto ta1_y_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 256);

    auto ta1_y_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 257);

    auto ta1_y_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 258);

    auto ta1_y_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 259);

    auto ta1_y_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 261);

    auto ta1_y_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 262);

    auto ta1_y_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 263);

    auto ta1_y_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 264);

    auto ta1_y_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 265);

    auto ta1_y_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 266);

    auto ta1_y_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 267);

    auto ta1_y_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 268);

    auto ta1_y_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 269);

    auto ta1_y_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 270);

    auto ta1_y_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 271);

    auto ta1_y_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 272);

    auto ta1_y_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 273);

    auto ta1_y_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 274);

    auto ta1_y_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 275);

    auto ta1_y_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 276);

    auto ta1_y_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 277);

    auto ta1_y_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 278);

    auto ta1_y_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 279);

    auto ta1_y_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 281);

    auto ta1_y_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 283);

    auto ta1_y_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 284);

    auto ta1_y_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 286);

    auto ta1_y_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 287);

    auto ta1_y_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 288);

    auto ta1_y_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 290);

    auto ta1_y_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 291);

    auto ta1_y_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 292);

    auto ta1_y_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 293);

    auto ta1_y_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 294);

    auto ta1_y_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 295);

    auto ta1_y_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 296);

    auto ta1_y_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 297);

    auto ta1_y_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 298);

    auto ta1_y_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 299);

    auto ta1_z_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 300);

    auto ta1_z_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 301);

    auto ta1_z_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 302);

    auto ta1_z_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 303);

    auto ta1_z_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 304);

    auto ta1_z_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 305);

    auto ta1_z_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 306);

    auto ta1_z_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 307);

    auto ta1_z_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 308);

    auto ta1_z_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 309);

    auto ta1_z_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 322);

    auto ta1_z_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 324);

    auto ta1_z_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 325);

    auto ta1_z_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 331);

    auto ta1_z_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 333);

    auto ta1_z_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 334);

    auto ta1_z_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 336);

    auto ta1_z_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 337);

    auto ta1_z_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 338);

    auto ta1_z_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 350);

    auto ta1_z_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 351);

    auto ta1_z_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 352);

    auto ta1_z_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 353);

    auto ta1_z_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 354);

    auto ta1_z_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 355);

    auto ta1_z_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 357);

    auto ta1_z_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 358);

    auto ta1_z_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 359);

    auto ta1_z_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 361);

    auto ta1_z_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 363);

    auto ta1_z_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 364);

    auto ta1_z_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 366);

    auto ta1_z_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 367);

    auto ta1_z_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 368);

    auto ta1_z_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 392);

    auto ta1_z_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 394);

    auto ta1_z_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 395);

    auto ta1_z_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 397);

    auto ta1_z_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 398);

    auto ta1_z_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 399);

    auto ta1_z_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 400);

    auto ta1_z_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 401);

    auto ta1_z_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 402);

    auto ta1_z_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 403);

    auto ta1_z_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 404);

    auto ta1_z_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 405);

    auto ta1_z_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 406);

    auto ta1_z_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 407);

    auto ta1_z_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 408);

    auto ta1_z_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 409);

    auto ta1_z_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 412);

    auto ta1_z_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 414);

    auto ta1_z_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 415);

    auto ta1_z_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 417);

    auto ta1_z_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 418);

    auto ta1_z_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 419);

    auto ta1_z_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 420);

    auto ta1_z_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 421);

    auto ta1_z_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 422);

    auto ta1_z_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 423);

    auto ta1_z_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 424);

    auto ta1_z_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 425);

    auto ta1_z_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 426);

    auto ta1_z_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 427);

    auto ta1_z_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 428);

    auto ta1_z_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 429);

    auto ta1_z_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 431);

    auto ta1_z_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 432);

    auto ta1_z_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 433);

    auto ta1_z_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 434);

    auto ta1_z_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 435);

    auto ta1_z_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 436);

    auto ta1_z_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 437);

    auto ta1_z_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 438);

    auto ta1_z_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 439);

    auto ta1_z_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 440);

    auto ta1_z_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 441);

    auto ta1_z_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 442);

    auto ta1_z_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 443);

    auto ta1_z_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 444);

    auto ta1_z_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 445);

    auto ta1_z_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 446);

    auto ta1_z_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 447);

    auto ta1_z_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 448);

    auto ta1_z_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 449);

    // Set up components of auxiliary buffer : GF

    auto ta1_x_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf);

    auto ta1_x_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 1);

    auto ta1_x_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 2);

    auto ta1_x_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 3);

    auto ta1_x_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 4);

    auto ta1_x_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 5);

    auto ta1_x_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 6);

    auto ta1_x_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 7);

    auto ta1_x_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 8);

    auto ta1_x_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 9);

    auto ta1_x_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 10);

    auto ta1_x_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 11);

    auto ta1_x_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 12);

    auto ta1_x_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 13);

    auto ta1_x_xxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 14);

    auto ta1_x_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 15);

    auto ta1_x_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 20);

    auto ta1_x_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 21);

    auto ta1_x_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 22);

    auto ta1_x_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 23);

    auto ta1_x_xxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 24);

    auto ta1_x_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 25);

    auto ta1_x_xxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 27);

    auto ta1_x_xxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 28);

    auto ta1_x_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 29);

    auto ta1_x_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 30);

    auto ta1_x_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 31);

    auto ta1_x_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 32);

    auto ta1_x_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 33);

    auto ta1_x_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 34);

    auto ta1_x_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 35);

    auto ta1_x_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 36);

    auto ta1_x_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 37);

    auto ta1_x_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 38);

    auto ta1_x_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 50);

    auto ta1_x_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 51);

    auto ta1_x_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 52);

    auto ta1_x_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 53);

    auto ta1_x_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 54);

    auto ta1_x_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 55);

    auto ta1_x_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 56);

    auto ta1_x_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 57);

    auto ta1_x_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 58);

    auto ta1_x_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 59);

    auto ta1_x_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 61);

    auto ta1_x_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 63);

    auto ta1_x_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 64);

    auto ta1_x_xzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 90);

    auto ta1_x_xzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 91);

    auto ta1_x_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 92);

    auto ta1_x_xzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 93);

    auto ta1_x_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 94);

    auto ta1_x_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 95);

    auto ta1_x_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 100);

    auto ta1_x_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 101);

    auto ta1_x_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 102);

    auto ta1_x_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 103);

    auto ta1_x_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 104);

    auto ta1_x_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 105);

    auto ta1_x_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 106);

    auto ta1_x_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 107);

    auto ta1_x_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 108);

    auto ta1_x_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 109);

    auto ta1_x_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 122);

    auto ta1_x_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 124);

    auto ta1_x_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 125);

    auto ta1_x_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 127);

    auto ta1_x_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 128);

    auto ta1_x_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 129);

    auto ta1_x_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 132);

    auto ta1_x_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 134);

    auto ta1_x_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 135);

    auto ta1_x_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 137);

    auto ta1_x_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 138);

    auto ta1_x_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 139);

    auto ta1_x_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 140);

    auto ta1_x_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 141);

    auto ta1_x_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 142);

    auto ta1_x_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 143);

    auto ta1_x_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 144);

    auto ta1_x_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 145);

    auto ta1_x_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 146);

    auto ta1_x_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 147);

    auto ta1_x_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 148);

    auto ta1_x_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 149);

    auto ta1_y_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 150);

    auto ta1_y_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 151);

    auto ta1_y_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 152);

    auto ta1_y_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 153);

    auto ta1_y_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 154);

    auto ta1_y_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 155);

    auto ta1_y_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 156);

    auto ta1_y_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 157);

    auto ta1_y_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 158);

    auto ta1_y_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 159);

    auto ta1_y_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 161);

    auto ta1_y_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 163);

    auto ta1_y_xxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 164);

    auto ta1_y_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 180);

    auto ta1_y_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 181);

    auto ta1_y_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 182);

    auto ta1_y_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 183);

    auto ta1_y_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 184);

    auto ta1_y_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 185);

    auto ta1_y_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 186);

    auto ta1_y_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 187);

    auto ta1_y_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 188);

    auto ta1_y_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 202);

    auto ta1_y_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 204);

    auto ta1_y_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 205);

    auto ta1_y_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 207);

    auto ta1_y_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 208);

    auto ta1_y_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 209);

    auto ta1_y_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 211);

    auto ta1_y_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 213);

    auto ta1_y_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 214);

    auto ta1_y_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 216);

    auto ta1_y_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 217);

    auto ta1_y_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 218);

    auto ta1_y_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 242);

    auto ta1_y_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 244);

    auto ta1_y_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 245);

    auto ta1_y_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 247);

    auto ta1_y_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 248);

    auto ta1_y_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 249);

    auto ta1_y_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 250);

    auto ta1_y_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 251);

    auto ta1_y_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 252);

    auto ta1_y_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 253);

    auto ta1_y_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 254);

    auto ta1_y_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 255);

    auto ta1_y_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 256);

    auto ta1_y_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 257);

    auto ta1_y_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 258);

    auto ta1_y_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 259);

    auto ta1_y_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 261);

    auto ta1_y_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 262);

    auto ta1_y_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 263);

    auto ta1_y_yyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 264);

    auto ta1_y_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 265);

    auto ta1_y_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 266);

    auto ta1_y_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 267);

    auto ta1_y_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 268);

    auto ta1_y_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 269);

    auto ta1_y_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 270);

    auto ta1_y_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 271);

    auto ta1_y_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 272);

    auto ta1_y_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 273);

    auto ta1_y_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 274);

    auto ta1_y_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 275);

    auto ta1_y_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 276);

    auto ta1_y_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 277);

    auto ta1_y_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 278);

    auto ta1_y_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 279);

    auto ta1_y_yzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 281);

    auto ta1_y_yzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 283);

    auto ta1_y_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 284);

    auto ta1_y_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 286);

    auto ta1_y_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 287);

    auto ta1_y_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 288);

    auto ta1_y_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 290);

    auto ta1_y_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 291);

    auto ta1_y_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 292);

    auto ta1_y_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 293);

    auto ta1_y_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 294);

    auto ta1_y_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 295);

    auto ta1_y_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 296);

    auto ta1_y_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 297);

    auto ta1_y_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 298);

    auto ta1_y_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 299);

    auto ta1_z_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 300);

    auto ta1_z_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 301);

    auto ta1_z_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 302);

    auto ta1_z_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 303);

    auto ta1_z_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 304);

    auto ta1_z_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 305);

    auto ta1_z_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 306);

    auto ta1_z_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 307);

    auto ta1_z_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 308);

    auto ta1_z_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 309);

    auto ta1_z_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 322);

    auto ta1_z_xxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 324);

    auto ta1_z_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 325);

    auto ta1_z_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 331);

    auto ta1_z_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 333);

    auto ta1_z_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 334);

    auto ta1_z_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 336);

    auto ta1_z_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 337);

    auto ta1_z_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 338);

    auto ta1_z_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 350);

    auto ta1_z_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 351);

    auto ta1_z_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 352);

    auto ta1_z_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 353);

    auto ta1_z_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 354);

    auto ta1_z_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 355);

    auto ta1_z_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 357);

    auto ta1_z_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 358);

    auto ta1_z_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 359);

    auto ta1_z_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 361);

    auto ta1_z_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 363);

    auto ta1_z_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 364);

    auto ta1_z_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 366);

    auto ta1_z_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 367);

    auto ta1_z_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 368);

    auto ta1_z_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 392);

    auto ta1_z_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 394);

    auto ta1_z_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 395);

    auto ta1_z_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 397);

    auto ta1_z_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 398);

    auto ta1_z_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 399);

    auto ta1_z_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 400);

    auto ta1_z_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 401);

    auto ta1_z_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 402);

    auto ta1_z_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 403);

    auto ta1_z_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 404);

    auto ta1_z_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 405);

    auto ta1_z_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 406);

    auto ta1_z_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 407);

    auto ta1_z_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 408);

    auto ta1_z_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 409);

    auto ta1_z_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 412);

    auto ta1_z_yyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 414);

    auto ta1_z_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 415);

    auto ta1_z_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 417);

    auto ta1_z_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 418);

    auto ta1_z_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 419);

    auto ta1_z_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 420);

    auto ta1_z_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 421);

    auto ta1_z_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 422);

    auto ta1_z_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 423);

    auto ta1_z_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 424);

    auto ta1_z_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 425);

    auto ta1_z_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 426);

    auto ta1_z_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 427);

    auto ta1_z_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 428);

    auto ta1_z_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 429);

    auto ta1_z_yzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 431);

    auto ta1_z_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 432);

    auto ta1_z_yzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 433);

    auto ta1_z_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 434);

    auto ta1_z_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 435);

    auto ta1_z_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 436);

    auto ta1_z_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 437);

    auto ta1_z_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 438);

    auto ta1_z_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 439);

    auto ta1_z_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 440);

    auto ta1_z_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 441);

    auto ta1_z_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 442);

    auto ta1_z_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 443);

    auto ta1_z_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 444);

    auto ta1_z_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 445);

    auto ta1_z_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 446);

    auto ta1_z_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 447);

    auto ta1_z_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 448);

    auto ta1_z_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 449);

    // Set up components of auxiliary buffer : GG

    auto ta_xxxx_xxxx_1 = pbuffer.data(idx_npot_1_gg);

    auto ta_xxxx_xxxy_1 = pbuffer.data(idx_npot_1_gg + 1);

    auto ta_xxxx_xxxz_1 = pbuffer.data(idx_npot_1_gg + 2);

    auto ta_xxxx_xxyy_1 = pbuffer.data(idx_npot_1_gg + 3);

    auto ta_xxxx_xxyz_1 = pbuffer.data(idx_npot_1_gg + 4);

    auto ta_xxxx_xxzz_1 = pbuffer.data(idx_npot_1_gg + 5);

    auto ta_xxxx_xyyy_1 = pbuffer.data(idx_npot_1_gg + 6);

    auto ta_xxxx_xyyz_1 = pbuffer.data(idx_npot_1_gg + 7);

    auto ta_xxxx_xyzz_1 = pbuffer.data(idx_npot_1_gg + 8);

    auto ta_xxxx_xzzz_1 = pbuffer.data(idx_npot_1_gg + 9);

    auto ta_xxxx_yyyy_1 = pbuffer.data(idx_npot_1_gg + 10);

    auto ta_xxxx_yyyz_1 = pbuffer.data(idx_npot_1_gg + 11);

    auto ta_xxxx_yyzz_1 = pbuffer.data(idx_npot_1_gg + 12);

    auto ta_xxxx_yzzz_1 = pbuffer.data(idx_npot_1_gg + 13);

    auto ta_xxxx_zzzz_1 = pbuffer.data(idx_npot_1_gg + 14);

    auto ta_xxxy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 15);

    auto ta_xxxy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 16);

    auto ta_xxxy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 17);

    auto ta_xxxy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 18);

    auto ta_xxxy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 20);

    auto ta_xxxy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 21);

    auto ta_xxxy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 24);

    auto ta_xxxy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 25);

    auto ta_xxxz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 30);

    auto ta_xxxz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 31);

    auto ta_xxxz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 32);

    auto ta_xxxz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 33);

    auto ta_xxxz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 35);

    auto ta_xxxz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 36);

    auto ta_xxxz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 39);

    auto ta_xxxz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 44);

    auto ta_xxyy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 45);

    auto ta_xxyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 46);

    auto ta_xxyy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 47);

    auto ta_xxyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 48);

    auto ta_xxyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 49);

    auto ta_xxyy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 50);

    auto ta_xxyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 51);

    auto ta_xxyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 52);

    auto ta_xxyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 53);

    auto ta_xxyy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 54);

    auto ta_xxyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 55);

    auto ta_xxyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 56);

    auto ta_xxyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 57);

    auto ta_xxyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 58);

    auto ta_xxzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 75);

    auto ta_xxzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 76);

    auto ta_xxzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 77);

    auto ta_xxzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 78);

    auto ta_xxzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 79);

    auto ta_xxzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 80);

    auto ta_xxzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 81);

    auto ta_xxzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 82);

    auto ta_xxzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 83);

    auto ta_xxzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 84);

    auto ta_xxzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 86);

    auto ta_xxzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 87);

    auto ta_xxzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 88);

    auto ta_xxzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 89);

    auto ta_xyyy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 90);

    auto ta_xyyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 91);

    auto ta_xyyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 93);

    auto ta_xyyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 96);

    auto ta_xyyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 100);

    auto ta_xyyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 101);

    auto ta_xyyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 102);

    auto ta_xyyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 103);

    auto ta_xzzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 135);

    auto ta_xzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 137);

    auto ta_xzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 140);

    auto ta_xzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 144);

    auto ta_xzzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 146);

    auto ta_xzzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 147);

    auto ta_xzzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 148);

    auto ta_xzzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 149);

    auto ta_yyyy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 150);

    auto ta_yyyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 151);

    auto ta_yyyy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 152);

    auto ta_yyyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 153);

    auto ta_yyyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 154);

    auto ta_yyyy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 155);

    auto ta_yyyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 156);

    auto ta_yyyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 157);

    auto ta_yyyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 158);

    auto ta_yyyy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 159);

    auto ta_yyyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 160);

    auto ta_yyyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 161);

    auto ta_yyyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 162);

    auto ta_yyyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 163);

    auto ta_yyyy_zzzz_1 = pbuffer.data(idx_npot_1_gg + 164);

    auto ta_yyyz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 166);

    auto ta_yyyz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 168);

    auto ta_yyyz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 171);

    auto ta_yyyz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 175);

    auto ta_yyyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 176);

    auto ta_yyyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 177);

    auto ta_yyyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 178);

    auto ta_yyyz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 179);

    auto ta_yyzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 181);

    auto ta_yyzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 182);

    auto ta_yyzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 183);

    auto ta_yyzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 184);

    auto ta_yyzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 185);

    auto ta_yyzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 186);

    auto ta_yyzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 187);

    auto ta_yyzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 188);

    auto ta_yyzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 189);

    auto ta_yyzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 190);

    auto ta_yyzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 191);

    auto ta_yyzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 192);

    auto ta_yyzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 193);

    auto ta_yyzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 194);

    auto ta_yzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 197);

    auto ta_yzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 200);

    auto ta_yzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 204);

    auto ta_yzzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 205);

    auto ta_yzzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 206);

    auto ta_yzzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 207);

    auto ta_yzzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 208);

    auto ta_yzzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 209);

    auto ta_zzzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 210);

    auto ta_zzzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 211);

    auto ta_zzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 212);

    auto ta_zzzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 213);

    auto ta_zzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 214);

    auto ta_zzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 215);

    auto ta_zzzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 216);

    auto ta_zzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 217);

    auto ta_zzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 218);

    auto ta_zzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 219);

    auto ta_zzzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 220);

    auto ta_zzzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 221);

    auto ta_zzzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 222);

    auto ta_zzzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 223);

    auto ta_zzzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 224);

    // Set up components of auxiliary buffer : GG

    auto ta1_x_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg);

    auto ta1_x_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 1);

    auto ta1_x_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 2);

    auto ta1_x_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 3);

    auto ta1_x_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 4);

    auto ta1_x_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 5);

    auto ta1_x_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 6);

    auto ta1_x_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 7);

    auto ta1_x_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 8);

    auto ta1_x_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 9);

    auto ta1_x_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 10);

    auto ta1_x_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 11);

    auto ta1_x_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 12);

    auto ta1_x_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 13);

    auto ta1_x_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 14);

    auto ta1_x_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 15);

    auto ta1_x_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 16);

    auto ta1_x_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 17);

    auto ta1_x_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 18);

    auto ta1_x_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 19);

    auto ta1_x_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 20);

    auto ta1_x_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 21);

    auto ta1_x_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 22);

    auto ta1_x_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 23);

    auto ta1_x_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 24);

    auto ta1_x_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 25);

    auto ta1_x_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 29);

    auto ta1_x_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 30);

    auto ta1_x_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 31);

    auto ta1_x_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 32);

    auto ta1_x_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 33);

    auto ta1_x_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 34);

    auto ta1_x_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 35);

    auto ta1_x_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 36);

    auto ta1_x_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 37);

    auto ta1_x_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 38);

    auto ta1_x_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 39);

    auto ta1_x_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 40);

    auto ta1_x_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 41);

    auto ta1_x_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 42);

    auto ta1_x_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 43);

    auto ta1_x_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 44);

    auto ta1_x_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 45);

    auto ta1_x_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 46);

    auto ta1_x_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 47);

    auto ta1_x_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 48);

    auto ta1_x_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 49);

    auto ta1_x_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 50);

    auto ta1_x_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 51);

    auto ta1_x_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 52);

    auto ta1_x_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 53);

    auto ta1_x_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 54);

    auto ta1_x_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 55);

    auto ta1_x_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 56);

    auto ta1_x_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 57);

    auto ta1_x_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 58);

    auto ta1_x_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 59);

    auto ta1_x_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 62);

    auto ta1_x_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 65);

    auto ta1_x_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 69);

    auto ta1_x_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 74);

    auto ta1_x_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 75);

    auto ta1_x_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 76);

    auto ta1_x_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 77);

    auto ta1_x_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 78);

    auto ta1_x_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 79);

    auto ta1_x_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 80);

    auto ta1_x_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 81);

    auto ta1_x_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 82);

    auto ta1_x_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 83);

    auto ta1_x_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 84);

    auto ta1_x_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 85);

    auto ta1_x_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 86);

    auto ta1_x_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 87);

    auto ta1_x_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 88);

    auto ta1_x_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 89);

    auto ta1_x_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 90);

    auto ta1_x_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 91);

    auto ta1_x_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 92);

    auto ta1_x_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 93);

    auto ta1_x_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 94);

    auto ta1_x_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 95);

    auto ta1_x_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 96);

    auto ta1_x_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 97);

    auto ta1_x_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 98);

    auto ta1_x_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 99);

    auto ta1_x_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 100);

    auto ta1_x_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 101);

    auto ta1_x_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 102);

    auto ta1_x_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 103);

    auto ta1_x_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 106);

    auto ta1_x_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 107);

    auto ta1_x_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 108);

    auto ta1_x_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 110);

    auto ta1_x_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 111);

    auto ta1_x_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 114);

    auto ta1_x_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 120);

    auto ta1_x_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 122);

    auto ta1_x_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 125);

    auto ta1_x_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 129);

    auto ta1_x_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 135);

    auto ta1_x_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 136);

    auto ta1_x_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 137);

    auto ta1_x_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 138);

    auto ta1_x_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 139);

    auto ta1_x_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 140);

    auto ta1_x_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 141);

    auto ta1_x_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 142);

    auto ta1_x_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 143);

    auto ta1_x_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 144);

    auto ta1_x_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 146);

    auto ta1_x_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 147);

    auto ta1_x_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 148);

    auto ta1_x_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 149);

    auto ta1_x_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 150);

    auto ta1_x_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 151);

    auto ta1_x_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 152);

    auto ta1_x_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 153);

    auto ta1_x_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 154);

    auto ta1_x_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 155);

    auto ta1_x_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 156);

    auto ta1_x_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 157);

    auto ta1_x_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 158);

    auto ta1_x_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 159);

    auto ta1_x_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 160);

    auto ta1_x_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 161);

    auto ta1_x_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 162);

    auto ta1_x_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 163);

    auto ta1_x_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 164);

    auto ta1_x_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 166);

    auto ta1_x_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 167);

    auto ta1_x_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 168);

    auto ta1_x_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 170);

    auto ta1_x_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 171);

    auto ta1_x_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 174);

    auto ta1_x_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 175);

    auto ta1_x_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 176);

    auto ta1_x_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 177);

    auto ta1_x_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 178);

    auto ta1_x_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 179);

    auto ta1_x_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 180);

    auto ta1_x_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 181);

    auto ta1_x_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 182);

    auto ta1_x_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 183);

    auto ta1_x_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 184);

    auto ta1_x_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 185);

    auto ta1_x_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 186);

    auto ta1_x_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 187);

    auto ta1_x_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 188);

    auto ta1_x_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 189);

    auto ta1_x_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 190);

    auto ta1_x_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 191);

    auto ta1_x_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 192);

    auto ta1_x_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 193);

    auto ta1_x_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 194);

    auto ta1_x_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 195);

    auto ta1_x_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 197);

    auto ta1_x_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 199);

    auto ta1_x_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 200);

    auto ta1_x_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 202);

    auto ta1_x_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 203);

    auto ta1_x_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 204);

    auto ta1_x_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 205);

    auto ta1_x_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 206);

    auto ta1_x_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 207);

    auto ta1_x_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 208);

    auto ta1_x_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 209);

    auto ta1_x_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 210);

    auto ta1_x_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 211);

    auto ta1_x_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 212);

    auto ta1_x_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 213);

    auto ta1_x_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 214);

    auto ta1_x_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 215);

    auto ta1_x_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 216);

    auto ta1_x_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 217);

    auto ta1_x_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 218);

    auto ta1_x_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 219);

    auto ta1_x_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 220);

    auto ta1_x_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 221);

    auto ta1_x_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 222);

    auto ta1_x_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 223);

    auto ta1_x_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 224);

    auto ta1_y_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 225);

    auto ta1_y_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 226);

    auto ta1_y_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 227);

    auto ta1_y_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 228);

    auto ta1_y_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 229);

    auto ta1_y_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 230);

    auto ta1_y_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 231);

    auto ta1_y_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 232);

    auto ta1_y_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 233);

    auto ta1_y_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 234);

    auto ta1_y_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 235);

    auto ta1_y_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 236);

    auto ta1_y_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 237);

    auto ta1_y_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 238);

    auto ta1_y_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 239);

    auto ta1_y_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 240);

    auto ta1_y_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 241);

    auto ta1_y_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 242);

    auto ta1_y_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 243);

    auto ta1_y_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 244);

    auto ta1_y_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 245);

    auto ta1_y_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 246);

    auto ta1_y_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 247);

    auto ta1_y_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 248);

    auto ta1_y_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 249);

    auto ta1_y_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 250);

    auto ta1_y_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 251);

    auto ta1_y_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 252);

    auto ta1_y_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 253);

    auto ta1_y_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 255);

    auto ta1_y_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 256);

    auto ta1_y_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 257);

    auto ta1_y_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 258);

    auto ta1_y_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 260);

    auto ta1_y_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 261);

    auto ta1_y_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 264);

    auto ta1_y_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 266);

    auto ta1_y_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 267);

    auto ta1_y_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 268);

    auto ta1_y_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 269);

    auto ta1_y_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 270);

    auto ta1_y_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 271);

    auto ta1_y_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 272);

    auto ta1_y_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 273);

    auto ta1_y_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 274);

    auto ta1_y_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 275);

    auto ta1_y_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 276);

    auto ta1_y_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 277);

    auto ta1_y_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 278);

    auto ta1_y_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 279);

    auto ta1_y_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 280);

    auto ta1_y_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 281);

    auto ta1_y_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 282);

    auto ta1_y_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 283);

    auto ta1_y_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 284);

    auto ta1_y_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 286);

    auto ta1_y_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 288);

    auto ta1_y_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 291);

    auto ta1_y_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 296);

    auto ta1_y_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 297);

    auto ta1_y_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 298);

    auto ta1_y_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 300);

    auto ta1_y_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 301);

    auto ta1_y_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 302);

    auto ta1_y_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 303);

    auto ta1_y_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 304);

    auto ta1_y_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 305);

    auto ta1_y_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 306);

    auto ta1_y_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 307);

    auto ta1_y_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 308);

    auto ta1_y_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 309);

    auto ta1_y_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 310);

    auto ta1_y_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 311);

    auto ta1_y_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 312);

    auto ta1_y_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 313);

    auto ta1_y_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 314);

    auto ta1_y_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 315);

    auto ta1_y_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 316);

    auto ta1_y_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 318);

    auto ta1_y_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 319);

    auto ta1_y_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 321);

    auto ta1_y_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 322);

    auto ta1_y_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 323);

    auto ta1_y_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 325);

    auto ta1_y_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 326);

    auto ta1_y_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 327);

    auto ta1_y_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 328);

    auto ta1_y_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 329);

    auto ta1_y_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 341);

    auto ta1_y_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 342);

    auto ta1_y_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 343);

    auto ta1_y_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 344);

    auto ta1_y_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 355);

    auto ta1_y_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 356);

    auto ta1_y_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 357);

    auto ta1_y_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 358);

    auto ta1_y_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 360);

    auto ta1_y_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 362);

    auto ta1_y_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 364);

    auto ta1_y_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 365);

    auto ta1_y_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 367);

    auto ta1_y_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 368);

    auto ta1_y_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 369);

    auto ta1_y_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 370);

    auto ta1_y_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 371);

    auto ta1_y_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 372);

    auto ta1_y_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 373);

    auto ta1_y_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 374);

    auto ta1_y_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 375);

    auto ta1_y_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 376);

    auto ta1_y_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 377);

    auto ta1_y_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 378);

    auto ta1_y_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 379);

    auto ta1_y_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 380);

    auto ta1_y_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 381);

    auto ta1_y_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 382);

    auto ta1_y_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 383);

    auto ta1_y_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 384);

    auto ta1_y_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 385);

    auto ta1_y_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 386);

    auto ta1_y_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 387);

    auto ta1_y_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 388);

    auto ta1_y_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 389);

    auto ta1_y_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 390);

    auto ta1_y_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 391);

    auto ta1_y_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 392);

    auto ta1_y_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 393);

    auto ta1_y_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 394);

    auto ta1_y_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 395);

    auto ta1_y_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 396);

    auto ta1_y_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 397);

    auto ta1_y_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 398);

    auto ta1_y_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 399);

    auto ta1_y_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 400);

    auto ta1_y_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 401);

    auto ta1_y_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 402);

    auto ta1_y_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 403);

    auto ta1_y_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 404);

    auto ta1_y_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 405);

    auto ta1_y_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 406);

    auto ta1_y_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 407);

    auto ta1_y_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 408);

    auto ta1_y_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 409);

    auto ta1_y_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 410);

    auto ta1_y_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 411);

    auto ta1_y_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 412);

    auto ta1_y_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 413);

    auto ta1_y_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 414);

    auto ta1_y_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 415);

    auto ta1_y_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 416);

    auto ta1_y_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 417);

    auto ta1_y_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 418);

    auto ta1_y_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 419);

    auto ta1_y_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 421);

    auto ta1_y_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 422);

    auto ta1_y_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 423);

    auto ta1_y_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 424);

    auto ta1_y_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 425);

    auto ta1_y_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 426);

    auto ta1_y_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 427);

    auto ta1_y_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 428);

    auto ta1_y_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 429);

    auto ta1_y_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 430);

    auto ta1_y_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 431);

    auto ta1_y_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 432);

    auto ta1_y_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 433);

    auto ta1_y_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 434);

    auto ta1_y_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 435);

    auto ta1_y_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 436);

    auto ta1_y_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 437);

    auto ta1_y_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 438);

    auto ta1_y_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 439);

    auto ta1_y_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 440);

    auto ta1_y_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 441);

    auto ta1_y_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 442);

    auto ta1_y_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 443);

    auto ta1_y_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 444);

    auto ta1_y_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 445);

    auto ta1_y_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 446);

    auto ta1_y_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 447);

    auto ta1_y_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 448);

    auto ta1_y_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 449);

    auto ta1_z_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 450);

    auto ta1_z_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 451);

    auto ta1_z_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 452);

    auto ta1_z_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 453);

    auto ta1_z_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 454);

    auto ta1_z_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 455);

    auto ta1_z_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 456);

    auto ta1_z_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 457);

    auto ta1_z_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 458);

    auto ta1_z_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 459);

    auto ta1_z_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 460);

    auto ta1_z_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 461);

    auto ta1_z_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 462);

    auto ta1_z_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 463);

    auto ta1_z_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 464);

    auto ta1_z_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 465);

    auto ta1_z_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 466);

    auto ta1_z_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 467);

    auto ta1_z_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 468);

    auto ta1_z_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 470);

    auto ta1_z_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 471);

    auto ta1_z_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 474);

    auto ta1_z_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 475);

    auto ta1_z_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 476);

    auto ta1_z_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 477);

    auto ta1_z_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 478);

    auto ta1_z_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 480);

    auto ta1_z_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 481);

    auto ta1_z_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 482);

    auto ta1_z_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 483);

    auto ta1_z_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 484);

    auto ta1_z_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 485);

    auto ta1_z_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 486);

    auto ta1_z_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 487);

    auto ta1_z_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 488);

    auto ta1_z_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 489);

    auto ta1_z_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 491);

    auto ta1_z_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 492);

    auto ta1_z_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 493);

    auto ta1_z_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 494);

    auto ta1_z_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 495);

    auto ta1_z_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 496);

    auto ta1_z_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 497);

    auto ta1_z_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 498);

    auto ta1_z_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 499);

    auto ta1_z_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 500);

    auto ta1_z_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 501);

    auto ta1_z_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 502);

    auto ta1_z_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 503);

    auto ta1_z_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 504);

    auto ta1_z_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 505);

    auto ta1_z_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 506);

    auto ta1_z_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 507);

    auto ta1_z_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 508);

    auto ta1_z_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 509);

    auto ta1_z_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 512);

    auto ta1_z_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 515);

    auto ta1_z_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 519);

    auto ta1_z_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 521);

    auto ta1_z_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 522);

    auto ta1_z_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 523);

    auto ta1_z_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 525);

    auto ta1_z_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 526);

    auto ta1_z_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 527);

    auto ta1_z_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 528);

    auto ta1_z_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 529);

    auto ta1_z_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 530);

    auto ta1_z_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 531);

    auto ta1_z_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 532);

    auto ta1_z_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 533);

    auto ta1_z_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 534);

    auto ta1_z_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 535);

    auto ta1_z_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 536);

    auto ta1_z_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 537);

    auto ta1_z_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 538);

    auto ta1_z_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 539);

    auto ta1_z_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 540);

    auto ta1_z_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 541);

    auto ta1_z_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 543);

    auto ta1_z_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 544);

    auto ta1_z_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 546);

    auto ta1_z_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 547);

    auto ta1_z_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 548);

    auto ta1_z_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 550);

    auto ta1_z_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 551);

    auto ta1_z_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 552);

    auto ta1_z_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 553);

    auto ta1_z_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 554);

    auto ta1_z_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 566);

    auto ta1_z_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 567);

    auto ta1_z_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 568);

    auto ta1_z_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 569);

    auto ta1_z_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 580);

    auto ta1_z_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 581);

    auto ta1_z_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 582);

    auto ta1_z_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 583);

    auto ta1_z_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 585);

    auto ta1_z_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 587);

    auto ta1_z_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 589);

    auto ta1_z_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 590);

    auto ta1_z_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 592);

    auto ta1_z_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 593);

    auto ta1_z_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 594);

    auto ta1_z_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 595);

    auto ta1_z_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 596);

    auto ta1_z_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 597);

    auto ta1_z_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 598);

    auto ta1_z_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 599);

    auto ta1_z_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 600);

    auto ta1_z_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 601);

    auto ta1_z_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 602);

    auto ta1_z_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 603);

    auto ta1_z_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 604);

    auto ta1_z_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 605);

    auto ta1_z_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 606);

    auto ta1_z_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 607);

    auto ta1_z_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 608);

    auto ta1_z_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 609);

    auto ta1_z_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 610);

    auto ta1_z_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 611);

    auto ta1_z_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 612);

    auto ta1_z_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 613);

    auto ta1_z_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 614);

    auto ta1_z_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 616);

    auto ta1_z_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 617);

    auto ta1_z_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 618);

    auto ta1_z_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 619);

    auto ta1_z_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 620);

    auto ta1_z_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 621);

    auto ta1_z_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 622);

    auto ta1_z_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 623);

    auto ta1_z_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 624);

    auto ta1_z_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 625);

    auto ta1_z_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 626);

    auto ta1_z_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 627);

    auto ta1_z_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 628);

    auto ta1_z_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 629);

    auto ta1_z_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 630);

    auto ta1_z_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 631);

    auto ta1_z_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 632);

    auto ta1_z_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 633);

    auto ta1_z_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 634);

    auto ta1_z_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 635);

    auto ta1_z_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 636);

    auto ta1_z_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 637);

    auto ta1_z_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 638);

    auto ta1_z_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 639);

    auto ta1_z_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 640);

    auto ta1_z_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 641);

    auto ta1_z_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 642);

    auto ta1_z_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 643);

    auto ta1_z_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 644);

    auto ta1_z_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 645);

    auto ta1_z_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 646);

    auto ta1_z_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 647);

    auto ta1_z_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 648);

    auto ta1_z_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 649);

    auto ta1_z_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 650);

    auto ta1_z_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 651);

    auto ta1_z_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 652);

    auto ta1_z_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 653);

    auto ta1_z_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 654);

    auto ta1_z_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 655);

    auto ta1_z_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 656);

    auto ta1_z_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 657);

    auto ta1_z_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 658);

    auto ta1_z_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 659);

    auto ta1_z_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 660);

    auto ta1_z_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 661);

    auto ta1_z_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 662);

    auto ta1_z_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 663);

    auto ta1_z_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 664);

    auto ta1_z_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 665);

    auto ta1_z_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 666);

    auto ta1_z_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 667);

    auto ta1_z_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 668);

    auto ta1_z_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 669);

    auto ta1_z_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 670);

    auto ta1_z_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 671);

    auto ta1_z_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 672);

    auto ta1_z_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 673);

    auto ta1_z_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 674);

    // Set up components of auxiliary buffer : GG

    auto ta1_x_xxxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg);

    auto ta1_x_xxxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 1);

    auto ta1_x_xxxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 2);

    auto ta1_x_xxxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 3);

    auto ta1_x_xxxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 4);

    auto ta1_x_xxxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 5);

    auto ta1_x_xxxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 6);

    auto ta1_x_xxxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 7);

    auto ta1_x_xxxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 8);

    auto ta1_x_xxxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 9);

    auto ta1_x_xxxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 10);

    auto ta1_x_xxxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 11);

    auto ta1_x_xxxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 12);

    auto ta1_x_xxxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 13);

    auto ta1_x_xxxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 14);

    auto ta1_x_xxxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 15);

    auto ta1_x_xxxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 16);

    auto ta1_x_xxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 17);

    auto ta1_x_xxxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 18);

    auto ta1_x_xxxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 19);

    auto ta1_x_xxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 20);

    auto ta1_x_xxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 21);

    auto ta1_x_xxxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 22);

    auto ta1_x_xxxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 23);

    auto ta1_x_xxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 24);

    auto ta1_x_xxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 25);

    auto ta1_x_xxxy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 29);

    auto ta1_x_xxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 30);

    auto ta1_x_xxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 31);

    auto ta1_x_xxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 32);

    auto ta1_x_xxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 33);

    auto ta1_x_xxxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 34);

    auto ta1_x_xxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 35);

    auto ta1_x_xxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 36);

    auto ta1_x_xxxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 37);

    auto ta1_x_xxxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 38);

    auto ta1_x_xxxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 39);

    auto ta1_x_xxxz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 40);

    auto ta1_x_xxxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 41);

    auto ta1_x_xxxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 42);

    auto ta1_x_xxxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 43);

    auto ta1_x_xxxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 44);

    auto ta1_x_xxyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 45);

    auto ta1_x_xxyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 46);

    auto ta1_x_xxyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 47);

    auto ta1_x_xxyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 48);

    auto ta1_x_xxyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 49);

    auto ta1_x_xxyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 50);

    auto ta1_x_xxyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 51);

    auto ta1_x_xxyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 52);

    auto ta1_x_xxyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 53);

    auto ta1_x_xxyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 54);

    auto ta1_x_xxyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 55);

    auto ta1_x_xxyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 56);

    auto ta1_x_xxyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 57);

    auto ta1_x_xxyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 58);

    auto ta1_x_xxyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 59);

    auto ta1_x_xxyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 62);

    auto ta1_x_xxyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 65);

    auto ta1_x_xxyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 69);

    auto ta1_x_xxyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 74);

    auto ta1_x_xxzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 75);

    auto ta1_x_xxzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 76);

    auto ta1_x_xxzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 77);

    auto ta1_x_xxzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 78);

    auto ta1_x_xxzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 79);

    auto ta1_x_xxzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 80);

    auto ta1_x_xxzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 81);

    auto ta1_x_xxzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 82);

    auto ta1_x_xxzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 83);

    auto ta1_x_xxzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 84);

    auto ta1_x_xxzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 85);

    auto ta1_x_xxzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 86);

    auto ta1_x_xxzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 87);

    auto ta1_x_xxzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 88);

    auto ta1_x_xxzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 89);

    auto ta1_x_xyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 90);

    auto ta1_x_xyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 91);

    auto ta1_x_xyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 92);

    auto ta1_x_xyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 93);

    auto ta1_x_xyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 94);

    auto ta1_x_xyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 95);

    auto ta1_x_xyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 96);

    auto ta1_x_xyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 97);

    auto ta1_x_xyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 98);

    auto ta1_x_xyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 99);

    auto ta1_x_xyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 100);

    auto ta1_x_xyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 101);

    auto ta1_x_xyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 102);

    auto ta1_x_xyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 103);

    auto ta1_x_xyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 106);

    auto ta1_x_xyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 107);

    auto ta1_x_xyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 108);

    auto ta1_x_xyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 110);

    auto ta1_x_xyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 111);

    auto ta1_x_xyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 114);

    auto ta1_x_xyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 120);

    auto ta1_x_xyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 122);

    auto ta1_x_xyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 125);

    auto ta1_x_xyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 129);

    auto ta1_x_xzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 135);

    auto ta1_x_xzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 136);

    auto ta1_x_xzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 137);

    auto ta1_x_xzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 138);

    auto ta1_x_xzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 139);

    auto ta1_x_xzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 140);

    auto ta1_x_xzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 141);

    auto ta1_x_xzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 142);

    auto ta1_x_xzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 143);

    auto ta1_x_xzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 144);

    auto ta1_x_xzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 146);

    auto ta1_x_xzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 147);

    auto ta1_x_xzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 148);

    auto ta1_x_xzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 149);

    auto ta1_x_yyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 150);

    auto ta1_x_yyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 151);

    auto ta1_x_yyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 152);

    auto ta1_x_yyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 153);

    auto ta1_x_yyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 154);

    auto ta1_x_yyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 155);

    auto ta1_x_yyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 156);

    auto ta1_x_yyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 157);

    auto ta1_x_yyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 158);

    auto ta1_x_yyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 159);

    auto ta1_x_yyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 160);

    auto ta1_x_yyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 161);

    auto ta1_x_yyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 162);

    auto ta1_x_yyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 163);

    auto ta1_x_yyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 164);

    auto ta1_x_yyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 166);

    auto ta1_x_yyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 167);

    auto ta1_x_yyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 168);

    auto ta1_x_yyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 170);

    auto ta1_x_yyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 171);

    auto ta1_x_yyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 174);

    auto ta1_x_yyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 175);

    auto ta1_x_yyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 176);

    auto ta1_x_yyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 177);

    auto ta1_x_yyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 178);

    auto ta1_x_yyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 179);

    auto ta1_x_yyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 180);

    auto ta1_x_yyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 181);

    auto ta1_x_yyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 182);

    auto ta1_x_yyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 183);

    auto ta1_x_yyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 184);

    auto ta1_x_yyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 185);

    auto ta1_x_yyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 186);

    auto ta1_x_yyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 187);

    auto ta1_x_yyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 188);

    auto ta1_x_yyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 189);

    auto ta1_x_yyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 190);

    auto ta1_x_yyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 191);

    auto ta1_x_yyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 192);

    auto ta1_x_yyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 193);

    auto ta1_x_yyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 194);

    auto ta1_x_yzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 195);

    auto ta1_x_yzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 197);

    auto ta1_x_yzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 199);

    auto ta1_x_yzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 200);

    auto ta1_x_yzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 202);

    auto ta1_x_yzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 203);

    auto ta1_x_yzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 204);

    auto ta1_x_yzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 205);

    auto ta1_x_yzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 206);

    auto ta1_x_yzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 207);

    auto ta1_x_yzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 208);

    auto ta1_x_yzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 209);

    auto ta1_x_zzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 210);

    auto ta1_x_zzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 211);

    auto ta1_x_zzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 212);

    auto ta1_x_zzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 213);

    auto ta1_x_zzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 214);

    auto ta1_x_zzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 215);

    auto ta1_x_zzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 216);

    auto ta1_x_zzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 217);

    auto ta1_x_zzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 218);

    auto ta1_x_zzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 219);

    auto ta1_x_zzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 220);

    auto ta1_x_zzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 221);

    auto ta1_x_zzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 222);

    auto ta1_x_zzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 223);

    auto ta1_x_zzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 224);

    auto ta1_y_xxxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 225);

    auto ta1_y_xxxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 226);

    auto ta1_y_xxxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 227);

    auto ta1_y_xxxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 228);

    auto ta1_y_xxxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 229);

    auto ta1_y_xxxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 230);

    auto ta1_y_xxxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 231);

    auto ta1_y_xxxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 232);

    auto ta1_y_xxxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 233);

    auto ta1_y_xxxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 234);

    auto ta1_y_xxxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 235);

    auto ta1_y_xxxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 236);

    auto ta1_y_xxxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 237);

    auto ta1_y_xxxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 238);

    auto ta1_y_xxxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 239);

    auto ta1_y_xxxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 240);

    auto ta1_y_xxxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 241);

    auto ta1_y_xxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 242);

    auto ta1_y_xxxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 243);

    auto ta1_y_xxxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 244);

    auto ta1_y_xxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 245);

    auto ta1_y_xxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 246);

    auto ta1_y_xxxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 247);

    auto ta1_y_xxxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 248);

    auto ta1_y_xxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 249);

    auto ta1_y_xxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 250);

    auto ta1_y_xxxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 251);

    auto ta1_y_xxxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 252);

    auto ta1_y_xxxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 253);

    auto ta1_y_xxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 255);

    auto ta1_y_xxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 256);

    auto ta1_y_xxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 257);

    auto ta1_y_xxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 258);

    auto ta1_y_xxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 260);

    auto ta1_y_xxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 261);

    auto ta1_y_xxxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 264);

    auto ta1_y_xxxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 266);

    auto ta1_y_xxxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 267);

    auto ta1_y_xxxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 268);

    auto ta1_y_xxxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 269);

    auto ta1_y_xxyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 270);

    auto ta1_y_xxyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 271);

    auto ta1_y_xxyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 272);

    auto ta1_y_xxyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 273);

    auto ta1_y_xxyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 274);

    auto ta1_y_xxyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 275);

    auto ta1_y_xxyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 276);

    auto ta1_y_xxyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 277);

    auto ta1_y_xxyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 278);

    auto ta1_y_xxyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 279);

    auto ta1_y_xxyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 280);

    auto ta1_y_xxyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 281);

    auto ta1_y_xxyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 282);

    auto ta1_y_xxyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 283);

    auto ta1_y_xxyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 284);

    auto ta1_y_xxyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 286);

    auto ta1_y_xxyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 288);

    auto ta1_y_xxyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 291);

    auto ta1_y_xxyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 296);

    auto ta1_y_xxyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 297);

    auto ta1_y_xxyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 298);

    auto ta1_y_xxzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 300);

    auto ta1_y_xxzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 301);

    auto ta1_y_xxzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 302);

    auto ta1_y_xxzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 303);

    auto ta1_y_xxzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 304);

    auto ta1_y_xxzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 305);

    auto ta1_y_xxzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 306);

    auto ta1_y_xxzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 307);

    auto ta1_y_xxzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 308);

    auto ta1_y_xxzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 309);

    auto ta1_y_xxzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 310);

    auto ta1_y_xxzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 311);

    auto ta1_y_xxzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 312);

    auto ta1_y_xxzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 313);

    auto ta1_y_xxzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 314);

    auto ta1_y_xyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 315);

    auto ta1_y_xyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 316);

    auto ta1_y_xyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 318);

    auto ta1_y_xyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 319);

    auto ta1_y_xyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 321);

    auto ta1_y_xyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 322);

    auto ta1_y_xyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 323);

    auto ta1_y_xyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 325);

    auto ta1_y_xyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 326);

    auto ta1_y_xyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 327);

    auto ta1_y_xyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 328);

    auto ta1_y_xyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 329);

    auto ta1_y_xyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 341);

    auto ta1_y_xyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 342);

    auto ta1_y_xyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 343);

    auto ta1_y_xyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 344);

    auto ta1_y_xyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 355);

    auto ta1_y_xyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 356);

    auto ta1_y_xyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 357);

    auto ta1_y_xyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 358);

    auto ta1_y_xzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 360);

    auto ta1_y_xzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 362);

    auto ta1_y_xzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 364);

    auto ta1_y_xzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 365);

    auto ta1_y_xzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 367);

    auto ta1_y_xzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 368);

    auto ta1_y_xzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 369);

    auto ta1_y_xzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 370);

    auto ta1_y_xzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 371);

    auto ta1_y_xzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 372);

    auto ta1_y_xzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 373);

    auto ta1_y_xzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 374);

    auto ta1_y_yyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 375);

    auto ta1_y_yyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 376);

    auto ta1_y_yyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 377);

    auto ta1_y_yyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 378);

    auto ta1_y_yyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 379);

    auto ta1_y_yyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 380);

    auto ta1_y_yyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 381);

    auto ta1_y_yyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 382);

    auto ta1_y_yyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 383);

    auto ta1_y_yyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 384);

    auto ta1_y_yyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 385);

    auto ta1_y_yyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 386);

    auto ta1_y_yyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 387);

    auto ta1_y_yyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 388);

    auto ta1_y_yyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 389);

    auto ta1_y_yyyz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 390);

    auto ta1_y_yyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 391);

    auto ta1_y_yyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 392);

    auto ta1_y_yyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 393);

    auto ta1_y_yyyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 394);

    auto ta1_y_yyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 395);

    auto ta1_y_yyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 396);

    auto ta1_y_yyyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 397);

    auto ta1_y_yyyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 398);

    auto ta1_y_yyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 399);

    auto ta1_y_yyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 400);

    auto ta1_y_yyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 401);

    auto ta1_y_yyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 402);

    auto ta1_y_yyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 403);

    auto ta1_y_yyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 404);

    auto ta1_y_yyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 405);

    auto ta1_y_yyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 406);

    auto ta1_y_yyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 407);

    auto ta1_y_yyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 408);

    auto ta1_y_yyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 409);

    auto ta1_y_yyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 410);

    auto ta1_y_yyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 411);

    auto ta1_y_yyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 412);

    auto ta1_y_yyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 413);

    auto ta1_y_yyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 414);

    auto ta1_y_yyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 415);

    auto ta1_y_yyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 416);

    auto ta1_y_yyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 417);

    auto ta1_y_yyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 418);

    auto ta1_y_yyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 419);

    auto ta1_y_yzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 421);

    auto ta1_y_yzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 422);

    auto ta1_y_yzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 423);

    auto ta1_y_yzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 424);

    auto ta1_y_yzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 425);

    auto ta1_y_yzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 426);

    auto ta1_y_yzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 427);

    auto ta1_y_yzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 428);

    auto ta1_y_yzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 429);

    auto ta1_y_yzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 430);

    auto ta1_y_yzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 431);

    auto ta1_y_yzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 432);

    auto ta1_y_yzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 433);

    auto ta1_y_yzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 434);

    auto ta1_y_zzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 435);

    auto ta1_y_zzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 436);

    auto ta1_y_zzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 437);

    auto ta1_y_zzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 438);

    auto ta1_y_zzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 439);

    auto ta1_y_zzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 440);

    auto ta1_y_zzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 441);

    auto ta1_y_zzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 442);

    auto ta1_y_zzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 443);

    auto ta1_y_zzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 444);

    auto ta1_y_zzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 445);

    auto ta1_y_zzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 446);

    auto ta1_y_zzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 447);

    auto ta1_y_zzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 448);

    auto ta1_y_zzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 449);

    auto ta1_z_xxxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 450);

    auto ta1_z_xxxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 451);

    auto ta1_z_xxxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 452);

    auto ta1_z_xxxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 453);

    auto ta1_z_xxxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 454);

    auto ta1_z_xxxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 455);

    auto ta1_z_xxxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 456);

    auto ta1_z_xxxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 457);

    auto ta1_z_xxxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 458);

    auto ta1_z_xxxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 459);

    auto ta1_z_xxxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 460);

    auto ta1_z_xxxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 461);

    auto ta1_z_xxxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 462);

    auto ta1_z_xxxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 463);

    auto ta1_z_xxxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 464);

    auto ta1_z_xxxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 465);

    auto ta1_z_xxxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 466);

    auto ta1_z_xxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 467);

    auto ta1_z_xxxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 468);

    auto ta1_z_xxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 470);

    auto ta1_z_xxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 471);

    auto ta1_z_xxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 474);

    auto ta1_z_xxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 475);

    auto ta1_z_xxxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 476);

    auto ta1_z_xxxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 477);

    auto ta1_z_xxxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 478);

    auto ta1_z_xxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 480);

    auto ta1_z_xxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 481);

    auto ta1_z_xxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 482);

    auto ta1_z_xxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 483);

    auto ta1_z_xxxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 484);

    auto ta1_z_xxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 485);

    auto ta1_z_xxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 486);

    auto ta1_z_xxxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 487);

    auto ta1_z_xxxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 488);

    auto ta1_z_xxxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 489);

    auto ta1_z_xxxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 491);

    auto ta1_z_xxxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 492);

    auto ta1_z_xxxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 493);

    auto ta1_z_xxxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 494);

    auto ta1_z_xxyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 495);

    auto ta1_z_xxyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 496);

    auto ta1_z_xxyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 497);

    auto ta1_z_xxyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 498);

    auto ta1_z_xxyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 499);

    auto ta1_z_xxyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 500);

    auto ta1_z_xxyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 501);

    auto ta1_z_xxyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 502);

    auto ta1_z_xxyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 503);

    auto ta1_z_xxyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 504);

    auto ta1_z_xxyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 505);

    auto ta1_z_xxyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 506);

    auto ta1_z_xxyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 507);

    auto ta1_z_xxyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 508);

    auto ta1_z_xxyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 509);

    auto ta1_z_xxyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 512);

    auto ta1_z_xxyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 515);

    auto ta1_z_xxyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 519);

    auto ta1_z_xxyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 521);

    auto ta1_z_xxyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 522);

    auto ta1_z_xxyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 523);

    auto ta1_z_xxzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 525);

    auto ta1_z_xxzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 526);

    auto ta1_z_xxzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 527);

    auto ta1_z_xxzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 528);

    auto ta1_z_xxzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 529);

    auto ta1_z_xxzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 530);

    auto ta1_z_xxzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 531);

    auto ta1_z_xxzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 532);

    auto ta1_z_xxzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 533);

    auto ta1_z_xxzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 534);

    auto ta1_z_xxzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 535);

    auto ta1_z_xxzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 536);

    auto ta1_z_xxzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 537);

    auto ta1_z_xxzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 538);

    auto ta1_z_xxzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 539);

    auto ta1_z_xyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 540);

    auto ta1_z_xyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 541);

    auto ta1_z_xyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 543);

    auto ta1_z_xyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 544);

    auto ta1_z_xyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 546);

    auto ta1_z_xyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 547);

    auto ta1_z_xyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 548);

    auto ta1_z_xyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 550);

    auto ta1_z_xyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 551);

    auto ta1_z_xyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 552);

    auto ta1_z_xyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 553);

    auto ta1_z_xyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 554);

    auto ta1_z_xyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 566);

    auto ta1_z_xyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 567);

    auto ta1_z_xyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 568);

    auto ta1_z_xyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 569);

    auto ta1_z_xyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 580);

    auto ta1_z_xyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 581);

    auto ta1_z_xyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 582);

    auto ta1_z_xyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 583);

    auto ta1_z_xzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 585);

    auto ta1_z_xzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 587);

    auto ta1_z_xzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 589);

    auto ta1_z_xzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 590);

    auto ta1_z_xzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 592);

    auto ta1_z_xzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 593);

    auto ta1_z_xzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 594);

    auto ta1_z_xzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 595);

    auto ta1_z_xzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 596);

    auto ta1_z_xzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 597);

    auto ta1_z_xzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 598);

    auto ta1_z_xzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 599);

    auto ta1_z_yyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 600);

    auto ta1_z_yyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 601);

    auto ta1_z_yyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 602);

    auto ta1_z_yyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 603);

    auto ta1_z_yyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 604);

    auto ta1_z_yyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 605);

    auto ta1_z_yyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 606);

    auto ta1_z_yyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 607);

    auto ta1_z_yyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 608);

    auto ta1_z_yyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 609);

    auto ta1_z_yyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 610);

    auto ta1_z_yyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 611);

    auto ta1_z_yyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 612);

    auto ta1_z_yyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 613);

    auto ta1_z_yyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 614);

    auto ta1_z_yyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 616);

    auto ta1_z_yyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 617);

    auto ta1_z_yyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 618);

    auto ta1_z_yyyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 619);

    auto ta1_z_yyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 620);

    auto ta1_z_yyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 621);

    auto ta1_z_yyyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 622);

    auto ta1_z_yyyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 623);

    auto ta1_z_yyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 624);

    auto ta1_z_yyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 625);

    auto ta1_z_yyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 626);

    auto ta1_z_yyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 627);

    auto ta1_z_yyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 628);

    auto ta1_z_yyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 629);

    auto ta1_z_yyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 630);

    auto ta1_z_yyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 631);

    auto ta1_z_yyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 632);

    auto ta1_z_yyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 633);

    auto ta1_z_yyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 634);

    auto ta1_z_yyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 635);

    auto ta1_z_yyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 636);

    auto ta1_z_yyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 637);

    auto ta1_z_yyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 638);

    auto ta1_z_yyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 639);

    auto ta1_z_yyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 640);

    auto ta1_z_yyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 641);

    auto ta1_z_yyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 642);

    auto ta1_z_yyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 643);

    auto ta1_z_yyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 644);

    auto ta1_z_yzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 645);

    auto ta1_z_yzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 646);

    auto ta1_z_yzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 647);

    auto ta1_z_yzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 648);

    auto ta1_z_yzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 649);

    auto ta1_z_yzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 650);

    auto ta1_z_yzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 651);

    auto ta1_z_yzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 652);

    auto ta1_z_yzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 653);

    auto ta1_z_yzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 654);

    auto ta1_z_yzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 655);

    auto ta1_z_yzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 656);

    auto ta1_z_yzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 657);

    auto ta1_z_yzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 658);

    auto ta1_z_yzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 659);

    auto ta1_z_zzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 660);

    auto ta1_z_zzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 661);

    auto ta1_z_zzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 662);

    auto ta1_z_zzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 663);

    auto ta1_z_zzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 664);

    auto ta1_z_zzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 665);

    auto ta1_z_zzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 666);

    auto ta1_z_zzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 667);

    auto ta1_z_zzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 668);

    auto ta1_z_zzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 669);

    auto ta1_z_zzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 670);

    auto ta1_z_zzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 671);

    auto ta1_z_zzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 672);

    auto ta1_z_zzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 673);

    auto ta1_z_zzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 674);

    // Set up 0-15 components of targeted buffer : HG

    auto ta1_x_xxxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg);

    auto ta1_x_xxxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 1);

    auto ta1_x_xxxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 2);

    auto ta1_x_xxxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 3);

    auto ta1_x_xxxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 4);

    auto ta1_x_xxxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 5);

    auto ta1_x_xxxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 6);

    auto ta1_x_xxxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 7);

    auto ta1_x_xxxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 8);

    auto ta1_x_xxxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 9);

    auto ta1_x_xxxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 10);

    auto ta1_x_xxxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 11);

    auto ta1_x_xxxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 12);

    auto ta1_x_xxxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 13);

    auto ta1_x_xxxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 14);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_x_xxx_xxxx_0,   \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxy_0,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxz_0,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxyy_0,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyz_0,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxzz_0,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xyyy_0,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyz_0,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyzz_0,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xzzz_0,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_yyyy_0,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_x_xxx_yyyz_0,   \
                             ta1_x_xxx_yyyz_1,   \
                             ta1_x_xxx_yyzz_0,   \
                             ta1_x_xxx_yyzz_1,   \
                             ta1_x_xxx_yzzz_0,   \
                             ta1_x_xxx_yzzz_1,   \
                             ta1_x_xxx_zzzz_0,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta1_x_xxxx_xxx_0,   \
                             ta1_x_xxxx_xxx_1,   \
                             ta1_x_xxxx_xxxx_0,  \
                             ta1_x_xxxx_xxxx_1,  \
                             ta1_x_xxxx_xxxy_0,  \
                             ta1_x_xxxx_xxxy_1,  \
                             ta1_x_xxxx_xxxz_0,  \
                             ta1_x_xxxx_xxxz_1,  \
                             ta1_x_xxxx_xxy_0,   \
                             ta1_x_xxxx_xxy_1,   \
                             ta1_x_xxxx_xxyy_0,  \
                             ta1_x_xxxx_xxyy_1,  \
                             ta1_x_xxxx_xxyz_0,  \
                             ta1_x_xxxx_xxyz_1,  \
                             ta1_x_xxxx_xxz_0,   \
                             ta1_x_xxxx_xxz_1,   \
                             ta1_x_xxxx_xxzz_0,  \
                             ta1_x_xxxx_xxzz_1,  \
                             ta1_x_xxxx_xyy_0,   \
                             ta1_x_xxxx_xyy_1,   \
                             ta1_x_xxxx_xyyy_0,  \
                             ta1_x_xxxx_xyyy_1,  \
                             ta1_x_xxxx_xyyz_0,  \
                             ta1_x_xxxx_xyyz_1,  \
                             ta1_x_xxxx_xyz_0,   \
                             ta1_x_xxxx_xyz_1,   \
                             ta1_x_xxxx_xyzz_0,  \
                             ta1_x_xxxx_xyzz_1,  \
                             ta1_x_xxxx_xzz_0,   \
                             ta1_x_xxxx_xzz_1,   \
                             ta1_x_xxxx_xzzz_0,  \
                             ta1_x_xxxx_xzzz_1,  \
                             ta1_x_xxxx_yyy_0,   \
                             ta1_x_xxxx_yyy_1,   \
                             ta1_x_xxxx_yyyy_0,  \
                             ta1_x_xxxx_yyyy_1,  \
                             ta1_x_xxxx_yyyz_0,  \
                             ta1_x_xxxx_yyyz_1,  \
                             ta1_x_xxxx_yyz_0,   \
                             ta1_x_xxxx_yyz_1,   \
                             ta1_x_xxxx_yyzz_0,  \
                             ta1_x_xxxx_yyzz_1,  \
                             ta1_x_xxxx_yzz_0,   \
                             ta1_x_xxxx_yzz_1,   \
                             ta1_x_xxxx_yzzz_0,  \
                             ta1_x_xxxx_yzzz_1,  \
                             ta1_x_xxxx_zzz_0,   \
                             ta1_x_xxxx_zzz_1,   \
                             ta1_x_xxxx_zzzz_0,  \
                             ta1_x_xxxx_zzzz_1,  \
                             ta1_x_xxxxx_xxxx_0, \
                             ta1_x_xxxxx_xxxy_0, \
                             ta1_x_xxxxx_xxxz_0, \
                             ta1_x_xxxxx_xxyy_0, \
                             ta1_x_xxxxx_xxyz_0, \
                             ta1_x_xxxxx_xxzz_0, \
                             ta1_x_xxxxx_xyyy_0, \
                             ta1_x_xxxxx_xyyz_0, \
                             ta1_x_xxxxx_xyzz_0, \
                             ta1_x_xxxxx_xzzz_0, \
                             ta1_x_xxxxx_yyyy_0, \
                             ta1_x_xxxxx_yyyz_0, \
                             ta1_x_xxxxx_yyzz_0, \
                             ta1_x_xxxxx_yzzz_0, \
                             ta1_x_xxxxx_zzzz_0, \
                             ta_xxxx_xxxx_1,     \
                             ta_xxxx_xxxy_1,     \
                             ta_xxxx_xxxz_1,     \
                             ta_xxxx_xxyy_1,     \
                             ta_xxxx_xxyz_1,     \
                             ta_xxxx_xxzz_1,     \
                             ta_xxxx_xyyy_1,     \
                             ta_xxxx_xyyz_1,     \
                             ta_xxxx_xyzz_1,     \
                             ta_xxxx_xzzz_1,     \
                             ta_xxxx_yyyy_1,     \
                             ta_xxxx_yyyz_1,     \
                             ta_xxxx_yyzz_1,     \
                             ta_xxxx_yzzz_1,     \
                             ta_xxxx_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxx_xxxx_0[i] = 4.0 * ta1_x_xxx_xxxx_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxxx_1[i] * fe_0 + 4.0 * ta1_x_xxxx_xxx_0[i] * fe_0 -
                                4.0 * ta1_x_xxxx_xxx_1[i] * fe_0 + ta_xxxx_xxxx_1[i] + ta1_x_xxxx_xxxx_0[i] * pa_x[i] -
                                ta1_x_xxxx_xxxx_1[i] * pc_x[i];

        ta1_x_xxxxx_xxxy_0[i] = 4.0 * ta1_x_xxx_xxxy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxxy_1[i] * fe_0 + 3.0 * ta1_x_xxxx_xxy_0[i] * fe_0 -
                                3.0 * ta1_x_xxxx_xxy_1[i] * fe_0 + ta_xxxx_xxxy_1[i] + ta1_x_xxxx_xxxy_0[i] * pa_x[i] -
                                ta1_x_xxxx_xxxy_1[i] * pc_x[i];

        ta1_x_xxxxx_xxxz_0[i] = 4.0 * ta1_x_xxx_xxxz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxxz_1[i] * fe_0 + 3.0 * ta1_x_xxxx_xxz_0[i] * fe_0 -
                                3.0 * ta1_x_xxxx_xxz_1[i] * fe_0 + ta_xxxx_xxxz_1[i] + ta1_x_xxxx_xxxz_0[i] * pa_x[i] -
                                ta1_x_xxxx_xxxz_1[i] * pc_x[i];

        ta1_x_xxxxx_xxyy_0[i] = 4.0 * ta1_x_xxx_xxyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxxx_xyy_0[i] * fe_0 -
                                2.0 * ta1_x_xxxx_xyy_1[i] * fe_0 + ta_xxxx_xxyy_1[i] + ta1_x_xxxx_xxyy_0[i] * pa_x[i] -
                                ta1_x_xxxx_xxyy_1[i] * pc_x[i];

        ta1_x_xxxxx_xxyz_0[i] = 4.0 * ta1_x_xxx_xxyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxyz_1[i] * fe_0 + 2.0 * ta1_x_xxxx_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxx_xyz_1[i] * fe_0 + ta_xxxx_xxyz_1[i] + ta1_x_xxxx_xxyz_0[i] * pa_x[i] -
                                ta1_x_xxxx_xxyz_1[i] * pc_x[i];

        ta1_x_xxxxx_xxzz_0[i] = 4.0 * ta1_x_xxx_xxzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxxx_xzz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxx_xzz_1[i] * fe_0 + ta_xxxx_xxzz_1[i] + ta1_x_xxxx_xxzz_0[i] * pa_x[i] -
                                ta1_x_xxxx_xxzz_1[i] * pc_x[i];

        ta1_x_xxxxx_xyyy_0[i] = 4.0 * ta1_x_xxx_xyyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyyy_1[i] * fe_0 + ta1_x_xxxx_yyy_0[i] * fe_0 -
                                ta1_x_xxxx_yyy_1[i] * fe_0 + ta_xxxx_xyyy_1[i] + ta1_x_xxxx_xyyy_0[i] * pa_x[i] - ta1_x_xxxx_xyyy_1[i] * pc_x[i];

        ta1_x_xxxxx_xyyz_0[i] = 4.0 * ta1_x_xxx_xyyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyyz_1[i] * fe_0 + ta1_x_xxxx_yyz_0[i] * fe_0 -
                                ta1_x_xxxx_yyz_1[i] * fe_0 + ta_xxxx_xyyz_1[i] + ta1_x_xxxx_xyyz_0[i] * pa_x[i] - ta1_x_xxxx_xyyz_1[i] * pc_x[i];

        ta1_x_xxxxx_xyzz_0[i] = 4.0 * ta1_x_xxx_xyzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyzz_1[i] * fe_0 + ta1_x_xxxx_yzz_0[i] * fe_0 -
                                ta1_x_xxxx_yzz_1[i] * fe_0 + ta_xxxx_xyzz_1[i] + ta1_x_xxxx_xyzz_0[i] * pa_x[i] - ta1_x_xxxx_xyzz_1[i] * pc_x[i];

        ta1_x_xxxxx_xzzz_0[i] = 4.0 * ta1_x_xxx_xzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xzzz_1[i] * fe_0 + ta1_x_xxxx_zzz_0[i] * fe_0 -
                                ta1_x_xxxx_zzz_1[i] * fe_0 + ta_xxxx_xzzz_1[i] + ta1_x_xxxx_xzzz_0[i] * pa_x[i] - ta1_x_xxxx_xzzz_1[i] * pc_x[i];

        ta1_x_xxxxx_yyyy_0[i] = 4.0 * ta1_x_xxx_yyyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyyy_1[i] * fe_0 + ta_xxxx_yyyy_1[i] +
                                ta1_x_xxxx_yyyy_0[i] * pa_x[i] - ta1_x_xxxx_yyyy_1[i] * pc_x[i];

        ta1_x_xxxxx_yyyz_0[i] = 4.0 * ta1_x_xxx_yyyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyyz_1[i] * fe_0 + ta_xxxx_yyyz_1[i] +
                                ta1_x_xxxx_yyyz_0[i] * pa_x[i] - ta1_x_xxxx_yyyz_1[i] * pc_x[i];

        ta1_x_xxxxx_yyzz_0[i] = 4.0 * ta1_x_xxx_yyzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyzz_1[i] * fe_0 + ta_xxxx_yyzz_1[i] +
                                ta1_x_xxxx_yyzz_0[i] * pa_x[i] - ta1_x_xxxx_yyzz_1[i] * pc_x[i];

        ta1_x_xxxxx_yzzz_0[i] = 4.0 * ta1_x_xxx_yzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yzzz_1[i] * fe_0 + ta_xxxx_yzzz_1[i] +
                                ta1_x_xxxx_yzzz_0[i] * pa_x[i] - ta1_x_xxxx_yzzz_1[i] * pc_x[i];

        ta1_x_xxxxx_zzzz_0[i] = 4.0 * ta1_x_xxx_zzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_zzzz_1[i] * fe_0 + ta_xxxx_zzzz_1[i] +
                                ta1_x_xxxx_zzzz_0[i] * pa_x[i] - ta1_x_xxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : HG

    auto ta1_x_xxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 15);

    auto ta1_x_xxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 16);

    auto ta1_x_xxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 17);

    auto ta1_x_xxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 18);

    auto ta1_x_xxxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 19);

    auto ta1_x_xxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 20);

    auto ta1_x_xxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 21);

    auto ta1_x_xxxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 22);

    auto ta1_x_xxxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 23);

    auto ta1_x_xxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 24);

    auto ta1_x_xxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 25);

    auto ta1_x_xxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 26);

    auto ta1_x_xxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 27);

    auto ta1_x_xxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 28);

    auto ta1_x_xxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 29);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xxxx_xxx_0,   \
                             ta1_x_xxxx_xxx_1,   \
                             ta1_x_xxxx_xxxx_0,  \
                             ta1_x_xxxx_xxxx_1,  \
                             ta1_x_xxxx_xxxy_0,  \
                             ta1_x_xxxx_xxxy_1,  \
                             ta1_x_xxxx_xxxz_0,  \
                             ta1_x_xxxx_xxxz_1,  \
                             ta1_x_xxxx_xxy_0,   \
                             ta1_x_xxxx_xxy_1,   \
                             ta1_x_xxxx_xxyy_0,  \
                             ta1_x_xxxx_xxyy_1,  \
                             ta1_x_xxxx_xxyz_0,  \
                             ta1_x_xxxx_xxyz_1,  \
                             ta1_x_xxxx_xxz_0,   \
                             ta1_x_xxxx_xxz_1,   \
                             ta1_x_xxxx_xxzz_0,  \
                             ta1_x_xxxx_xxzz_1,  \
                             ta1_x_xxxx_xyy_0,   \
                             ta1_x_xxxx_xyy_1,   \
                             ta1_x_xxxx_xyyy_0,  \
                             ta1_x_xxxx_xyyy_1,  \
                             ta1_x_xxxx_xyyz_0,  \
                             ta1_x_xxxx_xyyz_1,  \
                             ta1_x_xxxx_xyz_0,   \
                             ta1_x_xxxx_xyz_1,   \
                             ta1_x_xxxx_xyzz_0,  \
                             ta1_x_xxxx_xyzz_1,  \
                             ta1_x_xxxx_xzz_0,   \
                             ta1_x_xxxx_xzz_1,   \
                             ta1_x_xxxx_xzzz_0,  \
                             ta1_x_xxxx_xzzz_1,  \
                             ta1_x_xxxx_yyy_0,   \
                             ta1_x_xxxx_yyy_1,   \
                             ta1_x_xxxx_yyyy_0,  \
                             ta1_x_xxxx_yyyy_1,  \
                             ta1_x_xxxx_yyyz_0,  \
                             ta1_x_xxxx_yyyz_1,  \
                             ta1_x_xxxx_yyz_0,   \
                             ta1_x_xxxx_yyz_1,   \
                             ta1_x_xxxx_yyzz_0,  \
                             ta1_x_xxxx_yyzz_1,  \
                             ta1_x_xxxx_yzz_0,   \
                             ta1_x_xxxx_yzz_1,   \
                             ta1_x_xxxx_yzzz_0,  \
                             ta1_x_xxxx_yzzz_1,  \
                             ta1_x_xxxx_zzz_0,   \
                             ta1_x_xxxx_zzz_1,   \
                             ta1_x_xxxx_zzzz_0,  \
                             ta1_x_xxxx_zzzz_1,  \
                             ta1_x_xxxxy_xxxx_0, \
                             ta1_x_xxxxy_xxxy_0, \
                             ta1_x_xxxxy_xxxz_0, \
                             ta1_x_xxxxy_xxyy_0, \
                             ta1_x_xxxxy_xxyz_0, \
                             ta1_x_xxxxy_xxzz_0, \
                             ta1_x_xxxxy_xyyy_0, \
                             ta1_x_xxxxy_xyyz_0, \
                             ta1_x_xxxxy_xyzz_0, \
                             ta1_x_xxxxy_xzzz_0, \
                             ta1_x_xxxxy_yyyy_0, \
                             ta1_x_xxxxy_yyyz_0, \
                             ta1_x_xxxxy_yyzz_0, \
                             ta1_x_xxxxy_yzzz_0, \
                             ta1_x_xxxxy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxy_xxxx_0[i] = ta1_x_xxxx_xxxx_0[i] * pa_y[i] - ta1_x_xxxx_xxxx_1[i] * pc_y[i];

        ta1_x_xxxxy_xxxy_0[i] =
            ta1_x_xxxx_xxx_0[i] * fe_0 - ta1_x_xxxx_xxx_1[i] * fe_0 + ta1_x_xxxx_xxxy_0[i] * pa_y[i] - ta1_x_xxxx_xxxy_1[i] * pc_y[i];

        ta1_x_xxxxy_xxxz_0[i] = ta1_x_xxxx_xxxz_0[i] * pa_y[i] - ta1_x_xxxx_xxxz_1[i] * pc_y[i];

        ta1_x_xxxxy_xxyy_0[i] =
            2.0 * ta1_x_xxxx_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xxy_1[i] * fe_0 + ta1_x_xxxx_xxyy_0[i] * pa_y[i] - ta1_x_xxxx_xxyy_1[i] * pc_y[i];

        ta1_x_xxxxy_xxyz_0[i] =
            ta1_x_xxxx_xxz_0[i] * fe_0 - ta1_x_xxxx_xxz_1[i] * fe_0 + ta1_x_xxxx_xxyz_0[i] * pa_y[i] - ta1_x_xxxx_xxyz_1[i] * pc_y[i];

        ta1_x_xxxxy_xxzz_0[i] = ta1_x_xxxx_xxzz_0[i] * pa_y[i] - ta1_x_xxxx_xxzz_1[i] * pc_y[i];

        ta1_x_xxxxy_xyyy_0[i] =
            3.0 * ta1_x_xxxx_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxxx_xyy_1[i] * fe_0 + ta1_x_xxxx_xyyy_0[i] * pa_y[i] - ta1_x_xxxx_xyyy_1[i] * pc_y[i];

        ta1_x_xxxxy_xyyz_0[i] =
            2.0 * ta1_x_xxxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xyz_1[i] * fe_0 + ta1_x_xxxx_xyyz_0[i] * pa_y[i] - ta1_x_xxxx_xyyz_1[i] * pc_y[i];

        ta1_x_xxxxy_xyzz_0[i] =
            ta1_x_xxxx_xzz_0[i] * fe_0 - ta1_x_xxxx_xzz_1[i] * fe_0 + ta1_x_xxxx_xyzz_0[i] * pa_y[i] - ta1_x_xxxx_xyzz_1[i] * pc_y[i];

        ta1_x_xxxxy_xzzz_0[i] = ta1_x_xxxx_xzzz_0[i] * pa_y[i] - ta1_x_xxxx_xzzz_1[i] * pc_y[i];

        ta1_x_xxxxy_yyyy_0[i] =
            4.0 * ta1_x_xxxx_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxxx_yyy_1[i] * fe_0 + ta1_x_xxxx_yyyy_0[i] * pa_y[i] - ta1_x_xxxx_yyyy_1[i] * pc_y[i];

        ta1_x_xxxxy_yyyz_0[i] =
            3.0 * ta1_x_xxxx_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxxx_yyz_1[i] * fe_0 + ta1_x_xxxx_yyyz_0[i] * pa_y[i] - ta1_x_xxxx_yyyz_1[i] * pc_y[i];

        ta1_x_xxxxy_yyzz_0[i] =
            2.0 * ta1_x_xxxx_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_yzz_1[i] * fe_0 + ta1_x_xxxx_yyzz_0[i] * pa_y[i] - ta1_x_xxxx_yyzz_1[i] * pc_y[i];

        ta1_x_xxxxy_yzzz_0[i] =
            ta1_x_xxxx_zzz_0[i] * fe_0 - ta1_x_xxxx_zzz_1[i] * fe_0 + ta1_x_xxxx_yzzz_0[i] * pa_y[i] - ta1_x_xxxx_yzzz_1[i] * pc_y[i];

        ta1_x_xxxxy_zzzz_0[i] = ta1_x_xxxx_zzzz_0[i] * pa_y[i] - ta1_x_xxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : HG

    auto ta1_x_xxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 30);

    auto ta1_x_xxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 31);

    auto ta1_x_xxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 32);

    auto ta1_x_xxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 33);

    auto ta1_x_xxxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 34);

    auto ta1_x_xxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 35);

    auto ta1_x_xxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 36);

    auto ta1_x_xxxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 37);

    auto ta1_x_xxxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 38);

    auto ta1_x_xxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 39);

    auto ta1_x_xxxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 40);

    auto ta1_x_xxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 41);

    auto ta1_x_xxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 42);

    auto ta1_x_xxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 43);

    auto ta1_x_xxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 44);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_xxxx_xxx_0,   \
                             ta1_x_xxxx_xxx_1,   \
                             ta1_x_xxxx_xxxx_0,  \
                             ta1_x_xxxx_xxxx_1,  \
                             ta1_x_xxxx_xxxy_0,  \
                             ta1_x_xxxx_xxxy_1,  \
                             ta1_x_xxxx_xxxz_0,  \
                             ta1_x_xxxx_xxxz_1,  \
                             ta1_x_xxxx_xxy_0,   \
                             ta1_x_xxxx_xxy_1,   \
                             ta1_x_xxxx_xxyy_0,  \
                             ta1_x_xxxx_xxyy_1,  \
                             ta1_x_xxxx_xxyz_0,  \
                             ta1_x_xxxx_xxyz_1,  \
                             ta1_x_xxxx_xxz_0,   \
                             ta1_x_xxxx_xxz_1,   \
                             ta1_x_xxxx_xxzz_0,  \
                             ta1_x_xxxx_xxzz_1,  \
                             ta1_x_xxxx_xyy_0,   \
                             ta1_x_xxxx_xyy_1,   \
                             ta1_x_xxxx_xyyy_0,  \
                             ta1_x_xxxx_xyyy_1,  \
                             ta1_x_xxxx_xyyz_0,  \
                             ta1_x_xxxx_xyyz_1,  \
                             ta1_x_xxxx_xyz_0,   \
                             ta1_x_xxxx_xyz_1,   \
                             ta1_x_xxxx_xyzz_0,  \
                             ta1_x_xxxx_xyzz_1,  \
                             ta1_x_xxxx_xzz_0,   \
                             ta1_x_xxxx_xzz_1,   \
                             ta1_x_xxxx_xzzz_0,  \
                             ta1_x_xxxx_xzzz_1,  \
                             ta1_x_xxxx_yyy_0,   \
                             ta1_x_xxxx_yyy_1,   \
                             ta1_x_xxxx_yyyy_0,  \
                             ta1_x_xxxx_yyyy_1,  \
                             ta1_x_xxxx_yyyz_0,  \
                             ta1_x_xxxx_yyyz_1,  \
                             ta1_x_xxxx_yyz_0,   \
                             ta1_x_xxxx_yyz_1,   \
                             ta1_x_xxxx_yyzz_0,  \
                             ta1_x_xxxx_yyzz_1,  \
                             ta1_x_xxxx_yzz_0,   \
                             ta1_x_xxxx_yzz_1,   \
                             ta1_x_xxxx_yzzz_0,  \
                             ta1_x_xxxx_yzzz_1,  \
                             ta1_x_xxxx_zzz_0,   \
                             ta1_x_xxxx_zzz_1,   \
                             ta1_x_xxxx_zzzz_0,  \
                             ta1_x_xxxx_zzzz_1,  \
                             ta1_x_xxxxz_xxxx_0, \
                             ta1_x_xxxxz_xxxy_0, \
                             ta1_x_xxxxz_xxxz_0, \
                             ta1_x_xxxxz_xxyy_0, \
                             ta1_x_xxxxz_xxyz_0, \
                             ta1_x_xxxxz_xxzz_0, \
                             ta1_x_xxxxz_xyyy_0, \
                             ta1_x_xxxxz_xyyz_0, \
                             ta1_x_xxxxz_xyzz_0, \
                             ta1_x_xxxxz_xzzz_0, \
                             ta1_x_xxxxz_yyyy_0, \
                             ta1_x_xxxxz_yyyz_0, \
                             ta1_x_xxxxz_yyzz_0, \
                             ta1_x_xxxxz_yzzz_0, \
                             ta1_x_xxxxz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxz_xxxx_0[i] = ta1_x_xxxx_xxxx_0[i] * pa_z[i] - ta1_x_xxxx_xxxx_1[i] * pc_z[i];

        ta1_x_xxxxz_xxxy_0[i] = ta1_x_xxxx_xxxy_0[i] * pa_z[i] - ta1_x_xxxx_xxxy_1[i] * pc_z[i];

        ta1_x_xxxxz_xxxz_0[i] =
            ta1_x_xxxx_xxx_0[i] * fe_0 - ta1_x_xxxx_xxx_1[i] * fe_0 + ta1_x_xxxx_xxxz_0[i] * pa_z[i] - ta1_x_xxxx_xxxz_1[i] * pc_z[i];

        ta1_x_xxxxz_xxyy_0[i] = ta1_x_xxxx_xxyy_0[i] * pa_z[i] - ta1_x_xxxx_xxyy_1[i] * pc_z[i];

        ta1_x_xxxxz_xxyz_0[i] =
            ta1_x_xxxx_xxy_0[i] * fe_0 - ta1_x_xxxx_xxy_1[i] * fe_0 + ta1_x_xxxx_xxyz_0[i] * pa_z[i] - ta1_x_xxxx_xxyz_1[i] * pc_z[i];

        ta1_x_xxxxz_xxzz_0[i] =
            2.0 * ta1_x_xxxx_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xxz_1[i] * fe_0 + ta1_x_xxxx_xxzz_0[i] * pa_z[i] - ta1_x_xxxx_xxzz_1[i] * pc_z[i];

        ta1_x_xxxxz_xyyy_0[i] = ta1_x_xxxx_xyyy_0[i] * pa_z[i] - ta1_x_xxxx_xyyy_1[i] * pc_z[i];

        ta1_x_xxxxz_xyyz_0[i] =
            ta1_x_xxxx_xyy_0[i] * fe_0 - ta1_x_xxxx_xyy_1[i] * fe_0 + ta1_x_xxxx_xyyz_0[i] * pa_z[i] - ta1_x_xxxx_xyyz_1[i] * pc_z[i];

        ta1_x_xxxxz_xyzz_0[i] =
            2.0 * ta1_x_xxxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xyz_1[i] * fe_0 + ta1_x_xxxx_xyzz_0[i] * pa_z[i] - ta1_x_xxxx_xyzz_1[i] * pc_z[i];

        ta1_x_xxxxz_xzzz_0[i] =
            3.0 * ta1_x_xxxx_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxxx_xzz_1[i] * fe_0 + ta1_x_xxxx_xzzz_0[i] * pa_z[i] - ta1_x_xxxx_xzzz_1[i] * pc_z[i];

        ta1_x_xxxxz_yyyy_0[i] = ta1_x_xxxx_yyyy_0[i] * pa_z[i] - ta1_x_xxxx_yyyy_1[i] * pc_z[i];

        ta1_x_xxxxz_yyyz_0[i] =
            ta1_x_xxxx_yyy_0[i] * fe_0 - ta1_x_xxxx_yyy_1[i] * fe_0 + ta1_x_xxxx_yyyz_0[i] * pa_z[i] - ta1_x_xxxx_yyyz_1[i] * pc_z[i];

        ta1_x_xxxxz_yyzz_0[i] =
            2.0 * ta1_x_xxxx_yyz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_yyz_1[i] * fe_0 + ta1_x_xxxx_yyzz_0[i] * pa_z[i] - ta1_x_xxxx_yyzz_1[i] * pc_z[i];

        ta1_x_xxxxz_yzzz_0[i] =
            3.0 * ta1_x_xxxx_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxxx_yzz_1[i] * fe_0 + ta1_x_xxxx_yzzz_0[i] * pa_z[i] - ta1_x_xxxx_yzzz_1[i] * pc_z[i];

        ta1_x_xxxxz_zzzz_0[i] =
            4.0 * ta1_x_xxxx_zzz_0[i] * fe_0 - 4.0 * ta1_x_xxxx_zzz_1[i] * fe_0 + ta1_x_xxxx_zzzz_0[i] * pa_z[i] - ta1_x_xxxx_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : HG

    auto ta1_x_xxxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 45);

    auto ta1_x_xxxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 46);

    auto ta1_x_xxxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 47);

    auto ta1_x_xxxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 48);

    auto ta1_x_xxxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 49);

    auto ta1_x_xxxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 50);

    auto ta1_x_xxxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 51);

    auto ta1_x_xxxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 52);

    auto ta1_x_xxxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 53);

    auto ta1_x_xxxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 54);

    auto ta1_x_xxxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 55);

    auto ta1_x_xxxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 56);

    auto ta1_x_xxxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 57);

    auto ta1_x_xxxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 58);

    auto ta1_x_xxxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 59);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxx_xxxx_0,   \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxy_0,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxz_0,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxyy_0,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyz_0,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxzz_0,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xyyy_0,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyz_0,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyzz_0,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xzzz_0,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_zzzz_0,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta1_x_xxxy_xxx_0,   \
                             ta1_x_xxxy_xxx_1,   \
                             ta1_x_xxxy_xxxx_0,  \
                             ta1_x_xxxy_xxxx_1,  \
                             ta1_x_xxxy_xxxy_0,  \
                             ta1_x_xxxy_xxxy_1,  \
                             ta1_x_xxxy_xxxz_0,  \
                             ta1_x_xxxy_xxxz_1,  \
                             ta1_x_xxxy_xxy_0,   \
                             ta1_x_xxxy_xxy_1,   \
                             ta1_x_xxxy_xxyy_0,  \
                             ta1_x_xxxy_xxyy_1,  \
                             ta1_x_xxxy_xxyz_0,  \
                             ta1_x_xxxy_xxyz_1,  \
                             ta1_x_xxxy_xxz_0,   \
                             ta1_x_xxxy_xxz_1,   \
                             ta1_x_xxxy_xxzz_0,  \
                             ta1_x_xxxy_xxzz_1,  \
                             ta1_x_xxxy_xyy_0,   \
                             ta1_x_xxxy_xyy_1,   \
                             ta1_x_xxxy_xyyy_0,  \
                             ta1_x_xxxy_xyyy_1,  \
                             ta1_x_xxxy_xyyz_0,  \
                             ta1_x_xxxy_xyyz_1,  \
                             ta1_x_xxxy_xyz_0,   \
                             ta1_x_xxxy_xyz_1,   \
                             ta1_x_xxxy_xyzz_0,  \
                             ta1_x_xxxy_xyzz_1,  \
                             ta1_x_xxxy_xzz_0,   \
                             ta1_x_xxxy_xzz_1,   \
                             ta1_x_xxxy_xzzz_0,  \
                             ta1_x_xxxy_xzzz_1,  \
                             ta1_x_xxxy_zzzz_0,  \
                             ta1_x_xxxy_zzzz_1,  \
                             ta1_x_xxxyy_xxxx_0, \
                             ta1_x_xxxyy_xxxy_0, \
                             ta1_x_xxxyy_xxxz_0, \
                             ta1_x_xxxyy_xxyy_0, \
                             ta1_x_xxxyy_xxyz_0, \
                             ta1_x_xxxyy_xxzz_0, \
                             ta1_x_xxxyy_xyyy_0, \
                             ta1_x_xxxyy_xyyz_0, \
                             ta1_x_xxxyy_xyzz_0, \
                             ta1_x_xxxyy_xzzz_0, \
                             ta1_x_xxxyy_yyyy_0, \
                             ta1_x_xxxyy_yyyz_0, \
                             ta1_x_xxxyy_yyzz_0, \
                             ta1_x_xxxyy_yzzz_0, \
                             ta1_x_xxxyy_zzzz_0, \
                             ta1_x_xxyy_yyyy_0,  \
                             ta1_x_xxyy_yyyy_1,  \
                             ta1_x_xxyy_yyyz_0,  \
                             ta1_x_xxyy_yyyz_1,  \
                             ta1_x_xxyy_yyzz_0,  \
                             ta1_x_xxyy_yyzz_1,  \
                             ta1_x_xxyy_yzzz_0,  \
                             ta1_x_xxyy_yzzz_1,  \
                             ta1_x_xyy_yyyy_0,   \
                             ta1_x_xyy_yyyy_1,   \
                             ta1_x_xyy_yyyz_0,   \
                             ta1_x_xyy_yyyz_1,   \
                             ta1_x_xyy_yyzz_0,   \
                             ta1_x_xyy_yyzz_1,   \
                             ta1_x_xyy_yzzz_0,   \
                             ta1_x_xyy_yzzz_1,   \
                             ta_xxyy_yyyy_1,     \
                             ta_xxyy_yyyz_1,     \
                             ta_xxyy_yyzz_1,     \
                             ta_xxyy_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyy_xxxx_0[i] =
            ta1_x_xxx_xxxx_0[i] * fe_0 - ta1_x_xxx_xxxx_1[i] * fe_0 + ta1_x_xxxy_xxxx_0[i] * pa_y[i] - ta1_x_xxxy_xxxx_1[i] * pc_y[i];

        ta1_x_xxxyy_xxxy_0[i] = ta1_x_xxx_xxxy_0[i] * fe_0 - ta1_x_xxx_xxxy_1[i] * fe_0 + ta1_x_xxxy_xxx_0[i] * fe_0 - ta1_x_xxxy_xxx_1[i] * fe_0 +
                                ta1_x_xxxy_xxxy_0[i] * pa_y[i] - ta1_x_xxxy_xxxy_1[i] * pc_y[i];

        ta1_x_xxxyy_xxxz_0[i] =
            ta1_x_xxx_xxxz_0[i] * fe_0 - ta1_x_xxx_xxxz_1[i] * fe_0 + ta1_x_xxxy_xxxz_0[i] * pa_y[i] - ta1_x_xxxy_xxxz_1[i] * pc_y[i];

        ta1_x_xxxyy_xxyy_0[i] = ta1_x_xxx_xxyy_0[i] * fe_0 - ta1_x_xxx_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxxy_xxy_0[i] * fe_0 -
                                2.0 * ta1_x_xxxy_xxy_1[i] * fe_0 + ta1_x_xxxy_xxyy_0[i] * pa_y[i] - ta1_x_xxxy_xxyy_1[i] * pc_y[i];

        ta1_x_xxxyy_xxyz_0[i] = ta1_x_xxx_xxyz_0[i] * fe_0 - ta1_x_xxx_xxyz_1[i] * fe_0 + ta1_x_xxxy_xxz_0[i] * fe_0 - ta1_x_xxxy_xxz_1[i] * fe_0 +
                                ta1_x_xxxy_xxyz_0[i] * pa_y[i] - ta1_x_xxxy_xxyz_1[i] * pc_y[i];

        ta1_x_xxxyy_xxzz_0[i] =
            ta1_x_xxx_xxzz_0[i] * fe_0 - ta1_x_xxx_xxzz_1[i] * fe_0 + ta1_x_xxxy_xxzz_0[i] * pa_y[i] - ta1_x_xxxy_xxzz_1[i] * pc_y[i];

        ta1_x_xxxyy_xyyy_0[i] = ta1_x_xxx_xyyy_0[i] * fe_0 - ta1_x_xxx_xyyy_1[i] * fe_0 + 3.0 * ta1_x_xxxy_xyy_0[i] * fe_0 -
                                3.0 * ta1_x_xxxy_xyy_1[i] * fe_0 + ta1_x_xxxy_xyyy_0[i] * pa_y[i] - ta1_x_xxxy_xyyy_1[i] * pc_y[i];

        ta1_x_xxxyy_xyyz_0[i] = ta1_x_xxx_xyyz_0[i] * fe_0 - ta1_x_xxx_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxxy_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxy_xyz_1[i] * fe_0 + ta1_x_xxxy_xyyz_0[i] * pa_y[i] - ta1_x_xxxy_xyyz_1[i] * pc_y[i];

        ta1_x_xxxyy_xyzz_0[i] = ta1_x_xxx_xyzz_0[i] * fe_0 - ta1_x_xxx_xyzz_1[i] * fe_0 + ta1_x_xxxy_xzz_0[i] * fe_0 - ta1_x_xxxy_xzz_1[i] * fe_0 +
                                ta1_x_xxxy_xyzz_0[i] * pa_y[i] - ta1_x_xxxy_xyzz_1[i] * pc_y[i];

        ta1_x_xxxyy_xzzz_0[i] =
            ta1_x_xxx_xzzz_0[i] * fe_0 - ta1_x_xxx_xzzz_1[i] * fe_0 + ta1_x_xxxy_xzzz_0[i] * pa_y[i] - ta1_x_xxxy_xzzz_1[i] * pc_y[i];

        ta1_x_xxxyy_yyyy_0[i] = 2.0 * ta1_x_xyy_yyyy_0[i] * fe_0 - 2.0 * ta1_x_xyy_yyyy_1[i] * fe_0 + ta_xxyy_yyyy_1[i] +
                                ta1_x_xxyy_yyyy_0[i] * pa_x[i] - ta1_x_xxyy_yyyy_1[i] * pc_x[i];

        ta1_x_xxxyy_yyyz_0[i] = 2.0 * ta1_x_xyy_yyyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_yyyz_1[i] * fe_0 + ta_xxyy_yyyz_1[i] +
                                ta1_x_xxyy_yyyz_0[i] * pa_x[i] - ta1_x_xxyy_yyyz_1[i] * pc_x[i];

        ta1_x_xxxyy_yyzz_0[i] = 2.0 * ta1_x_xyy_yyzz_0[i] * fe_0 - 2.0 * ta1_x_xyy_yyzz_1[i] * fe_0 + ta_xxyy_yyzz_1[i] +
                                ta1_x_xxyy_yyzz_0[i] * pa_x[i] - ta1_x_xxyy_yyzz_1[i] * pc_x[i];

        ta1_x_xxxyy_yzzz_0[i] = 2.0 * ta1_x_xyy_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xyy_yzzz_1[i] * fe_0 + ta_xxyy_yzzz_1[i] +
                                ta1_x_xxyy_yzzz_0[i] * pa_x[i] - ta1_x_xxyy_yzzz_1[i] * pc_x[i];

        ta1_x_xxxyy_zzzz_0[i] =
            ta1_x_xxx_zzzz_0[i] * fe_0 - ta1_x_xxx_zzzz_1[i] * fe_0 + ta1_x_xxxy_zzzz_0[i] * pa_y[i] - ta1_x_xxxy_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : HG

    auto ta1_x_xxxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 60);

    auto ta1_x_xxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 61);

    auto ta1_x_xxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 62);

    auto ta1_x_xxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 63);

    auto ta1_x_xxxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 64);

    auto ta1_x_xxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 65);

    auto ta1_x_xxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 66);

    auto ta1_x_xxxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 67);

    auto ta1_x_xxxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 68);

    auto ta1_x_xxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 69);

    auto ta1_x_xxxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 70);

    auto ta1_x_xxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 71);

    auto ta1_x_xxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 72);

    auto ta1_x_xxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 73);

    auto ta1_x_xxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 74);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxxy_xxxy_0,  \
                             ta1_x_xxxy_xxxy_1,  \
                             ta1_x_xxxy_xxyy_0,  \
                             ta1_x_xxxy_xxyy_1,  \
                             ta1_x_xxxy_xyyy_0,  \
                             ta1_x_xxxy_xyyy_1,  \
                             ta1_x_xxxy_yyyy_0,  \
                             ta1_x_xxxy_yyyy_1,  \
                             ta1_x_xxxyz_xxxx_0, \
                             ta1_x_xxxyz_xxxy_0, \
                             ta1_x_xxxyz_xxxz_0, \
                             ta1_x_xxxyz_xxyy_0, \
                             ta1_x_xxxyz_xxyz_0, \
                             ta1_x_xxxyz_xxzz_0, \
                             ta1_x_xxxyz_xyyy_0, \
                             ta1_x_xxxyz_xyyz_0, \
                             ta1_x_xxxyz_xyzz_0, \
                             ta1_x_xxxyz_xzzz_0, \
                             ta1_x_xxxyz_yyyy_0, \
                             ta1_x_xxxyz_yyyz_0, \
                             ta1_x_xxxyz_yyzz_0, \
                             ta1_x_xxxyz_yzzz_0, \
                             ta1_x_xxxyz_zzzz_0, \
                             ta1_x_xxxz_xxxx_0,  \
                             ta1_x_xxxz_xxxx_1,  \
                             ta1_x_xxxz_xxxz_0,  \
                             ta1_x_xxxz_xxxz_1,  \
                             ta1_x_xxxz_xxyz_0,  \
                             ta1_x_xxxz_xxyz_1,  \
                             ta1_x_xxxz_xxz_0,   \
                             ta1_x_xxxz_xxz_1,   \
                             ta1_x_xxxz_xxzz_0,  \
                             ta1_x_xxxz_xxzz_1,  \
                             ta1_x_xxxz_xyyz_0,  \
                             ta1_x_xxxz_xyyz_1,  \
                             ta1_x_xxxz_xyz_0,   \
                             ta1_x_xxxz_xyz_1,   \
                             ta1_x_xxxz_xyzz_0,  \
                             ta1_x_xxxz_xyzz_1,  \
                             ta1_x_xxxz_xzz_0,   \
                             ta1_x_xxxz_xzz_1,   \
                             ta1_x_xxxz_xzzz_0,  \
                             ta1_x_xxxz_xzzz_1,  \
                             ta1_x_xxxz_yyyz_0,  \
                             ta1_x_xxxz_yyyz_1,  \
                             ta1_x_xxxz_yyz_0,   \
                             ta1_x_xxxz_yyz_1,   \
                             ta1_x_xxxz_yyzz_0,  \
                             ta1_x_xxxz_yyzz_1,  \
                             ta1_x_xxxz_yzz_0,   \
                             ta1_x_xxxz_yzz_1,   \
                             ta1_x_xxxz_yzzz_0,  \
                             ta1_x_xxxz_yzzz_1,  \
                             ta1_x_xxxz_zzz_0,   \
                             ta1_x_xxxz_zzz_1,   \
                             ta1_x_xxxz_zzzz_0,  \
                             ta1_x_xxxz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyz_xxxx_0[i] = ta1_x_xxxz_xxxx_0[i] * pa_y[i] - ta1_x_xxxz_xxxx_1[i] * pc_y[i];

        ta1_x_xxxyz_xxxy_0[i] = ta1_x_xxxy_xxxy_0[i] * pa_z[i] - ta1_x_xxxy_xxxy_1[i] * pc_z[i];

        ta1_x_xxxyz_xxxz_0[i] = ta1_x_xxxz_xxxz_0[i] * pa_y[i] - ta1_x_xxxz_xxxz_1[i] * pc_y[i];

        ta1_x_xxxyz_xxyy_0[i] = ta1_x_xxxy_xxyy_0[i] * pa_z[i] - ta1_x_xxxy_xxyy_1[i] * pc_z[i];

        ta1_x_xxxyz_xxyz_0[i] =
            ta1_x_xxxz_xxz_0[i] * fe_0 - ta1_x_xxxz_xxz_1[i] * fe_0 + ta1_x_xxxz_xxyz_0[i] * pa_y[i] - ta1_x_xxxz_xxyz_1[i] * pc_y[i];

        ta1_x_xxxyz_xxzz_0[i] = ta1_x_xxxz_xxzz_0[i] * pa_y[i] - ta1_x_xxxz_xxzz_1[i] * pc_y[i];

        ta1_x_xxxyz_xyyy_0[i] = ta1_x_xxxy_xyyy_0[i] * pa_z[i] - ta1_x_xxxy_xyyy_1[i] * pc_z[i];

        ta1_x_xxxyz_xyyz_0[i] =
            2.0 * ta1_x_xxxz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xyz_1[i] * fe_0 + ta1_x_xxxz_xyyz_0[i] * pa_y[i] - ta1_x_xxxz_xyyz_1[i] * pc_y[i];

        ta1_x_xxxyz_xyzz_0[i] =
            ta1_x_xxxz_xzz_0[i] * fe_0 - ta1_x_xxxz_xzz_1[i] * fe_0 + ta1_x_xxxz_xyzz_0[i] * pa_y[i] - ta1_x_xxxz_xyzz_1[i] * pc_y[i];

        ta1_x_xxxyz_xzzz_0[i] = ta1_x_xxxz_xzzz_0[i] * pa_y[i] - ta1_x_xxxz_xzzz_1[i] * pc_y[i];

        ta1_x_xxxyz_yyyy_0[i] = ta1_x_xxxy_yyyy_0[i] * pa_z[i] - ta1_x_xxxy_yyyy_1[i] * pc_z[i];

        ta1_x_xxxyz_yyyz_0[i] =
            3.0 * ta1_x_xxxz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxxz_yyz_1[i] * fe_0 + ta1_x_xxxz_yyyz_0[i] * pa_y[i] - ta1_x_xxxz_yyyz_1[i] * pc_y[i];

        ta1_x_xxxyz_yyzz_0[i] =
            2.0 * ta1_x_xxxz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_yzz_1[i] * fe_0 + ta1_x_xxxz_yyzz_0[i] * pa_y[i] - ta1_x_xxxz_yyzz_1[i] * pc_y[i];

        ta1_x_xxxyz_yzzz_0[i] =
            ta1_x_xxxz_zzz_0[i] * fe_0 - ta1_x_xxxz_zzz_1[i] * fe_0 + ta1_x_xxxz_yzzz_0[i] * pa_y[i] - ta1_x_xxxz_yzzz_1[i] * pc_y[i];

        ta1_x_xxxyz_zzzz_0[i] = ta1_x_xxxz_zzzz_0[i] * pa_y[i] - ta1_x_xxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : HG

    auto ta1_x_xxxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 75);

    auto ta1_x_xxxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 76);

    auto ta1_x_xxxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 77);

    auto ta1_x_xxxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 78);

    auto ta1_x_xxxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 79);

    auto ta1_x_xxxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 80);

    auto ta1_x_xxxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 81);

    auto ta1_x_xxxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 82);

    auto ta1_x_xxxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 83);

    auto ta1_x_xxxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 84);

    auto ta1_x_xxxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 85);

    auto ta1_x_xxxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 86);

    auto ta1_x_xxxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 87);

    auto ta1_x_xxxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 88);

    auto ta1_x_xxxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 89);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxx_xxxx_0,   \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxy_0,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxz_0,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxyy_0,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyz_0,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxzz_0,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xyyy_0,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyz_0,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyzz_0,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xzzz_0,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_yyyy_0,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_x_xxxz_xxx_0,   \
                             ta1_x_xxxz_xxx_1,   \
                             ta1_x_xxxz_xxxx_0,  \
                             ta1_x_xxxz_xxxx_1,  \
                             ta1_x_xxxz_xxxy_0,  \
                             ta1_x_xxxz_xxxy_1,  \
                             ta1_x_xxxz_xxxz_0,  \
                             ta1_x_xxxz_xxxz_1,  \
                             ta1_x_xxxz_xxy_0,   \
                             ta1_x_xxxz_xxy_1,   \
                             ta1_x_xxxz_xxyy_0,  \
                             ta1_x_xxxz_xxyy_1,  \
                             ta1_x_xxxz_xxyz_0,  \
                             ta1_x_xxxz_xxyz_1,  \
                             ta1_x_xxxz_xxz_0,   \
                             ta1_x_xxxz_xxz_1,   \
                             ta1_x_xxxz_xxzz_0,  \
                             ta1_x_xxxz_xxzz_1,  \
                             ta1_x_xxxz_xyy_0,   \
                             ta1_x_xxxz_xyy_1,   \
                             ta1_x_xxxz_xyyy_0,  \
                             ta1_x_xxxz_xyyy_1,  \
                             ta1_x_xxxz_xyyz_0,  \
                             ta1_x_xxxz_xyyz_1,  \
                             ta1_x_xxxz_xyz_0,   \
                             ta1_x_xxxz_xyz_1,   \
                             ta1_x_xxxz_xyzz_0,  \
                             ta1_x_xxxz_xyzz_1,  \
                             ta1_x_xxxz_xzz_0,   \
                             ta1_x_xxxz_xzz_1,   \
                             ta1_x_xxxz_xzzz_0,  \
                             ta1_x_xxxz_xzzz_1,  \
                             ta1_x_xxxz_yyyy_0,  \
                             ta1_x_xxxz_yyyy_1,  \
                             ta1_x_xxxzz_xxxx_0, \
                             ta1_x_xxxzz_xxxy_0, \
                             ta1_x_xxxzz_xxxz_0, \
                             ta1_x_xxxzz_xxyy_0, \
                             ta1_x_xxxzz_xxyz_0, \
                             ta1_x_xxxzz_xxzz_0, \
                             ta1_x_xxxzz_xyyy_0, \
                             ta1_x_xxxzz_xyyz_0, \
                             ta1_x_xxxzz_xyzz_0, \
                             ta1_x_xxxzz_xzzz_0, \
                             ta1_x_xxxzz_yyyy_0, \
                             ta1_x_xxxzz_yyyz_0, \
                             ta1_x_xxxzz_yyzz_0, \
                             ta1_x_xxxzz_yzzz_0, \
                             ta1_x_xxxzz_zzzz_0, \
                             ta1_x_xxzz_yyyz_0,  \
                             ta1_x_xxzz_yyyz_1,  \
                             ta1_x_xxzz_yyzz_0,  \
                             ta1_x_xxzz_yyzz_1,  \
                             ta1_x_xxzz_yzzz_0,  \
                             ta1_x_xxzz_yzzz_1,  \
                             ta1_x_xxzz_zzzz_0,  \
                             ta1_x_xxzz_zzzz_1,  \
                             ta1_x_xzz_yyyz_0,   \
                             ta1_x_xzz_yyyz_1,   \
                             ta1_x_xzz_yyzz_0,   \
                             ta1_x_xzz_yyzz_1,   \
                             ta1_x_xzz_yzzz_0,   \
                             ta1_x_xzz_yzzz_1,   \
                             ta1_x_xzz_zzzz_0,   \
                             ta1_x_xzz_zzzz_1,   \
                             ta_xxzz_yyyz_1,     \
                             ta_xxzz_yyzz_1,     \
                             ta_xxzz_yzzz_1,     \
                             ta_xxzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzz_xxxx_0[i] =
            ta1_x_xxx_xxxx_0[i] * fe_0 - ta1_x_xxx_xxxx_1[i] * fe_0 + ta1_x_xxxz_xxxx_0[i] * pa_z[i] - ta1_x_xxxz_xxxx_1[i] * pc_z[i];

        ta1_x_xxxzz_xxxy_0[i] =
            ta1_x_xxx_xxxy_0[i] * fe_0 - ta1_x_xxx_xxxy_1[i] * fe_0 + ta1_x_xxxz_xxxy_0[i] * pa_z[i] - ta1_x_xxxz_xxxy_1[i] * pc_z[i];

        ta1_x_xxxzz_xxxz_0[i] = ta1_x_xxx_xxxz_0[i] * fe_0 - ta1_x_xxx_xxxz_1[i] * fe_0 + ta1_x_xxxz_xxx_0[i] * fe_0 - ta1_x_xxxz_xxx_1[i] * fe_0 +
                                ta1_x_xxxz_xxxz_0[i] * pa_z[i] - ta1_x_xxxz_xxxz_1[i] * pc_z[i];

        ta1_x_xxxzz_xxyy_0[i] =
            ta1_x_xxx_xxyy_0[i] * fe_0 - ta1_x_xxx_xxyy_1[i] * fe_0 + ta1_x_xxxz_xxyy_0[i] * pa_z[i] - ta1_x_xxxz_xxyy_1[i] * pc_z[i];

        ta1_x_xxxzz_xxyz_0[i] = ta1_x_xxx_xxyz_0[i] * fe_0 - ta1_x_xxx_xxyz_1[i] * fe_0 + ta1_x_xxxz_xxy_0[i] * fe_0 - ta1_x_xxxz_xxy_1[i] * fe_0 +
                                ta1_x_xxxz_xxyz_0[i] * pa_z[i] - ta1_x_xxxz_xxyz_1[i] * pc_z[i];

        ta1_x_xxxzz_xxzz_0[i] = ta1_x_xxx_xxzz_0[i] * fe_0 - ta1_x_xxx_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxxz_xxz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxz_xxz_1[i] * fe_0 + ta1_x_xxxz_xxzz_0[i] * pa_z[i] - ta1_x_xxxz_xxzz_1[i] * pc_z[i];

        ta1_x_xxxzz_xyyy_0[i] =
            ta1_x_xxx_xyyy_0[i] * fe_0 - ta1_x_xxx_xyyy_1[i] * fe_0 + ta1_x_xxxz_xyyy_0[i] * pa_z[i] - ta1_x_xxxz_xyyy_1[i] * pc_z[i];

        ta1_x_xxxzz_xyyz_0[i] = ta1_x_xxx_xyyz_0[i] * fe_0 - ta1_x_xxx_xyyz_1[i] * fe_0 + ta1_x_xxxz_xyy_0[i] * fe_0 - ta1_x_xxxz_xyy_1[i] * fe_0 +
                                ta1_x_xxxz_xyyz_0[i] * pa_z[i] - ta1_x_xxxz_xyyz_1[i] * pc_z[i];

        ta1_x_xxxzz_xyzz_0[i] = ta1_x_xxx_xyzz_0[i] * fe_0 - ta1_x_xxx_xyzz_1[i] * fe_0 + 2.0 * ta1_x_xxxz_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxz_xyz_1[i] * fe_0 + ta1_x_xxxz_xyzz_0[i] * pa_z[i] - ta1_x_xxxz_xyzz_1[i] * pc_z[i];

        ta1_x_xxxzz_xzzz_0[i] = ta1_x_xxx_xzzz_0[i] * fe_0 - ta1_x_xxx_xzzz_1[i] * fe_0 + 3.0 * ta1_x_xxxz_xzz_0[i] * fe_0 -
                                3.0 * ta1_x_xxxz_xzz_1[i] * fe_0 + ta1_x_xxxz_xzzz_0[i] * pa_z[i] - ta1_x_xxxz_xzzz_1[i] * pc_z[i];

        ta1_x_xxxzz_yyyy_0[i] =
            ta1_x_xxx_yyyy_0[i] * fe_0 - ta1_x_xxx_yyyy_1[i] * fe_0 + ta1_x_xxxz_yyyy_0[i] * pa_z[i] - ta1_x_xxxz_yyyy_1[i] * pc_z[i];

        ta1_x_xxxzz_yyyz_0[i] = 2.0 * ta1_x_xzz_yyyz_0[i] * fe_0 - 2.0 * ta1_x_xzz_yyyz_1[i] * fe_0 + ta_xxzz_yyyz_1[i] +
                                ta1_x_xxzz_yyyz_0[i] * pa_x[i] - ta1_x_xxzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxxzz_yyzz_0[i] = 2.0 * ta1_x_xzz_yyzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_yyzz_1[i] * fe_0 + ta_xxzz_yyzz_1[i] +
                                ta1_x_xxzz_yyzz_0[i] * pa_x[i] - ta1_x_xxzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxxzz_yzzz_0[i] = 2.0 * ta1_x_xzz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_yzzz_1[i] * fe_0 + ta_xxzz_yzzz_1[i] +
                                ta1_x_xxzz_yzzz_0[i] * pa_x[i] - ta1_x_xxzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxxzz_zzzz_0[i] = 2.0 * ta1_x_xzz_zzzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_zzzz_1[i] * fe_0 + ta_xxzz_zzzz_1[i] +
                                ta1_x_xxzz_zzzz_0[i] * pa_x[i] - ta1_x_xxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : HG

    auto ta1_x_xxyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 90);

    auto ta1_x_xxyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 91);

    auto ta1_x_xxyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 92);

    auto ta1_x_xxyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 93);

    auto ta1_x_xxyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 94);

    auto ta1_x_xxyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 95);

    auto ta1_x_xxyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 96);

    auto ta1_x_xxyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 97);

    auto ta1_x_xxyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 98);

    auto ta1_x_xxyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 99);

    auto ta1_x_xxyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 100);

    auto ta1_x_xxyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 101);

    auto ta1_x_xxyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 102);

    auto ta1_x_xxyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 103);

    auto ta1_x_xxyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 104);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxy_xxxx_0,   \
                             ta1_x_xxy_xxxx_1,   \
                             ta1_x_xxy_xxxy_0,   \
                             ta1_x_xxy_xxxy_1,   \
                             ta1_x_xxy_xxxz_0,   \
                             ta1_x_xxy_xxxz_1,   \
                             ta1_x_xxy_xxyy_0,   \
                             ta1_x_xxy_xxyy_1,   \
                             ta1_x_xxy_xxyz_0,   \
                             ta1_x_xxy_xxyz_1,   \
                             ta1_x_xxy_xxzz_0,   \
                             ta1_x_xxy_xxzz_1,   \
                             ta1_x_xxy_xyyy_0,   \
                             ta1_x_xxy_xyyy_1,   \
                             ta1_x_xxy_xyyz_0,   \
                             ta1_x_xxy_xyyz_1,   \
                             ta1_x_xxy_xyzz_0,   \
                             ta1_x_xxy_xyzz_1,   \
                             ta1_x_xxy_xzzz_0,   \
                             ta1_x_xxy_xzzz_1,   \
                             ta1_x_xxy_zzzz_0,   \
                             ta1_x_xxy_zzzz_1,   \
                             ta1_x_xxyy_xxx_0,   \
                             ta1_x_xxyy_xxx_1,   \
                             ta1_x_xxyy_xxxx_0,  \
                             ta1_x_xxyy_xxxx_1,  \
                             ta1_x_xxyy_xxxy_0,  \
                             ta1_x_xxyy_xxxy_1,  \
                             ta1_x_xxyy_xxxz_0,  \
                             ta1_x_xxyy_xxxz_1,  \
                             ta1_x_xxyy_xxy_0,   \
                             ta1_x_xxyy_xxy_1,   \
                             ta1_x_xxyy_xxyy_0,  \
                             ta1_x_xxyy_xxyy_1,  \
                             ta1_x_xxyy_xxyz_0,  \
                             ta1_x_xxyy_xxyz_1,  \
                             ta1_x_xxyy_xxz_0,   \
                             ta1_x_xxyy_xxz_1,   \
                             ta1_x_xxyy_xxzz_0,  \
                             ta1_x_xxyy_xxzz_1,  \
                             ta1_x_xxyy_xyy_0,   \
                             ta1_x_xxyy_xyy_1,   \
                             ta1_x_xxyy_xyyy_0,  \
                             ta1_x_xxyy_xyyy_1,  \
                             ta1_x_xxyy_xyyz_0,  \
                             ta1_x_xxyy_xyyz_1,  \
                             ta1_x_xxyy_xyz_0,   \
                             ta1_x_xxyy_xyz_1,   \
                             ta1_x_xxyy_xyzz_0,  \
                             ta1_x_xxyy_xyzz_1,  \
                             ta1_x_xxyy_xzz_0,   \
                             ta1_x_xxyy_xzz_1,   \
                             ta1_x_xxyy_xzzz_0,  \
                             ta1_x_xxyy_xzzz_1,  \
                             ta1_x_xxyy_zzzz_0,  \
                             ta1_x_xxyy_zzzz_1,  \
                             ta1_x_xxyyy_xxxx_0, \
                             ta1_x_xxyyy_xxxy_0, \
                             ta1_x_xxyyy_xxxz_0, \
                             ta1_x_xxyyy_xxyy_0, \
                             ta1_x_xxyyy_xxyz_0, \
                             ta1_x_xxyyy_xxzz_0, \
                             ta1_x_xxyyy_xyyy_0, \
                             ta1_x_xxyyy_xyyz_0, \
                             ta1_x_xxyyy_xyzz_0, \
                             ta1_x_xxyyy_xzzz_0, \
                             ta1_x_xxyyy_yyyy_0, \
                             ta1_x_xxyyy_yyyz_0, \
                             ta1_x_xxyyy_yyzz_0, \
                             ta1_x_xxyyy_yzzz_0, \
                             ta1_x_xxyyy_zzzz_0, \
                             ta1_x_xyyy_yyyy_0,  \
                             ta1_x_xyyy_yyyy_1,  \
                             ta1_x_xyyy_yyyz_0,  \
                             ta1_x_xyyy_yyyz_1,  \
                             ta1_x_xyyy_yyzz_0,  \
                             ta1_x_xyyy_yyzz_1,  \
                             ta1_x_xyyy_yzzz_0,  \
                             ta1_x_xyyy_yzzz_1,  \
                             ta1_x_yyy_yyyy_0,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyz_0,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyzz_0,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yzzz_0,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta_xyyy_yyyy_1,     \
                             ta_xyyy_yyyz_1,     \
                             ta_xyyy_yyzz_1,     \
                             ta_xyyy_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyy_xxxx_0[i] =
            2.0 * ta1_x_xxy_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxxx_1[i] * fe_0 + ta1_x_xxyy_xxxx_0[i] * pa_y[i] - ta1_x_xxyy_xxxx_1[i] * pc_y[i];

        ta1_x_xxyyy_xxxy_0[i] = 2.0 * ta1_x_xxy_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxxy_1[i] * fe_0 + ta1_x_xxyy_xxx_0[i] * fe_0 -
                                ta1_x_xxyy_xxx_1[i] * fe_0 + ta1_x_xxyy_xxxy_0[i] * pa_y[i] - ta1_x_xxyy_xxxy_1[i] * pc_y[i];

        ta1_x_xxyyy_xxxz_0[i] =
            2.0 * ta1_x_xxy_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxxz_1[i] * fe_0 + ta1_x_xxyy_xxxz_0[i] * pa_y[i] - ta1_x_xxyy_xxxz_1[i] * pc_y[i];

        ta1_x_xxyyy_xxyy_0[i] = 2.0 * ta1_x_xxy_xxyy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxyy_xxy_0[i] * fe_0 -
                                2.0 * ta1_x_xxyy_xxy_1[i] * fe_0 + ta1_x_xxyy_xxyy_0[i] * pa_y[i] - ta1_x_xxyy_xxyy_1[i] * pc_y[i];

        ta1_x_xxyyy_xxyz_0[i] = 2.0 * ta1_x_xxy_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxyz_1[i] * fe_0 + ta1_x_xxyy_xxz_0[i] * fe_0 -
                                ta1_x_xxyy_xxz_1[i] * fe_0 + ta1_x_xxyy_xxyz_0[i] * pa_y[i] - ta1_x_xxyy_xxyz_1[i] * pc_y[i];

        ta1_x_xxyyy_xxzz_0[i] =
            2.0 * ta1_x_xxy_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxzz_1[i] * fe_0 + ta1_x_xxyy_xxzz_0[i] * pa_y[i] - ta1_x_xxyy_xxzz_1[i] * pc_y[i];

        ta1_x_xxyyy_xyyy_0[i] = 2.0 * ta1_x_xxy_xyyy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xyyy_1[i] * fe_0 + 3.0 * ta1_x_xxyy_xyy_0[i] * fe_0 -
                                3.0 * ta1_x_xxyy_xyy_1[i] * fe_0 + ta1_x_xxyy_xyyy_0[i] * pa_y[i] - ta1_x_xxyy_xyyy_1[i] * pc_y[i];

        ta1_x_xxyyy_xyyz_0[i] = 2.0 * ta1_x_xxy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxyy_xyz_1[i] * fe_0 + ta1_x_xxyy_xyyz_0[i] * pa_y[i] - ta1_x_xxyy_xyyz_1[i] * pc_y[i];

        ta1_x_xxyyy_xyzz_0[i] = 2.0 * ta1_x_xxy_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xyzz_1[i] * fe_0 + ta1_x_xxyy_xzz_0[i] * fe_0 -
                                ta1_x_xxyy_xzz_1[i] * fe_0 + ta1_x_xxyy_xyzz_0[i] * pa_y[i] - ta1_x_xxyy_xyzz_1[i] * pc_y[i];

        ta1_x_xxyyy_xzzz_0[i] =
            2.0 * ta1_x_xxy_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xzzz_1[i] * fe_0 + ta1_x_xxyy_xzzz_0[i] * pa_y[i] - ta1_x_xxyy_xzzz_1[i] * pc_y[i];

        ta1_x_xxyyy_yyyy_0[i] = ta1_x_yyy_yyyy_0[i] * fe_0 - ta1_x_yyy_yyyy_1[i] * fe_0 + ta_xyyy_yyyy_1[i] + ta1_x_xyyy_yyyy_0[i] * pa_x[i] -
                                ta1_x_xyyy_yyyy_1[i] * pc_x[i];

        ta1_x_xxyyy_yyyz_0[i] = ta1_x_yyy_yyyz_0[i] * fe_0 - ta1_x_yyy_yyyz_1[i] * fe_0 + ta_xyyy_yyyz_1[i] + ta1_x_xyyy_yyyz_0[i] * pa_x[i] -
                                ta1_x_xyyy_yyyz_1[i] * pc_x[i];

        ta1_x_xxyyy_yyzz_0[i] = ta1_x_yyy_yyzz_0[i] * fe_0 - ta1_x_yyy_yyzz_1[i] * fe_0 + ta_xyyy_yyzz_1[i] + ta1_x_xyyy_yyzz_0[i] * pa_x[i] -
                                ta1_x_xyyy_yyzz_1[i] * pc_x[i];

        ta1_x_xxyyy_yzzz_0[i] = ta1_x_yyy_yzzz_0[i] * fe_0 - ta1_x_yyy_yzzz_1[i] * fe_0 + ta_xyyy_yzzz_1[i] + ta1_x_xyyy_yzzz_0[i] * pa_x[i] -
                                ta1_x_xyyy_yzzz_1[i] * pc_x[i];

        ta1_x_xxyyy_zzzz_0[i] =
            2.0 * ta1_x_xxy_zzzz_0[i] * fe_0 - 2.0 * ta1_x_xxy_zzzz_1[i] * fe_0 + ta1_x_xxyy_zzzz_0[i] * pa_y[i] - ta1_x_xxyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 105-120 components of targeted buffer : HG

    auto ta1_x_xxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 105);

    auto ta1_x_xxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 106);

    auto ta1_x_xxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 107);

    auto ta1_x_xxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 108);

    auto ta1_x_xxyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 109);

    auto ta1_x_xxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 110);

    auto ta1_x_xxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 111);

    auto ta1_x_xxyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 112);

    auto ta1_x_xxyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 113);

    auto ta1_x_xxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 114);

    auto ta1_x_xxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 115);

    auto ta1_x_xxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 116);

    auto ta1_x_xxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 117);

    auto ta1_x_xxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 118);

    auto ta1_x_xxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 119);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxyy_xxxx_0,  \
                             ta1_x_xxyy_xxxx_1,  \
                             ta1_x_xxyy_xxxy_0,  \
                             ta1_x_xxyy_xxxy_1,  \
                             ta1_x_xxyy_xxy_0,   \
                             ta1_x_xxyy_xxy_1,   \
                             ta1_x_xxyy_xxyy_0,  \
                             ta1_x_xxyy_xxyy_1,  \
                             ta1_x_xxyy_xxyz_0,  \
                             ta1_x_xxyy_xxyz_1,  \
                             ta1_x_xxyy_xyy_0,   \
                             ta1_x_xxyy_xyy_1,   \
                             ta1_x_xxyy_xyyy_0,  \
                             ta1_x_xxyy_xyyy_1,  \
                             ta1_x_xxyy_xyyz_0,  \
                             ta1_x_xxyy_xyyz_1,  \
                             ta1_x_xxyy_xyz_0,   \
                             ta1_x_xxyy_xyz_1,   \
                             ta1_x_xxyy_xyzz_0,  \
                             ta1_x_xxyy_xyzz_1,  \
                             ta1_x_xxyy_yyy_0,   \
                             ta1_x_xxyy_yyy_1,   \
                             ta1_x_xxyy_yyyy_0,  \
                             ta1_x_xxyy_yyyy_1,  \
                             ta1_x_xxyy_yyyz_0,  \
                             ta1_x_xxyy_yyyz_1,  \
                             ta1_x_xxyy_yyz_0,   \
                             ta1_x_xxyy_yyz_1,   \
                             ta1_x_xxyy_yyzz_0,  \
                             ta1_x_xxyy_yyzz_1,  \
                             ta1_x_xxyy_yzz_0,   \
                             ta1_x_xxyy_yzz_1,   \
                             ta1_x_xxyy_yzzz_0,  \
                             ta1_x_xxyy_yzzz_1,  \
                             ta1_x_xxyyz_xxxx_0, \
                             ta1_x_xxyyz_xxxy_0, \
                             ta1_x_xxyyz_xxxz_0, \
                             ta1_x_xxyyz_xxyy_0, \
                             ta1_x_xxyyz_xxyz_0, \
                             ta1_x_xxyyz_xxzz_0, \
                             ta1_x_xxyyz_xyyy_0, \
                             ta1_x_xxyyz_xyyz_0, \
                             ta1_x_xxyyz_xyzz_0, \
                             ta1_x_xxyyz_xzzz_0, \
                             ta1_x_xxyyz_yyyy_0, \
                             ta1_x_xxyyz_yyyz_0, \
                             ta1_x_xxyyz_yyzz_0, \
                             ta1_x_xxyyz_yzzz_0, \
                             ta1_x_xxyyz_zzzz_0, \
                             ta1_x_xxyz_xxxz_0,  \
                             ta1_x_xxyz_xxxz_1,  \
                             ta1_x_xxyz_xxzz_0,  \
                             ta1_x_xxyz_xxzz_1,  \
                             ta1_x_xxyz_xzzz_0,  \
                             ta1_x_xxyz_xzzz_1,  \
                             ta1_x_xxyz_zzzz_0,  \
                             ta1_x_xxyz_zzzz_1,  \
                             ta1_x_xxz_xxxz_0,   \
                             ta1_x_xxz_xxxz_1,   \
                             ta1_x_xxz_xxzz_0,   \
                             ta1_x_xxz_xxzz_1,   \
                             ta1_x_xxz_xzzz_0,   \
                             ta1_x_xxz_xzzz_1,   \
                             ta1_x_xxz_zzzz_0,   \
                             ta1_x_xxz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyz_xxxx_0[i] = ta1_x_xxyy_xxxx_0[i] * pa_z[i] - ta1_x_xxyy_xxxx_1[i] * pc_z[i];

        ta1_x_xxyyz_xxxy_0[i] = ta1_x_xxyy_xxxy_0[i] * pa_z[i] - ta1_x_xxyy_xxxy_1[i] * pc_z[i];

        ta1_x_xxyyz_xxxz_0[i] =
            ta1_x_xxz_xxxz_0[i] * fe_0 - ta1_x_xxz_xxxz_1[i] * fe_0 + ta1_x_xxyz_xxxz_0[i] * pa_y[i] - ta1_x_xxyz_xxxz_1[i] * pc_y[i];

        ta1_x_xxyyz_xxyy_0[i] = ta1_x_xxyy_xxyy_0[i] * pa_z[i] - ta1_x_xxyy_xxyy_1[i] * pc_z[i];

        ta1_x_xxyyz_xxyz_0[i] =
            ta1_x_xxyy_xxy_0[i] * fe_0 - ta1_x_xxyy_xxy_1[i] * fe_0 + ta1_x_xxyy_xxyz_0[i] * pa_z[i] - ta1_x_xxyy_xxyz_1[i] * pc_z[i];

        ta1_x_xxyyz_xxzz_0[i] =
            ta1_x_xxz_xxzz_0[i] * fe_0 - ta1_x_xxz_xxzz_1[i] * fe_0 + ta1_x_xxyz_xxzz_0[i] * pa_y[i] - ta1_x_xxyz_xxzz_1[i] * pc_y[i];

        ta1_x_xxyyz_xyyy_0[i] = ta1_x_xxyy_xyyy_0[i] * pa_z[i] - ta1_x_xxyy_xyyy_1[i] * pc_z[i];

        ta1_x_xxyyz_xyyz_0[i] =
            ta1_x_xxyy_xyy_0[i] * fe_0 - ta1_x_xxyy_xyy_1[i] * fe_0 + ta1_x_xxyy_xyyz_0[i] * pa_z[i] - ta1_x_xxyy_xyyz_1[i] * pc_z[i];

        ta1_x_xxyyz_xyzz_0[i] =
            2.0 * ta1_x_xxyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxyy_xyz_1[i] * fe_0 + ta1_x_xxyy_xyzz_0[i] * pa_z[i] - ta1_x_xxyy_xyzz_1[i] * pc_z[i];

        ta1_x_xxyyz_xzzz_0[i] =
            ta1_x_xxz_xzzz_0[i] * fe_0 - ta1_x_xxz_xzzz_1[i] * fe_0 + ta1_x_xxyz_xzzz_0[i] * pa_y[i] - ta1_x_xxyz_xzzz_1[i] * pc_y[i];

        ta1_x_xxyyz_yyyy_0[i] = ta1_x_xxyy_yyyy_0[i] * pa_z[i] - ta1_x_xxyy_yyyy_1[i] * pc_z[i];

        ta1_x_xxyyz_yyyz_0[i] =
            ta1_x_xxyy_yyy_0[i] * fe_0 - ta1_x_xxyy_yyy_1[i] * fe_0 + ta1_x_xxyy_yyyz_0[i] * pa_z[i] - ta1_x_xxyy_yyyz_1[i] * pc_z[i];

        ta1_x_xxyyz_yyzz_0[i] =
            2.0 * ta1_x_xxyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_xxyy_yyz_1[i] * fe_0 + ta1_x_xxyy_yyzz_0[i] * pa_z[i] - ta1_x_xxyy_yyzz_1[i] * pc_z[i];

        ta1_x_xxyyz_yzzz_0[i] =
            3.0 * ta1_x_xxyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yzz_1[i] * fe_0 + ta1_x_xxyy_yzzz_0[i] * pa_z[i] - ta1_x_xxyy_yzzz_1[i] * pc_z[i];

        ta1_x_xxyyz_zzzz_0[i] =
            ta1_x_xxz_zzzz_0[i] * fe_0 - ta1_x_xxz_zzzz_1[i] * fe_0 + ta1_x_xxyz_zzzz_0[i] * pa_y[i] - ta1_x_xxyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : HG

    auto ta1_x_xxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 120);

    auto ta1_x_xxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 121);

    auto ta1_x_xxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 122);

    auto ta1_x_xxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 123);

    auto ta1_x_xxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 124);

    auto ta1_x_xxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 125);

    auto ta1_x_xxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 126);

    auto ta1_x_xxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 127);

    auto ta1_x_xxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 128);

    auto ta1_x_xxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 129);

    auto ta1_x_xxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 130);

    auto ta1_x_xxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 131);

    auto ta1_x_xxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 132);

    auto ta1_x_xxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 133);

    auto ta1_x_xxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 134);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xxyzz_xxxx_0, \
                             ta1_x_xxyzz_xxxy_0, \
                             ta1_x_xxyzz_xxxz_0, \
                             ta1_x_xxyzz_xxyy_0, \
                             ta1_x_xxyzz_xxyz_0, \
                             ta1_x_xxyzz_xxzz_0, \
                             ta1_x_xxyzz_xyyy_0, \
                             ta1_x_xxyzz_xyyz_0, \
                             ta1_x_xxyzz_xyzz_0, \
                             ta1_x_xxyzz_xzzz_0, \
                             ta1_x_xxyzz_yyyy_0, \
                             ta1_x_xxyzz_yyyz_0, \
                             ta1_x_xxyzz_yyzz_0, \
                             ta1_x_xxyzz_yzzz_0, \
                             ta1_x_xxyzz_zzzz_0, \
                             ta1_x_xxzz_xxx_0,   \
                             ta1_x_xxzz_xxx_1,   \
                             ta1_x_xxzz_xxxx_0,  \
                             ta1_x_xxzz_xxxx_1,  \
                             ta1_x_xxzz_xxxy_0,  \
                             ta1_x_xxzz_xxxy_1,  \
                             ta1_x_xxzz_xxxz_0,  \
                             ta1_x_xxzz_xxxz_1,  \
                             ta1_x_xxzz_xxy_0,   \
                             ta1_x_xxzz_xxy_1,   \
                             ta1_x_xxzz_xxyy_0,  \
                             ta1_x_xxzz_xxyy_1,  \
                             ta1_x_xxzz_xxyz_0,  \
                             ta1_x_xxzz_xxyz_1,  \
                             ta1_x_xxzz_xxz_0,   \
                             ta1_x_xxzz_xxz_1,   \
                             ta1_x_xxzz_xxzz_0,  \
                             ta1_x_xxzz_xxzz_1,  \
                             ta1_x_xxzz_xyy_0,   \
                             ta1_x_xxzz_xyy_1,   \
                             ta1_x_xxzz_xyyy_0,  \
                             ta1_x_xxzz_xyyy_1,  \
                             ta1_x_xxzz_xyyz_0,  \
                             ta1_x_xxzz_xyyz_1,  \
                             ta1_x_xxzz_xyz_0,   \
                             ta1_x_xxzz_xyz_1,   \
                             ta1_x_xxzz_xyzz_0,  \
                             ta1_x_xxzz_xyzz_1,  \
                             ta1_x_xxzz_xzz_0,   \
                             ta1_x_xxzz_xzz_1,   \
                             ta1_x_xxzz_xzzz_0,  \
                             ta1_x_xxzz_xzzz_1,  \
                             ta1_x_xxzz_yyy_0,   \
                             ta1_x_xxzz_yyy_1,   \
                             ta1_x_xxzz_yyyy_0,  \
                             ta1_x_xxzz_yyyy_1,  \
                             ta1_x_xxzz_yyyz_0,  \
                             ta1_x_xxzz_yyyz_1,  \
                             ta1_x_xxzz_yyz_0,   \
                             ta1_x_xxzz_yyz_1,   \
                             ta1_x_xxzz_yyzz_0,  \
                             ta1_x_xxzz_yyzz_1,  \
                             ta1_x_xxzz_yzz_0,   \
                             ta1_x_xxzz_yzz_1,   \
                             ta1_x_xxzz_yzzz_0,  \
                             ta1_x_xxzz_yzzz_1,  \
                             ta1_x_xxzz_zzz_0,   \
                             ta1_x_xxzz_zzz_1,   \
                             ta1_x_xxzz_zzzz_0,  \
                             ta1_x_xxzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzz_xxxx_0[i] = ta1_x_xxzz_xxxx_0[i] * pa_y[i] - ta1_x_xxzz_xxxx_1[i] * pc_y[i];

        ta1_x_xxyzz_xxxy_0[i] =
            ta1_x_xxzz_xxx_0[i] * fe_0 - ta1_x_xxzz_xxx_1[i] * fe_0 + ta1_x_xxzz_xxxy_0[i] * pa_y[i] - ta1_x_xxzz_xxxy_1[i] * pc_y[i];

        ta1_x_xxyzz_xxxz_0[i] = ta1_x_xxzz_xxxz_0[i] * pa_y[i] - ta1_x_xxzz_xxxz_1[i] * pc_y[i];

        ta1_x_xxyzz_xxyy_0[i] =
            2.0 * ta1_x_xxzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxzz_xxy_1[i] * fe_0 + ta1_x_xxzz_xxyy_0[i] * pa_y[i] - ta1_x_xxzz_xxyy_1[i] * pc_y[i];

        ta1_x_xxyzz_xxyz_0[i] =
            ta1_x_xxzz_xxz_0[i] * fe_0 - ta1_x_xxzz_xxz_1[i] * fe_0 + ta1_x_xxzz_xxyz_0[i] * pa_y[i] - ta1_x_xxzz_xxyz_1[i] * pc_y[i];

        ta1_x_xxyzz_xxzz_0[i] = ta1_x_xxzz_xxzz_0[i] * pa_y[i] - ta1_x_xxzz_xxzz_1[i] * pc_y[i];

        ta1_x_xxyzz_xyyy_0[i] =
            3.0 * ta1_x_xxzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xyy_1[i] * fe_0 + ta1_x_xxzz_xyyy_0[i] * pa_y[i] - ta1_x_xxzz_xyyy_1[i] * pc_y[i];

        ta1_x_xxyzz_xyyz_0[i] =
            2.0 * ta1_x_xxzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxzz_xyz_1[i] * fe_0 + ta1_x_xxzz_xyyz_0[i] * pa_y[i] - ta1_x_xxzz_xyyz_1[i] * pc_y[i];

        ta1_x_xxyzz_xyzz_0[i] =
            ta1_x_xxzz_xzz_0[i] * fe_0 - ta1_x_xxzz_xzz_1[i] * fe_0 + ta1_x_xxzz_xyzz_0[i] * pa_y[i] - ta1_x_xxzz_xyzz_1[i] * pc_y[i];

        ta1_x_xxyzz_xzzz_0[i] = ta1_x_xxzz_xzzz_0[i] * pa_y[i] - ta1_x_xxzz_xzzz_1[i] * pc_y[i];

        ta1_x_xxyzz_yyyy_0[i] =
            4.0 * ta1_x_xxzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxzz_yyy_1[i] * fe_0 + ta1_x_xxzz_yyyy_0[i] * pa_y[i] - ta1_x_xxzz_yyyy_1[i] * pc_y[i];

        ta1_x_xxyzz_yyyz_0[i] =
            3.0 * ta1_x_xxzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yyz_1[i] * fe_0 + ta1_x_xxzz_yyyz_0[i] * pa_y[i] - ta1_x_xxzz_yyyz_1[i] * pc_y[i];

        ta1_x_xxyzz_yyzz_0[i] =
            2.0 * ta1_x_xxzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxzz_yzz_1[i] * fe_0 + ta1_x_xxzz_yyzz_0[i] * pa_y[i] - ta1_x_xxzz_yyzz_1[i] * pc_y[i];

        ta1_x_xxyzz_yzzz_0[i] =
            ta1_x_xxzz_zzz_0[i] * fe_0 - ta1_x_xxzz_zzz_1[i] * fe_0 + ta1_x_xxzz_yzzz_0[i] * pa_y[i] - ta1_x_xxzz_yzzz_1[i] * pc_y[i];

        ta1_x_xxyzz_zzzz_0[i] = ta1_x_xxzz_zzzz_0[i] * pa_y[i] - ta1_x_xxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : HG

    auto ta1_x_xxzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 135);

    auto ta1_x_xxzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 136);

    auto ta1_x_xxzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 137);

    auto ta1_x_xxzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 138);

    auto ta1_x_xxzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 139);

    auto ta1_x_xxzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 140);

    auto ta1_x_xxzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 141);

    auto ta1_x_xxzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 142);

    auto ta1_x_xxzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 143);

    auto ta1_x_xxzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 144);

    auto ta1_x_xxzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 145);

    auto ta1_x_xxzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 146);

    auto ta1_x_xxzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 147);

    auto ta1_x_xxzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 148);

    auto ta1_x_xxzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 149);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxz_xxxx_0,   \
                             ta1_x_xxz_xxxx_1,   \
                             ta1_x_xxz_xxxy_0,   \
                             ta1_x_xxz_xxxy_1,   \
                             ta1_x_xxz_xxxz_0,   \
                             ta1_x_xxz_xxxz_1,   \
                             ta1_x_xxz_xxyy_0,   \
                             ta1_x_xxz_xxyy_1,   \
                             ta1_x_xxz_xxyz_0,   \
                             ta1_x_xxz_xxyz_1,   \
                             ta1_x_xxz_xxzz_0,   \
                             ta1_x_xxz_xxzz_1,   \
                             ta1_x_xxz_xyyy_0,   \
                             ta1_x_xxz_xyyy_1,   \
                             ta1_x_xxz_xyyz_0,   \
                             ta1_x_xxz_xyyz_1,   \
                             ta1_x_xxz_xyzz_0,   \
                             ta1_x_xxz_xyzz_1,   \
                             ta1_x_xxz_xzzz_0,   \
                             ta1_x_xxz_xzzz_1,   \
                             ta1_x_xxz_yyyy_0,   \
                             ta1_x_xxz_yyyy_1,   \
                             ta1_x_xxzz_xxx_0,   \
                             ta1_x_xxzz_xxx_1,   \
                             ta1_x_xxzz_xxxx_0,  \
                             ta1_x_xxzz_xxxx_1,  \
                             ta1_x_xxzz_xxxy_0,  \
                             ta1_x_xxzz_xxxy_1,  \
                             ta1_x_xxzz_xxxz_0,  \
                             ta1_x_xxzz_xxxz_1,  \
                             ta1_x_xxzz_xxy_0,   \
                             ta1_x_xxzz_xxy_1,   \
                             ta1_x_xxzz_xxyy_0,  \
                             ta1_x_xxzz_xxyy_1,  \
                             ta1_x_xxzz_xxyz_0,  \
                             ta1_x_xxzz_xxyz_1,  \
                             ta1_x_xxzz_xxz_0,   \
                             ta1_x_xxzz_xxz_1,   \
                             ta1_x_xxzz_xxzz_0,  \
                             ta1_x_xxzz_xxzz_1,  \
                             ta1_x_xxzz_xyy_0,   \
                             ta1_x_xxzz_xyy_1,   \
                             ta1_x_xxzz_xyyy_0,  \
                             ta1_x_xxzz_xyyy_1,  \
                             ta1_x_xxzz_xyyz_0,  \
                             ta1_x_xxzz_xyyz_1,  \
                             ta1_x_xxzz_xyz_0,   \
                             ta1_x_xxzz_xyz_1,   \
                             ta1_x_xxzz_xyzz_0,  \
                             ta1_x_xxzz_xyzz_1,  \
                             ta1_x_xxzz_xzz_0,   \
                             ta1_x_xxzz_xzz_1,   \
                             ta1_x_xxzz_xzzz_0,  \
                             ta1_x_xxzz_xzzz_1,  \
                             ta1_x_xxzz_yyyy_0,  \
                             ta1_x_xxzz_yyyy_1,  \
                             ta1_x_xxzzz_xxxx_0, \
                             ta1_x_xxzzz_xxxy_0, \
                             ta1_x_xxzzz_xxxz_0, \
                             ta1_x_xxzzz_xxyy_0, \
                             ta1_x_xxzzz_xxyz_0, \
                             ta1_x_xxzzz_xxzz_0, \
                             ta1_x_xxzzz_xyyy_0, \
                             ta1_x_xxzzz_xyyz_0, \
                             ta1_x_xxzzz_xyzz_0, \
                             ta1_x_xxzzz_xzzz_0, \
                             ta1_x_xxzzz_yyyy_0, \
                             ta1_x_xxzzz_yyyz_0, \
                             ta1_x_xxzzz_yyzz_0, \
                             ta1_x_xxzzz_yzzz_0, \
                             ta1_x_xxzzz_zzzz_0, \
                             ta1_x_xzzz_yyyz_0,  \
                             ta1_x_xzzz_yyyz_1,  \
                             ta1_x_xzzz_yyzz_0,  \
                             ta1_x_xzzz_yyzz_1,  \
                             ta1_x_xzzz_yzzz_0,  \
                             ta1_x_xzzz_yzzz_1,  \
                             ta1_x_xzzz_zzzz_0,  \
                             ta1_x_xzzz_zzzz_1,  \
                             ta1_x_zzz_yyyz_0,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyzz_0,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yzzz_0,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_zzzz_0,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta_xzzz_yyyz_1,     \
                             ta_xzzz_yyzz_1,     \
                             ta_xzzz_yzzz_1,     \
                             ta_xzzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzz_xxxx_0[i] =
            2.0 * ta1_x_xxz_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxxx_1[i] * fe_0 + ta1_x_xxzz_xxxx_0[i] * pa_z[i] - ta1_x_xxzz_xxxx_1[i] * pc_z[i];

        ta1_x_xxzzz_xxxy_0[i] =
            2.0 * ta1_x_xxz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxxy_1[i] * fe_0 + ta1_x_xxzz_xxxy_0[i] * pa_z[i] - ta1_x_xxzz_xxxy_1[i] * pc_z[i];

        ta1_x_xxzzz_xxxz_0[i] = 2.0 * ta1_x_xxz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxxz_1[i] * fe_0 + ta1_x_xxzz_xxx_0[i] * fe_0 -
                                ta1_x_xxzz_xxx_1[i] * fe_0 + ta1_x_xxzz_xxxz_0[i] * pa_z[i] - ta1_x_xxzz_xxxz_1[i] * pc_z[i];

        ta1_x_xxzzz_xxyy_0[i] =
            2.0 * ta1_x_xxz_xxyy_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxyy_1[i] * fe_0 + ta1_x_xxzz_xxyy_0[i] * pa_z[i] - ta1_x_xxzz_xxyy_1[i] * pc_z[i];

        ta1_x_xxzzz_xxyz_0[i] = 2.0 * ta1_x_xxz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxyz_1[i] * fe_0 + ta1_x_xxzz_xxy_0[i] * fe_0 -
                                ta1_x_xxzz_xxy_1[i] * fe_0 + ta1_x_xxzz_xxyz_0[i] * pa_z[i] - ta1_x_xxzz_xxyz_1[i] * pc_z[i];

        ta1_x_xxzzz_xxzz_0[i] = 2.0 * ta1_x_xxz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxzz_xxz_0[i] * fe_0 -
                                2.0 * ta1_x_xxzz_xxz_1[i] * fe_0 + ta1_x_xxzz_xxzz_0[i] * pa_z[i] - ta1_x_xxzz_xxzz_1[i] * pc_z[i];

        ta1_x_xxzzz_xyyy_0[i] =
            2.0 * ta1_x_xxz_xyyy_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyyy_1[i] * fe_0 + ta1_x_xxzz_xyyy_0[i] * pa_z[i] - ta1_x_xxzz_xyyy_1[i] * pc_z[i];

        ta1_x_xxzzz_xyyz_0[i] = 2.0 * ta1_x_xxz_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyyz_1[i] * fe_0 + ta1_x_xxzz_xyy_0[i] * fe_0 -
                                ta1_x_xxzz_xyy_1[i] * fe_0 + ta1_x_xxzz_xyyz_0[i] * pa_z[i] - ta1_x_xxzz_xyyz_1[i] * pc_z[i];

        ta1_x_xxzzz_xyzz_0[i] = 2.0 * ta1_x_xxz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyzz_1[i] * fe_0 + 2.0 * ta1_x_xxzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxzz_xyz_1[i] * fe_0 + ta1_x_xxzz_xyzz_0[i] * pa_z[i] - ta1_x_xxzz_xyzz_1[i] * pc_z[i];

        ta1_x_xxzzz_xzzz_0[i] = 2.0 * ta1_x_xxz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xzzz_1[i] * fe_0 + 3.0 * ta1_x_xxzz_xzz_0[i] * fe_0 -
                                3.0 * ta1_x_xxzz_xzz_1[i] * fe_0 + ta1_x_xxzz_xzzz_0[i] * pa_z[i] - ta1_x_xxzz_xzzz_1[i] * pc_z[i];

        ta1_x_xxzzz_yyyy_0[i] =
            2.0 * ta1_x_xxz_yyyy_0[i] * fe_0 - 2.0 * ta1_x_xxz_yyyy_1[i] * fe_0 + ta1_x_xxzz_yyyy_0[i] * pa_z[i] - ta1_x_xxzz_yyyy_1[i] * pc_z[i];

        ta1_x_xxzzz_yyyz_0[i] = ta1_x_zzz_yyyz_0[i] * fe_0 - ta1_x_zzz_yyyz_1[i] * fe_0 + ta_xzzz_yyyz_1[i] + ta1_x_xzzz_yyyz_0[i] * pa_x[i] -
                                ta1_x_xzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxzzz_yyzz_0[i] = ta1_x_zzz_yyzz_0[i] * fe_0 - ta1_x_zzz_yyzz_1[i] * fe_0 + ta_xzzz_yyzz_1[i] + ta1_x_xzzz_yyzz_0[i] * pa_x[i] -
                                ta1_x_xzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxzzz_yzzz_0[i] = ta1_x_zzz_yzzz_0[i] * fe_0 - ta1_x_zzz_yzzz_1[i] * fe_0 + ta_xzzz_yzzz_1[i] + ta1_x_xzzz_yzzz_0[i] * pa_x[i] -
                                ta1_x_xzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxzzz_zzzz_0[i] = ta1_x_zzz_zzzz_0[i] * fe_0 - ta1_x_zzz_zzzz_1[i] * fe_0 + ta_xzzz_zzzz_1[i] + ta1_x_xzzz_zzzz_0[i] * pa_x[i] -
                                ta1_x_xzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : HG

    auto ta1_x_xyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 150);

    auto ta1_x_xyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 151);

    auto ta1_x_xyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 152);

    auto ta1_x_xyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 153);

    auto ta1_x_xyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 154);

    auto ta1_x_xyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 155);

    auto ta1_x_xyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 156);

    auto ta1_x_xyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 157);

    auto ta1_x_xyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 158);

    auto ta1_x_xyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 159);

    auto ta1_x_xyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 160);

    auto ta1_x_xyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 161);

    auto ta1_x_xyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 162);

    auto ta1_x_xyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 163);

    auto ta1_x_xyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 164);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xyy_xxxx_0,   \
                             ta1_x_xyy_xxxx_1,   \
                             ta1_x_xyy_xxxz_0,   \
                             ta1_x_xyy_xxxz_1,   \
                             ta1_x_xyy_xxzz_0,   \
                             ta1_x_xyy_xxzz_1,   \
                             ta1_x_xyy_xzzz_0,   \
                             ta1_x_xyy_xzzz_1,   \
                             ta1_x_xyyy_xxxx_0,  \
                             ta1_x_xyyy_xxxx_1,  \
                             ta1_x_xyyy_xxxz_0,  \
                             ta1_x_xyyy_xxxz_1,  \
                             ta1_x_xyyy_xxzz_0,  \
                             ta1_x_xyyy_xxzz_1,  \
                             ta1_x_xyyy_xzzz_0,  \
                             ta1_x_xyyy_xzzz_1,  \
                             ta1_x_xyyyy_xxxx_0, \
                             ta1_x_xyyyy_xxxy_0, \
                             ta1_x_xyyyy_xxxz_0, \
                             ta1_x_xyyyy_xxyy_0, \
                             ta1_x_xyyyy_xxyz_0, \
                             ta1_x_xyyyy_xxzz_0, \
                             ta1_x_xyyyy_xyyy_0, \
                             ta1_x_xyyyy_xyyz_0, \
                             ta1_x_xyyyy_xyzz_0, \
                             ta1_x_xyyyy_xzzz_0, \
                             ta1_x_xyyyy_yyyy_0, \
                             ta1_x_xyyyy_yyyz_0, \
                             ta1_x_xyyyy_yyzz_0, \
                             ta1_x_xyyyy_yzzz_0, \
                             ta1_x_xyyyy_zzzz_0, \
                             ta1_x_yyyy_xxxy_0,  \
                             ta1_x_yyyy_xxxy_1,  \
                             ta1_x_yyyy_xxy_0,   \
                             ta1_x_yyyy_xxy_1,   \
                             ta1_x_yyyy_xxyy_0,  \
                             ta1_x_yyyy_xxyy_1,  \
                             ta1_x_yyyy_xxyz_0,  \
                             ta1_x_yyyy_xxyz_1,  \
                             ta1_x_yyyy_xyy_0,   \
                             ta1_x_yyyy_xyy_1,   \
                             ta1_x_yyyy_xyyy_0,  \
                             ta1_x_yyyy_xyyy_1,  \
                             ta1_x_yyyy_xyyz_0,  \
                             ta1_x_yyyy_xyyz_1,  \
                             ta1_x_yyyy_xyz_0,   \
                             ta1_x_yyyy_xyz_1,   \
                             ta1_x_yyyy_xyzz_0,  \
                             ta1_x_yyyy_xyzz_1,  \
                             ta1_x_yyyy_yyy_0,   \
                             ta1_x_yyyy_yyy_1,   \
                             ta1_x_yyyy_yyyy_0,  \
                             ta1_x_yyyy_yyyy_1,  \
                             ta1_x_yyyy_yyyz_0,  \
                             ta1_x_yyyy_yyyz_1,  \
                             ta1_x_yyyy_yyz_0,   \
                             ta1_x_yyyy_yyz_1,   \
                             ta1_x_yyyy_yyzz_0,  \
                             ta1_x_yyyy_yyzz_1,  \
                             ta1_x_yyyy_yzz_0,   \
                             ta1_x_yyyy_yzz_1,   \
                             ta1_x_yyyy_yzzz_0,  \
                             ta1_x_yyyy_yzzz_1,  \
                             ta1_x_yyyy_zzzz_0,  \
                             ta1_x_yyyy_zzzz_1,  \
                             ta_yyyy_xxxy_1,     \
                             ta_yyyy_xxyy_1,     \
                             ta_yyyy_xxyz_1,     \
                             ta_yyyy_xyyy_1,     \
                             ta_yyyy_xyyz_1,     \
                             ta_yyyy_xyzz_1,     \
                             ta_yyyy_yyyy_1,     \
                             ta_yyyy_yyyz_1,     \
                             ta_yyyy_yyzz_1,     \
                             ta_yyyy_yzzz_1,     \
                             ta_yyyy_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyy_xxxx_0[i] =
            3.0 * ta1_x_xyy_xxxx_0[i] * fe_0 - 3.0 * ta1_x_xyy_xxxx_1[i] * fe_0 + ta1_x_xyyy_xxxx_0[i] * pa_y[i] - ta1_x_xyyy_xxxx_1[i] * pc_y[i];

        ta1_x_xyyyy_xxxy_0[i] = 3.0 * ta1_x_yyyy_xxy_0[i] * fe_0 - 3.0 * ta1_x_yyyy_xxy_1[i] * fe_0 + ta_yyyy_xxxy_1[i] +
                                ta1_x_yyyy_xxxy_0[i] * pa_x[i] - ta1_x_yyyy_xxxy_1[i] * pc_x[i];

        ta1_x_xyyyy_xxxz_0[i] =
            3.0 * ta1_x_xyy_xxxz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xxxz_1[i] * fe_0 + ta1_x_xyyy_xxxz_0[i] * pa_y[i] - ta1_x_xyyy_xxxz_1[i] * pc_y[i];

        ta1_x_xyyyy_xxyy_0[i] = 2.0 * ta1_x_yyyy_xyy_0[i] * fe_0 - 2.0 * ta1_x_yyyy_xyy_1[i] * fe_0 + ta_yyyy_xxyy_1[i] +
                                ta1_x_yyyy_xxyy_0[i] * pa_x[i] - ta1_x_yyyy_xxyy_1[i] * pc_x[i];

        ta1_x_xyyyy_xxyz_0[i] = 2.0 * ta1_x_yyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyy_xyz_1[i] * fe_0 + ta_yyyy_xxyz_1[i] +
                                ta1_x_yyyy_xxyz_0[i] * pa_x[i] - ta1_x_yyyy_xxyz_1[i] * pc_x[i];

        ta1_x_xyyyy_xxzz_0[i] =
            3.0 * ta1_x_xyy_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xxzz_1[i] * fe_0 + ta1_x_xyyy_xxzz_0[i] * pa_y[i] - ta1_x_xyyy_xxzz_1[i] * pc_y[i];

        ta1_x_xyyyy_xyyy_0[i] = ta1_x_yyyy_yyy_0[i] * fe_0 - ta1_x_yyyy_yyy_1[i] * fe_0 + ta_yyyy_xyyy_1[i] + ta1_x_yyyy_xyyy_0[i] * pa_x[i] -
                                ta1_x_yyyy_xyyy_1[i] * pc_x[i];

        ta1_x_xyyyy_xyyz_0[i] = ta1_x_yyyy_yyz_0[i] * fe_0 - ta1_x_yyyy_yyz_1[i] * fe_0 + ta_yyyy_xyyz_1[i] + ta1_x_yyyy_xyyz_0[i] * pa_x[i] -
                                ta1_x_yyyy_xyyz_1[i] * pc_x[i];

        ta1_x_xyyyy_xyzz_0[i] = ta1_x_yyyy_yzz_0[i] * fe_0 - ta1_x_yyyy_yzz_1[i] * fe_0 + ta_yyyy_xyzz_1[i] + ta1_x_yyyy_xyzz_0[i] * pa_x[i] -
                                ta1_x_yyyy_xyzz_1[i] * pc_x[i];

        ta1_x_xyyyy_xzzz_0[i] =
            3.0 * ta1_x_xyy_xzzz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xzzz_1[i] * fe_0 + ta1_x_xyyy_xzzz_0[i] * pa_y[i] - ta1_x_xyyy_xzzz_1[i] * pc_y[i];

        ta1_x_xyyyy_yyyy_0[i] = ta_yyyy_yyyy_1[i] + ta1_x_yyyy_yyyy_0[i] * pa_x[i] - ta1_x_yyyy_yyyy_1[i] * pc_x[i];

        ta1_x_xyyyy_yyyz_0[i] = ta_yyyy_yyyz_1[i] + ta1_x_yyyy_yyyz_0[i] * pa_x[i] - ta1_x_yyyy_yyyz_1[i] * pc_x[i];

        ta1_x_xyyyy_yyzz_0[i] = ta_yyyy_yyzz_1[i] + ta1_x_yyyy_yyzz_0[i] * pa_x[i] - ta1_x_yyyy_yyzz_1[i] * pc_x[i];

        ta1_x_xyyyy_yzzz_0[i] = ta_yyyy_yzzz_1[i] + ta1_x_yyyy_yzzz_0[i] * pa_x[i] - ta1_x_yyyy_yzzz_1[i] * pc_x[i];

        ta1_x_xyyyy_zzzz_0[i] = ta_yyyy_zzzz_1[i] + ta1_x_yyyy_zzzz_0[i] * pa_x[i] - ta1_x_yyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 165-180 components of targeted buffer : HG

    auto ta1_x_xyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 165);

    auto ta1_x_xyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 166);

    auto ta1_x_xyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 167);

    auto ta1_x_xyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 168);

    auto ta1_x_xyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 169);

    auto ta1_x_xyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 170);

    auto ta1_x_xyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 171);

    auto ta1_x_xyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 172);

    auto ta1_x_xyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 173);

    auto ta1_x_xyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 174);

    auto ta1_x_xyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 175);

    auto ta1_x_xyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 176);

    auto ta1_x_xyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 177);

    auto ta1_x_xyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 178);

    auto ta1_x_xyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 179);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyyy_xxxx_0,  \
                             ta1_x_xyyy_xxxx_1,  \
                             ta1_x_xyyy_xxxy_0,  \
                             ta1_x_xyyy_xxxy_1,  \
                             ta1_x_xyyy_xxy_0,   \
                             ta1_x_xyyy_xxy_1,   \
                             ta1_x_xyyy_xxyy_0,  \
                             ta1_x_xyyy_xxyy_1,  \
                             ta1_x_xyyy_xxyz_0,  \
                             ta1_x_xyyy_xxyz_1,  \
                             ta1_x_xyyy_xyy_0,   \
                             ta1_x_xyyy_xyy_1,   \
                             ta1_x_xyyy_xyyy_0,  \
                             ta1_x_xyyy_xyyy_1,  \
                             ta1_x_xyyy_xyyz_0,  \
                             ta1_x_xyyy_xyyz_1,  \
                             ta1_x_xyyy_xyz_0,   \
                             ta1_x_xyyy_xyz_1,   \
                             ta1_x_xyyy_xyzz_0,  \
                             ta1_x_xyyy_xyzz_1,  \
                             ta1_x_xyyy_yyyy_0,  \
                             ta1_x_xyyy_yyyy_1,  \
                             ta1_x_xyyyz_xxxx_0, \
                             ta1_x_xyyyz_xxxy_0, \
                             ta1_x_xyyyz_xxxz_0, \
                             ta1_x_xyyyz_xxyy_0, \
                             ta1_x_xyyyz_xxyz_0, \
                             ta1_x_xyyyz_xxzz_0, \
                             ta1_x_xyyyz_xyyy_0, \
                             ta1_x_xyyyz_xyyz_0, \
                             ta1_x_xyyyz_xyzz_0, \
                             ta1_x_xyyyz_xzzz_0, \
                             ta1_x_xyyyz_yyyy_0, \
                             ta1_x_xyyyz_yyyz_0, \
                             ta1_x_xyyyz_yyzz_0, \
                             ta1_x_xyyyz_yzzz_0, \
                             ta1_x_xyyyz_zzzz_0, \
                             ta1_x_xyyz_xxxz_0,  \
                             ta1_x_xyyz_xxxz_1,  \
                             ta1_x_xyyz_xxzz_0,  \
                             ta1_x_xyyz_xxzz_1,  \
                             ta1_x_xyyz_xzzz_0,  \
                             ta1_x_xyyz_xzzz_1,  \
                             ta1_x_xyz_xxxz_0,   \
                             ta1_x_xyz_xxxz_1,   \
                             ta1_x_xyz_xxzz_0,   \
                             ta1_x_xyz_xxzz_1,   \
                             ta1_x_xyz_xzzz_0,   \
                             ta1_x_xyz_xzzz_1,   \
                             ta1_x_yyyz_yyyz_0,  \
                             ta1_x_yyyz_yyyz_1,  \
                             ta1_x_yyyz_yyzz_0,  \
                             ta1_x_yyyz_yyzz_1,  \
                             ta1_x_yyyz_yzzz_0,  \
                             ta1_x_yyyz_yzzz_1,  \
                             ta1_x_yyyz_zzzz_0,  \
                             ta1_x_yyyz_zzzz_1,  \
                             ta_yyyz_yyyz_1,     \
                             ta_yyyz_yyzz_1,     \
                             ta_yyyz_yzzz_1,     \
                             ta_yyyz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyz_xxxx_0[i] = ta1_x_xyyy_xxxx_0[i] * pa_z[i] - ta1_x_xyyy_xxxx_1[i] * pc_z[i];

        ta1_x_xyyyz_xxxy_0[i] = ta1_x_xyyy_xxxy_0[i] * pa_z[i] - ta1_x_xyyy_xxxy_1[i] * pc_z[i];

        ta1_x_xyyyz_xxxz_0[i] =
            2.0 * ta1_x_xyz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xyz_xxxz_1[i] * fe_0 + ta1_x_xyyz_xxxz_0[i] * pa_y[i] - ta1_x_xyyz_xxxz_1[i] * pc_y[i];

        ta1_x_xyyyz_xxyy_0[i] = ta1_x_xyyy_xxyy_0[i] * pa_z[i] - ta1_x_xyyy_xxyy_1[i] * pc_z[i];

        ta1_x_xyyyz_xxyz_0[i] =
            ta1_x_xyyy_xxy_0[i] * fe_0 - ta1_x_xyyy_xxy_1[i] * fe_0 + ta1_x_xyyy_xxyz_0[i] * pa_z[i] - ta1_x_xyyy_xxyz_1[i] * pc_z[i];

        ta1_x_xyyyz_xxzz_0[i] =
            2.0 * ta1_x_xyz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xyz_xxzz_1[i] * fe_0 + ta1_x_xyyz_xxzz_0[i] * pa_y[i] - ta1_x_xyyz_xxzz_1[i] * pc_y[i];

        ta1_x_xyyyz_xyyy_0[i] = ta1_x_xyyy_xyyy_0[i] * pa_z[i] - ta1_x_xyyy_xyyy_1[i] * pc_z[i];

        ta1_x_xyyyz_xyyz_0[i] =
            ta1_x_xyyy_xyy_0[i] * fe_0 - ta1_x_xyyy_xyy_1[i] * fe_0 + ta1_x_xyyy_xyyz_0[i] * pa_z[i] - ta1_x_xyyy_xyyz_1[i] * pc_z[i];

        ta1_x_xyyyz_xyzz_0[i] =
            2.0 * ta1_x_xyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_xyz_1[i] * fe_0 + ta1_x_xyyy_xyzz_0[i] * pa_z[i] - ta1_x_xyyy_xyzz_1[i] * pc_z[i];

        ta1_x_xyyyz_xzzz_0[i] =
            2.0 * ta1_x_xyz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xyz_xzzz_1[i] * fe_0 + ta1_x_xyyz_xzzz_0[i] * pa_y[i] - ta1_x_xyyz_xzzz_1[i] * pc_y[i];

        ta1_x_xyyyz_yyyy_0[i] = ta1_x_xyyy_yyyy_0[i] * pa_z[i] - ta1_x_xyyy_yyyy_1[i] * pc_z[i];

        ta1_x_xyyyz_yyyz_0[i] = ta_yyyz_yyyz_1[i] + ta1_x_yyyz_yyyz_0[i] * pa_x[i] - ta1_x_yyyz_yyyz_1[i] * pc_x[i];

        ta1_x_xyyyz_yyzz_0[i] = ta_yyyz_yyzz_1[i] + ta1_x_yyyz_yyzz_0[i] * pa_x[i] - ta1_x_yyyz_yyzz_1[i] * pc_x[i];

        ta1_x_xyyyz_yzzz_0[i] = ta_yyyz_yzzz_1[i] + ta1_x_yyyz_yzzz_0[i] * pa_x[i] - ta1_x_yyyz_yzzz_1[i] * pc_x[i];

        ta1_x_xyyyz_zzzz_0[i] = ta_yyyz_zzzz_1[i] + ta1_x_yyyz_zzzz_0[i] * pa_x[i] - ta1_x_yyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 180-195 components of targeted buffer : HG

    auto ta1_x_xyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 180);

    auto ta1_x_xyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 181);

    auto ta1_x_xyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 182);

    auto ta1_x_xyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 183);

    auto ta1_x_xyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 184);

    auto ta1_x_xyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 185);

    auto ta1_x_xyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 186);

    auto ta1_x_xyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 187);

    auto ta1_x_xyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 188);

    auto ta1_x_xyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 189);

    auto ta1_x_xyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 190);

    auto ta1_x_xyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 191);

    auto ta1_x_xyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 192);

    auto ta1_x_xyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 193);

    auto ta1_x_xyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 194);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyy_xxxy_0,   \
                             ta1_x_xyy_xxxy_1,   \
                             ta1_x_xyy_xxyy_0,   \
                             ta1_x_xyy_xxyy_1,   \
                             ta1_x_xyy_xyyy_0,   \
                             ta1_x_xyy_xyyy_1,   \
                             ta1_x_xyyz_xxxy_0,  \
                             ta1_x_xyyz_xxxy_1,  \
                             ta1_x_xyyz_xxyy_0,  \
                             ta1_x_xyyz_xxyy_1,  \
                             ta1_x_xyyz_xyyy_0,  \
                             ta1_x_xyyz_xyyy_1,  \
                             ta1_x_xyyzz_xxxx_0, \
                             ta1_x_xyyzz_xxxy_0, \
                             ta1_x_xyyzz_xxxz_0, \
                             ta1_x_xyyzz_xxyy_0, \
                             ta1_x_xyyzz_xxyz_0, \
                             ta1_x_xyyzz_xxzz_0, \
                             ta1_x_xyyzz_xyyy_0, \
                             ta1_x_xyyzz_xyyz_0, \
                             ta1_x_xyyzz_xyzz_0, \
                             ta1_x_xyyzz_xzzz_0, \
                             ta1_x_xyyzz_yyyy_0, \
                             ta1_x_xyyzz_yyyz_0, \
                             ta1_x_xyyzz_yyzz_0, \
                             ta1_x_xyyzz_yzzz_0, \
                             ta1_x_xyyzz_zzzz_0, \
                             ta1_x_xyzz_xxxx_0,  \
                             ta1_x_xyzz_xxxx_1,  \
                             ta1_x_xyzz_xxxz_0,  \
                             ta1_x_xyzz_xxxz_1,  \
                             ta1_x_xyzz_xxzz_0,  \
                             ta1_x_xyzz_xxzz_1,  \
                             ta1_x_xyzz_xzzz_0,  \
                             ta1_x_xyzz_xzzz_1,  \
                             ta1_x_xzz_xxxx_0,   \
                             ta1_x_xzz_xxxx_1,   \
                             ta1_x_xzz_xxxz_0,   \
                             ta1_x_xzz_xxxz_1,   \
                             ta1_x_xzz_xxzz_0,   \
                             ta1_x_xzz_xxzz_1,   \
                             ta1_x_xzz_xzzz_0,   \
                             ta1_x_xzz_xzzz_1,   \
                             ta1_x_yyzz_xxyz_0,  \
                             ta1_x_yyzz_xxyz_1,  \
                             ta1_x_yyzz_xyyz_0,  \
                             ta1_x_yyzz_xyyz_1,  \
                             ta1_x_yyzz_xyz_0,   \
                             ta1_x_yyzz_xyz_1,   \
                             ta1_x_yyzz_xyzz_0,  \
                             ta1_x_yyzz_xyzz_1,  \
                             ta1_x_yyzz_yyyy_0,  \
                             ta1_x_yyzz_yyyy_1,  \
                             ta1_x_yyzz_yyyz_0,  \
                             ta1_x_yyzz_yyyz_1,  \
                             ta1_x_yyzz_yyz_0,   \
                             ta1_x_yyzz_yyz_1,   \
                             ta1_x_yyzz_yyzz_0,  \
                             ta1_x_yyzz_yyzz_1,  \
                             ta1_x_yyzz_yzz_0,   \
                             ta1_x_yyzz_yzz_1,   \
                             ta1_x_yyzz_yzzz_0,  \
                             ta1_x_yyzz_yzzz_1,  \
                             ta1_x_yyzz_zzzz_0,  \
                             ta1_x_yyzz_zzzz_1,  \
                             ta_yyzz_xxyz_1,     \
                             ta_yyzz_xyyz_1,     \
                             ta_yyzz_xyzz_1,     \
                             ta_yyzz_yyyy_1,     \
                             ta_yyzz_yyyz_1,     \
                             ta_yyzz_yyzz_1,     \
                             ta_yyzz_yzzz_1,     \
                             ta_yyzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzz_xxxx_0[i] =
            ta1_x_xzz_xxxx_0[i] * fe_0 - ta1_x_xzz_xxxx_1[i] * fe_0 + ta1_x_xyzz_xxxx_0[i] * pa_y[i] - ta1_x_xyzz_xxxx_1[i] * pc_y[i];

        ta1_x_xyyzz_xxxy_0[i] =
            ta1_x_xyy_xxxy_0[i] * fe_0 - ta1_x_xyy_xxxy_1[i] * fe_0 + ta1_x_xyyz_xxxy_0[i] * pa_z[i] - ta1_x_xyyz_xxxy_1[i] * pc_z[i];

        ta1_x_xyyzz_xxxz_0[i] =
            ta1_x_xzz_xxxz_0[i] * fe_0 - ta1_x_xzz_xxxz_1[i] * fe_0 + ta1_x_xyzz_xxxz_0[i] * pa_y[i] - ta1_x_xyzz_xxxz_1[i] * pc_y[i];

        ta1_x_xyyzz_xxyy_0[i] =
            ta1_x_xyy_xxyy_0[i] * fe_0 - ta1_x_xyy_xxyy_1[i] * fe_0 + ta1_x_xyyz_xxyy_0[i] * pa_z[i] - ta1_x_xyyz_xxyy_1[i] * pc_z[i];

        ta1_x_xyyzz_xxyz_0[i] = 2.0 * ta1_x_yyzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyzz_xyz_1[i] * fe_0 + ta_yyzz_xxyz_1[i] +
                                ta1_x_yyzz_xxyz_0[i] * pa_x[i] - ta1_x_yyzz_xxyz_1[i] * pc_x[i];

        ta1_x_xyyzz_xxzz_0[i] =
            ta1_x_xzz_xxzz_0[i] * fe_0 - ta1_x_xzz_xxzz_1[i] * fe_0 + ta1_x_xyzz_xxzz_0[i] * pa_y[i] - ta1_x_xyzz_xxzz_1[i] * pc_y[i];

        ta1_x_xyyzz_xyyy_0[i] =
            ta1_x_xyy_xyyy_0[i] * fe_0 - ta1_x_xyy_xyyy_1[i] * fe_0 + ta1_x_xyyz_xyyy_0[i] * pa_z[i] - ta1_x_xyyz_xyyy_1[i] * pc_z[i];

        ta1_x_xyyzz_xyyz_0[i] = ta1_x_yyzz_yyz_0[i] * fe_0 - ta1_x_yyzz_yyz_1[i] * fe_0 + ta_yyzz_xyyz_1[i] + ta1_x_yyzz_xyyz_0[i] * pa_x[i] -
                                ta1_x_yyzz_xyyz_1[i] * pc_x[i];

        ta1_x_xyyzz_xyzz_0[i] = ta1_x_yyzz_yzz_0[i] * fe_0 - ta1_x_yyzz_yzz_1[i] * fe_0 + ta_yyzz_xyzz_1[i] + ta1_x_yyzz_xyzz_0[i] * pa_x[i] -
                                ta1_x_yyzz_xyzz_1[i] * pc_x[i];

        ta1_x_xyyzz_xzzz_0[i] =
            ta1_x_xzz_xzzz_0[i] * fe_0 - ta1_x_xzz_xzzz_1[i] * fe_0 + ta1_x_xyzz_xzzz_0[i] * pa_y[i] - ta1_x_xyzz_xzzz_1[i] * pc_y[i];

        ta1_x_xyyzz_yyyy_0[i] = ta_yyzz_yyyy_1[i] + ta1_x_yyzz_yyyy_0[i] * pa_x[i] - ta1_x_yyzz_yyyy_1[i] * pc_x[i];

        ta1_x_xyyzz_yyyz_0[i] = ta_yyzz_yyyz_1[i] + ta1_x_yyzz_yyyz_0[i] * pa_x[i] - ta1_x_yyzz_yyyz_1[i] * pc_x[i];

        ta1_x_xyyzz_yyzz_0[i] = ta_yyzz_yyzz_1[i] + ta1_x_yyzz_yyzz_0[i] * pa_x[i] - ta1_x_yyzz_yyzz_1[i] * pc_x[i];

        ta1_x_xyyzz_yzzz_0[i] = ta_yyzz_yzzz_1[i] + ta1_x_yyzz_yzzz_0[i] * pa_x[i] - ta1_x_yyzz_yzzz_1[i] * pc_x[i];

        ta1_x_xyyzz_zzzz_0[i] = ta_yyzz_zzzz_1[i] + ta1_x_yyzz_zzzz_0[i] * pa_x[i] - ta1_x_yyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : HG

    auto ta1_x_xyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 195);

    auto ta1_x_xyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 196);

    auto ta1_x_xyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 197);

    auto ta1_x_xyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 198);

    auto ta1_x_xyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 199);

    auto ta1_x_xyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 200);

    auto ta1_x_xyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 201);

    auto ta1_x_xyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 202);

    auto ta1_x_xyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 203);

    auto ta1_x_xyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 204);

    auto ta1_x_xyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 205);

    auto ta1_x_xyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 206);

    auto ta1_x_xyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 207);

    auto ta1_x_xyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 208);

    auto ta1_x_xyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 209);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xyzzz_xxxx_0, \
                             ta1_x_xyzzz_xxxy_0, \
                             ta1_x_xyzzz_xxxz_0, \
                             ta1_x_xyzzz_xxyy_0, \
                             ta1_x_xyzzz_xxyz_0, \
                             ta1_x_xyzzz_xxzz_0, \
                             ta1_x_xyzzz_xyyy_0, \
                             ta1_x_xyzzz_xyyz_0, \
                             ta1_x_xyzzz_xyzz_0, \
                             ta1_x_xyzzz_xzzz_0, \
                             ta1_x_xyzzz_yyyy_0, \
                             ta1_x_xyzzz_yyyz_0, \
                             ta1_x_xyzzz_yyzz_0, \
                             ta1_x_xyzzz_yzzz_0, \
                             ta1_x_xyzzz_zzzz_0, \
                             ta1_x_xzzz_xxx_0,   \
                             ta1_x_xzzz_xxx_1,   \
                             ta1_x_xzzz_xxxx_0,  \
                             ta1_x_xzzz_xxxx_1,  \
                             ta1_x_xzzz_xxxy_0,  \
                             ta1_x_xzzz_xxxy_1,  \
                             ta1_x_xzzz_xxxz_0,  \
                             ta1_x_xzzz_xxxz_1,  \
                             ta1_x_xzzz_xxy_0,   \
                             ta1_x_xzzz_xxy_1,   \
                             ta1_x_xzzz_xxyy_0,  \
                             ta1_x_xzzz_xxyy_1,  \
                             ta1_x_xzzz_xxyz_0,  \
                             ta1_x_xzzz_xxyz_1,  \
                             ta1_x_xzzz_xxz_0,   \
                             ta1_x_xzzz_xxz_1,   \
                             ta1_x_xzzz_xxzz_0,  \
                             ta1_x_xzzz_xxzz_1,  \
                             ta1_x_xzzz_xyy_0,   \
                             ta1_x_xzzz_xyy_1,   \
                             ta1_x_xzzz_xyyy_0,  \
                             ta1_x_xzzz_xyyy_1,  \
                             ta1_x_xzzz_xyyz_0,  \
                             ta1_x_xzzz_xyyz_1,  \
                             ta1_x_xzzz_xyz_0,   \
                             ta1_x_xzzz_xyz_1,   \
                             ta1_x_xzzz_xyzz_0,  \
                             ta1_x_xzzz_xyzz_1,  \
                             ta1_x_xzzz_xzz_0,   \
                             ta1_x_xzzz_xzz_1,   \
                             ta1_x_xzzz_xzzz_0,  \
                             ta1_x_xzzz_xzzz_1,  \
                             ta1_x_xzzz_zzzz_0,  \
                             ta1_x_xzzz_zzzz_1,  \
                             ta1_x_yzzz_yyyy_0,  \
                             ta1_x_yzzz_yyyy_1,  \
                             ta1_x_yzzz_yyyz_0,  \
                             ta1_x_yzzz_yyyz_1,  \
                             ta1_x_yzzz_yyzz_0,  \
                             ta1_x_yzzz_yyzz_1,  \
                             ta1_x_yzzz_yzzz_0,  \
                             ta1_x_yzzz_yzzz_1,  \
                             ta_yzzz_yyyy_1,     \
                             ta_yzzz_yyyz_1,     \
                             ta_yzzz_yyzz_1,     \
                             ta_yzzz_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzzz_xxxx_0[i] = ta1_x_xzzz_xxxx_0[i] * pa_y[i] - ta1_x_xzzz_xxxx_1[i] * pc_y[i];

        ta1_x_xyzzz_xxxy_0[i] =
            ta1_x_xzzz_xxx_0[i] * fe_0 - ta1_x_xzzz_xxx_1[i] * fe_0 + ta1_x_xzzz_xxxy_0[i] * pa_y[i] - ta1_x_xzzz_xxxy_1[i] * pc_y[i];

        ta1_x_xyzzz_xxxz_0[i] = ta1_x_xzzz_xxxz_0[i] * pa_y[i] - ta1_x_xzzz_xxxz_1[i] * pc_y[i];

        ta1_x_xyzzz_xxyy_0[i] =
            2.0 * ta1_x_xzzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xzzz_xxy_1[i] * fe_0 + ta1_x_xzzz_xxyy_0[i] * pa_y[i] - ta1_x_xzzz_xxyy_1[i] * pc_y[i];

        ta1_x_xyzzz_xxyz_0[i] =
            ta1_x_xzzz_xxz_0[i] * fe_0 - ta1_x_xzzz_xxz_1[i] * fe_0 + ta1_x_xzzz_xxyz_0[i] * pa_y[i] - ta1_x_xzzz_xxyz_1[i] * pc_y[i];

        ta1_x_xyzzz_xxzz_0[i] = ta1_x_xzzz_xxzz_0[i] * pa_y[i] - ta1_x_xzzz_xxzz_1[i] * pc_y[i];

        ta1_x_xyzzz_xyyy_0[i] =
            3.0 * ta1_x_xzzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xzzz_xyy_1[i] * fe_0 + ta1_x_xzzz_xyyy_0[i] * pa_y[i] - ta1_x_xzzz_xyyy_1[i] * pc_y[i];

        ta1_x_xyzzz_xyyz_0[i] =
            2.0 * ta1_x_xzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_xyz_1[i] * fe_0 + ta1_x_xzzz_xyyz_0[i] * pa_y[i] - ta1_x_xzzz_xyyz_1[i] * pc_y[i];

        ta1_x_xyzzz_xyzz_0[i] =
            ta1_x_xzzz_xzz_0[i] * fe_0 - ta1_x_xzzz_xzz_1[i] * fe_0 + ta1_x_xzzz_xyzz_0[i] * pa_y[i] - ta1_x_xzzz_xyzz_1[i] * pc_y[i];

        ta1_x_xyzzz_xzzz_0[i] = ta1_x_xzzz_xzzz_0[i] * pa_y[i] - ta1_x_xzzz_xzzz_1[i] * pc_y[i];

        ta1_x_xyzzz_yyyy_0[i] = ta_yzzz_yyyy_1[i] + ta1_x_yzzz_yyyy_0[i] * pa_x[i] - ta1_x_yzzz_yyyy_1[i] * pc_x[i];

        ta1_x_xyzzz_yyyz_0[i] = ta_yzzz_yyyz_1[i] + ta1_x_yzzz_yyyz_0[i] * pa_x[i] - ta1_x_yzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xyzzz_yyzz_0[i] = ta_yzzz_yyzz_1[i] + ta1_x_yzzz_yyzz_0[i] * pa_x[i] - ta1_x_yzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xyzzz_yzzz_0[i] = ta_yzzz_yzzz_1[i] + ta1_x_yzzz_yzzz_0[i] * pa_x[i] - ta1_x_yzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xyzzz_zzzz_0[i] = ta1_x_xzzz_zzzz_0[i] * pa_y[i] - ta1_x_xzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : HG

    auto ta1_x_xzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 210);

    auto ta1_x_xzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 211);

    auto ta1_x_xzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 212);

    auto ta1_x_xzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 213);

    auto ta1_x_xzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 214);

    auto ta1_x_xzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 215);

    auto ta1_x_xzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 216);

    auto ta1_x_xzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 217);

    auto ta1_x_xzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 218);

    auto ta1_x_xzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 219);

    auto ta1_x_xzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 220);

    auto ta1_x_xzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 221);

    auto ta1_x_xzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 222);

    auto ta1_x_xzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 223);

    auto ta1_x_xzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 224);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xzz_xxxx_0,   \
                             ta1_x_xzz_xxxx_1,   \
                             ta1_x_xzz_xxxy_0,   \
                             ta1_x_xzz_xxxy_1,   \
                             ta1_x_xzz_xxyy_0,   \
                             ta1_x_xzz_xxyy_1,   \
                             ta1_x_xzz_xyyy_0,   \
                             ta1_x_xzz_xyyy_1,   \
                             ta1_x_xzzz_xxxx_0,  \
                             ta1_x_xzzz_xxxx_1,  \
                             ta1_x_xzzz_xxxy_0,  \
                             ta1_x_xzzz_xxxy_1,  \
                             ta1_x_xzzz_xxyy_0,  \
                             ta1_x_xzzz_xxyy_1,  \
                             ta1_x_xzzz_xyyy_0,  \
                             ta1_x_xzzz_xyyy_1,  \
                             ta1_x_xzzzz_xxxx_0, \
                             ta1_x_xzzzz_xxxy_0, \
                             ta1_x_xzzzz_xxxz_0, \
                             ta1_x_xzzzz_xxyy_0, \
                             ta1_x_xzzzz_xxyz_0, \
                             ta1_x_xzzzz_xxzz_0, \
                             ta1_x_xzzzz_xyyy_0, \
                             ta1_x_xzzzz_xyyz_0, \
                             ta1_x_xzzzz_xyzz_0, \
                             ta1_x_xzzzz_xzzz_0, \
                             ta1_x_xzzzz_yyyy_0, \
                             ta1_x_xzzzz_yyyz_0, \
                             ta1_x_xzzzz_yyzz_0, \
                             ta1_x_xzzzz_yzzz_0, \
                             ta1_x_xzzzz_zzzz_0, \
                             ta1_x_zzzz_xxxz_0,  \
                             ta1_x_zzzz_xxxz_1,  \
                             ta1_x_zzzz_xxyz_0,  \
                             ta1_x_zzzz_xxyz_1,  \
                             ta1_x_zzzz_xxz_0,   \
                             ta1_x_zzzz_xxz_1,   \
                             ta1_x_zzzz_xxzz_0,  \
                             ta1_x_zzzz_xxzz_1,  \
                             ta1_x_zzzz_xyyz_0,  \
                             ta1_x_zzzz_xyyz_1,  \
                             ta1_x_zzzz_xyz_0,   \
                             ta1_x_zzzz_xyz_1,   \
                             ta1_x_zzzz_xyzz_0,  \
                             ta1_x_zzzz_xyzz_1,  \
                             ta1_x_zzzz_xzz_0,   \
                             ta1_x_zzzz_xzz_1,   \
                             ta1_x_zzzz_xzzz_0,  \
                             ta1_x_zzzz_xzzz_1,  \
                             ta1_x_zzzz_yyyy_0,  \
                             ta1_x_zzzz_yyyy_1,  \
                             ta1_x_zzzz_yyyz_0,  \
                             ta1_x_zzzz_yyyz_1,  \
                             ta1_x_zzzz_yyz_0,   \
                             ta1_x_zzzz_yyz_1,   \
                             ta1_x_zzzz_yyzz_0,  \
                             ta1_x_zzzz_yyzz_1,  \
                             ta1_x_zzzz_yzz_0,   \
                             ta1_x_zzzz_yzz_1,   \
                             ta1_x_zzzz_yzzz_0,  \
                             ta1_x_zzzz_yzzz_1,  \
                             ta1_x_zzzz_zzz_0,   \
                             ta1_x_zzzz_zzz_1,   \
                             ta1_x_zzzz_zzzz_0,  \
                             ta1_x_zzzz_zzzz_1,  \
                             ta_zzzz_xxxz_1,     \
                             ta_zzzz_xxyz_1,     \
                             ta_zzzz_xxzz_1,     \
                             ta_zzzz_xyyz_1,     \
                             ta_zzzz_xyzz_1,     \
                             ta_zzzz_xzzz_1,     \
                             ta_zzzz_yyyy_1,     \
                             ta_zzzz_yyyz_1,     \
                             ta_zzzz_yyzz_1,     \
                             ta_zzzz_yzzz_1,     \
                             ta_zzzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzz_xxxx_0[i] =
            3.0 * ta1_x_xzz_xxxx_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxxx_1[i] * fe_0 + ta1_x_xzzz_xxxx_0[i] * pa_z[i] - ta1_x_xzzz_xxxx_1[i] * pc_z[i];

        ta1_x_xzzzz_xxxy_0[i] =
            3.0 * ta1_x_xzz_xxxy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxxy_1[i] * fe_0 + ta1_x_xzzz_xxxy_0[i] * pa_z[i] - ta1_x_xzzz_xxxy_1[i] * pc_z[i];

        ta1_x_xzzzz_xxxz_0[i] = 3.0 * ta1_x_zzzz_xxz_0[i] * fe_0 - 3.0 * ta1_x_zzzz_xxz_1[i] * fe_0 + ta_zzzz_xxxz_1[i] +
                                ta1_x_zzzz_xxxz_0[i] * pa_x[i] - ta1_x_zzzz_xxxz_1[i] * pc_x[i];

        ta1_x_xzzzz_xxyy_0[i] =
            3.0 * ta1_x_xzz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxyy_1[i] * fe_0 + ta1_x_xzzz_xxyy_0[i] * pa_z[i] - ta1_x_xzzz_xxyy_1[i] * pc_z[i];

        ta1_x_xzzzz_xxyz_0[i] = 2.0 * ta1_x_zzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xyz_1[i] * fe_0 + ta_zzzz_xxyz_1[i] +
                                ta1_x_zzzz_xxyz_0[i] * pa_x[i] - ta1_x_zzzz_xxyz_1[i] * pc_x[i];

        ta1_x_xzzzz_xxzz_0[i] = 2.0 * ta1_x_zzzz_xzz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xzz_1[i] * fe_0 + ta_zzzz_xxzz_1[i] +
                                ta1_x_zzzz_xxzz_0[i] * pa_x[i] - ta1_x_zzzz_xxzz_1[i] * pc_x[i];

        ta1_x_xzzzz_xyyy_0[i] =
            3.0 * ta1_x_xzz_xyyy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xyyy_1[i] * fe_0 + ta1_x_xzzz_xyyy_0[i] * pa_z[i] - ta1_x_xzzz_xyyy_1[i] * pc_z[i];

        ta1_x_xzzzz_xyyz_0[i] = ta1_x_zzzz_yyz_0[i] * fe_0 - ta1_x_zzzz_yyz_1[i] * fe_0 + ta_zzzz_xyyz_1[i] + ta1_x_zzzz_xyyz_0[i] * pa_x[i] -
                                ta1_x_zzzz_xyyz_1[i] * pc_x[i];

        ta1_x_xzzzz_xyzz_0[i] = ta1_x_zzzz_yzz_0[i] * fe_0 - ta1_x_zzzz_yzz_1[i] * fe_0 + ta_zzzz_xyzz_1[i] + ta1_x_zzzz_xyzz_0[i] * pa_x[i] -
                                ta1_x_zzzz_xyzz_1[i] * pc_x[i];

        ta1_x_xzzzz_xzzz_0[i] = ta1_x_zzzz_zzz_0[i] * fe_0 - ta1_x_zzzz_zzz_1[i] * fe_0 + ta_zzzz_xzzz_1[i] + ta1_x_zzzz_xzzz_0[i] * pa_x[i] -
                                ta1_x_zzzz_xzzz_1[i] * pc_x[i];

        ta1_x_xzzzz_yyyy_0[i] = ta_zzzz_yyyy_1[i] + ta1_x_zzzz_yyyy_0[i] * pa_x[i] - ta1_x_zzzz_yyyy_1[i] * pc_x[i];

        ta1_x_xzzzz_yyyz_0[i] = ta_zzzz_yyyz_1[i] + ta1_x_zzzz_yyyz_0[i] * pa_x[i] - ta1_x_zzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xzzzz_yyzz_0[i] = ta_zzzz_yyzz_1[i] + ta1_x_zzzz_yyzz_0[i] * pa_x[i] - ta1_x_zzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xzzzz_yzzz_0[i] = ta_zzzz_yzzz_1[i] + ta1_x_zzzz_yzzz_0[i] * pa_x[i] - ta1_x_zzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xzzzz_zzzz_0[i] = ta_zzzz_zzzz_1[i] + ta1_x_zzzz_zzzz_0[i] * pa_x[i] - ta1_x_zzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : HG

    auto ta1_x_yyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 225);

    auto ta1_x_yyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 226);

    auto ta1_x_yyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 227);

    auto ta1_x_yyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 228);

    auto ta1_x_yyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 229);

    auto ta1_x_yyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 230);

    auto ta1_x_yyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 231);

    auto ta1_x_yyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 232);

    auto ta1_x_yyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 233);

    auto ta1_x_yyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 234);

    auto ta1_x_yyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 235);

    auto ta1_x_yyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 236);

    auto ta1_x_yyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 237);

    auto ta1_x_yyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 238);

    auto ta1_x_yyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 239);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yyy_xxxx_0,   \
                             ta1_x_yyy_xxxx_1,   \
                             ta1_x_yyy_xxxy_0,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxxz_0,   \
                             ta1_x_yyy_xxxz_1,   \
                             ta1_x_yyy_xxyy_0,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyz_0,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xxzz_0,   \
                             ta1_x_yyy_xxzz_1,   \
                             ta1_x_yyy_xyyy_0,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyz_0,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyzz_0,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_xzzz_0,   \
                             ta1_x_yyy_xzzz_1,   \
                             ta1_x_yyy_yyyy_0,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyz_0,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyzz_0,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yzzz_0,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta1_x_yyy_zzzz_0,   \
                             ta1_x_yyy_zzzz_1,   \
                             ta1_x_yyyy_xxx_0,   \
                             ta1_x_yyyy_xxx_1,   \
                             ta1_x_yyyy_xxxx_0,  \
                             ta1_x_yyyy_xxxx_1,  \
                             ta1_x_yyyy_xxxy_0,  \
                             ta1_x_yyyy_xxxy_1,  \
                             ta1_x_yyyy_xxxz_0,  \
                             ta1_x_yyyy_xxxz_1,  \
                             ta1_x_yyyy_xxy_0,   \
                             ta1_x_yyyy_xxy_1,   \
                             ta1_x_yyyy_xxyy_0,  \
                             ta1_x_yyyy_xxyy_1,  \
                             ta1_x_yyyy_xxyz_0,  \
                             ta1_x_yyyy_xxyz_1,  \
                             ta1_x_yyyy_xxz_0,   \
                             ta1_x_yyyy_xxz_1,   \
                             ta1_x_yyyy_xxzz_0,  \
                             ta1_x_yyyy_xxzz_1,  \
                             ta1_x_yyyy_xyy_0,   \
                             ta1_x_yyyy_xyy_1,   \
                             ta1_x_yyyy_xyyy_0,  \
                             ta1_x_yyyy_xyyy_1,  \
                             ta1_x_yyyy_xyyz_0,  \
                             ta1_x_yyyy_xyyz_1,  \
                             ta1_x_yyyy_xyz_0,   \
                             ta1_x_yyyy_xyz_1,   \
                             ta1_x_yyyy_xyzz_0,  \
                             ta1_x_yyyy_xyzz_1,  \
                             ta1_x_yyyy_xzz_0,   \
                             ta1_x_yyyy_xzz_1,   \
                             ta1_x_yyyy_xzzz_0,  \
                             ta1_x_yyyy_xzzz_1,  \
                             ta1_x_yyyy_yyy_0,   \
                             ta1_x_yyyy_yyy_1,   \
                             ta1_x_yyyy_yyyy_0,  \
                             ta1_x_yyyy_yyyy_1,  \
                             ta1_x_yyyy_yyyz_0,  \
                             ta1_x_yyyy_yyyz_1,  \
                             ta1_x_yyyy_yyz_0,   \
                             ta1_x_yyyy_yyz_1,   \
                             ta1_x_yyyy_yyzz_0,  \
                             ta1_x_yyyy_yyzz_1,  \
                             ta1_x_yyyy_yzz_0,   \
                             ta1_x_yyyy_yzz_1,   \
                             ta1_x_yyyy_yzzz_0,  \
                             ta1_x_yyyy_yzzz_1,  \
                             ta1_x_yyyy_zzz_0,   \
                             ta1_x_yyyy_zzz_1,   \
                             ta1_x_yyyy_zzzz_0,  \
                             ta1_x_yyyy_zzzz_1,  \
                             ta1_x_yyyyy_xxxx_0, \
                             ta1_x_yyyyy_xxxy_0, \
                             ta1_x_yyyyy_xxxz_0, \
                             ta1_x_yyyyy_xxyy_0, \
                             ta1_x_yyyyy_xxyz_0, \
                             ta1_x_yyyyy_xxzz_0, \
                             ta1_x_yyyyy_xyyy_0, \
                             ta1_x_yyyyy_xyyz_0, \
                             ta1_x_yyyyy_xyzz_0, \
                             ta1_x_yyyyy_xzzz_0, \
                             ta1_x_yyyyy_yyyy_0, \
                             ta1_x_yyyyy_yyyz_0, \
                             ta1_x_yyyyy_yyzz_0, \
                             ta1_x_yyyyy_yzzz_0, \
                             ta1_x_yyyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyy_xxxx_0[i] =
            4.0 * ta1_x_yyy_xxxx_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxxx_1[i] * fe_0 + ta1_x_yyyy_xxxx_0[i] * pa_y[i] - ta1_x_yyyy_xxxx_1[i] * pc_y[i];

        ta1_x_yyyyy_xxxy_0[i] = 4.0 * ta1_x_yyy_xxxy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxxy_1[i] * fe_0 + ta1_x_yyyy_xxx_0[i] * fe_0 -
                                ta1_x_yyyy_xxx_1[i] * fe_0 + ta1_x_yyyy_xxxy_0[i] * pa_y[i] - ta1_x_yyyy_xxxy_1[i] * pc_y[i];

        ta1_x_yyyyy_xxxz_0[i] =
            4.0 * ta1_x_yyy_xxxz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxxz_1[i] * fe_0 + ta1_x_yyyy_xxxz_0[i] * pa_y[i] - ta1_x_yyyy_xxxz_1[i] * pc_y[i];

        ta1_x_yyyyy_xxyy_0[i] = 4.0 * ta1_x_yyy_xxyy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxyy_1[i] * fe_0 + 2.0 * ta1_x_yyyy_xxy_0[i] * fe_0 -
                                2.0 * ta1_x_yyyy_xxy_1[i] * fe_0 + ta1_x_yyyy_xxyy_0[i] * pa_y[i] - ta1_x_yyyy_xxyy_1[i] * pc_y[i];

        ta1_x_yyyyy_xxyz_0[i] = 4.0 * ta1_x_yyy_xxyz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxyz_1[i] * fe_0 + ta1_x_yyyy_xxz_0[i] * fe_0 -
                                ta1_x_yyyy_xxz_1[i] * fe_0 + ta1_x_yyyy_xxyz_0[i] * pa_y[i] - ta1_x_yyyy_xxyz_1[i] * pc_y[i];

        ta1_x_yyyyy_xxzz_0[i] =
            4.0 * ta1_x_yyy_xxzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxzz_1[i] * fe_0 + ta1_x_yyyy_xxzz_0[i] * pa_y[i] - ta1_x_yyyy_xxzz_1[i] * pc_y[i];

        ta1_x_yyyyy_xyyy_0[i] = 4.0 * ta1_x_yyy_xyyy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xyyy_1[i] * fe_0 + 3.0 * ta1_x_yyyy_xyy_0[i] * fe_0 -
                                3.0 * ta1_x_yyyy_xyy_1[i] * fe_0 + ta1_x_yyyy_xyyy_0[i] * pa_y[i] - ta1_x_yyyy_xyyy_1[i] * pc_y[i];

        ta1_x_yyyyy_xyyz_0[i] = 4.0 * ta1_x_yyy_xyyz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yyyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_yyyy_xyz_1[i] * fe_0 + ta1_x_yyyy_xyyz_0[i] * pa_y[i] - ta1_x_yyyy_xyyz_1[i] * pc_y[i];

        ta1_x_yyyyy_xyzz_0[i] = 4.0 * ta1_x_yyy_xyzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xyzz_1[i] * fe_0 + ta1_x_yyyy_xzz_0[i] * fe_0 -
                                ta1_x_yyyy_xzz_1[i] * fe_0 + ta1_x_yyyy_xyzz_0[i] * pa_y[i] - ta1_x_yyyy_xyzz_1[i] * pc_y[i];

        ta1_x_yyyyy_xzzz_0[i] =
            4.0 * ta1_x_yyy_xzzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xzzz_1[i] * fe_0 + ta1_x_yyyy_xzzz_0[i] * pa_y[i] - ta1_x_yyyy_xzzz_1[i] * pc_y[i];

        ta1_x_yyyyy_yyyy_0[i] = 4.0 * ta1_x_yyy_yyyy_0[i] * fe_0 - 4.0 * ta1_x_yyy_yyyy_1[i] * fe_0 + 4.0 * ta1_x_yyyy_yyy_0[i] * fe_0 -
                                4.0 * ta1_x_yyyy_yyy_1[i] * fe_0 + ta1_x_yyyy_yyyy_0[i] * pa_y[i] - ta1_x_yyyy_yyyy_1[i] * pc_y[i];

        ta1_x_yyyyy_yyyz_0[i] = 4.0 * ta1_x_yyy_yyyz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yyyy_yyz_0[i] * fe_0 -
                                3.0 * ta1_x_yyyy_yyz_1[i] * fe_0 + ta1_x_yyyy_yyyz_0[i] * pa_y[i] - ta1_x_yyyy_yyyz_1[i] * pc_y[i];

        ta1_x_yyyyy_yyzz_0[i] = 4.0 * ta1_x_yyy_yyzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yyyy_yzz_0[i] * fe_0 -
                                2.0 * ta1_x_yyyy_yzz_1[i] * fe_0 + ta1_x_yyyy_yyzz_0[i] * pa_y[i] - ta1_x_yyyy_yyzz_1[i] * pc_y[i];

        ta1_x_yyyyy_yzzz_0[i] = 4.0 * ta1_x_yyy_yzzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yzzz_1[i] * fe_0 + ta1_x_yyyy_zzz_0[i] * fe_0 -
                                ta1_x_yyyy_zzz_1[i] * fe_0 + ta1_x_yyyy_yzzz_0[i] * pa_y[i] - ta1_x_yyyy_yzzz_1[i] * pc_y[i];

        ta1_x_yyyyy_zzzz_0[i] =
            4.0 * ta1_x_yyy_zzzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_zzzz_1[i] * fe_0 + ta1_x_yyyy_zzzz_0[i] * pa_y[i] - ta1_x_yyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 240-255 components of targeted buffer : HG

    auto ta1_x_yyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 240);

    auto ta1_x_yyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 241);

    auto ta1_x_yyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 242);

    auto ta1_x_yyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 243);

    auto ta1_x_yyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 244);

    auto ta1_x_yyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 245);

    auto ta1_x_yyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 246);

    auto ta1_x_yyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 247);

    auto ta1_x_yyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 248);

    auto ta1_x_yyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 249);

    auto ta1_x_yyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 250);

    auto ta1_x_yyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 251);

    auto ta1_x_yyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 252);

    auto ta1_x_yyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 253);

    auto ta1_x_yyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 254);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyyy_xxxx_0,  \
                             ta1_x_yyyy_xxxx_1,  \
                             ta1_x_yyyy_xxxy_0,  \
                             ta1_x_yyyy_xxxy_1,  \
                             ta1_x_yyyy_xxy_0,   \
                             ta1_x_yyyy_xxy_1,   \
                             ta1_x_yyyy_xxyy_0,  \
                             ta1_x_yyyy_xxyy_1,  \
                             ta1_x_yyyy_xxyz_0,  \
                             ta1_x_yyyy_xxyz_1,  \
                             ta1_x_yyyy_xyy_0,   \
                             ta1_x_yyyy_xyy_1,   \
                             ta1_x_yyyy_xyyy_0,  \
                             ta1_x_yyyy_xyyy_1,  \
                             ta1_x_yyyy_xyyz_0,  \
                             ta1_x_yyyy_xyyz_1,  \
                             ta1_x_yyyy_xyz_0,   \
                             ta1_x_yyyy_xyz_1,   \
                             ta1_x_yyyy_xyzz_0,  \
                             ta1_x_yyyy_xyzz_1,  \
                             ta1_x_yyyy_yyy_0,   \
                             ta1_x_yyyy_yyy_1,   \
                             ta1_x_yyyy_yyyy_0,  \
                             ta1_x_yyyy_yyyy_1,  \
                             ta1_x_yyyy_yyyz_0,  \
                             ta1_x_yyyy_yyyz_1,  \
                             ta1_x_yyyy_yyz_0,   \
                             ta1_x_yyyy_yyz_1,   \
                             ta1_x_yyyy_yyzz_0,  \
                             ta1_x_yyyy_yyzz_1,  \
                             ta1_x_yyyy_yzz_0,   \
                             ta1_x_yyyy_yzz_1,   \
                             ta1_x_yyyy_yzzz_0,  \
                             ta1_x_yyyy_yzzz_1,  \
                             ta1_x_yyyyz_xxxx_0, \
                             ta1_x_yyyyz_xxxy_0, \
                             ta1_x_yyyyz_xxxz_0, \
                             ta1_x_yyyyz_xxyy_0, \
                             ta1_x_yyyyz_xxyz_0, \
                             ta1_x_yyyyz_xxzz_0, \
                             ta1_x_yyyyz_xyyy_0, \
                             ta1_x_yyyyz_xyyz_0, \
                             ta1_x_yyyyz_xyzz_0, \
                             ta1_x_yyyyz_xzzz_0, \
                             ta1_x_yyyyz_yyyy_0, \
                             ta1_x_yyyyz_yyyz_0, \
                             ta1_x_yyyyz_yyzz_0, \
                             ta1_x_yyyyz_yzzz_0, \
                             ta1_x_yyyyz_zzzz_0, \
                             ta1_x_yyyz_xxxz_0,  \
                             ta1_x_yyyz_xxxz_1,  \
                             ta1_x_yyyz_xxzz_0,  \
                             ta1_x_yyyz_xxzz_1,  \
                             ta1_x_yyyz_xzzz_0,  \
                             ta1_x_yyyz_xzzz_1,  \
                             ta1_x_yyyz_zzzz_0,  \
                             ta1_x_yyyz_zzzz_1,  \
                             ta1_x_yyz_xxxz_0,   \
                             ta1_x_yyz_xxxz_1,   \
                             ta1_x_yyz_xxzz_0,   \
                             ta1_x_yyz_xxzz_1,   \
                             ta1_x_yyz_xzzz_0,   \
                             ta1_x_yyz_xzzz_1,   \
                             ta1_x_yyz_zzzz_0,   \
                             ta1_x_yyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyz_xxxx_0[i] = ta1_x_yyyy_xxxx_0[i] * pa_z[i] - ta1_x_yyyy_xxxx_1[i] * pc_z[i];

        ta1_x_yyyyz_xxxy_0[i] = ta1_x_yyyy_xxxy_0[i] * pa_z[i] - ta1_x_yyyy_xxxy_1[i] * pc_z[i];

        ta1_x_yyyyz_xxxz_0[i] =
            3.0 * ta1_x_yyz_xxxz_0[i] * fe_0 - 3.0 * ta1_x_yyz_xxxz_1[i] * fe_0 + ta1_x_yyyz_xxxz_0[i] * pa_y[i] - ta1_x_yyyz_xxxz_1[i] * pc_y[i];

        ta1_x_yyyyz_xxyy_0[i] = ta1_x_yyyy_xxyy_0[i] * pa_z[i] - ta1_x_yyyy_xxyy_1[i] * pc_z[i];

        ta1_x_yyyyz_xxyz_0[i] =
            ta1_x_yyyy_xxy_0[i] * fe_0 - ta1_x_yyyy_xxy_1[i] * fe_0 + ta1_x_yyyy_xxyz_0[i] * pa_z[i] - ta1_x_yyyy_xxyz_1[i] * pc_z[i];

        ta1_x_yyyyz_xxzz_0[i] =
            3.0 * ta1_x_yyz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_yyz_xxzz_1[i] * fe_0 + ta1_x_yyyz_xxzz_0[i] * pa_y[i] - ta1_x_yyyz_xxzz_1[i] * pc_y[i];

        ta1_x_yyyyz_xyyy_0[i] = ta1_x_yyyy_xyyy_0[i] * pa_z[i] - ta1_x_yyyy_xyyy_1[i] * pc_z[i];

        ta1_x_yyyyz_xyyz_0[i] =
            ta1_x_yyyy_xyy_0[i] * fe_0 - ta1_x_yyyy_xyy_1[i] * fe_0 + ta1_x_yyyy_xyyz_0[i] * pa_z[i] - ta1_x_yyyy_xyyz_1[i] * pc_z[i];

        ta1_x_yyyyz_xyzz_0[i] =
            2.0 * ta1_x_yyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyy_xyz_1[i] * fe_0 + ta1_x_yyyy_xyzz_0[i] * pa_z[i] - ta1_x_yyyy_xyzz_1[i] * pc_z[i];

        ta1_x_yyyyz_xzzz_0[i] =
            3.0 * ta1_x_yyz_xzzz_0[i] * fe_0 - 3.0 * ta1_x_yyz_xzzz_1[i] * fe_0 + ta1_x_yyyz_xzzz_0[i] * pa_y[i] - ta1_x_yyyz_xzzz_1[i] * pc_y[i];

        ta1_x_yyyyz_yyyy_0[i] = ta1_x_yyyy_yyyy_0[i] * pa_z[i] - ta1_x_yyyy_yyyy_1[i] * pc_z[i];

        ta1_x_yyyyz_yyyz_0[i] =
            ta1_x_yyyy_yyy_0[i] * fe_0 - ta1_x_yyyy_yyy_1[i] * fe_0 + ta1_x_yyyy_yyyz_0[i] * pa_z[i] - ta1_x_yyyy_yyyz_1[i] * pc_z[i];

        ta1_x_yyyyz_yyzz_0[i] =
            2.0 * ta1_x_yyyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_yyyy_yyz_1[i] * fe_0 + ta1_x_yyyy_yyzz_0[i] * pa_z[i] - ta1_x_yyyy_yyzz_1[i] * pc_z[i];

        ta1_x_yyyyz_yzzz_0[i] =
            3.0 * ta1_x_yyyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_yyyy_yzz_1[i] * fe_0 + ta1_x_yyyy_yzzz_0[i] * pa_z[i] - ta1_x_yyyy_yzzz_1[i] * pc_z[i];

        ta1_x_yyyyz_zzzz_0[i] =
            3.0 * ta1_x_yyz_zzzz_0[i] * fe_0 - 3.0 * ta1_x_yyz_zzzz_1[i] * fe_0 + ta1_x_yyyz_zzzz_0[i] * pa_y[i] - ta1_x_yyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : HG

    auto ta1_x_yyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 255);

    auto ta1_x_yyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 256);

    auto ta1_x_yyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 257);

    auto ta1_x_yyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 258);

    auto ta1_x_yyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 259);

    auto ta1_x_yyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 260);

    auto ta1_x_yyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 261);

    auto ta1_x_yyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 262);

    auto ta1_x_yyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 263);

    auto ta1_x_yyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 264);

    auto ta1_x_yyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 265);

    auto ta1_x_yyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 266);

    auto ta1_x_yyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 267);

    auto ta1_x_yyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 268);

    auto ta1_x_yyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 269);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyy_xxxy_0,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxyy_0,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xyyy_0,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_yyyy_0,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyyz_xxxy_0,  \
                             ta1_x_yyyz_xxxy_1,  \
                             ta1_x_yyyz_xxyy_0,  \
                             ta1_x_yyyz_xxyy_1,  \
                             ta1_x_yyyz_xyyy_0,  \
                             ta1_x_yyyz_xyyy_1,  \
                             ta1_x_yyyz_yyyy_0,  \
                             ta1_x_yyyz_yyyy_1,  \
                             ta1_x_yyyzz_xxxx_0, \
                             ta1_x_yyyzz_xxxy_0, \
                             ta1_x_yyyzz_xxxz_0, \
                             ta1_x_yyyzz_xxyy_0, \
                             ta1_x_yyyzz_xxyz_0, \
                             ta1_x_yyyzz_xxzz_0, \
                             ta1_x_yyyzz_xyyy_0, \
                             ta1_x_yyyzz_xyyz_0, \
                             ta1_x_yyyzz_xyzz_0, \
                             ta1_x_yyyzz_xzzz_0, \
                             ta1_x_yyyzz_yyyy_0, \
                             ta1_x_yyyzz_yyyz_0, \
                             ta1_x_yyyzz_yyzz_0, \
                             ta1_x_yyyzz_yzzz_0, \
                             ta1_x_yyyzz_zzzz_0, \
                             ta1_x_yyzz_xxxx_0,  \
                             ta1_x_yyzz_xxxx_1,  \
                             ta1_x_yyzz_xxxz_0,  \
                             ta1_x_yyzz_xxxz_1,  \
                             ta1_x_yyzz_xxyz_0,  \
                             ta1_x_yyzz_xxyz_1,  \
                             ta1_x_yyzz_xxz_0,   \
                             ta1_x_yyzz_xxz_1,   \
                             ta1_x_yyzz_xxzz_0,  \
                             ta1_x_yyzz_xxzz_1,  \
                             ta1_x_yyzz_xyyz_0,  \
                             ta1_x_yyzz_xyyz_1,  \
                             ta1_x_yyzz_xyz_0,   \
                             ta1_x_yyzz_xyz_1,   \
                             ta1_x_yyzz_xyzz_0,  \
                             ta1_x_yyzz_xyzz_1,  \
                             ta1_x_yyzz_xzz_0,   \
                             ta1_x_yyzz_xzz_1,   \
                             ta1_x_yyzz_xzzz_0,  \
                             ta1_x_yyzz_xzzz_1,  \
                             ta1_x_yyzz_yyyz_0,  \
                             ta1_x_yyzz_yyyz_1,  \
                             ta1_x_yyzz_yyz_0,   \
                             ta1_x_yyzz_yyz_1,   \
                             ta1_x_yyzz_yyzz_0,  \
                             ta1_x_yyzz_yyzz_1,  \
                             ta1_x_yyzz_yzz_0,   \
                             ta1_x_yyzz_yzz_1,   \
                             ta1_x_yyzz_yzzz_0,  \
                             ta1_x_yyzz_yzzz_1,  \
                             ta1_x_yyzz_zzz_0,   \
                             ta1_x_yyzz_zzz_1,   \
                             ta1_x_yyzz_zzzz_0,  \
                             ta1_x_yyzz_zzzz_1,  \
                             ta1_x_yzz_xxxx_0,   \
                             ta1_x_yzz_xxxx_1,   \
                             ta1_x_yzz_xxxz_0,   \
                             ta1_x_yzz_xxxz_1,   \
                             ta1_x_yzz_xxyz_0,   \
                             ta1_x_yzz_xxyz_1,   \
                             ta1_x_yzz_xxzz_0,   \
                             ta1_x_yzz_xxzz_1,   \
                             ta1_x_yzz_xyyz_0,   \
                             ta1_x_yzz_xyyz_1,   \
                             ta1_x_yzz_xyzz_0,   \
                             ta1_x_yzz_xyzz_1,   \
                             ta1_x_yzz_xzzz_0,   \
                             ta1_x_yzz_xzzz_1,   \
                             ta1_x_yzz_yyyz_0,   \
                             ta1_x_yzz_yyyz_1,   \
                             ta1_x_yzz_yyzz_0,   \
                             ta1_x_yzz_yyzz_1,   \
                             ta1_x_yzz_yzzz_0,   \
                             ta1_x_yzz_yzzz_1,   \
                             ta1_x_yzz_zzzz_0,   \
                             ta1_x_yzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzz_xxxx_0[i] =
            2.0 * ta1_x_yzz_xxxx_0[i] * fe_0 - 2.0 * ta1_x_yzz_xxxx_1[i] * fe_0 + ta1_x_yyzz_xxxx_0[i] * pa_y[i] - ta1_x_yyzz_xxxx_1[i] * pc_y[i];

        ta1_x_yyyzz_xxxy_0[i] =
            ta1_x_yyy_xxxy_0[i] * fe_0 - ta1_x_yyy_xxxy_1[i] * fe_0 + ta1_x_yyyz_xxxy_0[i] * pa_z[i] - ta1_x_yyyz_xxxy_1[i] * pc_z[i];

        ta1_x_yyyzz_xxxz_0[i] =
            2.0 * ta1_x_yzz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xxxz_1[i] * fe_0 + ta1_x_yyzz_xxxz_0[i] * pa_y[i] - ta1_x_yyzz_xxxz_1[i] * pc_y[i];

        ta1_x_yyyzz_xxyy_0[i] =
            ta1_x_yyy_xxyy_0[i] * fe_0 - ta1_x_yyy_xxyy_1[i] * fe_0 + ta1_x_yyyz_xxyy_0[i] * pa_z[i] - ta1_x_yyyz_xxyy_1[i] * pc_z[i];

        ta1_x_yyyzz_xxyz_0[i] = 2.0 * ta1_x_yzz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xxyz_1[i] * fe_0 + ta1_x_yyzz_xxz_0[i] * fe_0 -
                                ta1_x_yyzz_xxz_1[i] * fe_0 + ta1_x_yyzz_xxyz_0[i] * pa_y[i] - ta1_x_yyzz_xxyz_1[i] * pc_y[i];

        ta1_x_yyyzz_xxzz_0[i] =
            2.0 * ta1_x_yzz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xxzz_1[i] * fe_0 + ta1_x_yyzz_xxzz_0[i] * pa_y[i] - ta1_x_yyzz_xxzz_1[i] * pc_y[i];

        ta1_x_yyyzz_xyyy_0[i] =
            ta1_x_yyy_xyyy_0[i] * fe_0 - ta1_x_yyy_xyyy_1[i] * fe_0 + ta1_x_yyyz_xyyy_0[i] * pa_z[i] - ta1_x_yyyz_xyyy_1[i] * pc_z[i];

        ta1_x_yyyzz_xyyz_0[i] = 2.0 * ta1_x_yzz_xyyz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yyzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_yyzz_xyz_1[i] * fe_0 + ta1_x_yyzz_xyyz_0[i] * pa_y[i] - ta1_x_yyzz_xyyz_1[i] * pc_y[i];

        ta1_x_yyyzz_xyzz_0[i] = 2.0 * ta1_x_yzz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xyzz_1[i] * fe_0 + ta1_x_yyzz_xzz_0[i] * fe_0 -
                                ta1_x_yyzz_xzz_1[i] * fe_0 + ta1_x_yyzz_xyzz_0[i] * pa_y[i] - ta1_x_yyzz_xyzz_1[i] * pc_y[i];

        ta1_x_yyyzz_xzzz_0[i] =
            2.0 * ta1_x_yzz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xzzz_1[i] * fe_0 + ta1_x_yyzz_xzzz_0[i] * pa_y[i] - ta1_x_yyzz_xzzz_1[i] * pc_y[i];

        ta1_x_yyyzz_yyyy_0[i] =
            ta1_x_yyy_yyyy_0[i] * fe_0 - ta1_x_yyy_yyyy_1[i] * fe_0 + ta1_x_yyyz_yyyy_0[i] * pa_z[i] - ta1_x_yyyz_yyyy_1[i] * pc_z[i];

        ta1_x_yyyzz_yyyz_0[i] = 2.0 * ta1_x_yzz_yyyz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yyzz_yyz_0[i] * fe_0 -
                                3.0 * ta1_x_yyzz_yyz_1[i] * fe_0 + ta1_x_yyzz_yyyz_0[i] * pa_y[i] - ta1_x_yyzz_yyyz_1[i] * pc_y[i];

        ta1_x_yyyzz_yyzz_0[i] = 2.0 * ta1_x_yzz_yyzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yyzz_yzz_0[i] * fe_0 -
                                2.0 * ta1_x_yyzz_yzz_1[i] * fe_0 + ta1_x_yyzz_yyzz_0[i] * pa_y[i] - ta1_x_yyzz_yyzz_1[i] * pc_y[i];

        ta1_x_yyyzz_yzzz_0[i] = 2.0 * ta1_x_yzz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yzzz_1[i] * fe_0 + ta1_x_yyzz_zzz_0[i] * fe_0 -
                                ta1_x_yyzz_zzz_1[i] * fe_0 + ta1_x_yyzz_yzzz_0[i] * pa_y[i] - ta1_x_yyzz_yzzz_1[i] * pc_y[i];

        ta1_x_yyyzz_zzzz_0[i] =
            2.0 * ta1_x_yzz_zzzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_zzzz_1[i] * fe_0 + ta1_x_yyzz_zzzz_0[i] * pa_y[i] - ta1_x_yyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 270-285 components of targeted buffer : HG

    auto ta1_x_yyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 270);

    auto ta1_x_yyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 271);

    auto ta1_x_yyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 272);

    auto ta1_x_yyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 273);

    auto ta1_x_yyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 274);

    auto ta1_x_yyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 275);

    auto ta1_x_yyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 276);

    auto ta1_x_yyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 277);

    auto ta1_x_yyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 278);

    auto ta1_x_yyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 279);

    auto ta1_x_yyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 280);

    auto ta1_x_yyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 281);

    auto ta1_x_yyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 282);

    auto ta1_x_yyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 283);

    auto ta1_x_yyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 284);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyz_xxxy_0,   \
                             ta1_x_yyz_xxxy_1,   \
                             ta1_x_yyz_xxyy_0,   \
                             ta1_x_yyz_xxyy_1,   \
                             ta1_x_yyz_xyyy_0,   \
                             ta1_x_yyz_xyyy_1,   \
                             ta1_x_yyz_yyyy_0,   \
                             ta1_x_yyz_yyyy_1,   \
                             ta1_x_yyzz_xxxy_0,  \
                             ta1_x_yyzz_xxxy_1,  \
                             ta1_x_yyzz_xxyy_0,  \
                             ta1_x_yyzz_xxyy_1,  \
                             ta1_x_yyzz_xyyy_0,  \
                             ta1_x_yyzz_xyyy_1,  \
                             ta1_x_yyzz_yyyy_0,  \
                             ta1_x_yyzz_yyyy_1,  \
                             ta1_x_yyzzz_xxxx_0, \
                             ta1_x_yyzzz_xxxy_0, \
                             ta1_x_yyzzz_xxxz_0, \
                             ta1_x_yyzzz_xxyy_0, \
                             ta1_x_yyzzz_xxyz_0, \
                             ta1_x_yyzzz_xxzz_0, \
                             ta1_x_yyzzz_xyyy_0, \
                             ta1_x_yyzzz_xyyz_0, \
                             ta1_x_yyzzz_xyzz_0, \
                             ta1_x_yyzzz_xzzz_0, \
                             ta1_x_yyzzz_yyyy_0, \
                             ta1_x_yyzzz_yyyz_0, \
                             ta1_x_yyzzz_yyzz_0, \
                             ta1_x_yyzzz_yzzz_0, \
                             ta1_x_yyzzz_zzzz_0, \
                             ta1_x_yzzz_xxxx_0,  \
                             ta1_x_yzzz_xxxx_1,  \
                             ta1_x_yzzz_xxxz_0,  \
                             ta1_x_yzzz_xxxz_1,  \
                             ta1_x_yzzz_xxyz_0,  \
                             ta1_x_yzzz_xxyz_1,  \
                             ta1_x_yzzz_xxz_0,   \
                             ta1_x_yzzz_xxz_1,   \
                             ta1_x_yzzz_xxzz_0,  \
                             ta1_x_yzzz_xxzz_1,  \
                             ta1_x_yzzz_xyyz_0,  \
                             ta1_x_yzzz_xyyz_1,  \
                             ta1_x_yzzz_xyz_0,   \
                             ta1_x_yzzz_xyz_1,   \
                             ta1_x_yzzz_xyzz_0,  \
                             ta1_x_yzzz_xyzz_1,  \
                             ta1_x_yzzz_xzz_0,   \
                             ta1_x_yzzz_xzz_1,   \
                             ta1_x_yzzz_xzzz_0,  \
                             ta1_x_yzzz_xzzz_1,  \
                             ta1_x_yzzz_yyyz_0,  \
                             ta1_x_yzzz_yyyz_1,  \
                             ta1_x_yzzz_yyz_0,   \
                             ta1_x_yzzz_yyz_1,   \
                             ta1_x_yzzz_yyzz_0,  \
                             ta1_x_yzzz_yyzz_1,  \
                             ta1_x_yzzz_yzz_0,   \
                             ta1_x_yzzz_yzz_1,   \
                             ta1_x_yzzz_yzzz_0,  \
                             ta1_x_yzzz_yzzz_1,  \
                             ta1_x_yzzz_zzz_0,   \
                             ta1_x_yzzz_zzz_1,   \
                             ta1_x_yzzz_zzzz_0,  \
                             ta1_x_yzzz_zzzz_1,  \
                             ta1_x_zzz_xxxx_0,   \
                             ta1_x_zzz_xxxx_1,   \
                             ta1_x_zzz_xxxz_0,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxyz_0,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxzz_0,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xyyz_0,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyzz_0,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xzzz_0,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_yyyz_0,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyzz_0,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yzzz_0,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_zzzz_0,   \
                             ta1_x_zzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzz_xxxx_0[i] =
            ta1_x_zzz_xxxx_0[i] * fe_0 - ta1_x_zzz_xxxx_1[i] * fe_0 + ta1_x_yzzz_xxxx_0[i] * pa_y[i] - ta1_x_yzzz_xxxx_1[i] * pc_y[i];

        ta1_x_yyzzz_xxxy_0[i] =
            2.0 * ta1_x_yyz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_yyz_xxxy_1[i] * fe_0 + ta1_x_yyzz_xxxy_0[i] * pa_z[i] - ta1_x_yyzz_xxxy_1[i] * pc_z[i];

        ta1_x_yyzzz_xxxz_0[i] =
            ta1_x_zzz_xxxz_0[i] * fe_0 - ta1_x_zzz_xxxz_1[i] * fe_0 + ta1_x_yzzz_xxxz_0[i] * pa_y[i] - ta1_x_yzzz_xxxz_1[i] * pc_y[i];

        ta1_x_yyzzz_xxyy_0[i] =
            2.0 * ta1_x_yyz_xxyy_0[i] * fe_0 - 2.0 * ta1_x_yyz_xxyy_1[i] * fe_0 + ta1_x_yyzz_xxyy_0[i] * pa_z[i] - ta1_x_yyzz_xxyy_1[i] * pc_z[i];

        ta1_x_yyzzz_xxyz_0[i] = ta1_x_zzz_xxyz_0[i] * fe_0 - ta1_x_zzz_xxyz_1[i] * fe_0 + ta1_x_yzzz_xxz_0[i] * fe_0 - ta1_x_yzzz_xxz_1[i] * fe_0 +
                                ta1_x_yzzz_xxyz_0[i] * pa_y[i] - ta1_x_yzzz_xxyz_1[i] * pc_y[i];

        ta1_x_yyzzz_xxzz_0[i] =
            ta1_x_zzz_xxzz_0[i] * fe_0 - ta1_x_zzz_xxzz_1[i] * fe_0 + ta1_x_yzzz_xxzz_0[i] * pa_y[i] - ta1_x_yzzz_xxzz_1[i] * pc_y[i];

        ta1_x_yyzzz_xyyy_0[i] =
            2.0 * ta1_x_yyz_xyyy_0[i] * fe_0 - 2.0 * ta1_x_yyz_xyyy_1[i] * fe_0 + ta1_x_yyzz_xyyy_0[i] * pa_z[i] - ta1_x_yyzz_xyyy_1[i] * pc_z[i];

        ta1_x_yyzzz_xyyz_0[i] = ta1_x_zzz_xyyz_0[i] * fe_0 - ta1_x_zzz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_yzzz_xyz_1[i] * fe_0 + ta1_x_yzzz_xyyz_0[i] * pa_y[i] - ta1_x_yzzz_xyyz_1[i] * pc_y[i];

        ta1_x_yyzzz_xyzz_0[i] = ta1_x_zzz_xyzz_0[i] * fe_0 - ta1_x_zzz_xyzz_1[i] * fe_0 + ta1_x_yzzz_xzz_0[i] * fe_0 - ta1_x_yzzz_xzz_1[i] * fe_0 +
                                ta1_x_yzzz_xyzz_0[i] * pa_y[i] - ta1_x_yzzz_xyzz_1[i] * pc_y[i];

        ta1_x_yyzzz_xzzz_0[i] =
            ta1_x_zzz_xzzz_0[i] * fe_0 - ta1_x_zzz_xzzz_1[i] * fe_0 + ta1_x_yzzz_xzzz_0[i] * pa_y[i] - ta1_x_yzzz_xzzz_1[i] * pc_y[i];

        ta1_x_yyzzz_yyyy_0[i] =
            2.0 * ta1_x_yyz_yyyy_0[i] * fe_0 - 2.0 * ta1_x_yyz_yyyy_1[i] * fe_0 + ta1_x_yyzz_yyyy_0[i] * pa_z[i] - ta1_x_yyzz_yyyy_1[i] * pc_z[i];

        ta1_x_yyzzz_yyyz_0[i] = ta1_x_zzz_yyyz_0[i] * fe_0 - ta1_x_zzz_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yzzz_yyz_0[i] * fe_0 -
                                3.0 * ta1_x_yzzz_yyz_1[i] * fe_0 + ta1_x_yzzz_yyyz_0[i] * pa_y[i] - ta1_x_yzzz_yyyz_1[i] * pc_y[i];

        ta1_x_yyzzz_yyzz_0[i] = ta1_x_zzz_yyzz_0[i] * fe_0 - ta1_x_zzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yzzz_yzz_0[i] * fe_0 -
                                2.0 * ta1_x_yzzz_yzz_1[i] * fe_0 + ta1_x_yzzz_yyzz_0[i] * pa_y[i] - ta1_x_yzzz_yyzz_1[i] * pc_y[i];

        ta1_x_yyzzz_yzzz_0[i] = ta1_x_zzz_yzzz_0[i] * fe_0 - ta1_x_zzz_yzzz_1[i] * fe_0 + ta1_x_yzzz_zzz_0[i] * fe_0 - ta1_x_yzzz_zzz_1[i] * fe_0 +
                                ta1_x_yzzz_yzzz_0[i] * pa_y[i] - ta1_x_yzzz_yzzz_1[i] * pc_y[i];

        ta1_x_yyzzz_zzzz_0[i] =
            ta1_x_zzz_zzzz_0[i] * fe_0 - ta1_x_zzz_zzzz_1[i] * fe_0 + ta1_x_yzzz_zzzz_0[i] * pa_y[i] - ta1_x_yzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 285-300 components of targeted buffer : HG

    auto ta1_x_yzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 285);

    auto ta1_x_yzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 286);

    auto ta1_x_yzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 287);

    auto ta1_x_yzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 288);

    auto ta1_x_yzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 289);

    auto ta1_x_yzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 290);

    auto ta1_x_yzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 291);

    auto ta1_x_yzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 292);

    auto ta1_x_yzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 293);

    auto ta1_x_yzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 294);

    auto ta1_x_yzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 295);

    auto ta1_x_yzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 296);

    auto ta1_x_yzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 297);

    auto ta1_x_yzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 298);

    auto ta1_x_yzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 299);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yzzzz_xxxx_0, \
                             ta1_x_yzzzz_xxxy_0, \
                             ta1_x_yzzzz_xxxz_0, \
                             ta1_x_yzzzz_xxyy_0, \
                             ta1_x_yzzzz_xxyz_0, \
                             ta1_x_yzzzz_xxzz_0, \
                             ta1_x_yzzzz_xyyy_0, \
                             ta1_x_yzzzz_xyyz_0, \
                             ta1_x_yzzzz_xyzz_0, \
                             ta1_x_yzzzz_xzzz_0, \
                             ta1_x_yzzzz_yyyy_0, \
                             ta1_x_yzzzz_yyyz_0, \
                             ta1_x_yzzzz_yyzz_0, \
                             ta1_x_yzzzz_yzzz_0, \
                             ta1_x_yzzzz_zzzz_0, \
                             ta1_x_zzzz_xxx_0,   \
                             ta1_x_zzzz_xxx_1,   \
                             ta1_x_zzzz_xxxx_0,  \
                             ta1_x_zzzz_xxxx_1,  \
                             ta1_x_zzzz_xxxy_0,  \
                             ta1_x_zzzz_xxxy_1,  \
                             ta1_x_zzzz_xxxz_0,  \
                             ta1_x_zzzz_xxxz_1,  \
                             ta1_x_zzzz_xxy_0,   \
                             ta1_x_zzzz_xxy_1,   \
                             ta1_x_zzzz_xxyy_0,  \
                             ta1_x_zzzz_xxyy_1,  \
                             ta1_x_zzzz_xxyz_0,  \
                             ta1_x_zzzz_xxyz_1,  \
                             ta1_x_zzzz_xxz_0,   \
                             ta1_x_zzzz_xxz_1,   \
                             ta1_x_zzzz_xxzz_0,  \
                             ta1_x_zzzz_xxzz_1,  \
                             ta1_x_zzzz_xyy_0,   \
                             ta1_x_zzzz_xyy_1,   \
                             ta1_x_zzzz_xyyy_0,  \
                             ta1_x_zzzz_xyyy_1,  \
                             ta1_x_zzzz_xyyz_0,  \
                             ta1_x_zzzz_xyyz_1,  \
                             ta1_x_zzzz_xyz_0,   \
                             ta1_x_zzzz_xyz_1,   \
                             ta1_x_zzzz_xyzz_0,  \
                             ta1_x_zzzz_xyzz_1,  \
                             ta1_x_zzzz_xzz_0,   \
                             ta1_x_zzzz_xzz_1,   \
                             ta1_x_zzzz_xzzz_0,  \
                             ta1_x_zzzz_xzzz_1,  \
                             ta1_x_zzzz_yyy_0,   \
                             ta1_x_zzzz_yyy_1,   \
                             ta1_x_zzzz_yyyy_0,  \
                             ta1_x_zzzz_yyyy_1,  \
                             ta1_x_zzzz_yyyz_0,  \
                             ta1_x_zzzz_yyyz_1,  \
                             ta1_x_zzzz_yyz_0,   \
                             ta1_x_zzzz_yyz_1,   \
                             ta1_x_zzzz_yyzz_0,  \
                             ta1_x_zzzz_yyzz_1,  \
                             ta1_x_zzzz_yzz_0,   \
                             ta1_x_zzzz_yzz_1,   \
                             ta1_x_zzzz_yzzz_0,  \
                             ta1_x_zzzz_yzzz_1,  \
                             ta1_x_zzzz_zzz_0,   \
                             ta1_x_zzzz_zzz_1,   \
                             ta1_x_zzzz_zzzz_0,  \
                             ta1_x_zzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzz_xxxx_0[i] = ta1_x_zzzz_xxxx_0[i] * pa_y[i] - ta1_x_zzzz_xxxx_1[i] * pc_y[i];

        ta1_x_yzzzz_xxxy_0[i] =
            ta1_x_zzzz_xxx_0[i] * fe_0 - ta1_x_zzzz_xxx_1[i] * fe_0 + ta1_x_zzzz_xxxy_0[i] * pa_y[i] - ta1_x_zzzz_xxxy_1[i] * pc_y[i];

        ta1_x_yzzzz_xxxz_0[i] = ta1_x_zzzz_xxxz_0[i] * pa_y[i] - ta1_x_zzzz_xxxz_1[i] * pc_y[i];

        ta1_x_yzzzz_xxyy_0[i] =
            2.0 * ta1_x_zzzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xxy_1[i] * fe_0 + ta1_x_zzzz_xxyy_0[i] * pa_y[i] - ta1_x_zzzz_xxyy_1[i] * pc_y[i];

        ta1_x_yzzzz_xxyz_0[i] =
            ta1_x_zzzz_xxz_0[i] * fe_0 - ta1_x_zzzz_xxz_1[i] * fe_0 + ta1_x_zzzz_xxyz_0[i] * pa_y[i] - ta1_x_zzzz_xxyz_1[i] * pc_y[i];

        ta1_x_yzzzz_xxzz_0[i] = ta1_x_zzzz_xxzz_0[i] * pa_y[i] - ta1_x_zzzz_xxzz_1[i] * pc_y[i];

        ta1_x_yzzzz_xyyy_0[i] =
            3.0 * ta1_x_zzzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_zzzz_xyy_1[i] * fe_0 + ta1_x_zzzz_xyyy_0[i] * pa_y[i] - ta1_x_zzzz_xyyy_1[i] * pc_y[i];

        ta1_x_yzzzz_xyyz_0[i] =
            2.0 * ta1_x_zzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xyz_1[i] * fe_0 + ta1_x_zzzz_xyyz_0[i] * pa_y[i] - ta1_x_zzzz_xyyz_1[i] * pc_y[i];

        ta1_x_yzzzz_xyzz_0[i] =
            ta1_x_zzzz_xzz_0[i] * fe_0 - ta1_x_zzzz_xzz_1[i] * fe_0 + ta1_x_zzzz_xyzz_0[i] * pa_y[i] - ta1_x_zzzz_xyzz_1[i] * pc_y[i];

        ta1_x_yzzzz_xzzz_0[i] = ta1_x_zzzz_xzzz_0[i] * pa_y[i] - ta1_x_zzzz_xzzz_1[i] * pc_y[i];

        ta1_x_yzzzz_yyyy_0[i] =
            4.0 * ta1_x_zzzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_zzzz_yyy_1[i] * fe_0 + ta1_x_zzzz_yyyy_0[i] * pa_y[i] - ta1_x_zzzz_yyyy_1[i] * pc_y[i];

        ta1_x_yzzzz_yyyz_0[i] =
            3.0 * ta1_x_zzzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_zzzz_yyz_1[i] * fe_0 + ta1_x_zzzz_yyyz_0[i] * pa_y[i] - ta1_x_zzzz_yyyz_1[i] * pc_y[i];

        ta1_x_yzzzz_yyzz_0[i] =
            2.0 * ta1_x_zzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_yzz_1[i] * fe_0 + ta1_x_zzzz_yyzz_0[i] * pa_y[i] - ta1_x_zzzz_yyzz_1[i] * pc_y[i];

        ta1_x_yzzzz_yzzz_0[i] =
            ta1_x_zzzz_zzz_0[i] * fe_0 - ta1_x_zzzz_zzz_1[i] * fe_0 + ta1_x_zzzz_yzzz_0[i] * pa_y[i] - ta1_x_zzzz_yzzz_1[i] * pc_y[i];

        ta1_x_yzzzz_zzzz_0[i] = ta1_x_zzzz_zzzz_0[i] * pa_y[i] - ta1_x_zzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 300-315 components of targeted buffer : HG

    auto ta1_x_zzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 300);

    auto ta1_x_zzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 301);

    auto ta1_x_zzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 302);

    auto ta1_x_zzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 303);

    auto ta1_x_zzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 304);

    auto ta1_x_zzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 305);

    auto ta1_x_zzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 306);

    auto ta1_x_zzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 307);

    auto ta1_x_zzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 308);

    auto ta1_x_zzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 309);

    auto ta1_x_zzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 310);

    auto ta1_x_zzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 311);

    auto ta1_x_zzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 312);

    auto ta1_x_zzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 313);

    auto ta1_x_zzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 314);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_zzz_xxxx_0,   \
                             ta1_x_zzz_xxxx_1,   \
                             ta1_x_zzz_xxxy_0,   \
                             ta1_x_zzz_xxxy_1,   \
                             ta1_x_zzz_xxxz_0,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxyy_0,   \
                             ta1_x_zzz_xxyy_1,   \
                             ta1_x_zzz_xxyz_0,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxzz_0,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xyyy_0,   \
                             ta1_x_zzz_xyyy_1,   \
                             ta1_x_zzz_xyyz_0,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyzz_0,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xzzz_0,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_yyyy_0,   \
                             ta1_x_zzz_yyyy_1,   \
                             ta1_x_zzz_yyyz_0,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyzz_0,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yzzz_0,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_zzzz_0,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta1_x_zzzz_xxx_0,   \
                             ta1_x_zzzz_xxx_1,   \
                             ta1_x_zzzz_xxxx_0,  \
                             ta1_x_zzzz_xxxx_1,  \
                             ta1_x_zzzz_xxxy_0,  \
                             ta1_x_zzzz_xxxy_1,  \
                             ta1_x_zzzz_xxxz_0,  \
                             ta1_x_zzzz_xxxz_1,  \
                             ta1_x_zzzz_xxy_0,   \
                             ta1_x_zzzz_xxy_1,   \
                             ta1_x_zzzz_xxyy_0,  \
                             ta1_x_zzzz_xxyy_1,  \
                             ta1_x_zzzz_xxyz_0,  \
                             ta1_x_zzzz_xxyz_1,  \
                             ta1_x_zzzz_xxz_0,   \
                             ta1_x_zzzz_xxz_1,   \
                             ta1_x_zzzz_xxzz_0,  \
                             ta1_x_zzzz_xxzz_1,  \
                             ta1_x_zzzz_xyy_0,   \
                             ta1_x_zzzz_xyy_1,   \
                             ta1_x_zzzz_xyyy_0,  \
                             ta1_x_zzzz_xyyy_1,  \
                             ta1_x_zzzz_xyyz_0,  \
                             ta1_x_zzzz_xyyz_1,  \
                             ta1_x_zzzz_xyz_0,   \
                             ta1_x_zzzz_xyz_1,   \
                             ta1_x_zzzz_xyzz_0,  \
                             ta1_x_zzzz_xyzz_1,  \
                             ta1_x_zzzz_xzz_0,   \
                             ta1_x_zzzz_xzz_1,   \
                             ta1_x_zzzz_xzzz_0,  \
                             ta1_x_zzzz_xzzz_1,  \
                             ta1_x_zzzz_yyy_0,   \
                             ta1_x_zzzz_yyy_1,   \
                             ta1_x_zzzz_yyyy_0,  \
                             ta1_x_zzzz_yyyy_1,  \
                             ta1_x_zzzz_yyyz_0,  \
                             ta1_x_zzzz_yyyz_1,  \
                             ta1_x_zzzz_yyz_0,   \
                             ta1_x_zzzz_yyz_1,   \
                             ta1_x_zzzz_yyzz_0,  \
                             ta1_x_zzzz_yyzz_1,  \
                             ta1_x_zzzz_yzz_0,   \
                             ta1_x_zzzz_yzz_1,   \
                             ta1_x_zzzz_yzzz_0,  \
                             ta1_x_zzzz_yzzz_1,  \
                             ta1_x_zzzz_zzz_0,   \
                             ta1_x_zzzz_zzz_1,   \
                             ta1_x_zzzz_zzzz_0,  \
                             ta1_x_zzzz_zzzz_1,  \
                             ta1_x_zzzzz_xxxx_0, \
                             ta1_x_zzzzz_xxxy_0, \
                             ta1_x_zzzzz_xxxz_0, \
                             ta1_x_zzzzz_xxyy_0, \
                             ta1_x_zzzzz_xxyz_0, \
                             ta1_x_zzzzz_xxzz_0, \
                             ta1_x_zzzzz_xyyy_0, \
                             ta1_x_zzzzz_xyyz_0, \
                             ta1_x_zzzzz_xyzz_0, \
                             ta1_x_zzzzz_xzzz_0, \
                             ta1_x_zzzzz_yyyy_0, \
                             ta1_x_zzzzz_yyyz_0, \
                             ta1_x_zzzzz_yyzz_0, \
                             ta1_x_zzzzz_yzzz_0, \
                             ta1_x_zzzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzz_xxxx_0[i] =
            4.0 * ta1_x_zzz_xxxx_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxxx_1[i] * fe_0 + ta1_x_zzzz_xxxx_0[i] * pa_z[i] - ta1_x_zzzz_xxxx_1[i] * pc_z[i];

        ta1_x_zzzzz_xxxy_0[i] =
            4.0 * ta1_x_zzz_xxxy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxxy_1[i] * fe_0 + ta1_x_zzzz_xxxy_0[i] * pa_z[i] - ta1_x_zzzz_xxxy_1[i] * pc_z[i];

        ta1_x_zzzzz_xxxz_0[i] = 4.0 * ta1_x_zzz_xxxz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxxz_1[i] * fe_0 + ta1_x_zzzz_xxx_0[i] * fe_0 -
                                ta1_x_zzzz_xxx_1[i] * fe_0 + ta1_x_zzzz_xxxz_0[i] * pa_z[i] - ta1_x_zzzz_xxxz_1[i] * pc_z[i];

        ta1_x_zzzzz_xxyy_0[i] =
            4.0 * ta1_x_zzz_xxyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxyy_1[i] * fe_0 + ta1_x_zzzz_xxyy_0[i] * pa_z[i] - ta1_x_zzzz_xxyy_1[i] * pc_z[i];

        ta1_x_zzzzz_xxyz_0[i] = 4.0 * ta1_x_zzz_xxyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxyz_1[i] * fe_0 + ta1_x_zzzz_xxy_0[i] * fe_0 -
                                ta1_x_zzzz_xxy_1[i] * fe_0 + ta1_x_zzzz_xxyz_0[i] * pa_z[i] - ta1_x_zzzz_xxyz_1[i] * pc_z[i];

        ta1_x_zzzzz_xxzz_0[i] = 4.0 * ta1_x_zzz_xxzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxzz_1[i] * fe_0 + 2.0 * ta1_x_zzzz_xxz_0[i] * fe_0 -
                                2.0 * ta1_x_zzzz_xxz_1[i] * fe_0 + ta1_x_zzzz_xxzz_0[i] * pa_z[i] - ta1_x_zzzz_xxzz_1[i] * pc_z[i];

        ta1_x_zzzzz_xyyy_0[i] =
            4.0 * ta1_x_zzz_xyyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyyy_1[i] * fe_0 + ta1_x_zzzz_xyyy_0[i] * pa_z[i] - ta1_x_zzzz_xyyy_1[i] * pc_z[i];

        ta1_x_zzzzz_xyyz_0[i] = 4.0 * ta1_x_zzz_xyyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyyz_1[i] * fe_0 + ta1_x_zzzz_xyy_0[i] * fe_0 -
                                ta1_x_zzzz_xyy_1[i] * fe_0 + ta1_x_zzzz_xyyz_0[i] * pa_z[i] - ta1_x_zzzz_xyyz_1[i] * pc_z[i];

        ta1_x_zzzzz_xyzz_0[i] = 4.0 * ta1_x_zzz_xyzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyzz_1[i] * fe_0 + 2.0 * ta1_x_zzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_x_zzzz_xyz_1[i] * fe_0 + ta1_x_zzzz_xyzz_0[i] * pa_z[i] - ta1_x_zzzz_xyzz_1[i] * pc_z[i];

        ta1_x_zzzzz_xzzz_0[i] = 4.0 * ta1_x_zzz_xzzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xzzz_1[i] * fe_0 + 3.0 * ta1_x_zzzz_xzz_0[i] * fe_0 -
                                3.0 * ta1_x_zzzz_xzz_1[i] * fe_0 + ta1_x_zzzz_xzzz_0[i] * pa_z[i] - ta1_x_zzzz_xzzz_1[i] * pc_z[i];

        ta1_x_zzzzz_yyyy_0[i] =
            4.0 * ta1_x_zzz_yyyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyyy_1[i] * fe_0 + ta1_x_zzzz_yyyy_0[i] * pa_z[i] - ta1_x_zzzz_yyyy_1[i] * pc_z[i];

        ta1_x_zzzzz_yyyz_0[i] = 4.0 * ta1_x_zzz_yyyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyyz_1[i] * fe_0 + ta1_x_zzzz_yyy_0[i] * fe_0 -
                                ta1_x_zzzz_yyy_1[i] * fe_0 + ta1_x_zzzz_yyyz_0[i] * pa_z[i] - ta1_x_zzzz_yyyz_1[i] * pc_z[i];

        ta1_x_zzzzz_yyzz_0[i] = 4.0 * ta1_x_zzz_yyzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_zzzz_yyz_0[i] * fe_0 -
                                2.0 * ta1_x_zzzz_yyz_1[i] * fe_0 + ta1_x_zzzz_yyzz_0[i] * pa_z[i] - ta1_x_zzzz_yyzz_1[i] * pc_z[i];

        ta1_x_zzzzz_yzzz_0[i] = 4.0 * ta1_x_zzz_yzzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yzzz_1[i] * fe_0 + 3.0 * ta1_x_zzzz_yzz_0[i] * fe_0 -
                                3.0 * ta1_x_zzzz_yzz_1[i] * fe_0 + ta1_x_zzzz_yzzz_0[i] * pa_z[i] - ta1_x_zzzz_yzzz_1[i] * pc_z[i];

        ta1_x_zzzzz_zzzz_0[i] = 4.0 * ta1_x_zzz_zzzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_zzzz_1[i] * fe_0 + 4.0 * ta1_x_zzzz_zzz_0[i] * fe_0 -
                                4.0 * ta1_x_zzzz_zzz_1[i] * fe_0 + ta1_x_zzzz_zzzz_0[i] * pa_z[i] - ta1_x_zzzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 315-330 components of targeted buffer : HG

    auto ta1_y_xxxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 315);

    auto ta1_y_xxxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 316);

    auto ta1_y_xxxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 317);

    auto ta1_y_xxxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 318);

    auto ta1_y_xxxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 319);

    auto ta1_y_xxxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 320);

    auto ta1_y_xxxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 321);

    auto ta1_y_xxxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 322);

    auto ta1_y_xxxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 323);

    auto ta1_y_xxxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 324);

    auto ta1_y_xxxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 325);

    auto ta1_y_xxxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 326);

    auto ta1_y_xxxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 327);

    auto ta1_y_xxxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 328);

    auto ta1_y_xxxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 329);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xxx_xxxx_0,   \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxy_0,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxz_0,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxyy_0,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyz_0,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxzz_0,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xyyy_0,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyz_0,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyzz_0,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xzzz_0,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_yyyy_0,   \
                             ta1_y_xxx_yyyy_1,   \
                             ta1_y_xxx_yyyz_0,   \
                             ta1_y_xxx_yyyz_1,   \
                             ta1_y_xxx_yyzz_0,   \
                             ta1_y_xxx_yyzz_1,   \
                             ta1_y_xxx_yzzz_0,   \
                             ta1_y_xxx_yzzz_1,   \
                             ta1_y_xxx_zzzz_0,   \
                             ta1_y_xxx_zzzz_1,   \
                             ta1_y_xxxx_xxx_0,   \
                             ta1_y_xxxx_xxx_1,   \
                             ta1_y_xxxx_xxxx_0,  \
                             ta1_y_xxxx_xxxx_1,  \
                             ta1_y_xxxx_xxxy_0,  \
                             ta1_y_xxxx_xxxy_1,  \
                             ta1_y_xxxx_xxxz_0,  \
                             ta1_y_xxxx_xxxz_1,  \
                             ta1_y_xxxx_xxy_0,   \
                             ta1_y_xxxx_xxy_1,   \
                             ta1_y_xxxx_xxyy_0,  \
                             ta1_y_xxxx_xxyy_1,  \
                             ta1_y_xxxx_xxyz_0,  \
                             ta1_y_xxxx_xxyz_1,  \
                             ta1_y_xxxx_xxz_0,   \
                             ta1_y_xxxx_xxz_1,   \
                             ta1_y_xxxx_xxzz_0,  \
                             ta1_y_xxxx_xxzz_1,  \
                             ta1_y_xxxx_xyy_0,   \
                             ta1_y_xxxx_xyy_1,   \
                             ta1_y_xxxx_xyyy_0,  \
                             ta1_y_xxxx_xyyy_1,  \
                             ta1_y_xxxx_xyyz_0,  \
                             ta1_y_xxxx_xyyz_1,  \
                             ta1_y_xxxx_xyz_0,   \
                             ta1_y_xxxx_xyz_1,   \
                             ta1_y_xxxx_xyzz_0,  \
                             ta1_y_xxxx_xyzz_1,  \
                             ta1_y_xxxx_xzz_0,   \
                             ta1_y_xxxx_xzz_1,   \
                             ta1_y_xxxx_xzzz_0,  \
                             ta1_y_xxxx_xzzz_1,  \
                             ta1_y_xxxx_yyy_0,   \
                             ta1_y_xxxx_yyy_1,   \
                             ta1_y_xxxx_yyyy_0,  \
                             ta1_y_xxxx_yyyy_1,  \
                             ta1_y_xxxx_yyyz_0,  \
                             ta1_y_xxxx_yyyz_1,  \
                             ta1_y_xxxx_yyz_0,   \
                             ta1_y_xxxx_yyz_1,   \
                             ta1_y_xxxx_yyzz_0,  \
                             ta1_y_xxxx_yyzz_1,  \
                             ta1_y_xxxx_yzz_0,   \
                             ta1_y_xxxx_yzz_1,   \
                             ta1_y_xxxx_yzzz_0,  \
                             ta1_y_xxxx_yzzz_1,  \
                             ta1_y_xxxx_zzz_0,   \
                             ta1_y_xxxx_zzz_1,   \
                             ta1_y_xxxx_zzzz_0,  \
                             ta1_y_xxxx_zzzz_1,  \
                             ta1_y_xxxxx_xxxx_0, \
                             ta1_y_xxxxx_xxxy_0, \
                             ta1_y_xxxxx_xxxz_0, \
                             ta1_y_xxxxx_xxyy_0, \
                             ta1_y_xxxxx_xxyz_0, \
                             ta1_y_xxxxx_xxzz_0, \
                             ta1_y_xxxxx_xyyy_0, \
                             ta1_y_xxxxx_xyyz_0, \
                             ta1_y_xxxxx_xyzz_0, \
                             ta1_y_xxxxx_xzzz_0, \
                             ta1_y_xxxxx_yyyy_0, \
                             ta1_y_xxxxx_yyyz_0, \
                             ta1_y_xxxxx_yyzz_0, \
                             ta1_y_xxxxx_yzzz_0, \
                             ta1_y_xxxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxx_xxxx_0[i] = 4.0 * ta1_y_xxx_xxxx_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxxx_1[i] * fe_0 + 4.0 * ta1_y_xxxx_xxx_0[i] * fe_0 -
                                4.0 * ta1_y_xxxx_xxx_1[i] * fe_0 + ta1_y_xxxx_xxxx_0[i] * pa_x[i] - ta1_y_xxxx_xxxx_1[i] * pc_x[i];

        ta1_y_xxxxx_xxxy_0[i] = 4.0 * ta1_y_xxx_xxxy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xxxx_xxy_0[i] * fe_0 -
                                3.0 * ta1_y_xxxx_xxy_1[i] * fe_0 + ta1_y_xxxx_xxxy_0[i] * pa_x[i] - ta1_y_xxxx_xxxy_1[i] * pc_x[i];

        ta1_y_xxxxx_xxxz_0[i] = 4.0 * ta1_y_xxx_xxxz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xxxx_xxz_0[i] * fe_0 -
                                3.0 * ta1_y_xxxx_xxz_1[i] * fe_0 + ta1_y_xxxx_xxxz_0[i] * pa_x[i] - ta1_y_xxxx_xxxz_1[i] * pc_x[i];

        ta1_y_xxxxx_xxyy_0[i] = 4.0 * ta1_y_xxx_xxyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xxxx_xyy_0[i] * fe_0 -
                                2.0 * ta1_y_xxxx_xyy_1[i] * fe_0 + ta1_y_xxxx_xxyy_0[i] * pa_x[i] - ta1_y_xxxx_xxyy_1[i] * pc_x[i];

        ta1_y_xxxxx_xxyz_0[i] = 4.0 * ta1_y_xxx_xxyz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxxx_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_xxxx_xyz_1[i] * fe_0 + ta1_y_xxxx_xxyz_0[i] * pa_x[i] - ta1_y_xxxx_xxyz_1[i] * pc_x[i];

        ta1_y_xxxxx_xxzz_0[i] = 4.0 * ta1_y_xxx_xxzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxxx_xzz_0[i] * fe_0 -
                                2.0 * ta1_y_xxxx_xzz_1[i] * fe_0 + ta1_y_xxxx_xxzz_0[i] * pa_x[i] - ta1_y_xxxx_xxzz_1[i] * pc_x[i];

        ta1_y_xxxxx_xyyy_0[i] = 4.0 * ta1_y_xxx_xyyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyyy_1[i] * fe_0 + ta1_y_xxxx_yyy_0[i] * fe_0 -
                                ta1_y_xxxx_yyy_1[i] * fe_0 + ta1_y_xxxx_xyyy_0[i] * pa_x[i] - ta1_y_xxxx_xyyy_1[i] * pc_x[i];

        ta1_y_xxxxx_xyyz_0[i] = 4.0 * ta1_y_xxx_xyyz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyyz_1[i] * fe_0 + ta1_y_xxxx_yyz_0[i] * fe_0 -
                                ta1_y_xxxx_yyz_1[i] * fe_0 + ta1_y_xxxx_xyyz_0[i] * pa_x[i] - ta1_y_xxxx_xyyz_1[i] * pc_x[i];

        ta1_y_xxxxx_xyzz_0[i] = 4.0 * ta1_y_xxx_xyzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyzz_1[i] * fe_0 + ta1_y_xxxx_yzz_0[i] * fe_0 -
                                ta1_y_xxxx_yzz_1[i] * fe_0 + ta1_y_xxxx_xyzz_0[i] * pa_x[i] - ta1_y_xxxx_xyzz_1[i] * pc_x[i];

        ta1_y_xxxxx_xzzz_0[i] = 4.0 * ta1_y_xxx_xzzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xzzz_1[i] * fe_0 + ta1_y_xxxx_zzz_0[i] * fe_0 -
                                ta1_y_xxxx_zzz_1[i] * fe_0 + ta1_y_xxxx_xzzz_0[i] * pa_x[i] - ta1_y_xxxx_xzzz_1[i] * pc_x[i];

        ta1_y_xxxxx_yyyy_0[i] =
            4.0 * ta1_y_xxx_yyyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_yyyy_1[i] * fe_0 + ta1_y_xxxx_yyyy_0[i] * pa_x[i] - ta1_y_xxxx_yyyy_1[i] * pc_x[i];

        ta1_y_xxxxx_yyyz_0[i] =
            4.0 * ta1_y_xxx_yyyz_0[i] * fe_0 - 4.0 * ta1_y_xxx_yyyz_1[i] * fe_0 + ta1_y_xxxx_yyyz_0[i] * pa_x[i] - ta1_y_xxxx_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxx_yyzz_0[i] =
            4.0 * ta1_y_xxx_yyzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_yyzz_1[i] * fe_0 + ta1_y_xxxx_yyzz_0[i] * pa_x[i] - ta1_y_xxxx_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxx_yzzz_0[i] =
            4.0 * ta1_y_xxx_yzzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_yzzz_1[i] * fe_0 + ta1_y_xxxx_yzzz_0[i] * pa_x[i] - ta1_y_xxxx_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxx_zzzz_0[i] =
            4.0 * ta1_y_xxx_zzzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_zzzz_1[i] * fe_0 + ta1_y_xxxx_zzzz_0[i] * pa_x[i] - ta1_y_xxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 330-345 components of targeted buffer : HG

    auto ta1_y_xxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 330);

    auto ta1_y_xxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 331);

    auto ta1_y_xxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 332);

    auto ta1_y_xxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 333);

    auto ta1_y_xxxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 334);

    auto ta1_y_xxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 335);

    auto ta1_y_xxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 336);

    auto ta1_y_xxxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 337);

    auto ta1_y_xxxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 338);

    auto ta1_y_xxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 339);

    auto ta1_y_xxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 340);

    auto ta1_y_xxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 341);

    auto ta1_y_xxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 342);

    auto ta1_y_xxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 343);

    auto ta1_y_xxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 344);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxxx_xxx_0,   \
                             ta1_y_xxxx_xxx_1,   \
                             ta1_y_xxxx_xxxx_0,  \
                             ta1_y_xxxx_xxxx_1,  \
                             ta1_y_xxxx_xxxy_0,  \
                             ta1_y_xxxx_xxxy_1,  \
                             ta1_y_xxxx_xxxz_0,  \
                             ta1_y_xxxx_xxxz_1,  \
                             ta1_y_xxxx_xxy_0,   \
                             ta1_y_xxxx_xxy_1,   \
                             ta1_y_xxxx_xxyy_0,  \
                             ta1_y_xxxx_xxyy_1,  \
                             ta1_y_xxxx_xxyz_0,  \
                             ta1_y_xxxx_xxyz_1,  \
                             ta1_y_xxxx_xxz_0,   \
                             ta1_y_xxxx_xxz_1,   \
                             ta1_y_xxxx_xxzz_0,  \
                             ta1_y_xxxx_xxzz_1,  \
                             ta1_y_xxxx_xyy_0,   \
                             ta1_y_xxxx_xyy_1,   \
                             ta1_y_xxxx_xyyy_0,  \
                             ta1_y_xxxx_xyyy_1,  \
                             ta1_y_xxxx_xyyz_0,  \
                             ta1_y_xxxx_xyyz_1,  \
                             ta1_y_xxxx_xyz_0,   \
                             ta1_y_xxxx_xyz_1,   \
                             ta1_y_xxxx_xyzz_0,  \
                             ta1_y_xxxx_xyzz_1,  \
                             ta1_y_xxxx_xzz_0,   \
                             ta1_y_xxxx_xzz_1,   \
                             ta1_y_xxxx_xzzz_0,  \
                             ta1_y_xxxx_xzzz_1,  \
                             ta1_y_xxxx_zzzz_0,  \
                             ta1_y_xxxx_zzzz_1,  \
                             ta1_y_xxxxy_xxxx_0, \
                             ta1_y_xxxxy_xxxy_0, \
                             ta1_y_xxxxy_xxxz_0, \
                             ta1_y_xxxxy_xxyy_0, \
                             ta1_y_xxxxy_xxyz_0, \
                             ta1_y_xxxxy_xxzz_0, \
                             ta1_y_xxxxy_xyyy_0, \
                             ta1_y_xxxxy_xyyz_0, \
                             ta1_y_xxxxy_xyzz_0, \
                             ta1_y_xxxxy_xzzz_0, \
                             ta1_y_xxxxy_yyyy_0, \
                             ta1_y_xxxxy_yyyz_0, \
                             ta1_y_xxxxy_yyzz_0, \
                             ta1_y_xxxxy_yzzz_0, \
                             ta1_y_xxxxy_zzzz_0, \
                             ta1_y_xxxy_yyyy_0,  \
                             ta1_y_xxxy_yyyy_1,  \
                             ta1_y_xxxy_yyyz_0,  \
                             ta1_y_xxxy_yyyz_1,  \
                             ta1_y_xxxy_yyzz_0,  \
                             ta1_y_xxxy_yyzz_1,  \
                             ta1_y_xxxy_yzzz_0,  \
                             ta1_y_xxxy_yzzz_1,  \
                             ta1_y_xxy_yyyy_0,   \
                             ta1_y_xxy_yyyy_1,   \
                             ta1_y_xxy_yyyz_0,   \
                             ta1_y_xxy_yyyz_1,   \
                             ta1_y_xxy_yyzz_0,   \
                             ta1_y_xxy_yyzz_1,   \
                             ta1_y_xxy_yzzz_0,   \
                             ta1_y_xxy_yzzz_1,   \
                             ta_xxxx_xxxx_1,     \
                             ta_xxxx_xxxy_1,     \
                             ta_xxxx_xxxz_1,     \
                             ta_xxxx_xxyy_1,     \
                             ta_xxxx_xxyz_1,     \
                             ta_xxxx_xxzz_1,     \
                             ta_xxxx_xyyy_1,     \
                             ta_xxxx_xyyz_1,     \
                             ta_xxxx_xyzz_1,     \
                             ta_xxxx_xzzz_1,     \
                             ta_xxxx_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxy_xxxx_0[i] = ta_xxxx_xxxx_1[i] + ta1_y_xxxx_xxxx_0[i] * pa_y[i] - ta1_y_xxxx_xxxx_1[i] * pc_y[i];

        ta1_y_xxxxy_xxxy_0[i] = ta1_y_xxxx_xxx_0[i] * fe_0 - ta1_y_xxxx_xxx_1[i] * fe_0 + ta_xxxx_xxxy_1[i] + ta1_y_xxxx_xxxy_0[i] * pa_y[i] -
                                ta1_y_xxxx_xxxy_1[i] * pc_y[i];

        ta1_y_xxxxy_xxxz_0[i] = ta_xxxx_xxxz_1[i] + ta1_y_xxxx_xxxz_0[i] * pa_y[i] - ta1_y_xxxx_xxxz_1[i] * pc_y[i];

        ta1_y_xxxxy_xxyy_0[i] = 2.0 * ta1_y_xxxx_xxy_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xxy_1[i] * fe_0 + ta_xxxx_xxyy_1[i] +
                                ta1_y_xxxx_xxyy_0[i] * pa_y[i] - ta1_y_xxxx_xxyy_1[i] * pc_y[i];

        ta1_y_xxxxy_xxyz_0[i] = ta1_y_xxxx_xxz_0[i] * fe_0 - ta1_y_xxxx_xxz_1[i] * fe_0 + ta_xxxx_xxyz_1[i] + ta1_y_xxxx_xxyz_0[i] * pa_y[i] -
                                ta1_y_xxxx_xxyz_1[i] * pc_y[i];

        ta1_y_xxxxy_xxzz_0[i] = ta_xxxx_xxzz_1[i] + ta1_y_xxxx_xxzz_0[i] * pa_y[i] - ta1_y_xxxx_xxzz_1[i] * pc_y[i];

        ta1_y_xxxxy_xyyy_0[i] = 3.0 * ta1_y_xxxx_xyy_0[i] * fe_0 - 3.0 * ta1_y_xxxx_xyy_1[i] * fe_0 + ta_xxxx_xyyy_1[i] +
                                ta1_y_xxxx_xyyy_0[i] * pa_y[i] - ta1_y_xxxx_xyyy_1[i] * pc_y[i];

        ta1_y_xxxxy_xyyz_0[i] = 2.0 * ta1_y_xxxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xyz_1[i] * fe_0 + ta_xxxx_xyyz_1[i] +
                                ta1_y_xxxx_xyyz_0[i] * pa_y[i] - ta1_y_xxxx_xyyz_1[i] * pc_y[i];

        ta1_y_xxxxy_xyzz_0[i] = ta1_y_xxxx_xzz_0[i] * fe_0 - ta1_y_xxxx_xzz_1[i] * fe_0 + ta_xxxx_xyzz_1[i] + ta1_y_xxxx_xyzz_0[i] * pa_y[i] -
                                ta1_y_xxxx_xyzz_1[i] * pc_y[i];

        ta1_y_xxxxy_xzzz_0[i] = ta_xxxx_xzzz_1[i] + ta1_y_xxxx_xzzz_0[i] * pa_y[i] - ta1_y_xxxx_xzzz_1[i] * pc_y[i];

        ta1_y_xxxxy_yyyy_0[i] =
            3.0 * ta1_y_xxy_yyyy_0[i] * fe_0 - 3.0 * ta1_y_xxy_yyyy_1[i] * fe_0 + ta1_y_xxxy_yyyy_0[i] * pa_x[i] - ta1_y_xxxy_yyyy_1[i] * pc_x[i];

        ta1_y_xxxxy_yyyz_0[i] =
            3.0 * ta1_y_xxy_yyyz_0[i] * fe_0 - 3.0 * ta1_y_xxy_yyyz_1[i] * fe_0 + ta1_y_xxxy_yyyz_0[i] * pa_x[i] - ta1_y_xxxy_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxy_yyzz_0[i] =
            3.0 * ta1_y_xxy_yyzz_0[i] * fe_0 - 3.0 * ta1_y_xxy_yyzz_1[i] * fe_0 + ta1_y_xxxy_yyzz_0[i] * pa_x[i] - ta1_y_xxxy_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxy_yzzz_0[i] =
            3.0 * ta1_y_xxy_yzzz_0[i] * fe_0 - 3.0 * ta1_y_xxy_yzzz_1[i] * fe_0 + ta1_y_xxxy_yzzz_0[i] * pa_x[i] - ta1_y_xxxy_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxy_zzzz_0[i] = ta_xxxx_zzzz_1[i] + ta1_y_xxxx_zzzz_0[i] * pa_y[i] - ta1_y_xxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 345-360 components of targeted buffer : HG

    auto ta1_y_xxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 345);

    auto ta1_y_xxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 346);

    auto ta1_y_xxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 347);

    auto ta1_y_xxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 348);

    auto ta1_y_xxxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 349);

    auto ta1_y_xxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 350);

    auto ta1_y_xxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 351);

    auto ta1_y_xxxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 352);

    auto ta1_y_xxxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 353);

    auto ta1_y_xxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 354);

    auto ta1_y_xxxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 355);

    auto ta1_y_xxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 356);

    auto ta1_y_xxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 357);

    auto ta1_y_xxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 358);

    auto ta1_y_xxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 359);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxxx_xxx_0,   \
                             ta1_y_xxxx_xxx_1,   \
                             ta1_y_xxxx_xxxx_0,  \
                             ta1_y_xxxx_xxxx_1,  \
                             ta1_y_xxxx_xxxy_0,  \
                             ta1_y_xxxx_xxxy_1,  \
                             ta1_y_xxxx_xxxz_0,  \
                             ta1_y_xxxx_xxxz_1,  \
                             ta1_y_xxxx_xxy_0,   \
                             ta1_y_xxxx_xxy_1,   \
                             ta1_y_xxxx_xxyy_0,  \
                             ta1_y_xxxx_xxyy_1,  \
                             ta1_y_xxxx_xxyz_0,  \
                             ta1_y_xxxx_xxyz_1,  \
                             ta1_y_xxxx_xxz_0,   \
                             ta1_y_xxxx_xxz_1,   \
                             ta1_y_xxxx_xxzz_0,  \
                             ta1_y_xxxx_xxzz_1,  \
                             ta1_y_xxxx_xyy_0,   \
                             ta1_y_xxxx_xyy_1,   \
                             ta1_y_xxxx_xyyy_0,  \
                             ta1_y_xxxx_xyyy_1,  \
                             ta1_y_xxxx_xyyz_0,  \
                             ta1_y_xxxx_xyyz_1,  \
                             ta1_y_xxxx_xyz_0,   \
                             ta1_y_xxxx_xyz_1,   \
                             ta1_y_xxxx_xyzz_0,  \
                             ta1_y_xxxx_xyzz_1,  \
                             ta1_y_xxxx_xzz_0,   \
                             ta1_y_xxxx_xzz_1,   \
                             ta1_y_xxxx_xzzz_0,  \
                             ta1_y_xxxx_xzzz_1,  \
                             ta1_y_xxxx_yyyy_0,  \
                             ta1_y_xxxx_yyyy_1,  \
                             ta1_y_xxxxz_xxxx_0, \
                             ta1_y_xxxxz_xxxy_0, \
                             ta1_y_xxxxz_xxxz_0, \
                             ta1_y_xxxxz_xxyy_0, \
                             ta1_y_xxxxz_xxyz_0, \
                             ta1_y_xxxxz_xxzz_0, \
                             ta1_y_xxxxz_xyyy_0, \
                             ta1_y_xxxxz_xyyz_0, \
                             ta1_y_xxxxz_xyzz_0, \
                             ta1_y_xxxxz_xzzz_0, \
                             ta1_y_xxxxz_yyyy_0, \
                             ta1_y_xxxxz_yyyz_0, \
                             ta1_y_xxxxz_yyzz_0, \
                             ta1_y_xxxxz_yzzz_0, \
                             ta1_y_xxxxz_zzzz_0, \
                             ta1_y_xxxz_yyyz_0,  \
                             ta1_y_xxxz_yyyz_1,  \
                             ta1_y_xxxz_yyzz_0,  \
                             ta1_y_xxxz_yyzz_1,  \
                             ta1_y_xxxz_yzzz_0,  \
                             ta1_y_xxxz_yzzz_1,  \
                             ta1_y_xxxz_zzzz_0,  \
                             ta1_y_xxxz_zzzz_1,  \
                             ta1_y_xxz_yyyz_0,   \
                             ta1_y_xxz_yyyz_1,   \
                             ta1_y_xxz_yyzz_0,   \
                             ta1_y_xxz_yyzz_1,   \
                             ta1_y_xxz_yzzz_0,   \
                             ta1_y_xxz_yzzz_1,   \
                             ta1_y_xxz_zzzz_0,   \
                             ta1_y_xxz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxz_xxxx_0[i] = ta1_y_xxxx_xxxx_0[i] * pa_z[i] - ta1_y_xxxx_xxxx_1[i] * pc_z[i];

        ta1_y_xxxxz_xxxy_0[i] = ta1_y_xxxx_xxxy_0[i] * pa_z[i] - ta1_y_xxxx_xxxy_1[i] * pc_z[i];

        ta1_y_xxxxz_xxxz_0[i] =
            ta1_y_xxxx_xxx_0[i] * fe_0 - ta1_y_xxxx_xxx_1[i] * fe_0 + ta1_y_xxxx_xxxz_0[i] * pa_z[i] - ta1_y_xxxx_xxxz_1[i] * pc_z[i];

        ta1_y_xxxxz_xxyy_0[i] = ta1_y_xxxx_xxyy_0[i] * pa_z[i] - ta1_y_xxxx_xxyy_1[i] * pc_z[i];

        ta1_y_xxxxz_xxyz_0[i] =
            ta1_y_xxxx_xxy_0[i] * fe_0 - ta1_y_xxxx_xxy_1[i] * fe_0 + ta1_y_xxxx_xxyz_0[i] * pa_z[i] - ta1_y_xxxx_xxyz_1[i] * pc_z[i];

        ta1_y_xxxxz_xxzz_0[i] =
            2.0 * ta1_y_xxxx_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xxz_1[i] * fe_0 + ta1_y_xxxx_xxzz_0[i] * pa_z[i] - ta1_y_xxxx_xxzz_1[i] * pc_z[i];

        ta1_y_xxxxz_xyyy_0[i] = ta1_y_xxxx_xyyy_0[i] * pa_z[i] - ta1_y_xxxx_xyyy_1[i] * pc_z[i];

        ta1_y_xxxxz_xyyz_0[i] =
            ta1_y_xxxx_xyy_0[i] * fe_0 - ta1_y_xxxx_xyy_1[i] * fe_0 + ta1_y_xxxx_xyyz_0[i] * pa_z[i] - ta1_y_xxxx_xyyz_1[i] * pc_z[i];

        ta1_y_xxxxz_xyzz_0[i] =
            2.0 * ta1_y_xxxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xyz_1[i] * fe_0 + ta1_y_xxxx_xyzz_0[i] * pa_z[i] - ta1_y_xxxx_xyzz_1[i] * pc_z[i];

        ta1_y_xxxxz_xzzz_0[i] =
            3.0 * ta1_y_xxxx_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxxx_xzz_1[i] * fe_0 + ta1_y_xxxx_xzzz_0[i] * pa_z[i] - ta1_y_xxxx_xzzz_1[i] * pc_z[i];

        ta1_y_xxxxz_yyyy_0[i] = ta1_y_xxxx_yyyy_0[i] * pa_z[i] - ta1_y_xxxx_yyyy_1[i] * pc_z[i];

        ta1_y_xxxxz_yyyz_0[i] =
            3.0 * ta1_y_xxz_yyyz_0[i] * fe_0 - 3.0 * ta1_y_xxz_yyyz_1[i] * fe_0 + ta1_y_xxxz_yyyz_0[i] * pa_x[i] - ta1_y_xxxz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxz_yyzz_0[i] =
            3.0 * ta1_y_xxz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_xxz_yyzz_1[i] * fe_0 + ta1_y_xxxz_yyzz_0[i] * pa_x[i] - ta1_y_xxxz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxz_yzzz_0[i] =
            3.0 * ta1_y_xxz_yzzz_0[i] * fe_0 - 3.0 * ta1_y_xxz_yzzz_1[i] * fe_0 + ta1_y_xxxz_yzzz_0[i] * pa_x[i] - ta1_y_xxxz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxz_zzzz_0[i] =
            3.0 * ta1_y_xxz_zzzz_0[i] * fe_0 - 3.0 * ta1_y_xxz_zzzz_1[i] * fe_0 + ta1_y_xxxz_zzzz_0[i] * pa_x[i] - ta1_y_xxxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 360-375 components of targeted buffer : HG

    auto ta1_y_xxxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 360);

    auto ta1_y_xxxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 361);

    auto ta1_y_xxxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 362);

    auto ta1_y_xxxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 363);

    auto ta1_y_xxxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 364);

    auto ta1_y_xxxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 365);

    auto ta1_y_xxxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 366);

    auto ta1_y_xxxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 367);

    auto ta1_y_xxxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 368);

    auto ta1_y_xxxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 369);

    auto ta1_y_xxxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 370);

    auto ta1_y_xxxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 371);

    auto ta1_y_xxxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 372);

    auto ta1_y_xxxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 373);

    auto ta1_y_xxxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 374);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxx_xxxx_0,   \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxz_0,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxzz_0,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xzzz_0,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxxy_xxxx_0,  \
                             ta1_y_xxxy_xxxx_1,  \
                             ta1_y_xxxy_xxxz_0,  \
                             ta1_y_xxxy_xxxz_1,  \
                             ta1_y_xxxy_xxzz_0,  \
                             ta1_y_xxxy_xxzz_1,  \
                             ta1_y_xxxy_xzzz_0,  \
                             ta1_y_xxxy_xzzz_1,  \
                             ta1_y_xxxyy_xxxx_0, \
                             ta1_y_xxxyy_xxxy_0, \
                             ta1_y_xxxyy_xxxz_0, \
                             ta1_y_xxxyy_xxyy_0, \
                             ta1_y_xxxyy_xxyz_0, \
                             ta1_y_xxxyy_xxzz_0, \
                             ta1_y_xxxyy_xyyy_0, \
                             ta1_y_xxxyy_xyyz_0, \
                             ta1_y_xxxyy_xyzz_0, \
                             ta1_y_xxxyy_xzzz_0, \
                             ta1_y_xxxyy_yyyy_0, \
                             ta1_y_xxxyy_yyyz_0, \
                             ta1_y_xxxyy_yyzz_0, \
                             ta1_y_xxxyy_yzzz_0, \
                             ta1_y_xxxyy_zzzz_0, \
                             ta1_y_xxyy_xxxy_0,  \
                             ta1_y_xxyy_xxxy_1,  \
                             ta1_y_xxyy_xxy_0,   \
                             ta1_y_xxyy_xxy_1,   \
                             ta1_y_xxyy_xxyy_0,  \
                             ta1_y_xxyy_xxyy_1,  \
                             ta1_y_xxyy_xxyz_0,  \
                             ta1_y_xxyy_xxyz_1,  \
                             ta1_y_xxyy_xyy_0,   \
                             ta1_y_xxyy_xyy_1,   \
                             ta1_y_xxyy_xyyy_0,  \
                             ta1_y_xxyy_xyyy_1,  \
                             ta1_y_xxyy_xyyz_0,  \
                             ta1_y_xxyy_xyyz_1,  \
                             ta1_y_xxyy_xyz_0,   \
                             ta1_y_xxyy_xyz_1,   \
                             ta1_y_xxyy_xyzz_0,  \
                             ta1_y_xxyy_xyzz_1,  \
                             ta1_y_xxyy_yyy_0,   \
                             ta1_y_xxyy_yyy_1,   \
                             ta1_y_xxyy_yyyy_0,  \
                             ta1_y_xxyy_yyyy_1,  \
                             ta1_y_xxyy_yyyz_0,  \
                             ta1_y_xxyy_yyyz_1,  \
                             ta1_y_xxyy_yyz_0,   \
                             ta1_y_xxyy_yyz_1,   \
                             ta1_y_xxyy_yyzz_0,  \
                             ta1_y_xxyy_yyzz_1,  \
                             ta1_y_xxyy_yzz_0,   \
                             ta1_y_xxyy_yzz_1,   \
                             ta1_y_xxyy_yzzz_0,  \
                             ta1_y_xxyy_yzzz_1,  \
                             ta1_y_xxyy_zzzz_0,  \
                             ta1_y_xxyy_zzzz_1,  \
                             ta1_y_xyy_xxxy_0,   \
                             ta1_y_xyy_xxxy_1,   \
                             ta1_y_xyy_xxyy_0,   \
                             ta1_y_xyy_xxyy_1,   \
                             ta1_y_xyy_xxyz_0,   \
                             ta1_y_xyy_xxyz_1,   \
                             ta1_y_xyy_xyyy_0,   \
                             ta1_y_xyy_xyyy_1,   \
                             ta1_y_xyy_xyyz_0,   \
                             ta1_y_xyy_xyyz_1,   \
                             ta1_y_xyy_xyzz_0,   \
                             ta1_y_xyy_xyzz_1,   \
                             ta1_y_xyy_yyyy_0,   \
                             ta1_y_xyy_yyyy_1,   \
                             ta1_y_xyy_yyyz_0,   \
                             ta1_y_xyy_yyyz_1,   \
                             ta1_y_xyy_yyzz_0,   \
                             ta1_y_xyy_yyzz_1,   \
                             ta1_y_xyy_yzzz_0,   \
                             ta1_y_xyy_yzzz_1,   \
                             ta1_y_xyy_zzzz_0,   \
                             ta1_y_xyy_zzzz_1,   \
                             ta_xxxy_xxxx_1,     \
                             ta_xxxy_xxxz_1,     \
                             ta_xxxy_xxzz_1,     \
                             ta_xxxy_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyy_xxxx_0[i] = ta1_y_xxx_xxxx_0[i] * fe_0 - ta1_y_xxx_xxxx_1[i] * fe_0 + ta_xxxy_xxxx_1[i] + ta1_y_xxxy_xxxx_0[i] * pa_y[i] -
                                ta1_y_xxxy_xxxx_1[i] * pc_y[i];

        ta1_y_xxxyy_xxxy_0[i] = 2.0 * ta1_y_xyy_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xxyy_xxy_0[i] * fe_0 -
                                3.0 * ta1_y_xxyy_xxy_1[i] * fe_0 + ta1_y_xxyy_xxxy_0[i] * pa_x[i] - ta1_y_xxyy_xxxy_1[i] * pc_x[i];

        ta1_y_xxxyy_xxxz_0[i] = ta1_y_xxx_xxxz_0[i] * fe_0 - ta1_y_xxx_xxxz_1[i] * fe_0 + ta_xxxy_xxxz_1[i] + ta1_y_xxxy_xxxz_0[i] * pa_y[i] -
                                ta1_y_xxxy_xxxz_1[i] * pc_y[i];

        ta1_y_xxxyy_xxyy_0[i] = 2.0 * ta1_y_xyy_xxyy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xxyy_xyy_0[i] * fe_0 -
                                2.0 * ta1_y_xxyy_xyy_1[i] * fe_0 + ta1_y_xxyy_xxyy_0[i] * pa_x[i] - ta1_y_xxyy_xxyy_1[i] * pc_x[i];

        ta1_y_xxxyy_xxyz_0[i] = 2.0 * ta1_y_xyy_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xyy_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_xxyy_xyz_1[i] * fe_0 + ta1_y_xxyy_xxyz_0[i] * pa_x[i] - ta1_y_xxyy_xxyz_1[i] * pc_x[i];

        ta1_y_xxxyy_xxzz_0[i] = ta1_y_xxx_xxzz_0[i] * fe_0 - ta1_y_xxx_xxzz_1[i] * fe_0 + ta_xxxy_xxzz_1[i] + ta1_y_xxxy_xxzz_0[i] * pa_y[i] -
                                ta1_y_xxxy_xxzz_1[i] * pc_y[i];

        ta1_y_xxxyy_xyyy_0[i] = 2.0 * ta1_y_xyy_xyyy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xyyy_1[i] * fe_0 + ta1_y_xxyy_yyy_0[i] * fe_0 -
                                ta1_y_xxyy_yyy_1[i] * fe_0 + ta1_y_xxyy_xyyy_0[i] * pa_x[i] - ta1_y_xxyy_xyyy_1[i] * pc_x[i];

        ta1_y_xxxyy_xyyz_0[i] = 2.0 * ta1_y_xyy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xyy_xyyz_1[i] * fe_0 + ta1_y_xxyy_yyz_0[i] * fe_0 -
                                ta1_y_xxyy_yyz_1[i] * fe_0 + ta1_y_xxyy_xyyz_0[i] * pa_x[i] - ta1_y_xxyy_xyyz_1[i] * pc_x[i];

        ta1_y_xxxyy_xyzz_0[i] = 2.0 * ta1_y_xyy_xyzz_0[i] * fe_0 - 2.0 * ta1_y_xyy_xyzz_1[i] * fe_0 + ta1_y_xxyy_yzz_0[i] * fe_0 -
                                ta1_y_xxyy_yzz_1[i] * fe_0 + ta1_y_xxyy_xyzz_0[i] * pa_x[i] - ta1_y_xxyy_xyzz_1[i] * pc_x[i];

        ta1_y_xxxyy_xzzz_0[i] = ta1_y_xxx_xzzz_0[i] * fe_0 - ta1_y_xxx_xzzz_1[i] * fe_0 + ta_xxxy_xzzz_1[i] + ta1_y_xxxy_xzzz_0[i] * pa_y[i] -
                                ta1_y_xxxy_xzzz_1[i] * pc_y[i];

        ta1_y_xxxyy_yyyy_0[i] =
            2.0 * ta1_y_xyy_yyyy_0[i] * fe_0 - 2.0 * ta1_y_xyy_yyyy_1[i] * fe_0 + ta1_y_xxyy_yyyy_0[i] * pa_x[i] - ta1_y_xxyy_yyyy_1[i] * pc_x[i];

        ta1_y_xxxyy_yyyz_0[i] =
            2.0 * ta1_y_xyy_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xyy_yyyz_1[i] * fe_0 + ta1_y_xxyy_yyyz_0[i] * pa_x[i] - ta1_y_xxyy_yyyz_1[i] * pc_x[i];

        ta1_y_xxxyy_yyzz_0[i] =
            2.0 * ta1_y_xyy_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xyy_yyzz_1[i] * fe_0 + ta1_y_xxyy_yyzz_0[i] * pa_x[i] - ta1_y_xxyy_yyzz_1[i] * pc_x[i];

        ta1_y_xxxyy_yzzz_0[i] =
            2.0 * ta1_y_xyy_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xyy_yzzz_1[i] * fe_0 + ta1_y_xxyy_yzzz_0[i] * pa_x[i] - ta1_y_xxyy_yzzz_1[i] * pc_x[i];

        ta1_y_xxxyy_zzzz_0[i] =
            2.0 * ta1_y_xyy_zzzz_0[i] * fe_0 - 2.0 * ta1_y_xyy_zzzz_1[i] * fe_0 + ta1_y_xxyy_zzzz_0[i] * pa_x[i] - ta1_y_xxyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 375-390 components of targeted buffer : HG

    auto ta1_y_xxxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 375);

    auto ta1_y_xxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 376);

    auto ta1_y_xxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 377);

    auto ta1_y_xxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 378);

    auto ta1_y_xxxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 379);

    auto ta1_y_xxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 380);

    auto ta1_y_xxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 381);

    auto ta1_y_xxxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 382);

    auto ta1_y_xxxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 383);

    auto ta1_y_xxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 384);

    auto ta1_y_xxxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 385);

    auto ta1_y_xxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 386);

    auto ta1_y_xxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 387);

    auto ta1_y_xxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 388);

    auto ta1_y_xxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 389);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxxy_xxxx_0,  \
                             ta1_y_xxxy_xxxx_1,  \
                             ta1_y_xxxy_xxxy_0,  \
                             ta1_y_xxxy_xxxy_1,  \
                             ta1_y_xxxy_xxy_0,   \
                             ta1_y_xxxy_xxy_1,   \
                             ta1_y_xxxy_xxyy_0,  \
                             ta1_y_xxxy_xxyy_1,  \
                             ta1_y_xxxy_xxyz_0,  \
                             ta1_y_xxxy_xxyz_1,  \
                             ta1_y_xxxy_xyy_0,   \
                             ta1_y_xxxy_xyy_1,   \
                             ta1_y_xxxy_xyyy_0,  \
                             ta1_y_xxxy_xyyy_1,  \
                             ta1_y_xxxy_xyyz_0,  \
                             ta1_y_xxxy_xyyz_1,  \
                             ta1_y_xxxy_xyz_0,   \
                             ta1_y_xxxy_xyz_1,   \
                             ta1_y_xxxy_xyzz_0,  \
                             ta1_y_xxxy_xyzz_1,  \
                             ta1_y_xxxy_yyyy_0,  \
                             ta1_y_xxxy_yyyy_1,  \
                             ta1_y_xxxyz_xxxx_0, \
                             ta1_y_xxxyz_xxxy_0, \
                             ta1_y_xxxyz_xxxz_0, \
                             ta1_y_xxxyz_xxyy_0, \
                             ta1_y_xxxyz_xxyz_0, \
                             ta1_y_xxxyz_xxzz_0, \
                             ta1_y_xxxyz_xyyy_0, \
                             ta1_y_xxxyz_xyyz_0, \
                             ta1_y_xxxyz_xyzz_0, \
                             ta1_y_xxxyz_xzzz_0, \
                             ta1_y_xxxyz_yyyy_0, \
                             ta1_y_xxxyz_yyyz_0, \
                             ta1_y_xxxyz_yyzz_0, \
                             ta1_y_xxxyz_yzzz_0, \
                             ta1_y_xxxyz_zzzz_0, \
                             ta1_y_xxxz_xxxz_0,  \
                             ta1_y_xxxz_xxxz_1,  \
                             ta1_y_xxxz_xxzz_0,  \
                             ta1_y_xxxz_xxzz_1,  \
                             ta1_y_xxxz_xzzz_0,  \
                             ta1_y_xxxz_xzzz_1,  \
                             ta1_y_xxxz_zzzz_0,  \
                             ta1_y_xxxz_zzzz_1,  \
                             ta1_y_xxyz_yyyz_0,  \
                             ta1_y_xxyz_yyyz_1,  \
                             ta1_y_xxyz_yyzz_0,  \
                             ta1_y_xxyz_yyzz_1,  \
                             ta1_y_xxyz_yzzz_0,  \
                             ta1_y_xxyz_yzzz_1,  \
                             ta1_y_xyz_yyyz_0,   \
                             ta1_y_xyz_yyyz_1,   \
                             ta1_y_xyz_yyzz_0,   \
                             ta1_y_xyz_yyzz_1,   \
                             ta1_y_xyz_yzzz_0,   \
                             ta1_y_xyz_yzzz_1,   \
                             ta_xxxz_xxxz_1,     \
                             ta_xxxz_xxzz_1,     \
                             ta_xxxz_xzzz_1,     \
                             ta_xxxz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyz_xxxx_0[i] = ta1_y_xxxy_xxxx_0[i] * pa_z[i] - ta1_y_xxxy_xxxx_1[i] * pc_z[i];

        ta1_y_xxxyz_xxxy_0[i] = ta1_y_xxxy_xxxy_0[i] * pa_z[i] - ta1_y_xxxy_xxxy_1[i] * pc_z[i];

        ta1_y_xxxyz_xxxz_0[i] = ta_xxxz_xxxz_1[i] + ta1_y_xxxz_xxxz_0[i] * pa_y[i] - ta1_y_xxxz_xxxz_1[i] * pc_y[i];

        ta1_y_xxxyz_xxyy_0[i] = ta1_y_xxxy_xxyy_0[i] * pa_z[i] - ta1_y_xxxy_xxyy_1[i] * pc_z[i];

        ta1_y_xxxyz_xxyz_0[i] =
            ta1_y_xxxy_xxy_0[i] * fe_0 - ta1_y_xxxy_xxy_1[i] * fe_0 + ta1_y_xxxy_xxyz_0[i] * pa_z[i] - ta1_y_xxxy_xxyz_1[i] * pc_z[i];

        ta1_y_xxxyz_xxzz_0[i] = ta_xxxz_xxzz_1[i] + ta1_y_xxxz_xxzz_0[i] * pa_y[i] - ta1_y_xxxz_xxzz_1[i] * pc_y[i];

        ta1_y_xxxyz_xyyy_0[i] = ta1_y_xxxy_xyyy_0[i] * pa_z[i] - ta1_y_xxxy_xyyy_1[i] * pc_z[i];

        ta1_y_xxxyz_xyyz_0[i] =
            ta1_y_xxxy_xyy_0[i] * fe_0 - ta1_y_xxxy_xyy_1[i] * fe_0 + ta1_y_xxxy_xyyz_0[i] * pa_z[i] - ta1_y_xxxy_xyyz_1[i] * pc_z[i];

        ta1_y_xxxyz_xyzz_0[i] =
            2.0 * ta1_y_xxxy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xyz_1[i] * fe_0 + ta1_y_xxxy_xyzz_0[i] * pa_z[i] - ta1_y_xxxy_xyzz_1[i] * pc_z[i];

        ta1_y_xxxyz_xzzz_0[i] = ta_xxxz_xzzz_1[i] + ta1_y_xxxz_xzzz_0[i] * pa_y[i] - ta1_y_xxxz_xzzz_1[i] * pc_y[i];

        ta1_y_xxxyz_yyyy_0[i] = ta1_y_xxxy_yyyy_0[i] * pa_z[i] - ta1_y_xxxy_yyyy_1[i] * pc_z[i];

        ta1_y_xxxyz_yyyz_0[i] =
            2.0 * ta1_y_xyz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xyz_yyyz_1[i] * fe_0 + ta1_y_xxyz_yyyz_0[i] * pa_x[i] - ta1_y_xxyz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxyz_yyzz_0[i] =
            2.0 * ta1_y_xyz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xyz_yyzz_1[i] * fe_0 + ta1_y_xxyz_yyzz_0[i] * pa_x[i] - ta1_y_xxyz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxyz_yzzz_0[i] =
            2.0 * ta1_y_xyz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xyz_yzzz_1[i] * fe_0 + ta1_y_xxyz_yzzz_0[i] * pa_x[i] - ta1_y_xxyz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxyz_zzzz_0[i] = ta_xxxz_zzzz_1[i] + ta1_y_xxxz_zzzz_0[i] * pa_y[i] - ta1_y_xxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 390-405 components of targeted buffer : HG

    auto ta1_y_xxxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 390);

    auto ta1_y_xxxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 391);

    auto ta1_y_xxxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 392);

    auto ta1_y_xxxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 393);

    auto ta1_y_xxxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 394);

    auto ta1_y_xxxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 395);

    auto ta1_y_xxxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 396);

    auto ta1_y_xxxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 397);

    auto ta1_y_xxxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 398);

    auto ta1_y_xxxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 399);

    auto ta1_y_xxxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 400);

    auto ta1_y_xxxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 401);

    auto ta1_y_xxxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 402);

    auto ta1_y_xxxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 403);

    auto ta1_y_xxxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 404);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxx_xxxx_0,   \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxy_0,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxyy_0,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xyyy_0,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxxz_xxxx_0,  \
                             ta1_y_xxxz_xxxx_1,  \
                             ta1_y_xxxz_xxxy_0,  \
                             ta1_y_xxxz_xxxy_1,  \
                             ta1_y_xxxz_xxyy_0,  \
                             ta1_y_xxxz_xxyy_1,  \
                             ta1_y_xxxz_xyyy_0,  \
                             ta1_y_xxxz_xyyy_1,  \
                             ta1_y_xxxzz_xxxx_0, \
                             ta1_y_xxxzz_xxxy_0, \
                             ta1_y_xxxzz_xxxz_0, \
                             ta1_y_xxxzz_xxyy_0, \
                             ta1_y_xxxzz_xxyz_0, \
                             ta1_y_xxxzz_xxzz_0, \
                             ta1_y_xxxzz_xyyy_0, \
                             ta1_y_xxxzz_xyyz_0, \
                             ta1_y_xxxzz_xyzz_0, \
                             ta1_y_xxxzz_xzzz_0, \
                             ta1_y_xxxzz_yyyy_0, \
                             ta1_y_xxxzz_yyyz_0, \
                             ta1_y_xxxzz_yyzz_0, \
                             ta1_y_xxxzz_yzzz_0, \
                             ta1_y_xxxzz_zzzz_0, \
                             ta1_y_xxzz_xxxz_0,  \
                             ta1_y_xxzz_xxxz_1,  \
                             ta1_y_xxzz_xxyz_0,  \
                             ta1_y_xxzz_xxyz_1,  \
                             ta1_y_xxzz_xxz_0,   \
                             ta1_y_xxzz_xxz_1,   \
                             ta1_y_xxzz_xxzz_0,  \
                             ta1_y_xxzz_xxzz_1,  \
                             ta1_y_xxzz_xyyz_0,  \
                             ta1_y_xxzz_xyyz_1,  \
                             ta1_y_xxzz_xyz_0,   \
                             ta1_y_xxzz_xyz_1,   \
                             ta1_y_xxzz_xyzz_0,  \
                             ta1_y_xxzz_xyzz_1,  \
                             ta1_y_xxzz_xzz_0,   \
                             ta1_y_xxzz_xzz_1,   \
                             ta1_y_xxzz_xzzz_0,  \
                             ta1_y_xxzz_xzzz_1,  \
                             ta1_y_xxzz_yyyy_0,  \
                             ta1_y_xxzz_yyyy_1,  \
                             ta1_y_xxzz_yyyz_0,  \
                             ta1_y_xxzz_yyyz_1,  \
                             ta1_y_xxzz_yyz_0,   \
                             ta1_y_xxzz_yyz_1,   \
                             ta1_y_xxzz_yyzz_0,  \
                             ta1_y_xxzz_yyzz_1,  \
                             ta1_y_xxzz_yzz_0,   \
                             ta1_y_xxzz_yzz_1,   \
                             ta1_y_xxzz_yzzz_0,  \
                             ta1_y_xxzz_yzzz_1,  \
                             ta1_y_xxzz_zzz_0,   \
                             ta1_y_xxzz_zzz_1,   \
                             ta1_y_xxzz_zzzz_0,  \
                             ta1_y_xxzz_zzzz_1,  \
                             ta1_y_xzz_xxxz_0,   \
                             ta1_y_xzz_xxxz_1,   \
                             ta1_y_xzz_xxyz_0,   \
                             ta1_y_xzz_xxyz_1,   \
                             ta1_y_xzz_xxzz_0,   \
                             ta1_y_xzz_xxzz_1,   \
                             ta1_y_xzz_xyyz_0,   \
                             ta1_y_xzz_xyyz_1,   \
                             ta1_y_xzz_xyzz_0,   \
                             ta1_y_xzz_xyzz_1,   \
                             ta1_y_xzz_xzzz_0,   \
                             ta1_y_xzz_xzzz_1,   \
                             ta1_y_xzz_yyyy_0,   \
                             ta1_y_xzz_yyyy_1,   \
                             ta1_y_xzz_yyyz_0,   \
                             ta1_y_xzz_yyyz_1,   \
                             ta1_y_xzz_yyzz_0,   \
                             ta1_y_xzz_yyzz_1,   \
                             ta1_y_xzz_yzzz_0,   \
                             ta1_y_xzz_yzzz_1,   \
                             ta1_y_xzz_zzzz_0,   \
                             ta1_y_xzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzz_xxxx_0[i] =
            ta1_y_xxx_xxxx_0[i] * fe_0 - ta1_y_xxx_xxxx_1[i] * fe_0 + ta1_y_xxxz_xxxx_0[i] * pa_z[i] - ta1_y_xxxz_xxxx_1[i] * pc_z[i];

        ta1_y_xxxzz_xxxy_0[i] =
            ta1_y_xxx_xxxy_0[i] * fe_0 - ta1_y_xxx_xxxy_1[i] * fe_0 + ta1_y_xxxz_xxxy_0[i] * pa_z[i] - ta1_y_xxxz_xxxy_1[i] * pc_z[i];

        ta1_y_xxxzz_xxxz_0[i] = 2.0 * ta1_y_xzz_xxxz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xxzz_xxz_0[i] * fe_0 -
                                3.0 * ta1_y_xxzz_xxz_1[i] * fe_0 + ta1_y_xxzz_xxxz_0[i] * pa_x[i] - ta1_y_xxzz_xxxz_1[i] * pc_x[i];

        ta1_y_xxxzz_xxyy_0[i] =
            ta1_y_xxx_xxyy_0[i] * fe_0 - ta1_y_xxx_xxyy_1[i] * fe_0 + ta1_y_xxxz_xxyy_0[i] * pa_z[i] - ta1_y_xxxz_xxyy_1[i] * pc_z[i];

        ta1_y_xxxzz_xxyz_0[i] = 2.0 * ta1_y_xzz_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_xxzz_xyz_1[i] * fe_0 + ta1_y_xxzz_xxyz_0[i] * pa_x[i] - ta1_y_xxzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxxzz_xxzz_0[i] = 2.0 * ta1_y_xzz_xxzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxzz_xzz_0[i] * fe_0 -
                                2.0 * ta1_y_xxzz_xzz_1[i] * fe_0 + ta1_y_xxzz_xxzz_0[i] * pa_x[i] - ta1_y_xxzz_xxzz_1[i] * pc_x[i];

        ta1_y_xxxzz_xyyy_0[i] =
            ta1_y_xxx_xyyy_0[i] * fe_0 - ta1_y_xxx_xyyy_1[i] * fe_0 + ta1_y_xxxz_xyyy_0[i] * pa_z[i] - ta1_y_xxxz_xyyy_1[i] * pc_z[i];

        ta1_y_xxxzz_xyyz_0[i] = 2.0 * ta1_y_xzz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xyyz_1[i] * fe_0 + ta1_y_xxzz_yyz_0[i] * fe_0 -
                                ta1_y_xxzz_yyz_1[i] * fe_0 + ta1_y_xxzz_xyyz_0[i] * pa_x[i] - ta1_y_xxzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxxzz_xyzz_0[i] = 2.0 * ta1_y_xzz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xyzz_1[i] * fe_0 + ta1_y_xxzz_yzz_0[i] * fe_0 -
                                ta1_y_xxzz_yzz_1[i] * fe_0 + ta1_y_xxzz_xyzz_0[i] * pa_x[i] - ta1_y_xxzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxxzz_xzzz_0[i] = 2.0 * ta1_y_xzz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xzzz_1[i] * fe_0 + ta1_y_xxzz_zzz_0[i] * fe_0 -
                                ta1_y_xxzz_zzz_1[i] * fe_0 + ta1_y_xxzz_xzzz_0[i] * pa_x[i] - ta1_y_xxzz_xzzz_1[i] * pc_x[i];

        ta1_y_xxxzz_yyyy_0[i] =
            2.0 * ta1_y_xzz_yyyy_0[i] * fe_0 - 2.0 * ta1_y_xzz_yyyy_1[i] * fe_0 + ta1_y_xxzz_yyyy_0[i] * pa_x[i] - ta1_y_xxzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxxzz_yyyz_0[i] =
            2.0 * ta1_y_xzz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xzz_yyyz_1[i] * fe_0 + ta1_y_xxzz_yyyz_0[i] * pa_x[i] - ta1_y_xxzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxzz_yyzz_0[i] =
            2.0 * ta1_y_xzz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_yyzz_1[i] * fe_0 + ta1_y_xxzz_yyzz_0[i] * pa_x[i] - ta1_y_xxzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxzz_yzzz_0[i] =
            2.0 * ta1_y_xzz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_yzzz_1[i] * fe_0 + ta1_y_xxzz_yzzz_0[i] * pa_x[i] - ta1_y_xxzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxzz_zzzz_0[i] =
            2.0 * ta1_y_xzz_zzzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_zzzz_1[i] * fe_0 + ta1_y_xxzz_zzzz_0[i] * pa_x[i] - ta1_y_xxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 405-420 components of targeted buffer : HG

    auto ta1_y_xxyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 405);

    auto ta1_y_xxyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 406);

    auto ta1_y_xxyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 407);

    auto ta1_y_xxyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 408);

    auto ta1_y_xxyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 409);

    auto ta1_y_xxyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 410);

    auto ta1_y_xxyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 411);

    auto ta1_y_xxyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 412);

    auto ta1_y_xxyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 413);

    auto ta1_y_xxyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 414);

    auto ta1_y_xxyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 415);

    auto ta1_y_xxyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 416);

    auto ta1_y_xxyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 417);

    auto ta1_y_xxyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 418);

    auto ta1_y_xxyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 419);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxy_xxxx_0,   \
                             ta1_y_xxy_xxxx_1,   \
                             ta1_y_xxy_xxxz_0,   \
                             ta1_y_xxy_xxxz_1,   \
                             ta1_y_xxy_xxzz_0,   \
                             ta1_y_xxy_xxzz_1,   \
                             ta1_y_xxy_xzzz_0,   \
                             ta1_y_xxy_xzzz_1,   \
                             ta1_y_xxyy_xxxx_0,  \
                             ta1_y_xxyy_xxxx_1,  \
                             ta1_y_xxyy_xxxz_0,  \
                             ta1_y_xxyy_xxxz_1,  \
                             ta1_y_xxyy_xxzz_0,  \
                             ta1_y_xxyy_xxzz_1,  \
                             ta1_y_xxyy_xzzz_0,  \
                             ta1_y_xxyy_xzzz_1,  \
                             ta1_y_xxyyy_xxxx_0, \
                             ta1_y_xxyyy_xxxy_0, \
                             ta1_y_xxyyy_xxxz_0, \
                             ta1_y_xxyyy_xxyy_0, \
                             ta1_y_xxyyy_xxyz_0, \
                             ta1_y_xxyyy_xxzz_0, \
                             ta1_y_xxyyy_xyyy_0, \
                             ta1_y_xxyyy_xyyz_0, \
                             ta1_y_xxyyy_xyzz_0, \
                             ta1_y_xxyyy_xzzz_0, \
                             ta1_y_xxyyy_yyyy_0, \
                             ta1_y_xxyyy_yyyz_0, \
                             ta1_y_xxyyy_yyzz_0, \
                             ta1_y_xxyyy_yzzz_0, \
                             ta1_y_xxyyy_zzzz_0, \
                             ta1_y_xyyy_xxxy_0,  \
                             ta1_y_xyyy_xxxy_1,  \
                             ta1_y_xyyy_xxy_0,   \
                             ta1_y_xyyy_xxy_1,   \
                             ta1_y_xyyy_xxyy_0,  \
                             ta1_y_xyyy_xxyy_1,  \
                             ta1_y_xyyy_xxyz_0,  \
                             ta1_y_xyyy_xxyz_1,  \
                             ta1_y_xyyy_xyy_0,   \
                             ta1_y_xyyy_xyy_1,   \
                             ta1_y_xyyy_xyyy_0,  \
                             ta1_y_xyyy_xyyy_1,  \
                             ta1_y_xyyy_xyyz_0,  \
                             ta1_y_xyyy_xyyz_1,  \
                             ta1_y_xyyy_xyz_0,   \
                             ta1_y_xyyy_xyz_1,   \
                             ta1_y_xyyy_xyzz_0,  \
                             ta1_y_xyyy_xyzz_1,  \
                             ta1_y_xyyy_yyy_0,   \
                             ta1_y_xyyy_yyy_1,   \
                             ta1_y_xyyy_yyyy_0,  \
                             ta1_y_xyyy_yyyy_1,  \
                             ta1_y_xyyy_yyyz_0,  \
                             ta1_y_xyyy_yyyz_1,  \
                             ta1_y_xyyy_yyz_0,   \
                             ta1_y_xyyy_yyz_1,   \
                             ta1_y_xyyy_yyzz_0,  \
                             ta1_y_xyyy_yyzz_1,  \
                             ta1_y_xyyy_yzz_0,   \
                             ta1_y_xyyy_yzz_1,   \
                             ta1_y_xyyy_yzzz_0,  \
                             ta1_y_xyyy_yzzz_1,  \
                             ta1_y_xyyy_zzzz_0,  \
                             ta1_y_xyyy_zzzz_1,  \
                             ta1_y_yyy_xxxy_0,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxyy_0,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyz_0,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xyyy_0,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyz_0,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyzz_0,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_yyyy_0,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyz_0,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyzz_0,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yzzz_0,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_zzzz_0,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta_xxyy_xxxx_1,     \
                             ta_xxyy_xxxz_1,     \
                             ta_xxyy_xxzz_1,     \
                             ta_xxyy_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyy_xxxx_0[i] = 2.0 * ta1_y_xxy_xxxx_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxxx_1[i] * fe_0 + ta_xxyy_xxxx_1[i] +
                                ta1_y_xxyy_xxxx_0[i] * pa_y[i] - ta1_y_xxyy_xxxx_1[i] * pc_y[i];

        ta1_y_xxyyy_xxxy_0[i] = ta1_y_yyy_xxxy_0[i] * fe_0 - ta1_y_yyy_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xyyy_xxy_0[i] * fe_0 -
                                3.0 * ta1_y_xyyy_xxy_1[i] * fe_0 + ta1_y_xyyy_xxxy_0[i] * pa_x[i] - ta1_y_xyyy_xxxy_1[i] * pc_x[i];

        ta1_y_xxyyy_xxxz_0[i] = 2.0 * ta1_y_xxy_xxxz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxxz_1[i] * fe_0 + ta_xxyy_xxxz_1[i] +
                                ta1_y_xxyy_xxxz_0[i] * pa_y[i] - ta1_y_xxyy_xxxz_1[i] * pc_y[i];

        ta1_y_xxyyy_xxyy_0[i] = ta1_y_yyy_xxyy_0[i] * fe_0 - ta1_y_yyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xyyy_xyy_0[i] * fe_0 -
                                2.0 * ta1_y_xyyy_xyy_1[i] * fe_0 + ta1_y_xyyy_xxyy_0[i] * pa_x[i] - ta1_y_xyyy_xxyy_1[i] * pc_x[i];

        ta1_y_xxyyy_xxyz_0[i] = ta1_y_yyy_xxyz_0[i] * fe_0 - ta1_y_yyy_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xyyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_xyyy_xyz_1[i] * fe_0 + ta1_y_xyyy_xxyz_0[i] * pa_x[i] - ta1_y_xyyy_xxyz_1[i] * pc_x[i];

        ta1_y_xxyyy_xxzz_0[i] = 2.0 * ta1_y_xxy_xxzz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxzz_1[i] * fe_0 + ta_xxyy_xxzz_1[i] +
                                ta1_y_xxyy_xxzz_0[i] * pa_y[i] - ta1_y_xxyy_xxzz_1[i] * pc_y[i];

        ta1_y_xxyyy_xyyy_0[i] = ta1_y_yyy_xyyy_0[i] * fe_0 - ta1_y_yyy_xyyy_1[i] * fe_0 + ta1_y_xyyy_yyy_0[i] * fe_0 - ta1_y_xyyy_yyy_1[i] * fe_0 +
                                ta1_y_xyyy_xyyy_0[i] * pa_x[i] - ta1_y_xyyy_xyyy_1[i] * pc_x[i];

        ta1_y_xxyyy_xyyz_0[i] = ta1_y_yyy_xyyz_0[i] * fe_0 - ta1_y_yyy_xyyz_1[i] * fe_0 + ta1_y_xyyy_yyz_0[i] * fe_0 - ta1_y_xyyy_yyz_1[i] * fe_0 +
                                ta1_y_xyyy_xyyz_0[i] * pa_x[i] - ta1_y_xyyy_xyyz_1[i] * pc_x[i];

        ta1_y_xxyyy_xyzz_0[i] = ta1_y_yyy_xyzz_0[i] * fe_0 - ta1_y_yyy_xyzz_1[i] * fe_0 + ta1_y_xyyy_yzz_0[i] * fe_0 - ta1_y_xyyy_yzz_1[i] * fe_0 +
                                ta1_y_xyyy_xyzz_0[i] * pa_x[i] - ta1_y_xyyy_xyzz_1[i] * pc_x[i];

        ta1_y_xxyyy_xzzz_0[i] = 2.0 * ta1_y_xxy_xzzz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xzzz_1[i] * fe_0 + ta_xxyy_xzzz_1[i] +
                                ta1_y_xxyy_xzzz_0[i] * pa_y[i] - ta1_y_xxyy_xzzz_1[i] * pc_y[i];

        ta1_y_xxyyy_yyyy_0[i] =
            ta1_y_yyy_yyyy_0[i] * fe_0 - ta1_y_yyy_yyyy_1[i] * fe_0 + ta1_y_xyyy_yyyy_0[i] * pa_x[i] - ta1_y_xyyy_yyyy_1[i] * pc_x[i];

        ta1_y_xxyyy_yyyz_0[i] =
            ta1_y_yyy_yyyz_0[i] * fe_0 - ta1_y_yyy_yyyz_1[i] * fe_0 + ta1_y_xyyy_yyyz_0[i] * pa_x[i] - ta1_y_xyyy_yyyz_1[i] * pc_x[i];

        ta1_y_xxyyy_yyzz_0[i] =
            ta1_y_yyy_yyzz_0[i] * fe_0 - ta1_y_yyy_yyzz_1[i] * fe_0 + ta1_y_xyyy_yyzz_0[i] * pa_x[i] - ta1_y_xyyy_yyzz_1[i] * pc_x[i];

        ta1_y_xxyyy_yzzz_0[i] =
            ta1_y_yyy_yzzz_0[i] * fe_0 - ta1_y_yyy_yzzz_1[i] * fe_0 + ta1_y_xyyy_yzzz_0[i] * pa_x[i] - ta1_y_xyyy_yzzz_1[i] * pc_x[i];

        ta1_y_xxyyy_zzzz_0[i] =
            ta1_y_yyy_zzzz_0[i] * fe_0 - ta1_y_yyy_zzzz_1[i] * fe_0 + ta1_y_xyyy_zzzz_0[i] * pa_x[i] - ta1_y_xyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 420-435 components of targeted buffer : HG

    auto ta1_y_xxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 420);

    auto ta1_y_xxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 421);

    auto ta1_y_xxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 422);

    auto ta1_y_xxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 423);

    auto ta1_y_xxyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 424);

    auto ta1_y_xxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 425);

    auto ta1_y_xxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 426);

    auto ta1_y_xxyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 427);

    auto ta1_y_xxyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 428);

    auto ta1_y_xxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 429);

    auto ta1_y_xxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 430);

    auto ta1_y_xxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 431);

    auto ta1_y_xxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 432);

    auto ta1_y_xxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 433);

    auto ta1_y_xxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 434);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxyy_xxx_0,   \
                             ta1_y_xxyy_xxx_1,   \
                             ta1_y_xxyy_xxxx_0,  \
                             ta1_y_xxyy_xxxx_1,  \
                             ta1_y_xxyy_xxxy_0,  \
                             ta1_y_xxyy_xxxy_1,  \
                             ta1_y_xxyy_xxxz_0,  \
                             ta1_y_xxyy_xxxz_1,  \
                             ta1_y_xxyy_xxy_0,   \
                             ta1_y_xxyy_xxy_1,   \
                             ta1_y_xxyy_xxyy_0,  \
                             ta1_y_xxyy_xxyy_1,  \
                             ta1_y_xxyy_xxyz_0,  \
                             ta1_y_xxyy_xxyz_1,  \
                             ta1_y_xxyy_xxz_0,   \
                             ta1_y_xxyy_xxz_1,   \
                             ta1_y_xxyy_xxzz_0,  \
                             ta1_y_xxyy_xxzz_1,  \
                             ta1_y_xxyy_xyy_0,   \
                             ta1_y_xxyy_xyy_1,   \
                             ta1_y_xxyy_xyyy_0,  \
                             ta1_y_xxyy_xyyy_1,  \
                             ta1_y_xxyy_xyyz_0,  \
                             ta1_y_xxyy_xyyz_1,  \
                             ta1_y_xxyy_xyz_0,   \
                             ta1_y_xxyy_xyz_1,   \
                             ta1_y_xxyy_xyzz_0,  \
                             ta1_y_xxyy_xyzz_1,  \
                             ta1_y_xxyy_xzz_0,   \
                             ta1_y_xxyy_xzz_1,   \
                             ta1_y_xxyy_xzzz_0,  \
                             ta1_y_xxyy_xzzz_1,  \
                             ta1_y_xxyy_yyyy_0,  \
                             ta1_y_xxyy_yyyy_1,  \
                             ta1_y_xxyyz_xxxx_0, \
                             ta1_y_xxyyz_xxxy_0, \
                             ta1_y_xxyyz_xxxz_0, \
                             ta1_y_xxyyz_xxyy_0, \
                             ta1_y_xxyyz_xxyz_0, \
                             ta1_y_xxyyz_xxzz_0, \
                             ta1_y_xxyyz_xyyy_0, \
                             ta1_y_xxyyz_xyyz_0, \
                             ta1_y_xxyyz_xyzz_0, \
                             ta1_y_xxyyz_xzzz_0, \
                             ta1_y_xxyyz_yyyy_0, \
                             ta1_y_xxyyz_yyyz_0, \
                             ta1_y_xxyyz_yyzz_0, \
                             ta1_y_xxyyz_yzzz_0, \
                             ta1_y_xxyyz_zzzz_0, \
                             ta1_y_xyyz_yyyz_0,  \
                             ta1_y_xyyz_yyyz_1,  \
                             ta1_y_xyyz_yyzz_0,  \
                             ta1_y_xyyz_yyzz_1,  \
                             ta1_y_xyyz_yzzz_0,  \
                             ta1_y_xyyz_yzzz_1,  \
                             ta1_y_xyyz_zzzz_0,  \
                             ta1_y_xyyz_zzzz_1,  \
                             ta1_y_yyz_yyyz_0,   \
                             ta1_y_yyz_yyyz_1,   \
                             ta1_y_yyz_yyzz_0,   \
                             ta1_y_yyz_yyzz_1,   \
                             ta1_y_yyz_yzzz_0,   \
                             ta1_y_yyz_yzzz_1,   \
                             ta1_y_yyz_zzzz_0,   \
                             ta1_y_yyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyz_xxxx_0[i] = ta1_y_xxyy_xxxx_0[i] * pa_z[i] - ta1_y_xxyy_xxxx_1[i] * pc_z[i];

        ta1_y_xxyyz_xxxy_0[i] = ta1_y_xxyy_xxxy_0[i] * pa_z[i] - ta1_y_xxyy_xxxy_1[i] * pc_z[i];

        ta1_y_xxyyz_xxxz_0[i] =
            ta1_y_xxyy_xxx_0[i] * fe_0 - ta1_y_xxyy_xxx_1[i] * fe_0 + ta1_y_xxyy_xxxz_0[i] * pa_z[i] - ta1_y_xxyy_xxxz_1[i] * pc_z[i];

        ta1_y_xxyyz_xxyy_0[i] = ta1_y_xxyy_xxyy_0[i] * pa_z[i] - ta1_y_xxyy_xxyy_1[i] * pc_z[i];

        ta1_y_xxyyz_xxyz_0[i] =
            ta1_y_xxyy_xxy_0[i] * fe_0 - ta1_y_xxyy_xxy_1[i] * fe_0 + ta1_y_xxyy_xxyz_0[i] * pa_z[i] - ta1_y_xxyy_xxyz_1[i] * pc_z[i];

        ta1_y_xxyyz_xxzz_0[i] =
            2.0 * ta1_y_xxyy_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxyy_xxz_1[i] * fe_0 + ta1_y_xxyy_xxzz_0[i] * pa_z[i] - ta1_y_xxyy_xxzz_1[i] * pc_z[i];

        ta1_y_xxyyz_xyyy_0[i] = ta1_y_xxyy_xyyy_0[i] * pa_z[i] - ta1_y_xxyy_xyyy_1[i] * pc_z[i];

        ta1_y_xxyyz_xyyz_0[i] =
            ta1_y_xxyy_xyy_0[i] * fe_0 - ta1_y_xxyy_xyy_1[i] * fe_0 + ta1_y_xxyy_xyyz_0[i] * pa_z[i] - ta1_y_xxyy_xyyz_1[i] * pc_z[i];

        ta1_y_xxyyz_xyzz_0[i] =
            2.0 * ta1_y_xxyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxyy_xyz_1[i] * fe_0 + ta1_y_xxyy_xyzz_0[i] * pa_z[i] - ta1_y_xxyy_xyzz_1[i] * pc_z[i];

        ta1_y_xxyyz_xzzz_0[i] =
            3.0 * ta1_y_xxyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xzz_1[i] * fe_0 + ta1_y_xxyy_xzzz_0[i] * pa_z[i] - ta1_y_xxyy_xzzz_1[i] * pc_z[i];

        ta1_y_xxyyz_yyyy_0[i] = ta1_y_xxyy_yyyy_0[i] * pa_z[i] - ta1_y_xxyy_yyyy_1[i] * pc_z[i];

        ta1_y_xxyyz_yyyz_0[i] =
            ta1_y_yyz_yyyz_0[i] * fe_0 - ta1_y_yyz_yyyz_1[i] * fe_0 + ta1_y_xyyz_yyyz_0[i] * pa_x[i] - ta1_y_xyyz_yyyz_1[i] * pc_x[i];

        ta1_y_xxyyz_yyzz_0[i] =
            ta1_y_yyz_yyzz_0[i] * fe_0 - ta1_y_yyz_yyzz_1[i] * fe_0 + ta1_y_xyyz_yyzz_0[i] * pa_x[i] - ta1_y_xyyz_yyzz_1[i] * pc_x[i];

        ta1_y_xxyyz_yzzz_0[i] =
            ta1_y_yyz_yzzz_0[i] * fe_0 - ta1_y_yyz_yzzz_1[i] * fe_0 + ta1_y_xyyz_yzzz_0[i] * pa_x[i] - ta1_y_xyyz_yzzz_1[i] * pc_x[i];

        ta1_y_xxyyz_zzzz_0[i] =
            ta1_y_yyz_zzzz_0[i] * fe_0 - ta1_y_yyz_zzzz_1[i] * fe_0 + ta1_y_xyyz_zzzz_0[i] * pa_x[i] - ta1_y_xyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 435-450 components of targeted buffer : HG

    auto ta1_y_xxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 435);

    auto ta1_y_xxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 436);

    auto ta1_y_xxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 437);

    auto ta1_y_xxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 438);

    auto ta1_y_xxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 439);

    auto ta1_y_xxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 440);

    auto ta1_y_xxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 441);

    auto ta1_y_xxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 442);

    auto ta1_y_xxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 443);

    auto ta1_y_xxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 444);

    auto ta1_y_xxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 445);

    auto ta1_y_xxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 446);

    auto ta1_y_xxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 447);

    auto ta1_y_xxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 448);

    auto ta1_y_xxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 449);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxy_xxxy_0,   \
                             ta1_y_xxy_xxxy_1,   \
                             ta1_y_xxy_xxyy_0,   \
                             ta1_y_xxy_xxyy_1,   \
                             ta1_y_xxy_xyyy_0,   \
                             ta1_y_xxy_xyyy_1,   \
                             ta1_y_xxyz_xxxy_0,  \
                             ta1_y_xxyz_xxxy_1,  \
                             ta1_y_xxyz_xxyy_0,  \
                             ta1_y_xxyz_xxyy_1,  \
                             ta1_y_xxyz_xyyy_0,  \
                             ta1_y_xxyz_xyyy_1,  \
                             ta1_y_xxyzz_xxxx_0, \
                             ta1_y_xxyzz_xxxy_0, \
                             ta1_y_xxyzz_xxxz_0, \
                             ta1_y_xxyzz_xxyy_0, \
                             ta1_y_xxyzz_xxyz_0, \
                             ta1_y_xxyzz_xxzz_0, \
                             ta1_y_xxyzz_xyyy_0, \
                             ta1_y_xxyzz_xyyz_0, \
                             ta1_y_xxyzz_xyzz_0, \
                             ta1_y_xxyzz_xzzz_0, \
                             ta1_y_xxyzz_yyyy_0, \
                             ta1_y_xxyzz_yyyz_0, \
                             ta1_y_xxyzz_yyzz_0, \
                             ta1_y_xxyzz_yzzz_0, \
                             ta1_y_xxyzz_zzzz_0, \
                             ta1_y_xxzz_xxxx_0,  \
                             ta1_y_xxzz_xxxx_1,  \
                             ta1_y_xxzz_xxxz_0,  \
                             ta1_y_xxzz_xxxz_1,  \
                             ta1_y_xxzz_xxyz_0,  \
                             ta1_y_xxzz_xxyz_1,  \
                             ta1_y_xxzz_xxz_0,   \
                             ta1_y_xxzz_xxz_1,   \
                             ta1_y_xxzz_xxzz_0,  \
                             ta1_y_xxzz_xxzz_1,  \
                             ta1_y_xxzz_xyyz_0,  \
                             ta1_y_xxzz_xyyz_1,  \
                             ta1_y_xxzz_xyz_0,   \
                             ta1_y_xxzz_xyz_1,   \
                             ta1_y_xxzz_xyzz_0,  \
                             ta1_y_xxzz_xyzz_1,  \
                             ta1_y_xxzz_xzz_0,   \
                             ta1_y_xxzz_xzz_1,   \
                             ta1_y_xxzz_xzzz_0,  \
                             ta1_y_xxzz_xzzz_1,  \
                             ta1_y_xxzz_zzzz_0,  \
                             ta1_y_xxzz_zzzz_1,  \
                             ta1_y_xyzz_yyyy_0,  \
                             ta1_y_xyzz_yyyy_1,  \
                             ta1_y_xyzz_yyyz_0,  \
                             ta1_y_xyzz_yyyz_1,  \
                             ta1_y_xyzz_yyzz_0,  \
                             ta1_y_xyzz_yyzz_1,  \
                             ta1_y_xyzz_yzzz_0,  \
                             ta1_y_xyzz_yzzz_1,  \
                             ta1_y_yzz_yyyy_0,   \
                             ta1_y_yzz_yyyy_1,   \
                             ta1_y_yzz_yyyz_0,   \
                             ta1_y_yzz_yyyz_1,   \
                             ta1_y_yzz_yyzz_0,   \
                             ta1_y_yzz_yyzz_1,   \
                             ta1_y_yzz_yzzz_0,   \
                             ta1_y_yzz_yzzz_1,   \
                             ta_xxzz_xxxx_1,     \
                             ta_xxzz_xxxz_1,     \
                             ta_xxzz_xxyz_1,     \
                             ta_xxzz_xxzz_1,     \
                             ta_xxzz_xyyz_1,     \
                             ta_xxzz_xyzz_1,     \
                             ta_xxzz_xzzz_1,     \
                             ta_xxzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzz_xxxx_0[i] = ta_xxzz_xxxx_1[i] + ta1_y_xxzz_xxxx_0[i] * pa_y[i] - ta1_y_xxzz_xxxx_1[i] * pc_y[i];

        ta1_y_xxyzz_xxxy_0[i] =
            ta1_y_xxy_xxxy_0[i] * fe_0 - ta1_y_xxy_xxxy_1[i] * fe_0 + ta1_y_xxyz_xxxy_0[i] * pa_z[i] - ta1_y_xxyz_xxxy_1[i] * pc_z[i];

        ta1_y_xxyzz_xxxz_0[i] = ta_xxzz_xxxz_1[i] + ta1_y_xxzz_xxxz_0[i] * pa_y[i] - ta1_y_xxzz_xxxz_1[i] * pc_y[i];

        ta1_y_xxyzz_xxyy_0[i] =
            ta1_y_xxy_xxyy_0[i] * fe_0 - ta1_y_xxy_xxyy_1[i] * fe_0 + ta1_y_xxyz_xxyy_0[i] * pa_z[i] - ta1_y_xxyz_xxyy_1[i] * pc_z[i];

        ta1_y_xxyzz_xxyz_0[i] = ta1_y_xxzz_xxz_0[i] * fe_0 - ta1_y_xxzz_xxz_1[i] * fe_0 + ta_xxzz_xxyz_1[i] + ta1_y_xxzz_xxyz_0[i] * pa_y[i] -
                                ta1_y_xxzz_xxyz_1[i] * pc_y[i];

        ta1_y_xxyzz_xxzz_0[i] = ta_xxzz_xxzz_1[i] + ta1_y_xxzz_xxzz_0[i] * pa_y[i] - ta1_y_xxzz_xxzz_1[i] * pc_y[i];

        ta1_y_xxyzz_xyyy_0[i] =
            ta1_y_xxy_xyyy_0[i] * fe_0 - ta1_y_xxy_xyyy_1[i] * fe_0 + ta1_y_xxyz_xyyy_0[i] * pa_z[i] - ta1_y_xxyz_xyyy_1[i] * pc_z[i];

        ta1_y_xxyzz_xyyz_0[i] = 2.0 * ta1_y_xxzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxzz_xyz_1[i] * fe_0 + ta_xxzz_xyyz_1[i] +
                                ta1_y_xxzz_xyyz_0[i] * pa_y[i] - ta1_y_xxzz_xyyz_1[i] * pc_y[i];

        ta1_y_xxyzz_xyzz_0[i] = ta1_y_xxzz_xzz_0[i] * fe_0 - ta1_y_xxzz_xzz_1[i] * fe_0 + ta_xxzz_xyzz_1[i] + ta1_y_xxzz_xyzz_0[i] * pa_y[i] -
                                ta1_y_xxzz_xyzz_1[i] * pc_y[i];

        ta1_y_xxyzz_xzzz_0[i] = ta_xxzz_xzzz_1[i] + ta1_y_xxzz_xzzz_0[i] * pa_y[i] - ta1_y_xxzz_xzzz_1[i] * pc_y[i];

        ta1_y_xxyzz_yyyy_0[i] =
            ta1_y_yzz_yyyy_0[i] * fe_0 - ta1_y_yzz_yyyy_1[i] * fe_0 + ta1_y_xyzz_yyyy_0[i] * pa_x[i] - ta1_y_xyzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxyzz_yyyz_0[i] =
            ta1_y_yzz_yyyz_0[i] * fe_0 - ta1_y_yzz_yyyz_1[i] * fe_0 + ta1_y_xyzz_yyyz_0[i] * pa_x[i] - ta1_y_xyzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxyzz_yyzz_0[i] =
            ta1_y_yzz_yyzz_0[i] * fe_0 - ta1_y_yzz_yyzz_1[i] * fe_0 + ta1_y_xyzz_yyzz_0[i] * pa_x[i] - ta1_y_xyzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxyzz_yzzz_0[i] =
            ta1_y_yzz_yzzz_0[i] * fe_0 - ta1_y_yzz_yzzz_1[i] * fe_0 + ta1_y_xyzz_yzzz_0[i] * pa_x[i] - ta1_y_xyzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxyzz_zzzz_0[i] = ta_xxzz_zzzz_1[i] + ta1_y_xxzz_zzzz_0[i] * pa_y[i] - ta1_y_xxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 450-465 components of targeted buffer : HG

    auto ta1_y_xxzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 450);

    auto ta1_y_xxzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 451);

    auto ta1_y_xxzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 452);

    auto ta1_y_xxzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 453);

    auto ta1_y_xxzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 454);

    auto ta1_y_xxzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 455);

    auto ta1_y_xxzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 456);

    auto ta1_y_xxzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 457);

    auto ta1_y_xxzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 458);

    auto ta1_y_xxzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 459);

    auto ta1_y_xxzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 460);

    auto ta1_y_xxzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 461);

    auto ta1_y_xxzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 462);

    auto ta1_y_xxzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 463);

    auto ta1_y_xxzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 464);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxz_xxxx_0,   \
                             ta1_y_xxz_xxxx_1,   \
                             ta1_y_xxz_xxxy_0,   \
                             ta1_y_xxz_xxxy_1,   \
                             ta1_y_xxz_xxyy_0,   \
                             ta1_y_xxz_xxyy_1,   \
                             ta1_y_xxz_xyyy_0,   \
                             ta1_y_xxz_xyyy_1,   \
                             ta1_y_xxzz_xxxx_0,  \
                             ta1_y_xxzz_xxxx_1,  \
                             ta1_y_xxzz_xxxy_0,  \
                             ta1_y_xxzz_xxxy_1,  \
                             ta1_y_xxzz_xxyy_0,  \
                             ta1_y_xxzz_xxyy_1,  \
                             ta1_y_xxzz_xyyy_0,  \
                             ta1_y_xxzz_xyyy_1,  \
                             ta1_y_xxzzz_xxxx_0, \
                             ta1_y_xxzzz_xxxy_0, \
                             ta1_y_xxzzz_xxxz_0, \
                             ta1_y_xxzzz_xxyy_0, \
                             ta1_y_xxzzz_xxyz_0, \
                             ta1_y_xxzzz_xxzz_0, \
                             ta1_y_xxzzz_xyyy_0, \
                             ta1_y_xxzzz_xyyz_0, \
                             ta1_y_xxzzz_xyzz_0, \
                             ta1_y_xxzzz_xzzz_0, \
                             ta1_y_xxzzz_yyyy_0, \
                             ta1_y_xxzzz_yyyz_0, \
                             ta1_y_xxzzz_yyzz_0, \
                             ta1_y_xxzzz_yzzz_0, \
                             ta1_y_xxzzz_zzzz_0, \
                             ta1_y_xzzz_xxxz_0,  \
                             ta1_y_xzzz_xxxz_1,  \
                             ta1_y_xzzz_xxyz_0,  \
                             ta1_y_xzzz_xxyz_1,  \
                             ta1_y_xzzz_xxz_0,   \
                             ta1_y_xzzz_xxz_1,   \
                             ta1_y_xzzz_xxzz_0,  \
                             ta1_y_xzzz_xxzz_1,  \
                             ta1_y_xzzz_xyyz_0,  \
                             ta1_y_xzzz_xyyz_1,  \
                             ta1_y_xzzz_xyz_0,   \
                             ta1_y_xzzz_xyz_1,   \
                             ta1_y_xzzz_xyzz_0,  \
                             ta1_y_xzzz_xyzz_1,  \
                             ta1_y_xzzz_xzz_0,   \
                             ta1_y_xzzz_xzz_1,   \
                             ta1_y_xzzz_xzzz_0,  \
                             ta1_y_xzzz_xzzz_1,  \
                             ta1_y_xzzz_yyyy_0,  \
                             ta1_y_xzzz_yyyy_1,  \
                             ta1_y_xzzz_yyyz_0,  \
                             ta1_y_xzzz_yyyz_1,  \
                             ta1_y_xzzz_yyz_0,   \
                             ta1_y_xzzz_yyz_1,   \
                             ta1_y_xzzz_yyzz_0,  \
                             ta1_y_xzzz_yyzz_1,  \
                             ta1_y_xzzz_yzz_0,   \
                             ta1_y_xzzz_yzz_1,   \
                             ta1_y_xzzz_yzzz_0,  \
                             ta1_y_xzzz_yzzz_1,  \
                             ta1_y_xzzz_zzz_0,   \
                             ta1_y_xzzz_zzz_1,   \
                             ta1_y_xzzz_zzzz_0,  \
                             ta1_y_xzzz_zzzz_1,  \
                             ta1_y_zzz_xxxz_0,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxyz_0,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxzz_0,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xyyz_0,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyzz_0,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xzzz_0,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_yyyy_0,   \
                             ta1_y_zzz_yyyy_1,   \
                             ta1_y_zzz_yyyz_0,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyzz_0,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yzzz_0,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_zzzz_0,   \
                             ta1_y_zzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzz_xxxx_0[i] =
            2.0 * ta1_y_xxz_xxxx_0[i] * fe_0 - 2.0 * ta1_y_xxz_xxxx_1[i] * fe_0 + ta1_y_xxzz_xxxx_0[i] * pa_z[i] - ta1_y_xxzz_xxxx_1[i] * pc_z[i];

        ta1_y_xxzzz_xxxy_0[i] =
            2.0 * ta1_y_xxz_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xxz_xxxy_1[i] * fe_0 + ta1_y_xxzz_xxxy_0[i] * pa_z[i] - ta1_y_xxzz_xxxy_1[i] * pc_z[i];

        ta1_y_xxzzz_xxxz_0[i] = ta1_y_zzz_xxxz_0[i] * fe_0 - ta1_y_zzz_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xzzz_xxz_0[i] * fe_0 -
                                3.0 * ta1_y_xzzz_xxz_1[i] * fe_0 + ta1_y_xzzz_xxxz_0[i] * pa_x[i] - ta1_y_xzzz_xxxz_1[i] * pc_x[i];

        ta1_y_xxzzz_xxyy_0[i] =
            2.0 * ta1_y_xxz_xxyy_0[i] * fe_0 - 2.0 * ta1_y_xxz_xxyy_1[i] * fe_0 + ta1_y_xxzz_xxyy_0[i] * pa_z[i] - ta1_y_xxzz_xxyy_1[i] * pc_z[i];

        ta1_y_xxzzz_xxyz_0[i] = ta1_y_zzz_xxyz_0[i] * fe_0 - ta1_y_zzz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_xzzz_xyz_1[i] * fe_0 + ta1_y_xzzz_xxyz_0[i] * pa_x[i] - ta1_y_xzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxzzz_xxzz_0[i] = ta1_y_zzz_xxzz_0[i] * fe_0 - ta1_y_zzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xzzz_xzz_0[i] * fe_0 -
                                2.0 * ta1_y_xzzz_xzz_1[i] * fe_0 + ta1_y_xzzz_xxzz_0[i] * pa_x[i] - ta1_y_xzzz_xxzz_1[i] * pc_x[i];

        ta1_y_xxzzz_xyyy_0[i] =
            2.0 * ta1_y_xxz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_xxz_xyyy_1[i] * fe_0 + ta1_y_xxzz_xyyy_0[i] * pa_z[i] - ta1_y_xxzz_xyyy_1[i] * pc_z[i];

        ta1_y_xxzzz_xyyz_0[i] = ta1_y_zzz_xyyz_0[i] * fe_0 - ta1_y_zzz_xyyz_1[i] * fe_0 + ta1_y_xzzz_yyz_0[i] * fe_0 - ta1_y_xzzz_yyz_1[i] * fe_0 +
                                ta1_y_xzzz_xyyz_0[i] * pa_x[i] - ta1_y_xzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxzzz_xyzz_0[i] = ta1_y_zzz_xyzz_0[i] * fe_0 - ta1_y_zzz_xyzz_1[i] * fe_0 + ta1_y_xzzz_yzz_0[i] * fe_0 - ta1_y_xzzz_yzz_1[i] * fe_0 +
                                ta1_y_xzzz_xyzz_0[i] * pa_x[i] - ta1_y_xzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxzzz_xzzz_0[i] = ta1_y_zzz_xzzz_0[i] * fe_0 - ta1_y_zzz_xzzz_1[i] * fe_0 + ta1_y_xzzz_zzz_0[i] * fe_0 - ta1_y_xzzz_zzz_1[i] * fe_0 +
                                ta1_y_xzzz_xzzz_0[i] * pa_x[i] - ta1_y_xzzz_xzzz_1[i] * pc_x[i];

        ta1_y_xxzzz_yyyy_0[i] =
            ta1_y_zzz_yyyy_0[i] * fe_0 - ta1_y_zzz_yyyy_1[i] * fe_0 + ta1_y_xzzz_yyyy_0[i] * pa_x[i] - ta1_y_xzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxzzz_yyyz_0[i] =
            ta1_y_zzz_yyyz_0[i] * fe_0 - ta1_y_zzz_yyyz_1[i] * fe_0 + ta1_y_xzzz_yyyz_0[i] * pa_x[i] - ta1_y_xzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxzzz_yyzz_0[i] =
            ta1_y_zzz_yyzz_0[i] * fe_0 - ta1_y_zzz_yyzz_1[i] * fe_0 + ta1_y_xzzz_yyzz_0[i] * pa_x[i] - ta1_y_xzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxzzz_yzzz_0[i] =
            ta1_y_zzz_yzzz_0[i] * fe_0 - ta1_y_zzz_yzzz_1[i] * fe_0 + ta1_y_xzzz_yzzz_0[i] * pa_x[i] - ta1_y_xzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxzzz_zzzz_0[i] =
            ta1_y_zzz_zzzz_0[i] * fe_0 - ta1_y_zzz_zzzz_1[i] * fe_0 + ta1_y_xzzz_zzzz_0[i] * pa_x[i] - ta1_y_xzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 465-480 components of targeted buffer : HG

    auto ta1_y_xyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 465);

    auto ta1_y_xyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 466);

    auto ta1_y_xyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 467);

    auto ta1_y_xyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 468);

    auto ta1_y_xyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 469);

    auto ta1_y_xyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 470);

    auto ta1_y_xyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 471);

    auto ta1_y_xyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 472);

    auto ta1_y_xyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 473);

    auto ta1_y_xyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 474);

    auto ta1_y_xyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 475);

    auto ta1_y_xyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 476);

    auto ta1_y_xyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 477);

    auto ta1_y_xyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 478);

    auto ta1_y_xyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 479);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyyyy_xxxx_0, \
                             ta1_y_xyyyy_xxxy_0, \
                             ta1_y_xyyyy_xxxz_0, \
                             ta1_y_xyyyy_xxyy_0, \
                             ta1_y_xyyyy_xxyz_0, \
                             ta1_y_xyyyy_xxzz_0, \
                             ta1_y_xyyyy_xyyy_0, \
                             ta1_y_xyyyy_xyyz_0, \
                             ta1_y_xyyyy_xyzz_0, \
                             ta1_y_xyyyy_xzzz_0, \
                             ta1_y_xyyyy_yyyy_0, \
                             ta1_y_xyyyy_yyyz_0, \
                             ta1_y_xyyyy_yyzz_0, \
                             ta1_y_xyyyy_yzzz_0, \
                             ta1_y_xyyyy_zzzz_0, \
                             ta1_y_yyyy_xxx_0,   \
                             ta1_y_yyyy_xxx_1,   \
                             ta1_y_yyyy_xxxx_0,  \
                             ta1_y_yyyy_xxxx_1,  \
                             ta1_y_yyyy_xxxy_0,  \
                             ta1_y_yyyy_xxxy_1,  \
                             ta1_y_yyyy_xxxz_0,  \
                             ta1_y_yyyy_xxxz_1,  \
                             ta1_y_yyyy_xxy_0,   \
                             ta1_y_yyyy_xxy_1,   \
                             ta1_y_yyyy_xxyy_0,  \
                             ta1_y_yyyy_xxyy_1,  \
                             ta1_y_yyyy_xxyz_0,  \
                             ta1_y_yyyy_xxyz_1,  \
                             ta1_y_yyyy_xxz_0,   \
                             ta1_y_yyyy_xxz_1,   \
                             ta1_y_yyyy_xxzz_0,  \
                             ta1_y_yyyy_xxzz_1,  \
                             ta1_y_yyyy_xyy_0,   \
                             ta1_y_yyyy_xyy_1,   \
                             ta1_y_yyyy_xyyy_0,  \
                             ta1_y_yyyy_xyyy_1,  \
                             ta1_y_yyyy_xyyz_0,  \
                             ta1_y_yyyy_xyyz_1,  \
                             ta1_y_yyyy_xyz_0,   \
                             ta1_y_yyyy_xyz_1,   \
                             ta1_y_yyyy_xyzz_0,  \
                             ta1_y_yyyy_xyzz_1,  \
                             ta1_y_yyyy_xzz_0,   \
                             ta1_y_yyyy_xzz_1,   \
                             ta1_y_yyyy_xzzz_0,  \
                             ta1_y_yyyy_xzzz_1,  \
                             ta1_y_yyyy_yyy_0,   \
                             ta1_y_yyyy_yyy_1,   \
                             ta1_y_yyyy_yyyy_0,  \
                             ta1_y_yyyy_yyyy_1,  \
                             ta1_y_yyyy_yyyz_0,  \
                             ta1_y_yyyy_yyyz_1,  \
                             ta1_y_yyyy_yyz_0,   \
                             ta1_y_yyyy_yyz_1,   \
                             ta1_y_yyyy_yyzz_0,  \
                             ta1_y_yyyy_yyzz_1,  \
                             ta1_y_yyyy_yzz_0,   \
                             ta1_y_yyyy_yzz_1,   \
                             ta1_y_yyyy_yzzz_0,  \
                             ta1_y_yyyy_yzzz_1,  \
                             ta1_y_yyyy_zzz_0,   \
                             ta1_y_yyyy_zzz_1,   \
                             ta1_y_yyyy_zzzz_0,  \
                             ta1_y_yyyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyy_xxxx_0[i] =
            4.0 * ta1_y_yyyy_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyyy_xxx_1[i] * fe_0 + ta1_y_yyyy_xxxx_0[i] * pa_x[i] - ta1_y_yyyy_xxxx_1[i] * pc_x[i];

        ta1_y_xyyyy_xxxy_0[i] =
            3.0 * ta1_y_yyyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyyy_xxy_1[i] * fe_0 + ta1_y_yyyy_xxxy_0[i] * pa_x[i] - ta1_y_yyyy_xxxy_1[i] * pc_x[i];

        ta1_y_xyyyy_xxxz_0[i] =
            3.0 * ta1_y_yyyy_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyyy_xxz_1[i] * fe_0 + ta1_y_yyyy_xxxz_0[i] * pa_x[i] - ta1_y_yyyy_xxxz_1[i] * pc_x[i];

        ta1_y_xyyyy_xxyy_0[i] =
            2.0 * ta1_y_yyyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xyy_1[i] * fe_0 + ta1_y_yyyy_xxyy_0[i] * pa_x[i] - ta1_y_yyyy_xxyy_1[i] * pc_x[i];

        ta1_y_xyyyy_xxyz_0[i] =
            2.0 * ta1_y_yyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xyz_1[i] * fe_0 + ta1_y_yyyy_xxyz_0[i] * pa_x[i] - ta1_y_yyyy_xxyz_1[i] * pc_x[i];

        ta1_y_xyyyy_xxzz_0[i] =
            2.0 * ta1_y_yyyy_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xzz_1[i] * fe_0 + ta1_y_yyyy_xxzz_0[i] * pa_x[i] - ta1_y_yyyy_xxzz_1[i] * pc_x[i];

        ta1_y_xyyyy_xyyy_0[i] =
            ta1_y_yyyy_yyy_0[i] * fe_0 - ta1_y_yyyy_yyy_1[i] * fe_0 + ta1_y_yyyy_xyyy_0[i] * pa_x[i] - ta1_y_yyyy_xyyy_1[i] * pc_x[i];

        ta1_y_xyyyy_xyyz_0[i] =
            ta1_y_yyyy_yyz_0[i] * fe_0 - ta1_y_yyyy_yyz_1[i] * fe_0 + ta1_y_yyyy_xyyz_0[i] * pa_x[i] - ta1_y_yyyy_xyyz_1[i] * pc_x[i];

        ta1_y_xyyyy_xyzz_0[i] =
            ta1_y_yyyy_yzz_0[i] * fe_0 - ta1_y_yyyy_yzz_1[i] * fe_0 + ta1_y_yyyy_xyzz_0[i] * pa_x[i] - ta1_y_yyyy_xyzz_1[i] * pc_x[i];

        ta1_y_xyyyy_xzzz_0[i] =
            ta1_y_yyyy_zzz_0[i] * fe_0 - ta1_y_yyyy_zzz_1[i] * fe_0 + ta1_y_yyyy_xzzz_0[i] * pa_x[i] - ta1_y_yyyy_xzzz_1[i] * pc_x[i];

        ta1_y_xyyyy_yyyy_0[i] = ta1_y_yyyy_yyyy_0[i] * pa_x[i] - ta1_y_yyyy_yyyy_1[i] * pc_x[i];

        ta1_y_xyyyy_yyyz_0[i] = ta1_y_yyyy_yyyz_0[i] * pa_x[i] - ta1_y_yyyy_yyyz_1[i] * pc_x[i];

        ta1_y_xyyyy_yyzz_0[i] = ta1_y_yyyy_yyzz_0[i] * pa_x[i] - ta1_y_yyyy_yyzz_1[i] * pc_x[i];

        ta1_y_xyyyy_yzzz_0[i] = ta1_y_yyyy_yzzz_0[i] * pa_x[i] - ta1_y_yyyy_yzzz_1[i] * pc_x[i];

        ta1_y_xyyyy_zzzz_0[i] = ta1_y_yyyy_zzzz_0[i] * pa_x[i] - ta1_y_yyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 480-495 components of targeted buffer : HG

    auto ta1_y_xyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 480);

    auto ta1_y_xyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 481);

    auto ta1_y_xyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 482);

    auto ta1_y_xyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 483);

    auto ta1_y_xyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 484);

    auto ta1_y_xyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 485);

    auto ta1_y_xyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 486);

    auto ta1_y_xyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 487);

    auto ta1_y_xyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 488);

    auto ta1_y_xyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 489);

    auto ta1_y_xyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 490);

    auto ta1_y_xyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 491);

    auto ta1_y_xyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 492);

    auto ta1_y_xyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 493);

    auto ta1_y_xyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 494);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xyyy_xxxx_0,  \
                             ta1_y_xyyy_xxxx_1,  \
                             ta1_y_xyyy_xxxy_0,  \
                             ta1_y_xyyy_xxxy_1,  \
                             ta1_y_xyyy_xxyy_0,  \
                             ta1_y_xyyy_xxyy_1,  \
                             ta1_y_xyyy_xyyy_0,  \
                             ta1_y_xyyy_xyyy_1,  \
                             ta1_y_xyyyz_xxxx_0, \
                             ta1_y_xyyyz_xxxy_0, \
                             ta1_y_xyyyz_xxxz_0, \
                             ta1_y_xyyyz_xxyy_0, \
                             ta1_y_xyyyz_xxyz_0, \
                             ta1_y_xyyyz_xxzz_0, \
                             ta1_y_xyyyz_xyyy_0, \
                             ta1_y_xyyyz_xyyz_0, \
                             ta1_y_xyyyz_xyzz_0, \
                             ta1_y_xyyyz_xzzz_0, \
                             ta1_y_xyyyz_yyyy_0, \
                             ta1_y_xyyyz_yyyz_0, \
                             ta1_y_xyyyz_yyzz_0, \
                             ta1_y_xyyyz_yzzz_0, \
                             ta1_y_xyyyz_zzzz_0, \
                             ta1_y_yyyz_xxxz_0,  \
                             ta1_y_yyyz_xxxz_1,  \
                             ta1_y_yyyz_xxyz_0,  \
                             ta1_y_yyyz_xxyz_1,  \
                             ta1_y_yyyz_xxz_0,   \
                             ta1_y_yyyz_xxz_1,   \
                             ta1_y_yyyz_xxzz_0,  \
                             ta1_y_yyyz_xxzz_1,  \
                             ta1_y_yyyz_xyyz_0,  \
                             ta1_y_yyyz_xyyz_1,  \
                             ta1_y_yyyz_xyz_0,   \
                             ta1_y_yyyz_xyz_1,   \
                             ta1_y_yyyz_xyzz_0,  \
                             ta1_y_yyyz_xyzz_1,  \
                             ta1_y_yyyz_xzz_0,   \
                             ta1_y_yyyz_xzz_1,   \
                             ta1_y_yyyz_xzzz_0,  \
                             ta1_y_yyyz_xzzz_1,  \
                             ta1_y_yyyz_yyyy_0,  \
                             ta1_y_yyyz_yyyy_1,  \
                             ta1_y_yyyz_yyyz_0,  \
                             ta1_y_yyyz_yyyz_1,  \
                             ta1_y_yyyz_yyz_0,   \
                             ta1_y_yyyz_yyz_1,   \
                             ta1_y_yyyz_yyzz_0,  \
                             ta1_y_yyyz_yyzz_1,  \
                             ta1_y_yyyz_yzz_0,   \
                             ta1_y_yyyz_yzz_1,   \
                             ta1_y_yyyz_yzzz_0,  \
                             ta1_y_yyyz_yzzz_1,  \
                             ta1_y_yyyz_zzz_0,   \
                             ta1_y_yyyz_zzz_1,   \
                             ta1_y_yyyz_zzzz_0,  \
                             ta1_y_yyyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyz_xxxx_0[i] = ta1_y_xyyy_xxxx_0[i] * pa_z[i] - ta1_y_xyyy_xxxx_1[i] * pc_z[i];

        ta1_y_xyyyz_xxxy_0[i] = ta1_y_xyyy_xxxy_0[i] * pa_z[i] - ta1_y_xyyy_xxxy_1[i] * pc_z[i];

        ta1_y_xyyyz_xxxz_0[i] =
            3.0 * ta1_y_yyyz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyyz_xxz_1[i] * fe_0 + ta1_y_yyyz_xxxz_0[i] * pa_x[i] - ta1_y_yyyz_xxxz_1[i] * pc_x[i];

        ta1_y_xyyyz_xxyy_0[i] = ta1_y_xyyy_xxyy_0[i] * pa_z[i] - ta1_y_xyyy_xxyy_1[i] * pc_z[i];

        ta1_y_xyyyz_xxyz_0[i] =
            2.0 * ta1_y_yyyz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xyz_1[i] * fe_0 + ta1_y_yyyz_xxyz_0[i] * pa_x[i] - ta1_y_yyyz_xxyz_1[i] * pc_x[i];

        ta1_y_xyyyz_xxzz_0[i] =
            2.0 * ta1_y_yyyz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xzz_1[i] * fe_0 + ta1_y_yyyz_xxzz_0[i] * pa_x[i] - ta1_y_yyyz_xxzz_1[i] * pc_x[i];

        ta1_y_xyyyz_xyyy_0[i] = ta1_y_xyyy_xyyy_0[i] * pa_z[i] - ta1_y_xyyy_xyyy_1[i] * pc_z[i];

        ta1_y_xyyyz_xyyz_0[i] =
            ta1_y_yyyz_yyz_0[i] * fe_0 - ta1_y_yyyz_yyz_1[i] * fe_0 + ta1_y_yyyz_xyyz_0[i] * pa_x[i] - ta1_y_yyyz_xyyz_1[i] * pc_x[i];

        ta1_y_xyyyz_xyzz_0[i] =
            ta1_y_yyyz_yzz_0[i] * fe_0 - ta1_y_yyyz_yzz_1[i] * fe_0 + ta1_y_yyyz_xyzz_0[i] * pa_x[i] - ta1_y_yyyz_xyzz_1[i] * pc_x[i];

        ta1_y_xyyyz_xzzz_0[i] =
            ta1_y_yyyz_zzz_0[i] * fe_0 - ta1_y_yyyz_zzz_1[i] * fe_0 + ta1_y_yyyz_xzzz_0[i] * pa_x[i] - ta1_y_yyyz_xzzz_1[i] * pc_x[i];

        ta1_y_xyyyz_yyyy_0[i] = ta1_y_yyyz_yyyy_0[i] * pa_x[i] - ta1_y_yyyz_yyyy_1[i] * pc_x[i];

        ta1_y_xyyyz_yyyz_0[i] = ta1_y_yyyz_yyyz_0[i] * pa_x[i] - ta1_y_yyyz_yyyz_1[i] * pc_x[i];

        ta1_y_xyyyz_yyzz_0[i] = ta1_y_yyyz_yyzz_0[i] * pa_x[i] - ta1_y_yyyz_yyzz_1[i] * pc_x[i];

        ta1_y_xyyyz_yzzz_0[i] = ta1_y_yyyz_yzzz_0[i] * pa_x[i] - ta1_y_yyyz_yzzz_1[i] * pc_x[i];

        ta1_y_xyyyz_zzzz_0[i] = ta1_y_yyyz_zzzz_0[i] * pa_x[i] - ta1_y_yyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 495-510 components of targeted buffer : HG

    auto ta1_y_xyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 495);

    auto ta1_y_xyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 496);

    auto ta1_y_xyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 497);

    auto ta1_y_xyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 498);

    auto ta1_y_xyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 499);

    auto ta1_y_xyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 500);

    auto ta1_y_xyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 501);

    auto ta1_y_xyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 502);

    auto ta1_y_xyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 503);

    auto ta1_y_xyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 504);

    auto ta1_y_xyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 505);

    auto ta1_y_xyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 506);

    auto ta1_y_xyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 507);

    auto ta1_y_xyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 508);

    auto ta1_y_xyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 509);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyyzz_xxxx_0, \
                             ta1_y_xyyzz_xxxy_0, \
                             ta1_y_xyyzz_xxxz_0, \
                             ta1_y_xyyzz_xxyy_0, \
                             ta1_y_xyyzz_xxyz_0, \
                             ta1_y_xyyzz_xxzz_0, \
                             ta1_y_xyyzz_xyyy_0, \
                             ta1_y_xyyzz_xyyz_0, \
                             ta1_y_xyyzz_xyzz_0, \
                             ta1_y_xyyzz_xzzz_0, \
                             ta1_y_xyyzz_yyyy_0, \
                             ta1_y_xyyzz_yyyz_0, \
                             ta1_y_xyyzz_yyzz_0, \
                             ta1_y_xyyzz_yzzz_0, \
                             ta1_y_xyyzz_zzzz_0, \
                             ta1_y_yyzz_xxx_0,   \
                             ta1_y_yyzz_xxx_1,   \
                             ta1_y_yyzz_xxxx_0,  \
                             ta1_y_yyzz_xxxx_1,  \
                             ta1_y_yyzz_xxxy_0,  \
                             ta1_y_yyzz_xxxy_1,  \
                             ta1_y_yyzz_xxxz_0,  \
                             ta1_y_yyzz_xxxz_1,  \
                             ta1_y_yyzz_xxy_0,   \
                             ta1_y_yyzz_xxy_1,   \
                             ta1_y_yyzz_xxyy_0,  \
                             ta1_y_yyzz_xxyy_1,  \
                             ta1_y_yyzz_xxyz_0,  \
                             ta1_y_yyzz_xxyz_1,  \
                             ta1_y_yyzz_xxz_0,   \
                             ta1_y_yyzz_xxz_1,   \
                             ta1_y_yyzz_xxzz_0,  \
                             ta1_y_yyzz_xxzz_1,  \
                             ta1_y_yyzz_xyy_0,   \
                             ta1_y_yyzz_xyy_1,   \
                             ta1_y_yyzz_xyyy_0,  \
                             ta1_y_yyzz_xyyy_1,  \
                             ta1_y_yyzz_xyyz_0,  \
                             ta1_y_yyzz_xyyz_1,  \
                             ta1_y_yyzz_xyz_0,   \
                             ta1_y_yyzz_xyz_1,   \
                             ta1_y_yyzz_xyzz_0,  \
                             ta1_y_yyzz_xyzz_1,  \
                             ta1_y_yyzz_xzz_0,   \
                             ta1_y_yyzz_xzz_1,   \
                             ta1_y_yyzz_xzzz_0,  \
                             ta1_y_yyzz_xzzz_1,  \
                             ta1_y_yyzz_yyy_0,   \
                             ta1_y_yyzz_yyy_1,   \
                             ta1_y_yyzz_yyyy_0,  \
                             ta1_y_yyzz_yyyy_1,  \
                             ta1_y_yyzz_yyyz_0,  \
                             ta1_y_yyzz_yyyz_1,  \
                             ta1_y_yyzz_yyz_0,   \
                             ta1_y_yyzz_yyz_1,   \
                             ta1_y_yyzz_yyzz_0,  \
                             ta1_y_yyzz_yyzz_1,  \
                             ta1_y_yyzz_yzz_0,   \
                             ta1_y_yyzz_yzz_1,   \
                             ta1_y_yyzz_yzzz_0,  \
                             ta1_y_yyzz_yzzz_1,  \
                             ta1_y_yyzz_zzz_0,   \
                             ta1_y_yyzz_zzz_1,   \
                             ta1_y_yyzz_zzzz_0,  \
                             ta1_y_yyzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzz_xxxx_0[i] =
            4.0 * ta1_y_yyzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyzz_xxx_1[i] * fe_0 + ta1_y_yyzz_xxxx_0[i] * pa_x[i] - ta1_y_yyzz_xxxx_1[i] * pc_x[i];

        ta1_y_xyyzz_xxxy_0[i] =
            3.0 * ta1_y_yyzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxy_1[i] * fe_0 + ta1_y_yyzz_xxxy_0[i] * pa_x[i] - ta1_y_yyzz_xxxy_1[i] * pc_x[i];

        ta1_y_xyyzz_xxxz_0[i] =
            3.0 * ta1_y_yyzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxz_1[i] * fe_0 + ta1_y_yyzz_xxxz_0[i] * pa_x[i] - ta1_y_yyzz_xxxz_1[i] * pc_x[i];

        ta1_y_xyyzz_xxyy_0[i] =
            2.0 * ta1_y_yyzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyzz_xyy_1[i] * fe_0 + ta1_y_yyzz_xxyy_0[i] * pa_x[i] - ta1_y_yyzz_xxyy_1[i] * pc_x[i];

        ta1_y_xyyzz_xxyz_0[i] =
            2.0 * ta1_y_yyzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyzz_xyz_1[i] * fe_0 + ta1_y_yyzz_xxyz_0[i] * pa_x[i] - ta1_y_yyzz_xxyz_1[i] * pc_x[i];

        ta1_y_xyyzz_xxzz_0[i] =
            2.0 * ta1_y_yyzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyzz_xzz_1[i] * fe_0 + ta1_y_yyzz_xxzz_0[i] * pa_x[i] - ta1_y_yyzz_xxzz_1[i] * pc_x[i];

        ta1_y_xyyzz_xyyy_0[i] =
            ta1_y_yyzz_yyy_0[i] * fe_0 - ta1_y_yyzz_yyy_1[i] * fe_0 + ta1_y_yyzz_xyyy_0[i] * pa_x[i] - ta1_y_yyzz_xyyy_1[i] * pc_x[i];

        ta1_y_xyyzz_xyyz_0[i] =
            ta1_y_yyzz_yyz_0[i] * fe_0 - ta1_y_yyzz_yyz_1[i] * fe_0 + ta1_y_yyzz_xyyz_0[i] * pa_x[i] - ta1_y_yyzz_xyyz_1[i] * pc_x[i];

        ta1_y_xyyzz_xyzz_0[i] =
            ta1_y_yyzz_yzz_0[i] * fe_0 - ta1_y_yyzz_yzz_1[i] * fe_0 + ta1_y_yyzz_xyzz_0[i] * pa_x[i] - ta1_y_yyzz_xyzz_1[i] * pc_x[i];

        ta1_y_xyyzz_xzzz_0[i] =
            ta1_y_yyzz_zzz_0[i] * fe_0 - ta1_y_yyzz_zzz_1[i] * fe_0 + ta1_y_yyzz_xzzz_0[i] * pa_x[i] - ta1_y_yyzz_xzzz_1[i] * pc_x[i];

        ta1_y_xyyzz_yyyy_0[i] = ta1_y_yyzz_yyyy_0[i] * pa_x[i] - ta1_y_yyzz_yyyy_1[i] * pc_x[i];

        ta1_y_xyyzz_yyyz_0[i] = ta1_y_yyzz_yyyz_0[i] * pa_x[i] - ta1_y_yyzz_yyyz_1[i] * pc_x[i];

        ta1_y_xyyzz_yyzz_0[i] = ta1_y_yyzz_yyzz_0[i] * pa_x[i] - ta1_y_yyzz_yyzz_1[i] * pc_x[i];

        ta1_y_xyyzz_yzzz_0[i] = ta1_y_yyzz_yzzz_0[i] * pa_x[i] - ta1_y_yyzz_yzzz_1[i] * pc_x[i];

        ta1_y_xyyzz_zzzz_0[i] = ta1_y_yyzz_zzzz_0[i] * pa_x[i] - ta1_y_yyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 510-525 components of targeted buffer : HG

    auto ta1_y_xyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 510);

    auto ta1_y_xyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 511);

    auto ta1_y_xyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 512);

    auto ta1_y_xyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 513);

    auto ta1_y_xyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 514);

    auto ta1_y_xyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 515);

    auto ta1_y_xyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 516);

    auto ta1_y_xyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 517);

    auto ta1_y_xyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 518);

    auto ta1_y_xyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 519);

    auto ta1_y_xyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 520);

    auto ta1_y_xyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 521);

    auto ta1_y_xyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 522);

    auto ta1_y_xyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 523);

    auto ta1_y_xyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 524);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xyzzz_xxxx_0, \
                             ta1_y_xyzzz_xxxy_0, \
                             ta1_y_xyzzz_xxxz_0, \
                             ta1_y_xyzzz_xxyy_0, \
                             ta1_y_xyzzz_xxyz_0, \
                             ta1_y_xyzzz_xxzz_0, \
                             ta1_y_xyzzz_xyyy_0, \
                             ta1_y_xyzzz_xyyz_0, \
                             ta1_y_xyzzz_xyzz_0, \
                             ta1_y_xyzzz_xzzz_0, \
                             ta1_y_xyzzz_yyyy_0, \
                             ta1_y_xyzzz_yyyz_0, \
                             ta1_y_xyzzz_yyzz_0, \
                             ta1_y_xyzzz_yzzz_0, \
                             ta1_y_xyzzz_zzzz_0, \
                             ta1_y_xzzz_xxxx_0,  \
                             ta1_y_xzzz_xxxx_1,  \
                             ta1_y_xzzz_xxxz_0,  \
                             ta1_y_xzzz_xxxz_1,  \
                             ta1_y_xzzz_xxzz_0,  \
                             ta1_y_xzzz_xxzz_1,  \
                             ta1_y_xzzz_xzzz_0,  \
                             ta1_y_xzzz_xzzz_1,  \
                             ta1_y_yzzz_xxxy_0,  \
                             ta1_y_yzzz_xxxy_1,  \
                             ta1_y_yzzz_xxy_0,   \
                             ta1_y_yzzz_xxy_1,   \
                             ta1_y_yzzz_xxyy_0,  \
                             ta1_y_yzzz_xxyy_1,  \
                             ta1_y_yzzz_xxyz_0,  \
                             ta1_y_yzzz_xxyz_1,  \
                             ta1_y_yzzz_xyy_0,   \
                             ta1_y_yzzz_xyy_1,   \
                             ta1_y_yzzz_xyyy_0,  \
                             ta1_y_yzzz_xyyy_1,  \
                             ta1_y_yzzz_xyyz_0,  \
                             ta1_y_yzzz_xyyz_1,  \
                             ta1_y_yzzz_xyz_0,   \
                             ta1_y_yzzz_xyz_1,   \
                             ta1_y_yzzz_xyzz_0,  \
                             ta1_y_yzzz_xyzz_1,  \
                             ta1_y_yzzz_yyy_0,   \
                             ta1_y_yzzz_yyy_1,   \
                             ta1_y_yzzz_yyyy_0,  \
                             ta1_y_yzzz_yyyy_1,  \
                             ta1_y_yzzz_yyyz_0,  \
                             ta1_y_yzzz_yyyz_1,  \
                             ta1_y_yzzz_yyz_0,   \
                             ta1_y_yzzz_yyz_1,   \
                             ta1_y_yzzz_yyzz_0,  \
                             ta1_y_yzzz_yyzz_1,  \
                             ta1_y_yzzz_yzz_0,   \
                             ta1_y_yzzz_yzz_1,   \
                             ta1_y_yzzz_yzzz_0,  \
                             ta1_y_yzzz_yzzz_1,  \
                             ta1_y_yzzz_zzzz_0,  \
                             ta1_y_yzzz_zzzz_1,  \
                             ta_xzzz_xxxx_1,     \
                             ta_xzzz_xxxz_1,     \
                             ta_xzzz_xxzz_1,     \
                             ta_xzzz_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzzz_xxxx_0[i] = ta_xzzz_xxxx_1[i] + ta1_y_xzzz_xxxx_0[i] * pa_y[i] - ta1_y_xzzz_xxxx_1[i] * pc_y[i];

        ta1_y_xyzzz_xxxy_0[i] =
            3.0 * ta1_y_yzzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yzzz_xxy_1[i] * fe_0 + ta1_y_yzzz_xxxy_0[i] * pa_x[i] - ta1_y_yzzz_xxxy_1[i] * pc_x[i];

        ta1_y_xyzzz_xxxz_0[i] = ta_xzzz_xxxz_1[i] + ta1_y_xzzz_xxxz_0[i] * pa_y[i] - ta1_y_xzzz_xxxz_1[i] * pc_y[i];

        ta1_y_xyzzz_xxyy_0[i] =
            2.0 * ta1_y_yzzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xyy_1[i] * fe_0 + ta1_y_yzzz_xxyy_0[i] * pa_x[i] - ta1_y_yzzz_xxyy_1[i] * pc_x[i];

        ta1_y_xyzzz_xxyz_0[i] =
            2.0 * ta1_y_yzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xyz_1[i] * fe_0 + ta1_y_yzzz_xxyz_0[i] * pa_x[i] - ta1_y_yzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xyzzz_xxzz_0[i] = ta_xzzz_xxzz_1[i] + ta1_y_xzzz_xxzz_0[i] * pa_y[i] - ta1_y_xzzz_xxzz_1[i] * pc_y[i];

        ta1_y_xyzzz_xyyy_0[i] =
            ta1_y_yzzz_yyy_0[i] * fe_0 - ta1_y_yzzz_yyy_1[i] * fe_0 + ta1_y_yzzz_xyyy_0[i] * pa_x[i] - ta1_y_yzzz_xyyy_1[i] * pc_x[i];

        ta1_y_xyzzz_xyyz_0[i] =
            ta1_y_yzzz_yyz_0[i] * fe_0 - ta1_y_yzzz_yyz_1[i] * fe_0 + ta1_y_yzzz_xyyz_0[i] * pa_x[i] - ta1_y_yzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xyzzz_xyzz_0[i] =
            ta1_y_yzzz_yzz_0[i] * fe_0 - ta1_y_yzzz_yzz_1[i] * fe_0 + ta1_y_yzzz_xyzz_0[i] * pa_x[i] - ta1_y_yzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xyzzz_xzzz_0[i] = ta_xzzz_xzzz_1[i] + ta1_y_xzzz_xzzz_0[i] * pa_y[i] - ta1_y_xzzz_xzzz_1[i] * pc_y[i];

        ta1_y_xyzzz_yyyy_0[i] = ta1_y_yzzz_yyyy_0[i] * pa_x[i] - ta1_y_yzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xyzzz_yyyz_0[i] = ta1_y_yzzz_yyyz_0[i] * pa_x[i] - ta1_y_yzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xyzzz_yyzz_0[i] = ta1_y_yzzz_yyzz_0[i] * pa_x[i] - ta1_y_yzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xyzzz_yzzz_0[i] = ta1_y_yzzz_yzzz_0[i] * pa_x[i] - ta1_y_yzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xyzzz_zzzz_0[i] = ta1_y_yzzz_zzzz_0[i] * pa_x[i] - ta1_y_yzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 525-540 components of targeted buffer : HG

    auto ta1_y_xzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 525);

    auto ta1_y_xzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 526);

    auto ta1_y_xzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 527);

    auto ta1_y_xzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 528);

    auto ta1_y_xzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 529);

    auto ta1_y_xzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 530);

    auto ta1_y_xzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 531);

    auto ta1_y_xzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 532);

    auto ta1_y_xzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 533);

    auto ta1_y_xzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 534);

    auto ta1_y_xzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 535);

    auto ta1_y_xzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 536);

    auto ta1_y_xzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 537);

    auto ta1_y_xzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 538);

    auto ta1_y_xzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 539);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xzzzz_xxxx_0, \
                             ta1_y_xzzzz_xxxy_0, \
                             ta1_y_xzzzz_xxxz_0, \
                             ta1_y_xzzzz_xxyy_0, \
                             ta1_y_xzzzz_xxyz_0, \
                             ta1_y_xzzzz_xxzz_0, \
                             ta1_y_xzzzz_xyyy_0, \
                             ta1_y_xzzzz_xyyz_0, \
                             ta1_y_xzzzz_xyzz_0, \
                             ta1_y_xzzzz_xzzz_0, \
                             ta1_y_xzzzz_yyyy_0, \
                             ta1_y_xzzzz_yyyz_0, \
                             ta1_y_xzzzz_yyzz_0, \
                             ta1_y_xzzzz_yzzz_0, \
                             ta1_y_xzzzz_zzzz_0, \
                             ta1_y_zzzz_xxx_0,   \
                             ta1_y_zzzz_xxx_1,   \
                             ta1_y_zzzz_xxxx_0,  \
                             ta1_y_zzzz_xxxx_1,  \
                             ta1_y_zzzz_xxxy_0,  \
                             ta1_y_zzzz_xxxy_1,  \
                             ta1_y_zzzz_xxxz_0,  \
                             ta1_y_zzzz_xxxz_1,  \
                             ta1_y_zzzz_xxy_0,   \
                             ta1_y_zzzz_xxy_1,   \
                             ta1_y_zzzz_xxyy_0,  \
                             ta1_y_zzzz_xxyy_1,  \
                             ta1_y_zzzz_xxyz_0,  \
                             ta1_y_zzzz_xxyz_1,  \
                             ta1_y_zzzz_xxz_0,   \
                             ta1_y_zzzz_xxz_1,   \
                             ta1_y_zzzz_xxzz_0,  \
                             ta1_y_zzzz_xxzz_1,  \
                             ta1_y_zzzz_xyy_0,   \
                             ta1_y_zzzz_xyy_1,   \
                             ta1_y_zzzz_xyyy_0,  \
                             ta1_y_zzzz_xyyy_1,  \
                             ta1_y_zzzz_xyyz_0,  \
                             ta1_y_zzzz_xyyz_1,  \
                             ta1_y_zzzz_xyz_0,   \
                             ta1_y_zzzz_xyz_1,   \
                             ta1_y_zzzz_xyzz_0,  \
                             ta1_y_zzzz_xyzz_1,  \
                             ta1_y_zzzz_xzz_0,   \
                             ta1_y_zzzz_xzz_1,   \
                             ta1_y_zzzz_xzzz_0,  \
                             ta1_y_zzzz_xzzz_1,  \
                             ta1_y_zzzz_yyy_0,   \
                             ta1_y_zzzz_yyy_1,   \
                             ta1_y_zzzz_yyyy_0,  \
                             ta1_y_zzzz_yyyy_1,  \
                             ta1_y_zzzz_yyyz_0,  \
                             ta1_y_zzzz_yyyz_1,  \
                             ta1_y_zzzz_yyz_0,   \
                             ta1_y_zzzz_yyz_1,   \
                             ta1_y_zzzz_yyzz_0,  \
                             ta1_y_zzzz_yyzz_1,  \
                             ta1_y_zzzz_yzz_0,   \
                             ta1_y_zzzz_yzz_1,   \
                             ta1_y_zzzz_yzzz_0,  \
                             ta1_y_zzzz_yzzz_1,  \
                             ta1_y_zzzz_zzz_0,   \
                             ta1_y_zzzz_zzz_1,   \
                             ta1_y_zzzz_zzzz_0,  \
                             ta1_y_zzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzz_xxxx_0[i] =
            4.0 * ta1_y_zzzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_zzzz_xxx_1[i] * fe_0 + ta1_y_zzzz_xxxx_0[i] * pa_x[i] - ta1_y_zzzz_xxxx_1[i] * pc_x[i];

        ta1_y_xzzzz_xxxy_0[i] =
            3.0 * ta1_y_zzzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_zzzz_xxy_1[i] * fe_0 + ta1_y_zzzz_xxxy_0[i] * pa_x[i] - ta1_y_zzzz_xxxy_1[i] * pc_x[i];

        ta1_y_xzzzz_xxxz_0[i] =
            3.0 * ta1_y_zzzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_zzzz_xxz_1[i] * fe_0 + ta1_y_zzzz_xxxz_0[i] * pa_x[i] - ta1_y_zzzz_xxxz_1[i] * pc_x[i];

        ta1_y_xzzzz_xxyy_0[i] =
            2.0 * ta1_y_zzzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xyy_1[i] * fe_0 + ta1_y_zzzz_xxyy_0[i] * pa_x[i] - ta1_y_zzzz_xxyy_1[i] * pc_x[i];

        ta1_y_xzzzz_xxyz_0[i] =
            2.0 * ta1_y_zzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xyz_1[i] * fe_0 + ta1_y_zzzz_xxyz_0[i] * pa_x[i] - ta1_y_zzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xzzzz_xxzz_0[i] =
            2.0 * ta1_y_zzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xzz_1[i] * fe_0 + ta1_y_zzzz_xxzz_0[i] * pa_x[i] - ta1_y_zzzz_xxzz_1[i] * pc_x[i];

        ta1_y_xzzzz_xyyy_0[i] =
            ta1_y_zzzz_yyy_0[i] * fe_0 - ta1_y_zzzz_yyy_1[i] * fe_0 + ta1_y_zzzz_xyyy_0[i] * pa_x[i] - ta1_y_zzzz_xyyy_1[i] * pc_x[i];

        ta1_y_xzzzz_xyyz_0[i] =
            ta1_y_zzzz_yyz_0[i] * fe_0 - ta1_y_zzzz_yyz_1[i] * fe_0 + ta1_y_zzzz_xyyz_0[i] * pa_x[i] - ta1_y_zzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xzzzz_xyzz_0[i] =
            ta1_y_zzzz_yzz_0[i] * fe_0 - ta1_y_zzzz_yzz_1[i] * fe_0 + ta1_y_zzzz_xyzz_0[i] * pa_x[i] - ta1_y_zzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xzzzz_xzzz_0[i] =
            ta1_y_zzzz_zzz_0[i] * fe_0 - ta1_y_zzzz_zzz_1[i] * fe_0 + ta1_y_zzzz_xzzz_0[i] * pa_x[i] - ta1_y_zzzz_xzzz_1[i] * pc_x[i];

        ta1_y_xzzzz_yyyy_0[i] = ta1_y_zzzz_yyyy_0[i] * pa_x[i] - ta1_y_zzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xzzzz_yyyz_0[i] = ta1_y_zzzz_yyyz_0[i] * pa_x[i] - ta1_y_zzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xzzzz_yyzz_0[i] = ta1_y_zzzz_yyzz_0[i] * pa_x[i] - ta1_y_zzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xzzzz_yzzz_0[i] = ta1_y_zzzz_yzzz_0[i] * pa_x[i] - ta1_y_zzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xzzzz_zzzz_0[i] = ta1_y_zzzz_zzzz_0[i] * pa_x[i] - ta1_y_zzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 540-555 components of targeted buffer : HG

    auto ta1_y_yyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 540);

    auto ta1_y_yyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 541);

    auto ta1_y_yyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 542);

    auto ta1_y_yyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 543);

    auto ta1_y_yyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 544);

    auto ta1_y_yyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 545);

    auto ta1_y_yyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 546);

    auto ta1_y_yyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 547);

    auto ta1_y_yyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 548);

    auto ta1_y_yyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 549);

    auto ta1_y_yyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 550);

    auto ta1_y_yyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 551);

    auto ta1_y_yyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 552);

    auto ta1_y_yyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 553);

    auto ta1_y_yyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 554);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_y_yyy_xxxx_0,   \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxy_0,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxxz_0,   \
                             ta1_y_yyy_xxxz_1,   \
                             ta1_y_yyy_xxyy_0,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyz_0,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xxzz_0,   \
                             ta1_y_yyy_xxzz_1,   \
                             ta1_y_yyy_xyyy_0,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyz_0,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyzz_0,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_xzzz_0,   \
                             ta1_y_yyy_xzzz_1,   \
                             ta1_y_yyy_yyyy_0,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyz_0,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyzz_0,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yzzz_0,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_zzzz_0,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta1_y_yyyy_xxx_0,   \
                             ta1_y_yyyy_xxx_1,   \
                             ta1_y_yyyy_xxxx_0,  \
                             ta1_y_yyyy_xxxx_1,  \
                             ta1_y_yyyy_xxxy_0,  \
                             ta1_y_yyyy_xxxy_1,  \
                             ta1_y_yyyy_xxxz_0,  \
                             ta1_y_yyyy_xxxz_1,  \
                             ta1_y_yyyy_xxy_0,   \
                             ta1_y_yyyy_xxy_1,   \
                             ta1_y_yyyy_xxyy_0,  \
                             ta1_y_yyyy_xxyy_1,  \
                             ta1_y_yyyy_xxyz_0,  \
                             ta1_y_yyyy_xxyz_1,  \
                             ta1_y_yyyy_xxz_0,   \
                             ta1_y_yyyy_xxz_1,   \
                             ta1_y_yyyy_xxzz_0,  \
                             ta1_y_yyyy_xxzz_1,  \
                             ta1_y_yyyy_xyy_0,   \
                             ta1_y_yyyy_xyy_1,   \
                             ta1_y_yyyy_xyyy_0,  \
                             ta1_y_yyyy_xyyy_1,  \
                             ta1_y_yyyy_xyyz_0,  \
                             ta1_y_yyyy_xyyz_1,  \
                             ta1_y_yyyy_xyz_0,   \
                             ta1_y_yyyy_xyz_1,   \
                             ta1_y_yyyy_xyzz_0,  \
                             ta1_y_yyyy_xyzz_1,  \
                             ta1_y_yyyy_xzz_0,   \
                             ta1_y_yyyy_xzz_1,   \
                             ta1_y_yyyy_xzzz_0,  \
                             ta1_y_yyyy_xzzz_1,  \
                             ta1_y_yyyy_yyy_0,   \
                             ta1_y_yyyy_yyy_1,   \
                             ta1_y_yyyy_yyyy_0,  \
                             ta1_y_yyyy_yyyy_1,  \
                             ta1_y_yyyy_yyyz_0,  \
                             ta1_y_yyyy_yyyz_1,  \
                             ta1_y_yyyy_yyz_0,   \
                             ta1_y_yyyy_yyz_1,   \
                             ta1_y_yyyy_yyzz_0,  \
                             ta1_y_yyyy_yyzz_1,  \
                             ta1_y_yyyy_yzz_0,   \
                             ta1_y_yyyy_yzz_1,   \
                             ta1_y_yyyy_yzzz_0,  \
                             ta1_y_yyyy_yzzz_1,  \
                             ta1_y_yyyy_zzz_0,   \
                             ta1_y_yyyy_zzz_1,   \
                             ta1_y_yyyy_zzzz_0,  \
                             ta1_y_yyyy_zzzz_1,  \
                             ta1_y_yyyyy_xxxx_0, \
                             ta1_y_yyyyy_xxxy_0, \
                             ta1_y_yyyyy_xxxz_0, \
                             ta1_y_yyyyy_xxyy_0, \
                             ta1_y_yyyyy_xxyz_0, \
                             ta1_y_yyyyy_xxzz_0, \
                             ta1_y_yyyyy_xyyy_0, \
                             ta1_y_yyyyy_xyyz_0, \
                             ta1_y_yyyyy_xyzz_0, \
                             ta1_y_yyyyy_xzzz_0, \
                             ta1_y_yyyyy_yyyy_0, \
                             ta1_y_yyyyy_yyyz_0, \
                             ta1_y_yyyyy_yyzz_0, \
                             ta1_y_yyyyy_yzzz_0, \
                             ta1_y_yyyyy_zzzz_0, \
                             ta_yyyy_xxxx_1,     \
                             ta_yyyy_xxxy_1,     \
                             ta_yyyy_xxxz_1,     \
                             ta_yyyy_xxyy_1,     \
                             ta_yyyy_xxyz_1,     \
                             ta_yyyy_xxzz_1,     \
                             ta_yyyy_xyyy_1,     \
                             ta_yyyy_xyyz_1,     \
                             ta_yyyy_xyzz_1,     \
                             ta_yyyy_xzzz_1,     \
                             ta_yyyy_yyyy_1,     \
                             ta_yyyy_yyyz_1,     \
                             ta_yyyy_yyzz_1,     \
                             ta_yyyy_yzzz_1,     \
                             ta_yyyy_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyy_xxxx_0[i] = 4.0 * ta1_y_yyy_xxxx_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxx_1[i] * fe_0 + ta_yyyy_xxxx_1[i] +
                                ta1_y_yyyy_xxxx_0[i] * pa_y[i] - ta1_y_yyyy_xxxx_1[i] * pc_y[i];

        ta1_y_yyyyy_xxxy_0[i] = 4.0 * ta1_y_yyy_xxxy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxy_1[i] * fe_0 + ta1_y_yyyy_xxx_0[i] * fe_0 -
                                ta1_y_yyyy_xxx_1[i] * fe_0 + ta_yyyy_xxxy_1[i] + ta1_y_yyyy_xxxy_0[i] * pa_y[i] - ta1_y_yyyy_xxxy_1[i] * pc_y[i];

        ta1_y_yyyyy_xxxz_0[i] = 4.0 * ta1_y_yyy_xxxz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxz_1[i] * fe_0 + ta_yyyy_xxxz_1[i] +
                                ta1_y_yyyy_xxxz_0[i] * pa_y[i] - ta1_y_yyyy_xxxz_1[i] * pc_y[i];

        ta1_y_yyyyy_xxyy_0[i] = 4.0 * ta1_y_yyy_xxyy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_yyyy_xxy_0[i] * fe_0 -
                                2.0 * ta1_y_yyyy_xxy_1[i] * fe_0 + ta_yyyy_xxyy_1[i] + ta1_y_yyyy_xxyy_0[i] * pa_y[i] -
                                ta1_y_yyyy_xxyy_1[i] * pc_y[i];

        ta1_y_yyyyy_xxyz_0[i] = 4.0 * ta1_y_yyy_xxyz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxyz_1[i] * fe_0 + ta1_y_yyyy_xxz_0[i] * fe_0 -
                                ta1_y_yyyy_xxz_1[i] * fe_0 + ta_yyyy_xxyz_1[i] + ta1_y_yyyy_xxyz_0[i] * pa_y[i] - ta1_y_yyyy_xxyz_1[i] * pc_y[i];

        ta1_y_yyyyy_xxzz_0[i] = 4.0 * ta1_y_yyy_xxzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxzz_1[i] * fe_0 + ta_yyyy_xxzz_1[i] +
                                ta1_y_yyyy_xxzz_0[i] * pa_y[i] - ta1_y_yyyy_xxzz_1[i] * pc_y[i];

        ta1_y_yyyyy_xyyy_0[i] = 4.0 * ta1_y_yyy_xyyy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xyyy_1[i] * fe_0 + 3.0 * ta1_y_yyyy_xyy_0[i] * fe_0 -
                                3.0 * ta1_y_yyyy_xyy_1[i] * fe_0 + ta_yyyy_xyyy_1[i] + ta1_y_yyyy_xyyy_0[i] * pa_y[i] -
                                ta1_y_yyyy_xyyy_1[i] * pc_y[i];

        ta1_y_yyyyy_xyyz_0[i] = 4.0 * ta1_y_yyy_xyyz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xyyz_1[i] * fe_0 + 2.0 * ta1_y_yyyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyy_xyz_1[i] * fe_0 + ta_yyyy_xyyz_1[i] + ta1_y_yyyy_xyyz_0[i] * pa_y[i] -
                                ta1_y_yyyy_xyyz_1[i] * pc_y[i];

        ta1_y_yyyyy_xyzz_0[i] = 4.0 * ta1_y_yyy_xyzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xyzz_1[i] * fe_0 + ta1_y_yyyy_xzz_0[i] * fe_0 -
                                ta1_y_yyyy_xzz_1[i] * fe_0 + ta_yyyy_xyzz_1[i] + ta1_y_yyyy_xyzz_0[i] * pa_y[i] - ta1_y_yyyy_xyzz_1[i] * pc_y[i];

        ta1_y_yyyyy_xzzz_0[i] = 4.0 * ta1_y_yyy_xzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xzzz_1[i] * fe_0 + ta_yyyy_xzzz_1[i] +
                                ta1_y_yyyy_xzzz_0[i] * pa_y[i] - ta1_y_yyyy_xzzz_1[i] * pc_y[i];

        ta1_y_yyyyy_yyyy_0[i] = 4.0 * ta1_y_yyy_yyyy_0[i] * fe_0 - 4.0 * ta1_y_yyy_yyyy_1[i] * fe_0 + 4.0 * ta1_y_yyyy_yyy_0[i] * fe_0 -
                                4.0 * ta1_y_yyyy_yyy_1[i] * fe_0 + ta_yyyy_yyyy_1[i] + ta1_y_yyyy_yyyy_0[i] * pa_y[i] -
                                ta1_y_yyyy_yyyy_1[i] * pc_y[i];

        ta1_y_yyyyy_yyyz_0[i] = 4.0 * ta1_y_yyy_yyyz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yyyz_1[i] * fe_0 + 3.0 * ta1_y_yyyy_yyz_0[i] * fe_0 -
                                3.0 * ta1_y_yyyy_yyz_1[i] * fe_0 + ta_yyyy_yyyz_1[i] + ta1_y_yyyy_yyyz_0[i] * pa_y[i] -
                                ta1_y_yyyy_yyyz_1[i] * pc_y[i];

        ta1_y_yyyyy_yyzz_0[i] = 4.0 * ta1_y_yyy_yyzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyy_yzz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyy_yzz_1[i] * fe_0 + ta_yyyy_yyzz_1[i] + ta1_y_yyyy_yyzz_0[i] * pa_y[i] -
                                ta1_y_yyyy_yyzz_1[i] * pc_y[i];

        ta1_y_yyyyy_yzzz_0[i] = 4.0 * ta1_y_yyy_yzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yzzz_1[i] * fe_0 + ta1_y_yyyy_zzz_0[i] * fe_0 -
                                ta1_y_yyyy_zzz_1[i] * fe_0 + ta_yyyy_yzzz_1[i] + ta1_y_yyyy_yzzz_0[i] * pa_y[i] - ta1_y_yyyy_yzzz_1[i] * pc_y[i];

        ta1_y_yyyyy_zzzz_0[i] = 4.0 * ta1_y_yyy_zzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_zzzz_1[i] * fe_0 + ta_yyyy_zzzz_1[i] +
                                ta1_y_yyyy_zzzz_0[i] * pa_y[i] - ta1_y_yyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 555-570 components of targeted buffer : HG

    auto ta1_y_yyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 555);

    auto ta1_y_yyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 556);

    auto ta1_y_yyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 557);

    auto ta1_y_yyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 558);

    auto ta1_y_yyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 559);

    auto ta1_y_yyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 560);

    auto ta1_y_yyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 561);

    auto ta1_y_yyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 562);

    auto ta1_y_yyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 563);

    auto ta1_y_yyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 564);

    auto ta1_y_yyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 565);

    auto ta1_y_yyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 566);

    auto ta1_y_yyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 567);

    auto ta1_y_yyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 568);

    auto ta1_y_yyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 569);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_yyyy_xxx_0,   \
                             ta1_y_yyyy_xxx_1,   \
                             ta1_y_yyyy_xxxx_0,  \
                             ta1_y_yyyy_xxxx_1,  \
                             ta1_y_yyyy_xxxy_0,  \
                             ta1_y_yyyy_xxxy_1,  \
                             ta1_y_yyyy_xxxz_0,  \
                             ta1_y_yyyy_xxxz_1,  \
                             ta1_y_yyyy_xxy_0,   \
                             ta1_y_yyyy_xxy_1,   \
                             ta1_y_yyyy_xxyy_0,  \
                             ta1_y_yyyy_xxyy_1,  \
                             ta1_y_yyyy_xxyz_0,  \
                             ta1_y_yyyy_xxyz_1,  \
                             ta1_y_yyyy_xxz_0,   \
                             ta1_y_yyyy_xxz_1,   \
                             ta1_y_yyyy_xxzz_0,  \
                             ta1_y_yyyy_xxzz_1,  \
                             ta1_y_yyyy_xyy_0,   \
                             ta1_y_yyyy_xyy_1,   \
                             ta1_y_yyyy_xyyy_0,  \
                             ta1_y_yyyy_xyyy_1,  \
                             ta1_y_yyyy_xyyz_0,  \
                             ta1_y_yyyy_xyyz_1,  \
                             ta1_y_yyyy_xyz_0,   \
                             ta1_y_yyyy_xyz_1,   \
                             ta1_y_yyyy_xyzz_0,  \
                             ta1_y_yyyy_xyzz_1,  \
                             ta1_y_yyyy_xzz_0,   \
                             ta1_y_yyyy_xzz_1,   \
                             ta1_y_yyyy_xzzz_0,  \
                             ta1_y_yyyy_xzzz_1,  \
                             ta1_y_yyyy_yyy_0,   \
                             ta1_y_yyyy_yyy_1,   \
                             ta1_y_yyyy_yyyy_0,  \
                             ta1_y_yyyy_yyyy_1,  \
                             ta1_y_yyyy_yyyz_0,  \
                             ta1_y_yyyy_yyyz_1,  \
                             ta1_y_yyyy_yyz_0,   \
                             ta1_y_yyyy_yyz_1,   \
                             ta1_y_yyyy_yyzz_0,  \
                             ta1_y_yyyy_yyzz_1,  \
                             ta1_y_yyyy_yzz_0,   \
                             ta1_y_yyyy_yzz_1,   \
                             ta1_y_yyyy_yzzz_0,  \
                             ta1_y_yyyy_yzzz_1,  \
                             ta1_y_yyyy_zzz_0,   \
                             ta1_y_yyyy_zzz_1,   \
                             ta1_y_yyyy_zzzz_0,  \
                             ta1_y_yyyy_zzzz_1,  \
                             ta1_y_yyyyz_xxxx_0, \
                             ta1_y_yyyyz_xxxy_0, \
                             ta1_y_yyyyz_xxxz_0, \
                             ta1_y_yyyyz_xxyy_0, \
                             ta1_y_yyyyz_xxyz_0, \
                             ta1_y_yyyyz_xxzz_0, \
                             ta1_y_yyyyz_xyyy_0, \
                             ta1_y_yyyyz_xyyz_0, \
                             ta1_y_yyyyz_xyzz_0, \
                             ta1_y_yyyyz_xzzz_0, \
                             ta1_y_yyyyz_yyyy_0, \
                             ta1_y_yyyyz_yyyz_0, \
                             ta1_y_yyyyz_yyzz_0, \
                             ta1_y_yyyyz_yzzz_0, \
                             ta1_y_yyyyz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyz_xxxx_0[i] = ta1_y_yyyy_xxxx_0[i] * pa_z[i] - ta1_y_yyyy_xxxx_1[i] * pc_z[i];

        ta1_y_yyyyz_xxxy_0[i] = ta1_y_yyyy_xxxy_0[i] * pa_z[i] - ta1_y_yyyy_xxxy_1[i] * pc_z[i];

        ta1_y_yyyyz_xxxz_0[i] =
            ta1_y_yyyy_xxx_0[i] * fe_0 - ta1_y_yyyy_xxx_1[i] * fe_0 + ta1_y_yyyy_xxxz_0[i] * pa_z[i] - ta1_y_yyyy_xxxz_1[i] * pc_z[i];

        ta1_y_yyyyz_xxyy_0[i] = ta1_y_yyyy_xxyy_0[i] * pa_z[i] - ta1_y_yyyy_xxyy_1[i] * pc_z[i];

        ta1_y_yyyyz_xxyz_0[i] =
            ta1_y_yyyy_xxy_0[i] * fe_0 - ta1_y_yyyy_xxy_1[i] * fe_0 + ta1_y_yyyy_xxyz_0[i] * pa_z[i] - ta1_y_yyyy_xxyz_1[i] * pc_z[i];

        ta1_y_yyyyz_xxzz_0[i] =
            2.0 * ta1_y_yyyy_xxz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xxz_1[i] * fe_0 + ta1_y_yyyy_xxzz_0[i] * pa_z[i] - ta1_y_yyyy_xxzz_1[i] * pc_z[i];

        ta1_y_yyyyz_xyyy_0[i] = ta1_y_yyyy_xyyy_0[i] * pa_z[i] - ta1_y_yyyy_xyyy_1[i] * pc_z[i];

        ta1_y_yyyyz_xyyz_0[i] =
            ta1_y_yyyy_xyy_0[i] * fe_0 - ta1_y_yyyy_xyy_1[i] * fe_0 + ta1_y_yyyy_xyyz_0[i] * pa_z[i] - ta1_y_yyyy_xyyz_1[i] * pc_z[i];

        ta1_y_yyyyz_xyzz_0[i] =
            2.0 * ta1_y_yyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xyz_1[i] * fe_0 + ta1_y_yyyy_xyzz_0[i] * pa_z[i] - ta1_y_yyyy_xyzz_1[i] * pc_z[i];

        ta1_y_yyyyz_xzzz_0[i] =
            3.0 * ta1_y_yyyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_yyyy_xzz_1[i] * fe_0 + ta1_y_yyyy_xzzz_0[i] * pa_z[i] - ta1_y_yyyy_xzzz_1[i] * pc_z[i];

        ta1_y_yyyyz_yyyy_0[i] = ta1_y_yyyy_yyyy_0[i] * pa_z[i] - ta1_y_yyyy_yyyy_1[i] * pc_z[i];

        ta1_y_yyyyz_yyyz_0[i] =
            ta1_y_yyyy_yyy_0[i] * fe_0 - ta1_y_yyyy_yyy_1[i] * fe_0 + ta1_y_yyyy_yyyz_0[i] * pa_z[i] - ta1_y_yyyy_yyyz_1[i] * pc_z[i];

        ta1_y_yyyyz_yyzz_0[i] =
            2.0 * ta1_y_yyyy_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_yyz_1[i] * fe_0 + ta1_y_yyyy_yyzz_0[i] * pa_z[i] - ta1_y_yyyy_yyzz_1[i] * pc_z[i];

        ta1_y_yyyyz_yzzz_0[i] =
            3.0 * ta1_y_yyyy_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyyy_yzz_1[i] * fe_0 + ta1_y_yyyy_yzzz_0[i] * pa_z[i] - ta1_y_yyyy_yzzz_1[i] * pc_z[i];

        ta1_y_yyyyz_zzzz_0[i] =
            4.0 * ta1_y_yyyy_zzz_0[i] * fe_0 - 4.0 * ta1_y_yyyy_zzz_1[i] * fe_0 + ta1_y_yyyy_zzzz_0[i] * pa_z[i] - ta1_y_yyyy_zzzz_1[i] * pc_z[i];
    }

    // Set up 570-585 components of targeted buffer : HG

    auto ta1_y_yyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 570);

    auto ta1_y_yyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 571);

    auto ta1_y_yyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 572);

    auto ta1_y_yyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 573);

    auto ta1_y_yyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 574);

    auto ta1_y_yyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 575);

    auto ta1_y_yyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 576);

    auto ta1_y_yyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 577);

    auto ta1_y_yyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 578);

    auto ta1_y_yyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 579);

    auto ta1_y_yyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 580);

    auto ta1_y_yyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 581);

    auto ta1_y_yyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 582);

    auto ta1_y_yyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 583);

    auto ta1_y_yyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 584);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyy_xxxx_0,   \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxy_0,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxyy_0,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyz_0,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xyyy_0,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyz_0,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyzz_0,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_yyyy_0,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyz_0,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyzz_0,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yzzz_0,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyyz_xxxx_0,  \
                             ta1_y_yyyz_xxxx_1,  \
                             ta1_y_yyyz_xxxy_0,  \
                             ta1_y_yyyz_xxxy_1,  \
                             ta1_y_yyyz_xxy_0,   \
                             ta1_y_yyyz_xxy_1,   \
                             ta1_y_yyyz_xxyy_0,  \
                             ta1_y_yyyz_xxyy_1,  \
                             ta1_y_yyyz_xxyz_0,  \
                             ta1_y_yyyz_xxyz_1,  \
                             ta1_y_yyyz_xyy_0,   \
                             ta1_y_yyyz_xyy_1,   \
                             ta1_y_yyyz_xyyy_0,  \
                             ta1_y_yyyz_xyyy_1,  \
                             ta1_y_yyyz_xyyz_0,  \
                             ta1_y_yyyz_xyyz_1,  \
                             ta1_y_yyyz_xyz_0,   \
                             ta1_y_yyyz_xyz_1,   \
                             ta1_y_yyyz_xyzz_0,  \
                             ta1_y_yyyz_xyzz_1,  \
                             ta1_y_yyyz_yyy_0,   \
                             ta1_y_yyyz_yyy_1,   \
                             ta1_y_yyyz_yyyy_0,  \
                             ta1_y_yyyz_yyyy_1,  \
                             ta1_y_yyyz_yyyz_0,  \
                             ta1_y_yyyz_yyyz_1,  \
                             ta1_y_yyyz_yyz_0,   \
                             ta1_y_yyyz_yyz_1,   \
                             ta1_y_yyyz_yyzz_0,  \
                             ta1_y_yyyz_yyzz_1,  \
                             ta1_y_yyyz_yzz_0,   \
                             ta1_y_yyyz_yzz_1,   \
                             ta1_y_yyyz_yzzz_0,  \
                             ta1_y_yyyz_yzzz_1,  \
                             ta1_y_yyyzz_xxxx_0, \
                             ta1_y_yyyzz_xxxy_0, \
                             ta1_y_yyyzz_xxxz_0, \
                             ta1_y_yyyzz_xxyy_0, \
                             ta1_y_yyyzz_xxyz_0, \
                             ta1_y_yyyzz_xxzz_0, \
                             ta1_y_yyyzz_xyyy_0, \
                             ta1_y_yyyzz_xyyz_0, \
                             ta1_y_yyyzz_xyzz_0, \
                             ta1_y_yyyzz_xzzz_0, \
                             ta1_y_yyyzz_yyyy_0, \
                             ta1_y_yyyzz_yyyz_0, \
                             ta1_y_yyyzz_yyzz_0, \
                             ta1_y_yyyzz_yzzz_0, \
                             ta1_y_yyyzz_zzzz_0, \
                             ta1_y_yyzz_xxxz_0,  \
                             ta1_y_yyzz_xxxz_1,  \
                             ta1_y_yyzz_xxzz_0,  \
                             ta1_y_yyzz_xxzz_1,  \
                             ta1_y_yyzz_xzzz_0,  \
                             ta1_y_yyzz_xzzz_1,  \
                             ta1_y_yyzz_zzzz_0,  \
                             ta1_y_yyzz_zzzz_1,  \
                             ta1_y_yzz_xxxz_0,   \
                             ta1_y_yzz_xxxz_1,   \
                             ta1_y_yzz_xxzz_0,   \
                             ta1_y_yzz_xxzz_1,   \
                             ta1_y_yzz_xzzz_0,   \
                             ta1_y_yzz_xzzz_1,   \
                             ta1_y_yzz_zzzz_0,   \
                             ta1_y_yzz_zzzz_1,   \
                             ta_yyzz_xxxz_1,     \
                             ta_yyzz_xxzz_1,     \
                             ta_yyzz_xzzz_1,     \
                             ta_yyzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzz_xxxx_0[i] =
            ta1_y_yyy_xxxx_0[i] * fe_0 - ta1_y_yyy_xxxx_1[i] * fe_0 + ta1_y_yyyz_xxxx_0[i] * pa_z[i] - ta1_y_yyyz_xxxx_1[i] * pc_z[i];

        ta1_y_yyyzz_xxxy_0[i] =
            ta1_y_yyy_xxxy_0[i] * fe_0 - ta1_y_yyy_xxxy_1[i] * fe_0 + ta1_y_yyyz_xxxy_0[i] * pa_z[i] - ta1_y_yyyz_xxxy_1[i] * pc_z[i];

        ta1_y_yyyzz_xxxz_0[i] = 2.0 * ta1_y_yzz_xxxz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xxxz_1[i] * fe_0 + ta_yyzz_xxxz_1[i] +
                                ta1_y_yyzz_xxxz_0[i] * pa_y[i] - ta1_y_yyzz_xxxz_1[i] * pc_y[i];

        ta1_y_yyyzz_xxyy_0[i] =
            ta1_y_yyy_xxyy_0[i] * fe_0 - ta1_y_yyy_xxyy_1[i] * fe_0 + ta1_y_yyyz_xxyy_0[i] * pa_z[i] - ta1_y_yyyz_xxyy_1[i] * pc_z[i];

        ta1_y_yyyzz_xxyz_0[i] = ta1_y_yyy_xxyz_0[i] * fe_0 - ta1_y_yyy_xxyz_1[i] * fe_0 + ta1_y_yyyz_xxy_0[i] * fe_0 - ta1_y_yyyz_xxy_1[i] * fe_0 +
                                ta1_y_yyyz_xxyz_0[i] * pa_z[i] - ta1_y_yyyz_xxyz_1[i] * pc_z[i];

        ta1_y_yyyzz_xxzz_0[i] = 2.0 * ta1_y_yzz_xxzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xxzz_1[i] * fe_0 + ta_yyzz_xxzz_1[i] +
                                ta1_y_yyzz_xxzz_0[i] * pa_y[i] - ta1_y_yyzz_xxzz_1[i] * pc_y[i];

        ta1_y_yyyzz_xyyy_0[i] =
            ta1_y_yyy_xyyy_0[i] * fe_0 - ta1_y_yyy_xyyy_1[i] * fe_0 + ta1_y_yyyz_xyyy_0[i] * pa_z[i] - ta1_y_yyyz_xyyy_1[i] * pc_z[i];

        ta1_y_yyyzz_xyyz_0[i] = ta1_y_yyy_xyyz_0[i] * fe_0 - ta1_y_yyy_xyyz_1[i] * fe_0 + ta1_y_yyyz_xyy_0[i] * fe_0 - ta1_y_yyyz_xyy_1[i] * fe_0 +
                                ta1_y_yyyz_xyyz_0[i] * pa_z[i] - ta1_y_yyyz_xyyz_1[i] * pc_z[i];

        ta1_y_yyyzz_xyzz_0[i] = ta1_y_yyy_xyzz_0[i] * fe_0 - ta1_y_yyy_xyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyz_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyz_xyz_1[i] * fe_0 + ta1_y_yyyz_xyzz_0[i] * pa_z[i] - ta1_y_yyyz_xyzz_1[i] * pc_z[i];

        ta1_y_yyyzz_xzzz_0[i] = 2.0 * ta1_y_yzz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xzzz_1[i] * fe_0 + ta_yyzz_xzzz_1[i] +
                                ta1_y_yyzz_xzzz_0[i] * pa_y[i] - ta1_y_yyzz_xzzz_1[i] * pc_y[i];

        ta1_y_yyyzz_yyyy_0[i] =
            ta1_y_yyy_yyyy_0[i] * fe_0 - ta1_y_yyy_yyyy_1[i] * fe_0 + ta1_y_yyyz_yyyy_0[i] * pa_z[i] - ta1_y_yyyz_yyyy_1[i] * pc_z[i];

        ta1_y_yyyzz_yyyz_0[i] = ta1_y_yyy_yyyz_0[i] * fe_0 - ta1_y_yyy_yyyz_1[i] * fe_0 + ta1_y_yyyz_yyy_0[i] * fe_0 - ta1_y_yyyz_yyy_1[i] * fe_0 +
                                ta1_y_yyyz_yyyz_0[i] * pa_z[i] - ta1_y_yyyz_yyyz_1[i] * pc_z[i];

        ta1_y_yyyzz_yyzz_0[i] = ta1_y_yyy_yyzz_0[i] * fe_0 - ta1_y_yyy_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyz_yyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyz_yyz_1[i] * fe_0 + ta1_y_yyyz_yyzz_0[i] * pa_z[i] - ta1_y_yyyz_yyzz_1[i] * pc_z[i];

        ta1_y_yyyzz_yzzz_0[i] = ta1_y_yyy_yzzz_0[i] * fe_0 - ta1_y_yyy_yzzz_1[i] * fe_0 + 3.0 * ta1_y_yyyz_yzz_0[i] * fe_0 -
                                3.0 * ta1_y_yyyz_yzz_1[i] * fe_0 + ta1_y_yyyz_yzzz_0[i] * pa_z[i] - ta1_y_yyyz_yzzz_1[i] * pc_z[i];

        ta1_y_yyyzz_zzzz_0[i] = 2.0 * ta1_y_yzz_zzzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_zzzz_1[i] * fe_0 + ta_yyzz_zzzz_1[i] +
                                ta1_y_yyzz_zzzz_0[i] * pa_y[i] - ta1_y_yyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 585-600 components of targeted buffer : HG

    auto ta1_y_yyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 585);

    auto ta1_y_yyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 586);

    auto ta1_y_yyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 587);

    auto ta1_y_yyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 588);

    auto ta1_y_yyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 589);

    auto ta1_y_yyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 590);

    auto ta1_y_yyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 591);

    auto ta1_y_yyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 592);

    auto ta1_y_yyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 593);

    auto ta1_y_yyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 594);

    auto ta1_y_yyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 595);

    auto ta1_y_yyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 596);

    auto ta1_y_yyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 597);

    auto ta1_y_yyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 598);

    auto ta1_y_yyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 599);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyz_xxxx_0,   \
                             ta1_y_yyz_xxxx_1,   \
                             ta1_y_yyz_xxxy_0,   \
                             ta1_y_yyz_xxxy_1,   \
                             ta1_y_yyz_xxyy_0,   \
                             ta1_y_yyz_xxyy_1,   \
                             ta1_y_yyz_xxyz_0,   \
                             ta1_y_yyz_xxyz_1,   \
                             ta1_y_yyz_xyyy_0,   \
                             ta1_y_yyz_xyyy_1,   \
                             ta1_y_yyz_xyyz_0,   \
                             ta1_y_yyz_xyyz_1,   \
                             ta1_y_yyz_xyzz_0,   \
                             ta1_y_yyz_xyzz_1,   \
                             ta1_y_yyz_yyyy_0,   \
                             ta1_y_yyz_yyyy_1,   \
                             ta1_y_yyz_yyyz_0,   \
                             ta1_y_yyz_yyyz_1,   \
                             ta1_y_yyz_yyzz_0,   \
                             ta1_y_yyz_yyzz_1,   \
                             ta1_y_yyz_yzzz_0,   \
                             ta1_y_yyz_yzzz_1,   \
                             ta1_y_yyzz_xxxx_0,  \
                             ta1_y_yyzz_xxxx_1,  \
                             ta1_y_yyzz_xxxy_0,  \
                             ta1_y_yyzz_xxxy_1,  \
                             ta1_y_yyzz_xxy_0,   \
                             ta1_y_yyzz_xxy_1,   \
                             ta1_y_yyzz_xxyy_0,  \
                             ta1_y_yyzz_xxyy_1,  \
                             ta1_y_yyzz_xxyz_0,  \
                             ta1_y_yyzz_xxyz_1,  \
                             ta1_y_yyzz_xyy_0,   \
                             ta1_y_yyzz_xyy_1,   \
                             ta1_y_yyzz_xyyy_0,  \
                             ta1_y_yyzz_xyyy_1,  \
                             ta1_y_yyzz_xyyz_0,  \
                             ta1_y_yyzz_xyyz_1,  \
                             ta1_y_yyzz_xyz_0,   \
                             ta1_y_yyzz_xyz_1,   \
                             ta1_y_yyzz_xyzz_0,  \
                             ta1_y_yyzz_xyzz_1,  \
                             ta1_y_yyzz_yyy_0,   \
                             ta1_y_yyzz_yyy_1,   \
                             ta1_y_yyzz_yyyy_0,  \
                             ta1_y_yyzz_yyyy_1,  \
                             ta1_y_yyzz_yyyz_0,  \
                             ta1_y_yyzz_yyyz_1,  \
                             ta1_y_yyzz_yyz_0,   \
                             ta1_y_yyzz_yyz_1,   \
                             ta1_y_yyzz_yyzz_0,  \
                             ta1_y_yyzz_yyzz_1,  \
                             ta1_y_yyzz_yzz_0,   \
                             ta1_y_yyzz_yzz_1,   \
                             ta1_y_yyzz_yzzz_0,  \
                             ta1_y_yyzz_yzzz_1,  \
                             ta1_y_yyzzz_xxxx_0, \
                             ta1_y_yyzzz_xxxy_0, \
                             ta1_y_yyzzz_xxxz_0, \
                             ta1_y_yyzzz_xxyy_0, \
                             ta1_y_yyzzz_xxyz_0, \
                             ta1_y_yyzzz_xxzz_0, \
                             ta1_y_yyzzz_xyyy_0, \
                             ta1_y_yyzzz_xyyz_0, \
                             ta1_y_yyzzz_xyzz_0, \
                             ta1_y_yyzzz_xzzz_0, \
                             ta1_y_yyzzz_yyyy_0, \
                             ta1_y_yyzzz_yyyz_0, \
                             ta1_y_yyzzz_yyzz_0, \
                             ta1_y_yyzzz_yzzz_0, \
                             ta1_y_yyzzz_zzzz_0, \
                             ta1_y_yzzz_xxxz_0,  \
                             ta1_y_yzzz_xxxz_1,  \
                             ta1_y_yzzz_xxzz_0,  \
                             ta1_y_yzzz_xxzz_1,  \
                             ta1_y_yzzz_xzzz_0,  \
                             ta1_y_yzzz_xzzz_1,  \
                             ta1_y_yzzz_zzzz_0,  \
                             ta1_y_yzzz_zzzz_1,  \
                             ta1_y_zzz_xxxz_0,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxzz_0,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xzzz_0,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_zzzz_0,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta_yzzz_xxxz_1,     \
                             ta_yzzz_xxzz_1,     \
                             ta_yzzz_xzzz_1,     \
                             ta_yzzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzz_xxxx_0[i] =
            2.0 * ta1_y_yyz_xxxx_0[i] * fe_0 - 2.0 * ta1_y_yyz_xxxx_1[i] * fe_0 + ta1_y_yyzz_xxxx_0[i] * pa_z[i] - ta1_y_yyzz_xxxx_1[i] * pc_z[i];

        ta1_y_yyzzz_xxxy_0[i] =
            2.0 * ta1_y_yyz_xxxy_0[i] * fe_0 - 2.0 * ta1_y_yyz_xxxy_1[i] * fe_0 + ta1_y_yyzz_xxxy_0[i] * pa_z[i] - ta1_y_yyzz_xxxy_1[i] * pc_z[i];

        ta1_y_yyzzz_xxxz_0[i] = ta1_y_zzz_xxxz_0[i] * fe_0 - ta1_y_zzz_xxxz_1[i] * fe_0 + ta_yzzz_xxxz_1[i] + ta1_y_yzzz_xxxz_0[i] * pa_y[i] -
                                ta1_y_yzzz_xxxz_1[i] * pc_y[i];

        ta1_y_yyzzz_xxyy_0[i] =
            2.0 * ta1_y_yyz_xxyy_0[i] * fe_0 - 2.0 * ta1_y_yyz_xxyy_1[i] * fe_0 + ta1_y_yyzz_xxyy_0[i] * pa_z[i] - ta1_y_yyzz_xxyy_1[i] * pc_z[i];

        ta1_y_yyzzz_xxyz_0[i] = 2.0 * ta1_y_yyz_xxyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xxyz_1[i] * fe_0 + ta1_y_yyzz_xxy_0[i] * fe_0 -
                                ta1_y_yyzz_xxy_1[i] * fe_0 + ta1_y_yyzz_xxyz_0[i] * pa_z[i] - ta1_y_yyzz_xxyz_1[i] * pc_z[i];

        ta1_y_yyzzz_xxzz_0[i] = ta1_y_zzz_xxzz_0[i] * fe_0 - ta1_y_zzz_xxzz_1[i] * fe_0 + ta_yzzz_xxzz_1[i] + ta1_y_yzzz_xxzz_0[i] * pa_y[i] -
                                ta1_y_yzzz_xxzz_1[i] * pc_y[i];

        ta1_y_yyzzz_xyyy_0[i] =
            2.0 * ta1_y_yyz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyyy_1[i] * fe_0 + ta1_y_yyzz_xyyy_0[i] * pa_z[i] - ta1_y_yyzz_xyyy_1[i] * pc_z[i];

        ta1_y_yyzzz_xyyz_0[i] = 2.0 * ta1_y_yyz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyyz_1[i] * fe_0 + ta1_y_yyzz_xyy_0[i] * fe_0 -
                                ta1_y_yyzz_xyy_1[i] * fe_0 + ta1_y_yyzz_xyyz_0[i] * pa_z[i] - ta1_y_yyzz_xyyz_1[i] * pc_z[i];

        ta1_y_yyzzz_xyzz_0[i] = 2.0 * ta1_y_yyz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyzz_1[i] * fe_0 + 2.0 * ta1_y_yyzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyzz_xyz_1[i] * fe_0 + ta1_y_yyzz_xyzz_0[i] * pa_z[i] - ta1_y_yyzz_xyzz_1[i] * pc_z[i];

        ta1_y_yyzzz_xzzz_0[i] = ta1_y_zzz_xzzz_0[i] * fe_0 - ta1_y_zzz_xzzz_1[i] * fe_0 + ta_yzzz_xzzz_1[i] + ta1_y_yzzz_xzzz_0[i] * pa_y[i] -
                                ta1_y_yzzz_xzzz_1[i] * pc_y[i];

        ta1_y_yyzzz_yyyy_0[i] =
            2.0 * ta1_y_yyz_yyyy_0[i] * fe_0 - 2.0 * ta1_y_yyz_yyyy_1[i] * fe_0 + ta1_y_yyzz_yyyy_0[i] * pa_z[i] - ta1_y_yyzz_yyyy_1[i] * pc_z[i];

        ta1_y_yyzzz_yyyz_0[i] = 2.0 * ta1_y_yyz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yyyz_1[i] * fe_0 + ta1_y_yyzz_yyy_0[i] * fe_0 -
                                ta1_y_yyzz_yyy_1[i] * fe_0 + ta1_y_yyzz_yyyz_0[i] * pa_z[i] - ta1_y_yyzz_yyyz_1[i] * pc_z[i];

        ta1_y_yyzzz_yyzz_0[i] = 2.0 * ta1_y_yyz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyzz_yyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyzz_yyz_1[i] * fe_0 + ta1_y_yyzz_yyzz_0[i] * pa_z[i] - ta1_y_yyzz_yyzz_1[i] * pc_z[i];

        ta1_y_yyzzz_yzzz_0[i] = 2.0 * ta1_y_yyz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yzzz_1[i] * fe_0 + 3.0 * ta1_y_yyzz_yzz_0[i] * fe_0 -
                                3.0 * ta1_y_yyzz_yzz_1[i] * fe_0 + ta1_y_yyzz_yzzz_0[i] * pa_z[i] - ta1_y_yyzz_yzzz_1[i] * pc_z[i];

        ta1_y_yyzzz_zzzz_0[i] = ta1_y_zzz_zzzz_0[i] * fe_0 - ta1_y_zzz_zzzz_1[i] * fe_0 + ta_yzzz_zzzz_1[i] + ta1_y_yzzz_zzzz_0[i] * pa_y[i] -
                                ta1_y_yzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 600-615 components of targeted buffer : HG

    auto ta1_y_yzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 600);

    auto ta1_y_yzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 601);

    auto ta1_y_yzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 602);

    auto ta1_y_yzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 603);

    auto ta1_y_yzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 604);

    auto ta1_y_yzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 605);

    auto ta1_y_yzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 606);

    auto ta1_y_yzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 607);

    auto ta1_y_yzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 608);

    auto ta1_y_yzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 609);

    auto ta1_y_yzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 610);

    auto ta1_y_yzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 611);

    auto ta1_y_yzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 612);

    auto ta1_y_yzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 613);

    auto ta1_y_yzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 614);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yzz_xxxy_0,   \
                             ta1_y_yzz_xxxy_1,   \
                             ta1_y_yzz_xxyy_0,   \
                             ta1_y_yzz_xxyy_1,   \
                             ta1_y_yzz_xyyy_0,   \
                             ta1_y_yzz_xyyy_1,   \
                             ta1_y_yzz_yyyy_0,   \
                             ta1_y_yzz_yyyy_1,   \
                             ta1_y_yzzz_xxxy_0,  \
                             ta1_y_yzzz_xxxy_1,  \
                             ta1_y_yzzz_xxyy_0,  \
                             ta1_y_yzzz_xxyy_1,  \
                             ta1_y_yzzz_xyyy_0,  \
                             ta1_y_yzzz_xyyy_1,  \
                             ta1_y_yzzz_yyyy_0,  \
                             ta1_y_yzzz_yyyy_1,  \
                             ta1_y_yzzzz_xxxx_0, \
                             ta1_y_yzzzz_xxxy_0, \
                             ta1_y_yzzzz_xxxz_0, \
                             ta1_y_yzzzz_xxyy_0, \
                             ta1_y_yzzzz_xxyz_0, \
                             ta1_y_yzzzz_xxzz_0, \
                             ta1_y_yzzzz_xyyy_0, \
                             ta1_y_yzzzz_xyyz_0, \
                             ta1_y_yzzzz_xyzz_0, \
                             ta1_y_yzzzz_xzzz_0, \
                             ta1_y_yzzzz_yyyy_0, \
                             ta1_y_yzzzz_yyyz_0, \
                             ta1_y_yzzzz_yyzz_0, \
                             ta1_y_yzzzz_yzzz_0, \
                             ta1_y_yzzzz_zzzz_0, \
                             ta1_y_zzzz_xxxx_0,  \
                             ta1_y_zzzz_xxxx_1,  \
                             ta1_y_zzzz_xxxz_0,  \
                             ta1_y_zzzz_xxxz_1,  \
                             ta1_y_zzzz_xxyz_0,  \
                             ta1_y_zzzz_xxyz_1,  \
                             ta1_y_zzzz_xxz_0,   \
                             ta1_y_zzzz_xxz_1,   \
                             ta1_y_zzzz_xxzz_0,  \
                             ta1_y_zzzz_xxzz_1,  \
                             ta1_y_zzzz_xyyz_0,  \
                             ta1_y_zzzz_xyyz_1,  \
                             ta1_y_zzzz_xyz_0,   \
                             ta1_y_zzzz_xyz_1,   \
                             ta1_y_zzzz_xyzz_0,  \
                             ta1_y_zzzz_xyzz_1,  \
                             ta1_y_zzzz_xzz_0,   \
                             ta1_y_zzzz_xzz_1,   \
                             ta1_y_zzzz_xzzz_0,  \
                             ta1_y_zzzz_xzzz_1,  \
                             ta1_y_zzzz_yyyz_0,  \
                             ta1_y_zzzz_yyyz_1,  \
                             ta1_y_zzzz_yyz_0,   \
                             ta1_y_zzzz_yyz_1,   \
                             ta1_y_zzzz_yyzz_0,  \
                             ta1_y_zzzz_yyzz_1,  \
                             ta1_y_zzzz_yzz_0,   \
                             ta1_y_zzzz_yzz_1,   \
                             ta1_y_zzzz_yzzz_0,  \
                             ta1_y_zzzz_yzzz_1,  \
                             ta1_y_zzzz_zzz_0,   \
                             ta1_y_zzzz_zzz_1,   \
                             ta1_y_zzzz_zzzz_0,  \
                             ta1_y_zzzz_zzzz_1,  \
                             ta_zzzz_xxxx_1,     \
                             ta_zzzz_xxxz_1,     \
                             ta_zzzz_xxyz_1,     \
                             ta_zzzz_xxzz_1,     \
                             ta_zzzz_xyyz_1,     \
                             ta_zzzz_xyzz_1,     \
                             ta_zzzz_xzzz_1,     \
                             ta_zzzz_yyyz_1,     \
                             ta_zzzz_yyzz_1,     \
                             ta_zzzz_yzzz_1,     \
                             ta_zzzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzz_xxxx_0[i] = ta_zzzz_xxxx_1[i] + ta1_y_zzzz_xxxx_0[i] * pa_y[i] - ta1_y_zzzz_xxxx_1[i] * pc_y[i];

        ta1_y_yzzzz_xxxy_0[i] =
            3.0 * ta1_y_yzz_xxxy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxxy_1[i] * fe_0 + ta1_y_yzzz_xxxy_0[i] * pa_z[i] - ta1_y_yzzz_xxxy_1[i] * pc_z[i];

        ta1_y_yzzzz_xxxz_0[i] = ta_zzzz_xxxz_1[i] + ta1_y_zzzz_xxxz_0[i] * pa_y[i] - ta1_y_zzzz_xxxz_1[i] * pc_y[i];

        ta1_y_yzzzz_xxyy_0[i] =
            3.0 * ta1_y_yzz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxyy_1[i] * fe_0 + ta1_y_yzzz_xxyy_0[i] * pa_z[i] - ta1_y_yzzz_xxyy_1[i] * pc_z[i];

        ta1_y_yzzzz_xxyz_0[i] = ta1_y_zzzz_xxz_0[i] * fe_0 - ta1_y_zzzz_xxz_1[i] * fe_0 + ta_zzzz_xxyz_1[i] + ta1_y_zzzz_xxyz_0[i] * pa_y[i] -
                                ta1_y_zzzz_xxyz_1[i] * pc_y[i];

        ta1_y_yzzzz_xxzz_0[i] = ta_zzzz_xxzz_1[i] + ta1_y_zzzz_xxzz_0[i] * pa_y[i] - ta1_y_zzzz_xxzz_1[i] * pc_y[i];

        ta1_y_yzzzz_xyyy_0[i] =
            3.0 * ta1_y_yzz_xyyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xyyy_1[i] * fe_0 + ta1_y_yzzz_xyyy_0[i] * pa_z[i] - ta1_y_yzzz_xyyy_1[i] * pc_z[i];

        ta1_y_yzzzz_xyyz_0[i] = 2.0 * ta1_y_zzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xyz_1[i] * fe_0 + ta_zzzz_xyyz_1[i] +
                                ta1_y_zzzz_xyyz_0[i] * pa_y[i] - ta1_y_zzzz_xyyz_1[i] * pc_y[i];

        ta1_y_yzzzz_xyzz_0[i] = ta1_y_zzzz_xzz_0[i] * fe_0 - ta1_y_zzzz_xzz_1[i] * fe_0 + ta_zzzz_xyzz_1[i] + ta1_y_zzzz_xyzz_0[i] * pa_y[i] -
                                ta1_y_zzzz_xyzz_1[i] * pc_y[i];

        ta1_y_yzzzz_xzzz_0[i] = ta_zzzz_xzzz_1[i] + ta1_y_zzzz_xzzz_0[i] * pa_y[i] - ta1_y_zzzz_xzzz_1[i] * pc_y[i];

        ta1_y_yzzzz_yyyy_0[i] =
            3.0 * ta1_y_yzz_yyyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_yyyy_1[i] * fe_0 + ta1_y_yzzz_yyyy_0[i] * pa_z[i] - ta1_y_yzzz_yyyy_1[i] * pc_z[i];

        ta1_y_yzzzz_yyyz_0[i] = 3.0 * ta1_y_zzzz_yyz_0[i] * fe_0 - 3.0 * ta1_y_zzzz_yyz_1[i] * fe_0 + ta_zzzz_yyyz_1[i] +
                                ta1_y_zzzz_yyyz_0[i] * pa_y[i] - ta1_y_zzzz_yyyz_1[i] * pc_y[i];

        ta1_y_yzzzz_yyzz_0[i] = 2.0 * ta1_y_zzzz_yzz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_yzz_1[i] * fe_0 + ta_zzzz_yyzz_1[i] +
                                ta1_y_zzzz_yyzz_0[i] * pa_y[i] - ta1_y_zzzz_yyzz_1[i] * pc_y[i];

        ta1_y_yzzzz_yzzz_0[i] = ta1_y_zzzz_zzz_0[i] * fe_0 - ta1_y_zzzz_zzz_1[i] * fe_0 + ta_zzzz_yzzz_1[i] + ta1_y_zzzz_yzzz_0[i] * pa_y[i] -
                                ta1_y_zzzz_yzzz_1[i] * pc_y[i];

        ta1_y_yzzzz_zzzz_0[i] = ta_zzzz_zzzz_1[i] + ta1_y_zzzz_zzzz_0[i] * pa_y[i] - ta1_y_zzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 615-630 components of targeted buffer : HG

    auto ta1_y_zzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 615);

    auto ta1_y_zzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 616);

    auto ta1_y_zzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 617);

    auto ta1_y_zzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 618);

    auto ta1_y_zzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 619);

    auto ta1_y_zzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 620);

    auto ta1_y_zzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 621);

    auto ta1_y_zzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 622);

    auto ta1_y_zzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 623);

    auto ta1_y_zzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 624);

    auto ta1_y_zzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 625);

    auto ta1_y_zzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 626);

    auto ta1_y_zzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 627);

    auto ta1_y_zzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 628);

    auto ta1_y_zzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 629);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_zzz_xxxx_0,   \
                             ta1_y_zzz_xxxx_1,   \
                             ta1_y_zzz_xxxy_0,   \
                             ta1_y_zzz_xxxy_1,   \
                             ta1_y_zzz_xxxz_0,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxyy_0,   \
                             ta1_y_zzz_xxyy_1,   \
                             ta1_y_zzz_xxyz_0,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxzz_0,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xyyy_0,   \
                             ta1_y_zzz_xyyy_1,   \
                             ta1_y_zzz_xyyz_0,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyzz_0,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xzzz_0,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_yyyy_0,   \
                             ta1_y_zzz_yyyy_1,   \
                             ta1_y_zzz_yyyz_0,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyzz_0,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yzzz_0,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_zzzz_0,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta1_y_zzzz_xxx_0,   \
                             ta1_y_zzzz_xxx_1,   \
                             ta1_y_zzzz_xxxx_0,  \
                             ta1_y_zzzz_xxxx_1,  \
                             ta1_y_zzzz_xxxy_0,  \
                             ta1_y_zzzz_xxxy_1,  \
                             ta1_y_zzzz_xxxz_0,  \
                             ta1_y_zzzz_xxxz_1,  \
                             ta1_y_zzzz_xxy_0,   \
                             ta1_y_zzzz_xxy_1,   \
                             ta1_y_zzzz_xxyy_0,  \
                             ta1_y_zzzz_xxyy_1,  \
                             ta1_y_zzzz_xxyz_0,  \
                             ta1_y_zzzz_xxyz_1,  \
                             ta1_y_zzzz_xxz_0,   \
                             ta1_y_zzzz_xxz_1,   \
                             ta1_y_zzzz_xxzz_0,  \
                             ta1_y_zzzz_xxzz_1,  \
                             ta1_y_zzzz_xyy_0,   \
                             ta1_y_zzzz_xyy_1,   \
                             ta1_y_zzzz_xyyy_0,  \
                             ta1_y_zzzz_xyyy_1,  \
                             ta1_y_zzzz_xyyz_0,  \
                             ta1_y_zzzz_xyyz_1,  \
                             ta1_y_zzzz_xyz_0,   \
                             ta1_y_zzzz_xyz_1,   \
                             ta1_y_zzzz_xyzz_0,  \
                             ta1_y_zzzz_xyzz_1,  \
                             ta1_y_zzzz_xzz_0,   \
                             ta1_y_zzzz_xzz_1,   \
                             ta1_y_zzzz_xzzz_0,  \
                             ta1_y_zzzz_xzzz_1,  \
                             ta1_y_zzzz_yyy_0,   \
                             ta1_y_zzzz_yyy_1,   \
                             ta1_y_zzzz_yyyy_0,  \
                             ta1_y_zzzz_yyyy_1,  \
                             ta1_y_zzzz_yyyz_0,  \
                             ta1_y_zzzz_yyyz_1,  \
                             ta1_y_zzzz_yyz_0,   \
                             ta1_y_zzzz_yyz_1,   \
                             ta1_y_zzzz_yyzz_0,  \
                             ta1_y_zzzz_yyzz_1,  \
                             ta1_y_zzzz_yzz_0,   \
                             ta1_y_zzzz_yzz_1,   \
                             ta1_y_zzzz_yzzz_0,  \
                             ta1_y_zzzz_yzzz_1,  \
                             ta1_y_zzzz_zzz_0,   \
                             ta1_y_zzzz_zzz_1,   \
                             ta1_y_zzzz_zzzz_0,  \
                             ta1_y_zzzz_zzzz_1,  \
                             ta1_y_zzzzz_xxxx_0, \
                             ta1_y_zzzzz_xxxy_0, \
                             ta1_y_zzzzz_xxxz_0, \
                             ta1_y_zzzzz_xxyy_0, \
                             ta1_y_zzzzz_xxyz_0, \
                             ta1_y_zzzzz_xxzz_0, \
                             ta1_y_zzzzz_xyyy_0, \
                             ta1_y_zzzzz_xyyz_0, \
                             ta1_y_zzzzz_xyzz_0, \
                             ta1_y_zzzzz_xzzz_0, \
                             ta1_y_zzzzz_yyyy_0, \
                             ta1_y_zzzzz_yyyz_0, \
                             ta1_y_zzzzz_yyzz_0, \
                             ta1_y_zzzzz_yzzz_0, \
                             ta1_y_zzzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzz_xxxx_0[i] =
            4.0 * ta1_y_zzz_xxxx_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxx_1[i] * fe_0 + ta1_y_zzzz_xxxx_0[i] * pa_z[i] - ta1_y_zzzz_xxxx_1[i] * pc_z[i];

        ta1_y_zzzzz_xxxy_0[i] =
            4.0 * ta1_y_zzz_xxxy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxy_1[i] * fe_0 + ta1_y_zzzz_xxxy_0[i] * pa_z[i] - ta1_y_zzzz_xxxy_1[i] * pc_z[i];

        ta1_y_zzzzz_xxxz_0[i] = 4.0 * ta1_y_zzz_xxxz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxz_1[i] * fe_0 + ta1_y_zzzz_xxx_0[i] * fe_0 -
                                ta1_y_zzzz_xxx_1[i] * fe_0 + ta1_y_zzzz_xxxz_0[i] * pa_z[i] - ta1_y_zzzz_xxxz_1[i] * pc_z[i];

        ta1_y_zzzzz_xxyy_0[i] =
            4.0 * ta1_y_zzz_xxyy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxyy_1[i] * fe_0 + ta1_y_zzzz_xxyy_0[i] * pa_z[i] - ta1_y_zzzz_xxyy_1[i] * pc_z[i];

        ta1_y_zzzzz_xxyz_0[i] = 4.0 * ta1_y_zzz_xxyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxyz_1[i] * fe_0 + ta1_y_zzzz_xxy_0[i] * fe_0 -
                                ta1_y_zzzz_xxy_1[i] * fe_0 + ta1_y_zzzz_xxyz_0[i] * pa_z[i] - ta1_y_zzzz_xxyz_1[i] * pc_z[i];

        ta1_y_zzzzz_xxzz_0[i] = 4.0 * ta1_y_zzz_xxzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_zzzz_xxz_0[i] * fe_0 -
                                2.0 * ta1_y_zzzz_xxz_1[i] * fe_0 + ta1_y_zzzz_xxzz_0[i] * pa_z[i] - ta1_y_zzzz_xxzz_1[i] * pc_z[i];

        ta1_y_zzzzz_xyyy_0[i] =
            4.0 * ta1_y_zzz_xyyy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xyyy_1[i] * fe_0 + ta1_y_zzzz_xyyy_0[i] * pa_z[i] - ta1_y_zzzz_xyyy_1[i] * pc_z[i];

        ta1_y_zzzzz_xyyz_0[i] = 4.0 * ta1_y_zzz_xyyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xyyz_1[i] * fe_0 + ta1_y_zzzz_xyy_0[i] * fe_0 -
                                ta1_y_zzzz_xyy_1[i] * fe_0 + ta1_y_zzzz_xyyz_0[i] * pa_z[i] - ta1_y_zzzz_xyyz_1[i] * pc_z[i];

        ta1_y_zzzzz_xyzz_0[i] = 4.0 * ta1_y_zzz_xyzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xyzz_1[i] * fe_0 + 2.0 * ta1_y_zzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_y_zzzz_xyz_1[i] * fe_0 + ta1_y_zzzz_xyzz_0[i] * pa_z[i] - ta1_y_zzzz_xyzz_1[i] * pc_z[i];

        ta1_y_zzzzz_xzzz_0[i] = 4.0 * ta1_y_zzz_xzzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xzzz_1[i] * fe_0 + 3.0 * ta1_y_zzzz_xzz_0[i] * fe_0 -
                                3.0 * ta1_y_zzzz_xzz_1[i] * fe_0 + ta1_y_zzzz_xzzz_0[i] * pa_z[i] - ta1_y_zzzz_xzzz_1[i] * pc_z[i];

        ta1_y_zzzzz_yyyy_0[i] =
            4.0 * ta1_y_zzz_yyyy_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyyy_1[i] * fe_0 + ta1_y_zzzz_yyyy_0[i] * pa_z[i] - ta1_y_zzzz_yyyy_1[i] * pc_z[i];

        ta1_y_zzzzz_yyyz_0[i] = 4.0 * ta1_y_zzz_yyyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyyz_1[i] * fe_0 + ta1_y_zzzz_yyy_0[i] * fe_0 -
                                ta1_y_zzzz_yyy_1[i] * fe_0 + ta1_y_zzzz_yyyz_0[i] * pa_z[i] - ta1_y_zzzz_yyyz_1[i] * pc_z[i];

        ta1_y_zzzzz_yyzz_0[i] = 4.0 * ta1_y_zzz_yyzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyzz_1[i] * fe_0 + 2.0 * ta1_y_zzzz_yyz_0[i] * fe_0 -
                                2.0 * ta1_y_zzzz_yyz_1[i] * fe_0 + ta1_y_zzzz_yyzz_0[i] * pa_z[i] - ta1_y_zzzz_yyzz_1[i] * pc_z[i];

        ta1_y_zzzzz_yzzz_0[i] = 4.0 * ta1_y_zzz_yzzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yzzz_1[i] * fe_0 + 3.0 * ta1_y_zzzz_yzz_0[i] * fe_0 -
                                3.0 * ta1_y_zzzz_yzz_1[i] * fe_0 + ta1_y_zzzz_yzzz_0[i] * pa_z[i] - ta1_y_zzzz_yzzz_1[i] * pc_z[i];

        ta1_y_zzzzz_zzzz_0[i] = 4.0 * ta1_y_zzz_zzzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_zzzz_1[i] * fe_0 + 4.0 * ta1_y_zzzz_zzz_0[i] * fe_0 -
                                4.0 * ta1_y_zzzz_zzz_1[i] * fe_0 + ta1_y_zzzz_zzzz_0[i] * pa_z[i] - ta1_y_zzzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 630-645 components of targeted buffer : HG

    auto ta1_z_xxxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 630);

    auto ta1_z_xxxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 631);

    auto ta1_z_xxxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 632);

    auto ta1_z_xxxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 633);

    auto ta1_z_xxxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 634);

    auto ta1_z_xxxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 635);

    auto ta1_z_xxxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 636);

    auto ta1_z_xxxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 637);

    auto ta1_z_xxxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 638);

    auto ta1_z_xxxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 639);

    auto ta1_z_xxxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 640);

    auto ta1_z_xxxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 641);

    auto ta1_z_xxxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 642);

    auto ta1_z_xxxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 643);

    auto ta1_z_xxxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 644);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xxx_xxxx_0,   \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxy_0,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxz_0,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxyy_0,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyz_0,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxzz_0,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xyyy_0,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyz_0,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyzz_0,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xzzz_0,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_yyyy_0,   \
                             ta1_z_xxx_yyyy_1,   \
                             ta1_z_xxx_yyyz_0,   \
                             ta1_z_xxx_yyyz_1,   \
                             ta1_z_xxx_yyzz_0,   \
                             ta1_z_xxx_yyzz_1,   \
                             ta1_z_xxx_yzzz_0,   \
                             ta1_z_xxx_yzzz_1,   \
                             ta1_z_xxx_zzzz_0,   \
                             ta1_z_xxx_zzzz_1,   \
                             ta1_z_xxxx_xxx_0,   \
                             ta1_z_xxxx_xxx_1,   \
                             ta1_z_xxxx_xxxx_0,  \
                             ta1_z_xxxx_xxxx_1,  \
                             ta1_z_xxxx_xxxy_0,  \
                             ta1_z_xxxx_xxxy_1,  \
                             ta1_z_xxxx_xxxz_0,  \
                             ta1_z_xxxx_xxxz_1,  \
                             ta1_z_xxxx_xxy_0,   \
                             ta1_z_xxxx_xxy_1,   \
                             ta1_z_xxxx_xxyy_0,  \
                             ta1_z_xxxx_xxyy_1,  \
                             ta1_z_xxxx_xxyz_0,  \
                             ta1_z_xxxx_xxyz_1,  \
                             ta1_z_xxxx_xxz_0,   \
                             ta1_z_xxxx_xxz_1,   \
                             ta1_z_xxxx_xxzz_0,  \
                             ta1_z_xxxx_xxzz_1,  \
                             ta1_z_xxxx_xyy_0,   \
                             ta1_z_xxxx_xyy_1,   \
                             ta1_z_xxxx_xyyy_0,  \
                             ta1_z_xxxx_xyyy_1,  \
                             ta1_z_xxxx_xyyz_0,  \
                             ta1_z_xxxx_xyyz_1,  \
                             ta1_z_xxxx_xyz_0,   \
                             ta1_z_xxxx_xyz_1,   \
                             ta1_z_xxxx_xyzz_0,  \
                             ta1_z_xxxx_xyzz_1,  \
                             ta1_z_xxxx_xzz_0,   \
                             ta1_z_xxxx_xzz_1,   \
                             ta1_z_xxxx_xzzz_0,  \
                             ta1_z_xxxx_xzzz_1,  \
                             ta1_z_xxxx_yyy_0,   \
                             ta1_z_xxxx_yyy_1,   \
                             ta1_z_xxxx_yyyy_0,  \
                             ta1_z_xxxx_yyyy_1,  \
                             ta1_z_xxxx_yyyz_0,  \
                             ta1_z_xxxx_yyyz_1,  \
                             ta1_z_xxxx_yyz_0,   \
                             ta1_z_xxxx_yyz_1,   \
                             ta1_z_xxxx_yyzz_0,  \
                             ta1_z_xxxx_yyzz_1,  \
                             ta1_z_xxxx_yzz_0,   \
                             ta1_z_xxxx_yzz_1,   \
                             ta1_z_xxxx_yzzz_0,  \
                             ta1_z_xxxx_yzzz_1,  \
                             ta1_z_xxxx_zzz_0,   \
                             ta1_z_xxxx_zzz_1,   \
                             ta1_z_xxxx_zzzz_0,  \
                             ta1_z_xxxx_zzzz_1,  \
                             ta1_z_xxxxx_xxxx_0, \
                             ta1_z_xxxxx_xxxy_0, \
                             ta1_z_xxxxx_xxxz_0, \
                             ta1_z_xxxxx_xxyy_0, \
                             ta1_z_xxxxx_xxyz_0, \
                             ta1_z_xxxxx_xxzz_0, \
                             ta1_z_xxxxx_xyyy_0, \
                             ta1_z_xxxxx_xyyz_0, \
                             ta1_z_xxxxx_xyzz_0, \
                             ta1_z_xxxxx_xzzz_0, \
                             ta1_z_xxxxx_yyyy_0, \
                             ta1_z_xxxxx_yyyz_0, \
                             ta1_z_xxxxx_yyzz_0, \
                             ta1_z_xxxxx_yzzz_0, \
                             ta1_z_xxxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxx_xxxx_0[i] = 4.0 * ta1_z_xxx_xxxx_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxxx_1[i] * fe_0 + 4.0 * ta1_z_xxxx_xxx_0[i] * fe_0 -
                                4.0 * ta1_z_xxxx_xxx_1[i] * fe_0 + ta1_z_xxxx_xxxx_0[i] * pa_x[i] - ta1_z_xxxx_xxxx_1[i] * pc_x[i];

        ta1_z_xxxxx_xxxy_0[i] = 4.0 * ta1_z_xxx_xxxy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xxxx_xxy_0[i] * fe_0 -
                                3.0 * ta1_z_xxxx_xxy_1[i] * fe_0 + ta1_z_xxxx_xxxy_0[i] * pa_x[i] - ta1_z_xxxx_xxxy_1[i] * pc_x[i];

        ta1_z_xxxxx_xxxz_0[i] = 4.0 * ta1_z_xxx_xxxz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xxxx_xxz_0[i] * fe_0 -
                                3.0 * ta1_z_xxxx_xxz_1[i] * fe_0 + ta1_z_xxxx_xxxz_0[i] * pa_x[i] - ta1_z_xxxx_xxxz_1[i] * pc_x[i];

        ta1_z_xxxxx_xxyy_0[i] = 4.0 * ta1_z_xxx_xxyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxxx_xyy_0[i] * fe_0 -
                                2.0 * ta1_z_xxxx_xyy_1[i] * fe_0 + ta1_z_xxxx_xxyy_0[i] * pa_x[i] - ta1_z_xxxx_xxyy_1[i] * pc_x[i];

        ta1_z_xxxxx_xxyz_0[i] = 4.0 * ta1_z_xxx_xxyz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxxx_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_xxxx_xyz_1[i] * fe_0 + ta1_z_xxxx_xxyz_0[i] * pa_x[i] - ta1_z_xxxx_xxyz_1[i] * pc_x[i];

        ta1_z_xxxxx_xxzz_0[i] = 4.0 * ta1_z_xxx_xxzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xxxx_xzz_0[i] * fe_0 -
                                2.0 * ta1_z_xxxx_xzz_1[i] * fe_0 + ta1_z_xxxx_xxzz_0[i] * pa_x[i] - ta1_z_xxxx_xxzz_1[i] * pc_x[i];

        ta1_z_xxxxx_xyyy_0[i] = 4.0 * ta1_z_xxx_xyyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyyy_1[i] * fe_0 + ta1_z_xxxx_yyy_0[i] * fe_0 -
                                ta1_z_xxxx_yyy_1[i] * fe_0 + ta1_z_xxxx_xyyy_0[i] * pa_x[i] - ta1_z_xxxx_xyyy_1[i] * pc_x[i];

        ta1_z_xxxxx_xyyz_0[i] = 4.0 * ta1_z_xxx_xyyz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyyz_1[i] * fe_0 + ta1_z_xxxx_yyz_0[i] * fe_0 -
                                ta1_z_xxxx_yyz_1[i] * fe_0 + ta1_z_xxxx_xyyz_0[i] * pa_x[i] - ta1_z_xxxx_xyyz_1[i] * pc_x[i];

        ta1_z_xxxxx_xyzz_0[i] = 4.0 * ta1_z_xxx_xyzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyzz_1[i] * fe_0 + ta1_z_xxxx_yzz_0[i] * fe_0 -
                                ta1_z_xxxx_yzz_1[i] * fe_0 + ta1_z_xxxx_xyzz_0[i] * pa_x[i] - ta1_z_xxxx_xyzz_1[i] * pc_x[i];

        ta1_z_xxxxx_xzzz_0[i] = 4.0 * ta1_z_xxx_xzzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xzzz_1[i] * fe_0 + ta1_z_xxxx_zzz_0[i] * fe_0 -
                                ta1_z_xxxx_zzz_1[i] * fe_0 + ta1_z_xxxx_xzzz_0[i] * pa_x[i] - ta1_z_xxxx_xzzz_1[i] * pc_x[i];

        ta1_z_xxxxx_yyyy_0[i] =
            4.0 * ta1_z_xxx_yyyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_yyyy_1[i] * fe_0 + ta1_z_xxxx_yyyy_0[i] * pa_x[i] - ta1_z_xxxx_yyyy_1[i] * pc_x[i];

        ta1_z_xxxxx_yyyz_0[i] =
            4.0 * ta1_z_xxx_yyyz_0[i] * fe_0 - 4.0 * ta1_z_xxx_yyyz_1[i] * fe_0 + ta1_z_xxxx_yyyz_0[i] * pa_x[i] - ta1_z_xxxx_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxx_yyzz_0[i] =
            4.0 * ta1_z_xxx_yyzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_yyzz_1[i] * fe_0 + ta1_z_xxxx_yyzz_0[i] * pa_x[i] - ta1_z_xxxx_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxx_yzzz_0[i] =
            4.0 * ta1_z_xxx_yzzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_yzzz_1[i] * fe_0 + ta1_z_xxxx_yzzz_0[i] * pa_x[i] - ta1_z_xxxx_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxx_zzzz_0[i] =
            4.0 * ta1_z_xxx_zzzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_zzzz_1[i] * fe_0 + ta1_z_xxxx_zzzz_0[i] * pa_x[i] - ta1_z_xxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 645-660 components of targeted buffer : HG

    auto ta1_z_xxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 645);

    auto ta1_z_xxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 646);

    auto ta1_z_xxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 647);

    auto ta1_z_xxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 648);

    auto ta1_z_xxxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 649);

    auto ta1_z_xxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 650);

    auto ta1_z_xxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 651);

    auto ta1_z_xxxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 652);

    auto ta1_z_xxxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 653);

    auto ta1_z_xxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 654);

    auto ta1_z_xxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 655);

    auto ta1_z_xxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 656);

    auto ta1_z_xxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 657);

    auto ta1_z_xxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 658);

    auto ta1_z_xxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 659);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxxx_xxx_0,   \
                             ta1_z_xxxx_xxx_1,   \
                             ta1_z_xxxx_xxxx_0,  \
                             ta1_z_xxxx_xxxx_1,  \
                             ta1_z_xxxx_xxxy_0,  \
                             ta1_z_xxxx_xxxy_1,  \
                             ta1_z_xxxx_xxxz_0,  \
                             ta1_z_xxxx_xxxz_1,  \
                             ta1_z_xxxx_xxy_0,   \
                             ta1_z_xxxx_xxy_1,   \
                             ta1_z_xxxx_xxyy_0,  \
                             ta1_z_xxxx_xxyy_1,  \
                             ta1_z_xxxx_xxyz_0,  \
                             ta1_z_xxxx_xxyz_1,  \
                             ta1_z_xxxx_xxz_0,   \
                             ta1_z_xxxx_xxz_1,   \
                             ta1_z_xxxx_xxzz_0,  \
                             ta1_z_xxxx_xxzz_1,  \
                             ta1_z_xxxx_xyy_0,   \
                             ta1_z_xxxx_xyy_1,   \
                             ta1_z_xxxx_xyyy_0,  \
                             ta1_z_xxxx_xyyy_1,  \
                             ta1_z_xxxx_xyyz_0,  \
                             ta1_z_xxxx_xyyz_1,  \
                             ta1_z_xxxx_xyz_0,   \
                             ta1_z_xxxx_xyz_1,   \
                             ta1_z_xxxx_xyzz_0,  \
                             ta1_z_xxxx_xyzz_1,  \
                             ta1_z_xxxx_xzz_0,   \
                             ta1_z_xxxx_xzz_1,   \
                             ta1_z_xxxx_xzzz_0,  \
                             ta1_z_xxxx_xzzz_1,  \
                             ta1_z_xxxx_zzzz_0,  \
                             ta1_z_xxxx_zzzz_1,  \
                             ta1_z_xxxxy_xxxx_0, \
                             ta1_z_xxxxy_xxxy_0, \
                             ta1_z_xxxxy_xxxz_0, \
                             ta1_z_xxxxy_xxyy_0, \
                             ta1_z_xxxxy_xxyz_0, \
                             ta1_z_xxxxy_xxzz_0, \
                             ta1_z_xxxxy_xyyy_0, \
                             ta1_z_xxxxy_xyyz_0, \
                             ta1_z_xxxxy_xyzz_0, \
                             ta1_z_xxxxy_xzzz_0, \
                             ta1_z_xxxxy_yyyy_0, \
                             ta1_z_xxxxy_yyyz_0, \
                             ta1_z_xxxxy_yyzz_0, \
                             ta1_z_xxxxy_yzzz_0, \
                             ta1_z_xxxxy_zzzz_0, \
                             ta1_z_xxxy_yyyy_0,  \
                             ta1_z_xxxy_yyyy_1,  \
                             ta1_z_xxxy_yyyz_0,  \
                             ta1_z_xxxy_yyyz_1,  \
                             ta1_z_xxxy_yyzz_0,  \
                             ta1_z_xxxy_yyzz_1,  \
                             ta1_z_xxxy_yzzz_0,  \
                             ta1_z_xxxy_yzzz_1,  \
                             ta1_z_xxy_yyyy_0,   \
                             ta1_z_xxy_yyyy_1,   \
                             ta1_z_xxy_yyyz_0,   \
                             ta1_z_xxy_yyyz_1,   \
                             ta1_z_xxy_yyzz_0,   \
                             ta1_z_xxy_yyzz_1,   \
                             ta1_z_xxy_yzzz_0,   \
                             ta1_z_xxy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxy_xxxx_0[i] = ta1_z_xxxx_xxxx_0[i] * pa_y[i] - ta1_z_xxxx_xxxx_1[i] * pc_y[i];

        ta1_z_xxxxy_xxxy_0[i] =
            ta1_z_xxxx_xxx_0[i] * fe_0 - ta1_z_xxxx_xxx_1[i] * fe_0 + ta1_z_xxxx_xxxy_0[i] * pa_y[i] - ta1_z_xxxx_xxxy_1[i] * pc_y[i];

        ta1_z_xxxxy_xxxz_0[i] = ta1_z_xxxx_xxxz_0[i] * pa_y[i] - ta1_z_xxxx_xxxz_1[i] * pc_y[i];

        ta1_z_xxxxy_xxyy_0[i] =
            2.0 * ta1_z_xxxx_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xxy_1[i] * fe_0 + ta1_z_xxxx_xxyy_0[i] * pa_y[i] - ta1_z_xxxx_xxyy_1[i] * pc_y[i];

        ta1_z_xxxxy_xxyz_0[i] =
            ta1_z_xxxx_xxz_0[i] * fe_0 - ta1_z_xxxx_xxz_1[i] * fe_0 + ta1_z_xxxx_xxyz_0[i] * pa_y[i] - ta1_z_xxxx_xxyz_1[i] * pc_y[i];

        ta1_z_xxxxy_xxzz_0[i] = ta1_z_xxxx_xxzz_0[i] * pa_y[i] - ta1_z_xxxx_xxzz_1[i] * pc_y[i];

        ta1_z_xxxxy_xyyy_0[i] =
            3.0 * ta1_z_xxxx_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxxx_xyy_1[i] * fe_0 + ta1_z_xxxx_xyyy_0[i] * pa_y[i] - ta1_z_xxxx_xyyy_1[i] * pc_y[i];

        ta1_z_xxxxy_xyyz_0[i] =
            2.0 * ta1_z_xxxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xyz_1[i] * fe_0 + ta1_z_xxxx_xyyz_0[i] * pa_y[i] - ta1_z_xxxx_xyyz_1[i] * pc_y[i];

        ta1_z_xxxxy_xyzz_0[i] =
            ta1_z_xxxx_xzz_0[i] * fe_0 - ta1_z_xxxx_xzz_1[i] * fe_0 + ta1_z_xxxx_xyzz_0[i] * pa_y[i] - ta1_z_xxxx_xyzz_1[i] * pc_y[i];

        ta1_z_xxxxy_xzzz_0[i] = ta1_z_xxxx_xzzz_0[i] * pa_y[i] - ta1_z_xxxx_xzzz_1[i] * pc_y[i];

        ta1_z_xxxxy_yyyy_0[i] =
            3.0 * ta1_z_xxy_yyyy_0[i] * fe_0 - 3.0 * ta1_z_xxy_yyyy_1[i] * fe_0 + ta1_z_xxxy_yyyy_0[i] * pa_x[i] - ta1_z_xxxy_yyyy_1[i] * pc_x[i];

        ta1_z_xxxxy_yyyz_0[i] =
            3.0 * ta1_z_xxy_yyyz_0[i] * fe_0 - 3.0 * ta1_z_xxy_yyyz_1[i] * fe_0 + ta1_z_xxxy_yyyz_0[i] * pa_x[i] - ta1_z_xxxy_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxy_yyzz_0[i] =
            3.0 * ta1_z_xxy_yyzz_0[i] * fe_0 - 3.0 * ta1_z_xxy_yyzz_1[i] * fe_0 + ta1_z_xxxy_yyzz_0[i] * pa_x[i] - ta1_z_xxxy_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxy_yzzz_0[i] =
            3.0 * ta1_z_xxy_yzzz_0[i] * fe_0 - 3.0 * ta1_z_xxy_yzzz_1[i] * fe_0 + ta1_z_xxxy_yzzz_0[i] * pa_x[i] - ta1_z_xxxy_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxy_zzzz_0[i] = ta1_z_xxxx_zzzz_0[i] * pa_y[i] - ta1_z_xxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 660-675 components of targeted buffer : HG

    auto ta1_z_xxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 660);

    auto ta1_z_xxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 661);

    auto ta1_z_xxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 662);

    auto ta1_z_xxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 663);

    auto ta1_z_xxxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 664);

    auto ta1_z_xxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 665);

    auto ta1_z_xxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 666);

    auto ta1_z_xxxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 667);

    auto ta1_z_xxxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 668);

    auto ta1_z_xxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 669);

    auto ta1_z_xxxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 670);

    auto ta1_z_xxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 671);

    auto ta1_z_xxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 672);

    auto ta1_z_xxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 673);

    auto ta1_z_xxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 674);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxxx_xxx_0,   \
                             ta1_z_xxxx_xxx_1,   \
                             ta1_z_xxxx_xxxx_0,  \
                             ta1_z_xxxx_xxxx_1,  \
                             ta1_z_xxxx_xxxy_0,  \
                             ta1_z_xxxx_xxxy_1,  \
                             ta1_z_xxxx_xxxz_0,  \
                             ta1_z_xxxx_xxxz_1,  \
                             ta1_z_xxxx_xxy_0,   \
                             ta1_z_xxxx_xxy_1,   \
                             ta1_z_xxxx_xxyy_0,  \
                             ta1_z_xxxx_xxyy_1,  \
                             ta1_z_xxxx_xxyz_0,  \
                             ta1_z_xxxx_xxyz_1,  \
                             ta1_z_xxxx_xxz_0,   \
                             ta1_z_xxxx_xxz_1,   \
                             ta1_z_xxxx_xxzz_0,  \
                             ta1_z_xxxx_xxzz_1,  \
                             ta1_z_xxxx_xyy_0,   \
                             ta1_z_xxxx_xyy_1,   \
                             ta1_z_xxxx_xyyy_0,  \
                             ta1_z_xxxx_xyyy_1,  \
                             ta1_z_xxxx_xyyz_0,  \
                             ta1_z_xxxx_xyyz_1,  \
                             ta1_z_xxxx_xyz_0,   \
                             ta1_z_xxxx_xyz_1,   \
                             ta1_z_xxxx_xyzz_0,  \
                             ta1_z_xxxx_xyzz_1,  \
                             ta1_z_xxxx_xzz_0,   \
                             ta1_z_xxxx_xzz_1,   \
                             ta1_z_xxxx_xzzz_0,  \
                             ta1_z_xxxx_xzzz_1,  \
                             ta1_z_xxxx_yyyy_0,  \
                             ta1_z_xxxx_yyyy_1,  \
                             ta1_z_xxxxz_xxxx_0, \
                             ta1_z_xxxxz_xxxy_0, \
                             ta1_z_xxxxz_xxxz_0, \
                             ta1_z_xxxxz_xxyy_0, \
                             ta1_z_xxxxz_xxyz_0, \
                             ta1_z_xxxxz_xxzz_0, \
                             ta1_z_xxxxz_xyyy_0, \
                             ta1_z_xxxxz_xyyz_0, \
                             ta1_z_xxxxz_xyzz_0, \
                             ta1_z_xxxxz_xzzz_0, \
                             ta1_z_xxxxz_yyyy_0, \
                             ta1_z_xxxxz_yyyz_0, \
                             ta1_z_xxxxz_yyzz_0, \
                             ta1_z_xxxxz_yzzz_0, \
                             ta1_z_xxxxz_zzzz_0, \
                             ta1_z_xxxz_yyyz_0,  \
                             ta1_z_xxxz_yyyz_1,  \
                             ta1_z_xxxz_yyzz_0,  \
                             ta1_z_xxxz_yyzz_1,  \
                             ta1_z_xxxz_yzzz_0,  \
                             ta1_z_xxxz_yzzz_1,  \
                             ta1_z_xxxz_zzzz_0,  \
                             ta1_z_xxxz_zzzz_1,  \
                             ta1_z_xxz_yyyz_0,   \
                             ta1_z_xxz_yyyz_1,   \
                             ta1_z_xxz_yyzz_0,   \
                             ta1_z_xxz_yyzz_1,   \
                             ta1_z_xxz_yzzz_0,   \
                             ta1_z_xxz_yzzz_1,   \
                             ta1_z_xxz_zzzz_0,   \
                             ta1_z_xxz_zzzz_1,   \
                             ta_xxxx_xxxx_1,     \
                             ta_xxxx_xxxy_1,     \
                             ta_xxxx_xxxz_1,     \
                             ta_xxxx_xxyy_1,     \
                             ta_xxxx_xxyz_1,     \
                             ta_xxxx_xxzz_1,     \
                             ta_xxxx_xyyy_1,     \
                             ta_xxxx_xyyz_1,     \
                             ta_xxxx_xyzz_1,     \
                             ta_xxxx_xzzz_1,     \
                             ta_xxxx_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxz_xxxx_0[i] = ta_xxxx_xxxx_1[i] + ta1_z_xxxx_xxxx_0[i] * pa_z[i] - ta1_z_xxxx_xxxx_1[i] * pc_z[i];

        ta1_z_xxxxz_xxxy_0[i] = ta_xxxx_xxxy_1[i] + ta1_z_xxxx_xxxy_0[i] * pa_z[i] - ta1_z_xxxx_xxxy_1[i] * pc_z[i];

        ta1_z_xxxxz_xxxz_0[i] = ta1_z_xxxx_xxx_0[i] * fe_0 - ta1_z_xxxx_xxx_1[i] * fe_0 + ta_xxxx_xxxz_1[i] + ta1_z_xxxx_xxxz_0[i] * pa_z[i] -
                                ta1_z_xxxx_xxxz_1[i] * pc_z[i];

        ta1_z_xxxxz_xxyy_0[i] = ta_xxxx_xxyy_1[i] + ta1_z_xxxx_xxyy_0[i] * pa_z[i] - ta1_z_xxxx_xxyy_1[i] * pc_z[i];

        ta1_z_xxxxz_xxyz_0[i] = ta1_z_xxxx_xxy_0[i] * fe_0 - ta1_z_xxxx_xxy_1[i] * fe_0 + ta_xxxx_xxyz_1[i] + ta1_z_xxxx_xxyz_0[i] * pa_z[i] -
                                ta1_z_xxxx_xxyz_1[i] * pc_z[i];

        ta1_z_xxxxz_xxzz_0[i] = 2.0 * ta1_z_xxxx_xxz_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xxz_1[i] * fe_0 + ta_xxxx_xxzz_1[i] +
                                ta1_z_xxxx_xxzz_0[i] * pa_z[i] - ta1_z_xxxx_xxzz_1[i] * pc_z[i];

        ta1_z_xxxxz_xyyy_0[i] = ta_xxxx_xyyy_1[i] + ta1_z_xxxx_xyyy_0[i] * pa_z[i] - ta1_z_xxxx_xyyy_1[i] * pc_z[i];

        ta1_z_xxxxz_xyyz_0[i] = ta1_z_xxxx_xyy_0[i] * fe_0 - ta1_z_xxxx_xyy_1[i] * fe_0 + ta_xxxx_xyyz_1[i] + ta1_z_xxxx_xyyz_0[i] * pa_z[i] -
                                ta1_z_xxxx_xyyz_1[i] * pc_z[i];

        ta1_z_xxxxz_xyzz_0[i] = 2.0 * ta1_z_xxxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xyz_1[i] * fe_0 + ta_xxxx_xyzz_1[i] +
                                ta1_z_xxxx_xyzz_0[i] * pa_z[i] - ta1_z_xxxx_xyzz_1[i] * pc_z[i];

        ta1_z_xxxxz_xzzz_0[i] = 3.0 * ta1_z_xxxx_xzz_0[i] * fe_0 - 3.0 * ta1_z_xxxx_xzz_1[i] * fe_0 + ta_xxxx_xzzz_1[i] +
                                ta1_z_xxxx_xzzz_0[i] * pa_z[i] - ta1_z_xxxx_xzzz_1[i] * pc_z[i];

        ta1_z_xxxxz_yyyy_0[i] = ta_xxxx_yyyy_1[i] + ta1_z_xxxx_yyyy_0[i] * pa_z[i] - ta1_z_xxxx_yyyy_1[i] * pc_z[i];

        ta1_z_xxxxz_yyyz_0[i] =
            3.0 * ta1_z_xxz_yyyz_0[i] * fe_0 - 3.0 * ta1_z_xxz_yyyz_1[i] * fe_0 + ta1_z_xxxz_yyyz_0[i] * pa_x[i] - ta1_z_xxxz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxz_yyzz_0[i] =
            3.0 * ta1_z_xxz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_xxz_yyzz_1[i] * fe_0 + ta1_z_xxxz_yyzz_0[i] * pa_x[i] - ta1_z_xxxz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxz_yzzz_0[i] =
            3.0 * ta1_z_xxz_yzzz_0[i] * fe_0 - 3.0 * ta1_z_xxz_yzzz_1[i] * fe_0 + ta1_z_xxxz_yzzz_0[i] * pa_x[i] - ta1_z_xxxz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxz_zzzz_0[i] =
            3.0 * ta1_z_xxz_zzzz_0[i] * fe_0 - 3.0 * ta1_z_xxz_zzzz_1[i] * fe_0 + ta1_z_xxxz_zzzz_0[i] * pa_x[i] - ta1_z_xxxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 675-690 components of targeted buffer : HG

    auto ta1_z_xxxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 675);

    auto ta1_z_xxxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 676);

    auto ta1_z_xxxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 677);

    auto ta1_z_xxxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 678);

    auto ta1_z_xxxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 679);

    auto ta1_z_xxxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 680);

    auto ta1_z_xxxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 681);

    auto ta1_z_xxxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 682);

    auto ta1_z_xxxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 683);

    auto ta1_z_xxxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 684);

    auto ta1_z_xxxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 685);

    auto ta1_z_xxxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 686);

    auto ta1_z_xxxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 687);

    auto ta1_z_xxxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 688);

    auto ta1_z_xxxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 689);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxx_xxxx_0,   \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxz_0,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxzz_0,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xzzz_0,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxxy_xxxx_0,  \
                             ta1_z_xxxy_xxxx_1,  \
                             ta1_z_xxxy_xxxz_0,  \
                             ta1_z_xxxy_xxxz_1,  \
                             ta1_z_xxxy_xxzz_0,  \
                             ta1_z_xxxy_xxzz_1,  \
                             ta1_z_xxxy_xzzz_0,  \
                             ta1_z_xxxy_xzzz_1,  \
                             ta1_z_xxxyy_xxxx_0, \
                             ta1_z_xxxyy_xxxy_0, \
                             ta1_z_xxxyy_xxxz_0, \
                             ta1_z_xxxyy_xxyy_0, \
                             ta1_z_xxxyy_xxyz_0, \
                             ta1_z_xxxyy_xxzz_0, \
                             ta1_z_xxxyy_xyyy_0, \
                             ta1_z_xxxyy_xyyz_0, \
                             ta1_z_xxxyy_xyzz_0, \
                             ta1_z_xxxyy_xzzz_0, \
                             ta1_z_xxxyy_yyyy_0, \
                             ta1_z_xxxyy_yyyz_0, \
                             ta1_z_xxxyy_yyzz_0, \
                             ta1_z_xxxyy_yzzz_0, \
                             ta1_z_xxxyy_zzzz_0, \
                             ta1_z_xxyy_xxxy_0,  \
                             ta1_z_xxyy_xxxy_1,  \
                             ta1_z_xxyy_xxy_0,   \
                             ta1_z_xxyy_xxy_1,   \
                             ta1_z_xxyy_xxyy_0,  \
                             ta1_z_xxyy_xxyy_1,  \
                             ta1_z_xxyy_xxyz_0,  \
                             ta1_z_xxyy_xxyz_1,  \
                             ta1_z_xxyy_xyy_0,   \
                             ta1_z_xxyy_xyy_1,   \
                             ta1_z_xxyy_xyyy_0,  \
                             ta1_z_xxyy_xyyy_1,  \
                             ta1_z_xxyy_xyyz_0,  \
                             ta1_z_xxyy_xyyz_1,  \
                             ta1_z_xxyy_xyz_0,   \
                             ta1_z_xxyy_xyz_1,   \
                             ta1_z_xxyy_xyzz_0,  \
                             ta1_z_xxyy_xyzz_1,  \
                             ta1_z_xxyy_yyy_0,   \
                             ta1_z_xxyy_yyy_1,   \
                             ta1_z_xxyy_yyyy_0,  \
                             ta1_z_xxyy_yyyy_1,  \
                             ta1_z_xxyy_yyyz_0,  \
                             ta1_z_xxyy_yyyz_1,  \
                             ta1_z_xxyy_yyz_0,   \
                             ta1_z_xxyy_yyz_1,   \
                             ta1_z_xxyy_yyzz_0,  \
                             ta1_z_xxyy_yyzz_1,  \
                             ta1_z_xxyy_yzz_0,   \
                             ta1_z_xxyy_yzz_1,   \
                             ta1_z_xxyy_yzzz_0,  \
                             ta1_z_xxyy_yzzz_1,  \
                             ta1_z_xxyy_zzzz_0,  \
                             ta1_z_xxyy_zzzz_1,  \
                             ta1_z_xyy_xxxy_0,   \
                             ta1_z_xyy_xxxy_1,   \
                             ta1_z_xyy_xxyy_0,   \
                             ta1_z_xyy_xxyy_1,   \
                             ta1_z_xyy_xxyz_0,   \
                             ta1_z_xyy_xxyz_1,   \
                             ta1_z_xyy_xyyy_0,   \
                             ta1_z_xyy_xyyy_1,   \
                             ta1_z_xyy_xyyz_0,   \
                             ta1_z_xyy_xyyz_1,   \
                             ta1_z_xyy_xyzz_0,   \
                             ta1_z_xyy_xyzz_1,   \
                             ta1_z_xyy_yyyy_0,   \
                             ta1_z_xyy_yyyy_1,   \
                             ta1_z_xyy_yyyz_0,   \
                             ta1_z_xyy_yyyz_1,   \
                             ta1_z_xyy_yyzz_0,   \
                             ta1_z_xyy_yyzz_1,   \
                             ta1_z_xyy_yzzz_0,   \
                             ta1_z_xyy_yzzz_1,   \
                             ta1_z_xyy_zzzz_0,   \
                             ta1_z_xyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyy_xxxx_0[i] =
            ta1_z_xxx_xxxx_0[i] * fe_0 - ta1_z_xxx_xxxx_1[i] * fe_0 + ta1_z_xxxy_xxxx_0[i] * pa_y[i] - ta1_z_xxxy_xxxx_1[i] * pc_y[i];

        ta1_z_xxxyy_xxxy_0[i] = 2.0 * ta1_z_xyy_xxxy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xxyy_xxy_0[i] * fe_0 -
                                3.0 * ta1_z_xxyy_xxy_1[i] * fe_0 + ta1_z_xxyy_xxxy_0[i] * pa_x[i] - ta1_z_xxyy_xxxy_1[i] * pc_x[i];

        ta1_z_xxxyy_xxxz_0[i] =
            ta1_z_xxx_xxxz_0[i] * fe_0 - ta1_z_xxx_xxxz_1[i] * fe_0 + ta1_z_xxxy_xxxz_0[i] * pa_y[i] - ta1_z_xxxy_xxxz_1[i] * pc_y[i];

        ta1_z_xxxyy_xxyy_0[i] = 2.0 * ta1_z_xyy_xxyy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxyy_xyy_0[i] * fe_0 -
                                2.0 * ta1_z_xxyy_xyy_1[i] * fe_0 + ta1_z_xxyy_xxyy_0[i] * pa_x[i] - ta1_z_xxyy_xxyy_1[i] * pc_x[i];

        ta1_z_xxxyy_xxyz_0[i] = 2.0 * ta1_z_xyy_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xyy_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_xxyy_xyz_1[i] * fe_0 + ta1_z_xxyy_xxyz_0[i] * pa_x[i] - ta1_z_xxyy_xxyz_1[i] * pc_x[i];

        ta1_z_xxxyy_xxzz_0[i] =
            ta1_z_xxx_xxzz_0[i] * fe_0 - ta1_z_xxx_xxzz_1[i] * fe_0 + ta1_z_xxxy_xxzz_0[i] * pa_y[i] - ta1_z_xxxy_xxzz_1[i] * pc_y[i];

        ta1_z_xxxyy_xyyy_0[i] = 2.0 * ta1_z_xyy_xyyy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xyyy_1[i] * fe_0 + ta1_z_xxyy_yyy_0[i] * fe_0 -
                                ta1_z_xxyy_yyy_1[i] * fe_0 + ta1_z_xxyy_xyyy_0[i] * pa_x[i] - ta1_z_xxyy_xyyy_1[i] * pc_x[i];

        ta1_z_xxxyy_xyyz_0[i] = 2.0 * ta1_z_xyy_xyyz_0[i] * fe_0 - 2.0 * ta1_z_xyy_xyyz_1[i] * fe_0 + ta1_z_xxyy_yyz_0[i] * fe_0 -
                                ta1_z_xxyy_yyz_1[i] * fe_0 + ta1_z_xxyy_xyyz_0[i] * pa_x[i] - ta1_z_xxyy_xyyz_1[i] * pc_x[i];

        ta1_z_xxxyy_xyzz_0[i] = 2.0 * ta1_z_xyy_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xyy_xyzz_1[i] * fe_0 + ta1_z_xxyy_yzz_0[i] * fe_0 -
                                ta1_z_xxyy_yzz_1[i] * fe_0 + ta1_z_xxyy_xyzz_0[i] * pa_x[i] - ta1_z_xxyy_xyzz_1[i] * pc_x[i];

        ta1_z_xxxyy_xzzz_0[i] =
            ta1_z_xxx_xzzz_0[i] * fe_0 - ta1_z_xxx_xzzz_1[i] * fe_0 + ta1_z_xxxy_xzzz_0[i] * pa_y[i] - ta1_z_xxxy_xzzz_1[i] * pc_y[i];

        ta1_z_xxxyy_yyyy_0[i] =
            2.0 * ta1_z_xyy_yyyy_0[i] * fe_0 - 2.0 * ta1_z_xyy_yyyy_1[i] * fe_0 + ta1_z_xxyy_yyyy_0[i] * pa_x[i] - ta1_z_xxyy_yyyy_1[i] * pc_x[i];

        ta1_z_xxxyy_yyyz_0[i] =
            2.0 * ta1_z_xyy_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xyy_yyyz_1[i] * fe_0 + ta1_z_xxyy_yyyz_0[i] * pa_x[i] - ta1_z_xxyy_yyyz_1[i] * pc_x[i];

        ta1_z_xxxyy_yyzz_0[i] =
            2.0 * ta1_z_xyy_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xyy_yyzz_1[i] * fe_0 + ta1_z_xxyy_yyzz_0[i] * pa_x[i] - ta1_z_xxyy_yyzz_1[i] * pc_x[i];

        ta1_z_xxxyy_yzzz_0[i] =
            2.0 * ta1_z_xyy_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xyy_yzzz_1[i] * fe_0 + ta1_z_xxyy_yzzz_0[i] * pa_x[i] - ta1_z_xxyy_yzzz_1[i] * pc_x[i];

        ta1_z_xxxyy_zzzz_0[i] =
            2.0 * ta1_z_xyy_zzzz_0[i] * fe_0 - 2.0 * ta1_z_xyy_zzzz_1[i] * fe_0 + ta1_z_xxyy_zzzz_0[i] * pa_x[i] - ta1_z_xxyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 690-705 components of targeted buffer : HG

    auto ta1_z_xxxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 690);

    auto ta1_z_xxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 691);

    auto ta1_z_xxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 692);

    auto ta1_z_xxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 693);

    auto ta1_z_xxxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 694);

    auto ta1_z_xxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 695);

    auto ta1_z_xxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 696);

    auto ta1_z_xxxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 697);

    auto ta1_z_xxxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 698);

    auto ta1_z_xxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 699);

    auto ta1_z_xxxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 700);

    auto ta1_z_xxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 701);

    auto ta1_z_xxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 702);

    auto ta1_z_xxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 703);

    auto ta1_z_xxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 704);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxxy_xxxy_0,  \
                             ta1_z_xxxy_xxxy_1,  \
                             ta1_z_xxxy_xxyy_0,  \
                             ta1_z_xxxy_xxyy_1,  \
                             ta1_z_xxxy_xyyy_0,  \
                             ta1_z_xxxy_xyyy_1,  \
                             ta1_z_xxxy_yyyy_0,  \
                             ta1_z_xxxy_yyyy_1,  \
                             ta1_z_xxxyz_xxxx_0, \
                             ta1_z_xxxyz_xxxy_0, \
                             ta1_z_xxxyz_xxxz_0, \
                             ta1_z_xxxyz_xxyy_0, \
                             ta1_z_xxxyz_xxyz_0, \
                             ta1_z_xxxyz_xxzz_0, \
                             ta1_z_xxxyz_xyyy_0, \
                             ta1_z_xxxyz_xyyz_0, \
                             ta1_z_xxxyz_xyzz_0, \
                             ta1_z_xxxyz_xzzz_0, \
                             ta1_z_xxxyz_yyyy_0, \
                             ta1_z_xxxyz_yyyz_0, \
                             ta1_z_xxxyz_yyzz_0, \
                             ta1_z_xxxyz_yzzz_0, \
                             ta1_z_xxxyz_zzzz_0, \
                             ta1_z_xxxz_xxxx_0,  \
                             ta1_z_xxxz_xxxx_1,  \
                             ta1_z_xxxz_xxxz_0,  \
                             ta1_z_xxxz_xxxz_1,  \
                             ta1_z_xxxz_xxyz_0,  \
                             ta1_z_xxxz_xxyz_1,  \
                             ta1_z_xxxz_xxz_0,   \
                             ta1_z_xxxz_xxz_1,   \
                             ta1_z_xxxz_xxzz_0,  \
                             ta1_z_xxxz_xxzz_1,  \
                             ta1_z_xxxz_xyyz_0,  \
                             ta1_z_xxxz_xyyz_1,  \
                             ta1_z_xxxz_xyz_0,   \
                             ta1_z_xxxz_xyz_1,   \
                             ta1_z_xxxz_xyzz_0,  \
                             ta1_z_xxxz_xyzz_1,  \
                             ta1_z_xxxz_xzz_0,   \
                             ta1_z_xxxz_xzz_1,   \
                             ta1_z_xxxz_xzzz_0,  \
                             ta1_z_xxxz_xzzz_1,  \
                             ta1_z_xxxz_zzzz_0,  \
                             ta1_z_xxxz_zzzz_1,  \
                             ta1_z_xxyz_yyyz_0,  \
                             ta1_z_xxyz_yyyz_1,  \
                             ta1_z_xxyz_yyzz_0,  \
                             ta1_z_xxyz_yyzz_1,  \
                             ta1_z_xxyz_yzzz_0,  \
                             ta1_z_xxyz_yzzz_1,  \
                             ta1_z_xyz_yyyz_0,   \
                             ta1_z_xyz_yyyz_1,   \
                             ta1_z_xyz_yyzz_0,   \
                             ta1_z_xyz_yyzz_1,   \
                             ta1_z_xyz_yzzz_0,   \
                             ta1_z_xyz_yzzz_1,   \
                             ta_xxxy_xxxy_1,     \
                             ta_xxxy_xxyy_1,     \
                             ta_xxxy_xyyy_1,     \
                             ta_xxxy_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyz_xxxx_0[i] = ta1_z_xxxz_xxxx_0[i] * pa_y[i] - ta1_z_xxxz_xxxx_1[i] * pc_y[i];

        ta1_z_xxxyz_xxxy_0[i] = ta_xxxy_xxxy_1[i] + ta1_z_xxxy_xxxy_0[i] * pa_z[i] - ta1_z_xxxy_xxxy_1[i] * pc_z[i];

        ta1_z_xxxyz_xxxz_0[i] = ta1_z_xxxz_xxxz_0[i] * pa_y[i] - ta1_z_xxxz_xxxz_1[i] * pc_y[i];

        ta1_z_xxxyz_xxyy_0[i] = ta_xxxy_xxyy_1[i] + ta1_z_xxxy_xxyy_0[i] * pa_z[i] - ta1_z_xxxy_xxyy_1[i] * pc_z[i];

        ta1_z_xxxyz_xxyz_0[i] =
            ta1_z_xxxz_xxz_0[i] * fe_0 - ta1_z_xxxz_xxz_1[i] * fe_0 + ta1_z_xxxz_xxyz_0[i] * pa_y[i] - ta1_z_xxxz_xxyz_1[i] * pc_y[i];

        ta1_z_xxxyz_xxzz_0[i] = ta1_z_xxxz_xxzz_0[i] * pa_y[i] - ta1_z_xxxz_xxzz_1[i] * pc_y[i];

        ta1_z_xxxyz_xyyy_0[i] = ta_xxxy_xyyy_1[i] + ta1_z_xxxy_xyyy_0[i] * pa_z[i] - ta1_z_xxxy_xyyy_1[i] * pc_z[i];

        ta1_z_xxxyz_xyyz_0[i] =
            2.0 * ta1_z_xxxz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xyz_1[i] * fe_0 + ta1_z_xxxz_xyyz_0[i] * pa_y[i] - ta1_z_xxxz_xyyz_1[i] * pc_y[i];

        ta1_z_xxxyz_xyzz_0[i] =
            ta1_z_xxxz_xzz_0[i] * fe_0 - ta1_z_xxxz_xzz_1[i] * fe_0 + ta1_z_xxxz_xyzz_0[i] * pa_y[i] - ta1_z_xxxz_xyzz_1[i] * pc_y[i];

        ta1_z_xxxyz_xzzz_0[i] = ta1_z_xxxz_xzzz_0[i] * pa_y[i] - ta1_z_xxxz_xzzz_1[i] * pc_y[i];

        ta1_z_xxxyz_yyyy_0[i] = ta_xxxy_yyyy_1[i] + ta1_z_xxxy_yyyy_0[i] * pa_z[i] - ta1_z_xxxy_yyyy_1[i] * pc_z[i];

        ta1_z_xxxyz_yyyz_0[i] =
            2.0 * ta1_z_xyz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xyz_yyyz_1[i] * fe_0 + ta1_z_xxyz_yyyz_0[i] * pa_x[i] - ta1_z_xxyz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxyz_yyzz_0[i] =
            2.0 * ta1_z_xyz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xyz_yyzz_1[i] * fe_0 + ta1_z_xxyz_yyzz_0[i] * pa_x[i] - ta1_z_xxyz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxyz_yzzz_0[i] =
            2.0 * ta1_z_xyz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xyz_yzzz_1[i] * fe_0 + ta1_z_xxyz_yzzz_0[i] * pa_x[i] - ta1_z_xxyz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxyz_zzzz_0[i] = ta1_z_xxxz_zzzz_0[i] * pa_y[i] - ta1_z_xxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 705-720 components of targeted buffer : HG

    auto ta1_z_xxxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 705);

    auto ta1_z_xxxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 706);

    auto ta1_z_xxxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 707);

    auto ta1_z_xxxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 708);

    auto ta1_z_xxxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 709);

    auto ta1_z_xxxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 710);

    auto ta1_z_xxxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 711);

    auto ta1_z_xxxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 712);

    auto ta1_z_xxxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 713);

    auto ta1_z_xxxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 714);

    auto ta1_z_xxxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 715);

    auto ta1_z_xxxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 716);

    auto ta1_z_xxxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 717);

    auto ta1_z_xxxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 718);

    auto ta1_z_xxxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 719);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxx_xxxx_0,   \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxy_0,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxyy_0,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xyyy_0,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxxz_xxxx_0,  \
                             ta1_z_xxxz_xxxx_1,  \
                             ta1_z_xxxz_xxxy_0,  \
                             ta1_z_xxxz_xxxy_1,  \
                             ta1_z_xxxz_xxyy_0,  \
                             ta1_z_xxxz_xxyy_1,  \
                             ta1_z_xxxz_xyyy_0,  \
                             ta1_z_xxxz_xyyy_1,  \
                             ta1_z_xxxzz_xxxx_0, \
                             ta1_z_xxxzz_xxxy_0, \
                             ta1_z_xxxzz_xxxz_0, \
                             ta1_z_xxxzz_xxyy_0, \
                             ta1_z_xxxzz_xxyz_0, \
                             ta1_z_xxxzz_xxzz_0, \
                             ta1_z_xxxzz_xyyy_0, \
                             ta1_z_xxxzz_xyyz_0, \
                             ta1_z_xxxzz_xyzz_0, \
                             ta1_z_xxxzz_xzzz_0, \
                             ta1_z_xxxzz_yyyy_0, \
                             ta1_z_xxxzz_yyyz_0, \
                             ta1_z_xxxzz_yyzz_0, \
                             ta1_z_xxxzz_yzzz_0, \
                             ta1_z_xxxzz_zzzz_0, \
                             ta1_z_xxzz_xxxz_0,  \
                             ta1_z_xxzz_xxxz_1,  \
                             ta1_z_xxzz_xxyz_0,  \
                             ta1_z_xxzz_xxyz_1,  \
                             ta1_z_xxzz_xxz_0,   \
                             ta1_z_xxzz_xxz_1,   \
                             ta1_z_xxzz_xxzz_0,  \
                             ta1_z_xxzz_xxzz_1,  \
                             ta1_z_xxzz_xyyz_0,  \
                             ta1_z_xxzz_xyyz_1,  \
                             ta1_z_xxzz_xyz_0,   \
                             ta1_z_xxzz_xyz_1,   \
                             ta1_z_xxzz_xyzz_0,  \
                             ta1_z_xxzz_xyzz_1,  \
                             ta1_z_xxzz_xzz_0,   \
                             ta1_z_xxzz_xzz_1,   \
                             ta1_z_xxzz_xzzz_0,  \
                             ta1_z_xxzz_xzzz_1,  \
                             ta1_z_xxzz_yyyy_0,  \
                             ta1_z_xxzz_yyyy_1,  \
                             ta1_z_xxzz_yyyz_0,  \
                             ta1_z_xxzz_yyyz_1,  \
                             ta1_z_xxzz_yyz_0,   \
                             ta1_z_xxzz_yyz_1,   \
                             ta1_z_xxzz_yyzz_0,  \
                             ta1_z_xxzz_yyzz_1,  \
                             ta1_z_xxzz_yzz_0,   \
                             ta1_z_xxzz_yzz_1,   \
                             ta1_z_xxzz_yzzz_0,  \
                             ta1_z_xxzz_yzzz_1,  \
                             ta1_z_xxzz_zzz_0,   \
                             ta1_z_xxzz_zzz_1,   \
                             ta1_z_xxzz_zzzz_0,  \
                             ta1_z_xxzz_zzzz_1,  \
                             ta1_z_xzz_xxxz_0,   \
                             ta1_z_xzz_xxxz_1,   \
                             ta1_z_xzz_xxyz_0,   \
                             ta1_z_xzz_xxyz_1,   \
                             ta1_z_xzz_xxzz_0,   \
                             ta1_z_xzz_xxzz_1,   \
                             ta1_z_xzz_xyyz_0,   \
                             ta1_z_xzz_xyyz_1,   \
                             ta1_z_xzz_xyzz_0,   \
                             ta1_z_xzz_xyzz_1,   \
                             ta1_z_xzz_xzzz_0,   \
                             ta1_z_xzz_xzzz_1,   \
                             ta1_z_xzz_yyyy_0,   \
                             ta1_z_xzz_yyyy_1,   \
                             ta1_z_xzz_yyyz_0,   \
                             ta1_z_xzz_yyyz_1,   \
                             ta1_z_xzz_yyzz_0,   \
                             ta1_z_xzz_yyzz_1,   \
                             ta1_z_xzz_yzzz_0,   \
                             ta1_z_xzz_yzzz_1,   \
                             ta1_z_xzz_zzzz_0,   \
                             ta1_z_xzz_zzzz_1,   \
                             ta_xxxz_xxxx_1,     \
                             ta_xxxz_xxxy_1,     \
                             ta_xxxz_xxyy_1,     \
                             ta_xxxz_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzz_xxxx_0[i] = ta1_z_xxx_xxxx_0[i] * fe_0 - ta1_z_xxx_xxxx_1[i] * fe_0 + ta_xxxz_xxxx_1[i] + ta1_z_xxxz_xxxx_0[i] * pa_z[i] -
                                ta1_z_xxxz_xxxx_1[i] * pc_z[i];

        ta1_z_xxxzz_xxxy_0[i] = ta1_z_xxx_xxxy_0[i] * fe_0 - ta1_z_xxx_xxxy_1[i] * fe_0 + ta_xxxz_xxxy_1[i] + ta1_z_xxxz_xxxy_0[i] * pa_z[i] -
                                ta1_z_xxxz_xxxy_1[i] * pc_z[i];

        ta1_z_xxxzz_xxxz_0[i] = 2.0 * ta1_z_xzz_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xxzz_xxz_0[i] * fe_0 -
                                3.0 * ta1_z_xxzz_xxz_1[i] * fe_0 + ta1_z_xxzz_xxxz_0[i] * pa_x[i] - ta1_z_xxzz_xxxz_1[i] * pc_x[i];

        ta1_z_xxxzz_xxyy_0[i] = ta1_z_xxx_xxyy_0[i] * fe_0 - ta1_z_xxx_xxyy_1[i] * fe_0 + ta_xxxz_xxyy_1[i] + ta1_z_xxxz_xxyy_0[i] * pa_z[i] -
                                ta1_z_xxxz_xxyy_1[i] * pc_z[i];

        ta1_z_xxxzz_xxyz_0[i] = 2.0 * ta1_z_xzz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_xxzz_xyz_1[i] * fe_0 + ta1_z_xxzz_xxyz_0[i] * pa_x[i] - ta1_z_xxzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxxzz_xxzz_0[i] = 2.0 * ta1_z_xzz_xxzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xxzz_xzz_0[i] * fe_0 -
                                2.0 * ta1_z_xxzz_xzz_1[i] * fe_0 + ta1_z_xxzz_xxzz_0[i] * pa_x[i] - ta1_z_xxzz_xxzz_1[i] * pc_x[i];

        ta1_z_xxxzz_xyyy_0[i] = ta1_z_xxx_xyyy_0[i] * fe_0 - ta1_z_xxx_xyyy_1[i] * fe_0 + ta_xxxz_xyyy_1[i] + ta1_z_xxxz_xyyy_0[i] * pa_z[i] -
                                ta1_z_xxxz_xyyy_1[i] * pc_z[i];

        ta1_z_xxxzz_xyyz_0[i] = 2.0 * ta1_z_xzz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xyyz_1[i] * fe_0 + ta1_z_xxzz_yyz_0[i] * fe_0 -
                                ta1_z_xxzz_yyz_1[i] * fe_0 + ta1_z_xxzz_xyyz_0[i] * pa_x[i] - ta1_z_xxzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxxzz_xyzz_0[i] = 2.0 * ta1_z_xzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xyzz_1[i] * fe_0 + ta1_z_xxzz_yzz_0[i] * fe_0 -
                                ta1_z_xxzz_yzz_1[i] * fe_0 + ta1_z_xxzz_xyzz_0[i] * pa_x[i] - ta1_z_xxzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxxzz_xzzz_0[i] = 2.0 * ta1_z_xzz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xzzz_1[i] * fe_0 + ta1_z_xxzz_zzz_0[i] * fe_0 -
                                ta1_z_xxzz_zzz_1[i] * fe_0 + ta1_z_xxzz_xzzz_0[i] * pa_x[i] - ta1_z_xxzz_xzzz_1[i] * pc_x[i];

        ta1_z_xxxzz_yyyy_0[i] =
            2.0 * ta1_z_xzz_yyyy_0[i] * fe_0 - 2.0 * ta1_z_xzz_yyyy_1[i] * fe_0 + ta1_z_xxzz_yyyy_0[i] * pa_x[i] - ta1_z_xxzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxxzz_yyyz_0[i] =
            2.0 * ta1_z_xzz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xzz_yyyz_1[i] * fe_0 + ta1_z_xxzz_yyyz_0[i] * pa_x[i] - ta1_z_xxzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxzz_yyzz_0[i] =
            2.0 * ta1_z_xzz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_yyzz_1[i] * fe_0 + ta1_z_xxzz_yyzz_0[i] * pa_x[i] - ta1_z_xxzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxzz_yzzz_0[i] =
            2.0 * ta1_z_xzz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_yzzz_1[i] * fe_0 + ta1_z_xxzz_yzzz_0[i] * pa_x[i] - ta1_z_xxzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxzz_zzzz_0[i] =
            2.0 * ta1_z_xzz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_zzzz_1[i] * fe_0 + ta1_z_xxzz_zzzz_0[i] * pa_x[i] - ta1_z_xxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 720-735 components of targeted buffer : HG

    auto ta1_z_xxyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 720);

    auto ta1_z_xxyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 721);

    auto ta1_z_xxyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 722);

    auto ta1_z_xxyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 723);

    auto ta1_z_xxyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 724);

    auto ta1_z_xxyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 725);

    auto ta1_z_xxyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 726);

    auto ta1_z_xxyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 727);

    auto ta1_z_xxyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 728);

    auto ta1_z_xxyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 729);

    auto ta1_z_xxyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 730);

    auto ta1_z_xxyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 731);

    auto ta1_z_xxyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 732);

    auto ta1_z_xxyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 733);

    auto ta1_z_xxyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 734);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxy_xxxx_0,   \
                             ta1_z_xxy_xxxx_1,   \
                             ta1_z_xxy_xxxz_0,   \
                             ta1_z_xxy_xxxz_1,   \
                             ta1_z_xxy_xxzz_0,   \
                             ta1_z_xxy_xxzz_1,   \
                             ta1_z_xxy_xzzz_0,   \
                             ta1_z_xxy_xzzz_1,   \
                             ta1_z_xxyy_xxxx_0,  \
                             ta1_z_xxyy_xxxx_1,  \
                             ta1_z_xxyy_xxxz_0,  \
                             ta1_z_xxyy_xxxz_1,  \
                             ta1_z_xxyy_xxzz_0,  \
                             ta1_z_xxyy_xxzz_1,  \
                             ta1_z_xxyy_xzzz_0,  \
                             ta1_z_xxyy_xzzz_1,  \
                             ta1_z_xxyyy_xxxx_0, \
                             ta1_z_xxyyy_xxxy_0, \
                             ta1_z_xxyyy_xxxz_0, \
                             ta1_z_xxyyy_xxyy_0, \
                             ta1_z_xxyyy_xxyz_0, \
                             ta1_z_xxyyy_xxzz_0, \
                             ta1_z_xxyyy_xyyy_0, \
                             ta1_z_xxyyy_xyyz_0, \
                             ta1_z_xxyyy_xyzz_0, \
                             ta1_z_xxyyy_xzzz_0, \
                             ta1_z_xxyyy_yyyy_0, \
                             ta1_z_xxyyy_yyyz_0, \
                             ta1_z_xxyyy_yyzz_0, \
                             ta1_z_xxyyy_yzzz_0, \
                             ta1_z_xxyyy_zzzz_0, \
                             ta1_z_xyyy_xxxy_0,  \
                             ta1_z_xyyy_xxxy_1,  \
                             ta1_z_xyyy_xxy_0,   \
                             ta1_z_xyyy_xxy_1,   \
                             ta1_z_xyyy_xxyy_0,  \
                             ta1_z_xyyy_xxyy_1,  \
                             ta1_z_xyyy_xxyz_0,  \
                             ta1_z_xyyy_xxyz_1,  \
                             ta1_z_xyyy_xyy_0,   \
                             ta1_z_xyyy_xyy_1,   \
                             ta1_z_xyyy_xyyy_0,  \
                             ta1_z_xyyy_xyyy_1,  \
                             ta1_z_xyyy_xyyz_0,  \
                             ta1_z_xyyy_xyyz_1,  \
                             ta1_z_xyyy_xyz_0,   \
                             ta1_z_xyyy_xyz_1,   \
                             ta1_z_xyyy_xyzz_0,  \
                             ta1_z_xyyy_xyzz_1,  \
                             ta1_z_xyyy_yyy_0,   \
                             ta1_z_xyyy_yyy_1,   \
                             ta1_z_xyyy_yyyy_0,  \
                             ta1_z_xyyy_yyyy_1,  \
                             ta1_z_xyyy_yyyz_0,  \
                             ta1_z_xyyy_yyyz_1,  \
                             ta1_z_xyyy_yyz_0,   \
                             ta1_z_xyyy_yyz_1,   \
                             ta1_z_xyyy_yyzz_0,  \
                             ta1_z_xyyy_yyzz_1,  \
                             ta1_z_xyyy_yzz_0,   \
                             ta1_z_xyyy_yzz_1,   \
                             ta1_z_xyyy_yzzz_0,  \
                             ta1_z_xyyy_yzzz_1,  \
                             ta1_z_xyyy_zzzz_0,  \
                             ta1_z_xyyy_zzzz_1,  \
                             ta1_z_yyy_xxxy_0,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxyy_0,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyz_0,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xyyy_0,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyz_0,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyzz_0,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_yyyy_0,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyz_0,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyzz_0,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yzzz_0,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_zzzz_0,   \
                             ta1_z_yyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyy_xxxx_0[i] =
            2.0 * ta1_z_xxy_xxxx_0[i] * fe_0 - 2.0 * ta1_z_xxy_xxxx_1[i] * fe_0 + ta1_z_xxyy_xxxx_0[i] * pa_y[i] - ta1_z_xxyy_xxxx_1[i] * pc_y[i];

        ta1_z_xxyyy_xxxy_0[i] = ta1_z_yyy_xxxy_0[i] * fe_0 - ta1_z_yyy_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xyyy_xxy_0[i] * fe_0 -
                                3.0 * ta1_z_xyyy_xxy_1[i] * fe_0 + ta1_z_xyyy_xxxy_0[i] * pa_x[i] - ta1_z_xyyy_xxxy_1[i] * pc_x[i];

        ta1_z_xxyyy_xxxz_0[i] =
            2.0 * ta1_z_xxy_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xxy_xxxz_1[i] * fe_0 + ta1_z_xxyy_xxxz_0[i] * pa_y[i] - ta1_z_xxyy_xxxz_1[i] * pc_y[i];

        ta1_z_xxyyy_xxyy_0[i] = ta1_z_yyy_xxyy_0[i] * fe_0 - ta1_z_yyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xyyy_xyy_0[i] * fe_0 -
                                2.0 * ta1_z_xyyy_xyy_1[i] * fe_0 + ta1_z_xyyy_xxyy_0[i] * pa_x[i] - ta1_z_xyyy_xxyy_1[i] * pc_x[i];

        ta1_z_xxyyy_xxyz_0[i] = ta1_z_yyy_xxyz_0[i] * fe_0 - ta1_z_yyy_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xyyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_xyyy_xyz_1[i] * fe_0 + ta1_z_xyyy_xxyz_0[i] * pa_x[i] - ta1_z_xyyy_xxyz_1[i] * pc_x[i];

        ta1_z_xxyyy_xxzz_0[i] =
            2.0 * ta1_z_xxy_xxzz_0[i] * fe_0 - 2.0 * ta1_z_xxy_xxzz_1[i] * fe_0 + ta1_z_xxyy_xxzz_0[i] * pa_y[i] - ta1_z_xxyy_xxzz_1[i] * pc_y[i];

        ta1_z_xxyyy_xyyy_0[i] = ta1_z_yyy_xyyy_0[i] * fe_0 - ta1_z_yyy_xyyy_1[i] * fe_0 + ta1_z_xyyy_yyy_0[i] * fe_0 - ta1_z_xyyy_yyy_1[i] * fe_0 +
                                ta1_z_xyyy_xyyy_0[i] * pa_x[i] - ta1_z_xyyy_xyyy_1[i] * pc_x[i];

        ta1_z_xxyyy_xyyz_0[i] = ta1_z_yyy_xyyz_0[i] * fe_0 - ta1_z_yyy_xyyz_1[i] * fe_0 + ta1_z_xyyy_yyz_0[i] * fe_0 - ta1_z_xyyy_yyz_1[i] * fe_0 +
                                ta1_z_xyyy_xyyz_0[i] * pa_x[i] - ta1_z_xyyy_xyyz_1[i] * pc_x[i];

        ta1_z_xxyyy_xyzz_0[i] = ta1_z_yyy_xyzz_0[i] * fe_0 - ta1_z_yyy_xyzz_1[i] * fe_0 + ta1_z_xyyy_yzz_0[i] * fe_0 - ta1_z_xyyy_yzz_1[i] * fe_0 +
                                ta1_z_xyyy_xyzz_0[i] * pa_x[i] - ta1_z_xyyy_xyzz_1[i] * pc_x[i];

        ta1_z_xxyyy_xzzz_0[i] =
            2.0 * ta1_z_xxy_xzzz_0[i] * fe_0 - 2.0 * ta1_z_xxy_xzzz_1[i] * fe_0 + ta1_z_xxyy_xzzz_0[i] * pa_y[i] - ta1_z_xxyy_xzzz_1[i] * pc_y[i];

        ta1_z_xxyyy_yyyy_0[i] =
            ta1_z_yyy_yyyy_0[i] * fe_0 - ta1_z_yyy_yyyy_1[i] * fe_0 + ta1_z_xyyy_yyyy_0[i] * pa_x[i] - ta1_z_xyyy_yyyy_1[i] * pc_x[i];

        ta1_z_xxyyy_yyyz_0[i] =
            ta1_z_yyy_yyyz_0[i] * fe_0 - ta1_z_yyy_yyyz_1[i] * fe_0 + ta1_z_xyyy_yyyz_0[i] * pa_x[i] - ta1_z_xyyy_yyyz_1[i] * pc_x[i];

        ta1_z_xxyyy_yyzz_0[i] =
            ta1_z_yyy_yyzz_0[i] * fe_0 - ta1_z_yyy_yyzz_1[i] * fe_0 + ta1_z_xyyy_yyzz_0[i] * pa_x[i] - ta1_z_xyyy_yyzz_1[i] * pc_x[i];

        ta1_z_xxyyy_yzzz_0[i] =
            ta1_z_yyy_yzzz_0[i] * fe_0 - ta1_z_yyy_yzzz_1[i] * fe_0 + ta1_z_xyyy_yzzz_0[i] * pa_x[i] - ta1_z_xyyy_yzzz_1[i] * pc_x[i];

        ta1_z_xxyyy_zzzz_0[i] =
            ta1_z_yyy_zzzz_0[i] * fe_0 - ta1_z_yyy_zzzz_1[i] * fe_0 + ta1_z_xyyy_zzzz_0[i] * pa_x[i] - ta1_z_xyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 735-750 components of targeted buffer : HG

    auto ta1_z_xxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 735);

    auto ta1_z_xxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 736);

    auto ta1_z_xxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 737);

    auto ta1_z_xxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 738);

    auto ta1_z_xxyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 739);

    auto ta1_z_xxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 740);

    auto ta1_z_xxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 741);

    auto ta1_z_xxyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 742);

    auto ta1_z_xxyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 743);

    auto ta1_z_xxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 744);

    auto ta1_z_xxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 745);

    auto ta1_z_xxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 746);

    auto ta1_z_xxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 747);

    auto ta1_z_xxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 748);

    auto ta1_z_xxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 749);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxyy_xxxx_0,  \
                             ta1_z_xxyy_xxxx_1,  \
                             ta1_z_xxyy_xxxy_0,  \
                             ta1_z_xxyy_xxxy_1,  \
                             ta1_z_xxyy_xxy_0,   \
                             ta1_z_xxyy_xxy_1,   \
                             ta1_z_xxyy_xxyy_0,  \
                             ta1_z_xxyy_xxyy_1,  \
                             ta1_z_xxyy_xxyz_0,  \
                             ta1_z_xxyy_xxyz_1,  \
                             ta1_z_xxyy_xyy_0,   \
                             ta1_z_xxyy_xyy_1,   \
                             ta1_z_xxyy_xyyy_0,  \
                             ta1_z_xxyy_xyyy_1,  \
                             ta1_z_xxyy_xyyz_0,  \
                             ta1_z_xxyy_xyyz_1,  \
                             ta1_z_xxyy_xyz_0,   \
                             ta1_z_xxyy_xyz_1,   \
                             ta1_z_xxyy_xyzz_0,  \
                             ta1_z_xxyy_xyzz_1,  \
                             ta1_z_xxyy_yyyy_0,  \
                             ta1_z_xxyy_yyyy_1,  \
                             ta1_z_xxyyz_xxxx_0, \
                             ta1_z_xxyyz_xxxy_0, \
                             ta1_z_xxyyz_xxxz_0, \
                             ta1_z_xxyyz_xxyy_0, \
                             ta1_z_xxyyz_xxyz_0, \
                             ta1_z_xxyyz_xxzz_0, \
                             ta1_z_xxyyz_xyyy_0, \
                             ta1_z_xxyyz_xyyz_0, \
                             ta1_z_xxyyz_xyzz_0, \
                             ta1_z_xxyyz_xzzz_0, \
                             ta1_z_xxyyz_yyyy_0, \
                             ta1_z_xxyyz_yyyz_0, \
                             ta1_z_xxyyz_yyzz_0, \
                             ta1_z_xxyyz_yzzz_0, \
                             ta1_z_xxyyz_zzzz_0, \
                             ta1_z_xxyz_xxxz_0,  \
                             ta1_z_xxyz_xxxz_1,  \
                             ta1_z_xxyz_xxzz_0,  \
                             ta1_z_xxyz_xxzz_1,  \
                             ta1_z_xxyz_xzzz_0,  \
                             ta1_z_xxyz_xzzz_1,  \
                             ta1_z_xxz_xxxz_0,   \
                             ta1_z_xxz_xxxz_1,   \
                             ta1_z_xxz_xxzz_0,   \
                             ta1_z_xxz_xxzz_1,   \
                             ta1_z_xxz_xzzz_0,   \
                             ta1_z_xxz_xzzz_1,   \
                             ta1_z_xyyz_yyyz_0,  \
                             ta1_z_xyyz_yyyz_1,  \
                             ta1_z_xyyz_yyzz_0,  \
                             ta1_z_xyyz_yyzz_1,  \
                             ta1_z_xyyz_yzzz_0,  \
                             ta1_z_xyyz_yzzz_1,  \
                             ta1_z_xyyz_zzzz_0,  \
                             ta1_z_xyyz_zzzz_1,  \
                             ta1_z_yyz_yyyz_0,   \
                             ta1_z_yyz_yyyz_1,   \
                             ta1_z_yyz_yyzz_0,   \
                             ta1_z_yyz_yyzz_1,   \
                             ta1_z_yyz_yzzz_0,   \
                             ta1_z_yyz_yzzz_1,   \
                             ta1_z_yyz_zzzz_0,   \
                             ta1_z_yyz_zzzz_1,   \
                             ta_xxyy_xxxx_1,     \
                             ta_xxyy_xxxy_1,     \
                             ta_xxyy_xxyy_1,     \
                             ta_xxyy_xxyz_1,     \
                             ta_xxyy_xyyy_1,     \
                             ta_xxyy_xyyz_1,     \
                             ta_xxyy_xyzz_1,     \
                             ta_xxyy_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyz_xxxx_0[i] = ta_xxyy_xxxx_1[i] + ta1_z_xxyy_xxxx_0[i] * pa_z[i] - ta1_z_xxyy_xxxx_1[i] * pc_z[i];

        ta1_z_xxyyz_xxxy_0[i] = ta_xxyy_xxxy_1[i] + ta1_z_xxyy_xxxy_0[i] * pa_z[i] - ta1_z_xxyy_xxxy_1[i] * pc_z[i];

        ta1_z_xxyyz_xxxz_0[i] =
            ta1_z_xxz_xxxz_0[i] * fe_0 - ta1_z_xxz_xxxz_1[i] * fe_0 + ta1_z_xxyz_xxxz_0[i] * pa_y[i] - ta1_z_xxyz_xxxz_1[i] * pc_y[i];

        ta1_z_xxyyz_xxyy_0[i] = ta_xxyy_xxyy_1[i] + ta1_z_xxyy_xxyy_0[i] * pa_z[i] - ta1_z_xxyy_xxyy_1[i] * pc_z[i];

        ta1_z_xxyyz_xxyz_0[i] = ta1_z_xxyy_xxy_0[i] * fe_0 - ta1_z_xxyy_xxy_1[i] * fe_0 + ta_xxyy_xxyz_1[i] + ta1_z_xxyy_xxyz_0[i] * pa_z[i] -
                                ta1_z_xxyy_xxyz_1[i] * pc_z[i];

        ta1_z_xxyyz_xxzz_0[i] =
            ta1_z_xxz_xxzz_0[i] * fe_0 - ta1_z_xxz_xxzz_1[i] * fe_0 + ta1_z_xxyz_xxzz_0[i] * pa_y[i] - ta1_z_xxyz_xxzz_1[i] * pc_y[i];

        ta1_z_xxyyz_xyyy_0[i] = ta_xxyy_xyyy_1[i] + ta1_z_xxyy_xyyy_0[i] * pa_z[i] - ta1_z_xxyy_xyyy_1[i] * pc_z[i];

        ta1_z_xxyyz_xyyz_0[i] = ta1_z_xxyy_xyy_0[i] * fe_0 - ta1_z_xxyy_xyy_1[i] * fe_0 + ta_xxyy_xyyz_1[i] + ta1_z_xxyy_xyyz_0[i] * pa_z[i] -
                                ta1_z_xxyy_xyyz_1[i] * pc_z[i];

        ta1_z_xxyyz_xyzz_0[i] = 2.0 * ta1_z_xxyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxyy_xyz_1[i] * fe_0 + ta_xxyy_xyzz_1[i] +
                                ta1_z_xxyy_xyzz_0[i] * pa_z[i] - ta1_z_xxyy_xyzz_1[i] * pc_z[i];

        ta1_z_xxyyz_xzzz_0[i] =
            ta1_z_xxz_xzzz_0[i] * fe_0 - ta1_z_xxz_xzzz_1[i] * fe_0 + ta1_z_xxyz_xzzz_0[i] * pa_y[i] - ta1_z_xxyz_xzzz_1[i] * pc_y[i];

        ta1_z_xxyyz_yyyy_0[i] = ta_xxyy_yyyy_1[i] + ta1_z_xxyy_yyyy_0[i] * pa_z[i] - ta1_z_xxyy_yyyy_1[i] * pc_z[i];

        ta1_z_xxyyz_yyyz_0[i] =
            ta1_z_yyz_yyyz_0[i] * fe_0 - ta1_z_yyz_yyyz_1[i] * fe_0 + ta1_z_xyyz_yyyz_0[i] * pa_x[i] - ta1_z_xyyz_yyyz_1[i] * pc_x[i];

        ta1_z_xxyyz_yyzz_0[i] =
            ta1_z_yyz_yyzz_0[i] * fe_0 - ta1_z_yyz_yyzz_1[i] * fe_0 + ta1_z_xyyz_yyzz_0[i] * pa_x[i] - ta1_z_xyyz_yyzz_1[i] * pc_x[i];

        ta1_z_xxyyz_yzzz_0[i] =
            ta1_z_yyz_yzzz_0[i] * fe_0 - ta1_z_yyz_yzzz_1[i] * fe_0 + ta1_z_xyyz_yzzz_0[i] * pa_x[i] - ta1_z_xyyz_yzzz_1[i] * pc_x[i];

        ta1_z_xxyyz_zzzz_0[i] =
            ta1_z_yyz_zzzz_0[i] * fe_0 - ta1_z_yyz_zzzz_1[i] * fe_0 + ta1_z_xyyz_zzzz_0[i] * pa_x[i] - ta1_z_xyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 750-765 components of targeted buffer : HG

    auto ta1_z_xxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 750);

    auto ta1_z_xxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 751);

    auto ta1_z_xxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 752);

    auto ta1_z_xxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 753);

    auto ta1_z_xxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 754);

    auto ta1_z_xxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 755);

    auto ta1_z_xxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 756);

    auto ta1_z_xxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 757);

    auto ta1_z_xxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 758);

    auto ta1_z_xxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 759);

    auto ta1_z_xxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 760);

    auto ta1_z_xxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 761);

    auto ta1_z_xxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 762);

    auto ta1_z_xxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 763);

    auto ta1_z_xxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 764);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxyzz_xxxx_0, \
                             ta1_z_xxyzz_xxxy_0, \
                             ta1_z_xxyzz_xxxz_0, \
                             ta1_z_xxyzz_xxyy_0, \
                             ta1_z_xxyzz_xxyz_0, \
                             ta1_z_xxyzz_xxzz_0, \
                             ta1_z_xxyzz_xyyy_0, \
                             ta1_z_xxyzz_xyyz_0, \
                             ta1_z_xxyzz_xyzz_0, \
                             ta1_z_xxyzz_xzzz_0, \
                             ta1_z_xxyzz_yyyy_0, \
                             ta1_z_xxyzz_yyyz_0, \
                             ta1_z_xxyzz_yyzz_0, \
                             ta1_z_xxyzz_yzzz_0, \
                             ta1_z_xxyzz_zzzz_0, \
                             ta1_z_xxzz_xxx_0,   \
                             ta1_z_xxzz_xxx_1,   \
                             ta1_z_xxzz_xxxx_0,  \
                             ta1_z_xxzz_xxxx_1,  \
                             ta1_z_xxzz_xxxy_0,  \
                             ta1_z_xxzz_xxxy_1,  \
                             ta1_z_xxzz_xxxz_0,  \
                             ta1_z_xxzz_xxxz_1,  \
                             ta1_z_xxzz_xxy_0,   \
                             ta1_z_xxzz_xxy_1,   \
                             ta1_z_xxzz_xxyy_0,  \
                             ta1_z_xxzz_xxyy_1,  \
                             ta1_z_xxzz_xxyz_0,  \
                             ta1_z_xxzz_xxyz_1,  \
                             ta1_z_xxzz_xxz_0,   \
                             ta1_z_xxzz_xxz_1,   \
                             ta1_z_xxzz_xxzz_0,  \
                             ta1_z_xxzz_xxzz_1,  \
                             ta1_z_xxzz_xyy_0,   \
                             ta1_z_xxzz_xyy_1,   \
                             ta1_z_xxzz_xyyy_0,  \
                             ta1_z_xxzz_xyyy_1,  \
                             ta1_z_xxzz_xyyz_0,  \
                             ta1_z_xxzz_xyyz_1,  \
                             ta1_z_xxzz_xyz_0,   \
                             ta1_z_xxzz_xyz_1,   \
                             ta1_z_xxzz_xyzz_0,  \
                             ta1_z_xxzz_xyzz_1,  \
                             ta1_z_xxzz_xzz_0,   \
                             ta1_z_xxzz_xzz_1,   \
                             ta1_z_xxzz_xzzz_0,  \
                             ta1_z_xxzz_xzzz_1,  \
                             ta1_z_xxzz_zzzz_0,  \
                             ta1_z_xxzz_zzzz_1,  \
                             ta1_z_xyzz_yyyy_0,  \
                             ta1_z_xyzz_yyyy_1,  \
                             ta1_z_xyzz_yyyz_0,  \
                             ta1_z_xyzz_yyyz_1,  \
                             ta1_z_xyzz_yyzz_0,  \
                             ta1_z_xyzz_yyzz_1,  \
                             ta1_z_xyzz_yzzz_0,  \
                             ta1_z_xyzz_yzzz_1,  \
                             ta1_z_yzz_yyyy_0,   \
                             ta1_z_yzz_yyyy_1,   \
                             ta1_z_yzz_yyyz_0,   \
                             ta1_z_yzz_yyyz_1,   \
                             ta1_z_yzz_yyzz_0,   \
                             ta1_z_yzz_yyzz_1,   \
                             ta1_z_yzz_yzzz_0,   \
                             ta1_z_yzz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzz_xxxx_0[i] = ta1_z_xxzz_xxxx_0[i] * pa_y[i] - ta1_z_xxzz_xxxx_1[i] * pc_y[i];

        ta1_z_xxyzz_xxxy_0[i] =
            ta1_z_xxzz_xxx_0[i] * fe_0 - ta1_z_xxzz_xxx_1[i] * fe_0 + ta1_z_xxzz_xxxy_0[i] * pa_y[i] - ta1_z_xxzz_xxxy_1[i] * pc_y[i];

        ta1_z_xxyzz_xxxz_0[i] = ta1_z_xxzz_xxxz_0[i] * pa_y[i] - ta1_z_xxzz_xxxz_1[i] * pc_y[i];

        ta1_z_xxyzz_xxyy_0[i] =
            2.0 * ta1_z_xxzz_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxzz_xxy_1[i] * fe_0 + ta1_z_xxzz_xxyy_0[i] * pa_y[i] - ta1_z_xxzz_xxyy_1[i] * pc_y[i];

        ta1_z_xxyzz_xxyz_0[i] =
            ta1_z_xxzz_xxz_0[i] * fe_0 - ta1_z_xxzz_xxz_1[i] * fe_0 + ta1_z_xxzz_xxyz_0[i] * pa_y[i] - ta1_z_xxzz_xxyz_1[i] * pc_y[i];

        ta1_z_xxyzz_xxzz_0[i] = ta1_z_xxzz_xxzz_0[i] * pa_y[i] - ta1_z_xxzz_xxzz_1[i] * pc_y[i];

        ta1_z_xxyzz_xyyy_0[i] =
            3.0 * ta1_z_xxzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xyy_1[i] * fe_0 + ta1_z_xxzz_xyyy_0[i] * pa_y[i] - ta1_z_xxzz_xyyy_1[i] * pc_y[i];

        ta1_z_xxyzz_xyyz_0[i] =
            2.0 * ta1_z_xxzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxzz_xyz_1[i] * fe_0 + ta1_z_xxzz_xyyz_0[i] * pa_y[i] - ta1_z_xxzz_xyyz_1[i] * pc_y[i];

        ta1_z_xxyzz_xyzz_0[i] =
            ta1_z_xxzz_xzz_0[i] * fe_0 - ta1_z_xxzz_xzz_1[i] * fe_0 + ta1_z_xxzz_xyzz_0[i] * pa_y[i] - ta1_z_xxzz_xyzz_1[i] * pc_y[i];

        ta1_z_xxyzz_xzzz_0[i] = ta1_z_xxzz_xzzz_0[i] * pa_y[i] - ta1_z_xxzz_xzzz_1[i] * pc_y[i];

        ta1_z_xxyzz_yyyy_0[i] =
            ta1_z_yzz_yyyy_0[i] * fe_0 - ta1_z_yzz_yyyy_1[i] * fe_0 + ta1_z_xyzz_yyyy_0[i] * pa_x[i] - ta1_z_xyzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxyzz_yyyz_0[i] =
            ta1_z_yzz_yyyz_0[i] * fe_0 - ta1_z_yzz_yyyz_1[i] * fe_0 + ta1_z_xyzz_yyyz_0[i] * pa_x[i] - ta1_z_xyzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxyzz_yyzz_0[i] =
            ta1_z_yzz_yyzz_0[i] * fe_0 - ta1_z_yzz_yyzz_1[i] * fe_0 + ta1_z_xyzz_yyzz_0[i] * pa_x[i] - ta1_z_xyzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxyzz_yzzz_0[i] =
            ta1_z_yzz_yzzz_0[i] * fe_0 - ta1_z_yzz_yzzz_1[i] * fe_0 + ta1_z_xyzz_yzzz_0[i] * pa_x[i] - ta1_z_xyzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxyzz_zzzz_0[i] = ta1_z_xxzz_zzzz_0[i] * pa_y[i] - ta1_z_xxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 765-780 components of targeted buffer : HG

    auto ta1_z_xxzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 765);

    auto ta1_z_xxzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 766);

    auto ta1_z_xxzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 767);

    auto ta1_z_xxzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 768);

    auto ta1_z_xxzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 769);

    auto ta1_z_xxzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 770);

    auto ta1_z_xxzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 771);

    auto ta1_z_xxzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 772);

    auto ta1_z_xxzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 773);

    auto ta1_z_xxzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 774);

    auto ta1_z_xxzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 775);

    auto ta1_z_xxzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 776);

    auto ta1_z_xxzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 777);

    auto ta1_z_xxzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 778);

    auto ta1_z_xxzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 779);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxz_xxxx_0,   \
                             ta1_z_xxz_xxxx_1,   \
                             ta1_z_xxz_xxxy_0,   \
                             ta1_z_xxz_xxxy_1,   \
                             ta1_z_xxz_xxyy_0,   \
                             ta1_z_xxz_xxyy_1,   \
                             ta1_z_xxz_xyyy_0,   \
                             ta1_z_xxz_xyyy_1,   \
                             ta1_z_xxzz_xxxx_0,  \
                             ta1_z_xxzz_xxxx_1,  \
                             ta1_z_xxzz_xxxy_0,  \
                             ta1_z_xxzz_xxxy_1,  \
                             ta1_z_xxzz_xxyy_0,  \
                             ta1_z_xxzz_xxyy_1,  \
                             ta1_z_xxzz_xyyy_0,  \
                             ta1_z_xxzz_xyyy_1,  \
                             ta1_z_xxzzz_xxxx_0, \
                             ta1_z_xxzzz_xxxy_0, \
                             ta1_z_xxzzz_xxxz_0, \
                             ta1_z_xxzzz_xxyy_0, \
                             ta1_z_xxzzz_xxyz_0, \
                             ta1_z_xxzzz_xxzz_0, \
                             ta1_z_xxzzz_xyyy_0, \
                             ta1_z_xxzzz_xyyz_0, \
                             ta1_z_xxzzz_xyzz_0, \
                             ta1_z_xxzzz_xzzz_0, \
                             ta1_z_xxzzz_yyyy_0, \
                             ta1_z_xxzzz_yyyz_0, \
                             ta1_z_xxzzz_yyzz_0, \
                             ta1_z_xxzzz_yzzz_0, \
                             ta1_z_xxzzz_zzzz_0, \
                             ta1_z_xzzz_xxxz_0,  \
                             ta1_z_xzzz_xxxz_1,  \
                             ta1_z_xzzz_xxyz_0,  \
                             ta1_z_xzzz_xxyz_1,  \
                             ta1_z_xzzz_xxz_0,   \
                             ta1_z_xzzz_xxz_1,   \
                             ta1_z_xzzz_xxzz_0,  \
                             ta1_z_xzzz_xxzz_1,  \
                             ta1_z_xzzz_xyyz_0,  \
                             ta1_z_xzzz_xyyz_1,  \
                             ta1_z_xzzz_xyz_0,   \
                             ta1_z_xzzz_xyz_1,   \
                             ta1_z_xzzz_xyzz_0,  \
                             ta1_z_xzzz_xyzz_1,  \
                             ta1_z_xzzz_xzz_0,   \
                             ta1_z_xzzz_xzz_1,   \
                             ta1_z_xzzz_xzzz_0,  \
                             ta1_z_xzzz_xzzz_1,  \
                             ta1_z_xzzz_yyyy_0,  \
                             ta1_z_xzzz_yyyy_1,  \
                             ta1_z_xzzz_yyyz_0,  \
                             ta1_z_xzzz_yyyz_1,  \
                             ta1_z_xzzz_yyz_0,   \
                             ta1_z_xzzz_yyz_1,   \
                             ta1_z_xzzz_yyzz_0,  \
                             ta1_z_xzzz_yyzz_1,  \
                             ta1_z_xzzz_yzz_0,   \
                             ta1_z_xzzz_yzz_1,   \
                             ta1_z_xzzz_yzzz_0,  \
                             ta1_z_xzzz_yzzz_1,  \
                             ta1_z_xzzz_zzz_0,   \
                             ta1_z_xzzz_zzz_1,   \
                             ta1_z_xzzz_zzzz_0,  \
                             ta1_z_xzzz_zzzz_1,  \
                             ta1_z_zzz_xxxz_0,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxyz_0,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxzz_0,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xyyz_0,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyzz_0,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xzzz_0,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_yyyy_0,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyz_0,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyzz_0,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yzzz_0,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_zzzz_0,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta_xxzz_xxxx_1,     \
                             ta_xxzz_xxxy_1,     \
                             ta_xxzz_xxyy_1,     \
                             ta_xxzz_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzz_xxxx_0[i] = 2.0 * ta1_z_xxz_xxxx_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxxx_1[i] * fe_0 + ta_xxzz_xxxx_1[i] +
                                ta1_z_xxzz_xxxx_0[i] * pa_z[i] - ta1_z_xxzz_xxxx_1[i] * pc_z[i];

        ta1_z_xxzzz_xxxy_0[i] = 2.0 * ta1_z_xxz_xxxy_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxxy_1[i] * fe_0 + ta_xxzz_xxxy_1[i] +
                                ta1_z_xxzz_xxxy_0[i] * pa_z[i] - ta1_z_xxzz_xxxy_1[i] * pc_z[i];

        ta1_z_xxzzz_xxxz_0[i] = ta1_z_zzz_xxxz_0[i] * fe_0 - ta1_z_zzz_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xzzz_xxz_0[i] * fe_0 -
                                3.0 * ta1_z_xzzz_xxz_1[i] * fe_0 + ta1_z_xzzz_xxxz_0[i] * pa_x[i] - ta1_z_xzzz_xxxz_1[i] * pc_x[i];

        ta1_z_xxzzz_xxyy_0[i] = 2.0 * ta1_z_xxz_xxyy_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxyy_1[i] * fe_0 + ta_xxzz_xxyy_1[i] +
                                ta1_z_xxzz_xxyy_0[i] * pa_z[i] - ta1_z_xxzz_xxyy_1[i] * pc_z[i];

        ta1_z_xxzzz_xxyz_0[i] = ta1_z_zzz_xxyz_0[i] * fe_0 - ta1_z_zzz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_xzzz_xyz_1[i] * fe_0 + ta1_z_xzzz_xxyz_0[i] * pa_x[i] - ta1_z_xzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxzzz_xxzz_0[i] = ta1_z_zzz_xxzz_0[i] * fe_0 - ta1_z_zzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xzzz_xzz_0[i] * fe_0 -
                                2.0 * ta1_z_xzzz_xzz_1[i] * fe_0 + ta1_z_xzzz_xxzz_0[i] * pa_x[i] - ta1_z_xzzz_xxzz_1[i] * pc_x[i];

        ta1_z_xxzzz_xyyy_0[i] = 2.0 * ta1_z_xxz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_xxz_xyyy_1[i] * fe_0 + ta_xxzz_xyyy_1[i] +
                                ta1_z_xxzz_xyyy_0[i] * pa_z[i] - ta1_z_xxzz_xyyy_1[i] * pc_z[i];

        ta1_z_xxzzz_xyyz_0[i] = ta1_z_zzz_xyyz_0[i] * fe_0 - ta1_z_zzz_xyyz_1[i] * fe_0 + ta1_z_xzzz_yyz_0[i] * fe_0 - ta1_z_xzzz_yyz_1[i] * fe_0 +
                                ta1_z_xzzz_xyyz_0[i] * pa_x[i] - ta1_z_xzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxzzz_xyzz_0[i] = ta1_z_zzz_xyzz_0[i] * fe_0 - ta1_z_zzz_xyzz_1[i] * fe_0 + ta1_z_xzzz_yzz_0[i] * fe_0 - ta1_z_xzzz_yzz_1[i] * fe_0 +
                                ta1_z_xzzz_xyzz_0[i] * pa_x[i] - ta1_z_xzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxzzz_xzzz_0[i] = ta1_z_zzz_xzzz_0[i] * fe_0 - ta1_z_zzz_xzzz_1[i] * fe_0 + ta1_z_xzzz_zzz_0[i] * fe_0 - ta1_z_xzzz_zzz_1[i] * fe_0 +
                                ta1_z_xzzz_xzzz_0[i] * pa_x[i] - ta1_z_xzzz_xzzz_1[i] * pc_x[i];

        ta1_z_xxzzz_yyyy_0[i] =
            ta1_z_zzz_yyyy_0[i] * fe_0 - ta1_z_zzz_yyyy_1[i] * fe_0 + ta1_z_xzzz_yyyy_0[i] * pa_x[i] - ta1_z_xzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxzzz_yyyz_0[i] =
            ta1_z_zzz_yyyz_0[i] * fe_0 - ta1_z_zzz_yyyz_1[i] * fe_0 + ta1_z_xzzz_yyyz_0[i] * pa_x[i] - ta1_z_xzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxzzz_yyzz_0[i] =
            ta1_z_zzz_yyzz_0[i] * fe_0 - ta1_z_zzz_yyzz_1[i] * fe_0 + ta1_z_xzzz_yyzz_0[i] * pa_x[i] - ta1_z_xzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxzzz_yzzz_0[i] =
            ta1_z_zzz_yzzz_0[i] * fe_0 - ta1_z_zzz_yzzz_1[i] * fe_0 + ta1_z_xzzz_yzzz_0[i] * pa_x[i] - ta1_z_xzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxzzz_zzzz_0[i] =
            ta1_z_zzz_zzzz_0[i] * fe_0 - ta1_z_zzz_zzzz_1[i] * fe_0 + ta1_z_xzzz_zzzz_0[i] * pa_x[i] - ta1_z_xzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 780-795 components of targeted buffer : HG

    auto ta1_z_xyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 780);

    auto ta1_z_xyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 781);

    auto ta1_z_xyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 782);

    auto ta1_z_xyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 783);

    auto ta1_z_xyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 784);

    auto ta1_z_xyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 785);

    auto ta1_z_xyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 786);

    auto ta1_z_xyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 787);

    auto ta1_z_xyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 788);

    auto ta1_z_xyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 789);

    auto ta1_z_xyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 790);

    auto ta1_z_xyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 791);

    auto ta1_z_xyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 792);

    auto ta1_z_xyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 793);

    auto ta1_z_xyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 794);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyyyy_xxxx_0, \
                             ta1_z_xyyyy_xxxy_0, \
                             ta1_z_xyyyy_xxxz_0, \
                             ta1_z_xyyyy_xxyy_0, \
                             ta1_z_xyyyy_xxyz_0, \
                             ta1_z_xyyyy_xxzz_0, \
                             ta1_z_xyyyy_xyyy_0, \
                             ta1_z_xyyyy_xyyz_0, \
                             ta1_z_xyyyy_xyzz_0, \
                             ta1_z_xyyyy_xzzz_0, \
                             ta1_z_xyyyy_yyyy_0, \
                             ta1_z_xyyyy_yyyz_0, \
                             ta1_z_xyyyy_yyzz_0, \
                             ta1_z_xyyyy_yzzz_0, \
                             ta1_z_xyyyy_zzzz_0, \
                             ta1_z_yyyy_xxx_0,   \
                             ta1_z_yyyy_xxx_1,   \
                             ta1_z_yyyy_xxxx_0,  \
                             ta1_z_yyyy_xxxx_1,  \
                             ta1_z_yyyy_xxxy_0,  \
                             ta1_z_yyyy_xxxy_1,  \
                             ta1_z_yyyy_xxxz_0,  \
                             ta1_z_yyyy_xxxz_1,  \
                             ta1_z_yyyy_xxy_0,   \
                             ta1_z_yyyy_xxy_1,   \
                             ta1_z_yyyy_xxyy_0,  \
                             ta1_z_yyyy_xxyy_1,  \
                             ta1_z_yyyy_xxyz_0,  \
                             ta1_z_yyyy_xxyz_1,  \
                             ta1_z_yyyy_xxz_0,   \
                             ta1_z_yyyy_xxz_1,   \
                             ta1_z_yyyy_xxzz_0,  \
                             ta1_z_yyyy_xxzz_1,  \
                             ta1_z_yyyy_xyy_0,   \
                             ta1_z_yyyy_xyy_1,   \
                             ta1_z_yyyy_xyyy_0,  \
                             ta1_z_yyyy_xyyy_1,  \
                             ta1_z_yyyy_xyyz_0,  \
                             ta1_z_yyyy_xyyz_1,  \
                             ta1_z_yyyy_xyz_0,   \
                             ta1_z_yyyy_xyz_1,   \
                             ta1_z_yyyy_xyzz_0,  \
                             ta1_z_yyyy_xyzz_1,  \
                             ta1_z_yyyy_xzz_0,   \
                             ta1_z_yyyy_xzz_1,   \
                             ta1_z_yyyy_xzzz_0,  \
                             ta1_z_yyyy_xzzz_1,  \
                             ta1_z_yyyy_yyy_0,   \
                             ta1_z_yyyy_yyy_1,   \
                             ta1_z_yyyy_yyyy_0,  \
                             ta1_z_yyyy_yyyy_1,  \
                             ta1_z_yyyy_yyyz_0,  \
                             ta1_z_yyyy_yyyz_1,  \
                             ta1_z_yyyy_yyz_0,   \
                             ta1_z_yyyy_yyz_1,   \
                             ta1_z_yyyy_yyzz_0,  \
                             ta1_z_yyyy_yyzz_1,  \
                             ta1_z_yyyy_yzz_0,   \
                             ta1_z_yyyy_yzz_1,   \
                             ta1_z_yyyy_yzzz_0,  \
                             ta1_z_yyyy_yzzz_1,  \
                             ta1_z_yyyy_zzz_0,   \
                             ta1_z_yyyy_zzz_1,   \
                             ta1_z_yyyy_zzzz_0,  \
                             ta1_z_yyyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyy_xxxx_0[i] =
            4.0 * ta1_z_yyyy_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyyy_xxx_1[i] * fe_0 + ta1_z_yyyy_xxxx_0[i] * pa_x[i] - ta1_z_yyyy_xxxx_1[i] * pc_x[i];

        ta1_z_xyyyy_xxxy_0[i] =
            3.0 * ta1_z_yyyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyyy_xxy_1[i] * fe_0 + ta1_z_yyyy_xxxy_0[i] * pa_x[i] - ta1_z_yyyy_xxxy_1[i] * pc_x[i];

        ta1_z_xyyyy_xxxz_0[i] =
            3.0 * ta1_z_yyyy_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyyy_xxz_1[i] * fe_0 + ta1_z_yyyy_xxxz_0[i] * pa_x[i] - ta1_z_yyyy_xxxz_1[i] * pc_x[i];

        ta1_z_xyyyy_xxyy_0[i] =
            2.0 * ta1_z_yyyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xyy_1[i] * fe_0 + ta1_z_yyyy_xxyy_0[i] * pa_x[i] - ta1_z_yyyy_xxyy_1[i] * pc_x[i];

        ta1_z_xyyyy_xxyz_0[i] =
            2.0 * ta1_z_yyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xyz_1[i] * fe_0 + ta1_z_yyyy_xxyz_0[i] * pa_x[i] - ta1_z_yyyy_xxyz_1[i] * pc_x[i];

        ta1_z_xyyyy_xxzz_0[i] =
            2.0 * ta1_z_yyyy_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xzz_1[i] * fe_0 + ta1_z_yyyy_xxzz_0[i] * pa_x[i] - ta1_z_yyyy_xxzz_1[i] * pc_x[i];

        ta1_z_xyyyy_xyyy_0[i] =
            ta1_z_yyyy_yyy_0[i] * fe_0 - ta1_z_yyyy_yyy_1[i] * fe_0 + ta1_z_yyyy_xyyy_0[i] * pa_x[i] - ta1_z_yyyy_xyyy_1[i] * pc_x[i];

        ta1_z_xyyyy_xyyz_0[i] =
            ta1_z_yyyy_yyz_0[i] * fe_0 - ta1_z_yyyy_yyz_1[i] * fe_0 + ta1_z_yyyy_xyyz_0[i] * pa_x[i] - ta1_z_yyyy_xyyz_1[i] * pc_x[i];

        ta1_z_xyyyy_xyzz_0[i] =
            ta1_z_yyyy_yzz_0[i] * fe_0 - ta1_z_yyyy_yzz_1[i] * fe_0 + ta1_z_yyyy_xyzz_0[i] * pa_x[i] - ta1_z_yyyy_xyzz_1[i] * pc_x[i];

        ta1_z_xyyyy_xzzz_0[i] =
            ta1_z_yyyy_zzz_0[i] * fe_0 - ta1_z_yyyy_zzz_1[i] * fe_0 + ta1_z_yyyy_xzzz_0[i] * pa_x[i] - ta1_z_yyyy_xzzz_1[i] * pc_x[i];

        ta1_z_xyyyy_yyyy_0[i] = ta1_z_yyyy_yyyy_0[i] * pa_x[i] - ta1_z_yyyy_yyyy_1[i] * pc_x[i];

        ta1_z_xyyyy_yyyz_0[i] = ta1_z_yyyy_yyyz_0[i] * pa_x[i] - ta1_z_yyyy_yyyz_1[i] * pc_x[i];

        ta1_z_xyyyy_yyzz_0[i] = ta1_z_yyyy_yyzz_0[i] * pa_x[i] - ta1_z_yyyy_yyzz_1[i] * pc_x[i];

        ta1_z_xyyyy_yzzz_0[i] = ta1_z_yyyy_yzzz_0[i] * pa_x[i] - ta1_z_yyyy_yzzz_1[i] * pc_x[i];

        ta1_z_xyyyy_zzzz_0[i] = ta1_z_yyyy_zzzz_0[i] * pa_x[i] - ta1_z_yyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 795-810 components of targeted buffer : HG

    auto ta1_z_xyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 795);

    auto ta1_z_xyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 796);

    auto ta1_z_xyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 797);

    auto ta1_z_xyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 798);

    auto ta1_z_xyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 799);

    auto ta1_z_xyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 800);

    auto ta1_z_xyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 801);

    auto ta1_z_xyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 802);

    auto ta1_z_xyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 803);

    auto ta1_z_xyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 804);

    auto ta1_z_xyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 805);

    auto ta1_z_xyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 806);

    auto ta1_z_xyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 807);

    auto ta1_z_xyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 808);

    auto ta1_z_xyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 809);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xyyy_xxxx_0,  \
                             ta1_z_xyyy_xxxx_1,  \
                             ta1_z_xyyy_xxxy_0,  \
                             ta1_z_xyyy_xxxy_1,  \
                             ta1_z_xyyy_xxyy_0,  \
                             ta1_z_xyyy_xxyy_1,  \
                             ta1_z_xyyy_xyyy_0,  \
                             ta1_z_xyyy_xyyy_1,  \
                             ta1_z_xyyyz_xxxx_0, \
                             ta1_z_xyyyz_xxxy_0, \
                             ta1_z_xyyyz_xxxz_0, \
                             ta1_z_xyyyz_xxyy_0, \
                             ta1_z_xyyyz_xxyz_0, \
                             ta1_z_xyyyz_xxzz_0, \
                             ta1_z_xyyyz_xyyy_0, \
                             ta1_z_xyyyz_xyyz_0, \
                             ta1_z_xyyyz_xyzz_0, \
                             ta1_z_xyyyz_xzzz_0, \
                             ta1_z_xyyyz_yyyy_0, \
                             ta1_z_xyyyz_yyyz_0, \
                             ta1_z_xyyyz_yyzz_0, \
                             ta1_z_xyyyz_yzzz_0, \
                             ta1_z_xyyyz_zzzz_0, \
                             ta1_z_yyyz_xxxz_0,  \
                             ta1_z_yyyz_xxxz_1,  \
                             ta1_z_yyyz_xxyz_0,  \
                             ta1_z_yyyz_xxyz_1,  \
                             ta1_z_yyyz_xxz_0,   \
                             ta1_z_yyyz_xxz_1,   \
                             ta1_z_yyyz_xxzz_0,  \
                             ta1_z_yyyz_xxzz_1,  \
                             ta1_z_yyyz_xyyz_0,  \
                             ta1_z_yyyz_xyyz_1,  \
                             ta1_z_yyyz_xyz_0,   \
                             ta1_z_yyyz_xyz_1,   \
                             ta1_z_yyyz_xyzz_0,  \
                             ta1_z_yyyz_xyzz_1,  \
                             ta1_z_yyyz_xzz_0,   \
                             ta1_z_yyyz_xzz_1,   \
                             ta1_z_yyyz_xzzz_0,  \
                             ta1_z_yyyz_xzzz_1,  \
                             ta1_z_yyyz_yyyy_0,  \
                             ta1_z_yyyz_yyyy_1,  \
                             ta1_z_yyyz_yyyz_0,  \
                             ta1_z_yyyz_yyyz_1,  \
                             ta1_z_yyyz_yyz_0,   \
                             ta1_z_yyyz_yyz_1,   \
                             ta1_z_yyyz_yyzz_0,  \
                             ta1_z_yyyz_yyzz_1,  \
                             ta1_z_yyyz_yzz_0,   \
                             ta1_z_yyyz_yzz_1,   \
                             ta1_z_yyyz_yzzz_0,  \
                             ta1_z_yyyz_yzzz_1,  \
                             ta1_z_yyyz_zzz_0,   \
                             ta1_z_yyyz_zzz_1,   \
                             ta1_z_yyyz_zzzz_0,  \
                             ta1_z_yyyz_zzzz_1,  \
                             ta_xyyy_xxxx_1,     \
                             ta_xyyy_xxxy_1,     \
                             ta_xyyy_xxyy_1,     \
                             ta_xyyy_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyz_xxxx_0[i] = ta_xyyy_xxxx_1[i] + ta1_z_xyyy_xxxx_0[i] * pa_z[i] - ta1_z_xyyy_xxxx_1[i] * pc_z[i];

        ta1_z_xyyyz_xxxy_0[i] = ta_xyyy_xxxy_1[i] + ta1_z_xyyy_xxxy_0[i] * pa_z[i] - ta1_z_xyyy_xxxy_1[i] * pc_z[i];

        ta1_z_xyyyz_xxxz_0[i] =
            3.0 * ta1_z_yyyz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyyz_xxz_1[i] * fe_0 + ta1_z_yyyz_xxxz_0[i] * pa_x[i] - ta1_z_yyyz_xxxz_1[i] * pc_x[i];

        ta1_z_xyyyz_xxyy_0[i] = ta_xyyy_xxyy_1[i] + ta1_z_xyyy_xxyy_0[i] * pa_z[i] - ta1_z_xyyy_xxyy_1[i] * pc_z[i];

        ta1_z_xyyyz_xxyz_0[i] =
            2.0 * ta1_z_yyyz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xyz_1[i] * fe_0 + ta1_z_yyyz_xxyz_0[i] * pa_x[i] - ta1_z_yyyz_xxyz_1[i] * pc_x[i];

        ta1_z_xyyyz_xxzz_0[i] =
            2.0 * ta1_z_yyyz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xzz_1[i] * fe_0 + ta1_z_yyyz_xxzz_0[i] * pa_x[i] - ta1_z_yyyz_xxzz_1[i] * pc_x[i];

        ta1_z_xyyyz_xyyy_0[i] = ta_xyyy_xyyy_1[i] + ta1_z_xyyy_xyyy_0[i] * pa_z[i] - ta1_z_xyyy_xyyy_1[i] * pc_z[i];

        ta1_z_xyyyz_xyyz_0[i] =
            ta1_z_yyyz_yyz_0[i] * fe_0 - ta1_z_yyyz_yyz_1[i] * fe_0 + ta1_z_yyyz_xyyz_0[i] * pa_x[i] - ta1_z_yyyz_xyyz_1[i] * pc_x[i];

        ta1_z_xyyyz_xyzz_0[i] =
            ta1_z_yyyz_yzz_0[i] * fe_0 - ta1_z_yyyz_yzz_1[i] * fe_0 + ta1_z_yyyz_xyzz_0[i] * pa_x[i] - ta1_z_yyyz_xyzz_1[i] * pc_x[i];

        ta1_z_xyyyz_xzzz_0[i] =
            ta1_z_yyyz_zzz_0[i] * fe_0 - ta1_z_yyyz_zzz_1[i] * fe_0 + ta1_z_yyyz_xzzz_0[i] * pa_x[i] - ta1_z_yyyz_xzzz_1[i] * pc_x[i];

        ta1_z_xyyyz_yyyy_0[i] = ta1_z_yyyz_yyyy_0[i] * pa_x[i] - ta1_z_yyyz_yyyy_1[i] * pc_x[i];

        ta1_z_xyyyz_yyyz_0[i] = ta1_z_yyyz_yyyz_0[i] * pa_x[i] - ta1_z_yyyz_yyyz_1[i] * pc_x[i];

        ta1_z_xyyyz_yyzz_0[i] = ta1_z_yyyz_yyzz_0[i] * pa_x[i] - ta1_z_yyyz_yyzz_1[i] * pc_x[i];

        ta1_z_xyyyz_yzzz_0[i] = ta1_z_yyyz_yzzz_0[i] * pa_x[i] - ta1_z_yyyz_yzzz_1[i] * pc_x[i];

        ta1_z_xyyyz_zzzz_0[i] = ta1_z_yyyz_zzzz_0[i] * pa_x[i] - ta1_z_yyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 810-825 components of targeted buffer : HG

    auto ta1_z_xyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 810);

    auto ta1_z_xyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 811);

    auto ta1_z_xyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 812);

    auto ta1_z_xyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 813);

    auto ta1_z_xyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 814);

    auto ta1_z_xyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 815);

    auto ta1_z_xyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 816);

    auto ta1_z_xyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 817);

    auto ta1_z_xyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 818);

    auto ta1_z_xyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 819);

    auto ta1_z_xyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 820);

    auto ta1_z_xyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 821);

    auto ta1_z_xyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 822);

    auto ta1_z_xyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 823);

    auto ta1_z_xyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 824);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyyzz_xxxx_0, \
                             ta1_z_xyyzz_xxxy_0, \
                             ta1_z_xyyzz_xxxz_0, \
                             ta1_z_xyyzz_xxyy_0, \
                             ta1_z_xyyzz_xxyz_0, \
                             ta1_z_xyyzz_xxzz_0, \
                             ta1_z_xyyzz_xyyy_0, \
                             ta1_z_xyyzz_xyyz_0, \
                             ta1_z_xyyzz_xyzz_0, \
                             ta1_z_xyyzz_xzzz_0, \
                             ta1_z_xyyzz_yyyy_0, \
                             ta1_z_xyyzz_yyyz_0, \
                             ta1_z_xyyzz_yyzz_0, \
                             ta1_z_xyyzz_yzzz_0, \
                             ta1_z_xyyzz_zzzz_0, \
                             ta1_z_yyzz_xxx_0,   \
                             ta1_z_yyzz_xxx_1,   \
                             ta1_z_yyzz_xxxx_0,  \
                             ta1_z_yyzz_xxxx_1,  \
                             ta1_z_yyzz_xxxy_0,  \
                             ta1_z_yyzz_xxxy_1,  \
                             ta1_z_yyzz_xxxz_0,  \
                             ta1_z_yyzz_xxxz_1,  \
                             ta1_z_yyzz_xxy_0,   \
                             ta1_z_yyzz_xxy_1,   \
                             ta1_z_yyzz_xxyy_0,  \
                             ta1_z_yyzz_xxyy_1,  \
                             ta1_z_yyzz_xxyz_0,  \
                             ta1_z_yyzz_xxyz_1,  \
                             ta1_z_yyzz_xxz_0,   \
                             ta1_z_yyzz_xxz_1,   \
                             ta1_z_yyzz_xxzz_0,  \
                             ta1_z_yyzz_xxzz_1,  \
                             ta1_z_yyzz_xyy_0,   \
                             ta1_z_yyzz_xyy_1,   \
                             ta1_z_yyzz_xyyy_0,  \
                             ta1_z_yyzz_xyyy_1,  \
                             ta1_z_yyzz_xyyz_0,  \
                             ta1_z_yyzz_xyyz_1,  \
                             ta1_z_yyzz_xyz_0,   \
                             ta1_z_yyzz_xyz_1,   \
                             ta1_z_yyzz_xyzz_0,  \
                             ta1_z_yyzz_xyzz_1,  \
                             ta1_z_yyzz_xzz_0,   \
                             ta1_z_yyzz_xzz_1,   \
                             ta1_z_yyzz_xzzz_0,  \
                             ta1_z_yyzz_xzzz_1,  \
                             ta1_z_yyzz_yyy_0,   \
                             ta1_z_yyzz_yyy_1,   \
                             ta1_z_yyzz_yyyy_0,  \
                             ta1_z_yyzz_yyyy_1,  \
                             ta1_z_yyzz_yyyz_0,  \
                             ta1_z_yyzz_yyyz_1,  \
                             ta1_z_yyzz_yyz_0,   \
                             ta1_z_yyzz_yyz_1,   \
                             ta1_z_yyzz_yyzz_0,  \
                             ta1_z_yyzz_yyzz_1,  \
                             ta1_z_yyzz_yzz_0,   \
                             ta1_z_yyzz_yzz_1,   \
                             ta1_z_yyzz_yzzz_0,  \
                             ta1_z_yyzz_yzzz_1,  \
                             ta1_z_yyzz_zzz_0,   \
                             ta1_z_yyzz_zzz_1,   \
                             ta1_z_yyzz_zzzz_0,  \
                             ta1_z_yyzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzz_xxxx_0[i] =
            4.0 * ta1_z_yyzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyzz_xxx_1[i] * fe_0 + ta1_z_yyzz_xxxx_0[i] * pa_x[i] - ta1_z_yyzz_xxxx_1[i] * pc_x[i];

        ta1_z_xyyzz_xxxy_0[i] =
            3.0 * ta1_z_yyzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxy_1[i] * fe_0 + ta1_z_yyzz_xxxy_0[i] * pa_x[i] - ta1_z_yyzz_xxxy_1[i] * pc_x[i];

        ta1_z_xyyzz_xxxz_0[i] =
            3.0 * ta1_z_yyzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxz_1[i] * fe_0 + ta1_z_yyzz_xxxz_0[i] * pa_x[i] - ta1_z_yyzz_xxxz_1[i] * pc_x[i];

        ta1_z_xyyzz_xxyy_0[i] =
            2.0 * ta1_z_yyzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyzz_xyy_1[i] * fe_0 + ta1_z_yyzz_xxyy_0[i] * pa_x[i] - ta1_z_yyzz_xxyy_1[i] * pc_x[i];

        ta1_z_xyyzz_xxyz_0[i] =
            2.0 * ta1_z_yyzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyzz_xyz_1[i] * fe_0 + ta1_z_yyzz_xxyz_0[i] * pa_x[i] - ta1_z_yyzz_xxyz_1[i] * pc_x[i];

        ta1_z_xyyzz_xxzz_0[i] =
            2.0 * ta1_z_yyzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyzz_xzz_1[i] * fe_0 + ta1_z_yyzz_xxzz_0[i] * pa_x[i] - ta1_z_yyzz_xxzz_1[i] * pc_x[i];

        ta1_z_xyyzz_xyyy_0[i] =
            ta1_z_yyzz_yyy_0[i] * fe_0 - ta1_z_yyzz_yyy_1[i] * fe_0 + ta1_z_yyzz_xyyy_0[i] * pa_x[i] - ta1_z_yyzz_xyyy_1[i] * pc_x[i];

        ta1_z_xyyzz_xyyz_0[i] =
            ta1_z_yyzz_yyz_0[i] * fe_0 - ta1_z_yyzz_yyz_1[i] * fe_0 + ta1_z_yyzz_xyyz_0[i] * pa_x[i] - ta1_z_yyzz_xyyz_1[i] * pc_x[i];

        ta1_z_xyyzz_xyzz_0[i] =
            ta1_z_yyzz_yzz_0[i] * fe_0 - ta1_z_yyzz_yzz_1[i] * fe_0 + ta1_z_yyzz_xyzz_0[i] * pa_x[i] - ta1_z_yyzz_xyzz_1[i] * pc_x[i];

        ta1_z_xyyzz_xzzz_0[i] =
            ta1_z_yyzz_zzz_0[i] * fe_0 - ta1_z_yyzz_zzz_1[i] * fe_0 + ta1_z_yyzz_xzzz_0[i] * pa_x[i] - ta1_z_yyzz_xzzz_1[i] * pc_x[i];

        ta1_z_xyyzz_yyyy_0[i] = ta1_z_yyzz_yyyy_0[i] * pa_x[i] - ta1_z_yyzz_yyyy_1[i] * pc_x[i];

        ta1_z_xyyzz_yyyz_0[i] = ta1_z_yyzz_yyyz_0[i] * pa_x[i] - ta1_z_yyzz_yyyz_1[i] * pc_x[i];

        ta1_z_xyyzz_yyzz_0[i] = ta1_z_yyzz_yyzz_0[i] * pa_x[i] - ta1_z_yyzz_yyzz_1[i] * pc_x[i];

        ta1_z_xyyzz_yzzz_0[i] = ta1_z_yyzz_yzzz_0[i] * pa_x[i] - ta1_z_yyzz_yzzz_1[i] * pc_x[i];

        ta1_z_xyyzz_zzzz_0[i] = ta1_z_yyzz_zzzz_0[i] * pa_x[i] - ta1_z_yyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 825-840 components of targeted buffer : HG

    auto ta1_z_xyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 825);

    auto ta1_z_xyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 826);

    auto ta1_z_xyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 827);

    auto ta1_z_xyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 828);

    auto ta1_z_xyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 829);

    auto ta1_z_xyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 830);

    auto ta1_z_xyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 831);

    auto ta1_z_xyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 832);

    auto ta1_z_xyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 833);

    auto ta1_z_xyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 834);

    auto ta1_z_xyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 835);

    auto ta1_z_xyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 836);

    auto ta1_z_xyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 837);

    auto ta1_z_xyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 838);

    auto ta1_z_xyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 839);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xyzzz_xxxx_0, \
                             ta1_z_xyzzz_xxxy_0, \
                             ta1_z_xyzzz_xxxz_0, \
                             ta1_z_xyzzz_xxyy_0, \
                             ta1_z_xyzzz_xxyz_0, \
                             ta1_z_xyzzz_xxzz_0, \
                             ta1_z_xyzzz_xyyy_0, \
                             ta1_z_xyzzz_xyyz_0, \
                             ta1_z_xyzzz_xyzz_0, \
                             ta1_z_xyzzz_xzzz_0, \
                             ta1_z_xyzzz_yyyy_0, \
                             ta1_z_xyzzz_yyyz_0, \
                             ta1_z_xyzzz_yyzz_0, \
                             ta1_z_xyzzz_yzzz_0, \
                             ta1_z_xyzzz_zzzz_0, \
                             ta1_z_xzzz_xxxx_0,  \
                             ta1_z_xzzz_xxxx_1,  \
                             ta1_z_xzzz_xxxz_0,  \
                             ta1_z_xzzz_xxxz_1,  \
                             ta1_z_xzzz_xxzz_0,  \
                             ta1_z_xzzz_xxzz_1,  \
                             ta1_z_xzzz_xzzz_0,  \
                             ta1_z_xzzz_xzzz_1,  \
                             ta1_z_yzzz_xxxy_0,  \
                             ta1_z_yzzz_xxxy_1,  \
                             ta1_z_yzzz_xxy_0,   \
                             ta1_z_yzzz_xxy_1,   \
                             ta1_z_yzzz_xxyy_0,  \
                             ta1_z_yzzz_xxyy_1,  \
                             ta1_z_yzzz_xxyz_0,  \
                             ta1_z_yzzz_xxyz_1,  \
                             ta1_z_yzzz_xyy_0,   \
                             ta1_z_yzzz_xyy_1,   \
                             ta1_z_yzzz_xyyy_0,  \
                             ta1_z_yzzz_xyyy_1,  \
                             ta1_z_yzzz_xyyz_0,  \
                             ta1_z_yzzz_xyyz_1,  \
                             ta1_z_yzzz_xyz_0,   \
                             ta1_z_yzzz_xyz_1,   \
                             ta1_z_yzzz_xyzz_0,  \
                             ta1_z_yzzz_xyzz_1,  \
                             ta1_z_yzzz_yyy_0,   \
                             ta1_z_yzzz_yyy_1,   \
                             ta1_z_yzzz_yyyy_0,  \
                             ta1_z_yzzz_yyyy_1,  \
                             ta1_z_yzzz_yyyz_0,  \
                             ta1_z_yzzz_yyyz_1,  \
                             ta1_z_yzzz_yyz_0,   \
                             ta1_z_yzzz_yyz_1,   \
                             ta1_z_yzzz_yyzz_0,  \
                             ta1_z_yzzz_yyzz_1,  \
                             ta1_z_yzzz_yzz_0,   \
                             ta1_z_yzzz_yzz_1,   \
                             ta1_z_yzzz_yzzz_0,  \
                             ta1_z_yzzz_yzzz_1,  \
                             ta1_z_yzzz_zzzz_0,  \
                             ta1_z_yzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzzz_xxxx_0[i] = ta1_z_xzzz_xxxx_0[i] * pa_y[i] - ta1_z_xzzz_xxxx_1[i] * pc_y[i];

        ta1_z_xyzzz_xxxy_0[i] =
            3.0 * ta1_z_yzzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yzzz_xxy_1[i] * fe_0 + ta1_z_yzzz_xxxy_0[i] * pa_x[i] - ta1_z_yzzz_xxxy_1[i] * pc_x[i];

        ta1_z_xyzzz_xxxz_0[i] = ta1_z_xzzz_xxxz_0[i] * pa_y[i] - ta1_z_xzzz_xxxz_1[i] * pc_y[i];

        ta1_z_xyzzz_xxyy_0[i] =
            2.0 * ta1_z_yzzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xyy_1[i] * fe_0 + ta1_z_yzzz_xxyy_0[i] * pa_x[i] - ta1_z_yzzz_xxyy_1[i] * pc_x[i];

        ta1_z_xyzzz_xxyz_0[i] =
            2.0 * ta1_z_yzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xyz_1[i] * fe_0 + ta1_z_yzzz_xxyz_0[i] * pa_x[i] - ta1_z_yzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xyzzz_xxzz_0[i] = ta1_z_xzzz_xxzz_0[i] * pa_y[i] - ta1_z_xzzz_xxzz_1[i] * pc_y[i];

        ta1_z_xyzzz_xyyy_0[i] =
            ta1_z_yzzz_yyy_0[i] * fe_0 - ta1_z_yzzz_yyy_1[i] * fe_0 + ta1_z_yzzz_xyyy_0[i] * pa_x[i] - ta1_z_yzzz_xyyy_1[i] * pc_x[i];

        ta1_z_xyzzz_xyyz_0[i] =
            ta1_z_yzzz_yyz_0[i] * fe_0 - ta1_z_yzzz_yyz_1[i] * fe_0 + ta1_z_yzzz_xyyz_0[i] * pa_x[i] - ta1_z_yzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xyzzz_xyzz_0[i] =
            ta1_z_yzzz_yzz_0[i] * fe_0 - ta1_z_yzzz_yzz_1[i] * fe_0 + ta1_z_yzzz_xyzz_0[i] * pa_x[i] - ta1_z_yzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xyzzz_xzzz_0[i] = ta1_z_xzzz_xzzz_0[i] * pa_y[i] - ta1_z_xzzz_xzzz_1[i] * pc_y[i];

        ta1_z_xyzzz_yyyy_0[i] = ta1_z_yzzz_yyyy_0[i] * pa_x[i] - ta1_z_yzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xyzzz_yyyz_0[i] = ta1_z_yzzz_yyyz_0[i] * pa_x[i] - ta1_z_yzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xyzzz_yyzz_0[i] = ta1_z_yzzz_yyzz_0[i] * pa_x[i] - ta1_z_yzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xyzzz_yzzz_0[i] = ta1_z_yzzz_yzzz_0[i] * pa_x[i] - ta1_z_yzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xyzzz_zzzz_0[i] = ta1_z_yzzz_zzzz_0[i] * pa_x[i] - ta1_z_yzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 840-855 components of targeted buffer : HG

    auto ta1_z_xzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 840);

    auto ta1_z_xzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 841);

    auto ta1_z_xzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 842);

    auto ta1_z_xzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 843);

    auto ta1_z_xzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 844);

    auto ta1_z_xzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 845);

    auto ta1_z_xzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 846);

    auto ta1_z_xzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 847);

    auto ta1_z_xzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 848);

    auto ta1_z_xzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 849);

    auto ta1_z_xzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 850);

    auto ta1_z_xzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 851);

    auto ta1_z_xzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 852);

    auto ta1_z_xzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 853);

    auto ta1_z_xzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 854);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xzzzz_xxxx_0, \
                             ta1_z_xzzzz_xxxy_0, \
                             ta1_z_xzzzz_xxxz_0, \
                             ta1_z_xzzzz_xxyy_0, \
                             ta1_z_xzzzz_xxyz_0, \
                             ta1_z_xzzzz_xxzz_0, \
                             ta1_z_xzzzz_xyyy_0, \
                             ta1_z_xzzzz_xyyz_0, \
                             ta1_z_xzzzz_xyzz_0, \
                             ta1_z_xzzzz_xzzz_0, \
                             ta1_z_xzzzz_yyyy_0, \
                             ta1_z_xzzzz_yyyz_0, \
                             ta1_z_xzzzz_yyzz_0, \
                             ta1_z_xzzzz_yzzz_0, \
                             ta1_z_xzzzz_zzzz_0, \
                             ta1_z_zzzz_xxx_0,   \
                             ta1_z_zzzz_xxx_1,   \
                             ta1_z_zzzz_xxxx_0,  \
                             ta1_z_zzzz_xxxx_1,  \
                             ta1_z_zzzz_xxxy_0,  \
                             ta1_z_zzzz_xxxy_1,  \
                             ta1_z_zzzz_xxxz_0,  \
                             ta1_z_zzzz_xxxz_1,  \
                             ta1_z_zzzz_xxy_0,   \
                             ta1_z_zzzz_xxy_1,   \
                             ta1_z_zzzz_xxyy_0,  \
                             ta1_z_zzzz_xxyy_1,  \
                             ta1_z_zzzz_xxyz_0,  \
                             ta1_z_zzzz_xxyz_1,  \
                             ta1_z_zzzz_xxz_0,   \
                             ta1_z_zzzz_xxz_1,   \
                             ta1_z_zzzz_xxzz_0,  \
                             ta1_z_zzzz_xxzz_1,  \
                             ta1_z_zzzz_xyy_0,   \
                             ta1_z_zzzz_xyy_1,   \
                             ta1_z_zzzz_xyyy_0,  \
                             ta1_z_zzzz_xyyy_1,  \
                             ta1_z_zzzz_xyyz_0,  \
                             ta1_z_zzzz_xyyz_1,  \
                             ta1_z_zzzz_xyz_0,   \
                             ta1_z_zzzz_xyz_1,   \
                             ta1_z_zzzz_xyzz_0,  \
                             ta1_z_zzzz_xyzz_1,  \
                             ta1_z_zzzz_xzz_0,   \
                             ta1_z_zzzz_xzz_1,   \
                             ta1_z_zzzz_xzzz_0,  \
                             ta1_z_zzzz_xzzz_1,  \
                             ta1_z_zzzz_yyy_0,   \
                             ta1_z_zzzz_yyy_1,   \
                             ta1_z_zzzz_yyyy_0,  \
                             ta1_z_zzzz_yyyy_1,  \
                             ta1_z_zzzz_yyyz_0,  \
                             ta1_z_zzzz_yyyz_1,  \
                             ta1_z_zzzz_yyz_0,   \
                             ta1_z_zzzz_yyz_1,   \
                             ta1_z_zzzz_yyzz_0,  \
                             ta1_z_zzzz_yyzz_1,  \
                             ta1_z_zzzz_yzz_0,   \
                             ta1_z_zzzz_yzz_1,   \
                             ta1_z_zzzz_yzzz_0,  \
                             ta1_z_zzzz_yzzz_1,  \
                             ta1_z_zzzz_zzz_0,   \
                             ta1_z_zzzz_zzz_1,   \
                             ta1_z_zzzz_zzzz_0,  \
                             ta1_z_zzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzz_xxxx_0[i] =
            4.0 * ta1_z_zzzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_zzzz_xxx_1[i] * fe_0 + ta1_z_zzzz_xxxx_0[i] * pa_x[i] - ta1_z_zzzz_xxxx_1[i] * pc_x[i];

        ta1_z_xzzzz_xxxy_0[i] =
            3.0 * ta1_z_zzzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_zzzz_xxy_1[i] * fe_0 + ta1_z_zzzz_xxxy_0[i] * pa_x[i] - ta1_z_zzzz_xxxy_1[i] * pc_x[i];

        ta1_z_xzzzz_xxxz_0[i] =
            3.0 * ta1_z_zzzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_zzzz_xxz_1[i] * fe_0 + ta1_z_zzzz_xxxz_0[i] * pa_x[i] - ta1_z_zzzz_xxxz_1[i] * pc_x[i];

        ta1_z_xzzzz_xxyy_0[i] =
            2.0 * ta1_z_zzzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xyy_1[i] * fe_0 + ta1_z_zzzz_xxyy_0[i] * pa_x[i] - ta1_z_zzzz_xxyy_1[i] * pc_x[i];

        ta1_z_xzzzz_xxyz_0[i] =
            2.0 * ta1_z_zzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xyz_1[i] * fe_0 + ta1_z_zzzz_xxyz_0[i] * pa_x[i] - ta1_z_zzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xzzzz_xxzz_0[i] =
            2.0 * ta1_z_zzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xzz_1[i] * fe_0 + ta1_z_zzzz_xxzz_0[i] * pa_x[i] - ta1_z_zzzz_xxzz_1[i] * pc_x[i];

        ta1_z_xzzzz_xyyy_0[i] =
            ta1_z_zzzz_yyy_0[i] * fe_0 - ta1_z_zzzz_yyy_1[i] * fe_0 + ta1_z_zzzz_xyyy_0[i] * pa_x[i] - ta1_z_zzzz_xyyy_1[i] * pc_x[i];

        ta1_z_xzzzz_xyyz_0[i] =
            ta1_z_zzzz_yyz_0[i] * fe_0 - ta1_z_zzzz_yyz_1[i] * fe_0 + ta1_z_zzzz_xyyz_0[i] * pa_x[i] - ta1_z_zzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xzzzz_xyzz_0[i] =
            ta1_z_zzzz_yzz_0[i] * fe_0 - ta1_z_zzzz_yzz_1[i] * fe_0 + ta1_z_zzzz_xyzz_0[i] * pa_x[i] - ta1_z_zzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xzzzz_xzzz_0[i] =
            ta1_z_zzzz_zzz_0[i] * fe_0 - ta1_z_zzzz_zzz_1[i] * fe_0 + ta1_z_zzzz_xzzz_0[i] * pa_x[i] - ta1_z_zzzz_xzzz_1[i] * pc_x[i];

        ta1_z_xzzzz_yyyy_0[i] = ta1_z_zzzz_yyyy_0[i] * pa_x[i] - ta1_z_zzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xzzzz_yyyz_0[i] = ta1_z_zzzz_yyyz_0[i] * pa_x[i] - ta1_z_zzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xzzzz_yyzz_0[i] = ta1_z_zzzz_yyzz_0[i] * pa_x[i] - ta1_z_zzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xzzzz_yzzz_0[i] = ta1_z_zzzz_yzzz_0[i] * pa_x[i] - ta1_z_zzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xzzzz_zzzz_0[i] = ta1_z_zzzz_zzzz_0[i] * pa_x[i] - ta1_z_zzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 855-870 components of targeted buffer : HG

    auto ta1_z_yyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 855);

    auto ta1_z_yyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 856);

    auto ta1_z_yyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 857);

    auto ta1_z_yyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 858);

    auto ta1_z_yyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 859);

    auto ta1_z_yyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 860);

    auto ta1_z_yyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 861);

    auto ta1_z_yyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 862);

    auto ta1_z_yyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 863);

    auto ta1_z_yyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 864);

    auto ta1_z_yyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 865);

    auto ta1_z_yyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 866);

    auto ta1_z_yyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 867);

    auto ta1_z_yyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 868);

    auto ta1_z_yyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 869);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yyy_xxxx_0,   \
                             ta1_z_yyy_xxxx_1,   \
                             ta1_z_yyy_xxxy_0,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxxz_0,   \
                             ta1_z_yyy_xxxz_1,   \
                             ta1_z_yyy_xxyy_0,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyz_0,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xxzz_0,   \
                             ta1_z_yyy_xxzz_1,   \
                             ta1_z_yyy_xyyy_0,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyz_0,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyzz_0,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_xzzz_0,   \
                             ta1_z_yyy_xzzz_1,   \
                             ta1_z_yyy_yyyy_0,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyz_0,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyzz_0,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yzzz_0,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_zzzz_0,   \
                             ta1_z_yyy_zzzz_1,   \
                             ta1_z_yyyy_xxx_0,   \
                             ta1_z_yyyy_xxx_1,   \
                             ta1_z_yyyy_xxxx_0,  \
                             ta1_z_yyyy_xxxx_1,  \
                             ta1_z_yyyy_xxxy_0,  \
                             ta1_z_yyyy_xxxy_1,  \
                             ta1_z_yyyy_xxxz_0,  \
                             ta1_z_yyyy_xxxz_1,  \
                             ta1_z_yyyy_xxy_0,   \
                             ta1_z_yyyy_xxy_1,   \
                             ta1_z_yyyy_xxyy_0,  \
                             ta1_z_yyyy_xxyy_1,  \
                             ta1_z_yyyy_xxyz_0,  \
                             ta1_z_yyyy_xxyz_1,  \
                             ta1_z_yyyy_xxz_0,   \
                             ta1_z_yyyy_xxz_1,   \
                             ta1_z_yyyy_xxzz_0,  \
                             ta1_z_yyyy_xxzz_1,  \
                             ta1_z_yyyy_xyy_0,   \
                             ta1_z_yyyy_xyy_1,   \
                             ta1_z_yyyy_xyyy_0,  \
                             ta1_z_yyyy_xyyy_1,  \
                             ta1_z_yyyy_xyyz_0,  \
                             ta1_z_yyyy_xyyz_1,  \
                             ta1_z_yyyy_xyz_0,   \
                             ta1_z_yyyy_xyz_1,   \
                             ta1_z_yyyy_xyzz_0,  \
                             ta1_z_yyyy_xyzz_1,  \
                             ta1_z_yyyy_xzz_0,   \
                             ta1_z_yyyy_xzz_1,   \
                             ta1_z_yyyy_xzzz_0,  \
                             ta1_z_yyyy_xzzz_1,  \
                             ta1_z_yyyy_yyy_0,   \
                             ta1_z_yyyy_yyy_1,   \
                             ta1_z_yyyy_yyyy_0,  \
                             ta1_z_yyyy_yyyy_1,  \
                             ta1_z_yyyy_yyyz_0,  \
                             ta1_z_yyyy_yyyz_1,  \
                             ta1_z_yyyy_yyz_0,   \
                             ta1_z_yyyy_yyz_1,   \
                             ta1_z_yyyy_yyzz_0,  \
                             ta1_z_yyyy_yyzz_1,  \
                             ta1_z_yyyy_yzz_0,   \
                             ta1_z_yyyy_yzz_1,   \
                             ta1_z_yyyy_yzzz_0,  \
                             ta1_z_yyyy_yzzz_1,  \
                             ta1_z_yyyy_zzz_0,   \
                             ta1_z_yyyy_zzz_1,   \
                             ta1_z_yyyy_zzzz_0,  \
                             ta1_z_yyyy_zzzz_1,  \
                             ta1_z_yyyyy_xxxx_0, \
                             ta1_z_yyyyy_xxxy_0, \
                             ta1_z_yyyyy_xxxz_0, \
                             ta1_z_yyyyy_xxyy_0, \
                             ta1_z_yyyyy_xxyz_0, \
                             ta1_z_yyyyy_xxzz_0, \
                             ta1_z_yyyyy_xyyy_0, \
                             ta1_z_yyyyy_xyyz_0, \
                             ta1_z_yyyyy_xyzz_0, \
                             ta1_z_yyyyy_xzzz_0, \
                             ta1_z_yyyyy_yyyy_0, \
                             ta1_z_yyyyy_yyyz_0, \
                             ta1_z_yyyyy_yyzz_0, \
                             ta1_z_yyyyy_yzzz_0, \
                             ta1_z_yyyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyy_xxxx_0[i] =
            4.0 * ta1_z_yyy_xxxx_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxx_1[i] * fe_0 + ta1_z_yyyy_xxxx_0[i] * pa_y[i] - ta1_z_yyyy_xxxx_1[i] * pc_y[i];

        ta1_z_yyyyy_xxxy_0[i] = 4.0 * ta1_z_yyy_xxxy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxy_1[i] * fe_0 + ta1_z_yyyy_xxx_0[i] * fe_0 -
                                ta1_z_yyyy_xxx_1[i] * fe_0 + ta1_z_yyyy_xxxy_0[i] * pa_y[i] - ta1_z_yyyy_xxxy_1[i] * pc_y[i];

        ta1_z_yyyyy_xxxz_0[i] =
            4.0 * ta1_z_yyy_xxxz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxz_1[i] * fe_0 + ta1_z_yyyy_xxxz_0[i] * pa_y[i] - ta1_z_yyyy_xxxz_1[i] * pc_y[i];

        ta1_z_yyyyy_xxyy_0[i] = 4.0 * ta1_z_yyy_xxyy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_yyyy_xxy_0[i] * fe_0 -
                                2.0 * ta1_z_yyyy_xxy_1[i] * fe_0 + ta1_z_yyyy_xxyy_0[i] * pa_y[i] - ta1_z_yyyy_xxyy_1[i] * pc_y[i];

        ta1_z_yyyyy_xxyz_0[i] = 4.0 * ta1_z_yyy_xxyz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxyz_1[i] * fe_0 + ta1_z_yyyy_xxz_0[i] * fe_0 -
                                ta1_z_yyyy_xxz_1[i] * fe_0 + ta1_z_yyyy_xxyz_0[i] * pa_y[i] - ta1_z_yyyy_xxyz_1[i] * pc_y[i];

        ta1_z_yyyyy_xxzz_0[i] =
            4.0 * ta1_z_yyy_xxzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxzz_1[i] * fe_0 + ta1_z_yyyy_xxzz_0[i] * pa_y[i] - ta1_z_yyyy_xxzz_1[i] * pc_y[i];

        ta1_z_yyyyy_xyyy_0[i] = 4.0 * ta1_z_yyy_xyyy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xyyy_1[i] * fe_0 + 3.0 * ta1_z_yyyy_xyy_0[i] * fe_0 -
                                3.0 * ta1_z_yyyy_xyy_1[i] * fe_0 + ta1_z_yyyy_xyyy_0[i] * pa_y[i] - ta1_z_yyyy_xyyy_1[i] * pc_y[i];

        ta1_z_yyyyy_xyyz_0[i] = 4.0 * ta1_z_yyy_xyyz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yyyy_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_yyyy_xyz_1[i] * fe_0 + ta1_z_yyyy_xyyz_0[i] * pa_y[i] - ta1_z_yyyy_xyyz_1[i] * pc_y[i];

        ta1_z_yyyyy_xyzz_0[i] = 4.0 * ta1_z_yyy_xyzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xyzz_1[i] * fe_0 + ta1_z_yyyy_xzz_0[i] * fe_0 -
                                ta1_z_yyyy_xzz_1[i] * fe_0 + ta1_z_yyyy_xyzz_0[i] * pa_y[i] - ta1_z_yyyy_xyzz_1[i] * pc_y[i];

        ta1_z_yyyyy_xzzz_0[i] =
            4.0 * ta1_z_yyy_xzzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xzzz_1[i] * fe_0 + ta1_z_yyyy_xzzz_0[i] * pa_y[i] - ta1_z_yyyy_xzzz_1[i] * pc_y[i];

        ta1_z_yyyyy_yyyy_0[i] = 4.0 * ta1_z_yyy_yyyy_0[i] * fe_0 - 4.0 * ta1_z_yyy_yyyy_1[i] * fe_0 + 4.0 * ta1_z_yyyy_yyy_0[i] * fe_0 -
                                4.0 * ta1_z_yyyy_yyy_1[i] * fe_0 + ta1_z_yyyy_yyyy_0[i] * pa_y[i] - ta1_z_yyyy_yyyy_1[i] * pc_y[i];

        ta1_z_yyyyy_yyyz_0[i] = 4.0 * ta1_z_yyy_yyyz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yyyy_yyz_0[i] * fe_0 -
                                3.0 * ta1_z_yyyy_yyz_1[i] * fe_0 + ta1_z_yyyy_yyyz_0[i] * pa_y[i] - ta1_z_yyyy_yyyz_1[i] * pc_y[i];

        ta1_z_yyyyy_yyzz_0[i] = 4.0 * ta1_z_yyy_yyzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yyyy_yzz_0[i] * fe_0 -
                                2.0 * ta1_z_yyyy_yzz_1[i] * fe_0 + ta1_z_yyyy_yyzz_0[i] * pa_y[i] - ta1_z_yyyy_yyzz_1[i] * pc_y[i];

        ta1_z_yyyyy_yzzz_0[i] = 4.0 * ta1_z_yyy_yzzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yzzz_1[i] * fe_0 + ta1_z_yyyy_zzz_0[i] * fe_0 -
                                ta1_z_yyyy_zzz_1[i] * fe_0 + ta1_z_yyyy_yzzz_0[i] * pa_y[i] - ta1_z_yyyy_yzzz_1[i] * pc_y[i];

        ta1_z_yyyyy_zzzz_0[i] =
            4.0 * ta1_z_yyy_zzzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_zzzz_1[i] * fe_0 + ta1_z_yyyy_zzzz_0[i] * pa_y[i] - ta1_z_yyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 870-885 components of targeted buffer : HG

    auto ta1_z_yyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 870);

    auto ta1_z_yyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 871);

    auto ta1_z_yyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 872);

    auto ta1_z_yyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 873);

    auto ta1_z_yyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 874);

    auto ta1_z_yyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 875);

    auto ta1_z_yyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 876);

    auto ta1_z_yyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 877);

    auto ta1_z_yyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 878);

    auto ta1_z_yyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 879);

    auto ta1_z_yyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 880);

    auto ta1_z_yyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 881);

    auto ta1_z_yyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 882);

    auto ta1_z_yyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 883);

    auto ta1_z_yyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 884);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyyy_xxxx_0,  \
                             ta1_z_yyyy_xxxx_1,  \
                             ta1_z_yyyy_xxxy_0,  \
                             ta1_z_yyyy_xxxy_1,  \
                             ta1_z_yyyy_xxy_0,   \
                             ta1_z_yyyy_xxy_1,   \
                             ta1_z_yyyy_xxyy_0,  \
                             ta1_z_yyyy_xxyy_1,  \
                             ta1_z_yyyy_xxyz_0,  \
                             ta1_z_yyyy_xxyz_1,  \
                             ta1_z_yyyy_xyy_0,   \
                             ta1_z_yyyy_xyy_1,   \
                             ta1_z_yyyy_xyyy_0,  \
                             ta1_z_yyyy_xyyy_1,  \
                             ta1_z_yyyy_xyyz_0,  \
                             ta1_z_yyyy_xyyz_1,  \
                             ta1_z_yyyy_xyz_0,   \
                             ta1_z_yyyy_xyz_1,   \
                             ta1_z_yyyy_xyzz_0,  \
                             ta1_z_yyyy_xyzz_1,  \
                             ta1_z_yyyy_yyy_0,   \
                             ta1_z_yyyy_yyy_1,   \
                             ta1_z_yyyy_yyyy_0,  \
                             ta1_z_yyyy_yyyy_1,  \
                             ta1_z_yyyy_yyyz_0,  \
                             ta1_z_yyyy_yyyz_1,  \
                             ta1_z_yyyy_yyz_0,   \
                             ta1_z_yyyy_yyz_1,   \
                             ta1_z_yyyy_yyzz_0,  \
                             ta1_z_yyyy_yyzz_1,  \
                             ta1_z_yyyy_yzz_0,   \
                             ta1_z_yyyy_yzz_1,   \
                             ta1_z_yyyy_yzzz_0,  \
                             ta1_z_yyyy_yzzz_1,  \
                             ta1_z_yyyyz_xxxx_0, \
                             ta1_z_yyyyz_xxxy_0, \
                             ta1_z_yyyyz_xxxz_0, \
                             ta1_z_yyyyz_xxyy_0, \
                             ta1_z_yyyyz_xxyz_0, \
                             ta1_z_yyyyz_xxzz_0, \
                             ta1_z_yyyyz_xyyy_0, \
                             ta1_z_yyyyz_xyyz_0, \
                             ta1_z_yyyyz_xyzz_0, \
                             ta1_z_yyyyz_xzzz_0, \
                             ta1_z_yyyyz_yyyy_0, \
                             ta1_z_yyyyz_yyyz_0, \
                             ta1_z_yyyyz_yyzz_0, \
                             ta1_z_yyyyz_yzzz_0, \
                             ta1_z_yyyyz_zzzz_0, \
                             ta1_z_yyyz_xxxz_0,  \
                             ta1_z_yyyz_xxxz_1,  \
                             ta1_z_yyyz_xxzz_0,  \
                             ta1_z_yyyz_xxzz_1,  \
                             ta1_z_yyyz_xzzz_0,  \
                             ta1_z_yyyz_xzzz_1,  \
                             ta1_z_yyyz_zzzz_0,  \
                             ta1_z_yyyz_zzzz_1,  \
                             ta1_z_yyz_xxxz_0,   \
                             ta1_z_yyz_xxxz_1,   \
                             ta1_z_yyz_xxzz_0,   \
                             ta1_z_yyz_xxzz_1,   \
                             ta1_z_yyz_xzzz_0,   \
                             ta1_z_yyz_xzzz_1,   \
                             ta1_z_yyz_zzzz_0,   \
                             ta1_z_yyz_zzzz_1,   \
                             ta_yyyy_xxxx_1,     \
                             ta_yyyy_xxxy_1,     \
                             ta_yyyy_xxyy_1,     \
                             ta_yyyy_xxyz_1,     \
                             ta_yyyy_xyyy_1,     \
                             ta_yyyy_xyyz_1,     \
                             ta_yyyy_xyzz_1,     \
                             ta_yyyy_yyyy_1,     \
                             ta_yyyy_yyyz_1,     \
                             ta_yyyy_yyzz_1,     \
                             ta_yyyy_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyz_xxxx_0[i] = ta_yyyy_xxxx_1[i] + ta1_z_yyyy_xxxx_0[i] * pa_z[i] - ta1_z_yyyy_xxxx_1[i] * pc_z[i];

        ta1_z_yyyyz_xxxy_0[i] = ta_yyyy_xxxy_1[i] + ta1_z_yyyy_xxxy_0[i] * pa_z[i] - ta1_z_yyyy_xxxy_1[i] * pc_z[i];

        ta1_z_yyyyz_xxxz_0[i] =
            3.0 * ta1_z_yyz_xxxz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxxz_1[i] * fe_0 + ta1_z_yyyz_xxxz_0[i] * pa_y[i] - ta1_z_yyyz_xxxz_1[i] * pc_y[i];

        ta1_z_yyyyz_xxyy_0[i] = ta_yyyy_xxyy_1[i] + ta1_z_yyyy_xxyy_0[i] * pa_z[i] - ta1_z_yyyy_xxyy_1[i] * pc_z[i];

        ta1_z_yyyyz_xxyz_0[i] = ta1_z_yyyy_xxy_0[i] * fe_0 - ta1_z_yyyy_xxy_1[i] * fe_0 + ta_yyyy_xxyz_1[i] + ta1_z_yyyy_xxyz_0[i] * pa_z[i] -
                                ta1_z_yyyy_xxyz_1[i] * pc_z[i];

        ta1_z_yyyyz_xxzz_0[i] =
            3.0 * ta1_z_yyz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxzz_1[i] * fe_0 + ta1_z_yyyz_xxzz_0[i] * pa_y[i] - ta1_z_yyyz_xxzz_1[i] * pc_y[i];

        ta1_z_yyyyz_xyyy_0[i] = ta_yyyy_xyyy_1[i] + ta1_z_yyyy_xyyy_0[i] * pa_z[i] - ta1_z_yyyy_xyyy_1[i] * pc_z[i];

        ta1_z_yyyyz_xyyz_0[i] = ta1_z_yyyy_xyy_0[i] * fe_0 - ta1_z_yyyy_xyy_1[i] * fe_0 + ta_yyyy_xyyz_1[i] + ta1_z_yyyy_xyyz_0[i] * pa_z[i] -
                                ta1_z_yyyy_xyyz_1[i] * pc_z[i];

        ta1_z_yyyyz_xyzz_0[i] = 2.0 * ta1_z_yyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xyz_1[i] * fe_0 + ta_yyyy_xyzz_1[i] +
                                ta1_z_yyyy_xyzz_0[i] * pa_z[i] - ta1_z_yyyy_xyzz_1[i] * pc_z[i];

        ta1_z_yyyyz_xzzz_0[i] =
            3.0 * ta1_z_yyz_xzzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xzzz_1[i] * fe_0 + ta1_z_yyyz_xzzz_0[i] * pa_y[i] - ta1_z_yyyz_xzzz_1[i] * pc_y[i];

        ta1_z_yyyyz_yyyy_0[i] = ta_yyyy_yyyy_1[i] + ta1_z_yyyy_yyyy_0[i] * pa_z[i] - ta1_z_yyyy_yyyy_1[i] * pc_z[i];

        ta1_z_yyyyz_yyyz_0[i] = ta1_z_yyyy_yyy_0[i] * fe_0 - ta1_z_yyyy_yyy_1[i] * fe_0 + ta_yyyy_yyyz_1[i] + ta1_z_yyyy_yyyz_0[i] * pa_z[i] -
                                ta1_z_yyyy_yyyz_1[i] * pc_z[i];

        ta1_z_yyyyz_yyzz_0[i] = 2.0 * ta1_z_yyyy_yyz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_yyz_1[i] * fe_0 + ta_yyyy_yyzz_1[i] +
                                ta1_z_yyyy_yyzz_0[i] * pa_z[i] - ta1_z_yyyy_yyzz_1[i] * pc_z[i];

        ta1_z_yyyyz_yzzz_0[i] = 3.0 * ta1_z_yyyy_yzz_0[i] * fe_0 - 3.0 * ta1_z_yyyy_yzz_1[i] * fe_0 + ta_yyyy_yzzz_1[i] +
                                ta1_z_yyyy_yzzz_0[i] * pa_z[i] - ta1_z_yyyy_yzzz_1[i] * pc_z[i];

        ta1_z_yyyyz_zzzz_0[i] =
            3.0 * ta1_z_yyz_zzzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_zzzz_1[i] * fe_0 + ta1_z_yyyz_zzzz_0[i] * pa_y[i] - ta1_z_yyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 885-900 components of targeted buffer : HG

    auto ta1_z_yyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 885);

    auto ta1_z_yyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 886);

    auto ta1_z_yyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 887);

    auto ta1_z_yyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 888);

    auto ta1_z_yyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 889);

    auto ta1_z_yyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 890);

    auto ta1_z_yyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 891);

    auto ta1_z_yyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 892);

    auto ta1_z_yyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 893);

    auto ta1_z_yyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 894);

    auto ta1_z_yyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 895);

    auto ta1_z_yyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 896);

    auto ta1_z_yyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 897);

    auto ta1_z_yyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 898);

    auto ta1_z_yyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 899);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyy_xxxy_0,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxyy_0,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xyyy_0,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_yyyy_0,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyyz_xxxy_0,  \
                             ta1_z_yyyz_xxxy_1,  \
                             ta1_z_yyyz_xxyy_0,  \
                             ta1_z_yyyz_xxyy_1,  \
                             ta1_z_yyyz_xyyy_0,  \
                             ta1_z_yyyz_xyyy_1,  \
                             ta1_z_yyyz_yyyy_0,  \
                             ta1_z_yyyz_yyyy_1,  \
                             ta1_z_yyyzz_xxxx_0, \
                             ta1_z_yyyzz_xxxy_0, \
                             ta1_z_yyyzz_xxxz_0, \
                             ta1_z_yyyzz_xxyy_0, \
                             ta1_z_yyyzz_xxyz_0, \
                             ta1_z_yyyzz_xxzz_0, \
                             ta1_z_yyyzz_xyyy_0, \
                             ta1_z_yyyzz_xyyz_0, \
                             ta1_z_yyyzz_xyzz_0, \
                             ta1_z_yyyzz_xzzz_0, \
                             ta1_z_yyyzz_yyyy_0, \
                             ta1_z_yyyzz_yyyz_0, \
                             ta1_z_yyyzz_yyzz_0, \
                             ta1_z_yyyzz_yzzz_0, \
                             ta1_z_yyyzz_zzzz_0, \
                             ta1_z_yyzz_xxxx_0,  \
                             ta1_z_yyzz_xxxx_1,  \
                             ta1_z_yyzz_xxxz_0,  \
                             ta1_z_yyzz_xxxz_1,  \
                             ta1_z_yyzz_xxyz_0,  \
                             ta1_z_yyzz_xxyz_1,  \
                             ta1_z_yyzz_xxz_0,   \
                             ta1_z_yyzz_xxz_1,   \
                             ta1_z_yyzz_xxzz_0,  \
                             ta1_z_yyzz_xxzz_1,  \
                             ta1_z_yyzz_xyyz_0,  \
                             ta1_z_yyzz_xyyz_1,  \
                             ta1_z_yyzz_xyz_0,   \
                             ta1_z_yyzz_xyz_1,   \
                             ta1_z_yyzz_xyzz_0,  \
                             ta1_z_yyzz_xyzz_1,  \
                             ta1_z_yyzz_xzz_0,   \
                             ta1_z_yyzz_xzz_1,   \
                             ta1_z_yyzz_xzzz_0,  \
                             ta1_z_yyzz_xzzz_1,  \
                             ta1_z_yyzz_yyyz_0,  \
                             ta1_z_yyzz_yyyz_1,  \
                             ta1_z_yyzz_yyz_0,   \
                             ta1_z_yyzz_yyz_1,   \
                             ta1_z_yyzz_yyzz_0,  \
                             ta1_z_yyzz_yyzz_1,  \
                             ta1_z_yyzz_yzz_0,   \
                             ta1_z_yyzz_yzz_1,   \
                             ta1_z_yyzz_yzzz_0,  \
                             ta1_z_yyzz_yzzz_1,  \
                             ta1_z_yyzz_zzz_0,   \
                             ta1_z_yyzz_zzz_1,   \
                             ta1_z_yyzz_zzzz_0,  \
                             ta1_z_yyzz_zzzz_1,  \
                             ta1_z_yzz_xxxx_0,   \
                             ta1_z_yzz_xxxx_1,   \
                             ta1_z_yzz_xxxz_0,   \
                             ta1_z_yzz_xxxz_1,   \
                             ta1_z_yzz_xxyz_0,   \
                             ta1_z_yzz_xxyz_1,   \
                             ta1_z_yzz_xxzz_0,   \
                             ta1_z_yzz_xxzz_1,   \
                             ta1_z_yzz_xyyz_0,   \
                             ta1_z_yzz_xyyz_1,   \
                             ta1_z_yzz_xyzz_0,   \
                             ta1_z_yzz_xyzz_1,   \
                             ta1_z_yzz_xzzz_0,   \
                             ta1_z_yzz_xzzz_1,   \
                             ta1_z_yzz_yyyz_0,   \
                             ta1_z_yzz_yyyz_1,   \
                             ta1_z_yzz_yyzz_0,   \
                             ta1_z_yzz_yyzz_1,   \
                             ta1_z_yzz_yzzz_0,   \
                             ta1_z_yzz_yzzz_1,   \
                             ta1_z_yzz_zzzz_0,   \
                             ta1_z_yzz_zzzz_1,   \
                             ta_yyyz_xxxy_1,     \
                             ta_yyyz_xxyy_1,     \
                             ta_yyyz_xyyy_1,     \
                             ta_yyyz_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzz_xxxx_0[i] =
            2.0 * ta1_z_yzz_xxxx_0[i] * fe_0 - 2.0 * ta1_z_yzz_xxxx_1[i] * fe_0 + ta1_z_yyzz_xxxx_0[i] * pa_y[i] - ta1_z_yyzz_xxxx_1[i] * pc_y[i];

        ta1_z_yyyzz_xxxy_0[i] = ta1_z_yyy_xxxy_0[i] * fe_0 - ta1_z_yyy_xxxy_1[i] * fe_0 + ta_yyyz_xxxy_1[i] + ta1_z_yyyz_xxxy_0[i] * pa_z[i] -
                                ta1_z_yyyz_xxxy_1[i] * pc_z[i];

        ta1_z_yyyzz_xxxz_0[i] =
            2.0 * ta1_z_yzz_xxxz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xxxz_1[i] * fe_0 + ta1_z_yyzz_xxxz_0[i] * pa_y[i] - ta1_z_yyzz_xxxz_1[i] * pc_y[i];

        ta1_z_yyyzz_xxyy_0[i] = ta1_z_yyy_xxyy_0[i] * fe_0 - ta1_z_yyy_xxyy_1[i] * fe_0 + ta_yyyz_xxyy_1[i] + ta1_z_yyyz_xxyy_0[i] * pa_z[i] -
                                ta1_z_yyyz_xxyy_1[i] * pc_z[i];

        ta1_z_yyyzz_xxyz_0[i] = 2.0 * ta1_z_yzz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xxyz_1[i] * fe_0 + ta1_z_yyzz_xxz_0[i] * fe_0 -
                                ta1_z_yyzz_xxz_1[i] * fe_0 + ta1_z_yyzz_xxyz_0[i] * pa_y[i] - ta1_z_yyzz_xxyz_1[i] * pc_y[i];

        ta1_z_yyyzz_xxzz_0[i] =
            2.0 * ta1_z_yzz_xxzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xxzz_1[i] * fe_0 + ta1_z_yyzz_xxzz_0[i] * pa_y[i] - ta1_z_yyzz_xxzz_1[i] * pc_y[i];

        ta1_z_yyyzz_xyyy_0[i] = ta1_z_yyy_xyyy_0[i] * fe_0 - ta1_z_yyy_xyyy_1[i] * fe_0 + ta_yyyz_xyyy_1[i] + ta1_z_yyyz_xyyy_0[i] * pa_z[i] -
                                ta1_z_yyyz_xyyy_1[i] * pc_z[i];

        ta1_z_yyyzz_xyyz_0[i] = 2.0 * ta1_z_yzz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yyzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_yyzz_xyz_1[i] * fe_0 + ta1_z_yyzz_xyyz_0[i] * pa_y[i] - ta1_z_yyzz_xyyz_1[i] * pc_y[i];

        ta1_z_yyyzz_xyzz_0[i] = 2.0 * ta1_z_yzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyzz_1[i] * fe_0 + ta1_z_yyzz_xzz_0[i] * fe_0 -
                                ta1_z_yyzz_xzz_1[i] * fe_0 + ta1_z_yyzz_xyzz_0[i] * pa_y[i] - ta1_z_yyzz_xyzz_1[i] * pc_y[i];

        ta1_z_yyyzz_xzzz_0[i] =
            2.0 * ta1_z_yzz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xzzz_1[i] * fe_0 + ta1_z_yyzz_xzzz_0[i] * pa_y[i] - ta1_z_yyzz_xzzz_1[i] * pc_y[i];

        ta1_z_yyyzz_yyyy_0[i] = ta1_z_yyy_yyyy_0[i] * fe_0 - ta1_z_yyy_yyyy_1[i] * fe_0 + ta_yyyz_yyyy_1[i] + ta1_z_yyyz_yyyy_0[i] * pa_z[i] -
                                ta1_z_yyyz_yyyy_1[i] * pc_z[i];

        ta1_z_yyyzz_yyyz_0[i] = 2.0 * ta1_z_yzz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yyzz_yyz_0[i] * fe_0 -
                                3.0 * ta1_z_yyzz_yyz_1[i] * fe_0 + ta1_z_yyzz_yyyz_0[i] * pa_y[i] - ta1_z_yyzz_yyyz_1[i] * pc_y[i];

        ta1_z_yyyzz_yyzz_0[i] = 2.0 * ta1_z_yzz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yyzz_yzz_0[i] * fe_0 -
                                2.0 * ta1_z_yyzz_yzz_1[i] * fe_0 + ta1_z_yyzz_yyzz_0[i] * pa_y[i] - ta1_z_yyzz_yyzz_1[i] * pc_y[i];

        ta1_z_yyyzz_yzzz_0[i] = 2.0 * ta1_z_yzz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yzzz_1[i] * fe_0 + ta1_z_yyzz_zzz_0[i] * fe_0 -
                                ta1_z_yyzz_zzz_1[i] * fe_0 + ta1_z_yyzz_yzzz_0[i] * pa_y[i] - ta1_z_yyzz_yzzz_1[i] * pc_y[i];

        ta1_z_yyyzz_zzzz_0[i] =
            2.0 * ta1_z_yzz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_zzzz_1[i] * fe_0 + ta1_z_yyzz_zzzz_0[i] * pa_y[i] - ta1_z_yyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 900-915 components of targeted buffer : HG

    auto ta1_z_yyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 900);

    auto ta1_z_yyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 901);

    auto ta1_z_yyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 902);

    auto ta1_z_yyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 903);

    auto ta1_z_yyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 904);

    auto ta1_z_yyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 905);

    auto ta1_z_yyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 906);

    auto ta1_z_yyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 907);

    auto ta1_z_yyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 908);

    auto ta1_z_yyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 909);

    auto ta1_z_yyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 910);

    auto ta1_z_yyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 911);

    auto ta1_z_yyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 912);

    auto ta1_z_yyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 913);

    auto ta1_z_yyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 914);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyz_xxxy_0,   \
                             ta1_z_yyz_xxxy_1,   \
                             ta1_z_yyz_xxyy_0,   \
                             ta1_z_yyz_xxyy_1,   \
                             ta1_z_yyz_xyyy_0,   \
                             ta1_z_yyz_xyyy_1,   \
                             ta1_z_yyz_yyyy_0,   \
                             ta1_z_yyz_yyyy_1,   \
                             ta1_z_yyzz_xxxy_0,  \
                             ta1_z_yyzz_xxxy_1,  \
                             ta1_z_yyzz_xxyy_0,  \
                             ta1_z_yyzz_xxyy_1,  \
                             ta1_z_yyzz_xyyy_0,  \
                             ta1_z_yyzz_xyyy_1,  \
                             ta1_z_yyzz_yyyy_0,  \
                             ta1_z_yyzz_yyyy_1,  \
                             ta1_z_yyzzz_xxxx_0, \
                             ta1_z_yyzzz_xxxy_0, \
                             ta1_z_yyzzz_xxxz_0, \
                             ta1_z_yyzzz_xxyy_0, \
                             ta1_z_yyzzz_xxyz_0, \
                             ta1_z_yyzzz_xxzz_0, \
                             ta1_z_yyzzz_xyyy_0, \
                             ta1_z_yyzzz_xyyz_0, \
                             ta1_z_yyzzz_xyzz_0, \
                             ta1_z_yyzzz_xzzz_0, \
                             ta1_z_yyzzz_yyyy_0, \
                             ta1_z_yyzzz_yyyz_0, \
                             ta1_z_yyzzz_yyzz_0, \
                             ta1_z_yyzzz_yzzz_0, \
                             ta1_z_yyzzz_zzzz_0, \
                             ta1_z_yzzz_xxxx_0,  \
                             ta1_z_yzzz_xxxx_1,  \
                             ta1_z_yzzz_xxxz_0,  \
                             ta1_z_yzzz_xxxz_1,  \
                             ta1_z_yzzz_xxyz_0,  \
                             ta1_z_yzzz_xxyz_1,  \
                             ta1_z_yzzz_xxz_0,   \
                             ta1_z_yzzz_xxz_1,   \
                             ta1_z_yzzz_xxzz_0,  \
                             ta1_z_yzzz_xxzz_1,  \
                             ta1_z_yzzz_xyyz_0,  \
                             ta1_z_yzzz_xyyz_1,  \
                             ta1_z_yzzz_xyz_0,   \
                             ta1_z_yzzz_xyz_1,   \
                             ta1_z_yzzz_xyzz_0,  \
                             ta1_z_yzzz_xyzz_1,  \
                             ta1_z_yzzz_xzz_0,   \
                             ta1_z_yzzz_xzz_1,   \
                             ta1_z_yzzz_xzzz_0,  \
                             ta1_z_yzzz_xzzz_1,  \
                             ta1_z_yzzz_yyyz_0,  \
                             ta1_z_yzzz_yyyz_1,  \
                             ta1_z_yzzz_yyz_0,   \
                             ta1_z_yzzz_yyz_1,   \
                             ta1_z_yzzz_yyzz_0,  \
                             ta1_z_yzzz_yyzz_1,  \
                             ta1_z_yzzz_yzz_0,   \
                             ta1_z_yzzz_yzz_1,   \
                             ta1_z_yzzz_yzzz_0,  \
                             ta1_z_yzzz_yzzz_1,  \
                             ta1_z_yzzz_zzz_0,   \
                             ta1_z_yzzz_zzz_1,   \
                             ta1_z_yzzz_zzzz_0,  \
                             ta1_z_yzzz_zzzz_1,  \
                             ta1_z_zzz_xxxx_0,   \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxz_0,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxyz_0,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxzz_0,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xyyz_0,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyzz_0,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xzzz_0,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_yyyz_0,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyzz_0,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yzzz_0,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_zzzz_0,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta_yyzz_xxxy_1,     \
                             ta_yyzz_xxyy_1,     \
                             ta_yyzz_xyyy_1,     \
                             ta_yyzz_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzz_xxxx_0[i] =
            ta1_z_zzz_xxxx_0[i] * fe_0 - ta1_z_zzz_xxxx_1[i] * fe_0 + ta1_z_yzzz_xxxx_0[i] * pa_y[i] - ta1_z_yzzz_xxxx_1[i] * pc_y[i];

        ta1_z_yyzzz_xxxy_0[i] = 2.0 * ta1_z_yyz_xxxy_0[i] * fe_0 - 2.0 * ta1_z_yyz_xxxy_1[i] * fe_0 + ta_yyzz_xxxy_1[i] +
                                ta1_z_yyzz_xxxy_0[i] * pa_z[i] - ta1_z_yyzz_xxxy_1[i] * pc_z[i];

        ta1_z_yyzzz_xxxz_0[i] =
            ta1_z_zzz_xxxz_0[i] * fe_0 - ta1_z_zzz_xxxz_1[i] * fe_0 + ta1_z_yzzz_xxxz_0[i] * pa_y[i] - ta1_z_yzzz_xxxz_1[i] * pc_y[i];

        ta1_z_yyzzz_xxyy_0[i] = 2.0 * ta1_z_yyz_xxyy_0[i] * fe_0 - 2.0 * ta1_z_yyz_xxyy_1[i] * fe_0 + ta_yyzz_xxyy_1[i] +
                                ta1_z_yyzz_xxyy_0[i] * pa_z[i] - ta1_z_yyzz_xxyy_1[i] * pc_z[i];

        ta1_z_yyzzz_xxyz_0[i] = ta1_z_zzz_xxyz_0[i] * fe_0 - ta1_z_zzz_xxyz_1[i] * fe_0 + ta1_z_yzzz_xxz_0[i] * fe_0 - ta1_z_yzzz_xxz_1[i] * fe_0 +
                                ta1_z_yzzz_xxyz_0[i] * pa_y[i] - ta1_z_yzzz_xxyz_1[i] * pc_y[i];

        ta1_z_yyzzz_xxzz_0[i] =
            ta1_z_zzz_xxzz_0[i] * fe_0 - ta1_z_zzz_xxzz_1[i] * fe_0 + ta1_z_yzzz_xxzz_0[i] * pa_y[i] - ta1_z_yzzz_xxzz_1[i] * pc_y[i];

        ta1_z_yyzzz_xyyy_0[i] = 2.0 * ta1_z_yyz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyyy_1[i] * fe_0 + ta_yyzz_xyyy_1[i] +
                                ta1_z_yyzz_xyyy_0[i] * pa_z[i] - ta1_z_yyzz_xyyy_1[i] * pc_z[i];

        ta1_z_yyzzz_xyyz_0[i] = ta1_z_zzz_xyyz_0[i] * fe_0 - ta1_z_zzz_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_yzzz_xyz_1[i] * fe_0 + ta1_z_yzzz_xyyz_0[i] * pa_y[i] - ta1_z_yzzz_xyyz_1[i] * pc_y[i];

        ta1_z_yyzzz_xyzz_0[i] = ta1_z_zzz_xyzz_0[i] * fe_0 - ta1_z_zzz_xyzz_1[i] * fe_0 + ta1_z_yzzz_xzz_0[i] * fe_0 - ta1_z_yzzz_xzz_1[i] * fe_0 +
                                ta1_z_yzzz_xyzz_0[i] * pa_y[i] - ta1_z_yzzz_xyzz_1[i] * pc_y[i];

        ta1_z_yyzzz_xzzz_0[i] =
            ta1_z_zzz_xzzz_0[i] * fe_0 - ta1_z_zzz_xzzz_1[i] * fe_0 + ta1_z_yzzz_xzzz_0[i] * pa_y[i] - ta1_z_yzzz_xzzz_1[i] * pc_y[i];

        ta1_z_yyzzz_yyyy_0[i] = 2.0 * ta1_z_yyz_yyyy_0[i] * fe_0 - 2.0 * ta1_z_yyz_yyyy_1[i] * fe_0 + ta_yyzz_yyyy_1[i] +
                                ta1_z_yyzz_yyyy_0[i] * pa_z[i] - ta1_z_yyzz_yyyy_1[i] * pc_z[i];

        ta1_z_yyzzz_yyyz_0[i] = ta1_z_zzz_yyyz_0[i] * fe_0 - ta1_z_zzz_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yzzz_yyz_0[i] * fe_0 -
                                3.0 * ta1_z_yzzz_yyz_1[i] * fe_0 + ta1_z_yzzz_yyyz_0[i] * pa_y[i] - ta1_z_yzzz_yyyz_1[i] * pc_y[i];

        ta1_z_yyzzz_yyzz_0[i] = ta1_z_zzz_yyzz_0[i] * fe_0 - ta1_z_zzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yzzz_yzz_0[i] * fe_0 -
                                2.0 * ta1_z_yzzz_yzz_1[i] * fe_0 + ta1_z_yzzz_yyzz_0[i] * pa_y[i] - ta1_z_yzzz_yyzz_1[i] * pc_y[i];

        ta1_z_yyzzz_yzzz_0[i] = ta1_z_zzz_yzzz_0[i] * fe_0 - ta1_z_zzz_yzzz_1[i] * fe_0 + ta1_z_yzzz_zzz_0[i] * fe_0 - ta1_z_yzzz_zzz_1[i] * fe_0 +
                                ta1_z_yzzz_yzzz_0[i] * pa_y[i] - ta1_z_yzzz_yzzz_1[i] * pc_y[i];

        ta1_z_yyzzz_zzzz_0[i] =
            ta1_z_zzz_zzzz_0[i] * fe_0 - ta1_z_zzz_zzzz_1[i] * fe_0 + ta1_z_yzzz_zzzz_0[i] * pa_y[i] - ta1_z_yzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 915-930 components of targeted buffer : HG

    auto ta1_z_yzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 915);

    auto ta1_z_yzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 916);

    auto ta1_z_yzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 917);

    auto ta1_z_yzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 918);

    auto ta1_z_yzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 919);

    auto ta1_z_yzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 920);

    auto ta1_z_yzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 921);

    auto ta1_z_yzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 922);

    auto ta1_z_yzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 923);

    auto ta1_z_yzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 924);

    auto ta1_z_yzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 925);

    auto ta1_z_yzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 926);

    auto ta1_z_yzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 927);

    auto ta1_z_yzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 928);

    auto ta1_z_yzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 929);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yzzzz_xxxx_0, \
                             ta1_z_yzzzz_xxxy_0, \
                             ta1_z_yzzzz_xxxz_0, \
                             ta1_z_yzzzz_xxyy_0, \
                             ta1_z_yzzzz_xxyz_0, \
                             ta1_z_yzzzz_xxzz_0, \
                             ta1_z_yzzzz_xyyy_0, \
                             ta1_z_yzzzz_xyyz_0, \
                             ta1_z_yzzzz_xyzz_0, \
                             ta1_z_yzzzz_xzzz_0, \
                             ta1_z_yzzzz_yyyy_0, \
                             ta1_z_yzzzz_yyyz_0, \
                             ta1_z_yzzzz_yyzz_0, \
                             ta1_z_yzzzz_yzzz_0, \
                             ta1_z_yzzzz_zzzz_0, \
                             ta1_z_zzzz_xxx_0,   \
                             ta1_z_zzzz_xxx_1,   \
                             ta1_z_zzzz_xxxx_0,  \
                             ta1_z_zzzz_xxxx_1,  \
                             ta1_z_zzzz_xxxy_0,  \
                             ta1_z_zzzz_xxxy_1,  \
                             ta1_z_zzzz_xxxz_0,  \
                             ta1_z_zzzz_xxxz_1,  \
                             ta1_z_zzzz_xxy_0,   \
                             ta1_z_zzzz_xxy_1,   \
                             ta1_z_zzzz_xxyy_0,  \
                             ta1_z_zzzz_xxyy_1,  \
                             ta1_z_zzzz_xxyz_0,  \
                             ta1_z_zzzz_xxyz_1,  \
                             ta1_z_zzzz_xxz_0,   \
                             ta1_z_zzzz_xxz_1,   \
                             ta1_z_zzzz_xxzz_0,  \
                             ta1_z_zzzz_xxzz_1,  \
                             ta1_z_zzzz_xyy_0,   \
                             ta1_z_zzzz_xyy_1,   \
                             ta1_z_zzzz_xyyy_0,  \
                             ta1_z_zzzz_xyyy_1,  \
                             ta1_z_zzzz_xyyz_0,  \
                             ta1_z_zzzz_xyyz_1,  \
                             ta1_z_zzzz_xyz_0,   \
                             ta1_z_zzzz_xyz_1,   \
                             ta1_z_zzzz_xyzz_0,  \
                             ta1_z_zzzz_xyzz_1,  \
                             ta1_z_zzzz_xzz_0,   \
                             ta1_z_zzzz_xzz_1,   \
                             ta1_z_zzzz_xzzz_0,  \
                             ta1_z_zzzz_xzzz_1,  \
                             ta1_z_zzzz_yyy_0,   \
                             ta1_z_zzzz_yyy_1,   \
                             ta1_z_zzzz_yyyy_0,  \
                             ta1_z_zzzz_yyyy_1,  \
                             ta1_z_zzzz_yyyz_0,  \
                             ta1_z_zzzz_yyyz_1,  \
                             ta1_z_zzzz_yyz_0,   \
                             ta1_z_zzzz_yyz_1,   \
                             ta1_z_zzzz_yyzz_0,  \
                             ta1_z_zzzz_yyzz_1,  \
                             ta1_z_zzzz_yzz_0,   \
                             ta1_z_zzzz_yzz_1,   \
                             ta1_z_zzzz_yzzz_0,  \
                             ta1_z_zzzz_yzzz_1,  \
                             ta1_z_zzzz_zzz_0,   \
                             ta1_z_zzzz_zzz_1,   \
                             ta1_z_zzzz_zzzz_0,  \
                             ta1_z_zzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzz_xxxx_0[i] = ta1_z_zzzz_xxxx_0[i] * pa_y[i] - ta1_z_zzzz_xxxx_1[i] * pc_y[i];

        ta1_z_yzzzz_xxxy_0[i] =
            ta1_z_zzzz_xxx_0[i] * fe_0 - ta1_z_zzzz_xxx_1[i] * fe_0 + ta1_z_zzzz_xxxy_0[i] * pa_y[i] - ta1_z_zzzz_xxxy_1[i] * pc_y[i];

        ta1_z_yzzzz_xxxz_0[i] = ta1_z_zzzz_xxxz_0[i] * pa_y[i] - ta1_z_zzzz_xxxz_1[i] * pc_y[i];

        ta1_z_yzzzz_xxyy_0[i] =
            2.0 * ta1_z_zzzz_xxy_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xxy_1[i] * fe_0 + ta1_z_zzzz_xxyy_0[i] * pa_y[i] - ta1_z_zzzz_xxyy_1[i] * pc_y[i];

        ta1_z_yzzzz_xxyz_0[i] =
            ta1_z_zzzz_xxz_0[i] * fe_0 - ta1_z_zzzz_xxz_1[i] * fe_0 + ta1_z_zzzz_xxyz_0[i] * pa_y[i] - ta1_z_zzzz_xxyz_1[i] * pc_y[i];

        ta1_z_yzzzz_xxzz_0[i] = ta1_z_zzzz_xxzz_0[i] * pa_y[i] - ta1_z_zzzz_xxzz_1[i] * pc_y[i];

        ta1_z_yzzzz_xyyy_0[i] =
            3.0 * ta1_z_zzzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_zzzz_xyy_1[i] * fe_0 + ta1_z_zzzz_xyyy_0[i] * pa_y[i] - ta1_z_zzzz_xyyy_1[i] * pc_y[i];

        ta1_z_yzzzz_xyyz_0[i] =
            2.0 * ta1_z_zzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xyz_1[i] * fe_0 + ta1_z_zzzz_xyyz_0[i] * pa_y[i] - ta1_z_zzzz_xyyz_1[i] * pc_y[i];

        ta1_z_yzzzz_xyzz_0[i] =
            ta1_z_zzzz_xzz_0[i] * fe_0 - ta1_z_zzzz_xzz_1[i] * fe_0 + ta1_z_zzzz_xyzz_0[i] * pa_y[i] - ta1_z_zzzz_xyzz_1[i] * pc_y[i];

        ta1_z_yzzzz_xzzz_0[i] = ta1_z_zzzz_xzzz_0[i] * pa_y[i] - ta1_z_zzzz_xzzz_1[i] * pc_y[i];

        ta1_z_yzzzz_yyyy_0[i] =
            4.0 * ta1_z_zzzz_yyy_0[i] * fe_0 - 4.0 * ta1_z_zzzz_yyy_1[i] * fe_0 + ta1_z_zzzz_yyyy_0[i] * pa_y[i] - ta1_z_zzzz_yyyy_1[i] * pc_y[i];

        ta1_z_yzzzz_yyyz_0[i] =
            3.0 * ta1_z_zzzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_zzzz_yyz_1[i] * fe_0 + ta1_z_zzzz_yyyz_0[i] * pa_y[i] - ta1_z_zzzz_yyyz_1[i] * pc_y[i];

        ta1_z_yzzzz_yyzz_0[i] =
            2.0 * ta1_z_zzzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_yzz_1[i] * fe_0 + ta1_z_zzzz_yyzz_0[i] * pa_y[i] - ta1_z_zzzz_yyzz_1[i] * pc_y[i];

        ta1_z_yzzzz_yzzz_0[i] =
            ta1_z_zzzz_zzz_0[i] * fe_0 - ta1_z_zzzz_zzz_1[i] * fe_0 + ta1_z_zzzz_yzzz_0[i] * pa_y[i] - ta1_z_zzzz_yzzz_1[i] * pc_y[i];

        ta1_z_yzzzz_zzzz_0[i] = ta1_z_zzzz_zzzz_0[i] * pa_y[i] - ta1_z_zzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 930-945 components of targeted buffer : HG

    auto ta1_z_zzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 930);

    auto ta1_z_zzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 931);

    auto ta1_z_zzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 932);

    auto ta1_z_zzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 933);

    auto ta1_z_zzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 934);

    auto ta1_z_zzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 935);

    auto ta1_z_zzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 936);

    auto ta1_z_zzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 937);

    auto ta1_z_zzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 938);

    auto ta1_z_zzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 939);

    auto ta1_z_zzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 940);

    auto ta1_z_zzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 941);

    auto ta1_z_zzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 942);

    auto ta1_z_zzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 943);

    auto ta1_z_zzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 944);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_z_zzz_xxxx_0,   \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxy_0,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxz_0,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxyy_0,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyz_0,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxzz_0,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xyyy_0,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyz_0,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyzz_0,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xzzz_0,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_yyyy_0,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyz_0,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyzz_0,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yzzz_0,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_zzzz_0,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta1_z_zzzz_xxx_0,   \
                             ta1_z_zzzz_xxx_1,   \
                             ta1_z_zzzz_xxxx_0,  \
                             ta1_z_zzzz_xxxx_1,  \
                             ta1_z_zzzz_xxxy_0,  \
                             ta1_z_zzzz_xxxy_1,  \
                             ta1_z_zzzz_xxxz_0,  \
                             ta1_z_zzzz_xxxz_1,  \
                             ta1_z_zzzz_xxy_0,   \
                             ta1_z_zzzz_xxy_1,   \
                             ta1_z_zzzz_xxyy_0,  \
                             ta1_z_zzzz_xxyy_1,  \
                             ta1_z_zzzz_xxyz_0,  \
                             ta1_z_zzzz_xxyz_1,  \
                             ta1_z_zzzz_xxz_0,   \
                             ta1_z_zzzz_xxz_1,   \
                             ta1_z_zzzz_xxzz_0,  \
                             ta1_z_zzzz_xxzz_1,  \
                             ta1_z_zzzz_xyy_0,   \
                             ta1_z_zzzz_xyy_1,   \
                             ta1_z_zzzz_xyyy_0,  \
                             ta1_z_zzzz_xyyy_1,  \
                             ta1_z_zzzz_xyyz_0,  \
                             ta1_z_zzzz_xyyz_1,  \
                             ta1_z_zzzz_xyz_0,   \
                             ta1_z_zzzz_xyz_1,   \
                             ta1_z_zzzz_xyzz_0,  \
                             ta1_z_zzzz_xyzz_1,  \
                             ta1_z_zzzz_xzz_0,   \
                             ta1_z_zzzz_xzz_1,   \
                             ta1_z_zzzz_xzzz_0,  \
                             ta1_z_zzzz_xzzz_1,  \
                             ta1_z_zzzz_yyy_0,   \
                             ta1_z_zzzz_yyy_1,   \
                             ta1_z_zzzz_yyyy_0,  \
                             ta1_z_zzzz_yyyy_1,  \
                             ta1_z_zzzz_yyyz_0,  \
                             ta1_z_zzzz_yyyz_1,  \
                             ta1_z_zzzz_yyz_0,   \
                             ta1_z_zzzz_yyz_1,   \
                             ta1_z_zzzz_yyzz_0,  \
                             ta1_z_zzzz_yyzz_1,  \
                             ta1_z_zzzz_yzz_0,   \
                             ta1_z_zzzz_yzz_1,   \
                             ta1_z_zzzz_yzzz_0,  \
                             ta1_z_zzzz_yzzz_1,  \
                             ta1_z_zzzz_zzz_0,   \
                             ta1_z_zzzz_zzz_1,   \
                             ta1_z_zzzz_zzzz_0,  \
                             ta1_z_zzzz_zzzz_1,  \
                             ta1_z_zzzzz_xxxx_0, \
                             ta1_z_zzzzz_xxxy_0, \
                             ta1_z_zzzzz_xxxz_0, \
                             ta1_z_zzzzz_xxyy_0, \
                             ta1_z_zzzzz_xxyz_0, \
                             ta1_z_zzzzz_xxzz_0, \
                             ta1_z_zzzzz_xyyy_0, \
                             ta1_z_zzzzz_xyyz_0, \
                             ta1_z_zzzzz_xyzz_0, \
                             ta1_z_zzzzz_xzzz_0, \
                             ta1_z_zzzzz_yyyy_0, \
                             ta1_z_zzzzz_yyyz_0, \
                             ta1_z_zzzzz_yyzz_0, \
                             ta1_z_zzzzz_yzzz_0, \
                             ta1_z_zzzzz_zzzz_0, \
                             ta_zzzz_xxxx_1,     \
                             ta_zzzz_xxxy_1,     \
                             ta_zzzz_xxxz_1,     \
                             ta_zzzz_xxyy_1,     \
                             ta_zzzz_xxyz_1,     \
                             ta_zzzz_xxzz_1,     \
                             ta_zzzz_xyyy_1,     \
                             ta_zzzz_xyyz_1,     \
                             ta_zzzz_xyzz_1,     \
                             ta_zzzz_xzzz_1,     \
                             ta_zzzz_yyyy_1,     \
                             ta_zzzz_yyyz_1,     \
                             ta_zzzz_yyzz_1,     \
                             ta_zzzz_yzzz_1,     \
                             ta_zzzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzz_xxxx_0[i] = 4.0 * ta1_z_zzz_xxxx_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxx_1[i] * fe_0 + ta_zzzz_xxxx_1[i] +
                                ta1_z_zzzz_xxxx_0[i] * pa_z[i] - ta1_z_zzzz_xxxx_1[i] * pc_z[i];

        ta1_z_zzzzz_xxxy_0[i] = 4.0 * ta1_z_zzz_xxxy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxy_1[i] * fe_0 + ta_zzzz_xxxy_1[i] +
                                ta1_z_zzzz_xxxy_0[i] * pa_z[i] - ta1_z_zzzz_xxxy_1[i] * pc_z[i];

        ta1_z_zzzzz_xxxz_0[i] = 4.0 * ta1_z_zzz_xxxz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxz_1[i] * fe_0 + ta1_z_zzzz_xxx_0[i] * fe_0 -
                                ta1_z_zzzz_xxx_1[i] * fe_0 + ta_zzzz_xxxz_1[i] + ta1_z_zzzz_xxxz_0[i] * pa_z[i] - ta1_z_zzzz_xxxz_1[i] * pc_z[i];

        ta1_z_zzzzz_xxyy_0[i] = 4.0 * ta1_z_zzz_xxyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxyy_1[i] * fe_0 + ta_zzzz_xxyy_1[i] +
                                ta1_z_zzzz_xxyy_0[i] * pa_z[i] - ta1_z_zzzz_xxyy_1[i] * pc_z[i];

        ta1_z_zzzzz_xxyz_0[i] = 4.0 * ta1_z_zzz_xxyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxyz_1[i] * fe_0 + ta1_z_zzzz_xxy_0[i] * fe_0 -
                                ta1_z_zzzz_xxy_1[i] * fe_0 + ta_zzzz_xxyz_1[i] + ta1_z_zzzz_xxyz_0[i] * pa_z[i] - ta1_z_zzzz_xxyz_1[i] * pc_z[i];

        ta1_z_zzzzz_xxzz_0[i] = 4.0 * ta1_z_zzz_xxzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_zzzz_xxz_0[i] * fe_0 -
                                2.0 * ta1_z_zzzz_xxz_1[i] * fe_0 + ta_zzzz_xxzz_1[i] + ta1_z_zzzz_xxzz_0[i] * pa_z[i] -
                                ta1_z_zzzz_xxzz_1[i] * pc_z[i];

        ta1_z_zzzzz_xyyy_0[i] = 4.0 * ta1_z_zzz_xyyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyyy_1[i] * fe_0 + ta_zzzz_xyyy_1[i] +
                                ta1_z_zzzz_xyyy_0[i] * pa_z[i] - ta1_z_zzzz_xyyy_1[i] * pc_z[i];

        ta1_z_zzzzz_xyyz_0[i] = 4.0 * ta1_z_zzz_xyyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyyz_1[i] * fe_0 + ta1_z_zzzz_xyy_0[i] * fe_0 -
                                ta1_z_zzzz_xyy_1[i] * fe_0 + ta_zzzz_xyyz_1[i] + ta1_z_zzzz_xyyz_0[i] * pa_z[i] - ta1_z_zzzz_xyyz_1[i] * pc_z[i];

        ta1_z_zzzzz_xyzz_0[i] = 4.0 * ta1_z_zzz_xyzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyzz_1[i] * fe_0 + 2.0 * ta1_z_zzzz_xyz_0[i] * fe_0 -
                                2.0 * ta1_z_zzzz_xyz_1[i] * fe_0 + ta_zzzz_xyzz_1[i] + ta1_z_zzzz_xyzz_0[i] * pa_z[i] -
                                ta1_z_zzzz_xyzz_1[i] * pc_z[i];

        ta1_z_zzzzz_xzzz_0[i] = 4.0 * ta1_z_zzz_xzzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xzzz_1[i] * fe_0 + 3.0 * ta1_z_zzzz_xzz_0[i] * fe_0 -
                                3.0 * ta1_z_zzzz_xzz_1[i] * fe_0 + ta_zzzz_xzzz_1[i] + ta1_z_zzzz_xzzz_0[i] * pa_z[i] -
                                ta1_z_zzzz_xzzz_1[i] * pc_z[i];

        ta1_z_zzzzz_yyyy_0[i] = 4.0 * ta1_z_zzz_yyyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyyy_1[i] * fe_0 + ta_zzzz_yyyy_1[i] +
                                ta1_z_zzzz_yyyy_0[i] * pa_z[i] - ta1_z_zzzz_yyyy_1[i] * pc_z[i];

        ta1_z_zzzzz_yyyz_0[i] = 4.0 * ta1_z_zzz_yyyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyyz_1[i] * fe_0 + ta1_z_zzzz_yyy_0[i] * fe_0 -
                                ta1_z_zzzz_yyy_1[i] * fe_0 + ta_zzzz_yyyz_1[i] + ta1_z_zzzz_yyyz_0[i] * pa_z[i] - ta1_z_zzzz_yyyz_1[i] * pc_z[i];

        ta1_z_zzzzz_yyzz_0[i] = 4.0 * ta1_z_zzz_yyzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_zzzz_yyz_0[i] * fe_0 -
                                2.0 * ta1_z_zzzz_yyz_1[i] * fe_0 + ta_zzzz_yyzz_1[i] + ta1_z_zzzz_yyzz_0[i] * pa_z[i] -
                                ta1_z_zzzz_yyzz_1[i] * pc_z[i];

        ta1_z_zzzzz_yzzz_0[i] = 4.0 * ta1_z_zzz_yzzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yzzz_1[i] * fe_0 + 3.0 * ta1_z_zzzz_yzz_0[i] * fe_0 -
                                3.0 * ta1_z_zzzz_yzz_1[i] * fe_0 + ta_zzzz_yzzz_1[i] + ta1_z_zzzz_yzzz_0[i] * pa_z[i] -
                                ta1_z_zzzz_yzzz_1[i] * pc_z[i];

        ta1_z_zzzzz_zzzz_0[i] = 4.0 * ta1_z_zzz_zzzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_zzzz_1[i] * fe_0 + 4.0 * ta1_z_zzzz_zzz_0[i] * fe_0 -
                                4.0 * ta1_z_zzzz_zzz_1[i] * fe_0 + ta_zzzz_zzzz_1[i] + ta1_z_zzzz_zzzz_0[i] * pa_z[i] -
                                ta1_z_zzzz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
