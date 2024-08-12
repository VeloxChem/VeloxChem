#include "NuclearPotentialGeom010PrimRecGG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gg,
                                        const size_t              idx_npot_geom_010_0_dg,
                                        const size_t              idx_npot_geom_010_1_dg,
                                        const size_t              idx_npot_geom_010_0_ff,
                                        const size_t              idx_npot_geom_010_1_ff,
                                        const size_t              idx_npot_1_fg,
                                        const size_t              idx_npot_geom_010_0_fg,
                                        const size_t              idx_npot_geom_010_1_fg,
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

    // Set up components of auxiliary buffer : DG

    auto ta1_x_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg);

    auto ta1_x_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 1);

    auto ta1_x_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 2);

    auto ta1_x_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 3);

    auto ta1_x_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 4);

    auto ta1_x_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 5);

    auto ta1_x_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 6);

    auto ta1_x_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 7);

    auto ta1_x_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 8);

    auto ta1_x_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 9);

    auto ta1_x_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 10);

    auto ta1_x_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 11);

    auto ta1_x_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 12);

    auto ta1_x_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 13);

    auto ta1_x_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 14);

    auto ta1_x_xy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 15);

    auto ta1_x_xy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 17);

    auto ta1_x_xy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 20);

    auto ta1_x_xy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 24);

    auto ta1_x_xz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 30);

    auto ta1_x_xz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 31);

    auto ta1_x_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 32);

    auto ta1_x_xz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 33);

    auto ta1_x_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 35);

    auto ta1_x_xz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 36);

    auto ta1_x_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 39);

    auto ta1_x_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 45);

    auto ta1_x_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 46);

    auto ta1_x_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 47);

    auto ta1_x_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 48);

    auto ta1_x_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 49);

    auto ta1_x_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 50);

    auto ta1_x_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 51);

    auto ta1_x_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 52);

    auto ta1_x_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 53);

    auto ta1_x_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 54);

    auto ta1_x_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 55);

    auto ta1_x_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 56);

    auto ta1_x_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 57);

    auto ta1_x_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 58);

    auto ta1_x_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 59);

    auto ta1_x_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 62);

    auto ta1_x_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 65);

    auto ta1_x_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 69);

    auto ta1_x_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 74);

    auto ta1_x_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 75);

    auto ta1_x_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 76);

    auto ta1_x_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 77);

    auto ta1_x_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 78);

    auto ta1_x_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 79);

    auto ta1_x_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 80);

    auto ta1_x_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 81);

    auto ta1_x_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 82);

    auto ta1_x_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 83);

    auto ta1_x_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 84);

    auto ta1_x_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 85);

    auto ta1_x_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 86);

    auto ta1_x_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 87);

    auto ta1_x_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 88);

    auto ta1_x_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 89);

    auto ta1_y_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 90);

    auto ta1_y_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 91);

    auto ta1_y_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 92);

    auto ta1_y_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 93);

    auto ta1_y_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 94);

    auto ta1_y_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 95);

    auto ta1_y_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 96);

    auto ta1_y_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 97);

    auto ta1_y_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 98);

    auto ta1_y_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 99);

    auto ta1_y_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 100);

    auto ta1_y_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 101);

    auto ta1_y_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 102);

    auto ta1_y_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 103);

    auto ta1_y_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 104);

    auto ta1_y_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 115);

    auto ta1_y_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 116);

    auto ta1_y_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 117);

    auto ta1_y_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 118);

    auto ta1_y_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 131);

    auto ta1_y_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 132);

    auto ta1_y_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 133);

    auto ta1_y_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 134);

    auto ta1_y_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 135);

    auto ta1_y_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 136);

    auto ta1_y_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 137);

    auto ta1_y_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 138);

    auto ta1_y_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 139);

    auto ta1_y_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 140);

    auto ta1_y_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 141);

    auto ta1_y_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 142);

    auto ta1_y_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 143);

    auto ta1_y_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 144);

    auto ta1_y_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 145);

    auto ta1_y_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 146);

    auto ta1_y_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 147);

    auto ta1_y_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 148);

    auto ta1_y_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 149);

    auto ta1_y_yz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 151);

    auto ta1_y_yz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 153);

    auto ta1_y_yz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 156);

    auto ta1_y_yz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 160);

    auto ta1_y_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 161);

    auto ta1_y_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 162);

    auto ta1_y_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 163);

    auto ta1_y_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 165);

    auto ta1_y_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 166);

    auto ta1_y_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 167);

    auto ta1_y_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 168);

    auto ta1_y_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 169);

    auto ta1_y_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 170);

    auto ta1_y_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 171);

    auto ta1_y_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 172);

    auto ta1_y_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 173);

    auto ta1_y_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 174);

    auto ta1_y_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 175);

    auto ta1_y_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 176);

    auto ta1_y_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 177);

    auto ta1_y_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 178);

    auto ta1_y_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 179);

    auto ta1_z_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 180);

    auto ta1_z_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 181);

    auto ta1_z_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 182);

    auto ta1_z_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 183);

    auto ta1_z_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 184);

    auto ta1_z_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 185);

    auto ta1_z_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 186);

    auto ta1_z_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 187);

    auto ta1_z_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 188);

    auto ta1_z_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 189);

    auto ta1_z_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 190);

    auto ta1_z_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 191);

    auto ta1_z_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 192);

    auto ta1_z_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 193);

    auto ta1_z_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 194);

    auto ta1_z_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 205);

    auto ta1_z_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 206);

    auto ta1_z_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 207);

    auto ta1_z_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 208);

    auto ta1_z_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 221);

    auto ta1_z_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 222);

    auto ta1_z_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 223);

    auto ta1_z_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 224);

    auto ta1_z_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 225);

    auto ta1_z_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 226);

    auto ta1_z_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 227);

    auto ta1_z_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 228);

    auto ta1_z_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 229);

    auto ta1_z_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 230);

    auto ta1_z_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 231);

    auto ta1_z_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 232);

    auto ta1_z_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 233);

    auto ta1_z_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 234);

    auto ta1_z_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 235);

    auto ta1_z_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 236);

    auto ta1_z_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 237);

    auto ta1_z_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 238);

    auto ta1_z_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 239);

    auto ta1_z_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 242);

    auto ta1_z_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 245);

    auto ta1_z_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 249);

    auto ta1_z_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 251);

    auto ta1_z_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 252);

    auto ta1_z_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 253);

    auto ta1_z_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 254);

    auto ta1_z_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 255);

    auto ta1_z_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 256);

    auto ta1_z_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 257);

    auto ta1_z_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 258);

    auto ta1_z_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 259);

    auto ta1_z_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 260);

    auto ta1_z_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 261);

    auto ta1_z_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 262);

    auto ta1_z_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 263);

    auto ta1_z_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 264);

    auto ta1_z_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 265);

    auto ta1_z_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 266);

    auto ta1_z_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 267);

    auto ta1_z_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 268);

    auto ta1_z_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 269);

    // Set up components of auxiliary buffer : DG

    auto ta1_x_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg);

    auto ta1_x_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 1);

    auto ta1_x_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 2);

    auto ta1_x_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 3);

    auto ta1_x_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 4);

    auto ta1_x_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 5);

    auto ta1_x_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 6);

    auto ta1_x_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 7);

    auto ta1_x_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 8);

    auto ta1_x_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 9);

    auto ta1_x_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 10);

    auto ta1_x_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 11);

    auto ta1_x_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 12);

    auto ta1_x_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 13);

    auto ta1_x_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 14);

    auto ta1_x_xy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 15);

    auto ta1_x_xy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 17);

    auto ta1_x_xy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 20);

    auto ta1_x_xy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 24);

    auto ta1_x_xz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 30);

    auto ta1_x_xz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 31);

    auto ta1_x_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 32);

    auto ta1_x_xz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 33);

    auto ta1_x_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 35);

    auto ta1_x_xz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 36);

    auto ta1_x_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 39);

    auto ta1_x_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 45);

    auto ta1_x_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 46);

    auto ta1_x_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 47);

    auto ta1_x_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 48);

    auto ta1_x_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 49);

    auto ta1_x_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 50);

    auto ta1_x_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 51);

    auto ta1_x_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 52);

    auto ta1_x_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 53);

    auto ta1_x_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 54);

    auto ta1_x_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 55);

    auto ta1_x_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 56);

    auto ta1_x_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 57);

    auto ta1_x_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 58);

    auto ta1_x_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 59);

    auto ta1_x_yz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 62);

    auto ta1_x_yz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 65);

    auto ta1_x_yz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 69);

    auto ta1_x_yz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 74);

    auto ta1_x_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 75);

    auto ta1_x_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 76);

    auto ta1_x_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 77);

    auto ta1_x_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 78);

    auto ta1_x_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 79);

    auto ta1_x_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 80);

    auto ta1_x_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 81);

    auto ta1_x_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 82);

    auto ta1_x_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 83);

    auto ta1_x_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 84);

    auto ta1_x_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 85);

    auto ta1_x_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 86);

    auto ta1_x_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 87);

    auto ta1_x_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 88);

    auto ta1_x_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 89);

    auto ta1_y_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 90);

    auto ta1_y_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 91);

    auto ta1_y_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 92);

    auto ta1_y_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 93);

    auto ta1_y_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 94);

    auto ta1_y_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 95);

    auto ta1_y_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 96);

    auto ta1_y_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 97);

    auto ta1_y_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 98);

    auto ta1_y_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 99);

    auto ta1_y_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 100);

    auto ta1_y_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 101);

    auto ta1_y_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 102);

    auto ta1_y_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 103);

    auto ta1_y_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 104);

    auto ta1_y_xy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 115);

    auto ta1_y_xy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 116);

    auto ta1_y_xy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 117);

    auto ta1_y_xy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 118);

    auto ta1_y_xz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 131);

    auto ta1_y_xz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 132);

    auto ta1_y_xz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 133);

    auto ta1_y_xz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 134);

    auto ta1_y_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 135);

    auto ta1_y_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 136);

    auto ta1_y_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 137);

    auto ta1_y_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 138);

    auto ta1_y_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 139);

    auto ta1_y_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 140);

    auto ta1_y_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 141);

    auto ta1_y_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 142);

    auto ta1_y_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 143);

    auto ta1_y_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 144);

    auto ta1_y_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 145);

    auto ta1_y_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 146);

    auto ta1_y_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 147);

    auto ta1_y_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 148);

    auto ta1_y_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 149);

    auto ta1_y_yz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 151);

    auto ta1_y_yz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 153);

    auto ta1_y_yz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 156);

    auto ta1_y_yz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 160);

    auto ta1_y_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 161);

    auto ta1_y_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 162);

    auto ta1_y_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 163);

    auto ta1_y_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 165);

    auto ta1_y_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 166);

    auto ta1_y_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 167);

    auto ta1_y_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 168);

    auto ta1_y_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 169);

    auto ta1_y_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 170);

    auto ta1_y_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 171);

    auto ta1_y_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 172);

    auto ta1_y_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 173);

    auto ta1_y_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 174);

    auto ta1_y_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 175);

    auto ta1_y_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 176);

    auto ta1_y_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 177);

    auto ta1_y_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 178);

    auto ta1_y_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 179);

    auto ta1_z_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 180);

    auto ta1_z_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 181);

    auto ta1_z_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 182);

    auto ta1_z_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 183);

    auto ta1_z_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 184);

    auto ta1_z_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 185);

    auto ta1_z_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 186);

    auto ta1_z_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 187);

    auto ta1_z_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 188);

    auto ta1_z_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 189);

    auto ta1_z_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 190);

    auto ta1_z_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 191);

    auto ta1_z_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 192);

    auto ta1_z_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 193);

    auto ta1_z_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 194);

    auto ta1_z_xy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 205);

    auto ta1_z_xy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 206);

    auto ta1_z_xy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 207);

    auto ta1_z_xy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 208);

    auto ta1_z_xz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 221);

    auto ta1_z_xz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 222);

    auto ta1_z_xz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 223);

    auto ta1_z_xz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 224);

    auto ta1_z_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 225);

    auto ta1_z_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 226);

    auto ta1_z_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 227);

    auto ta1_z_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 228);

    auto ta1_z_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 229);

    auto ta1_z_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 230);

    auto ta1_z_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 231);

    auto ta1_z_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 232);

    auto ta1_z_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 233);

    auto ta1_z_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 234);

    auto ta1_z_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 235);

    auto ta1_z_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 236);

    auto ta1_z_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 237);

    auto ta1_z_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 238);

    auto ta1_z_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 239);

    auto ta1_z_yz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 242);

    auto ta1_z_yz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 245);

    auto ta1_z_yz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 249);

    auto ta1_z_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 251);

    auto ta1_z_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 252);

    auto ta1_z_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 253);

    auto ta1_z_yz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 254);

    auto ta1_z_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 255);

    auto ta1_z_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 256);

    auto ta1_z_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 257);

    auto ta1_z_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 258);

    auto ta1_z_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 259);

    auto ta1_z_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 260);

    auto ta1_z_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 261);

    auto ta1_z_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 262);

    auto ta1_z_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 263);

    auto ta1_z_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 264);

    auto ta1_z_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 265);

    auto ta1_z_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 266);

    auto ta1_z_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 267);

    auto ta1_z_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 268);

    auto ta1_z_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 269);

    // Set up components of auxiliary buffer : FF

    auto ta1_x_xxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff);

    auto ta1_x_xxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 1);

    auto ta1_x_xxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 2);

    auto ta1_x_xxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 3);

    auto ta1_x_xxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 4);

    auto ta1_x_xxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 5);

    auto ta1_x_xxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 6);

    auto ta1_x_xxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 7);

    auto ta1_x_xxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 8);

    auto ta1_x_xxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 9);

    auto ta1_x_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 10);

    auto ta1_x_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 11);

    auto ta1_x_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 12);

    auto ta1_x_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 13);

    auto ta1_x_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 14);

    auto ta1_x_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 15);

    auto ta1_x_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 20);

    auto ta1_x_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 21);

    auto ta1_x_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 22);

    auto ta1_x_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 23);

    auto ta1_x_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 24);

    auto ta1_x_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 25);

    auto ta1_x_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 27);

    auto ta1_x_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 28);

    auto ta1_x_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 29);

    auto ta1_x_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 31);

    auto ta1_x_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 33);

    auto ta1_x_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 34);

    auto ta1_x_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 50);

    auto ta1_x_xzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 51);

    auto ta1_x_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 52);

    auto ta1_x_xzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 53);

    auto ta1_x_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 54);

    auto ta1_x_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 55);

    auto ta1_x_yyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 60);

    auto ta1_x_yyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 61);

    auto ta1_x_yyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 62);

    auto ta1_x_yyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 63);

    auto ta1_x_yyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 64);

    auto ta1_x_yyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 65);

    auto ta1_x_yyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 66);

    auto ta1_x_yyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 67);

    auto ta1_x_yyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 68);

    auto ta1_x_yyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 69);

    auto ta1_x_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 82);

    auto ta1_x_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 84);

    auto ta1_x_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 85);

    auto ta1_x_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 87);

    auto ta1_x_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 88);

    auto ta1_x_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 89);

    auto ta1_x_zzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 90);

    auto ta1_x_zzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 91);

    auto ta1_x_zzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 92);

    auto ta1_x_zzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 93);

    auto ta1_x_zzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 94);

    auto ta1_x_zzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 95);

    auto ta1_x_zzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 96);

    auto ta1_x_zzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 97);

    auto ta1_x_zzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 98);

    auto ta1_x_zzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 99);

    auto ta1_y_xxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 100);

    auto ta1_y_xxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 101);

    auto ta1_y_xxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 102);

    auto ta1_y_xxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 103);

    auto ta1_y_xxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 104);

    auto ta1_y_xxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 105);

    auto ta1_y_xxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 106);

    auto ta1_y_xxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 107);

    auto ta1_y_xxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 108);

    auto ta1_y_xxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 109);

    auto ta1_y_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 111);

    auto ta1_y_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 113);

    auto ta1_y_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 114);

    auto ta1_y_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 131);

    auto ta1_y_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 133);

    auto ta1_y_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 134);

    auto ta1_y_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 136);

    auto ta1_y_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 137);

    auto ta1_y_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 138);

    auto ta1_y_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 152);

    auto ta1_y_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 154);

    auto ta1_y_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 155);

    auto ta1_y_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 157);

    auto ta1_y_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 158);

    auto ta1_y_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 159);

    auto ta1_y_yyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 160);

    auto ta1_y_yyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 161);

    auto ta1_y_yyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 162);

    auto ta1_y_yyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 163);

    auto ta1_y_yyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 164);

    auto ta1_y_yyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 165);

    auto ta1_y_yyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 166);

    auto ta1_y_yyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 167);

    auto ta1_y_yyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 168);

    auto ta1_y_yyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 169);

    auto ta1_y_yyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 171);

    auto ta1_y_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 172);

    auto ta1_y_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 173);

    auto ta1_y_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 174);

    auto ta1_y_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 175);

    auto ta1_y_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 176);

    auto ta1_y_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 177);

    auto ta1_y_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 178);

    auto ta1_y_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 179);

    auto ta1_y_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 181);

    auto ta1_y_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 183);

    auto ta1_y_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 184);

    auto ta1_y_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 186);

    auto ta1_y_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 187);

    auto ta1_y_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 188);

    auto ta1_y_zzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 190);

    auto ta1_y_zzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 191);

    auto ta1_y_zzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 192);

    auto ta1_y_zzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 193);

    auto ta1_y_zzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 194);

    auto ta1_y_zzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 195);

    auto ta1_y_zzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 196);

    auto ta1_y_zzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 197);

    auto ta1_y_zzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 198);

    auto ta1_y_zzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 199);

    auto ta1_z_xxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 200);

    auto ta1_z_xxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 201);

    auto ta1_z_xxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 202);

    auto ta1_z_xxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 203);

    auto ta1_z_xxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 204);

    auto ta1_z_xxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 205);

    auto ta1_z_xxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 206);

    auto ta1_z_xxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 207);

    auto ta1_z_xxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 208);

    auto ta1_z_xxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 209);

    auto ta1_z_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 222);

    auto ta1_z_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 224);

    auto ta1_z_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 225);

    auto ta1_z_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 231);

    auto ta1_z_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 233);

    auto ta1_z_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 234);

    auto ta1_z_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 236);

    auto ta1_z_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 237);

    auto ta1_z_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 238);

    auto ta1_z_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 252);

    auto ta1_z_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 254);

    auto ta1_z_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 255);

    auto ta1_z_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 257);

    auto ta1_z_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 258);

    auto ta1_z_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 259);

    auto ta1_z_yyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 260);

    auto ta1_z_yyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 261);

    auto ta1_z_yyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 262);

    auto ta1_z_yyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 263);

    auto ta1_z_yyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 264);

    auto ta1_z_yyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 265);

    auto ta1_z_yyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 266);

    auto ta1_z_yyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 267);

    auto ta1_z_yyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 268);

    auto ta1_z_yyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 269);

    auto ta1_z_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 272);

    auto ta1_z_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 274);

    auto ta1_z_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 275);

    auto ta1_z_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 277);

    auto ta1_z_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 278);

    auto ta1_z_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 279);

    auto ta1_z_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 281);

    auto ta1_z_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 282);

    auto ta1_z_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 283);

    auto ta1_z_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 284);

    auto ta1_z_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 285);

    auto ta1_z_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 286);

    auto ta1_z_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 287);

    auto ta1_z_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 288);

    auto ta1_z_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 289);

    auto ta1_z_zzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 290);

    auto ta1_z_zzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 291);

    auto ta1_z_zzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 292);

    auto ta1_z_zzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 293);

    auto ta1_z_zzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 294);

    auto ta1_z_zzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 295);

    auto ta1_z_zzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 296);

    auto ta1_z_zzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 297);

    auto ta1_z_zzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 298);

    auto ta1_z_zzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 299);

    // Set up components of auxiliary buffer : FF

    auto ta1_x_xxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff);

    auto ta1_x_xxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 1);

    auto ta1_x_xxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 2);

    auto ta1_x_xxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 3);

    auto ta1_x_xxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 4);

    auto ta1_x_xxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 5);

    auto ta1_x_xxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 6);

    auto ta1_x_xxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 7);

    auto ta1_x_xxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 8);

    auto ta1_x_xxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 9);

    auto ta1_x_xxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 10);

    auto ta1_x_xxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 11);

    auto ta1_x_xxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 12);

    auto ta1_x_xxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 13);

    auto ta1_x_xxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 14);

    auto ta1_x_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 15);

    auto ta1_x_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 20);

    auto ta1_x_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 21);

    auto ta1_x_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 22);

    auto ta1_x_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 23);

    auto ta1_x_xxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 24);

    auto ta1_x_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 25);

    auto ta1_x_xxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 27);

    auto ta1_x_xxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 28);

    auto ta1_x_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 29);

    auto ta1_x_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 31);

    auto ta1_x_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 33);

    auto ta1_x_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 34);

    auto ta1_x_xzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 50);

    auto ta1_x_xzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 51);

    auto ta1_x_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 52);

    auto ta1_x_xzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 53);

    auto ta1_x_xzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 54);

    auto ta1_x_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 55);

    auto ta1_x_yyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 60);

    auto ta1_x_yyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 61);

    auto ta1_x_yyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 62);

    auto ta1_x_yyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 63);

    auto ta1_x_yyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 64);

    auto ta1_x_yyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 65);

    auto ta1_x_yyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 66);

    auto ta1_x_yyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 67);

    auto ta1_x_yyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 68);

    auto ta1_x_yyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 69);

    auto ta1_x_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 82);

    auto ta1_x_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 84);

    auto ta1_x_yzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 85);

    auto ta1_x_yzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 87);

    auto ta1_x_yzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 88);

    auto ta1_x_yzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 89);

    auto ta1_x_zzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 90);

    auto ta1_x_zzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 91);

    auto ta1_x_zzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 92);

    auto ta1_x_zzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 93);

    auto ta1_x_zzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 94);

    auto ta1_x_zzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 95);

    auto ta1_x_zzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 96);

    auto ta1_x_zzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 97);

    auto ta1_x_zzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 98);

    auto ta1_x_zzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 99);

    auto ta1_y_xxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 100);

    auto ta1_y_xxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 101);

    auto ta1_y_xxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 102);

    auto ta1_y_xxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 103);

    auto ta1_y_xxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 104);

    auto ta1_y_xxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 105);

    auto ta1_y_xxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 106);

    auto ta1_y_xxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 107);

    auto ta1_y_xxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 108);

    auto ta1_y_xxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 109);

    auto ta1_y_xxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 111);

    auto ta1_y_xxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 113);

    auto ta1_y_xxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 114);

    auto ta1_y_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 131);

    auto ta1_y_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 133);

    auto ta1_y_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 134);

    auto ta1_y_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 136);

    auto ta1_y_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 137);

    auto ta1_y_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 138);

    auto ta1_y_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 152);

    auto ta1_y_xzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 154);

    auto ta1_y_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 155);

    auto ta1_y_xzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 157);

    auto ta1_y_xzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 158);

    auto ta1_y_xzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 159);

    auto ta1_y_yyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 160);

    auto ta1_y_yyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 161);

    auto ta1_y_yyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 162);

    auto ta1_y_yyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 163);

    auto ta1_y_yyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 164);

    auto ta1_y_yyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 165);

    auto ta1_y_yyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 166);

    auto ta1_y_yyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 167);

    auto ta1_y_yyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 168);

    auto ta1_y_yyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 169);

    auto ta1_y_yyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 171);

    auto ta1_y_yyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 172);

    auto ta1_y_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 173);

    auto ta1_y_yyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 174);

    auto ta1_y_yyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 175);

    auto ta1_y_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 176);

    auto ta1_y_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 177);

    auto ta1_y_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 178);

    auto ta1_y_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 179);

    auto ta1_y_yzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 181);

    auto ta1_y_yzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 183);

    auto ta1_y_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 184);

    auto ta1_y_yzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 186);

    auto ta1_y_yzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 187);

    auto ta1_y_yzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 188);

    auto ta1_y_zzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 190);

    auto ta1_y_zzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 191);

    auto ta1_y_zzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 192);

    auto ta1_y_zzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 193);

    auto ta1_y_zzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 194);

    auto ta1_y_zzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 195);

    auto ta1_y_zzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 196);

    auto ta1_y_zzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 197);

    auto ta1_y_zzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 198);

    auto ta1_y_zzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 199);

    auto ta1_z_xxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 200);

    auto ta1_z_xxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 201);

    auto ta1_z_xxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 202);

    auto ta1_z_xxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 203);

    auto ta1_z_xxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 204);

    auto ta1_z_xxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 205);

    auto ta1_z_xxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 206);

    auto ta1_z_xxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 207);

    auto ta1_z_xxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 208);

    auto ta1_z_xxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 209);

    auto ta1_z_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 222);

    auto ta1_z_xxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 224);

    auto ta1_z_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 225);

    auto ta1_z_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 231);

    auto ta1_z_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 233);

    auto ta1_z_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 234);

    auto ta1_z_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 236);

    auto ta1_z_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 237);

    auto ta1_z_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 238);

    auto ta1_z_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 252);

    auto ta1_z_xzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 254);

    auto ta1_z_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 255);

    auto ta1_z_xzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 257);

    auto ta1_z_xzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 258);

    auto ta1_z_xzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 259);

    auto ta1_z_yyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 260);

    auto ta1_z_yyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 261);

    auto ta1_z_yyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 262);

    auto ta1_z_yyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 263);

    auto ta1_z_yyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 264);

    auto ta1_z_yyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 265);

    auto ta1_z_yyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 266);

    auto ta1_z_yyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 267);

    auto ta1_z_yyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 268);

    auto ta1_z_yyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 269);

    auto ta1_z_yyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 272);

    auto ta1_z_yyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 274);

    auto ta1_z_yyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 275);

    auto ta1_z_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 277);

    auto ta1_z_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 278);

    auto ta1_z_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 279);

    auto ta1_z_yzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 281);

    auto ta1_z_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 282);

    auto ta1_z_yzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 283);

    auto ta1_z_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 284);

    auto ta1_z_yzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 285);

    auto ta1_z_yzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 286);

    auto ta1_z_yzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 287);

    auto ta1_z_yzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 288);

    auto ta1_z_yzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 289);

    auto ta1_z_zzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 290);

    auto ta1_z_zzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 291);

    auto ta1_z_zzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 292);

    auto ta1_z_zzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 293);

    auto ta1_z_zzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 294);

    auto ta1_z_zzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 295);

    auto ta1_z_zzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 296);

    auto ta1_z_zzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 297);

    auto ta1_z_zzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 298);

    auto ta1_z_zzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 299);

    // Set up components of auxiliary buffer : FG

    auto ta_xxx_xxxx_1 = pbuffer.data(idx_npot_1_fg);

    auto ta_xxx_xxxy_1 = pbuffer.data(idx_npot_1_fg + 1);

    auto ta_xxx_xxxz_1 = pbuffer.data(idx_npot_1_fg + 2);

    auto ta_xxx_xxyy_1 = pbuffer.data(idx_npot_1_fg + 3);

    auto ta_xxx_xxyz_1 = pbuffer.data(idx_npot_1_fg + 4);

    auto ta_xxx_xxzz_1 = pbuffer.data(idx_npot_1_fg + 5);

    auto ta_xxx_xyyy_1 = pbuffer.data(idx_npot_1_fg + 6);

    auto ta_xxx_xyyz_1 = pbuffer.data(idx_npot_1_fg + 7);

    auto ta_xxx_xyzz_1 = pbuffer.data(idx_npot_1_fg + 8);

    auto ta_xxx_xzzz_1 = pbuffer.data(idx_npot_1_fg + 9);

    auto ta_xxx_yyyy_1 = pbuffer.data(idx_npot_1_fg + 10);

    auto ta_xxx_yyyz_1 = pbuffer.data(idx_npot_1_fg + 11);

    auto ta_xxx_yyzz_1 = pbuffer.data(idx_npot_1_fg + 12);

    auto ta_xxx_yzzz_1 = pbuffer.data(idx_npot_1_fg + 13);

    auto ta_xxx_zzzz_1 = pbuffer.data(idx_npot_1_fg + 14);

    auto ta_xxy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 15);

    auto ta_xxy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 16);

    auto ta_xxy_xxxz_1 = pbuffer.data(idx_npot_1_fg + 17);

    auto ta_xxy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 18);

    auto ta_xxy_xxzz_1 = pbuffer.data(idx_npot_1_fg + 20);

    auto ta_xxy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 21);

    auto ta_xxy_xzzz_1 = pbuffer.data(idx_npot_1_fg + 24);

    auto ta_xxy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 25);

    auto ta_xxz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 30);

    auto ta_xxz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 31);

    auto ta_xxz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 32);

    auto ta_xxz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 33);

    auto ta_xxz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 35);

    auto ta_xxz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 36);

    auto ta_xxz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 39);

    auto ta_xxz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 44);

    auto ta_xyy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 45);

    auto ta_xyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 46);

    auto ta_xyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 48);

    auto ta_xyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 51);

    auto ta_xyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 55);

    auto ta_xyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 56);

    auto ta_xyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 57);

    auto ta_xyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 58);

    auto ta_xzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 75);

    auto ta_xzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 77);

    auto ta_xzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 80);

    auto ta_xzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 84);

    auto ta_xzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 86);

    auto ta_xzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 87);

    auto ta_xzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 88);

    auto ta_xzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 89);

    auto ta_yyy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 90);

    auto ta_yyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 91);

    auto ta_yyy_xxxz_1 = pbuffer.data(idx_npot_1_fg + 92);

    auto ta_yyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 93);

    auto ta_yyy_xxyz_1 = pbuffer.data(idx_npot_1_fg + 94);

    auto ta_yyy_xxzz_1 = pbuffer.data(idx_npot_1_fg + 95);

    auto ta_yyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 96);

    auto ta_yyy_xyyz_1 = pbuffer.data(idx_npot_1_fg + 97);

    auto ta_yyy_xyzz_1 = pbuffer.data(idx_npot_1_fg + 98);

    auto ta_yyy_xzzz_1 = pbuffer.data(idx_npot_1_fg + 99);

    auto ta_yyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 100);

    auto ta_yyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 101);

    auto ta_yyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 102);

    auto ta_yyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 103);

    auto ta_yyy_zzzz_1 = pbuffer.data(idx_npot_1_fg + 104);

    auto ta_yyz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 106);

    auto ta_yyz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 108);

    auto ta_yyz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 111);

    auto ta_yyz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 115);

    auto ta_yyz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 116);

    auto ta_yyz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 117);

    auto ta_yyz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 118);

    auto ta_yyz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 119);

    auto ta_yzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 122);

    auto ta_yzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 125);

    auto ta_yzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 129);

    auto ta_yzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 130);

    auto ta_yzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 131);

    auto ta_yzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 132);

    auto ta_yzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 133);

    auto ta_yzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 134);

    auto ta_zzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 135);

    auto ta_zzz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 136);

    auto ta_zzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 137);

    auto ta_zzz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 138);

    auto ta_zzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 139);

    auto ta_zzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 140);

    auto ta_zzz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 141);

    auto ta_zzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 142);

    auto ta_zzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 143);

    auto ta_zzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 144);

    auto ta_zzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 145);

    auto ta_zzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 146);

    auto ta_zzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 147);

    auto ta_zzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 148);

    auto ta_zzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 149);

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

    auto ta1_x_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 25);

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

    auto ta1_x_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 41);

    auto ta1_x_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 42);

    auto ta1_x_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 43);

    auto ta1_x_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 44);

    auto ta1_x_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 45);

    auto ta1_x_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 46);

    auto ta1_x_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 47);

    auto ta1_x_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 48);

    auto ta1_x_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 49);

    auto ta1_x_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 50);

    auto ta1_x_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 51);

    auto ta1_x_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 52);

    auto ta1_x_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 53);

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

    auto ta1_x_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 79);

    auto ta1_x_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 80);

    auto ta1_x_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 81);

    auto ta1_x_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 82);

    auto ta1_x_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 83);

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

    auto ta1_x_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 116);

    auto ta1_x_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 117);

    auto ta1_x_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 118);

    auto ta1_x_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 119);

    auto ta1_x_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 120);

    auto ta1_x_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 122);

    auto ta1_x_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 124);

    auto ta1_x_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 125);

    auto ta1_x_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 127);

    auto ta1_x_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 128);

    auto ta1_x_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 129);

    auto ta1_x_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 130);

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

    auto ta1_y_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 169);

    auto ta1_y_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 170);

    auto ta1_y_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 171);

    auto ta1_y_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 172);

    auto ta1_y_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 173);

    auto ta1_y_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 174);

    auto ta1_y_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 175);

    auto ta1_y_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 176);

    auto ta1_y_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 177);

    auto ta1_y_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 178);

    auto ta1_y_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 180);

    auto ta1_y_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 181);

    auto ta1_y_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 182);

    auto ta1_y_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 183);

    auto ta1_y_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 185);

    auto ta1_y_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 186);

    auto ta1_y_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 189);

    auto ta1_y_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 191);

    auto ta1_y_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 192);

    auto ta1_y_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 193);

    auto ta1_y_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 194);

    auto ta1_y_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 195);

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

    auto ta1_y_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 225);

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

    auto ta1_y_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 257);

    auto ta1_y_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 258);

    auto ta1_y_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 259);

    auto ta1_y_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 260);

    auto ta1_y_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 261);

    auto ta1_y_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 262);

    auto ta1_y_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 263);

    auto ta1_y_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 264);

    auto ta1_y_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 265);

    auto ta1_y_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 266);

    auto ta1_y_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 267);

    auto ta1_y_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 268);

    auto ta1_y_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 269);

    auto ta1_y_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 271);

    auto ta1_y_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 272);

    auto ta1_y_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 273);

    auto ta1_y_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 274);

    auto ta1_y_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 275);

    auto ta1_y_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 276);

    auto ta1_y_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 277);

    auto ta1_y_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 278);

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

    auto ta1_z_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 316);

    auto ta1_z_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 317);

    auto ta1_z_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 318);

    auto ta1_z_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 320);

    auto ta1_z_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 321);

    auto ta1_z_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 324);

    auto ta1_z_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 325);

    auto ta1_z_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 326);

    auto ta1_z_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 327);

    auto ta1_z_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 328);

    auto ta1_z_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 330);

    auto ta1_z_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 331);

    auto ta1_z_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 332);

    auto ta1_z_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 333);

    auto ta1_z_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 334);

    auto ta1_z_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 335);

    auto ta1_z_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 336);

    auto ta1_z_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 337);

    auto ta1_z_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 338);

    auto ta1_z_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 339);

    auto ta1_z_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 341);

    auto ta1_z_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 342);

    auto ta1_z_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 343);

    auto ta1_z_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 344);

    auto ta1_z_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 345);

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

    auto ta1_z_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 375);

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

    auto ta1_z_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 409);

    auto ta1_z_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 410);

    auto ta1_z_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 411);

    auto ta1_z_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 412);

    auto ta1_z_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 413);

    auto ta1_z_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 414);

    auto ta1_z_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 415);

    auto ta1_z_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 416);

    auto ta1_z_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 417);

    auto ta1_z_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 418);

    auto ta1_z_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 419);

    auto ta1_z_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 420);

    auto ta1_z_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 421);

    auto ta1_z_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 422);

    auto ta1_z_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 423);

    auto ta1_z_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 424);

    auto ta1_z_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 425);

    auto ta1_z_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 426);

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

    auto ta1_x_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 25);

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

    auto ta1_x_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 41);

    auto ta1_x_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 42);

    auto ta1_x_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 43);

    auto ta1_x_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 44);

    auto ta1_x_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 45);

    auto ta1_x_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 46);

    auto ta1_x_xyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 47);

    auto ta1_x_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 48);

    auto ta1_x_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 49);

    auto ta1_x_xyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 50);

    auto ta1_x_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 51);

    auto ta1_x_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 52);

    auto ta1_x_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 53);

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

    auto ta1_x_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 79);

    auto ta1_x_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 80);

    auto ta1_x_xzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 81);

    auto ta1_x_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 82);

    auto ta1_x_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 83);

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

    auto ta1_x_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 116);

    auto ta1_x_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 117);

    auto ta1_x_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 118);

    auto ta1_x_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 119);

    auto ta1_x_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 120);

    auto ta1_x_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 122);

    auto ta1_x_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 124);

    auto ta1_x_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 125);

    auto ta1_x_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 127);

    auto ta1_x_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 128);

    auto ta1_x_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 129);

    auto ta1_x_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 130);

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

    auto ta1_y_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 169);

    auto ta1_y_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 170);

    auto ta1_y_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 171);

    auto ta1_y_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 172);

    auto ta1_y_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 173);

    auto ta1_y_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 174);

    auto ta1_y_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 175);

    auto ta1_y_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 176);

    auto ta1_y_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 177);

    auto ta1_y_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 178);

    auto ta1_y_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 180);

    auto ta1_y_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 181);

    auto ta1_y_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 182);

    auto ta1_y_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 183);

    auto ta1_y_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 185);

    auto ta1_y_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 186);

    auto ta1_y_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 189);

    auto ta1_y_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 191);

    auto ta1_y_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 192);

    auto ta1_y_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 193);

    auto ta1_y_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 194);

    auto ta1_y_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 195);

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

    auto ta1_y_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 225);

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

    auto ta1_y_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 257);

    auto ta1_y_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 258);

    auto ta1_y_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 259);

    auto ta1_y_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 260);

    auto ta1_y_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 261);

    auto ta1_y_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 262);

    auto ta1_y_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 263);

    auto ta1_y_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 264);

    auto ta1_y_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 265);

    auto ta1_y_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 266);

    auto ta1_y_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 267);

    auto ta1_y_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 268);

    auto ta1_y_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 269);

    auto ta1_y_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 271);

    auto ta1_y_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 272);

    auto ta1_y_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 273);

    auto ta1_y_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 274);

    auto ta1_y_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 275);

    auto ta1_y_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 276);

    auto ta1_y_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 277);

    auto ta1_y_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 278);

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

    auto ta1_z_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 316);

    auto ta1_z_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 317);

    auto ta1_z_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 318);

    auto ta1_z_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 320);

    auto ta1_z_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 321);

    auto ta1_z_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 324);

    auto ta1_z_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 325);

    auto ta1_z_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 326);

    auto ta1_z_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 327);

    auto ta1_z_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 328);

    auto ta1_z_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 330);

    auto ta1_z_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 331);

    auto ta1_z_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 332);

    auto ta1_z_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 333);

    auto ta1_z_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 334);

    auto ta1_z_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 335);

    auto ta1_z_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 336);

    auto ta1_z_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 337);

    auto ta1_z_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 338);

    auto ta1_z_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 339);

    auto ta1_z_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 341);

    auto ta1_z_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 342);

    auto ta1_z_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 343);

    auto ta1_z_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 344);

    auto ta1_z_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 345);

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

    auto ta1_z_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 375);

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

    auto ta1_z_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 409);

    auto ta1_z_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 410);

    auto ta1_z_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 411);

    auto ta1_z_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 412);

    auto ta1_z_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 413);

    auto ta1_z_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 414);

    auto ta1_z_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 415);

    auto ta1_z_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 416);

    auto ta1_z_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 417);

    auto ta1_z_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 418);

    auto ta1_z_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 419);

    auto ta1_z_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 420);

    auto ta1_z_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 421);

    auto ta1_z_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 422);

    auto ta1_z_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 423);

    auto ta1_z_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 424);

    auto ta1_z_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 425);

    auto ta1_z_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 426);

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

    // Set up 0-15 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_x_xx_xxxx_0,   \
                             ta1_x_xx_xxxx_1,   \
                             ta1_x_xx_xxxy_0,   \
                             ta1_x_xx_xxxy_1,   \
                             ta1_x_xx_xxxz_0,   \
                             ta1_x_xx_xxxz_1,   \
                             ta1_x_xx_xxyy_0,   \
                             ta1_x_xx_xxyy_1,   \
                             ta1_x_xx_xxyz_0,   \
                             ta1_x_xx_xxyz_1,   \
                             ta1_x_xx_xxzz_0,   \
                             ta1_x_xx_xxzz_1,   \
                             ta1_x_xx_xyyy_0,   \
                             ta1_x_xx_xyyy_1,   \
                             ta1_x_xx_xyyz_0,   \
                             ta1_x_xx_xyyz_1,   \
                             ta1_x_xx_xyzz_0,   \
                             ta1_x_xx_xyzz_1,   \
                             ta1_x_xx_xzzz_0,   \
                             ta1_x_xx_xzzz_1,   \
                             ta1_x_xx_yyyy_0,   \
                             ta1_x_xx_yyyy_1,   \
                             ta1_x_xx_yyyz_0,   \
                             ta1_x_xx_yyyz_1,   \
                             ta1_x_xx_yyzz_0,   \
                             ta1_x_xx_yyzz_1,   \
                             ta1_x_xx_yzzz_0,   \
                             ta1_x_xx_yzzz_1,   \
                             ta1_x_xx_zzzz_0,   \
                             ta1_x_xx_zzzz_1,   \
                             ta1_x_xxx_xxx_0,   \
                             ta1_x_xxx_xxx_1,   \
                             ta1_x_xxx_xxxx_0,  \
                             ta1_x_xxx_xxxx_1,  \
                             ta1_x_xxx_xxxy_0,  \
                             ta1_x_xxx_xxxy_1,  \
                             ta1_x_xxx_xxxz_0,  \
                             ta1_x_xxx_xxxz_1,  \
                             ta1_x_xxx_xxy_0,   \
                             ta1_x_xxx_xxy_1,   \
                             ta1_x_xxx_xxyy_0,  \
                             ta1_x_xxx_xxyy_1,  \
                             ta1_x_xxx_xxyz_0,  \
                             ta1_x_xxx_xxyz_1,  \
                             ta1_x_xxx_xxz_0,   \
                             ta1_x_xxx_xxz_1,   \
                             ta1_x_xxx_xxzz_0,  \
                             ta1_x_xxx_xxzz_1,  \
                             ta1_x_xxx_xyy_0,   \
                             ta1_x_xxx_xyy_1,   \
                             ta1_x_xxx_xyyy_0,  \
                             ta1_x_xxx_xyyy_1,  \
                             ta1_x_xxx_xyyz_0,  \
                             ta1_x_xxx_xyyz_1,  \
                             ta1_x_xxx_xyz_0,   \
                             ta1_x_xxx_xyz_1,   \
                             ta1_x_xxx_xyzz_0,  \
                             ta1_x_xxx_xyzz_1,  \
                             ta1_x_xxx_xzz_0,   \
                             ta1_x_xxx_xzz_1,   \
                             ta1_x_xxx_xzzz_0,  \
                             ta1_x_xxx_xzzz_1,  \
                             ta1_x_xxx_yyy_0,   \
                             ta1_x_xxx_yyy_1,   \
                             ta1_x_xxx_yyyy_0,  \
                             ta1_x_xxx_yyyy_1,  \
                             ta1_x_xxx_yyyz_0,  \
                             ta1_x_xxx_yyyz_1,  \
                             ta1_x_xxx_yyz_0,   \
                             ta1_x_xxx_yyz_1,   \
                             ta1_x_xxx_yyzz_0,  \
                             ta1_x_xxx_yyzz_1,  \
                             ta1_x_xxx_yzz_0,   \
                             ta1_x_xxx_yzz_1,   \
                             ta1_x_xxx_yzzz_0,  \
                             ta1_x_xxx_yzzz_1,  \
                             ta1_x_xxx_zzz_0,   \
                             ta1_x_xxx_zzz_1,   \
                             ta1_x_xxx_zzzz_0,  \
                             ta1_x_xxx_zzzz_1,  \
                             ta1_x_xxxx_xxxx_0, \
                             ta1_x_xxxx_xxxy_0, \
                             ta1_x_xxxx_xxxz_0, \
                             ta1_x_xxxx_xxyy_0, \
                             ta1_x_xxxx_xxyz_0, \
                             ta1_x_xxxx_xxzz_0, \
                             ta1_x_xxxx_xyyy_0, \
                             ta1_x_xxxx_xyyz_0, \
                             ta1_x_xxxx_xyzz_0, \
                             ta1_x_xxxx_xzzz_0, \
                             ta1_x_xxxx_yyyy_0, \
                             ta1_x_xxxx_yyyz_0, \
                             ta1_x_xxxx_yyzz_0, \
                             ta1_x_xxxx_yzzz_0, \
                             ta1_x_xxxx_zzzz_0, \
                             ta_xxx_xxxx_1,     \
                             ta_xxx_xxxy_1,     \
                             ta_xxx_xxxz_1,     \
                             ta_xxx_xxyy_1,     \
                             ta_xxx_xxyz_1,     \
                             ta_xxx_xxzz_1,     \
                             ta_xxx_xyyy_1,     \
                             ta_xxx_xyyz_1,     \
                             ta_xxx_xyzz_1,     \
                             ta_xxx_xzzz_1,     \
                             ta_xxx_yyyy_1,     \
                             ta_xxx_yyyz_1,     \
                             ta_xxx_yyzz_1,     \
                             ta_xxx_yzzz_1,     \
                             ta_xxx_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_xxxx_0[i] = 3.0 * ta1_x_xx_xxxx_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxx_1[i] * fe_0 + 4.0 * ta1_x_xxx_xxx_0[i] * fe_0 -
                               4.0 * ta1_x_xxx_xxx_1[i] * fe_0 + ta_xxx_xxxx_1[i] + ta1_x_xxx_xxxx_0[i] * pa_x[i] - ta1_x_xxx_xxxx_1[i] * pc_x[i];

        ta1_x_xxxx_xxxy_0[i] = 3.0 * ta1_x_xx_xxxy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxy_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxy_0[i] * fe_0 -
                               3.0 * ta1_x_xxx_xxy_1[i] * fe_0 + ta_xxx_xxxy_1[i] + ta1_x_xxx_xxxy_0[i] * pa_x[i] - ta1_x_xxx_xxxy_1[i] * pc_x[i];

        ta1_x_xxxx_xxxz_0[i] = 3.0 * ta1_x_xx_xxxz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxz_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxz_0[i] * fe_0 -
                               3.0 * ta1_x_xxx_xxz_1[i] * fe_0 + ta_xxx_xxxz_1[i] + ta1_x_xxx_xxxz_0[i] * pa_x[i] - ta1_x_xxx_xxxz_1[i] * pc_x[i];

        ta1_x_xxxx_xxyy_0[i] = 3.0 * ta1_x_xx_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyy_0[i] * fe_0 -
                               2.0 * ta1_x_xxx_xyy_1[i] * fe_0 + ta_xxx_xxyy_1[i] + ta1_x_xxx_xxyy_0[i] * pa_x[i] - ta1_x_xxx_xxyy_1[i] * pc_x[i];

        ta1_x_xxxx_xxyz_0[i] = 3.0 * ta1_x_xx_xxyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyz_0[i] * fe_0 -
                               2.0 * ta1_x_xxx_xyz_1[i] * fe_0 + ta_xxx_xxyz_1[i] + ta1_x_xxx_xxyz_0[i] * pa_x[i] - ta1_x_xxx_xxyz_1[i] * pc_x[i];

        ta1_x_xxxx_xxzz_0[i] = 3.0 * ta1_x_xx_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xzz_0[i] * fe_0 -
                               2.0 * ta1_x_xxx_xzz_1[i] * fe_0 + ta_xxx_xxzz_1[i] + ta1_x_xxx_xxzz_0[i] * pa_x[i] - ta1_x_xxx_xxzz_1[i] * pc_x[i];

        ta1_x_xxxx_xyyy_0[i] = 3.0 * ta1_x_xx_xyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyy_1[i] * fe_0 + ta1_x_xxx_yyy_0[i] * fe_0 -
                               ta1_x_xxx_yyy_1[i] * fe_0 + ta_xxx_xyyy_1[i] + ta1_x_xxx_xyyy_0[i] * pa_x[i] - ta1_x_xxx_xyyy_1[i] * pc_x[i];

        ta1_x_xxxx_xyyz_0[i] = 3.0 * ta1_x_xx_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyz_1[i] * fe_0 + ta1_x_xxx_yyz_0[i] * fe_0 -
                               ta1_x_xxx_yyz_1[i] * fe_0 + ta_xxx_xyyz_1[i] + ta1_x_xxx_xyyz_0[i] * pa_x[i] - ta1_x_xxx_xyyz_1[i] * pc_x[i];

        ta1_x_xxxx_xyzz_0[i] = 3.0 * ta1_x_xx_xyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyzz_1[i] * fe_0 + ta1_x_xxx_yzz_0[i] * fe_0 -
                               ta1_x_xxx_yzz_1[i] * fe_0 + ta_xxx_xyzz_1[i] + ta1_x_xxx_xyzz_0[i] * pa_x[i] - ta1_x_xxx_xyzz_1[i] * pc_x[i];

        ta1_x_xxxx_xzzz_0[i] = 3.0 * ta1_x_xx_xzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xzzz_1[i] * fe_0 + ta1_x_xxx_zzz_0[i] * fe_0 -
                               ta1_x_xxx_zzz_1[i] * fe_0 + ta_xxx_xzzz_1[i] + ta1_x_xxx_xzzz_0[i] * pa_x[i] - ta1_x_xxx_xzzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyyy_0[i] = 3.0 * ta1_x_xx_yyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyy_1[i] * fe_0 + ta_xxx_yyyy_1[i] + ta1_x_xxx_yyyy_0[i] * pa_x[i] -
                               ta1_x_xxx_yyyy_1[i] * pc_x[i];

        ta1_x_xxxx_yyyz_0[i] = 3.0 * ta1_x_xx_yyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyz_1[i] * fe_0 + ta_xxx_yyyz_1[i] + ta1_x_xxx_yyyz_0[i] * pa_x[i] -
                               ta1_x_xxx_yyyz_1[i] * pc_x[i];

        ta1_x_xxxx_yyzz_0[i] = 3.0 * ta1_x_xx_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyzz_1[i] * fe_0 + ta_xxx_yyzz_1[i] + ta1_x_xxx_yyzz_0[i] * pa_x[i] -
                               ta1_x_xxx_yyzz_1[i] * pc_x[i];

        ta1_x_xxxx_yzzz_0[i] = 3.0 * ta1_x_xx_yzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yzzz_1[i] * fe_0 + ta_xxx_yzzz_1[i] + ta1_x_xxx_yzzz_0[i] * pa_x[i] -
                               ta1_x_xxx_yzzz_1[i] * pc_x[i];

        ta1_x_xxxx_zzzz_0[i] = 3.0 * ta1_x_xx_zzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_zzzz_1[i] * fe_0 + ta_xxx_zzzz_1[i] + ta1_x_xxx_zzzz_0[i] * pa_x[i] -
                               ta1_x_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : GG

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

    auto ta1_x_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 26);

    auto ta1_x_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 27);

    auto ta1_x_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 28);

    auto ta1_x_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 29);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_xxx_xxx_0,   \
                             ta1_x_xxx_xxx_1,   \
                             ta1_x_xxx_xxxx_0,  \
                             ta1_x_xxx_xxxx_1,  \
                             ta1_x_xxx_xxxy_0,  \
                             ta1_x_xxx_xxxy_1,  \
                             ta1_x_xxx_xxxz_0,  \
                             ta1_x_xxx_xxxz_1,  \
                             ta1_x_xxx_xxy_0,   \
                             ta1_x_xxx_xxy_1,   \
                             ta1_x_xxx_xxyy_0,  \
                             ta1_x_xxx_xxyy_1,  \
                             ta1_x_xxx_xxyz_0,  \
                             ta1_x_xxx_xxyz_1,  \
                             ta1_x_xxx_xxz_0,   \
                             ta1_x_xxx_xxz_1,   \
                             ta1_x_xxx_xxzz_0,  \
                             ta1_x_xxx_xxzz_1,  \
                             ta1_x_xxx_xyy_0,   \
                             ta1_x_xxx_xyy_1,   \
                             ta1_x_xxx_xyyy_0,  \
                             ta1_x_xxx_xyyy_1,  \
                             ta1_x_xxx_xyyz_0,  \
                             ta1_x_xxx_xyyz_1,  \
                             ta1_x_xxx_xyz_0,   \
                             ta1_x_xxx_xyz_1,   \
                             ta1_x_xxx_xyzz_0,  \
                             ta1_x_xxx_xyzz_1,  \
                             ta1_x_xxx_xzz_0,   \
                             ta1_x_xxx_xzz_1,   \
                             ta1_x_xxx_xzzz_0,  \
                             ta1_x_xxx_xzzz_1,  \
                             ta1_x_xxx_yyy_0,   \
                             ta1_x_xxx_yyy_1,   \
                             ta1_x_xxx_yyyy_0,  \
                             ta1_x_xxx_yyyy_1,  \
                             ta1_x_xxx_yyyz_0,  \
                             ta1_x_xxx_yyyz_1,  \
                             ta1_x_xxx_yyz_0,   \
                             ta1_x_xxx_yyz_1,   \
                             ta1_x_xxx_yyzz_0,  \
                             ta1_x_xxx_yyzz_1,  \
                             ta1_x_xxx_yzz_0,   \
                             ta1_x_xxx_yzz_1,   \
                             ta1_x_xxx_yzzz_0,  \
                             ta1_x_xxx_yzzz_1,  \
                             ta1_x_xxx_zzz_0,   \
                             ta1_x_xxx_zzz_1,   \
                             ta1_x_xxx_zzzz_0,  \
                             ta1_x_xxx_zzzz_1,  \
                             ta1_x_xxxy_xxxx_0, \
                             ta1_x_xxxy_xxxy_0, \
                             ta1_x_xxxy_xxxz_0, \
                             ta1_x_xxxy_xxyy_0, \
                             ta1_x_xxxy_xxyz_0, \
                             ta1_x_xxxy_xxzz_0, \
                             ta1_x_xxxy_xyyy_0, \
                             ta1_x_xxxy_xyyz_0, \
                             ta1_x_xxxy_xyzz_0, \
                             ta1_x_xxxy_xzzz_0, \
                             ta1_x_xxxy_yyyy_0, \
                             ta1_x_xxxy_yyyz_0, \
                             ta1_x_xxxy_yyzz_0, \
                             ta1_x_xxxy_yzzz_0, \
                             ta1_x_xxxy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxy_xxxx_0[i] = ta1_x_xxx_xxxx_0[i] * pa_y[i] - ta1_x_xxx_xxxx_1[i] * pc_y[i];

        ta1_x_xxxy_xxxy_0[i] = ta1_x_xxx_xxx_0[i] * fe_0 - ta1_x_xxx_xxx_1[i] * fe_0 + ta1_x_xxx_xxxy_0[i] * pa_y[i] - ta1_x_xxx_xxxy_1[i] * pc_y[i];

        ta1_x_xxxy_xxxz_0[i] = ta1_x_xxx_xxxz_0[i] * pa_y[i] - ta1_x_xxx_xxxz_1[i] * pc_y[i];

        ta1_x_xxxy_xxyy_0[i] =
            2.0 * ta1_x_xxx_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxy_1[i] * fe_0 + ta1_x_xxx_xxyy_0[i] * pa_y[i] - ta1_x_xxx_xxyy_1[i] * pc_y[i];

        ta1_x_xxxy_xxyz_0[i] = ta1_x_xxx_xxz_0[i] * fe_0 - ta1_x_xxx_xxz_1[i] * fe_0 + ta1_x_xxx_xxyz_0[i] * pa_y[i] - ta1_x_xxx_xxyz_1[i] * pc_y[i];

        ta1_x_xxxy_xxzz_0[i] = ta1_x_xxx_xxzz_0[i] * pa_y[i] - ta1_x_xxx_xxzz_1[i] * pc_y[i];

        ta1_x_xxxy_xyyy_0[i] =
            3.0 * ta1_x_xxx_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxx_xyy_1[i] * fe_0 + ta1_x_xxx_xyyy_0[i] * pa_y[i] - ta1_x_xxx_xyyy_1[i] * pc_y[i];

        ta1_x_xxxy_xyyz_0[i] =
            2.0 * ta1_x_xxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xyz_1[i] * fe_0 + ta1_x_xxx_xyyz_0[i] * pa_y[i] - ta1_x_xxx_xyyz_1[i] * pc_y[i];

        ta1_x_xxxy_xyzz_0[i] = ta1_x_xxx_xzz_0[i] * fe_0 - ta1_x_xxx_xzz_1[i] * fe_0 + ta1_x_xxx_xyzz_0[i] * pa_y[i] - ta1_x_xxx_xyzz_1[i] * pc_y[i];

        ta1_x_xxxy_xzzz_0[i] = ta1_x_xxx_xzzz_0[i] * pa_y[i] - ta1_x_xxx_xzzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyyy_0[i] =
            4.0 * ta1_x_xxx_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyy_1[i] * fe_0 + ta1_x_xxx_yyyy_0[i] * pa_y[i] - ta1_x_xxx_yyyy_1[i] * pc_y[i];

        ta1_x_xxxy_yyyz_0[i] =
            3.0 * ta1_x_xxx_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxx_yyz_1[i] * fe_0 + ta1_x_xxx_yyyz_0[i] * pa_y[i] - ta1_x_xxx_yyyz_1[i] * pc_y[i];

        ta1_x_xxxy_yyzz_0[i] =
            2.0 * ta1_x_xxx_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yzz_1[i] * fe_0 + ta1_x_xxx_yyzz_0[i] * pa_y[i] - ta1_x_xxx_yyzz_1[i] * pc_y[i];

        ta1_x_xxxy_yzzz_0[i] = ta1_x_xxx_zzz_0[i] * fe_0 - ta1_x_xxx_zzz_1[i] * fe_0 + ta1_x_xxx_yzzz_0[i] * pa_y[i] - ta1_x_xxx_yzzz_1[i] * pc_y[i];

        ta1_x_xxxy_zzzz_0[i] = ta1_x_xxx_zzzz_0[i] * pa_y[i] - ta1_x_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_x_xxx_xxx_0,   \
                             ta1_x_xxx_xxx_1,   \
                             ta1_x_xxx_xxxx_0,  \
                             ta1_x_xxx_xxxx_1,  \
                             ta1_x_xxx_xxxy_0,  \
                             ta1_x_xxx_xxxy_1,  \
                             ta1_x_xxx_xxxz_0,  \
                             ta1_x_xxx_xxxz_1,  \
                             ta1_x_xxx_xxy_0,   \
                             ta1_x_xxx_xxy_1,   \
                             ta1_x_xxx_xxyy_0,  \
                             ta1_x_xxx_xxyy_1,  \
                             ta1_x_xxx_xxyz_0,  \
                             ta1_x_xxx_xxyz_1,  \
                             ta1_x_xxx_xxz_0,   \
                             ta1_x_xxx_xxz_1,   \
                             ta1_x_xxx_xxzz_0,  \
                             ta1_x_xxx_xxzz_1,  \
                             ta1_x_xxx_xyy_0,   \
                             ta1_x_xxx_xyy_1,   \
                             ta1_x_xxx_xyyy_0,  \
                             ta1_x_xxx_xyyy_1,  \
                             ta1_x_xxx_xyyz_0,  \
                             ta1_x_xxx_xyyz_1,  \
                             ta1_x_xxx_xyz_0,   \
                             ta1_x_xxx_xyz_1,   \
                             ta1_x_xxx_xyzz_0,  \
                             ta1_x_xxx_xyzz_1,  \
                             ta1_x_xxx_xzz_0,   \
                             ta1_x_xxx_xzz_1,   \
                             ta1_x_xxx_xzzz_0,  \
                             ta1_x_xxx_xzzz_1,  \
                             ta1_x_xxx_yyy_0,   \
                             ta1_x_xxx_yyy_1,   \
                             ta1_x_xxx_yyyy_0,  \
                             ta1_x_xxx_yyyy_1,  \
                             ta1_x_xxx_yyyz_0,  \
                             ta1_x_xxx_yyyz_1,  \
                             ta1_x_xxx_yyz_0,   \
                             ta1_x_xxx_yyz_1,   \
                             ta1_x_xxx_yyzz_0,  \
                             ta1_x_xxx_yyzz_1,  \
                             ta1_x_xxx_yzz_0,   \
                             ta1_x_xxx_yzz_1,   \
                             ta1_x_xxx_yzzz_0,  \
                             ta1_x_xxx_yzzz_1,  \
                             ta1_x_xxx_zzz_0,   \
                             ta1_x_xxx_zzz_1,   \
                             ta1_x_xxx_zzzz_0,  \
                             ta1_x_xxx_zzzz_1,  \
                             ta1_x_xxxz_xxxx_0, \
                             ta1_x_xxxz_xxxy_0, \
                             ta1_x_xxxz_xxxz_0, \
                             ta1_x_xxxz_xxyy_0, \
                             ta1_x_xxxz_xxyz_0, \
                             ta1_x_xxxz_xxzz_0, \
                             ta1_x_xxxz_xyyy_0, \
                             ta1_x_xxxz_xyyz_0, \
                             ta1_x_xxxz_xyzz_0, \
                             ta1_x_xxxz_xzzz_0, \
                             ta1_x_xxxz_yyyy_0, \
                             ta1_x_xxxz_yyyz_0, \
                             ta1_x_xxxz_yyzz_0, \
                             ta1_x_xxxz_yzzz_0, \
                             ta1_x_xxxz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxz_xxxx_0[i] = ta1_x_xxx_xxxx_0[i] * pa_z[i] - ta1_x_xxx_xxxx_1[i] * pc_z[i];

        ta1_x_xxxz_xxxy_0[i] = ta1_x_xxx_xxxy_0[i] * pa_z[i] - ta1_x_xxx_xxxy_1[i] * pc_z[i];

        ta1_x_xxxz_xxxz_0[i] = ta1_x_xxx_xxx_0[i] * fe_0 - ta1_x_xxx_xxx_1[i] * fe_0 + ta1_x_xxx_xxxz_0[i] * pa_z[i] - ta1_x_xxx_xxxz_1[i] * pc_z[i];

        ta1_x_xxxz_xxyy_0[i] = ta1_x_xxx_xxyy_0[i] * pa_z[i] - ta1_x_xxx_xxyy_1[i] * pc_z[i];

        ta1_x_xxxz_xxyz_0[i] = ta1_x_xxx_xxy_0[i] * fe_0 - ta1_x_xxx_xxy_1[i] * fe_0 + ta1_x_xxx_xxyz_0[i] * pa_z[i] - ta1_x_xxx_xxyz_1[i] * pc_z[i];

        ta1_x_xxxz_xxzz_0[i] =
            2.0 * ta1_x_xxx_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxz_1[i] * fe_0 + ta1_x_xxx_xxzz_0[i] * pa_z[i] - ta1_x_xxx_xxzz_1[i] * pc_z[i];

        ta1_x_xxxz_xyyy_0[i] = ta1_x_xxx_xyyy_0[i] * pa_z[i] - ta1_x_xxx_xyyy_1[i] * pc_z[i];

        ta1_x_xxxz_xyyz_0[i] = ta1_x_xxx_xyy_0[i] * fe_0 - ta1_x_xxx_xyy_1[i] * fe_0 + ta1_x_xxx_xyyz_0[i] * pa_z[i] - ta1_x_xxx_xyyz_1[i] * pc_z[i];

        ta1_x_xxxz_xyzz_0[i] =
            2.0 * ta1_x_xxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xyz_1[i] * fe_0 + ta1_x_xxx_xyzz_0[i] * pa_z[i] - ta1_x_xxx_xyzz_1[i] * pc_z[i];

        ta1_x_xxxz_xzzz_0[i] =
            3.0 * ta1_x_xxx_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xzz_1[i] * fe_0 + ta1_x_xxx_xzzz_0[i] * pa_z[i] - ta1_x_xxx_xzzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyyy_0[i] = ta1_x_xxx_yyyy_0[i] * pa_z[i] - ta1_x_xxx_yyyy_1[i] * pc_z[i];

        ta1_x_xxxz_yyyz_0[i] = ta1_x_xxx_yyy_0[i] * fe_0 - ta1_x_xxx_yyy_1[i] * fe_0 + ta1_x_xxx_yyyz_0[i] * pa_z[i] - ta1_x_xxx_yyyz_1[i] * pc_z[i];

        ta1_x_xxxz_yyzz_0[i] =
            2.0 * ta1_x_xxx_yyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yyz_1[i] * fe_0 + ta1_x_xxx_yyzz_0[i] * pa_z[i] - ta1_x_xxx_yyzz_1[i] * pc_z[i];

        ta1_x_xxxz_yzzz_0[i] =
            3.0 * ta1_x_xxx_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_yzz_1[i] * fe_0 + ta1_x_xxx_yzzz_0[i] * pa_z[i] - ta1_x_xxx_yzzz_1[i] * pc_z[i];

        ta1_x_xxxz_zzzz_0[i] =
            4.0 * ta1_x_xxx_zzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_zzz_1[i] * fe_0 + ta1_x_xxx_zzzz_0[i] * pa_z[i] - ta1_x_xxx_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xx_xxxx_0,   \
                             ta1_x_xx_xxxx_1,   \
                             ta1_x_xx_xxxy_0,   \
                             ta1_x_xx_xxxy_1,   \
                             ta1_x_xx_xxxz_0,   \
                             ta1_x_xx_xxxz_1,   \
                             ta1_x_xx_xxyy_0,   \
                             ta1_x_xx_xxyy_1,   \
                             ta1_x_xx_xxyz_0,   \
                             ta1_x_xx_xxyz_1,   \
                             ta1_x_xx_xxzz_0,   \
                             ta1_x_xx_xxzz_1,   \
                             ta1_x_xx_xyyy_0,   \
                             ta1_x_xx_xyyy_1,   \
                             ta1_x_xx_xyyz_0,   \
                             ta1_x_xx_xyyz_1,   \
                             ta1_x_xx_xyzz_0,   \
                             ta1_x_xx_xyzz_1,   \
                             ta1_x_xx_xzzz_0,   \
                             ta1_x_xx_xzzz_1,   \
                             ta1_x_xx_zzzz_0,   \
                             ta1_x_xx_zzzz_1,   \
                             ta1_x_xxy_xxx_0,   \
                             ta1_x_xxy_xxx_1,   \
                             ta1_x_xxy_xxxx_0,  \
                             ta1_x_xxy_xxxx_1,  \
                             ta1_x_xxy_xxxy_0,  \
                             ta1_x_xxy_xxxy_1,  \
                             ta1_x_xxy_xxxz_0,  \
                             ta1_x_xxy_xxxz_1,  \
                             ta1_x_xxy_xxy_0,   \
                             ta1_x_xxy_xxy_1,   \
                             ta1_x_xxy_xxyy_0,  \
                             ta1_x_xxy_xxyy_1,  \
                             ta1_x_xxy_xxyz_0,  \
                             ta1_x_xxy_xxyz_1,  \
                             ta1_x_xxy_xxz_0,   \
                             ta1_x_xxy_xxz_1,   \
                             ta1_x_xxy_xxzz_0,  \
                             ta1_x_xxy_xxzz_1,  \
                             ta1_x_xxy_xyy_0,   \
                             ta1_x_xxy_xyy_1,   \
                             ta1_x_xxy_xyyy_0,  \
                             ta1_x_xxy_xyyy_1,  \
                             ta1_x_xxy_xyyz_0,  \
                             ta1_x_xxy_xyyz_1,  \
                             ta1_x_xxy_xyz_0,   \
                             ta1_x_xxy_xyz_1,   \
                             ta1_x_xxy_xyzz_0,  \
                             ta1_x_xxy_xyzz_1,  \
                             ta1_x_xxy_xzz_0,   \
                             ta1_x_xxy_xzz_1,   \
                             ta1_x_xxy_xzzz_0,  \
                             ta1_x_xxy_xzzz_1,  \
                             ta1_x_xxy_zzzz_0,  \
                             ta1_x_xxy_zzzz_1,  \
                             ta1_x_xxyy_xxxx_0, \
                             ta1_x_xxyy_xxxy_0, \
                             ta1_x_xxyy_xxxz_0, \
                             ta1_x_xxyy_xxyy_0, \
                             ta1_x_xxyy_xxyz_0, \
                             ta1_x_xxyy_xxzz_0, \
                             ta1_x_xxyy_xyyy_0, \
                             ta1_x_xxyy_xyyz_0, \
                             ta1_x_xxyy_xyzz_0, \
                             ta1_x_xxyy_xzzz_0, \
                             ta1_x_xxyy_yyyy_0, \
                             ta1_x_xxyy_yyyz_0, \
                             ta1_x_xxyy_yyzz_0, \
                             ta1_x_xxyy_yzzz_0, \
                             ta1_x_xxyy_zzzz_0, \
                             ta1_x_xyy_yyyy_0,  \
                             ta1_x_xyy_yyyy_1,  \
                             ta1_x_xyy_yyyz_0,  \
                             ta1_x_xyy_yyyz_1,  \
                             ta1_x_xyy_yyzz_0,  \
                             ta1_x_xyy_yyzz_1,  \
                             ta1_x_xyy_yzzz_0,  \
                             ta1_x_xyy_yzzz_1,  \
                             ta1_x_yy_yyyy_0,   \
                             ta1_x_yy_yyyy_1,   \
                             ta1_x_yy_yyyz_0,   \
                             ta1_x_yy_yyyz_1,   \
                             ta1_x_yy_yyzz_0,   \
                             ta1_x_yy_yyzz_1,   \
                             ta1_x_yy_yzzz_0,   \
                             ta1_x_yy_yzzz_1,   \
                             ta_xyy_yyyy_1,     \
                             ta_xyy_yyyz_1,     \
                             ta_xyy_yyzz_1,     \
                             ta_xyy_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyy_xxxx_0[i] = ta1_x_xx_xxxx_0[i] * fe_0 - ta1_x_xx_xxxx_1[i] * fe_0 + ta1_x_xxy_xxxx_0[i] * pa_y[i] - ta1_x_xxy_xxxx_1[i] * pc_y[i];

        ta1_x_xxyy_xxxy_0[i] = ta1_x_xx_xxxy_0[i] * fe_0 - ta1_x_xx_xxxy_1[i] * fe_0 + ta1_x_xxy_xxx_0[i] * fe_0 - ta1_x_xxy_xxx_1[i] * fe_0 +
                               ta1_x_xxy_xxxy_0[i] * pa_y[i] - ta1_x_xxy_xxxy_1[i] * pc_y[i];

        ta1_x_xxyy_xxxz_0[i] = ta1_x_xx_xxxz_0[i] * fe_0 - ta1_x_xx_xxxz_1[i] * fe_0 + ta1_x_xxy_xxxz_0[i] * pa_y[i] - ta1_x_xxy_xxxz_1[i] * pc_y[i];

        ta1_x_xxyy_xxyy_0[i] = ta1_x_xx_xxyy_0[i] * fe_0 - ta1_x_xx_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxy_xxy_0[i] * fe_0 -
                               2.0 * ta1_x_xxy_xxy_1[i] * fe_0 + ta1_x_xxy_xxyy_0[i] * pa_y[i] - ta1_x_xxy_xxyy_1[i] * pc_y[i];

        ta1_x_xxyy_xxyz_0[i] = ta1_x_xx_xxyz_0[i] * fe_0 - ta1_x_xx_xxyz_1[i] * fe_0 + ta1_x_xxy_xxz_0[i] * fe_0 - ta1_x_xxy_xxz_1[i] * fe_0 +
                               ta1_x_xxy_xxyz_0[i] * pa_y[i] - ta1_x_xxy_xxyz_1[i] * pc_y[i];

        ta1_x_xxyy_xxzz_0[i] = ta1_x_xx_xxzz_0[i] * fe_0 - ta1_x_xx_xxzz_1[i] * fe_0 + ta1_x_xxy_xxzz_0[i] * pa_y[i] - ta1_x_xxy_xxzz_1[i] * pc_y[i];

        ta1_x_xxyy_xyyy_0[i] = ta1_x_xx_xyyy_0[i] * fe_0 - ta1_x_xx_xyyy_1[i] * fe_0 + 3.0 * ta1_x_xxy_xyy_0[i] * fe_0 -
                               3.0 * ta1_x_xxy_xyy_1[i] * fe_0 + ta1_x_xxy_xyyy_0[i] * pa_y[i] - ta1_x_xxy_xyyy_1[i] * pc_y[i];

        ta1_x_xxyy_xyyz_0[i] = ta1_x_xx_xyyz_0[i] * fe_0 - ta1_x_xx_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxy_xyz_0[i] * fe_0 -
                               2.0 * ta1_x_xxy_xyz_1[i] * fe_0 + ta1_x_xxy_xyyz_0[i] * pa_y[i] - ta1_x_xxy_xyyz_1[i] * pc_y[i];

        ta1_x_xxyy_xyzz_0[i] = ta1_x_xx_xyzz_0[i] * fe_0 - ta1_x_xx_xyzz_1[i] * fe_0 + ta1_x_xxy_xzz_0[i] * fe_0 - ta1_x_xxy_xzz_1[i] * fe_0 +
                               ta1_x_xxy_xyzz_0[i] * pa_y[i] - ta1_x_xxy_xyzz_1[i] * pc_y[i];

        ta1_x_xxyy_xzzz_0[i] = ta1_x_xx_xzzz_0[i] * fe_0 - ta1_x_xx_xzzz_1[i] * fe_0 + ta1_x_xxy_xzzz_0[i] * pa_y[i] - ta1_x_xxy_xzzz_1[i] * pc_y[i];

        ta1_x_xxyy_yyyy_0[i] =
            ta1_x_yy_yyyy_0[i] * fe_0 - ta1_x_yy_yyyy_1[i] * fe_0 + ta_xyy_yyyy_1[i] + ta1_x_xyy_yyyy_0[i] * pa_x[i] - ta1_x_xyy_yyyy_1[i] * pc_x[i];

        ta1_x_xxyy_yyyz_0[i] =
            ta1_x_yy_yyyz_0[i] * fe_0 - ta1_x_yy_yyyz_1[i] * fe_0 + ta_xyy_yyyz_1[i] + ta1_x_xyy_yyyz_0[i] * pa_x[i] - ta1_x_xyy_yyyz_1[i] * pc_x[i];

        ta1_x_xxyy_yyzz_0[i] =
            ta1_x_yy_yyzz_0[i] * fe_0 - ta1_x_yy_yyzz_1[i] * fe_0 + ta_xyy_yyzz_1[i] + ta1_x_xyy_yyzz_0[i] * pa_x[i] - ta1_x_xyy_yyzz_1[i] * pc_x[i];

        ta1_x_xxyy_yzzz_0[i] =
            ta1_x_yy_yzzz_0[i] * fe_0 - ta1_x_yy_yzzz_1[i] * fe_0 + ta_xyy_yzzz_1[i] + ta1_x_xyy_yzzz_0[i] * pa_x[i] - ta1_x_xyy_yzzz_1[i] * pc_x[i];

        ta1_x_xxyy_zzzz_0[i] = ta1_x_xx_zzzz_0[i] * fe_0 - ta1_x_xx_zzzz_1[i] * fe_0 + ta1_x_xxy_zzzz_0[i] * pa_y[i] - ta1_x_xxy_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto ta1_x_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 60);

    auto ta1_x_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 61);

    auto ta1_x_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 62);

    auto ta1_x_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 63);

    auto ta1_x_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 64);

    auto ta1_x_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 65);

    auto ta1_x_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 66);

    auto ta1_x_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 67);

    auto ta1_x_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 68);

    auto ta1_x_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 69);

    auto ta1_x_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 70);

    auto ta1_x_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 71);

    auto ta1_x_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 72);

    auto ta1_x_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 73);

    auto ta1_x_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 74);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xxy_xxxy_0,  \
                             ta1_x_xxy_xxxy_1,  \
                             ta1_x_xxy_xxyy_0,  \
                             ta1_x_xxy_xxyy_1,  \
                             ta1_x_xxy_xyyy_0,  \
                             ta1_x_xxy_xyyy_1,  \
                             ta1_x_xxy_yyyy_0,  \
                             ta1_x_xxy_yyyy_1,  \
                             ta1_x_xxyz_xxxx_0, \
                             ta1_x_xxyz_xxxy_0, \
                             ta1_x_xxyz_xxxz_0, \
                             ta1_x_xxyz_xxyy_0, \
                             ta1_x_xxyz_xxyz_0, \
                             ta1_x_xxyz_xxzz_0, \
                             ta1_x_xxyz_xyyy_0, \
                             ta1_x_xxyz_xyyz_0, \
                             ta1_x_xxyz_xyzz_0, \
                             ta1_x_xxyz_xzzz_0, \
                             ta1_x_xxyz_yyyy_0, \
                             ta1_x_xxyz_yyyz_0, \
                             ta1_x_xxyz_yyzz_0, \
                             ta1_x_xxyz_yzzz_0, \
                             ta1_x_xxyz_zzzz_0, \
                             ta1_x_xxz_xxxx_0,  \
                             ta1_x_xxz_xxxx_1,  \
                             ta1_x_xxz_xxxz_0,  \
                             ta1_x_xxz_xxxz_1,  \
                             ta1_x_xxz_xxyz_0,  \
                             ta1_x_xxz_xxyz_1,  \
                             ta1_x_xxz_xxz_0,   \
                             ta1_x_xxz_xxz_1,   \
                             ta1_x_xxz_xxzz_0,  \
                             ta1_x_xxz_xxzz_1,  \
                             ta1_x_xxz_xyyz_0,  \
                             ta1_x_xxz_xyyz_1,  \
                             ta1_x_xxz_xyz_0,   \
                             ta1_x_xxz_xyz_1,   \
                             ta1_x_xxz_xyzz_0,  \
                             ta1_x_xxz_xyzz_1,  \
                             ta1_x_xxz_xzz_0,   \
                             ta1_x_xxz_xzz_1,   \
                             ta1_x_xxz_xzzz_0,  \
                             ta1_x_xxz_xzzz_1,  \
                             ta1_x_xxz_yyyz_0,  \
                             ta1_x_xxz_yyyz_1,  \
                             ta1_x_xxz_yyz_0,   \
                             ta1_x_xxz_yyz_1,   \
                             ta1_x_xxz_yyzz_0,  \
                             ta1_x_xxz_yyzz_1,  \
                             ta1_x_xxz_yzz_0,   \
                             ta1_x_xxz_yzz_1,   \
                             ta1_x_xxz_yzzz_0,  \
                             ta1_x_xxz_yzzz_1,  \
                             ta1_x_xxz_zzz_0,   \
                             ta1_x_xxz_zzz_1,   \
                             ta1_x_xxz_zzzz_0,  \
                             ta1_x_xxz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyz_xxxx_0[i] = ta1_x_xxz_xxxx_0[i] * pa_y[i] - ta1_x_xxz_xxxx_1[i] * pc_y[i];

        ta1_x_xxyz_xxxy_0[i] = ta1_x_xxy_xxxy_0[i] * pa_z[i] - ta1_x_xxy_xxxy_1[i] * pc_z[i];

        ta1_x_xxyz_xxxz_0[i] = ta1_x_xxz_xxxz_0[i] * pa_y[i] - ta1_x_xxz_xxxz_1[i] * pc_y[i];

        ta1_x_xxyz_xxyy_0[i] = ta1_x_xxy_xxyy_0[i] * pa_z[i] - ta1_x_xxy_xxyy_1[i] * pc_z[i];

        ta1_x_xxyz_xxyz_0[i] = ta1_x_xxz_xxz_0[i] * fe_0 - ta1_x_xxz_xxz_1[i] * fe_0 + ta1_x_xxz_xxyz_0[i] * pa_y[i] - ta1_x_xxz_xxyz_1[i] * pc_y[i];

        ta1_x_xxyz_xxzz_0[i] = ta1_x_xxz_xxzz_0[i] * pa_y[i] - ta1_x_xxz_xxzz_1[i] * pc_y[i];

        ta1_x_xxyz_xyyy_0[i] = ta1_x_xxy_xyyy_0[i] * pa_z[i] - ta1_x_xxy_xyyy_1[i] * pc_z[i];

        ta1_x_xxyz_xyyz_0[i] =
            2.0 * ta1_x_xxz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyz_1[i] * fe_0 + ta1_x_xxz_xyyz_0[i] * pa_y[i] - ta1_x_xxz_xyyz_1[i] * pc_y[i];

        ta1_x_xxyz_xyzz_0[i] = ta1_x_xxz_xzz_0[i] * fe_0 - ta1_x_xxz_xzz_1[i] * fe_0 + ta1_x_xxz_xyzz_0[i] * pa_y[i] - ta1_x_xxz_xyzz_1[i] * pc_y[i];

        ta1_x_xxyz_xzzz_0[i] = ta1_x_xxz_xzzz_0[i] * pa_y[i] - ta1_x_xxz_xzzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyyy_0[i] = ta1_x_xxy_yyyy_0[i] * pa_z[i] - ta1_x_xxy_yyyy_1[i] * pc_z[i];

        ta1_x_xxyz_yyyz_0[i] =
            3.0 * ta1_x_xxz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxz_yyz_1[i] * fe_0 + ta1_x_xxz_yyyz_0[i] * pa_y[i] - ta1_x_xxz_yyyz_1[i] * pc_y[i];

        ta1_x_xxyz_yyzz_0[i] =
            2.0 * ta1_x_xxz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_yzz_1[i] * fe_0 + ta1_x_xxz_yyzz_0[i] * pa_y[i] - ta1_x_xxz_yyzz_1[i] * pc_y[i];

        ta1_x_xxyz_yzzz_0[i] = ta1_x_xxz_zzz_0[i] * fe_0 - ta1_x_xxz_zzz_1[i] * fe_0 + ta1_x_xxz_yzzz_0[i] * pa_y[i] - ta1_x_xxz_yzzz_1[i] * pc_y[i];

        ta1_x_xxyz_zzzz_0[i] = ta1_x_xxz_zzzz_0[i] * pa_y[i] - ta1_x_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_xx_xxxx_0,   \
                             ta1_x_xx_xxxx_1,   \
                             ta1_x_xx_xxxy_0,   \
                             ta1_x_xx_xxxy_1,   \
                             ta1_x_xx_xxxz_0,   \
                             ta1_x_xx_xxxz_1,   \
                             ta1_x_xx_xxyy_0,   \
                             ta1_x_xx_xxyy_1,   \
                             ta1_x_xx_xxyz_0,   \
                             ta1_x_xx_xxyz_1,   \
                             ta1_x_xx_xxzz_0,   \
                             ta1_x_xx_xxzz_1,   \
                             ta1_x_xx_xyyy_0,   \
                             ta1_x_xx_xyyy_1,   \
                             ta1_x_xx_xyyz_0,   \
                             ta1_x_xx_xyyz_1,   \
                             ta1_x_xx_xyzz_0,   \
                             ta1_x_xx_xyzz_1,   \
                             ta1_x_xx_xzzz_0,   \
                             ta1_x_xx_xzzz_1,   \
                             ta1_x_xx_yyyy_0,   \
                             ta1_x_xx_yyyy_1,   \
                             ta1_x_xxz_xxx_0,   \
                             ta1_x_xxz_xxx_1,   \
                             ta1_x_xxz_xxxx_0,  \
                             ta1_x_xxz_xxxx_1,  \
                             ta1_x_xxz_xxxy_0,  \
                             ta1_x_xxz_xxxy_1,  \
                             ta1_x_xxz_xxxz_0,  \
                             ta1_x_xxz_xxxz_1,  \
                             ta1_x_xxz_xxy_0,   \
                             ta1_x_xxz_xxy_1,   \
                             ta1_x_xxz_xxyy_0,  \
                             ta1_x_xxz_xxyy_1,  \
                             ta1_x_xxz_xxyz_0,  \
                             ta1_x_xxz_xxyz_1,  \
                             ta1_x_xxz_xxz_0,   \
                             ta1_x_xxz_xxz_1,   \
                             ta1_x_xxz_xxzz_0,  \
                             ta1_x_xxz_xxzz_1,  \
                             ta1_x_xxz_xyy_0,   \
                             ta1_x_xxz_xyy_1,   \
                             ta1_x_xxz_xyyy_0,  \
                             ta1_x_xxz_xyyy_1,  \
                             ta1_x_xxz_xyyz_0,  \
                             ta1_x_xxz_xyyz_1,  \
                             ta1_x_xxz_xyz_0,   \
                             ta1_x_xxz_xyz_1,   \
                             ta1_x_xxz_xyzz_0,  \
                             ta1_x_xxz_xyzz_1,  \
                             ta1_x_xxz_xzz_0,   \
                             ta1_x_xxz_xzz_1,   \
                             ta1_x_xxz_xzzz_0,  \
                             ta1_x_xxz_xzzz_1,  \
                             ta1_x_xxz_yyyy_0,  \
                             ta1_x_xxz_yyyy_1,  \
                             ta1_x_xxzz_xxxx_0, \
                             ta1_x_xxzz_xxxy_0, \
                             ta1_x_xxzz_xxxz_0, \
                             ta1_x_xxzz_xxyy_0, \
                             ta1_x_xxzz_xxyz_0, \
                             ta1_x_xxzz_xxzz_0, \
                             ta1_x_xxzz_xyyy_0, \
                             ta1_x_xxzz_xyyz_0, \
                             ta1_x_xxzz_xyzz_0, \
                             ta1_x_xxzz_xzzz_0, \
                             ta1_x_xxzz_yyyy_0, \
                             ta1_x_xxzz_yyyz_0, \
                             ta1_x_xxzz_yyzz_0, \
                             ta1_x_xxzz_yzzz_0, \
                             ta1_x_xxzz_zzzz_0, \
                             ta1_x_xzz_yyyz_0,  \
                             ta1_x_xzz_yyyz_1,  \
                             ta1_x_xzz_yyzz_0,  \
                             ta1_x_xzz_yyzz_1,  \
                             ta1_x_xzz_yzzz_0,  \
                             ta1_x_xzz_yzzz_1,  \
                             ta1_x_xzz_zzzz_0,  \
                             ta1_x_xzz_zzzz_1,  \
                             ta1_x_zz_yyyz_0,   \
                             ta1_x_zz_yyyz_1,   \
                             ta1_x_zz_yyzz_0,   \
                             ta1_x_zz_yyzz_1,   \
                             ta1_x_zz_yzzz_0,   \
                             ta1_x_zz_yzzz_1,   \
                             ta1_x_zz_zzzz_0,   \
                             ta1_x_zz_zzzz_1,   \
                             ta_xzz_yyyz_1,     \
                             ta_xzz_yyzz_1,     \
                             ta_xzz_yzzz_1,     \
                             ta_xzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzz_xxxx_0[i] = ta1_x_xx_xxxx_0[i] * fe_0 - ta1_x_xx_xxxx_1[i] * fe_0 + ta1_x_xxz_xxxx_0[i] * pa_z[i] - ta1_x_xxz_xxxx_1[i] * pc_z[i];

        ta1_x_xxzz_xxxy_0[i] = ta1_x_xx_xxxy_0[i] * fe_0 - ta1_x_xx_xxxy_1[i] * fe_0 + ta1_x_xxz_xxxy_0[i] * pa_z[i] - ta1_x_xxz_xxxy_1[i] * pc_z[i];

        ta1_x_xxzz_xxxz_0[i] = ta1_x_xx_xxxz_0[i] * fe_0 - ta1_x_xx_xxxz_1[i] * fe_0 + ta1_x_xxz_xxx_0[i] * fe_0 - ta1_x_xxz_xxx_1[i] * fe_0 +
                               ta1_x_xxz_xxxz_0[i] * pa_z[i] - ta1_x_xxz_xxxz_1[i] * pc_z[i];

        ta1_x_xxzz_xxyy_0[i] = ta1_x_xx_xxyy_0[i] * fe_0 - ta1_x_xx_xxyy_1[i] * fe_0 + ta1_x_xxz_xxyy_0[i] * pa_z[i] - ta1_x_xxz_xxyy_1[i] * pc_z[i];

        ta1_x_xxzz_xxyz_0[i] = ta1_x_xx_xxyz_0[i] * fe_0 - ta1_x_xx_xxyz_1[i] * fe_0 + ta1_x_xxz_xxy_0[i] * fe_0 - ta1_x_xxz_xxy_1[i] * fe_0 +
                               ta1_x_xxz_xxyz_0[i] * pa_z[i] - ta1_x_xxz_xxyz_1[i] * pc_z[i];

        ta1_x_xxzz_xxzz_0[i] = ta1_x_xx_xxzz_0[i] * fe_0 - ta1_x_xx_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xxz_0[i] * fe_0 -
                               2.0 * ta1_x_xxz_xxz_1[i] * fe_0 + ta1_x_xxz_xxzz_0[i] * pa_z[i] - ta1_x_xxz_xxzz_1[i] * pc_z[i];

        ta1_x_xxzz_xyyy_0[i] = ta1_x_xx_xyyy_0[i] * fe_0 - ta1_x_xx_xyyy_1[i] * fe_0 + ta1_x_xxz_xyyy_0[i] * pa_z[i] - ta1_x_xxz_xyyy_1[i] * pc_z[i];

        ta1_x_xxzz_xyyz_0[i] = ta1_x_xx_xyyz_0[i] * fe_0 - ta1_x_xx_xyyz_1[i] * fe_0 + ta1_x_xxz_xyy_0[i] * fe_0 - ta1_x_xxz_xyy_1[i] * fe_0 +
                               ta1_x_xxz_xyyz_0[i] * pa_z[i] - ta1_x_xxz_xyyz_1[i] * pc_z[i];

        ta1_x_xxzz_xyzz_0[i] = ta1_x_xx_xyzz_0[i] * fe_0 - ta1_x_xx_xyzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xyz_0[i] * fe_0 -
                               2.0 * ta1_x_xxz_xyz_1[i] * fe_0 + ta1_x_xxz_xyzz_0[i] * pa_z[i] - ta1_x_xxz_xyzz_1[i] * pc_z[i];

        ta1_x_xxzz_xzzz_0[i] = ta1_x_xx_xzzz_0[i] * fe_0 - ta1_x_xx_xzzz_1[i] * fe_0 + 3.0 * ta1_x_xxz_xzz_0[i] * fe_0 -
                               3.0 * ta1_x_xxz_xzz_1[i] * fe_0 + ta1_x_xxz_xzzz_0[i] * pa_z[i] - ta1_x_xxz_xzzz_1[i] * pc_z[i];

        ta1_x_xxzz_yyyy_0[i] = ta1_x_xx_yyyy_0[i] * fe_0 - ta1_x_xx_yyyy_1[i] * fe_0 + ta1_x_xxz_yyyy_0[i] * pa_z[i] - ta1_x_xxz_yyyy_1[i] * pc_z[i];

        ta1_x_xxzz_yyyz_0[i] =
            ta1_x_zz_yyyz_0[i] * fe_0 - ta1_x_zz_yyyz_1[i] * fe_0 + ta_xzz_yyyz_1[i] + ta1_x_xzz_yyyz_0[i] * pa_x[i] - ta1_x_xzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxzz_yyzz_0[i] =
            ta1_x_zz_yyzz_0[i] * fe_0 - ta1_x_zz_yyzz_1[i] * fe_0 + ta_xzz_yyzz_1[i] + ta1_x_xzz_yyzz_0[i] * pa_x[i] - ta1_x_xzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxzz_yzzz_0[i] =
            ta1_x_zz_yzzz_0[i] * fe_0 - ta1_x_zz_yzzz_1[i] * fe_0 + ta_xzz_yzzz_1[i] + ta1_x_xzz_yzzz_0[i] * pa_x[i] - ta1_x_xzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxzz_zzzz_0[i] =
            ta1_x_zz_zzzz_0[i] * fe_0 - ta1_x_zz_zzzz_1[i] * fe_0 + ta_xzz_zzzz_1[i] + ta1_x_xzz_zzzz_0[i] * pa_x[i] - ta1_x_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : GG

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

    auto ta1_x_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 104);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xy_xxxx_0,   \
                             ta1_x_xy_xxxx_1,   \
                             ta1_x_xy_xxxz_0,   \
                             ta1_x_xy_xxxz_1,   \
                             ta1_x_xy_xxzz_0,   \
                             ta1_x_xy_xxzz_1,   \
                             ta1_x_xy_xzzz_0,   \
                             ta1_x_xy_xzzz_1,   \
                             ta1_x_xyy_xxxx_0,  \
                             ta1_x_xyy_xxxx_1,  \
                             ta1_x_xyy_xxxz_0,  \
                             ta1_x_xyy_xxxz_1,  \
                             ta1_x_xyy_xxzz_0,  \
                             ta1_x_xyy_xxzz_1,  \
                             ta1_x_xyy_xzzz_0,  \
                             ta1_x_xyy_xzzz_1,  \
                             ta1_x_xyyy_xxxx_0, \
                             ta1_x_xyyy_xxxy_0, \
                             ta1_x_xyyy_xxxz_0, \
                             ta1_x_xyyy_xxyy_0, \
                             ta1_x_xyyy_xxyz_0, \
                             ta1_x_xyyy_xxzz_0, \
                             ta1_x_xyyy_xyyy_0, \
                             ta1_x_xyyy_xyyz_0, \
                             ta1_x_xyyy_xyzz_0, \
                             ta1_x_xyyy_xzzz_0, \
                             ta1_x_xyyy_yyyy_0, \
                             ta1_x_xyyy_yyyz_0, \
                             ta1_x_xyyy_yyzz_0, \
                             ta1_x_xyyy_yzzz_0, \
                             ta1_x_xyyy_zzzz_0, \
                             ta1_x_yyy_xxxy_0,  \
                             ta1_x_yyy_xxxy_1,  \
                             ta1_x_yyy_xxy_0,   \
                             ta1_x_yyy_xxy_1,   \
                             ta1_x_yyy_xxyy_0,  \
                             ta1_x_yyy_xxyy_1,  \
                             ta1_x_yyy_xxyz_0,  \
                             ta1_x_yyy_xxyz_1,  \
                             ta1_x_yyy_xyy_0,   \
                             ta1_x_yyy_xyy_1,   \
                             ta1_x_yyy_xyyy_0,  \
                             ta1_x_yyy_xyyy_1,  \
                             ta1_x_yyy_xyyz_0,  \
                             ta1_x_yyy_xyyz_1,  \
                             ta1_x_yyy_xyz_0,   \
                             ta1_x_yyy_xyz_1,   \
                             ta1_x_yyy_xyzz_0,  \
                             ta1_x_yyy_xyzz_1,  \
                             ta1_x_yyy_yyy_0,   \
                             ta1_x_yyy_yyy_1,   \
                             ta1_x_yyy_yyyy_0,  \
                             ta1_x_yyy_yyyy_1,  \
                             ta1_x_yyy_yyyz_0,  \
                             ta1_x_yyy_yyyz_1,  \
                             ta1_x_yyy_yyz_0,   \
                             ta1_x_yyy_yyz_1,   \
                             ta1_x_yyy_yyzz_0,  \
                             ta1_x_yyy_yyzz_1,  \
                             ta1_x_yyy_yzz_0,   \
                             ta1_x_yyy_yzz_1,   \
                             ta1_x_yyy_yzzz_0,  \
                             ta1_x_yyy_yzzz_1,  \
                             ta1_x_yyy_zzzz_0,  \
                             ta1_x_yyy_zzzz_1,  \
                             ta_yyy_xxxy_1,     \
                             ta_yyy_xxyy_1,     \
                             ta_yyy_xxyz_1,     \
                             ta_yyy_xyyy_1,     \
                             ta_yyy_xyyz_1,     \
                             ta_yyy_xyzz_1,     \
                             ta_yyy_yyyy_1,     \
                             ta_yyy_yyyz_1,     \
                             ta_yyy_yyzz_1,     \
                             ta_yyy_yzzz_1,     \
                             ta_yyy_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyy_xxxx_0[i] =
            2.0 * ta1_x_xy_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxx_1[i] * fe_0 + ta1_x_xyy_xxxx_0[i] * pa_y[i] - ta1_x_xyy_xxxx_1[i] * pc_y[i];

        ta1_x_xyyy_xxxy_0[i] = 3.0 * ta1_x_yyy_xxy_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxy_1[i] * fe_0 + ta_yyy_xxxy_1[i] + ta1_x_yyy_xxxy_0[i] * pa_x[i] -
                               ta1_x_yyy_xxxy_1[i] * pc_x[i];

        ta1_x_xyyy_xxxz_0[i] =
            2.0 * ta1_x_xy_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxz_1[i] * fe_0 + ta1_x_xyy_xxxz_0[i] * pa_y[i] - ta1_x_xyy_xxxz_1[i] * pc_y[i];

        ta1_x_xyyy_xxyy_0[i] = 2.0 * ta1_x_yyy_xyy_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyy_1[i] * fe_0 + ta_yyy_xxyy_1[i] + ta1_x_yyy_xxyy_0[i] * pa_x[i] -
                               ta1_x_yyy_xxyy_1[i] * pc_x[i];

        ta1_x_xyyy_xxyz_0[i] = 2.0 * ta1_x_yyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyz_1[i] * fe_0 + ta_yyy_xxyz_1[i] + ta1_x_yyy_xxyz_0[i] * pa_x[i] -
                               ta1_x_yyy_xxyz_1[i] * pc_x[i];

        ta1_x_xyyy_xxzz_0[i] =
            2.0 * ta1_x_xy_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxzz_1[i] * fe_0 + ta1_x_xyy_xxzz_0[i] * pa_y[i] - ta1_x_xyy_xxzz_1[i] * pc_y[i];

        ta1_x_xyyy_xyyy_0[i] =
            ta1_x_yyy_yyy_0[i] * fe_0 - ta1_x_yyy_yyy_1[i] * fe_0 + ta_yyy_xyyy_1[i] + ta1_x_yyy_xyyy_0[i] * pa_x[i] - ta1_x_yyy_xyyy_1[i] * pc_x[i];

        ta1_x_xyyy_xyyz_0[i] =
            ta1_x_yyy_yyz_0[i] * fe_0 - ta1_x_yyy_yyz_1[i] * fe_0 + ta_yyy_xyyz_1[i] + ta1_x_yyy_xyyz_0[i] * pa_x[i] - ta1_x_yyy_xyyz_1[i] * pc_x[i];

        ta1_x_xyyy_xyzz_0[i] =
            ta1_x_yyy_yzz_0[i] * fe_0 - ta1_x_yyy_yzz_1[i] * fe_0 + ta_yyy_xyzz_1[i] + ta1_x_yyy_xyzz_0[i] * pa_x[i] - ta1_x_yyy_xyzz_1[i] * pc_x[i];

        ta1_x_xyyy_xzzz_0[i] =
            2.0 * ta1_x_xy_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xzzz_1[i] * fe_0 + ta1_x_xyy_xzzz_0[i] * pa_y[i] - ta1_x_xyy_xzzz_1[i] * pc_y[i];

        ta1_x_xyyy_yyyy_0[i] = ta_yyy_yyyy_1[i] + ta1_x_yyy_yyyy_0[i] * pa_x[i] - ta1_x_yyy_yyyy_1[i] * pc_x[i];

        ta1_x_xyyy_yyyz_0[i] = ta_yyy_yyyz_1[i] + ta1_x_yyy_yyyz_0[i] * pa_x[i] - ta1_x_yyy_yyyz_1[i] * pc_x[i];

        ta1_x_xyyy_yyzz_0[i] = ta_yyy_yyzz_1[i] + ta1_x_yyy_yyzz_0[i] * pa_x[i] - ta1_x_yyy_yyzz_1[i] * pc_x[i];

        ta1_x_xyyy_yzzz_0[i] = ta_yyy_yzzz_1[i] + ta1_x_yyy_yzzz_0[i] * pa_x[i] - ta1_x_yyy_yzzz_1[i] * pc_x[i];

        ta1_x_xyyy_zzzz_0[i] = ta_yyy_zzzz_1[i] + ta1_x_yyy_zzzz_0[i] * pa_x[i] - ta1_x_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto ta1_x_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 105);

    auto ta1_x_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 106);

    auto ta1_x_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 107);

    auto ta1_x_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 108);

    auto ta1_x_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 109);

    auto ta1_x_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 110);

    auto ta1_x_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 111);

    auto ta1_x_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 112);

    auto ta1_x_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 113);

    auto ta1_x_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 114);

    auto ta1_x_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 115);

    auto ta1_x_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 116);

    auto ta1_x_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 117);

    auto ta1_x_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 118);

    auto ta1_x_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 119);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xyy_xxxx_0,  \
                             ta1_x_xyy_xxxx_1,  \
                             ta1_x_xyy_xxxy_0,  \
                             ta1_x_xyy_xxxy_1,  \
                             ta1_x_xyy_xxy_0,   \
                             ta1_x_xyy_xxy_1,   \
                             ta1_x_xyy_xxyy_0,  \
                             ta1_x_xyy_xxyy_1,  \
                             ta1_x_xyy_xxyz_0,  \
                             ta1_x_xyy_xxyz_1,  \
                             ta1_x_xyy_xyy_0,   \
                             ta1_x_xyy_xyy_1,   \
                             ta1_x_xyy_xyyy_0,  \
                             ta1_x_xyy_xyyy_1,  \
                             ta1_x_xyy_xyyz_0,  \
                             ta1_x_xyy_xyyz_1,  \
                             ta1_x_xyy_xyz_0,   \
                             ta1_x_xyy_xyz_1,   \
                             ta1_x_xyy_xyzz_0,  \
                             ta1_x_xyy_xyzz_1,  \
                             ta1_x_xyy_yyyy_0,  \
                             ta1_x_xyy_yyyy_1,  \
                             ta1_x_xyyz_xxxx_0, \
                             ta1_x_xyyz_xxxy_0, \
                             ta1_x_xyyz_xxxz_0, \
                             ta1_x_xyyz_xxyy_0, \
                             ta1_x_xyyz_xxyz_0, \
                             ta1_x_xyyz_xxzz_0, \
                             ta1_x_xyyz_xyyy_0, \
                             ta1_x_xyyz_xyyz_0, \
                             ta1_x_xyyz_xyzz_0, \
                             ta1_x_xyyz_xzzz_0, \
                             ta1_x_xyyz_yyyy_0, \
                             ta1_x_xyyz_yyyz_0, \
                             ta1_x_xyyz_yyzz_0, \
                             ta1_x_xyyz_yzzz_0, \
                             ta1_x_xyyz_zzzz_0, \
                             ta1_x_xyz_xxxz_0,  \
                             ta1_x_xyz_xxxz_1,  \
                             ta1_x_xyz_xxzz_0,  \
                             ta1_x_xyz_xxzz_1,  \
                             ta1_x_xyz_xzzz_0,  \
                             ta1_x_xyz_xzzz_1,  \
                             ta1_x_xz_xxxz_0,   \
                             ta1_x_xz_xxxz_1,   \
                             ta1_x_xz_xxzz_0,   \
                             ta1_x_xz_xxzz_1,   \
                             ta1_x_xz_xzzz_0,   \
                             ta1_x_xz_xzzz_1,   \
                             ta1_x_yyz_yyyz_0,  \
                             ta1_x_yyz_yyyz_1,  \
                             ta1_x_yyz_yyzz_0,  \
                             ta1_x_yyz_yyzz_1,  \
                             ta1_x_yyz_yzzz_0,  \
                             ta1_x_yyz_yzzz_1,  \
                             ta1_x_yyz_zzzz_0,  \
                             ta1_x_yyz_zzzz_1,  \
                             ta_yyz_yyyz_1,     \
                             ta_yyz_yyzz_1,     \
                             ta_yyz_yzzz_1,     \
                             ta_yyz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyz_xxxx_0[i] = ta1_x_xyy_xxxx_0[i] * pa_z[i] - ta1_x_xyy_xxxx_1[i] * pc_z[i];

        ta1_x_xyyz_xxxy_0[i] = ta1_x_xyy_xxxy_0[i] * pa_z[i] - ta1_x_xyy_xxxy_1[i] * pc_z[i];

        ta1_x_xyyz_xxxz_0[i] = ta1_x_xz_xxxz_0[i] * fe_0 - ta1_x_xz_xxxz_1[i] * fe_0 + ta1_x_xyz_xxxz_0[i] * pa_y[i] - ta1_x_xyz_xxxz_1[i] * pc_y[i];

        ta1_x_xyyz_xxyy_0[i] = ta1_x_xyy_xxyy_0[i] * pa_z[i] - ta1_x_xyy_xxyy_1[i] * pc_z[i];

        ta1_x_xyyz_xxyz_0[i] = ta1_x_xyy_xxy_0[i] * fe_0 - ta1_x_xyy_xxy_1[i] * fe_0 + ta1_x_xyy_xxyz_0[i] * pa_z[i] - ta1_x_xyy_xxyz_1[i] * pc_z[i];

        ta1_x_xyyz_xxzz_0[i] = ta1_x_xz_xxzz_0[i] * fe_0 - ta1_x_xz_xxzz_1[i] * fe_0 + ta1_x_xyz_xxzz_0[i] * pa_y[i] - ta1_x_xyz_xxzz_1[i] * pc_y[i];

        ta1_x_xyyz_xyyy_0[i] = ta1_x_xyy_xyyy_0[i] * pa_z[i] - ta1_x_xyy_xyyy_1[i] * pc_z[i];

        ta1_x_xyyz_xyyz_0[i] = ta1_x_xyy_xyy_0[i] * fe_0 - ta1_x_xyy_xyy_1[i] * fe_0 + ta1_x_xyy_xyyz_0[i] * pa_z[i] - ta1_x_xyy_xyyz_1[i] * pc_z[i];

        ta1_x_xyyz_xyzz_0[i] =
            2.0 * ta1_x_xyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_xyz_1[i] * fe_0 + ta1_x_xyy_xyzz_0[i] * pa_z[i] - ta1_x_xyy_xyzz_1[i] * pc_z[i];

        ta1_x_xyyz_xzzz_0[i] = ta1_x_xz_xzzz_0[i] * fe_0 - ta1_x_xz_xzzz_1[i] * fe_0 + ta1_x_xyz_xzzz_0[i] * pa_y[i] - ta1_x_xyz_xzzz_1[i] * pc_y[i];

        ta1_x_xyyz_yyyy_0[i] = ta1_x_xyy_yyyy_0[i] * pa_z[i] - ta1_x_xyy_yyyy_1[i] * pc_z[i];

        ta1_x_xyyz_yyyz_0[i] = ta_yyz_yyyz_1[i] + ta1_x_yyz_yyyz_0[i] * pa_x[i] - ta1_x_yyz_yyyz_1[i] * pc_x[i];

        ta1_x_xyyz_yyzz_0[i] = ta_yyz_yyzz_1[i] + ta1_x_yyz_yyzz_0[i] * pa_x[i] - ta1_x_yyz_yyzz_1[i] * pc_x[i];

        ta1_x_xyyz_yzzz_0[i] = ta_yyz_yzzz_1[i] + ta1_x_yyz_yzzz_0[i] * pa_x[i] - ta1_x_yyz_yzzz_1[i] * pc_x[i];

        ta1_x_xyyz_zzzz_0[i] = ta_yyz_zzzz_1[i] + ta1_x_yyz_zzzz_0[i] * pa_x[i] - ta1_x_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto ta1_x_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 120);

    auto ta1_x_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 121);

    auto ta1_x_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 122);

    auto ta1_x_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 123);

    auto ta1_x_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 124);

    auto ta1_x_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 125);

    auto ta1_x_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 126);

    auto ta1_x_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 127);

    auto ta1_x_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 128);

    auto ta1_x_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 129);

    auto ta1_x_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 130);

    auto ta1_x_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 131);

    auto ta1_x_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 132);

    auto ta1_x_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 133);

    auto ta1_x_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 134);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xyzz_xxxx_0, \
                             ta1_x_xyzz_xxxy_0, \
                             ta1_x_xyzz_xxxz_0, \
                             ta1_x_xyzz_xxyy_0, \
                             ta1_x_xyzz_xxyz_0, \
                             ta1_x_xyzz_xxzz_0, \
                             ta1_x_xyzz_xyyy_0, \
                             ta1_x_xyzz_xyyz_0, \
                             ta1_x_xyzz_xyzz_0, \
                             ta1_x_xyzz_xzzz_0, \
                             ta1_x_xyzz_yyyy_0, \
                             ta1_x_xyzz_yyyz_0, \
                             ta1_x_xyzz_yyzz_0, \
                             ta1_x_xyzz_yzzz_0, \
                             ta1_x_xyzz_zzzz_0, \
                             ta1_x_xzz_xxx_0,   \
                             ta1_x_xzz_xxx_1,   \
                             ta1_x_xzz_xxxx_0,  \
                             ta1_x_xzz_xxxx_1,  \
                             ta1_x_xzz_xxxy_0,  \
                             ta1_x_xzz_xxxy_1,  \
                             ta1_x_xzz_xxxz_0,  \
                             ta1_x_xzz_xxxz_1,  \
                             ta1_x_xzz_xxy_0,   \
                             ta1_x_xzz_xxy_1,   \
                             ta1_x_xzz_xxyy_0,  \
                             ta1_x_xzz_xxyy_1,  \
                             ta1_x_xzz_xxyz_0,  \
                             ta1_x_xzz_xxyz_1,  \
                             ta1_x_xzz_xxz_0,   \
                             ta1_x_xzz_xxz_1,   \
                             ta1_x_xzz_xxzz_0,  \
                             ta1_x_xzz_xxzz_1,  \
                             ta1_x_xzz_xyy_0,   \
                             ta1_x_xzz_xyy_1,   \
                             ta1_x_xzz_xyyy_0,  \
                             ta1_x_xzz_xyyy_1,  \
                             ta1_x_xzz_xyyz_0,  \
                             ta1_x_xzz_xyyz_1,  \
                             ta1_x_xzz_xyz_0,   \
                             ta1_x_xzz_xyz_1,   \
                             ta1_x_xzz_xyzz_0,  \
                             ta1_x_xzz_xyzz_1,  \
                             ta1_x_xzz_xzz_0,   \
                             ta1_x_xzz_xzz_1,   \
                             ta1_x_xzz_xzzz_0,  \
                             ta1_x_xzz_xzzz_1,  \
                             ta1_x_xzz_zzzz_0,  \
                             ta1_x_xzz_zzzz_1,  \
                             ta1_x_yzz_yyyy_0,  \
                             ta1_x_yzz_yyyy_1,  \
                             ta1_x_yzz_yyyz_0,  \
                             ta1_x_yzz_yyyz_1,  \
                             ta1_x_yzz_yyzz_0,  \
                             ta1_x_yzz_yyzz_1,  \
                             ta1_x_yzz_yzzz_0,  \
                             ta1_x_yzz_yzzz_1,  \
                             ta_yzz_yyyy_1,     \
                             ta_yzz_yyyz_1,     \
                             ta_yzz_yyzz_1,     \
                             ta_yzz_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzz_xxxx_0[i] = ta1_x_xzz_xxxx_0[i] * pa_y[i] - ta1_x_xzz_xxxx_1[i] * pc_y[i];

        ta1_x_xyzz_xxxy_0[i] = ta1_x_xzz_xxx_0[i] * fe_0 - ta1_x_xzz_xxx_1[i] * fe_0 + ta1_x_xzz_xxxy_0[i] * pa_y[i] - ta1_x_xzz_xxxy_1[i] * pc_y[i];

        ta1_x_xyzz_xxxz_0[i] = ta1_x_xzz_xxxz_0[i] * pa_y[i] - ta1_x_xzz_xxxz_1[i] * pc_y[i];

        ta1_x_xyzz_xxyy_0[i] =
            2.0 * ta1_x_xzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xzz_xxy_1[i] * fe_0 + ta1_x_xzz_xxyy_0[i] * pa_y[i] - ta1_x_xzz_xxyy_1[i] * pc_y[i];

        ta1_x_xyzz_xxyz_0[i] = ta1_x_xzz_xxz_0[i] * fe_0 - ta1_x_xzz_xxz_1[i] * fe_0 + ta1_x_xzz_xxyz_0[i] * pa_y[i] - ta1_x_xzz_xxyz_1[i] * pc_y[i];

        ta1_x_xyzz_xxzz_0[i] = ta1_x_xzz_xxzz_0[i] * pa_y[i] - ta1_x_xzz_xxzz_1[i] * pc_y[i];

        ta1_x_xyzz_xyyy_0[i] =
            3.0 * ta1_x_xzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xyy_1[i] * fe_0 + ta1_x_xzz_xyyy_0[i] * pa_y[i] - ta1_x_xzz_xyyy_1[i] * pc_y[i];

        ta1_x_xyzz_xyyz_0[i] =
            2.0 * ta1_x_xzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xzz_xyz_1[i] * fe_0 + ta1_x_xzz_xyyz_0[i] * pa_y[i] - ta1_x_xzz_xyyz_1[i] * pc_y[i];

        ta1_x_xyzz_xyzz_0[i] = ta1_x_xzz_xzz_0[i] * fe_0 - ta1_x_xzz_xzz_1[i] * fe_0 + ta1_x_xzz_xyzz_0[i] * pa_y[i] - ta1_x_xzz_xyzz_1[i] * pc_y[i];

        ta1_x_xyzz_xzzz_0[i] = ta1_x_xzz_xzzz_0[i] * pa_y[i] - ta1_x_xzz_xzzz_1[i] * pc_y[i];

        ta1_x_xyzz_yyyy_0[i] = ta_yzz_yyyy_1[i] + ta1_x_yzz_yyyy_0[i] * pa_x[i] - ta1_x_yzz_yyyy_1[i] * pc_x[i];

        ta1_x_xyzz_yyyz_0[i] = ta_yzz_yyyz_1[i] + ta1_x_yzz_yyyz_0[i] * pa_x[i] - ta1_x_yzz_yyyz_1[i] * pc_x[i];

        ta1_x_xyzz_yyzz_0[i] = ta_yzz_yyzz_1[i] + ta1_x_yzz_yyzz_0[i] * pa_x[i] - ta1_x_yzz_yyzz_1[i] * pc_x[i];

        ta1_x_xyzz_yzzz_0[i] = ta_yzz_yzzz_1[i] + ta1_x_yzz_yzzz_0[i] * pa_x[i] - ta1_x_yzz_yzzz_1[i] * pc_x[i];

        ta1_x_xyzz_zzzz_0[i] = ta1_x_xzz_zzzz_0[i] * pa_y[i] - ta1_x_xzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : GG

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

    auto ta1_x_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 145);

    auto ta1_x_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 146);

    auto ta1_x_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 147);

    auto ta1_x_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 148);

    auto ta1_x_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 149);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_xz_xxxx_0,   \
                             ta1_x_xz_xxxx_1,   \
                             ta1_x_xz_xxxy_0,   \
                             ta1_x_xz_xxxy_1,   \
                             ta1_x_xz_xxyy_0,   \
                             ta1_x_xz_xxyy_1,   \
                             ta1_x_xz_xyyy_0,   \
                             ta1_x_xz_xyyy_1,   \
                             ta1_x_xzz_xxxx_0,  \
                             ta1_x_xzz_xxxx_1,  \
                             ta1_x_xzz_xxxy_0,  \
                             ta1_x_xzz_xxxy_1,  \
                             ta1_x_xzz_xxyy_0,  \
                             ta1_x_xzz_xxyy_1,  \
                             ta1_x_xzz_xyyy_0,  \
                             ta1_x_xzz_xyyy_1,  \
                             ta1_x_xzzz_xxxx_0, \
                             ta1_x_xzzz_xxxy_0, \
                             ta1_x_xzzz_xxxz_0, \
                             ta1_x_xzzz_xxyy_0, \
                             ta1_x_xzzz_xxyz_0, \
                             ta1_x_xzzz_xxzz_0, \
                             ta1_x_xzzz_xyyy_0, \
                             ta1_x_xzzz_xyyz_0, \
                             ta1_x_xzzz_xyzz_0, \
                             ta1_x_xzzz_xzzz_0, \
                             ta1_x_xzzz_yyyy_0, \
                             ta1_x_xzzz_yyyz_0, \
                             ta1_x_xzzz_yyzz_0, \
                             ta1_x_xzzz_yzzz_0, \
                             ta1_x_xzzz_zzzz_0, \
                             ta1_x_zzz_xxxz_0,  \
                             ta1_x_zzz_xxxz_1,  \
                             ta1_x_zzz_xxyz_0,  \
                             ta1_x_zzz_xxyz_1,  \
                             ta1_x_zzz_xxz_0,   \
                             ta1_x_zzz_xxz_1,   \
                             ta1_x_zzz_xxzz_0,  \
                             ta1_x_zzz_xxzz_1,  \
                             ta1_x_zzz_xyyz_0,  \
                             ta1_x_zzz_xyyz_1,  \
                             ta1_x_zzz_xyz_0,   \
                             ta1_x_zzz_xyz_1,   \
                             ta1_x_zzz_xyzz_0,  \
                             ta1_x_zzz_xyzz_1,  \
                             ta1_x_zzz_xzz_0,   \
                             ta1_x_zzz_xzz_1,   \
                             ta1_x_zzz_xzzz_0,  \
                             ta1_x_zzz_xzzz_1,  \
                             ta1_x_zzz_yyyy_0,  \
                             ta1_x_zzz_yyyy_1,  \
                             ta1_x_zzz_yyyz_0,  \
                             ta1_x_zzz_yyyz_1,  \
                             ta1_x_zzz_yyz_0,   \
                             ta1_x_zzz_yyz_1,   \
                             ta1_x_zzz_yyzz_0,  \
                             ta1_x_zzz_yyzz_1,  \
                             ta1_x_zzz_yzz_0,   \
                             ta1_x_zzz_yzz_1,   \
                             ta1_x_zzz_yzzz_0,  \
                             ta1_x_zzz_yzzz_1,  \
                             ta1_x_zzz_zzz_0,   \
                             ta1_x_zzz_zzz_1,   \
                             ta1_x_zzz_zzzz_0,  \
                             ta1_x_zzz_zzzz_1,  \
                             ta_zzz_xxxz_1,     \
                             ta_zzz_xxyz_1,     \
                             ta_zzz_xxzz_1,     \
                             ta_zzz_xyyz_1,     \
                             ta_zzz_xyzz_1,     \
                             ta_zzz_xzzz_1,     \
                             ta_zzz_yyyy_1,     \
                             ta_zzz_yyyz_1,     \
                             ta_zzz_yyzz_1,     \
                             ta_zzz_yzzz_1,     \
                             ta_zzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzz_xxxx_0[i] =
            2.0 * ta1_x_xz_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxx_1[i] * fe_0 + ta1_x_xzz_xxxx_0[i] * pa_z[i] - ta1_x_xzz_xxxx_1[i] * pc_z[i];

        ta1_x_xzzz_xxxy_0[i] =
            2.0 * ta1_x_xz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxy_1[i] * fe_0 + ta1_x_xzz_xxxy_0[i] * pa_z[i] - ta1_x_xzz_xxxy_1[i] * pc_z[i];

        ta1_x_xzzz_xxxz_0[i] = 3.0 * ta1_x_zzz_xxz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxz_1[i] * fe_0 + ta_zzz_xxxz_1[i] + ta1_x_zzz_xxxz_0[i] * pa_x[i] -
                               ta1_x_zzz_xxxz_1[i] * pc_x[i];

        ta1_x_xzzz_xxyy_0[i] =
            2.0 * ta1_x_xz_xxyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxyy_1[i] * fe_0 + ta1_x_xzz_xxyy_0[i] * pa_z[i] - ta1_x_xzz_xxyy_1[i] * pc_z[i];

        ta1_x_xzzz_xxyz_0[i] = 2.0 * ta1_x_zzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyz_1[i] * fe_0 + ta_zzz_xxyz_1[i] + ta1_x_zzz_xxyz_0[i] * pa_x[i] -
                               ta1_x_zzz_xxyz_1[i] * pc_x[i];

        ta1_x_xzzz_xxzz_0[i] = 2.0 * ta1_x_zzz_xzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xzz_1[i] * fe_0 + ta_zzz_xxzz_1[i] + ta1_x_zzz_xxzz_0[i] * pa_x[i] -
                               ta1_x_zzz_xxzz_1[i] * pc_x[i];

        ta1_x_xzzz_xyyy_0[i] =
            2.0 * ta1_x_xz_xyyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xyyy_1[i] * fe_0 + ta1_x_xzz_xyyy_0[i] * pa_z[i] - ta1_x_xzz_xyyy_1[i] * pc_z[i];

        ta1_x_xzzz_xyyz_0[i] =
            ta1_x_zzz_yyz_0[i] * fe_0 - ta1_x_zzz_yyz_1[i] * fe_0 + ta_zzz_xyyz_1[i] + ta1_x_zzz_xyyz_0[i] * pa_x[i] - ta1_x_zzz_xyyz_1[i] * pc_x[i];

        ta1_x_xzzz_xyzz_0[i] =
            ta1_x_zzz_yzz_0[i] * fe_0 - ta1_x_zzz_yzz_1[i] * fe_0 + ta_zzz_xyzz_1[i] + ta1_x_zzz_xyzz_0[i] * pa_x[i] - ta1_x_zzz_xyzz_1[i] * pc_x[i];

        ta1_x_xzzz_xzzz_0[i] =
            ta1_x_zzz_zzz_0[i] * fe_0 - ta1_x_zzz_zzz_1[i] * fe_0 + ta_zzz_xzzz_1[i] + ta1_x_zzz_xzzz_0[i] * pa_x[i] - ta1_x_zzz_xzzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyyy_0[i] = ta_zzz_yyyy_1[i] + ta1_x_zzz_yyyy_0[i] * pa_x[i] - ta1_x_zzz_yyyy_1[i] * pc_x[i];

        ta1_x_xzzz_yyyz_0[i] = ta_zzz_yyyz_1[i] + ta1_x_zzz_yyyz_0[i] * pa_x[i] - ta1_x_zzz_yyyz_1[i] * pc_x[i];

        ta1_x_xzzz_yyzz_0[i] = ta_zzz_yyzz_1[i] + ta1_x_zzz_yyzz_0[i] * pa_x[i] - ta1_x_zzz_yyzz_1[i] * pc_x[i];

        ta1_x_xzzz_yzzz_0[i] = ta_zzz_yzzz_1[i] + ta1_x_zzz_yzzz_0[i] * pa_x[i] - ta1_x_zzz_yzzz_1[i] * pc_x[i];

        ta1_x_xzzz_zzzz_0[i] = ta_zzz_zzzz_1[i] + ta1_x_zzz_zzzz_0[i] * pa_x[i] - ta1_x_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_yy_xxxx_0,   \
                             ta1_x_yy_xxxx_1,   \
                             ta1_x_yy_xxxy_0,   \
                             ta1_x_yy_xxxy_1,   \
                             ta1_x_yy_xxxz_0,   \
                             ta1_x_yy_xxxz_1,   \
                             ta1_x_yy_xxyy_0,   \
                             ta1_x_yy_xxyy_1,   \
                             ta1_x_yy_xxyz_0,   \
                             ta1_x_yy_xxyz_1,   \
                             ta1_x_yy_xxzz_0,   \
                             ta1_x_yy_xxzz_1,   \
                             ta1_x_yy_xyyy_0,   \
                             ta1_x_yy_xyyy_1,   \
                             ta1_x_yy_xyyz_0,   \
                             ta1_x_yy_xyyz_1,   \
                             ta1_x_yy_xyzz_0,   \
                             ta1_x_yy_xyzz_1,   \
                             ta1_x_yy_xzzz_0,   \
                             ta1_x_yy_xzzz_1,   \
                             ta1_x_yy_yyyy_0,   \
                             ta1_x_yy_yyyy_1,   \
                             ta1_x_yy_yyyz_0,   \
                             ta1_x_yy_yyyz_1,   \
                             ta1_x_yy_yyzz_0,   \
                             ta1_x_yy_yyzz_1,   \
                             ta1_x_yy_yzzz_0,   \
                             ta1_x_yy_yzzz_1,   \
                             ta1_x_yy_zzzz_0,   \
                             ta1_x_yy_zzzz_1,   \
                             ta1_x_yyy_xxx_0,   \
                             ta1_x_yyy_xxx_1,   \
                             ta1_x_yyy_xxxx_0,  \
                             ta1_x_yyy_xxxx_1,  \
                             ta1_x_yyy_xxxy_0,  \
                             ta1_x_yyy_xxxy_1,  \
                             ta1_x_yyy_xxxz_0,  \
                             ta1_x_yyy_xxxz_1,  \
                             ta1_x_yyy_xxy_0,   \
                             ta1_x_yyy_xxy_1,   \
                             ta1_x_yyy_xxyy_0,  \
                             ta1_x_yyy_xxyy_1,  \
                             ta1_x_yyy_xxyz_0,  \
                             ta1_x_yyy_xxyz_1,  \
                             ta1_x_yyy_xxz_0,   \
                             ta1_x_yyy_xxz_1,   \
                             ta1_x_yyy_xxzz_0,  \
                             ta1_x_yyy_xxzz_1,  \
                             ta1_x_yyy_xyy_0,   \
                             ta1_x_yyy_xyy_1,   \
                             ta1_x_yyy_xyyy_0,  \
                             ta1_x_yyy_xyyy_1,  \
                             ta1_x_yyy_xyyz_0,  \
                             ta1_x_yyy_xyyz_1,  \
                             ta1_x_yyy_xyz_0,   \
                             ta1_x_yyy_xyz_1,   \
                             ta1_x_yyy_xyzz_0,  \
                             ta1_x_yyy_xyzz_1,  \
                             ta1_x_yyy_xzz_0,   \
                             ta1_x_yyy_xzz_1,   \
                             ta1_x_yyy_xzzz_0,  \
                             ta1_x_yyy_xzzz_1,  \
                             ta1_x_yyy_yyy_0,   \
                             ta1_x_yyy_yyy_1,   \
                             ta1_x_yyy_yyyy_0,  \
                             ta1_x_yyy_yyyy_1,  \
                             ta1_x_yyy_yyyz_0,  \
                             ta1_x_yyy_yyyz_1,  \
                             ta1_x_yyy_yyz_0,   \
                             ta1_x_yyy_yyz_1,   \
                             ta1_x_yyy_yyzz_0,  \
                             ta1_x_yyy_yyzz_1,  \
                             ta1_x_yyy_yzz_0,   \
                             ta1_x_yyy_yzz_1,   \
                             ta1_x_yyy_yzzz_0,  \
                             ta1_x_yyy_yzzz_1,  \
                             ta1_x_yyy_zzz_0,   \
                             ta1_x_yyy_zzz_1,   \
                             ta1_x_yyy_zzzz_0,  \
                             ta1_x_yyy_zzzz_1,  \
                             ta1_x_yyyy_xxxx_0, \
                             ta1_x_yyyy_xxxy_0, \
                             ta1_x_yyyy_xxxz_0, \
                             ta1_x_yyyy_xxyy_0, \
                             ta1_x_yyyy_xxyz_0, \
                             ta1_x_yyyy_xxzz_0, \
                             ta1_x_yyyy_xyyy_0, \
                             ta1_x_yyyy_xyyz_0, \
                             ta1_x_yyyy_xyzz_0, \
                             ta1_x_yyyy_xzzz_0, \
                             ta1_x_yyyy_yyyy_0, \
                             ta1_x_yyyy_yyyz_0, \
                             ta1_x_yyyy_yyzz_0, \
                             ta1_x_yyyy_yzzz_0, \
                             ta1_x_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyy_xxxx_0[i] =
            3.0 * ta1_x_yy_xxxx_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxx_1[i] * fe_0 + ta1_x_yyy_xxxx_0[i] * pa_y[i] - ta1_x_yyy_xxxx_1[i] * pc_y[i];

        ta1_x_yyyy_xxxy_0[i] = 3.0 * ta1_x_yy_xxxy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxy_1[i] * fe_0 + ta1_x_yyy_xxx_0[i] * fe_0 -
                               ta1_x_yyy_xxx_1[i] * fe_0 + ta1_x_yyy_xxxy_0[i] * pa_y[i] - ta1_x_yyy_xxxy_1[i] * pc_y[i];

        ta1_x_yyyy_xxxz_0[i] =
            3.0 * ta1_x_yy_xxxz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxz_1[i] * fe_0 + ta1_x_yyy_xxxz_0[i] * pa_y[i] - ta1_x_yyy_xxxz_1[i] * pc_y[i];

        ta1_x_yyyy_xxyy_0[i] = 3.0 * ta1_x_yy_xxyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxy_0[i] * fe_0 -
                               2.0 * ta1_x_yyy_xxy_1[i] * fe_0 + ta1_x_yyy_xxyy_0[i] * pa_y[i] - ta1_x_yyy_xxyy_1[i] * pc_y[i];

        ta1_x_yyyy_xxyz_0[i] = 3.0 * ta1_x_yy_xxyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyz_1[i] * fe_0 + ta1_x_yyy_xxz_0[i] * fe_0 -
                               ta1_x_yyy_xxz_1[i] * fe_0 + ta1_x_yyy_xxyz_0[i] * pa_y[i] - ta1_x_yyy_xxyz_1[i] * pc_y[i];

        ta1_x_yyyy_xxzz_0[i] =
            3.0 * ta1_x_yy_xxzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxzz_1[i] * fe_0 + ta1_x_yyy_xxzz_0[i] * pa_y[i] - ta1_x_yyy_xxzz_1[i] * pc_y[i];

        ta1_x_yyyy_xyyy_0[i] = 3.0 * ta1_x_yy_xyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyy_1[i] * fe_0 + 3.0 * ta1_x_yyy_xyy_0[i] * fe_0 -
                               3.0 * ta1_x_yyy_xyy_1[i] * fe_0 + ta1_x_yyy_xyyy_0[i] * pa_y[i] - ta1_x_yyy_xyyy_1[i] * pc_y[i];

        ta1_x_yyyy_xyyz_0[i] = 3.0 * ta1_x_yy_xyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyz_0[i] * fe_0 -
                               2.0 * ta1_x_yyy_xyz_1[i] * fe_0 + ta1_x_yyy_xyyz_0[i] * pa_y[i] - ta1_x_yyy_xyyz_1[i] * pc_y[i];

        ta1_x_yyyy_xyzz_0[i] = 3.0 * ta1_x_yy_xyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyzz_1[i] * fe_0 + ta1_x_yyy_xzz_0[i] * fe_0 -
                               ta1_x_yyy_xzz_1[i] * fe_0 + ta1_x_yyy_xyzz_0[i] * pa_y[i] - ta1_x_yyy_xyzz_1[i] * pc_y[i];

        ta1_x_yyyy_xzzz_0[i] =
            3.0 * ta1_x_yy_xzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xzzz_1[i] * fe_0 + ta1_x_yyy_xzzz_0[i] * pa_y[i] - ta1_x_yyy_xzzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyyy_0[i] = 3.0 * ta1_x_yy_yyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyy_1[i] * fe_0 + 4.0 * ta1_x_yyy_yyy_0[i] * fe_0 -
                               4.0 * ta1_x_yyy_yyy_1[i] * fe_0 + ta1_x_yyy_yyyy_0[i] * pa_y[i] - ta1_x_yyy_yyyy_1[i] * pc_y[i];

        ta1_x_yyyy_yyyz_0[i] = 3.0 * ta1_x_yy_yyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yyy_yyz_0[i] * fe_0 -
                               3.0 * ta1_x_yyy_yyz_1[i] * fe_0 + ta1_x_yyy_yyyz_0[i] * pa_y[i] - ta1_x_yyy_yyyz_1[i] * pc_y[i];

        ta1_x_yyyy_yyzz_0[i] = 3.0 * ta1_x_yy_yyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_yzz_0[i] * fe_0 -
                               2.0 * ta1_x_yyy_yzz_1[i] * fe_0 + ta1_x_yyy_yyzz_0[i] * pa_y[i] - ta1_x_yyy_yyzz_1[i] * pc_y[i];

        ta1_x_yyyy_yzzz_0[i] = 3.0 * ta1_x_yy_yzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yzzz_1[i] * fe_0 + ta1_x_yyy_zzz_0[i] * fe_0 -
                               ta1_x_yyy_zzz_1[i] * fe_0 + ta1_x_yyy_yzzz_0[i] * pa_y[i] - ta1_x_yyy_yzzz_1[i] * pc_y[i];

        ta1_x_yyyy_zzzz_0[i] =
            3.0 * ta1_x_yy_zzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_zzzz_1[i] * fe_0 + ta1_x_yyy_zzzz_0[i] * pa_y[i] - ta1_x_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : GG

    auto ta1_x_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 165);

    auto ta1_x_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 166);

    auto ta1_x_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 167);

    auto ta1_x_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 168);

    auto ta1_x_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 169);

    auto ta1_x_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 170);

    auto ta1_x_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 171);

    auto ta1_x_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 172);

    auto ta1_x_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 173);

    auto ta1_x_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 174);

    auto ta1_x_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 175);

    auto ta1_x_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 176);

    auto ta1_x_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 177);

    auto ta1_x_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 178);

    auto ta1_x_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 179);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yyy_xxxx_0,  \
                             ta1_x_yyy_xxxx_1,  \
                             ta1_x_yyy_xxxy_0,  \
                             ta1_x_yyy_xxxy_1,  \
                             ta1_x_yyy_xxy_0,   \
                             ta1_x_yyy_xxy_1,   \
                             ta1_x_yyy_xxyy_0,  \
                             ta1_x_yyy_xxyy_1,  \
                             ta1_x_yyy_xxyz_0,  \
                             ta1_x_yyy_xxyz_1,  \
                             ta1_x_yyy_xyy_0,   \
                             ta1_x_yyy_xyy_1,   \
                             ta1_x_yyy_xyyy_0,  \
                             ta1_x_yyy_xyyy_1,  \
                             ta1_x_yyy_xyyz_0,  \
                             ta1_x_yyy_xyyz_1,  \
                             ta1_x_yyy_xyz_0,   \
                             ta1_x_yyy_xyz_1,   \
                             ta1_x_yyy_xyzz_0,  \
                             ta1_x_yyy_xyzz_1,  \
                             ta1_x_yyy_yyy_0,   \
                             ta1_x_yyy_yyy_1,   \
                             ta1_x_yyy_yyyy_0,  \
                             ta1_x_yyy_yyyy_1,  \
                             ta1_x_yyy_yyyz_0,  \
                             ta1_x_yyy_yyyz_1,  \
                             ta1_x_yyy_yyz_0,   \
                             ta1_x_yyy_yyz_1,   \
                             ta1_x_yyy_yyzz_0,  \
                             ta1_x_yyy_yyzz_1,  \
                             ta1_x_yyy_yzz_0,   \
                             ta1_x_yyy_yzz_1,   \
                             ta1_x_yyy_yzzz_0,  \
                             ta1_x_yyy_yzzz_1,  \
                             ta1_x_yyyz_xxxx_0, \
                             ta1_x_yyyz_xxxy_0, \
                             ta1_x_yyyz_xxxz_0, \
                             ta1_x_yyyz_xxyy_0, \
                             ta1_x_yyyz_xxyz_0, \
                             ta1_x_yyyz_xxzz_0, \
                             ta1_x_yyyz_xyyy_0, \
                             ta1_x_yyyz_xyyz_0, \
                             ta1_x_yyyz_xyzz_0, \
                             ta1_x_yyyz_xzzz_0, \
                             ta1_x_yyyz_yyyy_0, \
                             ta1_x_yyyz_yyyz_0, \
                             ta1_x_yyyz_yyzz_0, \
                             ta1_x_yyyz_yzzz_0, \
                             ta1_x_yyyz_zzzz_0, \
                             ta1_x_yyz_xxxz_0,  \
                             ta1_x_yyz_xxxz_1,  \
                             ta1_x_yyz_xxzz_0,  \
                             ta1_x_yyz_xxzz_1,  \
                             ta1_x_yyz_xzzz_0,  \
                             ta1_x_yyz_xzzz_1,  \
                             ta1_x_yyz_zzzz_0,  \
                             ta1_x_yyz_zzzz_1,  \
                             ta1_x_yz_xxxz_0,   \
                             ta1_x_yz_xxxz_1,   \
                             ta1_x_yz_xxzz_0,   \
                             ta1_x_yz_xxzz_1,   \
                             ta1_x_yz_xzzz_0,   \
                             ta1_x_yz_xzzz_1,   \
                             ta1_x_yz_zzzz_0,   \
                             ta1_x_yz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyz_xxxx_0[i] = ta1_x_yyy_xxxx_0[i] * pa_z[i] - ta1_x_yyy_xxxx_1[i] * pc_z[i];

        ta1_x_yyyz_xxxy_0[i] = ta1_x_yyy_xxxy_0[i] * pa_z[i] - ta1_x_yyy_xxxy_1[i] * pc_z[i];

        ta1_x_yyyz_xxxz_0[i] =
            2.0 * ta1_x_yz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxxz_1[i] * fe_0 + ta1_x_yyz_xxxz_0[i] * pa_y[i] - ta1_x_yyz_xxxz_1[i] * pc_y[i];

        ta1_x_yyyz_xxyy_0[i] = ta1_x_yyy_xxyy_0[i] * pa_z[i] - ta1_x_yyy_xxyy_1[i] * pc_z[i];

        ta1_x_yyyz_xxyz_0[i] = ta1_x_yyy_xxy_0[i] * fe_0 - ta1_x_yyy_xxy_1[i] * fe_0 + ta1_x_yyy_xxyz_0[i] * pa_z[i] - ta1_x_yyy_xxyz_1[i] * pc_z[i];

        ta1_x_yyyz_xxzz_0[i] =
            2.0 * ta1_x_yz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxzz_1[i] * fe_0 + ta1_x_yyz_xxzz_0[i] * pa_y[i] - ta1_x_yyz_xxzz_1[i] * pc_y[i];

        ta1_x_yyyz_xyyy_0[i] = ta1_x_yyy_xyyy_0[i] * pa_z[i] - ta1_x_yyy_xyyy_1[i] * pc_z[i];

        ta1_x_yyyz_xyyz_0[i] = ta1_x_yyy_xyy_0[i] * fe_0 - ta1_x_yyy_xyy_1[i] * fe_0 + ta1_x_yyy_xyyz_0[i] * pa_z[i] - ta1_x_yyy_xyyz_1[i] * pc_z[i];

        ta1_x_yyyz_xyzz_0[i] =
            2.0 * ta1_x_yyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyz_1[i] * fe_0 + ta1_x_yyy_xyzz_0[i] * pa_z[i] - ta1_x_yyy_xyzz_1[i] * pc_z[i];

        ta1_x_yyyz_xzzz_0[i] =
            2.0 * ta1_x_yz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xzzz_1[i] * fe_0 + ta1_x_yyz_xzzz_0[i] * pa_y[i] - ta1_x_yyz_xzzz_1[i] * pc_y[i];

        ta1_x_yyyz_yyyy_0[i] = ta1_x_yyy_yyyy_0[i] * pa_z[i] - ta1_x_yyy_yyyy_1[i] * pc_z[i];

        ta1_x_yyyz_yyyz_0[i] = ta1_x_yyy_yyy_0[i] * fe_0 - ta1_x_yyy_yyy_1[i] * fe_0 + ta1_x_yyy_yyyz_0[i] * pa_z[i] - ta1_x_yyy_yyyz_1[i] * pc_z[i];

        ta1_x_yyyz_yyzz_0[i] =
            2.0 * ta1_x_yyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_yyz_1[i] * fe_0 + ta1_x_yyy_yyzz_0[i] * pa_z[i] - ta1_x_yyy_yyzz_1[i] * pc_z[i];

        ta1_x_yyyz_yzzz_0[i] =
            3.0 * ta1_x_yyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_yzz_1[i] * fe_0 + ta1_x_yyy_yzzz_0[i] * pa_z[i] - ta1_x_yyy_yzzz_1[i] * pc_z[i];

        ta1_x_yyyz_zzzz_0[i] =
            2.0 * ta1_x_yz_zzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_zzzz_1[i] * fe_0 + ta1_x_yyz_zzzz_0[i] * pa_y[i] - ta1_x_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 180-195 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yy_xxxy_0,   \
                             ta1_x_yy_xxxy_1,   \
                             ta1_x_yy_xxyy_0,   \
                             ta1_x_yy_xxyy_1,   \
                             ta1_x_yy_xyyy_0,   \
                             ta1_x_yy_xyyy_1,   \
                             ta1_x_yy_yyyy_0,   \
                             ta1_x_yy_yyyy_1,   \
                             ta1_x_yyz_xxxy_0,  \
                             ta1_x_yyz_xxxy_1,  \
                             ta1_x_yyz_xxyy_0,  \
                             ta1_x_yyz_xxyy_1,  \
                             ta1_x_yyz_xyyy_0,  \
                             ta1_x_yyz_xyyy_1,  \
                             ta1_x_yyz_yyyy_0,  \
                             ta1_x_yyz_yyyy_1,  \
                             ta1_x_yyzz_xxxx_0, \
                             ta1_x_yyzz_xxxy_0, \
                             ta1_x_yyzz_xxxz_0, \
                             ta1_x_yyzz_xxyy_0, \
                             ta1_x_yyzz_xxyz_0, \
                             ta1_x_yyzz_xxzz_0, \
                             ta1_x_yyzz_xyyy_0, \
                             ta1_x_yyzz_xyyz_0, \
                             ta1_x_yyzz_xyzz_0, \
                             ta1_x_yyzz_xzzz_0, \
                             ta1_x_yyzz_yyyy_0, \
                             ta1_x_yyzz_yyyz_0, \
                             ta1_x_yyzz_yyzz_0, \
                             ta1_x_yyzz_yzzz_0, \
                             ta1_x_yyzz_zzzz_0, \
                             ta1_x_yzz_xxxx_0,  \
                             ta1_x_yzz_xxxx_1,  \
                             ta1_x_yzz_xxxz_0,  \
                             ta1_x_yzz_xxxz_1,  \
                             ta1_x_yzz_xxyz_0,  \
                             ta1_x_yzz_xxyz_1,  \
                             ta1_x_yzz_xxz_0,   \
                             ta1_x_yzz_xxz_1,   \
                             ta1_x_yzz_xxzz_0,  \
                             ta1_x_yzz_xxzz_1,  \
                             ta1_x_yzz_xyyz_0,  \
                             ta1_x_yzz_xyyz_1,  \
                             ta1_x_yzz_xyz_0,   \
                             ta1_x_yzz_xyz_1,   \
                             ta1_x_yzz_xyzz_0,  \
                             ta1_x_yzz_xyzz_1,  \
                             ta1_x_yzz_xzz_0,   \
                             ta1_x_yzz_xzz_1,   \
                             ta1_x_yzz_xzzz_0,  \
                             ta1_x_yzz_xzzz_1,  \
                             ta1_x_yzz_yyyz_0,  \
                             ta1_x_yzz_yyyz_1,  \
                             ta1_x_yzz_yyz_0,   \
                             ta1_x_yzz_yyz_1,   \
                             ta1_x_yzz_yyzz_0,  \
                             ta1_x_yzz_yyzz_1,  \
                             ta1_x_yzz_yzz_0,   \
                             ta1_x_yzz_yzz_1,   \
                             ta1_x_yzz_yzzz_0,  \
                             ta1_x_yzz_yzzz_1,  \
                             ta1_x_yzz_zzz_0,   \
                             ta1_x_yzz_zzz_1,   \
                             ta1_x_yzz_zzzz_0,  \
                             ta1_x_yzz_zzzz_1,  \
                             ta1_x_zz_xxxx_0,   \
                             ta1_x_zz_xxxx_1,   \
                             ta1_x_zz_xxxz_0,   \
                             ta1_x_zz_xxxz_1,   \
                             ta1_x_zz_xxyz_0,   \
                             ta1_x_zz_xxyz_1,   \
                             ta1_x_zz_xxzz_0,   \
                             ta1_x_zz_xxzz_1,   \
                             ta1_x_zz_xyyz_0,   \
                             ta1_x_zz_xyyz_1,   \
                             ta1_x_zz_xyzz_0,   \
                             ta1_x_zz_xyzz_1,   \
                             ta1_x_zz_xzzz_0,   \
                             ta1_x_zz_xzzz_1,   \
                             ta1_x_zz_yyyz_0,   \
                             ta1_x_zz_yyyz_1,   \
                             ta1_x_zz_yyzz_0,   \
                             ta1_x_zz_yyzz_1,   \
                             ta1_x_zz_yzzz_0,   \
                             ta1_x_zz_yzzz_1,   \
                             ta1_x_zz_zzzz_0,   \
                             ta1_x_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzz_xxxx_0[i] = ta1_x_zz_xxxx_0[i] * fe_0 - ta1_x_zz_xxxx_1[i] * fe_0 + ta1_x_yzz_xxxx_0[i] * pa_y[i] - ta1_x_yzz_xxxx_1[i] * pc_y[i];

        ta1_x_yyzz_xxxy_0[i] = ta1_x_yy_xxxy_0[i] * fe_0 - ta1_x_yy_xxxy_1[i] * fe_0 + ta1_x_yyz_xxxy_0[i] * pa_z[i] - ta1_x_yyz_xxxy_1[i] * pc_z[i];

        ta1_x_yyzz_xxxz_0[i] = ta1_x_zz_xxxz_0[i] * fe_0 - ta1_x_zz_xxxz_1[i] * fe_0 + ta1_x_yzz_xxxz_0[i] * pa_y[i] - ta1_x_yzz_xxxz_1[i] * pc_y[i];

        ta1_x_yyzz_xxyy_0[i] = ta1_x_yy_xxyy_0[i] * fe_0 - ta1_x_yy_xxyy_1[i] * fe_0 + ta1_x_yyz_xxyy_0[i] * pa_z[i] - ta1_x_yyz_xxyy_1[i] * pc_z[i];

        ta1_x_yyzz_xxyz_0[i] = ta1_x_zz_xxyz_0[i] * fe_0 - ta1_x_zz_xxyz_1[i] * fe_0 + ta1_x_yzz_xxz_0[i] * fe_0 - ta1_x_yzz_xxz_1[i] * fe_0 +
                               ta1_x_yzz_xxyz_0[i] * pa_y[i] - ta1_x_yzz_xxyz_1[i] * pc_y[i];

        ta1_x_yyzz_xxzz_0[i] = ta1_x_zz_xxzz_0[i] * fe_0 - ta1_x_zz_xxzz_1[i] * fe_0 + ta1_x_yzz_xxzz_0[i] * pa_y[i] - ta1_x_yzz_xxzz_1[i] * pc_y[i];

        ta1_x_yyzz_xyyy_0[i] = ta1_x_yy_xyyy_0[i] * fe_0 - ta1_x_yy_xyyy_1[i] * fe_0 + ta1_x_yyz_xyyy_0[i] * pa_z[i] - ta1_x_yyz_xyyy_1[i] * pc_z[i];

        ta1_x_yyzz_xyyz_0[i] = ta1_x_zz_xyyz_0[i] * fe_0 - ta1_x_zz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_x_yzz_xyz_1[i] * fe_0 + ta1_x_yzz_xyyz_0[i] * pa_y[i] - ta1_x_yzz_xyyz_1[i] * pc_y[i];

        ta1_x_yyzz_xyzz_0[i] = ta1_x_zz_xyzz_0[i] * fe_0 - ta1_x_zz_xyzz_1[i] * fe_0 + ta1_x_yzz_xzz_0[i] * fe_0 - ta1_x_yzz_xzz_1[i] * fe_0 +
                               ta1_x_yzz_xyzz_0[i] * pa_y[i] - ta1_x_yzz_xyzz_1[i] * pc_y[i];

        ta1_x_yyzz_xzzz_0[i] = ta1_x_zz_xzzz_0[i] * fe_0 - ta1_x_zz_xzzz_1[i] * fe_0 + ta1_x_yzz_xzzz_0[i] * pa_y[i] - ta1_x_yzz_xzzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyyy_0[i] = ta1_x_yy_yyyy_0[i] * fe_0 - ta1_x_yy_yyyy_1[i] * fe_0 + ta1_x_yyz_yyyy_0[i] * pa_z[i] - ta1_x_yyz_yyyy_1[i] * pc_z[i];

        ta1_x_yyzz_yyyz_0[i] = ta1_x_zz_yyyz_0[i] * fe_0 - ta1_x_zz_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yzz_yyz_0[i] * fe_0 -
                               3.0 * ta1_x_yzz_yyz_1[i] * fe_0 + ta1_x_yzz_yyyz_0[i] * pa_y[i] - ta1_x_yzz_yyyz_1[i] * pc_y[i];

        ta1_x_yyzz_yyzz_0[i] = ta1_x_zz_yyzz_0[i] * fe_0 - ta1_x_zz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yzz_yzz_0[i] * fe_0 -
                               2.0 * ta1_x_yzz_yzz_1[i] * fe_0 + ta1_x_yzz_yyzz_0[i] * pa_y[i] - ta1_x_yzz_yyzz_1[i] * pc_y[i];

        ta1_x_yyzz_yzzz_0[i] = ta1_x_zz_yzzz_0[i] * fe_0 - ta1_x_zz_yzzz_1[i] * fe_0 + ta1_x_yzz_zzz_0[i] * fe_0 - ta1_x_yzz_zzz_1[i] * fe_0 +
                               ta1_x_yzz_yzzz_0[i] * pa_y[i] - ta1_x_yzz_yzzz_1[i] * pc_y[i];

        ta1_x_yyzz_zzzz_0[i] = ta1_x_zz_zzzz_0[i] * fe_0 - ta1_x_zz_zzzz_1[i] * fe_0 + ta1_x_yzz_zzzz_0[i] * pa_y[i] - ta1_x_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 195-210 components of targeted buffer : GG

    auto ta1_x_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 195);

    auto ta1_x_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 196);

    auto ta1_x_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 197);

    auto ta1_x_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 198);

    auto ta1_x_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 199);

    auto ta1_x_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 200);

    auto ta1_x_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 201);

    auto ta1_x_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 202);

    auto ta1_x_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 203);

    auto ta1_x_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 204);

    auto ta1_x_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 205);

    auto ta1_x_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 206);

    auto ta1_x_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 207);

    auto ta1_x_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 208);

    auto ta1_x_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 209);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_yzzz_xxxx_0, \
                             ta1_x_yzzz_xxxy_0, \
                             ta1_x_yzzz_xxxz_0, \
                             ta1_x_yzzz_xxyy_0, \
                             ta1_x_yzzz_xxyz_0, \
                             ta1_x_yzzz_xxzz_0, \
                             ta1_x_yzzz_xyyy_0, \
                             ta1_x_yzzz_xyyz_0, \
                             ta1_x_yzzz_xyzz_0, \
                             ta1_x_yzzz_xzzz_0, \
                             ta1_x_yzzz_yyyy_0, \
                             ta1_x_yzzz_yyyz_0, \
                             ta1_x_yzzz_yyzz_0, \
                             ta1_x_yzzz_yzzz_0, \
                             ta1_x_yzzz_zzzz_0, \
                             ta1_x_zzz_xxx_0,   \
                             ta1_x_zzz_xxx_1,   \
                             ta1_x_zzz_xxxx_0,  \
                             ta1_x_zzz_xxxx_1,  \
                             ta1_x_zzz_xxxy_0,  \
                             ta1_x_zzz_xxxy_1,  \
                             ta1_x_zzz_xxxz_0,  \
                             ta1_x_zzz_xxxz_1,  \
                             ta1_x_zzz_xxy_0,   \
                             ta1_x_zzz_xxy_1,   \
                             ta1_x_zzz_xxyy_0,  \
                             ta1_x_zzz_xxyy_1,  \
                             ta1_x_zzz_xxyz_0,  \
                             ta1_x_zzz_xxyz_1,  \
                             ta1_x_zzz_xxz_0,   \
                             ta1_x_zzz_xxz_1,   \
                             ta1_x_zzz_xxzz_0,  \
                             ta1_x_zzz_xxzz_1,  \
                             ta1_x_zzz_xyy_0,   \
                             ta1_x_zzz_xyy_1,   \
                             ta1_x_zzz_xyyy_0,  \
                             ta1_x_zzz_xyyy_1,  \
                             ta1_x_zzz_xyyz_0,  \
                             ta1_x_zzz_xyyz_1,  \
                             ta1_x_zzz_xyz_0,   \
                             ta1_x_zzz_xyz_1,   \
                             ta1_x_zzz_xyzz_0,  \
                             ta1_x_zzz_xyzz_1,  \
                             ta1_x_zzz_xzz_0,   \
                             ta1_x_zzz_xzz_1,   \
                             ta1_x_zzz_xzzz_0,  \
                             ta1_x_zzz_xzzz_1,  \
                             ta1_x_zzz_yyy_0,   \
                             ta1_x_zzz_yyy_1,   \
                             ta1_x_zzz_yyyy_0,  \
                             ta1_x_zzz_yyyy_1,  \
                             ta1_x_zzz_yyyz_0,  \
                             ta1_x_zzz_yyyz_1,  \
                             ta1_x_zzz_yyz_0,   \
                             ta1_x_zzz_yyz_1,   \
                             ta1_x_zzz_yyzz_0,  \
                             ta1_x_zzz_yyzz_1,  \
                             ta1_x_zzz_yzz_0,   \
                             ta1_x_zzz_yzz_1,   \
                             ta1_x_zzz_yzzz_0,  \
                             ta1_x_zzz_yzzz_1,  \
                             ta1_x_zzz_zzz_0,   \
                             ta1_x_zzz_zzz_1,   \
                             ta1_x_zzz_zzzz_0,  \
                             ta1_x_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzz_xxxx_0[i] = ta1_x_zzz_xxxx_0[i] * pa_y[i] - ta1_x_zzz_xxxx_1[i] * pc_y[i];

        ta1_x_yzzz_xxxy_0[i] = ta1_x_zzz_xxx_0[i] * fe_0 - ta1_x_zzz_xxx_1[i] * fe_0 + ta1_x_zzz_xxxy_0[i] * pa_y[i] - ta1_x_zzz_xxxy_1[i] * pc_y[i];

        ta1_x_yzzz_xxxz_0[i] = ta1_x_zzz_xxxz_0[i] * pa_y[i] - ta1_x_zzz_xxxz_1[i] * pc_y[i];

        ta1_x_yzzz_xxyy_0[i] =
            2.0 * ta1_x_zzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_zzz_xxy_1[i] * fe_0 + ta1_x_zzz_xxyy_0[i] * pa_y[i] - ta1_x_zzz_xxyy_1[i] * pc_y[i];

        ta1_x_yzzz_xxyz_0[i] = ta1_x_zzz_xxz_0[i] * fe_0 - ta1_x_zzz_xxz_1[i] * fe_0 + ta1_x_zzz_xxyz_0[i] * pa_y[i] - ta1_x_zzz_xxyz_1[i] * pc_y[i];

        ta1_x_yzzz_xxzz_0[i] = ta1_x_zzz_xxzz_0[i] * pa_y[i] - ta1_x_zzz_xxzz_1[i] * pc_y[i];

        ta1_x_yzzz_xyyy_0[i] =
            3.0 * ta1_x_zzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_zzz_xyy_1[i] * fe_0 + ta1_x_zzz_xyyy_0[i] * pa_y[i] - ta1_x_zzz_xyyy_1[i] * pc_y[i];

        ta1_x_yzzz_xyyz_0[i] =
            2.0 * ta1_x_zzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyz_1[i] * fe_0 + ta1_x_zzz_xyyz_0[i] * pa_y[i] - ta1_x_zzz_xyyz_1[i] * pc_y[i];

        ta1_x_yzzz_xyzz_0[i] = ta1_x_zzz_xzz_0[i] * fe_0 - ta1_x_zzz_xzz_1[i] * fe_0 + ta1_x_zzz_xyzz_0[i] * pa_y[i] - ta1_x_zzz_xyzz_1[i] * pc_y[i];

        ta1_x_yzzz_xzzz_0[i] = ta1_x_zzz_xzzz_0[i] * pa_y[i] - ta1_x_zzz_xzzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyyy_0[i] =
            4.0 * ta1_x_zzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyy_1[i] * fe_0 + ta1_x_zzz_yyyy_0[i] * pa_y[i] - ta1_x_zzz_yyyy_1[i] * pc_y[i];

        ta1_x_yzzz_yyyz_0[i] =
            3.0 * ta1_x_zzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_zzz_yyz_1[i] * fe_0 + ta1_x_zzz_yyyz_0[i] * pa_y[i] - ta1_x_zzz_yyyz_1[i] * pc_y[i];

        ta1_x_yzzz_yyzz_0[i] =
            2.0 * ta1_x_zzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_yzz_1[i] * fe_0 + ta1_x_zzz_yyzz_0[i] * pa_y[i] - ta1_x_zzz_yyzz_1[i] * pc_y[i];

        ta1_x_yzzz_yzzz_0[i] = ta1_x_zzz_zzz_0[i] * fe_0 - ta1_x_zzz_zzz_1[i] * fe_0 + ta1_x_zzz_yzzz_0[i] * pa_y[i] - ta1_x_zzz_yzzz_1[i] * pc_y[i];

        ta1_x_yzzz_zzzz_0[i] = ta1_x_zzz_zzzz_0[i] * pa_y[i] - ta1_x_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_x_zz_xxxx_0,   \
                             ta1_x_zz_xxxx_1,   \
                             ta1_x_zz_xxxy_0,   \
                             ta1_x_zz_xxxy_1,   \
                             ta1_x_zz_xxxz_0,   \
                             ta1_x_zz_xxxz_1,   \
                             ta1_x_zz_xxyy_0,   \
                             ta1_x_zz_xxyy_1,   \
                             ta1_x_zz_xxyz_0,   \
                             ta1_x_zz_xxyz_1,   \
                             ta1_x_zz_xxzz_0,   \
                             ta1_x_zz_xxzz_1,   \
                             ta1_x_zz_xyyy_0,   \
                             ta1_x_zz_xyyy_1,   \
                             ta1_x_zz_xyyz_0,   \
                             ta1_x_zz_xyyz_1,   \
                             ta1_x_zz_xyzz_0,   \
                             ta1_x_zz_xyzz_1,   \
                             ta1_x_zz_xzzz_0,   \
                             ta1_x_zz_xzzz_1,   \
                             ta1_x_zz_yyyy_0,   \
                             ta1_x_zz_yyyy_1,   \
                             ta1_x_zz_yyyz_0,   \
                             ta1_x_zz_yyyz_1,   \
                             ta1_x_zz_yyzz_0,   \
                             ta1_x_zz_yyzz_1,   \
                             ta1_x_zz_yzzz_0,   \
                             ta1_x_zz_yzzz_1,   \
                             ta1_x_zz_zzzz_0,   \
                             ta1_x_zz_zzzz_1,   \
                             ta1_x_zzz_xxx_0,   \
                             ta1_x_zzz_xxx_1,   \
                             ta1_x_zzz_xxxx_0,  \
                             ta1_x_zzz_xxxx_1,  \
                             ta1_x_zzz_xxxy_0,  \
                             ta1_x_zzz_xxxy_1,  \
                             ta1_x_zzz_xxxz_0,  \
                             ta1_x_zzz_xxxz_1,  \
                             ta1_x_zzz_xxy_0,   \
                             ta1_x_zzz_xxy_1,   \
                             ta1_x_zzz_xxyy_0,  \
                             ta1_x_zzz_xxyy_1,  \
                             ta1_x_zzz_xxyz_0,  \
                             ta1_x_zzz_xxyz_1,  \
                             ta1_x_zzz_xxz_0,   \
                             ta1_x_zzz_xxz_1,   \
                             ta1_x_zzz_xxzz_0,  \
                             ta1_x_zzz_xxzz_1,  \
                             ta1_x_zzz_xyy_0,   \
                             ta1_x_zzz_xyy_1,   \
                             ta1_x_zzz_xyyy_0,  \
                             ta1_x_zzz_xyyy_1,  \
                             ta1_x_zzz_xyyz_0,  \
                             ta1_x_zzz_xyyz_1,  \
                             ta1_x_zzz_xyz_0,   \
                             ta1_x_zzz_xyz_1,   \
                             ta1_x_zzz_xyzz_0,  \
                             ta1_x_zzz_xyzz_1,  \
                             ta1_x_zzz_xzz_0,   \
                             ta1_x_zzz_xzz_1,   \
                             ta1_x_zzz_xzzz_0,  \
                             ta1_x_zzz_xzzz_1,  \
                             ta1_x_zzz_yyy_0,   \
                             ta1_x_zzz_yyy_1,   \
                             ta1_x_zzz_yyyy_0,  \
                             ta1_x_zzz_yyyy_1,  \
                             ta1_x_zzz_yyyz_0,  \
                             ta1_x_zzz_yyyz_1,  \
                             ta1_x_zzz_yyz_0,   \
                             ta1_x_zzz_yyz_1,   \
                             ta1_x_zzz_yyzz_0,  \
                             ta1_x_zzz_yyzz_1,  \
                             ta1_x_zzz_yzz_0,   \
                             ta1_x_zzz_yzz_1,   \
                             ta1_x_zzz_yzzz_0,  \
                             ta1_x_zzz_yzzz_1,  \
                             ta1_x_zzz_zzz_0,   \
                             ta1_x_zzz_zzz_1,   \
                             ta1_x_zzz_zzzz_0,  \
                             ta1_x_zzz_zzzz_1,  \
                             ta1_x_zzzz_xxxx_0, \
                             ta1_x_zzzz_xxxy_0, \
                             ta1_x_zzzz_xxxz_0, \
                             ta1_x_zzzz_xxyy_0, \
                             ta1_x_zzzz_xxyz_0, \
                             ta1_x_zzzz_xxzz_0, \
                             ta1_x_zzzz_xyyy_0, \
                             ta1_x_zzzz_xyyz_0, \
                             ta1_x_zzzz_xyzz_0, \
                             ta1_x_zzzz_xzzz_0, \
                             ta1_x_zzzz_yyyy_0, \
                             ta1_x_zzzz_yyyz_0, \
                             ta1_x_zzzz_yyzz_0, \
                             ta1_x_zzzz_yzzz_0, \
                             ta1_x_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzz_xxxx_0[i] =
            3.0 * ta1_x_zz_xxxx_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxx_1[i] * fe_0 + ta1_x_zzz_xxxx_0[i] * pa_z[i] - ta1_x_zzz_xxxx_1[i] * pc_z[i];

        ta1_x_zzzz_xxxy_0[i] =
            3.0 * ta1_x_zz_xxxy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxy_1[i] * fe_0 + ta1_x_zzz_xxxy_0[i] * pa_z[i] - ta1_x_zzz_xxxy_1[i] * pc_z[i];

        ta1_x_zzzz_xxxz_0[i] = 3.0 * ta1_x_zz_xxxz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxz_1[i] * fe_0 + ta1_x_zzz_xxx_0[i] * fe_0 -
                               ta1_x_zzz_xxx_1[i] * fe_0 + ta1_x_zzz_xxxz_0[i] * pa_z[i] - ta1_x_zzz_xxxz_1[i] * pc_z[i];

        ta1_x_zzzz_xxyy_0[i] =
            3.0 * ta1_x_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyy_1[i] * fe_0 + ta1_x_zzz_xxyy_0[i] * pa_z[i] - ta1_x_zzz_xxyy_1[i] * pc_z[i];

        ta1_x_zzzz_xxyz_0[i] = 3.0 * ta1_x_zz_xxyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyz_1[i] * fe_0 + ta1_x_zzz_xxy_0[i] * fe_0 -
                               ta1_x_zzz_xxy_1[i] * fe_0 + ta1_x_zzz_xxyz_0[i] * pa_z[i] - ta1_x_zzz_xxyz_1[i] * pc_z[i];

        ta1_x_zzzz_xxzz_0[i] = 3.0 * ta1_x_zz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxz_0[i] * fe_0 -
                               2.0 * ta1_x_zzz_xxz_1[i] * fe_0 + ta1_x_zzz_xxzz_0[i] * pa_z[i] - ta1_x_zzz_xxzz_1[i] * pc_z[i];

        ta1_x_zzzz_xyyy_0[i] =
            3.0 * ta1_x_zz_xyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyy_1[i] * fe_0 + ta1_x_zzz_xyyy_0[i] * pa_z[i] - ta1_x_zzz_xyyy_1[i] * pc_z[i];

        ta1_x_zzzz_xyyz_0[i] = 3.0 * ta1_x_zz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyz_1[i] * fe_0 + ta1_x_zzz_xyy_0[i] * fe_0 -
                               ta1_x_zzz_xyy_1[i] * fe_0 + ta1_x_zzz_xyyz_0[i] * pa_z[i] - ta1_x_zzz_xyyz_1[i] * pc_z[i];

        ta1_x_zzzz_xyzz_0[i] = 3.0 * ta1_x_zz_xyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_x_zzz_xyz_1[i] * fe_0 + ta1_x_zzz_xyzz_0[i] * pa_z[i] - ta1_x_zzz_xyzz_1[i] * pc_z[i];

        ta1_x_zzzz_xzzz_0[i] = 3.0 * ta1_x_zz_xzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_xzz_0[i] * fe_0 -
                               3.0 * ta1_x_zzz_xzz_1[i] * fe_0 + ta1_x_zzz_xzzz_0[i] * pa_z[i] - ta1_x_zzz_xzzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyyy_0[i] =
            3.0 * ta1_x_zz_yyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyy_1[i] * fe_0 + ta1_x_zzz_yyyy_0[i] * pa_z[i] - ta1_x_zzz_yyyy_1[i] * pc_z[i];

        ta1_x_zzzz_yyyz_0[i] = 3.0 * ta1_x_zz_yyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyz_1[i] * fe_0 + ta1_x_zzz_yyy_0[i] * fe_0 -
                               ta1_x_zzz_yyy_1[i] * fe_0 + ta1_x_zzz_yyyz_0[i] * pa_z[i] - ta1_x_zzz_yyyz_1[i] * pc_z[i];

        ta1_x_zzzz_yyzz_0[i] = 3.0 * ta1_x_zz_yyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_yyz_0[i] * fe_0 -
                               2.0 * ta1_x_zzz_yyz_1[i] * fe_0 + ta1_x_zzz_yyzz_0[i] * pa_z[i] - ta1_x_zzz_yyzz_1[i] * pc_z[i];

        ta1_x_zzzz_yzzz_0[i] = 3.0 * ta1_x_zz_yzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_yzz_0[i] * fe_0 -
                               3.0 * ta1_x_zzz_yzz_1[i] * fe_0 + ta1_x_zzz_yzzz_0[i] * pa_z[i] - ta1_x_zzz_yzzz_1[i] * pc_z[i];

        ta1_x_zzzz_zzzz_0[i] = 3.0 * ta1_x_zz_zzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_zzzz_1[i] * fe_0 + 4.0 * ta1_x_zzz_zzz_0[i] * fe_0 -
                               4.0 * ta1_x_zzz_zzz_1[i] * fe_0 + ta1_x_zzz_zzzz_0[i] * pa_z[i] - ta1_x_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 225-240 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xx_xxxx_0,   \
                             ta1_y_xx_xxxx_1,   \
                             ta1_y_xx_xxxy_0,   \
                             ta1_y_xx_xxxy_1,   \
                             ta1_y_xx_xxxz_0,   \
                             ta1_y_xx_xxxz_1,   \
                             ta1_y_xx_xxyy_0,   \
                             ta1_y_xx_xxyy_1,   \
                             ta1_y_xx_xxyz_0,   \
                             ta1_y_xx_xxyz_1,   \
                             ta1_y_xx_xxzz_0,   \
                             ta1_y_xx_xxzz_1,   \
                             ta1_y_xx_xyyy_0,   \
                             ta1_y_xx_xyyy_1,   \
                             ta1_y_xx_xyyz_0,   \
                             ta1_y_xx_xyyz_1,   \
                             ta1_y_xx_xyzz_0,   \
                             ta1_y_xx_xyzz_1,   \
                             ta1_y_xx_xzzz_0,   \
                             ta1_y_xx_xzzz_1,   \
                             ta1_y_xx_yyyy_0,   \
                             ta1_y_xx_yyyy_1,   \
                             ta1_y_xx_yyyz_0,   \
                             ta1_y_xx_yyyz_1,   \
                             ta1_y_xx_yyzz_0,   \
                             ta1_y_xx_yyzz_1,   \
                             ta1_y_xx_yzzz_0,   \
                             ta1_y_xx_yzzz_1,   \
                             ta1_y_xx_zzzz_0,   \
                             ta1_y_xx_zzzz_1,   \
                             ta1_y_xxx_xxx_0,   \
                             ta1_y_xxx_xxx_1,   \
                             ta1_y_xxx_xxxx_0,  \
                             ta1_y_xxx_xxxx_1,  \
                             ta1_y_xxx_xxxy_0,  \
                             ta1_y_xxx_xxxy_1,  \
                             ta1_y_xxx_xxxz_0,  \
                             ta1_y_xxx_xxxz_1,  \
                             ta1_y_xxx_xxy_0,   \
                             ta1_y_xxx_xxy_1,   \
                             ta1_y_xxx_xxyy_0,  \
                             ta1_y_xxx_xxyy_1,  \
                             ta1_y_xxx_xxyz_0,  \
                             ta1_y_xxx_xxyz_1,  \
                             ta1_y_xxx_xxz_0,   \
                             ta1_y_xxx_xxz_1,   \
                             ta1_y_xxx_xxzz_0,  \
                             ta1_y_xxx_xxzz_1,  \
                             ta1_y_xxx_xyy_0,   \
                             ta1_y_xxx_xyy_1,   \
                             ta1_y_xxx_xyyy_0,  \
                             ta1_y_xxx_xyyy_1,  \
                             ta1_y_xxx_xyyz_0,  \
                             ta1_y_xxx_xyyz_1,  \
                             ta1_y_xxx_xyz_0,   \
                             ta1_y_xxx_xyz_1,   \
                             ta1_y_xxx_xyzz_0,  \
                             ta1_y_xxx_xyzz_1,  \
                             ta1_y_xxx_xzz_0,   \
                             ta1_y_xxx_xzz_1,   \
                             ta1_y_xxx_xzzz_0,  \
                             ta1_y_xxx_xzzz_1,  \
                             ta1_y_xxx_yyy_0,   \
                             ta1_y_xxx_yyy_1,   \
                             ta1_y_xxx_yyyy_0,  \
                             ta1_y_xxx_yyyy_1,  \
                             ta1_y_xxx_yyyz_0,  \
                             ta1_y_xxx_yyyz_1,  \
                             ta1_y_xxx_yyz_0,   \
                             ta1_y_xxx_yyz_1,   \
                             ta1_y_xxx_yyzz_0,  \
                             ta1_y_xxx_yyzz_1,  \
                             ta1_y_xxx_yzz_0,   \
                             ta1_y_xxx_yzz_1,   \
                             ta1_y_xxx_yzzz_0,  \
                             ta1_y_xxx_yzzz_1,  \
                             ta1_y_xxx_zzz_0,   \
                             ta1_y_xxx_zzz_1,   \
                             ta1_y_xxx_zzzz_0,  \
                             ta1_y_xxx_zzzz_1,  \
                             ta1_y_xxxx_xxxx_0, \
                             ta1_y_xxxx_xxxy_0, \
                             ta1_y_xxxx_xxxz_0, \
                             ta1_y_xxxx_xxyy_0, \
                             ta1_y_xxxx_xxyz_0, \
                             ta1_y_xxxx_xxzz_0, \
                             ta1_y_xxxx_xyyy_0, \
                             ta1_y_xxxx_xyyz_0, \
                             ta1_y_xxxx_xyzz_0, \
                             ta1_y_xxxx_xzzz_0, \
                             ta1_y_xxxx_yyyy_0, \
                             ta1_y_xxxx_yyyz_0, \
                             ta1_y_xxxx_yyzz_0, \
                             ta1_y_xxxx_yzzz_0, \
                             ta1_y_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxx_xxxx_0[i] = 3.0 * ta1_y_xx_xxxx_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxx_1[i] * fe_0 + 4.0 * ta1_y_xxx_xxx_0[i] * fe_0 -
                               4.0 * ta1_y_xxx_xxx_1[i] * fe_0 + ta1_y_xxx_xxxx_0[i] * pa_x[i] - ta1_y_xxx_xxxx_1[i] * pc_x[i];

        ta1_y_xxxx_xxxy_0[i] = 3.0 * ta1_y_xx_xxxy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxy_0[i] * fe_0 -
                               3.0 * ta1_y_xxx_xxy_1[i] * fe_0 + ta1_y_xxx_xxxy_0[i] * pa_x[i] - ta1_y_xxx_xxxy_1[i] * pc_x[i];

        ta1_y_xxxx_xxxz_0[i] = 3.0 * ta1_y_xx_xxxz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxz_0[i] * fe_0 -
                               3.0 * ta1_y_xxx_xxz_1[i] * fe_0 + ta1_y_xxx_xxxz_0[i] * pa_x[i] - ta1_y_xxx_xxxz_1[i] * pc_x[i];

        ta1_y_xxxx_xxyy_0[i] = 3.0 * ta1_y_xx_xxyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyy_0[i] * fe_0 -
                               2.0 * ta1_y_xxx_xyy_1[i] * fe_0 + ta1_y_xxx_xxyy_0[i] * pa_x[i] - ta1_y_xxx_xxyy_1[i] * pc_x[i];

        ta1_y_xxxx_xxyz_0[i] = 3.0 * ta1_y_xx_xxyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyz_0[i] * fe_0 -
                               2.0 * ta1_y_xxx_xyz_1[i] * fe_0 + ta1_y_xxx_xxyz_0[i] * pa_x[i] - ta1_y_xxx_xxyz_1[i] * pc_x[i];

        ta1_y_xxxx_xxzz_0[i] = 3.0 * ta1_y_xx_xxzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xzz_0[i] * fe_0 -
                               2.0 * ta1_y_xxx_xzz_1[i] * fe_0 + ta1_y_xxx_xxzz_0[i] * pa_x[i] - ta1_y_xxx_xxzz_1[i] * pc_x[i];

        ta1_y_xxxx_xyyy_0[i] = 3.0 * ta1_y_xx_xyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyy_1[i] * fe_0 + ta1_y_xxx_yyy_0[i] * fe_0 -
                               ta1_y_xxx_yyy_1[i] * fe_0 + ta1_y_xxx_xyyy_0[i] * pa_x[i] - ta1_y_xxx_xyyy_1[i] * pc_x[i];

        ta1_y_xxxx_xyyz_0[i] = 3.0 * ta1_y_xx_xyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyz_1[i] * fe_0 + ta1_y_xxx_yyz_0[i] * fe_0 -
                               ta1_y_xxx_yyz_1[i] * fe_0 + ta1_y_xxx_xyyz_0[i] * pa_x[i] - ta1_y_xxx_xyyz_1[i] * pc_x[i];

        ta1_y_xxxx_xyzz_0[i] = 3.0 * ta1_y_xx_xyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyzz_1[i] * fe_0 + ta1_y_xxx_yzz_0[i] * fe_0 -
                               ta1_y_xxx_yzz_1[i] * fe_0 + ta1_y_xxx_xyzz_0[i] * pa_x[i] - ta1_y_xxx_xyzz_1[i] * pc_x[i];

        ta1_y_xxxx_xzzz_0[i] = 3.0 * ta1_y_xx_xzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xzzz_1[i] * fe_0 + ta1_y_xxx_zzz_0[i] * fe_0 -
                               ta1_y_xxx_zzz_1[i] * fe_0 + ta1_y_xxx_xzzz_0[i] * pa_x[i] - ta1_y_xxx_xzzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyyy_0[i] =
            3.0 * ta1_y_xx_yyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyy_1[i] * fe_0 + ta1_y_xxx_yyyy_0[i] * pa_x[i] - ta1_y_xxx_yyyy_1[i] * pc_x[i];

        ta1_y_xxxx_yyyz_0[i] =
            3.0 * ta1_y_xx_yyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyz_1[i] * fe_0 + ta1_y_xxx_yyyz_0[i] * pa_x[i] - ta1_y_xxx_yyyz_1[i] * pc_x[i];

        ta1_y_xxxx_yyzz_0[i] =
            3.0 * ta1_y_xx_yyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyzz_1[i] * fe_0 + ta1_y_xxx_yyzz_0[i] * pa_x[i] - ta1_y_xxx_yyzz_1[i] * pc_x[i];

        ta1_y_xxxx_yzzz_0[i] =
            3.0 * ta1_y_xx_yzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yzzz_1[i] * fe_0 + ta1_y_xxx_yzzz_0[i] * pa_x[i] - ta1_y_xxx_yzzz_1[i] * pc_x[i];

        ta1_y_xxxx_zzzz_0[i] =
            3.0 * ta1_y_xx_zzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_zzzz_1[i] * fe_0 + ta1_y_xxx_zzzz_0[i] * pa_x[i] - ta1_y_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : GG

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

    auto ta1_y_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 254);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xxx_xxx_0,   \
                             ta1_y_xxx_xxx_1,   \
                             ta1_y_xxx_xxxx_0,  \
                             ta1_y_xxx_xxxx_1,  \
                             ta1_y_xxx_xxxy_0,  \
                             ta1_y_xxx_xxxy_1,  \
                             ta1_y_xxx_xxxz_0,  \
                             ta1_y_xxx_xxxz_1,  \
                             ta1_y_xxx_xxy_0,   \
                             ta1_y_xxx_xxy_1,   \
                             ta1_y_xxx_xxyy_0,  \
                             ta1_y_xxx_xxyy_1,  \
                             ta1_y_xxx_xxyz_0,  \
                             ta1_y_xxx_xxyz_1,  \
                             ta1_y_xxx_xxz_0,   \
                             ta1_y_xxx_xxz_1,   \
                             ta1_y_xxx_xxzz_0,  \
                             ta1_y_xxx_xxzz_1,  \
                             ta1_y_xxx_xyy_0,   \
                             ta1_y_xxx_xyy_1,   \
                             ta1_y_xxx_xyyy_0,  \
                             ta1_y_xxx_xyyy_1,  \
                             ta1_y_xxx_xyyz_0,  \
                             ta1_y_xxx_xyyz_1,  \
                             ta1_y_xxx_xyz_0,   \
                             ta1_y_xxx_xyz_1,   \
                             ta1_y_xxx_xyzz_0,  \
                             ta1_y_xxx_xyzz_1,  \
                             ta1_y_xxx_xzz_0,   \
                             ta1_y_xxx_xzz_1,   \
                             ta1_y_xxx_xzzz_0,  \
                             ta1_y_xxx_xzzz_1,  \
                             ta1_y_xxx_zzzz_0,  \
                             ta1_y_xxx_zzzz_1,  \
                             ta1_y_xxxy_xxxx_0, \
                             ta1_y_xxxy_xxxy_0, \
                             ta1_y_xxxy_xxxz_0, \
                             ta1_y_xxxy_xxyy_0, \
                             ta1_y_xxxy_xxyz_0, \
                             ta1_y_xxxy_xxzz_0, \
                             ta1_y_xxxy_xyyy_0, \
                             ta1_y_xxxy_xyyz_0, \
                             ta1_y_xxxy_xyzz_0, \
                             ta1_y_xxxy_xzzz_0, \
                             ta1_y_xxxy_yyyy_0, \
                             ta1_y_xxxy_yyyz_0, \
                             ta1_y_xxxy_yyzz_0, \
                             ta1_y_xxxy_yzzz_0, \
                             ta1_y_xxxy_zzzz_0, \
                             ta1_y_xxy_yyyy_0,  \
                             ta1_y_xxy_yyyy_1,  \
                             ta1_y_xxy_yyyz_0,  \
                             ta1_y_xxy_yyyz_1,  \
                             ta1_y_xxy_yyzz_0,  \
                             ta1_y_xxy_yyzz_1,  \
                             ta1_y_xxy_yzzz_0,  \
                             ta1_y_xxy_yzzz_1,  \
                             ta1_y_xy_yyyy_0,   \
                             ta1_y_xy_yyyy_1,   \
                             ta1_y_xy_yyyz_0,   \
                             ta1_y_xy_yyyz_1,   \
                             ta1_y_xy_yyzz_0,   \
                             ta1_y_xy_yyzz_1,   \
                             ta1_y_xy_yzzz_0,   \
                             ta1_y_xy_yzzz_1,   \
                             ta_xxx_xxxx_1,     \
                             ta_xxx_xxxy_1,     \
                             ta_xxx_xxxz_1,     \
                             ta_xxx_xxyy_1,     \
                             ta_xxx_xxyz_1,     \
                             ta_xxx_xxzz_1,     \
                             ta_xxx_xyyy_1,     \
                             ta_xxx_xyyz_1,     \
                             ta_xxx_xyzz_1,     \
                             ta_xxx_xzzz_1,     \
                             ta_xxx_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxy_xxxx_0[i] = ta_xxx_xxxx_1[i] + ta1_y_xxx_xxxx_0[i] * pa_y[i] - ta1_y_xxx_xxxx_1[i] * pc_y[i];

        ta1_y_xxxy_xxxy_0[i] =
            ta1_y_xxx_xxx_0[i] * fe_0 - ta1_y_xxx_xxx_1[i] * fe_0 + ta_xxx_xxxy_1[i] + ta1_y_xxx_xxxy_0[i] * pa_y[i] - ta1_y_xxx_xxxy_1[i] * pc_y[i];

        ta1_y_xxxy_xxxz_0[i] = ta_xxx_xxxz_1[i] + ta1_y_xxx_xxxz_0[i] * pa_y[i] - ta1_y_xxx_xxxz_1[i] * pc_y[i];

        ta1_y_xxxy_xxyy_0[i] = 2.0 * ta1_y_xxx_xxy_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxy_1[i] * fe_0 + ta_xxx_xxyy_1[i] + ta1_y_xxx_xxyy_0[i] * pa_y[i] -
                               ta1_y_xxx_xxyy_1[i] * pc_y[i];

        ta1_y_xxxy_xxyz_0[i] =
            ta1_y_xxx_xxz_0[i] * fe_0 - ta1_y_xxx_xxz_1[i] * fe_0 + ta_xxx_xxyz_1[i] + ta1_y_xxx_xxyz_0[i] * pa_y[i] - ta1_y_xxx_xxyz_1[i] * pc_y[i];

        ta1_y_xxxy_xxzz_0[i] = ta_xxx_xxzz_1[i] + ta1_y_xxx_xxzz_0[i] * pa_y[i] - ta1_y_xxx_xxzz_1[i] * pc_y[i];

        ta1_y_xxxy_xyyy_0[i] = 3.0 * ta1_y_xxx_xyy_0[i] * fe_0 - 3.0 * ta1_y_xxx_xyy_1[i] * fe_0 + ta_xxx_xyyy_1[i] + ta1_y_xxx_xyyy_0[i] * pa_y[i] -
                               ta1_y_xxx_xyyy_1[i] * pc_y[i];

        ta1_y_xxxy_xyyz_0[i] = 2.0 * ta1_y_xxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xyz_1[i] * fe_0 + ta_xxx_xyyz_1[i] + ta1_y_xxx_xyyz_0[i] * pa_y[i] -
                               ta1_y_xxx_xyyz_1[i] * pc_y[i];

        ta1_y_xxxy_xyzz_0[i] =
            ta1_y_xxx_xzz_0[i] * fe_0 - ta1_y_xxx_xzz_1[i] * fe_0 + ta_xxx_xyzz_1[i] + ta1_y_xxx_xyzz_0[i] * pa_y[i] - ta1_y_xxx_xyzz_1[i] * pc_y[i];

        ta1_y_xxxy_xzzz_0[i] = ta_xxx_xzzz_1[i] + ta1_y_xxx_xzzz_0[i] * pa_y[i] - ta1_y_xxx_xzzz_1[i] * pc_y[i];

        ta1_y_xxxy_yyyy_0[i] =
            2.0 * ta1_y_xy_yyyy_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyy_1[i] * fe_0 + ta1_y_xxy_yyyy_0[i] * pa_x[i] - ta1_y_xxy_yyyy_1[i] * pc_x[i];

        ta1_y_xxxy_yyyz_0[i] =
            2.0 * ta1_y_xy_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyz_1[i] * fe_0 + ta1_y_xxy_yyyz_0[i] * pa_x[i] - ta1_y_xxy_yyyz_1[i] * pc_x[i];

        ta1_y_xxxy_yyzz_0[i] =
            2.0 * ta1_y_xy_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyzz_1[i] * fe_0 + ta1_y_xxy_yyzz_0[i] * pa_x[i] - ta1_y_xxy_yyzz_1[i] * pc_x[i];

        ta1_y_xxxy_yzzz_0[i] =
            2.0 * ta1_y_xy_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yzzz_1[i] * fe_0 + ta1_y_xxy_yzzz_0[i] * pa_x[i] - ta1_y_xxy_yzzz_1[i] * pc_x[i];

        ta1_y_xxxy_zzzz_0[i] = ta_xxx_zzzz_1[i] + ta1_y_xxx_zzzz_0[i] * pa_y[i] - ta1_y_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : GG

    auto ta1_y_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 255);

    auto ta1_y_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 256);

    auto ta1_y_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 257);

    auto ta1_y_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 258);

    auto ta1_y_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 259);

    auto ta1_y_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 260);

    auto ta1_y_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 261);

    auto ta1_y_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 262);

    auto ta1_y_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 263);

    auto ta1_y_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 264);

    auto ta1_y_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 265);

    auto ta1_y_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 266);

    auto ta1_y_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 267);

    auto ta1_y_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 268);

    auto ta1_y_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 269);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxx_xxx_0,   \
                             ta1_y_xxx_xxx_1,   \
                             ta1_y_xxx_xxxx_0,  \
                             ta1_y_xxx_xxxx_1,  \
                             ta1_y_xxx_xxxy_0,  \
                             ta1_y_xxx_xxxy_1,  \
                             ta1_y_xxx_xxxz_0,  \
                             ta1_y_xxx_xxxz_1,  \
                             ta1_y_xxx_xxy_0,   \
                             ta1_y_xxx_xxy_1,   \
                             ta1_y_xxx_xxyy_0,  \
                             ta1_y_xxx_xxyy_1,  \
                             ta1_y_xxx_xxyz_0,  \
                             ta1_y_xxx_xxyz_1,  \
                             ta1_y_xxx_xxz_0,   \
                             ta1_y_xxx_xxz_1,   \
                             ta1_y_xxx_xxzz_0,  \
                             ta1_y_xxx_xxzz_1,  \
                             ta1_y_xxx_xyy_0,   \
                             ta1_y_xxx_xyy_1,   \
                             ta1_y_xxx_xyyy_0,  \
                             ta1_y_xxx_xyyy_1,  \
                             ta1_y_xxx_xyyz_0,  \
                             ta1_y_xxx_xyyz_1,  \
                             ta1_y_xxx_xyz_0,   \
                             ta1_y_xxx_xyz_1,   \
                             ta1_y_xxx_xyzz_0,  \
                             ta1_y_xxx_xyzz_1,  \
                             ta1_y_xxx_xzz_0,   \
                             ta1_y_xxx_xzz_1,   \
                             ta1_y_xxx_xzzz_0,  \
                             ta1_y_xxx_xzzz_1,  \
                             ta1_y_xxx_yyyy_0,  \
                             ta1_y_xxx_yyyy_1,  \
                             ta1_y_xxxz_xxxx_0, \
                             ta1_y_xxxz_xxxy_0, \
                             ta1_y_xxxz_xxxz_0, \
                             ta1_y_xxxz_xxyy_0, \
                             ta1_y_xxxz_xxyz_0, \
                             ta1_y_xxxz_xxzz_0, \
                             ta1_y_xxxz_xyyy_0, \
                             ta1_y_xxxz_xyyz_0, \
                             ta1_y_xxxz_xyzz_0, \
                             ta1_y_xxxz_xzzz_0, \
                             ta1_y_xxxz_yyyy_0, \
                             ta1_y_xxxz_yyyz_0, \
                             ta1_y_xxxz_yyzz_0, \
                             ta1_y_xxxz_yzzz_0, \
                             ta1_y_xxxz_zzzz_0, \
                             ta1_y_xxz_yyyz_0,  \
                             ta1_y_xxz_yyyz_1,  \
                             ta1_y_xxz_yyzz_0,  \
                             ta1_y_xxz_yyzz_1,  \
                             ta1_y_xxz_yzzz_0,  \
                             ta1_y_xxz_yzzz_1,  \
                             ta1_y_xxz_zzzz_0,  \
                             ta1_y_xxz_zzzz_1,  \
                             ta1_y_xz_yyyz_0,   \
                             ta1_y_xz_yyyz_1,   \
                             ta1_y_xz_yyzz_0,   \
                             ta1_y_xz_yyzz_1,   \
                             ta1_y_xz_yzzz_0,   \
                             ta1_y_xz_yzzz_1,   \
                             ta1_y_xz_zzzz_0,   \
                             ta1_y_xz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxz_xxxx_0[i] = ta1_y_xxx_xxxx_0[i] * pa_z[i] - ta1_y_xxx_xxxx_1[i] * pc_z[i];

        ta1_y_xxxz_xxxy_0[i] = ta1_y_xxx_xxxy_0[i] * pa_z[i] - ta1_y_xxx_xxxy_1[i] * pc_z[i];

        ta1_y_xxxz_xxxz_0[i] = ta1_y_xxx_xxx_0[i] * fe_0 - ta1_y_xxx_xxx_1[i] * fe_0 + ta1_y_xxx_xxxz_0[i] * pa_z[i] - ta1_y_xxx_xxxz_1[i] * pc_z[i];

        ta1_y_xxxz_xxyy_0[i] = ta1_y_xxx_xxyy_0[i] * pa_z[i] - ta1_y_xxx_xxyy_1[i] * pc_z[i];

        ta1_y_xxxz_xxyz_0[i] = ta1_y_xxx_xxy_0[i] * fe_0 - ta1_y_xxx_xxy_1[i] * fe_0 + ta1_y_xxx_xxyz_0[i] * pa_z[i] - ta1_y_xxx_xxyz_1[i] * pc_z[i];

        ta1_y_xxxz_xxzz_0[i] =
            2.0 * ta1_y_xxx_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxz_1[i] * fe_0 + ta1_y_xxx_xxzz_0[i] * pa_z[i] - ta1_y_xxx_xxzz_1[i] * pc_z[i];

        ta1_y_xxxz_xyyy_0[i] = ta1_y_xxx_xyyy_0[i] * pa_z[i] - ta1_y_xxx_xyyy_1[i] * pc_z[i];

        ta1_y_xxxz_xyyz_0[i] = ta1_y_xxx_xyy_0[i] * fe_0 - ta1_y_xxx_xyy_1[i] * fe_0 + ta1_y_xxx_xyyz_0[i] * pa_z[i] - ta1_y_xxx_xyyz_1[i] * pc_z[i];

        ta1_y_xxxz_xyzz_0[i] =
            2.0 * ta1_y_xxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xyz_1[i] * fe_0 + ta1_y_xxx_xyzz_0[i] * pa_z[i] - ta1_y_xxx_xyzz_1[i] * pc_z[i];

        ta1_y_xxxz_xzzz_0[i] =
            3.0 * ta1_y_xxx_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xzz_1[i] * fe_0 + ta1_y_xxx_xzzz_0[i] * pa_z[i] - ta1_y_xxx_xzzz_1[i] * pc_z[i];

        ta1_y_xxxz_yyyy_0[i] = ta1_y_xxx_yyyy_0[i] * pa_z[i] - ta1_y_xxx_yyyy_1[i] * pc_z[i];

        ta1_y_xxxz_yyyz_0[i] =
            2.0 * ta1_y_xz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyyz_1[i] * fe_0 + ta1_y_xxz_yyyz_0[i] * pa_x[i] - ta1_y_xxz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxz_yyzz_0[i] =
            2.0 * ta1_y_xz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyzz_1[i] * fe_0 + ta1_y_xxz_yyzz_0[i] * pa_x[i] - ta1_y_xxz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxz_yzzz_0[i] =
            2.0 * ta1_y_xz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yzzz_1[i] * fe_0 + ta1_y_xxz_yzzz_0[i] * pa_x[i] - ta1_y_xxz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxz_zzzz_0[i] =
            2.0 * ta1_y_xz_zzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_zzzz_1[i] * fe_0 + ta1_y_xxz_zzzz_0[i] * pa_x[i] - ta1_y_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 270-285 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xx_xxxx_0,   \
                             ta1_y_xx_xxxx_1,   \
                             ta1_y_xx_xxxz_0,   \
                             ta1_y_xx_xxxz_1,   \
                             ta1_y_xx_xxzz_0,   \
                             ta1_y_xx_xxzz_1,   \
                             ta1_y_xx_xzzz_0,   \
                             ta1_y_xx_xzzz_1,   \
                             ta1_y_xxy_xxxx_0,  \
                             ta1_y_xxy_xxxx_1,  \
                             ta1_y_xxy_xxxz_0,  \
                             ta1_y_xxy_xxxz_1,  \
                             ta1_y_xxy_xxzz_0,  \
                             ta1_y_xxy_xxzz_1,  \
                             ta1_y_xxy_xzzz_0,  \
                             ta1_y_xxy_xzzz_1,  \
                             ta1_y_xxyy_xxxx_0, \
                             ta1_y_xxyy_xxxy_0, \
                             ta1_y_xxyy_xxxz_0, \
                             ta1_y_xxyy_xxyy_0, \
                             ta1_y_xxyy_xxyz_0, \
                             ta1_y_xxyy_xxzz_0, \
                             ta1_y_xxyy_xyyy_0, \
                             ta1_y_xxyy_xyyz_0, \
                             ta1_y_xxyy_xyzz_0, \
                             ta1_y_xxyy_xzzz_0, \
                             ta1_y_xxyy_yyyy_0, \
                             ta1_y_xxyy_yyyz_0, \
                             ta1_y_xxyy_yyzz_0, \
                             ta1_y_xxyy_yzzz_0, \
                             ta1_y_xxyy_zzzz_0, \
                             ta1_y_xyy_xxxy_0,  \
                             ta1_y_xyy_xxxy_1,  \
                             ta1_y_xyy_xxy_0,   \
                             ta1_y_xyy_xxy_1,   \
                             ta1_y_xyy_xxyy_0,  \
                             ta1_y_xyy_xxyy_1,  \
                             ta1_y_xyy_xxyz_0,  \
                             ta1_y_xyy_xxyz_1,  \
                             ta1_y_xyy_xyy_0,   \
                             ta1_y_xyy_xyy_1,   \
                             ta1_y_xyy_xyyy_0,  \
                             ta1_y_xyy_xyyy_1,  \
                             ta1_y_xyy_xyyz_0,  \
                             ta1_y_xyy_xyyz_1,  \
                             ta1_y_xyy_xyz_0,   \
                             ta1_y_xyy_xyz_1,   \
                             ta1_y_xyy_xyzz_0,  \
                             ta1_y_xyy_xyzz_1,  \
                             ta1_y_xyy_yyy_0,   \
                             ta1_y_xyy_yyy_1,   \
                             ta1_y_xyy_yyyy_0,  \
                             ta1_y_xyy_yyyy_1,  \
                             ta1_y_xyy_yyyz_0,  \
                             ta1_y_xyy_yyyz_1,  \
                             ta1_y_xyy_yyz_0,   \
                             ta1_y_xyy_yyz_1,   \
                             ta1_y_xyy_yyzz_0,  \
                             ta1_y_xyy_yyzz_1,  \
                             ta1_y_xyy_yzz_0,   \
                             ta1_y_xyy_yzz_1,   \
                             ta1_y_xyy_yzzz_0,  \
                             ta1_y_xyy_yzzz_1,  \
                             ta1_y_xyy_zzzz_0,  \
                             ta1_y_xyy_zzzz_1,  \
                             ta1_y_yy_xxxy_0,   \
                             ta1_y_yy_xxxy_1,   \
                             ta1_y_yy_xxyy_0,   \
                             ta1_y_yy_xxyy_1,   \
                             ta1_y_yy_xxyz_0,   \
                             ta1_y_yy_xxyz_1,   \
                             ta1_y_yy_xyyy_0,   \
                             ta1_y_yy_xyyy_1,   \
                             ta1_y_yy_xyyz_0,   \
                             ta1_y_yy_xyyz_1,   \
                             ta1_y_yy_xyzz_0,   \
                             ta1_y_yy_xyzz_1,   \
                             ta1_y_yy_yyyy_0,   \
                             ta1_y_yy_yyyy_1,   \
                             ta1_y_yy_yyyz_0,   \
                             ta1_y_yy_yyyz_1,   \
                             ta1_y_yy_yyzz_0,   \
                             ta1_y_yy_yyzz_1,   \
                             ta1_y_yy_yzzz_0,   \
                             ta1_y_yy_yzzz_1,   \
                             ta1_y_yy_zzzz_0,   \
                             ta1_y_yy_zzzz_1,   \
                             ta_xxy_xxxx_1,     \
                             ta_xxy_xxxz_1,     \
                             ta_xxy_xxzz_1,     \
                             ta_xxy_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyy_xxxx_0[i] =
            ta1_y_xx_xxxx_0[i] * fe_0 - ta1_y_xx_xxxx_1[i] * fe_0 + ta_xxy_xxxx_1[i] + ta1_y_xxy_xxxx_0[i] * pa_y[i] - ta1_y_xxy_xxxx_1[i] * pc_y[i];

        ta1_y_xxyy_xxxy_0[i] = ta1_y_yy_xxxy_0[i] * fe_0 - ta1_y_yy_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xyy_xxy_0[i] * fe_0 -
                               3.0 * ta1_y_xyy_xxy_1[i] * fe_0 + ta1_y_xyy_xxxy_0[i] * pa_x[i] - ta1_y_xyy_xxxy_1[i] * pc_x[i];

        ta1_y_xxyy_xxxz_0[i] =
            ta1_y_xx_xxxz_0[i] * fe_0 - ta1_y_xx_xxxz_1[i] * fe_0 + ta_xxy_xxxz_1[i] + ta1_y_xxy_xxxz_0[i] * pa_y[i] - ta1_y_xxy_xxxz_1[i] * pc_y[i];

        ta1_y_xxyy_xxyy_0[i] = ta1_y_yy_xxyy_0[i] * fe_0 - ta1_y_yy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyy_0[i] * fe_0 -
                               2.0 * ta1_y_xyy_xyy_1[i] * fe_0 + ta1_y_xyy_xxyy_0[i] * pa_x[i] - ta1_y_xyy_xxyy_1[i] * pc_x[i];

        ta1_y_xxyy_xxyz_0[i] = ta1_y_yy_xxyz_0[i] * fe_0 - ta1_y_yy_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyz_0[i] * fe_0 -
                               2.0 * ta1_y_xyy_xyz_1[i] * fe_0 + ta1_y_xyy_xxyz_0[i] * pa_x[i] - ta1_y_xyy_xxyz_1[i] * pc_x[i];

        ta1_y_xxyy_xxzz_0[i] =
            ta1_y_xx_xxzz_0[i] * fe_0 - ta1_y_xx_xxzz_1[i] * fe_0 + ta_xxy_xxzz_1[i] + ta1_y_xxy_xxzz_0[i] * pa_y[i] - ta1_y_xxy_xxzz_1[i] * pc_y[i];

        ta1_y_xxyy_xyyy_0[i] = ta1_y_yy_xyyy_0[i] * fe_0 - ta1_y_yy_xyyy_1[i] * fe_0 + ta1_y_xyy_yyy_0[i] * fe_0 - ta1_y_xyy_yyy_1[i] * fe_0 +
                               ta1_y_xyy_xyyy_0[i] * pa_x[i] - ta1_y_xyy_xyyy_1[i] * pc_x[i];

        ta1_y_xxyy_xyyz_0[i] = ta1_y_yy_xyyz_0[i] * fe_0 - ta1_y_yy_xyyz_1[i] * fe_0 + ta1_y_xyy_yyz_0[i] * fe_0 - ta1_y_xyy_yyz_1[i] * fe_0 +
                               ta1_y_xyy_xyyz_0[i] * pa_x[i] - ta1_y_xyy_xyyz_1[i] * pc_x[i];

        ta1_y_xxyy_xyzz_0[i] = ta1_y_yy_xyzz_0[i] * fe_0 - ta1_y_yy_xyzz_1[i] * fe_0 + ta1_y_xyy_yzz_0[i] * fe_0 - ta1_y_xyy_yzz_1[i] * fe_0 +
                               ta1_y_xyy_xyzz_0[i] * pa_x[i] - ta1_y_xyy_xyzz_1[i] * pc_x[i];

        ta1_y_xxyy_xzzz_0[i] =
            ta1_y_xx_xzzz_0[i] * fe_0 - ta1_y_xx_xzzz_1[i] * fe_0 + ta_xxy_xzzz_1[i] + ta1_y_xxy_xzzz_0[i] * pa_y[i] - ta1_y_xxy_xzzz_1[i] * pc_y[i];

        ta1_y_xxyy_yyyy_0[i] = ta1_y_yy_yyyy_0[i] * fe_0 - ta1_y_yy_yyyy_1[i] * fe_0 + ta1_y_xyy_yyyy_0[i] * pa_x[i] - ta1_y_xyy_yyyy_1[i] * pc_x[i];

        ta1_y_xxyy_yyyz_0[i] = ta1_y_yy_yyyz_0[i] * fe_0 - ta1_y_yy_yyyz_1[i] * fe_0 + ta1_y_xyy_yyyz_0[i] * pa_x[i] - ta1_y_xyy_yyyz_1[i] * pc_x[i];

        ta1_y_xxyy_yyzz_0[i] = ta1_y_yy_yyzz_0[i] * fe_0 - ta1_y_yy_yyzz_1[i] * fe_0 + ta1_y_xyy_yyzz_0[i] * pa_x[i] - ta1_y_xyy_yyzz_1[i] * pc_x[i];

        ta1_y_xxyy_yzzz_0[i] = ta1_y_yy_yzzz_0[i] * fe_0 - ta1_y_yy_yzzz_1[i] * fe_0 + ta1_y_xyy_yzzz_0[i] * pa_x[i] - ta1_y_xyy_yzzz_1[i] * pc_x[i];

        ta1_y_xxyy_zzzz_0[i] = ta1_y_yy_zzzz_0[i] * fe_0 - ta1_y_yy_zzzz_1[i] * fe_0 + ta1_y_xyy_zzzz_0[i] * pa_x[i] - ta1_y_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 285-300 components of targeted buffer : GG

    auto ta1_y_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 285);

    auto ta1_y_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 286);

    auto ta1_y_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 287);

    auto ta1_y_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 288);

    auto ta1_y_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 289);

    auto ta1_y_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 290);

    auto ta1_y_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 291);

    auto ta1_y_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 292);

    auto ta1_y_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 293);

    auto ta1_y_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 294);

    auto ta1_y_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 295);

    auto ta1_y_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 296);

    auto ta1_y_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 297);

    auto ta1_y_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 298);

    auto ta1_y_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 299);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_xxy_xxxx_0,  \
                             ta1_y_xxy_xxxx_1,  \
                             ta1_y_xxy_xxxy_0,  \
                             ta1_y_xxy_xxxy_1,  \
                             ta1_y_xxy_xxy_0,   \
                             ta1_y_xxy_xxy_1,   \
                             ta1_y_xxy_xxyy_0,  \
                             ta1_y_xxy_xxyy_1,  \
                             ta1_y_xxy_xxyz_0,  \
                             ta1_y_xxy_xxyz_1,  \
                             ta1_y_xxy_xyy_0,   \
                             ta1_y_xxy_xyy_1,   \
                             ta1_y_xxy_xyyy_0,  \
                             ta1_y_xxy_xyyy_1,  \
                             ta1_y_xxy_xyyz_0,  \
                             ta1_y_xxy_xyyz_1,  \
                             ta1_y_xxy_xyz_0,   \
                             ta1_y_xxy_xyz_1,   \
                             ta1_y_xxy_xyzz_0,  \
                             ta1_y_xxy_xyzz_1,  \
                             ta1_y_xxy_yyyy_0,  \
                             ta1_y_xxy_yyyy_1,  \
                             ta1_y_xxyz_xxxx_0, \
                             ta1_y_xxyz_xxxy_0, \
                             ta1_y_xxyz_xxxz_0, \
                             ta1_y_xxyz_xxyy_0, \
                             ta1_y_xxyz_xxyz_0, \
                             ta1_y_xxyz_xxzz_0, \
                             ta1_y_xxyz_xyyy_0, \
                             ta1_y_xxyz_xyyz_0, \
                             ta1_y_xxyz_xyzz_0, \
                             ta1_y_xxyz_xzzz_0, \
                             ta1_y_xxyz_yyyy_0, \
                             ta1_y_xxyz_yyyz_0, \
                             ta1_y_xxyz_yyzz_0, \
                             ta1_y_xxyz_yzzz_0, \
                             ta1_y_xxyz_zzzz_0, \
                             ta1_y_xxz_xxxz_0,  \
                             ta1_y_xxz_xxxz_1,  \
                             ta1_y_xxz_xxzz_0,  \
                             ta1_y_xxz_xxzz_1,  \
                             ta1_y_xxz_xzzz_0,  \
                             ta1_y_xxz_xzzz_1,  \
                             ta1_y_xxz_zzzz_0,  \
                             ta1_y_xxz_zzzz_1,  \
                             ta1_y_xyz_yyyz_0,  \
                             ta1_y_xyz_yyyz_1,  \
                             ta1_y_xyz_yyzz_0,  \
                             ta1_y_xyz_yyzz_1,  \
                             ta1_y_xyz_yzzz_0,  \
                             ta1_y_xyz_yzzz_1,  \
                             ta1_y_yz_yyyz_0,   \
                             ta1_y_yz_yyyz_1,   \
                             ta1_y_yz_yyzz_0,   \
                             ta1_y_yz_yyzz_1,   \
                             ta1_y_yz_yzzz_0,   \
                             ta1_y_yz_yzzz_1,   \
                             ta_xxz_xxxz_1,     \
                             ta_xxz_xxzz_1,     \
                             ta_xxz_xzzz_1,     \
                             ta_xxz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyz_xxxx_0[i] = ta1_y_xxy_xxxx_0[i] * pa_z[i] - ta1_y_xxy_xxxx_1[i] * pc_z[i];

        ta1_y_xxyz_xxxy_0[i] = ta1_y_xxy_xxxy_0[i] * pa_z[i] - ta1_y_xxy_xxxy_1[i] * pc_z[i];

        ta1_y_xxyz_xxxz_0[i] = ta_xxz_xxxz_1[i] + ta1_y_xxz_xxxz_0[i] * pa_y[i] - ta1_y_xxz_xxxz_1[i] * pc_y[i];

        ta1_y_xxyz_xxyy_0[i] = ta1_y_xxy_xxyy_0[i] * pa_z[i] - ta1_y_xxy_xxyy_1[i] * pc_z[i];

        ta1_y_xxyz_xxyz_0[i] = ta1_y_xxy_xxy_0[i] * fe_0 - ta1_y_xxy_xxy_1[i] * fe_0 + ta1_y_xxy_xxyz_0[i] * pa_z[i] - ta1_y_xxy_xxyz_1[i] * pc_z[i];

        ta1_y_xxyz_xxzz_0[i] = ta_xxz_xxzz_1[i] + ta1_y_xxz_xxzz_0[i] * pa_y[i] - ta1_y_xxz_xxzz_1[i] * pc_y[i];

        ta1_y_xxyz_xyyy_0[i] = ta1_y_xxy_xyyy_0[i] * pa_z[i] - ta1_y_xxy_xyyy_1[i] * pc_z[i];

        ta1_y_xxyz_xyyz_0[i] = ta1_y_xxy_xyy_0[i] * fe_0 - ta1_y_xxy_xyy_1[i] * fe_0 + ta1_y_xxy_xyyz_0[i] * pa_z[i] - ta1_y_xxy_xyyz_1[i] * pc_z[i];

        ta1_y_xxyz_xyzz_0[i] =
            2.0 * ta1_y_xxy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xyz_1[i] * fe_0 + ta1_y_xxy_xyzz_0[i] * pa_z[i] - ta1_y_xxy_xyzz_1[i] * pc_z[i];

        ta1_y_xxyz_xzzz_0[i] = ta_xxz_xzzz_1[i] + ta1_y_xxz_xzzz_0[i] * pa_y[i] - ta1_y_xxz_xzzz_1[i] * pc_y[i];

        ta1_y_xxyz_yyyy_0[i] = ta1_y_xxy_yyyy_0[i] * pa_z[i] - ta1_y_xxy_yyyy_1[i] * pc_z[i];

        ta1_y_xxyz_yyyz_0[i] = ta1_y_yz_yyyz_0[i] * fe_0 - ta1_y_yz_yyyz_1[i] * fe_0 + ta1_y_xyz_yyyz_0[i] * pa_x[i] - ta1_y_xyz_yyyz_1[i] * pc_x[i];

        ta1_y_xxyz_yyzz_0[i] = ta1_y_yz_yyzz_0[i] * fe_0 - ta1_y_yz_yyzz_1[i] * fe_0 + ta1_y_xyz_yyzz_0[i] * pa_x[i] - ta1_y_xyz_yyzz_1[i] * pc_x[i];

        ta1_y_xxyz_yzzz_0[i] = ta1_y_yz_yzzz_0[i] * fe_0 - ta1_y_yz_yzzz_1[i] * fe_0 + ta1_y_xyz_yzzz_0[i] * pa_x[i] - ta1_y_xyz_yzzz_1[i] * pc_x[i];

        ta1_y_xxyz_zzzz_0[i] = ta_xxz_zzzz_1[i] + ta1_y_xxz_zzzz_0[i] * pa_y[i] - ta1_y_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 300-315 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xx_xxxx_0,   \
                             ta1_y_xx_xxxx_1,   \
                             ta1_y_xx_xxxy_0,   \
                             ta1_y_xx_xxxy_1,   \
                             ta1_y_xx_xxyy_0,   \
                             ta1_y_xx_xxyy_1,   \
                             ta1_y_xx_xyyy_0,   \
                             ta1_y_xx_xyyy_1,   \
                             ta1_y_xxz_xxxx_0,  \
                             ta1_y_xxz_xxxx_1,  \
                             ta1_y_xxz_xxxy_0,  \
                             ta1_y_xxz_xxxy_1,  \
                             ta1_y_xxz_xxyy_0,  \
                             ta1_y_xxz_xxyy_1,  \
                             ta1_y_xxz_xyyy_0,  \
                             ta1_y_xxz_xyyy_1,  \
                             ta1_y_xxzz_xxxx_0, \
                             ta1_y_xxzz_xxxy_0, \
                             ta1_y_xxzz_xxxz_0, \
                             ta1_y_xxzz_xxyy_0, \
                             ta1_y_xxzz_xxyz_0, \
                             ta1_y_xxzz_xxzz_0, \
                             ta1_y_xxzz_xyyy_0, \
                             ta1_y_xxzz_xyyz_0, \
                             ta1_y_xxzz_xyzz_0, \
                             ta1_y_xxzz_xzzz_0, \
                             ta1_y_xxzz_yyyy_0, \
                             ta1_y_xxzz_yyyz_0, \
                             ta1_y_xxzz_yyzz_0, \
                             ta1_y_xxzz_yzzz_0, \
                             ta1_y_xxzz_zzzz_0, \
                             ta1_y_xzz_xxxz_0,  \
                             ta1_y_xzz_xxxz_1,  \
                             ta1_y_xzz_xxyz_0,  \
                             ta1_y_xzz_xxyz_1,  \
                             ta1_y_xzz_xxz_0,   \
                             ta1_y_xzz_xxz_1,   \
                             ta1_y_xzz_xxzz_0,  \
                             ta1_y_xzz_xxzz_1,  \
                             ta1_y_xzz_xyyz_0,  \
                             ta1_y_xzz_xyyz_1,  \
                             ta1_y_xzz_xyz_0,   \
                             ta1_y_xzz_xyz_1,   \
                             ta1_y_xzz_xyzz_0,  \
                             ta1_y_xzz_xyzz_1,  \
                             ta1_y_xzz_xzz_0,   \
                             ta1_y_xzz_xzz_1,   \
                             ta1_y_xzz_xzzz_0,  \
                             ta1_y_xzz_xzzz_1,  \
                             ta1_y_xzz_yyyy_0,  \
                             ta1_y_xzz_yyyy_1,  \
                             ta1_y_xzz_yyyz_0,  \
                             ta1_y_xzz_yyyz_1,  \
                             ta1_y_xzz_yyz_0,   \
                             ta1_y_xzz_yyz_1,   \
                             ta1_y_xzz_yyzz_0,  \
                             ta1_y_xzz_yyzz_1,  \
                             ta1_y_xzz_yzz_0,   \
                             ta1_y_xzz_yzz_1,   \
                             ta1_y_xzz_yzzz_0,  \
                             ta1_y_xzz_yzzz_1,  \
                             ta1_y_xzz_zzz_0,   \
                             ta1_y_xzz_zzz_1,   \
                             ta1_y_xzz_zzzz_0,  \
                             ta1_y_xzz_zzzz_1,  \
                             ta1_y_zz_xxxz_0,   \
                             ta1_y_zz_xxxz_1,   \
                             ta1_y_zz_xxyz_0,   \
                             ta1_y_zz_xxyz_1,   \
                             ta1_y_zz_xxzz_0,   \
                             ta1_y_zz_xxzz_1,   \
                             ta1_y_zz_xyyz_0,   \
                             ta1_y_zz_xyyz_1,   \
                             ta1_y_zz_xyzz_0,   \
                             ta1_y_zz_xyzz_1,   \
                             ta1_y_zz_xzzz_0,   \
                             ta1_y_zz_xzzz_1,   \
                             ta1_y_zz_yyyy_0,   \
                             ta1_y_zz_yyyy_1,   \
                             ta1_y_zz_yyyz_0,   \
                             ta1_y_zz_yyyz_1,   \
                             ta1_y_zz_yyzz_0,   \
                             ta1_y_zz_yyzz_1,   \
                             ta1_y_zz_yzzz_0,   \
                             ta1_y_zz_yzzz_1,   \
                             ta1_y_zz_zzzz_0,   \
                             ta1_y_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzz_xxxx_0[i] = ta1_y_xx_xxxx_0[i] * fe_0 - ta1_y_xx_xxxx_1[i] * fe_0 + ta1_y_xxz_xxxx_0[i] * pa_z[i] - ta1_y_xxz_xxxx_1[i] * pc_z[i];

        ta1_y_xxzz_xxxy_0[i] = ta1_y_xx_xxxy_0[i] * fe_0 - ta1_y_xx_xxxy_1[i] * fe_0 + ta1_y_xxz_xxxy_0[i] * pa_z[i] - ta1_y_xxz_xxxy_1[i] * pc_z[i];

        ta1_y_xxzz_xxxz_0[i] = ta1_y_zz_xxxz_0[i] * fe_0 - ta1_y_zz_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xzz_xxz_0[i] * fe_0 -
                               3.0 * ta1_y_xzz_xxz_1[i] * fe_0 + ta1_y_xzz_xxxz_0[i] * pa_x[i] - ta1_y_xzz_xxxz_1[i] * pc_x[i];

        ta1_y_xxzz_xxyy_0[i] = ta1_y_xx_xxyy_0[i] * fe_0 - ta1_y_xx_xxyy_1[i] * fe_0 + ta1_y_xxz_xxyy_0[i] * pa_z[i] - ta1_y_xxz_xxyy_1[i] * pc_z[i];

        ta1_y_xxzz_xxyz_0[i] = ta1_y_zz_xxyz_0[i] * fe_0 - ta1_y_zz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_y_xzz_xyz_1[i] * fe_0 + ta1_y_xzz_xxyz_0[i] * pa_x[i] - ta1_y_xzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxzz_xxzz_0[i] = ta1_y_zz_xxzz_0[i] * fe_0 - ta1_y_zz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xzz_0[i] * fe_0 -
                               2.0 * ta1_y_xzz_xzz_1[i] * fe_0 + ta1_y_xzz_xxzz_0[i] * pa_x[i] - ta1_y_xzz_xxzz_1[i] * pc_x[i];

        ta1_y_xxzz_xyyy_0[i] = ta1_y_xx_xyyy_0[i] * fe_0 - ta1_y_xx_xyyy_1[i] * fe_0 + ta1_y_xxz_xyyy_0[i] * pa_z[i] - ta1_y_xxz_xyyy_1[i] * pc_z[i];

        ta1_y_xxzz_xyyz_0[i] = ta1_y_zz_xyyz_0[i] * fe_0 - ta1_y_zz_xyyz_1[i] * fe_0 + ta1_y_xzz_yyz_0[i] * fe_0 - ta1_y_xzz_yyz_1[i] * fe_0 +
                               ta1_y_xzz_xyyz_0[i] * pa_x[i] - ta1_y_xzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxzz_xyzz_0[i] = ta1_y_zz_xyzz_0[i] * fe_0 - ta1_y_zz_xyzz_1[i] * fe_0 + ta1_y_xzz_yzz_0[i] * fe_0 - ta1_y_xzz_yzz_1[i] * fe_0 +
                               ta1_y_xzz_xyzz_0[i] * pa_x[i] - ta1_y_xzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxzz_xzzz_0[i] = ta1_y_zz_xzzz_0[i] * fe_0 - ta1_y_zz_xzzz_1[i] * fe_0 + ta1_y_xzz_zzz_0[i] * fe_0 - ta1_y_xzz_zzz_1[i] * fe_0 +
                               ta1_y_xzz_xzzz_0[i] * pa_x[i] - ta1_y_xzz_xzzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyyy_0[i] = ta1_y_zz_yyyy_0[i] * fe_0 - ta1_y_zz_yyyy_1[i] * fe_0 + ta1_y_xzz_yyyy_0[i] * pa_x[i] - ta1_y_xzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxzz_yyyz_0[i] = ta1_y_zz_yyyz_0[i] * fe_0 - ta1_y_zz_yyyz_1[i] * fe_0 + ta1_y_xzz_yyyz_0[i] * pa_x[i] - ta1_y_xzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxzz_yyzz_0[i] = ta1_y_zz_yyzz_0[i] * fe_0 - ta1_y_zz_yyzz_1[i] * fe_0 + ta1_y_xzz_yyzz_0[i] * pa_x[i] - ta1_y_xzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxzz_yzzz_0[i] = ta1_y_zz_yzzz_0[i] * fe_0 - ta1_y_zz_yzzz_1[i] * fe_0 + ta1_y_xzz_yzzz_0[i] * pa_x[i] - ta1_y_xzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxzz_zzzz_0[i] = ta1_y_zz_zzzz_0[i] * fe_0 - ta1_y_zz_zzzz_1[i] * fe_0 + ta1_y_xzz_zzzz_0[i] * pa_x[i] - ta1_y_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : GG

    auto ta1_y_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 315);

    auto ta1_y_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 316);

    auto ta1_y_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 317);

    auto ta1_y_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 318);

    auto ta1_y_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 319);

    auto ta1_y_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 320);

    auto ta1_y_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 321);

    auto ta1_y_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 322);

    auto ta1_y_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 323);

    auto ta1_y_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 324);

    auto ta1_y_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 325);

    auto ta1_y_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 326);

    auto ta1_y_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 327);

    auto ta1_y_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 328);

    auto ta1_y_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 329);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xyyy_xxxx_0, \
                             ta1_y_xyyy_xxxy_0, \
                             ta1_y_xyyy_xxxz_0, \
                             ta1_y_xyyy_xxyy_0, \
                             ta1_y_xyyy_xxyz_0, \
                             ta1_y_xyyy_xxzz_0, \
                             ta1_y_xyyy_xyyy_0, \
                             ta1_y_xyyy_xyyz_0, \
                             ta1_y_xyyy_xyzz_0, \
                             ta1_y_xyyy_xzzz_0, \
                             ta1_y_xyyy_yyyy_0, \
                             ta1_y_xyyy_yyyz_0, \
                             ta1_y_xyyy_yyzz_0, \
                             ta1_y_xyyy_yzzz_0, \
                             ta1_y_xyyy_zzzz_0, \
                             ta1_y_yyy_xxx_0,   \
                             ta1_y_yyy_xxx_1,   \
                             ta1_y_yyy_xxxx_0,  \
                             ta1_y_yyy_xxxx_1,  \
                             ta1_y_yyy_xxxy_0,  \
                             ta1_y_yyy_xxxy_1,  \
                             ta1_y_yyy_xxxz_0,  \
                             ta1_y_yyy_xxxz_1,  \
                             ta1_y_yyy_xxy_0,   \
                             ta1_y_yyy_xxy_1,   \
                             ta1_y_yyy_xxyy_0,  \
                             ta1_y_yyy_xxyy_1,  \
                             ta1_y_yyy_xxyz_0,  \
                             ta1_y_yyy_xxyz_1,  \
                             ta1_y_yyy_xxz_0,   \
                             ta1_y_yyy_xxz_1,   \
                             ta1_y_yyy_xxzz_0,  \
                             ta1_y_yyy_xxzz_1,  \
                             ta1_y_yyy_xyy_0,   \
                             ta1_y_yyy_xyy_1,   \
                             ta1_y_yyy_xyyy_0,  \
                             ta1_y_yyy_xyyy_1,  \
                             ta1_y_yyy_xyyz_0,  \
                             ta1_y_yyy_xyyz_1,  \
                             ta1_y_yyy_xyz_0,   \
                             ta1_y_yyy_xyz_1,   \
                             ta1_y_yyy_xyzz_0,  \
                             ta1_y_yyy_xyzz_1,  \
                             ta1_y_yyy_xzz_0,   \
                             ta1_y_yyy_xzz_1,   \
                             ta1_y_yyy_xzzz_0,  \
                             ta1_y_yyy_xzzz_1,  \
                             ta1_y_yyy_yyy_0,   \
                             ta1_y_yyy_yyy_1,   \
                             ta1_y_yyy_yyyy_0,  \
                             ta1_y_yyy_yyyy_1,  \
                             ta1_y_yyy_yyyz_0,  \
                             ta1_y_yyy_yyyz_1,  \
                             ta1_y_yyy_yyz_0,   \
                             ta1_y_yyy_yyz_1,   \
                             ta1_y_yyy_yyzz_0,  \
                             ta1_y_yyy_yyzz_1,  \
                             ta1_y_yyy_yzz_0,   \
                             ta1_y_yyy_yzz_1,   \
                             ta1_y_yyy_yzzz_0,  \
                             ta1_y_yyy_yzzz_1,  \
                             ta1_y_yyy_zzz_0,   \
                             ta1_y_yyy_zzz_1,   \
                             ta1_y_yyy_zzzz_0,  \
                             ta1_y_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyy_xxxx_0[i] =
            4.0 * ta1_y_yyy_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxx_1[i] * fe_0 + ta1_y_yyy_xxxx_0[i] * pa_x[i] - ta1_y_yyy_xxxx_1[i] * pc_x[i];

        ta1_y_xyyy_xxxy_0[i] =
            3.0 * ta1_y_yyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxy_1[i] * fe_0 + ta1_y_yyy_xxxy_0[i] * pa_x[i] - ta1_y_yyy_xxxy_1[i] * pc_x[i];

        ta1_y_xyyy_xxxz_0[i] =
            3.0 * ta1_y_yyy_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxz_1[i] * fe_0 + ta1_y_yyy_xxxz_0[i] * pa_x[i] - ta1_y_yyy_xxxz_1[i] * pc_x[i];

        ta1_y_xyyy_xxyy_0[i] =
            2.0 * ta1_y_yyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyy_1[i] * fe_0 + ta1_y_yyy_xxyy_0[i] * pa_x[i] - ta1_y_yyy_xxyy_1[i] * pc_x[i];

        ta1_y_xyyy_xxyz_0[i] =
            2.0 * ta1_y_yyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyz_1[i] * fe_0 + ta1_y_yyy_xxyz_0[i] * pa_x[i] - ta1_y_yyy_xxyz_1[i] * pc_x[i];

        ta1_y_xyyy_xxzz_0[i] =
            2.0 * ta1_y_yyy_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xzz_1[i] * fe_0 + ta1_y_yyy_xxzz_0[i] * pa_x[i] - ta1_y_yyy_xxzz_1[i] * pc_x[i];

        ta1_y_xyyy_xyyy_0[i] = ta1_y_yyy_yyy_0[i] * fe_0 - ta1_y_yyy_yyy_1[i] * fe_0 + ta1_y_yyy_xyyy_0[i] * pa_x[i] - ta1_y_yyy_xyyy_1[i] * pc_x[i];

        ta1_y_xyyy_xyyz_0[i] = ta1_y_yyy_yyz_0[i] * fe_0 - ta1_y_yyy_yyz_1[i] * fe_0 + ta1_y_yyy_xyyz_0[i] * pa_x[i] - ta1_y_yyy_xyyz_1[i] * pc_x[i];

        ta1_y_xyyy_xyzz_0[i] = ta1_y_yyy_yzz_0[i] * fe_0 - ta1_y_yyy_yzz_1[i] * fe_0 + ta1_y_yyy_xyzz_0[i] * pa_x[i] - ta1_y_yyy_xyzz_1[i] * pc_x[i];

        ta1_y_xyyy_xzzz_0[i] = ta1_y_yyy_zzz_0[i] * fe_0 - ta1_y_yyy_zzz_1[i] * fe_0 + ta1_y_yyy_xzzz_0[i] * pa_x[i] - ta1_y_yyy_xzzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyyy_0[i] = ta1_y_yyy_yyyy_0[i] * pa_x[i] - ta1_y_yyy_yyyy_1[i] * pc_x[i];

        ta1_y_xyyy_yyyz_0[i] = ta1_y_yyy_yyyz_0[i] * pa_x[i] - ta1_y_yyy_yyyz_1[i] * pc_x[i];

        ta1_y_xyyy_yyzz_0[i] = ta1_y_yyy_yyzz_0[i] * pa_x[i] - ta1_y_yyy_yyzz_1[i] * pc_x[i];

        ta1_y_xyyy_yzzz_0[i] = ta1_y_yyy_yzzz_0[i] * pa_x[i] - ta1_y_yyy_yzzz_1[i] * pc_x[i];

        ta1_y_xyyy_zzzz_0[i] = ta1_y_yyy_zzzz_0[i] * pa_x[i] - ta1_y_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 330-345 components of targeted buffer : GG

    auto ta1_y_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 330);

    auto ta1_y_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 331);

    auto ta1_y_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 332);

    auto ta1_y_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 333);

    auto ta1_y_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 334);

    auto ta1_y_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 335);

    auto ta1_y_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 336);

    auto ta1_y_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 337);

    auto ta1_y_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 338);

    auto ta1_y_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 339);

    auto ta1_y_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 340);

    auto ta1_y_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 341);

    auto ta1_y_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 342);

    auto ta1_y_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 343);

    auto ta1_y_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 344);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xyy_xxxx_0,  \
                             ta1_y_xyy_xxxx_1,  \
                             ta1_y_xyy_xxxy_0,  \
                             ta1_y_xyy_xxxy_1,  \
                             ta1_y_xyy_xxyy_0,  \
                             ta1_y_xyy_xxyy_1,  \
                             ta1_y_xyy_xyyy_0,  \
                             ta1_y_xyy_xyyy_1,  \
                             ta1_y_xyyz_xxxx_0, \
                             ta1_y_xyyz_xxxy_0, \
                             ta1_y_xyyz_xxxz_0, \
                             ta1_y_xyyz_xxyy_0, \
                             ta1_y_xyyz_xxyz_0, \
                             ta1_y_xyyz_xxzz_0, \
                             ta1_y_xyyz_xyyy_0, \
                             ta1_y_xyyz_xyyz_0, \
                             ta1_y_xyyz_xyzz_0, \
                             ta1_y_xyyz_xzzz_0, \
                             ta1_y_xyyz_yyyy_0, \
                             ta1_y_xyyz_yyyz_0, \
                             ta1_y_xyyz_yyzz_0, \
                             ta1_y_xyyz_yzzz_0, \
                             ta1_y_xyyz_zzzz_0, \
                             ta1_y_yyz_xxxz_0,  \
                             ta1_y_yyz_xxxz_1,  \
                             ta1_y_yyz_xxyz_0,  \
                             ta1_y_yyz_xxyz_1,  \
                             ta1_y_yyz_xxz_0,   \
                             ta1_y_yyz_xxz_1,   \
                             ta1_y_yyz_xxzz_0,  \
                             ta1_y_yyz_xxzz_1,  \
                             ta1_y_yyz_xyyz_0,  \
                             ta1_y_yyz_xyyz_1,  \
                             ta1_y_yyz_xyz_0,   \
                             ta1_y_yyz_xyz_1,   \
                             ta1_y_yyz_xyzz_0,  \
                             ta1_y_yyz_xyzz_1,  \
                             ta1_y_yyz_xzz_0,   \
                             ta1_y_yyz_xzz_1,   \
                             ta1_y_yyz_xzzz_0,  \
                             ta1_y_yyz_xzzz_1,  \
                             ta1_y_yyz_yyyy_0,  \
                             ta1_y_yyz_yyyy_1,  \
                             ta1_y_yyz_yyyz_0,  \
                             ta1_y_yyz_yyyz_1,  \
                             ta1_y_yyz_yyz_0,   \
                             ta1_y_yyz_yyz_1,   \
                             ta1_y_yyz_yyzz_0,  \
                             ta1_y_yyz_yyzz_1,  \
                             ta1_y_yyz_yzz_0,   \
                             ta1_y_yyz_yzz_1,   \
                             ta1_y_yyz_yzzz_0,  \
                             ta1_y_yyz_yzzz_1,  \
                             ta1_y_yyz_zzz_0,   \
                             ta1_y_yyz_zzz_1,   \
                             ta1_y_yyz_zzzz_0,  \
                             ta1_y_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyz_xxxx_0[i] = ta1_y_xyy_xxxx_0[i] * pa_z[i] - ta1_y_xyy_xxxx_1[i] * pc_z[i];

        ta1_y_xyyz_xxxy_0[i] = ta1_y_xyy_xxxy_0[i] * pa_z[i] - ta1_y_xyy_xxxy_1[i] * pc_z[i];

        ta1_y_xyyz_xxxz_0[i] =
            3.0 * ta1_y_yyz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyz_xxz_1[i] * fe_0 + ta1_y_yyz_xxxz_0[i] * pa_x[i] - ta1_y_yyz_xxxz_1[i] * pc_x[i];

        ta1_y_xyyz_xxyy_0[i] = ta1_y_xyy_xxyy_0[i] * pa_z[i] - ta1_y_xyy_xxyy_1[i] * pc_z[i];

        ta1_y_xyyz_xxyz_0[i] =
            2.0 * ta1_y_yyz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyz_1[i] * fe_0 + ta1_y_yyz_xxyz_0[i] * pa_x[i] - ta1_y_yyz_xxyz_1[i] * pc_x[i];

        ta1_y_xyyz_xxzz_0[i] =
            2.0 * ta1_y_yyz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xzz_1[i] * fe_0 + ta1_y_yyz_xxzz_0[i] * pa_x[i] - ta1_y_yyz_xxzz_1[i] * pc_x[i];

        ta1_y_xyyz_xyyy_0[i] = ta1_y_xyy_xyyy_0[i] * pa_z[i] - ta1_y_xyy_xyyy_1[i] * pc_z[i];

        ta1_y_xyyz_xyyz_0[i] = ta1_y_yyz_yyz_0[i] * fe_0 - ta1_y_yyz_yyz_1[i] * fe_0 + ta1_y_yyz_xyyz_0[i] * pa_x[i] - ta1_y_yyz_xyyz_1[i] * pc_x[i];

        ta1_y_xyyz_xyzz_0[i] = ta1_y_yyz_yzz_0[i] * fe_0 - ta1_y_yyz_yzz_1[i] * fe_0 + ta1_y_yyz_xyzz_0[i] * pa_x[i] - ta1_y_yyz_xyzz_1[i] * pc_x[i];

        ta1_y_xyyz_xzzz_0[i] = ta1_y_yyz_zzz_0[i] * fe_0 - ta1_y_yyz_zzz_1[i] * fe_0 + ta1_y_yyz_xzzz_0[i] * pa_x[i] - ta1_y_yyz_xzzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyyy_0[i] = ta1_y_yyz_yyyy_0[i] * pa_x[i] - ta1_y_yyz_yyyy_1[i] * pc_x[i];

        ta1_y_xyyz_yyyz_0[i] = ta1_y_yyz_yyyz_0[i] * pa_x[i] - ta1_y_yyz_yyyz_1[i] * pc_x[i];

        ta1_y_xyyz_yyzz_0[i] = ta1_y_yyz_yyzz_0[i] * pa_x[i] - ta1_y_yyz_yyzz_1[i] * pc_x[i];

        ta1_y_xyyz_yzzz_0[i] = ta1_y_yyz_yzzz_0[i] * pa_x[i] - ta1_y_yyz_yzzz_1[i] * pc_x[i];

        ta1_y_xyyz_zzzz_0[i] = ta1_y_yyz_zzzz_0[i] * pa_x[i] - ta1_y_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 345-360 components of targeted buffer : GG

    auto ta1_y_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 345);

    auto ta1_y_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 346);

    auto ta1_y_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 347);

    auto ta1_y_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 348);

    auto ta1_y_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 349);

    auto ta1_y_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 350);

    auto ta1_y_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 351);

    auto ta1_y_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 352);

    auto ta1_y_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 353);

    auto ta1_y_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 354);

    auto ta1_y_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 355);

    auto ta1_y_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 356);

    auto ta1_y_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 357);

    auto ta1_y_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 358);

    auto ta1_y_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 359);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xyzz_xxxx_0, \
                             ta1_y_xyzz_xxxy_0, \
                             ta1_y_xyzz_xxxz_0, \
                             ta1_y_xyzz_xxyy_0, \
                             ta1_y_xyzz_xxyz_0, \
                             ta1_y_xyzz_xxzz_0, \
                             ta1_y_xyzz_xyyy_0, \
                             ta1_y_xyzz_xyyz_0, \
                             ta1_y_xyzz_xyzz_0, \
                             ta1_y_xyzz_xzzz_0, \
                             ta1_y_xyzz_yyyy_0, \
                             ta1_y_xyzz_yyyz_0, \
                             ta1_y_xyzz_yyzz_0, \
                             ta1_y_xyzz_yzzz_0, \
                             ta1_y_xyzz_zzzz_0, \
                             ta1_y_xzz_xxxx_0,  \
                             ta1_y_xzz_xxxx_1,  \
                             ta1_y_xzz_xxxz_0,  \
                             ta1_y_xzz_xxxz_1,  \
                             ta1_y_xzz_xxzz_0,  \
                             ta1_y_xzz_xxzz_1,  \
                             ta1_y_xzz_xzzz_0,  \
                             ta1_y_xzz_xzzz_1,  \
                             ta1_y_yzz_xxxy_0,  \
                             ta1_y_yzz_xxxy_1,  \
                             ta1_y_yzz_xxy_0,   \
                             ta1_y_yzz_xxy_1,   \
                             ta1_y_yzz_xxyy_0,  \
                             ta1_y_yzz_xxyy_1,  \
                             ta1_y_yzz_xxyz_0,  \
                             ta1_y_yzz_xxyz_1,  \
                             ta1_y_yzz_xyy_0,   \
                             ta1_y_yzz_xyy_1,   \
                             ta1_y_yzz_xyyy_0,  \
                             ta1_y_yzz_xyyy_1,  \
                             ta1_y_yzz_xyyz_0,  \
                             ta1_y_yzz_xyyz_1,  \
                             ta1_y_yzz_xyz_0,   \
                             ta1_y_yzz_xyz_1,   \
                             ta1_y_yzz_xyzz_0,  \
                             ta1_y_yzz_xyzz_1,  \
                             ta1_y_yzz_yyy_0,   \
                             ta1_y_yzz_yyy_1,   \
                             ta1_y_yzz_yyyy_0,  \
                             ta1_y_yzz_yyyy_1,  \
                             ta1_y_yzz_yyyz_0,  \
                             ta1_y_yzz_yyyz_1,  \
                             ta1_y_yzz_yyz_0,   \
                             ta1_y_yzz_yyz_1,   \
                             ta1_y_yzz_yyzz_0,  \
                             ta1_y_yzz_yyzz_1,  \
                             ta1_y_yzz_yzz_0,   \
                             ta1_y_yzz_yzz_1,   \
                             ta1_y_yzz_yzzz_0,  \
                             ta1_y_yzz_yzzz_1,  \
                             ta1_y_yzz_zzzz_0,  \
                             ta1_y_yzz_zzzz_1,  \
                             ta_xzz_xxxx_1,     \
                             ta_xzz_xxxz_1,     \
                             ta_xzz_xxzz_1,     \
                             ta_xzz_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzz_xxxx_0[i] = ta_xzz_xxxx_1[i] + ta1_y_xzz_xxxx_0[i] * pa_y[i] - ta1_y_xzz_xxxx_1[i] * pc_y[i];

        ta1_y_xyzz_xxxy_0[i] =
            3.0 * ta1_y_yzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxy_1[i] * fe_0 + ta1_y_yzz_xxxy_0[i] * pa_x[i] - ta1_y_yzz_xxxy_1[i] * pc_x[i];

        ta1_y_xyzz_xxxz_0[i] = ta_xzz_xxxz_1[i] + ta1_y_xzz_xxxz_0[i] * pa_y[i] - ta1_y_xzz_xxxz_1[i] * pc_y[i];

        ta1_y_xyzz_xxyy_0[i] =
            2.0 * ta1_y_yzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyy_1[i] * fe_0 + ta1_y_yzz_xxyy_0[i] * pa_x[i] - ta1_y_yzz_xxyy_1[i] * pc_x[i];

        ta1_y_xyzz_xxyz_0[i] =
            2.0 * ta1_y_yzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyz_1[i] * fe_0 + ta1_y_yzz_xxyz_0[i] * pa_x[i] - ta1_y_yzz_xxyz_1[i] * pc_x[i];

        ta1_y_xyzz_xxzz_0[i] = ta_xzz_xxzz_1[i] + ta1_y_xzz_xxzz_0[i] * pa_y[i] - ta1_y_xzz_xxzz_1[i] * pc_y[i];

        ta1_y_xyzz_xyyy_0[i] = ta1_y_yzz_yyy_0[i] * fe_0 - ta1_y_yzz_yyy_1[i] * fe_0 + ta1_y_yzz_xyyy_0[i] * pa_x[i] - ta1_y_yzz_xyyy_1[i] * pc_x[i];

        ta1_y_xyzz_xyyz_0[i] = ta1_y_yzz_yyz_0[i] * fe_0 - ta1_y_yzz_yyz_1[i] * fe_0 + ta1_y_yzz_xyyz_0[i] * pa_x[i] - ta1_y_yzz_xyyz_1[i] * pc_x[i];

        ta1_y_xyzz_xyzz_0[i] = ta1_y_yzz_yzz_0[i] * fe_0 - ta1_y_yzz_yzz_1[i] * fe_0 + ta1_y_yzz_xyzz_0[i] * pa_x[i] - ta1_y_yzz_xyzz_1[i] * pc_x[i];

        ta1_y_xyzz_xzzz_0[i] = ta_xzz_xzzz_1[i] + ta1_y_xzz_xzzz_0[i] * pa_y[i] - ta1_y_xzz_xzzz_1[i] * pc_y[i];

        ta1_y_xyzz_yyyy_0[i] = ta1_y_yzz_yyyy_0[i] * pa_x[i] - ta1_y_yzz_yyyy_1[i] * pc_x[i];

        ta1_y_xyzz_yyyz_0[i] = ta1_y_yzz_yyyz_0[i] * pa_x[i] - ta1_y_yzz_yyyz_1[i] * pc_x[i];

        ta1_y_xyzz_yyzz_0[i] = ta1_y_yzz_yyzz_0[i] * pa_x[i] - ta1_y_yzz_yyzz_1[i] * pc_x[i];

        ta1_y_xyzz_yzzz_0[i] = ta1_y_yzz_yzzz_0[i] * pa_x[i] - ta1_y_yzz_yzzz_1[i] * pc_x[i];

        ta1_y_xyzz_zzzz_0[i] = ta1_y_yzz_zzzz_0[i] * pa_x[i] - ta1_y_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 360-375 components of targeted buffer : GG

    auto ta1_y_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 360);

    auto ta1_y_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 361);

    auto ta1_y_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 362);

    auto ta1_y_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 363);

    auto ta1_y_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 364);

    auto ta1_y_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 365);

    auto ta1_y_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 366);

    auto ta1_y_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 367);

    auto ta1_y_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 368);

    auto ta1_y_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 369);

    auto ta1_y_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 370);

    auto ta1_y_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 371);

    auto ta1_y_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 372);

    auto ta1_y_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 373);

    auto ta1_y_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 374);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xzzz_xxxx_0, \
                             ta1_y_xzzz_xxxy_0, \
                             ta1_y_xzzz_xxxz_0, \
                             ta1_y_xzzz_xxyy_0, \
                             ta1_y_xzzz_xxyz_0, \
                             ta1_y_xzzz_xxzz_0, \
                             ta1_y_xzzz_xyyy_0, \
                             ta1_y_xzzz_xyyz_0, \
                             ta1_y_xzzz_xyzz_0, \
                             ta1_y_xzzz_xzzz_0, \
                             ta1_y_xzzz_yyyy_0, \
                             ta1_y_xzzz_yyyz_0, \
                             ta1_y_xzzz_yyzz_0, \
                             ta1_y_xzzz_yzzz_0, \
                             ta1_y_xzzz_zzzz_0, \
                             ta1_y_zzz_xxx_0,   \
                             ta1_y_zzz_xxx_1,   \
                             ta1_y_zzz_xxxx_0,  \
                             ta1_y_zzz_xxxx_1,  \
                             ta1_y_zzz_xxxy_0,  \
                             ta1_y_zzz_xxxy_1,  \
                             ta1_y_zzz_xxxz_0,  \
                             ta1_y_zzz_xxxz_1,  \
                             ta1_y_zzz_xxy_0,   \
                             ta1_y_zzz_xxy_1,   \
                             ta1_y_zzz_xxyy_0,  \
                             ta1_y_zzz_xxyy_1,  \
                             ta1_y_zzz_xxyz_0,  \
                             ta1_y_zzz_xxyz_1,  \
                             ta1_y_zzz_xxz_0,   \
                             ta1_y_zzz_xxz_1,   \
                             ta1_y_zzz_xxzz_0,  \
                             ta1_y_zzz_xxzz_1,  \
                             ta1_y_zzz_xyy_0,   \
                             ta1_y_zzz_xyy_1,   \
                             ta1_y_zzz_xyyy_0,  \
                             ta1_y_zzz_xyyy_1,  \
                             ta1_y_zzz_xyyz_0,  \
                             ta1_y_zzz_xyyz_1,  \
                             ta1_y_zzz_xyz_0,   \
                             ta1_y_zzz_xyz_1,   \
                             ta1_y_zzz_xyzz_0,  \
                             ta1_y_zzz_xyzz_1,  \
                             ta1_y_zzz_xzz_0,   \
                             ta1_y_zzz_xzz_1,   \
                             ta1_y_zzz_xzzz_0,  \
                             ta1_y_zzz_xzzz_1,  \
                             ta1_y_zzz_yyy_0,   \
                             ta1_y_zzz_yyy_1,   \
                             ta1_y_zzz_yyyy_0,  \
                             ta1_y_zzz_yyyy_1,  \
                             ta1_y_zzz_yyyz_0,  \
                             ta1_y_zzz_yyyz_1,  \
                             ta1_y_zzz_yyz_0,   \
                             ta1_y_zzz_yyz_1,   \
                             ta1_y_zzz_yyzz_0,  \
                             ta1_y_zzz_yyzz_1,  \
                             ta1_y_zzz_yzz_0,   \
                             ta1_y_zzz_yzz_1,   \
                             ta1_y_zzz_yzzz_0,  \
                             ta1_y_zzz_yzzz_1,  \
                             ta1_y_zzz_zzz_0,   \
                             ta1_y_zzz_zzz_1,   \
                             ta1_y_zzz_zzzz_0,  \
                             ta1_y_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzz_xxxx_0[i] =
            4.0 * ta1_y_zzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxx_1[i] * fe_0 + ta1_y_zzz_xxxx_0[i] * pa_x[i] - ta1_y_zzz_xxxx_1[i] * pc_x[i];

        ta1_y_xzzz_xxxy_0[i] =
            3.0 * ta1_y_zzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxy_1[i] * fe_0 + ta1_y_zzz_xxxy_0[i] * pa_x[i] - ta1_y_zzz_xxxy_1[i] * pc_x[i];

        ta1_y_xzzz_xxxz_0[i] =
            3.0 * ta1_y_zzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxz_1[i] * fe_0 + ta1_y_zzz_xxxz_0[i] * pa_x[i] - ta1_y_zzz_xxxz_1[i] * pc_x[i];

        ta1_y_xzzz_xxyy_0[i] =
            2.0 * ta1_y_zzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyy_1[i] * fe_0 + ta1_y_zzz_xxyy_0[i] * pa_x[i] - ta1_y_zzz_xxyy_1[i] * pc_x[i];

        ta1_y_xzzz_xxyz_0[i] =
            2.0 * ta1_y_zzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyz_1[i] * fe_0 + ta1_y_zzz_xxyz_0[i] * pa_x[i] - ta1_y_zzz_xxyz_1[i] * pc_x[i];

        ta1_y_xzzz_xxzz_0[i] =
            2.0 * ta1_y_zzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xzz_1[i] * fe_0 + ta1_y_zzz_xxzz_0[i] * pa_x[i] - ta1_y_zzz_xxzz_1[i] * pc_x[i];

        ta1_y_xzzz_xyyy_0[i] = ta1_y_zzz_yyy_0[i] * fe_0 - ta1_y_zzz_yyy_1[i] * fe_0 + ta1_y_zzz_xyyy_0[i] * pa_x[i] - ta1_y_zzz_xyyy_1[i] * pc_x[i];

        ta1_y_xzzz_xyyz_0[i] = ta1_y_zzz_yyz_0[i] * fe_0 - ta1_y_zzz_yyz_1[i] * fe_0 + ta1_y_zzz_xyyz_0[i] * pa_x[i] - ta1_y_zzz_xyyz_1[i] * pc_x[i];

        ta1_y_xzzz_xyzz_0[i] = ta1_y_zzz_yzz_0[i] * fe_0 - ta1_y_zzz_yzz_1[i] * fe_0 + ta1_y_zzz_xyzz_0[i] * pa_x[i] - ta1_y_zzz_xyzz_1[i] * pc_x[i];

        ta1_y_xzzz_xzzz_0[i] = ta1_y_zzz_zzz_0[i] * fe_0 - ta1_y_zzz_zzz_1[i] * fe_0 + ta1_y_zzz_xzzz_0[i] * pa_x[i] - ta1_y_zzz_xzzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyyy_0[i] = ta1_y_zzz_yyyy_0[i] * pa_x[i] - ta1_y_zzz_yyyy_1[i] * pc_x[i];

        ta1_y_xzzz_yyyz_0[i] = ta1_y_zzz_yyyz_0[i] * pa_x[i] - ta1_y_zzz_yyyz_1[i] * pc_x[i];

        ta1_y_xzzz_yyzz_0[i] = ta1_y_zzz_yyzz_0[i] * pa_x[i] - ta1_y_zzz_yyzz_1[i] * pc_x[i];

        ta1_y_xzzz_yzzz_0[i] = ta1_y_zzz_yzzz_0[i] * pa_x[i] - ta1_y_zzz_yzzz_1[i] * pc_x[i];

        ta1_y_xzzz_zzzz_0[i] = ta1_y_zzz_zzzz_0[i] * pa_x[i] - ta1_y_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 375-390 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_y_yy_xxxx_0,   \
                             ta1_y_yy_xxxx_1,   \
                             ta1_y_yy_xxxy_0,   \
                             ta1_y_yy_xxxy_1,   \
                             ta1_y_yy_xxxz_0,   \
                             ta1_y_yy_xxxz_1,   \
                             ta1_y_yy_xxyy_0,   \
                             ta1_y_yy_xxyy_1,   \
                             ta1_y_yy_xxyz_0,   \
                             ta1_y_yy_xxyz_1,   \
                             ta1_y_yy_xxzz_0,   \
                             ta1_y_yy_xxzz_1,   \
                             ta1_y_yy_xyyy_0,   \
                             ta1_y_yy_xyyy_1,   \
                             ta1_y_yy_xyyz_0,   \
                             ta1_y_yy_xyyz_1,   \
                             ta1_y_yy_xyzz_0,   \
                             ta1_y_yy_xyzz_1,   \
                             ta1_y_yy_xzzz_0,   \
                             ta1_y_yy_xzzz_1,   \
                             ta1_y_yy_yyyy_0,   \
                             ta1_y_yy_yyyy_1,   \
                             ta1_y_yy_yyyz_0,   \
                             ta1_y_yy_yyyz_1,   \
                             ta1_y_yy_yyzz_0,   \
                             ta1_y_yy_yyzz_1,   \
                             ta1_y_yy_yzzz_0,   \
                             ta1_y_yy_yzzz_1,   \
                             ta1_y_yy_zzzz_0,   \
                             ta1_y_yy_zzzz_1,   \
                             ta1_y_yyy_xxx_0,   \
                             ta1_y_yyy_xxx_1,   \
                             ta1_y_yyy_xxxx_0,  \
                             ta1_y_yyy_xxxx_1,  \
                             ta1_y_yyy_xxxy_0,  \
                             ta1_y_yyy_xxxy_1,  \
                             ta1_y_yyy_xxxz_0,  \
                             ta1_y_yyy_xxxz_1,  \
                             ta1_y_yyy_xxy_0,   \
                             ta1_y_yyy_xxy_1,   \
                             ta1_y_yyy_xxyy_0,  \
                             ta1_y_yyy_xxyy_1,  \
                             ta1_y_yyy_xxyz_0,  \
                             ta1_y_yyy_xxyz_1,  \
                             ta1_y_yyy_xxz_0,   \
                             ta1_y_yyy_xxz_1,   \
                             ta1_y_yyy_xxzz_0,  \
                             ta1_y_yyy_xxzz_1,  \
                             ta1_y_yyy_xyy_0,   \
                             ta1_y_yyy_xyy_1,   \
                             ta1_y_yyy_xyyy_0,  \
                             ta1_y_yyy_xyyy_1,  \
                             ta1_y_yyy_xyyz_0,  \
                             ta1_y_yyy_xyyz_1,  \
                             ta1_y_yyy_xyz_0,   \
                             ta1_y_yyy_xyz_1,   \
                             ta1_y_yyy_xyzz_0,  \
                             ta1_y_yyy_xyzz_1,  \
                             ta1_y_yyy_xzz_0,   \
                             ta1_y_yyy_xzz_1,   \
                             ta1_y_yyy_xzzz_0,  \
                             ta1_y_yyy_xzzz_1,  \
                             ta1_y_yyy_yyy_0,   \
                             ta1_y_yyy_yyy_1,   \
                             ta1_y_yyy_yyyy_0,  \
                             ta1_y_yyy_yyyy_1,  \
                             ta1_y_yyy_yyyz_0,  \
                             ta1_y_yyy_yyyz_1,  \
                             ta1_y_yyy_yyz_0,   \
                             ta1_y_yyy_yyz_1,   \
                             ta1_y_yyy_yyzz_0,  \
                             ta1_y_yyy_yyzz_1,  \
                             ta1_y_yyy_yzz_0,   \
                             ta1_y_yyy_yzz_1,   \
                             ta1_y_yyy_yzzz_0,  \
                             ta1_y_yyy_yzzz_1,  \
                             ta1_y_yyy_zzz_0,   \
                             ta1_y_yyy_zzz_1,   \
                             ta1_y_yyy_zzzz_0,  \
                             ta1_y_yyy_zzzz_1,  \
                             ta1_y_yyyy_xxxx_0, \
                             ta1_y_yyyy_xxxy_0, \
                             ta1_y_yyyy_xxxz_0, \
                             ta1_y_yyyy_xxyy_0, \
                             ta1_y_yyyy_xxyz_0, \
                             ta1_y_yyyy_xxzz_0, \
                             ta1_y_yyyy_xyyy_0, \
                             ta1_y_yyyy_xyyz_0, \
                             ta1_y_yyyy_xyzz_0, \
                             ta1_y_yyyy_xzzz_0, \
                             ta1_y_yyyy_yyyy_0, \
                             ta1_y_yyyy_yyyz_0, \
                             ta1_y_yyyy_yyzz_0, \
                             ta1_y_yyyy_yzzz_0, \
                             ta1_y_yyyy_zzzz_0, \
                             ta_yyy_xxxx_1,     \
                             ta_yyy_xxxy_1,     \
                             ta_yyy_xxxz_1,     \
                             ta_yyy_xxyy_1,     \
                             ta_yyy_xxyz_1,     \
                             ta_yyy_xxzz_1,     \
                             ta_yyy_xyyy_1,     \
                             ta_yyy_xyyz_1,     \
                             ta_yyy_xyzz_1,     \
                             ta_yyy_xzzz_1,     \
                             ta_yyy_yyyy_1,     \
                             ta_yyy_yyyz_1,     \
                             ta_yyy_yyzz_1,     \
                             ta_yyy_yzzz_1,     \
                             ta_yyy_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyy_xxxx_0[i] = 3.0 * ta1_y_yy_xxxx_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxx_1[i] * fe_0 + ta_yyy_xxxx_1[i] + ta1_y_yyy_xxxx_0[i] * pa_y[i] -
                               ta1_y_yyy_xxxx_1[i] * pc_y[i];

        ta1_y_yyyy_xxxy_0[i] = 3.0 * ta1_y_yy_xxxy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxy_1[i] * fe_0 + ta1_y_yyy_xxx_0[i] * fe_0 -
                               ta1_y_yyy_xxx_1[i] * fe_0 + ta_yyy_xxxy_1[i] + ta1_y_yyy_xxxy_0[i] * pa_y[i] - ta1_y_yyy_xxxy_1[i] * pc_y[i];

        ta1_y_yyyy_xxxz_0[i] = 3.0 * ta1_y_yy_xxxz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxz_1[i] * fe_0 + ta_yyy_xxxz_1[i] + ta1_y_yyy_xxxz_0[i] * pa_y[i] -
                               ta1_y_yyy_xxxz_1[i] * pc_y[i];

        ta1_y_yyyy_xxyy_0[i] = 3.0 * ta1_y_yy_xxyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxy_0[i] * fe_0 -
                               2.0 * ta1_y_yyy_xxy_1[i] * fe_0 + ta_yyy_xxyy_1[i] + ta1_y_yyy_xxyy_0[i] * pa_y[i] - ta1_y_yyy_xxyy_1[i] * pc_y[i];

        ta1_y_yyyy_xxyz_0[i] = 3.0 * ta1_y_yy_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyz_1[i] * fe_0 + ta1_y_yyy_xxz_0[i] * fe_0 -
                               ta1_y_yyy_xxz_1[i] * fe_0 + ta_yyy_xxyz_1[i] + ta1_y_yyy_xxyz_0[i] * pa_y[i] - ta1_y_yyy_xxyz_1[i] * pc_y[i];

        ta1_y_yyyy_xxzz_0[i] = 3.0 * ta1_y_yy_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxzz_1[i] * fe_0 + ta_yyy_xxzz_1[i] + ta1_y_yyy_xxzz_0[i] * pa_y[i] -
                               ta1_y_yyy_xxzz_1[i] * pc_y[i];

        ta1_y_yyyy_xyyy_0[i] = 3.0 * ta1_y_yy_xyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyy_1[i] * fe_0 + 3.0 * ta1_y_yyy_xyy_0[i] * fe_0 -
                               3.0 * ta1_y_yyy_xyy_1[i] * fe_0 + ta_yyy_xyyy_1[i] + ta1_y_yyy_xyyy_0[i] * pa_y[i] - ta1_y_yyy_xyyy_1[i] * pc_y[i];

        ta1_y_yyyy_xyyz_0[i] = 3.0 * ta1_y_yy_xyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyz_0[i] * fe_0 -
                               2.0 * ta1_y_yyy_xyz_1[i] * fe_0 + ta_yyy_xyyz_1[i] + ta1_y_yyy_xyyz_0[i] * pa_y[i] - ta1_y_yyy_xyyz_1[i] * pc_y[i];

        ta1_y_yyyy_xyzz_0[i] = 3.0 * ta1_y_yy_xyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyzz_1[i] * fe_0 + ta1_y_yyy_xzz_0[i] * fe_0 -
                               ta1_y_yyy_xzz_1[i] * fe_0 + ta_yyy_xyzz_1[i] + ta1_y_yyy_xyzz_0[i] * pa_y[i] - ta1_y_yyy_xyzz_1[i] * pc_y[i];

        ta1_y_yyyy_xzzz_0[i] = 3.0 * ta1_y_yy_xzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xzzz_1[i] * fe_0 + ta_yyy_xzzz_1[i] + ta1_y_yyy_xzzz_0[i] * pa_y[i] -
                               ta1_y_yyy_xzzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyyy_0[i] = 3.0 * ta1_y_yy_yyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyy_1[i] * fe_0 + 4.0 * ta1_y_yyy_yyy_0[i] * fe_0 -
                               4.0 * ta1_y_yyy_yyy_1[i] * fe_0 + ta_yyy_yyyy_1[i] + ta1_y_yyy_yyyy_0[i] * pa_y[i] - ta1_y_yyy_yyyy_1[i] * pc_y[i];

        ta1_y_yyyy_yyyz_0[i] = 3.0 * ta1_y_yy_yyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyz_1[i] * fe_0 + 3.0 * ta1_y_yyy_yyz_0[i] * fe_0 -
                               3.0 * ta1_y_yyy_yyz_1[i] * fe_0 + ta_yyy_yyyz_1[i] + ta1_y_yyy_yyyz_0[i] * pa_y[i] - ta1_y_yyy_yyyz_1[i] * pc_y[i];

        ta1_y_yyyy_yyzz_0[i] = 3.0 * ta1_y_yy_yyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yzz_0[i] * fe_0 -
                               2.0 * ta1_y_yyy_yzz_1[i] * fe_0 + ta_yyy_yyzz_1[i] + ta1_y_yyy_yyzz_0[i] * pa_y[i] - ta1_y_yyy_yyzz_1[i] * pc_y[i];

        ta1_y_yyyy_yzzz_0[i] = 3.0 * ta1_y_yy_yzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yzzz_1[i] * fe_0 + ta1_y_yyy_zzz_0[i] * fe_0 -
                               ta1_y_yyy_zzz_1[i] * fe_0 + ta_yyy_yzzz_1[i] + ta1_y_yyy_yzzz_0[i] * pa_y[i] - ta1_y_yyy_yzzz_1[i] * pc_y[i];

        ta1_y_yyyy_zzzz_0[i] = 3.0 * ta1_y_yy_zzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_zzzz_1[i] * fe_0 + ta_yyy_zzzz_1[i] + ta1_y_yyy_zzzz_0[i] * pa_y[i] -
                               ta1_y_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 390-405 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_y_yyy_xxx_0,   \
                             ta1_y_yyy_xxx_1,   \
                             ta1_y_yyy_xxxx_0,  \
                             ta1_y_yyy_xxxx_1,  \
                             ta1_y_yyy_xxxy_0,  \
                             ta1_y_yyy_xxxy_1,  \
                             ta1_y_yyy_xxxz_0,  \
                             ta1_y_yyy_xxxz_1,  \
                             ta1_y_yyy_xxy_0,   \
                             ta1_y_yyy_xxy_1,   \
                             ta1_y_yyy_xxyy_0,  \
                             ta1_y_yyy_xxyy_1,  \
                             ta1_y_yyy_xxyz_0,  \
                             ta1_y_yyy_xxyz_1,  \
                             ta1_y_yyy_xxz_0,   \
                             ta1_y_yyy_xxz_1,   \
                             ta1_y_yyy_xxzz_0,  \
                             ta1_y_yyy_xxzz_1,  \
                             ta1_y_yyy_xyy_0,   \
                             ta1_y_yyy_xyy_1,   \
                             ta1_y_yyy_xyyy_0,  \
                             ta1_y_yyy_xyyy_1,  \
                             ta1_y_yyy_xyyz_0,  \
                             ta1_y_yyy_xyyz_1,  \
                             ta1_y_yyy_xyz_0,   \
                             ta1_y_yyy_xyz_1,   \
                             ta1_y_yyy_xyzz_0,  \
                             ta1_y_yyy_xyzz_1,  \
                             ta1_y_yyy_xzz_0,   \
                             ta1_y_yyy_xzz_1,   \
                             ta1_y_yyy_xzzz_0,  \
                             ta1_y_yyy_xzzz_1,  \
                             ta1_y_yyy_yyy_0,   \
                             ta1_y_yyy_yyy_1,   \
                             ta1_y_yyy_yyyy_0,  \
                             ta1_y_yyy_yyyy_1,  \
                             ta1_y_yyy_yyyz_0,  \
                             ta1_y_yyy_yyyz_1,  \
                             ta1_y_yyy_yyz_0,   \
                             ta1_y_yyy_yyz_1,   \
                             ta1_y_yyy_yyzz_0,  \
                             ta1_y_yyy_yyzz_1,  \
                             ta1_y_yyy_yzz_0,   \
                             ta1_y_yyy_yzz_1,   \
                             ta1_y_yyy_yzzz_0,  \
                             ta1_y_yyy_yzzz_1,  \
                             ta1_y_yyy_zzz_0,   \
                             ta1_y_yyy_zzz_1,   \
                             ta1_y_yyy_zzzz_0,  \
                             ta1_y_yyy_zzzz_1,  \
                             ta1_y_yyyz_xxxx_0, \
                             ta1_y_yyyz_xxxy_0, \
                             ta1_y_yyyz_xxxz_0, \
                             ta1_y_yyyz_xxyy_0, \
                             ta1_y_yyyz_xxyz_0, \
                             ta1_y_yyyz_xxzz_0, \
                             ta1_y_yyyz_xyyy_0, \
                             ta1_y_yyyz_xyyz_0, \
                             ta1_y_yyyz_xyzz_0, \
                             ta1_y_yyyz_xzzz_0, \
                             ta1_y_yyyz_yyyy_0, \
                             ta1_y_yyyz_yyyz_0, \
                             ta1_y_yyyz_yyzz_0, \
                             ta1_y_yyyz_yzzz_0, \
                             ta1_y_yyyz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyz_xxxx_0[i] = ta1_y_yyy_xxxx_0[i] * pa_z[i] - ta1_y_yyy_xxxx_1[i] * pc_z[i];

        ta1_y_yyyz_xxxy_0[i] = ta1_y_yyy_xxxy_0[i] * pa_z[i] - ta1_y_yyy_xxxy_1[i] * pc_z[i];

        ta1_y_yyyz_xxxz_0[i] = ta1_y_yyy_xxx_0[i] * fe_0 - ta1_y_yyy_xxx_1[i] * fe_0 + ta1_y_yyy_xxxz_0[i] * pa_z[i] - ta1_y_yyy_xxxz_1[i] * pc_z[i];

        ta1_y_yyyz_xxyy_0[i] = ta1_y_yyy_xxyy_0[i] * pa_z[i] - ta1_y_yyy_xxyy_1[i] * pc_z[i];

        ta1_y_yyyz_xxyz_0[i] = ta1_y_yyy_xxy_0[i] * fe_0 - ta1_y_yyy_xxy_1[i] * fe_0 + ta1_y_yyy_xxyz_0[i] * pa_z[i] - ta1_y_yyy_xxyz_1[i] * pc_z[i];

        ta1_y_yyyz_xxzz_0[i] =
            2.0 * ta1_y_yyy_xxz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xxz_1[i] * fe_0 + ta1_y_yyy_xxzz_0[i] * pa_z[i] - ta1_y_yyy_xxzz_1[i] * pc_z[i];

        ta1_y_yyyz_xyyy_0[i] = ta1_y_yyy_xyyy_0[i] * pa_z[i] - ta1_y_yyy_xyyy_1[i] * pc_z[i];

        ta1_y_yyyz_xyyz_0[i] = ta1_y_yyy_xyy_0[i] * fe_0 - ta1_y_yyy_xyy_1[i] * fe_0 + ta1_y_yyy_xyyz_0[i] * pa_z[i] - ta1_y_yyy_xyyz_1[i] * pc_z[i];

        ta1_y_yyyz_xyzz_0[i] =
            2.0 * ta1_y_yyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyz_1[i] * fe_0 + ta1_y_yyy_xyzz_0[i] * pa_z[i] - ta1_y_yyy_xyzz_1[i] * pc_z[i];

        ta1_y_yyyz_xzzz_0[i] =
            3.0 * ta1_y_yyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xzz_1[i] * fe_0 + ta1_y_yyy_xzzz_0[i] * pa_z[i] - ta1_y_yyy_xzzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyyy_0[i] = ta1_y_yyy_yyyy_0[i] * pa_z[i] - ta1_y_yyy_yyyy_1[i] * pc_z[i];

        ta1_y_yyyz_yyyz_0[i] = ta1_y_yyy_yyy_0[i] * fe_0 - ta1_y_yyy_yyy_1[i] * fe_0 + ta1_y_yyy_yyyz_0[i] * pa_z[i] - ta1_y_yyy_yyyz_1[i] * pc_z[i];

        ta1_y_yyyz_yyzz_0[i] =
            2.0 * ta1_y_yyy_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_yyz_1[i] * fe_0 + ta1_y_yyy_yyzz_0[i] * pa_z[i] - ta1_y_yyy_yyzz_1[i] * pc_z[i];

        ta1_y_yyyz_yzzz_0[i] =
            3.0 * ta1_y_yyy_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_yzz_1[i] * fe_0 + ta1_y_yyy_yzzz_0[i] * pa_z[i] - ta1_y_yyy_yzzz_1[i] * pc_z[i];

        ta1_y_yyyz_zzzz_0[i] =
            4.0 * ta1_y_yyy_zzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_zzz_1[i] * fe_0 + ta1_y_yyy_zzzz_0[i] * pa_z[i] - ta1_y_yyy_zzzz_1[i] * pc_z[i];
    }

    // Set up 405-420 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_yy_xxxx_0,   \
                             ta1_y_yy_xxxx_1,   \
                             ta1_y_yy_xxxy_0,   \
                             ta1_y_yy_xxxy_1,   \
                             ta1_y_yy_xxyy_0,   \
                             ta1_y_yy_xxyy_1,   \
                             ta1_y_yy_xxyz_0,   \
                             ta1_y_yy_xxyz_1,   \
                             ta1_y_yy_xyyy_0,   \
                             ta1_y_yy_xyyy_1,   \
                             ta1_y_yy_xyyz_0,   \
                             ta1_y_yy_xyyz_1,   \
                             ta1_y_yy_xyzz_0,   \
                             ta1_y_yy_xyzz_1,   \
                             ta1_y_yy_yyyy_0,   \
                             ta1_y_yy_yyyy_1,   \
                             ta1_y_yy_yyyz_0,   \
                             ta1_y_yy_yyyz_1,   \
                             ta1_y_yy_yyzz_0,   \
                             ta1_y_yy_yyzz_1,   \
                             ta1_y_yy_yzzz_0,   \
                             ta1_y_yy_yzzz_1,   \
                             ta1_y_yyz_xxxx_0,  \
                             ta1_y_yyz_xxxx_1,  \
                             ta1_y_yyz_xxxy_0,  \
                             ta1_y_yyz_xxxy_1,  \
                             ta1_y_yyz_xxy_0,   \
                             ta1_y_yyz_xxy_1,   \
                             ta1_y_yyz_xxyy_0,  \
                             ta1_y_yyz_xxyy_1,  \
                             ta1_y_yyz_xxyz_0,  \
                             ta1_y_yyz_xxyz_1,  \
                             ta1_y_yyz_xyy_0,   \
                             ta1_y_yyz_xyy_1,   \
                             ta1_y_yyz_xyyy_0,  \
                             ta1_y_yyz_xyyy_1,  \
                             ta1_y_yyz_xyyz_0,  \
                             ta1_y_yyz_xyyz_1,  \
                             ta1_y_yyz_xyz_0,   \
                             ta1_y_yyz_xyz_1,   \
                             ta1_y_yyz_xyzz_0,  \
                             ta1_y_yyz_xyzz_1,  \
                             ta1_y_yyz_yyy_0,   \
                             ta1_y_yyz_yyy_1,   \
                             ta1_y_yyz_yyyy_0,  \
                             ta1_y_yyz_yyyy_1,  \
                             ta1_y_yyz_yyyz_0,  \
                             ta1_y_yyz_yyyz_1,  \
                             ta1_y_yyz_yyz_0,   \
                             ta1_y_yyz_yyz_1,   \
                             ta1_y_yyz_yyzz_0,  \
                             ta1_y_yyz_yyzz_1,  \
                             ta1_y_yyz_yzz_0,   \
                             ta1_y_yyz_yzz_1,   \
                             ta1_y_yyz_yzzz_0,  \
                             ta1_y_yyz_yzzz_1,  \
                             ta1_y_yyzz_xxxx_0, \
                             ta1_y_yyzz_xxxy_0, \
                             ta1_y_yyzz_xxxz_0, \
                             ta1_y_yyzz_xxyy_0, \
                             ta1_y_yyzz_xxyz_0, \
                             ta1_y_yyzz_xxzz_0, \
                             ta1_y_yyzz_xyyy_0, \
                             ta1_y_yyzz_xyyz_0, \
                             ta1_y_yyzz_xyzz_0, \
                             ta1_y_yyzz_xzzz_0, \
                             ta1_y_yyzz_yyyy_0, \
                             ta1_y_yyzz_yyyz_0, \
                             ta1_y_yyzz_yyzz_0, \
                             ta1_y_yyzz_yzzz_0, \
                             ta1_y_yyzz_zzzz_0, \
                             ta1_y_yzz_xxxz_0,  \
                             ta1_y_yzz_xxxz_1,  \
                             ta1_y_yzz_xxzz_0,  \
                             ta1_y_yzz_xxzz_1,  \
                             ta1_y_yzz_xzzz_0,  \
                             ta1_y_yzz_xzzz_1,  \
                             ta1_y_yzz_zzzz_0,  \
                             ta1_y_yzz_zzzz_1,  \
                             ta1_y_zz_xxxz_0,   \
                             ta1_y_zz_xxxz_1,   \
                             ta1_y_zz_xxzz_0,   \
                             ta1_y_zz_xxzz_1,   \
                             ta1_y_zz_xzzz_0,   \
                             ta1_y_zz_xzzz_1,   \
                             ta1_y_zz_zzzz_0,   \
                             ta1_y_zz_zzzz_1,   \
                             ta_yzz_xxxz_1,     \
                             ta_yzz_xxzz_1,     \
                             ta_yzz_xzzz_1,     \
                             ta_yzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzz_xxxx_0[i] = ta1_y_yy_xxxx_0[i] * fe_0 - ta1_y_yy_xxxx_1[i] * fe_0 + ta1_y_yyz_xxxx_0[i] * pa_z[i] - ta1_y_yyz_xxxx_1[i] * pc_z[i];

        ta1_y_yyzz_xxxy_0[i] = ta1_y_yy_xxxy_0[i] * fe_0 - ta1_y_yy_xxxy_1[i] * fe_0 + ta1_y_yyz_xxxy_0[i] * pa_z[i] - ta1_y_yyz_xxxy_1[i] * pc_z[i];

        ta1_y_yyzz_xxxz_0[i] =
            ta1_y_zz_xxxz_0[i] * fe_0 - ta1_y_zz_xxxz_1[i] * fe_0 + ta_yzz_xxxz_1[i] + ta1_y_yzz_xxxz_0[i] * pa_y[i] - ta1_y_yzz_xxxz_1[i] * pc_y[i];

        ta1_y_yyzz_xxyy_0[i] = ta1_y_yy_xxyy_0[i] * fe_0 - ta1_y_yy_xxyy_1[i] * fe_0 + ta1_y_yyz_xxyy_0[i] * pa_z[i] - ta1_y_yyz_xxyy_1[i] * pc_z[i];

        ta1_y_yyzz_xxyz_0[i] = ta1_y_yy_xxyz_0[i] * fe_0 - ta1_y_yy_xxyz_1[i] * fe_0 + ta1_y_yyz_xxy_0[i] * fe_0 - ta1_y_yyz_xxy_1[i] * fe_0 +
                               ta1_y_yyz_xxyz_0[i] * pa_z[i] - ta1_y_yyz_xxyz_1[i] * pc_z[i];

        ta1_y_yyzz_xxzz_0[i] =
            ta1_y_zz_xxzz_0[i] * fe_0 - ta1_y_zz_xxzz_1[i] * fe_0 + ta_yzz_xxzz_1[i] + ta1_y_yzz_xxzz_0[i] * pa_y[i] - ta1_y_yzz_xxzz_1[i] * pc_y[i];

        ta1_y_yyzz_xyyy_0[i] = ta1_y_yy_xyyy_0[i] * fe_0 - ta1_y_yy_xyyy_1[i] * fe_0 + ta1_y_yyz_xyyy_0[i] * pa_z[i] - ta1_y_yyz_xyyy_1[i] * pc_z[i];

        ta1_y_yyzz_xyyz_0[i] = ta1_y_yy_xyyz_0[i] * fe_0 - ta1_y_yy_xyyz_1[i] * fe_0 + ta1_y_yyz_xyy_0[i] * fe_0 - ta1_y_yyz_xyy_1[i] * fe_0 +
                               ta1_y_yyz_xyyz_0[i] * pa_z[i] - ta1_y_yyz_xyyz_1[i] * pc_z[i];

        ta1_y_yyzz_xyzz_0[i] = ta1_y_yy_xyzz_0[i] * fe_0 - ta1_y_yy_xyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_xyz_0[i] * fe_0 -
                               2.0 * ta1_y_yyz_xyz_1[i] * fe_0 + ta1_y_yyz_xyzz_0[i] * pa_z[i] - ta1_y_yyz_xyzz_1[i] * pc_z[i];

        ta1_y_yyzz_xzzz_0[i] =
            ta1_y_zz_xzzz_0[i] * fe_0 - ta1_y_zz_xzzz_1[i] * fe_0 + ta_yzz_xzzz_1[i] + ta1_y_yzz_xzzz_0[i] * pa_y[i] - ta1_y_yzz_xzzz_1[i] * pc_y[i];

        ta1_y_yyzz_yyyy_0[i] = ta1_y_yy_yyyy_0[i] * fe_0 - ta1_y_yy_yyyy_1[i] * fe_0 + ta1_y_yyz_yyyy_0[i] * pa_z[i] - ta1_y_yyz_yyyy_1[i] * pc_z[i];

        ta1_y_yyzz_yyyz_0[i] = ta1_y_yy_yyyz_0[i] * fe_0 - ta1_y_yy_yyyz_1[i] * fe_0 + ta1_y_yyz_yyy_0[i] * fe_0 - ta1_y_yyz_yyy_1[i] * fe_0 +
                               ta1_y_yyz_yyyz_0[i] * pa_z[i] - ta1_y_yyz_yyyz_1[i] * pc_z[i];

        ta1_y_yyzz_yyzz_0[i] = ta1_y_yy_yyzz_0[i] * fe_0 - ta1_y_yy_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_yyz_0[i] * fe_0 -
                               2.0 * ta1_y_yyz_yyz_1[i] * fe_0 + ta1_y_yyz_yyzz_0[i] * pa_z[i] - ta1_y_yyz_yyzz_1[i] * pc_z[i];

        ta1_y_yyzz_yzzz_0[i] = ta1_y_yy_yzzz_0[i] * fe_0 - ta1_y_yy_yzzz_1[i] * fe_0 + 3.0 * ta1_y_yyz_yzz_0[i] * fe_0 -
                               3.0 * ta1_y_yyz_yzz_1[i] * fe_0 + ta1_y_yyz_yzzz_0[i] * pa_z[i] - ta1_y_yyz_yzzz_1[i] * pc_z[i];

        ta1_y_yyzz_zzzz_0[i] =
            ta1_y_zz_zzzz_0[i] * fe_0 - ta1_y_zz_zzzz_1[i] * fe_0 + ta_yzz_zzzz_1[i] + ta1_y_yzz_zzzz_0[i] * pa_y[i] - ta1_y_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 420-435 components of targeted buffer : GG

    auto ta1_y_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 420);

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_yz_xxxy_0,   \
                             ta1_y_yz_xxxy_1,   \
                             ta1_y_yz_xxyy_0,   \
                             ta1_y_yz_xxyy_1,   \
                             ta1_y_yz_xyyy_0,   \
                             ta1_y_yz_xyyy_1,   \
                             ta1_y_yz_yyyy_0,   \
                             ta1_y_yz_yyyy_1,   \
                             ta1_y_yzz_xxxy_0,  \
                             ta1_y_yzz_xxxy_1,  \
                             ta1_y_yzz_xxyy_0,  \
                             ta1_y_yzz_xxyy_1,  \
                             ta1_y_yzz_xyyy_0,  \
                             ta1_y_yzz_xyyy_1,  \
                             ta1_y_yzz_yyyy_0,  \
                             ta1_y_yzz_yyyy_1,  \
                             ta1_y_yzzz_xxxx_0, \
                             ta1_y_yzzz_xxxy_0, \
                             ta1_y_yzzz_xxxz_0, \
                             ta1_y_yzzz_xxyy_0, \
                             ta1_y_yzzz_xxyz_0, \
                             ta1_y_yzzz_xxzz_0, \
                             ta1_y_yzzz_xyyy_0, \
                             ta1_y_yzzz_xyyz_0, \
                             ta1_y_yzzz_xyzz_0, \
                             ta1_y_yzzz_xzzz_0, \
                             ta1_y_yzzz_yyyy_0, \
                             ta1_y_yzzz_yyyz_0, \
                             ta1_y_yzzz_yyzz_0, \
                             ta1_y_yzzz_yzzz_0, \
                             ta1_y_yzzz_zzzz_0, \
                             ta1_y_zzz_xxxx_0,  \
                             ta1_y_zzz_xxxx_1,  \
                             ta1_y_zzz_xxxz_0,  \
                             ta1_y_zzz_xxxz_1,  \
                             ta1_y_zzz_xxyz_0,  \
                             ta1_y_zzz_xxyz_1,  \
                             ta1_y_zzz_xxz_0,   \
                             ta1_y_zzz_xxz_1,   \
                             ta1_y_zzz_xxzz_0,  \
                             ta1_y_zzz_xxzz_1,  \
                             ta1_y_zzz_xyyz_0,  \
                             ta1_y_zzz_xyyz_1,  \
                             ta1_y_zzz_xyz_0,   \
                             ta1_y_zzz_xyz_1,   \
                             ta1_y_zzz_xyzz_0,  \
                             ta1_y_zzz_xyzz_1,  \
                             ta1_y_zzz_xzz_0,   \
                             ta1_y_zzz_xzz_1,   \
                             ta1_y_zzz_xzzz_0,  \
                             ta1_y_zzz_xzzz_1,  \
                             ta1_y_zzz_yyyz_0,  \
                             ta1_y_zzz_yyyz_1,  \
                             ta1_y_zzz_yyz_0,   \
                             ta1_y_zzz_yyz_1,   \
                             ta1_y_zzz_yyzz_0,  \
                             ta1_y_zzz_yyzz_1,  \
                             ta1_y_zzz_yzz_0,   \
                             ta1_y_zzz_yzz_1,   \
                             ta1_y_zzz_yzzz_0,  \
                             ta1_y_zzz_yzzz_1,  \
                             ta1_y_zzz_zzz_0,   \
                             ta1_y_zzz_zzz_1,   \
                             ta1_y_zzz_zzzz_0,  \
                             ta1_y_zzz_zzzz_1,  \
                             ta_zzz_xxxx_1,     \
                             ta_zzz_xxxz_1,     \
                             ta_zzz_xxyz_1,     \
                             ta_zzz_xxzz_1,     \
                             ta_zzz_xyyz_1,     \
                             ta_zzz_xyzz_1,     \
                             ta_zzz_xzzz_1,     \
                             ta_zzz_yyyz_1,     \
                             ta_zzz_yyzz_1,     \
                             ta_zzz_yzzz_1,     \
                             ta_zzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzz_xxxx_0[i] = ta_zzz_xxxx_1[i] + ta1_y_zzz_xxxx_0[i] * pa_y[i] - ta1_y_zzz_xxxx_1[i] * pc_y[i];

        ta1_y_yzzz_xxxy_0[i] =
            2.0 * ta1_y_yz_xxxy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxxy_1[i] * fe_0 + ta1_y_yzz_xxxy_0[i] * pa_z[i] - ta1_y_yzz_xxxy_1[i] * pc_z[i];

        ta1_y_yzzz_xxxz_0[i] = ta_zzz_xxxz_1[i] + ta1_y_zzz_xxxz_0[i] * pa_y[i] - ta1_y_zzz_xxxz_1[i] * pc_y[i];

        ta1_y_yzzz_xxyy_0[i] =
            2.0 * ta1_y_yz_xxyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxyy_1[i] * fe_0 + ta1_y_yzz_xxyy_0[i] * pa_z[i] - ta1_y_yzz_xxyy_1[i] * pc_z[i];

        ta1_y_yzzz_xxyz_0[i] =
            ta1_y_zzz_xxz_0[i] * fe_0 - ta1_y_zzz_xxz_1[i] * fe_0 + ta_zzz_xxyz_1[i] + ta1_y_zzz_xxyz_0[i] * pa_y[i] - ta1_y_zzz_xxyz_1[i] * pc_y[i];

        ta1_y_yzzz_xxzz_0[i] = ta_zzz_xxzz_1[i] + ta1_y_zzz_xxzz_0[i] * pa_y[i] - ta1_y_zzz_xxzz_1[i] * pc_y[i];

        ta1_y_yzzz_xyyy_0[i] =
            2.0 * ta1_y_yz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xyyy_1[i] * fe_0 + ta1_y_yzz_xyyy_0[i] * pa_z[i] - ta1_y_yzz_xyyy_1[i] * pc_z[i];

        ta1_y_yzzz_xyyz_0[i] = 2.0 * ta1_y_zzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyz_1[i] * fe_0 + ta_zzz_xyyz_1[i] + ta1_y_zzz_xyyz_0[i] * pa_y[i] -
                               ta1_y_zzz_xyyz_1[i] * pc_y[i];

        ta1_y_yzzz_xyzz_0[i] =
            ta1_y_zzz_xzz_0[i] * fe_0 - ta1_y_zzz_xzz_1[i] * fe_0 + ta_zzz_xyzz_1[i] + ta1_y_zzz_xyzz_0[i] * pa_y[i] - ta1_y_zzz_xyzz_1[i] * pc_y[i];

        ta1_y_yzzz_xzzz_0[i] = ta_zzz_xzzz_1[i] + ta1_y_zzz_xzzz_0[i] * pa_y[i] - ta1_y_zzz_xzzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyyy_0[i] =
            2.0 * ta1_y_yz_yyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_yyyy_1[i] * fe_0 + ta1_y_yzz_yyyy_0[i] * pa_z[i] - ta1_y_yzz_yyyy_1[i] * pc_z[i];

        ta1_y_yzzz_yyyz_0[i] = 3.0 * ta1_y_zzz_yyz_0[i] * fe_0 - 3.0 * ta1_y_zzz_yyz_1[i] * fe_0 + ta_zzz_yyyz_1[i] + ta1_y_zzz_yyyz_0[i] * pa_y[i] -
                               ta1_y_zzz_yyyz_1[i] * pc_y[i];

        ta1_y_yzzz_yyzz_0[i] = 2.0 * ta1_y_zzz_yzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_yzz_1[i] * fe_0 + ta_zzz_yyzz_1[i] + ta1_y_zzz_yyzz_0[i] * pa_y[i] -
                               ta1_y_zzz_yyzz_1[i] * pc_y[i];

        ta1_y_yzzz_yzzz_0[i] =
            ta1_y_zzz_zzz_0[i] * fe_0 - ta1_y_zzz_zzz_1[i] * fe_0 + ta_zzz_yzzz_1[i] + ta1_y_zzz_yzzz_0[i] * pa_y[i] - ta1_y_zzz_yzzz_1[i] * pc_y[i];

        ta1_y_yzzz_zzzz_0[i] = ta_zzz_zzzz_1[i] + ta1_y_zzz_zzzz_0[i] * pa_y[i] - ta1_y_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 435-450 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_y_zz_xxxx_0,   \
                             ta1_y_zz_xxxx_1,   \
                             ta1_y_zz_xxxy_0,   \
                             ta1_y_zz_xxxy_1,   \
                             ta1_y_zz_xxxz_0,   \
                             ta1_y_zz_xxxz_1,   \
                             ta1_y_zz_xxyy_0,   \
                             ta1_y_zz_xxyy_1,   \
                             ta1_y_zz_xxyz_0,   \
                             ta1_y_zz_xxyz_1,   \
                             ta1_y_zz_xxzz_0,   \
                             ta1_y_zz_xxzz_1,   \
                             ta1_y_zz_xyyy_0,   \
                             ta1_y_zz_xyyy_1,   \
                             ta1_y_zz_xyyz_0,   \
                             ta1_y_zz_xyyz_1,   \
                             ta1_y_zz_xyzz_0,   \
                             ta1_y_zz_xyzz_1,   \
                             ta1_y_zz_xzzz_0,   \
                             ta1_y_zz_xzzz_1,   \
                             ta1_y_zz_yyyy_0,   \
                             ta1_y_zz_yyyy_1,   \
                             ta1_y_zz_yyyz_0,   \
                             ta1_y_zz_yyyz_1,   \
                             ta1_y_zz_yyzz_0,   \
                             ta1_y_zz_yyzz_1,   \
                             ta1_y_zz_yzzz_0,   \
                             ta1_y_zz_yzzz_1,   \
                             ta1_y_zz_zzzz_0,   \
                             ta1_y_zz_zzzz_1,   \
                             ta1_y_zzz_xxx_0,   \
                             ta1_y_zzz_xxx_1,   \
                             ta1_y_zzz_xxxx_0,  \
                             ta1_y_zzz_xxxx_1,  \
                             ta1_y_zzz_xxxy_0,  \
                             ta1_y_zzz_xxxy_1,  \
                             ta1_y_zzz_xxxz_0,  \
                             ta1_y_zzz_xxxz_1,  \
                             ta1_y_zzz_xxy_0,   \
                             ta1_y_zzz_xxy_1,   \
                             ta1_y_zzz_xxyy_0,  \
                             ta1_y_zzz_xxyy_1,  \
                             ta1_y_zzz_xxyz_0,  \
                             ta1_y_zzz_xxyz_1,  \
                             ta1_y_zzz_xxz_0,   \
                             ta1_y_zzz_xxz_1,   \
                             ta1_y_zzz_xxzz_0,  \
                             ta1_y_zzz_xxzz_1,  \
                             ta1_y_zzz_xyy_0,   \
                             ta1_y_zzz_xyy_1,   \
                             ta1_y_zzz_xyyy_0,  \
                             ta1_y_zzz_xyyy_1,  \
                             ta1_y_zzz_xyyz_0,  \
                             ta1_y_zzz_xyyz_1,  \
                             ta1_y_zzz_xyz_0,   \
                             ta1_y_zzz_xyz_1,   \
                             ta1_y_zzz_xyzz_0,  \
                             ta1_y_zzz_xyzz_1,  \
                             ta1_y_zzz_xzz_0,   \
                             ta1_y_zzz_xzz_1,   \
                             ta1_y_zzz_xzzz_0,  \
                             ta1_y_zzz_xzzz_1,  \
                             ta1_y_zzz_yyy_0,   \
                             ta1_y_zzz_yyy_1,   \
                             ta1_y_zzz_yyyy_0,  \
                             ta1_y_zzz_yyyy_1,  \
                             ta1_y_zzz_yyyz_0,  \
                             ta1_y_zzz_yyyz_1,  \
                             ta1_y_zzz_yyz_0,   \
                             ta1_y_zzz_yyz_1,   \
                             ta1_y_zzz_yyzz_0,  \
                             ta1_y_zzz_yyzz_1,  \
                             ta1_y_zzz_yzz_0,   \
                             ta1_y_zzz_yzz_1,   \
                             ta1_y_zzz_yzzz_0,  \
                             ta1_y_zzz_yzzz_1,  \
                             ta1_y_zzz_zzz_0,   \
                             ta1_y_zzz_zzz_1,   \
                             ta1_y_zzz_zzzz_0,  \
                             ta1_y_zzz_zzzz_1,  \
                             ta1_y_zzzz_xxxx_0, \
                             ta1_y_zzzz_xxxy_0, \
                             ta1_y_zzzz_xxxz_0, \
                             ta1_y_zzzz_xxyy_0, \
                             ta1_y_zzzz_xxyz_0, \
                             ta1_y_zzzz_xxzz_0, \
                             ta1_y_zzzz_xyyy_0, \
                             ta1_y_zzzz_xyyz_0, \
                             ta1_y_zzzz_xyzz_0, \
                             ta1_y_zzzz_xzzz_0, \
                             ta1_y_zzzz_yyyy_0, \
                             ta1_y_zzzz_yyyz_0, \
                             ta1_y_zzzz_yyzz_0, \
                             ta1_y_zzzz_yzzz_0, \
                             ta1_y_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzz_xxxx_0[i] =
            3.0 * ta1_y_zz_xxxx_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxx_1[i] * fe_0 + ta1_y_zzz_xxxx_0[i] * pa_z[i] - ta1_y_zzz_xxxx_1[i] * pc_z[i];

        ta1_y_zzzz_xxxy_0[i] =
            3.0 * ta1_y_zz_xxxy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxy_1[i] * fe_0 + ta1_y_zzz_xxxy_0[i] * pa_z[i] - ta1_y_zzz_xxxy_1[i] * pc_z[i];

        ta1_y_zzzz_xxxz_0[i] = 3.0 * ta1_y_zz_xxxz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxz_1[i] * fe_0 + ta1_y_zzz_xxx_0[i] * fe_0 -
                               ta1_y_zzz_xxx_1[i] * fe_0 + ta1_y_zzz_xxxz_0[i] * pa_z[i] - ta1_y_zzz_xxxz_1[i] * pc_z[i];

        ta1_y_zzzz_xxyy_0[i] =
            3.0 * ta1_y_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyy_1[i] * fe_0 + ta1_y_zzz_xxyy_0[i] * pa_z[i] - ta1_y_zzz_xxyy_1[i] * pc_z[i];

        ta1_y_zzzz_xxyz_0[i] = 3.0 * ta1_y_zz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyz_1[i] * fe_0 + ta1_y_zzz_xxy_0[i] * fe_0 -
                               ta1_y_zzz_xxy_1[i] * fe_0 + ta1_y_zzz_xxyz_0[i] * pa_z[i] - ta1_y_zzz_xxyz_1[i] * pc_z[i];

        ta1_y_zzzz_xxzz_0[i] = 3.0 * ta1_y_zz_xxzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxz_0[i] * fe_0 -
                               2.0 * ta1_y_zzz_xxz_1[i] * fe_0 + ta1_y_zzz_xxzz_0[i] * pa_z[i] - ta1_y_zzz_xxzz_1[i] * pc_z[i];

        ta1_y_zzzz_xyyy_0[i] =
            3.0 * ta1_y_zz_xyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyy_1[i] * fe_0 + ta1_y_zzz_xyyy_0[i] * pa_z[i] - ta1_y_zzz_xyyy_1[i] * pc_z[i];

        ta1_y_zzzz_xyyz_0[i] = 3.0 * ta1_y_zz_xyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyz_1[i] * fe_0 + ta1_y_zzz_xyy_0[i] * fe_0 -
                               ta1_y_zzz_xyy_1[i] * fe_0 + ta1_y_zzz_xyyz_0[i] * pa_z[i] - ta1_y_zzz_xyyz_1[i] * pc_z[i];

        ta1_y_zzzz_xyzz_0[i] = 3.0 * ta1_y_zz_xyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_y_zzz_xyz_1[i] * fe_0 + ta1_y_zzz_xyzz_0[i] * pa_z[i] - ta1_y_zzz_xyzz_1[i] * pc_z[i];

        ta1_y_zzzz_xzzz_0[i] = 3.0 * ta1_y_zz_xzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_xzz_0[i] * fe_0 -
                               3.0 * ta1_y_zzz_xzz_1[i] * fe_0 + ta1_y_zzz_xzzz_0[i] * pa_z[i] - ta1_y_zzz_xzzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyyy_0[i] =
            3.0 * ta1_y_zz_yyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyy_1[i] * fe_0 + ta1_y_zzz_yyyy_0[i] * pa_z[i] - ta1_y_zzz_yyyy_1[i] * pc_z[i];

        ta1_y_zzzz_yyyz_0[i] = 3.0 * ta1_y_zz_yyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyz_1[i] * fe_0 + ta1_y_zzz_yyy_0[i] * fe_0 -
                               ta1_y_zzz_yyy_1[i] * fe_0 + ta1_y_zzz_yyyz_0[i] * pa_z[i] - ta1_y_zzz_yyyz_1[i] * pc_z[i];

        ta1_y_zzzz_yyzz_0[i] = 3.0 * ta1_y_zz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yyz_0[i] * fe_0 -
                               2.0 * ta1_y_zzz_yyz_1[i] * fe_0 + ta1_y_zzz_yyzz_0[i] * pa_z[i] - ta1_y_zzz_yyzz_1[i] * pc_z[i];

        ta1_y_zzzz_yzzz_0[i] = 3.0 * ta1_y_zz_yzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_yzz_0[i] * fe_0 -
                               3.0 * ta1_y_zzz_yzz_1[i] * fe_0 + ta1_y_zzz_yzzz_0[i] * pa_z[i] - ta1_y_zzz_yzzz_1[i] * pc_z[i];

        ta1_y_zzzz_zzzz_0[i] = 3.0 * ta1_y_zz_zzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_zzzz_1[i] * fe_0 + 4.0 * ta1_y_zzz_zzz_0[i] * fe_0 -
                               4.0 * ta1_y_zzz_zzz_1[i] * fe_0 + ta1_y_zzz_zzzz_0[i] * pa_z[i] - ta1_y_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 450-465 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xx_xxxx_0,   \
                             ta1_z_xx_xxxx_1,   \
                             ta1_z_xx_xxxy_0,   \
                             ta1_z_xx_xxxy_1,   \
                             ta1_z_xx_xxxz_0,   \
                             ta1_z_xx_xxxz_1,   \
                             ta1_z_xx_xxyy_0,   \
                             ta1_z_xx_xxyy_1,   \
                             ta1_z_xx_xxyz_0,   \
                             ta1_z_xx_xxyz_1,   \
                             ta1_z_xx_xxzz_0,   \
                             ta1_z_xx_xxzz_1,   \
                             ta1_z_xx_xyyy_0,   \
                             ta1_z_xx_xyyy_1,   \
                             ta1_z_xx_xyyz_0,   \
                             ta1_z_xx_xyyz_1,   \
                             ta1_z_xx_xyzz_0,   \
                             ta1_z_xx_xyzz_1,   \
                             ta1_z_xx_xzzz_0,   \
                             ta1_z_xx_xzzz_1,   \
                             ta1_z_xx_yyyy_0,   \
                             ta1_z_xx_yyyy_1,   \
                             ta1_z_xx_yyyz_0,   \
                             ta1_z_xx_yyyz_1,   \
                             ta1_z_xx_yyzz_0,   \
                             ta1_z_xx_yyzz_1,   \
                             ta1_z_xx_yzzz_0,   \
                             ta1_z_xx_yzzz_1,   \
                             ta1_z_xx_zzzz_0,   \
                             ta1_z_xx_zzzz_1,   \
                             ta1_z_xxx_xxx_0,   \
                             ta1_z_xxx_xxx_1,   \
                             ta1_z_xxx_xxxx_0,  \
                             ta1_z_xxx_xxxx_1,  \
                             ta1_z_xxx_xxxy_0,  \
                             ta1_z_xxx_xxxy_1,  \
                             ta1_z_xxx_xxxz_0,  \
                             ta1_z_xxx_xxxz_1,  \
                             ta1_z_xxx_xxy_0,   \
                             ta1_z_xxx_xxy_1,   \
                             ta1_z_xxx_xxyy_0,  \
                             ta1_z_xxx_xxyy_1,  \
                             ta1_z_xxx_xxyz_0,  \
                             ta1_z_xxx_xxyz_1,  \
                             ta1_z_xxx_xxz_0,   \
                             ta1_z_xxx_xxz_1,   \
                             ta1_z_xxx_xxzz_0,  \
                             ta1_z_xxx_xxzz_1,  \
                             ta1_z_xxx_xyy_0,   \
                             ta1_z_xxx_xyy_1,   \
                             ta1_z_xxx_xyyy_0,  \
                             ta1_z_xxx_xyyy_1,  \
                             ta1_z_xxx_xyyz_0,  \
                             ta1_z_xxx_xyyz_1,  \
                             ta1_z_xxx_xyz_0,   \
                             ta1_z_xxx_xyz_1,   \
                             ta1_z_xxx_xyzz_0,  \
                             ta1_z_xxx_xyzz_1,  \
                             ta1_z_xxx_xzz_0,   \
                             ta1_z_xxx_xzz_1,   \
                             ta1_z_xxx_xzzz_0,  \
                             ta1_z_xxx_xzzz_1,  \
                             ta1_z_xxx_yyy_0,   \
                             ta1_z_xxx_yyy_1,   \
                             ta1_z_xxx_yyyy_0,  \
                             ta1_z_xxx_yyyy_1,  \
                             ta1_z_xxx_yyyz_0,  \
                             ta1_z_xxx_yyyz_1,  \
                             ta1_z_xxx_yyz_0,   \
                             ta1_z_xxx_yyz_1,   \
                             ta1_z_xxx_yyzz_0,  \
                             ta1_z_xxx_yyzz_1,  \
                             ta1_z_xxx_yzz_0,   \
                             ta1_z_xxx_yzz_1,   \
                             ta1_z_xxx_yzzz_0,  \
                             ta1_z_xxx_yzzz_1,  \
                             ta1_z_xxx_zzz_0,   \
                             ta1_z_xxx_zzz_1,   \
                             ta1_z_xxx_zzzz_0,  \
                             ta1_z_xxx_zzzz_1,  \
                             ta1_z_xxxx_xxxx_0, \
                             ta1_z_xxxx_xxxy_0, \
                             ta1_z_xxxx_xxxz_0, \
                             ta1_z_xxxx_xxyy_0, \
                             ta1_z_xxxx_xxyz_0, \
                             ta1_z_xxxx_xxzz_0, \
                             ta1_z_xxxx_xyyy_0, \
                             ta1_z_xxxx_xyyz_0, \
                             ta1_z_xxxx_xyzz_0, \
                             ta1_z_xxxx_xzzz_0, \
                             ta1_z_xxxx_yyyy_0, \
                             ta1_z_xxxx_yyyz_0, \
                             ta1_z_xxxx_yyzz_0, \
                             ta1_z_xxxx_yzzz_0, \
                             ta1_z_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxx_xxxx_0[i] = 3.0 * ta1_z_xx_xxxx_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxx_1[i] * fe_0 + 4.0 * ta1_z_xxx_xxx_0[i] * fe_0 -
                               4.0 * ta1_z_xxx_xxx_1[i] * fe_0 + ta1_z_xxx_xxxx_0[i] * pa_x[i] - ta1_z_xxx_xxxx_1[i] * pc_x[i];

        ta1_z_xxxx_xxxy_0[i] = 3.0 * ta1_z_xx_xxxy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxy_0[i] * fe_0 -
                               3.0 * ta1_z_xxx_xxy_1[i] * fe_0 + ta1_z_xxx_xxxy_0[i] * pa_x[i] - ta1_z_xxx_xxxy_1[i] * pc_x[i];

        ta1_z_xxxx_xxxz_0[i] = 3.0 * ta1_z_xx_xxxz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxz_0[i] * fe_0 -
                               3.0 * ta1_z_xxx_xxz_1[i] * fe_0 + ta1_z_xxx_xxxz_0[i] * pa_x[i] - ta1_z_xxx_xxxz_1[i] * pc_x[i];

        ta1_z_xxxx_xxyy_0[i] = 3.0 * ta1_z_xx_xxyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyy_0[i] * fe_0 -
                               2.0 * ta1_z_xxx_xyy_1[i] * fe_0 + ta1_z_xxx_xxyy_0[i] * pa_x[i] - ta1_z_xxx_xxyy_1[i] * pc_x[i];

        ta1_z_xxxx_xxyz_0[i] = 3.0 * ta1_z_xx_xxyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyz_0[i] * fe_0 -
                               2.0 * ta1_z_xxx_xyz_1[i] * fe_0 + ta1_z_xxx_xxyz_0[i] * pa_x[i] - ta1_z_xxx_xxyz_1[i] * pc_x[i];

        ta1_z_xxxx_xxzz_0[i] = 3.0 * ta1_z_xx_xxzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xzz_0[i] * fe_0 -
                               2.0 * ta1_z_xxx_xzz_1[i] * fe_0 + ta1_z_xxx_xxzz_0[i] * pa_x[i] - ta1_z_xxx_xxzz_1[i] * pc_x[i];

        ta1_z_xxxx_xyyy_0[i] = 3.0 * ta1_z_xx_xyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyy_1[i] * fe_0 + ta1_z_xxx_yyy_0[i] * fe_0 -
                               ta1_z_xxx_yyy_1[i] * fe_0 + ta1_z_xxx_xyyy_0[i] * pa_x[i] - ta1_z_xxx_xyyy_1[i] * pc_x[i];

        ta1_z_xxxx_xyyz_0[i] = 3.0 * ta1_z_xx_xyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyz_1[i] * fe_0 + ta1_z_xxx_yyz_0[i] * fe_0 -
                               ta1_z_xxx_yyz_1[i] * fe_0 + ta1_z_xxx_xyyz_0[i] * pa_x[i] - ta1_z_xxx_xyyz_1[i] * pc_x[i];

        ta1_z_xxxx_xyzz_0[i] = 3.0 * ta1_z_xx_xyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyzz_1[i] * fe_0 + ta1_z_xxx_yzz_0[i] * fe_0 -
                               ta1_z_xxx_yzz_1[i] * fe_0 + ta1_z_xxx_xyzz_0[i] * pa_x[i] - ta1_z_xxx_xyzz_1[i] * pc_x[i];

        ta1_z_xxxx_xzzz_0[i] = 3.0 * ta1_z_xx_xzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xzzz_1[i] * fe_0 + ta1_z_xxx_zzz_0[i] * fe_0 -
                               ta1_z_xxx_zzz_1[i] * fe_0 + ta1_z_xxx_xzzz_0[i] * pa_x[i] - ta1_z_xxx_xzzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyyy_0[i] =
            3.0 * ta1_z_xx_yyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyy_1[i] * fe_0 + ta1_z_xxx_yyyy_0[i] * pa_x[i] - ta1_z_xxx_yyyy_1[i] * pc_x[i];

        ta1_z_xxxx_yyyz_0[i] =
            3.0 * ta1_z_xx_yyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyz_1[i] * fe_0 + ta1_z_xxx_yyyz_0[i] * pa_x[i] - ta1_z_xxx_yyyz_1[i] * pc_x[i];

        ta1_z_xxxx_yyzz_0[i] =
            3.0 * ta1_z_xx_yyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyzz_1[i] * fe_0 + ta1_z_xxx_yyzz_0[i] * pa_x[i] - ta1_z_xxx_yyzz_1[i] * pc_x[i];

        ta1_z_xxxx_yzzz_0[i] =
            3.0 * ta1_z_xx_yzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yzzz_1[i] * fe_0 + ta1_z_xxx_yzzz_0[i] * pa_x[i] - ta1_z_xxx_yzzz_1[i] * pc_x[i];

        ta1_z_xxxx_zzzz_0[i] =
            3.0 * ta1_z_xx_zzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_zzzz_1[i] * fe_0 + ta1_z_xxx_zzzz_0[i] * pa_x[i] - ta1_z_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 465-480 components of targeted buffer : GG

    auto ta1_z_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 465);

    auto ta1_z_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 466);

    auto ta1_z_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 467);

    auto ta1_z_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 468);

    auto ta1_z_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 469);

    auto ta1_z_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 470);

    auto ta1_z_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 471);

    auto ta1_z_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 472);

    auto ta1_z_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 473);

    auto ta1_z_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 474);

    auto ta1_z_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 475);

    auto ta1_z_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 476);

    auto ta1_z_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 477);

    auto ta1_z_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 478);

    auto ta1_z_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 479);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxx_xxx_0,   \
                             ta1_z_xxx_xxx_1,   \
                             ta1_z_xxx_xxxx_0,  \
                             ta1_z_xxx_xxxx_1,  \
                             ta1_z_xxx_xxxy_0,  \
                             ta1_z_xxx_xxxy_1,  \
                             ta1_z_xxx_xxxz_0,  \
                             ta1_z_xxx_xxxz_1,  \
                             ta1_z_xxx_xxy_0,   \
                             ta1_z_xxx_xxy_1,   \
                             ta1_z_xxx_xxyy_0,  \
                             ta1_z_xxx_xxyy_1,  \
                             ta1_z_xxx_xxyz_0,  \
                             ta1_z_xxx_xxyz_1,  \
                             ta1_z_xxx_xxz_0,   \
                             ta1_z_xxx_xxz_1,   \
                             ta1_z_xxx_xxzz_0,  \
                             ta1_z_xxx_xxzz_1,  \
                             ta1_z_xxx_xyy_0,   \
                             ta1_z_xxx_xyy_1,   \
                             ta1_z_xxx_xyyy_0,  \
                             ta1_z_xxx_xyyy_1,  \
                             ta1_z_xxx_xyyz_0,  \
                             ta1_z_xxx_xyyz_1,  \
                             ta1_z_xxx_xyz_0,   \
                             ta1_z_xxx_xyz_1,   \
                             ta1_z_xxx_xyzz_0,  \
                             ta1_z_xxx_xyzz_1,  \
                             ta1_z_xxx_xzz_0,   \
                             ta1_z_xxx_xzz_1,   \
                             ta1_z_xxx_xzzz_0,  \
                             ta1_z_xxx_xzzz_1,  \
                             ta1_z_xxx_zzzz_0,  \
                             ta1_z_xxx_zzzz_1,  \
                             ta1_z_xxxy_xxxx_0, \
                             ta1_z_xxxy_xxxy_0, \
                             ta1_z_xxxy_xxxz_0, \
                             ta1_z_xxxy_xxyy_0, \
                             ta1_z_xxxy_xxyz_0, \
                             ta1_z_xxxy_xxzz_0, \
                             ta1_z_xxxy_xyyy_0, \
                             ta1_z_xxxy_xyyz_0, \
                             ta1_z_xxxy_xyzz_0, \
                             ta1_z_xxxy_xzzz_0, \
                             ta1_z_xxxy_yyyy_0, \
                             ta1_z_xxxy_yyyz_0, \
                             ta1_z_xxxy_yyzz_0, \
                             ta1_z_xxxy_yzzz_0, \
                             ta1_z_xxxy_zzzz_0, \
                             ta1_z_xxy_yyyy_0,  \
                             ta1_z_xxy_yyyy_1,  \
                             ta1_z_xxy_yyyz_0,  \
                             ta1_z_xxy_yyyz_1,  \
                             ta1_z_xxy_yyzz_0,  \
                             ta1_z_xxy_yyzz_1,  \
                             ta1_z_xxy_yzzz_0,  \
                             ta1_z_xxy_yzzz_1,  \
                             ta1_z_xy_yyyy_0,   \
                             ta1_z_xy_yyyy_1,   \
                             ta1_z_xy_yyyz_0,   \
                             ta1_z_xy_yyyz_1,   \
                             ta1_z_xy_yyzz_0,   \
                             ta1_z_xy_yyzz_1,   \
                             ta1_z_xy_yzzz_0,   \
                             ta1_z_xy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxy_xxxx_0[i] = ta1_z_xxx_xxxx_0[i] * pa_y[i] - ta1_z_xxx_xxxx_1[i] * pc_y[i];

        ta1_z_xxxy_xxxy_0[i] = ta1_z_xxx_xxx_0[i] * fe_0 - ta1_z_xxx_xxx_1[i] * fe_0 + ta1_z_xxx_xxxy_0[i] * pa_y[i] - ta1_z_xxx_xxxy_1[i] * pc_y[i];

        ta1_z_xxxy_xxxz_0[i] = ta1_z_xxx_xxxz_0[i] * pa_y[i] - ta1_z_xxx_xxxz_1[i] * pc_y[i];

        ta1_z_xxxy_xxyy_0[i] =
            2.0 * ta1_z_xxx_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxy_1[i] * fe_0 + ta1_z_xxx_xxyy_0[i] * pa_y[i] - ta1_z_xxx_xxyy_1[i] * pc_y[i];

        ta1_z_xxxy_xxyz_0[i] = ta1_z_xxx_xxz_0[i] * fe_0 - ta1_z_xxx_xxz_1[i] * fe_0 + ta1_z_xxx_xxyz_0[i] * pa_y[i] - ta1_z_xxx_xxyz_1[i] * pc_y[i];

        ta1_z_xxxy_xxzz_0[i] = ta1_z_xxx_xxzz_0[i] * pa_y[i] - ta1_z_xxx_xxzz_1[i] * pc_y[i];

        ta1_z_xxxy_xyyy_0[i] =
            3.0 * ta1_z_xxx_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxx_xyy_1[i] * fe_0 + ta1_z_xxx_xyyy_0[i] * pa_y[i] - ta1_z_xxx_xyyy_1[i] * pc_y[i];

        ta1_z_xxxy_xyyz_0[i] =
            2.0 * ta1_z_xxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xyz_1[i] * fe_0 + ta1_z_xxx_xyyz_0[i] * pa_y[i] - ta1_z_xxx_xyyz_1[i] * pc_y[i];

        ta1_z_xxxy_xyzz_0[i] = ta1_z_xxx_xzz_0[i] * fe_0 - ta1_z_xxx_xzz_1[i] * fe_0 + ta1_z_xxx_xyzz_0[i] * pa_y[i] - ta1_z_xxx_xyzz_1[i] * pc_y[i];

        ta1_z_xxxy_xzzz_0[i] = ta1_z_xxx_xzzz_0[i] * pa_y[i] - ta1_z_xxx_xzzz_1[i] * pc_y[i];

        ta1_z_xxxy_yyyy_0[i] =
            2.0 * ta1_z_xy_yyyy_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyy_1[i] * fe_0 + ta1_z_xxy_yyyy_0[i] * pa_x[i] - ta1_z_xxy_yyyy_1[i] * pc_x[i];

        ta1_z_xxxy_yyyz_0[i] =
            2.0 * ta1_z_xy_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyz_1[i] * fe_0 + ta1_z_xxy_yyyz_0[i] * pa_x[i] - ta1_z_xxy_yyyz_1[i] * pc_x[i];

        ta1_z_xxxy_yyzz_0[i] =
            2.0 * ta1_z_xy_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyzz_1[i] * fe_0 + ta1_z_xxy_yyzz_0[i] * pa_x[i] - ta1_z_xxy_yyzz_1[i] * pc_x[i];

        ta1_z_xxxy_yzzz_0[i] =
            2.0 * ta1_z_xy_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yzzz_1[i] * fe_0 + ta1_z_xxy_yzzz_0[i] * pa_x[i] - ta1_z_xxy_yzzz_1[i] * pc_x[i];

        ta1_z_xxxy_zzzz_0[i] = ta1_z_xxx_zzzz_0[i] * pa_y[i] - ta1_z_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 480-495 components of targeted buffer : GG

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

    auto ta1_z_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 490);

    auto ta1_z_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 491);

    auto ta1_z_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 492);

    auto ta1_z_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 493);

    auto ta1_z_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 494);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xxx_xxx_0,   \
                             ta1_z_xxx_xxx_1,   \
                             ta1_z_xxx_xxxx_0,  \
                             ta1_z_xxx_xxxx_1,  \
                             ta1_z_xxx_xxxy_0,  \
                             ta1_z_xxx_xxxy_1,  \
                             ta1_z_xxx_xxxz_0,  \
                             ta1_z_xxx_xxxz_1,  \
                             ta1_z_xxx_xxy_0,   \
                             ta1_z_xxx_xxy_1,   \
                             ta1_z_xxx_xxyy_0,  \
                             ta1_z_xxx_xxyy_1,  \
                             ta1_z_xxx_xxyz_0,  \
                             ta1_z_xxx_xxyz_1,  \
                             ta1_z_xxx_xxz_0,   \
                             ta1_z_xxx_xxz_1,   \
                             ta1_z_xxx_xxzz_0,  \
                             ta1_z_xxx_xxzz_1,  \
                             ta1_z_xxx_xyy_0,   \
                             ta1_z_xxx_xyy_1,   \
                             ta1_z_xxx_xyyy_0,  \
                             ta1_z_xxx_xyyy_1,  \
                             ta1_z_xxx_xyyz_0,  \
                             ta1_z_xxx_xyyz_1,  \
                             ta1_z_xxx_xyz_0,   \
                             ta1_z_xxx_xyz_1,   \
                             ta1_z_xxx_xyzz_0,  \
                             ta1_z_xxx_xyzz_1,  \
                             ta1_z_xxx_xzz_0,   \
                             ta1_z_xxx_xzz_1,   \
                             ta1_z_xxx_xzzz_0,  \
                             ta1_z_xxx_xzzz_1,  \
                             ta1_z_xxx_yyyy_0,  \
                             ta1_z_xxx_yyyy_1,  \
                             ta1_z_xxxz_xxxx_0, \
                             ta1_z_xxxz_xxxy_0, \
                             ta1_z_xxxz_xxxz_0, \
                             ta1_z_xxxz_xxyy_0, \
                             ta1_z_xxxz_xxyz_0, \
                             ta1_z_xxxz_xxzz_0, \
                             ta1_z_xxxz_xyyy_0, \
                             ta1_z_xxxz_xyyz_0, \
                             ta1_z_xxxz_xyzz_0, \
                             ta1_z_xxxz_xzzz_0, \
                             ta1_z_xxxz_yyyy_0, \
                             ta1_z_xxxz_yyyz_0, \
                             ta1_z_xxxz_yyzz_0, \
                             ta1_z_xxxz_yzzz_0, \
                             ta1_z_xxxz_zzzz_0, \
                             ta1_z_xxz_yyyz_0,  \
                             ta1_z_xxz_yyyz_1,  \
                             ta1_z_xxz_yyzz_0,  \
                             ta1_z_xxz_yyzz_1,  \
                             ta1_z_xxz_yzzz_0,  \
                             ta1_z_xxz_yzzz_1,  \
                             ta1_z_xxz_zzzz_0,  \
                             ta1_z_xxz_zzzz_1,  \
                             ta1_z_xz_yyyz_0,   \
                             ta1_z_xz_yyyz_1,   \
                             ta1_z_xz_yyzz_0,   \
                             ta1_z_xz_yyzz_1,   \
                             ta1_z_xz_yzzz_0,   \
                             ta1_z_xz_yzzz_1,   \
                             ta1_z_xz_zzzz_0,   \
                             ta1_z_xz_zzzz_1,   \
                             ta_xxx_xxxx_1,     \
                             ta_xxx_xxxy_1,     \
                             ta_xxx_xxxz_1,     \
                             ta_xxx_xxyy_1,     \
                             ta_xxx_xxyz_1,     \
                             ta_xxx_xxzz_1,     \
                             ta_xxx_xyyy_1,     \
                             ta_xxx_xyyz_1,     \
                             ta_xxx_xyzz_1,     \
                             ta_xxx_xzzz_1,     \
                             ta_xxx_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxz_xxxx_0[i] = ta_xxx_xxxx_1[i] + ta1_z_xxx_xxxx_0[i] * pa_z[i] - ta1_z_xxx_xxxx_1[i] * pc_z[i];

        ta1_z_xxxz_xxxy_0[i] = ta_xxx_xxxy_1[i] + ta1_z_xxx_xxxy_0[i] * pa_z[i] - ta1_z_xxx_xxxy_1[i] * pc_z[i];

        ta1_z_xxxz_xxxz_0[i] =
            ta1_z_xxx_xxx_0[i] * fe_0 - ta1_z_xxx_xxx_1[i] * fe_0 + ta_xxx_xxxz_1[i] + ta1_z_xxx_xxxz_0[i] * pa_z[i] - ta1_z_xxx_xxxz_1[i] * pc_z[i];

        ta1_z_xxxz_xxyy_0[i] = ta_xxx_xxyy_1[i] + ta1_z_xxx_xxyy_0[i] * pa_z[i] - ta1_z_xxx_xxyy_1[i] * pc_z[i];

        ta1_z_xxxz_xxyz_0[i] =
            ta1_z_xxx_xxy_0[i] * fe_0 - ta1_z_xxx_xxy_1[i] * fe_0 + ta_xxx_xxyz_1[i] + ta1_z_xxx_xxyz_0[i] * pa_z[i] - ta1_z_xxx_xxyz_1[i] * pc_z[i];

        ta1_z_xxxz_xxzz_0[i] = 2.0 * ta1_z_xxx_xxz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxz_1[i] * fe_0 + ta_xxx_xxzz_1[i] + ta1_z_xxx_xxzz_0[i] * pa_z[i] -
                               ta1_z_xxx_xxzz_1[i] * pc_z[i];

        ta1_z_xxxz_xyyy_0[i] = ta_xxx_xyyy_1[i] + ta1_z_xxx_xyyy_0[i] * pa_z[i] - ta1_z_xxx_xyyy_1[i] * pc_z[i];

        ta1_z_xxxz_xyyz_0[i] =
            ta1_z_xxx_xyy_0[i] * fe_0 - ta1_z_xxx_xyy_1[i] * fe_0 + ta_xxx_xyyz_1[i] + ta1_z_xxx_xyyz_0[i] * pa_z[i] - ta1_z_xxx_xyyz_1[i] * pc_z[i];

        ta1_z_xxxz_xyzz_0[i] = 2.0 * ta1_z_xxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xyz_1[i] * fe_0 + ta_xxx_xyzz_1[i] + ta1_z_xxx_xyzz_0[i] * pa_z[i] -
                               ta1_z_xxx_xyzz_1[i] * pc_z[i];

        ta1_z_xxxz_xzzz_0[i] = 3.0 * ta1_z_xxx_xzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xzz_1[i] * fe_0 + ta_xxx_xzzz_1[i] + ta1_z_xxx_xzzz_0[i] * pa_z[i] -
                               ta1_z_xxx_xzzz_1[i] * pc_z[i];

        ta1_z_xxxz_yyyy_0[i] = ta_xxx_yyyy_1[i] + ta1_z_xxx_yyyy_0[i] * pa_z[i] - ta1_z_xxx_yyyy_1[i] * pc_z[i];

        ta1_z_xxxz_yyyz_0[i] =
            2.0 * ta1_z_xz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyyz_1[i] * fe_0 + ta1_z_xxz_yyyz_0[i] * pa_x[i] - ta1_z_xxz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxz_yyzz_0[i] =
            2.0 * ta1_z_xz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyzz_1[i] * fe_0 + ta1_z_xxz_yyzz_0[i] * pa_x[i] - ta1_z_xxz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxz_yzzz_0[i] =
            2.0 * ta1_z_xz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yzzz_1[i] * fe_0 + ta1_z_xxz_yzzz_0[i] * pa_x[i] - ta1_z_xxz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxz_zzzz_0[i] =
            2.0 * ta1_z_xz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_zzzz_1[i] * fe_0 + ta1_z_xxz_zzzz_0[i] * pa_x[i] - ta1_z_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 495-510 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xx_xxxx_0,   \
                             ta1_z_xx_xxxx_1,   \
                             ta1_z_xx_xxxz_0,   \
                             ta1_z_xx_xxxz_1,   \
                             ta1_z_xx_xxzz_0,   \
                             ta1_z_xx_xxzz_1,   \
                             ta1_z_xx_xzzz_0,   \
                             ta1_z_xx_xzzz_1,   \
                             ta1_z_xxy_xxxx_0,  \
                             ta1_z_xxy_xxxx_1,  \
                             ta1_z_xxy_xxxz_0,  \
                             ta1_z_xxy_xxxz_1,  \
                             ta1_z_xxy_xxzz_0,  \
                             ta1_z_xxy_xxzz_1,  \
                             ta1_z_xxy_xzzz_0,  \
                             ta1_z_xxy_xzzz_1,  \
                             ta1_z_xxyy_xxxx_0, \
                             ta1_z_xxyy_xxxy_0, \
                             ta1_z_xxyy_xxxz_0, \
                             ta1_z_xxyy_xxyy_0, \
                             ta1_z_xxyy_xxyz_0, \
                             ta1_z_xxyy_xxzz_0, \
                             ta1_z_xxyy_xyyy_0, \
                             ta1_z_xxyy_xyyz_0, \
                             ta1_z_xxyy_xyzz_0, \
                             ta1_z_xxyy_xzzz_0, \
                             ta1_z_xxyy_yyyy_0, \
                             ta1_z_xxyy_yyyz_0, \
                             ta1_z_xxyy_yyzz_0, \
                             ta1_z_xxyy_yzzz_0, \
                             ta1_z_xxyy_zzzz_0, \
                             ta1_z_xyy_xxxy_0,  \
                             ta1_z_xyy_xxxy_1,  \
                             ta1_z_xyy_xxy_0,   \
                             ta1_z_xyy_xxy_1,   \
                             ta1_z_xyy_xxyy_0,  \
                             ta1_z_xyy_xxyy_1,  \
                             ta1_z_xyy_xxyz_0,  \
                             ta1_z_xyy_xxyz_1,  \
                             ta1_z_xyy_xyy_0,   \
                             ta1_z_xyy_xyy_1,   \
                             ta1_z_xyy_xyyy_0,  \
                             ta1_z_xyy_xyyy_1,  \
                             ta1_z_xyy_xyyz_0,  \
                             ta1_z_xyy_xyyz_1,  \
                             ta1_z_xyy_xyz_0,   \
                             ta1_z_xyy_xyz_1,   \
                             ta1_z_xyy_xyzz_0,  \
                             ta1_z_xyy_xyzz_1,  \
                             ta1_z_xyy_yyy_0,   \
                             ta1_z_xyy_yyy_1,   \
                             ta1_z_xyy_yyyy_0,  \
                             ta1_z_xyy_yyyy_1,  \
                             ta1_z_xyy_yyyz_0,  \
                             ta1_z_xyy_yyyz_1,  \
                             ta1_z_xyy_yyz_0,   \
                             ta1_z_xyy_yyz_1,   \
                             ta1_z_xyy_yyzz_0,  \
                             ta1_z_xyy_yyzz_1,  \
                             ta1_z_xyy_yzz_0,   \
                             ta1_z_xyy_yzz_1,   \
                             ta1_z_xyy_yzzz_0,  \
                             ta1_z_xyy_yzzz_1,  \
                             ta1_z_xyy_zzzz_0,  \
                             ta1_z_xyy_zzzz_1,  \
                             ta1_z_yy_xxxy_0,   \
                             ta1_z_yy_xxxy_1,   \
                             ta1_z_yy_xxyy_0,   \
                             ta1_z_yy_xxyy_1,   \
                             ta1_z_yy_xxyz_0,   \
                             ta1_z_yy_xxyz_1,   \
                             ta1_z_yy_xyyy_0,   \
                             ta1_z_yy_xyyy_1,   \
                             ta1_z_yy_xyyz_0,   \
                             ta1_z_yy_xyyz_1,   \
                             ta1_z_yy_xyzz_0,   \
                             ta1_z_yy_xyzz_1,   \
                             ta1_z_yy_yyyy_0,   \
                             ta1_z_yy_yyyy_1,   \
                             ta1_z_yy_yyyz_0,   \
                             ta1_z_yy_yyyz_1,   \
                             ta1_z_yy_yyzz_0,   \
                             ta1_z_yy_yyzz_1,   \
                             ta1_z_yy_yzzz_0,   \
                             ta1_z_yy_yzzz_1,   \
                             ta1_z_yy_zzzz_0,   \
                             ta1_z_yy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyy_xxxx_0[i] = ta1_z_xx_xxxx_0[i] * fe_0 - ta1_z_xx_xxxx_1[i] * fe_0 + ta1_z_xxy_xxxx_0[i] * pa_y[i] - ta1_z_xxy_xxxx_1[i] * pc_y[i];

        ta1_z_xxyy_xxxy_0[i] = ta1_z_yy_xxxy_0[i] * fe_0 - ta1_z_yy_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xyy_xxy_0[i] * fe_0 -
                               3.0 * ta1_z_xyy_xxy_1[i] * fe_0 + ta1_z_xyy_xxxy_0[i] * pa_x[i] - ta1_z_xyy_xxxy_1[i] * pc_x[i];

        ta1_z_xxyy_xxxz_0[i] = ta1_z_xx_xxxz_0[i] * fe_0 - ta1_z_xx_xxxz_1[i] * fe_0 + ta1_z_xxy_xxxz_0[i] * pa_y[i] - ta1_z_xxy_xxxz_1[i] * pc_y[i];

        ta1_z_xxyy_xxyy_0[i] = ta1_z_yy_xxyy_0[i] * fe_0 - ta1_z_yy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyy_0[i] * fe_0 -
                               2.0 * ta1_z_xyy_xyy_1[i] * fe_0 + ta1_z_xyy_xxyy_0[i] * pa_x[i] - ta1_z_xyy_xxyy_1[i] * pc_x[i];

        ta1_z_xxyy_xxyz_0[i] = ta1_z_yy_xxyz_0[i] * fe_0 - ta1_z_yy_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyz_0[i] * fe_0 -
                               2.0 * ta1_z_xyy_xyz_1[i] * fe_0 + ta1_z_xyy_xxyz_0[i] * pa_x[i] - ta1_z_xyy_xxyz_1[i] * pc_x[i];

        ta1_z_xxyy_xxzz_0[i] = ta1_z_xx_xxzz_0[i] * fe_0 - ta1_z_xx_xxzz_1[i] * fe_0 + ta1_z_xxy_xxzz_0[i] * pa_y[i] - ta1_z_xxy_xxzz_1[i] * pc_y[i];

        ta1_z_xxyy_xyyy_0[i] = ta1_z_yy_xyyy_0[i] * fe_0 - ta1_z_yy_xyyy_1[i] * fe_0 + ta1_z_xyy_yyy_0[i] * fe_0 - ta1_z_xyy_yyy_1[i] * fe_0 +
                               ta1_z_xyy_xyyy_0[i] * pa_x[i] - ta1_z_xyy_xyyy_1[i] * pc_x[i];

        ta1_z_xxyy_xyyz_0[i] = ta1_z_yy_xyyz_0[i] * fe_0 - ta1_z_yy_xyyz_1[i] * fe_0 + ta1_z_xyy_yyz_0[i] * fe_0 - ta1_z_xyy_yyz_1[i] * fe_0 +
                               ta1_z_xyy_xyyz_0[i] * pa_x[i] - ta1_z_xyy_xyyz_1[i] * pc_x[i];

        ta1_z_xxyy_xyzz_0[i] = ta1_z_yy_xyzz_0[i] * fe_0 - ta1_z_yy_xyzz_1[i] * fe_0 + ta1_z_xyy_yzz_0[i] * fe_0 - ta1_z_xyy_yzz_1[i] * fe_0 +
                               ta1_z_xyy_xyzz_0[i] * pa_x[i] - ta1_z_xyy_xyzz_1[i] * pc_x[i];

        ta1_z_xxyy_xzzz_0[i] = ta1_z_xx_xzzz_0[i] * fe_0 - ta1_z_xx_xzzz_1[i] * fe_0 + ta1_z_xxy_xzzz_0[i] * pa_y[i] - ta1_z_xxy_xzzz_1[i] * pc_y[i];

        ta1_z_xxyy_yyyy_0[i] = ta1_z_yy_yyyy_0[i] * fe_0 - ta1_z_yy_yyyy_1[i] * fe_0 + ta1_z_xyy_yyyy_0[i] * pa_x[i] - ta1_z_xyy_yyyy_1[i] * pc_x[i];

        ta1_z_xxyy_yyyz_0[i] = ta1_z_yy_yyyz_0[i] * fe_0 - ta1_z_yy_yyyz_1[i] * fe_0 + ta1_z_xyy_yyyz_0[i] * pa_x[i] - ta1_z_xyy_yyyz_1[i] * pc_x[i];

        ta1_z_xxyy_yyzz_0[i] = ta1_z_yy_yyzz_0[i] * fe_0 - ta1_z_yy_yyzz_1[i] * fe_0 + ta1_z_xyy_yyzz_0[i] * pa_x[i] - ta1_z_xyy_yyzz_1[i] * pc_x[i];

        ta1_z_xxyy_yzzz_0[i] = ta1_z_yy_yzzz_0[i] * fe_0 - ta1_z_yy_yzzz_1[i] * fe_0 + ta1_z_xyy_yzzz_0[i] * pa_x[i] - ta1_z_xyy_yzzz_1[i] * pc_x[i];

        ta1_z_xxyy_zzzz_0[i] = ta1_z_yy_zzzz_0[i] * fe_0 - ta1_z_yy_zzzz_1[i] * fe_0 + ta1_z_xyy_zzzz_0[i] * pa_x[i] - ta1_z_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 510-525 components of targeted buffer : GG

    auto ta1_z_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 510);

    auto ta1_z_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 511);

    auto ta1_z_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 512);

    auto ta1_z_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 513);

    auto ta1_z_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 514);

    auto ta1_z_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 515);

    auto ta1_z_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 516);

    auto ta1_z_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 517);

    auto ta1_z_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 518);

    auto ta1_z_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 519);

    auto ta1_z_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 520);

    auto ta1_z_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 521);

    auto ta1_z_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 522);

    auto ta1_z_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 523);

    auto ta1_z_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 524);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_xxy_xxxy_0,  \
                             ta1_z_xxy_xxxy_1,  \
                             ta1_z_xxy_xxyy_0,  \
                             ta1_z_xxy_xxyy_1,  \
                             ta1_z_xxy_xyyy_0,  \
                             ta1_z_xxy_xyyy_1,  \
                             ta1_z_xxy_yyyy_0,  \
                             ta1_z_xxy_yyyy_1,  \
                             ta1_z_xxyz_xxxx_0, \
                             ta1_z_xxyz_xxxy_0, \
                             ta1_z_xxyz_xxxz_0, \
                             ta1_z_xxyz_xxyy_0, \
                             ta1_z_xxyz_xxyz_0, \
                             ta1_z_xxyz_xxzz_0, \
                             ta1_z_xxyz_xyyy_0, \
                             ta1_z_xxyz_xyyz_0, \
                             ta1_z_xxyz_xyzz_0, \
                             ta1_z_xxyz_xzzz_0, \
                             ta1_z_xxyz_yyyy_0, \
                             ta1_z_xxyz_yyyz_0, \
                             ta1_z_xxyz_yyzz_0, \
                             ta1_z_xxyz_yzzz_0, \
                             ta1_z_xxyz_zzzz_0, \
                             ta1_z_xxz_xxxx_0,  \
                             ta1_z_xxz_xxxx_1,  \
                             ta1_z_xxz_xxxz_0,  \
                             ta1_z_xxz_xxxz_1,  \
                             ta1_z_xxz_xxyz_0,  \
                             ta1_z_xxz_xxyz_1,  \
                             ta1_z_xxz_xxz_0,   \
                             ta1_z_xxz_xxz_1,   \
                             ta1_z_xxz_xxzz_0,  \
                             ta1_z_xxz_xxzz_1,  \
                             ta1_z_xxz_xyyz_0,  \
                             ta1_z_xxz_xyyz_1,  \
                             ta1_z_xxz_xyz_0,   \
                             ta1_z_xxz_xyz_1,   \
                             ta1_z_xxz_xyzz_0,  \
                             ta1_z_xxz_xyzz_1,  \
                             ta1_z_xxz_xzz_0,   \
                             ta1_z_xxz_xzz_1,   \
                             ta1_z_xxz_xzzz_0,  \
                             ta1_z_xxz_xzzz_1,  \
                             ta1_z_xxz_zzzz_0,  \
                             ta1_z_xxz_zzzz_1,  \
                             ta1_z_xyz_yyyz_0,  \
                             ta1_z_xyz_yyyz_1,  \
                             ta1_z_xyz_yyzz_0,  \
                             ta1_z_xyz_yyzz_1,  \
                             ta1_z_xyz_yzzz_0,  \
                             ta1_z_xyz_yzzz_1,  \
                             ta1_z_yz_yyyz_0,   \
                             ta1_z_yz_yyyz_1,   \
                             ta1_z_yz_yyzz_0,   \
                             ta1_z_yz_yyzz_1,   \
                             ta1_z_yz_yzzz_0,   \
                             ta1_z_yz_yzzz_1,   \
                             ta_xxy_xxxy_1,     \
                             ta_xxy_xxyy_1,     \
                             ta_xxy_xyyy_1,     \
                             ta_xxy_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyz_xxxx_0[i] = ta1_z_xxz_xxxx_0[i] * pa_y[i] - ta1_z_xxz_xxxx_1[i] * pc_y[i];

        ta1_z_xxyz_xxxy_0[i] = ta_xxy_xxxy_1[i] + ta1_z_xxy_xxxy_0[i] * pa_z[i] - ta1_z_xxy_xxxy_1[i] * pc_z[i];

        ta1_z_xxyz_xxxz_0[i] = ta1_z_xxz_xxxz_0[i] * pa_y[i] - ta1_z_xxz_xxxz_1[i] * pc_y[i];

        ta1_z_xxyz_xxyy_0[i] = ta_xxy_xxyy_1[i] + ta1_z_xxy_xxyy_0[i] * pa_z[i] - ta1_z_xxy_xxyy_1[i] * pc_z[i];

        ta1_z_xxyz_xxyz_0[i] = ta1_z_xxz_xxz_0[i] * fe_0 - ta1_z_xxz_xxz_1[i] * fe_0 + ta1_z_xxz_xxyz_0[i] * pa_y[i] - ta1_z_xxz_xxyz_1[i] * pc_y[i];

        ta1_z_xxyz_xxzz_0[i] = ta1_z_xxz_xxzz_0[i] * pa_y[i] - ta1_z_xxz_xxzz_1[i] * pc_y[i];

        ta1_z_xxyz_xyyy_0[i] = ta_xxy_xyyy_1[i] + ta1_z_xxy_xyyy_0[i] * pa_z[i] - ta1_z_xxy_xyyy_1[i] * pc_z[i];

        ta1_z_xxyz_xyyz_0[i] =
            2.0 * ta1_z_xxz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxz_xyz_1[i] * fe_0 + ta1_z_xxz_xyyz_0[i] * pa_y[i] - ta1_z_xxz_xyyz_1[i] * pc_y[i];

        ta1_z_xxyz_xyzz_0[i] = ta1_z_xxz_xzz_0[i] * fe_0 - ta1_z_xxz_xzz_1[i] * fe_0 + ta1_z_xxz_xyzz_0[i] * pa_y[i] - ta1_z_xxz_xyzz_1[i] * pc_y[i];

        ta1_z_xxyz_xzzz_0[i] = ta1_z_xxz_xzzz_0[i] * pa_y[i] - ta1_z_xxz_xzzz_1[i] * pc_y[i];

        ta1_z_xxyz_yyyy_0[i] = ta_xxy_yyyy_1[i] + ta1_z_xxy_yyyy_0[i] * pa_z[i] - ta1_z_xxy_yyyy_1[i] * pc_z[i];

        ta1_z_xxyz_yyyz_0[i] = ta1_z_yz_yyyz_0[i] * fe_0 - ta1_z_yz_yyyz_1[i] * fe_0 + ta1_z_xyz_yyyz_0[i] * pa_x[i] - ta1_z_xyz_yyyz_1[i] * pc_x[i];

        ta1_z_xxyz_yyzz_0[i] = ta1_z_yz_yyzz_0[i] * fe_0 - ta1_z_yz_yyzz_1[i] * fe_0 + ta1_z_xyz_yyzz_0[i] * pa_x[i] - ta1_z_xyz_yyzz_1[i] * pc_x[i];

        ta1_z_xxyz_yzzz_0[i] = ta1_z_yz_yzzz_0[i] * fe_0 - ta1_z_yz_yzzz_1[i] * fe_0 + ta1_z_xyz_yzzz_0[i] * pa_x[i] - ta1_z_xyz_yzzz_1[i] * pc_x[i];

        ta1_z_xxyz_zzzz_0[i] = ta1_z_xxz_zzzz_0[i] * pa_y[i] - ta1_z_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 525-540 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xx_xxxx_0,   \
                             ta1_z_xx_xxxx_1,   \
                             ta1_z_xx_xxxy_0,   \
                             ta1_z_xx_xxxy_1,   \
                             ta1_z_xx_xxyy_0,   \
                             ta1_z_xx_xxyy_1,   \
                             ta1_z_xx_xyyy_0,   \
                             ta1_z_xx_xyyy_1,   \
                             ta1_z_xxz_xxxx_0,  \
                             ta1_z_xxz_xxxx_1,  \
                             ta1_z_xxz_xxxy_0,  \
                             ta1_z_xxz_xxxy_1,  \
                             ta1_z_xxz_xxyy_0,  \
                             ta1_z_xxz_xxyy_1,  \
                             ta1_z_xxz_xyyy_0,  \
                             ta1_z_xxz_xyyy_1,  \
                             ta1_z_xxzz_xxxx_0, \
                             ta1_z_xxzz_xxxy_0, \
                             ta1_z_xxzz_xxxz_0, \
                             ta1_z_xxzz_xxyy_0, \
                             ta1_z_xxzz_xxyz_0, \
                             ta1_z_xxzz_xxzz_0, \
                             ta1_z_xxzz_xyyy_0, \
                             ta1_z_xxzz_xyyz_0, \
                             ta1_z_xxzz_xyzz_0, \
                             ta1_z_xxzz_xzzz_0, \
                             ta1_z_xxzz_yyyy_0, \
                             ta1_z_xxzz_yyyz_0, \
                             ta1_z_xxzz_yyzz_0, \
                             ta1_z_xxzz_yzzz_0, \
                             ta1_z_xxzz_zzzz_0, \
                             ta1_z_xzz_xxxz_0,  \
                             ta1_z_xzz_xxxz_1,  \
                             ta1_z_xzz_xxyz_0,  \
                             ta1_z_xzz_xxyz_1,  \
                             ta1_z_xzz_xxz_0,   \
                             ta1_z_xzz_xxz_1,   \
                             ta1_z_xzz_xxzz_0,  \
                             ta1_z_xzz_xxzz_1,  \
                             ta1_z_xzz_xyyz_0,  \
                             ta1_z_xzz_xyyz_1,  \
                             ta1_z_xzz_xyz_0,   \
                             ta1_z_xzz_xyz_1,   \
                             ta1_z_xzz_xyzz_0,  \
                             ta1_z_xzz_xyzz_1,  \
                             ta1_z_xzz_xzz_0,   \
                             ta1_z_xzz_xzz_1,   \
                             ta1_z_xzz_xzzz_0,  \
                             ta1_z_xzz_xzzz_1,  \
                             ta1_z_xzz_yyyy_0,  \
                             ta1_z_xzz_yyyy_1,  \
                             ta1_z_xzz_yyyz_0,  \
                             ta1_z_xzz_yyyz_1,  \
                             ta1_z_xzz_yyz_0,   \
                             ta1_z_xzz_yyz_1,   \
                             ta1_z_xzz_yyzz_0,  \
                             ta1_z_xzz_yyzz_1,  \
                             ta1_z_xzz_yzz_0,   \
                             ta1_z_xzz_yzz_1,   \
                             ta1_z_xzz_yzzz_0,  \
                             ta1_z_xzz_yzzz_1,  \
                             ta1_z_xzz_zzz_0,   \
                             ta1_z_xzz_zzz_1,   \
                             ta1_z_xzz_zzzz_0,  \
                             ta1_z_xzz_zzzz_1,  \
                             ta1_z_zz_xxxz_0,   \
                             ta1_z_zz_xxxz_1,   \
                             ta1_z_zz_xxyz_0,   \
                             ta1_z_zz_xxyz_1,   \
                             ta1_z_zz_xxzz_0,   \
                             ta1_z_zz_xxzz_1,   \
                             ta1_z_zz_xyyz_0,   \
                             ta1_z_zz_xyyz_1,   \
                             ta1_z_zz_xyzz_0,   \
                             ta1_z_zz_xyzz_1,   \
                             ta1_z_zz_xzzz_0,   \
                             ta1_z_zz_xzzz_1,   \
                             ta1_z_zz_yyyy_0,   \
                             ta1_z_zz_yyyy_1,   \
                             ta1_z_zz_yyyz_0,   \
                             ta1_z_zz_yyyz_1,   \
                             ta1_z_zz_yyzz_0,   \
                             ta1_z_zz_yyzz_1,   \
                             ta1_z_zz_yzzz_0,   \
                             ta1_z_zz_yzzz_1,   \
                             ta1_z_zz_zzzz_0,   \
                             ta1_z_zz_zzzz_1,   \
                             ta_xxz_xxxx_1,     \
                             ta_xxz_xxxy_1,     \
                             ta_xxz_xxyy_1,     \
                             ta_xxz_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzz_xxxx_0[i] =
            ta1_z_xx_xxxx_0[i] * fe_0 - ta1_z_xx_xxxx_1[i] * fe_0 + ta_xxz_xxxx_1[i] + ta1_z_xxz_xxxx_0[i] * pa_z[i] - ta1_z_xxz_xxxx_1[i] * pc_z[i];

        ta1_z_xxzz_xxxy_0[i] =
            ta1_z_xx_xxxy_0[i] * fe_0 - ta1_z_xx_xxxy_1[i] * fe_0 + ta_xxz_xxxy_1[i] + ta1_z_xxz_xxxy_0[i] * pa_z[i] - ta1_z_xxz_xxxy_1[i] * pc_z[i];

        ta1_z_xxzz_xxxz_0[i] = ta1_z_zz_xxxz_0[i] * fe_0 - ta1_z_zz_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xzz_xxz_0[i] * fe_0 -
                               3.0 * ta1_z_xzz_xxz_1[i] * fe_0 + ta1_z_xzz_xxxz_0[i] * pa_x[i] - ta1_z_xzz_xxxz_1[i] * pc_x[i];

        ta1_z_xxzz_xxyy_0[i] =
            ta1_z_xx_xxyy_0[i] * fe_0 - ta1_z_xx_xxyy_1[i] * fe_0 + ta_xxz_xxyy_1[i] + ta1_z_xxz_xxyy_0[i] * pa_z[i] - ta1_z_xxz_xxyy_1[i] * pc_z[i];

        ta1_z_xxzz_xxyz_0[i] = ta1_z_zz_xxyz_0[i] * fe_0 - ta1_z_zz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_z_xzz_xyz_1[i] * fe_0 + ta1_z_xzz_xxyz_0[i] * pa_x[i] - ta1_z_xzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxzz_xxzz_0[i] = ta1_z_zz_xxzz_0[i] * fe_0 - ta1_z_zz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xzz_0[i] * fe_0 -
                               2.0 * ta1_z_xzz_xzz_1[i] * fe_0 + ta1_z_xzz_xxzz_0[i] * pa_x[i] - ta1_z_xzz_xxzz_1[i] * pc_x[i];

        ta1_z_xxzz_xyyy_0[i] =
            ta1_z_xx_xyyy_0[i] * fe_0 - ta1_z_xx_xyyy_1[i] * fe_0 + ta_xxz_xyyy_1[i] + ta1_z_xxz_xyyy_0[i] * pa_z[i] - ta1_z_xxz_xyyy_1[i] * pc_z[i];

        ta1_z_xxzz_xyyz_0[i] = ta1_z_zz_xyyz_0[i] * fe_0 - ta1_z_zz_xyyz_1[i] * fe_0 + ta1_z_xzz_yyz_0[i] * fe_0 - ta1_z_xzz_yyz_1[i] * fe_0 +
                               ta1_z_xzz_xyyz_0[i] * pa_x[i] - ta1_z_xzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxzz_xyzz_0[i] = ta1_z_zz_xyzz_0[i] * fe_0 - ta1_z_zz_xyzz_1[i] * fe_0 + ta1_z_xzz_yzz_0[i] * fe_0 - ta1_z_xzz_yzz_1[i] * fe_0 +
                               ta1_z_xzz_xyzz_0[i] * pa_x[i] - ta1_z_xzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxzz_xzzz_0[i] = ta1_z_zz_xzzz_0[i] * fe_0 - ta1_z_zz_xzzz_1[i] * fe_0 + ta1_z_xzz_zzz_0[i] * fe_0 - ta1_z_xzz_zzz_1[i] * fe_0 +
                               ta1_z_xzz_xzzz_0[i] * pa_x[i] - ta1_z_xzz_xzzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyyy_0[i] = ta1_z_zz_yyyy_0[i] * fe_0 - ta1_z_zz_yyyy_1[i] * fe_0 + ta1_z_xzz_yyyy_0[i] * pa_x[i] - ta1_z_xzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxzz_yyyz_0[i] = ta1_z_zz_yyyz_0[i] * fe_0 - ta1_z_zz_yyyz_1[i] * fe_0 + ta1_z_xzz_yyyz_0[i] * pa_x[i] - ta1_z_xzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxzz_yyzz_0[i] = ta1_z_zz_yyzz_0[i] * fe_0 - ta1_z_zz_yyzz_1[i] * fe_0 + ta1_z_xzz_yyzz_0[i] * pa_x[i] - ta1_z_xzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxzz_yzzz_0[i] = ta1_z_zz_yzzz_0[i] * fe_0 - ta1_z_zz_yzzz_1[i] * fe_0 + ta1_z_xzz_yzzz_0[i] * pa_x[i] - ta1_z_xzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxzz_zzzz_0[i] = ta1_z_zz_zzzz_0[i] * fe_0 - ta1_z_zz_zzzz_1[i] * fe_0 + ta1_z_xzz_zzzz_0[i] * pa_x[i] - ta1_z_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 540-555 components of targeted buffer : GG

    auto ta1_z_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 540);

    auto ta1_z_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 541);

    auto ta1_z_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 542);

    auto ta1_z_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 543);

    auto ta1_z_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 544);

    auto ta1_z_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 545);

    auto ta1_z_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 546);

    auto ta1_z_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 547);

    auto ta1_z_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 548);

    auto ta1_z_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 549);

    auto ta1_z_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 550);

    auto ta1_z_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 551);

    auto ta1_z_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 552);

    auto ta1_z_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 553);

    auto ta1_z_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 554);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xyyy_xxxx_0, \
                             ta1_z_xyyy_xxxy_0, \
                             ta1_z_xyyy_xxxz_0, \
                             ta1_z_xyyy_xxyy_0, \
                             ta1_z_xyyy_xxyz_0, \
                             ta1_z_xyyy_xxzz_0, \
                             ta1_z_xyyy_xyyy_0, \
                             ta1_z_xyyy_xyyz_0, \
                             ta1_z_xyyy_xyzz_0, \
                             ta1_z_xyyy_xzzz_0, \
                             ta1_z_xyyy_yyyy_0, \
                             ta1_z_xyyy_yyyz_0, \
                             ta1_z_xyyy_yyzz_0, \
                             ta1_z_xyyy_yzzz_0, \
                             ta1_z_xyyy_zzzz_0, \
                             ta1_z_yyy_xxx_0,   \
                             ta1_z_yyy_xxx_1,   \
                             ta1_z_yyy_xxxx_0,  \
                             ta1_z_yyy_xxxx_1,  \
                             ta1_z_yyy_xxxy_0,  \
                             ta1_z_yyy_xxxy_1,  \
                             ta1_z_yyy_xxxz_0,  \
                             ta1_z_yyy_xxxz_1,  \
                             ta1_z_yyy_xxy_0,   \
                             ta1_z_yyy_xxy_1,   \
                             ta1_z_yyy_xxyy_0,  \
                             ta1_z_yyy_xxyy_1,  \
                             ta1_z_yyy_xxyz_0,  \
                             ta1_z_yyy_xxyz_1,  \
                             ta1_z_yyy_xxz_0,   \
                             ta1_z_yyy_xxz_1,   \
                             ta1_z_yyy_xxzz_0,  \
                             ta1_z_yyy_xxzz_1,  \
                             ta1_z_yyy_xyy_0,   \
                             ta1_z_yyy_xyy_1,   \
                             ta1_z_yyy_xyyy_0,  \
                             ta1_z_yyy_xyyy_1,  \
                             ta1_z_yyy_xyyz_0,  \
                             ta1_z_yyy_xyyz_1,  \
                             ta1_z_yyy_xyz_0,   \
                             ta1_z_yyy_xyz_1,   \
                             ta1_z_yyy_xyzz_0,  \
                             ta1_z_yyy_xyzz_1,  \
                             ta1_z_yyy_xzz_0,   \
                             ta1_z_yyy_xzz_1,   \
                             ta1_z_yyy_xzzz_0,  \
                             ta1_z_yyy_xzzz_1,  \
                             ta1_z_yyy_yyy_0,   \
                             ta1_z_yyy_yyy_1,   \
                             ta1_z_yyy_yyyy_0,  \
                             ta1_z_yyy_yyyy_1,  \
                             ta1_z_yyy_yyyz_0,  \
                             ta1_z_yyy_yyyz_1,  \
                             ta1_z_yyy_yyz_0,   \
                             ta1_z_yyy_yyz_1,   \
                             ta1_z_yyy_yyzz_0,  \
                             ta1_z_yyy_yyzz_1,  \
                             ta1_z_yyy_yzz_0,   \
                             ta1_z_yyy_yzz_1,   \
                             ta1_z_yyy_yzzz_0,  \
                             ta1_z_yyy_yzzz_1,  \
                             ta1_z_yyy_zzz_0,   \
                             ta1_z_yyy_zzz_1,   \
                             ta1_z_yyy_zzzz_0,  \
                             ta1_z_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyy_xxxx_0[i] =
            4.0 * ta1_z_yyy_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxx_1[i] * fe_0 + ta1_z_yyy_xxxx_0[i] * pa_x[i] - ta1_z_yyy_xxxx_1[i] * pc_x[i];

        ta1_z_xyyy_xxxy_0[i] =
            3.0 * ta1_z_yyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxy_1[i] * fe_0 + ta1_z_yyy_xxxy_0[i] * pa_x[i] - ta1_z_yyy_xxxy_1[i] * pc_x[i];

        ta1_z_xyyy_xxxz_0[i] =
            3.0 * ta1_z_yyy_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxz_1[i] * fe_0 + ta1_z_yyy_xxxz_0[i] * pa_x[i] - ta1_z_yyy_xxxz_1[i] * pc_x[i];

        ta1_z_xyyy_xxyy_0[i] =
            2.0 * ta1_z_yyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyy_1[i] * fe_0 + ta1_z_yyy_xxyy_0[i] * pa_x[i] - ta1_z_yyy_xxyy_1[i] * pc_x[i];

        ta1_z_xyyy_xxyz_0[i] =
            2.0 * ta1_z_yyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyz_1[i] * fe_0 + ta1_z_yyy_xxyz_0[i] * pa_x[i] - ta1_z_yyy_xxyz_1[i] * pc_x[i];

        ta1_z_xyyy_xxzz_0[i] =
            2.0 * ta1_z_yyy_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xzz_1[i] * fe_0 + ta1_z_yyy_xxzz_0[i] * pa_x[i] - ta1_z_yyy_xxzz_1[i] * pc_x[i];

        ta1_z_xyyy_xyyy_0[i] = ta1_z_yyy_yyy_0[i] * fe_0 - ta1_z_yyy_yyy_1[i] * fe_0 + ta1_z_yyy_xyyy_0[i] * pa_x[i] - ta1_z_yyy_xyyy_1[i] * pc_x[i];

        ta1_z_xyyy_xyyz_0[i] = ta1_z_yyy_yyz_0[i] * fe_0 - ta1_z_yyy_yyz_1[i] * fe_0 + ta1_z_yyy_xyyz_0[i] * pa_x[i] - ta1_z_yyy_xyyz_1[i] * pc_x[i];

        ta1_z_xyyy_xyzz_0[i] = ta1_z_yyy_yzz_0[i] * fe_0 - ta1_z_yyy_yzz_1[i] * fe_0 + ta1_z_yyy_xyzz_0[i] * pa_x[i] - ta1_z_yyy_xyzz_1[i] * pc_x[i];

        ta1_z_xyyy_xzzz_0[i] = ta1_z_yyy_zzz_0[i] * fe_0 - ta1_z_yyy_zzz_1[i] * fe_0 + ta1_z_yyy_xzzz_0[i] * pa_x[i] - ta1_z_yyy_xzzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyyy_0[i] = ta1_z_yyy_yyyy_0[i] * pa_x[i] - ta1_z_yyy_yyyy_1[i] * pc_x[i];

        ta1_z_xyyy_yyyz_0[i] = ta1_z_yyy_yyyz_0[i] * pa_x[i] - ta1_z_yyy_yyyz_1[i] * pc_x[i];

        ta1_z_xyyy_yyzz_0[i] = ta1_z_yyy_yyzz_0[i] * pa_x[i] - ta1_z_yyy_yyzz_1[i] * pc_x[i];

        ta1_z_xyyy_yzzz_0[i] = ta1_z_yyy_yzzz_0[i] * pa_x[i] - ta1_z_yyy_yzzz_1[i] * pc_x[i];

        ta1_z_xyyy_zzzz_0[i] = ta1_z_yyy_zzzz_0[i] * pa_x[i] - ta1_z_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 555-570 components of targeted buffer : GG

    auto ta1_z_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 555);

    auto ta1_z_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 556);

    auto ta1_z_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 557);

    auto ta1_z_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 558);

    auto ta1_z_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 559);

    auto ta1_z_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 560);

    auto ta1_z_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 561);

    auto ta1_z_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 562);

    auto ta1_z_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 563);

    auto ta1_z_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 564);

    auto ta1_z_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 565);

    auto ta1_z_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 566);

    auto ta1_z_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 567);

    auto ta1_z_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 568);

    auto ta1_z_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 569);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xyy_xxxx_0,  \
                             ta1_z_xyy_xxxx_1,  \
                             ta1_z_xyy_xxxy_0,  \
                             ta1_z_xyy_xxxy_1,  \
                             ta1_z_xyy_xxyy_0,  \
                             ta1_z_xyy_xxyy_1,  \
                             ta1_z_xyy_xyyy_0,  \
                             ta1_z_xyy_xyyy_1,  \
                             ta1_z_xyyz_xxxx_0, \
                             ta1_z_xyyz_xxxy_0, \
                             ta1_z_xyyz_xxxz_0, \
                             ta1_z_xyyz_xxyy_0, \
                             ta1_z_xyyz_xxyz_0, \
                             ta1_z_xyyz_xxzz_0, \
                             ta1_z_xyyz_xyyy_0, \
                             ta1_z_xyyz_xyyz_0, \
                             ta1_z_xyyz_xyzz_0, \
                             ta1_z_xyyz_xzzz_0, \
                             ta1_z_xyyz_yyyy_0, \
                             ta1_z_xyyz_yyyz_0, \
                             ta1_z_xyyz_yyzz_0, \
                             ta1_z_xyyz_yzzz_0, \
                             ta1_z_xyyz_zzzz_0, \
                             ta1_z_yyz_xxxz_0,  \
                             ta1_z_yyz_xxxz_1,  \
                             ta1_z_yyz_xxyz_0,  \
                             ta1_z_yyz_xxyz_1,  \
                             ta1_z_yyz_xxz_0,   \
                             ta1_z_yyz_xxz_1,   \
                             ta1_z_yyz_xxzz_0,  \
                             ta1_z_yyz_xxzz_1,  \
                             ta1_z_yyz_xyyz_0,  \
                             ta1_z_yyz_xyyz_1,  \
                             ta1_z_yyz_xyz_0,   \
                             ta1_z_yyz_xyz_1,   \
                             ta1_z_yyz_xyzz_0,  \
                             ta1_z_yyz_xyzz_1,  \
                             ta1_z_yyz_xzz_0,   \
                             ta1_z_yyz_xzz_1,   \
                             ta1_z_yyz_xzzz_0,  \
                             ta1_z_yyz_xzzz_1,  \
                             ta1_z_yyz_yyyy_0,  \
                             ta1_z_yyz_yyyy_1,  \
                             ta1_z_yyz_yyyz_0,  \
                             ta1_z_yyz_yyyz_1,  \
                             ta1_z_yyz_yyz_0,   \
                             ta1_z_yyz_yyz_1,   \
                             ta1_z_yyz_yyzz_0,  \
                             ta1_z_yyz_yyzz_1,  \
                             ta1_z_yyz_yzz_0,   \
                             ta1_z_yyz_yzz_1,   \
                             ta1_z_yyz_yzzz_0,  \
                             ta1_z_yyz_yzzz_1,  \
                             ta1_z_yyz_zzz_0,   \
                             ta1_z_yyz_zzz_1,   \
                             ta1_z_yyz_zzzz_0,  \
                             ta1_z_yyz_zzzz_1,  \
                             ta_xyy_xxxx_1,     \
                             ta_xyy_xxxy_1,     \
                             ta_xyy_xxyy_1,     \
                             ta_xyy_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyz_xxxx_0[i] = ta_xyy_xxxx_1[i] + ta1_z_xyy_xxxx_0[i] * pa_z[i] - ta1_z_xyy_xxxx_1[i] * pc_z[i];

        ta1_z_xyyz_xxxy_0[i] = ta_xyy_xxxy_1[i] + ta1_z_xyy_xxxy_0[i] * pa_z[i] - ta1_z_xyy_xxxy_1[i] * pc_z[i];

        ta1_z_xyyz_xxxz_0[i] =
            3.0 * ta1_z_yyz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxz_1[i] * fe_0 + ta1_z_yyz_xxxz_0[i] * pa_x[i] - ta1_z_yyz_xxxz_1[i] * pc_x[i];

        ta1_z_xyyz_xxyy_0[i] = ta_xyy_xxyy_1[i] + ta1_z_xyy_xxyy_0[i] * pa_z[i] - ta1_z_xyy_xxyy_1[i] * pc_z[i];

        ta1_z_xyyz_xxyz_0[i] =
            2.0 * ta1_z_yyz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyz_1[i] * fe_0 + ta1_z_yyz_xxyz_0[i] * pa_x[i] - ta1_z_yyz_xxyz_1[i] * pc_x[i];

        ta1_z_xyyz_xxzz_0[i] =
            2.0 * ta1_z_yyz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xzz_1[i] * fe_0 + ta1_z_yyz_xxzz_0[i] * pa_x[i] - ta1_z_yyz_xxzz_1[i] * pc_x[i];

        ta1_z_xyyz_xyyy_0[i] = ta_xyy_xyyy_1[i] + ta1_z_xyy_xyyy_0[i] * pa_z[i] - ta1_z_xyy_xyyy_1[i] * pc_z[i];

        ta1_z_xyyz_xyyz_0[i] = ta1_z_yyz_yyz_0[i] * fe_0 - ta1_z_yyz_yyz_1[i] * fe_0 + ta1_z_yyz_xyyz_0[i] * pa_x[i] - ta1_z_yyz_xyyz_1[i] * pc_x[i];

        ta1_z_xyyz_xyzz_0[i] = ta1_z_yyz_yzz_0[i] * fe_0 - ta1_z_yyz_yzz_1[i] * fe_0 + ta1_z_yyz_xyzz_0[i] * pa_x[i] - ta1_z_yyz_xyzz_1[i] * pc_x[i];

        ta1_z_xyyz_xzzz_0[i] = ta1_z_yyz_zzz_0[i] * fe_0 - ta1_z_yyz_zzz_1[i] * fe_0 + ta1_z_yyz_xzzz_0[i] * pa_x[i] - ta1_z_yyz_xzzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyyy_0[i] = ta1_z_yyz_yyyy_0[i] * pa_x[i] - ta1_z_yyz_yyyy_1[i] * pc_x[i];

        ta1_z_xyyz_yyyz_0[i] = ta1_z_yyz_yyyz_0[i] * pa_x[i] - ta1_z_yyz_yyyz_1[i] * pc_x[i];

        ta1_z_xyyz_yyzz_0[i] = ta1_z_yyz_yyzz_0[i] * pa_x[i] - ta1_z_yyz_yyzz_1[i] * pc_x[i];

        ta1_z_xyyz_yzzz_0[i] = ta1_z_yyz_yzzz_0[i] * pa_x[i] - ta1_z_yyz_yzzz_1[i] * pc_x[i];

        ta1_z_xyyz_zzzz_0[i] = ta1_z_yyz_zzzz_0[i] * pa_x[i] - ta1_z_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 570-585 components of targeted buffer : GG

    auto ta1_z_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 570);

    auto ta1_z_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 571);

    auto ta1_z_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 572);

    auto ta1_z_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 573);

    auto ta1_z_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 574);

    auto ta1_z_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 575);

    auto ta1_z_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 576);

    auto ta1_z_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 577);

    auto ta1_z_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 578);

    auto ta1_z_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 579);

    auto ta1_z_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 580);

    auto ta1_z_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 581);

    auto ta1_z_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 582);

    auto ta1_z_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 583);

    auto ta1_z_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 584);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xyzz_xxxx_0, \
                             ta1_z_xyzz_xxxy_0, \
                             ta1_z_xyzz_xxxz_0, \
                             ta1_z_xyzz_xxyy_0, \
                             ta1_z_xyzz_xxyz_0, \
                             ta1_z_xyzz_xxzz_0, \
                             ta1_z_xyzz_xyyy_0, \
                             ta1_z_xyzz_xyyz_0, \
                             ta1_z_xyzz_xyzz_0, \
                             ta1_z_xyzz_xzzz_0, \
                             ta1_z_xyzz_yyyy_0, \
                             ta1_z_xyzz_yyyz_0, \
                             ta1_z_xyzz_yyzz_0, \
                             ta1_z_xyzz_yzzz_0, \
                             ta1_z_xyzz_zzzz_0, \
                             ta1_z_xzz_xxxx_0,  \
                             ta1_z_xzz_xxxx_1,  \
                             ta1_z_xzz_xxxz_0,  \
                             ta1_z_xzz_xxxz_1,  \
                             ta1_z_xzz_xxzz_0,  \
                             ta1_z_xzz_xxzz_1,  \
                             ta1_z_xzz_xzzz_0,  \
                             ta1_z_xzz_xzzz_1,  \
                             ta1_z_yzz_xxxy_0,  \
                             ta1_z_yzz_xxxy_1,  \
                             ta1_z_yzz_xxy_0,   \
                             ta1_z_yzz_xxy_1,   \
                             ta1_z_yzz_xxyy_0,  \
                             ta1_z_yzz_xxyy_1,  \
                             ta1_z_yzz_xxyz_0,  \
                             ta1_z_yzz_xxyz_1,  \
                             ta1_z_yzz_xyy_0,   \
                             ta1_z_yzz_xyy_1,   \
                             ta1_z_yzz_xyyy_0,  \
                             ta1_z_yzz_xyyy_1,  \
                             ta1_z_yzz_xyyz_0,  \
                             ta1_z_yzz_xyyz_1,  \
                             ta1_z_yzz_xyz_0,   \
                             ta1_z_yzz_xyz_1,   \
                             ta1_z_yzz_xyzz_0,  \
                             ta1_z_yzz_xyzz_1,  \
                             ta1_z_yzz_yyy_0,   \
                             ta1_z_yzz_yyy_1,   \
                             ta1_z_yzz_yyyy_0,  \
                             ta1_z_yzz_yyyy_1,  \
                             ta1_z_yzz_yyyz_0,  \
                             ta1_z_yzz_yyyz_1,  \
                             ta1_z_yzz_yyz_0,   \
                             ta1_z_yzz_yyz_1,   \
                             ta1_z_yzz_yyzz_0,  \
                             ta1_z_yzz_yyzz_1,  \
                             ta1_z_yzz_yzz_0,   \
                             ta1_z_yzz_yzz_1,   \
                             ta1_z_yzz_yzzz_0,  \
                             ta1_z_yzz_yzzz_1,  \
                             ta1_z_yzz_zzzz_0,  \
                             ta1_z_yzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzz_xxxx_0[i] = ta1_z_xzz_xxxx_0[i] * pa_y[i] - ta1_z_xzz_xxxx_1[i] * pc_y[i];

        ta1_z_xyzz_xxxy_0[i] =
            3.0 * ta1_z_yzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yzz_xxy_1[i] * fe_0 + ta1_z_yzz_xxxy_0[i] * pa_x[i] - ta1_z_yzz_xxxy_1[i] * pc_x[i];

        ta1_z_xyzz_xxxz_0[i] = ta1_z_xzz_xxxz_0[i] * pa_y[i] - ta1_z_xzz_xxxz_1[i] * pc_y[i];

        ta1_z_xyzz_xxyy_0[i] =
            2.0 * ta1_z_yzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyy_1[i] * fe_0 + ta1_z_yzz_xxyy_0[i] * pa_x[i] - ta1_z_yzz_xxyy_1[i] * pc_x[i];

        ta1_z_xyzz_xxyz_0[i] =
            2.0 * ta1_z_yzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyz_1[i] * fe_0 + ta1_z_yzz_xxyz_0[i] * pa_x[i] - ta1_z_yzz_xxyz_1[i] * pc_x[i];

        ta1_z_xyzz_xxzz_0[i] = ta1_z_xzz_xxzz_0[i] * pa_y[i] - ta1_z_xzz_xxzz_1[i] * pc_y[i];

        ta1_z_xyzz_xyyy_0[i] = ta1_z_yzz_yyy_0[i] * fe_0 - ta1_z_yzz_yyy_1[i] * fe_0 + ta1_z_yzz_xyyy_0[i] * pa_x[i] - ta1_z_yzz_xyyy_1[i] * pc_x[i];

        ta1_z_xyzz_xyyz_0[i] = ta1_z_yzz_yyz_0[i] * fe_0 - ta1_z_yzz_yyz_1[i] * fe_0 + ta1_z_yzz_xyyz_0[i] * pa_x[i] - ta1_z_yzz_xyyz_1[i] * pc_x[i];

        ta1_z_xyzz_xyzz_0[i] = ta1_z_yzz_yzz_0[i] * fe_0 - ta1_z_yzz_yzz_1[i] * fe_0 + ta1_z_yzz_xyzz_0[i] * pa_x[i] - ta1_z_yzz_xyzz_1[i] * pc_x[i];

        ta1_z_xyzz_xzzz_0[i] = ta1_z_xzz_xzzz_0[i] * pa_y[i] - ta1_z_xzz_xzzz_1[i] * pc_y[i];

        ta1_z_xyzz_yyyy_0[i] = ta1_z_yzz_yyyy_0[i] * pa_x[i] - ta1_z_yzz_yyyy_1[i] * pc_x[i];

        ta1_z_xyzz_yyyz_0[i] = ta1_z_yzz_yyyz_0[i] * pa_x[i] - ta1_z_yzz_yyyz_1[i] * pc_x[i];

        ta1_z_xyzz_yyzz_0[i] = ta1_z_yzz_yyzz_0[i] * pa_x[i] - ta1_z_yzz_yyzz_1[i] * pc_x[i];

        ta1_z_xyzz_yzzz_0[i] = ta1_z_yzz_yzzz_0[i] * pa_x[i] - ta1_z_yzz_yzzz_1[i] * pc_x[i];

        ta1_z_xyzz_zzzz_0[i] = ta1_z_yzz_zzzz_0[i] * pa_x[i] - ta1_z_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 585-600 components of targeted buffer : GG

    auto ta1_z_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 585);

    auto ta1_z_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 586);

    auto ta1_z_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 587);

    auto ta1_z_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 588);

    auto ta1_z_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 589);

    auto ta1_z_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 590);

    auto ta1_z_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 591);

    auto ta1_z_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 592);

    auto ta1_z_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 593);

    auto ta1_z_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 594);

    auto ta1_z_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 595);

    auto ta1_z_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 596);

    auto ta1_z_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 597);

    auto ta1_z_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 598);

    auto ta1_z_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 599);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xzzz_xxxx_0, \
                             ta1_z_xzzz_xxxy_0, \
                             ta1_z_xzzz_xxxz_0, \
                             ta1_z_xzzz_xxyy_0, \
                             ta1_z_xzzz_xxyz_0, \
                             ta1_z_xzzz_xxzz_0, \
                             ta1_z_xzzz_xyyy_0, \
                             ta1_z_xzzz_xyyz_0, \
                             ta1_z_xzzz_xyzz_0, \
                             ta1_z_xzzz_xzzz_0, \
                             ta1_z_xzzz_yyyy_0, \
                             ta1_z_xzzz_yyyz_0, \
                             ta1_z_xzzz_yyzz_0, \
                             ta1_z_xzzz_yzzz_0, \
                             ta1_z_xzzz_zzzz_0, \
                             ta1_z_zzz_xxx_0,   \
                             ta1_z_zzz_xxx_1,   \
                             ta1_z_zzz_xxxx_0,  \
                             ta1_z_zzz_xxxx_1,  \
                             ta1_z_zzz_xxxy_0,  \
                             ta1_z_zzz_xxxy_1,  \
                             ta1_z_zzz_xxxz_0,  \
                             ta1_z_zzz_xxxz_1,  \
                             ta1_z_zzz_xxy_0,   \
                             ta1_z_zzz_xxy_1,   \
                             ta1_z_zzz_xxyy_0,  \
                             ta1_z_zzz_xxyy_1,  \
                             ta1_z_zzz_xxyz_0,  \
                             ta1_z_zzz_xxyz_1,  \
                             ta1_z_zzz_xxz_0,   \
                             ta1_z_zzz_xxz_1,   \
                             ta1_z_zzz_xxzz_0,  \
                             ta1_z_zzz_xxzz_1,  \
                             ta1_z_zzz_xyy_0,   \
                             ta1_z_zzz_xyy_1,   \
                             ta1_z_zzz_xyyy_0,  \
                             ta1_z_zzz_xyyy_1,  \
                             ta1_z_zzz_xyyz_0,  \
                             ta1_z_zzz_xyyz_1,  \
                             ta1_z_zzz_xyz_0,   \
                             ta1_z_zzz_xyz_1,   \
                             ta1_z_zzz_xyzz_0,  \
                             ta1_z_zzz_xyzz_1,  \
                             ta1_z_zzz_xzz_0,   \
                             ta1_z_zzz_xzz_1,   \
                             ta1_z_zzz_xzzz_0,  \
                             ta1_z_zzz_xzzz_1,  \
                             ta1_z_zzz_yyy_0,   \
                             ta1_z_zzz_yyy_1,   \
                             ta1_z_zzz_yyyy_0,  \
                             ta1_z_zzz_yyyy_1,  \
                             ta1_z_zzz_yyyz_0,  \
                             ta1_z_zzz_yyyz_1,  \
                             ta1_z_zzz_yyz_0,   \
                             ta1_z_zzz_yyz_1,   \
                             ta1_z_zzz_yyzz_0,  \
                             ta1_z_zzz_yyzz_1,  \
                             ta1_z_zzz_yzz_0,   \
                             ta1_z_zzz_yzz_1,   \
                             ta1_z_zzz_yzzz_0,  \
                             ta1_z_zzz_yzzz_1,  \
                             ta1_z_zzz_zzz_0,   \
                             ta1_z_zzz_zzz_1,   \
                             ta1_z_zzz_zzzz_0,  \
                             ta1_z_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzz_xxxx_0[i] =
            4.0 * ta1_z_zzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxx_1[i] * fe_0 + ta1_z_zzz_xxxx_0[i] * pa_x[i] - ta1_z_zzz_xxxx_1[i] * pc_x[i];

        ta1_z_xzzz_xxxy_0[i] =
            3.0 * ta1_z_zzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxy_1[i] * fe_0 + ta1_z_zzz_xxxy_0[i] * pa_x[i] - ta1_z_zzz_xxxy_1[i] * pc_x[i];

        ta1_z_xzzz_xxxz_0[i] =
            3.0 * ta1_z_zzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxz_1[i] * fe_0 + ta1_z_zzz_xxxz_0[i] * pa_x[i] - ta1_z_zzz_xxxz_1[i] * pc_x[i];

        ta1_z_xzzz_xxyy_0[i] =
            2.0 * ta1_z_zzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyy_1[i] * fe_0 + ta1_z_zzz_xxyy_0[i] * pa_x[i] - ta1_z_zzz_xxyy_1[i] * pc_x[i];

        ta1_z_xzzz_xxyz_0[i] =
            2.0 * ta1_z_zzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyz_1[i] * fe_0 + ta1_z_zzz_xxyz_0[i] * pa_x[i] - ta1_z_zzz_xxyz_1[i] * pc_x[i];

        ta1_z_xzzz_xxzz_0[i] =
            2.0 * ta1_z_zzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xzz_1[i] * fe_0 + ta1_z_zzz_xxzz_0[i] * pa_x[i] - ta1_z_zzz_xxzz_1[i] * pc_x[i];

        ta1_z_xzzz_xyyy_0[i] = ta1_z_zzz_yyy_0[i] * fe_0 - ta1_z_zzz_yyy_1[i] * fe_0 + ta1_z_zzz_xyyy_0[i] * pa_x[i] - ta1_z_zzz_xyyy_1[i] * pc_x[i];

        ta1_z_xzzz_xyyz_0[i] = ta1_z_zzz_yyz_0[i] * fe_0 - ta1_z_zzz_yyz_1[i] * fe_0 + ta1_z_zzz_xyyz_0[i] * pa_x[i] - ta1_z_zzz_xyyz_1[i] * pc_x[i];

        ta1_z_xzzz_xyzz_0[i] = ta1_z_zzz_yzz_0[i] * fe_0 - ta1_z_zzz_yzz_1[i] * fe_0 + ta1_z_zzz_xyzz_0[i] * pa_x[i] - ta1_z_zzz_xyzz_1[i] * pc_x[i];

        ta1_z_xzzz_xzzz_0[i] = ta1_z_zzz_zzz_0[i] * fe_0 - ta1_z_zzz_zzz_1[i] * fe_0 + ta1_z_zzz_xzzz_0[i] * pa_x[i] - ta1_z_zzz_xzzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyyy_0[i] = ta1_z_zzz_yyyy_0[i] * pa_x[i] - ta1_z_zzz_yyyy_1[i] * pc_x[i];

        ta1_z_xzzz_yyyz_0[i] = ta1_z_zzz_yyyz_0[i] * pa_x[i] - ta1_z_zzz_yyyz_1[i] * pc_x[i];

        ta1_z_xzzz_yyzz_0[i] = ta1_z_zzz_yyzz_0[i] * pa_x[i] - ta1_z_zzz_yyzz_1[i] * pc_x[i];

        ta1_z_xzzz_yzzz_0[i] = ta1_z_zzz_yzzz_0[i] * pa_x[i] - ta1_z_zzz_yzzz_1[i] * pc_x[i];

        ta1_z_xzzz_zzzz_0[i] = ta1_z_zzz_zzzz_0[i] * pa_x[i] - ta1_z_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 600-615 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_z_yy_xxxx_0,   \
                             ta1_z_yy_xxxx_1,   \
                             ta1_z_yy_xxxy_0,   \
                             ta1_z_yy_xxxy_1,   \
                             ta1_z_yy_xxxz_0,   \
                             ta1_z_yy_xxxz_1,   \
                             ta1_z_yy_xxyy_0,   \
                             ta1_z_yy_xxyy_1,   \
                             ta1_z_yy_xxyz_0,   \
                             ta1_z_yy_xxyz_1,   \
                             ta1_z_yy_xxzz_0,   \
                             ta1_z_yy_xxzz_1,   \
                             ta1_z_yy_xyyy_0,   \
                             ta1_z_yy_xyyy_1,   \
                             ta1_z_yy_xyyz_0,   \
                             ta1_z_yy_xyyz_1,   \
                             ta1_z_yy_xyzz_0,   \
                             ta1_z_yy_xyzz_1,   \
                             ta1_z_yy_xzzz_0,   \
                             ta1_z_yy_xzzz_1,   \
                             ta1_z_yy_yyyy_0,   \
                             ta1_z_yy_yyyy_1,   \
                             ta1_z_yy_yyyz_0,   \
                             ta1_z_yy_yyyz_1,   \
                             ta1_z_yy_yyzz_0,   \
                             ta1_z_yy_yyzz_1,   \
                             ta1_z_yy_yzzz_0,   \
                             ta1_z_yy_yzzz_1,   \
                             ta1_z_yy_zzzz_0,   \
                             ta1_z_yy_zzzz_1,   \
                             ta1_z_yyy_xxx_0,   \
                             ta1_z_yyy_xxx_1,   \
                             ta1_z_yyy_xxxx_0,  \
                             ta1_z_yyy_xxxx_1,  \
                             ta1_z_yyy_xxxy_0,  \
                             ta1_z_yyy_xxxy_1,  \
                             ta1_z_yyy_xxxz_0,  \
                             ta1_z_yyy_xxxz_1,  \
                             ta1_z_yyy_xxy_0,   \
                             ta1_z_yyy_xxy_1,   \
                             ta1_z_yyy_xxyy_0,  \
                             ta1_z_yyy_xxyy_1,  \
                             ta1_z_yyy_xxyz_0,  \
                             ta1_z_yyy_xxyz_1,  \
                             ta1_z_yyy_xxz_0,   \
                             ta1_z_yyy_xxz_1,   \
                             ta1_z_yyy_xxzz_0,  \
                             ta1_z_yyy_xxzz_1,  \
                             ta1_z_yyy_xyy_0,   \
                             ta1_z_yyy_xyy_1,   \
                             ta1_z_yyy_xyyy_0,  \
                             ta1_z_yyy_xyyy_1,  \
                             ta1_z_yyy_xyyz_0,  \
                             ta1_z_yyy_xyyz_1,  \
                             ta1_z_yyy_xyz_0,   \
                             ta1_z_yyy_xyz_1,   \
                             ta1_z_yyy_xyzz_0,  \
                             ta1_z_yyy_xyzz_1,  \
                             ta1_z_yyy_xzz_0,   \
                             ta1_z_yyy_xzz_1,   \
                             ta1_z_yyy_xzzz_0,  \
                             ta1_z_yyy_xzzz_1,  \
                             ta1_z_yyy_yyy_0,   \
                             ta1_z_yyy_yyy_1,   \
                             ta1_z_yyy_yyyy_0,  \
                             ta1_z_yyy_yyyy_1,  \
                             ta1_z_yyy_yyyz_0,  \
                             ta1_z_yyy_yyyz_1,  \
                             ta1_z_yyy_yyz_0,   \
                             ta1_z_yyy_yyz_1,   \
                             ta1_z_yyy_yyzz_0,  \
                             ta1_z_yyy_yyzz_1,  \
                             ta1_z_yyy_yzz_0,   \
                             ta1_z_yyy_yzz_1,   \
                             ta1_z_yyy_yzzz_0,  \
                             ta1_z_yyy_yzzz_1,  \
                             ta1_z_yyy_zzz_0,   \
                             ta1_z_yyy_zzz_1,   \
                             ta1_z_yyy_zzzz_0,  \
                             ta1_z_yyy_zzzz_1,  \
                             ta1_z_yyyy_xxxx_0, \
                             ta1_z_yyyy_xxxy_0, \
                             ta1_z_yyyy_xxxz_0, \
                             ta1_z_yyyy_xxyy_0, \
                             ta1_z_yyyy_xxyz_0, \
                             ta1_z_yyyy_xxzz_0, \
                             ta1_z_yyyy_xyyy_0, \
                             ta1_z_yyyy_xyyz_0, \
                             ta1_z_yyyy_xyzz_0, \
                             ta1_z_yyyy_xzzz_0, \
                             ta1_z_yyyy_yyyy_0, \
                             ta1_z_yyyy_yyyz_0, \
                             ta1_z_yyyy_yyzz_0, \
                             ta1_z_yyyy_yzzz_0, \
                             ta1_z_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyy_xxxx_0[i] =
            3.0 * ta1_z_yy_xxxx_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxx_1[i] * fe_0 + ta1_z_yyy_xxxx_0[i] * pa_y[i] - ta1_z_yyy_xxxx_1[i] * pc_y[i];

        ta1_z_yyyy_xxxy_0[i] = 3.0 * ta1_z_yy_xxxy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxy_1[i] * fe_0 + ta1_z_yyy_xxx_0[i] * fe_0 -
                               ta1_z_yyy_xxx_1[i] * fe_0 + ta1_z_yyy_xxxy_0[i] * pa_y[i] - ta1_z_yyy_xxxy_1[i] * pc_y[i];

        ta1_z_yyyy_xxxz_0[i] =
            3.0 * ta1_z_yy_xxxz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxz_1[i] * fe_0 + ta1_z_yyy_xxxz_0[i] * pa_y[i] - ta1_z_yyy_xxxz_1[i] * pc_y[i];

        ta1_z_yyyy_xxyy_0[i] = 3.0 * ta1_z_yy_xxyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxy_0[i] * fe_0 -
                               2.0 * ta1_z_yyy_xxy_1[i] * fe_0 + ta1_z_yyy_xxyy_0[i] * pa_y[i] - ta1_z_yyy_xxyy_1[i] * pc_y[i];

        ta1_z_yyyy_xxyz_0[i] = 3.0 * ta1_z_yy_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyz_1[i] * fe_0 + ta1_z_yyy_xxz_0[i] * fe_0 -
                               ta1_z_yyy_xxz_1[i] * fe_0 + ta1_z_yyy_xxyz_0[i] * pa_y[i] - ta1_z_yyy_xxyz_1[i] * pc_y[i];

        ta1_z_yyyy_xxzz_0[i] =
            3.0 * ta1_z_yy_xxzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxzz_1[i] * fe_0 + ta1_z_yyy_xxzz_0[i] * pa_y[i] - ta1_z_yyy_xxzz_1[i] * pc_y[i];

        ta1_z_yyyy_xyyy_0[i] = 3.0 * ta1_z_yy_xyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyy_1[i] * fe_0 + 3.0 * ta1_z_yyy_xyy_0[i] * fe_0 -
                               3.0 * ta1_z_yyy_xyy_1[i] * fe_0 + ta1_z_yyy_xyyy_0[i] * pa_y[i] - ta1_z_yyy_xyyy_1[i] * pc_y[i];

        ta1_z_yyyy_xyyz_0[i] = 3.0 * ta1_z_yy_xyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xyz_0[i] * fe_0 -
                               2.0 * ta1_z_yyy_xyz_1[i] * fe_0 + ta1_z_yyy_xyyz_0[i] * pa_y[i] - ta1_z_yyy_xyyz_1[i] * pc_y[i];

        ta1_z_yyyy_xyzz_0[i] = 3.0 * ta1_z_yy_xyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyzz_1[i] * fe_0 + ta1_z_yyy_xzz_0[i] * fe_0 -
                               ta1_z_yyy_xzz_1[i] * fe_0 + ta1_z_yyy_xyzz_0[i] * pa_y[i] - ta1_z_yyy_xyzz_1[i] * pc_y[i];

        ta1_z_yyyy_xzzz_0[i] =
            3.0 * ta1_z_yy_xzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xzzz_1[i] * fe_0 + ta1_z_yyy_xzzz_0[i] * pa_y[i] - ta1_z_yyy_xzzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyyy_0[i] = 3.0 * ta1_z_yy_yyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyy_1[i] * fe_0 + 4.0 * ta1_z_yyy_yyy_0[i] * fe_0 -
                               4.0 * ta1_z_yyy_yyy_1[i] * fe_0 + ta1_z_yyy_yyyy_0[i] * pa_y[i] - ta1_z_yyy_yyyy_1[i] * pc_y[i];

        ta1_z_yyyy_yyyz_0[i] = 3.0 * ta1_z_yy_yyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yyy_yyz_0[i] * fe_0 -
                               3.0 * ta1_z_yyy_yyz_1[i] * fe_0 + ta1_z_yyy_yyyz_0[i] * pa_y[i] - ta1_z_yyy_yyyz_1[i] * pc_y[i];

        ta1_z_yyyy_yyzz_0[i] = 3.0 * ta1_z_yy_yyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yzz_0[i] * fe_0 -
                               2.0 * ta1_z_yyy_yzz_1[i] * fe_0 + ta1_z_yyy_yyzz_0[i] * pa_y[i] - ta1_z_yyy_yyzz_1[i] * pc_y[i];

        ta1_z_yyyy_yzzz_0[i] = 3.0 * ta1_z_yy_yzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yzzz_1[i] * fe_0 + ta1_z_yyy_zzz_0[i] * fe_0 -
                               ta1_z_yyy_zzz_1[i] * fe_0 + ta1_z_yyy_yzzz_0[i] * pa_y[i] - ta1_z_yyy_yzzz_1[i] * pc_y[i];

        ta1_z_yyyy_zzzz_0[i] =
            3.0 * ta1_z_yy_zzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_zzzz_1[i] * fe_0 + ta1_z_yyy_zzzz_0[i] * pa_y[i] - ta1_z_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 615-630 components of targeted buffer : GG

    auto ta1_z_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 615);

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yyy_xxxx_0,  \
                             ta1_z_yyy_xxxx_1,  \
                             ta1_z_yyy_xxxy_0,  \
                             ta1_z_yyy_xxxy_1,  \
                             ta1_z_yyy_xxy_0,   \
                             ta1_z_yyy_xxy_1,   \
                             ta1_z_yyy_xxyy_0,  \
                             ta1_z_yyy_xxyy_1,  \
                             ta1_z_yyy_xxyz_0,  \
                             ta1_z_yyy_xxyz_1,  \
                             ta1_z_yyy_xyy_0,   \
                             ta1_z_yyy_xyy_1,   \
                             ta1_z_yyy_xyyy_0,  \
                             ta1_z_yyy_xyyy_1,  \
                             ta1_z_yyy_xyyz_0,  \
                             ta1_z_yyy_xyyz_1,  \
                             ta1_z_yyy_xyz_0,   \
                             ta1_z_yyy_xyz_1,   \
                             ta1_z_yyy_xyzz_0,  \
                             ta1_z_yyy_xyzz_1,  \
                             ta1_z_yyy_yyy_0,   \
                             ta1_z_yyy_yyy_1,   \
                             ta1_z_yyy_yyyy_0,  \
                             ta1_z_yyy_yyyy_1,  \
                             ta1_z_yyy_yyyz_0,  \
                             ta1_z_yyy_yyyz_1,  \
                             ta1_z_yyy_yyz_0,   \
                             ta1_z_yyy_yyz_1,   \
                             ta1_z_yyy_yyzz_0,  \
                             ta1_z_yyy_yyzz_1,  \
                             ta1_z_yyy_yzz_0,   \
                             ta1_z_yyy_yzz_1,   \
                             ta1_z_yyy_yzzz_0,  \
                             ta1_z_yyy_yzzz_1,  \
                             ta1_z_yyyz_xxxx_0, \
                             ta1_z_yyyz_xxxy_0, \
                             ta1_z_yyyz_xxxz_0, \
                             ta1_z_yyyz_xxyy_0, \
                             ta1_z_yyyz_xxyz_0, \
                             ta1_z_yyyz_xxzz_0, \
                             ta1_z_yyyz_xyyy_0, \
                             ta1_z_yyyz_xyyz_0, \
                             ta1_z_yyyz_xyzz_0, \
                             ta1_z_yyyz_xzzz_0, \
                             ta1_z_yyyz_yyyy_0, \
                             ta1_z_yyyz_yyyz_0, \
                             ta1_z_yyyz_yyzz_0, \
                             ta1_z_yyyz_yzzz_0, \
                             ta1_z_yyyz_zzzz_0, \
                             ta1_z_yyz_xxxz_0,  \
                             ta1_z_yyz_xxxz_1,  \
                             ta1_z_yyz_xxzz_0,  \
                             ta1_z_yyz_xxzz_1,  \
                             ta1_z_yyz_xzzz_0,  \
                             ta1_z_yyz_xzzz_1,  \
                             ta1_z_yyz_zzzz_0,  \
                             ta1_z_yyz_zzzz_1,  \
                             ta1_z_yz_xxxz_0,   \
                             ta1_z_yz_xxxz_1,   \
                             ta1_z_yz_xxzz_0,   \
                             ta1_z_yz_xxzz_1,   \
                             ta1_z_yz_xzzz_0,   \
                             ta1_z_yz_xzzz_1,   \
                             ta1_z_yz_zzzz_0,   \
                             ta1_z_yz_zzzz_1,   \
                             ta_yyy_xxxx_1,     \
                             ta_yyy_xxxy_1,     \
                             ta_yyy_xxyy_1,     \
                             ta_yyy_xxyz_1,     \
                             ta_yyy_xyyy_1,     \
                             ta_yyy_xyyz_1,     \
                             ta_yyy_xyzz_1,     \
                             ta_yyy_yyyy_1,     \
                             ta_yyy_yyyz_1,     \
                             ta_yyy_yyzz_1,     \
                             ta_yyy_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyz_xxxx_0[i] = ta_yyy_xxxx_1[i] + ta1_z_yyy_xxxx_0[i] * pa_z[i] - ta1_z_yyy_xxxx_1[i] * pc_z[i];

        ta1_z_yyyz_xxxy_0[i] = ta_yyy_xxxy_1[i] + ta1_z_yyy_xxxy_0[i] * pa_z[i] - ta1_z_yyy_xxxy_1[i] * pc_z[i];

        ta1_z_yyyz_xxxz_0[i] =
            2.0 * ta1_z_yz_xxxz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxxz_1[i] * fe_0 + ta1_z_yyz_xxxz_0[i] * pa_y[i] - ta1_z_yyz_xxxz_1[i] * pc_y[i];

        ta1_z_yyyz_xxyy_0[i] = ta_yyy_xxyy_1[i] + ta1_z_yyy_xxyy_0[i] * pa_z[i] - ta1_z_yyy_xxyy_1[i] * pc_z[i];

        ta1_z_yyyz_xxyz_0[i] =
            ta1_z_yyy_xxy_0[i] * fe_0 - ta1_z_yyy_xxy_1[i] * fe_0 + ta_yyy_xxyz_1[i] + ta1_z_yyy_xxyz_0[i] * pa_z[i] - ta1_z_yyy_xxyz_1[i] * pc_z[i];

        ta1_z_yyyz_xxzz_0[i] =
            2.0 * ta1_z_yz_xxzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxzz_1[i] * fe_0 + ta1_z_yyz_xxzz_0[i] * pa_y[i] - ta1_z_yyz_xxzz_1[i] * pc_y[i];

        ta1_z_yyyz_xyyy_0[i] = ta_yyy_xyyy_1[i] + ta1_z_yyy_xyyy_0[i] * pa_z[i] - ta1_z_yyy_xyyy_1[i] * pc_z[i];

        ta1_z_yyyz_xyyz_0[i] =
            ta1_z_yyy_xyy_0[i] * fe_0 - ta1_z_yyy_xyy_1[i] * fe_0 + ta_yyy_xyyz_1[i] + ta1_z_yyy_xyyz_0[i] * pa_z[i] - ta1_z_yyy_xyyz_1[i] * pc_z[i];

        ta1_z_yyyz_xyzz_0[i] = 2.0 * ta1_z_yyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyz_1[i] * fe_0 + ta_yyy_xyzz_1[i] + ta1_z_yyy_xyzz_0[i] * pa_z[i] -
                               ta1_z_yyy_xyzz_1[i] * pc_z[i];

        ta1_z_yyyz_xzzz_0[i] =
            2.0 * ta1_z_yz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xzzz_1[i] * fe_0 + ta1_z_yyz_xzzz_0[i] * pa_y[i] - ta1_z_yyz_xzzz_1[i] * pc_y[i];

        ta1_z_yyyz_yyyy_0[i] = ta_yyy_yyyy_1[i] + ta1_z_yyy_yyyy_0[i] * pa_z[i] - ta1_z_yyy_yyyy_1[i] * pc_z[i];

        ta1_z_yyyz_yyyz_0[i] =
            ta1_z_yyy_yyy_0[i] * fe_0 - ta1_z_yyy_yyy_1[i] * fe_0 + ta_yyy_yyyz_1[i] + ta1_z_yyy_yyyz_0[i] * pa_z[i] - ta1_z_yyy_yyyz_1[i] * pc_z[i];

        ta1_z_yyyz_yyzz_0[i] = 2.0 * ta1_z_yyy_yyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_yyz_1[i] * fe_0 + ta_yyy_yyzz_1[i] + ta1_z_yyy_yyzz_0[i] * pa_z[i] -
                               ta1_z_yyy_yyzz_1[i] * pc_z[i];

        ta1_z_yyyz_yzzz_0[i] = 3.0 * ta1_z_yyy_yzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_yzz_1[i] * fe_0 + ta_yyy_yzzz_1[i] + ta1_z_yyy_yzzz_0[i] * pa_z[i] -
                               ta1_z_yyy_yzzz_1[i] * pc_z[i];

        ta1_z_yyyz_zzzz_0[i] =
            2.0 * ta1_z_yz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_zzzz_1[i] * fe_0 + ta1_z_yyz_zzzz_0[i] * pa_y[i] - ta1_z_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 630-645 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yy_xxxy_0,   \
                             ta1_z_yy_xxxy_1,   \
                             ta1_z_yy_xxyy_0,   \
                             ta1_z_yy_xxyy_1,   \
                             ta1_z_yy_xyyy_0,   \
                             ta1_z_yy_xyyy_1,   \
                             ta1_z_yy_yyyy_0,   \
                             ta1_z_yy_yyyy_1,   \
                             ta1_z_yyz_xxxy_0,  \
                             ta1_z_yyz_xxxy_1,  \
                             ta1_z_yyz_xxyy_0,  \
                             ta1_z_yyz_xxyy_1,  \
                             ta1_z_yyz_xyyy_0,  \
                             ta1_z_yyz_xyyy_1,  \
                             ta1_z_yyz_yyyy_0,  \
                             ta1_z_yyz_yyyy_1,  \
                             ta1_z_yyzz_xxxx_0, \
                             ta1_z_yyzz_xxxy_0, \
                             ta1_z_yyzz_xxxz_0, \
                             ta1_z_yyzz_xxyy_0, \
                             ta1_z_yyzz_xxyz_0, \
                             ta1_z_yyzz_xxzz_0, \
                             ta1_z_yyzz_xyyy_0, \
                             ta1_z_yyzz_xyyz_0, \
                             ta1_z_yyzz_xyzz_0, \
                             ta1_z_yyzz_xzzz_0, \
                             ta1_z_yyzz_yyyy_0, \
                             ta1_z_yyzz_yyyz_0, \
                             ta1_z_yyzz_yyzz_0, \
                             ta1_z_yyzz_yzzz_0, \
                             ta1_z_yyzz_zzzz_0, \
                             ta1_z_yzz_xxxx_0,  \
                             ta1_z_yzz_xxxx_1,  \
                             ta1_z_yzz_xxxz_0,  \
                             ta1_z_yzz_xxxz_1,  \
                             ta1_z_yzz_xxyz_0,  \
                             ta1_z_yzz_xxyz_1,  \
                             ta1_z_yzz_xxz_0,   \
                             ta1_z_yzz_xxz_1,   \
                             ta1_z_yzz_xxzz_0,  \
                             ta1_z_yzz_xxzz_1,  \
                             ta1_z_yzz_xyyz_0,  \
                             ta1_z_yzz_xyyz_1,  \
                             ta1_z_yzz_xyz_0,   \
                             ta1_z_yzz_xyz_1,   \
                             ta1_z_yzz_xyzz_0,  \
                             ta1_z_yzz_xyzz_1,  \
                             ta1_z_yzz_xzz_0,   \
                             ta1_z_yzz_xzz_1,   \
                             ta1_z_yzz_xzzz_0,  \
                             ta1_z_yzz_xzzz_1,  \
                             ta1_z_yzz_yyyz_0,  \
                             ta1_z_yzz_yyyz_1,  \
                             ta1_z_yzz_yyz_0,   \
                             ta1_z_yzz_yyz_1,   \
                             ta1_z_yzz_yyzz_0,  \
                             ta1_z_yzz_yyzz_1,  \
                             ta1_z_yzz_yzz_0,   \
                             ta1_z_yzz_yzz_1,   \
                             ta1_z_yzz_yzzz_0,  \
                             ta1_z_yzz_yzzz_1,  \
                             ta1_z_yzz_zzz_0,   \
                             ta1_z_yzz_zzz_1,   \
                             ta1_z_yzz_zzzz_0,  \
                             ta1_z_yzz_zzzz_1,  \
                             ta1_z_zz_xxxx_0,   \
                             ta1_z_zz_xxxx_1,   \
                             ta1_z_zz_xxxz_0,   \
                             ta1_z_zz_xxxz_1,   \
                             ta1_z_zz_xxyz_0,   \
                             ta1_z_zz_xxyz_1,   \
                             ta1_z_zz_xxzz_0,   \
                             ta1_z_zz_xxzz_1,   \
                             ta1_z_zz_xyyz_0,   \
                             ta1_z_zz_xyyz_1,   \
                             ta1_z_zz_xyzz_0,   \
                             ta1_z_zz_xyzz_1,   \
                             ta1_z_zz_xzzz_0,   \
                             ta1_z_zz_xzzz_1,   \
                             ta1_z_zz_yyyz_0,   \
                             ta1_z_zz_yyyz_1,   \
                             ta1_z_zz_yyzz_0,   \
                             ta1_z_zz_yyzz_1,   \
                             ta1_z_zz_yzzz_0,   \
                             ta1_z_zz_yzzz_1,   \
                             ta1_z_zz_zzzz_0,   \
                             ta1_z_zz_zzzz_1,   \
                             ta_yyz_xxxy_1,     \
                             ta_yyz_xxyy_1,     \
                             ta_yyz_xyyy_1,     \
                             ta_yyz_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzz_xxxx_0[i] = ta1_z_zz_xxxx_0[i] * fe_0 - ta1_z_zz_xxxx_1[i] * fe_0 + ta1_z_yzz_xxxx_0[i] * pa_y[i] - ta1_z_yzz_xxxx_1[i] * pc_y[i];

        ta1_z_yyzz_xxxy_0[i] =
            ta1_z_yy_xxxy_0[i] * fe_0 - ta1_z_yy_xxxy_1[i] * fe_0 + ta_yyz_xxxy_1[i] + ta1_z_yyz_xxxy_0[i] * pa_z[i] - ta1_z_yyz_xxxy_1[i] * pc_z[i];

        ta1_z_yyzz_xxxz_0[i] = ta1_z_zz_xxxz_0[i] * fe_0 - ta1_z_zz_xxxz_1[i] * fe_0 + ta1_z_yzz_xxxz_0[i] * pa_y[i] - ta1_z_yzz_xxxz_1[i] * pc_y[i];

        ta1_z_yyzz_xxyy_0[i] =
            ta1_z_yy_xxyy_0[i] * fe_0 - ta1_z_yy_xxyy_1[i] * fe_0 + ta_yyz_xxyy_1[i] + ta1_z_yyz_xxyy_0[i] * pa_z[i] - ta1_z_yyz_xxyy_1[i] * pc_z[i];

        ta1_z_yyzz_xxyz_0[i] = ta1_z_zz_xxyz_0[i] * fe_0 - ta1_z_zz_xxyz_1[i] * fe_0 + ta1_z_yzz_xxz_0[i] * fe_0 - ta1_z_yzz_xxz_1[i] * fe_0 +
                               ta1_z_yzz_xxyz_0[i] * pa_y[i] - ta1_z_yzz_xxyz_1[i] * pc_y[i];

        ta1_z_yyzz_xxzz_0[i] = ta1_z_zz_xxzz_0[i] * fe_0 - ta1_z_zz_xxzz_1[i] * fe_0 + ta1_z_yzz_xxzz_0[i] * pa_y[i] - ta1_z_yzz_xxzz_1[i] * pc_y[i];

        ta1_z_yyzz_xyyy_0[i] =
            ta1_z_yy_xyyy_0[i] * fe_0 - ta1_z_yy_xyyy_1[i] * fe_0 + ta_yyz_xyyy_1[i] + ta1_z_yyz_xyyy_0[i] * pa_z[i] - ta1_z_yyz_xyyy_1[i] * pc_z[i];

        ta1_z_yyzz_xyyz_0[i] = ta1_z_zz_xyyz_0[i] * fe_0 - ta1_z_zz_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_z_yzz_xyz_1[i] * fe_0 + ta1_z_yzz_xyyz_0[i] * pa_y[i] - ta1_z_yzz_xyyz_1[i] * pc_y[i];

        ta1_z_yyzz_xyzz_0[i] = ta1_z_zz_xyzz_0[i] * fe_0 - ta1_z_zz_xyzz_1[i] * fe_0 + ta1_z_yzz_xzz_0[i] * fe_0 - ta1_z_yzz_xzz_1[i] * fe_0 +
                               ta1_z_yzz_xyzz_0[i] * pa_y[i] - ta1_z_yzz_xyzz_1[i] * pc_y[i];

        ta1_z_yyzz_xzzz_0[i] = ta1_z_zz_xzzz_0[i] * fe_0 - ta1_z_zz_xzzz_1[i] * fe_0 + ta1_z_yzz_xzzz_0[i] * pa_y[i] - ta1_z_yzz_xzzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyyy_0[i] =
            ta1_z_yy_yyyy_0[i] * fe_0 - ta1_z_yy_yyyy_1[i] * fe_0 + ta_yyz_yyyy_1[i] + ta1_z_yyz_yyyy_0[i] * pa_z[i] - ta1_z_yyz_yyyy_1[i] * pc_z[i];

        ta1_z_yyzz_yyyz_0[i] = ta1_z_zz_yyyz_0[i] * fe_0 - ta1_z_zz_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yzz_yyz_0[i] * fe_0 -
                               3.0 * ta1_z_yzz_yyz_1[i] * fe_0 + ta1_z_yzz_yyyz_0[i] * pa_y[i] - ta1_z_yzz_yyyz_1[i] * pc_y[i];

        ta1_z_yyzz_yyzz_0[i] = ta1_z_zz_yyzz_0[i] * fe_0 - ta1_z_zz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yzz_yzz_0[i] * fe_0 -
                               2.0 * ta1_z_yzz_yzz_1[i] * fe_0 + ta1_z_yzz_yyzz_0[i] * pa_y[i] - ta1_z_yzz_yyzz_1[i] * pc_y[i];

        ta1_z_yyzz_yzzz_0[i] = ta1_z_zz_yzzz_0[i] * fe_0 - ta1_z_zz_yzzz_1[i] * fe_0 + ta1_z_yzz_zzz_0[i] * fe_0 - ta1_z_yzz_zzz_1[i] * fe_0 +
                               ta1_z_yzz_yzzz_0[i] * pa_y[i] - ta1_z_yzz_yzzz_1[i] * pc_y[i];

        ta1_z_yyzz_zzzz_0[i] = ta1_z_zz_zzzz_0[i] * fe_0 - ta1_z_zz_zzzz_1[i] * fe_0 + ta1_z_yzz_zzzz_0[i] * pa_y[i] - ta1_z_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 645-660 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_z_yzzz_xxxx_0, \
                             ta1_z_yzzz_xxxy_0, \
                             ta1_z_yzzz_xxxz_0, \
                             ta1_z_yzzz_xxyy_0, \
                             ta1_z_yzzz_xxyz_0, \
                             ta1_z_yzzz_xxzz_0, \
                             ta1_z_yzzz_xyyy_0, \
                             ta1_z_yzzz_xyyz_0, \
                             ta1_z_yzzz_xyzz_0, \
                             ta1_z_yzzz_xzzz_0, \
                             ta1_z_yzzz_yyyy_0, \
                             ta1_z_yzzz_yyyz_0, \
                             ta1_z_yzzz_yyzz_0, \
                             ta1_z_yzzz_yzzz_0, \
                             ta1_z_yzzz_zzzz_0, \
                             ta1_z_zzz_xxx_0,   \
                             ta1_z_zzz_xxx_1,   \
                             ta1_z_zzz_xxxx_0,  \
                             ta1_z_zzz_xxxx_1,  \
                             ta1_z_zzz_xxxy_0,  \
                             ta1_z_zzz_xxxy_1,  \
                             ta1_z_zzz_xxxz_0,  \
                             ta1_z_zzz_xxxz_1,  \
                             ta1_z_zzz_xxy_0,   \
                             ta1_z_zzz_xxy_1,   \
                             ta1_z_zzz_xxyy_0,  \
                             ta1_z_zzz_xxyy_1,  \
                             ta1_z_zzz_xxyz_0,  \
                             ta1_z_zzz_xxyz_1,  \
                             ta1_z_zzz_xxz_0,   \
                             ta1_z_zzz_xxz_1,   \
                             ta1_z_zzz_xxzz_0,  \
                             ta1_z_zzz_xxzz_1,  \
                             ta1_z_zzz_xyy_0,   \
                             ta1_z_zzz_xyy_1,   \
                             ta1_z_zzz_xyyy_0,  \
                             ta1_z_zzz_xyyy_1,  \
                             ta1_z_zzz_xyyz_0,  \
                             ta1_z_zzz_xyyz_1,  \
                             ta1_z_zzz_xyz_0,   \
                             ta1_z_zzz_xyz_1,   \
                             ta1_z_zzz_xyzz_0,  \
                             ta1_z_zzz_xyzz_1,  \
                             ta1_z_zzz_xzz_0,   \
                             ta1_z_zzz_xzz_1,   \
                             ta1_z_zzz_xzzz_0,  \
                             ta1_z_zzz_xzzz_1,  \
                             ta1_z_zzz_yyy_0,   \
                             ta1_z_zzz_yyy_1,   \
                             ta1_z_zzz_yyyy_0,  \
                             ta1_z_zzz_yyyy_1,  \
                             ta1_z_zzz_yyyz_0,  \
                             ta1_z_zzz_yyyz_1,  \
                             ta1_z_zzz_yyz_0,   \
                             ta1_z_zzz_yyz_1,   \
                             ta1_z_zzz_yyzz_0,  \
                             ta1_z_zzz_yyzz_1,  \
                             ta1_z_zzz_yzz_0,   \
                             ta1_z_zzz_yzz_1,   \
                             ta1_z_zzz_yzzz_0,  \
                             ta1_z_zzz_yzzz_1,  \
                             ta1_z_zzz_zzz_0,   \
                             ta1_z_zzz_zzz_1,   \
                             ta1_z_zzz_zzzz_0,  \
                             ta1_z_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzz_xxxx_0[i] = ta1_z_zzz_xxxx_0[i] * pa_y[i] - ta1_z_zzz_xxxx_1[i] * pc_y[i];

        ta1_z_yzzz_xxxy_0[i] = ta1_z_zzz_xxx_0[i] * fe_0 - ta1_z_zzz_xxx_1[i] * fe_0 + ta1_z_zzz_xxxy_0[i] * pa_y[i] - ta1_z_zzz_xxxy_1[i] * pc_y[i];

        ta1_z_yzzz_xxxz_0[i] = ta1_z_zzz_xxxz_0[i] * pa_y[i] - ta1_z_zzz_xxxz_1[i] * pc_y[i];

        ta1_z_yzzz_xxyy_0[i] =
            2.0 * ta1_z_zzz_xxy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xxy_1[i] * fe_0 + ta1_z_zzz_xxyy_0[i] * pa_y[i] - ta1_z_zzz_xxyy_1[i] * pc_y[i];

        ta1_z_yzzz_xxyz_0[i] = ta1_z_zzz_xxz_0[i] * fe_0 - ta1_z_zzz_xxz_1[i] * fe_0 + ta1_z_zzz_xxyz_0[i] * pa_y[i] - ta1_z_zzz_xxyz_1[i] * pc_y[i];

        ta1_z_yzzz_xxzz_0[i] = ta1_z_zzz_xxzz_0[i] * pa_y[i] - ta1_z_zzz_xxzz_1[i] * pc_y[i];

        ta1_z_yzzz_xyyy_0[i] =
            3.0 * ta1_z_zzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_zzz_xyy_1[i] * fe_0 + ta1_z_zzz_xyyy_0[i] * pa_y[i] - ta1_z_zzz_xyyy_1[i] * pc_y[i];

        ta1_z_yzzz_xyyz_0[i] =
            2.0 * ta1_z_zzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyz_1[i] * fe_0 + ta1_z_zzz_xyyz_0[i] * pa_y[i] - ta1_z_zzz_xyyz_1[i] * pc_y[i];

        ta1_z_yzzz_xyzz_0[i] = ta1_z_zzz_xzz_0[i] * fe_0 - ta1_z_zzz_xzz_1[i] * fe_0 + ta1_z_zzz_xyzz_0[i] * pa_y[i] - ta1_z_zzz_xyzz_1[i] * pc_y[i];

        ta1_z_yzzz_xzzz_0[i] = ta1_z_zzz_xzzz_0[i] * pa_y[i] - ta1_z_zzz_xzzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyyy_0[i] =
            4.0 * ta1_z_zzz_yyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyy_1[i] * fe_0 + ta1_z_zzz_yyyy_0[i] * pa_y[i] - ta1_z_zzz_yyyy_1[i] * pc_y[i];

        ta1_z_yzzz_yyyz_0[i] =
            3.0 * ta1_z_zzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_zzz_yyz_1[i] * fe_0 + ta1_z_zzz_yyyz_0[i] * pa_y[i] - ta1_z_zzz_yyyz_1[i] * pc_y[i];

        ta1_z_yzzz_yyzz_0[i] =
            2.0 * ta1_z_zzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_yzz_1[i] * fe_0 + ta1_z_zzz_yyzz_0[i] * pa_y[i] - ta1_z_zzz_yyzz_1[i] * pc_y[i];

        ta1_z_yzzz_yzzz_0[i] = ta1_z_zzz_zzz_0[i] * fe_0 - ta1_z_zzz_zzz_1[i] * fe_0 + ta1_z_zzz_yzzz_0[i] * pa_y[i] - ta1_z_zzz_yzzz_1[i] * pc_y[i];

        ta1_z_yzzz_zzzz_0[i] = ta1_z_zzz_zzzz_0[i] * pa_y[i] - ta1_z_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 660-675 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_z_zz_xxxx_0,   \
                             ta1_z_zz_xxxx_1,   \
                             ta1_z_zz_xxxy_0,   \
                             ta1_z_zz_xxxy_1,   \
                             ta1_z_zz_xxxz_0,   \
                             ta1_z_zz_xxxz_1,   \
                             ta1_z_zz_xxyy_0,   \
                             ta1_z_zz_xxyy_1,   \
                             ta1_z_zz_xxyz_0,   \
                             ta1_z_zz_xxyz_1,   \
                             ta1_z_zz_xxzz_0,   \
                             ta1_z_zz_xxzz_1,   \
                             ta1_z_zz_xyyy_0,   \
                             ta1_z_zz_xyyy_1,   \
                             ta1_z_zz_xyyz_0,   \
                             ta1_z_zz_xyyz_1,   \
                             ta1_z_zz_xyzz_0,   \
                             ta1_z_zz_xyzz_1,   \
                             ta1_z_zz_xzzz_0,   \
                             ta1_z_zz_xzzz_1,   \
                             ta1_z_zz_yyyy_0,   \
                             ta1_z_zz_yyyy_1,   \
                             ta1_z_zz_yyyz_0,   \
                             ta1_z_zz_yyyz_1,   \
                             ta1_z_zz_yyzz_0,   \
                             ta1_z_zz_yyzz_1,   \
                             ta1_z_zz_yzzz_0,   \
                             ta1_z_zz_yzzz_1,   \
                             ta1_z_zz_zzzz_0,   \
                             ta1_z_zz_zzzz_1,   \
                             ta1_z_zzz_xxx_0,   \
                             ta1_z_zzz_xxx_1,   \
                             ta1_z_zzz_xxxx_0,  \
                             ta1_z_zzz_xxxx_1,  \
                             ta1_z_zzz_xxxy_0,  \
                             ta1_z_zzz_xxxy_1,  \
                             ta1_z_zzz_xxxz_0,  \
                             ta1_z_zzz_xxxz_1,  \
                             ta1_z_zzz_xxy_0,   \
                             ta1_z_zzz_xxy_1,   \
                             ta1_z_zzz_xxyy_0,  \
                             ta1_z_zzz_xxyy_1,  \
                             ta1_z_zzz_xxyz_0,  \
                             ta1_z_zzz_xxyz_1,  \
                             ta1_z_zzz_xxz_0,   \
                             ta1_z_zzz_xxz_1,   \
                             ta1_z_zzz_xxzz_0,  \
                             ta1_z_zzz_xxzz_1,  \
                             ta1_z_zzz_xyy_0,   \
                             ta1_z_zzz_xyy_1,   \
                             ta1_z_zzz_xyyy_0,  \
                             ta1_z_zzz_xyyy_1,  \
                             ta1_z_zzz_xyyz_0,  \
                             ta1_z_zzz_xyyz_1,  \
                             ta1_z_zzz_xyz_0,   \
                             ta1_z_zzz_xyz_1,   \
                             ta1_z_zzz_xyzz_0,  \
                             ta1_z_zzz_xyzz_1,  \
                             ta1_z_zzz_xzz_0,   \
                             ta1_z_zzz_xzz_1,   \
                             ta1_z_zzz_xzzz_0,  \
                             ta1_z_zzz_xzzz_1,  \
                             ta1_z_zzz_yyy_0,   \
                             ta1_z_zzz_yyy_1,   \
                             ta1_z_zzz_yyyy_0,  \
                             ta1_z_zzz_yyyy_1,  \
                             ta1_z_zzz_yyyz_0,  \
                             ta1_z_zzz_yyyz_1,  \
                             ta1_z_zzz_yyz_0,   \
                             ta1_z_zzz_yyz_1,   \
                             ta1_z_zzz_yyzz_0,  \
                             ta1_z_zzz_yyzz_1,  \
                             ta1_z_zzz_yzz_0,   \
                             ta1_z_zzz_yzz_1,   \
                             ta1_z_zzz_yzzz_0,  \
                             ta1_z_zzz_yzzz_1,  \
                             ta1_z_zzz_zzz_0,   \
                             ta1_z_zzz_zzz_1,   \
                             ta1_z_zzz_zzzz_0,  \
                             ta1_z_zzz_zzzz_1,  \
                             ta1_z_zzzz_xxxx_0, \
                             ta1_z_zzzz_xxxy_0, \
                             ta1_z_zzzz_xxxz_0, \
                             ta1_z_zzzz_xxyy_0, \
                             ta1_z_zzzz_xxyz_0, \
                             ta1_z_zzzz_xxzz_0, \
                             ta1_z_zzzz_xyyy_0, \
                             ta1_z_zzzz_xyyz_0, \
                             ta1_z_zzzz_xyzz_0, \
                             ta1_z_zzzz_xzzz_0, \
                             ta1_z_zzzz_yyyy_0, \
                             ta1_z_zzzz_yyyz_0, \
                             ta1_z_zzzz_yyzz_0, \
                             ta1_z_zzzz_yzzz_0, \
                             ta1_z_zzzz_zzzz_0, \
                             ta_zzz_xxxx_1,     \
                             ta_zzz_xxxy_1,     \
                             ta_zzz_xxxz_1,     \
                             ta_zzz_xxyy_1,     \
                             ta_zzz_xxyz_1,     \
                             ta_zzz_xxzz_1,     \
                             ta_zzz_xyyy_1,     \
                             ta_zzz_xyyz_1,     \
                             ta_zzz_xyzz_1,     \
                             ta_zzz_xzzz_1,     \
                             ta_zzz_yyyy_1,     \
                             ta_zzz_yyyz_1,     \
                             ta_zzz_yyzz_1,     \
                             ta_zzz_yzzz_1,     \
                             ta_zzz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzz_xxxx_0[i] = 3.0 * ta1_z_zz_xxxx_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxx_1[i] * fe_0 + ta_zzz_xxxx_1[i] + ta1_z_zzz_xxxx_0[i] * pa_z[i] -
                               ta1_z_zzz_xxxx_1[i] * pc_z[i];

        ta1_z_zzzz_xxxy_0[i] = 3.0 * ta1_z_zz_xxxy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxy_1[i] * fe_0 + ta_zzz_xxxy_1[i] + ta1_z_zzz_xxxy_0[i] * pa_z[i] -
                               ta1_z_zzz_xxxy_1[i] * pc_z[i];

        ta1_z_zzzz_xxxz_0[i] = 3.0 * ta1_z_zz_xxxz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxz_1[i] * fe_0 + ta1_z_zzz_xxx_0[i] * fe_0 -
                               ta1_z_zzz_xxx_1[i] * fe_0 + ta_zzz_xxxz_1[i] + ta1_z_zzz_xxxz_0[i] * pa_z[i] - ta1_z_zzz_xxxz_1[i] * pc_z[i];

        ta1_z_zzzz_xxyy_0[i] = 3.0 * ta1_z_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyy_1[i] * fe_0 + ta_zzz_xxyy_1[i] + ta1_z_zzz_xxyy_0[i] * pa_z[i] -
                               ta1_z_zzz_xxyy_1[i] * pc_z[i];

        ta1_z_zzzz_xxyz_0[i] = 3.0 * ta1_z_zz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyz_1[i] * fe_0 + ta1_z_zzz_xxy_0[i] * fe_0 -
                               ta1_z_zzz_xxy_1[i] * fe_0 + ta_zzz_xxyz_1[i] + ta1_z_zzz_xxyz_0[i] * pa_z[i] - ta1_z_zzz_xxyz_1[i] * pc_z[i];

        ta1_z_zzzz_xxzz_0[i] = 3.0 * ta1_z_zz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxz_0[i] * fe_0 -
                               2.0 * ta1_z_zzz_xxz_1[i] * fe_0 + ta_zzz_xxzz_1[i] + ta1_z_zzz_xxzz_0[i] * pa_z[i] - ta1_z_zzz_xxzz_1[i] * pc_z[i];

        ta1_z_zzzz_xyyy_0[i] = 3.0 * ta1_z_zz_xyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyy_1[i] * fe_0 + ta_zzz_xyyy_1[i] + ta1_z_zzz_xyyy_0[i] * pa_z[i] -
                               ta1_z_zzz_xyyy_1[i] * pc_z[i];

        ta1_z_zzzz_xyyz_0[i] = 3.0 * ta1_z_zz_xyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyz_1[i] * fe_0 + ta1_z_zzz_xyy_0[i] * fe_0 -
                               ta1_z_zzz_xyy_1[i] * fe_0 + ta_zzz_xyyz_1[i] + ta1_z_zzz_xyyz_0[i] * pa_z[i] - ta1_z_zzz_xyyz_1[i] * pc_z[i];

        ta1_z_zzzz_xyzz_0[i] = 3.0 * ta1_z_zz_xyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyz_0[i] * fe_0 -
                               2.0 * ta1_z_zzz_xyz_1[i] * fe_0 + ta_zzz_xyzz_1[i] + ta1_z_zzz_xyzz_0[i] * pa_z[i] - ta1_z_zzz_xyzz_1[i] * pc_z[i];

        ta1_z_zzzz_xzzz_0[i] = 3.0 * ta1_z_zz_xzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_xzz_0[i] * fe_0 -
                               3.0 * ta1_z_zzz_xzz_1[i] * fe_0 + ta_zzz_xzzz_1[i] + ta1_z_zzz_xzzz_0[i] * pa_z[i] - ta1_z_zzz_xzzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyyy_0[i] = 3.0 * ta1_z_zz_yyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyy_1[i] * fe_0 + ta_zzz_yyyy_1[i] + ta1_z_zzz_yyyy_0[i] * pa_z[i] -
                               ta1_z_zzz_yyyy_1[i] * pc_z[i];

        ta1_z_zzzz_yyyz_0[i] = 3.0 * ta1_z_zz_yyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyz_1[i] * fe_0 + ta1_z_zzz_yyy_0[i] * fe_0 -
                               ta1_z_zzz_yyy_1[i] * fe_0 + ta_zzz_yyyz_1[i] + ta1_z_zzz_yyyz_0[i] * pa_z[i] - ta1_z_zzz_yyyz_1[i] * pc_z[i];

        ta1_z_zzzz_yyzz_0[i] = 3.0 * ta1_z_zz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyz_0[i] * fe_0 -
                               2.0 * ta1_z_zzz_yyz_1[i] * fe_0 + ta_zzz_yyzz_1[i] + ta1_z_zzz_yyzz_0[i] * pa_z[i] - ta1_z_zzz_yyzz_1[i] * pc_z[i];

        ta1_z_zzzz_yzzz_0[i] = 3.0 * ta1_z_zz_yzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_yzz_0[i] * fe_0 -
                               3.0 * ta1_z_zzz_yzz_1[i] * fe_0 + ta_zzz_yzzz_1[i] + ta1_z_zzz_yzzz_0[i] * pa_z[i] - ta1_z_zzz_yzzz_1[i] * pc_z[i];

        ta1_z_zzzz_zzzz_0[i] = 3.0 * ta1_z_zz_zzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_zzzz_1[i] * fe_0 + 4.0 * ta1_z_zzz_zzz_0[i] * fe_0 -
                               4.0 * ta1_z_zzz_zzz_1[i] * fe_0 + ta_zzz_zzzz_1[i] + ta1_z_zzz_zzzz_0[i] * pa_z[i] - ta1_z_zzz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
