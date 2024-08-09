#include "NuclearPotentialGeom020PrimRecFG.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_fg(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_fg,
                                        const size_t idx_npot_geom_020_0_pg,
                                        const size_t idx_npot_geom_020_1_pg,
                                        const size_t idx_npot_geom_020_0_df,
                                        const size_t idx_npot_geom_020_1_df,
                                        const size_t idx_npot_geom_010_1_dg,
                                        const size_t idx_npot_geom_020_0_dg,
                                        const size_t idx_npot_geom_020_1_dg,
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

    // Set up components of auxiliary buffer : PG

    auto ta2_xx_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg);

    auto ta2_xx_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 1);

    auto ta2_xx_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 2);

    auto ta2_xx_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 3);

    auto ta2_xx_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 4);

    auto ta2_xx_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 5);

    auto ta2_xx_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 6);

    auto ta2_xx_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 7);

    auto ta2_xx_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 8);

    auto ta2_xx_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 9);

    auto ta2_xx_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 10);

    auto ta2_xx_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 11);

    auto ta2_xx_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 12);

    auto ta2_xx_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 13);

    auto ta2_xx_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 14);

    auto ta2_xx_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 15);

    auto ta2_xx_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 16);

    auto ta2_xx_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 17);

    auto ta2_xx_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 18);

    auto ta2_xx_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 19);

    auto ta2_xx_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 20);

    auto ta2_xx_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 21);

    auto ta2_xx_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 22);

    auto ta2_xx_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 23);

    auto ta2_xx_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 24);

    auto ta2_xx_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 25);

    auto ta2_xx_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 26);

    auto ta2_xx_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 27);

    auto ta2_xx_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 28);

    auto ta2_xx_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 29);

    auto ta2_xx_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 30);

    auto ta2_xx_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 31);

    auto ta2_xx_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 32);

    auto ta2_xx_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 33);

    auto ta2_xx_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 34);

    auto ta2_xx_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 35);

    auto ta2_xx_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 36);

    auto ta2_xx_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 37);

    auto ta2_xx_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 38);

    auto ta2_xx_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 39);

    auto ta2_xx_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 40);

    auto ta2_xx_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 41);

    auto ta2_xx_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 42);

    auto ta2_xx_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 43);

    auto ta2_xx_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 44);

    auto ta2_xy_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 45);

    auto ta2_xy_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 46);

    auto ta2_xy_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 47);

    auto ta2_xy_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 48);

    auto ta2_xy_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 49);

    auto ta2_xy_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 50);

    auto ta2_xy_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 51);

    auto ta2_xy_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 52);

    auto ta2_xy_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 53);

    auto ta2_xy_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 54);

    auto ta2_xy_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 55);

    auto ta2_xy_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 56);

    auto ta2_xy_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 57);

    auto ta2_xy_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 58);

    auto ta2_xy_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 59);

    auto ta2_xy_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 60);

    auto ta2_xy_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 61);

    auto ta2_xy_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 62);

    auto ta2_xy_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 63);

    auto ta2_xy_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 64);

    auto ta2_xy_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 65);

    auto ta2_xy_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 66);

    auto ta2_xy_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 67);

    auto ta2_xy_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 68);

    auto ta2_xy_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 69);

    auto ta2_xy_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 70);

    auto ta2_xy_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 71);

    auto ta2_xy_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 72);

    auto ta2_xy_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 73);

    auto ta2_xy_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 74);

    auto ta2_xy_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 75);

    auto ta2_xy_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 76);

    auto ta2_xy_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 77);

    auto ta2_xy_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 78);

    auto ta2_xy_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 79);

    auto ta2_xy_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 80);

    auto ta2_xy_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 81);

    auto ta2_xy_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 82);

    auto ta2_xy_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 83);

    auto ta2_xy_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 84);

    auto ta2_xy_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 85);

    auto ta2_xy_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 86);

    auto ta2_xy_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 87);

    auto ta2_xy_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 88);

    auto ta2_xy_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 89);

    auto ta2_xz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 90);

    auto ta2_xz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 91);

    auto ta2_xz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 92);

    auto ta2_xz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 93);

    auto ta2_xz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 94);

    auto ta2_xz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 95);

    auto ta2_xz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 96);

    auto ta2_xz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 97);

    auto ta2_xz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 98);

    auto ta2_xz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 99);

    auto ta2_xz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 100);

    auto ta2_xz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 101);

    auto ta2_xz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 102);

    auto ta2_xz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 103);

    auto ta2_xz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 104);

    auto ta2_xz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 105);

    auto ta2_xz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 106);

    auto ta2_xz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 107);

    auto ta2_xz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 108);

    auto ta2_xz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 109);

    auto ta2_xz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 110);

    auto ta2_xz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 111);

    auto ta2_xz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 112);

    auto ta2_xz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 113);

    auto ta2_xz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 114);

    auto ta2_xz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 115);

    auto ta2_xz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 116);

    auto ta2_xz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 117);

    auto ta2_xz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 118);

    auto ta2_xz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 119);

    auto ta2_xz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 120);

    auto ta2_xz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 121);

    auto ta2_xz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 122);

    auto ta2_xz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 123);

    auto ta2_xz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 124);

    auto ta2_xz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 125);

    auto ta2_xz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 126);

    auto ta2_xz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 127);

    auto ta2_xz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 128);

    auto ta2_xz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 129);

    auto ta2_xz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 130);

    auto ta2_xz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 131);

    auto ta2_xz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 132);

    auto ta2_xz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 133);

    auto ta2_xz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 134);

    auto ta2_yy_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 135);

    auto ta2_yy_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 136);

    auto ta2_yy_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 137);

    auto ta2_yy_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 138);

    auto ta2_yy_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 139);

    auto ta2_yy_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 140);

    auto ta2_yy_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 141);

    auto ta2_yy_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 142);

    auto ta2_yy_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 143);

    auto ta2_yy_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 144);

    auto ta2_yy_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 145);

    auto ta2_yy_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 146);

    auto ta2_yy_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 147);

    auto ta2_yy_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 148);

    auto ta2_yy_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 149);

    auto ta2_yy_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 150);

    auto ta2_yy_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 151);

    auto ta2_yy_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 152);

    auto ta2_yy_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 153);

    auto ta2_yy_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 154);

    auto ta2_yy_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 155);

    auto ta2_yy_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 156);

    auto ta2_yy_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 157);

    auto ta2_yy_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 158);

    auto ta2_yy_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 159);

    auto ta2_yy_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 160);

    auto ta2_yy_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 161);

    auto ta2_yy_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 162);

    auto ta2_yy_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 163);

    auto ta2_yy_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 164);

    auto ta2_yy_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 165);

    auto ta2_yy_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 166);

    auto ta2_yy_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 167);

    auto ta2_yy_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 168);

    auto ta2_yy_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 169);

    auto ta2_yy_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 170);

    auto ta2_yy_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 171);

    auto ta2_yy_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 172);

    auto ta2_yy_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 173);

    auto ta2_yy_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 174);

    auto ta2_yy_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 175);

    auto ta2_yy_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 176);

    auto ta2_yy_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 177);

    auto ta2_yy_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 178);

    auto ta2_yy_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 179);

    auto ta2_yz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 180);

    auto ta2_yz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 181);

    auto ta2_yz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 182);

    auto ta2_yz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 183);

    auto ta2_yz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 184);

    auto ta2_yz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 185);

    auto ta2_yz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 186);

    auto ta2_yz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 187);

    auto ta2_yz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 188);

    auto ta2_yz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 189);

    auto ta2_yz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 190);

    auto ta2_yz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 191);

    auto ta2_yz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 192);

    auto ta2_yz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 193);

    auto ta2_yz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 194);

    auto ta2_yz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 195);

    auto ta2_yz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 196);

    auto ta2_yz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 197);

    auto ta2_yz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 198);

    auto ta2_yz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 199);

    auto ta2_yz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 200);

    auto ta2_yz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 201);

    auto ta2_yz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 202);

    auto ta2_yz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 203);

    auto ta2_yz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 204);

    auto ta2_yz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 205);

    auto ta2_yz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 206);

    auto ta2_yz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 207);

    auto ta2_yz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 208);

    auto ta2_yz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 209);

    auto ta2_yz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 210);

    auto ta2_yz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 211);

    auto ta2_yz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 212);

    auto ta2_yz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 213);

    auto ta2_yz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 214);

    auto ta2_yz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 215);

    auto ta2_yz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 216);

    auto ta2_yz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 217);

    auto ta2_yz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 218);

    auto ta2_yz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 219);

    auto ta2_yz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 220);

    auto ta2_yz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 221);

    auto ta2_yz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 222);

    auto ta2_yz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 223);

    auto ta2_yz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 224);

    auto ta2_zz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 225);

    auto ta2_zz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 226);

    auto ta2_zz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 227);

    auto ta2_zz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 228);

    auto ta2_zz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 229);

    auto ta2_zz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 230);

    auto ta2_zz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 231);

    auto ta2_zz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 232);

    auto ta2_zz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 233);

    auto ta2_zz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 234);

    auto ta2_zz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 235);

    auto ta2_zz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 236);

    auto ta2_zz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 237);

    auto ta2_zz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 238);

    auto ta2_zz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 239);

    auto ta2_zz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 240);

    auto ta2_zz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 241);

    auto ta2_zz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 242);

    auto ta2_zz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 243);

    auto ta2_zz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 244);

    auto ta2_zz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 245);

    auto ta2_zz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 246);

    auto ta2_zz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 247);

    auto ta2_zz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 248);

    auto ta2_zz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 249);

    auto ta2_zz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 250);

    auto ta2_zz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 251);

    auto ta2_zz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 252);

    auto ta2_zz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 253);

    auto ta2_zz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 254);

    auto ta2_zz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 255);

    auto ta2_zz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 256);

    auto ta2_zz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 257);

    auto ta2_zz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 258);

    auto ta2_zz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 259);

    auto ta2_zz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 260);

    auto ta2_zz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 261);

    auto ta2_zz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 262);

    auto ta2_zz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 263);

    auto ta2_zz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 264);

    auto ta2_zz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 265);

    auto ta2_zz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 266);

    auto ta2_zz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 267);

    auto ta2_zz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 268);

    auto ta2_zz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 269);

    // Set up components of auxiliary buffer : PG

    auto ta2_xx_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg);

    auto ta2_xx_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 1);

    auto ta2_xx_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 2);

    auto ta2_xx_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 3);

    auto ta2_xx_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 4);

    auto ta2_xx_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 5);

    auto ta2_xx_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 6);

    auto ta2_xx_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 7);

    auto ta2_xx_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 8);

    auto ta2_xx_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 9);

    auto ta2_xx_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 10);

    auto ta2_xx_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 11);

    auto ta2_xx_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 12);

    auto ta2_xx_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 13);

    auto ta2_xx_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 14);

    auto ta2_xx_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 15);

    auto ta2_xx_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 16);

    auto ta2_xx_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 17);

    auto ta2_xx_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 18);

    auto ta2_xx_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 19);

    auto ta2_xx_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 20);

    auto ta2_xx_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 21);

    auto ta2_xx_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 22);

    auto ta2_xx_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 23);

    auto ta2_xx_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 24);

    auto ta2_xx_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 25);

    auto ta2_xx_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 26);

    auto ta2_xx_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 27);

    auto ta2_xx_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 28);

    auto ta2_xx_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 29);

    auto ta2_xx_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 30);

    auto ta2_xx_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 31);

    auto ta2_xx_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 32);

    auto ta2_xx_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 33);

    auto ta2_xx_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 34);

    auto ta2_xx_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 35);

    auto ta2_xx_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 36);

    auto ta2_xx_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 37);

    auto ta2_xx_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 38);

    auto ta2_xx_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 39);

    auto ta2_xx_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 40);

    auto ta2_xx_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 41);

    auto ta2_xx_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 42);

    auto ta2_xx_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 43);

    auto ta2_xx_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 44);

    auto ta2_xy_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 45);

    auto ta2_xy_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 46);

    auto ta2_xy_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 47);

    auto ta2_xy_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 48);

    auto ta2_xy_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 49);

    auto ta2_xy_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 50);

    auto ta2_xy_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 51);

    auto ta2_xy_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 52);

    auto ta2_xy_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 53);

    auto ta2_xy_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 54);

    auto ta2_xy_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 55);

    auto ta2_xy_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 56);

    auto ta2_xy_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 57);

    auto ta2_xy_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 58);

    auto ta2_xy_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 59);

    auto ta2_xy_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 60);

    auto ta2_xy_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 61);

    auto ta2_xy_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 62);

    auto ta2_xy_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 63);

    auto ta2_xy_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 64);

    auto ta2_xy_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 65);

    auto ta2_xy_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 66);

    auto ta2_xy_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 67);

    auto ta2_xy_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 68);

    auto ta2_xy_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 69);

    auto ta2_xy_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 70);

    auto ta2_xy_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 71);

    auto ta2_xy_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 72);

    auto ta2_xy_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 73);

    auto ta2_xy_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 74);

    auto ta2_xy_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 75);

    auto ta2_xy_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 76);

    auto ta2_xy_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 77);

    auto ta2_xy_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 78);

    auto ta2_xy_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 79);

    auto ta2_xy_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 80);

    auto ta2_xy_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 81);

    auto ta2_xy_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 82);

    auto ta2_xy_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 83);

    auto ta2_xy_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 84);

    auto ta2_xy_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 85);

    auto ta2_xy_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 86);

    auto ta2_xy_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 87);

    auto ta2_xy_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 88);

    auto ta2_xy_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 89);

    auto ta2_xz_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 90);

    auto ta2_xz_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 91);

    auto ta2_xz_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 92);

    auto ta2_xz_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 93);

    auto ta2_xz_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 94);

    auto ta2_xz_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 95);

    auto ta2_xz_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 96);

    auto ta2_xz_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 97);

    auto ta2_xz_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 98);

    auto ta2_xz_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 99);

    auto ta2_xz_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 100);

    auto ta2_xz_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 101);

    auto ta2_xz_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 102);

    auto ta2_xz_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 103);

    auto ta2_xz_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 104);

    auto ta2_xz_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 105);

    auto ta2_xz_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 106);

    auto ta2_xz_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 107);

    auto ta2_xz_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 108);

    auto ta2_xz_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 109);

    auto ta2_xz_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 110);

    auto ta2_xz_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 111);

    auto ta2_xz_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 112);

    auto ta2_xz_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 113);

    auto ta2_xz_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 114);

    auto ta2_xz_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 115);

    auto ta2_xz_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 116);

    auto ta2_xz_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 117);

    auto ta2_xz_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 118);

    auto ta2_xz_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 119);

    auto ta2_xz_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 120);

    auto ta2_xz_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 121);

    auto ta2_xz_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 122);

    auto ta2_xz_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 123);

    auto ta2_xz_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 124);

    auto ta2_xz_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 125);

    auto ta2_xz_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 126);

    auto ta2_xz_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 127);

    auto ta2_xz_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 128);

    auto ta2_xz_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 129);

    auto ta2_xz_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 130);

    auto ta2_xz_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 131);

    auto ta2_xz_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 132);

    auto ta2_xz_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 133);

    auto ta2_xz_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 134);

    auto ta2_yy_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 135);

    auto ta2_yy_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 136);

    auto ta2_yy_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 137);

    auto ta2_yy_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 138);

    auto ta2_yy_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 139);

    auto ta2_yy_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 140);

    auto ta2_yy_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 141);

    auto ta2_yy_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 142);

    auto ta2_yy_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 143);

    auto ta2_yy_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 144);

    auto ta2_yy_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 145);

    auto ta2_yy_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 146);

    auto ta2_yy_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 147);

    auto ta2_yy_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 148);

    auto ta2_yy_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 149);

    auto ta2_yy_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 150);

    auto ta2_yy_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 151);

    auto ta2_yy_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 152);

    auto ta2_yy_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 153);

    auto ta2_yy_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 154);

    auto ta2_yy_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 155);

    auto ta2_yy_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 156);

    auto ta2_yy_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 157);

    auto ta2_yy_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 158);

    auto ta2_yy_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 159);

    auto ta2_yy_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 160);

    auto ta2_yy_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 161);

    auto ta2_yy_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 162);

    auto ta2_yy_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 163);

    auto ta2_yy_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 164);

    auto ta2_yy_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 165);

    auto ta2_yy_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 166);

    auto ta2_yy_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 167);

    auto ta2_yy_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 168);

    auto ta2_yy_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 169);

    auto ta2_yy_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 170);

    auto ta2_yy_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 171);

    auto ta2_yy_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 172);

    auto ta2_yy_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 173);

    auto ta2_yy_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 174);

    auto ta2_yy_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 175);

    auto ta2_yy_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 176);

    auto ta2_yy_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 177);

    auto ta2_yy_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 178);

    auto ta2_yy_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 179);

    auto ta2_yz_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 180);

    auto ta2_yz_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 181);

    auto ta2_yz_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 182);

    auto ta2_yz_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 183);

    auto ta2_yz_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 184);

    auto ta2_yz_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 185);

    auto ta2_yz_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 186);

    auto ta2_yz_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 187);

    auto ta2_yz_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 188);

    auto ta2_yz_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 189);

    auto ta2_yz_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 190);

    auto ta2_yz_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 191);

    auto ta2_yz_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 192);

    auto ta2_yz_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 193);

    auto ta2_yz_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 194);

    auto ta2_yz_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 195);

    auto ta2_yz_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 196);

    auto ta2_yz_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 197);

    auto ta2_yz_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 198);

    auto ta2_yz_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 199);

    auto ta2_yz_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 200);

    auto ta2_yz_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 201);

    auto ta2_yz_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 202);

    auto ta2_yz_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 203);

    auto ta2_yz_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 204);

    auto ta2_yz_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 205);

    auto ta2_yz_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 206);

    auto ta2_yz_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 207);

    auto ta2_yz_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 208);

    auto ta2_yz_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 209);

    auto ta2_yz_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 210);

    auto ta2_yz_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 211);

    auto ta2_yz_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 212);

    auto ta2_yz_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 213);

    auto ta2_yz_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 214);

    auto ta2_yz_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 215);

    auto ta2_yz_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 216);

    auto ta2_yz_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 217);

    auto ta2_yz_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 218);

    auto ta2_yz_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 219);

    auto ta2_yz_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 220);

    auto ta2_yz_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 221);

    auto ta2_yz_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 222);

    auto ta2_yz_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 223);

    auto ta2_yz_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 224);

    auto ta2_zz_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 225);

    auto ta2_zz_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 226);

    auto ta2_zz_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 227);

    auto ta2_zz_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 228);

    auto ta2_zz_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 229);

    auto ta2_zz_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 230);

    auto ta2_zz_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 231);

    auto ta2_zz_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 232);

    auto ta2_zz_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 233);

    auto ta2_zz_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 234);

    auto ta2_zz_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 235);

    auto ta2_zz_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 236);

    auto ta2_zz_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 237);

    auto ta2_zz_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 238);

    auto ta2_zz_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 239);

    auto ta2_zz_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 240);

    auto ta2_zz_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 241);

    auto ta2_zz_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 242);

    auto ta2_zz_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 243);

    auto ta2_zz_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 244);

    auto ta2_zz_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 245);

    auto ta2_zz_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 246);

    auto ta2_zz_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 247);

    auto ta2_zz_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 248);

    auto ta2_zz_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 249);

    auto ta2_zz_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 250);

    auto ta2_zz_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 251);

    auto ta2_zz_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 252);

    auto ta2_zz_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 253);

    auto ta2_zz_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 254);

    auto ta2_zz_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 255);

    auto ta2_zz_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 256);

    auto ta2_zz_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 257);

    auto ta2_zz_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 258);

    auto ta2_zz_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 259);

    auto ta2_zz_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 260);

    auto ta2_zz_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 261);

    auto ta2_zz_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 262);

    auto ta2_zz_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 263);

    auto ta2_zz_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 264);

    auto ta2_zz_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 265);

    auto ta2_zz_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 266);

    auto ta2_zz_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 267);

    auto ta2_zz_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 268);

    auto ta2_zz_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 269);

    // Set up components of auxiliary buffer : DF

    auto ta2_xx_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df);

    auto ta2_xx_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 1);

    auto ta2_xx_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 2);

    auto ta2_xx_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 3);

    auto ta2_xx_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 4);

    auto ta2_xx_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 5);

    auto ta2_xx_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 6);

    auto ta2_xx_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 7);

    auto ta2_xx_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 8);

    auto ta2_xx_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 9);

    auto ta2_xx_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 22);

    auto ta2_xx_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 24);

    auto ta2_xx_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 25);

    auto ta2_xx_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 30);

    auto ta2_xx_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 31);

    auto ta2_xx_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 32);

    auto ta2_xx_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 33);

    auto ta2_xx_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 34);

    auto ta2_xx_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 35);

    auto ta2_xx_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 36);

    auto ta2_xx_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 37);

    auto ta2_xx_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 38);

    auto ta2_xx_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 39);

    auto ta2_xx_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 50);

    auto ta2_xx_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 51);

    auto ta2_xx_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 52);

    auto ta2_xx_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 53);

    auto ta2_xx_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 54);

    auto ta2_xx_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 55);

    auto ta2_xx_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 56);

    auto ta2_xx_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 57);

    auto ta2_xx_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 58);

    auto ta2_xx_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 59);

    auto ta2_xy_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 60);

    auto ta2_xy_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 61);

    auto ta2_xy_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 62);

    auto ta2_xy_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 63);

    auto ta2_xy_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 64);

    auto ta2_xy_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 65);

    auto ta2_xy_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 66);

    auto ta2_xy_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 67);

    auto ta2_xy_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 68);

    auto ta2_xy_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 69);

    auto ta2_xy_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 71);

    auto ta2_xy_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 73);

    auto ta2_xy_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 74);

    auto ta2_xy_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 90);

    auto ta2_xy_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 91);

    auto ta2_xy_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 92);

    auto ta2_xy_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 93);

    auto ta2_xy_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 94);

    auto ta2_xy_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 95);

    auto ta2_xy_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 96);

    auto ta2_xy_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 97);

    auto ta2_xy_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 98);

    auto ta2_xy_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 99);

    auto ta2_xy_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 110);

    auto ta2_xy_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 111);

    auto ta2_xy_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 112);

    auto ta2_xy_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 113);

    auto ta2_xy_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 114);

    auto ta2_xy_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 115);

    auto ta2_xy_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 116);

    auto ta2_xy_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 117);

    auto ta2_xy_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 118);

    auto ta2_xy_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 119);

    auto ta2_xz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 120);

    auto ta2_xz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 121);

    auto ta2_xz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 122);

    auto ta2_xz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 123);

    auto ta2_xz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 124);

    auto ta2_xz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 125);

    auto ta2_xz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 126);

    auto ta2_xz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 127);

    auto ta2_xz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 128);

    auto ta2_xz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 129);

    auto ta2_xz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 142);

    auto ta2_xz_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 144);

    auto ta2_xz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 145);

    auto ta2_xz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 150);

    auto ta2_xz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 151);

    auto ta2_xz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 152);

    auto ta2_xz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 153);

    auto ta2_xz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 154);

    auto ta2_xz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 155);

    auto ta2_xz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 156);

    auto ta2_xz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 157);

    auto ta2_xz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 158);

    auto ta2_xz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 159);

    auto ta2_xz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 170);

    auto ta2_xz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 171);

    auto ta2_xz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 172);

    auto ta2_xz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 173);

    auto ta2_xz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 174);

    auto ta2_xz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 175);

    auto ta2_xz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 176);

    auto ta2_xz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 177);

    auto ta2_xz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 178);

    auto ta2_xz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 179);

    auto ta2_yy_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 180);

    auto ta2_yy_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 181);

    auto ta2_yy_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 182);

    auto ta2_yy_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 183);

    auto ta2_yy_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 184);

    auto ta2_yy_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 185);

    auto ta2_yy_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 186);

    auto ta2_yy_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 187);

    auto ta2_yy_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 188);

    auto ta2_yy_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 189);

    auto ta2_yy_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 210);

    auto ta2_yy_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 211);

    auto ta2_yy_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 212);

    auto ta2_yy_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 213);

    auto ta2_yy_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 214);

    auto ta2_yy_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 215);

    auto ta2_yy_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 216);

    auto ta2_yy_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 217);

    auto ta2_yy_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 218);

    auto ta2_yy_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 219);

    auto ta2_yy_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 224);

    auto ta2_yy_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 227);

    auto ta2_yy_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 228);

    auto ta2_yy_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 230);

    auto ta2_yy_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 231);

    auto ta2_yy_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 232);

    auto ta2_yy_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 233);

    auto ta2_yy_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 234);

    auto ta2_yy_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 235);

    auto ta2_yy_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 236);

    auto ta2_yy_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 237);

    auto ta2_yy_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 238);

    auto ta2_yy_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 239);

    auto ta2_yz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 240);

    auto ta2_yz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 241);

    auto ta2_yz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 242);

    auto ta2_yz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 243);

    auto ta2_yz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 244);

    auto ta2_yz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 245);

    auto ta2_yz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 246);

    auto ta2_yz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 247);

    auto ta2_yz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 248);

    auto ta2_yz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 249);

    auto ta2_yz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 270);

    auto ta2_yz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 271);

    auto ta2_yz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 272);

    auto ta2_yz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 273);

    auto ta2_yz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 274);

    auto ta2_yz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 275);

    auto ta2_yz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 276);

    auto ta2_yz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 277);

    auto ta2_yz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 278);

    auto ta2_yz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 279);

    auto ta2_yz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 284);

    auto ta2_yz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 287);

    auto ta2_yz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 288);

    auto ta2_yz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 290);

    auto ta2_yz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 291);

    auto ta2_yz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 292);

    auto ta2_yz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 293);

    auto ta2_yz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 294);

    auto ta2_yz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 295);

    auto ta2_yz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 296);

    auto ta2_yz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 297);

    auto ta2_yz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 298);

    auto ta2_yz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 299);

    auto ta2_zz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 300);

    auto ta2_zz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 301);

    auto ta2_zz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 302);

    auto ta2_zz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 303);

    auto ta2_zz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 304);

    auto ta2_zz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 305);

    auto ta2_zz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 306);

    auto ta2_zz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 307);

    auto ta2_zz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 308);

    auto ta2_zz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 309);

    auto ta2_zz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 330);

    auto ta2_zz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 331);

    auto ta2_zz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 332);

    auto ta2_zz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 333);

    auto ta2_zz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 334);

    auto ta2_zz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 335);

    auto ta2_zz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 336);

    auto ta2_zz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 337);

    auto ta2_zz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 338);

    auto ta2_zz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 339);

    auto ta2_zz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 344);

    auto ta2_zz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 347);

    auto ta2_zz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 348);

    auto ta2_zz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 350);

    auto ta2_zz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 351);

    auto ta2_zz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 352);

    auto ta2_zz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 353);

    auto ta2_zz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 354);

    auto ta2_zz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 355);

    auto ta2_zz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 356);

    auto ta2_zz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 357);

    auto ta2_zz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 358);

    auto ta2_zz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 359);

    // Set up components of auxiliary buffer : DF

    auto ta2_xx_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df);

    auto ta2_xx_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 1);

    auto ta2_xx_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 2);

    auto ta2_xx_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 3);

    auto ta2_xx_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 4);

    auto ta2_xx_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 5);

    auto ta2_xx_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 6);

    auto ta2_xx_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 7);

    auto ta2_xx_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 8);

    auto ta2_xx_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 9);

    auto ta2_xx_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 22);

    auto ta2_xx_xz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 24);

    auto ta2_xx_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 25);

    auto ta2_xx_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 30);

    auto ta2_xx_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 31);

    auto ta2_xx_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 32);

    auto ta2_xx_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 33);

    auto ta2_xx_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 34);

    auto ta2_xx_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 35);

    auto ta2_xx_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 36);

    auto ta2_xx_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 37);

    auto ta2_xx_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 38);

    auto ta2_xx_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 39);

    auto ta2_xx_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 50);

    auto ta2_xx_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 51);

    auto ta2_xx_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 52);

    auto ta2_xx_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 53);

    auto ta2_xx_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 54);

    auto ta2_xx_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 55);

    auto ta2_xx_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 56);

    auto ta2_xx_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 57);

    auto ta2_xx_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 58);

    auto ta2_xx_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 59);

    auto ta2_xy_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 60);

    auto ta2_xy_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 61);

    auto ta2_xy_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 62);

    auto ta2_xy_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 63);

    auto ta2_xy_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 64);

    auto ta2_xy_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 65);

    auto ta2_xy_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 66);

    auto ta2_xy_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 67);

    auto ta2_xy_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 68);

    auto ta2_xy_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 69);

    auto ta2_xy_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 71);

    auto ta2_xy_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 73);

    auto ta2_xy_xy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 74);

    auto ta2_xy_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 90);

    auto ta2_xy_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 91);

    auto ta2_xy_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 92);

    auto ta2_xy_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 93);

    auto ta2_xy_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 94);

    auto ta2_xy_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 95);

    auto ta2_xy_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 96);

    auto ta2_xy_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 97);

    auto ta2_xy_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 98);

    auto ta2_xy_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 99);

    auto ta2_xy_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 110);

    auto ta2_xy_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 111);

    auto ta2_xy_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 112);

    auto ta2_xy_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 113);

    auto ta2_xy_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 114);

    auto ta2_xy_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 115);

    auto ta2_xy_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 116);

    auto ta2_xy_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 117);

    auto ta2_xy_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 118);

    auto ta2_xy_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 119);

    auto ta2_xz_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 120);

    auto ta2_xz_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 121);

    auto ta2_xz_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 122);

    auto ta2_xz_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 123);

    auto ta2_xz_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 124);

    auto ta2_xz_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 125);

    auto ta2_xz_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 126);

    auto ta2_xz_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 127);

    auto ta2_xz_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 128);

    auto ta2_xz_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 129);

    auto ta2_xz_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 142);

    auto ta2_xz_xz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 144);

    auto ta2_xz_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 145);

    auto ta2_xz_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 150);

    auto ta2_xz_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 151);

    auto ta2_xz_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 152);

    auto ta2_xz_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 153);

    auto ta2_xz_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 154);

    auto ta2_xz_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 155);

    auto ta2_xz_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 156);

    auto ta2_xz_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 157);

    auto ta2_xz_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 158);

    auto ta2_xz_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 159);

    auto ta2_xz_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 170);

    auto ta2_xz_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 171);

    auto ta2_xz_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 172);

    auto ta2_xz_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 173);

    auto ta2_xz_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 174);

    auto ta2_xz_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 175);

    auto ta2_xz_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 176);

    auto ta2_xz_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 177);

    auto ta2_xz_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 178);

    auto ta2_xz_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 179);

    auto ta2_yy_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 180);

    auto ta2_yy_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 181);

    auto ta2_yy_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 182);

    auto ta2_yy_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 183);

    auto ta2_yy_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 184);

    auto ta2_yy_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 185);

    auto ta2_yy_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 186);

    auto ta2_yy_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 187);

    auto ta2_yy_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 188);

    auto ta2_yy_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 189);

    auto ta2_yy_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 210);

    auto ta2_yy_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 211);

    auto ta2_yy_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 212);

    auto ta2_yy_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 213);

    auto ta2_yy_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 214);

    auto ta2_yy_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 215);

    auto ta2_yy_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 216);

    auto ta2_yy_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 217);

    auto ta2_yy_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 218);

    auto ta2_yy_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 219);

    auto ta2_yy_yz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 224);

    auto ta2_yy_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 227);

    auto ta2_yy_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 228);

    auto ta2_yy_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 230);

    auto ta2_yy_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 231);

    auto ta2_yy_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 232);

    auto ta2_yy_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 233);

    auto ta2_yy_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 234);

    auto ta2_yy_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 235);

    auto ta2_yy_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 236);

    auto ta2_yy_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 237);

    auto ta2_yy_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 238);

    auto ta2_yy_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 239);

    auto ta2_yz_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 240);

    auto ta2_yz_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 241);

    auto ta2_yz_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 242);

    auto ta2_yz_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 243);

    auto ta2_yz_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 244);

    auto ta2_yz_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 245);

    auto ta2_yz_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 246);

    auto ta2_yz_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 247);

    auto ta2_yz_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 248);

    auto ta2_yz_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 249);

    auto ta2_yz_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 270);

    auto ta2_yz_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 271);

    auto ta2_yz_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 272);

    auto ta2_yz_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 273);

    auto ta2_yz_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 274);

    auto ta2_yz_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 275);

    auto ta2_yz_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 276);

    auto ta2_yz_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 277);

    auto ta2_yz_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 278);

    auto ta2_yz_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 279);

    auto ta2_yz_yz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 284);

    auto ta2_yz_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 287);

    auto ta2_yz_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 288);

    auto ta2_yz_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 290);

    auto ta2_yz_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 291);

    auto ta2_yz_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 292);

    auto ta2_yz_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 293);

    auto ta2_yz_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 294);

    auto ta2_yz_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 295);

    auto ta2_yz_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 296);

    auto ta2_yz_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 297);

    auto ta2_yz_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 298);

    auto ta2_yz_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 299);

    auto ta2_zz_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 300);

    auto ta2_zz_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 301);

    auto ta2_zz_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 302);

    auto ta2_zz_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 303);

    auto ta2_zz_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 304);

    auto ta2_zz_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 305);

    auto ta2_zz_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 306);

    auto ta2_zz_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 307);

    auto ta2_zz_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 308);

    auto ta2_zz_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 309);

    auto ta2_zz_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 330);

    auto ta2_zz_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 331);

    auto ta2_zz_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 332);

    auto ta2_zz_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 333);

    auto ta2_zz_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 334);

    auto ta2_zz_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 335);

    auto ta2_zz_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 336);

    auto ta2_zz_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 337);

    auto ta2_zz_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 338);

    auto ta2_zz_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 339);

    auto ta2_zz_yz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 344);

    auto ta2_zz_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 347);

    auto ta2_zz_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 348);

    auto ta2_zz_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 350);

    auto ta2_zz_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 351);

    auto ta2_zz_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 352);

    auto ta2_zz_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 353);

    auto ta2_zz_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 354);

    auto ta2_zz_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 355);

    auto ta2_zz_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 356);

    auto ta2_zz_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 357);

    auto ta2_zz_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 358);

    auto ta2_zz_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 359);

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

    auto ta1_x_xy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 16);

    auto ta1_x_xy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 18);

    auto ta1_x_xy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 21);

    auto ta1_x_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 32);

    auto ta1_x_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 35);

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

    auto ta1_x_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 71);

    auto ta1_x_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 72);

    auto ta1_x_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 73);

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

    auto ta1_y_xy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 106);

    auto ta1_y_xy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 108);

    auto ta1_y_xy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 111);

    auto ta1_y_xy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 115);

    auto ta1_y_xy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 116);

    auto ta1_y_xy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 117);

    auto ta1_y_xy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 118);

    auto ta1_y_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 122);

    auto ta1_y_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 125);

    auto ta1_y_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 129);

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

    auto ta1_y_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 161);

    auto ta1_y_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 162);

    auto ta1_y_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 163);

    auto ta1_y_yz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 164);

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

    auto ta1_z_xy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 196);

    auto ta1_z_xy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 198);

    auto ta1_z_xy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 201);

    auto ta1_z_xz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 210);

    auto ta1_z_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 212);

    auto ta1_z_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 215);

    auto ta1_z_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 219);

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

    auto ta1_z_yz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 250);

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

    // Set up components of auxiliary buffer : DG

    auto ta2_xx_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg);

    auto ta2_xx_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 1);

    auto ta2_xx_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 2);

    auto ta2_xx_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 3);

    auto ta2_xx_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 4);

    auto ta2_xx_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 5);

    auto ta2_xx_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 6);

    auto ta2_xx_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 7);

    auto ta2_xx_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 8);

    auto ta2_xx_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 9);

    auto ta2_xx_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 10);

    auto ta2_xx_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 11);

    auto ta2_xx_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 12);

    auto ta2_xx_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 13);

    auto ta2_xx_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 14);

    auto ta2_xx_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 15);

    auto ta2_xx_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 16);

    auto ta2_xx_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 17);

    auto ta2_xx_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 18);

    auto ta2_xx_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 20);

    auto ta2_xx_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 21);

    auto ta2_xx_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 24);

    auto ta2_xx_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 25);

    auto ta2_xx_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 30);

    auto ta2_xx_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 31);

    auto ta2_xx_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 32);

    auto ta2_xx_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 33);

    auto ta2_xx_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 34);

    auto ta2_xx_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 35);

    auto ta2_xx_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 36);

    auto ta2_xx_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 37);

    auto ta2_xx_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 38);

    auto ta2_xx_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 39);

    auto ta2_xx_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 44);

    auto ta2_xx_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 45);

    auto ta2_xx_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 46);

    auto ta2_xx_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 47);

    auto ta2_xx_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 48);

    auto ta2_xx_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 49);

    auto ta2_xx_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 50);

    auto ta2_xx_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 51);

    auto ta2_xx_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 52);

    auto ta2_xx_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 53);

    auto ta2_xx_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 54);

    auto ta2_xx_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 55);

    auto ta2_xx_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 56);

    auto ta2_xx_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 57);

    auto ta2_xx_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 58);

    auto ta2_xx_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 59);

    auto ta2_xx_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 62);

    auto ta2_xx_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 65);

    auto ta2_xx_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 69);

    auto ta2_xx_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 71);

    auto ta2_xx_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 72);

    auto ta2_xx_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 73);

    auto ta2_xx_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 74);

    auto ta2_xx_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 75);

    auto ta2_xx_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 76);

    auto ta2_xx_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 77);

    auto ta2_xx_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 78);

    auto ta2_xx_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 79);

    auto ta2_xx_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 80);

    auto ta2_xx_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 81);

    auto ta2_xx_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 82);

    auto ta2_xx_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 83);

    auto ta2_xx_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 84);

    auto ta2_xx_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 85);

    auto ta2_xx_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 86);

    auto ta2_xx_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 87);

    auto ta2_xx_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 88);

    auto ta2_xx_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 89);

    auto ta2_xy_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 90);

    auto ta2_xy_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 91);

    auto ta2_xy_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 92);

    auto ta2_xy_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 93);

    auto ta2_xy_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 94);

    auto ta2_xy_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 95);

    auto ta2_xy_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 96);

    auto ta2_xy_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 97);

    auto ta2_xy_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 98);

    auto ta2_xy_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 99);

    auto ta2_xy_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 100);

    auto ta2_xy_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 101);

    auto ta2_xy_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 102);

    auto ta2_xy_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 103);

    auto ta2_xy_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 104);

    auto ta2_xy_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 105);

    auto ta2_xy_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 106);

    auto ta2_xy_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 108);

    auto ta2_xy_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 109);

    auto ta2_xy_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 111);

    auto ta2_xy_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 112);

    auto ta2_xy_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 113);

    auto ta2_xy_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 115);

    auto ta2_xy_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 116);

    auto ta2_xy_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 117);

    auto ta2_xy_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 118);

    auto ta2_xy_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 120);

    auto ta2_xy_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 121);

    auto ta2_xy_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 122);

    auto ta2_xy_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 123);

    auto ta2_xy_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 125);

    auto ta2_xy_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 126);

    auto ta2_xy_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 129);

    auto ta2_xy_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 135);

    auto ta2_xy_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 136);

    auto ta2_xy_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 137);

    auto ta2_xy_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 138);

    auto ta2_xy_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 139);

    auto ta2_xy_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 140);

    auto ta2_xy_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 141);

    auto ta2_xy_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 142);

    auto ta2_xy_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 143);

    auto ta2_xy_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 144);

    auto ta2_xy_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 145);

    auto ta2_xy_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 146);

    auto ta2_xy_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 147);

    auto ta2_xy_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 148);

    auto ta2_xy_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 149);

    auto ta2_xy_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 151);

    auto ta2_xy_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 153);

    auto ta2_xy_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 156);

    auto ta2_xy_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 160);

    auto ta2_xy_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 161);

    auto ta2_xy_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 162);

    auto ta2_xy_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 163);

    auto ta2_xy_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 164);

    auto ta2_xy_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 165);

    auto ta2_xy_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 166);

    auto ta2_xy_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 167);

    auto ta2_xy_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 168);

    auto ta2_xy_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 169);

    auto ta2_xy_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 170);

    auto ta2_xy_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 171);

    auto ta2_xy_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 172);

    auto ta2_xy_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 173);

    auto ta2_xy_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 174);

    auto ta2_xy_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 175);

    auto ta2_xy_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 176);

    auto ta2_xy_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 177);

    auto ta2_xy_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 178);

    auto ta2_xy_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 179);

    auto ta2_xz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 180);

    auto ta2_xz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 181);

    auto ta2_xz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 182);

    auto ta2_xz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 183);

    auto ta2_xz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 184);

    auto ta2_xz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 185);

    auto ta2_xz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 186);

    auto ta2_xz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 187);

    auto ta2_xz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 188);

    auto ta2_xz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 189);

    auto ta2_xz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 190);

    auto ta2_xz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 191);

    auto ta2_xz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 192);

    auto ta2_xz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 193);

    auto ta2_xz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 194);

    auto ta2_xz_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 195);

    auto ta2_xz_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 196);

    auto ta2_xz_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 197);

    auto ta2_xz_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 198);

    auto ta2_xz_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 200);

    auto ta2_xz_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 201);

    auto ta2_xz_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 204);

    auto ta2_xz_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 210);

    auto ta2_xz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 212);

    auto ta2_xz_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 214);

    auto ta2_xz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 215);

    auto ta2_xz_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 217);

    auto ta2_xz_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 218);

    auto ta2_xz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 219);

    auto ta2_xz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 221);

    auto ta2_xz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 222);

    auto ta2_xz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 223);

    auto ta2_xz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 224);

    auto ta2_xz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 225);

    auto ta2_xz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 226);

    auto ta2_xz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 227);

    auto ta2_xz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 228);

    auto ta2_xz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 229);

    auto ta2_xz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 230);

    auto ta2_xz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 231);

    auto ta2_xz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 232);

    auto ta2_xz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 233);

    auto ta2_xz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 234);

    auto ta2_xz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 235);

    auto ta2_xz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 236);

    auto ta2_xz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 237);

    auto ta2_xz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 238);

    auto ta2_xz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 239);

    auto ta2_xz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 242);

    auto ta2_xz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 245);

    auto ta2_xz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 249);

    auto ta2_xz_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 250);

    auto ta2_xz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 251);

    auto ta2_xz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 252);

    auto ta2_xz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 253);

    auto ta2_xz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 254);

    auto ta2_xz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 255);

    auto ta2_xz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 256);

    auto ta2_xz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 257);

    auto ta2_xz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 258);

    auto ta2_xz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 259);

    auto ta2_xz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 260);

    auto ta2_xz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 261);

    auto ta2_xz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 262);

    auto ta2_xz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 263);

    auto ta2_xz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 264);

    auto ta2_xz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 265);

    auto ta2_xz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 266);

    auto ta2_xz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 267);

    auto ta2_xz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 268);

    auto ta2_xz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 269);

    auto ta2_yy_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 270);

    auto ta2_yy_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 271);

    auto ta2_yy_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 272);

    auto ta2_yy_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 273);

    auto ta2_yy_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 274);

    auto ta2_yy_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 275);

    auto ta2_yy_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 276);

    auto ta2_yy_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 277);

    auto ta2_yy_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 278);

    auto ta2_yy_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 279);

    auto ta2_yy_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 280);

    auto ta2_yy_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 281);

    auto ta2_yy_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 282);

    auto ta2_yy_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 283);

    auto ta2_yy_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 284);

    auto ta2_yy_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 285);

    auto ta2_yy_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 286);

    auto ta2_yy_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 288);

    auto ta2_yy_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 291);

    auto ta2_yy_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 295);

    auto ta2_yy_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 296);

    auto ta2_yy_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 297);

    auto ta2_yy_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 298);

    auto ta2_yy_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 302);

    auto ta2_yy_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 305);

    auto ta2_yy_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 309);

    auto ta2_yy_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 311);

    auto ta2_yy_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 312);

    auto ta2_yy_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 313);

    auto ta2_yy_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 314);

    auto ta2_yy_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 315);

    auto ta2_yy_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 316);

    auto ta2_yy_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 317);

    auto ta2_yy_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 318);

    auto ta2_yy_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 319);

    auto ta2_yy_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 320);

    auto ta2_yy_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 321);

    auto ta2_yy_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 322);

    auto ta2_yy_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 323);

    auto ta2_yy_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 324);

    auto ta2_yy_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 325);

    auto ta2_yy_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 326);

    auto ta2_yy_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 327);

    auto ta2_yy_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 328);

    auto ta2_yy_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 329);

    auto ta2_yy_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 331);

    auto ta2_yy_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 333);

    auto ta2_yy_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 334);

    auto ta2_yy_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 336);

    auto ta2_yy_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 337);

    auto ta2_yy_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 338);

    auto ta2_yy_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 340);

    auto ta2_yy_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 341);

    auto ta2_yy_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 342);

    auto ta2_yy_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 343);

    auto ta2_yy_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 344);

    auto ta2_yy_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 345);

    auto ta2_yy_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 346);

    auto ta2_yy_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 347);

    auto ta2_yy_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 348);

    auto ta2_yy_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 349);

    auto ta2_yy_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 350);

    auto ta2_yy_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 351);

    auto ta2_yy_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 352);

    auto ta2_yy_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 353);

    auto ta2_yy_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 354);

    auto ta2_yy_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 355);

    auto ta2_yy_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 356);

    auto ta2_yy_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 357);

    auto ta2_yy_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 358);

    auto ta2_yy_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 359);

    auto ta2_yz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 360);

    auto ta2_yz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 361);

    auto ta2_yz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 362);

    auto ta2_yz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 363);

    auto ta2_yz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 364);

    auto ta2_yz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 365);

    auto ta2_yz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 366);

    auto ta2_yz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 367);

    auto ta2_yz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 368);

    auto ta2_yz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 369);

    auto ta2_yz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 370);

    auto ta2_yz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 371);

    auto ta2_yz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 372);

    auto ta2_yz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 373);

    auto ta2_yz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 374);

    auto ta2_yz_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 376);

    auto ta2_yz_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 378);

    auto ta2_yz_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 381);

    auto ta2_yz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 385);

    auto ta2_yz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 386);

    auto ta2_yz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 387);

    auto ta2_yz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 388);

    auto ta2_yz_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 390);

    auto ta2_yz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 392);

    auto ta2_yz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 395);

    auto ta2_yz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 399);

    auto ta2_yz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 401);

    auto ta2_yz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 402);

    auto ta2_yz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 403);

    auto ta2_yz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 404);

    auto ta2_yz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 405);

    auto ta2_yz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 406);

    auto ta2_yz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 407);

    auto ta2_yz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 408);

    auto ta2_yz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 409);

    auto ta2_yz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 410);

    auto ta2_yz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 411);

    auto ta2_yz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 412);

    auto ta2_yz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 413);

    auto ta2_yz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 414);

    auto ta2_yz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 415);

    auto ta2_yz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 416);

    auto ta2_yz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 417);

    auto ta2_yz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 418);

    auto ta2_yz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 419);

    auto ta2_yz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 422);

    auto ta2_yz_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 424);

    auto ta2_yz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 425);

    auto ta2_yz_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 427);

    auto ta2_yz_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 428);

    auto ta2_yz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 429);

    auto ta2_yz_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 430);

    auto ta2_yz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 431);

    auto ta2_yz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 432);

    auto ta2_yz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 433);

    auto ta2_yz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 434);

    auto ta2_yz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 435);

    auto ta2_yz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 436);

    auto ta2_yz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 437);

    auto ta2_yz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 438);

    auto ta2_yz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 439);

    auto ta2_yz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 440);

    auto ta2_yz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 441);

    auto ta2_yz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 442);

    auto ta2_yz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 443);

    auto ta2_yz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 444);

    auto ta2_yz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 445);

    auto ta2_yz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 446);

    auto ta2_yz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 447);

    auto ta2_yz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 448);

    auto ta2_yz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 449);

    auto ta2_zz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 450);

    auto ta2_zz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 451);

    auto ta2_zz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 452);

    auto ta2_zz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 453);

    auto ta2_zz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 454);

    auto ta2_zz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 455);

    auto ta2_zz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 456);

    auto ta2_zz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 457);

    auto ta2_zz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 458);

    auto ta2_zz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 459);

    auto ta2_zz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 460);

    auto ta2_zz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 461);

    auto ta2_zz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 462);

    auto ta2_zz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 463);

    auto ta2_zz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 464);

    auto ta2_zz_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 466);

    auto ta2_zz_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 468);

    auto ta2_zz_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 471);

    auto ta2_zz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 475);

    auto ta2_zz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 476);

    auto ta2_zz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 477);

    auto ta2_zz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 478);

    auto ta2_zz_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 480);

    auto ta2_zz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 482);

    auto ta2_zz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 485);

    auto ta2_zz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 489);

    auto ta2_zz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 491);

    auto ta2_zz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 492);

    auto ta2_zz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 493);

    auto ta2_zz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 494);

    auto ta2_zz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 495);

    auto ta2_zz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 496);

    auto ta2_zz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 497);

    auto ta2_zz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 498);

    auto ta2_zz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 499);

    auto ta2_zz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 500);

    auto ta2_zz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 501);

    auto ta2_zz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 502);

    auto ta2_zz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 503);

    auto ta2_zz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 504);

    auto ta2_zz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 505);

    auto ta2_zz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 506);

    auto ta2_zz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 507);

    auto ta2_zz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 508);

    auto ta2_zz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 509);

    auto ta2_zz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 512);

    auto ta2_zz_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 514);

    auto ta2_zz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 515);

    auto ta2_zz_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 517);

    auto ta2_zz_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 518);

    auto ta2_zz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 519);

    auto ta2_zz_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 520);

    auto ta2_zz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 521);

    auto ta2_zz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 522);

    auto ta2_zz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 523);

    auto ta2_zz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 524);

    auto ta2_zz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 525);

    auto ta2_zz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 526);

    auto ta2_zz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 527);

    auto ta2_zz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 528);

    auto ta2_zz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 529);

    auto ta2_zz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 530);

    auto ta2_zz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 531);

    auto ta2_zz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 532);

    auto ta2_zz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 533);

    auto ta2_zz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 534);

    auto ta2_zz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 535);

    auto ta2_zz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 536);

    auto ta2_zz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 537);

    auto ta2_zz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 538);

    auto ta2_zz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 539);

    // Set up components of auxiliary buffer : DG

    auto ta2_xx_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg);

    auto ta2_xx_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 1);

    auto ta2_xx_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 2);

    auto ta2_xx_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 3);

    auto ta2_xx_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 4);

    auto ta2_xx_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 5);

    auto ta2_xx_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 6);

    auto ta2_xx_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 7);

    auto ta2_xx_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 8);

    auto ta2_xx_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 9);

    auto ta2_xx_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 10);

    auto ta2_xx_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 11);

    auto ta2_xx_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 12);

    auto ta2_xx_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 13);

    auto ta2_xx_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 14);

    auto ta2_xx_xy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 15);

    auto ta2_xx_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 16);

    auto ta2_xx_xy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 17);

    auto ta2_xx_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 18);

    auto ta2_xx_xy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 20);

    auto ta2_xx_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 21);

    auto ta2_xx_xy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 24);

    auto ta2_xx_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 25);

    auto ta2_xx_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 30);

    auto ta2_xx_xz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 31);

    auto ta2_xx_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 32);

    auto ta2_xx_xz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 33);

    auto ta2_xx_xz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 34);

    auto ta2_xx_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 35);

    auto ta2_xx_xz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 36);

    auto ta2_xx_xz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 37);

    auto ta2_xx_xz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 38);

    auto ta2_xx_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 39);

    auto ta2_xx_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 44);

    auto ta2_xx_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 45);

    auto ta2_xx_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 46);

    auto ta2_xx_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 47);

    auto ta2_xx_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 48);

    auto ta2_xx_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 49);

    auto ta2_xx_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 50);

    auto ta2_xx_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 51);

    auto ta2_xx_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 52);

    auto ta2_xx_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 53);

    auto ta2_xx_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 54);

    auto ta2_xx_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 55);

    auto ta2_xx_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 56);

    auto ta2_xx_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 57);

    auto ta2_xx_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 58);

    auto ta2_xx_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 59);

    auto ta2_xx_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 62);

    auto ta2_xx_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 65);

    auto ta2_xx_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 69);

    auto ta2_xx_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 71);

    auto ta2_xx_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 72);

    auto ta2_xx_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 73);

    auto ta2_xx_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 74);

    auto ta2_xx_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 75);

    auto ta2_xx_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 76);

    auto ta2_xx_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 77);

    auto ta2_xx_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 78);

    auto ta2_xx_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 79);

    auto ta2_xx_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 80);

    auto ta2_xx_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 81);

    auto ta2_xx_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 82);

    auto ta2_xx_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 83);

    auto ta2_xx_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 84);

    auto ta2_xx_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 85);

    auto ta2_xx_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 86);

    auto ta2_xx_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 87);

    auto ta2_xx_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 88);

    auto ta2_xx_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 89);

    auto ta2_xy_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 90);

    auto ta2_xy_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 91);

    auto ta2_xy_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 92);

    auto ta2_xy_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 93);

    auto ta2_xy_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 94);

    auto ta2_xy_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 95);

    auto ta2_xy_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 96);

    auto ta2_xy_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 97);

    auto ta2_xy_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 98);

    auto ta2_xy_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 99);

    auto ta2_xy_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 100);

    auto ta2_xy_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 101);

    auto ta2_xy_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 102);

    auto ta2_xy_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 103);

    auto ta2_xy_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 104);

    auto ta2_xy_xy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 105);

    auto ta2_xy_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 106);

    auto ta2_xy_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 108);

    auto ta2_xy_xy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 109);

    auto ta2_xy_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 111);

    auto ta2_xy_xy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 112);

    auto ta2_xy_xy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 113);

    auto ta2_xy_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 115);

    auto ta2_xy_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 116);

    auto ta2_xy_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 117);

    auto ta2_xy_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 118);

    auto ta2_xy_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 120);

    auto ta2_xy_xz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 121);

    auto ta2_xy_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 122);

    auto ta2_xy_xz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 123);

    auto ta2_xy_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 125);

    auto ta2_xy_xz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 126);

    auto ta2_xy_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 129);

    auto ta2_xy_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 135);

    auto ta2_xy_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 136);

    auto ta2_xy_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 137);

    auto ta2_xy_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 138);

    auto ta2_xy_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 139);

    auto ta2_xy_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 140);

    auto ta2_xy_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 141);

    auto ta2_xy_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 142);

    auto ta2_xy_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 143);

    auto ta2_xy_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 144);

    auto ta2_xy_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 145);

    auto ta2_xy_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 146);

    auto ta2_xy_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 147);

    auto ta2_xy_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 148);

    auto ta2_xy_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 149);

    auto ta2_xy_yz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 151);

    auto ta2_xy_yz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 153);

    auto ta2_xy_yz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 156);

    auto ta2_xy_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 160);

    auto ta2_xy_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 161);

    auto ta2_xy_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 162);

    auto ta2_xy_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 163);

    auto ta2_xy_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 164);

    auto ta2_xy_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 165);

    auto ta2_xy_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 166);

    auto ta2_xy_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 167);

    auto ta2_xy_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 168);

    auto ta2_xy_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 169);

    auto ta2_xy_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 170);

    auto ta2_xy_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 171);

    auto ta2_xy_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 172);

    auto ta2_xy_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 173);

    auto ta2_xy_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 174);

    auto ta2_xy_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 175);

    auto ta2_xy_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 176);

    auto ta2_xy_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 177);

    auto ta2_xy_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 178);

    auto ta2_xy_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 179);

    auto ta2_xz_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 180);

    auto ta2_xz_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 181);

    auto ta2_xz_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 182);

    auto ta2_xz_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 183);

    auto ta2_xz_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 184);

    auto ta2_xz_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 185);

    auto ta2_xz_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 186);

    auto ta2_xz_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 187);

    auto ta2_xz_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 188);

    auto ta2_xz_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 189);

    auto ta2_xz_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 190);

    auto ta2_xz_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 191);

    auto ta2_xz_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 192);

    auto ta2_xz_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 193);

    auto ta2_xz_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 194);

    auto ta2_xz_xy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 195);

    auto ta2_xz_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 196);

    auto ta2_xz_xy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 197);

    auto ta2_xz_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 198);

    auto ta2_xz_xy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 200);

    auto ta2_xz_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 201);

    auto ta2_xz_xy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 204);

    auto ta2_xz_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 210);

    auto ta2_xz_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 212);

    auto ta2_xz_xz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 214);

    auto ta2_xz_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 215);

    auto ta2_xz_xz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 217);

    auto ta2_xz_xz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 218);

    auto ta2_xz_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 219);

    auto ta2_xz_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 221);

    auto ta2_xz_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 222);

    auto ta2_xz_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 223);

    auto ta2_xz_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 224);

    auto ta2_xz_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 225);

    auto ta2_xz_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 226);

    auto ta2_xz_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 227);

    auto ta2_xz_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 228);

    auto ta2_xz_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 229);

    auto ta2_xz_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 230);

    auto ta2_xz_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 231);

    auto ta2_xz_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 232);

    auto ta2_xz_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 233);

    auto ta2_xz_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 234);

    auto ta2_xz_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 235);

    auto ta2_xz_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 236);

    auto ta2_xz_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 237);

    auto ta2_xz_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 238);

    auto ta2_xz_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 239);

    auto ta2_xz_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 242);

    auto ta2_xz_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 245);

    auto ta2_xz_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 249);

    auto ta2_xz_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 250);

    auto ta2_xz_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 251);

    auto ta2_xz_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 252);

    auto ta2_xz_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 253);

    auto ta2_xz_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 254);

    auto ta2_xz_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 255);

    auto ta2_xz_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 256);

    auto ta2_xz_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 257);

    auto ta2_xz_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 258);

    auto ta2_xz_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 259);

    auto ta2_xz_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 260);

    auto ta2_xz_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 261);

    auto ta2_xz_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 262);

    auto ta2_xz_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 263);

    auto ta2_xz_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 264);

    auto ta2_xz_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 265);

    auto ta2_xz_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 266);

    auto ta2_xz_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 267);

    auto ta2_xz_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 268);

    auto ta2_xz_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 269);

    auto ta2_yy_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 270);

    auto ta2_yy_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 271);

    auto ta2_yy_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 272);

    auto ta2_yy_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 273);

    auto ta2_yy_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 274);

    auto ta2_yy_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 275);

    auto ta2_yy_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 276);

    auto ta2_yy_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 277);

    auto ta2_yy_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 278);

    auto ta2_yy_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 279);

    auto ta2_yy_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 280);

    auto ta2_yy_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 281);

    auto ta2_yy_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 282);

    auto ta2_yy_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 283);

    auto ta2_yy_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 284);

    auto ta2_yy_xy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 285);

    auto ta2_yy_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 286);

    auto ta2_yy_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 288);

    auto ta2_yy_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 291);

    auto ta2_yy_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 295);

    auto ta2_yy_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 296);

    auto ta2_yy_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 297);

    auto ta2_yy_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 298);

    auto ta2_yy_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 302);

    auto ta2_yy_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 305);

    auto ta2_yy_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 309);

    auto ta2_yy_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 311);

    auto ta2_yy_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 312);

    auto ta2_yy_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 313);

    auto ta2_yy_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 314);

    auto ta2_yy_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 315);

    auto ta2_yy_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 316);

    auto ta2_yy_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 317);

    auto ta2_yy_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 318);

    auto ta2_yy_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 319);

    auto ta2_yy_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 320);

    auto ta2_yy_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 321);

    auto ta2_yy_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 322);

    auto ta2_yy_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 323);

    auto ta2_yy_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 324);

    auto ta2_yy_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 325);

    auto ta2_yy_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 326);

    auto ta2_yy_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 327);

    auto ta2_yy_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 328);

    auto ta2_yy_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 329);

    auto ta2_yy_yz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 331);

    auto ta2_yy_yz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 333);

    auto ta2_yy_yz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 334);

    auto ta2_yy_yz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 336);

    auto ta2_yy_yz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 337);

    auto ta2_yy_yz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 338);

    auto ta2_yy_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 340);

    auto ta2_yy_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 341);

    auto ta2_yy_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 342);

    auto ta2_yy_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 343);

    auto ta2_yy_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 344);

    auto ta2_yy_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 345);

    auto ta2_yy_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 346);

    auto ta2_yy_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 347);

    auto ta2_yy_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 348);

    auto ta2_yy_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 349);

    auto ta2_yy_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 350);

    auto ta2_yy_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 351);

    auto ta2_yy_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 352);

    auto ta2_yy_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 353);

    auto ta2_yy_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 354);

    auto ta2_yy_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 355);

    auto ta2_yy_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 356);

    auto ta2_yy_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 357);

    auto ta2_yy_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 358);

    auto ta2_yy_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 359);

    auto ta2_yz_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 360);

    auto ta2_yz_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 361);

    auto ta2_yz_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 362);

    auto ta2_yz_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 363);

    auto ta2_yz_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 364);

    auto ta2_yz_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 365);

    auto ta2_yz_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 366);

    auto ta2_yz_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 367);

    auto ta2_yz_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 368);

    auto ta2_yz_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 369);

    auto ta2_yz_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 370);

    auto ta2_yz_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 371);

    auto ta2_yz_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 372);

    auto ta2_yz_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 373);

    auto ta2_yz_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 374);

    auto ta2_yz_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 376);

    auto ta2_yz_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 378);

    auto ta2_yz_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 381);

    auto ta2_yz_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 385);

    auto ta2_yz_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 386);

    auto ta2_yz_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 387);

    auto ta2_yz_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 388);

    auto ta2_yz_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 390);

    auto ta2_yz_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 392);

    auto ta2_yz_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 395);

    auto ta2_yz_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 399);

    auto ta2_yz_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 401);

    auto ta2_yz_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 402);

    auto ta2_yz_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 403);

    auto ta2_yz_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 404);

    auto ta2_yz_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 405);

    auto ta2_yz_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 406);

    auto ta2_yz_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 407);

    auto ta2_yz_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 408);

    auto ta2_yz_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 409);

    auto ta2_yz_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 410);

    auto ta2_yz_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 411);

    auto ta2_yz_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 412);

    auto ta2_yz_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 413);

    auto ta2_yz_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 414);

    auto ta2_yz_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 415);

    auto ta2_yz_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 416);

    auto ta2_yz_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 417);

    auto ta2_yz_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 418);

    auto ta2_yz_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 419);

    auto ta2_yz_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 422);

    auto ta2_yz_yz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 424);

    auto ta2_yz_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 425);

    auto ta2_yz_yz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 427);

    auto ta2_yz_yz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 428);

    auto ta2_yz_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 429);

    auto ta2_yz_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 430);

    auto ta2_yz_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 431);

    auto ta2_yz_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 432);

    auto ta2_yz_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 433);

    auto ta2_yz_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 434);

    auto ta2_yz_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 435);

    auto ta2_yz_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 436);

    auto ta2_yz_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 437);

    auto ta2_yz_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 438);

    auto ta2_yz_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 439);

    auto ta2_yz_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 440);

    auto ta2_yz_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 441);

    auto ta2_yz_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 442);

    auto ta2_yz_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 443);

    auto ta2_yz_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 444);

    auto ta2_yz_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 445);

    auto ta2_yz_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 446);

    auto ta2_yz_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 447);

    auto ta2_yz_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 448);

    auto ta2_yz_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 449);

    auto ta2_zz_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 450);

    auto ta2_zz_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 451);

    auto ta2_zz_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 452);

    auto ta2_zz_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 453);

    auto ta2_zz_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 454);

    auto ta2_zz_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 455);

    auto ta2_zz_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 456);

    auto ta2_zz_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 457);

    auto ta2_zz_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 458);

    auto ta2_zz_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 459);

    auto ta2_zz_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 460);

    auto ta2_zz_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 461);

    auto ta2_zz_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 462);

    auto ta2_zz_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 463);

    auto ta2_zz_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 464);

    auto ta2_zz_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 466);

    auto ta2_zz_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 468);

    auto ta2_zz_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 471);

    auto ta2_zz_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 475);

    auto ta2_zz_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 476);

    auto ta2_zz_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 477);

    auto ta2_zz_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 478);

    auto ta2_zz_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 480);

    auto ta2_zz_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 482);

    auto ta2_zz_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 485);

    auto ta2_zz_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 489);

    auto ta2_zz_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 491);

    auto ta2_zz_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 492);

    auto ta2_zz_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 493);

    auto ta2_zz_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 494);

    auto ta2_zz_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 495);

    auto ta2_zz_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 496);

    auto ta2_zz_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 497);

    auto ta2_zz_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 498);

    auto ta2_zz_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 499);

    auto ta2_zz_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 500);

    auto ta2_zz_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 501);

    auto ta2_zz_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 502);

    auto ta2_zz_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 503);

    auto ta2_zz_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 504);

    auto ta2_zz_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 505);

    auto ta2_zz_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 506);

    auto ta2_zz_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 507);

    auto ta2_zz_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 508);

    auto ta2_zz_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 509);

    auto ta2_zz_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 512);

    auto ta2_zz_yz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 514);

    auto ta2_zz_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 515);

    auto ta2_zz_yz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 517);

    auto ta2_zz_yz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 518);

    auto ta2_zz_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 519);

    auto ta2_zz_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 520);

    auto ta2_zz_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 521);

    auto ta2_zz_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 522);

    auto ta2_zz_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 523);

    auto ta2_zz_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 524);

    auto ta2_zz_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 525);

    auto ta2_zz_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 526);

    auto ta2_zz_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 527);

    auto ta2_zz_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 528);

    auto ta2_zz_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 529);

    auto ta2_zz_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 530);

    auto ta2_zz_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 531);

    auto ta2_zz_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 532);

    auto ta2_zz_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 533);

    auto ta2_zz_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 534);

    auto ta2_zz_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 535);

    auto ta2_zz_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 536);

    auto ta2_zz_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 537);

    auto ta2_zz_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 538);

    auto ta2_zz_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 539);

    // Set up 0-15 components of targeted buffer : FG

    auto ta2_xx_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg);

    auto ta2_xx_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 1);

    auto ta2_xx_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 2);

    auto ta2_xx_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 3);

    auto ta2_xx_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 4);

    auto ta2_xx_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 5);

    auto ta2_xx_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 6);

    auto ta2_xx_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 7);

    auto ta2_xx_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 8);

    auto ta2_xx_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 9);

    auto ta2_xx_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 10);

    auto ta2_xx_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 11);

    auto ta2_xx_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 12);

    auto ta2_xx_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 13);

    auto ta2_xx_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 14);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xx_xxxx_1, ta1_x_xx_xxxy_1, ta1_x_xx_xxxz_1, ta1_x_xx_xxyy_1, ta1_x_xx_xxyz_1, ta1_x_xx_xxzz_1, ta1_x_xx_xyyy_1, ta1_x_xx_xyyz_1, ta1_x_xx_xyzz_1, ta1_x_xx_xzzz_1, ta1_x_xx_yyyy_1, ta1_x_xx_yyyz_1, ta1_x_xx_yyzz_1, ta1_x_xx_yzzz_1, ta1_x_xx_zzzz_1, ta2_xx_x_xxxx_0, ta2_xx_x_xxxx_1, ta2_xx_x_xxxy_0, ta2_xx_x_xxxy_1, ta2_xx_x_xxxz_0, ta2_xx_x_xxxz_1, ta2_xx_x_xxyy_0, ta2_xx_x_xxyy_1, ta2_xx_x_xxyz_0, ta2_xx_x_xxyz_1, ta2_xx_x_xxzz_0, ta2_xx_x_xxzz_1, ta2_xx_x_xyyy_0, ta2_xx_x_xyyy_1, ta2_xx_x_xyyz_0, ta2_xx_x_xyyz_1, ta2_xx_x_xyzz_0, ta2_xx_x_xyzz_1, ta2_xx_x_xzzz_0, ta2_xx_x_xzzz_1, ta2_xx_x_yyyy_0, ta2_xx_x_yyyy_1, ta2_xx_x_yyyz_0, ta2_xx_x_yyyz_1, ta2_xx_x_yyzz_0, ta2_xx_x_yyzz_1, ta2_xx_x_yzzz_0, ta2_xx_x_yzzz_1, ta2_xx_x_zzzz_0, ta2_xx_x_zzzz_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxxx_0, ta2_xx_xx_xxxx_1, ta2_xx_xx_xxxy_0, ta2_xx_xx_xxxy_1, ta2_xx_xx_xxxz_0, ta2_xx_xx_xxxz_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxyy_0, ta2_xx_xx_xxyy_1, ta2_xx_xx_xxyz_0, ta2_xx_xx_xxyz_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xxzz_0, ta2_xx_xx_xxzz_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyyy_0, ta2_xx_xx_xyyy_1, ta2_xx_xx_xyyz_0, ta2_xx_xx_xyyz_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xyzz_0, ta2_xx_xx_xyzz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_xzzz_0, ta2_xx_xx_xzzz_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyyy_0, ta2_xx_xx_yyyy_1, ta2_xx_xx_yyyz_0, ta2_xx_xx_yyyz_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yyzz_0, ta2_xx_xx_yyzz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_yzzz_0, ta2_xx_xx_yzzz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xx_zzzz_0, ta2_xx_xx_zzzz_1, ta2_xx_xxx_xxxx_0, ta2_xx_xxx_xxxy_0, ta2_xx_xxx_xxxz_0, ta2_xx_xxx_xxyy_0, ta2_xx_xxx_xxyz_0, ta2_xx_xxx_xxzz_0, ta2_xx_xxx_xyyy_0, ta2_xx_xxx_xyyz_0, ta2_xx_xxx_xyzz_0, ta2_xx_xxx_xzzz_0, ta2_xx_xxx_yyyy_0, ta2_xx_xxx_yyyz_0, ta2_xx_xxx_yyzz_0, ta2_xx_xxx_yzzz_0, ta2_xx_xxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxx_xxxx_0[i] = 2.0 * ta2_xx_x_xxxx_0[i] * fe_0 - 2.0 * ta2_xx_x_xxxx_1[i] * fe_0 + 4.0 * ta2_xx_xx_xxx_0[i] * fe_0 - 4.0 * ta2_xx_xx_xxx_1[i] * fe_0 + 2.0 * ta1_x_xx_xxxx_1[i] + ta2_xx_xx_xxxx_0[i] * pa_x[i] - ta2_xx_xx_xxxx_1[i] * pc_x[i];

        ta2_xx_xxx_xxxy_0[i] = 2.0 * ta2_xx_x_xxxy_0[i] * fe_0 - 2.0 * ta2_xx_x_xxxy_1[i] * fe_0 + 3.0 * ta2_xx_xx_xxy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxy_1[i] * fe_0 + 2.0 * ta1_x_xx_xxxy_1[i] + ta2_xx_xx_xxxy_0[i] * pa_x[i] - ta2_xx_xx_xxxy_1[i] * pc_x[i];

        ta2_xx_xxx_xxxz_0[i] = 2.0 * ta2_xx_x_xxxz_0[i] * fe_0 - 2.0 * ta2_xx_x_xxxz_1[i] * fe_0 + 3.0 * ta2_xx_xx_xxz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxz_1[i] * fe_0 + 2.0 * ta1_x_xx_xxxz_1[i] + ta2_xx_xx_xxxz_0[i] * pa_x[i] - ta2_xx_xx_xxxz_1[i] * pc_x[i];

        ta2_xx_xxx_xxyy_0[i] = 2.0 * ta2_xx_x_xxyy_0[i] * fe_0 - 2.0 * ta2_xx_x_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_xx_xyy_0[i] * fe_0 - 2.0 * ta2_xx_xx_xyy_1[i] * fe_0 + 2.0 * ta1_x_xx_xxyy_1[i] + ta2_xx_xx_xxyy_0[i] * pa_x[i] - ta2_xx_xx_xxyy_1[i] * pc_x[i];

        ta2_xx_xxx_xxyz_0[i] = 2.0 * ta2_xx_x_xxyz_0[i] * fe_0 - 2.0 * ta2_xx_x_xxyz_1[i] * fe_0 + 2.0 * ta2_xx_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xyz_1[i] * fe_0 + 2.0 * ta1_x_xx_xxyz_1[i] + ta2_xx_xx_xxyz_0[i] * pa_x[i] - ta2_xx_xx_xxyz_1[i] * pc_x[i];

        ta2_xx_xxx_xxzz_0[i] = 2.0 * ta2_xx_x_xxzz_0[i] * fe_0 - 2.0 * ta2_xx_x_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_xx_xzz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xxzz_1[i] + ta2_xx_xx_xxzz_0[i] * pa_x[i] - ta2_xx_xx_xxzz_1[i] * pc_x[i];

        ta2_xx_xxx_xyyy_0[i] = 2.0 * ta2_xx_x_xyyy_0[i] * fe_0 - 2.0 * ta2_xx_x_xyyy_1[i] * fe_0 + ta2_xx_xx_yyy_0[i] * fe_0 - ta2_xx_xx_yyy_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyy_1[i] + ta2_xx_xx_xyyy_0[i] * pa_x[i] - ta2_xx_xx_xyyy_1[i] * pc_x[i];

        ta2_xx_xxx_xyyz_0[i] = 2.0 * ta2_xx_x_xyyz_0[i] * fe_0 - 2.0 * ta2_xx_x_xyyz_1[i] * fe_0 + ta2_xx_xx_yyz_0[i] * fe_0 - ta2_xx_xx_yyz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyz_1[i] + ta2_xx_xx_xyyz_0[i] * pa_x[i] - ta2_xx_xx_xyyz_1[i] * pc_x[i];

        ta2_xx_xxx_xyzz_0[i] = 2.0 * ta2_xx_x_xyzz_0[i] * fe_0 - 2.0 * ta2_xx_x_xyzz_1[i] * fe_0 + ta2_xx_xx_yzz_0[i] * fe_0 - ta2_xx_xx_yzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyzz_1[i] + ta2_xx_xx_xyzz_0[i] * pa_x[i] - ta2_xx_xx_xyzz_1[i] * pc_x[i];

        ta2_xx_xxx_xzzz_0[i] = 2.0 * ta2_xx_x_xzzz_0[i] * fe_0 - 2.0 * ta2_xx_x_xzzz_1[i] * fe_0 + ta2_xx_xx_zzz_0[i] * fe_0 - ta2_xx_xx_zzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xzzz_1[i] + ta2_xx_xx_xzzz_0[i] * pa_x[i] - ta2_xx_xx_xzzz_1[i] * pc_x[i];

        ta2_xx_xxx_yyyy_0[i] = 2.0 * ta2_xx_x_yyyy_0[i] * fe_0 - 2.0 * ta2_xx_x_yyyy_1[i] * fe_0 + 2.0 * ta1_x_xx_yyyy_1[i] + ta2_xx_xx_yyyy_0[i] * pa_x[i] - ta2_xx_xx_yyyy_1[i] * pc_x[i];

        ta2_xx_xxx_yyyz_0[i] = 2.0 * ta2_xx_x_yyyz_0[i] * fe_0 - 2.0 * ta2_xx_x_yyyz_1[i] * fe_0 + 2.0 * ta1_x_xx_yyyz_1[i] + ta2_xx_xx_yyyz_0[i] * pa_x[i] - ta2_xx_xx_yyyz_1[i] * pc_x[i];

        ta2_xx_xxx_yyzz_0[i] = 2.0 * ta2_xx_x_yyzz_0[i] * fe_0 - 2.0 * ta2_xx_x_yyzz_1[i] * fe_0 + 2.0 * ta1_x_xx_yyzz_1[i] + ta2_xx_xx_yyzz_0[i] * pa_x[i] - ta2_xx_xx_yyzz_1[i] * pc_x[i];

        ta2_xx_xxx_yzzz_0[i] = 2.0 * ta2_xx_x_yzzz_0[i] * fe_0 - 2.0 * ta2_xx_x_yzzz_1[i] * fe_0 + 2.0 * ta1_x_xx_yzzz_1[i] + ta2_xx_xx_yzzz_0[i] * pa_x[i] - ta2_xx_xx_yzzz_1[i] * pc_x[i];

        ta2_xx_xxx_zzzz_0[i] = 2.0 * ta2_xx_x_zzzz_0[i] * fe_0 - 2.0 * ta2_xx_x_zzzz_1[i] * fe_0 + 2.0 * ta1_x_xx_zzzz_1[i] + ta2_xx_xx_zzzz_0[i] * pa_x[i] - ta2_xx_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto ta2_xx_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 15);

    auto ta2_xx_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 16);

    auto ta2_xx_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 17);

    auto ta2_xx_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 18);

    auto ta2_xx_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 19);

    auto ta2_xx_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 20);

    auto ta2_xx_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 21);

    auto ta2_xx_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 22);

    auto ta2_xx_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 23);

    auto ta2_xx_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 24);

    auto ta2_xx_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 25);

    auto ta2_xx_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 26);

    auto ta2_xx_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 27);

    auto ta2_xx_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 28);

    auto ta2_xx_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 29);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxxx_0, ta2_xx_xx_xxxx_1, ta2_xx_xx_xxxy_0, ta2_xx_xx_xxxy_1, ta2_xx_xx_xxxz_0, ta2_xx_xx_xxxz_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxyy_0, ta2_xx_xx_xxyy_1, ta2_xx_xx_xxyz_0, ta2_xx_xx_xxyz_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xxzz_0, ta2_xx_xx_xxzz_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyyy_0, ta2_xx_xx_xyyy_1, ta2_xx_xx_xyyz_0, ta2_xx_xx_xyyz_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xyzz_0, ta2_xx_xx_xyzz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_xzzz_0, ta2_xx_xx_xzzz_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyyy_0, ta2_xx_xx_yyyy_1, ta2_xx_xx_yyyz_0, ta2_xx_xx_yyyz_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yyzz_0, ta2_xx_xx_yyzz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_yzzz_0, ta2_xx_xx_yzzz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xx_zzzz_0, ta2_xx_xx_zzzz_1, ta2_xx_xxy_xxxx_0, ta2_xx_xxy_xxxy_0, ta2_xx_xxy_xxxz_0, ta2_xx_xxy_xxyy_0, ta2_xx_xxy_xxyz_0, ta2_xx_xxy_xxzz_0, ta2_xx_xxy_xyyy_0, ta2_xx_xxy_xyyz_0, ta2_xx_xxy_xyzz_0, ta2_xx_xxy_xzzz_0, ta2_xx_xxy_yyyy_0, ta2_xx_xxy_yyyz_0, ta2_xx_xxy_yyzz_0, ta2_xx_xxy_yzzz_0, ta2_xx_xxy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxy_xxxx_0[i] = ta2_xx_xx_xxxx_0[i] * pa_y[i] - ta2_xx_xx_xxxx_1[i] * pc_y[i];

        ta2_xx_xxy_xxxy_0[i] = ta2_xx_xx_xxx_0[i] * fe_0 - ta2_xx_xx_xxx_1[i] * fe_0 + ta2_xx_xx_xxxy_0[i] * pa_y[i] - ta2_xx_xx_xxxy_1[i] * pc_y[i];

        ta2_xx_xxy_xxxz_0[i] = ta2_xx_xx_xxxz_0[i] * pa_y[i] - ta2_xx_xx_xxxz_1[i] * pc_y[i];

        ta2_xx_xxy_xxyy_0[i] = 2.0 * ta2_xx_xx_xxy_0[i] * fe_0 - 2.0 * ta2_xx_xx_xxy_1[i] * fe_0 + ta2_xx_xx_xxyy_0[i] * pa_y[i] - ta2_xx_xx_xxyy_1[i] * pc_y[i];

        ta2_xx_xxy_xxyz_0[i] = ta2_xx_xx_xxz_0[i] * fe_0 - ta2_xx_xx_xxz_1[i] * fe_0 + ta2_xx_xx_xxyz_0[i] * pa_y[i] - ta2_xx_xx_xxyz_1[i] * pc_y[i];

        ta2_xx_xxy_xxzz_0[i] = ta2_xx_xx_xxzz_0[i] * pa_y[i] - ta2_xx_xx_xxzz_1[i] * pc_y[i];

        ta2_xx_xxy_xyyy_0[i] = 3.0 * ta2_xx_xx_xyy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xyy_1[i] * fe_0 + ta2_xx_xx_xyyy_0[i] * pa_y[i] - ta2_xx_xx_xyyy_1[i] * pc_y[i];

        ta2_xx_xxy_xyyz_0[i] = 2.0 * ta2_xx_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xyz_1[i] * fe_0 + ta2_xx_xx_xyyz_0[i] * pa_y[i] - ta2_xx_xx_xyyz_1[i] * pc_y[i];

        ta2_xx_xxy_xyzz_0[i] = ta2_xx_xx_xzz_0[i] * fe_0 - ta2_xx_xx_xzz_1[i] * fe_0 + ta2_xx_xx_xyzz_0[i] * pa_y[i] - ta2_xx_xx_xyzz_1[i] * pc_y[i];

        ta2_xx_xxy_xzzz_0[i] = ta2_xx_xx_xzzz_0[i] * pa_y[i] - ta2_xx_xx_xzzz_1[i] * pc_y[i];

        ta2_xx_xxy_yyyy_0[i] = 4.0 * ta2_xx_xx_yyy_0[i] * fe_0 - 4.0 * ta2_xx_xx_yyy_1[i] * fe_0 + ta2_xx_xx_yyyy_0[i] * pa_y[i] - ta2_xx_xx_yyyy_1[i] * pc_y[i];

        ta2_xx_xxy_yyyz_0[i] = 3.0 * ta2_xx_xx_yyz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yyz_1[i] * fe_0 + ta2_xx_xx_yyyz_0[i] * pa_y[i] - ta2_xx_xx_yyyz_1[i] * pc_y[i];

        ta2_xx_xxy_yyzz_0[i] = 2.0 * ta2_xx_xx_yzz_0[i] * fe_0 - 2.0 * ta2_xx_xx_yzz_1[i] * fe_0 + ta2_xx_xx_yyzz_0[i] * pa_y[i] - ta2_xx_xx_yyzz_1[i] * pc_y[i];

        ta2_xx_xxy_yzzz_0[i] = ta2_xx_xx_zzz_0[i] * fe_0 - ta2_xx_xx_zzz_1[i] * fe_0 + ta2_xx_xx_yzzz_0[i] * pa_y[i] - ta2_xx_xx_yzzz_1[i] * pc_y[i];

        ta2_xx_xxy_zzzz_0[i] = ta2_xx_xx_zzzz_0[i] * pa_y[i] - ta2_xx_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto ta2_xx_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 30);

    auto ta2_xx_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 31);

    auto ta2_xx_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 32);

    auto ta2_xx_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 33);

    auto ta2_xx_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 34);

    auto ta2_xx_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 35);

    auto ta2_xx_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 36);

    auto ta2_xx_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 37);

    auto ta2_xx_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 38);

    auto ta2_xx_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 39);

    auto ta2_xx_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 40);

    auto ta2_xx_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 41);

    auto ta2_xx_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 42);

    auto ta2_xx_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 43);

    auto ta2_xx_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 44);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxxx_0, ta2_xx_xx_xxxx_1, ta2_xx_xx_xxxy_0, ta2_xx_xx_xxxy_1, ta2_xx_xx_xxxz_0, ta2_xx_xx_xxxz_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxyy_0, ta2_xx_xx_xxyy_1, ta2_xx_xx_xxyz_0, ta2_xx_xx_xxyz_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xxzz_0, ta2_xx_xx_xxzz_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyyy_0, ta2_xx_xx_xyyy_1, ta2_xx_xx_xyyz_0, ta2_xx_xx_xyyz_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xyzz_0, ta2_xx_xx_xyzz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_xzzz_0, ta2_xx_xx_xzzz_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyyy_0, ta2_xx_xx_yyyy_1, ta2_xx_xx_yyyz_0, ta2_xx_xx_yyyz_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yyzz_0, ta2_xx_xx_yyzz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_yzzz_0, ta2_xx_xx_yzzz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xx_zzzz_0, ta2_xx_xx_zzzz_1, ta2_xx_xxz_xxxx_0, ta2_xx_xxz_xxxy_0, ta2_xx_xxz_xxxz_0, ta2_xx_xxz_xxyy_0, ta2_xx_xxz_xxyz_0, ta2_xx_xxz_xxzz_0, ta2_xx_xxz_xyyy_0, ta2_xx_xxz_xyyz_0, ta2_xx_xxz_xyzz_0, ta2_xx_xxz_xzzz_0, ta2_xx_xxz_yyyy_0, ta2_xx_xxz_yyyz_0, ta2_xx_xxz_yyzz_0, ta2_xx_xxz_yzzz_0, ta2_xx_xxz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxz_xxxx_0[i] = ta2_xx_xx_xxxx_0[i] * pa_z[i] - ta2_xx_xx_xxxx_1[i] * pc_z[i];

        ta2_xx_xxz_xxxy_0[i] = ta2_xx_xx_xxxy_0[i] * pa_z[i] - ta2_xx_xx_xxxy_1[i] * pc_z[i];

        ta2_xx_xxz_xxxz_0[i] = ta2_xx_xx_xxx_0[i] * fe_0 - ta2_xx_xx_xxx_1[i] * fe_0 + ta2_xx_xx_xxxz_0[i] * pa_z[i] - ta2_xx_xx_xxxz_1[i] * pc_z[i];

        ta2_xx_xxz_xxyy_0[i] = ta2_xx_xx_xxyy_0[i] * pa_z[i] - ta2_xx_xx_xxyy_1[i] * pc_z[i];

        ta2_xx_xxz_xxyz_0[i] = ta2_xx_xx_xxy_0[i] * fe_0 - ta2_xx_xx_xxy_1[i] * fe_0 + ta2_xx_xx_xxyz_0[i] * pa_z[i] - ta2_xx_xx_xxyz_1[i] * pc_z[i];

        ta2_xx_xxz_xxzz_0[i] = 2.0 * ta2_xx_xx_xxz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xxz_1[i] * fe_0 + ta2_xx_xx_xxzz_0[i] * pa_z[i] - ta2_xx_xx_xxzz_1[i] * pc_z[i];

        ta2_xx_xxz_xyyy_0[i] = ta2_xx_xx_xyyy_0[i] * pa_z[i] - ta2_xx_xx_xyyy_1[i] * pc_z[i];

        ta2_xx_xxz_xyyz_0[i] = ta2_xx_xx_xyy_0[i] * fe_0 - ta2_xx_xx_xyy_1[i] * fe_0 + ta2_xx_xx_xyyz_0[i] * pa_z[i] - ta2_xx_xx_xyyz_1[i] * pc_z[i];

        ta2_xx_xxz_xyzz_0[i] = 2.0 * ta2_xx_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xyz_1[i] * fe_0 + ta2_xx_xx_xyzz_0[i] * pa_z[i] - ta2_xx_xx_xyzz_1[i] * pc_z[i];

        ta2_xx_xxz_xzzz_0[i] = 3.0 * ta2_xx_xx_xzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xzz_1[i] * fe_0 + ta2_xx_xx_xzzz_0[i] * pa_z[i] - ta2_xx_xx_xzzz_1[i] * pc_z[i];

        ta2_xx_xxz_yyyy_0[i] = ta2_xx_xx_yyyy_0[i] * pa_z[i] - ta2_xx_xx_yyyy_1[i] * pc_z[i];

        ta2_xx_xxz_yyyz_0[i] = ta2_xx_xx_yyy_0[i] * fe_0 - ta2_xx_xx_yyy_1[i] * fe_0 + ta2_xx_xx_yyyz_0[i] * pa_z[i] - ta2_xx_xx_yyyz_1[i] * pc_z[i];

        ta2_xx_xxz_yyzz_0[i] = 2.0 * ta2_xx_xx_yyz_0[i] * fe_0 - 2.0 * ta2_xx_xx_yyz_1[i] * fe_0 + ta2_xx_xx_yyzz_0[i] * pa_z[i] - ta2_xx_xx_yyzz_1[i] * pc_z[i];

        ta2_xx_xxz_yzzz_0[i] = 3.0 * ta2_xx_xx_yzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yzz_1[i] * fe_0 + ta2_xx_xx_yzzz_0[i] * pa_z[i] - ta2_xx_xx_yzzz_1[i] * pc_z[i];

        ta2_xx_xxz_zzzz_0[i] = 4.0 * ta2_xx_xx_zzz_0[i] * fe_0 - 4.0 * ta2_xx_xx_zzz_1[i] * fe_0 + ta2_xx_xx_zzzz_0[i] * pa_z[i] - ta2_xx_xx_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto ta2_xx_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 45);

    auto ta2_xx_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 46);

    auto ta2_xx_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 47);

    auto ta2_xx_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 48);

    auto ta2_xx_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 49);

    auto ta2_xx_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 50);

    auto ta2_xx_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 51);

    auto ta2_xx_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 52);

    auto ta2_xx_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 53);

    auto ta2_xx_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 54);

    auto ta2_xx_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 55);

    auto ta2_xx_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 56);

    auto ta2_xx_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 57);

    auto ta2_xx_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 58);

    auto ta2_xx_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 59);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_yy_xxxy_1, ta1_x_yy_xxyy_1, ta1_x_yy_xxyz_1, ta1_x_yy_xyyy_1, ta1_x_yy_xyyz_1, ta1_x_yy_xyzz_1, ta1_x_yy_yyyy_1, ta1_x_yy_yyyz_1, ta1_x_yy_yyzz_1, ta1_x_yy_yzzz_1, ta1_x_yy_zzzz_1, ta2_xx_x_xxxx_0, ta2_xx_x_xxxx_1, ta2_xx_x_xxxz_0, ta2_xx_x_xxxz_1, ta2_xx_x_xxzz_0, ta2_xx_x_xxzz_1, ta2_xx_x_xzzz_0, ta2_xx_x_xzzz_1, ta2_xx_xy_xxxx_0, ta2_xx_xy_xxxx_1, ta2_xx_xy_xxxz_0, ta2_xx_xy_xxxz_1, ta2_xx_xy_xxzz_0, ta2_xx_xy_xxzz_1, ta2_xx_xy_xzzz_0, ta2_xx_xy_xzzz_1, ta2_xx_xyy_xxxx_0, ta2_xx_xyy_xxxy_0, ta2_xx_xyy_xxxz_0, ta2_xx_xyy_xxyy_0, ta2_xx_xyy_xxyz_0, ta2_xx_xyy_xxzz_0, ta2_xx_xyy_xyyy_0, ta2_xx_xyy_xyyz_0, ta2_xx_xyy_xyzz_0, ta2_xx_xyy_xzzz_0, ta2_xx_xyy_yyyy_0, ta2_xx_xyy_yyyz_0, ta2_xx_xyy_yyzz_0, ta2_xx_xyy_yzzz_0, ta2_xx_xyy_zzzz_0, ta2_xx_yy_xxxy_0, ta2_xx_yy_xxxy_1, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xxyy_0, ta2_xx_yy_xxyy_1, ta2_xx_yy_xxyz_0, ta2_xx_yy_xxyz_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyyy_0, ta2_xx_yy_xyyy_1, ta2_xx_yy_xyyz_0, ta2_xx_yy_xyyz_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_xyzz_0, ta2_xx_yy_xyzz_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyyy_0, ta2_xx_yy_yyyy_1, ta2_xx_yy_yyyz_0, ta2_xx_yy_yyyz_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yyzz_0, ta2_xx_yy_yyzz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yy_yzzz_0, ta2_xx_yy_yzzz_1, ta2_xx_yy_zzzz_0, ta2_xx_yy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyy_xxxx_0[i] = ta2_xx_x_xxxx_0[i] * fe_0 - ta2_xx_x_xxxx_1[i] * fe_0 + ta2_xx_xy_xxxx_0[i] * pa_y[i] - ta2_xx_xy_xxxx_1[i] * pc_y[i];

        ta2_xx_xyy_xxxy_0[i] = 3.0 * ta2_xx_yy_xxy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxy_1[i] * fe_0 + 2.0 * ta1_x_yy_xxxy_1[i] + ta2_xx_yy_xxxy_0[i] * pa_x[i] - ta2_xx_yy_xxxy_1[i] * pc_x[i];

        ta2_xx_xyy_xxxz_0[i] = ta2_xx_x_xxxz_0[i] * fe_0 - ta2_xx_x_xxxz_1[i] * fe_0 + ta2_xx_xy_xxxz_0[i] * pa_y[i] - ta2_xx_xy_xxxz_1[i] * pc_y[i];

        ta2_xx_xyy_xxyy_0[i] = 2.0 * ta2_xx_yy_xyy_0[i] * fe_0 - 2.0 * ta2_xx_yy_xyy_1[i] * fe_0 + 2.0 * ta1_x_yy_xxyy_1[i] + ta2_xx_yy_xxyy_0[i] * pa_x[i] - ta2_xx_yy_xxyy_1[i] * pc_x[i];

        ta2_xx_xyy_xxyz_0[i] = 2.0 * ta2_xx_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xx_yy_xyz_1[i] * fe_0 + 2.0 * ta1_x_yy_xxyz_1[i] + ta2_xx_yy_xxyz_0[i] * pa_x[i] - ta2_xx_yy_xxyz_1[i] * pc_x[i];

        ta2_xx_xyy_xxzz_0[i] = ta2_xx_x_xxzz_0[i] * fe_0 - ta2_xx_x_xxzz_1[i] * fe_0 + ta2_xx_xy_xxzz_0[i] * pa_y[i] - ta2_xx_xy_xxzz_1[i] * pc_y[i];

        ta2_xx_xyy_xyyy_0[i] = ta2_xx_yy_yyy_0[i] * fe_0 - ta2_xx_yy_yyy_1[i] * fe_0 + 2.0 * ta1_x_yy_xyyy_1[i] + ta2_xx_yy_xyyy_0[i] * pa_x[i] - ta2_xx_yy_xyyy_1[i] * pc_x[i];

        ta2_xx_xyy_xyyz_0[i] = ta2_xx_yy_yyz_0[i] * fe_0 - ta2_xx_yy_yyz_1[i] * fe_0 + 2.0 * ta1_x_yy_xyyz_1[i] + ta2_xx_yy_xyyz_0[i] * pa_x[i] - ta2_xx_yy_xyyz_1[i] * pc_x[i];

        ta2_xx_xyy_xyzz_0[i] = ta2_xx_yy_yzz_0[i] * fe_0 - ta2_xx_yy_yzz_1[i] * fe_0 + 2.0 * ta1_x_yy_xyzz_1[i] + ta2_xx_yy_xyzz_0[i] * pa_x[i] - ta2_xx_yy_xyzz_1[i] * pc_x[i];

        ta2_xx_xyy_xzzz_0[i] = ta2_xx_x_xzzz_0[i] * fe_0 - ta2_xx_x_xzzz_1[i] * fe_0 + ta2_xx_xy_xzzz_0[i] * pa_y[i] - ta2_xx_xy_xzzz_1[i] * pc_y[i];

        ta2_xx_xyy_yyyy_0[i] = 2.0 * ta1_x_yy_yyyy_1[i] + ta2_xx_yy_yyyy_0[i] * pa_x[i] - ta2_xx_yy_yyyy_1[i] * pc_x[i];

        ta2_xx_xyy_yyyz_0[i] = 2.0 * ta1_x_yy_yyyz_1[i] + ta2_xx_yy_yyyz_0[i] * pa_x[i] - ta2_xx_yy_yyyz_1[i] * pc_x[i];

        ta2_xx_xyy_yyzz_0[i] = 2.0 * ta1_x_yy_yyzz_1[i] + ta2_xx_yy_yyzz_0[i] * pa_x[i] - ta2_xx_yy_yyzz_1[i] * pc_x[i];

        ta2_xx_xyy_yzzz_0[i] = 2.0 * ta1_x_yy_yzzz_1[i] + ta2_xx_yy_yzzz_0[i] * pa_x[i] - ta2_xx_yy_yzzz_1[i] * pc_x[i];

        ta2_xx_xyy_zzzz_0[i] = 2.0 * ta1_x_yy_zzzz_1[i] + ta2_xx_yy_zzzz_0[i] * pa_x[i] - ta2_xx_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto ta2_xx_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 60);

    auto ta2_xx_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 61);

    auto ta2_xx_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 62);

    auto ta2_xx_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 63);

    auto ta2_xx_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 64);

    auto ta2_xx_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 65);

    auto ta2_xx_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 66);

    auto ta2_xx_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 67);

    auto ta2_xx_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 68);

    auto ta2_xx_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 69);

    auto ta2_xx_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 70);

    auto ta2_xx_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 71);

    auto ta2_xx_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 72);

    auto ta2_xx_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 73);

    auto ta2_xx_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 74);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_yz_yyyz_1, ta1_x_yz_yyzz_1, ta1_x_yz_yzzz_1, ta2_xx_xy_xxxy_0, ta2_xx_xy_xxxy_1, ta2_xx_xy_xxyy_0, ta2_xx_xy_xxyy_1, ta2_xx_xy_xyyy_0, ta2_xx_xy_xyyy_1, ta2_xx_xy_yyyy_0, ta2_xx_xy_yyyy_1, ta2_xx_xyz_xxxx_0, ta2_xx_xyz_xxxy_0, ta2_xx_xyz_xxxz_0, ta2_xx_xyz_xxyy_0, ta2_xx_xyz_xxyz_0, ta2_xx_xyz_xxzz_0, ta2_xx_xyz_xyyy_0, ta2_xx_xyz_xyyz_0, ta2_xx_xyz_xyzz_0, ta2_xx_xyz_xzzz_0, ta2_xx_xyz_yyyy_0, ta2_xx_xyz_yyyz_0, ta2_xx_xyz_yyzz_0, ta2_xx_xyz_yzzz_0, ta2_xx_xyz_zzzz_0, ta2_xx_xz_xxxx_0, ta2_xx_xz_xxxx_1, ta2_xx_xz_xxxz_0, ta2_xx_xz_xxxz_1, ta2_xx_xz_xxyz_0, ta2_xx_xz_xxyz_1, ta2_xx_xz_xxz_0, ta2_xx_xz_xxz_1, ta2_xx_xz_xxzz_0, ta2_xx_xz_xxzz_1, ta2_xx_xz_xyyz_0, ta2_xx_xz_xyyz_1, ta2_xx_xz_xyz_0, ta2_xx_xz_xyz_1, ta2_xx_xz_xyzz_0, ta2_xx_xz_xyzz_1, ta2_xx_xz_xzz_0, ta2_xx_xz_xzz_1, ta2_xx_xz_xzzz_0, ta2_xx_xz_xzzz_1, ta2_xx_xz_zzzz_0, ta2_xx_xz_zzzz_1, ta2_xx_yz_yyyz_0, ta2_xx_yz_yyyz_1, ta2_xx_yz_yyzz_0, ta2_xx_yz_yyzz_1, ta2_xx_yz_yzzz_0, ta2_xx_yz_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyz_xxxx_0[i] = ta2_xx_xz_xxxx_0[i] * pa_y[i] - ta2_xx_xz_xxxx_1[i] * pc_y[i];

        ta2_xx_xyz_xxxy_0[i] = ta2_xx_xy_xxxy_0[i] * pa_z[i] - ta2_xx_xy_xxxy_1[i] * pc_z[i];

        ta2_xx_xyz_xxxz_0[i] = ta2_xx_xz_xxxz_0[i] * pa_y[i] - ta2_xx_xz_xxxz_1[i] * pc_y[i];

        ta2_xx_xyz_xxyy_0[i] = ta2_xx_xy_xxyy_0[i] * pa_z[i] - ta2_xx_xy_xxyy_1[i] * pc_z[i];

        ta2_xx_xyz_xxyz_0[i] = ta2_xx_xz_xxz_0[i] * fe_0 - ta2_xx_xz_xxz_1[i] * fe_0 + ta2_xx_xz_xxyz_0[i] * pa_y[i] - ta2_xx_xz_xxyz_1[i] * pc_y[i];

        ta2_xx_xyz_xxzz_0[i] = ta2_xx_xz_xxzz_0[i] * pa_y[i] - ta2_xx_xz_xxzz_1[i] * pc_y[i];

        ta2_xx_xyz_xyyy_0[i] = ta2_xx_xy_xyyy_0[i] * pa_z[i] - ta2_xx_xy_xyyy_1[i] * pc_z[i];

        ta2_xx_xyz_xyyz_0[i] = 2.0 * ta2_xx_xz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xz_xyz_1[i] * fe_0 + ta2_xx_xz_xyyz_0[i] * pa_y[i] - ta2_xx_xz_xyyz_1[i] * pc_y[i];

        ta2_xx_xyz_xyzz_0[i] = ta2_xx_xz_xzz_0[i] * fe_0 - ta2_xx_xz_xzz_1[i] * fe_0 + ta2_xx_xz_xyzz_0[i] * pa_y[i] - ta2_xx_xz_xyzz_1[i] * pc_y[i];

        ta2_xx_xyz_xzzz_0[i] = ta2_xx_xz_xzzz_0[i] * pa_y[i] - ta2_xx_xz_xzzz_1[i] * pc_y[i];

        ta2_xx_xyz_yyyy_0[i] = ta2_xx_xy_yyyy_0[i] * pa_z[i] - ta2_xx_xy_yyyy_1[i] * pc_z[i];

        ta2_xx_xyz_yyyz_0[i] = 2.0 * ta1_x_yz_yyyz_1[i] + ta2_xx_yz_yyyz_0[i] * pa_x[i] - ta2_xx_yz_yyyz_1[i] * pc_x[i];

        ta2_xx_xyz_yyzz_0[i] = 2.0 * ta1_x_yz_yyzz_1[i] + ta2_xx_yz_yyzz_0[i] * pa_x[i] - ta2_xx_yz_yyzz_1[i] * pc_x[i];

        ta2_xx_xyz_yzzz_0[i] = 2.0 * ta1_x_yz_yzzz_1[i] + ta2_xx_yz_yzzz_0[i] * pa_x[i] - ta2_xx_yz_yzzz_1[i] * pc_x[i];

        ta2_xx_xyz_zzzz_0[i] = ta2_xx_xz_zzzz_0[i] * pa_y[i] - ta2_xx_xz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto ta2_xx_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 75);

    auto ta2_xx_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 76);

    auto ta2_xx_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 77);

    auto ta2_xx_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 78);

    auto ta2_xx_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 79);

    auto ta2_xx_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 80);

    auto ta2_xx_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 81);

    auto ta2_xx_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 82);

    auto ta2_xx_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 83);

    auto ta2_xx_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 84);

    auto ta2_xx_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 85);

    auto ta2_xx_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 86);

    auto ta2_xx_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 87);

    auto ta2_xx_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 88);

    auto ta2_xx_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 89);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_zz_xxxz_1, ta1_x_zz_xxyz_1, ta1_x_zz_xxzz_1, ta1_x_zz_xyyz_1, ta1_x_zz_xyzz_1, ta1_x_zz_xzzz_1, ta1_x_zz_yyyy_1, ta1_x_zz_yyyz_1, ta1_x_zz_yyzz_1, ta1_x_zz_yzzz_1, ta1_x_zz_zzzz_1, ta2_xx_x_xxxx_0, ta2_xx_x_xxxx_1, ta2_xx_x_xxxy_0, ta2_xx_x_xxxy_1, ta2_xx_x_xxyy_0, ta2_xx_x_xxyy_1, ta2_xx_x_xyyy_0, ta2_xx_x_xyyy_1, ta2_xx_xz_xxxx_0, ta2_xx_xz_xxxx_1, ta2_xx_xz_xxxy_0, ta2_xx_xz_xxxy_1, ta2_xx_xz_xxyy_0, ta2_xx_xz_xxyy_1, ta2_xx_xz_xyyy_0, ta2_xx_xz_xyyy_1, ta2_xx_xzz_xxxx_0, ta2_xx_xzz_xxxy_0, ta2_xx_xzz_xxxz_0, ta2_xx_xzz_xxyy_0, ta2_xx_xzz_xxyz_0, ta2_xx_xzz_xxzz_0, ta2_xx_xzz_xyyy_0, ta2_xx_xzz_xyyz_0, ta2_xx_xzz_xyzz_0, ta2_xx_xzz_xzzz_0, ta2_xx_xzz_yyyy_0, ta2_xx_xzz_yyyz_0, ta2_xx_xzz_yyzz_0, ta2_xx_xzz_yzzz_0, ta2_xx_xzz_zzzz_0, ta2_xx_zz_xxxz_0, ta2_xx_zz_xxxz_1, ta2_xx_zz_xxyz_0, ta2_xx_zz_xxyz_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xxzz_0, ta2_xx_zz_xxzz_1, ta2_xx_zz_xyyz_0, ta2_xx_zz_xyyz_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xyzz_0, ta2_xx_zz_xyzz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_xzzz_0, ta2_xx_zz_xzzz_1, ta2_xx_zz_yyyy_0, ta2_xx_zz_yyyy_1, ta2_xx_zz_yyyz_0, ta2_xx_zz_yyyz_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yyzz_0, ta2_xx_zz_yyzz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_yzzz_0, ta2_xx_zz_yzzz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, ta2_xx_zz_zzzz_0, ta2_xx_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzz_xxxx_0[i] = ta2_xx_x_xxxx_0[i] * fe_0 - ta2_xx_x_xxxx_1[i] * fe_0 + ta2_xx_xz_xxxx_0[i] * pa_z[i] - ta2_xx_xz_xxxx_1[i] * pc_z[i];

        ta2_xx_xzz_xxxy_0[i] = ta2_xx_x_xxxy_0[i] * fe_0 - ta2_xx_x_xxxy_1[i] * fe_0 + ta2_xx_xz_xxxy_0[i] * pa_z[i] - ta2_xx_xz_xxxy_1[i] * pc_z[i];

        ta2_xx_xzz_xxxz_0[i] = 3.0 * ta2_xx_zz_xxz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxxz_1[i] + ta2_xx_zz_xxxz_0[i] * pa_x[i] - ta2_xx_zz_xxxz_1[i] * pc_x[i];

        ta2_xx_xzz_xxyy_0[i] = ta2_xx_x_xxyy_0[i] * fe_0 - ta2_xx_x_xxyy_1[i] * fe_0 + ta2_xx_xz_xxyy_0[i] * pa_z[i] - ta2_xx_xz_xxyy_1[i] * pc_z[i];

        ta2_xx_xzz_xxyz_0[i] = 2.0 * ta2_xx_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xyz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxyz_1[i] + ta2_xx_zz_xxyz_0[i] * pa_x[i] - ta2_xx_zz_xxyz_1[i] * pc_x[i];

        ta2_xx_xzz_xxzz_0[i] = 2.0 * ta2_xx_zz_xzz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxzz_1[i] + ta2_xx_zz_xxzz_0[i] * pa_x[i] - ta2_xx_zz_xxzz_1[i] * pc_x[i];

        ta2_xx_xzz_xyyy_0[i] = ta2_xx_x_xyyy_0[i] * fe_0 - ta2_xx_x_xyyy_1[i] * fe_0 + ta2_xx_xz_xyyy_0[i] * pa_z[i] - ta2_xx_xz_xyyy_1[i] * pc_z[i];

        ta2_xx_xzz_xyyz_0[i] = ta2_xx_zz_yyz_0[i] * fe_0 - ta2_xx_zz_yyz_1[i] * fe_0 + 2.0 * ta1_x_zz_xyyz_1[i] + ta2_xx_zz_xyyz_0[i] * pa_x[i] - ta2_xx_zz_xyyz_1[i] * pc_x[i];

        ta2_xx_xzz_xyzz_0[i] = ta2_xx_zz_yzz_0[i] * fe_0 - ta2_xx_zz_yzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xyzz_1[i] + ta2_xx_zz_xyzz_0[i] * pa_x[i] - ta2_xx_zz_xyzz_1[i] * pc_x[i];

        ta2_xx_xzz_xzzz_0[i] = ta2_xx_zz_zzz_0[i] * fe_0 - ta2_xx_zz_zzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xzzz_1[i] + ta2_xx_zz_xzzz_0[i] * pa_x[i] - ta2_xx_zz_xzzz_1[i] * pc_x[i];

        ta2_xx_xzz_yyyy_0[i] = 2.0 * ta1_x_zz_yyyy_1[i] + ta2_xx_zz_yyyy_0[i] * pa_x[i] - ta2_xx_zz_yyyy_1[i] * pc_x[i];

        ta2_xx_xzz_yyyz_0[i] = 2.0 * ta1_x_zz_yyyz_1[i] + ta2_xx_zz_yyyz_0[i] * pa_x[i] - ta2_xx_zz_yyyz_1[i] * pc_x[i];

        ta2_xx_xzz_yyzz_0[i] = 2.0 * ta1_x_zz_yyzz_1[i] + ta2_xx_zz_yyzz_0[i] * pa_x[i] - ta2_xx_zz_yyzz_1[i] * pc_x[i];

        ta2_xx_xzz_yzzz_0[i] = 2.0 * ta1_x_zz_yzzz_1[i] + ta2_xx_zz_yzzz_0[i] * pa_x[i] - ta2_xx_zz_yzzz_1[i] * pc_x[i];

        ta2_xx_xzz_zzzz_0[i] = 2.0 * ta1_x_zz_zzzz_1[i] + ta2_xx_zz_zzzz_0[i] * pa_x[i] - ta2_xx_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto ta2_xx_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 90);

    auto ta2_xx_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 91);

    auto ta2_xx_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 92);

    auto ta2_xx_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 93);

    auto ta2_xx_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 94);

    auto ta2_xx_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 95);

    auto ta2_xx_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 96);

    auto ta2_xx_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 97);

    auto ta2_xx_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 98);

    auto ta2_xx_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 99);

    auto ta2_xx_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 100);

    auto ta2_xx_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 101);

    auto ta2_xx_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 102);

    auto ta2_xx_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 103);

    auto ta2_xx_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 104);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_y_xxxx_0, ta2_xx_y_xxxx_1, ta2_xx_y_xxxy_0, ta2_xx_y_xxxy_1, ta2_xx_y_xxxz_0, ta2_xx_y_xxxz_1, ta2_xx_y_xxyy_0, ta2_xx_y_xxyy_1, ta2_xx_y_xxyz_0, ta2_xx_y_xxyz_1, ta2_xx_y_xxzz_0, ta2_xx_y_xxzz_1, ta2_xx_y_xyyy_0, ta2_xx_y_xyyy_1, ta2_xx_y_xyyz_0, ta2_xx_y_xyyz_1, ta2_xx_y_xyzz_0, ta2_xx_y_xyzz_1, ta2_xx_y_xzzz_0, ta2_xx_y_xzzz_1, ta2_xx_y_yyyy_0, ta2_xx_y_yyyy_1, ta2_xx_y_yyyz_0, ta2_xx_y_yyyz_1, ta2_xx_y_yyzz_0, ta2_xx_y_yyzz_1, ta2_xx_y_yzzz_0, ta2_xx_y_yzzz_1, ta2_xx_y_zzzz_0, ta2_xx_y_zzzz_1, ta2_xx_yy_xxx_0, ta2_xx_yy_xxx_1, ta2_xx_yy_xxxx_0, ta2_xx_yy_xxxx_1, ta2_xx_yy_xxxy_0, ta2_xx_yy_xxxy_1, ta2_xx_yy_xxxz_0, ta2_xx_yy_xxxz_1, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xxyy_0, ta2_xx_yy_xxyy_1, ta2_xx_yy_xxyz_0, ta2_xx_yy_xxyz_1, ta2_xx_yy_xxz_0, ta2_xx_yy_xxz_1, ta2_xx_yy_xxzz_0, ta2_xx_yy_xxzz_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyyy_0, ta2_xx_yy_xyyy_1, ta2_xx_yy_xyyz_0, ta2_xx_yy_xyyz_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_xyzz_0, ta2_xx_yy_xyzz_1, ta2_xx_yy_xzz_0, ta2_xx_yy_xzz_1, ta2_xx_yy_xzzz_0, ta2_xx_yy_xzzz_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyyy_0, ta2_xx_yy_yyyy_1, ta2_xx_yy_yyyz_0, ta2_xx_yy_yyyz_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yyzz_0, ta2_xx_yy_yyzz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yy_yzzz_0, ta2_xx_yy_yzzz_1, ta2_xx_yy_zzz_0, ta2_xx_yy_zzz_1, ta2_xx_yy_zzzz_0, ta2_xx_yy_zzzz_1, ta2_xx_yyy_xxxx_0, ta2_xx_yyy_xxxy_0, ta2_xx_yyy_xxxz_0, ta2_xx_yyy_xxyy_0, ta2_xx_yyy_xxyz_0, ta2_xx_yyy_xxzz_0, ta2_xx_yyy_xyyy_0, ta2_xx_yyy_xyyz_0, ta2_xx_yyy_xyzz_0, ta2_xx_yyy_xzzz_0, ta2_xx_yyy_yyyy_0, ta2_xx_yyy_yyyz_0, ta2_xx_yyy_yyzz_0, ta2_xx_yyy_yzzz_0, ta2_xx_yyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyy_xxxx_0[i] = 2.0 * ta2_xx_y_xxxx_0[i] * fe_0 - 2.0 * ta2_xx_y_xxxx_1[i] * fe_0 + ta2_xx_yy_xxxx_0[i] * pa_y[i] - ta2_xx_yy_xxxx_1[i] * pc_y[i];

        ta2_xx_yyy_xxxy_0[i] = 2.0 * ta2_xx_y_xxxy_0[i] * fe_0 - 2.0 * ta2_xx_y_xxxy_1[i] * fe_0 + ta2_xx_yy_xxx_0[i] * fe_0 - ta2_xx_yy_xxx_1[i] * fe_0 + ta2_xx_yy_xxxy_0[i] * pa_y[i] - ta2_xx_yy_xxxy_1[i] * pc_y[i];

        ta2_xx_yyy_xxxz_0[i] = 2.0 * ta2_xx_y_xxxz_0[i] * fe_0 - 2.0 * ta2_xx_y_xxxz_1[i] * fe_0 + ta2_xx_yy_xxxz_0[i] * pa_y[i] - ta2_xx_yy_xxxz_1[i] * pc_y[i];

        ta2_xx_yyy_xxyy_0[i] = 2.0 * ta2_xx_y_xxyy_0[i] * fe_0 - 2.0 * ta2_xx_y_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_yy_xxy_0[i] * fe_0 - 2.0 * ta2_xx_yy_xxy_1[i] * fe_0 + ta2_xx_yy_xxyy_0[i] * pa_y[i] - ta2_xx_yy_xxyy_1[i] * pc_y[i];

        ta2_xx_yyy_xxyz_0[i] = 2.0 * ta2_xx_y_xxyz_0[i] * fe_0 - 2.0 * ta2_xx_y_xxyz_1[i] * fe_0 + ta2_xx_yy_xxz_0[i] * fe_0 - ta2_xx_yy_xxz_1[i] * fe_0 + ta2_xx_yy_xxyz_0[i] * pa_y[i] - ta2_xx_yy_xxyz_1[i] * pc_y[i];

        ta2_xx_yyy_xxzz_0[i] = 2.0 * ta2_xx_y_xxzz_0[i] * fe_0 - 2.0 * ta2_xx_y_xxzz_1[i] * fe_0 + ta2_xx_yy_xxzz_0[i] * pa_y[i] - ta2_xx_yy_xxzz_1[i] * pc_y[i];

        ta2_xx_yyy_xyyy_0[i] = 2.0 * ta2_xx_y_xyyy_0[i] * fe_0 - 2.0 * ta2_xx_y_xyyy_1[i] * fe_0 + 3.0 * ta2_xx_yy_xyy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xyy_1[i] * fe_0 + ta2_xx_yy_xyyy_0[i] * pa_y[i] - ta2_xx_yy_xyyy_1[i] * pc_y[i];

        ta2_xx_yyy_xyyz_0[i] = 2.0 * ta2_xx_y_xyyz_0[i] * fe_0 - 2.0 * ta2_xx_y_xyyz_1[i] * fe_0 + 2.0 * ta2_xx_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xx_yy_xyz_1[i] * fe_0 + ta2_xx_yy_xyyz_0[i] * pa_y[i] - ta2_xx_yy_xyyz_1[i] * pc_y[i];

        ta2_xx_yyy_xyzz_0[i] = 2.0 * ta2_xx_y_xyzz_0[i] * fe_0 - 2.0 * ta2_xx_y_xyzz_1[i] * fe_0 + ta2_xx_yy_xzz_0[i] * fe_0 - ta2_xx_yy_xzz_1[i] * fe_0 + ta2_xx_yy_xyzz_0[i] * pa_y[i] - ta2_xx_yy_xyzz_1[i] * pc_y[i];

        ta2_xx_yyy_xzzz_0[i] = 2.0 * ta2_xx_y_xzzz_0[i] * fe_0 - 2.0 * ta2_xx_y_xzzz_1[i] * fe_0 + ta2_xx_yy_xzzz_0[i] * pa_y[i] - ta2_xx_yy_xzzz_1[i] * pc_y[i];

        ta2_xx_yyy_yyyy_0[i] = 2.0 * ta2_xx_y_yyyy_0[i] * fe_0 - 2.0 * ta2_xx_y_yyyy_1[i] * fe_0 + 4.0 * ta2_xx_yy_yyy_0[i] * fe_0 - 4.0 * ta2_xx_yy_yyy_1[i] * fe_0 + ta2_xx_yy_yyyy_0[i] * pa_y[i] - ta2_xx_yy_yyyy_1[i] * pc_y[i];

        ta2_xx_yyy_yyyz_0[i] = 2.0 * ta2_xx_y_yyyz_0[i] * fe_0 - 2.0 * ta2_xx_y_yyyz_1[i] * fe_0 + 3.0 * ta2_xx_yy_yyz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yyz_1[i] * fe_0 + ta2_xx_yy_yyyz_0[i] * pa_y[i] - ta2_xx_yy_yyyz_1[i] * pc_y[i];

        ta2_xx_yyy_yyzz_0[i] = 2.0 * ta2_xx_y_yyzz_0[i] * fe_0 - 2.0 * ta2_xx_y_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_yy_yzz_0[i] * fe_0 - 2.0 * ta2_xx_yy_yzz_1[i] * fe_0 + ta2_xx_yy_yyzz_0[i] * pa_y[i] - ta2_xx_yy_yyzz_1[i] * pc_y[i];

        ta2_xx_yyy_yzzz_0[i] = 2.0 * ta2_xx_y_yzzz_0[i] * fe_0 - 2.0 * ta2_xx_y_yzzz_1[i] * fe_0 + ta2_xx_yy_zzz_0[i] * fe_0 - ta2_xx_yy_zzz_1[i] * fe_0 + ta2_xx_yy_yzzz_0[i] * pa_y[i] - ta2_xx_yy_yzzz_1[i] * pc_y[i];

        ta2_xx_yyy_zzzz_0[i] = 2.0 * ta2_xx_y_zzzz_0[i] * fe_0 - 2.0 * ta2_xx_y_zzzz_1[i] * fe_0 + ta2_xx_yy_zzzz_0[i] * pa_y[i] - ta2_xx_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto ta2_xx_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 105);

    auto ta2_xx_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 106);

    auto ta2_xx_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 107);

    auto ta2_xx_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 108);

    auto ta2_xx_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 109);

    auto ta2_xx_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 110);

    auto ta2_xx_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 111);

    auto ta2_xx_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 112);

    auto ta2_xx_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 113);

    auto ta2_xx_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 114);

    auto ta2_xx_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 115);

    auto ta2_xx_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 116);

    auto ta2_xx_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 117);

    auto ta2_xx_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 118);

    auto ta2_xx_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 119);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_yy_xxxx_0, ta2_xx_yy_xxxx_1, ta2_xx_yy_xxxy_0, ta2_xx_yy_xxxy_1, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xxyy_0, ta2_xx_yy_xxyy_1, ta2_xx_yy_xxyz_0, ta2_xx_yy_xxyz_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyyy_0, ta2_xx_yy_xyyy_1, ta2_xx_yy_xyyz_0, ta2_xx_yy_xyyz_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_xyzz_0, ta2_xx_yy_xyzz_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyyy_0, ta2_xx_yy_yyyy_1, ta2_xx_yy_yyyz_0, ta2_xx_yy_yyyz_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yyzz_0, ta2_xx_yy_yyzz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yy_yzzz_0, ta2_xx_yy_yzzz_1, ta2_xx_yyz_xxxx_0, ta2_xx_yyz_xxxy_0, ta2_xx_yyz_xxxz_0, ta2_xx_yyz_xxyy_0, ta2_xx_yyz_xxyz_0, ta2_xx_yyz_xxzz_0, ta2_xx_yyz_xyyy_0, ta2_xx_yyz_xyyz_0, ta2_xx_yyz_xyzz_0, ta2_xx_yyz_xzzz_0, ta2_xx_yyz_yyyy_0, ta2_xx_yyz_yyyz_0, ta2_xx_yyz_yyzz_0, ta2_xx_yyz_yzzz_0, ta2_xx_yyz_zzzz_0, ta2_xx_yz_xxxz_0, ta2_xx_yz_xxxz_1, ta2_xx_yz_xxzz_0, ta2_xx_yz_xxzz_1, ta2_xx_yz_xzzz_0, ta2_xx_yz_xzzz_1, ta2_xx_yz_zzzz_0, ta2_xx_yz_zzzz_1, ta2_xx_z_xxxz_0, ta2_xx_z_xxxz_1, ta2_xx_z_xxzz_0, ta2_xx_z_xxzz_1, ta2_xx_z_xzzz_0, ta2_xx_z_xzzz_1, ta2_xx_z_zzzz_0, ta2_xx_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyz_xxxx_0[i] = ta2_xx_yy_xxxx_0[i] * pa_z[i] - ta2_xx_yy_xxxx_1[i] * pc_z[i];

        ta2_xx_yyz_xxxy_0[i] = ta2_xx_yy_xxxy_0[i] * pa_z[i] - ta2_xx_yy_xxxy_1[i] * pc_z[i];

        ta2_xx_yyz_xxxz_0[i] = ta2_xx_z_xxxz_0[i] * fe_0 - ta2_xx_z_xxxz_1[i] * fe_0 + ta2_xx_yz_xxxz_0[i] * pa_y[i] - ta2_xx_yz_xxxz_1[i] * pc_y[i];

        ta2_xx_yyz_xxyy_0[i] = ta2_xx_yy_xxyy_0[i] * pa_z[i] - ta2_xx_yy_xxyy_1[i] * pc_z[i];

        ta2_xx_yyz_xxyz_0[i] = ta2_xx_yy_xxy_0[i] * fe_0 - ta2_xx_yy_xxy_1[i] * fe_0 + ta2_xx_yy_xxyz_0[i] * pa_z[i] - ta2_xx_yy_xxyz_1[i] * pc_z[i];

        ta2_xx_yyz_xxzz_0[i] = ta2_xx_z_xxzz_0[i] * fe_0 - ta2_xx_z_xxzz_1[i] * fe_0 + ta2_xx_yz_xxzz_0[i] * pa_y[i] - ta2_xx_yz_xxzz_1[i] * pc_y[i];

        ta2_xx_yyz_xyyy_0[i] = ta2_xx_yy_xyyy_0[i] * pa_z[i] - ta2_xx_yy_xyyy_1[i] * pc_z[i];

        ta2_xx_yyz_xyyz_0[i] = ta2_xx_yy_xyy_0[i] * fe_0 - ta2_xx_yy_xyy_1[i] * fe_0 + ta2_xx_yy_xyyz_0[i] * pa_z[i] - ta2_xx_yy_xyyz_1[i] * pc_z[i];

        ta2_xx_yyz_xyzz_0[i] = 2.0 * ta2_xx_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xx_yy_xyz_1[i] * fe_0 + ta2_xx_yy_xyzz_0[i] * pa_z[i] - ta2_xx_yy_xyzz_1[i] * pc_z[i];

        ta2_xx_yyz_xzzz_0[i] = ta2_xx_z_xzzz_0[i] * fe_0 - ta2_xx_z_xzzz_1[i] * fe_0 + ta2_xx_yz_xzzz_0[i] * pa_y[i] - ta2_xx_yz_xzzz_1[i] * pc_y[i];

        ta2_xx_yyz_yyyy_0[i] = ta2_xx_yy_yyyy_0[i] * pa_z[i] - ta2_xx_yy_yyyy_1[i] * pc_z[i];

        ta2_xx_yyz_yyyz_0[i] = ta2_xx_yy_yyy_0[i] * fe_0 - ta2_xx_yy_yyy_1[i] * fe_0 + ta2_xx_yy_yyyz_0[i] * pa_z[i] - ta2_xx_yy_yyyz_1[i] * pc_z[i];

        ta2_xx_yyz_yyzz_0[i] = 2.0 * ta2_xx_yy_yyz_0[i] * fe_0 - 2.0 * ta2_xx_yy_yyz_1[i] * fe_0 + ta2_xx_yy_yyzz_0[i] * pa_z[i] - ta2_xx_yy_yyzz_1[i] * pc_z[i];

        ta2_xx_yyz_yzzz_0[i] = 3.0 * ta2_xx_yy_yzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yzz_1[i] * fe_0 + ta2_xx_yy_yzzz_0[i] * pa_z[i] - ta2_xx_yy_yzzz_1[i] * pc_z[i];

        ta2_xx_yyz_zzzz_0[i] = ta2_xx_z_zzzz_0[i] * fe_0 - ta2_xx_z_zzzz_1[i] * fe_0 + ta2_xx_yz_zzzz_0[i] * pa_y[i] - ta2_xx_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto ta2_xx_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 120);

    auto ta2_xx_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 121);

    auto ta2_xx_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 122);

    auto ta2_xx_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 123);

    auto ta2_xx_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 124);

    auto ta2_xx_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 125);

    auto ta2_xx_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 126);

    auto ta2_xx_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 127);

    auto ta2_xx_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 128);

    auto ta2_xx_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 129);

    auto ta2_xx_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 130);

    auto ta2_xx_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 131);

    auto ta2_xx_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 132);

    auto ta2_xx_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 133);

    auto ta2_xx_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 134);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_yzz_xxxx_0, ta2_xx_yzz_xxxy_0, ta2_xx_yzz_xxxz_0, ta2_xx_yzz_xxyy_0, ta2_xx_yzz_xxyz_0, ta2_xx_yzz_xxzz_0, ta2_xx_yzz_xyyy_0, ta2_xx_yzz_xyyz_0, ta2_xx_yzz_xyzz_0, ta2_xx_yzz_xzzz_0, ta2_xx_yzz_yyyy_0, ta2_xx_yzz_yyyz_0, ta2_xx_yzz_yyzz_0, ta2_xx_yzz_yzzz_0, ta2_xx_yzz_zzzz_0, ta2_xx_zz_xxx_0, ta2_xx_zz_xxx_1, ta2_xx_zz_xxxx_0, ta2_xx_zz_xxxx_1, ta2_xx_zz_xxxy_0, ta2_xx_zz_xxxy_1, ta2_xx_zz_xxxz_0, ta2_xx_zz_xxxz_1, ta2_xx_zz_xxy_0, ta2_xx_zz_xxy_1, ta2_xx_zz_xxyy_0, ta2_xx_zz_xxyy_1, ta2_xx_zz_xxyz_0, ta2_xx_zz_xxyz_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xxzz_0, ta2_xx_zz_xxzz_1, ta2_xx_zz_xyy_0, ta2_xx_zz_xyy_1, ta2_xx_zz_xyyy_0, ta2_xx_zz_xyyy_1, ta2_xx_zz_xyyz_0, ta2_xx_zz_xyyz_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xyzz_0, ta2_xx_zz_xyzz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_xzzz_0, ta2_xx_zz_xzzz_1, ta2_xx_zz_yyy_0, ta2_xx_zz_yyy_1, ta2_xx_zz_yyyy_0, ta2_xx_zz_yyyy_1, ta2_xx_zz_yyyz_0, ta2_xx_zz_yyyz_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yyzz_0, ta2_xx_zz_yyzz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_yzzz_0, ta2_xx_zz_yzzz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, ta2_xx_zz_zzzz_0, ta2_xx_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzz_xxxx_0[i] = ta2_xx_zz_xxxx_0[i] * pa_y[i] - ta2_xx_zz_xxxx_1[i] * pc_y[i];

        ta2_xx_yzz_xxxy_0[i] = ta2_xx_zz_xxx_0[i] * fe_0 - ta2_xx_zz_xxx_1[i] * fe_0 + ta2_xx_zz_xxxy_0[i] * pa_y[i] - ta2_xx_zz_xxxy_1[i] * pc_y[i];

        ta2_xx_yzz_xxxz_0[i] = ta2_xx_zz_xxxz_0[i] * pa_y[i] - ta2_xx_zz_xxxz_1[i] * pc_y[i];

        ta2_xx_yzz_xxyy_0[i] = 2.0 * ta2_xx_zz_xxy_0[i] * fe_0 - 2.0 * ta2_xx_zz_xxy_1[i] * fe_0 + ta2_xx_zz_xxyy_0[i] * pa_y[i] - ta2_xx_zz_xxyy_1[i] * pc_y[i];

        ta2_xx_yzz_xxyz_0[i] = ta2_xx_zz_xxz_0[i] * fe_0 - ta2_xx_zz_xxz_1[i] * fe_0 + ta2_xx_zz_xxyz_0[i] * pa_y[i] - ta2_xx_zz_xxyz_1[i] * pc_y[i];

        ta2_xx_yzz_xxzz_0[i] = ta2_xx_zz_xxzz_0[i] * pa_y[i] - ta2_xx_zz_xxzz_1[i] * pc_y[i];

        ta2_xx_yzz_xyyy_0[i] = 3.0 * ta2_xx_zz_xyy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xyy_1[i] * fe_0 + ta2_xx_zz_xyyy_0[i] * pa_y[i] - ta2_xx_zz_xyyy_1[i] * pc_y[i];

        ta2_xx_yzz_xyyz_0[i] = 2.0 * ta2_xx_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xyz_1[i] * fe_0 + ta2_xx_zz_xyyz_0[i] * pa_y[i] - ta2_xx_zz_xyyz_1[i] * pc_y[i];

        ta2_xx_yzz_xyzz_0[i] = ta2_xx_zz_xzz_0[i] * fe_0 - ta2_xx_zz_xzz_1[i] * fe_0 + ta2_xx_zz_xyzz_0[i] * pa_y[i] - ta2_xx_zz_xyzz_1[i] * pc_y[i];

        ta2_xx_yzz_xzzz_0[i] = ta2_xx_zz_xzzz_0[i] * pa_y[i] - ta2_xx_zz_xzzz_1[i] * pc_y[i];

        ta2_xx_yzz_yyyy_0[i] = 4.0 * ta2_xx_zz_yyy_0[i] * fe_0 - 4.0 * ta2_xx_zz_yyy_1[i] * fe_0 + ta2_xx_zz_yyyy_0[i] * pa_y[i] - ta2_xx_zz_yyyy_1[i] * pc_y[i];

        ta2_xx_yzz_yyyz_0[i] = 3.0 * ta2_xx_zz_yyz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yyz_1[i] * fe_0 + ta2_xx_zz_yyyz_0[i] * pa_y[i] - ta2_xx_zz_yyyz_1[i] * pc_y[i];

        ta2_xx_yzz_yyzz_0[i] = 2.0 * ta2_xx_zz_yzz_0[i] * fe_0 - 2.0 * ta2_xx_zz_yzz_1[i] * fe_0 + ta2_xx_zz_yyzz_0[i] * pa_y[i] - ta2_xx_zz_yyzz_1[i] * pc_y[i];

        ta2_xx_yzz_yzzz_0[i] = ta2_xx_zz_zzz_0[i] * fe_0 - ta2_xx_zz_zzz_1[i] * fe_0 + ta2_xx_zz_yzzz_0[i] * pa_y[i] - ta2_xx_zz_yzzz_1[i] * pc_y[i];

        ta2_xx_yzz_zzzz_0[i] = ta2_xx_zz_zzzz_0[i] * pa_y[i] - ta2_xx_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto ta2_xx_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 135);

    auto ta2_xx_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 136);

    auto ta2_xx_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 137);

    auto ta2_xx_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 138);

    auto ta2_xx_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 139);

    auto ta2_xx_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 140);

    auto ta2_xx_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 141);

    auto ta2_xx_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 142);

    auto ta2_xx_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 143);

    auto ta2_xx_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 144);

    auto ta2_xx_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 145);

    auto ta2_xx_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 146);

    auto ta2_xx_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 147);

    auto ta2_xx_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 148);

    auto ta2_xx_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 149);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_z_xxxx_0, ta2_xx_z_xxxx_1, ta2_xx_z_xxxy_0, ta2_xx_z_xxxy_1, ta2_xx_z_xxxz_0, ta2_xx_z_xxxz_1, ta2_xx_z_xxyy_0, ta2_xx_z_xxyy_1, ta2_xx_z_xxyz_0, ta2_xx_z_xxyz_1, ta2_xx_z_xxzz_0, ta2_xx_z_xxzz_1, ta2_xx_z_xyyy_0, ta2_xx_z_xyyy_1, ta2_xx_z_xyyz_0, ta2_xx_z_xyyz_1, ta2_xx_z_xyzz_0, ta2_xx_z_xyzz_1, ta2_xx_z_xzzz_0, ta2_xx_z_xzzz_1, ta2_xx_z_yyyy_0, ta2_xx_z_yyyy_1, ta2_xx_z_yyyz_0, ta2_xx_z_yyyz_1, ta2_xx_z_yyzz_0, ta2_xx_z_yyzz_1, ta2_xx_z_yzzz_0, ta2_xx_z_yzzz_1, ta2_xx_z_zzzz_0, ta2_xx_z_zzzz_1, ta2_xx_zz_xxx_0, ta2_xx_zz_xxx_1, ta2_xx_zz_xxxx_0, ta2_xx_zz_xxxx_1, ta2_xx_zz_xxxy_0, ta2_xx_zz_xxxy_1, ta2_xx_zz_xxxz_0, ta2_xx_zz_xxxz_1, ta2_xx_zz_xxy_0, ta2_xx_zz_xxy_1, ta2_xx_zz_xxyy_0, ta2_xx_zz_xxyy_1, ta2_xx_zz_xxyz_0, ta2_xx_zz_xxyz_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xxzz_0, ta2_xx_zz_xxzz_1, ta2_xx_zz_xyy_0, ta2_xx_zz_xyy_1, ta2_xx_zz_xyyy_0, ta2_xx_zz_xyyy_1, ta2_xx_zz_xyyz_0, ta2_xx_zz_xyyz_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xyzz_0, ta2_xx_zz_xyzz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_xzzz_0, ta2_xx_zz_xzzz_1, ta2_xx_zz_yyy_0, ta2_xx_zz_yyy_1, ta2_xx_zz_yyyy_0, ta2_xx_zz_yyyy_1, ta2_xx_zz_yyyz_0, ta2_xx_zz_yyyz_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yyzz_0, ta2_xx_zz_yyzz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_yzzz_0, ta2_xx_zz_yzzz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, ta2_xx_zz_zzzz_0, ta2_xx_zz_zzzz_1, ta2_xx_zzz_xxxx_0, ta2_xx_zzz_xxxy_0, ta2_xx_zzz_xxxz_0, ta2_xx_zzz_xxyy_0, ta2_xx_zzz_xxyz_0, ta2_xx_zzz_xxzz_0, ta2_xx_zzz_xyyy_0, ta2_xx_zzz_xyyz_0, ta2_xx_zzz_xyzz_0, ta2_xx_zzz_xzzz_0, ta2_xx_zzz_yyyy_0, ta2_xx_zzz_yyyz_0, ta2_xx_zzz_yyzz_0, ta2_xx_zzz_yzzz_0, ta2_xx_zzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzz_xxxx_0[i] = 2.0 * ta2_xx_z_xxxx_0[i] * fe_0 - 2.0 * ta2_xx_z_xxxx_1[i] * fe_0 + ta2_xx_zz_xxxx_0[i] * pa_z[i] - ta2_xx_zz_xxxx_1[i] * pc_z[i];

        ta2_xx_zzz_xxxy_0[i] = 2.0 * ta2_xx_z_xxxy_0[i] * fe_0 - 2.0 * ta2_xx_z_xxxy_1[i] * fe_0 + ta2_xx_zz_xxxy_0[i] * pa_z[i] - ta2_xx_zz_xxxy_1[i] * pc_z[i];

        ta2_xx_zzz_xxxz_0[i] = 2.0 * ta2_xx_z_xxxz_0[i] * fe_0 - 2.0 * ta2_xx_z_xxxz_1[i] * fe_0 + ta2_xx_zz_xxx_0[i] * fe_0 - ta2_xx_zz_xxx_1[i] * fe_0 + ta2_xx_zz_xxxz_0[i] * pa_z[i] - ta2_xx_zz_xxxz_1[i] * pc_z[i];

        ta2_xx_zzz_xxyy_0[i] = 2.0 * ta2_xx_z_xxyy_0[i] * fe_0 - 2.0 * ta2_xx_z_xxyy_1[i] * fe_0 + ta2_xx_zz_xxyy_0[i] * pa_z[i] - ta2_xx_zz_xxyy_1[i] * pc_z[i];

        ta2_xx_zzz_xxyz_0[i] = 2.0 * ta2_xx_z_xxyz_0[i] * fe_0 - 2.0 * ta2_xx_z_xxyz_1[i] * fe_0 + ta2_xx_zz_xxy_0[i] * fe_0 - ta2_xx_zz_xxy_1[i] * fe_0 + ta2_xx_zz_xxyz_0[i] * pa_z[i] - ta2_xx_zz_xxyz_1[i] * pc_z[i];

        ta2_xx_zzz_xxzz_0[i] = 2.0 * ta2_xx_z_xxzz_0[i] * fe_0 - 2.0 * ta2_xx_z_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_zz_xxz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xxz_1[i] * fe_0 + ta2_xx_zz_xxzz_0[i] * pa_z[i] - ta2_xx_zz_xxzz_1[i] * pc_z[i];

        ta2_xx_zzz_xyyy_0[i] = 2.0 * ta2_xx_z_xyyy_0[i] * fe_0 - 2.0 * ta2_xx_z_xyyy_1[i] * fe_0 + ta2_xx_zz_xyyy_0[i] * pa_z[i] - ta2_xx_zz_xyyy_1[i] * pc_z[i];

        ta2_xx_zzz_xyyz_0[i] = 2.0 * ta2_xx_z_xyyz_0[i] * fe_0 - 2.0 * ta2_xx_z_xyyz_1[i] * fe_0 + ta2_xx_zz_xyy_0[i] * fe_0 - ta2_xx_zz_xyy_1[i] * fe_0 + ta2_xx_zz_xyyz_0[i] * pa_z[i] - ta2_xx_zz_xyyz_1[i] * pc_z[i];

        ta2_xx_zzz_xyzz_0[i] = 2.0 * ta2_xx_z_xyzz_0[i] * fe_0 - 2.0 * ta2_xx_z_xyzz_1[i] * fe_0 + 2.0 * ta2_xx_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xyz_1[i] * fe_0 + ta2_xx_zz_xyzz_0[i] * pa_z[i] - ta2_xx_zz_xyzz_1[i] * pc_z[i];

        ta2_xx_zzz_xzzz_0[i] = 2.0 * ta2_xx_z_xzzz_0[i] * fe_0 - 2.0 * ta2_xx_z_xzzz_1[i] * fe_0 + 3.0 * ta2_xx_zz_xzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xzz_1[i] * fe_0 + ta2_xx_zz_xzzz_0[i] * pa_z[i] - ta2_xx_zz_xzzz_1[i] * pc_z[i];

        ta2_xx_zzz_yyyy_0[i] = 2.0 * ta2_xx_z_yyyy_0[i] * fe_0 - 2.0 * ta2_xx_z_yyyy_1[i] * fe_0 + ta2_xx_zz_yyyy_0[i] * pa_z[i] - ta2_xx_zz_yyyy_1[i] * pc_z[i];

        ta2_xx_zzz_yyyz_0[i] = 2.0 * ta2_xx_z_yyyz_0[i] * fe_0 - 2.0 * ta2_xx_z_yyyz_1[i] * fe_0 + ta2_xx_zz_yyy_0[i] * fe_0 - ta2_xx_zz_yyy_1[i] * fe_0 + ta2_xx_zz_yyyz_0[i] * pa_z[i] - ta2_xx_zz_yyyz_1[i] * pc_z[i];

        ta2_xx_zzz_yyzz_0[i] = 2.0 * ta2_xx_z_yyzz_0[i] * fe_0 - 2.0 * ta2_xx_z_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_zz_yyz_0[i] * fe_0 - 2.0 * ta2_xx_zz_yyz_1[i] * fe_0 + ta2_xx_zz_yyzz_0[i] * pa_z[i] - ta2_xx_zz_yyzz_1[i] * pc_z[i];

        ta2_xx_zzz_yzzz_0[i] = 2.0 * ta2_xx_z_yzzz_0[i] * fe_0 - 2.0 * ta2_xx_z_yzzz_1[i] * fe_0 + 3.0 * ta2_xx_zz_yzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yzz_1[i] * fe_0 + ta2_xx_zz_yzzz_0[i] * pa_z[i] - ta2_xx_zz_yzzz_1[i] * pc_z[i];

        ta2_xx_zzz_zzzz_0[i] = 2.0 * ta2_xx_z_zzzz_0[i] * fe_0 - 2.0 * ta2_xx_z_zzzz_1[i] * fe_0 + 4.0 * ta2_xx_zz_zzz_0[i] * fe_0 - 4.0 * ta2_xx_zz_zzz_1[i] * fe_0 + ta2_xx_zz_zzzz_0[i] * pa_z[i] - ta2_xx_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 150-165 components of targeted buffer : FG

    auto ta2_xy_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 150);

    auto ta2_xy_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 151);

    auto ta2_xy_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 152);

    auto ta2_xy_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 153);

    auto ta2_xy_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 154);

    auto ta2_xy_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 155);

    auto ta2_xy_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 156);

    auto ta2_xy_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 157);

    auto ta2_xy_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 158);

    auto ta2_xy_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 159);

    auto ta2_xy_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 160);

    auto ta2_xy_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 161);

    auto ta2_xy_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 162);

    auto ta2_xy_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 163);

    auto ta2_xy_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 164);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xx_xxxx_1, ta1_y_xx_xxxy_1, ta1_y_xx_xxxz_1, ta1_y_xx_xxyy_1, ta1_y_xx_xxyz_1, ta1_y_xx_xxzz_1, ta1_y_xx_xyyy_1, ta1_y_xx_xyyz_1, ta1_y_xx_xyzz_1, ta1_y_xx_xzzz_1, ta1_y_xx_yyyy_1, ta1_y_xx_yyyz_1, ta1_y_xx_yyzz_1, ta1_y_xx_yzzz_1, ta1_y_xx_zzzz_1, ta2_xy_x_xxxx_0, ta2_xy_x_xxxx_1, ta2_xy_x_xxxy_0, ta2_xy_x_xxxy_1, ta2_xy_x_xxxz_0, ta2_xy_x_xxxz_1, ta2_xy_x_xxyy_0, ta2_xy_x_xxyy_1, ta2_xy_x_xxyz_0, ta2_xy_x_xxyz_1, ta2_xy_x_xxzz_0, ta2_xy_x_xxzz_1, ta2_xy_x_xyyy_0, ta2_xy_x_xyyy_1, ta2_xy_x_xyyz_0, ta2_xy_x_xyyz_1, ta2_xy_x_xyzz_0, ta2_xy_x_xyzz_1, ta2_xy_x_xzzz_0, ta2_xy_x_xzzz_1, ta2_xy_x_yyyy_0, ta2_xy_x_yyyy_1, ta2_xy_x_yyyz_0, ta2_xy_x_yyyz_1, ta2_xy_x_yyzz_0, ta2_xy_x_yyzz_1, ta2_xy_x_yzzz_0, ta2_xy_x_yzzz_1, ta2_xy_x_zzzz_0, ta2_xy_x_zzzz_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxxx_0, ta2_xy_xx_xxxx_1, ta2_xy_xx_xxxy_0, ta2_xy_xx_xxxy_1, ta2_xy_xx_xxxz_0, ta2_xy_xx_xxxz_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxyy_0, ta2_xy_xx_xxyy_1, ta2_xy_xx_xxyz_0, ta2_xy_xx_xxyz_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xxzz_0, ta2_xy_xx_xxzz_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyyy_0, ta2_xy_xx_xyyy_1, ta2_xy_xx_xyyz_0, ta2_xy_xx_xyyz_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xyzz_0, ta2_xy_xx_xyzz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_xzzz_0, ta2_xy_xx_xzzz_1, ta2_xy_xx_yyy_0, ta2_xy_xx_yyy_1, ta2_xy_xx_yyyy_0, ta2_xy_xx_yyyy_1, ta2_xy_xx_yyyz_0, ta2_xy_xx_yyyz_1, ta2_xy_xx_yyz_0, ta2_xy_xx_yyz_1, ta2_xy_xx_yyzz_0, ta2_xy_xx_yyzz_1, ta2_xy_xx_yzz_0, ta2_xy_xx_yzz_1, ta2_xy_xx_yzzz_0, ta2_xy_xx_yzzz_1, ta2_xy_xx_zzz_0, ta2_xy_xx_zzz_1, ta2_xy_xx_zzzz_0, ta2_xy_xx_zzzz_1, ta2_xy_xxx_xxxx_0, ta2_xy_xxx_xxxy_0, ta2_xy_xxx_xxxz_0, ta2_xy_xxx_xxyy_0, ta2_xy_xxx_xxyz_0, ta2_xy_xxx_xxzz_0, ta2_xy_xxx_xyyy_0, ta2_xy_xxx_xyyz_0, ta2_xy_xxx_xyzz_0, ta2_xy_xxx_xzzz_0, ta2_xy_xxx_yyyy_0, ta2_xy_xxx_yyyz_0, ta2_xy_xxx_yyzz_0, ta2_xy_xxx_yzzz_0, ta2_xy_xxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxx_xxxx_0[i] = 2.0 * ta2_xy_x_xxxx_0[i] * fe_0 - 2.0 * ta2_xy_x_xxxx_1[i] * fe_0 + 4.0 * ta2_xy_xx_xxx_0[i] * fe_0 - 4.0 * ta2_xy_xx_xxx_1[i] * fe_0 + ta1_y_xx_xxxx_1[i] + ta2_xy_xx_xxxx_0[i] * pa_x[i] - ta2_xy_xx_xxxx_1[i] * pc_x[i];

        ta2_xy_xxx_xxxy_0[i] = 2.0 * ta2_xy_x_xxxy_0[i] * fe_0 - 2.0 * ta2_xy_x_xxxy_1[i] * fe_0 + 3.0 * ta2_xy_xx_xxy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxy_1[i] * fe_0 + ta1_y_xx_xxxy_1[i] + ta2_xy_xx_xxxy_0[i] * pa_x[i] - ta2_xy_xx_xxxy_1[i] * pc_x[i];

        ta2_xy_xxx_xxxz_0[i] = 2.0 * ta2_xy_x_xxxz_0[i] * fe_0 - 2.0 * ta2_xy_x_xxxz_1[i] * fe_0 + 3.0 * ta2_xy_xx_xxz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxz_1[i] * fe_0 + ta1_y_xx_xxxz_1[i] + ta2_xy_xx_xxxz_0[i] * pa_x[i] - ta2_xy_xx_xxxz_1[i] * pc_x[i];

        ta2_xy_xxx_xxyy_0[i] = 2.0 * ta2_xy_x_xxyy_0[i] * fe_0 - 2.0 * ta2_xy_x_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_xx_xyy_0[i] * fe_0 - 2.0 * ta2_xy_xx_xyy_1[i] * fe_0 + ta1_y_xx_xxyy_1[i] + ta2_xy_xx_xxyy_0[i] * pa_x[i] - ta2_xy_xx_xxyy_1[i] * pc_x[i];

        ta2_xy_xxx_xxyz_0[i] = 2.0 * ta2_xy_x_xxyz_0[i] * fe_0 - 2.0 * ta2_xy_x_xxyz_1[i] * fe_0 + 2.0 * ta2_xy_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xyz_1[i] * fe_0 + ta1_y_xx_xxyz_1[i] + ta2_xy_xx_xxyz_0[i] * pa_x[i] - ta2_xy_xx_xxyz_1[i] * pc_x[i];

        ta2_xy_xxx_xxzz_0[i] = 2.0 * ta2_xy_x_xxzz_0[i] * fe_0 - 2.0 * ta2_xy_x_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_xx_xzz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xzz_1[i] * fe_0 + ta1_y_xx_xxzz_1[i] + ta2_xy_xx_xxzz_0[i] * pa_x[i] - ta2_xy_xx_xxzz_1[i] * pc_x[i];

        ta2_xy_xxx_xyyy_0[i] = 2.0 * ta2_xy_x_xyyy_0[i] * fe_0 - 2.0 * ta2_xy_x_xyyy_1[i] * fe_0 + ta2_xy_xx_yyy_0[i] * fe_0 - ta2_xy_xx_yyy_1[i] * fe_0 + ta1_y_xx_xyyy_1[i] + ta2_xy_xx_xyyy_0[i] * pa_x[i] - ta2_xy_xx_xyyy_1[i] * pc_x[i];

        ta2_xy_xxx_xyyz_0[i] = 2.0 * ta2_xy_x_xyyz_0[i] * fe_0 - 2.0 * ta2_xy_x_xyyz_1[i] * fe_0 + ta2_xy_xx_yyz_0[i] * fe_0 - ta2_xy_xx_yyz_1[i] * fe_0 + ta1_y_xx_xyyz_1[i] + ta2_xy_xx_xyyz_0[i] * pa_x[i] - ta2_xy_xx_xyyz_1[i] * pc_x[i];

        ta2_xy_xxx_xyzz_0[i] = 2.0 * ta2_xy_x_xyzz_0[i] * fe_0 - 2.0 * ta2_xy_x_xyzz_1[i] * fe_0 + ta2_xy_xx_yzz_0[i] * fe_0 - ta2_xy_xx_yzz_1[i] * fe_0 + ta1_y_xx_xyzz_1[i] + ta2_xy_xx_xyzz_0[i] * pa_x[i] - ta2_xy_xx_xyzz_1[i] * pc_x[i];

        ta2_xy_xxx_xzzz_0[i] = 2.0 * ta2_xy_x_xzzz_0[i] * fe_0 - 2.0 * ta2_xy_x_xzzz_1[i] * fe_0 + ta2_xy_xx_zzz_0[i] * fe_0 - ta2_xy_xx_zzz_1[i] * fe_0 + ta1_y_xx_xzzz_1[i] + ta2_xy_xx_xzzz_0[i] * pa_x[i] - ta2_xy_xx_xzzz_1[i] * pc_x[i];

        ta2_xy_xxx_yyyy_0[i] = 2.0 * ta2_xy_x_yyyy_0[i] * fe_0 - 2.0 * ta2_xy_x_yyyy_1[i] * fe_0 + ta1_y_xx_yyyy_1[i] + ta2_xy_xx_yyyy_0[i] * pa_x[i] - ta2_xy_xx_yyyy_1[i] * pc_x[i];

        ta2_xy_xxx_yyyz_0[i] = 2.0 * ta2_xy_x_yyyz_0[i] * fe_0 - 2.0 * ta2_xy_x_yyyz_1[i] * fe_0 + ta1_y_xx_yyyz_1[i] + ta2_xy_xx_yyyz_0[i] * pa_x[i] - ta2_xy_xx_yyyz_1[i] * pc_x[i];

        ta2_xy_xxx_yyzz_0[i] = 2.0 * ta2_xy_x_yyzz_0[i] * fe_0 - 2.0 * ta2_xy_x_yyzz_1[i] * fe_0 + ta1_y_xx_yyzz_1[i] + ta2_xy_xx_yyzz_0[i] * pa_x[i] - ta2_xy_xx_yyzz_1[i] * pc_x[i];

        ta2_xy_xxx_yzzz_0[i] = 2.0 * ta2_xy_x_yzzz_0[i] * fe_0 - 2.0 * ta2_xy_x_yzzz_1[i] * fe_0 + ta1_y_xx_yzzz_1[i] + ta2_xy_xx_yzzz_0[i] * pa_x[i] - ta2_xy_xx_yzzz_1[i] * pc_x[i];

        ta2_xy_xxx_zzzz_0[i] = 2.0 * ta2_xy_x_zzzz_0[i] * fe_0 - 2.0 * ta2_xy_x_zzzz_1[i] * fe_0 + ta1_y_xx_zzzz_1[i] + ta2_xy_xx_zzzz_0[i] * pa_x[i] - ta2_xy_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 165-180 components of targeted buffer : FG

    auto ta2_xy_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 165);

    auto ta2_xy_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 166);

    auto ta2_xy_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 167);

    auto ta2_xy_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 168);

    auto ta2_xy_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 169);

    auto ta2_xy_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 170);

    auto ta2_xy_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 171);

    auto ta2_xy_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 172);

    auto ta2_xy_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 173);

    auto ta2_xy_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 174);

    auto ta2_xy_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 175);

    auto ta2_xy_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 176);

    auto ta2_xy_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 177);

    auto ta2_xy_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 178);

    auto ta2_xy_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 179);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xx_xxxx_1, ta1_x_xx_xxxy_1, ta1_x_xx_xxxz_1, ta1_x_xx_xxyy_1, ta1_x_xx_xxyz_1, ta1_x_xx_xxzz_1, ta1_x_xx_xyyy_1, ta1_x_xx_xyyz_1, ta1_x_xx_xyzz_1, ta1_x_xx_xzzz_1, ta1_x_xx_zzzz_1, ta1_y_xy_yyyy_1, ta1_y_xy_yyyz_1, ta1_y_xy_yyzz_1, ta1_y_xy_yzzz_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxxx_0, ta2_xy_xx_xxxx_1, ta2_xy_xx_xxxy_0, ta2_xy_xx_xxxy_1, ta2_xy_xx_xxxz_0, ta2_xy_xx_xxxz_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxyy_0, ta2_xy_xx_xxyy_1, ta2_xy_xx_xxyz_0, ta2_xy_xx_xxyz_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xxzz_0, ta2_xy_xx_xxzz_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyyy_0, ta2_xy_xx_xyyy_1, ta2_xy_xx_xyyz_0, ta2_xy_xx_xyyz_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xyzz_0, ta2_xy_xx_xyzz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_xzzz_0, ta2_xy_xx_xzzz_1, ta2_xy_xx_zzzz_0, ta2_xy_xx_zzzz_1, ta2_xy_xxy_xxxx_0, ta2_xy_xxy_xxxy_0, ta2_xy_xxy_xxxz_0, ta2_xy_xxy_xxyy_0, ta2_xy_xxy_xxyz_0, ta2_xy_xxy_xxzz_0, ta2_xy_xxy_xyyy_0, ta2_xy_xxy_xyyz_0, ta2_xy_xxy_xyzz_0, ta2_xy_xxy_xzzz_0, ta2_xy_xxy_yyyy_0, ta2_xy_xxy_yyyz_0, ta2_xy_xxy_yyzz_0, ta2_xy_xxy_yzzz_0, ta2_xy_xxy_zzzz_0, ta2_xy_xy_yyyy_0, ta2_xy_xy_yyyy_1, ta2_xy_xy_yyyz_0, ta2_xy_xy_yyyz_1, ta2_xy_xy_yyzz_0, ta2_xy_xy_yyzz_1, ta2_xy_xy_yzzz_0, ta2_xy_xy_yzzz_1, ta2_xy_y_yyyy_0, ta2_xy_y_yyyy_1, ta2_xy_y_yyyz_0, ta2_xy_y_yyyz_1, ta2_xy_y_yyzz_0, ta2_xy_y_yyzz_1, ta2_xy_y_yzzz_0, ta2_xy_y_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxy_xxxx_0[i] = ta1_x_xx_xxxx_1[i] + ta2_xy_xx_xxxx_0[i] * pa_y[i] - ta2_xy_xx_xxxx_1[i] * pc_y[i];

        ta2_xy_xxy_xxxy_0[i] = ta2_xy_xx_xxx_0[i] * fe_0 - ta2_xy_xx_xxx_1[i] * fe_0 + ta1_x_xx_xxxy_1[i] + ta2_xy_xx_xxxy_0[i] * pa_y[i] - ta2_xy_xx_xxxy_1[i] * pc_y[i];

        ta2_xy_xxy_xxxz_0[i] = ta1_x_xx_xxxz_1[i] + ta2_xy_xx_xxxz_0[i] * pa_y[i] - ta2_xy_xx_xxxz_1[i] * pc_y[i];

        ta2_xy_xxy_xxyy_0[i] = 2.0 * ta2_xy_xx_xxy_0[i] * fe_0 - 2.0 * ta2_xy_xx_xxy_1[i] * fe_0 + ta1_x_xx_xxyy_1[i] + ta2_xy_xx_xxyy_0[i] * pa_y[i] - ta2_xy_xx_xxyy_1[i] * pc_y[i];

        ta2_xy_xxy_xxyz_0[i] = ta2_xy_xx_xxz_0[i] * fe_0 - ta2_xy_xx_xxz_1[i] * fe_0 + ta1_x_xx_xxyz_1[i] + ta2_xy_xx_xxyz_0[i] * pa_y[i] - ta2_xy_xx_xxyz_1[i] * pc_y[i];

        ta2_xy_xxy_xxzz_0[i] = ta1_x_xx_xxzz_1[i] + ta2_xy_xx_xxzz_0[i] * pa_y[i] - ta2_xy_xx_xxzz_1[i] * pc_y[i];

        ta2_xy_xxy_xyyy_0[i] = 3.0 * ta2_xy_xx_xyy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xyy_1[i] * fe_0 + ta1_x_xx_xyyy_1[i] + ta2_xy_xx_xyyy_0[i] * pa_y[i] - ta2_xy_xx_xyyy_1[i] * pc_y[i];

        ta2_xy_xxy_xyyz_0[i] = 2.0 * ta2_xy_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xyz_1[i] * fe_0 + ta1_x_xx_xyyz_1[i] + ta2_xy_xx_xyyz_0[i] * pa_y[i] - ta2_xy_xx_xyyz_1[i] * pc_y[i];

        ta2_xy_xxy_xyzz_0[i] = ta2_xy_xx_xzz_0[i] * fe_0 - ta2_xy_xx_xzz_1[i] * fe_0 + ta1_x_xx_xyzz_1[i] + ta2_xy_xx_xyzz_0[i] * pa_y[i] - ta2_xy_xx_xyzz_1[i] * pc_y[i];

        ta2_xy_xxy_xzzz_0[i] = ta1_x_xx_xzzz_1[i] + ta2_xy_xx_xzzz_0[i] * pa_y[i] - ta2_xy_xx_xzzz_1[i] * pc_y[i];

        ta2_xy_xxy_yyyy_0[i] = ta2_xy_y_yyyy_0[i] * fe_0 - ta2_xy_y_yyyy_1[i] * fe_0 + ta1_y_xy_yyyy_1[i] + ta2_xy_xy_yyyy_0[i] * pa_x[i] - ta2_xy_xy_yyyy_1[i] * pc_x[i];

        ta2_xy_xxy_yyyz_0[i] = ta2_xy_y_yyyz_0[i] * fe_0 - ta2_xy_y_yyyz_1[i] * fe_0 + ta1_y_xy_yyyz_1[i] + ta2_xy_xy_yyyz_0[i] * pa_x[i] - ta2_xy_xy_yyyz_1[i] * pc_x[i];

        ta2_xy_xxy_yyzz_0[i] = ta2_xy_y_yyzz_0[i] * fe_0 - ta2_xy_y_yyzz_1[i] * fe_0 + ta1_y_xy_yyzz_1[i] + ta2_xy_xy_yyzz_0[i] * pa_x[i] - ta2_xy_xy_yyzz_1[i] * pc_x[i];

        ta2_xy_xxy_yzzz_0[i] = ta2_xy_y_yzzz_0[i] * fe_0 - ta2_xy_y_yzzz_1[i] * fe_0 + ta1_y_xy_yzzz_1[i] + ta2_xy_xy_yzzz_0[i] * pa_x[i] - ta2_xy_xy_yzzz_1[i] * pc_x[i];

        ta2_xy_xxy_zzzz_0[i] = ta1_x_xx_zzzz_1[i] + ta2_xy_xx_zzzz_0[i] * pa_y[i] - ta2_xy_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 180-195 components of targeted buffer : FG

    auto ta2_xy_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 180);

    auto ta2_xy_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 181);

    auto ta2_xy_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 182);

    auto ta2_xy_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 183);

    auto ta2_xy_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 184);

    auto ta2_xy_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 185);

    auto ta2_xy_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 186);

    auto ta2_xy_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 187);

    auto ta2_xy_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 188);

    auto ta2_xy_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 189);

    auto ta2_xy_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 190);

    auto ta2_xy_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 191);

    auto ta2_xy_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 192);

    auto ta2_xy_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 193);

    auto ta2_xy_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 194);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxxx_0, ta2_xy_xx_xxxx_1, ta2_xy_xx_xxxy_0, ta2_xy_xx_xxxy_1, ta2_xy_xx_xxxz_0, ta2_xy_xx_xxxz_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxyy_0, ta2_xy_xx_xxyy_1, ta2_xy_xx_xxyz_0, ta2_xy_xx_xxyz_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xxzz_0, ta2_xy_xx_xxzz_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyyy_0, ta2_xy_xx_xyyy_1, ta2_xy_xx_xyyz_0, ta2_xy_xx_xyyz_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xyzz_0, ta2_xy_xx_xyzz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_xzzz_0, ta2_xy_xx_xzzz_1, ta2_xy_xx_yyy_0, ta2_xy_xx_yyy_1, ta2_xy_xx_yyyy_0, ta2_xy_xx_yyyy_1, ta2_xy_xx_yyyz_0, ta2_xy_xx_yyyz_1, ta2_xy_xx_yyz_0, ta2_xy_xx_yyz_1, ta2_xy_xx_yyzz_0, ta2_xy_xx_yyzz_1, ta2_xy_xx_yzz_0, ta2_xy_xx_yzz_1, ta2_xy_xx_yzzz_0, ta2_xy_xx_yzzz_1, ta2_xy_xx_zzz_0, ta2_xy_xx_zzz_1, ta2_xy_xx_zzzz_0, ta2_xy_xx_zzzz_1, ta2_xy_xxz_xxxx_0, ta2_xy_xxz_xxxy_0, ta2_xy_xxz_xxxz_0, ta2_xy_xxz_xxyy_0, ta2_xy_xxz_xxyz_0, ta2_xy_xxz_xxzz_0, ta2_xy_xxz_xyyy_0, ta2_xy_xxz_xyyz_0, ta2_xy_xxz_xyzz_0, ta2_xy_xxz_xzzz_0, ta2_xy_xxz_yyyy_0, ta2_xy_xxz_yyyz_0, ta2_xy_xxz_yyzz_0, ta2_xy_xxz_yzzz_0, ta2_xy_xxz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxz_xxxx_0[i] = ta2_xy_xx_xxxx_0[i] * pa_z[i] - ta2_xy_xx_xxxx_1[i] * pc_z[i];

        ta2_xy_xxz_xxxy_0[i] = ta2_xy_xx_xxxy_0[i] * pa_z[i] - ta2_xy_xx_xxxy_1[i] * pc_z[i];

        ta2_xy_xxz_xxxz_0[i] = ta2_xy_xx_xxx_0[i] * fe_0 - ta2_xy_xx_xxx_1[i] * fe_0 + ta2_xy_xx_xxxz_0[i] * pa_z[i] - ta2_xy_xx_xxxz_1[i] * pc_z[i];

        ta2_xy_xxz_xxyy_0[i] = ta2_xy_xx_xxyy_0[i] * pa_z[i] - ta2_xy_xx_xxyy_1[i] * pc_z[i];

        ta2_xy_xxz_xxyz_0[i] = ta2_xy_xx_xxy_0[i] * fe_0 - ta2_xy_xx_xxy_1[i] * fe_0 + ta2_xy_xx_xxyz_0[i] * pa_z[i] - ta2_xy_xx_xxyz_1[i] * pc_z[i];

        ta2_xy_xxz_xxzz_0[i] = 2.0 * ta2_xy_xx_xxz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xxz_1[i] * fe_0 + ta2_xy_xx_xxzz_0[i] * pa_z[i] - ta2_xy_xx_xxzz_1[i] * pc_z[i];

        ta2_xy_xxz_xyyy_0[i] = ta2_xy_xx_xyyy_0[i] * pa_z[i] - ta2_xy_xx_xyyy_1[i] * pc_z[i];

        ta2_xy_xxz_xyyz_0[i] = ta2_xy_xx_xyy_0[i] * fe_0 - ta2_xy_xx_xyy_1[i] * fe_0 + ta2_xy_xx_xyyz_0[i] * pa_z[i] - ta2_xy_xx_xyyz_1[i] * pc_z[i];

        ta2_xy_xxz_xyzz_0[i] = 2.0 * ta2_xy_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xyz_1[i] * fe_0 + ta2_xy_xx_xyzz_0[i] * pa_z[i] - ta2_xy_xx_xyzz_1[i] * pc_z[i];

        ta2_xy_xxz_xzzz_0[i] = 3.0 * ta2_xy_xx_xzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xzz_1[i] * fe_0 + ta2_xy_xx_xzzz_0[i] * pa_z[i] - ta2_xy_xx_xzzz_1[i] * pc_z[i];

        ta2_xy_xxz_yyyy_0[i] = ta2_xy_xx_yyyy_0[i] * pa_z[i] - ta2_xy_xx_yyyy_1[i] * pc_z[i];

        ta2_xy_xxz_yyyz_0[i] = ta2_xy_xx_yyy_0[i] * fe_0 - ta2_xy_xx_yyy_1[i] * fe_0 + ta2_xy_xx_yyyz_0[i] * pa_z[i] - ta2_xy_xx_yyyz_1[i] * pc_z[i];

        ta2_xy_xxz_yyzz_0[i] = 2.0 * ta2_xy_xx_yyz_0[i] * fe_0 - 2.0 * ta2_xy_xx_yyz_1[i] * fe_0 + ta2_xy_xx_yyzz_0[i] * pa_z[i] - ta2_xy_xx_yyzz_1[i] * pc_z[i];

        ta2_xy_xxz_yzzz_0[i] = 3.0 * ta2_xy_xx_yzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yzz_1[i] * fe_0 + ta2_xy_xx_yzzz_0[i] * pa_z[i] - ta2_xy_xx_yzzz_1[i] * pc_z[i];

        ta2_xy_xxz_zzzz_0[i] = 4.0 * ta2_xy_xx_zzz_0[i] * fe_0 - 4.0 * ta2_xy_xx_zzz_1[i] * fe_0 + ta2_xy_xx_zzzz_0[i] * pa_z[i] - ta2_xy_xx_zzzz_1[i] * pc_z[i];
    }

    // Set up 195-210 components of targeted buffer : FG

    auto ta2_xy_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 195);

    auto ta2_xy_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 196);

    auto ta2_xy_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 197);

    auto ta2_xy_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 198);

    auto ta2_xy_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 199);

    auto ta2_xy_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 200);

    auto ta2_xy_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 201);

    auto ta2_xy_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 202);

    auto ta2_xy_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 203);

    auto ta2_xy_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 204);

    auto ta2_xy_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 205);

    auto ta2_xy_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 206);

    auto ta2_xy_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 207);

    auto ta2_xy_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 208);

    auto ta2_xy_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 209);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_yy_xxxx_1, ta1_y_yy_xxxy_1, ta1_y_yy_xxxz_1, ta1_y_yy_xxyy_1, ta1_y_yy_xxyz_1, ta1_y_yy_xxzz_1, ta1_y_yy_xyyy_1, ta1_y_yy_xyyz_1, ta1_y_yy_xyzz_1, ta1_y_yy_xzzz_1, ta1_y_yy_yyyy_1, ta1_y_yy_yyyz_1, ta1_y_yy_yyzz_1, ta1_y_yy_yzzz_1, ta1_y_yy_zzzz_1, ta2_xy_xyy_xxxx_0, ta2_xy_xyy_xxxy_0, ta2_xy_xyy_xxxz_0, ta2_xy_xyy_xxyy_0, ta2_xy_xyy_xxyz_0, ta2_xy_xyy_xxzz_0, ta2_xy_xyy_xyyy_0, ta2_xy_xyy_xyyz_0, ta2_xy_xyy_xyzz_0, ta2_xy_xyy_xzzz_0, ta2_xy_xyy_yyyy_0, ta2_xy_xyy_yyyz_0, ta2_xy_xyy_yyzz_0, ta2_xy_xyy_yzzz_0, ta2_xy_xyy_zzzz_0, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxxx_0, ta2_xy_yy_xxxx_1, ta2_xy_yy_xxxy_0, ta2_xy_yy_xxxy_1, ta2_xy_yy_xxxz_0, ta2_xy_yy_xxxz_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxyy_0, ta2_xy_yy_xxyy_1, ta2_xy_yy_xxyz_0, ta2_xy_yy_xxyz_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xxzz_0, ta2_xy_yy_xxzz_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyyy_0, ta2_xy_yy_xyyy_1, ta2_xy_yy_xyyz_0, ta2_xy_yy_xyyz_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xyzz_0, ta2_xy_yy_xyzz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_xzzz_0, ta2_xy_yy_xzzz_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyyy_0, ta2_xy_yy_yyyy_1, ta2_xy_yy_yyyz_0, ta2_xy_yy_yyyz_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yyzz_0, ta2_xy_yy_yyzz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_yzzz_0, ta2_xy_yy_yzzz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, ta2_xy_yy_zzzz_0, ta2_xy_yy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyy_xxxx_0[i] = 4.0 * ta2_xy_yy_xxx_0[i] * fe_0 - 4.0 * ta2_xy_yy_xxx_1[i] * fe_0 + ta1_y_yy_xxxx_1[i] + ta2_xy_yy_xxxx_0[i] * pa_x[i] - ta2_xy_yy_xxxx_1[i] * pc_x[i];

        ta2_xy_xyy_xxxy_0[i] = 3.0 * ta2_xy_yy_xxy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxy_1[i] * fe_0 + ta1_y_yy_xxxy_1[i] + ta2_xy_yy_xxxy_0[i] * pa_x[i] - ta2_xy_yy_xxxy_1[i] * pc_x[i];

        ta2_xy_xyy_xxxz_0[i] = 3.0 * ta2_xy_yy_xxz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxz_1[i] * fe_0 + ta1_y_yy_xxxz_1[i] + ta2_xy_yy_xxxz_0[i] * pa_x[i] - ta2_xy_yy_xxxz_1[i] * pc_x[i];

        ta2_xy_xyy_xxyy_0[i] = 2.0 * ta2_xy_yy_xyy_0[i] * fe_0 - 2.0 * ta2_xy_yy_xyy_1[i] * fe_0 + ta1_y_yy_xxyy_1[i] + ta2_xy_yy_xxyy_0[i] * pa_x[i] - ta2_xy_yy_xxyy_1[i] * pc_x[i];

        ta2_xy_xyy_xxyz_0[i] = 2.0 * ta2_xy_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xyz_1[i] * fe_0 + ta1_y_yy_xxyz_1[i] + ta2_xy_yy_xxyz_0[i] * pa_x[i] - ta2_xy_yy_xxyz_1[i] * pc_x[i];

        ta2_xy_xyy_xxzz_0[i] = 2.0 * ta2_xy_yy_xzz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xzz_1[i] * fe_0 + ta1_y_yy_xxzz_1[i] + ta2_xy_yy_xxzz_0[i] * pa_x[i] - ta2_xy_yy_xxzz_1[i] * pc_x[i];

        ta2_xy_xyy_xyyy_0[i] = ta2_xy_yy_yyy_0[i] * fe_0 - ta2_xy_yy_yyy_1[i] * fe_0 + ta1_y_yy_xyyy_1[i] + ta2_xy_yy_xyyy_0[i] * pa_x[i] - ta2_xy_yy_xyyy_1[i] * pc_x[i];

        ta2_xy_xyy_xyyz_0[i] = ta2_xy_yy_yyz_0[i] * fe_0 - ta2_xy_yy_yyz_1[i] * fe_0 + ta1_y_yy_xyyz_1[i] + ta2_xy_yy_xyyz_0[i] * pa_x[i] - ta2_xy_yy_xyyz_1[i] * pc_x[i];

        ta2_xy_xyy_xyzz_0[i] = ta2_xy_yy_yzz_0[i] * fe_0 - ta2_xy_yy_yzz_1[i] * fe_0 + ta1_y_yy_xyzz_1[i] + ta2_xy_yy_xyzz_0[i] * pa_x[i] - ta2_xy_yy_xyzz_1[i] * pc_x[i];

        ta2_xy_xyy_xzzz_0[i] = ta2_xy_yy_zzz_0[i] * fe_0 - ta2_xy_yy_zzz_1[i] * fe_0 + ta1_y_yy_xzzz_1[i] + ta2_xy_yy_xzzz_0[i] * pa_x[i] - ta2_xy_yy_xzzz_1[i] * pc_x[i];

        ta2_xy_xyy_yyyy_0[i] = ta1_y_yy_yyyy_1[i] + ta2_xy_yy_yyyy_0[i] * pa_x[i] - ta2_xy_yy_yyyy_1[i] * pc_x[i];

        ta2_xy_xyy_yyyz_0[i] = ta1_y_yy_yyyz_1[i] + ta2_xy_yy_yyyz_0[i] * pa_x[i] - ta2_xy_yy_yyyz_1[i] * pc_x[i];

        ta2_xy_xyy_yyzz_0[i] = ta1_y_yy_yyzz_1[i] + ta2_xy_yy_yyzz_0[i] * pa_x[i] - ta2_xy_yy_yyzz_1[i] * pc_x[i];

        ta2_xy_xyy_yzzz_0[i] = ta1_y_yy_yzzz_1[i] + ta2_xy_yy_yzzz_0[i] * pa_x[i] - ta2_xy_yy_yzzz_1[i] * pc_x[i];

        ta2_xy_xyy_zzzz_0[i] = ta1_y_yy_zzzz_1[i] + ta2_xy_yy_zzzz_0[i] * pa_x[i] - ta2_xy_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 210-225 components of targeted buffer : FG

    auto ta2_xy_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 210);

    auto ta2_xy_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 211);

    auto ta2_xy_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 212);

    auto ta2_xy_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 213);

    auto ta2_xy_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 214);

    auto ta2_xy_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 215);

    auto ta2_xy_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 216);

    auto ta2_xy_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 217);

    auto ta2_xy_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 218);

    auto ta2_xy_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 219);

    auto ta2_xy_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 220);

    auto ta2_xy_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 221);

    auto ta2_xy_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 222);

    auto ta2_xy_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 223);

    auto ta2_xy_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 224);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xz_xxxz_1, ta1_x_xz_xxzz_1, ta1_x_xz_xzzz_1, ta1_y_yz_yyyz_1, ta1_y_yz_yyzz_1, ta1_y_yz_yzzz_1, ta1_y_yz_zzzz_1, ta2_xy_xy_xxxx_0, ta2_xy_xy_xxxx_1, ta2_xy_xy_xxxy_0, ta2_xy_xy_xxxy_1, ta2_xy_xy_xxy_0, ta2_xy_xy_xxy_1, ta2_xy_xy_xxyy_0, ta2_xy_xy_xxyy_1, ta2_xy_xy_xxyz_0, ta2_xy_xy_xxyz_1, ta2_xy_xy_xyy_0, ta2_xy_xy_xyy_1, ta2_xy_xy_xyyy_0, ta2_xy_xy_xyyy_1, ta2_xy_xy_xyyz_0, ta2_xy_xy_xyyz_1, ta2_xy_xy_xyz_0, ta2_xy_xy_xyz_1, ta2_xy_xy_xyzz_0, ta2_xy_xy_xyzz_1, ta2_xy_xy_yyyy_0, ta2_xy_xy_yyyy_1, ta2_xy_xyz_xxxx_0, ta2_xy_xyz_xxxy_0, ta2_xy_xyz_xxxz_0, ta2_xy_xyz_xxyy_0, ta2_xy_xyz_xxyz_0, ta2_xy_xyz_xxzz_0, ta2_xy_xyz_xyyy_0, ta2_xy_xyz_xyyz_0, ta2_xy_xyz_xyzz_0, ta2_xy_xyz_xzzz_0, ta2_xy_xyz_yyyy_0, ta2_xy_xyz_yyyz_0, ta2_xy_xyz_yyzz_0, ta2_xy_xyz_yzzz_0, ta2_xy_xyz_zzzz_0, ta2_xy_xz_xxxz_0, ta2_xy_xz_xxxz_1, ta2_xy_xz_xxzz_0, ta2_xy_xz_xxzz_1, ta2_xy_xz_xzzz_0, ta2_xy_xz_xzzz_1, ta2_xy_yz_yyyz_0, ta2_xy_yz_yyyz_1, ta2_xy_yz_yyzz_0, ta2_xy_yz_yyzz_1, ta2_xy_yz_yzzz_0, ta2_xy_yz_yzzz_1, ta2_xy_yz_zzzz_0, ta2_xy_yz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyz_xxxx_0[i] = ta2_xy_xy_xxxx_0[i] * pa_z[i] - ta2_xy_xy_xxxx_1[i] * pc_z[i];

        ta2_xy_xyz_xxxy_0[i] = ta2_xy_xy_xxxy_0[i] * pa_z[i] - ta2_xy_xy_xxxy_1[i] * pc_z[i];

        ta2_xy_xyz_xxxz_0[i] = ta1_x_xz_xxxz_1[i] + ta2_xy_xz_xxxz_0[i] * pa_y[i] - ta2_xy_xz_xxxz_1[i] * pc_y[i];

        ta2_xy_xyz_xxyy_0[i] = ta2_xy_xy_xxyy_0[i] * pa_z[i] - ta2_xy_xy_xxyy_1[i] * pc_z[i];

        ta2_xy_xyz_xxyz_0[i] = ta2_xy_xy_xxy_0[i] * fe_0 - ta2_xy_xy_xxy_1[i] * fe_0 + ta2_xy_xy_xxyz_0[i] * pa_z[i] - ta2_xy_xy_xxyz_1[i] * pc_z[i];

        ta2_xy_xyz_xxzz_0[i] = ta1_x_xz_xxzz_1[i] + ta2_xy_xz_xxzz_0[i] * pa_y[i] - ta2_xy_xz_xxzz_1[i] * pc_y[i];

        ta2_xy_xyz_xyyy_0[i] = ta2_xy_xy_xyyy_0[i] * pa_z[i] - ta2_xy_xy_xyyy_1[i] * pc_z[i];

        ta2_xy_xyz_xyyz_0[i] = ta2_xy_xy_xyy_0[i] * fe_0 - ta2_xy_xy_xyy_1[i] * fe_0 + ta2_xy_xy_xyyz_0[i] * pa_z[i] - ta2_xy_xy_xyyz_1[i] * pc_z[i];

        ta2_xy_xyz_xyzz_0[i] = 2.0 * ta2_xy_xy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xy_xyz_1[i] * fe_0 + ta2_xy_xy_xyzz_0[i] * pa_z[i] - ta2_xy_xy_xyzz_1[i] * pc_z[i];

        ta2_xy_xyz_xzzz_0[i] = ta1_x_xz_xzzz_1[i] + ta2_xy_xz_xzzz_0[i] * pa_y[i] - ta2_xy_xz_xzzz_1[i] * pc_y[i];

        ta2_xy_xyz_yyyy_0[i] = ta2_xy_xy_yyyy_0[i] * pa_z[i] - ta2_xy_xy_yyyy_1[i] * pc_z[i];

        ta2_xy_xyz_yyyz_0[i] = ta1_y_yz_yyyz_1[i] + ta2_xy_yz_yyyz_0[i] * pa_x[i] - ta2_xy_yz_yyyz_1[i] * pc_x[i];

        ta2_xy_xyz_yyzz_0[i] = ta1_y_yz_yyzz_1[i] + ta2_xy_yz_yyzz_0[i] * pa_x[i] - ta2_xy_yz_yyzz_1[i] * pc_x[i];

        ta2_xy_xyz_yzzz_0[i] = ta1_y_yz_yzzz_1[i] + ta2_xy_yz_yzzz_0[i] * pa_x[i] - ta2_xy_yz_yzzz_1[i] * pc_x[i];

        ta2_xy_xyz_zzzz_0[i] = ta1_y_yz_zzzz_1[i] + ta2_xy_yz_zzzz_0[i] * pa_x[i] - ta2_xy_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : FG

    auto ta2_xy_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 225);

    auto ta2_xy_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 226);

    auto ta2_xy_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 227);

    auto ta2_xy_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 228);

    auto ta2_xy_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 229);

    auto ta2_xy_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 230);

    auto ta2_xy_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 231);

    auto ta2_xy_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 232);

    auto ta2_xy_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 233);

    auto ta2_xy_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 234);

    auto ta2_xy_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 235);

    auto ta2_xy_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 236);

    auto ta2_xy_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 237);

    auto ta2_xy_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 238);

    auto ta2_xy_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 239);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_zz_xxxz_1, ta1_y_zz_xxyz_1, ta1_y_zz_xxzz_1, ta1_y_zz_xyyz_1, ta1_y_zz_xyzz_1, ta1_y_zz_xzzz_1, ta1_y_zz_yyyy_1, ta1_y_zz_yyyz_1, ta1_y_zz_yyzz_1, ta1_y_zz_yzzz_1, ta1_y_zz_zzzz_1, ta2_xy_x_xxxx_0, ta2_xy_x_xxxx_1, ta2_xy_x_xxxy_0, ta2_xy_x_xxxy_1, ta2_xy_x_xxyy_0, ta2_xy_x_xxyy_1, ta2_xy_x_xyyy_0, ta2_xy_x_xyyy_1, ta2_xy_xz_xxxx_0, ta2_xy_xz_xxxx_1, ta2_xy_xz_xxxy_0, ta2_xy_xz_xxxy_1, ta2_xy_xz_xxyy_0, ta2_xy_xz_xxyy_1, ta2_xy_xz_xyyy_0, ta2_xy_xz_xyyy_1, ta2_xy_xzz_xxxx_0, ta2_xy_xzz_xxxy_0, ta2_xy_xzz_xxxz_0, ta2_xy_xzz_xxyy_0, ta2_xy_xzz_xxyz_0, ta2_xy_xzz_xxzz_0, ta2_xy_xzz_xyyy_0, ta2_xy_xzz_xyyz_0, ta2_xy_xzz_xyzz_0, ta2_xy_xzz_xzzz_0, ta2_xy_xzz_yyyy_0, ta2_xy_xzz_yyyz_0, ta2_xy_xzz_yyzz_0, ta2_xy_xzz_yzzz_0, ta2_xy_xzz_zzzz_0, ta2_xy_zz_xxxz_0, ta2_xy_zz_xxxz_1, ta2_xy_zz_xxyz_0, ta2_xy_zz_xxyz_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xxzz_0, ta2_xy_zz_xxzz_1, ta2_xy_zz_xyyz_0, ta2_xy_zz_xyyz_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xyzz_0, ta2_xy_zz_xyzz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_xzzz_0, ta2_xy_zz_xzzz_1, ta2_xy_zz_yyyy_0, ta2_xy_zz_yyyy_1, ta2_xy_zz_yyyz_0, ta2_xy_zz_yyyz_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yyzz_0, ta2_xy_zz_yyzz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_yzzz_0, ta2_xy_zz_yzzz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, ta2_xy_zz_zzzz_0, ta2_xy_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzz_xxxx_0[i] = ta2_xy_x_xxxx_0[i] * fe_0 - ta2_xy_x_xxxx_1[i] * fe_0 + ta2_xy_xz_xxxx_0[i] * pa_z[i] - ta2_xy_xz_xxxx_1[i] * pc_z[i];

        ta2_xy_xzz_xxxy_0[i] = ta2_xy_x_xxxy_0[i] * fe_0 - ta2_xy_x_xxxy_1[i] * fe_0 + ta2_xy_xz_xxxy_0[i] * pa_z[i] - ta2_xy_xz_xxxy_1[i] * pc_z[i];

        ta2_xy_xzz_xxxz_0[i] = 3.0 * ta2_xy_zz_xxz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxz_1[i] * fe_0 + ta1_y_zz_xxxz_1[i] + ta2_xy_zz_xxxz_0[i] * pa_x[i] - ta2_xy_zz_xxxz_1[i] * pc_x[i];

        ta2_xy_xzz_xxyy_0[i] = ta2_xy_x_xxyy_0[i] * fe_0 - ta2_xy_x_xxyy_1[i] * fe_0 + ta2_xy_xz_xxyy_0[i] * pa_z[i] - ta2_xy_xz_xxyy_1[i] * pc_z[i];

        ta2_xy_xzz_xxyz_0[i] = 2.0 * ta2_xy_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xyz_1[i] * fe_0 + ta1_y_zz_xxyz_1[i] + ta2_xy_zz_xxyz_0[i] * pa_x[i] - ta2_xy_zz_xxyz_1[i] * pc_x[i];

        ta2_xy_xzz_xxzz_0[i] = 2.0 * ta2_xy_zz_xzz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xzz_1[i] * fe_0 + ta1_y_zz_xxzz_1[i] + ta2_xy_zz_xxzz_0[i] * pa_x[i] - ta2_xy_zz_xxzz_1[i] * pc_x[i];

        ta2_xy_xzz_xyyy_0[i] = ta2_xy_x_xyyy_0[i] * fe_0 - ta2_xy_x_xyyy_1[i] * fe_0 + ta2_xy_xz_xyyy_0[i] * pa_z[i] - ta2_xy_xz_xyyy_1[i] * pc_z[i];

        ta2_xy_xzz_xyyz_0[i] = ta2_xy_zz_yyz_0[i] * fe_0 - ta2_xy_zz_yyz_1[i] * fe_0 + ta1_y_zz_xyyz_1[i] + ta2_xy_zz_xyyz_0[i] * pa_x[i] - ta2_xy_zz_xyyz_1[i] * pc_x[i];

        ta2_xy_xzz_xyzz_0[i] = ta2_xy_zz_yzz_0[i] * fe_0 - ta2_xy_zz_yzz_1[i] * fe_0 + ta1_y_zz_xyzz_1[i] + ta2_xy_zz_xyzz_0[i] * pa_x[i] - ta2_xy_zz_xyzz_1[i] * pc_x[i];

        ta2_xy_xzz_xzzz_0[i] = ta2_xy_zz_zzz_0[i] * fe_0 - ta2_xy_zz_zzz_1[i] * fe_0 + ta1_y_zz_xzzz_1[i] + ta2_xy_zz_xzzz_0[i] * pa_x[i] - ta2_xy_zz_xzzz_1[i] * pc_x[i];

        ta2_xy_xzz_yyyy_0[i] = ta1_y_zz_yyyy_1[i] + ta2_xy_zz_yyyy_0[i] * pa_x[i] - ta2_xy_zz_yyyy_1[i] * pc_x[i];

        ta2_xy_xzz_yyyz_0[i] = ta1_y_zz_yyyz_1[i] + ta2_xy_zz_yyyz_0[i] * pa_x[i] - ta2_xy_zz_yyyz_1[i] * pc_x[i];

        ta2_xy_xzz_yyzz_0[i] = ta1_y_zz_yyzz_1[i] + ta2_xy_zz_yyzz_0[i] * pa_x[i] - ta2_xy_zz_yyzz_1[i] * pc_x[i];

        ta2_xy_xzz_yzzz_0[i] = ta1_y_zz_yzzz_1[i] + ta2_xy_zz_yzzz_0[i] * pa_x[i] - ta2_xy_zz_yzzz_1[i] * pc_x[i];

        ta2_xy_xzz_zzzz_0[i] = ta1_y_zz_zzzz_1[i] + ta2_xy_zz_zzzz_0[i] * pa_x[i] - ta2_xy_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : FG

    auto ta2_xy_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 240);

    auto ta2_xy_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 241);

    auto ta2_xy_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 242);

    auto ta2_xy_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 243);

    auto ta2_xy_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 244);

    auto ta2_xy_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 245);

    auto ta2_xy_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 246);

    auto ta2_xy_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 247);

    auto ta2_xy_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 248);

    auto ta2_xy_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 249);

    auto ta2_xy_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 250);

    auto ta2_xy_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 251);

    auto ta2_xy_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 252);

    auto ta2_xy_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 253);

    auto ta2_xy_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 254);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yy_xxxx_1, ta1_x_yy_xxxy_1, ta1_x_yy_xxxz_1, ta1_x_yy_xxyy_1, ta1_x_yy_xxyz_1, ta1_x_yy_xxzz_1, ta1_x_yy_xyyy_1, ta1_x_yy_xyyz_1, ta1_x_yy_xyzz_1, ta1_x_yy_xzzz_1, ta1_x_yy_yyyy_1, ta1_x_yy_yyyz_1, ta1_x_yy_yyzz_1, ta1_x_yy_yzzz_1, ta1_x_yy_zzzz_1, ta2_xy_y_xxxx_0, ta2_xy_y_xxxx_1, ta2_xy_y_xxxy_0, ta2_xy_y_xxxy_1, ta2_xy_y_xxxz_0, ta2_xy_y_xxxz_1, ta2_xy_y_xxyy_0, ta2_xy_y_xxyy_1, ta2_xy_y_xxyz_0, ta2_xy_y_xxyz_1, ta2_xy_y_xxzz_0, ta2_xy_y_xxzz_1, ta2_xy_y_xyyy_0, ta2_xy_y_xyyy_1, ta2_xy_y_xyyz_0, ta2_xy_y_xyyz_1, ta2_xy_y_xyzz_0, ta2_xy_y_xyzz_1, ta2_xy_y_xzzz_0, ta2_xy_y_xzzz_1, ta2_xy_y_yyyy_0, ta2_xy_y_yyyy_1, ta2_xy_y_yyyz_0, ta2_xy_y_yyyz_1, ta2_xy_y_yyzz_0, ta2_xy_y_yyzz_1, ta2_xy_y_yzzz_0, ta2_xy_y_yzzz_1, ta2_xy_y_zzzz_0, ta2_xy_y_zzzz_1, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxxx_0, ta2_xy_yy_xxxx_1, ta2_xy_yy_xxxy_0, ta2_xy_yy_xxxy_1, ta2_xy_yy_xxxz_0, ta2_xy_yy_xxxz_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxyy_0, ta2_xy_yy_xxyy_1, ta2_xy_yy_xxyz_0, ta2_xy_yy_xxyz_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xxzz_0, ta2_xy_yy_xxzz_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyyy_0, ta2_xy_yy_xyyy_1, ta2_xy_yy_xyyz_0, ta2_xy_yy_xyyz_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xyzz_0, ta2_xy_yy_xyzz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_xzzz_0, ta2_xy_yy_xzzz_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyyy_0, ta2_xy_yy_yyyy_1, ta2_xy_yy_yyyz_0, ta2_xy_yy_yyyz_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yyzz_0, ta2_xy_yy_yyzz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_yzzz_0, ta2_xy_yy_yzzz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, ta2_xy_yy_zzzz_0, ta2_xy_yy_zzzz_1, ta2_xy_yyy_xxxx_0, ta2_xy_yyy_xxxy_0, ta2_xy_yyy_xxxz_0, ta2_xy_yyy_xxyy_0, ta2_xy_yyy_xxyz_0, ta2_xy_yyy_xxzz_0, ta2_xy_yyy_xyyy_0, ta2_xy_yyy_xyyz_0, ta2_xy_yyy_xyzz_0, ta2_xy_yyy_xzzz_0, ta2_xy_yyy_yyyy_0, ta2_xy_yyy_yyyz_0, ta2_xy_yyy_yyzz_0, ta2_xy_yyy_yzzz_0, ta2_xy_yyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyy_xxxx_0[i] = 2.0 * ta2_xy_y_xxxx_0[i] * fe_0 - 2.0 * ta2_xy_y_xxxx_1[i] * fe_0 + ta1_x_yy_xxxx_1[i] + ta2_xy_yy_xxxx_0[i] * pa_y[i] - ta2_xy_yy_xxxx_1[i] * pc_y[i];

        ta2_xy_yyy_xxxy_0[i] = 2.0 * ta2_xy_y_xxxy_0[i] * fe_0 - 2.0 * ta2_xy_y_xxxy_1[i] * fe_0 + ta2_xy_yy_xxx_0[i] * fe_0 - ta2_xy_yy_xxx_1[i] * fe_0 + ta1_x_yy_xxxy_1[i] + ta2_xy_yy_xxxy_0[i] * pa_y[i] - ta2_xy_yy_xxxy_1[i] * pc_y[i];

        ta2_xy_yyy_xxxz_0[i] = 2.0 * ta2_xy_y_xxxz_0[i] * fe_0 - 2.0 * ta2_xy_y_xxxz_1[i] * fe_0 + ta1_x_yy_xxxz_1[i] + ta2_xy_yy_xxxz_0[i] * pa_y[i] - ta2_xy_yy_xxxz_1[i] * pc_y[i];

        ta2_xy_yyy_xxyy_0[i] = 2.0 * ta2_xy_y_xxyy_0[i] * fe_0 - 2.0 * ta2_xy_y_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_yy_xxy_0[i] * fe_0 - 2.0 * ta2_xy_yy_xxy_1[i] * fe_0 + ta1_x_yy_xxyy_1[i] + ta2_xy_yy_xxyy_0[i] * pa_y[i] - ta2_xy_yy_xxyy_1[i] * pc_y[i];

        ta2_xy_yyy_xxyz_0[i] = 2.0 * ta2_xy_y_xxyz_0[i] * fe_0 - 2.0 * ta2_xy_y_xxyz_1[i] * fe_0 + ta2_xy_yy_xxz_0[i] * fe_0 - ta2_xy_yy_xxz_1[i] * fe_0 + ta1_x_yy_xxyz_1[i] + ta2_xy_yy_xxyz_0[i] * pa_y[i] - ta2_xy_yy_xxyz_1[i] * pc_y[i];

        ta2_xy_yyy_xxzz_0[i] = 2.0 * ta2_xy_y_xxzz_0[i] * fe_0 - 2.0 * ta2_xy_y_xxzz_1[i] * fe_0 + ta1_x_yy_xxzz_1[i] + ta2_xy_yy_xxzz_0[i] * pa_y[i] - ta2_xy_yy_xxzz_1[i] * pc_y[i];

        ta2_xy_yyy_xyyy_0[i] = 2.0 * ta2_xy_y_xyyy_0[i] * fe_0 - 2.0 * ta2_xy_y_xyyy_1[i] * fe_0 + 3.0 * ta2_xy_yy_xyy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xyy_1[i] * fe_0 + ta1_x_yy_xyyy_1[i] + ta2_xy_yy_xyyy_0[i] * pa_y[i] - ta2_xy_yy_xyyy_1[i] * pc_y[i];

        ta2_xy_yyy_xyyz_0[i] = 2.0 * ta2_xy_y_xyyz_0[i] * fe_0 - 2.0 * ta2_xy_y_xyyz_1[i] * fe_0 + 2.0 * ta2_xy_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xyz_1[i] * fe_0 + ta1_x_yy_xyyz_1[i] + ta2_xy_yy_xyyz_0[i] * pa_y[i] - ta2_xy_yy_xyyz_1[i] * pc_y[i];

        ta2_xy_yyy_xyzz_0[i] = 2.0 * ta2_xy_y_xyzz_0[i] * fe_0 - 2.0 * ta2_xy_y_xyzz_1[i] * fe_0 + ta2_xy_yy_xzz_0[i] * fe_0 - ta2_xy_yy_xzz_1[i] * fe_0 + ta1_x_yy_xyzz_1[i] + ta2_xy_yy_xyzz_0[i] * pa_y[i] - ta2_xy_yy_xyzz_1[i] * pc_y[i];

        ta2_xy_yyy_xzzz_0[i] = 2.0 * ta2_xy_y_xzzz_0[i] * fe_0 - 2.0 * ta2_xy_y_xzzz_1[i] * fe_0 + ta1_x_yy_xzzz_1[i] + ta2_xy_yy_xzzz_0[i] * pa_y[i] - ta2_xy_yy_xzzz_1[i] * pc_y[i];

        ta2_xy_yyy_yyyy_0[i] = 2.0 * ta2_xy_y_yyyy_0[i] * fe_0 - 2.0 * ta2_xy_y_yyyy_1[i] * fe_0 + 4.0 * ta2_xy_yy_yyy_0[i] * fe_0 - 4.0 * ta2_xy_yy_yyy_1[i] * fe_0 + ta1_x_yy_yyyy_1[i] + ta2_xy_yy_yyyy_0[i] * pa_y[i] - ta2_xy_yy_yyyy_1[i] * pc_y[i];

        ta2_xy_yyy_yyyz_0[i] = 2.0 * ta2_xy_y_yyyz_0[i] * fe_0 - 2.0 * ta2_xy_y_yyyz_1[i] * fe_0 + 3.0 * ta2_xy_yy_yyz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yyz_1[i] * fe_0 + ta1_x_yy_yyyz_1[i] + ta2_xy_yy_yyyz_0[i] * pa_y[i] - ta2_xy_yy_yyyz_1[i] * pc_y[i];

        ta2_xy_yyy_yyzz_0[i] = 2.0 * ta2_xy_y_yyzz_0[i] * fe_0 - 2.0 * ta2_xy_y_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_yy_yzz_0[i] * fe_0 - 2.0 * ta2_xy_yy_yzz_1[i] * fe_0 + ta1_x_yy_yyzz_1[i] + ta2_xy_yy_yyzz_0[i] * pa_y[i] - ta2_xy_yy_yyzz_1[i] * pc_y[i];

        ta2_xy_yyy_yzzz_0[i] = 2.0 * ta2_xy_y_yzzz_0[i] * fe_0 - 2.0 * ta2_xy_y_yzzz_1[i] * fe_0 + ta2_xy_yy_zzz_0[i] * fe_0 - ta2_xy_yy_zzz_1[i] * fe_0 + ta1_x_yy_yzzz_1[i] + ta2_xy_yy_yzzz_0[i] * pa_y[i] - ta2_xy_yy_yzzz_1[i] * pc_y[i];

        ta2_xy_yyy_zzzz_0[i] = 2.0 * ta2_xy_y_zzzz_0[i] * fe_0 - 2.0 * ta2_xy_y_zzzz_1[i] * fe_0 + ta1_x_yy_zzzz_1[i] + ta2_xy_yy_zzzz_0[i] * pa_y[i] - ta2_xy_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : FG

    auto ta2_xy_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 255);

    auto ta2_xy_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 256);

    auto ta2_xy_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 257);

    auto ta2_xy_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 258);

    auto ta2_xy_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 259);

    auto ta2_xy_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 260);

    auto ta2_xy_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 261);

    auto ta2_xy_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 262);

    auto ta2_xy_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 263);

    auto ta2_xy_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 264);

    auto ta2_xy_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 265);

    auto ta2_xy_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 266);

    auto ta2_xy_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 267);

    auto ta2_xy_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 268);

    auto ta2_xy_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 269);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxxx_0, ta2_xy_yy_xxxx_1, ta2_xy_yy_xxxy_0, ta2_xy_yy_xxxy_1, ta2_xy_yy_xxxz_0, ta2_xy_yy_xxxz_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxyy_0, ta2_xy_yy_xxyy_1, ta2_xy_yy_xxyz_0, ta2_xy_yy_xxyz_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xxzz_0, ta2_xy_yy_xxzz_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyyy_0, ta2_xy_yy_xyyy_1, ta2_xy_yy_xyyz_0, ta2_xy_yy_xyyz_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xyzz_0, ta2_xy_yy_xyzz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_xzzz_0, ta2_xy_yy_xzzz_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyyy_0, ta2_xy_yy_yyyy_1, ta2_xy_yy_yyyz_0, ta2_xy_yy_yyyz_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yyzz_0, ta2_xy_yy_yyzz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_yzzz_0, ta2_xy_yy_yzzz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, ta2_xy_yy_zzzz_0, ta2_xy_yy_zzzz_1, ta2_xy_yyz_xxxx_0, ta2_xy_yyz_xxxy_0, ta2_xy_yyz_xxxz_0, ta2_xy_yyz_xxyy_0, ta2_xy_yyz_xxyz_0, ta2_xy_yyz_xxzz_0, ta2_xy_yyz_xyyy_0, ta2_xy_yyz_xyyz_0, ta2_xy_yyz_xyzz_0, ta2_xy_yyz_xzzz_0, ta2_xy_yyz_yyyy_0, ta2_xy_yyz_yyyz_0, ta2_xy_yyz_yyzz_0, ta2_xy_yyz_yzzz_0, ta2_xy_yyz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyz_xxxx_0[i] = ta2_xy_yy_xxxx_0[i] * pa_z[i] - ta2_xy_yy_xxxx_1[i] * pc_z[i];

        ta2_xy_yyz_xxxy_0[i] = ta2_xy_yy_xxxy_0[i] * pa_z[i] - ta2_xy_yy_xxxy_1[i] * pc_z[i];

        ta2_xy_yyz_xxxz_0[i] = ta2_xy_yy_xxx_0[i] * fe_0 - ta2_xy_yy_xxx_1[i] * fe_0 + ta2_xy_yy_xxxz_0[i] * pa_z[i] - ta2_xy_yy_xxxz_1[i] * pc_z[i];

        ta2_xy_yyz_xxyy_0[i] = ta2_xy_yy_xxyy_0[i] * pa_z[i] - ta2_xy_yy_xxyy_1[i] * pc_z[i];

        ta2_xy_yyz_xxyz_0[i] = ta2_xy_yy_xxy_0[i] * fe_0 - ta2_xy_yy_xxy_1[i] * fe_0 + ta2_xy_yy_xxyz_0[i] * pa_z[i] - ta2_xy_yy_xxyz_1[i] * pc_z[i];

        ta2_xy_yyz_xxzz_0[i] = 2.0 * ta2_xy_yy_xxz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xxz_1[i] * fe_0 + ta2_xy_yy_xxzz_0[i] * pa_z[i] - ta2_xy_yy_xxzz_1[i] * pc_z[i];

        ta2_xy_yyz_xyyy_0[i] = ta2_xy_yy_xyyy_0[i] * pa_z[i] - ta2_xy_yy_xyyy_1[i] * pc_z[i];

        ta2_xy_yyz_xyyz_0[i] = ta2_xy_yy_xyy_0[i] * fe_0 - ta2_xy_yy_xyy_1[i] * fe_0 + ta2_xy_yy_xyyz_0[i] * pa_z[i] - ta2_xy_yy_xyyz_1[i] * pc_z[i];

        ta2_xy_yyz_xyzz_0[i] = 2.0 * ta2_xy_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xyz_1[i] * fe_0 + ta2_xy_yy_xyzz_0[i] * pa_z[i] - ta2_xy_yy_xyzz_1[i] * pc_z[i];

        ta2_xy_yyz_xzzz_0[i] = 3.0 * ta2_xy_yy_xzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xzz_1[i] * fe_0 + ta2_xy_yy_xzzz_0[i] * pa_z[i] - ta2_xy_yy_xzzz_1[i] * pc_z[i];

        ta2_xy_yyz_yyyy_0[i] = ta2_xy_yy_yyyy_0[i] * pa_z[i] - ta2_xy_yy_yyyy_1[i] * pc_z[i];

        ta2_xy_yyz_yyyz_0[i] = ta2_xy_yy_yyy_0[i] * fe_0 - ta2_xy_yy_yyy_1[i] * fe_0 + ta2_xy_yy_yyyz_0[i] * pa_z[i] - ta2_xy_yy_yyyz_1[i] * pc_z[i];

        ta2_xy_yyz_yyzz_0[i] = 2.0 * ta2_xy_yy_yyz_0[i] * fe_0 - 2.0 * ta2_xy_yy_yyz_1[i] * fe_0 + ta2_xy_yy_yyzz_0[i] * pa_z[i] - ta2_xy_yy_yyzz_1[i] * pc_z[i];

        ta2_xy_yyz_yzzz_0[i] = 3.0 * ta2_xy_yy_yzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yzz_1[i] * fe_0 + ta2_xy_yy_yzzz_0[i] * pa_z[i] - ta2_xy_yy_yzzz_1[i] * pc_z[i];

        ta2_xy_yyz_zzzz_0[i] = 4.0 * ta2_xy_yy_zzz_0[i] * fe_0 - 4.0 * ta2_xy_yy_zzz_1[i] * fe_0 + ta2_xy_yy_zzzz_0[i] * pa_z[i] - ta2_xy_yy_zzzz_1[i] * pc_z[i];
    }

    // Set up 270-285 components of targeted buffer : FG

    auto ta2_xy_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 270);

    auto ta2_xy_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 271);

    auto ta2_xy_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 272);

    auto ta2_xy_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 273);

    auto ta2_xy_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 274);

    auto ta2_xy_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 275);

    auto ta2_xy_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 276);

    auto ta2_xy_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 277);

    auto ta2_xy_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 278);

    auto ta2_xy_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 279);

    auto ta2_xy_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 280);

    auto ta2_xy_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 281);

    auto ta2_xy_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 282);

    auto ta2_xy_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 283);

    auto ta2_xy_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 284);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_zz_xxxx_1, ta1_x_zz_xxxz_1, ta1_x_zz_xxyz_1, ta1_x_zz_xxzz_1, ta1_x_zz_xyyz_1, ta1_x_zz_xyzz_1, ta1_x_zz_xzzz_1, ta1_x_zz_yyyz_1, ta1_x_zz_yyzz_1, ta1_x_zz_yzzz_1, ta1_x_zz_zzzz_1, ta2_xy_y_xxxy_0, ta2_xy_y_xxxy_1, ta2_xy_y_xxyy_0, ta2_xy_y_xxyy_1, ta2_xy_y_xyyy_0, ta2_xy_y_xyyy_1, ta2_xy_y_yyyy_0, ta2_xy_y_yyyy_1, ta2_xy_yz_xxxy_0, ta2_xy_yz_xxxy_1, ta2_xy_yz_xxyy_0, ta2_xy_yz_xxyy_1, ta2_xy_yz_xyyy_0, ta2_xy_yz_xyyy_1, ta2_xy_yz_yyyy_0, ta2_xy_yz_yyyy_1, ta2_xy_yzz_xxxx_0, ta2_xy_yzz_xxxy_0, ta2_xy_yzz_xxxz_0, ta2_xy_yzz_xxyy_0, ta2_xy_yzz_xxyz_0, ta2_xy_yzz_xxzz_0, ta2_xy_yzz_xyyy_0, ta2_xy_yzz_xyyz_0, ta2_xy_yzz_xyzz_0, ta2_xy_yzz_xzzz_0, ta2_xy_yzz_yyyy_0, ta2_xy_yzz_yyyz_0, ta2_xy_yzz_yyzz_0, ta2_xy_yzz_yzzz_0, ta2_xy_yzz_zzzz_0, ta2_xy_zz_xxxx_0, ta2_xy_zz_xxxx_1, ta2_xy_zz_xxxz_0, ta2_xy_zz_xxxz_1, ta2_xy_zz_xxyz_0, ta2_xy_zz_xxyz_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xxzz_0, ta2_xy_zz_xxzz_1, ta2_xy_zz_xyyz_0, ta2_xy_zz_xyyz_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xyzz_0, ta2_xy_zz_xyzz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_xzzz_0, ta2_xy_zz_xzzz_1, ta2_xy_zz_yyyz_0, ta2_xy_zz_yyyz_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yyzz_0, ta2_xy_zz_yyzz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_yzzz_0, ta2_xy_zz_yzzz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, ta2_xy_zz_zzzz_0, ta2_xy_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzz_xxxx_0[i] = ta1_x_zz_xxxx_1[i] + ta2_xy_zz_xxxx_0[i] * pa_y[i] - ta2_xy_zz_xxxx_1[i] * pc_y[i];

        ta2_xy_yzz_xxxy_0[i] = ta2_xy_y_xxxy_0[i] * fe_0 - ta2_xy_y_xxxy_1[i] * fe_0 + ta2_xy_yz_xxxy_0[i] * pa_z[i] - ta2_xy_yz_xxxy_1[i] * pc_z[i];

        ta2_xy_yzz_xxxz_0[i] = ta1_x_zz_xxxz_1[i] + ta2_xy_zz_xxxz_0[i] * pa_y[i] - ta2_xy_zz_xxxz_1[i] * pc_y[i];

        ta2_xy_yzz_xxyy_0[i] = ta2_xy_y_xxyy_0[i] * fe_0 - ta2_xy_y_xxyy_1[i] * fe_0 + ta2_xy_yz_xxyy_0[i] * pa_z[i] - ta2_xy_yz_xxyy_1[i] * pc_z[i];

        ta2_xy_yzz_xxyz_0[i] = ta2_xy_zz_xxz_0[i] * fe_0 - ta2_xy_zz_xxz_1[i] * fe_0 + ta1_x_zz_xxyz_1[i] + ta2_xy_zz_xxyz_0[i] * pa_y[i] - ta2_xy_zz_xxyz_1[i] * pc_y[i];

        ta2_xy_yzz_xxzz_0[i] = ta1_x_zz_xxzz_1[i] + ta2_xy_zz_xxzz_0[i] * pa_y[i] - ta2_xy_zz_xxzz_1[i] * pc_y[i];

        ta2_xy_yzz_xyyy_0[i] = ta2_xy_y_xyyy_0[i] * fe_0 - ta2_xy_y_xyyy_1[i] * fe_0 + ta2_xy_yz_xyyy_0[i] * pa_z[i] - ta2_xy_yz_xyyy_1[i] * pc_z[i];

        ta2_xy_yzz_xyyz_0[i] = 2.0 * ta2_xy_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xyz_1[i] * fe_0 + ta1_x_zz_xyyz_1[i] + ta2_xy_zz_xyyz_0[i] * pa_y[i] - ta2_xy_zz_xyyz_1[i] * pc_y[i];

        ta2_xy_yzz_xyzz_0[i] = ta2_xy_zz_xzz_0[i] * fe_0 - ta2_xy_zz_xzz_1[i] * fe_0 + ta1_x_zz_xyzz_1[i] + ta2_xy_zz_xyzz_0[i] * pa_y[i] - ta2_xy_zz_xyzz_1[i] * pc_y[i];

        ta2_xy_yzz_xzzz_0[i] = ta1_x_zz_xzzz_1[i] + ta2_xy_zz_xzzz_0[i] * pa_y[i] - ta2_xy_zz_xzzz_1[i] * pc_y[i];

        ta2_xy_yzz_yyyy_0[i] = ta2_xy_y_yyyy_0[i] * fe_0 - ta2_xy_y_yyyy_1[i] * fe_0 + ta2_xy_yz_yyyy_0[i] * pa_z[i] - ta2_xy_yz_yyyy_1[i] * pc_z[i];

        ta2_xy_yzz_yyyz_0[i] = 3.0 * ta2_xy_zz_yyz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yyz_1[i] * fe_0 + ta1_x_zz_yyyz_1[i] + ta2_xy_zz_yyyz_0[i] * pa_y[i] - ta2_xy_zz_yyyz_1[i] * pc_y[i];

        ta2_xy_yzz_yyzz_0[i] = 2.0 * ta2_xy_zz_yzz_0[i] * fe_0 - 2.0 * ta2_xy_zz_yzz_1[i] * fe_0 + ta1_x_zz_yyzz_1[i] + ta2_xy_zz_yyzz_0[i] * pa_y[i] - ta2_xy_zz_yyzz_1[i] * pc_y[i];

        ta2_xy_yzz_yzzz_0[i] = ta2_xy_zz_zzz_0[i] * fe_0 - ta2_xy_zz_zzz_1[i] * fe_0 + ta1_x_zz_yzzz_1[i] + ta2_xy_zz_yzzz_0[i] * pa_y[i] - ta2_xy_zz_yzzz_1[i] * pc_y[i];

        ta2_xy_yzz_zzzz_0[i] = ta1_x_zz_zzzz_1[i] + ta2_xy_zz_zzzz_0[i] * pa_y[i] - ta2_xy_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 285-300 components of targeted buffer : FG

    auto ta2_xy_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 285);

    auto ta2_xy_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 286);

    auto ta2_xy_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 287);

    auto ta2_xy_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 288);

    auto ta2_xy_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 289);

    auto ta2_xy_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 290);

    auto ta2_xy_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 291);

    auto ta2_xy_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 292);

    auto ta2_xy_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 293);

    auto ta2_xy_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 294);

    auto ta2_xy_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 295);

    auto ta2_xy_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 296);

    auto ta2_xy_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 297);

    auto ta2_xy_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 298);

    auto ta2_xy_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 299);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_z_xxxx_0, ta2_xy_z_xxxx_1, ta2_xy_z_xxxy_0, ta2_xy_z_xxxy_1, ta2_xy_z_xxxz_0, ta2_xy_z_xxxz_1, ta2_xy_z_xxyy_0, ta2_xy_z_xxyy_1, ta2_xy_z_xxyz_0, ta2_xy_z_xxyz_1, ta2_xy_z_xxzz_0, ta2_xy_z_xxzz_1, ta2_xy_z_xyyy_0, ta2_xy_z_xyyy_1, ta2_xy_z_xyyz_0, ta2_xy_z_xyyz_1, ta2_xy_z_xyzz_0, ta2_xy_z_xyzz_1, ta2_xy_z_xzzz_0, ta2_xy_z_xzzz_1, ta2_xy_z_yyyy_0, ta2_xy_z_yyyy_1, ta2_xy_z_yyyz_0, ta2_xy_z_yyyz_1, ta2_xy_z_yyzz_0, ta2_xy_z_yyzz_1, ta2_xy_z_yzzz_0, ta2_xy_z_yzzz_1, ta2_xy_z_zzzz_0, ta2_xy_z_zzzz_1, ta2_xy_zz_xxx_0, ta2_xy_zz_xxx_1, ta2_xy_zz_xxxx_0, ta2_xy_zz_xxxx_1, ta2_xy_zz_xxxy_0, ta2_xy_zz_xxxy_1, ta2_xy_zz_xxxz_0, ta2_xy_zz_xxxz_1, ta2_xy_zz_xxy_0, ta2_xy_zz_xxy_1, ta2_xy_zz_xxyy_0, ta2_xy_zz_xxyy_1, ta2_xy_zz_xxyz_0, ta2_xy_zz_xxyz_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xxzz_0, ta2_xy_zz_xxzz_1, ta2_xy_zz_xyy_0, ta2_xy_zz_xyy_1, ta2_xy_zz_xyyy_0, ta2_xy_zz_xyyy_1, ta2_xy_zz_xyyz_0, ta2_xy_zz_xyyz_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xyzz_0, ta2_xy_zz_xyzz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_xzzz_0, ta2_xy_zz_xzzz_1, ta2_xy_zz_yyy_0, ta2_xy_zz_yyy_1, ta2_xy_zz_yyyy_0, ta2_xy_zz_yyyy_1, ta2_xy_zz_yyyz_0, ta2_xy_zz_yyyz_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yyzz_0, ta2_xy_zz_yyzz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_yzzz_0, ta2_xy_zz_yzzz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, ta2_xy_zz_zzzz_0, ta2_xy_zz_zzzz_1, ta2_xy_zzz_xxxx_0, ta2_xy_zzz_xxxy_0, ta2_xy_zzz_xxxz_0, ta2_xy_zzz_xxyy_0, ta2_xy_zzz_xxyz_0, ta2_xy_zzz_xxzz_0, ta2_xy_zzz_xyyy_0, ta2_xy_zzz_xyyz_0, ta2_xy_zzz_xyzz_0, ta2_xy_zzz_xzzz_0, ta2_xy_zzz_yyyy_0, ta2_xy_zzz_yyyz_0, ta2_xy_zzz_yyzz_0, ta2_xy_zzz_yzzz_0, ta2_xy_zzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzz_xxxx_0[i] = 2.0 * ta2_xy_z_xxxx_0[i] * fe_0 - 2.0 * ta2_xy_z_xxxx_1[i] * fe_0 + ta2_xy_zz_xxxx_0[i] * pa_z[i] - ta2_xy_zz_xxxx_1[i] * pc_z[i];

        ta2_xy_zzz_xxxy_0[i] = 2.0 * ta2_xy_z_xxxy_0[i] * fe_0 - 2.0 * ta2_xy_z_xxxy_1[i] * fe_0 + ta2_xy_zz_xxxy_0[i] * pa_z[i] - ta2_xy_zz_xxxy_1[i] * pc_z[i];

        ta2_xy_zzz_xxxz_0[i] = 2.0 * ta2_xy_z_xxxz_0[i] * fe_0 - 2.0 * ta2_xy_z_xxxz_1[i] * fe_0 + ta2_xy_zz_xxx_0[i] * fe_0 - ta2_xy_zz_xxx_1[i] * fe_0 + ta2_xy_zz_xxxz_0[i] * pa_z[i] - ta2_xy_zz_xxxz_1[i] * pc_z[i];

        ta2_xy_zzz_xxyy_0[i] = 2.0 * ta2_xy_z_xxyy_0[i] * fe_0 - 2.0 * ta2_xy_z_xxyy_1[i] * fe_0 + ta2_xy_zz_xxyy_0[i] * pa_z[i] - ta2_xy_zz_xxyy_1[i] * pc_z[i];

        ta2_xy_zzz_xxyz_0[i] = 2.0 * ta2_xy_z_xxyz_0[i] * fe_0 - 2.0 * ta2_xy_z_xxyz_1[i] * fe_0 + ta2_xy_zz_xxy_0[i] * fe_0 - ta2_xy_zz_xxy_1[i] * fe_0 + ta2_xy_zz_xxyz_0[i] * pa_z[i] - ta2_xy_zz_xxyz_1[i] * pc_z[i];

        ta2_xy_zzz_xxzz_0[i] = 2.0 * ta2_xy_z_xxzz_0[i] * fe_0 - 2.0 * ta2_xy_z_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_zz_xxz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xxz_1[i] * fe_0 + ta2_xy_zz_xxzz_0[i] * pa_z[i] - ta2_xy_zz_xxzz_1[i] * pc_z[i];

        ta2_xy_zzz_xyyy_0[i] = 2.0 * ta2_xy_z_xyyy_0[i] * fe_0 - 2.0 * ta2_xy_z_xyyy_1[i] * fe_0 + ta2_xy_zz_xyyy_0[i] * pa_z[i] - ta2_xy_zz_xyyy_1[i] * pc_z[i];

        ta2_xy_zzz_xyyz_0[i] = 2.0 * ta2_xy_z_xyyz_0[i] * fe_0 - 2.0 * ta2_xy_z_xyyz_1[i] * fe_0 + ta2_xy_zz_xyy_0[i] * fe_0 - ta2_xy_zz_xyy_1[i] * fe_0 + ta2_xy_zz_xyyz_0[i] * pa_z[i] - ta2_xy_zz_xyyz_1[i] * pc_z[i];

        ta2_xy_zzz_xyzz_0[i] = 2.0 * ta2_xy_z_xyzz_0[i] * fe_0 - 2.0 * ta2_xy_z_xyzz_1[i] * fe_0 + 2.0 * ta2_xy_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xyz_1[i] * fe_0 + ta2_xy_zz_xyzz_0[i] * pa_z[i] - ta2_xy_zz_xyzz_1[i] * pc_z[i];

        ta2_xy_zzz_xzzz_0[i] = 2.0 * ta2_xy_z_xzzz_0[i] * fe_0 - 2.0 * ta2_xy_z_xzzz_1[i] * fe_0 + 3.0 * ta2_xy_zz_xzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xzz_1[i] * fe_0 + ta2_xy_zz_xzzz_0[i] * pa_z[i] - ta2_xy_zz_xzzz_1[i] * pc_z[i];

        ta2_xy_zzz_yyyy_0[i] = 2.0 * ta2_xy_z_yyyy_0[i] * fe_0 - 2.0 * ta2_xy_z_yyyy_1[i] * fe_0 + ta2_xy_zz_yyyy_0[i] * pa_z[i] - ta2_xy_zz_yyyy_1[i] * pc_z[i];

        ta2_xy_zzz_yyyz_0[i] = 2.0 * ta2_xy_z_yyyz_0[i] * fe_0 - 2.0 * ta2_xy_z_yyyz_1[i] * fe_0 + ta2_xy_zz_yyy_0[i] * fe_0 - ta2_xy_zz_yyy_1[i] * fe_0 + ta2_xy_zz_yyyz_0[i] * pa_z[i] - ta2_xy_zz_yyyz_1[i] * pc_z[i];

        ta2_xy_zzz_yyzz_0[i] = 2.0 * ta2_xy_z_yyzz_0[i] * fe_0 - 2.0 * ta2_xy_z_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_zz_yyz_0[i] * fe_0 - 2.0 * ta2_xy_zz_yyz_1[i] * fe_0 + ta2_xy_zz_yyzz_0[i] * pa_z[i] - ta2_xy_zz_yyzz_1[i] * pc_z[i];

        ta2_xy_zzz_yzzz_0[i] = 2.0 * ta2_xy_z_yzzz_0[i] * fe_0 - 2.0 * ta2_xy_z_yzzz_1[i] * fe_0 + 3.0 * ta2_xy_zz_yzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yzz_1[i] * fe_0 + ta2_xy_zz_yzzz_0[i] * pa_z[i] - ta2_xy_zz_yzzz_1[i] * pc_z[i];

        ta2_xy_zzz_zzzz_0[i] = 2.0 * ta2_xy_z_zzzz_0[i] * fe_0 - 2.0 * ta2_xy_z_zzzz_1[i] * fe_0 + 4.0 * ta2_xy_zz_zzz_0[i] * fe_0 - 4.0 * ta2_xy_zz_zzz_1[i] * fe_0 + ta2_xy_zz_zzzz_0[i] * pa_z[i] - ta2_xy_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 300-315 components of targeted buffer : FG

    auto ta2_xz_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 300);

    auto ta2_xz_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 301);

    auto ta2_xz_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 302);

    auto ta2_xz_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 303);

    auto ta2_xz_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 304);

    auto ta2_xz_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 305);

    auto ta2_xz_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 306);

    auto ta2_xz_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 307);

    auto ta2_xz_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 308);

    auto ta2_xz_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 309);

    auto ta2_xz_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 310);

    auto ta2_xz_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 311);

    auto ta2_xz_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 312);

    auto ta2_xz_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 313);

    auto ta2_xz_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 314);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xx_xxxx_1, ta1_z_xx_xxxy_1, ta1_z_xx_xxxz_1, ta1_z_xx_xxyy_1, ta1_z_xx_xxyz_1, ta1_z_xx_xxzz_1, ta1_z_xx_xyyy_1, ta1_z_xx_xyyz_1, ta1_z_xx_xyzz_1, ta1_z_xx_xzzz_1, ta1_z_xx_yyyy_1, ta1_z_xx_yyyz_1, ta1_z_xx_yyzz_1, ta1_z_xx_yzzz_1, ta1_z_xx_zzzz_1, ta2_xz_x_xxxx_0, ta2_xz_x_xxxx_1, ta2_xz_x_xxxy_0, ta2_xz_x_xxxy_1, ta2_xz_x_xxxz_0, ta2_xz_x_xxxz_1, ta2_xz_x_xxyy_0, ta2_xz_x_xxyy_1, ta2_xz_x_xxyz_0, ta2_xz_x_xxyz_1, ta2_xz_x_xxzz_0, ta2_xz_x_xxzz_1, ta2_xz_x_xyyy_0, ta2_xz_x_xyyy_1, ta2_xz_x_xyyz_0, ta2_xz_x_xyyz_1, ta2_xz_x_xyzz_0, ta2_xz_x_xyzz_1, ta2_xz_x_xzzz_0, ta2_xz_x_xzzz_1, ta2_xz_x_yyyy_0, ta2_xz_x_yyyy_1, ta2_xz_x_yyyz_0, ta2_xz_x_yyyz_1, ta2_xz_x_yyzz_0, ta2_xz_x_yyzz_1, ta2_xz_x_yzzz_0, ta2_xz_x_yzzz_1, ta2_xz_x_zzzz_0, ta2_xz_x_zzzz_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxxx_0, ta2_xz_xx_xxxx_1, ta2_xz_xx_xxxy_0, ta2_xz_xx_xxxy_1, ta2_xz_xx_xxxz_0, ta2_xz_xx_xxxz_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxyy_0, ta2_xz_xx_xxyy_1, ta2_xz_xx_xxyz_0, ta2_xz_xx_xxyz_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xxzz_0, ta2_xz_xx_xxzz_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyyy_0, ta2_xz_xx_xyyy_1, ta2_xz_xx_xyyz_0, ta2_xz_xx_xyyz_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xyzz_0, ta2_xz_xx_xyzz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_xzzz_0, ta2_xz_xx_xzzz_1, ta2_xz_xx_yyy_0, ta2_xz_xx_yyy_1, ta2_xz_xx_yyyy_0, ta2_xz_xx_yyyy_1, ta2_xz_xx_yyyz_0, ta2_xz_xx_yyyz_1, ta2_xz_xx_yyz_0, ta2_xz_xx_yyz_1, ta2_xz_xx_yyzz_0, ta2_xz_xx_yyzz_1, ta2_xz_xx_yzz_0, ta2_xz_xx_yzz_1, ta2_xz_xx_yzzz_0, ta2_xz_xx_yzzz_1, ta2_xz_xx_zzz_0, ta2_xz_xx_zzz_1, ta2_xz_xx_zzzz_0, ta2_xz_xx_zzzz_1, ta2_xz_xxx_xxxx_0, ta2_xz_xxx_xxxy_0, ta2_xz_xxx_xxxz_0, ta2_xz_xxx_xxyy_0, ta2_xz_xxx_xxyz_0, ta2_xz_xxx_xxzz_0, ta2_xz_xxx_xyyy_0, ta2_xz_xxx_xyyz_0, ta2_xz_xxx_xyzz_0, ta2_xz_xxx_xzzz_0, ta2_xz_xxx_yyyy_0, ta2_xz_xxx_yyyz_0, ta2_xz_xxx_yyzz_0, ta2_xz_xxx_yzzz_0, ta2_xz_xxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxx_xxxx_0[i] = 2.0 * ta2_xz_x_xxxx_0[i] * fe_0 - 2.0 * ta2_xz_x_xxxx_1[i] * fe_0 + 4.0 * ta2_xz_xx_xxx_0[i] * fe_0 - 4.0 * ta2_xz_xx_xxx_1[i] * fe_0 + ta1_z_xx_xxxx_1[i] + ta2_xz_xx_xxxx_0[i] * pa_x[i] - ta2_xz_xx_xxxx_1[i] * pc_x[i];

        ta2_xz_xxx_xxxy_0[i] = 2.0 * ta2_xz_x_xxxy_0[i] * fe_0 - 2.0 * ta2_xz_x_xxxy_1[i] * fe_0 + 3.0 * ta2_xz_xx_xxy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxy_1[i] * fe_0 + ta1_z_xx_xxxy_1[i] + ta2_xz_xx_xxxy_0[i] * pa_x[i] - ta2_xz_xx_xxxy_1[i] * pc_x[i];

        ta2_xz_xxx_xxxz_0[i] = 2.0 * ta2_xz_x_xxxz_0[i] * fe_0 - 2.0 * ta2_xz_x_xxxz_1[i] * fe_0 + 3.0 * ta2_xz_xx_xxz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxz_1[i] * fe_0 + ta1_z_xx_xxxz_1[i] + ta2_xz_xx_xxxz_0[i] * pa_x[i] - ta2_xz_xx_xxxz_1[i] * pc_x[i];

        ta2_xz_xxx_xxyy_0[i] = 2.0 * ta2_xz_x_xxyy_0[i] * fe_0 - 2.0 * ta2_xz_x_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_xx_xyy_0[i] * fe_0 - 2.0 * ta2_xz_xx_xyy_1[i] * fe_0 + ta1_z_xx_xxyy_1[i] + ta2_xz_xx_xxyy_0[i] * pa_x[i] - ta2_xz_xx_xxyy_1[i] * pc_x[i];

        ta2_xz_xxx_xxyz_0[i] = 2.0 * ta2_xz_x_xxyz_0[i] * fe_0 - 2.0 * ta2_xz_x_xxyz_1[i] * fe_0 + 2.0 * ta2_xz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xyz_1[i] * fe_0 + ta1_z_xx_xxyz_1[i] + ta2_xz_xx_xxyz_0[i] * pa_x[i] - ta2_xz_xx_xxyz_1[i] * pc_x[i];

        ta2_xz_xxx_xxzz_0[i] = 2.0 * ta2_xz_x_xxzz_0[i] * fe_0 - 2.0 * ta2_xz_x_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_xx_xzz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xzz_1[i] * fe_0 + ta1_z_xx_xxzz_1[i] + ta2_xz_xx_xxzz_0[i] * pa_x[i] - ta2_xz_xx_xxzz_1[i] * pc_x[i];

        ta2_xz_xxx_xyyy_0[i] = 2.0 * ta2_xz_x_xyyy_0[i] * fe_0 - 2.0 * ta2_xz_x_xyyy_1[i] * fe_0 + ta2_xz_xx_yyy_0[i] * fe_0 - ta2_xz_xx_yyy_1[i] * fe_0 + ta1_z_xx_xyyy_1[i] + ta2_xz_xx_xyyy_0[i] * pa_x[i] - ta2_xz_xx_xyyy_1[i] * pc_x[i];

        ta2_xz_xxx_xyyz_0[i] = 2.0 * ta2_xz_x_xyyz_0[i] * fe_0 - 2.0 * ta2_xz_x_xyyz_1[i] * fe_0 + ta2_xz_xx_yyz_0[i] * fe_0 - ta2_xz_xx_yyz_1[i] * fe_0 + ta1_z_xx_xyyz_1[i] + ta2_xz_xx_xyyz_0[i] * pa_x[i] - ta2_xz_xx_xyyz_1[i] * pc_x[i];

        ta2_xz_xxx_xyzz_0[i] = 2.0 * ta2_xz_x_xyzz_0[i] * fe_0 - 2.0 * ta2_xz_x_xyzz_1[i] * fe_0 + ta2_xz_xx_yzz_0[i] * fe_0 - ta2_xz_xx_yzz_1[i] * fe_0 + ta1_z_xx_xyzz_1[i] + ta2_xz_xx_xyzz_0[i] * pa_x[i] - ta2_xz_xx_xyzz_1[i] * pc_x[i];

        ta2_xz_xxx_xzzz_0[i] = 2.0 * ta2_xz_x_xzzz_0[i] * fe_0 - 2.0 * ta2_xz_x_xzzz_1[i] * fe_0 + ta2_xz_xx_zzz_0[i] * fe_0 - ta2_xz_xx_zzz_1[i] * fe_0 + ta1_z_xx_xzzz_1[i] + ta2_xz_xx_xzzz_0[i] * pa_x[i] - ta2_xz_xx_xzzz_1[i] * pc_x[i];

        ta2_xz_xxx_yyyy_0[i] = 2.0 * ta2_xz_x_yyyy_0[i] * fe_0 - 2.0 * ta2_xz_x_yyyy_1[i] * fe_0 + ta1_z_xx_yyyy_1[i] + ta2_xz_xx_yyyy_0[i] * pa_x[i] - ta2_xz_xx_yyyy_1[i] * pc_x[i];

        ta2_xz_xxx_yyyz_0[i] = 2.0 * ta2_xz_x_yyyz_0[i] * fe_0 - 2.0 * ta2_xz_x_yyyz_1[i] * fe_0 + ta1_z_xx_yyyz_1[i] + ta2_xz_xx_yyyz_0[i] * pa_x[i] - ta2_xz_xx_yyyz_1[i] * pc_x[i];

        ta2_xz_xxx_yyzz_0[i] = 2.0 * ta2_xz_x_yyzz_0[i] * fe_0 - 2.0 * ta2_xz_x_yyzz_1[i] * fe_0 + ta1_z_xx_yyzz_1[i] + ta2_xz_xx_yyzz_0[i] * pa_x[i] - ta2_xz_xx_yyzz_1[i] * pc_x[i];

        ta2_xz_xxx_yzzz_0[i] = 2.0 * ta2_xz_x_yzzz_0[i] * fe_0 - 2.0 * ta2_xz_x_yzzz_1[i] * fe_0 + ta1_z_xx_yzzz_1[i] + ta2_xz_xx_yzzz_0[i] * pa_x[i] - ta2_xz_xx_yzzz_1[i] * pc_x[i];

        ta2_xz_xxx_zzzz_0[i] = 2.0 * ta2_xz_x_zzzz_0[i] * fe_0 - 2.0 * ta2_xz_x_zzzz_1[i] * fe_0 + ta1_z_xx_zzzz_1[i] + ta2_xz_xx_zzzz_0[i] * pa_x[i] - ta2_xz_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : FG

    auto ta2_xz_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 315);

    auto ta2_xz_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 316);

    auto ta2_xz_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 317);

    auto ta2_xz_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 318);

    auto ta2_xz_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 319);

    auto ta2_xz_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 320);

    auto ta2_xz_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 321);

    auto ta2_xz_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 322);

    auto ta2_xz_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 323);

    auto ta2_xz_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 324);

    auto ta2_xz_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 325);

    auto ta2_xz_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 326);

    auto ta2_xz_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 327);

    auto ta2_xz_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 328);

    auto ta2_xz_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 329);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxxx_0, ta2_xz_xx_xxxx_1, ta2_xz_xx_xxxy_0, ta2_xz_xx_xxxy_1, ta2_xz_xx_xxxz_0, ta2_xz_xx_xxxz_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxyy_0, ta2_xz_xx_xxyy_1, ta2_xz_xx_xxyz_0, ta2_xz_xx_xxyz_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xxzz_0, ta2_xz_xx_xxzz_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyyy_0, ta2_xz_xx_xyyy_1, ta2_xz_xx_xyyz_0, ta2_xz_xx_xyyz_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xyzz_0, ta2_xz_xx_xyzz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_xzzz_0, ta2_xz_xx_xzzz_1, ta2_xz_xx_yyy_0, ta2_xz_xx_yyy_1, ta2_xz_xx_yyyy_0, ta2_xz_xx_yyyy_1, ta2_xz_xx_yyyz_0, ta2_xz_xx_yyyz_1, ta2_xz_xx_yyz_0, ta2_xz_xx_yyz_1, ta2_xz_xx_yyzz_0, ta2_xz_xx_yyzz_1, ta2_xz_xx_yzz_0, ta2_xz_xx_yzz_1, ta2_xz_xx_yzzz_0, ta2_xz_xx_yzzz_1, ta2_xz_xx_zzz_0, ta2_xz_xx_zzz_1, ta2_xz_xx_zzzz_0, ta2_xz_xx_zzzz_1, ta2_xz_xxy_xxxx_0, ta2_xz_xxy_xxxy_0, ta2_xz_xxy_xxxz_0, ta2_xz_xxy_xxyy_0, ta2_xz_xxy_xxyz_0, ta2_xz_xxy_xxzz_0, ta2_xz_xxy_xyyy_0, ta2_xz_xxy_xyyz_0, ta2_xz_xxy_xyzz_0, ta2_xz_xxy_xzzz_0, ta2_xz_xxy_yyyy_0, ta2_xz_xxy_yyyz_0, ta2_xz_xxy_yyzz_0, ta2_xz_xxy_yzzz_0, ta2_xz_xxy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxy_xxxx_0[i] = ta2_xz_xx_xxxx_0[i] * pa_y[i] - ta2_xz_xx_xxxx_1[i] * pc_y[i];

        ta2_xz_xxy_xxxy_0[i] = ta2_xz_xx_xxx_0[i] * fe_0 - ta2_xz_xx_xxx_1[i] * fe_0 + ta2_xz_xx_xxxy_0[i] * pa_y[i] - ta2_xz_xx_xxxy_1[i] * pc_y[i];

        ta2_xz_xxy_xxxz_0[i] = ta2_xz_xx_xxxz_0[i] * pa_y[i] - ta2_xz_xx_xxxz_1[i] * pc_y[i];

        ta2_xz_xxy_xxyy_0[i] = 2.0 * ta2_xz_xx_xxy_0[i] * fe_0 - 2.0 * ta2_xz_xx_xxy_1[i] * fe_0 + ta2_xz_xx_xxyy_0[i] * pa_y[i] - ta2_xz_xx_xxyy_1[i] * pc_y[i];

        ta2_xz_xxy_xxyz_0[i] = ta2_xz_xx_xxz_0[i] * fe_0 - ta2_xz_xx_xxz_1[i] * fe_0 + ta2_xz_xx_xxyz_0[i] * pa_y[i] - ta2_xz_xx_xxyz_1[i] * pc_y[i];

        ta2_xz_xxy_xxzz_0[i] = ta2_xz_xx_xxzz_0[i] * pa_y[i] - ta2_xz_xx_xxzz_1[i] * pc_y[i];

        ta2_xz_xxy_xyyy_0[i] = 3.0 * ta2_xz_xx_xyy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xyy_1[i] * fe_0 + ta2_xz_xx_xyyy_0[i] * pa_y[i] - ta2_xz_xx_xyyy_1[i] * pc_y[i];

        ta2_xz_xxy_xyyz_0[i] = 2.0 * ta2_xz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xyz_1[i] * fe_0 + ta2_xz_xx_xyyz_0[i] * pa_y[i] - ta2_xz_xx_xyyz_1[i] * pc_y[i];

        ta2_xz_xxy_xyzz_0[i] = ta2_xz_xx_xzz_0[i] * fe_0 - ta2_xz_xx_xzz_1[i] * fe_0 + ta2_xz_xx_xyzz_0[i] * pa_y[i] - ta2_xz_xx_xyzz_1[i] * pc_y[i];

        ta2_xz_xxy_xzzz_0[i] = ta2_xz_xx_xzzz_0[i] * pa_y[i] - ta2_xz_xx_xzzz_1[i] * pc_y[i];

        ta2_xz_xxy_yyyy_0[i] = 4.0 * ta2_xz_xx_yyy_0[i] * fe_0 - 4.0 * ta2_xz_xx_yyy_1[i] * fe_0 + ta2_xz_xx_yyyy_0[i] * pa_y[i] - ta2_xz_xx_yyyy_1[i] * pc_y[i];

        ta2_xz_xxy_yyyz_0[i] = 3.0 * ta2_xz_xx_yyz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yyz_1[i] * fe_0 + ta2_xz_xx_yyyz_0[i] * pa_y[i] - ta2_xz_xx_yyyz_1[i] * pc_y[i];

        ta2_xz_xxy_yyzz_0[i] = 2.0 * ta2_xz_xx_yzz_0[i] * fe_0 - 2.0 * ta2_xz_xx_yzz_1[i] * fe_0 + ta2_xz_xx_yyzz_0[i] * pa_y[i] - ta2_xz_xx_yyzz_1[i] * pc_y[i];

        ta2_xz_xxy_yzzz_0[i] = ta2_xz_xx_zzz_0[i] * fe_0 - ta2_xz_xx_zzz_1[i] * fe_0 + ta2_xz_xx_yzzz_0[i] * pa_y[i] - ta2_xz_xx_yzzz_1[i] * pc_y[i];

        ta2_xz_xxy_zzzz_0[i] = ta2_xz_xx_zzzz_0[i] * pa_y[i] - ta2_xz_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 330-345 components of targeted buffer : FG

    auto ta2_xz_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 330);

    auto ta2_xz_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 331);

    auto ta2_xz_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 332);

    auto ta2_xz_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 333);

    auto ta2_xz_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 334);

    auto ta2_xz_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 335);

    auto ta2_xz_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 336);

    auto ta2_xz_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 337);

    auto ta2_xz_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 338);

    auto ta2_xz_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 339);

    auto ta2_xz_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 340);

    auto ta2_xz_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 341);

    auto ta2_xz_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 342);

    auto ta2_xz_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 343);

    auto ta2_xz_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 344);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xx_xxxx_1, ta1_x_xx_xxxy_1, ta1_x_xx_xxxz_1, ta1_x_xx_xxyy_1, ta1_x_xx_xxyz_1, ta1_x_xx_xxzz_1, ta1_x_xx_xyyy_1, ta1_x_xx_xyyz_1, ta1_x_xx_xyzz_1, ta1_x_xx_xzzz_1, ta1_x_xx_yyyy_1, ta1_z_xz_yyyz_1, ta1_z_xz_yyzz_1, ta1_z_xz_yzzz_1, ta1_z_xz_zzzz_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxxx_0, ta2_xz_xx_xxxx_1, ta2_xz_xx_xxxy_0, ta2_xz_xx_xxxy_1, ta2_xz_xx_xxxz_0, ta2_xz_xx_xxxz_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxyy_0, ta2_xz_xx_xxyy_1, ta2_xz_xx_xxyz_0, ta2_xz_xx_xxyz_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xxzz_0, ta2_xz_xx_xxzz_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyyy_0, ta2_xz_xx_xyyy_1, ta2_xz_xx_xyyz_0, ta2_xz_xx_xyyz_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xyzz_0, ta2_xz_xx_xyzz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_xzzz_0, ta2_xz_xx_xzzz_1, ta2_xz_xx_yyyy_0, ta2_xz_xx_yyyy_1, ta2_xz_xxz_xxxx_0, ta2_xz_xxz_xxxy_0, ta2_xz_xxz_xxxz_0, ta2_xz_xxz_xxyy_0, ta2_xz_xxz_xxyz_0, ta2_xz_xxz_xxzz_0, ta2_xz_xxz_xyyy_0, ta2_xz_xxz_xyyz_0, ta2_xz_xxz_xyzz_0, ta2_xz_xxz_xzzz_0, ta2_xz_xxz_yyyy_0, ta2_xz_xxz_yyyz_0, ta2_xz_xxz_yyzz_0, ta2_xz_xxz_yzzz_0, ta2_xz_xxz_zzzz_0, ta2_xz_xz_yyyz_0, ta2_xz_xz_yyyz_1, ta2_xz_xz_yyzz_0, ta2_xz_xz_yyzz_1, ta2_xz_xz_yzzz_0, ta2_xz_xz_yzzz_1, ta2_xz_xz_zzzz_0, ta2_xz_xz_zzzz_1, ta2_xz_z_yyyz_0, ta2_xz_z_yyyz_1, ta2_xz_z_yyzz_0, ta2_xz_z_yyzz_1, ta2_xz_z_yzzz_0, ta2_xz_z_yzzz_1, ta2_xz_z_zzzz_0, ta2_xz_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxz_xxxx_0[i] = ta1_x_xx_xxxx_1[i] + ta2_xz_xx_xxxx_0[i] * pa_z[i] - ta2_xz_xx_xxxx_1[i] * pc_z[i];

        ta2_xz_xxz_xxxy_0[i] = ta1_x_xx_xxxy_1[i] + ta2_xz_xx_xxxy_0[i] * pa_z[i] - ta2_xz_xx_xxxy_1[i] * pc_z[i];

        ta2_xz_xxz_xxxz_0[i] = ta2_xz_xx_xxx_0[i] * fe_0 - ta2_xz_xx_xxx_1[i] * fe_0 + ta1_x_xx_xxxz_1[i] + ta2_xz_xx_xxxz_0[i] * pa_z[i] - ta2_xz_xx_xxxz_1[i] * pc_z[i];

        ta2_xz_xxz_xxyy_0[i] = ta1_x_xx_xxyy_1[i] + ta2_xz_xx_xxyy_0[i] * pa_z[i] - ta2_xz_xx_xxyy_1[i] * pc_z[i];

        ta2_xz_xxz_xxyz_0[i] = ta2_xz_xx_xxy_0[i] * fe_0 - ta2_xz_xx_xxy_1[i] * fe_0 + ta1_x_xx_xxyz_1[i] + ta2_xz_xx_xxyz_0[i] * pa_z[i] - ta2_xz_xx_xxyz_1[i] * pc_z[i];

        ta2_xz_xxz_xxzz_0[i] = 2.0 * ta2_xz_xx_xxz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xxz_1[i] * fe_0 + ta1_x_xx_xxzz_1[i] + ta2_xz_xx_xxzz_0[i] * pa_z[i] - ta2_xz_xx_xxzz_1[i] * pc_z[i];

        ta2_xz_xxz_xyyy_0[i] = ta1_x_xx_xyyy_1[i] + ta2_xz_xx_xyyy_0[i] * pa_z[i] - ta2_xz_xx_xyyy_1[i] * pc_z[i];

        ta2_xz_xxz_xyyz_0[i] = ta2_xz_xx_xyy_0[i] * fe_0 - ta2_xz_xx_xyy_1[i] * fe_0 + ta1_x_xx_xyyz_1[i] + ta2_xz_xx_xyyz_0[i] * pa_z[i] - ta2_xz_xx_xyyz_1[i] * pc_z[i];

        ta2_xz_xxz_xyzz_0[i] = 2.0 * ta2_xz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xyz_1[i] * fe_0 + ta1_x_xx_xyzz_1[i] + ta2_xz_xx_xyzz_0[i] * pa_z[i] - ta2_xz_xx_xyzz_1[i] * pc_z[i];

        ta2_xz_xxz_xzzz_0[i] = 3.0 * ta2_xz_xx_xzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xzz_1[i] * fe_0 + ta1_x_xx_xzzz_1[i] + ta2_xz_xx_xzzz_0[i] * pa_z[i] - ta2_xz_xx_xzzz_1[i] * pc_z[i];

        ta2_xz_xxz_yyyy_0[i] = ta1_x_xx_yyyy_1[i] + ta2_xz_xx_yyyy_0[i] * pa_z[i] - ta2_xz_xx_yyyy_1[i] * pc_z[i];

        ta2_xz_xxz_yyyz_0[i] = ta2_xz_z_yyyz_0[i] * fe_0 - ta2_xz_z_yyyz_1[i] * fe_0 + ta1_z_xz_yyyz_1[i] + ta2_xz_xz_yyyz_0[i] * pa_x[i] - ta2_xz_xz_yyyz_1[i] * pc_x[i];

        ta2_xz_xxz_yyzz_0[i] = ta2_xz_z_yyzz_0[i] * fe_0 - ta2_xz_z_yyzz_1[i] * fe_0 + ta1_z_xz_yyzz_1[i] + ta2_xz_xz_yyzz_0[i] * pa_x[i] - ta2_xz_xz_yyzz_1[i] * pc_x[i];

        ta2_xz_xxz_yzzz_0[i] = ta2_xz_z_yzzz_0[i] * fe_0 - ta2_xz_z_yzzz_1[i] * fe_0 + ta1_z_xz_yzzz_1[i] + ta2_xz_xz_yzzz_0[i] * pa_x[i] - ta2_xz_xz_yzzz_1[i] * pc_x[i];

        ta2_xz_xxz_zzzz_0[i] = ta2_xz_z_zzzz_0[i] * fe_0 - ta2_xz_z_zzzz_1[i] * fe_0 + ta1_z_xz_zzzz_1[i] + ta2_xz_xz_zzzz_0[i] * pa_x[i] - ta2_xz_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 345-360 components of targeted buffer : FG

    auto ta2_xz_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 345);

    auto ta2_xz_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 346);

    auto ta2_xz_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 347);

    auto ta2_xz_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 348);

    auto ta2_xz_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 349);

    auto ta2_xz_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 350);

    auto ta2_xz_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 351);

    auto ta2_xz_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 352);

    auto ta2_xz_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 353);

    auto ta2_xz_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 354);

    auto ta2_xz_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 355);

    auto ta2_xz_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 356);

    auto ta2_xz_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 357);

    auto ta2_xz_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 358);

    auto ta2_xz_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 359);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_yy_xxxy_1, ta1_z_yy_xxyy_1, ta1_z_yy_xxyz_1, ta1_z_yy_xyyy_1, ta1_z_yy_xyyz_1, ta1_z_yy_xyzz_1, ta1_z_yy_yyyy_1, ta1_z_yy_yyyz_1, ta1_z_yy_yyzz_1, ta1_z_yy_yzzz_1, ta1_z_yy_zzzz_1, ta2_xz_x_xxxx_0, ta2_xz_x_xxxx_1, ta2_xz_x_xxxz_0, ta2_xz_x_xxxz_1, ta2_xz_x_xxzz_0, ta2_xz_x_xxzz_1, ta2_xz_x_xzzz_0, ta2_xz_x_xzzz_1, ta2_xz_xy_xxxx_0, ta2_xz_xy_xxxx_1, ta2_xz_xy_xxxz_0, ta2_xz_xy_xxxz_1, ta2_xz_xy_xxzz_0, ta2_xz_xy_xxzz_1, ta2_xz_xy_xzzz_0, ta2_xz_xy_xzzz_1, ta2_xz_xyy_xxxx_0, ta2_xz_xyy_xxxy_0, ta2_xz_xyy_xxxz_0, ta2_xz_xyy_xxyy_0, ta2_xz_xyy_xxyz_0, ta2_xz_xyy_xxzz_0, ta2_xz_xyy_xyyy_0, ta2_xz_xyy_xyyz_0, ta2_xz_xyy_xyzz_0, ta2_xz_xyy_xzzz_0, ta2_xz_xyy_yyyy_0, ta2_xz_xyy_yyyz_0, ta2_xz_xyy_yyzz_0, ta2_xz_xyy_yzzz_0, ta2_xz_xyy_zzzz_0, ta2_xz_yy_xxxy_0, ta2_xz_yy_xxxy_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xxyy_0, ta2_xz_yy_xxyy_1, ta2_xz_yy_xxyz_0, ta2_xz_yy_xxyz_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyyy_0, ta2_xz_yy_xyyy_1, ta2_xz_yy_xyyz_0, ta2_xz_yy_xyyz_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_xyzz_0, ta2_xz_yy_xyzz_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyyy_0, ta2_xz_yy_yyyy_1, ta2_xz_yy_yyyz_0, ta2_xz_yy_yyyz_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yyzz_0, ta2_xz_yy_yyzz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yy_yzzz_0, ta2_xz_yy_yzzz_1, ta2_xz_yy_zzzz_0, ta2_xz_yy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyy_xxxx_0[i] = ta2_xz_x_xxxx_0[i] * fe_0 - ta2_xz_x_xxxx_1[i] * fe_0 + ta2_xz_xy_xxxx_0[i] * pa_y[i] - ta2_xz_xy_xxxx_1[i] * pc_y[i];

        ta2_xz_xyy_xxxy_0[i] = 3.0 * ta2_xz_yy_xxy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxy_1[i] * fe_0 + ta1_z_yy_xxxy_1[i] + ta2_xz_yy_xxxy_0[i] * pa_x[i] - ta2_xz_yy_xxxy_1[i] * pc_x[i];

        ta2_xz_xyy_xxxz_0[i] = ta2_xz_x_xxxz_0[i] * fe_0 - ta2_xz_x_xxxz_1[i] * fe_0 + ta2_xz_xy_xxxz_0[i] * pa_y[i] - ta2_xz_xy_xxxz_1[i] * pc_y[i];

        ta2_xz_xyy_xxyy_0[i] = 2.0 * ta2_xz_yy_xyy_0[i] * fe_0 - 2.0 * ta2_xz_yy_xyy_1[i] * fe_0 + ta1_z_yy_xxyy_1[i] + ta2_xz_yy_xxyy_0[i] * pa_x[i] - ta2_xz_yy_xxyy_1[i] * pc_x[i];

        ta2_xz_xyy_xxyz_0[i] = 2.0 * ta2_xz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xz_yy_xyz_1[i] * fe_0 + ta1_z_yy_xxyz_1[i] + ta2_xz_yy_xxyz_0[i] * pa_x[i] - ta2_xz_yy_xxyz_1[i] * pc_x[i];

        ta2_xz_xyy_xxzz_0[i] = ta2_xz_x_xxzz_0[i] * fe_0 - ta2_xz_x_xxzz_1[i] * fe_0 + ta2_xz_xy_xxzz_0[i] * pa_y[i] - ta2_xz_xy_xxzz_1[i] * pc_y[i];

        ta2_xz_xyy_xyyy_0[i] = ta2_xz_yy_yyy_0[i] * fe_0 - ta2_xz_yy_yyy_1[i] * fe_0 + ta1_z_yy_xyyy_1[i] + ta2_xz_yy_xyyy_0[i] * pa_x[i] - ta2_xz_yy_xyyy_1[i] * pc_x[i];

        ta2_xz_xyy_xyyz_0[i] = ta2_xz_yy_yyz_0[i] * fe_0 - ta2_xz_yy_yyz_1[i] * fe_0 + ta1_z_yy_xyyz_1[i] + ta2_xz_yy_xyyz_0[i] * pa_x[i] - ta2_xz_yy_xyyz_1[i] * pc_x[i];

        ta2_xz_xyy_xyzz_0[i] = ta2_xz_yy_yzz_0[i] * fe_0 - ta2_xz_yy_yzz_1[i] * fe_0 + ta1_z_yy_xyzz_1[i] + ta2_xz_yy_xyzz_0[i] * pa_x[i] - ta2_xz_yy_xyzz_1[i] * pc_x[i];

        ta2_xz_xyy_xzzz_0[i] = ta2_xz_x_xzzz_0[i] * fe_0 - ta2_xz_x_xzzz_1[i] * fe_0 + ta2_xz_xy_xzzz_0[i] * pa_y[i] - ta2_xz_xy_xzzz_1[i] * pc_y[i];

        ta2_xz_xyy_yyyy_0[i] = ta1_z_yy_yyyy_1[i] + ta2_xz_yy_yyyy_0[i] * pa_x[i] - ta2_xz_yy_yyyy_1[i] * pc_x[i];

        ta2_xz_xyy_yyyz_0[i] = ta1_z_yy_yyyz_1[i] + ta2_xz_yy_yyyz_0[i] * pa_x[i] - ta2_xz_yy_yyyz_1[i] * pc_x[i];

        ta2_xz_xyy_yyzz_0[i] = ta1_z_yy_yyzz_1[i] + ta2_xz_yy_yyzz_0[i] * pa_x[i] - ta2_xz_yy_yyzz_1[i] * pc_x[i];

        ta2_xz_xyy_yzzz_0[i] = ta1_z_yy_yzzz_1[i] + ta2_xz_yy_yzzz_0[i] * pa_x[i] - ta2_xz_yy_yzzz_1[i] * pc_x[i];

        ta2_xz_xyy_zzzz_0[i] = ta1_z_yy_zzzz_1[i] + ta2_xz_yy_zzzz_0[i] * pa_x[i] - ta2_xz_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 360-375 components of targeted buffer : FG

    auto ta2_xz_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 360);

    auto ta2_xz_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 361);

    auto ta2_xz_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 362);

    auto ta2_xz_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 363);

    auto ta2_xz_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 364);

    auto ta2_xz_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 365);

    auto ta2_xz_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 366);

    auto ta2_xz_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 367);

    auto ta2_xz_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 368);

    auto ta2_xz_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 369);

    auto ta2_xz_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 370);

    auto ta2_xz_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 371);

    auto ta2_xz_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 372);

    auto ta2_xz_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 373);

    auto ta2_xz_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 374);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xy_xxxy_1, ta1_x_xy_xxyy_1, ta1_x_xy_xyyy_1, ta1_z_yz_yyyy_1, ta1_z_yz_yyyz_1, ta1_z_yz_yyzz_1, ta1_z_yz_yzzz_1, ta2_xz_xy_xxxy_0, ta2_xz_xy_xxxy_1, ta2_xz_xy_xxyy_0, ta2_xz_xy_xxyy_1, ta2_xz_xy_xyyy_0, ta2_xz_xy_xyyy_1, ta2_xz_xyz_xxxx_0, ta2_xz_xyz_xxxy_0, ta2_xz_xyz_xxxz_0, ta2_xz_xyz_xxyy_0, ta2_xz_xyz_xxyz_0, ta2_xz_xyz_xxzz_0, ta2_xz_xyz_xyyy_0, ta2_xz_xyz_xyyz_0, ta2_xz_xyz_xyzz_0, ta2_xz_xyz_xzzz_0, ta2_xz_xyz_yyyy_0, ta2_xz_xyz_yyyz_0, ta2_xz_xyz_yyzz_0, ta2_xz_xyz_yzzz_0, ta2_xz_xyz_zzzz_0, ta2_xz_xz_xxxx_0, ta2_xz_xz_xxxx_1, ta2_xz_xz_xxxz_0, ta2_xz_xz_xxxz_1, ta2_xz_xz_xxyz_0, ta2_xz_xz_xxyz_1, ta2_xz_xz_xxz_0, ta2_xz_xz_xxz_1, ta2_xz_xz_xxzz_0, ta2_xz_xz_xxzz_1, ta2_xz_xz_xyyz_0, ta2_xz_xz_xyyz_1, ta2_xz_xz_xyz_0, ta2_xz_xz_xyz_1, ta2_xz_xz_xyzz_0, ta2_xz_xz_xyzz_1, ta2_xz_xz_xzz_0, ta2_xz_xz_xzz_1, ta2_xz_xz_xzzz_0, ta2_xz_xz_xzzz_1, ta2_xz_xz_zzzz_0, ta2_xz_xz_zzzz_1, ta2_xz_yz_yyyy_0, ta2_xz_yz_yyyy_1, ta2_xz_yz_yyyz_0, ta2_xz_yz_yyyz_1, ta2_xz_yz_yyzz_0, ta2_xz_yz_yyzz_1, ta2_xz_yz_yzzz_0, ta2_xz_yz_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyz_xxxx_0[i] = ta2_xz_xz_xxxx_0[i] * pa_y[i] - ta2_xz_xz_xxxx_1[i] * pc_y[i];

        ta2_xz_xyz_xxxy_0[i] = ta1_x_xy_xxxy_1[i] + ta2_xz_xy_xxxy_0[i] * pa_z[i] - ta2_xz_xy_xxxy_1[i] * pc_z[i];

        ta2_xz_xyz_xxxz_0[i] = ta2_xz_xz_xxxz_0[i] * pa_y[i] - ta2_xz_xz_xxxz_1[i] * pc_y[i];

        ta2_xz_xyz_xxyy_0[i] = ta1_x_xy_xxyy_1[i] + ta2_xz_xy_xxyy_0[i] * pa_z[i] - ta2_xz_xy_xxyy_1[i] * pc_z[i];

        ta2_xz_xyz_xxyz_0[i] = ta2_xz_xz_xxz_0[i] * fe_0 - ta2_xz_xz_xxz_1[i] * fe_0 + ta2_xz_xz_xxyz_0[i] * pa_y[i] - ta2_xz_xz_xxyz_1[i] * pc_y[i];

        ta2_xz_xyz_xxzz_0[i] = ta2_xz_xz_xxzz_0[i] * pa_y[i] - ta2_xz_xz_xxzz_1[i] * pc_y[i];

        ta2_xz_xyz_xyyy_0[i] = ta1_x_xy_xyyy_1[i] + ta2_xz_xy_xyyy_0[i] * pa_z[i] - ta2_xz_xy_xyyy_1[i] * pc_z[i];

        ta2_xz_xyz_xyyz_0[i] = 2.0 * ta2_xz_xz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xz_xyz_1[i] * fe_0 + ta2_xz_xz_xyyz_0[i] * pa_y[i] - ta2_xz_xz_xyyz_1[i] * pc_y[i];

        ta2_xz_xyz_xyzz_0[i] = ta2_xz_xz_xzz_0[i] * fe_0 - ta2_xz_xz_xzz_1[i] * fe_0 + ta2_xz_xz_xyzz_0[i] * pa_y[i] - ta2_xz_xz_xyzz_1[i] * pc_y[i];

        ta2_xz_xyz_xzzz_0[i] = ta2_xz_xz_xzzz_0[i] * pa_y[i] - ta2_xz_xz_xzzz_1[i] * pc_y[i];

        ta2_xz_xyz_yyyy_0[i] = ta1_z_yz_yyyy_1[i] + ta2_xz_yz_yyyy_0[i] * pa_x[i] - ta2_xz_yz_yyyy_1[i] * pc_x[i];

        ta2_xz_xyz_yyyz_0[i] = ta1_z_yz_yyyz_1[i] + ta2_xz_yz_yyyz_0[i] * pa_x[i] - ta2_xz_yz_yyyz_1[i] * pc_x[i];

        ta2_xz_xyz_yyzz_0[i] = ta1_z_yz_yyzz_1[i] + ta2_xz_yz_yyzz_0[i] * pa_x[i] - ta2_xz_yz_yyzz_1[i] * pc_x[i];

        ta2_xz_xyz_yzzz_0[i] = ta1_z_yz_yzzz_1[i] + ta2_xz_yz_yzzz_0[i] * pa_x[i] - ta2_xz_yz_yzzz_1[i] * pc_x[i];

        ta2_xz_xyz_zzzz_0[i] = ta2_xz_xz_zzzz_0[i] * pa_y[i] - ta2_xz_xz_zzzz_1[i] * pc_y[i];
    }

    // Set up 375-390 components of targeted buffer : FG

    auto ta2_xz_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 375);

    auto ta2_xz_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 376);

    auto ta2_xz_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 377);

    auto ta2_xz_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 378);

    auto ta2_xz_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 379);

    auto ta2_xz_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 380);

    auto ta2_xz_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 381);

    auto ta2_xz_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 382);

    auto ta2_xz_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 383);

    auto ta2_xz_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 384);

    auto ta2_xz_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 385);

    auto ta2_xz_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 386);

    auto ta2_xz_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 387);

    auto ta2_xz_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 388);

    auto ta2_xz_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 389);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_zz_xxxx_1, ta1_z_zz_xxxy_1, ta1_z_zz_xxxz_1, ta1_z_zz_xxyy_1, ta1_z_zz_xxyz_1, ta1_z_zz_xxzz_1, ta1_z_zz_xyyy_1, ta1_z_zz_xyyz_1, ta1_z_zz_xyzz_1, ta1_z_zz_xzzz_1, ta1_z_zz_yyyy_1, ta1_z_zz_yyyz_1, ta1_z_zz_yyzz_1, ta1_z_zz_yzzz_1, ta1_z_zz_zzzz_1, ta2_xz_xzz_xxxx_0, ta2_xz_xzz_xxxy_0, ta2_xz_xzz_xxxz_0, ta2_xz_xzz_xxyy_0, ta2_xz_xzz_xxyz_0, ta2_xz_xzz_xxzz_0, ta2_xz_xzz_xyyy_0, ta2_xz_xzz_xyyz_0, ta2_xz_xzz_xyzz_0, ta2_xz_xzz_xzzz_0, ta2_xz_xzz_yyyy_0, ta2_xz_xzz_yyyz_0, ta2_xz_xzz_yyzz_0, ta2_xz_xzz_yzzz_0, ta2_xz_xzz_zzzz_0, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxxx_0, ta2_xz_zz_xxxx_1, ta2_xz_zz_xxxy_0, ta2_xz_zz_xxxy_1, ta2_xz_zz_xxxz_0, ta2_xz_zz_xxxz_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxyy_0, ta2_xz_zz_xxyy_1, ta2_xz_zz_xxyz_0, ta2_xz_zz_xxyz_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xxzz_0, ta2_xz_zz_xxzz_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyyy_0, ta2_xz_zz_xyyy_1, ta2_xz_zz_xyyz_0, ta2_xz_zz_xyyz_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xyzz_0, ta2_xz_zz_xyzz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_xzzz_0, ta2_xz_zz_xzzz_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyyy_0, ta2_xz_zz_yyyy_1, ta2_xz_zz_yyyz_0, ta2_xz_zz_yyyz_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yyzz_0, ta2_xz_zz_yyzz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_yzzz_0, ta2_xz_zz_yzzz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, ta2_xz_zz_zzzz_0, ta2_xz_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzz_xxxx_0[i] = 4.0 * ta2_xz_zz_xxx_0[i] * fe_0 - 4.0 * ta2_xz_zz_xxx_1[i] * fe_0 + ta1_z_zz_xxxx_1[i] + ta2_xz_zz_xxxx_0[i] * pa_x[i] - ta2_xz_zz_xxxx_1[i] * pc_x[i];

        ta2_xz_xzz_xxxy_0[i] = 3.0 * ta2_xz_zz_xxy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxy_1[i] * fe_0 + ta1_z_zz_xxxy_1[i] + ta2_xz_zz_xxxy_0[i] * pa_x[i] - ta2_xz_zz_xxxy_1[i] * pc_x[i];

        ta2_xz_xzz_xxxz_0[i] = 3.0 * ta2_xz_zz_xxz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxz_1[i] * fe_0 + ta1_z_zz_xxxz_1[i] + ta2_xz_zz_xxxz_0[i] * pa_x[i] - ta2_xz_zz_xxxz_1[i] * pc_x[i];

        ta2_xz_xzz_xxyy_0[i] = 2.0 * ta2_xz_zz_xyy_0[i] * fe_0 - 2.0 * ta2_xz_zz_xyy_1[i] * fe_0 + ta1_z_zz_xxyy_1[i] + ta2_xz_zz_xxyy_0[i] * pa_x[i] - ta2_xz_zz_xxyy_1[i] * pc_x[i];

        ta2_xz_xzz_xxyz_0[i] = 2.0 * ta2_xz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xyz_1[i] * fe_0 + ta1_z_zz_xxyz_1[i] + ta2_xz_zz_xxyz_0[i] * pa_x[i] - ta2_xz_zz_xxyz_1[i] * pc_x[i];

        ta2_xz_xzz_xxzz_0[i] = 2.0 * ta2_xz_zz_xzz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xzz_1[i] * fe_0 + ta1_z_zz_xxzz_1[i] + ta2_xz_zz_xxzz_0[i] * pa_x[i] - ta2_xz_zz_xxzz_1[i] * pc_x[i];

        ta2_xz_xzz_xyyy_0[i] = ta2_xz_zz_yyy_0[i] * fe_0 - ta2_xz_zz_yyy_1[i] * fe_0 + ta1_z_zz_xyyy_1[i] + ta2_xz_zz_xyyy_0[i] * pa_x[i] - ta2_xz_zz_xyyy_1[i] * pc_x[i];

        ta2_xz_xzz_xyyz_0[i] = ta2_xz_zz_yyz_0[i] * fe_0 - ta2_xz_zz_yyz_1[i] * fe_0 + ta1_z_zz_xyyz_1[i] + ta2_xz_zz_xyyz_0[i] * pa_x[i] - ta2_xz_zz_xyyz_1[i] * pc_x[i];

        ta2_xz_xzz_xyzz_0[i] = ta2_xz_zz_yzz_0[i] * fe_0 - ta2_xz_zz_yzz_1[i] * fe_0 + ta1_z_zz_xyzz_1[i] + ta2_xz_zz_xyzz_0[i] * pa_x[i] - ta2_xz_zz_xyzz_1[i] * pc_x[i];

        ta2_xz_xzz_xzzz_0[i] = ta2_xz_zz_zzz_0[i] * fe_0 - ta2_xz_zz_zzz_1[i] * fe_0 + ta1_z_zz_xzzz_1[i] + ta2_xz_zz_xzzz_0[i] * pa_x[i] - ta2_xz_zz_xzzz_1[i] * pc_x[i];

        ta2_xz_xzz_yyyy_0[i] = ta1_z_zz_yyyy_1[i] + ta2_xz_zz_yyyy_0[i] * pa_x[i] - ta2_xz_zz_yyyy_1[i] * pc_x[i];

        ta2_xz_xzz_yyyz_0[i] = ta1_z_zz_yyyz_1[i] + ta2_xz_zz_yyyz_0[i] * pa_x[i] - ta2_xz_zz_yyyz_1[i] * pc_x[i];

        ta2_xz_xzz_yyzz_0[i] = ta1_z_zz_yyzz_1[i] + ta2_xz_zz_yyzz_0[i] * pa_x[i] - ta2_xz_zz_yyzz_1[i] * pc_x[i];

        ta2_xz_xzz_yzzz_0[i] = ta1_z_zz_yzzz_1[i] + ta2_xz_zz_yzzz_0[i] * pa_x[i] - ta2_xz_zz_yzzz_1[i] * pc_x[i];

        ta2_xz_xzz_zzzz_0[i] = ta1_z_zz_zzzz_1[i] + ta2_xz_zz_zzzz_0[i] * pa_x[i] - ta2_xz_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 390-405 components of targeted buffer : FG

    auto ta2_xz_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 390);

    auto ta2_xz_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 391);

    auto ta2_xz_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 392);

    auto ta2_xz_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 393);

    auto ta2_xz_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 394);

    auto ta2_xz_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 395);

    auto ta2_xz_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 396);

    auto ta2_xz_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 397);

    auto ta2_xz_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 398);

    auto ta2_xz_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 399);

    auto ta2_xz_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 400);

    auto ta2_xz_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 401);

    auto ta2_xz_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 402);

    auto ta2_xz_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 403);

    auto ta2_xz_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 404);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_y_xxxx_0, ta2_xz_y_xxxx_1, ta2_xz_y_xxxy_0, ta2_xz_y_xxxy_1, ta2_xz_y_xxxz_0, ta2_xz_y_xxxz_1, ta2_xz_y_xxyy_0, ta2_xz_y_xxyy_1, ta2_xz_y_xxyz_0, ta2_xz_y_xxyz_1, ta2_xz_y_xxzz_0, ta2_xz_y_xxzz_1, ta2_xz_y_xyyy_0, ta2_xz_y_xyyy_1, ta2_xz_y_xyyz_0, ta2_xz_y_xyyz_1, ta2_xz_y_xyzz_0, ta2_xz_y_xyzz_1, ta2_xz_y_xzzz_0, ta2_xz_y_xzzz_1, ta2_xz_y_yyyy_0, ta2_xz_y_yyyy_1, ta2_xz_y_yyyz_0, ta2_xz_y_yyyz_1, ta2_xz_y_yyzz_0, ta2_xz_y_yyzz_1, ta2_xz_y_yzzz_0, ta2_xz_y_yzzz_1, ta2_xz_y_zzzz_0, ta2_xz_y_zzzz_1, ta2_xz_yy_xxx_0, ta2_xz_yy_xxx_1, ta2_xz_yy_xxxx_0, ta2_xz_yy_xxxx_1, ta2_xz_yy_xxxy_0, ta2_xz_yy_xxxy_1, ta2_xz_yy_xxxz_0, ta2_xz_yy_xxxz_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xxyy_0, ta2_xz_yy_xxyy_1, ta2_xz_yy_xxyz_0, ta2_xz_yy_xxyz_1, ta2_xz_yy_xxz_0, ta2_xz_yy_xxz_1, ta2_xz_yy_xxzz_0, ta2_xz_yy_xxzz_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyyy_0, ta2_xz_yy_xyyy_1, ta2_xz_yy_xyyz_0, ta2_xz_yy_xyyz_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_xyzz_0, ta2_xz_yy_xyzz_1, ta2_xz_yy_xzz_0, ta2_xz_yy_xzz_1, ta2_xz_yy_xzzz_0, ta2_xz_yy_xzzz_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyyy_0, ta2_xz_yy_yyyy_1, ta2_xz_yy_yyyz_0, ta2_xz_yy_yyyz_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yyzz_0, ta2_xz_yy_yyzz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yy_yzzz_0, ta2_xz_yy_yzzz_1, ta2_xz_yy_zzz_0, ta2_xz_yy_zzz_1, ta2_xz_yy_zzzz_0, ta2_xz_yy_zzzz_1, ta2_xz_yyy_xxxx_0, ta2_xz_yyy_xxxy_0, ta2_xz_yyy_xxxz_0, ta2_xz_yyy_xxyy_0, ta2_xz_yyy_xxyz_0, ta2_xz_yyy_xxzz_0, ta2_xz_yyy_xyyy_0, ta2_xz_yyy_xyyz_0, ta2_xz_yyy_xyzz_0, ta2_xz_yyy_xzzz_0, ta2_xz_yyy_yyyy_0, ta2_xz_yyy_yyyz_0, ta2_xz_yyy_yyzz_0, ta2_xz_yyy_yzzz_0, ta2_xz_yyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyy_xxxx_0[i] = 2.0 * ta2_xz_y_xxxx_0[i] * fe_0 - 2.0 * ta2_xz_y_xxxx_1[i] * fe_0 + ta2_xz_yy_xxxx_0[i] * pa_y[i] - ta2_xz_yy_xxxx_1[i] * pc_y[i];

        ta2_xz_yyy_xxxy_0[i] = 2.0 * ta2_xz_y_xxxy_0[i] * fe_0 - 2.0 * ta2_xz_y_xxxy_1[i] * fe_0 + ta2_xz_yy_xxx_0[i] * fe_0 - ta2_xz_yy_xxx_1[i] * fe_0 + ta2_xz_yy_xxxy_0[i] * pa_y[i] - ta2_xz_yy_xxxy_1[i] * pc_y[i];

        ta2_xz_yyy_xxxz_0[i] = 2.0 * ta2_xz_y_xxxz_0[i] * fe_0 - 2.0 * ta2_xz_y_xxxz_1[i] * fe_0 + ta2_xz_yy_xxxz_0[i] * pa_y[i] - ta2_xz_yy_xxxz_1[i] * pc_y[i];

        ta2_xz_yyy_xxyy_0[i] = 2.0 * ta2_xz_y_xxyy_0[i] * fe_0 - 2.0 * ta2_xz_y_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_yy_xxy_0[i] * fe_0 - 2.0 * ta2_xz_yy_xxy_1[i] * fe_0 + ta2_xz_yy_xxyy_0[i] * pa_y[i] - ta2_xz_yy_xxyy_1[i] * pc_y[i];

        ta2_xz_yyy_xxyz_0[i] = 2.0 * ta2_xz_y_xxyz_0[i] * fe_0 - 2.0 * ta2_xz_y_xxyz_1[i] * fe_0 + ta2_xz_yy_xxz_0[i] * fe_0 - ta2_xz_yy_xxz_1[i] * fe_0 + ta2_xz_yy_xxyz_0[i] * pa_y[i] - ta2_xz_yy_xxyz_1[i] * pc_y[i];

        ta2_xz_yyy_xxzz_0[i] = 2.0 * ta2_xz_y_xxzz_0[i] * fe_0 - 2.0 * ta2_xz_y_xxzz_1[i] * fe_0 + ta2_xz_yy_xxzz_0[i] * pa_y[i] - ta2_xz_yy_xxzz_1[i] * pc_y[i];

        ta2_xz_yyy_xyyy_0[i] = 2.0 * ta2_xz_y_xyyy_0[i] * fe_0 - 2.0 * ta2_xz_y_xyyy_1[i] * fe_0 + 3.0 * ta2_xz_yy_xyy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xyy_1[i] * fe_0 + ta2_xz_yy_xyyy_0[i] * pa_y[i] - ta2_xz_yy_xyyy_1[i] * pc_y[i];

        ta2_xz_yyy_xyyz_0[i] = 2.0 * ta2_xz_y_xyyz_0[i] * fe_0 - 2.0 * ta2_xz_y_xyyz_1[i] * fe_0 + 2.0 * ta2_xz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xz_yy_xyz_1[i] * fe_0 + ta2_xz_yy_xyyz_0[i] * pa_y[i] - ta2_xz_yy_xyyz_1[i] * pc_y[i];

        ta2_xz_yyy_xyzz_0[i] = 2.0 * ta2_xz_y_xyzz_0[i] * fe_0 - 2.0 * ta2_xz_y_xyzz_1[i] * fe_0 + ta2_xz_yy_xzz_0[i] * fe_0 - ta2_xz_yy_xzz_1[i] * fe_0 + ta2_xz_yy_xyzz_0[i] * pa_y[i] - ta2_xz_yy_xyzz_1[i] * pc_y[i];

        ta2_xz_yyy_xzzz_0[i] = 2.0 * ta2_xz_y_xzzz_0[i] * fe_0 - 2.0 * ta2_xz_y_xzzz_1[i] * fe_0 + ta2_xz_yy_xzzz_0[i] * pa_y[i] - ta2_xz_yy_xzzz_1[i] * pc_y[i];

        ta2_xz_yyy_yyyy_0[i] = 2.0 * ta2_xz_y_yyyy_0[i] * fe_0 - 2.0 * ta2_xz_y_yyyy_1[i] * fe_0 + 4.0 * ta2_xz_yy_yyy_0[i] * fe_0 - 4.0 * ta2_xz_yy_yyy_1[i] * fe_0 + ta2_xz_yy_yyyy_0[i] * pa_y[i] - ta2_xz_yy_yyyy_1[i] * pc_y[i];

        ta2_xz_yyy_yyyz_0[i] = 2.0 * ta2_xz_y_yyyz_0[i] * fe_0 - 2.0 * ta2_xz_y_yyyz_1[i] * fe_0 + 3.0 * ta2_xz_yy_yyz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yyz_1[i] * fe_0 + ta2_xz_yy_yyyz_0[i] * pa_y[i] - ta2_xz_yy_yyyz_1[i] * pc_y[i];

        ta2_xz_yyy_yyzz_0[i] = 2.0 * ta2_xz_y_yyzz_0[i] * fe_0 - 2.0 * ta2_xz_y_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_yy_yzz_0[i] * fe_0 - 2.0 * ta2_xz_yy_yzz_1[i] * fe_0 + ta2_xz_yy_yyzz_0[i] * pa_y[i] - ta2_xz_yy_yyzz_1[i] * pc_y[i];

        ta2_xz_yyy_yzzz_0[i] = 2.0 * ta2_xz_y_yzzz_0[i] * fe_0 - 2.0 * ta2_xz_y_yzzz_1[i] * fe_0 + ta2_xz_yy_zzz_0[i] * fe_0 - ta2_xz_yy_zzz_1[i] * fe_0 + ta2_xz_yy_yzzz_0[i] * pa_y[i] - ta2_xz_yy_yzzz_1[i] * pc_y[i];

        ta2_xz_yyy_zzzz_0[i] = 2.0 * ta2_xz_y_zzzz_0[i] * fe_0 - 2.0 * ta2_xz_y_zzzz_1[i] * fe_0 + ta2_xz_yy_zzzz_0[i] * pa_y[i] - ta2_xz_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 405-420 components of targeted buffer : FG

    auto ta2_xz_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 405);

    auto ta2_xz_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 406);

    auto ta2_xz_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 407);

    auto ta2_xz_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 408);

    auto ta2_xz_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 409);

    auto ta2_xz_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 410);

    auto ta2_xz_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 411);

    auto ta2_xz_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 412);

    auto ta2_xz_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 413);

    auto ta2_xz_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 414);

    auto ta2_xz_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 415);

    auto ta2_xz_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 416);

    auto ta2_xz_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 417);

    auto ta2_xz_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 418);

    auto ta2_xz_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 419);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yy_xxxx_1, ta1_x_yy_xxxy_1, ta1_x_yy_xxyy_1, ta1_x_yy_xxyz_1, ta1_x_yy_xyyy_1, ta1_x_yy_xyyz_1, ta1_x_yy_xyzz_1, ta1_x_yy_yyyy_1, ta1_x_yy_yyyz_1, ta1_x_yy_yyzz_1, ta1_x_yy_yzzz_1, ta2_xz_yy_xxxx_0, ta2_xz_yy_xxxx_1, ta2_xz_yy_xxxy_0, ta2_xz_yy_xxxy_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xxyy_0, ta2_xz_yy_xxyy_1, ta2_xz_yy_xxyz_0, ta2_xz_yy_xxyz_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyyy_0, ta2_xz_yy_xyyy_1, ta2_xz_yy_xyyz_0, ta2_xz_yy_xyyz_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_xyzz_0, ta2_xz_yy_xyzz_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyyy_0, ta2_xz_yy_yyyy_1, ta2_xz_yy_yyyz_0, ta2_xz_yy_yyyz_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yyzz_0, ta2_xz_yy_yyzz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yy_yzzz_0, ta2_xz_yy_yzzz_1, ta2_xz_yyz_xxxx_0, ta2_xz_yyz_xxxy_0, ta2_xz_yyz_xxxz_0, ta2_xz_yyz_xxyy_0, ta2_xz_yyz_xxyz_0, ta2_xz_yyz_xxzz_0, ta2_xz_yyz_xyyy_0, ta2_xz_yyz_xyyz_0, ta2_xz_yyz_xyzz_0, ta2_xz_yyz_xzzz_0, ta2_xz_yyz_yyyy_0, ta2_xz_yyz_yyyz_0, ta2_xz_yyz_yyzz_0, ta2_xz_yyz_yzzz_0, ta2_xz_yyz_zzzz_0, ta2_xz_yz_xxxz_0, ta2_xz_yz_xxxz_1, ta2_xz_yz_xxzz_0, ta2_xz_yz_xxzz_1, ta2_xz_yz_xzzz_0, ta2_xz_yz_xzzz_1, ta2_xz_yz_zzzz_0, ta2_xz_yz_zzzz_1, ta2_xz_z_xxxz_0, ta2_xz_z_xxxz_1, ta2_xz_z_xxzz_0, ta2_xz_z_xxzz_1, ta2_xz_z_xzzz_0, ta2_xz_z_xzzz_1, ta2_xz_z_zzzz_0, ta2_xz_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyz_xxxx_0[i] = ta1_x_yy_xxxx_1[i] + ta2_xz_yy_xxxx_0[i] * pa_z[i] - ta2_xz_yy_xxxx_1[i] * pc_z[i];

        ta2_xz_yyz_xxxy_0[i] = ta1_x_yy_xxxy_1[i] + ta2_xz_yy_xxxy_0[i] * pa_z[i] - ta2_xz_yy_xxxy_1[i] * pc_z[i];

        ta2_xz_yyz_xxxz_0[i] = ta2_xz_z_xxxz_0[i] * fe_0 - ta2_xz_z_xxxz_1[i] * fe_0 + ta2_xz_yz_xxxz_0[i] * pa_y[i] - ta2_xz_yz_xxxz_1[i] * pc_y[i];

        ta2_xz_yyz_xxyy_0[i] = ta1_x_yy_xxyy_1[i] + ta2_xz_yy_xxyy_0[i] * pa_z[i] - ta2_xz_yy_xxyy_1[i] * pc_z[i];

        ta2_xz_yyz_xxyz_0[i] = ta2_xz_yy_xxy_0[i] * fe_0 - ta2_xz_yy_xxy_1[i] * fe_0 + ta1_x_yy_xxyz_1[i] + ta2_xz_yy_xxyz_0[i] * pa_z[i] - ta2_xz_yy_xxyz_1[i] * pc_z[i];

        ta2_xz_yyz_xxzz_0[i] = ta2_xz_z_xxzz_0[i] * fe_0 - ta2_xz_z_xxzz_1[i] * fe_0 + ta2_xz_yz_xxzz_0[i] * pa_y[i] - ta2_xz_yz_xxzz_1[i] * pc_y[i];

        ta2_xz_yyz_xyyy_0[i] = ta1_x_yy_xyyy_1[i] + ta2_xz_yy_xyyy_0[i] * pa_z[i] - ta2_xz_yy_xyyy_1[i] * pc_z[i];

        ta2_xz_yyz_xyyz_0[i] = ta2_xz_yy_xyy_0[i] * fe_0 - ta2_xz_yy_xyy_1[i] * fe_0 + ta1_x_yy_xyyz_1[i] + ta2_xz_yy_xyyz_0[i] * pa_z[i] - ta2_xz_yy_xyyz_1[i] * pc_z[i];

        ta2_xz_yyz_xyzz_0[i] = 2.0 * ta2_xz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_xz_yy_xyz_1[i] * fe_0 + ta1_x_yy_xyzz_1[i] + ta2_xz_yy_xyzz_0[i] * pa_z[i] - ta2_xz_yy_xyzz_1[i] * pc_z[i];

        ta2_xz_yyz_xzzz_0[i] = ta2_xz_z_xzzz_0[i] * fe_0 - ta2_xz_z_xzzz_1[i] * fe_0 + ta2_xz_yz_xzzz_0[i] * pa_y[i] - ta2_xz_yz_xzzz_1[i] * pc_y[i];

        ta2_xz_yyz_yyyy_0[i] = ta1_x_yy_yyyy_1[i] + ta2_xz_yy_yyyy_0[i] * pa_z[i] - ta2_xz_yy_yyyy_1[i] * pc_z[i];

        ta2_xz_yyz_yyyz_0[i] = ta2_xz_yy_yyy_0[i] * fe_0 - ta2_xz_yy_yyy_1[i] * fe_0 + ta1_x_yy_yyyz_1[i] + ta2_xz_yy_yyyz_0[i] * pa_z[i] - ta2_xz_yy_yyyz_1[i] * pc_z[i];

        ta2_xz_yyz_yyzz_0[i] = 2.0 * ta2_xz_yy_yyz_0[i] * fe_0 - 2.0 * ta2_xz_yy_yyz_1[i] * fe_0 + ta1_x_yy_yyzz_1[i] + ta2_xz_yy_yyzz_0[i] * pa_z[i] - ta2_xz_yy_yyzz_1[i] * pc_z[i];

        ta2_xz_yyz_yzzz_0[i] = 3.0 * ta2_xz_yy_yzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yzz_1[i] * fe_0 + ta1_x_yy_yzzz_1[i] + ta2_xz_yy_yzzz_0[i] * pa_z[i] - ta2_xz_yy_yzzz_1[i] * pc_z[i];

        ta2_xz_yyz_zzzz_0[i] = ta2_xz_z_zzzz_0[i] * fe_0 - ta2_xz_z_zzzz_1[i] * fe_0 + ta2_xz_yz_zzzz_0[i] * pa_y[i] - ta2_xz_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 420-435 components of targeted buffer : FG

    auto ta2_xz_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 420);

    auto ta2_xz_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 421);

    auto ta2_xz_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 422);

    auto ta2_xz_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 423);

    auto ta2_xz_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 424);

    auto ta2_xz_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 425);

    auto ta2_xz_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 426);

    auto ta2_xz_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 427);

    auto ta2_xz_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 428);

    auto ta2_xz_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 429);

    auto ta2_xz_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 430);

    auto ta2_xz_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 431);

    auto ta2_xz_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 432);

    auto ta2_xz_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 433);

    auto ta2_xz_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 434);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_yzz_xxxx_0, ta2_xz_yzz_xxxy_0, ta2_xz_yzz_xxxz_0, ta2_xz_yzz_xxyy_0, ta2_xz_yzz_xxyz_0, ta2_xz_yzz_xxzz_0, ta2_xz_yzz_xyyy_0, ta2_xz_yzz_xyyz_0, ta2_xz_yzz_xyzz_0, ta2_xz_yzz_xzzz_0, ta2_xz_yzz_yyyy_0, ta2_xz_yzz_yyyz_0, ta2_xz_yzz_yyzz_0, ta2_xz_yzz_yzzz_0, ta2_xz_yzz_zzzz_0, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxxx_0, ta2_xz_zz_xxxx_1, ta2_xz_zz_xxxy_0, ta2_xz_zz_xxxy_1, ta2_xz_zz_xxxz_0, ta2_xz_zz_xxxz_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxyy_0, ta2_xz_zz_xxyy_1, ta2_xz_zz_xxyz_0, ta2_xz_zz_xxyz_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xxzz_0, ta2_xz_zz_xxzz_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyyy_0, ta2_xz_zz_xyyy_1, ta2_xz_zz_xyyz_0, ta2_xz_zz_xyyz_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xyzz_0, ta2_xz_zz_xyzz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_xzzz_0, ta2_xz_zz_xzzz_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyyy_0, ta2_xz_zz_yyyy_1, ta2_xz_zz_yyyz_0, ta2_xz_zz_yyyz_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yyzz_0, ta2_xz_zz_yyzz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_yzzz_0, ta2_xz_zz_yzzz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, ta2_xz_zz_zzzz_0, ta2_xz_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzz_xxxx_0[i] = ta2_xz_zz_xxxx_0[i] * pa_y[i] - ta2_xz_zz_xxxx_1[i] * pc_y[i];

        ta2_xz_yzz_xxxy_0[i] = ta2_xz_zz_xxx_0[i] * fe_0 - ta2_xz_zz_xxx_1[i] * fe_0 + ta2_xz_zz_xxxy_0[i] * pa_y[i] - ta2_xz_zz_xxxy_1[i] * pc_y[i];

        ta2_xz_yzz_xxxz_0[i] = ta2_xz_zz_xxxz_0[i] * pa_y[i] - ta2_xz_zz_xxxz_1[i] * pc_y[i];

        ta2_xz_yzz_xxyy_0[i] = 2.0 * ta2_xz_zz_xxy_0[i] * fe_0 - 2.0 * ta2_xz_zz_xxy_1[i] * fe_0 + ta2_xz_zz_xxyy_0[i] * pa_y[i] - ta2_xz_zz_xxyy_1[i] * pc_y[i];

        ta2_xz_yzz_xxyz_0[i] = ta2_xz_zz_xxz_0[i] * fe_0 - ta2_xz_zz_xxz_1[i] * fe_0 + ta2_xz_zz_xxyz_0[i] * pa_y[i] - ta2_xz_zz_xxyz_1[i] * pc_y[i];

        ta2_xz_yzz_xxzz_0[i] = ta2_xz_zz_xxzz_0[i] * pa_y[i] - ta2_xz_zz_xxzz_1[i] * pc_y[i];

        ta2_xz_yzz_xyyy_0[i] = 3.0 * ta2_xz_zz_xyy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xyy_1[i] * fe_0 + ta2_xz_zz_xyyy_0[i] * pa_y[i] - ta2_xz_zz_xyyy_1[i] * pc_y[i];

        ta2_xz_yzz_xyyz_0[i] = 2.0 * ta2_xz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xyz_1[i] * fe_0 + ta2_xz_zz_xyyz_0[i] * pa_y[i] - ta2_xz_zz_xyyz_1[i] * pc_y[i];

        ta2_xz_yzz_xyzz_0[i] = ta2_xz_zz_xzz_0[i] * fe_0 - ta2_xz_zz_xzz_1[i] * fe_0 + ta2_xz_zz_xyzz_0[i] * pa_y[i] - ta2_xz_zz_xyzz_1[i] * pc_y[i];

        ta2_xz_yzz_xzzz_0[i] = ta2_xz_zz_xzzz_0[i] * pa_y[i] - ta2_xz_zz_xzzz_1[i] * pc_y[i];

        ta2_xz_yzz_yyyy_0[i] = 4.0 * ta2_xz_zz_yyy_0[i] * fe_0 - 4.0 * ta2_xz_zz_yyy_1[i] * fe_0 + ta2_xz_zz_yyyy_0[i] * pa_y[i] - ta2_xz_zz_yyyy_1[i] * pc_y[i];

        ta2_xz_yzz_yyyz_0[i] = 3.0 * ta2_xz_zz_yyz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yyz_1[i] * fe_0 + ta2_xz_zz_yyyz_0[i] * pa_y[i] - ta2_xz_zz_yyyz_1[i] * pc_y[i];

        ta2_xz_yzz_yyzz_0[i] = 2.0 * ta2_xz_zz_yzz_0[i] * fe_0 - 2.0 * ta2_xz_zz_yzz_1[i] * fe_0 + ta2_xz_zz_yyzz_0[i] * pa_y[i] - ta2_xz_zz_yyzz_1[i] * pc_y[i];

        ta2_xz_yzz_yzzz_0[i] = ta2_xz_zz_zzz_0[i] * fe_0 - ta2_xz_zz_zzz_1[i] * fe_0 + ta2_xz_zz_yzzz_0[i] * pa_y[i] - ta2_xz_zz_yzzz_1[i] * pc_y[i];

        ta2_xz_yzz_zzzz_0[i] = ta2_xz_zz_zzzz_0[i] * pa_y[i] - ta2_xz_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 435-450 components of targeted buffer : FG

    auto ta2_xz_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 435);

    auto ta2_xz_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 436);

    auto ta2_xz_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 437);

    auto ta2_xz_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 438);

    auto ta2_xz_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 439);

    auto ta2_xz_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 440);

    auto ta2_xz_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 441);

    auto ta2_xz_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 442);

    auto ta2_xz_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 443);

    auto ta2_xz_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 444);

    auto ta2_xz_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 445);

    auto ta2_xz_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 446);

    auto ta2_xz_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 447);

    auto ta2_xz_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 448);

    auto ta2_xz_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 449);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zz_xxxx_1, ta1_x_zz_xxxy_1, ta1_x_zz_xxxz_1, ta1_x_zz_xxyy_1, ta1_x_zz_xxyz_1, ta1_x_zz_xxzz_1, ta1_x_zz_xyyy_1, ta1_x_zz_xyyz_1, ta1_x_zz_xyzz_1, ta1_x_zz_xzzz_1, ta1_x_zz_yyyy_1, ta1_x_zz_yyyz_1, ta1_x_zz_yyzz_1, ta1_x_zz_yzzz_1, ta1_x_zz_zzzz_1, ta2_xz_z_xxxx_0, ta2_xz_z_xxxx_1, ta2_xz_z_xxxy_0, ta2_xz_z_xxxy_1, ta2_xz_z_xxxz_0, ta2_xz_z_xxxz_1, ta2_xz_z_xxyy_0, ta2_xz_z_xxyy_1, ta2_xz_z_xxyz_0, ta2_xz_z_xxyz_1, ta2_xz_z_xxzz_0, ta2_xz_z_xxzz_1, ta2_xz_z_xyyy_0, ta2_xz_z_xyyy_1, ta2_xz_z_xyyz_0, ta2_xz_z_xyyz_1, ta2_xz_z_xyzz_0, ta2_xz_z_xyzz_1, ta2_xz_z_xzzz_0, ta2_xz_z_xzzz_1, ta2_xz_z_yyyy_0, ta2_xz_z_yyyy_1, ta2_xz_z_yyyz_0, ta2_xz_z_yyyz_1, ta2_xz_z_yyzz_0, ta2_xz_z_yyzz_1, ta2_xz_z_yzzz_0, ta2_xz_z_yzzz_1, ta2_xz_z_zzzz_0, ta2_xz_z_zzzz_1, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxxx_0, ta2_xz_zz_xxxx_1, ta2_xz_zz_xxxy_0, ta2_xz_zz_xxxy_1, ta2_xz_zz_xxxz_0, ta2_xz_zz_xxxz_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxyy_0, ta2_xz_zz_xxyy_1, ta2_xz_zz_xxyz_0, ta2_xz_zz_xxyz_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xxzz_0, ta2_xz_zz_xxzz_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyyy_0, ta2_xz_zz_xyyy_1, ta2_xz_zz_xyyz_0, ta2_xz_zz_xyyz_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xyzz_0, ta2_xz_zz_xyzz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_xzzz_0, ta2_xz_zz_xzzz_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyyy_0, ta2_xz_zz_yyyy_1, ta2_xz_zz_yyyz_0, ta2_xz_zz_yyyz_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yyzz_0, ta2_xz_zz_yyzz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_yzzz_0, ta2_xz_zz_yzzz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, ta2_xz_zz_zzzz_0, ta2_xz_zz_zzzz_1, ta2_xz_zzz_xxxx_0, ta2_xz_zzz_xxxy_0, ta2_xz_zzz_xxxz_0, ta2_xz_zzz_xxyy_0, ta2_xz_zzz_xxyz_0, ta2_xz_zzz_xxzz_0, ta2_xz_zzz_xyyy_0, ta2_xz_zzz_xyyz_0, ta2_xz_zzz_xyzz_0, ta2_xz_zzz_xzzz_0, ta2_xz_zzz_yyyy_0, ta2_xz_zzz_yyyz_0, ta2_xz_zzz_yyzz_0, ta2_xz_zzz_yzzz_0, ta2_xz_zzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzz_xxxx_0[i] = 2.0 * ta2_xz_z_xxxx_0[i] * fe_0 - 2.0 * ta2_xz_z_xxxx_1[i] * fe_0 + ta1_x_zz_xxxx_1[i] + ta2_xz_zz_xxxx_0[i] * pa_z[i] - ta2_xz_zz_xxxx_1[i] * pc_z[i];

        ta2_xz_zzz_xxxy_0[i] = 2.0 * ta2_xz_z_xxxy_0[i] * fe_0 - 2.0 * ta2_xz_z_xxxy_1[i] * fe_0 + ta1_x_zz_xxxy_1[i] + ta2_xz_zz_xxxy_0[i] * pa_z[i] - ta2_xz_zz_xxxy_1[i] * pc_z[i];

        ta2_xz_zzz_xxxz_0[i] = 2.0 * ta2_xz_z_xxxz_0[i] * fe_0 - 2.0 * ta2_xz_z_xxxz_1[i] * fe_0 + ta2_xz_zz_xxx_0[i] * fe_0 - ta2_xz_zz_xxx_1[i] * fe_0 + ta1_x_zz_xxxz_1[i] + ta2_xz_zz_xxxz_0[i] * pa_z[i] - ta2_xz_zz_xxxz_1[i] * pc_z[i];

        ta2_xz_zzz_xxyy_0[i] = 2.0 * ta2_xz_z_xxyy_0[i] * fe_0 - 2.0 * ta2_xz_z_xxyy_1[i] * fe_0 + ta1_x_zz_xxyy_1[i] + ta2_xz_zz_xxyy_0[i] * pa_z[i] - ta2_xz_zz_xxyy_1[i] * pc_z[i];

        ta2_xz_zzz_xxyz_0[i] = 2.0 * ta2_xz_z_xxyz_0[i] * fe_0 - 2.0 * ta2_xz_z_xxyz_1[i] * fe_0 + ta2_xz_zz_xxy_0[i] * fe_0 - ta2_xz_zz_xxy_1[i] * fe_0 + ta1_x_zz_xxyz_1[i] + ta2_xz_zz_xxyz_0[i] * pa_z[i] - ta2_xz_zz_xxyz_1[i] * pc_z[i];

        ta2_xz_zzz_xxzz_0[i] = 2.0 * ta2_xz_z_xxzz_0[i] * fe_0 - 2.0 * ta2_xz_z_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_zz_xxz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xxz_1[i] * fe_0 + ta1_x_zz_xxzz_1[i] + ta2_xz_zz_xxzz_0[i] * pa_z[i] - ta2_xz_zz_xxzz_1[i] * pc_z[i];

        ta2_xz_zzz_xyyy_0[i] = 2.0 * ta2_xz_z_xyyy_0[i] * fe_0 - 2.0 * ta2_xz_z_xyyy_1[i] * fe_0 + ta1_x_zz_xyyy_1[i] + ta2_xz_zz_xyyy_0[i] * pa_z[i] - ta2_xz_zz_xyyy_1[i] * pc_z[i];

        ta2_xz_zzz_xyyz_0[i] = 2.0 * ta2_xz_z_xyyz_0[i] * fe_0 - 2.0 * ta2_xz_z_xyyz_1[i] * fe_0 + ta2_xz_zz_xyy_0[i] * fe_0 - ta2_xz_zz_xyy_1[i] * fe_0 + ta1_x_zz_xyyz_1[i] + ta2_xz_zz_xyyz_0[i] * pa_z[i] - ta2_xz_zz_xyyz_1[i] * pc_z[i];

        ta2_xz_zzz_xyzz_0[i] = 2.0 * ta2_xz_z_xyzz_0[i] * fe_0 - 2.0 * ta2_xz_z_xyzz_1[i] * fe_0 + 2.0 * ta2_xz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xyz_1[i] * fe_0 + ta1_x_zz_xyzz_1[i] + ta2_xz_zz_xyzz_0[i] * pa_z[i] - ta2_xz_zz_xyzz_1[i] * pc_z[i];

        ta2_xz_zzz_xzzz_0[i] = 2.0 * ta2_xz_z_xzzz_0[i] * fe_0 - 2.0 * ta2_xz_z_xzzz_1[i] * fe_0 + 3.0 * ta2_xz_zz_xzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xzz_1[i] * fe_0 + ta1_x_zz_xzzz_1[i] + ta2_xz_zz_xzzz_0[i] * pa_z[i] - ta2_xz_zz_xzzz_1[i] * pc_z[i];

        ta2_xz_zzz_yyyy_0[i] = 2.0 * ta2_xz_z_yyyy_0[i] * fe_0 - 2.0 * ta2_xz_z_yyyy_1[i] * fe_0 + ta1_x_zz_yyyy_1[i] + ta2_xz_zz_yyyy_0[i] * pa_z[i] - ta2_xz_zz_yyyy_1[i] * pc_z[i];

        ta2_xz_zzz_yyyz_0[i] = 2.0 * ta2_xz_z_yyyz_0[i] * fe_0 - 2.0 * ta2_xz_z_yyyz_1[i] * fe_0 + ta2_xz_zz_yyy_0[i] * fe_0 - ta2_xz_zz_yyy_1[i] * fe_0 + ta1_x_zz_yyyz_1[i] + ta2_xz_zz_yyyz_0[i] * pa_z[i] - ta2_xz_zz_yyyz_1[i] * pc_z[i];

        ta2_xz_zzz_yyzz_0[i] = 2.0 * ta2_xz_z_yyzz_0[i] * fe_0 - 2.0 * ta2_xz_z_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_zz_yyz_0[i] * fe_0 - 2.0 * ta2_xz_zz_yyz_1[i] * fe_0 + ta1_x_zz_yyzz_1[i] + ta2_xz_zz_yyzz_0[i] * pa_z[i] - ta2_xz_zz_yyzz_1[i] * pc_z[i];

        ta2_xz_zzz_yzzz_0[i] = 2.0 * ta2_xz_z_yzzz_0[i] * fe_0 - 2.0 * ta2_xz_z_yzzz_1[i] * fe_0 + 3.0 * ta2_xz_zz_yzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yzz_1[i] * fe_0 + ta1_x_zz_yzzz_1[i] + ta2_xz_zz_yzzz_0[i] * pa_z[i] - ta2_xz_zz_yzzz_1[i] * pc_z[i];

        ta2_xz_zzz_zzzz_0[i] = 2.0 * ta2_xz_z_zzzz_0[i] * fe_0 - 2.0 * ta2_xz_z_zzzz_1[i] * fe_0 + 4.0 * ta2_xz_zz_zzz_0[i] * fe_0 - 4.0 * ta2_xz_zz_zzz_1[i] * fe_0 + ta1_x_zz_zzzz_1[i] + ta2_xz_zz_zzzz_0[i] * pa_z[i] - ta2_xz_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 450-465 components of targeted buffer : FG

    auto ta2_yy_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 450);

    auto ta2_yy_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 451);

    auto ta2_yy_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 452);

    auto ta2_yy_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 453);

    auto ta2_yy_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 454);

    auto ta2_yy_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 455);

    auto ta2_yy_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 456);

    auto ta2_yy_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 457);

    auto ta2_yy_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 458);

    auto ta2_yy_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 459);

    auto ta2_yy_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 460);

    auto ta2_yy_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 461);

    auto ta2_yy_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 462);

    auto ta2_yy_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 463);

    auto ta2_yy_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 464);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_x_xxxx_0, ta2_yy_x_xxxx_1, ta2_yy_x_xxxy_0, ta2_yy_x_xxxy_1, ta2_yy_x_xxxz_0, ta2_yy_x_xxxz_1, ta2_yy_x_xxyy_0, ta2_yy_x_xxyy_1, ta2_yy_x_xxyz_0, ta2_yy_x_xxyz_1, ta2_yy_x_xxzz_0, ta2_yy_x_xxzz_1, ta2_yy_x_xyyy_0, ta2_yy_x_xyyy_1, ta2_yy_x_xyyz_0, ta2_yy_x_xyyz_1, ta2_yy_x_xyzz_0, ta2_yy_x_xyzz_1, ta2_yy_x_xzzz_0, ta2_yy_x_xzzz_1, ta2_yy_x_yyyy_0, ta2_yy_x_yyyy_1, ta2_yy_x_yyyz_0, ta2_yy_x_yyyz_1, ta2_yy_x_yyzz_0, ta2_yy_x_yyzz_1, ta2_yy_x_yzzz_0, ta2_yy_x_yzzz_1, ta2_yy_x_zzzz_0, ta2_yy_x_zzzz_1, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxxx_0, ta2_yy_xx_xxxx_1, ta2_yy_xx_xxxy_0, ta2_yy_xx_xxxy_1, ta2_yy_xx_xxxz_0, ta2_yy_xx_xxxz_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxyy_0, ta2_yy_xx_xxyy_1, ta2_yy_xx_xxyz_0, ta2_yy_xx_xxyz_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xxzz_0, ta2_yy_xx_xxzz_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyyy_0, ta2_yy_xx_xyyy_1, ta2_yy_xx_xyyz_0, ta2_yy_xx_xyyz_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xyzz_0, ta2_yy_xx_xyzz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_xzzz_0, ta2_yy_xx_xzzz_1, ta2_yy_xx_yyy_0, ta2_yy_xx_yyy_1, ta2_yy_xx_yyyy_0, ta2_yy_xx_yyyy_1, ta2_yy_xx_yyyz_0, ta2_yy_xx_yyyz_1, ta2_yy_xx_yyz_0, ta2_yy_xx_yyz_1, ta2_yy_xx_yyzz_0, ta2_yy_xx_yyzz_1, ta2_yy_xx_yzz_0, ta2_yy_xx_yzz_1, ta2_yy_xx_yzzz_0, ta2_yy_xx_yzzz_1, ta2_yy_xx_zzz_0, ta2_yy_xx_zzz_1, ta2_yy_xx_zzzz_0, ta2_yy_xx_zzzz_1, ta2_yy_xxx_xxxx_0, ta2_yy_xxx_xxxy_0, ta2_yy_xxx_xxxz_0, ta2_yy_xxx_xxyy_0, ta2_yy_xxx_xxyz_0, ta2_yy_xxx_xxzz_0, ta2_yy_xxx_xyyy_0, ta2_yy_xxx_xyyz_0, ta2_yy_xxx_xyzz_0, ta2_yy_xxx_xzzz_0, ta2_yy_xxx_yyyy_0, ta2_yy_xxx_yyyz_0, ta2_yy_xxx_yyzz_0, ta2_yy_xxx_yzzz_0, ta2_yy_xxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxx_xxxx_0[i] = 2.0 * ta2_yy_x_xxxx_0[i] * fe_0 - 2.0 * ta2_yy_x_xxxx_1[i] * fe_0 + 4.0 * ta2_yy_xx_xxx_0[i] * fe_0 - 4.0 * ta2_yy_xx_xxx_1[i] * fe_0 + ta2_yy_xx_xxxx_0[i] * pa_x[i] - ta2_yy_xx_xxxx_1[i] * pc_x[i];

        ta2_yy_xxx_xxxy_0[i] = 2.0 * ta2_yy_x_xxxy_0[i] * fe_0 - 2.0 * ta2_yy_x_xxxy_1[i] * fe_0 + 3.0 * ta2_yy_xx_xxy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxy_1[i] * fe_0 + ta2_yy_xx_xxxy_0[i] * pa_x[i] - ta2_yy_xx_xxxy_1[i] * pc_x[i];

        ta2_yy_xxx_xxxz_0[i] = 2.0 * ta2_yy_x_xxxz_0[i] * fe_0 - 2.0 * ta2_yy_x_xxxz_1[i] * fe_0 + 3.0 * ta2_yy_xx_xxz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxz_1[i] * fe_0 + ta2_yy_xx_xxxz_0[i] * pa_x[i] - ta2_yy_xx_xxxz_1[i] * pc_x[i];

        ta2_yy_xxx_xxyy_0[i] = 2.0 * ta2_yy_x_xxyy_0[i] * fe_0 - 2.0 * ta2_yy_x_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_xx_xyy_0[i] * fe_0 - 2.0 * ta2_yy_xx_xyy_1[i] * fe_0 + ta2_yy_xx_xxyy_0[i] * pa_x[i] - ta2_yy_xx_xxyy_1[i] * pc_x[i];

        ta2_yy_xxx_xxyz_0[i] = 2.0 * ta2_yy_x_xxyz_0[i] * fe_0 - 2.0 * ta2_yy_x_xxyz_1[i] * fe_0 + 2.0 * ta2_yy_xx_xyz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xyz_1[i] * fe_0 + ta2_yy_xx_xxyz_0[i] * pa_x[i] - ta2_yy_xx_xxyz_1[i] * pc_x[i];

        ta2_yy_xxx_xxzz_0[i] = 2.0 * ta2_yy_x_xxzz_0[i] * fe_0 - 2.0 * ta2_yy_x_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_xx_xzz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xzz_1[i] * fe_0 + ta2_yy_xx_xxzz_0[i] * pa_x[i] - ta2_yy_xx_xxzz_1[i] * pc_x[i];

        ta2_yy_xxx_xyyy_0[i] = 2.0 * ta2_yy_x_xyyy_0[i] * fe_0 - 2.0 * ta2_yy_x_xyyy_1[i] * fe_0 + ta2_yy_xx_yyy_0[i] * fe_0 - ta2_yy_xx_yyy_1[i] * fe_0 + ta2_yy_xx_xyyy_0[i] * pa_x[i] - ta2_yy_xx_xyyy_1[i] * pc_x[i];

        ta2_yy_xxx_xyyz_0[i] = 2.0 * ta2_yy_x_xyyz_0[i] * fe_0 - 2.0 * ta2_yy_x_xyyz_1[i] * fe_0 + ta2_yy_xx_yyz_0[i] * fe_0 - ta2_yy_xx_yyz_1[i] * fe_0 + ta2_yy_xx_xyyz_0[i] * pa_x[i] - ta2_yy_xx_xyyz_1[i] * pc_x[i];

        ta2_yy_xxx_xyzz_0[i] = 2.0 * ta2_yy_x_xyzz_0[i] * fe_0 - 2.0 * ta2_yy_x_xyzz_1[i] * fe_0 + ta2_yy_xx_yzz_0[i] * fe_0 - ta2_yy_xx_yzz_1[i] * fe_0 + ta2_yy_xx_xyzz_0[i] * pa_x[i] - ta2_yy_xx_xyzz_1[i] * pc_x[i];

        ta2_yy_xxx_xzzz_0[i] = 2.0 * ta2_yy_x_xzzz_0[i] * fe_0 - 2.0 * ta2_yy_x_xzzz_1[i] * fe_0 + ta2_yy_xx_zzz_0[i] * fe_0 - ta2_yy_xx_zzz_1[i] * fe_0 + ta2_yy_xx_xzzz_0[i] * pa_x[i] - ta2_yy_xx_xzzz_1[i] * pc_x[i];

        ta2_yy_xxx_yyyy_0[i] = 2.0 * ta2_yy_x_yyyy_0[i] * fe_0 - 2.0 * ta2_yy_x_yyyy_1[i] * fe_0 + ta2_yy_xx_yyyy_0[i] * pa_x[i] - ta2_yy_xx_yyyy_1[i] * pc_x[i];

        ta2_yy_xxx_yyyz_0[i] = 2.0 * ta2_yy_x_yyyz_0[i] * fe_0 - 2.0 * ta2_yy_x_yyyz_1[i] * fe_0 + ta2_yy_xx_yyyz_0[i] * pa_x[i] - ta2_yy_xx_yyyz_1[i] * pc_x[i];

        ta2_yy_xxx_yyzz_0[i] = 2.0 * ta2_yy_x_yyzz_0[i] * fe_0 - 2.0 * ta2_yy_x_yyzz_1[i] * fe_0 + ta2_yy_xx_yyzz_0[i] * pa_x[i] - ta2_yy_xx_yyzz_1[i] * pc_x[i];

        ta2_yy_xxx_yzzz_0[i] = 2.0 * ta2_yy_x_yzzz_0[i] * fe_0 - 2.0 * ta2_yy_x_yzzz_1[i] * fe_0 + ta2_yy_xx_yzzz_0[i] * pa_x[i] - ta2_yy_xx_yzzz_1[i] * pc_x[i];

        ta2_yy_xxx_zzzz_0[i] = 2.0 * ta2_yy_x_zzzz_0[i] * fe_0 - 2.0 * ta2_yy_x_zzzz_1[i] * fe_0 + ta2_yy_xx_zzzz_0[i] * pa_x[i] - ta2_yy_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 465-480 components of targeted buffer : FG

    auto ta2_yy_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 465);

    auto ta2_yy_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 466);

    auto ta2_yy_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 467);

    auto ta2_yy_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 468);

    auto ta2_yy_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 469);

    auto ta2_yy_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 470);

    auto ta2_yy_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 471);

    auto ta2_yy_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 472);

    auto ta2_yy_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 473);

    auto ta2_yy_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 474);

    auto ta2_yy_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 475);

    auto ta2_yy_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 476);

    auto ta2_yy_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 477);

    auto ta2_yy_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 478);

    auto ta2_yy_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 479);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xx_xxxx_1, ta1_y_xx_xxxy_1, ta1_y_xx_xxxz_1, ta1_y_xx_xxyy_1, ta1_y_xx_xxyz_1, ta1_y_xx_xxzz_1, ta1_y_xx_xyyy_1, ta1_y_xx_xyyz_1, ta1_y_xx_xyzz_1, ta1_y_xx_xzzz_1, ta1_y_xx_zzzz_1, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxxx_0, ta2_yy_xx_xxxx_1, ta2_yy_xx_xxxy_0, ta2_yy_xx_xxxy_1, ta2_yy_xx_xxxz_0, ta2_yy_xx_xxxz_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxyy_0, ta2_yy_xx_xxyy_1, ta2_yy_xx_xxyz_0, ta2_yy_xx_xxyz_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xxzz_0, ta2_yy_xx_xxzz_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyyy_0, ta2_yy_xx_xyyy_1, ta2_yy_xx_xyyz_0, ta2_yy_xx_xyyz_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xyzz_0, ta2_yy_xx_xyzz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_xzzz_0, ta2_yy_xx_xzzz_1, ta2_yy_xx_zzzz_0, ta2_yy_xx_zzzz_1, ta2_yy_xxy_xxxx_0, ta2_yy_xxy_xxxy_0, ta2_yy_xxy_xxxz_0, ta2_yy_xxy_xxyy_0, ta2_yy_xxy_xxyz_0, ta2_yy_xxy_xxzz_0, ta2_yy_xxy_xyyy_0, ta2_yy_xxy_xyyz_0, ta2_yy_xxy_xyzz_0, ta2_yy_xxy_xzzz_0, ta2_yy_xxy_yyyy_0, ta2_yy_xxy_yyyz_0, ta2_yy_xxy_yyzz_0, ta2_yy_xxy_yzzz_0, ta2_yy_xxy_zzzz_0, ta2_yy_xy_yyyy_0, ta2_yy_xy_yyyy_1, ta2_yy_xy_yyyz_0, ta2_yy_xy_yyyz_1, ta2_yy_xy_yyzz_0, ta2_yy_xy_yyzz_1, ta2_yy_xy_yzzz_0, ta2_yy_xy_yzzz_1, ta2_yy_y_yyyy_0, ta2_yy_y_yyyy_1, ta2_yy_y_yyyz_0, ta2_yy_y_yyyz_1, ta2_yy_y_yyzz_0, ta2_yy_y_yyzz_1, ta2_yy_y_yzzz_0, ta2_yy_y_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxy_xxxx_0[i] = 2.0 * ta1_y_xx_xxxx_1[i] + ta2_yy_xx_xxxx_0[i] * pa_y[i] - ta2_yy_xx_xxxx_1[i] * pc_y[i];

        ta2_yy_xxy_xxxy_0[i] = ta2_yy_xx_xxx_0[i] * fe_0 - ta2_yy_xx_xxx_1[i] * fe_0 + 2.0 * ta1_y_xx_xxxy_1[i] + ta2_yy_xx_xxxy_0[i] * pa_y[i] - ta2_yy_xx_xxxy_1[i] * pc_y[i];

        ta2_yy_xxy_xxxz_0[i] = 2.0 * ta1_y_xx_xxxz_1[i] + ta2_yy_xx_xxxz_0[i] * pa_y[i] - ta2_yy_xx_xxxz_1[i] * pc_y[i];

        ta2_yy_xxy_xxyy_0[i] = 2.0 * ta2_yy_xx_xxy_0[i] * fe_0 - 2.0 * ta2_yy_xx_xxy_1[i] * fe_0 + 2.0 * ta1_y_xx_xxyy_1[i] + ta2_yy_xx_xxyy_0[i] * pa_y[i] - ta2_yy_xx_xxyy_1[i] * pc_y[i];

        ta2_yy_xxy_xxyz_0[i] = ta2_yy_xx_xxz_0[i] * fe_0 - ta2_yy_xx_xxz_1[i] * fe_0 + 2.0 * ta1_y_xx_xxyz_1[i] + ta2_yy_xx_xxyz_0[i] * pa_y[i] - ta2_yy_xx_xxyz_1[i] * pc_y[i];

        ta2_yy_xxy_xxzz_0[i] = 2.0 * ta1_y_xx_xxzz_1[i] + ta2_yy_xx_xxzz_0[i] * pa_y[i] - ta2_yy_xx_xxzz_1[i] * pc_y[i];

        ta2_yy_xxy_xyyy_0[i] = 3.0 * ta2_yy_xx_xyy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xyy_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyy_1[i] + ta2_yy_xx_xyyy_0[i] * pa_y[i] - ta2_yy_xx_xyyy_1[i] * pc_y[i];

        ta2_yy_xxy_xyyz_0[i] = 2.0 * ta2_yy_xx_xyz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xyz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyz_1[i] + ta2_yy_xx_xyyz_0[i] * pa_y[i] - ta2_yy_xx_xyyz_1[i] * pc_y[i];

        ta2_yy_xxy_xyzz_0[i] = ta2_yy_xx_xzz_0[i] * fe_0 - ta2_yy_xx_xzz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyzz_1[i] + ta2_yy_xx_xyzz_0[i] * pa_y[i] - ta2_yy_xx_xyzz_1[i] * pc_y[i];

        ta2_yy_xxy_xzzz_0[i] = 2.0 * ta1_y_xx_xzzz_1[i] + ta2_yy_xx_xzzz_0[i] * pa_y[i] - ta2_yy_xx_xzzz_1[i] * pc_y[i];

        ta2_yy_xxy_yyyy_0[i] = ta2_yy_y_yyyy_0[i] * fe_0 - ta2_yy_y_yyyy_1[i] * fe_0 + ta2_yy_xy_yyyy_0[i] * pa_x[i] - ta2_yy_xy_yyyy_1[i] * pc_x[i];

        ta2_yy_xxy_yyyz_0[i] = ta2_yy_y_yyyz_0[i] * fe_0 - ta2_yy_y_yyyz_1[i] * fe_0 + ta2_yy_xy_yyyz_0[i] * pa_x[i] - ta2_yy_xy_yyyz_1[i] * pc_x[i];

        ta2_yy_xxy_yyzz_0[i] = ta2_yy_y_yyzz_0[i] * fe_0 - ta2_yy_y_yyzz_1[i] * fe_0 + ta2_yy_xy_yyzz_0[i] * pa_x[i] - ta2_yy_xy_yyzz_1[i] * pc_x[i];

        ta2_yy_xxy_yzzz_0[i] = ta2_yy_y_yzzz_0[i] * fe_0 - ta2_yy_y_yzzz_1[i] * fe_0 + ta2_yy_xy_yzzz_0[i] * pa_x[i] - ta2_yy_xy_yzzz_1[i] * pc_x[i];

        ta2_yy_xxy_zzzz_0[i] = 2.0 * ta1_y_xx_zzzz_1[i] + ta2_yy_xx_zzzz_0[i] * pa_y[i] - ta2_yy_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 480-495 components of targeted buffer : FG

    auto ta2_yy_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 480);

    auto ta2_yy_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 481);

    auto ta2_yy_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 482);

    auto ta2_yy_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 483);

    auto ta2_yy_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 484);

    auto ta2_yy_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 485);

    auto ta2_yy_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 486);

    auto ta2_yy_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 487);

    auto ta2_yy_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 488);

    auto ta2_yy_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 489);

    auto ta2_yy_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 490);

    auto ta2_yy_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 491);

    auto ta2_yy_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 492);

    auto ta2_yy_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 493);

    auto ta2_yy_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 494);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxxx_0, ta2_yy_xx_xxxx_1, ta2_yy_xx_xxxy_0, ta2_yy_xx_xxxy_1, ta2_yy_xx_xxxz_0, ta2_yy_xx_xxxz_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxyy_0, ta2_yy_xx_xxyy_1, ta2_yy_xx_xxyz_0, ta2_yy_xx_xxyz_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xxzz_0, ta2_yy_xx_xxzz_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyyy_0, ta2_yy_xx_xyyy_1, ta2_yy_xx_xyyz_0, ta2_yy_xx_xyyz_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xyzz_0, ta2_yy_xx_xyzz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_xzzz_0, ta2_yy_xx_xzzz_1, ta2_yy_xx_yyyy_0, ta2_yy_xx_yyyy_1, ta2_yy_xxz_xxxx_0, ta2_yy_xxz_xxxy_0, ta2_yy_xxz_xxxz_0, ta2_yy_xxz_xxyy_0, ta2_yy_xxz_xxyz_0, ta2_yy_xxz_xxzz_0, ta2_yy_xxz_xyyy_0, ta2_yy_xxz_xyyz_0, ta2_yy_xxz_xyzz_0, ta2_yy_xxz_xzzz_0, ta2_yy_xxz_yyyy_0, ta2_yy_xxz_yyyz_0, ta2_yy_xxz_yyzz_0, ta2_yy_xxz_yzzz_0, ta2_yy_xxz_zzzz_0, ta2_yy_xz_yyyz_0, ta2_yy_xz_yyyz_1, ta2_yy_xz_yyzz_0, ta2_yy_xz_yyzz_1, ta2_yy_xz_yzzz_0, ta2_yy_xz_yzzz_1, ta2_yy_xz_zzzz_0, ta2_yy_xz_zzzz_1, ta2_yy_z_yyyz_0, ta2_yy_z_yyyz_1, ta2_yy_z_yyzz_0, ta2_yy_z_yyzz_1, ta2_yy_z_yzzz_0, ta2_yy_z_yzzz_1, ta2_yy_z_zzzz_0, ta2_yy_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxz_xxxx_0[i] = ta2_yy_xx_xxxx_0[i] * pa_z[i] - ta2_yy_xx_xxxx_1[i] * pc_z[i];

        ta2_yy_xxz_xxxy_0[i] = ta2_yy_xx_xxxy_0[i] * pa_z[i] - ta2_yy_xx_xxxy_1[i] * pc_z[i];

        ta2_yy_xxz_xxxz_0[i] = ta2_yy_xx_xxx_0[i] * fe_0 - ta2_yy_xx_xxx_1[i] * fe_0 + ta2_yy_xx_xxxz_0[i] * pa_z[i] - ta2_yy_xx_xxxz_1[i] * pc_z[i];

        ta2_yy_xxz_xxyy_0[i] = ta2_yy_xx_xxyy_0[i] * pa_z[i] - ta2_yy_xx_xxyy_1[i] * pc_z[i];

        ta2_yy_xxz_xxyz_0[i] = ta2_yy_xx_xxy_0[i] * fe_0 - ta2_yy_xx_xxy_1[i] * fe_0 + ta2_yy_xx_xxyz_0[i] * pa_z[i] - ta2_yy_xx_xxyz_1[i] * pc_z[i];

        ta2_yy_xxz_xxzz_0[i] = 2.0 * ta2_yy_xx_xxz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xxz_1[i] * fe_0 + ta2_yy_xx_xxzz_0[i] * pa_z[i] - ta2_yy_xx_xxzz_1[i] * pc_z[i];

        ta2_yy_xxz_xyyy_0[i] = ta2_yy_xx_xyyy_0[i] * pa_z[i] - ta2_yy_xx_xyyy_1[i] * pc_z[i];

        ta2_yy_xxz_xyyz_0[i] = ta2_yy_xx_xyy_0[i] * fe_0 - ta2_yy_xx_xyy_1[i] * fe_0 + ta2_yy_xx_xyyz_0[i] * pa_z[i] - ta2_yy_xx_xyyz_1[i] * pc_z[i];

        ta2_yy_xxz_xyzz_0[i] = 2.0 * ta2_yy_xx_xyz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xyz_1[i] * fe_0 + ta2_yy_xx_xyzz_0[i] * pa_z[i] - ta2_yy_xx_xyzz_1[i] * pc_z[i];

        ta2_yy_xxz_xzzz_0[i] = 3.0 * ta2_yy_xx_xzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xzz_1[i] * fe_0 + ta2_yy_xx_xzzz_0[i] * pa_z[i] - ta2_yy_xx_xzzz_1[i] * pc_z[i];

        ta2_yy_xxz_yyyy_0[i] = ta2_yy_xx_yyyy_0[i] * pa_z[i] - ta2_yy_xx_yyyy_1[i] * pc_z[i];

        ta2_yy_xxz_yyyz_0[i] = ta2_yy_z_yyyz_0[i] * fe_0 - ta2_yy_z_yyyz_1[i] * fe_0 + ta2_yy_xz_yyyz_0[i] * pa_x[i] - ta2_yy_xz_yyyz_1[i] * pc_x[i];

        ta2_yy_xxz_yyzz_0[i] = ta2_yy_z_yyzz_0[i] * fe_0 - ta2_yy_z_yyzz_1[i] * fe_0 + ta2_yy_xz_yyzz_0[i] * pa_x[i] - ta2_yy_xz_yyzz_1[i] * pc_x[i];

        ta2_yy_xxz_yzzz_0[i] = ta2_yy_z_yzzz_0[i] * fe_0 - ta2_yy_z_yzzz_1[i] * fe_0 + ta2_yy_xz_yzzz_0[i] * pa_x[i] - ta2_yy_xz_yzzz_1[i] * pc_x[i];

        ta2_yy_xxz_zzzz_0[i] = ta2_yy_z_zzzz_0[i] * fe_0 - ta2_yy_z_zzzz_1[i] * fe_0 + ta2_yy_xz_zzzz_0[i] * pa_x[i] - ta2_yy_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 495-510 components of targeted buffer : FG

    auto ta2_yy_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 495);

    auto ta2_yy_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 496);

    auto ta2_yy_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 497);

    auto ta2_yy_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 498);

    auto ta2_yy_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 499);

    auto ta2_yy_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 500);

    auto ta2_yy_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 501);

    auto ta2_yy_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 502);

    auto ta2_yy_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 503);

    auto ta2_yy_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 504);

    auto ta2_yy_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 505);

    auto ta2_yy_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 506);

    auto ta2_yy_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 507);

    auto ta2_yy_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 508);

    auto ta2_yy_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 509);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xyy_xxxx_0, ta2_yy_xyy_xxxy_0, ta2_yy_xyy_xxxz_0, ta2_yy_xyy_xxyy_0, ta2_yy_xyy_xxyz_0, ta2_yy_xyy_xxzz_0, ta2_yy_xyy_xyyy_0, ta2_yy_xyy_xyyz_0, ta2_yy_xyy_xyzz_0, ta2_yy_xyy_xzzz_0, ta2_yy_xyy_yyyy_0, ta2_yy_xyy_yyyz_0, ta2_yy_xyy_yyzz_0, ta2_yy_xyy_yzzz_0, ta2_yy_xyy_zzzz_0, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxxx_0, ta2_yy_yy_xxxx_1, ta2_yy_yy_xxxy_0, ta2_yy_yy_xxxy_1, ta2_yy_yy_xxxz_0, ta2_yy_yy_xxxz_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxyy_0, ta2_yy_yy_xxyy_1, ta2_yy_yy_xxyz_0, ta2_yy_yy_xxyz_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xxzz_0, ta2_yy_yy_xxzz_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyyy_0, ta2_yy_yy_xyyy_1, ta2_yy_yy_xyyz_0, ta2_yy_yy_xyyz_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xyzz_0, ta2_yy_yy_xyzz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_xzzz_0, ta2_yy_yy_xzzz_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyyy_0, ta2_yy_yy_yyyy_1, ta2_yy_yy_yyyz_0, ta2_yy_yy_yyyz_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yyzz_0, ta2_yy_yy_yyzz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_yzzz_0, ta2_yy_yy_yzzz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, ta2_yy_yy_zzzz_0, ta2_yy_yy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyy_xxxx_0[i] = 4.0 * ta2_yy_yy_xxx_0[i] * fe_0 - 4.0 * ta2_yy_yy_xxx_1[i] * fe_0 + ta2_yy_yy_xxxx_0[i] * pa_x[i] - ta2_yy_yy_xxxx_1[i] * pc_x[i];

        ta2_yy_xyy_xxxy_0[i] = 3.0 * ta2_yy_yy_xxy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxy_1[i] * fe_0 + ta2_yy_yy_xxxy_0[i] * pa_x[i] - ta2_yy_yy_xxxy_1[i] * pc_x[i];

        ta2_yy_xyy_xxxz_0[i] = 3.0 * ta2_yy_yy_xxz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxz_1[i] * fe_0 + ta2_yy_yy_xxxz_0[i] * pa_x[i] - ta2_yy_yy_xxxz_1[i] * pc_x[i];

        ta2_yy_xyy_xxyy_0[i] = 2.0 * ta2_yy_yy_xyy_0[i] * fe_0 - 2.0 * ta2_yy_yy_xyy_1[i] * fe_0 + ta2_yy_yy_xxyy_0[i] * pa_x[i] - ta2_yy_yy_xxyy_1[i] * pc_x[i];

        ta2_yy_xyy_xxyz_0[i] = 2.0 * ta2_yy_yy_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xyz_1[i] * fe_0 + ta2_yy_yy_xxyz_0[i] * pa_x[i] - ta2_yy_yy_xxyz_1[i] * pc_x[i];

        ta2_yy_xyy_xxzz_0[i] = 2.0 * ta2_yy_yy_xzz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xzz_1[i] * fe_0 + ta2_yy_yy_xxzz_0[i] * pa_x[i] - ta2_yy_yy_xxzz_1[i] * pc_x[i];

        ta2_yy_xyy_xyyy_0[i] = ta2_yy_yy_yyy_0[i] * fe_0 - ta2_yy_yy_yyy_1[i] * fe_0 + ta2_yy_yy_xyyy_0[i] * pa_x[i] - ta2_yy_yy_xyyy_1[i] * pc_x[i];

        ta2_yy_xyy_xyyz_0[i] = ta2_yy_yy_yyz_0[i] * fe_0 - ta2_yy_yy_yyz_1[i] * fe_0 + ta2_yy_yy_xyyz_0[i] * pa_x[i] - ta2_yy_yy_xyyz_1[i] * pc_x[i];

        ta2_yy_xyy_xyzz_0[i] = ta2_yy_yy_yzz_0[i] * fe_0 - ta2_yy_yy_yzz_1[i] * fe_0 + ta2_yy_yy_xyzz_0[i] * pa_x[i] - ta2_yy_yy_xyzz_1[i] * pc_x[i];

        ta2_yy_xyy_xzzz_0[i] = ta2_yy_yy_zzz_0[i] * fe_0 - ta2_yy_yy_zzz_1[i] * fe_0 + ta2_yy_yy_xzzz_0[i] * pa_x[i] - ta2_yy_yy_xzzz_1[i] * pc_x[i];

        ta2_yy_xyy_yyyy_0[i] = ta2_yy_yy_yyyy_0[i] * pa_x[i] - ta2_yy_yy_yyyy_1[i] * pc_x[i];

        ta2_yy_xyy_yyyz_0[i] = ta2_yy_yy_yyyz_0[i] * pa_x[i] - ta2_yy_yy_yyyz_1[i] * pc_x[i];

        ta2_yy_xyy_yyzz_0[i] = ta2_yy_yy_yyzz_0[i] * pa_x[i] - ta2_yy_yy_yyzz_1[i] * pc_x[i];

        ta2_yy_xyy_yzzz_0[i] = ta2_yy_yy_yzzz_0[i] * pa_x[i] - ta2_yy_yy_yzzz_1[i] * pc_x[i];

        ta2_yy_xyy_zzzz_0[i] = ta2_yy_yy_zzzz_0[i] * pa_x[i] - ta2_yy_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 510-525 components of targeted buffer : FG

    auto ta2_yy_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 510);

    auto ta2_yy_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 511);

    auto ta2_yy_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 512);

    auto ta2_yy_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 513);

    auto ta2_yy_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 514);

    auto ta2_yy_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 515);

    auto ta2_yy_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 516);

    auto ta2_yy_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 517);

    auto ta2_yy_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 518);

    auto ta2_yy_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 519);

    auto ta2_yy_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 520);

    auto ta2_yy_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 521);

    auto ta2_yy_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 522);

    auto ta2_yy_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 523);

    auto ta2_yy_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 524);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xz_xxxz_1, ta1_y_xz_xxzz_1, ta1_y_xz_xzzz_1, ta2_yy_xy_xxxx_0, ta2_yy_xy_xxxx_1, ta2_yy_xy_xxxy_0, ta2_yy_xy_xxxy_1, ta2_yy_xy_xxyy_0, ta2_yy_xy_xxyy_1, ta2_yy_xy_xyyy_0, ta2_yy_xy_xyyy_1, ta2_yy_xyz_xxxx_0, ta2_yy_xyz_xxxy_0, ta2_yy_xyz_xxxz_0, ta2_yy_xyz_xxyy_0, ta2_yy_xyz_xxyz_0, ta2_yy_xyz_xxzz_0, ta2_yy_xyz_xyyy_0, ta2_yy_xyz_xyyz_0, ta2_yy_xyz_xyzz_0, ta2_yy_xyz_xzzz_0, ta2_yy_xyz_yyyy_0, ta2_yy_xyz_yyyz_0, ta2_yy_xyz_yyzz_0, ta2_yy_xyz_yzzz_0, ta2_yy_xyz_zzzz_0, ta2_yy_xz_xxxz_0, ta2_yy_xz_xxxz_1, ta2_yy_xz_xxzz_0, ta2_yy_xz_xxzz_1, ta2_yy_xz_xzzz_0, ta2_yy_xz_xzzz_1, ta2_yy_yz_xxyz_0, ta2_yy_yz_xxyz_1, ta2_yy_yz_xyyz_0, ta2_yy_yz_xyyz_1, ta2_yy_yz_xyz_0, ta2_yy_yz_xyz_1, ta2_yy_yz_xyzz_0, ta2_yy_yz_xyzz_1, ta2_yy_yz_yyyy_0, ta2_yy_yz_yyyy_1, ta2_yy_yz_yyyz_0, ta2_yy_yz_yyyz_1, ta2_yy_yz_yyz_0, ta2_yy_yz_yyz_1, ta2_yy_yz_yyzz_0, ta2_yy_yz_yyzz_1, ta2_yy_yz_yzz_0, ta2_yy_yz_yzz_1, ta2_yy_yz_yzzz_0, ta2_yy_yz_yzzz_1, ta2_yy_yz_zzzz_0, ta2_yy_yz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyz_xxxx_0[i] = ta2_yy_xy_xxxx_0[i] * pa_z[i] - ta2_yy_xy_xxxx_1[i] * pc_z[i];

        ta2_yy_xyz_xxxy_0[i] = ta2_yy_xy_xxxy_0[i] * pa_z[i] - ta2_yy_xy_xxxy_1[i] * pc_z[i];

        ta2_yy_xyz_xxxz_0[i] = 2.0 * ta1_y_xz_xxxz_1[i] + ta2_yy_xz_xxxz_0[i] * pa_y[i] - ta2_yy_xz_xxxz_1[i] * pc_y[i];

        ta2_yy_xyz_xxyy_0[i] = ta2_yy_xy_xxyy_0[i] * pa_z[i] - ta2_yy_xy_xxyy_1[i] * pc_z[i];

        ta2_yy_xyz_xxyz_0[i] = 2.0 * ta2_yy_yz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yz_xyz_1[i] * fe_0 + ta2_yy_yz_xxyz_0[i] * pa_x[i] - ta2_yy_yz_xxyz_1[i] * pc_x[i];

        ta2_yy_xyz_xxzz_0[i] = 2.0 * ta1_y_xz_xxzz_1[i] + ta2_yy_xz_xxzz_0[i] * pa_y[i] - ta2_yy_xz_xxzz_1[i] * pc_y[i];

        ta2_yy_xyz_xyyy_0[i] = ta2_yy_xy_xyyy_0[i] * pa_z[i] - ta2_yy_xy_xyyy_1[i] * pc_z[i];

        ta2_yy_xyz_xyyz_0[i] = ta2_yy_yz_yyz_0[i] * fe_0 - ta2_yy_yz_yyz_1[i] * fe_0 + ta2_yy_yz_xyyz_0[i] * pa_x[i] - ta2_yy_yz_xyyz_1[i] * pc_x[i];

        ta2_yy_xyz_xyzz_0[i] = ta2_yy_yz_yzz_0[i] * fe_0 - ta2_yy_yz_yzz_1[i] * fe_0 + ta2_yy_yz_xyzz_0[i] * pa_x[i] - ta2_yy_yz_xyzz_1[i] * pc_x[i];

        ta2_yy_xyz_xzzz_0[i] = 2.0 * ta1_y_xz_xzzz_1[i] + ta2_yy_xz_xzzz_0[i] * pa_y[i] - ta2_yy_xz_xzzz_1[i] * pc_y[i];

        ta2_yy_xyz_yyyy_0[i] = ta2_yy_yz_yyyy_0[i] * pa_x[i] - ta2_yy_yz_yyyy_1[i] * pc_x[i];

        ta2_yy_xyz_yyyz_0[i] = ta2_yy_yz_yyyz_0[i] * pa_x[i] - ta2_yy_yz_yyyz_1[i] * pc_x[i];

        ta2_yy_xyz_yyzz_0[i] = ta2_yy_yz_yyzz_0[i] * pa_x[i] - ta2_yy_yz_yyzz_1[i] * pc_x[i];

        ta2_yy_xyz_yzzz_0[i] = ta2_yy_yz_yzzz_0[i] * pa_x[i] - ta2_yy_yz_yzzz_1[i] * pc_x[i];

        ta2_yy_xyz_zzzz_0[i] = ta2_yy_yz_zzzz_0[i] * pa_x[i] - ta2_yy_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 525-540 components of targeted buffer : FG

    auto ta2_yy_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 525);

    auto ta2_yy_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 526);

    auto ta2_yy_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 527);

    auto ta2_yy_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 528);

    auto ta2_yy_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 529);

    auto ta2_yy_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 530);

    auto ta2_yy_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 531);

    auto ta2_yy_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 532);

    auto ta2_yy_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 533);

    auto ta2_yy_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 534);

    auto ta2_yy_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 535);

    auto ta2_yy_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 536);

    auto ta2_yy_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 537);

    auto ta2_yy_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 538);

    auto ta2_yy_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 539);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xzz_xxxx_0, ta2_yy_xzz_xxxy_0, ta2_yy_xzz_xxxz_0, ta2_yy_xzz_xxyy_0, ta2_yy_xzz_xxyz_0, ta2_yy_xzz_xxzz_0, ta2_yy_xzz_xyyy_0, ta2_yy_xzz_xyyz_0, ta2_yy_xzz_xyzz_0, ta2_yy_xzz_xzzz_0, ta2_yy_xzz_yyyy_0, ta2_yy_xzz_yyyz_0, ta2_yy_xzz_yyzz_0, ta2_yy_xzz_yzzz_0, ta2_yy_xzz_zzzz_0, ta2_yy_zz_xxx_0, ta2_yy_zz_xxx_1, ta2_yy_zz_xxxx_0, ta2_yy_zz_xxxx_1, ta2_yy_zz_xxxy_0, ta2_yy_zz_xxxy_1, ta2_yy_zz_xxxz_0, ta2_yy_zz_xxxz_1, ta2_yy_zz_xxy_0, ta2_yy_zz_xxy_1, ta2_yy_zz_xxyy_0, ta2_yy_zz_xxyy_1, ta2_yy_zz_xxyz_0, ta2_yy_zz_xxyz_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xxzz_0, ta2_yy_zz_xxzz_1, ta2_yy_zz_xyy_0, ta2_yy_zz_xyy_1, ta2_yy_zz_xyyy_0, ta2_yy_zz_xyyy_1, ta2_yy_zz_xyyz_0, ta2_yy_zz_xyyz_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xyzz_0, ta2_yy_zz_xyzz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_xzzz_0, ta2_yy_zz_xzzz_1, ta2_yy_zz_yyy_0, ta2_yy_zz_yyy_1, ta2_yy_zz_yyyy_0, ta2_yy_zz_yyyy_1, ta2_yy_zz_yyyz_0, ta2_yy_zz_yyyz_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yyzz_0, ta2_yy_zz_yyzz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_yzzz_0, ta2_yy_zz_yzzz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, ta2_yy_zz_zzzz_0, ta2_yy_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzz_xxxx_0[i] = 4.0 * ta2_yy_zz_xxx_0[i] * fe_0 - 4.0 * ta2_yy_zz_xxx_1[i] * fe_0 + ta2_yy_zz_xxxx_0[i] * pa_x[i] - ta2_yy_zz_xxxx_1[i] * pc_x[i];

        ta2_yy_xzz_xxxy_0[i] = 3.0 * ta2_yy_zz_xxy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxy_1[i] * fe_0 + ta2_yy_zz_xxxy_0[i] * pa_x[i] - ta2_yy_zz_xxxy_1[i] * pc_x[i];

        ta2_yy_xzz_xxxz_0[i] = 3.0 * ta2_yy_zz_xxz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxz_1[i] * fe_0 + ta2_yy_zz_xxxz_0[i] * pa_x[i] - ta2_yy_zz_xxxz_1[i] * pc_x[i];

        ta2_yy_xzz_xxyy_0[i] = 2.0 * ta2_yy_zz_xyy_0[i] * fe_0 - 2.0 * ta2_yy_zz_xyy_1[i] * fe_0 + ta2_yy_zz_xxyy_0[i] * pa_x[i] - ta2_yy_zz_xxyy_1[i] * pc_x[i];

        ta2_yy_xzz_xxyz_0[i] = 2.0 * ta2_yy_zz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xyz_1[i] * fe_0 + ta2_yy_zz_xxyz_0[i] * pa_x[i] - ta2_yy_zz_xxyz_1[i] * pc_x[i];

        ta2_yy_xzz_xxzz_0[i] = 2.0 * ta2_yy_zz_xzz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xzz_1[i] * fe_0 + ta2_yy_zz_xxzz_0[i] * pa_x[i] - ta2_yy_zz_xxzz_1[i] * pc_x[i];

        ta2_yy_xzz_xyyy_0[i] = ta2_yy_zz_yyy_0[i] * fe_0 - ta2_yy_zz_yyy_1[i] * fe_0 + ta2_yy_zz_xyyy_0[i] * pa_x[i] - ta2_yy_zz_xyyy_1[i] * pc_x[i];

        ta2_yy_xzz_xyyz_0[i] = ta2_yy_zz_yyz_0[i] * fe_0 - ta2_yy_zz_yyz_1[i] * fe_0 + ta2_yy_zz_xyyz_0[i] * pa_x[i] - ta2_yy_zz_xyyz_1[i] * pc_x[i];

        ta2_yy_xzz_xyzz_0[i] = ta2_yy_zz_yzz_0[i] * fe_0 - ta2_yy_zz_yzz_1[i] * fe_0 + ta2_yy_zz_xyzz_0[i] * pa_x[i] - ta2_yy_zz_xyzz_1[i] * pc_x[i];

        ta2_yy_xzz_xzzz_0[i] = ta2_yy_zz_zzz_0[i] * fe_0 - ta2_yy_zz_zzz_1[i] * fe_0 + ta2_yy_zz_xzzz_0[i] * pa_x[i] - ta2_yy_zz_xzzz_1[i] * pc_x[i];

        ta2_yy_xzz_yyyy_0[i] = ta2_yy_zz_yyyy_0[i] * pa_x[i] - ta2_yy_zz_yyyy_1[i] * pc_x[i];

        ta2_yy_xzz_yyyz_0[i] = ta2_yy_zz_yyyz_0[i] * pa_x[i] - ta2_yy_zz_yyyz_1[i] * pc_x[i];

        ta2_yy_xzz_yyzz_0[i] = ta2_yy_zz_yyzz_0[i] * pa_x[i] - ta2_yy_zz_yyzz_1[i] * pc_x[i];

        ta2_yy_xzz_yzzz_0[i] = ta2_yy_zz_yzzz_0[i] * pa_x[i] - ta2_yy_zz_yzzz_1[i] * pc_x[i];

        ta2_yy_xzz_zzzz_0[i] = ta2_yy_zz_zzzz_0[i] * pa_x[i] - ta2_yy_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 540-555 components of targeted buffer : FG

    auto ta2_yy_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 540);

    auto ta2_yy_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 541);

    auto ta2_yy_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 542);

    auto ta2_yy_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 543);

    auto ta2_yy_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 544);

    auto ta2_yy_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 545);

    auto ta2_yy_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 546);

    auto ta2_yy_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 547);

    auto ta2_yy_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 548);

    auto ta2_yy_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 549);

    auto ta2_yy_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 550);

    auto ta2_yy_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 551);

    auto ta2_yy_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 552);

    auto ta2_yy_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 553);

    auto ta2_yy_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 554);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yy_xxxx_1, ta1_y_yy_xxxy_1, ta1_y_yy_xxxz_1, ta1_y_yy_xxyy_1, ta1_y_yy_xxyz_1, ta1_y_yy_xxzz_1, ta1_y_yy_xyyy_1, ta1_y_yy_xyyz_1, ta1_y_yy_xyzz_1, ta1_y_yy_xzzz_1, ta1_y_yy_yyyy_1, ta1_y_yy_yyyz_1, ta1_y_yy_yyzz_1, ta1_y_yy_yzzz_1, ta1_y_yy_zzzz_1, ta2_yy_y_xxxx_0, ta2_yy_y_xxxx_1, ta2_yy_y_xxxy_0, ta2_yy_y_xxxy_1, ta2_yy_y_xxxz_0, ta2_yy_y_xxxz_1, ta2_yy_y_xxyy_0, ta2_yy_y_xxyy_1, ta2_yy_y_xxyz_0, ta2_yy_y_xxyz_1, ta2_yy_y_xxzz_0, ta2_yy_y_xxzz_1, ta2_yy_y_xyyy_0, ta2_yy_y_xyyy_1, ta2_yy_y_xyyz_0, ta2_yy_y_xyyz_1, ta2_yy_y_xyzz_0, ta2_yy_y_xyzz_1, ta2_yy_y_xzzz_0, ta2_yy_y_xzzz_1, ta2_yy_y_yyyy_0, ta2_yy_y_yyyy_1, ta2_yy_y_yyyz_0, ta2_yy_y_yyyz_1, ta2_yy_y_yyzz_0, ta2_yy_y_yyzz_1, ta2_yy_y_yzzz_0, ta2_yy_y_yzzz_1, ta2_yy_y_zzzz_0, ta2_yy_y_zzzz_1, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxxx_0, ta2_yy_yy_xxxx_1, ta2_yy_yy_xxxy_0, ta2_yy_yy_xxxy_1, ta2_yy_yy_xxxz_0, ta2_yy_yy_xxxz_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxyy_0, ta2_yy_yy_xxyy_1, ta2_yy_yy_xxyz_0, ta2_yy_yy_xxyz_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xxzz_0, ta2_yy_yy_xxzz_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyyy_0, ta2_yy_yy_xyyy_1, ta2_yy_yy_xyyz_0, ta2_yy_yy_xyyz_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xyzz_0, ta2_yy_yy_xyzz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_xzzz_0, ta2_yy_yy_xzzz_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyyy_0, ta2_yy_yy_yyyy_1, ta2_yy_yy_yyyz_0, ta2_yy_yy_yyyz_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yyzz_0, ta2_yy_yy_yyzz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_yzzz_0, ta2_yy_yy_yzzz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, ta2_yy_yy_zzzz_0, ta2_yy_yy_zzzz_1, ta2_yy_yyy_xxxx_0, ta2_yy_yyy_xxxy_0, ta2_yy_yyy_xxxz_0, ta2_yy_yyy_xxyy_0, ta2_yy_yyy_xxyz_0, ta2_yy_yyy_xxzz_0, ta2_yy_yyy_xyyy_0, ta2_yy_yyy_xyyz_0, ta2_yy_yyy_xyzz_0, ta2_yy_yyy_xzzz_0, ta2_yy_yyy_yyyy_0, ta2_yy_yyy_yyyz_0, ta2_yy_yyy_yyzz_0, ta2_yy_yyy_yzzz_0, ta2_yy_yyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyy_xxxx_0[i] = 2.0 * ta2_yy_y_xxxx_0[i] * fe_0 - 2.0 * ta2_yy_y_xxxx_1[i] * fe_0 + 2.0 * ta1_y_yy_xxxx_1[i] + ta2_yy_yy_xxxx_0[i] * pa_y[i] - ta2_yy_yy_xxxx_1[i] * pc_y[i];

        ta2_yy_yyy_xxxy_0[i] = 2.0 * ta2_yy_y_xxxy_0[i] * fe_0 - 2.0 * ta2_yy_y_xxxy_1[i] * fe_0 + ta2_yy_yy_xxx_0[i] * fe_0 - ta2_yy_yy_xxx_1[i] * fe_0 + 2.0 * ta1_y_yy_xxxy_1[i] + ta2_yy_yy_xxxy_0[i] * pa_y[i] - ta2_yy_yy_xxxy_1[i] * pc_y[i];

        ta2_yy_yyy_xxxz_0[i] = 2.0 * ta2_yy_y_xxxz_0[i] * fe_0 - 2.0 * ta2_yy_y_xxxz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxxz_1[i] + ta2_yy_yy_xxxz_0[i] * pa_y[i] - ta2_yy_yy_xxxz_1[i] * pc_y[i];

        ta2_yy_yyy_xxyy_0[i] = 2.0 * ta2_yy_y_xxyy_0[i] * fe_0 - 2.0 * ta2_yy_y_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_yy_xxy_0[i] * fe_0 - 2.0 * ta2_yy_yy_xxy_1[i] * fe_0 + 2.0 * ta1_y_yy_xxyy_1[i] + ta2_yy_yy_xxyy_0[i] * pa_y[i] - ta2_yy_yy_xxyy_1[i] * pc_y[i];

        ta2_yy_yyy_xxyz_0[i] = 2.0 * ta2_yy_y_xxyz_0[i] * fe_0 - 2.0 * ta2_yy_y_xxyz_1[i] * fe_0 + ta2_yy_yy_xxz_0[i] * fe_0 - ta2_yy_yy_xxz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxyz_1[i] + ta2_yy_yy_xxyz_0[i] * pa_y[i] - ta2_yy_yy_xxyz_1[i] * pc_y[i];

        ta2_yy_yyy_xxzz_0[i] = 2.0 * ta2_yy_y_xxzz_0[i] * fe_0 - 2.0 * ta2_yy_y_xxzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxzz_1[i] + ta2_yy_yy_xxzz_0[i] * pa_y[i] - ta2_yy_yy_xxzz_1[i] * pc_y[i];

        ta2_yy_yyy_xyyy_0[i] = 2.0 * ta2_yy_y_xyyy_0[i] * fe_0 - 2.0 * ta2_yy_y_xyyy_1[i] * fe_0 + 3.0 * ta2_yy_yy_xyy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xyy_1[i] * fe_0 + 2.0 * ta1_y_yy_xyyy_1[i] + ta2_yy_yy_xyyy_0[i] * pa_y[i] - ta2_yy_yy_xyyy_1[i] * pc_y[i];

        ta2_yy_yyy_xyyz_0[i] = 2.0 * ta2_yy_y_xyyz_0[i] * fe_0 - 2.0 * ta2_yy_y_xyyz_1[i] * fe_0 + 2.0 * ta2_yy_yy_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xyz_1[i] * fe_0 + 2.0 * ta1_y_yy_xyyz_1[i] + ta2_yy_yy_xyyz_0[i] * pa_y[i] - ta2_yy_yy_xyyz_1[i] * pc_y[i];

        ta2_yy_yyy_xyzz_0[i] = 2.0 * ta2_yy_y_xyzz_0[i] * fe_0 - 2.0 * ta2_yy_y_xyzz_1[i] * fe_0 + ta2_yy_yy_xzz_0[i] * fe_0 - ta2_yy_yy_xzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xyzz_1[i] + ta2_yy_yy_xyzz_0[i] * pa_y[i] - ta2_yy_yy_xyzz_1[i] * pc_y[i];

        ta2_yy_yyy_xzzz_0[i] = 2.0 * ta2_yy_y_xzzz_0[i] * fe_0 - 2.0 * ta2_yy_y_xzzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xzzz_1[i] + ta2_yy_yy_xzzz_0[i] * pa_y[i] - ta2_yy_yy_xzzz_1[i] * pc_y[i];

        ta2_yy_yyy_yyyy_0[i] = 2.0 * ta2_yy_y_yyyy_0[i] * fe_0 - 2.0 * ta2_yy_y_yyyy_1[i] * fe_0 + 4.0 * ta2_yy_yy_yyy_0[i] * fe_0 - 4.0 * ta2_yy_yy_yyy_1[i] * fe_0 + 2.0 * ta1_y_yy_yyyy_1[i] + ta2_yy_yy_yyyy_0[i] * pa_y[i] - ta2_yy_yy_yyyy_1[i] * pc_y[i];

        ta2_yy_yyy_yyyz_0[i] = 2.0 * ta2_yy_y_yyyz_0[i] * fe_0 - 2.0 * ta2_yy_y_yyyz_1[i] * fe_0 + 3.0 * ta2_yy_yy_yyz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yyz_1[i] * fe_0 + 2.0 * ta1_y_yy_yyyz_1[i] + ta2_yy_yy_yyyz_0[i] * pa_y[i] - ta2_yy_yy_yyyz_1[i] * pc_y[i];

        ta2_yy_yyy_yyzz_0[i] = 2.0 * ta2_yy_y_yyzz_0[i] * fe_0 - 2.0 * ta2_yy_y_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_yy_yzz_0[i] * fe_0 - 2.0 * ta2_yy_yy_yzz_1[i] * fe_0 + 2.0 * ta1_y_yy_yyzz_1[i] + ta2_yy_yy_yyzz_0[i] * pa_y[i] - ta2_yy_yy_yyzz_1[i] * pc_y[i];

        ta2_yy_yyy_yzzz_0[i] = 2.0 * ta2_yy_y_yzzz_0[i] * fe_0 - 2.0 * ta2_yy_y_yzzz_1[i] * fe_0 + ta2_yy_yy_zzz_0[i] * fe_0 - ta2_yy_yy_zzz_1[i] * fe_0 + 2.0 * ta1_y_yy_yzzz_1[i] + ta2_yy_yy_yzzz_0[i] * pa_y[i] - ta2_yy_yy_yzzz_1[i] * pc_y[i];

        ta2_yy_yyy_zzzz_0[i] = 2.0 * ta2_yy_y_zzzz_0[i] * fe_0 - 2.0 * ta2_yy_y_zzzz_1[i] * fe_0 + 2.0 * ta1_y_yy_zzzz_1[i] + ta2_yy_yy_zzzz_0[i] * pa_y[i] - ta2_yy_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 555-570 components of targeted buffer : FG

    auto ta2_yy_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 555);

    auto ta2_yy_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 556);

    auto ta2_yy_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 557);

    auto ta2_yy_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 558);

    auto ta2_yy_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 559);

    auto ta2_yy_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 560);

    auto ta2_yy_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 561);

    auto ta2_yy_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 562);

    auto ta2_yy_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 563);

    auto ta2_yy_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 564);

    auto ta2_yy_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 565);

    auto ta2_yy_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 566);

    auto ta2_yy_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 567);

    auto ta2_yy_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 568);

    auto ta2_yy_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 569);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxxx_0, ta2_yy_yy_xxxx_1, ta2_yy_yy_xxxy_0, ta2_yy_yy_xxxy_1, ta2_yy_yy_xxxz_0, ta2_yy_yy_xxxz_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxyy_0, ta2_yy_yy_xxyy_1, ta2_yy_yy_xxyz_0, ta2_yy_yy_xxyz_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xxzz_0, ta2_yy_yy_xxzz_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyyy_0, ta2_yy_yy_xyyy_1, ta2_yy_yy_xyyz_0, ta2_yy_yy_xyyz_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xyzz_0, ta2_yy_yy_xyzz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_xzzz_0, ta2_yy_yy_xzzz_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyyy_0, ta2_yy_yy_yyyy_1, ta2_yy_yy_yyyz_0, ta2_yy_yy_yyyz_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yyzz_0, ta2_yy_yy_yyzz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_yzzz_0, ta2_yy_yy_yzzz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, ta2_yy_yy_zzzz_0, ta2_yy_yy_zzzz_1, ta2_yy_yyz_xxxx_0, ta2_yy_yyz_xxxy_0, ta2_yy_yyz_xxxz_0, ta2_yy_yyz_xxyy_0, ta2_yy_yyz_xxyz_0, ta2_yy_yyz_xxzz_0, ta2_yy_yyz_xyyy_0, ta2_yy_yyz_xyyz_0, ta2_yy_yyz_xyzz_0, ta2_yy_yyz_xzzz_0, ta2_yy_yyz_yyyy_0, ta2_yy_yyz_yyyz_0, ta2_yy_yyz_yyzz_0, ta2_yy_yyz_yzzz_0, ta2_yy_yyz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyz_xxxx_0[i] = ta2_yy_yy_xxxx_0[i] * pa_z[i] - ta2_yy_yy_xxxx_1[i] * pc_z[i];

        ta2_yy_yyz_xxxy_0[i] = ta2_yy_yy_xxxy_0[i] * pa_z[i] - ta2_yy_yy_xxxy_1[i] * pc_z[i];

        ta2_yy_yyz_xxxz_0[i] = ta2_yy_yy_xxx_0[i] * fe_0 - ta2_yy_yy_xxx_1[i] * fe_0 + ta2_yy_yy_xxxz_0[i] * pa_z[i] - ta2_yy_yy_xxxz_1[i] * pc_z[i];

        ta2_yy_yyz_xxyy_0[i] = ta2_yy_yy_xxyy_0[i] * pa_z[i] - ta2_yy_yy_xxyy_1[i] * pc_z[i];

        ta2_yy_yyz_xxyz_0[i] = ta2_yy_yy_xxy_0[i] * fe_0 - ta2_yy_yy_xxy_1[i] * fe_0 + ta2_yy_yy_xxyz_0[i] * pa_z[i] - ta2_yy_yy_xxyz_1[i] * pc_z[i];

        ta2_yy_yyz_xxzz_0[i] = 2.0 * ta2_yy_yy_xxz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xxz_1[i] * fe_0 + ta2_yy_yy_xxzz_0[i] * pa_z[i] - ta2_yy_yy_xxzz_1[i] * pc_z[i];

        ta2_yy_yyz_xyyy_0[i] = ta2_yy_yy_xyyy_0[i] * pa_z[i] - ta2_yy_yy_xyyy_1[i] * pc_z[i];

        ta2_yy_yyz_xyyz_0[i] = ta2_yy_yy_xyy_0[i] * fe_0 - ta2_yy_yy_xyy_1[i] * fe_0 + ta2_yy_yy_xyyz_0[i] * pa_z[i] - ta2_yy_yy_xyyz_1[i] * pc_z[i];

        ta2_yy_yyz_xyzz_0[i] = 2.0 * ta2_yy_yy_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xyz_1[i] * fe_0 + ta2_yy_yy_xyzz_0[i] * pa_z[i] - ta2_yy_yy_xyzz_1[i] * pc_z[i];

        ta2_yy_yyz_xzzz_0[i] = 3.0 * ta2_yy_yy_xzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xzz_1[i] * fe_0 + ta2_yy_yy_xzzz_0[i] * pa_z[i] - ta2_yy_yy_xzzz_1[i] * pc_z[i];

        ta2_yy_yyz_yyyy_0[i] = ta2_yy_yy_yyyy_0[i] * pa_z[i] - ta2_yy_yy_yyyy_1[i] * pc_z[i];

        ta2_yy_yyz_yyyz_0[i] = ta2_yy_yy_yyy_0[i] * fe_0 - ta2_yy_yy_yyy_1[i] * fe_0 + ta2_yy_yy_yyyz_0[i] * pa_z[i] - ta2_yy_yy_yyyz_1[i] * pc_z[i];

        ta2_yy_yyz_yyzz_0[i] = 2.0 * ta2_yy_yy_yyz_0[i] * fe_0 - 2.0 * ta2_yy_yy_yyz_1[i] * fe_0 + ta2_yy_yy_yyzz_0[i] * pa_z[i] - ta2_yy_yy_yyzz_1[i] * pc_z[i];

        ta2_yy_yyz_yzzz_0[i] = 3.0 * ta2_yy_yy_yzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yzz_1[i] * fe_0 + ta2_yy_yy_yzzz_0[i] * pa_z[i] - ta2_yy_yy_yzzz_1[i] * pc_z[i];

        ta2_yy_yyz_zzzz_0[i] = 4.0 * ta2_yy_yy_zzz_0[i] * fe_0 - 4.0 * ta2_yy_yy_zzz_1[i] * fe_0 + ta2_yy_yy_zzzz_0[i] * pa_z[i] - ta2_yy_yy_zzzz_1[i] * pc_z[i];
    }

    // Set up 570-585 components of targeted buffer : FG

    auto ta2_yy_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 570);

    auto ta2_yy_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 571);

    auto ta2_yy_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 572);

    auto ta2_yy_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 573);

    auto ta2_yy_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 574);

    auto ta2_yy_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 575);

    auto ta2_yy_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 576);

    auto ta2_yy_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 577);

    auto ta2_yy_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 578);

    auto ta2_yy_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 579);

    auto ta2_yy_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 580);

    auto ta2_yy_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 581);

    auto ta2_yy_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 582);

    auto ta2_yy_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 583);

    auto ta2_yy_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 584);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_zz_xxxx_1, ta1_y_zz_xxxz_1, ta1_y_zz_xxyz_1, ta1_y_zz_xxzz_1, ta1_y_zz_xyyz_1, ta1_y_zz_xyzz_1, ta1_y_zz_xzzz_1, ta1_y_zz_yyyz_1, ta1_y_zz_yyzz_1, ta1_y_zz_yzzz_1, ta1_y_zz_zzzz_1, ta2_yy_y_xxxy_0, ta2_yy_y_xxxy_1, ta2_yy_y_xxyy_0, ta2_yy_y_xxyy_1, ta2_yy_y_xyyy_0, ta2_yy_y_xyyy_1, ta2_yy_y_yyyy_0, ta2_yy_y_yyyy_1, ta2_yy_yz_xxxy_0, ta2_yy_yz_xxxy_1, ta2_yy_yz_xxyy_0, ta2_yy_yz_xxyy_1, ta2_yy_yz_xyyy_0, ta2_yy_yz_xyyy_1, ta2_yy_yz_yyyy_0, ta2_yy_yz_yyyy_1, ta2_yy_yzz_xxxx_0, ta2_yy_yzz_xxxy_0, ta2_yy_yzz_xxxz_0, ta2_yy_yzz_xxyy_0, ta2_yy_yzz_xxyz_0, ta2_yy_yzz_xxzz_0, ta2_yy_yzz_xyyy_0, ta2_yy_yzz_xyyz_0, ta2_yy_yzz_xyzz_0, ta2_yy_yzz_xzzz_0, ta2_yy_yzz_yyyy_0, ta2_yy_yzz_yyyz_0, ta2_yy_yzz_yyzz_0, ta2_yy_yzz_yzzz_0, ta2_yy_yzz_zzzz_0, ta2_yy_zz_xxxx_0, ta2_yy_zz_xxxx_1, ta2_yy_zz_xxxz_0, ta2_yy_zz_xxxz_1, ta2_yy_zz_xxyz_0, ta2_yy_zz_xxyz_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xxzz_0, ta2_yy_zz_xxzz_1, ta2_yy_zz_xyyz_0, ta2_yy_zz_xyyz_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xyzz_0, ta2_yy_zz_xyzz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_xzzz_0, ta2_yy_zz_xzzz_1, ta2_yy_zz_yyyz_0, ta2_yy_zz_yyyz_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yyzz_0, ta2_yy_zz_yyzz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_yzzz_0, ta2_yy_zz_yzzz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, ta2_yy_zz_zzzz_0, ta2_yy_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzz_xxxx_0[i] = 2.0 * ta1_y_zz_xxxx_1[i] + ta2_yy_zz_xxxx_0[i] * pa_y[i] - ta2_yy_zz_xxxx_1[i] * pc_y[i];

        ta2_yy_yzz_xxxy_0[i] = ta2_yy_y_xxxy_0[i] * fe_0 - ta2_yy_y_xxxy_1[i] * fe_0 + ta2_yy_yz_xxxy_0[i] * pa_z[i] - ta2_yy_yz_xxxy_1[i] * pc_z[i];

        ta2_yy_yzz_xxxz_0[i] = 2.0 * ta1_y_zz_xxxz_1[i] + ta2_yy_zz_xxxz_0[i] * pa_y[i] - ta2_yy_zz_xxxz_1[i] * pc_y[i];

        ta2_yy_yzz_xxyy_0[i] = ta2_yy_y_xxyy_0[i] * fe_0 - ta2_yy_y_xxyy_1[i] * fe_0 + ta2_yy_yz_xxyy_0[i] * pa_z[i] - ta2_yy_yz_xxyy_1[i] * pc_z[i];

        ta2_yy_yzz_xxyz_0[i] = ta2_yy_zz_xxz_0[i] * fe_0 - ta2_yy_zz_xxz_1[i] * fe_0 + 2.0 * ta1_y_zz_xxyz_1[i] + ta2_yy_zz_xxyz_0[i] * pa_y[i] - ta2_yy_zz_xxyz_1[i] * pc_y[i];

        ta2_yy_yzz_xxzz_0[i] = 2.0 * ta1_y_zz_xxzz_1[i] + ta2_yy_zz_xxzz_0[i] * pa_y[i] - ta2_yy_zz_xxzz_1[i] * pc_y[i];

        ta2_yy_yzz_xyyy_0[i] = ta2_yy_y_xyyy_0[i] * fe_0 - ta2_yy_y_xyyy_1[i] * fe_0 + ta2_yy_yz_xyyy_0[i] * pa_z[i] - ta2_yy_yz_xyyy_1[i] * pc_z[i];

        ta2_yy_yzz_xyyz_0[i] = 2.0 * ta2_yy_zz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xyz_1[i] * fe_0 + 2.0 * ta1_y_zz_xyyz_1[i] + ta2_yy_zz_xyyz_0[i] * pa_y[i] - ta2_yy_zz_xyyz_1[i] * pc_y[i];

        ta2_yy_yzz_xyzz_0[i] = ta2_yy_zz_xzz_0[i] * fe_0 - ta2_yy_zz_xzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xyzz_1[i] + ta2_yy_zz_xyzz_0[i] * pa_y[i] - ta2_yy_zz_xyzz_1[i] * pc_y[i];

        ta2_yy_yzz_xzzz_0[i] = 2.0 * ta1_y_zz_xzzz_1[i] + ta2_yy_zz_xzzz_0[i] * pa_y[i] - ta2_yy_zz_xzzz_1[i] * pc_y[i];

        ta2_yy_yzz_yyyy_0[i] = ta2_yy_y_yyyy_0[i] * fe_0 - ta2_yy_y_yyyy_1[i] * fe_0 + ta2_yy_yz_yyyy_0[i] * pa_z[i] - ta2_yy_yz_yyyy_1[i] * pc_z[i];

        ta2_yy_yzz_yyyz_0[i] = 3.0 * ta2_yy_zz_yyz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yyz_1[i] * fe_0 + 2.0 * ta1_y_zz_yyyz_1[i] + ta2_yy_zz_yyyz_0[i] * pa_y[i] - ta2_yy_zz_yyyz_1[i] * pc_y[i];

        ta2_yy_yzz_yyzz_0[i] = 2.0 * ta2_yy_zz_yzz_0[i] * fe_0 - 2.0 * ta2_yy_zz_yzz_1[i] * fe_0 + 2.0 * ta1_y_zz_yyzz_1[i] + ta2_yy_zz_yyzz_0[i] * pa_y[i] - ta2_yy_zz_yyzz_1[i] * pc_y[i];

        ta2_yy_yzz_yzzz_0[i] = ta2_yy_zz_zzz_0[i] * fe_0 - ta2_yy_zz_zzz_1[i] * fe_0 + 2.0 * ta1_y_zz_yzzz_1[i] + ta2_yy_zz_yzzz_0[i] * pa_y[i] - ta2_yy_zz_yzzz_1[i] * pc_y[i];

        ta2_yy_yzz_zzzz_0[i] = 2.0 * ta1_y_zz_zzzz_1[i] + ta2_yy_zz_zzzz_0[i] * pa_y[i] - ta2_yy_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 585-600 components of targeted buffer : FG

    auto ta2_yy_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 585);

    auto ta2_yy_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 586);

    auto ta2_yy_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 587);

    auto ta2_yy_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 588);

    auto ta2_yy_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 589);

    auto ta2_yy_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 590);

    auto ta2_yy_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 591);

    auto ta2_yy_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 592);

    auto ta2_yy_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 593);

    auto ta2_yy_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 594);

    auto ta2_yy_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 595);

    auto ta2_yy_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 596);

    auto ta2_yy_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 597);

    auto ta2_yy_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 598);

    auto ta2_yy_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 599);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_z_xxxx_0, ta2_yy_z_xxxx_1, ta2_yy_z_xxxy_0, ta2_yy_z_xxxy_1, ta2_yy_z_xxxz_0, ta2_yy_z_xxxz_1, ta2_yy_z_xxyy_0, ta2_yy_z_xxyy_1, ta2_yy_z_xxyz_0, ta2_yy_z_xxyz_1, ta2_yy_z_xxzz_0, ta2_yy_z_xxzz_1, ta2_yy_z_xyyy_0, ta2_yy_z_xyyy_1, ta2_yy_z_xyyz_0, ta2_yy_z_xyyz_1, ta2_yy_z_xyzz_0, ta2_yy_z_xyzz_1, ta2_yy_z_xzzz_0, ta2_yy_z_xzzz_1, ta2_yy_z_yyyy_0, ta2_yy_z_yyyy_1, ta2_yy_z_yyyz_0, ta2_yy_z_yyyz_1, ta2_yy_z_yyzz_0, ta2_yy_z_yyzz_1, ta2_yy_z_yzzz_0, ta2_yy_z_yzzz_1, ta2_yy_z_zzzz_0, ta2_yy_z_zzzz_1, ta2_yy_zz_xxx_0, ta2_yy_zz_xxx_1, ta2_yy_zz_xxxx_0, ta2_yy_zz_xxxx_1, ta2_yy_zz_xxxy_0, ta2_yy_zz_xxxy_1, ta2_yy_zz_xxxz_0, ta2_yy_zz_xxxz_1, ta2_yy_zz_xxy_0, ta2_yy_zz_xxy_1, ta2_yy_zz_xxyy_0, ta2_yy_zz_xxyy_1, ta2_yy_zz_xxyz_0, ta2_yy_zz_xxyz_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xxzz_0, ta2_yy_zz_xxzz_1, ta2_yy_zz_xyy_0, ta2_yy_zz_xyy_1, ta2_yy_zz_xyyy_0, ta2_yy_zz_xyyy_1, ta2_yy_zz_xyyz_0, ta2_yy_zz_xyyz_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xyzz_0, ta2_yy_zz_xyzz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_xzzz_0, ta2_yy_zz_xzzz_1, ta2_yy_zz_yyy_0, ta2_yy_zz_yyy_1, ta2_yy_zz_yyyy_0, ta2_yy_zz_yyyy_1, ta2_yy_zz_yyyz_0, ta2_yy_zz_yyyz_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yyzz_0, ta2_yy_zz_yyzz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_yzzz_0, ta2_yy_zz_yzzz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, ta2_yy_zz_zzzz_0, ta2_yy_zz_zzzz_1, ta2_yy_zzz_xxxx_0, ta2_yy_zzz_xxxy_0, ta2_yy_zzz_xxxz_0, ta2_yy_zzz_xxyy_0, ta2_yy_zzz_xxyz_0, ta2_yy_zzz_xxzz_0, ta2_yy_zzz_xyyy_0, ta2_yy_zzz_xyyz_0, ta2_yy_zzz_xyzz_0, ta2_yy_zzz_xzzz_0, ta2_yy_zzz_yyyy_0, ta2_yy_zzz_yyyz_0, ta2_yy_zzz_yyzz_0, ta2_yy_zzz_yzzz_0, ta2_yy_zzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzz_xxxx_0[i] = 2.0 * ta2_yy_z_xxxx_0[i] * fe_0 - 2.0 * ta2_yy_z_xxxx_1[i] * fe_0 + ta2_yy_zz_xxxx_0[i] * pa_z[i] - ta2_yy_zz_xxxx_1[i] * pc_z[i];

        ta2_yy_zzz_xxxy_0[i] = 2.0 * ta2_yy_z_xxxy_0[i] * fe_0 - 2.0 * ta2_yy_z_xxxy_1[i] * fe_0 + ta2_yy_zz_xxxy_0[i] * pa_z[i] - ta2_yy_zz_xxxy_1[i] * pc_z[i];

        ta2_yy_zzz_xxxz_0[i] = 2.0 * ta2_yy_z_xxxz_0[i] * fe_0 - 2.0 * ta2_yy_z_xxxz_1[i] * fe_0 + ta2_yy_zz_xxx_0[i] * fe_0 - ta2_yy_zz_xxx_1[i] * fe_0 + ta2_yy_zz_xxxz_0[i] * pa_z[i] - ta2_yy_zz_xxxz_1[i] * pc_z[i];

        ta2_yy_zzz_xxyy_0[i] = 2.0 * ta2_yy_z_xxyy_0[i] * fe_0 - 2.0 * ta2_yy_z_xxyy_1[i] * fe_0 + ta2_yy_zz_xxyy_0[i] * pa_z[i] - ta2_yy_zz_xxyy_1[i] * pc_z[i];

        ta2_yy_zzz_xxyz_0[i] = 2.0 * ta2_yy_z_xxyz_0[i] * fe_0 - 2.0 * ta2_yy_z_xxyz_1[i] * fe_0 + ta2_yy_zz_xxy_0[i] * fe_0 - ta2_yy_zz_xxy_1[i] * fe_0 + ta2_yy_zz_xxyz_0[i] * pa_z[i] - ta2_yy_zz_xxyz_1[i] * pc_z[i];

        ta2_yy_zzz_xxzz_0[i] = 2.0 * ta2_yy_z_xxzz_0[i] * fe_0 - 2.0 * ta2_yy_z_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_zz_xxz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xxz_1[i] * fe_0 + ta2_yy_zz_xxzz_0[i] * pa_z[i] - ta2_yy_zz_xxzz_1[i] * pc_z[i];

        ta2_yy_zzz_xyyy_0[i] = 2.0 * ta2_yy_z_xyyy_0[i] * fe_0 - 2.0 * ta2_yy_z_xyyy_1[i] * fe_0 + ta2_yy_zz_xyyy_0[i] * pa_z[i] - ta2_yy_zz_xyyy_1[i] * pc_z[i];

        ta2_yy_zzz_xyyz_0[i] = 2.0 * ta2_yy_z_xyyz_0[i] * fe_0 - 2.0 * ta2_yy_z_xyyz_1[i] * fe_0 + ta2_yy_zz_xyy_0[i] * fe_0 - ta2_yy_zz_xyy_1[i] * fe_0 + ta2_yy_zz_xyyz_0[i] * pa_z[i] - ta2_yy_zz_xyyz_1[i] * pc_z[i];

        ta2_yy_zzz_xyzz_0[i] = 2.0 * ta2_yy_z_xyzz_0[i] * fe_0 - 2.0 * ta2_yy_z_xyzz_1[i] * fe_0 + 2.0 * ta2_yy_zz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xyz_1[i] * fe_0 + ta2_yy_zz_xyzz_0[i] * pa_z[i] - ta2_yy_zz_xyzz_1[i] * pc_z[i];

        ta2_yy_zzz_xzzz_0[i] = 2.0 * ta2_yy_z_xzzz_0[i] * fe_0 - 2.0 * ta2_yy_z_xzzz_1[i] * fe_0 + 3.0 * ta2_yy_zz_xzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xzz_1[i] * fe_0 + ta2_yy_zz_xzzz_0[i] * pa_z[i] - ta2_yy_zz_xzzz_1[i] * pc_z[i];

        ta2_yy_zzz_yyyy_0[i] = 2.0 * ta2_yy_z_yyyy_0[i] * fe_0 - 2.0 * ta2_yy_z_yyyy_1[i] * fe_0 + ta2_yy_zz_yyyy_0[i] * pa_z[i] - ta2_yy_zz_yyyy_1[i] * pc_z[i];

        ta2_yy_zzz_yyyz_0[i] = 2.0 * ta2_yy_z_yyyz_0[i] * fe_0 - 2.0 * ta2_yy_z_yyyz_1[i] * fe_0 + ta2_yy_zz_yyy_0[i] * fe_0 - ta2_yy_zz_yyy_1[i] * fe_0 + ta2_yy_zz_yyyz_0[i] * pa_z[i] - ta2_yy_zz_yyyz_1[i] * pc_z[i];

        ta2_yy_zzz_yyzz_0[i] = 2.0 * ta2_yy_z_yyzz_0[i] * fe_0 - 2.0 * ta2_yy_z_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_zz_yyz_0[i] * fe_0 - 2.0 * ta2_yy_zz_yyz_1[i] * fe_0 + ta2_yy_zz_yyzz_0[i] * pa_z[i] - ta2_yy_zz_yyzz_1[i] * pc_z[i];

        ta2_yy_zzz_yzzz_0[i] = 2.0 * ta2_yy_z_yzzz_0[i] * fe_0 - 2.0 * ta2_yy_z_yzzz_1[i] * fe_0 + 3.0 * ta2_yy_zz_yzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yzz_1[i] * fe_0 + ta2_yy_zz_yzzz_0[i] * pa_z[i] - ta2_yy_zz_yzzz_1[i] * pc_z[i];

        ta2_yy_zzz_zzzz_0[i] = 2.0 * ta2_yy_z_zzzz_0[i] * fe_0 - 2.0 * ta2_yy_z_zzzz_1[i] * fe_0 + 4.0 * ta2_yy_zz_zzz_0[i] * fe_0 - 4.0 * ta2_yy_zz_zzz_1[i] * fe_0 + ta2_yy_zz_zzzz_0[i] * pa_z[i] - ta2_yy_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 600-615 components of targeted buffer : FG

    auto ta2_yz_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 600);

    auto ta2_yz_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 601);

    auto ta2_yz_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 602);

    auto ta2_yz_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 603);

    auto ta2_yz_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 604);

    auto ta2_yz_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 605);

    auto ta2_yz_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 606);

    auto ta2_yz_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 607);

    auto ta2_yz_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 608);

    auto ta2_yz_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 609);

    auto ta2_yz_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 610);

    auto ta2_yz_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 611);

    auto ta2_yz_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 612);

    auto ta2_yz_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 613);

    auto ta2_yz_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 614);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_x_xxxx_0, ta2_yz_x_xxxx_1, ta2_yz_x_xxxy_0, ta2_yz_x_xxxy_1, ta2_yz_x_xxxz_0, ta2_yz_x_xxxz_1, ta2_yz_x_xxyy_0, ta2_yz_x_xxyy_1, ta2_yz_x_xxyz_0, ta2_yz_x_xxyz_1, ta2_yz_x_xxzz_0, ta2_yz_x_xxzz_1, ta2_yz_x_xyyy_0, ta2_yz_x_xyyy_1, ta2_yz_x_xyyz_0, ta2_yz_x_xyyz_1, ta2_yz_x_xyzz_0, ta2_yz_x_xyzz_1, ta2_yz_x_xzzz_0, ta2_yz_x_xzzz_1, ta2_yz_x_yyyy_0, ta2_yz_x_yyyy_1, ta2_yz_x_yyyz_0, ta2_yz_x_yyyz_1, ta2_yz_x_yyzz_0, ta2_yz_x_yyzz_1, ta2_yz_x_yzzz_0, ta2_yz_x_yzzz_1, ta2_yz_x_zzzz_0, ta2_yz_x_zzzz_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxxx_0, ta2_yz_xx_xxxx_1, ta2_yz_xx_xxxy_0, ta2_yz_xx_xxxy_1, ta2_yz_xx_xxxz_0, ta2_yz_xx_xxxz_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxyy_0, ta2_yz_xx_xxyy_1, ta2_yz_xx_xxyz_0, ta2_yz_xx_xxyz_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xxzz_0, ta2_yz_xx_xxzz_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyyy_0, ta2_yz_xx_xyyy_1, ta2_yz_xx_xyyz_0, ta2_yz_xx_xyyz_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xyzz_0, ta2_yz_xx_xyzz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_xzzz_0, ta2_yz_xx_xzzz_1, ta2_yz_xx_yyy_0, ta2_yz_xx_yyy_1, ta2_yz_xx_yyyy_0, ta2_yz_xx_yyyy_1, ta2_yz_xx_yyyz_0, ta2_yz_xx_yyyz_1, ta2_yz_xx_yyz_0, ta2_yz_xx_yyz_1, ta2_yz_xx_yyzz_0, ta2_yz_xx_yyzz_1, ta2_yz_xx_yzz_0, ta2_yz_xx_yzz_1, ta2_yz_xx_yzzz_0, ta2_yz_xx_yzzz_1, ta2_yz_xx_zzz_0, ta2_yz_xx_zzz_1, ta2_yz_xx_zzzz_0, ta2_yz_xx_zzzz_1, ta2_yz_xxx_xxxx_0, ta2_yz_xxx_xxxy_0, ta2_yz_xxx_xxxz_0, ta2_yz_xxx_xxyy_0, ta2_yz_xxx_xxyz_0, ta2_yz_xxx_xxzz_0, ta2_yz_xxx_xyyy_0, ta2_yz_xxx_xyyz_0, ta2_yz_xxx_xyzz_0, ta2_yz_xxx_xzzz_0, ta2_yz_xxx_yyyy_0, ta2_yz_xxx_yyyz_0, ta2_yz_xxx_yyzz_0, ta2_yz_xxx_yzzz_0, ta2_yz_xxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxx_xxxx_0[i] = 2.0 * ta2_yz_x_xxxx_0[i] * fe_0 - 2.0 * ta2_yz_x_xxxx_1[i] * fe_0 + 4.0 * ta2_yz_xx_xxx_0[i] * fe_0 - 4.0 * ta2_yz_xx_xxx_1[i] * fe_0 + ta2_yz_xx_xxxx_0[i] * pa_x[i] - ta2_yz_xx_xxxx_1[i] * pc_x[i];

        ta2_yz_xxx_xxxy_0[i] = 2.0 * ta2_yz_x_xxxy_0[i] * fe_0 - 2.0 * ta2_yz_x_xxxy_1[i] * fe_0 + 3.0 * ta2_yz_xx_xxy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxy_1[i] * fe_0 + ta2_yz_xx_xxxy_0[i] * pa_x[i] - ta2_yz_xx_xxxy_1[i] * pc_x[i];

        ta2_yz_xxx_xxxz_0[i] = 2.0 * ta2_yz_x_xxxz_0[i] * fe_0 - 2.0 * ta2_yz_x_xxxz_1[i] * fe_0 + 3.0 * ta2_yz_xx_xxz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxz_1[i] * fe_0 + ta2_yz_xx_xxxz_0[i] * pa_x[i] - ta2_yz_xx_xxxz_1[i] * pc_x[i];

        ta2_yz_xxx_xxyy_0[i] = 2.0 * ta2_yz_x_xxyy_0[i] * fe_0 - 2.0 * ta2_yz_x_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_xx_xyy_0[i] * fe_0 - 2.0 * ta2_yz_xx_xyy_1[i] * fe_0 + ta2_yz_xx_xxyy_0[i] * pa_x[i] - ta2_yz_xx_xxyy_1[i] * pc_x[i];

        ta2_yz_xxx_xxyz_0[i] = 2.0 * ta2_yz_x_xxyz_0[i] * fe_0 - 2.0 * ta2_yz_x_xxyz_1[i] * fe_0 + 2.0 * ta2_yz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xyz_1[i] * fe_0 + ta2_yz_xx_xxyz_0[i] * pa_x[i] - ta2_yz_xx_xxyz_1[i] * pc_x[i];

        ta2_yz_xxx_xxzz_0[i] = 2.0 * ta2_yz_x_xxzz_0[i] * fe_0 - 2.0 * ta2_yz_x_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_xx_xzz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xzz_1[i] * fe_0 + ta2_yz_xx_xxzz_0[i] * pa_x[i] - ta2_yz_xx_xxzz_1[i] * pc_x[i];

        ta2_yz_xxx_xyyy_0[i] = 2.0 * ta2_yz_x_xyyy_0[i] * fe_0 - 2.0 * ta2_yz_x_xyyy_1[i] * fe_0 + ta2_yz_xx_yyy_0[i] * fe_0 - ta2_yz_xx_yyy_1[i] * fe_0 + ta2_yz_xx_xyyy_0[i] * pa_x[i] - ta2_yz_xx_xyyy_1[i] * pc_x[i];

        ta2_yz_xxx_xyyz_0[i] = 2.0 * ta2_yz_x_xyyz_0[i] * fe_0 - 2.0 * ta2_yz_x_xyyz_1[i] * fe_0 + ta2_yz_xx_yyz_0[i] * fe_0 - ta2_yz_xx_yyz_1[i] * fe_0 + ta2_yz_xx_xyyz_0[i] * pa_x[i] - ta2_yz_xx_xyyz_1[i] * pc_x[i];

        ta2_yz_xxx_xyzz_0[i] = 2.0 * ta2_yz_x_xyzz_0[i] * fe_0 - 2.0 * ta2_yz_x_xyzz_1[i] * fe_0 + ta2_yz_xx_yzz_0[i] * fe_0 - ta2_yz_xx_yzz_1[i] * fe_0 + ta2_yz_xx_xyzz_0[i] * pa_x[i] - ta2_yz_xx_xyzz_1[i] * pc_x[i];

        ta2_yz_xxx_xzzz_0[i] = 2.0 * ta2_yz_x_xzzz_0[i] * fe_0 - 2.0 * ta2_yz_x_xzzz_1[i] * fe_0 + ta2_yz_xx_zzz_0[i] * fe_0 - ta2_yz_xx_zzz_1[i] * fe_0 + ta2_yz_xx_xzzz_0[i] * pa_x[i] - ta2_yz_xx_xzzz_1[i] * pc_x[i];

        ta2_yz_xxx_yyyy_0[i] = 2.0 * ta2_yz_x_yyyy_0[i] * fe_0 - 2.0 * ta2_yz_x_yyyy_1[i] * fe_0 + ta2_yz_xx_yyyy_0[i] * pa_x[i] - ta2_yz_xx_yyyy_1[i] * pc_x[i];

        ta2_yz_xxx_yyyz_0[i] = 2.0 * ta2_yz_x_yyyz_0[i] * fe_0 - 2.0 * ta2_yz_x_yyyz_1[i] * fe_0 + ta2_yz_xx_yyyz_0[i] * pa_x[i] - ta2_yz_xx_yyyz_1[i] * pc_x[i];

        ta2_yz_xxx_yyzz_0[i] = 2.0 * ta2_yz_x_yyzz_0[i] * fe_0 - 2.0 * ta2_yz_x_yyzz_1[i] * fe_0 + ta2_yz_xx_yyzz_0[i] * pa_x[i] - ta2_yz_xx_yyzz_1[i] * pc_x[i];

        ta2_yz_xxx_yzzz_0[i] = 2.0 * ta2_yz_x_yzzz_0[i] * fe_0 - 2.0 * ta2_yz_x_yzzz_1[i] * fe_0 + ta2_yz_xx_yzzz_0[i] * pa_x[i] - ta2_yz_xx_yzzz_1[i] * pc_x[i];

        ta2_yz_xxx_zzzz_0[i] = 2.0 * ta2_yz_x_zzzz_0[i] * fe_0 - 2.0 * ta2_yz_x_zzzz_1[i] * fe_0 + ta2_yz_xx_zzzz_0[i] * pa_x[i] - ta2_yz_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 615-630 components of targeted buffer : FG

    auto ta2_yz_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 615);

    auto ta2_yz_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 616);

    auto ta2_yz_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 617);

    auto ta2_yz_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 618);

    auto ta2_yz_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 619);

    auto ta2_yz_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 620);

    auto ta2_yz_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 621);

    auto ta2_yz_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 622);

    auto ta2_yz_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 623);

    auto ta2_yz_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 624);

    auto ta2_yz_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 625);

    auto ta2_yz_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 626);

    auto ta2_yz_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 627);

    auto ta2_yz_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 628);

    auto ta2_yz_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 629);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xx_xxxx_1, ta1_z_xx_xxxy_1, ta1_z_xx_xxxz_1, ta1_z_xx_xxyy_1, ta1_z_xx_xxyz_1, ta1_z_xx_xxzz_1, ta1_z_xx_xyyy_1, ta1_z_xx_xyyz_1, ta1_z_xx_xyzz_1, ta1_z_xx_xzzz_1, ta1_z_xx_zzzz_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxxx_0, ta2_yz_xx_xxxx_1, ta2_yz_xx_xxxy_0, ta2_yz_xx_xxxy_1, ta2_yz_xx_xxxz_0, ta2_yz_xx_xxxz_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxyy_0, ta2_yz_xx_xxyy_1, ta2_yz_xx_xxyz_0, ta2_yz_xx_xxyz_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xxzz_0, ta2_yz_xx_xxzz_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyyy_0, ta2_yz_xx_xyyy_1, ta2_yz_xx_xyyz_0, ta2_yz_xx_xyyz_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xyzz_0, ta2_yz_xx_xyzz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_xzzz_0, ta2_yz_xx_xzzz_1, ta2_yz_xx_zzzz_0, ta2_yz_xx_zzzz_1, ta2_yz_xxy_xxxx_0, ta2_yz_xxy_xxxy_0, ta2_yz_xxy_xxxz_0, ta2_yz_xxy_xxyy_0, ta2_yz_xxy_xxyz_0, ta2_yz_xxy_xxzz_0, ta2_yz_xxy_xyyy_0, ta2_yz_xxy_xyyz_0, ta2_yz_xxy_xyzz_0, ta2_yz_xxy_xzzz_0, ta2_yz_xxy_yyyy_0, ta2_yz_xxy_yyyz_0, ta2_yz_xxy_yyzz_0, ta2_yz_xxy_yzzz_0, ta2_yz_xxy_zzzz_0, ta2_yz_xy_yyyy_0, ta2_yz_xy_yyyy_1, ta2_yz_xy_yyyz_0, ta2_yz_xy_yyyz_1, ta2_yz_xy_yyzz_0, ta2_yz_xy_yyzz_1, ta2_yz_xy_yzzz_0, ta2_yz_xy_yzzz_1, ta2_yz_y_yyyy_0, ta2_yz_y_yyyy_1, ta2_yz_y_yyyz_0, ta2_yz_y_yyyz_1, ta2_yz_y_yyzz_0, ta2_yz_y_yyzz_1, ta2_yz_y_yzzz_0, ta2_yz_y_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxy_xxxx_0[i] = ta1_z_xx_xxxx_1[i] + ta2_yz_xx_xxxx_0[i] * pa_y[i] - ta2_yz_xx_xxxx_1[i] * pc_y[i];

        ta2_yz_xxy_xxxy_0[i] = ta2_yz_xx_xxx_0[i] * fe_0 - ta2_yz_xx_xxx_1[i] * fe_0 + ta1_z_xx_xxxy_1[i] + ta2_yz_xx_xxxy_0[i] * pa_y[i] - ta2_yz_xx_xxxy_1[i] * pc_y[i];

        ta2_yz_xxy_xxxz_0[i] = ta1_z_xx_xxxz_1[i] + ta2_yz_xx_xxxz_0[i] * pa_y[i] - ta2_yz_xx_xxxz_1[i] * pc_y[i];

        ta2_yz_xxy_xxyy_0[i] = 2.0 * ta2_yz_xx_xxy_0[i] * fe_0 - 2.0 * ta2_yz_xx_xxy_1[i] * fe_0 + ta1_z_xx_xxyy_1[i] + ta2_yz_xx_xxyy_0[i] * pa_y[i] - ta2_yz_xx_xxyy_1[i] * pc_y[i];

        ta2_yz_xxy_xxyz_0[i] = ta2_yz_xx_xxz_0[i] * fe_0 - ta2_yz_xx_xxz_1[i] * fe_0 + ta1_z_xx_xxyz_1[i] + ta2_yz_xx_xxyz_0[i] * pa_y[i] - ta2_yz_xx_xxyz_1[i] * pc_y[i];

        ta2_yz_xxy_xxzz_0[i] = ta1_z_xx_xxzz_1[i] + ta2_yz_xx_xxzz_0[i] * pa_y[i] - ta2_yz_xx_xxzz_1[i] * pc_y[i];

        ta2_yz_xxy_xyyy_0[i] = 3.0 * ta2_yz_xx_xyy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xyy_1[i] * fe_0 + ta1_z_xx_xyyy_1[i] + ta2_yz_xx_xyyy_0[i] * pa_y[i] - ta2_yz_xx_xyyy_1[i] * pc_y[i];

        ta2_yz_xxy_xyyz_0[i] = 2.0 * ta2_yz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xyz_1[i] * fe_0 + ta1_z_xx_xyyz_1[i] + ta2_yz_xx_xyyz_0[i] * pa_y[i] - ta2_yz_xx_xyyz_1[i] * pc_y[i];

        ta2_yz_xxy_xyzz_0[i] = ta2_yz_xx_xzz_0[i] * fe_0 - ta2_yz_xx_xzz_1[i] * fe_0 + ta1_z_xx_xyzz_1[i] + ta2_yz_xx_xyzz_0[i] * pa_y[i] - ta2_yz_xx_xyzz_1[i] * pc_y[i];

        ta2_yz_xxy_xzzz_0[i] = ta1_z_xx_xzzz_1[i] + ta2_yz_xx_xzzz_0[i] * pa_y[i] - ta2_yz_xx_xzzz_1[i] * pc_y[i];

        ta2_yz_xxy_yyyy_0[i] = ta2_yz_y_yyyy_0[i] * fe_0 - ta2_yz_y_yyyy_1[i] * fe_0 + ta2_yz_xy_yyyy_0[i] * pa_x[i] - ta2_yz_xy_yyyy_1[i] * pc_x[i];

        ta2_yz_xxy_yyyz_0[i] = ta2_yz_y_yyyz_0[i] * fe_0 - ta2_yz_y_yyyz_1[i] * fe_0 + ta2_yz_xy_yyyz_0[i] * pa_x[i] - ta2_yz_xy_yyyz_1[i] * pc_x[i];

        ta2_yz_xxy_yyzz_0[i] = ta2_yz_y_yyzz_0[i] * fe_0 - ta2_yz_y_yyzz_1[i] * fe_0 + ta2_yz_xy_yyzz_0[i] * pa_x[i] - ta2_yz_xy_yyzz_1[i] * pc_x[i];

        ta2_yz_xxy_yzzz_0[i] = ta2_yz_y_yzzz_0[i] * fe_0 - ta2_yz_y_yzzz_1[i] * fe_0 + ta2_yz_xy_yzzz_0[i] * pa_x[i] - ta2_yz_xy_yzzz_1[i] * pc_x[i];

        ta2_yz_xxy_zzzz_0[i] = ta1_z_xx_zzzz_1[i] + ta2_yz_xx_zzzz_0[i] * pa_y[i] - ta2_yz_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 630-645 components of targeted buffer : FG

    auto ta2_yz_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 630);

    auto ta2_yz_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 631);

    auto ta2_yz_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 632);

    auto ta2_yz_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 633);

    auto ta2_yz_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 634);

    auto ta2_yz_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 635);

    auto ta2_yz_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 636);

    auto ta2_yz_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 637);

    auto ta2_yz_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 638);

    auto ta2_yz_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 639);

    auto ta2_yz_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 640);

    auto ta2_yz_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 641);

    auto ta2_yz_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 642);

    auto ta2_yz_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 643);

    auto ta2_yz_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 644);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xx_xxxx_1, ta1_y_xx_xxxy_1, ta1_y_xx_xxxz_1, ta1_y_xx_xxyy_1, ta1_y_xx_xxyz_1, ta1_y_xx_xxzz_1, ta1_y_xx_xyyy_1, ta1_y_xx_xyyz_1, ta1_y_xx_xyzz_1, ta1_y_xx_xzzz_1, ta1_y_xx_yyyy_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxxx_0, ta2_yz_xx_xxxx_1, ta2_yz_xx_xxxy_0, ta2_yz_xx_xxxy_1, ta2_yz_xx_xxxz_0, ta2_yz_xx_xxxz_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxyy_0, ta2_yz_xx_xxyy_1, ta2_yz_xx_xxyz_0, ta2_yz_xx_xxyz_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xxzz_0, ta2_yz_xx_xxzz_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyyy_0, ta2_yz_xx_xyyy_1, ta2_yz_xx_xyyz_0, ta2_yz_xx_xyyz_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xyzz_0, ta2_yz_xx_xyzz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_xzzz_0, ta2_yz_xx_xzzz_1, ta2_yz_xx_yyyy_0, ta2_yz_xx_yyyy_1, ta2_yz_xxz_xxxx_0, ta2_yz_xxz_xxxy_0, ta2_yz_xxz_xxxz_0, ta2_yz_xxz_xxyy_0, ta2_yz_xxz_xxyz_0, ta2_yz_xxz_xxzz_0, ta2_yz_xxz_xyyy_0, ta2_yz_xxz_xyyz_0, ta2_yz_xxz_xyzz_0, ta2_yz_xxz_xzzz_0, ta2_yz_xxz_yyyy_0, ta2_yz_xxz_yyyz_0, ta2_yz_xxz_yyzz_0, ta2_yz_xxz_yzzz_0, ta2_yz_xxz_zzzz_0, ta2_yz_xz_yyyz_0, ta2_yz_xz_yyyz_1, ta2_yz_xz_yyzz_0, ta2_yz_xz_yyzz_1, ta2_yz_xz_yzzz_0, ta2_yz_xz_yzzz_1, ta2_yz_xz_zzzz_0, ta2_yz_xz_zzzz_1, ta2_yz_z_yyyz_0, ta2_yz_z_yyyz_1, ta2_yz_z_yyzz_0, ta2_yz_z_yyzz_1, ta2_yz_z_yzzz_0, ta2_yz_z_yzzz_1, ta2_yz_z_zzzz_0, ta2_yz_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxz_xxxx_0[i] = ta1_y_xx_xxxx_1[i] + ta2_yz_xx_xxxx_0[i] * pa_z[i] - ta2_yz_xx_xxxx_1[i] * pc_z[i];

        ta2_yz_xxz_xxxy_0[i] = ta1_y_xx_xxxy_1[i] + ta2_yz_xx_xxxy_0[i] * pa_z[i] - ta2_yz_xx_xxxy_1[i] * pc_z[i];

        ta2_yz_xxz_xxxz_0[i] = ta2_yz_xx_xxx_0[i] * fe_0 - ta2_yz_xx_xxx_1[i] * fe_0 + ta1_y_xx_xxxz_1[i] + ta2_yz_xx_xxxz_0[i] * pa_z[i] - ta2_yz_xx_xxxz_1[i] * pc_z[i];

        ta2_yz_xxz_xxyy_0[i] = ta1_y_xx_xxyy_1[i] + ta2_yz_xx_xxyy_0[i] * pa_z[i] - ta2_yz_xx_xxyy_1[i] * pc_z[i];

        ta2_yz_xxz_xxyz_0[i] = ta2_yz_xx_xxy_0[i] * fe_0 - ta2_yz_xx_xxy_1[i] * fe_0 + ta1_y_xx_xxyz_1[i] + ta2_yz_xx_xxyz_0[i] * pa_z[i] - ta2_yz_xx_xxyz_1[i] * pc_z[i];

        ta2_yz_xxz_xxzz_0[i] = 2.0 * ta2_yz_xx_xxz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xxz_1[i] * fe_0 + ta1_y_xx_xxzz_1[i] + ta2_yz_xx_xxzz_0[i] * pa_z[i] - ta2_yz_xx_xxzz_1[i] * pc_z[i];

        ta2_yz_xxz_xyyy_0[i] = ta1_y_xx_xyyy_1[i] + ta2_yz_xx_xyyy_0[i] * pa_z[i] - ta2_yz_xx_xyyy_1[i] * pc_z[i];

        ta2_yz_xxz_xyyz_0[i] = ta2_yz_xx_xyy_0[i] * fe_0 - ta2_yz_xx_xyy_1[i] * fe_0 + ta1_y_xx_xyyz_1[i] + ta2_yz_xx_xyyz_0[i] * pa_z[i] - ta2_yz_xx_xyyz_1[i] * pc_z[i];

        ta2_yz_xxz_xyzz_0[i] = 2.0 * ta2_yz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xyz_1[i] * fe_0 + ta1_y_xx_xyzz_1[i] + ta2_yz_xx_xyzz_0[i] * pa_z[i] - ta2_yz_xx_xyzz_1[i] * pc_z[i];

        ta2_yz_xxz_xzzz_0[i] = 3.0 * ta2_yz_xx_xzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xzz_1[i] * fe_0 + ta1_y_xx_xzzz_1[i] + ta2_yz_xx_xzzz_0[i] * pa_z[i] - ta2_yz_xx_xzzz_1[i] * pc_z[i];

        ta2_yz_xxz_yyyy_0[i] = ta1_y_xx_yyyy_1[i] + ta2_yz_xx_yyyy_0[i] * pa_z[i] - ta2_yz_xx_yyyy_1[i] * pc_z[i];

        ta2_yz_xxz_yyyz_0[i] = ta2_yz_z_yyyz_0[i] * fe_0 - ta2_yz_z_yyyz_1[i] * fe_0 + ta2_yz_xz_yyyz_0[i] * pa_x[i] - ta2_yz_xz_yyyz_1[i] * pc_x[i];

        ta2_yz_xxz_yyzz_0[i] = ta2_yz_z_yyzz_0[i] * fe_0 - ta2_yz_z_yyzz_1[i] * fe_0 + ta2_yz_xz_yyzz_0[i] * pa_x[i] - ta2_yz_xz_yyzz_1[i] * pc_x[i];

        ta2_yz_xxz_yzzz_0[i] = ta2_yz_z_yzzz_0[i] * fe_0 - ta2_yz_z_yzzz_1[i] * fe_0 + ta2_yz_xz_yzzz_0[i] * pa_x[i] - ta2_yz_xz_yzzz_1[i] * pc_x[i];

        ta2_yz_xxz_zzzz_0[i] = ta2_yz_z_zzzz_0[i] * fe_0 - ta2_yz_z_zzzz_1[i] * fe_0 + ta2_yz_xz_zzzz_0[i] * pa_x[i] - ta2_yz_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 645-660 components of targeted buffer : FG

    auto ta2_yz_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 645);

    auto ta2_yz_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 646);

    auto ta2_yz_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 647);

    auto ta2_yz_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 648);

    auto ta2_yz_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 649);

    auto ta2_yz_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 650);

    auto ta2_yz_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 651);

    auto ta2_yz_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 652);

    auto ta2_yz_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 653);

    auto ta2_yz_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 654);

    auto ta2_yz_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 655);

    auto ta2_yz_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 656);

    auto ta2_yz_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 657);

    auto ta2_yz_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 658);

    auto ta2_yz_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 659);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xyy_xxxx_0, ta2_yz_xyy_xxxy_0, ta2_yz_xyy_xxxz_0, ta2_yz_xyy_xxyy_0, ta2_yz_xyy_xxyz_0, ta2_yz_xyy_xxzz_0, ta2_yz_xyy_xyyy_0, ta2_yz_xyy_xyyz_0, ta2_yz_xyy_xyzz_0, ta2_yz_xyy_xzzz_0, ta2_yz_xyy_yyyy_0, ta2_yz_xyy_yyyz_0, ta2_yz_xyy_yyzz_0, ta2_yz_xyy_yzzz_0, ta2_yz_xyy_zzzz_0, ta2_yz_yy_xxx_0, ta2_yz_yy_xxx_1, ta2_yz_yy_xxxx_0, ta2_yz_yy_xxxx_1, ta2_yz_yy_xxxy_0, ta2_yz_yy_xxxy_1, ta2_yz_yy_xxxz_0, ta2_yz_yy_xxxz_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xxyy_0, ta2_yz_yy_xxyy_1, ta2_yz_yy_xxyz_0, ta2_yz_yy_xxyz_1, ta2_yz_yy_xxz_0, ta2_yz_yy_xxz_1, ta2_yz_yy_xxzz_0, ta2_yz_yy_xxzz_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyyy_0, ta2_yz_yy_xyyy_1, ta2_yz_yy_xyyz_0, ta2_yz_yy_xyyz_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_xyzz_0, ta2_yz_yy_xyzz_1, ta2_yz_yy_xzz_0, ta2_yz_yy_xzz_1, ta2_yz_yy_xzzz_0, ta2_yz_yy_xzzz_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyyy_0, ta2_yz_yy_yyyy_1, ta2_yz_yy_yyyz_0, ta2_yz_yy_yyyz_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yyzz_0, ta2_yz_yy_yyzz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_yzzz_0, ta2_yz_yy_yzzz_1, ta2_yz_yy_zzz_0, ta2_yz_yy_zzz_1, ta2_yz_yy_zzzz_0, ta2_yz_yy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyy_xxxx_0[i] = 4.0 * ta2_yz_yy_xxx_0[i] * fe_0 - 4.0 * ta2_yz_yy_xxx_1[i] * fe_0 + ta2_yz_yy_xxxx_0[i] * pa_x[i] - ta2_yz_yy_xxxx_1[i] * pc_x[i];

        ta2_yz_xyy_xxxy_0[i] = 3.0 * ta2_yz_yy_xxy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxy_1[i] * fe_0 + ta2_yz_yy_xxxy_0[i] * pa_x[i] - ta2_yz_yy_xxxy_1[i] * pc_x[i];

        ta2_yz_xyy_xxxz_0[i] = 3.0 * ta2_yz_yy_xxz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxz_1[i] * fe_0 + ta2_yz_yy_xxxz_0[i] * pa_x[i] - ta2_yz_yy_xxxz_1[i] * pc_x[i];

        ta2_yz_xyy_xxyy_0[i] = 2.0 * ta2_yz_yy_xyy_0[i] * fe_0 - 2.0 * ta2_yz_yy_xyy_1[i] * fe_0 + ta2_yz_yy_xxyy_0[i] * pa_x[i] - ta2_yz_yy_xxyy_1[i] * pc_x[i];

        ta2_yz_xyy_xxyz_0[i] = 2.0 * ta2_yz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yy_xyz_1[i] * fe_0 + ta2_yz_yy_xxyz_0[i] * pa_x[i] - ta2_yz_yy_xxyz_1[i] * pc_x[i];

        ta2_yz_xyy_xxzz_0[i] = 2.0 * ta2_yz_yy_xzz_0[i] * fe_0 - 2.0 * ta2_yz_yy_xzz_1[i] * fe_0 + ta2_yz_yy_xxzz_0[i] * pa_x[i] - ta2_yz_yy_xxzz_1[i] * pc_x[i];

        ta2_yz_xyy_xyyy_0[i] = ta2_yz_yy_yyy_0[i] * fe_0 - ta2_yz_yy_yyy_1[i] * fe_0 + ta2_yz_yy_xyyy_0[i] * pa_x[i] - ta2_yz_yy_xyyy_1[i] * pc_x[i];

        ta2_yz_xyy_xyyz_0[i] = ta2_yz_yy_yyz_0[i] * fe_0 - ta2_yz_yy_yyz_1[i] * fe_0 + ta2_yz_yy_xyyz_0[i] * pa_x[i] - ta2_yz_yy_xyyz_1[i] * pc_x[i];

        ta2_yz_xyy_xyzz_0[i] = ta2_yz_yy_yzz_0[i] * fe_0 - ta2_yz_yy_yzz_1[i] * fe_0 + ta2_yz_yy_xyzz_0[i] * pa_x[i] - ta2_yz_yy_xyzz_1[i] * pc_x[i];

        ta2_yz_xyy_xzzz_0[i] = ta2_yz_yy_zzz_0[i] * fe_0 - ta2_yz_yy_zzz_1[i] * fe_0 + ta2_yz_yy_xzzz_0[i] * pa_x[i] - ta2_yz_yy_xzzz_1[i] * pc_x[i];

        ta2_yz_xyy_yyyy_0[i] = ta2_yz_yy_yyyy_0[i] * pa_x[i] - ta2_yz_yy_yyyy_1[i] * pc_x[i];

        ta2_yz_xyy_yyyz_0[i] = ta2_yz_yy_yyyz_0[i] * pa_x[i] - ta2_yz_yy_yyyz_1[i] * pc_x[i];

        ta2_yz_xyy_yyzz_0[i] = ta2_yz_yy_yyzz_0[i] * pa_x[i] - ta2_yz_yy_yyzz_1[i] * pc_x[i];

        ta2_yz_xyy_yzzz_0[i] = ta2_yz_yy_yzzz_0[i] * pa_x[i] - ta2_yz_yy_yzzz_1[i] * pc_x[i];

        ta2_yz_xyy_zzzz_0[i] = ta2_yz_yy_zzzz_0[i] * pa_x[i] - ta2_yz_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 660-675 components of targeted buffer : FG

    auto ta2_yz_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 660);

    auto ta2_yz_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 661);

    auto ta2_yz_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 662);

    auto ta2_yz_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 663);

    auto ta2_yz_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 664);

    auto ta2_yz_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 665);

    auto ta2_yz_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 666);

    auto ta2_yz_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 667);

    auto ta2_yz_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 668);

    auto ta2_yz_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 669);

    auto ta2_yz_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 670);

    auto ta2_yz_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 671);

    auto ta2_yz_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 672);

    auto ta2_yz_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 673);

    auto ta2_yz_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 674);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xy_xxxy_1, ta1_y_xy_xxyy_1, ta1_y_xy_xyyy_1, ta1_z_xz_xxxx_1, ta1_z_xz_xxxz_1, ta1_z_xz_xxzz_1, ta1_z_xz_xzzz_1, ta2_yz_xy_xxxy_0, ta2_yz_xy_xxxy_1, ta2_yz_xy_xxyy_0, ta2_yz_xy_xxyy_1, ta2_yz_xy_xyyy_0, ta2_yz_xy_xyyy_1, ta2_yz_xyz_xxxx_0, ta2_yz_xyz_xxxy_0, ta2_yz_xyz_xxxz_0, ta2_yz_xyz_xxyy_0, ta2_yz_xyz_xxyz_0, ta2_yz_xyz_xxzz_0, ta2_yz_xyz_xyyy_0, ta2_yz_xyz_xyyz_0, ta2_yz_xyz_xyzz_0, ta2_yz_xyz_xzzz_0, ta2_yz_xyz_yyyy_0, ta2_yz_xyz_yyyz_0, ta2_yz_xyz_yyzz_0, ta2_yz_xyz_yzzz_0, ta2_yz_xyz_zzzz_0, ta2_yz_xz_xxxx_0, ta2_yz_xz_xxxx_1, ta2_yz_xz_xxxz_0, ta2_yz_xz_xxxz_1, ta2_yz_xz_xxzz_0, ta2_yz_xz_xxzz_1, ta2_yz_xz_xzzz_0, ta2_yz_xz_xzzz_1, ta2_yz_yz_xxyz_0, ta2_yz_yz_xxyz_1, ta2_yz_yz_xyyz_0, ta2_yz_yz_xyyz_1, ta2_yz_yz_xyz_0, ta2_yz_yz_xyz_1, ta2_yz_yz_xyzz_0, ta2_yz_yz_xyzz_1, ta2_yz_yz_yyyy_0, ta2_yz_yz_yyyy_1, ta2_yz_yz_yyyz_0, ta2_yz_yz_yyyz_1, ta2_yz_yz_yyz_0, ta2_yz_yz_yyz_1, ta2_yz_yz_yyzz_0, ta2_yz_yz_yyzz_1, ta2_yz_yz_yzz_0, ta2_yz_yz_yzz_1, ta2_yz_yz_yzzz_0, ta2_yz_yz_yzzz_1, ta2_yz_yz_zzzz_0, ta2_yz_yz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyz_xxxx_0[i] = ta1_z_xz_xxxx_1[i] + ta2_yz_xz_xxxx_0[i] * pa_y[i] - ta2_yz_xz_xxxx_1[i] * pc_y[i];

        ta2_yz_xyz_xxxy_0[i] = ta1_y_xy_xxxy_1[i] + ta2_yz_xy_xxxy_0[i] * pa_z[i] - ta2_yz_xy_xxxy_1[i] * pc_z[i];

        ta2_yz_xyz_xxxz_0[i] = ta1_z_xz_xxxz_1[i] + ta2_yz_xz_xxxz_0[i] * pa_y[i] - ta2_yz_xz_xxxz_1[i] * pc_y[i];

        ta2_yz_xyz_xxyy_0[i] = ta1_y_xy_xxyy_1[i] + ta2_yz_xy_xxyy_0[i] * pa_z[i] - ta2_yz_xy_xxyy_1[i] * pc_z[i];

        ta2_yz_xyz_xxyz_0[i] = 2.0 * ta2_yz_yz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xyz_1[i] * fe_0 + ta2_yz_yz_xxyz_0[i] * pa_x[i] - ta2_yz_yz_xxyz_1[i] * pc_x[i];

        ta2_yz_xyz_xxzz_0[i] = ta1_z_xz_xxzz_1[i] + ta2_yz_xz_xxzz_0[i] * pa_y[i] - ta2_yz_xz_xxzz_1[i] * pc_y[i];

        ta2_yz_xyz_xyyy_0[i] = ta1_y_xy_xyyy_1[i] + ta2_yz_xy_xyyy_0[i] * pa_z[i] - ta2_yz_xy_xyyy_1[i] * pc_z[i];

        ta2_yz_xyz_xyyz_0[i] = ta2_yz_yz_yyz_0[i] * fe_0 - ta2_yz_yz_yyz_1[i] * fe_0 + ta2_yz_yz_xyyz_0[i] * pa_x[i] - ta2_yz_yz_xyyz_1[i] * pc_x[i];

        ta2_yz_xyz_xyzz_0[i] = ta2_yz_yz_yzz_0[i] * fe_0 - ta2_yz_yz_yzz_1[i] * fe_0 + ta2_yz_yz_xyzz_0[i] * pa_x[i] - ta2_yz_yz_xyzz_1[i] * pc_x[i];

        ta2_yz_xyz_xzzz_0[i] = ta1_z_xz_xzzz_1[i] + ta2_yz_xz_xzzz_0[i] * pa_y[i] - ta2_yz_xz_xzzz_1[i] * pc_y[i];

        ta2_yz_xyz_yyyy_0[i] = ta2_yz_yz_yyyy_0[i] * pa_x[i] - ta2_yz_yz_yyyy_1[i] * pc_x[i];

        ta2_yz_xyz_yyyz_0[i] = ta2_yz_yz_yyyz_0[i] * pa_x[i] - ta2_yz_yz_yyyz_1[i] * pc_x[i];

        ta2_yz_xyz_yyzz_0[i] = ta2_yz_yz_yyzz_0[i] * pa_x[i] - ta2_yz_yz_yyzz_1[i] * pc_x[i];

        ta2_yz_xyz_yzzz_0[i] = ta2_yz_yz_yzzz_0[i] * pa_x[i] - ta2_yz_yz_yzzz_1[i] * pc_x[i];

        ta2_yz_xyz_zzzz_0[i] = ta2_yz_yz_zzzz_0[i] * pa_x[i] - ta2_yz_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 675-690 components of targeted buffer : FG

    auto ta2_yz_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 675);

    auto ta2_yz_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 676);

    auto ta2_yz_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 677);

    auto ta2_yz_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 678);

    auto ta2_yz_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 679);

    auto ta2_yz_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 680);

    auto ta2_yz_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 681);

    auto ta2_yz_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 682);

    auto ta2_yz_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 683);

    auto ta2_yz_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 684);

    auto ta2_yz_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 685);

    auto ta2_yz_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 686);

    auto ta2_yz_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 687);

    auto ta2_yz_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 688);

    auto ta2_yz_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 689);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xzz_xxxx_0, ta2_yz_xzz_xxxy_0, ta2_yz_xzz_xxxz_0, ta2_yz_xzz_xxyy_0, ta2_yz_xzz_xxyz_0, ta2_yz_xzz_xxzz_0, ta2_yz_xzz_xyyy_0, ta2_yz_xzz_xyyz_0, ta2_yz_xzz_xyzz_0, ta2_yz_xzz_xzzz_0, ta2_yz_xzz_yyyy_0, ta2_yz_xzz_yyyz_0, ta2_yz_xzz_yyzz_0, ta2_yz_xzz_yzzz_0, ta2_yz_xzz_zzzz_0, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxxx_0, ta2_yz_zz_xxxx_1, ta2_yz_zz_xxxy_0, ta2_yz_zz_xxxy_1, ta2_yz_zz_xxxz_0, ta2_yz_zz_xxxz_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxyy_0, ta2_yz_zz_xxyy_1, ta2_yz_zz_xxyz_0, ta2_yz_zz_xxyz_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xxzz_0, ta2_yz_zz_xxzz_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyyy_0, ta2_yz_zz_xyyy_1, ta2_yz_zz_xyyz_0, ta2_yz_zz_xyyz_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xyzz_0, ta2_yz_zz_xyzz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_xzzz_0, ta2_yz_zz_xzzz_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyyy_0, ta2_yz_zz_yyyy_1, ta2_yz_zz_yyyz_0, ta2_yz_zz_yyyz_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yyzz_0, ta2_yz_zz_yyzz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_yzzz_0, ta2_yz_zz_yzzz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, ta2_yz_zz_zzzz_0, ta2_yz_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzz_xxxx_0[i] = 4.0 * ta2_yz_zz_xxx_0[i] * fe_0 - 4.0 * ta2_yz_zz_xxx_1[i] * fe_0 + ta2_yz_zz_xxxx_0[i] * pa_x[i] - ta2_yz_zz_xxxx_1[i] * pc_x[i];

        ta2_yz_xzz_xxxy_0[i] = 3.0 * ta2_yz_zz_xxy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxy_1[i] * fe_0 + ta2_yz_zz_xxxy_0[i] * pa_x[i] - ta2_yz_zz_xxxy_1[i] * pc_x[i];

        ta2_yz_xzz_xxxz_0[i] = 3.0 * ta2_yz_zz_xxz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxz_1[i] * fe_0 + ta2_yz_zz_xxxz_0[i] * pa_x[i] - ta2_yz_zz_xxxz_1[i] * pc_x[i];

        ta2_yz_xzz_xxyy_0[i] = 2.0 * ta2_yz_zz_xyy_0[i] * fe_0 - 2.0 * ta2_yz_zz_xyy_1[i] * fe_0 + ta2_yz_zz_xxyy_0[i] * pa_x[i] - ta2_yz_zz_xxyy_1[i] * pc_x[i];

        ta2_yz_xzz_xxyz_0[i] = 2.0 * ta2_yz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xyz_1[i] * fe_0 + ta2_yz_zz_xxyz_0[i] * pa_x[i] - ta2_yz_zz_xxyz_1[i] * pc_x[i];

        ta2_yz_xzz_xxzz_0[i] = 2.0 * ta2_yz_zz_xzz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xzz_1[i] * fe_0 + ta2_yz_zz_xxzz_0[i] * pa_x[i] - ta2_yz_zz_xxzz_1[i] * pc_x[i];

        ta2_yz_xzz_xyyy_0[i] = ta2_yz_zz_yyy_0[i] * fe_0 - ta2_yz_zz_yyy_1[i] * fe_0 + ta2_yz_zz_xyyy_0[i] * pa_x[i] - ta2_yz_zz_xyyy_1[i] * pc_x[i];

        ta2_yz_xzz_xyyz_0[i] = ta2_yz_zz_yyz_0[i] * fe_0 - ta2_yz_zz_yyz_1[i] * fe_0 + ta2_yz_zz_xyyz_0[i] * pa_x[i] - ta2_yz_zz_xyyz_1[i] * pc_x[i];

        ta2_yz_xzz_xyzz_0[i] = ta2_yz_zz_yzz_0[i] * fe_0 - ta2_yz_zz_yzz_1[i] * fe_0 + ta2_yz_zz_xyzz_0[i] * pa_x[i] - ta2_yz_zz_xyzz_1[i] * pc_x[i];

        ta2_yz_xzz_xzzz_0[i] = ta2_yz_zz_zzz_0[i] * fe_0 - ta2_yz_zz_zzz_1[i] * fe_0 + ta2_yz_zz_xzzz_0[i] * pa_x[i] - ta2_yz_zz_xzzz_1[i] * pc_x[i];

        ta2_yz_xzz_yyyy_0[i] = ta2_yz_zz_yyyy_0[i] * pa_x[i] - ta2_yz_zz_yyyy_1[i] * pc_x[i];

        ta2_yz_xzz_yyyz_0[i] = ta2_yz_zz_yyyz_0[i] * pa_x[i] - ta2_yz_zz_yyyz_1[i] * pc_x[i];

        ta2_yz_xzz_yyzz_0[i] = ta2_yz_zz_yyzz_0[i] * pa_x[i] - ta2_yz_zz_yyzz_1[i] * pc_x[i];

        ta2_yz_xzz_yzzz_0[i] = ta2_yz_zz_yzzz_0[i] * pa_x[i] - ta2_yz_zz_yzzz_1[i] * pc_x[i];

        ta2_yz_xzz_zzzz_0[i] = ta2_yz_zz_zzzz_0[i] * pa_x[i] - ta2_yz_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 690-705 components of targeted buffer : FG

    auto ta2_yz_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 690);

    auto ta2_yz_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 691);

    auto ta2_yz_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 692);

    auto ta2_yz_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 693);

    auto ta2_yz_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 694);

    auto ta2_yz_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 695);

    auto ta2_yz_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 696);

    auto ta2_yz_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 697);

    auto ta2_yz_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 698);

    auto ta2_yz_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 699);

    auto ta2_yz_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 700);

    auto ta2_yz_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 701);

    auto ta2_yz_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 702);

    auto ta2_yz_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 703);

    auto ta2_yz_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 704);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yy_xxxx_1, ta1_z_yy_xxxy_1, ta1_z_yy_xxxz_1, ta1_z_yy_xxyy_1, ta1_z_yy_xxyz_1, ta1_z_yy_xxzz_1, ta1_z_yy_xyyy_1, ta1_z_yy_xyyz_1, ta1_z_yy_xyzz_1, ta1_z_yy_xzzz_1, ta1_z_yy_yyyy_1, ta1_z_yy_yyyz_1, ta1_z_yy_yyzz_1, ta1_z_yy_yzzz_1, ta1_z_yy_zzzz_1, ta2_yz_y_xxxx_0, ta2_yz_y_xxxx_1, ta2_yz_y_xxxy_0, ta2_yz_y_xxxy_1, ta2_yz_y_xxxz_0, ta2_yz_y_xxxz_1, ta2_yz_y_xxyy_0, ta2_yz_y_xxyy_1, ta2_yz_y_xxyz_0, ta2_yz_y_xxyz_1, ta2_yz_y_xxzz_0, ta2_yz_y_xxzz_1, ta2_yz_y_xyyy_0, ta2_yz_y_xyyy_1, ta2_yz_y_xyyz_0, ta2_yz_y_xyyz_1, ta2_yz_y_xyzz_0, ta2_yz_y_xyzz_1, ta2_yz_y_xzzz_0, ta2_yz_y_xzzz_1, ta2_yz_y_yyyy_0, ta2_yz_y_yyyy_1, ta2_yz_y_yyyz_0, ta2_yz_y_yyyz_1, ta2_yz_y_yyzz_0, ta2_yz_y_yyzz_1, ta2_yz_y_yzzz_0, ta2_yz_y_yzzz_1, ta2_yz_y_zzzz_0, ta2_yz_y_zzzz_1, ta2_yz_yy_xxx_0, ta2_yz_yy_xxx_1, ta2_yz_yy_xxxx_0, ta2_yz_yy_xxxx_1, ta2_yz_yy_xxxy_0, ta2_yz_yy_xxxy_1, ta2_yz_yy_xxxz_0, ta2_yz_yy_xxxz_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xxyy_0, ta2_yz_yy_xxyy_1, ta2_yz_yy_xxyz_0, ta2_yz_yy_xxyz_1, ta2_yz_yy_xxz_0, ta2_yz_yy_xxz_1, ta2_yz_yy_xxzz_0, ta2_yz_yy_xxzz_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyyy_0, ta2_yz_yy_xyyy_1, ta2_yz_yy_xyyz_0, ta2_yz_yy_xyyz_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_xyzz_0, ta2_yz_yy_xyzz_1, ta2_yz_yy_xzz_0, ta2_yz_yy_xzz_1, ta2_yz_yy_xzzz_0, ta2_yz_yy_xzzz_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyyy_0, ta2_yz_yy_yyyy_1, ta2_yz_yy_yyyz_0, ta2_yz_yy_yyyz_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yyzz_0, ta2_yz_yy_yyzz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_yzzz_0, ta2_yz_yy_yzzz_1, ta2_yz_yy_zzz_0, ta2_yz_yy_zzz_1, ta2_yz_yy_zzzz_0, ta2_yz_yy_zzzz_1, ta2_yz_yyy_xxxx_0, ta2_yz_yyy_xxxy_0, ta2_yz_yyy_xxxz_0, ta2_yz_yyy_xxyy_0, ta2_yz_yyy_xxyz_0, ta2_yz_yyy_xxzz_0, ta2_yz_yyy_xyyy_0, ta2_yz_yyy_xyyz_0, ta2_yz_yyy_xyzz_0, ta2_yz_yyy_xzzz_0, ta2_yz_yyy_yyyy_0, ta2_yz_yyy_yyyz_0, ta2_yz_yyy_yyzz_0, ta2_yz_yyy_yzzz_0, ta2_yz_yyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyy_xxxx_0[i] = 2.0 * ta2_yz_y_xxxx_0[i] * fe_0 - 2.0 * ta2_yz_y_xxxx_1[i] * fe_0 + ta1_z_yy_xxxx_1[i] + ta2_yz_yy_xxxx_0[i] * pa_y[i] - ta2_yz_yy_xxxx_1[i] * pc_y[i];

        ta2_yz_yyy_xxxy_0[i] = 2.0 * ta2_yz_y_xxxy_0[i] * fe_0 - 2.0 * ta2_yz_y_xxxy_1[i] * fe_0 + ta2_yz_yy_xxx_0[i] * fe_0 - ta2_yz_yy_xxx_1[i] * fe_0 + ta1_z_yy_xxxy_1[i] + ta2_yz_yy_xxxy_0[i] * pa_y[i] - ta2_yz_yy_xxxy_1[i] * pc_y[i];

        ta2_yz_yyy_xxxz_0[i] = 2.0 * ta2_yz_y_xxxz_0[i] * fe_0 - 2.0 * ta2_yz_y_xxxz_1[i] * fe_0 + ta1_z_yy_xxxz_1[i] + ta2_yz_yy_xxxz_0[i] * pa_y[i] - ta2_yz_yy_xxxz_1[i] * pc_y[i];

        ta2_yz_yyy_xxyy_0[i] = 2.0 * ta2_yz_y_xxyy_0[i] * fe_0 - 2.0 * ta2_yz_y_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_yy_xxy_0[i] * fe_0 - 2.0 * ta2_yz_yy_xxy_1[i] * fe_0 + ta1_z_yy_xxyy_1[i] + ta2_yz_yy_xxyy_0[i] * pa_y[i] - ta2_yz_yy_xxyy_1[i] * pc_y[i];

        ta2_yz_yyy_xxyz_0[i] = 2.0 * ta2_yz_y_xxyz_0[i] * fe_0 - 2.0 * ta2_yz_y_xxyz_1[i] * fe_0 + ta2_yz_yy_xxz_0[i] * fe_0 - ta2_yz_yy_xxz_1[i] * fe_0 + ta1_z_yy_xxyz_1[i] + ta2_yz_yy_xxyz_0[i] * pa_y[i] - ta2_yz_yy_xxyz_1[i] * pc_y[i];

        ta2_yz_yyy_xxzz_0[i] = 2.0 * ta2_yz_y_xxzz_0[i] * fe_0 - 2.0 * ta2_yz_y_xxzz_1[i] * fe_0 + ta1_z_yy_xxzz_1[i] + ta2_yz_yy_xxzz_0[i] * pa_y[i] - ta2_yz_yy_xxzz_1[i] * pc_y[i];

        ta2_yz_yyy_xyyy_0[i] = 2.0 * ta2_yz_y_xyyy_0[i] * fe_0 - 2.0 * ta2_yz_y_xyyy_1[i] * fe_0 + 3.0 * ta2_yz_yy_xyy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xyy_1[i] * fe_0 + ta1_z_yy_xyyy_1[i] + ta2_yz_yy_xyyy_0[i] * pa_y[i] - ta2_yz_yy_xyyy_1[i] * pc_y[i];

        ta2_yz_yyy_xyyz_0[i] = 2.0 * ta2_yz_y_xyyz_0[i] * fe_0 - 2.0 * ta2_yz_y_xyyz_1[i] * fe_0 + 2.0 * ta2_yz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yy_xyz_1[i] * fe_0 + ta1_z_yy_xyyz_1[i] + ta2_yz_yy_xyyz_0[i] * pa_y[i] - ta2_yz_yy_xyyz_1[i] * pc_y[i];

        ta2_yz_yyy_xyzz_0[i] = 2.0 * ta2_yz_y_xyzz_0[i] * fe_0 - 2.0 * ta2_yz_y_xyzz_1[i] * fe_0 + ta2_yz_yy_xzz_0[i] * fe_0 - ta2_yz_yy_xzz_1[i] * fe_0 + ta1_z_yy_xyzz_1[i] + ta2_yz_yy_xyzz_0[i] * pa_y[i] - ta2_yz_yy_xyzz_1[i] * pc_y[i];

        ta2_yz_yyy_xzzz_0[i] = 2.0 * ta2_yz_y_xzzz_0[i] * fe_0 - 2.0 * ta2_yz_y_xzzz_1[i] * fe_0 + ta1_z_yy_xzzz_1[i] + ta2_yz_yy_xzzz_0[i] * pa_y[i] - ta2_yz_yy_xzzz_1[i] * pc_y[i];

        ta2_yz_yyy_yyyy_0[i] = 2.0 * ta2_yz_y_yyyy_0[i] * fe_0 - 2.0 * ta2_yz_y_yyyy_1[i] * fe_0 + 4.0 * ta2_yz_yy_yyy_0[i] * fe_0 - 4.0 * ta2_yz_yy_yyy_1[i] * fe_0 + ta1_z_yy_yyyy_1[i] + ta2_yz_yy_yyyy_0[i] * pa_y[i] - ta2_yz_yy_yyyy_1[i] * pc_y[i];

        ta2_yz_yyy_yyyz_0[i] = 2.0 * ta2_yz_y_yyyz_0[i] * fe_0 - 2.0 * ta2_yz_y_yyyz_1[i] * fe_0 + 3.0 * ta2_yz_yy_yyz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yyz_1[i] * fe_0 + ta1_z_yy_yyyz_1[i] + ta2_yz_yy_yyyz_0[i] * pa_y[i] - ta2_yz_yy_yyyz_1[i] * pc_y[i];

        ta2_yz_yyy_yyzz_0[i] = 2.0 * ta2_yz_y_yyzz_0[i] * fe_0 - 2.0 * ta2_yz_y_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_yy_yzz_0[i] * fe_0 - 2.0 * ta2_yz_yy_yzz_1[i] * fe_0 + ta1_z_yy_yyzz_1[i] + ta2_yz_yy_yyzz_0[i] * pa_y[i] - ta2_yz_yy_yyzz_1[i] * pc_y[i];

        ta2_yz_yyy_yzzz_0[i] = 2.0 * ta2_yz_y_yzzz_0[i] * fe_0 - 2.0 * ta2_yz_y_yzzz_1[i] * fe_0 + ta2_yz_yy_zzz_0[i] * fe_0 - ta2_yz_yy_zzz_1[i] * fe_0 + ta1_z_yy_yzzz_1[i] + ta2_yz_yy_yzzz_0[i] * pa_y[i] - ta2_yz_yy_yzzz_1[i] * pc_y[i];

        ta2_yz_yyy_zzzz_0[i] = 2.0 * ta2_yz_y_zzzz_0[i] * fe_0 - 2.0 * ta2_yz_y_zzzz_1[i] * fe_0 + ta1_z_yy_zzzz_1[i] + ta2_yz_yy_zzzz_0[i] * pa_y[i] - ta2_yz_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 705-720 components of targeted buffer : FG

    auto ta2_yz_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 705);

    auto ta2_yz_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 706);

    auto ta2_yz_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 707);

    auto ta2_yz_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 708);

    auto ta2_yz_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 709);

    auto ta2_yz_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 710);

    auto ta2_yz_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 711);

    auto ta2_yz_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 712);

    auto ta2_yz_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 713);

    auto ta2_yz_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 714);

    auto ta2_yz_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 715);

    auto ta2_yz_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 716);

    auto ta2_yz_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 717);

    auto ta2_yz_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 718);

    auto ta2_yz_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 719);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yy_xxxx_1, ta1_y_yy_xxxy_1, ta1_y_yy_xxyy_1, ta1_y_yy_xxyz_1, ta1_y_yy_xyyy_1, ta1_y_yy_xyyz_1, ta1_y_yy_xyzz_1, ta1_y_yy_yyyy_1, ta1_y_yy_yyyz_1, ta1_y_yy_yyzz_1, ta1_y_yy_yzzz_1, ta1_z_yz_xxxz_1, ta1_z_yz_xxzz_1, ta1_z_yz_xzzz_1, ta1_z_yz_zzzz_1, ta2_yz_yy_xxxx_0, ta2_yz_yy_xxxx_1, ta2_yz_yy_xxxy_0, ta2_yz_yy_xxxy_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xxyy_0, ta2_yz_yy_xxyy_1, ta2_yz_yy_xxyz_0, ta2_yz_yy_xxyz_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyyy_0, ta2_yz_yy_xyyy_1, ta2_yz_yy_xyyz_0, ta2_yz_yy_xyyz_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_xyzz_0, ta2_yz_yy_xyzz_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyyy_0, ta2_yz_yy_yyyy_1, ta2_yz_yy_yyyz_0, ta2_yz_yy_yyyz_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yyzz_0, ta2_yz_yy_yyzz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_yzzz_0, ta2_yz_yy_yzzz_1, ta2_yz_yyz_xxxx_0, ta2_yz_yyz_xxxy_0, ta2_yz_yyz_xxxz_0, ta2_yz_yyz_xxyy_0, ta2_yz_yyz_xxyz_0, ta2_yz_yyz_xxzz_0, ta2_yz_yyz_xyyy_0, ta2_yz_yyz_xyyz_0, ta2_yz_yyz_xyzz_0, ta2_yz_yyz_xzzz_0, ta2_yz_yyz_yyyy_0, ta2_yz_yyz_yyyz_0, ta2_yz_yyz_yyzz_0, ta2_yz_yyz_yzzz_0, ta2_yz_yyz_zzzz_0, ta2_yz_yz_xxxz_0, ta2_yz_yz_xxxz_1, ta2_yz_yz_xxzz_0, ta2_yz_yz_xxzz_1, ta2_yz_yz_xzzz_0, ta2_yz_yz_xzzz_1, ta2_yz_yz_zzzz_0, ta2_yz_yz_zzzz_1, ta2_yz_z_xxxz_0, ta2_yz_z_xxxz_1, ta2_yz_z_xxzz_0, ta2_yz_z_xxzz_1, ta2_yz_z_xzzz_0, ta2_yz_z_xzzz_1, ta2_yz_z_zzzz_0, ta2_yz_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyz_xxxx_0[i] = ta1_y_yy_xxxx_1[i] + ta2_yz_yy_xxxx_0[i] * pa_z[i] - ta2_yz_yy_xxxx_1[i] * pc_z[i];

        ta2_yz_yyz_xxxy_0[i] = ta1_y_yy_xxxy_1[i] + ta2_yz_yy_xxxy_0[i] * pa_z[i] - ta2_yz_yy_xxxy_1[i] * pc_z[i];

        ta2_yz_yyz_xxxz_0[i] = ta2_yz_z_xxxz_0[i] * fe_0 - ta2_yz_z_xxxz_1[i] * fe_0 + ta1_z_yz_xxxz_1[i] + ta2_yz_yz_xxxz_0[i] * pa_y[i] - ta2_yz_yz_xxxz_1[i] * pc_y[i];

        ta2_yz_yyz_xxyy_0[i] = ta1_y_yy_xxyy_1[i] + ta2_yz_yy_xxyy_0[i] * pa_z[i] - ta2_yz_yy_xxyy_1[i] * pc_z[i];

        ta2_yz_yyz_xxyz_0[i] = ta2_yz_yy_xxy_0[i] * fe_0 - ta2_yz_yy_xxy_1[i] * fe_0 + ta1_y_yy_xxyz_1[i] + ta2_yz_yy_xxyz_0[i] * pa_z[i] - ta2_yz_yy_xxyz_1[i] * pc_z[i];

        ta2_yz_yyz_xxzz_0[i] = ta2_yz_z_xxzz_0[i] * fe_0 - ta2_yz_z_xxzz_1[i] * fe_0 + ta1_z_yz_xxzz_1[i] + ta2_yz_yz_xxzz_0[i] * pa_y[i] - ta2_yz_yz_xxzz_1[i] * pc_y[i];

        ta2_yz_yyz_xyyy_0[i] = ta1_y_yy_xyyy_1[i] + ta2_yz_yy_xyyy_0[i] * pa_z[i] - ta2_yz_yy_xyyy_1[i] * pc_z[i];

        ta2_yz_yyz_xyyz_0[i] = ta2_yz_yy_xyy_0[i] * fe_0 - ta2_yz_yy_xyy_1[i] * fe_0 + ta1_y_yy_xyyz_1[i] + ta2_yz_yy_xyyz_0[i] * pa_z[i] - ta2_yz_yy_xyyz_1[i] * pc_z[i];

        ta2_yz_yyz_xyzz_0[i] = 2.0 * ta2_yz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yy_xyz_1[i] * fe_0 + ta1_y_yy_xyzz_1[i] + ta2_yz_yy_xyzz_0[i] * pa_z[i] - ta2_yz_yy_xyzz_1[i] * pc_z[i];

        ta2_yz_yyz_xzzz_0[i] = ta2_yz_z_xzzz_0[i] * fe_0 - ta2_yz_z_xzzz_1[i] * fe_0 + ta1_z_yz_xzzz_1[i] + ta2_yz_yz_xzzz_0[i] * pa_y[i] - ta2_yz_yz_xzzz_1[i] * pc_y[i];

        ta2_yz_yyz_yyyy_0[i] = ta1_y_yy_yyyy_1[i] + ta2_yz_yy_yyyy_0[i] * pa_z[i] - ta2_yz_yy_yyyy_1[i] * pc_z[i];

        ta2_yz_yyz_yyyz_0[i] = ta2_yz_yy_yyy_0[i] * fe_0 - ta2_yz_yy_yyy_1[i] * fe_0 + ta1_y_yy_yyyz_1[i] + ta2_yz_yy_yyyz_0[i] * pa_z[i] - ta2_yz_yy_yyyz_1[i] * pc_z[i];

        ta2_yz_yyz_yyzz_0[i] = 2.0 * ta2_yz_yy_yyz_0[i] * fe_0 - 2.0 * ta2_yz_yy_yyz_1[i] * fe_0 + ta1_y_yy_yyzz_1[i] + ta2_yz_yy_yyzz_0[i] * pa_z[i] - ta2_yz_yy_yyzz_1[i] * pc_z[i];

        ta2_yz_yyz_yzzz_0[i] = 3.0 * ta2_yz_yy_yzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yzz_1[i] * fe_0 + ta1_y_yy_yzzz_1[i] + ta2_yz_yy_yzzz_0[i] * pa_z[i] - ta2_yz_yy_yzzz_1[i] * pc_z[i];

        ta2_yz_yyz_zzzz_0[i] = ta2_yz_z_zzzz_0[i] * fe_0 - ta2_yz_z_zzzz_1[i] * fe_0 + ta1_z_yz_zzzz_1[i] + ta2_yz_yz_zzzz_0[i] * pa_y[i] - ta2_yz_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 720-735 components of targeted buffer : FG

    auto ta2_yz_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 720);

    auto ta2_yz_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 721);

    auto ta2_yz_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 722);

    auto ta2_yz_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 723);

    auto ta2_yz_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 724);

    auto ta2_yz_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 725);

    auto ta2_yz_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 726);

    auto ta2_yz_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 727);

    auto ta2_yz_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 728);

    auto ta2_yz_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 729);

    auto ta2_yz_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 730);

    auto ta2_yz_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 731);

    auto ta2_yz_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 732);

    auto ta2_yz_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 733);

    auto ta2_yz_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 734);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_zz_xxxx_1, ta1_z_zz_xxxy_1, ta1_z_zz_xxxz_1, ta1_z_zz_xxyy_1, ta1_z_zz_xxyz_1, ta1_z_zz_xxzz_1, ta1_z_zz_xyyy_1, ta1_z_zz_xyyz_1, ta1_z_zz_xyzz_1, ta1_z_zz_xzzz_1, ta1_z_zz_yyyy_1, ta1_z_zz_yyyz_1, ta1_z_zz_yyzz_1, ta1_z_zz_yzzz_1, ta1_z_zz_zzzz_1, ta2_yz_yzz_xxxx_0, ta2_yz_yzz_xxxy_0, ta2_yz_yzz_xxxz_0, ta2_yz_yzz_xxyy_0, ta2_yz_yzz_xxyz_0, ta2_yz_yzz_xxzz_0, ta2_yz_yzz_xyyy_0, ta2_yz_yzz_xyyz_0, ta2_yz_yzz_xyzz_0, ta2_yz_yzz_xzzz_0, ta2_yz_yzz_yyyy_0, ta2_yz_yzz_yyyz_0, ta2_yz_yzz_yyzz_0, ta2_yz_yzz_yzzz_0, ta2_yz_yzz_zzzz_0, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxxx_0, ta2_yz_zz_xxxx_1, ta2_yz_zz_xxxy_0, ta2_yz_zz_xxxy_1, ta2_yz_zz_xxxz_0, ta2_yz_zz_xxxz_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxyy_0, ta2_yz_zz_xxyy_1, ta2_yz_zz_xxyz_0, ta2_yz_zz_xxyz_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xxzz_0, ta2_yz_zz_xxzz_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyyy_0, ta2_yz_zz_xyyy_1, ta2_yz_zz_xyyz_0, ta2_yz_zz_xyyz_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xyzz_0, ta2_yz_zz_xyzz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_xzzz_0, ta2_yz_zz_xzzz_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyyy_0, ta2_yz_zz_yyyy_1, ta2_yz_zz_yyyz_0, ta2_yz_zz_yyyz_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yyzz_0, ta2_yz_zz_yyzz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_yzzz_0, ta2_yz_zz_yzzz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, ta2_yz_zz_zzzz_0, ta2_yz_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzz_xxxx_0[i] = ta1_z_zz_xxxx_1[i] + ta2_yz_zz_xxxx_0[i] * pa_y[i] - ta2_yz_zz_xxxx_1[i] * pc_y[i];

        ta2_yz_yzz_xxxy_0[i] = ta2_yz_zz_xxx_0[i] * fe_0 - ta2_yz_zz_xxx_1[i] * fe_0 + ta1_z_zz_xxxy_1[i] + ta2_yz_zz_xxxy_0[i] * pa_y[i] - ta2_yz_zz_xxxy_1[i] * pc_y[i];

        ta2_yz_yzz_xxxz_0[i] = ta1_z_zz_xxxz_1[i] + ta2_yz_zz_xxxz_0[i] * pa_y[i] - ta2_yz_zz_xxxz_1[i] * pc_y[i];

        ta2_yz_yzz_xxyy_0[i] = 2.0 * ta2_yz_zz_xxy_0[i] * fe_0 - 2.0 * ta2_yz_zz_xxy_1[i] * fe_0 + ta1_z_zz_xxyy_1[i] + ta2_yz_zz_xxyy_0[i] * pa_y[i] - ta2_yz_zz_xxyy_1[i] * pc_y[i];

        ta2_yz_yzz_xxyz_0[i] = ta2_yz_zz_xxz_0[i] * fe_0 - ta2_yz_zz_xxz_1[i] * fe_0 + ta1_z_zz_xxyz_1[i] + ta2_yz_zz_xxyz_0[i] * pa_y[i] - ta2_yz_zz_xxyz_1[i] * pc_y[i];

        ta2_yz_yzz_xxzz_0[i] = ta1_z_zz_xxzz_1[i] + ta2_yz_zz_xxzz_0[i] * pa_y[i] - ta2_yz_zz_xxzz_1[i] * pc_y[i];

        ta2_yz_yzz_xyyy_0[i] = 3.0 * ta2_yz_zz_xyy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xyy_1[i] * fe_0 + ta1_z_zz_xyyy_1[i] + ta2_yz_zz_xyyy_0[i] * pa_y[i] - ta2_yz_zz_xyyy_1[i] * pc_y[i];

        ta2_yz_yzz_xyyz_0[i] = 2.0 * ta2_yz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xyz_1[i] * fe_0 + ta1_z_zz_xyyz_1[i] + ta2_yz_zz_xyyz_0[i] * pa_y[i] - ta2_yz_zz_xyyz_1[i] * pc_y[i];

        ta2_yz_yzz_xyzz_0[i] = ta2_yz_zz_xzz_0[i] * fe_0 - ta2_yz_zz_xzz_1[i] * fe_0 + ta1_z_zz_xyzz_1[i] + ta2_yz_zz_xyzz_0[i] * pa_y[i] - ta2_yz_zz_xyzz_1[i] * pc_y[i];

        ta2_yz_yzz_xzzz_0[i] = ta1_z_zz_xzzz_1[i] + ta2_yz_zz_xzzz_0[i] * pa_y[i] - ta2_yz_zz_xzzz_1[i] * pc_y[i];

        ta2_yz_yzz_yyyy_0[i] = 4.0 * ta2_yz_zz_yyy_0[i] * fe_0 - 4.0 * ta2_yz_zz_yyy_1[i] * fe_0 + ta1_z_zz_yyyy_1[i] + ta2_yz_zz_yyyy_0[i] * pa_y[i] - ta2_yz_zz_yyyy_1[i] * pc_y[i];

        ta2_yz_yzz_yyyz_0[i] = 3.0 * ta2_yz_zz_yyz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yyz_1[i] * fe_0 + ta1_z_zz_yyyz_1[i] + ta2_yz_zz_yyyz_0[i] * pa_y[i] - ta2_yz_zz_yyyz_1[i] * pc_y[i];

        ta2_yz_yzz_yyzz_0[i] = 2.0 * ta2_yz_zz_yzz_0[i] * fe_0 - 2.0 * ta2_yz_zz_yzz_1[i] * fe_0 + ta1_z_zz_yyzz_1[i] + ta2_yz_zz_yyzz_0[i] * pa_y[i] - ta2_yz_zz_yyzz_1[i] * pc_y[i];

        ta2_yz_yzz_yzzz_0[i] = ta2_yz_zz_zzz_0[i] * fe_0 - ta2_yz_zz_zzz_1[i] * fe_0 + ta1_z_zz_yzzz_1[i] + ta2_yz_zz_yzzz_0[i] * pa_y[i] - ta2_yz_zz_yzzz_1[i] * pc_y[i];

        ta2_yz_yzz_zzzz_0[i] = ta1_z_zz_zzzz_1[i] + ta2_yz_zz_zzzz_0[i] * pa_y[i] - ta2_yz_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 735-750 components of targeted buffer : FG

    auto ta2_yz_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 735);

    auto ta2_yz_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 736);

    auto ta2_yz_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 737);

    auto ta2_yz_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 738);

    auto ta2_yz_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 739);

    auto ta2_yz_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 740);

    auto ta2_yz_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 741);

    auto ta2_yz_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 742);

    auto ta2_yz_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 743);

    auto ta2_yz_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 744);

    auto ta2_yz_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 745);

    auto ta2_yz_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 746);

    auto ta2_yz_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 747);

    auto ta2_yz_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 748);

    auto ta2_yz_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 749);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zz_xxxx_1, ta1_y_zz_xxxy_1, ta1_y_zz_xxxz_1, ta1_y_zz_xxyy_1, ta1_y_zz_xxyz_1, ta1_y_zz_xxzz_1, ta1_y_zz_xyyy_1, ta1_y_zz_xyyz_1, ta1_y_zz_xyzz_1, ta1_y_zz_xzzz_1, ta1_y_zz_yyyy_1, ta1_y_zz_yyyz_1, ta1_y_zz_yyzz_1, ta1_y_zz_yzzz_1, ta1_y_zz_zzzz_1, ta2_yz_z_xxxx_0, ta2_yz_z_xxxx_1, ta2_yz_z_xxxy_0, ta2_yz_z_xxxy_1, ta2_yz_z_xxxz_0, ta2_yz_z_xxxz_1, ta2_yz_z_xxyy_0, ta2_yz_z_xxyy_1, ta2_yz_z_xxyz_0, ta2_yz_z_xxyz_1, ta2_yz_z_xxzz_0, ta2_yz_z_xxzz_1, ta2_yz_z_xyyy_0, ta2_yz_z_xyyy_1, ta2_yz_z_xyyz_0, ta2_yz_z_xyyz_1, ta2_yz_z_xyzz_0, ta2_yz_z_xyzz_1, ta2_yz_z_xzzz_0, ta2_yz_z_xzzz_1, ta2_yz_z_yyyy_0, ta2_yz_z_yyyy_1, ta2_yz_z_yyyz_0, ta2_yz_z_yyyz_1, ta2_yz_z_yyzz_0, ta2_yz_z_yyzz_1, ta2_yz_z_yzzz_0, ta2_yz_z_yzzz_1, ta2_yz_z_zzzz_0, ta2_yz_z_zzzz_1, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxxx_0, ta2_yz_zz_xxxx_1, ta2_yz_zz_xxxy_0, ta2_yz_zz_xxxy_1, ta2_yz_zz_xxxz_0, ta2_yz_zz_xxxz_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxyy_0, ta2_yz_zz_xxyy_1, ta2_yz_zz_xxyz_0, ta2_yz_zz_xxyz_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xxzz_0, ta2_yz_zz_xxzz_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyyy_0, ta2_yz_zz_xyyy_1, ta2_yz_zz_xyyz_0, ta2_yz_zz_xyyz_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xyzz_0, ta2_yz_zz_xyzz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_xzzz_0, ta2_yz_zz_xzzz_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyyy_0, ta2_yz_zz_yyyy_1, ta2_yz_zz_yyyz_0, ta2_yz_zz_yyyz_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yyzz_0, ta2_yz_zz_yyzz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_yzzz_0, ta2_yz_zz_yzzz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, ta2_yz_zz_zzzz_0, ta2_yz_zz_zzzz_1, ta2_yz_zzz_xxxx_0, ta2_yz_zzz_xxxy_0, ta2_yz_zzz_xxxz_0, ta2_yz_zzz_xxyy_0, ta2_yz_zzz_xxyz_0, ta2_yz_zzz_xxzz_0, ta2_yz_zzz_xyyy_0, ta2_yz_zzz_xyyz_0, ta2_yz_zzz_xyzz_0, ta2_yz_zzz_xzzz_0, ta2_yz_zzz_yyyy_0, ta2_yz_zzz_yyyz_0, ta2_yz_zzz_yyzz_0, ta2_yz_zzz_yzzz_0, ta2_yz_zzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzz_xxxx_0[i] = 2.0 * ta2_yz_z_xxxx_0[i] * fe_0 - 2.0 * ta2_yz_z_xxxx_1[i] * fe_0 + ta1_y_zz_xxxx_1[i] + ta2_yz_zz_xxxx_0[i] * pa_z[i] - ta2_yz_zz_xxxx_1[i] * pc_z[i];

        ta2_yz_zzz_xxxy_0[i] = 2.0 * ta2_yz_z_xxxy_0[i] * fe_0 - 2.0 * ta2_yz_z_xxxy_1[i] * fe_0 + ta1_y_zz_xxxy_1[i] + ta2_yz_zz_xxxy_0[i] * pa_z[i] - ta2_yz_zz_xxxy_1[i] * pc_z[i];

        ta2_yz_zzz_xxxz_0[i] = 2.0 * ta2_yz_z_xxxz_0[i] * fe_0 - 2.0 * ta2_yz_z_xxxz_1[i] * fe_0 + ta2_yz_zz_xxx_0[i] * fe_0 - ta2_yz_zz_xxx_1[i] * fe_0 + ta1_y_zz_xxxz_1[i] + ta2_yz_zz_xxxz_0[i] * pa_z[i] - ta2_yz_zz_xxxz_1[i] * pc_z[i];

        ta2_yz_zzz_xxyy_0[i] = 2.0 * ta2_yz_z_xxyy_0[i] * fe_0 - 2.0 * ta2_yz_z_xxyy_1[i] * fe_0 + ta1_y_zz_xxyy_1[i] + ta2_yz_zz_xxyy_0[i] * pa_z[i] - ta2_yz_zz_xxyy_1[i] * pc_z[i];

        ta2_yz_zzz_xxyz_0[i] = 2.0 * ta2_yz_z_xxyz_0[i] * fe_0 - 2.0 * ta2_yz_z_xxyz_1[i] * fe_0 + ta2_yz_zz_xxy_0[i] * fe_0 - ta2_yz_zz_xxy_1[i] * fe_0 + ta1_y_zz_xxyz_1[i] + ta2_yz_zz_xxyz_0[i] * pa_z[i] - ta2_yz_zz_xxyz_1[i] * pc_z[i];

        ta2_yz_zzz_xxzz_0[i] = 2.0 * ta2_yz_z_xxzz_0[i] * fe_0 - 2.0 * ta2_yz_z_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_zz_xxz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xxz_1[i] * fe_0 + ta1_y_zz_xxzz_1[i] + ta2_yz_zz_xxzz_0[i] * pa_z[i] - ta2_yz_zz_xxzz_1[i] * pc_z[i];

        ta2_yz_zzz_xyyy_0[i] = 2.0 * ta2_yz_z_xyyy_0[i] * fe_0 - 2.0 * ta2_yz_z_xyyy_1[i] * fe_0 + ta1_y_zz_xyyy_1[i] + ta2_yz_zz_xyyy_0[i] * pa_z[i] - ta2_yz_zz_xyyy_1[i] * pc_z[i];

        ta2_yz_zzz_xyyz_0[i] = 2.0 * ta2_yz_z_xyyz_0[i] * fe_0 - 2.0 * ta2_yz_z_xyyz_1[i] * fe_0 + ta2_yz_zz_xyy_0[i] * fe_0 - ta2_yz_zz_xyy_1[i] * fe_0 + ta1_y_zz_xyyz_1[i] + ta2_yz_zz_xyyz_0[i] * pa_z[i] - ta2_yz_zz_xyyz_1[i] * pc_z[i];

        ta2_yz_zzz_xyzz_0[i] = 2.0 * ta2_yz_z_xyzz_0[i] * fe_0 - 2.0 * ta2_yz_z_xyzz_1[i] * fe_0 + 2.0 * ta2_yz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xyz_1[i] * fe_0 + ta1_y_zz_xyzz_1[i] + ta2_yz_zz_xyzz_0[i] * pa_z[i] - ta2_yz_zz_xyzz_1[i] * pc_z[i];

        ta2_yz_zzz_xzzz_0[i] = 2.0 * ta2_yz_z_xzzz_0[i] * fe_0 - 2.0 * ta2_yz_z_xzzz_1[i] * fe_0 + 3.0 * ta2_yz_zz_xzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xzz_1[i] * fe_0 + ta1_y_zz_xzzz_1[i] + ta2_yz_zz_xzzz_0[i] * pa_z[i] - ta2_yz_zz_xzzz_1[i] * pc_z[i];

        ta2_yz_zzz_yyyy_0[i] = 2.0 * ta2_yz_z_yyyy_0[i] * fe_0 - 2.0 * ta2_yz_z_yyyy_1[i] * fe_0 + ta1_y_zz_yyyy_1[i] + ta2_yz_zz_yyyy_0[i] * pa_z[i] - ta2_yz_zz_yyyy_1[i] * pc_z[i];

        ta2_yz_zzz_yyyz_0[i] = 2.0 * ta2_yz_z_yyyz_0[i] * fe_0 - 2.0 * ta2_yz_z_yyyz_1[i] * fe_0 + ta2_yz_zz_yyy_0[i] * fe_0 - ta2_yz_zz_yyy_1[i] * fe_0 + ta1_y_zz_yyyz_1[i] + ta2_yz_zz_yyyz_0[i] * pa_z[i] - ta2_yz_zz_yyyz_1[i] * pc_z[i];

        ta2_yz_zzz_yyzz_0[i] = 2.0 * ta2_yz_z_yyzz_0[i] * fe_0 - 2.0 * ta2_yz_z_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_zz_yyz_0[i] * fe_0 - 2.0 * ta2_yz_zz_yyz_1[i] * fe_0 + ta1_y_zz_yyzz_1[i] + ta2_yz_zz_yyzz_0[i] * pa_z[i] - ta2_yz_zz_yyzz_1[i] * pc_z[i];

        ta2_yz_zzz_yzzz_0[i] = 2.0 * ta2_yz_z_yzzz_0[i] * fe_0 - 2.0 * ta2_yz_z_yzzz_1[i] * fe_0 + 3.0 * ta2_yz_zz_yzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yzz_1[i] * fe_0 + ta1_y_zz_yzzz_1[i] + ta2_yz_zz_yzzz_0[i] * pa_z[i] - ta2_yz_zz_yzzz_1[i] * pc_z[i];

        ta2_yz_zzz_zzzz_0[i] = 2.0 * ta2_yz_z_zzzz_0[i] * fe_0 - 2.0 * ta2_yz_z_zzzz_1[i] * fe_0 + 4.0 * ta2_yz_zz_zzz_0[i] * fe_0 - 4.0 * ta2_yz_zz_zzz_1[i] * fe_0 + ta1_y_zz_zzzz_1[i] + ta2_yz_zz_zzzz_0[i] * pa_z[i] - ta2_yz_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 750-765 components of targeted buffer : FG

    auto ta2_zz_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 750);

    auto ta2_zz_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 751);

    auto ta2_zz_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 752);

    auto ta2_zz_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 753);

    auto ta2_zz_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 754);

    auto ta2_zz_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 755);

    auto ta2_zz_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 756);

    auto ta2_zz_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 757);

    auto ta2_zz_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 758);

    auto ta2_zz_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 759);

    auto ta2_zz_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 760);

    auto ta2_zz_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 761);

    auto ta2_zz_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 762);

    auto ta2_zz_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 763);

    auto ta2_zz_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 764);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_x_xxxx_0, ta2_zz_x_xxxx_1, ta2_zz_x_xxxy_0, ta2_zz_x_xxxy_1, ta2_zz_x_xxxz_0, ta2_zz_x_xxxz_1, ta2_zz_x_xxyy_0, ta2_zz_x_xxyy_1, ta2_zz_x_xxyz_0, ta2_zz_x_xxyz_1, ta2_zz_x_xxzz_0, ta2_zz_x_xxzz_1, ta2_zz_x_xyyy_0, ta2_zz_x_xyyy_1, ta2_zz_x_xyyz_0, ta2_zz_x_xyyz_1, ta2_zz_x_xyzz_0, ta2_zz_x_xyzz_1, ta2_zz_x_xzzz_0, ta2_zz_x_xzzz_1, ta2_zz_x_yyyy_0, ta2_zz_x_yyyy_1, ta2_zz_x_yyyz_0, ta2_zz_x_yyyz_1, ta2_zz_x_yyzz_0, ta2_zz_x_yyzz_1, ta2_zz_x_yzzz_0, ta2_zz_x_yzzz_1, ta2_zz_x_zzzz_0, ta2_zz_x_zzzz_1, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxxx_0, ta2_zz_xx_xxxx_1, ta2_zz_xx_xxxy_0, ta2_zz_xx_xxxy_1, ta2_zz_xx_xxxz_0, ta2_zz_xx_xxxz_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxyy_0, ta2_zz_xx_xxyy_1, ta2_zz_xx_xxyz_0, ta2_zz_xx_xxyz_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xxzz_0, ta2_zz_xx_xxzz_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyyy_0, ta2_zz_xx_xyyy_1, ta2_zz_xx_xyyz_0, ta2_zz_xx_xyyz_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xyzz_0, ta2_zz_xx_xyzz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_xzzz_0, ta2_zz_xx_xzzz_1, ta2_zz_xx_yyy_0, ta2_zz_xx_yyy_1, ta2_zz_xx_yyyy_0, ta2_zz_xx_yyyy_1, ta2_zz_xx_yyyz_0, ta2_zz_xx_yyyz_1, ta2_zz_xx_yyz_0, ta2_zz_xx_yyz_1, ta2_zz_xx_yyzz_0, ta2_zz_xx_yyzz_1, ta2_zz_xx_yzz_0, ta2_zz_xx_yzz_1, ta2_zz_xx_yzzz_0, ta2_zz_xx_yzzz_1, ta2_zz_xx_zzz_0, ta2_zz_xx_zzz_1, ta2_zz_xx_zzzz_0, ta2_zz_xx_zzzz_1, ta2_zz_xxx_xxxx_0, ta2_zz_xxx_xxxy_0, ta2_zz_xxx_xxxz_0, ta2_zz_xxx_xxyy_0, ta2_zz_xxx_xxyz_0, ta2_zz_xxx_xxzz_0, ta2_zz_xxx_xyyy_0, ta2_zz_xxx_xyyz_0, ta2_zz_xxx_xyzz_0, ta2_zz_xxx_xzzz_0, ta2_zz_xxx_yyyy_0, ta2_zz_xxx_yyyz_0, ta2_zz_xxx_yyzz_0, ta2_zz_xxx_yzzz_0, ta2_zz_xxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxx_xxxx_0[i] = 2.0 * ta2_zz_x_xxxx_0[i] * fe_0 - 2.0 * ta2_zz_x_xxxx_1[i] * fe_0 + 4.0 * ta2_zz_xx_xxx_0[i] * fe_0 - 4.0 * ta2_zz_xx_xxx_1[i] * fe_0 + ta2_zz_xx_xxxx_0[i] * pa_x[i] - ta2_zz_xx_xxxx_1[i] * pc_x[i];

        ta2_zz_xxx_xxxy_0[i] = 2.0 * ta2_zz_x_xxxy_0[i] * fe_0 - 2.0 * ta2_zz_x_xxxy_1[i] * fe_0 + 3.0 * ta2_zz_xx_xxy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxy_1[i] * fe_0 + ta2_zz_xx_xxxy_0[i] * pa_x[i] - ta2_zz_xx_xxxy_1[i] * pc_x[i];

        ta2_zz_xxx_xxxz_0[i] = 2.0 * ta2_zz_x_xxxz_0[i] * fe_0 - 2.0 * ta2_zz_x_xxxz_1[i] * fe_0 + 3.0 * ta2_zz_xx_xxz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxz_1[i] * fe_0 + ta2_zz_xx_xxxz_0[i] * pa_x[i] - ta2_zz_xx_xxxz_1[i] * pc_x[i];

        ta2_zz_xxx_xxyy_0[i] = 2.0 * ta2_zz_x_xxyy_0[i] * fe_0 - 2.0 * ta2_zz_x_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_xx_xyy_0[i] * fe_0 - 2.0 * ta2_zz_xx_xyy_1[i] * fe_0 + ta2_zz_xx_xxyy_0[i] * pa_x[i] - ta2_zz_xx_xxyy_1[i] * pc_x[i];

        ta2_zz_xxx_xxyz_0[i] = 2.0 * ta2_zz_x_xxyz_0[i] * fe_0 - 2.0 * ta2_zz_x_xxyz_1[i] * fe_0 + 2.0 * ta2_zz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xyz_1[i] * fe_0 + ta2_zz_xx_xxyz_0[i] * pa_x[i] - ta2_zz_xx_xxyz_1[i] * pc_x[i];

        ta2_zz_xxx_xxzz_0[i] = 2.0 * ta2_zz_x_xxzz_0[i] * fe_0 - 2.0 * ta2_zz_x_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_xx_xzz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xzz_1[i] * fe_0 + ta2_zz_xx_xxzz_0[i] * pa_x[i] - ta2_zz_xx_xxzz_1[i] * pc_x[i];

        ta2_zz_xxx_xyyy_0[i] = 2.0 * ta2_zz_x_xyyy_0[i] * fe_0 - 2.0 * ta2_zz_x_xyyy_1[i] * fe_0 + ta2_zz_xx_yyy_0[i] * fe_0 - ta2_zz_xx_yyy_1[i] * fe_0 + ta2_zz_xx_xyyy_0[i] * pa_x[i] - ta2_zz_xx_xyyy_1[i] * pc_x[i];

        ta2_zz_xxx_xyyz_0[i] = 2.0 * ta2_zz_x_xyyz_0[i] * fe_0 - 2.0 * ta2_zz_x_xyyz_1[i] * fe_0 + ta2_zz_xx_yyz_0[i] * fe_0 - ta2_zz_xx_yyz_1[i] * fe_0 + ta2_zz_xx_xyyz_0[i] * pa_x[i] - ta2_zz_xx_xyyz_1[i] * pc_x[i];

        ta2_zz_xxx_xyzz_0[i] = 2.0 * ta2_zz_x_xyzz_0[i] * fe_0 - 2.0 * ta2_zz_x_xyzz_1[i] * fe_0 + ta2_zz_xx_yzz_0[i] * fe_0 - ta2_zz_xx_yzz_1[i] * fe_0 + ta2_zz_xx_xyzz_0[i] * pa_x[i] - ta2_zz_xx_xyzz_1[i] * pc_x[i];

        ta2_zz_xxx_xzzz_0[i] = 2.0 * ta2_zz_x_xzzz_0[i] * fe_0 - 2.0 * ta2_zz_x_xzzz_1[i] * fe_0 + ta2_zz_xx_zzz_0[i] * fe_0 - ta2_zz_xx_zzz_1[i] * fe_0 + ta2_zz_xx_xzzz_0[i] * pa_x[i] - ta2_zz_xx_xzzz_1[i] * pc_x[i];

        ta2_zz_xxx_yyyy_0[i] = 2.0 * ta2_zz_x_yyyy_0[i] * fe_0 - 2.0 * ta2_zz_x_yyyy_1[i] * fe_0 + ta2_zz_xx_yyyy_0[i] * pa_x[i] - ta2_zz_xx_yyyy_1[i] * pc_x[i];

        ta2_zz_xxx_yyyz_0[i] = 2.0 * ta2_zz_x_yyyz_0[i] * fe_0 - 2.0 * ta2_zz_x_yyyz_1[i] * fe_0 + ta2_zz_xx_yyyz_0[i] * pa_x[i] - ta2_zz_xx_yyyz_1[i] * pc_x[i];

        ta2_zz_xxx_yyzz_0[i] = 2.0 * ta2_zz_x_yyzz_0[i] * fe_0 - 2.0 * ta2_zz_x_yyzz_1[i] * fe_0 + ta2_zz_xx_yyzz_0[i] * pa_x[i] - ta2_zz_xx_yyzz_1[i] * pc_x[i];

        ta2_zz_xxx_yzzz_0[i] = 2.0 * ta2_zz_x_yzzz_0[i] * fe_0 - 2.0 * ta2_zz_x_yzzz_1[i] * fe_0 + ta2_zz_xx_yzzz_0[i] * pa_x[i] - ta2_zz_xx_yzzz_1[i] * pc_x[i];

        ta2_zz_xxx_zzzz_0[i] = 2.0 * ta2_zz_x_zzzz_0[i] * fe_0 - 2.0 * ta2_zz_x_zzzz_1[i] * fe_0 + ta2_zz_xx_zzzz_0[i] * pa_x[i] - ta2_zz_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 765-780 components of targeted buffer : FG

    auto ta2_zz_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 765);

    auto ta2_zz_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 766);

    auto ta2_zz_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 767);

    auto ta2_zz_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 768);

    auto ta2_zz_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 769);

    auto ta2_zz_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 770);

    auto ta2_zz_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 771);

    auto ta2_zz_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 772);

    auto ta2_zz_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 773);

    auto ta2_zz_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 774);

    auto ta2_zz_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 775);

    auto ta2_zz_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 776);

    auto ta2_zz_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 777);

    auto ta2_zz_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 778);

    auto ta2_zz_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 779);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxxx_0, ta2_zz_xx_xxxx_1, ta2_zz_xx_xxxy_0, ta2_zz_xx_xxxy_1, ta2_zz_xx_xxxz_0, ta2_zz_xx_xxxz_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxyy_0, ta2_zz_xx_xxyy_1, ta2_zz_xx_xxyz_0, ta2_zz_xx_xxyz_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xxzz_0, ta2_zz_xx_xxzz_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyyy_0, ta2_zz_xx_xyyy_1, ta2_zz_xx_xyyz_0, ta2_zz_xx_xyyz_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xyzz_0, ta2_zz_xx_xyzz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_xzzz_0, ta2_zz_xx_xzzz_1, ta2_zz_xx_zzzz_0, ta2_zz_xx_zzzz_1, ta2_zz_xxy_xxxx_0, ta2_zz_xxy_xxxy_0, ta2_zz_xxy_xxxz_0, ta2_zz_xxy_xxyy_0, ta2_zz_xxy_xxyz_0, ta2_zz_xxy_xxzz_0, ta2_zz_xxy_xyyy_0, ta2_zz_xxy_xyyz_0, ta2_zz_xxy_xyzz_0, ta2_zz_xxy_xzzz_0, ta2_zz_xxy_yyyy_0, ta2_zz_xxy_yyyz_0, ta2_zz_xxy_yyzz_0, ta2_zz_xxy_yzzz_0, ta2_zz_xxy_zzzz_0, ta2_zz_xy_yyyy_0, ta2_zz_xy_yyyy_1, ta2_zz_xy_yyyz_0, ta2_zz_xy_yyyz_1, ta2_zz_xy_yyzz_0, ta2_zz_xy_yyzz_1, ta2_zz_xy_yzzz_0, ta2_zz_xy_yzzz_1, ta2_zz_y_yyyy_0, ta2_zz_y_yyyy_1, ta2_zz_y_yyyz_0, ta2_zz_y_yyyz_1, ta2_zz_y_yyzz_0, ta2_zz_y_yyzz_1, ta2_zz_y_yzzz_0, ta2_zz_y_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxy_xxxx_0[i] = ta2_zz_xx_xxxx_0[i] * pa_y[i] - ta2_zz_xx_xxxx_1[i] * pc_y[i];

        ta2_zz_xxy_xxxy_0[i] = ta2_zz_xx_xxx_0[i] * fe_0 - ta2_zz_xx_xxx_1[i] * fe_0 + ta2_zz_xx_xxxy_0[i] * pa_y[i] - ta2_zz_xx_xxxy_1[i] * pc_y[i];

        ta2_zz_xxy_xxxz_0[i] = ta2_zz_xx_xxxz_0[i] * pa_y[i] - ta2_zz_xx_xxxz_1[i] * pc_y[i];

        ta2_zz_xxy_xxyy_0[i] = 2.0 * ta2_zz_xx_xxy_0[i] * fe_0 - 2.0 * ta2_zz_xx_xxy_1[i] * fe_0 + ta2_zz_xx_xxyy_0[i] * pa_y[i] - ta2_zz_xx_xxyy_1[i] * pc_y[i];

        ta2_zz_xxy_xxyz_0[i] = ta2_zz_xx_xxz_0[i] * fe_0 - ta2_zz_xx_xxz_1[i] * fe_0 + ta2_zz_xx_xxyz_0[i] * pa_y[i] - ta2_zz_xx_xxyz_1[i] * pc_y[i];

        ta2_zz_xxy_xxzz_0[i] = ta2_zz_xx_xxzz_0[i] * pa_y[i] - ta2_zz_xx_xxzz_1[i] * pc_y[i];

        ta2_zz_xxy_xyyy_0[i] = 3.0 * ta2_zz_xx_xyy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xyy_1[i] * fe_0 + ta2_zz_xx_xyyy_0[i] * pa_y[i] - ta2_zz_xx_xyyy_1[i] * pc_y[i];

        ta2_zz_xxy_xyyz_0[i] = 2.0 * ta2_zz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xyz_1[i] * fe_0 + ta2_zz_xx_xyyz_0[i] * pa_y[i] - ta2_zz_xx_xyyz_1[i] * pc_y[i];

        ta2_zz_xxy_xyzz_0[i] = ta2_zz_xx_xzz_0[i] * fe_0 - ta2_zz_xx_xzz_1[i] * fe_0 + ta2_zz_xx_xyzz_0[i] * pa_y[i] - ta2_zz_xx_xyzz_1[i] * pc_y[i];

        ta2_zz_xxy_xzzz_0[i] = ta2_zz_xx_xzzz_0[i] * pa_y[i] - ta2_zz_xx_xzzz_1[i] * pc_y[i];

        ta2_zz_xxy_yyyy_0[i] = ta2_zz_y_yyyy_0[i] * fe_0 - ta2_zz_y_yyyy_1[i] * fe_0 + ta2_zz_xy_yyyy_0[i] * pa_x[i] - ta2_zz_xy_yyyy_1[i] * pc_x[i];

        ta2_zz_xxy_yyyz_0[i] = ta2_zz_y_yyyz_0[i] * fe_0 - ta2_zz_y_yyyz_1[i] * fe_0 + ta2_zz_xy_yyyz_0[i] * pa_x[i] - ta2_zz_xy_yyyz_1[i] * pc_x[i];

        ta2_zz_xxy_yyzz_0[i] = ta2_zz_y_yyzz_0[i] * fe_0 - ta2_zz_y_yyzz_1[i] * fe_0 + ta2_zz_xy_yyzz_0[i] * pa_x[i] - ta2_zz_xy_yyzz_1[i] * pc_x[i];

        ta2_zz_xxy_yzzz_0[i] = ta2_zz_y_yzzz_0[i] * fe_0 - ta2_zz_y_yzzz_1[i] * fe_0 + ta2_zz_xy_yzzz_0[i] * pa_x[i] - ta2_zz_xy_yzzz_1[i] * pc_x[i];

        ta2_zz_xxy_zzzz_0[i] = ta2_zz_xx_zzzz_0[i] * pa_y[i] - ta2_zz_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 780-795 components of targeted buffer : FG

    auto ta2_zz_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 780);

    auto ta2_zz_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 781);

    auto ta2_zz_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 782);

    auto ta2_zz_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 783);

    auto ta2_zz_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 784);

    auto ta2_zz_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 785);

    auto ta2_zz_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 786);

    auto ta2_zz_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 787);

    auto ta2_zz_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 788);

    auto ta2_zz_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 789);

    auto ta2_zz_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 790);

    auto ta2_zz_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 791);

    auto ta2_zz_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 792);

    auto ta2_zz_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 793);

    auto ta2_zz_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 794);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xx_xxxx_1, ta1_z_xx_xxxy_1, ta1_z_xx_xxxz_1, ta1_z_xx_xxyy_1, ta1_z_xx_xxyz_1, ta1_z_xx_xxzz_1, ta1_z_xx_xyyy_1, ta1_z_xx_xyyz_1, ta1_z_xx_xyzz_1, ta1_z_xx_xzzz_1, ta1_z_xx_yyyy_1, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxxx_0, ta2_zz_xx_xxxx_1, ta2_zz_xx_xxxy_0, ta2_zz_xx_xxxy_1, ta2_zz_xx_xxxz_0, ta2_zz_xx_xxxz_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxyy_0, ta2_zz_xx_xxyy_1, ta2_zz_xx_xxyz_0, ta2_zz_xx_xxyz_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xxzz_0, ta2_zz_xx_xxzz_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyyy_0, ta2_zz_xx_xyyy_1, ta2_zz_xx_xyyz_0, ta2_zz_xx_xyyz_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xyzz_0, ta2_zz_xx_xyzz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_xzzz_0, ta2_zz_xx_xzzz_1, ta2_zz_xx_yyyy_0, ta2_zz_xx_yyyy_1, ta2_zz_xxz_xxxx_0, ta2_zz_xxz_xxxy_0, ta2_zz_xxz_xxxz_0, ta2_zz_xxz_xxyy_0, ta2_zz_xxz_xxyz_0, ta2_zz_xxz_xxzz_0, ta2_zz_xxz_xyyy_0, ta2_zz_xxz_xyyz_0, ta2_zz_xxz_xyzz_0, ta2_zz_xxz_xzzz_0, ta2_zz_xxz_yyyy_0, ta2_zz_xxz_yyyz_0, ta2_zz_xxz_yyzz_0, ta2_zz_xxz_yzzz_0, ta2_zz_xxz_zzzz_0, ta2_zz_xz_yyyz_0, ta2_zz_xz_yyyz_1, ta2_zz_xz_yyzz_0, ta2_zz_xz_yyzz_1, ta2_zz_xz_yzzz_0, ta2_zz_xz_yzzz_1, ta2_zz_xz_zzzz_0, ta2_zz_xz_zzzz_1, ta2_zz_z_yyyz_0, ta2_zz_z_yyyz_1, ta2_zz_z_yyzz_0, ta2_zz_z_yyzz_1, ta2_zz_z_yzzz_0, ta2_zz_z_yzzz_1, ta2_zz_z_zzzz_0, ta2_zz_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxz_xxxx_0[i] = 2.0 * ta1_z_xx_xxxx_1[i] + ta2_zz_xx_xxxx_0[i] * pa_z[i] - ta2_zz_xx_xxxx_1[i] * pc_z[i];

        ta2_zz_xxz_xxxy_0[i] = 2.0 * ta1_z_xx_xxxy_1[i] + ta2_zz_xx_xxxy_0[i] * pa_z[i] - ta2_zz_xx_xxxy_1[i] * pc_z[i];

        ta2_zz_xxz_xxxz_0[i] = ta2_zz_xx_xxx_0[i] * fe_0 - ta2_zz_xx_xxx_1[i] * fe_0 + 2.0 * ta1_z_xx_xxxz_1[i] + ta2_zz_xx_xxxz_0[i] * pa_z[i] - ta2_zz_xx_xxxz_1[i] * pc_z[i];

        ta2_zz_xxz_xxyy_0[i] = 2.0 * ta1_z_xx_xxyy_1[i] + ta2_zz_xx_xxyy_0[i] * pa_z[i] - ta2_zz_xx_xxyy_1[i] * pc_z[i];

        ta2_zz_xxz_xxyz_0[i] = ta2_zz_xx_xxy_0[i] * fe_0 - ta2_zz_xx_xxy_1[i] * fe_0 + 2.0 * ta1_z_xx_xxyz_1[i] + ta2_zz_xx_xxyz_0[i] * pa_z[i] - ta2_zz_xx_xxyz_1[i] * pc_z[i];

        ta2_zz_xxz_xxzz_0[i] = 2.0 * ta2_zz_xx_xxz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xxz_1[i] * fe_0 + 2.0 * ta1_z_xx_xxzz_1[i] + ta2_zz_xx_xxzz_0[i] * pa_z[i] - ta2_zz_xx_xxzz_1[i] * pc_z[i];

        ta2_zz_xxz_xyyy_0[i] = 2.0 * ta1_z_xx_xyyy_1[i] + ta2_zz_xx_xyyy_0[i] * pa_z[i] - ta2_zz_xx_xyyy_1[i] * pc_z[i];

        ta2_zz_xxz_xyyz_0[i] = ta2_zz_xx_xyy_0[i] * fe_0 - ta2_zz_xx_xyy_1[i] * fe_0 + 2.0 * ta1_z_xx_xyyz_1[i] + ta2_zz_xx_xyyz_0[i] * pa_z[i] - ta2_zz_xx_xyyz_1[i] * pc_z[i];

        ta2_zz_xxz_xyzz_0[i] = 2.0 * ta2_zz_xx_xyz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xyz_1[i] * fe_0 + 2.0 * ta1_z_xx_xyzz_1[i] + ta2_zz_xx_xyzz_0[i] * pa_z[i] - ta2_zz_xx_xyzz_1[i] * pc_z[i];

        ta2_zz_xxz_xzzz_0[i] = 3.0 * ta2_zz_xx_xzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xzz_1[i] * fe_0 + 2.0 * ta1_z_xx_xzzz_1[i] + ta2_zz_xx_xzzz_0[i] * pa_z[i] - ta2_zz_xx_xzzz_1[i] * pc_z[i];

        ta2_zz_xxz_yyyy_0[i] = 2.0 * ta1_z_xx_yyyy_1[i] + ta2_zz_xx_yyyy_0[i] * pa_z[i] - ta2_zz_xx_yyyy_1[i] * pc_z[i];

        ta2_zz_xxz_yyyz_0[i] = ta2_zz_z_yyyz_0[i] * fe_0 - ta2_zz_z_yyyz_1[i] * fe_0 + ta2_zz_xz_yyyz_0[i] * pa_x[i] - ta2_zz_xz_yyyz_1[i] * pc_x[i];

        ta2_zz_xxz_yyzz_0[i] = ta2_zz_z_yyzz_0[i] * fe_0 - ta2_zz_z_yyzz_1[i] * fe_0 + ta2_zz_xz_yyzz_0[i] * pa_x[i] - ta2_zz_xz_yyzz_1[i] * pc_x[i];

        ta2_zz_xxz_yzzz_0[i] = ta2_zz_z_yzzz_0[i] * fe_0 - ta2_zz_z_yzzz_1[i] * fe_0 + ta2_zz_xz_yzzz_0[i] * pa_x[i] - ta2_zz_xz_yzzz_1[i] * pc_x[i];

        ta2_zz_xxz_zzzz_0[i] = ta2_zz_z_zzzz_0[i] * fe_0 - ta2_zz_z_zzzz_1[i] * fe_0 + ta2_zz_xz_zzzz_0[i] * pa_x[i] - ta2_zz_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 795-810 components of targeted buffer : FG

    auto ta2_zz_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 795);

    auto ta2_zz_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 796);

    auto ta2_zz_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 797);

    auto ta2_zz_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 798);

    auto ta2_zz_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 799);

    auto ta2_zz_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 800);

    auto ta2_zz_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 801);

    auto ta2_zz_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 802);

    auto ta2_zz_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 803);

    auto ta2_zz_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 804);

    auto ta2_zz_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 805);

    auto ta2_zz_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 806);

    auto ta2_zz_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 807);

    auto ta2_zz_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 808);

    auto ta2_zz_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 809);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xyy_xxxx_0, ta2_zz_xyy_xxxy_0, ta2_zz_xyy_xxxz_0, ta2_zz_xyy_xxyy_0, ta2_zz_xyy_xxyz_0, ta2_zz_xyy_xxzz_0, ta2_zz_xyy_xyyy_0, ta2_zz_xyy_xyyz_0, ta2_zz_xyy_xyzz_0, ta2_zz_xyy_xzzz_0, ta2_zz_xyy_yyyy_0, ta2_zz_xyy_yyyz_0, ta2_zz_xyy_yyzz_0, ta2_zz_xyy_yzzz_0, ta2_zz_xyy_zzzz_0, ta2_zz_yy_xxx_0, ta2_zz_yy_xxx_1, ta2_zz_yy_xxxx_0, ta2_zz_yy_xxxx_1, ta2_zz_yy_xxxy_0, ta2_zz_yy_xxxy_1, ta2_zz_yy_xxxz_0, ta2_zz_yy_xxxz_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xxyy_0, ta2_zz_yy_xxyy_1, ta2_zz_yy_xxyz_0, ta2_zz_yy_xxyz_1, ta2_zz_yy_xxz_0, ta2_zz_yy_xxz_1, ta2_zz_yy_xxzz_0, ta2_zz_yy_xxzz_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyyy_0, ta2_zz_yy_xyyy_1, ta2_zz_yy_xyyz_0, ta2_zz_yy_xyyz_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_xyzz_0, ta2_zz_yy_xyzz_1, ta2_zz_yy_xzz_0, ta2_zz_yy_xzz_1, ta2_zz_yy_xzzz_0, ta2_zz_yy_xzzz_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyyy_0, ta2_zz_yy_yyyy_1, ta2_zz_yy_yyyz_0, ta2_zz_yy_yyyz_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yyzz_0, ta2_zz_yy_yyzz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_yzzz_0, ta2_zz_yy_yzzz_1, ta2_zz_yy_zzz_0, ta2_zz_yy_zzz_1, ta2_zz_yy_zzzz_0, ta2_zz_yy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyy_xxxx_0[i] = 4.0 * ta2_zz_yy_xxx_0[i] * fe_0 - 4.0 * ta2_zz_yy_xxx_1[i] * fe_0 + ta2_zz_yy_xxxx_0[i] * pa_x[i] - ta2_zz_yy_xxxx_1[i] * pc_x[i];

        ta2_zz_xyy_xxxy_0[i] = 3.0 * ta2_zz_yy_xxy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxy_1[i] * fe_0 + ta2_zz_yy_xxxy_0[i] * pa_x[i] - ta2_zz_yy_xxxy_1[i] * pc_x[i];

        ta2_zz_xyy_xxxz_0[i] = 3.0 * ta2_zz_yy_xxz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxz_1[i] * fe_0 + ta2_zz_yy_xxxz_0[i] * pa_x[i] - ta2_zz_yy_xxxz_1[i] * pc_x[i];

        ta2_zz_xyy_xxyy_0[i] = 2.0 * ta2_zz_yy_xyy_0[i] * fe_0 - 2.0 * ta2_zz_yy_xyy_1[i] * fe_0 + ta2_zz_yy_xxyy_0[i] * pa_x[i] - ta2_zz_yy_xxyy_1[i] * pc_x[i];

        ta2_zz_xyy_xxyz_0[i] = 2.0 * ta2_zz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yy_xyz_1[i] * fe_0 + ta2_zz_yy_xxyz_0[i] * pa_x[i] - ta2_zz_yy_xxyz_1[i] * pc_x[i];

        ta2_zz_xyy_xxzz_0[i] = 2.0 * ta2_zz_yy_xzz_0[i] * fe_0 - 2.0 * ta2_zz_yy_xzz_1[i] * fe_0 + ta2_zz_yy_xxzz_0[i] * pa_x[i] - ta2_zz_yy_xxzz_1[i] * pc_x[i];

        ta2_zz_xyy_xyyy_0[i] = ta2_zz_yy_yyy_0[i] * fe_0 - ta2_zz_yy_yyy_1[i] * fe_0 + ta2_zz_yy_xyyy_0[i] * pa_x[i] - ta2_zz_yy_xyyy_1[i] * pc_x[i];

        ta2_zz_xyy_xyyz_0[i] = ta2_zz_yy_yyz_0[i] * fe_0 - ta2_zz_yy_yyz_1[i] * fe_0 + ta2_zz_yy_xyyz_0[i] * pa_x[i] - ta2_zz_yy_xyyz_1[i] * pc_x[i];

        ta2_zz_xyy_xyzz_0[i] = ta2_zz_yy_yzz_0[i] * fe_0 - ta2_zz_yy_yzz_1[i] * fe_0 + ta2_zz_yy_xyzz_0[i] * pa_x[i] - ta2_zz_yy_xyzz_1[i] * pc_x[i];

        ta2_zz_xyy_xzzz_0[i] = ta2_zz_yy_zzz_0[i] * fe_0 - ta2_zz_yy_zzz_1[i] * fe_0 + ta2_zz_yy_xzzz_0[i] * pa_x[i] - ta2_zz_yy_xzzz_1[i] * pc_x[i];

        ta2_zz_xyy_yyyy_0[i] = ta2_zz_yy_yyyy_0[i] * pa_x[i] - ta2_zz_yy_yyyy_1[i] * pc_x[i];

        ta2_zz_xyy_yyyz_0[i] = ta2_zz_yy_yyyz_0[i] * pa_x[i] - ta2_zz_yy_yyyz_1[i] * pc_x[i];

        ta2_zz_xyy_yyzz_0[i] = ta2_zz_yy_yyzz_0[i] * pa_x[i] - ta2_zz_yy_yyzz_1[i] * pc_x[i];

        ta2_zz_xyy_yzzz_0[i] = ta2_zz_yy_yzzz_0[i] * pa_x[i] - ta2_zz_yy_yzzz_1[i] * pc_x[i];

        ta2_zz_xyy_zzzz_0[i] = ta2_zz_yy_zzzz_0[i] * pa_x[i] - ta2_zz_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 810-825 components of targeted buffer : FG

    auto ta2_zz_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 810);

    auto ta2_zz_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 811);

    auto ta2_zz_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 812);

    auto ta2_zz_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 813);

    auto ta2_zz_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 814);

    auto ta2_zz_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 815);

    auto ta2_zz_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 816);

    auto ta2_zz_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 817);

    auto ta2_zz_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 818);

    auto ta2_zz_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 819);

    auto ta2_zz_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 820);

    auto ta2_zz_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 821);

    auto ta2_zz_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 822);

    auto ta2_zz_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 823);

    auto ta2_zz_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 824);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xy_xxxy_1, ta1_z_xy_xxyy_1, ta1_z_xy_xyyy_1, ta2_zz_xy_xxxy_0, ta2_zz_xy_xxxy_1, ta2_zz_xy_xxyy_0, ta2_zz_xy_xxyy_1, ta2_zz_xy_xyyy_0, ta2_zz_xy_xyyy_1, ta2_zz_xyz_xxxx_0, ta2_zz_xyz_xxxy_0, ta2_zz_xyz_xxxz_0, ta2_zz_xyz_xxyy_0, ta2_zz_xyz_xxyz_0, ta2_zz_xyz_xxzz_0, ta2_zz_xyz_xyyy_0, ta2_zz_xyz_xyyz_0, ta2_zz_xyz_xyzz_0, ta2_zz_xyz_xzzz_0, ta2_zz_xyz_yyyy_0, ta2_zz_xyz_yyyz_0, ta2_zz_xyz_yyzz_0, ta2_zz_xyz_yzzz_0, ta2_zz_xyz_zzzz_0, ta2_zz_xz_xxxx_0, ta2_zz_xz_xxxx_1, ta2_zz_xz_xxxz_0, ta2_zz_xz_xxxz_1, ta2_zz_xz_xxzz_0, ta2_zz_xz_xxzz_1, ta2_zz_xz_xzzz_0, ta2_zz_xz_xzzz_1, ta2_zz_yz_xxyz_0, ta2_zz_yz_xxyz_1, ta2_zz_yz_xyyz_0, ta2_zz_yz_xyyz_1, ta2_zz_yz_xyz_0, ta2_zz_yz_xyz_1, ta2_zz_yz_xyzz_0, ta2_zz_yz_xyzz_1, ta2_zz_yz_yyyy_0, ta2_zz_yz_yyyy_1, ta2_zz_yz_yyyz_0, ta2_zz_yz_yyyz_1, ta2_zz_yz_yyz_0, ta2_zz_yz_yyz_1, ta2_zz_yz_yyzz_0, ta2_zz_yz_yyzz_1, ta2_zz_yz_yzz_0, ta2_zz_yz_yzz_1, ta2_zz_yz_yzzz_0, ta2_zz_yz_yzzz_1, ta2_zz_yz_zzzz_0, ta2_zz_yz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyz_xxxx_0[i] = ta2_zz_xz_xxxx_0[i] * pa_y[i] - ta2_zz_xz_xxxx_1[i] * pc_y[i];

        ta2_zz_xyz_xxxy_0[i] = 2.0 * ta1_z_xy_xxxy_1[i] + ta2_zz_xy_xxxy_0[i] * pa_z[i] - ta2_zz_xy_xxxy_1[i] * pc_z[i];

        ta2_zz_xyz_xxxz_0[i] = ta2_zz_xz_xxxz_0[i] * pa_y[i] - ta2_zz_xz_xxxz_1[i] * pc_y[i];

        ta2_zz_xyz_xxyy_0[i] = 2.0 * ta1_z_xy_xxyy_1[i] + ta2_zz_xy_xxyy_0[i] * pa_z[i] - ta2_zz_xy_xxyy_1[i] * pc_z[i];

        ta2_zz_xyz_xxyz_0[i] = 2.0 * ta2_zz_yz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xyz_1[i] * fe_0 + ta2_zz_yz_xxyz_0[i] * pa_x[i] - ta2_zz_yz_xxyz_1[i] * pc_x[i];

        ta2_zz_xyz_xxzz_0[i] = ta2_zz_xz_xxzz_0[i] * pa_y[i] - ta2_zz_xz_xxzz_1[i] * pc_y[i];

        ta2_zz_xyz_xyyy_0[i] = 2.0 * ta1_z_xy_xyyy_1[i] + ta2_zz_xy_xyyy_0[i] * pa_z[i] - ta2_zz_xy_xyyy_1[i] * pc_z[i];

        ta2_zz_xyz_xyyz_0[i] = ta2_zz_yz_yyz_0[i] * fe_0 - ta2_zz_yz_yyz_1[i] * fe_0 + ta2_zz_yz_xyyz_0[i] * pa_x[i] - ta2_zz_yz_xyyz_1[i] * pc_x[i];

        ta2_zz_xyz_xyzz_0[i] = ta2_zz_yz_yzz_0[i] * fe_0 - ta2_zz_yz_yzz_1[i] * fe_0 + ta2_zz_yz_xyzz_0[i] * pa_x[i] - ta2_zz_yz_xyzz_1[i] * pc_x[i];

        ta2_zz_xyz_xzzz_0[i] = ta2_zz_xz_xzzz_0[i] * pa_y[i] - ta2_zz_xz_xzzz_1[i] * pc_y[i];

        ta2_zz_xyz_yyyy_0[i] = ta2_zz_yz_yyyy_0[i] * pa_x[i] - ta2_zz_yz_yyyy_1[i] * pc_x[i];

        ta2_zz_xyz_yyyz_0[i] = ta2_zz_yz_yyyz_0[i] * pa_x[i] - ta2_zz_yz_yyyz_1[i] * pc_x[i];

        ta2_zz_xyz_yyzz_0[i] = ta2_zz_yz_yyzz_0[i] * pa_x[i] - ta2_zz_yz_yyzz_1[i] * pc_x[i];

        ta2_zz_xyz_yzzz_0[i] = ta2_zz_yz_yzzz_0[i] * pa_x[i] - ta2_zz_yz_yzzz_1[i] * pc_x[i];

        ta2_zz_xyz_zzzz_0[i] = ta2_zz_yz_zzzz_0[i] * pa_x[i] - ta2_zz_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 825-840 components of targeted buffer : FG

    auto ta2_zz_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 825);

    auto ta2_zz_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 826);

    auto ta2_zz_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 827);

    auto ta2_zz_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 828);

    auto ta2_zz_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 829);

    auto ta2_zz_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 830);

    auto ta2_zz_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 831);

    auto ta2_zz_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 832);

    auto ta2_zz_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 833);

    auto ta2_zz_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 834);

    auto ta2_zz_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 835);

    auto ta2_zz_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 836);

    auto ta2_zz_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 837);

    auto ta2_zz_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 838);

    auto ta2_zz_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 839);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xzz_xxxx_0, ta2_zz_xzz_xxxy_0, ta2_zz_xzz_xxxz_0, ta2_zz_xzz_xxyy_0, ta2_zz_xzz_xxyz_0, ta2_zz_xzz_xxzz_0, ta2_zz_xzz_xyyy_0, ta2_zz_xzz_xyyz_0, ta2_zz_xzz_xyzz_0, ta2_zz_xzz_xzzz_0, ta2_zz_xzz_yyyy_0, ta2_zz_xzz_yyyz_0, ta2_zz_xzz_yyzz_0, ta2_zz_xzz_yzzz_0, ta2_zz_xzz_zzzz_0, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxxx_0, ta2_zz_zz_xxxx_1, ta2_zz_zz_xxxy_0, ta2_zz_zz_xxxy_1, ta2_zz_zz_xxxz_0, ta2_zz_zz_xxxz_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxyy_0, ta2_zz_zz_xxyy_1, ta2_zz_zz_xxyz_0, ta2_zz_zz_xxyz_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xxzz_0, ta2_zz_zz_xxzz_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyyy_0, ta2_zz_zz_xyyy_1, ta2_zz_zz_xyyz_0, ta2_zz_zz_xyyz_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xyzz_0, ta2_zz_zz_xyzz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_xzzz_0, ta2_zz_zz_xzzz_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyyy_0, ta2_zz_zz_yyyy_1, ta2_zz_zz_yyyz_0, ta2_zz_zz_yyyz_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yyzz_0, ta2_zz_zz_yyzz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_yzzz_0, ta2_zz_zz_yzzz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, ta2_zz_zz_zzzz_0, ta2_zz_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzz_xxxx_0[i] = 4.0 * ta2_zz_zz_xxx_0[i] * fe_0 - 4.0 * ta2_zz_zz_xxx_1[i] * fe_0 + ta2_zz_zz_xxxx_0[i] * pa_x[i] - ta2_zz_zz_xxxx_1[i] * pc_x[i];

        ta2_zz_xzz_xxxy_0[i] = 3.0 * ta2_zz_zz_xxy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxy_1[i] * fe_0 + ta2_zz_zz_xxxy_0[i] * pa_x[i] - ta2_zz_zz_xxxy_1[i] * pc_x[i];

        ta2_zz_xzz_xxxz_0[i] = 3.0 * ta2_zz_zz_xxz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxz_1[i] * fe_0 + ta2_zz_zz_xxxz_0[i] * pa_x[i] - ta2_zz_zz_xxxz_1[i] * pc_x[i];

        ta2_zz_xzz_xxyy_0[i] = 2.0 * ta2_zz_zz_xyy_0[i] * fe_0 - 2.0 * ta2_zz_zz_xyy_1[i] * fe_0 + ta2_zz_zz_xxyy_0[i] * pa_x[i] - ta2_zz_zz_xxyy_1[i] * pc_x[i];

        ta2_zz_xzz_xxyz_0[i] = 2.0 * ta2_zz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xyz_1[i] * fe_0 + ta2_zz_zz_xxyz_0[i] * pa_x[i] - ta2_zz_zz_xxyz_1[i] * pc_x[i];

        ta2_zz_xzz_xxzz_0[i] = 2.0 * ta2_zz_zz_xzz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xzz_1[i] * fe_0 + ta2_zz_zz_xxzz_0[i] * pa_x[i] - ta2_zz_zz_xxzz_1[i] * pc_x[i];

        ta2_zz_xzz_xyyy_0[i] = ta2_zz_zz_yyy_0[i] * fe_0 - ta2_zz_zz_yyy_1[i] * fe_0 + ta2_zz_zz_xyyy_0[i] * pa_x[i] - ta2_zz_zz_xyyy_1[i] * pc_x[i];

        ta2_zz_xzz_xyyz_0[i] = ta2_zz_zz_yyz_0[i] * fe_0 - ta2_zz_zz_yyz_1[i] * fe_0 + ta2_zz_zz_xyyz_0[i] * pa_x[i] - ta2_zz_zz_xyyz_1[i] * pc_x[i];

        ta2_zz_xzz_xyzz_0[i] = ta2_zz_zz_yzz_0[i] * fe_0 - ta2_zz_zz_yzz_1[i] * fe_0 + ta2_zz_zz_xyzz_0[i] * pa_x[i] - ta2_zz_zz_xyzz_1[i] * pc_x[i];

        ta2_zz_xzz_xzzz_0[i] = ta2_zz_zz_zzz_0[i] * fe_0 - ta2_zz_zz_zzz_1[i] * fe_0 + ta2_zz_zz_xzzz_0[i] * pa_x[i] - ta2_zz_zz_xzzz_1[i] * pc_x[i];

        ta2_zz_xzz_yyyy_0[i] = ta2_zz_zz_yyyy_0[i] * pa_x[i] - ta2_zz_zz_yyyy_1[i] * pc_x[i];

        ta2_zz_xzz_yyyz_0[i] = ta2_zz_zz_yyyz_0[i] * pa_x[i] - ta2_zz_zz_yyyz_1[i] * pc_x[i];

        ta2_zz_xzz_yyzz_0[i] = ta2_zz_zz_yyzz_0[i] * pa_x[i] - ta2_zz_zz_yyzz_1[i] * pc_x[i];

        ta2_zz_xzz_yzzz_0[i] = ta2_zz_zz_yzzz_0[i] * pa_x[i] - ta2_zz_zz_yzzz_1[i] * pc_x[i];

        ta2_zz_xzz_zzzz_0[i] = ta2_zz_zz_zzzz_0[i] * pa_x[i] - ta2_zz_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 840-855 components of targeted buffer : FG

    auto ta2_zz_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 840);

    auto ta2_zz_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 841);

    auto ta2_zz_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 842);

    auto ta2_zz_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 843);

    auto ta2_zz_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 844);

    auto ta2_zz_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 845);

    auto ta2_zz_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 846);

    auto ta2_zz_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 847);

    auto ta2_zz_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 848);

    auto ta2_zz_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 849);

    auto ta2_zz_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 850);

    auto ta2_zz_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 851);

    auto ta2_zz_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 852);

    auto ta2_zz_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 853);

    auto ta2_zz_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 854);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_y_xxxx_0, ta2_zz_y_xxxx_1, ta2_zz_y_xxxy_0, ta2_zz_y_xxxy_1, ta2_zz_y_xxxz_0, ta2_zz_y_xxxz_1, ta2_zz_y_xxyy_0, ta2_zz_y_xxyy_1, ta2_zz_y_xxyz_0, ta2_zz_y_xxyz_1, ta2_zz_y_xxzz_0, ta2_zz_y_xxzz_1, ta2_zz_y_xyyy_0, ta2_zz_y_xyyy_1, ta2_zz_y_xyyz_0, ta2_zz_y_xyyz_1, ta2_zz_y_xyzz_0, ta2_zz_y_xyzz_1, ta2_zz_y_xzzz_0, ta2_zz_y_xzzz_1, ta2_zz_y_yyyy_0, ta2_zz_y_yyyy_1, ta2_zz_y_yyyz_0, ta2_zz_y_yyyz_1, ta2_zz_y_yyzz_0, ta2_zz_y_yyzz_1, ta2_zz_y_yzzz_0, ta2_zz_y_yzzz_1, ta2_zz_y_zzzz_0, ta2_zz_y_zzzz_1, ta2_zz_yy_xxx_0, ta2_zz_yy_xxx_1, ta2_zz_yy_xxxx_0, ta2_zz_yy_xxxx_1, ta2_zz_yy_xxxy_0, ta2_zz_yy_xxxy_1, ta2_zz_yy_xxxz_0, ta2_zz_yy_xxxz_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xxyy_0, ta2_zz_yy_xxyy_1, ta2_zz_yy_xxyz_0, ta2_zz_yy_xxyz_1, ta2_zz_yy_xxz_0, ta2_zz_yy_xxz_1, ta2_zz_yy_xxzz_0, ta2_zz_yy_xxzz_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyyy_0, ta2_zz_yy_xyyy_1, ta2_zz_yy_xyyz_0, ta2_zz_yy_xyyz_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_xyzz_0, ta2_zz_yy_xyzz_1, ta2_zz_yy_xzz_0, ta2_zz_yy_xzz_1, ta2_zz_yy_xzzz_0, ta2_zz_yy_xzzz_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyyy_0, ta2_zz_yy_yyyy_1, ta2_zz_yy_yyyz_0, ta2_zz_yy_yyyz_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yyzz_0, ta2_zz_yy_yyzz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_yzzz_0, ta2_zz_yy_yzzz_1, ta2_zz_yy_zzz_0, ta2_zz_yy_zzz_1, ta2_zz_yy_zzzz_0, ta2_zz_yy_zzzz_1, ta2_zz_yyy_xxxx_0, ta2_zz_yyy_xxxy_0, ta2_zz_yyy_xxxz_0, ta2_zz_yyy_xxyy_0, ta2_zz_yyy_xxyz_0, ta2_zz_yyy_xxzz_0, ta2_zz_yyy_xyyy_0, ta2_zz_yyy_xyyz_0, ta2_zz_yyy_xyzz_0, ta2_zz_yyy_xzzz_0, ta2_zz_yyy_yyyy_0, ta2_zz_yyy_yyyz_0, ta2_zz_yyy_yyzz_0, ta2_zz_yyy_yzzz_0, ta2_zz_yyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyy_xxxx_0[i] = 2.0 * ta2_zz_y_xxxx_0[i] * fe_0 - 2.0 * ta2_zz_y_xxxx_1[i] * fe_0 + ta2_zz_yy_xxxx_0[i] * pa_y[i] - ta2_zz_yy_xxxx_1[i] * pc_y[i];

        ta2_zz_yyy_xxxy_0[i] = 2.0 * ta2_zz_y_xxxy_0[i] * fe_0 - 2.0 * ta2_zz_y_xxxy_1[i] * fe_0 + ta2_zz_yy_xxx_0[i] * fe_0 - ta2_zz_yy_xxx_1[i] * fe_0 + ta2_zz_yy_xxxy_0[i] * pa_y[i] - ta2_zz_yy_xxxy_1[i] * pc_y[i];

        ta2_zz_yyy_xxxz_0[i] = 2.0 * ta2_zz_y_xxxz_0[i] * fe_0 - 2.0 * ta2_zz_y_xxxz_1[i] * fe_0 + ta2_zz_yy_xxxz_0[i] * pa_y[i] - ta2_zz_yy_xxxz_1[i] * pc_y[i];

        ta2_zz_yyy_xxyy_0[i] = 2.0 * ta2_zz_y_xxyy_0[i] * fe_0 - 2.0 * ta2_zz_y_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_yy_xxy_0[i] * fe_0 - 2.0 * ta2_zz_yy_xxy_1[i] * fe_0 + ta2_zz_yy_xxyy_0[i] * pa_y[i] - ta2_zz_yy_xxyy_1[i] * pc_y[i];

        ta2_zz_yyy_xxyz_0[i] = 2.0 * ta2_zz_y_xxyz_0[i] * fe_0 - 2.0 * ta2_zz_y_xxyz_1[i] * fe_0 + ta2_zz_yy_xxz_0[i] * fe_0 - ta2_zz_yy_xxz_1[i] * fe_0 + ta2_zz_yy_xxyz_0[i] * pa_y[i] - ta2_zz_yy_xxyz_1[i] * pc_y[i];

        ta2_zz_yyy_xxzz_0[i] = 2.0 * ta2_zz_y_xxzz_0[i] * fe_0 - 2.0 * ta2_zz_y_xxzz_1[i] * fe_0 + ta2_zz_yy_xxzz_0[i] * pa_y[i] - ta2_zz_yy_xxzz_1[i] * pc_y[i];

        ta2_zz_yyy_xyyy_0[i] = 2.0 * ta2_zz_y_xyyy_0[i] * fe_0 - 2.0 * ta2_zz_y_xyyy_1[i] * fe_0 + 3.0 * ta2_zz_yy_xyy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xyy_1[i] * fe_0 + ta2_zz_yy_xyyy_0[i] * pa_y[i] - ta2_zz_yy_xyyy_1[i] * pc_y[i];

        ta2_zz_yyy_xyyz_0[i] = 2.0 * ta2_zz_y_xyyz_0[i] * fe_0 - 2.0 * ta2_zz_y_xyyz_1[i] * fe_0 + 2.0 * ta2_zz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yy_xyz_1[i] * fe_0 + ta2_zz_yy_xyyz_0[i] * pa_y[i] - ta2_zz_yy_xyyz_1[i] * pc_y[i];

        ta2_zz_yyy_xyzz_0[i] = 2.0 * ta2_zz_y_xyzz_0[i] * fe_0 - 2.0 * ta2_zz_y_xyzz_1[i] * fe_0 + ta2_zz_yy_xzz_0[i] * fe_0 - ta2_zz_yy_xzz_1[i] * fe_0 + ta2_zz_yy_xyzz_0[i] * pa_y[i] - ta2_zz_yy_xyzz_1[i] * pc_y[i];

        ta2_zz_yyy_xzzz_0[i] = 2.0 * ta2_zz_y_xzzz_0[i] * fe_0 - 2.0 * ta2_zz_y_xzzz_1[i] * fe_0 + ta2_zz_yy_xzzz_0[i] * pa_y[i] - ta2_zz_yy_xzzz_1[i] * pc_y[i];

        ta2_zz_yyy_yyyy_0[i] = 2.0 * ta2_zz_y_yyyy_0[i] * fe_0 - 2.0 * ta2_zz_y_yyyy_1[i] * fe_0 + 4.0 * ta2_zz_yy_yyy_0[i] * fe_0 - 4.0 * ta2_zz_yy_yyy_1[i] * fe_0 + ta2_zz_yy_yyyy_0[i] * pa_y[i] - ta2_zz_yy_yyyy_1[i] * pc_y[i];

        ta2_zz_yyy_yyyz_0[i] = 2.0 * ta2_zz_y_yyyz_0[i] * fe_0 - 2.0 * ta2_zz_y_yyyz_1[i] * fe_0 + 3.0 * ta2_zz_yy_yyz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yyz_1[i] * fe_0 + ta2_zz_yy_yyyz_0[i] * pa_y[i] - ta2_zz_yy_yyyz_1[i] * pc_y[i];

        ta2_zz_yyy_yyzz_0[i] = 2.0 * ta2_zz_y_yyzz_0[i] * fe_0 - 2.0 * ta2_zz_y_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_yy_yzz_0[i] * fe_0 - 2.0 * ta2_zz_yy_yzz_1[i] * fe_0 + ta2_zz_yy_yyzz_0[i] * pa_y[i] - ta2_zz_yy_yyzz_1[i] * pc_y[i];

        ta2_zz_yyy_yzzz_0[i] = 2.0 * ta2_zz_y_yzzz_0[i] * fe_0 - 2.0 * ta2_zz_y_yzzz_1[i] * fe_0 + ta2_zz_yy_zzz_0[i] * fe_0 - ta2_zz_yy_zzz_1[i] * fe_0 + ta2_zz_yy_yzzz_0[i] * pa_y[i] - ta2_zz_yy_yzzz_1[i] * pc_y[i];

        ta2_zz_yyy_zzzz_0[i] = 2.0 * ta2_zz_y_zzzz_0[i] * fe_0 - 2.0 * ta2_zz_y_zzzz_1[i] * fe_0 + ta2_zz_yy_zzzz_0[i] * pa_y[i] - ta2_zz_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 855-870 components of targeted buffer : FG

    auto ta2_zz_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 855);

    auto ta2_zz_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 856);

    auto ta2_zz_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 857);

    auto ta2_zz_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 858);

    auto ta2_zz_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 859);

    auto ta2_zz_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 860);

    auto ta2_zz_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 861);

    auto ta2_zz_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 862);

    auto ta2_zz_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 863);

    auto ta2_zz_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 864);

    auto ta2_zz_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 865);

    auto ta2_zz_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 866);

    auto ta2_zz_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 867);

    auto ta2_zz_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 868);

    auto ta2_zz_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 869);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yy_xxxx_1, ta1_z_yy_xxxy_1, ta1_z_yy_xxyy_1, ta1_z_yy_xxyz_1, ta1_z_yy_xyyy_1, ta1_z_yy_xyyz_1, ta1_z_yy_xyzz_1, ta1_z_yy_yyyy_1, ta1_z_yy_yyyz_1, ta1_z_yy_yyzz_1, ta1_z_yy_yzzz_1, ta2_zz_yy_xxxx_0, ta2_zz_yy_xxxx_1, ta2_zz_yy_xxxy_0, ta2_zz_yy_xxxy_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xxyy_0, ta2_zz_yy_xxyy_1, ta2_zz_yy_xxyz_0, ta2_zz_yy_xxyz_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyyy_0, ta2_zz_yy_xyyy_1, ta2_zz_yy_xyyz_0, ta2_zz_yy_xyyz_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_xyzz_0, ta2_zz_yy_xyzz_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyyy_0, ta2_zz_yy_yyyy_1, ta2_zz_yy_yyyz_0, ta2_zz_yy_yyyz_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yyzz_0, ta2_zz_yy_yyzz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_yzzz_0, ta2_zz_yy_yzzz_1, ta2_zz_yyz_xxxx_0, ta2_zz_yyz_xxxy_0, ta2_zz_yyz_xxxz_0, ta2_zz_yyz_xxyy_0, ta2_zz_yyz_xxyz_0, ta2_zz_yyz_xxzz_0, ta2_zz_yyz_xyyy_0, ta2_zz_yyz_xyyz_0, ta2_zz_yyz_xyzz_0, ta2_zz_yyz_xzzz_0, ta2_zz_yyz_yyyy_0, ta2_zz_yyz_yyyz_0, ta2_zz_yyz_yyzz_0, ta2_zz_yyz_yzzz_0, ta2_zz_yyz_zzzz_0, ta2_zz_yz_xxxz_0, ta2_zz_yz_xxxz_1, ta2_zz_yz_xxzz_0, ta2_zz_yz_xxzz_1, ta2_zz_yz_xzzz_0, ta2_zz_yz_xzzz_1, ta2_zz_yz_zzzz_0, ta2_zz_yz_zzzz_1, ta2_zz_z_xxxz_0, ta2_zz_z_xxxz_1, ta2_zz_z_xxzz_0, ta2_zz_z_xxzz_1, ta2_zz_z_xzzz_0, ta2_zz_z_xzzz_1, ta2_zz_z_zzzz_0, ta2_zz_z_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyz_xxxx_0[i] = 2.0 * ta1_z_yy_xxxx_1[i] + ta2_zz_yy_xxxx_0[i] * pa_z[i] - ta2_zz_yy_xxxx_1[i] * pc_z[i];

        ta2_zz_yyz_xxxy_0[i] = 2.0 * ta1_z_yy_xxxy_1[i] + ta2_zz_yy_xxxy_0[i] * pa_z[i] - ta2_zz_yy_xxxy_1[i] * pc_z[i];

        ta2_zz_yyz_xxxz_0[i] = ta2_zz_z_xxxz_0[i] * fe_0 - ta2_zz_z_xxxz_1[i] * fe_0 + ta2_zz_yz_xxxz_0[i] * pa_y[i] - ta2_zz_yz_xxxz_1[i] * pc_y[i];

        ta2_zz_yyz_xxyy_0[i] = 2.0 * ta1_z_yy_xxyy_1[i] + ta2_zz_yy_xxyy_0[i] * pa_z[i] - ta2_zz_yy_xxyy_1[i] * pc_z[i];

        ta2_zz_yyz_xxyz_0[i] = ta2_zz_yy_xxy_0[i] * fe_0 - ta2_zz_yy_xxy_1[i] * fe_0 + 2.0 * ta1_z_yy_xxyz_1[i] + ta2_zz_yy_xxyz_0[i] * pa_z[i] - ta2_zz_yy_xxyz_1[i] * pc_z[i];

        ta2_zz_yyz_xxzz_0[i] = ta2_zz_z_xxzz_0[i] * fe_0 - ta2_zz_z_xxzz_1[i] * fe_0 + ta2_zz_yz_xxzz_0[i] * pa_y[i] - ta2_zz_yz_xxzz_1[i] * pc_y[i];

        ta2_zz_yyz_xyyy_0[i] = 2.0 * ta1_z_yy_xyyy_1[i] + ta2_zz_yy_xyyy_0[i] * pa_z[i] - ta2_zz_yy_xyyy_1[i] * pc_z[i];

        ta2_zz_yyz_xyyz_0[i] = ta2_zz_yy_xyy_0[i] * fe_0 - ta2_zz_yy_xyy_1[i] * fe_0 + 2.0 * ta1_z_yy_xyyz_1[i] + ta2_zz_yy_xyyz_0[i] * pa_z[i] - ta2_zz_yy_xyyz_1[i] * pc_z[i];

        ta2_zz_yyz_xyzz_0[i] = 2.0 * ta2_zz_yy_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yy_xyz_1[i] * fe_0 + 2.0 * ta1_z_yy_xyzz_1[i] + ta2_zz_yy_xyzz_0[i] * pa_z[i] - ta2_zz_yy_xyzz_1[i] * pc_z[i];

        ta2_zz_yyz_xzzz_0[i] = ta2_zz_z_xzzz_0[i] * fe_0 - ta2_zz_z_xzzz_1[i] * fe_0 + ta2_zz_yz_xzzz_0[i] * pa_y[i] - ta2_zz_yz_xzzz_1[i] * pc_y[i];

        ta2_zz_yyz_yyyy_0[i] = 2.0 * ta1_z_yy_yyyy_1[i] + ta2_zz_yy_yyyy_0[i] * pa_z[i] - ta2_zz_yy_yyyy_1[i] * pc_z[i];

        ta2_zz_yyz_yyyz_0[i] = ta2_zz_yy_yyy_0[i] * fe_0 - ta2_zz_yy_yyy_1[i] * fe_0 + 2.0 * ta1_z_yy_yyyz_1[i] + ta2_zz_yy_yyyz_0[i] * pa_z[i] - ta2_zz_yy_yyyz_1[i] * pc_z[i];

        ta2_zz_yyz_yyzz_0[i] = 2.0 * ta2_zz_yy_yyz_0[i] * fe_0 - 2.0 * ta2_zz_yy_yyz_1[i] * fe_0 + 2.0 * ta1_z_yy_yyzz_1[i] + ta2_zz_yy_yyzz_0[i] * pa_z[i] - ta2_zz_yy_yyzz_1[i] * pc_z[i];

        ta2_zz_yyz_yzzz_0[i] = 3.0 * ta2_zz_yy_yzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yzz_1[i] * fe_0 + 2.0 * ta1_z_yy_yzzz_1[i] + ta2_zz_yy_yzzz_0[i] * pa_z[i] - ta2_zz_yy_yzzz_1[i] * pc_z[i];

        ta2_zz_yyz_zzzz_0[i] = ta2_zz_z_zzzz_0[i] * fe_0 - ta2_zz_z_zzzz_1[i] * fe_0 + ta2_zz_yz_zzzz_0[i] * pa_y[i] - ta2_zz_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 870-885 components of targeted buffer : FG

    auto ta2_zz_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 870);

    auto ta2_zz_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 871);

    auto ta2_zz_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 872);

    auto ta2_zz_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 873);

    auto ta2_zz_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 874);

    auto ta2_zz_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 875);

    auto ta2_zz_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 876);

    auto ta2_zz_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 877);

    auto ta2_zz_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 878);

    auto ta2_zz_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 879);

    auto ta2_zz_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 880);

    auto ta2_zz_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 881);

    auto ta2_zz_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 882);

    auto ta2_zz_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 883);

    auto ta2_zz_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 884);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_yzz_xxxx_0, ta2_zz_yzz_xxxy_0, ta2_zz_yzz_xxxz_0, ta2_zz_yzz_xxyy_0, ta2_zz_yzz_xxyz_0, ta2_zz_yzz_xxzz_0, ta2_zz_yzz_xyyy_0, ta2_zz_yzz_xyyz_0, ta2_zz_yzz_xyzz_0, ta2_zz_yzz_xzzz_0, ta2_zz_yzz_yyyy_0, ta2_zz_yzz_yyyz_0, ta2_zz_yzz_yyzz_0, ta2_zz_yzz_yzzz_0, ta2_zz_yzz_zzzz_0, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxxx_0, ta2_zz_zz_xxxx_1, ta2_zz_zz_xxxy_0, ta2_zz_zz_xxxy_1, ta2_zz_zz_xxxz_0, ta2_zz_zz_xxxz_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxyy_0, ta2_zz_zz_xxyy_1, ta2_zz_zz_xxyz_0, ta2_zz_zz_xxyz_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xxzz_0, ta2_zz_zz_xxzz_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyyy_0, ta2_zz_zz_xyyy_1, ta2_zz_zz_xyyz_0, ta2_zz_zz_xyyz_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xyzz_0, ta2_zz_zz_xyzz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_xzzz_0, ta2_zz_zz_xzzz_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyyy_0, ta2_zz_zz_yyyy_1, ta2_zz_zz_yyyz_0, ta2_zz_zz_yyyz_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yyzz_0, ta2_zz_zz_yyzz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_yzzz_0, ta2_zz_zz_yzzz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, ta2_zz_zz_zzzz_0, ta2_zz_zz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzz_xxxx_0[i] = ta2_zz_zz_xxxx_0[i] * pa_y[i] - ta2_zz_zz_xxxx_1[i] * pc_y[i];

        ta2_zz_yzz_xxxy_0[i] = ta2_zz_zz_xxx_0[i] * fe_0 - ta2_zz_zz_xxx_1[i] * fe_0 + ta2_zz_zz_xxxy_0[i] * pa_y[i] - ta2_zz_zz_xxxy_1[i] * pc_y[i];

        ta2_zz_yzz_xxxz_0[i] = ta2_zz_zz_xxxz_0[i] * pa_y[i] - ta2_zz_zz_xxxz_1[i] * pc_y[i];

        ta2_zz_yzz_xxyy_0[i] = 2.0 * ta2_zz_zz_xxy_0[i] * fe_0 - 2.0 * ta2_zz_zz_xxy_1[i] * fe_0 + ta2_zz_zz_xxyy_0[i] * pa_y[i] - ta2_zz_zz_xxyy_1[i] * pc_y[i];

        ta2_zz_yzz_xxyz_0[i] = ta2_zz_zz_xxz_0[i] * fe_0 - ta2_zz_zz_xxz_1[i] * fe_0 + ta2_zz_zz_xxyz_0[i] * pa_y[i] - ta2_zz_zz_xxyz_1[i] * pc_y[i];

        ta2_zz_yzz_xxzz_0[i] = ta2_zz_zz_xxzz_0[i] * pa_y[i] - ta2_zz_zz_xxzz_1[i] * pc_y[i];

        ta2_zz_yzz_xyyy_0[i] = 3.0 * ta2_zz_zz_xyy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xyy_1[i] * fe_0 + ta2_zz_zz_xyyy_0[i] * pa_y[i] - ta2_zz_zz_xyyy_1[i] * pc_y[i];

        ta2_zz_yzz_xyyz_0[i] = 2.0 * ta2_zz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xyz_1[i] * fe_0 + ta2_zz_zz_xyyz_0[i] * pa_y[i] - ta2_zz_zz_xyyz_1[i] * pc_y[i];

        ta2_zz_yzz_xyzz_0[i] = ta2_zz_zz_xzz_0[i] * fe_0 - ta2_zz_zz_xzz_1[i] * fe_0 + ta2_zz_zz_xyzz_0[i] * pa_y[i] - ta2_zz_zz_xyzz_1[i] * pc_y[i];

        ta2_zz_yzz_xzzz_0[i] = ta2_zz_zz_xzzz_0[i] * pa_y[i] - ta2_zz_zz_xzzz_1[i] * pc_y[i];

        ta2_zz_yzz_yyyy_0[i] = 4.0 * ta2_zz_zz_yyy_0[i] * fe_0 - 4.0 * ta2_zz_zz_yyy_1[i] * fe_0 + ta2_zz_zz_yyyy_0[i] * pa_y[i] - ta2_zz_zz_yyyy_1[i] * pc_y[i];

        ta2_zz_yzz_yyyz_0[i] = 3.0 * ta2_zz_zz_yyz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yyz_1[i] * fe_0 + ta2_zz_zz_yyyz_0[i] * pa_y[i] - ta2_zz_zz_yyyz_1[i] * pc_y[i];

        ta2_zz_yzz_yyzz_0[i] = 2.0 * ta2_zz_zz_yzz_0[i] * fe_0 - 2.0 * ta2_zz_zz_yzz_1[i] * fe_0 + ta2_zz_zz_yyzz_0[i] * pa_y[i] - ta2_zz_zz_yyzz_1[i] * pc_y[i];

        ta2_zz_yzz_yzzz_0[i] = ta2_zz_zz_zzz_0[i] * fe_0 - ta2_zz_zz_zzz_1[i] * fe_0 + ta2_zz_zz_yzzz_0[i] * pa_y[i] - ta2_zz_zz_yzzz_1[i] * pc_y[i];

        ta2_zz_yzz_zzzz_0[i] = ta2_zz_zz_zzzz_0[i] * pa_y[i] - ta2_zz_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 885-900 components of targeted buffer : FG

    auto ta2_zz_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 885);

    auto ta2_zz_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 886);

    auto ta2_zz_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 887);

    auto ta2_zz_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 888);

    auto ta2_zz_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 889);

    auto ta2_zz_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 890);

    auto ta2_zz_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 891);

    auto ta2_zz_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 892);

    auto ta2_zz_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 893);

    auto ta2_zz_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 894);

    auto ta2_zz_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 895);

    auto ta2_zz_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 896);

    auto ta2_zz_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 897);

    auto ta2_zz_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 898);

    auto ta2_zz_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 899);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zz_xxxx_1, ta1_z_zz_xxxy_1, ta1_z_zz_xxxz_1, ta1_z_zz_xxyy_1, ta1_z_zz_xxyz_1, ta1_z_zz_xxzz_1, ta1_z_zz_xyyy_1, ta1_z_zz_xyyz_1, ta1_z_zz_xyzz_1, ta1_z_zz_xzzz_1, ta1_z_zz_yyyy_1, ta1_z_zz_yyyz_1, ta1_z_zz_yyzz_1, ta1_z_zz_yzzz_1, ta1_z_zz_zzzz_1, ta2_zz_z_xxxx_0, ta2_zz_z_xxxx_1, ta2_zz_z_xxxy_0, ta2_zz_z_xxxy_1, ta2_zz_z_xxxz_0, ta2_zz_z_xxxz_1, ta2_zz_z_xxyy_0, ta2_zz_z_xxyy_1, ta2_zz_z_xxyz_0, ta2_zz_z_xxyz_1, ta2_zz_z_xxzz_0, ta2_zz_z_xxzz_1, ta2_zz_z_xyyy_0, ta2_zz_z_xyyy_1, ta2_zz_z_xyyz_0, ta2_zz_z_xyyz_1, ta2_zz_z_xyzz_0, ta2_zz_z_xyzz_1, ta2_zz_z_xzzz_0, ta2_zz_z_xzzz_1, ta2_zz_z_yyyy_0, ta2_zz_z_yyyy_1, ta2_zz_z_yyyz_0, ta2_zz_z_yyyz_1, ta2_zz_z_yyzz_0, ta2_zz_z_yyzz_1, ta2_zz_z_yzzz_0, ta2_zz_z_yzzz_1, ta2_zz_z_zzzz_0, ta2_zz_z_zzzz_1, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxxx_0, ta2_zz_zz_xxxx_1, ta2_zz_zz_xxxy_0, ta2_zz_zz_xxxy_1, ta2_zz_zz_xxxz_0, ta2_zz_zz_xxxz_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxyy_0, ta2_zz_zz_xxyy_1, ta2_zz_zz_xxyz_0, ta2_zz_zz_xxyz_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xxzz_0, ta2_zz_zz_xxzz_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyyy_0, ta2_zz_zz_xyyy_1, ta2_zz_zz_xyyz_0, ta2_zz_zz_xyyz_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xyzz_0, ta2_zz_zz_xyzz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_xzzz_0, ta2_zz_zz_xzzz_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyyy_0, ta2_zz_zz_yyyy_1, ta2_zz_zz_yyyz_0, ta2_zz_zz_yyyz_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yyzz_0, ta2_zz_zz_yyzz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_yzzz_0, ta2_zz_zz_yzzz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, ta2_zz_zz_zzzz_0, ta2_zz_zz_zzzz_1, ta2_zz_zzz_xxxx_0, ta2_zz_zzz_xxxy_0, ta2_zz_zzz_xxxz_0, ta2_zz_zzz_xxyy_0, ta2_zz_zzz_xxyz_0, ta2_zz_zzz_xxzz_0, ta2_zz_zzz_xyyy_0, ta2_zz_zzz_xyyz_0, ta2_zz_zzz_xyzz_0, ta2_zz_zzz_xzzz_0, ta2_zz_zzz_yyyy_0, ta2_zz_zzz_yyyz_0, ta2_zz_zzz_yyzz_0, ta2_zz_zzz_yzzz_0, ta2_zz_zzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzz_xxxx_0[i] = 2.0 * ta2_zz_z_xxxx_0[i] * fe_0 - 2.0 * ta2_zz_z_xxxx_1[i] * fe_0 + 2.0 * ta1_z_zz_xxxx_1[i] + ta2_zz_zz_xxxx_0[i] * pa_z[i] - ta2_zz_zz_xxxx_1[i] * pc_z[i];

        ta2_zz_zzz_xxxy_0[i] = 2.0 * ta2_zz_z_xxxy_0[i] * fe_0 - 2.0 * ta2_zz_z_xxxy_1[i] * fe_0 + 2.0 * ta1_z_zz_xxxy_1[i] + ta2_zz_zz_xxxy_0[i] * pa_z[i] - ta2_zz_zz_xxxy_1[i] * pc_z[i];

        ta2_zz_zzz_xxxz_0[i] = 2.0 * ta2_zz_z_xxxz_0[i] * fe_0 - 2.0 * ta2_zz_z_xxxz_1[i] * fe_0 + ta2_zz_zz_xxx_0[i] * fe_0 - ta2_zz_zz_xxx_1[i] * fe_0 + 2.0 * ta1_z_zz_xxxz_1[i] + ta2_zz_zz_xxxz_0[i] * pa_z[i] - ta2_zz_zz_xxxz_1[i] * pc_z[i];

        ta2_zz_zzz_xxyy_0[i] = 2.0 * ta2_zz_z_xxyy_0[i] * fe_0 - 2.0 * ta2_zz_z_xxyy_1[i] * fe_0 + 2.0 * ta1_z_zz_xxyy_1[i] + ta2_zz_zz_xxyy_0[i] * pa_z[i] - ta2_zz_zz_xxyy_1[i] * pc_z[i];

        ta2_zz_zzz_xxyz_0[i] = 2.0 * ta2_zz_z_xxyz_0[i] * fe_0 - 2.0 * ta2_zz_z_xxyz_1[i] * fe_0 + ta2_zz_zz_xxy_0[i] * fe_0 - ta2_zz_zz_xxy_1[i] * fe_0 + 2.0 * ta1_z_zz_xxyz_1[i] + ta2_zz_zz_xxyz_0[i] * pa_z[i] - ta2_zz_zz_xxyz_1[i] * pc_z[i];

        ta2_zz_zzz_xxzz_0[i] = 2.0 * ta2_zz_z_xxzz_0[i] * fe_0 - 2.0 * ta2_zz_z_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_zz_xxz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xxz_1[i] * fe_0 + 2.0 * ta1_z_zz_xxzz_1[i] + ta2_zz_zz_xxzz_0[i] * pa_z[i] - ta2_zz_zz_xxzz_1[i] * pc_z[i];

        ta2_zz_zzz_xyyy_0[i] = 2.0 * ta2_zz_z_xyyy_0[i] * fe_0 - 2.0 * ta2_zz_z_xyyy_1[i] * fe_0 + 2.0 * ta1_z_zz_xyyy_1[i] + ta2_zz_zz_xyyy_0[i] * pa_z[i] - ta2_zz_zz_xyyy_1[i] * pc_z[i];

        ta2_zz_zzz_xyyz_0[i] = 2.0 * ta2_zz_z_xyyz_0[i] * fe_0 - 2.0 * ta2_zz_z_xyyz_1[i] * fe_0 + ta2_zz_zz_xyy_0[i] * fe_0 - ta2_zz_zz_xyy_1[i] * fe_0 + 2.0 * ta1_z_zz_xyyz_1[i] + ta2_zz_zz_xyyz_0[i] * pa_z[i] - ta2_zz_zz_xyyz_1[i] * pc_z[i];

        ta2_zz_zzz_xyzz_0[i] = 2.0 * ta2_zz_z_xyzz_0[i] * fe_0 - 2.0 * ta2_zz_z_xyzz_1[i] * fe_0 + 2.0 * ta2_zz_zz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xyz_1[i] * fe_0 + 2.0 * ta1_z_zz_xyzz_1[i] + ta2_zz_zz_xyzz_0[i] * pa_z[i] - ta2_zz_zz_xyzz_1[i] * pc_z[i];

        ta2_zz_zzz_xzzz_0[i] = 2.0 * ta2_zz_z_xzzz_0[i] * fe_0 - 2.0 * ta2_zz_z_xzzz_1[i] * fe_0 + 3.0 * ta2_zz_zz_xzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xzzz_1[i] + ta2_zz_zz_xzzz_0[i] * pa_z[i] - ta2_zz_zz_xzzz_1[i] * pc_z[i];

        ta2_zz_zzz_yyyy_0[i] = 2.0 * ta2_zz_z_yyyy_0[i] * fe_0 - 2.0 * ta2_zz_z_yyyy_1[i] * fe_0 + 2.0 * ta1_z_zz_yyyy_1[i] + ta2_zz_zz_yyyy_0[i] * pa_z[i] - ta2_zz_zz_yyyy_1[i] * pc_z[i];

        ta2_zz_zzz_yyyz_0[i] = 2.0 * ta2_zz_z_yyyz_0[i] * fe_0 - 2.0 * ta2_zz_z_yyyz_1[i] * fe_0 + ta2_zz_zz_yyy_0[i] * fe_0 - ta2_zz_zz_yyy_1[i] * fe_0 + 2.0 * ta1_z_zz_yyyz_1[i] + ta2_zz_zz_yyyz_0[i] * pa_z[i] - ta2_zz_zz_yyyz_1[i] * pc_z[i];

        ta2_zz_zzz_yyzz_0[i] = 2.0 * ta2_zz_z_yyzz_0[i] * fe_0 - 2.0 * ta2_zz_z_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_zz_yyz_0[i] * fe_0 - 2.0 * ta2_zz_zz_yyz_1[i] * fe_0 + 2.0 * ta1_z_zz_yyzz_1[i] + ta2_zz_zz_yyzz_0[i] * pa_z[i] - ta2_zz_zz_yyzz_1[i] * pc_z[i];

        ta2_zz_zzz_yzzz_0[i] = 2.0 * ta2_zz_z_yzzz_0[i] * fe_0 - 2.0 * ta2_zz_z_yzzz_1[i] * fe_0 + 3.0 * ta2_zz_zz_yzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yzz_1[i] * fe_0 + 2.0 * ta1_z_zz_yzzz_1[i] + ta2_zz_zz_yzzz_0[i] * pa_z[i] - ta2_zz_zz_yzzz_1[i] * pc_z[i];

        ta2_zz_zzz_zzzz_0[i] = 2.0 * ta2_zz_z_zzzz_0[i] * fe_0 - 2.0 * ta2_zz_z_zzzz_1[i] * fe_0 + 4.0 * ta2_zz_zz_zzz_0[i] * fe_0 - 4.0 * ta2_zz_zz_zzz_1[i] * fe_0 + 2.0 * ta1_z_zz_zzzz_1[i] + ta2_zz_zz_zzzz_0[i] * pa_z[i] - ta2_zz_zz_zzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

