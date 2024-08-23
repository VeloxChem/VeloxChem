#include "NuclearPotentialGeom010PrimRecIG.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_ig(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_ig,
                                        const size_t idx_npot_geom_010_0_gg,
                                        const size_t idx_npot_geom_010_1_gg,
                                        const size_t idx_npot_geom_010_0_hf,
                                        const size_t idx_npot_geom_010_1_hf,
                                        const size_t idx_npot_1_hg,
                                        const size_t idx_npot_geom_010_0_hg,
                                        const size_t idx_npot_geom_010_1_hg,
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

    auto ta1_x_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 95);

    auto ta1_x_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 96);

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

    auto ta1_x_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 140);

    auto ta1_x_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 141);

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

    auto ta1_y_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 245);

    auto ta1_y_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 246);

    auto ta1_y_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 249);

    auto ta1_y_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 250);

    auto ta1_y_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 251);

    auto ta1_y_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 252);

    auto ta1_y_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 253);

    auto ta1_y_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 255);

    auto ta1_y_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 256);

    auto ta1_y_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 258);

    auto ta1_y_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 261);

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

    auto ta1_y_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 393);

    auto ta1_y_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 394);

    auto ta1_y_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 396);

    auto ta1_y_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 397);

    auto ta1_y_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 398);

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

    auto ta1_y_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 425);

    auto ta1_y_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 426);

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

    auto ta1_z_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 467);

    auto ta1_z_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 470);

    auto ta1_z_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 474);

    auto ta1_z_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 475);

    auto ta1_z_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 476);

    auto ta1_z_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 477);

    auto ta1_z_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 478);

    auto ta1_z_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_gg + 480);

    auto ta1_z_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 481);

    auto ta1_z_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 482);

    auto ta1_z_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 483);

    auto ta1_z_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 485);

    auto ta1_z_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 486);

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

    auto ta1_z_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 620);

    auto ta1_z_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_gg + 621);

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

    auto ta1_z_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 647);

    auto ta1_z_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 649);

    auto ta1_z_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_gg + 650);

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

    auto ta1_x_xyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 95);

    auto ta1_x_xyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 96);

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

    auto ta1_x_xzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 140);

    auto ta1_x_xzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 141);

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

    auto ta1_y_xxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 245);

    auto ta1_y_xxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 246);

    auto ta1_y_xxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 249);

    auto ta1_y_xxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 250);

    auto ta1_y_xxxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 251);

    auto ta1_y_xxxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 252);

    auto ta1_y_xxxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 253);

    auto ta1_y_xxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 255);

    auto ta1_y_xxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 256);

    auto ta1_y_xxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 258);

    auto ta1_y_xxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 261);

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

    auto ta1_y_yyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 393);

    auto ta1_y_yyyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 394);

    auto ta1_y_yyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 396);

    auto ta1_y_yyyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 397);

    auto ta1_y_yyyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 398);

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

    auto ta1_y_yzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 425);

    auto ta1_y_yzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 426);

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

    auto ta1_z_xxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 467);

    auto ta1_z_xxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 470);

    auto ta1_z_xxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 474);

    auto ta1_z_xxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 475);

    auto ta1_z_xxxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 476);

    auto ta1_z_xxxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 477);

    auto ta1_z_xxxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 478);

    auto ta1_z_xxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_gg + 480);

    auto ta1_z_xxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 481);

    auto ta1_z_xxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 482);

    auto ta1_z_xxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 483);

    auto ta1_z_xxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 485);

    auto ta1_z_xxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 486);

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

    auto ta1_z_yyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 620);

    auto ta1_z_yyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_gg + 621);

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

    auto ta1_z_yzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 647);

    auto ta1_z_yzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 649);

    auto ta1_z_yzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_gg + 650);

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

    // Set up components of auxiliary buffer : HF

    auto ta1_x_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf);

    auto ta1_x_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 1);

    auto ta1_x_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 2);

    auto ta1_x_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 3);

    auto ta1_x_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 4);

    auto ta1_x_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 5);

    auto ta1_x_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 6);

    auto ta1_x_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 7);

    auto ta1_x_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 8);

    auto ta1_x_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 9);

    auto ta1_x_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 10);

    auto ta1_x_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 11);

    auto ta1_x_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 12);

    auto ta1_x_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 13);

    auto ta1_x_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 14);

    auto ta1_x_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 15);

    auto ta1_x_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 20);

    auto ta1_x_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 21);

    auto ta1_x_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 22);

    auto ta1_x_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 23);

    auto ta1_x_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 24);

    auto ta1_x_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 25);

    auto ta1_x_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 27);

    auto ta1_x_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 28);

    auto ta1_x_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 29);

    auto ta1_x_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 30);

    auto ta1_x_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 31);

    auto ta1_x_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 32);

    auto ta1_x_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 33);

    auto ta1_x_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 34);

    auto ta1_x_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 35);

    auto ta1_x_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 36);

    auto ta1_x_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 37);

    auto ta1_x_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 38);

    auto ta1_x_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 50);

    auto ta1_x_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 51);

    auto ta1_x_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 52);

    auto ta1_x_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 53);

    auto ta1_x_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 54);

    auto ta1_x_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 55);

    auto ta1_x_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 56);

    auto ta1_x_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 57);

    auto ta1_x_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 58);

    auto ta1_x_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 59);

    auto ta1_x_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 60);

    auto ta1_x_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 61);

    auto ta1_x_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 62);

    auto ta1_x_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 63);

    auto ta1_x_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 64);

    auto ta1_x_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 65);

    auto ta1_x_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 66);

    auto ta1_x_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 67);

    auto ta1_x_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 68);

    auto ta1_x_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 82);

    auto ta1_x_xxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 84);

    auto ta1_x_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 85);

    auto ta1_x_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 90);

    auto ta1_x_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 91);

    auto ta1_x_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 92);

    auto ta1_x_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 93);

    auto ta1_x_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 94);

    auto ta1_x_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 95);

    auto ta1_x_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 96);

    auto ta1_x_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 97);

    auto ta1_x_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 98);

    auto ta1_x_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 99);

    auto ta1_x_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 101);

    auto ta1_x_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 103);

    auto ta1_x_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 104);

    auto ta1_x_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 140);

    auto ta1_x_xzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 141);

    auto ta1_x_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 142);

    auto ta1_x_xzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 143);

    auto ta1_x_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 144);

    auto ta1_x_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 145);

    auto ta1_x_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 150);

    auto ta1_x_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 151);

    auto ta1_x_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 152);

    auto ta1_x_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 153);

    auto ta1_x_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 154);

    auto ta1_x_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 155);

    auto ta1_x_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 156);

    auto ta1_x_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 157);

    auto ta1_x_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 158);

    auto ta1_x_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 159);

    auto ta1_x_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 172);

    auto ta1_x_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 174);

    auto ta1_x_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 175);

    auto ta1_x_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 177);

    auto ta1_x_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 178);

    auto ta1_x_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 179);

    auto ta1_x_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 182);

    auto ta1_x_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 184);

    auto ta1_x_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 185);

    auto ta1_x_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 187);

    auto ta1_x_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 188);

    auto ta1_x_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 189);

    auto ta1_x_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 192);

    auto ta1_x_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 194);

    auto ta1_x_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 195);

    auto ta1_x_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 197);

    auto ta1_x_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 198);

    auto ta1_x_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 199);

    auto ta1_x_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 200);

    auto ta1_x_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 201);

    auto ta1_x_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 202);

    auto ta1_x_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 203);

    auto ta1_x_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 204);

    auto ta1_x_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 205);

    auto ta1_x_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 206);

    auto ta1_x_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 207);

    auto ta1_x_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 208);

    auto ta1_x_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 209);

    auto ta1_y_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 210);

    auto ta1_y_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 211);

    auto ta1_y_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 212);

    auto ta1_y_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 213);

    auto ta1_y_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 214);

    auto ta1_y_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 215);

    auto ta1_y_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 216);

    auto ta1_y_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 217);

    auto ta1_y_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 218);

    auto ta1_y_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 219);

    auto ta1_y_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 221);

    auto ta1_y_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 223);

    auto ta1_y_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 224);

    auto ta1_y_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 240);

    auto ta1_y_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 241);

    auto ta1_y_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 242);

    auto ta1_y_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 243);

    auto ta1_y_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 244);

    auto ta1_y_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 245);

    auto ta1_y_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 246);

    auto ta1_y_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 247);

    auto ta1_y_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 248);

    auto ta1_y_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 262);

    auto ta1_y_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 264);

    auto ta1_y_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 265);

    auto ta1_y_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 267);

    auto ta1_y_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 268);

    auto ta1_y_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 269);

    auto ta1_y_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 270);

    auto ta1_y_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 271);

    auto ta1_y_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 272);

    auto ta1_y_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 273);

    auto ta1_y_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 274);

    auto ta1_y_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 275);

    auto ta1_y_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 276);

    auto ta1_y_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 277);

    auto ta1_y_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 278);

    auto ta1_y_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 302);

    auto ta1_y_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 304);

    auto ta1_y_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 305);

    auto ta1_y_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 307);

    auto ta1_y_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 308);

    auto ta1_y_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 309);

    auto ta1_y_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 311);

    auto ta1_y_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 313);

    auto ta1_y_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 314);

    auto ta1_y_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 316);

    auto ta1_y_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 317);

    auto ta1_y_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 318);

    auto ta1_y_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 334);

    auto ta1_y_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 337);

    auto ta1_y_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 338);

    auto ta1_y_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 352);

    auto ta1_y_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 354);

    auto ta1_y_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 355);

    auto ta1_y_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 357);

    auto ta1_y_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 358);

    auto ta1_y_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 359);

    auto ta1_y_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 360);

    auto ta1_y_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 361);

    auto ta1_y_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 362);

    auto ta1_y_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 363);

    auto ta1_y_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 364);

    auto ta1_y_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 365);

    auto ta1_y_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 366);

    auto ta1_y_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 367);

    auto ta1_y_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 368);

    auto ta1_y_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 369);

    auto ta1_y_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 371);

    auto ta1_y_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 372);

    auto ta1_y_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 373);

    auto ta1_y_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 374);

    auto ta1_y_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 375);

    auto ta1_y_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 376);

    auto ta1_y_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 377);

    auto ta1_y_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 378);

    auto ta1_y_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 379);

    auto ta1_y_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 380);

    auto ta1_y_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 381);

    auto ta1_y_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 382);

    auto ta1_y_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 383);

    auto ta1_y_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 384);

    auto ta1_y_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 385);

    auto ta1_y_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 386);

    auto ta1_y_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 387);

    auto ta1_y_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 388);

    auto ta1_y_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 389);

    auto ta1_y_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 390);

    auto ta1_y_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 391);

    auto ta1_y_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 392);

    auto ta1_y_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 393);

    auto ta1_y_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 394);

    auto ta1_y_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 395);

    auto ta1_y_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 396);

    auto ta1_y_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 397);

    auto ta1_y_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 398);

    auto ta1_y_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 399);

    auto ta1_y_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 401);

    auto ta1_y_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 403);

    auto ta1_y_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 404);

    auto ta1_y_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 406);

    auto ta1_y_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 407);

    auto ta1_y_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 408);

    auto ta1_y_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 410);

    auto ta1_y_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 411);

    auto ta1_y_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 412);

    auto ta1_y_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 413);

    auto ta1_y_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 414);

    auto ta1_y_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 415);

    auto ta1_y_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 416);

    auto ta1_y_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 417);

    auto ta1_y_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 418);

    auto ta1_y_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 419);

    auto ta1_z_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 420);

    auto ta1_z_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 421);

    auto ta1_z_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 422);

    auto ta1_z_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 423);

    auto ta1_z_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 424);

    auto ta1_z_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 425);

    auto ta1_z_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 426);

    auto ta1_z_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 427);

    auto ta1_z_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 428);

    auto ta1_z_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 429);

    auto ta1_z_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 442);

    auto ta1_z_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 444);

    auto ta1_z_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 445);

    auto ta1_z_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 451);

    auto ta1_z_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 453);

    auto ta1_z_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 454);

    auto ta1_z_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 456);

    auto ta1_z_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 457);

    auto ta1_z_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 458);

    auto ta1_z_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 470);

    auto ta1_z_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 471);

    auto ta1_z_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 472);

    auto ta1_z_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 473);

    auto ta1_z_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 474);

    auto ta1_z_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 475);

    auto ta1_z_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 477);

    auto ta1_z_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 478);

    auto ta1_z_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 479);

    auto ta1_z_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 481);

    auto ta1_z_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 483);

    auto ta1_z_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 484);

    auto ta1_z_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 486);

    auto ta1_z_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 487);

    auto ta1_z_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 488);

    auto ta1_z_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 510);

    auto ta1_z_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 511);

    auto ta1_z_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 512);

    auto ta1_z_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 513);

    auto ta1_z_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 514);

    auto ta1_z_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 515);

    auto ta1_z_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 517);

    auto ta1_z_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 518);

    auto ta1_z_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 519);

    auto ta1_z_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 521);

    auto ta1_z_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 523);

    auto ta1_z_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 524);

    auto ta1_z_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 526);

    auto ta1_z_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 527);

    auto ta1_z_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 528);

    auto ta1_z_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 544);

    auto ta1_z_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 547);

    auto ta1_z_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 548);

    auto ta1_z_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 562);

    auto ta1_z_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 564);

    auto ta1_z_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 565);

    auto ta1_z_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 567);

    auto ta1_z_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 568);

    auto ta1_z_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 569);

    auto ta1_z_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 570);

    auto ta1_z_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 571);

    auto ta1_z_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 572);

    auto ta1_z_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 573);

    auto ta1_z_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 574);

    auto ta1_z_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 575);

    auto ta1_z_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 576);

    auto ta1_z_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 577);

    auto ta1_z_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 578);

    auto ta1_z_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 579);

    auto ta1_z_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 582);

    auto ta1_z_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 584);

    auto ta1_z_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 585);

    auto ta1_z_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 587);

    auto ta1_z_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 588);

    auto ta1_z_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 589);

    auto ta1_z_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 590);

    auto ta1_z_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 591);

    auto ta1_z_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 592);

    auto ta1_z_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 593);

    auto ta1_z_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 594);

    auto ta1_z_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 595);

    auto ta1_z_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 596);

    auto ta1_z_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 597);

    auto ta1_z_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 598);

    auto ta1_z_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 599);

    auto ta1_z_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 600);

    auto ta1_z_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 601);

    auto ta1_z_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 602);

    auto ta1_z_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 603);

    auto ta1_z_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 604);

    auto ta1_z_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 605);

    auto ta1_z_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 606);

    auto ta1_z_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 607);

    auto ta1_z_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 608);

    auto ta1_z_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 609);

    auto ta1_z_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 611);

    auto ta1_z_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 612);

    auto ta1_z_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 613);

    auto ta1_z_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 614);

    auto ta1_z_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 615);

    auto ta1_z_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 616);

    auto ta1_z_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 617);

    auto ta1_z_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 618);

    auto ta1_z_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 619);

    auto ta1_z_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 620);

    auto ta1_z_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 621);

    auto ta1_z_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 622);

    auto ta1_z_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 623);

    auto ta1_z_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 624);

    auto ta1_z_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 625);

    auto ta1_z_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 626);

    auto ta1_z_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 627);

    auto ta1_z_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 628);

    auto ta1_z_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 629);

    // Set up components of auxiliary buffer : HF

    auto ta1_x_xxxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf);

    auto ta1_x_xxxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 1);

    auto ta1_x_xxxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 2);

    auto ta1_x_xxxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 3);

    auto ta1_x_xxxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 4);

    auto ta1_x_xxxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 5);

    auto ta1_x_xxxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 6);

    auto ta1_x_xxxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 7);

    auto ta1_x_xxxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 8);

    auto ta1_x_xxxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 9);

    auto ta1_x_xxxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 10);

    auto ta1_x_xxxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 11);

    auto ta1_x_xxxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 12);

    auto ta1_x_xxxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 13);

    auto ta1_x_xxxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 14);

    auto ta1_x_xxxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 15);

    auto ta1_x_xxxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 20);

    auto ta1_x_xxxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 21);

    auto ta1_x_xxxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 22);

    auto ta1_x_xxxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 23);

    auto ta1_x_xxxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 24);

    auto ta1_x_xxxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 25);

    auto ta1_x_xxxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 27);

    auto ta1_x_xxxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 28);

    auto ta1_x_xxxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 29);

    auto ta1_x_xxxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 30);

    auto ta1_x_xxxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 31);

    auto ta1_x_xxxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 32);

    auto ta1_x_xxxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 33);

    auto ta1_x_xxxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 34);

    auto ta1_x_xxxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 35);

    auto ta1_x_xxxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 36);

    auto ta1_x_xxxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 37);

    auto ta1_x_xxxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 38);

    auto ta1_x_xxxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 50);

    auto ta1_x_xxxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 51);

    auto ta1_x_xxxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 52);

    auto ta1_x_xxxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 53);

    auto ta1_x_xxxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 54);

    auto ta1_x_xxxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 55);

    auto ta1_x_xxxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 56);

    auto ta1_x_xxxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 57);

    auto ta1_x_xxxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 58);

    auto ta1_x_xxxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 59);

    auto ta1_x_xxyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 60);

    auto ta1_x_xxyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 61);

    auto ta1_x_xxyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 62);

    auto ta1_x_xxyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 63);

    auto ta1_x_xxyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 64);

    auto ta1_x_xxyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 65);

    auto ta1_x_xxyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 66);

    auto ta1_x_xxyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 67);

    auto ta1_x_xxyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 68);

    auto ta1_x_xxyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 82);

    auto ta1_x_xxyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 84);

    auto ta1_x_xxyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 85);

    auto ta1_x_xxzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 90);

    auto ta1_x_xxzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 91);

    auto ta1_x_xxzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 92);

    auto ta1_x_xxzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 93);

    auto ta1_x_xxzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 94);

    auto ta1_x_xxzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 95);

    auto ta1_x_xxzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 96);

    auto ta1_x_xxzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 97);

    auto ta1_x_xxzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 98);

    auto ta1_x_xxzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 99);

    auto ta1_x_xyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 101);

    auto ta1_x_xyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 103);

    auto ta1_x_xyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 104);

    auto ta1_x_xzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 140);

    auto ta1_x_xzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 141);

    auto ta1_x_xzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 142);

    auto ta1_x_xzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 143);

    auto ta1_x_xzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 144);

    auto ta1_x_xzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 145);

    auto ta1_x_yyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 150);

    auto ta1_x_yyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 151);

    auto ta1_x_yyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 152);

    auto ta1_x_yyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 153);

    auto ta1_x_yyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 154);

    auto ta1_x_yyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 155);

    auto ta1_x_yyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 156);

    auto ta1_x_yyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 157);

    auto ta1_x_yyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 158);

    auto ta1_x_yyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 159);

    auto ta1_x_yyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 172);

    auto ta1_x_yyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 174);

    auto ta1_x_yyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 175);

    auto ta1_x_yyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 177);

    auto ta1_x_yyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 178);

    auto ta1_x_yyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 179);

    auto ta1_x_yyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 182);

    auto ta1_x_yyzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 184);

    auto ta1_x_yyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 185);

    auto ta1_x_yyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 187);

    auto ta1_x_yyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 188);

    auto ta1_x_yyzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 189);

    auto ta1_x_yzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 192);

    auto ta1_x_yzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 194);

    auto ta1_x_yzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 195);

    auto ta1_x_yzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 197);

    auto ta1_x_yzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 198);

    auto ta1_x_yzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 199);

    auto ta1_x_zzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 200);

    auto ta1_x_zzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 201);

    auto ta1_x_zzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 202);

    auto ta1_x_zzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 203);

    auto ta1_x_zzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 204);

    auto ta1_x_zzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 205);

    auto ta1_x_zzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 206);

    auto ta1_x_zzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 207);

    auto ta1_x_zzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 208);

    auto ta1_x_zzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 209);

    auto ta1_y_xxxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 210);

    auto ta1_y_xxxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 211);

    auto ta1_y_xxxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 212);

    auto ta1_y_xxxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 213);

    auto ta1_y_xxxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 214);

    auto ta1_y_xxxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 215);

    auto ta1_y_xxxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 216);

    auto ta1_y_xxxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 217);

    auto ta1_y_xxxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 218);

    auto ta1_y_xxxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 219);

    auto ta1_y_xxxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 221);

    auto ta1_y_xxxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 223);

    auto ta1_y_xxxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 224);

    auto ta1_y_xxxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 240);

    auto ta1_y_xxxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 241);

    auto ta1_y_xxxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 242);

    auto ta1_y_xxxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 243);

    auto ta1_y_xxxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 244);

    auto ta1_y_xxxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 245);

    auto ta1_y_xxxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 246);

    auto ta1_y_xxxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 247);

    auto ta1_y_xxxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 248);

    auto ta1_y_xxxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 262);

    auto ta1_y_xxxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 264);

    auto ta1_y_xxxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 265);

    auto ta1_y_xxxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 267);

    auto ta1_y_xxxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 268);

    auto ta1_y_xxxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 269);

    auto ta1_y_xxyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 270);

    auto ta1_y_xxyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 271);

    auto ta1_y_xxyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 272);

    auto ta1_y_xxyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 273);

    auto ta1_y_xxyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 274);

    auto ta1_y_xxyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 275);

    auto ta1_y_xxyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 276);

    auto ta1_y_xxyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 277);

    auto ta1_y_xxyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 278);

    auto ta1_y_xxzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 302);

    auto ta1_y_xxzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 304);

    auto ta1_y_xxzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 305);

    auto ta1_y_xxzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 307);

    auto ta1_y_xxzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 308);

    auto ta1_y_xxzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 309);

    auto ta1_y_xyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 311);

    auto ta1_y_xyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 313);

    auto ta1_y_xyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 314);

    auto ta1_y_xyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 316);

    auto ta1_y_xyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 317);

    auto ta1_y_xyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 318);

    auto ta1_y_xyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 334);

    auto ta1_y_xyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 337);

    auto ta1_y_xyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 338);

    auto ta1_y_xzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 352);

    auto ta1_y_xzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 354);

    auto ta1_y_xzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 355);

    auto ta1_y_xzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 357);

    auto ta1_y_xzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 358);

    auto ta1_y_xzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 359);

    auto ta1_y_yyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 360);

    auto ta1_y_yyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 361);

    auto ta1_y_yyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 362);

    auto ta1_y_yyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 363);

    auto ta1_y_yyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 364);

    auto ta1_y_yyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 365);

    auto ta1_y_yyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 366);

    auto ta1_y_yyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 367);

    auto ta1_y_yyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 368);

    auto ta1_y_yyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 369);

    auto ta1_y_yyyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 371);

    auto ta1_y_yyyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 372);

    auto ta1_y_yyyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 373);

    auto ta1_y_yyyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 374);

    auto ta1_y_yyyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 375);

    auto ta1_y_yyyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 376);

    auto ta1_y_yyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 377);

    auto ta1_y_yyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 378);

    auto ta1_y_yyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 379);

    auto ta1_y_yyyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 380);

    auto ta1_y_yyyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 381);

    auto ta1_y_yyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 382);

    auto ta1_y_yyyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 383);

    auto ta1_y_yyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 384);

    auto ta1_y_yyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 385);

    auto ta1_y_yyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 386);

    auto ta1_y_yyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 387);

    auto ta1_y_yyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 388);

    auto ta1_y_yyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 389);

    auto ta1_y_yyzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 390);

    auto ta1_y_yyzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 391);

    auto ta1_y_yyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 392);

    auto ta1_y_yyzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 393);

    auto ta1_y_yyzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 394);

    auto ta1_y_yyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 395);

    auto ta1_y_yyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 396);

    auto ta1_y_yyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 397);

    auto ta1_y_yyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 398);

    auto ta1_y_yyzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 399);

    auto ta1_y_yzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 401);

    auto ta1_y_yzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 403);

    auto ta1_y_yzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 404);

    auto ta1_y_yzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 406);

    auto ta1_y_yzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 407);

    auto ta1_y_yzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 408);

    auto ta1_y_zzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 410);

    auto ta1_y_zzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 411);

    auto ta1_y_zzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 412);

    auto ta1_y_zzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 413);

    auto ta1_y_zzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 414);

    auto ta1_y_zzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 415);

    auto ta1_y_zzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 416);

    auto ta1_y_zzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 417);

    auto ta1_y_zzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 418);

    auto ta1_y_zzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 419);

    auto ta1_z_xxxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 420);

    auto ta1_z_xxxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 421);

    auto ta1_z_xxxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 422);

    auto ta1_z_xxxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 423);

    auto ta1_z_xxxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 424);

    auto ta1_z_xxxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 425);

    auto ta1_z_xxxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 426);

    auto ta1_z_xxxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 427);

    auto ta1_z_xxxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 428);

    auto ta1_z_xxxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 429);

    auto ta1_z_xxxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 442);

    auto ta1_z_xxxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 444);

    auto ta1_z_xxxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 445);

    auto ta1_z_xxxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 451);

    auto ta1_z_xxxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 453);

    auto ta1_z_xxxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 454);

    auto ta1_z_xxxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 456);

    auto ta1_z_xxxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 457);

    auto ta1_z_xxxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 458);

    auto ta1_z_xxxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 470);

    auto ta1_z_xxxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 471);

    auto ta1_z_xxxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 472);

    auto ta1_z_xxxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 473);

    auto ta1_z_xxxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 474);

    auto ta1_z_xxxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 475);

    auto ta1_z_xxxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 477);

    auto ta1_z_xxxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 478);

    auto ta1_z_xxxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 479);

    auto ta1_z_xxyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 481);

    auto ta1_z_xxyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 483);

    auto ta1_z_xxyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 484);

    auto ta1_z_xxyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 486);

    auto ta1_z_xxyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 487);

    auto ta1_z_xxyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 488);

    auto ta1_z_xxzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 510);

    auto ta1_z_xxzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 511);

    auto ta1_z_xxzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 512);

    auto ta1_z_xxzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 513);

    auto ta1_z_xxzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 514);

    auto ta1_z_xxzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 515);

    auto ta1_z_xxzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 517);

    auto ta1_z_xxzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 518);

    auto ta1_z_xxzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 519);

    auto ta1_z_xyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 521);

    auto ta1_z_xyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 523);

    auto ta1_z_xyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 524);

    auto ta1_z_xyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 526);

    auto ta1_z_xyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 527);

    auto ta1_z_xyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 528);

    auto ta1_z_xyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 544);

    auto ta1_z_xyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 547);

    auto ta1_z_xyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 548);

    auto ta1_z_xzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 562);

    auto ta1_z_xzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 564);

    auto ta1_z_xzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 565);

    auto ta1_z_xzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 567);

    auto ta1_z_xzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 568);

    auto ta1_z_xzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 569);

    auto ta1_z_yyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 570);

    auto ta1_z_yyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 571);

    auto ta1_z_yyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 572);

    auto ta1_z_yyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 573);

    auto ta1_z_yyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 574);

    auto ta1_z_yyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 575);

    auto ta1_z_yyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 576);

    auto ta1_z_yyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 577);

    auto ta1_z_yyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 578);

    auto ta1_z_yyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 579);

    auto ta1_z_yyyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 582);

    auto ta1_z_yyyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 584);

    auto ta1_z_yyyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 585);

    auto ta1_z_yyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 587);

    auto ta1_z_yyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 588);

    auto ta1_z_yyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 589);

    auto ta1_z_yyyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 590);

    auto ta1_z_yyyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 591);

    auto ta1_z_yyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 592);

    auto ta1_z_yyyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 593);

    auto ta1_z_yyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 594);

    auto ta1_z_yyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 595);

    auto ta1_z_yyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 596);

    auto ta1_z_yyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 597);

    auto ta1_z_yyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 598);

    auto ta1_z_yyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 599);

    auto ta1_z_yyzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 600);

    auto ta1_z_yyzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 601);

    auto ta1_z_yyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 602);

    auto ta1_z_yyzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 603);

    auto ta1_z_yyzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 604);

    auto ta1_z_yyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 605);

    auto ta1_z_yyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 606);

    auto ta1_z_yyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 607);

    auto ta1_z_yyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 608);

    auto ta1_z_yyzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 609);

    auto ta1_z_yzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 611);

    auto ta1_z_yzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 612);

    auto ta1_z_yzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 613);

    auto ta1_z_yzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 614);

    auto ta1_z_yzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 615);

    auto ta1_z_yzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 616);

    auto ta1_z_yzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 617);

    auto ta1_z_yzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 618);

    auto ta1_z_yzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 619);

    auto ta1_z_zzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 620);

    auto ta1_z_zzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 621);

    auto ta1_z_zzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 622);

    auto ta1_z_zzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 623);

    auto ta1_z_zzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 624);

    auto ta1_z_zzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 625);

    auto ta1_z_zzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 626);

    auto ta1_z_zzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 627);

    auto ta1_z_zzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 628);

    auto ta1_z_zzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 629);

    // Set up components of auxiliary buffer : HG

    auto ta_xxxxx_xxxx_1 = pbuffer.data(idx_npot_1_hg);

    auto ta_xxxxx_xxxy_1 = pbuffer.data(idx_npot_1_hg + 1);

    auto ta_xxxxx_xxxz_1 = pbuffer.data(idx_npot_1_hg + 2);

    auto ta_xxxxx_xxyy_1 = pbuffer.data(idx_npot_1_hg + 3);

    auto ta_xxxxx_xxyz_1 = pbuffer.data(idx_npot_1_hg + 4);

    auto ta_xxxxx_xxzz_1 = pbuffer.data(idx_npot_1_hg + 5);

    auto ta_xxxxx_xyyy_1 = pbuffer.data(idx_npot_1_hg + 6);

    auto ta_xxxxx_xyyz_1 = pbuffer.data(idx_npot_1_hg + 7);

    auto ta_xxxxx_xyzz_1 = pbuffer.data(idx_npot_1_hg + 8);

    auto ta_xxxxx_xzzz_1 = pbuffer.data(idx_npot_1_hg + 9);

    auto ta_xxxxx_yyyy_1 = pbuffer.data(idx_npot_1_hg + 10);

    auto ta_xxxxx_yyyz_1 = pbuffer.data(idx_npot_1_hg + 11);

    auto ta_xxxxx_yyzz_1 = pbuffer.data(idx_npot_1_hg + 12);

    auto ta_xxxxx_yzzz_1 = pbuffer.data(idx_npot_1_hg + 13);

    auto ta_xxxxx_zzzz_1 = pbuffer.data(idx_npot_1_hg + 14);

    auto ta_xxxxy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 15);

    auto ta_xxxxy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 16);

    auto ta_xxxxy_xxxz_1 = pbuffer.data(idx_npot_1_hg + 17);

    auto ta_xxxxy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 18);

    auto ta_xxxxy_xxzz_1 = pbuffer.data(idx_npot_1_hg + 20);

    auto ta_xxxxy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 21);

    auto ta_xxxxy_xzzz_1 = pbuffer.data(idx_npot_1_hg + 24);

    auto ta_xxxxy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 25);

    auto ta_xxxxz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 30);

    auto ta_xxxxz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 31);

    auto ta_xxxxz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 32);

    auto ta_xxxxz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 33);

    auto ta_xxxxz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 35);

    auto ta_xxxxz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 36);

    auto ta_xxxxz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 39);

    auto ta_xxxxz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 44);

    auto ta_xxxyy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 45);

    auto ta_xxxyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 46);

    auto ta_xxxyy_xxxz_1 = pbuffer.data(idx_npot_1_hg + 47);

    auto ta_xxxyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 48);

    auto ta_xxxyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 49);

    auto ta_xxxyy_xxzz_1 = pbuffer.data(idx_npot_1_hg + 50);

    auto ta_xxxyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 51);

    auto ta_xxxyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 52);

    auto ta_xxxyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 53);

    auto ta_xxxyy_xzzz_1 = pbuffer.data(idx_npot_1_hg + 54);

    auto ta_xxxyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 55);

    auto ta_xxxyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 56);

    auto ta_xxxyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 57);

    auto ta_xxxyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 58);

    auto ta_xxxzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 75);

    auto ta_xxxzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 76);

    auto ta_xxxzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 77);

    auto ta_xxxzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 78);

    auto ta_xxxzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 79);

    auto ta_xxxzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 80);

    auto ta_xxxzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 81);

    auto ta_xxxzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 82);

    auto ta_xxxzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 83);

    auto ta_xxxzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 84);

    auto ta_xxxzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 86);

    auto ta_xxxzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 87);

    auto ta_xxxzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 88);

    auto ta_xxxzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 89);

    auto ta_xxyyy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 90);

    auto ta_xxyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 91);

    auto ta_xxyyy_xxxz_1 = pbuffer.data(idx_npot_1_hg + 92);

    auto ta_xxyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 93);

    auto ta_xxyyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 94);

    auto ta_xxyyy_xxzz_1 = pbuffer.data(idx_npot_1_hg + 95);

    auto ta_xxyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 96);

    auto ta_xxyyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 97);

    auto ta_xxyyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 98);

    auto ta_xxyyy_xzzz_1 = pbuffer.data(idx_npot_1_hg + 99);

    auto ta_xxyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 100);

    auto ta_xxyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 101);

    auto ta_xxyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 102);

    auto ta_xxyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 103);

    auto ta_xxyyz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 106);

    auto ta_xxyyz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 108);

    auto ta_xxyyz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 111);

    auto ta_xxyzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 122);

    auto ta_xxyzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 125);

    auto ta_xxyzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 129);

    auto ta_xxzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 135);

    auto ta_xxzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 136);

    auto ta_xxzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 137);

    auto ta_xxzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 138);

    auto ta_xxzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 139);

    auto ta_xxzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 140);

    auto ta_xxzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 141);

    auto ta_xxzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 142);

    auto ta_xxzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 143);

    auto ta_xxzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 144);

    auto ta_xxzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 146);

    auto ta_xxzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 147);

    auto ta_xxzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 148);

    auto ta_xxzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 149);

    auto ta_xyyyy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 150);

    auto ta_xyyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 151);

    auto ta_xyyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 153);

    auto ta_xyyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 156);

    auto ta_xyyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 160);

    auto ta_xyyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 161);

    auto ta_xyyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 162);

    auto ta_xyyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 163);

    auto ta_xyyzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 191);

    auto ta_xyyzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 192);

    auto ta_xyyzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 193);

    auto ta_xzzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 210);

    auto ta_xzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 212);

    auto ta_xzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 215);

    auto ta_xzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 219);

    auto ta_xzzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 221);

    auto ta_xzzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 222);

    auto ta_xzzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 223);

    auto ta_xzzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 224);

    auto ta_yyyyy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 225);

    auto ta_yyyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 226);

    auto ta_yyyyy_xxxz_1 = pbuffer.data(idx_npot_1_hg + 227);

    auto ta_yyyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 228);

    auto ta_yyyyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 229);

    auto ta_yyyyy_xxzz_1 = pbuffer.data(idx_npot_1_hg + 230);

    auto ta_yyyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 231);

    auto ta_yyyyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 232);

    auto ta_yyyyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 233);

    auto ta_yyyyy_xzzz_1 = pbuffer.data(idx_npot_1_hg + 234);

    auto ta_yyyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 235);

    auto ta_yyyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 236);

    auto ta_yyyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 237);

    auto ta_yyyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 238);

    auto ta_yyyyy_zzzz_1 = pbuffer.data(idx_npot_1_hg + 239);

    auto ta_yyyyz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 241);

    auto ta_yyyyz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 243);

    auto ta_yyyyz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 246);

    auto ta_yyyyz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 250);

    auto ta_yyyyz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 251);

    auto ta_yyyyz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 252);

    auto ta_yyyyz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 253);

    auto ta_yyyyz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 254);

    auto ta_yyyzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 256);

    auto ta_yyyzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 257);

    auto ta_yyyzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 258);

    auto ta_yyyzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 259);

    auto ta_yyyzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 260);

    auto ta_yyyzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 261);

    auto ta_yyyzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 262);

    auto ta_yyyzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 263);

    auto ta_yyyzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 264);

    auto ta_yyyzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 265);

    auto ta_yyyzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 266);

    auto ta_yyyzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 267);

    auto ta_yyyzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 268);

    auto ta_yyyzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 269);

    auto ta_yyzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 271);

    auto ta_yyzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 272);

    auto ta_yyzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 273);

    auto ta_yyzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 274);

    auto ta_yyzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 275);

    auto ta_yyzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 276);

    auto ta_yyzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 277);

    auto ta_yyzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 278);

    auto ta_yyzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 279);

    auto ta_yyzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 280);

    auto ta_yyzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 281);

    auto ta_yyzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 282);

    auto ta_yyzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 283);

    auto ta_yyzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 284);

    auto ta_yzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 287);

    auto ta_yzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 290);

    auto ta_yzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 294);

    auto ta_yzzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 295);

    auto ta_yzzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 296);

    auto ta_yzzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 297);

    auto ta_yzzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 298);

    auto ta_yzzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 299);

    auto ta_zzzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 300);

    auto ta_zzzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 301);

    auto ta_zzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 302);

    auto ta_zzzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 303);

    auto ta_zzzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 304);

    auto ta_zzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 305);

    auto ta_zzzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 306);

    auto ta_zzzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 307);

    auto ta_zzzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 308);

    auto ta_zzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 309);

    auto ta_zzzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 310);

    auto ta_zzzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 311);

    auto ta_zzzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 312);

    auto ta_zzzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 313);

    auto ta_zzzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 314);

    // Set up components of auxiliary buffer : HG

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

    auto ta1_x_xxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 29);

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

    auto ta1_x_xxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 62);

    auto ta1_x_xxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 65);

    auto ta1_x_xxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 69);

    auto ta1_x_xxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 74);

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

    auto ta1_x_xxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 106);

    auto ta1_x_xxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 107);

    auto ta1_x_xxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 108);

    auto ta1_x_xxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 110);

    auto ta1_x_xxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 111);

    auto ta1_x_xxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 114);

    auto ta1_x_xxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 115);

    auto ta1_x_xxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 119);

    auto ta1_x_xxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 120);

    auto ta1_x_xxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 122);

    auto ta1_x_xxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 124);

    auto ta1_x_xxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 125);

    auto ta1_x_xxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 127);

    auto ta1_x_xxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 128);

    auto ta1_x_xxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 129);

    auto ta1_x_xxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 134);

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

    auto ta1_x_xyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 166);

    auto ta1_x_xyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 167);

    auto ta1_x_xyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 168);

    auto ta1_x_xyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 170);

    auto ta1_x_xyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 171);

    auto ta1_x_xyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 174);

    auto ta1_x_xyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 180);

    auto ta1_x_xyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 181);

    auto ta1_x_xyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 182);

    auto ta1_x_xyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 183);

    auto ta1_x_xyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 185);

    auto ta1_x_xyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 186);

    auto ta1_x_xyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 189);

    auto ta1_x_xyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 191);

    auto ta1_x_xyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 192);

    auto ta1_x_xyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 193);

    auto ta1_x_xyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 195);

    auto ta1_x_xyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 197);

    auto ta1_x_xyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 200);

    auto ta1_x_xyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 204);

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

    auto ta1_x_xzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 221);

    auto ta1_x_xzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 222);

    auto ta1_x_xzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 223);

    auto ta1_x_xzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 224);

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

    auto ta1_x_yyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 241);

    auto ta1_x_yyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 242);

    auto ta1_x_yyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 243);

    auto ta1_x_yyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 245);

    auto ta1_x_yyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 246);

    auto ta1_x_yyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 249);

    auto ta1_x_yyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 250);

    auto ta1_x_yyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 251);

    auto ta1_x_yyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 252);

    auto ta1_x_yyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 253);

    auto ta1_x_yyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 254);

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

    auto ta1_x_yzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 285);

    auto ta1_x_yzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 287);

    auto ta1_x_yzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 289);

    auto ta1_x_yzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 290);

    auto ta1_x_yzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 292);

    auto ta1_x_yzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 293);

    auto ta1_x_yzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 294);

    auto ta1_x_yzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 295);

    auto ta1_x_yzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 296);

    auto ta1_x_yzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 297);

    auto ta1_x_yzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 298);

    auto ta1_x_yzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 299);

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

    auto ta1_y_xxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 345);

    auto ta1_y_xxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 346);

    auto ta1_y_xxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 347);

    auto ta1_y_xxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 348);

    auto ta1_y_xxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 350);

    auto ta1_y_xxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 351);

    auto ta1_y_xxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 354);

    auto ta1_y_xxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 356);

    auto ta1_y_xxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 357);

    auto ta1_y_xxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 358);

    auto ta1_y_xxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 359);

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

    auto ta1_y_xxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 376);

    auto ta1_y_xxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 378);

    auto ta1_y_xxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 381);

    auto ta1_y_xxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 386);

    auto ta1_y_xxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 387);

    auto ta1_y_xxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 388);

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

    auto ta1_y_xxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 420);

    auto ta1_y_xxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 421);

    auto ta1_y_xxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 423);

    auto ta1_y_xxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 426);

    auto ta1_y_xxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 431);

    auto ta1_y_xxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 432);

    auto ta1_y_xxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 433);

    auto ta1_y_xxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 434);

    auto ta1_y_xxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 436);

    auto ta1_y_xxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 437);

    auto ta1_y_xxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 438);

    auto ta1_y_xxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 440);

    auto ta1_y_xxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 441);

    auto ta1_y_xxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 444);

    auto ta1_y_xxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 445);

    auto ta1_y_xxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 446);

    auto ta1_y_xxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 447);

    auto ta1_y_xxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 448);

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

    auto ta1_y_xyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 465);

    auto ta1_y_xyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 466);

    auto ta1_y_xyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 468);

    auto ta1_y_xyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 469);

    auto ta1_y_xyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 471);

    auto ta1_y_xyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 472);

    auto ta1_y_xyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 473);

    auto ta1_y_xyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 475);

    auto ta1_y_xyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 476);

    auto ta1_y_xyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 477);

    auto ta1_y_xyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 478);

    auto ta1_y_xyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 479);

    auto ta1_y_xyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 491);

    auto ta1_y_xyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 492);

    auto ta1_y_xyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 493);

    auto ta1_y_xyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 494);

    auto ta1_y_xyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 499);

    auto ta1_y_xyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 502);

    auto ta1_y_xyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 503);

    auto ta1_y_xyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 505);

    auto ta1_y_xyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 506);

    auto ta1_y_xyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 507);

    auto ta1_y_xyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 508);

    auto ta1_y_xyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 509);

    auto ta1_y_xyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 520);

    auto ta1_y_xyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 521);

    auto ta1_y_xyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 522);

    auto ta1_y_xyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 523);

    auto ta1_y_xzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 525);

    auto ta1_y_xzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 527);

    auto ta1_y_xzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 529);

    auto ta1_y_xzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 530);

    auto ta1_y_xzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 532);

    auto ta1_y_xzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 533);

    auto ta1_y_xzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 534);

    auto ta1_y_xzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 535);

    auto ta1_y_xzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 536);

    auto ta1_y_xzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 537);

    auto ta1_y_xzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 538);

    auto ta1_y_xzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 539);

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

    auto ta1_z_xxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 645);

    auto ta1_z_xxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 646);

    auto ta1_z_xxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 647);

    auto ta1_z_xxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 648);

    auto ta1_z_xxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 650);

    auto ta1_z_xxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 651);

    auto ta1_z_xxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 654);

    auto ta1_z_xxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 655);

    auto ta1_z_xxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 656);

    auto ta1_z_xxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 657);

    auto ta1_z_xxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 658);

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

    auto ta1_z_xxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 671);

    auto ta1_z_xxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 672);

    auto ta1_z_xxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 673);

    auto ta1_z_xxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 674);

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

    auto ta1_z_xxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 692);

    auto ta1_z_xxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 695);

    auto ta1_z_xxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 699);

    auto ta1_z_xxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 701);

    auto ta1_z_xxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 702);

    auto ta1_z_xxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 703);

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

    auto ta1_z_xxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 736);

    auto ta1_z_xxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 737);

    auto ta1_z_xxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 738);

    auto ta1_z_xxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 740);

    auto ta1_z_xxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 741);

    auto ta1_z_xxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 744);

    auto ta1_z_xxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 746);

    auto ta1_z_xxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 747);

    auto ta1_z_xxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 748);

    auto ta1_z_xxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 749);

    auto ta1_z_xxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 750);

    auto ta1_z_xxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 752);

    auto ta1_z_xxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 755);

    auto ta1_z_xxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 759);

    auto ta1_z_xxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 760);

    auto ta1_z_xxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 761);

    auto ta1_z_xxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 762);

    auto ta1_z_xxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 763);

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

    auto ta1_z_xyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 780);

    auto ta1_z_xyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 781);

    auto ta1_z_xyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 783);

    auto ta1_z_xyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 784);

    auto ta1_z_xyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 786);

    auto ta1_z_xyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 787);

    auto ta1_z_xyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 788);

    auto ta1_z_xyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 790);

    auto ta1_z_xyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 791);

    auto ta1_z_xyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 792);

    auto ta1_z_xyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 793);

    auto ta1_z_xyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 794);

    auto ta1_z_xyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 806);

    auto ta1_z_xyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 807);

    auto ta1_z_xyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 808);

    auto ta1_z_xyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 809);

    auto ta1_z_xyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 814);

    auto ta1_z_xyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 817);

    auto ta1_z_xyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 818);

    auto ta1_z_xyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 820);

    auto ta1_z_xyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 821);

    auto ta1_z_xyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 822);

    auto ta1_z_xyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 823);

    auto ta1_z_xyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 824);

    auto ta1_z_xyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 835);

    auto ta1_z_xyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 836);

    auto ta1_z_xyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 837);

    auto ta1_z_xyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 838);

    auto ta1_z_xzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_hg + 840);

    auto ta1_z_xzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 842);

    auto ta1_z_xzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 844);

    auto ta1_z_xzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 845);

    auto ta1_z_xzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 847);

    auto ta1_z_xzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 848);

    auto ta1_z_xzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 849);

    auto ta1_z_xzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_hg + 850);

    auto ta1_z_xzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 851);

    auto ta1_z_xzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 852);

    auto ta1_z_xzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 853);

    auto ta1_z_xzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_hg + 854);

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

    // Set up components of auxiliary buffer : HG

    auto ta1_x_xxxxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg);

    auto ta1_x_xxxxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 1);

    auto ta1_x_xxxxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 2);

    auto ta1_x_xxxxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 3);

    auto ta1_x_xxxxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 4);

    auto ta1_x_xxxxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 5);

    auto ta1_x_xxxxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 6);

    auto ta1_x_xxxxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 7);

    auto ta1_x_xxxxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 8);

    auto ta1_x_xxxxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 9);

    auto ta1_x_xxxxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 10);

    auto ta1_x_xxxxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 11);

    auto ta1_x_xxxxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 12);

    auto ta1_x_xxxxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 13);

    auto ta1_x_xxxxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 14);

    auto ta1_x_xxxxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 15);

    auto ta1_x_xxxxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 16);

    auto ta1_x_xxxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 17);

    auto ta1_x_xxxxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 18);

    auto ta1_x_xxxxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 19);

    auto ta1_x_xxxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 20);

    auto ta1_x_xxxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 21);

    auto ta1_x_xxxxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 22);

    auto ta1_x_xxxxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 23);

    auto ta1_x_xxxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 24);

    auto ta1_x_xxxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 25);

    auto ta1_x_xxxxy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 29);

    auto ta1_x_xxxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 30);

    auto ta1_x_xxxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 31);

    auto ta1_x_xxxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 32);

    auto ta1_x_xxxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 33);

    auto ta1_x_xxxxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 34);

    auto ta1_x_xxxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 35);

    auto ta1_x_xxxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 36);

    auto ta1_x_xxxxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 37);

    auto ta1_x_xxxxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 38);

    auto ta1_x_xxxxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 39);

    auto ta1_x_xxxxz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 40);

    auto ta1_x_xxxxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 41);

    auto ta1_x_xxxxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 42);

    auto ta1_x_xxxxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 43);

    auto ta1_x_xxxxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 44);

    auto ta1_x_xxxyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 45);

    auto ta1_x_xxxyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 46);

    auto ta1_x_xxxyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 47);

    auto ta1_x_xxxyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 48);

    auto ta1_x_xxxyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 49);

    auto ta1_x_xxxyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 50);

    auto ta1_x_xxxyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 51);

    auto ta1_x_xxxyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 52);

    auto ta1_x_xxxyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 53);

    auto ta1_x_xxxyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 54);

    auto ta1_x_xxxyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 55);

    auto ta1_x_xxxyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 56);

    auto ta1_x_xxxyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 57);

    auto ta1_x_xxxyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 58);

    auto ta1_x_xxxyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 59);

    auto ta1_x_xxxyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 62);

    auto ta1_x_xxxyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 65);

    auto ta1_x_xxxyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 69);

    auto ta1_x_xxxyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 74);

    auto ta1_x_xxxzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 75);

    auto ta1_x_xxxzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 76);

    auto ta1_x_xxxzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 77);

    auto ta1_x_xxxzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 78);

    auto ta1_x_xxxzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 79);

    auto ta1_x_xxxzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 80);

    auto ta1_x_xxxzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 81);

    auto ta1_x_xxxzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 82);

    auto ta1_x_xxxzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 83);

    auto ta1_x_xxxzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 84);

    auto ta1_x_xxxzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 85);

    auto ta1_x_xxxzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 86);

    auto ta1_x_xxxzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 87);

    auto ta1_x_xxxzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 88);

    auto ta1_x_xxxzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 89);

    auto ta1_x_xxyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 90);

    auto ta1_x_xxyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 91);

    auto ta1_x_xxyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 92);

    auto ta1_x_xxyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 93);

    auto ta1_x_xxyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 94);

    auto ta1_x_xxyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 95);

    auto ta1_x_xxyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 96);

    auto ta1_x_xxyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 97);

    auto ta1_x_xxyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 98);

    auto ta1_x_xxyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 99);

    auto ta1_x_xxyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 100);

    auto ta1_x_xxyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 101);

    auto ta1_x_xxyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 102);

    auto ta1_x_xxyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 103);

    auto ta1_x_xxyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 104);

    auto ta1_x_xxyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 106);

    auto ta1_x_xxyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 107);

    auto ta1_x_xxyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 108);

    auto ta1_x_xxyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 110);

    auto ta1_x_xxyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 111);

    auto ta1_x_xxyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 114);

    auto ta1_x_xxyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 115);

    auto ta1_x_xxyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 119);

    auto ta1_x_xxyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 120);

    auto ta1_x_xxyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 122);

    auto ta1_x_xxyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 124);

    auto ta1_x_xxyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 125);

    auto ta1_x_xxyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 127);

    auto ta1_x_xxyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 128);

    auto ta1_x_xxyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 129);

    auto ta1_x_xxyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 134);

    auto ta1_x_xxzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 135);

    auto ta1_x_xxzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 136);

    auto ta1_x_xxzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 137);

    auto ta1_x_xxzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 138);

    auto ta1_x_xxzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 139);

    auto ta1_x_xxzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 140);

    auto ta1_x_xxzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 141);

    auto ta1_x_xxzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 142);

    auto ta1_x_xxzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 143);

    auto ta1_x_xxzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 144);

    auto ta1_x_xxzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 145);

    auto ta1_x_xxzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 146);

    auto ta1_x_xxzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 147);

    auto ta1_x_xxzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 148);

    auto ta1_x_xxzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 149);

    auto ta1_x_xyyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 150);

    auto ta1_x_xyyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 151);

    auto ta1_x_xyyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 152);

    auto ta1_x_xyyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 153);

    auto ta1_x_xyyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 154);

    auto ta1_x_xyyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 155);

    auto ta1_x_xyyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 156);

    auto ta1_x_xyyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 157);

    auto ta1_x_xyyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 158);

    auto ta1_x_xyyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 159);

    auto ta1_x_xyyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 160);

    auto ta1_x_xyyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 161);

    auto ta1_x_xyyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 162);

    auto ta1_x_xyyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 163);

    auto ta1_x_xyyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 166);

    auto ta1_x_xyyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 167);

    auto ta1_x_xyyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 168);

    auto ta1_x_xyyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 170);

    auto ta1_x_xyyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 171);

    auto ta1_x_xyyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 174);

    auto ta1_x_xyyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 180);

    auto ta1_x_xyyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 181);

    auto ta1_x_xyyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 182);

    auto ta1_x_xyyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 183);

    auto ta1_x_xyyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 185);

    auto ta1_x_xyyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 186);

    auto ta1_x_xyyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 189);

    auto ta1_x_xyyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 191);

    auto ta1_x_xyyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 192);

    auto ta1_x_xyyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 193);

    auto ta1_x_xyzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 195);

    auto ta1_x_xyzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 197);

    auto ta1_x_xyzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 200);

    auto ta1_x_xyzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 204);

    auto ta1_x_xzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 210);

    auto ta1_x_xzzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 211);

    auto ta1_x_xzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 212);

    auto ta1_x_xzzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 213);

    auto ta1_x_xzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 214);

    auto ta1_x_xzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 215);

    auto ta1_x_xzzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 216);

    auto ta1_x_xzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 217);

    auto ta1_x_xzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 218);

    auto ta1_x_xzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 219);

    auto ta1_x_xzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 221);

    auto ta1_x_xzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 222);

    auto ta1_x_xzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 223);

    auto ta1_x_xzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 224);

    auto ta1_x_yyyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 225);

    auto ta1_x_yyyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 226);

    auto ta1_x_yyyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 227);

    auto ta1_x_yyyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 228);

    auto ta1_x_yyyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 229);

    auto ta1_x_yyyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 230);

    auto ta1_x_yyyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 231);

    auto ta1_x_yyyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 232);

    auto ta1_x_yyyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 233);

    auto ta1_x_yyyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 234);

    auto ta1_x_yyyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 235);

    auto ta1_x_yyyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 236);

    auto ta1_x_yyyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 237);

    auto ta1_x_yyyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 238);

    auto ta1_x_yyyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 239);

    auto ta1_x_yyyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 241);

    auto ta1_x_yyyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 242);

    auto ta1_x_yyyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 243);

    auto ta1_x_yyyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 245);

    auto ta1_x_yyyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 246);

    auto ta1_x_yyyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 249);

    auto ta1_x_yyyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 250);

    auto ta1_x_yyyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 251);

    auto ta1_x_yyyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 252);

    auto ta1_x_yyyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 253);

    auto ta1_x_yyyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 254);

    auto ta1_x_yyyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 255);

    auto ta1_x_yyyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 256);

    auto ta1_x_yyyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 257);

    auto ta1_x_yyyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 258);

    auto ta1_x_yyyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 259);

    auto ta1_x_yyyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 260);

    auto ta1_x_yyyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 261);

    auto ta1_x_yyyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 262);

    auto ta1_x_yyyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 263);

    auto ta1_x_yyyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 264);

    auto ta1_x_yyyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 265);

    auto ta1_x_yyyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 266);

    auto ta1_x_yyyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 267);

    auto ta1_x_yyyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 268);

    auto ta1_x_yyyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 269);

    auto ta1_x_yyzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 270);

    auto ta1_x_yyzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 271);

    auto ta1_x_yyzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 272);

    auto ta1_x_yyzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 273);

    auto ta1_x_yyzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 274);

    auto ta1_x_yyzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 275);

    auto ta1_x_yyzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 276);

    auto ta1_x_yyzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 277);

    auto ta1_x_yyzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 278);

    auto ta1_x_yyzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 279);

    auto ta1_x_yyzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 280);

    auto ta1_x_yyzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 281);

    auto ta1_x_yyzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 282);

    auto ta1_x_yyzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 283);

    auto ta1_x_yyzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 284);

    auto ta1_x_yzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 285);

    auto ta1_x_yzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 287);

    auto ta1_x_yzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 289);

    auto ta1_x_yzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 290);

    auto ta1_x_yzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 292);

    auto ta1_x_yzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 293);

    auto ta1_x_yzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 294);

    auto ta1_x_yzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 295);

    auto ta1_x_yzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 296);

    auto ta1_x_yzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 297);

    auto ta1_x_yzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 298);

    auto ta1_x_yzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 299);

    auto ta1_x_zzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 300);

    auto ta1_x_zzzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 301);

    auto ta1_x_zzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 302);

    auto ta1_x_zzzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 303);

    auto ta1_x_zzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 304);

    auto ta1_x_zzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 305);

    auto ta1_x_zzzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 306);

    auto ta1_x_zzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 307);

    auto ta1_x_zzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 308);

    auto ta1_x_zzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 309);

    auto ta1_x_zzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 310);

    auto ta1_x_zzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 311);

    auto ta1_x_zzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 312);

    auto ta1_x_zzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 313);

    auto ta1_x_zzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 314);

    auto ta1_y_xxxxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 315);

    auto ta1_y_xxxxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 316);

    auto ta1_y_xxxxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 317);

    auto ta1_y_xxxxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 318);

    auto ta1_y_xxxxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 319);

    auto ta1_y_xxxxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 320);

    auto ta1_y_xxxxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 321);

    auto ta1_y_xxxxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 322);

    auto ta1_y_xxxxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 323);

    auto ta1_y_xxxxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 324);

    auto ta1_y_xxxxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 325);

    auto ta1_y_xxxxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 326);

    auto ta1_y_xxxxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 327);

    auto ta1_y_xxxxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 328);

    auto ta1_y_xxxxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 329);

    auto ta1_y_xxxxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 330);

    auto ta1_y_xxxxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 331);

    auto ta1_y_xxxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 332);

    auto ta1_y_xxxxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 333);

    auto ta1_y_xxxxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 334);

    auto ta1_y_xxxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 335);

    auto ta1_y_xxxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 336);

    auto ta1_y_xxxxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 337);

    auto ta1_y_xxxxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 338);

    auto ta1_y_xxxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 339);

    auto ta1_y_xxxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 340);

    auto ta1_y_xxxxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 341);

    auto ta1_y_xxxxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 342);

    auto ta1_y_xxxxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 343);

    auto ta1_y_xxxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 345);

    auto ta1_y_xxxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 346);

    auto ta1_y_xxxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 347);

    auto ta1_y_xxxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 348);

    auto ta1_y_xxxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 350);

    auto ta1_y_xxxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 351);

    auto ta1_y_xxxxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 354);

    auto ta1_y_xxxxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 356);

    auto ta1_y_xxxxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 357);

    auto ta1_y_xxxxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 358);

    auto ta1_y_xxxxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 359);

    auto ta1_y_xxxyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 360);

    auto ta1_y_xxxyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 361);

    auto ta1_y_xxxyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 362);

    auto ta1_y_xxxyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 363);

    auto ta1_y_xxxyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 364);

    auto ta1_y_xxxyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 365);

    auto ta1_y_xxxyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 366);

    auto ta1_y_xxxyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 367);

    auto ta1_y_xxxyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 368);

    auto ta1_y_xxxyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 369);

    auto ta1_y_xxxyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 370);

    auto ta1_y_xxxyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 371);

    auto ta1_y_xxxyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 372);

    auto ta1_y_xxxyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 373);

    auto ta1_y_xxxyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 374);

    auto ta1_y_xxxyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 376);

    auto ta1_y_xxxyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 378);

    auto ta1_y_xxxyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 381);

    auto ta1_y_xxxyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 386);

    auto ta1_y_xxxyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 387);

    auto ta1_y_xxxyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 388);

    auto ta1_y_xxxzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 390);

    auto ta1_y_xxxzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 391);

    auto ta1_y_xxxzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 392);

    auto ta1_y_xxxzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 393);

    auto ta1_y_xxxzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 394);

    auto ta1_y_xxxzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 395);

    auto ta1_y_xxxzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 396);

    auto ta1_y_xxxzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 397);

    auto ta1_y_xxxzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 398);

    auto ta1_y_xxxzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 399);

    auto ta1_y_xxxzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 400);

    auto ta1_y_xxxzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 401);

    auto ta1_y_xxxzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 402);

    auto ta1_y_xxxzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 403);

    auto ta1_y_xxxzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 404);

    auto ta1_y_xxyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 405);

    auto ta1_y_xxyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 406);

    auto ta1_y_xxyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 407);

    auto ta1_y_xxyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 408);

    auto ta1_y_xxyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 409);

    auto ta1_y_xxyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 410);

    auto ta1_y_xxyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 411);

    auto ta1_y_xxyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 412);

    auto ta1_y_xxyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 413);

    auto ta1_y_xxyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 414);

    auto ta1_y_xxyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 415);

    auto ta1_y_xxyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 416);

    auto ta1_y_xxyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 417);

    auto ta1_y_xxyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 418);

    auto ta1_y_xxyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 419);

    auto ta1_y_xxyyz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 420);

    auto ta1_y_xxyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 421);

    auto ta1_y_xxyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 423);

    auto ta1_y_xxyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 426);

    auto ta1_y_xxyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 431);

    auto ta1_y_xxyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 432);

    auto ta1_y_xxyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 433);

    auto ta1_y_xxyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 434);

    auto ta1_y_xxyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 436);

    auto ta1_y_xxyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 437);

    auto ta1_y_xxyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 438);

    auto ta1_y_xxyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 440);

    auto ta1_y_xxyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 441);

    auto ta1_y_xxyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 444);

    auto ta1_y_xxyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 445);

    auto ta1_y_xxyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 446);

    auto ta1_y_xxyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 447);

    auto ta1_y_xxyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 448);

    auto ta1_y_xxzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 450);

    auto ta1_y_xxzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 451);

    auto ta1_y_xxzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 452);

    auto ta1_y_xxzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 453);

    auto ta1_y_xxzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 454);

    auto ta1_y_xxzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 455);

    auto ta1_y_xxzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 456);

    auto ta1_y_xxzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 457);

    auto ta1_y_xxzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 458);

    auto ta1_y_xxzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 459);

    auto ta1_y_xxzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 460);

    auto ta1_y_xxzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 461);

    auto ta1_y_xxzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 462);

    auto ta1_y_xxzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 463);

    auto ta1_y_xxzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 464);

    auto ta1_y_xyyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 465);

    auto ta1_y_xyyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 466);

    auto ta1_y_xyyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 468);

    auto ta1_y_xyyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 469);

    auto ta1_y_xyyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 471);

    auto ta1_y_xyyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 472);

    auto ta1_y_xyyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 473);

    auto ta1_y_xyyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 475);

    auto ta1_y_xyyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 476);

    auto ta1_y_xyyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 477);

    auto ta1_y_xyyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 478);

    auto ta1_y_xyyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 479);

    auto ta1_y_xyyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 491);

    auto ta1_y_xyyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 492);

    auto ta1_y_xyyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 493);

    auto ta1_y_xyyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 494);

    auto ta1_y_xyyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 499);

    auto ta1_y_xyyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 502);

    auto ta1_y_xyyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 503);

    auto ta1_y_xyyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 505);

    auto ta1_y_xyyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 506);

    auto ta1_y_xyyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 507);

    auto ta1_y_xyyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 508);

    auto ta1_y_xyyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 509);

    auto ta1_y_xyzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 520);

    auto ta1_y_xyzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 521);

    auto ta1_y_xyzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 522);

    auto ta1_y_xyzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 523);

    auto ta1_y_xzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 525);

    auto ta1_y_xzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 527);

    auto ta1_y_xzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 529);

    auto ta1_y_xzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 530);

    auto ta1_y_xzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 532);

    auto ta1_y_xzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 533);

    auto ta1_y_xzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 534);

    auto ta1_y_xzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 535);

    auto ta1_y_xzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 536);

    auto ta1_y_xzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 537);

    auto ta1_y_xzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 538);

    auto ta1_y_xzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 539);

    auto ta1_y_yyyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 540);

    auto ta1_y_yyyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 541);

    auto ta1_y_yyyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 542);

    auto ta1_y_yyyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 543);

    auto ta1_y_yyyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 544);

    auto ta1_y_yyyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 545);

    auto ta1_y_yyyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 546);

    auto ta1_y_yyyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 547);

    auto ta1_y_yyyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 548);

    auto ta1_y_yyyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 549);

    auto ta1_y_yyyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 550);

    auto ta1_y_yyyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 551);

    auto ta1_y_yyyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 552);

    auto ta1_y_yyyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 553);

    auto ta1_y_yyyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 554);

    auto ta1_y_yyyyz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 555);

    auto ta1_y_yyyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 556);

    auto ta1_y_yyyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 557);

    auto ta1_y_yyyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 558);

    auto ta1_y_yyyyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 559);

    auto ta1_y_yyyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 560);

    auto ta1_y_yyyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 561);

    auto ta1_y_yyyyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 562);

    auto ta1_y_yyyyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 563);

    auto ta1_y_yyyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 564);

    auto ta1_y_yyyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 565);

    auto ta1_y_yyyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 566);

    auto ta1_y_yyyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 567);

    auto ta1_y_yyyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 568);

    auto ta1_y_yyyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 569);

    auto ta1_y_yyyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 570);

    auto ta1_y_yyyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 571);

    auto ta1_y_yyyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 572);

    auto ta1_y_yyyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 573);

    auto ta1_y_yyyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 574);

    auto ta1_y_yyyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 575);

    auto ta1_y_yyyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 576);

    auto ta1_y_yyyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 577);

    auto ta1_y_yyyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 578);

    auto ta1_y_yyyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 579);

    auto ta1_y_yyyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 580);

    auto ta1_y_yyyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 581);

    auto ta1_y_yyyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 582);

    auto ta1_y_yyyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 583);

    auto ta1_y_yyyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 584);

    auto ta1_y_yyzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 585);

    auto ta1_y_yyzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 586);

    auto ta1_y_yyzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 587);

    auto ta1_y_yyzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 588);

    auto ta1_y_yyzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 589);

    auto ta1_y_yyzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 590);

    auto ta1_y_yyzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 591);

    auto ta1_y_yyzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 592);

    auto ta1_y_yyzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 593);

    auto ta1_y_yyzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 594);

    auto ta1_y_yyzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 595);

    auto ta1_y_yyzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 596);

    auto ta1_y_yyzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 597);

    auto ta1_y_yyzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 598);

    auto ta1_y_yyzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 599);

    auto ta1_y_yzzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 601);

    auto ta1_y_yzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 602);

    auto ta1_y_yzzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 603);

    auto ta1_y_yzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 604);

    auto ta1_y_yzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 605);

    auto ta1_y_yzzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 606);

    auto ta1_y_yzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 607);

    auto ta1_y_yzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 608);

    auto ta1_y_yzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 609);

    auto ta1_y_yzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 610);

    auto ta1_y_yzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 611);

    auto ta1_y_yzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 612);

    auto ta1_y_yzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 613);

    auto ta1_y_yzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 614);

    auto ta1_y_zzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 615);

    auto ta1_y_zzzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 616);

    auto ta1_y_zzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 617);

    auto ta1_y_zzzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 618);

    auto ta1_y_zzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 619);

    auto ta1_y_zzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 620);

    auto ta1_y_zzzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 621);

    auto ta1_y_zzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 622);

    auto ta1_y_zzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 623);

    auto ta1_y_zzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 624);

    auto ta1_y_zzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 625);

    auto ta1_y_zzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 626);

    auto ta1_y_zzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 627);

    auto ta1_y_zzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 628);

    auto ta1_y_zzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 629);

    auto ta1_z_xxxxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 630);

    auto ta1_z_xxxxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 631);

    auto ta1_z_xxxxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 632);

    auto ta1_z_xxxxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 633);

    auto ta1_z_xxxxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 634);

    auto ta1_z_xxxxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 635);

    auto ta1_z_xxxxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 636);

    auto ta1_z_xxxxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 637);

    auto ta1_z_xxxxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 638);

    auto ta1_z_xxxxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 639);

    auto ta1_z_xxxxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 640);

    auto ta1_z_xxxxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 641);

    auto ta1_z_xxxxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 642);

    auto ta1_z_xxxxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 643);

    auto ta1_z_xxxxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 644);

    auto ta1_z_xxxxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 645);

    auto ta1_z_xxxxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 646);

    auto ta1_z_xxxxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 647);

    auto ta1_z_xxxxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 648);

    auto ta1_z_xxxxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 650);

    auto ta1_z_xxxxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 651);

    auto ta1_z_xxxxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 654);

    auto ta1_z_xxxxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 655);

    auto ta1_z_xxxxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 656);

    auto ta1_z_xxxxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 657);

    auto ta1_z_xxxxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 658);

    auto ta1_z_xxxxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 660);

    auto ta1_z_xxxxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 661);

    auto ta1_z_xxxxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 662);

    auto ta1_z_xxxxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 663);

    auto ta1_z_xxxxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 664);

    auto ta1_z_xxxxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 665);

    auto ta1_z_xxxxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 666);

    auto ta1_z_xxxxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 667);

    auto ta1_z_xxxxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 668);

    auto ta1_z_xxxxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 669);

    auto ta1_z_xxxxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 671);

    auto ta1_z_xxxxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 672);

    auto ta1_z_xxxxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 673);

    auto ta1_z_xxxxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 674);

    auto ta1_z_xxxyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 675);

    auto ta1_z_xxxyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 676);

    auto ta1_z_xxxyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 677);

    auto ta1_z_xxxyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 678);

    auto ta1_z_xxxyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 679);

    auto ta1_z_xxxyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 680);

    auto ta1_z_xxxyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 681);

    auto ta1_z_xxxyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 682);

    auto ta1_z_xxxyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 683);

    auto ta1_z_xxxyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 684);

    auto ta1_z_xxxyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 685);

    auto ta1_z_xxxyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 686);

    auto ta1_z_xxxyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 687);

    auto ta1_z_xxxyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 688);

    auto ta1_z_xxxyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 689);

    auto ta1_z_xxxyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 692);

    auto ta1_z_xxxyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 695);

    auto ta1_z_xxxyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 699);

    auto ta1_z_xxxyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 701);

    auto ta1_z_xxxyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 702);

    auto ta1_z_xxxyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 703);

    auto ta1_z_xxxzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 705);

    auto ta1_z_xxxzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 706);

    auto ta1_z_xxxzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 707);

    auto ta1_z_xxxzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 708);

    auto ta1_z_xxxzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 709);

    auto ta1_z_xxxzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 710);

    auto ta1_z_xxxzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 711);

    auto ta1_z_xxxzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 712);

    auto ta1_z_xxxzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 713);

    auto ta1_z_xxxzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 714);

    auto ta1_z_xxxzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 715);

    auto ta1_z_xxxzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 716);

    auto ta1_z_xxxzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 717);

    auto ta1_z_xxxzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 718);

    auto ta1_z_xxxzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 719);

    auto ta1_z_xxyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 720);

    auto ta1_z_xxyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 721);

    auto ta1_z_xxyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 722);

    auto ta1_z_xxyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 723);

    auto ta1_z_xxyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 724);

    auto ta1_z_xxyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 725);

    auto ta1_z_xxyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 726);

    auto ta1_z_xxyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 727);

    auto ta1_z_xxyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 728);

    auto ta1_z_xxyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 729);

    auto ta1_z_xxyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 730);

    auto ta1_z_xxyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 731);

    auto ta1_z_xxyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 732);

    auto ta1_z_xxyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 733);

    auto ta1_z_xxyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 734);

    auto ta1_z_xxyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 736);

    auto ta1_z_xxyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 737);

    auto ta1_z_xxyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 738);

    auto ta1_z_xxyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 740);

    auto ta1_z_xxyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 741);

    auto ta1_z_xxyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 744);

    auto ta1_z_xxyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 746);

    auto ta1_z_xxyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 747);

    auto ta1_z_xxyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 748);

    auto ta1_z_xxyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 749);

    auto ta1_z_xxyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 750);

    auto ta1_z_xxyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 752);

    auto ta1_z_xxyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 755);

    auto ta1_z_xxyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 759);

    auto ta1_z_xxyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 760);

    auto ta1_z_xxyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 761);

    auto ta1_z_xxyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 762);

    auto ta1_z_xxyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 763);

    auto ta1_z_xxzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 765);

    auto ta1_z_xxzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 766);

    auto ta1_z_xxzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 767);

    auto ta1_z_xxzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 768);

    auto ta1_z_xxzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 769);

    auto ta1_z_xxzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 770);

    auto ta1_z_xxzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 771);

    auto ta1_z_xxzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 772);

    auto ta1_z_xxzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 773);

    auto ta1_z_xxzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 774);

    auto ta1_z_xxzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 775);

    auto ta1_z_xxzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 776);

    auto ta1_z_xxzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 777);

    auto ta1_z_xxzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 778);

    auto ta1_z_xxzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 779);

    auto ta1_z_xyyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 780);

    auto ta1_z_xyyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 781);

    auto ta1_z_xyyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 783);

    auto ta1_z_xyyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 784);

    auto ta1_z_xyyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 786);

    auto ta1_z_xyyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 787);

    auto ta1_z_xyyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 788);

    auto ta1_z_xyyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 790);

    auto ta1_z_xyyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 791);

    auto ta1_z_xyyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 792);

    auto ta1_z_xyyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 793);

    auto ta1_z_xyyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 794);

    auto ta1_z_xyyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 806);

    auto ta1_z_xyyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 807);

    auto ta1_z_xyyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 808);

    auto ta1_z_xyyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 809);

    auto ta1_z_xyyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 814);

    auto ta1_z_xyyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 817);

    auto ta1_z_xyyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 818);

    auto ta1_z_xyyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 820);

    auto ta1_z_xyyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 821);

    auto ta1_z_xyyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 822);

    auto ta1_z_xyyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 823);

    auto ta1_z_xyyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 824);

    auto ta1_z_xyzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 835);

    auto ta1_z_xyzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 836);

    auto ta1_z_xyzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 837);

    auto ta1_z_xyzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 838);

    auto ta1_z_xzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 840);

    auto ta1_z_xzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 842);

    auto ta1_z_xzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 844);

    auto ta1_z_xzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 845);

    auto ta1_z_xzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 847);

    auto ta1_z_xzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 848);

    auto ta1_z_xzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 849);

    auto ta1_z_xzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 850);

    auto ta1_z_xzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 851);

    auto ta1_z_xzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 852);

    auto ta1_z_xzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 853);

    auto ta1_z_xzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 854);

    auto ta1_z_yyyyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 855);

    auto ta1_z_yyyyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 856);

    auto ta1_z_yyyyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 857);

    auto ta1_z_yyyyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 858);

    auto ta1_z_yyyyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 859);

    auto ta1_z_yyyyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 860);

    auto ta1_z_yyyyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 861);

    auto ta1_z_yyyyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 862);

    auto ta1_z_yyyyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 863);

    auto ta1_z_yyyyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 864);

    auto ta1_z_yyyyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 865);

    auto ta1_z_yyyyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 866);

    auto ta1_z_yyyyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 867);

    auto ta1_z_yyyyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 868);

    auto ta1_z_yyyyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 869);

    auto ta1_z_yyyyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 871);

    auto ta1_z_yyyyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 872);

    auto ta1_z_yyyyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 873);

    auto ta1_z_yyyyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 874);

    auto ta1_z_yyyyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 875);

    auto ta1_z_yyyyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 876);

    auto ta1_z_yyyyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 877);

    auto ta1_z_yyyyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 878);

    auto ta1_z_yyyyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 879);

    auto ta1_z_yyyyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 880);

    auto ta1_z_yyyyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 881);

    auto ta1_z_yyyyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 882);

    auto ta1_z_yyyyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 883);

    auto ta1_z_yyyyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 884);

    auto ta1_z_yyyzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 885);

    auto ta1_z_yyyzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 886);

    auto ta1_z_yyyzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 887);

    auto ta1_z_yyyzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 888);

    auto ta1_z_yyyzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 889);

    auto ta1_z_yyyzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 890);

    auto ta1_z_yyyzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 891);

    auto ta1_z_yyyzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 892);

    auto ta1_z_yyyzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 893);

    auto ta1_z_yyyzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 894);

    auto ta1_z_yyyzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 895);

    auto ta1_z_yyyzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 896);

    auto ta1_z_yyyzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 897);

    auto ta1_z_yyyzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 898);

    auto ta1_z_yyyzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 899);

    auto ta1_z_yyzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 900);

    auto ta1_z_yyzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 901);

    auto ta1_z_yyzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 902);

    auto ta1_z_yyzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 903);

    auto ta1_z_yyzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 904);

    auto ta1_z_yyzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 905);

    auto ta1_z_yyzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 906);

    auto ta1_z_yyzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 907);

    auto ta1_z_yyzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 908);

    auto ta1_z_yyzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 909);

    auto ta1_z_yyzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 910);

    auto ta1_z_yyzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 911);

    auto ta1_z_yyzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 912);

    auto ta1_z_yyzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 913);

    auto ta1_z_yyzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 914);

    auto ta1_z_yzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 915);

    auto ta1_z_yzzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 916);

    auto ta1_z_yzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 917);

    auto ta1_z_yzzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 918);

    auto ta1_z_yzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 919);

    auto ta1_z_yzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 920);

    auto ta1_z_yzzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 921);

    auto ta1_z_yzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 922);

    auto ta1_z_yzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 923);

    auto ta1_z_yzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 924);

    auto ta1_z_yzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 925);

    auto ta1_z_yzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 926);

    auto ta1_z_yzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 927);

    auto ta1_z_yzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 928);

    auto ta1_z_yzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 929);

    auto ta1_z_zzzzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_hg + 930);

    auto ta1_z_zzzzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 931);

    auto ta1_z_zzzzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 932);

    auto ta1_z_zzzzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 933);

    auto ta1_z_zzzzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 934);

    auto ta1_z_zzzzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 935);

    auto ta1_z_zzzzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 936);

    auto ta1_z_zzzzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 937);

    auto ta1_z_zzzzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 938);

    auto ta1_z_zzzzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 939);

    auto ta1_z_zzzzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_hg + 940);

    auto ta1_z_zzzzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 941);

    auto ta1_z_zzzzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 942);

    auto ta1_z_zzzzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 943);

    auto ta1_z_zzzzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_hg + 944);

    // Set up 0-15 components of targeted buffer : IG

    auto ta1_x_xxxxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig);

    auto ta1_x_xxxxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1);

    auto ta1_x_xxxxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 2);

    auto ta1_x_xxxxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 3);

    auto ta1_x_xxxxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 4);

    auto ta1_x_xxxxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 5);

    auto ta1_x_xxxxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 6);

    auto ta1_x_xxxxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 7);

    auto ta1_x_xxxxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 8);

    auto ta1_x_xxxxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 9);

    auto ta1_x_xxxxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 10);

    auto ta1_x_xxxxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 11);

    auto ta1_x_xxxxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 12);

    auto ta1_x_xxxxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 13);

    auto ta1_x_xxxxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 14);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xxxx_xxxx_0, ta1_x_xxxx_xxxx_1, ta1_x_xxxx_xxxy_0, ta1_x_xxxx_xxxy_1, ta1_x_xxxx_xxxz_0, ta1_x_xxxx_xxxz_1, ta1_x_xxxx_xxyy_0, ta1_x_xxxx_xxyy_1, ta1_x_xxxx_xxyz_0, ta1_x_xxxx_xxyz_1, ta1_x_xxxx_xxzz_0, ta1_x_xxxx_xxzz_1, ta1_x_xxxx_xyyy_0, ta1_x_xxxx_xyyy_1, ta1_x_xxxx_xyyz_0, ta1_x_xxxx_xyyz_1, ta1_x_xxxx_xyzz_0, ta1_x_xxxx_xyzz_1, ta1_x_xxxx_xzzz_0, ta1_x_xxxx_xzzz_1, ta1_x_xxxx_yyyy_0, ta1_x_xxxx_yyyy_1, ta1_x_xxxx_yyyz_0, ta1_x_xxxx_yyyz_1, ta1_x_xxxx_yyzz_0, ta1_x_xxxx_yyzz_1, ta1_x_xxxx_yzzz_0, ta1_x_xxxx_yzzz_1, ta1_x_xxxx_zzzz_0, ta1_x_xxxx_zzzz_1, ta1_x_xxxxx_xxx_0, ta1_x_xxxxx_xxx_1, ta1_x_xxxxx_xxxx_0, ta1_x_xxxxx_xxxx_1, ta1_x_xxxxx_xxxy_0, ta1_x_xxxxx_xxxy_1, ta1_x_xxxxx_xxxz_0, ta1_x_xxxxx_xxxz_1, ta1_x_xxxxx_xxy_0, ta1_x_xxxxx_xxy_1, ta1_x_xxxxx_xxyy_0, ta1_x_xxxxx_xxyy_1, ta1_x_xxxxx_xxyz_0, ta1_x_xxxxx_xxyz_1, ta1_x_xxxxx_xxz_0, ta1_x_xxxxx_xxz_1, ta1_x_xxxxx_xxzz_0, ta1_x_xxxxx_xxzz_1, ta1_x_xxxxx_xyy_0, ta1_x_xxxxx_xyy_1, ta1_x_xxxxx_xyyy_0, ta1_x_xxxxx_xyyy_1, ta1_x_xxxxx_xyyz_0, ta1_x_xxxxx_xyyz_1, ta1_x_xxxxx_xyz_0, ta1_x_xxxxx_xyz_1, ta1_x_xxxxx_xyzz_0, ta1_x_xxxxx_xyzz_1, ta1_x_xxxxx_xzz_0, ta1_x_xxxxx_xzz_1, ta1_x_xxxxx_xzzz_0, ta1_x_xxxxx_xzzz_1, ta1_x_xxxxx_yyy_0, ta1_x_xxxxx_yyy_1, ta1_x_xxxxx_yyyy_0, ta1_x_xxxxx_yyyy_1, ta1_x_xxxxx_yyyz_0, ta1_x_xxxxx_yyyz_1, ta1_x_xxxxx_yyz_0, ta1_x_xxxxx_yyz_1, ta1_x_xxxxx_yyzz_0, ta1_x_xxxxx_yyzz_1, ta1_x_xxxxx_yzz_0, ta1_x_xxxxx_yzz_1, ta1_x_xxxxx_yzzz_0, ta1_x_xxxxx_yzzz_1, ta1_x_xxxxx_zzz_0, ta1_x_xxxxx_zzz_1, ta1_x_xxxxx_zzzz_0, ta1_x_xxxxx_zzzz_1, ta1_x_xxxxxx_xxxx_0, ta1_x_xxxxxx_xxxy_0, ta1_x_xxxxxx_xxxz_0, ta1_x_xxxxxx_xxyy_0, ta1_x_xxxxxx_xxyz_0, ta1_x_xxxxxx_xxzz_0, ta1_x_xxxxxx_xyyy_0, ta1_x_xxxxxx_xyyz_0, ta1_x_xxxxxx_xyzz_0, ta1_x_xxxxxx_xzzz_0, ta1_x_xxxxxx_yyyy_0, ta1_x_xxxxxx_yyyz_0, ta1_x_xxxxxx_yyzz_0, ta1_x_xxxxxx_yzzz_0, ta1_x_xxxxxx_zzzz_0, ta_xxxxx_xxxx_1, ta_xxxxx_xxxy_1, ta_xxxxx_xxxz_1, ta_xxxxx_xxyy_1, ta_xxxxx_xxyz_1, ta_xxxxx_xxzz_1, ta_xxxxx_xyyy_1, ta_xxxxx_xyyz_1, ta_xxxxx_xyzz_1, ta_xxxxx_xzzz_1, ta_xxxxx_yyyy_1, ta_xxxxx_yyyz_1, ta_xxxxx_yyzz_1, ta_xxxxx_yzzz_1, ta_xxxxx_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxx_xxxx_0[i] = 5.0 * ta1_x_xxxx_xxxx_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxxx_1[i] * fe_0 + 4.0 * ta1_x_xxxxx_xxx_0[i] * fe_0 - 4.0 * ta1_x_xxxxx_xxx_1[i] * fe_0 + ta_xxxxx_xxxx_1[i] + ta1_x_xxxxx_xxxx_0[i] * pa_x[i] - ta1_x_xxxxx_xxxx_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxxy_0[i] = 5.0 * ta1_x_xxxx_xxxy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxxy_1[i] * fe_0 + 3.0 * ta1_x_xxxxx_xxy_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_xxy_1[i] * fe_0 + ta_xxxxx_xxxy_1[i] + ta1_x_xxxxx_xxxy_0[i] * pa_x[i] - ta1_x_xxxxx_xxxy_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxxz_0[i] = 5.0 * ta1_x_xxxx_xxxz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxxz_1[i] * fe_0 + 3.0 * ta1_x_xxxxx_xxz_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_xxz_1[i] * fe_0 + ta_xxxxx_xxxz_1[i] + ta1_x_xxxxx_xxxz_0[i] * pa_x[i] - ta1_x_xxxxx_xxxz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxyy_0[i] = 5.0 * ta1_x_xxxx_xxyy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxxxx_xyy_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xyy_1[i] * fe_0 + ta_xxxxx_xxyy_1[i] + ta1_x_xxxxx_xxyy_0[i] * pa_x[i] - ta1_x_xxxxx_xxyy_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxyz_0[i] = 5.0 * ta1_x_xxxx_xxyz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxyz_1[i] * fe_0 + 2.0 * ta1_x_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xyz_1[i] * fe_0 + ta_xxxxx_xxyz_1[i] + ta1_x_xxxxx_xxyz_0[i] * pa_x[i] - ta1_x_xxxxx_xxyz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxzz_0[i] = 5.0 * ta1_x_xxxx_xxzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxxxx_xzz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xzz_1[i] * fe_0 + ta_xxxxx_xxzz_1[i] + ta1_x_xxxxx_xxzz_0[i] * pa_x[i] - ta1_x_xxxxx_xxzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xyyy_0[i] = 5.0 * ta1_x_xxxx_xyyy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xyyy_1[i] * fe_0 + ta1_x_xxxxx_yyy_0[i] * fe_0 - ta1_x_xxxxx_yyy_1[i] * fe_0 + ta_xxxxx_xyyy_1[i] + ta1_x_xxxxx_xyyy_0[i] * pa_x[i] - ta1_x_xxxxx_xyyy_1[i] * pc_x[i];

        ta1_x_xxxxxx_xyyz_0[i] = 5.0 * ta1_x_xxxx_xyyz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xyyz_1[i] * fe_0 + ta1_x_xxxxx_yyz_0[i] * fe_0 - ta1_x_xxxxx_yyz_1[i] * fe_0 + ta_xxxxx_xyyz_1[i] + ta1_x_xxxxx_xyyz_0[i] * pa_x[i] - ta1_x_xxxxx_xyyz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xyzz_0[i] = 5.0 * ta1_x_xxxx_xyzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xyzz_1[i] * fe_0 + ta1_x_xxxxx_yzz_0[i] * fe_0 - ta1_x_xxxxx_yzz_1[i] * fe_0 + ta_xxxxx_xyzz_1[i] + ta1_x_xxxxx_xyzz_0[i] * pa_x[i] - ta1_x_xxxxx_xyzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xzzz_0[i] = 5.0 * ta1_x_xxxx_xzzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xzzz_1[i] * fe_0 + ta1_x_xxxxx_zzz_0[i] * fe_0 - ta1_x_xxxxx_zzz_1[i] * fe_0 + ta_xxxxx_xzzz_1[i] + ta1_x_xxxxx_xzzz_0[i] * pa_x[i] - ta1_x_xxxxx_xzzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_yyyy_0[i] = 5.0 * ta1_x_xxxx_yyyy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yyyy_1[i] * fe_0 + ta_xxxxx_yyyy_1[i] + ta1_x_xxxxx_yyyy_0[i] * pa_x[i] - ta1_x_xxxxx_yyyy_1[i] * pc_x[i];

        ta1_x_xxxxxx_yyyz_0[i] = 5.0 * ta1_x_xxxx_yyyz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yyyz_1[i] * fe_0 + ta_xxxxx_yyyz_1[i] + ta1_x_xxxxx_yyyz_0[i] * pa_x[i] - ta1_x_xxxxx_yyyz_1[i] * pc_x[i];

        ta1_x_xxxxxx_yyzz_0[i] = 5.0 * ta1_x_xxxx_yyzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yyzz_1[i] * fe_0 + ta_xxxxx_yyzz_1[i] + ta1_x_xxxxx_yyzz_0[i] * pa_x[i] - ta1_x_xxxxx_yyzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_yzzz_0[i] = 5.0 * ta1_x_xxxx_yzzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yzzz_1[i] * fe_0 + ta_xxxxx_yzzz_1[i] + ta1_x_xxxxx_yzzz_0[i] * pa_x[i] - ta1_x_xxxxx_yzzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_zzzz_0[i] = 5.0 * ta1_x_xxxx_zzzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_zzzz_1[i] * fe_0 + ta_xxxxx_zzzz_1[i] + ta1_x_xxxxx_zzzz_0[i] * pa_x[i] - ta1_x_xxxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : IG

    auto ta1_x_xxxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 15);

    auto ta1_x_xxxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 16);

    auto ta1_x_xxxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 17);

    auto ta1_x_xxxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 18);

    auto ta1_x_xxxxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 19);

    auto ta1_x_xxxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 20);

    auto ta1_x_xxxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 21);

    auto ta1_x_xxxxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 22);

    auto ta1_x_xxxxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 23);

    auto ta1_x_xxxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 24);

    auto ta1_x_xxxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 25);

    auto ta1_x_xxxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 26);

    auto ta1_x_xxxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 27);

    auto ta1_x_xxxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 28);

    auto ta1_x_xxxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 29);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxxxx_xxx_0, ta1_x_xxxxx_xxx_1, ta1_x_xxxxx_xxxx_0, ta1_x_xxxxx_xxxx_1, ta1_x_xxxxx_xxxy_0, ta1_x_xxxxx_xxxy_1, ta1_x_xxxxx_xxxz_0, ta1_x_xxxxx_xxxz_1, ta1_x_xxxxx_xxy_0, ta1_x_xxxxx_xxy_1, ta1_x_xxxxx_xxyy_0, ta1_x_xxxxx_xxyy_1, ta1_x_xxxxx_xxyz_0, ta1_x_xxxxx_xxyz_1, ta1_x_xxxxx_xxz_0, ta1_x_xxxxx_xxz_1, ta1_x_xxxxx_xxzz_0, ta1_x_xxxxx_xxzz_1, ta1_x_xxxxx_xyy_0, ta1_x_xxxxx_xyy_1, ta1_x_xxxxx_xyyy_0, ta1_x_xxxxx_xyyy_1, ta1_x_xxxxx_xyyz_0, ta1_x_xxxxx_xyyz_1, ta1_x_xxxxx_xyz_0, ta1_x_xxxxx_xyz_1, ta1_x_xxxxx_xyzz_0, ta1_x_xxxxx_xyzz_1, ta1_x_xxxxx_xzz_0, ta1_x_xxxxx_xzz_1, ta1_x_xxxxx_xzzz_0, ta1_x_xxxxx_xzzz_1, ta1_x_xxxxx_yyy_0, ta1_x_xxxxx_yyy_1, ta1_x_xxxxx_yyyy_0, ta1_x_xxxxx_yyyy_1, ta1_x_xxxxx_yyyz_0, ta1_x_xxxxx_yyyz_1, ta1_x_xxxxx_yyz_0, ta1_x_xxxxx_yyz_1, ta1_x_xxxxx_yyzz_0, ta1_x_xxxxx_yyzz_1, ta1_x_xxxxx_yzz_0, ta1_x_xxxxx_yzz_1, ta1_x_xxxxx_yzzz_0, ta1_x_xxxxx_yzzz_1, ta1_x_xxxxx_zzz_0, ta1_x_xxxxx_zzz_1, ta1_x_xxxxx_zzzz_0, ta1_x_xxxxx_zzzz_1, ta1_x_xxxxxy_xxxx_0, ta1_x_xxxxxy_xxxy_0, ta1_x_xxxxxy_xxxz_0, ta1_x_xxxxxy_xxyy_0, ta1_x_xxxxxy_xxyz_0, ta1_x_xxxxxy_xxzz_0, ta1_x_xxxxxy_xyyy_0, ta1_x_xxxxxy_xyyz_0, ta1_x_xxxxxy_xyzz_0, ta1_x_xxxxxy_xzzz_0, ta1_x_xxxxxy_yyyy_0, ta1_x_xxxxxy_yyyz_0, ta1_x_xxxxxy_yyzz_0, ta1_x_xxxxxy_yzzz_0, ta1_x_xxxxxy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxy_xxxx_0[i] = ta1_x_xxxxx_xxxx_0[i] * pa_y[i] - ta1_x_xxxxx_xxxx_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxxy_0[i] = ta1_x_xxxxx_xxx_0[i] * fe_0 - ta1_x_xxxxx_xxx_1[i] * fe_0 + ta1_x_xxxxx_xxxy_0[i] * pa_y[i] - ta1_x_xxxxx_xxxy_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxxz_0[i] = ta1_x_xxxxx_xxxz_0[i] * pa_y[i] - ta1_x_xxxxx_xxxz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxyy_0[i] = 2.0 * ta1_x_xxxxx_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xxy_1[i] * fe_0 + ta1_x_xxxxx_xxyy_0[i] * pa_y[i] - ta1_x_xxxxx_xxyy_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxyz_0[i] = ta1_x_xxxxx_xxz_0[i] * fe_0 - ta1_x_xxxxx_xxz_1[i] * fe_0 + ta1_x_xxxxx_xxyz_0[i] * pa_y[i] - ta1_x_xxxxx_xxyz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxzz_0[i] = ta1_x_xxxxx_xxzz_0[i] * pa_y[i] - ta1_x_xxxxx_xxzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xyyy_0[i] = 3.0 * ta1_x_xxxxx_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_xyy_1[i] * fe_0 + ta1_x_xxxxx_xyyy_0[i] * pa_y[i] - ta1_x_xxxxx_xyyy_1[i] * pc_y[i];

        ta1_x_xxxxxy_xyyz_0[i] = 2.0 * ta1_x_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xyz_1[i] * fe_0 + ta1_x_xxxxx_xyyz_0[i] * pa_y[i] - ta1_x_xxxxx_xyyz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xyzz_0[i] = ta1_x_xxxxx_xzz_0[i] * fe_0 - ta1_x_xxxxx_xzz_1[i] * fe_0 + ta1_x_xxxxx_xyzz_0[i] * pa_y[i] - ta1_x_xxxxx_xyzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xzzz_0[i] = ta1_x_xxxxx_xzzz_0[i] * pa_y[i] - ta1_x_xxxxx_xzzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_yyyy_0[i] = 4.0 * ta1_x_xxxxx_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxxxx_yyy_1[i] * fe_0 + ta1_x_xxxxx_yyyy_0[i] * pa_y[i] - ta1_x_xxxxx_yyyy_1[i] * pc_y[i];

        ta1_x_xxxxxy_yyyz_0[i] = 3.0 * ta1_x_xxxxx_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_yyz_1[i] * fe_0 + ta1_x_xxxxx_yyyz_0[i] * pa_y[i] - ta1_x_xxxxx_yyyz_1[i] * pc_y[i];

        ta1_x_xxxxxy_yyzz_0[i] = 2.0 * ta1_x_xxxxx_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_yzz_1[i] * fe_0 + ta1_x_xxxxx_yyzz_0[i] * pa_y[i] - ta1_x_xxxxx_yyzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_yzzz_0[i] = ta1_x_xxxxx_zzz_0[i] * fe_0 - ta1_x_xxxxx_zzz_1[i] * fe_0 + ta1_x_xxxxx_yzzz_0[i] * pa_y[i] - ta1_x_xxxxx_yzzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_zzzz_0[i] = ta1_x_xxxxx_zzzz_0[i] * pa_y[i] - ta1_x_xxxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : IG

    auto ta1_x_xxxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 30);

    auto ta1_x_xxxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 31);

    auto ta1_x_xxxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 32);

    auto ta1_x_xxxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 33);

    auto ta1_x_xxxxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 34);

    auto ta1_x_xxxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 35);

    auto ta1_x_xxxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 36);

    auto ta1_x_xxxxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 37);

    auto ta1_x_xxxxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 38);

    auto ta1_x_xxxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 39);

    auto ta1_x_xxxxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 40);

    auto ta1_x_xxxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 41);

    auto ta1_x_xxxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 42);

    auto ta1_x_xxxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 43);

    auto ta1_x_xxxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 44);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_xxxxx_xxx_0, ta1_x_xxxxx_xxx_1, ta1_x_xxxxx_xxxx_0, ta1_x_xxxxx_xxxx_1, ta1_x_xxxxx_xxxy_0, ta1_x_xxxxx_xxxy_1, ta1_x_xxxxx_xxxz_0, ta1_x_xxxxx_xxxz_1, ta1_x_xxxxx_xxy_0, ta1_x_xxxxx_xxy_1, ta1_x_xxxxx_xxyy_0, ta1_x_xxxxx_xxyy_1, ta1_x_xxxxx_xxyz_0, ta1_x_xxxxx_xxyz_1, ta1_x_xxxxx_xxz_0, ta1_x_xxxxx_xxz_1, ta1_x_xxxxx_xxzz_0, ta1_x_xxxxx_xxzz_1, ta1_x_xxxxx_xyy_0, ta1_x_xxxxx_xyy_1, ta1_x_xxxxx_xyyy_0, ta1_x_xxxxx_xyyy_1, ta1_x_xxxxx_xyyz_0, ta1_x_xxxxx_xyyz_1, ta1_x_xxxxx_xyz_0, ta1_x_xxxxx_xyz_1, ta1_x_xxxxx_xyzz_0, ta1_x_xxxxx_xyzz_1, ta1_x_xxxxx_xzz_0, ta1_x_xxxxx_xzz_1, ta1_x_xxxxx_xzzz_0, ta1_x_xxxxx_xzzz_1, ta1_x_xxxxx_yyy_0, ta1_x_xxxxx_yyy_1, ta1_x_xxxxx_yyyy_0, ta1_x_xxxxx_yyyy_1, ta1_x_xxxxx_yyyz_0, ta1_x_xxxxx_yyyz_1, ta1_x_xxxxx_yyz_0, ta1_x_xxxxx_yyz_1, ta1_x_xxxxx_yyzz_0, ta1_x_xxxxx_yyzz_1, ta1_x_xxxxx_yzz_0, ta1_x_xxxxx_yzz_1, ta1_x_xxxxx_yzzz_0, ta1_x_xxxxx_yzzz_1, ta1_x_xxxxx_zzz_0, ta1_x_xxxxx_zzz_1, ta1_x_xxxxx_zzzz_0, ta1_x_xxxxx_zzzz_1, ta1_x_xxxxxz_xxxx_0, ta1_x_xxxxxz_xxxy_0, ta1_x_xxxxxz_xxxz_0, ta1_x_xxxxxz_xxyy_0, ta1_x_xxxxxz_xxyz_0, ta1_x_xxxxxz_xxzz_0, ta1_x_xxxxxz_xyyy_0, ta1_x_xxxxxz_xyyz_0, ta1_x_xxxxxz_xyzz_0, ta1_x_xxxxxz_xzzz_0, ta1_x_xxxxxz_yyyy_0, ta1_x_xxxxxz_yyyz_0, ta1_x_xxxxxz_yyzz_0, ta1_x_xxxxxz_yzzz_0, ta1_x_xxxxxz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxz_xxxx_0[i] = ta1_x_xxxxx_xxxx_0[i] * pa_z[i] - ta1_x_xxxxx_xxxx_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxxy_0[i] = ta1_x_xxxxx_xxxy_0[i] * pa_z[i] - ta1_x_xxxxx_xxxy_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxxz_0[i] = ta1_x_xxxxx_xxx_0[i] * fe_0 - ta1_x_xxxxx_xxx_1[i] * fe_0 + ta1_x_xxxxx_xxxz_0[i] * pa_z[i] - ta1_x_xxxxx_xxxz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxyy_0[i] = ta1_x_xxxxx_xxyy_0[i] * pa_z[i] - ta1_x_xxxxx_xxyy_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxyz_0[i] = ta1_x_xxxxx_xxy_0[i] * fe_0 - ta1_x_xxxxx_xxy_1[i] * fe_0 + ta1_x_xxxxx_xxyz_0[i] * pa_z[i] - ta1_x_xxxxx_xxyz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxzz_0[i] = 2.0 * ta1_x_xxxxx_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xxz_1[i] * fe_0 + ta1_x_xxxxx_xxzz_0[i] * pa_z[i] - ta1_x_xxxxx_xxzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xyyy_0[i] = ta1_x_xxxxx_xyyy_0[i] * pa_z[i] - ta1_x_xxxxx_xyyy_1[i] * pc_z[i];

        ta1_x_xxxxxz_xyyz_0[i] = ta1_x_xxxxx_xyy_0[i] * fe_0 - ta1_x_xxxxx_xyy_1[i] * fe_0 + ta1_x_xxxxx_xyyz_0[i] * pa_z[i] - ta1_x_xxxxx_xyyz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xyzz_0[i] = 2.0 * ta1_x_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xyz_1[i] * fe_0 + ta1_x_xxxxx_xyzz_0[i] * pa_z[i] - ta1_x_xxxxx_xyzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xzzz_0[i] = 3.0 * ta1_x_xxxxx_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_xzz_1[i] * fe_0 + ta1_x_xxxxx_xzzz_0[i] * pa_z[i] - ta1_x_xxxxx_xzzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_yyyy_0[i] = ta1_x_xxxxx_yyyy_0[i] * pa_z[i] - ta1_x_xxxxx_yyyy_1[i] * pc_z[i];

        ta1_x_xxxxxz_yyyz_0[i] = ta1_x_xxxxx_yyy_0[i] * fe_0 - ta1_x_xxxxx_yyy_1[i] * fe_0 + ta1_x_xxxxx_yyyz_0[i] * pa_z[i] - ta1_x_xxxxx_yyyz_1[i] * pc_z[i];

        ta1_x_xxxxxz_yyzz_0[i] = 2.0 * ta1_x_xxxxx_yyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_yyz_1[i] * fe_0 + ta1_x_xxxxx_yyzz_0[i] * pa_z[i] - ta1_x_xxxxx_yyzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_yzzz_0[i] = 3.0 * ta1_x_xxxxx_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_yzz_1[i] * fe_0 + ta1_x_xxxxx_yzzz_0[i] * pa_z[i] - ta1_x_xxxxx_yzzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_zzzz_0[i] = 4.0 * ta1_x_xxxxx_zzz_0[i] * fe_0 - 4.0 * ta1_x_xxxxx_zzz_1[i] * fe_0 + ta1_x_xxxxx_zzzz_0[i] * pa_z[i] - ta1_x_xxxxx_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : IG

    auto ta1_x_xxxxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 45);

    auto ta1_x_xxxxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 46);

    auto ta1_x_xxxxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 47);

    auto ta1_x_xxxxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 48);

    auto ta1_x_xxxxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 49);

    auto ta1_x_xxxxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 50);

    auto ta1_x_xxxxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 51);

    auto ta1_x_xxxxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 52);

    auto ta1_x_xxxxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 53);

    auto ta1_x_xxxxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 54);

    auto ta1_x_xxxxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 55);

    auto ta1_x_xxxxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 56);

    auto ta1_x_xxxxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 57);

    auto ta1_x_xxxxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 58);

    auto ta1_x_xxxxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 59);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxxx_xxxx_0, ta1_x_xxxx_xxxx_1, ta1_x_xxxx_xxxy_0, ta1_x_xxxx_xxxy_1, ta1_x_xxxx_xxxz_0, ta1_x_xxxx_xxxz_1, ta1_x_xxxx_xxyy_0, ta1_x_xxxx_xxyy_1, ta1_x_xxxx_xxyz_0, ta1_x_xxxx_xxyz_1, ta1_x_xxxx_xxzz_0, ta1_x_xxxx_xxzz_1, ta1_x_xxxx_xyyy_0, ta1_x_xxxx_xyyy_1, ta1_x_xxxx_xyyz_0, ta1_x_xxxx_xyyz_1, ta1_x_xxxx_xyzz_0, ta1_x_xxxx_xyzz_1, ta1_x_xxxx_xzzz_0, ta1_x_xxxx_xzzz_1, ta1_x_xxxx_zzzz_0, ta1_x_xxxx_zzzz_1, ta1_x_xxxxy_xxx_0, ta1_x_xxxxy_xxx_1, ta1_x_xxxxy_xxxx_0, ta1_x_xxxxy_xxxx_1, ta1_x_xxxxy_xxxy_0, ta1_x_xxxxy_xxxy_1, ta1_x_xxxxy_xxxz_0, ta1_x_xxxxy_xxxz_1, ta1_x_xxxxy_xxy_0, ta1_x_xxxxy_xxy_1, ta1_x_xxxxy_xxyy_0, ta1_x_xxxxy_xxyy_1, ta1_x_xxxxy_xxyz_0, ta1_x_xxxxy_xxyz_1, ta1_x_xxxxy_xxz_0, ta1_x_xxxxy_xxz_1, ta1_x_xxxxy_xxzz_0, ta1_x_xxxxy_xxzz_1, ta1_x_xxxxy_xyy_0, ta1_x_xxxxy_xyy_1, ta1_x_xxxxy_xyyy_0, ta1_x_xxxxy_xyyy_1, ta1_x_xxxxy_xyyz_0, ta1_x_xxxxy_xyyz_1, ta1_x_xxxxy_xyz_0, ta1_x_xxxxy_xyz_1, ta1_x_xxxxy_xyzz_0, ta1_x_xxxxy_xyzz_1, ta1_x_xxxxy_xzz_0, ta1_x_xxxxy_xzz_1, ta1_x_xxxxy_xzzz_0, ta1_x_xxxxy_xzzz_1, ta1_x_xxxxy_zzzz_0, ta1_x_xxxxy_zzzz_1, ta1_x_xxxxyy_xxxx_0, ta1_x_xxxxyy_xxxy_0, ta1_x_xxxxyy_xxxz_0, ta1_x_xxxxyy_xxyy_0, ta1_x_xxxxyy_xxyz_0, ta1_x_xxxxyy_xxzz_0, ta1_x_xxxxyy_xyyy_0, ta1_x_xxxxyy_xyyz_0, ta1_x_xxxxyy_xyzz_0, ta1_x_xxxxyy_xzzz_0, ta1_x_xxxxyy_yyyy_0, ta1_x_xxxxyy_yyyz_0, ta1_x_xxxxyy_yyzz_0, ta1_x_xxxxyy_yzzz_0, ta1_x_xxxxyy_zzzz_0, ta1_x_xxxyy_yyyy_0, ta1_x_xxxyy_yyyy_1, ta1_x_xxxyy_yyyz_0, ta1_x_xxxyy_yyyz_1, ta1_x_xxxyy_yyzz_0, ta1_x_xxxyy_yyzz_1, ta1_x_xxxyy_yzzz_0, ta1_x_xxxyy_yzzz_1, ta1_x_xxyy_yyyy_0, ta1_x_xxyy_yyyy_1, ta1_x_xxyy_yyyz_0, ta1_x_xxyy_yyyz_1, ta1_x_xxyy_yyzz_0, ta1_x_xxyy_yyzz_1, ta1_x_xxyy_yzzz_0, ta1_x_xxyy_yzzz_1, ta_xxxyy_yyyy_1, ta_xxxyy_yyyz_1, ta_xxxyy_yyzz_1, ta_xxxyy_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyy_xxxx_0[i] = ta1_x_xxxx_xxxx_0[i] * fe_0 - ta1_x_xxxx_xxxx_1[i] * fe_0 + ta1_x_xxxxy_xxxx_0[i] * pa_y[i] - ta1_x_xxxxy_xxxx_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxxy_0[i] = ta1_x_xxxx_xxxy_0[i] * fe_0 - ta1_x_xxxx_xxxy_1[i] * fe_0 + ta1_x_xxxxy_xxx_0[i] * fe_0 - ta1_x_xxxxy_xxx_1[i] * fe_0 + ta1_x_xxxxy_xxxy_0[i] * pa_y[i] - ta1_x_xxxxy_xxxy_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxxz_0[i] = ta1_x_xxxx_xxxz_0[i] * fe_0 - ta1_x_xxxx_xxxz_1[i] * fe_0 + ta1_x_xxxxy_xxxz_0[i] * pa_y[i] - ta1_x_xxxxy_xxxz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxyy_0[i] = ta1_x_xxxx_xxyy_0[i] * fe_0 - ta1_x_xxxx_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxxxy_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxxy_xxy_1[i] * fe_0 + ta1_x_xxxxy_xxyy_0[i] * pa_y[i] - ta1_x_xxxxy_xxyy_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxyz_0[i] = ta1_x_xxxx_xxyz_0[i] * fe_0 - ta1_x_xxxx_xxyz_1[i] * fe_0 + ta1_x_xxxxy_xxz_0[i] * fe_0 - ta1_x_xxxxy_xxz_1[i] * fe_0 + ta1_x_xxxxy_xxyz_0[i] * pa_y[i] - ta1_x_xxxxy_xxyz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxzz_0[i] = ta1_x_xxxx_xxzz_0[i] * fe_0 - ta1_x_xxxx_xxzz_1[i] * fe_0 + ta1_x_xxxxy_xxzz_0[i] * pa_y[i] - ta1_x_xxxxy_xxzz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xyyy_0[i] = ta1_x_xxxx_xyyy_0[i] * fe_0 - ta1_x_xxxx_xyyy_1[i] * fe_0 + 3.0 * ta1_x_xxxxy_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxxxy_xyy_1[i] * fe_0 + ta1_x_xxxxy_xyyy_0[i] * pa_y[i] - ta1_x_xxxxy_xyyy_1[i] * pc_y[i];

        ta1_x_xxxxyy_xyyz_0[i] = ta1_x_xxxx_xyyz_0[i] * fe_0 - ta1_x_xxxx_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxxxy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxy_xyz_1[i] * fe_0 + ta1_x_xxxxy_xyyz_0[i] * pa_y[i] - ta1_x_xxxxy_xyyz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xyzz_0[i] = ta1_x_xxxx_xyzz_0[i] * fe_0 - ta1_x_xxxx_xyzz_1[i] * fe_0 + ta1_x_xxxxy_xzz_0[i] * fe_0 - ta1_x_xxxxy_xzz_1[i] * fe_0 + ta1_x_xxxxy_xyzz_0[i] * pa_y[i] - ta1_x_xxxxy_xyzz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xzzz_0[i] = ta1_x_xxxx_xzzz_0[i] * fe_0 - ta1_x_xxxx_xzzz_1[i] * fe_0 + ta1_x_xxxxy_xzzz_0[i] * pa_y[i] - ta1_x_xxxxy_xzzz_1[i] * pc_y[i];

        ta1_x_xxxxyy_yyyy_0[i] = 3.0 * ta1_x_xxyy_yyyy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yyyy_1[i] * fe_0 + ta_xxxyy_yyyy_1[i] + ta1_x_xxxyy_yyyy_0[i] * pa_x[i] - ta1_x_xxxyy_yyyy_1[i] * pc_x[i];

        ta1_x_xxxxyy_yyyz_0[i] = 3.0 * ta1_x_xxyy_yyyz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yyyz_1[i] * fe_0 + ta_xxxyy_yyyz_1[i] + ta1_x_xxxyy_yyyz_0[i] * pa_x[i] - ta1_x_xxxyy_yyyz_1[i] * pc_x[i];

        ta1_x_xxxxyy_yyzz_0[i] = 3.0 * ta1_x_xxyy_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yyzz_1[i] * fe_0 + ta_xxxyy_yyzz_1[i] + ta1_x_xxxyy_yyzz_0[i] * pa_x[i] - ta1_x_xxxyy_yyzz_1[i] * pc_x[i];

        ta1_x_xxxxyy_yzzz_0[i] = 3.0 * ta1_x_xxyy_yzzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yzzz_1[i] * fe_0 + ta_xxxyy_yzzz_1[i] + ta1_x_xxxyy_yzzz_0[i] * pa_x[i] - ta1_x_xxxyy_yzzz_1[i] * pc_x[i];

        ta1_x_xxxxyy_zzzz_0[i] = ta1_x_xxxx_zzzz_0[i] * fe_0 - ta1_x_xxxx_zzzz_1[i] * fe_0 + ta1_x_xxxxy_zzzz_0[i] * pa_y[i] - ta1_x_xxxxy_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : IG

    auto ta1_x_xxxxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 60);

    auto ta1_x_xxxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 61);

    auto ta1_x_xxxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 62);

    auto ta1_x_xxxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 63);

    auto ta1_x_xxxxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 64);

    auto ta1_x_xxxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 65);

    auto ta1_x_xxxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 66);

    auto ta1_x_xxxxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 67);

    auto ta1_x_xxxxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 68);

    auto ta1_x_xxxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 69);

    auto ta1_x_xxxxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 70);

    auto ta1_x_xxxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 71);

    auto ta1_x_xxxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 72);

    auto ta1_x_xxxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 73);

    auto ta1_x_xxxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 74);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxxxy_xxxy_0, ta1_x_xxxxy_xxxy_1, ta1_x_xxxxy_xxyy_0, ta1_x_xxxxy_xxyy_1, ta1_x_xxxxy_xyyy_0, ta1_x_xxxxy_xyyy_1, ta1_x_xxxxy_yyyy_0, ta1_x_xxxxy_yyyy_1, ta1_x_xxxxyz_xxxx_0, ta1_x_xxxxyz_xxxy_0, ta1_x_xxxxyz_xxxz_0, ta1_x_xxxxyz_xxyy_0, ta1_x_xxxxyz_xxyz_0, ta1_x_xxxxyz_xxzz_0, ta1_x_xxxxyz_xyyy_0, ta1_x_xxxxyz_xyyz_0, ta1_x_xxxxyz_xyzz_0, ta1_x_xxxxyz_xzzz_0, ta1_x_xxxxyz_yyyy_0, ta1_x_xxxxyz_yyyz_0, ta1_x_xxxxyz_yyzz_0, ta1_x_xxxxyz_yzzz_0, ta1_x_xxxxyz_zzzz_0, ta1_x_xxxxz_xxxx_0, ta1_x_xxxxz_xxxx_1, ta1_x_xxxxz_xxxz_0, ta1_x_xxxxz_xxxz_1, ta1_x_xxxxz_xxyz_0, ta1_x_xxxxz_xxyz_1, ta1_x_xxxxz_xxz_0, ta1_x_xxxxz_xxz_1, ta1_x_xxxxz_xxzz_0, ta1_x_xxxxz_xxzz_1, ta1_x_xxxxz_xyyz_0, ta1_x_xxxxz_xyyz_1, ta1_x_xxxxz_xyz_0, ta1_x_xxxxz_xyz_1, ta1_x_xxxxz_xyzz_0, ta1_x_xxxxz_xyzz_1, ta1_x_xxxxz_xzz_0, ta1_x_xxxxz_xzz_1, ta1_x_xxxxz_xzzz_0, ta1_x_xxxxz_xzzz_1, ta1_x_xxxxz_yyyz_0, ta1_x_xxxxz_yyyz_1, ta1_x_xxxxz_yyz_0, ta1_x_xxxxz_yyz_1, ta1_x_xxxxz_yyzz_0, ta1_x_xxxxz_yyzz_1, ta1_x_xxxxz_yzz_0, ta1_x_xxxxz_yzz_1, ta1_x_xxxxz_yzzz_0, ta1_x_xxxxz_yzzz_1, ta1_x_xxxxz_zzz_0, ta1_x_xxxxz_zzz_1, ta1_x_xxxxz_zzzz_0, ta1_x_xxxxz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyz_xxxx_0[i] = ta1_x_xxxxz_xxxx_0[i] * pa_y[i] - ta1_x_xxxxz_xxxx_1[i] * pc_y[i];

        ta1_x_xxxxyz_xxxy_0[i] = ta1_x_xxxxy_xxxy_0[i] * pa_z[i] - ta1_x_xxxxy_xxxy_1[i] * pc_z[i];

        ta1_x_xxxxyz_xxxz_0[i] = ta1_x_xxxxz_xxxz_0[i] * pa_y[i] - ta1_x_xxxxz_xxxz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xxyy_0[i] = ta1_x_xxxxy_xxyy_0[i] * pa_z[i] - ta1_x_xxxxy_xxyy_1[i] * pc_z[i];

        ta1_x_xxxxyz_xxyz_0[i] = ta1_x_xxxxz_xxz_0[i] * fe_0 - ta1_x_xxxxz_xxz_1[i] * fe_0 + ta1_x_xxxxz_xxyz_0[i] * pa_y[i] - ta1_x_xxxxz_xxyz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xxzz_0[i] = ta1_x_xxxxz_xxzz_0[i] * pa_y[i] - ta1_x_xxxxz_xxzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xyyy_0[i] = ta1_x_xxxxy_xyyy_0[i] * pa_z[i] - ta1_x_xxxxy_xyyy_1[i] * pc_z[i];

        ta1_x_xxxxyz_xyyz_0[i] = 2.0 * ta1_x_xxxxz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxz_xyz_1[i] * fe_0 + ta1_x_xxxxz_xyyz_0[i] * pa_y[i] - ta1_x_xxxxz_xyyz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xyzz_0[i] = ta1_x_xxxxz_xzz_0[i] * fe_0 - ta1_x_xxxxz_xzz_1[i] * fe_0 + ta1_x_xxxxz_xyzz_0[i] * pa_y[i] - ta1_x_xxxxz_xyzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xzzz_0[i] = ta1_x_xxxxz_xzzz_0[i] * pa_y[i] - ta1_x_xxxxz_xzzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_yyyy_0[i] = ta1_x_xxxxy_yyyy_0[i] * pa_z[i] - ta1_x_xxxxy_yyyy_1[i] * pc_z[i];

        ta1_x_xxxxyz_yyyz_0[i] = 3.0 * ta1_x_xxxxz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxxxz_yyz_1[i] * fe_0 + ta1_x_xxxxz_yyyz_0[i] * pa_y[i] - ta1_x_xxxxz_yyyz_1[i] * pc_y[i];

        ta1_x_xxxxyz_yyzz_0[i] = 2.0 * ta1_x_xxxxz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxxxz_yzz_1[i] * fe_0 + ta1_x_xxxxz_yyzz_0[i] * pa_y[i] - ta1_x_xxxxz_yyzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_yzzz_0[i] = ta1_x_xxxxz_zzz_0[i] * fe_0 - ta1_x_xxxxz_zzz_1[i] * fe_0 + ta1_x_xxxxz_yzzz_0[i] * pa_y[i] - ta1_x_xxxxz_yzzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_zzzz_0[i] = ta1_x_xxxxz_zzzz_0[i] * pa_y[i] - ta1_x_xxxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : IG

    auto ta1_x_xxxxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 75);

    auto ta1_x_xxxxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 76);

    auto ta1_x_xxxxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 77);

    auto ta1_x_xxxxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 78);

    auto ta1_x_xxxxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 79);

    auto ta1_x_xxxxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 80);

    auto ta1_x_xxxxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 81);

    auto ta1_x_xxxxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 82);

    auto ta1_x_xxxxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 83);

    auto ta1_x_xxxxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 84);

    auto ta1_x_xxxxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 85);

    auto ta1_x_xxxxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 86);

    auto ta1_x_xxxxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 87);

    auto ta1_x_xxxxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 88);

    auto ta1_x_xxxxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 89);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxxx_xxxx_0, ta1_x_xxxx_xxxx_1, ta1_x_xxxx_xxxy_0, ta1_x_xxxx_xxxy_1, ta1_x_xxxx_xxxz_0, ta1_x_xxxx_xxxz_1, ta1_x_xxxx_xxyy_0, ta1_x_xxxx_xxyy_1, ta1_x_xxxx_xxyz_0, ta1_x_xxxx_xxyz_1, ta1_x_xxxx_xxzz_0, ta1_x_xxxx_xxzz_1, ta1_x_xxxx_xyyy_0, ta1_x_xxxx_xyyy_1, ta1_x_xxxx_xyyz_0, ta1_x_xxxx_xyyz_1, ta1_x_xxxx_xyzz_0, ta1_x_xxxx_xyzz_1, ta1_x_xxxx_xzzz_0, ta1_x_xxxx_xzzz_1, ta1_x_xxxx_yyyy_0, ta1_x_xxxx_yyyy_1, ta1_x_xxxxz_xxx_0, ta1_x_xxxxz_xxx_1, ta1_x_xxxxz_xxxx_0, ta1_x_xxxxz_xxxx_1, ta1_x_xxxxz_xxxy_0, ta1_x_xxxxz_xxxy_1, ta1_x_xxxxz_xxxz_0, ta1_x_xxxxz_xxxz_1, ta1_x_xxxxz_xxy_0, ta1_x_xxxxz_xxy_1, ta1_x_xxxxz_xxyy_0, ta1_x_xxxxz_xxyy_1, ta1_x_xxxxz_xxyz_0, ta1_x_xxxxz_xxyz_1, ta1_x_xxxxz_xxz_0, ta1_x_xxxxz_xxz_1, ta1_x_xxxxz_xxzz_0, ta1_x_xxxxz_xxzz_1, ta1_x_xxxxz_xyy_0, ta1_x_xxxxz_xyy_1, ta1_x_xxxxz_xyyy_0, ta1_x_xxxxz_xyyy_1, ta1_x_xxxxz_xyyz_0, ta1_x_xxxxz_xyyz_1, ta1_x_xxxxz_xyz_0, ta1_x_xxxxz_xyz_1, ta1_x_xxxxz_xyzz_0, ta1_x_xxxxz_xyzz_1, ta1_x_xxxxz_xzz_0, ta1_x_xxxxz_xzz_1, ta1_x_xxxxz_xzzz_0, ta1_x_xxxxz_xzzz_1, ta1_x_xxxxz_yyyy_0, ta1_x_xxxxz_yyyy_1, ta1_x_xxxxzz_xxxx_0, ta1_x_xxxxzz_xxxy_0, ta1_x_xxxxzz_xxxz_0, ta1_x_xxxxzz_xxyy_0, ta1_x_xxxxzz_xxyz_0, ta1_x_xxxxzz_xxzz_0, ta1_x_xxxxzz_xyyy_0, ta1_x_xxxxzz_xyyz_0, ta1_x_xxxxzz_xyzz_0, ta1_x_xxxxzz_xzzz_0, ta1_x_xxxxzz_yyyy_0, ta1_x_xxxxzz_yyyz_0, ta1_x_xxxxzz_yyzz_0, ta1_x_xxxxzz_yzzz_0, ta1_x_xxxxzz_zzzz_0, ta1_x_xxxzz_yyyz_0, ta1_x_xxxzz_yyyz_1, ta1_x_xxxzz_yyzz_0, ta1_x_xxxzz_yyzz_1, ta1_x_xxxzz_yzzz_0, ta1_x_xxxzz_yzzz_1, ta1_x_xxxzz_zzzz_0, ta1_x_xxxzz_zzzz_1, ta1_x_xxzz_yyyz_0, ta1_x_xxzz_yyyz_1, ta1_x_xxzz_yyzz_0, ta1_x_xxzz_yyzz_1, ta1_x_xxzz_yzzz_0, ta1_x_xxzz_yzzz_1, ta1_x_xxzz_zzzz_0, ta1_x_xxzz_zzzz_1, ta_xxxzz_yyyz_1, ta_xxxzz_yyzz_1, ta_xxxzz_yzzz_1, ta_xxxzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxzz_xxxx_0[i] = ta1_x_xxxx_xxxx_0[i] * fe_0 - ta1_x_xxxx_xxxx_1[i] * fe_0 + ta1_x_xxxxz_xxxx_0[i] * pa_z[i] - ta1_x_xxxxz_xxxx_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxxy_0[i] = ta1_x_xxxx_xxxy_0[i] * fe_0 - ta1_x_xxxx_xxxy_1[i] * fe_0 + ta1_x_xxxxz_xxxy_0[i] * pa_z[i] - ta1_x_xxxxz_xxxy_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxxz_0[i] = ta1_x_xxxx_xxxz_0[i] * fe_0 - ta1_x_xxxx_xxxz_1[i] * fe_0 + ta1_x_xxxxz_xxx_0[i] * fe_0 - ta1_x_xxxxz_xxx_1[i] * fe_0 + ta1_x_xxxxz_xxxz_0[i] * pa_z[i] - ta1_x_xxxxz_xxxz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxyy_0[i] = ta1_x_xxxx_xxyy_0[i] * fe_0 - ta1_x_xxxx_xxyy_1[i] * fe_0 + ta1_x_xxxxz_xxyy_0[i] * pa_z[i] - ta1_x_xxxxz_xxyy_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxyz_0[i] = ta1_x_xxxx_xxyz_0[i] * fe_0 - ta1_x_xxxx_xxyz_1[i] * fe_0 + ta1_x_xxxxz_xxy_0[i] * fe_0 - ta1_x_xxxxz_xxy_1[i] * fe_0 + ta1_x_xxxxz_xxyz_0[i] * pa_z[i] - ta1_x_xxxxz_xxyz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxzz_0[i] = ta1_x_xxxx_xxzz_0[i] * fe_0 - ta1_x_xxxx_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxxxz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxxxz_xxz_1[i] * fe_0 + ta1_x_xxxxz_xxzz_0[i] * pa_z[i] - ta1_x_xxxxz_xxzz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xyyy_0[i] = ta1_x_xxxx_xyyy_0[i] * fe_0 - ta1_x_xxxx_xyyy_1[i] * fe_0 + ta1_x_xxxxz_xyyy_0[i] * pa_z[i] - ta1_x_xxxxz_xyyy_1[i] * pc_z[i];

        ta1_x_xxxxzz_xyyz_0[i] = ta1_x_xxxx_xyyz_0[i] * fe_0 - ta1_x_xxxx_xyyz_1[i] * fe_0 + ta1_x_xxxxz_xyy_0[i] * fe_0 - ta1_x_xxxxz_xyy_1[i] * fe_0 + ta1_x_xxxxz_xyyz_0[i] * pa_z[i] - ta1_x_xxxxz_xyyz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xyzz_0[i] = ta1_x_xxxx_xyzz_0[i] * fe_0 - ta1_x_xxxx_xyzz_1[i] * fe_0 + 2.0 * ta1_x_xxxxz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxxz_xyz_1[i] * fe_0 + ta1_x_xxxxz_xyzz_0[i] * pa_z[i] - ta1_x_xxxxz_xyzz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xzzz_0[i] = ta1_x_xxxx_xzzz_0[i] * fe_0 - ta1_x_xxxx_xzzz_1[i] * fe_0 + 3.0 * ta1_x_xxxxz_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxxxz_xzz_1[i] * fe_0 + ta1_x_xxxxz_xzzz_0[i] * pa_z[i] - ta1_x_xxxxz_xzzz_1[i] * pc_z[i];

        ta1_x_xxxxzz_yyyy_0[i] = ta1_x_xxxx_yyyy_0[i] * fe_0 - ta1_x_xxxx_yyyy_1[i] * fe_0 + ta1_x_xxxxz_yyyy_0[i] * pa_z[i] - ta1_x_xxxxz_yyyy_1[i] * pc_z[i];

        ta1_x_xxxxzz_yyyz_0[i] = 3.0 * ta1_x_xxzz_yyyz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yyyz_1[i] * fe_0 + ta_xxxzz_yyyz_1[i] + ta1_x_xxxzz_yyyz_0[i] * pa_x[i] - ta1_x_xxxzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxxxzz_yyzz_0[i] = 3.0 * ta1_x_xxzz_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yyzz_1[i] * fe_0 + ta_xxxzz_yyzz_1[i] + ta1_x_xxxzz_yyzz_0[i] * pa_x[i] - ta1_x_xxxzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxxxzz_yzzz_0[i] = 3.0 * ta1_x_xxzz_yzzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yzzz_1[i] * fe_0 + ta_xxxzz_yzzz_1[i] + ta1_x_xxxzz_yzzz_0[i] * pa_x[i] - ta1_x_xxxzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxxxzz_zzzz_0[i] = 3.0 * ta1_x_xxzz_zzzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_zzzz_1[i] * fe_0 + ta_xxxzz_zzzz_1[i] + ta1_x_xxxzz_zzzz_0[i] * pa_x[i] - ta1_x_xxxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : IG

    auto ta1_x_xxxyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 90);

    auto ta1_x_xxxyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 91);

    auto ta1_x_xxxyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 92);

    auto ta1_x_xxxyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 93);

    auto ta1_x_xxxyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 94);

    auto ta1_x_xxxyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 95);

    auto ta1_x_xxxyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 96);

    auto ta1_x_xxxyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 97);

    auto ta1_x_xxxyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 98);

    auto ta1_x_xxxyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 99);

    auto ta1_x_xxxyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 100);

    auto ta1_x_xxxyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 101);

    auto ta1_x_xxxyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 102);

    auto ta1_x_xxxyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 103);

    auto ta1_x_xxxyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 104);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxxy_xxxx_0, ta1_x_xxxy_xxxx_1, ta1_x_xxxy_xxxy_0, ta1_x_xxxy_xxxy_1, ta1_x_xxxy_xxxz_0, ta1_x_xxxy_xxxz_1, ta1_x_xxxy_xxyy_0, ta1_x_xxxy_xxyy_1, ta1_x_xxxy_xxyz_0, ta1_x_xxxy_xxyz_1, ta1_x_xxxy_xxzz_0, ta1_x_xxxy_xxzz_1, ta1_x_xxxy_xyyy_0, ta1_x_xxxy_xyyy_1, ta1_x_xxxy_xyyz_0, ta1_x_xxxy_xyyz_1, ta1_x_xxxy_xyzz_0, ta1_x_xxxy_xyzz_1, ta1_x_xxxy_xzzz_0, ta1_x_xxxy_xzzz_1, ta1_x_xxxy_zzzz_0, ta1_x_xxxy_zzzz_1, ta1_x_xxxyy_xxx_0, ta1_x_xxxyy_xxx_1, ta1_x_xxxyy_xxxx_0, ta1_x_xxxyy_xxxx_1, ta1_x_xxxyy_xxxy_0, ta1_x_xxxyy_xxxy_1, ta1_x_xxxyy_xxxz_0, ta1_x_xxxyy_xxxz_1, ta1_x_xxxyy_xxy_0, ta1_x_xxxyy_xxy_1, ta1_x_xxxyy_xxyy_0, ta1_x_xxxyy_xxyy_1, ta1_x_xxxyy_xxyz_0, ta1_x_xxxyy_xxyz_1, ta1_x_xxxyy_xxz_0, ta1_x_xxxyy_xxz_1, ta1_x_xxxyy_xxzz_0, ta1_x_xxxyy_xxzz_1, ta1_x_xxxyy_xyy_0, ta1_x_xxxyy_xyy_1, ta1_x_xxxyy_xyyy_0, ta1_x_xxxyy_xyyy_1, ta1_x_xxxyy_xyyz_0, ta1_x_xxxyy_xyyz_1, ta1_x_xxxyy_xyz_0, ta1_x_xxxyy_xyz_1, ta1_x_xxxyy_xyzz_0, ta1_x_xxxyy_xyzz_1, ta1_x_xxxyy_xzz_0, ta1_x_xxxyy_xzz_1, ta1_x_xxxyy_xzzz_0, ta1_x_xxxyy_xzzz_1, ta1_x_xxxyy_zzzz_0, ta1_x_xxxyy_zzzz_1, ta1_x_xxxyyy_xxxx_0, ta1_x_xxxyyy_xxxy_0, ta1_x_xxxyyy_xxxz_0, ta1_x_xxxyyy_xxyy_0, ta1_x_xxxyyy_xxyz_0, ta1_x_xxxyyy_xxzz_0, ta1_x_xxxyyy_xyyy_0, ta1_x_xxxyyy_xyyz_0, ta1_x_xxxyyy_xyzz_0, ta1_x_xxxyyy_xzzz_0, ta1_x_xxxyyy_yyyy_0, ta1_x_xxxyyy_yyyz_0, ta1_x_xxxyyy_yyzz_0, ta1_x_xxxyyy_yzzz_0, ta1_x_xxxyyy_zzzz_0, ta1_x_xxyyy_yyyy_0, ta1_x_xxyyy_yyyy_1, ta1_x_xxyyy_yyyz_0, ta1_x_xxyyy_yyyz_1, ta1_x_xxyyy_yyzz_0, ta1_x_xxyyy_yyzz_1, ta1_x_xxyyy_yzzz_0, ta1_x_xxyyy_yzzz_1, ta1_x_xyyy_yyyy_0, ta1_x_xyyy_yyyy_1, ta1_x_xyyy_yyyz_0, ta1_x_xyyy_yyyz_1, ta1_x_xyyy_yyzz_0, ta1_x_xyyy_yyzz_1, ta1_x_xyyy_yzzz_0, ta1_x_xyyy_yzzz_1, ta_xxyyy_yyyy_1, ta_xxyyy_yyyz_1, ta_xxyyy_yyzz_1, ta_xxyyy_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyy_xxxx_0[i] = 2.0 * ta1_x_xxxy_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxxx_1[i] * fe_0 + ta1_x_xxxyy_xxxx_0[i] * pa_y[i] - ta1_x_xxxyy_xxxx_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxxy_0[i] = 2.0 * ta1_x_xxxy_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxxy_1[i] * fe_0 + ta1_x_xxxyy_xxx_0[i] * fe_0 - ta1_x_xxxyy_xxx_1[i] * fe_0 + ta1_x_xxxyy_xxxy_0[i] * pa_y[i] - ta1_x_xxxyy_xxxy_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxxz_0[i] = 2.0 * ta1_x_xxxy_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxxz_1[i] * fe_0 + ta1_x_xxxyy_xxxz_0[i] * pa_y[i] - ta1_x_xxxyy_xxxz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxyy_0[i] = 2.0 * ta1_x_xxxy_xxyy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxxyy_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxyy_xxy_1[i] * fe_0 + ta1_x_xxxyy_xxyy_0[i] * pa_y[i] - ta1_x_xxxyy_xxyy_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxyz_0[i] = 2.0 * ta1_x_xxxy_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxyz_1[i] * fe_0 + ta1_x_xxxyy_xxz_0[i] * fe_0 - ta1_x_xxxyy_xxz_1[i] * fe_0 + ta1_x_xxxyy_xxyz_0[i] * pa_y[i] - ta1_x_xxxyy_xxyz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxzz_0[i] = 2.0 * ta1_x_xxxy_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxzz_1[i] * fe_0 + ta1_x_xxxyy_xxzz_0[i] * pa_y[i] - ta1_x_xxxyy_xxzz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xyyy_0[i] = 2.0 * ta1_x_xxxy_xyyy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xyyy_1[i] * fe_0 + 3.0 * ta1_x_xxxyy_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxxyy_xyy_1[i] * fe_0 + ta1_x_xxxyy_xyyy_0[i] * pa_y[i] - ta1_x_xxxyy_xyyy_1[i] * pc_y[i];

        ta1_x_xxxyyy_xyyz_0[i] = 2.0 * ta1_x_xxxy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxyy_xyz_1[i] * fe_0 + ta1_x_xxxyy_xyyz_0[i] * pa_y[i] - ta1_x_xxxyy_xyyz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xyzz_0[i] = 2.0 * ta1_x_xxxy_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xyzz_1[i] * fe_0 + ta1_x_xxxyy_xzz_0[i] * fe_0 - ta1_x_xxxyy_xzz_1[i] * fe_0 + ta1_x_xxxyy_xyzz_0[i] * pa_y[i] - ta1_x_xxxyy_xyzz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xzzz_0[i] = 2.0 * ta1_x_xxxy_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xzzz_1[i] * fe_0 + ta1_x_xxxyy_xzzz_0[i] * pa_y[i] - ta1_x_xxxyy_xzzz_1[i] * pc_y[i];

        ta1_x_xxxyyy_yyyy_0[i] = 2.0 * ta1_x_xyyy_yyyy_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yyyy_1[i] * fe_0 + ta_xxyyy_yyyy_1[i] + ta1_x_xxyyy_yyyy_0[i] * pa_x[i] - ta1_x_xxyyy_yyyy_1[i] * pc_x[i];

        ta1_x_xxxyyy_yyyz_0[i] = 2.0 * ta1_x_xyyy_yyyz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yyyz_1[i] * fe_0 + ta_xxyyy_yyyz_1[i] + ta1_x_xxyyy_yyyz_0[i] * pa_x[i] - ta1_x_xxyyy_yyyz_1[i] * pc_x[i];

        ta1_x_xxxyyy_yyzz_0[i] = 2.0 * ta1_x_xyyy_yyzz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yyzz_1[i] * fe_0 + ta_xxyyy_yyzz_1[i] + ta1_x_xxyyy_yyzz_0[i] * pa_x[i] - ta1_x_xxyyy_yyzz_1[i] * pc_x[i];

        ta1_x_xxxyyy_yzzz_0[i] = 2.0 * ta1_x_xyyy_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yzzz_1[i] * fe_0 + ta_xxyyy_yzzz_1[i] + ta1_x_xxyyy_yzzz_0[i] * pa_x[i] - ta1_x_xxyyy_yzzz_1[i] * pc_x[i];

        ta1_x_xxxyyy_zzzz_0[i] = 2.0 * ta1_x_xxxy_zzzz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_zzzz_1[i] * fe_0 + ta1_x_xxxyy_zzzz_0[i] * pa_y[i] - ta1_x_xxxyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 105-120 components of targeted buffer : IG

    auto ta1_x_xxxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 105);

    auto ta1_x_xxxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 106);

    auto ta1_x_xxxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 107);

    auto ta1_x_xxxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 108);

    auto ta1_x_xxxyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 109);

    auto ta1_x_xxxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 110);

    auto ta1_x_xxxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 111);

    auto ta1_x_xxxyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 112);

    auto ta1_x_xxxyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 113);

    auto ta1_x_xxxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 114);

    auto ta1_x_xxxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 115);

    auto ta1_x_xxxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 116);

    auto ta1_x_xxxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 117);

    auto ta1_x_xxxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 118);

    auto ta1_x_xxxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 119);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxxyy_xxxx_0, ta1_x_xxxyy_xxxx_1, ta1_x_xxxyy_xxxy_0, ta1_x_xxxyy_xxxy_1, ta1_x_xxxyy_xxy_0, ta1_x_xxxyy_xxy_1, ta1_x_xxxyy_xxyy_0, ta1_x_xxxyy_xxyy_1, ta1_x_xxxyy_xxyz_0, ta1_x_xxxyy_xxyz_1, ta1_x_xxxyy_xyy_0, ta1_x_xxxyy_xyy_1, ta1_x_xxxyy_xyyy_0, ta1_x_xxxyy_xyyy_1, ta1_x_xxxyy_xyyz_0, ta1_x_xxxyy_xyyz_1, ta1_x_xxxyy_xyz_0, ta1_x_xxxyy_xyz_1, ta1_x_xxxyy_xyzz_0, ta1_x_xxxyy_xyzz_1, ta1_x_xxxyy_yyy_0, ta1_x_xxxyy_yyy_1, ta1_x_xxxyy_yyyy_0, ta1_x_xxxyy_yyyy_1, ta1_x_xxxyy_yyyz_0, ta1_x_xxxyy_yyyz_1, ta1_x_xxxyy_yyz_0, ta1_x_xxxyy_yyz_1, ta1_x_xxxyy_yyzz_0, ta1_x_xxxyy_yyzz_1, ta1_x_xxxyy_yzz_0, ta1_x_xxxyy_yzz_1, ta1_x_xxxyy_yzzz_0, ta1_x_xxxyy_yzzz_1, ta1_x_xxxyyz_xxxx_0, ta1_x_xxxyyz_xxxy_0, ta1_x_xxxyyz_xxxz_0, ta1_x_xxxyyz_xxyy_0, ta1_x_xxxyyz_xxyz_0, ta1_x_xxxyyz_xxzz_0, ta1_x_xxxyyz_xyyy_0, ta1_x_xxxyyz_xyyz_0, ta1_x_xxxyyz_xyzz_0, ta1_x_xxxyyz_xzzz_0, ta1_x_xxxyyz_yyyy_0, ta1_x_xxxyyz_yyyz_0, ta1_x_xxxyyz_yyzz_0, ta1_x_xxxyyz_yzzz_0, ta1_x_xxxyyz_zzzz_0, ta1_x_xxxyz_xxxz_0, ta1_x_xxxyz_xxxz_1, ta1_x_xxxyz_xxzz_0, ta1_x_xxxyz_xxzz_1, ta1_x_xxxyz_xzzz_0, ta1_x_xxxyz_xzzz_1, ta1_x_xxxyz_zzzz_0, ta1_x_xxxyz_zzzz_1, ta1_x_xxxz_xxxz_0, ta1_x_xxxz_xxxz_1, ta1_x_xxxz_xxzz_0, ta1_x_xxxz_xxzz_1, ta1_x_xxxz_xzzz_0, ta1_x_xxxz_xzzz_1, ta1_x_xxxz_zzzz_0, ta1_x_xxxz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyz_xxxx_0[i] = ta1_x_xxxyy_xxxx_0[i] * pa_z[i] - ta1_x_xxxyy_xxxx_1[i] * pc_z[i];

        ta1_x_xxxyyz_xxxy_0[i] = ta1_x_xxxyy_xxxy_0[i] * pa_z[i] - ta1_x_xxxyy_xxxy_1[i] * pc_z[i];

        ta1_x_xxxyyz_xxxz_0[i] = ta1_x_xxxz_xxxz_0[i] * fe_0 - ta1_x_xxxz_xxxz_1[i] * fe_0 + ta1_x_xxxyz_xxxz_0[i] * pa_y[i] - ta1_x_xxxyz_xxxz_1[i] * pc_y[i];

        ta1_x_xxxyyz_xxyy_0[i] = ta1_x_xxxyy_xxyy_0[i] * pa_z[i] - ta1_x_xxxyy_xxyy_1[i] * pc_z[i];

        ta1_x_xxxyyz_xxyz_0[i] = ta1_x_xxxyy_xxy_0[i] * fe_0 - ta1_x_xxxyy_xxy_1[i] * fe_0 + ta1_x_xxxyy_xxyz_0[i] * pa_z[i] - ta1_x_xxxyy_xxyz_1[i] * pc_z[i];

        ta1_x_xxxyyz_xxzz_0[i] = ta1_x_xxxz_xxzz_0[i] * fe_0 - ta1_x_xxxz_xxzz_1[i] * fe_0 + ta1_x_xxxyz_xxzz_0[i] * pa_y[i] - ta1_x_xxxyz_xxzz_1[i] * pc_y[i];

        ta1_x_xxxyyz_xyyy_0[i] = ta1_x_xxxyy_xyyy_0[i] * pa_z[i] - ta1_x_xxxyy_xyyy_1[i] * pc_z[i];

        ta1_x_xxxyyz_xyyz_0[i] = ta1_x_xxxyy_xyy_0[i] * fe_0 - ta1_x_xxxyy_xyy_1[i] * fe_0 + ta1_x_xxxyy_xyyz_0[i] * pa_z[i] - ta1_x_xxxyy_xyyz_1[i] * pc_z[i];

        ta1_x_xxxyyz_xyzz_0[i] = 2.0 * ta1_x_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxyy_xyz_1[i] * fe_0 + ta1_x_xxxyy_xyzz_0[i] * pa_z[i] - ta1_x_xxxyy_xyzz_1[i] * pc_z[i];

        ta1_x_xxxyyz_xzzz_0[i] = ta1_x_xxxz_xzzz_0[i] * fe_0 - ta1_x_xxxz_xzzz_1[i] * fe_0 + ta1_x_xxxyz_xzzz_0[i] * pa_y[i] - ta1_x_xxxyz_xzzz_1[i] * pc_y[i];

        ta1_x_xxxyyz_yyyy_0[i] = ta1_x_xxxyy_yyyy_0[i] * pa_z[i] - ta1_x_xxxyy_yyyy_1[i] * pc_z[i];

        ta1_x_xxxyyz_yyyz_0[i] = ta1_x_xxxyy_yyy_0[i] * fe_0 - ta1_x_xxxyy_yyy_1[i] * fe_0 + ta1_x_xxxyy_yyyz_0[i] * pa_z[i] - ta1_x_xxxyy_yyyz_1[i] * pc_z[i];

        ta1_x_xxxyyz_yyzz_0[i] = 2.0 * ta1_x_xxxyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_xxxyy_yyz_1[i] * fe_0 + ta1_x_xxxyy_yyzz_0[i] * pa_z[i] - ta1_x_xxxyy_yyzz_1[i] * pc_z[i];

        ta1_x_xxxyyz_yzzz_0[i] = 3.0 * ta1_x_xxxyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxxyy_yzz_1[i] * fe_0 + ta1_x_xxxyy_yzzz_0[i] * pa_z[i] - ta1_x_xxxyy_yzzz_1[i] * pc_z[i];

        ta1_x_xxxyyz_zzzz_0[i] = ta1_x_xxxz_zzzz_0[i] * fe_0 - ta1_x_xxxz_zzzz_1[i] * fe_0 + ta1_x_xxxyz_zzzz_0[i] * pa_y[i] - ta1_x_xxxyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : IG

    auto ta1_x_xxxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 120);

    auto ta1_x_xxxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 121);

    auto ta1_x_xxxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 122);

    auto ta1_x_xxxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 123);

    auto ta1_x_xxxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 124);

    auto ta1_x_xxxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 125);

    auto ta1_x_xxxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 126);

    auto ta1_x_xxxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 127);

    auto ta1_x_xxxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 128);

    auto ta1_x_xxxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 129);

    auto ta1_x_xxxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 130);

    auto ta1_x_xxxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 131);

    auto ta1_x_xxxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 132);

    auto ta1_x_xxxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 133);

    auto ta1_x_xxxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 134);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxxyzz_xxxx_0, ta1_x_xxxyzz_xxxy_0, ta1_x_xxxyzz_xxxz_0, ta1_x_xxxyzz_xxyy_0, ta1_x_xxxyzz_xxyz_0, ta1_x_xxxyzz_xxzz_0, ta1_x_xxxyzz_xyyy_0, ta1_x_xxxyzz_xyyz_0, ta1_x_xxxyzz_xyzz_0, ta1_x_xxxyzz_xzzz_0, ta1_x_xxxyzz_yyyy_0, ta1_x_xxxyzz_yyyz_0, ta1_x_xxxyzz_yyzz_0, ta1_x_xxxyzz_yzzz_0, ta1_x_xxxyzz_zzzz_0, ta1_x_xxxzz_xxx_0, ta1_x_xxxzz_xxx_1, ta1_x_xxxzz_xxxx_0, ta1_x_xxxzz_xxxx_1, ta1_x_xxxzz_xxxy_0, ta1_x_xxxzz_xxxy_1, ta1_x_xxxzz_xxxz_0, ta1_x_xxxzz_xxxz_1, ta1_x_xxxzz_xxy_0, ta1_x_xxxzz_xxy_1, ta1_x_xxxzz_xxyy_0, ta1_x_xxxzz_xxyy_1, ta1_x_xxxzz_xxyz_0, ta1_x_xxxzz_xxyz_1, ta1_x_xxxzz_xxz_0, ta1_x_xxxzz_xxz_1, ta1_x_xxxzz_xxzz_0, ta1_x_xxxzz_xxzz_1, ta1_x_xxxzz_xyy_0, ta1_x_xxxzz_xyy_1, ta1_x_xxxzz_xyyy_0, ta1_x_xxxzz_xyyy_1, ta1_x_xxxzz_xyyz_0, ta1_x_xxxzz_xyyz_1, ta1_x_xxxzz_xyz_0, ta1_x_xxxzz_xyz_1, ta1_x_xxxzz_xyzz_0, ta1_x_xxxzz_xyzz_1, ta1_x_xxxzz_xzz_0, ta1_x_xxxzz_xzz_1, ta1_x_xxxzz_xzzz_0, ta1_x_xxxzz_xzzz_1, ta1_x_xxxzz_yyy_0, ta1_x_xxxzz_yyy_1, ta1_x_xxxzz_yyyy_0, ta1_x_xxxzz_yyyy_1, ta1_x_xxxzz_yyyz_0, ta1_x_xxxzz_yyyz_1, ta1_x_xxxzz_yyz_0, ta1_x_xxxzz_yyz_1, ta1_x_xxxzz_yyzz_0, ta1_x_xxxzz_yyzz_1, ta1_x_xxxzz_yzz_0, ta1_x_xxxzz_yzz_1, ta1_x_xxxzz_yzzz_0, ta1_x_xxxzz_yzzz_1, ta1_x_xxxzz_zzz_0, ta1_x_xxxzz_zzz_1, ta1_x_xxxzz_zzzz_0, ta1_x_xxxzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyzz_xxxx_0[i] = ta1_x_xxxzz_xxxx_0[i] * pa_y[i] - ta1_x_xxxzz_xxxx_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxxy_0[i] = ta1_x_xxxzz_xxx_0[i] * fe_0 - ta1_x_xxxzz_xxx_1[i] * fe_0 + ta1_x_xxxzz_xxxy_0[i] * pa_y[i] - ta1_x_xxxzz_xxxy_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxxz_0[i] = ta1_x_xxxzz_xxxz_0[i] * pa_y[i] - ta1_x_xxxzz_xxxz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxyy_0[i] = 2.0 * ta1_x_xxxzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_xxy_1[i] * fe_0 + ta1_x_xxxzz_xxyy_0[i] * pa_y[i] - ta1_x_xxxzz_xxyy_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxyz_0[i] = ta1_x_xxxzz_xxz_0[i] * fe_0 - ta1_x_xxxzz_xxz_1[i] * fe_0 + ta1_x_xxxzz_xxyz_0[i] * pa_y[i] - ta1_x_xxxzz_xxyz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxzz_0[i] = ta1_x_xxxzz_xxzz_0[i] * pa_y[i] - ta1_x_xxxzz_xxzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xyyy_0[i] = 3.0 * ta1_x_xxxzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxxzz_xyy_1[i] * fe_0 + ta1_x_xxxzz_xyyy_0[i] * pa_y[i] - ta1_x_xxxzz_xyyy_1[i] * pc_y[i];

        ta1_x_xxxyzz_xyyz_0[i] = 2.0 * ta1_x_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_xyz_1[i] * fe_0 + ta1_x_xxxzz_xyyz_0[i] * pa_y[i] - ta1_x_xxxzz_xyyz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xyzz_0[i] = ta1_x_xxxzz_xzz_0[i] * fe_0 - ta1_x_xxxzz_xzz_1[i] * fe_0 + ta1_x_xxxzz_xyzz_0[i] * pa_y[i] - ta1_x_xxxzz_xyzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xzzz_0[i] = ta1_x_xxxzz_xzzz_0[i] * pa_y[i] - ta1_x_xxxzz_xzzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_yyyy_0[i] = 4.0 * ta1_x_xxxzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxxzz_yyy_1[i] * fe_0 + ta1_x_xxxzz_yyyy_0[i] * pa_y[i] - ta1_x_xxxzz_yyyy_1[i] * pc_y[i];

        ta1_x_xxxyzz_yyyz_0[i] = 3.0 * ta1_x_xxxzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxxzz_yyz_1[i] * fe_0 + ta1_x_xxxzz_yyyz_0[i] * pa_y[i] - ta1_x_xxxzz_yyyz_1[i] * pc_y[i];

        ta1_x_xxxyzz_yyzz_0[i] = 2.0 * ta1_x_xxxzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_yzz_1[i] * fe_0 + ta1_x_xxxzz_yyzz_0[i] * pa_y[i] - ta1_x_xxxzz_yyzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_yzzz_0[i] = ta1_x_xxxzz_zzz_0[i] * fe_0 - ta1_x_xxxzz_zzz_1[i] * fe_0 + ta1_x_xxxzz_yzzz_0[i] * pa_y[i] - ta1_x_xxxzz_yzzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_zzzz_0[i] = ta1_x_xxxzz_zzzz_0[i] * pa_y[i] - ta1_x_xxxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : IG

    auto ta1_x_xxxzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 135);

    auto ta1_x_xxxzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 136);

    auto ta1_x_xxxzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 137);

    auto ta1_x_xxxzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 138);

    auto ta1_x_xxxzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 139);

    auto ta1_x_xxxzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 140);

    auto ta1_x_xxxzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 141);

    auto ta1_x_xxxzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 142);

    auto ta1_x_xxxzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 143);

    auto ta1_x_xxxzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 144);

    auto ta1_x_xxxzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 145);

    auto ta1_x_xxxzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 146);

    auto ta1_x_xxxzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 147);

    auto ta1_x_xxxzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 148);

    auto ta1_x_xxxzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 149);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxxz_xxxx_0, ta1_x_xxxz_xxxx_1, ta1_x_xxxz_xxxy_0, ta1_x_xxxz_xxxy_1, ta1_x_xxxz_xxxz_0, ta1_x_xxxz_xxxz_1, ta1_x_xxxz_xxyy_0, ta1_x_xxxz_xxyy_1, ta1_x_xxxz_xxyz_0, ta1_x_xxxz_xxyz_1, ta1_x_xxxz_xxzz_0, ta1_x_xxxz_xxzz_1, ta1_x_xxxz_xyyy_0, ta1_x_xxxz_xyyy_1, ta1_x_xxxz_xyyz_0, ta1_x_xxxz_xyyz_1, ta1_x_xxxz_xyzz_0, ta1_x_xxxz_xyzz_1, ta1_x_xxxz_xzzz_0, ta1_x_xxxz_xzzz_1, ta1_x_xxxz_yyyy_0, ta1_x_xxxz_yyyy_1, ta1_x_xxxzz_xxx_0, ta1_x_xxxzz_xxx_1, ta1_x_xxxzz_xxxx_0, ta1_x_xxxzz_xxxx_1, ta1_x_xxxzz_xxxy_0, ta1_x_xxxzz_xxxy_1, ta1_x_xxxzz_xxxz_0, ta1_x_xxxzz_xxxz_1, ta1_x_xxxzz_xxy_0, ta1_x_xxxzz_xxy_1, ta1_x_xxxzz_xxyy_0, ta1_x_xxxzz_xxyy_1, ta1_x_xxxzz_xxyz_0, ta1_x_xxxzz_xxyz_1, ta1_x_xxxzz_xxz_0, ta1_x_xxxzz_xxz_1, ta1_x_xxxzz_xxzz_0, ta1_x_xxxzz_xxzz_1, ta1_x_xxxzz_xyy_0, ta1_x_xxxzz_xyy_1, ta1_x_xxxzz_xyyy_0, ta1_x_xxxzz_xyyy_1, ta1_x_xxxzz_xyyz_0, ta1_x_xxxzz_xyyz_1, ta1_x_xxxzz_xyz_0, ta1_x_xxxzz_xyz_1, ta1_x_xxxzz_xyzz_0, ta1_x_xxxzz_xyzz_1, ta1_x_xxxzz_xzz_0, ta1_x_xxxzz_xzz_1, ta1_x_xxxzz_xzzz_0, ta1_x_xxxzz_xzzz_1, ta1_x_xxxzz_yyyy_0, ta1_x_xxxzz_yyyy_1, ta1_x_xxxzzz_xxxx_0, ta1_x_xxxzzz_xxxy_0, ta1_x_xxxzzz_xxxz_0, ta1_x_xxxzzz_xxyy_0, ta1_x_xxxzzz_xxyz_0, ta1_x_xxxzzz_xxzz_0, ta1_x_xxxzzz_xyyy_0, ta1_x_xxxzzz_xyyz_0, ta1_x_xxxzzz_xyzz_0, ta1_x_xxxzzz_xzzz_0, ta1_x_xxxzzz_yyyy_0, ta1_x_xxxzzz_yyyz_0, ta1_x_xxxzzz_yyzz_0, ta1_x_xxxzzz_yzzz_0, ta1_x_xxxzzz_zzzz_0, ta1_x_xxzzz_yyyz_0, ta1_x_xxzzz_yyyz_1, ta1_x_xxzzz_yyzz_0, ta1_x_xxzzz_yyzz_1, ta1_x_xxzzz_yzzz_0, ta1_x_xxzzz_yzzz_1, ta1_x_xxzzz_zzzz_0, ta1_x_xxzzz_zzzz_1, ta1_x_xzzz_yyyz_0, ta1_x_xzzz_yyyz_1, ta1_x_xzzz_yyzz_0, ta1_x_xzzz_yyzz_1, ta1_x_xzzz_yzzz_0, ta1_x_xzzz_yzzz_1, ta1_x_xzzz_zzzz_0, ta1_x_xzzz_zzzz_1, ta_xxzzz_yyyz_1, ta_xxzzz_yyzz_1, ta_xxzzz_yzzz_1, ta_xxzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzzz_xxxx_0[i] = 2.0 * ta1_x_xxxz_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxxx_1[i] * fe_0 + ta1_x_xxxzz_xxxx_0[i] * pa_z[i] - ta1_x_xxxzz_xxxx_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxxy_0[i] = 2.0 * ta1_x_xxxz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxxy_1[i] * fe_0 + ta1_x_xxxzz_xxxy_0[i] * pa_z[i] - ta1_x_xxxzz_xxxy_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxxz_0[i] = 2.0 * ta1_x_xxxz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxxz_1[i] * fe_0 + ta1_x_xxxzz_xxx_0[i] * fe_0 - ta1_x_xxxzz_xxx_1[i] * fe_0 + ta1_x_xxxzz_xxxz_0[i] * pa_z[i] - ta1_x_xxxzz_xxxz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxyy_0[i] = 2.0 * ta1_x_xxxz_xxyy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxyy_1[i] * fe_0 + ta1_x_xxxzz_xxyy_0[i] * pa_z[i] - ta1_x_xxxzz_xxyy_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxyz_0[i] = 2.0 * ta1_x_xxxz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxyz_1[i] * fe_0 + ta1_x_xxxzz_xxy_0[i] * fe_0 - ta1_x_xxxzz_xxy_1[i] * fe_0 + ta1_x_xxxzz_xxyz_0[i] * pa_z[i] - ta1_x_xxxzz_xxyz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxzz_0[i] = 2.0 * ta1_x_xxxz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxxzz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_xxz_1[i] * fe_0 + ta1_x_xxxzz_xxzz_0[i] * pa_z[i] - ta1_x_xxxzz_xxzz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xyyy_0[i] = 2.0 * ta1_x_xxxz_xyyy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xyyy_1[i] * fe_0 + ta1_x_xxxzz_xyyy_0[i] * pa_z[i] - ta1_x_xxxzz_xyyy_1[i] * pc_z[i];

        ta1_x_xxxzzz_xyyz_0[i] = 2.0 * ta1_x_xxxz_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xyyz_1[i] * fe_0 + ta1_x_xxxzz_xyy_0[i] * fe_0 - ta1_x_xxxzz_xyy_1[i] * fe_0 + ta1_x_xxxzz_xyyz_0[i] * pa_z[i] - ta1_x_xxxzz_xyyz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xyzz_0[i] = 2.0 * ta1_x_xxxz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xyzz_1[i] * fe_0 + 2.0 * ta1_x_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_xyz_1[i] * fe_0 + ta1_x_xxxzz_xyzz_0[i] * pa_z[i] - ta1_x_xxxzz_xyzz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xzzz_0[i] = 2.0 * ta1_x_xxxz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xzzz_1[i] * fe_0 + 3.0 * ta1_x_xxxzz_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxxzz_xzz_1[i] * fe_0 + ta1_x_xxxzz_xzzz_0[i] * pa_z[i] - ta1_x_xxxzz_xzzz_1[i] * pc_z[i];

        ta1_x_xxxzzz_yyyy_0[i] = 2.0 * ta1_x_xxxz_yyyy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_yyyy_1[i] * fe_0 + ta1_x_xxxzz_yyyy_0[i] * pa_z[i] - ta1_x_xxxzz_yyyy_1[i] * pc_z[i];

        ta1_x_xxxzzz_yyyz_0[i] = 2.0 * ta1_x_xzzz_yyyz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_yyyz_1[i] * fe_0 + ta_xxzzz_yyyz_1[i] + ta1_x_xxzzz_yyyz_0[i] * pa_x[i] - ta1_x_xxzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxxzzz_yyzz_0[i] = 2.0 * ta1_x_xzzz_yyzz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_yyzz_1[i] * fe_0 + ta_xxzzz_yyzz_1[i] + ta1_x_xxzzz_yyzz_0[i] * pa_x[i] - ta1_x_xxzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxxzzz_yzzz_0[i] = 2.0 * ta1_x_xzzz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_yzzz_1[i] * fe_0 + ta_xxzzz_yzzz_1[i] + ta1_x_xxzzz_yzzz_0[i] * pa_x[i] - ta1_x_xxzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxxzzz_zzzz_0[i] = 2.0 * ta1_x_xzzz_zzzz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_zzzz_1[i] * fe_0 + ta_xxzzz_zzzz_1[i] + ta1_x_xxzzz_zzzz_0[i] * pa_x[i] - ta1_x_xxzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : IG

    auto ta1_x_xxyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 150);

    auto ta1_x_xxyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 151);

    auto ta1_x_xxyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 152);

    auto ta1_x_xxyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 153);

    auto ta1_x_xxyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 154);

    auto ta1_x_xxyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 155);

    auto ta1_x_xxyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 156);

    auto ta1_x_xxyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 157);

    auto ta1_x_xxyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 158);

    auto ta1_x_xxyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 159);

    auto ta1_x_xxyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 160);

    auto ta1_x_xxyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 161);

    auto ta1_x_xxyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 162);

    auto ta1_x_xxyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 163);

    auto ta1_x_xxyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 164);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxyy_xxxx_0, ta1_x_xxyy_xxxx_1, ta1_x_xxyy_xxxy_0, ta1_x_xxyy_xxxy_1, ta1_x_xxyy_xxxz_0, ta1_x_xxyy_xxxz_1, ta1_x_xxyy_xxyy_0, ta1_x_xxyy_xxyy_1, ta1_x_xxyy_xxyz_0, ta1_x_xxyy_xxyz_1, ta1_x_xxyy_xxzz_0, ta1_x_xxyy_xxzz_1, ta1_x_xxyy_xyyy_0, ta1_x_xxyy_xyyy_1, ta1_x_xxyy_xyyz_0, ta1_x_xxyy_xyyz_1, ta1_x_xxyy_xyzz_0, ta1_x_xxyy_xyzz_1, ta1_x_xxyy_xzzz_0, ta1_x_xxyy_xzzz_1, ta1_x_xxyy_zzzz_0, ta1_x_xxyy_zzzz_1, ta1_x_xxyyy_xxx_0, ta1_x_xxyyy_xxx_1, ta1_x_xxyyy_xxxx_0, ta1_x_xxyyy_xxxx_1, ta1_x_xxyyy_xxxy_0, ta1_x_xxyyy_xxxy_1, ta1_x_xxyyy_xxxz_0, ta1_x_xxyyy_xxxz_1, ta1_x_xxyyy_xxy_0, ta1_x_xxyyy_xxy_1, ta1_x_xxyyy_xxyy_0, ta1_x_xxyyy_xxyy_1, ta1_x_xxyyy_xxyz_0, ta1_x_xxyyy_xxyz_1, ta1_x_xxyyy_xxz_0, ta1_x_xxyyy_xxz_1, ta1_x_xxyyy_xxzz_0, ta1_x_xxyyy_xxzz_1, ta1_x_xxyyy_xyy_0, ta1_x_xxyyy_xyy_1, ta1_x_xxyyy_xyyy_0, ta1_x_xxyyy_xyyy_1, ta1_x_xxyyy_xyyz_0, ta1_x_xxyyy_xyyz_1, ta1_x_xxyyy_xyz_0, ta1_x_xxyyy_xyz_1, ta1_x_xxyyy_xyzz_0, ta1_x_xxyyy_xyzz_1, ta1_x_xxyyy_xzz_0, ta1_x_xxyyy_xzz_1, ta1_x_xxyyy_xzzz_0, ta1_x_xxyyy_xzzz_1, ta1_x_xxyyy_zzzz_0, ta1_x_xxyyy_zzzz_1, ta1_x_xxyyyy_xxxx_0, ta1_x_xxyyyy_xxxy_0, ta1_x_xxyyyy_xxxz_0, ta1_x_xxyyyy_xxyy_0, ta1_x_xxyyyy_xxyz_0, ta1_x_xxyyyy_xxzz_0, ta1_x_xxyyyy_xyyy_0, ta1_x_xxyyyy_xyyz_0, ta1_x_xxyyyy_xyzz_0, ta1_x_xxyyyy_xzzz_0, ta1_x_xxyyyy_yyyy_0, ta1_x_xxyyyy_yyyz_0, ta1_x_xxyyyy_yyzz_0, ta1_x_xxyyyy_yzzz_0, ta1_x_xxyyyy_zzzz_0, ta1_x_xyyyy_yyyy_0, ta1_x_xyyyy_yyyy_1, ta1_x_xyyyy_yyyz_0, ta1_x_xyyyy_yyyz_1, ta1_x_xyyyy_yyzz_0, ta1_x_xyyyy_yyzz_1, ta1_x_xyyyy_yzzz_0, ta1_x_xyyyy_yzzz_1, ta1_x_yyyy_yyyy_0, ta1_x_yyyy_yyyy_1, ta1_x_yyyy_yyyz_0, ta1_x_yyyy_yyyz_1, ta1_x_yyyy_yyzz_0, ta1_x_yyyy_yyzz_1, ta1_x_yyyy_yzzz_0, ta1_x_yyyy_yzzz_1, ta_xyyyy_yyyy_1, ta_xyyyy_yyyz_1, ta_xyyyy_yyzz_1, ta_xyyyy_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyy_xxxx_0[i] = 3.0 * ta1_x_xxyy_xxxx_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxxx_1[i] * fe_0 + ta1_x_xxyyy_xxxx_0[i] * pa_y[i] - ta1_x_xxyyy_xxxx_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxxy_0[i] = 3.0 * ta1_x_xxyy_xxxy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxxy_1[i] * fe_0 + ta1_x_xxyyy_xxx_0[i] * fe_0 - ta1_x_xxyyy_xxx_1[i] * fe_0 + ta1_x_xxyyy_xxxy_0[i] * pa_y[i] - ta1_x_xxyyy_xxxy_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxxz_0[i] = 3.0 * ta1_x_xxyy_xxxz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxxz_1[i] * fe_0 + ta1_x_xxyyy_xxxz_0[i] * pa_y[i] - ta1_x_xxyyy_xxxz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxyy_0[i] = 3.0 * ta1_x_xxyy_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxyy_1[i] * fe_0 + 2.0 * ta1_x_xxyyy_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxyyy_xxy_1[i] * fe_0 + ta1_x_xxyyy_xxyy_0[i] * pa_y[i] - ta1_x_xxyyy_xxyy_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxyz_0[i] = 3.0 * ta1_x_xxyy_xxyz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxyz_1[i] * fe_0 + ta1_x_xxyyy_xxz_0[i] * fe_0 - ta1_x_xxyyy_xxz_1[i] * fe_0 + ta1_x_xxyyy_xxyz_0[i] * pa_y[i] - ta1_x_xxyyy_xxyz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxzz_0[i] = 3.0 * ta1_x_xxyy_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxzz_1[i] * fe_0 + ta1_x_xxyyy_xxzz_0[i] * pa_y[i] - ta1_x_xxyyy_xxzz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xyyy_0[i] = 3.0 * ta1_x_xxyy_xyyy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xyyy_1[i] * fe_0 + 3.0 * ta1_x_xxyyy_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxyyy_xyy_1[i] * fe_0 + ta1_x_xxyyy_xyyy_0[i] * pa_y[i] - ta1_x_xxyyy_xyyy_1[i] * pc_y[i];

        ta1_x_xxyyyy_xyyz_0[i] = 3.0 * ta1_x_xxyy_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxyyy_xyz_1[i] * fe_0 + ta1_x_xxyyy_xyyz_0[i] * pa_y[i] - ta1_x_xxyyy_xyyz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xyzz_0[i] = 3.0 * ta1_x_xxyy_xyzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xyzz_1[i] * fe_0 + ta1_x_xxyyy_xzz_0[i] * fe_0 - ta1_x_xxyyy_xzz_1[i] * fe_0 + ta1_x_xxyyy_xyzz_0[i] * pa_y[i] - ta1_x_xxyyy_xyzz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xzzz_0[i] = 3.0 * ta1_x_xxyy_xzzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xzzz_1[i] * fe_0 + ta1_x_xxyyy_xzzz_0[i] * pa_y[i] - ta1_x_xxyyy_xzzz_1[i] * pc_y[i];

        ta1_x_xxyyyy_yyyy_0[i] = ta1_x_yyyy_yyyy_0[i] * fe_0 - ta1_x_yyyy_yyyy_1[i] * fe_0 + ta_xyyyy_yyyy_1[i] + ta1_x_xyyyy_yyyy_0[i] * pa_x[i] - ta1_x_xyyyy_yyyy_1[i] * pc_x[i];

        ta1_x_xxyyyy_yyyz_0[i] = ta1_x_yyyy_yyyz_0[i] * fe_0 - ta1_x_yyyy_yyyz_1[i] * fe_0 + ta_xyyyy_yyyz_1[i] + ta1_x_xyyyy_yyyz_0[i] * pa_x[i] - ta1_x_xyyyy_yyyz_1[i] * pc_x[i];

        ta1_x_xxyyyy_yyzz_0[i] = ta1_x_yyyy_yyzz_0[i] * fe_0 - ta1_x_yyyy_yyzz_1[i] * fe_0 + ta_xyyyy_yyzz_1[i] + ta1_x_xyyyy_yyzz_0[i] * pa_x[i] - ta1_x_xyyyy_yyzz_1[i] * pc_x[i];

        ta1_x_xxyyyy_yzzz_0[i] = ta1_x_yyyy_yzzz_0[i] * fe_0 - ta1_x_yyyy_yzzz_1[i] * fe_0 + ta_xyyyy_yzzz_1[i] + ta1_x_xyyyy_yzzz_0[i] * pa_x[i] - ta1_x_xyyyy_yzzz_1[i] * pc_x[i];

        ta1_x_xxyyyy_zzzz_0[i] = 3.0 * ta1_x_xxyy_zzzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_zzzz_1[i] * fe_0 + ta1_x_xxyyy_zzzz_0[i] * pa_y[i] - ta1_x_xxyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : IG

    auto ta1_x_xxyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 165);

    auto ta1_x_xxyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 166);

    auto ta1_x_xxyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 167);

    auto ta1_x_xxyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 168);

    auto ta1_x_xxyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 169);

    auto ta1_x_xxyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 170);

    auto ta1_x_xxyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 171);

    auto ta1_x_xxyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 172);

    auto ta1_x_xxyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 173);

    auto ta1_x_xxyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 174);

    auto ta1_x_xxyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 175);

    auto ta1_x_xxyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 176);

    auto ta1_x_xxyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 177);

    auto ta1_x_xxyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 178);

    auto ta1_x_xxyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 179);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxyyy_xxxx_0, ta1_x_xxyyy_xxxx_1, ta1_x_xxyyy_xxxy_0, ta1_x_xxyyy_xxxy_1, ta1_x_xxyyy_xxy_0, ta1_x_xxyyy_xxy_1, ta1_x_xxyyy_xxyy_0, ta1_x_xxyyy_xxyy_1, ta1_x_xxyyy_xxyz_0, ta1_x_xxyyy_xxyz_1, ta1_x_xxyyy_xyy_0, ta1_x_xxyyy_xyy_1, ta1_x_xxyyy_xyyy_0, ta1_x_xxyyy_xyyy_1, ta1_x_xxyyy_xyyz_0, ta1_x_xxyyy_xyyz_1, ta1_x_xxyyy_xyz_0, ta1_x_xxyyy_xyz_1, ta1_x_xxyyy_xyzz_0, ta1_x_xxyyy_xyzz_1, ta1_x_xxyyy_yyy_0, ta1_x_xxyyy_yyy_1, ta1_x_xxyyy_yyyy_0, ta1_x_xxyyy_yyyy_1, ta1_x_xxyyy_yyyz_0, ta1_x_xxyyy_yyyz_1, ta1_x_xxyyy_yyz_0, ta1_x_xxyyy_yyz_1, ta1_x_xxyyy_yyzz_0, ta1_x_xxyyy_yyzz_1, ta1_x_xxyyy_yzz_0, ta1_x_xxyyy_yzz_1, ta1_x_xxyyy_yzzz_0, ta1_x_xxyyy_yzzz_1, ta1_x_xxyyyz_xxxx_0, ta1_x_xxyyyz_xxxy_0, ta1_x_xxyyyz_xxxz_0, ta1_x_xxyyyz_xxyy_0, ta1_x_xxyyyz_xxyz_0, ta1_x_xxyyyz_xxzz_0, ta1_x_xxyyyz_xyyy_0, ta1_x_xxyyyz_xyyz_0, ta1_x_xxyyyz_xyzz_0, ta1_x_xxyyyz_xzzz_0, ta1_x_xxyyyz_yyyy_0, ta1_x_xxyyyz_yyyz_0, ta1_x_xxyyyz_yyzz_0, ta1_x_xxyyyz_yzzz_0, ta1_x_xxyyyz_zzzz_0, ta1_x_xxyyz_xxxz_0, ta1_x_xxyyz_xxxz_1, ta1_x_xxyyz_xxzz_0, ta1_x_xxyyz_xxzz_1, ta1_x_xxyyz_xzzz_0, ta1_x_xxyyz_xzzz_1, ta1_x_xxyyz_zzzz_0, ta1_x_xxyyz_zzzz_1, ta1_x_xxyz_xxxz_0, ta1_x_xxyz_xxxz_1, ta1_x_xxyz_xxzz_0, ta1_x_xxyz_xxzz_1, ta1_x_xxyz_xzzz_0, ta1_x_xxyz_xzzz_1, ta1_x_xxyz_zzzz_0, ta1_x_xxyz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyz_xxxx_0[i] = ta1_x_xxyyy_xxxx_0[i] * pa_z[i] - ta1_x_xxyyy_xxxx_1[i] * pc_z[i];

        ta1_x_xxyyyz_xxxy_0[i] = ta1_x_xxyyy_xxxy_0[i] * pa_z[i] - ta1_x_xxyyy_xxxy_1[i] * pc_z[i];

        ta1_x_xxyyyz_xxxz_0[i] = 2.0 * ta1_x_xxyz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_xxxz_1[i] * fe_0 + ta1_x_xxyyz_xxxz_0[i] * pa_y[i] - ta1_x_xxyyz_xxxz_1[i] * pc_y[i];

        ta1_x_xxyyyz_xxyy_0[i] = ta1_x_xxyyy_xxyy_0[i] * pa_z[i] - ta1_x_xxyyy_xxyy_1[i] * pc_z[i];

        ta1_x_xxyyyz_xxyz_0[i] = ta1_x_xxyyy_xxy_0[i] * fe_0 - ta1_x_xxyyy_xxy_1[i] * fe_0 + ta1_x_xxyyy_xxyz_0[i] * pa_z[i] - ta1_x_xxyyy_xxyz_1[i] * pc_z[i];

        ta1_x_xxyyyz_xxzz_0[i] = 2.0 * ta1_x_xxyz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_xxzz_1[i] * fe_0 + ta1_x_xxyyz_xxzz_0[i] * pa_y[i] - ta1_x_xxyyz_xxzz_1[i] * pc_y[i];

        ta1_x_xxyyyz_xyyy_0[i] = ta1_x_xxyyy_xyyy_0[i] * pa_z[i] - ta1_x_xxyyy_xyyy_1[i] * pc_z[i];

        ta1_x_xxyyyz_xyyz_0[i] = ta1_x_xxyyy_xyy_0[i] * fe_0 - ta1_x_xxyyy_xyy_1[i] * fe_0 + ta1_x_xxyyy_xyyz_0[i] * pa_z[i] - ta1_x_xxyyy_xyyz_1[i] * pc_z[i];

        ta1_x_xxyyyz_xyzz_0[i] = 2.0 * ta1_x_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxyyy_xyz_1[i] * fe_0 + ta1_x_xxyyy_xyzz_0[i] * pa_z[i] - ta1_x_xxyyy_xyzz_1[i] * pc_z[i];

        ta1_x_xxyyyz_xzzz_0[i] = 2.0 * ta1_x_xxyz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_xzzz_1[i] * fe_0 + ta1_x_xxyyz_xzzz_0[i] * pa_y[i] - ta1_x_xxyyz_xzzz_1[i] * pc_y[i];

        ta1_x_xxyyyz_yyyy_0[i] = ta1_x_xxyyy_yyyy_0[i] * pa_z[i] - ta1_x_xxyyy_yyyy_1[i] * pc_z[i];

        ta1_x_xxyyyz_yyyz_0[i] = ta1_x_xxyyy_yyy_0[i] * fe_0 - ta1_x_xxyyy_yyy_1[i] * fe_0 + ta1_x_xxyyy_yyyz_0[i] * pa_z[i] - ta1_x_xxyyy_yyyz_1[i] * pc_z[i];

        ta1_x_xxyyyz_yyzz_0[i] = 2.0 * ta1_x_xxyyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_xxyyy_yyz_1[i] * fe_0 + ta1_x_xxyyy_yyzz_0[i] * pa_z[i] - ta1_x_xxyyy_yyzz_1[i] * pc_z[i];

        ta1_x_xxyyyz_yzzz_0[i] = 3.0 * ta1_x_xxyyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxyyy_yzz_1[i] * fe_0 + ta1_x_xxyyy_yzzz_0[i] * pa_z[i] - ta1_x_xxyyy_yzzz_1[i] * pc_z[i];

        ta1_x_xxyyyz_zzzz_0[i] = 2.0 * ta1_x_xxyz_zzzz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_zzzz_1[i] * fe_0 + ta1_x_xxyyz_zzzz_0[i] * pa_y[i] - ta1_x_xxyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 180-195 components of targeted buffer : IG

    auto ta1_x_xxyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 180);

    auto ta1_x_xxyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 181);

    auto ta1_x_xxyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 182);

    auto ta1_x_xxyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 183);

    auto ta1_x_xxyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 184);

    auto ta1_x_xxyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 185);

    auto ta1_x_xxyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 186);

    auto ta1_x_xxyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 187);

    auto ta1_x_xxyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 188);

    auto ta1_x_xxyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 189);

    auto ta1_x_xxyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 190);

    auto ta1_x_xxyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 191);

    auto ta1_x_xxyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 192);

    auto ta1_x_xxyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 193);

    auto ta1_x_xxyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 194);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xxyy_xxxy_0, ta1_x_xxyy_xxxy_1, ta1_x_xxyy_xxyy_0, ta1_x_xxyy_xxyy_1, ta1_x_xxyy_xyyy_0, ta1_x_xxyy_xyyy_1, ta1_x_xxyy_yyyy_0, ta1_x_xxyy_yyyy_1, ta1_x_xxyyz_xxxy_0, ta1_x_xxyyz_xxxy_1, ta1_x_xxyyz_xxyy_0, ta1_x_xxyyz_xxyy_1, ta1_x_xxyyz_xyyy_0, ta1_x_xxyyz_xyyy_1, ta1_x_xxyyz_yyyy_0, ta1_x_xxyyz_yyyy_1, ta1_x_xxyyzz_xxxx_0, ta1_x_xxyyzz_xxxy_0, ta1_x_xxyyzz_xxxz_0, ta1_x_xxyyzz_xxyy_0, ta1_x_xxyyzz_xxyz_0, ta1_x_xxyyzz_xxzz_0, ta1_x_xxyyzz_xyyy_0, ta1_x_xxyyzz_xyyz_0, ta1_x_xxyyzz_xyzz_0, ta1_x_xxyyzz_xzzz_0, ta1_x_xxyyzz_yyyy_0, ta1_x_xxyyzz_yyyz_0, ta1_x_xxyyzz_yyzz_0, ta1_x_xxyyzz_yzzz_0, ta1_x_xxyyzz_zzzz_0, ta1_x_xxyzz_xxxx_0, ta1_x_xxyzz_xxxx_1, ta1_x_xxyzz_xxxz_0, ta1_x_xxyzz_xxxz_1, ta1_x_xxyzz_xxyz_0, ta1_x_xxyzz_xxyz_1, ta1_x_xxyzz_xxz_0, ta1_x_xxyzz_xxz_1, ta1_x_xxyzz_xxzz_0, ta1_x_xxyzz_xxzz_1, ta1_x_xxyzz_xyyz_0, ta1_x_xxyzz_xyyz_1, ta1_x_xxyzz_xyz_0, ta1_x_xxyzz_xyz_1, ta1_x_xxyzz_xyzz_0, ta1_x_xxyzz_xyzz_1, ta1_x_xxyzz_xzz_0, ta1_x_xxyzz_xzz_1, ta1_x_xxyzz_xzzz_0, ta1_x_xxyzz_xzzz_1, ta1_x_xxyzz_zzzz_0, ta1_x_xxyzz_zzzz_1, ta1_x_xxzz_xxxx_0, ta1_x_xxzz_xxxx_1, ta1_x_xxzz_xxxz_0, ta1_x_xxzz_xxxz_1, ta1_x_xxzz_xxyz_0, ta1_x_xxzz_xxyz_1, ta1_x_xxzz_xxzz_0, ta1_x_xxzz_xxzz_1, ta1_x_xxzz_xyyz_0, ta1_x_xxzz_xyyz_1, ta1_x_xxzz_xyzz_0, ta1_x_xxzz_xyzz_1, ta1_x_xxzz_xzzz_0, ta1_x_xxzz_xzzz_1, ta1_x_xxzz_zzzz_0, ta1_x_xxzz_zzzz_1, ta1_x_xyyzz_yyyz_0, ta1_x_xyyzz_yyyz_1, ta1_x_xyyzz_yyzz_0, ta1_x_xyyzz_yyzz_1, ta1_x_xyyzz_yzzz_0, ta1_x_xyyzz_yzzz_1, ta1_x_yyzz_yyyz_0, ta1_x_yyzz_yyyz_1, ta1_x_yyzz_yyzz_0, ta1_x_yyzz_yyzz_1, ta1_x_yyzz_yzzz_0, ta1_x_yyzz_yzzz_1, ta_xyyzz_yyyz_1, ta_xyyzz_yyzz_1, ta_xyyzz_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyzz_xxxx_0[i] = ta1_x_xxzz_xxxx_0[i] * fe_0 - ta1_x_xxzz_xxxx_1[i] * fe_0 + ta1_x_xxyzz_xxxx_0[i] * pa_y[i] - ta1_x_xxyzz_xxxx_1[i] * pc_y[i];

        ta1_x_xxyyzz_xxxy_0[i] = ta1_x_xxyy_xxxy_0[i] * fe_0 - ta1_x_xxyy_xxxy_1[i] * fe_0 + ta1_x_xxyyz_xxxy_0[i] * pa_z[i] - ta1_x_xxyyz_xxxy_1[i] * pc_z[i];

        ta1_x_xxyyzz_xxxz_0[i] = ta1_x_xxzz_xxxz_0[i] * fe_0 - ta1_x_xxzz_xxxz_1[i] * fe_0 + ta1_x_xxyzz_xxxz_0[i] * pa_y[i] - ta1_x_xxyzz_xxxz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xxyy_0[i] = ta1_x_xxyy_xxyy_0[i] * fe_0 - ta1_x_xxyy_xxyy_1[i] * fe_0 + ta1_x_xxyyz_xxyy_0[i] * pa_z[i] - ta1_x_xxyyz_xxyy_1[i] * pc_z[i];

        ta1_x_xxyyzz_xxyz_0[i] = ta1_x_xxzz_xxyz_0[i] * fe_0 - ta1_x_xxzz_xxyz_1[i] * fe_0 + ta1_x_xxyzz_xxz_0[i] * fe_0 - ta1_x_xxyzz_xxz_1[i] * fe_0 + ta1_x_xxyzz_xxyz_0[i] * pa_y[i] - ta1_x_xxyzz_xxyz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xxzz_0[i] = ta1_x_xxzz_xxzz_0[i] * fe_0 - ta1_x_xxzz_xxzz_1[i] * fe_0 + ta1_x_xxyzz_xxzz_0[i] * pa_y[i] - ta1_x_xxyzz_xxzz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xyyy_0[i] = ta1_x_xxyy_xyyy_0[i] * fe_0 - ta1_x_xxyy_xyyy_1[i] * fe_0 + ta1_x_xxyyz_xyyy_0[i] * pa_z[i] - ta1_x_xxyyz_xyyy_1[i] * pc_z[i];

        ta1_x_xxyyzz_xyyz_0[i] = ta1_x_xxzz_xyyz_0[i] * fe_0 - ta1_x_xxzz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_xxyzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxyzz_xyz_1[i] * fe_0 + ta1_x_xxyzz_xyyz_0[i] * pa_y[i] - ta1_x_xxyzz_xyyz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xyzz_0[i] = ta1_x_xxzz_xyzz_0[i] * fe_0 - ta1_x_xxzz_xyzz_1[i] * fe_0 + ta1_x_xxyzz_xzz_0[i] * fe_0 - ta1_x_xxyzz_xzz_1[i] * fe_0 + ta1_x_xxyzz_xyzz_0[i] * pa_y[i] - ta1_x_xxyzz_xyzz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xzzz_0[i] = ta1_x_xxzz_xzzz_0[i] * fe_0 - ta1_x_xxzz_xzzz_1[i] * fe_0 + ta1_x_xxyzz_xzzz_0[i] * pa_y[i] - ta1_x_xxyzz_xzzz_1[i] * pc_y[i];

        ta1_x_xxyyzz_yyyy_0[i] = ta1_x_xxyy_yyyy_0[i] * fe_0 - ta1_x_xxyy_yyyy_1[i] * fe_0 + ta1_x_xxyyz_yyyy_0[i] * pa_z[i] - ta1_x_xxyyz_yyyy_1[i] * pc_z[i];

        ta1_x_xxyyzz_yyyz_0[i] = ta1_x_yyzz_yyyz_0[i] * fe_0 - ta1_x_yyzz_yyyz_1[i] * fe_0 + ta_xyyzz_yyyz_1[i] + ta1_x_xyyzz_yyyz_0[i] * pa_x[i] - ta1_x_xyyzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxyyzz_yyzz_0[i] = ta1_x_yyzz_yyzz_0[i] * fe_0 - ta1_x_yyzz_yyzz_1[i] * fe_0 + ta_xyyzz_yyzz_1[i] + ta1_x_xyyzz_yyzz_0[i] * pa_x[i] - ta1_x_xyyzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxyyzz_yzzz_0[i] = ta1_x_yyzz_yzzz_0[i] * fe_0 - ta1_x_yyzz_yzzz_1[i] * fe_0 + ta_xyyzz_yzzz_1[i] + ta1_x_xyyzz_yzzz_0[i] * pa_x[i] - ta1_x_xyyzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxyyzz_zzzz_0[i] = ta1_x_xxzz_zzzz_0[i] * fe_0 - ta1_x_xxzz_zzzz_1[i] * fe_0 + ta1_x_xxyzz_zzzz_0[i] * pa_y[i] - ta1_x_xxyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 195-210 components of targeted buffer : IG

    auto ta1_x_xxyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 195);

    auto ta1_x_xxyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 196);

    auto ta1_x_xxyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 197);

    auto ta1_x_xxyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 198);

    auto ta1_x_xxyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 199);

    auto ta1_x_xxyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 200);

    auto ta1_x_xxyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 201);

    auto ta1_x_xxyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 202);

    auto ta1_x_xxyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 203);

    auto ta1_x_xxyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 204);

    auto ta1_x_xxyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 205);

    auto ta1_x_xxyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 206);

    auto ta1_x_xxyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 207);

    auto ta1_x_xxyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 208);

    auto ta1_x_xxyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 209);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxyzzz_xxxx_0, ta1_x_xxyzzz_xxxy_0, ta1_x_xxyzzz_xxxz_0, ta1_x_xxyzzz_xxyy_0, ta1_x_xxyzzz_xxyz_0, ta1_x_xxyzzz_xxzz_0, ta1_x_xxyzzz_xyyy_0, ta1_x_xxyzzz_xyyz_0, ta1_x_xxyzzz_xyzz_0, ta1_x_xxyzzz_xzzz_0, ta1_x_xxyzzz_yyyy_0, ta1_x_xxyzzz_yyyz_0, ta1_x_xxyzzz_yyzz_0, ta1_x_xxyzzz_yzzz_0, ta1_x_xxyzzz_zzzz_0, ta1_x_xxzzz_xxx_0, ta1_x_xxzzz_xxx_1, ta1_x_xxzzz_xxxx_0, ta1_x_xxzzz_xxxx_1, ta1_x_xxzzz_xxxy_0, ta1_x_xxzzz_xxxy_1, ta1_x_xxzzz_xxxz_0, ta1_x_xxzzz_xxxz_1, ta1_x_xxzzz_xxy_0, ta1_x_xxzzz_xxy_1, ta1_x_xxzzz_xxyy_0, ta1_x_xxzzz_xxyy_1, ta1_x_xxzzz_xxyz_0, ta1_x_xxzzz_xxyz_1, ta1_x_xxzzz_xxz_0, ta1_x_xxzzz_xxz_1, ta1_x_xxzzz_xxzz_0, ta1_x_xxzzz_xxzz_1, ta1_x_xxzzz_xyy_0, ta1_x_xxzzz_xyy_1, ta1_x_xxzzz_xyyy_0, ta1_x_xxzzz_xyyy_1, ta1_x_xxzzz_xyyz_0, ta1_x_xxzzz_xyyz_1, ta1_x_xxzzz_xyz_0, ta1_x_xxzzz_xyz_1, ta1_x_xxzzz_xyzz_0, ta1_x_xxzzz_xyzz_1, ta1_x_xxzzz_xzz_0, ta1_x_xxzzz_xzz_1, ta1_x_xxzzz_xzzz_0, ta1_x_xxzzz_xzzz_1, ta1_x_xxzzz_yyy_0, ta1_x_xxzzz_yyy_1, ta1_x_xxzzz_yyyy_0, ta1_x_xxzzz_yyyy_1, ta1_x_xxzzz_yyyz_0, ta1_x_xxzzz_yyyz_1, ta1_x_xxzzz_yyz_0, ta1_x_xxzzz_yyz_1, ta1_x_xxzzz_yyzz_0, ta1_x_xxzzz_yyzz_1, ta1_x_xxzzz_yzz_0, ta1_x_xxzzz_yzz_1, ta1_x_xxzzz_yzzz_0, ta1_x_xxzzz_yzzz_1, ta1_x_xxzzz_zzz_0, ta1_x_xxzzz_zzz_1, ta1_x_xxzzz_zzzz_0, ta1_x_xxzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzzz_xxxx_0[i] = ta1_x_xxzzz_xxxx_0[i] * pa_y[i] - ta1_x_xxzzz_xxxx_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxxy_0[i] = ta1_x_xxzzz_xxx_0[i] * fe_0 - ta1_x_xxzzz_xxx_1[i] * fe_0 + ta1_x_xxzzz_xxxy_0[i] * pa_y[i] - ta1_x_xxzzz_xxxy_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxxz_0[i] = ta1_x_xxzzz_xxxz_0[i] * pa_y[i] - ta1_x_xxzzz_xxxz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxyy_0[i] = 2.0 * ta1_x_xxzzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_xxy_1[i] * fe_0 + ta1_x_xxzzz_xxyy_0[i] * pa_y[i] - ta1_x_xxzzz_xxyy_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxyz_0[i] = ta1_x_xxzzz_xxz_0[i] * fe_0 - ta1_x_xxzzz_xxz_1[i] * fe_0 + ta1_x_xxzzz_xxyz_0[i] * pa_y[i] - ta1_x_xxzzz_xxyz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxzz_0[i] = ta1_x_xxzzz_xxzz_0[i] * pa_y[i] - ta1_x_xxzzz_xxzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xyyy_0[i] = 3.0 * ta1_x_xxzzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxzzz_xyy_1[i] * fe_0 + ta1_x_xxzzz_xyyy_0[i] * pa_y[i] - ta1_x_xxzzz_xyyy_1[i] * pc_y[i];

        ta1_x_xxyzzz_xyyz_0[i] = 2.0 * ta1_x_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_xyz_1[i] * fe_0 + ta1_x_xxzzz_xyyz_0[i] * pa_y[i] - ta1_x_xxzzz_xyyz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xyzz_0[i] = ta1_x_xxzzz_xzz_0[i] * fe_0 - ta1_x_xxzzz_xzz_1[i] * fe_0 + ta1_x_xxzzz_xyzz_0[i] * pa_y[i] - ta1_x_xxzzz_xyzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xzzz_0[i] = ta1_x_xxzzz_xzzz_0[i] * pa_y[i] - ta1_x_xxzzz_xzzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_yyyy_0[i] = 4.0 * ta1_x_xxzzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxzzz_yyy_1[i] * fe_0 + ta1_x_xxzzz_yyyy_0[i] * pa_y[i] - ta1_x_xxzzz_yyyy_1[i] * pc_y[i];

        ta1_x_xxyzzz_yyyz_0[i] = 3.0 * ta1_x_xxzzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxzzz_yyz_1[i] * fe_0 + ta1_x_xxzzz_yyyz_0[i] * pa_y[i] - ta1_x_xxzzz_yyyz_1[i] * pc_y[i];

        ta1_x_xxyzzz_yyzz_0[i] = 2.0 * ta1_x_xxzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_yzz_1[i] * fe_0 + ta1_x_xxzzz_yyzz_0[i] * pa_y[i] - ta1_x_xxzzz_yyzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_yzzz_0[i] = ta1_x_xxzzz_zzz_0[i] * fe_0 - ta1_x_xxzzz_zzz_1[i] * fe_0 + ta1_x_xxzzz_yzzz_0[i] * pa_y[i] - ta1_x_xxzzz_yzzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_zzzz_0[i] = ta1_x_xxzzz_zzzz_0[i] * pa_y[i] - ta1_x_xxzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : IG

    auto ta1_x_xxzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 210);

    auto ta1_x_xxzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 211);

    auto ta1_x_xxzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 212);

    auto ta1_x_xxzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 213);

    auto ta1_x_xxzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 214);

    auto ta1_x_xxzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 215);

    auto ta1_x_xxzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 216);

    auto ta1_x_xxzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 217);

    auto ta1_x_xxzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 218);

    auto ta1_x_xxzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 219);

    auto ta1_x_xxzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 220);

    auto ta1_x_xxzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 221);

    auto ta1_x_xxzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 222);

    auto ta1_x_xxzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 223);

    auto ta1_x_xxzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 224);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxzz_xxxx_0, ta1_x_xxzz_xxxx_1, ta1_x_xxzz_xxxy_0, ta1_x_xxzz_xxxy_1, ta1_x_xxzz_xxxz_0, ta1_x_xxzz_xxxz_1, ta1_x_xxzz_xxyy_0, ta1_x_xxzz_xxyy_1, ta1_x_xxzz_xxyz_0, ta1_x_xxzz_xxyz_1, ta1_x_xxzz_xxzz_0, ta1_x_xxzz_xxzz_1, ta1_x_xxzz_xyyy_0, ta1_x_xxzz_xyyy_1, ta1_x_xxzz_xyyz_0, ta1_x_xxzz_xyyz_1, ta1_x_xxzz_xyzz_0, ta1_x_xxzz_xyzz_1, ta1_x_xxzz_xzzz_0, ta1_x_xxzz_xzzz_1, ta1_x_xxzz_yyyy_0, ta1_x_xxzz_yyyy_1, ta1_x_xxzzz_xxx_0, ta1_x_xxzzz_xxx_1, ta1_x_xxzzz_xxxx_0, ta1_x_xxzzz_xxxx_1, ta1_x_xxzzz_xxxy_0, ta1_x_xxzzz_xxxy_1, ta1_x_xxzzz_xxxz_0, ta1_x_xxzzz_xxxz_1, ta1_x_xxzzz_xxy_0, ta1_x_xxzzz_xxy_1, ta1_x_xxzzz_xxyy_0, ta1_x_xxzzz_xxyy_1, ta1_x_xxzzz_xxyz_0, ta1_x_xxzzz_xxyz_1, ta1_x_xxzzz_xxz_0, ta1_x_xxzzz_xxz_1, ta1_x_xxzzz_xxzz_0, ta1_x_xxzzz_xxzz_1, ta1_x_xxzzz_xyy_0, ta1_x_xxzzz_xyy_1, ta1_x_xxzzz_xyyy_0, ta1_x_xxzzz_xyyy_1, ta1_x_xxzzz_xyyz_0, ta1_x_xxzzz_xyyz_1, ta1_x_xxzzz_xyz_0, ta1_x_xxzzz_xyz_1, ta1_x_xxzzz_xyzz_0, ta1_x_xxzzz_xyzz_1, ta1_x_xxzzz_xzz_0, ta1_x_xxzzz_xzz_1, ta1_x_xxzzz_xzzz_0, ta1_x_xxzzz_xzzz_1, ta1_x_xxzzz_yyyy_0, ta1_x_xxzzz_yyyy_1, ta1_x_xxzzzz_xxxx_0, ta1_x_xxzzzz_xxxy_0, ta1_x_xxzzzz_xxxz_0, ta1_x_xxzzzz_xxyy_0, ta1_x_xxzzzz_xxyz_0, ta1_x_xxzzzz_xxzz_0, ta1_x_xxzzzz_xyyy_0, ta1_x_xxzzzz_xyyz_0, ta1_x_xxzzzz_xyzz_0, ta1_x_xxzzzz_xzzz_0, ta1_x_xxzzzz_yyyy_0, ta1_x_xxzzzz_yyyz_0, ta1_x_xxzzzz_yyzz_0, ta1_x_xxzzzz_yzzz_0, ta1_x_xxzzzz_zzzz_0, ta1_x_xzzzz_yyyz_0, ta1_x_xzzzz_yyyz_1, ta1_x_xzzzz_yyzz_0, ta1_x_xzzzz_yyzz_1, ta1_x_xzzzz_yzzz_0, ta1_x_xzzzz_yzzz_1, ta1_x_xzzzz_zzzz_0, ta1_x_xzzzz_zzzz_1, ta1_x_zzzz_yyyz_0, ta1_x_zzzz_yyyz_1, ta1_x_zzzz_yyzz_0, ta1_x_zzzz_yyzz_1, ta1_x_zzzz_yzzz_0, ta1_x_zzzz_yzzz_1, ta1_x_zzzz_zzzz_0, ta1_x_zzzz_zzzz_1, ta_xzzzz_yyyz_1, ta_xzzzz_yyzz_1, ta_xzzzz_yzzz_1, ta_xzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzzz_xxxx_0[i] = 3.0 * ta1_x_xxzz_xxxx_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxxx_1[i] * fe_0 + ta1_x_xxzzz_xxxx_0[i] * pa_z[i] - ta1_x_xxzzz_xxxx_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxxy_0[i] = 3.0 * ta1_x_xxzz_xxxy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxxy_1[i] * fe_0 + ta1_x_xxzzz_xxxy_0[i] * pa_z[i] - ta1_x_xxzzz_xxxy_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxxz_0[i] = 3.0 * ta1_x_xxzz_xxxz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxxz_1[i] * fe_0 + ta1_x_xxzzz_xxx_0[i] * fe_0 - ta1_x_xxzzz_xxx_1[i] * fe_0 + ta1_x_xxzzz_xxxz_0[i] * pa_z[i] - ta1_x_xxzzz_xxxz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxyy_0[i] = 3.0 * ta1_x_xxzz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxyy_1[i] * fe_0 + ta1_x_xxzzz_xxyy_0[i] * pa_z[i] - ta1_x_xxzzz_xxyy_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxyz_0[i] = 3.0 * ta1_x_xxzz_xxyz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxyz_1[i] * fe_0 + ta1_x_xxzzz_xxy_0[i] * fe_0 - ta1_x_xxzzz_xxy_1[i] * fe_0 + ta1_x_xxzzz_xxyz_0[i] * pa_z[i] - ta1_x_xxzzz_xxyz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxzz_0[i] = 3.0 * ta1_x_xxzz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxzz_1[i] * fe_0 + 2.0 * ta1_x_xxzzz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_xxz_1[i] * fe_0 + ta1_x_xxzzz_xxzz_0[i] * pa_z[i] - ta1_x_xxzzz_xxzz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xyyy_0[i] = 3.0 * ta1_x_xxzz_xyyy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xyyy_1[i] * fe_0 + ta1_x_xxzzz_xyyy_0[i] * pa_z[i] - ta1_x_xxzzz_xyyy_1[i] * pc_z[i];

        ta1_x_xxzzzz_xyyz_0[i] = 3.0 * ta1_x_xxzz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xyyz_1[i] * fe_0 + ta1_x_xxzzz_xyy_0[i] * fe_0 - ta1_x_xxzzz_xyy_1[i] * fe_0 + ta1_x_xxzzz_xyyz_0[i] * pa_z[i] - ta1_x_xxzzz_xyyz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xyzz_0[i] = 3.0 * ta1_x_xxzz_xyzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xyzz_1[i] * fe_0 + 2.0 * ta1_x_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_xyz_1[i] * fe_0 + ta1_x_xxzzz_xyzz_0[i] * pa_z[i] - ta1_x_xxzzz_xyzz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xzzz_0[i] = 3.0 * ta1_x_xxzz_xzzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xzzz_1[i] * fe_0 + 3.0 * ta1_x_xxzzz_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxzzz_xzz_1[i] * fe_0 + ta1_x_xxzzz_xzzz_0[i] * pa_z[i] - ta1_x_xxzzz_xzzz_1[i] * pc_z[i];

        ta1_x_xxzzzz_yyyy_0[i] = 3.0 * ta1_x_xxzz_yyyy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yyyy_1[i] * fe_0 + ta1_x_xxzzz_yyyy_0[i] * pa_z[i] - ta1_x_xxzzz_yyyy_1[i] * pc_z[i];

        ta1_x_xxzzzz_yyyz_0[i] = ta1_x_zzzz_yyyz_0[i] * fe_0 - ta1_x_zzzz_yyyz_1[i] * fe_0 + ta_xzzzz_yyyz_1[i] + ta1_x_xzzzz_yyyz_0[i] * pa_x[i] - ta1_x_xzzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xxzzzz_yyzz_0[i] = ta1_x_zzzz_yyzz_0[i] * fe_0 - ta1_x_zzzz_yyzz_1[i] * fe_0 + ta_xzzzz_yyzz_1[i] + ta1_x_xzzzz_yyzz_0[i] * pa_x[i] - ta1_x_xzzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xxzzzz_yzzz_0[i] = ta1_x_zzzz_yzzz_0[i] * fe_0 - ta1_x_zzzz_yzzz_1[i] * fe_0 + ta_xzzzz_yzzz_1[i] + ta1_x_xzzzz_yzzz_0[i] * pa_x[i] - ta1_x_xzzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xxzzzz_zzzz_0[i] = ta1_x_zzzz_zzzz_0[i] * fe_0 - ta1_x_zzzz_zzzz_1[i] * fe_0 + ta_xzzzz_zzzz_1[i] + ta1_x_xzzzz_zzzz_0[i] * pa_x[i] - ta1_x_xzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : IG

    auto ta1_x_xyyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 225);

    auto ta1_x_xyyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 226);

    auto ta1_x_xyyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 227);

    auto ta1_x_xyyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 228);

    auto ta1_x_xyyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 229);

    auto ta1_x_xyyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 230);

    auto ta1_x_xyyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 231);

    auto ta1_x_xyyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 232);

    auto ta1_x_xyyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 233);

    auto ta1_x_xyyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 234);

    auto ta1_x_xyyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 235);

    auto ta1_x_xyyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 236);

    auto ta1_x_xyyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 237);

    auto ta1_x_xyyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 238);

    auto ta1_x_xyyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 239);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyyy_xxxx_0, ta1_x_xyyy_xxxx_1, ta1_x_xyyy_xxxz_0, ta1_x_xyyy_xxxz_1, ta1_x_xyyy_xxzz_0, ta1_x_xyyy_xxzz_1, ta1_x_xyyy_xzzz_0, ta1_x_xyyy_xzzz_1, ta1_x_xyyyy_xxxx_0, ta1_x_xyyyy_xxxx_1, ta1_x_xyyyy_xxxz_0, ta1_x_xyyyy_xxxz_1, ta1_x_xyyyy_xxzz_0, ta1_x_xyyyy_xxzz_1, ta1_x_xyyyy_xzzz_0, ta1_x_xyyyy_xzzz_1, ta1_x_xyyyyy_xxxx_0, ta1_x_xyyyyy_xxxy_0, ta1_x_xyyyyy_xxxz_0, ta1_x_xyyyyy_xxyy_0, ta1_x_xyyyyy_xxyz_0, ta1_x_xyyyyy_xxzz_0, ta1_x_xyyyyy_xyyy_0, ta1_x_xyyyyy_xyyz_0, ta1_x_xyyyyy_xyzz_0, ta1_x_xyyyyy_xzzz_0, ta1_x_xyyyyy_yyyy_0, ta1_x_xyyyyy_yyyz_0, ta1_x_xyyyyy_yyzz_0, ta1_x_xyyyyy_yzzz_0, ta1_x_xyyyyy_zzzz_0, ta1_x_yyyyy_xxxy_0, ta1_x_yyyyy_xxxy_1, ta1_x_yyyyy_xxy_0, ta1_x_yyyyy_xxy_1, ta1_x_yyyyy_xxyy_0, ta1_x_yyyyy_xxyy_1, ta1_x_yyyyy_xxyz_0, ta1_x_yyyyy_xxyz_1, ta1_x_yyyyy_xyy_0, ta1_x_yyyyy_xyy_1, ta1_x_yyyyy_xyyy_0, ta1_x_yyyyy_xyyy_1, ta1_x_yyyyy_xyyz_0, ta1_x_yyyyy_xyyz_1, ta1_x_yyyyy_xyz_0, ta1_x_yyyyy_xyz_1, ta1_x_yyyyy_xyzz_0, ta1_x_yyyyy_xyzz_1, ta1_x_yyyyy_yyy_0, ta1_x_yyyyy_yyy_1, ta1_x_yyyyy_yyyy_0, ta1_x_yyyyy_yyyy_1, ta1_x_yyyyy_yyyz_0, ta1_x_yyyyy_yyyz_1, ta1_x_yyyyy_yyz_0, ta1_x_yyyyy_yyz_1, ta1_x_yyyyy_yyzz_0, ta1_x_yyyyy_yyzz_1, ta1_x_yyyyy_yzz_0, ta1_x_yyyyy_yzz_1, ta1_x_yyyyy_yzzz_0, ta1_x_yyyyy_yzzz_1, ta1_x_yyyyy_zzzz_0, ta1_x_yyyyy_zzzz_1, ta_yyyyy_xxxy_1, ta_yyyyy_xxyy_1, ta_yyyyy_xxyz_1, ta_yyyyy_xyyy_1, ta_yyyyy_xyyz_1, ta_yyyyy_xyzz_1, ta_yyyyy_yyyy_1, ta_yyyyy_yyyz_1, ta_yyyyy_yyzz_1, ta_yyyyy_yzzz_1, ta_yyyyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyy_xxxx_0[i] = 4.0 * ta1_x_xyyy_xxxx_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xxxx_1[i] * fe_0 + ta1_x_xyyyy_xxxx_0[i] * pa_y[i] - ta1_x_xyyyy_xxxx_1[i] * pc_y[i];

        ta1_x_xyyyyy_xxxy_0[i] = 3.0 * ta1_x_yyyyy_xxy_0[i] * fe_0 - 3.0 * ta1_x_yyyyy_xxy_1[i] * fe_0 + ta_yyyyy_xxxy_1[i] + ta1_x_yyyyy_xxxy_0[i] * pa_x[i] - ta1_x_yyyyy_xxxy_1[i] * pc_x[i];

        ta1_x_xyyyyy_xxxz_0[i] = 4.0 * ta1_x_xyyy_xxxz_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xxxz_1[i] * fe_0 + ta1_x_xyyyy_xxxz_0[i] * pa_y[i] - ta1_x_xyyyy_xxxz_1[i] * pc_y[i];

        ta1_x_xyyyyy_xxyy_0[i] = 2.0 * ta1_x_yyyyy_xyy_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_xyy_1[i] * fe_0 + ta_yyyyy_xxyy_1[i] + ta1_x_yyyyy_xxyy_0[i] * pa_x[i] - ta1_x_yyyyy_xxyy_1[i] * pc_x[i];

        ta1_x_xyyyyy_xxyz_0[i] = 2.0 * ta1_x_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_xyz_1[i] * fe_0 + ta_yyyyy_xxyz_1[i] + ta1_x_yyyyy_xxyz_0[i] * pa_x[i] - ta1_x_yyyyy_xxyz_1[i] * pc_x[i];

        ta1_x_xyyyyy_xxzz_0[i] = 4.0 * ta1_x_xyyy_xxzz_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xxzz_1[i] * fe_0 + ta1_x_xyyyy_xxzz_0[i] * pa_y[i] - ta1_x_xyyyy_xxzz_1[i] * pc_y[i];

        ta1_x_xyyyyy_xyyy_0[i] = ta1_x_yyyyy_yyy_0[i] * fe_0 - ta1_x_yyyyy_yyy_1[i] * fe_0 + ta_yyyyy_xyyy_1[i] + ta1_x_yyyyy_xyyy_0[i] * pa_x[i] - ta1_x_yyyyy_xyyy_1[i] * pc_x[i];

        ta1_x_xyyyyy_xyyz_0[i] = ta1_x_yyyyy_yyz_0[i] * fe_0 - ta1_x_yyyyy_yyz_1[i] * fe_0 + ta_yyyyy_xyyz_1[i] + ta1_x_yyyyy_xyyz_0[i] * pa_x[i] - ta1_x_yyyyy_xyyz_1[i] * pc_x[i];

        ta1_x_xyyyyy_xyzz_0[i] = ta1_x_yyyyy_yzz_0[i] * fe_0 - ta1_x_yyyyy_yzz_1[i] * fe_0 + ta_yyyyy_xyzz_1[i] + ta1_x_yyyyy_xyzz_0[i] * pa_x[i] - ta1_x_yyyyy_xyzz_1[i] * pc_x[i];

        ta1_x_xyyyyy_xzzz_0[i] = 4.0 * ta1_x_xyyy_xzzz_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xzzz_1[i] * fe_0 + ta1_x_xyyyy_xzzz_0[i] * pa_y[i] - ta1_x_xyyyy_xzzz_1[i] * pc_y[i];

        ta1_x_xyyyyy_yyyy_0[i] = ta_yyyyy_yyyy_1[i] + ta1_x_yyyyy_yyyy_0[i] * pa_x[i] - ta1_x_yyyyy_yyyy_1[i] * pc_x[i];

        ta1_x_xyyyyy_yyyz_0[i] = ta_yyyyy_yyyz_1[i] + ta1_x_yyyyy_yyyz_0[i] * pa_x[i] - ta1_x_yyyyy_yyyz_1[i] * pc_x[i];

        ta1_x_xyyyyy_yyzz_0[i] = ta_yyyyy_yyzz_1[i] + ta1_x_yyyyy_yyzz_0[i] * pa_x[i] - ta1_x_yyyyy_yyzz_1[i] * pc_x[i];

        ta1_x_xyyyyy_yzzz_0[i] = ta_yyyyy_yzzz_1[i] + ta1_x_yyyyy_yzzz_0[i] * pa_x[i] - ta1_x_yyyyy_yzzz_1[i] * pc_x[i];

        ta1_x_xyyyyy_zzzz_0[i] = ta_yyyyy_zzzz_1[i] + ta1_x_yyyyy_zzzz_0[i] * pa_x[i] - ta1_x_yyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : IG

    auto ta1_x_xyyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 240);

    auto ta1_x_xyyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 241);

    auto ta1_x_xyyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 242);

    auto ta1_x_xyyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 243);

    auto ta1_x_xyyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 244);

    auto ta1_x_xyyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 245);

    auto ta1_x_xyyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 246);

    auto ta1_x_xyyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 247);

    auto ta1_x_xyyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 248);

    auto ta1_x_xyyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 249);

    auto ta1_x_xyyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 250);

    auto ta1_x_xyyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 251);

    auto ta1_x_xyyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 252);

    auto ta1_x_xyyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 253);

    auto ta1_x_xyyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 254);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xyyyy_xxxx_0, ta1_x_xyyyy_xxxx_1, ta1_x_xyyyy_xxxy_0, ta1_x_xyyyy_xxxy_1, ta1_x_xyyyy_xxy_0, ta1_x_xyyyy_xxy_1, ta1_x_xyyyy_xxyy_0, ta1_x_xyyyy_xxyy_1, ta1_x_xyyyy_xxyz_0, ta1_x_xyyyy_xxyz_1, ta1_x_xyyyy_xyy_0, ta1_x_xyyyy_xyy_1, ta1_x_xyyyy_xyyy_0, ta1_x_xyyyy_xyyy_1, ta1_x_xyyyy_xyyz_0, ta1_x_xyyyy_xyyz_1, ta1_x_xyyyy_xyz_0, ta1_x_xyyyy_xyz_1, ta1_x_xyyyy_xyzz_0, ta1_x_xyyyy_xyzz_1, ta1_x_xyyyy_yyyy_0, ta1_x_xyyyy_yyyy_1, ta1_x_xyyyyz_xxxx_0, ta1_x_xyyyyz_xxxy_0, ta1_x_xyyyyz_xxxz_0, ta1_x_xyyyyz_xxyy_0, ta1_x_xyyyyz_xxyz_0, ta1_x_xyyyyz_xxzz_0, ta1_x_xyyyyz_xyyy_0, ta1_x_xyyyyz_xyyz_0, ta1_x_xyyyyz_xyzz_0, ta1_x_xyyyyz_xzzz_0, ta1_x_xyyyyz_yyyy_0, ta1_x_xyyyyz_yyyz_0, ta1_x_xyyyyz_yyzz_0, ta1_x_xyyyyz_yzzz_0, ta1_x_xyyyyz_zzzz_0, ta1_x_xyyyz_xxxz_0, ta1_x_xyyyz_xxxz_1, ta1_x_xyyyz_xxzz_0, ta1_x_xyyyz_xxzz_1, ta1_x_xyyyz_xzzz_0, ta1_x_xyyyz_xzzz_1, ta1_x_xyyz_xxxz_0, ta1_x_xyyz_xxxz_1, ta1_x_xyyz_xxzz_0, ta1_x_xyyz_xxzz_1, ta1_x_xyyz_xzzz_0, ta1_x_xyyz_xzzz_1, ta1_x_yyyyz_yyyz_0, ta1_x_yyyyz_yyyz_1, ta1_x_yyyyz_yyzz_0, ta1_x_yyyyz_yyzz_1, ta1_x_yyyyz_yzzz_0, ta1_x_yyyyz_yzzz_1, ta1_x_yyyyz_zzzz_0, ta1_x_yyyyz_zzzz_1, ta_yyyyz_yyyz_1, ta_yyyyz_yyzz_1, ta_yyyyz_yzzz_1, ta_yyyyz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyz_xxxx_0[i] = ta1_x_xyyyy_xxxx_0[i] * pa_z[i] - ta1_x_xyyyy_xxxx_1[i] * pc_z[i];

        ta1_x_xyyyyz_xxxy_0[i] = ta1_x_xyyyy_xxxy_0[i] * pa_z[i] - ta1_x_xyyyy_xxxy_1[i] * pc_z[i];

        ta1_x_xyyyyz_xxxz_0[i] = 3.0 * ta1_x_xyyz_xxxz_0[i] * fe_0 - 3.0 * ta1_x_xyyz_xxxz_1[i] * fe_0 + ta1_x_xyyyz_xxxz_0[i] * pa_y[i] - ta1_x_xyyyz_xxxz_1[i] * pc_y[i];

        ta1_x_xyyyyz_xxyy_0[i] = ta1_x_xyyyy_xxyy_0[i] * pa_z[i] - ta1_x_xyyyy_xxyy_1[i] * pc_z[i];

        ta1_x_xyyyyz_xxyz_0[i] = ta1_x_xyyyy_xxy_0[i] * fe_0 - ta1_x_xyyyy_xxy_1[i] * fe_0 + ta1_x_xyyyy_xxyz_0[i] * pa_z[i] - ta1_x_xyyyy_xxyz_1[i] * pc_z[i];

        ta1_x_xyyyyz_xxzz_0[i] = 3.0 * ta1_x_xyyz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xyyz_xxzz_1[i] * fe_0 + ta1_x_xyyyz_xxzz_0[i] * pa_y[i] - ta1_x_xyyyz_xxzz_1[i] * pc_y[i];

        ta1_x_xyyyyz_xyyy_0[i] = ta1_x_xyyyy_xyyy_0[i] * pa_z[i] - ta1_x_xyyyy_xyyy_1[i] * pc_z[i];

        ta1_x_xyyyyz_xyyz_0[i] = ta1_x_xyyyy_xyy_0[i] * fe_0 - ta1_x_xyyyy_xyy_1[i] * fe_0 + ta1_x_xyyyy_xyyz_0[i] * pa_z[i] - ta1_x_xyyyy_xyyz_1[i] * pc_z[i];

        ta1_x_xyyyyz_xyzz_0[i] = 2.0 * ta1_x_xyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xyyyy_xyz_1[i] * fe_0 + ta1_x_xyyyy_xyzz_0[i] * pa_z[i] - ta1_x_xyyyy_xyzz_1[i] * pc_z[i];

        ta1_x_xyyyyz_xzzz_0[i] = 3.0 * ta1_x_xyyz_xzzz_0[i] * fe_0 - 3.0 * ta1_x_xyyz_xzzz_1[i] * fe_0 + ta1_x_xyyyz_xzzz_0[i] * pa_y[i] - ta1_x_xyyyz_xzzz_1[i] * pc_y[i];

        ta1_x_xyyyyz_yyyy_0[i] = ta1_x_xyyyy_yyyy_0[i] * pa_z[i] - ta1_x_xyyyy_yyyy_1[i] * pc_z[i];

        ta1_x_xyyyyz_yyyz_0[i] = ta_yyyyz_yyyz_1[i] + ta1_x_yyyyz_yyyz_0[i] * pa_x[i] - ta1_x_yyyyz_yyyz_1[i] * pc_x[i];

        ta1_x_xyyyyz_yyzz_0[i] = ta_yyyyz_yyzz_1[i] + ta1_x_yyyyz_yyzz_0[i] * pa_x[i] - ta1_x_yyyyz_yyzz_1[i] * pc_x[i];

        ta1_x_xyyyyz_yzzz_0[i] = ta_yyyyz_yzzz_1[i] + ta1_x_yyyyz_yzzz_0[i] * pa_x[i] - ta1_x_yyyyz_yzzz_1[i] * pc_x[i];

        ta1_x_xyyyyz_zzzz_0[i] = ta_yyyyz_zzzz_1[i] + ta1_x_yyyyz_zzzz_0[i] * pa_x[i] - ta1_x_yyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 255-270 components of targeted buffer : IG

    auto ta1_x_xyyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 255);

    auto ta1_x_xyyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 256);

    auto ta1_x_xyyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 257);

    auto ta1_x_xyyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 258);

    auto ta1_x_xyyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 259);

    auto ta1_x_xyyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 260);

    auto ta1_x_xyyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 261);

    auto ta1_x_xyyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 262);

    auto ta1_x_xyyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 263);

    auto ta1_x_xyyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 264);

    auto ta1_x_xyyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 265);

    auto ta1_x_xyyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 266);

    auto ta1_x_xyyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 267);

    auto ta1_x_xyyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 268);

    auto ta1_x_xyyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 269);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xyyy_xxxy_0, ta1_x_xyyy_xxxy_1, ta1_x_xyyy_xxyy_0, ta1_x_xyyy_xxyy_1, ta1_x_xyyy_xyyy_0, ta1_x_xyyy_xyyy_1, ta1_x_xyyyz_xxxy_0, ta1_x_xyyyz_xxxy_1, ta1_x_xyyyz_xxyy_0, ta1_x_xyyyz_xxyy_1, ta1_x_xyyyz_xyyy_0, ta1_x_xyyyz_xyyy_1, ta1_x_xyyyzz_xxxx_0, ta1_x_xyyyzz_xxxy_0, ta1_x_xyyyzz_xxxz_0, ta1_x_xyyyzz_xxyy_0, ta1_x_xyyyzz_xxyz_0, ta1_x_xyyyzz_xxzz_0, ta1_x_xyyyzz_xyyy_0, ta1_x_xyyyzz_xyyz_0, ta1_x_xyyyzz_xyzz_0, ta1_x_xyyyzz_xzzz_0, ta1_x_xyyyzz_yyyy_0, ta1_x_xyyyzz_yyyz_0, ta1_x_xyyyzz_yyzz_0, ta1_x_xyyyzz_yzzz_0, ta1_x_xyyyzz_zzzz_0, ta1_x_xyyzz_xxxx_0, ta1_x_xyyzz_xxxx_1, ta1_x_xyyzz_xxxz_0, ta1_x_xyyzz_xxxz_1, ta1_x_xyyzz_xxzz_0, ta1_x_xyyzz_xxzz_1, ta1_x_xyyzz_xzzz_0, ta1_x_xyyzz_xzzz_1, ta1_x_xyzz_xxxx_0, ta1_x_xyzz_xxxx_1, ta1_x_xyzz_xxxz_0, ta1_x_xyzz_xxxz_1, ta1_x_xyzz_xxzz_0, ta1_x_xyzz_xxzz_1, ta1_x_xyzz_xzzz_0, ta1_x_xyzz_xzzz_1, ta1_x_yyyzz_xxyz_0, ta1_x_yyyzz_xxyz_1, ta1_x_yyyzz_xyyz_0, ta1_x_yyyzz_xyyz_1, ta1_x_yyyzz_xyz_0, ta1_x_yyyzz_xyz_1, ta1_x_yyyzz_xyzz_0, ta1_x_yyyzz_xyzz_1, ta1_x_yyyzz_yyyy_0, ta1_x_yyyzz_yyyy_1, ta1_x_yyyzz_yyyz_0, ta1_x_yyyzz_yyyz_1, ta1_x_yyyzz_yyz_0, ta1_x_yyyzz_yyz_1, ta1_x_yyyzz_yyzz_0, ta1_x_yyyzz_yyzz_1, ta1_x_yyyzz_yzz_0, ta1_x_yyyzz_yzz_1, ta1_x_yyyzz_yzzz_0, ta1_x_yyyzz_yzzz_1, ta1_x_yyyzz_zzzz_0, ta1_x_yyyzz_zzzz_1, ta_yyyzz_xxyz_1, ta_yyyzz_xyyz_1, ta_yyyzz_xyzz_1, ta_yyyzz_yyyy_1, ta_yyyzz_yyyz_1, ta_yyyzz_yyzz_1, ta_yyyzz_yzzz_1, ta_yyyzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyzz_xxxx_0[i] = 2.0 * ta1_x_xyzz_xxxx_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xxxx_1[i] * fe_0 + ta1_x_xyyzz_xxxx_0[i] * pa_y[i] - ta1_x_xyyzz_xxxx_1[i] * pc_y[i];

        ta1_x_xyyyzz_xxxy_0[i] = ta1_x_xyyy_xxxy_0[i] * fe_0 - ta1_x_xyyy_xxxy_1[i] * fe_0 + ta1_x_xyyyz_xxxy_0[i] * pa_z[i] - ta1_x_xyyyz_xxxy_1[i] * pc_z[i];

        ta1_x_xyyyzz_xxxz_0[i] = 2.0 * ta1_x_xyzz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xxxz_1[i] * fe_0 + ta1_x_xyyzz_xxxz_0[i] * pa_y[i] - ta1_x_xyyzz_xxxz_1[i] * pc_y[i];

        ta1_x_xyyyzz_xxyy_0[i] = ta1_x_xyyy_xxyy_0[i] * fe_0 - ta1_x_xyyy_xxyy_1[i] * fe_0 + ta1_x_xyyyz_xxyy_0[i] * pa_z[i] - ta1_x_xyyyz_xxyy_1[i] * pc_z[i];

        ta1_x_xyyyzz_xxyz_0[i] = 2.0 * ta1_x_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyzz_xyz_1[i] * fe_0 + ta_yyyzz_xxyz_1[i] + ta1_x_yyyzz_xxyz_0[i] * pa_x[i] - ta1_x_yyyzz_xxyz_1[i] * pc_x[i];

        ta1_x_xyyyzz_xxzz_0[i] = 2.0 * ta1_x_xyzz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xxzz_1[i] * fe_0 + ta1_x_xyyzz_xxzz_0[i] * pa_y[i] - ta1_x_xyyzz_xxzz_1[i] * pc_y[i];

        ta1_x_xyyyzz_xyyy_0[i] = ta1_x_xyyy_xyyy_0[i] * fe_0 - ta1_x_xyyy_xyyy_1[i] * fe_0 + ta1_x_xyyyz_xyyy_0[i] * pa_z[i] - ta1_x_xyyyz_xyyy_1[i] * pc_z[i];

        ta1_x_xyyyzz_xyyz_0[i] = ta1_x_yyyzz_yyz_0[i] * fe_0 - ta1_x_yyyzz_yyz_1[i] * fe_0 + ta_yyyzz_xyyz_1[i] + ta1_x_yyyzz_xyyz_0[i] * pa_x[i] - ta1_x_yyyzz_xyyz_1[i] * pc_x[i];

        ta1_x_xyyyzz_xyzz_0[i] = ta1_x_yyyzz_yzz_0[i] * fe_0 - ta1_x_yyyzz_yzz_1[i] * fe_0 + ta_yyyzz_xyzz_1[i] + ta1_x_yyyzz_xyzz_0[i] * pa_x[i] - ta1_x_yyyzz_xyzz_1[i] * pc_x[i];

        ta1_x_xyyyzz_xzzz_0[i] = 2.0 * ta1_x_xyzz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xzzz_1[i] * fe_0 + ta1_x_xyyzz_xzzz_0[i] * pa_y[i] - ta1_x_xyyzz_xzzz_1[i] * pc_y[i];

        ta1_x_xyyyzz_yyyy_0[i] = ta_yyyzz_yyyy_1[i] + ta1_x_yyyzz_yyyy_0[i] * pa_x[i] - ta1_x_yyyzz_yyyy_1[i] * pc_x[i];

        ta1_x_xyyyzz_yyyz_0[i] = ta_yyyzz_yyyz_1[i] + ta1_x_yyyzz_yyyz_0[i] * pa_x[i] - ta1_x_yyyzz_yyyz_1[i] * pc_x[i];

        ta1_x_xyyyzz_yyzz_0[i] = ta_yyyzz_yyzz_1[i] + ta1_x_yyyzz_yyzz_0[i] * pa_x[i] - ta1_x_yyyzz_yyzz_1[i] * pc_x[i];

        ta1_x_xyyyzz_yzzz_0[i] = ta_yyyzz_yzzz_1[i] + ta1_x_yyyzz_yzzz_0[i] * pa_x[i] - ta1_x_yyyzz_yzzz_1[i] * pc_x[i];

        ta1_x_xyyyzz_zzzz_0[i] = ta_yyyzz_zzzz_1[i] + ta1_x_yyyzz_zzzz_0[i] * pa_x[i] - ta1_x_yyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 270-285 components of targeted buffer : IG

    auto ta1_x_xyyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 270);

    auto ta1_x_xyyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 271);

    auto ta1_x_xyyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 272);

    auto ta1_x_xyyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 273);

    auto ta1_x_xyyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 274);

    auto ta1_x_xyyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 275);

    auto ta1_x_xyyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 276);

    auto ta1_x_xyyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 277);

    auto ta1_x_xyyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 278);

    auto ta1_x_xyyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 279);

    auto ta1_x_xyyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 280);

    auto ta1_x_xyyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 281);

    auto ta1_x_xyyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 282);

    auto ta1_x_xyyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 283);

    auto ta1_x_xyyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 284);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xyyz_xxxy_0, ta1_x_xyyz_xxxy_1, ta1_x_xyyz_xxyy_0, ta1_x_xyyz_xxyy_1, ta1_x_xyyz_xyyy_0, ta1_x_xyyz_xyyy_1, ta1_x_xyyzz_xxxy_0, ta1_x_xyyzz_xxxy_1, ta1_x_xyyzz_xxyy_0, ta1_x_xyyzz_xxyy_1, ta1_x_xyyzz_xyyy_0, ta1_x_xyyzz_xyyy_1, ta1_x_xyyzzz_xxxx_0, ta1_x_xyyzzz_xxxy_0, ta1_x_xyyzzz_xxxz_0, ta1_x_xyyzzz_xxyy_0, ta1_x_xyyzzz_xxyz_0, ta1_x_xyyzzz_xxzz_0, ta1_x_xyyzzz_xyyy_0, ta1_x_xyyzzz_xyyz_0, ta1_x_xyyzzz_xyzz_0, ta1_x_xyyzzz_xzzz_0, ta1_x_xyyzzz_yyyy_0, ta1_x_xyyzzz_yyyz_0, ta1_x_xyyzzz_yyzz_0, ta1_x_xyyzzz_yzzz_0, ta1_x_xyyzzz_zzzz_0, ta1_x_xyzzz_xxxx_0, ta1_x_xyzzz_xxxx_1, ta1_x_xyzzz_xxxz_0, ta1_x_xyzzz_xxxz_1, ta1_x_xyzzz_xxzz_0, ta1_x_xyzzz_xxzz_1, ta1_x_xyzzz_xzzz_0, ta1_x_xyzzz_xzzz_1, ta1_x_xzzz_xxxx_0, ta1_x_xzzz_xxxx_1, ta1_x_xzzz_xxxz_0, ta1_x_xzzz_xxxz_1, ta1_x_xzzz_xxzz_0, ta1_x_xzzz_xxzz_1, ta1_x_xzzz_xzzz_0, ta1_x_xzzz_xzzz_1, ta1_x_yyzzz_xxyz_0, ta1_x_yyzzz_xxyz_1, ta1_x_yyzzz_xyyz_0, ta1_x_yyzzz_xyyz_1, ta1_x_yyzzz_xyz_0, ta1_x_yyzzz_xyz_1, ta1_x_yyzzz_xyzz_0, ta1_x_yyzzz_xyzz_1, ta1_x_yyzzz_yyyy_0, ta1_x_yyzzz_yyyy_1, ta1_x_yyzzz_yyyz_0, ta1_x_yyzzz_yyyz_1, ta1_x_yyzzz_yyz_0, ta1_x_yyzzz_yyz_1, ta1_x_yyzzz_yyzz_0, ta1_x_yyzzz_yyzz_1, ta1_x_yyzzz_yzz_0, ta1_x_yyzzz_yzz_1, ta1_x_yyzzz_yzzz_0, ta1_x_yyzzz_yzzz_1, ta1_x_yyzzz_zzzz_0, ta1_x_yyzzz_zzzz_1, ta_yyzzz_xxyz_1, ta_yyzzz_xyyz_1, ta_yyzzz_xyzz_1, ta_yyzzz_yyyy_1, ta_yyzzz_yyyz_1, ta_yyzzz_yyzz_1, ta_yyzzz_yzzz_1, ta_yyzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzzz_xxxx_0[i] = ta1_x_xzzz_xxxx_0[i] * fe_0 - ta1_x_xzzz_xxxx_1[i] * fe_0 + ta1_x_xyzzz_xxxx_0[i] * pa_y[i] - ta1_x_xyzzz_xxxx_1[i] * pc_y[i];

        ta1_x_xyyzzz_xxxy_0[i] = 2.0 * ta1_x_xyyz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xyyz_xxxy_1[i] * fe_0 + ta1_x_xyyzz_xxxy_0[i] * pa_z[i] - ta1_x_xyyzz_xxxy_1[i] * pc_z[i];

        ta1_x_xyyzzz_xxxz_0[i] = ta1_x_xzzz_xxxz_0[i] * fe_0 - ta1_x_xzzz_xxxz_1[i] * fe_0 + ta1_x_xyzzz_xxxz_0[i] * pa_y[i] - ta1_x_xyzzz_xxxz_1[i] * pc_y[i];

        ta1_x_xyyzzz_xxyy_0[i] = 2.0 * ta1_x_xyyz_xxyy_0[i] * fe_0 - 2.0 * ta1_x_xyyz_xxyy_1[i] * fe_0 + ta1_x_xyyzz_xxyy_0[i] * pa_z[i] - ta1_x_xyyzz_xxyy_1[i] * pc_z[i];

        ta1_x_xyyzzz_xxyz_0[i] = 2.0 * ta1_x_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyzzz_xyz_1[i] * fe_0 + ta_yyzzz_xxyz_1[i] + ta1_x_yyzzz_xxyz_0[i] * pa_x[i] - ta1_x_yyzzz_xxyz_1[i] * pc_x[i];

        ta1_x_xyyzzz_xxzz_0[i] = ta1_x_xzzz_xxzz_0[i] * fe_0 - ta1_x_xzzz_xxzz_1[i] * fe_0 + ta1_x_xyzzz_xxzz_0[i] * pa_y[i] - ta1_x_xyzzz_xxzz_1[i] * pc_y[i];

        ta1_x_xyyzzz_xyyy_0[i] = 2.0 * ta1_x_xyyz_xyyy_0[i] * fe_0 - 2.0 * ta1_x_xyyz_xyyy_1[i] * fe_0 + ta1_x_xyyzz_xyyy_0[i] * pa_z[i] - ta1_x_xyyzz_xyyy_1[i] * pc_z[i];

        ta1_x_xyyzzz_xyyz_0[i] = ta1_x_yyzzz_yyz_0[i] * fe_0 - ta1_x_yyzzz_yyz_1[i] * fe_0 + ta_yyzzz_xyyz_1[i] + ta1_x_yyzzz_xyyz_0[i] * pa_x[i] - ta1_x_yyzzz_xyyz_1[i] * pc_x[i];

        ta1_x_xyyzzz_xyzz_0[i] = ta1_x_yyzzz_yzz_0[i] * fe_0 - ta1_x_yyzzz_yzz_1[i] * fe_0 + ta_yyzzz_xyzz_1[i] + ta1_x_yyzzz_xyzz_0[i] * pa_x[i] - ta1_x_yyzzz_xyzz_1[i] * pc_x[i];

        ta1_x_xyyzzz_xzzz_0[i] = ta1_x_xzzz_xzzz_0[i] * fe_0 - ta1_x_xzzz_xzzz_1[i] * fe_0 + ta1_x_xyzzz_xzzz_0[i] * pa_y[i] - ta1_x_xyzzz_xzzz_1[i] * pc_y[i];

        ta1_x_xyyzzz_yyyy_0[i] = ta_yyzzz_yyyy_1[i] + ta1_x_yyzzz_yyyy_0[i] * pa_x[i] - ta1_x_yyzzz_yyyy_1[i] * pc_x[i];

        ta1_x_xyyzzz_yyyz_0[i] = ta_yyzzz_yyyz_1[i] + ta1_x_yyzzz_yyyz_0[i] * pa_x[i] - ta1_x_yyzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xyyzzz_yyzz_0[i] = ta_yyzzz_yyzz_1[i] + ta1_x_yyzzz_yyzz_0[i] * pa_x[i] - ta1_x_yyzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xyyzzz_yzzz_0[i] = ta_yyzzz_yzzz_1[i] + ta1_x_yyzzz_yzzz_0[i] * pa_x[i] - ta1_x_yyzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xyyzzz_zzzz_0[i] = ta_yyzzz_zzzz_1[i] + ta1_x_yyzzz_zzzz_0[i] * pa_x[i] - ta1_x_yyzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 285-300 components of targeted buffer : IG

    auto ta1_x_xyzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 285);

    auto ta1_x_xyzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 286);

    auto ta1_x_xyzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 287);

    auto ta1_x_xyzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 288);

    auto ta1_x_xyzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 289);

    auto ta1_x_xyzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 290);

    auto ta1_x_xyzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 291);

    auto ta1_x_xyzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 292);

    auto ta1_x_xyzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 293);

    auto ta1_x_xyzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 294);

    auto ta1_x_xyzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 295);

    auto ta1_x_xyzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 296);

    auto ta1_x_xyzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 297);

    auto ta1_x_xyzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 298);

    auto ta1_x_xyzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 299);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyzzzz_xxxx_0, ta1_x_xyzzzz_xxxy_0, ta1_x_xyzzzz_xxxz_0, ta1_x_xyzzzz_xxyy_0, ta1_x_xyzzzz_xxyz_0, ta1_x_xyzzzz_xxzz_0, ta1_x_xyzzzz_xyyy_0, ta1_x_xyzzzz_xyyz_0, ta1_x_xyzzzz_xyzz_0, ta1_x_xyzzzz_xzzz_0, ta1_x_xyzzzz_yyyy_0, ta1_x_xyzzzz_yyyz_0, ta1_x_xyzzzz_yyzz_0, ta1_x_xyzzzz_yzzz_0, ta1_x_xyzzzz_zzzz_0, ta1_x_xzzzz_xxx_0, ta1_x_xzzzz_xxx_1, ta1_x_xzzzz_xxxx_0, ta1_x_xzzzz_xxxx_1, ta1_x_xzzzz_xxxy_0, ta1_x_xzzzz_xxxy_1, ta1_x_xzzzz_xxxz_0, ta1_x_xzzzz_xxxz_1, ta1_x_xzzzz_xxy_0, ta1_x_xzzzz_xxy_1, ta1_x_xzzzz_xxyy_0, ta1_x_xzzzz_xxyy_1, ta1_x_xzzzz_xxyz_0, ta1_x_xzzzz_xxyz_1, ta1_x_xzzzz_xxz_0, ta1_x_xzzzz_xxz_1, ta1_x_xzzzz_xxzz_0, ta1_x_xzzzz_xxzz_1, ta1_x_xzzzz_xyy_0, ta1_x_xzzzz_xyy_1, ta1_x_xzzzz_xyyy_0, ta1_x_xzzzz_xyyy_1, ta1_x_xzzzz_xyyz_0, ta1_x_xzzzz_xyyz_1, ta1_x_xzzzz_xyz_0, ta1_x_xzzzz_xyz_1, ta1_x_xzzzz_xyzz_0, ta1_x_xzzzz_xyzz_1, ta1_x_xzzzz_xzz_0, ta1_x_xzzzz_xzz_1, ta1_x_xzzzz_xzzz_0, ta1_x_xzzzz_xzzz_1, ta1_x_xzzzz_zzzz_0, ta1_x_xzzzz_zzzz_1, ta1_x_yzzzz_yyyy_0, ta1_x_yzzzz_yyyy_1, ta1_x_yzzzz_yyyz_0, ta1_x_yzzzz_yyyz_1, ta1_x_yzzzz_yyzz_0, ta1_x_yzzzz_yyzz_1, ta1_x_yzzzz_yzzz_0, ta1_x_yzzzz_yzzz_1, ta_yzzzz_yyyy_1, ta_yzzzz_yyyz_1, ta_yzzzz_yyzz_1, ta_yzzzz_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzzzz_xxxx_0[i] = ta1_x_xzzzz_xxxx_0[i] * pa_y[i] - ta1_x_xzzzz_xxxx_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxxy_0[i] = ta1_x_xzzzz_xxx_0[i] * fe_0 - ta1_x_xzzzz_xxx_1[i] * fe_0 + ta1_x_xzzzz_xxxy_0[i] * pa_y[i] - ta1_x_xzzzz_xxxy_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxxz_0[i] = ta1_x_xzzzz_xxxz_0[i] * pa_y[i] - ta1_x_xzzzz_xxxz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxyy_0[i] = 2.0 * ta1_x_xzzzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xzzzz_xxy_1[i] * fe_0 + ta1_x_xzzzz_xxyy_0[i] * pa_y[i] - ta1_x_xzzzz_xxyy_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxyz_0[i] = ta1_x_xzzzz_xxz_0[i] * fe_0 - ta1_x_xzzzz_xxz_1[i] * fe_0 + ta1_x_xzzzz_xxyz_0[i] * pa_y[i] - ta1_x_xzzzz_xxyz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxzz_0[i] = ta1_x_xzzzz_xxzz_0[i] * pa_y[i] - ta1_x_xzzzz_xxzz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xyyy_0[i] = 3.0 * ta1_x_xzzzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xzzzz_xyy_1[i] * fe_0 + ta1_x_xzzzz_xyyy_0[i] * pa_y[i] - ta1_x_xzzzz_xyyy_1[i] * pc_y[i];

        ta1_x_xyzzzz_xyyz_0[i] = 2.0 * ta1_x_xzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xzzzz_xyz_1[i] * fe_0 + ta1_x_xzzzz_xyyz_0[i] * pa_y[i] - ta1_x_xzzzz_xyyz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xyzz_0[i] = ta1_x_xzzzz_xzz_0[i] * fe_0 - ta1_x_xzzzz_xzz_1[i] * fe_0 + ta1_x_xzzzz_xyzz_0[i] * pa_y[i] - ta1_x_xzzzz_xyzz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xzzz_0[i] = ta1_x_xzzzz_xzzz_0[i] * pa_y[i] - ta1_x_xzzzz_xzzz_1[i] * pc_y[i];

        ta1_x_xyzzzz_yyyy_0[i] = ta_yzzzz_yyyy_1[i] + ta1_x_yzzzz_yyyy_0[i] * pa_x[i] - ta1_x_yzzzz_yyyy_1[i] * pc_x[i];

        ta1_x_xyzzzz_yyyz_0[i] = ta_yzzzz_yyyz_1[i] + ta1_x_yzzzz_yyyz_0[i] * pa_x[i] - ta1_x_yzzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xyzzzz_yyzz_0[i] = ta_yzzzz_yyzz_1[i] + ta1_x_yzzzz_yyzz_0[i] * pa_x[i] - ta1_x_yzzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xyzzzz_yzzz_0[i] = ta_yzzzz_yzzz_1[i] + ta1_x_yzzzz_yzzz_0[i] * pa_x[i] - ta1_x_yzzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xyzzzz_zzzz_0[i] = ta1_x_xzzzz_zzzz_0[i] * pa_y[i] - ta1_x_xzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 300-315 components of targeted buffer : IG

    auto ta1_x_xzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 300);

    auto ta1_x_xzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 301);

    auto ta1_x_xzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 302);

    auto ta1_x_xzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 303);

    auto ta1_x_xzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 304);

    auto ta1_x_xzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 305);

    auto ta1_x_xzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 306);

    auto ta1_x_xzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 307);

    auto ta1_x_xzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 308);

    auto ta1_x_xzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 309);

    auto ta1_x_xzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 310);

    auto ta1_x_xzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 311);

    auto ta1_x_xzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 312);

    auto ta1_x_xzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 313);

    auto ta1_x_xzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 314);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xzzz_xxxx_0, ta1_x_xzzz_xxxx_1, ta1_x_xzzz_xxxy_0, ta1_x_xzzz_xxxy_1, ta1_x_xzzz_xxyy_0, ta1_x_xzzz_xxyy_1, ta1_x_xzzz_xyyy_0, ta1_x_xzzz_xyyy_1, ta1_x_xzzzz_xxxx_0, ta1_x_xzzzz_xxxx_1, ta1_x_xzzzz_xxxy_0, ta1_x_xzzzz_xxxy_1, ta1_x_xzzzz_xxyy_0, ta1_x_xzzzz_xxyy_1, ta1_x_xzzzz_xyyy_0, ta1_x_xzzzz_xyyy_1, ta1_x_xzzzzz_xxxx_0, ta1_x_xzzzzz_xxxy_0, ta1_x_xzzzzz_xxxz_0, ta1_x_xzzzzz_xxyy_0, ta1_x_xzzzzz_xxyz_0, ta1_x_xzzzzz_xxzz_0, ta1_x_xzzzzz_xyyy_0, ta1_x_xzzzzz_xyyz_0, ta1_x_xzzzzz_xyzz_0, ta1_x_xzzzzz_xzzz_0, ta1_x_xzzzzz_yyyy_0, ta1_x_xzzzzz_yyyz_0, ta1_x_xzzzzz_yyzz_0, ta1_x_xzzzzz_yzzz_0, ta1_x_xzzzzz_zzzz_0, ta1_x_zzzzz_xxxz_0, ta1_x_zzzzz_xxxz_1, ta1_x_zzzzz_xxyz_0, ta1_x_zzzzz_xxyz_1, ta1_x_zzzzz_xxz_0, ta1_x_zzzzz_xxz_1, ta1_x_zzzzz_xxzz_0, ta1_x_zzzzz_xxzz_1, ta1_x_zzzzz_xyyz_0, ta1_x_zzzzz_xyyz_1, ta1_x_zzzzz_xyz_0, ta1_x_zzzzz_xyz_1, ta1_x_zzzzz_xyzz_0, ta1_x_zzzzz_xyzz_1, ta1_x_zzzzz_xzz_0, ta1_x_zzzzz_xzz_1, ta1_x_zzzzz_xzzz_0, ta1_x_zzzzz_xzzz_1, ta1_x_zzzzz_yyyy_0, ta1_x_zzzzz_yyyy_1, ta1_x_zzzzz_yyyz_0, ta1_x_zzzzz_yyyz_1, ta1_x_zzzzz_yyz_0, ta1_x_zzzzz_yyz_1, ta1_x_zzzzz_yyzz_0, ta1_x_zzzzz_yyzz_1, ta1_x_zzzzz_yzz_0, ta1_x_zzzzz_yzz_1, ta1_x_zzzzz_yzzz_0, ta1_x_zzzzz_yzzz_1, ta1_x_zzzzz_zzz_0, ta1_x_zzzzz_zzz_1, ta1_x_zzzzz_zzzz_0, ta1_x_zzzzz_zzzz_1, ta_zzzzz_xxxz_1, ta_zzzzz_xxyz_1, ta_zzzzz_xxzz_1, ta_zzzzz_xyyz_1, ta_zzzzz_xyzz_1, ta_zzzzz_xzzz_1, ta_zzzzz_yyyy_1, ta_zzzzz_yyyz_1, ta_zzzzz_yyzz_1, ta_zzzzz_yzzz_1, ta_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzzz_xxxx_0[i] = 4.0 * ta1_x_xzzz_xxxx_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xxxx_1[i] * fe_0 + ta1_x_xzzzz_xxxx_0[i] * pa_z[i] - ta1_x_xzzzz_xxxx_1[i] * pc_z[i];

        ta1_x_xzzzzz_xxxy_0[i] = 4.0 * ta1_x_xzzz_xxxy_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xxxy_1[i] * fe_0 + ta1_x_xzzzz_xxxy_0[i] * pa_z[i] - ta1_x_xzzzz_xxxy_1[i] * pc_z[i];

        ta1_x_xzzzzz_xxxz_0[i] = 3.0 * ta1_x_zzzzz_xxz_0[i] * fe_0 - 3.0 * ta1_x_zzzzz_xxz_1[i] * fe_0 + ta_zzzzz_xxxz_1[i] + ta1_x_zzzzz_xxxz_0[i] * pa_x[i] - ta1_x_zzzzz_xxxz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xxyy_0[i] = 4.0 * ta1_x_xzzz_xxyy_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xxyy_1[i] * fe_0 + ta1_x_xzzzz_xxyy_0[i] * pa_z[i] - ta1_x_xzzzz_xxyy_1[i] * pc_z[i];

        ta1_x_xzzzzz_xxyz_0[i] = 2.0 * ta1_x_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xyz_1[i] * fe_0 + ta_zzzzz_xxyz_1[i] + ta1_x_zzzzz_xxyz_0[i] * pa_x[i] - ta1_x_zzzzz_xxyz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xxzz_0[i] = 2.0 * ta1_x_zzzzz_xzz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xzz_1[i] * fe_0 + ta_zzzzz_xxzz_1[i] + ta1_x_zzzzz_xxzz_0[i] * pa_x[i] - ta1_x_zzzzz_xxzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xyyy_0[i] = 4.0 * ta1_x_xzzz_xyyy_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xyyy_1[i] * fe_0 + ta1_x_xzzzz_xyyy_0[i] * pa_z[i] - ta1_x_xzzzz_xyyy_1[i] * pc_z[i];

        ta1_x_xzzzzz_xyyz_0[i] = ta1_x_zzzzz_yyz_0[i] * fe_0 - ta1_x_zzzzz_yyz_1[i] * fe_0 + ta_zzzzz_xyyz_1[i] + ta1_x_zzzzz_xyyz_0[i] * pa_x[i] - ta1_x_zzzzz_xyyz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xyzz_0[i] = ta1_x_zzzzz_yzz_0[i] * fe_0 - ta1_x_zzzzz_yzz_1[i] * fe_0 + ta_zzzzz_xyzz_1[i] + ta1_x_zzzzz_xyzz_0[i] * pa_x[i] - ta1_x_zzzzz_xyzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xzzz_0[i] = ta1_x_zzzzz_zzz_0[i] * fe_0 - ta1_x_zzzzz_zzz_1[i] * fe_0 + ta_zzzzz_xzzz_1[i] + ta1_x_zzzzz_xzzz_0[i] * pa_x[i] - ta1_x_zzzzz_xzzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_yyyy_0[i] = ta_zzzzz_yyyy_1[i] + ta1_x_zzzzz_yyyy_0[i] * pa_x[i] - ta1_x_zzzzz_yyyy_1[i] * pc_x[i];

        ta1_x_xzzzzz_yyyz_0[i] = ta_zzzzz_yyyz_1[i] + ta1_x_zzzzz_yyyz_0[i] * pa_x[i] - ta1_x_zzzzz_yyyz_1[i] * pc_x[i];

        ta1_x_xzzzzz_yyzz_0[i] = ta_zzzzz_yyzz_1[i] + ta1_x_zzzzz_yyzz_0[i] * pa_x[i] - ta1_x_zzzzz_yyzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_yzzz_0[i] = ta_zzzzz_yzzz_1[i] + ta1_x_zzzzz_yzzz_0[i] * pa_x[i] - ta1_x_zzzzz_yzzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_zzzz_0[i] = ta_zzzzz_zzzz_1[i] + ta1_x_zzzzz_zzzz_0[i] * pa_x[i] - ta1_x_zzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : IG

    auto ta1_x_yyyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 315);

    auto ta1_x_yyyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 316);

    auto ta1_x_yyyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 317);

    auto ta1_x_yyyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 318);

    auto ta1_x_yyyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 319);

    auto ta1_x_yyyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 320);

    auto ta1_x_yyyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 321);

    auto ta1_x_yyyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 322);

    auto ta1_x_yyyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 323);

    auto ta1_x_yyyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 324);

    auto ta1_x_yyyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 325);

    auto ta1_x_yyyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 326);

    auto ta1_x_yyyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 327);

    auto ta1_x_yyyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 328);

    auto ta1_x_yyyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 329);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yyyy_xxxx_0, ta1_x_yyyy_xxxx_1, ta1_x_yyyy_xxxy_0, ta1_x_yyyy_xxxy_1, ta1_x_yyyy_xxxz_0, ta1_x_yyyy_xxxz_1, ta1_x_yyyy_xxyy_0, ta1_x_yyyy_xxyy_1, ta1_x_yyyy_xxyz_0, ta1_x_yyyy_xxyz_1, ta1_x_yyyy_xxzz_0, ta1_x_yyyy_xxzz_1, ta1_x_yyyy_xyyy_0, ta1_x_yyyy_xyyy_1, ta1_x_yyyy_xyyz_0, ta1_x_yyyy_xyyz_1, ta1_x_yyyy_xyzz_0, ta1_x_yyyy_xyzz_1, ta1_x_yyyy_xzzz_0, ta1_x_yyyy_xzzz_1, ta1_x_yyyy_yyyy_0, ta1_x_yyyy_yyyy_1, ta1_x_yyyy_yyyz_0, ta1_x_yyyy_yyyz_1, ta1_x_yyyy_yyzz_0, ta1_x_yyyy_yyzz_1, ta1_x_yyyy_yzzz_0, ta1_x_yyyy_yzzz_1, ta1_x_yyyy_zzzz_0, ta1_x_yyyy_zzzz_1, ta1_x_yyyyy_xxx_0, ta1_x_yyyyy_xxx_1, ta1_x_yyyyy_xxxx_0, ta1_x_yyyyy_xxxx_1, ta1_x_yyyyy_xxxy_0, ta1_x_yyyyy_xxxy_1, ta1_x_yyyyy_xxxz_0, ta1_x_yyyyy_xxxz_1, ta1_x_yyyyy_xxy_0, ta1_x_yyyyy_xxy_1, ta1_x_yyyyy_xxyy_0, ta1_x_yyyyy_xxyy_1, ta1_x_yyyyy_xxyz_0, ta1_x_yyyyy_xxyz_1, ta1_x_yyyyy_xxz_0, ta1_x_yyyyy_xxz_1, ta1_x_yyyyy_xxzz_0, ta1_x_yyyyy_xxzz_1, ta1_x_yyyyy_xyy_0, ta1_x_yyyyy_xyy_1, ta1_x_yyyyy_xyyy_0, ta1_x_yyyyy_xyyy_1, ta1_x_yyyyy_xyyz_0, ta1_x_yyyyy_xyyz_1, ta1_x_yyyyy_xyz_0, ta1_x_yyyyy_xyz_1, ta1_x_yyyyy_xyzz_0, ta1_x_yyyyy_xyzz_1, ta1_x_yyyyy_xzz_0, ta1_x_yyyyy_xzz_1, ta1_x_yyyyy_xzzz_0, ta1_x_yyyyy_xzzz_1, ta1_x_yyyyy_yyy_0, ta1_x_yyyyy_yyy_1, ta1_x_yyyyy_yyyy_0, ta1_x_yyyyy_yyyy_1, ta1_x_yyyyy_yyyz_0, ta1_x_yyyyy_yyyz_1, ta1_x_yyyyy_yyz_0, ta1_x_yyyyy_yyz_1, ta1_x_yyyyy_yyzz_0, ta1_x_yyyyy_yyzz_1, ta1_x_yyyyy_yzz_0, ta1_x_yyyyy_yzz_1, ta1_x_yyyyy_yzzz_0, ta1_x_yyyyy_yzzz_1, ta1_x_yyyyy_zzz_0, ta1_x_yyyyy_zzz_1, ta1_x_yyyyy_zzzz_0, ta1_x_yyyyy_zzzz_1, ta1_x_yyyyyy_xxxx_0, ta1_x_yyyyyy_xxxy_0, ta1_x_yyyyyy_xxxz_0, ta1_x_yyyyyy_xxyy_0, ta1_x_yyyyyy_xxyz_0, ta1_x_yyyyyy_xxzz_0, ta1_x_yyyyyy_xyyy_0, ta1_x_yyyyyy_xyyz_0, ta1_x_yyyyyy_xyzz_0, ta1_x_yyyyyy_xzzz_0, ta1_x_yyyyyy_yyyy_0, ta1_x_yyyyyy_yyyz_0, ta1_x_yyyyyy_yyzz_0, ta1_x_yyyyyy_yzzz_0, ta1_x_yyyyyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyy_xxxx_0[i] = 5.0 * ta1_x_yyyy_xxxx_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxxx_1[i] * fe_0 + ta1_x_yyyyy_xxxx_0[i] * pa_y[i] - ta1_x_yyyyy_xxxx_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxxy_0[i] = 5.0 * ta1_x_yyyy_xxxy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxxy_1[i] * fe_0 + ta1_x_yyyyy_xxx_0[i] * fe_0 - ta1_x_yyyyy_xxx_1[i] * fe_0 + ta1_x_yyyyy_xxxy_0[i] * pa_y[i] - ta1_x_yyyyy_xxxy_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxxz_0[i] = 5.0 * ta1_x_yyyy_xxxz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxxz_1[i] * fe_0 + ta1_x_yyyyy_xxxz_0[i] * pa_y[i] - ta1_x_yyyyy_xxxz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxyy_0[i] = 5.0 * ta1_x_yyyy_xxyy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_x_yyyyy_xxy_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_xxy_1[i] * fe_0 + ta1_x_yyyyy_xxyy_0[i] * pa_y[i] - ta1_x_yyyyy_xxyy_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxyz_0[i] = 5.0 * ta1_x_yyyy_xxyz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxyz_1[i] * fe_0 + ta1_x_yyyyy_xxz_0[i] * fe_0 - ta1_x_yyyyy_xxz_1[i] * fe_0 + ta1_x_yyyyy_xxyz_0[i] * pa_y[i] - ta1_x_yyyyy_xxyz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxzz_0[i] = 5.0 * ta1_x_yyyy_xxzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxzz_1[i] * fe_0 + ta1_x_yyyyy_xxzz_0[i] * pa_y[i] - ta1_x_yyyyy_xxzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xyyy_0[i] = 5.0 * ta1_x_yyyy_xyyy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xyyy_1[i] * fe_0 + 3.0 * ta1_x_yyyyy_xyy_0[i] * fe_0 - 3.0 * ta1_x_yyyyy_xyy_1[i] * fe_0 + ta1_x_yyyyy_xyyy_0[i] * pa_y[i] - ta1_x_yyyyy_xyyy_1[i] * pc_y[i];

        ta1_x_yyyyyy_xyyz_0[i] = 5.0 * ta1_x_yyyy_xyyz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_xyz_1[i] * fe_0 + ta1_x_yyyyy_xyyz_0[i] * pa_y[i] - ta1_x_yyyyy_xyyz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xyzz_0[i] = 5.0 * ta1_x_yyyy_xyzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xyzz_1[i] * fe_0 + ta1_x_yyyyy_xzz_0[i] * fe_0 - ta1_x_yyyyy_xzz_1[i] * fe_0 + ta1_x_yyyyy_xyzz_0[i] * pa_y[i] - ta1_x_yyyyy_xyzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xzzz_0[i] = 5.0 * ta1_x_yyyy_xzzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xzzz_1[i] * fe_0 + ta1_x_yyyyy_xzzz_0[i] * pa_y[i] - ta1_x_yyyyy_xzzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_yyyy_0[i] = 5.0 * ta1_x_yyyy_yyyy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yyyy_1[i] * fe_0 + 4.0 * ta1_x_yyyyy_yyy_0[i] * fe_0 - 4.0 * ta1_x_yyyyy_yyy_1[i] * fe_0 + ta1_x_yyyyy_yyyy_0[i] * pa_y[i] - ta1_x_yyyyy_yyyy_1[i] * pc_y[i];

        ta1_x_yyyyyy_yyyz_0[i] = 5.0 * ta1_x_yyyy_yyyz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yyyyy_yyz_0[i] * fe_0 - 3.0 * ta1_x_yyyyy_yyz_1[i] * fe_0 + ta1_x_yyyyy_yyyz_0[i] * pa_y[i] - ta1_x_yyyyy_yyyz_1[i] * pc_y[i];

        ta1_x_yyyyyy_yyzz_0[i] = 5.0 * ta1_x_yyyy_yyzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yyyyy_yzz_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_yzz_1[i] * fe_0 + ta1_x_yyyyy_yyzz_0[i] * pa_y[i] - ta1_x_yyyyy_yyzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_yzzz_0[i] = 5.0 * ta1_x_yyyy_yzzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yzzz_1[i] * fe_0 + ta1_x_yyyyy_zzz_0[i] * fe_0 - ta1_x_yyyyy_zzz_1[i] * fe_0 + ta1_x_yyyyy_yzzz_0[i] * pa_y[i] - ta1_x_yyyyy_yzzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_zzzz_0[i] = 5.0 * ta1_x_yyyy_zzzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_zzzz_1[i] * fe_0 + ta1_x_yyyyy_zzzz_0[i] * pa_y[i] - ta1_x_yyyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 330-345 components of targeted buffer : IG

    auto ta1_x_yyyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 330);

    auto ta1_x_yyyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 331);

    auto ta1_x_yyyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 332);

    auto ta1_x_yyyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 333);

    auto ta1_x_yyyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 334);

    auto ta1_x_yyyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 335);

    auto ta1_x_yyyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 336);

    auto ta1_x_yyyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 337);

    auto ta1_x_yyyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 338);

    auto ta1_x_yyyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 339);

    auto ta1_x_yyyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 340);

    auto ta1_x_yyyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 341);

    auto ta1_x_yyyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 342);

    auto ta1_x_yyyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 343);

    auto ta1_x_yyyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 344);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyyy_xxxx_0, ta1_x_yyyyy_xxxx_1, ta1_x_yyyyy_xxxy_0, ta1_x_yyyyy_xxxy_1, ta1_x_yyyyy_xxy_0, ta1_x_yyyyy_xxy_1, ta1_x_yyyyy_xxyy_0, ta1_x_yyyyy_xxyy_1, ta1_x_yyyyy_xxyz_0, ta1_x_yyyyy_xxyz_1, ta1_x_yyyyy_xyy_0, ta1_x_yyyyy_xyy_1, ta1_x_yyyyy_xyyy_0, ta1_x_yyyyy_xyyy_1, ta1_x_yyyyy_xyyz_0, ta1_x_yyyyy_xyyz_1, ta1_x_yyyyy_xyz_0, ta1_x_yyyyy_xyz_1, ta1_x_yyyyy_xyzz_0, ta1_x_yyyyy_xyzz_1, ta1_x_yyyyy_yyy_0, ta1_x_yyyyy_yyy_1, ta1_x_yyyyy_yyyy_0, ta1_x_yyyyy_yyyy_1, ta1_x_yyyyy_yyyz_0, ta1_x_yyyyy_yyyz_1, ta1_x_yyyyy_yyz_0, ta1_x_yyyyy_yyz_1, ta1_x_yyyyy_yyzz_0, ta1_x_yyyyy_yyzz_1, ta1_x_yyyyy_yzz_0, ta1_x_yyyyy_yzz_1, ta1_x_yyyyy_yzzz_0, ta1_x_yyyyy_yzzz_1, ta1_x_yyyyyz_xxxx_0, ta1_x_yyyyyz_xxxy_0, ta1_x_yyyyyz_xxxz_0, ta1_x_yyyyyz_xxyy_0, ta1_x_yyyyyz_xxyz_0, ta1_x_yyyyyz_xxzz_0, ta1_x_yyyyyz_xyyy_0, ta1_x_yyyyyz_xyyz_0, ta1_x_yyyyyz_xyzz_0, ta1_x_yyyyyz_xzzz_0, ta1_x_yyyyyz_yyyy_0, ta1_x_yyyyyz_yyyz_0, ta1_x_yyyyyz_yyzz_0, ta1_x_yyyyyz_yzzz_0, ta1_x_yyyyyz_zzzz_0, ta1_x_yyyyz_xxxz_0, ta1_x_yyyyz_xxxz_1, ta1_x_yyyyz_xxzz_0, ta1_x_yyyyz_xxzz_1, ta1_x_yyyyz_xzzz_0, ta1_x_yyyyz_xzzz_1, ta1_x_yyyyz_zzzz_0, ta1_x_yyyyz_zzzz_1, ta1_x_yyyz_xxxz_0, ta1_x_yyyz_xxxz_1, ta1_x_yyyz_xxzz_0, ta1_x_yyyz_xxzz_1, ta1_x_yyyz_xzzz_0, ta1_x_yyyz_xzzz_1, ta1_x_yyyz_zzzz_0, ta1_x_yyyz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyz_xxxx_0[i] = ta1_x_yyyyy_xxxx_0[i] * pa_z[i] - ta1_x_yyyyy_xxxx_1[i] * pc_z[i];

        ta1_x_yyyyyz_xxxy_0[i] = ta1_x_yyyyy_xxxy_0[i] * pa_z[i] - ta1_x_yyyyy_xxxy_1[i] * pc_z[i];

        ta1_x_yyyyyz_xxxz_0[i] = 4.0 * ta1_x_yyyz_xxxz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_xxxz_1[i] * fe_0 + ta1_x_yyyyz_xxxz_0[i] * pa_y[i] - ta1_x_yyyyz_xxxz_1[i] * pc_y[i];

        ta1_x_yyyyyz_xxyy_0[i] = ta1_x_yyyyy_xxyy_0[i] * pa_z[i] - ta1_x_yyyyy_xxyy_1[i] * pc_z[i];

        ta1_x_yyyyyz_xxyz_0[i] = ta1_x_yyyyy_xxy_0[i] * fe_0 - ta1_x_yyyyy_xxy_1[i] * fe_0 + ta1_x_yyyyy_xxyz_0[i] * pa_z[i] - ta1_x_yyyyy_xxyz_1[i] * pc_z[i];

        ta1_x_yyyyyz_xxzz_0[i] = 4.0 * ta1_x_yyyz_xxzz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_xxzz_1[i] * fe_0 + ta1_x_yyyyz_xxzz_0[i] * pa_y[i] - ta1_x_yyyyz_xxzz_1[i] * pc_y[i];

        ta1_x_yyyyyz_xyyy_0[i] = ta1_x_yyyyy_xyyy_0[i] * pa_z[i] - ta1_x_yyyyy_xyyy_1[i] * pc_z[i];

        ta1_x_yyyyyz_xyyz_0[i] = ta1_x_yyyyy_xyy_0[i] * fe_0 - ta1_x_yyyyy_xyy_1[i] * fe_0 + ta1_x_yyyyy_xyyz_0[i] * pa_z[i] - ta1_x_yyyyy_xyyz_1[i] * pc_z[i];

        ta1_x_yyyyyz_xyzz_0[i] = 2.0 * ta1_x_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_xyz_1[i] * fe_0 + ta1_x_yyyyy_xyzz_0[i] * pa_z[i] - ta1_x_yyyyy_xyzz_1[i] * pc_z[i];

        ta1_x_yyyyyz_xzzz_0[i] = 4.0 * ta1_x_yyyz_xzzz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_xzzz_1[i] * fe_0 + ta1_x_yyyyz_xzzz_0[i] * pa_y[i] - ta1_x_yyyyz_xzzz_1[i] * pc_y[i];

        ta1_x_yyyyyz_yyyy_0[i] = ta1_x_yyyyy_yyyy_0[i] * pa_z[i] - ta1_x_yyyyy_yyyy_1[i] * pc_z[i];

        ta1_x_yyyyyz_yyyz_0[i] = ta1_x_yyyyy_yyy_0[i] * fe_0 - ta1_x_yyyyy_yyy_1[i] * fe_0 + ta1_x_yyyyy_yyyz_0[i] * pa_z[i] - ta1_x_yyyyy_yyyz_1[i] * pc_z[i];

        ta1_x_yyyyyz_yyzz_0[i] = 2.0 * ta1_x_yyyyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_yyz_1[i] * fe_0 + ta1_x_yyyyy_yyzz_0[i] * pa_z[i] - ta1_x_yyyyy_yyzz_1[i] * pc_z[i];

        ta1_x_yyyyyz_yzzz_0[i] = 3.0 * ta1_x_yyyyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_yyyyy_yzz_1[i] * fe_0 + ta1_x_yyyyy_yzzz_0[i] * pa_z[i] - ta1_x_yyyyy_yzzz_1[i] * pc_z[i];

        ta1_x_yyyyyz_zzzz_0[i] = 4.0 * ta1_x_yyyz_zzzz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_zzzz_1[i] * fe_0 + ta1_x_yyyyz_zzzz_0[i] * pa_y[i] - ta1_x_yyyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 345-360 components of targeted buffer : IG

    auto ta1_x_yyyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 345);

    auto ta1_x_yyyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 346);

    auto ta1_x_yyyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 347);

    auto ta1_x_yyyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 348);

    auto ta1_x_yyyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 349);

    auto ta1_x_yyyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 350);

    auto ta1_x_yyyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 351);

    auto ta1_x_yyyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 352);

    auto ta1_x_yyyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 353);

    auto ta1_x_yyyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 354);

    auto ta1_x_yyyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 355);

    auto ta1_x_yyyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 356);

    auto ta1_x_yyyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 357);

    auto ta1_x_yyyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 358);

    auto ta1_x_yyyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 359);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyy_xxxy_0, ta1_x_yyyy_xxxy_1, ta1_x_yyyy_xxyy_0, ta1_x_yyyy_xxyy_1, ta1_x_yyyy_xyyy_0, ta1_x_yyyy_xyyy_1, ta1_x_yyyy_yyyy_0, ta1_x_yyyy_yyyy_1, ta1_x_yyyyz_xxxy_0, ta1_x_yyyyz_xxxy_1, ta1_x_yyyyz_xxyy_0, ta1_x_yyyyz_xxyy_1, ta1_x_yyyyz_xyyy_0, ta1_x_yyyyz_xyyy_1, ta1_x_yyyyz_yyyy_0, ta1_x_yyyyz_yyyy_1, ta1_x_yyyyzz_xxxx_0, ta1_x_yyyyzz_xxxy_0, ta1_x_yyyyzz_xxxz_0, ta1_x_yyyyzz_xxyy_0, ta1_x_yyyyzz_xxyz_0, ta1_x_yyyyzz_xxzz_0, ta1_x_yyyyzz_xyyy_0, ta1_x_yyyyzz_xyyz_0, ta1_x_yyyyzz_xyzz_0, ta1_x_yyyyzz_xzzz_0, ta1_x_yyyyzz_yyyy_0, ta1_x_yyyyzz_yyyz_0, ta1_x_yyyyzz_yyzz_0, ta1_x_yyyyzz_yzzz_0, ta1_x_yyyyzz_zzzz_0, ta1_x_yyyzz_xxxx_0, ta1_x_yyyzz_xxxx_1, ta1_x_yyyzz_xxxz_0, ta1_x_yyyzz_xxxz_1, ta1_x_yyyzz_xxyz_0, ta1_x_yyyzz_xxyz_1, ta1_x_yyyzz_xxz_0, ta1_x_yyyzz_xxz_1, ta1_x_yyyzz_xxzz_0, ta1_x_yyyzz_xxzz_1, ta1_x_yyyzz_xyyz_0, ta1_x_yyyzz_xyyz_1, ta1_x_yyyzz_xyz_0, ta1_x_yyyzz_xyz_1, ta1_x_yyyzz_xyzz_0, ta1_x_yyyzz_xyzz_1, ta1_x_yyyzz_xzz_0, ta1_x_yyyzz_xzz_1, ta1_x_yyyzz_xzzz_0, ta1_x_yyyzz_xzzz_1, ta1_x_yyyzz_yyyz_0, ta1_x_yyyzz_yyyz_1, ta1_x_yyyzz_yyz_0, ta1_x_yyyzz_yyz_1, ta1_x_yyyzz_yyzz_0, ta1_x_yyyzz_yyzz_1, ta1_x_yyyzz_yzz_0, ta1_x_yyyzz_yzz_1, ta1_x_yyyzz_yzzz_0, ta1_x_yyyzz_yzzz_1, ta1_x_yyyzz_zzz_0, ta1_x_yyyzz_zzz_1, ta1_x_yyyzz_zzzz_0, ta1_x_yyyzz_zzzz_1, ta1_x_yyzz_xxxx_0, ta1_x_yyzz_xxxx_1, ta1_x_yyzz_xxxz_0, ta1_x_yyzz_xxxz_1, ta1_x_yyzz_xxyz_0, ta1_x_yyzz_xxyz_1, ta1_x_yyzz_xxzz_0, ta1_x_yyzz_xxzz_1, ta1_x_yyzz_xyyz_0, ta1_x_yyzz_xyyz_1, ta1_x_yyzz_xyzz_0, ta1_x_yyzz_xyzz_1, ta1_x_yyzz_xzzz_0, ta1_x_yyzz_xzzz_1, ta1_x_yyzz_yyyz_0, ta1_x_yyzz_yyyz_1, ta1_x_yyzz_yyzz_0, ta1_x_yyzz_yyzz_1, ta1_x_yyzz_yzzz_0, ta1_x_yyzz_yzzz_1, ta1_x_yyzz_zzzz_0, ta1_x_yyzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyzz_xxxx_0[i] = 3.0 * ta1_x_yyzz_xxxx_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxxx_1[i] * fe_0 + ta1_x_yyyzz_xxxx_0[i] * pa_y[i] - ta1_x_yyyzz_xxxx_1[i] * pc_y[i];

        ta1_x_yyyyzz_xxxy_0[i] = ta1_x_yyyy_xxxy_0[i] * fe_0 - ta1_x_yyyy_xxxy_1[i] * fe_0 + ta1_x_yyyyz_xxxy_0[i] * pa_z[i] - ta1_x_yyyyz_xxxy_1[i] * pc_z[i];

        ta1_x_yyyyzz_xxxz_0[i] = 3.0 * ta1_x_yyzz_xxxz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxxz_1[i] * fe_0 + ta1_x_yyyzz_xxxz_0[i] * pa_y[i] - ta1_x_yyyzz_xxxz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xxyy_0[i] = ta1_x_yyyy_xxyy_0[i] * fe_0 - ta1_x_yyyy_xxyy_1[i] * fe_0 + ta1_x_yyyyz_xxyy_0[i] * pa_z[i] - ta1_x_yyyyz_xxyy_1[i] * pc_z[i];

        ta1_x_yyyyzz_xxyz_0[i] = 3.0 * ta1_x_yyzz_xxyz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxyz_1[i] * fe_0 + ta1_x_yyyzz_xxz_0[i] * fe_0 - ta1_x_yyyzz_xxz_1[i] * fe_0 + ta1_x_yyyzz_xxyz_0[i] * pa_y[i] - ta1_x_yyyzz_xxyz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xxzz_0[i] = 3.0 * ta1_x_yyzz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxzz_1[i] * fe_0 + ta1_x_yyyzz_xxzz_0[i] * pa_y[i] - ta1_x_yyyzz_xxzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xyyy_0[i] = ta1_x_yyyy_xyyy_0[i] * fe_0 - ta1_x_yyyy_xyyy_1[i] * fe_0 + ta1_x_yyyyz_xyyy_0[i] * pa_z[i] - ta1_x_yyyyz_xyyy_1[i] * pc_z[i];

        ta1_x_yyyyzz_xyyz_0[i] = 3.0 * ta1_x_yyzz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyyzz_xyz_1[i] * fe_0 + ta1_x_yyyzz_xyyz_0[i] * pa_y[i] - ta1_x_yyyzz_xyyz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xyzz_0[i] = 3.0 * ta1_x_yyzz_xyzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xyzz_1[i] * fe_0 + ta1_x_yyyzz_xzz_0[i] * fe_0 - ta1_x_yyyzz_xzz_1[i] * fe_0 + ta1_x_yyyzz_xyzz_0[i] * pa_y[i] - ta1_x_yyyzz_xyzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xzzz_0[i] = 3.0 * ta1_x_yyzz_xzzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xzzz_1[i] * fe_0 + ta1_x_yyyzz_xzzz_0[i] * pa_y[i] - ta1_x_yyyzz_xzzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_yyyy_0[i] = ta1_x_yyyy_yyyy_0[i] * fe_0 - ta1_x_yyyy_yyyy_1[i] * fe_0 + ta1_x_yyyyz_yyyy_0[i] * pa_z[i] - ta1_x_yyyyz_yyyy_1[i] * pc_z[i];

        ta1_x_yyyyzz_yyyz_0[i] = 3.0 * ta1_x_yyzz_yyyz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yyyzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_yyyzz_yyz_1[i] * fe_0 + ta1_x_yyyzz_yyyz_0[i] * pa_y[i] - ta1_x_yyyzz_yyyz_1[i] * pc_y[i];

        ta1_x_yyyyzz_yyzz_0[i] = 3.0 * ta1_x_yyzz_yyzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yyyzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_yyyzz_yzz_1[i] * fe_0 + ta1_x_yyyzz_yyzz_0[i] * pa_y[i] - ta1_x_yyyzz_yyzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_yzzz_0[i] = 3.0 * ta1_x_yyzz_yzzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yzzz_1[i] * fe_0 + ta1_x_yyyzz_zzz_0[i] * fe_0 - ta1_x_yyyzz_zzz_1[i] * fe_0 + ta1_x_yyyzz_yzzz_0[i] * pa_y[i] - ta1_x_yyyzz_yzzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_zzzz_0[i] = 3.0 * ta1_x_yyzz_zzzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_zzzz_1[i] * fe_0 + ta1_x_yyyzz_zzzz_0[i] * pa_y[i] - ta1_x_yyyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 360-375 components of targeted buffer : IG

    auto ta1_x_yyyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 360);

    auto ta1_x_yyyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 361);

    auto ta1_x_yyyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 362);

    auto ta1_x_yyyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 363);

    auto ta1_x_yyyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 364);

    auto ta1_x_yyyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 365);

    auto ta1_x_yyyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 366);

    auto ta1_x_yyyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 367);

    auto ta1_x_yyyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 368);

    auto ta1_x_yyyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 369);

    auto ta1_x_yyyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 370);

    auto ta1_x_yyyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 371);

    auto ta1_x_yyyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 372);

    auto ta1_x_yyyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 373);

    auto ta1_x_yyyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 374);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyz_xxxy_0, ta1_x_yyyz_xxxy_1, ta1_x_yyyz_xxyy_0, ta1_x_yyyz_xxyy_1, ta1_x_yyyz_xyyy_0, ta1_x_yyyz_xyyy_1, ta1_x_yyyz_yyyy_0, ta1_x_yyyz_yyyy_1, ta1_x_yyyzz_xxxy_0, ta1_x_yyyzz_xxxy_1, ta1_x_yyyzz_xxyy_0, ta1_x_yyyzz_xxyy_1, ta1_x_yyyzz_xyyy_0, ta1_x_yyyzz_xyyy_1, ta1_x_yyyzz_yyyy_0, ta1_x_yyyzz_yyyy_1, ta1_x_yyyzzz_xxxx_0, ta1_x_yyyzzz_xxxy_0, ta1_x_yyyzzz_xxxz_0, ta1_x_yyyzzz_xxyy_0, ta1_x_yyyzzz_xxyz_0, ta1_x_yyyzzz_xxzz_0, ta1_x_yyyzzz_xyyy_0, ta1_x_yyyzzz_xyyz_0, ta1_x_yyyzzz_xyzz_0, ta1_x_yyyzzz_xzzz_0, ta1_x_yyyzzz_yyyy_0, ta1_x_yyyzzz_yyyz_0, ta1_x_yyyzzz_yyzz_0, ta1_x_yyyzzz_yzzz_0, ta1_x_yyyzzz_zzzz_0, ta1_x_yyzzz_xxxx_0, ta1_x_yyzzz_xxxx_1, ta1_x_yyzzz_xxxz_0, ta1_x_yyzzz_xxxz_1, ta1_x_yyzzz_xxyz_0, ta1_x_yyzzz_xxyz_1, ta1_x_yyzzz_xxz_0, ta1_x_yyzzz_xxz_1, ta1_x_yyzzz_xxzz_0, ta1_x_yyzzz_xxzz_1, ta1_x_yyzzz_xyyz_0, ta1_x_yyzzz_xyyz_1, ta1_x_yyzzz_xyz_0, ta1_x_yyzzz_xyz_1, ta1_x_yyzzz_xyzz_0, ta1_x_yyzzz_xyzz_1, ta1_x_yyzzz_xzz_0, ta1_x_yyzzz_xzz_1, ta1_x_yyzzz_xzzz_0, ta1_x_yyzzz_xzzz_1, ta1_x_yyzzz_yyyz_0, ta1_x_yyzzz_yyyz_1, ta1_x_yyzzz_yyz_0, ta1_x_yyzzz_yyz_1, ta1_x_yyzzz_yyzz_0, ta1_x_yyzzz_yyzz_1, ta1_x_yyzzz_yzz_0, ta1_x_yyzzz_yzz_1, ta1_x_yyzzz_yzzz_0, ta1_x_yyzzz_yzzz_1, ta1_x_yyzzz_zzz_0, ta1_x_yyzzz_zzz_1, ta1_x_yyzzz_zzzz_0, ta1_x_yyzzz_zzzz_1, ta1_x_yzzz_xxxx_0, ta1_x_yzzz_xxxx_1, ta1_x_yzzz_xxxz_0, ta1_x_yzzz_xxxz_1, ta1_x_yzzz_xxyz_0, ta1_x_yzzz_xxyz_1, ta1_x_yzzz_xxzz_0, ta1_x_yzzz_xxzz_1, ta1_x_yzzz_xyyz_0, ta1_x_yzzz_xyyz_1, ta1_x_yzzz_xyzz_0, ta1_x_yzzz_xyzz_1, ta1_x_yzzz_xzzz_0, ta1_x_yzzz_xzzz_1, ta1_x_yzzz_yyyz_0, ta1_x_yzzz_yyyz_1, ta1_x_yzzz_yyzz_0, ta1_x_yzzz_yyzz_1, ta1_x_yzzz_yzzz_0, ta1_x_yzzz_yzzz_1, ta1_x_yzzz_zzzz_0, ta1_x_yzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzzz_xxxx_0[i] = 2.0 * ta1_x_yzzz_xxxx_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xxxx_1[i] * fe_0 + ta1_x_yyzzz_xxxx_0[i] * pa_y[i] - ta1_x_yyzzz_xxxx_1[i] * pc_y[i];

        ta1_x_yyyzzz_xxxy_0[i] = 2.0 * ta1_x_yyyz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_xxxy_1[i] * fe_0 + ta1_x_yyyzz_xxxy_0[i] * pa_z[i] - ta1_x_yyyzz_xxxy_1[i] * pc_z[i];

        ta1_x_yyyzzz_xxxz_0[i] = 2.0 * ta1_x_yzzz_xxxz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xxxz_1[i] * fe_0 + ta1_x_yyzzz_xxxz_0[i] * pa_y[i] - ta1_x_yyzzz_xxxz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xxyy_0[i] = 2.0 * ta1_x_yyyz_xxyy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_xxyy_1[i] * fe_0 + ta1_x_yyyzz_xxyy_0[i] * pa_z[i] - ta1_x_yyyzz_xxyy_1[i] * pc_z[i];

        ta1_x_yyyzzz_xxyz_0[i] = 2.0 * ta1_x_yzzz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xxyz_1[i] * fe_0 + ta1_x_yyzzz_xxz_0[i] * fe_0 - ta1_x_yyzzz_xxz_1[i] * fe_0 + ta1_x_yyzzz_xxyz_0[i] * pa_y[i] - ta1_x_yyzzz_xxyz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xxzz_0[i] = 2.0 * ta1_x_yzzz_xxzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xxzz_1[i] * fe_0 + ta1_x_yyzzz_xxzz_0[i] * pa_y[i] - ta1_x_yyzzz_xxzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xyyy_0[i] = 2.0 * ta1_x_yyyz_xyyy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_xyyy_1[i] * fe_0 + ta1_x_yyyzz_xyyy_0[i] * pa_z[i] - ta1_x_yyyzz_xyyy_1[i] * pc_z[i];

        ta1_x_yyyzzz_xyyz_0[i] = 2.0 * ta1_x_yzzz_xyyz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yyzzz_xyz_1[i] * fe_0 + ta1_x_yyzzz_xyyz_0[i] * pa_y[i] - ta1_x_yyzzz_xyyz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xyzz_0[i] = 2.0 * ta1_x_yzzz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xyzz_1[i] * fe_0 + ta1_x_yyzzz_xzz_0[i] * fe_0 - ta1_x_yyzzz_xzz_1[i] * fe_0 + ta1_x_yyzzz_xyzz_0[i] * pa_y[i] - ta1_x_yyzzz_xyzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xzzz_0[i] = 2.0 * ta1_x_yzzz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xzzz_1[i] * fe_0 + ta1_x_yyzzz_xzzz_0[i] * pa_y[i] - ta1_x_yyzzz_xzzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_yyyy_0[i] = 2.0 * ta1_x_yyyz_yyyy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_yyyy_1[i] * fe_0 + ta1_x_yyyzz_yyyy_0[i] * pa_z[i] - ta1_x_yyyzz_yyyy_1[i] * pc_z[i];

        ta1_x_yyyzzz_yyyz_0[i] = 2.0 * ta1_x_yzzz_yyyz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yyzzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_yyzzz_yyz_1[i] * fe_0 + ta1_x_yyzzz_yyyz_0[i] * pa_y[i] - ta1_x_yyzzz_yyyz_1[i] * pc_y[i];

        ta1_x_yyyzzz_yyzz_0[i] = 2.0 * ta1_x_yzzz_yyzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yyzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_yyzzz_yzz_1[i] * fe_0 + ta1_x_yyzzz_yyzz_0[i] * pa_y[i] - ta1_x_yyzzz_yyzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_yzzz_0[i] = 2.0 * ta1_x_yzzz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yzzz_1[i] * fe_0 + ta1_x_yyzzz_zzz_0[i] * fe_0 - ta1_x_yyzzz_zzz_1[i] * fe_0 + ta1_x_yyzzz_yzzz_0[i] * pa_y[i] - ta1_x_yyzzz_yzzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_zzzz_0[i] = 2.0 * ta1_x_yzzz_zzzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_zzzz_1[i] * fe_0 + ta1_x_yyzzz_zzzz_0[i] * pa_y[i] - ta1_x_yyzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 375-390 components of targeted buffer : IG

    auto ta1_x_yyzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 375);

    auto ta1_x_yyzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 376);

    auto ta1_x_yyzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 377);

    auto ta1_x_yyzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 378);

    auto ta1_x_yyzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 379);

    auto ta1_x_yyzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 380);

    auto ta1_x_yyzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 381);

    auto ta1_x_yyzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 382);

    auto ta1_x_yyzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 383);

    auto ta1_x_yyzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 384);

    auto ta1_x_yyzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 385);

    auto ta1_x_yyzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 386);

    auto ta1_x_yyzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 387);

    auto ta1_x_yyzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 388);

    auto ta1_x_yyzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 389);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyzz_xxxy_0, ta1_x_yyzz_xxxy_1, ta1_x_yyzz_xxyy_0, ta1_x_yyzz_xxyy_1, ta1_x_yyzz_xyyy_0, ta1_x_yyzz_xyyy_1, ta1_x_yyzz_yyyy_0, ta1_x_yyzz_yyyy_1, ta1_x_yyzzz_xxxy_0, ta1_x_yyzzz_xxxy_1, ta1_x_yyzzz_xxyy_0, ta1_x_yyzzz_xxyy_1, ta1_x_yyzzz_xyyy_0, ta1_x_yyzzz_xyyy_1, ta1_x_yyzzz_yyyy_0, ta1_x_yyzzz_yyyy_1, ta1_x_yyzzzz_xxxx_0, ta1_x_yyzzzz_xxxy_0, ta1_x_yyzzzz_xxxz_0, ta1_x_yyzzzz_xxyy_0, ta1_x_yyzzzz_xxyz_0, ta1_x_yyzzzz_xxzz_0, ta1_x_yyzzzz_xyyy_0, ta1_x_yyzzzz_xyyz_0, ta1_x_yyzzzz_xyzz_0, ta1_x_yyzzzz_xzzz_0, ta1_x_yyzzzz_yyyy_0, ta1_x_yyzzzz_yyyz_0, ta1_x_yyzzzz_yyzz_0, ta1_x_yyzzzz_yzzz_0, ta1_x_yyzzzz_zzzz_0, ta1_x_yzzzz_xxxx_0, ta1_x_yzzzz_xxxx_1, ta1_x_yzzzz_xxxz_0, ta1_x_yzzzz_xxxz_1, ta1_x_yzzzz_xxyz_0, ta1_x_yzzzz_xxyz_1, ta1_x_yzzzz_xxz_0, ta1_x_yzzzz_xxz_1, ta1_x_yzzzz_xxzz_0, ta1_x_yzzzz_xxzz_1, ta1_x_yzzzz_xyyz_0, ta1_x_yzzzz_xyyz_1, ta1_x_yzzzz_xyz_0, ta1_x_yzzzz_xyz_1, ta1_x_yzzzz_xyzz_0, ta1_x_yzzzz_xyzz_1, ta1_x_yzzzz_xzz_0, ta1_x_yzzzz_xzz_1, ta1_x_yzzzz_xzzz_0, ta1_x_yzzzz_xzzz_1, ta1_x_yzzzz_yyyz_0, ta1_x_yzzzz_yyyz_1, ta1_x_yzzzz_yyz_0, ta1_x_yzzzz_yyz_1, ta1_x_yzzzz_yyzz_0, ta1_x_yzzzz_yyzz_1, ta1_x_yzzzz_yzz_0, ta1_x_yzzzz_yzz_1, ta1_x_yzzzz_yzzz_0, ta1_x_yzzzz_yzzz_1, ta1_x_yzzzz_zzz_0, ta1_x_yzzzz_zzz_1, ta1_x_yzzzz_zzzz_0, ta1_x_yzzzz_zzzz_1, ta1_x_zzzz_xxxx_0, ta1_x_zzzz_xxxx_1, ta1_x_zzzz_xxxz_0, ta1_x_zzzz_xxxz_1, ta1_x_zzzz_xxyz_0, ta1_x_zzzz_xxyz_1, ta1_x_zzzz_xxzz_0, ta1_x_zzzz_xxzz_1, ta1_x_zzzz_xyyz_0, ta1_x_zzzz_xyyz_1, ta1_x_zzzz_xyzz_0, ta1_x_zzzz_xyzz_1, ta1_x_zzzz_xzzz_0, ta1_x_zzzz_xzzz_1, ta1_x_zzzz_yyyz_0, ta1_x_zzzz_yyyz_1, ta1_x_zzzz_yyzz_0, ta1_x_zzzz_yyzz_1, ta1_x_zzzz_yzzz_0, ta1_x_zzzz_yzzz_1, ta1_x_zzzz_zzzz_0, ta1_x_zzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzzz_xxxx_0[i] = ta1_x_zzzz_xxxx_0[i] * fe_0 - ta1_x_zzzz_xxxx_1[i] * fe_0 + ta1_x_yzzzz_xxxx_0[i] * pa_y[i] - ta1_x_yzzzz_xxxx_1[i] * pc_y[i];

        ta1_x_yyzzzz_xxxy_0[i] = 3.0 * ta1_x_yyzz_xxxy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxxy_1[i] * fe_0 + ta1_x_yyzzz_xxxy_0[i] * pa_z[i] - ta1_x_yyzzz_xxxy_1[i] * pc_z[i];

        ta1_x_yyzzzz_xxxz_0[i] = ta1_x_zzzz_xxxz_0[i] * fe_0 - ta1_x_zzzz_xxxz_1[i] * fe_0 + ta1_x_yzzzz_xxxz_0[i] * pa_y[i] - ta1_x_yzzzz_xxxz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xxyy_0[i] = 3.0 * ta1_x_yyzz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxyy_1[i] * fe_0 + ta1_x_yyzzz_xxyy_0[i] * pa_z[i] - ta1_x_yyzzz_xxyy_1[i] * pc_z[i];

        ta1_x_yyzzzz_xxyz_0[i] = ta1_x_zzzz_xxyz_0[i] * fe_0 - ta1_x_zzzz_xxyz_1[i] * fe_0 + ta1_x_yzzzz_xxz_0[i] * fe_0 - ta1_x_yzzzz_xxz_1[i] * fe_0 + ta1_x_yzzzz_xxyz_0[i] * pa_y[i] - ta1_x_yzzzz_xxyz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xxzz_0[i] = ta1_x_zzzz_xxzz_0[i] * fe_0 - ta1_x_zzzz_xxzz_1[i] * fe_0 + ta1_x_yzzzz_xxzz_0[i] * pa_y[i] - ta1_x_yzzzz_xxzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xyyy_0[i] = 3.0 * ta1_x_yyzz_xyyy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xyyy_1[i] * fe_0 + ta1_x_yyzzz_xyyy_0[i] * pa_z[i] - ta1_x_yyzzz_xyyy_1[i] * pc_z[i];

        ta1_x_yyzzzz_xyyz_0[i] = ta1_x_zzzz_xyyz_0[i] * fe_0 - ta1_x_zzzz_xyyz_1[i] * fe_0 + 2.0 * ta1_x_yzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yzzzz_xyz_1[i] * fe_0 + ta1_x_yzzzz_xyyz_0[i] * pa_y[i] - ta1_x_yzzzz_xyyz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xyzz_0[i] = ta1_x_zzzz_xyzz_0[i] * fe_0 - ta1_x_zzzz_xyzz_1[i] * fe_0 + ta1_x_yzzzz_xzz_0[i] * fe_0 - ta1_x_yzzzz_xzz_1[i] * fe_0 + ta1_x_yzzzz_xyzz_0[i] * pa_y[i] - ta1_x_yzzzz_xyzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xzzz_0[i] = ta1_x_zzzz_xzzz_0[i] * fe_0 - ta1_x_zzzz_xzzz_1[i] * fe_0 + ta1_x_yzzzz_xzzz_0[i] * pa_y[i] - ta1_x_yzzzz_xzzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_yyyy_0[i] = 3.0 * ta1_x_yyzz_yyyy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yyyy_1[i] * fe_0 + ta1_x_yyzzz_yyyy_0[i] * pa_z[i] - ta1_x_yyzzz_yyyy_1[i] * pc_z[i];

        ta1_x_yyzzzz_yyyz_0[i] = ta1_x_zzzz_yyyz_0[i] * fe_0 - ta1_x_zzzz_yyyz_1[i] * fe_0 + 3.0 * ta1_x_yzzzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_yzzzz_yyz_1[i] * fe_0 + ta1_x_yzzzz_yyyz_0[i] * pa_y[i] - ta1_x_yzzzz_yyyz_1[i] * pc_y[i];

        ta1_x_yyzzzz_yyzz_0[i] = ta1_x_zzzz_yyzz_0[i] * fe_0 - ta1_x_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_yzzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_yzzzz_yzz_1[i] * fe_0 + ta1_x_yzzzz_yyzz_0[i] * pa_y[i] - ta1_x_yzzzz_yyzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_yzzz_0[i] = ta1_x_zzzz_yzzz_0[i] * fe_0 - ta1_x_zzzz_yzzz_1[i] * fe_0 + ta1_x_yzzzz_zzz_0[i] * fe_0 - ta1_x_yzzzz_zzz_1[i] * fe_0 + ta1_x_yzzzz_yzzz_0[i] * pa_y[i] - ta1_x_yzzzz_yzzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_zzzz_0[i] = ta1_x_zzzz_zzzz_0[i] * fe_0 - ta1_x_zzzz_zzzz_1[i] * fe_0 + ta1_x_yzzzz_zzzz_0[i] * pa_y[i] - ta1_x_yzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 390-405 components of targeted buffer : IG

    auto ta1_x_yzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 390);

    auto ta1_x_yzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 391);

    auto ta1_x_yzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 392);

    auto ta1_x_yzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 393);

    auto ta1_x_yzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 394);

    auto ta1_x_yzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 395);

    auto ta1_x_yzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 396);

    auto ta1_x_yzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 397);

    auto ta1_x_yzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 398);

    auto ta1_x_yzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 399);

    auto ta1_x_yzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 400);

    auto ta1_x_yzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 401);

    auto ta1_x_yzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 402);

    auto ta1_x_yzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 403);

    auto ta1_x_yzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 404);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yzzzzz_xxxx_0, ta1_x_yzzzzz_xxxy_0, ta1_x_yzzzzz_xxxz_0, ta1_x_yzzzzz_xxyy_0, ta1_x_yzzzzz_xxyz_0, ta1_x_yzzzzz_xxzz_0, ta1_x_yzzzzz_xyyy_0, ta1_x_yzzzzz_xyyz_0, ta1_x_yzzzzz_xyzz_0, ta1_x_yzzzzz_xzzz_0, ta1_x_yzzzzz_yyyy_0, ta1_x_yzzzzz_yyyz_0, ta1_x_yzzzzz_yyzz_0, ta1_x_yzzzzz_yzzz_0, ta1_x_yzzzzz_zzzz_0, ta1_x_zzzzz_xxx_0, ta1_x_zzzzz_xxx_1, ta1_x_zzzzz_xxxx_0, ta1_x_zzzzz_xxxx_1, ta1_x_zzzzz_xxxy_0, ta1_x_zzzzz_xxxy_1, ta1_x_zzzzz_xxxz_0, ta1_x_zzzzz_xxxz_1, ta1_x_zzzzz_xxy_0, ta1_x_zzzzz_xxy_1, ta1_x_zzzzz_xxyy_0, ta1_x_zzzzz_xxyy_1, ta1_x_zzzzz_xxyz_0, ta1_x_zzzzz_xxyz_1, ta1_x_zzzzz_xxz_0, ta1_x_zzzzz_xxz_1, ta1_x_zzzzz_xxzz_0, ta1_x_zzzzz_xxzz_1, ta1_x_zzzzz_xyy_0, ta1_x_zzzzz_xyy_1, ta1_x_zzzzz_xyyy_0, ta1_x_zzzzz_xyyy_1, ta1_x_zzzzz_xyyz_0, ta1_x_zzzzz_xyyz_1, ta1_x_zzzzz_xyz_0, ta1_x_zzzzz_xyz_1, ta1_x_zzzzz_xyzz_0, ta1_x_zzzzz_xyzz_1, ta1_x_zzzzz_xzz_0, ta1_x_zzzzz_xzz_1, ta1_x_zzzzz_xzzz_0, ta1_x_zzzzz_xzzz_1, ta1_x_zzzzz_yyy_0, ta1_x_zzzzz_yyy_1, ta1_x_zzzzz_yyyy_0, ta1_x_zzzzz_yyyy_1, ta1_x_zzzzz_yyyz_0, ta1_x_zzzzz_yyyz_1, ta1_x_zzzzz_yyz_0, ta1_x_zzzzz_yyz_1, ta1_x_zzzzz_yyzz_0, ta1_x_zzzzz_yyzz_1, ta1_x_zzzzz_yzz_0, ta1_x_zzzzz_yzz_1, ta1_x_zzzzz_yzzz_0, ta1_x_zzzzz_yzzz_1, ta1_x_zzzzz_zzz_0, ta1_x_zzzzz_zzz_1, ta1_x_zzzzz_zzzz_0, ta1_x_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzzz_xxxx_0[i] = ta1_x_zzzzz_xxxx_0[i] * pa_y[i] - ta1_x_zzzzz_xxxx_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxxy_0[i] = ta1_x_zzzzz_xxx_0[i] * fe_0 - ta1_x_zzzzz_xxx_1[i] * fe_0 + ta1_x_zzzzz_xxxy_0[i] * pa_y[i] - ta1_x_zzzzz_xxxy_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxxz_0[i] = ta1_x_zzzzz_xxxz_0[i] * pa_y[i] - ta1_x_zzzzz_xxxz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxyy_0[i] = 2.0 * ta1_x_zzzzz_xxy_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xxy_1[i] * fe_0 + ta1_x_zzzzz_xxyy_0[i] * pa_y[i] - ta1_x_zzzzz_xxyy_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxyz_0[i] = ta1_x_zzzzz_xxz_0[i] * fe_0 - ta1_x_zzzzz_xxz_1[i] * fe_0 + ta1_x_zzzzz_xxyz_0[i] * pa_y[i] - ta1_x_zzzzz_xxyz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxzz_0[i] = ta1_x_zzzzz_xxzz_0[i] * pa_y[i] - ta1_x_zzzzz_xxzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xyyy_0[i] = 3.0 * ta1_x_zzzzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_zzzzz_xyy_1[i] * fe_0 + ta1_x_zzzzz_xyyy_0[i] * pa_y[i] - ta1_x_zzzzz_xyyy_1[i] * pc_y[i];

        ta1_x_yzzzzz_xyyz_0[i] = 2.0 * ta1_x_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xyz_1[i] * fe_0 + ta1_x_zzzzz_xyyz_0[i] * pa_y[i] - ta1_x_zzzzz_xyyz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xyzz_0[i] = ta1_x_zzzzz_xzz_0[i] * fe_0 - ta1_x_zzzzz_xzz_1[i] * fe_0 + ta1_x_zzzzz_xyzz_0[i] * pa_y[i] - ta1_x_zzzzz_xyzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xzzz_0[i] = ta1_x_zzzzz_xzzz_0[i] * pa_y[i] - ta1_x_zzzzz_xzzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_yyyy_0[i] = 4.0 * ta1_x_zzzzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_zzzzz_yyy_1[i] * fe_0 + ta1_x_zzzzz_yyyy_0[i] * pa_y[i] - ta1_x_zzzzz_yyyy_1[i] * pc_y[i];

        ta1_x_yzzzzz_yyyz_0[i] = 3.0 * ta1_x_zzzzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_zzzzz_yyz_1[i] * fe_0 + ta1_x_zzzzz_yyyz_0[i] * pa_y[i] - ta1_x_zzzzz_yyyz_1[i] * pc_y[i];

        ta1_x_yzzzzz_yyzz_0[i] = 2.0 * ta1_x_zzzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_yzz_1[i] * fe_0 + ta1_x_zzzzz_yyzz_0[i] * pa_y[i] - ta1_x_zzzzz_yyzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_yzzz_0[i] = ta1_x_zzzzz_zzz_0[i] * fe_0 - ta1_x_zzzzz_zzz_1[i] * fe_0 + ta1_x_zzzzz_yzzz_0[i] * pa_y[i] - ta1_x_zzzzz_yzzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_zzzz_0[i] = ta1_x_zzzzz_zzzz_0[i] * pa_y[i] - ta1_x_zzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 405-420 components of targeted buffer : IG

    auto ta1_x_zzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 405);

    auto ta1_x_zzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 406);

    auto ta1_x_zzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 407);

    auto ta1_x_zzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 408);

    auto ta1_x_zzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 409);

    auto ta1_x_zzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 410);

    auto ta1_x_zzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 411);

    auto ta1_x_zzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 412);

    auto ta1_x_zzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 413);

    auto ta1_x_zzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 414);

    auto ta1_x_zzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 415);

    auto ta1_x_zzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 416);

    auto ta1_x_zzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 417);

    auto ta1_x_zzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 418);

    auto ta1_x_zzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 419);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zzzz_xxxx_0, ta1_x_zzzz_xxxx_1, ta1_x_zzzz_xxxy_0, ta1_x_zzzz_xxxy_1, ta1_x_zzzz_xxxz_0, ta1_x_zzzz_xxxz_1, ta1_x_zzzz_xxyy_0, ta1_x_zzzz_xxyy_1, ta1_x_zzzz_xxyz_0, ta1_x_zzzz_xxyz_1, ta1_x_zzzz_xxzz_0, ta1_x_zzzz_xxzz_1, ta1_x_zzzz_xyyy_0, ta1_x_zzzz_xyyy_1, ta1_x_zzzz_xyyz_0, ta1_x_zzzz_xyyz_1, ta1_x_zzzz_xyzz_0, ta1_x_zzzz_xyzz_1, ta1_x_zzzz_xzzz_0, ta1_x_zzzz_xzzz_1, ta1_x_zzzz_yyyy_0, ta1_x_zzzz_yyyy_1, ta1_x_zzzz_yyyz_0, ta1_x_zzzz_yyyz_1, ta1_x_zzzz_yyzz_0, ta1_x_zzzz_yyzz_1, ta1_x_zzzz_yzzz_0, ta1_x_zzzz_yzzz_1, ta1_x_zzzz_zzzz_0, ta1_x_zzzz_zzzz_1, ta1_x_zzzzz_xxx_0, ta1_x_zzzzz_xxx_1, ta1_x_zzzzz_xxxx_0, ta1_x_zzzzz_xxxx_1, ta1_x_zzzzz_xxxy_0, ta1_x_zzzzz_xxxy_1, ta1_x_zzzzz_xxxz_0, ta1_x_zzzzz_xxxz_1, ta1_x_zzzzz_xxy_0, ta1_x_zzzzz_xxy_1, ta1_x_zzzzz_xxyy_0, ta1_x_zzzzz_xxyy_1, ta1_x_zzzzz_xxyz_0, ta1_x_zzzzz_xxyz_1, ta1_x_zzzzz_xxz_0, ta1_x_zzzzz_xxz_1, ta1_x_zzzzz_xxzz_0, ta1_x_zzzzz_xxzz_1, ta1_x_zzzzz_xyy_0, ta1_x_zzzzz_xyy_1, ta1_x_zzzzz_xyyy_0, ta1_x_zzzzz_xyyy_1, ta1_x_zzzzz_xyyz_0, ta1_x_zzzzz_xyyz_1, ta1_x_zzzzz_xyz_0, ta1_x_zzzzz_xyz_1, ta1_x_zzzzz_xyzz_0, ta1_x_zzzzz_xyzz_1, ta1_x_zzzzz_xzz_0, ta1_x_zzzzz_xzz_1, ta1_x_zzzzz_xzzz_0, ta1_x_zzzzz_xzzz_1, ta1_x_zzzzz_yyy_0, ta1_x_zzzzz_yyy_1, ta1_x_zzzzz_yyyy_0, ta1_x_zzzzz_yyyy_1, ta1_x_zzzzz_yyyz_0, ta1_x_zzzzz_yyyz_1, ta1_x_zzzzz_yyz_0, ta1_x_zzzzz_yyz_1, ta1_x_zzzzz_yyzz_0, ta1_x_zzzzz_yyzz_1, ta1_x_zzzzz_yzz_0, ta1_x_zzzzz_yzz_1, ta1_x_zzzzz_yzzz_0, ta1_x_zzzzz_yzzz_1, ta1_x_zzzzz_zzz_0, ta1_x_zzzzz_zzz_1, ta1_x_zzzzz_zzzz_0, ta1_x_zzzzz_zzzz_1, ta1_x_zzzzzz_xxxx_0, ta1_x_zzzzzz_xxxy_0, ta1_x_zzzzzz_xxxz_0, ta1_x_zzzzzz_xxyy_0, ta1_x_zzzzzz_xxyz_0, ta1_x_zzzzzz_xxzz_0, ta1_x_zzzzzz_xyyy_0, ta1_x_zzzzzz_xyyz_0, ta1_x_zzzzzz_xyzz_0, ta1_x_zzzzzz_xzzz_0, ta1_x_zzzzzz_yyyy_0, ta1_x_zzzzzz_yyyz_0, ta1_x_zzzzzz_yyzz_0, ta1_x_zzzzzz_yzzz_0, ta1_x_zzzzzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzzz_xxxx_0[i] = 5.0 * ta1_x_zzzz_xxxx_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxxx_1[i] * fe_0 + ta1_x_zzzzz_xxxx_0[i] * pa_z[i] - ta1_x_zzzzz_xxxx_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxxy_0[i] = 5.0 * ta1_x_zzzz_xxxy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxxy_1[i] * fe_0 + ta1_x_zzzzz_xxxy_0[i] * pa_z[i] - ta1_x_zzzzz_xxxy_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxxz_0[i] = 5.0 * ta1_x_zzzz_xxxz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxxz_1[i] * fe_0 + ta1_x_zzzzz_xxx_0[i] * fe_0 - ta1_x_zzzzz_xxx_1[i] * fe_0 + ta1_x_zzzzz_xxxz_0[i] * pa_z[i] - ta1_x_zzzzz_xxxz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxyy_0[i] = 5.0 * ta1_x_zzzz_xxyy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxyy_1[i] * fe_0 + ta1_x_zzzzz_xxyy_0[i] * pa_z[i] - ta1_x_zzzzz_xxyy_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxyz_0[i] = 5.0 * ta1_x_zzzz_xxyz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxyz_1[i] * fe_0 + ta1_x_zzzzz_xxy_0[i] * fe_0 - ta1_x_zzzzz_xxy_1[i] * fe_0 + ta1_x_zzzzz_xxyz_0[i] * pa_z[i] - ta1_x_zzzzz_xxyz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxzz_0[i] = 5.0 * ta1_x_zzzz_xxzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_x_zzzzz_xxz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xxz_1[i] * fe_0 + ta1_x_zzzzz_xxzz_0[i] * pa_z[i] - ta1_x_zzzzz_xxzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xyyy_0[i] = 5.0 * ta1_x_zzzz_xyyy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xyyy_1[i] * fe_0 + ta1_x_zzzzz_xyyy_0[i] * pa_z[i] - ta1_x_zzzzz_xyyy_1[i] * pc_z[i];

        ta1_x_zzzzzz_xyyz_0[i] = 5.0 * ta1_x_zzzz_xyyz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xyyz_1[i] * fe_0 + ta1_x_zzzzz_xyy_0[i] * fe_0 - ta1_x_zzzzz_xyy_1[i] * fe_0 + ta1_x_zzzzz_xyyz_0[i] * pa_z[i] - ta1_x_zzzzz_xyyz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xyzz_0[i] = 5.0 * ta1_x_zzzz_xyzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xyzz_1[i] * fe_0 + 2.0 * ta1_x_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xyz_1[i] * fe_0 + ta1_x_zzzzz_xyzz_0[i] * pa_z[i] - ta1_x_zzzzz_xyzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xzzz_0[i] = 5.0 * ta1_x_zzzz_xzzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xzzz_1[i] * fe_0 + 3.0 * ta1_x_zzzzz_xzz_0[i] * fe_0 - 3.0 * ta1_x_zzzzz_xzz_1[i] * fe_0 + ta1_x_zzzzz_xzzz_0[i] * pa_z[i] - ta1_x_zzzzz_xzzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_yyyy_0[i] = 5.0 * ta1_x_zzzz_yyyy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yyyy_1[i] * fe_0 + ta1_x_zzzzz_yyyy_0[i] * pa_z[i] - ta1_x_zzzzz_yyyy_1[i] * pc_z[i];

        ta1_x_zzzzzz_yyyz_0[i] = 5.0 * ta1_x_zzzz_yyyz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yyyz_1[i] * fe_0 + ta1_x_zzzzz_yyy_0[i] * fe_0 - ta1_x_zzzzz_yyy_1[i] * fe_0 + ta1_x_zzzzz_yyyz_0[i] * pa_z[i] - ta1_x_zzzzz_yyyz_1[i] * pc_z[i];

        ta1_x_zzzzzz_yyzz_0[i] = 5.0 * ta1_x_zzzz_yyzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_zzzzz_yyz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_yyz_1[i] * fe_0 + ta1_x_zzzzz_yyzz_0[i] * pa_z[i] - ta1_x_zzzzz_yyzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_yzzz_0[i] = 5.0 * ta1_x_zzzz_yzzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yzzz_1[i] * fe_0 + 3.0 * ta1_x_zzzzz_yzz_0[i] * fe_0 - 3.0 * ta1_x_zzzzz_yzz_1[i] * fe_0 + ta1_x_zzzzz_yzzz_0[i] * pa_z[i] - ta1_x_zzzzz_yzzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_zzzz_0[i] = 5.0 * ta1_x_zzzz_zzzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_zzzz_1[i] * fe_0 + 4.0 * ta1_x_zzzzz_zzz_0[i] * fe_0 - 4.0 * ta1_x_zzzzz_zzz_1[i] * fe_0 + ta1_x_zzzzz_zzzz_0[i] * pa_z[i] - ta1_x_zzzzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 420-435 components of targeted buffer : IG

    auto ta1_y_xxxxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 420);

    auto ta1_y_xxxxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 421);

    auto ta1_y_xxxxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 422);

    auto ta1_y_xxxxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 423);

    auto ta1_y_xxxxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 424);

    auto ta1_y_xxxxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 425);

    auto ta1_y_xxxxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 426);

    auto ta1_y_xxxxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 427);

    auto ta1_y_xxxxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 428);

    auto ta1_y_xxxxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 429);

    auto ta1_y_xxxxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 430);

    auto ta1_y_xxxxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 431);

    auto ta1_y_xxxxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 432);

    auto ta1_y_xxxxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 433);

    auto ta1_y_xxxxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 434);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xxxx_xxxx_0, ta1_y_xxxx_xxxx_1, ta1_y_xxxx_xxxy_0, ta1_y_xxxx_xxxy_1, ta1_y_xxxx_xxxz_0, ta1_y_xxxx_xxxz_1, ta1_y_xxxx_xxyy_0, ta1_y_xxxx_xxyy_1, ta1_y_xxxx_xxyz_0, ta1_y_xxxx_xxyz_1, ta1_y_xxxx_xxzz_0, ta1_y_xxxx_xxzz_1, ta1_y_xxxx_xyyy_0, ta1_y_xxxx_xyyy_1, ta1_y_xxxx_xyyz_0, ta1_y_xxxx_xyyz_1, ta1_y_xxxx_xyzz_0, ta1_y_xxxx_xyzz_1, ta1_y_xxxx_xzzz_0, ta1_y_xxxx_xzzz_1, ta1_y_xxxx_yyyy_0, ta1_y_xxxx_yyyy_1, ta1_y_xxxx_yyyz_0, ta1_y_xxxx_yyyz_1, ta1_y_xxxx_yyzz_0, ta1_y_xxxx_yyzz_1, ta1_y_xxxx_yzzz_0, ta1_y_xxxx_yzzz_1, ta1_y_xxxx_zzzz_0, ta1_y_xxxx_zzzz_1, ta1_y_xxxxx_xxx_0, ta1_y_xxxxx_xxx_1, ta1_y_xxxxx_xxxx_0, ta1_y_xxxxx_xxxx_1, ta1_y_xxxxx_xxxy_0, ta1_y_xxxxx_xxxy_1, ta1_y_xxxxx_xxxz_0, ta1_y_xxxxx_xxxz_1, ta1_y_xxxxx_xxy_0, ta1_y_xxxxx_xxy_1, ta1_y_xxxxx_xxyy_0, ta1_y_xxxxx_xxyy_1, ta1_y_xxxxx_xxyz_0, ta1_y_xxxxx_xxyz_1, ta1_y_xxxxx_xxz_0, ta1_y_xxxxx_xxz_1, ta1_y_xxxxx_xxzz_0, ta1_y_xxxxx_xxzz_1, ta1_y_xxxxx_xyy_0, ta1_y_xxxxx_xyy_1, ta1_y_xxxxx_xyyy_0, ta1_y_xxxxx_xyyy_1, ta1_y_xxxxx_xyyz_0, ta1_y_xxxxx_xyyz_1, ta1_y_xxxxx_xyz_0, ta1_y_xxxxx_xyz_1, ta1_y_xxxxx_xyzz_0, ta1_y_xxxxx_xyzz_1, ta1_y_xxxxx_xzz_0, ta1_y_xxxxx_xzz_1, ta1_y_xxxxx_xzzz_0, ta1_y_xxxxx_xzzz_1, ta1_y_xxxxx_yyy_0, ta1_y_xxxxx_yyy_1, ta1_y_xxxxx_yyyy_0, ta1_y_xxxxx_yyyy_1, ta1_y_xxxxx_yyyz_0, ta1_y_xxxxx_yyyz_1, ta1_y_xxxxx_yyz_0, ta1_y_xxxxx_yyz_1, ta1_y_xxxxx_yyzz_0, ta1_y_xxxxx_yyzz_1, ta1_y_xxxxx_yzz_0, ta1_y_xxxxx_yzz_1, ta1_y_xxxxx_yzzz_0, ta1_y_xxxxx_yzzz_1, ta1_y_xxxxx_zzz_0, ta1_y_xxxxx_zzz_1, ta1_y_xxxxx_zzzz_0, ta1_y_xxxxx_zzzz_1, ta1_y_xxxxxx_xxxx_0, ta1_y_xxxxxx_xxxy_0, ta1_y_xxxxxx_xxxz_0, ta1_y_xxxxxx_xxyy_0, ta1_y_xxxxxx_xxyz_0, ta1_y_xxxxxx_xxzz_0, ta1_y_xxxxxx_xyyy_0, ta1_y_xxxxxx_xyyz_0, ta1_y_xxxxxx_xyzz_0, ta1_y_xxxxxx_xzzz_0, ta1_y_xxxxxx_yyyy_0, ta1_y_xxxxxx_yyyz_0, ta1_y_xxxxxx_yyzz_0, ta1_y_xxxxxx_yzzz_0, ta1_y_xxxxxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxx_xxxx_0[i] = 5.0 * ta1_y_xxxx_xxxx_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxxx_1[i] * fe_0 + 4.0 * ta1_y_xxxxx_xxx_0[i] * fe_0 - 4.0 * ta1_y_xxxxx_xxx_1[i] * fe_0 + ta1_y_xxxxx_xxxx_0[i] * pa_x[i] - ta1_y_xxxxx_xxxx_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxxy_0[i] = 5.0 * ta1_y_xxxx_xxxy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xxxxx_xxy_0[i] * fe_0 - 3.0 * ta1_y_xxxxx_xxy_1[i] * fe_0 + ta1_y_xxxxx_xxxy_0[i] * pa_x[i] - ta1_y_xxxxx_xxxy_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxxz_0[i] = 5.0 * ta1_y_xxxx_xxxz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xxxxx_xxz_0[i] * fe_0 - 3.0 * ta1_y_xxxxx_xxz_1[i] * fe_0 + ta1_y_xxxxx_xxxz_0[i] * pa_x[i] - ta1_y_xxxxx_xxxz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxyy_0[i] = 5.0 * ta1_y_xxxx_xxyy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xxxxx_xyy_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xyy_1[i] * fe_0 + ta1_y_xxxxx_xxyy_0[i] * pa_x[i] - ta1_y_xxxxx_xxyy_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxyz_0[i] = 5.0 * ta1_y_xxxx_xxyz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xyz_1[i] * fe_0 + ta1_y_xxxxx_xxyz_0[i] * pa_x[i] - ta1_y_xxxxx_xxyz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxzz_0[i] = 5.0 * ta1_y_xxxx_xxzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxxxx_xzz_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xzz_1[i] * fe_0 + ta1_y_xxxxx_xxzz_0[i] * pa_x[i] - ta1_y_xxxxx_xxzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xyyy_0[i] = 5.0 * ta1_y_xxxx_xyyy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xyyy_1[i] * fe_0 + ta1_y_xxxxx_yyy_0[i] * fe_0 - ta1_y_xxxxx_yyy_1[i] * fe_0 + ta1_y_xxxxx_xyyy_0[i] * pa_x[i] - ta1_y_xxxxx_xyyy_1[i] * pc_x[i];

        ta1_y_xxxxxx_xyyz_0[i] = 5.0 * ta1_y_xxxx_xyyz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xyyz_1[i] * fe_0 + ta1_y_xxxxx_yyz_0[i] * fe_0 - ta1_y_xxxxx_yyz_1[i] * fe_0 + ta1_y_xxxxx_xyyz_0[i] * pa_x[i] - ta1_y_xxxxx_xyyz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xyzz_0[i] = 5.0 * ta1_y_xxxx_xyzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xyzz_1[i] * fe_0 + ta1_y_xxxxx_yzz_0[i] * fe_0 - ta1_y_xxxxx_yzz_1[i] * fe_0 + ta1_y_xxxxx_xyzz_0[i] * pa_x[i] - ta1_y_xxxxx_xyzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xzzz_0[i] = 5.0 * ta1_y_xxxx_xzzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xzzz_1[i] * fe_0 + ta1_y_xxxxx_zzz_0[i] * fe_0 - ta1_y_xxxxx_zzz_1[i] * fe_0 + ta1_y_xxxxx_xzzz_0[i] * pa_x[i] - ta1_y_xxxxx_xzzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_yyyy_0[i] = 5.0 * ta1_y_xxxx_yyyy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yyyy_1[i] * fe_0 + ta1_y_xxxxx_yyyy_0[i] * pa_x[i] - ta1_y_xxxxx_yyyy_1[i] * pc_x[i];

        ta1_y_xxxxxx_yyyz_0[i] = 5.0 * ta1_y_xxxx_yyyz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yyyz_1[i] * fe_0 + ta1_y_xxxxx_yyyz_0[i] * pa_x[i] - ta1_y_xxxxx_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxxx_yyzz_0[i] = 5.0 * ta1_y_xxxx_yyzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yyzz_1[i] * fe_0 + ta1_y_xxxxx_yyzz_0[i] * pa_x[i] - ta1_y_xxxxx_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_yzzz_0[i] = 5.0 * ta1_y_xxxx_yzzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yzzz_1[i] * fe_0 + ta1_y_xxxxx_yzzz_0[i] * pa_x[i] - ta1_y_xxxxx_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_zzzz_0[i] = 5.0 * ta1_y_xxxx_zzzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_zzzz_1[i] * fe_0 + ta1_y_xxxxx_zzzz_0[i] * pa_x[i] - ta1_y_xxxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 435-450 components of targeted buffer : IG

    auto ta1_y_xxxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 435);

    auto ta1_y_xxxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 436);

    auto ta1_y_xxxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 437);

    auto ta1_y_xxxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 438);

    auto ta1_y_xxxxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 439);

    auto ta1_y_xxxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 440);

    auto ta1_y_xxxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 441);

    auto ta1_y_xxxxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 442);

    auto ta1_y_xxxxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 443);

    auto ta1_y_xxxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 444);

    auto ta1_y_xxxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 445);

    auto ta1_y_xxxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 446);

    auto ta1_y_xxxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 447);

    auto ta1_y_xxxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 448);

    auto ta1_y_xxxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 449);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxxx_xxx_0, ta1_y_xxxxx_xxx_1, ta1_y_xxxxx_xxxx_0, ta1_y_xxxxx_xxxx_1, ta1_y_xxxxx_xxxy_0, ta1_y_xxxxx_xxxy_1, ta1_y_xxxxx_xxxz_0, ta1_y_xxxxx_xxxz_1, ta1_y_xxxxx_xxy_0, ta1_y_xxxxx_xxy_1, ta1_y_xxxxx_xxyy_0, ta1_y_xxxxx_xxyy_1, ta1_y_xxxxx_xxyz_0, ta1_y_xxxxx_xxyz_1, ta1_y_xxxxx_xxz_0, ta1_y_xxxxx_xxz_1, ta1_y_xxxxx_xxzz_0, ta1_y_xxxxx_xxzz_1, ta1_y_xxxxx_xyy_0, ta1_y_xxxxx_xyy_1, ta1_y_xxxxx_xyyy_0, ta1_y_xxxxx_xyyy_1, ta1_y_xxxxx_xyyz_0, ta1_y_xxxxx_xyyz_1, ta1_y_xxxxx_xyz_0, ta1_y_xxxxx_xyz_1, ta1_y_xxxxx_xyzz_0, ta1_y_xxxxx_xyzz_1, ta1_y_xxxxx_xzz_0, ta1_y_xxxxx_xzz_1, ta1_y_xxxxx_xzzz_0, ta1_y_xxxxx_xzzz_1, ta1_y_xxxxx_zzzz_0, ta1_y_xxxxx_zzzz_1, ta1_y_xxxxxy_xxxx_0, ta1_y_xxxxxy_xxxy_0, ta1_y_xxxxxy_xxxz_0, ta1_y_xxxxxy_xxyy_0, ta1_y_xxxxxy_xxyz_0, ta1_y_xxxxxy_xxzz_0, ta1_y_xxxxxy_xyyy_0, ta1_y_xxxxxy_xyyz_0, ta1_y_xxxxxy_xyzz_0, ta1_y_xxxxxy_xzzz_0, ta1_y_xxxxxy_yyyy_0, ta1_y_xxxxxy_yyyz_0, ta1_y_xxxxxy_yyzz_0, ta1_y_xxxxxy_yzzz_0, ta1_y_xxxxxy_zzzz_0, ta1_y_xxxxy_yyyy_0, ta1_y_xxxxy_yyyy_1, ta1_y_xxxxy_yyyz_0, ta1_y_xxxxy_yyyz_1, ta1_y_xxxxy_yyzz_0, ta1_y_xxxxy_yyzz_1, ta1_y_xxxxy_yzzz_0, ta1_y_xxxxy_yzzz_1, ta1_y_xxxy_yyyy_0, ta1_y_xxxy_yyyy_1, ta1_y_xxxy_yyyz_0, ta1_y_xxxy_yyyz_1, ta1_y_xxxy_yyzz_0, ta1_y_xxxy_yyzz_1, ta1_y_xxxy_yzzz_0, ta1_y_xxxy_yzzz_1, ta_xxxxx_xxxx_1, ta_xxxxx_xxxy_1, ta_xxxxx_xxxz_1, ta_xxxxx_xxyy_1, ta_xxxxx_xxyz_1, ta_xxxxx_xxzz_1, ta_xxxxx_xyyy_1, ta_xxxxx_xyyz_1, ta_xxxxx_xyzz_1, ta_xxxxx_xzzz_1, ta_xxxxx_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxy_xxxx_0[i] = ta_xxxxx_xxxx_1[i] + ta1_y_xxxxx_xxxx_0[i] * pa_y[i] - ta1_y_xxxxx_xxxx_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxxy_0[i] = ta1_y_xxxxx_xxx_0[i] * fe_0 - ta1_y_xxxxx_xxx_1[i] * fe_0 + ta_xxxxx_xxxy_1[i] + ta1_y_xxxxx_xxxy_0[i] * pa_y[i] - ta1_y_xxxxx_xxxy_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxxz_0[i] = ta_xxxxx_xxxz_1[i] + ta1_y_xxxxx_xxxz_0[i] * pa_y[i] - ta1_y_xxxxx_xxxz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxyy_0[i] = 2.0 * ta1_y_xxxxx_xxy_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xxy_1[i] * fe_0 + ta_xxxxx_xxyy_1[i] + ta1_y_xxxxx_xxyy_0[i] * pa_y[i] - ta1_y_xxxxx_xxyy_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxyz_0[i] = ta1_y_xxxxx_xxz_0[i] * fe_0 - ta1_y_xxxxx_xxz_1[i] * fe_0 + ta_xxxxx_xxyz_1[i] + ta1_y_xxxxx_xxyz_0[i] * pa_y[i] - ta1_y_xxxxx_xxyz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxzz_0[i] = ta_xxxxx_xxzz_1[i] + ta1_y_xxxxx_xxzz_0[i] * pa_y[i] - ta1_y_xxxxx_xxzz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xyyy_0[i] = 3.0 * ta1_y_xxxxx_xyy_0[i] * fe_0 - 3.0 * ta1_y_xxxxx_xyy_1[i] * fe_0 + ta_xxxxx_xyyy_1[i] + ta1_y_xxxxx_xyyy_0[i] * pa_y[i] - ta1_y_xxxxx_xyyy_1[i] * pc_y[i];

        ta1_y_xxxxxy_xyyz_0[i] = 2.0 * ta1_y_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xyz_1[i] * fe_0 + ta_xxxxx_xyyz_1[i] + ta1_y_xxxxx_xyyz_0[i] * pa_y[i] - ta1_y_xxxxx_xyyz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xyzz_0[i] = ta1_y_xxxxx_xzz_0[i] * fe_0 - ta1_y_xxxxx_xzz_1[i] * fe_0 + ta_xxxxx_xyzz_1[i] + ta1_y_xxxxx_xyzz_0[i] * pa_y[i] - ta1_y_xxxxx_xyzz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xzzz_0[i] = ta_xxxxx_xzzz_1[i] + ta1_y_xxxxx_xzzz_0[i] * pa_y[i] - ta1_y_xxxxx_xzzz_1[i] * pc_y[i];

        ta1_y_xxxxxy_yyyy_0[i] = 4.0 * ta1_y_xxxy_yyyy_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yyyy_1[i] * fe_0 + ta1_y_xxxxy_yyyy_0[i] * pa_x[i] - ta1_y_xxxxy_yyyy_1[i] * pc_x[i];

        ta1_y_xxxxxy_yyyz_0[i] = 4.0 * ta1_y_xxxy_yyyz_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yyyz_1[i] * fe_0 + ta1_y_xxxxy_yyyz_0[i] * pa_x[i] - ta1_y_xxxxy_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxxy_yyzz_0[i] = 4.0 * ta1_y_xxxy_yyzz_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yyzz_1[i] * fe_0 + ta1_y_xxxxy_yyzz_0[i] * pa_x[i] - ta1_y_xxxxy_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxxy_yzzz_0[i] = 4.0 * ta1_y_xxxy_yzzz_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yzzz_1[i] * fe_0 + ta1_y_xxxxy_yzzz_0[i] * pa_x[i] - ta1_y_xxxxy_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxxy_zzzz_0[i] = ta_xxxxx_zzzz_1[i] + ta1_y_xxxxx_zzzz_0[i] * pa_y[i] - ta1_y_xxxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 450-465 components of targeted buffer : IG

    auto ta1_y_xxxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 450);

    auto ta1_y_xxxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 451);

    auto ta1_y_xxxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 452);

    auto ta1_y_xxxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 453);

    auto ta1_y_xxxxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 454);

    auto ta1_y_xxxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 455);

    auto ta1_y_xxxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 456);

    auto ta1_y_xxxxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 457);

    auto ta1_y_xxxxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 458);

    auto ta1_y_xxxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 459);

    auto ta1_y_xxxxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 460);

    auto ta1_y_xxxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 461);

    auto ta1_y_xxxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 462);

    auto ta1_y_xxxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 463);

    auto ta1_y_xxxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 464);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxxx_xxx_0, ta1_y_xxxxx_xxx_1, ta1_y_xxxxx_xxxx_0, ta1_y_xxxxx_xxxx_1, ta1_y_xxxxx_xxxy_0, ta1_y_xxxxx_xxxy_1, ta1_y_xxxxx_xxxz_0, ta1_y_xxxxx_xxxz_1, ta1_y_xxxxx_xxy_0, ta1_y_xxxxx_xxy_1, ta1_y_xxxxx_xxyy_0, ta1_y_xxxxx_xxyy_1, ta1_y_xxxxx_xxyz_0, ta1_y_xxxxx_xxyz_1, ta1_y_xxxxx_xxz_0, ta1_y_xxxxx_xxz_1, ta1_y_xxxxx_xxzz_0, ta1_y_xxxxx_xxzz_1, ta1_y_xxxxx_xyy_0, ta1_y_xxxxx_xyy_1, ta1_y_xxxxx_xyyy_0, ta1_y_xxxxx_xyyy_1, ta1_y_xxxxx_xyyz_0, ta1_y_xxxxx_xyyz_1, ta1_y_xxxxx_xyz_0, ta1_y_xxxxx_xyz_1, ta1_y_xxxxx_xyzz_0, ta1_y_xxxxx_xyzz_1, ta1_y_xxxxx_xzz_0, ta1_y_xxxxx_xzz_1, ta1_y_xxxxx_xzzz_0, ta1_y_xxxxx_xzzz_1, ta1_y_xxxxx_yyyy_0, ta1_y_xxxxx_yyyy_1, ta1_y_xxxxxz_xxxx_0, ta1_y_xxxxxz_xxxy_0, ta1_y_xxxxxz_xxxz_0, ta1_y_xxxxxz_xxyy_0, ta1_y_xxxxxz_xxyz_0, ta1_y_xxxxxz_xxzz_0, ta1_y_xxxxxz_xyyy_0, ta1_y_xxxxxz_xyyz_0, ta1_y_xxxxxz_xyzz_0, ta1_y_xxxxxz_xzzz_0, ta1_y_xxxxxz_yyyy_0, ta1_y_xxxxxz_yyyz_0, ta1_y_xxxxxz_yyzz_0, ta1_y_xxxxxz_yzzz_0, ta1_y_xxxxxz_zzzz_0, ta1_y_xxxxz_yyyz_0, ta1_y_xxxxz_yyyz_1, ta1_y_xxxxz_yyzz_0, ta1_y_xxxxz_yyzz_1, ta1_y_xxxxz_yzzz_0, ta1_y_xxxxz_yzzz_1, ta1_y_xxxxz_zzzz_0, ta1_y_xxxxz_zzzz_1, ta1_y_xxxz_yyyz_0, ta1_y_xxxz_yyyz_1, ta1_y_xxxz_yyzz_0, ta1_y_xxxz_yyzz_1, ta1_y_xxxz_yzzz_0, ta1_y_xxxz_yzzz_1, ta1_y_xxxz_zzzz_0, ta1_y_xxxz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxz_xxxx_0[i] = ta1_y_xxxxx_xxxx_0[i] * pa_z[i] - ta1_y_xxxxx_xxxx_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxxy_0[i] = ta1_y_xxxxx_xxxy_0[i] * pa_z[i] - ta1_y_xxxxx_xxxy_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxxz_0[i] = ta1_y_xxxxx_xxx_0[i] * fe_0 - ta1_y_xxxxx_xxx_1[i] * fe_0 + ta1_y_xxxxx_xxxz_0[i] * pa_z[i] - ta1_y_xxxxx_xxxz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxyy_0[i] = ta1_y_xxxxx_xxyy_0[i] * pa_z[i] - ta1_y_xxxxx_xxyy_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxyz_0[i] = ta1_y_xxxxx_xxy_0[i] * fe_0 - ta1_y_xxxxx_xxy_1[i] * fe_0 + ta1_y_xxxxx_xxyz_0[i] * pa_z[i] - ta1_y_xxxxx_xxyz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxzz_0[i] = 2.0 * ta1_y_xxxxx_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xxz_1[i] * fe_0 + ta1_y_xxxxx_xxzz_0[i] * pa_z[i] - ta1_y_xxxxx_xxzz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xyyy_0[i] = ta1_y_xxxxx_xyyy_0[i] * pa_z[i] - ta1_y_xxxxx_xyyy_1[i] * pc_z[i];

        ta1_y_xxxxxz_xyyz_0[i] = ta1_y_xxxxx_xyy_0[i] * fe_0 - ta1_y_xxxxx_xyy_1[i] * fe_0 + ta1_y_xxxxx_xyyz_0[i] * pa_z[i] - ta1_y_xxxxx_xyyz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xyzz_0[i] = 2.0 * ta1_y_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xyz_1[i] * fe_0 + ta1_y_xxxxx_xyzz_0[i] * pa_z[i] - ta1_y_xxxxx_xyzz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xzzz_0[i] = 3.0 * ta1_y_xxxxx_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxxxx_xzz_1[i] * fe_0 + ta1_y_xxxxx_xzzz_0[i] * pa_z[i] - ta1_y_xxxxx_xzzz_1[i] * pc_z[i];

        ta1_y_xxxxxz_yyyy_0[i] = ta1_y_xxxxx_yyyy_0[i] * pa_z[i] - ta1_y_xxxxx_yyyy_1[i] * pc_z[i];

        ta1_y_xxxxxz_yyyz_0[i] = 4.0 * ta1_y_xxxz_yyyz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_yyyz_1[i] * fe_0 + ta1_y_xxxxz_yyyz_0[i] * pa_x[i] - ta1_y_xxxxz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxxz_yyzz_0[i] = 4.0 * ta1_y_xxxz_yyzz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_yyzz_1[i] * fe_0 + ta1_y_xxxxz_yyzz_0[i] * pa_x[i] - ta1_y_xxxxz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxxz_yzzz_0[i] = 4.0 * ta1_y_xxxz_yzzz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_yzzz_1[i] * fe_0 + ta1_y_xxxxz_yzzz_0[i] * pa_x[i] - ta1_y_xxxxz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxxz_zzzz_0[i] = 4.0 * ta1_y_xxxz_zzzz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_zzzz_1[i] * fe_0 + ta1_y_xxxxz_zzzz_0[i] * pa_x[i] - ta1_y_xxxxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 465-480 components of targeted buffer : IG

    auto ta1_y_xxxxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 465);

    auto ta1_y_xxxxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 466);

    auto ta1_y_xxxxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 467);

    auto ta1_y_xxxxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 468);

    auto ta1_y_xxxxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 469);

    auto ta1_y_xxxxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 470);

    auto ta1_y_xxxxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 471);

    auto ta1_y_xxxxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 472);

    auto ta1_y_xxxxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 473);

    auto ta1_y_xxxxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 474);

    auto ta1_y_xxxxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 475);

    auto ta1_y_xxxxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 476);

    auto ta1_y_xxxxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 477);

    auto ta1_y_xxxxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 478);

    auto ta1_y_xxxxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 479);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxx_xxxx_0, ta1_y_xxxx_xxxx_1, ta1_y_xxxx_xxxz_0, ta1_y_xxxx_xxxz_1, ta1_y_xxxx_xxzz_0, ta1_y_xxxx_xxzz_1, ta1_y_xxxx_xzzz_0, ta1_y_xxxx_xzzz_1, ta1_y_xxxxy_xxxx_0, ta1_y_xxxxy_xxxx_1, ta1_y_xxxxy_xxxz_0, ta1_y_xxxxy_xxxz_1, ta1_y_xxxxy_xxzz_0, ta1_y_xxxxy_xxzz_1, ta1_y_xxxxy_xzzz_0, ta1_y_xxxxy_xzzz_1, ta1_y_xxxxyy_xxxx_0, ta1_y_xxxxyy_xxxy_0, ta1_y_xxxxyy_xxxz_0, ta1_y_xxxxyy_xxyy_0, ta1_y_xxxxyy_xxyz_0, ta1_y_xxxxyy_xxzz_0, ta1_y_xxxxyy_xyyy_0, ta1_y_xxxxyy_xyyz_0, ta1_y_xxxxyy_xyzz_0, ta1_y_xxxxyy_xzzz_0, ta1_y_xxxxyy_yyyy_0, ta1_y_xxxxyy_yyyz_0, ta1_y_xxxxyy_yyzz_0, ta1_y_xxxxyy_yzzz_0, ta1_y_xxxxyy_zzzz_0, ta1_y_xxxyy_xxxy_0, ta1_y_xxxyy_xxxy_1, ta1_y_xxxyy_xxy_0, ta1_y_xxxyy_xxy_1, ta1_y_xxxyy_xxyy_0, ta1_y_xxxyy_xxyy_1, ta1_y_xxxyy_xxyz_0, ta1_y_xxxyy_xxyz_1, ta1_y_xxxyy_xyy_0, ta1_y_xxxyy_xyy_1, ta1_y_xxxyy_xyyy_0, ta1_y_xxxyy_xyyy_1, ta1_y_xxxyy_xyyz_0, ta1_y_xxxyy_xyyz_1, ta1_y_xxxyy_xyz_0, ta1_y_xxxyy_xyz_1, ta1_y_xxxyy_xyzz_0, ta1_y_xxxyy_xyzz_1, ta1_y_xxxyy_yyy_0, ta1_y_xxxyy_yyy_1, ta1_y_xxxyy_yyyy_0, ta1_y_xxxyy_yyyy_1, ta1_y_xxxyy_yyyz_0, ta1_y_xxxyy_yyyz_1, ta1_y_xxxyy_yyz_0, ta1_y_xxxyy_yyz_1, ta1_y_xxxyy_yyzz_0, ta1_y_xxxyy_yyzz_1, ta1_y_xxxyy_yzz_0, ta1_y_xxxyy_yzz_1, ta1_y_xxxyy_yzzz_0, ta1_y_xxxyy_yzzz_1, ta1_y_xxxyy_zzzz_0, ta1_y_xxxyy_zzzz_1, ta1_y_xxyy_xxxy_0, ta1_y_xxyy_xxxy_1, ta1_y_xxyy_xxyy_0, ta1_y_xxyy_xxyy_1, ta1_y_xxyy_xxyz_0, ta1_y_xxyy_xxyz_1, ta1_y_xxyy_xyyy_0, ta1_y_xxyy_xyyy_1, ta1_y_xxyy_xyyz_0, ta1_y_xxyy_xyyz_1, ta1_y_xxyy_xyzz_0, ta1_y_xxyy_xyzz_1, ta1_y_xxyy_yyyy_0, ta1_y_xxyy_yyyy_1, ta1_y_xxyy_yyyz_0, ta1_y_xxyy_yyyz_1, ta1_y_xxyy_yyzz_0, ta1_y_xxyy_yyzz_1, ta1_y_xxyy_yzzz_0, ta1_y_xxyy_yzzz_1, ta1_y_xxyy_zzzz_0, ta1_y_xxyy_zzzz_1, ta_xxxxy_xxxx_1, ta_xxxxy_xxxz_1, ta_xxxxy_xxzz_1, ta_xxxxy_xzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyy_xxxx_0[i] = ta1_y_xxxx_xxxx_0[i] * fe_0 - ta1_y_xxxx_xxxx_1[i] * fe_0 + ta_xxxxy_xxxx_1[i] + ta1_y_xxxxy_xxxx_0[i] * pa_y[i] - ta1_y_xxxxy_xxxx_1[i] * pc_y[i];

        ta1_y_xxxxyy_xxxy_0[i] = 3.0 * ta1_y_xxyy_xxxy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xxxyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_xxxyy_xxy_1[i] * fe_0 + ta1_y_xxxyy_xxxy_0[i] * pa_x[i] - ta1_y_xxxyy_xxxy_1[i] * pc_x[i];

        ta1_y_xxxxyy_xxxz_0[i] = ta1_y_xxxx_xxxz_0[i] * fe_0 - ta1_y_xxxx_xxxz_1[i] * fe_0 + ta_xxxxy_xxxz_1[i] + ta1_y_xxxxy_xxxz_0[i] * pa_y[i] - ta1_y_xxxxy_xxxz_1[i] * pc_y[i];

        ta1_y_xxxxyy_xxyy_0[i] = 3.0 * ta1_y_xxyy_xxyy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xxxyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_xxxyy_xyy_1[i] * fe_0 + ta1_y_xxxyy_xxyy_0[i] * pa_x[i] - ta1_y_xxxyy_xxyy_1[i] * pc_x[i];

        ta1_y_xxxxyy_xxyz_0[i] = 3.0 * ta1_y_xxyy_xxyz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxyy_xyz_1[i] * fe_0 + ta1_y_xxxyy_xxyz_0[i] * pa_x[i] - ta1_y_xxxyy_xxyz_1[i] * pc_x[i];

        ta1_y_xxxxyy_xxzz_0[i] = ta1_y_xxxx_xxzz_0[i] * fe_0 - ta1_y_xxxx_xxzz_1[i] * fe_0 + ta_xxxxy_xxzz_1[i] + ta1_y_xxxxy_xxzz_0[i] * pa_y[i] - ta1_y_xxxxy_xxzz_1[i] * pc_y[i];

        ta1_y_xxxxyy_xyyy_0[i] = 3.0 * ta1_y_xxyy_xyyy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xyyy_1[i] * fe_0 + ta1_y_xxxyy_yyy_0[i] * fe_0 - ta1_y_xxxyy_yyy_1[i] * fe_0 + ta1_y_xxxyy_xyyy_0[i] * pa_x[i] - ta1_y_xxxyy_xyyy_1[i] * pc_x[i];

        ta1_y_xxxxyy_xyyz_0[i] = 3.0 * ta1_y_xxyy_xyyz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xyyz_1[i] * fe_0 + ta1_y_xxxyy_yyz_0[i] * fe_0 - ta1_y_xxxyy_yyz_1[i] * fe_0 + ta1_y_xxxyy_xyyz_0[i] * pa_x[i] - ta1_y_xxxyy_xyyz_1[i] * pc_x[i];

        ta1_y_xxxxyy_xyzz_0[i] = 3.0 * ta1_y_xxyy_xyzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xyzz_1[i] * fe_0 + ta1_y_xxxyy_yzz_0[i] * fe_0 - ta1_y_xxxyy_yzz_1[i] * fe_0 + ta1_y_xxxyy_xyzz_0[i] * pa_x[i] - ta1_y_xxxyy_xyzz_1[i] * pc_x[i];

        ta1_y_xxxxyy_xzzz_0[i] = ta1_y_xxxx_xzzz_0[i] * fe_0 - ta1_y_xxxx_xzzz_1[i] * fe_0 + ta_xxxxy_xzzz_1[i] + ta1_y_xxxxy_xzzz_0[i] * pa_y[i] - ta1_y_xxxxy_xzzz_1[i] * pc_y[i];

        ta1_y_xxxxyy_yyyy_0[i] = 3.0 * ta1_y_xxyy_yyyy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yyyy_1[i] * fe_0 + ta1_y_xxxyy_yyyy_0[i] * pa_x[i] - ta1_y_xxxyy_yyyy_1[i] * pc_x[i];

        ta1_y_xxxxyy_yyyz_0[i] = 3.0 * ta1_y_xxyy_yyyz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yyyz_1[i] * fe_0 + ta1_y_xxxyy_yyyz_0[i] * pa_x[i] - ta1_y_xxxyy_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxyy_yyzz_0[i] = 3.0 * ta1_y_xxyy_yyzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yyzz_1[i] * fe_0 + ta1_y_xxxyy_yyzz_0[i] * pa_x[i] - ta1_y_xxxyy_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxyy_yzzz_0[i] = 3.0 * ta1_y_xxyy_yzzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yzzz_1[i] * fe_0 + ta1_y_xxxyy_yzzz_0[i] * pa_x[i] - ta1_y_xxxyy_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxyy_zzzz_0[i] = 3.0 * ta1_y_xxyy_zzzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_zzzz_1[i] * fe_0 + ta1_y_xxxyy_zzzz_0[i] * pa_x[i] - ta1_y_xxxyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 480-495 components of targeted buffer : IG

    auto ta1_y_xxxxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 480);

    auto ta1_y_xxxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 481);

    auto ta1_y_xxxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 482);

    auto ta1_y_xxxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 483);

    auto ta1_y_xxxxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 484);

    auto ta1_y_xxxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 485);

    auto ta1_y_xxxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 486);

    auto ta1_y_xxxxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 487);

    auto ta1_y_xxxxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 488);

    auto ta1_y_xxxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 489);

    auto ta1_y_xxxxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 490);

    auto ta1_y_xxxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 491);

    auto ta1_y_xxxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 492);

    auto ta1_y_xxxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 493);

    auto ta1_y_xxxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 494);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxxxy_xxxx_0, ta1_y_xxxxy_xxxx_1, ta1_y_xxxxy_xxxy_0, ta1_y_xxxxy_xxxy_1, ta1_y_xxxxy_xxy_0, ta1_y_xxxxy_xxy_1, ta1_y_xxxxy_xxyy_0, ta1_y_xxxxy_xxyy_1, ta1_y_xxxxy_xxyz_0, ta1_y_xxxxy_xxyz_1, ta1_y_xxxxy_xyy_0, ta1_y_xxxxy_xyy_1, ta1_y_xxxxy_xyyy_0, ta1_y_xxxxy_xyyy_1, ta1_y_xxxxy_xyyz_0, ta1_y_xxxxy_xyyz_1, ta1_y_xxxxy_xyz_0, ta1_y_xxxxy_xyz_1, ta1_y_xxxxy_xyzz_0, ta1_y_xxxxy_xyzz_1, ta1_y_xxxxy_yyyy_0, ta1_y_xxxxy_yyyy_1, ta1_y_xxxxyz_xxxx_0, ta1_y_xxxxyz_xxxy_0, ta1_y_xxxxyz_xxxz_0, ta1_y_xxxxyz_xxyy_0, ta1_y_xxxxyz_xxyz_0, ta1_y_xxxxyz_xxzz_0, ta1_y_xxxxyz_xyyy_0, ta1_y_xxxxyz_xyyz_0, ta1_y_xxxxyz_xyzz_0, ta1_y_xxxxyz_xzzz_0, ta1_y_xxxxyz_yyyy_0, ta1_y_xxxxyz_yyyz_0, ta1_y_xxxxyz_yyzz_0, ta1_y_xxxxyz_yzzz_0, ta1_y_xxxxyz_zzzz_0, ta1_y_xxxxz_xxxz_0, ta1_y_xxxxz_xxxz_1, ta1_y_xxxxz_xxzz_0, ta1_y_xxxxz_xxzz_1, ta1_y_xxxxz_xzzz_0, ta1_y_xxxxz_xzzz_1, ta1_y_xxxxz_zzzz_0, ta1_y_xxxxz_zzzz_1, ta1_y_xxxyz_yyyz_0, ta1_y_xxxyz_yyyz_1, ta1_y_xxxyz_yyzz_0, ta1_y_xxxyz_yyzz_1, ta1_y_xxxyz_yzzz_0, ta1_y_xxxyz_yzzz_1, ta1_y_xxyz_yyyz_0, ta1_y_xxyz_yyyz_1, ta1_y_xxyz_yyzz_0, ta1_y_xxyz_yyzz_1, ta1_y_xxyz_yzzz_0, ta1_y_xxyz_yzzz_1, ta_xxxxz_xxxz_1, ta_xxxxz_xxzz_1, ta_xxxxz_xzzz_1, ta_xxxxz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyz_xxxx_0[i] = ta1_y_xxxxy_xxxx_0[i] * pa_z[i] - ta1_y_xxxxy_xxxx_1[i] * pc_z[i];

        ta1_y_xxxxyz_xxxy_0[i] = ta1_y_xxxxy_xxxy_0[i] * pa_z[i] - ta1_y_xxxxy_xxxy_1[i] * pc_z[i];

        ta1_y_xxxxyz_xxxz_0[i] = ta_xxxxz_xxxz_1[i] + ta1_y_xxxxz_xxxz_0[i] * pa_y[i] - ta1_y_xxxxz_xxxz_1[i] * pc_y[i];

        ta1_y_xxxxyz_xxyy_0[i] = ta1_y_xxxxy_xxyy_0[i] * pa_z[i] - ta1_y_xxxxy_xxyy_1[i] * pc_z[i];

        ta1_y_xxxxyz_xxyz_0[i] = ta1_y_xxxxy_xxy_0[i] * fe_0 - ta1_y_xxxxy_xxy_1[i] * fe_0 + ta1_y_xxxxy_xxyz_0[i] * pa_z[i] - ta1_y_xxxxy_xxyz_1[i] * pc_z[i];

        ta1_y_xxxxyz_xxzz_0[i] = ta_xxxxz_xxzz_1[i] + ta1_y_xxxxz_xxzz_0[i] * pa_y[i] - ta1_y_xxxxz_xxzz_1[i] * pc_y[i];

        ta1_y_xxxxyz_xyyy_0[i] = ta1_y_xxxxy_xyyy_0[i] * pa_z[i] - ta1_y_xxxxy_xyyy_1[i] * pc_z[i];

        ta1_y_xxxxyz_xyyz_0[i] = ta1_y_xxxxy_xyy_0[i] * fe_0 - ta1_y_xxxxy_xyy_1[i] * fe_0 + ta1_y_xxxxy_xyyz_0[i] * pa_z[i] - ta1_y_xxxxy_xyyz_1[i] * pc_z[i];

        ta1_y_xxxxyz_xyzz_0[i] = 2.0 * ta1_y_xxxxy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxxy_xyz_1[i] * fe_0 + ta1_y_xxxxy_xyzz_0[i] * pa_z[i] - ta1_y_xxxxy_xyzz_1[i] * pc_z[i];

        ta1_y_xxxxyz_xzzz_0[i] = ta_xxxxz_xzzz_1[i] + ta1_y_xxxxz_xzzz_0[i] * pa_y[i] - ta1_y_xxxxz_xzzz_1[i] * pc_y[i];

        ta1_y_xxxxyz_yyyy_0[i] = ta1_y_xxxxy_yyyy_0[i] * pa_z[i] - ta1_y_xxxxy_yyyy_1[i] * pc_z[i];

        ta1_y_xxxxyz_yyyz_0[i] = 3.0 * ta1_y_xxyz_yyyz_0[i] * fe_0 - 3.0 * ta1_y_xxyz_yyyz_1[i] * fe_0 + ta1_y_xxxyz_yyyz_0[i] * pa_x[i] - ta1_y_xxxyz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxyz_yyzz_0[i] = 3.0 * ta1_y_xxyz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_xxyz_yyzz_1[i] * fe_0 + ta1_y_xxxyz_yyzz_0[i] * pa_x[i] - ta1_y_xxxyz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxyz_yzzz_0[i] = 3.0 * ta1_y_xxyz_yzzz_0[i] * fe_0 - 3.0 * ta1_y_xxyz_yzzz_1[i] * fe_0 + ta1_y_xxxyz_yzzz_0[i] * pa_x[i] - ta1_y_xxxyz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxyz_zzzz_0[i] = ta_xxxxz_zzzz_1[i] + ta1_y_xxxxz_zzzz_0[i] * pa_y[i] - ta1_y_xxxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 495-510 components of targeted buffer : IG

    auto ta1_y_xxxxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 495);

    auto ta1_y_xxxxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 496);

    auto ta1_y_xxxxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 497);

    auto ta1_y_xxxxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 498);

    auto ta1_y_xxxxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 499);

    auto ta1_y_xxxxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 500);

    auto ta1_y_xxxxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 501);

    auto ta1_y_xxxxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 502);

    auto ta1_y_xxxxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 503);

    auto ta1_y_xxxxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 504);

    auto ta1_y_xxxxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 505);

    auto ta1_y_xxxxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 506);

    auto ta1_y_xxxxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 507);

    auto ta1_y_xxxxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 508);

    auto ta1_y_xxxxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 509);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxx_xxxx_0, ta1_y_xxxx_xxxx_1, ta1_y_xxxx_xxxy_0, ta1_y_xxxx_xxxy_1, ta1_y_xxxx_xxyy_0, ta1_y_xxxx_xxyy_1, ta1_y_xxxx_xyyy_0, ta1_y_xxxx_xyyy_1, ta1_y_xxxxz_xxxx_0, ta1_y_xxxxz_xxxx_1, ta1_y_xxxxz_xxxy_0, ta1_y_xxxxz_xxxy_1, ta1_y_xxxxz_xxyy_0, ta1_y_xxxxz_xxyy_1, ta1_y_xxxxz_xyyy_0, ta1_y_xxxxz_xyyy_1, ta1_y_xxxxzz_xxxx_0, ta1_y_xxxxzz_xxxy_0, ta1_y_xxxxzz_xxxz_0, ta1_y_xxxxzz_xxyy_0, ta1_y_xxxxzz_xxyz_0, ta1_y_xxxxzz_xxzz_0, ta1_y_xxxxzz_xyyy_0, ta1_y_xxxxzz_xyyz_0, ta1_y_xxxxzz_xyzz_0, ta1_y_xxxxzz_xzzz_0, ta1_y_xxxxzz_yyyy_0, ta1_y_xxxxzz_yyyz_0, ta1_y_xxxxzz_yyzz_0, ta1_y_xxxxzz_yzzz_0, ta1_y_xxxxzz_zzzz_0, ta1_y_xxxzz_xxxz_0, ta1_y_xxxzz_xxxz_1, ta1_y_xxxzz_xxyz_0, ta1_y_xxxzz_xxyz_1, ta1_y_xxxzz_xxz_0, ta1_y_xxxzz_xxz_1, ta1_y_xxxzz_xxzz_0, ta1_y_xxxzz_xxzz_1, ta1_y_xxxzz_xyyz_0, ta1_y_xxxzz_xyyz_1, ta1_y_xxxzz_xyz_0, ta1_y_xxxzz_xyz_1, ta1_y_xxxzz_xyzz_0, ta1_y_xxxzz_xyzz_1, ta1_y_xxxzz_xzz_0, ta1_y_xxxzz_xzz_1, ta1_y_xxxzz_xzzz_0, ta1_y_xxxzz_xzzz_1, ta1_y_xxxzz_yyyy_0, ta1_y_xxxzz_yyyy_1, ta1_y_xxxzz_yyyz_0, ta1_y_xxxzz_yyyz_1, ta1_y_xxxzz_yyz_0, ta1_y_xxxzz_yyz_1, ta1_y_xxxzz_yyzz_0, ta1_y_xxxzz_yyzz_1, ta1_y_xxxzz_yzz_0, ta1_y_xxxzz_yzz_1, ta1_y_xxxzz_yzzz_0, ta1_y_xxxzz_yzzz_1, ta1_y_xxxzz_zzz_0, ta1_y_xxxzz_zzz_1, ta1_y_xxxzz_zzzz_0, ta1_y_xxxzz_zzzz_1, ta1_y_xxzz_xxxz_0, ta1_y_xxzz_xxxz_1, ta1_y_xxzz_xxyz_0, ta1_y_xxzz_xxyz_1, ta1_y_xxzz_xxzz_0, ta1_y_xxzz_xxzz_1, ta1_y_xxzz_xyyz_0, ta1_y_xxzz_xyyz_1, ta1_y_xxzz_xyzz_0, ta1_y_xxzz_xyzz_1, ta1_y_xxzz_xzzz_0, ta1_y_xxzz_xzzz_1, ta1_y_xxzz_yyyy_0, ta1_y_xxzz_yyyy_1, ta1_y_xxzz_yyyz_0, ta1_y_xxzz_yyyz_1, ta1_y_xxzz_yyzz_0, ta1_y_xxzz_yyzz_1, ta1_y_xxzz_yzzz_0, ta1_y_xxzz_yzzz_1, ta1_y_xxzz_zzzz_0, ta1_y_xxzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxzz_xxxx_0[i] = ta1_y_xxxx_xxxx_0[i] * fe_0 - ta1_y_xxxx_xxxx_1[i] * fe_0 + ta1_y_xxxxz_xxxx_0[i] * pa_z[i] - ta1_y_xxxxz_xxxx_1[i] * pc_z[i];

        ta1_y_xxxxzz_xxxy_0[i] = ta1_y_xxxx_xxxy_0[i] * fe_0 - ta1_y_xxxx_xxxy_1[i] * fe_0 + ta1_y_xxxxz_xxxy_0[i] * pa_z[i] - ta1_y_xxxxz_xxxy_1[i] * pc_z[i];

        ta1_y_xxxxzz_xxxz_0[i] = 3.0 * ta1_y_xxzz_xxxz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xxxzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_xxxzz_xxz_1[i] * fe_0 + ta1_y_xxxzz_xxxz_0[i] * pa_x[i] - ta1_y_xxxzz_xxxz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xxyy_0[i] = ta1_y_xxxx_xxyy_0[i] * fe_0 - ta1_y_xxxx_xxyy_1[i] * fe_0 + ta1_y_xxxxz_xxyy_0[i] * pa_z[i] - ta1_y_xxxxz_xxyy_1[i] * pc_z[i];

        ta1_y_xxxxzz_xxyz_0[i] = 3.0 * ta1_y_xxzz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxzz_xyz_1[i] * fe_0 + ta1_y_xxxzz_xxyz_0[i] * pa_x[i] - ta1_y_xxxzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xxzz_0[i] = 3.0 * ta1_y_xxzz_xxzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxxzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_xxxzz_xzz_1[i] * fe_0 + ta1_y_xxxzz_xxzz_0[i] * pa_x[i] - ta1_y_xxxzz_xxzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xyyy_0[i] = ta1_y_xxxx_xyyy_0[i] * fe_0 - ta1_y_xxxx_xyyy_1[i] * fe_0 + ta1_y_xxxxz_xyyy_0[i] * pa_z[i] - ta1_y_xxxxz_xyyy_1[i] * pc_z[i];

        ta1_y_xxxxzz_xyyz_0[i] = 3.0 * ta1_y_xxzz_xyyz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xyyz_1[i] * fe_0 + ta1_y_xxxzz_yyz_0[i] * fe_0 - ta1_y_xxxzz_yyz_1[i] * fe_0 + ta1_y_xxxzz_xyyz_0[i] * pa_x[i] - ta1_y_xxxzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xyzz_0[i] = 3.0 * ta1_y_xxzz_xyzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xyzz_1[i] * fe_0 + ta1_y_xxxzz_yzz_0[i] * fe_0 - ta1_y_xxxzz_yzz_1[i] * fe_0 + ta1_y_xxxzz_xyzz_0[i] * pa_x[i] - ta1_y_xxxzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xzzz_0[i] = 3.0 * ta1_y_xxzz_xzzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xzzz_1[i] * fe_0 + ta1_y_xxxzz_zzz_0[i] * fe_0 - ta1_y_xxxzz_zzz_1[i] * fe_0 + ta1_y_xxxzz_xzzz_0[i] * pa_x[i] - ta1_y_xxxzz_xzzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_yyyy_0[i] = 3.0 * ta1_y_xxzz_yyyy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yyyy_1[i] * fe_0 + ta1_y_xxxzz_yyyy_0[i] * pa_x[i] - ta1_y_xxxzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxxxzz_yyyz_0[i] = 3.0 * ta1_y_xxzz_yyyz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yyyz_1[i] * fe_0 + ta1_y_xxxzz_yyyz_0[i] * pa_x[i] - ta1_y_xxxzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxxzz_yyzz_0[i] = 3.0 * ta1_y_xxzz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yyzz_1[i] * fe_0 + ta1_y_xxxzz_yyzz_0[i] * pa_x[i] - ta1_y_xxxzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_yzzz_0[i] = 3.0 * ta1_y_xxzz_yzzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yzzz_1[i] * fe_0 + ta1_y_xxxzz_yzzz_0[i] * pa_x[i] - ta1_y_xxxzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_zzzz_0[i] = 3.0 * ta1_y_xxzz_zzzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_zzzz_1[i] * fe_0 + ta1_y_xxxzz_zzzz_0[i] * pa_x[i] - ta1_y_xxxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 510-525 components of targeted buffer : IG

    auto ta1_y_xxxyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 510);

    auto ta1_y_xxxyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 511);

    auto ta1_y_xxxyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 512);

    auto ta1_y_xxxyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 513);

    auto ta1_y_xxxyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 514);

    auto ta1_y_xxxyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 515);

    auto ta1_y_xxxyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 516);

    auto ta1_y_xxxyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 517);

    auto ta1_y_xxxyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 518);

    auto ta1_y_xxxyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 519);

    auto ta1_y_xxxyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 520);

    auto ta1_y_xxxyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 521);

    auto ta1_y_xxxyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 522);

    auto ta1_y_xxxyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 523);

    auto ta1_y_xxxyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 524);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxy_xxxx_0, ta1_y_xxxy_xxxx_1, ta1_y_xxxy_xxxz_0, ta1_y_xxxy_xxxz_1, ta1_y_xxxy_xxzz_0, ta1_y_xxxy_xxzz_1, ta1_y_xxxy_xzzz_0, ta1_y_xxxy_xzzz_1, ta1_y_xxxyy_xxxx_0, ta1_y_xxxyy_xxxx_1, ta1_y_xxxyy_xxxz_0, ta1_y_xxxyy_xxxz_1, ta1_y_xxxyy_xxzz_0, ta1_y_xxxyy_xxzz_1, ta1_y_xxxyy_xzzz_0, ta1_y_xxxyy_xzzz_1, ta1_y_xxxyyy_xxxx_0, ta1_y_xxxyyy_xxxy_0, ta1_y_xxxyyy_xxxz_0, ta1_y_xxxyyy_xxyy_0, ta1_y_xxxyyy_xxyz_0, ta1_y_xxxyyy_xxzz_0, ta1_y_xxxyyy_xyyy_0, ta1_y_xxxyyy_xyyz_0, ta1_y_xxxyyy_xyzz_0, ta1_y_xxxyyy_xzzz_0, ta1_y_xxxyyy_yyyy_0, ta1_y_xxxyyy_yyyz_0, ta1_y_xxxyyy_yyzz_0, ta1_y_xxxyyy_yzzz_0, ta1_y_xxxyyy_zzzz_0, ta1_y_xxyyy_xxxy_0, ta1_y_xxyyy_xxxy_1, ta1_y_xxyyy_xxy_0, ta1_y_xxyyy_xxy_1, ta1_y_xxyyy_xxyy_0, ta1_y_xxyyy_xxyy_1, ta1_y_xxyyy_xxyz_0, ta1_y_xxyyy_xxyz_1, ta1_y_xxyyy_xyy_0, ta1_y_xxyyy_xyy_1, ta1_y_xxyyy_xyyy_0, ta1_y_xxyyy_xyyy_1, ta1_y_xxyyy_xyyz_0, ta1_y_xxyyy_xyyz_1, ta1_y_xxyyy_xyz_0, ta1_y_xxyyy_xyz_1, ta1_y_xxyyy_xyzz_0, ta1_y_xxyyy_xyzz_1, ta1_y_xxyyy_yyy_0, ta1_y_xxyyy_yyy_1, ta1_y_xxyyy_yyyy_0, ta1_y_xxyyy_yyyy_1, ta1_y_xxyyy_yyyz_0, ta1_y_xxyyy_yyyz_1, ta1_y_xxyyy_yyz_0, ta1_y_xxyyy_yyz_1, ta1_y_xxyyy_yyzz_0, ta1_y_xxyyy_yyzz_1, ta1_y_xxyyy_yzz_0, ta1_y_xxyyy_yzz_1, ta1_y_xxyyy_yzzz_0, ta1_y_xxyyy_yzzz_1, ta1_y_xxyyy_zzzz_0, ta1_y_xxyyy_zzzz_1, ta1_y_xyyy_xxxy_0, ta1_y_xyyy_xxxy_1, ta1_y_xyyy_xxyy_0, ta1_y_xyyy_xxyy_1, ta1_y_xyyy_xxyz_0, ta1_y_xyyy_xxyz_1, ta1_y_xyyy_xyyy_0, ta1_y_xyyy_xyyy_1, ta1_y_xyyy_xyyz_0, ta1_y_xyyy_xyyz_1, ta1_y_xyyy_xyzz_0, ta1_y_xyyy_xyzz_1, ta1_y_xyyy_yyyy_0, ta1_y_xyyy_yyyy_1, ta1_y_xyyy_yyyz_0, ta1_y_xyyy_yyyz_1, ta1_y_xyyy_yyzz_0, ta1_y_xyyy_yyzz_1, ta1_y_xyyy_yzzz_0, ta1_y_xyyy_yzzz_1, ta1_y_xyyy_zzzz_0, ta1_y_xyyy_zzzz_1, ta_xxxyy_xxxx_1, ta_xxxyy_xxxz_1, ta_xxxyy_xxzz_1, ta_xxxyy_xzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyy_xxxx_0[i] = 2.0 * ta1_y_xxxy_xxxx_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xxxx_1[i] * fe_0 + ta_xxxyy_xxxx_1[i] + ta1_y_xxxyy_xxxx_0[i] * pa_y[i] - ta1_y_xxxyy_xxxx_1[i] * pc_y[i];

        ta1_y_xxxyyy_xxxy_0[i] = 2.0 * ta1_y_xyyy_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xxyyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_xxyyy_xxy_1[i] * fe_0 + ta1_y_xxyyy_xxxy_0[i] * pa_x[i] - ta1_y_xxyyy_xxxy_1[i] * pc_x[i];

        ta1_y_xxxyyy_xxxz_0[i] = 2.0 * ta1_y_xxxy_xxxz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xxxz_1[i] * fe_0 + ta_xxxyy_xxxz_1[i] + ta1_y_xxxyy_xxxz_0[i] * pa_y[i] - ta1_y_xxxyy_xxxz_1[i] * pc_y[i];

        ta1_y_xxxyyy_xxyy_0[i] = 2.0 * ta1_y_xyyy_xxyy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xxyyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_xxyyy_xyy_1[i] * fe_0 + ta1_y_xxyyy_xxyy_0[i] * pa_x[i] - ta1_y_xxyyy_xxyy_1[i] * pc_x[i];

        ta1_y_xxxyyy_xxyz_0[i] = 2.0 * ta1_y_xyyy_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxyyy_xyz_1[i] * fe_0 + ta1_y_xxyyy_xxyz_0[i] * pa_x[i] - ta1_y_xxyyy_xxyz_1[i] * pc_x[i];

        ta1_y_xxxyyy_xxzz_0[i] = 2.0 * ta1_y_xxxy_xxzz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xxzz_1[i] * fe_0 + ta_xxxyy_xxzz_1[i] + ta1_y_xxxyy_xxzz_0[i] * pa_y[i] - ta1_y_xxxyy_xxzz_1[i] * pc_y[i];

        ta1_y_xxxyyy_xyyy_0[i] = 2.0 * ta1_y_xyyy_xyyy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xyyy_1[i] * fe_0 + ta1_y_xxyyy_yyy_0[i] * fe_0 - ta1_y_xxyyy_yyy_1[i] * fe_0 + ta1_y_xxyyy_xyyy_0[i] * pa_x[i] - ta1_y_xxyyy_xyyy_1[i] * pc_x[i];

        ta1_y_xxxyyy_xyyz_0[i] = 2.0 * ta1_y_xyyy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xyyz_1[i] * fe_0 + ta1_y_xxyyy_yyz_0[i] * fe_0 - ta1_y_xxyyy_yyz_1[i] * fe_0 + ta1_y_xxyyy_xyyz_0[i] * pa_x[i] - ta1_y_xxyyy_xyyz_1[i] * pc_x[i];

        ta1_y_xxxyyy_xyzz_0[i] = 2.0 * ta1_y_xyyy_xyzz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xyzz_1[i] * fe_0 + ta1_y_xxyyy_yzz_0[i] * fe_0 - ta1_y_xxyyy_yzz_1[i] * fe_0 + ta1_y_xxyyy_xyzz_0[i] * pa_x[i] - ta1_y_xxyyy_xyzz_1[i] * pc_x[i];

        ta1_y_xxxyyy_xzzz_0[i] = 2.0 * ta1_y_xxxy_xzzz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xzzz_1[i] * fe_0 + ta_xxxyy_xzzz_1[i] + ta1_y_xxxyy_xzzz_0[i] * pa_y[i] - ta1_y_xxxyy_xzzz_1[i] * pc_y[i];

        ta1_y_xxxyyy_yyyy_0[i] = 2.0 * ta1_y_xyyy_yyyy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yyyy_1[i] * fe_0 + ta1_y_xxyyy_yyyy_0[i] * pa_x[i] - ta1_y_xxyyy_yyyy_1[i] * pc_x[i];

        ta1_y_xxxyyy_yyyz_0[i] = 2.0 * ta1_y_xyyy_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yyyz_1[i] * fe_0 + ta1_y_xxyyy_yyyz_0[i] * pa_x[i] - ta1_y_xxyyy_yyyz_1[i] * pc_x[i];

        ta1_y_xxxyyy_yyzz_0[i] = 2.0 * ta1_y_xyyy_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yyzz_1[i] * fe_0 + ta1_y_xxyyy_yyzz_0[i] * pa_x[i] - ta1_y_xxyyy_yyzz_1[i] * pc_x[i];

        ta1_y_xxxyyy_yzzz_0[i] = 2.0 * ta1_y_xyyy_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yzzz_1[i] * fe_0 + ta1_y_xxyyy_yzzz_0[i] * pa_x[i] - ta1_y_xxyyy_yzzz_1[i] * pc_x[i];

        ta1_y_xxxyyy_zzzz_0[i] = 2.0 * ta1_y_xyyy_zzzz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_zzzz_1[i] * fe_0 + ta1_y_xxyyy_zzzz_0[i] * pa_x[i] - ta1_y_xxyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 525-540 components of targeted buffer : IG

    auto ta1_y_xxxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 525);

    auto ta1_y_xxxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 526);

    auto ta1_y_xxxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 527);

    auto ta1_y_xxxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 528);

    auto ta1_y_xxxyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 529);

    auto ta1_y_xxxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 530);

    auto ta1_y_xxxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 531);

    auto ta1_y_xxxyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 532);

    auto ta1_y_xxxyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 533);

    auto ta1_y_xxxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 534);

    auto ta1_y_xxxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 535);

    auto ta1_y_xxxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 536);

    auto ta1_y_xxxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 537);

    auto ta1_y_xxxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 538);

    auto ta1_y_xxxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 539);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxyy_xxx_0, ta1_y_xxxyy_xxx_1, ta1_y_xxxyy_xxxx_0, ta1_y_xxxyy_xxxx_1, ta1_y_xxxyy_xxxy_0, ta1_y_xxxyy_xxxy_1, ta1_y_xxxyy_xxxz_0, ta1_y_xxxyy_xxxz_1, ta1_y_xxxyy_xxy_0, ta1_y_xxxyy_xxy_1, ta1_y_xxxyy_xxyy_0, ta1_y_xxxyy_xxyy_1, ta1_y_xxxyy_xxyz_0, ta1_y_xxxyy_xxyz_1, ta1_y_xxxyy_xxz_0, ta1_y_xxxyy_xxz_1, ta1_y_xxxyy_xxzz_0, ta1_y_xxxyy_xxzz_1, ta1_y_xxxyy_xyy_0, ta1_y_xxxyy_xyy_1, ta1_y_xxxyy_xyyy_0, ta1_y_xxxyy_xyyy_1, ta1_y_xxxyy_xyyz_0, ta1_y_xxxyy_xyyz_1, ta1_y_xxxyy_xyz_0, ta1_y_xxxyy_xyz_1, ta1_y_xxxyy_xyzz_0, ta1_y_xxxyy_xyzz_1, ta1_y_xxxyy_xzz_0, ta1_y_xxxyy_xzz_1, ta1_y_xxxyy_xzzz_0, ta1_y_xxxyy_xzzz_1, ta1_y_xxxyy_yyyy_0, ta1_y_xxxyy_yyyy_1, ta1_y_xxxyyz_xxxx_0, ta1_y_xxxyyz_xxxy_0, ta1_y_xxxyyz_xxxz_0, ta1_y_xxxyyz_xxyy_0, ta1_y_xxxyyz_xxyz_0, ta1_y_xxxyyz_xxzz_0, ta1_y_xxxyyz_xyyy_0, ta1_y_xxxyyz_xyyz_0, ta1_y_xxxyyz_xyzz_0, ta1_y_xxxyyz_xzzz_0, ta1_y_xxxyyz_yyyy_0, ta1_y_xxxyyz_yyyz_0, ta1_y_xxxyyz_yyzz_0, ta1_y_xxxyyz_yzzz_0, ta1_y_xxxyyz_zzzz_0, ta1_y_xxyyz_yyyz_0, ta1_y_xxyyz_yyyz_1, ta1_y_xxyyz_yyzz_0, ta1_y_xxyyz_yyzz_1, ta1_y_xxyyz_yzzz_0, ta1_y_xxyyz_yzzz_1, ta1_y_xxyyz_zzzz_0, ta1_y_xxyyz_zzzz_1, ta1_y_xyyz_yyyz_0, ta1_y_xyyz_yyyz_1, ta1_y_xyyz_yyzz_0, ta1_y_xyyz_yyzz_1, ta1_y_xyyz_yzzz_0, ta1_y_xyyz_yzzz_1, ta1_y_xyyz_zzzz_0, ta1_y_xyyz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyz_xxxx_0[i] = ta1_y_xxxyy_xxxx_0[i] * pa_z[i] - ta1_y_xxxyy_xxxx_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxxy_0[i] = ta1_y_xxxyy_xxxy_0[i] * pa_z[i] - ta1_y_xxxyy_xxxy_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxxz_0[i] = ta1_y_xxxyy_xxx_0[i] * fe_0 - ta1_y_xxxyy_xxx_1[i] * fe_0 + ta1_y_xxxyy_xxxz_0[i] * pa_z[i] - ta1_y_xxxyy_xxxz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxyy_0[i] = ta1_y_xxxyy_xxyy_0[i] * pa_z[i] - ta1_y_xxxyy_xxyy_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxyz_0[i] = ta1_y_xxxyy_xxy_0[i] * fe_0 - ta1_y_xxxyy_xxy_1[i] * fe_0 + ta1_y_xxxyy_xxyz_0[i] * pa_z[i] - ta1_y_xxxyy_xxyz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxzz_0[i] = 2.0 * ta1_y_xxxyy_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxxyy_xxz_1[i] * fe_0 + ta1_y_xxxyy_xxzz_0[i] * pa_z[i] - ta1_y_xxxyy_xxzz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xyyy_0[i] = ta1_y_xxxyy_xyyy_0[i] * pa_z[i] - ta1_y_xxxyy_xyyy_1[i] * pc_z[i];

        ta1_y_xxxyyz_xyyz_0[i] = ta1_y_xxxyy_xyy_0[i] * fe_0 - ta1_y_xxxyy_xyy_1[i] * fe_0 + ta1_y_xxxyy_xyyz_0[i] * pa_z[i] - ta1_y_xxxyy_xyyz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xyzz_0[i] = 2.0 * ta1_y_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxyy_xyz_1[i] * fe_0 + ta1_y_xxxyy_xyzz_0[i] * pa_z[i] - ta1_y_xxxyy_xyzz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xzzz_0[i] = 3.0 * ta1_y_xxxyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxxyy_xzz_1[i] * fe_0 + ta1_y_xxxyy_xzzz_0[i] * pa_z[i] - ta1_y_xxxyy_xzzz_1[i] * pc_z[i];

        ta1_y_xxxyyz_yyyy_0[i] = ta1_y_xxxyy_yyyy_0[i] * pa_z[i] - ta1_y_xxxyy_yyyy_1[i] * pc_z[i];

        ta1_y_xxxyyz_yyyz_0[i] = 2.0 * ta1_y_xyyz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_yyyz_1[i] * fe_0 + ta1_y_xxyyz_yyyz_0[i] * pa_x[i] - ta1_y_xxyyz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxyyz_yyzz_0[i] = 2.0 * ta1_y_xyyz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_yyzz_1[i] * fe_0 + ta1_y_xxyyz_yyzz_0[i] * pa_x[i] - ta1_y_xxyyz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxyyz_yzzz_0[i] = 2.0 * ta1_y_xyyz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_yzzz_1[i] * fe_0 + ta1_y_xxyyz_yzzz_0[i] * pa_x[i] - ta1_y_xxyyz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxyyz_zzzz_0[i] = 2.0 * ta1_y_xyyz_zzzz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_zzzz_1[i] * fe_0 + ta1_y_xxyyz_zzzz_0[i] * pa_x[i] - ta1_y_xxyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 540-555 components of targeted buffer : IG

    auto ta1_y_xxxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 540);

    auto ta1_y_xxxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 541);

    auto ta1_y_xxxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 542);

    auto ta1_y_xxxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 543);

    auto ta1_y_xxxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 544);

    auto ta1_y_xxxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 545);

    auto ta1_y_xxxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 546);

    auto ta1_y_xxxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 547);

    auto ta1_y_xxxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 548);

    auto ta1_y_xxxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 549);

    auto ta1_y_xxxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 550);

    auto ta1_y_xxxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 551);

    auto ta1_y_xxxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 552);

    auto ta1_y_xxxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 553);

    auto ta1_y_xxxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 554);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxxy_xxxy_0, ta1_y_xxxy_xxxy_1, ta1_y_xxxy_xxyy_0, ta1_y_xxxy_xxyy_1, ta1_y_xxxy_xyyy_0, ta1_y_xxxy_xyyy_1, ta1_y_xxxyz_xxxy_0, ta1_y_xxxyz_xxxy_1, ta1_y_xxxyz_xxyy_0, ta1_y_xxxyz_xxyy_1, ta1_y_xxxyz_xyyy_0, ta1_y_xxxyz_xyyy_1, ta1_y_xxxyzz_xxxx_0, ta1_y_xxxyzz_xxxy_0, ta1_y_xxxyzz_xxxz_0, ta1_y_xxxyzz_xxyy_0, ta1_y_xxxyzz_xxyz_0, ta1_y_xxxyzz_xxzz_0, ta1_y_xxxyzz_xyyy_0, ta1_y_xxxyzz_xyyz_0, ta1_y_xxxyzz_xyzz_0, ta1_y_xxxyzz_xzzz_0, ta1_y_xxxyzz_yyyy_0, ta1_y_xxxyzz_yyyz_0, ta1_y_xxxyzz_yyzz_0, ta1_y_xxxyzz_yzzz_0, ta1_y_xxxyzz_zzzz_0, ta1_y_xxxzz_xxxx_0, ta1_y_xxxzz_xxxx_1, ta1_y_xxxzz_xxxz_0, ta1_y_xxxzz_xxxz_1, ta1_y_xxxzz_xxyz_0, ta1_y_xxxzz_xxyz_1, ta1_y_xxxzz_xxz_0, ta1_y_xxxzz_xxz_1, ta1_y_xxxzz_xxzz_0, ta1_y_xxxzz_xxzz_1, ta1_y_xxxzz_xyyz_0, ta1_y_xxxzz_xyyz_1, ta1_y_xxxzz_xyz_0, ta1_y_xxxzz_xyz_1, ta1_y_xxxzz_xyzz_0, ta1_y_xxxzz_xyzz_1, ta1_y_xxxzz_xzz_0, ta1_y_xxxzz_xzz_1, ta1_y_xxxzz_xzzz_0, ta1_y_xxxzz_xzzz_1, ta1_y_xxxzz_zzzz_0, ta1_y_xxxzz_zzzz_1, ta1_y_xxyzz_yyyy_0, ta1_y_xxyzz_yyyy_1, ta1_y_xxyzz_yyyz_0, ta1_y_xxyzz_yyyz_1, ta1_y_xxyzz_yyzz_0, ta1_y_xxyzz_yyzz_1, ta1_y_xxyzz_yzzz_0, ta1_y_xxyzz_yzzz_1, ta1_y_xyzz_yyyy_0, ta1_y_xyzz_yyyy_1, ta1_y_xyzz_yyyz_0, ta1_y_xyzz_yyyz_1, ta1_y_xyzz_yyzz_0, ta1_y_xyzz_yyzz_1, ta1_y_xyzz_yzzz_0, ta1_y_xyzz_yzzz_1, ta_xxxzz_xxxx_1, ta_xxxzz_xxxz_1, ta_xxxzz_xxyz_1, ta_xxxzz_xxzz_1, ta_xxxzz_xyyz_1, ta_xxxzz_xyzz_1, ta_xxxzz_xzzz_1, ta_xxxzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyzz_xxxx_0[i] = ta_xxxzz_xxxx_1[i] + ta1_y_xxxzz_xxxx_0[i] * pa_y[i] - ta1_y_xxxzz_xxxx_1[i] * pc_y[i];

        ta1_y_xxxyzz_xxxy_0[i] = ta1_y_xxxy_xxxy_0[i] * fe_0 - ta1_y_xxxy_xxxy_1[i] * fe_0 + ta1_y_xxxyz_xxxy_0[i] * pa_z[i] - ta1_y_xxxyz_xxxy_1[i] * pc_z[i];

        ta1_y_xxxyzz_xxxz_0[i] = ta_xxxzz_xxxz_1[i] + ta1_y_xxxzz_xxxz_0[i] * pa_y[i] - ta1_y_xxxzz_xxxz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xxyy_0[i] = ta1_y_xxxy_xxyy_0[i] * fe_0 - ta1_y_xxxy_xxyy_1[i] * fe_0 + ta1_y_xxxyz_xxyy_0[i] * pa_z[i] - ta1_y_xxxyz_xxyy_1[i] * pc_z[i];

        ta1_y_xxxyzz_xxyz_0[i] = ta1_y_xxxzz_xxz_0[i] * fe_0 - ta1_y_xxxzz_xxz_1[i] * fe_0 + ta_xxxzz_xxyz_1[i] + ta1_y_xxxzz_xxyz_0[i] * pa_y[i] - ta1_y_xxxzz_xxyz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xxzz_0[i] = ta_xxxzz_xxzz_1[i] + ta1_y_xxxzz_xxzz_0[i] * pa_y[i] - ta1_y_xxxzz_xxzz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xyyy_0[i] = ta1_y_xxxy_xyyy_0[i] * fe_0 - ta1_y_xxxy_xyyy_1[i] * fe_0 + ta1_y_xxxyz_xyyy_0[i] * pa_z[i] - ta1_y_xxxyz_xyyy_1[i] * pc_z[i];

        ta1_y_xxxyzz_xyyz_0[i] = 2.0 * ta1_y_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxxzz_xyz_1[i] * fe_0 + ta_xxxzz_xyyz_1[i] + ta1_y_xxxzz_xyyz_0[i] * pa_y[i] - ta1_y_xxxzz_xyyz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xyzz_0[i] = ta1_y_xxxzz_xzz_0[i] * fe_0 - ta1_y_xxxzz_xzz_1[i] * fe_0 + ta_xxxzz_xyzz_1[i] + ta1_y_xxxzz_xyzz_0[i] * pa_y[i] - ta1_y_xxxzz_xyzz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xzzz_0[i] = ta_xxxzz_xzzz_1[i] + ta1_y_xxxzz_xzzz_0[i] * pa_y[i] - ta1_y_xxxzz_xzzz_1[i] * pc_y[i];

        ta1_y_xxxyzz_yyyy_0[i] = 2.0 * ta1_y_xyzz_yyyy_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yyyy_1[i] * fe_0 + ta1_y_xxyzz_yyyy_0[i] * pa_x[i] - ta1_y_xxyzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxxyzz_yyyz_0[i] = 2.0 * ta1_y_xyzz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yyyz_1[i] * fe_0 + ta1_y_xxyzz_yyyz_0[i] * pa_x[i] - ta1_y_xxyzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxyzz_yyzz_0[i] = 2.0 * ta1_y_xyzz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yyzz_1[i] * fe_0 + ta1_y_xxyzz_yyzz_0[i] * pa_x[i] - ta1_y_xxyzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxyzz_yzzz_0[i] = 2.0 * ta1_y_xyzz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yzzz_1[i] * fe_0 + ta1_y_xxyzz_yzzz_0[i] * pa_x[i] - ta1_y_xxyzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxyzz_zzzz_0[i] = ta_xxxzz_zzzz_1[i] + ta1_y_xxxzz_zzzz_0[i] * pa_y[i] - ta1_y_xxxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 555-570 components of targeted buffer : IG

    auto ta1_y_xxxzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 555);

    auto ta1_y_xxxzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 556);

    auto ta1_y_xxxzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 557);

    auto ta1_y_xxxzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 558);

    auto ta1_y_xxxzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 559);

    auto ta1_y_xxxzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 560);

    auto ta1_y_xxxzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 561);

    auto ta1_y_xxxzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 562);

    auto ta1_y_xxxzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 563);

    auto ta1_y_xxxzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 564);

    auto ta1_y_xxxzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 565);

    auto ta1_y_xxxzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 566);

    auto ta1_y_xxxzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 567);

    auto ta1_y_xxxzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 568);

    auto ta1_y_xxxzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 569);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxz_xxxx_0, ta1_y_xxxz_xxxx_1, ta1_y_xxxz_xxxy_0, ta1_y_xxxz_xxxy_1, ta1_y_xxxz_xxyy_0, ta1_y_xxxz_xxyy_1, ta1_y_xxxz_xyyy_0, ta1_y_xxxz_xyyy_1, ta1_y_xxxzz_xxxx_0, ta1_y_xxxzz_xxxx_1, ta1_y_xxxzz_xxxy_0, ta1_y_xxxzz_xxxy_1, ta1_y_xxxzz_xxyy_0, ta1_y_xxxzz_xxyy_1, ta1_y_xxxzz_xyyy_0, ta1_y_xxxzz_xyyy_1, ta1_y_xxxzzz_xxxx_0, ta1_y_xxxzzz_xxxy_0, ta1_y_xxxzzz_xxxz_0, ta1_y_xxxzzz_xxyy_0, ta1_y_xxxzzz_xxyz_0, ta1_y_xxxzzz_xxzz_0, ta1_y_xxxzzz_xyyy_0, ta1_y_xxxzzz_xyyz_0, ta1_y_xxxzzz_xyzz_0, ta1_y_xxxzzz_xzzz_0, ta1_y_xxxzzz_yyyy_0, ta1_y_xxxzzz_yyyz_0, ta1_y_xxxzzz_yyzz_0, ta1_y_xxxzzz_yzzz_0, ta1_y_xxxzzz_zzzz_0, ta1_y_xxzzz_xxxz_0, ta1_y_xxzzz_xxxz_1, ta1_y_xxzzz_xxyz_0, ta1_y_xxzzz_xxyz_1, ta1_y_xxzzz_xxz_0, ta1_y_xxzzz_xxz_1, ta1_y_xxzzz_xxzz_0, ta1_y_xxzzz_xxzz_1, ta1_y_xxzzz_xyyz_0, ta1_y_xxzzz_xyyz_1, ta1_y_xxzzz_xyz_0, ta1_y_xxzzz_xyz_1, ta1_y_xxzzz_xyzz_0, ta1_y_xxzzz_xyzz_1, ta1_y_xxzzz_xzz_0, ta1_y_xxzzz_xzz_1, ta1_y_xxzzz_xzzz_0, ta1_y_xxzzz_xzzz_1, ta1_y_xxzzz_yyyy_0, ta1_y_xxzzz_yyyy_1, ta1_y_xxzzz_yyyz_0, ta1_y_xxzzz_yyyz_1, ta1_y_xxzzz_yyz_0, ta1_y_xxzzz_yyz_1, ta1_y_xxzzz_yyzz_0, ta1_y_xxzzz_yyzz_1, ta1_y_xxzzz_yzz_0, ta1_y_xxzzz_yzz_1, ta1_y_xxzzz_yzzz_0, ta1_y_xxzzz_yzzz_1, ta1_y_xxzzz_zzz_0, ta1_y_xxzzz_zzz_1, ta1_y_xxzzz_zzzz_0, ta1_y_xxzzz_zzzz_1, ta1_y_xzzz_xxxz_0, ta1_y_xzzz_xxxz_1, ta1_y_xzzz_xxyz_0, ta1_y_xzzz_xxyz_1, ta1_y_xzzz_xxzz_0, ta1_y_xzzz_xxzz_1, ta1_y_xzzz_xyyz_0, ta1_y_xzzz_xyyz_1, ta1_y_xzzz_xyzz_0, ta1_y_xzzz_xyzz_1, ta1_y_xzzz_xzzz_0, ta1_y_xzzz_xzzz_1, ta1_y_xzzz_yyyy_0, ta1_y_xzzz_yyyy_1, ta1_y_xzzz_yyyz_0, ta1_y_xzzz_yyyz_1, ta1_y_xzzz_yyzz_0, ta1_y_xzzz_yyzz_1, ta1_y_xzzz_yzzz_0, ta1_y_xzzz_yzzz_1, ta1_y_xzzz_zzzz_0, ta1_y_xzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzzz_xxxx_0[i] = 2.0 * ta1_y_xxxz_xxxx_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xxxx_1[i] * fe_0 + ta1_y_xxxzz_xxxx_0[i] * pa_z[i] - ta1_y_xxxzz_xxxx_1[i] * pc_z[i];

        ta1_y_xxxzzz_xxxy_0[i] = 2.0 * ta1_y_xxxz_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xxxy_1[i] * fe_0 + ta1_y_xxxzz_xxxy_0[i] * pa_z[i] - ta1_y_xxxzz_xxxy_1[i] * pc_z[i];

        ta1_y_xxxzzz_xxxz_0[i] = 2.0 * ta1_y_xzzz_xxxz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xxzzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_xxzzz_xxz_1[i] * fe_0 + ta1_y_xxzzz_xxxz_0[i] * pa_x[i] - ta1_y_xxzzz_xxxz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xxyy_0[i] = 2.0 * ta1_y_xxxz_xxyy_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xxyy_1[i] * fe_0 + ta1_y_xxxzz_xxyy_0[i] * pa_z[i] - ta1_y_xxxzz_xxyy_1[i] * pc_z[i];

        ta1_y_xxxzzz_xxyz_0[i] = 2.0 * ta1_y_xzzz_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxzzz_xyz_1[i] * fe_0 + ta1_y_xxzzz_xxyz_0[i] * pa_x[i] - ta1_y_xxzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xxzz_0[i] = 2.0 * ta1_y_xzzz_xxzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_xxzzz_xzz_1[i] * fe_0 + ta1_y_xxzzz_xxzz_0[i] * pa_x[i] - ta1_y_xxzzz_xxzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xyyy_0[i] = 2.0 * ta1_y_xxxz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xyyy_1[i] * fe_0 + ta1_y_xxxzz_xyyy_0[i] * pa_z[i] - ta1_y_xxxzz_xyyy_1[i] * pc_z[i];

        ta1_y_xxxzzz_xyyz_0[i] = 2.0 * ta1_y_xzzz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xyyz_1[i] * fe_0 + ta1_y_xxzzz_yyz_0[i] * fe_0 - ta1_y_xxzzz_yyz_1[i] * fe_0 + ta1_y_xxzzz_xyyz_0[i] * pa_x[i] - ta1_y_xxzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xyzz_0[i] = 2.0 * ta1_y_xzzz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xyzz_1[i] * fe_0 + ta1_y_xxzzz_yzz_0[i] * fe_0 - ta1_y_xxzzz_yzz_1[i] * fe_0 + ta1_y_xxzzz_xyzz_0[i] * pa_x[i] - ta1_y_xxzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xzzz_0[i] = 2.0 * ta1_y_xzzz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xzzz_1[i] * fe_0 + ta1_y_xxzzz_zzz_0[i] * fe_0 - ta1_y_xxzzz_zzz_1[i] * fe_0 + ta1_y_xxzzz_xzzz_0[i] * pa_x[i] - ta1_y_xxzzz_xzzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_yyyy_0[i] = 2.0 * ta1_y_xzzz_yyyy_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yyyy_1[i] * fe_0 + ta1_y_xxzzz_yyyy_0[i] * pa_x[i] - ta1_y_xxzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxxzzz_yyyz_0[i] = 2.0 * ta1_y_xzzz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yyyz_1[i] * fe_0 + ta1_y_xxzzz_yyyz_0[i] * pa_x[i] - ta1_y_xxzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxxzzz_yyzz_0[i] = 2.0 * ta1_y_xzzz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yyzz_1[i] * fe_0 + ta1_y_xxzzz_yyzz_0[i] * pa_x[i] - ta1_y_xxzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_yzzz_0[i] = 2.0 * ta1_y_xzzz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yzzz_1[i] * fe_0 + ta1_y_xxzzz_yzzz_0[i] * pa_x[i] - ta1_y_xxzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_zzzz_0[i] = 2.0 * ta1_y_xzzz_zzzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_zzzz_1[i] * fe_0 + ta1_y_xxzzz_zzzz_0[i] * pa_x[i] - ta1_y_xxzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 570-585 components of targeted buffer : IG

    auto ta1_y_xxyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 570);

    auto ta1_y_xxyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 571);

    auto ta1_y_xxyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 572);

    auto ta1_y_xxyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 573);

    auto ta1_y_xxyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 574);

    auto ta1_y_xxyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 575);

    auto ta1_y_xxyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 576);

    auto ta1_y_xxyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 577);

    auto ta1_y_xxyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 578);

    auto ta1_y_xxyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 579);

    auto ta1_y_xxyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 580);

    auto ta1_y_xxyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 581);

    auto ta1_y_xxyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 582);

    auto ta1_y_xxyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 583);

    auto ta1_y_xxyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 584);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxyy_xxxx_0, ta1_y_xxyy_xxxx_1, ta1_y_xxyy_xxxz_0, ta1_y_xxyy_xxxz_1, ta1_y_xxyy_xxzz_0, ta1_y_xxyy_xxzz_1, ta1_y_xxyy_xzzz_0, ta1_y_xxyy_xzzz_1, ta1_y_xxyyy_xxxx_0, ta1_y_xxyyy_xxxx_1, ta1_y_xxyyy_xxxz_0, ta1_y_xxyyy_xxxz_1, ta1_y_xxyyy_xxzz_0, ta1_y_xxyyy_xxzz_1, ta1_y_xxyyy_xzzz_0, ta1_y_xxyyy_xzzz_1, ta1_y_xxyyyy_xxxx_0, ta1_y_xxyyyy_xxxy_0, ta1_y_xxyyyy_xxxz_0, ta1_y_xxyyyy_xxyy_0, ta1_y_xxyyyy_xxyz_0, ta1_y_xxyyyy_xxzz_0, ta1_y_xxyyyy_xyyy_0, ta1_y_xxyyyy_xyyz_0, ta1_y_xxyyyy_xyzz_0, ta1_y_xxyyyy_xzzz_0, ta1_y_xxyyyy_yyyy_0, ta1_y_xxyyyy_yyyz_0, ta1_y_xxyyyy_yyzz_0, ta1_y_xxyyyy_yzzz_0, ta1_y_xxyyyy_zzzz_0, ta1_y_xyyyy_xxxy_0, ta1_y_xyyyy_xxxy_1, ta1_y_xyyyy_xxy_0, ta1_y_xyyyy_xxy_1, ta1_y_xyyyy_xxyy_0, ta1_y_xyyyy_xxyy_1, ta1_y_xyyyy_xxyz_0, ta1_y_xyyyy_xxyz_1, ta1_y_xyyyy_xyy_0, ta1_y_xyyyy_xyy_1, ta1_y_xyyyy_xyyy_0, ta1_y_xyyyy_xyyy_1, ta1_y_xyyyy_xyyz_0, ta1_y_xyyyy_xyyz_1, ta1_y_xyyyy_xyz_0, ta1_y_xyyyy_xyz_1, ta1_y_xyyyy_xyzz_0, ta1_y_xyyyy_xyzz_1, ta1_y_xyyyy_yyy_0, ta1_y_xyyyy_yyy_1, ta1_y_xyyyy_yyyy_0, ta1_y_xyyyy_yyyy_1, ta1_y_xyyyy_yyyz_0, ta1_y_xyyyy_yyyz_1, ta1_y_xyyyy_yyz_0, ta1_y_xyyyy_yyz_1, ta1_y_xyyyy_yyzz_0, ta1_y_xyyyy_yyzz_1, ta1_y_xyyyy_yzz_0, ta1_y_xyyyy_yzz_1, ta1_y_xyyyy_yzzz_0, ta1_y_xyyyy_yzzz_1, ta1_y_xyyyy_zzzz_0, ta1_y_xyyyy_zzzz_1, ta1_y_yyyy_xxxy_0, ta1_y_yyyy_xxxy_1, ta1_y_yyyy_xxyy_0, ta1_y_yyyy_xxyy_1, ta1_y_yyyy_xxyz_0, ta1_y_yyyy_xxyz_1, ta1_y_yyyy_xyyy_0, ta1_y_yyyy_xyyy_1, ta1_y_yyyy_xyyz_0, ta1_y_yyyy_xyyz_1, ta1_y_yyyy_xyzz_0, ta1_y_yyyy_xyzz_1, ta1_y_yyyy_yyyy_0, ta1_y_yyyy_yyyy_1, ta1_y_yyyy_yyyz_0, ta1_y_yyyy_yyyz_1, ta1_y_yyyy_yyzz_0, ta1_y_yyyy_yyzz_1, ta1_y_yyyy_yzzz_0, ta1_y_yyyy_yzzz_1, ta1_y_yyyy_zzzz_0, ta1_y_yyyy_zzzz_1, ta_xxyyy_xxxx_1, ta_xxyyy_xxxz_1, ta_xxyyy_xxzz_1, ta_xxyyy_xzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyy_xxxx_0[i] = 3.0 * ta1_y_xxyy_xxxx_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxxx_1[i] * fe_0 + ta_xxyyy_xxxx_1[i] + ta1_y_xxyyy_xxxx_0[i] * pa_y[i] - ta1_y_xxyyy_xxxx_1[i] * pc_y[i];

        ta1_y_xxyyyy_xxxy_0[i] = ta1_y_yyyy_xxxy_0[i] * fe_0 - ta1_y_yyyy_xxxy_1[i] * fe_0 + 3.0 * ta1_y_xyyyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_xyyyy_xxy_1[i] * fe_0 + ta1_y_xyyyy_xxxy_0[i] * pa_x[i] - ta1_y_xyyyy_xxxy_1[i] * pc_x[i];

        ta1_y_xxyyyy_xxxz_0[i] = 3.0 * ta1_y_xxyy_xxxz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxxz_1[i] * fe_0 + ta_xxyyy_xxxz_1[i] + ta1_y_xxyyy_xxxz_0[i] * pa_y[i] - ta1_y_xxyyy_xxxz_1[i] * pc_y[i];

        ta1_y_xxyyyy_xxyy_0[i] = ta1_y_yyyy_xxyy_0[i] * fe_0 - ta1_y_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_xyyyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_xyyyy_xyy_1[i] * fe_0 + ta1_y_xyyyy_xxyy_0[i] * pa_x[i] - ta1_y_xyyyy_xxyy_1[i] * pc_x[i];

        ta1_y_xxyyyy_xxyz_0[i] = ta1_y_yyyy_xxyz_0[i] * fe_0 - ta1_y_yyyy_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xyyyy_xyz_1[i] * fe_0 + ta1_y_xyyyy_xxyz_0[i] * pa_x[i] - ta1_y_xyyyy_xxyz_1[i] * pc_x[i];

        ta1_y_xxyyyy_xxzz_0[i] = 3.0 * ta1_y_xxyy_xxzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxzz_1[i] * fe_0 + ta_xxyyy_xxzz_1[i] + ta1_y_xxyyy_xxzz_0[i] * pa_y[i] - ta1_y_xxyyy_xxzz_1[i] * pc_y[i];

        ta1_y_xxyyyy_xyyy_0[i] = ta1_y_yyyy_xyyy_0[i] * fe_0 - ta1_y_yyyy_xyyy_1[i] * fe_0 + ta1_y_xyyyy_yyy_0[i] * fe_0 - ta1_y_xyyyy_yyy_1[i] * fe_0 + ta1_y_xyyyy_xyyy_0[i] * pa_x[i] - ta1_y_xyyyy_xyyy_1[i] * pc_x[i];

        ta1_y_xxyyyy_xyyz_0[i] = ta1_y_yyyy_xyyz_0[i] * fe_0 - ta1_y_yyyy_xyyz_1[i] * fe_0 + ta1_y_xyyyy_yyz_0[i] * fe_0 - ta1_y_xyyyy_yyz_1[i] * fe_0 + ta1_y_xyyyy_xyyz_0[i] * pa_x[i] - ta1_y_xyyyy_xyyz_1[i] * pc_x[i];

        ta1_y_xxyyyy_xyzz_0[i] = ta1_y_yyyy_xyzz_0[i] * fe_0 - ta1_y_yyyy_xyzz_1[i] * fe_0 + ta1_y_xyyyy_yzz_0[i] * fe_0 - ta1_y_xyyyy_yzz_1[i] * fe_0 + ta1_y_xyyyy_xyzz_0[i] * pa_x[i] - ta1_y_xyyyy_xyzz_1[i] * pc_x[i];

        ta1_y_xxyyyy_xzzz_0[i] = 3.0 * ta1_y_xxyy_xzzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xzzz_1[i] * fe_0 + ta_xxyyy_xzzz_1[i] + ta1_y_xxyyy_xzzz_0[i] * pa_y[i] - ta1_y_xxyyy_xzzz_1[i] * pc_y[i];

        ta1_y_xxyyyy_yyyy_0[i] = ta1_y_yyyy_yyyy_0[i] * fe_0 - ta1_y_yyyy_yyyy_1[i] * fe_0 + ta1_y_xyyyy_yyyy_0[i] * pa_x[i] - ta1_y_xyyyy_yyyy_1[i] * pc_x[i];

        ta1_y_xxyyyy_yyyz_0[i] = ta1_y_yyyy_yyyz_0[i] * fe_0 - ta1_y_yyyy_yyyz_1[i] * fe_0 + ta1_y_xyyyy_yyyz_0[i] * pa_x[i] - ta1_y_xyyyy_yyyz_1[i] * pc_x[i];

        ta1_y_xxyyyy_yyzz_0[i] = ta1_y_yyyy_yyzz_0[i] * fe_0 - ta1_y_yyyy_yyzz_1[i] * fe_0 + ta1_y_xyyyy_yyzz_0[i] * pa_x[i] - ta1_y_xyyyy_yyzz_1[i] * pc_x[i];

        ta1_y_xxyyyy_yzzz_0[i] = ta1_y_yyyy_yzzz_0[i] * fe_0 - ta1_y_yyyy_yzzz_1[i] * fe_0 + ta1_y_xyyyy_yzzz_0[i] * pa_x[i] - ta1_y_xyyyy_yzzz_1[i] * pc_x[i];

        ta1_y_xxyyyy_zzzz_0[i] = ta1_y_yyyy_zzzz_0[i] * fe_0 - ta1_y_yyyy_zzzz_1[i] * fe_0 + ta1_y_xyyyy_zzzz_0[i] * pa_x[i] - ta1_y_xyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 585-600 components of targeted buffer : IG

    auto ta1_y_xxyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 585);

    auto ta1_y_xxyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 586);

    auto ta1_y_xxyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 587);

    auto ta1_y_xxyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 588);

    auto ta1_y_xxyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 589);

    auto ta1_y_xxyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 590);

    auto ta1_y_xxyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 591);

    auto ta1_y_xxyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 592);

    auto ta1_y_xxyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 593);

    auto ta1_y_xxyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 594);

    auto ta1_y_xxyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 595);

    auto ta1_y_xxyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 596);

    auto ta1_y_xxyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 597);

    auto ta1_y_xxyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 598);

    auto ta1_y_xxyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 599);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxyyy_xxx_0, ta1_y_xxyyy_xxx_1, ta1_y_xxyyy_xxxx_0, ta1_y_xxyyy_xxxx_1, ta1_y_xxyyy_xxxy_0, ta1_y_xxyyy_xxxy_1, ta1_y_xxyyy_xxxz_0, ta1_y_xxyyy_xxxz_1, ta1_y_xxyyy_xxy_0, ta1_y_xxyyy_xxy_1, ta1_y_xxyyy_xxyy_0, ta1_y_xxyyy_xxyy_1, ta1_y_xxyyy_xxyz_0, ta1_y_xxyyy_xxyz_1, ta1_y_xxyyy_xxz_0, ta1_y_xxyyy_xxz_1, ta1_y_xxyyy_xxzz_0, ta1_y_xxyyy_xxzz_1, ta1_y_xxyyy_xyy_0, ta1_y_xxyyy_xyy_1, ta1_y_xxyyy_xyyy_0, ta1_y_xxyyy_xyyy_1, ta1_y_xxyyy_xyyz_0, ta1_y_xxyyy_xyyz_1, ta1_y_xxyyy_xyz_0, ta1_y_xxyyy_xyz_1, ta1_y_xxyyy_xyzz_0, ta1_y_xxyyy_xyzz_1, ta1_y_xxyyy_xzz_0, ta1_y_xxyyy_xzz_1, ta1_y_xxyyy_xzzz_0, ta1_y_xxyyy_xzzz_1, ta1_y_xxyyy_yyyy_0, ta1_y_xxyyy_yyyy_1, ta1_y_xxyyyz_xxxx_0, ta1_y_xxyyyz_xxxy_0, ta1_y_xxyyyz_xxxz_0, ta1_y_xxyyyz_xxyy_0, ta1_y_xxyyyz_xxyz_0, ta1_y_xxyyyz_xxzz_0, ta1_y_xxyyyz_xyyy_0, ta1_y_xxyyyz_xyyz_0, ta1_y_xxyyyz_xyzz_0, ta1_y_xxyyyz_xzzz_0, ta1_y_xxyyyz_yyyy_0, ta1_y_xxyyyz_yyyz_0, ta1_y_xxyyyz_yyzz_0, ta1_y_xxyyyz_yzzz_0, ta1_y_xxyyyz_zzzz_0, ta1_y_xyyyz_yyyz_0, ta1_y_xyyyz_yyyz_1, ta1_y_xyyyz_yyzz_0, ta1_y_xyyyz_yyzz_1, ta1_y_xyyyz_yzzz_0, ta1_y_xyyyz_yzzz_1, ta1_y_xyyyz_zzzz_0, ta1_y_xyyyz_zzzz_1, ta1_y_yyyz_yyyz_0, ta1_y_yyyz_yyyz_1, ta1_y_yyyz_yyzz_0, ta1_y_yyyz_yyzz_1, ta1_y_yyyz_yzzz_0, ta1_y_yyyz_yzzz_1, ta1_y_yyyz_zzzz_0, ta1_y_yyyz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyz_xxxx_0[i] = ta1_y_xxyyy_xxxx_0[i] * pa_z[i] - ta1_y_xxyyy_xxxx_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxxy_0[i] = ta1_y_xxyyy_xxxy_0[i] * pa_z[i] - ta1_y_xxyyy_xxxy_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxxz_0[i] = ta1_y_xxyyy_xxx_0[i] * fe_0 - ta1_y_xxyyy_xxx_1[i] * fe_0 + ta1_y_xxyyy_xxxz_0[i] * pa_z[i] - ta1_y_xxyyy_xxxz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxyy_0[i] = ta1_y_xxyyy_xxyy_0[i] * pa_z[i] - ta1_y_xxyyy_xxyy_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxyz_0[i] = ta1_y_xxyyy_xxy_0[i] * fe_0 - ta1_y_xxyyy_xxy_1[i] * fe_0 + ta1_y_xxyyy_xxyz_0[i] * pa_z[i] - ta1_y_xxyyy_xxyz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxzz_0[i] = 2.0 * ta1_y_xxyyy_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxyyy_xxz_1[i] * fe_0 + ta1_y_xxyyy_xxzz_0[i] * pa_z[i] - ta1_y_xxyyy_xxzz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xyyy_0[i] = ta1_y_xxyyy_xyyy_0[i] * pa_z[i] - ta1_y_xxyyy_xyyy_1[i] * pc_z[i];

        ta1_y_xxyyyz_xyyz_0[i] = ta1_y_xxyyy_xyy_0[i] * fe_0 - ta1_y_xxyyy_xyy_1[i] * fe_0 + ta1_y_xxyyy_xyyz_0[i] * pa_z[i] - ta1_y_xxyyy_xyyz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xyzz_0[i] = 2.0 * ta1_y_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxyyy_xyz_1[i] * fe_0 + ta1_y_xxyyy_xyzz_0[i] * pa_z[i] - ta1_y_xxyyy_xyzz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xzzz_0[i] = 3.0 * ta1_y_xxyyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxyyy_xzz_1[i] * fe_0 + ta1_y_xxyyy_xzzz_0[i] * pa_z[i] - ta1_y_xxyyy_xzzz_1[i] * pc_z[i];

        ta1_y_xxyyyz_yyyy_0[i] = ta1_y_xxyyy_yyyy_0[i] * pa_z[i] - ta1_y_xxyyy_yyyy_1[i] * pc_z[i];

        ta1_y_xxyyyz_yyyz_0[i] = ta1_y_yyyz_yyyz_0[i] * fe_0 - ta1_y_yyyz_yyyz_1[i] * fe_0 + ta1_y_xyyyz_yyyz_0[i] * pa_x[i] - ta1_y_xyyyz_yyyz_1[i] * pc_x[i];

        ta1_y_xxyyyz_yyzz_0[i] = ta1_y_yyyz_yyzz_0[i] * fe_0 - ta1_y_yyyz_yyzz_1[i] * fe_0 + ta1_y_xyyyz_yyzz_0[i] * pa_x[i] - ta1_y_xyyyz_yyzz_1[i] * pc_x[i];

        ta1_y_xxyyyz_yzzz_0[i] = ta1_y_yyyz_yzzz_0[i] * fe_0 - ta1_y_yyyz_yzzz_1[i] * fe_0 + ta1_y_xyyyz_yzzz_0[i] * pa_x[i] - ta1_y_xyyyz_yzzz_1[i] * pc_x[i];

        ta1_y_xxyyyz_zzzz_0[i] = ta1_y_yyyz_zzzz_0[i] * fe_0 - ta1_y_yyyz_zzzz_1[i] * fe_0 + ta1_y_xyyyz_zzzz_0[i] * pa_x[i] - ta1_y_xyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 600-615 components of targeted buffer : IG

    auto ta1_y_xxyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 600);

    auto ta1_y_xxyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 601);

    auto ta1_y_xxyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 602);

    auto ta1_y_xxyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 603);

    auto ta1_y_xxyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 604);

    auto ta1_y_xxyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 605);

    auto ta1_y_xxyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 606);

    auto ta1_y_xxyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 607);

    auto ta1_y_xxyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 608);

    auto ta1_y_xxyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 609);

    auto ta1_y_xxyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 610);

    auto ta1_y_xxyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 611);

    auto ta1_y_xxyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 612);

    auto ta1_y_xxyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 613);

    auto ta1_y_xxyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 614);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxyy_xxxx_0, ta1_y_xxyy_xxxx_1, ta1_y_xxyy_xxxy_0, ta1_y_xxyy_xxxy_1, ta1_y_xxyy_xxyy_0, ta1_y_xxyy_xxyy_1, ta1_y_xxyy_xyyy_0, ta1_y_xxyy_xyyy_1, ta1_y_xxyyz_xxxx_0, ta1_y_xxyyz_xxxx_1, ta1_y_xxyyz_xxxy_0, ta1_y_xxyyz_xxxy_1, ta1_y_xxyyz_xxyy_0, ta1_y_xxyyz_xxyy_1, ta1_y_xxyyz_xyyy_0, ta1_y_xxyyz_xyyy_1, ta1_y_xxyyzz_xxxx_0, ta1_y_xxyyzz_xxxy_0, ta1_y_xxyyzz_xxxz_0, ta1_y_xxyyzz_xxyy_0, ta1_y_xxyyzz_xxyz_0, ta1_y_xxyyzz_xxzz_0, ta1_y_xxyyzz_xyyy_0, ta1_y_xxyyzz_xyyz_0, ta1_y_xxyyzz_xyzz_0, ta1_y_xxyyzz_xzzz_0, ta1_y_xxyyzz_yyyy_0, ta1_y_xxyyzz_yyyz_0, ta1_y_xxyyzz_yyzz_0, ta1_y_xxyyzz_yzzz_0, ta1_y_xxyyzz_zzzz_0, ta1_y_xxyzz_xxxz_0, ta1_y_xxyzz_xxxz_1, ta1_y_xxyzz_xxzz_0, ta1_y_xxyzz_xxzz_1, ta1_y_xxyzz_xzzz_0, ta1_y_xxyzz_xzzz_1, ta1_y_xxzz_xxxz_0, ta1_y_xxzz_xxxz_1, ta1_y_xxzz_xxzz_0, ta1_y_xxzz_xxzz_1, ta1_y_xxzz_xzzz_0, ta1_y_xxzz_xzzz_1, ta1_y_xyyzz_xxyz_0, ta1_y_xyyzz_xxyz_1, ta1_y_xyyzz_xyyz_0, ta1_y_xyyzz_xyyz_1, ta1_y_xyyzz_xyz_0, ta1_y_xyyzz_xyz_1, ta1_y_xyyzz_xyzz_0, ta1_y_xyyzz_xyzz_1, ta1_y_xyyzz_yyyy_0, ta1_y_xyyzz_yyyy_1, ta1_y_xyyzz_yyyz_0, ta1_y_xyyzz_yyyz_1, ta1_y_xyyzz_yyz_0, ta1_y_xyyzz_yyz_1, ta1_y_xyyzz_yyzz_0, ta1_y_xyyzz_yyzz_1, ta1_y_xyyzz_yzz_0, ta1_y_xyyzz_yzz_1, ta1_y_xyyzz_yzzz_0, ta1_y_xyyzz_yzzz_1, ta1_y_xyyzz_zzzz_0, ta1_y_xyyzz_zzzz_1, ta1_y_yyzz_xxyz_0, ta1_y_yyzz_xxyz_1, ta1_y_yyzz_xyyz_0, ta1_y_yyzz_xyyz_1, ta1_y_yyzz_xyzz_0, ta1_y_yyzz_xyzz_1, ta1_y_yyzz_yyyy_0, ta1_y_yyzz_yyyy_1, ta1_y_yyzz_yyyz_0, ta1_y_yyzz_yyyz_1, ta1_y_yyzz_yyzz_0, ta1_y_yyzz_yyzz_1, ta1_y_yyzz_yzzz_0, ta1_y_yyzz_yzzz_1, ta1_y_yyzz_zzzz_0, ta1_y_yyzz_zzzz_1, ta_xxyzz_xxxz_1, ta_xxyzz_xxzz_1, ta_xxyzz_xzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyzz_xxxx_0[i] = ta1_y_xxyy_xxxx_0[i] * fe_0 - ta1_y_xxyy_xxxx_1[i] * fe_0 + ta1_y_xxyyz_xxxx_0[i] * pa_z[i] - ta1_y_xxyyz_xxxx_1[i] * pc_z[i];

        ta1_y_xxyyzz_xxxy_0[i] = ta1_y_xxyy_xxxy_0[i] * fe_0 - ta1_y_xxyy_xxxy_1[i] * fe_0 + ta1_y_xxyyz_xxxy_0[i] * pa_z[i] - ta1_y_xxyyz_xxxy_1[i] * pc_z[i];

        ta1_y_xxyyzz_xxxz_0[i] = ta1_y_xxzz_xxxz_0[i] * fe_0 - ta1_y_xxzz_xxxz_1[i] * fe_0 + ta_xxyzz_xxxz_1[i] + ta1_y_xxyzz_xxxz_0[i] * pa_y[i] - ta1_y_xxyzz_xxxz_1[i] * pc_y[i];

        ta1_y_xxyyzz_xxyy_0[i] = ta1_y_xxyy_xxyy_0[i] * fe_0 - ta1_y_xxyy_xxyy_1[i] * fe_0 + ta1_y_xxyyz_xxyy_0[i] * pa_z[i] - ta1_y_xxyyz_xxyy_1[i] * pc_z[i];

        ta1_y_xxyyzz_xxyz_0[i] = ta1_y_yyzz_xxyz_0[i] * fe_0 - ta1_y_yyzz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xyyzz_xyz_1[i] * fe_0 + ta1_y_xyyzz_xxyz_0[i] * pa_x[i] - ta1_y_xyyzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxyyzz_xxzz_0[i] = ta1_y_xxzz_xxzz_0[i] * fe_0 - ta1_y_xxzz_xxzz_1[i] * fe_0 + ta_xxyzz_xxzz_1[i] + ta1_y_xxyzz_xxzz_0[i] * pa_y[i] - ta1_y_xxyzz_xxzz_1[i] * pc_y[i];

        ta1_y_xxyyzz_xyyy_0[i] = ta1_y_xxyy_xyyy_0[i] * fe_0 - ta1_y_xxyy_xyyy_1[i] * fe_0 + ta1_y_xxyyz_xyyy_0[i] * pa_z[i] - ta1_y_xxyyz_xyyy_1[i] * pc_z[i];

        ta1_y_xxyyzz_xyyz_0[i] = ta1_y_yyzz_xyyz_0[i] * fe_0 - ta1_y_yyzz_xyyz_1[i] * fe_0 + ta1_y_xyyzz_yyz_0[i] * fe_0 - ta1_y_xyyzz_yyz_1[i] * fe_0 + ta1_y_xyyzz_xyyz_0[i] * pa_x[i] - ta1_y_xyyzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxyyzz_xyzz_0[i] = ta1_y_yyzz_xyzz_0[i] * fe_0 - ta1_y_yyzz_xyzz_1[i] * fe_0 + ta1_y_xyyzz_yzz_0[i] * fe_0 - ta1_y_xyyzz_yzz_1[i] * fe_0 + ta1_y_xyyzz_xyzz_0[i] * pa_x[i] - ta1_y_xyyzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxyyzz_xzzz_0[i] = ta1_y_xxzz_xzzz_0[i] * fe_0 - ta1_y_xxzz_xzzz_1[i] * fe_0 + ta_xxyzz_xzzz_1[i] + ta1_y_xxyzz_xzzz_0[i] * pa_y[i] - ta1_y_xxyzz_xzzz_1[i] * pc_y[i];

        ta1_y_xxyyzz_yyyy_0[i] = ta1_y_yyzz_yyyy_0[i] * fe_0 - ta1_y_yyzz_yyyy_1[i] * fe_0 + ta1_y_xyyzz_yyyy_0[i] * pa_x[i] - ta1_y_xyyzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxyyzz_yyyz_0[i] = ta1_y_yyzz_yyyz_0[i] * fe_0 - ta1_y_yyzz_yyyz_1[i] * fe_0 + ta1_y_xyyzz_yyyz_0[i] * pa_x[i] - ta1_y_xyyzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxyyzz_yyzz_0[i] = ta1_y_yyzz_yyzz_0[i] * fe_0 - ta1_y_yyzz_yyzz_1[i] * fe_0 + ta1_y_xyyzz_yyzz_0[i] * pa_x[i] - ta1_y_xyyzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxyyzz_yzzz_0[i] = ta1_y_yyzz_yzzz_0[i] * fe_0 - ta1_y_yyzz_yzzz_1[i] * fe_0 + ta1_y_xyyzz_yzzz_0[i] * pa_x[i] - ta1_y_xyyzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxyyzz_zzzz_0[i] = ta1_y_yyzz_zzzz_0[i] * fe_0 - ta1_y_yyzz_zzzz_1[i] * fe_0 + ta1_y_xyyzz_zzzz_0[i] * pa_x[i] - ta1_y_xyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 615-630 components of targeted buffer : IG

    auto ta1_y_xxyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 615);

    auto ta1_y_xxyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 616);

    auto ta1_y_xxyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 617);

    auto ta1_y_xxyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 618);

    auto ta1_y_xxyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 619);

    auto ta1_y_xxyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 620);

    auto ta1_y_xxyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 621);

    auto ta1_y_xxyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 622);

    auto ta1_y_xxyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 623);

    auto ta1_y_xxyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 624);

    auto ta1_y_xxyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 625);

    auto ta1_y_xxyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 626);

    auto ta1_y_xxyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 627);

    auto ta1_y_xxyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 628);

    auto ta1_y_xxyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 629);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxyz_xxxy_0, ta1_y_xxyz_xxxy_1, ta1_y_xxyz_xxyy_0, ta1_y_xxyz_xxyy_1, ta1_y_xxyz_xyyy_0, ta1_y_xxyz_xyyy_1, ta1_y_xxyzz_xxxy_0, ta1_y_xxyzz_xxxy_1, ta1_y_xxyzz_xxyy_0, ta1_y_xxyzz_xxyy_1, ta1_y_xxyzz_xyyy_0, ta1_y_xxyzz_xyyy_1, ta1_y_xxyzzz_xxxx_0, ta1_y_xxyzzz_xxxy_0, ta1_y_xxyzzz_xxxz_0, ta1_y_xxyzzz_xxyy_0, ta1_y_xxyzzz_xxyz_0, ta1_y_xxyzzz_xxzz_0, ta1_y_xxyzzz_xyyy_0, ta1_y_xxyzzz_xyyz_0, ta1_y_xxyzzz_xyzz_0, ta1_y_xxyzzz_xzzz_0, ta1_y_xxyzzz_yyyy_0, ta1_y_xxyzzz_yyyz_0, ta1_y_xxyzzz_yyzz_0, ta1_y_xxyzzz_yzzz_0, ta1_y_xxyzzz_zzzz_0, ta1_y_xxzzz_xxxx_0, ta1_y_xxzzz_xxxx_1, ta1_y_xxzzz_xxxz_0, ta1_y_xxzzz_xxxz_1, ta1_y_xxzzz_xxyz_0, ta1_y_xxzzz_xxyz_1, ta1_y_xxzzz_xxz_0, ta1_y_xxzzz_xxz_1, ta1_y_xxzzz_xxzz_0, ta1_y_xxzzz_xxzz_1, ta1_y_xxzzz_xyyz_0, ta1_y_xxzzz_xyyz_1, ta1_y_xxzzz_xyz_0, ta1_y_xxzzz_xyz_1, ta1_y_xxzzz_xyzz_0, ta1_y_xxzzz_xyzz_1, ta1_y_xxzzz_xzz_0, ta1_y_xxzzz_xzz_1, ta1_y_xxzzz_xzzz_0, ta1_y_xxzzz_xzzz_1, ta1_y_xxzzz_zzzz_0, ta1_y_xxzzz_zzzz_1, ta1_y_xyzzz_yyyy_0, ta1_y_xyzzz_yyyy_1, ta1_y_xyzzz_yyyz_0, ta1_y_xyzzz_yyyz_1, ta1_y_xyzzz_yyzz_0, ta1_y_xyzzz_yyzz_1, ta1_y_xyzzz_yzzz_0, ta1_y_xyzzz_yzzz_1, ta1_y_yzzz_yyyy_0, ta1_y_yzzz_yyyy_1, ta1_y_yzzz_yyyz_0, ta1_y_yzzz_yyyz_1, ta1_y_yzzz_yyzz_0, ta1_y_yzzz_yyzz_1, ta1_y_yzzz_yzzz_0, ta1_y_yzzz_yzzz_1, ta_xxzzz_xxxx_1, ta_xxzzz_xxxz_1, ta_xxzzz_xxyz_1, ta_xxzzz_xxzz_1, ta_xxzzz_xyyz_1, ta_xxzzz_xyzz_1, ta_xxzzz_xzzz_1, ta_xxzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzzz_xxxx_0[i] = ta_xxzzz_xxxx_1[i] + ta1_y_xxzzz_xxxx_0[i] * pa_y[i] - ta1_y_xxzzz_xxxx_1[i] * pc_y[i];

        ta1_y_xxyzzz_xxxy_0[i] = 2.0 * ta1_y_xxyz_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xxyz_xxxy_1[i] * fe_0 + ta1_y_xxyzz_xxxy_0[i] * pa_z[i] - ta1_y_xxyzz_xxxy_1[i] * pc_z[i];

        ta1_y_xxyzzz_xxxz_0[i] = ta_xxzzz_xxxz_1[i] + ta1_y_xxzzz_xxxz_0[i] * pa_y[i] - ta1_y_xxzzz_xxxz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xxyy_0[i] = 2.0 * ta1_y_xxyz_xxyy_0[i] * fe_0 - 2.0 * ta1_y_xxyz_xxyy_1[i] * fe_0 + ta1_y_xxyzz_xxyy_0[i] * pa_z[i] - ta1_y_xxyzz_xxyy_1[i] * pc_z[i];

        ta1_y_xxyzzz_xxyz_0[i] = ta1_y_xxzzz_xxz_0[i] * fe_0 - ta1_y_xxzzz_xxz_1[i] * fe_0 + ta_xxzzz_xxyz_1[i] + ta1_y_xxzzz_xxyz_0[i] * pa_y[i] - ta1_y_xxzzz_xxyz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xxzz_0[i] = ta_xxzzz_xxzz_1[i] + ta1_y_xxzzz_xxzz_0[i] * pa_y[i] - ta1_y_xxzzz_xxzz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xyyy_0[i] = 2.0 * ta1_y_xxyz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_xxyz_xyyy_1[i] * fe_0 + ta1_y_xxyzz_xyyy_0[i] * pa_z[i] - ta1_y_xxyzz_xyyy_1[i] * pc_z[i];

        ta1_y_xxyzzz_xyyz_0[i] = 2.0 * ta1_y_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xxzzz_xyz_1[i] * fe_0 + ta_xxzzz_xyyz_1[i] + ta1_y_xxzzz_xyyz_0[i] * pa_y[i] - ta1_y_xxzzz_xyyz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xyzz_0[i] = ta1_y_xxzzz_xzz_0[i] * fe_0 - ta1_y_xxzzz_xzz_1[i] * fe_0 + ta_xxzzz_xyzz_1[i] + ta1_y_xxzzz_xyzz_0[i] * pa_y[i] - ta1_y_xxzzz_xyzz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xzzz_0[i] = ta_xxzzz_xzzz_1[i] + ta1_y_xxzzz_xzzz_0[i] * pa_y[i] - ta1_y_xxzzz_xzzz_1[i] * pc_y[i];

        ta1_y_xxyzzz_yyyy_0[i] = ta1_y_yzzz_yyyy_0[i] * fe_0 - ta1_y_yzzz_yyyy_1[i] * fe_0 + ta1_y_xyzzz_yyyy_0[i] * pa_x[i] - ta1_y_xyzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxyzzz_yyyz_0[i] = ta1_y_yzzz_yyyz_0[i] * fe_0 - ta1_y_yzzz_yyyz_1[i] * fe_0 + ta1_y_xyzzz_yyyz_0[i] * pa_x[i] - ta1_y_xyzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxyzzz_yyzz_0[i] = ta1_y_yzzz_yyzz_0[i] * fe_0 - ta1_y_yzzz_yyzz_1[i] * fe_0 + ta1_y_xyzzz_yyzz_0[i] * pa_x[i] - ta1_y_xyzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxyzzz_yzzz_0[i] = ta1_y_yzzz_yzzz_0[i] * fe_0 - ta1_y_yzzz_yzzz_1[i] * fe_0 + ta1_y_xyzzz_yzzz_0[i] * pa_x[i] - ta1_y_xyzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxyzzz_zzzz_0[i] = ta_xxzzz_zzzz_1[i] + ta1_y_xxzzz_zzzz_0[i] * pa_y[i] - ta1_y_xxzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 630-645 components of targeted buffer : IG

    auto ta1_y_xxzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 630);

    auto ta1_y_xxzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 631);

    auto ta1_y_xxzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 632);

    auto ta1_y_xxzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 633);

    auto ta1_y_xxzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 634);

    auto ta1_y_xxzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 635);

    auto ta1_y_xxzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 636);

    auto ta1_y_xxzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 637);

    auto ta1_y_xxzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 638);

    auto ta1_y_xxzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 639);

    auto ta1_y_xxzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 640);

    auto ta1_y_xxzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 641);

    auto ta1_y_xxzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 642);

    auto ta1_y_xxzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 643);

    auto ta1_y_xxzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 644);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxzz_xxxx_0, ta1_y_xxzz_xxxx_1, ta1_y_xxzz_xxxy_0, ta1_y_xxzz_xxxy_1, ta1_y_xxzz_xxyy_0, ta1_y_xxzz_xxyy_1, ta1_y_xxzz_xyyy_0, ta1_y_xxzz_xyyy_1, ta1_y_xxzzz_xxxx_0, ta1_y_xxzzz_xxxx_1, ta1_y_xxzzz_xxxy_0, ta1_y_xxzzz_xxxy_1, ta1_y_xxzzz_xxyy_0, ta1_y_xxzzz_xxyy_1, ta1_y_xxzzz_xyyy_0, ta1_y_xxzzz_xyyy_1, ta1_y_xxzzzz_xxxx_0, ta1_y_xxzzzz_xxxy_0, ta1_y_xxzzzz_xxxz_0, ta1_y_xxzzzz_xxyy_0, ta1_y_xxzzzz_xxyz_0, ta1_y_xxzzzz_xxzz_0, ta1_y_xxzzzz_xyyy_0, ta1_y_xxzzzz_xyyz_0, ta1_y_xxzzzz_xyzz_0, ta1_y_xxzzzz_xzzz_0, ta1_y_xxzzzz_yyyy_0, ta1_y_xxzzzz_yyyz_0, ta1_y_xxzzzz_yyzz_0, ta1_y_xxzzzz_yzzz_0, ta1_y_xxzzzz_zzzz_0, ta1_y_xzzzz_xxxz_0, ta1_y_xzzzz_xxxz_1, ta1_y_xzzzz_xxyz_0, ta1_y_xzzzz_xxyz_1, ta1_y_xzzzz_xxz_0, ta1_y_xzzzz_xxz_1, ta1_y_xzzzz_xxzz_0, ta1_y_xzzzz_xxzz_1, ta1_y_xzzzz_xyyz_0, ta1_y_xzzzz_xyyz_1, ta1_y_xzzzz_xyz_0, ta1_y_xzzzz_xyz_1, ta1_y_xzzzz_xyzz_0, ta1_y_xzzzz_xyzz_1, ta1_y_xzzzz_xzz_0, ta1_y_xzzzz_xzz_1, ta1_y_xzzzz_xzzz_0, ta1_y_xzzzz_xzzz_1, ta1_y_xzzzz_yyyy_0, ta1_y_xzzzz_yyyy_1, ta1_y_xzzzz_yyyz_0, ta1_y_xzzzz_yyyz_1, ta1_y_xzzzz_yyz_0, ta1_y_xzzzz_yyz_1, ta1_y_xzzzz_yyzz_0, ta1_y_xzzzz_yyzz_1, ta1_y_xzzzz_yzz_0, ta1_y_xzzzz_yzz_1, ta1_y_xzzzz_yzzz_0, ta1_y_xzzzz_yzzz_1, ta1_y_xzzzz_zzz_0, ta1_y_xzzzz_zzz_1, ta1_y_xzzzz_zzzz_0, ta1_y_xzzzz_zzzz_1, ta1_y_zzzz_xxxz_0, ta1_y_zzzz_xxxz_1, ta1_y_zzzz_xxyz_0, ta1_y_zzzz_xxyz_1, ta1_y_zzzz_xxzz_0, ta1_y_zzzz_xxzz_1, ta1_y_zzzz_xyyz_0, ta1_y_zzzz_xyyz_1, ta1_y_zzzz_xyzz_0, ta1_y_zzzz_xyzz_1, ta1_y_zzzz_xzzz_0, ta1_y_zzzz_xzzz_1, ta1_y_zzzz_yyyy_0, ta1_y_zzzz_yyyy_1, ta1_y_zzzz_yyyz_0, ta1_y_zzzz_yyyz_1, ta1_y_zzzz_yyzz_0, ta1_y_zzzz_yyzz_1, ta1_y_zzzz_yzzz_0, ta1_y_zzzz_yzzz_1, ta1_y_zzzz_zzzz_0, ta1_y_zzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzzz_xxxx_0[i] = 3.0 * ta1_y_xxzz_xxxx_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxxx_1[i] * fe_0 + ta1_y_xxzzz_xxxx_0[i] * pa_z[i] - ta1_y_xxzzz_xxxx_1[i] * pc_z[i];

        ta1_y_xxzzzz_xxxy_0[i] = 3.0 * ta1_y_xxzz_xxxy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxxy_1[i] * fe_0 + ta1_y_xxzzz_xxxy_0[i] * pa_z[i] - ta1_y_xxzzz_xxxy_1[i] * pc_z[i];

        ta1_y_xxzzzz_xxxz_0[i] = ta1_y_zzzz_xxxz_0[i] * fe_0 - ta1_y_zzzz_xxxz_1[i] * fe_0 + 3.0 * ta1_y_xzzzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_xzzzz_xxz_1[i] * fe_0 + ta1_y_xzzzz_xxxz_0[i] * pa_x[i] - ta1_y_xzzzz_xxxz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xxyy_0[i] = 3.0 * ta1_y_xxzz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxyy_1[i] * fe_0 + ta1_y_xxzzz_xxyy_0[i] * pa_z[i] - ta1_y_xxzzz_xxyy_1[i] * pc_z[i];

        ta1_y_xxzzzz_xxyz_0[i] = ta1_y_zzzz_xxyz_0[i] * fe_0 - ta1_y_zzzz_xxyz_1[i] * fe_0 + 2.0 * ta1_y_xzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xzzzz_xyz_1[i] * fe_0 + ta1_y_xzzzz_xxyz_0[i] * pa_x[i] - ta1_y_xzzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xxzz_0[i] = ta1_y_zzzz_xxzz_0[i] * fe_0 - ta1_y_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xzzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_xzzzz_xzz_1[i] * fe_0 + ta1_y_xzzzz_xxzz_0[i] * pa_x[i] - ta1_y_xzzzz_xxzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xyyy_0[i] = 3.0 * ta1_y_xxzz_xyyy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xyyy_1[i] * fe_0 + ta1_y_xxzzz_xyyy_0[i] * pa_z[i] - ta1_y_xxzzz_xyyy_1[i] * pc_z[i];

        ta1_y_xxzzzz_xyyz_0[i] = ta1_y_zzzz_xyyz_0[i] * fe_0 - ta1_y_zzzz_xyyz_1[i] * fe_0 + ta1_y_xzzzz_yyz_0[i] * fe_0 - ta1_y_xzzzz_yyz_1[i] * fe_0 + ta1_y_xzzzz_xyyz_0[i] * pa_x[i] - ta1_y_xzzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xyzz_0[i] = ta1_y_zzzz_xyzz_0[i] * fe_0 - ta1_y_zzzz_xyzz_1[i] * fe_0 + ta1_y_xzzzz_yzz_0[i] * fe_0 - ta1_y_xzzzz_yzz_1[i] * fe_0 + ta1_y_xzzzz_xyzz_0[i] * pa_x[i] - ta1_y_xzzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xzzz_0[i] = ta1_y_zzzz_xzzz_0[i] * fe_0 - ta1_y_zzzz_xzzz_1[i] * fe_0 + ta1_y_xzzzz_zzz_0[i] * fe_0 - ta1_y_xzzzz_zzz_1[i] * fe_0 + ta1_y_xzzzz_xzzz_0[i] * pa_x[i] - ta1_y_xzzzz_xzzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_yyyy_0[i] = ta1_y_zzzz_yyyy_0[i] * fe_0 - ta1_y_zzzz_yyyy_1[i] * fe_0 + ta1_y_xzzzz_yyyy_0[i] * pa_x[i] - ta1_y_xzzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xxzzzz_yyyz_0[i] = ta1_y_zzzz_yyyz_0[i] * fe_0 - ta1_y_zzzz_yyyz_1[i] * fe_0 + ta1_y_xzzzz_yyyz_0[i] * pa_x[i] - ta1_y_xzzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xxzzzz_yyzz_0[i] = ta1_y_zzzz_yyzz_0[i] * fe_0 - ta1_y_zzzz_yyzz_1[i] * fe_0 + ta1_y_xzzzz_yyzz_0[i] * pa_x[i] - ta1_y_xzzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_yzzz_0[i] = ta1_y_zzzz_yzzz_0[i] * fe_0 - ta1_y_zzzz_yzzz_1[i] * fe_0 + ta1_y_xzzzz_yzzz_0[i] * pa_x[i] - ta1_y_xzzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_zzzz_0[i] = ta1_y_zzzz_zzzz_0[i] * fe_0 - ta1_y_zzzz_zzzz_1[i] * fe_0 + ta1_y_xzzzz_zzzz_0[i] * pa_x[i] - ta1_y_xzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 645-660 components of targeted buffer : IG

    auto ta1_y_xyyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 645);

    auto ta1_y_xyyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 646);

    auto ta1_y_xyyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 647);

    auto ta1_y_xyyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 648);

    auto ta1_y_xyyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 649);

    auto ta1_y_xyyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 650);

    auto ta1_y_xyyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 651);

    auto ta1_y_xyyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 652);

    auto ta1_y_xyyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 653);

    auto ta1_y_xyyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 654);

    auto ta1_y_xyyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 655);

    auto ta1_y_xyyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 656);

    auto ta1_y_xyyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 657);

    auto ta1_y_xyyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 658);

    auto ta1_y_xyyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 659);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyyyy_xxxx_0, ta1_y_xyyyyy_xxxy_0, ta1_y_xyyyyy_xxxz_0, ta1_y_xyyyyy_xxyy_0, ta1_y_xyyyyy_xxyz_0, ta1_y_xyyyyy_xxzz_0, ta1_y_xyyyyy_xyyy_0, ta1_y_xyyyyy_xyyz_0, ta1_y_xyyyyy_xyzz_0, ta1_y_xyyyyy_xzzz_0, ta1_y_xyyyyy_yyyy_0, ta1_y_xyyyyy_yyyz_0, ta1_y_xyyyyy_yyzz_0, ta1_y_xyyyyy_yzzz_0, ta1_y_xyyyyy_zzzz_0, ta1_y_yyyyy_xxx_0, ta1_y_yyyyy_xxx_1, ta1_y_yyyyy_xxxx_0, ta1_y_yyyyy_xxxx_1, ta1_y_yyyyy_xxxy_0, ta1_y_yyyyy_xxxy_1, ta1_y_yyyyy_xxxz_0, ta1_y_yyyyy_xxxz_1, ta1_y_yyyyy_xxy_0, ta1_y_yyyyy_xxy_1, ta1_y_yyyyy_xxyy_0, ta1_y_yyyyy_xxyy_1, ta1_y_yyyyy_xxyz_0, ta1_y_yyyyy_xxyz_1, ta1_y_yyyyy_xxz_0, ta1_y_yyyyy_xxz_1, ta1_y_yyyyy_xxzz_0, ta1_y_yyyyy_xxzz_1, ta1_y_yyyyy_xyy_0, ta1_y_yyyyy_xyy_1, ta1_y_yyyyy_xyyy_0, ta1_y_yyyyy_xyyy_1, ta1_y_yyyyy_xyyz_0, ta1_y_yyyyy_xyyz_1, ta1_y_yyyyy_xyz_0, ta1_y_yyyyy_xyz_1, ta1_y_yyyyy_xyzz_0, ta1_y_yyyyy_xyzz_1, ta1_y_yyyyy_xzz_0, ta1_y_yyyyy_xzz_1, ta1_y_yyyyy_xzzz_0, ta1_y_yyyyy_xzzz_1, ta1_y_yyyyy_yyy_0, ta1_y_yyyyy_yyy_1, ta1_y_yyyyy_yyyy_0, ta1_y_yyyyy_yyyy_1, ta1_y_yyyyy_yyyz_0, ta1_y_yyyyy_yyyz_1, ta1_y_yyyyy_yyz_0, ta1_y_yyyyy_yyz_1, ta1_y_yyyyy_yyzz_0, ta1_y_yyyyy_yyzz_1, ta1_y_yyyyy_yzz_0, ta1_y_yyyyy_yzz_1, ta1_y_yyyyy_yzzz_0, ta1_y_yyyyy_yzzz_1, ta1_y_yyyyy_zzz_0, ta1_y_yyyyy_zzz_1, ta1_y_yyyyy_zzzz_0, ta1_y_yyyyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyy_xxxx_0[i] = 4.0 * ta1_y_yyyyy_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyyyy_xxx_1[i] * fe_0 + ta1_y_yyyyy_xxxx_0[i] * pa_x[i] - ta1_y_yyyyy_xxxx_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxxy_0[i] = 3.0 * ta1_y_yyyyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_xxy_1[i] * fe_0 + ta1_y_yyyyy_xxxy_0[i] * pa_x[i] - ta1_y_yyyyy_xxxy_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxxz_0[i] = 3.0 * ta1_y_yyyyy_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_xxz_1[i] * fe_0 + ta1_y_yyyyy_xxxz_0[i] * pa_x[i] - ta1_y_yyyyy_xxxz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxyy_0[i] = 2.0 * ta1_y_yyyyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xyy_1[i] * fe_0 + ta1_y_yyyyy_xxyy_0[i] * pa_x[i] - ta1_y_yyyyy_xxyy_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxyz_0[i] = 2.0 * ta1_y_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xyz_1[i] * fe_0 + ta1_y_yyyyy_xxyz_0[i] * pa_x[i] - ta1_y_yyyyy_xxyz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxzz_0[i] = 2.0 * ta1_y_yyyyy_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xzz_1[i] * fe_0 + ta1_y_yyyyy_xxzz_0[i] * pa_x[i] - ta1_y_yyyyy_xxzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xyyy_0[i] = ta1_y_yyyyy_yyy_0[i] * fe_0 - ta1_y_yyyyy_yyy_1[i] * fe_0 + ta1_y_yyyyy_xyyy_0[i] * pa_x[i] - ta1_y_yyyyy_xyyy_1[i] * pc_x[i];

        ta1_y_xyyyyy_xyyz_0[i] = ta1_y_yyyyy_yyz_0[i] * fe_0 - ta1_y_yyyyy_yyz_1[i] * fe_0 + ta1_y_yyyyy_xyyz_0[i] * pa_x[i] - ta1_y_yyyyy_xyyz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xyzz_0[i] = ta1_y_yyyyy_yzz_0[i] * fe_0 - ta1_y_yyyyy_yzz_1[i] * fe_0 + ta1_y_yyyyy_xyzz_0[i] * pa_x[i] - ta1_y_yyyyy_xyzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xzzz_0[i] = ta1_y_yyyyy_zzz_0[i] * fe_0 - ta1_y_yyyyy_zzz_1[i] * fe_0 + ta1_y_yyyyy_xzzz_0[i] * pa_x[i] - ta1_y_yyyyy_xzzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_yyyy_0[i] = ta1_y_yyyyy_yyyy_0[i] * pa_x[i] - ta1_y_yyyyy_yyyy_1[i] * pc_x[i];

        ta1_y_xyyyyy_yyyz_0[i] = ta1_y_yyyyy_yyyz_0[i] * pa_x[i] - ta1_y_yyyyy_yyyz_1[i] * pc_x[i];

        ta1_y_xyyyyy_yyzz_0[i] = ta1_y_yyyyy_yyzz_0[i] * pa_x[i] - ta1_y_yyyyy_yyzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_yzzz_0[i] = ta1_y_yyyyy_yzzz_0[i] * pa_x[i] - ta1_y_yyyyy_yzzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_zzzz_0[i] = ta1_y_yyyyy_zzzz_0[i] * pa_x[i] - ta1_y_yyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 660-675 components of targeted buffer : IG

    auto ta1_y_xyyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 660);

    auto ta1_y_xyyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 661);

    auto ta1_y_xyyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 662);

    auto ta1_y_xyyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 663);

    auto ta1_y_xyyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 664);

    auto ta1_y_xyyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 665);

    auto ta1_y_xyyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 666);

    auto ta1_y_xyyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 667);

    auto ta1_y_xyyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 668);

    auto ta1_y_xyyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 669);

    auto ta1_y_xyyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 670);

    auto ta1_y_xyyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 671);

    auto ta1_y_xyyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 672);

    auto ta1_y_xyyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 673);

    auto ta1_y_xyyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 674);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xyyyy_xxxx_0, ta1_y_xyyyy_xxxx_1, ta1_y_xyyyy_xxxy_0, ta1_y_xyyyy_xxxy_1, ta1_y_xyyyy_xxyy_0, ta1_y_xyyyy_xxyy_1, ta1_y_xyyyy_xyyy_0, ta1_y_xyyyy_xyyy_1, ta1_y_xyyyyz_xxxx_0, ta1_y_xyyyyz_xxxy_0, ta1_y_xyyyyz_xxxz_0, ta1_y_xyyyyz_xxyy_0, ta1_y_xyyyyz_xxyz_0, ta1_y_xyyyyz_xxzz_0, ta1_y_xyyyyz_xyyy_0, ta1_y_xyyyyz_xyyz_0, ta1_y_xyyyyz_xyzz_0, ta1_y_xyyyyz_xzzz_0, ta1_y_xyyyyz_yyyy_0, ta1_y_xyyyyz_yyyz_0, ta1_y_xyyyyz_yyzz_0, ta1_y_xyyyyz_yzzz_0, ta1_y_xyyyyz_zzzz_0, ta1_y_yyyyz_xxxz_0, ta1_y_yyyyz_xxxz_1, ta1_y_yyyyz_xxyz_0, ta1_y_yyyyz_xxyz_1, ta1_y_yyyyz_xxz_0, ta1_y_yyyyz_xxz_1, ta1_y_yyyyz_xxzz_0, ta1_y_yyyyz_xxzz_1, ta1_y_yyyyz_xyyz_0, ta1_y_yyyyz_xyyz_1, ta1_y_yyyyz_xyz_0, ta1_y_yyyyz_xyz_1, ta1_y_yyyyz_xyzz_0, ta1_y_yyyyz_xyzz_1, ta1_y_yyyyz_xzz_0, ta1_y_yyyyz_xzz_1, ta1_y_yyyyz_xzzz_0, ta1_y_yyyyz_xzzz_1, ta1_y_yyyyz_yyyy_0, ta1_y_yyyyz_yyyy_1, ta1_y_yyyyz_yyyz_0, ta1_y_yyyyz_yyyz_1, ta1_y_yyyyz_yyz_0, ta1_y_yyyyz_yyz_1, ta1_y_yyyyz_yyzz_0, ta1_y_yyyyz_yyzz_1, ta1_y_yyyyz_yzz_0, ta1_y_yyyyz_yzz_1, ta1_y_yyyyz_yzzz_0, ta1_y_yyyyz_yzzz_1, ta1_y_yyyyz_zzz_0, ta1_y_yyyyz_zzz_1, ta1_y_yyyyz_zzzz_0, ta1_y_yyyyz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyz_xxxx_0[i] = ta1_y_xyyyy_xxxx_0[i] * pa_z[i] - ta1_y_xyyyy_xxxx_1[i] * pc_z[i];

        ta1_y_xyyyyz_xxxy_0[i] = ta1_y_xyyyy_xxxy_0[i] * pa_z[i] - ta1_y_xyyyy_xxxy_1[i] * pc_z[i];

        ta1_y_xyyyyz_xxxz_0[i] = 3.0 * ta1_y_yyyyz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyyyz_xxz_1[i] * fe_0 + ta1_y_yyyyz_xxxz_0[i] * pa_x[i] - ta1_y_yyyyz_xxxz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xxyy_0[i] = ta1_y_xyyyy_xxyy_0[i] * pa_z[i] - ta1_y_xyyyy_xxyy_1[i] * pc_z[i];

        ta1_y_xyyyyz_xxyz_0[i] = 2.0 * ta1_y_yyyyz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyz_xyz_1[i] * fe_0 + ta1_y_yyyyz_xxyz_0[i] * pa_x[i] - ta1_y_yyyyz_xxyz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xxzz_0[i] = 2.0 * ta1_y_yyyyz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyyyz_xzz_1[i] * fe_0 + ta1_y_yyyyz_xxzz_0[i] * pa_x[i] - ta1_y_yyyyz_xxzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xyyy_0[i] = ta1_y_xyyyy_xyyy_0[i] * pa_z[i] - ta1_y_xyyyy_xyyy_1[i] * pc_z[i];

        ta1_y_xyyyyz_xyyz_0[i] = ta1_y_yyyyz_yyz_0[i] * fe_0 - ta1_y_yyyyz_yyz_1[i] * fe_0 + ta1_y_yyyyz_xyyz_0[i] * pa_x[i] - ta1_y_yyyyz_xyyz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xyzz_0[i] = ta1_y_yyyyz_yzz_0[i] * fe_0 - ta1_y_yyyyz_yzz_1[i] * fe_0 + ta1_y_yyyyz_xyzz_0[i] * pa_x[i] - ta1_y_yyyyz_xyzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xzzz_0[i] = ta1_y_yyyyz_zzz_0[i] * fe_0 - ta1_y_yyyyz_zzz_1[i] * fe_0 + ta1_y_yyyyz_xzzz_0[i] * pa_x[i] - ta1_y_yyyyz_xzzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_yyyy_0[i] = ta1_y_yyyyz_yyyy_0[i] * pa_x[i] - ta1_y_yyyyz_yyyy_1[i] * pc_x[i];

        ta1_y_xyyyyz_yyyz_0[i] = ta1_y_yyyyz_yyyz_0[i] * pa_x[i] - ta1_y_yyyyz_yyyz_1[i] * pc_x[i];

        ta1_y_xyyyyz_yyzz_0[i] = ta1_y_yyyyz_yyzz_0[i] * pa_x[i] - ta1_y_yyyyz_yyzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_yzzz_0[i] = ta1_y_yyyyz_yzzz_0[i] * pa_x[i] - ta1_y_yyyyz_yzzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_zzzz_0[i] = ta1_y_yyyyz_zzzz_0[i] * pa_x[i] - ta1_y_yyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 675-690 components of targeted buffer : IG

    auto ta1_y_xyyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 675);

    auto ta1_y_xyyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 676);

    auto ta1_y_xyyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 677);

    auto ta1_y_xyyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 678);

    auto ta1_y_xyyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 679);

    auto ta1_y_xyyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 680);

    auto ta1_y_xyyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 681);

    auto ta1_y_xyyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 682);

    auto ta1_y_xyyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 683);

    auto ta1_y_xyyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 684);

    auto ta1_y_xyyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 685);

    auto ta1_y_xyyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 686);

    auto ta1_y_xyyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 687);

    auto ta1_y_xyyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 688);

    auto ta1_y_xyyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 689);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyyzz_xxxx_0, ta1_y_xyyyzz_xxxy_0, ta1_y_xyyyzz_xxxz_0, ta1_y_xyyyzz_xxyy_0, ta1_y_xyyyzz_xxyz_0, ta1_y_xyyyzz_xxzz_0, ta1_y_xyyyzz_xyyy_0, ta1_y_xyyyzz_xyyz_0, ta1_y_xyyyzz_xyzz_0, ta1_y_xyyyzz_xzzz_0, ta1_y_xyyyzz_yyyy_0, ta1_y_xyyyzz_yyyz_0, ta1_y_xyyyzz_yyzz_0, ta1_y_xyyyzz_yzzz_0, ta1_y_xyyyzz_zzzz_0, ta1_y_yyyzz_xxx_0, ta1_y_yyyzz_xxx_1, ta1_y_yyyzz_xxxx_0, ta1_y_yyyzz_xxxx_1, ta1_y_yyyzz_xxxy_0, ta1_y_yyyzz_xxxy_1, ta1_y_yyyzz_xxxz_0, ta1_y_yyyzz_xxxz_1, ta1_y_yyyzz_xxy_0, ta1_y_yyyzz_xxy_1, ta1_y_yyyzz_xxyy_0, ta1_y_yyyzz_xxyy_1, ta1_y_yyyzz_xxyz_0, ta1_y_yyyzz_xxyz_1, ta1_y_yyyzz_xxz_0, ta1_y_yyyzz_xxz_1, ta1_y_yyyzz_xxzz_0, ta1_y_yyyzz_xxzz_1, ta1_y_yyyzz_xyy_0, ta1_y_yyyzz_xyy_1, ta1_y_yyyzz_xyyy_0, ta1_y_yyyzz_xyyy_1, ta1_y_yyyzz_xyyz_0, ta1_y_yyyzz_xyyz_1, ta1_y_yyyzz_xyz_0, ta1_y_yyyzz_xyz_1, ta1_y_yyyzz_xyzz_0, ta1_y_yyyzz_xyzz_1, ta1_y_yyyzz_xzz_0, ta1_y_yyyzz_xzz_1, ta1_y_yyyzz_xzzz_0, ta1_y_yyyzz_xzzz_1, ta1_y_yyyzz_yyy_0, ta1_y_yyyzz_yyy_1, ta1_y_yyyzz_yyyy_0, ta1_y_yyyzz_yyyy_1, ta1_y_yyyzz_yyyz_0, ta1_y_yyyzz_yyyz_1, ta1_y_yyyzz_yyz_0, ta1_y_yyyzz_yyz_1, ta1_y_yyyzz_yyzz_0, ta1_y_yyyzz_yyzz_1, ta1_y_yyyzz_yzz_0, ta1_y_yyyzz_yzz_1, ta1_y_yyyzz_yzzz_0, ta1_y_yyyzz_yzzz_1, ta1_y_yyyzz_zzz_0, ta1_y_yyyzz_zzz_1, ta1_y_yyyzz_zzzz_0, ta1_y_yyyzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyzz_xxxx_0[i] = 4.0 * ta1_y_yyyzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyyzz_xxx_1[i] * fe_0 + ta1_y_yyyzz_xxxx_0[i] * pa_x[i] - ta1_y_yyyzz_xxxx_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxxy_0[i] = 3.0 * ta1_y_yyyzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyyzz_xxy_1[i] * fe_0 + ta1_y_yyyzz_xxxy_0[i] * pa_x[i] - ta1_y_yyyzz_xxxy_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxxz_0[i] = 3.0 * ta1_y_yyyzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyyzz_xxz_1[i] * fe_0 + ta1_y_yyyzz_xxxz_0[i] * pa_x[i] - ta1_y_yyyzz_xxxz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxyy_0[i] = 2.0 * ta1_y_yyyzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_xyy_1[i] * fe_0 + ta1_y_yyyzz_xxyy_0[i] * pa_x[i] - ta1_y_yyyzz_xxyy_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxyz_0[i] = 2.0 * ta1_y_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_xyz_1[i] * fe_0 + ta1_y_yyyzz_xxyz_0[i] * pa_x[i] - ta1_y_yyyzz_xxyz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxzz_0[i] = 2.0 * ta1_y_yyyzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_xzz_1[i] * fe_0 + ta1_y_yyyzz_xxzz_0[i] * pa_x[i] - ta1_y_yyyzz_xxzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xyyy_0[i] = ta1_y_yyyzz_yyy_0[i] * fe_0 - ta1_y_yyyzz_yyy_1[i] * fe_0 + ta1_y_yyyzz_xyyy_0[i] * pa_x[i] - ta1_y_yyyzz_xyyy_1[i] * pc_x[i];

        ta1_y_xyyyzz_xyyz_0[i] = ta1_y_yyyzz_yyz_0[i] * fe_0 - ta1_y_yyyzz_yyz_1[i] * fe_0 + ta1_y_yyyzz_xyyz_0[i] * pa_x[i] - ta1_y_yyyzz_xyyz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xyzz_0[i] = ta1_y_yyyzz_yzz_0[i] * fe_0 - ta1_y_yyyzz_yzz_1[i] * fe_0 + ta1_y_yyyzz_xyzz_0[i] * pa_x[i] - ta1_y_yyyzz_xyzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xzzz_0[i] = ta1_y_yyyzz_zzz_0[i] * fe_0 - ta1_y_yyyzz_zzz_1[i] * fe_0 + ta1_y_yyyzz_xzzz_0[i] * pa_x[i] - ta1_y_yyyzz_xzzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_yyyy_0[i] = ta1_y_yyyzz_yyyy_0[i] * pa_x[i] - ta1_y_yyyzz_yyyy_1[i] * pc_x[i];

        ta1_y_xyyyzz_yyyz_0[i] = ta1_y_yyyzz_yyyz_0[i] * pa_x[i] - ta1_y_yyyzz_yyyz_1[i] * pc_x[i];

        ta1_y_xyyyzz_yyzz_0[i] = ta1_y_yyyzz_yyzz_0[i] * pa_x[i] - ta1_y_yyyzz_yyzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_yzzz_0[i] = ta1_y_yyyzz_yzzz_0[i] * pa_x[i] - ta1_y_yyyzz_yzzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_zzzz_0[i] = ta1_y_yyyzz_zzzz_0[i] * pa_x[i] - ta1_y_yyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 690-705 components of targeted buffer : IG

    auto ta1_y_xyyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 690);

    auto ta1_y_xyyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 691);

    auto ta1_y_xyyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 692);

    auto ta1_y_xyyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 693);

    auto ta1_y_xyyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 694);

    auto ta1_y_xyyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 695);

    auto ta1_y_xyyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 696);

    auto ta1_y_xyyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 697);

    auto ta1_y_xyyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 698);

    auto ta1_y_xyyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 699);

    auto ta1_y_xyyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 700);

    auto ta1_y_xyyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 701);

    auto ta1_y_xyyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 702);

    auto ta1_y_xyyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 703);

    auto ta1_y_xyyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 704);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyzzz_xxxx_0, ta1_y_xyyzzz_xxxy_0, ta1_y_xyyzzz_xxxz_0, ta1_y_xyyzzz_xxyy_0, ta1_y_xyyzzz_xxyz_0, ta1_y_xyyzzz_xxzz_0, ta1_y_xyyzzz_xyyy_0, ta1_y_xyyzzz_xyyz_0, ta1_y_xyyzzz_xyzz_0, ta1_y_xyyzzz_xzzz_0, ta1_y_xyyzzz_yyyy_0, ta1_y_xyyzzz_yyyz_0, ta1_y_xyyzzz_yyzz_0, ta1_y_xyyzzz_yzzz_0, ta1_y_xyyzzz_zzzz_0, ta1_y_yyzzz_xxx_0, ta1_y_yyzzz_xxx_1, ta1_y_yyzzz_xxxx_0, ta1_y_yyzzz_xxxx_1, ta1_y_yyzzz_xxxy_0, ta1_y_yyzzz_xxxy_1, ta1_y_yyzzz_xxxz_0, ta1_y_yyzzz_xxxz_1, ta1_y_yyzzz_xxy_0, ta1_y_yyzzz_xxy_1, ta1_y_yyzzz_xxyy_0, ta1_y_yyzzz_xxyy_1, ta1_y_yyzzz_xxyz_0, ta1_y_yyzzz_xxyz_1, ta1_y_yyzzz_xxz_0, ta1_y_yyzzz_xxz_1, ta1_y_yyzzz_xxzz_0, ta1_y_yyzzz_xxzz_1, ta1_y_yyzzz_xyy_0, ta1_y_yyzzz_xyy_1, ta1_y_yyzzz_xyyy_0, ta1_y_yyzzz_xyyy_1, ta1_y_yyzzz_xyyz_0, ta1_y_yyzzz_xyyz_1, ta1_y_yyzzz_xyz_0, ta1_y_yyzzz_xyz_1, ta1_y_yyzzz_xyzz_0, ta1_y_yyzzz_xyzz_1, ta1_y_yyzzz_xzz_0, ta1_y_yyzzz_xzz_1, ta1_y_yyzzz_xzzz_0, ta1_y_yyzzz_xzzz_1, ta1_y_yyzzz_yyy_0, ta1_y_yyzzz_yyy_1, ta1_y_yyzzz_yyyy_0, ta1_y_yyzzz_yyyy_1, ta1_y_yyzzz_yyyz_0, ta1_y_yyzzz_yyyz_1, ta1_y_yyzzz_yyz_0, ta1_y_yyzzz_yyz_1, ta1_y_yyzzz_yyzz_0, ta1_y_yyzzz_yyzz_1, ta1_y_yyzzz_yzz_0, ta1_y_yyzzz_yzz_1, ta1_y_yyzzz_yzzz_0, ta1_y_yyzzz_yzzz_1, ta1_y_yyzzz_zzz_0, ta1_y_yyzzz_zzz_1, ta1_y_yyzzz_zzzz_0, ta1_y_yyzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzzz_xxxx_0[i] = 4.0 * ta1_y_yyzzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyzzz_xxx_1[i] * fe_0 + ta1_y_yyzzz_xxxx_0[i] * pa_x[i] - ta1_y_yyzzz_xxxx_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxxy_0[i] = 3.0 * ta1_y_yyzzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyzzz_xxy_1[i] * fe_0 + ta1_y_yyzzz_xxxy_0[i] * pa_x[i] - ta1_y_yyzzz_xxxy_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxxz_0[i] = 3.0 * ta1_y_yyzzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyzzz_xxz_1[i] * fe_0 + ta1_y_yyzzz_xxxz_0[i] * pa_x[i] - ta1_y_yyzzz_xxxz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxyy_0[i] = 2.0 * ta1_y_yyzzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_xyy_1[i] * fe_0 + ta1_y_yyzzz_xxyy_0[i] * pa_x[i] - ta1_y_yyzzz_xxyy_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxyz_0[i] = 2.0 * ta1_y_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_xyz_1[i] * fe_0 + ta1_y_yyzzz_xxyz_0[i] * pa_x[i] - ta1_y_yyzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxzz_0[i] = 2.0 * ta1_y_yyzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_xzz_1[i] * fe_0 + ta1_y_yyzzz_xxzz_0[i] * pa_x[i] - ta1_y_yyzzz_xxzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xyyy_0[i] = ta1_y_yyzzz_yyy_0[i] * fe_0 - ta1_y_yyzzz_yyy_1[i] * fe_0 + ta1_y_yyzzz_xyyy_0[i] * pa_x[i] - ta1_y_yyzzz_xyyy_1[i] * pc_x[i];

        ta1_y_xyyzzz_xyyz_0[i] = ta1_y_yyzzz_yyz_0[i] * fe_0 - ta1_y_yyzzz_yyz_1[i] * fe_0 + ta1_y_yyzzz_xyyz_0[i] * pa_x[i] - ta1_y_yyzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xyzz_0[i] = ta1_y_yyzzz_yzz_0[i] * fe_0 - ta1_y_yyzzz_yzz_1[i] * fe_0 + ta1_y_yyzzz_xyzz_0[i] * pa_x[i] - ta1_y_yyzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xzzz_0[i] = ta1_y_yyzzz_zzz_0[i] * fe_0 - ta1_y_yyzzz_zzz_1[i] * fe_0 + ta1_y_yyzzz_xzzz_0[i] * pa_x[i] - ta1_y_yyzzz_xzzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_yyyy_0[i] = ta1_y_yyzzz_yyyy_0[i] * pa_x[i] - ta1_y_yyzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xyyzzz_yyyz_0[i] = ta1_y_yyzzz_yyyz_0[i] * pa_x[i] - ta1_y_yyzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xyyzzz_yyzz_0[i] = ta1_y_yyzzz_yyzz_0[i] * pa_x[i] - ta1_y_yyzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_yzzz_0[i] = ta1_y_yyzzz_yzzz_0[i] * pa_x[i] - ta1_y_yyzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_zzzz_0[i] = ta1_y_yyzzz_zzzz_0[i] * pa_x[i] - ta1_y_yyzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 705-720 components of targeted buffer : IG

    auto ta1_y_xyzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 705);

    auto ta1_y_xyzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 706);

    auto ta1_y_xyzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 707);

    auto ta1_y_xyzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 708);

    auto ta1_y_xyzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 709);

    auto ta1_y_xyzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 710);

    auto ta1_y_xyzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 711);

    auto ta1_y_xyzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 712);

    auto ta1_y_xyzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 713);

    auto ta1_y_xyzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 714);

    auto ta1_y_xyzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 715);

    auto ta1_y_xyzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 716);

    auto ta1_y_xyzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 717);

    auto ta1_y_xyzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 718);

    auto ta1_y_xyzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 719);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xyzzzz_xxxx_0, ta1_y_xyzzzz_xxxy_0, ta1_y_xyzzzz_xxxz_0, ta1_y_xyzzzz_xxyy_0, ta1_y_xyzzzz_xxyz_0, ta1_y_xyzzzz_xxzz_0, ta1_y_xyzzzz_xyyy_0, ta1_y_xyzzzz_xyyz_0, ta1_y_xyzzzz_xyzz_0, ta1_y_xyzzzz_xzzz_0, ta1_y_xyzzzz_yyyy_0, ta1_y_xyzzzz_yyyz_0, ta1_y_xyzzzz_yyzz_0, ta1_y_xyzzzz_yzzz_0, ta1_y_xyzzzz_zzzz_0, ta1_y_xzzzz_xxxx_0, ta1_y_xzzzz_xxxx_1, ta1_y_xzzzz_xxxz_0, ta1_y_xzzzz_xxxz_1, ta1_y_xzzzz_xxzz_0, ta1_y_xzzzz_xxzz_1, ta1_y_xzzzz_xzzz_0, ta1_y_xzzzz_xzzz_1, ta1_y_yzzzz_xxxy_0, ta1_y_yzzzz_xxxy_1, ta1_y_yzzzz_xxy_0, ta1_y_yzzzz_xxy_1, ta1_y_yzzzz_xxyy_0, ta1_y_yzzzz_xxyy_1, ta1_y_yzzzz_xxyz_0, ta1_y_yzzzz_xxyz_1, ta1_y_yzzzz_xyy_0, ta1_y_yzzzz_xyy_1, ta1_y_yzzzz_xyyy_0, ta1_y_yzzzz_xyyy_1, ta1_y_yzzzz_xyyz_0, ta1_y_yzzzz_xyyz_1, ta1_y_yzzzz_xyz_0, ta1_y_yzzzz_xyz_1, ta1_y_yzzzz_xyzz_0, ta1_y_yzzzz_xyzz_1, ta1_y_yzzzz_yyy_0, ta1_y_yzzzz_yyy_1, ta1_y_yzzzz_yyyy_0, ta1_y_yzzzz_yyyy_1, ta1_y_yzzzz_yyyz_0, ta1_y_yzzzz_yyyz_1, ta1_y_yzzzz_yyz_0, ta1_y_yzzzz_yyz_1, ta1_y_yzzzz_yyzz_0, ta1_y_yzzzz_yyzz_1, ta1_y_yzzzz_yzz_0, ta1_y_yzzzz_yzz_1, ta1_y_yzzzz_yzzz_0, ta1_y_yzzzz_yzzz_1, ta1_y_yzzzz_zzzz_0, ta1_y_yzzzz_zzzz_1, ta_xzzzz_xxxx_1, ta_xzzzz_xxxz_1, ta_xzzzz_xxzz_1, ta_xzzzz_xzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzzzz_xxxx_0[i] = ta_xzzzz_xxxx_1[i] + ta1_y_xzzzz_xxxx_0[i] * pa_y[i] - ta1_y_xzzzz_xxxx_1[i] * pc_y[i];

        ta1_y_xyzzzz_xxxy_0[i] = 3.0 * ta1_y_yzzzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yzzzz_xxy_1[i] * fe_0 + ta1_y_yzzzz_xxxy_0[i] * pa_x[i] - ta1_y_yzzzz_xxxy_1[i] * pc_x[i];

        ta1_y_xyzzzz_xxxz_0[i] = ta_xzzzz_xxxz_1[i] + ta1_y_xzzzz_xxxz_0[i] * pa_y[i] - ta1_y_xzzzz_xxxz_1[i] * pc_y[i];

        ta1_y_xyzzzz_xxyy_0[i] = 2.0 * ta1_y_yzzzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yzzzz_xyy_1[i] * fe_0 + ta1_y_yzzzz_xxyy_0[i] * pa_x[i] - ta1_y_yzzzz_xxyy_1[i] * pc_x[i];

        ta1_y_xyzzzz_xxyz_0[i] = 2.0 * ta1_y_yzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yzzzz_xyz_1[i] * fe_0 + ta1_y_yzzzz_xxyz_0[i] * pa_x[i] - ta1_y_yzzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xyzzzz_xxzz_0[i] = ta_xzzzz_xxzz_1[i] + ta1_y_xzzzz_xxzz_0[i] * pa_y[i] - ta1_y_xzzzz_xxzz_1[i] * pc_y[i];

        ta1_y_xyzzzz_xyyy_0[i] = ta1_y_yzzzz_yyy_0[i] * fe_0 - ta1_y_yzzzz_yyy_1[i] * fe_0 + ta1_y_yzzzz_xyyy_0[i] * pa_x[i] - ta1_y_yzzzz_xyyy_1[i] * pc_x[i];

        ta1_y_xyzzzz_xyyz_0[i] = ta1_y_yzzzz_yyz_0[i] * fe_0 - ta1_y_yzzzz_yyz_1[i] * fe_0 + ta1_y_yzzzz_xyyz_0[i] * pa_x[i] - ta1_y_yzzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xyzzzz_xyzz_0[i] = ta1_y_yzzzz_yzz_0[i] * fe_0 - ta1_y_yzzzz_yzz_1[i] * fe_0 + ta1_y_yzzzz_xyzz_0[i] * pa_x[i] - ta1_y_yzzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xyzzzz_xzzz_0[i] = ta_xzzzz_xzzz_1[i] + ta1_y_xzzzz_xzzz_0[i] * pa_y[i] - ta1_y_xzzzz_xzzz_1[i] * pc_y[i];

        ta1_y_xyzzzz_yyyy_0[i] = ta1_y_yzzzz_yyyy_0[i] * pa_x[i] - ta1_y_yzzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xyzzzz_yyyz_0[i] = ta1_y_yzzzz_yyyz_0[i] * pa_x[i] - ta1_y_yzzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xyzzzz_yyzz_0[i] = ta1_y_yzzzz_yyzz_0[i] * pa_x[i] - ta1_y_yzzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xyzzzz_yzzz_0[i] = ta1_y_yzzzz_yzzz_0[i] * pa_x[i] - ta1_y_yzzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xyzzzz_zzzz_0[i] = ta1_y_yzzzz_zzzz_0[i] * pa_x[i] - ta1_y_yzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 720-735 components of targeted buffer : IG

    auto ta1_y_xzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 720);

    auto ta1_y_xzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 721);

    auto ta1_y_xzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 722);

    auto ta1_y_xzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 723);

    auto ta1_y_xzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 724);

    auto ta1_y_xzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 725);

    auto ta1_y_xzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 726);

    auto ta1_y_xzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 727);

    auto ta1_y_xzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 728);

    auto ta1_y_xzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 729);

    auto ta1_y_xzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 730);

    auto ta1_y_xzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 731);

    auto ta1_y_xzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 732);

    auto ta1_y_xzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 733);

    auto ta1_y_xzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 734);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xzzzzz_xxxx_0, ta1_y_xzzzzz_xxxy_0, ta1_y_xzzzzz_xxxz_0, ta1_y_xzzzzz_xxyy_0, ta1_y_xzzzzz_xxyz_0, ta1_y_xzzzzz_xxzz_0, ta1_y_xzzzzz_xyyy_0, ta1_y_xzzzzz_xyyz_0, ta1_y_xzzzzz_xyzz_0, ta1_y_xzzzzz_xzzz_0, ta1_y_xzzzzz_yyyy_0, ta1_y_xzzzzz_yyyz_0, ta1_y_xzzzzz_yyzz_0, ta1_y_xzzzzz_yzzz_0, ta1_y_xzzzzz_zzzz_0, ta1_y_zzzzz_xxx_0, ta1_y_zzzzz_xxx_1, ta1_y_zzzzz_xxxx_0, ta1_y_zzzzz_xxxx_1, ta1_y_zzzzz_xxxy_0, ta1_y_zzzzz_xxxy_1, ta1_y_zzzzz_xxxz_0, ta1_y_zzzzz_xxxz_1, ta1_y_zzzzz_xxy_0, ta1_y_zzzzz_xxy_1, ta1_y_zzzzz_xxyy_0, ta1_y_zzzzz_xxyy_1, ta1_y_zzzzz_xxyz_0, ta1_y_zzzzz_xxyz_1, ta1_y_zzzzz_xxz_0, ta1_y_zzzzz_xxz_1, ta1_y_zzzzz_xxzz_0, ta1_y_zzzzz_xxzz_1, ta1_y_zzzzz_xyy_0, ta1_y_zzzzz_xyy_1, ta1_y_zzzzz_xyyy_0, ta1_y_zzzzz_xyyy_1, ta1_y_zzzzz_xyyz_0, ta1_y_zzzzz_xyyz_1, ta1_y_zzzzz_xyz_0, ta1_y_zzzzz_xyz_1, ta1_y_zzzzz_xyzz_0, ta1_y_zzzzz_xyzz_1, ta1_y_zzzzz_xzz_0, ta1_y_zzzzz_xzz_1, ta1_y_zzzzz_xzzz_0, ta1_y_zzzzz_xzzz_1, ta1_y_zzzzz_yyy_0, ta1_y_zzzzz_yyy_1, ta1_y_zzzzz_yyyy_0, ta1_y_zzzzz_yyyy_1, ta1_y_zzzzz_yyyz_0, ta1_y_zzzzz_yyyz_1, ta1_y_zzzzz_yyz_0, ta1_y_zzzzz_yyz_1, ta1_y_zzzzz_yyzz_0, ta1_y_zzzzz_yyzz_1, ta1_y_zzzzz_yzz_0, ta1_y_zzzzz_yzz_1, ta1_y_zzzzz_yzzz_0, ta1_y_zzzzz_yzzz_1, ta1_y_zzzzz_zzz_0, ta1_y_zzzzz_zzz_1, ta1_y_zzzzz_zzzz_0, ta1_y_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzzz_xxxx_0[i] = 4.0 * ta1_y_zzzzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_zzzzz_xxx_1[i] * fe_0 + ta1_y_zzzzz_xxxx_0[i] * pa_x[i] - ta1_y_zzzzz_xxxx_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxxy_0[i] = 3.0 * ta1_y_zzzzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_zzzzz_xxy_1[i] * fe_0 + ta1_y_zzzzz_xxxy_0[i] * pa_x[i] - ta1_y_zzzzz_xxxy_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxxz_0[i] = 3.0 * ta1_y_zzzzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_zzzzz_xxz_1[i] * fe_0 + ta1_y_zzzzz_xxxz_0[i] * pa_x[i] - ta1_y_zzzzz_xxxz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxyy_0[i] = 2.0 * ta1_y_zzzzz_xyy_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xyy_1[i] * fe_0 + ta1_y_zzzzz_xxyy_0[i] * pa_x[i] - ta1_y_zzzzz_xxyy_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxyz_0[i] = 2.0 * ta1_y_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xyz_1[i] * fe_0 + ta1_y_zzzzz_xxyz_0[i] * pa_x[i] - ta1_y_zzzzz_xxyz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxzz_0[i] = 2.0 * ta1_y_zzzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xzz_1[i] * fe_0 + ta1_y_zzzzz_xxzz_0[i] * pa_x[i] - ta1_y_zzzzz_xxzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xyyy_0[i] = ta1_y_zzzzz_yyy_0[i] * fe_0 - ta1_y_zzzzz_yyy_1[i] * fe_0 + ta1_y_zzzzz_xyyy_0[i] * pa_x[i] - ta1_y_zzzzz_xyyy_1[i] * pc_x[i];

        ta1_y_xzzzzz_xyyz_0[i] = ta1_y_zzzzz_yyz_0[i] * fe_0 - ta1_y_zzzzz_yyz_1[i] * fe_0 + ta1_y_zzzzz_xyyz_0[i] * pa_x[i] - ta1_y_zzzzz_xyyz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xyzz_0[i] = ta1_y_zzzzz_yzz_0[i] * fe_0 - ta1_y_zzzzz_yzz_1[i] * fe_0 + ta1_y_zzzzz_xyzz_0[i] * pa_x[i] - ta1_y_zzzzz_xyzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xzzz_0[i] = ta1_y_zzzzz_zzz_0[i] * fe_0 - ta1_y_zzzzz_zzz_1[i] * fe_0 + ta1_y_zzzzz_xzzz_0[i] * pa_x[i] - ta1_y_zzzzz_xzzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_yyyy_0[i] = ta1_y_zzzzz_yyyy_0[i] * pa_x[i] - ta1_y_zzzzz_yyyy_1[i] * pc_x[i];

        ta1_y_xzzzzz_yyyz_0[i] = ta1_y_zzzzz_yyyz_0[i] * pa_x[i] - ta1_y_zzzzz_yyyz_1[i] * pc_x[i];

        ta1_y_xzzzzz_yyzz_0[i] = ta1_y_zzzzz_yyzz_0[i] * pa_x[i] - ta1_y_zzzzz_yyzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_yzzz_0[i] = ta1_y_zzzzz_yzzz_0[i] * pa_x[i] - ta1_y_zzzzz_yzzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_zzzz_0[i] = ta1_y_zzzzz_zzzz_0[i] * pa_x[i] - ta1_y_zzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 735-750 components of targeted buffer : IG

    auto ta1_y_yyyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 735);

    auto ta1_y_yyyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 736);

    auto ta1_y_yyyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 737);

    auto ta1_y_yyyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 738);

    auto ta1_y_yyyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 739);

    auto ta1_y_yyyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 740);

    auto ta1_y_yyyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 741);

    auto ta1_y_yyyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 742);

    auto ta1_y_yyyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 743);

    auto ta1_y_yyyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 744);

    auto ta1_y_yyyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 745);

    auto ta1_y_yyyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 746);

    auto ta1_y_yyyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 747);

    auto ta1_y_yyyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 748);

    auto ta1_y_yyyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 749);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yyyy_xxxx_0, ta1_y_yyyy_xxxx_1, ta1_y_yyyy_xxxy_0, ta1_y_yyyy_xxxy_1, ta1_y_yyyy_xxxz_0, ta1_y_yyyy_xxxz_1, ta1_y_yyyy_xxyy_0, ta1_y_yyyy_xxyy_1, ta1_y_yyyy_xxyz_0, ta1_y_yyyy_xxyz_1, ta1_y_yyyy_xxzz_0, ta1_y_yyyy_xxzz_1, ta1_y_yyyy_xyyy_0, ta1_y_yyyy_xyyy_1, ta1_y_yyyy_xyyz_0, ta1_y_yyyy_xyyz_1, ta1_y_yyyy_xyzz_0, ta1_y_yyyy_xyzz_1, ta1_y_yyyy_xzzz_0, ta1_y_yyyy_xzzz_1, ta1_y_yyyy_yyyy_0, ta1_y_yyyy_yyyy_1, ta1_y_yyyy_yyyz_0, ta1_y_yyyy_yyyz_1, ta1_y_yyyy_yyzz_0, ta1_y_yyyy_yyzz_1, ta1_y_yyyy_yzzz_0, ta1_y_yyyy_yzzz_1, ta1_y_yyyy_zzzz_0, ta1_y_yyyy_zzzz_1, ta1_y_yyyyy_xxx_0, ta1_y_yyyyy_xxx_1, ta1_y_yyyyy_xxxx_0, ta1_y_yyyyy_xxxx_1, ta1_y_yyyyy_xxxy_0, ta1_y_yyyyy_xxxy_1, ta1_y_yyyyy_xxxz_0, ta1_y_yyyyy_xxxz_1, ta1_y_yyyyy_xxy_0, ta1_y_yyyyy_xxy_1, ta1_y_yyyyy_xxyy_0, ta1_y_yyyyy_xxyy_1, ta1_y_yyyyy_xxyz_0, ta1_y_yyyyy_xxyz_1, ta1_y_yyyyy_xxz_0, ta1_y_yyyyy_xxz_1, ta1_y_yyyyy_xxzz_0, ta1_y_yyyyy_xxzz_1, ta1_y_yyyyy_xyy_0, ta1_y_yyyyy_xyy_1, ta1_y_yyyyy_xyyy_0, ta1_y_yyyyy_xyyy_1, ta1_y_yyyyy_xyyz_0, ta1_y_yyyyy_xyyz_1, ta1_y_yyyyy_xyz_0, ta1_y_yyyyy_xyz_1, ta1_y_yyyyy_xyzz_0, ta1_y_yyyyy_xyzz_1, ta1_y_yyyyy_xzz_0, ta1_y_yyyyy_xzz_1, ta1_y_yyyyy_xzzz_0, ta1_y_yyyyy_xzzz_1, ta1_y_yyyyy_yyy_0, ta1_y_yyyyy_yyy_1, ta1_y_yyyyy_yyyy_0, ta1_y_yyyyy_yyyy_1, ta1_y_yyyyy_yyyz_0, ta1_y_yyyyy_yyyz_1, ta1_y_yyyyy_yyz_0, ta1_y_yyyyy_yyz_1, ta1_y_yyyyy_yyzz_0, ta1_y_yyyyy_yyzz_1, ta1_y_yyyyy_yzz_0, ta1_y_yyyyy_yzz_1, ta1_y_yyyyy_yzzz_0, ta1_y_yyyyy_yzzz_1, ta1_y_yyyyy_zzz_0, ta1_y_yyyyy_zzz_1, ta1_y_yyyyy_zzzz_0, ta1_y_yyyyy_zzzz_1, ta1_y_yyyyyy_xxxx_0, ta1_y_yyyyyy_xxxy_0, ta1_y_yyyyyy_xxxz_0, ta1_y_yyyyyy_xxyy_0, ta1_y_yyyyyy_xxyz_0, ta1_y_yyyyyy_xxzz_0, ta1_y_yyyyyy_xyyy_0, ta1_y_yyyyyy_xyyz_0, ta1_y_yyyyyy_xyzz_0, ta1_y_yyyyyy_xzzz_0, ta1_y_yyyyyy_yyyy_0, ta1_y_yyyyyy_yyyz_0, ta1_y_yyyyyy_yyzz_0, ta1_y_yyyyyy_yzzz_0, ta1_y_yyyyyy_zzzz_0, ta_yyyyy_xxxx_1, ta_yyyyy_xxxy_1, ta_yyyyy_xxxz_1, ta_yyyyy_xxyy_1, ta_yyyyy_xxyz_1, ta_yyyyy_xxzz_1, ta_yyyyy_xyyy_1, ta_yyyyy_xyyz_1, ta_yyyyy_xyzz_1, ta_yyyyy_xzzz_1, ta_yyyyy_yyyy_1, ta_yyyyy_yyyz_1, ta_yyyyy_yyzz_1, ta_yyyyy_yzzz_1, ta_yyyyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyy_xxxx_0[i] = 5.0 * ta1_y_yyyy_xxxx_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxxx_1[i] * fe_0 + ta_yyyyy_xxxx_1[i] + ta1_y_yyyyy_xxxx_0[i] * pa_y[i] - ta1_y_yyyyy_xxxx_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxxy_0[i] = 5.0 * ta1_y_yyyy_xxxy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxxy_1[i] * fe_0 + ta1_y_yyyyy_xxx_0[i] * fe_0 - ta1_y_yyyyy_xxx_1[i] * fe_0 + ta_yyyyy_xxxy_1[i] + ta1_y_yyyyy_xxxy_0[i] * pa_y[i] - ta1_y_yyyyy_xxxy_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxxz_0[i] = 5.0 * ta1_y_yyyy_xxxz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxxz_1[i] * fe_0 + ta_yyyyy_xxxz_1[i] + ta1_y_yyyyy_xxxz_0[i] * pa_y[i] - ta1_y_yyyyy_xxxz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxyy_0[i] = 5.0 * ta1_y_yyyy_xxyy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_y_yyyyy_xxy_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xxy_1[i] * fe_0 + ta_yyyyy_xxyy_1[i] + ta1_y_yyyyy_xxyy_0[i] * pa_y[i] - ta1_y_yyyyy_xxyy_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxyz_0[i] = 5.0 * ta1_y_yyyy_xxyz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxyz_1[i] * fe_0 + ta1_y_yyyyy_xxz_0[i] * fe_0 - ta1_y_yyyyy_xxz_1[i] * fe_0 + ta_yyyyy_xxyz_1[i] + ta1_y_yyyyy_xxyz_0[i] * pa_y[i] - ta1_y_yyyyy_xxyz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxzz_0[i] = 5.0 * ta1_y_yyyy_xxzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxzz_1[i] * fe_0 + ta_yyyyy_xxzz_1[i] + ta1_y_yyyyy_xxzz_0[i] * pa_y[i] - ta1_y_yyyyy_xxzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xyyy_0[i] = 5.0 * ta1_y_yyyy_xyyy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xyyy_1[i] * fe_0 + 3.0 * ta1_y_yyyyy_xyy_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_xyy_1[i] * fe_0 + ta_yyyyy_xyyy_1[i] + ta1_y_yyyyy_xyyy_0[i] * pa_y[i] - ta1_y_yyyyy_xyyy_1[i] * pc_y[i];

        ta1_y_yyyyyy_xyyz_0[i] = 5.0 * ta1_y_yyyy_xyyz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xyyz_1[i] * fe_0 + 2.0 * ta1_y_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xyz_1[i] * fe_0 + ta_yyyyy_xyyz_1[i] + ta1_y_yyyyy_xyyz_0[i] * pa_y[i] - ta1_y_yyyyy_xyyz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xyzz_0[i] = 5.0 * ta1_y_yyyy_xyzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xyzz_1[i] * fe_0 + ta1_y_yyyyy_xzz_0[i] * fe_0 - ta1_y_yyyyy_xzz_1[i] * fe_0 + ta_yyyyy_xyzz_1[i] + ta1_y_yyyyy_xyzz_0[i] * pa_y[i] - ta1_y_yyyyy_xyzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xzzz_0[i] = 5.0 * ta1_y_yyyy_xzzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xzzz_1[i] * fe_0 + ta_yyyyy_xzzz_1[i] + ta1_y_yyyyy_xzzz_0[i] * pa_y[i] - ta1_y_yyyyy_xzzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_yyyy_0[i] = 5.0 * ta1_y_yyyy_yyyy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yyyy_1[i] * fe_0 + 4.0 * ta1_y_yyyyy_yyy_0[i] * fe_0 - 4.0 * ta1_y_yyyyy_yyy_1[i] * fe_0 + ta_yyyyy_yyyy_1[i] + ta1_y_yyyyy_yyyy_0[i] * pa_y[i] - ta1_y_yyyyy_yyyy_1[i] * pc_y[i];

        ta1_y_yyyyyy_yyyz_0[i] = 5.0 * ta1_y_yyyy_yyyz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yyyz_1[i] * fe_0 + 3.0 * ta1_y_yyyyy_yyz_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_yyz_1[i] * fe_0 + ta_yyyyy_yyyz_1[i] + ta1_y_yyyyy_yyyz_0[i] * pa_y[i] - ta1_y_yyyyy_yyyz_1[i] * pc_y[i];

        ta1_y_yyyyyy_yyzz_0[i] = 5.0 * ta1_y_yyyy_yyzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyyy_yzz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_yzz_1[i] * fe_0 + ta_yyyyy_yyzz_1[i] + ta1_y_yyyyy_yyzz_0[i] * pa_y[i] - ta1_y_yyyyy_yyzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_yzzz_0[i] = 5.0 * ta1_y_yyyy_yzzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yzzz_1[i] * fe_0 + ta1_y_yyyyy_zzz_0[i] * fe_0 - ta1_y_yyyyy_zzz_1[i] * fe_0 + ta_yyyyy_yzzz_1[i] + ta1_y_yyyyy_yzzz_0[i] * pa_y[i] - ta1_y_yyyyy_yzzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_zzzz_0[i] = 5.0 * ta1_y_yyyy_zzzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_zzzz_1[i] * fe_0 + ta_yyyyy_zzzz_1[i] + ta1_y_yyyyy_zzzz_0[i] * pa_y[i] - ta1_y_yyyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 750-765 components of targeted buffer : IG

    auto ta1_y_yyyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 750);

    auto ta1_y_yyyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 751);

    auto ta1_y_yyyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 752);

    auto ta1_y_yyyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 753);

    auto ta1_y_yyyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 754);

    auto ta1_y_yyyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 755);

    auto ta1_y_yyyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 756);

    auto ta1_y_yyyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 757);

    auto ta1_y_yyyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 758);

    auto ta1_y_yyyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 759);

    auto ta1_y_yyyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 760);

    auto ta1_y_yyyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 761);

    auto ta1_y_yyyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 762);

    auto ta1_y_yyyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 763);

    auto ta1_y_yyyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 764);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_yyyyy_xxx_0, ta1_y_yyyyy_xxx_1, ta1_y_yyyyy_xxxx_0, ta1_y_yyyyy_xxxx_1, ta1_y_yyyyy_xxxy_0, ta1_y_yyyyy_xxxy_1, ta1_y_yyyyy_xxxz_0, ta1_y_yyyyy_xxxz_1, ta1_y_yyyyy_xxy_0, ta1_y_yyyyy_xxy_1, ta1_y_yyyyy_xxyy_0, ta1_y_yyyyy_xxyy_1, ta1_y_yyyyy_xxyz_0, ta1_y_yyyyy_xxyz_1, ta1_y_yyyyy_xxz_0, ta1_y_yyyyy_xxz_1, ta1_y_yyyyy_xxzz_0, ta1_y_yyyyy_xxzz_1, ta1_y_yyyyy_xyy_0, ta1_y_yyyyy_xyy_1, ta1_y_yyyyy_xyyy_0, ta1_y_yyyyy_xyyy_1, ta1_y_yyyyy_xyyz_0, ta1_y_yyyyy_xyyz_1, ta1_y_yyyyy_xyz_0, ta1_y_yyyyy_xyz_1, ta1_y_yyyyy_xyzz_0, ta1_y_yyyyy_xyzz_1, ta1_y_yyyyy_xzz_0, ta1_y_yyyyy_xzz_1, ta1_y_yyyyy_xzzz_0, ta1_y_yyyyy_xzzz_1, ta1_y_yyyyy_yyy_0, ta1_y_yyyyy_yyy_1, ta1_y_yyyyy_yyyy_0, ta1_y_yyyyy_yyyy_1, ta1_y_yyyyy_yyyz_0, ta1_y_yyyyy_yyyz_1, ta1_y_yyyyy_yyz_0, ta1_y_yyyyy_yyz_1, ta1_y_yyyyy_yyzz_0, ta1_y_yyyyy_yyzz_1, ta1_y_yyyyy_yzz_0, ta1_y_yyyyy_yzz_1, ta1_y_yyyyy_yzzz_0, ta1_y_yyyyy_yzzz_1, ta1_y_yyyyy_zzz_0, ta1_y_yyyyy_zzz_1, ta1_y_yyyyy_zzzz_0, ta1_y_yyyyy_zzzz_1, ta1_y_yyyyyz_xxxx_0, ta1_y_yyyyyz_xxxy_0, ta1_y_yyyyyz_xxxz_0, ta1_y_yyyyyz_xxyy_0, ta1_y_yyyyyz_xxyz_0, ta1_y_yyyyyz_xxzz_0, ta1_y_yyyyyz_xyyy_0, ta1_y_yyyyyz_xyyz_0, ta1_y_yyyyyz_xyzz_0, ta1_y_yyyyyz_xzzz_0, ta1_y_yyyyyz_yyyy_0, ta1_y_yyyyyz_yyyz_0, ta1_y_yyyyyz_yyzz_0, ta1_y_yyyyyz_yzzz_0, ta1_y_yyyyyz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyz_xxxx_0[i] = ta1_y_yyyyy_xxxx_0[i] * pa_z[i] - ta1_y_yyyyy_xxxx_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxxy_0[i] = ta1_y_yyyyy_xxxy_0[i] * pa_z[i] - ta1_y_yyyyy_xxxy_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxxz_0[i] = ta1_y_yyyyy_xxx_0[i] * fe_0 - ta1_y_yyyyy_xxx_1[i] * fe_0 + ta1_y_yyyyy_xxxz_0[i] * pa_z[i] - ta1_y_yyyyy_xxxz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxyy_0[i] = ta1_y_yyyyy_xxyy_0[i] * pa_z[i] - ta1_y_yyyyy_xxyy_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxyz_0[i] = ta1_y_yyyyy_xxy_0[i] * fe_0 - ta1_y_yyyyy_xxy_1[i] * fe_0 + ta1_y_yyyyy_xxyz_0[i] * pa_z[i] - ta1_y_yyyyy_xxyz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxzz_0[i] = 2.0 * ta1_y_yyyyy_xxz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xxz_1[i] * fe_0 + ta1_y_yyyyy_xxzz_0[i] * pa_z[i] - ta1_y_yyyyy_xxzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xyyy_0[i] = ta1_y_yyyyy_xyyy_0[i] * pa_z[i] - ta1_y_yyyyy_xyyy_1[i] * pc_z[i];

        ta1_y_yyyyyz_xyyz_0[i] = ta1_y_yyyyy_xyy_0[i] * fe_0 - ta1_y_yyyyy_xyy_1[i] * fe_0 + ta1_y_yyyyy_xyyz_0[i] * pa_z[i] - ta1_y_yyyyy_xyyz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xyzz_0[i] = 2.0 * ta1_y_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xyz_1[i] * fe_0 + ta1_y_yyyyy_xyzz_0[i] * pa_z[i] - ta1_y_yyyyy_xyzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xzzz_0[i] = 3.0 * ta1_y_yyyyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_xzz_1[i] * fe_0 + ta1_y_yyyyy_xzzz_0[i] * pa_z[i] - ta1_y_yyyyy_xzzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_yyyy_0[i] = ta1_y_yyyyy_yyyy_0[i] * pa_z[i] - ta1_y_yyyyy_yyyy_1[i] * pc_z[i];

        ta1_y_yyyyyz_yyyz_0[i] = ta1_y_yyyyy_yyy_0[i] * fe_0 - ta1_y_yyyyy_yyy_1[i] * fe_0 + ta1_y_yyyyy_yyyz_0[i] * pa_z[i] - ta1_y_yyyyy_yyyz_1[i] * pc_z[i];

        ta1_y_yyyyyz_yyzz_0[i] = 2.0 * ta1_y_yyyyy_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_yyz_1[i] * fe_0 + ta1_y_yyyyy_yyzz_0[i] * pa_z[i] - ta1_y_yyyyy_yyzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_yzzz_0[i] = 3.0 * ta1_y_yyyyy_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_yzz_1[i] * fe_0 + ta1_y_yyyyy_yzzz_0[i] * pa_z[i] - ta1_y_yyyyy_yzzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_zzzz_0[i] = 4.0 * ta1_y_yyyyy_zzz_0[i] * fe_0 - 4.0 * ta1_y_yyyyy_zzz_1[i] * fe_0 + ta1_y_yyyyy_zzzz_0[i] * pa_z[i] - ta1_y_yyyyy_zzzz_1[i] * pc_z[i];
    }

    // Set up 765-780 components of targeted buffer : IG

    auto ta1_y_yyyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 765);

    auto ta1_y_yyyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 766);

    auto ta1_y_yyyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 767);

    auto ta1_y_yyyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 768);

    auto ta1_y_yyyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 769);

    auto ta1_y_yyyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 770);

    auto ta1_y_yyyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 771);

    auto ta1_y_yyyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 772);

    auto ta1_y_yyyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 773);

    auto ta1_y_yyyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 774);

    auto ta1_y_yyyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 775);

    auto ta1_y_yyyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 776);

    auto ta1_y_yyyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 777);

    auto ta1_y_yyyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 778);

    auto ta1_y_yyyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 779);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyyy_xxxx_0, ta1_y_yyyy_xxxx_1, ta1_y_yyyy_xxxy_0, ta1_y_yyyy_xxxy_1, ta1_y_yyyy_xxyy_0, ta1_y_yyyy_xxyy_1, ta1_y_yyyy_xxyz_0, ta1_y_yyyy_xxyz_1, ta1_y_yyyy_xyyy_0, ta1_y_yyyy_xyyy_1, ta1_y_yyyy_xyyz_0, ta1_y_yyyy_xyyz_1, ta1_y_yyyy_xyzz_0, ta1_y_yyyy_xyzz_1, ta1_y_yyyy_yyyy_0, ta1_y_yyyy_yyyy_1, ta1_y_yyyy_yyyz_0, ta1_y_yyyy_yyyz_1, ta1_y_yyyy_yyzz_0, ta1_y_yyyy_yyzz_1, ta1_y_yyyy_yzzz_0, ta1_y_yyyy_yzzz_1, ta1_y_yyyyz_xxxx_0, ta1_y_yyyyz_xxxx_1, ta1_y_yyyyz_xxxy_0, ta1_y_yyyyz_xxxy_1, ta1_y_yyyyz_xxy_0, ta1_y_yyyyz_xxy_1, ta1_y_yyyyz_xxyy_0, ta1_y_yyyyz_xxyy_1, ta1_y_yyyyz_xxyz_0, ta1_y_yyyyz_xxyz_1, ta1_y_yyyyz_xyy_0, ta1_y_yyyyz_xyy_1, ta1_y_yyyyz_xyyy_0, ta1_y_yyyyz_xyyy_1, ta1_y_yyyyz_xyyz_0, ta1_y_yyyyz_xyyz_1, ta1_y_yyyyz_xyz_0, ta1_y_yyyyz_xyz_1, ta1_y_yyyyz_xyzz_0, ta1_y_yyyyz_xyzz_1, ta1_y_yyyyz_yyy_0, ta1_y_yyyyz_yyy_1, ta1_y_yyyyz_yyyy_0, ta1_y_yyyyz_yyyy_1, ta1_y_yyyyz_yyyz_0, ta1_y_yyyyz_yyyz_1, ta1_y_yyyyz_yyz_0, ta1_y_yyyyz_yyz_1, ta1_y_yyyyz_yyzz_0, ta1_y_yyyyz_yyzz_1, ta1_y_yyyyz_yzz_0, ta1_y_yyyyz_yzz_1, ta1_y_yyyyz_yzzz_0, ta1_y_yyyyz_yzzz_1, ta1_y_yyyyzz_xxxx_0, ta1_y_yyyyzz_xxxy_0, ta1_y_yyyyzz_xxxz_0, ta1_y_yyyyzz_xxyy_0, ta1_y_yyyyzz_xxyz_0, ta1_y_yyyyzz_xxzz_0, ta1_y_yyyyzz_xyyy_0, ta1_y_yyyyzz_xyyz_0, ta1_y_yyyyzz_xyzz_0, ta1_y_yyyyzz_xzzz_0, ta1_y_yyyyzz_yyyy_0, ta1_y_yyyyzz_yyyz_0, ta1_y_yyyyzz_yyzz_0, ta1_y_yyyyzz_yzzz_0, ta1_y_yyyyzz_zzzz_0, ta1_y_yyyzz_xxxz_0, ta1_y_yyyzz_xxxz_1, ta1_y_yyyzz_xxzz_0, ta1_y_yyyzz_xxzz_1, ta1_y_yyyzz_xzzz_0, ta1_y_yyyzz_xzzz_1, ta1_y_yyyzz_zzzz_0, ta1_y_yyyzz_zzzz_1, ta1_y_yyzz_xxxz_0, ta1_y_yyzz_xxxz_1, ta1_y_yyzz_xxzz_0, ta1_y_yyzz_xxzz_1, ta1_y_yyzz_xzzz_0, ta1_y_yyzz_xzzz_1, ta1_y_yyzz_zzzz_0, ta1_y_yyzz_zzzz_1, ta_yyyzz_xxxz_1, ta_yyyzz_xxzz_1, ta_yyyzz_xzzz_1, ta_yyyzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyzz_xxxx_0[i] = ta1_y_yyyy_xxxx_0[i] * fe_0 - ta1_y_yyyy_xxxx_1[i] * fe_0 + ta1_y_yyyyz_xxxx_0[i] * pa_z[i] - ta1_y_yyyyz_xxxx_1[i] * pc_z[i];

        ta1_y_yyyyzz_xxxy_0[i] = ta1_y_yyyy_xxxy_0[i] * fe_0 - ta1_y_yyyy_xxxy_1[i] * fe_0 + ta1_y_yyyyz_xxxy_0[i] * pa_z[i] - ta1_y_yyyyz_xxxy_1[i] * pc_z[i];

        ta1_y_yyyyzz_xxxz_0[i] = 3.0 * ta1_y_yyzz_xxxz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxxz_1[i] * fe_0 + ta_yyyzz_xxxz_1[i] + ta1_y_yyyzz_xxxz_0[i] * pa_y[i] - ta1_y_yyyzz_xxxz_1[i] * pc_y[i];

        ta1_y_yyyyzz_xxyy_0[i] = ta1_y_yyyy_xxyy_0[i] * fe_0 - ta1_y_yyyy_xxyy_1[i] * fe_0 + ta1_y_yyyyz_xxyy_0[i] * pa_z[i] - ta1_y_yyyyz_xxyy_1[i] * pc_z[i];

        ta1_y_yyyyzz_xxyz_0[i] = ta1_y_yyyy_xxyz_0[i] * fe_0 - ta1_y_yyyy_xxyz_1[i] * fe_0 + ta1_y_yyyyz_xxy_0[i] * fe_0 - ta1_y_yyyyz_xxy_1[i] * fe_0 + ta1_y_yyyyz_xxyz_0[i] * pa_z[i] - ta1_y_yyyyz_xxyz_1[i] * pc_z[i];

        ta1_y_yyyyzz_xxzz_0[i] = 3.0 * ta1_y_yyzz_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxzz_1[i] * fe_0 + ta_yyyzz_xxzz_1[i] + ta1_y_yyyzz_xxzz_0[i] * pa_y[i] - ta1_y_yyyzz_xxzz_1[i] * pc_y[i];

        ta1_y_yyyyzz_xyyy_0[i] = ta1_y_yyyy_xyyy_0[i] * fe_0 - ta1_y_yyyy_xyyy_1[i] * fe_0 + ta1_y_yyyyz_xyyy_0[i] * pa_z[i] - ta1_y_yyyyz_xyyy_1[i] * pc_z[i];

        ta1_y_yyyyzz_xyyz_0[i] = ta1_y_yyyy_xyyz_0[i] * fe_0 - ta1_y_yyyy_xyyz_1[i] * fe_0 + ta1_y_yyyyz_xyy_0[i] * fe_0 - ta1_y_yyyyz_xyy_1[i] * fe_0 + ta1_y_yyyyz_xyyz_0[i] * pa_z[i] - ta1_y_yyyyz_xyyz_1[i] * pc_z[i];

        ta1_y_yyyyzz_xyzz_0[i] = ta1_y_yyyy_xyzz_0[i] * fe_0 - ta1_y_yyyy_xyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyyz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyz_xyz_1[i] * fe_0 + ta1_y_yyyyz_xyzz_0[i] * pa_z[i] - ta1_y_yyyyz_xyzz_1[i] * pc_z[i];

        ta1_y_yyyyzz_xzzz_0[i] = 3.0 * ta1_y_yyzz_xzzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xzzz_1[i] * fe_0 + ta_yyyzz_xzzz_1[i] + ta1_y_yyyzz_xzzz_0[i] * pa_y[i] - ta1_y_yyyzz_xzzz_1[i] * pc_y[i];

        ta1_y_yyyyzz_yyyy_0[i] = ta1_y_yyyy_yyyy_0[i] * fe_0 - ta1_y_yyyy_yyyy_1[i] * fe_0 + ta1_y_yyyyz_yyyy_0[i] * pa_z[i] - ta1_y_yyyyz_yyyy_1[i] * pc_z[i];

        ta1_y_yyyyzz_yyyz_0[i] = ta1_y_yyyy_yyyz_0[i] * fe_0 - ta1_y_yyyy_yyyz_1[i] * fe_0 + ta1_y_yyyyz_yyy_0[i] * fe_0 - ta1_y_yyyyz_yyy_1[i] * fe_0 + ta1_y_yyyyz_yyyz_0[i] * pa_z[i] - ta1_y_yyyyz_yyyz_1[i] * pc_z[i];

        ta1_y_yyyyzz_yyzz_0[i] = ta1_y_yyyy_yyzz_0[i] * fe_0 - ta1_y_yyyy_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyyz_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyyyz_yyz_1[i] * fe_0 + ta1_y_yyyyz_yyzz_0[i] * pa_z[i] - ta1_y_yyyyz_yyzz_1[i] * pc_z[i];

        ta1_y_yyyyzz_yzzz_0[i] = ta1_y_yyyy_yzzz_0[i] * fe_0 - ta1_y_yyyy_yzzz_1[i] * fe_0 + 3.0 * ta1_y_yyyyz_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyyyz_yzz_1[i] * fe_0 + ta1_y_yyyyz_yzzz_0[i] * pa_z[i] - ta1_y_yyyyz_yzzz_1[i] * pc_z[i];

        ta1_y_yyyyzz_zzzz_0[i] = 3.0 * ta1_y_yyzz_zzzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_zzzz_1[i] * fe_0 + ta_yyyzz_zzzz_1[i] + ta1_y_yyyzz_zzzz_0[i] * pa_y[i] - ta1_y_yyyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 780-795 components of targeted buffer : IG

    auto ta1_y_yyyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 780);

    auto ta1_y_yyyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 781);

    auto ta1_y_yyyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 782);

    auto ta1_y_yyyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 783);

    auto ta1_y_yyyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 784);

    auto ta1_y_yyyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 785);

    auto ta1_y_yyyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 786);

    auto ta1_y_yyyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 787);

    auto ta1_y_yyyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 788);

    auto ta1_y_yyyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 789);

    auto ta1_y_yyyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 790);

    auto ta1_y_yyyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 791);

    auto ta1_y_yyyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 792);

    auto ta1_y_yyyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 793);

    auto ta1_y_yyyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 794);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyyz_xxxx_0, ta1_y_yyyz_xxxx_1, ta1_y_yyyz_xxxy_0, ta1_y_yyyz_xxxy_1, ta1_y_yyyz_xxyy_0, ta1_y_yyyz_xxyy_1, ta1_y_yyyz_xxyz_0, ta1_y_yyyz_xxyz_1, ta1_y_yyyz_xyyy_0, ta1_y_yyyz_xyyy_1, ta1_y_yyyz_xyyz_0, ta1_y_yyyz_xyyz_1, ta1_y_yyyz_xyzz_0, ta1_y_yyyz_xyzz_1, ta1_y_yyyz_yyyy_0, ta1_y_yyyz_yyyy_1, ta1_y_yyyz_yyyz_0, ta1_y_yyyz_yyyz_1, ta1_y_yyyz_yyzz_0, ta1_y_yyyz_yyzz_1, ta1_y_yyyz_yzzz_0, ta1_y_yyyz_yzzz_1, ta1_y_yyyzz_xxxx_0, ta1_y_yyyzz_xxxx_1, ta1_y_yyyzz_xxxy_0, ta1_y_yyyzz_xxxy_1, ta1_y_yyyzz_xxy_0, ta1_y_yyyzz_xxy_1, ta1_y_yyyzz_xxyy_0, ta1_y_yyyzz_xxyy_1, ta1_y_yyyzz_xxyz_0, ta1_y_yyyzz_xxyz_1, ta1_y_yyyzz_xyy_0, ta1_y_yyyzz_xyy_1, ta1_y_yyyzz_xyyy_0, ta1_y_yyyzz_xyyy_1, ta1_y_yyyzz_xyyz_0, ta1_y_yyyzz_xyyz_1, ta1_y_yyyzz_xyz_0, ta1_y_yyyzz_xyz_1, ta1_y_yyyzz_xyzz_0, ta1_y_yyyzz_xyzz_1, ta1_y_yyyzz_yyy_0, ta1_y_yyyzz_yyy_1, ta1_y_yyyzz_yyyy_0, ta1_y_yyyzz_yyyy_1, ta1_y_yyyzz_yyyz_0, ta1_y_yyyzz_yyyz_1, ta1_y_yyyzz_yyz_0, ta1_y_yyyzz_yyz_1, ta1_y_yyyzz_yyzz_0, ta1_y_yyyzz_yyzz_1, ta1_y_yyyzz_yzz_0, ta1_y_yyyzz_yzz_1, ta1_y_yyyzz_yzzz_0, ta1_y_yyyzz_yzzz_1, ta1_y_yyyzzz_xxxx_0, ta1_y_yyyzzz_xxxy_0, ta1_y_yyyzzz_xxxz_0, ta1_y_yyyzzz_xxyy_0, ta1_y_yyyzzz_xxyz_0, ta1_y_yyyzzz_xxzz_0, ta1_y_yyyzzz_xyyy_0, ta1_y_yyyzzz_xyyz_0, ta1_y_yyyzzz_xyzz_0, ta1_y_yyyzzz_xzzz_0, ta1_y_yyyzzz_yyyy_0, ta1_y_yyyzzz_yyyz_0, ta1_y_yyyzzz_yyzz_0, ta1_y_yyyzzz_yzzz_0, ta1_y_yyyzzz_zzzz_0, ta1_y_yyzzz_xxxz_0, ta1_y_yyzzz_xxxz_1, ta1_y_yyzzz_xxzz_0, ta1_y_yyzzz_xxzz_1, ta1_y_yyzzz_xzzz_0, ta1_y_yyzzz_xzzz_1, ta1_y_yyzzz_zzzz_0, ta1_y_yyzzz_zzzz_1, ta1_y_yzzz_xxxz_0, ta1_y_yzzz_xxxz_1, ta1_y_yzzz_xxzz_0, ta1_y_yzzz_xxzz_1, ta1_y_yzzz_xzzz_0, ta1_y_yzzz_xzzz_1, ta1_y_yzzz_zzzz_0, ta1_y_yzzz_zzzz_1, ta_yyzzz_xxxz_1, ta_yyzzz_xxzz_1, ta_yyzzz_xzzz_1, ta_yyzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzzz_xxxx_0[i] = 2.0 * ta1_y_yyyz_xxxx_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xxxx_1[i] * fe_0 + ta1_y_yyyzz_xxxx_0[i] * pa_z[i] - ta1_y_yyyzz_xxxx_1[i] * pc_z[i];

        ta1_y_yyyzzz_xxxy_0[i] = 2.0 * ta1_y_yyyz_xxxy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xxxy_1[i] * fe_0 + ta1_y_yyyzz_xxxy_0[i] * pa_z[i] - ta1_y_yyyzz_xxxy_1[i] * pc_z[i];

        ta1_y_yyyzzz_xxxz_0[i] = 2.0 * ta1_y_yzzz_xxxz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xxxz_1[i] * fe_0 + ta_yyzzz_xxxz_1[i] + ta1_y_yyzzz_xxxz_0[i] * pa_y[i] - ta1_y_yyzzz_xxxz_1[i] * pc_y[i];

        ta1_y_yyyzzz_xxyy_0[i] = 2.0 * ta1_y_yyyz_xxyy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xxyy_1[i] * fe_0 + ta1_y_yyyzz_xxyy_0[i] * pa_z[i] - ta1_y_yyyzz_xxyy_1[i] * pc_z[i];

        ta1_y_yyyzzz_xxyz_0[i] = 2.0 * ta1_y_yyyz_xxyz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xxyz_1[i] * fe_0 + ta1_y_yyyzz_xxy_0[i] * fe_0 - ta1_y_yyyzz_xxy_1[i] * fe_0 + ta1_y_yyyzz_xxyz_0[i] * pa_z[i] - ta1_y_yyyzz_xxyz_1[i] * pc_z[i];

        ta1_y_yyyzzz_xxzz_0[i] = 2.0 * ta1_y_yzzz_xxzz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xxzz_1[i] * fe_0 + ta_yyzzz_xxzz_1[i] + ta1_y_yyzzz_xxzz_0[i] * pa_y[i] - ta1_y_yyzzz_xxzz_1[i] * pc_y[i];

        ta1_y_yyyzzz_xyyy_0[i] = 2.0 * ta1_y_yyyz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xyyy_1[i] * fe_0 + ta1_y_yyyzz_xyyy_0[i] * pa_z[i] - ta1_y_yyyzz_xyyy_1[i] * pc_z[i];

        ta1_y_yyyzzz_xyyz_0[i] = 2.0 * ta1_y_yyyz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xyyz_1[i] * fe_0 + ta1_y_yyyzz_xyy_0[i] * fe_0 - ta1_y_yyyzz_xyy_1[i] * fe_0 + ta1_y_yyyzz_xyyz_0[i] * pa_z[i] - ta1_y_yyyzz_xyyz_1[i] * pc_z[i];

        ta1_y_yyyzzz_xyzz_0[i] = 2.0 * ta1_y_yyyz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_xyz_1[i] * fe_0 + ta1_y_yyyzz_xyzz_0[i] * pa_z[i] - ta1_y_yyyzz_xyzz_1[i] * pc_z[i];

        ta1_y_yyyzzz_xzzz_0[i] = 2.0 * ta1_y_yzzz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xzzz_1[i] * fe_0 + ta_yyzzz_xzzz_1[i] + ta1_y_yyzzz_xzzz_0[i] * pa_y[i] - ta1_y_yyzzz_xzzz_1[i] * pc_y[i];

        ta1_y_yyyzzz_yyyy_0[i] = 2.0 * ta1_y_yyyz_yyyy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yyyy_1[i] * fe_0 + ta1_y_yyyzz_yyyy_0[i] * pa_z[i] - ta1_y_yyyzz_yyyy_1[i] * pc_z[i];

        ta1_y_yyyzzz_yyyz_0[i] = 2.0 * ta1_y_yyyz_yyyz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yyyz_1[i] * fe_0 + ta1_y_yyyzz_yyy_0[i] * fe_0 - ta1_y_yyyzz_yyy_1[i] * fe_0 + ta1_y_yyyzz_yyyz_0[i] * pa_z[i] - ta1_y_yyyzz_yyyz_1[i] * pc_z[i];

        ta1_y_yyyzzz_yyzz_0[i] = 2.0 * ta1_y_yyyz_yyzz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyyzz_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_yyz_1[i] * fe_0 + ta1_y_yyyzz_yyzz_0[i] * pa_z[i] - ta1_y_yyyzz_yyzz_1[i] * pc_z[i];

        ta1_y_yyyzzz_yzzz_0[i] = 2.0 * ta1_y_yyyz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yzzz_1[i] * fe_0 + 3.0 * ta1_y_yyyzz_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyyzz_yzz_1[i] * fe_0 + ta1_y_yyyzz_yzzz_0[i] * pa_z[i] - ta1_y_yyyzz_yzzz_1[i] * pc_z[i];

        ta1_y_yyyzzz_zzzz_0[i] = 2.0 * ta1_y_yzzz_zzzz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_zzzz_1[i] * fe_0 + ta_yyzzz_zzzz_1[i] + ta1_y_yyzzz_zzzz_0[i] * pa_y[i] - ta1_y_yyzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 795-810 components of targeted buffer : IG

    auto ta1_y_yyzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 795);

    auto ta1_y_yyzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 796);

    auto ta1_y_yyzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 797);

    auto ta1_y_yyzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 798);

    auto ta1_y_yyzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 799);

    auto ta1_y_yyzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 800);

    auto ta1_y_yyzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 801);

    auto ta1_y_yyzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 802);

    auto ta1_y_yyzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 803);

    auto ta1_y_yyzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 804);

    auto ta1_y_yyzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 805);

    auto ta1_y_yyzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 806);

    auto ta1_y_yyzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 807);

    auto ta1_y_yyzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 808);

    auto ta1_y_yyzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 809);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyzz_xxxx_0, ta1_y_yyzz_xxxx_1, ta1_y_yyzz_xxxy_0, ta1_y_yyzz_xxxy_1, ta1_y_yyzz_xxyy_0, ta1_y_yyzz_xxyy_1, ta1_y_yyzz_xxyz_0, ta1_y_yyzz_xxyz_1, ta1_y_yyzz_xyyy_0, ta1_y_yyzz_xyyy_1, ta1_y_yyzz_xyyz_0, ta1_y_yyzz_xyyz_1, ta1_y_yyzz_xyzz_0, ta1_y_yyzz_xyzz_1, ta1_y_yyzz_yyyy_0, ta1_y_yyzz_yyyy_1, ta1_y_yyzz_yyyz_0, ta1_y_yyzz_yyyz_1, ta1_y_yyzz_yyzz_0, ta1_y_yyzz_yyzz_1, ta1_y_yyzz_yzzz_0, ta1_y_yyzz_yzzz_1, ta1_y_yyzzz_xxxx_0, ta1_y_yyzzz_xxxx_1, ta1_y_yyzzz_xxxy_0, ta1_y_yyzzz_xxxy_1, ta1_y_yyzzz_xxy_0, ta1_y_yyzzz_xxy_1, ta1_y_yyzzz_xxyy_0, ta1_y_yyzzz_xxyy_1, ta1_y_yyzzz_xxyz_0, ta1_y_yyzzz_xxyz_1, ta1_y_yyzzz_xyy_0, ta1_y_yyzzz_xyy_1, ta1_y_yyzzz_xyyy_0, ta1_y_yyzzz_xyyy_1, ta1_y_yyzzz_xyyz_0, ta1_y_yyzzz_xyyz_1, ta1_y_yyzzz_xyz_0, ta1_y_yyzzz_xyz_1, ta1_y_yyzzz_xyzz_0, ta1_y_yyzzz_xyzz_1, ta1_y_yyzzz_yyy_0, ta1_y_yyzzz_yyy_1, ta1_y_yyzzz_yyyy_0, ta1_y_yyzzz_yyyy_1, ta1_y_yyzzz_yyyz_0, ta1_y_yyzzz_yyyz_1, ta1_y_yyzzz_yyz_0, ta1_y_yyzzz_yyz_1, ta1_y_yyzzz_yyzz_0, ta1_y_yyzzz_yyzz_1, ta1_y_yyzzz_yzz_0, ta1_y_yyzzz_yzz_1, ta1_y_yyzzz_yzzz_0, ta1_y_yyzzz_yzzz_1, ta1_y_yyzzzz_xxxx_0, ta1_y_yyzzzz_xxxy_0, ta1_y_yyzzzz_xxxz_0, ta1_y_yyzzzz_xxyy_0, ta1_y_yyzzzz_xxyz_0, ta1_y_yyzzzz_xxzz_0, ta1_y_yyzzzz_xyyy_0, ta1_y_yyzzzz_xyyz_0, ta1_y_yyzzzz_xyzz_0, ta1_y_yyzzzz_xzzz_0, ta1_y_yyzzzz_yyyy_0, ta1_y_yyzzzz_yyyz_0, ta1_y_yyzzzz_yyzz_0, ta1_y_yyzzzz_yzzz_0, ta1_y_yyzzzz_zzzz_0, ta1_y_yzzzz_xxxz_0, ta1_y_yzzzz_xxxz_1, ta1_y_yzzzz_xxzz_0, ta1_y_yzzzz_xxzz_1, ta1_y_yzzzz_xzzz_0, ta1_y_yzzzz_xzzz_1, ta1_y_yzzzz_zzzz_0, ta1_y_yzzzz_zzzz_1, ta1_y_zzzz_xxxz_0, ta1_y_zzzz_xxxz_1, ta1_y_zzzz_xxzz_0, ta1_y_zzzz_xxzz_1, ta1_y_zzzz_xzzz_0, ta1_y_zzzz_xzzz_1, ta1_y_zzzz_zzzz_0, ta1_y_zzzz_zzzz_1, ta_yzzzz_xxxz_1, ta_yzzzz_xxzz_1, ta_yzzzz_xzzz_1, ta_yzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzzz_xxxx_0[i] = 3.0 * ta1_y_yyzz_xxxx_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxxx_1[i] * fe_0 + ta1_y_yyzzz_xxxx_0[i] * pa_z[i] - ta1_y_yyzzz_xxxx_1[i] * pc_z[i];

        ta1_y_yyzzzz_xxxy_0[i] = 3.0 * ta1_y_yyzz_xxxy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxxy_1[i] * fe_0 + ta1_y_yyzzz_xxxy_0[i] * pa_z[i] - ta1_y_yyzzz_xxxy_1[i] * pc_z[i];

        ta1_y_yyzzzz_xxxz_0[i] = ta1_y_zzzz_xxxz_0[i] * fe_0 - ta1_y_zzzz_xxxz_1[i] * fe_0 + ta_yzzzz_xxxz_1[i] + ta1_y_yzzzz_xxxz_0[i] * pa_y[i] - ta1_y_yzzzz_xxxz_1[i] * pc_y[i];

        ta1_y_yyzzzz_xxyy_0[i] = 3.0 * ta1_y_yyzz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxyy_1[i] * fe_0 + ta1_y_yyzzz_xxyy_0[i] * pa_z[i] - ta1_y_yyzzz_xxyy_1[i] * pc_z[i];

        ta1_y_yyzzzz_xxyz_0[i] = 3.0 * ta1_y_yyzz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxyz_1[i] * fe_0 + ta1_y_yyzzz_xxy_0[i] * fe_0 - ta1_y_yyzzz_xxy_1[i] * fe_0 + ta1_y_yyzzz_xxyz_0[i] * pa_z[i] - ta1_y_yyzzz_xxyz_1[i] * pc_z[i];

        ta1_y_yyzzzz_xxzz_0[i] = ta1_y_zzzz_xxzz_0[i] * fe_0 - ta1_y_zzzz_xxzz_1[i] * fe_0 + ta_yzzzz_xxzz_1[i] + ta1_y_yzzzz_xxzz_0[i] * pa_y[i] - ta1_y_yzzzz_xxzz_1[i] * pc_y[i];

        ta1_y_yyzzzz_xyyy_0[i] = 3.0 * ta1_y_yyzz_xyyy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xyyy_1[i] * fe_0 + ta1_y_yyzzz_xyyy_0[i] * pa_z[i] - ta1_y_yyzzz_xyyy_1[i] * pc_z[i];

        ta1_y_yyzzzz_xyyz_0[i] = 3.0 * ta1_y_yyzz_xyyz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xyyz_1[i] * fe_0 + ta1_y_yyzzz_xyy_0[i] * fe_0 - ta1_y_yyzzz_xyy_1[i] * fe_0 + ta1_y_yyzzz_xyyz_0[i] * pa_z[i] - ta1_y_yyzzz_xyyz_1[i] * pc_z[i];

        ta1_y_yyzzzz_xyzz_0[i] = 3.0 * ta1_y_yyzz_xyzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xyzz_1[i] * fe_0 + 2.0 * ta1_y_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_xyz_1[i] * fe_0 + ta1_y_yyzzz_xyzz_0[i] * pa_z[i] - ta1_y_yyzzz_xyzz_1[i] * pc_z[i];

        ta1_y_yyzzzz_xzzz_0[i] = ta1_y_zzzz_xzzz_0[i] * fe_0 - ta1_y_zzzz_xzzz_1[i] * fe_0 + ta_yzzzz_xzzz_1[i] + ta1_y_yzzzz_xzzz_0[i] * pa_y[i] - ta1_y_yzzzz_xzzz_1[i] * pc_y[i];

        ta1_y_yyzzzz_yyyy_0[i] = 3.0 * ta1_y_yyzz_yyyy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yyyy_1[i] * fe_0 + ta1_y_yyzzz_yyyy_0[i] * pa_z[i] - ta1_y_yyzzz_yyyy_1[i] * pc_z[i];

        ta1_y_yyzzzz_yyyz_0[i] = 3.0 * ta1_y_yyzz_yyyz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yyyz_1[i] * fe_0 + ta1_y_yyzzz_yyy_0[i] * fe_0 - ta1_y_yyzzz_yyy_1[i] * fe_0 + ta1_y_yyzzz_yyyz_0[i] * pa_z[i] - ta1_y_yyzzz_yyyz_1[i] * pc_z[i];

        ta1_y_yyzzzz_yyzz_0[i] = 3.0 * ta1_y_yyzz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yyzz_1[i] * fe_0 + 2.0 * ta1_y_yyzzz_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_yyz_1[i] * fe_0 + ta1_y_yyzzz_yyzz_0[i] * pa_z[i] - ta1_y_yyzzz_yyzz_1[i] * pc_z[i];

        ta1_y_yyzzzz_yzzz_0[i] = 3.0 * ta1_y_yyzz_yzzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yzzz_1[i] * fe_0 + 3.0 * ta1_y_yyzzz_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyzzz_yzz_1[i] * fe_0 + ta1_y_yyzzz_yzzz_0[i] * pa_z[i] - ta1_y_yyzzz_yzzz_1[i] * pc_z[i];

        ta1_y_yyzzzz_zzzz_0[i] = ta1_y_zzzz_zzzz_0[i] * fe_0 - ta1_y_zzzz_zzzz_1[i] * fe_0 + ta_yzzzz_zzzz_1[i] + ta1_y_yzzzz_zzzz_0[i] * pa_y[i] - ta1_y_yzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 810-825 components of targeted buffer : IG

    auto ta1_y_yzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 810);

    auto ta1_y_yzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 811);

    auto ta1_y_yzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 812);

    auto ta1_y_yzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 813);

    auto ta1_y_yzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 814);

    auto ta1_y_yzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 815);

    auto ta1_y_yzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 816);

    auto ta1_y_yzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 817);

    auto ta1_y_yzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 818);

    auto ta1_y_yzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 819);

    auto ta1_y_yzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 820);

    auto ta1_y_yzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 821);

    auto ta1_y_yzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 822);

    auto ta1_y_yzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 823);

    auto ta1_y_yzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 824);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yzzz_xxxy_0, ta1_y_yzzz_xxxy_1, ta1_y_yzzz_xxyy_0, ta1_y_yzzz_xxyy_1, ta1_y_yzzz_xyyy_0, ta1_y_yzzz_xyyy_1, ta1_y_yzzz_yyyy_0, ta1_y_yzzz_yyyy_1, ta1_y_yzzzz_xxxy_0, ta1_y_yzzzz_xxxy_1, ta1_y_yzzzz_xxyy_0, ta1_y_yzzzz_xxyy_1, ta1_y_yzzzz_xyyy_0, ta1_y_yzzzz_xyyy_1, ta1_y_yzzzz_yyyy_0, ta1_y_yzzzz_yyyy_1, ta1_y_yzzzzz_xxxx_0, ta1_y_yzzzzz_xxxy_0, ta1_y_yzzzzz_xxxz_0, ta1_y_yzzzzz_xxyy_0, ta1_y_yzzzzz_xxyz_0, ta1_y_yzzzzz_xxzz_0, ta1_y_yzzzzz_xyyy_0, ta1_y_yzzzzz_xyyz_0, ta1_y_yzzzzz_xyzz_0, ta1_y_yzzzzz_xzzz_0, ta1_y_yzzzzz_yyyy_0, ta1_y_yzzzzz_yyyz_0, ta1_y_yzzzzz_yyzz_0, ta1_y_yzzzzz_yzzz_0, ta1_y_yzzzzz_zzzz_0, ta1_y_zzzzz_xxxx_0, ta1_y_zzzzz_xxxx_1, ta1_y_zzzzz_xxxz_0, ta1_y_zzzzz_xxxz_1, ta1_y_zzzzz_xxyz_0, ta1_y_zzzzz_xxyz_1, ta1_y_zzzzz_xxz_0, ta1_y_zzzzz_xxz_1, ta1_y_zzzzz_xxzz_0, ta1_y_zzzzz_xxzz_1, ta1_y_zzzzz_xyyz_0, ta1_y_zzzzz_xyyz_1, ta1_y_zzzzz_xyz_0, ta1_y_zzzzz_xyz_1, ta1_y_zzzzz_xyzz_0, ta1_y_zzzzz_xyzz_1, ta1_y_zzzzz_xzz_0, ta1_y_zzzzz_xzz_1, ta1_y_zzzzz_xzzz_0, ta1_y_zzzzz_xzzz_1, ta1_y_zzzzz_yyyz_0, ta1_y_zzzzz_yyyz_1, ta1_y_zzzzz_yyz_0, ta1_y_zzzzz_yyz_1, ta1_y_zzzzz_yyzz_0, ta1_y_zzzzz_yyzz_1, ta1_y_zzzzz_yzz_0, ta1_y_zzzzz_yzz_1, ta1_y_zzzzz_yzzz_0, ta1_y_zzzzz_yzzz_1, ta1_y_zzzzz_zzz_0, ta1_y_zzzzz_zzz_1, ta1_y_zzzzz_zzzz_0, ta1_y_zzzzz_zzzz_1, ta_zzzzz_xxxx_1, ta_zzzzz_xxxz_1, ta_zzzzz_xxyz_1, ta_zzzzz_xxzz_1, ta_zzzzz_xyyz_1, ta_zzzzz_xyzz_1, ta_zzzzz_xzzz_1, ta_zzzzz_yyyz_1, ta_zzzzz_yyzz_1, ta_zzzzz_yzzz_1, ta_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzzz_xxxx_0[i] = ta_zzzzz_xxxx_1[i] + ta1_y_zzzzz_xxxx_0[i] * pa_y[i] - ta1_y_zzzzz_xxxx_1[i] * pc_y[i];

        ta1_y_yzzzzz_xxxy_0[i] = 4.0 * ta1_y_yzzz_xxxy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_xxxy_1[i] * fe_0 + ta1_y_yzzzz_xxxy_0[i] * pa_z[i] - ta1_y_yzzzz_xxxy_1[i] * pc_z[i];

        ta1_y_yzzzzz_xxxz_0[i] = ta_zzzzz_xxxz_1[i] + ta1_y_zzzzz_xxxz_0[i] * pa_y[i] - ta1_y_zzzzz_xxxz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xxyy_0[i] = 4.0 * ta1_y_yzzz_xxyy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_xxyy_1[i] * fe_0 + ta1_y_yzzzz_xxyy_0[i] * pa_z[i] - ta1_y_yzzzz_xxyy_1[i] * pc_z[i];

        ta1_y_yzzzzz_xxyz_0[i] = ta1_y_zzzzz_xxz_0[i] * fe_0 - ta1_y_zzzzz_xxz_1[i] * fe_0 + ta_zzzzz_xxyz_1[i] + ta1_y_zzzzz_xxyz_0[i] * pa_y[i] - ta1_y_zzzzz_xxyz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xxzz_0[i] = ta_zzzzz_xxzz_1[i] + ta1_y_zzzzz_xxzz_0[i] * pa_y[i] - ta1_y_zzzzz_xxzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xyyy_0[i] = 4.0 * ta1_y_yzzz_xyyy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_xyyy_1[i] * fe_0 + ta1_y_yzzzz_xyyy_0[i] * pa_z[i] - ta1_y_yzzzz_xyyy_1[i] * pc_z[i];

        ta1_y_yzzzzz_xyyz_0[i] = 2.0 * ta1_y_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xyz_1[i] * fe_0 + ta_zzzzz_xyyz_1[i] + ta1_y_zzzzz_xyyz_0[i] * pa_y[i] - ta1_y_zzzzz_xyyz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xyzz_0[i] = ta1_y_zzzzz_xzz_0[i] * fe_0 - ta1_y_zzzzz_xzz_1[i] * fe_0 + ta_zzzzz_xyzz_1[i] + ta1_y_zzzzz_xyzz_0[i] * pa_y[i] - ta1_y_zzzzz_xyzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xzzz_0[i] = ta_zzzzz_xzzz_1[i] + ta1_y_zzzzz_xzzz_0[i] * pa_y[i] - ta1_y_zzzzz_xzzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_yyyy_0[i] = 4.0 * ta1_y_yzzz_yyyy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_yyyy_1[i] * fe_0 + ta1_y_yzzzz_yyyy_0[i] * pa_z[i] - ta1_y_yzzzz_yyyy_1[i] * pc_z[i];

        ta1_y_yzzzzz_yyyz_0[i] = 3.0 * ta1_y_zzzzz_yyz_0[i] * fe_0 - 3.0 * ta1_y_zzzzz_yyz_1[i] * fe_0 + ta_zzzzz_yyyz_1[i] + ta1_y_zzzzz_yyyz_0[i] * pa_y[i] - ta1_y_zzzzz_yyyz_1[i] * pc_y[i];

        ta1_y_yzzzzz_yyzz_0[i] = 2.0 * ta1_y_zzzzz_yzz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_yzz_1[i] * fe_0 + ta_zzzzz_yyzz_1[i] + ta1_y_zzzzz_yyzz_0[i] * pa_y[i] - ta1_y_zzzzz_yyzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_yzzz_0[i] = ta1_y_zzzzz_zzz_0[i] * fe_0 - ta1_y_zzzzz_zzz_1[i] * fe_0 + ta_zzzzz_yzzz_1[i] + ta1_y_zzzzz_yzzz_0[i] * pa_y[i] - ta1_y_zzzzz_yzzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_zzzz_0[i] = ta_zzzzz_zzzz_1[i] + ta1_y_zzzzz_zzzz_0[i] * pa_y[i] - ta1_y_zzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 825-840 components of targeted buffer : IG

    auto ta1_y_zzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 825);

    auto ta1_y_zzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 826);

    auto ta1_y_zzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 827);

    auto ta1_y_zzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 828);

    auto ta1_y_zzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 829);

    auto ta1_y_zzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 830);

    auto ta1_y_zzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 831);

    auto ta1_y_zzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 832);

    auto ta1_y_zzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 833);

    auto ta1_y_zzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 834);

    auto ta1_y_zzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 835);

    auto ta1_y_zzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 836);

    auto ta1_y_zzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 837);

    auto ta1_y_zzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 838);

    auto ta1_y_zzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 839);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zzzz_xxxx_0, ta1_y_zzzz_xxxx_1, ta1_y_zzzz_xxxy_0, ta1_y_zzzz_xxxy_1, ta1_y_zzzz_xxxz_0, ta1_y_zzzz_xxxz_1, ta1_y_zzzz_xxyy_0, ta1_y_zzzz_xxyy_1, ta1_y_zzzz_xxyz_0, ta1_y_zzzz_xxyz_1, ta1_y_zzzz_xxzz_0, ta1_y_zzzz_xxzz_1, ta1_y_zzzz_xyyy_0, ta1_y_zzzz_xyyy_1, ta1_y_zzzz_xyyz_0, ta1_y_zzzz_xyyz_1, ta1_y_zzzz_xyzz_0, ta1_y_zzzz_xyzz_1, ta1_y_zzzz_xzzz_0, ta1_y_zzzz_xzzz_1, ta1_y_zzzz_yyyy_0, ta1_y_zzzz_yyyy_1, ta1_y_zzzz_yyyz_0, ta1_y_zzzz_yyyz_1, ta1_y_zzzz_yyzz_0, ta1_y_zzzz_yyzz_1, ta1_y_zzzz_yzzz_0, ta1_y_zzzz_yzzz_1, ta1_y_zzzz_zzzz_0, ta1_y_zzzz_zzzz_1, ta1_y_zzzzz_xxx_0, ta1_y_zzzzz_xxx_1, ta1_y_zzzzz_xxxx_0, ta1_y_zzzzz_xxxx_1, ta1_y_zzzzz_xxxy_0, ta1_y_zzzzz_xxxy_1, ta1_y_zzzzz_xxxz_0, ta1_y_zzzzz_xxxz_1, ta1_y_zzzzz_xxy_0, ta1_y_zzzzz_xxy_1, ta1_y_zzzzz_xxyy_0, ta1_y_zzzzz_xxyy_1, ta1_y_zzzzz_xxyz_0, ta1_y_zzzzz_xxyz_1, ta1_y_zzzzz_xxz_0, ta1_y_zzzzz_xxz_1, ta1_y_zzzzz_xxzz_0, ta1_y_zzzzz_xxzz_1, ta1_y_zzzzz_xyy_0, ta1_y_zzzzz_xyy_1, ta1_y_zzzzz_xyyy_0, ta1_y_zzzzz_xyyy_1, ta1_y_zzzzz_xyyz_0, ta1_y_zzzzz_xyyz_1, ta1_y_zzzzz_xyz_0, ta1_y_zzzzz_xyz_1, ta1_y_zzzzz_xyzz_0, ta1_y_zzzzz_xyzz_1, ta1_y_zzzzz_xzz_0, ta1_y_zzzzz_xzz_1, ta1_y_zzzzz_xzzz_0, ta1_y_zzzzz_xzzz_1, ta1_y_zzzzz_yyy_0, ta1_y_zzzzz_yyy_1, ta1_y_zzzzz_yyyy_0, ta1_y_zzzzz_yyyy_1, ta1_y_zzzzz_yyyz_0, ta1_y_zzzzz_yyyz_1, ta1_y_zzzzz_yyz_0, ta1_y_zzzzz_yyz_1, ta1_y_zzzzz_yyzz_0, ta1_y_zzzzz_yyzz_1, ta1_y_zzzzz_yzz_0, ta1_y_zzzzz_yzz_1, ta1_y_zzzzz_yzzz_0, ta1_y_zzzzz_yzzz_1, ta1_y_zzzzz_zzz_0, ta1_y_zzzzz_zzz_1, ta1_y_zzzzz_zzzz_0, ta1_y_zzzzz_zzzz_1, ta1_y_zzzzzz_xxxx_0, ta1_y_zzzzzz_xxxy_0, ta1_y_zzzzzz_xxxz_0, ta1_y_zzzzzz_xxyy_0, ta1_y_zzzzzz_xxyz_0, ta1_y_zzzzzz_xxzz_0, ta1_y_zzzzzz_xyyy_0, ta1_y_zzzzzz_xyyz_0, ta1_y_zzzzzz_xyzz_0, ta1_y_zzzzzz_xzzz_0, ta1_y_zzzzzz_yyyy_0, ta1_y_zzzzzz_yyyz_0, ta1_y_zzzzzz_yyzz_0, ta1_y_zzzzzz_yzzz_0, ta1_y_zzzzzz_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzzz_xxxx_0[i] = 5.0 * ta1_y_zzzz_xxxx_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxxx_1[i] * fe_0 + ta1_y_zzzzz_xxxx_0[i] * pa_z[i] - ta1_y_zzzzz_xxxx_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxxy_0[i] = 5.0 * ta1_y_zzzz_xxxy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxxy_1[i] * fe_0 + ta1_y_zzzzz_xxxy_0[i] * pa_z[i] - ta1_y_zzzzz_xxxy_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxxz_0[i] = 5.0 * ta1_y_zzzz_xxxz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxxz_1[i] * fe_0 + ta1_y_zzzzz_xxx_0[i] * fe_0 - ta1_y_zzzzz_xxx_1[i] * fe_0 + ta1_y_zzzzz_xxxz_0[i] * pa_z[i] - ta1_y_zzzzz_xxxz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxyy_0[i] = 5.0 * ta1_y_zzzz_xxyy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxyy_1[i] * fe_0 + ta1_y_zzzzz_xxyy_0[i] * pa_z[i] - ta1_y_zzzzz_xxyy_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxyz_0[i] = 5.0 * ta1_y_zzzz_xxyz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxyz_1[i] * fe_0 + ta1_y_zzzzz_xxy_0[i] * fe_0 - ta1_y_zzzzz_xxy_1[i] * fe_0 + ta1_y_zzzzz_xxyz_0[i] * pa_z[i] - ta1_y_zzzzz_xxyz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxzz_0[i] = 5.0 * ta1_y_zzzz_xxzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_zzzzz_xxz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xxz_1[i] * fe_0 + ta1_y_zzzzz_xxzz_0[i] * pa_z[i] - ta1_y_zzzzz_xxzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xyyy_0[i] = 5.0 * ta1_y_zzzz_xyyy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xyyy_1[i] * fe_0 + ta1_y_zzzzz_xyyy_0[i] * pa_z[i] - ta1_y_zzzzz_xyyy_1[i] * pc_z[i];

        ta1_y_zzzzzz_xyyz_0[i] = 5.0 * ta1_y_zzzz_xyyz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xyyz_1[i] * fe_0 + ta1_y_zzzzz_xyy_0[i] * fe_0 - ta1_y_zzzzz_xyy_1[i] * fe_0 + ta1_y_zzzzz_xyyz_0[i] * pa_z[i] - ta1_y_zzzzz_xyyz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xyzz_0[i] = 5.0 * ta1_y_zzzz_xyzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xyzz_1[i] * fe_0 + 2.0 * ta1_y_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xyz_1[i] * fe_0 + ta1_y_zzzzz_xyzz_0[i] * pa_z[i] - ta1_y_zzzzz_xyzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xzzz_0[i] = 5.0 * ta1_y_zzzz_xzzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xzzz_1[i] * fe_0 + 3.0 * ta1_y_zzzzz_xzz_0[i] * fe_0 - 3.0 * ta1_y_zzzzz_xzz_1[i] * fe_0 + ta1_y_zzzzz_xzzz_0[i] * pa_z[i] - ta1_y_zzzzz_xzzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_yyyy_0[i] = 5.0 * ta1_y_zzzz_yyyy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yyyy_1[i] * fe_0 + ta1_y_zzzzz_yyyy_0[i] * pa_z[i] - ta1_y_zzzzz_yyyy_1[i] * pc_z[i];

        ta1_y_zzzzzz_yyyz_0[i] = 5.0 * ta1_y_zzzz_yyyz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yyyz_1[i] * fe_0 + ta1_y_zzzzz_yyy_0[i] * fe_0 - ta1_y_zzzzz_yyy_1[i] * fe_0 + ta1_y_zzzzz_yyyz_0[i] * pa_z[i] - ta1_y_zzzzz_yyyz_1[i] * pc_z[i];

        ta1_y_zzzzzz_yyzz_0[i] = 5.0 * ta1_y_zzzz_yyzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_y_zzzzz_yyz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_yyz_1[i] * fe_0 + ta1_y_zzzzz_yyzz_0[i] * pa_z[i] - ta1_y_zzzzz_yyzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_yzzz_0[i] = 5.0 * ta1_y_zzzz_yzzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yzzz_1[i] * fe_0 + 3.0 * ta1_y_zzzzz_yzz_0[i] * fe_0 - 3.0 * ta1_y_zzzzz_yzz_1[i] * fe_0 + ta1_y_zzzzz_yzzz_0[i] * pa_z[i] - ta1_y_zzzzz_yzzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_zzzz_0[i] = 5.0 * ta1_y_zzzz_zzzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_zzzz_1[i] * fe_0 + 4.0 * ta1_y_zzzzz_zzz_0[i] * fe_0 - 4.0 * ta1_y_zzzzz_zzz_1[i] * fe_0 + ta1_y_zzzzz_zzzz_0[i] * pa_z[i] - ta1_y_zzzzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 840-855 components of targeted buffer : IG

    auto ta1_z_xxxxxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 840);

    auto ta1_z_xxxxxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 841);

    auto ta1_z_xxxxxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 842);

    auto ta1_z_xxxxxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 843);

    auto ta1_z_xxxxxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 844);

    auto ta1_z_xxxxxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 845);

    auto ta1_z_xxxxxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 846);

    auto ta1_z_xxxxxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 847);

    auto ta1_z_xxxxxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 848);

    auto ta1_z_xxxxxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 849);

    auto ta1_z_xxxxxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 850);

    auto ta1_z_xxxxxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 851);

    auto ta1_z_xxxxxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 852);

    auto ta1_z_xxxxxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 853);

    auto ta1_z_xxxxxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 854);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xxxx_xxxx_0, ta1_z_xxxx_xxxx_1, ta1_z_xxxx_xxxy_0, ta1_z_xxxx_xxxy_1, ta1_z_xxxx_xxxz_0, ta1_z_xxxx_xxxz_1, ta1_z_xxxx_xxyy_0, ta1_z_xxxx_xxyy_1, ta1_z_xxxx_xxyz_0, ta1_z_xxxx_xxyz_1, ta1_z_xxxx_xxzz_0, ta1_z_xxxx_xxzz_1, ta1_z_xxxx_xyyy_0, ta1_z_xxxx_xyyy_1, ta1_z_xxxx_xyyz_0, ta1_z_xxxx_xyyz_1, ta1_z_xxxx_xyzz_0, ta1_z_xxxx_xyzz_1, ta1_z_xxxx_xzzz_0, ta1_z_xxxx_xzzz_1, ta1_z_xxxx_yyyy_0, ta1_z_xxxx_yyyy_1, ta1_z_xxxx_yyyz_0, ta1_z_xxxx_yyyz_1, ta1_z_xxxx_yyzz_0, ta1_z_xxxx_yyzz_1, ta1_z_xxxx_yzzz_0, ta1_z_xxxx_yzzz_1, ta1_z_xxxx_zzzz_0, ta1_z_xxxx_zzzz_1, ta1_z_xxxxx_xxx_0, ta1_z_xxxxx_xxx_1, ta1_z_xxxxx_xxxx_0, ta1_z_xxxxx_xxxx_1, ta1_z_xxxxx_xxxy_0, ta1_z_xxxxx_xxxy_1, ta1_z_xxxxx_xxxz_0, ta1_z_xxxxx_xxxz_1, ta1_z_xxxxx_xxy_0, ta1_z_xxxxx_xxy_1, ta1_z_xxxxx_xxyy_0, ta1_z_xxxxx_xxyy_1, ta1_z_xxxxx_xxyz_0, ta1_z_xxxxx_xxyz_1, ta1_z_xxxxx_xxz_0, ta1_z_xxxxx_xxz_1, ta1_z_xxxxx_xxzz_0, ta1_z_xxxxx_xxzz_1, ta1_z_xxxxx_xyy_0, ta1_z_xxxxx_xyy_1, ta1_z_xxxxx_xyyy_0, ta1_z_xxxxx_xyyy_1, ta1_z_xxxxx_xyyz_0, ta1_z_xxxxx_xyyz_1, ta1_z_xxxxx_xyz_0, ta1_z_xxxxx_xyz_1, ta1_z_xxxxx_xyzz_0, ta1_z_xxxxx_xyzz_1, ta1_z_xxxxx_xzz_0, ta1_z_xxxxx_xzz_1, ta1_z_xxxxx_xzzz_0, ta1_z_xxxxx_xzzz_1, ta1_z_xxxxx_yyy_0, ta1_z_xxxxx_yyy_1, ta1_z_xxxxx_yyyy_0, ta1_z_xxxxx_yyyy_1, ta1_z_xxxxx_yyyz_0, ta1_z_xxxxx_yyyz_1, ta1_z_xxxxx_yyz_0, ta1_z_xxxxx_yyz_1, ta1_z_xxxxx_yyzz_0, ta1_z_xxxxx_yyzz_1, ta1_z_xxxxx_yzz_0, ta1_z_xxxxx_yzz_1, ta1_z_xxxxx_yzzz_0, ta1_z_xxxxx_yzzz_1, ta1_z_xxxxx_zzz_0, ta1_z_xxxxx_zzz_1, ta1_z_xxxxx_zzzz_0, ta1_z_xxxxx_zzzz_1, ta1_z_xxxxxx_xxxx_0, ta1_z_xxxxxx_xxxy_0, ta1_z_xxxxxx_xxxz_0, ta1_z_xxxxxx_xxyy_0, ta1_z_xxxxxx_xxyz_0, ta1_z_xxxxxx_xxzz_0, ta1_z_xxxxxx_xyyy_0, ta1_z_xxxxxx_xyyz_0, ta1_z_xxxxxx_xyzz_0, ta1_z_xxxxxx_xzzz_0, ta1_z_xxxxxx_yyyy_0, ta1_z_xxxxxx_yyyz_0, ta1_z_xxxxxx_yyzz_0, ta1_z_xxxxxx_yzzz_0, ta1_z_xxxxxx_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxx_xxxx_0[i] = 5.0 * ta1_z_xxxx_xxxx_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxxx_1[i] * fe_0 + 4.0 * ta1_z_xxxxx_xxx_0[i] * fe_0 - 4.0 * ta1_z_xxxxx_xxx_1[i] * fe_0 + ta1_z_xxxxx_xxxx_0[i] * pa_x[i] - ta1_z_xxxxx_xxxx_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxxy_0[i] = 5.0 * ta1_z_xxxx_xxxy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xxxxx_xxy_0[i] * fe_0 - 3.0 * ta1_z_xxxxx_xxy_1[i] * fe_0 + ta1_z_xxxxx_xxxy_0[i] * pa_x[i] - ta1_z_xxxxx_xxxy_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxxz_0[i] = 5.0 * ta1_z_xxxx_xxxz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xxxxx_xxz_0[i] * fe_0 - 3.0 * ta1_z_xxxxx_xxz_1[i] * fe_0 + ta1_z_xxxxx_xxxz_0[i] * pa_x[i] - ta1_z_xxxxx_xxxz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxyy_0[i] = 5.0 * ta1_z_xxxx_xxyy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxxxx_xyy_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xyy_1[i] * fe_0 + ta1_z_xxxxx_xxyy_0[i] * pa_x[i] - ta1_z_xxxxx_xxyy_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxyz_0[i] = 5.0 * ta1_z_xxxx_xxyz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xyz_1[i] * fe_0 + ta1_z_xxxxx_xxyz_0[i] * pa_x[i] - ta1_z_xxxxx_xxyz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxzz_0[i] = 5.0 * ta1_z_xxxx_xxzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xxxxx_xzz_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xzz_1[i] * fe_0 + ta1_z_xxxxx_xxzz_0[i] * pa_x[i] - ta1_z_xxxxx_xxzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xyyy_0[i] = 5.0 * ta1_z_xxxx_xyyy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xyyy_1[i] * fe_0 + ta1_z_xxxxx_yyy_0[i] * fe_0 - ta1_z_xxxxx_yyy_1[i] * fe_0 + ta1_z_xxxxx_xyyy_0[i] * pa_x[i] - ta1_z_xxxxx_xyyy_1[i] * pc_x[i];

        ta1_z_xxxxxx_xyyz_0[i] = 5.0 * ta1_z_xxxx_xyyz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xyyz_1[i] * fe_0 + ta1_z_xxxxx_yyz_0[i] * fe_0 - ta1_z_xxxxx_yyz_1[i] * fe_0 + ta1_z_xxxxx_xyyz_0[i] * pa_x[i] - ta1_z_xxxxx_xyyz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xyzz_0[i] = 5.0 * ta1_z_xxxx_xyzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xyzz_1[i] * fe_0 + ta1_z_xxxxx_yzz_0[i] * fe_0 - ta1_z_xxxxx_yzz_1[i] * fe_0 + ta1_z_xxxxx_xyzz_0[i] * pa_x[i] - ta1_z_xxxxx_xyzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xzzz_0[i] = 5.0 * ta1_z_xxxx_xzzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xzzz_1[i] * fe_0 + ta1_z_xxxxx_zzz_0[i] * fe_0 - ta1_z_xxxxx_zzz_1[i] * fe_0 + ta1_z_xxxxx_xzzz_0[i] * pa_x[i] - ta1_z_xxxxx_xzzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_yyyy_0[i] = 5.0 * ta1_z_xxxx_yyyy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yyyy_1[i] * fe_0 + ta1_z_xxxxx_yyyy_0[i] * pa_x[i] - ta1_z_xxxxx_yyyy_1[i] * pc_x[i];

        ta1_z_xxxxxx_yyyz_0[i] = 5.0 * ta1_z_xxxx_yyyz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yyyz_1[i] * fe_0 + ta1_z_xxxxx_yyyz_0[i] * pa_x[i] - ta1_z_xxxxx_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxxx_yyzz_0[i] = 5.0 * ta1_z_xxxx_yyzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yyzz_1[i] * fe_0 + ta1_z_xxxxx_yyzz_0[i] * pa_x[i] - ta1_z_xxxxx_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_yzzz_0[i] = 5.0 * ta1_z_xxxx_yzzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yzzz_1[i] * fe_0 + ta1_z_xxxxx_yzzz_0[i] * pa_x[i] - ta1_z_xxxxx_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_zzzz_0[i] = 5.0 * ta1_z_xxxx_zzzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_zzzz_1[i] * fe_0 + ta1_z_xxxxx_zzzz_0[i] * pa_x[i] - ta1_z_xxxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 855-870 components of targeted buffer : IG

    auto ta1_z_xxxxxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 855);

    auto ta1_z_xxxxxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 856);

    auto ta1_z_xxxxxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 857);

    auto ta1_z_xxxxxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 858);

    auto ta1_z_xxxxxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 859);

    auto ta1_z_xxxxxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 860);

    auto ta1_z_xxxxxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 861);

    auto ta1_z_xxxxxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 862);

    auto ta1_z_xxxxxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 863);

    auto ta1_z_xxxxxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 864);

    auto ta1_z_xxxxxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 865);

    auto ta1_z_xxxxxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 866);

    auto ta1_z_xxxxxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 867);

    auto ta1_z_xxxxxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 868);

    auto ta1_z_xxxxxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 869);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxxx_xxx_0, ta1_z_xxxxx_xxx_1, ta1_z_xxxxx_xxxx_0, ta1_z_xxxxx_xxxx_1, ta1_z_xxxxx_xxxy_0, ta1_z_xxxxx_xxxy_1, ta1_z_xxxxx_xxxz_0, ta1_z_xxxxx_xxxz_1, ta1_z_xxxxx_xxy_0, ta1_z_xxxxx_xxy_1, ta1_z_xxxxx_xxyy_0, ta1_z_xxxxx_xxyy_1, ta1_z_xxxxx_xxyz_0, ta1_z_xxxxx_xxyz_1, ta1_z_xxxxx_xxz_0, ta1_z_xxxxx_xxz_1, ta1_z_xxxxx_xxzz_0, ta1_z_xxxxx_xxzz_1, ta1_z_xxxxx_xyy_0, ta1_z_xxxxx_xyy_1, ta1_z_xxxxx_xyyy_0, ta1_z_xxxxx_xyyy_1, ta1_z_xxxxx_xyyz_0, ta1_z_xxxxx_xyyz_1, ta1_z_xxxxx_xyz_0, ta1_z_xxxxx_xyz_1, ta1_z_xxxxx_xyzz_0, ta1_z_xxxxx_xyzz_1, ta1_z_xxxxx_xzz_0, ta1_z_xxxxx_xzz_1, ta1_z_xxxxx_xzzz_0, ta1_z_xxxxx_xzzz_1, ta1_z_xxxxx_zzzz_0, ta1_z_xxxxx_zzzz_1, ta1_z_xxxxxy_xxxx_0, ta1_z_xxxxxy_xxxy_0, ta1_z_xxxxxy_xxxz_0, ta1_z_xxxxxy_xxyy_0, ta1_z_xxxxxy_xxyz_0, ta1_z_xxxxxy_xxzz_0, ta1_z_xxxxxy_xyyy_0, ta1_z_xxxxxy_xyyz_0, ta1_z_xxxxxy_xyzz_0, ta1_z_xxxxxy_xzzz_0, ta1_z_xxxxxy_yyyy_0, ta1_z_xxxxxy_yyyz_0, ta1_z_xxxxxy_yyzz_0, ta1_z_xxxxxy_yzzz_0, ta1_z_xxxxxy_zzzz_0, ta1_z_xxxxy_yyyy_0, ta1_z_xxxxy_yyyy_1, ta1_z_xxxxy_yyyz_0, ta1_z_xxxxy_yyyz_1, ta1_z_xxxxy_yyzz_0, ta1_z_xxxxy_yyzz_1, ta1_z_xxxxy_yzzz_0, ta1_z_xxxxy_yzzz_1, ta1_z_xxxy_yyyy_0, ta1_z_xxxy_yyyy_1, ta1_z_xxxy_yyyz_0, ta1_z_xxxy_yyyz_1, ta1_z_xxxy_yyzz_0, ta1_z_xxxy_yyzz_1, ta1_z_xxxy_yzzz_0, ta1_z_xxxy_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxy_xxxx_0[i] = ta1_z_xxxxx_xxxx_0[i] * pa_y[i] - ta1_z_xxxxx_xxxx_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxxy_0[i] = ta1_z_xxxxx_xxx_0[i] * fe_0 - ta1_z_xxxxx_xxx_1[i] * fe_0 + ta1_z_xxxxx_xxxy_0[i] * pa_y[i] - ta1_z_xxxxx_xxxy_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxxz_0[i] = ta1_z_xxxxx_xxxz_0[i] * pa_y[i] - ta1_z_xxxxx_xxxz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxyy_0[i] = 2.0 * ta1_z_xxxxx_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xxy_1[i] * fe_0 + ta1_z_xxxxx_xxyy_0[i] * pa_y[i] - ta1_z_xxxxx_xxyy_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxyz_0[i] = ta1_z_xxxxx_xxz_0[i] * fe_0 - ta1_z_xxxxx_xxz_1[i] * fe_0 + ta1_z_xxxxx_xxyz_0[i] * pa_y[i] - ta1_z_xxxxx_xxyz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxzz_0[i] = ta1_z_xxxxx_xxzz_0[i] * pa_y[i] - ta1_z_xxxxx_xxzz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xyyy_0[i] = 3.0 * ta1_z_xxxxx_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxxxx_xyy_1[i] * fe_0 + ta1_z_xxxxx_xyyy_0[i] * pa_y[i] - ta1_z_xxxxx_xyyy_1[i] * pc_y[i];

        ta1_z_xxxxxy_xyyz_0[i] = 2.0 * ta1_z_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xyz_1[i] * fe_0 + ta1_z_xxxxx_xyyz_0[i] * pa_y[i] - ta1_z_xxxxx_xyyz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xyzz_0[i] = ta1_z_xxxxx_xzz_0[i] * fe_0 - ta1_z_xxxxx_xzz_1[i] * fe_0 + ta1_z_xxxxx_xyzz_0[i] * pa_y[i] - ta1_z_xxxxx_xyzz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xzzz_0[i] = ta1_z_xxxxx_xzzz_0[i] * pa_y[i] - ta1_z_xxxxx_xzzz_1[i] * pc_y[i];

        ta1_z_xxxxxy_yyyy_0[i] = 4.0 * ta1_z_xxxy_yyyy_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yyyy_1[i] * fe_0 + ta1_z_xxxxy_yyyy_0[i] * pa_x[i] - ta1_z_xxxxy_yyyy_1[i] * pc_x[i];

        ta1_z_xxxxxy_yyyz_0[i] = 4.0 * ta1_z_xxxy_yyyz_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yyyz_1[i] * fe_0 + ta1_z_xxxxy_yyyz_0[i] * pa_x[i] - ta1_z_xxxxy_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxxy_yyzz_0[i] = 4.0 * ta1_z_xxxy_yyzz_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yyzz_1[i] * fe_0 + ta1_z_xxxxy_yyzz_0[i] * pa_x[i] - ta1_z_xxxxy_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxxy_yzzz_0[i] = 4.0 * ta1_z_xxxy_yzzz_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yzzz_1[i] * fe_0 + ta1_z_xxxxy_yzzz_0[i] * pa_x[i] - ta1_z_xxxxy_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxxy_zzzz_0[i] = ta1_z_xxxxx_zzzz_0[i] * pa_y[i] - ta1_z_xxxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 870-885 components of targeted buffer : IG

    auto ta1_z_xxxxxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 870);

    auto ta1_z_xxxxxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 871);

    auto ta1_z_xxxxxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 872);

    auto ta1_z_xxxxxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 873);

    auto ta1_z_xxxxxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 874);

    auto ta1_z_xxxxxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 875);

    auto ta1_z_xxxxxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 876);

    auto ta1_z_xxxxxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 877);

    auto ta1_z_xxxxxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 878);

    auto ta1_z_xxxxxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 879);

    auto ta1_z_xxxxxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 880);

    auto ta1_z_xxxxxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 881);

    auto ta1_z_xxxxxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 882);

    auto ta1_z_xxxxxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 883);

    auto ta1_z_xxxxxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 884);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxxx_xxx_0, ta1_z_xxxxx_xxx_1, ta1_z_xxxxx_xxxx_0, ta1_z_xxxxx_xxxx_1, ta1_z_xxxxx_xxxy_0, ta1_z_xxxxx_xxxy_1, ta1_z_xxxxx_xxxz_0, ta1_z_xxxxx_xxxz_1, ta1_z_xxxxx_xxy_0, ta1_z_xxxxx_xxy_1, ta1_z_xxxxx_xxyy_0, ta1_z_xxxxx_xxyy_1, ta1_z_xxxxx_xxyz_0, ta1_z_xxxxx_xxyz_1, ta1_z_xxxxx_xxz_0, ta1_z_xxxxx_xxz_1, ta1_z_xxxxx_xxzz_0, ta1_z_xxxxx_xxzz_1, ta1_z_xxxxx_xyy_0, ta1_z_xxxxx_xyy_1, ta1_z_xxxxx_xyyy_0, ta1_z_xxxxx_xyyy_1, ta1_z_xxxxx_xyyz_0, ta1_z_xxxxx_xyyz_1, ta1_z_xxxxx_xyz_0, ta1_z_xxxxx_xyz_1, ta1_z_xxxxx_xyzz_0, ta1_z_xxxxx_xyzz_1, ta1_z_xxxxx_xzz_0, ta1_z_xxxxx_xzz_1, ta1_z_xxxxx_xzzz_0, ta1_z_xxxxx_xzzz_1, ta1_z_xxxxx_yyyy_0, ta1_z_xxxxx_yyyy_1, ta1_z_xxxxxz_xxxx_0, ta1_z_xxxxxz_xxxy_0, ta1_z_xxxxxz_xxxz_0, ta1_z_xxxxxz_xxyy_0, ta1_z_xxxxxz_xxyz_0, ta1_z_xxxxxz_xxzz_0, ta1_z_xxxxxz_xyyy_0, ta1_z_xxxxxz_xyyz_0, ta1_z_xxxxxz_xyzz_0, ta1_z_xxxxxz_xzzz_0, ta1_z_xxxxxz_yyyy_0, ta1_z_xxxxxz_yyyz_0, ta1_z_xxxxxz_yyzz_0, ta1_z_xxxxxz_yzzz_0, ta1_z_xxxxxz_zzzz_0, ta1_z_xxxxz_yyyz_0, ta1_z_xxxxz_yyyz_1, ta1_z_xxxxz_yyzz_0, ta1_z_xxxxz_yyzz_1, ta1_z_xxxxz_yzzz_0, ta1_z_xxxxz_yzzz_1, ta1_z_xxxxz_zzzz_0, ta1_z_xxxxz_zzzz_1, ta1_z_xxxz_yyyz_0, ta1_z_xxxz_yyyz_1, ta1_z_xxxz_yyzz_0, ta1_z_xxxz_yyzz_1, ta1_z_xxxz_yzzz_0, ta1_z_xxxz_yzzz_1, ta1_z_xxxz_zzzz_0, ta1_z_xxxz_zzzz_1, ta_xxxxx_xxxx_1, ta_xxxxx_xxxy_1, ta_xxxxx_xxxz_1, ta_xxxxx_xxyy_1, ta_xxxxx_xxyz_1, ta_xxxxx_xxzz_1, ta_xxxxx_xyyy_1, ta_xxxxx_xyyz_1, ta_xxxxx_xyzz_1, ta_xxxxx_xzzz_1, ta_xxxxx_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxz_xxxx_0[i] = ta_xxxxx_xxxx_1[i] + ta1_z_xxxxx_xxxx_0[i] * pa_z[i] - ta1_z_xxxxx_xxxx_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxxy_0[i] = ta_xxxxx_xxxy_1[i] + ta1_z_xxxxx_xxxy_0[i] * pa_z[i] - ta1_z_xxxxx_xxxy_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxxz_0[i] = ta1_z_xxxxx_xxx_0[i] * fe_0 - ta1_z_xxxxx_xxx_1[i] * fe_0 + ta_xxxxx_xxxz_1[i] + ta1_z_xxxxx_xxxz_0[i] * pa_z[i] - ta1_z_xxxxx_xxxz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxyy_0[i] = ta_xxxxx_xxyy_1[i] + ta1_z_xxxxx_xxyy_0[i] * pa_z[i] - ta1_z_xxxxx_xxyy_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxyz_0[i] = ta1_z_xxxxx_xxy_0[i] * fe_0 - ta1_z_xxxxx_xxy_1[i] * fe_0 + ta_xxxxx_xxyz_1[i] + ta1_z_xxxxx_xxyz_0[i] * pa_z[i] - ta1_z_xxxxx_xxyz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxzz_0[i] = 2.0 * ta1_z_xxxxx_xxz_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xxz_1[i] * fe_0 + ta_xxxxx_xxzz_1[i] + ta1_z_xxxxx_xxzz_0[i] * pa_z[i] - ta1_z_xxxxx_xxzz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xyyy_0[i] = ta_xxxxx_xyyy_1[i] + ta1_z_xxxxx_xyyy_0[i] * pa_z[i] - ta1_z_xxxxx_xyyy_1[i] * pc_z[i];

        ta1_z_xxxxxz_xyyz_0[i] = ta1_z_xxxxx_xyy_0[i] * fe_0 - ta1_z_xxxxx_xyy_1[i] * fe_0 + ta_xxxxx_xyyz_1[i] + ta1_z_xxxxx_xyyz_0[i] * pa_z[i] - ta1_z_xxxxx_xyyz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xyzz_0[i] = 2.0 * ta1_z_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xyz_1[i] * fe_0 + ta_xxxxx_xyzz_1[i] + ta1_z_xxxxx_xyzz_0[i] * pa_z[i] - ta1_z_xxxxx_xyzz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xzzz_0[i] = 3.0 * ta1_z_xxxxx_xzz_0[i] * fe_0 - 3.0 * ta1_z_xxxxx_xzz_1[i] * fe_0 + ta_xxxxx_xzzz_1[i] + ta1_z_xxxxx_xzzz_0[i] * pa_z[i] - ta1_z_xxxxx_xzzz_1[i] * pc_z[i];

        ta1_z_xxxxxz_yyyy_0[i] = ta_xxxxx_yyyy_1[i] + ta1_z_xxxxx_yyyy_0[i] * pa_z[i] - ta1_z_xxxxx_yyyy_1[i] * pc_z[i];

        ta1_z_xxxxxz_yyyz_0[i] = 4.0 * ta1_z_xxxz_yyyz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_yyyz_1[i] * fe_0 + ta1_z_xxxxz_yyyz_0[i] * pa_x[i] - ta1_z_xxxxz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxxz_yyzz_0[i] = 4.0 * ta1_z_xxxz_yyzz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_yyzz_1[i] * fe_0 + ta1_z_xxxxz_yyzz_0[i] * pa_x[i] - ta1_z_xxxxz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxxz_yzzz_0[i] = 4.0 * ta1_z_xxxz_yzzz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_yzzz_1[i] * fe_0 + ta1_z_xxxxz_yzzz_0[i] * pa_x[i] - ta1_z_xxxxz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxxz_zzzz_0[i] = 4.0 * ta1_z_xxxz_zzzz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_zzzz_1[i] * fe_0 + ta1_z_xxxxz_zzzz_0[i] * pa_x[i] - ta1_z_xxxxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 885-900 components of targeted buffer : IG

    auto ta1_z_xxxxyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 885);

    auto ta1_z_xxxxyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 886);

    auto ta1_z_xxxxyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 887);

    auto ta1_z_xxxxyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 888);

    auto ta1_z_xxxxyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 889);

    auto ta1_z_xxxxyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 890);

    auto ta1_z_xxxxyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 891);

    auto ta1_z_xxxxyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 892);

    auto ta1_z_xxxxyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 893);

    auto ta1_z_xxxxyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 894);

    auto ta1_z_xxxxyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 895);

    auto ta1_z_xxxxyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 896);

    auto ta1_z_xxxxyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 897);

    auto ta1_z_xxxxyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 898);

    auto ta1_z_xxxxyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 899);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxx_xxxx_0, ta1_z_xxxx_xxxx_1, ta1_z_xxxx_xxxz_0, ta1_z_xxxx_xxxz_1, ta1_z_xxxx_xxzz_0, ta1_z_xxxx_xxzz_1, ta1_z_xxxx_xzzz_0, ta1_z_xxxx_xzzz_1, ta1_z_xxxxy_xxxx_0, ta1_z_xxxxy_xxxx_1, ta1_z_xxxxy_xxxz_0, ta1_z_xxxxy_xxxz_1, ta1_z_xxxxy_xxzz_0, ta1_z_xxxxy_xxzz_1, ta1_z_xxxxy_xzzz_0, ta1_z_xxxxy_xzzz_1, ta1_z_xxxxyy_xxxx_0, ta1_z_xxxxyy_xxxy_0, ta1_z_xxxxyy_xxxz_0, ta1_z_xxxxyy_xxyy_0, ta1_z_xxxxyy_xxyz_0, ta1_z_xxxxyy_xxzz_0, ta1_z_xxxxyy_xyyy_0, ta1_z_xxxxyy_xyyz_0, ta1_z_xxxxyy_xyzz_0, ta1_z_xxxxyy_xzzz_0, ta1_z_xxxxyy_yyyy_0, ta1_z_xxxxyy_yyyz_0, ta1_z_xxxxyy_yyzz_0, ta1_z_xxxxyy_yzzz_0, ta1_z_xxxxyy_zzzz_0, ta1_z_xxxyy_xxxy_0, ta1_z_xxxyy_xxxy_1, ta1_z_xxxyy_xxy_0, ta1_z_xxxyy_xxy_1, ta1_z_xxxyy_xxyy_0, ta1_z_xxxyy_xxyy_1, ta1_z_xxxyy_xxyz_0, ta1_z_xxxyy_xxyz_1, ta1_z_xxxyy_xyy_0, ta1_z_xxxyy_xyy_1, ta1_z_xxxyy_xyyy_0, ta1_z_xxxyy_xyyy_1, ta1_z_xxxyy_xyyz_0, ta1_z_xxxyy_xyyz_1, ta1_z_xxxyy_xyz_0, ta1_z_xxxyy_xyz_1, ta1_z_xxxyy_xyzz_0, ta1_z_xxxyy_xyzz_1, ta1_z_xxxyy_yyy_0, ta1_z_xxxyy_yyy_1, ta1_z_xxxyy_yyyy_0, ta1_z_xxxyy_yyyy_1, ta1_z_xxxyy_yyyz_0, ta1_z_xxxyy_yyyz_1, ta1_z_xxxyy_yyz_0, ta1_z_xxxyy_yyz_1, ta1_z_xxxyy_yyzz_0, ta1_z_xxxyy_yyzz_1, ta1_z_xxxyy_yzz_0, ta1_z_xxxyy_yzz_1, ta1_z_xxxyy_yzzz_0, ta1_z_xxxyy_yzzz_1, ta1_z_xxxyy_zzzz_0, ta1_z_xxxyy_zzzz_1, ta1_z_xxyy_xxxy_0, ta1_z_xxyy_xxxy_1, ta1_z_xxyy_xxyy_0, ta1_z_xxyy_xxyy_1, ta1_z_xxyy_xxyz_0, ta1_z_xxyy_xxyz_1, ta1_z_xxyy_xyyy_0, ta1_z_xxyy_xyyy_1, ta1_z_xxyy_xyyz_0, ta1_z_xxyy_xyyz_1, ta1_z_xxyy_xyzz_0, ta1_z_xxyy_xyzz_1, ta1_z_xxyy_yyyy_0, ta1_z_xxyy_yyyy_1, ta1_z_xxyy_yyyz_0, ta1_z_xxyy_yyyz_1, ta1_z_xxyy_yyzz_0, ta1_z_xxyy_yyzz_1, ta1_z_xxyy_yzzz_0, ta1_z_xxyy_yzzz_1, ta1_z_xxyy_zzzz_0, ta1_z_xxyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyy_xxxx_0[i] = ta1_z_xxxx_xxxx_0[i] * fe_0 - ta1_z_xxxx_xxxx_1[i] * fe_0 + ta1_z_xxxxy_xxxx_0[i] * pa_y[i] - ta1_z_xxxxy_xxxx_1[i] * pc_y[i];

        ta1_z_xxxxyy_xxxy_0[i] = 3.0 * ta1_z_xxyy_xxxy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xxxyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_xxxyy_xxy_1[i] * fe_0 + ta1_z_xxxyy_xxxy_0[i] * pa_x[i] - ta1_z_xxxyy_xxxy_1[i] * pc_x[i];

        ta1_z_xxxxyy_xxxz_0[i] = ta1_z_xxxx_xxxz_0[i] * fe_0 - ta1_z_xxxx_xxxz_1[i] * fe_0 + ta1_z_xxxxy_xxxz_0[i] * pa_y[i] - ta1_z_xxxxy_xxxz_1[i] * pc_y[i];

        ta1_z_xxxxyy_xxyy_0[i] = 3.0 * ta1_z_xxyy_xxyy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxxyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_xxxyy_xyy_1[i] * fe_0 + ta1_z_xxxyy_xxyy_0[i] * pa_x[i] - ta1_z_xxxyy_xxyy_1[i] * pc_x[i];

        ta1_z_xxxxyy_xxyz_0[i] = 3.0 * ta1_z_xxyy_xxyz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxyy_xyz_1[i] * fe_0 + ta1_z_xxxyy_xxyz_0[i] * pa_x[i] - ta1_z_xxxyy_xxyz_1[i] * pc_x[i];

        ta1_z_xxxxyy_xxzz_0[i] = ta1_z_xxxx_xxzz_0[i] * fe_0 - ta1_z_xxxx_xxzz_1[i] * fe_0 + ta1_z_xxxxy_xxzz_0[i] * pa_y[i] - ta1_z_xxxxy_xxzz_1[i] * pc_y[i];

        ta1_z_xxxxyy_xyyy_0[i] = 3.0 * ta1_z_xxyy_xyyy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xyyy_1[i] * fe_0 + ta1_z_xxxyy_yyy_0[i] * fe_0 - ta1_z_xxxyy_yyy_1[i] * fe_0 + ta1_z_xxxyy_xyyy_0[i] * pa_x[i] - ta1_z_xxxyy_xyyy_1[i] * pc_x[i];

        ta1_z_xxxxyy_xyyz_0[i] = 3.0 * ta1_z_xxyy_xyyz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xyyz_1[i] * fe_0 + ta1_z_xxxyy_yyz_0[i] * fe_0 - ta1_z_xxxyy_yyz_1[i] * fe_0 + ta1_z_xxxyy_xyyz_0[i] * pa_x[i] - ta1_z_xxxyy_xyyz_1[i] * pc_x[i];

        ta1_z_xxxxyy_xyzz_0[i] = 3.0 * ta1_z_xxyy_xyzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xyzz_1[i] * fe_0 + ta1_z_xxxyy_yzz_0[i] * fe_0 - ta1_z_xxxyy_yzz_1[i] * fe_0 + ta1_z_xxxyy_xyzz_0[i] * pa_x[i] - ta1_z_xxxyy_xyzz_1[i] * pc_x[i];

        ta1_z_xxxxyy_xzzz_0[i] = ta1_z_xxxx_xzzz_0[i] * fe_0 - ta1_z_xxxx_xzzz_1[i] * fe_0 + ta1_z_xxxxy_xzzz_0[i] * pa_y[i] - ta1_z_xxxxy_xzzz_1[i] * pc_y[i];

        ta1_z_xxxxyy_yyyy_0[i] = 3.0 * ta1_z_xxyy_yyyy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yyyy_1[i] * fe_0 + ta1_z_xxxyy_yyyy_0[i] * pa_x[i] - ta1_z_xxxyy_yyyy_1[i] * pc_x[i];

        ta1_z_xxxxyy_yyyz_0[i] = 3.0 * ta1_z_xxyy_yyyz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yyyz_1[i] * fe_0 + ta1_z_xxxyy_yyyz_0[i] * pa_x[i] - ta1_z_xxxyy_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxyy_yyzz_0[i] = 3.0 * ta1_z_xxyy_yyzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yyzz_1[i] * fe_0 + ta1_z_xxxyy_yyzz_0[i] * pa_x[i] - ta1_z_xxxyy_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxyy_yzzz_0[i] = 3.0 * ta1_z_xxyy_yzzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yzzz_1[i] * fe_0 + ta1_z_xxxyy_yzzz_0[i] * pa_x[i] - ta1_z_xxxyy_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxyy_zzzz_0[i] = 3.0 * ta1_z_xxyy_zzzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_zzzz_1[i] * fe_0 + ta1_z_xxxyy_zzzz_0[i] * pa_x[i] - ta1_z_xxxyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 900-915 components of targeted buffer : IG

    auto ta1_z_xxxxyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 900);

    auto ta1_z_xxxxyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 901);

    auto ta1_z_xxxxyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 902);

    auto ta1_z_xxxxyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 903);

    auto ta1_z_xxxxyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 904);

    auto ta1_z_xxxxyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 905);

    auto ta1_z_xxxxyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 906);

    auto ta1_z_xxxxyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 907);

    auto ta1_z_xxxxyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 908);

    auto ta1_z_xxxxyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 909);

    auto ta1_z_xxxxyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 910);

    auto ta1_z_xxxxyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 911);

    auto ta1_z_xxxxyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 912);

    auto ta1_z_xxxxyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 913);

    auto ta1_z_xxxxyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 914);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxxxy_xxxy_0, ta1_z_xxxxy_xxxy_1, ta1_z_xxxxy_xxyy_0, ta1_z_xxxxy_xxyy_1, ta1_z_xxxxy_xyyy_0, ta1_z_xxxxy_xyyy_1, ta1_z_xxxxy_yyyy_0, ta1_z_xxxxy_yyyy_1, ta1_z_xxxxyz_xxxx_0, ta1_z_xxxxyz_xxxy_0, ta1_z_xxxxyz_xxxz_0, ta1_z_xxxxyz_xxyy_0, ta1_z_xxxxyz_xxyz_0, ta1_z_xxxxyz_xxzz_0, ta1_z_xxxxyz_xyyy_0, ta1_z_xxxxyz_xyyz_0, ta1_z_xxxxyz_xyzz_0, ta1_z_xxxxyz_xzzz_0, ta1_z_xxxxyz_yyyy_0, ta1_z_xxxxyz_yyyz_0, ta1_z_xxxxyz_yyzz_0, ta1_z_xxxxyz_yzzz_0, ta1_z_xxxxyz_zzzz_0, ta1_z_xxxxz_xxxx_0, ta1_z_xxxxz_xxxx_1, ta1_z_xxxxz_xxxz_0, ta1_z_xxxxz_xxxz_1, ta1_z_xxxxz_xxyz_0, ta1_z_xxxxz_xxyz_1, ta1_z_xxxxz_xxz_0, ta1_z_xxxxz_xxz_1, ta1_z_xxxxz_xxzz_0, ta1_z_xxxxz_xxzz_1, ta1_z_xxxxz_xyyz_0, ta1_z_xxxxz_xyyz_1, ta1_z_xxxxz_xyz_0, ta1_z_xxxxz_xyz_1, ta1_z_xxxxz_xyzz_0, ta1_z_xxxxz_xyzz_1, ta1_z_xxxxz_xzz_0, ta1_z_xxxxz_xzz_1, ta1_z_xxxxz_xzzz_0, ta1_z_xxxxz_xzzz_1, ta1_z_xxxxz_zzzz_0, ta1_z_xxxxz_zzzz_1, ta1_z_xxxyz_yyyz_0, ta1_z_xxxyz_yyyz_1, ta1_z_xxxyz_yyzz_0, ta1_z_xxxyz_yyzz_1, ta1_z_xxxyz_yzzz_0, ta1_z_xxxyz_yzzz_1, ta1_z_xxyz_yyyz_0, ta1_z_xxyz_yyyz_1, ta1_z_xxyz_yyzz_0, ta1_z_xxyz_yyzz_1, ta1_z_xxyz_yzzz_0, ta1_z_xxyz_yzzz_1, ta_xxxxy_xxxy_1, ta_xxxxy_xxyy_1, ta_xxxxy_xyyy_1, ta_xxxxy_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyz_xxxx_0[i] = ta1_z_xxxxz_xxxx_0[i] * pa_y[i] - ta1_z_xxxxz_xxxx_1[i] * pc_y[i];

        ta1_z_xxxxyz_xxxy_0[i] = ta_xxxxy_xxxy_1[i] + ta1_z_xxxxy_xxxy_0[i] * pa_z[i] - ta1_z_xxxxy_xxxy_1[i] * pc_z[i];

        ta1_z_xxxxyz_xxxz_0[i] = ta1_z_xxxxz_xxxz_0[i] * pa_y[i] - ta1_z_xxxxz_xxxz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xxyy_0[i] = ta_xxxxy_xxyy_1[i] + ta1_z_xxxxy_xxyy_0[i] * pa_z[i] - ta1_z_xxxxy_xxyy_1[i] * pc_z[i];

        ta1_z_xxxxyz_xxyz_0[i] = ta1_z_xxxxz_xxz_0[i] * fe_0 - ta1_z_xxxxz_xxz_1[i] * fe_0 + ta1_z_xxxxz_xxyz_0[i] * pa_y[i] - ta1_z_xxxxz_xxyz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xxzz_0[i] = ta1_z_xxxxz_xxzz_0[i] * pa_y[i] - ta1_z_xxxxz_xxzz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xyyy_0[i] = ta_xxxxy_xyyy_1[i] + ta1_z_xxxxy_xyyy_0[i] * pa_z[i] - ta1_z_xxxxy_xyyy_1[i] * pc_z[i];

        ta1_z_xxxxyz_xyyz_0[i] = 2.0 * ta1_z_xxxxz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxxz_xyz_1[i] * fe_0 + ta1_z_xxxxz_xyyz_0[i] * pa_y[i] - ta1_z_xxxxz_xyyz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xyzz_0[i] = ta1_z_xxxxz_xzz_0[i] * fe_0 - ta1_z_xxxxz_xzz_1[i] * fe_0 + ta1_z_xxxxz_xyzz_0[i] * pa_y[i] - ta1_z_xxxxz_xyzz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xzzz_0[i] = ta1_z_xxxxz_xzzz_0[i] * pa_y[i] - ta1_z_xxxxz_xzzz_1[i] * pc_y[i];

        ta1_z_xxxxyz_yyyy_0[i] = ta_xxxxy_yyyy_1[i] + ta1_z_xxxxy_yyyy_0[i] * pa_z[i] - ta1_z_xxxxy_yyyy_1[i] * pc_z[i];

        ta1_z_xxxxyz_yyyz_0[i] = 3.0 * ta1_z_xxyz_yyyz_0[i] * fe_0 - 3.0 * ta1_z_xxyz_yyyz_1[i] * fe_0 + ta1_z_xxxyz_yyyz_0[i] * pa_x[i] - ta1_z_xxxyz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxyz_yyzz_0[i] = 3.0 * ta1_z_xxyz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_xxyz_yyzz_1[i] * fe_0 + ta1_z_xxxyz_yyzz_0[i] * pa_x[i] - ta1_z_xxxyz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxyz_yzzz_0[i] = 3.0 * ta1_z_xxyz_yzzz_0[i] * fe_0 - 3.0 * ta1_z_xxyz_yzzz_1[i] * fe_0 + ta1_z_xxxyz_yzzz_0[i] * pa_x[i] - ta1_z_xxxyz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxyz_zzzz_0[i] = ta1_z_xxxxz_zzzz_0[i] * pa_y[i] - ta1_z_xxxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 915-930 components of targeted buffer : IG

    auto ta1_z_xxxxzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 915);

    auto ta1_z_xxxxzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 916);

    auto ta1_z_xxxxzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 917);

    auto ta1_z_xxxxzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 918);

    auto ta1_z_xxxxzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 919);

    auto ta1_z_xxxxzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 920);

    auto ta1_z_xxxxzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 921);

    auto ta1_z_xxxxzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 922);

    auto ta1_z_xxxxzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 923);

    auto ta1_z_xxxxzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 924);

    auto ta1_z_xxxxzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 925);

    auto ta1_z_xxxxzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 926);

    auto ta1_z_xxxxzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 927);

    auto ta1_z_xxxxzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 928);

    auto ta1_z_xxxxzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 929);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxx_xxxx_0, ta1_z_xxxx_xxxx_1, ta1_z_xxxx_xxxy_0, ta1_z_xxxx_xxxy_1, ta1_z_xxxx_xxyy_0, ta1_z_xxxx_xxyy_1, ta1_z_xxxx_xyyy_0, ta1_z_xxxx_xyyy_1, ta1_z_xxxxz_xxxx_0, ta1_z_xxxxz_xxxx_1, ta1_z_xxxxz_xxxy_0, ta1_z_xxxxz_xxxy_1, ta1_z_xxxxz_xxyy_0, ta1_z_xxxxz_xxyy_1, ta1_z_xxxxz_xyyy_0, ta1_z_xxxxz_xyyy_1, ta1_z_xxxxzz_xxxx_0, ta1_z_xxxxzz_xxxy_0, ta1_z_xxxxzz_xxxz_0, ta1_z_xxxxzz_xxyy_0, ta1_z_xxxxzz_xxyz_0, ta1_z_xxxxzz_xxzz_0, ta1_z_xxxxzz_xyyy_0, ta1_z_xxxxzz_xyyz_0, ta1_z_xxxxzz_xyzz_0, ta1_z_xxxxzz_xzzz_0, ta1_z_xxxxzz_yyyy_0, ta1_z_xxxxzz_yyyz_0, ta1_z_xxxxzz_yyzz_0, ta1_z_xxxxzz_yzzz_0, ta1_z_xxxxzz_zzzz_0, ta1_z_xxxzz_xxxz_0, ta1_z_xxxzz_xxxz_1, ta1_z_xxxzz_xxyz_0, ta1_z_xxxzz_xxyz_1, ta1_z_xxxzz_xxz_0, ta1_z_xxxzz_xxz_1, ta1_z_xxxzz_xxzz_0, ta1_z_xxxzz_xxzz_1, ta1_z_xxxzz_xyyz_0, ta1_z_xxxzz_xyyz_1, ta1_z_xxxzz_xyz_0, ta1_z_xxxzz_xyz_1, ta1_z_xxxzz_xyzz_0, ta1_z_xxxzz_xyzz_1, ta1_z_xxxzz_xzz_0, ta1_z_xxxzz_xzz_1, ta1_z_xxxzz_xzzz_0, ta1_z_xxxzz_xzzz_1, ta1_z_xxxzz_yyyy_0, ta1_z_xxxzz_yyyy_1, ta1_z_xxxzz_yyyz_0, ta1_z_xxxzz_yyyz_1, ta1_z_xxxzz_yyz_0, ta1_z_xxxzz_yyz_1, ta1_z_xxxzz_yyzz_0, ta1_z_xxxzz_yyzz_1, ta1_z_xxxzz_yzz_0, ta1_z_xxxzz_yzz_1, ta1_z_xxxzz_yzzz_0, ta1_z_xxxzz_yzzz_1, ta1_z_xxxzz_zzz_0, ta1_z_xxxzz_zzz_1, ta1_z_xxxzz_zzzz_0, ta1_z_xxxzz_zzzz_1, ta1_z_xxzz_xxxz_0, ta1_z_xxzz_xxxz_1, ta1_z_xxzz_xxyz_0, ta1_z_xxzz_xxyz_1, ta1_z_xxzz_xxzz_0, ta1_z_xxzz_xxzz_1, ta1_z_xxzz_xyyz_0, ta1_z_xxzz_xyyz_1, ta1_z_xxzz_xyzz_0, ta1_z_xxzz_xyzz_1, ta1_z_xxzz_xzzz_0, ta1_z_xxzz_xzzz_1, ta1_z_xxzz_yyyy_0, ta1_z_xxzz_yyyy_1, ta1_z_xxzz_yyyz_0, ta1_z_xxzz_yyyz_1, ta1_z_xxzz_yyzz_0, ta1_z_xxzz_yyzz_1, ta1_z_xxzz_yzzz_0, ta1_z_xxzz_yzzz_1, ta1_z_xxzz_zzzz_0, ta1_z_xxzz_zzzz_1, ta_xxxxz_xxxx_1, ta_xxxxz_xxxy_1, ta_xxxxz_xxyy_1, ta_xxxxz_xyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxzz_xxxx_0[i] = ta1_z_xxxx_xxxx_0[i] * fe_0 - ta1_z_xxxx_xxxx_1[i] * fe_0 + ta_xxxxz_xxxx_1[i] + ta1_z_xxxxz_xxxx_0[i] * pa_z[i] - ta1_z_xxxxz_xxxx_1[i] * pc_z[i];

        ta1_z_xxxxzz_xxxy_0[i] = ta1_z_xxxx_xxxy_0[i] * fe_0 - ta1_z_xxxx_xxxy_1[i] * fe_0 + ta_xxxxz_xxxy_1[i] + ta1_z_xxxxz_xxxy_0[i] * pa_z[i] - ta1_z_xxxxz_xxxy_1[i] * pc_z[i];

        ta1_z_xxxxzz_xxxz_0[i] = 3.0 * ta1_z_xxzz_xxxz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xxxzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_xxxzz_xxz_1[i] * fe_0 + ta1_z_xxxzz_xxxz_0[i] * pa_x[i] - ta1_z_xxxzz_xxxz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xxyy_0[i] = ta1_z_xxxx_xxyy_0[i] * fe_0 - ta1_z_xxxx_xxyy_1[i] * fe_0 + ta_xxxxz_xxyy_1[i] + ta1_z_xxxxz_xxyy_0[i] * pa_z[i] - ta1_z_xxxxz_xxyy_1[i] * pc_z[i];

        ta1_z_xxxxzz_xxyz_0[i] = 3.0 * ta1_z_xxzz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxzz_xyz_1[i] * fe_0 + ta1_z_xxxzz_xxyz_0[i] * pa_x[i] - ta1_z_xxxzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xxzz_0[i] = 3.0 * ta1_z_xxzz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xxxzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_xxxzz_xzz_1[i] * fe_0 + ta1_z_xxxzz_xxzz_0[i] * pa_x[i] - ta1_z_xxxzz_xxzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xyyy_0[i] = ta1_z_xxxx_xyyy_0[i] * fe_0 - ta1_z_xxxx_xyyy_1[i] * fe_0 + ta_xxxxz_xyyy_1[i] + ta1_z_xxxxz_xyyy_0[i] * pa_z[i] - ta1_z_xxxxz_xyyy_1[i] * pc_z[i];

        ta1_z_xxxxzz_xyyz_0[i] = 3.0 * ta1_z_xxzz_xyyz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xyyz_1[i] * fe_0 + ta1_z_xxxzz_yyz_0[i] * fe_0 - ta1_z_xxxzz_yyz_1[i] * fe_0 + ta1_z_xxxzz_xyyz_0[i] * pa_x[i] - ta1_z_xxxzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xyzz_0[i] = 3.0 * ta1_z_xxzz_xyzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xyzz_1[i] * fe_0 + ta1_z_xxxzz_yzz_0[i] * fe_0 - ta1_z_xxxzz_yzz_1[i] * fe_0 + ta1_z_xxxzz_xyzz_0[i] * pa_x[i] - ta1_z_xxxzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xzzz_0[i] = 3.0 * ta1_z_xxzz_xzzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xzzz_1[i] * fe_0 + ta1_z_xxxzz_zzz_0[i] * fe_0 - ta1_z_xxxzz_zzz_1[i] * fe_0 + ta1_z_xxxzz_xzzz_0[i] * pa_x[i] - ta1_z_xxxzz_xzzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_yyyy_0[i] = 3.0 * ta1_z_xxzz_yyyy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yyyy_1[i] * fe_0 + ta1_z_xxxzz_yyyy_0[i] * pa_x[i] - ta1_z_xxxzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxxxzz_yyyz_0[i] = 3.0 * ta1_z_xxzz_yyyz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yyyz_1[i] * fe_0 + ta1_z_xxxzz_yyyz_0[i] * pa_x[i] - ta1_z_xxxzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxxzz_yyzz_0[i] = 3.0 * ta1_z_xxzz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yyzz_1[i] * fe_0 + ta1_z_xxxzz_yyzz_0[i] * pa_x[i] - ta1_z_xxxzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_yzzz_0[i] = 3.0 * ta1_z_xxzz_yzzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yzzz_1[i] * fe_0 + ta1_z_xxxzz_yzzz_0[i] * pa_x[i] - ta1_z_xxxzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_zzzz_0[i] = 3.0 * ta1_z_xxzz_zzzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_zzzz_1[i] * fe_0 + ta1_z_xxxzz_zzzz_0[i] * pa_x[i] - ta1_z_xxxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 930-945 components of targeted buffer : IG

    auto ta1_z_xxxyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 930);

    auto ta1_z_xxxyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 931);

    auto ta1_z_xxxyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 932);

    auto ta1_z_xxxyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 933);

    auto ta1_z_xxxyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 934);

    auto ta1_z_xxxyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 935);

    auto ta1_z_xxxyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 936);

    auto ta1_z_xxxyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 937);

    auto ta1_z_xxxyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 938);

    auto ta1_z_xxxyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 939);

    auto ta1_z_xxxyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 940);

    auto ta1_z_xxxyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 941);

    auto ta1_z_xxxyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 942);

    auto ta1_z_xxxyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 943);

    auto ta1_z_xxxyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 944);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxy_xxxx_0, ta1_z_xxxy_xxxx_1, ta1_z_xxxy_xxxz_0, ta1_z_xxxy_xxxz_1, ta1_z_xxxy_xxzz_0, ta1_z_xxxy_xxzz_1, ta1_z_xxxy_xzzz_0, ta1_z_xxxy_xzzz_1, ta1_z_xxxyy_xxxx_0, ta1_z_xxxyy_xxxx_1, ta1_z_xxxyy_xxxz_0, ta1_z_xxxyy_xxxz_1, ta1_z_xxxyy_xxzz_0, ta1_z_xxxyy_xxzz_1, ta1_z_xxxyy_xzzz_0, ta1_z_xxxyy_xzzz_1, ta1_z_xxxyyy_xxxx_0, ta1_z_xxxyyy_xxxy_0, ta1_z_xxxyyy_xxxz_0, ta1_z_xxxyyy_xxyy_0, ta1_z_xxxyyy_xxyz_0, ta1_z_xxxyyy_xxzz_0, ta1_z_xxxyyy_xyyy_0, ta1_z_xxxyyy_xyyz_0, ta1_z_xxxyyy_xyzz_0, ta1_z_xxxyyy_xzzz_0, ta1_z_xxxyyy_yyyy_0, ta1_z_xxxyyy_yyyz_0, ta1_z_xxxyyy_yyzz_0, ta1_z_xxxyyy_yzzz_0, ta1_z_xxxyyy_zzzz_0, ta1_z_xxyyy_xxxy_0, ta1_z_xxyyy_xxxy_1, ta1_z_xxyyy_xxy_0, ta1_z_xxyyy_xxy_1, ta1_z_xxyyy_xxyy_0, ta1_z_xxyyy_xxyy_1, ta1_z_xxyyy_xxyz_0, ta1_z_xxyyy_xxyz_1, ta1_z_xxyyy_xyy_0, ta1_z_xxyyy_xyy_1, ta1_z_xxyyy_xyyy_0, ta1_z_xxyyy_xyyy_1, ta1_z_xxyyy_xyyz_0, ta1_z_xxyyy_xyyz_1, ta1_z_xxyyy_xyz_0, ta1_z_xxyyy_xyz_1, ta1_z_xxyyy_xyzz_0, ta1_z_xxyyy_xyzz_1, ta1_z_xxyyy_yyy_0, ta1_z_xxyyy_yyy_1, ta1_z_xxyyy_yyyy_0, ta1_z_xxyyy_yyyy_1, ta1_z_xxyyy_yyyz_0, ta1_z_xxyyy_yyyz_1, ta1_z_xxyyy_yyz_0, ta1_z_xxyyy_yyz_1, ta1_z_xxyyy_yyzz_0, ta1_z_xxyyy_yyzz_1, ta1_z_xxyyy_yzz_0, ta1_z_xxyyy_yzz_1, ta1_z_xxyyy_yzzz_0, ta1_z_xxyyy_yzzz_1, ta1_z_xxyyy_zzzz_0, ta1_z_xxyyy_zzzz_1, ta1_z_xyyy_xxxy_0, ta1_z_xyyy_xxxy_1, ta1_z_xyyy_xxyy_0, ta1_z_xyyy_xxyy_1, ta1_z_xyyy_xxyz_0, ta1_z_xyyy_xxyz_1, ta1_z_xyyy_xyyy_0, ta1_z_xyyy_xyyy_1, ta1_z_xyyy_xyyz_0, ta1_z_xyyy_xyyz_1, ta1_z_xyyy_xyzz_0, ta1_z_xyyy_xyzz_1, ta1_z_xyyy_yyyy_0, ta1_z_xyyy_yyyy_1, ta1_z_xyyy_yyyz_0, ta1_z_xyyy_yyyz_1, ta1_z_xyyy_yyzz_0, ta1_z_xyyy_yyzz_1, ta1_z_xyyy_yzzz_0, ta1_z_xyyy_yzzz_1, ta1_z_xyyy_zzzz_0, ta1_z_xyyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyy_xxxx_0[i] = 2.0 * ta1_z_xxxy_xxxx_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xxxx_1[i] * fe_0 + ta1_z_xxxyy_xxxx_0[i] * pa_y[i] - ta1_z_xxxyy_xxxx_1[i] * pc_y[i];

        ta1_z_xxxyyy_xxxy_0[i] = 2.0 * ta1_z_xyyy_xxxy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xxyyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_xxyyy_xxy_1[i] * fe_0 + ta1_z_xxyyy_xxxy_0[i] * pa_x[i] - ta1_z_xxyyy_xxxy_1[i] * pc_x[i];

        ta1_z_xxxyyy_xxxz_0[i] = 2.0 * ta1_z_xxxy_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xxxz_1[i] * fe_0 + ta1_z_xxxyy_xxxz_0[i] * pa_y[i] - ta1_z_xxxyy_xxxz_1[i] * pc_y[i];

        ta1_z_xxxyyy_xxyy_0[i] = 2.0 * ta1_z_xyyy_xxyy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxyyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_xxyyy_xyy_1[i] * fe_0 + ta1_z_xxyyy_xxyy_0[i] * pa_x[i] - ta1_z_xxyyy_xxyy_1[i] * pc_x[i];

        ta1_z_xxxyyy_xxyz_0[i] = 2.0 * ta1_z_xyyy_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxyyy_xyz_1[i] * fe_0 + ta1_z_xxyyy_xxyz_0[i] * pa_x[i] - ta1_z_xxyyy_xxyz_1[i] * pc_x[i];

        ta1_z_xxxyyy_xxzz_0[i] = 2.0 * ta1_z_xxxy_xxzz_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xxzz_1[i] * fe_0 + ta1_z_xxxyy_xxzz_0[i] * pa_y[i] - ta1_z_xxxyy_xxzz_1[i] * pc_y[i];

        ta1_z_xxxyyy_xyyy_0[i] = 2.0 * ta1_z_xyyy_xyyy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xyyy_1[i] * fe_0 + ta1_z_xxyyy_yyy_0[i] * fe_0 - ta1_z_xxyyy_yyy_1[i] * fe_0 + ta1_z_xxyyy_xyyy_0[i] * pa_x[i] - ta1_z_xxyyy_xyyy_1[i] * pc_x[i];

        ta1_z_xxxyyy_xyyz_0[i] = 2.0 * ta1_z_xyyy_xyyz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xyyz_1[i] * fe_0 + ta1_z_xxyyy_yyz_0[i] * fe_0 - ta1_z_xxyyy_yyz_1[i] * fe_0 + ta1_z_xxyyy_xyyz_0[i] * pa_x[i] - ta1_z_xxyyy_xyyz_1[i] * pc_x[i];

        ta1_z_xxxyyy_xyzz_0[i] = 2.0 * ta1_z_xyyy_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xyzz_1[i] * fe_0 + ta1_z_xxyyy_yzz_0[i] * fe_0 - ta1_z_xxyyy_yzz_1[i] * fe_0 + ta1_z_xxyyy_xyzz_0[i] * pa_x[i] - ta1_z_xxyyy_xyzz_1[i] * pc_x[i];

        ta1_z_xxxyyy_xzzz_0[i] = 2.0 * ta1_z_xxxy_xzzz_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xzzz_1[i] * fe_0 + ta1_z_xxxyy_xzzz_0[i] * pa_y[i] - ta1_z_xxxyy_xzzz_1[i] * pc_y[i];

        ta1_z_xxxyyy_yyyy_0[i] = 2.0 * ta1_z_xyyy_yyyy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yyyy_1[i] * fe_0 + ta1_z_xxyyy_yyyy_0[i] * pa_x[i] - ta1_z_xxyyy_yyyy_1[i] * pc_x[i];

        ta1_z_xxxyyy_yyyz_0[i] = 2.0 * ta1_z_xyyy_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yyyz_1[i] * fe_0 + ta1_z_xxyyy_yyyz_0[i] * pa_x[i] - ta1_z_xxyyy_yyyz_1[i] * pc_x[i];

        ta1_z_xxxyyy_yyzz_0[i] = 2.0 * ta1_z_xyyy_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yyzz_1[i] * fe_0 + ta1_z_xxyyy_yyzz_0[i] * pa_x[i] - ta1_z_xxyyy_yyzz_1[i] * pc_x[i];

        ta1_z_xxxyyy_yzzz_0[i] = 2.0 * ta1_z_xyyy_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yzzz_1[i] * fe_0 + ta1_z_xxyyy_yzzz_0[i] * pa_x[i] - ta1_z_xxyyy_yzzz_1[i] * pc_x[i];

        ta1_z_xxxyyy_zzzz_0[i] = 2.0 * ta1_z_xyyy_zzzz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_zzzz_1[i] * fe_0 + ta1_z_xxyyy_zzzz_0[i] * pa_x[i] - ta1_z_xxyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 945-960 components of targeted buffer : IG

    auto ta1_z_xxxyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 945);

    auto ta1_z_xxxyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 946);

    auto ta1_z_xxxyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 947);

    auto ta1_z_xxxyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 948);

    auto ta1_z_xxxyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 949);

    auto ta1_z_xxxyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 950);

    auto ta1_z_xxxyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 951);

    auto ta1_z_xxxyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 952);

    auto ta1_z_xxxyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 953);

    auto ta1_z_xxxyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 954);

    auto ta1_z_xxxyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 955);

    auto ta1_z_xxxyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 956);

    auto ta1_z_xxxyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 957);

    auto ta1_z_xxxyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 958);

    auto ta1_z_xxxyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 959);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxxyy_xxxx_0, ta1_z_xxxyy_xxxx_1, ta1_z_xxxyy_xxxy_0, ta1_z_xxxyy_xxxy_1, ta1_z_xxxyy_xxy_0, ta1_z_xxxyy_xxy_1, ta1_z_xxxyy_xxyy_0, ta1_z_xxxyy_xxyy_1, ta1_z_xxxyy_xxyz_0, ta1_z_xxxyy_xxyz_1, ta1_z_xxxyy_xyy_0, ta1_z_xxxyy_xyy_1, ta1_z_xxxyy_xyyy_0, ta1_z_xxxyy_xyyy_1, ta1_z_xxxyy_xyyz_0, ta1_z_xxxyy_xyyz_1, ta1_z_xxxyy_xyz_0, ta1_z_xxxyy_xyz_1, ta1_z_xxxyy_xyzz_0, ta1_z_xxxyy_xyzz_1, ta1_z_xxxyy_yyyy_0, ta1_z_xxxyy_yyyy_1, ta1_z_xxxyyz_xxxx_0, ta1_z_xxxyyz_xxxy_0, ta1_z_xxxyyz_xxxz_0, ta1_z_xxxyyz_xxyy_0, ta1_z_xxxyyz_xxyz_0, ta1_z_xxxyyz_xxzz_0, ta1_z_xxxyyz_xyyy_0, ta1_z_xxxyyz_xyyz_0, ta1_z_xxxyyz_xyzz_0, ta1_z_xxxyyz_xzzz_0, ta1_z_xxxyyz_yyyy_0, ta1_z_xxxyyz_yyyz_0, ta1_z_xxxyyz_yyzz_0, ta1_z_xxxyyz_yzzz_0, ta1_z_xxxyyz_zzzz_0, ta1_z_xxxyz_xxxz_0, ta1_z_xxxyz_xxxz_1, ta1_z_xxxyz_xxzz_0, ta1_z_xxxyz_xxzz_1, ta1_z_xxxyz_xzzz_0, ta1_z_xxxyz_xzzz_1, ta1_z_xxxz_xxxz_0, ta1_z_xxxz_xxxz_1, ta1_z_xxxz_xxzz_0, ta1_z_xxxz_xxzz_1, ta1_z_xxxz_xzzz_0, ta1_z_xxxz_xzzz_1, ta1_z_xxyyz_yyyz_0, ta1_z_xxyyz_yyyz_1, ta1_z_xxyyz_yyzz_0, ta1_z_xxyyz_yyzz_1, ta1_z_xxyyz_yzzz_0, ta1_z_xxyyz_yzzz_1, ta1_z_xxyyz_zzzz_0, ta1_z_xxyyz_zzzz_1, ta1_z_xyyz_yyyz_0, ta1_z_xyyz_yyyz_1, ta1_z_xyyz_yyzz_0, ta1_z_xyyz_yyzz_1, ta1_z_xyyz_yzzz_0, ta1_z_xyyz_yzzz_1, ta1_z_xyyz_zzzz_0, ta1_z_xyyz_zzzz_1, ta_xxxyy_xxxx_1, ta_xxxyy_xxxy_1, ta_xxxyy_xxyy_1, ta_xxxyy_xxyz_1, ta_xxxyy_xyyy_1, ta_xxxyy_xyyz_1, ta_xxxyy_xyzz_1, ta_xxxyy_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyz_xxxx_0[i] = ta_xxxyy_xxxx_1[i] + ta1_z_xxxyy_xxxx_0[i] * pa_z[i] - ta1_z_xxxyy_xxxx_1[i] * pc_z[i];

        ta1_z_xxxyyz_xxxy_0[i] = ta_xxxyy_xxxy_1[i] + ta1_z_xxxyy_xxxy_0[i] * pa_z[i] - ta1_z_xxxyy_xxxy_1[i] * pc_z[i];

        ta1_z_xxxyyz_xxxz_0[i] = ta1_z_xxxz_xxxz_0[i] * fe_0 - ta1_z_xxxz_xxxz_1[i] * fe_0 + ta1_z_xxxyz_xxxz_0[i] * pa_y[i] - ta1_z_xxxyz_xxxz_1[i] * pc_y[i];

        ta1_z_xxxyyz_xxyy_0[i] = ta_xxxyy_xxyy_1[i] + ta1_z_xxxyy_xxyy_0[i] * pa_z[i] - ta1_z_xxxyy_xxyy_1[i] * pc_z[i];

        ta1_z_xxxyyz_xxyz_0[i] = ta1_z_xxxyy_xxy_0[i] * fe_0 - ta1_z_xxxyy_xxy_1[i] * fe_0 + ta_xxxyy_xxyz_1[i] + ta1_z_xxxyy_xxyz_0[i] * pa_z[i] - ta1_z_xxxyy_xxyz_1[i] * pc_z[i];

        ta1_z_xxxyyz_xxzz_0[i] = ta1_z_xxxz_xxzz_0[i] * fe_0 - ta1_z_xxxz_xxzz_1[i] * fe_0 + ta1_z_xxxyz_xxzz_0[i] * pa_y[i] - ta1_z_xxxyz_xxzz_1[i] * pc_y[i];

        ta1_z_xxxyyz_xyyy_0[i] = ta_xxxyy_xyyy_1[i] + ta1_z_xxxyy_xyyy_0[i] * pa_z[i] - ta1_z_xxxyy_xyyy_1[i] * pc_z[i];

        ta1_z_xxxyyz_xyyz_0[i] = ta1_z_xxxyy_xyy_0[i] * fe_0 - ta1_z_xxxyy_xyy_1[i] * fe_0 + ta_xxxyy_xyyz_1[i] + ta1_z_xxxyy_xyyz_0[i] * pa_z[i] - ta1_z_xxxyy_xyyz_1[i] * pc_z[i];

        ta1_z_xxxyyz_xyzz_0[i] = 2.0 * ta1_z_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxyy_xyz_1[i] * fe_0 + ta_xxxyy_xyzz_1[i] + ta1_z_xxxyy_xyzz_0[i] * pa_z[i] - ta1_z_xxxyy_xyzz_1[i] * pc_z[i];

        ta1_z_xxxyyz_xzzz_0[i] = ta1_z_xxxz_xzzz_0[i] * fe_0 - ta1_z_xxxz_xzzz_1[i] * fe_0 + ta1_z_xxxyz_xzzz_0[i] * pa_y[i] - ta1_z_xxxyz_xzzz_1[i] * pc_y[i];

        ta1_z_xxxyyz_yyyy_0[i] = ta_xxxyy_yyyy_1[i] + ta1_z_xxxyy_yyyy_0[i] * pa_z[i] - ta1_z_xxxyy_yyyy_1[i] * pc_z[i];

        ta1_z_xxxyyz_yyyz_0[i] = 2.0 * ta1_z_xyyz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_yyyz_1[i] * fe_0 + ta1_z_xxyyz_yyyz_0[i] * pa_x[i] - ta1_z_xxyyz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxyyz_yyzz_0[i] = 2.0 * ta1_z_xyyz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_yyzz_1[i] * fe_0 + ta1_z_xxyyz_yyzz_0[i] * pa_x[i] - ta1_z_xxyyz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxyyz_yzzz_0[i] = 2.0 * ta1_z_xyyz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_yzzz_1[i] * fe_0 + ta1_z_xxyyz_yzzz_0[i] * pa_x[i] - ta1_z_xxyyz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxyyz_zzzz_0[i] = 2.0 * ta1_z_xyyz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_zzzz_1[i] * fe_0 + ta1_z_xxyyz_zzzz_0[i] * pa_x[i] - ta1_z_xxyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 960-975 components of targeted buffer : IG

    auto ta1_z_xxxyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 960);

    auto ta1_z_xxxyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 961);

    auto ta1_z_xxxyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 962);

    auto ta1_z_xxxyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 963);

    auto ta1_z_xxxyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 964);

    auto ta1_z_xxxyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 965);

    auto ta1_z_xxxyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 966);

    auto ta1_z_xxxyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 967);

    auto ta1_z_xxxyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 968);

    auto ta1_z_xxxyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 969);

    auto ta1_z_xxxyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 970);

    auto ta1_z_xxxyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 971);

    auto ta1_z_xxxyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 972);

    auto ta1_z_xxxyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 973);

    auto ta1_z_xxxyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 974);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxyzz_xxxx_0, ta1_z_xxxyzz_xxxy_0, ta1_z_xxxyzz_xxxz_0, ta1_z_xxxyzz_xxyy_0, ta1_z_xxxyzz_xxyz_0, ta1_z_xxxyzz_xxzz_0, ta1_z_xxxyzz_xyyy_0, ta1_z_xxxyzz_xyyz_0, ta1_z_xxxyzz_xyzz_0, ta1_z_xxxyzz_xzzz_0, ta1_z_xxxyzz_yyyy_0, ta1_z_xxxyzz_yyyz_0, ta1_z_xxxyzz_yyzz_0, ta1_z_xxxyzz_yzzz_0, ta1_z_xxxyzz_zzzz_0, ta1_z_xxxzz_xxx_0, ta1_z_xxxzz_xxx_1, ta1_z_xxxzz_xxxx_0, ta1_z_xxxzz_xxxx_1, ta1_z_xxxzz_xxxy_0, ta1_z_xxxzz_xxxy_1, ta1_z_xxxzz_xxxz_0, ta1_z_xxxzz_xxxz_1, ta1_z_xxxzz_xxy_0, ta1_z_xxxzz_xxy_1, ta1_z_xxxzz_xxyy_0, ta1_z_xxxzz_xxyy_1, ta1_z_xxxzz_xxyz_0, ta1_z_xxxzz_xxyz_1, ta1_z_xxxzz_xxz_0, ta1_z_xxxzz_xxz_1, ta1_z_xxxzz_xxzz_0, ta1_z_xxxzz_xxzz_1, ta1_z_xxxzz_xyy_0, ta1_z_xxxzz_xyy_1, ta1_z_xxxzz_xyyy_0, ta1_z_xxxzz_xyyy_1, ta1_z_xxxzz_xyyz_0, ta1_z_xxxzz_xyyz_1, ta1_z_xxxzz_xyz_0, ta1_z_xxxzz_xyz_1, ta1_z_xxxzz_xyzz_0, ta1_z_xxxzz_xyzz_1, ta1_z_xxxzz_xzz_0, ta1_z_xxxzz_xzz_1, ta1_z_xxxzz_xzzz_0, ta1_z_xxxzz_xzzz_1, ta1_z_xxxzz_zzzz_0, ta1_z_xxxzz_zzzz_1, ta1_z_xxyzz_yyyy_0, ta1_z_xxyzz_yyyy_1, ta1_z_xxyzz_yyyz_0, ta1_z_xxyzz_yyyz_1, ta1_z_xxyzz_yyzz_0, ta1_z_xxyzz_yyzz_1, ta1_z_xxyzz_yzzz_0, ta1_z_xxyzz_yzzz_1, ta1_z_xyzz_yyyy_0, ta1_z_xyzz_yyyy_1, ta1_z_xyzz_yyyz_0, ta1_z_xyzz_yyyz_1, ta1_z_xyzz_yyzz_0, ta1_z_xyzz_yyzz_1, ta1_z_xyzz_yzzz_0, ta1_z_xyzz_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyzz_xxxx_0[i] = ta1_z_xxxzz_xxxx_0[i] * pa_y[i] - ta1_z_xxxzz_xxxx_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxxy_0[i] = ta1_z_xxxzz_xxx_0[i] * fe_0 - ta1_z_xxxzz_xxx_1[i] * fe_0 + ta1_z_xxxzz_xxxy_0[i] * pa_y[i] - ta1_z_xxxzz_xxxy_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxxz_0[i] = ta1_z_xxxzz_xxxz_0[i] * pa_y[i] - ta1_z_xxxzz_xxxz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxyy_0[i] = 2.0 * ta1_z_xxxzz_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxxzz_xxy_1[i] * fe_0 + ta1_z_xxxzz_xxyy_0[i] * pa_y[i] - ta1_z_xxxzz_xxyy_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxyz_0[i] = ta1_z_xxxzz_xxz_0[i] * fe_0 - ta1_z_xxxzz_xxz_1[i] * fe_0 + ta1_z_xxxzz_xxyz_0[i] * pa_y[i] - ta1_z_xxxzz_xxyz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxzz_0[i] = ta1_z_xxxzz_xxzz_0[i] * pa_y[i] - ta1_z_xxxzz_xxzz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xyyy_0[i] = 3.0 * ta1_z_xxxzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxxzz_xyy_1[i] * fe_0 + ta1_z_xxxzz_xyyy_0[i] * pa_y[i] - ta1_z_xxxzz_xyyy_1[i] * pc_y[i];

        ta1_z_xxxyzz_xyyz_0[i] = 2.0 * ta1_z_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxxzz_xyz_1[i] * fe_0 + ta1_z_xxxzz_xyyz_0[i] * pa_y[i] - ta1_z_xxxzz_xyyz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xyzz_0[i] = ta1_z_xxxzz_xzz_0[i] * fe_0 - ta1_z_xxxzz_xzz_1[i] * fe_0 + ta1_z_xxxzz_xyzz_0[i] * pa_y[i] - ta1_z_xxxzz_xyzz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xzzz_0[i] = ta1_z_xxxzz_xzzz_0[i] * pa_y[i] - ta1_z_xxxzz_xzzz_1[i] * pc_y[i];

        ta1_z_xxxyzz_yyyy_0[i] = 2.0 * ta1_z_xyzz_yyyy_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yyyy_1[i] * fe_0 + ta1_z_xxyzz_yyyy_0[i] * pa_x[i] - ta1_z_xxyzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxxyzz_yyyz_0[i] = 2.0 * ta1_z_xyzz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yyyz_1[i] * fe_0 + ta1_z_xxyzz_yyyz_0[i] * pa_x[i] - ta1_z_xxyzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxyzz_yyzz_0[i] = 2.0 * ta1_z_xyzz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yyzz_1[i] * fe_0 + ta1_z_xxyzz_yyzz_0[i] * pa_x[i] - ta1_z_xxyzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxyzz_yzzz_0[i] = 2.0 * ta1_z_xyzz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yzzz_1[i] * fe_0 + ta1_z_xxyzz_yzzz_0[i] * pa_x[i] - ta1_z_xxyzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxyzz_zzzz_0[i] = ta1_z_xxxzz_zzzz_0[i] * pa_y[i] - ta1_z_xxxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 975-990 components of targeted buffer : IG

    auto ta1_z_xxxzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 975);

    auto ta1_z_xxxzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 976);

    auto ta1_z_xxxzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 977);

    auto ta1_z_xxxzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 978);

    auto ta1_z_xxxzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 979);

    auto ta1_z_xxxzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 980);

    auto ta1_z_xxxzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 981);

    auto ta1_z_xxxzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 982);

    auto ta1_z_xxxzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 983);

    auto ta1_z_xxxzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 984);

    auto ta1_z_xxxzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 985);

    auto ta1_z_xxxzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 986);

    auto ta1_z_xxxzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 987);

    auto ta1_z_xxxzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 988);

    auto ta1_z_xxxzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 989);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxz_xxxx_0, ta1_z_xxxz_xxxx_1, ta1_z_xxxz_xxxy_0, ta1_z_xxxz_xxxy_1, ta1_z_xxxz_xxyy_0, ta1_z_xxxz_xxyy_1, ta1_z_xxxz_xyyy_0, ta1_z_xxxz_xyyy_1, ta1_z_xxxzz_xxxx_0, ta1_z_xxxzz_xxxx_1, ta1_z_xxxzz_xxxy_0, ta1_z_xxxzz_xxxy_1, ta1_z_xxxzz_xxyy_0, ta1_z_xxxzz_xxyy_1, ta1_z_xxxzz_xyyy_0, ta1_z_xxxzz_xyyy_1, ta1_z_xxxzzz_xxxx_0, ta1_z_xxxzzz_xxxy_0, ta1_z_xxxzzz_xxxz_0, ta1_z_xxxzzz_xxyy_0, ta1_z_xxxzzz_xxyz_0, ta1_z_xxxzzz_xxzz_0, ta1_z_xxxzzz_xyyy_0, ta1_z_xxxzzz_xyyz_0, ta1_z_xxxzzz_xyzz_0, ta1_z_xxxzzz_xzzz_0, ta1_z_xxxzzz_yyyy_0, ta1_z_xxxzzz_yyyz_0, ta1_z_xxxzzz_yyzz_0, ta1_z_xxxzzz_yzzz_0, ta1_z_xxxzzz_zzzz_0, ta1_z_xxzzz_xxxz_0, ta1_z_xxzzz_xxxz_1, ta1_z_xxzzz_xxyz_0, ta1_z_xxzzz_xxyz_1, ta1_z_xxzzz_xxz_0, ta1_z_xxzzz_xxz_1, ta1_z_xxzzz_xxzz_0, ta1_z_xxzzz_xxzz_1, ta1_z_xxzzz_xyyz_0, ta1_z_xxzzz_xyyz_1, ta1_z_xxzzz_xyz_0, ta1_z_xxzzz_xyz_1, ta1_z_xxzzz_xyzz_0, ta1_z_xxzzz_xyzz_1, ta1_z_xxzzz_xzz_0, ta1_z_xxzzz_xzz_1, ta1_z_xxzzz_xzzz_0, ta1_z_xxzzz_xzzz_1, ta1_z_xxzzz_yyyy_0, ta1_z_xxzzz_yyyy_1, ta1_z_xxzzz_yyyz_0, ta1_z_xxzzz_yyyz_1, ta1_z_xxzzz_yyz_0, ta1_z_xxzzz_yyz_1, ta1_z_xxzzz_yyzz_0, ta1_z_xxzzz_yyzz_1, ta1_z_xxzzz_yzz_0, ta1_z_xxzzz_yzz_1, ta1_z_xxzzz_yzzz_0, ta1_z_xxzzz_yzzz_1, ta1_z_xxzzz_zzz_0, ta1_z_xxzzz_zzz_1, ta1_z_xxzzz_zzzz_0, ta1_z_xxzzz_zzzz_1, ta1_z_xzzz_xxxz_0, ta1_z_xzzz_xxxz_1, ta1_z_xzzz_xxyz_0, ta1_z_xzzz_xxyz_1, ta1_z_xzzz_xxzz_0, ta1_z_xzzz_xxzz_1, ta1_z_xzzz_xyyz_0, ta1_z_xzzz_xyyz_1, ta1_z_xzzz_xyzz_0, ta1_z_xzzz_xyzz_1, ta1_z_xzzz_xzzz_0, ta1_z_xzzz_xzzz_1, ta1_z_xzzz_yyyy_0, ta1_z_xzzz_yyyy_1, ta1_z_xzzz_yyyz_0, ta1_z_xzzz_yyyz_1, ta1_z_xzzz_yyzz_0, ta1_z_xzzz_yyzz_1, ta1_z_xzzz_yzzz_0, ta1_z_xzzz_yzzz_1, ta1_z_xzzz_zzzz_0, ta1_z_xzzz_zzzz_1, ta_xxxzz_xxxx_1, ta_xxxzz_xxxy_1, ta_xxxzz_xxyy_1, ta_xxxzz_xyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzzz_xxxx_0[i] = 2.0 * ta1_z_xxxz_xxxx_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xxxx_1[i] * fe_0 + ta_xxxzz_xxxx_1[i] + ta1_z_xxxzz_xxxx_0[i] * pa_z[i] - ta1_z_xxxzz_xxxx_1[i] * pc_z[i];

        ta1_z_xxxzzz_xxxy_0[i] = 2.0 * ta1_z_xxxz_xxxy_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xxxy_1[i] * fe_0 + ta_xxxzz_xxxy_1[i] + ta1_z_xxxzz_xxxy_0[i] * pa_z[i] - ta1_z_xxxzz_xxxy_1[i] * pc_z[i];

        ta1_z_xxxzzz_xxxz_0[i] = 2.0 * ta1_z_xzzz_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xxzzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_xxzzz_xxz_1[i] * fe_0 + ta1_z_xxzzz_xxxz_0[i] * pa_x[i] - ta1_z_xxzzz_xxxz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xxyy_0[i] = 2.0 * ta1_z_xxxz_xxyy_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xxyy_1[i] * fe_0 + ta_xxxzz_xxyy_1[i] + ta1_z_xxxzz_xxyy_0[i] * pa_z[i] - ta1_z_xxxzz_xxyy_1[i] * pc_z[i];

        ta1_z_xxxzzz_xxyz_0[i] = 2.0 * ta1_z_xzzz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxzzz_xyz_1[i] * fe_0 + ta1_z_xxzzz_xxyz_0[i] * pa_x[i] - ta1_z_xxzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xxzz_0[i] = 2.0 * ta1_z_xzzz_xxzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xxzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_xxzzz_xzz_1[i] * fe_0 + ta1_z_xxzzz_xxzz_0[i] * pa_x[i] - ta1_z_xxzzz_xxzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xyyy_0[i] = 2.0 * ta1_z_xxxz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xyyy_1[i] * fe_0 + ta_xxxzz_xyyy_1[i] + ta1_z_xxxzz_xyyy_0[i] * pa_z[i] - ta1_z_xxxzz_xyyy_1[i] * pc_z[i];

        ta1_z_xxxzzz_xyyz_0[i] = 2.0 * ta1_z_xzzz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xyyz_1[i] * fe_0 + ta1_z_xxzzz_yyz_0[i] * fe_0 - ta1_z_xxzzz_yyz_1[i] * fe_0 + ta1_z_xxzzz_xyyz_0[i] * pa_x[i] - ta1_z_xxzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xyzz_0[i] = 2.0 * ta1_z_xzzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xyzz_1[i] * fe_0 + ta1_z_xxzzz_yzz_0[i] * fe_0 - ta1_z_xxzzz_yzz_1[i] * fe_0 + ta1_z_xxzzz_xyzz_0[i] * pa_x[i] - ta1_z_xxzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xzzz_0[i] = 2.0 * ta1_z_xzzz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xzzz_1[i] * fe_0 + ta1_z_xxzzz_zzz_0[i] * fe_0 - ta1_z_xxzzz_zzz_1[i] * fe_0 + ta1_z_xxzzz_xzzz_0[i] * pa_x[i] - ta1_z_xxzzz_xzzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_yyyy_0[i] = 2.0 * ta1_z_xzzz_yyyy_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yyyy_1[i] * fe_0 + ta1_z_xxzzz_yyyy_0[i] * pa_x[i] - ta1_z_xxzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxxzzz_yyyz_0[i] = 2.0 * ta1_z_xzzz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yyyz_1[i] * fe_0 + ta1_z_xxzzz_yyyz_0[i] * pa_x[i] - ta1_z_xxzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxxzzz_yyzz_0[i] = 2.0 * ta1_z_xzzz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yyzz_1[i] * fe_0 + ta1_z_xxzzz_yyzz_0[i] * pa_x[i] - ta1_z_xxzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_yzzz_0[i] = 2.0 * ta1_z_xzzz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yzzz_1[i] * fe_0 + ta1_z_xxzzz_yzzz_0[i] * pa_x[i] - ta1_z_xxzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_zzzz_0[i] = 2.0 * ta1_z_xzzz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_zzzz_1[i] * fe_0 + ta1_z_xxzzz_zzzz_0[i] * pa_x[i] - ta1_z_xxzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 990-1005 components of targeted buffer : IG

    auto ta1_z_xxyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 990);

    auto ta1_z_xxyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 991);

    auto ta1_z_xxyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 992);

    auto ta1_z_xxyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 993);

    auto ta1_z_xxyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 994);

    auto ta1_z_xxyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 995);

    auto ta1_z_xxyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 996);

    auto ta1_z_xxyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 997);

    auto ta1_z_xxyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 998);

    auto ta1_z_xxyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 999);

    auto ta1_z_xxyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1000);

    auto ta1_z_xxyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1001);

    auto ta1_z_xxyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1002);

    auto ta1_z_xxyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1003);

    auto ta1_z_xxyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1004);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxyy_xxxx_0, ta1_z_xxyy_xxxx_1, ta1_z_xxyy_xxxz_0, ta1_z_xxyy_xxxz_1, ta1_z_xxyy_xxzz_0, ta1_z_xxyy_xxzz_1, ta1_z_xxyy_xzzz_0, ta1_z_xxyy_xzzz_1, ta1_z_xxyyy_xxxx_0, ta1_z_xxyyy_xxxx_1, ta1_z_xxyyy_xxxz_0, ta1_z_xxyyy_xxxz_1, ta1_z_xxyyy_xxzz_0, ta1_z_xxyyy_xxzz_1, ta1_z_xxyyy_xzzz_0, ta1_z_xxyyy_xzzz_1, ta1_z_xxyyyy_xxxx_0, ta1_z_xxyyyy_xxxy_0, ta1_z_xxyyyy_xxxz_0, ta1_z_xxyyyy_xxyy_0, ta1_z_xxyyyy_xxyz_0, ta1_z_xxyyyy_xxzz_0, ta1_z_xxyyyy_xyyy_0, ta1_z_xxyyyy_xyyz_0, ta1_z_xxyyyy_xyzz_0, ta1_z_xxyyyy_xzzz_0, ta1_z_xxyyyy_yyyy_0, ta1_z_xxyyyy_yyyz_0, ta1_z_xxyyyy_yyzz_0, ta1_z_xxyyyy_yzzz_0, ta1_z_xxyyyy_zzzz_0, ta1_z_xyyyy_xxxy_0, ta1_z_xyyyy_xxxy_1, ta1_z_xyyyy_xxy_0, ta1_z_xyyyy_xxy_1, ta1_z_xyyyy_xxyy_0, ta1_z_xyyyy_xxyy_1, ta1_z_xyyyy_xxyz_0, ta1_z_xyyyy_xxyz_1, ta1_z_xyyyy_xyy_0, ta1_z_xyyyy_xyy_1, ta1_z_xyyyy_xyyy_0, ta1_z_xyyyy_xyyy_1, ta1_z_xyyyy_xyyz_0, ta1_z_xyyyy_xyyz_1, ta1_z_xyyyy_xyz_0, ta1_z_xyyyy_xyz_1, ta1_z_xyyyy_xyzz_0, ta1_z_xyyyy_xyzz_1, ta1_z_xyyyy_yyy_0, ta1_z_xyyyy_yyy_1, ta1_z_xyyyy_yyyy_0, ta1_z_xyyyy_yyyy_1, ta1_z_xyyyy_yyyz_0, ta1_z_xyyyy_yyyz_1, ta1_z_xyyyy_yyz_0, ta1_z_xyyyy_yyz_1, ta1_z_xyyyy_yyzz_0, ta1_z_xyyyy_yyzz_1, ta1_z_xyyyy_yzz_0, ta1_z_xyyyy_yzz_1, ta1_z_xyyyy_yzzz_0, ta1_z_xyyyy_yzzz_1, ta1_z_xyyyy_zzzz_0, ta1_z_xyyyy_zzzz_1, ta1_z_yyyy_xxxy_0, ta1_z_yyyy_xxxy_1, ta1_z_yyyy_xxyy_0, ta1_z_yyyy_xxyy_1, ta1_z_yyyy_xxyz_0, ta1_z_yyyy_xxyz_1, ta1_z_yyyy_xyyy_0, ta1_z_yyyy_xyyy_1, ta1_z_yyyy_xyyz_0, ta1_z_yyyy_xyyz_1, ta1_z_yyyy_xyzz_0, ta1_z_yyyy_xyzz_1, ta1_z_yyyy_yyyy_0, ta1_z_yyyy_yyyy_1, ta1_z_yyyy_yyyz_0, ta1_z_yyyy_yyyz_1, ta1_z_yyyy_yyzz_0, ta1_z_yyyy_yyzz_1, ta1_z_yyyy_yzzz_0, ta1_z_yyyy_yzzz_1, ta1_z_yyyy_zzzz_0, ta1_z_yyyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyy_xxxx_0[i] = 3.0 * ta1_z_xxyy_xxxx_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxxx_1[i] * fe_0 + ta1_z_xxyyy_xxxx_0[i] * pa_y[i] - ta1_z_xxyyy_xxxx_1[i] * pc_y[i];

        ta1_z_xxyyyy_xxxy_0[i] = ta1_z_yyyy_xxxy_0[i] * fe_0 - ta1_z_yyyy_xxxy_1[i] * fe_0 + 3.0 * ta1_z_xyyyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_xyyyy_xxy_1[i] * fe_0 + ta1_z_xyyyy_xxxy_0[i] * pa_x[i] - ta1_z_xyyyy_xxxy_1[i] * pc_x[i];

        ta1_z_xxyyyy_xxxz_0[i] = 3.0 * ta1_z_xxyy_xxxz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxxz_1[i] * fe_0 + ta1_z_xxyyy_xxxz_0[i] * pa_y[i] - ta1_z_xxyyy_xxxz_1[i] * pc_y[i];

        ta1_z_xxyyyy_xxyy_0[i] = ta1_z_yyyy_xxyy_0[i] * fe_0 - ta1_z_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xyyyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_xyyyy_xyy_1[i] * fe_0 + ta1_z_xyyyy_xxyy_0[i] * pa_x[i] - ta1_z_xyyyy_xxyy_1[i] * pc_x[i];

        ta1_z_xxyyyy_xxyz_0[i] = ta1_z_yyyy_xxyz_0[i] * fe_0 - ta1_z_yyyy_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xyyyy_xyz_1[i] * fe_0 + ta1_z_xyyyy_xxyz_0[i] * pa_x[i] - ta1_z_xyyyy_xxyz_1[i] * pc_x[i];

        ta1_z_xxyyyy_xxzz_0[i] = 3.0 * ta1_z_xxyy_xxzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxzz_1[i] * fe_0 + ta1_z_xxyyy_xxzz_0[i] * pa_y[i] - ta1_z_xxyyy_xxzz_1[i] * pc_y[i];

        ta1_z_xxyyyy_xyyy_0[i] = ta1_z_yyyy_xyyy_0[i] * fe_0 - ta1_z_yyyy_xyyy_1[i] * fe_0 + ta1_z_xyyyy_yyy_0[i] * fe_0 - ta1_z_xyyyy_yyy_1[i] * fe_0 + ta1_z_xyyyy_xyyy_0[i] * pa_x[i] - ta1_z_xyyyy_xyyy_1[i] * pc_x[i];

        ta1_z_xxyyyy_xyyz_0[i] = ta1_z_yyyy_xyyz_0[i] * fe_0 - ta1_z_yyyy_xyyz_1[i] * fe_0 + ta1_z_xyyyy_yyz_0[i] * fe_0 - ta1_z_xyyyy_yyz_1[i] * fe_0 + ta1_z_xyyyy_xyyz_0[i] * pa_x[i] - ta1_z_xyyyy_xyyz_1[i] * pc_x[i];

        ta1_z_xxyyyy_xyzz_0[i] = ta1_z_yyyy_xyzz_0[i] * fe_0 - ta1_z_yyyy_xyzz_1[i] * fe_0 + ta1_z_xyyyy_yzz_0[i] * fe_0 - ta1_z_xyyyy_yzz_1[i] * fe_0 + ta1_z_xyyyy_xyzz_0[i] * pa_x[i] - ta1_z_xyyyy_xyzz_1[i] * pc_x[i];

        ta1_z_xxyyyy_xzzz_0[i] = 3.0 * ta1_z_xxyy_xzzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xzzz_1[i] * fe_0 + ta1_z_xxyyy_xzzz_0[i] * pa_y[i] - ta1_z_xxyyy_xzzz_1[i] * pc_y[i];

        ta1_z_xxyyyy_yyyy_0[i] = ta1_z_yyyy_yyyy_0[i] * fe_0 - ta1_z_yyyy_yyyy_1[i] * fe_0 + ta1_z_xyyyy_yyyy_0[i] * pa_x[i] - ta1_z_xyyyy_yyyy_1[i] * pc_x[i];

        ta1_z_xxyyyy_yyyz_0[i] = ta1_z_yyyy_yyyz_0[i] * fe_0 - ta1_z_yyyy_yyyz_1[i] * fe_0 + ta1_z_xyyyy_yyyz_0[i] * pa_x[i] - ta1_z_xyyyy_yyyz_1[i] * pc_x[i];

        ta1_z_xxyyyy_yyzz_0[i] = ta1_z_yyyy_yyzz_0[i] * fe_0 - ta1_z_yyyy_yyzz_1[i] * fe_0 + ta1_z_xyyyy_yyzz_0[i] * pa_x[i] - ta1_z_xyyyy_yyzz_1[i] * pc_x[i];

        ta1_z_xxyyyy_yzzz_0[i] = ta1_z_yyyy_yzzz_0[i] * fe_0 - ta1_z_yyyy_yzzz_1[i] * fe_0 + ta1_z_xyyyy_yzzz_0[i] * pa_x[i] - ta1_z_xyyyy_yzzz_1[i] * pc_x[i];

        ta1_z_xxyyyy_zzzz_0[i] = ta1_z_yyyy_zzzz_0[i] * fe_0 - ta1_z_yyyy_zzzz_1[i] * fe_0 + ta1_z_xyyyy_zzzz_0[i] * pa_x[i] - ta1_z_xyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 1005-1020 components of targeted buffer : IG

    auto ta1_z_xxyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1005);

    auto ta1_z_xxyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1006);

    auto ta1_z_xxyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1007);

    auto ta1_z_xxyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1008);

    auto ta1_z_xxyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1009);

    auto ta1_z_xxyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1010);

    auto ta1_z_xxyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1011);

    auto ta1_z_xxyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1012);

    auto ta1_z_xxyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1013);

    auto ta1_z_xxyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1014);

    auto ta1_z_xxyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1015);

    auto ta1_z_xxyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1016);

    auto ta1_z_xxyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1017);

    auto ta1_z_xxyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1018);

    auto ta1_z_xxyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1019);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxyyy_xxxx_0, ta1_z_xxyyy_xxxx_1, ta1_z_xxyyy_xxxy_0, ta1_z_xxyyy_xxxy_1, ta1_z_xxyyy_xxy_0, ta1_z_xxyyy_xxy_1, ta1_z_xxyyy_xxyy_0, ta1_z_xxyyy_xxyy_1, ta1_z_xxyyy_xxyz_0, ta1_z_xxyyy_xxyz_1, ta1_z_xxyyy_xyy_0, ta1_z_xxyyy_xyy_1, ta1_z_xxyyy_xyyy_0, ta1_z_xxyyy_xyyy_1, ta1_z_xxyyy_xyyz_0, ta1_z_xxyyy_xyyz_1, ta1_z_xxyyy_xyz_0, ta1_z_xxyyy_xyz_1, ta1_z_xxyyy_xyzz_0, ta1_z_xxyyy_xyzz_1, ta1_z_xxyyy_yyyy_0, ta1_z_xxyyy_yyyy_1, ta1_z_xxyyyz_xxxx_0, ta1_z_xxyyyz_xxxy_0, ta1_z_xxyyyz_xxxz_0, ta1_z_xxyyyz_xxyy_0, ta1_z_xxyyyz_xxyz_0, ta1_z_xxyyyz_xxzz_0, ta1_z_xxyyyz_xyyy_0, ta1_z_xxyyyz_xyyz_0, ta1_z_xxyyyz_xyzz_0, ta1_z_xxyyyz_xzzz_0, ta1_z_xxyyyz_yyyy_0, ta1_z_xxyyyz_yyyz_0, ta1_z_xxyyyz_yyzz_0, ta1_z_xxyyyz_yzzz_0, ta1_z_xxyyyz_zzzz_0, ta1_z_xxyyz_xxxz_0, ta1_z_xxyyz_xxxz_1, ta1_z_xxyyz_xxzz_0, ta1_z_xxyyz_xxzz_1, ta1_z_xxyyz_xzzz_0, ta1_z_xxyyz_xzzz_1, ta1_z_xxyz_xxxz_0, ta1_z_xxyz_xxxz_1, ta1_z_xxyz_xxzz_0, ta1_z_xxyz_xxzz_1, ta1_z_xxyz_xzzz_0, ta1_z_xxyz_xzzz_1, ta1_z_xyyyz_yyyz_0, ta1_z_xyyyz_yyyz_1, ta1_z_xyyyz_yyzz_0, ta1_z_xyyyz_yyzz_1, ta1_z_xyyyz_yzzz_0, ta1_z_xyyyz_yzzz_1, ta1_z_xyyyz_zzzz_0, ta1_z_xyyyz_zzzz_1, ta1_z_yyyz_yyyz_0, ta1_z_yyyz_yyyz_1, ta1_z_yyyz_yyzz_0, ta1_z_yyyz_yyzz_1, ta1_z_yyyz_yzzz_0, ta1_z_yyyz_yzzz_1, ta1_z_yyyz_zzzz_0, ta1_z_yyyz_zzzz_1, ta_xxyyy_xxxx_1, ta_xxyyy_xxxy_1, ta_xxyyy_xxyy_1, ta_xxyyy_xxyz_1, ta_xxyyy_xyyy_1, ta_xxyyy_xyyz_1, ta_xxyyy_xyzz_1, ta_xxyyy_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyz_xxxx_0[i] = ta_xxyyy_xxxx_1[i] + ta1_z_xxyyy_xxxx_0[i] * pa_z[i] - ta1_z_xxyyy_xxxx_1[i] * pc_z[i];

        ta1_z_xxyyyz_xxxy_0[i] = ta_xxyyy_xxxy_1[i] + ta1_z_xxyyy_xxxy_0[i] * pa_z[i] - ta1_z_xxyyy_xxxy_1[i] * pc_z[i];

        ta1_z_xxyyyz_xxxz_0[i] = 2.0 * ta1_z_xxyz_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xxyz_xxxz_1[i] * fe_0 + ta1_z_xxyyz_xxxz_0[i] * pa_y[i] - ta1_z_xxyyz_xxxz_1[i] * pc_y[i];

        ta1_z_xxyyyz_xxyy_0[i] = ta_xxyyy_xxyy_1[i] + ta1_z_xxyyy_xxyy_0[i] * pa_z[i] - ta1_z_xxyyy_xxyy_1[i] * pc_z[i];

        ta1_z_xxyyyz_xxyz_0[i] = ta1_z_xxyyy_xxy_0[i] * fe_0 - ta1_z_xxyyy_xxy_1[i] * fe_0 + ta_xxyyy_xxyz_1[i] + ta1_z_xxyyy_xxyz_0[i] * pa_z[i] - ta1_z_xxyyy_xxyz_1[i] * pc_z[i];

        ta1_z_xxyyyz_xxzz_0[i] = 2.0 * ta1_z_xxyz_xxzz_0[i] * fe_0 - 2.0 * ta1_z_xxyz_xxzz_1[i] * fe_0 + ta1_z_xxyyz_xxzz_0[i] * pa_y[i] - ta1_z_xxyyz_xxzz_1[i] * pc_y[i];

        ta1_z_xxyyyz_xyyy_0[i] = ta_xxyyy_xyyy_1[i] + ta1_z_xxyyy_xyyy_0[i] * pa_z[i] - ta1_z_xxyyy_xyyy_1[i] * pc_z[i];

        ta1_z_xxyyyz_xyyz_0[i] = ta1_z_xxyyy_xyy_0[i] * fe_0 - ta1_z_xxyyy_xyy_1[i] * fe_0 + ta_xxyyy_xyyz_1[i] + ta1_z_xxyyy_xyyz_0[i] * pa_z[i] - ta1_z_xxyyy_xyyz_1[i] * pc_z[i];

        ta1_z_xxyyyz_xyzz_0[i] = 2.0 * ta1_z_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxyyy_xyz_1[i] * fe_0 + ta_xxyyy_xyzz_1[i] + ta1_z_xxyyy_xyzz_0[i] * pa_z[i] - ta1_z_xxyyy_xyzz_1[i] * pc_z[i];

        ta1_z_xxyyyz_xzzz_0[i] = 2.0 * ta1_z_xxyz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_xxyz_xzzz_1[i] * fe_0 + ta1_z_xxyyz_xzzz_0[i] * pa_y[i] - ta1_z_xxyyz_xzzz_1[i] * pc_y[i];

        ta1_z_xxyyyz_yyyy_0[i] = ta_xxyyy_yyyy_1[i] + ta1_z_xxyyy_yyyy_0[i] * pa_z[i] - ta1_z_xxyyy_yyyy_1[i] * pc_z[i];

        ta1_z_xxyyyz_yyyz_0[i] = ta1_z_yyyz_yyyz_0[i] * fe_0 - ta1_z_yyyz_yyyz_1[i] * fe_0 + ta1_z_xyyyz_yyyz_0[i] * pa_x[i] - ta1_z_xyyyz_yyyz_1[i] * pc_x[i];

        ta1_z_xxyyyz_yyzz_0[i] = ta1_z_yyyz_yyzz_0[i] * fe_0 - ta1_z_yyyz_yyzz_1[i] * fe_0 + ta1_z_xyyyz_yyzz_0[i] * pa_x[i] - ta1_z_xyyyz_yyzz_1[i] * pc_x[i];

        ta1_z_xxyyyz_yzzz_0[i] = ta1_z_yyyz_yzzz_0[i] * fe_0 - ta1_z_yyyz_yzzz_1[i] * fe_0 + ta1_z_xyyyz_yzzz_0[i] * pa_x[i] - ta1_z_xyyyz_yzzz_1[i] * pc_x[i];

        ta1_z_xxyyyz_zzzz_0[i] = ta1_z_yyyz_zzzz_0[i] * fe_0 - ta1_z_yyyz_zzzz_1[i] * fe_0 + ta1_z_xyyyz_zzzz_0[i] * pa_x[i] - ta1_z_xyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1020-1035 components of targeted buffer : IG

    auto ta1_z_xxyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1020);

    auto ta1_z_xxyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1021);

    auto ta1_z_xxyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1022);

    auto ta1_z_xxyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1023);

    auto ta1_z_xxyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1024);

    auto ta1_z_xxyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1025);

    auto ta1_z_xxyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1026);

    auto ta1_z_xxyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1027);

    auto ta1_z_xxyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1028);

    auto ta1_z_xxyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1029);

    auto ta1_z_xxyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1030);

    auto ta1_z_xxyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1031);

    auto ta1_z_xxyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1032);

    auto ta1_z_xxyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1033);

    auto ta1_z_xxyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1034);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxyy_xxxy_0, ta1_z_xxyy_xxxy_1, ta1_z_xxyy_xxyy_0, ta1_z_xxyy_xxyy_1, ta1_z_xxyy_xyyy_0, ta1_z_xxyy_xyyy_1, ta1_z_xxyyz_xxxy_0, ta1_z_xxyyz_xxxy_1, ta1_z_xxyyz_xxyy_0, ta1_z_xxyyz_xxyy_1, ta1_z_xxyyz_xyyy_0, ta1_z_xxyyz_xyyy_1, ta1_z_xxyyzz_xxxx_0, ta1_z_xxyyzz_xxxy_0, ta1_z_xxyyzz_xxxz_0, ta1_z_xxyyzz_xxyy_0, ta1_z_xxyyzz_xxyz_0, ta1_z_xxyyzz_xxzz_0, ta1_z_xxyyzz_xyyy_0, ta1_z_xxyyzz_xyyz_0, ta1_z_xxyyzz_xyzz_0, ta1_z_xxyyzz_xzzz_0, ta1_z_xxyyzz_yyyy_0, ta1_z_xxyyzz_yyyz_0, ta1_z_xxyyzz_yyzz_0, ta1_z_xxyyzz_yzzz_0, ta1_z_xxyyzz_zzzz_0, ta1_z_xxyzz_xxxx_0, ta1_z_xxyzz_xxxx_1, ta1_z_xxyzz_xxxz_0, ta1_z_xxyzz_xxxz_1, ta1_z_xxyzz_xxzz_0, ta1_z_xxyzz_xxzz_1, ta1_z_xxyzz_xzzz_0, ta1_z_xxyzz_xzzz_1, ta1_z_xxzz_xxxx_0, ta1_z_xxzz_xxxx_1, ta1_z_xxzz_xxxz_0, ta1_z_xxzz_xxxz_1, ta1_z_xxzz_xxzz_0, ta1_z_xxzz_xxzz_1, ta1_z_xxzz_xzzz_0, ta1_z_xxzz_xzzz_1, ta1_z_xyyzz_xxyz_0, ta1_z_xyyzz_xxyz_1, ta1_z_xyyzz_xyyz_0, ta1_z_xyyzz_xyyz_1, ta1_z_xyyzz_xyz_0, ta1_z_xyyzz_xyz_1, ta1_z_xyyzz_xyzz_0, ta1_z_xyyzz_xyzz_1, ta1_z_xyyzz_yyyy_0, ta1_z_xyyzz_yyyy_1, ta1_z_xyyzz_yyyz_0, ta1_z_xyyzz_yyyz_1, ta1_z_xyyzz_yyz_0, ta1_z_xyyzz_yyz_1, ta1_z_xyyzz_yyzz_0, ta1_z_xyyzz_yyzz_1, ta1_z_xyyzz_yzz_0, ta1_z_xyyzz_yzz_1, ta1_z_xyyzz_yzzz_0, ta1_z_xyyzz_yzzz_1, ta1_z_xyyzz_zzzz_0, ta1_z_xyyzz_zzzz_1, ta1_z_yyzz_xxyz_0, ta1_z_yyzz_xxyz_1, ta1_z_yyzz_xyyz_0, ta1_z_yyzz_xyyz_1, ta1_z_yyzz_xyzz_0, ta1_z_yyzz_xyzz_1, ta1_z_yyzz_yyyy_0, ta1_z_yyzz_yyyy_1, ta1_z_yyzz_yyyz_0, ta1_z_yyzz_yyyz_1, ta1_z_yyzz_yyzz_0, ta1_z_yyzz_yyzz_1, ta1_z_yyzz_yzzz_0, ta1_z_yyzz_yzzz_1, ta1_z_yyzz_zzzz_0, ta1_z_yyzz_zzzz_1, ta_xxyyz_xxxy_1, ta_xxyyz_xxyy_1, ta_xxyyz_xyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyzz_xxxx_0[i] = ta1_z_xxzz_xxxx_0[i] * fe_0 - ta1_z_xxzz_xxxx_1[i] * fe_0 + ta1_z_xxyzz_xxxx_0[i] * pa_y[i] - ta1_z_xxyzz_xxxx_1[i] * pc_y[i];

        ta1_z_xxyyzz_xxxy_0[i] = ta1_z_xxyy_xxxy_0[i] * fe_0 - ta1_z_xxyy_xxxy_1[i] * fe_0 + ta_xxyyz_xxxy_1[i] + ta1_z_xxyyz_xxxy_0[i] * pa_z[i] - ta1_z_xxyyz_xxxy_1[i] * pc_z[i];

        ta1_z_xxyyzz_xxxz_0[i] = ta1_z_xxzz_xxxz_0[i] * fe_0 - ta1_z_xxzz_xxxz_1[i] * fe_0 + ta1_z_xxyzz_xxxz_0[i] * pa_y[i] - ta1_z_xxyzz_xxxz_1[i] * pc_y[i];

        ta1_z_xxyyzz_xxyy_0[i] = ta1_z_xxyy_xxyy_0[i] * fe_0 - ta1_z_xxyy_xxyy_1[i] * fe_0 + ta_xxyyz_xxyy_1[i] + ta1_z_xxyyz_xxyy_0[i] * pa_z[i] - ta1_z_xxyyz_xxyy_1[i] * pc_z[i];

        ta1_z_xxyyzz_xxyz_0[i] = ta1_z_yyzz_xxyz_0[i] * fe_0 - ta1_z_yyzz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xyyzz_xyz_1[i] * fe_0 + ta1_z_xyyzz_xxyz_0[i] * pa_x[i] - ta1_z_xyyzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxyyzz_xxzz_0[i] = ta1_z_xxzz_xxzz_0[i] * fe_0 - ta1_z_xxzz_xxzz_1[i] * fe_0 + ta1_z_xxyzz_xxzz_0[i] * pa_y[i] - ta1_z_xxyzz_xxzz_1[i] * pc_y[i];

        ta1_z_xxyyzz_xyyy_0[i] = ta1_z_xxyy_xyyy_0[i] * fe_0 - ta1_z_xxyy_xyyy_1[i] * fe_0 + ta_xxyyz_xyyy_1[i] + ta1_z_xxyyz_xyyy_0[i] * pa_z[i] - ta1_z_xxyyz_xyyy_1[i] * pc_z[i];

        ta1_z_xxyyzz_xyyz_0[i] = ta1_z_yyzz_xyyz_0[i] * fe_0 - ta1_z_yyzz_xyyz_1[i] * fe_0 + ta1_z_xyyzz_yyz_0[i] * fe_0 - ta1_z_xyyzz_yyz_1[i] * fe_0 + ta1_z_xyyzz_xyyz_0[i] * pa_x[i] - ta1_z_xyyzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxyyzz_xyzz_0[i] = ta1_z_yyzz_xyzz_0[i] * fe_0 - ta1_z_yyzz_xyzz_1[i] * fe_0 + ta1_z_xyyzz_yzz_0[i] * fe_0 - ta1_z_xyyzz_yzz_1[i] * fe_0 + ta1_z_xyyzz_xyzz_0[i] * pa_x[i] - ta1_z_xyyzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxyyzz_xzzz_0[i] = ta1_z_xxzz_xzzz_0[i] * fe_0 - ta1_z_xxzz_xzzz_1[i] * fe_0 + ta1_z_xxyzz_xzzz_0[i] * pa_y[i] - ta1_z_xxyzz_xzzz_1[i] * pc_y[i];

        ta1_z_xxyyzz_yyyy_0[i] = ta1_z_yyzz_yyyy_0[i] * fe_0 - ta1_z_yyzz_yyyy_1[i] * fe_0 + ta1_z_xyyzz_yyyy_0[i] * pa_x[i] - ta1_z_xyyzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxyyzz_yyyz_0[i] = ta1_z_yyzz_yyyz_0[i] * fe_0 - ta1_z_yyzz_yyyz_1[i] * fe_0 + ta1_z_xyyzz_yyyz_0[i] * pa_x[i] - ta1_z_xyyzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxyyzz_yyzz_0[i] = ta1_z_yyzz_yyzz_0[i] * fe_0 - ta1_z_yyzz_yyzz_1[i] * fe_0 + ta1_z_xyyzz_yyzz_0[i] * pa_x[i] - ta1_z_xyyzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxyyzz_yzzz_0[i] = ta1_z_yyzz_yzzz_0[i] * fe_0 - ta1_z_yyzz_yzzz_1[i] * fe_0 + ta1_z_xyyzz_yzzz_0[i] * pa_x[i] - ta1_z_xyyzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxyyzz_zzzz_0[i] = ta1_z_yyzz_zzzz_0[i] * fe_0 - ta1_z_yyzz_zzzz_1[i] * fe_0 + ta1_z_xyyzz_zzzz_0[i] * pa_x[i] - ta1_z_xyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1035-1050 components of targeted buffer : IG

    auto ta1_z_xxyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1035);

    auto ta1_z_xxyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1036);

    auto ta1_z_xxyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1037);

    auto ta1_z_xxyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1038);

    auto ta1_z_xxyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1039);

    auto ta1_z_xxyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1040);

    auto ta1_z_xxyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1041);

    auto ta1_z_xxyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1042);

    auto ta1_z_xxyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1043);

    auto ta1_z_xxyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1044);

    auto ta1_z_xxyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1045);

    auto ta1_z_xxyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1046);

    auto ta1_z_xxyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1047);

    auto ta1_z_xxyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1048);

    auto ta1_z_xxyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1049);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxyzzz_xxxx_0, ta1_z_xxyzzz_xxxy_0, ta1_z_xxyzzz_xxxz_0, ta1_z_xxyzzz_xxyy_0, ta1_z_xxyzzz_xxyz_0, ta1_z_xxyzzz_xxzz_0, ta1_z_xxyzzz_xyyy_0, ta1_z_xxyzzz_xyyz_0, ta1_z_xxyzzz_xyzz_0, ta1_z_xxyzzz_xzzz_0, ta1_z_xxyzzz_yyyy_0, ta1_z_xxyzzz_yyyz_0, ta1_z_xxyzzz_yyzz_0, ta1_z_xxyzzz_yzzz_0, ta1_z_xxyzzz_zzzz_0, ta1_z_xxzzz_xxx_0, ta1_z_xxzzz_xxx_1, ta1_z_xxzzz_xxxx_0, ta1_z_xxzzz_xxxx_1, ta1_z_xxzzz_xxxy_0, ta1_z_xxzzz_xxxy_1, ta1_z_xxzzz_xxxz_0, ta1_z_xxzzz_xxxz_1, ta1_z_xxzzz_xxy_0, ta1_z_xxzzz_xxy_1, ta1_z_xxzzz_xxyy_0, ta1_z_xxzzz_xxyy_1, ta1_z_xxzzz_xxyz_0, ta1_z_xxzzz_xxyz_1, ta1_z_xxzzz_xxz_0, ta1_z_xxzzz_xxz_1, ta1_z_xxzzz_xxzz_0, ta1_z_xxzzz_xxzz_1, ta1_z_xxzzz_xyy_0, ta1_z_xxzzz_xyy_1, ta1_z_xxzzz_xyyy_0, ta1_z_xxzzz_xyyy_1, ta1_z_xxzzz_xyyz_0, ta1_z_xxzzz_xyyz_1, ta1_z_xxzzz_xyz_0, ta1_z_xxzzz_xyz_1, ta1_z_xxzzz_xyzz_0, ta1_z_xxzzz_xyzz_1, ta1_z_xxzzz_xzz_0, ta1_z_xxzzz_xzz_1, ta1_z_xxzzz_xzzz_0, ta1_z_xxzzz_xzzz_1, ta1_z_xxzzz_zzzz_0, ta1_z_xxzzz_zzzz_1, ta1_z_xyzzz_yyyy_0, ta1_z_xyzzz_yyyy_1, ta1_z_xyzzz_yyyz_0, ta1_z_xyzzz_yyyz_1, ta1_z_xyzzz_yyzz_0, ta1_z_xyzzz_yyzz_1, ta1_z_xyzzz_yzzz_0, ta1_z_xyzzz_yzzz_1, ta1_z_yzzz_yyyy_0, ta1_z_yzzz_yyyy_1, ta1_z_yzzz_yyyz_0, ta1_z_yzzz_yyyz_1, ta1_z_yzzz_yyzz_0, ta1_z_yzzz_yyzz_1, ta1_z_yzzz_yzzz_0, ta1_z_yzzz_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzzz_xxxx_0[i] = ta1_z_xxzzz_xxxx_0[i] * pa_y[i] - ta1_z_xxzzz_xxxx_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxxy_0[i] = ta1_z_xxzzz_xxx_0[i] * fe_0 - ta1_z_xxzzz_xxx_1[i] * fe_0 + ta1_z_xxzzz_xxxy_0[i] * pa_y[i] - ta1_z_xxzzz_xxxy_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxxz_0[i] = ta1_z_xxzzz_xxxz_0[i] * pa_y[i] - ta1_z_xxzzz_xxxz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxyy_0[i] = 2.0 * ta1_z_xxzzz_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxzzz_xxy_1[i] * fe_0 + ta1_z_xxzzz_xxyy_0[i] * pa_y[i] - ta1_z_xxzzz_xxyy_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxyz_0[i] = ta1_z_xxzzz_xxz_0[i] * fe_0 - ta1_z_xxzzz_xxz_1[i] * fe_0 + ta1_z_xxzzz_xxyz_0[i] * pa_y[i] - ta1_z_xxzzz_xxyz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxzz_0[i] = ta1_z_xxzzz_xxzz_0[i] * pa_y[i] - ta1_z_xxzzz_xxzz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xyyy_0[i] = 3.0 * ta1_z_xxzzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxzzz_xyy_1[i] * fe_0 + ta1_z_xxzzz_xyyy_0[i] * pa_y[i] - ta1_z_xxzzz_xyyy_1[i] * pc_y[i];

        ta1_z_xxyzzz_xyyz_0[i] = 2.0 * ta1_z_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xxzzz_xyz_1[i] * fe_0 + ta1_z_xxzzz_xyyz_0[i] * pa_y[i] - ta1_z_xxzzz_xyyz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xyzz_0[i] = ta1_z_xxzzz_xzz_0[i] * fe_0 - ta1_z_xxzzz_xzz_1[i] * fe_0 + ta1_z_xxzzz_xyzz_0[i] * pa_y[i] - ta1_z_xxzzz_xyzz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xzzz_0[i] = ta1_z_xxzzz_xzzz_0[i] * pa_y[i] - ta1_z_xxzzz_xzzz_1[i] * pc_y[i];

        ta1_z_xxyzzz_yyyy_0[i] = ta1_z_yzzz_yyyy_0[i] * fe_0 - ta1_z_yzzz_yyyy_1[i] * fe_0 + ta1_z_xyzzz_yyyy_0[i] * pa_x[i] - ta1_z_xyzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxyzzz_yyyz_0[i] = ta1_z_yzzz_yyyz_0[i] * fe_0 - ta1_z_yzzz_yyyz_1[i] * fe_0 + ta1_z_xyzzz_yyyz_0[i] * pa_x[i] - ta1_z_xyzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxyzzz_yyzz_0[i] = ta1_z_yzzz_yyzz_0[i] * fe_0 - ta1_z_yzzz_yyzz_1[i] * fe_0 + ta1_z_xyzzz_yyzz_0[i] * pa_x[i] - ta1_z_xyzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxyzzz_yzzz_0[i] = ta1_z_yzzz_yzzz_0[i] * fe_0 - ta1_z_yzzz_yzzz_1[i] * fe_0 + ta1_z_xyzzz_yzzz_0[i] * pa_x[i] - ta1_z_xyzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxyzzz_zzzz_0[i] = ta1_z_xxzzz_zzzz_0[i] * pa_y[i] - ta1_z_xxzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1050-1065 components of targeted buffer : IG

    auto ta1_z_xxzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1050);

    auto ta1_z_xxzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1051);

    auto ta1_z_xxzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1052);

    auto ta1_z_xxzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1053);

    auto ta1_z_xxzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1054);

    auto ta1_z_xxzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1055);

    auto ta1_z_xxzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1056);

    auto ta1_z_xxzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1057);

    auto ta1_z_xxzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1058);

    auto ta1_z_xxzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1059);

    auto ta1_z_xxzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1060);

    auto ta1_z_xxzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1061);

    auto ta1_z_xxzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1062);

    auto ta1_z_xxzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1063);

    auto ta1_z_xxzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1064);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxzz_xxxx_0, ta1_z_xxzz_xxxx_1, ta1_z_xxzz_xxxy_0, ta1_z_xxzz_xxxy_1, ta1_z_xxzz_xxyy_0, ta1_z_xxzz_xxyy_1, ta1_z_xxzz_xyyy_0, ta1_z_xxzz_xyyy_1, ta1_z_xxzzz_xxxx_0, ta1_z_xxzzz_xxxx_1, ta1_z_xxzzz_xxxy_0, ta1_z_xxzzz_xxxy_1, ta1_z_xxzzz_xxyy_0, ta1_z_xxzzz_xxyy_1, ta1_z_xxzzz_xyyy_0, ta1_z_xxzzz_xyyy_1, ta1_z_xxzzzz_xxxx_0, ta1_z_xxzzzz_xxxy_0, ta1_z_xxzzzz_xxxz_0, ta1_z_xxzzzz_xxyy_0, ta1_z_xxzzzz_xxyz_0, ta1_z_xxzzzz_xxzz_0, ta1_z_xxzzzz_xyyy_0, ta1_z_xxzzzz_xyyz_0, ta1_z_xxzzzz_xyzz_0, ta1_z_xxzzzz_xzzz_0, ta1_z_xxzzzz_yyyy_0, ta1_z_xxzzzz_yyyz_0, ta1_z_xxzzzz_yyzz_0, ta1_z_xxzzzz_yzzz_0, ta1_z_xxzzzz_zzzz_0, ta1_z_xzzzz_xxxz_0, ta1_z_xzzzz_xxxz_1, ta1_z_xzzzz_xxyz_0, ta1_z_xzzzz_xxyz_1, ta1_z_xzzzz_xxz_0, ta1_z_xzzzz_xxz_1, ta1_z_xzzzz_xxzz_0, ta1_z_xzzzz_xxzz_1, ta1_z_xzzzz_xyyz_0, ta1_z_xzzzz_xyyz_1, ta1_z_xzzzz_xyz_0, ta1_z_xzzzz_xyz_1, ta1_z_xzzzz_xyzz_0, ta1_z_xzzzz_xyzz_1, ta1_z_xzzzz_xzz_0, ta1_z_xzzzz_xzz_1, ta1_z_xzzzz_xzzz_0, ta1_z_xzzzz_xzzz_1, ta1_z_xzzzz_yyyy_0, ta1_z_xzzzz_yyyy_1, ta1_z_xzzzz_yyyz_0, ta1_z_xzzzz_yyyz_1, ta1_z_xzzzz_yyz_0, ta1_z_xzzzz_yyz_1, ta1_z_xzzzz_yyzz_0, ta1_z_xzzzz_yyzz_1, ta1_z_xzzzz_yzz_0, ta1_z_xzzzz_yzz_1, ta1_z_xzzzz_yzzz_0, ta1_z_xzzzz_yzzz_1, ta1_z_xzzzz_zzz_0, ta1_z_xzzzz_zzz_1, ta1_z_xzzzz_zzzz_0, ta1_z_xzzzz_zzzz_1, ta1_z_zzzz_xxxz_0, ta1_z_zzzz_xxxz_1, ta1_z_zzzz_xxyz_0, ta1_z_zzzz_xxyz_1, ta1_z_zzzz_xxzz_0, ta1_z_zzzz_xxzz_1, ta1_z_zzzz_xyyz_0, ta1_z_zzzz_xyyz_1, ta1_z_zzzz_xyzz_0, ta1_z_zzzz_xyzz_1, ta1_z_zzzz_xzzz_0, ta1_z_zzzz_xzzz_1, ta1_z_zzzz_yyyy_0, ta1_z_zzzz_yyyy_1, ta1_z_zzzz_yyyz_0, ta1_z_zzzz_yyyz_1, ta1_z_zzzz_yyzz_0, ta1_z_zzzz_yyzz_1, ta1_z_zzzz_yzzz_0, ta1_z_zzzz_yzzz_1, ta1_z_zzzz_zzzz_0, ta1_z_zzzz_zzzz_1, ta_xxzzz_xxxx_1, ta_xxzzz_xxxy_1, ta_xxzzz_xxyy_1, ta_xxzzz_xyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzzz_xxxx_0[i] = 3.0 * ta1_z_xxzz_xxxx_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxxx_1[i] * fe_0 + ta_xxzzz_xxxx_1[i] + ta1_z_xxzzz_xxxx_0[i] * pa_z[i] - ta1_z_xxzzz_xxxx_1[i] * pc_z[i];

        ta1_z_xxzzzz_xxxy_0[i] = 3.0 * ta1_z_xxzz_xxxy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxxy_1[i] * fe_0 + ta_xxzzz_xxxy_1[i] + ta1_z_xxzzz_xxxy_0[i] * pa_z[i] - ta1_z_xxzzz_xxxy_1[i] * pc_z[i];

        ta1_z_xxzzzz_xxxz_0[i] = ta1_z_zzzz_xxxz_0[i] * fe_0 - ta1_z_zzzz_xxxz_1[i] * fe_0 + 3.0 * ta1_z_xzzzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_xzzzz_xxz_1[i] * fe_0 + ta1_z_xzzzz_xxxz_0[i] * pa_x[i] - ta1_z_xzzzz_xxxz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xxyy_0[i] = 3.0 * ta1_z_xxzz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxyy_1[i] * fe_0 + ta_xxzzz_xxyy_1[i] + ta1_z_xxzzz_xxyy_0[i] * pa_z[i] - ta1_z_xxzzz_xxyy_1[i] * pc_z[i];

        ta1_z_xxzzzz_xxyz_0[i] = ta1_z_zzzz_xxyz_0[i] * fe_0 - ta1_z_zzzz_xxyz_1[i] * fe_0 + 2.0 * ta1_z_xzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xzzzz_xyz_1[i] * fe_0 + ta1_z_xzzzz_xxyz_0[i] * pa_x[i] - ta1_z_xzzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xxzz_0[i] = ta1_z_zzzz_xxzz_0[i] * fe_0 - ta1_z_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_xzzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_xzzzz_xzz_1[i] * fe_0 + ta1_z_xzzzz_xxzz_0[i] * pa_x[i] - ta1_z_xzzzz_xxzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xyyy_0[i] = 3.0 * ta1_z_xxzz_xyyy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xyyy_1[i] * fe_0 + ta_xxzzz_xyyy_1[i] + ta1_z_xxzzz_xyyy_0[i] * pa_z[i] - ta1_z_xxzzz_xyyy_1[i] * pc_z[i];

        ta1_z_xxzzzz_xyyz_0[i] = ta1_z_zzzz_xyyz_0[i] * fe_0 - ta1_z_zzzz_xyyz_1[i] * fe_0 + ta1_z_xzzzz_yyz_0[i] * fe_0 - ta1_z_xzzzz_yyz_1[i] * fe_0 + ta1_z_xzzzz_xyyz_0[i] * pa_x[i] - ta1_z_xzzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xyzz_0[i] = ta1_z_zzzz_xyzz_0[i] * fe_0 - ta1_z_zzzz_xyzz_1[i] * fe_0 + ta1_z_xzzzz_yzz_0[i] * fe_0 - ta1_z_xzzzz_yzz_1[i] * fe_0 + ta1_z_xzzzz_xyzz_0[i] * pa_x[i] - ta1_z_xzzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xzzz_0[i] = ta1_z_zzzz_xzzz_0[i] * fe_0 - ta1_z_zzzz_xzzz_1[i] * fe_0 + ta1_z_xzzzz_zzz_0[i] * fe_0 - ta1_z_xzzzz_zzz_1[i] * fe_0 + ta1_z_xzzzz_xzzz_0[i] * pa_x[i] - ta1_z_xzzzz_xzzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_yyyy_0[i] = ta1_z_zzzz_yyyy_0[i] * fe_0 - ta1_z_zzzz_yyyy_1[i] * fe_0 + ta1_z_xzzzz_yyyy_0[i] * pa_x[i] - ta1_z_xzzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xxzzzz_yyyz_0[i] = ta1_z_zzzz_yyyz_0[i] * fe_0 - ta1_z_zzzz_yyyz_1[i] * fe_0 + ta1_z_xzzzz_yyyz_0[i] * pa_x[i] - ta1_z_xzzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xxzzzz_yyzz_0[i] = ta1_z_zzzz_yyzz_0[i] * fe_0 - ta1_z_zzzz_yyzz_1[i] * fe_0 + ta1_z_xzzzz_yyzz_0[i] * pa_x[i] - ta1_z_xzzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_yzzz_0[i] = ta1_z_zzzz_yzzz_0[i] * fe_0 - ta1_z_zzzz_yzzz_1[i] * fe_0 + ta1_z_xzzzz_yzzz_0[i] * pa_x[i] - ta1_z_xzzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_zzzz_0[i] = ta1_z_zzzz_zzzz_0[i] * fe_0 - ta1_z_zzzz_zzzz_1[i] * fe_0 + ta1_z_xzzzz_zzzz_0[i] * pa_x[i] - ta1_z_xzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1065-1080 components of targeted buffer : IG

    auto ta1_z_xyyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1065);

    auto ta1_z_xyyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1066);

    auto ta1_z_xyyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1067);

    auto ta1_z_xyyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1068);

    auto ta1_z_xyyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1069);

    auto ta1_z_xyyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1070);

    auto ta1_z_xyyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1071);

    auto ta1_z_xyyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1072);

    auto ta1_z_xyyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1073);

    auto ta1_z_xyyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1074);

    auto ta1_z_xyyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1075);

    auto ta1_z_xyyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1076);

    auto ta1_z_xyyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1077);

    auto ta1_z_xyyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1078);

    auto ta1_z_xyyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1079);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyyyy_xxxx_0, ta1_z_xyyyyy_xxxy_0, ta1_z_xyyyyy_xxxz_0, ta1_z_xyyyyy_xxyy_0, ta1_z_xyyyyy_xxyz_0, ta1_z_xyyyyy_xxzz_0, ta1_z_xyyyyy_xyyy_0, ta1_z_xyyyyy_xyyz_0, ta1_z_xyyyyy_xyzz_0, ta1_z_xyyyyy_xzzz_0, ta1_z_xyyyyy_yyyy_0, ta1_z_xyyyyy_yyyz_0, ta1_z_xyyyyy_yyzz_0, ta1_z_xyyyyy_yzzz_0, ta1_z_xyyyyy_zzzz_0, ta1_z_yyyyy_xxx_0, ta1_z_yyyyy_xxx_1, ta1_z_yyyyy_xxxx_0, ta1_z_yyyyy_xxxx_1, ta1_z_yyyyy_xxxy_0, ta1_z_yyyyy_xxxy_1, ta1_z_yyyyy_xxxz_0, ta1_z_yyyyy_xxxz_1, ta1_z_yyyyy_xxy_0, ta1_z_yyyyy_xxy_1, ta1_z_yyyyy_xxyy_0, ta1_z_yyyyy_xxyy_1, ta1_z_yyyyy_xxyz_0, ta1_z_yyyyy_xxyz_1, ta1_z_yyyyy_xxz_0, ta1_z_yyyyy_xxz_1, ta1_z_yyyyy_xxzz_0, ta1_z_yyyyy_xxzz_1, ta1_z_yyyyy_xyy_0, ta1_z_yyyyy_xyy_1, ta1_z_yyyyy_xyyy_0, ta1_z_yyyyy_xyyy_1, ta1_z_yyyyy_xyyz_0, ta1_z_yyyyy_xyyz_1, ta1_z_yyyyy_xyz_0, ta1_z_yyyyy_xyz_1, ta1_z_yyyyy_xyzz_0, ta1_z_yyyyy_xyzz_1, ta1_z_yyyyy_xzz_0, ta1_z_yyyyy_xzz_1, ta1_z_yyyyy_xzzz_0, ta1_z_yyyyy_xzzz_1, ta1_z_yyyyy_yyy_0, ta1_z_yyyyy_yyy_1, ta1_z_yyyyy_yyyy_0, ta1_z_yyyyy_yyyy_1, ta1_z_yyyyy_yyyz_0, ta1_z_yyyyy_yyyz_1, ta1_z_yyyyy_yyz_0, ta1_z_yyyyy_yyz_1, ta1_z_yyyyy_yyzz_0, ta1_z_yyyyy_yyzz_1, ta1_z_yyyyy_yzz_0, ta1_z_yyyyy_yzz_1, ta1_z_yyyyy_yzzz_0, ta1_z_yyyyy_yzzz_1, ta1_z_yyyyy_zzz_0, ta1_z_yyyyy_zzz_1, ta1_z_yyyyy_zzzz_0, ta1_z_yyyyy_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyy_xxxx_0[i] = 4.0 * ta1_z_yyyyy_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyyyy_xxx_1[i] * fe_0 + ta1_z_yyyyy_xxxx_0[i] * pa_x[i] - ta1_z_yyyyy_xxxx_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxxy_0[i] = 3.0 * ta1_z_yyyyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyyyy_xxy_1[i] * fe_0 + ta1_z_yyyyy_xxxy_0[i] * pa_x[i] - ta1_z_yyyyy_xxxy_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxxz_0[i] = 3.0 * ta1_z_yyyyy_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyyyy_xxz_1[i] * fe_0 + ta1_z_yyyyy_xxxz_0[i] * pa_x[i] - ta1_z_yyyyy_xxxz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxyy_0[i] = 2.0 * ta1_z_yyyyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xyy_1[i] * fe_0 + ta1_z_yyyyy_xxyy_0[i] * pa_x[i] - ta1_z_yyyyy_xxyy_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxyz_0[i] = 2.0 * ta1_z_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xyz_1[i] * fe_0 + ta1_z_yyyyy_xxyz_0[i] * pa_x[i] - ta1_z_yyyyy_xxyz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxzz_0[i] = 2.0 * ta1_z_yyyyy_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xzz_1[i] * fe_0 + ta1_z_yyyyy_xxzz_0[i] * pa_x[i] - ta1_z_yyyyy_xxzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xyyy_0[i] = ta1_z_yyyyy_yyy_0[i] * fe_0 - ta1_z_yyyyy_yyy_1[i] * fe_0 + ta1_z_yyyyy_xyyy_0[i] * pa_x[i] - ta1_z_yyyyy_xyyy_1[i] * pc_x[i];

        ta1_z_xyyyyy_xyyz_0[i] = ta1_z_yyyyy_yyz_0[i] * fe_0 - ta1_z_yyyyy_yyz_1[i] * fe_0 + ta1_z_yyyyy_xyyz_0[i] * pa_x[i] - ta1_z_yyyyy_xyyz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xyzz_0[i] = ta1_z_yyyyy_yzz_0[i] * fe_0 - ta1_z_yyyyy_yzz_1[i] * fe_0 + ta1_z_yyyyy_xyzz_0[i] * pa_x[i] - ta1_z_yyyyy_xyzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xzzz_0[i] = ta1_z_yyyyy_zzz_0[i] * fe_0 - ta1_z_yyyyy_zzz_1[i] * fe_0 + ta1_z_yyyyy_xzzz_0[i] * pa_x[i] - ta1_z_yyyyy_xzzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_yyyy_0[i] = ta1_z_yyyyy_yyyy_0[i] * pa_x[i] - ta1_z_yyyyy_yyyy_1[i] * pc_x[i];

        ta1_z_xyyyyy_yyyz_0[i] = ta1_z_yyyyy_yyyz_0[i] * pa_x[i] - ta1_z_yyyyy_yyyz_1[i] * pc_x[i];

        ta1_z_xyyyyy_yyzz_0[i] = ta1_z_yyyyy_yyzz_0[i] * pa_x[i] - ta1_z_yyyyy_yyzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_yzzz_0[i] = ta1_z_yyyyy_yzzz_0[i] * pa_x[i] - ta1_z_yyyyy_yzzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_zzzz_0[i] = ta1_z_yyyyy_zzzz_0[i] * pa_x[i] - ta1_z_yyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 1080-1095 components of targeted buffer : IG

    auto ta1_z_xyyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1080);

    auto ta1_z_xyyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1081);

    auto ta1_z_xyyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1082);

    auto ta1_z_xyyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1083);

    auto ta1_z_xyyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1084);

    auto ta1_z_xyyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1085);

    auto ta1_z_xyyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1086);

    auto ta1_z_xyyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1087);

    auto ta1_z_xyyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1088);

    auto ta1_z_xyyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1089);

    auto ta1_z_xyyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1090);

    auto ta1_z_xyyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1091);

    auto ta1_z_xyyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1092);

    auto ta1_z_xyyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1093);

    auto ta1_z_xyyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1094);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xyyyy_xxxx_0, ta1_z_xyyyy_xxxx_1, ta1_z_xyyyy_xxxy_0, ta1_z_xyyyy_xxxy_1, ta1_z_xyyyy_xxyy_0, ta1_z_xyyyy_xxyy_1, ta1_z_xyyyy_xyyy_0, ta1_z_xyyyy_xyyy_1, ta1_z_xyyyyz_xxxx_0, ta1_z_xyyyyz_xxxy_0, ta1_z_xyyyyz_xxxz_0, ta1_z_xyyyyz_xxyy_0, ta1_z_xyyyyz_xxyz_0, ta1_z_xyyyyz_xxzz_0, ta1_z_xyyyyz_xyyy_0, ta1_z_xyyyyz_xyyz_0, ta1_z_xyyyyz_xyzz_0, ta1_z_xyyyyz_xzzz_0, ta1_z_xyyyyz_yyyy_0, ta1_z_xyyyyz_yyyz_0, ta1_z_xyyyyz_yyzz_0, ta1_z_xyyyyz_yzzz_0, ta1_z_xyyyyz_zzzz_0, ta1_z_yyyyz_xxxz_0, ta1_z_yyyyz_xxxz_1, ta1_z_yyyyz_xxyz_0, ta1_z_yyyyz_xxyz_1, ta1_z_yyyyz_xxz_0, ta1_z_yyyyz_xxz_1, ta1_z_yyyyz_xxzz_0, ta1_z_yyyyz_xxzz_1, ta1_z_yyyyz_xyyz_0, ta1_z_yyyyz_xyyz_1, ta1_z_yyyyz_xyz_0, ta1_z_yyyyz_xyz_1, ta1_z_yyyyz_xyzz_0, ta1_z_yyyyz_xyzz_1, ta1_z_yyyyz_xzz_0, ta1_z_yyyyz_xzz_1, ta1_z_yyyyz_xzzz_0, ta1_z_yyyyz_xzzz_1, ta1_z_yyyyz_yyyy_0, ta1_z_yyyyz_yyyy_1, ta1_z_yyyyz_yyyz_0, ta1_z_yyyyz_yyyz_1, ta1_z_yyyyz_yyz_0, ta1_z_yyyyz_yyz_1, ta1_z_yyyyz_yyzz_0, ta1_z_yyyyz_yyzz_1, ta1_z_yyyyz_yzz_0, ta1_z_yyyyz_yzz_1, ta1_z_yyyyz_yzzz_0, ta1_z_yyyyz_yzzz_1, ta1_z_yyyyz_zzz_0, ta1_z_yyyyz_zzz_1, ta1_z_yyyyz_zzzz_0, ta1_z_yyyyz_zzzz_1, ta_xyyyy_xxxx_1, ta_xyyyy_xxxy_1, ta_xyyyy_xxyy_1, ta_xyyyy_xyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyz_xxxx_0[i] = ta_xyyyy_xxxx_1[i] + ta1_z_xyyyy_xxxx_0[i] * pa_z[i] - ta1_z_xyyyy_xxxx_1[i] * pc_z[i];

        ta1_z_xyyyyz_xxxy_0[i] = ta_xyyyy_xxxy_1[i] + ta1_z_xyyyy_xxxy_0[i] * pa_z[i] - ta1_z_xyyyy_xxxy_1[i] * pc_z[i];

        ta1_z_xyyyyz_xxxz_0[i] = 3.0 * ta1_z_yyyyz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyyyz_xxz_1[i] * fe_0 + ta1_z_yyyyz_xxxz_0[i] * pa_x[i] - ta1_z_yyyyz_xxxz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xxyy_0[i] = ta_xyyyy_xxyy_1[i] + ta1_z_xyyyy_xxyy_0[i] * pa_z[i] - ta1_z_xyyyy_xxyy_1[i] * pc_z[i];

        ta1_z_xyyyyz_xxyz_0[i] = 2.0 * ta1_z_yyyyz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyyz_xyz_1[i] * fe_0 + ta1_z_yyyyz_xxyz_0[i] * pa_x[i] - ta1_z_yyyyz_xxyz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xxzz_0[i] = 2.0 * ta1_z_yyyyz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyyyz_xzz_1[i] * fe_0 + ta1_z_yyyyz_xxzz_0[i] * pa_x[i] - ta1_z_yyyyz_xxzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xyyy_0[i] = ta_xyyyy_xyyy_1[i] + ta1_z_xyyyy_xyyy_0[i] * pa_z[i] - ta1_z_xyyyy_xyyy_1[i] * pc_z[i];

        ta1_z_xyyyyz_xyyz_0[i] = ta1_z_yyyyz_yyz_0[i] * fe_0 - ta1_z_yyyyz_yyz_1[i] * fe_0 + ta1_z_yyyyz_xyyz_0[i] * pa_x[i] - ta1_z_yyyyz_xyyz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xyzz_0[i] = ta1_z_yyyyz_yzz_0[i] * fe_0 - ta1_z_yyyyz_yzz_1[i] * fe_0 + ta1_z_yyyyz_xyzz_0[i] * pa_x[i] - ta1_z_yyyyz_xyzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xzzz_0[i] = ta1_z_yyyyz_zzz_0[i] * fe_0 - ta1_z_yyyyz_zzz_1[i] * fe_0 + ta1_z_yyyyz_xzzz_0[i] * pa_x[i] - ta1_z_yyyyz_xzzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_yyyy_0[i] = ta1_z_yyyyz_yyyy_0[i] * pa_x[i] - ta1_z_yyyyz_yyyy_1[i] * pc_x[i];

        ta1_z_xyyyyz_yyyz_0[i] = ta1_z_yyyyz_yyyz_0[i] * pa_x[i] - ta1_z_yyyyz_yyyz_1[i] * pc_x[i];

        ta1_z_xyyyyz_yyzz_0[i] = ta1_z_yyyyz_yyzz_0[i] * pa_x[i] - ta1_z_yyyyz_yyzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_yzzz_0[i] = ta1_z_yyyyz_yzzz_0[i] * pa_x[i] - ta1_z_yyyyz_yzzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_zzzz_0[i] = ta1_z_yyyyz_zzzz_0[i] * pa_x[i] - ta1_z_yyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1095-1110 components of targeted buffer : IG

    auto ta1_z_xyyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1095);

    auto ta1_z_xyyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1096);

    auto ta1_z_xyyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1097);

    auto ta1_z_xyyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1098);

    auto ta1_z_xyyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1099);

    auto ta1_z_xyyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1100);

    auto ta1_z_xyyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1101);

    auto ta1_z_xyyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1102);

    auto ta1_z_xyyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1103);

    auto ta1_z_xyyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1104);

    auto ta1_z_xyyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1105);

    auto ta1_z_xyyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1106);

    auto ta1_z_xyyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1107);

    auto ta1_z_xyyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1108);

    auto ta1_z_xyyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1109);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyyzz_xxxx_0, ta1_z_xyyyzz_xxxy_0, ta1_z_xyyyzz_xxxz_0, ta1_z_xyyyzz_xxyy_0, ta1_z_xyyyzz_xxyz_0, ta1_z_xyyyzz_xxzz_0, ta1_z_xyyyzz_xyyy_0, ta1_z_xyyyzz_xyyz_0, ta1_z_xyyyzz_xyzz_0, ta1_z_xyyyzz_xzzz_0, ta1_z_xyyyzz_yyyy_0, ta1_z_xyyyzz_yyyz_0, ta1_z_xyyyzz_yyzz_0, ta1_z_xyyyzz_yzzz_0, ta1_z_xyyyzz_zzzz_0, ta1_z_yyyzz_xxx_0, ta1_z_yyyzz_xxx_1, ta1_z_yyyzz_xxxx_0, ta1_z_yyyzz_xxxx_1, ta1_z_yyyzz_xxxy_0, ta1_z_yyyzz_xxxy_1, ta1_z_yyyzz_xxxz_0, ta1_z_yyyzz_xxxz_1, ta1_z_yyyzz_xxy_0, ta1_z_yyyzz_xxy_1, ta1_z_yyyzz_xxyy_0, ta1_z_yyyzz_xxyy_1, ta1_z_yyyzz_xxyz_0, ta1_z_yyyzz_xxyz_1, ta1_z_yyyzz_xxz_0, ta1_z_yyyzz_xxz_1, ta1_z_yyyzz_xxzz_0, ta1_z_yyyzz_xxzz_1, ta1_z_yyyzz_xyy_0, ta1_z_yyyzz_xyy_1, ta1_z_yyyzz_xyyy_0, ta1_z_yyyzz_xyyy_1, ta1_z_yyyzz_xyyz_0, ta1_z_yyyzz_xyyz_1, ta1_z_yyyzz_xyz_0, ta1_z_yyyzz_xyz_1, ta1_z_yyyzz_xyzz_0, ta1_z_yyyzz_xyzz_1, ta1_z_yyyzz_xzz_0, ta1_z_yyyzz_xzz_1, ta1_z_yyyzz_xzzz_0, ta1_z_yyyzz_xzzz_1, ta1_z_yyyzz_yyy_0, ta1_z_yyyzz_yyy_1, ta1_z_yyyzz_yyyy_0, ta1_z_yyyzz_yyyy_1, ta1_z_yyyzz_yyyz_0, ta1_z_yyyzz_yyyz_1, ta1_z_yyyzz_yyz_0, ta1_z_yyyzz_yyz_1, ta1_z_yyyzz_yyzz_0, ta1_z_yyyzz_yyzz_1, ta1_z_yyyzz_yzz_0, ta1_z_yyyzz_yzz_1, ta1_z_yyyzz_yzzz_0, ta1_z_yyyzz_yzzz_1, ta1_z_yyyzz_zzz_0, ta1_z_yyyzz_zzz_1, ta1_z_yyyzz_zzzz_0, ta1_z_yyyzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyzz_xxxx_0[i] = 4.0 * ta1_z_yyyzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyyzz_xxx_1[i] * fe_0 + ta1_z_yyyzz_xxxx_0[i] * pa_x[i] - ta1_z_yyyzz_xxxx_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxxy_0[i] = 3.0 * ta1_z_yyyzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyyzz_xxy_1[i] * fe_0 + ta1_z_yyyzz_xxxy_0[i] * pa_x[i] - ta1_z_yyyzz_xxxy_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxxz_0[i] = 3.0 * ta1_z_yyyzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyyzz_xxz_1[i] * fe_0 + ta1_z_yyyzz_xxxz_0[i] * pa_x[i] - ta1_z_yyyzz_xxxz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxyy_0[i] = 2.0 * ta1_z_yyyzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_xyy_1[i] * fe_0 + ta1_z_yyyzz_xxyy_0[i] * pa_x[i] - ta1_z_yyyzz_xxyy_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxyz_0[i] = 2.0 * ta1_z_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_xyz_1[i] * fe_0 + ta1_z_yyyzz_xxyz_0[i] * pa_x[i] - ta1_z_yyyzz_xxyz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxzz_0[i] = 2.0 * ta1_z_yyyzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_xzz_1[i] * fe_0 + ta1_z_yyyzz_xxzz_0[i] * pa_x[i] - ta1_z_yyyzz_xxzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xyyy_0[i] = ta1_z_yyyzz_yyy_0[i] * fe_0 - ta1_z_yyyzz_yyy_1[i] * fe_0 + ta1_z_yyyzz_xyyy_0[i] * pa_x[i] - ta1_z_yyyzz_xyyy_1[i] * pc_x[i];

        ta1_z_xyyyzz_xyyz_0[i] = ta1_z_yyyzz_yyz_0[i] * fe_0 - ta1_z_yyyzz_yyz_1[i] * fe_0 + ta1_z_yyyzz_xyyz_0[i] * pa_x[i] - ta1_z_yyyzz_xyyz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xyzz_0[i] = ta1_z_yyyzz_yzz_0[i] * fe_0 - ta1_z_yyyzz_yzz_1[i] * fe_0 + ta1_z_yyyzz_xyzz_0[i] * pa_x[i] - ta1_z_yyyzz_xyzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xzzz_0[i] = ta1_z_yyyzz_zzz_0[i] * fe_0 - ta1_z_yyyzz_zzz_1[i] * fe_0 + ta1_z_yyyzz_xzzz_0[i] * pa_x[i] - ta1_z_yyyzz_xzzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_yyyy_0[i] = ta1_z_yyyzz_yyyy_0[i] * pa_x[i] - ta1_z_yyyzz_yyyy_1[i] * pc_x[i];

        ta1_z_xyyyzz_yyyz_0[i] = ta1_z_yyyzz_yyyz_0[i] * pa_x[i] - ta1_z_yyyzz_yyyz_1[i] * pc_x[i];

        ta1_z_xyyyzz_yyzz_0[i] = ta1_z_yyyzz_yyzz_0[i] * pa_x[i] - ta1_z_yyyzz_yyzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_yzzz_0[i] = ta1_z_yyyzz_yzzz_0[i] * pa_x[i] - ta1_z_yyyzz_yzzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_zzzz_0[i] = ta1_z_yyyzz_zzzz_0[i] * pa_x[i] - ta1_z_yyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1110-1125 components of targeted buffer : IG

    auto ta1_z_xyyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1110);

    auto ta1_z_xyyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1111);

    auto ta1_z_xyyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1112);

    auto ta1_z_xyyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1113);

    auto ta1_z_xyyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1114);

    auto ta1_z_xyyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1115);

    auto ta1_z_xyyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1116);

    auto ta1_z_xyyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1117);

    auto ta1_z_xyyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1118);

    auto ta1_z_xyyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1119);

    auto ta1_z_xyyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1120);

    auto ta1_z_xyyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1121);

    auto ta1_z_xyyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1122);

    auto ta1_z_xyyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1123);

    auto ta1_z_xyyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1124);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyzzz_xxxx_0, ta1_z_xyyzzz_xxxy_0, ta1_z_xyyzzz_xxxz_0, ta1_z_xyyzzz_xxyy_0, ta1_z_xyyzzz_xxyz_0, ta1_z_xyyzzz_xxzz_0, ta1_z_xyyzzz_xyyy_0, ta1_z_xyyzzz_xyyz_0, ta1_z_xyyzzz_xyzz_0, ta1_z_xyyzzz_xzzz_0, ta1_z_xyyzzz_yyyy_0, ta1_z_xyyzzz_yyyz_0, ta1_z_xyyzzz_yyzz_0, ta1_z_xyyzzz_yzzz_0, ta1_z_xyyzzz_zzzz_0, ta1_z_yyzzz_xxx_0, ta1_z_yyzzz_xxx_1, ta1_z_yyzzz_xxxx_0, ta1_z_yyzzz_xxxx_1, ta1_z_yyzzz_xxxy_0, ta1_z_yyzzz_xxxy_1, ta1_z_yyzzz_xxxz_0, ta1_z_yyzzz_xxxz_1, ta1_z_yyzzz_xxy_0, ta1_z_yyzzz_xxy_1, ta1_z_yyzzz_xxyy_0, ta1_z_yyzzz_xxyy_1, ta1_z_yyzzz_xxyz_0, ta1_z_yyzzz_xxyz_1, ta1_z_yyzzz_xxz_0, ta1_z_yyzzz_xxz_1, ta1_z_yyzzz_xxzz_0, ta1_z_yyzzz_xxzz_1, ta1_z_yyzzz_xyy_0, ta1_z_yyzzz_xyy_1, ta1_z_yyzzz_xyyy_0, ta1_z_yyzzz_xyyy_1, ta1_z_yyzzz_xyyz_0, ta1_z_yyzzz_xyyz_1, ta1_z_yyzzz_xyz_0, ta1_z_yyzzz_xyz_1, ta1_z_yyzzz_xyzz_0, ta1_z_yyzzz_xyzz_1, ta1_z_yyzzz_xzz_0, ta1_z_yyzzz_xzz_1, ta1_z_yyzzz_xzzz_0, ta1_z_yyzzz_xzzz_1, ta1_z_yyzzz_yyy_0, ta1_z_yyzzz_yyy_1, ta1_z_yyzzz_yyyy_0, ta1_z_yyzzz_yyyy_1, ta1_z_yyzzz_yyyz_0, ta1_z_yyzzz_yyyz_1, ta1_z_yyzzz_yyz_0, ta1_z_yyzzz_yyz_1, ta1_z_yyzzz_yyzz_0, ta1_z_yyzzz_yyzz_1, ta1_z_yyzzz_yzz_0, ta1_z_yyzzz_yzz_1, ta1_z_yyzzz_yzzz_0, ta1_z_yyzzz_yzzz_1, ta1_z_yyzzz_zzz_0, ta1_z_yyzzz_zzz_1, ta1_z_yyzzz_zzzz_0, ta1_z_yyzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzzz_xxxx_0[i] = 4.0 * ta1_z_yyzzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyzzz_xxx_1[i] * fe_0 + ta1_z_yyzzz_xxxx_0[i] * pa_x[i] - ta1_z_yyzzz_xxxx_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxxy_0[i] = 3.0 * ta1_z_yyzzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyzzz_xxy_1[i] * fe_0 + ta1_z_yyzzz_xxxy_0[i] * pa_x[i] - ta1_z_yyzzz_xxxy_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxxz_0[i] = 3.0 * ta1_z_yyzzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyzzz_xxz_1[i] * fe_0 + ta1_z_yyzzz_xxxz_0[i] * pa_x[i] - ta1_z_yyzzz_xxxz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxyy_0[i] = 2.0 * ta1_z_yyzzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_xyy_1[i] * fe_0 + ta1_z_yyzzz_xxyy_0[i] * pa_x[i] - ta1_z_yyzzz_xxyy_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxyz_0[i] = 2.0 * ta1_z_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_xyz_1[i] * fe_0 + ta1_z_yyzzz_xxyz_0[i] * pa_x[i] - ta1_z_yyzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxzz_0[i] = 2.0 * ta1_z_yyzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_xzz_1[i] * fe_0 + ta1_z_yyzzz_xxzz_0[i] * pa_x[i] - ta1_z_yyzzz_xxzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xyyy_0[i] = ta1_z_yyzzz_yyy_0[i] * fe_0 - ta1_z_yyzzz_yyy_1[i] * fe_0 + ta1_z_yyzzz_xyyy_0[i] * pa_x[i] - ta1_z_yyzzz_xyyy_1[i] * pc_x[i];

        ta1_z_xyyzzz_xyyz_0[i] = ta1_z_yyzzz_yyz_0[i] * fe_0 - ta1_z_yyzzz_yyz_1[i] * fe_0 + ta1_z_yyzzz_xyyz_0[i] * pa_x[i] - ta1_z_yyzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xyzz_0[i] = ta1_z_yyzzz_yzz_0[i] * fe_0 - ta1_z_yyzzz_yzz_1[i] * fe_0 + ta1_z_yyzzz_xyzz_0[i] * pa_x[i] - ta1_z_yyzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xzzz_0[i] = ta1_z_yyzzz_zzz_0[i] * fe_0 - ta1_z_yyzzz_zzz_1[i] * fe_0 + ta1_z_yyzzz_xzzz_0[i] * pa_x[i] - ta1_z_yyzzz_xzzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_yyyy_0[i] = ta1_z_yyzzz_yyyy_0[i] * pa_x[i] - ta1_z_yyzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xyyzzz_yyyz_0[i] = ta1_z_yyzzz_yyyz_0[i] * pa_x[i] - ta1_z_yyzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xyyzzz_yyzz_0[i] = ta1_z_yyzzz_yyzz_0[i] * pa_x[i] - ta1_z_yyzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_yzzz_0[i] = ta1_z_yyzzz_yzzz_0[i] * pa_x[i] - ta1_z_yyzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_zzzz_0[i] = ta1_z_yyzzz_zzzz_0[i] * pa_x[i] - ta1_z_yyzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1125-1140 components of targeted buffer : IG

    auto ta1_z_xyzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1125);

    auto ta1_z_xyzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1126);

    auto ta1_z_xyzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1127);

    auto ta1_z_xyzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1128);

    auto ta1_z_xyzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1129);

    auto ta1_z_xyzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1130);

    auto ta1_z_xyzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1131);

    auto ta1_z_xyzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1132);

    auto ta1_z_xyzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1133);

    auto ta1_z_xyzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1134);

    auto ta1_z_xyzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1135);

    auto ta1_z_xyzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1136);

    auto ta1_z_xyzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1137);

    auto ta1_z_xyzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1138);

    auto ta1_z_xyzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1139);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xyzzzz_xxxx_0, ta1_z_xyzzzz_xxxy_0, ta1_z_xyzzzz_xxxz_0, ta1_z_xyzzzz_xxyy_0, ta1_z_xyzzzz_xxyz_0, ta1_z_xyzzzz_xxzz_0, ta1_z_xyzzzz_xyyy_0, ta1_z_xyzzzz_xyyz_0, ta1_z_xyzzzz_xyzz_0, ta1_z_xyzzzz_xzzz_0, ta1_z_xyzzzz_yyyy_0, ta1_z_xyzzzz_yyyz_0, ta1_z_xyzzzz_yyzz_0, ta1_z_xyzzzz_yzzz_0, ta1_z_xyzzzz_zzzz_0, ta1_z_xzzzz_xxxx_0, ta1_z_xzzzz_xxxx_1, ta1_z_xzzzz_xxxz_0, ta1_z_xzzzz_xxxz_1, ta1_z_xzzzz_xxzz_0, ta1_z_xzzzz_xxzz_1, ta1_z_xzzzz_xzzz_0, ta1_z_xzzzz_xzzz_1, ta1_z_yzzzz_xxxy_0, ta1_z_yzzzz_xxxy_1, ta1_z_yzzzz_xxy_0, ta1_z_yzzzz_xxy_1, ta1_z_yzzzz_xxyy_0, ta1_z_yzzzz_xxyy_1, ta1_z_yzzzz_xxyz_0, ta1_z_yzzzz_xxyz_1, ta1_z_yzzzz_xyy_0, ta1_z_yzzzz_xyy_1, ta1_z_yzzzz_xyyy_0, ta1_z_yzzzz_xyyy_1, ta1_z_yzzzz_xyyz_0, ta1_z_yzzzz_xyyz_1, ta1_z_yzzzz_xyz_0, ta1_z_yzzzz_xyz_1, ta1_z_yzzzz_xyzz_0, ta1_z_yzzzz_xyzz_1, ta1_z_yzzzz_yyy_0, ta1_z_yzzzz_yyy_1, ta1_z_yzzzz_yyyy_0, ta1_z_yzzzz_yyyy_1, ta1_z_yzzzz_yyyz_0, ta1_z_yzzzz_yyyz_1, ta1_z_yzzzz_yyz_0, ta1_z_yzzzz_yyz_1, ta1_z_yzzzz_yyzz_0, ta1_z_yzzzz_yyzz_1, ta1_z_yzzzz_yzz_0, ta1_z_yzzzz_yzz_1, ta1_z_yzzzz_yzzz_0, ta1_z_yzzzz_yzzz_1, ta1_z_yzzzz_zzzz_0, ta1_z_yzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzzzz_xxxx_0[i] = ta1_z_xzzzz_xxxx_0[i] * pa_y[i] - ta1_z_xzzzz_xxxx_1[i] * pc_y[i];

        ta1_z_xyzzzz_xxxy_0[i] = 3.0 * ta1_z_yzzzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yzzzz_xxy_1[i] * fe_0 + ta1_z_yzzzz_xxxy_0[i] * pa_x[i] - ta1_z_yzzzz_xxxy_1[i] * pc_x[i];

        ta1_z_xyzzzz_xxxz_0[i] = ta1_z_xzzzz_xxxz_0[i] * pa_y[i] - ta1_z_xzzzz_xxxz_1[i] * pc_y[i];

        ta1_z_xyzzzz_xxyy_0[i] = 2.0 * ta1_z_yzzzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yzzzz_xyy_1[i] * fe_0 + ta1_z_yzzzz_xxyy_0[i] * pa_x[i] - ta1_z_yzzzz_xxyy_1[i] * pc_x[i];

        ta1_z_xyzzzz_xxyz_0[i] = 2.0 * ta1_z_yzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yzzzz_xyz_1[i] * fe_0 + ta1_z_yzzzz_xxyz_0[i] * pa_x[i] - ta1_z_yzzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xyzzzz_xxzz_0[i] = ta1_z_xzzzz_xxzz_0[i] * pa_y[i] - ta1_z_xzzzz_xxzz_1[i] * pc_y[i];

        ta1_z_xyzzzz_xyyy_0[i] = ta1_z_yzzzz_yyy_0[i] * fe_0 - ta1_z_yzzzz_yyy_1[i] * fe_0 + ta1_z_yzzzz_xyyy_0[i] * pa_x[i] - ta1_z_yzzzz_xyyy_1[i] * pc_x[i];

        ta1_z_xyzzzz_xyyz_0[i] = ta1_z_yzzzz_yyz_0[i] * fe_0 - ta1_z_yzzzz_yyz_1[i] * fe_0 + ta1_z_yzzzz_xyyz_0[i] * pa_x[i] - ta1_z_yzzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xyzzzz_xyzz_0[i] = ta1_z_yzzzz_yzz_0[i] * fe_0 - ta1_z_yzzzz_yzz_1[i] * fe_0 + ta1_z_yzzzz_xyzz_0[i] * pa_x[i] - ta1_z_yzzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xyzzzz_xzzz_0[i] = ta1_z_xzzzz_xzzz_0[i] * pa_y[i] - ta1_z_xzzzz_xzzz_1[i] * pc_y[i];

        ta1_z_xyzzzz_yyyy_0[i] = ta1_z_yzzzz_yyyy_0[i] * pa_x[i] - ta1_z_yzzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xyzzzz_yyyz_0[i] = ta1_z_yzzzz_yyyz_0[i] * pa_x[i] - ta1_z_yzzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xyzzzz_yyzz_0[i] = ta1_z_yzzzz_yyzz_0[i] * pa_x[i] - ta1_z_yzzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xyzzzz_yzzz_0[i] = ta1_z_yzzzz_yzzz_0[i] * pa_x[i] - ta1_z_yzzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xyzzzz_zzzz_0[i] = ta1_z_yzzzz_zzzz_0[i] * pa_x[i] - ta1_z_yzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1140-1155 components of targeted buffer : IG

    auto ta1_z_xzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1140);

    auto ta1_z_xzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1141);

    auto ta1_z_xzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1142);

    auto ta1_z_xzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1143);

    auto ta1_z_xzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1144);

    auto ta1_z_xzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1145);

    auto ta1_z_xzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1146);

    auto ta1_z_xzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1147);

    auto ta1_z_xzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1148);

    auto ta1_z_xzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1149);

    auto ta1_z_xzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1150);

    auto ta1_z_xzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1151);

    auto ta1_z_xzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1152);

    auto ta1_z_xzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1153);

    auto ta1_z_xzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1154);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xzzzzz_xxxx_0, ta1_z_xzzzzz_xxxy_0, ta1_z_xzzzzz_xxxz_0, ta1_z_xzzzzz_xxyy_0, ta1_z_xzzzzz_xxyz_0, ta1_z_xzzzzz_xxzz_0, ta1_z_xzzzzz_xyyy_0, ta1_z_xzzzzz_xyyz_0, ta1_z_xzzzzz_xyzz_0, ta1_z_xzzzzz_xzzz_0, ta1_z_xzzzzz_yyyy_0, ta1_z_xzzzzz_yyyz_0, ta1_z_xzzzzz_yyzz_0, ta1_z_xzzzzz_yzzz_0, ta1_z_xzzzzz_zzzz_0, ta1_z_zzzzz_xxx_0, ta1_z_zzzzz_xxx_1, ta1_z_zzzzz_xxxx_0, ta1_z_zzzzz_xxxx_1, ta1_z_zzzzz_xxxy_0, ta1_z_zzzzz_xxxy_1, ta1_z_zzzzz_xxxz_0, ta1_z_zzzzz_xxxz_1, ta1_z_zzzzz_xxy_0, ta1_z_zzzzz_xxy_1, ta1_z_zzzzz_xxyy_0, ta1_z_zzzzz_xxyy_1, ta1_z_zzzzz_xxyz_0, ta1_z_zzzzz_xxyz_1, ta1_z_zzzzz_xxz_0, ta1_z_zzzzz_xxz_1, ta1_z_zzzzz_xxzz_0, ta1_z_zzzzz_xxzz_1, ta1_z_zzzzz_xyy_0, ta1_z_zzzzz_xyy_1, ta1_z_zzzzz_xyyy_0, ta1_z_zzzzz_xyyy_1, ta1_z_zzzzz_xyyz_0, ta1_z_zzzzz_xyyz_1, ta1_z_zzzzz_xyz_0, ta1_z_zzzzz_xyz_1, ta1_z_zzzzz_xyzz_0, ta1_z_zzzzz_xyzz_1, ta1_z_zzzzz_xzz_0, ta1_z_zzzzz_xzz_1, ta1_z_zzzzz_xzzz_0, ta1_z_zzzzz_xzzz_1, ta1_z_zzzzz_yyy_0, ta1_z_zzzzz_yyy_1, ta1_z_zzzzz_yyyy_0, ta1_z_zzzzz_yyyy_1, ta1_z_zzzzz_yyyz_0, ta1_z_zzzzz_yyyz_1, ta1_z_zzzzz_yyz_0, ta1_z_zzzzz_yyz_1, ta1_z_zzzzz_yyzz_0, ta1_z_zzzzz_yyzz_1, ta1_z_zzzzz_yzz_0, ta1_z_zzzzz_yzz_1, ta1_z_zzzzz_yzzz_0, ta1_z_zzzzz_yzzz_1, ta1_z_zzzzz_zzz_0, ta1_z_zzzzz_zzz_1, ta1_z_zzzzz_zzzz_0, ta1_z_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzzz_xxxx_0[i] = 4.0 * ta1_z_zzzzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_zzzzz_xxx_1[i] * fe_0 + ta1_z_zzzzz_xxxx_0[i] * pa_x[i] - ta1_z_zzzzz_xxxx_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxxy_0[i] = 3.0 * ta1_z_zzzzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_xxy_1[i] * fe_0 + ta1_z_zzzzz_xxxy_0[i] * pa_x[i] - ta1_z_zzzzz_xxxy_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxxz_0[i] = 3.0 * ta1_z_zzzzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_xxz_1[i] * fe_0 + ta1_z_zzzzz_xxxz_0[i] * pa_x[i] - ta1_z_zzzzz_xxxz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxyy_0[i] = 2.0 * ta1_z_zzzzz_xyy_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xyy_1[i] * fe_0 + ta1_z_zzzzz_xxyy_0[i] * pa_x[i] - ta1_z_zzzzz_xxyy_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxyz_0[i] = 2.0 * ta1_z_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xyz_1[i] * fe_0 + ta1_z_zzzzz_xxyz_0[i] * pa_x[i] - ta1_z_zzzzz_xxyz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxzz_0[i] = 2.0 * ta1_z_zzzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xzz_1[i] * fe_0 + ta1_z_zzzzz_xxzz_0[i] * pa_x[i] - ta1_z_zzzzz_xxzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xyyy_0[i] = ta1_z_zzzzz_yyy_0[i] * fe_0 - ta1_z_zzzzz_yyy_1[i] * fe_0 + ta1_z_zzzzz_xyyy_0[i] * pa_x[i] - ta1_z_zzzzz_xyyy_1[i] * pc_x[i];

        ta1_z_xzzzzz_xyyz_0[i] = ta1_z_zzzzz_yyz_0[i] * fe_0 - ta1_z_zzzzz_yyz_1[i] * fe_0 + ta1_z_zzzzz_xyyz_0[i] * pa_x[i] - ta1_z_zzzzz_xyyz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xyzz_0[i] = ta1_z_zzzzz_yzz_0[i] * fe_0 - ta1_z_zzzzz_yzz_1[i] * fe_0 + ta1_z_zzzzz_xyzz_0[i] * pa_x[i] - ta1_z_zzzzz_xyzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xzzz_0[i] = ta1_z_zzzzz_zzz_0[i] * fe_0 - ta1_z_zzzzz_zzz_1[i] * fe_0 + ta1_z_zzzzz_xzzz_0[i] * pa_x[i] - ta1_z_zzzzz_xzzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_yyyy_0[i] = ta1_z_zzzzz_yyyy_0[i] * pa_x[i] - ta1_z_zzzzz_yyyy_1[i] * pc_x[i];

        ta1_z_xzzzzz_yyyz_0[i] = ta1_z_zzzzz_yyyz_0[i] * pa_x[i] - ta1_z_zzzzz_yyyz_1[i] * pc_x[i];

        ta1_z_xzzzzz_yyzz_0[i] = ta1_z_zzzzz_yyzz_0[i] * pa_x[i] - ta1_z_zzzzz_yyzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_yzzz_0[i] = ta1_z_zzzzz_yzzz_0[i] * pa_x[i] - ta1_z_zzzzz_yzzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_zzzz_0[i] = ta1_z_zzzzz_zzzz_0[i] * pa_x[i] - ta1_z_zzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1155-1170 components of targeted buffer : IG

    auto ta1_z_yyyyyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1155);

    auto ta1_z_yyyyyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1156);

    auto ta1_z_yyyyyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1157);

    auto ta1_z_yyyyyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1158);

    auto ta1_z_yyyyyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1159);

    auto ta1_z_yyyyyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1160);

    auto ta1_z_yyyyyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1161);

    auto ta1_z_yyyyyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1162);

    auto ta1_z_yyyyyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1163);

    auto ta1_z_yyyyyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1164);

    auto ta1_z_yyyyyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1165);

    auto ta1_z_yyyyyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1166);

    auto ta1_z_yyyyyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1167);

    auto ta1_z_yyyyyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1168);

    auto ta1_z_yyyyyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1169);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yyyy_xxxx_0, ta1_z_yyyy_xxxx_1, ta1_z_yyyy_xxxy_0, ta1_z_yyyy_xxxy_1, ta1_z_yyyy_xxxz_0, ta1_z_yyyy_xxxz_1, ta1_z_yyyy_xxyy_0, ta1_z_yyyy_xxyy_1, ta1_z_yyyy_xxyz_0, ta1_z_yyyy_xxyz_1, ta1_z_yyyy_xxzz_0, ta1_z_yyyy_xxzz_1, ta1_z_yyyy_xyyy_0, ta1_z_yyyy_xyyy_1, ta1_z_yyyy_xyyz_0, ta1_z_yyyy_xyyz_1, ta1_z_yyyy_xyzz_0, ta1_z_yyyy_xyzz_1, ta1_z_yyyy_xzzz_0, ta1_z_yyyy_xzzz_1, ta1_z_yyyy_yyyy_0, ta1_z_yyyy_yyyy_1, ta1_z_yyyy_yyyz_0, ta1_z_yyyy_yyyz_1, ta1_z_yyyy_yyzz_0, ta1_z_yyyy_yyzz_1, ta1_z_yyyy_yzzz_0, ta1_z_yyyy_yzzz_1, ta1_z_yyyy_zzzz_0, ta1_z_yyyy_zzzz_1, ta1_z_yyyyy_xxx_0, ta1_z_yyyyy_xxx_1, ta1_z_yyyyy_xxxx_0, ta1_z_yyyyy_xxxx_1, ta1_z_yyyyy_xxxy_0, ta1_z_yyyyy_xxxy_1, ta1_z_yyyyy_xxxz_0, ta1_z_yyyyy_xxxz_1, ta1_z_yyyyy_xxy_0, ta1_z_yyyyy_xxy_1, ta1_z_yyyyy_xxyy_0, ta1_z_yyyyy_xxyy_1, ta1_z_yyyyy_xxyz_0, ta1_z_yyyyy_xxyz_1, ta1_z_yyyyy_xxz_0, ta1_z_yyyyy_xxz_1, ta1_z_yyyyy_xxzz_0, ta1_z_yyyyy_xxzz_1, ta1_z_yyyyy_xyy_0, ta1_z_yyyyy_xyy_1, ta1_z_yyyyy_xyyy_0, ta1_z_yyyyy_xyyy_1, ta1_z_yyyyy_xyyz_0, ta1_z_yyyyy_xyyz_1, ta1_z_yyyyy_xyz_0, ta1_z_yyyyy_xyz_1, ta1_z_yyyyy_xyzz_0, ta1_z_yyyyy_xyzz_1, ta1_z_yyyyy_xzz_0, ta1_z_yyyyy_xzz_1, ta1_z_yyyyy_xzzz_0, ta1_z_yyyyy_xzzz_1, ta1_z_yyyyy_yyy_0, ta1_z_yyyyy_yyy_1, ta1_z_yyyyy_yyyy_0, ta1_z_yyyyy_yyyy_1, ta1_z_yyyyy_yyyz_0, ta1_z_yyyyy_yyyz_1, ta1_z_yyyyy_yyz_0, ta1_z_yyyyy_yyz_1, ta1_z_yyyyy_yyzz_0, ta1_z_yyyyy_yyzz_1, ta1_z_yyyyy_yzz_0, ta1_z_yyyyy_yzz_1, ta1_z_yyyyy_yzzz_0, ta1_z_yyyyy_yzzz_1, ta1_z_yyyyy_zzz_0, ta1_z_yyyyy_zzz_1, ta1_z_yyyyy_zzzz_0, ta1_z_yyyyy_zzzz_1, ta1_z_yyyyyy_xxxx_0, ta1_z_yyyyyy_xxxy_0, ta1_z_yyyyyy_xxxz_0, ta1_z_yyyyyy_xxyy_0, ta1_z_yyyyyy_xxyz_0, ta1_z_yyyyyy_xxzz_0, ta1_z_yyyyyy_xyyy_0, ta1_z_yyyyyy_xyyz_0, ta1_z_yyyyyy_xyzz_0, ta1_z_yyyyyy_xzzz_0, ta1_z_yyyyyy_yyyy_0, ta1_z_yyyyyy_yyyz_0, ta1_z_yyyyyy_yyzz_0, ta1_z_yyyyyy_yzzz_0, ta1_z_yyyyyy_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyy_xxxx_0[i] = 5.0 * ta1_z_yyyy_xxxx_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxxx_1[i] * fe_0 + ta1_z_yyyyy_xxxx_0[i] * pa_y[i] - ta1_z_yyyyy_xxxx_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxxy_0[i] = 5.0 * ta1_z_yyyy_xxxy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxxy_1[i] * fe_0 + ta1_z_yyyyy_xxx_0[i] * fe_0 - ta1_z_yyyyy_xxx_1[i] * fe_0 + ta1_z_yyyyy_xxxy_0[i] * pa_y[i] - ta1_z_yyyyy_xxxy_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxxz_0[i] = 5.0 * ta1_z_yyyy_xxxz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxxz_1[i] * fe_0 + ta1_z_yyyyy_xxxz_0[i] * pa_y[i] - ta1_z_yyyyy_xxxz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxyy_0[i] = 5.0 * ta1_z_yyyy_xxyy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_yyyyy_xxy_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xxy_1[i] * fe_0 + ta1_z_yyyyy_xxyy_0[i] * pa_y[i] - ta1_z_yyyyy_xxyy_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxyz_0[i] = 5.0 * ta1_z_yyyy_xxyz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxyz_1[i] * fe_0 + ta1_z_yyyyy_xxz_0[i] * fe_0 - ta1_z_yyyyy_xxz_1[i] * fe_0 + ta1_z_yyyyy_xxyz_0[i] * pa_y[i] - ta1_z_yyyyy_xxyz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxzz_0[i] = 5.0 * ta1_z_yyyy_xxzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxzz_1[i] * fe_0 + ta1_z_yyyyy_xxzz_0[i] * pa_y[i] - ta1_z_yyyyy_xxzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xyyy_0[i] = 5.0 * ta1_z_yyyy_xyyy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xyyy_1[i] * fe_0 + 3.0 * ta1_z_yyyyy_xyy_0[i] * fe_0 - 3.0 * ta1_z_yyyyy_xyy_1[i] * fe_0 + ta1_z_yyyyy_xyyy_0[i] * pa_y[i] - ta1_z_yyyyy_xyyy_1[i] * pc_y[i];

        ta1_z_yyyyyy_xyyz_0[i] = 5.0 * ta1_z_yyyy_xyyz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xyz_1[i] * fe_0 + ta1_z_yyyyy_xyyz_0[i] * pa_y[i] - ta1_z_yyyyy_xyyz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xyzz_0[i] = 5.0 * ta1_z_yyyy_xyzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xyzz_1[i] * fe_0 + ta1_z_yyyyy_xzz_0[i] * fe_0 - ta1_z_yyyyy_xzz_1[i] * fe_0 + ta1_z_yyyyy_xyzz_0[i] * pa_y[i] - ta1_z_yyyyy_xyzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xzzz_0[i] = 5.0 * ta1_z_yyyy_xzzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xzzz_1[i] * fe_0 + ta1_z_yyyyy_xzzz_0[i] * pa_y[i] - ta1_z_yyyyy_xzzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_yyyy_0[i] = 5.0 * ta1_z_yyyy_yyyy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yyyy_1[i] * fe_0 + 4.0 * ta1_z_yyyyy_yyy_0[i] * fe_0 - 4.0 * ta1_z_yyyyy_yyy_1[i] * fe_0 + ta1_z_yyyyy_yyyy_0[i] * pa_y[i] - ta1_z_yyyyy_yyyy_1[i] * pc_y[i];

        ta1_z_yyyyyy_yyyz_0[i] = 5.0 * ta1_z_yyyy_yyyz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yyyyy_yyz_0[i] * fe_0 - 3.0 * ta1_z_yyyyy_yyz_1[i] * fe_0 + ta1_z_yyyyy_yyyz_0[i] * pa_y[i] - ta1_z_yyyyy_yyyz_1[i] * pc_y[i];

        ta1_z_yyyyyy_yyzz_0[i] = 5.0 * ta1_z_yyyy_yyzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yyyyy_yzz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_yzz_1[i] * fe_0 + ta1_z_yyyyy_yyzz_0[i] * pa_y[i] - ta1_z_yyyyy_yyzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_yzzz_0[i] = 5.0 * ta1_z_yyyy_yzzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yzzz_1[i] * fe_0 + ta1_z_yyyyy_zzz_0[i] * fe_0 - ta1_z_yyyyy_zzz_1[i] * fe_0 + ta1_z_yyyyy_yzzz_0[i] * pa_y[i] - ta1_z_yyyyy_yzzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_zzzz_0[i] = 5.0 * ta1_z_yyyy_zzzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_zzzz_1[i] * fe_0 + ta1_z_yyyyy_zzzz_0[i] * pa_y[i] - ta1_z_yyyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 1170-1185 components of targeted buffer : IG

    auto ta1_z_yyyyyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1170);

    auto ta1_z_yyyyyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1171);

    auto ta1_z_yyyyyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1172);

    auto ta1_z_yyyyyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1173);

    auto ta1_z_yyyyyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1174);

    auto ta1_z_yyyyyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1175);

    auto ta1_z_yyyyyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1176);

    auto ta1_z_yyyyyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1177);

    auto ta1_z_yyyyyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1178);

    auto ta1_z_yyyyyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1179);

    auto ta1_z_yyyyyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1180);

    auto ta1_z_yyyyyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1181);

    auto ta1_z_yyyyyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1182);

    auto ta1_z_yyyyyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1183);

    auto ta1_z_yyyyyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1184);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyyy_xxxx_0, ta1_z_yyyyy_xxxx_1, ta1_z_yyyyy_xxxy_0, ta1_z_yyyyy_xxxy_1, ta1_z_yyyyy_xxy_0, ta1_z_yyyyy_xxy_1, ta1_z_yyyyy_xxyy_0, ta1_z_yyyyy_xxyy_1, ta1_z_yyyyy_xxyz_0, ta1_z_yyyyy_xxyz_1, ta1_z_yyyyy_xyy_0, ta1_z_yyyyy_xyy_1, ta1_z_yyyyy_xyyy_0, ta1_z_yyyyy_xyyy_1, ta1_z_yyyyy_xyyz_0, ta1_z_yyyyy_xyyz_1, ta1_z_yyyyy_xyz_0, ta1_z_yyyyy_xyz_1, ta1_z_yyyyy_xyzz_0, ta1_z_yyyyy_xyzz_1, ta1_z_yyyyy_yyy_0, ta1_z_yyyyy_yyy_1, ta1_z_yyyyy_yyyy_0, ta1_z_yyyyy_yyyy_1, ta1_z_yyyyy_yyyz_0, ta1_z_yyyyy_yyyz_1, ta1_z_yyyyy_yyz_0, ta1_z_yyyyy_yyz_1, ta1_z_yyyyy_yyzz_0, ta1_z_yyyyy_yyzz_1, ta1_z_yyyyy_yzz_0, ta1_z_yyyyy_yzz_1, ta1_z_yyyyy_yzzz_0, ta1_z_yyyyy_yzzz_1, ta1_z_yyyyyz_xxxx_0, ta1_z_yyyyyz_xxxy_0, ta1_z_yyyyyz_xxxz_0, ta1_z_yyyyyz_xxyy_0, ta1_z_yyyyyz_xxyz_0, ta1_z_yyyyyz_xxzz_0, ta1_z_yyyyyz_xyyy_0, ta1_z_yyyyyz_xyyz_0, ta1_z_yyyyyz_xyzz_0, ta1_z_yyyyyz_xzzz_0, ta1_z_yyyyyz_yyyy_0, ta1_z_yyyyyz_yyyz_0, ta1_z_yyyyyz_yyzz_0, ta1_z_yyyyyz_yzzz_0, ta1_z_yyyyyz_zzzz_0, ta1_z_yyyyz_xxxz_0, ta1_z_yyyyz_xxxz_1, ta1_z_yyyyz_xxzz_0, ta1_z_yyyyz_xxzz_1, ta1_z_yyyyz_xzzz_0, ta1_z_yyyyz_xzzz_1, ta1_z_yyyyz_zzzz_0, ta1_z_yyyyz_zzzz_1, ta1_z_yyyz_xxxz_0, ta1_z_yyyz_xxxz_1, ta1_z_yyyz_xxzz_0, ta1_z_yyyz_xxzz_1, ta1_z_yyyz_xzzz_0, ta1_z_yyyz_xzzz_1, ta1_z_yyyz_zzzz_0, ta1_z_yyyz_zzzz_1, ta_yyyyy_xxxx_1, ta_yyyyy_xxxy_1, ta_yyyyy_xxyy_1, ta_yyyyy_xxyz_1, ta_yyyyy_xyyy_1, ta_yyyyy_xyyz_1, ta_yyyyy_xyzz_1, ta_yyyyy_yyyy_1, ta_yyyyy_yyyz_1, ta_yyyyy_yyzz_1, ta_yyyyy_yzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyz_xxxx_0[i] = ta_yyyyy_xxxx_1[i] + ta1_z_yyyyy_xxxx_0[i] * pa_z[i] - ta1_z_yyyyy_xxxx_1[i] * pc_z[i];

        ta1_z_yyyyyz_xxxy_0[i] = ta_yyyyy_xxxy_1[i] + ta1_z_yyyyy_xxxy_0[i] * pa_z[i] - ta1_z_yyyyy_xxxy_1[i] * pc_z[i];

        ta1_z_yyyyyz_xxxz_0[i] = 4.0 * ta1_z_yyyz_xxxz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_xxxz_1[i] * fe_0 + ta1_z_yyyyz_xxxz_0[i] * pa_y[i] - ta1_z_yyyyz_xxxz_1[i] * pc_y[i];

        ta1_z_yyyyyz_xxyy_0[i] = ta_yyyyy_xxyy_1[i] + ta1_z_yyyyy_xxyy_0[i] * pa_z[i] - ta1_z_yyyyy_xxyy_1[i] * pc_z[i];

        ta1_z_yyyyyz_xxyz_0[i] = ta1_z_yyyyy_xxy_0[i] * fe_0 - ta1_z_yyyyy_xxy_1[i] * fe_0 + ta_yyyyy_xxyz_1[i] + ta1_z_yyyyy_xxyz_0[i] * pa_z[i] - ta1_z_yyyyy_xxyz_1[i] * pc_z[i];

        ta1_z_yyyyyz_xxzz_0[i] = 4.0 * ta1_z_yyyz_xxzz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_xxzz_1[i] * fe_0 + ta1_z_yyyyz_xxzz_0[i] * pa_y[i] - ta1_z_yyyyz_xxzz_1[i] * pc_y[i];

        ta1_z_yyyyyz_xyyy_0[i] = ta_yyyyy_xyyy_1[i] + ta1_z_yyyyy_xyyy_0[i] * pa_z[i] - ta1_z_yyyyy_xyyy_1[i] * pc_z[i];

        ta1_z_yyyyyz_xyyz_0[i] = ta1_z_yyyyy_xyy_0[i] * fe_0 - ta1_z_yyyyy_xyy_1[i] * fe_0 + ta_yyyyy_xyyz_1[i] + ta1_z_yyyyy_xyyz_0[i] * pa_z[i] - ta1_z_yyyyy_xyyz_1[i] * pc_z[i];

        ta1_z_yyyyyz_xyzz_0[i] = 2.0 * ta1_z_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xyz_1[i] * fe_0 + ta_yyyyy_xyzz_1[i] + ta1_z_yyyyy_xyzz_0[i] * pa_z[i] - ta1_z_yyyyy_xyzz_1[i] * pc_z[i];

        ta1_z_yyyyyz_xzzz_0[i] = 4.0 * ta1_z_yyyz_xzzz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_xzzz_1[i] * fe_0 + ta1_z_yyyyz_xzzz_0[i] * pa_y[i] - ta1_z_yyyyz_xzzz_1[i] * pc_y[i];

        ta1_z_yyyyyz_yyyy_0[i] = ta_yyyyy_yyyy_1[i] + ta1_z_yyyyy_yyyy_0[i] * pa_z[i] - ta1_z_yyyyy_yyyy_1[i] * pc_z[i];

        ta1_z_yyyyyz_yyyz_0[i] = ta1_z_yyyyy_yyy_0[i] * fe_0 - ta1_z_yyyyy_yyy_1[i] * fe_0 + ta_yyyyy_yyyz_1[i] + ta1_z_yyyyy_yyyz_0[i] * pa_z[i] - ta1_z_yyyyy_yyyz_1[i] * pc_z[i];

        ta1_z_yyyyyz_yyzz_0[i] = 2.0 * ta1_z_yyyyy_yyz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_yyz_1[i] * fe_0 + ta_yyyyy_yyzz_1[i] + ta1_z_yyyyy_yyzz_0[i] * pa_z[i] - ta1_z_yyyyy_yyzz_1[i] * pc_z[i];

        ta1_z_yyyyyz_yzzz_0[i] = 3.0 * ta1_z_yyyyy_yzz_0[i] * fe_0 - 3.0 * ta1_z_yyyyy_yzz_1[i] * fe_0 + ta_yyyyy_yzzz_1[i] + ta1_z_yyyyy_yzzz_0[i] * pa_z[i] - ta1_z_yyyyy_yzzz_1[i] * pc_z[i];

        ta1_z_yyyyyz_zzzz_0[i] = 4.0 * ta1_z_yyyz_zzzz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_zzzz_1[i] * fe_0 + ta1_z_yyyyz_zzzz_0[i] * pa_y[i] - ta1_z_yyyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1185-1200 components of targeted buffer : IG

    auto ta1_z_yyyyzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1185);

    auto ta1_z_yyyyzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1186);

    auto ta1_z_yyyyzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1187);

    auto ta1_z_yyyyzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1188);

    auto ta1_z_yyyyzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1189);

    auto ta1_z_yyyyzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1190);

    auto ta1_z_yyyyzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1191);

    auto ta1_z_yyyyzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1192);

    auto ta1_z_yyyyzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1193);

    auto ta1_z_yyyyzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1194);

    auto ta1_z_yyyyzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1195);

    auto ta1_z_yyyyzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1196);

    auto ta1_z_yyyyzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1197);

    auto ta1_z_yyyyzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1198);

    auto ta1_z_yyyyzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1199);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyy_xxxy_0, ta1_z_yyyy_xxxy_1, ta1_z_yyyy_xxyy_0, ta1_z_yyyy_xxyy_1, ta1_z_yyyy_xyyy_0, ta1_z_yyyy_xyyy_1, ta1_z_yyyy_yyyy_0, ta1_z_yyyy_yyyy_1, ta1_z_yyyyz_xxxy_0, ta1_z_yyyyz_xxxy_1, ta1_z_yyyyz_xxyy_0, ta1_z_yyyyz_xxyy_1, ta1_z_yyyyz_xyyy_0, ta1_z_yyyyz_xyyy_1, ta1_z_yyyyz_yyyy_0, ta1_z_yyyyz_yyyy_1, ta1_z_yyyyzz_xxxx_0, ta1_z_yyyyzz_xxxy_0, ta1_z_yyyyzz_xxxz_0, ta1_z_yyyyzz_xxyy_0, ta1_z_yyyyzz_xxyz_0, ta1_z_yyyyzz_xxzz_0, ta1_z_yyyyzz_xyyy_0, ta1_z_yyyyzz_xyyz_0, ta1_z_yyyyzz_xyzz_0, ta1_z_yyyyzz_xzzz_0, ta1_z_yyyyzz_yyyy_0, ta1_z_yyyyzz_yyyz_0, ta1_z_yyyyzz_yyzz_0, ta1_z_yyyyzz_yzzz_0, ta1_z_yyyyzz_zzzz_0, ta1_z_yyyzz_xxxx_0, ta1_z_yyyzz_xxxx_1, ta1_z_yyyzz_xxxz_0, ta1_z_yyyzz_xxxz_1, ta1_z_yyyzz_xxyz_0, ta1_z_yyyzz_xxyz_1, ta1_z_yyyzz_xxz_0, ta1_z_yyyzz_xxz_1, ta1_z_yyyzz_xxzz_0, ta1_z_yyyzz_xxzz_1, ta1_z_yyyzz_xyyz_0, ta1_z_yyyzz_xyyz_1, ta1_z_yyyzz_xyz_0, ta1_z_yyyzz_xyz_1, ta1_z_yyyzz_xyzz_0, ta1_z_yyyzz_xyzz_1, ta1_z_yyyzz_xzz_0, ta1_z_yyyzz_xzz_1, ta1_z_yyyzz_xzzz_0, ta1_z_yyyzz_xzzz_1, ta1_z_yyyzz_yyyz_0, ta1_z_yyyzz_yyyz_1, ta1_z_yyyzz_yyz_0, ta1_z_yyyzz_yyz_1, ta1_z_yyyzz_yyzz_0, ta1_z_yyyzz_yyzz_1, ta1_z_yyyzz_yzz_0, ta1_z_yyyzz_yzz_1, ta1_z_yyyzz_yzzz_0, ta1_z_yyyzz_yzzz_1, ta1_z_yyyzz_zzz_0, ta1_z_yyyzz_zzz_1, ta1_z_yyyzz_zzzz_0, ta1_z_yyyzz_zzzz_1, ta1_z_yyzz_xxxx_0, ta1_z_yyzz_xxxx_1, ta1_z_yyzz_xxxz_0, ta1_z_yyzz_xxxz_1, ta1_z_yyzz_xxyz_0, ta1_z_yyzz_xxyz_1, ta1_z_yyzz_xxzz_0, ta1_z_yyzz_xxzz_1, ta1_z_yyzz_xyyz_0, ta1_z_yyzz_xyyz_1, ta1_z_yyzz_xyzz_0, ta1_z_yyzz_xyzz_1, ta1_z_yyzz_xzzz_0, ta1_z_yyzz_xzzz_1, ta1_z_yyzz_yyyz_0, ta1_z_yyzz_yyyz_1, ta1_z_yyzz_yyzz_0, ta1_z_yyzz_yyzz_1, ta1_z_yyzz_yzzz_0, ta1_z_yyzz_yzzz_1, ta1_z_yyzz_zzzz_0, ta1_z_yyzz_zzzz_1, ta_yyyyz_xxxy_1, ta_yyyyz_xxyy_1, ta_yyyyz_xyyy_1, ta_yyyyz_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyzz_xxxx_0[i] = 3.0 * ta1_z_yyzz_xxxx_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxxx_1[i] * fe_0 + ta1_z_yyyzz_xxxx_0[i] * pa_y[i] - ta1_z_yyyzz_xxxx_1[i] * pc_y[i];

        ta1_z_yyyyzz_xxxy_0[i] = ta1_z_yyyy_xxxy_0[i] * fe_0 - ta1_z_yyyy_xxxy_1[i] * fe_0 + ta_yyyyz_xxxy_1[i] + ta1_z_yyyyz_xxxy_0[i] * pa_z[i] - ta1_z_yyyyz_xxxy_1[i] * pc_z[i];

        ta1_z_yyyyzz_xxxz_0[i] = 3.0 * ta1_z_yyzz_xxxz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxxz_1[i] * fe_0 + ta1_z_yyyzz_xxxz_0[i] * pa_y[i] - ta1_z_yyyzz_xxxz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xxyy_0[i] = ta1_z_yyyy_xxyy_0[i] * fe_0 - ta1_z_yyyy_xxyy_1[i] * fe_0 + ta_yyyyz_xxyy_1[i] + ta1_z_yyyyz_xxyy_0[i] * pa_z[i] - ta1_z_yyyyz_xxyy_1[i] * pc_z[i];

        ta1_z_yyyyzz_xxyz_0[i] = 3.0 * ta1_z_yyzz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxyz_1[i] * fe_0 + ta1_z_yyyzz_xxz_0[i] * fe_0 - ta1_z_yyyzz_xxz_1[i] * fe_0 + ta1_z_yyyzz_xxyz_0[i] * pa_y[i] - ta1_z_yyyzz_xxyz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xxzz_0[i] = 3.0 * ta1_z_yyzz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxzz_1[i] * fe_0 + ta1_z_yyyzz_xxzz_0[i] * pa_y[i] - ta1_z_yyyzz_xxzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xyyy_0[i] = ta1_z_yyyy_xyyy_0[i] * fe_0 - ta1_z_yyyy_xyyy_1[i] * fe_0 + ta_yyyyz_xyyy_1[i] + ta1_z_yyyyz_xyyy_0[i] * pa_z[i] - ta1_z_yyyyz_xyyy_1[i] * pc_z[i];

        ta1_z_yyyyzz_xyyz_0[i] = 3.0 * ta1_z_yyzz_xyyz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_xyz_1[i] * fe_0 + ta1_z_yyyzz_xyyz_0[i] * pa_y[i] - ta1_z_yyyzz_xyyz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xyzz_0[i] = 3.0 * ta1_z_yyzz_xyzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xyzz_1[i] * fe_0 + ta1_z_yyyzz_xzz_0[i] * fe_0 - ta1_z_yyyzz_xzz_1[i] * fe_0 + ta1_z_yyyzz_xyzz_0[i] * pa_y[i] - ta1_z_yyyzz_xyzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xzzz_0[i] = 3.0 * ta1_z_yyzz_xzzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xzzz_1[i] * fe_0 + ta1_z_yyyzz_xzzz_0[i] * pa_y[i] - ta1_z_yyyzz_xzzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_yyyy_0[i] = ta1_z_yyyy_yyyy_0[i] * fe_0 - ta1_z_yyyy_yyyy_1[i] * fe_0 + ta_yyyyz_yyyy_1[i] + ta1_z_yyyyz_yyyy_0[i] * pa_z[i] - ta1_z_yyyyz_yyyy_1[i] * pc_z[i];

        ta1_z_yyyyzz_yyyz_0[i] = 3.0 * ta1_z_yyzz_yyyz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yyyzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_yyyzz_yyz_1[i] * fe_0 + ta1_z_yyyzz_yyyz_0[i] * pa_y[i] - ta1_z_yyyzz_yyyz_1[i] * pc_y[i];

        ta1_z_yyyyzz_yyzz_0[i] = 3.0 * ta1_z_yyzz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yyyzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_yzz_1[i] * fe_0 + ta1_z_yyyzz_yyzz_0[i] * pa_y[i] - ta1_z_yyyzz_yyzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_yzzz_0[i] = 3.0 * ta1_z_yyzz_yzzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yzzz_1[i] * fe_0 + ta1_z_yyyzz_zzz_0[i] * fe_0 - ta1_z_yyyzz_zzz_1[i] * fe_0 + ta1_z_yyyzz_yzzz_0[i] * pa_y[i] - ta1_z_yyyzz_yzzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_zzzz_0[i] = 3.0 * ta1_z_yyzz_zzzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_zzzz_1[i] * fe_0 + ta1_z_yyyzz_zzzz_0[i] * pa_y[i] - ta1_z_yyyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1200-1215 components of targeted buffer : IG

    auto ta1_z_yyyzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1200);

    auto ta1_z_yyyzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1201);

    auto ta1_z_yyyzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1202);

    auto ta1_z_yyyzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1203);

    auto ta1_z_yyyzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1204);

    auto ta1_z_yyyzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1205);

    auto ta1_z_yyyzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1206);

    auto ta1_z_yyyzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1207);

    auto ta1_z_yyyzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1208);

    auto ta1_z_yyyzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1209);

    auto ta1_z_yyyzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1210);

    auto ta1_z_yyyzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1211);

    auto ta1_z_yyyzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1212);

    auto ta1_z_yyyzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1213);

    auto ta1_z_yyyzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1214);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyz_xxxy_0, ta1_z_yyyz_xxxy_1, ta1_z_yyyz_xxyy_0, ta1_z_yyyz_xxyy_1, ta1_z_yyyz_xyyy_0, ta1_z_yyyz_xyyy_1, ta1_z_yyyz_yyyy_0, ta1_z_yyyz_yyyy_1, ta1_z_yyyzz_xxxy_0, ta1_z_yyyzz_xxxy_1, ta1_z_yyyzz_xxyy_0, ta1_z_yyyzz_xxyy_1, ta1_z_yyyzz_xyyy_0, ta1_z_yyyzz_xyyy_1, ta1_z_yyyzz_yyyy_0, ta1_z_yyyzz_yyyy_1, ta1_z_yyyzzz_xxxx_0, ta1_z_yyyzzz_xxxy_0, ta1_z_yyyzzz_xxxz_0, ta1_z_yyyzzz_xxyy_0, ta1_z_yyyzzz_xxyz_0, ta1_z_yyyzzz_xxzz_0, ta1_z_yyyzzz_xyyy_0, ta1_z_yyyzzz_xyyz_0, ta1_z_yyyzzz_xyzz_0, ta1_z_yyyzzz_xzzz_0, ta1_z_yyyzzz_yyyy_0, ta1_z_yyyzzz_yyyz_0, ta1_z_yyyzzz_yyzz_0, ta1_z_yyyzzz_yzzz_0, ta1_z_yyyzzz_zzzz_0, ta1_z_yyzzz_xxxx_0, ta1_z_yyzzz_xxxx_1, ta1_z_yyzzz_xxxz_0, ta1_z_yyzzz_xxxz_1, ta1_z_yyzzz_xxyz_0, ta1_z_yyzzz_xxyz_1, ta1_z_yyzzz_xxz_0, ta1_z_yyzzz_xxz_1, ta1_z_yyzzz_xxzz_0, ta1_z_yyzzz_xxzz_1, ta1_z_yyzzz_xyyz_0, ta1_z_yyzzz_xyyz_1, ta1_z_yyzzz_xyz_0, ta1_z_yyzzz_xyz_1, ta1_z_yyzzz_xyzz_0, ta1_z_yyzzz_xyzz_1, ta1_z_yyzzz_xzz_0, ta1_z_yyzzz_xzz_1, ta1_z_yyzzz_xzzz_0, ta1_z_yyzzz_xzzz_1, ta1_z_yyzzz_yyyz_0, ta1_z_yyzzz_yyyz_1, ta1_z_yyzzz_yyz_0, ta1_z_yyzzz_yyz_1, ta1_z_yyzzz_yyzz_0, ta1_z_yyzzz_yyzz_1, ta1_z_yyzzz_yzz_0, ta1_z_yyzzz_yzz_1, ta1_z_yyzzz_yzzz_0, ta1_z_yyzzz_yzzz_1, ta1_z_yyzzz_zzz_0, ta1_z_yyzzz_zzz_1, ta1_z_yyzzz_zzzz_0, ta1_z_yyzzz_zzzz_1, ta1_z_yzzz_xxxx_0, ta1_z_yzzz_xxxx_1, ta1_z_yzzz_xxxz_0, ta1_z_yzzz_xxxz_1, ta1_z_yzzz_xxyz_0, ta1_z_yzzz_xxyz_1, ta1_z_yzzz_xxzz_0, ta1_z_yzzz_xxzz_1, ta1_z_yzzz_xyyz_0, ta1_z_yzzz_xyyz_1, ta1_z_yzzz_xyzz_0, ta1_z_yzzz_xyzz_1, ta1_z_yzzz_xzzz_0, ta1_z_yzzz_xzzz_1, ta1_z_yzzz_yyyz_0, ta1_z_yzzz_yyyz_1, ta1_z_yzzz_yyzz_0, ta1_z_yzzz_yyzz_1, ta1_z_yzzz_yzzz_0, ta1_z_yzzz_yzzz_1, ta1_z_yzzz_zzzz_0, ta1_z_yzzz_zzzz_1, ta_yyyzz_xxxy_1, ta_yyyzz_xxyy_1, ta_yyyzz_xyyy_1, ta_yyyzz_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzzz_xxxx_0[i] = 2.0 * ta1_z_yzzz_xxxx_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xxxx_1[i] * fe_0 + ta1_z_yyzzz_xxxx_0[i] * pa_y[i] - ta1_z_yyzzz_xxxx_1[i] * pc_y[i];

        ta1_z_yyyzzz_xxxy_0[i] = 2.0 * ta1_z_yyyz_xxxy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xxxy_1[i] * fe_0 + ta_yyyzz_xxxy_1[i] + ta1_z_yyyzz_xxxy_0[i] * pa_z[i] - ta1_z_yyyzz_xxxy_1[i] * pc_z[i];

        ta1_z_yyyzzz_xxxz_0[i] = 2.0 * ta1_z_yzzz_xxxz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xxxz_1[i] * fe_0 + ta1_z_yyzzz_xxxz_0[i] * pa_y[i] - ta1_z_yyzzz_xxxz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xxyy_0[i] = 2.0 * ta1_z_yyyz_xxyy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xxyy_1[i] * fe_0 + ta_yyyzz_xxyy_1[i] + ta1_z_yyyzz_xxyy_0[i] * pa_z[i] - ta1_z_yyyzz_xxyy_1[i] * pc_z[i];

        ta1_z_yyyzzz_xxyz_0[i] = 2.0 * ta1_z_yzzz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xxyz_1[i] * fe_0 + ta1_z_yyzzz_xxz_0[i] * fe_0 - ta1_z_yyzzz_xxz_1[i] * fe_0 + ta1_z_yyzzz_xxyz_0[i] * pa_y[i] - ta1_z_yyzzz_xxyz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xxzz_0[i] = 2.0 * ta1_z_yzzz_xxzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xxzz_1[i] * fe_0 + ta1_z_yyzzz_xxzz_0[i] * pa_y[i] - ta1_z_yyzzz_xxzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xyyy_0[i] = 2.0 * ta1_z_yyyz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xyyy_1[i] * fe_0 + ta_yyyzz_xyyy_1[i] + ta1_z_yyyzz_xyyy_0[i] * pa_z[i] - ta1_z_yyyzz_xyyy_1[i] * pc_z[i];

        ta1_z_yyyzzz_xyyz_0[i] = 2.0 * ta1_z_yzzz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_xyz_1[i] * fe_0 + ta1_z_yyzzz_xyyz_0[i] * pa_y[i] - ta1_z_yyzzz_xyyz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xyzz_0[i] = 2.0 * ta1_z_yzzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xyzz_1[i] * fe_0 + ta1_z_yyzzz_xzz_0[i] * fe_0 - ta1_z_yyzzz_xzz_1[i] * fe_0 + ta1_z_yyzzz_xyzz_0[i] * pa_y[i] - ta1_z_yyzzz_xyzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xzzz_0[i] = 2.0 * ta1_z_yzzz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xzzz_1[i] * fe_0 + ta1_z_yyzzz_xzzz_0[i] * pa_y[i] - ta1_z_yyzzz_xzzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_yyyy_0[i] = 2.0 * ta1_z_yyyz_yyyy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_yyyy_1[i] * fe_0 + ta_yyyzz_yyyy_1[i] + ta1_z_yyyzz_yyyy_0[i] * pa_z[i] - ta1_z_yyyzz_yyyy_1[i] * pc_z[i];

        ta1_z_yyyzzz_yyyz_0[i] = 2.0 * ta1_z_yzzz_yyyz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yyzzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_yyzzz_yyz_1[i] * fe_0 + ta1_z_yyzzz_yyyz_0[i] * pa_y[i] - ta1_z_yyzzz_yyyz_1[i] * pc_y[i];

        ta1_z_yyyzzz_yyzz_0[i] = 2.0 * ta1_z_yzzz_yyzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yyzzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_yzz_1[i] * fe_0 + ta1_z_yyzzz_yyzz_0[i] * pa_y[i] - ta1_z_yyzzz_yyzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_yzzz_0[i] = 2.0 * ta1_z_yzzz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yzzz_1[i] * fe_0 + ta1_z_yyzzz_zzz_0[i] * fe_0 - ta1_z_yyzzz_zzz_1[i] * fe_0 + ta1_z_yyzzz_yzzz_0[i] * pa_y[i] - ta1_z_yyzzz_yzzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_zzzz_0[i] = 2.0 * ta1_z_yzzz_zzzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_zzzz_1[i] * fe_0 + ta1_z_yyzzz_zzzz_0[i] * pa_y[i] - ta1_z_yyzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1215-1230 components of targeted buffer : IG

    auto ta1_z_yyzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1215);

    auto ta1_z_yyzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1216);

    auto ta1_z_yyzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1217);

    auto ta1_z_yyzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1218);

    auto ta1_z_yyzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1219);

    auto ta1_z_yyzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1220);

    auto ta1_z_yyzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1221);

    auto ta1_z_yyzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1222);

    auto ta1_z_yyzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1223);

    auto ta1_z_yyzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1224);

    auto ta1_z_yyzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1225);

    auto ta1_z_yyzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1226);

    auto ta1_z_yyzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1227);

    auto ta1_z_yyzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1228);

    auto ta1_z_yyzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1229);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyzz_xxxy_0, ta1_z_yyzz_xxxy_1, ta1_z_yyzz_xxyy_0, ta1_z_yyzz_xxyy_1, ta1_z_yyzz_xyyy_0, ta1_z_yyzz_xyyy_1, ta1_z_yyzz_yyyy_0, ta1_z_yyzz_yyyy_1, ta1_z_yyzzz_xxxy_0, ta1_z_yyzzz_xxxy_1, ta1_z_yyzzz_xxyy_0, ta1_z_yyzzz_xxyy_1, ta1_z_yyzzz_xyyy_0, ta1_z_yyzzz_xyyy_1, ta1_z_yyzzz_yyyy_0, ta1_z_yyzzz_yyyy_1, ta1_z_yyzzzz_xxxx_0, ta1_z_yyzzzz_xxxy_0, ta1_z_yyzzzz_xxxz_0, ta1_z_yyzzzz_xxyy_0, ta1_z_yyzzzz_xxyz_0, ta1_z_yyzzzz_xxzz_0, ta1_z_yyzzzz_xyyy_0, ta1_z_yyzzzz_xyyz_0, ta1_z_yyzzzz_xyzz_0, ta1_z_yyzzzz_xzzz_0, ta1_z_yyzzzz_yyyy_0, ta1_z_yyzzzz_yyyz_0, ta1_z_yyzzzz_yyzz_0, ta1_z_yyzzzz_yzzz_0, ta1_z_yyzzzz_zzzz_0, ta1_z_yzzzz_xxxx_0, ta1_z_yzzzz_xxxx_1, ta1_z_yzzzz_xxxz_0, ta1_z_yzzzz_xxxz_1, ta1_z_yzzzz_xxyz_0, ta1_z_yzzzz_xxyz_1, ta1_z_yzzzz_xxz_0, ta1_z_yzzzz_xxz_1, ta1_z_yzzzz_xxzz_0, ta1_z_yzzzz_xxzz_1, ta1_z_yzzzz_xyyz_0, ta1_z_yzzzz_xyyz_1, ta1_z_yzzzz_xyz_0, ta1_z_yzzzz_xyz_1, ta1_z_yzzzz_xyzz_0, ta1_z_yzzzz_xyzz_1, ta1_z_yzzzz_xzz_0, ta1_z_yzzzz_xzz_1, ta1_z_yzzzz_xzzz_0, ta1_z_yzzzz_xzzz_1, ta1_z_yzzzz_yyyz_0, ta1_z_yzzzz_yyyz_1, ta1_z_yzzzz_yyz_0, ta1_z_yzzzz_yyz_1, ta1_z_yzzzz_yyzz_0, ta1_z_yzzzz_yyzz_1, ta1_z_yzzzz_yzz_0, ta1_z_yzzzz_yzz_1, ta1_z_yzzzz_yzzz_0, ta1_z_yzzzz_yzzz_1, ta1_z_yzzzz_zzz_0, ta1_z_yzzzz_zzz_1, ta1_z_yzzzz_zzzz_0, ta1_z_yzzzz_zzzz_1, ta1_z_zzzz_xxxx_0, ta1_z_zzzz_xxxx_1, ta1_z_zzzz_xxxz_0, ta1_z_zzzz_xxxz_1, ta1_z_zzzz_xxyz_0, ta1_z_zzzz_xxyz_1, ta1_z_zzzz_xxzz_0, ta1_z_zzzz_xxzz_1, ta1_z_zzzz_xyyz_0, ta1_z_zzzz_xyyz_1, ta1_z_zzzz_xyzz_0, ta1_z_zzzz_xyzz_1, ta1_z_zzzz_xzzz_0, ta1_z_zzzz_xzzz_1, ta1_z_zzzz_yyyz_0, ta1_z_zzzz_yyyz_1, ta1_z_zzzz_yyzz_0, ta1_z_zzzz_yyzz_1, ta1_z_zzzz_yzzz_0, ta1_z_zzzz_yzzz_1, ta1_z_zzzz_zzzz_0, ta1_z_zzzz_zzzz_1, ta_yyzzz_xxxy_1, ta_yyzzz_xxyy_1, ta_yyzzz_xyyy_1, ta_yyzzz_yyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzzz_xxxx_0[i] = ta1_z_zzzz_xxxx_0[i] * fe_0 - ta1_z_zzzz_xxxx_1[i] * fe_0 + ta1_z_yzzzz_xxxx_0[i] * pa_y[i] - ta1_z_yzzzz_xxxx_1[i] * pc_y[i];

        ta1_z_yyzzzz_xxxy_0[i] = 3.0 * ta1_z_yyzz_xxxy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxxy_1[i] * fe_0 + ta_yyzzz_xxxy_1[i] + ta1_z_yyzzz_xxxy_0[i] * pa_z[i] - ta1_z_yyzzz_xxxy_1[i] * pc_z[i];

        ta1_z_yyzzzz_xxxz_0[i] = ta1_z_zzzz_xxxz_0[i] * fe_0 - ta1_z_zzzz_xxxz_1[i] * fe_0 + ta1_z_yzzzz_xxxz_0[i] * pa_y[i] - ta1_z_yzzzz_xxxz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xxyy_0[i] = 3.0 * ta1_z_yyzz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxyy_1[i] * fe_0 + ta_yyzzz_xxyy_1[i] + ta1_z_yyzzz_xxyy_0[i] * pa_z[i] - ta1_z_yyzzz_xxyy_1[i] * pc_z[i];

        ta1_z_yyzzzz_xxyz_0[i] = ta1_z_zzzz_xxyz_0[i] * fe_0 - ta1_z_zzzz_xxyz_1[i] * fe_0 + ta1_z_yzzzz_xxz_0[i] * fe_0 - ta1_z_yzzzz_xxz_1[i] * fe_0 + ta1_z_yzzzz_xxyz_0[i] * pa_y[i] - ta1_z_yzzzz_xxyz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xxzz_0[i] = ta1_z_zzzz_xxzz_0[i] * fe_0 - ta1_z_zzzz_xxzz_1[i] * fe_0 + ta1_z_yzzzz_xxzz_0[i] * pa_y[i] - ta1_z_yzzzz_xxzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xyyy_0[i] = 3.0 * ta1_z_yyzz_xyyy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xyyy_1[i] * fe_0 + ta_yyzzz_xyyy_1[i] + ta1_z_yyzzz_xyyy_0[i] * pa_z[i] - ta1_z_yyzzz_xyyy_1[i] * pc_z[i];

        ta1_z_yyzzzz_xyyz_0[i] = ta1_z_zzzz_xyyz_0[i] * fe_0 - ta1_z_zzzz_xyyz_1[i] * fe_0 + 2.0 * ta1_z_yzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yzzzz_xyz_1[i] * fe_0 + ta1_z_yzzzz_xyyz_0[i] * pa_y[i] - ta1_z_yzzzz_xyyz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xyzz_0[i] = ta1_z_zzzz_xyzz_0[i] * fe_0 - ta1_z_zzzz_xyzz_1[i] * fe_0 + ta1_z_yzzzz_xzz_0[i] * fe_0 - ta1_z_yzzzz_xzz_1[i] * fe_0 + ta1_z_yzzzz_xyzz_0[i] * pa_y[i] - ta1_z_yzzzz_xyzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xzzz_0[i] = ta1_z_zzzz_xzzz_0[i] * fe_0 - ta1_z_zzzz_xzzz_1[i] * fe_0 + ta1_z_yzzzz_xzzz_0[i] * pa_y[i] - ta1_z_yzzzz_xzzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_yyyy_0[i] = 3.0 * ta1_z_yyzz_yyyy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yyyy_1[i] * fe_0 + ta_yyzzz_yyyy_1[i] + ta1_z_yyzzz_yyyy_0[i] * pa_z[i] - ta1_z_yyzzz_yyyy_1[i] * pc_z[i];

        ta1_z_yyzzzz_yyyz_0[i] = ta1_z_zzzz_yyyz_0[i] * fe_0 - ta1_z_zzzz_yyyz_1[i] * fe_0 + 3.0 * ta1_z_yzzzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_yzzzz_yyz_1[i] * fe_0 + ta1_z_yzzzz_yyyz_0[i] * pa_y[i] - ta1_z_yzzzz_yyyz_1[i] * pc_y[i];

        ta1_z_yyzzzz_yyzz_0[i] = ta1_z_zzzz_yyzz_0[i] * fe_0 - ta1_z_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_yzzzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_yzzzz_yzz_1[i] * fe_0 + ta1_z_yzzzz_yyzz_0[i] * pa_y[i] - ta1_z_yzzzz_yyzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_yzzz_0[i] = ta1_z_zzzz_yzzz_0[i] * fe_0 - ta1_z_zzzz_yzzz_1[i] * fe_0 + ta1_z_yzzzz_zzz_0[i] * fe_0 - ta1_z_yzzzz_zzz_1[i] * fe_0 + ta1_z_yzzzz_yzzz_0[i] * pa_y[i] - ta1_z_yzzzz_yzzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_zzzz_0[i] = ta1_z_zzzz_zzzz_0[i] * fe_0 - ta1_z_zzzz_zzzz_1[i] * fe_0 + ta1_z_yzzzz_zzzz_0[i] * pa_y[i] - ta1_z_yzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1230-1245 components of targeted buffer : IG

    auto ta1_z_yzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1230);

    auto ta1_z_yzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1231);

    auto ta1_z_yzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1232);

    auto ta1_z_yzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1233);

    auto ta1_z_yzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1234);

    auto ta1_z_yzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1235);

    auto ta1_z_yzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1236);

    auto ta1_z_yzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1237);

    auto ta1_z_yzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1238);

    auto ta1_z_yzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1239);

    auto ta1_z_yzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1240);

    auto ta1_z_yzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1241);

    auto ta1_z_yzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1242);

    auto ta1_z_yzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1243);

    auto ta1_z_yzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1244);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yzzzzz_xxxx_0, ta1_z_yzzzzz_xxxy_0, ta1_z_yzzzzz_xxxz_0, ta1_z_yzzzzz_xxyy_0, ta1_z_yzzzzz_xxyz_0, ta1_z_yzzzzz_xxzz_0, ta1_z_yzzzzz_xyyy_0, ta1_z_yzzzzz_xyyz_0, ta1_z_yzzzzz_xyzz_0, ta1_z_yzzzzz_xzzz_0, ta1_z_yzzzzz_yyyy_0, ta1_z_yzzzzz_yyyz_0, ta1_z_yzzzzz_yyzz_0, ta1_z_yzzzzz_yzzz_0, ta1_z_yzzzzz_zzzz_0, ta1_z_zzzzz_xxx_0, ta1_z_zzzzz_xxx_1, ta1_z_zzzzz_xxxx_0, ta1_z_zzzzz_xxxx_1, ta1_z_zzzzz_xxxy_0, ta1_z_zzzzz_xxxy_1, ta1_z_zzzzz_xxxz_0, ta1_z_zzzzz_xxxz_1, ta1_z_zzzzz_xxy_0, ta1_z_zzzzz_xxy_1, ta1_z_zzzzz_xxyy_0, ta1_z_zzzzz_xxyy_1, ta1_z_zzzzz_xxyz_0, ta1_z_zzzzz_xxyz_1, ta1_z_zzzzz_xxz_0, ta1_z_zzzzz_xxz_1, ta1_z_zzzzz_xxzz_0, ta1_z_zzzzz_xxzz_1, ta1_z_zzzzz_xyy_0, ta1_z_zzzzz_xyy_1, ta1_z_zzzzz_xyyy_0, ta1_z_zzzzz_xyyy_1, ta1_z_zzzzz_xyyz_0, ta1_z_zzzzz_xyyz_1, ta1_z_zzzzz_xyz_0, ta1_z_zzzzz_xyz_1, ta1_z_zzzzz_xyzz_0, ta1_z_zzzzz_xyzz_1, ta1_z_zzzzz_xzz_0, ta1_z_zzzzz_xzz_1, ta1_z_zzzzz_xzzz_0, ta1_z_zzzzz_xzzz_1, ta1_z_zzzzz_yyy_0, ta1_z_zzzzz_yyy_1, ta1_z_zzzzz_yyyy_0, ta1_z_zzzzz_yyyy_1, ta1_z_zzzzz_yyyz_0, ta1_z_zzzzz_yyyz_1, ta1_z_zzzzz_yyz_0, ta1_z_zzzzz_yyz_1, ta1_z_zzzzz_yyzz_0, ta1_z_zzzzz_yyzz_1, ta1_z_zzzzz_yzz_0, ta1_z_zzzzz_yzz_1, ta1_z_zzzzz_yzzz_0, ta1_z_zzzzz_yzzz_1, ta1_z_zzzzz_zzz_0, ta1_z_zzzzz_zzz_1, ta1_z_zzzzz_zzzz_0, ta1_z_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzzz_xxxx_0[i] = ta1_z_zzzzz_xxxx_0[i] * pa_y[i] - ta1_z_zzzzz_xxxx_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxxy_0[i] = ta1_z_zzzzz_xxx_0[i] * fe_0 - ta1_z_zzzzz_xxx_1[i] * fe_0 + ta1_z_zzzzz_xxxy_0[i] * pa_y[i] - ta1_z_zzzzz_xxxy_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxxz_0[i] = ta1_z_zzzzz_xxxz_0[i] * pa_y[i] - ta1_z_zzzzz_xxxz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxyy_0[i] = 2.0 * ta1_z_zzzzz_xxy_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xxy_1[i] * fe_0 + ta1_z_zzzzz_xxyy_0[i] * pa_y[i] - ta1_z_zzzzz_xxyy_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxyz_0[i] = ta1_z_zzzzz_xxz_0[i] * fe_0 - ta1_z_zzzzz_xxz_1[i] * fe_0 + ta1_z_zzzzz_xxyz_0[i] * pa_y[i] - ta1_z_zzzzz_xxyz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxzz_0[i] = ta1_z_zzzzz_xxzz_0[i] * pa_y[i] - ta1_z_zzzzz_xxzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xyyy_0[i] = 3.0 * ta1_z_zzzzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_xyy_1[i] * fe_0 + ta1_z_zzzzz_xyyy_0[i] * pa_y[i] - ta1_z_zzzzz_xyyy_1[i] * pc_y[i];

        ta1_z_yzzzzz_xyyz_0[i] = 2.0 * ta1_z_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xyz_1[i] * fe_0 + ta1_z_zzzzz_xyyz_0[i] * pa_y[i] - ta1_z_zzzzz_xyyz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xyzz_0[i] = ta1_z_zzzzz_xzz_0[i] * fe_0 - ta1_z_zzzzz_xzz_1[i] * fe_0 + ta1_z_zzzzz_xyzz_0[i] * pa_y[i] - ta1_z_zzzzz_xyzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xzzz_0[i] = ta1_z_zzzzz_xzzz_0[i] * pa_y[i] - ta1_z_zzzzz_xzzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_yyyy_0[i] = 4.0 * ta1_z_zzzzz_yyy_0[i] * fe_0 - 4.0 * ta1_z_zzzzz_yyy_1[i] * fe_0 + ta1_z_zzzzz_yyyy_0[i] * pa_y[i] - ta1_z_zzzzz_yyyy_1[i] * pc_y[i];

        ta1_z_yzzzzz_yyyz_0[i] = 3.0 * ta1_z_zzzzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_yyz_1[i] * fe_0 + ta1_z_zzzzz_yyyz_0[i] * pa_y[i] - ta1_z_zzzzz_yyyz_1[i] * pc_y[i];

        ta1_z_yzzzzz_yyzz_0[i] = 2.0 * ta1_z_zzzzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_yzz_1[i] * fe_0 + ta1_z_zzzzz_yyzz_0[i] * pa_y[i] - ta1_z_zzzzz_yyzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_yzzz_0[i] = ta1_z_zzzzz_zzz_0[i] * fe_0 - ta1_z_zzzzz_zzz_1[i] * fe_0 + ta1_z_zzzzz_yzzz_0[i] * pa_y[i] - ta1_z_zzzzz_yzzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_zzzz_0[i] = ta1_z_zzzzz_zzzz_0[i] * pa_y[i] - ta1_z_zzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1245-1260 components of targeted buffer : IG

    auto ta1_z_zzzzzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1245);

    auto ta1_z_zzzzzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1246);

    auto ta1_z_zzzzzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1247);

    auto ta1_z_zzzzzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1248);

    auto ta1_z_zzzzzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1249);

    auto ta1_z_zzzzzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1250);

    auto ta1_z_zzzzzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1251);

    auto ta1_z_zzzzzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1252);

    auto ta1_z_zzzzzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1253);

    auto ta1_z_zzzzzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1254);

    auto ta1_z_zzzzzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1255);

    auto ta1_z_zzzzzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1256);

    auto ta1_z_zzzzzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1257);

    auto ta1_z_zzzzzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1258);

    auto ta1_z_zzzzzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_ig + 1259);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zzzz_xxxx_0, ta1_z_zzzz_xxxx_1, ta1_z_zzzz_xxxy_0, ta1_z_zzzz_xxxy_1, ta1_z_zzzz_xxxz_0, ta1_z_zzzz_xxxz_1, ta1_z_zzzz_xxyy_0, ta1_z_zzzz_xxyy_1, ta1_z_zzzz_xxyz_0, ta1_z_zzzz_xxyz_1, ta1_z_zzzz_xxzz_0, ta1_z_zzzz_xxzz_1, ta1_z_zzzz_xyyy_0, ta1_z_zzzz_xyyy_1, ta1_z_zzzz_xyyz_0, ta1_z_zzzz_xyyz_1, ta1_z_zzzz_xyzz_0, ta1_z_zzzz_xyzz_1, ta1_z_zzzz_xzzz_0, ta1_z_zzzz_xzzz_1, ta1_z_zzzz_yyyy_0, ta1_z_zzzz_yyyy_1, ta1_z_zzzz_yyyz_0, ta1_z_zzzz_yyyz_1, ta1_z_zzzz_yyzz_0, ta1_z_zzzz_yyzz_1, ta1_z_zzzz_yzzz_0, ta1_z_zzzz_yzzz_1, ta1_z_zzzz_zzzz_0, ta1_z_zzzz_zzzz_1, ta1_z_zzzzz_xxx_0, ta1_z_zzzzz_xxx_1, ta1_z_zzzzz_xxxx_0, ta1_z_zzzzz_xxxx_1, ta1_z_zzzzz_xxxy_0, ta1_z_zzzzz_xxxy_1, ta1_z_zzzzz_xxxz_0, ta1_z_zzzzz_xxxz_1, ta1_z_zzzzz_xxy_0, ta1_z_zzzzz_xxy_1, ta1_z_zzzzz_xxyy_0, ta1_z_zzzzz_xxyy_1, ta1_z_zzzzz_xxyz_0, ta1_z_zzzzz_xxyz_1, ta1_z_zzzzz_xxz_0, ta1_z_zzzzz_xxz_1, ta1_z_zzzzz_xxzz_0, ta1_z_zzzzz_xxzz_1, ta1_z_zzzzz_xyy_0, ta1_z_zzzzz_xyy_1, ta1_z_zzzzz_xyyy_0, ta1_z_zzzzz_xyyy_1, ta1_z_zzzzz_xyyz_0, ta1_z_zzzzz_xyyz_1, ta1_z_zzzzz_xyz_0, ta1_z_zzzzz_xyz_1, ta1_z_zzzzz_xyzz_0, ta1_z_zzzzz_xyzz_1, ta1_z_zzzzz_xzz_0, ta1_z_zzzzz_xzz_1, ta1_z_zzzzz_xzzz_0, ta1_z_zzzzz_xzzz_1, ta1_z_zzzzz_yyy_0, ta1_z_zzzzz_yyy_1, ta1_z_zzzzz_yyyy_0, ta1_z_zzzzz_yyyy_1, ta1_z_zzzzz_yyyz_0, ta1_z_zzzzz_yyyz_1, ta1_z_zzzzz_yyz_0, ta1_z_zzzzz_yyz_1, ta1_z_zzzzz_yyzz_0, ta1_z_zzzzz_yyzz_1, ta1_z_zzzzz_yzz_0, ta1_z_zzzzz_yzz_1, ta1_z_zzzzz_yzzz_0, ta1_z_zzzzz_yzzz_1, ta1_z_zzzzz_zzz_0, ta1_z_zzzzz_zzz_1, ta1_z_zzzzz_zzzz_0, ta1_z_zzzzz_zzzz_1, ta1_z_zzzzzz_xxxx_0, ta1_z_zzzzzz_xxxy_0, ta1_z_zzzzzz_xxxz_0, ta1_z_zzzzzz_xxyy_0, ta1_z_zzzzzz_xxyz_0, ta1_z_zzzzzz_xxzz_0, ta1_z_zzzzzz_xyyy_0, ta1_z_zzzzzz_xyyz_0, ta1_z_zzzzzz_xyzz_0, ta1_z_zzzzzz_xzzz_0, ta1_z_zzzzzz_yyyy_0, ta1_z_zzzzzz_yyyz_0, ta1_z_zzzzzz_yyzz_0, ta1_z_zzzzzz_yzzz_0, ta1_z_zzzzzz_zzzz_0, ta_zzzzz_xxxx_1, ta_zzzzz_xxxy_1, ta_zzzzz_xxxz_1, ta_zzzzz_xxyy_1, ta_zzzzz_xxyz_1, ta_zzzzz_xxzz_1, ta_zzzzz_xyyy_1, ta_zzzzz_xyyz_1, ta_zzzzz_xyzz_1, ta_zzzzz_xzzz_1, ta_zzzzz_yyyy_1, ta_zzzzz_yyyz_1, ta_zzzzz_yyzz_1, ta_zzzzz_yzzz_1, ta_zzzzz_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzzz_xxxx_0[i] = 5.0 * ta1_z_zzzz_xxxx_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxxx_1[i] * fe_0 + ta_zzzzz_xxxx_1[i] + ta1_z_zzzzz_xxxx_0[i] * pa_z[i] - ta1_z_zzzzz_xxxx_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxxy_0[i] = 5.0 * ta1_z_zzzz_xxxy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxxy_1[i] * fe_0 + ta_zzzzz_xxxy_1[i] + ta1_z_zzzzz_xxxy_0[i] * pa_z[i] - ta1_z_zzzzz_xxxy_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxxz_0[i] = 5.0 * ta1_z_zzzz_xxxz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxxz_1[i] * fe_0 + ta1_z_zzzzz_xxx_0[i] * fe_0 - ta1_z_zzzzz_xxx_1[i] * fe_0 + ta_zzzzz_xxxz_1[i] + ta1_z_zzzzz_xxxz_0[i] * pa_z[i] - ta1_z_zzzzz_xxxz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxyy_0[i] = 5.0 * ta1_z_zzzz_xxyy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxyy_1[i] * fe_0 + ta_zzzzz_xxyy_1[i] + ta1_z_zzzzz_xxyy_0[i] * pa_z[i] - ta1_z_zzzzz_xxyy_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxyz_0[i] = 5.0 * ta1_z_zzzz_xxyz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxyz_1[i] * fe_0 + ta1_z_zzzzz_xxy_0[i] * fe_0 - ta1_z_zzzzz_xxy_1[i] * fe_0 + ta_zzzzz_xxyz_1[i] + ta1_z_zzzzz_xxyz_0[i] * pa_z[i] - ta1_z_zzzzz_xxyz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxzz_0[i] = 5.0 * ta1_z_zzzz_xxzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta1_z_zzzzz_xxz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xxz_1[i] * fe_0 + ta_zzzzz_xxzz_1[i] + ta1_z_zzzzz_xxzz_0[i] * pa_z[i] - ta1_z_zzzzz_xxzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xyyy_0[i] = 5.0 * ta1_z_zzzz_xyyy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xyyy_1[i] * fe_0 + ta_zzzzz_xyyy_1[i] + ta1_z_zzzzz_xyyy_0[i] * pa_z[i] - ta1_z_zzzzz_xyyy_1[i] * pc_z[i];

        ta1_z_zzzzzz_xyyz_0[i] = 5.0 * ta1_z_zzzz_xyyz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xyyz_1[i] * fe_0 + ta1_z_zzzzz_xyy_0[i] * fe_0 - ta1_z_zzzzz_xyy_1[i] * fe_0 + ta_zzzzz_xyyz_1[i] + ta1_z_zzzzz_xyyz_0[i] * pa_z[i] - ta1_z_zzzzz_xyyz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xyzz_0[i] = 5.0 * ta1_z_zzzz_xyzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xyzz_1[i] * fe_0 + 2.0 * ta1_z_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xyz_1[i] * fe_0 + ta_zzzzz_xyzz_1[i] + ta1_z_zzzzz_xyzz_0[i] * pa_z[i] - ta1_z_zzzzz_xyzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xzzz_0[i] = 5.0 * ta1_z_zzzz_xzzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xzzz_1[i] * fe_0 + 3.0 * ta1_z_zzzzz_xzz_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_xzz_1[i] * fe_0 + ta_zzzzz_xzzz_1[i] + ta1_z_zzzzz_xzzz_0[i] * pa_z[i] - ta1_z_zzzzz_xzzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_yyyy_0[i] = 5.0 * ta1_z_zzzz_yyyy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yyyy_1[i] * fe_0 + ta_zzzzz_yyyy_1[i] + ta1_z_zzzzz_yyyy_0[i] * pa_z[i] - ta1_z_zzzzz_yyyy_1[i] * pc_z[i];

        ta1_z_zzzzzz_yyyz_0[i] = 5.0 * ta1_z_zzzz_yyyz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yyyz_1[i] * fe_0 + ta1_z_zzzzz_yyy_0[i] * fe_0 - ta1_z_zzzzz_yyy_1[i] * fe_0 + ta_zzzzz_yyyz_1[i] + ta1_z_zzzzz_yyyz_0[i] * pa_z[i] - ta1_z_zzzzz_yyyz_1[i] * pc_z[i];

        ta1_z_zzzzzz_yyzz_0[i] = 5.0 * ta1_z_zzzz_yyzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta1_z_zzzzz_yyz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_yyz_1[i] * fe_0 + ta_zzzzz_yyzz_1[i] + ta1_z_zzzzz_yyzz_0[i] * pa_z[i] - ta1_z_zzzzz_yyzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_yzzz_0[i] = 5.0 * ta1_z_zzzz_yzzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yzzz_1[i] * fe_0 + 3.0 * ta1_z_zzzzz_yzz_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_yzz_1[i] * fe_0 + ta_zzzzz_yzzz_1[i] + ta1_z_zzzzz_yzzz_0[i] * pa_z[i] - ta1_z_zzzzz_yzzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_zzzz_0[i] = 5.0 * ta1_z_zzzz_zzzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_zzzz_1[i] * fe_0 + 4.0 * ta1_z_zzzzz_zzz_0[i] * fe_0 - 4.0 * ta1_z_zzzzz_zzz_1[i] * fe_0 + ta_zzzzz_zzzz_1[i] + ta1_z_zzzzz_zzzz_0[i] * pa_z[i] - ta1_z_zzzzz_zzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

