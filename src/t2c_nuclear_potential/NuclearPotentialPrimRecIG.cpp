#include "NuclearPotentialPrimRecIG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_ig(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_ig,
                               const size_t              idx_npot_0_gg,
                               const size_t              idx_npot_1_gg,
                               const size_t              idx_npot_0_hf,
                               const size_t              idx_npot_1_hf,
                               const size_t              idx_npot_0_hg,
                               const size_t              idx_npot_1_hg,
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

    // Set up components of auxiliary buffer : GG

    auto ta_xxxx_xxxx_0 = pbuffer.data(idx_npot_0_gg);

    auto ta_xxxx_xxxy_0 = pbuffer.data(idx_npot_0_gg + 1);

    auto ta_xxxx_xxxz_0 = pbuffer.data(idx_npot_0_gg + 2);

    auto ta_xxxx_xxyy_0 = pbuffer.data(idx_npot_0_gg + 3);

    auto ta_xxxx_xxyz_0 = pbuffer.data(idx_npot_0_gg + 4);

    auto ta_xxxx_xxzz_0 = pbuffer.data(idx_npot_0_gg + 5);

    auto ta_xxxx_xyyy_0 = pbuffer.data(idx_npot_0_gg + 6);

    auto ta_xxxx_xyyz_0 = pbuffer.data(idx_npot_0_gg + 7);

    auto ta_xxxx_xyzz_0 = pbuffer.data(idx_npot_0_gg + 8);

    auto ta_xxxx_xzzz_0 = pbuffer.data(idx_npot_0_gg + 9);

    auto ta_xxxx_yyyy_0 = pbuffer.data(idx_npot_0_gg + 10);

    auto ta_xxxx_yyyz_0 = pbuffer.data(idx_npot_0_gg + 11);

    auto ta_xxxx_yyzz_0 = pbuffer.data(idx_npot_0_gg + 12);

    auto ta_xxxx_yzzz_0 = pbuffer.data(idx_npot_0_gg + 13);

    auto ta_xxxx_zzzz_0 = pbuffer.data(idx_npot_0_gg + 14);

    auto ta_xxxy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 15);

    auto ta_xxxy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 17);

    auto ta_xxxy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 20);

    auto ta_xxxy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 24);

    auto ta_xxxy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 25);

    auto ta_xxxy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 26);

    auto ta_xxxy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 27);

    auto ta_xxxy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 28);

    auto ta_xxxz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 30);

    auto ta_xxxz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 31);

    auto ta_xxxz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 32);

    auto ta_xxxz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 33);

    auto ta_xxxz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 35);

    auto ta_xxxz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 36);

    auto ta_xxxz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 39);

    auto ta_xxxz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 41);

    auto ta_xxxz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 42);

    auto ta_xxxz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 43);

    auto ta_xxxz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 44);

    auto ta_xxyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 45);

    auto ta_xxyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 46);

    auto ta_xxyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 47);

    auto ta_xxyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 48);

    auto ta_xxyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 49);

    auto ta_xxyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 50);

    auto ta_xxyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 51);

    auto ta_xxyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 52);

    auto ta_xxyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 53);

    auto ta_xxyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 54);

    auto ta_xxyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 55);

    auto ta_xxyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 56);

    auto ta_xxyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 57);

    auto ta_xxyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 58);

    auto ta_xxyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 59);

    auto ta_xxyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 62);

    auto ta_xxyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 65);

    auto ta_xxyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 69);

    auto ta_xxyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 71);

    auto ta_xxyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 72);

    auto ta_xxyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 73);

    auto ta_xxzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 75);

    auto ta_xxzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 76);

    auto ta_xxzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 77);

    auto ta_xxzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 78);

    auto ta_xxzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 79);

    auto ta_xxzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 80);

    auto ta_xxzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 81);

    auto ta_xxzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 82);

    auto ta_xxzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 83);

    auto ta_xxzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 84);

    auto ta_xxzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 85);

    auto ta_xxzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 86);

    auto ta_xxzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 87);

    auto ta_xxzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 88);

    auto ta_xxzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 89);

    auto ta_xyyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 91);

    auto ta_xyyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 93);

    auto ta_xyyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 94);

    auto ta_xyyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 96);

    auto ta_xyyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 97);

    auto ta_xyyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 98);

    auto ta_xyyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 100);

    auto ta_xyyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 101);

    auto ta_xyyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 102);

    auto ta_xyyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 103);

    auto ta_xyyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 104);

    auto ta_xyyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 116);

    auto ta_xyyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 117);

    auto ta_xyyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 118);

    auto ta_xyyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 119);

    auto ta_xyzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 130);

    auto ta_xyzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 131);

    auto ta_xyzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 132);

    auto ta_xyzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 133);

    auto ta_xzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 137);

    auto ta_xzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 139);

    auto ta_xzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 140);

    auto ta_xzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 142);

    auto ta_xzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 143);

    auto ta_xzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 144);

    auto ta_xzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 145);

    auto ta_xzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 146);

    auto ta_xzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 147);

    auto ta_xzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 148);

    auto ta_xzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 149);

    auto ta_yyyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 150);

    auto ta_yyyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 151);

    auto ta_yyyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 152);

    auto ta_yyyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 153);

    auto ta_yyyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 154);

    auto ta_yyyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 155);

    auto ta_yyyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 156);

    auto ta_yyyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 157);

    auto ta_yyyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 158);

    auto ta_yyyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 159);

    auto ta_yyyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 160);

    auto ta_yyyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 161);

    auto ta_yyyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 162);

    auto ta_yyyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 163);

    auto ta_yyyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 164);

    auto ta_yyyz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 166);

    auto ta_yyyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 167);

    auto ta_yyyz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 168);

    auto ta_yyyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 170);

    auto ta_yyyz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 171);

    auto ta_yyyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 174);

    auto ta_yyyz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 175);

    auto ta_yyyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 176);

    auto ta_yyyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 177);

    auto ta_yyyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 178);

    auto ta_yyyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 179);

    auto ta_yyzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 180);

    auto ta_yyzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 181);

    auto ta_yyzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 182);

    auto ta_yyzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 183);

    auto ta_yyzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 184);

    auto ta_yyzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 185);

    auto ta_yyzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 186);

    auto ta_yyzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 187);

    auto ta_yyzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 188);

    auto ta_yyzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 189);

    auto ta_yyzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 190);

    auto ta_yyzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 191);

    auto ta_yyzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 192);

    auto ta_yyzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 193);

    auto ta_yyzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 194);

    auto ta_yzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 195);

    auto ta_yzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 197);

    auto ta_yzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 199);

    auto ta_yzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 200);

    auto ta_yzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 202);

    auto ta_yzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 203);

    auto ta_yzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 204);

    auto ta_yzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 205);

    auto ta_yzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 206);

    auto ta_yzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 207);

    auto ta_yzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 208);

    auto ta_yzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 209);

    auto ta_zzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 210);

    auto ta_zzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 211);

    auto ta_zzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 212);

    auto ta_zzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 213);

    auto ta_zzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 214);

    auto ta_zzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 215);

    auto ta_zzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 216);

    auto ta_zzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 217);

    auto ta_zzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 218);

    auto ta_zzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 219);

    auto ta_zzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 220);

    auto ta_zzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 221);

    auto ta_zzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 222);

    auto ta_zzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 223);

    auto ta_zzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 224);

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

    auto ta_xxxy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 17);

    auto ta_xxxy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 20);

    auto ta_xxxy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 24);

    auto ta_xxxy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 25);

    auto ta_xxxy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 26);

    auto ta_xxxy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 27);

    auto ta_xxxy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 28);

    auto ta_xxxz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 30);

    auto ta_xxxz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 31);

    auto ta_xxxz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 32);

    auto ta_xxxz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 33);

    auto ta_xxxz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 35);

    auto ta_xxxz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 36);

    auto ta_xxxz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 39);

    auto ta_xxxz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 41);

    auto ta_xxxz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 42);

    auto ta_xxxz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 43);

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

    auto ta_xxyy_zzzz_1 = pbuffer.data(idx_npot_1_gg + 59);

    auto ta_xxyz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 62);

    auto ta_xxyz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 65);

    auto ta_xxyz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 69);

    auto ta_xxyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 71);

    auto ta_xxyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 72);

    auto ta_xxyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 73);

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

    auto ta_xxzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 85);

    auto ta_xxzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 86);

    auto ta_xxzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 87);

    auto ta_xxzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 88);

    auto ta_xxzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 89);

    auto ta_xyyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 91);

    auto ta_xyyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 93);

    auto ta_xyyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 94);

    auto ta_xyyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 96);

    auto ta_xyyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 97);

    auto ta_xyyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 98);

    auto ta_xyyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 100);

    auto ta_xyyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 101);

    auto ta_xyyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 102);

    auto ta_xyyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 103);

    auto ta_xyyy_zzzz_1 = pbuffer.data(idx_npot_1_gg + 104);

    auto ta_xyyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 116);

    auto ta_xyyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 117);

    auto ta_xyyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 118);

    auto ta_xyyz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 119);

    auto ta_xyzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 130);

    auto ta_xyzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 131);

    auto ta_xyzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 132);

    auto ta_xyzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 133);

    auto ta_xzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 137);

    auto ta_xzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 139);

    auto ta_xzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 140);

    auto ta_xzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 142);

    auto ta_xzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 143);

    auto ta_xzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 144);

    auto ta_xzzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 145);

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

    auto ta_yyyz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 167);

    auto ta_yyyz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 168);

    auto ta_yyyz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 170);

    auto ta_yyyz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 171);

    auto ta_yyyz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 174);

    auto ta_yyyz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 175);

    auto ta_yyyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 176);

    auto ta_yyyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 177);

    auto ta_yyyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 178);

    auto ta_yyyz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 179);

    auto ta_yyzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 180);

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

    auto ta_yzzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 195);

    auto ta_yzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 197);

    auto ta_yzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 199);

    auto ta_yzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 200);

    auto ta_yzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 202);

    auto ta_yzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 203);

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

    // Set up components of auxiliary buffer : HF

    auto ta_xxxxx_xxx_0 = pbuffer.data(idx_npot_0_hf);

    auto ta_xxxxx_xxy_0 = pbuffer.data(idx_npot_0_hf + 1);

    auto ta_xxxxx_xxz_0 = pbuffer.data(idx_npot_0_hf + 2);

    auto ta_xxxxx_xyy_0 = pbuffer.data(idx_npot_0_hf + 3);

    auto ta_xxxxx_xyz_0 = pbuffer.data(idx_npot_0_hf + 4);

    auto ta_xxxxx_xzz_0 = pbuffer.data(idx_npot_0_hf + 5);

    auto ta_xxxxx_yyy_0 = pbuffer.data(idx_npot_0_hf + 6);

    auto ta_xxxxx_yyz_0 = pbuffer.data(idx_npot_0_hf + 7);

    auto ta_xxxxx_yzz_0 = pbuffer.data(idx_npot_0_hf + 8);

    auto ta_xxxxx_zzz_0 = pbuffer.data(idx_npot_0_hf + 9);

    auto ta_xxxxz_xxz_0 = pbuffer.data(idx_npot_0_hf + 22);

    auto ta_xxxxz_xyz_0 = pbuffer.data(idx_npot_0_hf + 24);

    auto ta_xxxxz_xzz_0 = pbuffer.data(idx_npot_0_hf + 25);

    auto ta_xxxyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 31);

    auto ta_xxxyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 33);

    auto ta_xxxyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 34);

    auto ta_xxxyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 36);

    auto ta_xxxyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 37);

    auto ta_xxxyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 38);

    auto ta_xxxzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 50);

    auto ta_xxxzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 51);

    auto ta_xxxzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 52);

    auto ta_xxxzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 53);

    auto ta_xxxzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 54);

    auto ta_xxxzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 55);

    auto ta_xxxzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 57);

    auto ta_xxxzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 58);

    auto ta_xxxzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 59);

    auto ta_xxyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 61);

    auto ta_xxyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 63);

    auto ta_xxyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 64);

    auto ta_xxyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 66);

    auto ta_xxyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 67);

    auto ta_xxyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 68);

    auto ta_xxzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 90);

    auto ta_xxzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 91);

    auto ta_xxzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 92);

    auto ta_xxzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 93);

    auto ta_xxzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 94);

    auto ta_xxzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 95);

    auto ta_xxzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 97);

    auto ta_xxzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 98);

    auto ta_xxzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 99);

    auto ta_xyyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 101);

    auto ta_xyyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 103);

    auto ta_xyyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 104);

    auto ta_xyyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 106);

    auto ta_xyyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 107);

    auto ta_xyyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 108);

    auto ta_xyyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 124);

    auto ta_xyyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 127);

    auto ta_xyyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 128);

    auto ta_xzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 142);

    auto ta_xzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 144);

    auto ta_xzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 145);

    auto ta_xzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 147);

    auto ta_xzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 148);

    auto ta_xzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 149);

    auto ta_yyyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 150);

    auto ta_yyyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 151);

    auto ta_yyyyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 152);

    auto ta_yyyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 153);

    auto ta_yyyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 154);

    auto ta_yyyyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 155);

    auto ta_yyyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 156);

    auto ta_yyyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 157);

    auto ta_yyyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 158);

    auto ta_yyyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 159);

    auto ta_yyyyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 162);

    auto ta_yyyyz_xyz_0 = pbuffer.data(idx_npot_0_hf + 164);

    auto ta_yyyyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 165);

    auto ta_yyyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 167);

    auto ta_yyyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 168);

    auto ta_yyyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 169);

    auto ta_yyyzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 170);

    auto ta_yyyzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 171);

    auto ta_yyyzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 172);

    auto ta_yyyzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 173);

    auto ta_yyyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 174);

    auto ta_yyyzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 175);

    auto ta_yyyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 176);

    auto ta_yyyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 177);

    auto ta_yyyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 178);

    auto ta_yyyzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 179);

    auto ta_yyzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 180);

    auto ta_yyzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 181);

    auto ta_yyzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 182);

    auto ta_yyzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 183);

    auto ta_yyzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 184);

    auto ta_yyzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 185);

    auto ta_yyzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 186);

    auto ta_yyzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 187);

    auto ta_yyzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 188);

    auto ta_yyzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 189);

    auto ta_yzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 191);

    auto ta_yzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 192);

    auto ta_yzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 193);

    auto ta_yzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 194);

    auto ta_yzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 195);

    auto ta_yzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 196);

    auto ta_yzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 197);

    auto ta_yzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 198);

    auto ta_yzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 199);

    auto ta_zzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 200);

    auto ta_zzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 201);

    auto ta_zzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 202);

    auto ta_zzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 203);

    auto ta_zzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 204);

    auto ta_zzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 205);

    auto ta_zzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 206);

    auto ta_zzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 207);

    auto ta_zzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 208);

    auto ta_zzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 209);

    // Set up components of auxiliary buffer : HF

    auto ta_xxxxx_xxx_1 = pbuffer.data(idx_npot_1_hf);

    auto ta_xxxxx_xxy_1 = pbuffer.data(idx_npot_1_hf + 1);

    auto ta_xxxxx_xxz_1 = pbuffer.data(idx_npot_1_hf + 2);

    auto ta_xxxxx_xyy_1 = pbuffer.data(idx_npot_1_hf + 3);

    auto ta_xxxxx_xyz_1 = pbuffer.data(idx_npot_1_hf + 4);

    auto ta_xxxxx_xzz_1 = pbuffer.data(idx_npot_1_hf + 5);

    auto ta_xxxxx_yyy_1 = pbuffer.data(idx_npot_1_hf + 6);

    auto ta_xxxxx_yyz_1 = pbuffer.data(idx_npot_1_hf + 7);

    auto ta_xxxxx_yzz_1 = pbuffer.data(idx_npot_1_hf + 8);

    auto ta_xxxxx_zzz_1 = pbuffer.data(idx_npot_1_hf + 9);

    auto ta_xxxxz_xxz_1 = pbuffer.data(idx_npot_1_hf + 22);

    auto ta_xxxxz_xyz_1 = pbuffer.data(idx_npot_1_hf + 24);

    auto ta_xxxxz_xzz_1 = pbuffer.data(idx_npot_1_hf + 25);

    auto ta_xxxyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 31);

    auto ta_xxxyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 33);

    auto ta_xxxyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 34);

    auto ta_xxxyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 36);

    auto ta_xxxyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 37);

    auto ta_xxxyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 38);

    auto ta_xxxzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 50);

    auto ta_xxxzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 51);

    auto ta_xxxzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 52);

    auto ta_xxxzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 53);

    auto ta_xxxzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 54);

    auto ta_xxxzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 55);

    auto ta_xxxzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 57);

    auto ta_xxxzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 58);

    auto ta_xxxzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 59);

    auto ta_xxyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 61);

    auto ta_xxyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 63);

    auto ta_xxyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 64);

    auto ta_xxyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 66);

    auto ta_xxyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 67);

    auto ta_xxyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 68);

    auto ta_xxzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 90);

    auto ta_xxzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 91);

    auto ta_xxzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 92);

    auto ta_xxzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 93);

    auto ta_xxzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 94);

    auto ta_xxzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 95);

    auto ta_xxzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 97);

    auto ta_xxzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 98);

    auto ta_xxzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 99);

    auto ta_xyyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 101);

    auto ta_xyyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 103);

    auto ta_xyyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 104);

    auto ta_xyyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 106);

    auto ta_xyyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 107);

    auto ta_xyyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 108);

    auto ta_xyyzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 124);

    auto ta_xyyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 127);

    auto ta_xyyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 128);

    auto ta_xzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 142);

    auto ta_xzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 144);

    auto ta_xzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 145);

    auto ta_xzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 147);

    auto ta_xzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 148);

    auto ta_xzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 149);

    auto ta_yyyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 150);

    auto ta_yyyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 151);

    auto ta_yyyyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 152);

    auto ta_yyyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 153);

    auto ta_yyyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 154);

    auto ta_yyyyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 155);

    auto ta_yyyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 156);

    auto ta_yyyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 157);

    auto ta_yyyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 158);

    auto ta_yyyyy_zzz_1 = pbuffer.data(idx_npot_1_hf + 159);

    auto ta_yyyyz_xxz_1 = pbuffer.data(idx_npot_1_hf + 162);

    auto ta_yyyyz_xyz_1 = pbuffer.data(idx_npot_1_hf + 164);

    auto ta_yyyyz_xzz_1 = pbuffer.data(idx_npot_1_hf + 165);

    auto ta_yyyyz_yyz_1 = pbuffer.data(idx_npot_1_hf + 167);

    auto ta_yyyyz_yzz_1 = pbuffer.data(idx_npot_1_hf + 168);

    auto ta_yyyyz_zzz_1 = pbuffer.data(idx_npot_1_hf + 169);

    auto ta_yyyzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 170);

    auto ta_yyyzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 171);

    auto ta_yyyzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 172);

    auto ta_yyyzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 173);

    auto ta_yyyzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 174);

    auto ta_yyyzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 175);

    auto ta_yyyzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 176);

    auto ta_yyyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 177);

    auto ta_yyyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 178);

    auto ta_yyyzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 179);

    auto ta_yyzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 180);

    auto ta_yyzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 181);

    auto ta_yyzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 182);

    auto ta_yyzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 183);

    auto ta_yyzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 184);

    auto ta_yyzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 185);

    auto ta_yyzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 186);

    auto ta_yyzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 187);

    auto ta_yyzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 188);

    auto ta_yyzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 189);

    auto ta_yzzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 191);

    auto ta_yzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 192);

    auto ta_yzzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 193);

    auto ta_yzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 194);

    auto ta_yzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 195);

    auto ta_yzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 196);

    auto ta_yzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 197);

    auto ta_yzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 198);

    auto ta_yzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 199);

    auto ta_zzzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 200);

    auto ta_zzzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 201);

    auto ta_zzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 202);

    auto ta_zzzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 203);

    auto ta_zzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 204);

    auto ta_zzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 205);

    auto ta_zzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 206);

    auto ta_zzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 207);

    auto ta_zzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 208);

    auto ta_zzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 209);

    // Set up components of auxiliary buffer : HG

    auto ta_xxxxx_xxxx_0 = pbuffer.data(idx_npot_0_hg);

    auto ta_xxxxx_xxxy_0 = pbuffer.data(idx_npot_0_hg + 1);

    auto ta_xxxxx_xxxz_0 = pbuffer.data(idx_npot_0_hg + 2);

    auto ta_xxxxx_xxyy_0 = pbuffer.data(idx_npot_0_hg + 3);

    auto ta_xxxxx_xxyz_0 = pbuffer.data(idx_npot_0_hg + 4);

    auto ta_xxxxx_xxzz_0 = pbuffer.data(idx_npot_0_hg + 5);

    auto ta_xxxxx_xyyy_0 = pbuffer.data(idx_npot_0_hg + 6);

    auto ta_xxxxx_xyyz_0 = pbuffer.data(idx_npot_0_hg + 7);

    auto ta_xxxxx_xyzz_0 = pbuffer.data(idx_npot_0_hg + 8);

    auto ta_xxxxx_xzzz_0 = pbuffer.data(idx_npot_0_hg + 9);

    auto ta_xxxxx_yyyy_0 = pbuffer.data(idx_npot_0_hg + 10);

    auto ta_xxxxx_yyyz_0 = pbuffer.data(idx_npot_0_hg + 11);

    auto ta_xxxxx_yyzz_0 = pbuffer.data(idx_npot_0_hg + 12);

    auto ta_xxxxx_yzzz_0 = pbuffer.data(idx_npot_0_hg + 13);

    auto ta_xxxxx_zzzz_0 = pbuffer.data(idx_npot_0_hg + 14);

    auto ta_xxxxy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 15);

    auto ta_xxxxy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 16);

    auto ta_xxxxy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 17);

    auto ta_xxxxy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 18);

    auto ta_xxxxy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 20);

    auto ta_xxxxy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 21);

    auto ta_xxxxy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 24);

    auto ta_xxxxy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 25);

    auto ta_xxxxy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 26);

    auto ta_xxxxy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 27);

    auto ta_xxxxy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 28);

    auto ta_xxxxz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 30);

    auto ta_xxxxz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 31);

    auto ta_xxxxz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 32);

    auto ta_xxxxz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 33);

    auto ta_xxxxz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 34);

    auto ta_xxxxz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 35);

    auto ta_xxxxz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 36);

    auto ta_xxxxz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 37);

    auto ta_xxxxz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 38);

    auto ta_xxxxz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 39);

    auto ta_xxxxz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 41);

    auto ta_xxxxz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 42);

    auto ta_xxxxz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 43);

    auto ta_xxxxz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 44);

    auto ta_xxxyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 45);

    auto ta_xxxyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 46);

    auto ta_xxxyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 47);

    auto ta_xxxyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 48);

    auto ta_xxxyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 49);

    auto ta_xxxyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 50);

    auto ta_xxxyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 51);

    auto ta_xxxyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 52);

    auto ta_xxxyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 53);

    auto ta_xxxyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 54);

    auto ta_xxxyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 55);

    auto ta_xxxyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 56);

    auto ta_xxxyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 57);

    auto ta_xxxyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 58);

    auto ta_xxxyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 59);

    auto ta_xxxyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 62);

    auto ta_xxxyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 65);

    auto ta_xxxyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 69);

    auto ta_xxxyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 71);

    auto ta_xxxyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 72);

    auto ta_xxxyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 73);

    auto ta_xxxzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 75);

    auto ta_xxxzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 76);

    auto ta_xxxzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 77);

    auto ta_xxxzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 78);

    auto ta_xxxzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 79);

    auto ta_xxxzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 80);

    auto ta_xxxzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 81);

    auto ta_xxxzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 82);

    auto ta_xxxzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 83);

    auto ta_xxxzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 84);

    auto ta_xxxzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 85);

    auto ta_xxxzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 86);

    auto ta_xxxzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 87);

    auto ta_xxxzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 88);

    auto ta_xxxzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 89);

    auto ta_xxyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 90);

    auto ta_xxyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 91);

    auto ta_xxyyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 92);

    auto ta_xxyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 93);

    auto ta_xxyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 94);

    auto ta_xxyyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 95);

    auto ta_xxyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 96);

    auto ta_xxyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 97);

    auto ta_xxyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 98);

    auto ta_xxyyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 99);

    auto ta_xxyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 100);

    auto ta_xxyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 101);

    auto ta_xxyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 102);

    auto ta_xxyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 103);

    auto ta_xxyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 104);

    auto ta_xxyyz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 106);

    auto ta_xxyyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 107);

    auto ta_xxyyz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 108);

    auto ta_xxyyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 110);

    auto ta_xxyyz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 111);

    auto ta_xxyyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 114);

    auto ta_xxyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 116);

    auto ta_xxyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 117);

    auto ta_xxyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 118);

    auto ta_xxyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 119);

    auto ta_xxyzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 120);

    auto ta_xxyzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 122);

    auto ta_xxyzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 125);

    auto ta_xxyzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 129);

    auto ta_xxyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 130);

    auto ta_xxyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 131);

    auto ta_xxyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 132);

    auto ta_xxyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 133);

    auto ta_xxzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 135);

    auto ta_xxzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 136);

    auto ta_xxzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 137);

    auto ta_xxzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 138);

    auto ta_xxzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 139);

    auto ta_xxzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 140);

    auto ta_xxzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 141);

    auto ta_xxzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 142);

    auto ta_xxzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 143);

    auto ta_xxzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 144);

    auto ta_xxzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 145);

    auto ta_xxzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 146);

    auto ta_xxzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 147);

    auto ta_xxzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 148);

    auto ta_xxzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 149);

    auto ta_xyyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 150);

    auto ta_xyyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 151);

    auto ta_xyyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 153);

    auto ta_xyyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 154);

    auto ta_xyyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 156);

    auto ta_xyyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 157);

    auto ta_xyyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 158);

    auto ta_xyyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 160);

    auto ta_xyyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 161);

    auto ta_xyyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 162);

    auto ta_xyyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 163);

    auto ta_xyyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 164);

    auto ta_xyyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 176);

    auto ta_xyyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 177);

    auto ta_xyyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 178);

    auto ta_xyyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 179);

    auto ta_xyyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 184);

    auto ta_xyyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 187);

    auto ta_xyyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 188);

    auto ta_xyyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 190);

    auto ta_xyyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 191);

    auto ta_xyyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 192);

    auto ta_xyyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 193);

    auto ta_xyyzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 194);

    auto ta_xyzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 205);

    auto ta_xyzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 206);

    auto ta_xyzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 207);

    auto ta_xyzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 208);

    auto ta_xzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 210);

    auto ta_xzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 212);

    auto ta_xzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 214);

    auto ta_xzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 215);

    auto ta_xzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 217);

    auto ta_xzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 218);

    auto ta_xzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 219);

    auto ta_xzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 220);

    auto ta_xzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 221);

    auto ta_xzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 222);

    auto ta_xzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 223);

    auto ta_xzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 224);

    auto ta_yyyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 225);

    auto ta_yyyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 226);

    auto ta_yyyyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 227);

    auto ta_yyyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 228);

    auto ta_yyyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 229);

    auto ta_yyyyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 230);

    auto ta_yyyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 231);

    auto ta_yyyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 232);

    auto ta_yyyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 233);

    auto ta_yyyyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 234);

    auto ta_yyyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 235);

    auto ta_yyyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 236);

    auto ta_yyyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 237);

    auto ta_yyyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 238);

    auto ta_yyyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 239);

    auto ta_yyyyz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 241);

    auto ta_yyyyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 242);

    auto ta_yyyyz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 243);

    auto ta_yyyyz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 244);

    auto ta_yyyyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 245);

    auto ta_yyyyz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 246);

    auto ta_yyyyz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 247);

    auto ta_yyyyz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 248);

    auto ta_yyyyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 249);

    auto ta_yyyyz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 250);

    auto ta_yyyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 251);

    auto ta_yyyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 252);

    auto ta_yyyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 253);

    auto ta_yyyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 254);

    auto ta_yyyzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 255);

    auto ta_yyyzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 256);

    auto ta_yyyzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 257);

    auto ta_yyyzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 258);

    auto ta_yyyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 259);

    auto ta_yyyzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 260);

    auto ta_yyyzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 261);

    auto ta_yyyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 262);

    auto ta_yyyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 263);

    auto ta_yyyzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 264);

    auto ta_yyyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 265);

    auto ta_yyyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 266);

    auto ta_yyyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 267);

    auto ta_yyyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 268);

    auto ta_yyyzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 269);

    auto ta_yyzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 270);

    auto ta_yyzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 271);

    auto ta_yyzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 272);

    auto ta_yyzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 273);

    auto ta_yyzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 274);

    auto ta_yyzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 275);

    auto ta_yyzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 276);

    auto ta_yyzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 277);

    auto ta_yyzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 278);

    auto ta_yyzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 279);

    auto ta_yyzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 280);

    auto ta_yyzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 281);

    auto ta_yyzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 282);

    auto ta_yyzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 283);

    auto ta_yyzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 284);

    auto ta_yzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 285);

    auto ta_yzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 286);

    auto ta_yzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 287);

    auto ta_yzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 288);

    auto ta_yzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 289);

    auto ta_yzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 290);

    auto ta_yzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 291);

    auto ta_yzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 292);

    auto ta_yzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 293);

    auto ta_yzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 294);

    auto ta_yzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 295);

    auto ta_yzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 296);

    auto ta_yzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 297);

    auto ta_yzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 298);

    auto ta_yzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 299);

    auto ta_zzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 300);

    auto ta_zzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 301);

    auto ta_zzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 302);

    auto ta_zzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 303);

    auto ta_zzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 304);

    auto ta_zzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 305);

    auto ta_zzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 306);

    auto ta_zzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 307);

    auto ta_zzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 308);

    auto ta_zzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 309);

    auto ta_zzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 310);

    auto ta_zzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 311);

    auto ta_zzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 312);

    auto ta_zzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 313);

    auto ta_zzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 314);

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

    auto ta_xxxxy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 26);

    auto ta_xxxxy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 27);

    auto ta_xxxxy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 28);

    auto ta_xxxxz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 30);

    auto ta_xxxxz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 31);

    auto ta_xxxxz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 32);

    auto ta_xxxxz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 33);

    auto ta_xxxxz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 34);

    auto ta_xxxxz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 35);

    auto ta_xxxxz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 36);

    auto ta_xxxxz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 37);

    auto ta_xxxxz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 38);

    auto ta_xxxxz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 39);

    auto ta_xxxxz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 41);

    auto ta_xxxxz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 42);

    auto ta_xxxxz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 43);

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

    auto ta_xxxyy_zzzz_1 = pbuffer.data(idx_npot_1_hg + 59);

    auto ta_xxxyz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 62);

    auto ta_xxxyz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 65);

    auto ta_xxxyz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 69);

    auto ta_xxxyz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 71);

    auto ta_xxxyz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 72);

    auto ta_xxxyz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 73);

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

    auto ta_xxxzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 85);

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

    auto ta_xxyyy_zzzz_1 = pbuffer.data(idx_npot_1_hg + 104);

    auto ta_xxyyz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 106);

    auto ta_xxyyz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 107);

    auto ta_xxyyz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 108);

    auto ta_xxyyz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 110);

    auto ta_xxyyz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 111);

    auto ta_xxyyz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 114);

    auto ta_xxyyz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 116);

    auto ta_xxyyz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 117);

    auto ta_xxyyz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 118);

    auto ta_xxyyz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 119);

    auto ta_xxyzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 120);

    auto ta_xxyzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 122);

    auto ta_xxyzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 125);

    auto ta_xxyzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 129);

    auto ta_xxyzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 130);

    auto ta_xxyzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 131);

    auto ta_xxyzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 132);

    auto ta_xxyzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 133);

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

    auto ta_xxzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 145);

    auto ta_xxzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 146);

    auto ta_xxzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 147);

    auto ta_xxzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 148);

    auto ta_xxzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 149);

    auto ta_xyyyy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 150);

    auto ta_xyyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 151);

    auto ta_xyyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 153);

    auto ta_xyyyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 154);

    auto ta_xyyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 156);

    auto ta_xyyyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 157);

    auto ta_xyyyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 158);

    auto ta_xyyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 160);

    auto ta_xyyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 161);

    auto ta_xyyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 162);

    auto ta_xyyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 163);

    auto ta_xyyyy_zzzz_1 = pbuffer.data(idx_npot_1_hg + 164);

    auto ta_xyyyz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 176);

    auto ta_xyyyz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 177);

    auto ta_xyyyz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 178);

    auto ta_xyyyz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 179);

    auto ta_xyyzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 184);

    auto ta_xyyzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 187);

    auto ta_xyyzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 188);

    auto ta_xyyzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 190);

    auto ta_xyyzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 191);

    auto ta_xyyzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 192);

    auto ta_xyyzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 193);

    auto ta_xyyzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 194);

    auto ta_xyzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 205);

    auto ta_xyzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 206);

    auto ta_xyzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 207);

    auto ta_xyzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 208);

    auto ta_xzzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 210);

    auto ta_xzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 212);

    auto ta_xzzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 214);

    auto ta_xzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 215);

    auto ta_xzzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 217);

    auto ta_xzzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 218);

    auto ta_xzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 219);

    auto ta_xzzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 220);

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

    auto ta_yyyyz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 242);

    auto ta_yyyyz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 243);

    auto ta_yyyyz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 244);

    auto ta_yyyyz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 245);

    auto ta_yyyyz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 246);

    auto ta_yyyyz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 247);

    auto ta_yyyyz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 248);

    auto ta_yyyyz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 249);

    auto ta_yyyyz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 250);

    auto ta_yyyyz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 251);

    auto ta_yyyyz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 252);

    auto ta_yyyyz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 253);

    auto ta_yyyyz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 254);

    auto ta_yyyzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 255);

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

    auto ta_yyzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 270);

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

    auto ta_yzzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 285);

    auto ta_yzzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 286);

    auto ta_yzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 287);

    auto ta_yzzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 288);

    auto ta_yzzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 289);

    auto ta_yzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 290);

    auto ta_yzzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 291);

    auto ta_yzzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 292);

    auto ta_yzzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 293);

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

    // Set up 0-15 components of targeted buffer : IG

    auto ta_xxxxxx_xxxx_0 = pbuffer.data(idx_npot_0_ig);

    auto ta_xxxxxx_xxxy_0 = pbuffer.data(idx_npot_0_ig + 1);

    auto ta_xxxxxx_xxxz_0 = pbuffer.data(idx_npot_0_ig + 2);

    auto ta_xxxxxx_xxyy_0 = pbuffer.data(idx_npot_0_ig + 3);

    auto ta_xxxxxx_xxyz_0 = pbuffer.data(idx_npot_0_ig + 4);

    auto ta_xxxxxx_xxzz_0 = pbuffer.data(idx_npot_0_ig + 5);

    auto ta_xxxxxx_xyyy_0 = pbuffer.data(idx_npot_0_ig + 6);

    auto ta_xxxxxx_xyyz_0 = pbuffer.data(idx_npot_0_ig + 7);

    auto ta_xxxxxx_xyzz_0 = pbuffer.data(idx_npot_0_ig + 8);

    auto ta_xxxxxx_xzzz_0 = pbuffer.data(idx_npot_0_ig + 9);

    auto ta_xxxxxx_yyyy_0 = pbuffer.data(idx_npot_0_ig + 10);

    auto ta_xxxxxx_yyyz_0 = pbuffer.data(idx_npot_0_ig + 11);

    auto ta_xxxxxx_yyzz_0 = pbuffer.data(idx_npot_0_ig + 12);

    auto ta_xxxxxx_yzzz_0 = pbuffer.data(idx_npot_0_ig + 13);

    auto ta_xxxxxx_zzzz_0 = pbuffer.data(idx_npot_0_ig + 14);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xxxx_xxxx_0,   \
                             ta_xxxx_xxxx_1,   \
                             ta_xxxx_xxxy_0,   \
                             ta_xxxx_xxxy_1,   \
                             ta_xxxx_xxxz_0,   \
                             ta_xxxx_xxxz_1,   \
                             ta_xxxx_xxyy_0,   \
                             ta_xxxx_xxyy_1,   \
                             ta_xxxx_xxyz_0,   \
                             ta_xxxx_xxyz_1,   \
                             ta_xxxx_xxzz_0,   \
                             ta_xxxx_xxzz_1,   \
                             ta_xxxx_xyyy_0,   \
                             ta_xxxx_xyyy_1,   \
                             ta_xxxx_xyyz_0,   \
                             ta_xxxx_xyyz_1,   \
                             ta_xxxx_xyzz_0,   \
                             ta_xxxx_xyzz_1,   \
                             ta_xxxx_xzzz_0,   \
                             ta_xxxx_xzzz_1,   \
                             ta_xxxx_yyyy_0,   \
                             ta_xxxx_yyyy_1,   \
                             ta_xxxx_yyyz_0,   \
                             ta_xxxx_yyyz_1,   \
                             ta_xxxx_yyzz_0,   \
                             ta_xxxx_yyzz_1,   \
                             ta_xxxx_yzzz_0,   \
                             ta_xxxx_yzzz_1,   \
                             ta_xxxx_zzzz_0,   \
                             ta_xxxx_zzzz_1,   \
                             ta_xxxxx_xxx_0,   \
                             ta_xxxxx_xxx_1,   \
                             ta_xxxxx_xxxx_0,  \
                             ta_xxxxx_xxxx_1,  \
                             ta_xxxxx_xxxy_0,  \
                             ta_xxxxx_xxxy_1,  \
                             ta_xxxxx_xxxz_0,  \
                             ta_xxxxx_xxxz_1,  \
                             ta_xxxxx_xxy_0,   \
                             ta_xxxxx_xxy_1,   \
                             ta_xxxxx_xxyy_0,  \
                             ta_xxxxx_xxyy_1,  \
                             ta_xxxxx_xxyz_0,  \
                             ta_xxxxx_xxyz_1,  \
                             ta_xxxxx_xxz_0,   \
                             ta_xxxxx_xxz_1,   \
                             ta_xxxxx_xxzz_0,  \
                             ta_xxxxx_xxzz_1,  \
                             ta_xxxxx_xyy_0,   \
                             ta_xxxxx_xyy_1,   \
                             ta_xxxxx_xyyy_0,  \
                             ta_xxxxx_xyyy_1,  \
                             ta_xxxxx_xyyz_0,  \
                             ta_xxxxx_xyyz_1,  \
                             ta_xxxxx_xyz_0,   \
                             ta_xxxxx_xyz_1,   \
                             ta_xxxxx_xyzz_0,  \
                             ta_xxxxx_xyzz_1,  \
                             ta_xxxxx_xzz_0,   \
                             ta_xxxxx_xzz_1,   \
                             ta_xxxxx_xzzz_0,  \
                             ta_xxxxx_xzzz_1,  \
                             ta_xxxxx_yyy_0,   \
                             ta_xxxxx_yyy_1,   \
                             ta_xxxxx_yyyy_0,  \
                             ta_xxxxx_yyyy_1,  \
                             ta_xxxxx_yyyz_0,  \
                             ta_xxxxx_yyyz_1,  \
                             ta_xxxxx_yyz_0,   \
                             ta_xxxxx_yyz_1,   \
                             ta_xxxxx_yyzz_0,  \
                             ta_xxxxx_yyzz_1,  \
                             ta_xxxxx_yzz_0,   \
                             ta_xxxxx_yzz_1,   \
                             ta_xxxxx_yzzz_0,  \
                             ta_xxxxx_yzzz_1,  \
                             ta_xxxxx_zzz_0,   \
                             ta_xxxxx_zzz_1,   \
                             ta_xxxxx_zzzz_0,  \
                             ta_xxxxx_zzzz_1,  \
                             ta_xxxxxx_xxxx_0, \
                             ta_xxxxxx_xxxy_0, \
                             ta_xxxxxx_xxxz_0, \
                             ta_xxxxxx_xxyy_0, \
                             ta_xxxxxx_xxyz_0, \
                             ta_xxxxxx_xxzz_0, \
                             ta_xxxxxx_xyyy_0, \
                             ta_xxxxxx_xyyz_0, \
                             ta_xxxxxx_xyzz_0, \
                             ta_xxxxxx_xzzz_0, \
                             ta_xxxxxx_yyyy_0, \
                             ta_xxxxxx_yyyz_0, \
                             ta_xxxxxx_yyzz_0, \
                             ta_xxxxxx_yzzz_0, \
                             ta_xxxxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_xxxx_0[i] = 5.0 * ta_xxxx_xxxx_0[i] * fe_0 - 5.0 * ta_xxxx_xxxx_1[i] * fe_0 + 4.0 * ta_xxxxx_xxx_0[i] * fe_0 -
                              4.0 * ta_xxxxx_xxx_1[i] * fe_0 + ta_xxxxx_xxxx_0[i] * pa_x[i] - ta_xxxxx_xxxx_1[i] * pc_x[i];

        ta_xxxxxx_xxxy_0[i] = 5.0 * ta_xxxx_xxxy_0[i] * fe_0 - 5.0 * ta_xxxx_xxxy_1[i] * fe_0 + 3.0 * ta_xxxxx_xxy_0[i] * fe_0 -
                              3.0 * ta_xxxxx_xxy_1[i] * fe_0 + ta_xxxxx_xxxy_0[i] * pa_x[i] - ta_xxxxx_xxxy_1[i] * pc_x[i];

        ta_xxxxxx_xxxz_0[i] = 5.0 * ta_xxxx_xxxz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxz_1[i] * fe_0 + 3.0 * ta_xxxxx_xxz_0[i] * fe_0 -
                              3.0 * ta_xxxxx_xxz_1[i] * fe_0 + ta_xxxxx_xxxz_0[i] * pa_x[i] - ta_xxxxx_xxxz_1[i] * pc_x[i];

        ta_xxxxxx_xxyy_0[i] = 5.0 * ta_xxxx_xxyy_0[i] * fe_0 - 5.0 * ta_xxxx_xxyy_1[i] * fe_0 + 2.0 * ta_xxxxx_xyy_0[i] * fe_0 -
                              2.0 * ta_xxxxx_xyy_1[i] * fe_0 + ta_xxxxx_xxyy_0[i] * pa_x[i] - ta_xxxxx_xxyy_1[i] * pc_x[i];

        ta_xxxxxx_xxyz_0[i] = 5.0 * ta_xxxx_xxyz_0[i] * fe_0 - 5.0 * ta_xxxx_xxyz_1[i] * fe_0 + 2.0 * ta_xxxxx_xyz_0[i] * fe_0 -
                              2.0 * ta_xxxxx_xyz_1[i] * fe_0 + ta_xxxxx_xxyz_0[i] * pa_x[i] - ta_xxxxx_xxyz_1[i] * pc_x[i];

        ta_xxxxxx_xxzz_0[i] = 5.0 * ta_xxxx_xxzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxzz_1[i] * fe_0 + 2.0 * ta_xxxxx_xzz_0[i] * fe_0 -
                              2.0 * ta_xxxxx_xzz_1[i] * fe_0 + ta_xxxxx_xxzz_0[i] * pa_x[i] - ta_xxxxx_xxzz_1[i] * pc_x[i];

        ta_xxxxxx_xyyy_0[i] = 5.0 * ta_xxxx_xyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xyyy_1[i] * fe_0 + ta_xxxxx_yyy_0[i] * fe_0 - ta_xxxxx_yyy_1[i] * fe_0 +
                              ta_xxxxx_xyyy_0[i] * pa_x[i] - ta_xxxxx_xyyy_1[i] * pc_x[i];

        ta_xxxxxx_xyyz_0[i] = 5.0 * ta_xxxx_xyyz_0[i] * fe_0 - 5.0 * ta_xxxx_xyyz_1[i] * fe_0 + ta_xxxxx_yyz_0[i] * fe_0 - ta_xxxxx_yyz_1[i] * fe_0 +
                              ta_xxxxx_xyyz_0[i] * pa_x[i] - ta_xxxxx_xyyz_1[i] * pc_x[i];

        ta_xxxxxx_xyzz_0[i] = 5.0 * ta_xxxx_xyzz_0[i] * fe_0 - 5.0 * ta_xxxx_xyzz_1[i] * fe_0 + ta_xxxxx_yzz_0[i] * fe_0 - ta_xxxxx_yzz_1[i] * fe_0 +
                              ta_xxxxx_xyzz_0[i] * pa_x[i] - ta_xxxxx_xyzz_1[i] * pc_x[i];

        ta_xxxxxx_xzzz_0[i] = 5.0 * ta_xxxx_xzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xzzz_1[i] * fe_0 + ta_xxxxx_zzz_0[i] * fe_0 - ta_xxxxx_zzz_1[i] * fe_0 +
                              ta_xxxxx_xzzz_0[i] * pa_x[i] - ta_xxxxx_xzzz_1[i] * pc_x[i];

        ta_xxxxxx_yyyy_0[i] =
            5.0 * ta_xxxx_yyyy_0[i] * fe_0 - 5.0 * ta_xxxx_yyyy_1[i] * fe_0 + ta_xxxxx_yyyy_0[i] * pa_x[i] - ta_xxxxx_yyyy_1[i] * pc_x[i];

        ta_xxxxxx_yyyz_0[i] =
            5.0 * ta_xxxx_yyyz_0[i] * fe_0 - 5.0 * ta_xxxx_yyyz_1[i] * fe_0 + ta_xxxxx_yyyz_0[i] * pa_x[i] - ta_xxxxx_yyyz_1[i] * pc_x[i];

        ta_xxxxxx_yyzz_0[i] =
            5.0 * ta_xxxx_yyzz_0[i] * fe_0 - 5.0 * ta_xxxx_yyzz_1[i] * fe_0 + ta_xxxxx_yyzz_0[i] * pa_x[i] - ta_xxxxx_yyzz_1[i] * pc_x[i];

        ta_xxxxxx_yzzz_0[i] =
            5.0 * ta_xxxx_yzzz_0[i] * fe_0 - 5.0 * ta_xxxx_yzzz_1[i] * fe_0 + ta_xxxxx_yzzz_0[i] * pa_x[i] - ta_xxxxx_yzzz_1[i] * pc_x[i];

        ta_xxxxxx_zzzz_0[i] =
            5.0 * ta_xxxx_zzzz_0[i] * fe_0 - 5.0 * ta_xxxx_zzzz_1[i] * fe_0 + ta_xxxxx_zzzz_0[i] * pa_x[i] - ta_xxxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : IG

    auto ta_xxxxxy_xxxx_0 = pbuffer.data(idx_npot_0_ig + 15);

    auto ta_xxxxxy_xxxy_0 = pbuffer.data(idx_npot_0_ig + 16);

    auto ta_xxxxxy_xxxz_0 = pbuffer.data(idx_npot_0_ig + 17);

    auto ta_xxxxxy_xxyy_0 = pbuffer.data(idx_npot_0_ig + 18);

    auto ta_xxxxxy_xxyz_0 = pbuffer.data(idx_npot_0_ig + 19);

    auto ta_xxxxxy_xxzz_0 = pbuffer.data(idx_npot_0_ig + 20);

    auto ta_xxxxxy_xyyy_0 = pbuffer.data(idx_npot_0_ig + 21);

    auto ta_xxxxxy_xyyz_0 = pbuffer.data(idx_npot_0_ig + 22);

    auto ta_xxxxxy_xyzz_0 = pbuffer.data(idx_npot_0_ig + 23);

    auto ta_xxxxxy_xzzz_0 = pbuffer.data(idx_npot_0_ig + 24);

    auto ta_xxxxxy_yyyy_0 = pbuffer.data(idx_npot_0_ig + 25);

    auto ta_xxxxxy_yyyz_0 = pbuffer.data(idx_npot_0_ig + 26);

    auto ta_xxxxxy_yyzz_0 = pbuffer.data(idx_npot_0_ig + 27);

    auto ta_xxxxxy_yzzz_0 = pbuffer.data(idx_npot_0_ig + 28);

    auto ta_xxxxxy_zzzz_0 = pbuffer.data(idx_npot_0_ig + 29);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxxxx_xxx_0,   \
                             ta_xxxxx_xxx_1,   \
                             ta_xxxxx_xxxx_0,  \
                             ta_xxxxx_xxxx_1,  \
                             ta_xxxxx_xxxy_0,  \
                             ta_xxxxx_xxxy_1,  \
                             ta_xxxxx_xxxz_0,  \
                             ta_xxxxx_xxxz_1,  \
                             ta_xxxxx_xxy_0,   \
                             ta_xxxxx_xxy_1,   \
                             ta_xxxxx_xxyy_0,  \
                             ta_xxxxx_xxyy_1,  \
                             ta_xxxxx_xxyz_0,  \
                             ta_xxxxx_xxyz_1,  \
                             ta_xxxxx_xxz_0,   \
                             ta_xxxxx_xxz_1,   \
                             ta_xxxxx_xxzz_0,  \
                             ta_xxxxx_xxzz_1,  \
                             ta_xxxxx_xyy_0,   \
                             ta_xxxxx_xyy_1,   \
                             ta_xxxxx_xyyy_0,  \
                             ta_xxxxx_xyyy_1,  \
                             ta_xxxxx_xyyz_0,  \
                             ta_xxxxx_xyyz_1,  \
                             ta_xxxxx_xyz_0,   \
                             ta_xxxxx_xyz_1,   \
                             ta_xxxxx_xyzz_0,  \
                             ta_xxxxx_xyzz_1,  \
                             ta_xxxxx_xzz_0,   \
                             ta_xxxxx_xzz_1,   \
                             ta_xxxxx_xzzz_0,  \
                             ta_xxxxx_xzzz_1,  \
                             ta_xxxxx_zzzz_0,  \
                             ta_xxxxx_zzzz_1,  \
                             ta_xxxxxy_xxxx_0, \
                             ta_xxxxxy_xxxy_0, \
                             ta_xxxxxy_xxxz_0, \
                             ta_xxxxxy_xxyy_0, \
                             ta_xxxxxy_xxyz_0, \
                             ta_xxxxxy_xxzz_0, \
                             ta_xxxxxy_xyyy_0, \
                             ta_xxxxxy_xyyz_0, \
                             ta_xxxxxy_xyzz_0, \
                             ta_xxxxxy_xzzz_0, \
                             ta_xxxxxy_yyyy_0, \
                             ta_xxxxxy_yyyz_0, \
                             ta_xxxxxy_yyzz_0, \
                             ta_xxxxxy_yzzz_0, \
                             ta_xxxxxy_zzzz_0, \
                             ta_xxxxy_yyyy_0,  \
                             ta_xxxxy_yyyy_1,  \
                             ta_xxxxy_yyyz_0,  \
                             ta_xxxxy_yyyz_1,  \
                             ta_xxxxy_yyzz_0,  \
                             ta_xxxxy_yyzz_1,  \
                             ta_xxxxy_yzzz_0,  \
                             ta_xxxxy_yzzz_1,  \
                             ta_xxxy_yyyy_0,   \
                             ta_xxxy_yyyy_1,   \
                             ta_xxxy_yyyz_0,   \
                             ta_xxxy_yyyz_1,   \
                             ta_xxxy_yyzz_0,   \
                             ta_xxxy_yyzz_1,   \
                             ta_xxxy_yzzz_0,   \
                             ta_xxxy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxy_xxxx_0[i] = ta_xxxxx_xxxx_0[i] * pa_y[i] - ta_xxxxx_xxxx_1[i] * pc_y[i];

        ta_xxxxxy_xxxy_0[i] = ta_xxxxx_xxx_0[i] * fe_0 - ta_xxxxx_xxx_1[i] * fe_0 + ta_xxxxx_xxxy_0[i] * pa_y[i] - ta_xxxxx_xxxy_1[i] * pc_y[i];

        ta_xxxxxy_xxxz_0[i] = ta_xxxxx_xxxz_0[i] * pa_y[i] - ta_xxxxx_xxxz_1[i] * pc_y[i];

        ta_xxxxxy_xxyy_0[i] =
            2.0 * ta_xxxxx_xxy_0[i] * fe_0 - 2.0 * ta_xxxxx_xxy_1[i] * fe_0 + ta_xxxxx_xxyy_0[i] * pa_y[i] - ta_xxxxx_xxyy_1[i] * pc_y[i];

        ta_xxxxxy_xxyz_0[i] = ta_xxxxx_xxz_0[i] * fe_0 - ta_xxxxx_xxz_1[i] * fe_0 + ta_xxxxx_xxyz_0[i] * pa_y[i] - ta_xxxxx_xxyz_1[i] * pc_y[i];

        ta_xxxxxy_xxzz_0[i] = ta_xxxxx_xxzz_0[i] * pa_y[i] - ta_xxxxx_xxzz_1[i] * pc_y[i];

        ta_xxxxxy_xyyy_0[i] =
            3.0 * ta_xxxxx_xyy_0[i] * fe_0 - 3.0 * ta_xxxxx_xyy_1[i] * fe_0 + ta_xxxxx_xyyy_0[i] * pa_y[i] - ta_xxxxx_xyyy_1[i] * pc_y[i];

        ta_xxxxxy_xyyz_0[i] =
            2.0 * ta_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyz_1[i] * fe_0 + ta_xxxxx_xyyz_0[i] * pa_y[i] - ta_xxxxx_xyyz_1[i] * pc_y[i];

        ta_xxxxxy_xyzz_0[i] = ta_xxxxx_xzz_0[i] * fe_0 - ta_xxxxx_xzz_1[i] * fe_0 + ta_xxxxx_xyzz_0[i] * pa_y[i] - ta_xxxxx_xyzz_1[i] * pc_y[i];

        ta_xxxxxy_xzzz_0[i] = ta_xxxxx_xzzz_0[i] * pa_y[i] - ta_xxxxx_xzzz_1[i] * pc_y[i];

        ta_xxxxxy_yyyy_0[i] =
            4.0 * ta_xxxy_yyyy_0[i] * fe_0 - 4.0 * ta_xxxy_yyyy_1[i] * fe_0 + ta_xxxxy_yyyy_0[i] * pa_x[i] - ta_xxxxy_yyyy_1[i] * pc_x[i];

        ta_xxxxxy_yyyz_0[i] =
            4.0 * ta_xxxy_yyyz_0[i] * fe_0 - 4.0 * ta_xxxy_yyyz_1[i] * fe_0 + ta_xxxxy_yyyz_0[i] * pa_x[i] - ta_xxxxy_yyyz_1[i] * pc_x[i];

        ta_xxxxxy_yyzz_0[i] =
            4.0 * ta_xxxy_yyzz_0[i] * fe_0 - 4.0 * ta_xxxy_yyzz_1[i] * fe_0 + ta_xxxxy_yyzz_0[i] * pa_x[i] - ta_xxxxy_yyzz_1[i] * pc_x[i];

        ta_xxxxxy_yzzz_0[i] =
            4.0 * ta_xxxy_yzzz_0[i] * fe_0 - 4.0 * ta_xxxy_yzzz_1[i] * fe_0 + ta_xxxxy_yzzz_0[i] * pa_x[i] - ta_xxxxy_yzzz_1[i] * pc_x[i];

        ta_xxxxxy_zzzz_0[i] = ta_xxxxx_zzzz_0[i] * pa_y[i] - ta_xxxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : IG

    auto ta_xxxxxz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 30);

    auto ta_xxxxxz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 31);

    auto ta_xxxxxz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 32);

    auto ta_xxxxxz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 33);

    auto ta_xxxxxz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 34);

    auto ta_xxxxxz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 35);

    auto ta_xxxxxz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 36);

    auto ta_xxxxxz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 37);

    auto ta_xxxxxz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 38);

    auto ta_xxxxxz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 39);

    auto ta_xxxxxz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 40);

    auto ta_xxxxxz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 41);

    auto ta_xxxxxz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 42);

    auto ta_xxxxxz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 43);

    auto ta_xxxxxz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 44);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxxxx_xxx_0,   \
                             ta_xxxxx_xxx_1,   \
                             ta_xxxxx_xxxx_0,  \
                             ta_xxxxx_xxxx_1,  \
                             ta_xxxxx_xxxy_0,  \
                             ta_xxxxx_xxxy_1,  \
                             ta_xxxxx_xxxz_0,  \
                             ta_xxxxx_xxxz_1,  \
                             ta_xxxxx_xxy_0,   \
                             ta_xxxxx_xxy_1,   \
                             ta_xxxxx_xxyy_0,  \
                             ta_xxxxx_xxyy_1,  \
                             ta_xxxxx_xxyz_0,  \
                             ta_xxxxx_xxyz_1,  \
                             ta_xxxxx_xxz_0,   \
                             ta_xxxxx_xxz_1,   \
                             ta_xxxxx_xxzz_0,  \
                             ta_xxxxx_xxzz_1,  \
                             ta_xxxxx_xyy_0,   \
                             ta_xxxxx_xyy_1,   \
                             ta_xxxxx_xyyy_0,  \
                             ta_xxxxx_xyyy_1,  \
                             ta_xxxxx_xyyz_0,  \
                             ta_xxxxx_xyyz_1,  \
                             ta_xxxxx_xyz_0,   \
                             ta_xxxxx_xyz_1,   \
                             ta_xxxxx_xyzz_0,  \
                             ta_xxxxx_xyzz_1,  \
                             ta_xxxxx_xzz_0,   \
                             ta_xxxxx_xzz_1,   \
                             ta_xxxxx_xzzz_0,  \
                             ta_xxxxx_xzzz_1,  \
                             ta_xxxxx_yyyy_0,  \
                             ta_xxxxx_yyyy_1,  \
                             ta_xxxxxz_xxxx_0, \
                             ta_xxxxxz_xxxy_0, \
                             ta_xxxxxz_xxxz_0, \
                             ta_xxxxxz_xxyy_0, \
                             ta_xxxxxz_xxyz_0, \
                             ta_xxxxxz_xxzz_0, \
                             ta_xxxxxz_xyyy_0, \
                             ta_xxxxxz_xyyz_0, \
                             ta_xxxxxz_xyzz_0, \
                             ta_xxxxxz_xzzz_0, \
                             ta_xxxxxz_yyyy_0, \
                             ta_xxxxxz_yyyz_0, \
                             ta_xxxxxz_yyzz_0, \
                             ta_xxxxxz_yzzz_0, \
                             ta_xxxxxz_zzzz_0, \
                             ta_xxxxz_yyyz_0,  \
                             ta_xxxxz_yyyz_1,  \
                             ta_xxxxz_yyzz_0,  \
                             ta_xxxxz_yyzz_1,  \
                             ta_xxxxz_yzzz_0,  \
                             ta_xxxxz_yzzz_1,  \
                             ta_xxxxz_zzzz_0,  \
                             ta_xxxxz_zzzz_1,  \
                             ta_xxxz_yyyz_0,   \
                             ta_xxxz_yyyz_1,   \
                             ta_xxxz_yyzz_0,   \
                             ta_xxxz_yyzz_1,   \
                             ta_xxxz_yzzz_0,   \
                             ta_xxxz_yzzz_1,   \
                             ta_xxxz_zzzz_0,   \
                             ta_xxxz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxz_xxxx_0[i] = ta_xxxxx_xxxx_0[i] * pa_z[i] - ta_xxxxx_xxxx_1[i] * pc_z[i];

        ta_xxxxxz_xxxy_0[i] = ta_xxxxx_xxxy_0[i] * pa_z[i] - ta_xxxxx_xxxy_1[i] * pc_z[i];

        ta_xxxxxz_xxxz_0[i] = ta_xxxxx_xxx_0[i] * fe_0 - ta_xxxxx_xxx_1[i] * fe_0 + ta_xxxxx_xxxz_0[i] * pa_z[i] - ta_xxxxx_xxxz_1[i] * pc_z[i];

        ta_xxxxxz_xxyy_0[i] = ta_xxxxx_xxyy_0[i] * pa_z[i] - ta_xxxxx_xxyy_1[i] * pc_z[i];

        ta_xxxxxz_xxyz_0[i] = ta_xxxxx_xxy_0[i] * fe_0 - ta_xxxxx_xxy_1[i] * fe_0 + ta_xxxxx_xxyz_0[i] * pa_z[i] - ta_xxxxx_xxyz_1[i] * pc_z[i];

        ta_xxxxxz_xxzz_0[i] =
            2.0 * ta_xxxxx_xxz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxz_1[i] * fe_0 + ta_xxxxx_xxzz_0[i] * pa_z[i] - ta_xxxxx_xxzz_1[i] * pc_z[i];

        ta_xxxxxz_xyyy_0[i] = ta_xxxxx_xyyy_0[i] * pa_z[i] - ta_xxxxx_xyyy_1[i] * pc_z[i];

        ta_xxxxxz_xyyz_0[i] = ta_xxxxx_xyy_0[i] * fe_0 - ta_xxxxx_xyy_1[i] * fe_0 + ta_xxxxx_xyyz_0[i] * pa_z[i] - ta_xxxxx_xyyz_1[i] * pc_z[i];

        ta_xxxxxz_xyzz_0[i] =
            2.0 * ta_xxxxx_xyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyz_1[i] * fe_0 + ta_xxxxx_xyzz_0[i] * pa_z[i] - ta_xxxxx_xyzz_1[i] * pc_z[i];

        ta_xxxxxz_xzzz_0[i] =
            3.0 * ta_xxxxx_xzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xzz_1[i] * fe_0 + ta_xxxxx_xzzz_0[i] * pa_z[i] - ta_xxxxx_xzzz_1[i] * pc_z[i];

        ta_xxxxxz_yyyy_0[i] = ta_xxxxx_yyyy_0[i] * pa_z[i] - ta_xxxxx_yyyy_1[i] * pc_z[i];

        ta_xxxxxz_yyyz_0[i] =
            4.0 * ta_xxxz_yyyz_0[i] * fe_0 - 4.0 * ta_xxxz_yyyz_1[i] * fe_0 + ta_xxxxz_yyyz_0[i] * pa_x[i] - ta_xxxxz_yyyz_1[i] * pc_x[i];

        ta_xxxxxz_yyzz_0[i] =
            4.0 * ta_xxxz_yyzz_0[i] * fe_0 - 4.0 * ta_xxxz_yyzz_1[i] * fe_0 + ta_xxxxz_yyzz_0[i] * pa_x[i] - ta_xxxxz_yyzz_1[i] * pc_x[i];

        ta_xxxxxz_yzzz_0[i] =
            4.0 * ta_xxxz_yzzz_0[i] * fe_0 - 4.0 * ta_xxxz_yzzz_1[i] * fe_0 + ta_xxxxz_yzzz_0[i] * pa_x[i] - ta_xxxxz_yzzz_1[i] * pc_x[i];

        ta_xxxxxz_zzzz_0[i] =
            4.0 * ta_xxxz_zzzz_0[i] * fe_0 - 4.0 * ta_xxxz_zzzz_1[i] * fe_0 + ta_xxxxz_zzzz_0[i] * pa_x[i] - ta_xxxxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : IG

    auto ta_xxxxyy_xxxx_0 = pbuffer.data(idx_npot_0_ig + 45);

    auto ta_xxxxyy_xxxy_0 = pbuffer.data(idx_npot_0_ig + 46);

    auto ta_xxxxyy_xxxz_0 = pbuffer.data(idx_npot_0_ig + 47);

    auto ta_xxxxyy_xxyy_0 = pbuffer.data(idx_npot_0_ig + 48);

    auto ta_xxxxyy_xxyz_0 = pbuffer.data(idx_npot_0_ig + 49);

    auto ta_xxxxyy_xxzz_0 = pbuffer.data(idx_npot_0_ig + 50);

    auto ta_xxxxyy_xyyy_0 = pbuffer.data(idx_npot_0_ig + 51);

    auto ta_xxxxyy_xyyz_0 = pbuffer.data(idx_npot_0_ig + 52);

    auto ta_xxxxyy_xyzz_0 = pbuffer.data(idx_npot_0_ig + 53);

    auto ta_xxxxyy_xzzz_0 = pbuffer.data(idx_npot_0_ig + 54);

    auto ta_xxxxyy_yyyy_0 = pbuffer.data(idx_npot_0_ig + 55);

    auto ta_xxxxyy_yyyz_0 = pbuffer.data(idx_npot_0_ig + 56);

    auto ta_xxxxyy_yyzz_0 = pbuffer.data(idx_npot_0_ig + 57);

    auto ta_xxxxyy_yzzz_0 = pbuffer.data(idx_npot_0_ig + 58);

    auto ta_xxxxyy_zzzz_0 = pbuffer.data(idx_npot_0_ig + 59);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxxx_xxxx_0,   \
                             ta_xxxx_xxxx_1,   \
                             ta_xxxx_xxxz_0,   \
                             ta_xxxx_xxxz_1,   \
                             ta_xxxx_xxzz_0,   \
                             ta_xxxx_xxzz_1,   \
                             ta_xxxx_xzzz_0,   \
                             ta_xxxx_xzzz_1,   \
                             ta_xxxxy_xxxx_0,  \
                             ta_xxxxy_xxxx_1,  \
                             ta_xxxxy_xxxz_0,  \
                             ta_xxxxy_xxxz_1,  \
                             ta_xxxxy_xxzz_0,  \
                             ta_xxxxy_xxzz_1,  \
                             ta_xxxxy_xzzz_0,  \
                             ta_xxxxy_xzzz_1,  \
                             ta_xxxxyy_xxxx_0, \
                             ta_xxxxyy_xxxy_0, \
                             ta_xxxxyy_xxxz_0, \
                             ta_xxxxyy_xxyy_0, \
                             ta_xxxxyy_xxyz_0, \
                             ta_xxxxyy_xxzz_0, \
                             ta_xxxxyy_xyyy_0, \
                             ta_xxxxyy_xyyz_0, \
                             ta_xxxxyy_xyzz_0, \
                             ta_xxxxyy_xzzz_0, \
                             ta_xxxxyy_yyyy_0, \
                             ta_xxxxyy_yyyz_0, \
                             ta_xxxxyy_yyzz_0, \
                             ta_xxxxyy_yzzz_0, \
                             ta_xxxxyy_zzzz_0, \
                             ta_xxxyy_xxxy_0,  \
                             ta_xxxyy_xxxy_1,  \
                             ta_xxxyy_xxy_0,   \
                             ta_xxxyy_xxy_1,   \
                             ta_xxxyy_xxyy_0,  \
                             ta_xxxyy_xxyy_1,  \
                             ta_xxxyy_xxyz_0,  \
                             ta_xxxyy_xxyz_1,  \
                             ta_xxxyy_xyy_0,   \
                             ta_xxxyy_xyy_1,   \
                             ta_xxxyy_xyyy_0,  \
                             ta_xxxyy_xyyy_1,  \
                             ta_xxxyy_xyyz_0,  \
                             ta_xxxyy_xyyz_1,  \
                             ta_xxxyy_xyz_0,   \
                             ta_xxxyy_xyz_1,   \
                             ta_xxxyy_xyzz_0,  \
                             ta_xxxyy_xyzz_1,  \
                             ta_xxxyy_yyy_0,   \
                             ta_xxxyy_yyy_1,   \
                             ta_xxxyy_yyyy_0,  \
                             ta_xxxyy_yyyy_1,  \
                             ta_xxxyy_yyyz_0,  \
                             ta_xxxyy_yyyz_1,  \
                             ta_xxxyy_yyz_0,   \
                             ta_xxxyy_yyz_1,   \
                             ta_xxxyy_yyzz_0,  \
                             ta_xxxyy_yyzz_1,  \
                             ta_xxxyy_yzz_0,   \
                             ta_xxxyy_yzz_1,   \
                             ta_xxxyy_yzzz_0,  \
                             ta_xxxyy_yzzz_1,  \
                             ta_xxxyy_zzzz_0,  \
                             ta_xxxyy_zzzz_1,  \
                             ta_xxyy_xxxy_0,   \
                             ta_xxyy_xxxy_1,   \
                             ta_xxyy_xxyy_0,   \
                             ta_xxyy_xxyy_1,   \
                             ta_xxyy_xxyz_0,   \
                             ta_xxyy_xxyz_1,   \
                             ta_xxyy_xyyy_0,   \
                             ta_xxyy_xyyy_1,   \
                             ta_xxyy_xyyz_0,   \
                             ta_xxyy_xyyz_1,   \
                             ta_xxyy_xyzz_0,   \
                             ta_xxyy_xyzz_1,   \
                             ta_xxyy_yyyy_0,   \
                             ta_xxyy_yyyy_1,   \
                             ta_xxyy_yyyz_0,   \
                             ta_xxyy_yyyz_1,   \
                             ta_xxyy_yyzz_0,   \
                             ta_xxyy_yyzz_1,   \
                             ta_xxyy_yzzz_0,   \
                             ta_xxyy_yzzz_1,   \
                             ta_xxyy_zzzz_0,   \
                             ta_xxyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyy_xxxx_0[i] = ta_xxxx_xxxx_0[i] * fe_0 - ta_xxxx_xxxx_1[i] * fe_0 + ta_xxxxy_xxxx_0[i] * pa_y[i] - ta_xxxxy_xxxx_1[i] * pc_y[i];

        ta_xxxxyy_xxxy_0[i] = 3.0 * ta_xxyy_xxxy_0[i] * fe_0 - 3.0 * ta_xxyy_xxxy_1[i] * fe_0 + 3.0 * ta_xxxyy_xxy_0[i] * fe_0 -
                              3.0 * ta_xxxyy_xxy_1[i] * fe_0 + ta_xxxyy_xxxy_0[i] * pa_x[i] - ta_xxxyy_xxxy_1[i] * pc_x[i];

        ta_xxxxyy_xxxz_0[i] = ta_xxxx_xxxz_0[i] * fe_0 - ta_xxxx_xxxz_1[i] * fe_0 + ta_xxxxy_xxxz_0[i] * pa_y[i] - ta_xxxxy_xxxz_1[i] * pc_y[i];

        ta_xxxxyy_xxyy_0[i] = 3.0 * ta_xxyy_xxyy_0[i] * fe_0 - 3.0 * ta_xxyy_xxyy_1[i] * fe_0 + 2.0 * ta_xxxyy_xyy_0[i] * fe_0 -
                              2.0 * ta_xxxyy_xyy_1[i] * fe_0 + ta_xxxyy_xxyy_0[i] * pa_x[i] - ta_xxxyy_xxyy_1[i] * pc_x[i];

        ta_xxxxyy_xxyz_0[i] = 3.0 * ta_xxyy_xxyz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyz_1[i] * fe_0 + 2.0 * ta_xxxyy_xyz_0[i] * fe_0 -
                              2.0 * ta_xxxyy_xyz_1[i] * fe_0 + ta_xxxyy_xxyz_0[i] * pa_x[i] - ta_xxxyy_xxyz_1[i] * pc_x[i];

        ta_xxxxyy_xxzz_0[i] = ta_xxxx_xxzz_0[i] * fe_0 - ta_xxxx_xxzz_1[i] * fe_0 + ta_xxxxy_xxzz_0[i] * pa_y[i] - ta_xxxxy_xxzz_1[i] * pc_y[i];

        ta_xxxxyy_xyyy_0[i] = 3.0 * ta_xxyy_xyyy_0[i] * fe_0 - 3.0 * ta_xxyy_xyyy_1[i] * fe_0 + ta_xxxyy_yyy_0[i] * fe_0 - ta_xxxyy_yyy_1[i] * fe_0 +
                              ta_xxxyy_xyyy_0[i] * pa_x[i] - ta_xxxyy_xyyy_1[i] * pc_x[i];

        ta_xxxxyy_xyyz_0[i] = 3.0 * ta_xxyy_xyyz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyz_1[i] * fe_0 + ta_xxxyy_yyz_0[i] * fe_0 - ta_xxxyy_yyz_1[i] * fe_0 +
                              ta_xxxyy_xyyz_0[i] * pa_x[i] - ta_xxxyy_xyyz_1[i] * pc_x[i];

        ta_xxxxyy_xyzz_0[i] = 3.0 * ta_xxyy_xyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyzz_1[i] * fe_0 + ta_xxxyy_yzz_0[i] * fe_0 - ta_xxxyy_yzz_1[i] * fe_0 +
                              ta_xxxyy_xyzz_0[i] * pa_x[i] - ta_xxxyy_xyzz_1[i] * pc_x[i];

        ta_xxxxyy_xzzz_0[i] = ta_xxxx_xzzz_0[i] * fe_0 - ta_xxxx_xzzz_1[i] * fe_0 + ta_xxxxy_xzzz_0[i] * pa_y[i] - ta_xxxxy_xzzz_1[i] * pc_y[i];

        ta_xxxxyy_yyyy_0[i] =
            3.0 * ta_xxyy_yyyy_0[i] * fe_0 - 3.0 * ta_xxyy_yyyy_1[i] * fe_0 + ta_xxxyy_yyyy_0[i] * pa_x[i] - ta_xxxyy_yyyy_1[i] * pc_x[i];

        ta_xxxxyy_yyyz_0[i] =
            3.0 * ta_xxyy_yyyz_0[i] * fe_0 - 3.0 * ta_xxyy_yyyz_1[i] * fe_0 + ta_xxxyy_yyyz_0[i] * pa_x[i] - ta_xxxyy_yyyz_1[i] * pc_x[i];

        ta_xxxxyy_yyzz_0[i] =
            3.0 * ta_xxyy_yyzz_0[i] * fe_0 - 3.0 * ta_xxyy_yyzz_1[i] * fe_0 + ta_xxxyy_yyzz_0[i] * pa_x[i] - ta_xxxyy_yyzz_1[i] * pc_x[i];

        ta_xxxxyy_yzzz_0[i] =
            3.0 * ta_xxyy_yzzz_0[i] * fe_0 - 3.0 * ta_xxyy_yzzz_1[i] * fe_0 + ta_xxxyy_yzzz_0[i] * pa_x[i] - ta_xxxyy_yzzz_1[i] * pc_x[i];

        ta_xxxxyy_zzzz_0[i] =
            3.0 * ta_xxyy_zzzz_0[i] * fe_0 - 3.0 * ta_xxyy_zzzz_1[i] * fe_0 + ta_xxxyy_zzzz_0[i] * pa_x[i] - ta_xxxyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : IG

    auto ta_xxxxyz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 60);

    auto ta_xxxxyz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 61);

    auto ta_xxxxyz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 62);

    auto ta_xxxxyz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 63);

    auto ta_xxxxyz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 64);

    auto ta_xxxxyz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 65);

    auto ta_xxxxyz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 66);

    auto ta_xxxxyz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 67);

    auto ta_xxxxyz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 68);

    auto ta_xxxxyz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 69);

    auto ta_xxxxyz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 70);

    auto ta_xxxxyz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 71);

    auto ta_xxxxyz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 72);

    auto ta_xxxxyz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 73);

    auto ta_xxxxyz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 74);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta_xxxxy_xxxy_0,  \
                             ta_xxxxy_xxxy_1,  \
                             ta_xxxxy_xxyy_0,  \
                             ta_xxxxy_xxyy_1,  \
                             ta_xxxxy_xyyy_0,  \
                             ta_xxxxy_xyyy_1,  \
                             ta_xxxxy_yyyy_0,  \
                             ta_xxxxy_yyyy_1,  \
                             ta_xxxxyz_xxxx_0, \
                             ta_xxxxyz_xxxy_0, \
                             ta_xxxxyz_xxxz_0, \
                             ta_xxxxyz_xxyy_0, \
                             ta_xxxxyz_xxyz_0, \
                             ta_xxxxyz_xxzz_0, \
                             ta_xxxxyz_xyyy_0, \
                             ta_xxxxyz_xyyz_0, \
                             ta_xxxxyz_xyzz_0, \
                             ta_xxxxyz_xzzz_0, \
                             ta_xxxxyz_yyyy_0, \
                             ta_xxxxyz_yyyz_0, \
                             ta_xxxxyz_yyzz_0, \
                             ta_xxxxyz_yzzz_0, \
                             ta_xxxxyz_zzzz_0, \
                             ta_xxxxz_xxxx_0,  \
                             ta_xxxxz_xxxx_1,  \
                             ta_xxxxz_xxxz_0,  \
                             ta_xxxxz_xxxz_1,  \
                             ta_xxxxz_xxyz_0,  \
                             ta_xxxxz_xxyz_1,  \
                             ta_xxxxz_xxz_0,   \
                             ta_xxxxz_xxz_1,   \
                             ta_xxxxz_xxzz_0,  \
                             ta_xxxxz_xxzz_1,  \
                             ta_xxxxz_xyyz_0,  \
                             ta_xxxxz_xyyz_1,  \
                             ta_xxxxz_xyz_0,   \
                             ta_xxxxz_xyz_1,   \
                             ta_xxxxz_xyzz_0,  \
                             ta_xxxxz_xyzz_1,  \
                             ta_xxxxz_xzz_0,   \
                             ta_xxxxz_xzz_1,   \
                             ta_xxxxz_xzzz_0,  \
                             ta_xxxxz_xzzz_1,  \
                             ta_xxxxz_zzzz_0,  \
                             ta_xxxxz_zzzz_1,  \
                             ta_xxxyz_yyyz_0,  \
                             ta_xxxyz_yyyz_1,  \
                             ta_xxxyz_yyzz_0,  \
                             ta_xxxyz_yyzz_1,  \
                             ta_xxxyz_yzzz_0,  \
                             ta_xxxyz_yzzz_1,  \
                             ta_xxyz_yyyz_0,   \
                             ta_xxyz_yyyz_1,   \
                             ta_xxyz_yyzz_0,   \
                             ta_xxyz_yyzz_1,   \
                             ta_xxyz_yzzz_0,   \
                             ta_xxyz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyz_xxxx_0[i] = ta_xxxxz_xxxx_0[i] * pa_y[i] - ta_xxxxz_xxxx_1[i] * pc_y[i];

        ta_xxxxyz_xxxy_0[i] = ta_xxxxy_xxxy_0[i] * pa_z[i] - ta_xxxxy_xxxy_1[i] * pc_z[i];

        ta_xxxxyz_xxxz_0[i] = ta_xxxxz_xxxz_0[i] * pa_y[i] - ta_xxxxz_xxxz_1[i] * pc_y[i];

        ta_xxxxyz_xxyy_0[i] = ta_xxxxy_xxyy_0[i] * pa_z[i] - ta_xxxxy_xxyy_1[i] * pc_z[i];

        ta_xxxxyz_xxyz_0[i] = ta_xxxxz_xxz_0[i] * fe_0 - ta_xxxxz_xxz_1[i] * fe_0 + ta_xxxxz_xxyz_0[i] * pa_y[i] - ta_xxxxz_xxyz_1[i] * pc_y[i];

        ta_xxxxyz_xxzz_0[i] = ta_xxxxz_xxzz_0[i] * pa_y[i] - ta_xxxxz_xxzz_1[i] * pc_y[i];

        ta_xxxxyz_xyyy_0[i] = ta_xxxxy_xyyy_0[i] * pa_z[i] - ta_xxxxy_xyyy_1[i] * pc_z[i];

        ta_xxxxyz_xyyz_0[i] =
            2.0 * ta_xxxxz_xyz_0[i] * fe_0 - 2.0 * ta_xxxxz_xyz_1[i] * fe_0 + ta_xxxxz_xyyz_0[i] * pa_y[i] - ta_xxxxz_xyyz_1[i] * pc_y[i];

        ta_xxxxyz_xyzz_0[i] = ta_xxxxz_xzz_0[i] * fe_0 - ta_xxxxz_xzz_1[i] * fe_0 + ta_xxxxz_xyzz_0[i] * pa_y[i] - ta_xxxxz_xyzz_1[i] * pc_y[i];

        ta_xxxxyz_xzzz_0[i] = ta_xxxxz_xzzz_0[i] * pa_y[i] - ta_xxxxz_xzzz_1[i] * pc_y[i];

        ta_xxxxyz_yyyy_0[i] = ta_xxxxy_yyyy_0[i] * pa_z[i] - ta_xxxxy_yyyy_1[i] * pc_z[i];

        ta_xxxxyz_yyyz_0[i] =
            3.0 * ta_xxyz_yyyz_0[i] * fe_0 - 3.0 * ta_xxyz_yyyz_1[i] * fe_0 + ta_xxxyz_yyyz_0[i] * pa_x[i] - ta_xxxyz_yyyz_1[i] * pc_x[i];

        ta_xxxxyz_yyzz_0[i] =
            3.0 * ta_xxyz_yyzz_0[i] * fe_0 - 3.0 * ta_xxyz_yyzz_1[i] * fe_0 + ta_xxxyz_yyzz_0[i] * pa_x[i] - ta_xxxyz_yyzz_1[i] * pc_x[i];

        ta_xxxxyz_yzzz_0[i] =
            3.0 * ta_xxyz_yzzz_0[i] * fe_0 - 3.0 * ta_xxyz_yzzz_1[i] * fe_0 + ta_xxxyz_yzzz_0[i] * pa_x[i] - ta_xxxyz_yzzz_1[i] * pc_x[i];

        ta_xxxxyz_zzzz_0[i] = ta_xxxxz_zzzz_0[i] * pa_y[i] - ta_xxxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : IG

    auto ta_xxxxzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 75);

    auto ta_xxxxzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 76);

    auto ta_xxxxzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 77);

    auto ta_xxxxzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 78);

    auto ta_xxxxzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 79);

    auto ta_xxxxzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 80);

    auto ta_xxxxzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 81);

    auto ta_xxxxzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 82);

    auto ta_xxxxzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 83);

    auto ta_xxxxzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 84);

    auto ta_xxxxzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 85);

    auto ta_xxxxzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 86);

    auto ta_xxxxzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 87);

    auto ta_xxxxzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 88);

    auto ta_xxxxzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 89);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxxx_xxxx_0,   \
                             ta_xxxx_xxxx_1,   \
                             ta_xxxx_xxxy_0,   \
                             ta_xxxx_xxxy_1,   \
                             ta_xxxx_xxyy_0,   \
                             ta_xxxx_xxyy_1,   \
                             ta_xxxx_xyyy_0,   \
                             ta_xxxx_xyyy_1,   \
                             ta_xxxxz_xxxx_0,  \
                             ta_xxxxz_xxxx_1,  \
                             ta_xxxxz_xxxy_0,  \
                             ta_xxxxz_xxxy_1,  \
                             ta_xxxxz_xxyy_0,  \
                             ta_xxxxz_xxyy_1,  \
                             ta_xxxxz_xyyy_0,  \
                             ta_xxxxz_xyyy_1,  \
                             ta_xxxxzz_xxxx_0, \
                             ta_xxxxzz_xxxy_0, \
                             ta_xxxxzz_xxxz_0, \
                             ta_xxxxzz_xxyy_0, \
                             ta_xxxxzz_xxyz_0, \
                             ta_xxxxzz_xxzz_0, \
                             ta_xxxxzz_xyyy_0, \
                             ta_xxxxzz_xyyz_0, \
                             ta_xxxxzz_xyzz_0, \
                             ta_xxxxzz_xzzz_0, \
                             ta_xxxxzz_yyyy_0, \
                             ta_xxxxzz_yyyz_0, \
                             ta_xxxxzz_yyzz_0, \
                             ta_xxxxzz_yzzz_0, \
                             ta_xxxxzz_zzzz_0, \
                             ta_xxxzz_xxxz_0,  \
                             ta_xxxzz_xxxz_1,  \
                             ta_xxxzz_xxyz_0,  \
                             ta_xxxzz_xxyz_1,  \
                             ta_xxxzz_xxz_0,   \
                             ta_xxxzz_xxz_1,   \
                             ta_xxxzz_xxzz_0,  \
                             ta_xxxzz_xxzz_1,  \
                             ta_xxxzz_xyyz_0,  \
                             ta_xxxzz_xyyz_1,  \
                             ta_xxxzz_xyz_0,   \
                             ta_xxxzz_xyz_1,   \
                             ta_xxxzz_xyzz_0,  \
                             ta_xxxzz_xyzz_1,  \
                             ta_xxxzz_xzz_0,   \
                             ta_xxxzz_xzz_1,   \
                             ta_xxxzz_xzzz_0,  \
                             ta_xxxzz_xzzz_1,  \
                             ta_xxxzz_yyyy_0,  \
                             ta_xxxzz_yyyy_1,  \
                             ta_xxxzz_yyyz_0,  \
                             ta_xxxzz_yyyz_1,  \
                             ta_xxxzz_yyz_0,   \
                             ta_xxxzz_yyz_1,   \
                             ta_xxxzz_yyzz_0,  \
                             ta_xxxzz_yyzz_1,  \
                             ta_xxxzz_yzz_0,   \
                             ta_xxxzz_yzz_1,   \
                             ta_xxxzz_yzzz_0,  \
                             ta_xxxzz_yzzz_1,  \
                             ta_xxxzz_zzz_0,   \
                             ta_xxxzz_zzz_1,   \
                             ta_xxxzz_zzzz_0,  \
                             ta_xxxzz_zzzz_1,  \
                             ta_xxzz_xxxz_0,   \
                             ta_xxzz_xxxz_1,   \
                             ta_xxzz_xxyz_0,   \
                             ta_xxzz_xxyz_1,   \
                             ta_xxzz_xxzz_0,   \
                             ta_xxzz_xxzz_1,   \
                             ta_xxzz_xyyz_0,   \
                             ta_xxzz_xyyz_1,   \
                             ta_xxzz_xyzz_0,   \
                             ta_xxzz_xyzz_1,   \
                             ta_xxzz_xzzz_0,   \
                             ta_xxzz_xzzz_1,   \
                             ta_xxzz_yyyy_0,   \
                             ta_xxzz_yyyy_1,   \
                             ta_xxzz_yyyz_0,   \
                             ta_xxzz_yyyz_1,   \
                             ta_xxzz_yyzz_0,   \
                             ta_xxzz_yyzz_1,   \
                             ta_xxzz_yzzz_0,   \
                             ta_xxzz_yzzz_1,   \
                             ta_xxzz_zzzz_0,   \
                             ta_xxzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxzz_xxxx_0[i] = ta_xxxx_xxxx_0[i] * fe_0 - ta_xxxx_xxxx_1[i] * fe_0 + ta_xxxxz_xxxx_0[i] * pa_z[i] - ta_xxxxz_xxxx_1[i] * pc_z[i];

        ta_xxxxzz_xxxy_0[i] = ta_xxxx_xxxy_0[i] * fe_0 - ta_xxxx_xxxy_1[i] * fe_0 + ta_xxxxz_xxxy_0[i] * pa_z[i] - ta_xxxxz_xxxy_1[i] * pc_z[i];

        ta_xxxxzz_xxxz_0[i] = 3.0 * ta_xxzz_xxxz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxz_1[i] * fe_0 + 3.0 * ta_xxxzz_xxz_0[i] * fe_0 -
                              3.0 * ta_xxxzz_xxz_1[i] * fe_0 + ta_xxxzz_xxxz_0[i] * pa_x[i] - ta_xxxzz_xxxz_1[i] * pc_x[i];

        ta_xxxxzz_xxyy_0[i] = ta_xxxx_xxyy_0[i] * fe_0 - ta_xxxx_xxyy_1[i] * fe_0 + ta_xxxxz_xxyy_0[i] * pa_z[i] - ta_xxxxz_xxyy_1[i] * pc_z[i];

        ta_xxxxzz_xxyz_0[i] = 3.0 * ta_xxzz_xxyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyz_1[i] * fe_0 + 2.0 * ta_xxxzz_xyz_0[i] * fe_0 -
                              2.0 * ta_xxxzz_xyz_1[i] * fe_0 + ta_xxxzz_xxyz_0[i] * pa_x[i] - ta_xxxzz_xxyz_1[i] * pc_x[i];

        ta_xxxxzz_xxzz_0[i] = 3.0 * ta_xxzz_xxzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxzz_1[i] * fe_0 + 2.0 * ta_xxxzz_xzz_0[i] * fe_0 -
                              2.0 * ta_xxxzz_xzz_1[i] * fe_0 + ta_xxxzz_xxzz_0[i] * pa_x[i] - ta_xxxzz_xxzz_1[i] * pc_x[i];

        ta_xxxxzz_xyyy_0[i] = ta_xxxx_xyyy_0[i] * fe_0 - ta_xxxx_xyyy_1[i] * fe_0 + ta_xxxxz_xyyy_0[i] * pa_z[i] - ta_xxxxz_xyyy_1[i] * pc_z[i];

        ta_xxxxzz_xyyz_0[i] = 3.0 * ta_xxzz_xyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyz_1[i] * fe_0 + ta_xxxzz_yyz_0[i] * fe_0 - ta_xxxzz_yyz_1[i] * fe_0 +
                              ta_xxxzz_xyyz_0[i] * pa_x[i] - ta_xxxzz_xyyz_1[i] * pc_x[i];

        ta_xxxxzz_xyzz_0[i] = 3.0 * ta_xxzz_xyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyzz_1[i] * fe_0 + ta_xxxzz_yzz_0[i] * fe_0 - ta_xxxzz_yzz_1[i] * fe_0 +
                              ta_xxxzz_xyzz_0[i] * pa_x[i] - ta_xxxzz_xyzz_1[i] * pc_x[i];

        ta_xxxxzz_xzzz_0[i] = 3.0 * ta_xxzz_xzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xzzz_1[i] * fe_0 + ta_xxxzz_zzz_0[i] * fe_0 - ta_xxxzz_zzz_1[i] * fe_0 +
                              ta_xxxzz_xzzz_0[i] * pa_x[i] - ta_xxxzz_xzzz_1[i] * pc_x[i];

        ta_xxxxzz_yyyy_0[i] =
            3.0 * ta_xxzz_yyyy_0[i] * fe_0 - 3.0 * ta_xxzz_yyyy_1[i] * fe_0 + ta_xxxzz_yyyy_0[i] * pa_x[i] - ta_xxxzz_yyyy_1[i] * pc_x[i];

        ta_xxxxzz_yyyz_0[i] =
            3.0 * ta_xxzz_yyyz_0[i] * fe_0 - 3.0 * ta_xxzz_yyyz_1[i] * fe_0 + ta_xxxzz_yyyz_0[i] * pa_x[i] - ta_xxxzz_yyyz_1[i] * pc_x[i];

        ta_xxxxzz_yyzz_0[i] =
            3.0 * ta_xxzz_yyzz_0[i] * fe_0 - 3.0 * ta_xxzz_yyzz_1[i] * fe_0 + ta_xxxzz_yyzz_0[i] * pa_x[i] - ta_xxxzz_yyzz_1[i] * pc_x[i];

        ta_xxxxzz_yzzz_0[i] =
            3.0 * ta_xxzz_yzzz_0[i] * fe_0 - 3.0 * ta_xxzz_yzzz_1[i] * fe_0 + ta_xxxzz_yzzz_0[i] * pa_x[i] - ta_xxxzz_yzzz_1[i] * pc_x[i];

        ta_xxxxzz_zzzz_0[i] =
            3.0 * ta_xxzz_zzzz_0[i] * fe_0 - 3.0 * ta_xxzz_zzzz_1[i] * fe_0 + ta_xxxzz_zzzz_0[i] * pa_x[i] - ta_xxxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : IG

    auto ta_xxxyyy_xxxx_0 = pbuffer.data(idx_npot_0_ig + 90);

    auto ta_xxxyyy_xxxy_0 = pbuffer.data(idx_npot_0_ig + 91);

    auto ta_xxxyyy_xxxz_0 = pbuffer.data(idx_npot_0_ig + 92);

    auto ta_xxxyyy_xxyy_0 = pbuffer.data(idx_npot_0_ig + 93);

    auto ta_xxxyyy_xxyz_0 = pbuffer.data(idx_npot_0_ig + 94);

    auto ta_xxxyyy_xxzz_0 = pbuffer.data(idx_npot_0_ig + 95);

    auto ta_xxxyyy_xyyy_0 = pbuffer.data(idx_npot_0_ig + 96);

    auto ta_xxxyyy_xyyz_0 = pbuffer.data(idx_npot_0_ig + 97);

    auto ta_xxxyyy_xyzz_0 = pbuffer.data(idx_npot_0_ig + 98);

    auto ta_xxxyyy_xzzz_0 = pbuffer.data(idx_npot_0_ig + 99);

    auto ta_xxxyyy_yyyy_0 = pbuffer.data(idx_npot_0_ig + 100);

    auto ta_xxxyyy_yyyz_0 = pbuffer.data(idx_npot_0_ig + 101);

    auto ta_xxxyyy_yyzz_0 = pbuffer.data(idx_npot_0_ig + 102);

    auto ta_xxxyyy_yzzz_0 = pbuffer.data(idx_npot_0_ig + 103);

    auto ta_xxxyyy_zzzz_0 = pbuffer.data(idx_npot_0_ig + 104);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxxy_xxxx_0,   \
                             ta_xxxy_xxxx_1,   \
                             ta_xxxy_xxxz_0,   \
                             ta_xxxy_xxxz_1,   \
                             ta_xxxy_xxzz_0,   \
                             ta_xxxy_xxzz_1,   \
                             ta_xxxy_xzzz_0,   \
                             ta_xxxy_xzzz_1,   \
                             ta_xxxyy_xxxx_0,  \
                             ta_xxxyy_xxxx_1,  \
                             ta_xxxyy_xxxz_0,  \
                             ta_xxxyy_xxxz_1,  \
                             ta_xxxyy_xxzz_0,  \
                             ta_xxxyy_xxzz_1,  \
                             ta_xxxyy_xzzz_0,  \
                             ta_xxxyy_xzzz_1,  \
                             ta_xxxyyy_xxxx_0, \
                             ta_xxxyyy_xxxy_0, \
                             ta_xxxyyy_xxxz_0, \
                             ta_xxxyyy_xxyy_0, \
                             ta_xxxyyy_xxyz_0, \
                             ta_xxxyyy_xxzz_0, \
                             ta_xxxyyy_xyyy_0, \
                             ta_xxxyyy_xyyz_0, \
                             ta_xxxyyy_xyzz_0, \
                             ta_xxxyyy_xzzz_0, \
                             ta_xxxyyy_yyyy_0, \
                             ta_xxxyyy_yyyz_0, \
                             ta_xxxyyy_yyzz_0, \
                             ta_xxxyyy_yzzz_0, \
                             ta_xxxyyy_zzzz_0, \
                             ta_xxyyy_xxxy_0,  \
                             ta_xxyyy_xxxy_1,  \
                             ta_xxyyy_xxy_0,   \
                             ta_xxyyy_xxy_1,   \
                             ta_xxyyy_xxyy_0,  \
                             ta_xxyyy_xxyy_1,  \
                             ta_xxyyy_xxyz_0,  \
                             ta_xxyyy_xxyz_1,  \
                             ta_xxyyy_xyy_0,   \
                             ta_xxyyy_xyy_1,   \
                             ta_xxyyy_xyyy_0,  \
                             ta_xxyyy_xyyy_1,  \
                             ta_xxyyy_xyyz_0,  \
                             ta_xxyyy_xyyz_1,  \
                             ta_xxyyy_xyz_0,   \
                             ta_xxyyy_xyz_1,   \
                             ta_xxyyy_xyzz_0,  \
                             ta_xxyyy_xyzz_1,  \
                             ta_xxyyy_yyy_0,   \
                             ta_xxyyy_yyy_1,   \
                             ta_xxyyy_yyyy_0,  \
                             ta_xxyyy_yyyy_1,  \
                             ta_xxyyy_yyyz_0,  \
                             ta_xxyyy_yyyz_1,  \
                             ta_xxyyy_yyz_0,   \
                             ta_xxyyy_yyz_1,   \
                             ta_xxyyy_yyzz_0,  \
                             ta_xxyyy_yyzz_1,  \
                             ta_xxyyy_yzz_0,   \
                             ta_xxyyy_yzz_1,   \
                             ta_xxyyy_yzzz_0,  \
                             ta_xxyyy_yzzz_1,  \
                             ta_xxyyy_zzzz_0,  \
                             ta_xxyyy_zzzz_1,  \
                             ta_xyyy_xxxy_0,   \
                             ta_xyyy_xxxy_1,   \
                             ta_xyyy_xxyy_0,   \
                             ta_xyyy_xxyy_1,   \
                             ta_xyyy_xxyz_0,   \
                             ta_xyyy_xxyz_1,   \
                             ta_xyyy_xyyy_0,   \
                             ta_xyyy_xyyy_1,   \
                             ta_xyyy_xyyz_0,   \
                             ta_xyyy_xyyz_1,   \
                             ta_xyyy_xyzz_0,   \
                             ta_xyyy_xyzz_1,   \
                             ta_xyyy_yyyy_0,   \
                             ta_xyyy_yyyy_1,   \
                             ta_xyyy_yyyz_0,   \
                             ta_xyyy_yyyz_1,   \
                             ta_xyyy_yyzz_0,   \
                             ta_xyyy_yyzz_1,   \
                             ta_xyyy_yzzz_0,   \
                             ta_xyyy_yzzz_1,   \
                             ta_xyyy_zzzz_0,   \
                             ta_xyyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyy_xxxx_0[i] =
            2.0 * ta_xxxy_xxxx_0[i] * fe_0 - 2.0 * ta_xxxy_xxxx_1[i] * fe_0 + ta_xxxyy_xxxx_0[i] * pa_y[i] - ta_xxxyy_xxxx_1[i] * pc_y[i];

        ta_xxxyyy_xxxy_0[i] = 2.0 * ta_xyyy_xxxy_0[i] * fe_0 - 2.0 * ta_xyyy_xxxy_1[i] * fe_0 + 3.0 * ta_xxyyy_xxy_0[i] * fe_0 -
                              3.0 * ta_xxyyy_xxy_1[i] * fe_0 + ta_xxyyy_xxxy_0[i] * pa_x[i] - ta_xxyyy_xxxy_1[i] * pc_x[i];

        ta_xxxyyy_xxxz_0[i] =
            2.0 * ta_xxxy_xxxz_0[i] * fe_0 - 2.0 * ta_xxxy_xxxz_1[i] * fe_0 + ta_xxxyy_xxxz_0[i] * pa_y[i] - ta_xxxyy_xxxz_1[i] * pc_y[i];

        ta_xxxyyy_xxyy_0[i] = 2.0 * ta_xyyy_xxyy_0[i] * fe_0 - 2.0 * ta_xyyy_xxyy_1[i] * fe_0 + 2.0 * ta_xxyyy_xyy_0[i] * fe_0 -
                              2.0 * ta_xxyyy_xyy_1[i] * fe_0 + ta_xxyyy_xxyy_0[i] * pa_x[i] - ta_xxyyy_xxyy_1[i] * pc_x[i];

        ta_xxxyyy_xxyz_0[i] = 2.0 * ta_xyyy_xxyz_0[i] * fe_0 - 2.0 * ta_xyyy_xxyz_1[i] * fe_0 + 2.0 * ta_xxyyy_xyz_0[i] * fe_0 -
                              2.0 * ta_xxyyy_xyz_1[i] * fe_0 + ta_xxyyy_xxyz_0[i] * pa_x[i] - ta_xxyyy_xxyz_1[i] * pc_x[i];

        ta_xxxyyy_xxzz_0[i] =
            2.0 * ta_xxxy_xxzz_0[i] * fe_0 - 2.0 * ta_xxxy_xxzz_1[i] * fe_0 + ta_xxxyy_xxzz_0[i] * pa_y[i] - ta_xxxyy_xxzz_1[i] * pc_y[i];

        ta_xxxyyy_xyyy_0[i] = 2.0 * ta_xyyy_xyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xyyy_1[i] * fe_0 + ta_xxyyy_yyy_0[i] * fe_0 - ta_xxyyy_yyy_1[i] * fe_0 +
                              ta_xxyyy_xyyy_0[i] * pa_x[i] - ta_xxyyy_xyyy_1[i] * pc_x[i];

        ta_xxxyyy_xyyz_0[i] = 2.0 * ta_xyyy_xyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyz_1[i] * fe_0 + ta_xxyyy_yyz_0[i] * fe_0 - ta_xxyyy_yyz_1[i] * fe_0 +
                              ta_xxyyy_xyyz_0[i] * pa_x[i] - ta_xxyyy_xyyz_1[i] * pc_x[i];

        ta_xxxyyy_xyzz_0[i] = 2.0 * ta_xyyy_xyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyzz_1[i] * fe_0 + ta_xxyyy_yzz_0[i] * fe_0 - ta_xxyyy_yzz_1[i] * fe_0 +
                              ta_xxyyy_xyzz_0[i] * pa_x[i] - ta_xxyyy_xyzz_1[i] * pc_x[i];

        ta_xxxyyy_xzzz_0[i] =
            2.0 * ta_xxxy_xzzz_0[i] * fe_0 - 2.0 * ta_xxxy_xzzz_1[i] * fe_0 + ta_xxxyy_xzzz_0[i] * pa_y[i] - ta_xxxyy_xzzz_1[i] * pc_y[i];

        ta_xxxyyy_yyyy_0[i] =
            2.0 * ta_xyyy_yyyy_0[i] * fe_0 - 2.0 * ta_xyyy_yyyy_1[i] * fe_0 + ta_xxyyy_yyyy_0[i] * pa_x[i] - ta_xxyyy_yyyy_1[i] * pc_x[i];

        ta_xxxyyy_yyyz_0[i] =
            2.0 * ta_xyyy_yyyz_0[i] * fe_0 - 2.0 * ta_xyyy_yyyz_1[i] * fe_0 + ta_xxyyy_yyyz_0[i] * pa_x[i] - ta_xxyyy_yyyz_1[i] * pc_x[i];

        ta_xxxyyy_yyzz_0[i] =
            2.0 * ta_xyyy_yyzz_0[i] * fe_0 - 2.0 * ta_xyyy_yyzz_1[i] * fe_0 + ta_xxyyy_yyzz_0[i] * pa_x[i] - ta_xxyyy_yyzz_1[i] * pc_x[i];

        ta_xxxyyy_yzzz_0[i] =
            2.0 * ta_xyyy_yzzz_0[i] * fe_0 - 2.0 * ta_xyyy_yzzz_1[i] * fe_0 + ta_xxyyy_yzzz_0[i] * pa_x[i] - ta_xxyyy_yzzz_1[i] * pc_x[i];

        ta_xxxyyy_zzzz_0[i] =
            2.0 * ta_xyyy_zzzz_0[i] * fe_0 - 2.0 * ta_xyyy_zzzz_1[i] * fe_0 + ta_xxyyy_zzzz_0[i] * pa_x[i] - ta_xxyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : IG

    auto ta_xxxyyz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 105);

    auto ta_xxxyyz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 106);

    auto ta_xxxyyz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 107);

    auto ta_xxxyyz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 108);

    auto ta_xxxyyz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 109);

    auto ta_xxxyyz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 110);

    auto ta_xxxyyz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 111);

    auto ta_xxxyyz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 112);

    auto ta_xxxyyz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 113);

    auto ta_xxxyyz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 114);

    auto ta_xxxyyz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 115);

    auto ta_xxxyyz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 116);

    auto ta_xxxyyz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 117);

    auto ta_xxxyyz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 118);

    auto ta_xxxyyz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 119);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta_xxxyy_xxxx_0,  \
                             ta_xxxyy_xxxx_1,  \
                             ta_xxxyy_xxxy_0,  \
                             ta_xxxyy_xxxy_1,  \
                             ta_xxxyy_xxy_0,   \
                             ta_xxxyy_xxy_1,   \
                             ta_xxxyy_xxyy_0,  \
                             ta_xxxyy_xxyy_1,  \
                             ta_xxxyy_xxyz_0,  \
                             ta_xxxyy_xxyz_1,  \
                             ta_xxxyy_xyy_0,   \
                             ta_xxxyy_xyy_1,   \
                             ta_xxxyy_xyyy_0,  \
                             ta_xxxyy_xyyy_1,  \
                             ta_xxxyy_xyyz_0,  \
                             ta_xxxyy_xyyz_1,  \
                             ta_xxxyy_xyz_0,   \
                             ta_xxxyy_xyz_1,   \
                             ta_xxxyy_xyzz_0,  \
                             ta_xxxyy_xyzz_1,  \
                             ta_xxxyy_yyyy_0,  \
                             ta_xxxyy_yyyy_1,  \
                             ta_xxxyyz_xxxx_0, \
                             ta_xxxyyz_xxxy_0, \
                             ta_xxxyyz_xxxz_0, \
                             ta_xxxyyz_xxyy_0, \
                             ta_xxxyyz_xxyz_0, \
                             ta_xxxyyz_xxzz_0, \
                             ta_xxxyyz_xyyy_0, \
                             ta_xxxyyz_xyyz_0, \
                             ta_xxxyyz_xyzz_0, \
                             ta_xxxyyz_xzzz_0, \
                             ta_xxxyyz_yyyy_0, \
                             ta_xxxyyz_yyyz_0, \
                             ta_xxxyyz_yyzz_0, \
                             ta_xxxyyz_yzzz_0, \
                             ta_xxxyyz_zzzz_0, \
                             ta_xxxyz_xxxz_0,  \
                             ta_xxxyz_xxxz_1,  \
                             ta_xxxyz_xxzz_0,  \
                             ta_xxxyz_xxzz_1,  \
                             ta_xxxyz_xzzz_0,  \
                             ta_xxxyz_xzzz_1,  \
                             ta_xxxz_xxxz_0,   \
                             ta_xxxz_xxxz_1,   \
                             ta_xxxz_xxzz_0,   \
                             ta_xxxz_xxzz_1,   \
                             ta_xxxz_xzzz_0,   \
                             ta_xxxz_xzzz_1,   \
                             ta_xxyyz_yyyz_0,  \
                             ta_xxyyz_yyyz_1,  \
                             ta_xxyyz_yyzz_0,  \
                             ta_xxyyz_yyzz_1,  \
                             ta_xxyyz_yzzz_0,  \
                             ta_xxyyz_yzzz_1,  \
                             ta_xxyyz_zzzz_0,  \
                             ta_xxyyz_zzzz_1,  \
                             ta_xyyz_yyyz_0,   \
                             ta_xyyz_yyyz_1,   \
                             ta_xyyz_yyzz_0,   \
                             ta_xyyz_yyzz_1,   \
                             ta_xyyz_yzzz_0,   \
                             ta_xyyz_yzzz_1,   \
                             ta_xyyz_zzzz_0,   \
                             ta_xyyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyz_xxxx_0[i] = ta_xxxyy_xxxx_0[i] * pa_z[i] - ta_xxxyy_xxxx_1[i] * pc_z[i];

        ta_xxxyyz_xxxy_0[i] = ta_xxxyy_xxxy_0[i] * pa_z[i] - ta_xxxyy_xxxy_1[i] * pc_z[i];

        ta_xxxyyz_xxxz_0[i] = ta_xxxz_xxxz_0[i] * fe_0 - ta_xxxz_xxxz_1[i] * fe_0 + ta_xxxyz_xxxz_0[i] * pa_y[i] - ta_xxxyz_xxxz_1[i] * pc_y[i];

        ta_xxxyyz_xxyy_0[i] = ta_xxxyy_xxyy_0[i] * pa_z[i] - ta_xxxyy_xxyy_1[i] * pc_z[i];

        ta_xxxyyz_xxyz_0[i] = ta_xxxyy_xxy_0[i] * fe_0 - ta_xxxyy_xxy_1[i] * fe_0 + ta_xxxyy_xxyz_0[i] * pa_z[i] - ta_xxxyy_xxyz_1[i] * pc_z[i];

        ta_xxxyyz_xxzz_0[i] = ta_xxxz_xxzz_0[i] * fe_0 - ta_xxxz_xxzz_1[i] * fe_0 + ta_xxxyz_xxzz_0[i] * pa_y[i] - ta_xxxyz_xxzz_1[i] * pc_y[i];

        ta_xxxyyz_xyyy_0[i] = ta_xxxyy_xyyy_0[i] * pa_z[i] - ta_xxxyy_xyyy_1[i] * pc_z[i];

        ta_xxxyyz_xyyz_0[i] = ta_xxxyy_xyy_0[i] * fe_0 - ta_xxxyy_xyy_1[i] * fe_0 + ta_xxxyy_xyyz_0[i] * pa_z[i] - ta_xxxyy_xyyz_1[i] * pc_z[i];

        ta_xxxyyz_xyzz_0[i] =
            2.0 * ta_xxxyy_xyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xyz_1[i] * fe_0 + ta_xxxyy_xyzz_0[i] * pa_z[i] - ta_xxxyy_xyzz_1[i] * pc_z[i];

        ta_xxxyyz_xzzz_0[i] = ta_xxxz_xzzz_0[i] * fe_0 - ta_xxxz_xzzz_1[i] * fe_0 + ta_xxxyz_xzzz_0[i] * pa_y[i] - ta_xxxyz_xzzz_1[i] * pc_y[i];

        ta_xxxyyz_yyyy_0[i] = ta_xxxyy_yyyy_0[i] * pa_z[i] - ta_xxxyy_yyyy_1[i] * pc_z[i];

        ta_xxxyyz_yyyz_0[i] =
            2.0 * ta_xyyz_yyyz_0[i] * fe_0 - 2.0 * ta_xyyz_yyyz_1[i] * fe_0 + ta_xxyyz_yyyz_0[i] * pa_x[i] - ta_xxyyz_yyyz_1[i] * pc_x[i];

        ta_xxxyyz_yyzz_0[i] =
            2.0 * ta_xyyz_yyzz_0[i] * fe_0 - 2.0 * ta_xyyz_yyzz_1[i] * fe_0 + ta_xxyyz_yyzz_0[i] * pa_x[i] - ta_xxyyz_yyzz_1[i] * pc_x[i];

        ta_xxxyyz_yzzz_0[i] =
            2.0 * ta_xyyz_yzzz_0[i] * fe_0 - 2.0 * ta_xyyz_yzzz_1[i] * fe_0 + ta_xxyyz_yzzz_0[i] * pa_x[i] - ta_xxyyz_yzzz_1[i] * pc_x[i];

        ta_xxxyyz_zzzz_0[i] =
            2.0 * ta_xyyz_zzzz_0[i] * fe_0 - 2.0 * ta_xyyz_zzzz_1[i] * fe_0 + ta_xxyyz_zzzz_0[i] * pa_x[i] - ta_xxyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : IG

    auto ta_xxxyzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 120);

    auto ta_xxxyzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 121);

    auto ta_xxxyzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 122);

    auto ta_xxxyzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 123);

    auto ta_xxxyzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 124);

    auto ta_xxxyzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 125);

    auto ta_xxxyzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 126);

    auto ta_xxxyzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 127);

    auto ta_xxxyzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 128);

    auto ta_xxxyzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 129);

    auto ta_xxxyzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 130);

    auto ta_xxxyzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 131);

    auto ta_xxxyzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 132);

    auto ta_xxxyzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 133);

    auto ta_xxxyzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 134);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxxyzz_xxxx_0, \
                             ta_xxxyzz_xxxy_0, \
                             ta_xxxyzz_xxxz_0, \
                             ta_xxxyzz_xxyy_0, \
                             ta_xxxyzz_xxyz_0, \
                             ta_xxxyzz_xxzz_0, \
                             ta_xxxyzz_xyyy_0, \
                             ta_xxxyzz_xyyz_0, \
                             ta_xxxyzz_xyzz_0, \
                             ta_xxxyzz_xzzz_0, \
                             ta_xxxyzz_yyyy_0, \
                             ta_xxxyzz_yyyz_0, \
                             ta_xxxyzz_yyzz_0, \
                             ta_xxxyzz_yzzz_0, \
                             ta_xxxyzz_zzzz_0, \
                             ta_xxxzz_xxx_0,   \
                             ta_xxxzz_xxx_1,   \
                             ta_xxxzz_xxxx_0,  \
                             ta_xxxzz_xxxx_1,  \
                             ta_xxxzz_xxxy_0,  \
                             ta_xxxzz_xxxy_1,  \
                             ta_xxxzz_xxxz_0,  \
                             ta_xxxzz_xxxz_1,  \
                             ta_xxxzz_xxy_0,   \
                             ta_xxxzz_xxy_1,   \
                             ta_xxxzz_xxyy_0,  \
                             ta_xxxzz_xxyy_1,  \
                             ta_xxxzz_xxyz_0,  \
                             ta_xxxzz_xxyz_1,  \
                             ta_xxxzz_xxz_0,   \
                             ta_xxxzz_xxz_1,   \
                             ta_xxxzz_xxzz_0,  \
                             ta_xxxzz_xxzz_1,  \
                             ta_xxxzz_xyy_0,   \
                             ta_xxxzz_xyy_1,   \
                             ta_xxxzz_xyyy_0,  \
                             ta_xxxzz_xyyy_1,  \
                             ta_xxxzz_xyyz_0,  \
                             ta_xxxzz_xyyz_1,  \
                             ta_xxxzz_xyz_0,   \
                             ta_xxxzz_xyz_1,   \
                             ta_xxxzz_xyzz_0,  \
                             ta_xxxzz_xyzz_1,  \
                             ta_xxxzz_xzz_0,   \
                             ta_xxxzz_xzz_1,   \
                             ta_xxxzz_xzzz_0,  \
                             ta_xxxzz_xzzz_1,  \
                             ta_xxxzz_zzzz_0,  \
                             ta_xxxzz_zzzz_1,  \
                             ta_xxyzz_yyyy_0,  \
                             ta_xxyzz_yyyy_1,  \
                             ta_xxyzz_yyyz_0,  \
                             ta_xxyzz_yyyz_1,  \
                             ta_xxyzz_yyzz_0,  \
                             ta_xxyzz_yyzz_1,  \
                             ta_xxyzz_yzzz_0,  \
                             ta_xxyzz_yzzz_1,  \
                             ta_xyzz_yyyy_0,   \
                             ta_xyzz_yyyy_1,   \
                             ta_xyzz_yyyz_0,   \
                             ta_xyzz_yyyz_1,   \
                             ta_xyzz_yyzz_0,   \
                             ta_xyzz_yyzz_1,   \
                             ta_xyzz_yzzz_0,   \
                             ta_xyzz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyzz_xxxx_0[i] = ta_xxxzz_xxxx_0[i] * pa_y[i] - ta_xxxzz_xxxx_1[i] * pc_y[i];

        ta_xxxyzz_xxxy_0[i] = ta_xxxzz_xxx_0[i] * fe_0 - ta_xxxzz_xxx_1[i] * fe_0 + ta_xxxzz_xxxy_0[i] * pa_y[i] - ta_xxxzz_xxxy_1[i] * pc_y[i];

        ta_xxxyzz_xxxz_0[i] = ta_xxxzz_xxxz_0[i] * pa_y[i] - ta_xxxzz_xxxz_1[i] * pc_y[i];

        ta_xxxyzz_xxyy_0[i] =
            2.0 * ta_xxxzz_xxy_0[i] * fe_0 - 2.0 * ta_xxxzz_xxy_1[i] * fe_0 + ta_xxxzz_xxyy_0[i] * pa_y[i] - ta_xxxzz_xxyy_1[i] * pc_y[i];

        ta_xxxyzz_xxyz_0[i] = ta_xxxzz_xxz_0[i] * fe_0 - ta_xxxzz_xxz_1[i] * fe_0 + ta_xxxzz_xxyz_0[i] * pa_y[i] - ta_xxxzz_xxyz_1[i] * pc_y[i];

        ta_xxxyzz_xxzz_0[i] = ta_xxxzz_xxzz_0[i] * pa_y[i] - ta_xxxzz_xxzz_1[i] * pc_y[i];

        ta_xxxyzz_xyyy_0[i] =
            3.0 * ta_xxxzz_xyy_0[i] * fe_0 - 3.0 * ta_xxxzz_xyy_1[i] * fe_0 + ta_xxxzz_xyyy_0[i] * pa_y[i] - ta_xxxzz_xyyy_1[i] * pc_y[i];

        ta_xxxyzz_xyyz_0[i] =
            2.0 * ta_xxxzz_xyz_0[i] * fe_0 - 2.0 * ta_xxxzz_xyz_1[i] * fe_0 + ta_xxxzz_xyyz_0[i] * pa_y[i] - ta_xxxzz_xyyz_1[i] * pc_y[i];

        ta_xxxyzz_xyzz_0[i] = ta_xxxzz_xzz_0[i] * fe_0 - ta_xxxzz_xzz_1[i] * fe_0 + ta_xxxzz_xyzz_0[i] * pa_y[i] - ta_xxxzz_xyzz_1[i] * pc_y[i];

        ta_xxxyzz_xzzz_0[i] = ta_xxxzz_xzzz_0[i] * pa_y[i] - ta_xxxzz_xzzz_1[i] * pc_y[i];

        ta_xxxyzz_yyyy_0[i] =
            2.0 * ta_xyzz_yyyy_0[i] * fe_0 - 2.0 * ta_xyzz_yyyy_1[i] * fe_0 + ta_xxyzz_yyyy_0[i] * pa_x[i] - ta_xxyzz_yyyy_1[i] * pc_x[i];

        ta_xxxyzz_yyyz_0[i] =
            2.0 * ta_xyzz_yyyz_0[i] * fe_0 - 2.0 * ta_xyzz_yyyz_1[i] * fe_0 + ta_xxyzz_yyyz_0[i] * pa_x[i] - ta_xxyzz_yyyz_1[i] * pc_x[i];

        ta_xxxyzz_yyzz_0[i] =
            2.0 * ta_xyzz_yyzz_0[i] * fe_0 - 2.0 * ta_xyzz_yyzz_1[i] * fe_0 + ta_xxyzz_yyzz_0[i] * pa_x[i] - ta_xxyzz_yyzz_1[i] * pc_x[i];

        ta_xxxyzz_yzzz_0[i] =
            2.0 * ta_xyzz_yzzz_0[i] * fe_0 - 2.0 * ta_xyzz_yzzz_1[i] * fe_0 + ta_xxyzz_yzzz_0[i] * pa_x[i] - ta_xxyzz_yzzz_1[i] * pc_x[i];

        ta_xxxyzz_zzzz_0[i] = ta_xxxzz_zzzz_0[i] * pa_y[i] - ta_xxxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : IG

    auto ta_xxxzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 135);

    auto ta_xxxzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 136);

    auto ta_xxxzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 137);

    auto ta_xxxzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 138);

    auto ta_xxxzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 139);

    auto ta_xxxzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 140);

    auto ta_xxxzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 141);

    auto ta_xxxzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 142);

    auto ta_xxxzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 143);

    auto ta_xxxzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 144);

    auto ta_xxxzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 145);

    auto ta_xxxzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 146);

    auto ta_xxxzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 147);

    auto ta_xxxzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 148);

    auto ta_xxxzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 149);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxxz_xxxx_0,   \
                             ta_xxxz_xxxx_1,   \
                             ta_xxxz_xxxy_0,   \
                             ta_xxxz_xxxy_1,   \
                             ta_xxxz_xxyy_0,   \
                             ta_xxxz_xxyy_1,   \
                             ta_xxxz_xyyy_0,   \
                             ta_xxxz_xyyy_1,   \
                             ta_xxxzz_xxxx_0,  \
                             ta_xxxzz_xxxx_1,  \
                             ta_xxxzz_xxxy_0,  \
                             ta_xxxzz_xxxy_1,  \
                             ta_xxxzz_xxyy_0,  \
                             ta_xxxzz_xxyy_1,  \
                             ta_xxxzz_xyyy_0,  \
                             ta_xxxzz_xyyy_1,  \
                             ta_xxxzzz_xxxx_0, \
                             ta_xxxzzz_xxxy_0, \
                             ta_xxxzzz_xxxz_0, \
                             ta_xxxzzz_xxyy_0, \
                             ta_xxxzzz_xxyz_0, \
                             ta_xxxzzz_xxzz_0, \
                             ta_xxxzzz_xyyy_0, \
                             ta_xxxzzz_xyyz_0, \
                             ta_xxxzzz_xyzz_0, \
                             ta_xxxzzz_xzzz_0, \
                             ta_xxxzzz_yyyy_0, \
                             ta_xxxzzz_yyyz_0, \
                             ta_xxxzzz_yyzz_0, \
                             ta_xxxzzz_yzzz_0, \
                             ta_xxxzzz_zzzz_0, \
                             ta_xxzzz_xxxz_0,  \
                             ta_xxzzz_xxxz_1,  \
                             ta_xxzzz_xxyz_0,  \
                             ta_xxzzz_xxyz_1,  \
                             ta_xxzzz_xxz_0,   \
                             ta_xxzzz_xxz_1,   \
                             ta_xxzzz_xxzz_0,  \
                             ta_xxzzz_xxzz_1,  \
                             ta_xxzzz_xyyz_0,  \
                             ta_xxzzz_xyyz_1,  \
                             ta_xxzzz_xyz_0,   \
                             ta_xxzzz_xyz_1,   \
                             ta_xxzzz_xyzz_0,  \
                             ta_xxzzz_xyzz_1,  \
                             ta_xxzzz_xzz_0,   \
                             ta_xxzzz_xzz_1,   \
                             ta_xxzzz_xzzz_0,  \
                             ta_xxzzz_xzzz_1,  \
                             ta_xxzzz_yyyy_0,  \
                             ta_xxzzz_yyyy_1,  \
                             ta_xxzzz_yyyz_0,  \
                             ta_xxzzz_yyyz_1,  \
                             ta_xxzzz_yyz_0,   \
                             ta_xxzzz_yyz_1,   \
                             ta_xxzzz_yyzz_0,  \
                             ta_xxzzz_yyzz_1,  \
                             ta_xxzzz_yzz_0,   \
                             ta_xxzzz_yzz_1,   \
                             ta_xxzzz_yzzz_0,  \
                             ta_xxzzz_yzzz_1,  \
                             ta_xxzzz_zzz_0,   \
                             ta_xxzzz_zzz_1,   \
                             ta_xxzzz_zzzz_0,  \
                             ta_xxzzz_zzzz_1,  \
                             ta_xzzz_xxxz_0,   \
                             ta_xzzz_xxxz_1,   \
                             ta_xzzz_xxyz_0,   \
                             ta_xzzz_xxyz_1,   \
                             ta_xzzz_xxzz_0,   \
                             ta_xzzz_xxzz_1,   \
                             ta_xzzz_xyyz_0,   \
                             ta_xzzz_xyyz_1,   \
                             ta_xzzz_xyzz_0,   \
                             ta_xzzz_xyzz_1,   \
                             ta_xzzz_xzzz_0,   \
                             ta_xzzz_xzzz_1,   \
                             ta_xzzz_yyyy_0,   \
                             ta_xzzz_yyyy_1,   \
                             ta_xzzz_yyyz_0,   \
                             ta_xzzz_yyyz_1,   \
                             ta_xzzz_yyzz_0,   \
                             ta_xzzz_yyzz_1,   \
                             ta_xzzz_yzzz_0,   \
                             ta_xzzz_yzzz_1,   \
                             ta_xzzz_zzzz_0,   \
                             ta_xzzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzzz_xxxx_0[i] =
            2.0 * ta_xxxz_xxxx_0[i] * fe_0 - 2.0 * ta_xxxz_xxxx_1[i] * fe_0 + ta_xxxzz_xxxx_0[i] * pa_z[i] - ta_xxxzz_xxxx_1[i] * pc_z[i];

        ta_xxxzzz_xxxy_0[i] =
            2.0 * ta_xxxz_xxxy_0[i] * fe_0 - 2.0 * ta_xxxz_xxxy_1[i] * fe_0 + ta_xxxzz_xxxy_0[i] * pa_z[i] - ta_xxxzz_xxxy_1[i] * pc_z[i];

        ta_xxxzzz_xxxz_0[i] = 2.0 * ta_xzzz_xxxz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxz_1[i] * fe_0 + 3.0 * ta_xxzzz_xxz_0[i] * fe_0 -
                              3.0 * ta_xxzzz_xxz_1[i] * fe_0 + ta_xxzzz_xxxz_0[i] * pa_x[i] - ta_xxzzz_xxxz_1[i] * pc_x[i];

        ta_xxxzzz_xxyy_0[i] =
            2.0 * ta_xxxz_xxyy_0[i] * fe_0 - 2.0 * ta_xxxz_xxyy_1[i] * fe_0 + ta_xxxzz_xxyy_0[i] * pa_z[i] - ta_xxxzz_xxyy_1[i] * pc_z[i];

        ta_xxxzzz_xxyz_0[i] = 2.0 * ta_xzzz_xxyz_0[i] * fe_0 - 2.0 * ta_xzzz_xxyz_1[i] * fe_0 + 2.0 * ta_xxzzz_xyz_0[i] * fe_0 -
                              2.0 * ta_xxzzz_xyz_1[i] * fe_0 + ta_xxzzz_xxyz_0[i] * pa_x[i] - ta_xxzzz_xxyz_1[i] * pc_x[i];

        ta_xxxzzz_xxzz_0[i] = 2.0 * ta_xzzz_xxzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxzz_1[i] * fe_0 + 2.0 * ta_xxzzz_xzz_0[i] * fe_0 -
                              2.0 * ta_xxzzz_xzz_1[i] * fe_0 + ta_xxzzz_xxzz_0[i] * pa_x[i] - ta_xxzzz_xxzz_1[i] * pc_x[i];

        ta_xxxzzz_xyyy_0[i] =
            2.0 * ta_xxxz_xyyy_0[i] * fe_0 - 2.0 * ta_xxxz_xyyy_1[i] * fe_0 + ta_xxxzz_xyyy_0[i] * pa_z[i] - ta_xxxzz_xyyy_1[i] * pc_z[i];

        ta_xxxzzz_xyyz_0[i] = 2.0 * ta_xzzz_xyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyz_1[i] * fe_0 + ta_xxzzz_yyz_0[i] * fe_0 - ta_xxzzz_yyz_1[i] * fe_0 +
                              ta_xxzzz_xyyz_0[i] * pa_x[i] - ta_xxzzz_xyyz_1[i] * pc_x[i];

        ta_xxxzzz_xyzz_0[i] = 2.0 * ta_xzzz_xyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyzz_1[i] * fe_0 + ta_xxzzz_yzz_0[i] * fe_0 - ta_xxzzz_yzz_1[i] * fe_0 +
                              ta_xxzzz_xyzz_0[i] * pa_x[i] - ta_xxzzz_xyzz_1[i] * pc_x[i];

        ta_xxxzzz_xzzz_0[i] = 2.0 * ta_xzzz_xzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xzzz_1[i] * fe_0 + ta_xxzzz_zzz_0[i] * fe_0 - ta_xxzzz_zzz_1[i] * fe_0 +
                              ta_xxzzz_xzzz_0[i] * pa_x[i] - ta_xxzzz_xzzz_1[i] * pc_x[i];

        ta_xxxzzz_yyyy_0[i] =
            2.0 * ta_xzzz_yyyy_0[i] * fe_0 - 2.0 * ta_xzzz_yyyy_1[i] * fe_0 + ta_xxzzz_yyyy_0[i] * pa_x[i] - ta_xxzzz_yyyy_1[i] * pc_x[i];

        ta_xxxzzz_yyyz_0[i] =
            2.0 * ta_xzzz_yyyz_0[i] * fe_0 - 2.0 * ta_xzzz_yyyz_1[i] * fe_0 + ta_xxzzz_yyyz_0[i] * pa_x[i] - ta_xxzzz_yyyz_1[i] * pc_x[i];

        ta_xxxzzz_yyzz_0[i] =
            2.0 * ta_xzzz_yyzz_0[i] * fe_0 - 2.0 * ta_xzzz_yyzz_1[i] * fe_0 + ta_xxzzz_yyzz_0[i] * pa_x[i] - ta_xxzzz_yyzz_1[i] * pc_x[i];

        ta_xxxzzz_yzzz_0[i] =
            2.0 * ta_xzzz_yzzz_0[i] * fe_0 - 2.0 * ta_xzzz_yzzz_1[i] * fe_0 + ta_xxzzz_yzzz_0[i] * pa_x[i] - ta_xxzzz_yzzz_1[i] * pc_x[i];

        ta_xxxzzz_zzzz_0[i] =
            2.0 * ta_xzzz_zzzz_0[i] * fe_0 - 2.0 * ta_xzzz_zzzz_1[i] * fe_0 + ta_xxzzz_zzzz_0[i] * pa_x[i] - ta_xxzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : IG

    auto ta_xxyyyy_xxxx_0 = pbuffer.data(idx_npot_0_ig + 150);

    auto ta_xxyyyy_xxxy_0 = pbuffer.data(idx_npot_0_ig + 151);

    auto ta_xxyyyy_xxxz_0 = pbuffer.data(idx_npot_0_ig + 152);

    auto ta_xxyyyy_xxyy_0 = pbuffer.data(idx_npot_0_ig + 153);

    auto ta_xxyyyy_xxyz_0 = pbuffer.data(idx_npot_0_ig + 154);

    auto ta_xxyyyy_xxzz_0 = pbuffer.data(idx_npot_0_ig + 155);

    auto ta_xxyyyy_xyyy_0 = pbuffer.data(idx_npot_0_ig + 156);

    auto ta_xxyyyy_xyyz_0 = pbuffer.data(idx_npot_0_ig + 157);

    auto ta_xxyyyy_xyzz_0 = pbuffer.data(idx_npot_0_ig + 158);

    auto ta_xxyyyy_xzzz_0 = pbuffer.data(idx_npot_0_ig + 159);

    auto ta_xxyyyy_yyyy_0 = pbuffer.data(idx_npot_0_ig + 160);

    auto ta_xxyyyy_yyyz_0 = pbuffer.data(idx_npot_0_ig + 161);

    auto ta_xxyyyy_yyzz_0 = pbuffer.data(idx_npot_0_ig + 162);

    auto ta_xxyyyy_yzzz_0 = pbuffer.data(idx_npot_0_ig + 163);

    auto ta_xxyyyy_zzzz_0 = pbuffer.data(idx_npot_0_ig + 164);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxyy_xxxx_0,   \
                             ta_xxyy_xxxx_1,   \
                             ta_xxyy_xxxz_0,   \
                             ta_xxyy_xxxz_1,   \
                             ta_xxyy_xxzz_0,   \
                             ta_xxyy_xxzz_1,   \
                             ta_xxyy_xzzz_0,   \
                             ta_xxyy_xzzz_1,   \
                             ta_xxyyy_xxxx_0,  \
                             ta_xxyyy_xxxx_1,  \
                             ta_xxyyy_xxxz_0,  \
                             ta_xxyyy_xxxz_1,  \
                             ta_xxyyy_xxzz_0,  \
                             ta_xxyyy_xxzz_1,  \
                             ta_xxyyy_xzzz_0,  \
                             ta_xxyyy_xzzz_1,  \
                             ta_xxyyyy_xxxx_0, \
                             ta_xxyyyy_xxxy_0, \
                             ta_xxyyyy_xxxz_0, \
                             ta_xxyyyy_xxyy_0, \
                             ta_xxyyyy_xxyz_0, \
                             ta_xxyyyy_xxzz_0, \
                             ta_xxyyyy_xyyy_0, \
                             ta_xxyyyy_xyyz_0, \
                             ta_xxyyyy_xyzz_0, \
                             ta_xxyyyy_xzzz_0, \
                             ta_xxyyyy_yyyy_0, \
                             ta_xxyyyy_yyyz_0, \
                             ta_xxyyyy_yyzz_0, \
                             ta_xxyyyy_yzzz_0, \
                             ta_xxyyyy_zzzz_0, \
                             ta_xyyyy_xxxy_0,  \
                             ta_xyyyy_xxxy_1,  \
                             ta_xyyyy_xxy_0,   \
                             ta_xyyyy_xxy_1,   \
                             ta_xyyyy_xxyy_0,  \
                             ta_xyyyy_xxyy_1,  \
                             ta_xyyyy_xxyz_0,  \
                             ta_xyyyy_xxyz_1,  \
                             ta_xyyyy_xyy_0,   \
                             ta_xyyyy_xyy_1,   \
                             ta_xyyyy_xyyy_0,  \
                             ta_xyyyy_xyyy_1,  \
                             ta_xyyyy_xyyz_0,  \
                             ta_xyyyy_xyyz_1,  \
                             ta_xyyyy_xyz_0,   \
                             ta_xyyyy_xyz_1,   \
                             ta_xyyyy_xyzz_0,  \
                             ta_xyyyy_xyzz_1,  \
                             ta_xyyyy_yyy_0,   \
                             ta_xyyyy_yyy_1,   \
                             ta_xyyyy_yyyy_0,  \
                             ta_xyyyy_yyyy_1,  \
                             ta_xyyyy_yyyz_0,  \
                             ta_xyyyy_yyyz_1,  \
                             ta_xyyyy_yyz_0,   \
                             ta_xyyyy_yyz_1,   \
                             ta_xyyyy_yyzz_0,  \
                             ta_xyyyy_yyzz_1,  \
                             ta_xyyyy_yzz_0,   \
                             ta_xyyyy_yzz_1,   \
                             ta_xyyyy_yzzz_0,  \
                             ta_xyyyy_yzzz_1,  \
                             ta_xyyyy_zzzz_0,  \
                             ta_xyyyy_zzzz_1,  \
                             ta_yyyy_xxxy_0,   \
                             ta_yyyy_xxxy_1,   \
                             ta_yyyy_xxyy_0,   \
                             ta_yyyy_xxyy_1,   \
                             ta_yyyy_xxyz_0,   \
                             ta_yyyy_xxyz_1,   \
                             ta_yyyy_xyyy_0,   \
                             ta_yyyy_xyyy_1,   \
                             ta_yyyy_xyyz_0,   \
                             ta_yyyy_xyyz_1,   \
                             ta_yyyy_xyzz_0,   \
                             ta_yyyy_xyzz_1,   \
                             ta_yyyy_yyyy_0,   \
                             ta_yyyy_yyyy_1,   \
                             ta_yyyy_yyyz_0,   \
                             ta_yyyy_yyyz_1,   \
                             ta_yyyy_yyzz_0,   \
                             ta_yyyy_yyzz_1,   \
                             ta_yyyy_yzzz_0,   \
                             ta_yyyy_yzzz_1,   \
                             ta_yyyy_zzzz_0,   \
                             ta_yyyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyy_xxxx_0[i] =
            3.0 * ta_xxyy_xxxx_0[i] * fe_0 - 3.0 * ta_xxyy_xxxx_1[i] * fe_0 + ta_xxyyy_xxxx_0[i] * pa_y[i] - ta_xxyyy_xxxx_1[i] * pc_y[i];

        ta_xxyyyy_xxxy_0[i] = ta_yyyy_xxxy_0[i] * fe_0 - ta_yyyy_xxxy_1[i] * fe_0 + 3.0 * ta_xyyyy_xxy_0[i] * fe_0 - 3.0 * ta_xyyyy_xxy_1[i] * fe_0 +
                              ta_xyyyy_xxxy_0[i] * pa_x[i] - ta_xyyyy_xxxy_1[i] * pc_x[i];

        ta_xxyyyy_xxxz_0[i] =
            3.0 * ta_xxyy_xxxz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxz_1[i] * fe_0 + ta_xxyyy_xxxz_0[i] * pa_y[i] - ta_xxyyy_xxxz_1[i] * pc_y[i];

        ta_xxyyyy_xxyy_0[i] = ta_yyyy_xxyy_0[i] * fe_0 - ta_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta_xyyyy_xyy_0[i] * fe_0 - 2.0 * ta_xyyyy_xyy_1[i] * fe_0 +
                              ta_xyyyy_xxyy_0[i] * pa_x[i] - ta_xyyyy_xxyy_1[i] * pc_x[i];

        ta_xxyyyy_xxyz_0[i] = ta_yyyy_xxyz_0[i] * fe_0 - ta_yyyy_xxyz_1[i] * fe_0 + 2.0 * ta_xyyyy_xyz_0[i] * fe_0 - 2.0 * ta_xyyyy_xyz_1[i] * fe_0 +
                              ta_xyyyy_xxyz_0[i] * pa_x[i] - ta_xyyyy_xxyz_1[i] * pc_x[i];

        ta_xxyyyy_xxzz_0[i] =
            3.0 * ta_xxyy_xxzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxzz_1[i] * fe_0 + ta_xxyyy_xxzz_0[i] * pa_y[i] - ta_xxyyy_xxzz_1[i] * pc_y[i];

        ta_xxyyyy_xyyy_0[i] = ta_yyyy_xyyy_0[i] * fe_0 - ta_yyyy_xyyy_1[i] * fe_0 + ta_xyyyy_yyy_0[i] * fe_0 - ta_xyyyy_yyy_1[i] * fe_0 +
                              ta_xyyyy_xyyy_0[i] * pa_x[i] - ta_xyyyy_xyyy_1[i] * pc_x[i];

        ta_xxyyyy_xyyz_0[i] = ta_yyyy_xyyz_0[i] * fe_0 - ta_yyyy_xyyz_1[i] * fe_0 + ta_xyyyy_yyz_0[i] * fe_0 - ta_xyyyy_yyz_1[i] * fe_0 +
                              ta_xyyyy_xyyz_0[i] * pa_x[i] - ta_xyyyy_xyyz_1[i] * pc_x[i];

        ta_xxyyyy_xyzz_0[i] = ta_yyyy_xyzz_0[i] * fe_0 - ta_yyyy_xyzz_1[i] * fe_0 + ta_xyyyy_yzz_0[i] * fe_0 - ta_xyyyy_yzz_1[i] * fe_0 +
                              ta_xyyyy_xyzz_0[i] * pa_x[i] - ta_xyyyy_xyzz_1[i] * pc_x[i];

        ta_xxyyyy_xzzz_0[i] =
            3.0 * ta_xxyy_xzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xzzz_1[i] * fe_0 + ta_xxyyy_xzzz_0[i] * pa_y[i] - ta_xxyyy_xzzz_1[i] * pc_y[i];

        ta_xxyyyy_yyyy_0[i] = ta_yyyy_yyyy_0[i] * fe_0 - ta_yyyy_yyyy_1[i] * fe_0 + ta_xyyyy_yyyy_0[i] * pa_x[i] - ta_xyyyy_yyyy_1[i] * pc_x[i];

        ta_xxyyyy_yyyz_0[i] = ta_yyyy_yyyz_0[i] * fe_0 - ta_yyyy_yyyz_1[i] * fe_0 + ta_xyyyy_yyyz_0[i] * pa_x[i] - ta_xyyyy_yyyz_1[i] * pc_x[i];

        ta_xxyyyy_yyzz_0[i] = ta_yyyy_yyzz_0[i] * fe_0 - ta_yyyy_yyzz_1[i] * fe_0 + ta_xyyyy_yyzz_0[i] * pa_x[i] - ta_xyyyy_yyzz_1[i] * pc_x[i];

        ta_xxyyyy_yzzz_0[i] = ta_yyyy_yzzz_0[i] * fe_0 - ta_yyyy_yzzz_1[i] * fe_0 + ta_xyyyy_yzzz_0[i] * pa_x[i] - ta_xyyyy_yzzz_1[i] * pc_x[i];

        ta_xxyyyy_zzzz_0[i] = ta_yyyy_zzzz_0[i] * fe_0 - ta_yyyy_zzzz_1[i] * fe_0 + ta_xyyyy_zzzz_0[i] * pa_x[i] - ta_xyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 165-180 components of targeted buffer : IG

    auto ta_xxyyyz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 165);

    auto ta_xxyyyz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 166);

    auto ta_xxyyyz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 167);

    auto ta_xxyyyz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 168);

    auto ta_xxyyyz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 169);

    auto ta_xxyyyz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 170);

    auto ta_xxyyyz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 171);

    auto ta_xxyyyz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 172);

    auto ta_xxyyyz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 173);

    auto ta_xxyyyz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 174);

    auto ta_xxyyyz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 175);

    auto ta_xxyyyz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 176);

    auto ta_xxyyyz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 177);

    auto ta_xxyyyz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 178);

    auto ta_xxyyyz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 179);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta_xxyyy_xxxx_0,  \
                             ta_xxyyy_xxxx_1,  \
                             ta_xxyyy_xxxy_0,  \
                             ta_xxyyy_xxxy_1,  \
                             ta_xxyyy_xxy_0,   \
                             ta_xxyyy_xxy_1,   \
                             ta_xxyyy_xxyy_0,  \
                             ta_xxyyy_xxyy_1,  \
                             ta_xxyyy_xxyz_0,  \
                             ta_xxyyy_xxyz_1,  \
                             ta_xxyyy_xyy_0,   \
                             ta_xxyyy_xyy_1,   \
                             ta_xxyyy_xyyy_0,  \
                             ta_xxyyy_xyyy_1,  \
                             ta_xxyyy_xyyz_0,  \
                             ta_xxyyy_xyyz_1,  \
                             ta_xxyyy_xyz_0,   \
                             ta_xxyyy_xyz_1,   \
                             ta_xxyyy_xyzz_0,  \
                             ta_xxyyy_xyzz_1,  \
                             ta_xxyyy_yyyy_0,  \
                             ta_xxyyy_yyyy_1,  \
                             ta_xxyyyz_xxxx_0, \
                             ta_xxyyyz_xxxy_0, \
                             ta_xxyyyz_xxxz_0, \
                             ta_xxyyyz_xxyy_0, \
                             ta_xxyyyz_xxyz_0, \
                             ta_xxyyyz_xxzz_0, \
                             ta_xxyyyz_xyyy_0, \
                             ta_xxyyyz_xyyz_0, \
                             ta_xxyyyz_xyzz_0, \
                             ta_xxyyyz_xzzz_0, \
                             ta_xxyyyz_yyyy_0, \
                             ta_xxyyyz_yyyz_0, \
                             ta_xxyyyz_yyzz_0, \
                             ta_xxyyyz_yzzz_0, \
                             ta_xxyyyz_zzzz_0, \
                             ta_xxyyz_xxxz_0,  \
                             ta_xxyyz_xxxz_1,  \
                             ta_xxyyz_xxzz_0,  \
                             ta_xxyyz_xxzz_1,  \
                             ta_xxyyz_xzzz_0,  \
                             ta_xxyyz_xzzz_1,  \
                             ta_xxyz_xxxz_0,   \
                             ta_xxyz_xxxz_1,   \
                             ta_xxyz_xxzz_0,   \
                             ta_xxyz_xxzz_1,   \
                             ta_xxyz_xzzz_0,   \
                             ta_xxyz_xzzz_1,   \
                             ta_xyyyz_yyyz_0,  \
                             ta_xyyyz_yyyz_1,  \
                             ta_xyyyz_yyzz_0,  \
                             ta_xyyyz_yyzz_1,  \
                             ta_xyyyz_yzzz_0,  \
                             ta_xyyyz_yzzz_1,  \
                             ta_xyyyz_zzzz_0,  \
                             ta_xyyyz_zzzz_1,  \
                             ta_yyyz_yyyz_0,   \
                             ta_yyyz_yyyz_1,   \
                             ta_yyyz_yyzz_0,   \
                             ta_yyyz_yyzz_1,   \
                             ta_yyyz_yzzz_0,   \
                             ta_yyyz_yzzz_1,   \
                             ta_yyyz_zzzz_0,   \
                             ta_yyyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyz_xxxx_0[i] = ta_xxyyy_xxxx_0[i] * pa_z[i] - ta_xxyyy_xxxx_1[i] * pc_z[i];

        ta_xxyyyz_xxxy_0[i] = ta_xxyyy_xxxy_0[i] * pa_z[i] - ta_xxyyy_xxxy_1[i] * pc_z[i];

        ta_xxyyyz_xxxz_0[i] =
            2.0 * ta_xxyz_xxxz_0[i] * fe_0 - 2.0 * ta_xxyz_xxxz_1[i] * fe_0 + ta_xxyyz_xxxz_0[i] * pa_y[i] - ta_xxyyz_xxxz_1[i] * pc_y[i];

        ta_xxyyyz_xxyy_0[i] = ta_xxyyy_xxyy_0[i] * pa_z[i] - ta_xxyyy_xxyy_1[i] * pc_z[i];

        ta_xxyyyz_xxyz_0[i] = ta_xxyyy_xxy_0[i] * fe_0 - ta_xxyyy_xxy_1[i] * fe_0 + ta_xxyyy_xxyz_0[i] * pa_z[i] - ta_xxyyy_xxyz_1[i] * pc_z[i];

        ta_xxyyyz_xxzz_0[i] =
            2.0 * ta_xxyz_xxzz_0[i] * fe_0 - 2.0 * ta_xxyz_xxzz_1[i] * fe_0 + ta_xxyyz_xxzz_0[i] * pa_y[i] - ta_xxyyz_xxzz_1[i] * pc_y[i];

        ta_xxyyyz_xyyy_0[i] = ta_xxyyy_xyyy_0[i] * pa_z[i] - ta_xxyyy_xyyy_1[i] * pc_z[i];

        ta_xxyyyz_xyyz_0[i] = ta_xxyyy_xyy_0[i] * fe_0 - ta_xxyyy_xyy_1[i] * fe_0 + ta_xxyyy_xyyz_0[i] * pa_z[i] - ta_xxyyy_xyyz_1[i] * pc_z[i];

        ta_xxyyyz_xyzz_0[i] =
            2.0 * ta_xxyyy_xyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xyz_1[i] * fe_0 + ta_xxyyy_xyzz_0[i] * pa_z[i] - ta_xxyyy_xyzz_1[i] * pc_z[i];

        ta_xxyyyz_xzzz_0[i] =
            2.0 * ta_xxyz_xzzz_0[i] * fe_0 - 2.0 * ta_xxyz_xzzz_1[i] * fe_0 + ta_xxyyz_xzzz_0[i] * pa_y[i] - ta_xxyyz_xzzz_1[i] * pc_y[i];

        ta_xxyyyz_yyyy_0[i] = ta_xxyyy_yyyy_0[i] * pa_z[i] - ta_xxyyy_yyyy_1[i] * pc_z[i];

        ta_xxyyyz_yyyz_0[i] = ta_yyyz_yyyz_0[i] * fe_0 - ta_yyyz_yyyz_1[i] * fe_0 + ta_xyyyz_yyyz_0[i] * pa_x[i] - ta_xyyyz_yyyz_1[i] * pc_x[i];

        ta_xxyyyz_yyzz_0[i] = ta_yyyz_yyzz_0[i] * fe_0 - ta_yyyz_yyzz_1[i] * fe_0 + ta_xyyyz_yyzz_0[i] * pa_x[i] - ta_xyyyz_yyzz_1[i] * pc_x[i];

        ta_xxyyyz_yzzz_0[i] = ta_yyyz_yzzz_0[i] * fe_0 - ta_yyyz_yzzz_1[i] * fe_0 + ta_xyyyz_yzzz_0[i] * pa_x[i] - ta_xyyyz_yzzz_1[i] * pc_x[i];

        ta_xxyyyz_zzzz_0[i] = ta_yyyz_zzzz_0[i] * fe_0 - ta_yyyz_zzzz_1[i] * fe_0 + ta_xyyyz_zzzz_0[i] * pa_x[i] - ta_xyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 180-195 components of targeted buffer : IG

    auto ta_xxyyzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 180);

    auto ta_xxyyzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 181);

    auto ta_xxyyzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 182);

    auto ta_xxyyzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 183);

    auto ta_xxyyzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 184);

    auto ta_xxyyzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 185);

    auto ta_xxyyzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 186);

    auto ta_xxyyzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 187);

    auto ta_xxyyzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 188);

    auto ta_xxyyzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 189);

    auto ta_xxyyzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 190);

    auto ta_xxyyzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 191);

    auto ta_xxyyzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 192);

    auto ta_xxyyzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 193);

    auto ta_xxyyzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 194);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta_xxyy_xxxy_0,   \
                             ta_xxyy_xxxy_1,   \
                             ta_xxyy_xxyy_0,   \
                             ta_xxyy_xxyy_1,   \
                             ta_xxyy_xyyy_0,   \
                             ta_xxyy_xyyy_1,   \
                             ta_xxyyz_xxxy_0,  \
                             ta_xxyyz_xxxy_1,  \
                             ta_xxyyz_xxyy_0,  \
                             ta_xxyyz_xxyy_1,  \
                             ta_xxyyz_xyyy_0,  \
                             ta_xxyyz_xyyy_1,  \
                             ta_xxyyzz_xxxx_0, \
                             ta_xxyyzz_xxxy_0, \
                             ta_xxyyzz_xxxz_0, \
                             ta_xxyyzz_xxyy_0, \
                             ta_xxyyzz_xxyz_0, \
                             ta_xxyyzz_xxzz_0, \
                             ta_xxyyzz_xyyy_0, \
                             ta_xxyyzz_xyyz_0, \
                             ta_xxyyzz_xyzz_0, \
                             ta_xxyyzz_xzzz_0, \
                             ta_xxyyzz_yyyy_0, \
                             ta_xxyyzz_yyyz_0, \
                             ta_xxyyzz_yyzz_0, \
                             ta_xxyyzz_yzzz_0, \
                             ta_xxyyzz_zzzz_0, \
                             ta_xxyzz_xxxx_0,  \
                             ta_xxyzz_xxxx_1,  \
                             ta_xxyzz_xxxz_0,  \
                             ta_xxyzz_xxxz_1,  \
                             ta_xxyzz_xxzz_0,  \
                             ta_xxyzz_xxzz_1,  \
                             ta_xxyzz_xzzz_0,  \
                             ta_xxyzz_xzzz_1,  \
                             ta_xxzz_xxxx_0,   \
                             ta_xxzz_xxxx_1,   \
                             ta_xxzz_xxxz_0,   \
                             ta_xxzz_xxxz_1,   \
                             ta_xxzz_xxzz_0,   \
                             ta_xxzz_xxzz_1,   \
                             ta_xxzz_xzzz_0,   \
                             ta_xxzz_xzzz_1,   \
                             ta_xyyzz_xxyz_0,  \
                             ta_xyyzz_xxyz_1,  \
                             ta_xyyzz_xyyz_0,  \
                             ta_xyyzz_xyyz_1,  \
                             ta_xyyzz_xyz_0,   \
                             ta_xyyzz_xyz_1,   \
                             ta_xyyzz_xyzz_0,  \
                             ta_xyyzz_xyzz_1,  \
                             ta_xyyzz_yyyy_0,  \
                             ta_xyyzz_yyyy_1,  \
                             ta_xyyzz_yyyz_0,  \
                             ta_xyyzz_yyyz_1,  \
                             ta_xyyzz_yyz_0,   \
                             ta_xyyzz_yyz_1,   \
                             ta_xyyzz_yyzz_0,  \
                             ta_xyyzz_yyzz_1,  \
                             ta_xyyzz_yzz_0,   \
                             ta_xyyzz_yzz_1,   \
                             ta_xyyzz_yzzz_0,  \
                             ta_xyyzz_yzzz_1,  \
                             ta_xyyzz_zzzz_0,  \
                             ta_xyyzz_zzzz_1,  \
                             ta_yyzz_xxyz_0,   \
                             ta_yyzz_xxyz_1,   \
                             ta_yyzz_xyyz_0,   \
                             ta_yyzz_xyyz_1,   \
                             ta_yyzz_xyzz_0,   \
                             ta_yyzz_xyzz_1,   \
                             ta_yyzz_yyyy_0,   \
                             ta_yyzz_yyyy_1,   \
                             ta_yyzz_yyyz_0,   \
                             ta_yyzz_yyyz_1,   \
                             ta_yyzz_yyzz_0,   \
                             ta_yyzz_yyzz_1,   \
                             ta_yyzz_yzzz_0,   \
                             ta_yyzz_yzzz_1,   \
                             ta_yyzz_zzzz_0,   \
                             ta_yyzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyzz_xxxx_0[i] = ta_xxzz_xxxx_0[i] * fe_0 - ta_xxzz_xxxx_1[i] * fe_0 + ta_xxyzz_xxxx_0[i] * pa_y[i] - ta_xxyzz_xxxx_1[i] * pc_y[i];

        ta_xxyyzz_xxxy_0[i] = ta_xxyy_xxxy_0[i] * fe_0 - ta_xxyy_xxxy_1[i] * fe_0 + ta_xxyyz_xxxy_0[i] * pa_z[i] - ta_xxyyz_xxxy_1[i] * pc_z[i];

        ta_xxyyzz_xxxz_0[i] = ta_xxzz_xxxz_0[i] * fe_0 - ta_xxzz_xxxz_1[i] * fe_0 + ta_xxyzz_xxxz_0[i] * pa_y[i] - ta_xxyzz_xxxz_1[i] * pc_y[i];

        ta_xxyyzz_xxyy_0[i] = ta_xxyy_xxyy_0[i] * fe_0 - ta_xxyy_xxyy_1[i] * fe_0 + ta_xxyyz_xxyy_0[i] * pa_z[i] - ta_xxyyz_xxyy_1[i] * pc_z[i];

        ta_xxyyzz_xxyz_0[i] = ta_yyzz_xxyz_0[i] * fe_0 - ta_yyzz_xxyz_1[i] * fe_0 + 2.0 * ta_xyyzz_xyz_0[i] * fe_0 - 2.0 * ta_xyyzz_xyz_1[i] * fe_0 +
                              ta_xyyzz_xxyz_0[i] * pa_x[i] - ta_xyyzz_xxyz_1[i] * pc_x[i];

        ta_xxyyzz_xxzz_0[i] = ta_xxzz_xxzz_0[i] * fe_0 - ta_xxzz_xxzz_1[i] * fe_0 + ta_xxyzz_xxzz_0[i] * pa_y[i] - ta_xxyzz_xxzz_1[i] * pc_y[i];

        ta_xxyyzz_xyyy_0[i] = ta_xxyy_xyyy_0[i] * fe_0 - ta_xxyy_xyyy_1[i] * fe_0 + ta_xxyyz_xyyy_0[i] * pa_z[i] - ta_xxyyz_xyyy_1[i] * pc_z[i];

        ta_xxyyzz_xyyz_0[i] = ta_yyzz_xyyz_0[i] * fe_0 - ta_yyzz_xyyz_1[i] * fe_0 + ta_xyyzz_yyz_0[i] * fe_0 - ta_xyyzz_yyz_1[i] * fe_0 +
                              ta_xyyzz_xyyz_0[i] * pa_x[i] - ta_xyyzz_xyyz_1[i] * pc_x[i];

        ta_xxyyzz_xyzz_0[i] = ta_yyzz_xyzz_0[i] * fe_0 - ta_yyzz_xyzz_1[i] * fe_0 + ta_xyyzz_yzz_0[i] * fe_0 - ta_xyyzz_yzz_1[i] * fe_0 +
                              ta_xyyzz_xyzz_0[i] * pa_x[i] - ta_xyyzz_xyzz_1[i] * pc_x[i];

        ta_xxyyzz_xzzz_0[i] = ta_xxzz_xzzz_0[i] * fe_0 - ta_xxzz_xzzz_1[i] * fe_0 + ta_xxyzz_xzzz_0[i] * pa_y[i] - ta_xxyzz_xzzz_1[i] * pc_y[i];

        ta_xxyyzz_yyyy_0[i] = ta_yyzz_yyyy_0[i] * fe_0 - ta_yyzz_yyyy_1[i] * fe_0 + ta_xyyzz_yyyy_0[i] * pa_x[i] - ta_xyyzz_yyyy_1[i] * pc_x[i];

        ta_xxyyzz_yyyz_0[i] = ta_yyzz_yyyz_0[i] * fe_0 - ta_yyzz_yyyz_1[i] * fe_0 + ta_xyyzz_yyyz_0[i] * pa_x[i] - ta_xyyzz_yyyz_1[i] * pc_x[i];

        ta_xxyyzz_yyzz_0[i] = ta_yyzz_yyzz_0[i] * fe_0 - ta_yyzz_yyzz_1[i] * fe_0 + ta_xyyzz_yyzz_0[i] * pa_x[i] - ta_xyyzz_yyzz_1[i] * pc_x[i];

        ta_xxyyzz_yzzz_0[i] = ta_yyzz_yzzz_0[i] * fe_0 - ta_yyzz_yzzz_1[i] * fe_0 + ta_xyyzz_yzzz_0[i] * pa_x[i] - ta_xyyzz_yzzz_1[i] * pc_x[i];

        ta_xxyyzz_zzzz_0[i] = ta_yyzz_zzzz_0[i] * fe_0 - ta_yyzz_zzzz_1[i] * fe_0 + ta_xyyzz_zzzz_0[i] * pa_x[i] - ta_xyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : IG

    auto ta_xxyzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 195);

    auto ta_xxyzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 196);

    auto ta_xxyzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 197);

    auto ta_xxyzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 198);

    auto ta_xxyzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 199);

    auto ta_xxyzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 200);

    auto ta_xxyzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 201);

    auto ta_xxyzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 202);

    auto ta_xxyzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 203);

    auto ta_xxyzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 204);

    auto ta_xxyzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 205);

    auto ta_xxyzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 206);

    auto ta_xxyzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 207);

    auto ta_xxyzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 208);

    auto ta_xxyzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 209);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxyzzz_xxxx_0, \
                             ta_xxyzzz_xxxy_0, \
                             ta_xxyzzz_xxxz_0, \
                             ta_xxyzzz_xxyy_0, \
                             ta_xxyzzz_xxyz_0, \
                             ta_xxyzzz_xxzz_0, \
                             ta_xxyzzz_xyyy_0, \
                             ta_xxyzzz_xyyz_0, \
                             ta_xxyzzz_xyzz_0, \
                             ta_xxyzzz_xzzz_0, \
                             ta_xxyzzz_yyyy_0, \
                             ta_xxyzzz_yyyz_0, \
                             ta_xxyzzz_yyzz_0, \
                             ta_xxyzzz_yzzz_0, \
                             ta_xxyzzz_zzzz_0, \
                             ta_xxzzz_xxx_0,   \
                             ta_xxzzz_xxx_1,   \
                             ta_xxzzz_xxxx_0,  \
                             ta_xxzzz_xxxx_1,  \
                             ta_xxzzz_xxxy_0,  \
                             ta_xxzzz_xxxy_1,  \
                             ta_xxzzz_xxxz_0,  \
                             ta_xxzzz_xxxz_1,  \
                             ta_xxzzz_xxy_0,   \
                             ta_xxzzz_xxy_1,   \
                             ta_xxzzz_xxyy_0,  \
                             ta_xxzzz_xxyy_1,  \
                             ta_xxzzz_xxyz_0,  \
                             ta_xxzzz_xxyz_1,  \
                             ta_xxzzz_xxz_0,   \
                             ta_xxzzz_xxz_1,   \
                             ta_xxzzz_xxzz_0,  \
                             ta_xxzzz_xxzz_1,  \
                             ta_xxzzz_xyy_0,   \
                             ta_xxzzz_xyy_1,   \
                             ta_xxzzz_xyyy_0,  \
                             ta_xxzzz_xyyy_1,  \
                             ta_xxzzz_xyyz_0,  \
                             ta_xxzzz_xyyz_1,  \
                             ta_xxzzz_xyz_0,   \
                             ta_xxzzz_xyz_1,   \
                             ta_xxzzz_xyzz_0,  \
                             ta_xxzzz_xyzz_1,  \
                             ta_xxzzz_xzz_0,   \
                             ta_xxzzz_xzz_1,   \
                             ta_xxzzz_xzzz_0,  \
                             ta_xxzzz_xzzz_1,  \
                             ta_xxzzz_zzzz_0,  \
                             ta_xxzzz_zzzz_1,  \
                             ta_xyzzz_yyyy_0,  \
                             ta_xyzzz_yyyy_1,  \
                             ta_xyzzz_yyyz_0,  \
                             ta_xyzzz_yyyz_1,  \
                             ta_xyzzz_yyzz_0,  \
                             ta_xyzzz_yyzz_1,  \
                             ta_xyzzz_yzzz_0,  \
                             ta_xyzzz_yzzz_1,  \
                             ta_yzzz_yyyy_0,   \
                             ta_yzzz_yyyy_1,   \
                             ta_yzzz_yyyz_0,   \
                             ta_yzzz_yyyz_1,   \
                             ta_yzzz_yyzz_0,   \
                             ta_yzzz_yyzz_1,   \
                             ta_yzzz_yzzz_0,   \
                             ta_yzzz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzzz_xxxx_0[i] = ta_xxzzz_xxxx_0[i] * pa_y[i] - ta_xxzzz_xxxx_1[i] * pc_y[i];

        ta_xxyzzz_xxxy_0[i] = ta_xxzzz_xxx_0[i] * fe_0 - ta_xxzzz_xxx_1[i] * fe_0 + ta_xxzzz_xxxy_0[i] * pa_y[i] - ta_xxzzz_xxxy_1[i] * pc_y[i];

        ta_xxyzzz_xxxz_0[i] = ta_xxzzz_xxxz_0[i] * pa_y[i] - ta_xxzzz_xxxz_1[i] * pc_y[i];

        ta_xxyzzz_xxyy_0[i] =
            2.0 * ta_xxzzz_xxy_0[i] * fe_0 - 2.0 * ta_xxzzz_xxy_1[i] * fe_0 + ta_xxzzz_xxyy_0[i] * pa_y[i] - ta_xxzzz_xxyy_1[i] * pc_y[i];

        ta_xxyzzz_xxyz_0[i] = ta_xxzzz_xxz_0[i] * fe_0 - ta_xxzzz_xxz_1[i] * fe_0 + ta_xxzzz_xxyz_0[i] * pa_y[i] - ta_xxzzz_xxyz_1[i] * pc_y[i];

        ta_xxyzzz_xxzz_0[i] = ta_xxzzz_xxzz_0[i] * pa_y[i] - ta_xxzzz_xxzz_1[i] * pc_y[i];

        ta_xxyzzz_xyyy_0[i] =
            3.0 * ta_xxzzz_xyy_0[i] * fe_0 - 3.0 * ta_xxzzz_xyy_1[i] * fe_0 + ta_xxzzz_xyyy_0[i] * pa_y[i] - ta_xxzzz_xyyy_1[i] * pc_y[i];

        ta_xxyzzz_xyyz_0[i] =
            2.0 * ta_xxzzz_xyz_0[i] * fe_0 - 2.0 * ta_xxzzz_xyz_1[i] * fe_0 + ta_xxzzz_xyyz_0[i] * pa_y[i] - ta_xxzzz_xyyz_1[i] * pc_y[i];

        ta_xxyzzz_xyzz_0[i] = ta_xxzzz_xzz_0[i] * fe_0 - ta_xxzzz_xzz_1[i] * fe_0 + ta_xxzzz_xyzz_0[i] * pa_y[i] - ta_xxzzz_xyzz_1[i] * pc_y[i];

        ta_xxyzzz_xzzz_0[i] = ta_xxzzz_xzzz_0[i] * pa_y[i] - ta_xxzzz_xzzz_1[i] * pc_y[i];

        ta_xxyzzz_yyyy_0[i] = ta_yzzz_yyyy_0[i] * fe_0 - ta_yzzz_yyyy_1[i] * fe_0 + ta_xyzzz_yyyy_0[i] * pa_x[i] - ta_xyzzz_yyyy_1[i] * pc_x[i];

        ta_xxyzzz_yyyz_0[i] = ta_yzzz_yyyz_0[i] * fe_0 - ta_yzzz_yyyz_1[i] * fe_0 + ta_xyzzz_yyyz_0[i] * pa_x[i] - ta_xyzzz_yyyz_1[i] * pc_x[i];

        ta_xxyzzz_yyzz_0[i] = ta_yzzz_yyzz_0[i] * fe_0 - ta_yzzz_yyzz_1[i] * fe_0 + ta_xyzzz_yyzz_0[i] * pa_x[i] - ta_xyzzz_yyzz_1[i] * pc_x[i];

        ta_xxyzzz_yzzz_0[i] = ta_yzzz_yzzz_0[i] * fe_0 - ta_yzzz_yzzz_1[i] * fe_0 + ta_xyzzz_yzzz_0[i] * pa_x[i] - ta_xyzzz_yzzz_1[i] * pc_x[i];

        ta_xxyzzz_zzzz_0[i] = ta_xxzzz_zzzz_0[i] * pa_y[i] - ta_xxzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : IG

    auto ta_xxzzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 210);

    auto ta_xxzzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 211);

    auto ta_xxzzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 212);

    auto ta_xxzzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 213);

    auto ta_xxzzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 214);

    auto ta_xxzzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 215);

    auto ta_xxzzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 216);

    auto ta_xxzzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 217);

    auto ta_xxzzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 218);

    auto ta_xxzzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 219);

    auto ta_xxzzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 220);

    auto ta_xxzzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 221);

    auto ta_xxzzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 222);

    auto ta_xxzzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 223);

    auto ta_xxzzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 224);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxzz_xxxx_0,   \
                             ta_xxzz_xxxx_1,   \
                             ta_xxzz_xxxy_0,   \
                             ta_xxzz_xxxy_1,   \
                             ta_xxzz_xxyy_0,   \
                             ta_xxzz_xxyy_1,   \
                             ta_xxzz_xyyy_0,   \
                             ta_xxzz_xyyy_1,   \
                             ta_xxzzz_xxxx_0,  \
                             ta_xxzzz_xxxx_1,  \
                             ta_xxzzz_xxxy_0,  \
                             ta_xxzzz_xxxy_1,  \
                             ta_xxzzz_xxyy_0,  \
                             ta_xxzzz_xxyy_1,  \
                             ta_xxzzz_xyyy_0,  \
                             ta_xxzzz_xyyy_1,  \
                             ta_xxzzzz_xxxx_0, \
                             ta_xxzzzz_xxxy_0, \
                             ta_xxzzzz_xxxz_0, \
                             ta_xxzzzz_xxyy_0, \
                             ta_xxzzzz_xxyz_0, \
                             ta_xxzzzz_xxzz_0, \
                             ta_xxzzzz_xyyy_0, \
                             ta_xxzzzz_xyyz_0, \
                             ta_xxzzzz_xyzz_0, \
                             ta_xxzzzz_xzzz_0, \
                             ta_xxzzzz_yyyy_0, \
                             ta_xxzzzz_yyyz_0, \
                             ta_xxzzzz_yyzz_0, \
                             ta_xxzzzz_yzzz_0, \
                             ta_xxzzzz_zzzz_0, \
                             ta_xzzzz_xxxz_0,  \
                             ta_xzzzz_xxxz_1,  \
                             ta_xzzzz_xxyz_0,  \
                             ta_xzzzz_xxyz_1,  \
                             ta_xzzzz_xxz_0,   \
                             ta_xzzzz_xxz_1,   \
                             ta_xzzzz_xxzz_0,  \
                             ta_xzzzz_xxzz_1,  \
                             ta_xzzzz_xyyz_0,  \
                             ta_xzzzz_xyyz_1,  \
                             ta_xzzzz_xyz_0,   \
                             ta_xzzzz_xyz_1,   \
                             ta_xzzzz_xyzz_0,  \
                             ta_xzzzz_xyzz_1,  \
                             ta_xzzzz_xzz_0,   \
                             ta_xzzzz_xzz_1,   \
                             ta_xzzzz_xzzz_0,  \
                             ta_xzzzz_xzzz_1,  \
                             ta_xzzzz_yyyy_0,  \
                             ta_xzzzz_yyyy_1,  \
                             ta_xzzzz_yyyz_0,  \
                             ta_xzzzz_yyyz_1,  \
                             ta_xzzzz_yyz_0,   \
                             ta_xzzzz_yyz_1,   \
                             ta_xzzzz_yyzz_0,  \
                             ta_xzzzz_yyzz_1,  \
                             ta_xzzzz_yzz_0,   \
                             ta_xzzzz_yzz_1,   \
                             ta_xzzzz_yzzz_0,  \
                             ta_xzzzz_yzzz_1,  \
                             ta_xzzzz_zzz_0,   \
                             ta_xzzzz_zzz_1,   \
                             ta_xzzzz_zzzz_0,  \
                             ta_xzzzz_zzzz_1,  \
                             ta_zzzz_xxxz_0,   \
                             ta_zzzz_xxxz_1,   \
                             ta_zzzz_xxyz_0,   \
                             ta_zzzz_xxyz_1,   \
                             ta_zzzz_xxzz_0,   \
                             ta_zzzz_xxzz_1,   \
                             ta_zzzz_xyyz_0,   \
                             ta_zzzz_xyyz_1,   \
                             ta_zzzz_xyzz_0,   \
                             ta_zzzz_xyzz_1,   \
                             ta_zzzz_xzzz_0,   \
                             ta_zzzz_xzzz_1,   \
                             ta_zzzz_yyyy_0,   \
                             ta_zzzz_yyyy_1,   \
                             ta_zzzz_yyyz_0,   \
                             ta_zzzz_yyyz_1,   \
                             ta_zzzz_yyzz_0,   \
                             ta_zzzz_yyzz_1,   \
                             ta_zzzz_yzzz_0,   \
                             ta_zzzz_yzzz_1,   \
                             ta_zzzz_zzzz_0,   \
                             ta_zzzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzzz_xxxx_0[i] =
            3.0 * ta_xxzz_xxxx_0[i] * fe_0 - 3.0 * ta_xxzz_xxxx_1[i] * fe_0 + ta_xxzzz_xxxx_0[i] * pa_z[i] - ta_xxzzz_xxxx_1[i] * pc_z[i];

        ta_xxzzzz_xxxy_0[i] =
            3.0 * ta_xxzz_xxxy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxy_1[i] * fe_0 + ta_xxzzz_xxxy_0[i] * pa_z[i] - ta_xxzzz_xxxy_1[i] * pc_z[i];

        ta_xxzzzz_xxxz_0[i] = ta_zzzz_xxxz_0[i] * fe_0 - ta_zzzz_xxxz_1[i] * fe_0 + 3.0 * ta_xzzzz_xxz_0[i] * fe_0 - 3.0 * ta_xzzzz_xxz_1[i] * fe_0 +
                              ta_xzzzz_xxxz_0[i] * pa_x[i] - ta_xzzzz_xxxz_1[i] * pc_x[i];

        ta_xxzzzz_xxyy_0[i] =
            3.0 * ta_xxzz_xxyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxyy_1[i] * fe_0 + ta_xxzzz_xxyy_0[i] * pa_z[i] - ta_xxzzz_xxyy_1[i] * pc_z[i];

        ta_xxzzzz_xxyz_0[i] = ta_zzzz_xxyz_0[i] * fe_0 - ta_zzzz_xxyz_1[i] * fe_0 + 2.0 * ta_xzzzz_xyz_0[i] * fe_0 - 2.0 * ta_xzzzz_xyz_1[i] * fe_0 +
                              ta_xzzzz_xxyz_0[i] * pa_x[i] - ta_xzzzz_xxyz_1[i] * pc_x[i];

        ta_xxzzzz_xxzz_0[i] = ta_zzzz_xxzz_0[i] * fe_0 - ta_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta_xzzzz_xzz_0[i] * fe_0 - 2.0 * ta_xzzzz_xzz_1[i] * fe_0 +
                              ta_xzzzz_xxzz_0[i] * pa_x[i] - ta_xzzzz_xxzz_1[i] * pc_x[i];

        ta_xxzzzz_xyyy_0[i] =
            3.0 * ta_xxzz_xyyy_0[i] * fe_0 - 3.0 * ta_xxzz_xyyy_1[i] * fe_0 + ta_xxzzz_xyyy_0[i] * pa_z[i] - ta_xxzzz_xyyy_1[i] * pc_z[i];

        ta_xxzzzz_xyyz_0[i] = ta_zzzz_xyyz_0[i] * fe_0 - ta_zzzz_xyyz_1[i] * fe_0 + ta_xzzzz_yyz_0[i] * fe_0 - ta_xzzzz_yyz_1[i] * fe_0 +
                              ta_xzzzz_xyyz_0[i] * pa_x[i] - ta_xzzzz_xyyz_1[i] * pc_x[i];

        ta_xxzzzz_xyzz_0[i] = ta_zzzz_xyzz_0[i] * fe_0 - ta_zzzz_xyzz_1[i] * fe_0 + ta_xzzzz_yzz_0[i] * fe_0 - ta_xzzzz_yzz_1[i] * fe_0 +
                              ta_xzzzz_xyzz_0[i] * pa_x[i] - ta_xzzzz_xyzz_1[i] * pc_x[i];

        ta_xxzzzz_xzzz_0[i] = ta_zzzz_xzzz_0[i] * fe_0 - ta_zzzz_xzzz_1[i] * fe_0 + ta_xzzzz_zzz_0[i] * fe_0 - ta_xzzzz_zzz_1[i] * fe_0 +
                              ta_xzzzz_xzzz_0[i] * pa_x[i] - ta_xzzzz_xzzz_1[i] * pc_x[i];

        ta_xxzzzz_yyyy_0[i] = ta_zzzz_yyyy_0[i] * fe_0 - ta_zzzz_yyyy_1[i] * fe_0 + ta_xzzzz_yyyy_0[i] * pa_x[i] - ta_xzzzz_yyyy_1[i] * pc_x[i];

        ta_xxzzzz_yyyz_0[i] = ta_zzzz_yyyz_0[i] * fe_0 - ta_zzzz_yyyz_1[i] * fe_0 + ta_xzzzz_yyyz_0[i] * pa_x[i] - ta_xzzzz_yyyz_1[i] * pc_x[i];

        ta_xxzzzz_yyzz_0[i] = ta_zzzz_yyzz_0[i] * fe_0 - ta_zzzz_yyzz_1[i] * fe_0 + ta_xzzzz_yyzz_0[i] * pa_x[i] - ta_xzzzz_yyzz_1[i] * pc_x[i];

        ta_xxzzzz_yzzz_0[i] = ta_zzzz_yzzz_0[i] * fe_0 - ta_zzzz_yzzz_1[i] * fe_0 + ta_xzzzz_yzzz_0[i] * pa_x[i] - ta_xzzzz_yzzz_1[i] * pc_x[i];

        ta_xxzzzz_zzzz_0[i] = ta_zzzz_zzzz_0[i] * fe_0 - ta_zzzz_zzzz_1[i] * fe_0 + ta_xzzzz_zzzz_0[i] * pa_x[i] - ta_xzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : IG

    auto ta_xyyyyy_xxxx_0 = pbuffer.data(idx_npot_0_ig + 225);

    auto ta_xyyyyy_xxxy_0 = pbuffer.data(idx_npot_0_ig + 226);

    auto ta_xyyyyy_xxxz_0 = pbuffer.data(idx_npot_0_ig + 227);

    auto ta_xyyyyy_xxyy_0 = pbuffer.data(idx_npot_0_ig + 228);

    auto ta_xyyyyy_xxyz_0 = pbuffer.data(idx_npot_0_ig + 229);

    auto ta_xyyyyy_xxzz_0 = pbuffer.data(idx_npot_0_ig + 230);

    auto ta_xyyyyy_xyyy_0 = pbuffer.data(idx_npot_0_ig + 231);

    auto ta_xyyyyy_xyyz_0 = pbuffer.data(idx_npot_0_ig + 232);

    auto ta_xyyyyy_xyzz_0 = pbuffer.data(idx_npot_0_ig + 233);

    auto ta_xyyyyy_xzzz_0 = pbuffer.data(idx_npot_0_ig + 234);

    auto ta_xyyyyy_yyyy_0 = pbuffer.data(idx_npot_0_ig + 235);

    auto ta_xyyyyy_yyyz_0 = pbuffer.data(idx_npot_0_ig + 236);

    auto ta_xyyyyy_yyzz_0 = pbuffer.data(idx_npot_0_ig + 237);

    auto ta_xyyyyy_yzzz_0 = pbuffer.data(idx_npot_0_ig + 238);

    auto ta_xyyyyy_zzzz_0 = pbuffer.data(idx_npot_0_ig + 239);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xyyyyy_xxxx_0, \
                             ta_xyyyyy_xxxy_0, \
                             ta_xyyyyy_xxxz_0, \
                             ta_xyyyyy_xxyy_0, \
                             ta_xyyyyy_xxyz_0, \
                             ta_xyyyyy_xxzz_0, \
                             ta_xyyyyy_xyyy_0, \
                             ta_xyyyyy_xyyz_0, \
                             ta_xyyyyy_xyzz_0, \
                             ta_xyyyyy_xzzz_0, \
                             ta_xyyyyy_yyyy_0, \
                             ta_xyyyyy_yyyz_0, \
                             ta_xyyyyy_yyzz_0, \
                             ta_xyyyyy_yzzz_0, \
                             ta_xyyyyy_zzzz_0, \
                             ta_yyyyy_xxx_0,   \
                             ta_yyyyy_xxx_1,   \
                             ta_yyyyy_xxxx_0,  \
                             ta_yyyyy_xxxx_1,  \
                             ta_yyyyy_xxxy_0,  \
                             ta_yyyyy_xxxy_1,  \
                             ta_yyyyy_xxxz_0,  \
                             ta_yyyyy_xxxz_1,  \
                             ta_yyyyy_xxy_0,   \
                             ta_yyyyy_xxy_1,   \
                             ta_yyyyy_xxyy_0,  \
                             ta_yyyyy_xxyy_1,  \
                             ta_yyyyy_xxyz_0,  \
                             ta_yyyyy_xxyz_1,  \
                             ta_yyyyy_xxz_0,   \
                             ta_yyyyy_xxz_1,   \
                             ta_yyyyy_xxzz_0,  \
                             ta_yyyyy_xxzz_1,  \
                             ta_yyyyy_xyy_0,   \
                             ta_yyyyy_xyy_1,   \
                             ta_yyyyy_xyyy_0,  \
                             ta_yyyyy_xyyy_1,  \
                             ta_yyyyy_xyyz_0,  \
                             ta_yyyyy_xyyz_1,  \
                             ta_yyyyy_xyz_0,   \
                             ta_yyyyy_xyz_1,   \
                             ta_yyyyy_xyzz_0,  \
                             ta_yyyyy_xyzz_1,  \
                             ta_yyyyy_xzz_0,   \
                             ta_yyyyy_xzz_1,   \
                             ta_yyyyy_xzzz_0,  \
                             ta_yyyyy_xzzz_1,  \
                             ta_yyyyy_yyy_0,   \
                             ta_yyyyy_yyy_1,   \
                             ta_yyyyy_yyyy_0,  \
                             ta_yyyyy_yyyy_1,  \
                             ta_yyyyy_yyyz_0,  \
                             ta_yyyyy_yyyz_1,  \
                             ta_yyyyy_yyz_0,   \
                             ta_yyyyy_yyz_1,   \
                             ta_yyyyy_yyzz_0,  \
                             ta_yyyyy_yyzz_1,  \
                             ta_yyyyy_yzz_0,   \
                             ta_yyyyy_yzz_1,   \
                             ta_yyyyy_yzzz_0,  \
                             ta_yyyyy_yzzz_1,  \
                             ta_yyyyy_zzz_0,   \
                             ta_yyyyy_zzz_1,   \
                             ta_yyyyy_zzzz_0,  \
                             ta_yyyyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyy_xxxx_0[i] =
            4.0 * ta_yyyyy_xxx_0[i] * fe_0 - 4.0 * ta_yyyyy_xxx_1[i] * fe_0 + ta_yyyyy_xxxx_0[i] * pa_x[i] - ta_yyyyy_xxxx_1[i] * pc_x[i];

        ta_xyyyyy_xxxy_0[i] =
            3.0 * ta_yyyyy_xxy_0[i] * fe_0 - 3.0 * ta_yyyyy_xxy_1[i] * fe_0 + ta_yyyyy_xxxy_0[i] * pa_x[i] - ta_yyyyy_xxxy_1[i] * pc_x[i];

        ta_xyyyyy_xxxz_0[i] =
            3.0 * ta_yyyyy_xxz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxz_1[i] * fe_0 + ta_yyyyy_xxxz_0[i] * pa_x[i] - ta_yyyyy_xxxz_1[i] * pc_x[i];

        ta_xyyyyy_xxyy_0[i] =
            2.0 * ta_yyyyy_xyy_0[i] * fe_0 - 2.0 * ta_yyyyy_xyy_1[i] * fe_0 + ta_yyyyy_xxyy_0[i] * pa_x[i] - ta_yyyyy_xxyy_1[i] * pc_x[i];

        ta_xyyyyy_xxyz_0[i] =
            2.0 * ta_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyz_1[i] * fe_0 + ta_yyyyy_xxyz_0[i] * pa_x[i] - ta_yyyyy_xxyz_1[i] * pc_x[i];

        ta_xyyyyy_xxzz_0[i] =
            2.0 * ta_yyyyy_xzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xzz_1[i] * fe_0 + ta_yyyyy_xxzz_0[i] * pa_x[i] - ta_yyyyy_xxzz_1[i] * pc_x[i];

        ta_xyyyyy_xyyy_0[i] = ta_yyyyy_yyy_0[i] * fe_0 - ta_yyyyy_yyy_1[i] * fe_0 + ta_yyyyy_xyyy_0[i] * pa_x[i] - ta_yyyyy_xyyy_1[i] * pc_x[i];

        ta_xyyyyy_xyyz_0[i] = ta_yyyyy_yyz_0[i] * fe_0 - ta_yyyyy_yyz_1[i] * fe_0 + ta_yyyyy_xyyz_0[i] * pa_x[i] - ta_yyyyy_xyyz_1[i] * pc_x[i];

        ta_xyyyyy_xyzz_0[i] = ta_yyyyy_yzz_0[i] * fe_0 - ta_yyyyy_yzz_1[i] * fe_0 + ta_yyyyy_xyzz_0[i] * pa_x[i] - ta_yyyyy_xyzz_1[i] * pc_x[i];

        ta_xyyyyy_xzzz_0[i] = ta_yyyyy_zzz_0[i] * fe_0 - ta_yyyyy_zzz_1[i] * fe_0 + ta_yyyyy_xzzz_0[i] * pa_x[i] - ta_yyyyy_xzzz_1[i] * pc_x[i];

        ta_xyyyyy_yyyy_0[i] = ta_yyyyy_yyyy_0[i] * pa_x[i] - ta_yyyyy_yyyy_1[i] * pc_x[i];

        ta_xyyyyy_yyyz_0[i] = ta_yyyyy_yyyz_0[i] * pa_x[i] - ta_yyyyy_yyyz_1[i] * pc_x[i];

        ta_xyyyyy_yyzz_0[i] = ta_yyyyy_yyzz_0[i] * pa_x[i] - ta_yyyyy_yyzz_1[i] * pc_x[i];

        ta_xyyyyy_yzzz_0[i] = ta_yyyyy_yzzz_0[i] * pa_x[i] - ta_yyyyy_yzzz_1[i] * pc_x[i];

        ta_xyyyyy_zzzz_0[i] = ta_yyyyy_zzzz_0[i] * pa_x[i] - ta_yyyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : IG

    auto ta_xyyyyz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 240);

    auto ta_xyyyyz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 241);

    auto ta_xyyyyz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 242);

    auto ta_xyyyyz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 243);

    auto ta_xyyyyz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 244);

    auto ta_xyyyyz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 245);

    auto ta_xyyyyz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 246);

    auto ta_xyyyyz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 247);

    auto ta_xyyyyz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 248);

    auto ta_xyyyyz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 249);

    auto ta_xyyyyz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 250);

    auto ta_xyyyyz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 251);

    auto ta_xyyyyz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 252);

    auto ta_xyyyyz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 253);

    auto ta_xyyyyz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 254);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xyyyy_xxxx_0,  \
                             ta_xyyyy_xxxx_1,  \
                             ta_xyyyy_xxxy_0,  \
                             ta_xyyyy_xxxy_1,  \
                             ta_xyyyy_xxyy_0,  \
                             ta_xyyyy_xxyy_1,  \
                             ta_xyyyy_xyyy_0,  \
                             ta_xyyyy_xyyy_1,  \
                             ta_xyyyyz_xxxx_0, \
                             ta_xyyyyz_xxxy_0, \
                             ta_xyyyyz_xxxz_0, \
                             ta_xyyyyz_xxyy_0, \
                             ta_xyyyyz_xxyz_0, \
                             ta_xyyyyz_xxzz_0, \
                             ta_xyyyyz_xyyy_0, \
                             ta_xyyyyz_xyyz_0, \
                             ta_xyyyyz_xyzz_0, \
                             ta_xyyyyz_xzzz_0, \
                             ta_xyyyyz_yyyy_0, \
                             ta_xyyyyz_yyyz_0, \
                             ta_xyyyyz_yyzz_0, \
                             ta_xyyyyz_yzzz_0, \
                             ta_xyyyyz_zzzz_0, \
                             ta_yyyyz_xxxz_0,  \
                             ta_yyyyz_xxxz_1,  \
                             ta_yyyyz_xxyz_0,  \
                             ta_yyyyz_xxyz_1,  \
                             ta_yyyyz_xxz_0,   \
                             ta_yyyyz_xxz_1,   \
                             ta_yyyyz_xxzz_0,  \
                             ta_yyyyz_xxzz_1,  \
                             ta_yyyyz_xyyz_0,  \
                             ta_yyyyz_xyyz_1,  \
                             ta_yyyyz_xyz_0,   \
                             ta_yyyyz_xyz_1,   \
                             ta_yyyyz_xyzz_0,  \
                             ta_yyyyz_xyzz_1,  \
                             ta_yyyyz_xzz_0,   \
                             ta_yyyyz_xzz_1,   \
                             ta_yyyyz_xzzz_0,  \
                             ta_yyyyz_xzzz_1,  \
                             ta_yyyyz_yyyy_0,  \
                             ta_yyyyz_yyyy_1,  \
                             ta_yyyyz_yyyz_0,  \
                             ta_yyyyz_yyyz_1,  \
                             ta_yyyyz_yyz_0,   \
                             ta_yyyyz_yyz_1,   \
                             ta_yyyyz_yyzz_0,  \
                             ta_yyyyz_yyzz_1,  \
                             ta_yyyyz_yzz_0,   \
                             ta_yyyyz_yzz_1,   \
                             ta_yyyyz_yzzz_0,  \
                             ta_yyyyz_yzzz_1,  \
                             ta_yyyyz_zzz_0,   \
                             ta_yyyyz_zzz_1,   \
                             ta_yyyyz_zzzz_0,  \
                             ta_yyyyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyz_xxxx_0[i] = ta_xyyyy_xxxx_0[i] * pa_z[i] - ta_xyyyy_xxxx_1[i] * pc_z[i];

        ta_xyyyyz_xxxy_0[i] = ta_xyyyy_xxxy_0[i] * pa_z[i] - ta_xyyyy_xxxy_1[i] * pc_z[i];

        ta_xyyyyz_xxxz_0[i] =
            3.0 * ta_yyyyz_xxz_0[i] * fe_0 - 3.0 * ta_yyyyz_xxz_1[i] * fe_0 + ta_yyyyz_xxxz_0[i] * pa_x[i] - ta_yyyyz_xxxz_1[i] * pc_x[i];

        ta_xyyyyz_xxyy_0[i] = ta_xyyyy_xxyy_0[i] * pa_z[i] - ta_xyyyy_xxyy_1[i] * pc_z[i];

        ta_xyyyyz_xxyz_0[i] =
            2.0 * ta_yyyyz_xyz_0[i] * fe_0 - 2.0 * ta_yyyyz_xyz_1[i] * fe_0 + ta_yyyyz_xxyz_0[i] * pa_x[i] - ta_yyyyz_xxyz_1[i] * pc_x[i];

        ta_xyyyyz_xxzz_0[i] =
            2.0 * ta_yyyyz_xzz_0[i] * fe_0 - 2.0 * ta_yyyyz_xzz_1[i] * fe_0 + ta_yyyyz_xxzz_0[i] * pa_x[i] - ta_yyyyz_xxzz_1[i] * pc_x[i];

        ta_xyyyyz_xyyy_0[i] = ta_xyyyy_xyyy_0[i] * pa_z[i] - ta_xyyyy_xyyy_1[i] * pc_z[i];

        ta_xyyyyz_xyyz_0[i] = ta_yyyyz_yyz_0[i] * fe_0 - ta_yyyyz_yyz_1[i] * fe_0 + ta_yyyyz_xyyz_0[i] * pa_x[i] - ta_yyyyz_xyyz_1[i] * pc_x[i];

        ta_xyyyyz_xyzz_0[i] = ta_yyyyz_yzz_0[i] * fe_0 - ta_yyyyz_yzz_1[i] * fe_0 + ta_yyyyz_xyzz_0[i] * pa_x[i] - ta_yyyyz_xyzz_1[i] * pc_x[i];

        ta_xyyyyz_xzzz_0[i] = ta_yyyyz_zzz_0[i] * fe_0 - ta_yyyyz_zzz_1[i] * fe_0 + ta_yyyyz_xzzz_0[i] * pa_x[i] - ta_yyyyz_xzzz_1[i] * pc_x[i];

        ta_xyyyyz_yyyy_0[i] = ta_yyyyz_yyyy_0[i] * pa_x[i] - ta_yyyyz_yyyy_1[i] * pc_x[i];

        ta_xyyyyz_yyyz_0[i] = ta_yyyyz_yyyz_0[i] * pa_x[i] - ta_yyyyz_yyyz_1[i] * pc_x[i];

        ta_xyyyyz_yyzz_0[i] = ta_yyyyz_yyzz_0[i] * pa_x[i] - ta_yyyyz_yyzz_1[i] * pc_x[i];

        ta_xyyyyz_yzzz_0[i] = ta_yyyyz_yzzz_0[i] * pa_x[i] - ta_yyyyz_yzzz_1[i] * pc_x[i];

        ta_xyyyyz_zzzz_0[i] = ta_yyyyz_zzzz_0[i] * pa_x[i] - ta_yyyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 255-270 components of targeted buffer : IG

    auto ta_xyyyzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 255);

    auto ta_xyyyzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 256);

    auto ta_xyyyzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 257);

    auto ta_xyyyzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 258);

    auto ta_xyyyzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 259);

    auto ta_xyyyzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 260);

    auto ta_xyyyzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 261);

    auto ta_xyyyzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 262);

    auto ta_xyyyzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 263);

    auto ta_xyyyzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 264);

    auto ta_xyyyzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 265);

    auto ta_xyyyzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 266);

    auto ta_xyyyzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 267);

    auto ta_xyyyzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 268);

    auto ta_xyyyzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 269);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xyyyzz_xxxx_0, \
                             ta_xyyyzz_xxxy_0, \
                             ta_xyyyzz_xxxz_0, \
                             ta_xyyyzz_xxyy_0, \
                             ta_xyyyzz_xxyz_0, \
                             ta_xyyyzz_xxzz_0, \
                             ta_xyyyzz_xyyy_0, \
                             ta_xyyyzz_xyyz_0, \
                             ta_xyyyzz_xyzz_0, \
                             ta_xyyyzz_xzzz_0, \
                             ta_xyyyzz_yyyy_0, \
                             ta_xyyyzz_yyyz_0, \
                             ta_xyyyzz_yyzz_0, \
                             ta_xyyyzz_yzzz_0, \
                             ta_xyyyzz_zzzz_0, \
                             ta_yyyzz_xxx_0,   \
                             ta_yyyzz_xxx_1,   \
                             ta_yyyzz_xxxx_0,  \
                             ta_yyyzz_xxxx_1,  \
                             ta_yyyzz_xxxy_0,  \
                             ta_yyyzz_xxxy_1,  \
                             ta_yyyzz_xxxz_0,  \
                             ta_yyyzz_xxxz_1,  \
                             ta_yyyzz_xxy_0,   \
                             ta_yyyzz_xxy_1,   \
                             ta_yyyzz_xxyy_0,  \
                             ta_yyyzz_xxyy_1,  \
                             ta_yyyzz_xxyz_0,  \
                             ta_yyyzz_xxyz_1,  \
                             ta_yyyzz_xxz_0,   \
                             ta_yyyzz_xxz_1,   \
                             ta_yyyzz_xxzz_0,  \
                             ta_yyyzz_xxzz_1,  \
                             ta_yyyzz_xyy_0,   \
                             ta_yyyzz_xyy_1,   \
                             ta_yyyzz_xyyy_0,  \
                             ta_yyyzz_xyyy_1,  \
                             ta_yyyzz_xyyz_0,  \
                             ta_yyyzz_xyyz_1,  \
                             ta_yyyzz_xyz_0,   \
                             ta_yyyzz_xyz_1,   \
                             ta_yyyzz_xyzz_0,  \
                             ta_yyyzz_xyzz_1,  \
                             ta_yyyzz_xzz_0,   \
                             ta_yyyzz_xzz_1,   \
                             ta_yyyzz_xzzz_0,  \
                             ta_yyyzz_xzzz_1,  \
                             ta_yyyzz_yyy_0,   \
                             ta_yyyzz_yyy_1,   \
                             ta_yyyzz_yyyy_0,  \
                             ta_yyyzz_yyyy_1,  \
                             ta_yyyzz_yyyz_0,  \
                             ta_yyyzz_yyyz_1,  \
                             ta_yyyzz_yyz_0,   \
                             ta_yyyzz_yyz_1,   \
                             ta_yyyzz_yyzz_0,  \
                             ta_yyyzz_yyzz_1,  \
                             ta_yyyzz_yzz_0,   \
                             ta_yyyzz_yzz_1,   \
                             ta_yyyzz_yzzz_0,  \
                             ta_yyyzz_yzzz_1,  \
                             ta_yyyzz_zzz_0,   \
                             ta_yyyzz_zzz_1,   \
                             ta_yyyzz_zzzz_0,  \
                             ta_yyyzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyzz_xxxx_0[i] =
            4.0 * ta_yyyzz_xxx_0[i] * fe_0 - 4.0 * ta_yyyzz_xxx_1[i] * fe_0 + ta_yyyzz_xxxx_0[i] * pa_x[i] - ta_yyyzz_xxxx_1[i] * pc_x[i];

        ta_xyyyzz_xxxy_0[i] =
            3.0 * ta_yyyzz_xxy_0[i] * fe_0 - 3.0 * ta_yyyzz_xxy_1[i] * fe_0 + ta_yyyzz_xxxy_0[i] * pa_x[i] - ta_yyyzz_xxxy_1[i] * pc_x[i];

        ta_xyyyzz_xxxz_0[i] =
            3.0 * ta_yyyzz_xxz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxz_1[i] * fe_0 + ta_yyyzz_xxxz_0[i] * pa_x[i] - ta_yyyzz_xxxz_1[i] * pc_x[i];

        ta_xyyyzz_xxyy_0[i] =
            2.0 * ta_yyyzz_xyy_0[i] * fe_0 - 2.0 * ta_yyyzz_xyy_1[i] * fe_0 + ta_yyyzz_xxyy_0[i] * pa_x[i] - ta_yyyzz_xxyy_1[i] * pc_x[i];

        ta_xyyyzz_xxyz_0[i] =
            2.0 * ta_yyyzz_xyz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyz_1[i] * fe_0 + ta_yyyzz_xxyz_0[i] * pa_x[i] - ta_yyyzz_xxyz_1[i] * pc_x[i];

        ta_xyyyzz_xxzz_0[i] =
            2.0 * ta_yyyzz_xzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xzz_1[i] * fe_0 + ta_yyyzz_xxzz_0[i] * pa_x[i] - ta_yyyzz_xxzz_1[i] * pc_x[i];

        ta_xyyyzz_xyyy_0[i] = ta_yyyzz_yyy_0[i] * fe_0 - ta_yyyzz_yyy_1[i] * fe_0 + ta_yyyzz_xyyy_0[i] * pa_x[i] - ta_yyyzz_xyyy_1[i] * pc_x[i];

        ta_xyyyzz_xyyz_0[i] = ta_yyyzz_yyz_0[i] * fe_0 - ta_yyyzz_yyz_1[i] * fe_0 + ta_yyyzz_xyyz_0[i] * pa_x[i] - ta_yyyzz_xyyz_1[i] * pc_x[i];

        ta_xyyyzz_xyzz_0[i] = ta_yyyzz_yzz_0[i] * fe_0 - ta_yyyzz_yzz_1[i] * fe_0 + ta_yyyzz_xyzz_0[i] * pa_x[i] - ta_yyyzz_xyzz_1[i] * pc_x[i];

        ta_xyyyzz_xzzz_0[i] = ta_yyyzz_zzz_0[i] * fe_0 - ta_yyyzz_zzz_1[i] * fe_0 + ta_yyyzz_xzzz_0[i] * pa_x[i] - ta_yyyzz_xzzz_1[i] * pc_x[i];

        ta_xyyyzz_yyyy_0[i] = ta_yyyzz_yyyy_0[i] * pa_x[i] - ta_yyyzz_yyyy_1[i] * pc_x[i];

        ta_xyyyzz_yyyz_0[i] = ta_yyyzz_yyyz_0[i] * pa_x[i] - ta_yyyzz_yyyz_1[i] * pc_x[i];

        ta_xyyyzz_yyzz_0[i] = ta_yyyzz_yyzz_0[i] * pa_x[i] - ta_yyyzz_yyzz_1[i] * pc_x[i];

        ta_xyyyzz_yzzz_0[i] = ta_yyyzz_yzzz_0[i] * pa_x[i] - ta_yyyzz_yzzz_1[i] * pc_x[i];

        ta_xyyyzz_zzzz_0[i] = ta_yyyzz_zzzz_0[i] * pa_x[i] - ta_yyyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 270-285 components of targeted buffer : IG

    auto ta_xyyzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 270);

    auto ta_xyyzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 271);

    auto ta_xyyzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 272);

    auto ta_xyyzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 273);

    auto ta_xyyzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 274);

    auto ta_xyyzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 275);

    auto ta_xyyzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 276);

    auto ta_xyyzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 277);

    auto ta_xyyzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 278);

    auto ta_xyyzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 279);

    auto ta_xyyzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 280);

    auto ta_xyyzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 281);

    auto ta_xyyzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 282);

    auto ta_xyyzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 283);

    auto ta_xyyzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 284);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xyyzzz_xxxx_0, \
                             ta_xyyzzz_xxxy_0, \
                             ta_xyyzzz_xxxz_0, \
                             ta_xyyzzz_xxyy_0, \
                             ta_xyyzzz_xxyz_0, \
                             ta_xyyzzz_xxzz_0, \
                             ta_xyyzzz_xyyy_0, \
                             ta_xyyzzz_xyyz_0, \
                             ta_xyyzzz_xyzz_0, \
                             ta_xyyzzz_xzzz_0, \
                             ta_xyyzzz_yyyy_0, \
                             ta_xyyzzz_yyyz_0, \
                             ta_xyyzzz_yyzz_0, \
                             ta_xyyzzz_yzzz_0, \
                             ta_xyyzzz_zzzz_0, \
                             ta_yyzzz_xxx_0,   \
                             ta_yyzzz_xxx_1,   \
                             ta_yyzzz_xxxx_0,  \
                             ta_yyzzz_xxxx_1,  \
                             ta_yyzzz_xxxy_0,  \
                             ta_yyzzz_xxxy_1,  \
                             ta_yyzzz_xxxz_0,  \
                             ta_yyzzz_xxxz_1,  \
                             ta_yyzzz_xxy_0,   \
                             ta_yyzzz_xxy_1,   \
                             ta_yyzzz_xxyy_0,  \
                             ta_yyzzz_xxyy_1,  \
                             ta_yyzzz_xxyz_0,  \
                             ta_yyzzz_xxyz_1,  \
                             ta_yyzzz_xxz_0,   \
                             ta_yyzzz_xxz_1,   \
                             ta_yyzzz_xxzz_0,  \
                             ta_yyzzz_xxzz_1,  \
                             ta_yyzzz_xyy_0,   \
                             ta_yyzzz_xyy_1,   \
                             ta_yyzzz_xyyy_0,  \
                             ta_yyzzz_xyyy_1,  \
                             ta_yyzzz_xyyz_0,  \
                             ta_yyzzz_xyyz_1,  \
                             ta_yyzzz_xyz_0,   \
                             ta_yyzzz_xyz_1,   \
                             ta_yyzzz_xyzz_0,  \
                             ta_yyzzz_xyzz_1,  \
                             ta_yyzzz_xzz_0,   \
                             ta_yyzzz_xzz_1,   \
                             ta_yyzzz_xzzz_0,  \
                             ta_yyzzz_xzzz_1,  \
                             ta_yyzzz_yyy_0,   \
                             ta_yyzzz_yyy_1,   \
                             ta_yyzzz_yyyy_0,  \
                             ta_yyzzz_yyyy_1,  \
                             ta_yyzzz_yyyz_0,  \
                             ta_yyzzz_yyyz_1,  \
                             ta_yyzzz_yyz_0,   \
                             ta_yyzzz_yyz_1,   \
                             ta_yyzzz_yyzz_0,  \
                             ta_yyzzz_yyzz_1,  \
                             ta_yyzzz_yzz_0,   \
                             ta_yyzzz_yzz_1,   \
                             ta_yyzzz_yzzz_0,  \
                             ta_yyzzz_yzzz_1,  \
                             ta_yyzzz_zzz_0,   \
                             ta_yyzzz_zzz_1,   \
                             ta_yyzzz_zzzz_0,  \
                             ta_yyzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzzz_xxxx_0[i] =
            4.0 * ta_yyzzz_xxx_0[i] * fe_0 - 4.0 * ta_yyzzz_xxx_1[i] * fe_0 + ta_yyzzz_xxxx_0[i] * pa_x[i] - ta_yyzzz_xxxx_1[i] * pc_x[i];

        ta_xyyzzz_xxxy_0[i] =
            3.0 * ta_yyzzz_xxy_0[i] * fe_0 - 3.0 * ta_yyzzz_xxy_1[i] * fe_0 + ta_yyzzz_xxxy_0[i] * pa_x[i] - ta_yyzzz_xxxy_1[i] * pc_x[i];

        ta_xyyzzz_xxxz_0[i] =
            3.0 * ta_yyzzz_xxz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxz_1[i] * fe_0 + ta_yyzzz_xxxz_0[i] * pa_x[i] - ta_yyzzz_xxxz_1[i] * pc_x[i];

        ta_xyyzzz_xxyy_0[i] =
            2.0 * ta_yyzzz_xyy_0[i] * fe_0 - 2.0 * ta_yyzzz_xyy_1[i] * fe_0 + ta_yyzzz_xxyy_0[i] * pa_x[i] - ta_yyzzz_xxyy_1[i] * pc_x[i];

        ta_xyyzzz_xxyz_0[i] =
            2.0 * ta_yyzzz_xyz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyz_1[i] * fe_0 + ta_yyzzz_xxyz_0[i] * pa_x[i] - ta_yyzzz_xxyz_1[i] * pc_x[i];

        ta_xyyzzz_xxzz_0[i] =
            2.0 * ta_yyzzz_xzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xzz_1[i] * fe_0 + ta_yyzzz_xxzz_0[i] * pa_x[i] - ta_yyzzz_xxzz_1[i] * pc_x[i];

        ta_xyyzzz_xyyy_0[i] = ta_yyzzz_yyy_0[i] * fe_0 - ta_yyzzz_yyy_1[i] * fe_0 + ta_yyzzz_xyyy_0[i] * pa_x[i] - ta_yyzzz_xyyy_1[i] * pc_x[i];

        ta_xyyzzz_xyyz_0[i] = ta_yyzzz_yyz_0[i] * fe_0 - ta_yyzzz_yyz_1[i] * fe_0 + ta_yyzzz_xyyz_0[i] * pa_x[i] - ta_yyzzz_xyyz_1[i] * pc_x[i];

        ta_xyyzzz_xyzz_0[i] = ta_yyzzz_yzz_0[i] * fe_0 - ta_yyzzz_yzz_1[i] * fe_0 + ta_yyzzz_xyzz_0[i] * pa_x[i] - ta_yyzzz_xyzz_1[i] * pc_x[i];

        ta_xyyzzz_xzzz_0[i] = ta_yyzzz_zzz_0[i] * fe_0 - ta_yyzzz_zzz_1[i] * fe_0 + ta_yyzzz_xzzz_0[i] * pa_x[i] - ta_yyzzz_xzzz_1[i] * pc_x[i];

        ta_xyyzzz_yyyy_0[i] = ta_yyzzz_yyyy_0[i] * pa_x[i] - ta_yyzzz_yyyy_1[i] * pc_x[i];

        ta_xyyzzz_yyyz_0[i] = ta_yyzzz_yyyz_0[i] * pa_x[i] - ta_yyzzz_yyyz_1[i] * pc_x[i];

        ta_xyyzzz_yyzz_0[i] = ta_yyzzz_yyzz_0[i] * pa_x[i] - ta_yyzzz_yyzz_1[i] * pc_x[i];

        ta_xyyzzz_yzzz_0[i] = ta_yyzzz_yzzz_0[i] * pa_x[i] - ta_yyzzz_yzzz_1[i] * pc_x[i];

        ta_xyyzzz_zzzz_0[i] = ta_yyzzz_zzzz_0[i] * pa_x[i] - ta_yyzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 285-300 components of targeted buffer : IG

    auto ta_xyzzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 285);

    auto ta_xyzzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 286);

    auto ta_xyzzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 287);

    auto ta_xyzzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 288);

    auto ta_xyzzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 289);

    auto ta_xyzzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 290);

    auto ta_xyzzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 291);

    auto ta_xyzzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 292);

    auto ta_xyzzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 293);

    auto ta_xyzzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 294);

    auto ta_xyzzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 295);

    auto ta_xyzzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 296);

    auto ta_xyzzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 297);

    auto ta_xyzzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 298);

    auto ta_xyzzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 299);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xyzzzz_xxxx_0, \
                             ta_xyzzzz_xxxy_0, \
                             ta_xyzzzz_xxxz_0, \
                             ta_xyzzzz_xxyy_0, \
                             ta_xyzzzz_xxyz_0, \
                             ta_xyzzzz_xxzz_0, \
                             ta_xyzzzz_xyyy_0, \
                             ta_xyzzzz_xyyz_0, \
                             ta_xyzzzz_xyzz_0, \
                             ta_xyzzzz_xzzz_0, \
                             ta_xyzzzz_yyyy_0, \
                             ta_xyzzzz_yyyz_0, \
                             ta_xyzzzz_yyzz_0, \
                             ta_xyzzzz_yzzz_0, \
                             ta_xyzzzz_zzzz_0, \
                             ta_xzzzz_xxxx_0,  \
                             ta_xzzzz_xxxx_1,  \
                             ta_xzzzz_xxxz_0,  \
                             ta_xzzzz_xxxz_1,  \
                             ta_xzzzz_xxzz_0,  \
                             ta_xzzzz_xxzz_1,  \
                             ta_xzzzz_xzzz_0,  \
                             ta_xzzzz_xzzz_1,  \
                             ta_yzzzz_xxxy_0,  \
                             ta_yzzzz_xxxy_1,  \
                             ta_yzzzz_xxy_0,   \
                             ta_yzzzz_xxy_1,   \
                             ta_yzzzz_xxyy_0,  \
                             ta_yzzzz_xxyy_1,  \
                             ta_yzzzz_xxyz_0,  \
                             ta_yzzzz_xxyz_1,  \
                             ta_yzzzz_xyy_0,   \
                             ta_yzzzz_xyy_1,   \
                             ta_yzzzz_xyyy_0,  \
                             ta_yzzzz_xyyy_1,  \
                             ta_yzzzz_xyyz_0,  \
                             ta_yzzzz_xyyz_1,  \
                             ta_yzzzz_xyz_0,   \
                             ta_yzzzz_xyz_1,   \
                             ta_yzzzz_xyzz_0,  \
                             ta_yzzzz_xyzz_1,  \
                             ta_yzzzz_yyy_0,   \
                             ta_yzzzz_yyy_1,   \
                             ta_yzzzz_yyyy_0,  \
                             ta_yzzzz_yyyy_1,  \
                             ta_yzzzz_yyyz_0,  \
                             ta_yzzzz_yyyz_1,  \
                             ta_yzzzz_yyz_0,   \
                             ta_yzzzz_yyz_1,   \
                             ta_yzzzz_yyzz_0,  \
                             ta_yzzzz_yyzz_1,  \
                             ta_yzzzz_yzz_0,   \
                             ta_yzzzz_yzz_1,   \
                             ta_yzzzz_yzzz_0,  \
                             ta_yzzzz_yzzz_1,  \
                             ta_yzzzz_zzzz_0,  \
                             ta_yzzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzzz_xxxx_0[i] = ta_xzzzz_xxxx_0[i] * pa_y[i] - ta_xzzzz_xxxx_1[i] * pc_y[i];

        ta_xyzzzz_xxxy_0[i] =
            3.0 * ta_yzzzz_xxy_0[i] * fe_0 - 3.0 * ta_yzzzz_xxy_1[i] * fe_0 + ta_yzzzz_xxxy_0[i] * pa_x[i] - ta_yzzzz_xxxy_1[i] * pc_x[i];

        ta_xyzzzz_xxxz_0[i] = ta_xzzzz_xxxz_0[i] * pa_y[i] - ta_xzzzz_xxxz_1[i] * pc_y[i];

        ta_xyzzzz_xxyy_0[i] =
            2.0 * ta_yzzzz_xyy_0[i] * fe_0 - 2.0 * ta_yzzzz_xyy_1[i] * fe_0 + ta_yzzzz_xxyy_0[i] * pa_x[i] - ta_yzzzz_xxyy_1[i] * pc_x[i];

        ta_xyzzzz_xxyz_0[i] =
            2.0 * ta_yzzzz_xyz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyz_1[i] * fe_0 + ta_yzzzz_xxyz_0[i] * pa_x[i] - ta_yzzzz_xxyz_1[i] * pc_x[i];

        ta_xyzzzz_xxzz_0[i] = ta_xzzzz_xxzz_0[i] * pa_y[i] - ta_xzzzz_xxzz_1[i] * pc_y[i];

        ta_xyzzzz_xyyy_0[i] = ta_yzzzz_yyy_0[i] * fe_0 - ta_yzzzz_yyy_1[i] * fe_0 + ta_yzzzz_xyyy_0[i] * pa_x[i] - ta_yzzzz_xyyy_1[i] * pc_x[i];

        ta_xyzzzz_xyyz_0[i] = ta_yzzzz_yyz_0[i] * fe_0 - ta_yzzzz_yyz_1[i] * fe_0 + ta_yzzzz_xyyz_0[i] * pa_x[i] - ta_yzzzz_xyyz_1[i] * pc_x[i];

        ta_xyzzzz_xyzz_0[i] = ta_yzzzz_yzz_0[i] * fe_0 - ta_yzzzz_yzz_1[i] * fe_0 + ta_yzzzz_xyzz_0[i] * pa_x[i] - ta_yzzzz_xyzz_1[i] * pc_x[i];

        ta_xyzzzz_xzzz_0[i] = ta_xzzzz_xzzz_0[i] * pa_y[i] - ta_xzzzz_xzzz_1[i] * pc_y[i];

        ta_xyzzzz_yyyy_0[i] = ta_yzzzz_yyyy_0[i] * pa_x[i] - ta_yzzzz_yyyy_1[i] * pc_x[i];

        ta_xyzzzz_yyyz_0[i] = ta_yzzzz_yyyz_0[i] * pa_x[i] - ta_yzzzz_yyyz_1[i] * pc_x[i];

        ta_xyzzzz_yyzz_0[i] = ta_yzzzz_yyzz_0[i] * pa_x[i] - ta_yzzzz_yyzz_1[i] * pc_x[i];

        ta_xyzzzz_yzzz_0[i] = ta_yzzzz_yzzz_0[i] * pa_x[i] - ta_yzzzz_yzzz_1[i] * pc_x[i];

        ta_xyzzzz_zzzz_0[i] = ta_yzzzz_zzzz_0[i] * pa_x[i] - ta_yzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 300-315 components of targeted buffer : IG

    auto ta_xzzzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 300);

    auto ta_xzzzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 301);

    auto ta_xzzzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 302);

    auto ta_xzzzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 303);

    auto ta_xzzzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 304);

    auto ta_xzzzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 305);

    auto ta_xzzzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 306);

    auto ta_xzzzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 307);

    auto ta_xzzzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 308);

    auto ta_xzzzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 309);

    auto ta_xzzzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 310);

    auto ta_xzzzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 311);

    auto ta_xzzzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 312);

    auto ta_xzzzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 313);

    auto ta_xzzzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 314);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xzzzzz_xxxx_0, \
                             ta_xzzzzz_xxxy_0, \
                             ta_xzzzzz_xxxz_0, \
                             ta_xzzzzz_xxyy_0, \
                             ta_xzzzzz_xxyz_0, \
                             ta_xzzzzz_xxzz_0, \
                             ta_xzzzzz_xyyy_0, \
                             ta_xzzzzz_xyyz_0, \
                             ta_xzzzzz_xyzz_0, \
                             ta_xzzzzz_xzzz_0, \
                             ta_xzzzzz_yyyy_0, \
                             ta_xzzzzz_yyyz_0, \
                             ta_xzzzzz_yyzz_0, \
                             ta_xzzzzz_yzzz_0, \
                             ta_xzzzzz_zzzz_0, \
                             ta_zzzzz_xxx_0,   \
                             ta_zzzzz_xxx_1,   \
                             ta_zzzzz_xxxx_0,  \
                             ta_zzzzz_xxxx_1,  \
                             ta_zzzzz_xxxy_0,  \
                             ta_zzzzz_xxxy_1,  \
                             ta_zzzzz_xxxz_0,  \
                             ta_zzzzz_xxxz_1,  \
                             ta_zzzzz_xxy_0,   \
                             ta_zzzzz_xxy_1,   \
                             ta_zzzzz_xxyy_0,  \
                             ta_zzzzz_xxyy_1,  \
                             ta_zzzzz_xxyz_0,  \
                             ta_zzzzz_xxyz_1,  \
                             ta_zzzzz_xxz_0,   \
                             ta_zzzzz_xxz_1,   \
                             ta_zzzzz_xxzz_0,  \
                             ta_zzzzz_xxzz_1,  \
                             ta_zzzzz_xyy_0,   \
                             ta_zzzzz_xyy_1,   \
                             ta_zzzzz_xyyy_0,  \
                             ta_zzzzz_xyyy_1,  \
                             ta_zzzzz_xyyz_0,  \
                             ta_zzzzz_xyyz_1,  \
                             ta_zzzzz_xyz_0,   \
                             ta_zzzzz_xyz_1,   \
                             ta_zzzzz_xyzz_0,  \
                             ta_zzzzz_xyzz_1,  \
                             ta_zzzzz_xzz_0,   \
                             ta_zzzzz_xzz_1,   \
                             ta_zzzzz_xzzz_0,  \
                             ta_zzzzz_xzzz_1,  \
                             ta_zzzzz_yyy_0,   \
                             ta_zzzzz_yyy_1,   \
                             ta_zzzzz_yyyy_0,  \
                             ta_zzzzz_yyyy_1,  \
                             ta_zzzzz_yyyz_0,  \
                             ta_zzzzz_yyyz_1,  \
                             ta_zzzzz_yyz_0,   \
                             ta_zzzzz_yyz_1,   \
                             ta_zzzzz_yyzz_0,  \
                             ta_zzzzz_yyzz_1,  \
                             ta_zzzzz_yzz_0,   \
                             ta_zzzzz_yzz_1,   \
                             ta_zzzzz_yzzz_0,  \
                             ta_zzzzz_yzzz_1,  \
                             ta_zzzzz_zzz_0,   \
                             ta_zzzzz_zzz_1,   \
                             ta_zzzzz_zzzz_0,  \
                             ta_zzzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzzz_xxxx_0[i] =
            4.0 * ta_zzzzz_xxx_0[i] * fe_0 - 4.0 * ta_zzzzz_xxx_1[i] * fe_0 + ta_zzzzz_xxxx_0[i] * pa_x[i] - ta_zzzzz_xxxx_1[i] * pc_x[i];

        ta_xzzzzz_xxxy_0[i] =
            3.0 * ta_zzzzz_xxy_0[i] * fe_0 - 3.0 * ta_zzzzz_xxy_1[i] * fe_0 + ta_zzzzz_xxxy_0[i] * pa_x[i] - ta_zzzzz_xxxy_1[i] * pc_x[i];

        ta_xzzzzz_xxxz_0[i] =
            3.0 * ta_zzzzz_xxz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxz_1[i] * fe_0 + ta_zzzzz_xxxz_0[i] * pa_x[i] - ta_zzzzz_xxxz_1[i] * pc_x[i];

        ta_xzzzzz_xxyy_0[i] =
            2.0 * ta_zzzzz_xyy_0[i] * fe_0 - 2.0 * ta_zzzzz_xyy_1[i] * fe_0 + ta_zzzzz_xxyy_0[i] * pa_x[i] - ta_zzzzz_xxyy_1[i] * pc_x[i];

        ta_xzzzzz_xxyz_0[i] =
            2.0 * ta_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyz_1[i] * fe_0 + ta_zzzzz_xxyz_0[i] * pa_x[i] - ta_zzzzz_xxyz_1[i] * pc_x[i];

        ta_xzzzzz_xxzz_0[i] =
            2.0 * ta_zzzzz_xzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xzz_1[i] * fe_0 + ta_zzzzz_xxzz_0[i] * pa_x[i] - ta_zzzzz_xxzz_1[i] * pc_x[i];

        ta_xzzzzz_xyyy_0[i] = ta_zzzzz_yyy_0[i] * fe_0 - ta_zzzzz_yyy_1[i] * fe_0 + ta_zzzzz_xyyy_0[i] * pa_x[i] - ta_zzzzz_xyyy_1[i] * pc_x[i];

        ta_xzzzzz_xyyz_0[i] = ta_zzzzz_yyz_0[i] * fe_0 - ta_zzzzz_yyz_1[i] * fe_0 + ta_zzzzz_xyyz_0[i] * pa_x[i] - ta_zzzzz_xyyz_1[i] * pc_x[i];

        ta_xzzzzz_xyzz_0[i] = ta_zzzzz_yzz_0[i] * fe_0 - ta_zzzzz_yzz_1[i] * fe_0 + ta_zzzzz_xyzz_0[i] * pa_x[i] - ta_zzzzz_xyzz_1[i] * pc_x[i];

        ta_xzzzzz_xzzz_0[i] = ta_zzzzz_zzz_0[i] * fe_0 - ta_zzzzz_zzz_1[i] * fe_0 + ta_zzzzz_xzzz_0[i] * pa_x[i] - ta_zzzzz_xzzz_1[i] * pc_x[i];

        ta_xzzzzz_yyyy_0[i] = ta_zzzzz_yyyy_0[i] * pa_x[i] - ta_zzzzz_yyyy_1[i] * pc_x[i];

        ta_xzzzzz_yyyz_0[i] = ta_zzzzz_yyyz_0[i] * pa_x[i] - ta_zzzzz_yyyz_1[i] * pc_x[i];

        ta_xzzzzz_yyzz_0[i] = ta_zzzzz_yyzz_0[i] * pa_x[i] - ta_zzzzz_yyzz_1[i] * pc_x[i];

        ta_xzzzzz_yzzz_0[i] = ta_zzzzz_yzzz_0[i] * pa_x[i] - ta_zzzzz_yzzz_1[i] * pc_x[i];

        ta_xzzzzz_zzzz_0[i] = ta_zzzzz_zzzz_0[i] * pa_x[i] - ta_zzzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : IG

    auto ta_yyyyyy_xxxx_0 = pbuffer.data(idx_npot_0_ig + 315);

    auto ta_yyyyyy_xxxy_0 = pbuffer.data(idx_npot_0_ig + 316);

    auto ta_yyyyyy_xxxz_0 = pbuffer.data(idx_npot_0_ig + 317);

    auto ta_yyyyyy_xxyy_0 = pbuffer.data(idx_npot_0_ig + 318);

    auto ta_yyyyyy_xxyz_0 = pbuffer.data(idx_npot_0_ig + 319);

    auto ta_yyyyyy_xxzz_0 = pbuffer.data(idx_npot_0_ig + 320);

    auto ta_yyyyyy_xyyy_0 = pbuffer.data(idx_npot_0_ig + 321);

    auto ta_yyyyyy_xyyz_0 = pbuffer.data(idx_npot_0_ig + 322);

    auto ta_yyyyyy_xyzz_0 = pbuffer.data(idx_npot_0_ig + 323);

    auto ta_yyyyyy_xzzz_0 = pbuffer.data(idx_npot_0_ig + 324);

    auto ta_yyyyyy_yyyy_0 = pbuffer.data(idx_npot_0_ig + 325);

    auto ta_yyyyyy_yyyz_0 = pbuffer.data(idx_npot_0_ig + 326);

    auto ta_yyyyyy_yyzz_0 = pbuffer.data(idx_npot_0_ig + 327);

    auto ta_yyyyyy_yzzz_0 = pbuffer.data(idx_npot_0_ig + 328);

    auto ta_yyyyyy_zzzz_0 = pbuffer.data(idx_npot_0_ig + 329);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta_yyyy_xxxx_0,   \
                             ta_yyyy_xxxx_1,   \
                             ta_yyyy_xxxy_0,   \
                             ta_yyyy_xxxy_1,   \
                             ta_yyyy_xxxz_0,   \
                             ta_yyyy_xxxz_1,   \
                             ta_yyyy_xxyy_0,   \
                             ta_yyyy_xxyy_1,   \
                             ta_yyyy_xxyz_0,   \
                             ta_yyyy_xxyz_1,   \
                             ta_yyyy_xxzz_0,   \
                             ta_yyyy_xxzz_1,   \
                             ta_yyyy_xyyy_0,   \
                             ta_yyyy_xyyy_1,   \
                             ta_yyyy_xyyz_0,   \
                             ta_yyyy_xyyz_1,   \
                             ta_yyyy_xyzz_0,   \
                             ta_yyyy_xyzz_1,   \
                             ta_yyyy_xzzz_0,   \
                             ta_yyyy_xzzz_1,   \
                             ta_yyyy_yyyy_0,   \
                             ta_yyyy_yyyy_1,   \
                             ta_yyyy_yyyz_0,   \
                             ta_yyyy_yyyz_1,   \
                             ta_yyyy_yyzz_0,   \
                             ta_yyyy_yyzz_1,   \
                             ta_yyyy_yzzz_0,   \
                             ta_yyyy_yzzz_1,   \
                             ta_yyyy_zzzz_0,   \
                             ta_yyyy_zzzz_1,   \
                             ta_yyyyy_xxx_0,   \
                             ta_yyyyy_xxx_1,   \
                             ta_yyyyy_xxxx_0,  \
                             ta_yyyyy_xxxx_1,  \
                             ta_yyyyy_xxxy_0,  \
                             ta_yyyyy_xxxy_1,  \
                             ta_yyyyy_xxxz_0,  \
                             ta_yyyyy_xxxz_1,  \
                             ta_yyyyy_xxy_0,   \
                             ta_yyyyy_xxy_1,   \
                             ta_yyyyy_xxyy_0,  \
                             ta_yyyyy_xxyy_1,  \
                             ta_yyyyy_xxyz_0,  \
                             ta_yyyyy_xxyz_1,  \
                             ta_yyyyy_xxz_0,   \
                             ta_yyyyy_xxz_1,   \
                             ta_yyyyy_xxzz_0,  \
                             ta_yyyyy_xxzz_1,  \
                             ta_yyyyy_xyy_0,   \
                             ta_yyyyy_xyy_1,   \
                             ta_yyyyy_xyyy_0,  \
                             ta_yyyyy_xyyy_1,  \
                             ta_yyyyy_xyyz_0,  \
                             ta_yyyyy_xyyz_1,  \
                             ta_yyyyy_xyz_0,   \
                             ta_yyyyy_xyz_1,   \
                             ta_yyyyy_xyzz_0,  \
                             ta_yyyyy_xyzz_1,  \
                             ta_yyyyy_xzz_0,   \
                             ta_yyyyy_xzz_1,   \
                             ta_yyyyy_xzzz_0,  \
                             ta_yyyyy_xzzz_1,  \
                             ta_yyyyy_yyy_0,   \
                             ta_yyyyy_yyy_1,   \
                             ta_yyyyy_yyyy_0,  \
                             ta_yyyyy_yyyy_1,  \
                             ta_yyyyy_yyyz_0,  \
                             ta_yyyyy_yyyz_1,  \
                             ta_yyyyy_yyz_0,   \
                             ta_yyyyy_yyz_1,   \
                             ta_yyyyy_yyzz_0,  \
                             ta_yyyyy_yyzz_1,  \
                             ta_yyyyy_yzz_0,   \
                             ta_yyyyy_yzz_1,   \
                             ta_yyyyy_yzzz_0,  \
                             ta_yyyyy_yzzz_1,  \
                             ta_yyyyy_zzz_0,   \
                             ta_yyyyy_zzz_1,   \
                             ta_yyyyy_zzzz_0,  \
                             ta_yyyyy_zzzz_1,  \
                             ta_yyyyyy_xxxx_0, \
                             ta_yyyyyy_xxxy_0, \
                             ta_yyyyyy_xxxz_0, \
                             ta_yyyyyy_xxyy_0, \
                             ta_yyyyyy_xxyz_0, \
                             ta_yyyyyy_xxzz_0, \
                             ta_yyyyyy_xyyy_0, \
                             ta_yyyyyy_xyyz_0, \
                             ta_yyyyyy_xyzz_0, \
                             ta_yyyyyy_xzzz_0, \
                             ta_yyyyyy_yyyy_0, \
                             ta_yyyyyy_yyyz_0, \
                             ta_yyyyyy_yyzz_0, \
                             ta_yyyyyy_yzzz_0, \
                             ta_yyyyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyy_xxxx_0[i] =
            5.0 * ta_yyyy_xxxx_0[i] * fe_0 - 5.0 * ta_yyyy_xxxx_1[i] * fe_0 + ta_yyyyy_xxxx_0[i] * pa_y[i] - ta_yyyyy_xxxx_1[i] * pc_y[i];

        ta_yyyyyy_xxxy_0[i] = 5.0 * ta_yyyy_xxxy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxy_1[i] * fe_0 + ta_yyyyy_xxx_0[i] * fe_0 - ta_yyyyy_xxx_1[i] * fe_0 +
                              ta_yyyyy_xxxy_0[i] * pa_y[i] - ta_yyyyy_xxxy_1[i] * pc_y[i];

        ta_yyyyyy_xxxz_0[i] =
            5.0 * ta_yyyy_xxxz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxz_1[i] * fe_0 + ta_yyyyy_xxxz_0[i] * pa_y[i] - ta_yyyyy_xxxz_1[i] * pc_y[i];

        ta_yyyyyy_xxyy_0[i] = 5.0 * ta_yyyy_xxyy_0[i] * fe_0 - 5.0 * ta_yyyy_xxyy_1[i] * fe_0 + 2.0 * ta_yyyyy_xxy_0[i] * fe_0 -
                              2.0 * ta_yyyyy_xxy_1[i] * fe_0 + ta_yyyyy_xxyy_0[i] * pa_y[i] - ta_yyyyy_xxyy_1[i] * pc_y[i];

        ta_yyyyyy_xxyz_0[i] = 5.0 * ta_yyyy_xxyz_0[i] * fe_0 - 5.0 * ta_yyyy_xxyz_1[i] * fe_0 + ta_yyyyy_xxz_0[i] * fe_0 - ta_yyyyy_xxz_1[i] * fe_0 +
                              ta_yyyyy_xxyz_0[i] * pa_y[i] - ta_yyyyy_xxyz_1[i] * pc_y[i];

        ta_yyyyyy_xxzz_0[i] =
            5.0 * ta_yyyy_xxzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxzz_1[i] * fe_0 + ta_yyyyy_xxzz_0[i] * pa_y[i] - ta_yyyyy_xxzz_1[i] * pc_y[i];

        ta_yyyyyy_xyyy_0[i] = 5.0 * ta_yyyy_xyyy_0[i] * fe_0 - 5.0 * ta_yyyy_xyyy_1[i] * fe_0 + 3.0 * ta_yyyyy_xyy_0[i] * fe_0 -
                              3.0 * ta_yyyyy_xyy_1[i] * fe_0 + ta_yyyyy_xyyy_0[i] * pa_y[i] - ta_yyyyy_xyyy_1[i] * pc_y[i];

        ta_yyyyyy_xyyz_0[i] = 5.0 * ta_yyyy_xyyz_0[i] * fe_0 - 5.0 * ta_yyyy_xyyz_1[i] * fe_0 + 2.0 * ta_yyyyy_xyz_0[i] * fe_0 -
                              2.0 * ta_yyyyy_xyz_1[i] * fe_0 + ta_yyyyy_xyyz_0[i] * pa_y[i] - ta_yyyyy_xyyz_1[i] * pc_y[i];

        ta_yyyyyy_xyzz_0[i] = 5.0 * ta_yyyy_xyzz_0[i] * fe_0 - 5.0 * ta_yyyy_xyzz_1[i] * fe_0 + ta_yyyyy_xzz_0[i] * fe_0 - ta_yyyyy_xzz_1[i] * fe_0 +
                              ta_yyyyy_xyzz_0[i] * pa_y[i] - ta_yyyyy_xyzz_1[i] * pc_y[i];

        ta_yyyyyy_xzzz_0[i] =
            5.0 * ta_yyyy_xzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xzzz_1[i] * fe_0 + ta_yyyyy_xzzz_0[i] * pa_y[i] - ta_yyyyy_xzzz_1[i] * pc_y[i];

        ta_yyyyyy_yyyy_0[i] = 5.0 * ta_yyyy_yyyy_0[i] * fe_0 - 5.0 * ta_yyyy_yyyy_1[i] * fe_0 + 4.0 * ta_yyyyy_yyy_0[i] * fe_0 -
                              4.0 * ta_yyyyy_yyy_1[i] * fe_0 + ta_yyyyy_yyyy_0[i] * pa_y[i] - ta_yyyyy_yyyy_1[i] * pc_y[i];

        ta_yyyyyy_yyyz_0[i] = 5.0 * ta_yyyy_yyyz_0[i] * fe_0 - 5.0 * ta_yyyy_yyyz_1[i] * fe_0 + 3.0 * ta_yyyyy_yyz_0[i] * fe_0 -
                              3.0 * ta_yyyyy_yyz_1[i] * fe_0 + ta_yyyyy_yyyz_0[i] * pa_y[i] - ta_yyyyy_yyyz_1[i] * pc_y[i];

        ta_yyyyyy_yyzz_0[i] = 5.0 * ta_yyyy_yyzz_0[i] * fe_0 - 5.0 * ta_yyyy_yyzz_1[i] * fe_0 + 2.0 * ta_yyyyy_yzz_0[i] * fe_0 -
                              2.0 * ta_yyyyy_yzz_1[i] * fe_0 + ta_yyyyy_yyzz_0[i] * pa_y[i] - ta_yyyyy_yyzz_1[i] * pc_y[i];

        ta_yyyyyy_yzzz_0[i] = 5.0 * ta_yyyy_yzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yzzz_1[i] * fe_0 + ta_yyyyy_zzz_0[i] * fe_0 - ta_yyyyy_zzz_1[i] * fe_0 +
                              ta_yyyyy_yzzz_0[i] * pa_y[i] - ta_yyyyy_yzzz_1[i] * pc_y[i];

        ta_yyyyyy_zzzz_0[i] =
            5.0 * ta_yyyy_zzzz_0[i] * fe_0 - 5.0 * ta_yyyy_zzzz_1[i] * fe_0 + ta_yyyyy_zzzz_0[i] * pa_y[i] - ta_yyyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 330-345 components of targeted buffer : IG

    auto ta_yyyyyz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 330);

    auto ta_yyyyyz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 331);

    auto ta_yyyyyz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 332);

    auto ta_yyyyyz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 333);

    auto ta_yyyyyz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 334);

    auto ta_yyyyyz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 335);

    auto ta_yyyyyz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 336);

    auto ta_yyyyyz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 337);

    auto ta_yyyyyz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 338);

    auto ta_yyyyyz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 339);

    auto ta_yyyyyz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 340);

    auto ta_yyyyyz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 341);

    auto ta_yyyyyz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 342);

    auto ta_yyyyyz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 343);

    auto ta_yyyyyz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 344);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyyyy_xxxx_0,  \
                             ta_yyyyy_xxxx_1,  \
                             ta_yyyyy_xxxy_0,  \
                             ta_yyyyy_xxxy_1,  \
                             ta_yyyyy_xxy_0,   \
                             ta_yyyyy_xxy_1,   \
                             ta_yyyyy_xxyy_0,  \
                             ta_yyyyy_xxyy_1,  \
                             ta_yyyyy_xxyz_0,  \
                             ta_yyyyy_xxyz_1,  \
                             ta_yyyyy_xyy_0,   \
                             ta_yyyyy_xyy_1,   \
                             ta_yyyyy_xyyy_0,  \
                             ta_yyyyy_xyyy_1,  \
                             ta_yyyyy_xyyz_0,  \
                             ta_yyyyy_xyyz_1,  \
                             ta_yyyyy_xyz_0,   \
                             ta_yyyyy_xyz_1,   \
                             ta_yyyyy_xyzz_0,  \
                             ta_yyyyy_xyzz_1,  \
                             ta_yyyyy_yyy_0,   \
                             ta_yyyyy_yyy_1,   \
                             ta_yyyyy_yyyy_0,  \
                             ta_yyyyy_yyyy_1,  \
                             ta_yyyyy_yyyz_0,  \
                             ta_yyyyy_yyyz_1,  \
                             ta_yyyyy_yyz_0,   \
                             ta_yyyyy_yyz_1,   \
                             ta_yyyyy_yyzz_0,  \
                             ta_yyyyy_yyzz_1,  \
                             ta_yyyyy_yzz_0,   \
                             ta_yyyyy_yzz_1,   \
                             ta_yyyyy_yzzz_0,  \
                             ta_yyyyy_yzzz_1,  \
                             ta_yyyyyz_xxxx_0, \
                             ta_yyyyyz_xxxy_0, \
                             ta_yyyyyz_xxxz_0, \
                             ta_yyyyyz_xxyy_0, \
                             ta_yyyyyz_xxyz_0, \
                             ta_yyyyyz_xxzz_0, \
                             ta_yyyyyz_xyyy_0, \
                             ta_yyyyyz_xyyz_0, \
                             ta_yyyyyz_xyzz_0, \
                             ta_yyyyyz_xzzz_0, \
                             ta_yyyyyz_yyyy_0, \
                             ta_yyyyyz_yyyz_0, \
                             ta_yyyyyz_yyzz_0, \
                             ta_yyyyyz_yzzz_0, \
                             ta_yyyyyz_zzzz_0, \
                             ta_yyyyz_xxxz_0,  \
                             ta_yyyyz_xxxz_1,  \
                             ta_yyyyz_xxzz_0,  \
                             ta_yyyyz_xxzz_1,  \
                             ta_yyyyz_xzzz_0,  \
                             ta_yyyyz_xzzz_1,  \
                             ta_yyyyz_zzzz_0,  \
                             ta_yyyyz_zzzz_1,  \
                             ta_yyyz_xxxz_0,   \
                             ta_yyyz_xxxz_1,   \
                             ta_yyyz_xxzz_0,   \
                             ta_yyyz_xxzz_1,   \
                             ta_yyyz_xzzz_0,   \
                             ta_yyyz_xzzz_1,   \
                             ta_yyyz_zzzz_0,   \
                             ta_yyyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyz_xxxx_0[i] = ta_yyyyy_xxxx_0[i] * pa_z[i] - ta_yyyyy_xxxx_1[i] * pc_z[i];

        ta_yyyyyz_xxxy_0[i] = ta_yyyyy_xxxy_0[i] * pa_z[i] - ta_yyyyy_xxxy_1[i] * pc_z[i];

        ta_yyyyyz_xxxz_0[i] =
            4.0 * ta_yyyz_xxxz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxz_1[i] * fe_0 + ta_yyyyz_xxxz_0[i] * pa_y[i] - ta_yyyyz_xxxz_1[i] * pc_y[i];

        ta_yyyyyz_xxyy_0[i] = ta_yyyyy_xxyy_0[i] * pa_z[i] - ta_yyyyy_xxyy_1[i] * pc_z[i];

        ta_yyyyyz_xxyz_0[i] = ta_yyyyy_xxy_0[i] * fe_0 - ta_yyyyy_xxy_1[i] * fe_0 + ta_yyyyy_xxyz_0[i] * pa_z[i] - ta_yyyyy_xxyz_1[i] * pc_z[i];

        ta_yyyyyz_xxzz_0[i] =
            4.0 * ta_yyyz_xxzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxzz_1[i] * fe_0 + ta_yyyyz_xxzz_0[i] * pa_y[i] - ta_yyyyz_xxzz_1[i] * pc_y[i];

        ta_yyyyyz_xyyy_0[i] = ta_yyyyy_xyyy_0[i] * pa_z[i] - ta_yyyyy_xyyy_1[i] * pc_z[i];

        ta_yyyyyz_xyyz_0[i] = ta_yyyyy_xyy_0[i] * fe_0 - ta_yyyyy_xyy_1[i] * fe_0 + ta_yyyyy_xyyz_0[i] * pa_z[i] - ta_yyyyy_xyyz_1[i] * pc_z[i];

        ta_yyyyyz_xyzz_0[i] =
            2.0 * ta_yyyyy_xyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyz_1[i] * fe_0 + ta_yyyyy_xyzz_0[i] * pa_z[i] - ta_yyyyy_xyzz_1[i] * pc_z[i];

        ta_yyyyyz_xzzz_0[i] =
            4.0 * ta_yyyz_xzzz_0[i] * fe_0 - 4.0 * ta_yyyz_xzzz_1[i] * fe_0 + ta_yyyyz_xzzz_0[i] * pa_y[i] - ta_yyyyz_xzzz_1[i] * pc_y[i];

        ta_yyyyyz_yyyy_0[i] = ta_yyyyy_yyyy_0[i] * pa_z[i] - ta_yyyyy_yyyy_1[i] * pc_z[i];

        ta_yyyyyz_yyyz_0[i] = ta_yyyyy_yyy_0[i] * fe_0 - ta_yyyyy_yyy_1[i] * fe_0 + ta_yyyyy_yyyz_0[i] * pa_z[i] - ta_yyyyy_yyyz_1[i] * pc_z[i];

        ta_yyyyyz_yyzz_0[i] =
            2.0 * ta_yyyyy_yyz_0[i] * fe_0 - 2.0 * ta_yyyyy_yyz_1[i] * fe_0 + ta_yyyyy_yyzz_0[i] * pa_z[i] - ta_yyyyy_yyzz_1[i] * pc_z[i];

        ta_yyyyyz_yzzz_0[i] =
            3.0 * ta_yyyyy_yzz_0[i] * fe_0 - 3.0 * ta_yyyyy_yzz_1[i] * fe_0 + ta_yyyyy_yzzz_0[i] * pa_z[i] - ta_yyyyy_yzzz_1[i] * pc_z[i];

        ta_yyyyyz_zzzz_0[i] =
            4.0 * ta_yyyz_zzzz_0[i] * fe_0 - 4.0 * ta_yyyz_zzzz_1[i] * fe_0 + ta_yyyyz_zzzz_0[i] * pa_y[i] - ta_yyyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 345-360 components of targeted buffer : IG

    auto ta_yyyyzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 345);

    auto ta_yyyyzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 346);

    auto ta_yyyyzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 347);

    auto ta_yyyyzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 348);

    auto ta_yyyyzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 349);

    auto ta_yyyyzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 350);

    auto ta_yyyyzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 351);

    auto ta_yyyyzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 352);

    auto ta_yyyyzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 353);

    auto ta_yyyyzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 354);

    auto ta_yyyyzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 355);

    auto ta_yyyyzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 356);

    auto ta_yyyyzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 357);

    auto ta_yyyyzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 358);

    auto ta_yyyyzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 359);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyyy_xxxy_0,   \
                             ta_yyyy_xxxy_1,   \
                             ta_yyyy_xxyy_0,   \
                             ta_yyyy_xxyy_1,   \
                             ta_yyyy_xyyy_0,   \
                             ta_yyyy_xyyy_1,   \
                             ta_yyyy_yyyy_0,   \
                             ta_yyyy_yyyy_1,   \
                             ta_yyyyz_xxxy_0,  \
                             ta_yyyyz_xxxy_1,  \
                             ta_yyyyz_xxyy_0,  \
                             ta_yyyyz_xxyy_1,  \
                             ta_yyyyz_xyyy_0,  \
                             ta_yyyyz_xyyy_1,  \
                             ta_yyyyz_yyyy_0,  \
                             ta_yyyyz_yyyy_1,  \
                             ta_yyyyzz_xxxx_0, \
                             ta_yyyyzz_xxxy_0, \
                             ta_yyyyzz_xxxz_0, \
                             ta_yyyyzz_xxyy_0, \
                             ta_yyyyzz_xxyz_0, \
                             ta_yyyyzz_xxzz_0, \
                             ta_yyyyzz_xyyy_0, \
                             ta_yyyyzz_xyyz_0, \
                             ta_yyyyzz_xyzz_0, \
                             ta_yyyyzz_xzzz_0, \
                             ta_yyyyzz_yyyy_0, \
                             ta_yyyyzz_yyyz_0, \
                             ta_yyyyzz_yyzz_0, \
                             ta_yyyyzz_yzzz_0, \
                             ta_yyyyzz_zzzz_0, \
                             ta_yyyzz_xxxx_0,  \
                             ta_yyyzz_xxxx_1,  \
                             ta_yyyzz_xxxz_0,  \
                             ta_yyyzz_xxxz_1,  \
                             ta_yyyzz_xxyz_0,  \
                             ta_yyyzz_xxyz_1,  \
                             ta_yyyzz_xxz_0,   \
                             ta_yyyzz_xxz_1,   \
                             ta_yyyzz_xxzz_0,  \
                             ta_yyyzz_xxzz_1,  \
                             ta_yyyzz_xyyz_0,  \
                             ta_yyyzz_xyyz_1,  \
                             ta_yyyzz_xyz_0,   \
                             ta_yyyzz_xyz_1,   \
                             ta_yyyzz_xyzz_0,  \
                             ta_yyyzz_xyzz_1,  \
                             ta_yyyzz_xzz_0,   \
                             ta_yyyzz_xzz_1,   \
                             ta_yyyzz_xzzz_0,  \
                             ta_yyyzz_xzzz_1,  \
                             ta_yyyzz_yyyz_0,  \
                             ta_yyyzz_yyyz_1,  \
                             ta_yyyzz_yyz_0,   \
                             ta_yyyzz_yyz_1,   \
                             ta_yyyzz_yyzz_0,  \
                             ta_yyyzz_yyzz_1,  \
                             ta_yyyzz_yzz_0,   \
                             ta_yyyzz_yzz_1,   \
                             ta_yyyzz_yzzz_0,  \
                             ta_yyyzz_yzzz_1,  \
                             ta_yyyzz_zzz_0,   \
                             ta_yyyzz_zzz_1,   \
                             ta_yyyzz_zzzz_0,  \
                             ta_yyyzz_zzzz_1,  \
                             ta_yyzz_xxxx_0,   \
                             ta_yyzz_xxxx_1,   \
                             ta_yyzz_xxxz_0,   \
                             ta_yyzz_xxxz_1,   \
                             ta_yyzz_xxyz_0,   \
                             ta_yyzz_xxyz_1,   \
                             ta_yyzz_xxzz_0,   \
                             ta_yyzz_xxzz_1,   \
                             ta_yyzz_xyyz_0,   \
                             ta_yyzz_xyyz_1,   \
                             ta_yyzz_xyzz_0,   \
                             ta_yyzz_xyzz_1,   \
                             ta_yyzz_xzzz_0,   \
                             ta_yyzz_xzzz_1,   \
                             ta_yyzz_yyyz_0,   \
                             ta_yyzz_yyyz_1,   \
                             ta_yyzz_yyzz_0,   \
                             ta_yyzz_yyzz_1,   \
                             ta_yyzz_yzzz_0,   \
                             ta_yyzz_yzzz_1,   \
                             ta_yyzz_zzzz_0,   \
                             ta_yyzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyzz_xxxx_0[i] =
            3.0 * ta_yyzz_xxxx_0[i] * fe_0 - 3.0 * ta_yyzz_xxxx_1[i] * fe_0 + ta_yyyzz_xxxx_0[i] * pa_y[i] - ta_yyyzz_xxxx_1[i] * pc_y[i];

        ta_yyyyzz_xxxy_0[i] = ta_yyyy_xxxy_0[i] * fe_0 - ta_yyyy_xxxy_1[i] * fe_0 + ta_yyyyz_xxxy_0[i] * pa_z[i] - ta_yyyyz_xxxy_1[i] * pc_z[i];

        ta_yyyyzz_xxxz_0[i] =
            3.0 * ta_yyzz_xxxz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxz_1[i] * fe_0 + ta_yyyzz_xxxz_0[i] * pa_y[i] - ta_yyyzz_xxxz_1[i] * pc_y[i];

        ta_yyyyzz_xxyy_0[i] = ta_yyyy_xxyy_0[i] * fe_0 - ta_yyyy_xxyy_1[i] * fe_0 + ta_yyyyz_xxyy_0[i] * pa_z[i] - ta_yyyyz_xxyy_1[i] * pc_z[i];

        ta_yyyyzz_xxyz_0[i] = 3.0 * ta_yyzz_xxyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyz_1[i] * fe_0 + ta_yyyzz_xxz_0[i] * fe_0 - ta_yyyzz_xxz_1[i] * fe_0 +
                              ta_yyyzz_xxyz_0[i] * pa_y[i] - ta_yyyzz_xxyz_1[i] * pc_y[i];

        ta_yyyyzz_xxzz_0[i] =
            3.0 * ta_yyzz_xxzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxzz_1[i] * fe_0 + ta_yyyzz_xxzz_0[i] * pa_y[i] - ta_yyyzz_xxzz_1[i] * pc_y[i];

        ta_yyyyzz_xyyy_0[i] = ta_yyyy_xyyy_0[i] * fe_0 - ta_yyyy_xyyy_1[i] * fe_0 + ta_yyyyz_xyyy_0[i] * pa_z[i] - ta_yyyyz_xyyy_1[i] * pc_z[i];

        ta_yyyyzz_xyyz_0[i] = 3.0 * ta_yyzz_xyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xyyz_1[i] * fe_0 + 2.0 * ta_yyyzz_xyz_0[i] * fe_0 -
                              2.0 * ta_yyyzz_xyz_1[i] * fe_0 + ta_yyyzz_xyyz_0[i] * pa_y[i] - ta_yyyzz_xyyz_1[i] * pc_y[i];

        ta_yyyyzz_xyzz_0[i] = 3.0 * ta_yyzz_xyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xyzz_1[i] * fe_0 + ta_yyyzz_xzz_0[i] * fe_0 - ta_yyyzz_xzz_1[i] * fe_0 +
                              ta_yyyzz_xyzz_0[i] * pa_y[i] - ta_yyyzz_xyzz_1[i] * pc_y[i];

        ta_yyyyzz_xzzz_0[i] =
            3.0 * ta_yyzz_xzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xzzz_1[i] * fe_0 + ta_yyyzz_xzzz_0[i] * pa_y[i] - ta_yyyzz_xzzz_1[i] * pc_y[i];

        ta_yyyyzz_yyyy_0[i] = ta_yyyy_yyyy_0[i] * fe_0 - ta_yyyy_yyyy_1[i] * fe_0 + ta_yyyyz_yyyy_0[i] * pa_z[i] - ta_yyyyz_yyyy_1[i] * pc_z[i];

        ta_yyyyzz_yyyz_0[i] = 3.0 * ta_yyzz_yyyz_0[i] * fe_0 - 3.0 * ta_yyzz_yyyz_1[i] * fe_0 + 3.0 * ta_yyyzz_yyz_0[i] * fe_0 -
                              3.0 * ta_yyyzz_yyz_1[i] * fe_0 + ta_yyyzz_yyyz_0[i] * pa_y[i] - ta_yyyzz_yyyz_1[i] * pc_y[i];

        ta_yyyyzz_yyzz_0[i] = 3.0 * ta_yyzz_yyzz_0[i] * fe_0 - 3.0 * ta_yyzz_yyzz_1[i] * fe_0 + 2.0 * ta_yyyzz_yzz_0[i] * fe_0 -
                              2.0 * ta_yyyzz_yzz_1[i] * fe_0 + ta_yyyzz_yyzz_0[i] * pa_y[i] - ta_yyyzz_yyzz_1[i] * pc_y[i];

        ta_yyyyzz_yzzz_0[i] = 3.0 * ta_yyzz_yzzz_0[i] * fe_0 - 3.0 * ta_yyzz_yzzz_1[i] * fe_0 + ta_yyyzz_zzz_0[i] * fe_0 - ta_yyyzz_zzz_1[i] * fe_0 +
                              ta_yyyzz_yzzz_0[i] * pa_y[i] - ta_yyyzz_yzzz_1[i] * pc_y[i];

        ta_yyyyzz_zzzz_0[i] =
            3.0 * ta_yyzz_zzzz_0[i] * fe_0 - 3.0 * ta_yyzz_zzzz_1[i] * fe_0 + ta_yyyzz_zzzz_0[i] * pa_y[i] - ta_yyyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 360-375 components of targeted buffer : IG

    auto ta_yyyzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 360);

    auto ta_yyyzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 361);

    auto ta_yyyzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 362);

    auto ta_yyyzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 363);

    auto ta_yyyzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 364);

    auto ta_yyyzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 365);

    auto ta_yyyzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 366);

    auto ta_yyyzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 367);

    auto ta_yyyzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 368);

    auto ta_yyyzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 369);

    auto ta_yyyzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 370);

    auto ta_yyyzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 371);

    auto ta_yyyzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 372);

    auto ta_yyyzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 373);

    auto ta_yyyzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 374);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyyz_xxxy_0,   \
                             ta_yyyz_xxxy_1,   \
                             ta_yyyz_xxyy_0,   \
                             ta_yyyz_xxyy_1,   \
                             ta_yyyz_xyyy_0,   \
                             ta_yyyz_xyyy_1,   \
                             ta_yyyz_yyyy_0,   \
                             ta_yyyz_yyyy_1,   \
                             ta_yyyzz_xxxy_0,  \
                             ta_yyyzz_xxxy_1,  \
                             ta_yyyzz_xxyy_0,  \
                             ta_yyyzz_xxyy_1,  \
                             ta_yyyzz_xyyy_0,  \
                             ta_yyyzz_xyyy_1,  \
                             ta_yyyzz_yyyy_0,  \
                             ta_yyyzz_yyyy_1,  \
                             ta_yyyzzz_xxxx_0, \
                             ta_yyyzzz_xxxy_0, \
                             ta_yyyzzz_xxxz_0, \
                             ta_yyyzzz_xxyy_0, \
                             ta_yyyzzz_xxyz_0, \
                             ta_yyyzzz_xxzz_0, \
                             ta_yyyzzz_xyyy_0, \
                             ta_yyyzzz_xyyz_0, \
                             ta_yyyzzz_xyzz_0, \
                             ta_yyyzzz_xzzz_0, \
                             ta_yyyzzz_yyyy_0, \
                             ta_yyyzzz_yyyz_0, \
                             ta_yyyzzz_yyzz_0, \
                             ta_yyyzzz_yzzz_0, \
                             ta_yyyzzz_zzzz_0, \
                             ta_yyzzz_xxxx_0,  \
                             ta_yyzzz_xxxx_1,  \
                             ta_yyzzz_xxxz_0,  \
                             ta_yyzzz_xxxz_1,  \
                             ta_yyzzz_xxyz_0,  \
                             ta_yyzzz_xxyz_1,  \
                             ta_yyzzz_xxz_0,   \
                             ta_yyzzz_xxz_1,   \
                             ta_yyzzz_xxzz_0,  \
                             ta_yyzzz_xxzz_1,  \
                             ta_yyzzz_xyyz_0,  \
                             ta_yyzzz_xyyz_1,  \
                             ta_yyzzz_xyz_0,   \
                             ta_yyzzz_xyz_1,   \
                             ta_yyzzz_xyzz_0,  \
                             ta_yyzzz_xyzz_1,  \
                             ta_yyzzz_xzz_0,   \
                             ta_yyzzz_xzz_1,   \
                             ta_yyzzz_xzzz_0,  \
                             ta_yyzzz_xzzz_1,  \
                             ta_yyzzz_yyyz_0,  \
                             ta_yyzzz_yyyz_1,  \
                             ta_yyzzz_yyz_0,   \
                             ta_yyzzz_yyz_1,   \
                             ta_yyzzz_yyzz_0,  \
                             ta_yyzzz_yyzz_1,  \
                             ta_yyzzz_yzz_0,   \
                             ta_yyzzz_yzz_1,   \
                             ta_yyzzz_yzzz_0,  \
                             ta_yyzzz_yzzz_1,  \
                             ta_yyzzz_zzz_0,   \
                             ta_yyzzz_zzz_1,   \
                             ta_yyzzz_zzzz_0,  \
                             ta_yyzzz_zzzz_1,  \
                             ta_yzzz_xxxx_0,   \
                             ta_yzzz_xxxx_1,   \
                             ta_yzzz_xxxz_0,   \
                             ta_yzzz_xxxz_1,   \
                             ta_yzzz_xxyz_0,   \
                             ta_yzzz_xxyz_1,   \
                             ta_yzzz_xxzz_0,   \
                             ta_yzzz_xxzz_1,   \
                             ta_yzzz_xyyz_0,   \
                             ta_yzzz_xyyz_1,   \
                             ta_yzzz_xyzz_0,   \
                             ta_yzzz_xyzz_1,   \
                             ta_yzzz_xzzz_0,   \
                             ta_yzzz_xzzz_1,   \
                             ta_yzzz_yyyz_0,   \
                             ta_yzzz_yyyz_1,   \
                             ta_yzzz_yyzz_0,   \
                             ta_yzzz_yyzz_1,   \
                             ta_yzzz_yzzz_0,   \
                             ta_yzzz_yzzz_1,   \
                             ta_yzzz_zzzz_0,   \
                             ta_yzzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzzz_xxxx_0[i] =
            2.0 * ta_yzzz_xxxx_0[i] * fe_0 - 2.0 * ta_yzzz_xxxx_1[i] * fe_0 + ta_yyzzz_xxxx_0[i] * pa_y[i] - ta_yyzzz_xxxx_1[i] * pc_y[i];

        ta_yyyzzz_xxxy_0[i] =
            2.0 * ta_yyyz_xxxy_0[i] * fe_0 - 2.0 * ta_yyyz_xxxy_1[i] * fe_0 + ta_yyyzz_xxxy_0[i] * pa_z[i] - ta_yyyzz_xxxy_1[i] * pc_z[i];

        ta_yyyzzz_xxxz_0[i] =
            2.0 * ta_yzzz_xxxz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxz_1[i] * fe_0 + ta_yyzzz_xxxz_0[i] * pa_y[i] - ta_yyzzz_xxxz_1[i] * pc_y[i];

        ta_yyyzzz_xxyy_0[i] =
            2.0 * ta_yyyz_xxyy_0[i] * fe_0 - 2.0 * ta_yyyz_xxyy_1[i] * fe_0 + ta_yyyzz_xxyy_0[i] * pa_z[i] - ta_yyyzz_xxyy_1[i] * pc_z[i];

        ta_yyyzzz_xxyz_0[i] = 2.0 * ta_yzzz_xxyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyz_1[i] * fe_0 + ta_yyzzz_xxz_0[i] * fe_0 - ta_yyzzz_xxz_1[i] * fe_0 +
                              ta_yyzzz_xxyz_0[i] * pa_y[i] - ta_yyzzz_xxyz_1[i] * pc_y[i];

        ta_yyyzzz_xxzz_0[i] =
            2.0 * ta_yzzz_xxzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxzz_1[i] * fe_0 + ta_yyzzz_xxzz_0[i] * pa_y[i] - ta_yyzzz_xxzz_1[i] * pc_y[i];

        ta_yyyzzz_xyyy_0[i] =
            2.0 * ta_yyyz_xyyy_0[i] * fe_0 - 2.0 * ta_yyyz_xyyy_1[i] * fe_0 + ta_yyyzz_xyyy_0[i] * pa_z[i] - ta_yyyzz_xyyy_1[i] * pc_z[i];

        ta_yyyzzz_xyyz_0[i] = 2.0 * ta_yzzz_xyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyz_1[i] * fe_0 + 2.0 * ta_yyzzz_xyz_0[i] * fe_0 -
                              2.0 * ta_yyzzz_xyz_1[i] * fe_0 + ta_yyzzz_xyyz_0[i] * pa_y[i] - ta_yyzzz_xyyz_1[i] * pc_y[i];

        ta_yyyzzz_xyzz_0[i] = 2.0 * ta_yzzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyzz_1[i] * fe_0 + ta_yyzzz_xzz_0[i] * fe_0 - ta_yyzzz_xzz_1[i] * fe_0 +
                              ta_yyzzz_xyzz_0[i] * pa_y[i] - ta_yyzzz_xyzz_1[i] * pc_y[i];

        ta_yyyzzz_xzzz_0[i] =
            2.0 * ta_yzzz_xzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xzzz_1[i] * fe_0 + ta_yyzzz_xzzz_0[i] * pa_y[i] - ta_yyzzz_xzzz_1[i] * pc_y[i];

        ta_yyyzzz_yyyy_0[i] =
            2.0 * ta_yyyz_yyyy_0[i] * fe_0 - 2.0 * ta_yyyz_yyyy_1[i] * fe_0 + ta_yyyzz_yyyy_0[i] * pa_z[i] - ta_yyyzz_yyyy_1[i] * pc_z[i];

        ta_yyyzzz_yyyz_0[i] = 2.0 * ta_yzzz_yyyz_0[i] * fe_0 - 2.0 * ta_yzzz_yyyz_1[i] * fe_0 + 3.0 * ta_yyzzz_yyz_0[i] * fe_0 -
                              3.0 * ta_yyzzz_yyz_1[i] * fe_0 + ta_yyzzz_yyyz_0[i] * pa_y[i] - ta_yyzzz_yyyz_1[i] * pc_y[i];

        ta_yyyzzz_yyzz_0[i] = 2.0 * ta_yzzz_yyzz_0[i] * fe_0 - 2.0 * ta_yzzz_yyzz_1[i] * fe_0 + 2.0 * ta_yyzzz_yzz_0[i] * fe_0 -
                              2.0 * ta_yyzzz_yzz_1[i] * fe_0 + ta_yyzzz_yyzz_0[i] * pa_y[i] - ta_yyzzz_yyzz_1[i] * pc_y[i];

        ta_yyyzzz_yzzz_0[i] = 2.0 * ta_yzzz_yzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yzzz_1[i] * fe_0 + ta_yyzzz_zzz_0[i] * fe_0 - ta_yyzzz_zzz_1[i] * fe_0 +
                              ta_yyzzz_yzzz_0[i] * pa_y[i] - ta_yyzzz_yzzz_1[i] * pc_y[i];

        ta_yyyzzz_zzzz_0[i] =
            2.0 * ta_yzzz_zzzz_0[i] * fe_0 - 2.0 * ta_yzzz_zzzz_1[i] * fe_0 + ta_yyzzz_zzzz_0[i] * pa_y[i] - ta_yyzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 375-390 components of targeted buffer : IG

    auto ta_yyzzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 375);

    auto ta_yyzzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 376);

    auto ta_yyzzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 377);

    auto ta_yyzzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 378);

    auto ta_yyzzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 379);

    auto ta_yyzzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 380);

    auto ta_yyzzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 381);

    auto ta_yyzzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 382);

    auto ta_yyzzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 383);

    auto ta_yyzzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 384);

    auto ta_yyzzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 385);

    auto ta_yyzzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 386);

    auto ta_yyzzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 387);

    auto ta_yyzzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 388);

    auto ta_yyzzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 389);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyzz_xxxy_0,   \
                             ta_yyzz_xxxy_1,   \
                             ta_yyzz_xxyy_0,   \
                             ta_yyzz_xxyy_1,   \
                             ta_yyzz_xyyy_0,   \
                             ta_yyzz_xyyy_1,   \
                             ta_yyzz_yyyy_0,   \
                             ta_yyzz_yyyy_1,   \
                             ta_yyzzz_xxxy_0,  \
                             ta_yyzzz_xxxy_1,  \
                             ta_yyzzz_xxyy_0,  \
                             ta_yyzzz_xxyy_1,  \
                             ta_yyzzz_xyyy_0,  \
                             ta_yyzzz_xyyy_1,  \
                             ta_yyzzz_yyyy_0,  \
                             ta_yyzzz_yyyy_1,  \
                             ta_yyzzzz_xxxx_0, \
                             ta_yyzzzz_xxxy_0, \
                             ta_yyzzzz_xxxz_0, \
                             ta_yyzzzz_xxyy_0, \
                             ta_yyzzzz_xxyz_0, \
                             ta_yyzzzz_xxzz_0, \
                             ta_yyzzzz_xyyy_0, \
                             ta_yyzzzz_xyyz_0, \
                             ta_yyzzzz_xyzz_0, \
                             ta_yyzzzz_xzzz_0, \
                             ta_yyzzzz_yyyy_0, \
                             ta_yyzzzz_yyyz_0, \
                             ta_yyzzzz_yyzz_0, \
                             ta_yyzzzz_yzzz_0, \
                             ta_yyzzzz_zzzz_0, \
                             ta_yzzzz_xxxx_0,  \
                             ta_yzzzz_xxxx_1,  \
                             ta_yzzzz_xxxz_0,  \
                             ta_yzzzz_xxxz_1,  \
                             ta_yzzzz_xxyz_0,  \
                             ta_yzzzz_xxyz_1,  \
                             ta_yzzzz_xxz_0,   \
                             ta_yzzzz_xxz_1,   \
                             ta_yzzzz_xxzz_0,  \
                             ta_yzzzz_xxzz_1,  \
                             ta_yzzzz_xyyz_0,  \
                             ta_yzzzz_xyyz_1,  \
                             ta_yzzzz_xyz_0,   \
                             ta_yzzzz_xyz_1,   \
                             ta_yzzzz_xyzz_0,  \
                             ta_yzzzz_xyzz_1,  \
                             ta_yzzzz_xzz_0,   \
                             ta_yzzzz_xzz_1,   \
                             ta_yzzzz_xzzz_0,  \
                             ta_yzzzz_xzzz_1,  \
                             ta_yzzzz_yyyz_0,  \
                             ta_yzzzz_yyyz_1,  \
                             ta_yzzzz_yyz_0,   \
                             ta_yzzzz_yyz_1,   \
                             ta_yzzzz_yyzz_0,  \
                             ta_yzzzz_yyzz_1,  \
                             ta_yzzzz_yzz_0,   \
                             ta_yzzzz_yzz_1,   \
                             ta_yzzzz_yzzz_0,  \
                             ta_yzzzz_yzzz_1,  \
                             ta_yzzzz_zzz_0,   \
                             ta_yzzzz_zzz_1,   \
                             ta_yzzzz_zzzz_0,  \
                             ta_yzzzz_zzzz_1,  \
                             ta_zzzz_xxxx_0,   \
                             ta_zzzz_xxxx_1,   \
                             ta_zzzz_xxxz_0,   \
                             ta_zzzz_xxxz_1,   \
                             ta_zzzz_xxyz_0,   \
                             ta_zzzz_xxyz_1,   \
                             ta_zzzz_xxzz_0,   \
                             ta_zzzz_xxzz_1,   \
                             ta_zzzz_xyyz_0,   \
                             ta_zzzz_xyyz_1,   \
                             ta_zzzz_xyzz_0,   \
                             ta_zzzz_xyzz_1,   \
                             ta_zzzz_xzzz_0,   \
                             ta_zzzz_xzzz_1,   \
                             ta_zzzz_yyyz_0,   \
                             ta_zzzz_yyyz_1,   \
                             ta_zzzz_yyzz_0,   \
                             ta_zzzz_yyzz_1,   \
                             ta_zzzz_yzzz_0,   \
                             ta_zzzz_yzzz_1,   \
                             ta_zzzz_zzzz_0,   \
                             ta_zzzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzzz_xxxx_0[i] = ta_zzzz_xxxx_0[i] * fe_0 - ta_zzzz_xxxx_1[i] * fe_0 + ta_yzzzz_xxxx_0[i] * pa_y[i] - ta_yzzzz_xxxx_1[i] * pc_y[i];

        ta_yyzzzz_xxxy_0[i] =
            3.0 * ta_yyzz_xxxy_0[i] * fe_0 - 3.0 * ta_yyzz_xxxy_1[i] * fe_0 + ta_yyzzz_xxxy_0[i] * pa_z[i] - ta_yyzzz_xxxy_1[i] * pc_z[i];

        ta_yyzzzz_xxxz_0[i] = ta_zzzz_xxxz_0[i] * fe_0 - ta_zzzz_xxxz_1[i] * fe_0 + ta_yzzzz_xxxz_0[i] * pa_y[i] - ta_yzzzz_xxxz_1[i] * pc_y[i];

        ta_yyzzzz_xxyy_0[i] =
            3.0 * ta_yyzz_xxyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxyy_1[i] * fe_0 + ta_yyzzz_xxyy_0[i] * pa_z[i] - ta_yyzzz_xxyy_1[i] * pc_z[i];

        ta_yyzzzz_xxyz_0[i] = ta_zzzz_xxyz_0[i] * fe_0 - ta_zzzz_xxyz_1[i] * fe_0 + ta_yzzzz_xxz_0[i] * fe_0 - ta_yzzzz_xxz_1[i] * fe_0 +
                              ta_yzzzz_xxyz_0[i] * pa_y[i] - ta_yzzzz_xxyz_1[i] * pc_y[i];

        ta_yyzzzz_xxzz_0[i] = ta_zzzz_xxzz_0[i] * fe_0 - ta_zzzz_xxzz_1[i] * fe_0 + ta_yzzzz_xxzz_0[i] * pa_y[i] - ta_yzzzz_xxzz_1[i] * pc_y[i];

        ta_yyzzzz_xyyy_0[i] =
            3.0 * ta_yyzz_xyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xyyy_1[i] * fe_0 + ta_yyzzz_xyyy_0[i] * pa_z[i] - ta_yyzzz_xyyy_1[i] * pc_z[i];

        ta_yyzzzz_xyyz_0[i] = ta_zzzz_xyyz_0[i] * fe_0 - ta_zzzz_xyyz_1[i] * fe_0 + 2.0 * ta_yzzzz_xyz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyz_1[i] * fe_0 +
                              ta_yzzzz_xyyz_0[i] * pa_y[i] - ta_yzzzz_xyyz_1[i] * pc_y[i];

        ta_yyzzzz_xyzz_0[i] = ta_zzzz_xyzz_0[i] * fe_0 - ta_zzzz_xyzz_1[i] * fe_0 + ta_yzzzz_xzz_0[i] * fe_0 - ta_yzzzz_xzz_1[i] * fe_0 +
                              ta_yzzzz_xyzz_0[i] * pa_y[i] - ta_yzzzz_xyzz_1[i] * pc_y[i];

        ta_yyzzzz_xzzz_0[i] = ta_zzzz_xzzz_0[i] * fe_0 - ta_zzzz_xzzz_1[i] * fe_0 + ta_yzzzz_xzzz_0[i] * pa_y[i] - ta_yzzzz_xzzz_1[i] * pc_y[i];

        ta_yyzzzz_yyyy_0[i] =
            3.0 * ta_yyzz_yyyy_0[i] * fe_0 - 3.0 * ta_yyzz_yyyy_1[i] * fe_0 + ta_yyzzz_yyyy_0[i] * pa_z[i] - ta_yyzzz_yyyy_1[i] * pc_z[i];

        ta_yyzzzz_yyyz_0[i] = ta_zzzz_yyyz_0[i] * fe_0 - ta_zzzz_yyyz_1[i] * fe_0 + 3.0 * ta_yzzzz_yyz_0[i] * fe_0 - 3.0 * ta_yzzzz_yyz_1[i] * fe_0 +
                              ta_yzzzz_yyyz_0[i] * pa_y[i] - ta_yzzzz_yyyz_1[i] * pc_y[i];

        ta_yyzzzz_yyzz_0[i] = ta_zzzz_yyzz_0[i] * fe_0 - ta_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta_yzzzz_yzz_0[i] * fe_0 - 2.0 * ta_yzzzz_yzz_1[i] * fe_0 +
                              ta_yzzzz_yyzz_0[i] * pa_y[i] - ta_yzzzz_yyzz_1[i] * pc_y[i];

        ta_yyzzzz_yzzz_0[i] = ta_zzzz_yzzz_0[i] * fe_0 - ta_zzzz_yzzz_1[i] * fe_0 + ta_yzzzz_zzz_0[i] * fe_0 - ta_yzzzz_zzz_1[i] * fe_0 +
                              ta_yzzzz_yzzz_0[i] * pa_y[i] - ta_yzzzz_yzzz_1[i] * pc_y[i];

        ta_yyzzzz_zzzz_0[i] = ta_zzzz_zzzz_0[i] * fe_0 - ta_zzzz_zzzz_1[i] * fe_0 + ta_yzzzz_zzzz_0[i] * pa_y[i] - ta_yzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 390-405 components of targeted buffer : IG

    auto ta_yzzzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 390);

    auto ta_yzzzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 391);

    auto ta_yzzzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 392);

    auto ta_yzzzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 393);

    auto ta_yzzzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 394);

    auto ta_yzzzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 395);

    auto ta_yzzzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 396);

    auto ta_yzzzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 397);

    auto ta_yzzzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 398);

    auto ta_yzzzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 399);

    auto ta_yzzzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 400);

    auto ta_yzzzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 401);

    auto ta_yzzzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 402);

    auto ta_yzzzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 403);

    auto ta_yzzzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 404);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta_yzzzzz_xxxx_0, \
                             ta_yzzzzz_xxxy_0, \
                             ta_yzzzzz_xxxz_0, \
                             ta_yzzzzz_xxyy_0, \
                             ta_yzzzzz_xxyz_0, \
                             ta_yzzzzz_xxzz_0, \
                             ta_yzzzzz_xyyy_0, \
                             ta_yzzzzz_xyyz_0, \
                             ta_yzzzzz_xyzz_0, \
                             ta_yzzzzz_xzzz_0, \
                             ta_yzzzzz_yyyy_0, \
                             ta_yzzzzz_yyyz_0, \
                             ta_yzzzzz_yyzz_0, \
                             ta_yzzzzz_yzzz_0, \
                             ta_yzzzzz_zzzz_0, \
                             ta_zzzzz_xxx_0,   \
                             ta_zzzzz_xxx_1,   \
                             ta_zzzzz_xxxx_0,  \
                             ta_zzzzz_xxxx_1,  \
                             ta_zzzzz_xxxy_0,  \
                             ta_zzzzz_xxxy_1,  \
                             ta_zzzzz_xxxz_0,  \
                             ta_zzzzz_xxxz_1,  \
                             ta_zzzzz_xxy_0,   \
                             ta_zzzzz_xxy_1,   \
                             ta_zzzzz_xxyy_0,  \
                             ta_zzzzz_xxyy_1,  \
                             ta_zzzzz_xxyz_0,  \
                             ta_zzzzz_xxyz_1,  \
                             ta_zzzzz_xxz_0,   \
                             ta_zzzzz_xxz_1,   \
                             ta_zzzzz_xxzz_0,  \
                             ta_zzzzz_xxzz_1,  \
                             ta_zzzzz_xyy_0,   \
                             ta_zzzzz_xyy_1,   \
                             ta_zzzzz_xyyy_0,  \
                             ta_zzzzz_xyyy_1,  \
                             ta_zzzzz_xyyz_0,  \
                             ta_zzzzz_xyyz_1,  \
                             ta_zzzzz_xyz_0,   \
                             ta_zzzzz_xyz_1,   \
                             ta_zzzzz_xyzz_0,  \
                             ta_zzzzz_xyzz_1,  \
                             ta_zzzzz_xzz_0,   \
                             ta_zzzzz_xzz_1,   \
                             ta_zzzzz_xzzz_0,  \
                             ta_zzzzz_xzzz_1,  \
                             ta_zzzzz_yyy_0,   \
                             ta_zzzzz_yyy_1,   \
                             ta_zzzzz_yyyy_0,  \
                             ta_zzzzz_yyyy_1,  \
                             ta_zzzzz_yyyz_0,  \
                             ta_zzzzz_yyyz_1,  \
                             ta_zzzzz_yyz_0,   \
                             ta_zzzzz_yyz_1,   \
                             ta_zzzzz_yyzz_0,  \
                             ta_zzzzz_yyzz_1,  \
                             ta_zzzzz_yzz_0,   \
                             ta_zzzzz_yzz_1,   \
                             ta_zzzzz_yzzz_0,  \
                             ta_zzzzz_yzzz_1,  \
                             ta_zzzzz_zzz_0,   \
                             ta_zzzzz_zzz_1,   \
                             ta_zzzzz_zzzz_0,  \
                             ta_zzzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzzz_xxxx_0[i] = ta_zzzzz_xxxx_0[i] * pa_y[i] - ta_zzzzz_xxxx_1[i] * pc_y[i];

        ta_yzzzzz_xxxy_0[i] = ta_zzzzz_xxx_0[i] * fe_0 - ta_zzzzz_xxx_1[i] * fe_0 + ta_zzzzz_xxxy_0[i] * pa_y[i] - ta_zzzzz_xxxy_1[i] * pc_y[i];

        ta_yzzzzz_xxxz_0[i] = ta_zzzzz_xxxz_0[i] * pa_y[i] - ta_zzzzz_xxxz_1[i] * pc_y[i];

        ta_yzzzzz_xxyy_0[i] =
            2.0 * ta_zzzzz_xxy_0[i] * fe_0 - 2.0 * ta_zzzzz_xxy_1[i] * fe_0 + ta_zzzzz_xxyy_0[i] * pa_y[i] - ta_zzzzz_xxyy_1[i] * pc_y[i];

        ta_yzzzzz_xxyz_0[i] = ta_zzzzz_xxz_0[i] * fe_0 - ta_zzzzz_xxz_1[i] * fe_0 + ta_zzzzz_xxyz_0[i] * pa_y[i] - ta_zzzzz_xxyz_1[i] * pc_y[i];

        ta_yzzzzz_xxzz_0[i] = ta_zzzzz_xxzz_0[i] * pa_y[i] - ta_zzzzz_xxzz_1[i] * pc_y[i];

        ta_yzzzzz_xyyy_0[i] =
            3.0 * ta_zzzzz_xyy_0[i] * fe_0 - 3.0 * ta_zzzzz_xyy_1[i] * fe_0 + ta_zzzzz_xyyy_0[i] * pa_y[i] - ta_zzzzz_xyyy_1[i] * pc_y[i];

        ta_yzzzzz_xyyz_0[i] =
            2.0 * ta_zzzzz_xyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyz_1[i] * fe_0 + ta_zzzzz_xyyz_0[i] * pa_y[i] - ta_zzzzz_xyyz_1[i] * pc_y[i];

        ta_yzzzzz_xyzz_0[i] = ta_zzzzz_xzz_0[i] * fe_0 - ta_zzzzz_xzz_1[i] * fe_0 + ta_zzzzz_xyzz_0[i] * pa_y[i] - ta_zzzzz_xyzz_1[i] * pc_y[i];

        ta_yzzzzz_xzzz_0[i] = ta_zzzzz_xzzz_0[i] * pa_y[i] - ta_zzzzz_xzzz_1[i] * pc_y[i];

        ta_yzzzzz_yyyy_0[i] =
            4.0 * ta_zzzzz_yyy_0[i] * fe_0 - 4.0 * ta_zzzzz_yyy_1[i] * fe_0 + ta_zzzzz_yyyy_0[i] * pa_y[i] - ta_zzzzz_yyyy_1[i] * pc_y[i];

        ta_yzzzzz_yyyz_0[i] =
            3.0 * ta_zzzzz_yyz_0[i] * fe_0 - 3.0 * ta_zzzzz_yyz_1[i] * fe_0 + ta_zzzzz_yyyz_0[i] * pa_y[i] - ta_zzzzz_yyyz_1[i] * pc_y[i];

        ta_yzzzzz_yyzz_0[i] =
            2.0 * ta_zzzzz_yzz_0[i] * fe_0 - 2.0 * ta_zzzzz_yzz_1[i] * fe_0 + ta_zzzzz_yyzz_0[i] * pa_y[i] - ta_zzzzz_yyzz_1[i] * pc_y[i];

        ta_yzzzzz_yzzz_0[i] = ta_zzzzz_zzz_0[i] * fe_0 - ta_zzzzz_zzz_1[i] * fe_0 + ta_zzzzz_yzzz_0[i] * pa_y[i] - ta_zzzzz_yzzz_1[i] * pc_y[i];

        ta_yzzzzz_zzzz_0[i] = ta_zzzzz_zzzz_0[i] * pa_y[i] - ta_zzzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 405-420 components of targeted buffer : IG

    auto ta_zzzzzz_xxxx_0 = pbuffer.data(idx_npot_0_ig + 405);

    auto ta_zzzzzz_xxxy_0 = pbuffer.data(idx_npot_0_ig + 406);

    auto ta_zzzzzz_xxxz_0 = pbuffer.data(idx_npot_0_ig + 407);

    auto ta_zzzzzz_xxyy_0 = pbuffer.data(idx_npot_0_ig + 408);

    auto ta_zzzzzz_xxyz_0 = pbuffer.data(idx_npot_0_ig + 409);

    auto ta_zzzzzz_xxzz_0 = pbuffer.data(idx_npot_0_ig + 410);

    auto ta_zzzzzz_xyyy_0 = pbuffer.data(idx_npot_0_ig + 411);

    auto ta_zzzzzz_xyyz_0 = pbuffer.data(idx_npot_0_ig + 412);

    auto ta_zzzzzz_xyzz_0 = pbuffer.data(idx_npot_0_ig + 413);

    auto ta_zzzzzz_xzzz_0 = pbuffer.data(idx_npot_0_ig + 414);

    auto ta_zzzzzz_yyyy_0 = pbuffer.data(idx_npot_0_ig + 415);

    auto ta_zzzzzz_yyyz_0 = pbuffer.data(idx_npot_0_ig + 416);

    auto ta_zzzzzz_yyzz_0 = pbuffer.data(idx_npot_0_ig + 417);

    auto ta_zzzzzz_yzzz_0 = pbuffer.data(idx_npot_0_ig + 418);

    auto ta_zzzzzz_zzzz_0 = pbuffer.data(idx_npot_0_ig + 419);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta_zzzz_xxxx_0,   \
                             ta_zzzz_xxxx_1,   \
                             ta_zzzz_xxxy_0,   \
                             ta_zzzz_xxxy_1,   \
                             ta_zzzz_xxxz_0,   \
                             ta_zzzz_xxxz_1,   \
                             ta_zzzz_xxyy_0,   \
                             ta_zzzz_xxyy_1,   \
                             ta_zzzz_xxyz_0,   \
                             ta_zzzz_xxyz_1,   \
                             ta_zzzz_xxzz_0,   \
                             ta_zzzz_xxzz_1,   \
                             ta_zzzz_xyyy_0,   \
                             ta_zzzz_xyyy_1,   \
                             ta_zzzz_xyyz_0,   \
                             ta_zzzz_xyyz_1,   \
                             ta_zzzz_xyzz_0,   \
                             ta_zzzz_xyzz_1,   \
                             ta_zzzz_xzzz_0,   \
                             ta_zzzz_xzzz_1,   \
                             ta_zzzz_yyyy_0,   \
                             ta_zzzz_yyyy_1,   \
                             ta_zzzz_yyyz_0,   \
                             ta_zzzz_yyyz_1,   \
                             ta_zzzz_yyzz_0,   \
                             ta_zzzz_yyzz_1,   \
                             ta_zzzz_yzzz_0,   \
                             ta_zzzz_yzzz_1,   \
                             ta_zzzz_zzzz_0,   \
                             ta_zzzz_zzzz_1,   \
                             ta_zzzzz_xxx_0,   \
                             ta_zzzzz_xxx_1,   \
                             ta_zzzzz_xxxx_0,  \
                             ta_zzzzz_xxxx_1,  \
                             ta_zzzzz_xxxy_0,  \
                             ta_zzzzz_xxxy_1,  \
                             ta_zzzzz_xxxz_0,  \
                             ta_zzzzz_xxxz_1,  \
                             ta_zzzzz_xxy_0,   \
                             ta_zzzzz_xxy_1,   \
                             ta_zzzzz_xxyy_0,  \
                             ta_zzzzz_xxyy_1,  \
                             ta_zzzzz_xxyz_0,  \
                             ta_zzzzz_xxyz_1,  \
                             ta_zzzzz_xxz_0,   \
                             ta_zzzzz_xxz_1,   \
                             ta_zzzzz_xxzz_0,  \
                             ta_zzzzz_xxzz_1,  \
                             ta_zzzzz_xyy_0,   \
                             ta_zzzzz_xyy_1,   \
                             ta_zzzzz_xyyy_0,  \
                             ta_zzzzz_xyyy_1,  \
                             ta_zzzzz_xyyz_0,  \
                             ta_zzzzz_xyyz_1,  \
                             ta_zzzzz_xyz_0,   \
                             ta_zzzzz_xyz_1,   \
                             ta_zzzzz_xyzz_0,  \
                             ta_zzzzz_xyzz_1,  \
                             ta_zzzzz_xzz_0,   \
                             ta_zzzzz_xzz_1,   \
                             ta_zzzzz_xzzz_0,  \
                             ta_zzzzz_xzzz_1,  \
                             ta_zzzzz_yyy_0,   \
                             ta_zzzzz_yyy_1,   \
                             ta_zzzzz_yyyy_0,  \
                             ta_zzzzz_yyyy_1,  \
                             ta_zzzzz_yyyz_0,  \
                             ta_zzzzz_yyyz_1,  \
                             ta_zzzzz_yyz_0,   \
                             ta_zzzzz_yyz_1,   \
                             ta_zzzzz_yyzz_0,  \
                             ta_zzzzz_yyzz_1,  \
                             ta_zzzzz_yzz_0,   \
                             ta_zzzzz_yzz_1,   \
                             ta_zzzzz_yzzz_0,  \
                             ta_zzzzz_yzzz_1,  \
                             ta_zzzzz_zzz_0,   \
                             ta_zzzzz_zzz_1,   \
                             ta_zzzzz_zzzz_0,  \
                             ta_zzzzz_zzzz_1,  \
                             ta_zzzzzz_xxxx_0, \
                             ta_zzzzzz_xxxy_0, \
                             ta_zzzzzz_xxxz_0, \
                             ta_zzzzzz_xxyy_0, \
                             ta_zzzzzz_xxyz_0, \
                             ta_zzzzzz_xxzz_0, \
                             ta_zzzzzz_xyyy_0, \
                             ta_zzzzzz_xyyz_0, \
                             ta_zzzzzz_xyzz_0, \
                             ta_zzzzzz_xzzz_0, \
                             ta_zzzzzz_yyyy_0, \
                             ta_zzzzzz_yyyz_0, \
                             ta_zzzzzz_yyzz_0, \
                             ta_zzzzzz_yzzz_0, \
                             ta_zzzzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzzz_xxxx_0[i] =
            5.0 * ta_zzzz_xxxx_0[i] * fe_0 - 5.0 * ta_zzzz_xxxx_1[i] * fe_0 + ta_zzzzz_xxxx_0[i] * pa_z[i] - ta_zzzzz_xxxx_1[i] * pc_z[i];

        ta_zzzzzz_xxxy_0[i] =
            5.0 * ta_zzzz_xxxy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxy_1[i] * fe_0 + ta_zzzzz_xxxy_0[i] * pa_z[i] - ta_zzzzz_xxxy_1[i] * pc_z[i];

        ta_zzzzzz_xxxz_0[i] = 5.0 * ta_zzzz_xxxz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxz_1[i] * fe_0 + ta_zzzzz_xxx_0[i] * fe_0 - ta_zzzzz_xxx_1[i] * fe_0 +
                              ta_zzzzz_xxxz_0[i] * pa_z[i] - ta_zzzzz_xxxz_1[i] * pc_z[i];

        ta_zzzzzz_xxyy_0[i] =
            5.0 * ta_zzzz_xxyy_0[i] * fe_0 - 5.0 * ta_zzzz_xxyy_1[i] * fe_0 + ta_zzzzz_xxyy_0[i] * pa_z[i] - ta_zzzzz_xxyy_1[i] * pc_z[i];

        ta_zzzzzz_xxyz_0[i] = 5.0 * ta_zzzz_xxyz_0[i] * fe_0 - 5.0 * ta_zzzz_xxyz_1[i] * fe_0 + ta_zzzzz_xxy_0[i] * fe_0 - ta_zzzzz_xxy_1[i] * fe_0 +
                              ta_zzzzz_xxyz_0[i] * pa_z[i] - ta_zzzzz_xxyz_1[i] * pc_z[i];

        ta_zzzzzz_xxzz_0[i] = 5.0 * ta_zzzz_xxzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xxz_0[i] * fe_0 -
                              2.0 * ta_zzzzz_xxz_1[i] * fe_0 + ta_zzzzz_xxzz_0[i] * pa_z[i] - ta_zzzzz_xxzz_1[i] * pc_z[i];

        ta_zzzzzz_xyyy_0[i] =
            5.0 * ta_zzzz_xyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xyyy_1[i] * fe_0 + ta_zzzzz_xyyy_0[i] * pa_z[i] - ta_zzzzz_xyyy_1[i] * pc_z[i];

        ta_zzzzzz_xyyz_0[i] = 5.0 * ta_zzzz_xyyz_0[i] * fe_0 - 5.0 * ta_zzzz_xyyz_1[i] * fe_0 + ta_zzzzz_xyy_0[i] * fe_0 - ta_zzzzz_xyy_1[i] * fe_0 +
                              ta_zzzzz_xyyz_0[i] * pa_z[i] - ta_zzzzz_xyyz_1[i] * pc_z[i];

        ta_zzzzzz_xyzz_0[i] = 5.0 * ta_zzzz_xyzz_0[i] * fe_0 - 5.0 * ta_zzzz_xyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xyz_0[i] * fe_0 -
                              2.0 * ta_zzzzz_xyz_1[i] * fe_0 + ta_zzzzz_xyzz_0[i] * pa_z[i] - ta_zzzzz_xyzz_1[i] * pc_z[i];

        ta_zzzzzz_xzzz_0[i] = 5.0 * ta_zzzz_xzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_xzz_0[i] * fe_0 -
                              3.0 * ta_zzzzz_xzz_1[i] * fe_0 + ta_zzzzz_xzzz_0[i] * pa_z[i] - ta_zzzzz_xzzz_1[i] * pc_z[i];

        ta_zzzzzz_yyyy_0[i] =
            5.0 * ta_zzzz_yyyy_0[i] * fe_0 - 5.0 * ta_zzzz_yyyy_1[i] * fe_0 + ta_zzzzz_yyyy_0[i] * pa_z[i] - ta_zzzzz_yyyy_1[i] * pc_z[i];

        ta_zzzzzz_yyyz_0[i] = 5.0 * ta_zzzz_yyyz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyz_1[i] * fe_0 + ta_zzzzz_yyy_0[i] * fe_0 - ta_zzzzz_yyy_1[i] * fe_0 +
                              ta_zzzzz_yyyz_0[i] * pa_z[i] - ta_zzzzz_yyyz_1[i] * pc_z[i];

        ta_zzzzzz_yyzz_0[i] = 5.0 * ta_zzzz_yyzz_0[i] * fe_0 - 5.0 * ta_zzzz_yyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_yyz_0[i] * fe_0 -
                              2.0 * ta_zzzzz_yyz_1[i] * fe_0 + ta_zzzzz_yyzz_0[i] * pa_z[i] - ta_zzzzz_yyzz_1[i] * pc_z[i];

        ta_zzzzzz_yzzz_0[i] = 5.0 * ta_zzzz_yzzz_0[i] * fe_0 - 5.0 * ta_zzzz_yzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_yzz_0[i] * fe_0 -
                              3.0 * ta_zzzzz_yzz_1[i] * fe_0 + ta_zzzzz_yzzz_0[i] * pa_z[i] - ta_zzzzz_yzzz_1[i] * pc_z[i];

        ta_zzzzzz_zzzz_0[i] = 5.0 * ta_zzzz_zzzz_0[i] * fe_0 - 5.0 * ta_zzzz_zzzz_1[i] * fe_0 + 4.0 * ta_zzzzz_zzz_0[i] * fe_0 -
                              4.0 * ta_zzzzz_zzz_1[i] * fe_0 + ta_zzzzz_zzzz_0[i] * pa_z[i] - ta_zzzzz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
