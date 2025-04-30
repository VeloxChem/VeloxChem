#include "ThreeCenterOverlapGradientPrimRecHG.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_hg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hg,
                              const size_t idx_gg,
                              const size_t idx_hf,
                              const size_t idx_hg,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : GG

    auto ts_xxxx_xxxx = pbuffer.data(idx_gg);

    auto ts_xxxx_xxxy = pbuffer.data(idx_gg + 1);

    auto ts_xxxx_xxxz = pbuffer.data(idx_gg + 2);

    auto ts_xxxx_xxyy = pbuffer.data(idx_gg + 3);

    auto ts_xxxx_xxyz = pbuffer.data(idx_gg + 4);

    auto ts_xxxx_xxzz = pbuffer.data(idx_gg + 5);

    auto ts_xxxx_xyyy = pbuffer.data(idx_gg + 6);

    auto ts_xxxx_xyyz = pbuffer.data(idx_gg + 7);

    auto ts_xxxx_xyzz = pbuffer.data(idx_gg + 8);

    auto ts_xxxx_xzzz = pbuffer.data(idx_gg + 9);

    auto ts_xxxx_yyyy = pbuffer.data(idx_gg + 10);

    auto ts_xxxx_yyyz = pbuffer.data(idx_gg + 11);

    auto ts_xxxx_yyzz = pbuffer.data(idx_gg + 12);

    auto ts_xxxx_yzzz = pbuffer.data(idx_gg + 13);

    auto ts_xxxx_zzzz = pbuffer.data(idx_gg + 14);

    auto ts_xxxy_xxxx = pbuffer.data(idx_gg + 15);

    auto ts_xxxy_xxxy = pbuffer.data(idx_gg + 16);

    auto ts_xxxy_xxxz = pbuffer.data(idx_gg + 17);

    auto ts_xxxy_xxyy = pbuffer.data(idx_gg + 18);

    auto ts_xxxy_xxyz = pbuffer.data(idx_gg + 19);

    auto ts_xxxy_xxzz = pbuffer.data(idx_gg + 20);

    auto ts_xxxy_xyyy = pbuffer.data(idx_gg + 21);

    auto ts_xxxy_xyyz = pbuffer.data(idx_gg + 22);

    auto ts_xxxy_xyzz = pbuffer.data(idx_gg + 23);

    auto ts_xxxy_xzzz = pbuffer.data(idx_gg + 24);

    auto ts_xxxy_yyyy = pbuffer.data(idx_gg + 25);

    auto ts_xxxy_yyyz = pbuffer.data(idx_gg + 26);

    auto ts_xxxy_yyzz = pbuffer.data(idx_gg + 27);

    auto ts_xxxy_yzzz = pbuffer.data(idx_gg + 28);

    auto ts_xxxy_zzzz = pbuffer.data(idx_gg + 29);

    auto ts_xxxz_xxxx = pbuffer.data(idx_gg + 30);

    auto ts_xxxz_xxxy = pbuffer.data(idx_gg + 31);

    auto ts_xxxz_xxxz = pbuffer.data(idx_gg + 32);

    auto ts_xxxz_xxyy = pbuffer.data(idx_gg + 33);

    auto ts_xxxz_xxyz = pbuffer.data(idx_gg + 34);

    auto ts_xxxz_xxzz = pbuffer.data(idx_gg + 35);

    auto ts_xxxz_xyyy = pbuffer.data(idx_gg + 36);

    auto ts_xxxz_xyyz = pbuffer.data(idx_gg + 37);

    auto ts_xxxz_xyzz = pbuffer.data(idx_gg + 38);

    auto ts_xxxz_xzzz = pbuffer.data(idx_gg + 39);

    auto ts_xxxz_yyyy = pbuffer.data(idx_gg + 40);

    auto ts_xxxz_yyyz = pbuffer.data(idx_gg + 41);

    auto ts_xxxz_yyzz = pbuffer.data(idx_gg + 42);

    auto ts_xxxz_yzzz = pbuffer.data(idx_gg + 43);

    auto ts_xxxz_zzzz = pbuffer.data(idx_gg + 44);

    auto ts_xxyy_xxxx = pbuffer.data(idx_gg + 45);

    auto ts_xxyy_xxxy = pbuffer.data(idx_gg + 46);

    auto ts_xxyy_xxxz = pbuffer.data(idx_gg + 47);

    auto ts_xxyy_xxyy = pbuffer.data(idx_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_gg + 49);

    auto ts_xxyy_xxzz = pbuffer.data(idx_gg + 50);

    auto ts_xxyy_xyyy = pbuffer.data(idx_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_gg + 53);

    auto ts_xxyy_xzzz = pbuffer.data(idx_gg + 54);

    auto ts_xxyy_yyyy = pbuffer.data(idx_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_gg + 58);

    auto ts_xxyy_zzzz = pbuffer.data(idx_gg + 59);

    auto ts_xxyz_xxxx = pbuffer.data(idx_gg + 60);

    auto ts_xxyz_xxxy = pbuffer.data(idx_gg + 61);

    auto ts_xxyz_xxxz = pbuffer.data(idx_gg + 62);

    auto ts_xxyz_xxyy = pbuffer.data(idx_gg + 63);

    auto ts_xxyz_xxyz = pbuffer.data(idx_gg + 64);

    auto ts_xxyz_xxzz = pbuffer.data(idx_gg + 65);

    auto ts_xxyz_xyyy = pbuffer.data(idx_gg + 66);

    auto ts_xxyz_xyyz = pbuffer.data(idx_gg + 67);

    auto ts_xxyz_xyzz = pbuffer.data(idx_gg + 68);

    auto ts_xxyz_xzzz = pbuffer.data(idx_gg + 69);

    auto ts_xxyz_yyyy = pbuffer.data(idx_gg + 70);

    auto ts_xxyz_yyyz = pbuffer.data(idx_gg + 71);

    auto ts_xxyz_yyzz = pbuffer.data(idx_gg + 72);

    auto ts_xxyz_yzzz = pbuffer.data(idx_gg + 73);

    auto ts_xxyz_zzzz = pbuffer.data(idx_gg + 74);

    auto ts_xxzz_xxxx = pbuffer.data(idx_gg + 75);

    auto ts_xxzz_xxxy = pbuffer.data(idx_gg + 76);

    auto ts_xxzz_xxxz = pbuffer.data(idx_gg + 77);

    auto ts_xxzz_xxyy = pbuffer.data(idx_gg + 78);

    auto ts_xxzz_xxyz = pbuffer.data(idx_gg + 79);

    auto ts_xxzz_xxzz = pbuffer.data(idx_gg + 80);

    auto ts_xxzz_xyyy = pbuffer.data(idx_gg + 81);

    auto ts_xxzz_xyyz = pbuffer.data(idx_gg + 82);

    auto ts_xxzz_xyzz = pbuffer.data(idx_gg + 83);

    auto ts_xxzz_xzzz = pbuffer.data(idx_gg + 84);

    auto ts_xxzz_yyyy = pbuffer.data(idx_gg + 85);

    auto ts_xxzz_yyyz = pbuffer.data(idx_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_gg + 89);

    auto ts_xyyy_xxxx = pbuffer.data(idx_gg + 90);

    auto ts_xyyy_xxxy = pbuffer.data(idx_gg + 91);

    auto ts_xyyy_xxxz = pbuffer.data(idx_gg + 92);

    auto ts_xyyy_xxyy = pbuffer.data(idx_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_gg + 94);

    auto ts_xyyy_xxzz = pbuffer.data(idx_gg + 95);

    auto ts_xyyy_xyyy = pbuffer.data(idx_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_gg + 98);

    auto ts_xyyy_xzzz = pbuffer.data(idx_gg + 99);

    auto ts_xyyy_yyyy = pbuffer.data(idx_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_gg + 104);

    auto ts_xyyz_xxxx = pbuffer.data(idx_gg + 105);

    auto ts_xyyz_xxxy = pbuffer.data(idx_gg + 106);

    auto ts_xyyz_xxxz = pbuffer.data(idx_gg + 107);

    auto ts_xyyz_xxyy = pbuffer.data(idx_gg + 108);

    auto ts_xyyz_xxyz = pbuffer.data(idx_gg + 109);

    auto ts_xyyz_xxzz = pbuffer.data(idx_gg + 110);

    auto ts_xyyz_xyyy = pbuffer.data(idx_gg + 111);

    auto ts_xyyz_xyyz = pbuffer.data(idx_gg + 112);

    auto ts_xyyz_xyzz = pbuffer.data(idx_gg + 113);

    auto ts_xyyz_xzzz = pbuffer.data(idx_gg + 114);

    auto ts_xyyz_yyyy = pbuffer.data(idx_gg + 115);

    auto ts_xyyz_yyyz = pbuffer.data(idx_gg + 116);

    auto ts_xyyz_yyzz = pbuffer.data(idx_gg + 117);

    auto ts_xyyz_yzzz = pbuffer.data(idx_gg + 118);

    auto ts_xyyz_zzzz = pbuffer.data(idx_gg + 119);

    auto ts_xyzz_xxxx = pbuffer.data(idx_gg + 120);

    auto ts_xyzz_xxxy = pbuffer.data(idx_gg + 121);

    auto ts_xyzz_xxxz = pbuffer.data(idx_gg + 122);

    auto ts_xyzz_xxyy = pbuffer.data(idx_gg + 123);

    auto ts_xyzz_xxyz = pbuffer.data(idx_gg + 124);

    auto ts_xyzz_xxzz = pbuffer.data(idx_gg + 125);

    auto ts_xyzz_xyyy = pbuffer.data(idx_gg + 126);

    auto ts_xyzz_xyyz = pbuffer.data(idx_gg + 127);

    auto ts_xyzz_xyzz = pbuffer.data(idx_gg + 128);

    auto ts_xyzz_xzzz = pbuffer.data(idx_gg + 129);

    auto ts_xyzz_yyyy = pbuffer.data(idx_gg + 130);

    auto ts_xyzz_yyyz = pbuffer.data(idx_gg + 131);

    auto ts_xyzz_yyzz = pbuffer.data(idx_gg + 132);

    auto ts_xyzz_yzzz = pbuffer.data(idx_gg + 133);

    auto ts_xyzz_zzzz = pbuffer.data(idx_gg + 134);

    auto ts_xzzz_xxxx = pbuffer.data(idx_gg + 135);

    auto ts_xzzz_xxxy = pbuffer.data(idx_gg + 136);

    auto ts_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto ts_xzzz_xxyy = pbuffer.data(idx_gg + 138);

    auto ts_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_gg + 140);

    auto ts_xzzz_xyyy = pbuffer.data(idx_gg + 141);

    auto ts_xzzz_xyyz = pbuffer.data(idx_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_gg + 144);

    auto ts_xzzz_yyyy = pbuffer.data(idx_gg + 145);

    auto ts_xzzz_yyyz = pbuffer.data(idx_gg + 146);

    auto ts_xzzz_yyzz = pbuffer.data(idx_gg + 147);

    auto ts_xzzz_yzzz = pbuffer.data(idx_gg + 148);

    auto ts_xzzz_zzzz = pbuffer.data(idx_gg + 149);

    auto ts_yyyy_xxxx = pbuffer.data(idx_gg + 150);

    auto ts_yyyy_xxxy = pbuffer.data(idx_gg + 151);

    auto ts_yyyy_xxxz = pbuffer.data(idx_gg + 152);

    auto ts_yyyy_xxyy = pbuffer.data(idx_gg + 153);

    auto ts_yyyy_xxyz = pbuffer.data(idx_gg + 154);

    auto ts_yyyy_xxzz = pbuffer.data(idx_gg + 155);

    auto ts_yyyy_xyyy = pbuffer.data(idx_gg + 156);

    auto ts_yyyy_xyyz = pbuffer.data(idx_gg + 157);

    auto ts_yyyy_xyzz = pbuffer.data(idx_gg + 158);

    auto ts_yyyy_xzzz = pbuffer.data(idx_gg + 159);

    auto ts_yyyy_yyyy = pbuffer.data(idx_gg + 160);

    auto ts_yyyy_yyyz = pbuffer.data(idx_gg + 161);

    auto ts_yyyy_yyzz = pbuffer.data(idx_gg + 162);

    auto ts_yyyy_yzzz = pbuffer.data(idx_gg + 163);

    auto ts_yyyy_zzzz = pbuffer.data(idx_gg + 164);

    auto ts_yyyz_xxxx = pbuffer.data(idx_gg + 165);

    auto ts_yyyz_xxxy = pbuffer.data(idx_gg + 166);

    auto ts_yyyz_xxxz = pbuffer.data(idx_gg + 167);

    auto ts_yyyz_xxyy = pbuffer.data(idx_gg + 168);

    auto ts_yyyz_xxyz = pbuffer.data(idx_gg + 169);

    auto ts_yyyz_xxzz = pbuffer.data(idx_gg + 170);

    auto ts_yyyz_xyyy = pbuffer.data(idx_gg + 171);

    auto ts_yyyz_xyyz = pbuffer.data(idx_gg + 172);

    auto ts_yyyz_xyzz = pbuffer.data(idx_gg + 173);

    auto ts_yyyz_xzzz = pbuffer.data(idx_gg + 174);

    auto ts_yyyz_yyyy = pbuffer.data(idx_gg + 175);

    auto ts_yyyz_yyyz = pbuffer.data(idx_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_gg + 179);

    auto ts_yyzz_xxxx = pbuffer.data(idx_gg + 180);

    auto ts_yyzz_xxxy = pbuffer.data(idx_gg + 181);

    auto ts_yyzz_xxxz = pbuffer.data(idx_gg + 182);

    auto ts_yyzz_xxyy = pbuffer.data(idx_gg + 183);

    auto ts_yyzz_xxyz = pbuffer.data(idx_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_gg + 185);

    auto ts_yyzz_xyyy = pbuffer.data(idx_gg + 186);

    auto ts_yyzz_xyyz = pbuffer.data(idx_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_gg + 194);

    auto ts_yzzz_xxxx = pbuffer.data(idx_gg + 195);

    auto ts_yzzz_xxxy = pbuffer.data(idx_gg + 196);

    auto ts_yzzz_xxxz = pbuffer.data(idx_gg + 197);

    auto ts_yzzz_xxyy = pbuffer.data(idx_gg + 198);

    auto ts_yzzz_xxyz = pbuffer.data(idx_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_gg + 200);

    auto ts_yzzz_xyyy = pbuffer.data(idx_gg + 201);

    auto ts_yzzz_xyyz = pbuffer.data(idx_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_gg + 203);

    auto ts_yzzz_xzzz = pbuffer.data(idx_gg + 204);

    auto ts_yzzz_yyyy = pbuffer.data(idx_gg + 205);

    auto ts_yzzz_yyyz = pbuffer.data(idx_gg + 206);

    auto ts_yzzz_yyzz = pbuffer.data(idx_gg + 207);

    auto ts_yzzz_yzzz = pbuffer.data(idx_gg + 208);

    auto ts_yzzz_zzzz = pbuffer.data(idx_gg + 209);

    auto ts_zzzz_xxxx = pbuffer.data(idx_gg + 210);

    auto ts_zzzz_xxxy = pbuffer.data(idx_gg + 211);

    auto ts_zzzz_xxxz = pbuffer.data(idx_gg + 212);

    auto ts_zzzz_xxyy = pbuffer.data(idx_gg + 213);

    auto ts_zzzz_xxyz = pbuffer.data(idx_gg + 214);

    auto ts_zzzz_xxzz = pbuffer.data(idx_gg + 215);

    auto ts_zzzz_xyyy = pbuffer.data(idx_gg + 216);

    auto ts_zzzz_xyyz = pbuffer.data(idx_gg + 217);

    auto ts_zzzz_xyzz = pbuffer.data(idx_gg + 218);

    auto ts_zzzz_xzzz = pbuffer.data(idx_gg + 219);

    auto ts_zzzz_yyyy = pbuffer.data(idx_gg + 220);

    auto ts_zzzz_yyyz = pbuffer.data(idx_gg + 221);

    auto ts_zzzz_yyzz = pbuffer.data(idx_gg + 222);

    auto ts_zzzz_yzzz = pbuffer.data(idx_gg + 223);

    auto ts_zzzz_zzzz = pbuffer.data(idx_gg + 224);

    // Set up components of auxiliary buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_hf + 9);

    auto ts_xxxxy_xxx = pbuffer.data(idx_hf + 10);

    auto ts_xxxxy_xxy = pbuffer.data(idx_hf + 11);

    auto ts_xxxxy_xxz = pbuffer.data(idx_hf + 12);

    auto ts_xxxxy_xyy = pbuffer.data(idx_hf + 13);

    auto ts_xxxxy_xyz = pbuffer.data(idx_hf + 14);

    auto ts_xxxxy_xzz = pbuffer.data(idx_hf + 15);

    auto ts_xxxxy_yyy = pbuffer.data(idx_hf + 16);

    auto ts_xxxxy_yyz = pbuffer.data(idx_hf + 17);

    auto ts_xxxxy_yzz = pbuffer.data(idx_hf + 18);

    auto ts_xxxxy_zzz = pbuffer.data(idx_hf + 19);

    auto ts_xxxxz_xxx = pbuffer.data(idx_hf + 20);

    auto ts_xxxxz_xxy = pbuffer.data(idx_hf + 21);

    auto ts_xxxxz_xxz = pbuffer.data(idx_hf + 22);

    auto ts_xxxxz_xyy = pbuffer.data(idx_hf + 23);

    auto ts_xxxxz_xyz = pbuffer.data(idx_hf + 24);

    auto ts_xxxxz_xzz = pbuffer.data(idx_hf + 25);

    auto ts_xxxxz_yyy = pbuffer.data(idx_hf + 26);

    auto ts_xxxxz_yyz = pbuffer.data(idx_hf + 27);

    auto ts_xxxxz_yzz = pbuffer.data(idx_hf + 28);

    auto ts_xxxxz_zzz = pbuffer.data(idx_hf + 29);

    auto ts_xxxyy_xxx = pbuffer.data(idx_hf + 30);

    auto ts_xxxyy_xxy = pbuffer.data(idx_hf + 31);

    auto ts_xxxyy_xxz = pbuffer.data(idx_hf + 32);

    auto ts_xxxyy_xyy = pbuffer.data(idx_hf + 33);

    auto ts_xxxyy_xyz = pbuffer.data(idx_hf + 34);

    auto ts_xxxyy_xzz = pbuffer.data(idx_hf + 35);

    auto ts_xxxyy_yyy = pbuffer.data(idx_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_hf + 38);

    auto ts_xxxyy_zzz = pbuffer.data(idx_hf + 39);

    auto ts_xxxyz_xxx = pbuffer.data(idx_hf + 40);

    auto ts_xxxyz_xxy = pbuffer.data(idx_hf + 41);

    auto ts_xxxyz_xxz = pbuffer.data(idx_hf + 42);

    auto ts_xxxyz_xyy = pbuffer.data(idx_hf + 43);

    auto ts_xxxyz_xyz = pbuffer.data(idx_hf + 44);

    auto ts_xxxyz_xzz = pbuffer.data(idx_hf + 45);

    auto ts_xxxyz_yyy = pbuffer.data(idx_hf + 46);

    auto ts_xxxyz_yyz = pbuffer.data(idx_hf + 47);

    auto ts_xxxyz_yzz = pbuffer.data(idx_hf + 48);

    auto ts_xxxyz_zzz = pbuffer.data(idx_hf + 49);

    auto ts_xxxzz_xxx = pbuffer.data(idx_hf + 50);

    auto ts_xxxzz_xxy = pbuffer.data(idx_hf + 51);

    auto ts_xxxzz_xxz = pbuffer.data(idx_hf + 52);

    auto ts_xxxzz_xyy = pbuffer.data(idx_hf + 53);

    auto ts_xxxzz_xyz = pbuffer.data(idx_hf + 54);

    auto ts_xxxzz_xzz = pbuffer.data(idx_hf + 55);

    auto ts_xxxzz_yyy = pbuffer.data(idx_hf + 56);

    auto ts_xxxzz_yyz = pbuffer.data(idx_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_hf + 59);

    auto ts_xxyyy_xxx = pbuffer.data(idx_hf + 60);

    auto ts_xxyyy_xxy = pbuffer.data(idx_hf + 61);

    auto ts_xxyyy_xxz = pbuffer.data(idx_hf + 62);

    auto ts_xxyyy_xyy = pbuffer.data(idx_hf + 63);

    auto ts_xxyyy_xyz = pbuffer.data(idx_hf + 64);

    auto ts_xxyyy_xzz = pbuffer.data(idx_hf + 65);

    auto ts_xxyyy_yyy = pbuffer.data(idx_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_hf + 68);

    auto ts_xxyyy_zzz = pbuffer.data(idx_hf + 69);

    auto ts_xxyyz_xxx = pbuffer.data(idx_hf + 70);

    auto ts_xxyyz_xxy = pbuffer.data(idx_hf + 71);

    auto ts_xxyyz_xxz = pbuffer.data(idx_hf + 72);

    auto ts_xxyyz_xyy = pbuffer.data(idx_hf + 73);

    auto ts_xxyyz_xyz = pbuffer.data(idx_hf + 74);

    auto ts_xxyyz_xzz = pbuffer.data(idx_hf + 75);

    auto ts_xxyyz_yyy = pbuffer.data(idx_hf + 76);

    auto ts_xxyyz_yyz = pbuffer.data(idx_hf + 77);

    auto ts_xxyyz_yzz = pbuffer.data(idx_hf + 78);

    auto ts_xxyyz_zzz = pbuffer.data(idx_hf + 79);

    auto ts_xxyzz_xxx = pbuffer.data(idx_hf + 80);

    auto ts_xxyzz_xxy = pbuffer.data(idx_hf + 81);

    auto ts_xxyzz_xxz = pbuffer.data(idx_hf + 82);

    auto ts_xxyzz_xyy = pbuffer.data(idx_hf + 83);

    auto ts_xxyzz_xyz = pbuffer.data(idx_hf + 84);

    auto ts_xxyzz_xzz = pbuffer.data(idx_hf + 85);

    auto ts_xxyzz_yyy = pbuffer.data(idx_hf + 86);

    auto ts_xxyzz_yyz = pbuffer.data(idx_hf + 87);

    auto ts_xxyzz_yzz = pbuffer.data(idx_hf + 88);

    auto ts_xxyzz_zzz = pbuffer.data(idx_hf + 89);

    auto ts_xxzzz_xxx = pbuffer.data(idx_hf + 90);

    auto ts_xxzzz_xxy = pbuffer.data(idx_hf + 91);

    auto ts_xxzzz_xxz = pbuffer.data(idx_hf + 92);

    auto ts_xxzzz_xyy = pbuffer.data(idx_hf + 93);

    auto ts_xxzzz_xyz = pbuffer.data(idx_hf + 94);

    auto ts_xxzzz_xzz = pbuffer.data(idx_hf + 95);

    auto ts_xxzzz_yyy = pbuffer.data(idx_hf + 96);

    auto ts_xxzzz_yyz = pbuffer.data(idx_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_hf + 99);

    auto ts_xyyyy_xxx = pbuffer.data(idx_hf + 100);

    auto ts_xyyyy_xxy = pbuffer.data(idx_hf + 101);

    auto ts_xyyyy_xxz = pbuffer.data(idx_hf + 102);

    auto ts_xyyyy_xyy = pbuffer.data(idx_hf + 103);

    auto ts_xyyyy_xyz = pbuffer.data(idx_hf + 104);

    auto ts_xyyyy_xzz = pbuffer.data(idx_hf + 105);

    auto ts_xyyyy_yyy = pbuffer.data(idx_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_hf + 108);

    auto ts_xyyyy_zzz = pbuffer.data(idx_hf + 109);

    auto ts_xyyyz_xxx = pbuffer.data(idx_hf + 110);

    auto ts_xyyyz_xxy = pbuffer.data(idx_hf + 111);

    auto ts_xyyyz_xxz = pbuffer.data(idx_hf + 112);

    auto ts_xyyyz_xyy = pbuffer.data(idx_hf + 113);

    auto ts_xyyyz_xyz = pbuffer.data(idx_hf + 114);

    auto ts_xyyyz_xzz = pbuffer.data(idx_hf + 115);

    auto ts_xyyyz_yyy = pbuffer.data(idx_hf + 116);

    auto ts_xyyyz_yyz = pbuffer.data(idx_hf + 117);

    auto ts_xyyyz_yzz = pbuffer.data(idx_hf + 118);

    auto ts_xyyyz_zzz = pbuffer.data(idx_hf + 119);

    auto ts_xyyzz_xxx = pbuffer.data(idx_hf + 120);

    auto ts_xyyzz_xxy = pbuffer.data(idx_hf + 121);

    auto ts_xyyzz_xxz = pbuffer.data(idx_hf + 122);

    auto ts_xyyzz_xyy = pbuffer.data(idx_hf + 123);

    auto ts_xyyzz_xyz = pbuffer.data(idx_hf + 124);

    auto ts_xyyzz_xzz = pbuffer.data(idx_hf + 125);

    auto ts_xyyzz_yyy = pbuffer.data(idx_hf + 126);

    auto ts_xyyzz_yyz = pbuffer.data(idx_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_hf + 128);

    auto ts_xyyzz_zzz = pbuffer.data(idx_hf + 129);

    auto ts_xyzzz_xxx = pbuffer.data(idx_hf + 130);

    auto ts_xyzzz_xxy = pbuffer.data(idx_hf + 131);

    auto ts_xyzzz_xxz = pbuffer.data(idx_hf + 132);

    auto ts_xyzzz_xyy = pbuffer.data(idx_hf + 133);

    auto ts_xyzzz_xyz = pbuffer.data(idx_hf + 134);

    auto ts_xyzzz_xzz = pbuffer.data(idx_hf + 135);

    auto ts_xyzzz_yyy = pbuffer.data(idx_hf + 136);

    auto ts_xyzzz_yyz = pbuffer.data(idx_hf + 137);

    auto ts_xyzzz_yzz = pbuffer.data(idx_hf + 138);

    auto ts_xyzzz_zzz = pbuffer.data(idx_hf + 139);

    auto ts_xzzzz_xxx = pbuffer.data(idx_hf + 140);

    auto ts_xzzzz_xxy = pbuffer.data(idx_hf + 141);

    auto ts_xzzzz_xxz = pbuffer.data(idx_hf + 142);

    auto ts_xzzzz_xyy = pbuffer.data(idx_hf + 143);

    auto ts_xzzzz_xyz = pbuffer.data(idx_hf + 144);

    auto ts_xzzzz_xzz = pbuffer.data(idx_hf + 145);

    auto ts_xzzzz_yyy = pbuffer.data(idx_hf + 146);

    auto ts_xzzzz_yyz = pbuffer.data(idx_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_hf + 149);

    auto ts_yyyyy_xxx = pbuffer.data(idx_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_hf + 159);

    auto ts_yyyyz_xxx = pbuffer.data(idx_hf + 160);

    auto ts_yyyyz_xxy = pbuffer.data(idx_hf + 161);

    auto ts_yyyyz_xxz = pbuffer.data(idx_hf + 162);

    auto ts_yyyyz_xyy = pbuffer.data(idx_hf + 163);

    auto ts_yyyyz_xyz = pbuffer.data(idx_hf + 164);

    auto ts_yyyyz_xzz = pbuffer.data(idx_hf + 165);

    auto ts_yyyyz_yyy = pbuffer.data(idx_hf + 166);

    auto ts_yyyyz_yyz = pbuffer.data(idx_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_hf + 169);

    auto ts_yyyzz_xxx = pbuffer.data(idx_hf + 170);

    auto ts_yyyzz_xxy = pbuffer.data(idx_hf + 171);

    auto ts_yyyzz_xxz = pbuffer.data(idx_hf + 172);

    auto ts_yyyzz_xyy = pbuffer.data(idx_hf + 173);

    auto ts_yyyzz_xyz = pbuffer.data(idx_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_hf + 179);

    auto ts_yyzzz_xxx = pbuffer.data(idx_hf + 180);

    auto ts_yyzzz_xxy = pbuffer.data(idx_hf + 181);

    auto ts_yyzzz_xxz = pbuffer.data(idx_hf + 182);

    auto ts_yyzzz_xyy = pbuffer.data(idx_hf + 183);

    auto ts_yyzzz_xyz = pbuffer.data(idx_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_hf + 189);

    auto ts_yzzzz_xxx = pbuffer.data(idx_hf + 190);

    auto ts_yzzzz_xxy = pbuffer.data(idx_hf + 191);

    auto ts_yzzzz_xxz = pbuffer.data(idx_hf + 192);

    auto ts_yzzzz_xyy = pbuffer.data(idx_hf + 193);

    auto ts_yzzzz_xyz = pbuffer.data(idx_hf + 194);

    auto ts_yzzzz_xzz = pbuffer.data(idx_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_hf + 199);

    auto ts_zzzzz_xxx = pbuffer.data(idx_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_hf + 209);

    // Set up components of auxiliary buffer : HG

    auto ts_xxxxx_xxxx = pbuffer.data(idx_hg);

    auto ts_xxxxx_xxxy = pbuffer.data(idx_hg + 1);

    auto ts_xxxxx_xxxz = pbuffer.data(idx_hg + 2);

    auto ts_xxxxx_xxyy = pbuffer.data(idx_hg + 3);

    auto ts_xxxxx_xxyz = pbuffer.data(idx_hg + 4);

    auto ts_xxxxx_xxzz = pbuffer.data(idx_hg + 5);

    auto ts_xxxxx_xyyy = pbuffer.data(idx_hg + 6);

    auto ts_xxxxx_xyyz = pbuffer.data(idx_hg + 7);

    auto ts_xxxxx_xyzz = pbuffer.data(idx_hg + 8);

    auto ts_xxxxx_xzzz = pbuffer.data(idx_hg + 9);

    auto ts_xxxxx_yyyy = pbuffer.data(idx_hg + 10);

    auto ts_xxxxx_yyyz = pbuffer.data(idx_hg + 11);

    auto ts_xxxxx_yyzz = pbuffer.data(idx_hg + 12);

    auto ts_xxxxx_yzzz = pbuffer.data(idx_hg + 13);

    auto ts_xxxxx_zzzz = pbuffer.data(idx_hg + 14);

    auto ts_xxxxy_xxxx = pbuffer.data(idx_hg + 15);

    auto ts_xxxxy_xxxy = pbuffer.data(idx_hg + 16);

    auto ts_xxxxy_xxxz = pbuffer.data(idx_hg + 17);

    auto ts_xxxxy_xxyy = pbuffer.data(idx_hg + 18);

    auto ts_xxxxy_xxyz = pbuffer.data(idx_hg + 19);

    auto ts_xxxxy_xxzz = pbuffer.data(idx_hg + 20);

    auto ts_xxxxy_xyyy = pbuffer.data(idx_hg + 21);

    auto ts_xxxxy_xyyz = pbuffer.data(idx_hg + 22);

    auto ts_xxxxy_xyzz = pbuffer.data(idx_hg + 23);

    auto ts_xxxxy_xzzz = pbuffer.data(idx_hg + 24);

    auto ts_xxxxy_yyyy = pbuffer.data(idx_hg + 25);

    auto ts_xxxxy_yyyz = pbuffer.data(idx_hg + 26);

    auto ts_xxxxy_yyzz = pbuffer.data(idx_hg + 27);

    auto ts_xxxxy_yzzz = pbuffer.data(idx_hg + 28);

    auto ts_xxxxy_zzzz = pbuffer.data(idx_hg + 29);

    auto ts_xxxxz_xxxx = pbuffer.data(idx_hg + 30);

    auto ts_xxxxz_xxxy = pbuffer.data(idx_hg + 31);

    auto ts_xxxxz_xxxz = pbuffer.data(idx_hg + 32);

    auto ts_xxxxz_xxyy = pbuffer.data(idx_hg + 33);

    auto ts_xxxxz_xxyz = pbuffer.data(idx_hg + 34);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_hg + 35);

    auto ts_xxxxz_xyyy = pbuffer.data(idx_hg + 36);

    auto ts_xxxxz_xyyz = pbuffer.data(idx_hg + 37);

    auto ts_xxxxz_xyzz = pbuffer.data(idx_hg + 38);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_hg + 39);

    auto ts_xxxxz_yyyy = pbuffer.data(idx_hg + 40);

    auto ts_xxxxz_yyyz = pbuffer.data(idx_hg + 41);

    auto ts_xxxxz_yyzz = pbuffer.data(idx_hg + 42);

    auto ts_xxxxz_yzzz = pbuffer.data(idx_hg + 43);

    auto ts_xxxxz_zzzz = pbuffer.data(idx_hg + 44);

    auto ts_xxxyy_xxxx = pbuffer.data(idx_hg + 45);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_hg + 46);

    auto ts_xxxyy_xxxz = pbuffer.data(idx_hg + 47);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_hg + 48);

    auto ts_xxxyy_xxyz = pbuffer.data(idx_hg + 49);

    auto ts_xxxyy_xxzz = pbuffer.data(idx_hg + 50);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_hg + 51);

    auto ts_xxxyy_xyyz = pbuffer.data(idx_hg + 52);

    auto ts_xxxyy_xyzz = pbuffer.data(idx_hg + 53);

    auto ts_xxxyy_xzzz = pbuffer.data(idx_hg + 54);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_hg + 58);

    auto ts_xxxyy_zzzz = pbuffer.data(idx_hg + 59);

    auto ts_xxxyz_xxxx = pbuffer.data(idx_hg + 60);

    auto ts_xxxyz_xxxy = pbuffer.data(idx_hg + 61);

    auto ts_xxxyz_xxxz = pbuffer.data(idx_hg + 62);

    auto ts_xxxyz_xxyy = pbuffer.data(idx_hg + 63);

    auto ts_xxxyz_xxyz = pbuffer.data(idx_hg + 64);

    auto ts_xxxyz_xxzz = pbuffer.data(idx_hg + 65);

    auto ts_xxxyz_xyyy = pbuffer.data(idx_hg + 66);

    auto ts_xxxyz_xyyz = pbuffer.data(idx_hg + 67);

    auto ts_xxxyz_xyzz = pbuffer.data(idx_hg + 68);

    auto ts_xxxyz_xzzz = pbuffer.data(idx_hg + 69);

    auto ts_xxxyz_yyyy = pbuffer.data(idx_hg + 70);

    auto ts_xxxyz_yyyz = pbuffer.data(idx_hg + 71);

    auto ts_xxxyz_yyzz = pbuffer.data(idx_hg + 72);

    auto ts_xxxyz_yzzz = pbuffer.data(idx_hg + 73);

    auto ts_xxxyz_zzzz = pbuffer.data(idx_hg + 74);

    auto ts_xxxzz_xxxx = pbuffer.data(idx_hg + 75);

    auto ts_xxxzz_xxxy = pbuffer.data(idx_hg + 76);

    auto ts_xxxzz_xxxz = pbuffer.data(idx_hg + 77);

    auto ts_xxxzz_xxyy = pbuffer.data(idx_hg + 78);

    auto ts_xxxzz_xxyz = pbuffer.data(idx_hg + 79);

    auto ts_xxxzz_xxzz = pbuffer.data(idx_hg + 80);

    auto ts_xxxzz_xyyy = pbuffer.data(idx_hg + 81);

    auto ts_xxxzz_xyyz = pbuffer.data(idx_hg + 82);

    auto ts_xxxzz_xyzz = pbuffer.data(idx_hg + 83);

    auto ts_xxxzz_xzzz = pbuffer.data(idx_hg + 84);

    auto ts_xxxzz_yyyy = pbuffer.data(idx_hg + 85);

    auto ts_xxxzz_yyyz = pbuffer.data(idx_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_hg + 89);

    auto ts_xxyyy_xxxx = pbuffer.data(idx_hg + 90);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_hg + 91);

    auto ts_xxyyy_xxxz = pbuffer.data(idx_hg + 92);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_hg + 93);

    auto ts_xxyyy_xxyz = pbuffer.data(idx_hg + 94);

    auto ts_xxyyy_xxzz = pbuffer.data(idx_hg + 95);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_hg + 96);

    auto ts_xxyyy_xyyz = pbuffer.data(idx_hg + 97);

    auto ts_xxyyy_xyzz = pbuffer.data(idx_hg + 98);

    auto ts_xxyyy_xzzz = pbuffer.data(idx_hg + 99);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_hg + 103);

    auto ts_xxyyy_zzzz = pbuffer.data(idx_hg + 104);

    auto ts_xxyyz_xxxx = pbuffer.data(idx_hg + 105);

    auto ts_xxyyz_xxxy = pbuffer.data(idx_hg + 106);

    auto ts_xxyyz_xxxz = pbuffer.data(idx_hg + 107);

    auto ts_xxyyz_xxyy = pbuffer.data(idx_hg + 108);

    auto ts_xxyyz_xxyz = pbuffer.data(idx_hg + 109);

    auto ts_xxyyz_xxzz = pbuffer.data(idx_hg + 110);

    auto ts_xxyyz_xyyy = pbuffer.data(idx_hg + 111);

    auto ts_xxyyz_xyyz = pbuffer.data(idx_hg + 112);

    auto ts_xxyyz_xyzz = pbuffer.data(idx_hg + 113);

    auto ts_xxyyz_xzzz = pbuffer.data(idx_hg + 114);

    auto ts_xxyyz_yyyy = pbuffer.data(idx_hg + 115);

    auto ts_xxyyz_yyyz = pbuffer.data(idx_hg + 116);

    auto ts_xxyyz_yyzz = pbuffer.data(idx_hg + 117);

    auto ts_xxyyz_yzzz = pbuffer.data(idx_hg + 118);

    auto ts_xxyyz_zzzz = pbuffer.data(idx_hg + 119);

    auto ts_xxyzz_xxxx = pbuffer.data(idx_hg + 120);

    auto ts_xxyzz_xxxy = pbuffer.data(idx_hg + 121);

    auto ts_xxyzz_xxxz = pbuffer.data(idx_hg + 122);

    auto ts_xxyzz_xxyy = pbuffer.data(idx_hg + 123);

    auto ts_xxyzz_xxyz = pbuffer.data(idx_hg + 124);

    auto ts_xxyzz_xxzz = pbuffer.data(idx_hg + 125);

    auto ts_xxyzz_xyyy = pbuffer.data(idx_hg + 126);

    auto ts_xxyzz_xyyz = pbuffer.data(idx_hg + 127);

    auto ts_xxyzz_xyzz = pbuffer.data(idx_hg + 128);

    auto ts_xxyzz_xzzz = pbuffer.data(idx_hg + 129);

    auto ts_xxyzz_yyyy = pbuffer.data(idx_hg + 130);

    auto ts_xxyzz_yyyz = pbuffer.data(idx_hg + 131);

    auto ts_xxyzz_yyzz = pbuffer.data(idx_hg + 132);

    auto ts_xxyzz_yzzz = pbuffer.data(idx_hg + 133);

    auto ts_xxyzz_zzzz = pbuffer.data(idx_hg + 134);

    auto ts_xxzzz_xxxx = pbuffer.data(idx_hg + 135);

    auto ts_xxzzz_xxxy = pbuffer.data(idx_hg + 136);

    auto ts_xxzzz_xxxz = pbuffer.data(idx_hg + 137);

    auto ts_xxzzz_xxyy = pbuffer.data(idx_hg + 138);

    auto ts_xxzzz_xxyz = pbuffer.data(idx_hg + 139);

    auto ts_xxzzz_xxzz = pbuffer.data(idx_hg + 140);

    auto ts_xxzzz_xyyy = pbuffer.data(idx_hg + 141);

    auto ts_xxzzz_xyyz = pbuffer.data(idx_hg + 142);

    auto ts_xxzzz_xyzz = pbuffer.data(idx_hg + 143);

    auto ts_xxzzz_xzzz = pbuffer.data(idx_hg + 144);

    auto ts_xxzzz_yyyy = pbuffer.data(idx_hg + 145);

    auto ts_xxzzz_yyyz = pbuffer.data(idx_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_hg + 149);

    auto ts_xyyyy_xxxx = pbuffer.data(idx_hg + 150);

    auto ts_xyyyy_xxxy = pbuffer.data(idx_hg + 151);

    auto ts_xyyyy_xxxz = pbuffer.data(idx_hg + 152);

    auto ts_xyyyy_xxyy = pbuffer.data(idx_hg + 153);

    auto ts_xyyyy_xxyz = pbuffer.data(idx_hg + 154);

    auto ts_xyyyy_xxzz = pbuffer.data(idx_hg + 155);

    auto ts_xyyyy_xyyy = pbuffer.data(idx_hg + 156);

    auto ts_xyyyy_xyyz = pbuffer.data(idx_hg + 157);

    auto ts_xyyyy_xyzz = pbuffer.data(idx_hg + 158);

    auto ts_xyyyy_xzzz = pbuffer.data(idx_hg + 159);

    auto ts_xyyyy_yyyy = pbuffer.data(idx_hg + 160);

    auto ts_xyyyy_yyyz = pbuffer.data(idx_hg + 161);

    auto ts_xyyyy_yyzz = pbuffer.data(idx_hg + 162);

    auto ts_xyyyy_yzzz = pbuffer.data(idx_hg + 163);

    auto ts_xyyyy_zzzz = pbuffer.data(idx_hg + 164);

    auto ts_xyyyz_xxxx = pbuffer.data(idx_hg + 165);

    auto ts_xyyyz_xxxy = pbuffer.data(idx_hg + 166);

    auto ts_xyyyz_xxxz = pbuffer.data(idx_hg + 167);

    auto ts_xyyyz_xxyy = pbuffer.data(idx_hg + 168);

    auto ts_xyyyz_xxyz = pbuffer.data(idx_hg + 169);

    auto ts_xyyyz_xxzz = pbuffer.data(idx_hg + 170);

    auto ts_xyyyz_xyyy = pbuffer.data(idx_hg + 171);

    auto ts_xyyyz_xyyz = pbuffer.data(idx_hg + 172);

    auto ts_xyyyz_xyzz = pbuffer.data(idx_hg + 173);

    auto ts_xyyyz_xzzz = pbuffer.data(idx_hg + 174);

    auto ts_xyyyz_yyyy = pbuffer.data(idx_hg + 175);

    auto ts_xyyyz_yyyz = pbuffer.data(idx_hg + 176);

    auto ts_xyyyz_yyzz = pbuffer.data(idx_hg + 177);

    auto ts_xyyyz_yzzz = pbuffer.data(idx_hg + 178);

    auto ts_xyyyz_zzzz = pbuffer.data(idx_hg + 179);

    auto ts_xyyzz_xxxx = pbuffer.data(idx_hg + 180);

    auto ts_xyyzz_xxxy = pbuffer.data(idx_hg + 181);

    auto ts_xyyzz_xxxz = pbuffer.data(idx_hg + 182);

    auto ts_xyyzz_xxyy = pbuffer.data(idx_hg + 183);

    auto ts_xyyzz_xxyz = pbuffer.data(idx_hg + 184);

    auto ts_xyyzz_xxzz = pbuffer.data(idx_hg + 185);

    auto ts_xyyzz_xyyy = pbuffer.data(idx_hg + 186);

    auto ts_xyyzz_xyyz = pbuffer.data(idx_hg + 187);

    auto ts_xyyzz_xyzz = pbuffer.data(idx_hg + 188);

    auto ts_xyyzz_xzzz = pbuffer.data(idx_hg + 189);

    auto ts_xyyzz_yyyy = pbuffer.data(idx_hg + 190);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_hg + 193);

    auto ts_xyyzz_zzzz = pbuffer.data(idx_hg + 194);

    auto ts_xyzzz_xxxx = pbuffer.data(idx_hg + 195);

    auto ts_xyzzz_xxxy = pbuffer.data(idx_hg + 196);

    auto ts_xyzzz_xxxz = pbuffer.data(idx_hg + 197);

    auto ts_xyzzz_xxyy = pbuffer.data(idx_hg + 198);

    auto ts_xyzzz_xxyz = pbuffer.data(idx_hg + 199);

    auto ts_xyzzz_xxzz = pbuffer.data(idx_hg + 200);

    auto ts_xyzzz_xyyy = pbuffer.data(idx_hg + 201);

    auto ts_xyzzz_xyyz = pbuffer.data(idx_hg + 202);

    auto ts_xyzzz_xyzz = pbuffer.data(idx_hg + 203);

    auto ts_xyzzz_xzzz = pbuffer.data(idx_hg + 204);

    auto ts_xyzzz_yyyy = pbuffer.data(idx_hg + 205);

    auto ts_xyzzz_yyyz = pbuffer.data(idx_hg + 206);

    auto ts_xyzzz_yyzz = pbuffer.data(idx_hg + 207);

    auto ts_xyzzz_yzzz = pbuffer.data(idx_hg + 208);

    auto ts_xyzzz_zzzz = pbuffer.data(idx_hg + 209);

    auto ts_xzzzz_xxxx = pbuffer.data(idx_hg + 210);

    auto ts_xzzzz_xxxy = pbuffer.data(idx_hg + 211);

    auto ts_xzzzz_xxxz = pbuffer.data(idx_hg + 212);

    auto ts_xzzzz_xxyy = pbuffer.data(idx_hg + 213);

    auto ts_xzzzz_xxyz = pbuffer.data(idx_hg + 214);

    auto ts_xzzzz_xxzz = pbuffer.data(idx_hg + 215);

    auto ts_xzzzz_xyyy = pbuffer.data(idx_hg + 216);

    auto ts_xzzzz_xyyz = pbuffer.data(idx_hg + 217);

    auto ts_xzzzz_xyzz = pbuffer.data(idx_hg + 218);

    auto ts_xzzzz_xzzz = pbuffer.data(idx_hg + 219);

    auto ts_xzzzz_yyyy = pbuffer.data(idx_hg + 220);

    auto ts_xzzzz_yyyz = pbuffer.data(idx_hg + 221);

    auto ts_xzzzz_yyzz = pbuffer.data(idx_hg + 222);

    auto ts_xzzzz_yzzz = pbuffer.data(idx_hg + 223);

    auto ts_xzzzz_zzzz = pbuffer.data(idx_hg + 224);

    auto ts_yyyyy_xxxx = pbuffer.data(idx_hg + 225);

    auto ts_yyyyy_xxxy = pbuffer.data(idx_hg + 226);

    auto ts_yyyyy_xxxz = pbuffer.data(idx_hg + 227);

    auto ts_yyyyy_xxyy = pbuffer.data(idx_hg + 228);

    auto ts_yyyyy_xxyz = pbuffer.data(idx_hg + 229);

    auto ts_yyyyy_xxzz = pbuffer.data(idx_hg + 230);

    auto ts_yyyyy_xyyy = pbuffer.data(idx_hg + 231);

    auto ts_yyyyy_xyyz = pbuffer.data(idx_hg + 232);

    auto ts_yyyyy_xyzz = pbuffer.data(idx_hg + 233);

    auto ts_yyyyy_xzzz = pbuffer.data(idx_hg + 234);

    auto ts_yyyyy_yyyy = pbuffer.data(idx_hg + 235);

    auto ts_yyyyy_yyyz = pbuffer.data(idx_hg + 236);

    auto ts_yyyyy_yyzz = pbuffer.data(idx_hg + 237);

    auto ts_yyyyy_yzzz = pbuffer.data(idx_hg + 238);

    auto ts_yyyyy_zzzz = pbuffer.data(idx_hg + 239);

    auto ts_yyyyz_xxxx = pbuffer.data(idx_hg + 240);

    auto ts_yyyyz_xxxy = pbuffer.data(idx_hg + 241);

    auto ts_yyyyz_xxxz = pbuffer.data(idx_hg + 242);

    auto ts_yyyyz_xxyy = pbuffer.data(idx_hg + 243);

    auto ts_yyyyz_xxyz = pbuffer.data(idx_hg + 244);

    auto ts_yyyyz_xxzz = pbuffer.data(idx_hg + 245);

    auto ts_yyyyz_xyyy = pbuffer.data(idx_hg + 246);

    auto ts_yyyyz_xyyz = pbuffer.data(idx_hg + 247);

    auto ts_yyyyz_xyzz = pbuffer.data(idx_hg + 248);

    auto ts_yyyyz_xzzz = pbuffer.data(idx_hg + 249);

    auto ts_yyyyz_yyyy = pbuffer.data(idx_hg + 250);

    auto ts_yyyyz_yyyz = pbuffer.data(idx_hg + 251);

    auto ts_yyyyz_yyzz = pbuffer.data(idx_hg + 252);

    auto ts_yyyyz_yzzz = pbuffer.data(idx_hg + 253);

    auto ts_yyyyz_zzzz = pbuffer.data(idx_hg + 254);

    auto ts_yyyzz_xxxx = pbuffer.data(idx_hg + 255);

    auto ts_yyyzz_xxxy = pbuffer.data(idx_hg + 256);

    auto ts_yyyzz_xxxz = pbuffer.data(idx_hg + 257);

    auto ts_yyyzz_xxyy = pbuffer.data(idx_hg + 258);

    auto ts_yyyzz_xxyz = pbuffer.data(idx_hg + 259);

    auto ts_yyyzz_xxzz = pbuffer.data(idx_hg + 260);

    auto ts_yyyzz_xyyy = pbuffer.data(idx_hg + 261);

    auto ts_yyyzz_xyyz = pbuffer.data(idx_hg + 262);

    auto ts_yyyzz_xyzz = pbuffer.data(idx_hg + 263);

    auto ts_yyyzz_xzzz = pbuffer.data(idx_hg + 264);

    auto ts_yyyzz_yyyy = pbuffer.data(idx_hg + 265);

    auto ts_yyyzz_yyyz = pbuffer.data(idx_hg + 266);

    auto ts_yyyzz_yyzz = pbuffer.data(idx_hg + 267);

    auto ts_yyyzz_yzzz = pbuffer.data(idx_hg + 268);

    auto ts_yyyzz_zzzz = pbuffer.data(idx_hg + 269);

    auto ts_yyzzz_xxxx = pbuffer.data(idx_hg + 270);

    auto ts_yyzzz_xxxy = pbuffer.data(idx_hg + 271);

    auto ts_yyzzz_xxxz = pbuffer.data(idx_hg + 272);

    auto ts_yyzzz_xxyy = pbuffer.data(idx_hg + 273);

    auto ts_yyzzz_xxyz = pbuffer.data(idx_hg + 274);

    auto ts_yyzzz_xxzz = pbuffer.data(idx_hg + 275);

    auto ts_yyzzz_xyyy = pbuffer.data(idx_hg + 276);

    auto ts_yyzzz_xyyz = pbuffer.data(idx_hg + 277);

    auto ts_yyzzz_xyzz = pbuffer.data(idx_hg + 278);

    auto ts_yyzzz_xzzz = pbuffer.data(idx_hg + 279);

    auto ts_yyzzz_yyyy = pbuffer.data(idx_hg + 280);

    auto ts_yyzzz_yyyz = pbuffer.data(idx_hg + 281);

    auto ts_yyzzz_yyzz = pbuffer.data(idx_hg + 282);

    auto ts_yyzzz_yzzz = pbuffer.data(idx_hg + 283);

    auto ts_yyzzz_zzzz = pbuffer.data(idx_hg + 284);

    auto ts_yzzzz_xxxx = pbuffer.data(idx_hg + 285);

    auto ts_yzzzz_xxxy = pbuffer.data(idx_hg + 286);

    auto ts_yzzzz_xxxz = pbuffer.data(idx_hg + 287);

    auto ts_yzzzz_xxyy = pbuffer.data(idx_hg + 288);

    auto ts_yzzzz_xxyz = pbuffer.data(idx_hg + 289);

    auto ts_yzzzz_xxzz = pbuffer.data(idx_hg + 290);

    auto ts_yzzzz_xyyy = pbuffer.data(idx_hg + 291);

    auto ts_yzzzz_xyyz = pbuffer.data(idx_hg + 292);

    auto ts_yzzzz_xyzz = pbuffer.data(idx_hg + 293);

    auto ts_yzzzz_xzzz = pbuffer.data(idx_hg + 294);

    auto ts_yzzzz_yyyy = pbuffer.data(idx_hg + 295);

    auto ts_yzzzz_yyyz = pbuffer.data(idx_hg + 296);

    auto ts_yzzzz_yyzz = pbuffer.data(idx_hg + 297);

    auto ts_yzzzz_yzzz = pbuffer.data(idx_hg + 298);

    auto ts_yzzzz_zzzz = pbuffer.data(idx_hg + 299);

    auto ts_zzzzz_xxxx = pbuffer.data(idx_hg + 300);

    auto ts_zzzzz_xxxy = pbuffer.data(idx_hg + 301);

    auto ts_zzzzz_xxxz = pbuffer.data(idx_hg + 302);

    auto ts_zzzzz_xxyy = pbuffer.data(idx_hg + 303);

    auto ts_zzzzz_xxyz = pbuffer.data(idx_hg + 304);

    auto ts_zzzzz_xxzz = pbuffer.data(idx_hg + 305);

    auto ts_zzzzz_xyyy = pbuffer.data(idx_hg + 306);

    auto ts_zzzzz_xyyz = pbuffer.data(idx_hg + 307);

    auto ts_zzzzz_xyzz = pbuffer.data(idx_hg + 308);

    auto ts_zzzzz_xzzz = pbuffer.data(idx_hg + 309);

    auto ts_zzzzz_yyyy = pbuffer.data(idx_hg + 310);

    auto ts_zzzzz_yyyz = pbuffer.data(idx_hg + 311);

    auto ts_zzzzz_yyzz = pbuffer.data(idx_hg + 312);

    auto ts_zzzzz_yzzz = pbuffer.data(idx_hg + 313);

    auto ts_zzzzz_zzzz = pbuffer.data(idx_hg + 314);

    // Set up 0-15 components of targeted buffer : HG

    auto gs_x_xxxxx_xxxx = pbuffer.data(idx_g_hg);

    auto gs_x_xxxxx_xxxy = pbuffer.data(idx_g_hg + 1);

    auto gs_x_xxxxx_xxxz = pbuffer.data(idx_g_hg + 2);

    auto gs_x_xxxxx_xxyy = pbuffer.data(idx_g_hg + 3);

    auto gs_x_xxxxx_xxyz = pbuffer.data(idx_g_hg + 4);

    auto gs_x_xxxxx_xxzz = pbuffer.data(idx_g_hg + 5);

    auto gs_x_xxxxx_xyyy = pbuffer.data(idx_g_hg + 6);

    auto gs_x_xxxxx_xyyz = pbuffer.data(idx_g_hg + 7);

    auto gs_x_xxxxx_xyzz = pbuffer.data(idx_g_hg + 8);

    auto gs_x_xxxxx_xzzz = pbuffer.data(idx_g_hg + 9);

    auto gs_x_xxxxx_yyyy = pbuffer.data(idx_g_hg + 10);

    auto gs_x_xxxxx_yyyz = pbuffer.data(idx_g_hg + 11);

    auto gs_x_xxxxx_yyzz = pbuffer.data(idx_g_hg + 12);

    auto gs_x_xxxxx_yzzz = pbuffer.data(idx_g_hg + 13);

    auto gs_x_xxxxx_zzzz = pbuffer.data(idx_g_hg + 14);

    #pragma omp simd aligned(gc_x, gs_x_xxxxx_xxxx, gs_x_xxxxx_xxxy, gs_x_xxxxx_xxxz, gs_x_xxxxx_xxyy, gs_x_xxxxx_xxyz, gs_x_xxxxx_xxzz, gs_x_xxxxx_xyyy, gs_x_xxxxx_xyyz, gs_x_xxxxx_xyzz, gs_x_xxxxx_xzzz, gs_x_xxxxx_yyyy, gs_x_xxxxx_yyyz, gs_x_xxxxx_yyzz, gs_x_xxxxx_yzzz, gs_x_xxxxx_zzzz, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxzz, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyzz, ts_xxxx_xzzz, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyzz, ts_xxxx_yzzz, ts_xxxx_zzzz, ts_xxxxx_xxx, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxy, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxz, ts_xxxxx_xxzz, ts_xxxxx_xyy, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyz, ts_xxxxx_xyzz, ts_xxxxx_xzz, ts_xxxxx_xzzz, ts_xxxxx_yyy, ts_xxxxx_yyyy, ts_xxxxx_yyyz, ts_xxxxx_yyz, ts_xxxxx_yyzz, ts_xxxxx_yzz, ts_xxxxx_yzzz, ts_xxxxx_zzz, ts_xxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxx_xxxx[i] = 10.0 * ts_xxxx_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxy[i] = 10.0 * ts_xxxx_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxz[i] = 10.0 * ts_xxxx_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxyy[i] = 10.0 * ts_xxxx_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxyz[i] = 10.0 * ts_xxxx_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxzz[i] = 10.0 * ts_xxxx_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyyy[i] = 10.0 * ts_xxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyyz[i] = 10.0 * ts_xxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyzz[i] = 10.0 * ts_xxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xzzz[i] = 10.0 * ts_xxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyyy[i] = 10.0 * ts_xxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyyz[i] = 10.0 * ts_xxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyzz[i] = 10.0 * ts_xxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yzzz[i] = 10.0 * ts_xxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_zzzz[i] = 10.0 * ts_xxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 15-30 components of targeted buffer : HG

    auto gs_x_xxxxy_xxxx = pbuffer.data(idx_g_hg + 15);

    auto gs_x_xxxxy_xxxy = pbuffer.data(idx_g_hg + 16);

    auto gs_x_xxxxy_xxxz = pbuffer.data(idx_g_hg + 17);

    auto gs_x_xxxxy_xxyy = pbuffer.data(idx_g_hg + 18);

    auto gs_x_xxxxy_xxyz = pbuffer.data(idx_g_hg + 19);

    auto gs_x_xxxxy_xxzz = pbuffer.data(idx_g_hg + 20);

    auto gs_x_xxxxy_xyyy = pbuffer.data(idx_g_hg + 21);

    auto gs_x_xxxxy_xyyz = pbuffer.data(idx_g_hg + 22);

    auto gs_x_xxxxy_xyzz = pbuffer.data(idx_g_hg + 23);

    auto gs_x_xxxxy_xzzz = pbuffer.data(idx_g_hg + 24);

    auto gs_x_xxxxy_yyyy = pbuffer.data(idx_g_hg + 25);

    auto gs_x_xxxxy_yyyz = pbuffer.data(idx_g_hg + 26);

    auto gs_x_xxxxy_yyzz = pbuffer.data(idx_g_hg + 27);

    auto gs_x_xxxxy_yzzz = pbuffer.data(idx_g_hg + 28);

    auto gs_x_xxxxy_zzzz = pbuffer.data(idx_g_hg + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxxxy_xxxx, gs_x_xxxxy_xxxy, gs_x_xxxxy_xxxz, gs_x_xxxxy_xxyy, gs_x_xxxxy_xxyz, gs_x_xxxxy_xxzz, gs_x_xxxxy_xyyy, gs_x_xxxxy_xyyz, gs_x_xxxxy_xyzz, gs_x_xxxxy_xzzz, gs_x_xxxxy_yyyy, gs_x_xxxxy_yyyz, gs_x_xxxxy_yyzz, gs_x_xxxxy_yzzz, gs_x_xxxxy_zzzz, ts_xxxxy_xxx, ts_xxxxy_xxxx, ts_xxxxy_xxxy, ts_xxxxy_xxxz, ts_xxxxy_xxy, ts_xxxxy_xxyy, ts_xxxxy_xxyz, ts_xxxxy_xxz, ts_xxxxy_xxzz, ts_xxxxy_xyy, ts_xxxxy_xyyy, ts_xxxxy_xyyz, ts_xxxxy_xyz, ts_xxxxy_xyzz, ts_xxxxy_xzz, ts_xxxxy_xzzz, ts_xxxxy_yyy, ts_xxxxy_yyyy, ts_xxxxy_yyyz, ts_xxxxy_yyz, ts_xxxxy_yyzz, ts_xxxxy_yzz, ts_xxxxy_yzzz, ts_xxxxy_zzz, ts_xxxxy_zzzz, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxzz, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyzz, ts_xxxy_xzzz, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyzz, ts_xxxy_yzzz, ts_xxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxy_xxxx[i] = 8.0 * ts_xxxy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxy[i] = 8.0 * ts_xxxy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxz[i] = 8.0 * ts_xxxy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxyy[i] = 8.0 * ts_xxxy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxyz[i] = 8.0 * ts_xxxy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxzz[i] = 8.0 * ts_xxxy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyyy[i] = 8.0 * ts_xxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyyz[i] = 8.0 * ts_xxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyzz[i] = 8.0 * ts_xxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xzzz[i] = 8.0 * ts_xxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyyy[i] = 8.0 * ts_xxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyyz[i] = 8.0 * ts_xxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyzz[i] = 8.0 * ts_xxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yzzz[i] = 8.0 * ts_xxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_zzzz[i] = 8.0 * ts_xxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-45 components of targeted buffer : HG

    auto gs_x_xxxxz_xxxx = pbuffer.data(idx_g_hg + 30);

    auto gs_x_xxxxz_xxxy = pbuffer.data(idx_g_hg + 31);

    auto gs_x_xxxxz_xxxz = pbuffer.data(idx_g_hg + 32);

    auto gs_x_xxxxz_xxyy = pbuffer.data(idx_g_hg + 33);

    auto gs_x_xxxxz_xxyz = pbuffer.data(idx_g_hg + 34);

    auto gs_x_xxxxz_xxzz = pbuffer.data(idx_g_hg + 35);

    auto gs_x_xxxxz_xyyy = pbuffer.data(idx_g_hg + 36);

    auto gs_x_xxxxz_xyyz = pbuffer.data(idx_g_hg + 37);

    auto gs_x_xxxxz_xyzz = pbuffer.data(idx_g_hg + 38);

    auto gs_x_xxxxz_xzzz = pbuffer.data(idx_g_hg + 39);

    auto gs_x_xxxxz_yyyy = pbuffer.data(idx_g_hg + 40);

    auto gs_x_xxxxz_yyyz = pbuffer.data(idx_g_hg + 41);

    auto gs_x_xxxxz_yyzz = pbuffer.data(idx_g_hg + 42);

    auto gs_x_xxxxz_yzzz = pbuffer.data(idx_g_hg + 43);

    auto gs_x_xxxxz_zzzz = pbuffer.data(idx_g_hg + 44);

    #pragma omp simd aligned(gc_x, gs_x_xxxxz_xxxx, gs_x_xxxxz_xxxy, gs_x_xxxxz_xxxz, gs_x_xxxxz_xxyy, gs_x_xxxxz_xxyz, gs_x_xxxxz_xxzz, gs_x_xxxxz_xyyy, gs_x_xxxxz_xyyz, gs_x_xxxxz_xyzz, gs_x_xxxxz_xzzz, gs_x_xxxxz_yyyy, gs_x_xxxxz_yyyz, gs_x_xxxxz_yyzz, gs_x_xxxxz_yzzz, gs_x_xxxxz_zzzz, ts_xxxxz_xxx, ts_xxxxz_xxxx, ts_xxxxz_xxxy, ts_xxxxz_xxxz, ts_xxxxz_xxy, ts_xxxxz_xxyy, ts_xxxxz_xxyz, ts_xxxxz_xxz, ts_xxxxz_xxzz, ts_xxxxz_xyy, ts_xxxxz_xyyy, ts_xxxxz_xyyz, ts_xxxxz_xyz, ts_xxxxz_xyzz, ts_xxxxz_xzz, ts_xxxxz_xzzz, ts_xxxxz_yyy, ts_xxxxz_yyyy, ts_xxxxz_yyyz, ts_xxxxz_yyz, ts_xxxxz_yyzz, ts_xxxxz_yzz, ts_xxxxz_yzzz, ts_xxxxz_zzz, ts_xxxxz_zzzz, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxzz, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyzz, ts_xxxz_xzzz, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyzz, ts_xxxz_yzzz, ts_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxz_xxxx[i] = 8.0 * ts_xxxz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxy[i] = 8.0 * ts_xxxz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxz[i] = 8.0 * ts_xxxz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxyy[i] = 8.0 * ts_xxxz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxyz[i] = 8.0 * ts_xxxz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxzz[i] = 8.0 * ts_xxxz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyyy[i] = 8.0 * ts_xxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyyz[i] = 8.0 * ts_xxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyzz[i] = 8.0 * ts_xxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xzzz[i] = 8.0 * ts_xxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyyy[i] = 8.0 * ts_xxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyyz[i] = 8.0 * ts_xxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyzz[i] = 8.0 * ts_xxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yzzz[i] = 8.0 * ts_xxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_zzzz[i] = 8.0 * ts_xxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 45-60 components of targeted buffer : HG

    auto gs_x_xxxyy_xxxx = pbuffer.data(idx_g_hg + 45);

    auto gs_x_xxxyy_xxxy = pbuffer.data(idx_g_hg + 46);

    auto gs_x_xxxyy_xxxz = pbuffer.data(idx_g_hg + 47);

    auto gs_x_xxxyy_xxyy = pbuffer.data(idx_g_hg + 48);

    auto gs_x_xxxyy_xxyz = pbuffer.data(idx_g_hg + 49);

    auto gs_x_xxxyy_xxzz = pbuffer.data(idx_g_hg + 50);

    auto gs_x_xxxyy_xyyy = pbuffer.data(idx_g_hg + 51);

    auto gs_x_xxxyy_xyyz = pbuffer.data(idx_g_hg + 52);

    auto gs_x_xxxyy_xyzz = pbuffer.data(idx_g_hg + 53);

    auto gs_x_xxxyy_xzzz = pbuffer.data(idx_g_hg + 54);

    auto gs_x_xxxyy_yyyy = pbuffer.data(idx_g_hg + 55);

    auto gs_x_xxxyy_yyyz = pbuffer.data(idx_g_hg + 56);

    auto gs_x_xxxyy_yyzz = pbuffer.data(idx_g_hg + 57);

    auto gs_x_xxxyy_yzzz = pbuffer.data(idx_g_hg + 58);

    auto gs_x_xxxyy_zzzz = pbuffer.data(idx_g_hg + 59);

    #pragma omp simd aligned(gc_x, gs_x_xxxyy_xxxx, gs_x_xxxyy_xxxy, gs_x_xxxyy_xxxz, gs_x_xxxyy_xxyy, gs_x_xxxyy_xxyz, gs_x_xxxyy_xxzz, gs_x_xxxyy_xyyy, gs_x_xxxyy_xyyz, gs_x_xxxyy_xyzz, gs_x_xxxyy_xzzz, gs_x_xxxyy_yyyy, gs_x_xxxyy_yyyz, gs_x_xxxyy_yyzz, gs_x_xxxyy_yzzz, gs_x_xxxyy_zzzz, ts_xxxyy_xxx, ts_xxxyy_xxxx, ts_xxxyy_xxxy, ts_xxxyy_xxxz, ts_xxxyy_xxy, ts_xxxyy_xxyy, ts_xxxyy_xxyz, ts_xxxyy_xxz, ts_xxxyy_xxzz, ts_xxxyy_xyy, ts_xxxyy_xyyy, ts_xxxyy_xyyz, ts_xxxyy_xyz, ts_xxxyy_xyzz, ts_xxxyy_xzz, ts_xxxyy_xzzz, ts_xxxyy_yyy, ts_xxxyy_yyyy, ts_xxxyy_yyyz, ts_xxxyy_yyz, ts_xxxyy_yyzz, ts_xxxyy_yzz, ts_xxxyy_yzzz, ts_xxxyy_zzz, ts_xxxyy_zzzz, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxzz, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyzz, ts_xxyy_xzzz, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyzz, ts_xxyy_yzzz, ts_xxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyy_xxxx[i] = 6.0 * ts_xxyy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxy[i] = 6.0 * ts_xxyy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxz[i] = 6.0 * ts_xxyy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxyy[i] = 6.0 * ts_xxyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxyz[i] = 6.0 * ts_xxyy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxzz[i] = 6.0 * ts_xxyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyyy[i] = 6.0 * ts_xxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyyz[i] = 6.0 * ts_xxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyzz[i] = 6.0 * ts_xxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xzzz[i] = 6.0 * ts_xxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyyy[i] = 6.0 * ts_xxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyyz[i] = 6.0 * ts_xxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyzz[i] = 6.0 * ts_xxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yzzz[i] = 6.0 * ts_xxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_zzzz[i] = 6.0 * ts_xxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-75 components of targeted buffer : HG

    auto gs_x_xxxyz_xxxx = pbuffer.data(idx_g_hg + 60);

    auto gs_x_xxxyz_xxxy = pbuffer.data(idx_g_hg + 61);

    auto gs_x_xxxyz_xxxz = pbuffer.data(idx_g_hg + 62);

    auto gs_x_xxxyz_xxyy = pbuffer.data(idx_g_hg + 63);

    auto gs_x_xxxyz_xxyz = pbuffer.data(idx_g_hg + 64);

    auto gs_x_xxxyz_xxzz = pbuffer.data(idx_g_hg + 65);

    auto gs_x_xxxyz_xyyy = pbuffer.data(idx_g_hg + 66);

    auto gs_x_xxxyz_xyyz = pbuffer.data(idx_g_hg + 67);

    auto gs_x_xxxyz_xyzz = pbuffer.data(idx_g_hg + 68);

    auto gs_x_xxxyz_xzzz = pbuffer.data(idx_g_hg + 69);

    auto gs_x_xxxyz_yyyy = pbuffer.data(idx_g_hg + 70);

    auto gs_x_xxxyz_yyyz = pbuffer.data(idx_g_hg + 71);

    auto gs_x_xxxyz_yyzz = pbuffer.data(idx_g_hg + 72);

    auto gs_x_xxxyz_yzzz = pbuffer.data(idx_g_hg + 73);

    auto gs_x_xxxyz_zzzz = pbuffer.data(idx_g_hg + 74);

    #pragma omp simd aligned(gc_x, gs_x_xxxyz_xxxx, gs_x_xxxyz_xxxy, gs_x_xxxyz_xxxz, gs_x_xxxyz_xxyy, gs_x_xxxyz_xxyz, gs_x_xxxyz_xxzz, gs_x_xxxyz_xyyy, gs_x_xxxyz_xyyz, gs_x_xxxyz_xyzz, gs_x_xxxyz_xzzz, gs_x_xxxyz_yyyy, gs_x_xxxyz_yyyz, gs_x_xxxyz_yyzz, gs_x_xxxyz_yzzz, gs_x_xxxyz_zzzz, ts_xxxyz_xxx, ts_xxxyz_xxxx, ts_xxxyz_xxxy, ts_xxxyz_xxxz, ts_xxxyz_xxy, ts_xxxyz_xxyy, ts_xxxyz_xxyz, ts_xxxyz_xxz, ts_xxxyz_xxzz, ts_xxxyz_xyy, ts_xxxyz_xyyy, ts_xxxyz_xyyz, ts_xxxyz_xyz, ts_xxxyz_xyzz, ts_xxxyz_xzz, ts_xxxyz_xzzz, ts_xxxyz_yyy, ts_xxxyz_yyyy, ts_xxxyz_yyyz, ts_xxxyz_yyz, ts_xxxyz_yyzz, ts_xxxyz_yzz, ts_xxxyz_yzzz, ts_xxxyz_zzz, ts_xxxyz_zzzz, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxzz, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyzz, ts_xxyz_xzzz, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyzz, ts_xxyz_yzzz, ts_xxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyz_xxxx[i] = 6.0 * ts_xxyz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxy[i] = 6.0 * ts_xxyz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxz[i] = 6.0 * ts_xxyz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxyy[i] = 6.0 * ts_xxyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxyz[i] = 6.0 * ts_xxyz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxzz[i] = 6.0 * ts_xxyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyyy[i] = 6.0 * ts_xxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyyz[i] = 6.0 * ts_xxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyzz[i] = 6.0 * ts_xxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xzzz[i] = 6.0 * ts_xxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyyy[i] = 6.0 * ts_xxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyyz[i] = 6.0 * ts_xxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyzz[i] = 6.0 * ts_xxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yzzz[i] = 6.0 * ts_xxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_zzzz[i] = 6.0 * ts_xxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 75-90 components of targeted buffer : HG

    auto gs_x_xxxzz_xxxx = pbuffer.data(idx_g_hg + 75);

    auto gs_x_xxxzz_xxxy = pbuffer.data(idx_g_hg + 76);

    auto gs_x_xxxzz_xxxz = pbuffer.data(idx_g_hg + 77);

    auto gs_x_xxxzz_xxyy = pbuffer.data(idx_g_hg + 78);

    auto gs_x_xxxzz_xxyz = pbuffer.data(idx_g_hg + 79);

    auto gs_x_xxxzz_xxzz = pbuffer.data(idx_g_hg + 80);

    auto gs_x_xxxzz_xyyy = pbuffer.data(idx_g_hg + 81);

    auto gs_x_xxxzz_xyyz = pbuffer.data(idx_g_hg + 82);

    auto gs_x_xxxzz_xyzz = pbuffer.data(idx_g_hg + 83);

    auto gs_x_xxxzz_xzzz = pbuffer.data(idx_g_hg + 84);

    auto gs_x_xxxzz_yyyy = pbuffer.data(idx_g_hg + 85);

    auto gs_x_xxxzz_yyyz = pbuffer.data(idx_g_hg + 86);

    auto gs_x_xxxzz_yyzz = pbuffer.data(idx_g_hg + 87);

    auto gs_x_xxxzz_yzzz = pbuffer.data(idx_g_hg + 88);

    auto gs_x_xxxzz_zzzz = pbuffer.data(idx_g_hg + 89);

    #pragma omp simd aligned(gc_x, gs_x_xxxzz_xxxx, gs_x_xxxzz_xxxy, gs_x_xxxzz_xxxz, gs_x_xxxzz_xxyy, gs_x_xxxzz_xxyz, gs_x_xxxzz_xxzz, gs_x_xxxzz_xyyy, gs_x_xxxzz_xyyz, gs_x_xxxzz_xyzz, gs_x_xxxzz_xzzz, gs_x_xxxzz_yyyy, gs_x_xxxzz_yyyz, gs_x_xxxzz_yyzz, gs_x_xxxzz_yzzz, gs_x_xxxzz_zzzz, ts_xxxzz_xxx, ts_xxxzz_xxxx, ts_xxxzz_xxxy, ts_xxxzz_xxxz, ts_xxxzz_xxy, ts_xxxzz_xxyy, ts_xxxzz_xxyz, ts_xxxzz_xxz, ts_xxxzz_xxzz, ts_xxxzz_xyy, ts_xxxzz_xyyy, ts_xxxzz_xyyz, ts_xxxzz_xyz, ts_xxxzz_xyzz, ts_xxxzz_xzz, ts_xxxzz_xzzz, ts_xxxzz_yyy, ts_xxxzz_yyyy, ts_xxxzz_yyyz, ts_xxxzz_yyz, ts_xxxzz_yyzz, ts_xxxzz_yzz, ts_xxxzz_yzzz, ts_xxxzz_zzz, ts_xxxzz_zzzz, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxzz, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyzz, ts_xxzz_xzzz, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyzz, ts_xxzz_yzzz, ts_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxzz_xxxx[i] = 6.0 * ts_xxzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxy[i] = 6.0 * ts_xxzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxz[i] = 6.0 * ts_xxzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxyy[i] = 6.0 * ts_xxzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxyz[i] = 6.0 * ts_xxzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxzz[i] = 6.0 * ts_xxzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyyy[i] = 6.0 * ts_xxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyyz[i] = 6.0 * ts_xxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyzz[i] = 6.0 * ts_xxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xzzz[i] = 6.0 * ts_xxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyyy[i] = 6.0 * ts_xxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyyz[i] = 6.0 * ts_xxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyzz[i] = 6.0 * ts_xxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yzzz[i] = 6.0 * ts_xxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_zzzz[i] = 6.0 * ts_xxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-105 components of targeted buffer : HG

    auto gs_x_xxyyy_xxxx = pbuffer.data(idx_g_hg + 90);

    auto gs_x_xxyyy_xxxy = pbuffer.data(idx_g_hg + 91);

    auto gs_x_xxyyy_xxxz = pbuffer.data(idx_g_hg + 92);

    auto gs_x_xxyyy_xxyy = pbuffer.data(idx_g_hg + 93);

    auto gs_x_xxyyy_xxyz = pbuffer.data(idx_g_hg + 94);

    auto gs_x_xxyyy_xxzz = pbuffer.data(idx_g_hg + 95);

    auto gs_x_xxyyy_xyyy = pbuffer.data(idx_g_hg + 96);

    auto gs_x_xxyyy_xyyz = pbuffer.data(idx_g_hg + 97);

    auto gs_x_xxyyy_xyzz = pbuffer.data(idx_g_hg + 98);

    auto gs_x_xxyyy_xzzz = pbuffer.data(idx_g_hg + 99);

    auto gs_x_xxyyy_yyyy = pbuffer.data(idx_g_hg + 100);

    auto gs_x_xxyyy_yyyz = pbuffer.data(idx_g_hg + 101);

    auto gs_x_xxyyy_yyzz = pbuffer.data(idx_g_hg + 102);

    auto gs_x_xxyyy_yzzz = pbuffer.data(idx_g_hg + 103);

    auto gs_x_xxyyy_zzzz = pbuffer.data(idx_g_hg + 104);

    #pragma omp simd aligned(gc_x, gs_x_xxyyy_xxxx, gs_x_xxyyy_xxxy, gs_x_xxyyy_xxxz, gs_x_xxyyy_xxyy, gs_x_xxyyy_xxyz, gs_x_xxyyy_xxzz, gs_x_xxyyy_xyyy, gs_x_xxyyy_xyyz, gs_x_xxyyy_xyzz, gs_x_xxyyy_xzzz, gs_x_xxyyy_yyyy, gs_x_xxyyy_yyyz, gs_x_xxyyy_yyzz, gs_x_xxyyy_yzzz, gs_x_xxyyy_zzzz, ts_xxyyy_xxx, ts_xxyyy_xxxx, ts_xxyyy_xxxy, ts_xxyyy_xxxz, ts_xxyyy_xxy, ts_xxyyy_xxyy, ts_xxyyy_xxyz, ts_xxyyy_xxz, ts_xxyyy_xxzz, ts_xxyyy_xyy, ts_xxyyy_xyyy, ts_xxyyy_xyyz, ts_xxyyy_xyz, ts_xxyyy_xyzz, ts_xxyyy_xzz, ts_xxyyy_xzzz, ts_xxyyy_yyy, ts_xxyyy_yyyy, ts_xxyyy_yyyz, ts_xxyyy_yyz, ts_xxyyy_yyzz, ts_xxyyy_yzz, ts_xxyyy_yzzz, ts_xxyyy_zzz, ts_xxyyy_zzzz, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxzz, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyzz, ts_xyyy_xzzz, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyzz, ts_xyyy_yzzz, ts_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyy_xxxx[i] = 4.0 * ts_xyyy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxy[i] = 4.0 * ts_xyyy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxz[i] = 4.0 * ts_xyyy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxyy[i] = 4.0 * ts_xyyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxyz[i] = 4.0 * ts_xyyy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxzz[i] = 4.0 * ts_xyyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyyy[i] = 4.0 * ts_xyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyyz[i] = 4.0 * ts_xyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyzz[i] = 4.0 * ts_xyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xzzz[i] = 4.0 * ts_xyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyyy[i] = 4.0 * ts_xyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyyz[i] = 4.0 * ts_xyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyzz[i] = 4.0 * ts_xyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yzzz[i] = 4.0 * ts_xyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_zzzz[i] = 4.0 * ts_xyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-120 components of targeted buffer : HG

    auto gs_x_xxyyz_xxxx = pbuffer.data(idx_g_hg + 105);

    auto gs_x_xxyyz_xxxy = pbuffer.data(idx_g_hg + 106);

    auto gs_x_xxyyz_xxxz = pbuffer.data(idx_g_hg + 107);

    auto gs_x_xxyyz_xxyy = pbuffer.data(idx_g_hg + 108);

    auto gs_x_xxyyz_xxyz = pbuffer.data(idx_g_hg + 109);

    auto gs_x_xxyyz_xxzz = pbuffer.data(idx_g_hg + 110);

    auto gs_x_xxyyz_xyyy = pbuffer.data(idx_g_hg + 111);

    auto gs_x_xxyyz_xyyz = pbuffer.data(idx_g_hg + 112);

    auto gs_x_xxyyz_xyzz = pbuffer.data(idx_g_hg + 113);

    auto gs_x_xxyyz_xzzz = pbuffer.data(idx_g_hg + 114);

    auto gs_x_xxyyz_yyyy = pbuffer.data(idx_g_hg + 115);

    auto gs_x_xxyyz_yyyz = pbuffer.data(idx_g_hg + 116);

    auto gs_x_xxyyz_yyzz = pbuffer.data(idx_g_hg + 117);

    auto gs_x_xxyyz_yzzz = pbuffer.data(idx_g_hg + 118);

    auto gs_x_xxyyz_zzzz = pbuffer.data(idx_g_hg + 119);

    #pragma omp simd aligned(gc_x, gs_x_xxyyz_xxxx, gs_x_xxyyz_xxxy, gs_x_xxyyz_xxxz, gs_x_xxyyz_xxyy, gs_x_xxyyz_xxyz, gs_x_xxyyz_xxzz, gs_x_xxyyz_xyyy, gs_x_xxyyz_xyyz, gs_x_xxyyz_xyzz, gs_x_xxyyz_xzzz, gs_x_xxyyz_yyyy, gs_x_xxyyz_yyyz, gs_x_xxyyz_yyzz, gs_x_xxyyz_yzzz, gs_x_xxyyz_zzzz, ts_xxyyz_xxx, ts_xxyyz_xxxx, ts_xxyyz_xxxy, ts_xxyyz_xxxz, ts_xxyyz_xxy, ts_xxyyz_xxyy, ts_xxyyz_xxyz, ts_xxyyz_xxz, ts_xxyyz_xxzz, ts_xxyyz_xyy, ts_xxyyz_xyyy, ts_xxyyz_xyyz, ts_xxyyz_xyz, ts_xxyyz_xyzz, ts_xxyyz_xzz, ts_xxyyz_xzzz, ts_xxyyz_yyy, ts_xxyyz_yyyy, ts_xxyyz_yyyz, ts_xxyyz_yyz, ts_xxyyz_yyzz, ts_xxyyz_yzz, ts_xxyyz_yzzz, ts_xxyyz_zzz, ts_xxyyz_zzzz, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxzz, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyzz, ts_xyyz_xzzz, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyzz, ts_xyyz_yzzz, ts_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyz_xxxx[i] = 4.0 * ts_xyyz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxy[i] = 4.0 * ts_xyyz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxz[i] = 4.0 * ts_xyyz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxyy[i] = 4.0 * ts_xyyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxyz[i] = 4.0 * ts_xyyz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxzz[i] = 4.0 * ts_xyyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyyy[i] = 4.0 * ts_xyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyyz[i] = 4.0 * ts_xyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyzz[i] = 4.0 * ts_xyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xzzz[i] = 4.0 * ts_xyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyyy[i] = 4.0 * ts_xyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyyz[i] = 4.0 * ts_xyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyzz[i] = 4.0 * ts_xyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yzzz[i] = 4.0 * ts_xyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_zzzz[i] = 4.0 * ts_xyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 120-135 components of targeted buffer : HG

    auto gs_x_xxyzz_xxxx = pbuffer.data(idx_g_hg + 120);

    auto gs_x_xxyzz_xxxy = pbuffer.data(idx_g_hg + 121);

    auto gs_x_xxyzz_xxxz = pbuffer.data(idx_g_hg + 122);

    auto gs_x_xxyzz_xxyy = pbuffer.data(idx_g_hg + 123);

    auto gs_x_xxyzz_xxyz = pbuffer.data(idx_g_hg + 124);

    auto gs_x_xxyzz_xxzz = pbuffer.data(idx_g_hg + 125);

    auto gs_x_xxyzz_xyyy = pbuffer.data(idx_g_hg + 126);

    auto gs_x_xxyzz_xyyz = pbuffer.data(idx_g_hg + 127);

    auto gs_x_xxyzz_xyzz = pbuffer.data(idx_g_hg + 128);

    auto gs_x_xxyzz_xzzz = pbuffer.data(idx_g_hg + 129);

    auto gs_x_xxyzz_yyyy = pbuffer.data(idx_g_hg + 130);

    auto gs_x_xxyzz_yyyz = pbuffer.data(idx_g_hg + 131);

    auto gs_x_xxyzz_yyzz = pbuffer.data(idx_g_hg + 132);

    auto gs_x_xxyzz_yzzz = pbuffer.data(idx_g_hg + 133);

    auto gs_x_xxyzz_zzzz = pbuffer.data(idx_g_hg + 134);

    #pragma omp simd aligned(gc_x, gs_x_xxyzz_xxxx, gs_x_xxyzz_xxxy, gs_x_xxyzz_xxxz, gs_x_xxyzz_xxyy, gs_x_xxyzz_xxyz, gs_x_xxyzz_xxzz, gs_x_xxyzz_xyyy, gs_x_xxyzz_xyyz, gs_x_xxyzz_xyzz, gs_x_xxyzz_xzzz, gs_x_xxyzz_yyyy, gs_x_xxyzz_yyyz, gs_x_xxyzz_yyzz, gs_x_xxyzz_yzzz, gs_x_xxyzz_zzzz, ts_xxyzz_xxx, ts_xxyzz_xxxx, ts_xxyzz_xxxy, ts_xxyzz_xxxz, ts_xxyzz_xxy, ts_xxyzz_xxyy, ts_xxyzz_xxyz, ts_xxyzz_xxz, ts_xxyzz_xxzz, ts_xxyzz_xyy, ts_xxyzz_xyyy, ts_xxyzz_xyyz, ts_xxyzz_xyz, ts_xxyzz_xyzz, ts_xxyzz_xzz, ts_xxyzz_xzzz, ts_xxyzz_yyy, ts_xxyzz_yyyy, ts_xxyzz_yyyz, ts_xxyzz_yyz, ts_xxyzz_yyzz, ts_xxyzz_yzz, ts_xxyzz_yzzz, ts_xxyzz_zzz, ts_xxyzz_zzzz, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxzz, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyzz, ts_xyzz_xzzz, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyzz, ts_xyzz_yzzz, ts_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyzz_xxxx[i] = 4.0 * ts_xyzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxy[i] = 4.0 * ts_xyzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxz[i] = 4.0 * ts_xyzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxyy[i] = 4.0 * ts_xyzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxyz[i] = 4.0 * ts_xyzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxzz[i] = 4.0 * ts_xyzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyyy[i] = 4.0 * ts_xyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyyz[i] = 4.0 * ts_xyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyzz[i] = 4.0 * ts_xyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xzzz[i] = 4.0 * ts_xyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyyy[i] = 4.0 * ts_xyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyyz[i] = 4.0 * ts_xyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyzz[i] = 4.0 * ts_xyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yzzz[i] = 4.0 * ts_xyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_zzzz[i] = 4.0 * ts_xyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 135-150 components of targeted buffer : HG

    auto gs_x_xxzzz_xxxx = pbuffer.data(idx_g_hg + 135);

    auto gs_x_xxzzz_xxxy = pbuffer.data(idx_g_hg + 136);

    auto gs_x_xxzzz_xxxz = pbuffer.data(idx_g_hg + 137);

    auto gs_x_xxzzz_xxyy = pbuffer.data(idx_g_hg + 138);

    auto gs_x_xxzzz_xxyz = pbuffer.data(idx_g_hg + 139);

    auto gs_x_xxzzz_xxzz = pbuffer.data(idx_g_hg + 140);

    auto gs_x_xxzzz_xyyy = pbuffer.data(idx_g_hg + 141);

    auto gs_x_xxzzz_xyyz = pbuffer.data(idx_g_hg + 142);

    auto gs_x_xxzzz_xyzz = pbuffer.data(idx_g_hg + 143);

    auto gs_x_xxzzz_xzzz = pbuffer.data(idx_g_hg + 144);

    auto gs_x_xxzzz_yyyy = pbuffer.data(idx_g_hg + 145);

    auto gs_x_xxzzz_yyyz = pbuffer.data(idx_g_hg + 146);

    auto gs_x_xxzzz_yyzz = pbuffer.data(idx_g_hg + 147);

    auto gs_x_xxzzz_yzzz = pbuffer.data(idx_g_hg + 148);

    auto gs_x_xxzzz_zzzz = pbuffer.data(idx_g_hg + 149);

    #pragma omp simd aligned(gc_x, gs_x_xxzzz_xxxx, gs_x_xxzzz_xxxy, gs_x_xxzzz_xxxz, gs_x_xxzzz_xxyy, gs_x_xxzzz_xxyz, gs_x_xxzzz_xxzz, gs_x_xxzzz_xyyy, gs_x_xxzzz_xyyz, gs_x_xxzzz_xyzz, gs_x_xxzzz_xzzz, gs_x_xxzzz_yyyy, gs_x_xxzzz_yyyz, gs_x_xxzzz_yyzz, gs_x_xxzzz_yzzz, gs_x_xxzzz_zzzz, ts_xxzzz_xxx, ts_xxzzz_xxxx, ts_xxzzz_xxxy, ts_xxzzz_xxxz, ts_xxzzz_xxy, ts_xxzzz_xxyy, ts_xxzzz_xxyz, ts_xxzzz_xxz, ts_xxzzz_xxzz, ts_xxzzz_xyy, ts_xxzzz_xyyy, ts_xxzzz_xyyz, ts_xxzzz_xyz, ts_xxzzz_xyzz, ts_xxzzz_xzz, ts_xxzzz_xzzz, ts_xxzzz_yyy, ts_xxzzz_yyyy, ts_xxzzz_yyyz, ts_xxzzz_yyz, ts_xxzzz_yyzz, ts_xxzzz_yzz, ts_xxzzz_yzzz, ts_xxzzz_zzz, ts_xxzzz_zzzz, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxzz, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyzz, ts_xzzz_xzzz, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyzz, ts_xzzz_yzzz, ts_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzzz_xxxx[i] = 4.0 * ts_xzzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxy[i] = 4.0 * ts_xzzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxz[i] = 4.0 * ts_xzzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxyy[i] = 4.0 * ts_xzzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxyz[i] = 4.0 * ts_xzzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxzz[i] = 4.0 * ts_xzzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyyy[i] = 4.0 * ts_xzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyyz[i] = 4.0 * ts_xzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyzz[i] = 4.0 * ts_xzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xzzz[i] = 4.0 * ts_xzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyyy[i] = 4.0 * ts_xzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyyz[i] = 4.0 * ts_xzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyzz[i] = 4.0 * ts_xzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yzzz[i] = 4.0 * ts_xzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_zzzz[i] = 4.0 * ts_xzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 150-165 components of targeted buffer : HG

    auto gs_x_xyyyy_xxxx = pbuffer.data(idx_g_hg + 150);

    auto gs_x_xyyyy_xxxy = pbuffer.data(idx_g_hg + 151);

    auto gs_x_xyyyy_xxxz = pbuffer.data(idx_g_hg + 152);

    auto gs_x_xyyyy_xxyy = pbuffer.data(idx_g_hg + 153);

    auto gs_x_xyyyy_xxyz = pbuffer.data(idx_g_hg + 154);

    auto gs_x_xyyyy_xxzz = pbuffer.data(idx_g_hg + 155);

    auto gs_x_xyyyy_xyyy = pbuffer.data(idx_g_hg + 156);

    auto gs_x_xyyyy_xyyz = pbuffer.data(idx_g_hg + 157);

    auto gs_x_xyyyy_xyzz = pbuffer.data(idx_g_hg + 158);

    auto gs_x_xyyyy_xzzz = pbuffer.data(idx_g_hg + 159);

    auto gs_x_xyyyy_yyyy = pbuffer.data(idx_g_hg + 160);

    auto gs_x_xyyyy_yyyz = pbuffer.data(idx_g_hg + 161);

    auto gs_x_xyyyy_yyzz = pbuffer.data(idx_g_hg + 162);

    auto gs_x_xyyyy_yzzz = pbuffer.data(idx_g_hg + 163);

    auto gs_x_xyyyy_zzzz = pbuffer.data(idx_g_hg + 164);

    #pragma omp simd aligned(gc_x, gs_x_xyyyy_xxxx, gs_x_xyyyy_xxxy, gs_x_xyyyy_xxxz, gs_x_xyyyy_xxyy, gs_x_xyyyy_xxyz, gs_x_xyyyy_xxzz, gs_x_xyyyy_xyyy, gs_x_xyyyy_xyyz, gs_x_xyyyy_xyzz, gs_x_xyyyy_xzzz, gs_x_xyyyy_yyyy, gs_x_xyyyy_yyyz, gs_x_xyyyy_yyzz, gs_x_xyyyy_yzzz, gs_x_xyyyy_zzzz, ts_xyyyy_xxx, ts_xyyyy_xxxx, ts_xyyyy_xxxy, ts_xyyyy_xxxz, ts_xyyyy_xxy, ts_xyyyy_xxyy, ts_xyyyy_xxyz, ts_xyyyy_xxz, ts_xyyyy_xxzz, ts_xyyyy_xyy, ts_xyyyy_xyyy, ts_xyyyy_xyyz, ts_xyyyy_xyz, ts_xyyyy_xyzz, ts_xyyyy_xzz, ts_xyyyy_xzzz, ts_xyyyy_yyy, ts_xyyyy_yyyy, ts_xyyyy_yyyz, ts_xyyyy_yyz, ts_xyyyy_yyzz, ts_xyyyy_yzz, ts_xyyyy_yzzz, ts_xyyyy_zzz, ts_xyyyy_zzzz, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxzz, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyzz, ts_yyyy_xzzz, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyzz, ts_yyyy_yzzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyy_xxxx[i] = 2.0 * ts_yyyy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxy[i] = 2.0 * ts_yyyy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxz[i] = 2.0 * ts_yyyy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxyy[i] = 2.0 * ts_yyyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxyz[i] = 2.0 * ts_yyyy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxzz[i] = 2.0 * ts_yyyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyyy[i] = 2.0 * ts_yyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyyz[i] = 2.0 * ts_yyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyzz[i] = 2.0 * ts_yyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xzzz[i] = 2.0 * ts_yyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyyy[i] = 2.0 * ts_yyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyyz[i] = 2.0 * ts_yyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyzz[i] = 2.0 * ts_yyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yzzz[i] = 2.0 * ts_yyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_zzzz[i] = 2.0 * ts_yyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 165-180 components of targeted buffer : HG

    auto gs_x_xyyyz_xxxx = pbuffer.data(idx_g_hg + 165);

    auto gs_x_xyyyz_xxxy = pbuffer.data(idx_g_hg + 166);

    auto gs_x_xyyyz_xxxz = pbuffer.data(idx_g_hg + 167);

    auto gs_x_xyyyz_xxyy = pbuffer.data(idx_g_hg + 168);

    auto gs_x_xyyyz_xxyz = pbuffer.data(idx_g_hg + 169);

    auto gs_x_xyyyz_xxzz = pbuffer.data(idx_g_hg + 170);

    auto gs_x_xyyyz_xyyy = pbuffer.data(idx_g_hg + 171);

    auto gs_x_xyyyz_xyyz = pbuffer.data(idx_g_hg + 172);

    auto gs_x_xyyyz_xyzz = pbuffer.data(idx_g_hg + 173);

    auto gs_x_xyyyz_xzzz = pbuffer.data(idx_g_hg + 174);

    auto gs_x_xyyyz_yyyy = pbuffer.data(idx_g_hg + 175);

    auto gs_x_xyyyz_yyyz = pbuffer.data(idx_g_hg + 176);

    auto gs_x_xyyyz_yyzz = pbuffer.data(idx_g_hg + 177);

    auto gs_x_xyyyz_yzzz = pbuffer.data(idx_g_hg + 178);

    auto gs_x_xyyyz_zzzz = pbuffer.data(idx_g_hg + 179);

    #pragma omp simd aligned(gc_x, gs_x_xyyyz_xxxx, gs_x_xyyyz_xxxy, gs_x_xyyyz_xxxz, gs_x_xyyyz_xxyy, gs_x_xyyyz_xxyz, gs_x_xyyyz_xxzz, gs_x_xyyyz_xyyy, gs_x_xyyyz_xyyz, gs_x_xyyyz_xyzz, gs_x_xyyyz_xzzz, gs_x_xyyyz_yyyy, gs_x_xyyyz_yyyz, gs_x_xyyyz_yyzz, gs_x_xyyyz_yzzz, gs_x_xyyyz_zzzz, ts_xyyyz_xxx, ts_xyyyz_xxxx, ts_xyyyz_xxxy, ts_xyyyz_xxxz, ts_xyyyz_xxy, ts_xyyyz_xxyy, ts_xyyyz_xxyz, ts_xyyyz_xxz, ts_xyyyz_xxzz, ts_xyyyz_xyy, ts_xyyyz_xyyy, ts_xyyyz_xyyz, ts_xyyyz_xyz, ts_xyyyz_xyzz, ts_xyyyz_xzz, ts_xyyyz_xzzz, ts_xyyyz_yyy, ts_xyyyz_yyyy, ts_xyyyz_yyyz, ts_xyyyz_yyz, ts_xyyyz_yyzz, ts_xyyyz_yzz, ts_xyyyz_yzzz, ts_xyyyz_zzz, ts_xyyyz_zzzz, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxzz, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyzz, ts_yyyz_xzzz, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyzz, ts_yyyz_yzzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyz_xxxx[i] = 2.0 * ts_yyyz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxy[i] = 2.0 * ts_yyyz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxz[i] = 2.0 * ts_yyyz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxyy[i] = 2.0 * ts_yyyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxyz[i] = 2.0 * ts_yyyz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxzz[i] = 2.0 * ts_yyyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyyy[i] = 2.0 * ts_yyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyyz[i] = 2.0 * ts_yyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyzz[i] = 2.0 * ts_yyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xzzz[i] = 2.0 * ts_yyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyyy[i] = 2.0 * ts_yyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyyz[i] = 2.0 * ts_yyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyzz[i] = 2.0 * ts_yyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yzzz[i] = 2.0 * ts_yyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_zzzz[i] = 2.0 * ts_yyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 180-195 components of targeted buffer : HG

    auto gs_x_xyyzz_xxxx = pbuffer.data(idx_g_hg + 180);

    auto gs_x_xyyzz_xxxy = pbuffer.data(idx_g_hg + 181);

    auto gs_x_xyyzz_xxxz = pbuffer.data(idx_g_hg + 182);

    auto gs_x_xyyzz_xxyy = pbuffer.data(idx_g_hg + 183);

    auto gs_x_xyyzz_xxyz = pbuffer.data(idx_g_hg + 184);

    auto gs_x_xyyzz_xxzz = pbuffer.data(idx_g_hg + 185);

    auto gs_x_xyyzz_xyyy = pbuffer.data(idx_g_hg + 186);

    auto gs_x_xyyzz_xyyz = pbuffer.data(idx_g_hg + 187);

    auto gs_x_xyyzz_xyzz = pbuffer.data(idx_g_hg + 188);

    auto gs_x_xyyzz_xzzz = pbuffer.data(idx_g_hg + 189);

    auto gs_x_xyyzz_yyyy = pbuffer.data(idx_g_hg + 190);

    auto gs_x_xyyzz_yyyz = pbuffer.data(idx_g_hg + 191);

    auto gs_x_xyyzz_yyzz = pbuffer.data(idx_g_hg + 192);

    auto gs_x_xyyzz_yzzz = pbuffer.data(idx_g_hg + 193);

    auto gs_x_xyyzz_zzzz = pbuffer.data(idx_g_hg + 194);

    #pragma omp simd aligned(gc_x, gs_x_xyyzz_xxxx, gs_x_xyyzz_xxxy, gs_x_xyyzz_xxxz, gs_x_xyyzz_xxyy, gs_x_xyyzz_xxyz, gs_x_xyyzz_xxzz, gs_x_xyyzz_xyyy, gs_x_xyyzz_xyyz, gs_x_xyyzz_xyzz, gs_x_xyyzz_xzzz, gs_x_xyyzz_yyyy, gs_x_xyyzz_yyyz, gs_x_xyyzz_yyzz, gs_x_xyyzz_yzzz, gs_x_xyyzz_zzzz, ts_xyyzz_xxx, ts_xyyzz_xxxx, ts_xyyzz_xxxy, ts_xyyzz_xxxz, ts_xyyzz_xxy, ts_xyyzz_xxyy, ts_xyyzz_xxyz, ts_xyyzz_xxz, ts_xyyzz_xxzz, ts_xyyzz_xyy, ts_xyyzz_xyyy, ts_xyyzz_xyyz, ts_xyyzz_xyz, ts_xyyzz_xyzz, ts_xyyzz_xzz, ts_xyyzz_xzzz, ts_xyyzz_yyy, ts_xyyzz_yyyy, ts_xyyzz_yyyz, ts_xyyzz_yyz, ts_xyyzz_yyzz, ts_xyyzz_yzz, ts_xyyzz_yzzz, ts_xyyzz_zzz, ts_xyyzz_zzzz, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxzz, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyzz, ts_yyzz_xzzz, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyzz, ts_yyzz_yzzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyzz_xxxx[i] = 2.0 * ts_yyzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxy[i] = 2.0 * ts_yyzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxz[i] = 2.0 * ts_yyzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxyy[i] = 2.0 * ts_yyzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxyz[i] = 2.0 * ts_yyzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxzz[i] = 2.0 * ts_yyzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyyy[i] = 2.0 * ts_yyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyyz[i] = 2.0 * ts_yyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyzz[i] = 2.0 * ts_yyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xzzz[i] = 2.0 * ts_yyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyyy[i] = 2.0 * ts_yyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyyz[i] = 2.0 * ts_yyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyzz[i] = 2.0 * ts_yyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yzzz[i] = 2.0 * ts_yyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_zzzz[i] = 2.0 * ts_yyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 195-210 components of targeted buffer : HG

    auto gs_x_xyzzz_xxxx = pbuffer.data(idx_g_hg + 195);

    auto gs_x_xyzzz_xxxy = pbuffer.data(idx_g_hg + 196);

    auto gs_x_xyzzz_xxxz = pbuffer.data(idx_g_hg + 197);

    auto gs_x_xyzzz_xxyy = pbuffer.data(idx_g_hg + 198);

    auto gs_x_xyzzz_xxyz = pbuffer.data(idx_g_hg + 199);

    auto gs_x_xyzzz_xxzz = pbuffer.data(idx_g_hg + 200);

    auto gs_x_xyzzz_xyyy = pbuffer.data(idx_g_hg + 201);

    auto gs_x_xyzzz_xyyz = pbuffer.data(idx_g_hg + 202);

    auto gs_x_xyzzz_xyzz = pbuffer.data(idx_g_hg + 203);

    auto gs_x_xyzzz_xzzz = pbuffer.data(idx_g_hg + 204);

    auto gs_x_xyzzz_yyyy = pbuffer.data(idx_g_hg + 205);

    auto gs_x_xyzzz_yyyz = pbuffer.data(idx_g_hg + 206);

    auto gs_x_xyzzz_yyzz = pbuffer.data(idx_g_hg + 207);

    auto gs_x_xyzzz_yzzz = pbuffer.data(idx_g_hg + 208);

    auto gs_x_xyzzz_zzzz = pbuffer.data(idx_g_hg + 209);

    #pragma omp simd aligned(gc_x, gs_x_xyzzz_xxxx, gs_x_xyzzz_xxxy, gs_x_xyzzz_xxxz, gs_x_xyzzz_xxyy, gs_x_xyzzz_xxyz, gs_x_xyzzz_xxzz, gs_x_xyzzz_xyyy, gs_x_xyzzz_xyyz, gs_x_xyzzz_xyzz, gs_x_xyzzz_xzzz, gs_x_xyzzz_yyyy, gs_x_xyzzz_yyyz, gs_x_xyzzz_yyzz, gs_x_xyzzz_yzzz, gs_x_xyzzz_zzzz, ts_xyzzz_xxx, ts_xyzzz_xxxx, ts_xyzzz_xxxy, ts_xyzzz_xxxz, ts_xyzzz_xxy, ts_xyzzz_xxyy, ts_xyzzz_xxyz, ts_xyzzz_xxz, ts_xyzzz_xxzz, ts_xyzzz_xyy, ts_xyzzz_xyyy, ts_xyzzz_xyyz, ts_xyzzz_xyz, ts_xyzzz_xyzz, ts_xyzzz_xzz, ts_xyzzz_xzzz, ts_xyzzz_yyy, ts_xyzzz_yyyy, ts_xyzzz_yyyz, ts_xyzzz_yyz, ts_xyzzz_yyzz, ts_xyzzz_yzz, ts_xyzzz_yzzz, ts_xyzzz_zzz, ts_xyzzz_zzzz, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxzz, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyzz, ts_yzzz_xzzz, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyzz, ts_yzzz_yzzz, ts_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzzz_xxxx[i] = 2.0 * ts_yzzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxy[i] = 2.0 * ts_yzzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxz[i] = 2.0 * ts_yzzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxyy[i] = 2.0 * ts_yzzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxyz[i] = 2.0 * ts_yzzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxzz[i] = 2.0 * ts_yzzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyyy[i] = 2.0 * ts_yzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyyz[i] = 2.0 * ts_yzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyzz[i] = 2.0 * ts_yzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xzzz[i] = 2.0 * ts_yzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyyy[i] = 2.0 * ts_yzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyyz[i] = 2.0 * ts_yzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyzz[i] = 2.0 * ts_yzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yzzz[i] = 2.0 * ts_yzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_zzzz[i] = 2.0 * ts_yzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 210-225 components of targeted buffer : HG

    auto gs_x_xzzzz_xxxx = pbuffer.data(idx_g_hg + 210);

    auto gs_x_xzzzz_xxxy = pbuffer.data(idx_g_hg + 211);

    auto gs_x_xzzzz_xxxz = pbuffer.data(idx_g_hg + 212);

    auto gs_x_xzzzz_xxyy = pbuffer.data(idx_g_hg + 213);

    auto gs_x_xzzzz_xxyz = pbuffer.data(idx_g_hg + 214);

    auto gs_x_xzzzz_xxzz = pbuffer.data(idx_g_hg + 215);

    auto gs_x_xzzzz_xyyy = pbuffer.data(idx_g_hg + 216);

    auto gs_x_xzzzz_xyyz = pbuffer.data(idx_g_hg + 217);

    auto gs_x_xzzzz_xyzz = pbuffer.data(idx_g_hg + 218);

    auto gs_x_xzzzz_xzzz = pbuffer.data(idx_g_hg + 219);

    auto gs_x_xzzzz_yyyy = pbuffer.data(idx_g_hg + 220);

    auto gs_x_xzzzz_yyyz = pbuffer.data(idx_g_hg + 221);

    auto gs_x_xzzzz_yyzz = pbuffer.data(idx_g_hg + 222);

    auto gs_x_xzzzz_yzzz = pbuffer.data(idx_g_hg + 223);

    auto gs_x_xzzzz_zzzz = pbuffer.data(idx_g_hg + 224);

    #pragma omp simd aligned(gc_x, gs_x_xzzzz_xxxx, gs_x_xzzzz_xxxy, gs_x_xzzzz_xxxz, gs_x_xzzzz_xxyy, gs_x_xzzzz_xxyz, gs_x_xzzzz_xxzz, gs_x_xzzzz_xyyy, gs_x_xzzzz_xyyz, gs_x_xzzzz_xyzz, gs_x_xzzzz_xzzz, gs_x_xzzzz_yyyy, gs_x_xzzzz_yyyz, gs_x_xzzzz_yyzz, gs_x_xzzzz_yzzz, gs_x_xzzzz_zzzz, ts_xzzzz_xxx, ts_xzzzz_xxxx, ts_xzzzz_xxxy, ts_xzzzz_xxxz, ts_xzzzz_xxy, ts_xzzzz_xxyy, ts_xzzzz_xxyz, ts_xzzzz_xxz, ts_xzzzz_xxzz, ts_xzzzz_xyy, ts_xzzzz_xyyy, ts_xzzzz_xyyz, ts_xzzzz_xyz, ts_xzzzz_xyzz, ts_xzzzz_xzz, ts_xzzzz_xzzz, ts_xzzzz_yyy, ts_xzzzz_yyyy, ts_xzzzz_yyyz, ts_xzzzz_yyz, ts_xzzzz_yyzz, ts_xzzzz_yzz, ts_xzzzz_yzzz, ts_xzzzz_zzz, ts_xzzzz_zzzz, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzzz_xxxx[i] = 2.0 * ts_zzzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxy[i] = 2.0 * ts_zzzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxz[i] = 2.0 * ts_zzzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxyy[i] = 2.0 * ts_zzzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxyz[i] = 2.0 * ts_zzzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxzz[i] = 2.0 * ts_zzzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyyy[i] = 2.0 * ts_zzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyyz[i] = 2.0 * ts_zzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyzz[i] = 2.0 * ts_zzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xzzz[i] = 2.0 * ts_zzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyyy[i] = 2.0 * ts_zzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyyz[i] = 2.0 * ts_zzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyzz[i] = 2.0 * ts_zzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yzzz[i] = 2.0 * ts_zzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_zzzz[i] = 2.0 * ts_zzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 225-240 components of targeted buffer : HG

    auto gs_x_yyyyy_xxxx = pbuffer.data(idx_g_hg + 225);

    auto gs_x_yyyyy_xxxy = pbuffer.data(idx_g_hg + 226);

    auto gs_x_yyyyy_xxxz = pbuffer.data(idx_g_hg + 227);

    auto gs_x_yyyyy_xxyy = pbuffer.data(idx_g_hg + 228);

    auto gs_x_yyyyy_xxyz = pbuffer.data(idx_g_hg + 229);

    auto gs_x_yyyyy_xxzz = pbuffer.data(idx_g_hg + 230);

    auto gs_x_yyyyy_xyyy = pbuffer.data(idx_g_hg + 231);

    auto gs_x_yyyyy_xyyz = pbuffer.data(idx_g_hg + 232);

    auto gs_x_yyyyy_xyzz = pbuffer.data(idx_g_hg + 233);

    auto gs_x_yyyyy_xzzz = pbuffer.data(idx_g_hg + 234);

    auto gs_x_yyyyy_yyyy = pbuffer.data(idx_g_hg + 235);

    auto gs_x_yyyyy_yyyz = pbuffer.data(idx_g_hg + 236);

    auto gs_x_yyyyy_yyzz = pbuffer.data(idx_g_hg + 237);

    auto gs_x_yyyyy_yzzz = pbuffer.data(idx_g_hg + 238);

    auto gs_x_yyyyy_zzzz = pbuffer.data(idx_g_hg + 239);

    #pragma omp simd aligned(gc_x, gs_x_yyyyy_xxxx, gs_x_yyyyy_xxxy, gs_x_yyyyy_xxxz, gs_x_yyyyy_xxyy, gs_x_yyyyy_xxyz, gs_x_yyyyy_xxzz, gs_x_yyyyy_xyyy, gs_x_yyyyy_xyyz, gs_x_yyyyy_xyzz, gs_x_yyyyy_xzzz, gs_x_yyyyy_yyyy, gs_x_yyyyy_yyyz, gs_x_yyyyy_yyzz, gs_x_yyyyy_yzzz, gs_x_yyyyy_zzzz, ts_yyyyy_xxx, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxxz, ts_yyyyy_xxy, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xxz, ts_yyyyy_xxzz, ts_yyyyy_xyy, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyz, ts_yyyyy_xyzz, ts_yyyyy_xzz, ts_yyyyy_xzzz, ts_yyyyy_yyy, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyz, ts_yyyyy_yyzz, ts_yyyyy_yzz, ts_yyyyy_yzzz, ts_yyyyy_zzz, ts_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyy_xxxx[i] = 8.0 * ts_yyyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxy[i] = 6.0 * ts_yyyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxz[i] = 6.0 * ts_yyyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxyy[i] = 4.0 * ts_yyyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxyz[i] = 4.0 * ts_yyyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxzz[i] = 4.0 * ts_yyyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyyy[i] = 2.0 * ts_yyyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyyz[i] = 2.0 * ts_yyyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyzz[i] = 2.0 * ts_yyyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xzzz[i] = 2.0 * ts_yyyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyyy[i] = 2.0 * ts_yyyyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyyz[i] = 2.0 * ts_yyyyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyzz[i] = 2.0 * ts_yyyyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yzzz[i] = 2.0 * ts_yyyyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_zzzz[i] = 2.0 * ts_yyyyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 240-255 components of targeted buffer : HG

    auto gs_x_yyyyz_xxxx = pbuffer.data(idx_g_hg + 240);

    auto gs_x_yyyyz_xxxy = pbuffer.data(idx_g_hg + 241);

    auto gs_x_yyyyz_xxxz = pbuffer.data(idx_g_hg + 242);

    auto gs_x_yyyyz_xxyy = pbuffer.data(idx_g_hg + 243);

    auto gs_x_yyyyz_xxyz = pbuffer.data(idx_g_hg + 244);

    auto gs_x_yyyyz_xxzz = pbuffer.data(idx_g_hg + 245);

    auto gs_x_yyyyz_xyyy = pbuffer.data(idx_g_hg + 246);

    auto gs_x_yyyyz_xyyz = pbuffer.data(idx_g_hg + 247);

    auto gs_x_yyyyz_xyzz = pbuffer.data(idx_g_hg + 248);

    auto gs_x_yyyyz_xzzz = pbuffer.data(idx_g_hg + 249);

    auto gs_x_yyyyz_yyyy = pbuffer.data(idx_g_hg + 250);

    auto gs_x_yyyyz_yyyz = pbuffer.data(idx_g_hg + 251);

    auto gs_x_yyyyz_yyzz = pbuffer.data(idx_g_hg + 252);

    auto gs_x_yyyyz_yzzz = pbuffer.data(idx_g_hg + 253);

    auto gs_x_yyyyz_zzzz = pbuffer.data(idx_g_hg + 254);

    #pragma omp simd aligned(gc_x, gs_x_yyyyz_xxxx, gs_x_yyyyz_xxxy, gs_x_yyyyz_xxxz, gs_x_yyyyz_xxyy, gs_x_yyyyz_xxyz, gs_x_yyyyz_xxzz, gs_x_yyyyz_xyyy, gs_x_yyyyz_xyyz, gs_x_yyyyz_xyzz, gs_x_yyyyz_xzzz, gs_x_yyyyz_yyyy, gs_x_yyyyz_yyyz, gs_x_yyyyz_yyzz, gs_x_yyyyz_yzzz, gs_x_yyyyz_zzzz, ts_yyyyz_xxx, ts_yyyyz_xxxx, ts_yyyyz_xxxy, ts_yyyyz_xxxz, ts_yyyyz_xxy, ts_yyyyz_xxyy, ts_yyyyz_xxyz, ts_yyyyz_xxz, ts_yyyyz_xxzz, ts_yyyyz_xyy, ts_yyyyz_xyyy, ts_yyyyz_xyyz, ts_yyyyz_xyz, ts_yyyyz_xyzz, ts_yyyyz_xzz, ts_yyyyz_xzzz, ts_yyyyz_yyy, ts_yyyyz_yyyy, ts_yyyyz_yyyz, ts_yyyyz_yyz, ts_yyyyz_yyzz, ts_yyyyz_yzz, ts_yyyyz_yzzz, ts_yyyyz_zzz, ts_yyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyz_xxxx[i] = 8.0 * ts_yyyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxy[i] = 6.0 * ts_yyyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxz[i] = 6.0 * ts_yyyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxyy[i] = 4.0 * ts_yyyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxyz[i] = 4.0 * ts_yyyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxzz[i] = 4.0 * ts_yyyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyyy[i] = 2.0 * ts_yyyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyyz[i] = 2.0 * ts_yyyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyzz[i] = 2.0 * ts_yyyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xzzz[i] = 2.0 * ts_yyyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyyy[i] = 2.0 * ts_yyyyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyyz[i] = 2.0 * ts_yyyyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyzz[i] = 2.0 * ts_yyyyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yzzz[i] = 2.0 * ts_yyyyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_zzzz[i] = 2.0 * ts_yyyyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 255-270 components of targeted buffer : HG

    auto gs_x_yyyzz_xxxx = pbuffer.data(idx_g_hg + 255);

    auto gs_x_yyyzz_xxxy = pbuffer.data(idx_g_hg + 256);

    auto gs_x_yyyzz_xxxz = pbuffer.data(idx_g_hg + 257);

    auto gs_x_yyyzz_xxyy = pbuffer.data(idx_g_hg + 258);

    auto gs_x_yyyzz_xxyz = pbuffer.data(idx_g_hg + 259);

    auto gs_x_yyyzz_xxzz = pbuffer.data(idx_g_hg + 260);

    auto gs_x_yyyzz_xyyy = pbuffer.data(idx_g_hg + 261);

    auto gs_x_yyyzz_xyyz = pbuffer.data(idx_g_hg + 262);

    auto gs_x_yyyzz_xyzz = pbuffer.data(idx_g_hg + 263);

    auto gs_x_yyyzz_xzzz = pbuffer.data(idx_g_hg + 264);

    auto gs_x_yyyzz_yyyy = pbuffer.data(idx_g_hg + 265);

    auto gs_x_yyyzz_yyyz = pbuffer.data(idx_g_hg + 266);

    auto gs_x_yyyzz_yyzz = pbuffer.data(idx_g_hg + 267);

    auto gs_x_yyyzz_yzzz = pbuffer.data(idx_g_hg + 268);

    auto gs_x_yyyzz_zzzz = pbuffer.data(idx_g_hg + 269);

    #pragma omp simd aligned(gc_x, gs_x_yyyzz_xxxx, gs_x_yyyzz_xxxy, gs_x_yyyzz_xxxz, gs_x_yyyzz_xxyy, gs_x_yyyzz_xxyz, gs_x_yyyzz_xxzz, gs_x_yyyzz_xyyy, gs_x_yyyzz_xyyz, gs_x_yyyzz_xyzz, gs_x_yyyzz_xzzz, gs_x_yyyzz_yyyy, gs_x_yyyzz_yyyz, gs_x_yyyzz_yyzz, gs_x_yyyzz_yzzz, gs_x_yyyzz_zzzz, ts_yyyzz_xxx, ts_yyyzz_xxxx, ts_yyyzz_xxxy, ts_yyyzz_xxxz, ts_yyyzz_xxy, ts_yyyzz_xxyy, ts_yyyzz_xxyz, ts_yyyzz_xxz, ts_yyyzz_xxzz, ts_yyyzz_xyy, ts_yyyzz_xyyy, ts_yyyzz_xyyz, ts_yyyzz_xyz, ts_yyyzz_xyzz, ts_yyyzz_xzz, ts_yyyzz_xzzz, ts_yyyzz_yyy, ts_yyyzz_yyyy, ts_yyyzz_yyyz, ts_yyyzz_yyz, ts_yyyzz_yyzz, ts_yyyzz_yzz, ts_yyyzz_yzzz, ts_yyyzz_zzz, ts_yyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyzz_xxxx[i] = 8.0 * ts_yyyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxy[i] = 6.0 * ts_yyyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxz[i] = 6.0 * ts_yyyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxyy[i] = 4.0 * ts_yyyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxyz[i] = 4.0 * ts_yyyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxzz[i] = 4.0 * ts_yyyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyyy[i] = 2.0 * ts_yyyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyyz[i] = 2.0 * ts_yyyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyzz[i] = 2.0 * ts_yyyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xzzz[i] = 2.0 * ts_yyyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyyy[i] = 2.0 * ts_yyyzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyyz[i] = 2.0 * ts_yyyzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyzz[i] = 2.0 * ts_yyyzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yzzz[i] = 2.0 * ts_yyyzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_zzzz[i] = 2.0 * ts_yyyzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 270-285 components of targeted buffer : HG

    auto gs_x_yyzzz_xxxx = pbuffer.data(idx_g_hg + 270);

    auto gs_x_yyzzz_xxxy = pbuffer.data(idx_g_hg + 271);

    auto gs_x_yyzzz_xxxz = pbuffer.data(idx_g_hg + 272);

    auto gs_x_yyzzz_xxyy = pbuffer.data(idx_g_hg + 273);

    auto gs_x_yyzzz_xxyz = pbuffer.data(idx_g_hg + 274);

    auto gs_x_yyzzz_xxzz = pbuffer.data(idx_g_hg + 275);

    auto gs_x_yyzzz_xyyy = pbuffer.data(idx_g_hg + 276);

    auto gs_x_yyzzz_xyyz = pbuffer.data(idx_g_hg + 277);

    auto gs_x_yyzzz_xyzz = pbuffer.data(idx_g_hg + 278);

    auto gs_x_yyzzz_xzzz = pbuffer.data(idx_g_hg + 279);

    auto gs_x_yyzzz_yyyy = pbuffer.data(idx_g_hg + 280);

    auto gs_x_yyzzz_yyyz = pbuffer.data(idx_g_hg + 281);

    auto gs_x_yyzzz_yyzz = pbuffer.data(idx_g_hg + 282);

    auto gs_x_yyzzz_yzzz = pbuffer.data(idx_g_hg + 283);

    auto gs_x_yyzzz_zzzz = pbuffer.data(idx_g_hg + 284);

    #pragma omp simd aligned(gc_x, gs_x_yyzzz_xxxx, gs_x_yyzzz_xxxy, gs_x_yyzzz_xxxz, gs_x_yyzzz_xxyy, gs_x_yyzzz_xxyz, gs_x_yyzzz_xxzz, gs_x_yyzzz_xyyy, gs_x_yyzzz_xyyz, gs_x_yyzzz_xyzz, gs_x_yyzzz_xzzz, gs_x_yyzzz_yyyy, gs_x_yyzzz_yyyz, gs_x_yyzzz_yyzz, gs_x_yyzzz_yzzz, gs_x_yyzzz_zzzz, ts_yyzzz_xxx, ts_yyzzz_xxxx, ts_yyzzz_xxxy, ts_yyzzz_xxxz, ts_yyzzz_xxy, ts_yyzzz_xxyy, ts_yyzzz_xxyz, ts_yyzzz_xxz, ts_yyzzz_xxzz, ts_yyzzz_xyy, ts_yyzzz_xyyy, ts_yyzzz_xyyz, ts_yyzzz_xyz, ts_yyzzz_xyzz, ts_yyzzz_xzz, ts_yyzzz_xzzz, ts_yyzzz_yyy, ts_yyzzz_yyyy, ts_yyzzz_yyyz, ts_yyzzz_yyz, ts_yyzzz_yyzz, ts_yyzzz_yzz, ts_yyzzz_yzzz, ts_yyzzz_zzz, ts_yyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzzz_xxxx[i] = 8.0 * ts_yyzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxy[i] = 6.0 * ts_yyzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxz[i] = 6.0 * ts_yyzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxyy[i] = 4.0 * ts_yyzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxyz[i] = 4.0 * ts_yyzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxzz[i] = 4.0 * ts_yyzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyyy[i] = 2.0 * ts_yyzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyyz[i] = 2.0 * ts_yyzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyzz[i] = 2.0 * ts_yyzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xzzz[i] = 2.0 * ts_yyzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyyy[i] = 2.0 * ts_yyzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyyz[i] = 2.0 * ts_yyzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyzz[i] = 2.0 * ts_yyzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yzzz[i] = 2.0 * ts_yyzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_zzzz[i] = 2.0 * ts_yyzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 285-300 components of targeted buffer : HG

    auto gs_x_yzzzz_xxxx = pbuffer.data(idx_g_hg + 285);

    auto gs_x_yzzzz_xxxy = pbuffer.data(idx_g_hg + 286);

    auto gs_x_yzzzz_xxxz = pbuffer.data(idx_g_hg + 287);

    auto gs_x_yzzzz_xxyy = pbuffer.data(idx_g_hg + 288);

    auto gs_x_yzzzz_xxyz = pbuffer.data(idx_g_hg + 289);

    auto gs_x_yzzzz_xxzz = pbuffer.data(idx_g_hg + 290);

    auto gs_x_yzzzz_xyyy = pbuffer.data(idx_g_hg + 291);

    auto gs_x_yzzzz_xyyz = pbuffer.data(idx_g_hg + 292);

    auto gs_x_yzzzz_xyzz = pbuffer.data(idx_g_hg + 293);

    auto gs_x_yzzzz_xzzz = pbuffer.data(idx_g_hg + 294);

    auto gs_x_yzzzz_yyyy = pbuffer.data(idx_g_hg + 295);

    auto gs_x_yzzzz_yyyz = pbuffer.data(idx_g_hg + 296);

    auto gs_x_yzzzz_yyzz = pbuffer.data(idx_g_hg + 297);

    auto gs_x_yzzzz_yzzz = pbuffer.data(idx_g_hg + 298);

    auto gs_x_yzzzz_zzzz = pbuffer.data(idx_g_hg + 299);

    #pragma omp simd aligned(gc_x, gs_x_yzzzz_xxxx, gs_x_yzzzz_xxxy, gs_x_yzzzz_xxxz, gs_x_yzzzz_xxyy, gs_x_yzzzz_xxyz, gs_x_yzzzz_xxzz, gs_x_yzzzz_xyyy, gs_x_yzzzz_xyyz, gs_x_yzzzz_xyzz, gs_x_yzzzz_xzzz, gs_x_yzzzz_yyyy, gs_x_yzzzz_yyyz, gs_x_yzzzz_yyzz, gs_x_yzzzz_yzzz, gs_x_yzzzz_zzzz, ts_yzzzz_xxx, ts_yzzzz_xxxx, ts_yzzzz_xxxy, ts_yzzzz_xxxz, ts_yzzzz_xxy, ts_yzzzz_xxyy, ts_yzzzz_xxyz, ts_yzzzz_xxz, ts_yzzzz_xxzz, ts_yzzzz_xyy, ts_yzzzz_xyyy, ts_yzzzz_xyyz, ts_yzzzz_xyz, ts_yzzzz_xyzz, ts_yzzzz_xzz, ts_yzzzz_xzzz, ts_yzzzz_yyy, ts_yzzzz_yyyy, ts_yzzzz_yyyz, ts_yzzzz_yyz, ts_yzzzz_yyzz, ts_yzzzz_yzz, ts_yzzzz_yzzz, ts_yzzzz_zzz, ts_yzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzzz_xxxx[i] = 8.0 * ts_yzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxy[i] = 6.0 * ts_yzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxz[i] = 6.0 * ts_yzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxyy[i] = 4.0 * ts_yzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxyz[i] = 4.0 * ts_yzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxzz[i] = 4.0 * ts_yzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyyy[i] = 2.0 * ts_yzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyyz[i] = 2.0 * ts_yzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyzz[i] = 2.0 * ts_yzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xzzz[i] = 2.0 * ts_yzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyyy[i] = 2.0 * ts_yzzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyyz[i] = 2.0 * ts_yzzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyzz[i] = 2.0 * ts_yzzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yzzz[i] = 2.0 * ts_yzzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_zzzz[i] = 2.0 * ts_yzzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 300-315 components of targeted buffer : HG

    auto gs_x_zzzzz_xxxx = pbuffer.data(idx_g_hg + 300);

    auto gs_x_zzzzz_xxxy = pbuffer.data(idx_g_hg + 301);

    auto gs_x_zzzzz_xxxz = pbuffer.data(idx_g_hg + 302);

    auto gs_x_zzzzz_xxyy = pbuffer.data(idx_g_hg + 303);

    auto gs_x_zzzzz_xxyz = pbuffer.data(idx_g_hg + 304);

    auto gs_x_zzzzz_xxzz = pbuffer.data(idx_g_hg + 305);

    auto gs_x_zzzzz_xyyy = pbuffer.data(idx_g_hg + 306);

    auto gs_x_zzzzz_xyyz = pbuffer.data(idx_g_hg + 307);

    auto gs_x_zzzzz_xyzz = pbuffer.data(idx_g_hg + 308);

    auto gs_x_zzzzz_xzzz = pbuffer.data(idx_g_hg + 309);

    auto gs_x_zzzzz_yyyy = pbuffer.data(idx_g_hg + 310);

    auto gs_x_zzzzz_yyyz = pbuffer.data(idx_g_hg + 311);

    auto gs_x_zzzzz_yyzz = pbuffer.data(idx_g_hg + 312);

    auto gs_x_zzzzz_yzzz = pbuffer.data(idx_g_hg + 313);

    auto gs_x_zzzzz_zzzz = pbuffer.data(idx_g_hg + 314);

    #pragma omp simd aligned(gc_x, gs_x_zzzzz_xxxx, gs_x_zzzzz_xxxy, gs_x_zzzzz_xxxz, gs_x_zzzzz_xxyy, gs_x_zzzzz_xxyz, gs_x_zzzzz_xxzz, gs_x_zzzzz_xyyy, gs_x_zzzzz_xyyz, gs_x_zzzzz_xyzz, gs_x_zzzzz_xzzz, gs_x_zzzzz_yyyy, gs_x_zzzzz_yyyz, gs_x_zzzzz_yyzz, gs_x_zzzzz_yzzz, gs_x_zzzzz_zzzz, ts_zzzzz_xxx, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxy, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxz, ts_zzzzz_xxzz, ts_zzzzz_xyy, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyz, ts_zzzzz_xyzz, ts_zzzzz_xzz, ts_zzzzz_xzzz, ts_zzzzz_yyy, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyz, ts_zzzzz_yyzz, ts_zzzzz_yzz, ts_zzzzz_yzzz, ts_zzzzz_zzz, ts_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzzz_xxxx[i] = 8.0 * ts_zzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxy[i] = 6.0 * ts_zzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxz[i] = 6.0 * ts_zzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxyy[i] = 4.0 * ts_zzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxyz[i] = 4.0 * ts_zzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxzz[i] = 4.0 * ts_zzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyyy[i] = 2.0 * ts_zzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyyz[i] = 2.0 * ts_zzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyzz[i] = 2.0 * ts_zzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xzzz[i] = 2.0 * ts_zzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyyy[i] = 2.0 * ts_zzzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyyz[i] = 2.0 * ts_zzzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyzz[i] = 2.0 * ts_zzzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yzzz[i] = 2.0 * ts_zzzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_zzzz[i] = 2.0 * ts_zzzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 315-330 components of targeted buffer : HG

    auto gs_y_xxxxx_xxxx = pbuffer.data(idx_g_hg + 315);

    auto gs_y_xxxxx_xxxy = pbuffer.data(idx_g_hg + 316);

    auto gs_y_xxxxx_xxxz = pbuffer.data(idx_g_hg + 317);

    auto gs_y_xxxxx_xxyy = pbuffer.data(idx_g_hg + 318);

    auto gs_y_xxxxx_xxyz = pbuffer.data(idx_g_hg + 319);

    auto gs_y_xxxxx_xxzz = pbuffer.data(idx_g_hg + 320);

    auto gs_y_xxxxx_xyyy = pbuffer.data(idx_g_hg + 321);

    auto gs_y_xxxxx_xyyz = pbuffer.data(idx_g_hg + 322);

    auto gs_y_xxxxx_xyzz = pbuffer.data(idx_g_hg + 323);

    auto gs_y_xxxxx_xzzz = pbuffer.data(idx_g_hg + 324);

    auto gs_y_xxxxx_yyyy = pbuffer.data(idx_g_hg + 325);

    auto gs_y_xxxxx_yyyz = pbuffer.data(idx_g_hg + 326);

    auto gs_y_xxxxx_yyzz = pbuffer.data(idx_g_hg + 327);

    auto gs_y_xxxxx_yzzz = pbuffer.data(idx_g_hg + 328);

    auto gs_y_xxxxx_zzzz = pbuffer.data(idx_g_hg + 329);

    #pragma omp simd aligned(gc_y, gs_y_xxxxx_xxxx, gs_y_xxxxx_xxxy, gs_y_xxxxx_xxxz, gs_y_xxxxx_xxyy, gs_y_xxxxx_xxyz, gs_y_xxxxx_xxzz, gs_y_xxxxx_xyyy, gs_y_xxxxx_xyyz, gs_y_xxxxx_xyzz, gs_y_xxxxx_xzzz, gs_y_xxxxx_yyyy, gs_y_xxxxx_yyyz, gs_y_xxxxx_yyzz, gs_y_xxxxx_yzzz, gs_y_xxxxx_zzzz, ts_xxxxx_xxx, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxy, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxz, ts_xxxxx_xxzz, ts_xxxxx_xyy, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyz, ts_xxxxx_xyzz, ts_xxxxx_xzz, ts_xxxxx_xzzz, ts_xxxxx_yyy, ts_xxxxx_yyyy, ts_xxxxx_yyyz, ts_xxxxx_yyz, ts_xxxxx_yyzz, ts_xxxxx_yzz, ts_xxxxx_yzzz, ts_xxxxx_zzz, ts_xxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxx_xxxx[i] = 2.0 * ts_xxxxx_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxy[i] = 2.0 * ts_xxxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxz[i] = 2.0 * ts_xxxxx_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxyy[i] = 4.0 * ts_xxxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxyz[i] = 2.0 * ts_xxxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxzz[i] = 2.0 * ts_xxxxx_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyyy[i] = 6.0 * ts_xxxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyyz[i] = 4.0 * ts_xxxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyzz[i] = 2.0 * ts_xxxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xzzz[i] = 2.0 * ts_xxxxx_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyyy[i] = 8.0 * ts_xxxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyyz[i] = 6.0 * ts_xxxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyzz[i] = 4.0 * ts_xxxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yzzz[i] = 2.0 * ts_xxxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_zzzz[i] = 2.0 * ts_xxxxx_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 330-345 components of targeted buffer : HG

    auto gs_y_xxxxy_xxxx = pbuffer.data(idx_g_hg + 330);

    auto gs_y_xxxxy_xxxy = pbuffer.data(idx_g_hg + 331);

    auto gs_y_xxxxy_xxxz = pbuffer.data(idx_g_hg + 332);

    auto gs_y_xxxxy_xxyy = pbuffer.data(idx_g_hg + 333);

    auto gs_y_xxxxy_xxyz = pbuffer.data(idx_g_hg + 334);

    auto gs_y_xxxxy_xxzz = pbuffer.data(idx_g_hg + 335);

    auto gs_y_xxxxy_xyyy = pbuffer.data(idx_g_hg + 336);

    auto gs_y_xxxxy_xyyz = pbuffer.data(idx_g_hg + 337);

    auto gs_y_xxxxy_xyzz = pbuffer.data(idx_g_hg + 338);

    auto gs_y_xxxxy_xzzz = pbuffer.data(idx_g_hg + 339);

    auto gs_y_xxxxy_yyyy = pbuffer.data(idx_g_hg + 340);

    auto gs_y_xxxxy_yyyz = pbuffer.data(idx_g_hg + 341);

    auto gs_y_xxxxy_yyzz = pbuffer.data(idx_g_hg + 342);

    auto gs_y_xxxxy_yzzz = pbuffer.data(idx_g_hg + 343);

    auto gs_y_xxxxy_zzzz = pbuffer.data(idx_g_hg + 344);

    #pragma omp simd aligned(gc_y, gs_y_xxxxy_xxxx, gs_y_xxxxy_xxxy, gs_y_xxxxy_xxxz, gs_y_xxxxy_xxyy, gs_y_xxxxy_xxyz, gs_y_xxxxy_xxzz, gs_y_xxxxy_xyyy, gs_y_xxxxy_xyyz, gs_y_xxxxy_xyzz, gs_y_xxxxy_xzzz, gs_y_xxxxy_yyyy, gs_y_xxxxy_yyyz, gs_y_xxxxy_yyzz, gs_y_xxxxy_yzzz, gs_y_xxxxy_zzzz, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxzz, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyzz, ts_xxxx_xzzz, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyzz, ts_xxxx_yzzz, ts_xxxx_zzzz, ts_xxxxy_xxx, ts_xxxxy_xxxx, ts_xxxxy_xxxy, ts_xxxxy_xxxz, ts_xxxxy_xxy, ts_xxxxy_xxyy, ts_xxxxy_xxyz, ts_xxxxy_xxz, ts_xxxxy_xxzz, ts_xxxxy_xyy, ts_xxxxy_xyyy, ts_xxxxy_xyyz, ts_xxxxy_xyz, ts_xxxxy_xyzz, ts_xxxxy_xzz, ts_xxxxy_xzzz, ts_xxxxy_yyy, ts_xxxxy_yyyy, ts_xxxxy_yyyz, ts_xxxxy_yyz, ts_xxxxy_yyzz, ts_xxxxy_yzz, ts_xxxxy_yzzz, ts_xxxxy_zzz, ts_xxxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxy_xxxx[i] = 2.0 * ts_xxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxy[i] = 2.0 * ts_xxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxz[i] = 2.0 * ts_xxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxyy[i] = 2.0 * ts_xxxx_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxyz[i] = 2.0 * ts_xxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxzz[i] = 2.0 * ts_xxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyyy[i] = 2.0 * ts_xxxx_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyyz[i] = 2.0 * ts_xxxx_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyzz[i] = 2.0 * ts_xxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xzzz[i] = 2.0 * ts_xxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyyy[i] = 2.0 * ts_xxxx_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyyz[i] = 2.0 * ts_xxxx_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyzz[i] = 2.0 * ts_xxxx_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yzzz[i] = 2.0 * ts_xxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_zzzz[i] = 2.0 * ts_xxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 345-360 components of targeted buffer : HG

    auto gs_y_xxxxz_xxxx = pbuffer.data(idx_g_hg + 345);

    auto gs_y_xxxxz_xxxy = pbuffer.data(idx_g_hg + 346);

    auto gs_y_xxxxz_xxxz = pbuffer.data(idx_g_hg + 347);

    auto gs_y_xxxxz_xxyy = pbuffer.data(idx_g_hg + 348);

    auto gs_y_xxxxz_xxyz = pbuffer.data(idx_g_hg + 349);

    auto gs_y_xxxxz_xxzz = pbuffer.data(idx_g_hg + 350);

    auto gs_y_xxxxz_xyyy = pbuffer.data(idx_g_hg + 351);

    auto gs_y_xxxxz_xyyz = pbuffer.data(idx_g_hg + 352);

    auto gs_y_xxxxz_xyzz = pbuffer.data(idx_g_hg + 353);

    auto gs_y_xxxxz_xzzz = pbuffer.data(idx_g_hg + 354);

    auto gs_y_xxxxz_yyyy = pbuffer.data(idx_g_hg + 355);

    auto gs_y_xxxxz_yyyz = pbuffer.data(idx_g_hg + 356);

    auto gs_y_xxxxz_yyzz = pbuffer.data(idx_g_hg + 357);

    auto gs_y_xxxxz_yzzz = pbuffer.data(idx_g_hg + 358);

    auto gs_y_xxxxz_zzzz = pbuffer.data(idx_g_hg + 359);

    #pragma omp simd aligned(gc_y, gs_y_xxxxz_xxxx, gs_y_xxxxz_xxxy, gs_y_xxxxz_xxxz, gs_y_xxxxz_xxyy, gs_y_xxxxz_xxyz, gs_y_xxxxz_xxzz, gs_y_xxxxz_xyyy, gs_y_xxxxz_xyyz, gs_y_xxxxz_xyzz, gs_y_xxxxz_xzzz, gs_y_xxxxz_yyyy, gs_y_xxxxz_yyyz, gs_y_xxxxz_yyzz, gs_y_xxxxz_yzzz, gs_y_xxxxz_zzzz, ts_xxxxz_xxx, ts_xxxxz_xxxx, ts_xxxxz_xxxy, ts_xxxxz_xxxz, ts_xxxxz_xxy, ts_xxxxz_xxyy, ts_xxxxz_xxyz, ts_xxxxz_xxz, ts_xxxxz_xxzz, ts_xxxxz_xyy, ts_xxxxz_xyyy, ts_xxxxz_xyyz, ts_xxxxz_xyz, ts_xxxxz_xyzz, ts_xxxxz_xzz, ts_xxxxz_xzzz, ts_xxxxz_yyy, ts_xxxxz_yyyy, ts_xxxxz_yyyz, ts_xxxxz_yyz, ts_xxxxz_yyzz, ts_xxxxz_yzz, ts_xxxxz_yzzz, ts_xxxxz_zzz, ts_xxxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxz_xxxx[i] = 2.0 * ts_xxxxz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxy[i] = 2.0 * ts_xxxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxz[i] = 2.0 * ts_xxxxz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxyy[i] = 4.0 * ts_xxxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxyz[i] = 2.0 * ts_xxxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxzz[i] = 2.0 * ts_xxxxz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyyy[i] = 6.0 * ts_xxxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyyz[i] = 4.0 * ts_xxxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyzz[i] = 2.0 * ts_xxxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xzzz[i] = 2.0 * ts_xxxxz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyyy[i] = 8.0 * ts_xxxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyyz[i] = 6.0 * ts_xxxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyzz[i] = 4.0 * ts_xxxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yzzz[i] = 2.0 * ts_xxxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_zzzz[i] = 2.0 * ts_xxxxz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 360-375 components of targeted buffer : HG

    auto gs_y_xxxyy_xxxx = pbuffer.data(idx_g_hg + 360);

    auto gs_y_xxxyy_xxxy = pbuffer.data(idx_g_hg + 361);

    auto gs_y_xxxyy_xxxz = pbuffer.data(idx_g_hg + 362);

    auto gs_y_xxxyy_xxyy = pbuffer.data(idx_g_hg + 363);

    auto gs_y_xxxyy_xxyz = pbuffer.data(idx_g_hg + 364);

    auto gs_y_xxxyy_xxzz = pbuffer.data(idx_g_hg + 365);

    auto gs_y_xxxyy_xyyy = pbuffer.data(idx_g_hg + 366);

    auto gs_y_xxxyy_xyyz = pbuffer.data(idx_g_hg + 367);

    auto gs_y_xxxyy_xyzz = pbuffer.data(idx_g_hg + 368);

    auto gs_y_xxxyy_xzzz = pbuffer.data(idx_g_hg + 369);

    auto gs_y_xxxyy_yyyy = pbuffer.data(idx_g_hg + 370);

    auto gs_y_xxxyy_yyyz = pbuffer.data(idx_g_hg + 371);

    auto gs_y_xxxyy_yyzz = pbuffer.data(idx_g_hg + 372);

    auto gs_y_xxxyy_yzzz = pbuffer.data(idx_g_hg + 373);

    auto gs_y_xxxyy_zzzz = pbuffer.data(idx_g_hg + 374);

    #pragma omp simd aligned(gc_y, gs_y_xxxyy_xxxx, gs_y_xxxyy_xxxy, gs_y_xxxyy_xxxz, gs_y_xxxyy_xxyy, gs_y_xxxyy_xxyz, gs_y_xxxyy_xxzz, gs_y_xxxyy_xyyy, gs_y_xxxyy_xyyz, gs_y_xxxyy_xyzz, gs_y_xxxyy_xzzz, gs_y_xxxyy_yyyy, gs_y_xxxyy_yyyz, gs_y_xxxyy_yyzz, gs_y_xxxyy_yzzz, gs_y_xxxyy_zzzz, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxzz, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyzz, ts_xxxy_xzzz, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyzz, ts_xxxy_yzzz, ts_xxxy_zzzz, ts_xxxyy_xxx, ts_xxxyy_xxxx, ts_xxxyy_xxxy, ts_xxxyy_xxxz, ts_xxxyy_xxy, ts_xxxyy_xxyy, ts_xxxyy_xxyz, ts_xxxyy_xxz, ts_xxxyy_xxzz, ts_xxxyy_xyy, ts_xxxyy_xyyy, ts_xxxyy_xyyz, ts_xxxyy_xyz, ts_xxxyy_xyzz, ts_xxxyy_xzz, ts_xxxyy_xzzz, ts_xxxyy_yyy, ts_xxxyy_yyyy, ts_xxxyy_yyyz, ts_xxxyy_yyz, ts_xxxyy_yyzz, ts_xxxyy_yzz, ts_xxxyy_yzzz, ts_xxxyy_zzz, ts_xxxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyy_xxxx[i] = 4.0 * ts_xxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxy[i] = 4.0 * ts_xxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxz[i] = 4.0 * ts_xxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxyy[i] = 4.0 * ts_xxxy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxyz[i] = 4.0 * ts_xxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxzz[i] = 4.0 * ts_xxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyyy[i] = 4.0 * ts_xxxy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyyz[i] = 4.0 * ts_xxxy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyzz[i] = 4.0 * ts_xxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xzzz[i] = 4.0 * ts_xxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyyy[i] = 4.0 * ts_xxxy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyyz[i] = 4.0 * ts_xxxy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyzz[i] = 4.0 * ts_xxxy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yzzz[i] = 4.0 * ts_xxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_zzzz[i] = 4.0 * ts_xxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 375-390 components of targeted buffer : HG

    auto gs_y_xxxyz_xxxx = pbuffer.data(idx_g_hg + 375);

    auto gs_y_xxxyz_xxxy = pbuffer.data(idx_g_hg + 376);

    auto gs_y_xxxyz_xxxz = pbuffer.data(idx_g_hg + 377);

    auto gs_y_xxxyz_xxyy = pbuffer.data(idx_g_hg + 378);

    auto gs_y_xxxyz_xxyz = pbuffer.data(idx_g_hg + 379);

    auto gs_y_xxxyz_xxzz = pbuffer.data(idx_g_hg + 380);

    auto gs_y_xxxyz_xyyy = pbuffer.data(idx_g_hg + 381);

    auto gs_y_xxxyz_xyyz = pbuffer.data(idx_g_hg + 382);

    auto gs_y_xxxyz_xyzz = pbuffer.data(idx_g_hg + 383);

    auto gs_y_xxxyz_xzzz = pbuffer.data(idx_g_hg + 384);

    auto gs_y_xxxyz_yyyy = pbuffer.data(idx_g_hg + 385);

    auto gs_y_xxxyz_yyyz = pbuffer.data(idx_g_hg + 386);

    auto gs_y_xxxyz_yyzz = pbuffer.data(idx_g_hg + 387);

    auto gs_y_xxxyz_yzzz = pbuffer.data(idx_g_hg + 388);

    auto gs_y_xxxyz_zzzz = pbuffer.data(idx_g_hg + 389);

    #pragma omp simd aligned(gc_y, gs_y_xxxyz_xxxx, gs_y_xxxyz_xxxy, gs_y_xxxyz_xxxz, gs_y_xxxyz_xxyy, gs_y_xxxyz_xxyz, gs_y_xxxyz_xxzz, gs_y_xxxyz_xyyy, gs_y_xxxyz_xyyz, gs_y_xxxyz_xyzz, gs_y_xxxyz_xzzz, gs_y_xxxyz_yyyy, gs_y_xxxyz_yyyz, gs_y_xxxyz_yyzz, gs_y_xxxyz_yzzz, gs_y_xxxyz_zzzz, ts_xxxyz_xxx, ts_xxxyz_xxxx, ts_xxxyz_xxxy, ts_xxxyz_xxxz, ts_xxxyz_xxy, ts_xxxyz_xxyy, ts_xxxyz_xxyz, ts_xxxyz_xxz, ts_xxxyz_xxzz, ts_xxxyz_xyy, ts_xxxyz_xyyy, ts_xxxyz_xyyz, ts_xxxyz_xyz, ts_xxxyz_xyzz, ts_xxxyz_xzz, ts_xxxyz_xzzz, ts_xxxyz_yyy, ts_xxxyz_yyyy, ts_xxxyz_yyyz, ts_xxxyz_yyz, ts_xxxyz_yyzz, ts_xxxyz_yzz, ts_xxxyz_yzzz, ts_xxxyz_zzz, ts_xxxyz_zzzz, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxzz, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyzz, ts_xxxz_xzzz, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyzz, ts_xxxz_yzzz, ts_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyz_xxxx[i] = 2.0 * ts_xxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxy[i] = 2.0 * ts_xxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxz[i] = 2.0 * ts_xxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxyy[i] = 2.0 * ts_xxxz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxyz[i] = 2.0 * ts_xxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxzz[i] = 2.0 * ts_xxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyyy[i] = 2.0 * ts_xxxz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyyz[i] = 2.0 * ts_xxxz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyzz[i] = 2.0 * ts_xxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xzzz[i] = 2.0 * ts_xxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyyy[i] = 2.0 * ts_xxxz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyyz[i] = 2.0 * ts_xxxz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyzz[i] = 2.0 * ts_xxxz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yzzz[i] = 2.0 * ts_xxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_zzzz[i] = 2.0 * ts_xxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 390-405 components of targeted buffer : HG

    auto gs_y_xxxzz_xxxx = pbuffer.data(idx_g_hg + 390);

    auto gs_y_xxxzz_xxxy = pbuffer.data(idx_g_hg + 391);

    auto gs_y_xxxzz_xxxz = pbuffer.data(idx_g_hg + 392);

    auto gs_y_xxxzz_xxyy = pbuffer.data(idx_g_hg + 393);

    auto gs_y_xxxzz_xxyz = pbuffer.data(idx_g_hg + 394);

    auto gs_y_xxxzz_xxzz = pbuffer.data(idx_g_hg + 395);

    auto gs_y_xxxzz_xyyy = pbuffer.data(idx_g_hg + 396);

    auto gs_y_xxxzz_xyyz = pbuffer.data(idx_g_hg + 397);

    auto gs_y_xxxzz_xyzz = pbuffer.data(idx_g_hg + 398);

    auto gs_y_xxxzz_xzzz = pbuffer.data(idx_g_hg + 399);

    auto gs_y_xxxzz_yyyy = pbuffer.data(idx_g_hg + 400);

    auto gs_y_xxxzz_yyyz = pbuffer.data(idx_g_hg + 401);

    auto gs_y_xxxzz_yyzz = pbuffer.data(idx_g_hg + 402);

    auto gs_y_xxxzz_yzzz = pbuffer.data(idx_g_hg + 403);

    auto gs_y_xxxzz_zzzz = pbuffer.data(idx_g_hg + 404);

    #pragma omp simd aligned(gc_y, gs_y_xxxzz_xxxx, gs_y_xxxzz_xxxy, gs_y_xxxzz_xxxz, gs_y_xxxzz_xxyy, gs_y_xxxzz_xxyz, gs_y_xxxzz_xxzz, gs_y_xxxzz_xyyy, gs_y_xxxzz_xyyz, gs_y_xxxzz_xyzz, gs_y_xxxzz_xzzz, gs_y_xxxzz_yyyy, gs_y_xxxzz_yyyz, gs_y_xxxzz_yyzz, gs_y_xxxzz_yzzz, gs_y_xxxzz_zzzz, ts_xxxzz_xxx, ts_xxxzz_xxxx, ts_xxxzz_xxxy, ts_xxxzz_xxxz, ts_xxxzz_xxy, ts_xxxzz_xxyy, ts_xxxzz_xxyz, ts_xxxzz_xxz, ts_xxxzz_xxzz, ts_xxxzz_xyy, ts_xxxzz_xyyy, ts_xxxzz_xyyz, ts_xxxzz_xyz, ts_xxxzz_xyzz, ts_xxxzz_xzz, ts_xxxzz_xzzz, ts_xxxzz_yyy, ts_xxxzz_yyyy, ts_xxxzz_yyyz, ts_xxxzz_yyz, ts_xxxzz_yyzz, ts_xxxzz_yzz, ts_xxxzz_yzzz, ts_xxxzz_zzz, ts_xxxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxzz_xxxx[i] = 2.0 * ts_xxxzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxy[i] = 2.0 * ts_xxxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxz[i] = 2.0 * ts_xxxzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxyy[i] = 4.0 * ts_xxxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxyz[i] = 2.0 * ts_xxxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxzz[i] = 2.0 * ts_xxxzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyyy[i] = 6.0 * ts_xxxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyyz[i] = 4.0 * ts_xxxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyzz[i] = 2.0 * ts_xxxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xzzz[i] = 2.0 * ts_xxxzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyyy[i] = 8.0 * ts_xxxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyyz[i] = 6.0 * ts_xxxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyzz[i] = 4.0 * ts_xxxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yzzz[i] = 2.0 * ts_xxxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_zzzz[i] = 2.0 * ts_xxxzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 405-420 components of targeted buffer : HG

    auto gs_y_xxyyy_xxxx = pbuffer.data(idx_g_hg + 405);

    auto gs_y_xxyyy_xxxy = pbuffer.data(idx_g_hg + 406);

    auto gs_y_xxyyy_xxxz = pbuffer.data(idx_g_hg + 407);

    auto gs_y_xxyyy_xxyy = pbuffer.data(idx_g_hg + 408);

    auto gs_y_xxyyy_xxyz = pbuffer.data(idx_g_hg + 409);

    auto gs_y_xxyyy_xxzz = pbuffer.data(idx_g_hg + 410);

    auto gs_y_xxyyy_xyyy = pbuffer.data(idx_g_hg + 411);

    auto gs_y_xxyyy_xyyz = pbuffer.data(idx_g_hg + 412);

    auto gs_y_xxyyy_xyzz = pbuffer.data(idx_g_hg + 413);

    auto gs_y_xxyyy_xzzz = pbuffer.data(idx_g_hg + 414);

    auto gs_y_xxyyy_yyyy = pbuffer.data(idx_g_hg + 415);

    auto gs_y_xxyyy_yyyz = pbuffer.data(idx_g_hg + 416);

    auto gs_y_xxyyy_yyzz = pbuffer.data(idx_g_hg + 417);

    auto gs_y_xxyyy_yzzz = pbuffer.data(idx_g_hg + 418);

    auto gs_y_xxyyy_zzzz = pbuffer.data(idx_g_hg + 419);

    #pragma omp simd aligned(gc_y, gs_y_xxyyy_xxxx, gs_y_xxyyy_xxxy, gs_y_xxyyy_xxxz, gs_y_xxyyy_xxyy, gs_y_xxyyy_xxyz, gs_y_xxyyy_xxzz, gs_y_xxyyy_xyyy, gs_y_xxyyy_xyyz, gs_y_xxyyy_xyzz, gs_y_xxyyy_xzzz, gs_y_xxyyy_yyyy, gs_y_xxyyy_yyyz, gs_y_xxyyy_yyzz, gs_y_xxyyy_yzzz, gs_y_xxyyy_zzzz, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxzz, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyzz, ts_xxyy_xzzz, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyzz, ts_xxyy_yzzz, ts_xxyy_zzzz, ts_xxyyy_xxx, ts_xxyyy_xxxx, ts_xxyyy_xxxy, ts_xxyyy_xxxz, ts_xxyyy_xxy, ts_xxyyy_xxyy, ts_xxyyy_xxyz, ts_xxyyy_xxz, ts_xxyyy_xxzz, ts_xxyyy_xyy, ts_xxyyy_xyyy, ts_xxyyy_xyyz, ts_xxyyy_xyz, ts_xxyyy_xyzz, ts_xxyyy_xzz, ts_xxyyy_xzzz, ts_xxyyy_yyy, ts_xxyyy_yyyy, ts_xxyyy_yyyz, ts_xxyyy_yyz, ts_xxyyy_yyzz, ts_xxyyy_yzz, ts_xxyyy_yzzz, ts_xxyyy_zzz, ts_xxyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyy_xxxx[i] = 6.0 * ts_xxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxy[i] = 6.0 * ts_xxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxz[i] = 6.0 * ts_xxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxyy[i] = 6.0 * ts_xxyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxyz[i] = 6.0 * ts_xxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxzz[i] = 6.0 * ts_xxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyyy[i] = 6.0 * ts_xxyy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyyz[i] = 6.0 * ts_xxyy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyzz[i] = 6.0 * ts_xxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xzzz[i] = 6.0 * ts_xxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyyy[i] = 6.0 * ts_xxyy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyyz[i] = 6.0 * ts_xxyy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyzz[i] = 6.0 * ts_xxyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yzzz[i] = 6.0 * ts_xxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_zzzz[i] = 6.0 * ts_xxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 420-435 components of targeted buffer : HG

    auto gs_y_xxyyz_xxxx = pbuffer.data(idx_g_hg + 420);

    auto gs_y_xxyyz_xxxy = pbuffer.data(idx_g_hg + 421);

    auto gs_y_xxyyz_xxxz = pbuffer.data(idx_g_hg + 422);

    auto gs_y_xxyyz_xxyy = pbuffer.data(idx_g_hg + 423);

    auto gs_y_xxyyz_xxyz = pbuffer.data(idx_g_hg + 424);

    auto gs_y_xxyyz_xxzz = pbuffer.data(idx_g_hg + 425);

    auto gs_y_xxyyz_xyyy = pbuffer.data(idx_g_hg + 426);

    auto gs_y_xxyyz_xyyz = pbuffer.data(idx_g_hg + 427);

    auto gs_y_xxyyz_xyzz = pbuffer.data(idx_g_hg + 428);

    auto gs_y_xxyyz_xzzz = pbuffer.data(idx_g_hg + 429);

    auto gs_y_xxyyz_yyyy = pbuffer.data(idx_g_hg + 430);

    auto gs_y_xxyyz_yyyz = pbuffer.data(idx_g_hg + 431);

    auto gs_y_xxyyz_yyzz = pbuffer.data(idx_g_hg + 432);

    auto gs_y_xxyyz_yzzz = pbuffer.data(idx_g_hg + 433);

    auto gs_y_xxyyz_zzzz = pbuffer.data(idx_g_hg + 434);

    #pragma omp simd aligned(gc_y, gs_y_xxyyz_xxxx, gs_y_xxyyz_xxxy, gs_y_xxyyz_xxxz, gs_y_xxyyz_xxyy, gs_y_xxyyz_xxyz, gs_y_xxyyz_xxzz, gs_y_xxyyz_xyyy, gs_y_xxyyz_xyyz, gs_y_xxyyz_xyzz, gs_y_xxyyz_xzzz, gs_y_xxyyz_yyyy, gs_y_xxyyz_yyyz, gs_y_xxyyz_yyzz, gs_y_xxyyz_yzzz, gs_y_xxyyz_zzzz, ts_xxyyz_xxx, ts_xxyyz_xxxx, ts_xxyyz_xxxy, ts_xxyyz_xxxz, ts_xxyyz_xxy, ts_xxyyz_xxyy, ts_xxyyz_xxyz, ts_xxyyz_xxz, ts_xxyyz_xxzz, ts_xxyyz_xyy, ts_xxyyz_xyyy, ts_xxyyz_xyyz, ts_xxyyz_xyz, ts_xxyyz_xyzz, ts_xxyyz_xzz, ts_xxyyz_xzzz, ts_xxyyz_yyy, ts_xxyyz_yyyy, ts_xxyyz_yyyz, ts_xxyyz_yyz, ts_xxyyz_yyzz, ts_xxyyz_yzz, ts_xxyyz_yzzz, ts_xxyyz_zzz, ts_xxyyz_zzzz, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxzz, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyzz, ts_xxyz_xzzz, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyzz, ts_xxyz_yzzz, ts_xxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyz_xxxx[i] = 4.0 * ts_xxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxy[i] = 4.0 * ts_xxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxz[i] = 4.0 * ts_xxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxyy[i] = 4.0 * ts_xxyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxyz[i] = 4.0 * ts_xxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxzz[i] = 4.0 * ts_xxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyyy[i] = 4.0 * ts_xxyz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyyz[i] = 4.0 * ts_xxyz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyzz[i] = 4.0 * ts_xxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xzzz[i] = 4.0 * ts_xxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyyy[i] = 4.0 * ts_xxyz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyyz[i] = 4.0 * ts_xxyz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyzz[i] = 4.0 * ts_xxyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yzzz[i] = 4.0 * ts_xxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_zzzz[i] = 4.0 * ts_xxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 435-450 components of targeted buffer : HG

    auto gs_y_xxyzz_xxxx = pbuffer.data(idx_g_hg + 435);

    auto gs_y_xxyzz_xxxy = pbuffer.data(idx_g_hg + 436);

    auto gs_y_xxyzz_xxxz = pbuffer.data(idx_g_hg + 437);

    auto gs_y_xxyzz_xxyy = pbuffer.data(idx_g_hg + 438);

    auto gs_y_xxyzz_xxyz = pbuffer.data(idx_g_hg + 439);

    auto gs_y_xxyzz_xxzz = pbuffer.data(idx_g_hg + 440);

    auto gs_y_xxyzz_xyyy = pbuffer.data(idx_g_hg + 441);

    auto gs_y_xxyzz_xyyz = pbuffer.data(idx_g_hg + 442);

    auto gs_y_xxyzz_xyzz = pbuffer.data(idx_g_hg + 443);

    auto gs_y_xxyzz_xzzz = pbuffer.data(idx_g_hg + 444);

    auto gs_y_xxyzz_yyyy = pbuffer.data(idx_g_hg + 445);

    auto gs_y_xxyzz_yyyz = pbuffer.data(idx_g_hg + 446);

    auto gs_y_xxyzz_yyzz = pbuffer.data(idx_g_hg + 447);

    auto gs_y_xxyzz_yzzz = pbuffer.data(idx_g_hg + 448);

    auto gs_y_xxyzz_zzzz = pbuffer.data(idx_g_hg + 449);

    #pragma omp simd aligned(gc_y, gs_y_xxyzz_xxxx, gs_y_xxyzz_xxxy, gs_y_xxyzz_xxxz, gs_y_xxyzz_xxyy, gs_y_xxyzz_xxyz, gs_y_xxyzz_xxzz, gs_y_xxyzz_xyyy, gs_y_xxyzz_xyyz, gs_y_xxyzz_xyzz, gs_y_xxyzz_xzzz, gs_y_xxyzz_yyyy, gs_y_xxyzz_yyyz, gs_y_xxyzz_yyzz, gs_y_xxyzz_yzzz, gs_y_xxyzz_zzzz, ts_xxyzz_xxx, ts_xxyzz_xxxx, ts_xxyzz_xxxy, ts_xxyzz_xxxz, ts_xxyzz_xxy, ts_xxyzz_xxyy, ts_xxyzz_xxyz, ts_xxyzz_xxz, ts_xxyzz_xxzz, ts_xxyzz_xyy, ts_xxyzz_xyyy, ts_xxyzz_xyyz, ts_xxyzz_xyz, ts_xxyzz_xyzz, ts_xxyzz_xzz, ts_xxyzz_xzzz, ts_xxyzz_yyy, ts_xxyzz_yyyy, ts_xxyzz_yyyz, ts_xxyzz_yyz, ts_xxyzz_yyzz, ts_xxyzz_yzz, ts_xxyzz_yzzz, ts_xxyzz_zzz, ts_xxyzz_zzzz, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxzz, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyzz, ts_xxzz_xzzz, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyzz, ts_xxzz_yzzz, ts_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyzz_xxxx[i] = 2.0 * ts_xxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxy[i] = 2.0 * ts_xxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxz[i] = 2.0 * ts_xxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxyy[i] = 2.0 * ts_xxzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxyz[i] = 2.0 * ts_xxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxzz[i] = 2.0 * ts_xxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyyy[i] = 2.0 * ts_xxzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyyz[i] = 2.0 * ts_xxzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyzz[i] = 2.0 * ts_xxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xzzz[i] = 2.0 * ts_xxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyyy[i] = 2.0 * ts_xxzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyyz[i] = 2.0 * ts_xxzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyzz[i] = 2.0 * ts_xxzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yzzz[i] = 2.0 * ts_xxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_zzzz[i] = 2.0 * ts_xxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 450-465 components of targeted buffer : HG

    auto gs_y_xxzzz_xxxx = pbuffer.data(idx_g_hg + 450);

    auto gs_y_xxzzz_xxxy = pbuffer.data(idx_g_hg + 451);

    auto gs_y_xxzzz_xxxz = pbuffer.data(idx_g_hg + 452);

    auto gs_y_xxzzz_xxyy = pbuffer.data(idx_g_hg + 453);

    auto gs_y_xxzzz_xxyz = pbuffer.data(idx_g_hg + 454);

    auto gs_y_xxzzz_xxzz = pbuffer.data(idx_g_hg + 455);

    auto gs_y_xxzzz_xyyy = pbuffer.data(idx_g_hg + 456);

    auto gs_y_xxzzz_xyyz = pbuffer.data(idx_g_hg + 457);

    auto gs_y_xxzzz_xyzz = pbuffer.data(idx_g_hg + 458);

    auto gs_y_xxzzz_xzzz = pbuffer.data(idx_g_hg + 459);

    auto gs_y_xxzzz_yyyy = pbuffer.data(idx_g_hg + 460);

    auto gs_y_xxzzz_yyyz = pbuffer.data(idx_g_hg + 461);

    auto gs_y_xxzzz_yyzz = pbuffer.data(idx_g_hg + 462);

    auto gs_y_xxzzz_yzzz = pbuffer.data(idx_g_hg + 463);

    auto gs_y_xxzzz_zzzz = pbuffer.data(idx_g_hg + 464);

    #pragma omp simd aligned(gc_y, gs_y_xxzzz_xxxx, gs_y_xxzzz_xxxy, gs_y_xxzzz_xxxz, gs_y_xxzzz_xxyy, gs_y_xxzzz_xxyz, gs_y_xxzzz_xxzz, gs_y_xxzzz_xyyy, gs_y_xxzzz_xyyz, gs_y_xxzzz_xyzz, gs_y_xxzzz_xzzz, gs_y_xxzzz_yyyy, gs_y_xxzzz_yyyz, gs_y_xxzzz_yyzz, gs_y_xxzzz_yzzz, gs_y_xxzzz_zzzz, ts_xxzzz_xxx, ts_xxzzz_xxxx, ts_xxzzz_xxxy, ts_xxzzz_xxxz, ts_xxzzz_xxy, ts_xxzzz_xxyy, ts_xxzzz_xxyz, ts_xxzzz_xxz, ts_xxzzz_xxzz, ts_xxzzz_xyy, ts_xxzzz_xyyy, ts_xxzzz_xyyz, ts_xxzzz_xyz, ts_xxzzz_xyzz, ts_xxzzz_xzz, ts_xxzzz_xzzz, ts_xxzzz_yyy, ts_xxzzz_yyyy, ts_xxzzz_yyyz, ts_xxzzz_yyz, ts_xxzzz_yyzz, ts_xxzzz_yzz, ts_xxzzz_yzzz, ts_xxzzz_zzz, ts_xxzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzzz_xxxx[i] = 2.0 * ts_xxzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxy[i] = 2.0 * ts_xxzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxz[i] = 2.0 * ts_xxzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxyy[i] = 4.0 * ts_xxzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxyz[i] = 2.0 * ts_xxzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxzz[i] = 2.0 * ts_xxzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyyy[i] = 6.0 * ts_xxzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyyz[i] = 4.0 * ts_xxzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyzz[i] = 2.0 * ts_xxzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xzzz[i] = 2.0 * ts_xxzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyyy[i] = 8.0 * ts_xxzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyyz[i] = 6.0 * ts_xxzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyzz[i] = 4.0 * ts_xxzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yzzz[i] = 2.0 * ts_xxzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_zzzz[i] = 2.0 * ts_xxzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 465-480 components of targeted buffer : HG

    auto gs_y_xyyyy_xxxx = pbuffer.data(idx_g_hg + 465);

    auto gs_y_xyyyy_xxxy = pbuffer.data(idx_g_hg + 466);

    auto gs_y_xyyyy_xxxz = pbuffer.data(idx_g_hg + 467);

    auto gs_y_xyyyy_xxyy = pbuffer.data(idx_g_hg + 468);

    auto gs_y_xyyyy_xxyz = pbuffer.data(idx_g_hg + 469);

    auto gs_y_xyyyy_xxzz = pbuffer.data(idx_g_hg + 470);

    auto gs_y_xyyyy_xyyy = pbuffer.data(idx_g_hg + 471);

    auto gs_y_xyyyy_xyyz = pbuffer.data(idx_g_hg + 472);

    auto gs_y_xyyyy_xyzz = pbuffer.data(idx_g_hg + 473);

    auto gs_y_xyyyy_xzzz = pbuffer.data(idx_g_hg + 474);

    auto gs_y_xyyyy_yyyy = pbuffer.data(idx_g_hg + 475);

    auto gs_y_xyyyy_yyyz = pbuffer.data(idx_g_hg + 476);

    auto gs_y_xyyyy_yyzz = pbuffer.data(idx_g_hg + 477);

    auto gs_y_xyyyy_yzzz = pbuffer.data(idx_g_hg + 478);

    auto gs_y_xyyyy_zzzz = pbuffer.data(idx_g_hg + 479);

    #pragma omp simd aligned(gc_y, gs_y_xyyyy_xxxx, gs_y_xyyyy_xxxy, gs_y_xyyyy_xxxz, gs_y_xyyyy_xxyy, gs_y_xyyyy_xxyz, gs_y_xyyyy_xxzz, gs_y_xyyyy_xyyy, gs_y_xyyyy_xyyz, gs_y_xyyyy_xyzz, gs_y_xyyyy_xzzz, gs_y_xyyyy_yyyy, gs_y_xyyyy_yyyz, gs_y_xyyyy_yyzz, gs_y_xyyyy_yzzz, gs_y_xyyyy_zzzz, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxzz, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyzz, ts_xyyy_xzzz, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyzz, ts_xyyy_yzzz, ts_xyyy_zzzz, ts_xyyyy_xxx, ts_xyyyy_xxxx, ts_xyyyy_xxxy, ts_xyyyy_xxxz, ts_xyyyy_xxy, ts_xyyyy_xxyy, ts_xyyyy_xxyz, ts_xyyyy_xxz, ts_xyyyy_xxzz, ts_xyyyy_xyy, ts_xyyyy_xyyy, ts_xyyyy_xyyz, ts_xyyyy_xyz, ts_xyyyy_xyzz, ts_xyyyy_xzz, ts_xyyyy_xzzz, ts_xyyyy_yyy, ts_xyyyy_yyyy, ts_xyyyy_yyyz, ts_xyyyy_yyz, ts_xyyyy_yyzz, ts_xyyyy_yzz, ts_xyyyy_yzzz, ts_xyyyy_zzz, ts_xyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyy_xxxx[i] = 8.0 * ts_xyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxy[i] = 8.0 * ts_xyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxz[i] = 8.0 * ts_xyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxyy[i] = 8.0 * ts_xyyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxyz[i] = 8.0 * ts_xyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxzz[i] = 8.0 * ts_xyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyyy[i] = 8.0 * ts_xyyy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyyz[i] = 8.0 * ts_xyyy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyzz[i] = 8.0 * ts_xyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xzzz[i] = 8.0 * ts_xyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyyy[i] = 8.0 * ts_xyyy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyyz[i] = 8.0 * ts_xyyy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyzz[i] = 8.0 * ts_xyyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yzzz[i] = 8.0 * ts_xyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_zzzz[i] = 8.0 * ts_xyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 480-495 components of targeted buffer : HG

    auto gs_y_xyyyz_xxxx = pbuffer.data(idx_g_hg + 480);

    auto gs_y_xyyyz_xxxy = pbuffer.data(idx_g_hg + 481);

    auto gs_y_xyyyz_xxxz = pbuffer.data(idx_g_hg + 482);

    auto gs_y_xyyyz_xxyy = pbuffer.data(idx_g_hg + 483);

    auto gs_y_xyyyz_xxyz = pbuffer.data(idx_g_hg + 484);

    auto gs_y_xyyyz_xxzz = pbuffer.data(idx_g_hg + 485);

    auto gs_y_xyyyz_xyyy = pbuffer.data(idx_g_hg + 486);

    auto gs_y_xyyyz_xyyz = pbuffer.data(idx_g_hg + 487);

    auto gs_y_xyyyz_xyzz = pbuffer.data(idx_g_hg + 488);

    auto gs_y_xyyyz_xzzz = pbuffer.data(idx_g_hg + 489);

    auto gs_y_xyyyz_yyyy = pbuffer.data(idx_g_hg + 490);

    auto gs_y_xyyyz_yyyz = pbuffer.data(idx_g_hg + 491);

    auto gs_y_xyyyz_yyzz = pbuffer.data(idx_g_hg + 492);

    auto gs_y_xyyyz_yzzz = pbuffer.data(idx_g_hg + 493);

    auto gs_y_xyyyz_zzzz = pbuffer.data(idx_g_hg + 494);

    #pragma omp simd aligned(gc_y, gs_y_xyyyz_xxxx, gs_y_xyyyz_xxxy, gs_y_xyyyz_xxxz, gs_y_xyyyz_xxyy, gs_y_xyyyz_xxyz, gs_y_xyyyz_xxzz, gs_y_xyyyz_xyyy, gs_y_xyyyz_xyyz, gs_y_xyyyz_xyzz, gs_y_xyyyz_xzzz, gs_y_xyyyz_yyyy, gs_y_xyyyz_yyyz, gs_y_xyyyz_yyzz, gs_y_xyyyz_yzzz, gs_y_xyyyz_zzzz, ts_xyyyz_xxx, ts_xyyyz_xxxx, ts_xyyyz_xxxy, ts_xyyyz_xxxz, ts_xyyyz_xxy, ts_xyyyz_xxyy, ts_xyyyz_xxyz, ts_xyyyz_xxz, ts_xyyyz_xxzz, ts_xyyyz_xyy, ts_xyyyz_xyyy, ts_xyyyz_xyyz, ts_xyyyz_xyz, ts_xyyyz_xyzz, ts_xyyyz_xzz, ts_xyyyz_xzzz, ts_xyyyz_yyy, ts_xyyyz_yyyy, ts_xyyyz_yyyz, ts_xyyyz_yyz, ts_xyyyz_yyzz, ts_xyyyz_yzz, ts_xyyyz_yzzz, ts_xyyyz_zzz, ts_xyyyz_zzzz, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxzz, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyzz, ts_xyyz_xzzz, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyzz, ts_xyyz_yzzz, ts_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyz_xxxx[i] = 6.0 * ts_xyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxy[i] = 6.0 * ts_xyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxz[i] = 6.0 * ts_xyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxyy[i] = 6.0 * ts_xyyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxyz[i] = 6.0 * ts_xyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxzz[i] = 6.0 * ts_xyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyyy[i] = 6.0 * ts_xyyz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyyz[i] = 6.0 * ts_xyyz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyzz[i] = 6.0 * ts_xyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xzzz[i] = 6.0 * ts_xyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyyy[i] = 6.0 * ts_xyyz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyyz[i] = 6.0 * ts_xyyz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyzz[i] = 6.0 * ts_xyyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yzzz[i] = 6.0 * ts_xyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_zzzz[i] = 6.0 * ts_xyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 495-510 components of targeted buffer : HG

    auto gs_y_xyyzz_xxxx = pbuffer.data(idx_g_hg + 495);

    auto gs_y_xyyzz_xxxy = pbuffer.data(idx_g_hg + 496);

    auto gs_y_xyyzz_xxxz = pbuffer.data(idx_g_hg + 497);

    auto gs_y_xyyzz_xxyy = pbuffer.data(idx_g_hg + 498);

    auto gs_y_xyyzz_xxyz = pbuffer.data(idx_g_hg + 499);

    auto gs_y_xyyzz_xxzz = pbuffer.data(idx_g_hg + 500);

    auto gs_y_xyyzz_xyyy = pbuffer.data(idx_g_hg + 501);

    auto gs_y_xyyzz_xyyz = pbuffer.data(idx_g_hg + 502);

    auto gs_y_xyyzz_xyzz = pbuffer.data(idx_g_hg + 503);

    auto gs_y_xyyzz_xzzz = pbuffer.data(idx_g_hg + 504);

    auto gs_y_xyyzz_yyyy = pbuffer.data(idx_g_hg + 505);

    auto gs_y_xyyzz_yyyz = pbuffer.data(idx_g_hg + 506);

    auto gs_y_xyyzz_yyzz = pbuffer.data(idx_g_hg + 507);

    auto gs_y_xyyzz_yzzz = pbuffer.data(idx_g_hg + 508);

    auto gs_y_xyyzz_zzzz = pbuffer.data(idx_g_hg + 509);

    #pragma omp simd aligned(gc_y, gs_y_xyyzz_xxxx, gs_y_xyyzz_xxxy, gs_y_xyyzz_xxxz, gs_y_xyyzz_xxyy, gs_y_xyyzz_xxyz, gs_y_xyyzz_xxzz, gs_y_xyyzz_xyyy, gs_y_xyyzz_xyyz, gs_y_xyyzz_xyzz, gs_y_xyyzz_xzzz, gs_y_xyyzz_yyyy, gs_y_xyyzz_yyyz, gs_y_xyyzz_yyzz, gs_y_xyyzz_yzzz, gs_y_xyyzz_zzzz, ts_xyyzz_xxx, ts_xyyzz_xxxx, ts_xyyzz_xxxy, ts_xyyzz_xxxz, ts_xyyzz_xxy, ts_xyyzz_xxyy, ts_xyyzz_xxyz, ts_xyyzz_xxz, ts_xyyzz_xxzz, ts_xyyzz_xyy, ts_xyyzz_xyyy, ts_xyyzz_xyyz, ts_xyyzz_xyz, ts_xyyzz_xyzz, ts_xyyzz_xzz, ts_xyyzz_xzzz, ts_xyyzz_yyy, ts_xyyzz_yyyy, ts_xyyzz_yyyz, ts_xyyzz_yyz, ts_xyyzz_yyzz, ts_xyyzz_yzz, ts_xyyzz_yzzz, ts_xyyzz_zzz, ts_xyyzz_zzzz, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxzz, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyzz, ts_xyzz_xzzz, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyzz, ts_xyzz_yzzz, ts_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyzz_xxxx[i] = 4.0 * ts_xyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxy[i] = 4.0 * ts_xyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxz[i] = 4.0 * ts_xyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxyy[i] = 4.0 * ts_xyzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxyz[i] = 4.0 * ts_xyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxzz[i] = 4.0 * ts_xyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyyy[i] = 4.0 * ts_xyzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyyz[i] = 4.0 * ts_xyzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyzz[i] = 4.0 * ts_xyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xzzz[i] = 4.0 * ts_xyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyyy[i] = 4.0 * ts_xyzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyyz[i] = 4.0 * ts_xyzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyzz[i] = 4.0 * ts_xyzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yzzz[i] = 4.0 * ts_xyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_zzzz[i] = 4.0 * ts_xyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 510-525 components of targeted buffer : HG

    auto gs_y_xyzzz_xxxx = pbuffer.data(idx_g_hg + 510);

    auto gs_y_xyzzz_xxxy = pbuffer.data(idx_g_hg + 511);

    auto gs_y_xyzzz_xxxz = pbuffer.data(idx_g_hg + 512);

    auto gs_y_xyzzz_xxyy = pbuffer.data(idx_g_hg + 513);

    auto gs_y_xyzzz_xxyz = pbuffer.data(idx_g_hg + 514);

    auto gs_y_xyzzz_xxzz = pbuffer.data(idx_g_hg + 515);

    auto gs_y_xyzzz_xyyy = pbuffer.data(idx_g_hg + 516);

    auto gs_y_xyzzz_xyyz = pbuffer.data(idx_g_hg + 517);

    auto gs_y_xyzzz_xyzz = pbuffer.data(idx_g_hg + 518);

    auto gs_y_xyzzz_xzzz = pbuffer.data(idx_g_hg + 519);

    auto gs_y_xyzzz_yyyy = pbuffer.data(idx_g_hg + 520);

    auto gs_y_xyzzz_yyyz = pbuffer.data(idx_g_hg + 521);

    auto gs_y_xyzzz_yyzz = pbuffer.data(idx_g_hg + 522);

    auto gs_y_xyzzz_yzzz = pbuffer.data(idx_g_hg + 523);

    auto gs_y_xyzzz_zzzz = pbuffer.data(idx_g_hg + 524);

    #pragma omp simd aligned(gc_y, gs_y_xyzzz_xxxx, gs_y_xyzzz_xxxy, gs_y_xyzzz_xxxz, gs_y_xyzzz_xxyy, gs_y_xyzzz_xxyz, gs_y_xyzzz_xxzz, gs_y_xyzzz_xyyy, gs_y_xyzzz_xyyz, gs_y_xyzzz_xyzz, gs_y_xyzzz_xzzz, gs_y_xyzzz_yyyy, gs_y_xyzzz_yyyz, gs_y_xyzzz_yyzz, gs_y_xyzzz_yzzz, gs_y_xyzzz_zzzz, ts_xyzzz_xxx, ts_xyzzz_xxxx, ts_xyzzz_xxxy, ts_xyzzz_xxxz, ts_xyzzz_xxy, ts_xyzzz_xxyy, ts_xyzzz_xxyz, ts_xyzzz_xxz, ts_xyzzz_xxzz, ts_xyzzz_xyy, ts_xyzzz_xyyy, ts_xyzzz_xyyz, ts_xyzzz_xyz, ts_xyzzz_xyzz, ts_xyzzz_xzz, ts_xyzzz_xzzz, ts_xyzzz_yyy, ts_xyzzz_yyyy, ts_xyzzz_yyyz, ts_xyzzz_yyz, ts_xyzzz_yyzz, ts_xyzzz_yzz, ts_xyzzz_yzzz, ts_xyzzz_zzz, ts_xyzzz_zzzz, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxzz, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyzz, ts_xzzz_xzzz, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyzz, ts_xzzz_yzzz, ts_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzzz_xxxx[i] = 2.0 * ts_xzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxy[i] = 2.0 * ts_xzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxz[i] = 2.0 * ts_xzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxyy[i] = 2.0 * ts_xzzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxyz[i] = 2.0 * ts_xzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxzz[i] = 2.0 * ts_xzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyyy[i] = 2.0 * ts_xzzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyyz[i] = 2.0 * ts_xzzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyzz[i] = 2.0 * ts_xzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xzzz[i] = 2.0 * ts_xzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyyy[i] = 2.0 * ts_xzzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyyz[i] = 2.0 * ts_xzzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyzz[i] = 2.0 * ts_xzzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yzzz[i] = 2.0 * ts_xzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_zzzz[i] = 2.0 * ts_xzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 525-540 components of targeted buffer : HG

    auto gs_y_xzzzz_xxxx = pbuffer.data(idx_g_hg + 525);

    auto gs_y_xzzzz_xxxy = pbuffer.data(idx_g_hg + 526);

    auto gs_y_xzzzz_xxxz = pbuffer.data(idx_g_hg + 527);

    auto gs_y_xzzzz_xxyy = pbuffer.data(idx_g_hg + 528);

    auto gs_y_xzzzz_xxyz = pbuffer.data(idx_g_hg + 529);

    auto gs_y_xzzzz_xxzz = pbuffer.data(idx_g_hg + 530);

    auto gs_y_xzzzz_xyyy = pbuffer.data(idx_g_hg + 531);

    auto gs_y_xzzzz_xyyz = pbuffer.data(idx_g_hg + 532);

    auto gs_y_xzzzz_xyzz = pbuffer.data(idx_g_hg + 533);

    auto gs_y_xzzzz_xzzz = pbuffer.data(idx_g_hg + 534);

    auto gs_y_xzzzz_yyyy = pbuffer.data(idx_g_hg + 535);

    auto gs_y_xzzzz_yyyz = pbuffer.data(idx_g_hg + 536);

    auto gs_y_xzzzz_yyzz = pbuffer.data(idx_g_hg + 537);

    auto gs_y_xzzzz_yzzz = pbuffer.data(idx_g_hg + 538);

    auto gs_y_xzzzz_zzzz = pbuffer.data(idx_g_hg + 539);

    #pragma omp simd aligned(gc_y, gs_y_xzzzz_xxxx, gs_y_xzzzz_xxxy, gs_y_xzzzz_xxxz, gs_y_xzzzz_xxyy, gs_y_xzzzz_xxyz, gs_y_xzzzz_xxzz, gs_y_xzzzz_xyyy, gs_y_xzzzz_xyyz, gs_y_xzzzz_xyzz, gs_y_xzzzz_xzzz, gs_y_xzzzz_yyyy, gs_y_xzzzz_yyyz, gs_y_xzzzz_yyzz, gs_y_xzzzz_yzzz, gs_y_xzzzz_zzzz, ts_xzzzz_xxx, ts_xzzzz_xxxx, ts_xzzzz_xxxy, ts_xzzzz_xxxz, ts_xzzzz_xxy, ts_xzzzz_xxyy, ts_xzzzz_xxyz, ts_xzzzz_xxz, ts_xzzzz_xxzz, ts_xzzzz_xyy, ts_xzzzz_xyyy, ts_xzzzz_xyyz, ts_xzzzz_xyz, ts_xzzzz_xyzz, ts_xzzzz_xzz, ts_xzzzz_xzzz, ts_xzzzz_yyy, ts_xzzzz_yyyy, ts_xzzzz_yyyz, ts_xzzzz_yyz, ts_xzzzz_yyzz, ts_xzzzz_yzz, ts_xzzzz_yzzz, ts_xzzzz_zzz, ts_xzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzzz_xxxx[i] = 2.0 * ts_xzzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxy[i] = 2.0 * ts_xzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxz[i] = 2.0 * ts_xzzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxyy[i] = 4.0 * ts_xzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxyz[i] = 2.0 * ts_xzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxzz[i] = 2.0 * ts_xzzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyyy[i] = 6.0 * ts_xzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyyz[i] = 4.0 * ts_xzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyzz[i] = 2.0 * ts_xzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xzzz[i] = 2.0 * ts_xzzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyyy[i] = 8.0 * ts_xzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyyz[i] = 6.0 * ts_xzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyzz[i] = 4.0 * ts_xzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yzzz[i] = 2.0 * ts_xzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_zzzz[i] = 2.0 * ts_xzzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 540-555 components of targeted buffer : HG

    auto gs_y_yyyyy_xxxx = pbuffer.data(idx_g_hg + 540);

    auto gs_y_yyyyy_xxxy = pbuffer.data(idx_g_hg + 541);

    auto gs_y_yyyyy_xxxz = pbuffer.data(idx_g_hg + 542);

    auto gs_y_yyyyy_xxyy = pbuffer.data(idx_g_hg + 543);

    auto gs_y_yyyyy_xxyz = pbuffer.data(idx_g_hg + 544);

    auto gs_y_yyyyy_xxzz = pbuffer.data(idx_g_hg + 545);

    auto gs_y_yyyyy_xyyy = pbuffer.data(idx_g_hg + 546);

    auto gs_y_yyyyy_xyyz = pbuffer.data(idx_g_hg + 547);

    auto gs_y_yyyyy_xyzz = pbuffer.data(idx_g_hg + 548);

    auto gs_y_yyyyy_xzzz = pbuffer.data(idx_g_hg + 549);

    auto gs_y_yyyyy_yyyy = pbuffer.data(idx_g_hg + 550);

    auto gs_y_yyyyy_yyyz = pbuffer.data(idx_g_hg + 551);

    auto gs_y_yyyyy_yyzz = pbuffer.data(idx_g_hg + 552);

    auto gs_y_yyyyy_yzzz = pbuffer.data(idx_g_hg + 553);

    auto gs_y_yyyyy_zzzz = pbuffer.data(idx_g_hg + 554);

    #pragma omp simd aligned(gc_y, gs_y_yyyyy_xxxx, gs_y_yyyyy_xxxy, gs_y_yyyyy_xxxz, gs_y_yyyyy_xxyy, gs_y_yyyyy_xxyz, gs_y_yyyyy_xxzz, gs_y_yyyyy_xyyy, gs_y_yyyyy_xyyz, gs_y_yyyyy_xyzz, gs_y_yyyyy_xzzz, gs_y_yyyyy_yyyy, gs_y_yyyyy_yyyz, gs_y_yyyyy_yyzz, gs_y_yyyyy_yzzz, gs_y_yyyyy_zzzz, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxzz, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyzz, ts_yyyy_xzzz, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyzz, ts_yyyy_yzzz, ts_yyyy_zzzz, ts_yyyyy_xxx, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxxz, ts_yyyyy_xxy, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xxz, ts_yyyyy_xxzz, ts_yyyyy_xyy, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyz, ts_yyyyy_xyzz, ts_yyyyy_xzz, ts_yyyyy_xzzz, ts_yyyyy_yyy, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyz, ts_yyyyy_yyzz, ts_yyyyy_yzz, ts_yyyyy_yzzz, ts_yyyyy_zzz, ts_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyy_xxxx[i] = 10.0 * ts_yyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxy[i] = 10.0 * ts_yyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxz[i] = 10.0 * ts_yyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxyy[i] = 10.0 * ts_yyyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxyz[i] = 10.0 * ts_yyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxzz[i] = 10.0 * ts_yyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyyy[i] = 10.0 * ts_yyyy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyyz[i] = 10.0 * ts_yyyy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyzz[i] = 10.0 * ts_yyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xzzz[i] = 10.0 * ts_yyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyyy[i] = 10.0 * ts_yyyy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyyz[i] = 10.0 * ts_yyyy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyzz[i] = 10.0 * ts_yyyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yzzz[i] = 10.0 * ts_yyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_zzzz[i] = 10.0 * ts_yyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 555-570 components of targeted buffer : HG

    auto gs_y_yyyyz_xxxx = pbuffer.data(idx_g_hg + 555);

    auto gs_y_yyyyz_xxxy = pbuffer.data(idx_g_hg + 556);

    auto gs_y_yyyyz_xxxz = pbuffer.data(idx_g_hg + 557);

    auto gs_y_yyyyz_xxyy = pbuffer.data(idx_g_hg + 558);

    auto gs_y_yyyyz_xxyz = pbuffer.data(idx_g_hg + 559);

    auto gs_y_yyyyz_xxzz = pbuffer.data(idx_g_hg + 560);

    auto gs_y_yyyyz_xyyy = pbuffer.data(idx_g_hg + 561);

    auto gs_y_yyyyz_xyyz = pbuffer.data(idx_g_hg + 562);

    auto gs_y_yyyyz_xyzz = pbuffer.data(idx_g_hg + 563);

    auto gs_y_yyyyz_xzzz = pbuffer.data(idx_g_hg + 564);

    auto gs_y_yyyyz_yyyy = pbuffer.data(idx_g_hg + 565);

    auto gs_y_yyyyz_yyyz = pbuffer.data(idx_g_hg + 566);

    auto gs_y_yyyyz_yyzz = pbuffer.data(idx_g_hg + 567);

    auto gs_y_yyyyz_yzzz = pbuffer.data(idx_g_hg + 568);

    auto gs_y_yyyyz_zzzz = pbuffer.data(idx_g_hg + 569);

    #pragma omp simd aligned(gc_y, gs_y_yyyyz_xxxx, gs_y_yyyyz_xxxy, gs_y_yyyyz_xxxz, gs_y_yyyyz_xxyy, gs_y_yyyyz_xxyz, gs_y_yyyyz_xxzz, gs_y_yyyyz_xyyy, gs_y_yyyyz_xyyz, gs_y_yyyyz_xyzz, gs_y_yyyyz_xzzz, gs_y_yyyyz_yyyy, gs_y_yyyyz_yyyz, gs_y_yyyyz_yyzz, gs_y_yyyyz_yzzz, gs_y_yyyyz_zzzz, ts_yyyyz_xxx, ts_yyyyz_xxxx, ts_yyyyz_xxxy, ts_yyyyz_xxxz, ts_yyyyz_xxy, ts_yyyyz_xxyy, ts_yyyyz_xxyz, ts_yyyyz_xxz, ts_yyyyz_xxzz, ts_yyyyz_xyy, ts_yyyyz_xyyy, ts_yyyyz_xyyz, ts_yyyyz_xyz, ts_yyyyz_xyzz, ts_yyyyz_xzz, ts_yyyyz_xzzz, ts_yyyyz_yyy, ts_yyyyz_yyyy, ts_yyyyz_yyyz, ts_yyyyz_yyz, ts_yyyyz_yyzz, ts_yyyyz_yzz, ts_yyyyz_yzzz, ts_yyyyz_zzz, ts_yyyyz_zzzz, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxzz, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyzz, ts_yyyz_xzzz, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyzz, ts_yyyz_yzzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyz_xxxx[i] = 8.0 * ts_yyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxy[i] = 8.0 * ts_yyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxz[i] = 8.0 * ts_yyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxyy[i] = 8.0 * ts_yyyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxyz[i] = 8.0 * ts_yyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxzz[i] = 8.0 * ts_yyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyyy[i] = 8.0 * ts_yyyz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyyz[i] = 8.0 * ts_yyyz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyzz[i] = 8.0 * ts_yyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xzzz[i] = 8.0 * ts_yyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyyy[i] = 8.0 * ts_yyyz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyyz[i] = 8.0 * ts_yyyz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyzz[i] = 8.0 * ts_yyyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yzzz[i] = 8.0 * ts_yyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_zzzz[i] = 8.0 * ts_yyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 570-585 components of targeted buffer : HG

    auto gs_y_yyyzz_xxxx = pbuffer.data(idx_g_hg + 570);

    auto gs_y_yyyzz_xxxy = pbuffer.data(idx_g_hg + 571);

    auto gs_y_yyyzz_xxxz = pbuffer.data(idx_g_hg + 572);

    auto gs_y_yyyzz_xxyy = pbuffer.data(idx_g_hg + 573);

    auto gs_y_yyyzz_xxyz = pbuffer.data(idx_g_hg + 574);

    auto gs_y_yyyzz_xxzz = pbuffer.data(idx_g_hg + 575);

    auto gs_y_yyyzz_xyyy = pbuffer.data(idx_g_hg + 576);

    auto gs_y_yyyzz_xyyz = pbuffer.data(idx_g_hg + 577);

    auto gs_y_yyyzz_xyzz = pbuffer.data(idx_g_hg + 578);

    auto gs_y_yyyzz_xzzz = pbuffer.data(idx_g_hg + 579);

    auto gs_y_yyyzz_yyyy = pbuffer.data(idx_g_hg + 580);

    auto gs_y_yyyzz_yyyz = pbuffer.data(idx_g_hg + 581);

    auto gs_y_yyyzz_yyzz = pbuffer.data(idx_g_hg + 582);

    auto gs_y_yyyzz_yzzz = pbuffer.data(idx_g_hg + 583);

    auto gs_y_yyyzz_zzzz = pbuffer.data(idx_g_hg + 584);

    #pragma omp simd aligned(gc_y, gs_y_yyyzz_xxxx, gs_y_yyyzz_xxxy, gs_y_yyyzz_xxxz, gs_y_yyyzz_xxyy, gs_y_yyyzz_xxyz, gs_y_yyyzz_xxzz, gs_y_yyyzz_xyyy, gs_y_yyyzz_xyyz, gs_y_yyyzz_xyzz, gs_y_yyyzz_xzzz, gs_y_yyyzz_yyyy, gs_y_yyyzz_yyyz, gs_y_yyyzz_yyzz, gs_y_yyyzz_yzzz, gs_y_yyyzz_zzzz, ts_yyyzz_xxx, ts_yyyzz_xxxx, ts_yyyzz_xxxy, ts_yyyzz_xxxz, ts_yyyzz_xxy, ts_yyyzz_xxyy, ts_yyyzz_xxyz, ts_yyyzz_xxz, ts_yyyzz_xxzz, ts_yyyzz_xyy, ts_yyyzz_xyyy, ts_yyyzz_xyyz, ts_yyyzz_xyz, ts_yyyzz_xyzz, ts_yyyzz_xzz, ts_yyyzz_xzzz, ts_yyyzz_yyy, ts_yyyzz_yyyy, ts_yyyzz_yyyz, ts_yyyzz_yyz, ts_yyyzz_yyzz, ts_yyyzz_yzz, ts_yyyzz_yzzz, ts_yyyzz_zzz, ts_yyyzz_zzzz, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxzz, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyzz, ts_yyzz_xzzz, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyzz, ts_yyzz_yzzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyzz_xxxx[i] = 6.0 * ts_yyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxy[i] = 6.0 * ts_yyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxz[i] = 6.0 * ts_yyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxyy[i] = 6.0 * ts_yyzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxyz[i] = 6.0 * ts_yyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxzz[i] = 6.0 * ts_yyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyyy[i] = 6.0 * ts_yyzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyyz[i] = 6.0 * ts_yyzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyzz[i] = 6.0 * ts_yyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xzzz[i] = 6.0 * ts_yyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyyy[i] = 6.0 * ts_yyzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyyz[i] = 6.0 * ts_yyzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyzz[i] = 6.0 * ts_yyzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yzzz[i] = 6.0 * ts_yyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_zzzz[i] = 6.0 * ts_yyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 585-600 components of targeted buffer : HG

    auto gs_y_yyzzz_xxxx = pbuffer.data(idx_g_hg + 585);

    auto gs_y_yyzzz_xxxy = pbuffer.data(idx_g_hg + 586);

    auto gs_y_yyzzz_xxxz = pbuffer.data(idx_g_hg + 587);

    auto gs_y_yyzzz_xxyy = pbuffer.data(idx_g_hg + 588);

    auto gs_y_yyzzz_xxyz = pbuffer.data(idx_g_hg + 589);

    auto gs_y_yyzzz_xxzz = pbuffer.data(idx_g_hg + 590);

    auto gs_y_yyzzz_xyyy = pbuffer.data(idx_g_hg + 591);

    auto gs_y_yyzzz_xyyz = pbuffer.data(idx_g_hg + 592);

    auto gs_y_yyzzz_xyzz = pbuffer.data(idx_g_hg + 593);

    auto gs_y_yyzzz_xzzz = pbuffer.data(idx_g_hg + 594);

    auto gs_y_yyzzz_yyyy = pbuffer.data(idx_g_hg + 595);

    auto gs_y_yyzzz_yyyz = pbuffer.data(idx_g_hg + 596);

    auto gs_y_yyzzz_yyzz = pbuffer.data(idx_g_hg + 597);

    auto gs_y_yyzzz_yzzz = pbuffer.data(idx_g_hg + 598);

    auto gs_y_yyzzz_zzzz = pbuffer.data(idx_g_hg + 599);

    #pragma omp simd aligned(gc_y, gs_y_yyzzz_xxxx, gs_y_yyzzz_xxxy, gs_y_yyzzz_xxxz, gs_y_yyzzz_xxyy, gs_y_yyzzz_xxyz, gs_y_yyzzz_xxzz, gs_y_yyzzz_xyyy, gs_y_yyzzz_xyyz, gs_y_yyzzz_xyzz, gs_y_yyzzz_xzzz, gs_y_yyzzz_yyyy, gs_y_yyzzz_yyyz, gs_y_yyzzz_yyzz, gs_y_yyzzz_yzzz, gs_y_yyzzz_zzzz, ts_yyzzz_xxx, ts_yyzzz_xxxx, ts_yyzzz_xxxy, ts_yyzzz_xxxz, ts_yyzzz_xxy, ts_yyzzz_xxyy, ts_yyzzz_xxyz, ts_yyzzz_xxz, ts_yyzzz_xxzz, ts_yyzzz_xyy, ts_yyzzz_xyyy, ts_yyzzz_xyyz, ts_yyzzz_xyz, ts_yyzzz_xyzz, ts_yyzzz_xzz, ts_yyzzz_xzzz, ts_yyzzz_yyy, ts_yyzzz_yyyy, ts_yyzzz_yyyz, ts_yyzzz_yyz, ts_yyzzz_yyzz, ts_yyzzz_yzz, ts_yyzzz_yzzz, ts_yyzzz_zzz, ts_yyzzz_zzzz, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxzz, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyzz, ts_yzzz_xzzz, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyzz, ts_yzzz_yzzz, ts_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzzz_xxxx[i] = 4.0 * ts_yzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxy[i] = 4.0 * ts_yzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxz[i] = 4.0 * ts_yzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxyy[i] = 4.0 * ts_yzzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxyz[i] = 4.0 * ts_yzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxzz[i] = 4.0 * ts_yzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyyy[i] = 4.0 * ts_yzzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyyz[i] = 4.0 * ts_yzzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyzz[i] = 4.0 * ts_yzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xzzz[i] = 4.0 * ts_yzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyyy[i] = 4.0 * ts_yzzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyyz[i] = 4.0 * ts_yzzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyzz[i] = 4.0 * ts_yzzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yzzz[i] = 4.0 * ts_yzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_zzzz[i] = 4.0 * ts_yzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 600-615 components of targeted buffer : HG

    auto gs_y_yzzzz_xxxx = pbuffer.data(idx_g_hg + 600);

    auto gs_y_yzzzz_xxxy = pbuffer.data(idx_g_hg + 601);

    auto gs_y_yzzzz_xxxz = pbuffer.data(idx_g_hg + 602);

    auto gs_y_yzzzz_xxyy = pbuffer.data(idx_g_hg + 603);

    auto gs_y_yzzzz_xxyz = pbuffer.data(idx_g_hg + 604);

    auto gs_y_yzzzz_xxzz = pbuffer.data(idx_g_hg + 605);

    auto gs_y_yzzzz_xyyy = pbuffer.data(idx_g_hg + 606);

    auto gs_y_yzzzz_xyyz = pbuffer.data(idx_g_hg + 607);

    auto gs_y_yzzzz_xyzz = pbuffer.data(idx_g_hg + 608);

    auto gs_y_yzzzz_xzzz = pbuffer.data(idx_g_hg + 609);

    auto gs_y_yzzzz_yyyy = pbuffer.data(idx_g_hg + 610);

    auto gs_y_yzzzz_yyyz = pbuffer.data(idx_g_hg + 611);

    auto gs_y_yzzzz_yyzz = pbuffer.data(idx_g_hg + 612);

    auto gs_y_yzzzz_yzzz = pbuffer.data(idx_g_hg + 613);

    auto gs_y_yzzzz_zzzz = pbuffer.data(idx_g_hg + 614);

    #pragma omp simd aligned(gc_y, gs_y_yzzzz_xxxx, gs_y_yzzzz_xxxy, gs_y_yzzzz_xxxz, gs_y_yzzzz_xxyy, gs_y_yzzzz_xxyz, gs_y_yzzzz_xxzz, gs_y_yzzzz_xyyy, gs_y_yzzzz_xyyz, gs_y_yzzzz_xyzz, gs_y_yzzzz_xzzz, gs_y_yzzzz_yyyy, gs_y_yzzzz_yyyz, gs_y_yzzzz_yyzz, gs_y_yzzzz_yzzz, gs_y_yzzzz_zzzz, ts_yzzzz_xxx, ts_yzzzz_xxxx, ts_yzzzz_xxxy, ts_yzzzz_xxxz, ts_yzzzz_xxy, ts_yzzzz_xxyy, ts_yzzzz_xxyz, ts_yzzzz_xxz, ts_yzzzz_xxzz, ts_yzzzz_xyy, ts_yzzzz_xyyy, ts_yzzzz_xyyz, ts_yzzzz_xyz, ts_yzzzz_xyzz, ts_yzzzz_xzz, ts_yzzzz_xzzz, ts_yzzzz_yyy, ts_yzzzz_yyyy, ts_yzzzz_yyyz, ts_yzzzz_yyz, ts_yzzzz_yyzz, ts_yzzzz_yzz, ts_yzzzz_yzzz, ts_yzzzz_zzz, ts_yzzzz_zzzz, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzzz_xxxx[i] = 2.0 * ts_zzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxy[i] = 2.0 * ts_zzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxz[i] = 2.0 * ts_zzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxyy[i] = 2.0 * ts_zzzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxyz[i] = 2.0 * ts_zzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxzz[i] = 2.0 * ts_zzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyyy[i] = 2.0 * ts_zzzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyyz[i] = 2.0 * ts_zzzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyzz[i] = 2.0 * ts_zzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xzzz[i] = 2.0 * ts_zzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyyy[i] = 2.0 * ts_zzzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyyz[i] = 2.0 * ts_zzzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyzz[i] = 2.0 * ts_zzzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yzzz[i] = 2.0 * ts_zzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_zzzz[i] = 2.0 * ts_zzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 615-630 components of targeted buffer : HG

    auto gs_y_zzzzz_xxxx = pbuffer.data(idx_g_hg + 615);

    auto gs_y_zzzzz_xxxy = pbuffer.data(idx_g_hg + 616);

    auto gs_y_zzzzz_xxxz = pbuffer.data(idx_g_hg + 617);

    auto gs_y_zzzzz_xxyy = pbuffer.data(idx_g_hg + 618);

    auto gs_y_zzzzz_xxyz = pbuffer.data(idx_g_hg + 619);

    auto gs_y_zzzzz_xxzz = pbuffer.data(idx_g_hg + 620);

    auto gs_y_zzzzz_xyyy = pbuffer.data(idx_g_hg + 621);

    auto gs_y_zzzzz_xyyz = pbuffer.data(idx_g_hg + 622);

    auto gs_y_zzzzz_xyzz = pbuffer.data(idx_g_hg + 623);

    auto gs_y_zzzzz_xzzz = pbuffer.data(idx_g_hg + 624);

    auto gs_y_zzzzz_yyyy = pbuffer.data(idx_g_hg + 625);

    auto gs_y_zzzzz_yyyz = pbuffer.data(idx_g_hg + 626);

    auto gs_y_zzzzz_yyzz = pbuffer.data(idx_g_hg + 627);

    auto gs_y_zzzzz_yzzz = pbuffer.data(idx_g_hg + 628);

    auto gs_y_zzzzz_zzzz = pbuffer.data(idx_g_hg + 629);

    #pragma omp simd aligned(gc_y, gs_y_zzzzz_xxxx, gs_y_zzzzz_xxxy, gs_y_zzzzz_xxxz, gs_y_zzzzz_xxyy, gs_y_zzzzz_xxyz, gs_y_zzzzz_xxzz, gs_y_zzzzz_xyyy, gs_y_zzzzz_xyyz, gs_y_zzzzz_xyzz, gs_y_zzzzz_xzzz, gs_y_zzzzz_yyyy, gs_y_zzzzz_yyyz, gs_y_zzzzz_yyzz, gs_y_zzzzz_yzzz, gs_y_zzzzz_zzzz, ts_zzzzz_xxx, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxy, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxz, ts_zzzzz_xxzz, ts_zzzzz_xyy, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyz, ts_zzzzz_xyzz, ts_zzzzz_xzz, ts_zzzzz_xzzz, ts_zzzzz_yyy, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyz, ts_zzzzz_yyzz, ts_zzzzz_yzz, ts_zzzzz_yzzz, ts_zzzzz_zzz, ts_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzzz_xxxx[i] = 2.0 * ts_zzzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxy[i] = 2.0 * ts_zzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxz[i] = 2.0 * ts_zzzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxyy[i] = 4.0 * ts_zzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxyz[i] = 2.0 * ts_zzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxzz[i] = 2.0 * ts_zzzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyyy[i] = 6.0 * ts_zzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyyz[i] = 4.0 * ts_zzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyzz[i] = 2.0 * ts_zzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xzzz[i] = 2.0 * ts_zzzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyyy[i] = 8.0 * ts_zzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyyz[i] = 6.0 * ts_zzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyzz[i] = 4.0 * ts_zzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yzzz[i] = 2.0 * ts_zzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_zzzz[i] = 2.0 * ts_zzzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 630-645 components of targeted buffer : HG

    auto gs_z_xxxxx_xxxx = pbuffer.data(idx_g_hg + 630);

    auto gs_z_xxxxx_xxxy = pbuffer.data(idx_g_hg + 631);

    auto gs_z_xxxxx_xxxz = pbuffer.data(idx_g_hg + 632);

    auto gs_z_xxxxx_xxyy = pbuffer.data(idx_g_hg + 633);

    auto gs_z_xxxxx_xxyz = pbuffer.data(idx_g_hg + 634);

    auto gs_z_xxxxx_xxzz = pbuffer.data(idx_g_hg + 635);

    auto gs_z_xxxxx_xyyy = pbuffer.data(idx_g_hg + 636);

    auto gs_z_xxxxx_xyyz = pbuffer.data(idx_g_hg + 637);

    auto gs_z_xxxxx_xyzz = pbuffer.data(idx_g_hg + 638);

    auto gs_z_xxxxx_xzzz = pbuffer.data(idx_g_hg + 639);

    auto gs_z_xxxxx_yyyy = pbuffer.data(idx_g_hg + 640);

    auto gs_z_xxxxx_yyyz = pbuffer.data(idx_g_hg + 641);

    auto gs_z_xxxxx_yyzz = pbuffer.data(idx_g_hg + 642);

    auto gs_z_xxxxx_yzzz = pbuffer.data(idx_g_hg + 643);

    auto gs_z_xxxxx_zzzz = pbuffer.data(idx_g_hg + 644);

    #pragma omp simd aligned(gc_z, gs_z_xxxxx_xxxx, gs_z_xxxxx_xxxy, gs_z_xxxxx_xxxz, gs_z_xxxxx_xxyy, gs_z_xxxxx_xxyz, gs_z_xxxxx_xxzz, gs_z_xxxxx_xyyy, gs_z_xxxxx_xyyz, gs_z_xxxxx_xyzz, gs_z_xxxxx_xzzz, gs_z_xxxxx_yyyy, gs_z_xxxxx_yyyz, gs_z_xxxxx_yyzz, gs_z_xxxxx_yzzz, gs_z_xxxxx_zzzz, ts_xxxxx_xxx, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxy, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxz, ts_xxxxx_xxzz, ts_xxxxx_xyy, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyz, ts_xxxxx_xyzz, ts_xxxxx_xzz, ts_xxxxx_xzzz, ts_xxxxx_yyy, ts_xxxxx_yyyy, ts_xxxxx_yyyz, ts_xxxxx_yyz, ts_xxxxx_yyzz, ts_xxxxx_yzz, ts_xxxxx_yzzz, ts_xxxxx_zzz, ts_xxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxx_xxxx[i] = 2.0 * ts_xxxxx_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxy[i] = 2.0 * ts_xxxxx_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxz[i] = 2.0 * ts_xxxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxyy[i] = 2.0 * ts_xxxxx_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxyz[i] = 2.0 * ts_xxxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxzz[i] = 4.0 * ts_xxxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyyy[i] = 2.0 * ts_xxxxx_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyyz[i] = 2.0 * ts_xxxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyzz[i] = 4.0 * ts_xxxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xzzz[i] = 6.0 * ts_xxxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyyy[i] = 2.0 * ts_xxxxx_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyyz[i] = 2.0 * ts_xxxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyzz[i] = 4.0 * ts_xxxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yzzz[i] = 6.0 * ts_xxxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_zzzz[i] = 8.0 * ts_xxxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 645-660 components of targeted buffer : HG

    auto gs_z_xxxxy_xxxx = pbuffer.data(idx_g_hg + 645);

    auto gs_z_xxxxy_xxxy = pbuffer.data(idx_g_hg + 646);

    auto gs_z_xxxxy_xxxz = pbuffer.data(idx_g_hg + 647);

    auto gs_z_xxxxy_xxyy = pbuffer.data(idx_g_hg + 648);

    auto gs_z_xxxxy_xxyz = pbuffer.data(idx_g_hg + 649);

    auto gs_z_xxxxy_xxzz = pbuffer.data(idx_g_hg + 650);

    auto gs_z_xxxxy_xyyy = pbuffer.data(idx_g_hg + 651);

    auto gs_z_xxxxy_xyyz = pbuffer.data(idx_g_hg + 652);

    auto gs_z_xxxxy_xyzz = pbuffer.data(idx_g_hg + 653);

    auto gs_z_xxxxy_xzzz = pbuffer.data(idx_g_hg + 654);

    auto gs_z_xxxxy_yyyy = pbuffer.data(idx_g_hg + 655);

    auto gs_z_xxxxy_yyyz = pbuffer.data(idx_g_hg + 656);

    auto gs_z_xxxxy_yyzz = pbuffer.data(idx_g_hg + 657);

    auto gs_z_xxxxy_yzzz = pbuffer.data(idx_g_hg + 658);

    auto gs_z_xxxxy_zzzz = pbuffer.data(idx_g_hg + 659);

    #pragma omp simd aligned(gc_z, gs_z_xxxxy_xxxx, gs_z_xxxxy_xxxy, gs_z_xxxxy_xxxz, gs_z_xxxxy_xxyy, gs_z_xxxxy_xxyz, gs_z_xxxxy_xxzz, gs_z_xxxxy_xyyy, gs_z_xxxxy_xyyz, gs_z_xxxxy_xyzz, gs_z_xxxxy_xzzz, gs_z_xxxxy_yyyy, gs_z_xxxxy_yyyz, gs_z_xxxxy_yyzz, gs_z_xxxxy_yzzz, gs_z_xxxxy_zzzz, ts_xxxxy_xxx, ts_xxxxy_xxxx, ts_xxxxy_xxxy, ts_xxxxy_xxxz, ts_xxxxy_xxy, ts_xxxxy_xxyy, ts_xxxxy_xxyz, ts_xxxxy_xxz, ts_xxxxy_xxzz, ts_xxxxy_xyy, ts_xxxxy_xyyy, ts_xxxxy_xyyz, ts_xxxxy_xyz, ts_xxxxy_xyzz, ts_xxxxy_xzz, ts_xxxxy_xzzz, ts_xxxxy_yyy, ts_xxxxy_yyyy, ts_xxxxy_yyyz, ts_xxxxy_yyz, ts_xxxxy_yyzz, ts_xxxxy_yzz, ts_xxxxy_yzzz, ts_xxxxy_zzz, ts_xxxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxy_xxxx[i] = 2.0 * ts_xxxxy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxy[i] = 2.0 * ts_xxxxy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxz[i] = 2.0 * ts_xxxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxyy[i] = 2.0 * ts_xxxxy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxyz[i] = 2.0 * ts_xxxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxzz[i] = 4.0 * ts_xxxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyyy[i] = 2.0 * ts_xxxxy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyyz[i] = 2.0 * ts_xxxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyzz[i] = 4.0 * ts_xxxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xzzz[i] = 6.0 * ts_xxxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyyy[i] = 2.0 * ts_xxxxy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyyz[i] = 2.0 * ts_xxxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyzz[i] = 4.0 * ts_xxxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yzzz[i] = 6.0 * ts_xxxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_zzzz[i] = 8.0 * ts_xxxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 660-675 components of targeted buffer : HG

    auto gs_z_xxxxz_xxxx = pbuffer.data(idx_g_hg + 660);

    auto gs_z_xxxxz_xxxy = pbuffer.data(idx_g_hg + 661);

    auto gs_z_xxxxz_xxxz = pbuffer.data(idx_g_hg + 662);

    auto gs_z_xxxxz_xxyy = pbuffer.data(idx_g_hg + 663);

    auto gs_z_xxxxz_xxyz = pbuffer.data(idx_g_hg + 664);

    auto gs_z_xxxxz_xxzz = pbuffer.data(idx_g_hg + 665);

    auto gs_z_xxxxz_xyyy = pbuffer.data(idx_g_hg + 666);

    auto gs_z_xxxxz_xyyz = pbuffer.data(idx_g_hg + 667);

    auto gs_z_xxxxz_xyzz = pbuffer.data(idx_g_hg + 668);

    auto gs_z_xxxxz_xzzz = pbuffer.data(idx_g_hg + 669);

    auto gs_z_xxxxz_yyyy = pbuffer.data(idx_g_hg + 670);

    auto gs_z_xxxxz_yyyz = pbuffer.data(idx_g_hg + 671);

    auto gs_z_xxxxz_yyzz = pbuffer.data(idx_g_hg + 672);

    auto gs_z_xxxxz_yzzz = pbuffer.data(idx_g_hg + 673);

    auto gs_z_xxxxz_zzzz = pbuffer.data(idx_g_hg + 674);

    #pragma omp simd aligned(gc_z, gs_z_xxxxz_xxxx, gs_z_xxxxz_xxxy, gs_z_xxxxz_xxxz, gs_z_xxxxz_xxyy, gs_z_xxxxz_xxyz, gs_z_xxxxz_xxzz, gs_z_xxxxz_xyyy, gs_z_xxxxz_xyyz, gs_z_xxxxz_xyzz, gs_z_xxxxz_xzzz, gs_z_xxxxz_yyyy, gs_z_xxxxz_yyyz, gs_z_xxxxz_yyzz, gs_z_xxxxz_yzzz, gs_z_xxxxz_zzzz, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxzz, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyzz, ts_xxxx_xzzz, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyzz, ts_xxxx_yzzz, ts_xxxx_zzzz, ts_xxxxz_xxx, ts_xxxxz_xxxx, ts_xxxxz_xxxy, ts_xxxxz_xxxz, ts_xxxxz_xxy, ts_xxxxz_xxyy, ts_xxxxz_xxyz, ts_xxxxz_xxz, ts_xxxxz_xxzz, ts_xxxxz_xyy, ts_xxxxz_xyyy, ts_xxxxz_xyyz, ts_xxxxz_xyz, ts_xxxxz_xyzz, ts_xxxxz_xzz, ts_xxxxz_xzzz, ts_xxxxz_yyy, ts_xxxxz_yyyy, ts_xxxxz_yyyz, ts_xxxxz_yyz, ts_xxxxz_yyzz, ts_xxxxz_yzz, ts_xxxxz_yzzz, ts_xxxxz_zzz, ts_xxxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxz_xxxx[i] = 2.0 * ts_xxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxy[i] = 2.0 * ts_xxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxz[i] = 2.0 * ts_xxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxyy[i] = 2.0 * ts_xxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxyz[i] = 2.0 * ts_xxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxzz[i] = 2.0 * ts_xxxx_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyyy[i] = 2.0 * ts_xxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyyz[i] = 2.0 * ts_xxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyzz[i] = 2.0 * ts_xxxx_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xzzz[i] = 2.0 * ts_xxxx_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyyy[i] = 2.0 * ts_xxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyyz[i] = 2.0 * ts_xxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyzz[i] = 2.0 * ts_xxxx_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yzzz[i] = 2.0 * ts_xxxx_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_zzzz[i] = 2.0 * ts_xxxx_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 675-690 components of targeted buffer : HG

    auto gs_z_xxxyy_xxxx = pbuffer.data(idx_g_hg + 675);

    auto gs_z_xxxyy_xxxy = pbuffer.data(idx_g_hg + 676);

    auto gs_z_xxxyy_xxxz = pbuffer.data(idx_g_hg + 677);

    auto gs_z_xxxyy_xxyy = pbuffer.data(idx_g_hg + 678);

    auto gs_z_xxxyy_xxyz = pbuffer.data(idx_g_hg + 679);

    auto gs_z_xxxyy_xxzz = pbuffer.data(idx_g_hg + 680);

    auto gs_z_xxxyy_xyyy = pbuffer.data(idx_g_hg + 681);

    auto gs_z_xxxyy_xyyz = pbuffer.data(idx_g_hg + 682);

    auto gs_z_xxxyy_xyzz = pbuffer.data(idx_g_hg + 683);

    auto gs_z_xxxyy_xzzz = pbuffer.data(idx_g_hg + 684);

    auto gs_z_xxxyy_yyyy = pbuffer.data(idx_g_hg + 685);

    auto gs_z_xxxyy_yyyz = pbuffer.data(idx_g_hg + 686);

    auto gs_z_xxxyy_yyzz = pbuffer.data(idx_g_hg + 687);

    auto gs_z_xxxyy_yzzz = pbuffer.data(idx_g_hg + 688);

    auto gs_z_xxxyy_zzzz = pbuffer.data(idx_g_hg + 689);

    #pragma omp simd aligned(gc_z, gs_z_xxxyy_xxxx, gs_z_xxxyy_xxxy, gs_z_xxxyy_xxxz, gs_z_xxxyy_xxyy, gs_z_xxxyy_xxyz, gs_z_xxxyy_xxzz, gs_z_xxxyy_xyyy, gs_z_xxxyy_xyyz, gs_z_xxxyy_xyzz, gs_z_xxxyy_xzzz, gs_z_xxxyy_yyyy, gs_z_xxxyy_yyyz, gs_z_xxxyy_yyzz, gs_z_xxxyy_yzzz, gs_z_xxxyy_zzzz, ts_xxxyy_xxx, ts_xxxyy_xxxx, ts_xxxyy_xxxy, ts_xxxyy_xxxz, ts_xxxyy_xxy, ts_xxxyy_xxyy, ts_xxxyy_xxyz, ts_xxxyy_xxz, ts_xxxyy_xxzz, ts_xxxyy_xyy, ts_xxxyy_xyyy, ts_xxxyy_xyyz, ts_xxxyy_xyz, ts_xxxyy_xyzz, ts_xxxyy_xzz, ts_xxxyy_xzzz, ts_xxxyy_yyy, ts_xxxyy_yyyy, ts_xxxyy_yyyz, ts_xxxyy_yyz, ts_xxxyy_yyzz, ts_xxxyy_yzz, ts_xxxyy_yzzz, ts_xxxyy_zzz, ts_xxxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyy_xxxx[i] = 2.0 * ts_xxxyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxy[i] = 2.0 * ts_xxxyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxz[i] = 2.0 * ts_xxxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxyy[i] = 2.0 * ts_xxxyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxyz[i] = 2.0 * ts_xxxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxzz[i] = 4.0 * ts_xxxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyyy[i] = 2.0 * ts_xxxyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyyz[i] = 2.0 * ts_xxxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyzz[i] = 4.0 * ts_xxxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xzzz[i] = 6.0 * ts_xxxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyyy[i] = 2.0 * ts_xxxyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyyz[i] = 2.0 * ts_xxxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyzz[i] = 4.0 * ts_xxxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yzzz[i] = 6.0 * ts_xxxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_zzzz[i] = 8.0 * ts_xxxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 690-705 components of targeted buffer : HG

    auto gs_z_xxxyz_xxxx = pbuffer.data(idx_g_hg + 690);

    auto gs_z_xxxyz_xxxy = pbuffer.data(idx_g_hg + 691);

    auto gs_z_xxxyz_xxxz = pbuffer.data(idx_g_hg + 692);

    auto gs_z_xxxyz_xxyy = pbuffer.data(idx_g_hg + 693);

    auto gs_z_xxxyz_xxyz = pbuffer.data(idx_g_hg + 694);

    auto gs_z_xxxyz_xxzz = pbuffer.data(idx_g_hg + 695);

    auto gs_z_xxxyz_xyyy = pbuffer.data(idx_g_hg + 696);

    auto gs_z_xxxyz_xyyz = pbuffer.data(idx_g_hg + 697);

    auto gs_z_xxxyz_xyzz = pbuffer.data(idx_g_hg + 698);

    auto gs_z_xxxyz_xzzz = pbuffer.data(idx_g_hg + 699);

    auto gs_z_xxxyz_yyyy = pbuffer.data(idx_g_hg + 700);

    auto gs_z_xxxyz_yyyz = pbuffer.data(idx_g_hg + 701);

    auto gs_z_xxxyz_yyzz = pbuffer.data(idx_g_hg + 702);

    auto gs_z_xxxyz_yzzz = pbuffer.data(idx_g_hg + 703);

    auto gs_z_xxxyz_zzzz = pbuffer.data(idx_g_hg + 704);

    #pragma omp simd aligned(gc_z, gs_z_xxxyz_xxxx, gs_z_xxxyz_xxxy, gs_z_xxxyz_xxxz, gs_z_xxxyz_xxyy, gs_z_xxxyz_xxyz, gs_z_xxxyz_xxzz, gs_z_xxxyz_xyyy, gs_z_xxxyz_xyyz, gs_z_xxxyz_xyzz, gs_z_xxxyz_xzzz, gs_z_xxxyz_yyyy, gs_z_xxxyz_yyyz, gs_z_xxxyz_yyzz, gs_z_xxxyz_yzzz, gs_z_xxxyz_zzzz, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxzz, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyzz, ts_xxxy_xzzz, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyzz, ts_xxxy_yzzz, ts_xxxy_zzzz, ts_xxxyz_xxx, ts_xxxyz_xxxx, ts_xxxyz_xxxy, ts_xxxyz_xxxz, ts_xxxyz_xxy, ts_xxxyz_xxyy, ts_xxxyz_xxyz, ts_xxxyz_xxz, ts_xxxyz_xxzz, ts_xxxyz_xyy, ts_xxxyz_xyyy, ts_xxxyz_xyyz, ts_xxxyz_xyz, ts_xxxyz_xyzz, ts_xxxyz_xzz, ts_xxxyz_xzzz, ts_xxxyz_yyy, ts_xxxyz_yyyy, ts_xxxyz_yyyz, ts_xxxyz_yyz, ts_xxxyz_yyzz, ts_xxxyz_yzz, ts_xxxyz_yzzz, ts_xxxyz_zzz, ts_xxxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyz_xxxx[i] = 2.0 * ts_xxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxy[i] = 2.0 * ts_xxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxz[i] = 2.0 * ts_xxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxyy[i] = 2.0 * ts_xxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxyz[i] = 2.0 * ts_xxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxzz[i] = 2.0 * ts_xxxy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyyy[i] = 2.0 * ts_xxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyyz[i] = 2.0 * ts_xxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyzz[i] = 2.0 * ts_xxxy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xzzz[i] = 2.0 * ts_xxxy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyyy[i] = 2.0 * ts_xxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyyz[i] = 2.0 * ts_xxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyzz[i] = 2.0 * ts_xxxy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yzzz[i] = 2.0 * ts_xxxy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_zzzz[i] = 2.0 * ts_xxxy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 705-720 components of targeted buffer : HG

    auto gs_z_xxxzz_xxxx = pbuffer.data(idx_g_hg + 705);

    auto gs_z_xxxzz_xxxy = pbuffer.data(idx_g_hg + 706);

    auto gs_z_xxxzz_xxxz = pbuffer.data(idx_g_hg + 707);

    auto gs_z_xxxzz_xxyy = pbuffer.data(idx_g_hg + 708);

    auto gs_z_xxxzz_xxyz = pbuffer.data(idx_g_hg + 709);

    auto gs_z_xxxzz_xxzz = pbuffer.data(idx_g_hg + 710);

    auto gs_z_xxxzz_xyyy = pbuffer.data(idx_g_hg + 711);

    auto gs_z_xxxzz_xyyz = pbuffer.data(idx_g_hg + 712);

    auto gs_z_xxxzz_xyzz = pbuffer.data(idx_g_hg + 713);

    auto gs_z_xxxzz_xzzz = pbuffer.data(idx_g_hg + 714);

    auto gs_z_xxxzz_yyyy = pbuffer.data(idx_g_hg + 715);

    auto gs_z_xxxzz_yyyz = pbuffer.data(idx_g_hg + 716);

    auto gs_z_xxxzz_yyzz = pbuffer.data(idx_g_hg + 717);

    auto gs_z_xxxzz_yzzz = pbuffer.data(idx_g_hg + 718);

    auto gs_z_xxxzz_zzzz = pbuffer.data(idx_g_hg + 719);

    #pragma omp simd aligned(gc_z, gs_z_xxxzz_xxxx, gs_z_xxxzz_xxxy, gs_z_xxxzz_xxxz, gs_z_xxxzz_xxyy, gs_z_xxxzz_xxyz, gs_z_xxxzz_xxzz, gs_z_xxxzz_xyyy, gs_z_xxxzz_xyyz, gs_z_xxxzz_xyzz, gs_z_xxxzz_xzzz, gs_z_xxxzz_yyyy, gs_z_xxxzz_yyyz, gs_z_xxxzz_yyzz, gs_z_xxxzz_yzzz, gs_z_xxxzz_zzzz, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxzz, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyzz, ts_xxxz_xzzz, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyzz, ts_xxxz_yzzz, ts_xxxz_zzzz, ts_xxxzz_xxx, ts_xxxzz_xxxx, ts_xxxzz_xxxy, ts_xxxzz_xxxz, ts_xxxzz_xxy, ts_xxxzz_xxyy, ts_xxxzz_xxyz, ts_xxxzz_xxz, ts_xxxzz_xxzz, ts_xxxzz_xyy, ts_xxxzz_xyyy, ts_xxxzz_xyyz, ts_xxxzz_xyz, ts_xxxzz_xyzz, ts_xxxzz_xzz, ts_xxxzz_xzzz, ts_xxxzz_yyy, ts_xxxzz_yyyy, ts_xxxzz_yyyz, ts_xxxzz_yyz, ts_xxxzz_yyzz, ts_xxxzz_yzz, ts_xxxzz_yzzz, ts_xxxzz_zzz, ts_xxxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxzz_xxxx[i] = 4.0 * ts_xxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxy[i] = 4.0 * ts_xxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxz[i] = 4.0 * ts_xxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxyy[i] = 4.0 * ts_xxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxyz[i] = 4.0 * ts_xxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxzz[i] = 4.0 * ts_xxxz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyyy[i] = 4.0 * ts_xxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyyz[i] = 4.0 * ts_xxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyzz[i] = 4.0 * ts_xxxz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xzzz[i] = 4.0 * ts_xxxz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyyy[i] = 4.0 * ts_xxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyyz[i] = 4.0 * ts_xxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyzz[i] = 4.0 * ts_xxxz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yzzz[i] = 4.0 * ts_xxxz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_zzzz[i] = 4.0 * ts_xxxz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 720-735 components of targeted buffer : HG

    auto gs_z_xxyyy_xxxx = pbuffer.data(idx_g_hg + 720);

    auto gs_z_xxyyy_xxxy = pbuffer.data(idx_g_hg + 721);

    auto gs_z_xxyyy_xxxz = pbuffer.data(idx_g_hg + 722);

    auto gs_z_xxyyy_xxyy = pbuffer.data(idx_g_hg + 723);

    auto gs_z_xxyyy_xxyz = pbuffer.data(idx_g_hg + 724);

    auto gs_z_xxyyy_xxzz = pbuffer.data(idx_g_hg + 725);

    auto gs_z_xxyyy_xyyy = pbuffer.data(idx_g_hg + 726);

    auto gs_z_xxyyy_xyyz = pbuffer.data(idx_g_hg + 727);

    auto gs_z_xxyyy_xyzz = pbuffer.data(idx_g_hg + 728);

    auto gs_z_xxyyy_xzzz = pbuffer.data(idx_g_hg + 729);

    auto gs_z_xxyyy_yyyy = pbuffer.data(idx_g_hg + 730);

    auto gs_z_xxyyy_yyyz = pbuffer.data(idx_g_hg + 731);

    auto gs_z_xxyyy_yyzz = pbuffer.data(idx_g_hg + 732);

    auto gs_z_xxyyy_yzzz = pbuffer.data(idx_g_hg + 733);

    auto gs_z_xxyyy_zzzz = pbuffer.data(idx_g_hg + 734);

    #pragma omp simd aligned(gc_z, gs_z_xxyyy_xxxx, gs_z_xxyyy_xxxy, gs_z_xxyyy_xxxz, gs_z_xxyyy_xxyy, gs_z_xxyyy_xxyz, gs_z_xxyyy_xxzz, gs_z_xxyyy_xyyy, gs_z_xxyyy_xyyz, gs_z_xxyyy_xyzz, gs_z_xxyyy_xzzz, gs_z_xxyyy_yyyy, gs_z_xxyyy_yyyz, gs_z_xxyyy_yyzz, gs_z_xxyyy_yzzz, gs_z_xxyyy_zzzz, ts_xxyyy_xxx, ts_xxyyy_xxxx, ts_xxyyy_xxxy, ts_xxyyy_xxxz, ts_xxyyy_xxy, ts_xxyyy_xxyy, ts_xxyyy_xxyz, ts_xxyyy_xxz, ts_xxyyy_xxzz, ts_xxyyy_xyy, ts_xxyyy_xyyy, ts_xxyyy_xyyz, ts_xxyyy_xyz, ts_xxyyy_xyzz, ts_xxyyy_xzz, ts_xxyyy_xzzz, ts_xxyyy_yyy, ts_xxyyy_yyyy, ts_xxyyy_yyyz, ts_xxyyy_yyz, ts_xxyyy_yyzz, ts_xxyyy_yzz, ts_xxyyy_yzzz, ts_xxyyy_zzz, ts_xxyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyy_xxxx[i] = 2.0 * ts_xxyyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxy[i] = 2.0 * ts_xxyyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxz[i] = 2.0 * ts_xxyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxyy[i] = 2.0 * ts_xxyyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxyz[i] = 2.0 * ts_xxyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxzz[i] = 4.0 * ts_xxyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyyy[i] = 2.0 * ts_xxyyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyyz[i] = 2.0 * ts_xxyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyzz[i] = 4.0 * ts_xxyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xzzz[i] = 6.0 * ts_xxyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyyy[i] = 2.0 * ts_xxyyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyyz[i] = 2.0 * ts_xxyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyzz[i] = 4.0 * ts_xxyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yzzz[i] = 6.0 * ts_xxyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_zzzz[i] = 8.0 * ts_xxyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 735-750 components of targeted buffer : HG

    auto gs_z_xxyyz_xxxx = pbuffer.data(idx_g_hg + 735);

    auto gs_z_xxyyz_xxxy = pbuffer.data(idx_g_hg + 736);

    auto gs_z_xxyyz_xxxz = pbuffer.data(idx_g_hg + 737);

    auto gs_z_xxyyz_xxyy = pbuffer.data(idx_g_hg + 738);

    auto gs_z_xxyyz_xxyz = pbuffer.data(idx_g_hg + 739);

    auto gs_z_xxyyz_xxzz = pbuffer.data(idx_g_hg + 740);

    auto gs_z_xxyyz_xyyy = pbuffer.data(idx_g_hg + 741);

    auto gs_z_xxyyz_xyyz = pbuffer.data(idx_g_hg + 742);

    auto gs_z_xxyyz_xyzz = pbuffer.data(idx_g_hg + 743);

    auto gs_z_xxyyz_xzzz = pbuffer.data(idx_g_hg + 744);

    auto gs_z_xxyyz_yyyy = pbuffer.data(idx_g_hg + 745);

    auto gs_z_xxyyz_yyyz = pbuffer.data(idx_g_hg + 746);

    auto gs_z_xxyyz_yyzz = pbuffer.data(idx_g_hg + 747);

    auto gs_z_xxyyz_yzzz = pbuffer.data(idx_g_hg + 748);

    auto gs_z_xxyyz_zzzz = pbuffer.data(idx_g_hg + 749);

    #pragma omp simd aligned(gc_z, gs_z_xxyyz_xxxx, gs_z_xxyyz_xxxy, gs_z_xxyyz_xxxz, gs_z_xxyyz_xxyy, gs_z_xxyyz_xxyz, gs_z_xxyyz_xxzz, gs_z_xxyyz_xyyy, gs_z_xxyyz_xyyz, gs_z_xxyyz_xyzz, gs_z_xxyyz_xzzz, gs_z_xxyyz_yyyy, gs_z_xxyyz_yyyz, gs_z_xxyyz_yyzz, gs_z_xxyyz_yzzz, gs_z_xxyyz_zzzz, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxzz, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyzz, ts_xxyy_xzzz, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyzz, ts_xxyy_yzzz, ts_xxyy_zzzz, ts_xxyyz_xxx, ts_xxyyz_xxxx, ts_xxyyz_xxxy, ts_xxyyz_xxxz, ts_xxyyz_xxy, ts_xxyyz_xxyy, ts_xxyyz_xxyz, ts_xxyyz_xxz, ts_xxyyz_xxzz, ts_xxyyz_xyy, ts_xxyyz_xyyy, ts_xxyyz_xyyz, ts_xxyyz_xyz, ts_xxyyz_xyzz, ts_xxyyz_xzz, ts_xxyyz_xzzz, ts_xxyyz_yyy, ts_xxyyz_yyyy, ts_xxyyz_yyyz, ts_xxyyz_yyz, ts_xxyyz_yyzz, ts_xxyyz_yzz, ts_xxyyz_yzzz, ts_xxyyz_zzz, ts_xxyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyz_xxxx[i] = 2.0 * ts_xxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxy[i] = 2.0 * ts_xxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxz[i] = 2.0 * ts_xxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxyy[i] = 2.0 * ts_xxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxyz[i] = 2.0 * ts_xxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxzz[i] = 2.0 * ts_xxyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyyy[i] = 2.0 * ts_xxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyyz[i] = 2.0 * ts_xxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyzz[i] = 2.0 * ts_xxyy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xzzz[i] = 2.0 * ts_xxyy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyyy[i] = 2.0 * ts_xxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyyz[i] = 2.0 * ts_xxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyzz[i] = 2.0 * ts_xxyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yzzz[i] = 2.0 * ts_xxyy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_zzzz[i] = 2.0 * ts_xxyy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 750-765 components of targeted buffer : HG

    auto gs_z_xxyzz_xxxx = pbuffer.data(idx_g_hg + 750);

    auto gs_z_xxyzz_xxxy = pbuffer.data(idx_g_hg + 751);

    auto gs_z_xxyzz_xxxz = pbuffer.data(idx_g_hg + 752);

    auto gs_z_xxyzz_xxyy = pbuffer.data(idx_g_hg + 753);

    auto gs_z_xxyzz_xxyz = pbuffer.data(idx_g_hg + 754);

    auto gs_z_xxyzz_xxzz = pbuffer.data(idx_g_hg + 755);

    auto gs_z_xxyzz_xyyy = pbuffer.data(idx_g_hg + 756);

    auto gs_z_xxyzz_xyyz = pbuffer.data(idx_g_hg + 757);

    auto gs_z_xxyzz_xyzz = pbuffer.data(idx_g_hg + 758);

    auto gs_z_xxyzz_xzzz = pbuffer.data(idx_g_hg + 759);

    auto gs_z_xxyzz_yyyy = pbuffer.data(idx_g_hg + 760);

    auto gs_z_xxyzz_yyyz = pbuffer.data(idx_g_hg + 761);

    auto gs_z_xxyzz_yyzz = pbuffer.data(idx_g_hg + 762);

    auto gs_z_xxyzz_yzzz = pbuffer.data(idx_g_hg + 763);

    auto gs_z_xxyzz_zzzz = pbuffer.data(idx_g_hg + 764);

    #pragma omp simd aligned(gc_z, gs_z_xxyzz_xxxx, gs_z_xxyzz_xxxy, gs_z_xxyzz_xxxz, gs_z_xxyzz_xxyy, gs_z_xxyzz_xxyz, gs_z_xxyzz_xxzz, gs_z_xxyzz_xyyy, gs_z_xxyzz_xyyz, gs_z_xxyzz_xyzz, gs_z_xxyzz_xzzz, gs_z_xxyzz_yyyy, gs_z_xxyzz_yyyz, gs_z_xxyzz_yyzz, gs_z_xxyzz_yzzz, gs_z_xxyzz_zzzz, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxzz, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyzz, ts_xxyz_xzzz, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyzz, ts_xxyz_yzzz, ts_xxyz_zzzz, ts_xxyzz_xxx, ts_xxyzz_xxxx, ts_xxyzz_xxxy, ts_xxyzz_xxxz, ts_xxyzz_xxy, ts_xxyzz_xxyy, ts_xxyzz_xxyz, ts_xxyzz_xxz, ts_xxyzz_xxzz, ts_xxyzz_xyy, ts_xxyzz_xyyy, ts_xxyzz_xyyz, ts_xxyzz_xyz, ts_xxyzz_xyzz, ts_xxyzz_xzz, ts_xxyzz_xzzz, ts_xxyzz_yyy, ts_xxyzz_yyyy, ts_xxyzz_yyyz, ts_xxyzz_yyz, ts_xxyzz_yyzz, ts_xxyzz_yzz, ts_xxyzz_yzzz, ts_xxyzz_zzz, ts_xxyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyzz_xxxx[i] = 4.0 * ts_xxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxy[i] = 4.0 * ts_xxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxz[i] = 4.0 * ts_xxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxyy[i] = 4.0 * ts_xxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxyz[i] = 4.0 * ts_xxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxzz[i] = 4.0 * ts_xxyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyyy[i] = 4.0 * ts_xxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyyz[i] = 4.0 * ts_xxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyzz[i] = 4.0 * ts_xxyz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xzzz[i] = 4.0 * ts_xxyz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyyy[i] = 4.0 * ts_xxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyyz[i] = 4.0 * ts_xxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyzz[i] = 4.0 * ts_xxyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yzzz[i] = 4.0 * ts_xxyz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_zzzz[i] = 4.0 * ts_xxyz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 765-780 components of targeted buffer : HG

    auto gs_z_xxzzz_xxxx = pbuffer.data(idx_g_hg + 765);

    auto gs_z_xxzzz_xxxy = pbuffer.data(idx_g_hg + 766);

    auto gs_z_xxzzz_xxxz = pbuffer.data(idx_g_hg + 767);

    auto gs_z_xxzzz_xxyy = pbuffer.data(idx_g_hg + 768);

    auto gs_z_xxzzz_xxyz = pbuffer.data(idx_g_hg + 769);

    auto gs_z_xxzzz_xxzz = pbuffer.data(idx_g_hg + 770);

    auto gs_z_xxzzz_xyyy = pbuffer.data(idx_g_hg + 771);

    auto gs_z_xxzzz_xyyz = pbuffer.data(idx_g_hg + 772);

    auto gs_z_xxzzz_xyzz = pbuffer.data(idx_g_hg + 773);

    auto gs_z_xxzzz_xzzz = pbuffer.data(idx_g_hg + 774);

    auto gs_z_xxzzz_yyyy = pbuffer.data(idx_g_hg + 775);

    auto gs_z_xxzzz_yyyz = pbuffer.data(idx_g_hg + 776);

    auto gs_z_xxzzz_yyzz = pbuffer.data(idx_g_hg + 777);

    auto gs_z_xxzzz_yzzz = pbuffer.data(idx_g_hg + 778);

    auto gs_z_xxzzz_zzzz = pbuffer.data(idx_g_hg + 779);

    #pragma omp simd aligned(gc_z, gs_z_xxzzz_xxxx, gs_z_xxzzz_xxxy, gs_z_xxzzz_xxxz, gs_z_xxzzz_xxyy, gs_z_xxzzz_xxyz, gs_z_xxzzz_xxzz, gs_z_xxzzz_xyyy, gs_z_xxzzz_xyyz, gs_z_xxzzz_xyzz, gs_z_xxzzz_xzzz, gs_z_xxzzz_yyyy, gs_z_xxzzz_yyyz, gs_z_xxzzz_yyzz, gs_z_xxzzz_yzzz, gs_z_xxzzz_zzzz, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxzz, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyzz, ts_xxzz_xzzz, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyzz, ts_xxzz_yzzz, ts_xxzz_zzzz, ts_xxzzz_xxx, ts_xxzzz_xxxx, ts_xxzzz_xxxy, ts_xxzzz_xxxz, ts_xxzzz_xxy, ts_xxzzz_xxyy, ts_xxzzz_xxyz, ts_xxzzz_xxz, ts_xxzzz_xxzz, ts_xxzzz_xyy, ts_xxzzz_xyyy, ts_xxzzz_xyyz, ts_xxzzz_xyz, ts_xxzzz_xyzz, ts_xxzzz_xzz, ts_xxzzz_xzzz, ts_xxzzz_yyy, ts_xxzzz_yyyy, ts_xxzzz_yyyz, ts_xxzzz_yyz, ts_xxzzz_yyzz, ts_xxzzz_yzz, ts_xxzzz_yzzz, ts_xxzzz_zzz, ts_xxzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzzz_xxxx[i] = 6.0 * ts_xxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxy[i] = 6.0 * ts_xxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxz[i] = 6.0 * ts_xxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxyy[i] = 6.0 * ts_xxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxyz[i] = 6.0 * ts_xxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxzz[i] = 6.0 * ts_xxzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyyy[i] = 6.0 * ts_xxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyyz[i] = 6.0 * ts_xxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyzz[i] = 6.0 * ts_xxzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xzzz[i] = 6.0 * ts_xxzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyyy[i] = 6.0 * ts_xxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyyz[i] = 6.0 * ts_xxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyzz[i] = 6.0 * ts_xxzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yzzz[i] = 6.0 * ts_xxzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_zzzz[i] = 6.0 * ts_xxzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 780-795 components of targeted buffer : HG

    auto gs_z_xyyyy_xxxx = pbuffer.data(idx_g_hg + 780);

    auto gs_z_xyyyy_xxxy = pbuffer.data(idx_g_hg + 781);

    auto gs_z_xyyyy_xxxz = pbuffer.data(idx_g_hg + 782);

    auto gs_z_xyyyy_xxyy = pbuffer.data(idx_g_hg + 783);

    auto gs_z_xyyyy_xxyz = pbuffer.data(idx_g_hg + 784);

    auto gs_z_xyyyy_xxzz = pbuffer.data(idx_g_hg + 785);

    auto gs_z_xyyyy_xyyy = pbuffer.data(idx_g_hg + 786);

    auto gs_z_xyyyy_xyyz = pbuffer.data(idx_g_hg + 787);

    auto gs_z_xyyyy_xyzz = pbuffer.data(idx_g_hg + 788);

    auto gs_z_xyyyy_xzzz = pbuffer.data(idx_g_hg + 789);

    auto gs_z_xyyyy_yyyy = pbuffer.data(idx_g_hg + 790);

    auto gs_z_xyyyy_yyyz = pbuffer.data(idx_g_hg + 791);

    auto gs_z_xyyyy_yyzz = pbuffer.data(idx_g_hg + 792);

    auto gs_z_xyyyy_yzzz = pbuffer.data(idx_g_hg + 793);

    auto gs_z_xyyyy_zzzz = pbuffer.data(idx_g_hg + 794);

    #pragma omp simd aligned(gc_z, gs_z_xyyyy_xxxx, gs_z_xyyyy_xxxy, gs_z_xyyyy_xxxz, gs_z_xyyyy_xxyy, gs_z_xyyyy_xxyz, gs_z_xyyyy_xxzz, gs_z_xyyyy_xyyy, gs_z_xyyyy_xyyz, gs_z_xyyyy_xyzz, gs_z_xyyyy_xzzz, gs_z_xyyyy_yyyy, gs_z_xyyyy_yyyz, gs_z_xyyyy_yyzz, gs_z_xyyyy_yzzz, gs_z_xyyyy_zzzz, ts_xyyyy_xxx, ts_xyyyy_xxxx, ts_xyyyy_xxxy, ts_xyyyy_xxxz, ts_xyyyy_xxy, ts_xyyyy_xxyy, ts_xyyyy_xxyz, ts_xyyyy_xxz, ts_xyyyy_xxzz, ts_xyyyy_xyy, ts_xyyyy_xyyy, ts_xyyyy_xyyz, ts_xyyyy_xyz, ts_xyyyy_xyzz, ts_xyyyy_xzz, ts_xyyyy_xzzz, ts_xyyyy_yyy, ts_xyyyy_yyyy, ts_xyyyy_yyyz, ts_xyyyy_yyz, ts_xyyyy_yyzz, ts_xyyyy_yzz, ts_xyyyy_yzzz, ts_xyyyy_zzz, ts_xyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyy_xxxx[i] = 2.0 * ts_xyyyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxy[i] = 2.0 * ts_xyyyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxz[i] = 2.0 * ts_xyyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxyy[i] = 2.0 * ts_xyyyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxyz[i] = 2.0 * ts_xyyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxzz[i] = 4.0 * ts_xyyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyyy[i] = 2.0 * ts_xyyyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyyz[i] = 2.0 * ts_xyyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyzz[i] = 4.0 * ts_xyyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xzzz[i] = 6.0 * ts_xyyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyyy[i] = 2.0 * ts_xyyyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyyz[i] = 2.0 * ts_xyyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyzz[i] = 4.0 * ts_xyyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yzzz[i] = 6.0 * ts_xyyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_zzzz[i] = 8.0 * ts_xyyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 795-810 components of targeted buffer : HG

    auto gs_z_xyyyz_xxxx = pbuffer.data(idx_g_hg + 795);

    auto gs_z_xyyyz_xxxy = pbuffer.data(idx_g_hg + 796);

    auto gs_z_xyyyz_xxxz = pbuffer.data(idx_g_hg + 797);

    auto gs_z_xyyyz_xxyy = pbuffer.data(idx_g_hg + 798);

    auto gs_z_xyyyz_xxyz = pbuffer.data(idx_g_hg + 799);

    auto gs_z_xyyyz_xxzz = pbuffer.data(idx_g_hg + 800);

    auto gs_z_xyyyz_xyyy = pbuffer.data(idx_g_hg + 801);

    auto gs_z_xyyyz_xyyz = pbuffer.data(idx_g_hg + 802);

    auto gs_z_xyyyz_xyzz = pbuffer.data(idx_g_hg + 803);

    auto gs_z_xyyyz_xzzz = pbuffer.data(idx_g_hg + 804);

    auto gs_z_xyyyz_yyyy = pbuffer.data(idx_g_hg + 805);

    auto gs_z_xyyyz_yyyz = pbuffer.data(idx_g_hg + 806);

    auto gs_z_xyyyz_yyzz = pbuffer.data(idx_g_hg + 807);

    auto gs_z_xyyyz_yzzz = pbuffer.data(idx_g_hg + 808);

    auto gs_z_xyyyz_zzzz = pbuffer.data(idx_g_hg + 809);

    #pragma omp simd aligned(gc_z, gs_z_xyyyz_xxxx, gs_z_xyyyz_xxxy, gs_z_xyyyz_xxxz, gs_z_xyyyz_xxyy, gs_z_xyyyz_xxyz, gs_z_xyyyz_xxzz, gs_z_xyyyz_xyyy, gs_z_xyyyz_xyyz, gs_z_xyyyz_xyzz, gs_z_xyyyz_xzzz, gs_z_xyyyz_yyyy, gs_z_xyyyz_yyyz, gs_z_xyyyz_yyzz, gs_z_xyyyz_yzzz, gs_z_xyyyz_zzzz, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxzz, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyzz, ts_xyyy_xzzz, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyzz, ts_xyyy_yzzz, ts_xyyy_zzzz, ts_xyyyz_xxx, ts_xyyyz_xxxx, ts_xyyyz_xxxy, ts_xyyyz_xxxz, ts_xyyyz_xxy, ts_xyyyz_xxyy, ts_xyyyz_xxyz, ts_xyyyz_xxz, ts_xyyyz_xxzz, ts_xyyyz_xyy, ts_xyyyz_xyyy, ts_xyyyz_xyyz, ts_xyyyz_xyz, ts_xyyyz_xyzz, ts_xyyyz_xzz, ts_xyyyz_xzzz, ts_xyyyz_yyy, ts_xyyyz_yyyy, ts_xyyyz_yyyz, ts_xyyyz_yyz, ts_xyyyz_yyzz, ts_xyyyz_yzz, ts_xyyyz_yzzz, ts_xyyyz_zzz, ts_xyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyz_xxxx[i] = 2.0 * ts_xyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxy[i] = 2.0 * ts_xyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxz[i] = 2.0 * ts_xyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxyy[i] = 2.0 * ts_xyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxyz[i] = 2.0 * ts_xyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxzz[i] = 2.0 * ts_xyyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyyy[i] = 2.0 * ts_xyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyyz[i] = 2.0 * ts_xyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyzz[i] = 2.0 * ts_xyyy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xzzz[i] = 2.0 * ts_xyyy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyyy[i] = 2.0 * ts_xyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyyz[i] = 2.0 * ts_xyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyzz[i] = 2.0 * ts_xyyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yzzz[i] = 2.0 * ts_xyyy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_zzzz[i] = 2.0 * ts_xyyy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 810-825 components of targeted buffer : HG

    auto gs_z_xyyzz_xxxx = pbuffer.data(idx_g_hg + 810);

    auto gs_z_xyyzz_xxxy = pbuffer.data(idx_g_hg + 811);

    auto gs_z_xyyzz_xxxz = pbuffer.data(idx_g_hg + 812);

    auto gs_z_xyyzz_xxyy = pbuffer.data(idx_g_hg + 813);

    auto gs_z_xyyzz_xxyz = pbuffer.data(idx_g_hg + 814);

    auto gs_z_xyyzz_xxzz = pbuffer.data(idx_g_hg + 815);

    auto gs_z_xyyzz_xyyy = pbuffer.data(idx_g_hg + 816);

    auto gs_z_xyyzz_xyyz = pbuffer.data(idx_g_hg + 817);

    auto gs_z_xyyzz_xyzz = pbuffer.data(idx_g_hg + 818);

    auto gs_z_xyyzz_xzzz = pbuffer.data(idx_g_hg + 819);

    auto gs_z_xyyzz_yyyy = pbuffer.data(idx_g_hg + 820);

    auto gs_z_xyyzz_yyyz = pbuffer.data(idx_g_hg + 821);

    auto gs_z_xyyzz_yyzz = pbuffer.data(idx_g_hg + 822);

    auto gs_z_xyyzz_yzzz = pbuffer.data(idx_g_hg + 823);

    auto gs_z_xyyzz_zzzz = pbuffer.data(idx_g_hg + 824);

    #pragma omp simd aligned(gc_z, gs_z_xyyzz_xxxx, gs_z_xyyzz_xxxy, gs_z_xyyzz_xxxz, gs_z_xyyzz_xxyy, gs_z_xyyzz_xxyz, gs_z_xyyzz_xxzz, gs_z_xyyzz_xyyy, gs_z_xyyzz_xyyz, gs_z_xyyzz_xyzz, gs_z_xyyzz_xzzz, gs_z_xyyzz_yyyy, gs_z_xyyzz_yyyz, gs_z_xyyzz_yyzz, gs_z_xyyzz_yzzz, gs_z_xyyzz_zzzz, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxzz, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyzz, ts_xyyz_xzzz, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyzz, ts_xyyz_yzzz, ts_xyyz_zzzz, ts_xyyzz_xxx, ts_xyyzz_xxxx, ts_xyyzz_xxxy, ts_xyyzz_xxxz, ts_xyyzz_xxy, ts_xyyzz_xxyy, ts_xyyzz_xxyz, ts_xyyzz_xxz, ts_xyyzz_xxzz, ts_xyyzz_xyy, ts_xyyzz_xyyy, ts_xyyzz_xyyz, ts_xyyzz_xyz, ts_xyyzz_xyzz, ts_xyyzz_xzz, ts_xyyzz_xzzz, ts_xyyzz_yyy, ts_xyyzz_yyyy, ts_xyyzz_yyyz, ts_xyyzz_yyz, ts_xyyzz_yyzz, ts_xyyzz_yzz, ts_xyyzz_yzzz, ts_xyyzz_zzz, ts_xyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyzz_xxxx[i] = 4.0 * ts_xyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxy[i] = 4.0 * ts_xyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxz[i] = 4.0 * ts_xyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxyy[i] = 4.0 * ts_xyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxyz[i] = 4.0 * ts_xyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxzz[i] = 4.0 * ts_xyyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyyy[i] = 4.0 * ts_xyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyyz[i] = 4.0 * ts_xyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyzz[i] = 4.0 * ts_xyyz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xzzz[i] = 4.0 * ts_xyyz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyyy[i] = 4.0 * ts_xyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyyz[i] = 4.0 * ts_xyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyzz[i] = 4.0 * ts_xyyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yzzz[i] = 4.0 * ts_xyyz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_zzzz[i] = 4.0 * ts_xyyz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 825-840 components of targeted buffer : HG

    auto gs_z_xyzzz_xxxx = pbuffer.data(idx_g_hg + 825);

    auto gs_z_xyzzz_xxxy = pbuffer.data(idx_g_hg + 826);

    auto gs_z_xyzzz_xxxz = pbuffer.data(idx_g_hg + 827);

    auto gs_z_xyzzz_xxyy = pbuffer.data(idx_g_hg + 828);

    auto gs_z_xyzzz_xxyz = pbuffer.data(idx_g_hg + 829);

    auto gs_z_xyzzz_xxzz = pbuffer.data(idx_g_hg + 830);

    auto gs_z_xyzzz_xyyy = pbuffer.data(idx_g_hg + 831);

    auto gs_z_xyzzz_xyyz = pbuffer.data(idx_g_hg + 832);

    auto gs_z_xyzzz_xyzz = pbuffer.data(idx_g_hg + 833);

    auto gs_z_xyzzz_xzzz = pbuffer.data(idx_g_hg + 834);

    auto gs_z_xyzzz_yyyy = pbuffer.data(idx_g_hg + 835);

    auto gs_z_xyzzz_yyyz = pbuffer.data(idx_g_hg + 836);

    auto gs_z_xyzzz_yyzz = pbuffer.data(idx_g_hg + 837);

    auto gs_z_xyzzz_yzzz = pbuffer.data(idx_g_hg + 838);

    auto gs_z_xyzzz_zzzz = pbuffer.data(idx_g_hg + 839);

    #pragma omp simd aligned(gc_z, gs_z_xyzzz_xxxx, gs_z_xyzzz_xxxy, gs_z_xyzzz_xxxz, gs_z_xyzzz_xxyy, gs_z_xyzzz_xxyz, gs_z_xyzzz_xxzz, gs_z_xyzzz_xyyy, gs_z_xyzzz_xyyz, gs_z_xyzzz_xyzz, gs_z_xyzzz_xzzz, gs_z_xyzzz_yyyy, gs_z_xyzzz_yyyz, gs_z_xyzzz_yyzz, gs_z_xyzzz_yzzz, gs_z_xyzzz_zzzz, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxzz, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyzz, ts_xyzz_xzzz, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyzz, ts_xyzz_yzzz, ts_xyzz_zzzz, ts_xyzzz_xxx, ts_xyzzz_xxxx, ts_xyzzz_xxxy, ts_xyzzz_xxxz, ts_xyzzz_xxy, ts_xyzzz_xxyy, ts_xyzzz_xxyz, ts_xyzzz_xxz, ts_xyzzz_xxzz, ts_xyzzz_xyy, ts_xyzzz_xyyy, ts_xyzzz_xyyz, ts_xyzzz_xyz, ts_xyzzz_xyzz, ts_xyzzz_xzz, ts_xyzzz_xzzz, ts_xyzzz_yyy, ts_xyzzz_yyyy, ts_xyzzz_yyyz, ts_xyzzz_yyz, ts_xyzzz_yyzz, ts_xyzzz_yzz, ts_xyzzz_yzzz, ts_xyzzz_zzz, ts_xyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzzz_xxxx[i] = 6.0 * ts_xyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxy[i] = 6.0 * ts_xyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxz[i] = 6.0 * ts_xyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxyy[i] = 6.0 * ts_xyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxyz[i] = 6.0 * ts_xyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxzz[i] = 6.0 * ts_xyzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyyy[i] = 6.0 * ts_xyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyyz[i] = 6.0 * ts_xyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyzz[i] = 6.0 * ts_xyzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xzzz[i] = 6.0 * ts_xyzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyyy[i] = 6.0 * ts_xyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyyz[i] = 6.0 * ts_xyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyzz[i] = 6.0 * ts_xyzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yzzz[i] = 6.0 * ts_xyzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_zzzz[i] = 6.0 * ts_xyzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 840-855 components of targeted buffer : HG

    auto gs_z_xzzzz_xxxx = pbuffer.data(idx_g_hg + 840);

    auto gs_z_xzzzz_xxxy = pbuffer.data(idx_g_hg + 841);

    auto gs_z_xzzzz_xxxz = pbuffer.data(idx_g_hg + 842);

    auto gs_z_xzzzz_xxyy = pbuffer.data(idx_g_hg + 843);

    auto gs_z_xzzzz_xxyz = pbuffer.data(idx_g_hg + 844);

    auto gs_z_xzzzz_xxzz = pbuffer.data(idx_g_hg + 845);

    auto gs_z_xzzzz_xyyy = pbuffer.data(idx_g_hg + 846);

    auto gs_z_xzzzz_xyyz = pbuffer.data(idx_g_hg + 847);

    auto gs_z_xzzzz_xyzz = pbuffer.data(idx_g_hg + 848);

    auto gs_z_xzzzz_xzzz = pbuffer.data(idx_g_hg + 849);

    auto gs_z_xzzzz_yyyy = pbuffer.data(idx_g_hg + 850);

    auto gs_z_xzzzz_yyyz = pbuffer.data(idx_g_hg + 851);

    auto gs_z_xzzzz_yyzz = pbuffer.data(idx_g_hg + 852);

    auto gs_z_xzzzz_yzzz = pbuffer.data(idx_g_hg + 853);

    auto gs_z_xzzzz_zzzz = pbuffer.data(idx_g_hg + 854);

    #pragma omp simd aligned(gc_z, gs_z_xzzzz_xxxx, gs_z_xzzzz_xxxy, gs_z_xzzzz_xxxz, gs_z_xzzzz_xxyy, gs_z_xzzzz_xxyz, gs_z_xzzzz_xxzz, gs_z_xzzzz_xyyy, gs_z_xzzzz_xyyz, gs_z_xzzzz_xyzz, gs_z_xzzzz_xzzz, gs_z_xzzzz_yyyy, gs_z_xzzzz_yyyz, gs_z_xzzzz_yyzz, gs_z_xzzzz_yzzz, gs_z_xzzzz_zzzz, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxzz, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyzz, ts_xzzz_xzzz, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyzz, ts_xzzz_yzzz, ts_xzzz_zzzz, ts_xzzzz_xxx, ts_xzzzz_xxxx, ts_xzzzz_xxxy, ts_xzzzz_xxxz, ts_xzzzz_xxy, ts_xzzzz_xxyy, ts_xzzzz_xxyz, ts_xzzzz_xxz, ts_xzzzz_xxzz, ts_xzzzz_xyy, ts_xzzzz_xyyy, ts_xzzzz_xyyz, ts_xzzzz_xyz, ts_xzzzz_xyzz, ts_xzzzz_xzz, ts_xzzzz_xzzz, ts_xzzzz_yyy, ts_xzzzz_yyyy, ts_xzzzz_yyyz, ts_xzzzz_yyz, ts_xzzzz_yyzz, ts_xzzzz_yzz, ts_xzzzz_yzzz, ts_xzzzz_zzz, ts_xzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzzz_xxxx[i] = 8.0 * ts_xzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxy[i] = 8.0 * ts_xzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxz[i] = 8.0 * ts_xzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxyy[i] = 8.0 * ts_xzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxyz[i] = 8.0 * ts_xzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxzz[i] = 8.0 * ts_xzzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyyy[i] = 8.0 * ts_xzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyyz[i] = 8.0 * ts_xzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyzz[i] = 8.0 * ts_xzzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xzzz[i] = 8.0 * ts_xzzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyyy[i] = 8.0 * ts_xzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyyz[i] = 8.0 * ts_xzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyzz[i] = 8.0 * ts_xzzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yzzz[i] = 8.0 * ts_xzzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_zzzz[i] = 8.0 * ts_xzzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 855-870 components of targeted buffer : HG

    auto gs_z_yyyyy_xxxx = pbuffer.data(idx_g_hg + 855);

    auto gs_z_yyyyy_xxxy = pbuffer.data(idx_g_hg + 856);

    auto gs_z_yyyyy_xxxz = pbuffer.data(idx_g_hg + 857);

    auto gs_z_yyyyy_xxyy = pbuffer.data(idx_g_hg + 858);

    auto gs_z_yyyyy_xxyz = pbuffer.data(idx_g_hg + 859);

    auto gs_z_yyyyy_xxzz = pbuffer.data(idx_g_hg + 860);

    auto gs_z_yyyyy_xyyy = pbuffer.data(idx_g_hg + 861);

    auto gs_z_yyyyy_xyyz = pbuffer.data(idx_g_hg + 862);

    auto gs_z_yyyyy_xyzz = pbuffer.data(idx_g_hg + 863);

    auto gs_z_yyyyy_xzzz = pbuffer.data(idx_g_hg + 864);

    auto gs_z_yyyyy_yyyy = pbuffer.data(idx_g_hg + 865);

    auto gs_z_yyyyy_yyyz = pbuffer.data(idx_g_hg + 866);

    auto gs_z_yyyyy_yyzz = pbuffer.data(idx_g_hg + 867);

    auto gs_z_yyyyy_yzzz = pbuffer.data(idx_g_hg + 868);

    auto gs_z_yyyyy_zzzz = pbuffer.data(idx_g_hg + 869);

    #pragma omp simd aligned(gc_z, gs_z_yyyyy_xxxx, gs_z_yyyyy_xxxy, gs_z_yyyyy_xxxz, gs_z_yyyyy_xxyy, gs_z_yyyyy_xxyz, gs_z_yyyyy_xxzz, gs_z_yyyyy_xyyy, gs_z_yyyyy_xyyz, gs_z_yyyyy_xyzz, gs_z_yyyyy_xzzz, gs_z_yyyyy_yyyy, gs_z_yyyyy_yyyz, gs_z_yyyyy_yyzz, gs_z_yyyyy_yzzz, gs_z_yyyyy_zzzz, ts_yyyyy_xxx, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxxz, ts_yyyyy_xxy, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xxz, ts_yyyyy_xxzz, ts_yyyyy_xyy, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyz, ts_yyyyy_xyzz, ts_yyyyy_xzz, ts_yyyyy_xzzz, ts_yyyyy_yyy, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyz, ts_yyyyy_yyzz, ts_yyyyy_yzz, ts_yyyyy_yzzz, ts_yyyyy_zzz, ts_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyy_xxxx[i] = 2.0 * ts_yyyyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxy[i] = 2.0 * ts_yyyyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxz[i] = 2.0 * ts_yyyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxyy[i] = 2.0 * ts_yyyyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxyz[i] = 2.0 * ts_yyyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxzz[i] = 4.0 * ts_yyyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyyy[i] = 2.0 * ts_yyyyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyyz[i] = 2.0 * ts_yyyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyzz[i] = 4.0 * ts_yyyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xzzz[i] = 6.0 * ts_yyyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyyy[i] = 2.0 * ts_yyyyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyyz[i] = 2.0 * ts_yyyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyzz[i] = 4.0 * ts_yyyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yzzz[i] = 6.0 * ts_yyyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_zzzz[i] = 8.0 * ts_yyyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 870-885 components of targeted buffer : HG

    auto gs_z_yyyyz_xxxx = pbuffer.data(idx_g_hg + 870);

    auto gs_z_yyyyz_xxxy = pbuffer.data(idx_g_hg + 871);

    auto gs_z_yyyyz_xxxz = pbuffer.data(idx_g_hg + 872);

    auto gs_z_yyyyz_xxyy = pbuffer.data(idx_g_hg + 873);

    auto gs_z_yyyyz_xxyz = pbuffer.data(idx_g_hg + 874);

    auto gs_z_yyyyz_xxzz = pbuffer.data(idx_g_hg + 875);

    auto gs_z_yyyyz_xyyy = pbuffer.data(idx_g_hg + 876);

    auto gs_z_yyyyz_xyyz = pbuffer.data(idx_g_hg + 877);

    auto gs_z_yyyyz_xyzz = pbuffer.data(idx_g_hg + 878);

    auto gs_z_yyyyz_xzzz = pbuffer.data(idx_g_hg + 879);

    auto gs_z_yyyyz_yyyy = pbuffer.data(idx_g_hg + 880);

    auto gs_z_yyyyz_yyyz = pbuffer.data(idx_g_hg + 881);

    auto gs_z_yyyyz_yyzz = pbuffer.data(idx_g_hg + 882);

    auto gs_z_yyyyz_yzzz = pbuffer.data(idx_g_hg + 883);

    auto gs_z_yyyyz_zzzz = pbuffer.data(idx_g_hg + 884);

    #pragma omp simd aligned(gc_z, gs_z_yyyyz_xxxx, gs_z_yyyyz_xxxy, gs_z_yyyyz_xxxz, gs_z_yyyyz_xxyy, gs_z_yyyyz_xxyz, gs_z_yyyyz_xxzz, gs_z_yyyyz_xyyy, gs_z_yyyyz_xyyz, gs_z_yyyyz_xyzz, gs_z_yyyyz_xzzz, gs_z_yyyyz_yyyy, gs_z_yyyyz_yyyz, gs_z_yyyyz_yyzz, gs_z_yyyyz_yzzz, gs_z_yyyyz_zzzz, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxzz, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyzz, ts_yyyy_xzzz, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyzz, ts_yyyy_yzzz, ts_yyyy_zzzz, ts_yyyyz_xxx, ts_yyyyz_xxxx, ts_yyyyz_xxxy, ts_yyyyz_xxxz, ts_yyyyz_xxy, ts_yyyyz_xxyy, ts_yyyyz_xxyz, ts_yyyyz_xxz, ts_yyyyz_xxzz, ts_yyyyz_xyy, ts_yyyyz_xyyy, ts_yyyyz_xyyz, ts_yyyyz_xyz, ts_yyyyz_xyzz, ts_yyyyz_xzz, ts_yyyyz_xzzz, ts_yyyyz_yyy, ts_yyyyz_yyyy, ts_yyyyz_yyyz, ts_yyyyz_yyz, ts_yyyyz_yyzz, ts_yyyyz_yzz, ts_yyyyz_yzzz, ts_yyyyz_zzz, ts_yyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyz_xxxx[i] = 2.0 * ts_yyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxy[i] = 2.0 * ts_yyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxz[i] = 2.0 * ts_yyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxyy[i] = 2.0 * ts_yyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxyz[i] = 2.0 * ts_yyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxzz[i] = 2.0 * ts_yyyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyyy[i] = 2.0 * ts_yyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyyz[i] = 2.0 * ts_yyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyzz[i] = 2.0 * ts_yyyy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xzzz[i] = 2.0 * ts_yyyy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyyy[i] = 2.0 * ts_yyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyyz[i] = 2.0 * ts_yyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyzz[i] = 2.0 * ts_yyyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yzzz[i] = 2.0 * ts_yyyy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_zzzz[i] = 2.0 * ts_yyyy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 885-900 components of targeted buffer : HG

    auto gs_z_yyyzz_xxxx = pbuffer.data(idx_g_hg + 885);

    auto gs_z_yyyzz_xxxy = pbuffer.data(idx_g_hg + 886);

    auto gs_z_yyyzz_xxxz = pbuffer.data(idx_g_hg + 887);

    auto gs_z_yyyzz_xxyy = pbuffer.data(idx_g_hg + 888);

    auto gs_z_yyyzz_xxyz = pbuffer.data(idx_g_hg + 889);

    auto gs_z_yyyzz_xxzz = pbuffer.data(idx_g_hg + 890);

    auto gs_z_yyyzz_xyyy = pbuffer.data(idx_g_hg + 891);

    auto gs_z_yyyzz_xyyz = pbuffer.data(idx_g_hg + 892);

    auto gs_z_yyyzz_xyzz = pbuffer.data(idx_g_hg + 893);

    auto gs_z_yyyzz_xzzz = pbuffer.data(idx_g_hg + 894);

    auto gs_z_yyyzz_yyyy = pbuffer.data(idx_g_hg + 895);

    auto gs_z_yyyzz_yyyz = pbuffer.data(idx_g_hg + 896);

    auto gs_z_yyyzz_yyzz = pbuffer.data(idx_g_hg + 897);

    auto gs_z_yyyzz_yzzz = pbuffer.data(idx_g_hg + 898);

    auto gs_z_yyyzz_zzzz = pbuffer.data(idx_g_hg + 899);

    #pragma omp simd aligned(gc_z, gs_z_yyyzz_xxxx, gs_z_yyyzz_xxxy, gs_z_yyyzz_xxxz, gs_z_yyyzz_xxyy, gs_z_yyyzz_xxyz, gs_z_yyyzz_xxzz, gs_z_yyyzz_xyyy, gs_z_yyyzz_xyyz, gs_z_yyyzz_xyzz, gs_z_yyyzz_xzzz, gs_z_yyyzz_yyyy, gs_z_yyyzz_yyyz, gs_z_yyyzz_yyzz, gs_z_yyyzz_yzzz, gs_z_yyyzz_zzzz, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxzz, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyzz, ts_yyyz_xzzz, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyzz, ts_yyyz_yzzz, ts_yyyz_zzzz, ts_yyyzz_xxx, ts_yyyzz_xxxx, ts_yyyzz_xxxy, ts_yyyzz_xxxz, ts_yyyzz_xxy, ts_yyyzz_xxyy, ts_yyyzz_xxyz, ts_yyyzz_xxz, ts_yyyzz_xxzz, ts_yyyzz_xyy, ts_yyyzz_xyyy, ts_yyyzz_xyyz, ts_yyyzz_xyz, ts_yyyzz_xyzz, ts_yyyzz_xzz, ts_yyyzz_xzzz, ts_yyyzz_yyy, ts_yyyzz_yyyy, ts_yyyzz_yyyz, ts_yyyzz_yyz, ts_yyyzz_yyzz, ts_yyyzz_yzz, ts_yyyzz_yzzz, ts_yyyzz_zzz, ts_yyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyzz_xxxx[i] = 4.0 * ts_yyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxy[i] = 4.0 * ts_yyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxz[i] = 4.0 * ts_yyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxyy[i] = 4.0 * ts_yyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxyz[i] = 4.0 * ts_yyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxzz[i] = 4.0 * ts_yyyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyyy[i] = 4.0 * ts_yyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyyz[i] = 4.0 * ts_yyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyzz[i] = 4.0 * ts_yyyz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xzzz[i] = 4.0 * ts_yyyz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyyy[i] = 4.0 * ts_yyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyyz[i] = 4.0 * ts_yyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyzz[i] = 4.0 * ts_yyyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yzzz[i] = 4.0 * ts_yyyz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_zzzz[i] = 4.0 * ts_yyyz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 900-915 components of targeted buffer : HG

    auto gs_z_yyzzz_xxxx = pbuffer.data(idx_g_hg + 900);

    auto gs_z_yyzzz_xxxy = pbuffer.data(idx_g_hg + 901);

    auto gs_z_yyzzz_xxxz = pbuffer.data(idx_g_hg + 902);

    auto gs_z_yyzzz_xxyy = pbuffer.data(idx_g_hg + 903);

    auto gs_z_yyzzz_xxyz = pbuffer.data(idx_g_hg + 904);

    auto gs_z_yyzzz_xxzz = pbuffer.data(idx_g_hg + 905);

    auto gs_z_yyzzz_xyyy = pbuffer.data(idx_g_hg + 906);

    auto gs_z_yyzzz_xyyz = pbuffer.data(idx_g_hg + 907);

    auto gs_z_yyzzz_xyzz = pbuffer.data(idx_g_hg + 908);

    auto gs_z_yyzzz_xzzz = pbuffer.data(idx_g_hg + 909);

    auto gs_z_yyzzz_yyyy = pbuffer.data(idx_g_hg + 910);

    auto gs_z_yyzzz_yyyz = pbuffer.data(idx_g_hg + 911);

    auto gs_z_yyzzz_yyzz = pbuffer.data(idx_g_hg + 912);

    auto gs_z_yyzzz_yzzz = pbuffer.data(idx_g_hg + 913);

    auto gs_z_yyzzz_zzzz = pbuffer.data(idx_g_hg + 914);

    #pragma omp simd aligned(gc_z, gs_z_yyzzz_xxxx, gs_z_yyzzz_xxxy, gs_z_yyzzz_xxxz, gs_z_yyzzz_xxyy, gs_z_yyzzz_xxyz, gs_z_yyzzz_xxzz, gs_z_yyzzz_xyyy, gs_z_yyzzz_xyyz, gs_z_yyzzz_xyzz, gs_z_yyzzz_xzzz, gs_z_yyzzz_yyyy, gs_z_yyzzz_yyyz, gs_z_yyzzz_yyzz, gs_z_yyzzz_yzzz, gs_z_yyzzz_zzzz, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxzz, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyzz, ts_yyzz_xzzz, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyzz, ts_yyzz_yzzz, ts_yyzz_zzzz, ts_yyzzz_xxx, ts_yyzzz_xxxx, ts_yyzzz_xxxy, ts_yyzzz_xxxz, ts_yyzzz_xxy, ts_yyzzz_xxyy, ts_yyzzz_xxyz, ts_yyzzz_xxz, ts_yyzzz_xxzz, ts_yyzzz_xyy, ts_yyzzz_xyyy, ts_yyzzz_xyyz, ts_yyzzz_xyz, ts_yyzzz_xyzz, ts_yyzzz_xzz, ts_yyzzz_xzzz, ts_yyzzz_yyy, ts_yyzzz_yyyy, ts_yyzzz_yyyz, ts_yyzzz_yyz, ts_yyzzz_yyzz, ts_yyzzz_yzz, ts_yyzzz_yzzz, ts_yyzzz_zzz, ts_yyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzzz_xxxx[i] = 6.0 * ts_yyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxy[i] = 6.0 * ts_yyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxz[i] = 6.0 * ts_yyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxyy[i] = 6.0 * ts_yyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxyz[i] = 6.0 * ts_yyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxzz[i] = 6.0 * ts_yyzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyyy[i] = 6.0 * ts_yyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyyz[i] = 6.0 * ts_yyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyzz[i] = 6.0 * ts_yyzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xzzz[i] = 6.0 * ts_yyzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyyy[i] = 6.0 * ts_yyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyyz[i] = 6.0 * ts_yyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyzz[i] = 6.0 * ts_yyzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yzzz[i] = 6.0 * ts_yyzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_zzzz[i] = 6.0 * ts_yyzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 915-930 components of targeted buffer : HG

    auto gs_z_yzzzz_xxxx = pbuffer.data(idx_g_hg + 915);

    auto gs_z_yzzzz_xxxy = pbuffer.data(idx_g_hg + 916);

    auto gs_z_yzzzz_xxxz = pbuffer.data(idx_g_hg + 917);

    auto gs_z_yzzzz_xxyy = pbuffer.data(idx_g_hg + 918);

    auto gs_z_yzzzz_xxyz = pbuffer.data(idx_g_hg + 919);

    auto gs_z_yzzzz_xxzz = pbuffer.data(idx_g_hg + 920);

    auto gs_z_yzzzz_xyyy = pbuffer.data(idx_g_hg + 921);

    auto gs_z_yzzzz_xyyz = pbuffer.data(idx_g_hg + 922);

    auto gs_z_yzzzz_xyzz = pbuffer.data(idx_g_hg + 923);

    auto gs_z_yzzzz_xzzz = pbuffer.data(idx_g_hg + 924);

    auto gs_z_yzzzz_yyyy = pbuffer.data(idx_g_hg + 925);

    auto gs_z_yzzzz_yyyz = pbuffer.data(idx_g_hg + 926);

    auto gs_z_yzzzz_yyzz = pbuffer.data(idx_g_hg + 927);

    auto gs_z_yzzzz_yzzz = pbuffer.data(idx_g_hg + 928);

    auto gs_z_yzzzz_zzzz = pbuffer.data(idx_g_hg + 929);

    #pragma omp simd aligned(gc_z, gs_z_yzzzz_xxxx, gs_z_yzzzz_xxxy, gs_z_yzzzz_xxxz, gs_z_yzzzz_xxyy, gs_z_yzzzz_xxyz, gs_z_yzzzz_xxzz, gs_z_yzzzz_xyyy, gs_z_yzzzz_xyyz, gs_z_yzzzz_xyzz, gs_z_yzzzz_xzzz, gs_z_yzzzz_yyyy, gs_z_yzzzz_yyyz, gs_z_yzzzz_yyzz, gs_z_yzzzz_yzzz, gs_z_yzzzz_zzzz, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxzz, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyzz, ts_yzzz_xzzz, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyzz, ts_yzzz_yzzz, ts_yzzz_zzzz, ts_yzzzz_xxx, ts_yzzzz_xxxx, ts_yzzzz_xxxy, ts_yzzzz_xxxz, ts_yzzzz_xxy, ts_yzzzz_xxyy, ts_yzzzz_xxyz, ts_yzzzz_xxz, ts_yzzzz_xxzz, ts_yzzzz_xyy, ts_yzzzz_xyyy, ts_yzzzz_xyyz, ts_yzzzz_xyz, ts_yzzzz_xyzz, ts_yzzzz_xzz, ts_yzzzz_xzzz, ts_yzzzz_yyy, ts_yzzzz_yyyy, ts_yzzzz_yyyz, ts_yzzzz_yyz, ts_yzzzz_yyzz, ts_yzzzz_yzz, ts_yzzzz_yzzz, ts_yzzzz_zzz, ts_yzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzzz_xxxx[i] = 8.0 * ts_yzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxy[i] = 8.0 * ts_yzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxz[i] = 8.0 * ts_yzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxyy[i] = 8.0 * ts_yzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxyz[i] = 8.0 * ts_yzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxzz[i] = 8.0 * ts_yzzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyyy[i] = 8.0 * ts_yzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyyz[i] = 8.0 * ts_yzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyzz[i] = 8.0 * ts_yzzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xzzz[i] = 8.0 * ts_yzzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyyy[i] = 8.0 * ts_yzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyyz[i] = 8.0 * ts_yzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyzz[i] = 8.0 * ts_yzzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yzzz[i] = 8.0 * ts_yzzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_zzzz[i] = 8.0 * ts_yzzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 930-945 components of targeted buffer : HG

    auto gs_z_zzzzz_xxxx = pbuffer.data(idx_g_hg + 930);

    auto gs_z_zzzzz_xxxy = pbuffer.data(idx_g_hg + 931);

    auto gs_z_zzzzz_xxxz = pbuffer.data(idx_g_hg + 932);

    auto gs_z_zzzzz_xxyy = pbuffer.data(idx_g_hg + 933);

    auto gs_z_zzzzz_xxyz = pbuffer.data(idx_g_hg + 934);

    auto gs_z_zzzzz_xxzz = pbuffer.data(idx_g_hg + 935);

    auto gs_z_zzzzz_xyyy = pbuffer.data(idx_g_hg + 936);

    auto gs_z_zzzzz_xyyz = pbuffer.data(idx_g_hg + 937);

    auto gs_z_zzzzz_xyzz = pbuffer.data(idx_g_hg + 938);

    auto gs_z_zzzzz_xzzz = pbuffer.data(idx_g_hg + 939);

    auto gs_z_zzzzz_yyyy = pbuffer.data(idx_g_hg + 940);

    auto gs_z_zzzzz_yyyz = pbuffer.data(idx_g_hg + 941);

    auto gs_z_zzzzz_yyzz = pbuffer.data(idx_g_hg + 942);

    auto gs_z_zzzzz_yzzz = pbuffer.data(idx_g_hg + 943);

    auto gs_z_zzzzz_zzzz = pbuffer.data(idx_g_hg + 944);

    #pragma omp simd aligned(gc_z, gs_z_zzzzz_xxxx, gs_z_zzzzz_xxxy, gs_z_zzzzz_xxxz, gs_z_zzzzz_xxyy, gs_z_zzzzz_xxyz, gs_z_zzzzz_xxzz, gs_z_zzzzz_xyyy, gs_z_zzzzz_xyyz, gs_z_zzzzz_xyzz, gs_z_zzzzz_xzzz, gs_z_zzzzz_yyyy, gs_z_zzzzz_yyyz, gs_z_zzzzz_yyzz, gs_z_zzzzz_yzzz, gs_z_zzzzz_zzzz, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, ts_zzzzz_xxx, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxy, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxz, ts_zzzzz_xxzz, ts_zzzzz_xyy, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyz, ts_zzzzz_xyzz, ts_zzzzz_xzz, ts_zzzzz_xzzz, ts_zzzzz_yyy, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyz, ts_zzzzz_yyzz, ts_zzzzz_yzz, ts_zzzzz_yzzz, ts_zzzzz_zzz, ts_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzzz_xxxx[i] = 10.0 * ts_zzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxy[i] = 10.0 * ts_zzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxz[i] = 10.0 * ts_zzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxyy[i] = 10.0 * ts_zzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxyz[i] = 10.0 * ts_zzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxzz[i] = 10.0 * ts_zzzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyyy[i] = 10.0 * ts_zzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyyz[i] = 10.0 * ts_zzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyzz[i] = 10.0 * ts_zzzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xzzz[i] = 10.0 * ts_zzzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyyy[i] = 10.0 * ts_zzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyyz[i] = 10.0 * ts_zzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyzz[i] = 10.0 * ts_zzzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yzzz[i] = 10.0 * ts_zzzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_zzzz[i] = 10.0 * ts_zzzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_zzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

