#include "KineticEnergyPrimRecIG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ig(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ig,
                            const size_t              idx_ovl_gg,
                            const size_t              idx_kin_gg,
                            const size_t              idx_kin_hf,
                            const size_t              idx_kin_hg,
                            const size_t              idx_ovl_ig,
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

    // Set up components of auxiliary buffer : GG

    auto ts_xxxx_xxxx = pbuffer.data(idx_ovl_gg);

    auto ts_xxxx_xxxy = pbuffer.data(idx_ovl_gg + 1);

    auto ts_xxxx_xxxz = pbuffer.data(idx_ovl_gg + 2);

    auto ts_xxxx_xxyy = pbuffer.data(idx_ovl_gg + 3);

    auto ts_xxxx_xxyz = pbuffer.data(idx_ovl_gg + 4);

    auto ts_xxxx_xxzz = pbuffer.data(idx_ovl_gg + 5);

    auto ts_xxxx_xyyy = pbuffer.data(idx_ovl_gg + 6);

    auto ts_xxxx_xyyz = pbuffer.data(idx_ovl_gg + 7);

    auto ts_xxxx_xyzz = pbuffer.data(idx_ovl_gg + 8);

    auto ts_xxxx_xzzz = pbuffer.data(idx_ovl_gg + 9);

    auto ts_xxxx_yyyy = pbuffer.data(idx_ovl_gg + 10);

    auto ts_xxxx_yyyz = pbuffer.data(idx_ovl_gg + 11);

    auto ts_xxxx_yyzz = pbuffer.data(idx_ovl_gg + 12);

    auto ts_xxxx_yzzz = pbuffer.data(idx_ovl_gg + 13);

    auto ts_xxxx_zzzz = pbuffer.data(idx_ovl_gg + 14);

    auto ts_xxxy_xxxx = pbuffer.data(idx_ovl_gg + 15);

    auto ts_xxxy_xxxz = pbuffer.data(idx_ovl_gg + 17);

    auto ts_xxxy_xxzz = pbuffer.data(idx_ovl_gg + 20);

    auto ts_xxxy_xzzz = pbuffer.data(idx_ovl_gg + 24);

    auto ts_xxxz_xxxx = pbuffer.data(idx_ovl_gg + 30);

    auto ts_xxxz_xxxy = pbuffer.data(idx_ovl_gg + 31);

    auto ts_xxxz_xxyy = pbuffer.data(idx_ovl_gg + 33);

    auto ts_xxxz_xyyy = pbuffer.data(idx_ovl_gg + 36);

    auto ts_xxyy_xxxx = pbuffer.data(idx_ovl_gg + 45);

    auto ts_xxyy_xxxy = pbuffer.data(idx_ovl_gg + 46);

    auto ts_xxyy_xxxz = pbuffer.data(idx_ovl_gg + 47);

    auto ts_xxyy_xxyy = pbuffer.data(idx_ovl_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_ovl_gg + 49);

    auto ts_xxyy_xxzz = pbuffer.data(idx_ovl_gg + 50);

    auto ts_xxyy_xyyy = pbuffer.data(idx_ovl_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_ovl_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_ovl_gg + 53);

    auto ts_xxyy_xzzz = pbuffer.data(idx_ovl_gg + 54);

    auto ts_xxyy_yyyy = pbuffer.data(idx_ovl_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_ovl_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_ovl_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_ovl_gg + 58);

    auto ts_xxyy_zzzz = pbuffer.data(idx_ovl_gg + 59);

    auto ts_xxzz_xxxx = pbuffer.data(idx_ovl_gg + 75);

    auto ts_xxzz_xxxy = pbuffer.data(idx_ovl_gg + 76);

    auto ts_xxzz_xxxz = pbuffer.data(idx_ovl_gg + 77);

    auto ts_xxzz_xxyy = pbuffer.data(idx_ovl_gg + 78);

    auto ts_xxzz_xxyz = pbuffer.data(idx_ovl_gg + 79);

    auto ts_xxzz_xxzz = pbuffer.data(idx_ovl_gg + 80);

    auto ts_xxzz_xyyy = pbuffer.data(idx_ovl_gg + 81);

    auto ts_xxzz_xyyz = pbuffer.data(idx_ovl_gg + 82);

    auto ts_xxzz_xyzz = pbuffer.data(idx_ovl_gg + 83);

    auto ts_xxzz_xzzz = pbuffer.data(idx_ovl_gg + 84);

    auto ts_xxzz_yyyy = pbuffer.data(idx_ovl_gg + 85);

    auto ts_xxzz_yyyz = pbuffer.data(idx_ovl_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_ovl_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_ovl_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_ovl_gg + 89);

    auto ts_xyyy_xxxy = pbuffer.data(idx_ovl_gg + 91);

    auto ts_xyyy_xxyy = pbuffer.data(idx_ovl_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_ovl_gg + 94);

    auto ts_xyyy_xyyy = pbuffer.data(idx_ovl_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_ovl_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_ovl_gg + 98);

    auto ts_xyyy_yyyy = pbuffer.data(idx_ovl_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_ovl_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_ovl_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_ovl_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_ovl_gg + 104);

    auto ts_xzzz_xxxz = pbuffer.data(idx_ovl_gg + 137);

    auto ts_xzzz_xxyz = pbuffer.data(idx_ovl_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_ovl_gg + 140);

    auto ts_xzzz_xyyz = pbuffer.data(idx_ovl_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_ovl_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_ovl_gg + 144);

    auto ts_xzzz_yyyy = pbuffer.data(idx_ovl_gg + 145);

    auto ts_xzzz_yyyz = pbuffer.data(idx_ovl_gg + 146);

    auto ts_xzzz_yyzz = pbuffer.data(idx_ovl_gg + 147);

    auto ts_xzzz_yzzz = pbuffer.data(idx_ovl_gg + 148);

    auto ts_xzzz_zzzz = pbuffer.data(idx_ovl_gg + 149);

    auto ts_yyyy_xxxx = pbuffer.data(idx_ovl_gg + 150);

    auto ts_yyyy_xxxy = pbuffer.data(idx_ovl_gg + 151);

    auto ts_yyyy_xxxz = pbuffer.data(idx_ovl_gg + 152);

    auto ts_yyyy_xxyy = pbuffer.data(idx_ovl_gg + 153);

    auto ts_yyyy_xxyz = pbuffer.data(idx_ovl_gg + 154);

    auto ts_yyyy_xxzz = pbuffer.data(idx_ovl_gg + 155);

    auto ts_yyyy_xyyy = pbuffer.data(idx_ovl_gg + 156);

    auto ts_yyyy_xyyz = pbuffer.data(idx_ovl_gg + 157);

    auto ts_yyyy_xyzz = pbuffer.data(idx_ovl_gg + 158);

    auto ts_yyyy_xzzz = pbuffer.data(idx_ovl_gg + 159);

    auto ts_yyyy_yyyy = pbuffer.data(idx_ovl_gg + 160);

    auto ts_yyyy_yyyz = pbuffer.data(idx_ovl_gg + 161);

    auto ts_yyyy_yyzz = pbuffer.data(idx_ovl_gg + 162);

    auto ts_yyyy_yzzz = pbuffer.data(idx_ovl_gg + 163);

    auto ts_yyyy_zzzz = pbuffer.data(idx_ovl_gg + 164);

    auto ts_yyyz_xxxy = pbuffer.data(idx_ovl_gg + 166);

    auto ts_yyyz_xxyy = pbuffer.data(idx_ovl_gg + 168);

    auto ts_yyyz_xyyy = pbuffer.data(idx_ovl_gg + 171);

    auto ts_yyyz_yyyy = pbuffer.data(idx_ovl_gg + 175);

    auto ts_yyzz_xxxx = pbuffer.data(idx_ovl_gg + 180);

    auto ts_yyzz_xxxy = pbuffer.data(idx_ovl_gg + 181);

    auto ts_yyzz_xxxz = pbuffer.data(idx_ovl_gg + 182);

    auto ts_yyzz_xxyy = pbuffer.data(idx_ovl_gg + 183);

    auto ts_yyzz_xxyz = pbuffer.data(idx_ovl_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_ovl_gg + 185);

    auto ts_yyzz_xyyy = pbuffer.data(idx_ovl_gg + 186);

    auto ts_yyzz_xyyz = pbuffer.data(idx_ovl_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_ovl_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_ovl_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_ovl_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_ovl_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_ovl_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_ovl_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_ovl_gg + 194);

    auto ts_yzzz_xxxx = pbuffer.data(idx_ovl_gg + 195);

    auto ts_yzzz_xxxz = pbuffer.data(idx_ovl_gg + 197);

    auto ts_yzzz_xxyz = pbuffer.data(idx_ovl_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_ovl_gg + 200);

    auto ts_yzzz_xyyz = pbuffer.data(idx_ovl_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_ovl_gg + 203);

    auto ts_yzzz_xzzz = pbuffer.data(idx_ovl_gg + 204);

    auto ts_yzzz_yyyz = pbuffer.data(idx_ovl_gg + 206);

    auto ts_yzzz_yyzz = pbuffer.data(idx_ovl_gg + 207);

    auto ts_yzzz_yzzz = pbuffer.data(idx_ovl_gg + 208);

    auto ts_yzzz_zzzz = pbuffer.data(idx_ovl_gg + 209);

    auto ts_zzzz_xxxx = pbuffer.data(idx_ovl_gg + 210);

    auto ts_zzzz_xxxy = pbuffer.data(idx_ovl_gg + 211);

    auto ts_zzzz_xxxz = pbuffer.data(idx_ovl_gg + 212);

    auto ts_zzzz_xxyy = pbuffer.data(idx_ovl_gg + 213);

    auto ts_zzzz_xxyz = pbuffer.data(idx_ovl_gg + 214);

    auto ts_zzzz_xxzz = pbuffer.data(idx_ovl_gg + 215);

    auto ts_zzzz_xyyy = pbuffer.data(idx_ovl_gg + 216);

    auto ts_zzzz_xyyz = pbuffer.data(idx_ovl_gg + 217);

    auto ts_zzzz_xyzz = pbuffer.data(idx_ovl_gg + 218);

    auto ts_zzzz_xzzz = pbuffer.data(idx_ovl_gg + 219);

    auto ts_zzzz_yyyy = pbuffer.data(idx_ovl_gg + 220);

    auto ts_zzzz_yyyz = pbuffer.data(idx_ovl_gg + 221);

    auto ts_zzzz_yyzz = pbuffer.data(idx_ovl_gg + 222);

    auto ts_zzzz_yzzz = pbuffer.data(idx_ovl_gg + 223);

    auto ts_zzzz_zzzz = pbuffer.data(idx_ovl_gg + 224);

    // Set up components of auxiliary buffer : GG

    auto tk_xxxx_xxxx = pbuffer.data(idx_kin_gg);

    auto tk_xxxx_xxxy = pbuffer.data(idx_kin_gg + 1);

    auto tk_xxxx_xxxz = pbuffer.data(idx_kin_gg + 2);

    auto tk_xxxx_xxyy = pbuffer.data(idx_kin_gg + 3);

    auto tk_xxxx_xxyz = pbuffer.data(idx_kin_gg + 4);

    auto tk_xxxx_xxzz = pbuffer.data(idx_kin_gg + 5);

    auto tk_xxxx_xyyy = pbuffer.data(idx_kin_gg + 6);

    auto tk_xxxx_xyyz = pbuffer.data(idx_kin_gg + 7);

    auto tk_xxxx_xyzz = pbuffer.data(idx_kin_gg + 8);

    auto tk_xxxx_xzzz = pbuffer.data(idx_kin_gg + 9);

    auto tk_xxxx_yyyy = pbuffer.data(idx_kin_gg + 10);

    auto tk_xxxx_yyyz = pbuffer.data(idx_kin_gg + 11);

    auto tk_xxxx_yyzz = pbuffer.data(idx_kin_gg + 12);

    auto tk_xxxx_yzzz = pbuffer.data(idx_kin_gg + 13);

    auto tk_xxxx_zzzz = pbuffer.data(idx_kin_gg + 14);

    auto tk_xxxy_xxxx = pbuffer.data(idx_kin_gg + 15);

    auto tk_xxxy_xxxz = pbuffer.data(idx_kin_gg + 17);

    auto tk_xxxy_xxzz = pbuffer.data(idx_kin_gg + 20);

    auto tk_xxxy_xzzz = pbuffer.data(idx_kin_gg + 24);

    auto tk_xxxz_xxxx = pbuffer.data(idx_kin_gg + 30);

    auto tk_xxxz_xxxy = pbuffer.data(idx_kin_gg + 31);

    auto tk_xxxz_xxyy = pbuffer.data(idx_kin_gg + 33);

    auto tk_xxxz_xyyy = pbuffer.data(idx_kin_gg + 36);

    auto tk_xxyy_xxxx = pbuffer.data(idx_kin_gg + 45);

    auto tk_xxyy_xxxy = pbuffer.data(idx_kin_gg + 46);

    auto tk_xxyy_xxxz = pbuffer.data(idx_kin_gg + 47);

    auto tk_xxyy_xxyy = pbuffer.data(idx_kin_gg + 48);

    auto tk_xxyy_xxyz = pbuffer.data(idx_kin_gg + 49);

    auto tk_xxyy_xxzz = pbuffer.data(idx_kin_gg + 50);

    auto tk_xxyy_xyyy = pbuffer.data(idx_kin_gg + 51);

    auto tk_xxyy_xyyz = pbuffer.data(idx_kin_gg + 52);

    auto tk_xxyy_xyzz = pbuffer.data(idx_kin_gg + 53);

    auto tk_xxyy_xzzz = pbuffer.data(idx_kin_gg + 54);

    auto tk_xxyy_yyyy = pbuffer.data(idx_kin_gg + 55);

    auto tk_xxyy_yyyz = pbuffer.data(idx_kin_gg + 56);

    auto tk_xxyy_yyzz = pbuffer.data(idx_kin_gg + 57);

    auto tk_xxyy_yzzz = pbuffer.data(idx_kin_gg + 58);

    auto tk_xxyy_zzzz = pbuffer.data(idx_kin_gg + 59);

    auto tk_xxzz_xxxx = pbuffer.data(idx_kin_gg + 75);

    auto tk_xxzz_xxxy = pbuffer.data(idx_kin_gg + 76);

    auto tk_xxzz_xxxz = pbuffer.data(idx_kin_gg + 77);

    auto tk_xxzz_xxyy = pbuffer.data(idx_kin_gg + 78);

    auto tk_xxzz_xxyz = pbuffer.data(idx_kin_gg + 79);

    auto tk_xxzz_xxzz = pbuffer.data(idx_kin_gg + 80);

    auto tk_xxzz_xyyy = pbuffer.data(idx_kin_gg + 81);

    auto tk_xxzz_xyyz = pbuffer.data(idx_kin_gg + 82);

    auto tk_xxzz_xyzz = pbuffer.data(idx_kin_gg + 83);

    auto tk_xxzz_xzzz = pbuffer.data(idx_kin_gg + 84);

    auto tk_xxzz_yyyy = pbuffer.data(idx_kin_gg + 85);

    auto tk_xxzz_yyyz = pbuffer.data(idx_kin_gg + 86);

    auto tk_xxzz_yyzz = pbuffer.data(idx_kin_gg + 87);

    auto tk_xxzz_yzzz = pbuffer.data(idx_kin_gg + 88);

    auto tk_xxzz_zzzz = pbuffer.data(idx_kin_gg + 89);

    auto tk_xyyy_xxxy = pbuffer.data(idx_kin_gg + 91);

    auto tk_xyyy_xxyy = pbuffer.data(idx_kin_gg + 93);

    auto tk_xyyy_xxyz = pbuffer.data(idx_kin_gg + 94);

    auto tk_xyyy_xyyy = pbuffer.data(idx_kin_gg + 96);

    auto tk_xyyy_xyyz = pbuffer.data(idx_kin_gg + 97);

    auto tk_xyyy_xyzz = pbuffer.data(idx_kin_gg + 98);

    auto tk_xyyy_yyyy = pbuffer.data(idx_kin_gg + 100);

    auto tk_xyyy_yyyz = pbuffer.data(idx_kin_gg + 101);

    auto tk_xyyy_yyzz = pbuffer.data(idx_kin_gg + 102);

    auto tk_xyyy_yzzz = pbuffer.data(idx_kin_gg + 103);

    auto tk_xyyy_zzzz = pbuffer.data(idx_kin_gg + 104);

    auto tk_xzzz_xxxz = pbuffer.data(idx_kin_gg + 137);

    auto tk_xzzz_xxyz = pbuffer.data(idx_kin_gg + 139);

    auto tk_xzzz_xxzz = pbuffer.data(idx_kin_gg + 140);

    auto tk_xzzz_xyyz = pbuffer.data(idx_kin_gg + 142);

    auto tk_xzzz_xyzz = pbuffer.data(idx_kin_gg + 143);

    auto tk_xzzz_xzzz = pbuffer.data(idx_kin_gg + 144);

    auto tk_xzzz_yyyy = pbuffer.data(idx_kin_gg + 145);

    auto tk_xzzz_yyyz = pbuffer.data(idx_kin_gg + 146);

    auto tk_xzzz_yyzz = pbuffer.data(idx_kin_gg + 147);

    auto tk_xzzz_yzzz = pbuffer.data(idx_kin_gg + 148);

    auto tk_xzzz_zzzz = pbuffer.data(idx_kin_gg + 149);

    auto tk_yyyy_xxxx = pbuffer.data(idx_kin_gg + 150);

    auto tk_yyyy_xxxy = pbuffer.data(idx_kin_gg + 151);

    auto tk_yyyy_xxxz = pbuffer.data(idx_kin_gg + 152);

    auto tk_yyyy_xxyy = pbuffer.data(idx_kin_gg + 153);

    auto tk_yyyy_xxyz = pbuffer.data(idx_kin_gg + 154);

    auto tk_yyyy_xxzz = pbuffer.data(idx_kin_gg + 155);

    auto tk_yyyy_xyyy = pbuffer.data(idx_kin_gg + 156);

    auto tk_yyyy_xyyz = pbuffer.data(idx_kin_gg + 157);

    auto tk_yyyy_xyzz = pbuffer.data(idx_kin_gg + 158);

    auto tk_yyyy_xzzz = pbuffer.data(idx_kin_gg + 159);

    auto tk_yyyy_yyyy = pbuffer.data(idx_kin_gg + 160);

    auto tk_yyyy_yyyz = pbuffer.data(idx_kin_gg + 161);

    auto tk_yyyy_yyzz = pbuffer.data(idx_kin_gg + 162);

    auto tk_yyyy_yzzz = pbuffer.data(idx_kin_gg + 163);

    auto tk_yyyy_zzzz = pbuffer.data(idx_kin_gg + 164);

    auto tk_yyyz_xxxy = pbuffer.data(idx_kin_gg + 166);

    auto tk_yyyz_xxyy = pbuffer.data(idx_kin_gg + 168);

    auto tk_yyyz_xyyy = pbuffer.data(idx_kin_gg + 171);

    auto tk_yyyz_yyyy = pbuffer.data(idx_kin_gg + 175);

    auto tk_yyzz_xxxx = pbuffer.data(idx_kin_gg + 180);

    auto tk_yyzz_xxxy = pbuffer.data(idx_kin_gg + 181);

    auto tk_yyzz_xxxz = pbuffer.data(idx_kin_gg + 182);

    auto tk_yyzz_xxyy = pbuffer.data(idx_kin_gg + 183);

    auto tk_yyzz_xxyz = pbuffer.data(idx_kin_gg + 184);

    auto tk_yyzz_xxzz = pbuffer.data(idx_kin_gg + 185);

    auto tk_yyzz_xyyy = pbuffer.data(idx_kin_gg + 186);

    auto tk_yyzz_xyyz = pbuffer.data(idx_kin_gg + 187);

    auto tk_yyzz_xyzz = pbuffer.data(idx_kin_gg + 188);

    auto tk_yyzz_xzzz = pbuffer.data(idx_kin_gg + 189);

    auto tk_yyzz_yyyy = pbuffer.data(idx_kin_gg + 190);

    auto tk_yyzz_yyyz = pbuffer.data(idx_kin_gg + 191);

    auto tk_yyzz_yyzz = pbuffer.data(idx_kin_gg + 192);

    auto tk_yyzz_yzzz = pbuffer.data(idx_kin_gg + 193);

    auto tk_yyzz_zzzz = pbuffer.data(idx_kin_gg + 194);

    auto tk_yzzz_xxxx = pbuffer.data(idx_kin_gg + 195);

    auto tk_yzzz_xxxz = pbuffer.data(idx_kin_gg + 197);

    auto tk_yzzz_xxyz = pbuffer.data(idx_kin_gg + 199);

    auto tk_yzzz_xxzz = pbuffer.data(idx_kin_gg + 200);

    auto tk_yzzz_xyyz = pbuffer.data(idx_kin_gg + 202);

    auto tk_yzzz_xyzz = pbuffer.data(idx_kin_gg + 203);

    auto tk_yzzz_xzzz = pbuffer.data(idx_kin_gg + 204);

    auto tk_yzzz_yyyz = pbuffer.data(idx_kin_gg + 206);

    auto tk_yzzz_yyzz = pbuffer.data(idx_kin_gg + 207);

    auto tk_yzzz_yzzz = pbuffer.data(idx_kin_gg + 208);

    auto tk_yzzz_zzzz = pbuffer.data(idx_kin_gg + 209);

    auto tk_zzzz_xxxx = pbuffer.data(idx_kin_gg + 210);

    auto tk_zzzz_xxxy = pbuffer.data(idx_kin_gg + 211);

    auto tk_zzzz_xxxz = pbuffer.data(idx_kin_gg + 212);

    auto tk_zzzz_xxyy = pbuffer.data(idx_kin_gg + 213);

    auto tk_zzzz_xxyz = pbuffer.data(idx_kin_gg + 214);

    auto tk_zzzz_xxzz = pbuffer.data(idx_kin_gg + 215);

    auto tk_zzzz_xyyy = pbuffer.data(idx_kin_gg + 216);

    auto tk_zzzz_xyyz = pbuffer.data(idx_kin_gg + 217);

    auto tk_zzzz_xyzz = pbuffer.data(idx_kin_gg + 218);

    auto tk_zzzz_xzzz = pbuffer.data(idx_kin_gg + 219);

    auto tk_zzzz_yyyy = pbuffer.data(idx_kin_gg + 220);

    auto tk_zzzz_yyyz = pbuffer.data(idx_kin_gg + 221);

    auto tk_zzzz_yyzz = pbuffer.data(idx_kin_gg + 222);

    auto tk_zzzz_yzzz = pbuffer.data(idx_kin_gg + 223);

    auto tk_zzzz_zzzz = pbuffer.data(idx_kin_gg + 224);

    // Set up components of auxiliary buffer : HF

    auto tk_xxxxx_xxx = pbuffer.data(idx_kin_hf);

    auto tk_xxxxx_xxy = pbuffer.data(idx_kin_hf + 1);

    auto tk_xxxxx_xxz = pbuffer.data(idx_kin_hf + 2);

    auto tk_xxxxx_xyy = pbuffer.data(idx_kin_hf + 3);

    auto tk_xxxxx_xyz = pbuffer.data(idx_kin_hf + 4);

    auto tk_xxxxx_xzz = pbuffer.data(idx_kin_hf + 5);

    auto tk_xxxxx_yyy = pbuffer.data(idx_kin_hf + 6);

    auto tk_xxxxx_yyz = pbuffer.data(idx_kin_hf + 7);

    auto tk_xxxxx_yzz = pbuffer.data(idx_kin_hf + 8);

    auto tk_xxxxx_zzz = pbuffer.data(idx_kin_hf + 9);

    auto tk_xxxxz_xxz = pbuffer.data(idx_kin_hf + 22);

    auto tk_xxxxz_xyz = pbuffer.data(idx_kin_hf + 24);

    auto tk_xxxxz_xzz = pbuffer.data(idx_kin_hf + 25);

    auto tk_xxxxz_yyz = pbuffer.data(idx_kin_hf + 27);

    auto tk_xxxxz_yzz = pbuffer.data(idx_kin_hf + 28);

    auto tk_xxxxz_zzz = pbuffer.data(idx_kin_hf + 29);

    auto tk_xxxyy_xxx = pbuffer.data(idx_kin_hf + 30);

    auto tk_xxxyy_xxy = pbuffer.data(idx_kin_hf + 31);

    auto tk_xxxyy_xxz = pbuffer.data(idx_kin_hf + 32);

    auto tk_xxxyy_xyy = pbuffer.data(idx_kin_hf + 33);

    auto tk_xxxyy_xyz = pbuffer.data(idx_kin_hf + 34);

    auto tk_xxxyy_xzz = pbuffer.data(idx_kin_hf + 35);

    auto tk_xxxyy_yyy = pbuffer.data(idx_kin_hf + 36);

    auto tk_xxxyy_yyz = pbuffer.data(idx_kin_hf + 37);

    auto tk_xxxyy_yzz = pbuffer.data(idx_kin_hf + 38);

    auto tk_xxxyy_zzz = pbuffer.data(idx_kin_hf + 39);

    auto tk_xxxzz_xxx = pbuffer.data(idx_kin_hf + 50);

    auto tk_xxxzz_xxy = pbuffer.data(idx_kin_hf + 51);

    auto tk_xxxzz_xxz = pbuffer.data(idx_kin_hf + 52);

    auto tk_xxxzz_xyy = pbuffer.data(idx_kin_hf + 53);

    auto tk_xxxzz_xyz = pbuffer.data(idx_kin_hf + 54);

    auto tk_xxxzz_xzz = pbuffer.data(idx_kin_hf + 55);

    auto tk_xxxzz_yyy = pbuffer.data(idx_kin_hf + 56);

    auto tk_xxxzz_yyz = pbuffer.data(idx_kin_hf + 57);

    auto tk_xxxzz_yzz = pbuffer.data(idx_kin_hf + 58);

    auto tk_xxxzz_zzz = pbuffer.data(idx_kin_hf + 59);

    auto tk_xxyyy_xxx = pbuffer.data(idx_kin_hf + 60);

    auto tk_xxyyy_xxy = pbuffer.data(idx_kin_hf + 61);

    auto tk_xxyyy_xxz = pbuffer.data(idx_kin_hf + 62);

    auto tk_xxyyy_xyy = pbuffer.data(idx_kin_hf + 63);

    auto tk_xxyyy_xyz = pbuffer.data(idx_kin_hf + 64);

    auto tk_xxyyy_xzz = pbuffer.data(idx_kin_hf + 65);

    auto tk_xxyyy_yyy = pbuffer.data(idx_kin_hf + 66);

    auto tk_xxyyy_yyz = pbuffer.data(idx_kin_hf + 67);

    auto tk_xxyyy_yzz = pbuffer.data(idx_kin_hf + 68);

    auto tk_xxyyy_zzz = pbuffer.data(idx_kin_hf + 69);

    auto tk_xxzzz_xxx = pbuffer.data(idx_kin_hf + 90);

    auto tk_xxzzz_xxy = pbuffer.data(idx_kin_hf + 91);

    auto tk_xxzzz_xxz = pbuffer.data(idx_kin_hf + 92);

    auto tk_xxzzz_xyy = pbuffer.data(idx_kin_hf + 93);

    auto tk_xxzzz_xyz = pbuffer.data(idx_kin_hf + 94);

    auto tk_xxzzz_xzz = pbuffer.data(idx_kin_hf + 95);

    auto tk_xxzzz_yyy = pbuffer.data(idx_kin_hf + 96);

    auto tk_xxzzz_yyz = pbuffer.data(idx_kin_hf + 97);

    auto tk_xxzzz_yzz = pbuffer.data(idx_kin_hf + 98);

    auto tk_xxzzz_zzz = pbuffer.data(idx_kin_hf + 99);

    auto tk_xyyyy_xxy = pbuffer.data(idx_kin_hf + 101);

    auto tk_xyyyy_xyy = pbuffer.data(idx_kin_hf + 103);

    auto tk_xyyyy_xyz = pbuffer.data(idx_kin_hf + 104);

    auto tk_xyyyy_yyy = pbuffer.data(idx_kin_hf + 106);

    auto tk_xyyyy_yyz = pbuffer.data(idx_kin_hf + 107);

    auto tk_xyyyy_yzz = pbuffer.data(idx_kin_hf + 108);

    auto tk_xyyzz_xyz = pbuffer.data(idx_kin_hf + 124);

    auto tk_xyyzz_yyz = pbuffer.data(idx_kin_hf + 127);

    auto tk_xyyzz_yzz = pbuffer.data(idx_kin_hf + 128);

    auto tk_xzzzz_xxz = pbuffer.data(idx_kin_hf + 142);

    auto tk_xzzzz_xyz = pbuffer.data(idx_kin_hf + 144);

    auto tk_xzzzz_xzz = pbuffer.data(idx_kin_hf + 145);

    auto tk_xzzzz_yyz = pbuffer.data(idx_kin_hf + 147);

    auto tk_xzzzz_yzz = pbuffer.data(idx_kin_hf + 148);

    auto tk_xzzzz_zzz = pbuffer.data(idx_kin_hf + 149);

    auto tk_yyyyy_xxx = pbuffer.data(idx_kin_hf + 150);

    auto tk_yyyyy_xxy = pbuffer.data(idx_kin_hf + 151);

    auto tk_yyyyy_xxz = pbuffer.data(idx_kin_hf + 152);

    auto tk_yyyyy_xyy = pbuffer.data(idx_kin_hf + 153);

    auto tk_yyyyy_xyz = pbuffer.data(idx_kin_hf + 154);

    auto tk_yyyyy_xzz = pbuffer.data(idx_kin_hf + 155);

    auto tk_yyyyy_yyy = pbuffer.data(idx_kin_hf + 156);

    auto tk_yyyyy_yyz = pbuffer.data(idx_kin_hf + 157);

    auto tk_yyyyy_yzz = pbuffer.data(idx_kin_hf + 158);

    auto tk_yyyyy_zzz = pbuffer.data(idx_kin_hf + 159);

    auto tk_yyyyz_xxz = pbuffer.data(idx_kin_hf + 162);

    auto tk_yyyyz_xyz = pbuffer.data(idx_kin_hf + 164);

    auto tk_yyyyz_xzz = pbuffer.data(idx_kin_hf + 165);

    auto tk_yyyyz_yyz = pbuffer.data(idx_kin_hf + 167);

    auto tk_yyyyz_yzz = pbuffer.data(idx_kin_hf + 168);

    auto tk_yyyyz_zzz = pbuffer.data(idx_kin_hf + 169);

    auto tk_yyyzz_xxx = pbuffer.data(idx_kin_hf + 170);

    auto tk_yyyzz_xxy = pbuffer.data(idx_kin_hf + 171);

    auto tk_yyyzz_xxz = pbuffer.data(idx_kin_hf + 172);

    auto tk_yyyzz_xyy = pbuffer.data(idx_kin_hf + 173);

    auto tk_yyyzz_xyz = pbuffer.data(idx_kin_hf + 174);

    auto tk_yyyzz_xzz = pbuffer.data(idx_kin_hf + 175);

    auto tk_yyyzz_yyy = pbuffer.data(idx_kin_hf + 176);

    auto tk_yyyzz_yyz = pbuffer.data(idx_kin_hf + 177);

    auto tk_yyyzz_yzz = pbuffer.data(idx_kin_hf + 178);

    auto tk_yyyzz_zzz = pbuffer.data(idx_kin_hf + 179);

    auto tk_yyzzz_xxx = pbuffer.data(idx_kin_hf + 180);

    auto tk_yyzzz_xxy = pbuffer.data(idx_kin_hf + 181);

    auto tk_yyzzz_xxz = pbuffer.data(idx_kin_hf + 182);

    auto tk_yyzzz_xyy = pbuffer.data(idx_kin_hf + 183);

    auto tk_yyzzz_xyz = pbuffer.data(idx_kin_hf + 184);

    auto tk_yyzzz_xzz = pbuffer.data(idx_kin_hf + 185);

    auto tk_yyzzz_yyy = pbuffer.data(idx_kin_hf + 186);

    auto tk_yyzzz_yyz = pbuffer.data(idx_kin_hf + 187);

    auto tk_yyzzz_yzz = pbuffer.data(idx_kin_hf + 188);

    auto tk_yyzzz_zzz = pbuffer.data(idx_kin_hf + 189);

    auto tk_yzzzz_xxy = pbuffer.data(idx_kin_hf + 191);

    auto tk_yzzzz_xxz = pbuffer.data(idx_kin_hf + 192);

    auto tk_yzzzz_xyy = pbuffer.data(idx_kin_hf + 193);

    auto tk_yzzzz_xyz = pbuffer.data(idx_kin_hf + 194);

    auto tk_yzzzz_xzz = pbuffer.data(idx_kin_hf + 195);

    auto tk_yzzzz_yyy = pbuffer.data(idx_kin_hf + 196);

    auto tk_yzzzz_yyz = pbuffer.data(idx_kin_hf + 197);

    auto tk_yzzzz_yzz = pbuffer.data(idx_kin_hf + 198);

    auto tk_yzzzz_zzz = pbuffer.data(idx_kin_hf + 199);

    auto tk_zzzzz_xxx = pbuffer.data(idx_kin_hf + 200);

    auto tk_zzzzz_xxy = pbuffer.data(idx_kin_hf + 201);

    auto tk_zzzzz_xxz = pbuffer.data(idx_kin_hf + 202);

    auto tk_zzzzz_xyy = pbuffer.data(idx_kin_hf + 203);

    auto tk_zzzzz_xyz = pbuffer.data(idx_kin_hf + 204);

    auto tk_zzzzz_xzz = pbuffer.data(idx_kin_hf + 205);

    auto tk_zzzzz_yyy = pbuffer.data(idx_kin_hf + 206);

    auto tk_zzzzz_yyz = pbuffer.data(idx_kin_hf + 207);

    auto tk_zzzzz_yzz = pbuffer.data(idx_kin_hf + 208);

    auto tk_zzzzz_zzz = pbuffer.data(idx_kin_hf + 209);

    // Set up components of auxiliary buffer : HG

    auto tk_xxxxx_xxxx = pbuffer.data(idx_kin_hg);

    auto tk_xxxxx_xxxy = pbuffer.data(idx_kin_hg + 1);

    auto tk_xxxxx_xxxz = pbuffer.data(idx_kin_hg + 2);

    auto tk_xxxxx_xxyy = pbuffer.data(idx_kin_hg + 3);

    auto tk_xxxxx_xxyz = pbuffer.data(idx_kin_hg + 4);

    auto tk_xxxxx_xxzz = pbuffer.data(idx_kin_hg + 5);

    auto tk_xxxxx_xyyy = pbuffer.data(idx_kin_hg + 6);

    auto tk_xxxxx_xyyz = pbuffer.data(idx_kin_hg + 7);

    auto tk_xxxxx_xyzz = pbuffer.data(idx_kin_hg + 8);

    auto tk_xxxxx_xzzz = pbuffer.data(idx_kin_hg + 9);

    auto tk_xxxxx_yyyy = pbuffer.data(idx_kin_hg + 10);

    auto tk_xxxxx_yyyz = pbuffer.data(idx_kin_hg + 11);

    auto tk_xxxxx_yyzz = pbuffer.data(idx_kin_hg + 12);

    auto tk_xxxxx_yzzz = pbuffer.data(idx_kin_hg + 13);

    auto tk_xxxxx_zzzz = pbuffer.data(idx_kin_hg + 14);

    auto tk_xxxxy_xxxx = pbuffer.data(idx_kin_hg + 15);

    auto tk_xxxxy_xxxy = pbuffer.data(idx_kin_hg + 16);

    auto tk_xxxxy_xxxz = pbuffer.data(idx_kin_hg + 17);

    auto tk_xxxxy_xxyy = pbuffer.data(idx_kin_hg + 18);

    auto tk_xxxxy_xxzz = pbuffer.data(idx_kin_hg + 20);

    auto tk_xxxxy_xyyy = pbuffer.data(idx_kin_hg + 21);

    auto tk_xxxxy_xzzz = pbuffer.data(idx_kin_hg + 24);

    auto tk_xxxxy_yyyy = pbuffer.data(idx_kin_hg + 25);

    auto tk_xxxxz_xxxx = pbuffer.data(idx_kin_hg + 30);

    auto tk_xxxxz_xxxy = pbuffer.data(idx_kin_hg + 31);

    auto tk_xxxxz_xxxz = pbuffer.data(idx_kin_hg + 32);

    auto tk_xxxxz_xxyy = pbuffer.data(idx_kin_hg + 33);

    auto tk_xxxxz_xxyz = pbuffer.data(idx_kin_hg + 34);

    auto tk_xxxxz_xxzz = pbuffer.data(idx_kin_hg + 35);

    auto tk_xxxxz_xyyy = pbuffer.data(idx_kin_hg + 36);

    auto tk_xxxxz_xyyz = pbuffer.data(idx_kin_hg + 37);

    auto tk_xxxxz_xyzz = pbuffer.data(idx_kin_hg + 38);

    auto tk_xxxxz_xzzz = pbuffer.data(idx_kin_hg + 39);

    auto tk_xxxxz_yyyz = pbuffer.data(idx_kin_hg + 41);

    auto tk_xxxxz_yyzz = pbuffer.data(idx_kin_hg + 42);

    auto tk_xxxxz_yzzz = pbuffer.data(idx_kin_hg + 43);

    auto tk_xxxxz_zzzz = pbuffer.data(idx_kin_hg + 44);

    auto tk_xxxyy_xxxx = pbuffer.data(idx_kin_hg + 45);

    auto tk_xxxyy_xxxy = pbuffer.data(idx_kin_hg + 46);

    auto tk_xxxyy_xxxz = pbuffer.data(idx_kin_hg + 47);

    auto tk_xxxyy_xxyy = pbuffer.data(idx_kin_hg + 48);

    auto tk_xxxyy_xxyz = pbuffer.data(idx_kin_hg + 49);

    auto tk_xxxyy_xxzz = pbuffer.data(idx_kin_hg + 50);

    auto tk_xxxyy_xyyy = pbuffer.data(idx_kin_hg + 51);

    auto tk_xxxyy_xyyz = pbuffer.data(idx_kin_hg + 52);

    auto tk_xxxyy_xyzz = pbuffer.data(idx_kin_hg + 53);

    auto tk_xxxyy_xzzz = pbuffer.data(idx_kin_hg + 54);

    auto tk_xxxyy_yyyy = pbuffer.data(idx_kin_hg + 55);

    auto tk_xxxyy_yyyz = pbuffer.data(idx_kin_hg + 56);

    auto tk_xxxyy_yyzz = pbuffer.data(idx_kin_hg + 57);

    auto tk_xxxyy_yzzz = pbuffer.data(idx_kin_hg + 58);

    auto tk_xxxyy_zzzz = pbuffer.data(idx_kin_hg + 59);

    auto tk_xxxzz_xxxx = pbuffer.data(idx_kin_hg + 75);

    auto tk_xxxzz_xxxy = pbuffer.data(idx_kin_hg + 76);

    auto tk_xxxzz_xxxz = pbuffer.data(idx_kin_hg + 77);

    auto tk_xxxzz_xxyy = pbuffer.data(idx_kin_hg + 78);

    auto tk_xxxzz_xxyz = pbuffer.data(idx_kin_hg + 79);

    auto tk_xxxzz_xxzz = pbuffer.data(idx_kin_hg + 80);

    auto tk_xxxzz_xyyy = pbuffer.data(idx_kin_hg + 81);

    auto tk_xxxzz_xyyz = pbuffer.data(idx_kin_hg + 82);

    auto tk_xxxzz_xyzz = pbuffer.data(idx_kin_hg + 83);

    auto tk_xxxzz_xzzz = pbuffer.data(idx_kin_hg + 84);

    auto tk_xxxzz_yyyy = pbuffer.data(idx_kin_hg + 85);

    auto tk_xxxzz_yyyz = pbuffer.data(idx_kin_hg + 86);

    auto tk_xxxzz_yyzz = pbuffer.data(idx_kin_hg + 87);

    auto tk_xxxzz_yzzz = pbuffer.data(idx_kin_hg + 88);

    auto tk_xxxzz_zzzz = pbuffer.data(idx_kin_hg + 89);

    auto tk_xxyyy_xxxx = pbuffer.data(idx_kin_hg + 90);

    auto tk_xxyyy_xxxy = pbuffer.data(idx_kin_hg + 91);

    auto tk_xxyyy_xxxz = pbuffer.data(idx_kin_hg + 92);

    auto tk_xxyyy_xxyy = pbuffer.data(idx_kin_hg + 93);

    auto tk_xxyyy_xxyz = pbuffer.data(idx_kin_hg + 94);

    auto tk_xxyyy_xxzz = pbuffer.data(idx_kin_hg + 95);

    auto tk_xxyyy_xyyy = pbuffer.data(idx_kin_hg + 96);

    auto tk_xxyyy_xyyz = pbuffer.data(idx_kin_hg + 97);

    auto tk_xxyyy_xyzz = pbuffer.data(idx_kin_hg + 98);

    auto tk_xxyyy_xzzz = pbuffer.data(idx_kin_hg + 99);

    auto tk_xxyyy_yyyy = pbuffer.data(idx_kin_hg + 100);

    auto tk_xxyyy_yyyz = pbuffer.data(idx_kin_hg + 101);

    auto tk_xxyyy_yyzz = pbuffer.data(idx_kin_hg + 102);

    auto tk_xxyyy_yzzz = pbuffer.data(idx_kin_hg + 103);

    auto tk_xxyyy_zzzz = pbuffer.data(idx_kin_hg + 104);

    auto tk_xxyyz_xxxy = pbuffer.data(idx_kin_hg + 106);

    auto tk_xxyyz_xxyy = pbuffer.data(idx_kin_hg + 108);

    auto tk_xxyyz_xyyy = pbuffer.data(idx_kin_hg + 111);

    auto tk_xxyzz_xxxx = pbuffer.data(idx_kin_hg + 120);

    auto tk_xxyzz_xxxz = pbuffer.data(idx_kin_hg + 122);

    auto tk_xxyzz_xxzz = pbuffer.data(idx_kin_hg + 125);

    auto tk_xxyzz_xzzz = pbuffer.data(idx_kin_hg + 129);

    auto tk_xxzzz_xxxx = pbuffer.data(idx_kin_hg + 135);

    auto tk_xxzzz_xxxy = pbuffer.data(idx_kin_hg + 136);

    auto tk_xxzzz_xxxz = pbuffer.data(idx_kin_hg + 137);

    auto tk_xxzzz_xxyy = pbuffer.data(idx_kin_hg + 138);

    auto tk_xxzzz_xxyz = pbuffer.data(idx_kin_hg + 139);

    auto tk_xxzzz_xxzz = pbuffer.data(idx_kin_hg + 140);

    auto tk_xxzzz_xyyy = pbuffer.data(idx_kin_hg + 141);

    auto tk_xxzzz_xyyz = pbuffer.data(idx_kin_hg + 142);

    auto tk_xxzzz_xyzz = pbuffer.data(idx_kin_hg + 143);

    auto tk_xxzzz_xzzz = pbuffer.data(idx_kin_hg + 144);

    auto tk_xxzzz_yyyy = pbuffer.data(idx_kin_hg + 145);

    auto tk_xxzzz_yyyz = pbuffer.data(idx_kin_hg + 146);

    auto tk_xxzzz_yyzz = pbuffer.data(idx_kin_hg + 147);

    auto tk_xxzzz_yzzz = pbuffer.data(idx_kin_hg + 148);

    auto tk_xxzzz_zzzz = pbuffer.data(idx_kin_hg + 149);

    auto tk_xyyyy_xxxx = pbuffer.data(idx_kin_hg + 150);

    auto tk_xyyyy_xxxy = pbuffer.data(idx_kin_hg + 151);

    auto tk_xyyyy_xxyy = pbuffer.data(idx_kin_hg + 153);

    auto tk_xyyyy_xxyz = pbuffer.data(idx_kin_hg + 154);

    auto tk_xyyyy_xyyy = pbuffer.data(idx_kin_hg + 156);

    auto tk_xyyyy_xyyz = pbuffer.data(idx_kin_hg + 157);

    auto tk_xyyyy_xyzz = pbuffer.data(idx_kin_hg + 158);

    auto tk_xyyyy_yyyy = pbuffer.data(idx_kin_hg + 160);

    auto tk_xyyyy_yyyz = pbuffer.data(idx_kin_hg + 161);

    auto tk_xyyyy_yyzz = pbuffer.data(idx_kin_hg + 162);

    auto tk_xyyyy_yzzz = pbuffer.data(idx_kin_hg + 163);

    auto tk_xyyyy_zzzz = pbuffer.data(idx_kin_hg + 164);

    auto tk_xyyzz_xxyz = pbuffer.data(idx_kin_hg + 184);

    auto tk_xyyzz_xyyz = pbuffer.data(idx_kin_hg + 187);

    auto tk_xyyzz_xyzz = pbuffer.data(idx_kin_hg + 188);

    auto tk_xyyzz_yyyy = pbuffer.data(idx_kin_hg + 190);

    auto tk_xyyzz_yyyz = pbuffer.data(idx_kin_hg + 191);

    auto tk_xyyzz_yyzz = pbuffer.data(idx_kin_hg + 192);

    auto tk_xyyzz_yzzz = pbuffer.data(idx_kin_hg + 193);

    auto tk_xyyzz_zzzz = pbuffer.data(idx_kin_hg + 194);

    auto tk_xzzzz_xxxx = pbuffer.data(idx_kin_hg + 210);

    auto tk_xzzzz_xxxz = pbuffer.data(idx_kin_hg + 212);

    auto tk_xzzzz_xxyz = pbuffer.data(idx_kin_hg + 214);

    auto tk_xzzzz_xxzz = pbuffer.data(idx_kin_hg + 215);

    auto tk_xzzzz_xyyz = pbuffer.data(idx_kin_hg + 217);

    auto tk_xzzzz_xyzz = pbuffer.data(idx_kin_hg + 218);

    auto tk_xzzzz_xzzz = pbuffer.data(idx_kin_hg + 219);

    auto tk_xzzzz_yyyy = pbuffer.data(idx_kin_hg + 220);

    auto tk_xzzzz_yyyz = pbuffer.data(idx_kin_hg + 221);

    auto tk_xzzzz_yyzz = pbuffer.data(idx_kin_hg + 222);

    auto tk_xzzzz_yzzz = pbuffer.data(idx_kin_hg + 223);

    auto tk_xzzzz_zzzz = pbuffer.data(idx_kin_hg + 224);

    auto tk_yyyyy_xxxx = pbuffer.data(idx_kin_hg + 225);

    auto tk_yyyyy_xxxy = pbuffer.data(idx_kin_hg + 226);

    auto tk_yyyyy_xxxz = pbuffer.data(idx_kin_hg + 227);

    auto tk_yyyyy_xxyy = pbuffer.data(idx_kin_hg + 228);

    auto tk_yyyyy_xxyz = pbuffer.data(idx_kin_hg + 229);

    auto tk_yyyyy_xxzz = pbuffer.data(idx_kin_hg + 230);

    auto tk_yyyyy_xyyy = pbuffer.data(idx_kin_hg + 231);

    auto tk_yyyyy_xyyz = pbuffer.data(idx_kin_hg + 232);

    auto tk_yyyyy_xyzz = pbuffer.data(idx_kin_hg + 233);

    auto tk_yyyyy_xzzz = pbuffer.data(idx_kin_hg + 234);

    auto tk_yyyyy_yyyy = pbuffer.data(idx_kin_hg + 235);

    auto tk_yyyyy_yyyz = pbuffer.data(idx_kin_hg + 236);

    auto tk_yyyyy_yyzz = pbuffer.data(idx_kin_hg + 237);

    auto tk_yyyyy_yzzz = pbuffer.data(idx_kin_hg + 238);

    auto tk_yyyyy_zzzz = pbuffer.data(idx_kin_hg + 239);

    auto tk_yyyyz_xxxy = pbuffer.data(idx_kin_hg + 241);

    auto tk_yyyyz_xxxz = pbuffer.data(idx_kin_hg + 242);

    auto tk_yyyyz_xxyy = pbuffer.data(idx_kin_hg + 243);

    auto tk_yyyyz_xxyz = pbuffer.data(idx_kin_hg + 244);

    auto tk_yyyyz_xxzz = pbuffer.data(idx_kin_hg + 245);

    auto tk_yyyyz_xyyy = pbuffer.data(idx_kin_hg + 246);

    auto tk_yyyyz_xyyz = pbuffer.data(idx_kin_hg + 247);

    auto tk_yyyyz_xyzz = pbuffer.data(idx_kin_hg + 248);

    auto tk_yyyyz_xzzz = pbuffer.data(idx_kin_hg + 249);

    auto tk_yyyyz_yyyy = pbuffer.data(idx_kin_hg + 250);

    auto tk_yyyyz_yyyz = pbuffer.data(idx_kin_hg + 251);

    auto tk_yyyyz_yyzz = pbuffer.data(idx_kin_hg + 252);

    auto tk_yyyyz_yzzz = pbuffer.data(idx_kin_hg + 253);

    auto tk_yyyyz_zzzz = pbuffer.data(idx_kin_hg + 254);

    auto tk_yyyzz_xxxx = pbuffer.data(idx_kin_hg + 255);

    auto tk_yyyzz_xxxy = pbuffer.data(idx_kin_hg + 256);

    auto tk_yyyzz_xxxz = pbuffer.data(idx_kin_hg + 257);

    auto tk_yyyzz_xxyy = pbuffer.data(idx_kin_hg + 258);

    auto tk_yyyzz_xxyz = pbuffer.data(idx_kin_hg + 259);

    auto tk_yyyzz_xxzz = pbuffer.data(idx_kin_hg + 260);

    auto tk_yyyzz_xyyy = pbuffer.data(idx_kin_hg + 261);

    auto tk_yyyzz_xyyz = pbuffer.data(idx_kin_hg + 262);

    auto tk_yyyzz_xyzz = pbuffer.data(idx_kin_hg + 263);

    auto tk_yyyzz_xzzz = pbuffer.data(idx_kin_hg + 264);

    auto tk_yyyzz_yyyy = pbuffer.data(idx_kin_hg + 265);

    auto tk_yyyzz_yyyz = pbuffer.data(idx_kin_hg + 266);

    auto tk_yyyzz_yyzz = pbuffer.data(idx_kin_hg + 267);

    auto tk_yyyzz_yzzz = pbuffer.data(idx_kin_hg + 268);

    auto tk_yyyzz_zzzz = pbuffer.data(idx_kin_hg + 269);

    auto tk_yyzzz_xxxx = pbuffer.data(idx_kin_hg + 270);

    auto tk_yyzzz_xxxy = pbuffer.data(idx_kin_hg + 271);

    auto tk_yyzzz_xxxz = pbuffer.data(idx_kin_hg + 272);

    auto tk_yyzzz_xxyy = pbuffer.data(idx_kin_hg + 273);

    auto tk_yyzzz_xxyz = pbuffer.data(idx_kin_hg + 274);

    auto tk_yyzzz_xxzz = pbuffer.data(idx_kin_hg + 275);

    auto tk_yyzzz_xyyy = pbuffer.data(idx_kin_hg + 276);

    auto tk_yyzzz_xyyz = pbuffer.data(idx_kin_hg + 277);

    auto tk_yyzzz_xyzz = pbuffer.data(idx_kin_hg + 278);

    auto tk_yyzzz_xzzz = pbuffer.data(idx_kin_hg + 279);

    auto tk_yyzzz_yyyy = pbuffer.data(idx_kin_hg + 280);

    auto tk_yyzzz_yyyz = pbuffer.data(idx_kin_hg + 281);

    auto tk_yyzzz_yyzz = pbuffer.data(idx_kin_hg + 282);

    auto tk_yyzzz_yzzz = pbuffer.data(idx_kin_hg + 283);

    auto tk_yyzzz_zzzz = pbuffer.data(idx_kin_hg + 284);

    auto tk_yzzzz_xxxx = pbuffer.data(idx_kin_hg + 285);

    auto tk_yzzzz_xxxy = pbuffer.data(idx_kin_hg + 286);

    auto tk_yzzzz_xxxz = pbuffer.data(idx_kin_hg + 287);

    auto tk_yzzzz_xxyy = pbuffer.data(idx_kin_hg + 288);

    auto tk_yzzzz_xxyz = pbuffer.data(idx_kin_hg + 289);

    auto tk_yzzzz_xxzz = pbuffer.data(idx_kin_hg + 290);

    auto tk_yzzzz_xyyy = pbuffer.data(idx_kin_hg + 291);

    auto tk_yzzzz_xyyz = pbuffer.data(idx_kin_hg + 292);

    auto tk_yzzzz_xyzz = pbuffer.data(idx_kin_hg + 293);

    auto tk_yzzzz_xzzz = pbuffer.data(idx_kin_hg + 294);

    auto tk_yzzzz_yyyy = pbuffer.data(idx_kin_hg + 295);

    auto tk_yzzzz_yyyz = pbuffer.data(idx_kin_hg + 296);

    auto tk_yzzzz_yyzz = pbuffer.data(idx_kin_hg + 297);

    auto tk_yzzzz_yzzz = pbuffer.data(idx_kin_hg + 298);

    auto tk_yzzzz_zzzz = pbuffer.data(idx_kin_hg + 299);

    auto tk_zzzzz_xxxx = pbuffer.data(idx_kin_hg + 300);

    auto tk_zzzzz_xxxy = pbuffer.data(idx_kin_hg + 301);

    auto tk_zzzzz_xxxz = pbuffer.data(idx_kin_hg + 302);

    auto tk_zzzzz_xxyy = pbuffer.data(idx_kin_hg + 303);

    auto tk_zzzzz_xxyz = pbuffer.data(idx_kin_hg + 304);

    auto tk_zzzzz_xxzz = pbuffer.data(idx_kin_hg + 305);

    auto tk_zzzzz_xyyy = pbuffer.data(idx_kin_hg + 306);

    auto tk_zzzzz_xyyz = pbuffer.data(idx_kin_hg + 307);

    auto tk_zzzzz_xyzz = pbuffer.data(idx_kin_hg + 308);

    auto tk_zzzzz_xzzz = pbuffer.data(idx_kin_hg + 309);

    auto tk_zzzzz_yyyy = pbuffer.data(idx_kin_hg + 310);

    auto tk_zzzzz_yyyz = pbuffer.data(idx_kin_hg + 311);

    auto tk_zzzzz_yyzz = pbuffer.data(idx_kin_hg + 312);

    auto tk_zzzzz_yzzz = pbuffer.data(idx_kin_hg + 313);

    auto tk_zzzzz_zzzz = pbuffer.data(idx_kin_hg + 314);

    // Set up components of auxiliary buffer : IG

    auto ts_xxxxxx_xxxx = pbuffer.data(idx_ovl_ig);

    auto ts_xxxxxx_xxxy = pbuffer.data(idx_ovl_ig + 1);

    auto ts_xxxxxx_xxxz = pbuffer.data(idx_ovl_ig + 2);

    auto ts_xxxxxx_xxyy = pbuffer.data(idx_ovl_ig + 3);

    auto ts_xxxxxx_xxyz = pbuffer.data(idx_ovl_ig + 4);

    auto ts_xxxxxx_xxzz = pbuffer.data(idx_ovl_ig + 5);

    auto ts_xxxxxx_xyyy = pbuffer.data(idx_ovl_ig + 6);

    auto ts_xxxxxx_xyyz = pbuffer.data(idx_ovl_ig + 7);

    auto ts_xxxxxx_xyzz = pbuffer.data(idx_ovl_ig + 8);

    auto ts_xxxxxx_xzzz = pbuffer.data(idx_ovl_ig + 9);

    auto ts_xxxxxx_yyyy = pbuffer.data(idx_ovl_ig + 10);

    auto ts_xxxxxx_yyyz = pbuffer.data(idx_ovl_ig + 11);

    auto ts_xxxxxx_yyzz = pbuffer.data(idx_ovl_ig + 12);

    auto ts_xxxxxx_yzzz = pbuffer.data(idx_ovl_ig + 13);

    auto ts_xxxxxx_zzzz = pbuffer.data(idx_ovl_ig + 14);

    auto ts_xxxxxy_xxxx = pbuffer.data(idx_ovl_ig + 15);

    auto ts_xxxxxy_xxxy = pbuffer.data(idx_ovl_ig + 16);

    auto ts_xxxxxy_xxxz = pbuffer.data(idx_ovl_ig + 17);

    auto ts_xxxxxy_xxyy = pbuffer.data(idx_ovl_ig + 18);

    auto ts_xxxxxy_xxyz = pbuffer.data(idx_ovl_ig + 19);

    auto ts_xxxxxy_xxzz = pbuffer.data(idx_ovl_ig + 20);

    auto ts_xxxxxy_xyyy = pbuffer.data(idx_ovl_ig + 21);

    auto ts_xxxxxy_xyyz = pbuffer.data(idx_ovl_ig + 22);

    auto ts_xxxxxy_xyzz = pbuffer.data(idx_ovl_ig + 23);

    auto ts_xxxxxy_xzzz = pbuffer.data(idx_ovl_ig + 24);

    auto ts_xxxxxy_yyyy = pbuffer.data(idx_ovl_ig + 25);

    auto ts_xxxxxy_yyyz = pbuffer.data(idx_ovl_ig + 26);

    auto ts_xxxxxy_yyzz = pbuffer.data(idx_ovl_ig + 27);

    auto ts_xxxxxy_yzzz = pbuffer.data(idx_ovl_ig + 28);

    auto ts_xxxxxy_zzzz = pbuffer.data(idx_ovl_ig + 29);

    auto ts_xxxxxz_xxxx = pbuffer.data(idx_ovl_ig + 30);

    auto ts_xxxxxz_xxxy = pbuffer.data(idx_ovl_ig + 31);

    auto ts_xxxxxz_xxxz = pbuffer.data(idx_ovl_ig + 32);

    auto ts_xxxxxz_xxyy = pbuffer.data(idx_ovl_ig + 33);

    auto ts_xxxxxz_xxyz = pbuffer.data(idx_ovl_ig + 34);

    auto ts_xxxxxz_xxzz = pbuffer.data(idx_ovl_ig + 35);

    auto ts_xxxxxz_xyyy = pbuffer.data(idx_ovl_ig + 36);

    auto ts_xxxxxz_xyyz = pbuffer.data(idx_ovl_ig + 37);

    auto ts_xxxxxz_xyzz = pbuffer.data(idx_ovl_ig + 38);

    auto ts_xxxxxz_xzzz = pbuffer.data(idx_ovl_ig + 39);

    auto ts_xxxxxz_yyyy = pbuffer.data(idx_ovl_ig + 40);

    auto ts_xxxxxz_yyyz = pbuffer.data(idx_ovl_ig + 41);

    auto ts_xxxxxz_yyzz = pbuffer.data(idx_ovl_ig + 42);

    auto ts_xxxxxz_yzzz = pbuffer.data(idx_ovl_ig + 43);

    auto ts_xxxxxz_zzzz = pbuffer.data(idx_ovl_ig + 44);

    auto ts_xxxxyy_xxxx = pbuffer.data(idx_ovl_ig + 45);

    auto ts_xxxxyy_xxxy = pbuffer.data(idx_ovl_ig + 46);

    auto ts_xxxxyy_xxxz = pbuffer.data(idx_ovl_ig + 47);

    auto ts_xxxxyy_xxyy = pbuffer.data(idx_ovl_ig + 48);

    auto ts_xxxxyy_xxyz = pbuffer.data(idx_ovl_ig + 49);

    auto ts_xxxxyy_xxzz = pbuffer.data(idx_ovl_ig + 50);

    auto ts_xxxxyy_xyyy = pbuffer.data(idx_ovl_ig + 51);

    auto ts_xxxxyy_xyyz = pbuffer.data(idx_ovl_ig + 52);

    auto ts_xxxxyy_xyzz = pbuffer.data(idx_ovl_ig + 53);

    auto ts_xxxxyy_xzzz = pbuffer.data(idx_ovl_ig + 54);

    auto ts_xxxxyy_yyyy = pbuffer.data(idx_ovl_ig + 55);

    auto ts_xxxxyy_yyyz = pbuffer.data(idx_ovl_ig + 56);

    auto ts_xxxxyy_yyzz = pbuffer.data(idx_ovl_ig + 57);

    auto ts_xxxxyy_yzzz = pbuffer.data(idx_ovl_ig + 58);

    auto ts_xxxxyy_zzzz = pbuffer.data(idx_ovl_ig + 59);

    auto ts_xxxxyz_xxxx = pbuffer.data(idx_ovl_ig + 60);

    auto ts_xxxxyz_xxxy = pbuffer.data(idx_ovl_ig + 61);

    auto ts_xxxxyz_xxxz = pbuffer.data(idx_ovl_ig + 62);

    auto ts_xxxxyz_xxyy = pbuffer.data(idx_ovl_ig + 63);

    auto ts_xxxxyz_xxyz = pbuffer.data(idx_ovl_ig + 64);

    auto ts_xxxxyz_xxzz = pbuffer.data(idx_ovl_ig + 65);

    auto ts_xxxxyz_xyyy = pbuffer.data(idx_ovl_ig + 66);

    auto ts_xxxxyz_xyyz = pbuffer.data(idx_ovl_ig + 67);

    auto ts_xxxxyz_xyzz = pbuffer.data(idx_ovl_ig + 68);

    auto ts_xxxxyz_xzzz = pbuffer.data(idx_ovl_ig + 69);

    auto ts_xxxxyz_yyyy = pbuffer.data(idx_ovl_ig + 70);

    auto ts_xxxxyz_yyyz = pbuffer.data(idx_ovl_ig + 71);

    auto ts_xxxxyz_yyzz = pbuffer.data(idx_ovl_ig + 72);

    auto ts_xxxxyz_yzzz = pbuffer.data(idx_ovl_ig + 73);

    auto ts_xxxxyz_zzzz = pbuffer.data(idx_ovl_ig + 74);

    auto ts_xxxxzz_xxxx = pbuffer.data(idx_ovl_ig + 75);

    auto ts_xxxxzz_xxxy = pbuffer.data(idx_ovl_ig + 76);

    auto ts_xxxxzz_xxxz = pbuffer.data(idx_ovl_ig + 77);

    auto ts_xxxxzz_xxyy = pbuffer.data(idx_ovl_ig + 78);

    auto ts_xxxxzz_xxyz = pbuffer.data(idx_ovl_ig + 79);

    auto ts_xxxxzz_xxzz = pbuffer.data(idx_ovl_ig + 80);

    auto ts_xxxxzz_xyyy = pbuffer.data(idx_ovl_ig + 81);

    auto ts_xxxxzz_xyyz = pbuffer.data(idx_ovl_ig + 82);

    auto ts_xxxxzz_xyzz = pbuffer.data(idx_ovl_ig + 83);

    auto ts_xxxxzz_xzzz = pbuffer.data(idx_ovl_ig + 84);

    auto ts_xxxxzz_yyyy = pbuffer.data(idx_ovl_ig + 85);

    auto ts_xxxxzz_yyyz = pbuffer.data(idx_ovl_ig + 86);

    auto ts_xxxxzz_yyzz = pbuffer.data(idx_ovl_ig + 87);

    auto ts_xxxxzz_yzzz = pbuffer.data(idx_ovl_ig + 88);

    auto ts_xxxxzz_zzzz = pbuffer.data(idx_ovl_ig + 89);

    auto ts_xxxyyy_xxxx = pbuffer.data(idx_ovl_ig + 90);

    auto ts_xxxyyy_xxxy = pbuffer.data(idx_ovl_ig + 91);

    auto ts_xxxyyy_xxxz = pbuffer.data(idx_ovl_ig + 92);

    auto ts_xxxyyy_xxyy = pbuffer.data(idx_ovl_ig + 93);

    auto ts_xxxyyy_xxyz = pbuffer.data(idx_ovl_ig + 94);

    auto ts_xxxyyy_xxzz = pbuffer.data(idx_ovl_ig + 95);

    auto ts_xxxyyy_xyyy = pbuffer.data(idx_ovl_ig + 96);

    auto ts_xxxyyy_xyyz = pbuffer.data(idx_ovl_ig + 97);

    auto ts_xxxyyy_xyzz = pbuffer.data(idx_ovl_ig + 98);

    auto ts_xxxyyy_xzzz = pbuffer.data(idx_ovl_ig + 99);

    auto ts_xxxyyy_yyyy = pbuffer.data(idx_ovl_ig + 100);

    auto ts_xxxyyy_yyyz = pbuffer.data(idx_ovl_ig + 101);

    auto ts_xxxyyy_yyzz = pbuffer.data(idx_ovl_ig + 102);

    auto ts_xxxyyy_yzzz = pbuffer.data(idx_ovl_ig + 103);

    auto ts_xxxyyy_zzzz = pbuffer.data(idx_ovl_ig + 104);

    auto ts_xxxyyz_xxxx = pbuffer.data(idx_ovl_ig + 105);

    auto ts_xxxyyz_xxxy = pbuffer.data(idx_ovl_ig + 106);

    auto ts_xxxyyz_xxxz = pbuffer.data(idx_ovl_ig + 107);

    auto ts_xxxyyz_xxyy = pbuffer.data(idx_ovl_ig + 108);

    auto ts_xxxyyz_xxyz = pbuffer.data(idx_ovl_ig + 109);

    auto ts_xxxyyz_xxzz = pbuffer.data(idx_ovl_ig + 110);

    auto ts_xxxyyz_xyyy = pbuffer.data(idx_ovl_ig + 111);

    auto ts_xxxyyz_xyyz = pbuffer.data(idx_ovl_ig + 112);

    auto ts_xxxyyz_xyzz = pbuffer.data(idx_ovl_ig + 113);

    auto ts_xxxyyz_xzzz = pbuffer.data(idx_ovl_ig + 114);

    auto ts_xxxyyz_yyyy = pbuffer.data(idx_ovl_ig + 115);

    auto ts_xxxyyz_yyyz = pbuffer.data(idx_ovl_ig + 116);

    auto ts_xxxyyz_yyzz = pbuffer.data(idx_ovl_ig + 117);

    auto ts_xxxyyz_yzzz = pbuffer.data(idx_ovl_ig + 118);

    auto ts_xxxyyz_zzzz = pbuffer.data(idx_ovl_ig + 119);

    auto ts_xxxyzz_xxxx = pbuffer.data(idx_ovl_ig + 120);

    auto ts_xxxyzz_xxxy = pbuffer.data(idx_ovl_ig + 121);

    auto ts_xxxyzz_xxxz = pbuffer.data(idx_ovl_ig + 122);

    auto ts_xxxyzz_xxyy = pbuffer.data(idx_ovl_ig + 123);

    auto ts_xxxyzz_xxyz = pbuffer.data(idx_ovl_ig + 124);

    auto ts_xxxyzz_xxzz = pbuffer.data(idx_ovl_ig + 125);

    auto ts_xxxyzz_xyyy = pbuffer.data(idx_ovl_ig + 126);

    auto ts_xxxyzz_xyyz = pbuffer.data(idx_ovl_ig + 127);

    auto ts_xxxyzz_xyzz = pbuffer.data(idx_ovl_ig + 128);

    auto ts_xxxyzz_xzzz = pbuffer.data(idx_ovl_ig + 129);

    auto ts_xxxyzz_yyyy = pbuffer.data(idx_ovl_ig + 130);

    auto ts_xxxyzz_yyyz = pbuffer.data(idx_ovl_ig + 131);

    auto ts_xxxyzz_yyzz = pbuffer.data(idx_ovl_ig + 132);

    auto ts_xxxyzz_yzzz = pbuffer.data(idx_ovl_ig + 133);

    auto ts_xxxyzz_zzzz = pbuffer.data(idx_ovl_ig + 134);

    auto ts_xxxzzz_xxxx = pbuffer.data(idx_ovl_ig + 135);

    auto ts_xxxzzz_xxxy = pbuffer.data(idx_ovl_ig + 136);

    auto ts_xxxzzz_xxxz = pbuffer.data(idx_ovl_ig + 137);

    auto ts_xxxzzz_xxyy = pbuffer.data(idx_ovl_ig + 138);

    auto ts_xxxzzz_xxyz = pbuffer.data(idx_ovl_ig + 139);

    auto ts_xxxzzz_xxzz = pbuffer.data(idx_ovl_ig + 140);

    auto ts_xxxzzz_xyyy = pbuffer.data(idx_ovl_ig + 141);

    auto ts_xxxzzz_xyyz = pbuffer.data(idx_ovl_ig + 142);

    auto ts_xxxzzz_xyzz = pbuffer.data(idx_ovl_ig + 143);

    auto ts_xxxzzz_xzzz = pbuffer.data(idx_ovl_ig + 144);

    auto ts_xxxzzz_yyyy = pbuffer.data(idx_ovl_ig + 145);

    auto ts_xxxzzz_yyyz = pbuffer.data(idx_ovl_ig + 146);

    auto ts_xxxzzz_yyzz = pbuffer.data(idx_ovl_ig + 147);

    auto ts_xxxzzz_yzzz = pbuffer.data(idx_ovl_ig + 148);

    auto ts_xxxzzz_zzzz = pbuffer.data(idx_ovl_ig + 149);

    auto ts_xxyyyy_xxxx = pbuffer.data(idx_ovl_ig + 150);

    auto ts_xxyyyy_xxxy = pbuffer.data(idx_ovl_ig + 151);

    auto ts_xxyyyy_xxxz = pbuffer.data(idx_ovl_ig + 152);

    auto ts_xxyyyy_xxyy = pbuffer.data(idx_ovl_ig + 153);

    auto ts_xxyyyy_xxyz = pbuffer.data(idx_ovl_ig + 154);

    auto ts_xxyyyy_xxzz = pbuffer.data(idx_ovl_ig + 155);

    auto ts_xxyyyy_xyyy = pbuffer.data(idx_ovl_ig + 156);

    auto ts_xxyyyy_xyyz = pbuffer.data(idx_ovl_ig + 157);

    auto ts_xxyyyy_xyzz = pbuffer.data(idx_ovl_ig + 158);

    auto ts_xxyyyy_xzzz = pbuffer.data(idx_ovl_ig + 159);

    auto ts_xxyyyy_yyyy = pbuffer.data(idx_ovl_ig + 160);

    auto ts_xxyyyy_yyyz = pbuffer.data(idx_ovl_ig + 161);

    auto ts_xxyyyy_yyzz = pbuffer.data(idx_ovl_ig + 162);

    auto ts_xxyyyy_yzzz = pbuffer.data(idx_ovl_ig + 163);

    auto ts_xxyyyy_zzzz = pbuffer.data(idx_ovl_ig + 164);

    auto ts_xxyyyz_xxxx = pbuffer.data(idx_ovl_ig + 165);

    auto ts_xxyyyz_xxxy = pbuffer.data(idx_ovl_ig + 166);

    auto ts_xxyyyz_xxxz = pbuffer.data(idx_ovl_ig + 167);

    auto ts_xxyyyz_xxyy = pbuffer.data(idx_ovl_ig + 168);

    auto ts_xxyyyz_xxyz = pbuffer.data(idx_ovl_ig + 169);

    auto ts_xxyyyz_xxzz = pbuffer.data(idx_ovl_ig + 170);

    auto ts_xxyyyz_xyyy = pbuffer.data(idx_ovl_ig + 171);

    auto ts_xxyyyz_xyyz = pbuffer.data(idx_ovl_ig + 172);

    auto ts_xxyyyz_xyzz = pbuffer.data(idx_ovl_ig + 173);

    auto ts_xxyyyz_xzzz = pbuffer.data(idx_ovl_ig + 174);

    auto ts_xxyyyz_yyyy = pbuffer.data(idx_ovl_ig + 175);

    auto ts_xxyyyz_yyyz = pbuffer.data(idx_ovl_ig + 176);

    auto ts_xxyyyz_yyzz = pbuffer.data(idx_ovl_ig + 177);

    auto ts_xxyyyz_yzzz = pbuffer.data(idx_ovl_ig + 178);

    auto ts_xxyyyz_zzzz = pbuffer.data(idx_ovl_ig + 179);

    auto ts_xxyyzz_xxxx = pbuffer.data(idx_ovl_ig + 180);

    auto ts_xxyyzz_xxxy = pbuffer.data(idx_ovl_ig + 181);

    auto ts_xxyyzz_xxxz = pbuffer.data(idx_ovl_ig + 182);

    auto ts_xxyyzz_xxyy = pbuffer.data(idx_ovl_ig + 183);

    auto ts_xxyyzz_xxyz = pbuffer.data(idx_ovl_ig + 184);

    auto ts_xxyyzz_xxzz = pbuffer.data(idx_ovl_ig + 185);

    auto ts_xxyyzz_xyyy = pbuffer.data(idx_ovl_ig + 186);

    auto ts_xxyyzz_xyyz = pbuffer.data(idx_ovl_ig + 187);

    auto ts_xxyyzz_xyzz = pbuffer.data(idx_ovl_ig + 188);

    auto ts_xxyyzz_xzzz = pbuffer.data(idx_ovl_ig + 189);

    auto ts_xxyyzz_yyyy = pbuffer.data(idx_ovl_ig + 190);

    auto ts_xxyyzz_yyyz = pbuffer.data(idx_ovl_ig + 191);

    auto ts_xxyyzz_yyzz = pbuffer.data(idx_ovl_ig + 192);

    auto ts_xxyyzz_yzzz = pbuffer.data(idx_ovl_ig + 193);

    auto ts_xxyyzz_zzzz = pbuffer.data(idx_ovl_ig + 194);

    auto ts_xxyzzz_xxxx = pbuffer.data(idx_ovl_ig + 195);

    auto ts_xxyzzz_xxxy = pbuffer.data(idx_ovl_ig + 196);

    auto ts_xxyzzz_xxxz = pbuffer.data(idx_ovl_ig + 197);

    auto ts_xxyzzz_xxyy = pbuffer.data(idx_ovl_ig + 198);

    auto ts_xxyzzz_xxyz = pbuffer.data(idx_ovl_ig + 199);

    auto ts_xxyzzz_xxzz = pbuffer.data(idx_ovl_ig + 200);

    auto ts_xxyzzz_xyyy = pbuffer.data(idx_ovl_ig + 201);

    auto ts_xxyzzz_xyyz = pbuffer.data(idx_ovl_ig + 202);

    auto ts_xxyzzz_xyzz = pbuffer.data(idx_ovl_ig + 203);

    auto ts_xxyzzz_xzzz = pbuffer.data(idx_ovl_ig + 204);

    auto ts_xxyzzz_yyyy = pbuffer.data(idx_ovl_ig + 205);

    auto ts_xxyzzz_yyyz = pbuffer.data(idx_ovl_ig + 206);

    auto ts_xxyzzz_yyzz = pbuffer.data(idx_ovl_ig + 207);

    auto ts_xxyzzz_yzzz = pbuffer.data(idx_ovl_ig + 208);

    auto ts_xxyzzz_zzzz = pbuffer.data(idx_ovl_ig + 209);

    auto ts_xxzzzz_xxxx = pbuffer.data(idx_ovl_ig + 210);

    auto ts_xxzzzz_xxxy = pbuffer.data(idx_ovl_ig + 211);

    auto ts_xxzzzz_xxxz = pbuffer.data(idx_ovl_ig + 212);

    auto ts_xxzzzz_xxyy = pbuffer.data(idx_ovl_ig + 213);

    auto ts_xxzzzz_xxyz = pbuffer.data(idx_ovl_ig + 214);

    auto ts_xxzzzz_xxzz = pbuffer.data(idx_ovl_ig + 215);

    auto ts_xxzzzz_xyyy = pbuffer.data(idx_ovl_ig + 216);

    auto ts_xxzzzz_xyyz = pbuffer.data(idx_ovl_ig + 217);

    auto ts_xxzzzz_xyzz = pbuffer.data(idx_ovl_ig + 218);

    auto ts_xxzzzz_xzzz = pbuffer.data(idx_ovl_ig + 219);

    auto ts_xxzzzz_yyyy = pbuffer.data(idx_ovl_ig + 220);

    auto ts_xxzzzz_yyyz = pbuffer.data(idx_ovl_ig + 221);

    auto ts_xxzzzz_yyzz = pbuffer.data(idx_ovl_ig + 222);

    auto ts_xxzzzz_yzzz = pbuffer.data(idx_ovl_ig + 223);

    auto ts_xxzzzz_zzzz = pbuffer.data(idx_ovl_ig + 224);

    auto ts_xyyyyy_xxxx = pbuffer.data(idx_ovl_ig + 225);

    auto ts_xyyyyy_xxxy = pbuffer.data(idx_ovl_ig + 226);

    auto ts_xyyyyy_xxxz = pbuffer.data(idx_ovl_ig + 227);

    auto ts_xyyyyy_xxyy = pbuffer.data(idx_ovl_ig + 228);

    auto ts_xyyyyy_xxyz = pbuffer.data(idx_ovl_ig + 229);

    auto ts_xyyyyy_xxzz = pbuffer.data(idx_ovl_ig + 230);

    auto ts_xyyyyy_xyyy = pbuffer.data(idx_ovl_ig + 231);

    auto ts_xyyyyy_xyyz = pbuffer.data(idx_ovl_ig + 232);

    auto ts_xyyyyy_xyzz = pbuffer.data(idx_ovl_ig + 233);

    auto ts_xyyyyy_xzzz = pbuffer.data(idx_ovl_ig + 234);

    auto ts_xyyyyy_yyyy = pbuffer.data(idx_ovl_ig + 235);

    auto ts_xyyyyy_yyyz = pbuffer.data(idx_ovl_ig + 236);

    auto ts_xyyyyy_yyzz = pbuffer.data(idx_ovl_ig + 237);

    auto ts_xyyyyy_yzzz = pbuffer.data(idx_ovl_ig + 238);

    auto ts_xyyyyy_zzzz = pbuffer.data(idx_ovl_ig + 239);

    auto ts_xyyyyz_xxxx = pbuffer.data(idx_ovl_ig + 240);

    auto ts_xyyyyz_xxxy = pbuffer.data(idx_ovl_ig + 241);

    auto ts_xyyyyz_xxxz = pbuffer.data(idx_ovl_ig + 242);

    auto ts_xyyyyz_xxyy = pbuffer.data(idx_ovl_ig + 243);

    auto ts_xyyyyz_xxyz = pbuffer.data(idx_ovl_ig + 244);

    auto ts_xyyyyz_xxzz = pbuffer.data(idx_ovl_ig + 245);

    auto ts_xyyyyz_xyyy = pbuffer.data(idx_ovl_ig + 246);

    auto ts_xyyyyz_xyyz = pbuffer.data(idx_ovl_ig + 247);

    auto ts_xyyyyz_xyzz = pbuffer.data(idx_ovl_ig + 248);

    auto ts_xyyyyz_xzzz = pbuffer.data(idx_ovl_ig + 249);

    auto ts_xyyyyz_yyyy = pbuffer.data(idx_ovl_ig + 250);

    auto ts_xyyyyz_yyyz = pbuffer.data(idx_ovl_ig + 251);

    auto ts_xyyyyz_yyzz = pbuffer.data(idx_ovl_ig + 252);

    auto ts_xyyyyz_yzzz = pbuffer.data(idx_ovl_ig + 253);

    auto ts_xyyyyz_zzzz = pbuffer.data(idx_ovl_ig + 254);

    auto ts_xyyyzz_xxxx = pbuffer.data(idx_ovl_ig + 255);

    auto ts_xyyyzz_xxxy = pbuffer.data(idx_ovl_ig + 256);

    auto ts_xyyyzz_xxxz = pbuffer.data(idx_ovl_ig + 257);

    auto ts_xyyyzz_xxyy = pbuffer.data(idx_ovl_ig + 258);

    auto ts_xyyyzz_xxyz = pbuffer.data(idx_ovl_ig + 259);

    auto ts_xyyyzz_xxzz = pbuffer.data(idx_ovl_ig + 260);

    auto ts_xyyyzz_xyyy = pbuffer.data(idx_ovl_ig + 261);

    auto ts_xyyyzz_xyyz = pbuffer.data(idx_ovl_ig + 262);

    auto ts_xyyyzz_xyzz = pbuffer.data(idx_ovl_ig + 263);

    auto ts_xyyyzz_xzzz = pbuffer.data(idx_ovl_ig + 264);

    auto ts_xyyyzz_yyyy = pbuffer.data(idx_ovl_ig + 265);

    auto ts_xyyyzz_yyyz = pbuffer.data(idx_ovl_ig + 266);

    auto ts_xyyyzz_yyzz = pbuffer.data(idx_ovl_ig + 267);

    auto ts_xyyyzz_yzzz = pbuffer.data(idx_ovl_ig + 268);

    auto ts_xyyyzz_zzzz = pbuffer.data(idx_ovl_ig + 269);

    auto ts_xyyzzz_xxxx = pbuffer.data(idx_ovl_ig + 270);

    auto ts_xyyzzz_xxxy = pbuffer.data(idx_ovl_ig + 271);

    auto ts_xyyzzz_xxxz = pbuffer.data(idx_ovl_ig + 272);

    auto ts_xyyzzz_xxyy = pbuffer.data(idx_ovl_ig + 273);

    auto ts_xyyzzz_xxyz = pbuffer.data(idx_ovl_ig + 274);

    auto ts_xyyzzz_xxzz = pbuffer.data(idx_ovl_ig + 275);

    auto ts_xyyzzz_xyyy = pbuffer.data(idx_ovl_ig + 276);

    auto ts_xyyzzz_xyyz = pbuffer.data(idx_ovl_ig + 277);

    auto ts_xyyzzz_xyzz = pbuffer.data(idx_ovl_ig + 278);

    auto ts_xyyzzz_xzzz = pbuffer.data(idx_ovl_ig + 279);

    auto ts_xyyzzz_yyyy = pbuffer.data(idx_ovl_ig + 280);

    auto ts_xyyzzz_yyyz = pbuffer.data(idx_ovl_ig + 281);

    auto ts_xyyzzz_yyzz = pbuffer.data(idx_ovl_ig + 282);

    auto ts_xyyzzz_yzzz = pbuffer.data(idx_ovl_ig + 283);

    auto ts_xyyzzz_zzzz = pbuffer.data(idx_ovl_ig + 284);

    auto ts_xyzzzz_xxxx = pbuffer.data(idx_ovl_ig + 285);

    auto ts_xyzzzz_xxxy = pbuffer.data(idx_ovl_ig + 286);

    auto ts_xyzzzz_xxxz = pbuffer.data(idx_ovl_ig + 287);

    auto ts_xyzzzz_xxyy = pbuffer.data(idx_ovl_ig + 288);

    auto ts_xyzzzz_xxyz = pbuffer.data(idx_ovl_ig + 289);

    auto ts_xyzzzz_xxzz = pbuffer.data(idx_ovl_ig + 290);

    auto ts_xyzzzz_xyyy = pbuffer.data(idx_ovl_ig + 291);

    auto ts_xyzzzz_xyyz = pbuffer.data(idx_ovl_ig + 292);

    auto ts_xyzzzz_xyzz = pbuffer.data(idx_ovl_ig + 293);

    auto ts_xyzzzz_xzzz = pbuffer.data(idx_ovl_ig + 294);

    auto ts_xyzzzz_yyyy = pbuffer.data(idx_ovl_ig + 295);

    auto ts_xyzzzz_yyyz = pbuffer.data(idx_ovl_ig + 296);

    auto ts_xyzzzz_yyzz = pbuffer.data(idx_ovl_ig + 297);

    auto ts_xyzzzz_yzzz = pbuffer.data(idx_ovl_ig + 298);

    auto ts_xyzzzz_zzzz = pbuffer.data(idx_ovl_ig + 299);

    auto ts_xzzzzz_xxxx = pbuffer.data(idx_ovl_ig + 300);

    auto ts_xzzzzz_xxxy = pbuffer.data(idx_ovl_ig + 301);

    auto ts_xzzzzz_xxxz = pbuffer.data(idx_ovl_ig + 302);

    auto ts_xzzzzz_xxyy = pbuffer.data(idx_ovl_ig + 303);

    auto ts_xzzzzz_xxyz = pbuffer.data(idx_ovl_ig + 304);

    auto ts_xzzzzz_xxzz = pbuffer.data(idx_ovl_ig + 305);

    auto ts_xzzzzz_xyyy = pbuffer.data(idx_ovl_ig + 306);

    auto ts_xzzzzz_xyyz = pbuffer.data(idx_ovl_ig + 307);

    auto ts_xzzzzz_xyzz = pbuffer.data(idx_ovl_ig + 308);

    auto ts_xzzzzz_xzzz = pbuffer.data(idx_ovl_ig + 309);

    auto ts_xzzzzz_yyyy = pbuffer.data(idx_ovl_ig + 310);

    auto ts_xzzzzz_yyyz = pbuffer.data(idx_ovl_ig + 311);

    auto ts_xzzzzz_yyzz = pbuffer.data(idx_ovl_ig + 312);

    auto ts_xzzzzz_yzzz = pbuffer.data(idx_ovl_ig + 313);

    auto ts_xzzzzz_zzzz = pbuffer.data(idx_ovl_ig + 314);

    auto ts_yyyyyy_xxxx = pbuffer.data(idx_ovl_ig + 315);

    auto ts_yyyyyy_xxxy = pbuffer.data(idx_ovl_ig + 316);

    auto ts_yyyyyy_xxxz = pbuffer.data(idx_ovl_ig + 317);

    auto ts_yyyyyy_xxyy = pbuffer.data(idx_ovl_ig + 318);

    auto ts_yyyyyy_xxyz = pbuffer.data(idx_ovl_ig + 319);

    auto ts_yyyyyy_xxzz = pbuffer.data(idx_ovl_ig + 320);

    auto ts_yyyyyy_xyyy = pbuffer.data(idx_ovl_ig + 321);

    auto ts_yyyyyy_xyyz = pbuffer.data(idx_ovl_ig + 322);

    auto ts_yyyyyy_xyzz = pbuffer.data(idx_ovl_ig + 323);

    auto ts_yyyyyy_xzzz = pbuffer.data(idx_ovl_ig + 324);

    auto ts_yyyyyy_yyyy = pbuffer.data(idx_ovl_ig + 325);

    auto ts_yyyyyy_yyyz = pbuffer.data(idx_ovl_ig + 326);

    auto ts_yyyyyy_yyzz = pbuffer.data(idx_ovl_ig + 327);

    auto ts_yyyyyy_yzzz = pbuffer.data(idx_ovl_ig + 328);

    auto ts_yyyyyy_zzzz = pbuffer.data(idx_ovl_ig + 329);

    auto ts_yyyyyz_xxxx = pbuffer.data(idx_ovl_ig + 330);

    auto ts_yyyyyz_xxxy = pbuffer.data(idx_ovl_ig + 331);

    auto ts_yyyyyz_xxxz = pbuffer.data(idx_ovl_ig + 332);

    auto ts_yyyyyz_xxyy = pbuffer.data(idx_ovl_ig + 333);

    auto ts_yyyyyz_xxyz = pbuffer.data(idx_ovl_ig + 334);

    auto ts_yyyyyz_xxzz = pbuffer.data(idx_ovl_ig + 335);

    auto ts_yyyyyz_xyyy = pbuffer.data(idx_ovl_ig + 336);

    auto ts_yyyyyz_xyyz = pbuffer.data(idx_ovl_ig + 337);

    auto ts_yyyyyz_xyzz = pbuffer.data(idx_ovl_ig + 338);

    auto ts_yyyyyz_xzzz = pbuffer.data(idx_ovl_ig + 339);

    auto ts_yyyyyz_yyyy = pbuffer.data(idx_ovl_ig + 340);

    auto ts_yyyyyz_yyyz = pbuffer.data(idx_ovl_ig + 341);

    auto ts_yyyyyz_yyzz = pbuffer.data(idx_ovl_ig + 342);

    auto ts_yyyyyz_yzzz = pbuffer.data(idx_ovl_ig + 343);

    auto ts_yyyyyz_zzzz = pbuffer.data(idx_ovl_ig + 344);

    auto ts_yyyyzz_xxxx = pbuffer.data(idx_ovl_ig + 345);

    auto ts_yyyyzz_xxxy = pbuffer.data(idx_ovl_ig + 346);

    auto ts_yyyyzz_xxxz = pbuffer.data(idx_ovl_ig + 347);

    auto ts_yyyyzz_xxyy = pbuffer.data(idx_ovl_ig + 348);

    auto ts_yyyyzz_xxyz = pbuffer.data(idx_ovl_ig + 349);

    auto ts_yyyyzz_xxzz = pbuffer.data(idx_ovl_ig + 350);

    auto ts_yyyyzz_xyyy = pbuffer.data(idx_ovl_ig + 351);

    auto ts_yyyyzz_xyyz = pbuffer.data(idx_ovl_ig + 352);

    auto ts_yyyyzz_xyzz = pbuffer.data(idx_ovl_ig + 353);

    auto ts_yyyyzz_xzzz = pbuffer.data(idx_ovl_ig + 354);

    auto ts_yyyyzz_yyyy = pbuffer.data(idx_ovl_ig + 355);

    auto ts_yyyyzz_yyyz = pbuffer.data(idx_ovl_ig + 356);

    auto ts_yyyyzz_yyzz = pbuffer.data(idx_ovl_ig + 357);

    auto ts_yyyyzz_yzzz = pbuffer.data(idx_ovl_ig + 358);

    auto ts_yyyyzz_zzzz = pbuffer.data(idx_ovl_ig + 359);

    auto ts_yyyzzz_xxxx = pbuffer.data(idx_ovl_ig + 360);

    auto ts_yyyzzz_xxxy = pbuffer.data(idx_ovl_ig + 361);

    auto ts_yyyzzz_xxxz = pbuffer.data(idx_ovl_ig + 362);

    auto ts_yyyzzz_xxyy = pbuffer.data(idx_ovl_ig + 363);

    auto ts_yyyzzz_xxyz = pbuffer.data(idx_ovl_ig + 364);

    auto ts_yyyzzz_xxzz = pbuffer.data(idx_ovl_ig + 365);

    auto ts_yyyzzz_xyyy = pbuffer.data(idx_ovl_ig + 366);

    auto ts_yyyzzz_xyyz = pbuffer.data(idx_ovl_ig + 367);

    auto ts_yyyzzz_xyzz = pbuffer.data(idx_ovl_ig + 368);

    auto ts_yyyzzz_xzzz = pbuffer.data(idx_ovl_ig + 369);

    auto ts_yyyzzz_yyyy = pbuffer.data(idx_ovl_ig + 370);

    auto ts_yyyzzz_yyyz = pbuffer.data(idx_ovl_ig + 371);

    auto ts_yyyzzz_yyzz = pbuffer.data(idx_ovl_ig + 372);

    auto ts_yyyzzz_yzzz = pbuffer.data(idx_ovl_ig + 373);

    auto ts_yyyzzz_zzzz = pbuffer.data(idx_ovl_ig + 374);

    auto ts_yyzzzz_xxxx = pbuffer.data(idx_ovl_ig + 375);

    auto ts_yyzzzz_xxxy = pbuffer.data(idx_ovl_ig + 376);

    auto ts_yyzzzz_xxxz = pbuffer.data(idx_ovl_ig + 377);

    auto ts_yyzzzz_xxyy = pbuffer.data(idx_ovl_ig + 378);

    auto ts_yyzzzz_xxyz = pbuffer.data(idx_ovl_ig + 379);

    auto ts_yyzzzz_xxzz = pbuffer.data(idx_ovl_ig + 380);

    auto ts_yyzzzz_xyyy = pbuffer.data(idx_ovl_ig + 381);

    auto ts_yyzzzz_xyyz = pbuffer.data(idx_ovl_ig + 382);

    auto ts_yyzzzz_xyzz = pbuffer.data(idx_ovl_ig + 383);

    auto ts_yyzzzz_xzzz = pbuffer.data(idx_ovl_ig + 384);

    auto ts_yyzzzz_yyyy = pbuffer.data(idx_ovl_ig + 385);

    auto ts_yyzzzz_yyyz = pbuffer.data(idx_ovl_ig + 386);

    auto ts_yyzzzz_yyzz = pbuffer.data(idx_ovl_ig + 387);

    auto ts_yyzzzz_yzzz = pbuffer.data(idx_ovl_ig + 388);

    auto ts_yyzzzz_zzzz = pbuffer.data(idx_ovl_ig + 389);

    auto ts_yzzzzz_xxxx = pbuffer.data(idx_ovl_ig + 390);

    auto ts_yzzzzz_xxxy = pbuffer.data(idx_ovl_ig + 391);

    auto ts_yzzzzz_xxxz = pbuffer.data(idx_ovl_ig + 392);

    auto ts_yzzzzz_xxyy = pbuffer.data(idx_ovl_ig + 393);

    auto ts_yzzzzz_xxyz = pbuffer.data(idx_ovl_ig + 394);

    auto ts_yzzzzz_xxzz = pbuffer.data(idx_ovl_ig + 395);

    auto ts_yzzzzz_xyyy = pbuffer.data(idx_ovl_ig + 396);

    auto ts_yzzzzz_xyyz = pbuffer.data(idx_ovl_ig + 397);

    auto ts_yzzzzz_xyzz = pbuffer.data(idx_ovl_ig + 398);

    auto ts_yzzzzz_xzzz = pbuffer.data(idx_ovl_ig + 399);

    auto ts_yzzzzz_yyyy = pbuffer.data(idx_ovl_ig + 400);

    auto ts_yzzzzz_yyyz = pbuffer.data(idx_ovl_ig + 401);

    auto ts_yzzzzz_yyzz = pbuffer.data(idx_ovl_ig + 402);

    auto ts_yzzzzz_yzzz = pbuffer.data(idx_ovl_ig + 403);

    auto ts_yzzzzz_zzzz = pbuffer.data(idx_ovl_ig + 404);

    auto ts_zzzzzz_xxxx = pbuffer.data(idx_ovl_ig + 405);

    auto ts_zzzzzz_xxxy = pbuffer.data(idx_ovl_ig + 406);

    auto ts_zzzzzz_xxxz = pbuffer.data(idx_ovl_ig + 407);

    auto ts_zzzzzz_xxyy = pbuffer.data(idx_ovl_ig + 408);

    auto ts_zzzzzz_xxyz = pbuffer.data(idx_ovl_ig + 409);

    auto ts_zzzzzz_xxzz = pbuffer.data(idx_ovl_ig + 410);

    auto ts_zzzzzz_xyyy = pbuffer.data(idx_ovl_ig + 411);

    auto ts_zzzzzz_xyyz = pbuffer.data(idx_ovl_ig + 412);

    auto ts_zzzzzz_xyzz = pbuffer.data(idx_ovl_ig + 413);

    auto ts_zzzzzz_xzzz = pbuffer.data(idx_ovl_ig + 414);

    auto ts_zzzzzz_yyyy = pbuffer.data(idx_ovl_ig + 415);

    auto ts_zzzzzz_yyyz = pbuffer.data(idx_ovl_ig + 416);

    auto ts_zzzzzz_yyzz = pbuffer.data(idx_ovl_ig + 417);

    auto ts_zzzzzz_yzzz = pbuffer.data(idx_ovl_ig + 418);

    auto ts_zzzzzz_zzzz = pbuffer.data(idx_ovl_ig + 419);

    // Set up 0-15 components of targeted buffer : IG

    auto tk_xxxxxx_xxxx = pbuffer.data(idx_kin_ig);

    auto tk_xxxxxx_xxxy = pbuffer.data(idx_kin_ig + 1);

    auto tk_xxxxxx_xxxz = pbuffer.data(idx_kin_ig + 2);

    auto tk_xxxxxx_xxyy = pbuffer.data(idx_kin_ig + 3);

    auto tk_xxxxxx_xxyz = pbuffer.data(idx_kin_ig + 4);

    auto tk_xxxxxx_xxzz = pbuffer.data(idx_kin_ig + 5);

    auto tk_xxxxxx_xyyy = pbuffer.data(idx_kin_ig + 6);

    auto tk_xxxxxx_xyyz = pbuffer.data(idx_kin_ig + 7);

    auto tk_xxxxxx_xyzz = pbuffer.data(idx_kin_ig + 8);

    auto tk_xxxxxx_xzzz = pbuffer.data(idx_kin_ig + 9);

    auto tk_xxxxxx_yyyy = pbuffer.data(idx_kin_ig + 10);

    auto tk_xxxxxx_yyyz = pbuffer.data(idx_kin_ig + 11);

    auto tk_xxxxxx_yyzz = pbuffer.data(idx_kin_ig + 12);

    auto tk_xxxxxx_yzzz = pbuffer.data(idx_kin_ig + 13);

    auto tk_xxxxxx_zzzz = pbuffer.data(idx_kin_ig + 14);

#pragma omp simd aligned(pa_x,               \
                             tk_xxxx_xxxx,   \
                             tk_xxxx_xxxy,   \
                             tk_xxxx_xxxz,   \
                             tk_xxxx_xxyy,   \
                             tk_xxxx_xxyz,   \
                             tk_xxxx_xxzz,   \
                             tk_xxxx_xyyy,   \
                             tk_xxxx_xyyz,   \
                             tk_xxxx_xyzz,   \
                             tk_xxxx_xzzz,   \
                             tk_xxxx_yyyy,   \
                             tk_xxxx_yyyz,   \
                             tk_xxxx_yyzz,   \
                             tk_xxxx_yzzz,   \
                             tk_xxxx_zzzz,   \
                             tk_xxxxx_xxx,   \
                             tk_xxxxx_xxxx,  \
                             tk_xxxxx_xxxy,  \
                             tk_xxxxx_xxxz,  \
                             tk_xxxxx_xxy,   \
                             tk_xxxxx_xxyy,  \
                             tk_xxxxx_xxyz,  \
                             tk_xxxxx_xxz,   \
                             tk_xxxxx_xxzz,  \
                             tk_xxxxx_xyy,   \
                             tk_xxxxx_xyyy,  \
                             tk_xxxxx_xyyz,  \
                             tk_xxxxx_xyz,   \
                             tk_xxxxx_xyzz,  \
                             tk_xxxxx_xzz,   \
                             tk_xxxxx_xzzz,  \
                             tk_xxxxx_yyy,   \
                             tk_xxxxx_yyyy,  \
                             tk_xxxxx_yyyz,  \
                             tk_xxxxx_yyz,   \
                             tk_xxxxx_yyzz,  \
                             tk_xxxxx_yzz,   \
                             tk_xxxxx_yzzz,  \
                             tk_xxxxx_zzz,   \
                             tk_xxxxx_zzzz,  \
                             tk_xxxxxx_xxxx, \
                             tk_xxxxxx_xxxy, \
                             tk_xxxxxx_xxxz, \
                             tk_xxxxxx_xxyy, \
                             tk_xxxxxx_xxyz, \
                             tk_xxxxxx_xxzz, \
                             tk_xxxxxx_xyyy, \
                             tk_xxxxxx_xyyz, \
                             tk_xxxxxx_xyzz, \
                             tk_xxxxxx_xzzz, \
                             tk_xxxxxx_yyyy, \
                             tk_xxxxxx_yyyz, \
                             tk_xxxxxx_yyzz, \
                             tk_xxxxxx_yzzz, \
                             tk_xxxxxx_zzzz, \
                             ts_xxxx_xxxx,   \
                             ts_xxxx_xxxy,   \
                             ts_xxxx_xxxz,   \
                             ts_xxxx_xxyy,   \
                             ts_xxxx_xxyz,   \
                             ts_xxxx_xxzz,   \
                             ts_xxxx_xyyy,   \
                             ts_xxxx_xyyz,   \
                             ts_xxxx_xyzz,   \
                             ts_xxxx_xzzz,   \
                             ts_xxxx_yyyy,   \
                             ts_xxxx_yyyz,   \
                             ts_xxxx_yyzz,   \
                             ts_xxxx_yzzz,   \
                             ts_xxxx_zzzz,   \
                             ts_xxxxxx_xxxx, \
                             ts_xxxxxx_xxxy, \
                             ts_xxxxxx_xxxz, \
                             ts_xxxxxx_xxyy, \
                             ts_xxxxxx_xxyz, \
                             ts_xxxxxx_xxzz, \
                             ts_xxxxxx_xyyy, \
                             ts_xxxxxx_xyyz, \
                             ts_xxxxxx_xyzz, \
                             ts_xxxxxx_xzzz, \
                             ts_xxxxxx_yyyy, \
                             ts_xxxxxx_yyyz, \
                             ts_xxxxxx_yyzz, \
                             ts_xxxxxx_yzzz, \
                             ts_xxxxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_xxxx[i] = -10.0 * ts_xxxx_xxxx[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxx[i] * fe_0 + 4.0 * tk_xxxxx_xxx[i] * fe_0 +
                            tk_xxxxx_xxxx[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxx[i] * fz_0;

        tk_xxxxxx_xxxy[i] = -10.0 * ts_xxxx_xxxy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxy[i] * fe_0 + 3.0 * tk_xxxxx_xxy[i] * fe_0 +
                            tk_xxxxx_xxxy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxy[i] * fz_0;

        tk_xxxxxx_xxxz[i] = -10.0 * ts_xxxx_xxxz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxz[i] * fe_0 + 3.0 * tk_xxxxx_xxz[i] * fe_0 +
                            tk_xxxxx_xxxz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxz[i] * fz_0;

        tk_xxxxxx_xxyy[i] = -10.0 * ts_xxxx_xxyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyy[i] * fe_0 + 2.0 * tk_xxxxx_xyy[i] * fe_0 +
                            tk_xxxxx_xxyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyy[i] * fz_0;

        tk_xxxxxx_xxyz[i] = -10.0 * ts_xxxx_xxyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyz[i] * fe_0 + 2.0 * tk_xxxxx_xyz[i] * fe_0 +
                            tk_xxxxx_xxyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyz[i] * fz_0;

        tk_xxxxxx_xxzz[i] = -10.0 * ts_xxxx_xxzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxzz[i] * fe_0 + 2.0 * tk_xxxxx_xzz[i] * fe_0 +
                            tk_xxxxx_xxzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxzz[i] * fz_0;

        tk_xxxxxx_xyyy[i] = -10.0 * ts_xxxx_xyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyy[i] * fe_0 + tk_xxxxx_yyy[i] * fe_0 +
                            tk_xxxxx_xyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyy[i] * fz_0;

        tk_xxxxxx_xyyz[i] = -10.0 * ts_xxxx_xyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyz[i] * fe_0 + tk_xxxxx_yyz[i] * fe_0 +
                            tk_xxxxx_xyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyz[i] * fz_0;

        tk_xxxxxx_xyzz[i] = -10.0 * ts_xxxx_xyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyzz[i] * fe_0 + tk_xxxxx_yzz[i] * fe_0 +
                            tk_xxxxx_xyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyzz[i] * fz_0;

        tk_xxxxxx_xzzz[i] = -10.0 * ts_xxxx_xzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xzzz[i] * fe_0 + tk_xxxxx_zzz[i] * fe_0 +
                            tk_xxxxx_xzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xzzz[i] * fz_0;

        tk_xxxxxx_yyyy[i] =
            -10.0 * ts_xxxx_yyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyy[i] * fe_0 + tk_xxxxx_yyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyyy[i] * fz_0;

        tk_xxxxxx_yyyz[i] =
            -10.0 * ts_xxxx_yyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyz[i] * fe_0 + tk_xxxxx_yyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyyz[i] * fz_0;

        tk_xxxxxx_yyzz[i] =
            -10.0 * ts_xxxx_yyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyzz[i] * fe_0 + tk_xxxxx_yyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyzz[i] * fz_0;

        tk_xxxxxx_yzzz[i] =
            -10.0 * ts_xxxx_yzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yzzz[i] * fe_0 + tk_xxxxx_yzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yzzz[i] * fz_0;

        tk_xxxxxx_zzzz[i] =
            -10.0 * ts_xxxx_zzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_zzzz[i] * fe_0 + tk_xxxxx_zzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_zzzz[i] * fz_0;
    }

    // Set up 15-30 components of targeted buffer : IG

    auto tk_xxxxxy_xxxx = pbuffer.data(idx_kin_ig + 15);

    auto tk_xxxxxy_xxxy = pbuffer.data(idx_kin_ig + 16);

    auto tk_xxxxxy_xxxz = pbuffer.data(idx_kin_ig + 17);

    auto tk_xxxxxy_xxyy = pbuffer.data(idx_kin_ig + 18);

    auto tk_xxxxxy_xxyz = pbuffer.data(idx_kin_ig + 19);

    auto tk_xxxxxy_xxzz = pbuffer.data(idx_kin_ig + 20);

    auto tk_xxxxxy_xyyy = pbuffer.data(idx_kin_ig + 21);

    auto tk_xxxxxy_xyyz = pbuffer.data(idx_kin_ig + 22);

    auto tk_xxxxxy_xyzz = pbuffer.data(idx_kin_ig + 23);

    auto tk_xxxxxy_xzzz = pbuffer.data(idx_kin_ig + 24);

    auto tk_xxxxxy_yyyy = pbuffer.data(idx_kin_ig + 25);

    auto tk_xxxxxy_yyyz = pbuffer.data(idx_kin_ig + 26);

    auto tk_xxxxxy_yyzz = pbuffer.data(idx_kin_ig + 27);

    auto tk_xxxxxy_yzzz = pbuffer.data(idx_kin_ig + 28);

    auto tk_xxxxxy_zzzz = pbuffer.data(idx_kin_ig + 29);

#pragma omp simd aligned(pa_y,               \
                             tk_xxxxx_xxx,   \
                             tk_xxxxx_xxxx,  \
                             tk_xxxxx_xxxy,  \
                             tk_xxxxx_xxxz,  \
                             tk_xxxxx_xxy,   \
                             tk_xxxxx_xxyy,  \
                             tk_xxxxx_xxyz,  \
                             tk_xxxxx_xxz,   \
                             tk_xxxxx_xxzz,  \
                             tk_xxxxx_xyy,   \
                             tk_xxxxx_xyyy,  \
                             tk_xxxxx_xyyz,  \
                             tk_xxxxx_xyz,   \
                             tk_xxxxx_xyzz,  \
                             tk_xxxxx_xzz,   \
                             tk_xxxxx_xzzz,  \
                             tk_xxxxx_yyy,   \
                             tk_xxxxx_yyyy,  \
                             tk_xxxxx_yyyz,  \
                             tk_xxxxx_yyz,   \
                             tk_xxxxx_yyzz,  \
                             tk_xxxxx_yzz,   \
                             tk_xxxxx_yzzz,  \
                             tk_xxxxx_zzz,   \
                             tk_xxxxx_zzzz,  \
                             tk_xxxxxy_xxxx, \
                             tk_xxxxxy_xxxy, \
                             tk_xxxxxy_xxxz, \
                             tk_xxxxxy_xxyy, \
                             tk_xxxxxy_xxyz, \
                             tk_xxxxxy_xxzz, \
                             tk_xxxxxy_xyyy, \
                             tk_xxxxxy_xyyz, \
                             tk_xxxxxy_xyzz, \
                             tk_xxxxxy_xzzz, \
                             tk_xxxxxy_yyyy, \
                             tk_xxxxxy_yyyz, \
                             tk_xxxxxy_yyzz, \
                             tk_xxxxxy_yzzz, \
                             tk_xxxxxy_zzzz, \
                             ts_xxxxxy_xxxx, \
                             ts_xxxxxy_xxxy, \
                             ts_xxxxxy_xxxz, \
                             ts_xxxxxy_xxyy, \
                             ts_xxxxxy_xxyz, \
                             ts_xxxxxy_xxzz, \
                             ts_xxxxxy_xyyy, \
                             ts_xxxxxy_xyyz, \
                             ts_xxxxxy_xyzz, \
                             ts_xxxxxy_xzzz, \
                             ts_xxxxxy_yyyy, \
                             ts_xxxxxy_yyyz, \
                             ts_xxxxxy_yyzz, \
                             ts_xxxxxy_yzzz, \
                             ts_xxxxxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxy_xxxx[i] = tk_xxxxx_xxxx[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxx[i] * fz_0;

        tk_xxxxxy_xxxy[i] = tk_xxxxx_xxx[i] * fe_0 + tk_xxxxx_xxxy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxy[i] * fz_0;

        tk_xxxxxy_xxxz[i] = tk_xxxxx_xxxz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxz[i] * fz_0;

        tk_xxxxxy_xxyy[i] = 2.0 * tk_xxxxx_xxy[i] * fe_0 + tk_xxxxx_xxyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyy[i] * fz_0;

        tk_xxxxxy_xxyz[i] = tk_xxxxx_xxz[i] * fe_0 + tk_xxxxx_xxyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyz[i] * fz_0;

        tk_xxxxxy_xxzz[i] = tk_xxxxx_xxzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxzz[i] * fz_0;

        tk_xxxxxy_xyyy[i] = 3.0 * tk_xxxxx_xyy[i] * fe_0 + tk_xxxxx_xyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyy[i] * fz_0;

        tk_xxxxxy_xyyz[i] = 2.0 * tk_xxxxx_xyz[i] * fe_0 + tk_xxxxx_xyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyz[i] * fz_0;

        tk_xxxxxy_xyzz[i] = tk_xxxxx_xzz[i] * fe_0 + tk_xxxxx_xyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyzz[i] * fz_0;

        tk_xxxxxy_xzzz[i] = tk_xxxxx_xzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xzzz[i] * fz_0;

        tk_xxxxxy_yyyy[i] = 4.0 * tk_xxxxx_yyy[i] * fe_0 + tk_xxxxx_yyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyy[i] * fz_0;

        tk_xxxxxy_yyyz[i] = 3.0 * tk_xxxxx_yyz[i] * fe_0 + tk_xxxxx_yyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyz[i] * fz_0;

        tk_xxxxxy_yyzz[i] = 2.0 * tk_xxxxx_yzz[i] * fe_0 + tk_xxxxx_yyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyzz[i] * fz_0;

        tk_xxxxxy_yzzz[i] = tk_xxxxx_zzz[i] * fe_0 + tk_xxxxx_yzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yzzz[i] * fz_0;

        tk_xxxxxy_zzzz[i] = tk_xxxxx_zzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_zzzz[i] * fz_0;
    }

    // Set up 30-45 components of targeted buffer : IG

    auto tk_xxxxxz_xxxx = pbuffer.data(idx_kin_ig + 30);

    auto tk_xxxxxz_xxxy = pbuffer.data(idx_kin_ig + 31);

    auto tk_xxxxxz_xxxz = pbuffer.data(idx_kin_ig + 32);

    auto tk_xxxxxz_xxyy = pbuffer.data(idx_kin_ig + 33);

    auto tk_xxxxxz_xxyz = pbuffer.data(idx_kin_ig + 34);

    auto tk_xxxxxz_xxzz = pbuffer.data(idx_kin_ig + 35);

    auto tk_xxxxxz_xyyy = pbuffer.data(idx_kin_ig + 36);

    auto tk_xxxxxz_xyyz = pbuffer.data(idx_kin_ig + 37);

    auto tk_xxxxxz_xyzz = pbuffer.data(idx_kin_ig + 38);

    auto tk_xxxxxz_xzzz = pbuffer.data(idx_kin_ig + 39);

    auto tk_xxxxxz_yyyy = pbuffer.data(idx_kin_ig + 40);

    auto tk_xxxxxz_yyyz = pbuffer.data(idx_kin_ig + 41);

    auto tk_xxxxxz_yyzz = pbuffer.data(idx_kin_ig + 42);

    auto tk_xxxxxz_yzzz = pbuffer.data(idx_kin_ig + 43);

    auto tk_xxxxxz_zzzz = pbuffer.data(idx_kin_ig + 44);

#pragma omp simd aligned(pa_z,               \
                             tk_xxxxx_xxx,   \
                             tk_xxxxx_xxxx,  \
                             tk_xxxxx_xxxy,  \
                             tk_xxxxx_xxxz,  \
                             tk_xxxxx_xxy,   \
                             tk_xxxxx_xxyy,  \
                             tk_xxxxx_xxyz,  \
                             tk_xxxxx_xxz,   \
                             tk_xxxxx_xxzz,  \
                             tk_xxxxx_xyy,   \
                             tk_xxxxx_xyyy,  \
                             tk_xxxxx_xyyz,  \
                             tk_xxxxx_xyz,   \
                             tk_xxxxx_xyzz,  \
                             tk_xxxxx_xzz,   \
                             tk_xxxxx_xzzz,  \
                             tk_xxxxx_yyy,   \
                             tk_xxxxx_yyyy,  \
                             tk_xxxxx_yyyz,  \
                             tk_xxxxx_yyz,   \
                             tk_xxxxx_yyzz,  \
                             tk_xxxxx_yzz,   \
                             tk_xxxxx_yzzz,  \
                             tk_xxxxx_zzz,   \
                             tk_xxxxx_zzzz,  \
                             tk_xxxxxz_xxxx, \
                             tk_xxxxxz_xxxy, \
                             tk_xxxxxz_xxxz, \
                             tk_xxxxxz_xxyy, \
                             tk_xxxxxz_xxyz, \
                             tk_xxxxxz_xxzz, \
                             tk_xxxxxz_xyyy, \
                             tk_xxxxxz_xyyz, \
                             tk_xxxxxz_xyzz, \
                             tk_xxxxxz_xzzz, \
                             tk_xxxxxz_yyyy, \
                             tk_xxxxxz_yyyz, \
                             tk_xxxxxz_yyzz, \
                             tk_xxxxxz_yzzz, \
                             tk_xxxxxz_zzzz, \
                             ts_xxxxxz_xxxx, \
                             ts_xxxxxz_xxxy, \
                             ts_xxxxxz_xxxz, \
                             ts_xxxxxz_xxyy, \
                             ts_xxxxxz_xxyz, \
                             ts_xxxxxz_xxzz, \
                             ts_xxxxxz_xyyy, \
                             ts_xxxxxz_xyyz, \
                             ts_xxxxxz_xyzz, \
                             ts_xxxxxz_xzzz, \
                             ts_xxxxxz_yyyy, \
                             ts_xxxxxz_yyyz, \
                             ts_xxxxxz_yyzz, \
                             ts_xxxxxz_yzzz, \
                             ts_xxxxxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxz_xxxx[i] = tk_xxxxx_xxxx[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxx[i] * fz_0;

        tk_xxxxxz_xxxy[i] = tk_xxxxx_xxxy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxy[i] * fz_0;

        tk_xxxxxz_xxxz[i] = tk_xxxxx_xxx[i] * fe_0 + tk_xxxxx_xxxz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxz[i] * fz_0;

        tk_xxxxxz_xxyy[i] = tk_xxxxx_xxyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyy[i] * fz_0;

        tk_xxxxxz_xxyz[i] = tk_xxxxx_xxy[i] * fe_0 + tk_xxxxx_xxyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyz[i] * fz_0;

        tk_xxxxxz_xxzz[i] = 2.0 * tk_xxxxx_xxz[i] * fe_0 + tk_xxxxx_xxzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxzz[i] * fz_0;

        tk_xxxxxz_xyyy[i] = tk_xxxxx_xyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyy[i] * fz_0;

        tk_xxxxxz_xyyz[i] = tk_xxxxx_xyy[i] * fe_0 + tk_xxxxx_xyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyz[i] * fz_0;

        tk_xxxxxz_xyzz[i] = 2.0 * tk_xxxxx_xyz[i] * fe_0 + tk_xxxxx_xyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyzz[i] * fz_0;

        tk_xxxxxz_xzzz[i] = 3.0 * tk_xxxxx_xzz[i] * fe_0 + tk_xxxxx_xzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xzzz[i] * fz_0;

        tk_xxxxxz_yyyy[i] = tk_xxxxx_yyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyy[i] * fz_0;

        tk_xxxxxz_yyyz[i] = tk_xxxxx_yyy[i] * fe_0 + tk_xxxxx_yyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyz[i] * fz_0;

        tk_xxxxxz_yyzz[i] = 2.0 * tk_xxxxx_yyz[i] * fe_0 + tk_xxxxx_yyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyzz[i] * fz_0;

        tk_xxxxxz_yzzz[i] = 3.0 * tk_xxxxx_yzz[i] * fe_0 + tk_xxxxx_yzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yzzz[i] * fz_0;

        tk_xxxxxz_zzzz[i] = 4.0 * tk_xxxxx_zzz[i] * fe_0 + tk_xxxxx_zzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_zzzz[i] * fz_0;
    }

    // Set up 45-60 components of targeted buffer : IG

    auto tk_xxxxyy_xxxx = pbuffer.data(idx_kin_ig + 45);

    auto tk_xxxxyy_xxxy = pbuffer.data(idx_kin_ig + 46);

    auto tk_xxxxyy_xxxz = pbuffer.data(idx_kin_ig + 47);

    auto tk_xxxxyy_xxyy = pbuffer.data(idx_kin_ig + 48);

    auto tk_xxxxyy_xxyz = pbuffer.data(idx_kin_ig + 49);

    auto tk_xxxxyy_xxzz = pbuffer.data(idx_kin_ig + 50);

    auto tk_xxxxyy_xyyy = pbuffer.data(idx_kin_ig + 51);

    auto tk_xxxxyy_xyyz = pbuffer.data(idx_kin_ig + 52);

    auto tk_xxxxyy_xyzz = pbuffer.data(idx_kin_ig + 53);

    auto tk_xxxxyy_xzzz = pbuffer.data(idx_kin_ig + 54);

    auto tk_xxxxyy_yyyy = pbuffer.data(idx_kin_ig + 55);

    auto tk_xxxxyy_yyyz = pbuffer.data(idx_kin_ig + 56);

    auto tk_xxxxyy_yyzz = pbuffer.data(idx_kin_ig + 57);

    auto tk_xxxxyy_yzzz = pbuffer.data(idx_kin_ig + 58);

    auto tk_xxxxyy_zzzz = pbuffer.data(idx_kin_ig + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xxxx_xxxx,   \
                             tk_xxxx_xxxz,   \
                             tk_xxxx_xxzz,   \
                             tk_xxxx_xzzz,   \
                             tk_xxxxy_xxxx,  \
                             tk_xxxxy_xxxz,  \
                             tk_xxxxy_xxzz,  \
                             tk_xxxxy_xzzz,  \
                             tk_xxxxyy_xxxx, \
                             tk_xxxxyy_xxxy, \
                             tk_xxxxyy_xxxz, \
                             tk_xxxxyy_xxyy, \
                             tk_xxxxyy_xxyz, \
                             tk_xxxxyy_xxzz, \
                             tk_xxxxyy_xyyy, \
                             tk_xxxxyy_xyyz, \
                             tk_xxxxyy_xyzz, \
                             tk_xxxxyy_xzzz, \
                             tk_xxxxyy_yyyy, \
                             tk_xxxxyy_yyyz, \
                             tk_xxxxyy_yyzz, \
                             tk_xxxxyy_yzzz, \
                             tk_xxxxyy_zzzz, \
                             tk_xxxyy_xxxy,  \
                             tk_xxxyy_xxy,   \
                             tk_xxxyy_xxyy,  \
                             tk_xxxyy_xxyz,  \
                             tk_xxxyy_xyy,   \
                             tk_xxxyy_xyyy,  \
                             tk_xxxyy_xyyz,  \
                             tk_xxxyy_xyz,   \
                             tk_xxxyy_xyzz,  \
                             tk_xxxyy_yyy,   \
                             tk_xxxyy_yyyy,  \
                             tk_xxxyy_yyyz,  \
                             tk_xxxyy_yyz,   \
                             tk_xxxyy_yyzz,  \
                             tk_xxxyy_yzz,   \
                             tk_xxxyy_yzzz,  \
                             tk_xxxyy_zzzz,  \
                             tk_xxyy_xxxy,   \
                             tk_xxyy_xxyy,   \
                             tk_xxyy_xxyz,   \
                             tk_xxyy_xyyy,   \
                             tk_xxyy_xyyz,   \
                             tk_xxyy_xyzz,   \
                             tk_xxyy_yyyy,   \
                             tk_xxyy_yyyz,   \
                             tk_xxyy_yyzz,   \
                             tk_xxyy_yzzz,   \
                             tk_xxyy_zzzz,   \
                             ts_xxxx_xxxx,   \
                             ts_xxxx_xxxz,   \
                             ts_xxxx_xxzz,   \
                             ts_xxxx_xzzz,   \
                             ts_xxxxyy_xxxx, \
                             ts_xxxxyy_xxxy, \
                             ts_xxxxyy_xxxz, \
                             ts_xxxxyy_xxyy, \
                             ts_xxxxyy_xxyz, \
                             ts_xxxxyy_xxzz, \
                             ts_xxxxyy_xyyy, \
                             ts_xxxxyy_xyyz, \
                             ts_xxxxyy_xyzz, \
                             ts_xxxxyy_xzzz, \
                             ts_xxxxyy_yyyy, \
                             ts_xxxxyy_yyyz, \
                             ts_xxxxyy_yyzz, \
                             ts_xxxxyy_yzzz, \
                             ts_xxxxyy_zzzz, \
                             ts_xxyy_xxxy,   \
                             ts_xxyy_xxyy,   \
                             ts_xxyy_xxyz,   \
                             ts_xxyy_xyyy,   \
                             ts_xxyy_xyyz,   \
                             ts_xxyy_xyzz,   \
                             ts_xxyy_yyyy,   \
                             ts_xxyy_yyyz,   \
                             ts_xxyy_yyzz,   \
                             ts_xxyy_yzzz,   \
                             ts_xxyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxyy_xxxx[i] =
            -2.0 * ts_xxxx_xxxx[i] * fbe_0 * fz_0 + tk_xxxx_xxxx[i] * fe_0 + tk_xxxxy_xxxx[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxx[i] * fz_0;

        tk_xxxxyy_xxxy[i] = -6.0 * ts_xxyy_xxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxy[i] * fe_0 + 3.0 * tk_xxxyy_xxy[i] * fe_0 +
                            tk_xxxyy_xxxy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxy[i] * fz_0;

        tk_xxxxyy_xxxz[i] =
            -2.0 * ts_xxxx_xxxz[i] * fbe_0 * fz_0 + tk_xxxx_xxxz[i] * fe_0 + tk_xxxxy_xxxz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxz[i] * fz_0;

        tk_xxxxyy_xxyy[i] = -6.0 * ts_xxyy_xxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyy[i] * fe_0 + 2.0 * tk_xxxyy_xyy[i] * fe_0 +
                            tk_xxxyy_xxyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyy[i] * fz_0;

        tk_xxxxyy_xxyz[i] = -6.0 * ts_xxyy_xxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyz[i] * fe_0 + 2.0 * tk_xxxyy_xyz[i] * fe_0 +
                            tk_xxxyy_xxyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyz[i] * fz_0;

        tk_xxxxyy_xxzz[i] =
            -2.0 * ts_xxxx_xxzz[i] * fbe_0 * fz_0 + tk_xxxx_xxzz[i] * fe_0 + tk_xxxxy_xxzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxzz[i] * fz_0;

        tk_xxxxyy_xyyy[i] = -6.0 * ts_xxyy_xyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyy[i] * fe_0 + tk_xxxyy_yyy[i] * fe_0 +
                            tk_xxxyy_xyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyy[i] * fz_0;

        tk_xxxxyy_xyyz[i] = -6.0 * ts_xxyy_xyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyz[i] * fe_0 + tk_xxxyy_yyz[i] * fe_0 +
                            tk_xxxyy_xyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyz[i] * fz_0;

        tk_xxxxyy_xyzz[i] = -6.0 * ts_xxyy_xyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyzz[i] * fe_0 + tk_xxxyy_yzz[i] * fe_0 +
                            tk_xxxyy_xyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyzz[i] * fz_0;

        tk_xxxxyy_xzzz[i] =
            -2.0 * ts_xxxx_xzzz[i] * fbe_0 * fz_0 + tk_xxxx_xzzz[i] * fe_0 + tk_xxxxy_xzzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xzzz[i] * fz_0;

        tk_xxxxyy_yyyy[i] =
            -6.0 * ts_xxyy_yyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyy[i] * fe_0 + tk_xxxyy_yyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyyy[i] * fz_0;

        tk_xxxxyy_yyyz[i] =
            -6.0 * ts_xxyy_yyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyz[i] * fe_0 + tk_xxxyy_yyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyyz[i] * fz_0;

        tk_xxxxyy_yyzz[i] =
            -6.0 * ts_xxyy_yyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyzz[i] * fe_0 + tk_xxxyy_yyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyzz[i] * fz_0;

        tk_xxxxyy_yzzz[i] =
            -6.0 * ts_xxyy_yzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yzzz[i] * fe_0 + tk_xxxyy_yzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yzzz[i] * fz_0;

        tk_xxxxyy_zzzz[i] =
            -6.0 * ts_xxyy_zzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_zzzz[i] * fe_0 + tk_xxxyy_zzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_zzzz[i] * fz_0;
    }

    // Set up 60-75 components of targeted buffer : IG

    auto tk_xxxxyz_xxxx = pbuffer.data(idx_kin_ig + 60);

    auto tk_xxxxyz_xxxy = pbuffer.data(idx_kin_ig + 61);

    auto tk_xxxxyz_xxxz = pbuffer.data(idx_kin_ig + 62);

    auto tk_xxxxyz_xxyy = pbuffer.data(idx_kin_ig + 63);

    auto tk_xxxxyz_xxyz = pbuffer.data(idx_kin_ig + 64);

    auto tk_xxxxyz_xxzz = pbuffer.data(idx_kin_ig + 65);

    auto tk_xxxxyz_xyyy = pbuffer.data(idx_kin_ig + 66);

    auto tk_xxxxyz_xyyz = pbuffer.data(idx_kin_ig + 67);

    auto tk_xxxxyz_xyzz = pbuffer.data(idx_kin_ig + 68);

    auto tk_xxxxyz_xzzz = pbuffer.data(idx_kin_ig + 69);

    auto tk_xxxxyz_yyyy = pbuffer.data(idx_kin_ig + 70);

    auto tk_xxxxyz_yyyz = pbuffer.data(idx_kin_ig + 71);

    auto tk_xxxxyz_yyzz = pbuffer.data(idx_kin_ig + 72);

    auto tk_xxxxyz_yzzz = pbuffer.data(idx_kin_ig + 73);

    auto tk_xxxxyz_zzzz = pbuffer.data(idx_kin_ig + 74);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_xxxxy_xxxy,  \
                             tk_xxxxy_xxyy,  \
                             tk_xxxxy_xyyy,  \
                             tk_xxxxy_yyyy,  \
                             tk_xxxxyz_xxxx, \
                             tk_xxxxyz_xxxy, \
                             tk_xxxxyz_xxxz, \
                             tk_xxxxyz_xxyy, \
                             tk_xxxxyz_xxyz, \
                             tk_xxxxyz_xxzz, \
                             tk_xxxxyz_xyyy, \
                             tk_xxxxyz_xyyz, \
                             tk_xxxxyz_xyzz, \
                             tk_xxxxyz_xzzz, \
                             tk_xxxxyz_yyyy, \
                             tk_xxxxyz_yyyz, \
                             tk_xxxxyz_yyzz, \
                             tk_xxxxyz_yzzz, \
                             tk_xxxxyz_zzzz, \
                             tk_xxxxz_xxxx,  \
                             tk_xxxxz_xxxz,  \
                             tk_xxxxz_xxyz,  \
                             tk_xxxxz_xxz,   \
                             tk_xxxxz_xxzz,  \
                             tk_xxxxz_xyyz,  \
                             tk_xxxxz_xyz,   \
                             tk_xxxxz_xyzz,  \
                             tk_xxxxz_xzz,   \
                             tk_xxxxz_xzzz,  \
                             tk_xxxxz_yyyz,  \
                             tk_xxxxz_yyz,   \
                             tk_xxxxz_yyzz,  \
                             tk_xxxxz_yzz,   \
                             tk_xxxxz_yzzz,  \
                             tk_xxxxz_zzz,   \
                             tk_xxxxz_zzzz,  \
                             ts_xxxxyz_xxxx, \
                             ts_xxxxyz_xxxy, \
                             ts_xxxxyz_xxxz, \
                             ts_xxxxyz_xxyy, \
                             ts_xxxxyz_xxyz, \
                             ts_xxxxyz_xxzz, \
                             ts_xxxxyz_xyyy, \
                             ts_xxxxyz_xyyz, \
                             ts_xxxxyz_xyzz, \
                             ts_xxxxyz_xzzz, \
                             ts_xxxxyz_yyyy, \
                             ts_xxxxyz_yyyz, \
                             ts_xxxxyz_yyzz, \
                             ts_xxxxyz_yzzz, \
                             ts_xxxxyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxyz_xxxx[i] = tk_xxxxz_xxxx[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxx[i] * fz_0;

        tk_xxxxyz_xxxy[i] = tk_xxxxy_xxxy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxxy[i] * fz_0;

        tk_xxxxyz_xxxz[i] = tk_xxxxz_xxxz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxz[i] * fz_0;

        tk_xxxxyz_xxyy[i] = tk_xxxxy_xxyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxyy[i] * fz_0;

        tk_xxxxyz_xxyz[i] = tk_xxxxz_xxz[i] * fe_0 + tk_xxxxz_xxyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxyz[i] * fz_0;

        tk_xxxxyz_xxzz[i] = tk_xxxxz_xxzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxzz[i] * fz_0;

        tk_xxxxyz_xyyy[i] = tk_xxxxy_xyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xyyy[i] * fz_0;

        tk_xxxxyz_xyyz[i] = 2.0 * tk_xxxxz_xyz[i] * fe_0 + tk_xxxxz_xyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyyz[i] * fz_0;

        tk_xxxxyz_xyzz[i] = tk_xxxxz_xzz[i] * fe_0 + tk_xxxxz_xyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyzz[i] * fz_0;

        tk_xxxxyz_xzzz[i] = tk_xxxxz_xzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xzzz[i] * fz_0;

        tk_xxxxyz_yyyy[i] = tk_xxxxy_yyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_yyyy[i] * fz_0;

        tk_xxxxyz_yyyz[i] = 3.0 * tk_xxxxz_yyz[i] * fe_0 + tk_xxxxz_yyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyyz[i] * fz_0;

        tk_xxxxyz_yyzz[i] = 2.0 * tk_xxxxz_yzz[i] * fe_0 + tk_xxxxz_yyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyzz[i] * fz_0;

        tk_xxxxyz_yzzz[i] = tk_xxxxz_zzz[i] * fe_0 + tk_xxxxz_yzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yzzz[i] * fz_0;

        tk_xxxxyz_zzzz[i] = tk_xxxxz_zzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_zzzz[i] * fz_0;
    }

    // Set up 75-90 components of targeted buffer : IG

    auto tk_xxxxzz_xxxx = pbuffer.data(idx_kin_ig + 75);

    auto tk_xxxxzz_xxxy = pbuffer.data(idx_kin_ig + 76);

    auto tk_xxxxzz_xxxz = pbuffer.data(idx_kin_ig + 77);

    auto tk_xxxxzz_xxyy = pbuffer.data(idx_kin_ig + 78);

    auto tk_xxxxzz_xxyz = pbuffer.data(idx_kin_ig + 79);

    auto tk_xxxxzz_xxzz = pbuffer.data(idx_kin_ig + 80);

    auto tk_xxxxzz_xyyy = pbuffer.data(idx_kin_ig + 81);

    auto tk_xxxxzz_xyyz = pbuffer.data(idx_kin_ig + 82);

    auto tk_xxxxzz_xyzz = pbuffer.data(idx_kin_ig + 83);

    auto tk_xxxxzz_xzzz = pbuffer.data(idx_kin_ig + 84);

    auto tk_xxxxzz_yyyy = pbuffer.data(idx_kin_ig + 85);

    auto tk_xxxxzz_yyyz = pbuffer.data(idx_kin_ig + 86);

    auto tk_xxxxzz_yyzz = pbuffer.data(idx_kin_ig + 87);

    auto tk_xxxxzz_yzzz = pbuffer.data(idx_kin_ig + 88);

    auto tk_xxxxzz_zzzz = pbuffer.data(idx_kin_ig + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xxxx_xxxx,   \
                             tk_xxxx_xxxy,   \
                             tk_xxxx_xxyy,   \
                             tk_xxxx_xyyy,   \
                             tk_xxxxz_xxxx,  \
                             tk_xxxxz_xxxy,  \
                             tk_xxxxz_xxyy,  \
                             tk_xxxxz_xyyy,  \
                             tk_xxxxzz_xxxx, \
                             tk_xxxxzz_xxxy, \
                             tk_xxxxzz_xxxz, \
                             tk_xxxxzz_xxyy, \
                             tk_xxxxzz_xxyz, \
                             tk_xxxxzz_xxzz, \
                             tk_xxxxzz_xyyy, \
                             tk_xxxxzz_xyyz, \
                             tk_xxxxzz_xyzz, \
                             tk_xxxxzz_xzzz, \
                             tk_xxxxzz_yyyy, \
                             tk_xxxxzz_yyyz, \
                             tk_xxxxzz_yyzz, \
                             tk_xxxxzz_yzzz, \
                             tk_xxxxzz_zzzz, \
                             tk_xxxzz_xxxz,  \
                             tk_xxxzz_xxyz,  \
                             tk_xxxzz_xxz,   \
                             tk_xxxzz_xxzz,  \
                             tk_xxxzz_xyyz,  \
                             tk_xxxzz_xyz,   \
                             tk_xxxzz_xyzz,  \
                             tk_xxxzz_xzz,   \
                             tk_xxxzz_xzzz,  \
                             tk_xxxzz_yyyy,  \
                             tk_xxxzz_yyyz,  \
                             tk_xxxzz_yyz,   \
                             tk_xxxzz_yyzz,  \
                             tk_xxxzz_yzz,   \
                             tk_xxxzz_yzzz,  \
                             tk_xxxzz_zzz,   \
                             tk_xxxzz_zzzz,  \
                             tk_xxzz_xxxz,   \
                             tk_xxzz_xxyz,   \
                             tk_xxzz_xxzz,   \
                             tk_xxzz_xyyz,   \
                             tk_xxzz_xyzz,   \
                             tk_xxzz_xzzz,   \
                             tk_xxzz_yyyy,   \
                             tk_xxzz_yyyz,   \
                             tk_xxzz_yyzz,   \
                             tk_xxzz_yzzz,   \
                             tk_xxzz_zzzz,   \
                             ts_xxxx_xxxx,   \
                             ts_xxxx_xxxy,   \
                             ts_xxxx_xxyy,   \
                             ts_xxxx_xyyy,   \
                             ts_xxxxzz_xxxx, \
                             ts_xxxxzz_xxxy, \
                             ts_xxxxzz_xxxz, \
                             ts_xxxxzz_xxyy, \
                             ts_xxxxzz_xxyz, \
                             ts_xxxxzz_xxzz, \
                             ts_xxxxzz_xyyy, \
                             ts_xxxxzz_xyyz, \
                             ts_xxxxzz_xyzz, \
                             ts_xxxxzz_xzzz, \
                             ts_xxxxzz_yyyy, \
                             ts_xxxxzz_yyyz, \
                             ts_xxxxzz_yyzz, \
                             ts_xxxxzz_yzzz, \
                             ts_xxxxzz_zzzz, \
                             ts_xxzz_xxxz,   \
                             ts_xxzz_xxyz,   \
                             ts_xxzz_xxzz,   \
                             ts_xxzz_xyyz,   \
                             ts_xxzz_xyzz,   \
                             ts_xxzz_xzzz,   \
                             ts_xxzz_yyyy,   \
                             ts_xxzz_yyyz,   \
                             ts_xxzz_yyzz,   \
                             ts_xxzz_yzzz,   \
                             ts_xxzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxzz_xxxx[i] =
            -2.0 * ts_xxxx_xxxx[i] * fbe_0 * fz_0 + tk_xxxx_xxxx[i] * fe_0 + tk_xxxxz_xxxx[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxx[i] * fz_0;

        tk_xxxxzz_xxxy[i] =
            -2.0 * ts_xxxx_xxxy[i] * fbe_0 * fz_0 + tk_xxxx_xxxy[i] * fe_0 + tk_xxxxz_xxxy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxy[i] * fz_0;

        tk_xxxxzz_xxxz[i] = -6.0 * ts_xxzz_xxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxz[i] * fe_0 + 3.0 * tk_xxxzz_xxz[i] * fe_0 +
                            tk_xxxzz_xxxz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxz[i] * fz_0;

        tk_xxxxzz_xxyy[i] =
            -2.0 * ts_xxxx_xxyy[i] * fbe_0 * fz_0 + tk_xxxx_xxyy[i] * fe_0 + tk_xxxxz_xxyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxyy[i] * fz_0;

        tk_xxxxzz_xxyz[i] = -6.0 * ts_xxzz_xxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyz[i] * fe_0 + 2.0 * tk_xxxzz_xyz[i] * fe_0 +
                            tk_xxxzz_xxyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxyz[i] * fz_0;

        tk_xxxxzz_xxzz[i] = -6.0 * ts_xxzz_xxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxzz[i] * fe_0 + 2.0 * tk_xxxzz_xzz[i] * fe_0 +
                            tk_xxxzz_xxzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxzz[i] * fz_0;

        tk_xxxxzz_xyyy[i] =
            -2.0 * ts_xxxx_xyyy[i] * fbe_0 * fz_0 + tk_xxxx_xyyy[i] * fe_0 + tk_xxxxz_xyyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xyyy[i] * fz_0;

        tk_xxxxzz_xyyz[i] = -6.0 * ts_xxzz_xyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyz[i] * fe_0 + tk_xxxzz_yyz[i] * fe_0 +
                            tk_xxxzz_xyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyyz[i] * fz_0;

        tk_xxxxzz_xyzz[i] = -6.0 * ts_xxzz_xyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyzz[i] * fe_0 + tk_xxxzz_yzz[i] * fe_0 +
                            tk_xxxzz_xyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyzz[i] * fz_0;

        tk_xxxxzz_xzzz[i] = -6.0 * ts_xxzz_xzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xzzz[i] * fe_0 + tk_xxxzz_zzz[i] * fe_0 +
                            tk_xxxzz_xzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xzzz[i] * fz_0;

        tk_xxxxzz_yyyy[i] =
            -6.0 * ts_xxzz_yyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyy[i] * fe_0 + tk_xxxzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyyy[i] * fz_0;

        tk_xxxxzz_yyyz[i] =
            -6.0 * ts_xxzz_yyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyz[i] * fe_0 + tk_xxxzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyyz[i] * fz_0;

        tk_xxxxzz_yyzz[i] =
            -6.0 * ts_xxzz_yyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyzz[i] * fe_0 + tk_xxxzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyzz[i] * fz_0;

        tk_xxxxzz_yzzz[i] =
            -6.0 * ts_xxzz_yzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yzzz[i] * fe_0 + tk_xxxzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yzzz[i] * fz_0;

        tk_xxxxzz_zzzz[i] =
            -6.0 * ts_xxzz_zzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_zzzz[i] * fe_0 + tk_xxxzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_zzzz[i] * fz_0;
    }

    // Set up 90-105 components of targeted buffer : IG

    auto tk_xxxyyy_xxxx = pbuffer.data(idx_kin_ig + 90);

    auto tk_xxxyyy_xxxy = pbuffer.data(idx_kin_ig + 91);

    auto tk_xxxyyy_xxxz = pbuffer.data(idx_kin_ig + 92);

    auto tk_xxxyyy_xxyy = pbuffer.data(idx_kin_ig + 93);

    auto tk_xxxyyy_xxyz = pbuffer.data(idx_kin_ig + 94);

    auto tk_xxxyyy_xxzz = pbuffer.data(idx_kin_ig + 95);

    auto tk_xxxyyy_xyyy = pbuffer.data(idx_kin_ig + 96);

    auto tk_xxxyyy_xyyz = pbuffer.data(idx_kin_ig + 97);

    auto tk_xxxyyy_xyzz = pbuffer.data(idx_kin_ig + 98);

    auto tk_xxxyyy_xzzz = pbuffer.data(idx_kin_ig + 99);

    auto tk_xxxyyy_yyyy = pbuffer.data(idx_kin_ig + 100);

    auto tk_xxxyyy_yyyz = pbuffer.data(idx_kin_ig + 101);

    auto tk_xxxyyy_yyzz = pbuffer.data(idx_kin_ig + 102);

    auto tk_xxxyyy_yzzz = pbuffer.data(idx_kin_ig + 103);

    auto tk_xxxyyy_zzzz = pbuffer.data(idx_kin_ig + 104);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xxxy_xxxx,   \
                             tk_xxxy_xxxz,   \
                             tk_xxxy_xxzz,   \
                             tk_xxxy_xzzz,   \
                             tk_xxxyy_xxxx,  \
                             tk_xxxyy_xxxz,  \
                             tk_xxxyy_xxzz,  \
                             tk_xxxyy_xzzz,  \
                             tk_xxxyyy_xxxx, \
                             tk_xxxyyy_xxxy, \
                             tk_xxxyyy_xxxz, \
                             tk_xxxyyy_xxyy, \
                             tk_xxxyyy_xxyz, \
                             tk_xxxyyy_xxzz, \
                             tk_xxxyyy_xyyy, \
                             tk_xxxyyy_xyyz, \
                             tk_xxxyyy_xyzz, \
                             tk_xxxyyy_xzzz, \
                             tk_xxxyyy_yyyy, \
                             tk_xxxyyy_yyyz, \
                             tk_xxxyyy_yyzz, \
                             tk_xxxyyy_yzzz, \
                             tk_xxxyyy_zzzz, \
                             tk_xxyyy_xxxy,  \
                             tk_xxyyy_xxy,   \
                             tk_xxyyy_xxyy,  \
                             tk_xxyyy_xxyz,  \
                             tk_xxyyy_xyy,   \
                             tk_xxyyy_xyyy,  \
                             tk_xxyyy_xyyz,  \
                             tk_xxyyy_xyz,   \
                             tk_xxyyy_xyzz,  \
                             tk_xxyyy_yyy,   \
                             tk_xxyyy_yyyy,  \
                             tk_xxyyy_yyyz,  \
                             tk_xxyyy_yyz,   \
                             tk_xxyyy_yyzz,  \
                             tk_xxyyy_yzz,   \
                             tk_xxyyy_yzzz,  \
                             tk_xxyyy_zzzz,  \
                             tk_xyyy_xxxy,   \
                             tk_xyyy_xxyy,   \
                             tk_xyyy_xxyz,   \
                             tk_xyyy_xyyy,   \
                             tk_xyyy_xyyz,   \
                             tk_xyyy_xyzz,   \
                             tk_xyyy_yyyy,   \
                             tk_xyyy_yyyz,   \
                             tk_xyyy_yyzz,   \
                             tk_xyyy_yzzz,   \
                             tk_xyyy_zzzz,   \
                             ts_xxxy_xxxx,   \
                             ts_xxxy_xxxz,   \
                             ts_xxxy_xxzz,   \
                             ts_xxxy_xzzz,   \
                             ts_xxxyyy_xxxx, \
                             ts_xxxyyy_xxxy, \
                             ts_xxxyyy_xxxz, \
                             ts_xxxyyy_xxyy, \
                             ts_xxxyyy_xxyz, \
                             ts_xxxyyy_xxzz, \
                             ts_xxxyyy_xyyy, \
                             ts_xxxyyy_xyyz, \
                             ts_xxxyyy_xyzz, \
                             ts_xxxyyy_xzzz, \
                             ts_xxxyyy_yyyy, \
                             ts_xxxyyy_yyyz, \
                             ts_xxxyyy_yyzz, \
                             ts_xxxyyy_yzzz, \
                             ts_xxxyyy_zzzz, \
                             ts_xyyy_xxxy,   \
                             ts_xyyy_xxyy,   \
                             ts_xyyy_xxyz,   \
                             ts_xyyy_xyyy,   \
                             ts_xyyy_xyyz,   \
                             ts_xyyy_xyzz,   \
                             ts_xyyy_yyyy,   \
                             ts_xyyy_yyyz,   \
                             ts_xyyy_yyzz,   \
                             ts_xyyy_yzzz,   \
                             ts_xyyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyyy_xxxx[i] =
            -4.0 * ts_xxxy_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxx[i] * fe_0 + tk_xxxyy_xxxx[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxxx[i] * fz_0;

        tk_xxxyyy_xxxy[i] = -4.0 * ts_xyyy_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxy[i] * fe_0 + 3.0 * tk_xxyyy_xxy[i] * fe_0 +
                            tk_xxyyy_xxxy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxy[i] * fz_0;

        tk_xxxyyy_xxxz[i] =
            -4.0 * ts_xxxy_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxz[i] * fe_0 + tk_xxxyy_xxxz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxxz[i] * fz_0;

        tk_xxxyyy_xxyy[i] = -4.0 * ts_xyyy_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyy[i] * fe_0 + 2.0 * tk_xxyyy_xyy[i] * fe_0 +
                            tk_xxyyy_xxyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyy[i] * fz_0;

        tk_xxxyyy_xxyz[i] = -4.0 * ts_xyyy_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyz[i] * fe_0 + 2.0 * tk_xxyyy_xyz[i] * fe_0 +
                            tk_xxyyy_xxyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyz[i] * fz_0;

        tk_xxxyyy_xxzz[i] =
            -4.0 * ts_xxxy_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxzz[i] * fe_0 + tk_xxxyy_xxzz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxzz[i] * fz_0;

        tk_xxxyyy_xyyy[i] = -4.0 * ts_xyyy_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyy[i] * fe_0 + tk_xxyyy_yyy[i] * fe_0 +
                            tk_xxyyy_xyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyy[i] * fz_0;

        tk_xxxyyy_xyyz[i] = -4.0 * ts_xyyy_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyz[i] * fe_0 + tk_xxyyy_yyz[i] * fe_0 +
                            tk_xxyyy_xyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyz[i] * fz_0;

        tk_xxxyyy_xyzz[i] = -4.0 * ts_xyyy_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyzz[i] * fe_0 + tk_xxyyy_yzz[i] * fe_0 +
                            tk_xxyyy_xyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyzz[i] * fz_0;

        tk_xxxyyy_xzzz[i] =
            -4.0 * ts_xxxy_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xzzz[i] * fe_0 + tk_xxxyy_xzzz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xzzz[i] * fz_0;

        tk_xxxyyy_yyyy[i] =
            -4.0 * ts_xyyy_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyy[i] * fe_0 + tk_xxyyy_yyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyyy[i] * fz_0;

        tk_xxxyyy_yyyz[i] =
            -4.0 * ts_xyyy_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyz[i] * fe_0 + tk_xxyyy_yyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyyz[i] * fz_0;

        tk_xxxyyy_yyzz[i] =
            -4.0 * ts_xyyy_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyzz[i] * fe_0 + tk_xxyyy_yyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyzz[i] * fz_0;

        tk_xxxyyy_yzzz[i] =
            -4.0 * ts_xyyy_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yzzz[i] * fe_0 + tk_xxyyy_yzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yzzz[i] * fz_0;

        tk_xxxyyy_zzzz[i] =
            -4.0 * ts_xyyy_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_zzzz[i] * fe_0 + tk_xxyyy_zzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_zzzz[i] * fz_0;
    }

    // Set up 105-120 components of targeted buffer : IG

    auto tk_xxxyyz_xxxx = pbuffer.data(idx_kin_ig + 105);

    auto tk_xxxyyz_xxxy = pbuffer.data(idx_kin_ig + 106);

    auto tk_xxxyyz_xxxz = pbuffer.data(idx_kin_ig + 107);

    auto tk_xxxyyz_xxyy = pbuffer.data(idx_kin_ig + 108);

    auto tk_xxxyyz_xxyz = pbuffer.data(idx_kin_ig + 109);

    auto tk_xxxyyz_xxzz = pbuffer.data(idx_kin_ig + 110);

    auto tk_xxxyyz_xyyy = pbuffer.data(idx_kin_ig + 111);

    auto tk_xxxyyz_xyyz = pbuffer.data(idx_kin_ig + 112);

    auto tk_xxxyyz_xyzz = pbuffer.data(idx_kin_ig + 113);

    auto tk_xxxyyz_xzzz = pbuffer.data(idx_kin_ig + 114);

    auto tk_xxxyyz_yyyy = pbuffer.data(idx_kin_ig + 115);

    auto tk_xxxyyz_yyyz = pbuffer.data(idx_kin_ig + 116);

    auto tk_xxxyyz_yyzz = pbuffer.data(idx_kin_ig + 117);

    auto tk_xxxyyz_yzzz = pbuffer.data(idx_kin_ig + 118);

    auto tk_xxxyyz_zzzz = pbuffer.data(idx_kin_ig + 119);

#pragma omp simd aligned(pa_z,               \
                             tk_xxxyy_xxx,   \
                             tk_xxxyy_xxxx,  \
                             tk_xxxyy_xxxy,  \
                             tk_xxxyy_xxxz,  \
                             tk_xxxyy_xxy,   \
                             tk_xxxyy_xxyy,  \
                             tk_xxxyy_xxyz,  \
                             tk_xxxyy_xxz,   \
                             tk_xxxyy_xxzz,  \
                             tk_xxxyy_xyy,   \
                             tk_xxxyy_xyyy,  \
                             tk_xxxyy_xyyz,  \
                             tk_xxxyy_xyz,   \
                             tk_xxxyy_xyzz,  \
                             tk_xxxyy_xzz,   \
                             tk_xxxyy_xzzz,  \
                             tk_xxxyy_yyy,   \
                             tk_xxxyy_yyyy,  \
                             tk_xxxyy_yyyz,  \
                             tk_xxxyy_yyz,   \
                             tk_xxxyy_yyzz,  \
                             tk_xxxyy_yzz,   \
                             tk_xxxyy_yzzz,  \
                             tk_xxxyy_zzz,   \
                             tk_xxxyy_zzzz,  \
                             tk_xxxyyz_xxxx, \
                             tk_xxxyyz_xxxy, \
                             tk_xxxyyz_xxxz, \
                             tk_xxxyyz_xxyy, \
                             tk_xxxyyz_xxyz, \
                             tk_xxxyyz_xxzz, \
                             tk_xxxyyz_xyyy, \
                             tk_xxxyyz_xyyz, \
                             tk_xxxyyz_xyzz, \
                             tk_xxxyyz_xzzz, \
                             tk_xxxyyz_yyyy, \
                             tk_xxxyyz_yyyz, \
                             tk_xxxyyz_yyzz, \
                             tk_xxxyyz_yzzz, \
                             tk_xxxyyz_zzzz, \
                             ts_xxxyyz_xxxx, \
                             ts_xxxyyz_xxxy, \
                             ts_xxxyyz_xxxz, \
                             ts_xxxyyz_xxyy, \
                             ts_xxxyyz_xxyz, \
                             ts_xxxyyz_xxzz, \
                             ts_xxxyyz_xyyy, \
                             ts_xxxyyz_xyyz, \
                             ts_xxxyyz_xyzz, \
                             ts_xxxyyz_xzzz, \
                             ts_xxxyyz_yyyy, \
                             ts_xxxyyz_yyyz, \
                             ts_xxxyyz_yyzz, \
                             ts_xxxyyz_yzzz, \
                             ts_xxxyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyyz_xxxx[i] = tk_xxxyy_xxxx[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxx[i] * fz_0;

        tk_xxxyyz_xxxy[i] = tk_xxxyy_xxxy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxy[i] * fz_0;

        tk_xxxyyz_xxxz[i] = tk_xxxyy_xxx[i] * fe_0 + tk_xxxyy_xxxz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxz[i] * fz_0;

        tk_xxxyyz_xxyy[i] = tk_xxxyy_xxyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyy[i] * fz_0;

        tk_xxxyyz_xxyz[i] = tk_xxxyy_xxy[i] * fe_0 + tk_xxxyy_xxyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyz[i] * fz_0;

        tk_xxxyyz_xxzz[i] = 2.0 * tk_xxxyy_xxz[i] * fe_0 + tk_xxxyy_xxzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxzz[i] * fz_0;

        tk_xxxyyz_xyyy[i] = tk_xxxyy_xyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyy[i] * fz_0;

        tk_xxxyyz_xyyz[i] = tk_xxxyy_xyy[i] * fe_0 + tk_xxxyy_xyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyz[i] * fz_0;

        tk_xxxyyz_xyzz[i] = 2.0 * tk_xxxyy_xyz[i] * fe_0 + tk_xxxyy_xyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyzz[i] * fz_0;

        tk_xxxyyz_xzzz[i] = 3.0 * tk_xxxyy_xzz[i] * fe_0 + tk_xxxyy_xzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xzzz[i] * fz_0;

        tk_xxxyyz_yyyy[i] = tk_xxxyy_yyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyy[i] * fz_0;

        tk_xxxyyz_yyyz[i] = tk_xxxyy_yyy[i] * fe_0 + tk_xxxyy_yyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyz[i] * fz_0;

        tk_xxxyyz_yyzz[i] = 2.0 * tk_xxxyy_yyz[i] * fe_0 + tk_xxxyy_yyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyzz[i] * fz_0;

        tk_xxxyyz_yzzz[i] = 3.0 * tk_xxxyy_yzz[i] * fe_0 + tk_xxxyy_yzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yzzz[i] * fz_0;

        tk_xxxyyz_zzzz[i] = 4.0 * tk_xxxyy_zzz[i] * fe_0 + tk_xxxyy_zzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_zzzz[i] * fz_0;
    }

    // Set up 120-135 components of targeted buffer : IG

    auto tk_xxxyzz_xxxx = pbuffer.data(idx_kin_ig + 120);

    auto tk_xxxyzz_xxxy = pbuffer.data(idx_kin_ig + 121);

    auto tk_xxxyzz_xxxz = pbuffer.data(idx_kin_ig + 122);

    auto tk_xxxyzz_xxyy = pbuffer.data(idx_kin_ig + 123);

    auto tk_xxxyzz_xxyz = pbuffer.data(idx_kin_ig + 124);

    auto tk_xxxyzz_xxzz = pbuffer.data(idx_kin_ig + 125);

    auto tk_xxxyzz_xyyy = pbuffer.data(idx_kin_ig + 126);

    auto tk_xxxyzz_xyyz = pbuffer.data(idx_kin_ig + 127);

    auto tk_xxxyzz_xyzz = pbuffer.data(idx_kin_ig + 128);

    auto tk_xxxyzz_xzzz = pbuffer.data(idx_kin_ig + 129);

    auto tk_xxxyzz_yyyy = pbuffer.data(idx_kin_ig + 130);

    auto tk_xxxyzz_yyyz = pbuffer.data(idx_kin_ig + 131);

    auto tk_xxxyzz_yyzz = pbuffer.data(idx_kin_ig + 132);

    auto tk_xxxyzz_yzzz = pbuffer.data(idx_kin_ig + 133);

    auto tk_xxxyzz_zzzz = pbuffer.data(idx_kin_ig + 134);

#pragma omp simd aligned(pa_y,               \
                             tk_xxxyzz_xxxx, \
                             tk_xxxyzz_xxxy, \
                             tk_xxxyzz_xxxz, \
                             tk_xxxyzz_xxyy, \
                             tk_xxxyzz_xxyz, \
                             tk_xxxyzz_xxzz, \
                             tk_xxxyzz_xyyy, \
                             tk_xxxyzz_xyyz, \
                             tk_xxxyzz_xyzz, \
                             tk_xxxyzz_xzzz, \
                             tk_xxxyzz_yyyy, \
                             tk_xxxyzz_yyyz, \
                             tk_xxxyzz_yyzz, \
                             tk_xxxyzz_yzzz, \
                             tk_xxxyzz_zzzz, \
                             tk_xxxzz_xxx,   \
                             tk_xxxzz_xxxx,  \
                             tk_xxxzz_xxxy,  \
                             tk_xxxzz_xxxz,  \
                             tk_xxxzz_xxy,   \
                             tk_xxxzz_xxyy,  \
                             tk_xxxzz_xxyz,  \
                             tk_xxxzz_xxz,   \
                             tk_xxxzz_xxzz,  \
                             tk_xxxzz_xyy,   \
                             tk_xxxzz_xyyy,  \
                             tk_xxxzz_xyyz,  \
                             tk_xxxzz_xyz,   \
                             tk_xxxzz_xyzz,  \
                             tk_xxxzz_xzz,   \
                             tk_xxxzz_xzzz,  \
                             tk_xxxzz_yyy,   \
                             tk_xxxzz_yyyy,  \
                             tk_xxxzz_yyyz,  \
                             tk_xxxzz_yyz,   \
                             tk_xxxzz_yyzz,  \
                             tk_xxxzz_yzz,   \
                             tk_xxxzz_yzzz,  \
                             tk_xxxzz_zzz,   \
                             tk_xxxzz_zzzz,  \
                             ts_xxxyzz_xxxx, \
                             ts_xxxyzz_xxxy, \
                             ts_xxxyzz_xxxz, \
                             ts_xxxyzz_xxyy, \
                             ts_xxxyzz_xxyz, \
                             ts_xxxyzz_xxzz, \
                             ts_xxxyzz_xyyy, \
                             ts_xxxyzz_xyyz, \
                             ts_xxxyzz_xyzz, \
                             ts_xxxyzz_xzzz, \
                             ts_xxxyzz_yyyy, \
                             ts_xxxyzz_yyyz, \
                             ts_xxxyzz_yyzz, \
                             ts_xxxyzz_yzzz, \
                             ts_xxxyzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyzz_xxxx[i] = tk_xxxzz_xxxx[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxx[i] * fz_0;

        tk_xxxyzz_xxxy[i] = tk_xxxzz_xxx[i] * fe_0 + tk_xxxzz_xxxy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxy[i] * fz_0;

        tk_xxxyzz_xxxz[i] = tk_xxxzz_xxxz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxz[i] * fz_0;

        tk_xxxyzz_xxyy[i] = 2.0 * tk_xxxzz_xxy[i] * fe_0 + tk_xxxzz_xxyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyy[i] * fz_0;

        tk_xxxyzz_xxyz[i] = tk_xxxzz_xxz[i] * fe_0 + tk_xxxzz_xxyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyz[i] * fz_0;

        tk_xxxyzz_xxzz[i] = tk_xxxzz_xxzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxzz[i] * fz_0;

        tk_xxxyzz_xyyy[i] = 3.0 * tk_xxxzz_xyy[i] * fe_0 + tk_xxxzz_xyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyy[i] * fz_0;

        tk_xxxyzz_xyyz[i] = 2.0 * tk_xxxzz_xyz[i] * fe_0 + tk_xxxzz_xyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyz[i] * fz_0;

        tk_xxxyzz_xyzz[i] = tk_xxxzz_xzz[i] * fe_0 + tk_xxxzz_xyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyzz[i] * fz_0;

        tk_xxxyzz_xzzz[i] = tk_xxxzz_xzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xzzz[i] * fz_0;

        tk_xxxyzz_yyyy[i] = 4.0 * tk_xxxzz_yyy[i] * fe_0 + tk_xxxzz_yyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyy[i] * fz_0;

        tk_xxxyzz_yyyz[i] = 3.0 * tk_xxxzz_yyz[i] * fe_0 + tk_xxxzz_yyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyz[i] * fz_0;

        tk_xxxyzz_yyzz[i] = 2.0 * tk_xxxzz_yzz[i] * fe_0 + tk_xxxzz_yyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyzz[i] * fz_0;

        tk_xxxyzz_yzzz[i] = tk_xxxzz_zzz[i] * fe_0 + tk_xxxzz_yzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yzzz[i] * fz_0;

        tk_xxxyzz_zzzz[i] = tk_xxxzz_zzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_zzzz[i] * fz_0;
    }

    // Set up 135-150 components of targeted buffer : IG

    auto tk_xxxzzz_xxxx = pbuffer.data(idx_kin_ig + 135);

    auto tk_xxxzzz_xxxy = pbuffer.data(idx_kin_ig + 136);

    auto tk_xxxzzz_xxxz = pbuffer.data(idx_kin_ig + 137);

    auto tk_xxxzzz_xxyy = pbuffer.data(idx_kin_ig + 138);

    auto tk_xxxzzz_xxyz = pbuffer.data(idx_kin_ig + 139);

    auto tk_xxxzzz_xxzz = pbuffer.data(idx_kin_ig + 140);

    auto tk_xxxzzz_xyyy = pbuffer.data(idx_kin_ig + 141);

    auto tk_xxxzzz_xyyz = pbuffer.data(idx_kin_ig + 142);

    auto tk_xxxzzz_xyzz = pbuffer.data(idx_kin_ig + 143);

    auto tk_xxxzzz_xzzz = pbuffer.data(idx_kin_ig + 144);

    auto tk_xxxzzz_yyyy = pbuffer.data(idx_kin_ig + 145);

    auto tk_xxxzzz_yyyz = pbuffer.data(idx_kin_ig + 146);

    auto tk_xxxzzz_yyzz = pbuffer.data(idx_kin_ig + 147);

    auto tk_xxxzzz_yzzz = pbuffer.data(idx_kin_ig + 148);

    auto tk_xxxzzz_zzzz = pbuffer.data(idx_kin_ig + 149);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xxxz_xxxx,   \
                             tk_xxxz_xxxy,   \
                             tk_xxxz_xxyy,   \
                             tk_xxxz_xyyy,   \
                             tk_xxxzz_xxxx,  \
                             tk_xxxzz_xxxy,  \
                             tk_xxxzz_xxyy,  \
                             tk_xxxzz_xyyy,  \
                             tk_xxxzzz_xxxx, \
                             tk_xxxzzz_xxxy, \
                             tk_xxxzzz_xxxz, \
                             tk_xxxzzz_xxyy, \
                             tk_xxxzzz_xxyz, \
                             tk_xxxzzz_xxzz, \
                             tk_xxxzzz_xyyy, \
                             tk_xxxzzz_xyyz, \
                             tk_xxxzzz_xyzz, \
                             tk_xxxzzz_xzzz, \
                             tk_xxxzzz_yyyy, \
                             tk_xxxzzz_yyyz, \
                             tk_xxxzzz_yyzz, \
                             tk_xxxzzz_yzzz, \
                             tk_xxxzzz_zzzz, \
                             tk_xxzzz_xxxz,  \
                             tk_xxzzz_xxyz,  \
                             tk_xxzzz_xxz,   \
                             tk_xxzzz_xxzz,  \
                             tk_xxzzz_xyyz,  \
                             tk_xxzzz_xyz,   \
                             tk_xxzzz_xyzz,  \
                             tk_xxzzz_xzz,   \
                             tk_xxzzz_xzzz,  \
                             tk_xxzzz_yyyy,  \
                             tk_xxzzz_yyyz,  \
                             tk_xxzzz_yyz,   \
                             tk_xxzzz_yyzz,  \
                             tk_xxzzz_yzz,   \
                             tk_xxzzz_yzzz,  \
                             tk_xxzzz_zzz,   \
                             tk_xxzzz_zzzz,  \
                             tk_xzzz_xxxz,   \
                             tk_xzzz_xxyz,   \
                             tk_xzzz_xxzz,   \
                             tk_xzzz_xyyz,   \
                             tk_xzzz_xyzz,   \
                             tk_xzzz_xzzz,   \
                             tk_xzzz_yyyy,   \
                             tk_xzzz_yyyz,   \
                             tk_xzzz_yyzz,   \
                             tk_xzzz_yzzz,   \
                             tk_xzzz_zzzz,   \
                             ts_xxxz_xxxx,   \
                             ts_xxxz_xxxy,   \
                             ts_xxxz_xxyy,   \
                             ts_xxxz_xyyy,   \
                             ts_xxxzzz_xxxx, \
                             ts_xxxzzz_xxxy, \
                             ts_xxxzzz_xxxz, \
                             ts_xxxzzz_xxyy, \
                             ts_xxxzzz_xxyz, \
                             ts_xxxzzz_xxzz, \
                             ts_xxxzzz_xyyy, \
                             ts_xxxzzz_xyyz, \
                             ts_xxxzzz_xyzz, \
                             ts_xxxzzz_xzzz, \
                             ts_xxxzzz_yyyy, \
                             ts_xxxzzz_yyyz, \
                             ts_xxxzzz_yyzz, \
                             ts_xxxzzz_yzzz, \
                             ts_xxxzzz_zzzz, \
                             ts_xzzz_xxxz,   \
                             ts_xzzz_xxyz,   \
                             ts_xzzz_xxzz,   \
                             ts_xzzz_xyyz,   \
                             ts_xzzz_xyzz,   \
                             ts_xzzz_xzzz,   \
                             ts_xzzz_yyyy,   \
                             ts_xzzz_yyyz,   \
                             ts_xzzz_yyzz,   \
                             ts_xzzz_yzzz,   \
                             ts_xzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzzz_xxxx[i] =
            -4.0 * ts_xxxz_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxx[i] * fe_0 + tk_xxxzz_xxxx[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxxx[i] * fz_0;

        tk_xxxzzz_xxxy[i] =
            -4.0 * ts_xxxz_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxy[i] * fe_0 + tk_xxxzz_xxxy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxxy[i] * fz_0;

        tk_xxxzzz_xxxz[i] = -4.0 * ts_xzzz_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxz[i] * fe_0 + 3.0 * tk_xxzzz_xxz[i] * fe_0 +
                            tk_xxzzz_xxxz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxz[i] * fz_0;

        tk_xxxzzz_xxyy[i] =
            -4.0 * ts_xxxz_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxyy[i] * fe_0 + tk_xxxzz_xxyy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxyy[i] * fz_0;

        tk_xxxzzz_xxyz[i] = -4.0 * ts_xzzz_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxyz[i] * fe_0 + 2.0 * tk_xxzzz_xyz[i] * fe_0 +
                            tk_xxzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxyz[i] * fz_0;

        tk_xxxzzz_xxzz[i] = -4.0 * ts_xzzz_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxzz[i] * fe_0 + 2.0 * tk_xxzzz_xzz[i] * fe_0 +
                            tk_xxzzz_xxzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxzz[i] * fz_0;

        tk_xxxzzz_xyyy[i] =
            -4.0 * ts_xxxz_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xyyy[i] * fe_0 + tk_xxxzz_xyyy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xyyy[i] * fz_0;

        tk_xxxzzz_xyyz[i] = -4.0 * ts_xzzz_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyyz[i] * fe_0 + tk_xxzzz_yyz[i] * fe_0 +
                            tk_xxzzz_xyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyyz[i] * fz_0;

        tk_xxxzzz_xyzz[i] = -4.0 * ts_xzzz_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyzz[i] * fe_0 + tk_xxzzz_yzz[i] * fe_0 +
                            tk_xxzzz_xyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyzz[i] * fz_0;

        tk_xxxzzz_xzzz[i] = -4.0 * ts_xzzz_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xzzz[i] * fe_0 + tk_xxzzz_zzz[i] * fe_0 +
                            tk_xxzzz_xzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xzzz[i] * fz_0;

        tk_xxxzzz_yyyy[i] =
            -4.0 * ts_xzzz_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyy[i] * fe_0 + tk_xxzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyyy[i] * fz_0;

        tk_xxxzzz_yyyz[i] =
            -4.0 * ts_xzzz_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyz[i] * fe_0 + tk_xxzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyyz[i] * fz_0;

        tk_xxxzzz_yyzz[i] =
            -4.0 * ts_xzzz_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyzz[i] * fe_0 + tk_xxzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyzz[i] * fz_0;

        tk_xxxzzz_yzzz[i] =
            -4.0 * ts_xzzz_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yzzz[i] * fe_0 + tk_xxzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yzzz[i] * fz_0;

        tk_xxxzzz_zzzz[i] =
            -4.0 * ts_xzzz_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_zzzz[i] * fe_0 + tk_xxzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_zzzz[i] * fz_0;
    }

    // Set up 150-165 components of targeted buffer : IG

    auto tk_xxyyyy_xxxx = pbuffer.data(idx_kin_ig + 150);

    auto tk_xxyyyy_xxxy = pbuffer.data(idx_kin_ig + 151);

    auto tk_xxyyyy_xxxz = pbuffer.data(idx_kin_ig + 152);

    auto tk_xxyyyy_xxyy = pbuffer.data(idx_kin_ig + 153);

    auto tk_xxyyyy_xxyz = pbuffer.data(idx_kin_ig + 154);

    auto tk_xxyyyy_xxzz = pbuffer.data(idx_kin_ig + 155);

    auto tk_xxyyyy_xyyy = pbuffer.data(idx_kin_ig + 156);

    auto tk_xxyyyy_xyyz = pbuffer.data(idx_kin_ig + 157);

    auto tk_xxyyyy_xyzz = pbuffer.data(idx_kin_ig + 158);

    auto tk_xxyyyy_xzzz = pbuffer.data(idx_kin_ig + 159);

    auto tk_xxyyyy_yyyy = pbuffer.data(idx_kin_ig + 160);

    auto tk_xxyyyy_yyyz = pbuffer.data(idx_kin_ig + 161);

    auto tk_xxyyyy_yyzz = pbuffer.data(idx_kin_ig + 162);

    auto tk_xxyyyy_yzzz = pbuffer.data(idx_kin_ig + 163);

    auto tk_xxyyyy_zzzz = pbuffer.data(idx_kin_ig + 164);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xxyy_xxxx,   \
                             tk_xxyy_xxxz,   \
                             tk_xxyy_xxzz,   \
                             tk_xxyy_xzzz,   \
                             tk_xxyyy_xxxx,  \
                             tk_xxyyy_xxxz,  \
                             tk_xxyyy_xxzz,  \
                             tk_xxyyy_xzzz,  \
                             tk_xxyyyy_xxxx, \
                             tk_xxyyyy_xxxy, \
                             tk_xxyyyy_xxxz, \
                             tk_xxyyyy_xxyy, \
                             tk_xxyyyy_xxyz, \
                             tk_xxyyyy_xxzz, \
                             tk_xxyyyy_xyyy, \
                             tk_xxyyyy_xyyz, \
                             tk_xxyyyy_xyzz, \
                             tk_xxyyyy_xzzz, \
                             tk_xxyyyy_yyyy, \
                             tk_xxyyyy_yyyz, \
                             tk_xxyyyy_yyzz, \
                             tk_xxyyyy_yzzz, \
                             tk_xxyyyy_zzzz, \
                             tk_xyyyy_xxxy,  \
                             tk_xyyyy_xxy,   \
                             tk_xyyyy_xxyy,  \
                             tk_xyyyy_xxyz,  \
                             tk_xyyyy_xyy,   \
                             tk_xyyyy_xyyy,  \
                             tk_xyyyy_xyyz,  \
                             tk_xyyyy_xyz,   \
                             tk_xyyyy_xyzz,  \
                             tk_xyyyy_yyy,   \
                             tk_xyyyy_yyyy,  \
                             tk_xyyyy_yyyz,  \
                             tk_xyyyy_yyz,   \
                             tk_xyyyy_yyzz,  \
                             tk_xyyyy_yzz,   \
                             tk_xyyyy_yzzz,  \
                             tk_xyyyy_zzzz,  \
                             tk_yyyy_xxxy,   \
                             tk_yyyy_xxyy,   \
                             tk_yyyy_xxyz,   \
                             tk_yyyy_xyyy,   \
                             tk_yyyy_xyyz,   \
                             tk_yyyy_xyzz,   \
                             tk_yyyy_yyyy,   \
                             tk_yyyy_yyyz,   \
                             tk_yyyy_yyzz,   \
                             tk_yyyy_yzzz,   \
                             tk_yyyy_zzzz,   \
                             ts_xxyy_xxxx,   \
                             ts_xxyy_xxxz,   \
                             ts_xxyy_xxzz,   \
                             ts_xxyy_xzzz,   \
                             ts_xxyyyy_xxxx, \
                             ts_xxyyyy_xxxy, \
                             ts_xxyyyy_xxxz, \
                             ts_xxyyyy_xxyy, \
                             ts_xxyyyy_xxyz, \
                             ts_xxyyyy_xxzz, \
                             ts_xxyyyy_xyyy, \
                             ts_xxyyyy_xyyz, \
                             ts_xxyyyy_xyzz, \
                             ts_xxyyyy_xzzz, \
                             ts_xxyyyy_yyyy, \
                             ts_xxyyyy_yyyz, \
                             ts_xxyyyy_yyzz, \
                             ts_xxyyyy_yzzz, \
                             ts_xxyyyy_zzzz, \
                             ts_yyyy_xxxy,   \
                             ts_yyyy_xxyy,   \
                             ts_yyyy_xxyz,   \
                             ts_yyyy_xyyy,   \
                             ts_yyyy_xyyz,   \
                             ts_yyyy_xyzz,   \
                             ts_yyyy_yyyy,   \
                             ts_yyyy_yyyz,   \
                             ts_yyyy_yyzz,   \
                             ts_yyyy_yzzz,   \
                             ts_yyyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyyy_xxxx[i] =
            -6.0 * ts_xxyy_xxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxx[i] * fe_0 + tk_xxyyy_xxxx[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxxx[i] * fz_0;

        tk_xxyyyy_xxxy[i] = -2.0 * ts_yyyy_xxxy[i] * fbe_0 * fz_0 + tk_yyyy_xxxy[i] * fe_0 + 3.0 * tk_xyyyy_xxy[i] * fe_0 +
                            tk_xyyyy_xxxy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxy[i] * fz_0;

        tk_xxyyyy_xxxz[i] =
            -6.0 * ts_xxyy_xxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxz[i] * fe_0 + tk_xxyyy_xxxz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxxz[i] * fz_0;

        tk_xxyyyy_xxyy[i] = -2.0 * ts_yyyy_xxyy[i] * fbe_0 * fz_0 + tk_yyyy_xxyy[i] * fe_0 + 2.0 * tk_xyyyy_xyy[i] * fe_0 +
                            tk_xyyyy_xxyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyy[i] * fz_0;

        tk_xxyyyy_xxyz[i] = -2.0 * ts_yyyy_xxyz[i] * fbe_0 * fz_0 + tk_yyyy_xxyz[i] * fe_0 + 2.0 * tk_xyyyy_xyz[i] * fe_0 +
                            tk_xyyyy_xxyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyz[i] * fz_0;

        tk_xxyyyy_xxzz[i] =
            -6.0 * ts_xxyy_xxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxzz[i] * fe_0 + tk_xxyyy_xxzz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxzz[i] * fz_0;

        tk_xxyyyy_xyyy[i] = -2.0 * ts_yyyy_xyyy[i] * fbe_0 * fz_0 + tk_yyyy_xyyy[i] * fe_0 + tk_xyyyy_yyy[i] * fe_0 + tk_xyyyy_xyyy[i] * pa_x[i] +
                            2.0 * ts_xxyyyy_xyyy[i] * fz_0;

        tk_xxyyyy_xyyz[i] = -2.0 * ts_yyyy_xyyz[i] * fbe_0 * fz_0 + tk_yyyy_xyyz[i] * fe_0 + tk_xyyyy_yyz[i] * fe_0 + tk_xyyyy_xyyz[i] * pa_x[i] +
                            2.0 * ts_xxyyyy_xyyz[i] * fz_0;

        tk_xxyyyy_xyzz[i] = -2.0 * ts_yyyy_xyzz[i] * fbe_0 * fz_0 + tk_yyyy_xyzz[i] * fe_0 + tk_xyyyy_yzz[i] * fe_0 + tk_xyyyy_xyzz[i] * pa_x[i] +
                            2.0 * ts_xxyyyy_xyzz[i] * fz_0;

        tk_xxyyyy_xzzz[i] =
            -6.0 * ts_xxyy_xzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xzzz[i] * fe_0 + tk_xxyyy_xzzz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xzzz[i] * fz_0;

        tk_xxyyyy_yyyy[i] =
            -2.0 * ts_yyyy_yyyy[i] * fbe_0 * fz_0 + tk_yyyy_yyyy[i] * fe_0 + tk_xyyyy_yyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyy[i] * fz_0;

        tk_xxyyyy_yyyz[i] =
            -2.0 * ts_yyyy_yyyz[i] * fbe_0 * fz_0 + tk_yyyy_yyyz[i] * fe_0 + tk_xyyyy_yyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyz[i] * fz_0;

        tk_xxyyyy_yyzz[i] =
            -2.0 * ts_yyyy_yyzz[i] * fbe_0 * fz_0 + tk_yyyy_yyzz[i] * fe_0 + tk_xyyyy_yyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyzz[i] * fz_0;

        tk_xxyyyy_yzzz[i] =
            -2.0 * ts_yyyy_yzzz[i] * fbe_0 * fz_0 + tk_yyyy_yzzz[i] * fe_0 + tk_xyyyy_yzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yzzz[i] * fz_0;

        tk_xxyyyy_zzzz[i] =
            -2.0 * ts_yyyy_zzzz[i] * fbe_0 * fz_0 + tk_yyyy_zzzz[i] * fe_0 + tk_xyyyy_zzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_zzzz[i] * fz_0;
    }

    // Set up 165-180 components of targeted buffer : IG

    auto tk_xxyyyz_xxxx = pbuffer.data(idx_kin_ig + 165);

    auto tk_xxyyyz_xxxy = pbuffer.data(idx_kin_ig + 166);

    auto tk_xxyyyz_xxxz = pbuffer.data(idx_kin_ig + 167);

    auto tk_xxyyyz_xxyy = pbuffer.data(idx_kin_ig + 168);

    auto tk_xxyyyz_xxyz = pbuffer.data(idx_kin_ig + 169);

    auto tk_xxyyyz_xxzz = pbuffer.data(idx_kin_ig + 170);

    auto tk_xxyyyz_xyyy = pbuffer.data(idx_kin_ig + 171);

    auto tk_xxyyyz_xyyz = pbuffer.data(idx_kin_ig + 172);

    auto tk_xxyyyz_xyzz = pbuffer.data(idx_kin_ig + 173);

    auto tk_xxyyyz_xzzz = pbuffer.data(idx_kin_ig + 174);

    auto tk_xxyyyz_yyyy = pbuffer.data(idx_kin_ig + 175);

    auto tk_xxyyyz_yyyz = pbuffer.data(idx_kin_ig + 176);

    auto tk_xxyyyz_yyzz = pbuffer.data(idx_kin_ig + 177);

    auto tk_xxyyyz_yzzz = pbuffer.data(idx_kin_ig + 178);

    auto tk_xxyyyz_zzzz = pbuffer.data(idx_kin_ig + 179);

#pragma omp simd aligned(pa_z,               \
                             tk_xxyyy_xxx,   \
                             tk_xxyyy_xxxx,  \
                             tk_xxyyy_xxxy,  \
                             tk_xxyyy_xxxz,  \
                             tk_xxyyy_xxy,   \
                             tk_xxyyy_xxyy,  \
                             tk_xxyyy_xxyz,  \
                             tk_xxyyy_xxz,   \
                             tk_xxyyy_xxzz,  \
                             tk_xxyyy_xyy,   \
                             tk_xxyyy_xyyy,  \
                             tk_xxyyy_xyyz,  \
                             tk_xxyyy_xyz,   \
                             tk_xxyyy_xyzz,  \
                             tk_xxyyy_xzz,   \
                             tk_xxyyy_xzzz,  \
                             tk_xxyyy_yyy,   \
                             tk_xxyyy_yyyy,  \
                             tk_xxyyy_yyyz,  \
                             tk_xxyyy_yyz,   \
                             tk_xxyyy_yyzz,  \
                             tk_xxyyy_yzz,   \
                             tk_xxyyy_yzzz,  \
                             tk_xxyyy_zzz,   \
                             tk_xxyyy_zzzz,  \
                             tk_xxyyyz_xxxx, \
                             tk_xxyyyz_xxxy, \
                             tk_xxyyyz_xxxz, \
                             tk_xxyyyz_xxyy, \
                             tk_xxyyyz_xxyz, \
                             tk_xxyyyz_xxzz, \
                             tk_xxyyyz_xyyy, \
                             tk_xxyyyz_xyyz, \
                             tk_xxyyyz_xyzz, \
                             tk_xxyyyz_xzzz, \
                             tk_xxyyyz_yyyy, \
                             tk_xxyyyz_yyyz, \
                             tk_xxyyyz_yyzz, \
                             tk_xxyyyz_yzzz, \
                             tk_xxyyyz_zzzz, \
                             ts_xxyyyz_xxxx, \
                             ts_xxyyyz_xxxy, \
                             ts_xxyyyz_xxxz, \
                             ts_xxyyyz_xxyy, \
                             ts_xxyyyz_xxyz, \
                             ts_xxyyyz_xxzz, \
                             ts_xxyyyz_xyyy, \
                             ts_xxyyyz_xyyz, \
                             ts_xxyyyz_xyzz, \
                             ts_xxyyyz_xzzz, \
                             ts_xxyyyz_yyyy, \
                             ts_xxyyyz_yyyz, \
                             ts_xxyyyz_yyzz, \
                             ts_xxyyyz_yzzz, \
                             ts_xxyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyyz_xxxx[i] = tk_xxyyy_xxxx[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxx[i] * fz_0;

        tk_xxyyyz_xxxy[i] = tk_xxyyy_xxxy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxy[i] * fz_0;

        tk_xxyyyz_xxxz[i] = tk_xxyyy_xxx[i] * fe_0 + tk_xxyyy_xxxz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxz[i] * fz_0;

        tk_xxyyyz_xxyy[i] = tk_xxyyy_xxyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyy[i] * fz_0;

        tk_xxyyyz_xxyz[i] = tk_xxyyy_xxy[i] * fe_0 + tk_xxyyy_xxyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyz[i] * fz_0;

        tk_xxyyyz_xxzz[i] = 2.0 * tk_xxyyy_xxz[i] * fe_0 + tk_xxyyy_xxzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxzz[i] * fz_0;

        tk_xxyyyz_xyyy[i] = tk_xxyyy_xyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyy[i] * fz_0;

        tk_xxyyyz_xyyz[i] = tk_xxyyy_xyy[i] * fe_0 + tk_xxyyy_xyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyz[i] * fz_0;

        tk_xxyyyz_xyzz[i] = 2.0 * tk_xxyyy_xyz[i] * fe_0 + tk_xxyyy_xyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyzz[i] * fz_0;

        tk_xxyyyz_xzzz[i] = 3.0 * tk_xxyyy_xzz[i] * fe_0 + tk_xxyyy_xzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xzzz[i] * fz_0;

        tk_xxyyyz_yyyy[i] = tk_xxyyy_yyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyy[i] * fz_0;

        tk_xxyyyz_yyyz[i] = tk_xxyyy_yyy[i] * fe_0 + tk_xxyyy_yyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyz[i] * fz_0;

        tk_xxyyyz_yyzz[i] = 2.0 * tk_xxyyy_yyz[i] * fe_0 + tk_xxyyy_yyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyzz[i] * fz_0;

        tk_xxyyyz_yzzz[i] = 3.0 * tk_xxyyy_yzz[i] * fe_0 + tk_xxyyy_yzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yzzz[i] * fz_0;

        tk_xxyyyz_zzzz[i] = 4.0 * tk_xxyyy_zzz[i] * fe_0 + tk_xxyyy_zzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_zzzz[i] * fz_0;
    }

    // Set up 180-195 components of targeted buffer : IG

    auto tk_xxyyzz_xxxx = pbuffer.data(idx_kin_ig + 180);

    auto tk_xxyyzz_xxxy = pbuffer.data(idx_kin_ig + 181);

    auto tk_xxyyzz_xxxz = pbuffer.data(idx_kin_ig + 182);

    auto tk_xxyyzz_xxyy = pbuffer.data(idx_kin_ig + 183);

    auto tk_xxyyzz_xxyz = pbuffer.data(idx_kin_ig + 184);

    auto tk_xxyyzz_xxzz = pbuffer.data(idx_kin_ig + 185);

    auto tk_xxyyzz_xyyy = pbuffer.data(idx_kin_ig + 186);

    auto tk_xxyyzz_xyyz = pbuffer.data(idx_kin_ig + 187);

    auto tk_xxyyzz_xyzz = pbuffer.data(idx_kin_ig + 188);

    auto tk_xxyyzz_xzzz = pbuffer.data(idx_kin_ig + 189);

    auto tk_xxyyzz_yyyy = pbuffer.data(idx_kin_ig + 190);

    auto tk_xxyyzz_yyyz = pbuffer.data(idx_kin_ig + 191);

    auto tk_xxyyzz_yyzz = pbuffer.data(idx_kin_ig + 192);

    auto tk_xxyyzz_yzzz = pbuffer.data(idx_kin_ig + 193);

    auto tk_xxyyzz_zzzz = pbuffer.data(idx_kin_ig + 194);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tk_xxyy_xxxy,   \
                             tk_xxyy_xxyy,   \
                             tk_xxyy_xyyy,   \
                             tk_xxyyz_xxxy,  \
                             tk_xxyyz_xxyy,  \
                             tk_xxyyz_xyyy,  \
                             tk_xxyyzz_xxxx, \
                             tk_xxyyzz_xxxy, \
                             tk_xxyyzz_xxxz, \
                             tk_xxyyzz_xxyy, \
                             tk_xxyyzz_xxyz, \
                             tk_xxyyzz_xxzz, \
                             tk_xxyyzz_xyyy, \
                             tk_xxyyzz_xyyz, \
                             tk_xxyyzz_xyzz, \
                             tk_xxyyzz_xzzz, \
                             tk_xxyyzz_yyyy, \
                             tk_xxyyzz_yyyz, \
                             tk_xxyyzz_yyzz, \
                             tk_xxyyzz_yzzz, \
                             tk_xxyyzz_zzzz, \
                             tk_xxyzz_xxxx,  \
                             tk_xxyzz_xxxz,  \
                             tk_xxyzz_xxzz,  \
                             tk_xxyzz_xzzz,  \
                             tk_xxzz_xxxx,   \
                             tk_xxzz_xxxz,   \
                             tk_xxzz_xxzz,   \
                             tk_xxzz_xzzz,   \
                             tk_xyyzz_xxyz,  \
                             tk_xyyzz_xyyz,  \
                             tk_xyyzz_xyz,   \
                             tk_xyyzz_xyzz,  \
                             tk_xyyzz_yyyy,  \
                             tk_xyyzz_yyyz,  \
                             tk_xyyzz_yyz,   \
                             tk_xyyzz_yyzz,  \
                             tk_xyyzz_yzz,   \
                             tk_xyyzz_yzzz,  \
                             tk_xyyzz_zzzz,  \
                             tk_yyzz_xxyz,   \
                             tk_yyzz_xyyz,   \
                             tk_yyzz_xyzz,   \
                             tk_yyzz_yyyy,   \
                             tk_yyzz_yyyz,   \
                             tk_yyzz_yyzz,   \
                             tk_yyzz_yzzz,   \
                             tk_yyzz_zzzz,   \
                             ts_xxyy_xxxy,   \
                             ts_xxyy_xxyy,   \
                             ts_xxyy_xyyy,   \
                             ts_xxyyzz_xxxx, \
                             ts_xxyyzz_xxxy, \
                             ts_xxyyzz_xxxz, \
                             ts_xxyyzz_xxyy, \
                             ts_xxyyzz_xxyz, \
                             ts_xxyyzz_xxzz, \
                             ts_xxyyzz_xyyy, \
                             ts_xxyyzz_xyyz, \
                             ts_xxyyzz_xyzz, \
                             ts_xxyyzz_xzzz, \
                             ts_xxyyzz_yyyy, \
                             ts_xxyyzz_yyyz, \
                             ts_xxyyzz_yyzz, \
                             ts_xxyyzz_yzzz, \
                             ts_xxyyzz_zzzz, \
                             ts_xxzz_xxxx,   \
                             ts_xxzz_xxxz,   \
                             ts_xxzz_xxzz,   \
                             ts_xxzz_xzzz,   \
                             ts_yyzz_xxyz,   \
                             ts_yyzz_xyyz,   \
                             ts_yyzz_xyzz,   \
                             ts_yyzz_yyyy,   \
                             ts_yyzz_yyyz,   \
                             ts_yyzz_yyzz,   \
                             ts_yyzz_yzzz,   \
                             ts_yyzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyzz_xxxx[i] =
            -2.0 * ts_xxzz_xxxx[i] * fbe_0 * fz_0 + tk_xxzz_xxxx[i] * fe_0 + tk_xxyzz_xxxx[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxx[i] * fz_0;

        tk_xxyyzz_xxxy[i] =
            -2.0 * ts_xxyy_xxxy[i] * fbe_0 * fz_0 + tk_xxyy_xxxy[i] * fe_0 + tk_xxyyz_xxxy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxxy[i] * fz_0;

        tk_xxyyzz_xxxz[i] =
            -2.0 * ts_xxzz_xxxz[i] * fbe_0 * fz_0 + tk_xxzz_xxxz[i] * fe_0 + tk_xxyzz_xxxz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxz[i] * fz_0;

        tk_xxyyzz_xxyy[i] =
            -2.0 * ts_xxyy_xxyy[i] * fbe_0 * fz_0 + tk_xxyy_xxyy[i] * fe_0 + tk_xxyyz_xxyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxyy[i] * fz_0;

        tk_xxyyzz_xxyz[i] = -2.0 * ts_yyzz_xxyz[i] * fbe_0 * fz_0 + tk_yyzz_xxyz[i] * fe_0 + 2.0 * tk_xyyzz_xyz[i] * fe_0 +
                            tk_xyyzz_xxyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxyz[i] * fz_0;

        tk_xxyyzz_xxzz[i] =
            -2.0 * ts_xxzz_xxzz[i] * fbe_0 * fz_0 + tk_xxzz_xxzz[i] * fe_0 + tk_xxyzz_xxzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxzz[i] * fz_0;

        tk_xxyyzz_xyyy[i] =
            -2.0 * ts_xxyy_xyyy[i] * fbe_0 * fz_0 + tk_xxyy_xyyy[i] * fe_0 + tk_xxyyz_xyyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xyyy[i] * fz_0;

        tk_xxyyzz_xyyz[i] = -2.0 * ts_yyzz_xyyz[i] * fbe_0 * fz_0 + tk_yyzz_xyyz[i] * fe_0 + tk_xyyzz_yyz[i] * fe_0 + tk_xyyzz_xyyz[i] * pa_x[i] +
                            2.0 * ts_xxyyzz_xyyz[i] * fz_0;

        tk_xxyyzz_xyzz[i] = -2.0 * ts_yyzz_xyzz[i] * fbe_0 * fz_0 + tk_yyzz_xyzz[i] * fe_0 + tk_xyyzz_yzz[i] * fe_0 + tk_xyyzz_xyzz[i] * pa_x[i] +
                            2.0 * ts_xxyyzz_xyzz[i] * fz_0;

        tk_xxyyzz_xzzz[i] =
            -2.0 * ts_xxzz_xzzz[i] * fbe_0 * fz_0 + tk_xxzz_xzzz[i] * fe_0 + tk_xxyzz_xzzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xzzz[i] * fz_0;

        tk_xxyyzz_yyyy[i] =
            -2.0 * ts_yyzz_yyyy[i] * fbe_0 * fz_0 + tk_yyzz_yyyy[i] * fe_0 + tk_xyyzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyy[i] * fz_0;

        tk_xxyyzz_yyyz[i] =
            -2.0 * ts_yyzz_yyyz[i] * fbe_0 * fz_0 + tk_yyzz_yyyz[i] * fe_0 + tk_xyyzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyz[i] * fz_0;

        tk_xxyyzz_yyzz[i] =
            -2.0 * ts_yyzz_yyzz[i] * fbe_0 * fz_0 + tk_yyzz_yyzz[i] * fe_0 + tk_xyyzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyzz[i] * fz_0;

        tk_xxyyzz_yzzz[i] =
            -2.0 * ts_yyzz_yzzz[i] * fbe_0 * fz_0 + tk_yyzz_yzzz[i] * fe_0 + tk_xyyzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yzzz[i] * fz_0;

        tk_xxyyzz_zzzz[i] =
            -2.0 * ts_yyzz_zzzz[i] * fbe_0 * fz_0 + tk_yyzz_zzzz[i] * fe_0 + tk_xyyzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_zzzz[i] * fz_0;
    }

    // Set up 195-210 components of targeted buffer : IG

    auto tk_xxyzzz_xxxx = pbuffer.data(idx_kin_ig + 195);

    auto tk_xxyzzz_xxxy = pbuffer.data(idx_kin_ig + 196);

    auto tk_xxyzzz_xxxz = pbuffer.data(idx_kin_ig + 197);

    auto tk_xxyzzz_xxyy = pbuffer.data(idx_kin_ig + 198);

    auto tk_xxyzzz_xxyz = pbuffer.data(idx_kin_ig + 199);

    auto tk_xxyzzz_xxzz = pbuffer.data(idx_kin_ig + 200);

    auto tk_xxyzzz_xyyy = pbuffer.data(idx_kin_ig + 201);

    auto tk_xxyzzz_xyyz = pbuffer.data(idx_kin_ig + 202);

    auto tk_xxyzzz_xyzz = pbuffer.data(idx_kin_ig + 203);

    auto tk_xxyzzz_xzzz = pbuffer.data(idx_kin_ig + 204);

    auto tk_xxyzzz_yyyy = pbuffer.data(idx_kin_ig + 205);

    auto tk_xxyzzz_yyyz = pbuffer.data(idx_kin_ig + 206);

    auto tk_xxyzzz_yyzz = pbuffer.data(idx_kin_ig + 207);

    auto tk_xxyzzz_yzzz = pbuffer.data(idx_kin_ig + 208);

    auto tk_xxyzzz_zzzz = pbuffer.data(idx_kin_ig + 209);

#pragma omp simd aligned(pa_y,               \
                             tk_xxyzzz_xxxx, \
                             tk_xxyzzz_xxxy, \
                             tk_xxyzzz_xxxz, \
                             tk_xxyzzz_xxyy, \
                             tk_xxyzzz_xxyz, \
                             tk_xxyzzz_xxzz, \
                             tk_xxyzzz_xyyy, \
                             tk_xxyzzz_xyyz, \
                             tk_xxyzzz_xyzz, \
                             tk_xxyzzz_xzzz, \
                             tk_xxyzzz_yyyy, \
                             tk_xxyzzz_yyyz, \
                             tk_xxyzzz_yyzz, \
                             tk_xxyzzz_yzzz, \
                             tk_xxyzzz_zzzz, \
                             tk_xxzzz_xxx,   \
                             tk_xxzzz_xxxx,  \
                             tk_xxzzz_xxxy,  \
                             tk_xxzzz_xxxz,  \
                             tk_xxzzz_xxy,   \
                             tk_xxzzz_xxyy,  \
                             tk_xxzzz_xxyz,  \
                             tk_xxzzz_xxz,   \
                             tk_xxzzz_xxzz,  \
                             tk_xxzzz_xyy,   \
                             tk_xxzzz_xyyy,  \
                             tk_xxzzz_xyyz,  \
                             tk_xxzzz_xyz,   \
                             tk_xxzzz_xyzz,  \
                             tk_xxzzz_xzz,   \
                             tk_xxzzz_xzzz,  \
                             tk_xxzzz_yyy,   \
                             tk_xxzzz_yyyy,  \
                             tk_xxzzz_yyyz,  \
                             tk_xxzzz_yyz,   \
                             tk_xxzzz_yyzz,  \
                             tk_xxzzz_yzz,   \
                             tk_xxzzz_yzzz,  \
                             tk_xxzzz_zzz,   \
                             tk_xxzzz_zzzz,  \
                             ts_xxyzzz_xxxx, \
                             ts_xxyzzz_xxxy, \
                             ts_xxyzzz_xxxz, \
                             ts_xxyzzz_xxyy, \
                             ts_xxyzzz_xxyz, \
                             ts_xxyzzz_xxzz, \
                             ts_xxyzzz_xyyy, \
                             ts_xxyzzz_xyyz, \
                             ts_xxyzzz_xyzz, \
                             ts_xxyzzz_xzzz, \
                             ts_xxyzzz_yyyy, \
                             ts_xxyzzz_yyyz, \
                             ts_xxyzzz_yyzz, \
                             ts_xxyzzz_yzzz, \
                             ts_xxyzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzzz_xxxx[i] = tk_xxzzz_xxxx[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxx[i] * fz_0;

        tk_xxyzzz_xxxy[i] = tk_xxzzz_xxx[i] * fe_0 + tk_xxzzz_xxxy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxy[i] * fz_0;

        tk_xxyzzz_xxxz[i] = tk_xxzzz_xxxz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxz[i] * fz_0;

        tk_xxyzzz_xxyy[i] = 2.0 * tk_xxzzz_xxy[i] * fe_0 + tk_xxzzz_xxyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyy[i] * fz_0;

        tk_xxyzzz_xxyz[i] = tk_xxzzz_xxz[i] * fe_0 + tk_xxzzz_xxyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyz[i] * fz_0;

        tk_xxyzzz_xxzz[i] = tk_xxzzz_xxzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxzz[i] * fz_0;

        tk_xxyzzz_xyyy[i] = 3.0 * tk_xxzzz_xyy[i] * fe_0 + tk_xxzzz_xyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyy[i] * fz_0;

        tk_xxyzzz_xyyz[i] = 2.0 * tk_xxzzz_xyz[i] * fe_0 + tk_xxzzz_xyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyz[i] * fz_0;

        tk_xxyzzz_xyzz[i] = tk_xxzzz_xzz[i] * fe_0 + tk_xxzzz_xyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyzz[i] * fz_0;

        tk_xxyzzz_xzzz[i] = tk_xxzzz_xzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xzzz[i] * fz_0;

        tk_xxyzzz_yyyy[i] = 4.0 * tk_xxzzz_yyy[i] * fe_0 + tk_xxzzz_yyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyy[i] * fz_0;

        tk_xxyzzz_yyyz[i] = 3.0 * tk_xxzzz_yyz[i] * fe_0 + tk_xxzzz_yyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyz[i] * fz_0;

        tk_xxyzzz_yyzz[i] = 2.0 * tk_xxzzz_yzz[i] * fe_0 + tk_xxzzz_yyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyzz[i] * fz_0;

        tk_xxyzzz_yzzz[i] = tk_xxzzz_zzz[i] * fe_0 + tk_xxzzz_yzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yzzz[i] * fz_0;

        tk_xxyzzz_zzzz[i] = tk_xxzzz_zzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_zzzz[i] * fz_0;
    }

    // Set up 210-225 components of targeted buffer : IG

    auto tk_xxzzzz_xxxx = pbuffer.data(idx_kin_ig + 210);

    auto tk_xxzzzz_xxxy = pbuffer.data(idx_kin_ig + 211);

    auto tk_xxzzzz_xxxz = pbuffer.data(idx_kin_ig + 212);

    auto tk_xxzzzz_xxyy = pbuffer.data(idx_kin_ig + 213);

    auto tk_xxzzzz_xxyz = pbuffer.data(idx_kin_ig + 214);

    auto tk_xxzzzz_xxzz = pbuffer.data(idx_kin_ig + 215);

    auto tk_xxzzzz_xyyy = pbuffer.data(idx_kin_ig + 216);

    auto tk_xxzzzz_xyyz = pbuffer.data(idx_kin_ig + 217);

    auto tk_xxzzzz_xyzz = pbuffer.data(idx_kin_ig + 218);

    auto tk_xxzzzz_xzzz = pbuffer.data(idx_kin_ig + 219);

    auto tk_xxzzzz_yyyy = pbuffer.data(idx_kin_ig + 220);

    auto tk_xxzzzz_yyyz = pbuffer.data(idx_kin_ig + 221);

    auto tk_xxzzzz_yyzz = pbuffer.data(idx_kin_ig + 222);

    auto tk_xxzzzz_yzzz = pbuffer.data(idx_kin_ig + 223);

    auto tk_xxzzzz_zzzz = pbuffer.data(idx_kin_ig + 224);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xxzz_xxxx,   \
                             tk_xxzz_xxxy,   \
                             tk_xxzz_xxyy,   \
                             tk_xxzz_xyyy,   \
                             tk_xxzzz_xxxx,  \
                             tk_xxzzz_xxxy,  \
                             tk_xxzzz_xxyy,  \
                             tk_xxzzz_xyyy,  \
                             tk_xxzzzz_xxxx, \
                             tk_xxzzzz_xxxy, \
                             tk_xxzzzz_xxxz, \
                             tk_xxzzzz_xxyy, \
                             tk_xxzzzz_xxyz, \
                             tk_xxzzzz_xxzz, \
                             tk_xxzzzz_xyyy, \
                             tk_xxzzzz_xyyz, \
                             tk_xxzzzz_xyzz, \
                             tk_xxzzzz_xzzz, \
                             tk_xxzzzz_yyyy, \
                             tk_xxzzzz_yyyz, \
                             tk_xxzzzz_yyzz, \
                             tk_xxzzzz_yzzz, \
                             tk_xxzzzz_zzzz, \
                             tk_xzzzz_xxxz,  \
                             tk_xzzzz_xxyz,  \
                             tk_xzzzz_xxz,   \
                             tk_xzzzz_xxzz,  \
                             tk_xzzzz_xyyz,  \
                             tk_xzzzz_xyz,   \
                             tk_xzzzz_xyzz,  \
                             tk_xzzzz_xzz,   \
                             tk_xzzzz_xzzz,  \
                             tk_xzzzz_yyyy,  \
                             tk_xzzzz_yyyz,  \
                             tk_xzzzz_yyz,   \
                             tk_xzzzz_yyzz,  \
                             tk_xzzzz_yzz,   \
                             tk_xzzzz_yzzz,  \
                             tk_xzzzz_zzz,   \
                             tk_xzzzz_zzzz,  \
                             tk_zzzz_xxxz,   \
                             tk_zzzz_xxyz,   \
                             tk_zzzz_xxzz,   \
                             tk_zzzz_xyyz,   \
                             tk_zzzz_xyzz,   \
                             tk_zzzz_xzzz,   \
                             tk_zzzz_yyyy,   \
                             tk_zzzz_yyyz,   \
                             tk_zzzz_yyzz,   \
                             tk_zzzz_yzzz,   \
                             tk_zzzz_zzzz,   \
                             ts_xxzz_xxxx,   \
                             ts_xxzz_xxxy,   \
                             ts_xxzz_xxyy,   \
                             ts_xxzz_xyyy,   \
                             ts_xxzzzz_xxxx, \
                             ts_xxzzzz_xxxy, \
                             ts_xxzzzz_xxxz, \
                             ts_xxzzzz_xxyy, \
                             ts_xxzzzz_xxyz, \
                             ts_xxzzzz_xxzz, \
                             ts_xxzzzz_xyyy, \
                             ts_xxzzzz_xyyz, \
                             ts_xxzzzz_xyzz, \
                             ts_xxzzzz_xzzz, \
                             ts_xxzzzz_yyyy, \
                             ts_xxzzzz_yyyz, \
                             ts_xxzzzz_yyzz, \
                             ts_xxzzzz_yzzz, \
                             ts_xxzzzz_zzzz, \
                             ts_zzzz_xxxz,   \
                             ts_zzzz_xxyz,   \
                             ts_zzzz_xxzz,   \
                             ts_zzzz_xyyz,   \
                             ts_zzzz_xyzz,   \
                             ts_zzzz_xzzz,   \
                             ts_zzzz_yyyy,   \
                             ts_zzzz_yyyz,   \
                             ts_zzzz_yyzz,   \
                             ts_zzzz_yzzz,   \
                             ts_zzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzzz_xxxx[i] =
            -6.0 * ts_xxzz_xxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxx[i] * fe_0 + tk_xxzzz_xxxx[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxxx[i] * fz_0;

        tk_xxzzzz_xxxy[i] =
            -6.0 * ts_xxzz_xxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxy[i] * fe_0 + tk_xxzzz_xxxy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxxy[i] * fz_0;

        tk_xxzzzz_xxxz[i] = -2.0 * ts_zzzz_xxxz[i] * fbe_0 * fz_0 + tk_zzzz_xxxz[i] * fe_0 + 3.0 * tk_xzzzz_xxz[i] * fe_0 +
                            tk_xzzzz_xxxz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxz[i] * fz_0;

        tk_xxzzzz_xxyy[i] =
            -6.0 * ts_xxzz_xxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyy[i] * fe_0 + tk_xxzzz_xxyy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxyy[i] * fz_0;

        tk_xxzzzz_xxyz[i] = -2.0 * ts_zzzz_xxyz[i] * fbe_0 * fz_0 + tk_zzzz_xxyz[i] * fe_0 + 2.0 * tk_xzzzz_xyz[i] * fe_0 +
                            tk_xzzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxyz[i] * fz_0;

        tk_xxzzzz_xxzz[i] = -2.0 * ts_zzzz_xxzz[i] * fbe_0 * fz_0 + tk_zzzz_xxzz[i] * fe_0 + 2.0 * tk_xzzzz_xzz[i] * fe_0 +
                            tk_xzzzz_xxzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxzz[i] * fz_0;

        tk_xxzzzz_xyyy[i] =
            -6.0 * ts_xxzz_xyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyy[i] * fe_0 + tk_xxzzz_xyyy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xyyy[i] * fz_0;

        tk_xxzzzz_xyyz[i] = -2.0 * ts_zzzz_xyyz[i] * fbe_0 * fz_0 + tk_zzzz_xyyz[i] * fe_0 + tk_xzzzz_yyz[i] * fe_0 + tk_xzzzz_xyyz[i] * pa_x[i] +
                            2.0 * ts_xxzzzz_xyyz[i] * fz_0;

        tk_xxzzzz_xyzz[i] = -2.0 * ts_zzzz_xyzz[i] * fbe_0 * fz_0 + tk_zzzz_xyzz[i] * fe_0 + tk_xzzzz_yzz[i] * fe_0 + tk_xzzzz_xyzz[i] * pa_x[i] +
                            2.0 * ts_xxzzzz_xyzz[i] * fz_0;

        tk_xxzzzz_xzzz[i] = -2.0 * ts_zzzz_xzzz[i] * fbe_0 * fz_0 + tk_zzzz_xzzz[i] * fe_0 + tk_xzzzz_zzz[i] * fe_0 + tk_xzzzz_xzzz[i] * pa_x[i] +
                            2.0 * ts_xxzzzz_xzzz[i] * fz_0;

        tk_xxzzzz_yyyy[i] =
            -2.0 * ts_zzzz_yyyy[i] * fbe_0 * fz_0 + tk_zzzz_yyyy[i] * fe_0 + tk_xzzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyy[i] * fz_0;

        tk_xxzzzz_yyyz[i] =
            -2.0 * ts_zzzz_yyyz[i] * fbe_0 * fz_0 + tk_zzzz_yyyz[i] * fe_0 + tk_xzzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyz[i] * fz_0;

        tk_xxzzzz_yyzz[i] =
            -2.0 * ts_zzzz_yyzz[i] * fbe_0 * fz_0 + tk_zzzz_yyzz[i] * fe_0 + tk_xzzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyzz[i] * fz_0;

        tk_xxzzzz_yzzz[i] =
            -2.0 * ts_zzzz_yzzz[i] * fbe_0 * fz_0 + tk_zzzz_yzzz[i] * fe_0 + tk_xzzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yzzz[i] * fz_0;

        tk_xxzzzz_zzzz[i] =
            -2.0 * ts_zzzz_zzzz[i] * fbe_0 * fz_0 + tk_zzzz_zzzz[i] * fe_0 + tk_xzzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_zzzz[i] * fz_0;
    }

    // Set up 225-240 components of targeted buffer : IG

    auto tk_xyyyyy_xxxx = pbuffer.data(idx_kin_ig + 225);

    auto tk_xyyyyy_xxxy = pbuffer.data(idx_kin_ig + 226);

    auto tk_xyyyyy_xxxz = pbuffer.data(idx_kin_ig + 227);

    auto tk_xyyyyy_xxyy = pbuffer.data(idx_kin_ig + 228);

    auto tk_xyyyyy_xxyz = pbuffer.data(idx_kin_ig + 229);

    auto tk_xyyyyy_xxzz = pbuffer.data(idx_kin_ig + 230);

    auto tk_xyyyyy_xyyy = pbuffer.data(idx_kin_ig + 231);

    auto tk_xyyyyy_xyyz = pbuffer.data(idx_kin_ig + 232);

    auto tk_xyyyyy_xyzz = pbuffer.data(idx_kin_ig + 233);

    auto tk_xyyyyy_xzzz = pbuffer.data(idx_kin_ig + 234);

    auto tk_xyyyyy_yyyy = pbuffer.data(idx_kin_ig + 235);

    auto tk_xyyyyy_yyyz = pbuffer.data(idx_kin_ig + 236);

    auto tk_xyyyyy_yyzz = pbuffer.data(idx_kin_ig + 237);

    auto tk_xyyyyy_yzzz = pbuffer.data(idx_kin_ig + 238);

    auto tk_xyyyyy_zzzz = pbuffer.data(idx_kin_ig + 239);

#pragma omp simd aligned(pa_x,               \
                             tk_xyyyyy_xxxx, \
                             tk_xyyyyy_xxxy, \
                             tk_xyyyyy_xxxz, \
                             tk_xyyyyy_xxyy, \
                             tk_xyyyyy_xxyz, \
                             tk_xyyyyy_xxzz, \
                             tk_xyyyyy_xyyy, \
                             tk_xyyyyy_xyyz, \
                             tk_xyyyyy_xyzz, \
                             tk_xyyyyy_xzzz, \
                             tk_xyyyyy_yyyy, \
                             tk_xyyyyy_yyyz, \
                             tk_xyyyyy_yyzz, \
                             tk_xyyyyy_yzzz, \
                             tk_xyyyyy_zzzz, \
                             tk_yyyyy_xxx,   \
                             tk_yyyyy_xxxx,  \
                             tk_yyyyy_xxxy,  \
                             tk_yyyyy_xxxz,  \
                             tk_yyyyy_xxy,   \
                             tk_yyyyy_xxyy,  \
                             tk_yyyyy_xxyz,  \
                             tk_yyyyy_xxz,   \
                             tk_yyyyy_xxzz,  \
                             tk_yyyyy_xyy,   \
                             tk_yyyyy_xyyy,  \
                             tk_yyyyy_xyyz,  \
                             tk_yyyyy_xyz,   \
                             tk_yyyyy_xyzz,  \
                             tk_yyyyy_xzz,   \
                             tk_yyyyy_xzzz,  \
                             tk_yyyyy_yyy,   \
                             tk_yyyyy_yyyy,  \
                             tk_yyyyy_yyyz,  \
                             tk_yyyyy_yyz,   \
                             tk_yyyyy_yyzz,  \
                             tk_yyyyy_yzz,   \
                             tk_yyyyy_yzzz,  \
                             tk_yyyyy_zzz,   \
                             tk_yyyyy_zzzz,  \
                             ts_xyyyyy_xxxx, \
                             ts_xyyyyy_xxxy, \
                             ts_xyyyyy_xxxz, \
                             ts_xyyyyy_xxyy, \
                             ts_xyyyyy_xxyz, \
                             ts_xyyyyy_xxzz, \
                             ts_xyyyyy_xyyy, \
                             ts_xyyyyy_xyyz, \
                             ts_xyyyyy_xyzz, \
                             ts_xyyyyy_xzzz, \
                             ts_xyyyyy_yyyy, \
                             ts_xyyyyy_yyyz, \
                             ts_xyyyyy_yyzz, \
                             ts_xyyyyy_yzzz, \
                             ts_xyyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyy_xxxx[i] = 4.0 * tk_yyyyy_xxx[i] * fe_0 + tk_yyyyy_xxxx[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxx[i] * fz_0;

        tk_xyyyyy_xxxy[i] = 3.0 * tk_yyyyy_xxy[i] * fe_0 + tk_yyyyy_xxxy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxy[i] * fz_0;

        tk_xyyyyy_xxxz[i] = 3.0 * tk_yyyyy_xxz[i] * fe_0 + tk_yyyyy_xxxz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxz[i] * fz_0;

        tk_xyyyyy_xxyy[i] = 2.0 * tk_yyyyy_xyy[i] * fe_0 + tk_yyyyy_xxyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyy[i] * fz_0;

        tk_xyyyyy_xxyz[i] = 2.0 * tk_yyyyy_xyz[i] * fe_0 + tk_yyyyy_xxyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyz[i] * fz_0;

        tk_xyyyyy_xxzz[i] = 2.0 * tk_yyyyy_xzz[i] * fe_0 + tk_yyyyy_xxzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxzz[i] * fz_0;

        tk_xyyyyy_xyyy[i] = tk_yyyyy_yyy[i] * fe_0 + tk_yyyyy_xyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyy[i] * fz_0;

        tk_xyyyyy_xyyz[i] = tk_yyyyy_yyz[i] * fe_0 + tk_yyyyy_xyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyz[i] * fz_0;

        tk_xyyyyy_xyzz[i] = tk_yyyyy_yzz[i] * fe_0 + tk_yyyyy_xyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyzz[i] * fz_0;

        tk_xyyyyy_xzzz[i] = tk_yyyyy_zzz[i] * fe_0 + tk_yyyyy_xzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xzzz[i] * fz_0;

        tk_xyyyyy_yyyy[i] = tk_yyyyy_yyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyy[i] * fz_0;

        tk_xyyyyy_yyyz[i] = tk_yyyyy_yyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyz[i] * fz_0;

        tk_xyyyyy_yyzz[i] = tk_yyyyy_yyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyzz[i] * fz_0;

        tk_xyyyyy_yzzz[i] = tk_yyyyy_yzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yzzz[i] * fz_0;

        tk_xyyyyy_zzzz[i] = tk_yyyyy_zzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_zzzz[i] * fz_0;
    }

    // Set up 240-255 components of targeted buffer : IG

    auto tk_xyyyyz_xxxx = pbuffer.data(idx_kin_ig + 240);

    auto tk_xyyyyz_xxxy = pbuffer.data(idx_kin_ig + 241);

    auto tk_xyyyyz_xxxz = pbuffer.data(idx_kin_ig + 242);

    auto tk_xyyyyz_xxyy = pbuffer.data(idx_kin_ig + 243);

    auto tk_xyyyyz_xxyz = pbuffer.data(idx_kin_ig + 244);

    auto tk_xyyyyz_xxzz = pbuffer.data(idx_kin_ig + 245);

    auto tk_xyyyyz_xyyy = pbuffer.data(idx_kin_ig + 246);

    auto tk_xyyyyz_xyyz = pbuffer.data(idx_kin_ig + 247);

    auto tk_xyyyyz_xyzz = pbuffer.data(idx_kin_ig + 248);

    auto tk_xyyyyz_xzzz = pbuffer.data(idx_kin_ig + 249);

    auto tk_xyyyyz_yyyy = pbuffer.data(idx_kin_ig + 250);

    auto tk_xyyyyz_yyyz = pbuffer.data(idx_kin_ig + 251);

    auto tk_xyyyyz_yyzz = pbuffer.data(idx_kin_ig + 252);

    auto tk_xyyyyz_yzzz = pbuffer.data(idx_kin_ig + 253);

    auto tk_xyyyyz_zzzz = pbuffer.data(idx_kin_ig + 254);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xyyyy_xxxx,  \
                             tk_xyyyy_xxxy,  \
                             tk_xyyyy_xxyy,  \
                             tk_xyyyy_xyyy,  \
                             tk_xyyyyz_xxxx, \
                             tk_xyyyyz_xxxy, \
                             tk_xyyyyz_xxxz, \
                             tk_xyyyyz_xxyy, \
                             tk_xyyyyz_xxyz, \
                             tk_xyyyyz_xxzz, \
                             tk_xyyyyz_xyyy, \
                             tk_xyyyyz_xyyz, \
                             tk_xyyyyz_xyzz, \
                             tk_xyyyyz_xzzz, \
                             tk_xyyyyz_yyyy, \
                             tk_xyyyyz_yyyz, \
                             tk_xyyyyz_yyzz, \
                             tk_xyyyyz_yzzz, \
                             tk_xyyyyz_zzzz, \
                             tk_yyyyz_xxxz,  \
                             tk_yyyyz_xxyz,  \
                             tk_yyyyz_xxz,   \
                             tk_yyyyz_xxzz,  \
                             tk_yyyyz_xyyz,  \
                             tk_yyyyz_xyz,   \
                             tk_yyyyz_xyzz,  \
                             tk_yyyyz_xzz,   \
                             tk_yyyyz_xzzz,  \
                             tk_yyyyz_yyyy,  \
                             tk_yyyyz_yyyz,  \
                             tk_yyyyz_yyz,   \
                             tk_yyyyz_yyzz,  \
                             tk_yyyyz_yzz,   \
                             tk_yyyyz_yzzz,  \
                             tk_yyyyz_zzz,   \
                             tk_yyyyz_zzzz,  \
                             ts_xyyyyz_xxxx, \
                             ts_xyyyyz_xxxy, \
                             ts_xyyyyz_xxxz, \
                             ts_xyyyyz_xxyy, \
                             ts_xyyyyz_xxyz, \
                             ts_xyyyyz_xxzz, \
                             ts_xyyyyz_xyyy, \
                             ts_xyyyyz_xyyz, \
                             ts_xyyyyz_xyzz, \
                             ts_xyyyyz_xzzz, \
                             ts_xyyyyz_yyyy, \
                             ts_xyyyyz_yyyz, \
                             ts_xyyyyz_yyzz, \
                             ts_xyyyyz_yzzz, \
                             ts_xyyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyz_xxxx[i] = tk_xyyyy_xxxx[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxx[i] * fz_0;

        tk_xyyyyz_xxxy[i] = tk_xyyyy_xxxy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxy[i] * fz_0;

        tk_xyyyyz_xxxz[i] = 3.0 * tk_yyyyz_xxz[i] * fe_0 + tk_yyyyz_xxxz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxz[i] * fz_0;

        tk_xyyyyz_xxyy[i] = tk_xyyyy_xxyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxyy[i] * fz_0;

        tk_xyyyyz_xxyz[i] = 2.0 * tk_yyyyz_xyz[i] * fe_0 + tk_yyyyz_xxyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxyz[i] * fz_0;

        tk_xyyyyz_xxzz[i] = 2.0 * tk_yyyyz_xzz[i] * fe_0 + tk_yyyyz_xxzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxzz[i] * fz_0;

        tk_xyyyyz_xyyy[i] = tk_xyyyy_xyyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xyyy[i] * fz_0;

        tk_xyyyyz_xyyz[i] = tk_yyyyz_yyz[i] * fe_0 + tk_yyyyz_xyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyyz[i] * fz_0;

        tk_xyyyyz_xyzz[i] = tk_yyyyz_yzz[i] * fe_0 + tk_yyyyz_xyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyzz[i] * fz_0;

        tk_xyyyyz_xzzz[i] = tk_yyyyz_zzz[i] * fe_0 + tk_yyyyz_xzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xzzz[i] * fz_0;

        tk_xyyyyz_yyyy[i] = tk_yyyyz_yyyy[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyy[i] * fz_0;

        tk_xyyyyz_yyyz[i] = tk_yyyyz_yyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyz[i] * fz_0;

        tk_xyyyyz_yyzz[i] = tk_yyyyz_yyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyzz[i] * fz_0;

        tk_xyyyyz_yzzz[i] = tk_yyyyz_yzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yzzz[i] * fz_0;

        tk_xyyyyz_zzzz[i] = tk_yyyyz_zzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_zzzz[i] * fz_0;
    }

    // Set up 255-270 components of targeted buffer : IG

    auto tk_xyyyzz_xxxx = pbuffer.data(idx_kin_ig + 255);

    auto tk_xyyyzz_xxxy = pbuffer.data(idx_kin_ig + 256);

    auto tk_xyyyzz_xxxz = pbuffer.data(idx_kin_ig + 257);

    auto tk_xyyyzz_xxyy = pbuffer.data(idx_kin_ig + 258);

    auto tk_xyyyzz_xxyz = pbuffer.data(idx_kin_ig + 259);

    auto tk_xyyyzz_xxzz = pbuffer.data(idx_kin_ig + 260);

    auto tk_xyyyzz_xyyy = pbuffer.data(idx_kin_ig + 261);

    auto tk_xyyyzz_xyyz = pbuffer.data(idx_kin_ig + 262);

    auto tk_xyyyzz_xyzz = pbuffer.data(idx_kin_ig + 263);

    auto tk_xyyyzz_xzzz = pbuffer.data(idx_kin_ig + 264);

    auto tk_xyyyzz_yyyy = pbuffer.data(idx_kin_ig + 265);

    auto tk_xyyyzz_yyyz = pbuffer.data(idx_kin_ig + 266);

    auto tk_xyyyzz_yyzz = pbuffer.data(idx_kin_ig + 267);

    auto tk_xyyyzz_yzzz = pbuffer.data(idx_kin_ig + 268);

    auto tk_xyyyzz_zzzz = pbuffer.data(idx_kin_ig + 269);

#pragma omp simd aligned(pa_x,               \
                             tk_xyyyzz_xxxx, \
                             tk_xyyyzz_xxxy, \
                             tk_xyyyzz_xxxz, \
                             tk_xyyyzz_xxyy, \
                             tk_xyyyzz_xxyz, \
                             tk_xyyyzz_xxzz, \
                             tk_xyyyzz_xyyy, \
                             tk_xyyyzz_xyyz, \
                             tk_xyyyzz_xyzz, \
                             tk_xyyyzz_xzzz, \
                             tk_xyyyzz_yyyy, \
                             tk_xyyyzz_yyyz, \
                             tk_xyyyzz_yyzz, \
                             tk_xyyyzz_yzzz, \
                             tk_xyyyzz_zzzz, \
                             tk_yyyzz_xxx,   \
                             tk_yyyzz_xxxx,  \
                             tk_yyyzz_xxxy,  \
                             tk_yyyzz_xxxz,  \
                             tk_yyyzz_xxy,   \
                             tk_yyyzz_xxyy,  \
                             tk_yyyzz_xxyz,  \
                             tk_yyyzz_xxz,   \
                             tk_yyyzz_xxzz,  \
                             tk_yyyzz_xyy,   \
                             tk_yyyzz_xyyy,  \
                             tk_yyyzz_xyyz,  \
                             tk_yyyzz_xyz,   \
                             tk_yyyzz_xyzz,  \
                             tk_yyyzz_xzz,   \
                             tk_yyyzz_xzzz,  \
                             tk_yyyzz_yyy,   \
                             tk_yyyzz_yyyy,  \
                             tk_yyyzz_yyyz,  \
                             tk_yyyzz_yyz,   \
                             tk_yyyzz_yyzz,  \
                             tk_yyyzz_yzz,   \
                             tk_yyyzz_yzzz,  \
                             tk_yyyzz_zzz,   \
                             tk_yyyzz_zzzz,  \
                             ts_xyyyzz_xxxx, \
                             ts_xyyyzz_xxxy, \
                             ts_xyyyzz_xxxz, \
                             ts_xyyyzz_xxyy, \
                             ts_xyyyzz_xxyz, \
                             ts_xyyyzz_xxzz, \
                             ts_xyyyzz_xyyy, \
                             ts_xyyyzz_xyyz, \
                             ts_xyyyzz_xyzz, \
                             ts_xyyyzz_xzzz, \
                             ts_xyyyzz_yyyy, \
                             ts_xyyyzz_yyyz, \
                             ts_xyyyzz_yyzz, \
                             ts_xyyyzz_yzzz, \
                             ts_xyyyzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyzz_xxxx[i] = 4.0 * tk_yyyzz_xxx[i] * fe_0 + tk_yyyzz_xxxx[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxx[i] * fz_0;

        tk_xyyyzz_xxxy[i] = 3.0 * tk_yyyzz_xxy[i] * fe_0 + tk_yyyzz_xxxy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxy[i] * fz_0;

        tk_xyyyzz_xxxz[i] = 3.0 * tk_yyyzz_xxz[i] * fe_0 + tk_yyyzz_xxxz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxz[i] * fz_0;

        tk_xyyyzz_xxyy[i] = 2.0 * tk_yyyzz_xyy[i] * fe_0 + tk_yyyzz_xxyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyy[i] * fz_0;

        tk_xyyyzz_xxyz[i] = 2.0 * tk_yyyzz_xyz[i] * fe_0 + tk_yyyzz_xxyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyz[i] * fz_0;

        tk_xyyyzz_xxzz[i] = 2.0 * tk_yyyzz_xzz[i] * fe_0 + tk_yyyzz_xxzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxzz[i] * fz_0;

        tk_xyyyzz_xyyy[i] = tk_yyyzz_yyy[i] * fe_0 + tk_yyyzz_xyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyy[i] * fz_0;

        tk_xyyyzz_xyyz[i] = tk_yyyzz_yyz[i] * fe_0 + tk_yyyzz_xyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyz[i] * fz_0;

        tk_xyyyzz_xyzz[i] = tk_yyyzz_yzz[i] * fe_0 + tk_yyyzz_xyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyzz[i] * fz_0;

        tk_xyyyzz_xzzz[i] = tk_yyyzz_zzz[i] * fe_0 + tk_yyyzz_xzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xzzz[i] * fz_0;

        tk_xyyyzz_yyyy[i] = tk_yyyzz_yyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyy[i] * fz_0;

        tk_xyyyzz_yyyz[i] = tk_yyyzz_yyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyz[i] * fz_0;

        tk_xyyyzz_yyzz[i] = tk_yyyzz_yyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyzz[i] * fz_0;

        tk_xyyyzz_yzzz[i] = tk_yyyzz_yzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yzzz[i] * fz_0;

        tk_xyyyzz_zzzz[i] = tk_yyyzz_zzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_zzzz[i] * fz_0;
    }

    // Set up 270-285 components of targeted buffer : IG

    auto tk_xyyzzz_xxxx = pbuffer.data(idx_kin_ig + 270);

    auto tk_xyyzzz_xxxy = pbuffer.data(idx_kin_ig + 271);

    auto tk_xyyzzz_xxxz = pbuffer.data(idx_kin_ig + 272);

    auto tk_xyyzzz_xxyy = pbuffer.data(idx_kin_ig + 273);

    auto tk_xyyzzz_xxyz = pbuffer.data(idx_kin_ig + 274);

    auto tk_xyyzzz_xxzz = pbuffer.data(idx_kin_ig + 275);

    auto tk_xyyzzz_xyyy = pbuffer.data(idx_kin_ig + 276);

    auto tk_xyyzzz_xyyz = pbuffer.data(idx_kin_ig + 277);

    auto tk_xyyzzz_xyzz = pbuffer.data(idx_kin_ig + 278);

    auto tk_xyyzzz_xzzz = pbuffer.data(idx_kin_ig + 279);

    auto tk_xyyzzz_yyyy = pbuffer.data(idx_kin_ig + 280);

    auto tk_xyyzzz_yyyz = pbuffer.data(idx_kin_ig + 281);

    auto tk_xyyzzz_yyzz = pbuffer.data(idx_kin_ig + 282);

    auto tk_xyyzzz_yzzz = pbuffer.data(idx_kin_ig + 283);

    auto tk_xyyzzz_zzzz = pbuffer.data(idx_kin_ig + 284);

#pragma omp simd aligned(pa_x,               \
                             tk_xyyzzz_xxxx, \
                             tk_xyyzzz_xxxy, \
                             tk_xyyzzz_xxxz, \
                             tk_xyyzzz_xxyy, \
                             tk_xyyzzz_xxyz, \
                             tk_xyyzzz_xxzz, \
                             tk_xyyzzz_xyyy, \
                             tk_xyyzzz_xyyz, \
                             tk_xyyzzz_xyzz, \
                             tk_xyyzzz_xzzz, \
                             tk_xyyzzz_yyyy, \
                             tk_xyyzzz_yyyz, \
                             tk_xyyzzz_yyzz, \
                             tk_xyyzzz_yzzz, \
                             tk_xyyzzz_zzzz, \
                             tk_yyzzz_xxx,   \
                             tk_yyzzz_xxxx,  \
                             tk_yyzzz_xxxy,  \
                             tk_yyzzz_xxxz,  \
                             tk_yyzzz_xxy,   \
                             tk_yyzzz_xxyy,  \
                             tk_yyzzz_xxyz,  \
                             tk_yyzzz_xxz,   \
                             tk_yyzzz_xxzz,  \
                             tk_yyzzz_xyy,   \
                             tk_yyzzz_xyyy,  \
                             tk_yyzzz_xyyz,  \
                             tk_yyzzz_xyz,   \
                             tk_yyzzz_xyzz,  \
                             tk_yyzzz_xzz,   \
                             tk_yyzzz_xzzz,  \
                             tk_yyzzz_yyy,   \
                             tk_yyzzz_yyyy,  \
                             tk_yyzzz_yyyz,  \
                             tk_yyzzz_yyz,   \
                             tk_yyzzz_yyzz,  \
                             tk_yyzzz_yzz,   \
                             tk_yyzzz_yzzz,  \
                             tk_yyzzz_zzz,   \
                             tk_yyzzz_zzzz,  \
                             ts_xyyzzz_xxxx, \
                             ts_xyyzzz_xxxy, \
                             ts_xyyzzz_xxxz, \
                             ts_xyyzzz_xxyy, \
                             ts_xyyzzz_xxyz, \
                             ts_xyyzzz_xxzz, \
                             ts_xyyzzz_xyyy, \
                             ts_xyyzzz_xyyz, \
                             ts_xyyzzz_xyzz, \
                             ts_xyyzzz_xzzz, \
                             ts_xyyzzz_yyyy, \
                             ts_xyyzzz_yyyz, \
                             ts_xyyzzz_yyzz, \
                             ts_xyyzzz_yzzz, \
                             ts_xyyzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzzz_xxxx[i] = 4.0 * tk_yyzzz_xxx[i] * fe_0 + tk_yyzzz_xxxx[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxx[i] * fz_0;

        tk_xyyzzz_xxxy[i] = 3.0 * tk_yyzzz_xxy[i] * fe_0 + tk_yyzzz_xxxy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxy[i] * fz_0;

        tk_xyyzzz_xxxz[i] = 3.0 * tk_yyzzz_xxz[i] * fe_0 + tk_yyzzz_xxxz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxz[i] * fz_0;

        tk_xyyzzz_xxyy[i] = 2.0 * tk_yyzzz_xyy[i] * fe_0 + tk_yyzzz_xxyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyy[i] * fz_0;

        tk_xyyzzz_xxyz[i] = 2.0 * tk_yyzzz_xyz[i] * fe_0 + tk_yyzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyz[i] * fz_0;

        tk_xyyzzz_xxzz[i] = 2.0 * tk_yyzzz_xzz[i] * fe_0 + tk_yyzzz_xxzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxzz[i] * fz_0;

        tk_xyyzzz_xyyy[i] = tk_yyzzz_yyy[i] * fe_0 + tk_yyzzz_xyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyy[i] * fz_0;

        tk_xyyzzz_xyyz[i] = tk_yyzzz_yyz[i] * fe_0 + tk_yyzzz_xyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyz[i] * fz_0;

        tk_xyyzzz_xyzz[i] = tk_yyzzz_yzz[i] * fe_0 + tk_yyzzz_xyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyzz[i] * fz_0;

        tk_xyyzzz_xzzz[i] = tk_yyzzz_zzz[i] * fe_0 + tk_yyzzz_xzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xzzz[i] * fz_0;

        tk_xyyzzz_yyyy[i] = tk_yyzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyy[i] * fz_0;

        tk_xyyzzz_yyyz[i] = tk_yyzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyz[i] * fz_0;

        tk_xyyzzz_yyzz[i] = tk_yyzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyzz[i] * fz_0;

        tk_xyyzzz_yzzz[i] = tk_yyzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yzzz[i] * fz_0;

        tk_xyyzzz_zzzz[i] = tk_yyzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_zzzz[i] * fz_0;
    }

    // Set up 285-300 components of targeted buffer : IG

    auto tk_xyzzzz_xxxx = pbuffer.data(idx_kin_ig + 285);

    auto tk_xyzzzz_xxxy = pbuffer.data(idx_kin_ig + 286);

    auto tk_xyzzzz_xxxz = pbuffer.data(idx_kin_ig + 287);

    auto tk_xyzzzz_xxyy = pbuffer.data(idx_kin_ig + 288);

    auto tk_xyzzzz_xxyz = pbuffer.data(idx_kin_ig + 289);

    auto tk_xyzzzz_xxzz = pbuffer.data(idx_kin_ig + 290);

    auto tk_xyzzzz_xyyy = pbuffer.data(idx_kin_ig + 291);

    auto tk_xyzzzz_xyyz = pbuffer.data(idx_kin_ig + 292);

    auto tk_xyzzzz_xyzz = pbuffer.data(idx_kin_ig + 293);

    auto tk_xyzzzz_xzzz = pbuffer.data(idx_kin_ig + 294);

    auto tk_xyzzzz_yyyy = pbuffer.data(idx_kin_ig + 295);

    auto tk_xyzzzz_yyyz = pbuffer.data(idx_kin_ig + 296);

    auto tk_xyzzzz_yyzz = pbuffer.data(idx_kin_ig + 297);

    auto tk_xyzzzz_yzzz = pbuffer.data(idx_kin_ig + 298);

    auto tk_xyzzzz_zzzz = pbuffer.data(idx_kin_ig + 299);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xyzzzz_xxxx, \
                             tk_xyzzzz_xxxy, \
                             tk_xyzzzz_xxxz, \
                             tk_xyzzzz_xxyy, \
                             tk_xyzzzz_xxyz, \
                             tk_xyzzzz_xxzz, \
                             tk_xyzzzz_xyyy, \
                             tk_xyzzzz_xyyz, \
                             tk_xyzzzz_xyzz, \
                             tk_xyzzzz_xzzz, \
                             tk_xyzzzz_yyyy, \
                             tk_xyzzzz_yyyz, \
                             tk_xyzzzz_yyzz, \
                             tk_xyzzzz_yzzz, \
                             tk_xyzzzz_zzzz, \
                             tk_xzzzz_xxxx,  \
                             tk_xzzzz_xxxz,  \
                             tk_xzzzz_xxzz,  \
                             tk_xzzzz_xzzz,  \
                             tk_yzzzz_xxxy,  \
                             tk_yzzzz_xxy,   \
                             tk_yzzzz_xxyy,  \
                             tk_yzzzz_xxyz,  \
                             tk_yzzzz_xyy,   \
                             tk_yzzzz_xyyy,  \
                             tk_yzzzz_xyyz,  \
                             tk_yzzzz_xyz,   \
                             tk_yzzzz_xyzz,  \
                             tk_yzzzz_yyy,   \
                             tk_yzzzz_yyyy,  \
                             tk_yzzzz_yyyz,  \
                             tk_yzzzz_yyz,   \
                             tk_yzzzz_yyzz,  \
                             tk_yzzzz_yzz,   \
                             tk_yzzzz_yzzz,  \
                             tk_yzzzz_zzzz,  \
                             ts_xyzzzz_xxxx, \
                             ts_xyzzzz_xxxy, \
                             ts_xyzzzz_xxxz, \
                             ts_xyzzzz_xxyy, \
                             ts_xyzzzz_xxyz, \
                             ts_xyzzzz_xxzz, \
                             ts_xyzzzz_xyyy, \
                             ts_xyzzzz_xyyz, \
                             ts_xyzzzz_xyzz, \
                             ts_xyzzzz_xzzz, \
                             ts_xyzzzz_yyyy, \
                             ts_xyzzzz_yyyz, \
                             ts_xyzzzz_yyzz, \
                             ts_xyzzzz_yzzz, \
                             ts_xyzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzzz_xxxx[i] = tk_xzzzz_xxxx[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxx[i] * fz_0;

        tk_xyzzzz_xxxy[i] = 3.0 * tk_yzzzz_xxy[i] * fe_0 + tk_yzzzz_xxxy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxy[i] * fz_0;

        tk_xyzzzz_xxxz[i] = tk_xzzzz_xxxz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxz[i] * fz_0;

        tk_xyzzzz_xxyy[i] = 2.0 * tk_yzzzz_xyy[i] * fe_0 + tk_yzzzz_xxyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyy[i] * fz_0;

        tk_xyzzzz_xxyz[i] = 2.0 * tk_yzzzz_xyz[i] * fe_0 + tk_yzzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyz[i] * fz_0;

        tk_xyzzzz_xxzz[i] = tk_xzzzz_xxzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxzz[i] * fz_0;

        tk_xyzzzz_xyyy[i] = tk_yzzzz_yyy[i] * fe_0 + tk_yzzzz_xyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyy[i] * fz_0;

        tk_xyzzzz_xyyz[i] = tk_yzzzz_yyz[i] * fe_0 + tk_yzzzz_xyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyz[i] * fz_0;

        tk_xyzzzz_xyzz[i] = tk_yzzzz_yzz[i] * fe_0 + tk_yzzzz_xyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyzz[i] * fz_0;

        tk_xyzzzz_xzzz[i] = tk_xzzzz_xzzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xzzz[i] * fz_0;

        tk_xyzzzz_yyyy[i] = tk_yzzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyy[i] * fz_0;

        tk_xyzzzz_yyyz[i] = tk_yzzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyz[i] * fz_0;

        tk_xyzzzz_yyzz[i] = tk_yzzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyzz[i] * fz_0;

        tk_xyzzzz_yzzz[i] = tk_yzzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yzzz[i] * fz_0;

        tk_xyzzzz_zzzz[i] = tk_yzzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_zzzz[i] * fz_0;
    }

    // Set up 300-315 components of targeted buffer : IG

    auto tk_xzzzzz_xxxx = pbuffer.data(idx_kin_ig + 300);

    auto tk_xzzzzz_xxxy = pbuffer.data(idx_kin_ig + 301);

    auto tk_xzzzzz_xxxz = pbuffer.data(idx_kin_ig + 302);

    auto tk_xzzzzz_xxyy = pbuffer.data(idx_kin_ig + 303);

    auto tk_xzzzzz_xxyz = pbuffer.data(idx_kin_ig + 304);

    auto tk_xzzzzz_xxzz = pbuffer.data(idx_kin_ig + 305);

    auto tk_xzzzzz_xyyy = pbuffer.data(idx_kin_ig + 306);

    auto tk_xzzzzz_xyyz = pbuffer.data(idx_kin_ig + 307);

    auto tk_xzzzzz_xyzz = pbuffer.data(idx_kin_ig + 308);

    auto tk_xzzzzz_xzzz = pbuffer.data(idx_kin_ig + 309);

    auto tk_xzzzzz_yyyy = pbuffer.data(idx_kin_ig + 310);

    auto tk_xzzzzz_yyyz = pbuffer.data(idx_kin_ig + 311);

    auto tk_xzzzzz_yyzz = pbuffer.data(idx_kin_ig + 312);

    auto tk_xzzzzz_yzzz = pbuffer.data(idx_kin_ig + 313);

    auto tk_xzzzzz_zzzz = pbuffer.data(idx_kin_ig + 314);

#pragma omp simd aligned(pa_x,               \
                             tk_xzzzzz_xxxx, \
                             tk_xzzzzz_xxxy, \
                             tk_xzzzzz_xxxz, \
                             tk_xzzzzz_xxyy, \
                             tk_xzzzzz_xxyz, \
                             tk_xzzzzz_xxzz, \
                             tk_xzzzzz_xyyy, \
                             tk_xzzzzz_xyyz, \
                             tk_xzzzzz_xyzz, \
                             tk_xzzzzz_xzzz, \
                             tk_xzzzzz_yyyy, \
                             tk_xzzzzz_yyyz, \
                             tk_xzzzzz_yyzz, \
                             tk_xzzzzz_yzzz, \
                             tk_xzzzzz_zzzz, \
                             tk_zzzzz_xxx,   \
                             tk_zzzzz_xxxx,  \
                             tk_zzzzz_xxxy,  \
                             tk_zzzzz_xxxz,  \
                             tk_zzzzz_xxy,   \
                             tk_zzzzz_xxyy,  \
                             tk_zzzzz_xxyz,  \
                             tk_zzzzz_xxz,   \
                             tk_zzzzz_xxzz,  \
                             tk_zzzzz_xyy,   \
                             tk_zzzzz_xyyy,  \
                             tk_zzzzz_xyyz,  \
                             tk_zzzzz_xyz,   \
                             tk_zzzzz_xyzz,  \
                             tk_zzzzz_xzz,   \
                             tk_zzzzz_xzzz,  \
                             tk_zzzzz_yyy,   \
                             tk_zzzzz_yyyy,  \
                             tk_zzzzz_yyyz,  \
                             tk_zzzzz_yyz,   \
                             tk_zzzzz_yyzz,  \
                             tk_zzzzz_yzz,   \
                             tk_zzzzz_yzzz,  \
                             tk_zzzzz_zzz,   \
                             tk_zzzzz_zzzz,  \
                             ts_xzzzzz_xxxx, \
                             ts_xzzzzz_xxxy, \
                             ts_xzzzzz_xxxz, \
                             ts_xzzzzz_xxyy, \
                             ts_xzzzzz_xxyz, \
                             ts_xzzzzz_xxzz, \
                             ts_xzzzzz_xyyy, \
                             ts_xzzzzz_xyyz, \
                             ts_xzzzzz_xyzz, \
                             ts_xzzzzz_xzzz, \
                             ts_xzzzzz_yyyy, \
                             ts_xzzzzz_yyyz, \
                             ts_xzzzzz_yyzz, \
                             ts_xzzzzz_yzzz, \
                             ts_xzzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzzz_xxxx[i] = 4.0 * tk_zzzzz_xxx[i] * fe_0 + tk_zzzzz_xxxx[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxx[i] * fz_0;

        tk_xzzzzz_xxxy[i] = 3.0 * tk_zzzzz_xxy[i] * fe_0 + tk_zzzzz_xxxy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxy[i] * fz_0;

        tk_xzzzzz_xxxz[i] = 3.0 * tk_zzzzz_xxz[i] * fe_0 + tk_zzzzz_xxxz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxz[i] * fz_0;

        tk_xzzzzz_xxyy[i] = 2.0 * tk_zzzzz_xyy[i] * fe_0 + tk_zzzzz_xxyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyy[i] * fz_0;

        tk_xzzzzz_xxyz[i] = 2.0 * tk_zzzzz_xyz[i] * fe_0 + tk_zzzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyz[i] * fz_0;

        tk_xzzzzz_xxzz[i] = 2.0 * tk_zzzzz_xzz[i] * fe_0 + tk_zzzzz_xxzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxzz[i] * fz_0;

        tk_xzzzzz_xyyy[i] = tk_zzzzz_yyy[i] * fe_0 + tk_zzzzz_xyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyy[i] * fz_0;

        tk_xzzzzz_xyyz[i] = tk_zzzzz_yyz[i] * fe_0 + tk_zzzzz_xyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyz[i] * fz_0;

        tk_xzzzzz_xyzz[i] = tk_zzzzz_yzz[i] * fe_0 + tk_zzzzz_xyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyzz[i] * fz_0;

        tk_xzzzzz_xzzz[i] = tk_zzzzz_zzz[i] * fe_0 + tk_zzzzz_xzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xzzz[i] * fz_0;

        tk_xzzzzz_yyyy[i] = tk_zzzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyy[i] * fz_0;

        tk_xzzzzz_yyyz[i] = tk_zzzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyz[i] * fz_0;

        tk_xzzzzz_yyzz[i] = tk_zzzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyzz[i] * fz_0;

        tk_xzzzzz_yzzz[i] = tk_zzzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yzzz[i] * fz_0;

        tk_xzzzzz_zzzz[i] = tk_zzzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_zzzz[i] * fz_0;
    }

    // Set up 315-330 components of targeted buffer : IG

    auto tk_yyyyyy_xxxx = pbuffer.data(idx_kin_ig + 315);

    auto tk_yyyyyy_xxxy = pbuffer.data(idx_kin_ig + 316);

    auto tk_yyyyyy_xxxz = pbuffer.data(idx_kin_ig + 317);

    auto tk_yyyyyy_xxyy = pbuffer.data(idx_kin_ig + 318);

    auto tk_yyyyyy_xxyz = pbuffer.data(idx_kin_ig + 319);

    auto tk_yyyyyy_xxzz = pbuffer.data(idx_kin_ig + 320);

    auto tk_yyyyyy_xyyy = pbuffer.data(idx_kin_ig + 321);

    auto tk_yyyyyy_xyyz = pbuffer.data(idx_kin_ig + 322);

    auto tk_yyyyyy_xyzz = pbuffer.data(idx_kin_ig + 323);

    auto tk_yyyyyy_xzzz = pbuffer.data(idx_kin_ig + 324);

    auto tk_yyyyyy_yyyy = pbuffer.data(idx_kin_ig + 325);

    auto tk_yyyyyy_yyyz = pbuffer.data(idx_kin_ig + 326);

    auto tk_yyyyyy_yyzz = pbuffer.data(idx_kin_ig + 327);

    auto tk_yyyyyy_yzzz = pbuffer.data(idx_kin_ig + 328);

    auto tk_yyyyyy_zzzz = pbuffer.data(idx_kin_ig + 329);

#pragma omp simd aligned(pa_y,               \
                             tk_yyyy_xxxx,   \
                             tk_yyyy_xxxy,   \
                             tk_yyyy_xxxz,   \
                             tk_yyyy_xxyy,   \
                             tk_yyyy_xxyz,   \
                             tk_yyyy_xxzz,   \
                             tk_yyyy_xyyy,   \
                             tk_yyyy_xyyz,   \
                             tk_yyyy_xyzz,   \
                             tk_yyyy_xzzz,   \
                             tk_yyyy_yyyy,   \
                             tk_yyyy_yyyz,   \
                             tk_yyyy_yyzz,   \
                             tk_yyyy_yzzz,   \
                             tk_yyyy_zzzz,   \
                             tk_yyyyy_xxx,   \
                             tk_yyyyy_xxxx,  \
                             tk_yyyyy_xxxy,  \
                             tk_yyyyy_xxxz,  \
                             tk_yyyyy_xxy,   \
                             tk_yyyyy_xxyy,  \
                             tk_yyyyy_xxyz,  \
                             tk_yyyyy_xxz,   \
                             tk_yyyyy_xxzz,  \
                             tk_yyyyy_xyy,   \
                             tk_yyyyy_xyyy,  \
                             tk_yyyyy_xyyz,  \
                             tk_yyyyy_xyz,   \
                             tk_yyyyy_xyzz,  \
                             tk_yyyyy_xzz,   \
                             tk_yyyyy_xzzz,  \
                             tk_yyyyy_yyy,   \
                             tk_yyyyy_yyyy,  \
                             tk_yyyyy_yyyz,  \
                             tk_yyyyy_yyz,   \
                             tk_yyyyy_yyzz,  \
                             tk_yyyyy_yzz,   \
                             tk_yyyyy_yzzz,  \
                             tk_yyyyy_zzz,   \
                             tk_yyyyy_zzzz,  \
                             tk_yyyyyy_xxxx, \
                             tk_yyyyyy_xxxy, \
                             tk_yyyyyy_xxxz, \
                             tk_yyyyyy_xxyy, \
                             tk_yyyyyy_xxyz, \
                             tk_yyyyyy_xxzz, \
                             tk_yyyyyy_xyyy, \
                             tk_yyyyyy_xyyz, \
                             tk_yyyyyy_xyzz, \
                             tk_yyyyyy_xzzz, \
                             tk_yyyyyy_yyyy, \
                             tk_yyyyyy_yyyz, \
                             tk_yyyyyy_yyzz, \
                             tk_yyyyyy_yzzz, \
                             tk_yyyyyy_zzzz, \
                             ts_yyyy_xxxx,   \
                             ts_yyyy_xxxy,   \
                             ts_yyyy_xxxz,   \
                             ts_yyyy_xxyy,   \
                             ts_yyyy_xxyz,   \
                             ts_yyyy_xxzz,   \
                             ts_yyyy_xyyy,   \
                             ts_yyyy_xyyz,   \
                             ts_yyyy_xyzz,   \
                             ts_yyyy_xzzz,   \
                             ts_yyyy_yyyy,   \
                             ts_yyyy_yyyz,   \
                             ts_yyyy_yyzz,   \
                             ts_yyyy_yzzz,   \
                             ts_yyyy_zzzz,   \
                             ts_yyyyyy_xxxx, \
                             ts_yyyyyy_xxxy, \
                             ts_yyyyyy_xxxz, \
                             ts_yyyyyy_xxyy, \
                             ts_yyyyyy_xxyz, \
                             ts_yyyyyy_xxzz, \
                             ts_yyyyyy_xyyy, \
                             ts_yyyyyy_xyyz, \
                             ts_yyyyyy_xyzz, \
                             ts_yyyyyy_xzzz, \
                             ts_yyyyyy_yyyy, \
                             ts_yyyyyy_yyyz, \
                             ts_yyyyyy_yyzz, \
                             ts_yyyyyy_yzzz, \
                             ts_yyyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyyy_xxxx[i] =
            -10.0 * ts_yyyy_xxxx[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxx[i] * fe_0 + tk_yyyyy_xxxx[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxx[i] * fz_0;

        tk_yyyyyy_xxxy[i] = -10.0 * ts_yyyy_xxxy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxy[i] * fe_0 + tk_yyyyy_xxx[i] * fe_0 +
                            tk_yyyyy_xxxy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxy[i] * fz_0;

        tk_yyyyyy_xxxz[i] =
            -10.0 * ts_yyyy_xxxz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxz[i] * fe_0 + tk_yyyyy_xxxz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxz[i] * fz_0;

        tk_yyyyyy_xxyy[i] = -10.0 * ts_yyyy_xxyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyy[i] * fe_0 + 2.0 * tk_yyyyy_xxy[i] * fe_0 +
                            tk_yyyyy_xxyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyy[i] * fz_0;

        tk_yyyyyy_xxyz[i] = -10.0 * ts_yyyy_xxyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyz[i] * fe_0 + tk_yyyyy_xxz[i] * fe_0 +
                            tk_yyyyy_xxyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyz[i] * fz_0;

        tk_yyyyyy_xxzz[i] =
            -10.0 * ts_yyyy_xxzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxzz[i] * fe_0 + tk_yyyyy_xxzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxzz[i] * fz_0;

        tk_yyyyyy_xyyy[i] = -10.0 * ts_yyyy_xyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyy[i] * fe_0 + 3.0 * tk_yyyyy_xyy[i] * fe_0 +
                            tk_yyyyy_xyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyy[i] * fz_0;

        tk_yyyyyy_xyyz[i] = -10.0 * ts_yyyy_xyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyz[i] * fe_0 + 2.0 * tk_yyyyy_xyz[i] * fe_0 +
                            tk_yyyyy_xyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyz[i] * fz_0;

        tk_yyyyyy_xyzz[i] = -10.0 * ts_yyyy_xyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyzz[i] * fe_0 + tk_yyyyy_xzz[i] * fe_0 +
                            tk_yyyyy_xyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyzz[i] * fz_0;

        tk_yyyyyy_xzzz[i] =
            -10.0 * ts_yyyy_xzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xzzz[i] * fe_0 + tk_yyyyy_xzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xzzz[i] * fz_0;

        tk_yyyyyy_yyyy[i] = -10.0 * ts_yyyy_yyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyy[i] * fe_0 + 4.0 * tk_yyyyy_yyy[i] * fe_0 +
                            tk_yyyyy_yyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyy[i] * fz_0;

        tk_yyyyyy_yyyz[i] = -10.0 * ts_yyyy_yyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyz[i] * fe_0 + 3.0 * tk_yyyyy_yyz[i] * fe_0 +
                            tk_yyyyy_yyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyz[i] * fz_0;

        tk_yyyyyy_yyzz[i] = -10.0 * ts_yyyy_yyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyzz[i] * fe_0 + 2.0 * tk_yyyyy_yzz[i] * fe_0 +
                            tk_yyyyy_yyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyzz[i] * fz_0;

        tk_yyyyyy_yzzz[i] = -10.0 * ts_yyyy_yzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yzzz[i] * fe_0 + tk_yyyyy_zzz[i] * fe_0 +
                            tk_yyyyy_yzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yzzz[i] * fz_0;

        tk_yyyyyy_zzzz[i] =
            -10.0 * ts_yyyy_zzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_zzzz[i] * fe_0 + tk_yyyyy_zzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_zzzz[i] * fz_0;
    }

    // Set up 330-345 components of targeted buffer : IG

    auto tk_yyyyyz_xxxx = pbuffer.data(idx_kin_ig + 330);

    auto tk_yyyyyz_xxxy = pbuffer.data(idx_kin_ig + 331);

    auto tk_yyyyyz_xxxz = pbuffer.data(idx_kin_ig + 332);

    auto tk_yyyyyz_xxyy = pbuffer.data(idx_kin_ig + 333);

    auto tk_yyyyyz_xxyz = pbuffer.data(idx_kin_ig + 334);

    auto tk_yyyyyz_xxzz = pbuffer.data(idx_kin_ig + 335);

    auto tk_yyyyyz_xyyy = pbuffer.data(idx_kin_ig + 336);

    auto tk_yyyyyz_xyyz = pbuffer.data(idx_kin_ig + 337);

    auto tk_yyyyyz_xyzz = pbuffer.data(idx_kin_ig + 338);

    auto tk_yyyyyz_xzzz = pbuffer.data(idx_kin_ig + 339);

    auto tk_yyyyyz_yyyy = pbuffer.data(idx_kin_ig + 340);

    auto tk_yyyyyz_yyyz = pbuffer.data(idx_kin_ig + 341);

    auto tk_yyyyyz_yyzz = pbuffer.data(idx_kin_ig + 342);

    auto tk_yyyyyz_yzzz = pbuffer.data(idx_kin_ig + 343);

    auto tk_yyyyyz_zzzz = pbuffer.data(idx_kin_ig + 344);

#pragma omp simd aligned(pa_z,               \
                             tk_yyyyy_xxx,   \
                             tk_yyyyy_xxxx,  \
                             tk_yyyyy_xxxy,  \
                             tk_yyyyy_xxxz,  \
                             tk_yyyyy_xxy,   \
                             tk_yyyyy_xxyy,  \
                             tk_yyyyy_xxyz,  \
                             tk_yyyyy_xxz,   \
                             tk_yyyyy_xxzz,  \
                             tk_yyyyy_xyy,   \
                             tk_yyyyy_xyyy,  \
                             tk_yyyyy_xyyz,  \
                             tk_yyyyy_xyz,   \
                             tk_yyyyy_xyzz,  \
                             tk_yyyyy_xzz,   \
                             tk_yyyyy_xzzz,  \
                             tk_yyyyy_yyy,   \
                             tk_yyyyy_yyyy,  \
                             tk_yyyyy_yyyz,  \
                             tk_yyyyy_yyz,   \
                             tk_yyyyy_yyzz,  \
                             tk_yyyyy_yzz,   \
                             tk_yyyyy_yzzz,  \
                             tk_yyyyy_zzz,   \
                             tk_yyyyy_zzzz,  \
                             tk_yyyyyz_xxxx, \
                             tk_yyyyyz_xxxy, \
                             tk_yyyyyz_xxxz, \
                             tk_yyyyyz_xxyy, \
                             tk_yyyyyz_xxyz, \
                             tk_yyyyyz_xxzz, \
                             tk_yyyyyz_xyyy, \
                             tk_yyyyyz_xyyz, \
                             tk_yyyyyz_xyzz, \
                             tk_yyyyyz_xzzz, \
                             tk_yyyyyz_yyyy, \
                             tk_yyyyyz_yyyz, \
                             tk_yyyyyz_yyzz, \
                             tk_yyyyyz_yzzz, \
                             tk_yyyyyz_zzzz, \
                             ts_yyyyyz_xxxx, \
                             ts_yyyyyz_xxxy, \
                             ts_yyyyyz_xxxz, \
                             ts_yyyyyz_xxyy, \
                             ts_yyyyyz_xxyz, \
                             ts_yyyyyz_xxzz, \
                             ts_yyyyyz_xyyy, \
                             ts_yyyyyz_xyyz, \
                             ts_yyyyyz_xyzz, \
                             ts_yyyyyz_xzzz, \
                             ts_yyyyyz_yyyy, \
                             ts_yyyyyz_yyyz, \
                             ts_yyyyyz_yyzz, \
                             ts_yyyyyz_yzzz, \
                             ts_yyyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyyz_xxxx[i] = tk_yyyyy_xxxx[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxx[i] * fz_0;

        tk_yyyyyz_xxxy[i] = tk_yyyyy_xxxy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxy[i] * fz_0;

        tk_yyyyyz_xxxz[i] = tk_yyyyy_xxx[i] * fe_0 + tk_yyyyy_xxxz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxz[i] * fz_0;

        tk_yyyyyz_xxyy[i] = tk_yyyyy_xxyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyy[i] * fz_0;

        tk_yyyyyz_xxyz[i] = tk_yyyyy_xxy[i] * fe_0 + tk_yyyyy_xxyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyz[i] * fz_0;

        tk_yyyyyz_xxzz[i] = 2.0 * tk_yyyyy_xxz[i] * fe_0 + tk_yyyyy_xxzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxzz[i] * fz_0;

        tk_yyyyyz_xyyy[i] = tk_yyyyy_xyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyy[i] * fz_0;

        tk_yyyyyz_xyyz[i] = tk_yyyyy_xyy[i] * fe_0 + tk_yyyyy_xyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyz[i] * fz_0;

        tk_yyyyyz_xyzz[i] = 2.0 * tk_yyyyy_xyz[i] * fe_0 + tk_yyyyy_xyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyzz[i] * fz_0;

        tk_yyyyyz_xzzz[i] = 3.0 * tk_yyyyy_xzz[i] * fe_0 + tk_yyyyy_xzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xzzz[i] * fz_0;

        tk_yyyyyz_yyyy[i] = tk_yyyyy_yyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyy[i] * fz_0;

        tk_yyyyyz_yyyz[i] = tk_yyyyy_yyy[i] * fe_0 + tk_yyyyy_yyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyz[i] * fz_0;

        tk_yyyyyz_yyzz[i] = 2.0 * tk_yyyyy_yyz[i] * fe_0 + tk_yyyyy_yyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyzz[i] * fz_0;

        tk_yyyyyz_yzzz[i] = 3.0 * tk_yyyyy_yzz[i] * fe_0 + tk_yyyyy_yzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yzzz[i] * fz_0;

        tk_yyyyyz_zzzz[i] = 4.0 * tk_yyyyy_zzz[i] * fe_0 + tk_yyyyy_zzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_zzzz[i] * fz_0;
    }

    // Set up 345-360 components of targeted buffer : IG

    auto tk_yyyyzz_xxxx = pbuffer.data(idx_kin_ig + 345);

    auto tk_yyyyzz_xxxy = pbuffer.data(idx_kin_ig + 346);

    auto tk_yyyyzz_xxxz = pbuffer.data(idx_kin_ig + 347);

    auto tk_yyyyzz_xxyy = pbuffer.data(idx_kin_ig + 348);

    auto tk_yyyyzz_xxyz = pbuffer.data(idx_kin_ig + 349);

    auto tk_yyyyzz_xxzz = pbuffer.data(idx_kin_ig + 350);

    auto tk_yyyyzz_xyyy = pbuffer.data(idx_kin_ig + 351);

    auto tk_yyyyzz_xyyz = pbuffer.data(idx_kin_ig + 352);

    auto tk_yyyyzz_xyzz = pbuffer.data(idx_kin_ig + 353);

    auto tk_yyyyzz_xzzz = pbuffer.data(idx_kin_ig + 354);

    auto tk_yyyyzz_yyyy = pbuffer.data(idx_kin_ig + 355);

    auto tk_yyyyzz_yyyz = pbuffer.data(idx_kin_ig + 356);

    auto tk_yyyyzz_yyzz = pbuffer.data(idx_kin_ig + 357);

    auto tk_yyyyzz_yzzz = pbuffer.data(idx_kin_ig + 358);

    auto tk_yyyyzz_zzzz = pbuffer.data(idx_kin_ig + 359);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_yyyy_xxxy,   \
                             tk_yyyy_xxyy,   \
                             tk_yyyy_xyyy,   \
                             tk_yyyy_yyyy,   \
                             tk_yyyyz_xxxy,  \
                             tk_yyyyz_xxyy,  \
                             tk_yyyyz_xyyy,  \
                             tk_yyyyz_yyyy,  \
                             tk_yyyyzz_xxxx, \
                             tk_yyyyzz_xxxy, \
                             tk_yyyyzz_xxxz, \
                             tk_yyyyzz_xxyy, \
                             tk_yyyyzz_xxyz, \
                             tk_yyyyzz_xxzz, \
                             tk_yyyyzz_xyyy, \
                             tk_yyyyzz_xyyz, \
                             tk_yyyyzz_xyzz, \
                             tk_yyyyzz_xzzz, \
                             tk_yyyyzz_yyyy, \
                             tk_yyyyzz_yyyz, \
                             tk_yyyyzz_yyzz, \
                             tk_yyyyzz_yzzz, \
                             tk_yyyyzz_zzzz, \
                             tk_yyyzz_xxxx,  \
                             tk_yyyzz_xxxz,  \
                             tk_yyyzz_xxyz,  \
                             tk_yyyzz_xxz,   \
                             tk_yyyzz_xxzz,  \
                             tk_yyyzz_xyyz,  \
                             tk_yyyzz_xyz,   \
                             tk_yyyzz_xyzz,  \
                             tk_yyyzz_xzz,   \
                             tk_yyyzz_xzzz,  \
                             tk_yyyzz_yyyz,  \
                             tk_yyyzz_yyz,   \
                             tk_yyyzz_yyzz,  \
                             tk_yyyzz_yzz,   \
                             tk_yyyzz_yzzz,  \
                             tk_yyyzz_zzz,   \
                             tk_yyyzz_zzzz,  \
                             tk_yyzz_xxxx,   \
                             tk_yyzz_xxxz,   \
                             tk_yyzz_xxyz,   \
                             tk_yyzz_xxzz,   \
                             tk_yyzz_xyyz,   \
                             tk_yyzz_xyzz,   \
                             tk_yyzz_xzzz,   \
                             tk_yyzz_yyyz,   \
                             tk_yyzz_yyzz,   \
                             tk_yyzz_yzzz,   \
                             tk_yyzz_zzzz,   \
                             ts_yyyy_xxxy,   \
                             ts_yyyy_xxyy,   \
                             ts_yyyy_xyyy,   \
                             ts_yyyy_yyyy,   \
                             ts_yyyyzz_xxxx, \
                             ts_yyyyzz_xxxy, \
                             ts_yyyyzz_xxxz, \
                             ts_yyyyzz_xxyy, \
                             ts_yyyyzz_xxyz, \
                             ts_yyyyzz_xxzz, \
                             ts_yyyyzz_xyyy, \
                             ts_yyyyzz_xyyz, \
                             ts_yyyyzz_xyzz, \
                             ts_yyyyzz_xzzz, \
                             ts_yyyyzz_yyyy, \
                             ts_yyyyzz_yyyz, \
                             ts_yyyyzz_yyzz, \
                             ts_yyyyzz_yzzz, \
                             ts_yyyyzz_zzzz, \
                             ts_yyzz_xxxx,   \
                             ts_yyzz_xxxz,   \
                             ts_yyzz_xxyz,   \
                             ts_yyzz_xxzz,   \
                             ts_yyzz_xyyz,   \
                             ts_yyzz_xyzz,   \
                             ts_yyzz_xzzz,   \
                             ts_yyzz_yyyz,   \
                             ts_yyzz_yyzz,   \
                             ts_yyzz_yzzz,   \
                             ts_yyzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyzz_xxxx[i] =
            -6.0 * ts_yyzz_xxxx[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxx[i] * fe_0 + tk_yyyzz_xxxx[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxx[i] * fz_0;

        tk_yyyyzz_xxxy[i] =
            -2.0 * ts_yyyy_xxxy[i] * fbe_0 * fz_0 + tk_yyyy_xxxy[i] * fe_0 + tk_yyyyz_xxxy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxxy[i] * fz_0;

        tk_yyyyzz_xxxz[i] =
            -6.0 * ts_yyzz_xxxz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxz[i] * fe_0 + tk_yyyzz_xxxz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxz[i] * fz_0;

        tk_yyyyzz_xxyy[i] =
            -2.0 * ts_yyyy_xxyy[i] * fbe_0 * fz_0 + tk_yyyy_xxyy[i] * fe_0 + tk_yyyyz_xxyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxyy[i] * fz_0;

        tk_yyyyzz_xxyz[i] = -6.0 * ts_yyzz_xxyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyz[i] * fe_0 + tk_yyyzz_xxz[i] * fe_0 +
                            tk_yyyzz_xxyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxyz[i] * fz_0;

        tk_yyyyzz_xxzz[i] =
            -6.0 * ts_yyzz_xxzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxzz[i] * fe_0 + tk_yyyzz_xxzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxzz[i] * fz_0;

        tk_yyyyzz_xyyy[i] =
            -2.0 * ts_yyyy_xyyy[i] * fbe_0 * fz_0 + tk_yyyy_xyyy[i] * fe_0 + tk_yyyyz_xyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xyyy[i] * fz_0;

        tk_yyyyzz_xyyz[i] = -6.0 * ts_yyzz_xyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyz[i] * fe_0 + 2.0 * tk_yyyzz_xyz[i] * fe_0 +
                            tk_yyyzz_xyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyyz[i] * fz_0;

        tk_yyyyzz_xyzz[i] = -6.0 * ts_yyzz_xyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyzz[i] * fe_0 + tk_yyyzz_xzz[i] * fe_0 +
                            tk_yyyzz_xyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyzz[i] * fz_0;

        tk_yyyyzz_xzzz[i] =
            -6.0 * ts_yyzz_xzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xzzz[i] * fe_0 + tk_yyyzz_xzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xzzz[i] * fz_0;

        tk_yyyyzz_yyyy[i] =
            -2.0 * ts_yyyy_yyyy[i] * fbe_0 * fz_0 + tk_yyyy_yyyy[i] * fe_0 + tk_yyyyz_yyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_yyyy[i] * fz_0;

        tk_yyyyzz_yyyz[i] = -6.0 * ts_yyzz_yyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyz[i] * fe_0 + 3.0 * tk_yyyzz_yyz[i] * fe_0 +
                            tk_yyyzz_yyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyyz[i] * fz_0;

        tk_yyyyzz_yyzz[i] = -6.0 * ts_yyzz_yyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyzz[i] * fe_0 + 2.0 * tk_yyyzz_yzz[i] * fe_0 +
                            tk_yyyzz_yyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyzz[i] * fz_0;

        tk_yyyyzz_yzzz[i] = -6.0 * ts_yyzz_yzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yzzz[i] * fe_0 + tk_yyyzz_zzz[i] * fe_0 +
                            tk_yyyzz_yzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yzzz[i] * fz_0;

        tk_yyyyzz_zzzz[i] =
            -6.0 * ts_yyzz_zzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_zzzz[i] * fe_0 + tk_yyyzz_zzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_zzzz[i] * fz_0;
    }

    // Set up 360-375 components of targeted buffer : IG

    auto tk_yyyzzz_xxxx = pbuffer.data(idx_kin_ig + 360);

    auto tk_yyyzzz_xxxy = pbuffer.data(idx_kin_ig + 361);

    auto tk_yyyzzz_xxxz = pbuffer.data(idx_kin_ig + 362);

    auto tk_yyyzzz_xxyy = pbuffer.data(idx_kin_ig + 363);

    auto tk_yyyzzz_xxyz = pbuffer.data(idx_kin_ig + 364);

    auto tk_yyyzzz_xxzz = pbuffer.data(idx_kin_ig + 365);

    auto tk_yyyzzz_xyyy = pbuffer.data(idx_kin_ig + 366);

    auto tk_yyyzzz_xyyz = pbuffer.data(idx_kin_ig + 367);

    auto tk_yyyzzz_xyzz = pbuffer.data(idx_kin_ig + 368);

    auto tk_yyyzzz_xzzz = pbuffer.data(idx_kin_ig + 369);

    auto tk_yyyzzz_yyyy = pbuffer.data(idx_kin_ig + 370);

    auto tk_yyyzzz_yyyz = pbuffer.data(idx_kin_ig + 371);

    auto tk_yyyzzz_yyzz = pbuffer.data(idx_kin_ig + 372);

    auto tk_yyyzzz_yzzz = pbuffer.data(idx_kin_ig + 373);

    auto tk_yyyzzz_zzzz = pbuffer.data(idx_kin_ig + 374);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_yyyz_xxxy,   \
                             tk_yyyz_xxyy,   \
                             tk_yyyz_xyyy,   \
                             tk_yyyz_yyyy,   \
                             tk_yyyzz_xxxy,  \
                             tk_yyyzz_xxyy,  \
                             tk_yyyzz_xyyy,  \
                             tk_yyyzz_yyyy,  \
                             tk_yyyzzz_xxxx, \
                             tk_yyyzzz_xxxy, \
                             tk_yyyzzz_xxxz, \
                             tk_yyyzzz_xxyy, \
                             tk_yyyzzz_xxyz, \
                             tk_yyyzzz_xxzz, \
                             tk_yyyzzz_xyyy, \
                             tk_yyyzzz_xyyz, \
                             tk_yyyzzz_xyzz, \
                             tk_yyyzzz_xzzz, \
                             tk_yyyzzz_yyyy, \
                             tk_yyyzzz_yyyz, \
                             tk_yyyzzz_yyzz, \
                             tk_yyyzzz_yzzz, \
                             tk_yyyzzz_zzzz, \
                             tk_yyzzz_xxxx,  \
                             tk_yyzzz_xxxz,  \
                             tk_yyzzz_xxyz,  \
                             tk_yyzzz_xxz,   \
                             tk_yyzzz_xxzz,  \
                             tk_yyzzz_xyyz,  \
                             tk_yyzzz_xyz,   \
                             tk_yyzzz_xyzz,  \
                             tk_yyzzz_xzz,   \
                             tk_yyzzz_xzzz,  \
                             tk_yyzzz_yyyz,  \
                             tk_yyzzz_yyz,   \
                             tk_yyzzz_yyzz,  \
                             tk_yyzzz_yzz,   \
                             tk_yyzzz_yzzz,  \
                             tk_yyzzz_zzz,   \
                             tk_yyzzz_zzzz,  \
                             tk_yzzz_xxxx,   \
                             tk_yzzz_xxxz,   \
                             tk_yzzz_xxyz,   \
                             tk_yzzz_xxzz,   \
                             tk_yzzz_xyyz,   \
                             tk_yzzz_xyzz,   \
                             tk_yzzz_xzzz,   \
                             tk_yzzz_yyyz,   \
                             tk_yzzz_yyzz,   \
                             tk_yzzz_yzzz,   \
                             tk_yzzz_zzzz,   \
                             ts_yyyz_xxxy,   \
                             ts_yyyz_xxyy,   \
                             ts_yyyz_xyyy,   \
                             ts_yyyz_yyyy,   \
                             ts_yyyzzz_xxxx, \
                             ts_yyyzzz_xxxy, \
                             ts_yyyzzz_xxxz, \
                             ts_yyyzzz_xxyy, \
                             ts_yyyzzz_xxyz, \
                             ts_yyyzzz_xxzz, \
                             ts_yyyzzz_xyyy, \
                             ts_yyyzzz_xyyz, \
                             ts_yyyzzz_xyzz, \
                             ts_yyyzzz_xzzz, \
                             ts_yyyzzz_yyyy, \
                             ts_yyyzzz_yyyz, \
                             ts_yyyzzz_yyzz, \
                             ts_yyyzzz_yzzz, \
                             ts_yyyzzz_zzzz, \
                             ts_yzzz_xxxx,   \
                             ts_yzzz_xxxz,   \
                             ts_yzzz_xxyz,   \
                             ts_yzzz_xxzz,   \
                             ts_yzzz_xyyz,   \
                             ts_yzzz_xyzz,   \
                             ts_yzzz_xzzz,   \
                             ts_yzzz_yyyz,   \
                             ts_yzzz_yyzz,   \
                             ts_yzzz_yzzz,   \
                             ts_yzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzzz_xxxx[i] =
            -4.0 * ts_yzzz_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxx[i] * fe_0 + tk_yyzzz_xxxx[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxx[i] * fz_0;

        tk_yyyzzz_xxxy[i] =
            -4.0 * ts_yyyz_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxxy[i] * fe_0 + tk_yyyzz_xxxy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xxxy[i] * fz_0;

        tk_yyyzzz_xxxz[i] =
            -4.0 * ts_yzzz_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxz[i] * fe_0 + tk_yyzzz_xxxz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxz[i] * fz_0;

        tk_yyyzzz_xxyy[i] =
            -4.0 * ts_yyyz_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxyy[i] * fe_0 + tk_yyyzz_xxyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xxyy[i] * fz_0;

        tk_yyyzzz_xxyz[i] = -4.0 * ts_yzzz_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxyz[i] * fe_0 + tk_yyzzz_xxz[i] * fe_0 +
                            tk_yyzzz_xxyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxyz[i] * fz_0;

        tk_yyyzzz_xxzz[i] =
            -4.0 * ts_yzzz_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxzz[i] * fe_0 + tk_yyzzz_xxzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxzz[i] * fz_0;

        tk_yyyzzz_xyyy[i] =
            -4.0 * ts_yyyz_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xyyy[i] * fe_0 + tk_yyyzz_xyyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xyyy[i] * fz_0;

        tk_yyyzzz_xyyz[i] = -4.0 * ts_yzzz_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyyz[i] * fe_0 + 2.0 * tk_yyzzz_xyz[i] * fe_0 +
                            tk_yyzzz_xyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyyz[i] * fz_0;

        tk_yyyzzz_xyzz[i] = -4.0 * ts_yzzz_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyzz[i] * fe_0 + tk_yyzzz_xzz[i] * fe_0 +
                            tk_yyzzz_xyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyzz[i] * fz_0;

        tk_yyyzzz_xzzz[i] =
            -4.0 * ts_yzzz_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xzzz[i] * fe_0 + tk_yyzzz_xzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xzzz[i] * fz_0;

        tk_yyyzzz_yyyy[i] =
            -4.0 * ts_yyyz_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_yyyy[i] * fe_0 + tk_yyyzz_yyyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_yyyy[i] * fz_0;

        tk_yyyzzz_yyyz[i] = -4.0 * ts_yzzz_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyyz[i] * fe_0 + 3.0 * tk_yyzzz_yyz[i] * fe_0 +
                            tk_yyzzz_yyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyyz[i] * fz_0;

        tk_yyyzzz_yyzz[i] = -4.0 * ts_yzzz_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyzz[i] * fe_0 + 2.0 * tk_yyzzz_yzz[i] * fe_0 +
                            tk_yyzzz_yyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyzz[i] * fz_0;

        tk_yyyzzz_yzzz[i] = -4.0 * ts_yzzz_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yzzz[i] * fe_0 + tk_yyzzz_zzz[i] * fe_0 +
                            tk_yyzzz_yzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yzzz[i] * fz_0;

        tk_yyyzzz_zzzz[i] =
            -4.0 * ts_yzzz_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_zzzz[i] * fe_0 + tk_yyzzz_zzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_zzzz[i] * fz_0;
    }

    // Set up 375-390 components of targeted buffer : IG

    auto tk_yyzzzz_xxxx = pbuffer.data(idx_kin_ig + 375);

    auto tk_yyzzzz_xxxy = pbuffer.data(idx_kin_ig + 376);

    auto tk_yyzzzz_xxxz = pbuffer.data(idx_kin_ig + 377);

    auto tk_yyzzzz_xxyy = pbuffer.data(idx_kin_ig + 378);

    auto tk_yyzzzz_xxyz = pbuffer.data(idx_kin_ig + 379);

    auto tk_yyzzzz_xxzz = pbuffer.data(idx_kin_ig + 380);

    auto tk_yyzzzz_xyyy = pbuffer.data(idx_kin_ig + 381);

    auto tk_yyzzzz_xyyz = pbuffer.data(idx_kin_ig + 382);

    auto tk_yyzzzz_xyzz = pbuffer.data(idx_kin_ig + 383);

    auto tk_yyzzzz_xzzz = pbuffer.data(idx_kin_ig + 384);

    auto tk_yyzzzz_yyyy = pbuffer.data(idx_kin_ig + 385);

    auto tk_yyzzzz_yyyz = pbuffer.data(idx_kin_ig + 386);

    auto tk_yyzzzz_yyzz = pbuffer.data(idx_kin_ig + 387);

    auto tk_yyzzzz_yzzz = pbuffer.data(idx_kin_ig + 388);

    auto tk_yyzzzz_zzzz = pbuffer.data(idx_kin_ig + 389);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_yyzz_xxxy,   \
                             tk_yyzz_xxyy,   \
                             tk_yyzz_xyyy,   \
                             tk_yyzz_yyyy,   \
                             tk_yyzzz_xxxy,  \
                             tk_yyzzz_xxyy,  \
                             tk_yyzzz_xyyy,  \
                             tk_yyzzz_yyyy,  \
                             tk_yyzzzz_xxxx, \
                             tk_yyzzzz_xxxy, \
                             tk_yyzzzz_xxxz, \
                             tk_yyzzzz_xxyy, \
                             tk_yyzzzz_xxyz, \
                             tk_yyzzzz_xxzz, \
                             tk_yyzzzz_xyyy, \
                             tk_yyzzzz_xyyz, \
                             tk_yyzzzz_xyzz, \
                             tk_yyzzzz_xzzz, \
                             tk_yyzzzz_yyyy, \
                             tk_yyzzzz_yyyz, \
                             tk_yyzzzz_yyzz, \
                             tk_yyzzzz_yzzz, \
                             tk_yyzzzz_zzzz, \
                             tk_yzzzz_xxxx,  \
                             tk_yzzzz_xxxz,  \
                             tk_yzzzz_xxyz,  \
                             tk_yzzzz_xxz,   \
                             tk_yzzzz_xxzz,  \
                             tk_yzzzz_xyyz,  \
                             tk_yzzzz_xyz,   \
                             tk_yzzzz_xyzz,  \
                             tk_yzzzz_xzz,   \
                             tk_yzzzz_xzzz,  \
                             tk_yzzzz_yyyz,  \
                             tk_yzzzz_yyz,   \
                             tk_yzzzz_yyzz,  \
                             tk_yzzzz_yzz,   \
                             tk_yzzzz_yzzz,  \
                             tk_yzzzz_zzz,   \
                             tk_yzzzz_zzzz,  \
                             tk_zzzz_xxxx,   \
                             tk_zzzz_xxxz,   \
                             tk_zzzz_xxyz,   \
                             tk_zzzz_xxzz,   \
                             tk_zzzz_xyyz,   \
                             tk_zzzz_xyzz,   \
                             tk_zzzz_xzzz,   \
                             tk_zzzz_yyyz,   \
                             tk_zzzz_yyzz,   \
                             tk_zzzz_yzzz,   \
                             tk_zzzz_zzzz,   \
                             ts_yyzz_xxxy,   \
                             ts_yyzz_xxyy,   \
                             ts_yyzz_xyyy,   \
                             ts_yyzz_yyyy,   \
                             ts_yyzzzz_xxxx, \
                             ts_yyzzzz_xxxy, \
                             ts_yyzzzz_xxxz, \
                             ts_yyzzzz_xxyy, \
                             ts_yyzzzz_xxyz, \
                             ts_yyzzzz_xxzz, \
                             ts_yyzzzz_xyyy, \
                             ts_yyzzzz_xyyz, \
                             ts_yyzzzz_xyzz, \
                             ts_yyzzzz_xzzz, \
                             ts_yyzzzz_yyyy, \
                             ts_yyzzzz_yyyz, \
                             ts_yyzzzz_yyzz, \
                             ts_yyzzzz_yzzz, \
                             ts_yyzzzz_zzzz, \
                             ts_zzzz_xxxx,   \
                             ts_zzzz_xxxz,   \
                             ts_zzzz_xxyz,   \
                             ts_zzzz_xxzz,   \
                             ts_zzzz_xyyz,   \
                             ts_zzzz_xyzz,   \
                             ts_zzzz_xzzz,   \
                             ts_zzzz_yyyz,   \
                             ts_zzzz_yyzz,   \
                             ts_zzzz_yzzz,   \
                             ts_zzzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzzz_xxxx[i] =
            -2.0 * ts_zzzz_xxxx[i] * fbe_0 * fz_0 + tk_zzzz_xxxx[i] * fe_0 + tk_yzzzz_xxxx[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxx[i] * fz_0;

        tk_yyzzzz_xxxy[i] =
            -6.0 * ts_yyzz_xxxy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxy[i] * fe_0 + tk_yyzzz_xxxy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xxxy[i] * fz_0;

        tk_yyzzzz_xxxz[i] =
            -2.0 * ts_zzzz_xxxz[i] * fbe_0 * fz_0 + tk_zzzz_xxxz[i] * fe_0 + tk_yzzzz_xxxz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxz[i] * fz_0;

        tk_yyzzzz_xxyy[i] =
            -6.0 * ts_yyzz_xxyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyy[i] * fe_0 + tk_yyzzz_xxyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xxyy[i] * fz_0;

        tk_yyzzzz_xxyz[i] = -2.0 * ts_zzzz_xxyz[i] * fbe_0 * fz_0 + tk_zzzz_xxyz[i] * fe_0 + tk_yzzzz_xxz[i] * fe_0 + tk_yzzzz_xxyz[i] * pa_y[i] +
                            2.0 * ts_yyzzzz_xxyz[i] * fz_0;

        tk_yyzzzz_xxzz[i] =
            -2.0 * ts_zzzz_xxzz[i] * fbe_0 * fz_0 + tk_zzzz_xxzz[i] * fe_0 + tk_yzzzz_xxzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxzz[i] * fz_0;

        tk_yyzzzz_xyyy[i] =
            -6.0 * ts_yyzz_xyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyy[i] * fe_0 + tk_yyzzz_xyyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xyyy[i] * fz_0;

        tk_yyzzzz_xyyz[i] = -2.0 * ts_zzzz_xyyz[i] * fbe_0 * fz_0 + tk_zzzz_xyyz[i] * fe_0 + 2.0 * tk_yzzzz_xyz[i] * fe_0 +
                            tk_yzzzz_xyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyyz[i] * fz_0;

        tk_yyzzzz_xyzz[i] = -2.0 * ts_zzzz_xyzz[i] * fbe_0 * fz_0 + tk_zzzz_xyzz[i] * fe_0 + tk_yzzzz_xzz[i] * fe_0 + tk_yzzzz_xyzz[i] * pa_y[i] +
                            2.0 * ts_yyzzzz_xyzz[i] * fz_0;

        tk_yyzzzz_xzzz[i] =
            -2.0 * ts_zzzz_xzzz[i] * fbe_0 * fz_0 + tk_zzzz_xzzz[i] * fe_0 + tk_yzzzz_xzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xzzz[i] * fz_0;

        tk_yyzzzz_yyyy[i] =
            -6.0 * ts_yyzz_yyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyy[i] * fe_0 + tk_yyzzz_yyyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_yyyy[i] * fz_0;

        tk_yyzzzz_yyyz[i] = -2.0 * ts_zzzz_yyyz[i] * fbe_0 * fz_0 + tk_zzzz_yyyz[i] * fe_0 + 3.0 * tk_yzzzz_yyz[i] * fe_0 +
                            tk_yzzzz_yyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyyz[i] * fz_0;

        tk_yyzzzz_yyzz[i] = -2.0 * ts_zzzz_yyzz[i] * fbe_0 * fz_0 + tk_zzzz_yyzz[i] * fe_0 + 2.0 * tk_yzzzz_yzz[i] * fe_0 +
                            tk_yzzzz_yyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyzz[i] * fz_0;

        tk_yyzzzz_yzzz[i] = -2.0 * ts_zzzz_yzzz[i] * fbe_0 * fz_0 + tk_zzzz_yzzz[i] * fe_0 + tk_yzzzz_zzz[i] * fe_0 + tk_yzzzz_yzzz[i] * pa_y[i] +
                            2.0 * ts_yyzzzz_yzzz[i] * fz_0;

        tk_yyzzzz_zzzz[i] =
            -2.0 * ts_zzzz_zzzz[i] * fbe_0 * fz_0 + tk_zzzz_zzzz[i] * fe_0 + tk_yzzzz_zzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_zzzz[i] * fz_0;
    }

    // Set up 390-405 components of targeted buffer : IG

    auto tk_yzzzzz_xxxx = pbuffer.data(idx_kin_ig + 390);

    auto tk_yzzzzz_xxxy = pbuffer.data(idx_kin_ig + 391);

    auto tk_yzzzzz_xxxz = pbuffer.data(idx_kin_ig + 392);

    auto tk_yzzzzz_xxyy = pbuffer.data(idx_kin_ig + 393);

    auto tk_yzzzzz_xxyz = pbuffer.data(idx_kin_ig + 394);

    auto tk_yzzzzz_xxzz = pbuffer.data(idx_kin_ig + 395);

    auto tk_yzzzzz_xyyy = pbuffer.data(idx_kin_ig + 396);

    auto tk_yzzzzz_xyyz = pbuffer.data(idx_kin_ig + 397);

    auto tk_yzzzzz_xyzz = pbuffer.data(idx_kin_ig + 398);

    auto tk_yzzzzz_xzzz = pbuffer.data(idx_kin_ig + 399);

    auto tk_yzzzzz_yyyy = pbuffer.data(idx_kin_ig + 400);

    auto tk_yzzzzz_yyyz = pbuffer.data(idx_kin_ig + 401);

    auto tk_yzzzzz_yyzz = pbuffer.data(idx_kin_ig + 402);

    auto tk_yzzzzz_yzzz = pbuffer.data(idx_kin_ig + 403);

    auto tk_yzzzzz_zzzz = pbuffer.data(idx_kin_ig + 404);

#pragma omp simd aligned(pa_y,               \
                             tk_yzzzzz_xxxx, \
                             tk_yzzzzz_xxxy, \
                             tk_yzzzzz_xxxz, \
                             tk_yzzzzz_xxyy, \
                             tk_yzzzzz_xxyz, \
                             tk_yzzzzz_xxzz, \
                             tk_yzzzzz_xyyy, \
                             tk_yzzzzz_xyyz, \
                             tk_yzzzzz_xyzz, \
                             tk_yzzzzz_xzzz, \
                             tk_yzzzzz_yyyy, \
                             tk_yzzzzz_yyyz, \
                             tk_yzzzzz_yyzz, \
                             tk_yzzzzz_yzzz, \
                             tk_yzzzzz_zzzz, \
                             tk_zzzzz_xxx,   \
                             tk_zzzzz_xxxx,  \
                             tk_zzzzz_xxxy,  \
                             tk_zzzzz_xxxz,  \
                             tk_zzzzz_xxy,   \
                             tk_zzzzz_xxyy,  \
                             tk_zzzzz_xxyz,  \
                             tk_zzzzz_xxz,   \
                             tk_zzzzz_xxzz,  \
                             tk_zzzzz_xyy,   \
                             tk_zzzzz_xyyy,  \
                             tk_zzzzz_xyyz,  \
                             tk_zzzzz_xyz,   \
                             tk_zzzzz_xyzz,  \
                             tk_zzzzz_xzz,   \
                             tk_zzzzz_xzzz,  \
                             tk_zzzzz_yyy,   \
                             tk_zzzzz_yyyy,  \
                             tk_zzzzz_yyyz,  \
                             tk_zzzzz_yyz,   \
                             tk_zzzzz_yyzz,  \
                             tk_zzzzz_yzz,   \
                             tk_zzzzz_yzzz,  \
                             tk_zzzzz_zzz,   \
                             tk_zzzzz_zzzz,  \
                             ts_yzzzzz_xxxx, \
                             ts_yzzzzz_xxxy, \
                             ts_yzzzzz_xxxz, \
                             ts_yzzzzz_xxyy, \
                             ts_yzzzzz_xxyz, \
                             ts_yzzzzz_xxzz, \
                             ts_yzzzzz_xyyy, \
                             ts_yzzzzz_xyyz, \
                             ts_yzzzzz_xyzz, \
                             ts_yzzzzz_xzzz, \
                             ts_yzzzzz_yyyy, \
                             ts_yzzzzz_yyyz, \
                             ts_yzzzzz_yyzz, \
                             ts_yzzzzz_yzzz, \
                             ts_yzzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzzz_xxxx[i] = tk_zzzzz_xxxx[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxx[i] * fz_0;

        tk_yzzzzz_xxxy[i] = tk_zzzzz_xxx[i] * fe_0 + tk_zzzzz_xxxy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxy[i] * fz_0;

        tk_yzzzzz_xxxz[i] = tk_zzzzz_xxxz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxz[i] * fz_0;

        tk_yzzzzz_xxyy[i] = 2.0 * tk_zzzzz_xxy[i] * fe_0 + tk_zzzzz_xxyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyy[i] * fz_0;

        tk_yzzzzz_xxyz[i] = tk_zzzzz_xxz[i] * fe_0 + tk_zzzzz_xxyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyz[i] * fz_0;

        tk_yzzzzz_xxzz[i] = tk_zzzzz_xxzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxzz[i] * fz_0;

        tk_yzzzzz_xyyy[i] = 3.0 * tk_zzzzz_xyy[i] * fe_0 + tk_zzzzz_xyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyy[i] * fz_0;

        tk_yzzzzz_xyyz[i] = 2.0 * tk_zzzzz_xyz[i] * fe_0 + tk_zzzzz_xyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyz[i] * fz_0;

        tk_yzzzzz_xyzz[i] = tk_zzzzz_xzz[i] * fe_0 + tk_zzzzz_xyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyzz[i] * fz_0;

        tk_yzzzzz_xzzz[i] = tk_zzzzz_xzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xzzz[i] * fz_0;

        tk_yzzzzz_yyyy[i] = 4.0 * tk_zzzzz_yyy[i] * fe_0 + tk_zzzzz_yyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyy[i] * fz_0;

        tk_yzzzzz_yyyz[i] = 3.0 * tk_zzzzz_yyz[i] * fe_0 + tk_zzzzz_yyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyz[i] * fz_0;

        tk_yzzzzz_yyzz[i] = 2.0 * tk_zzzzz_yzz[i] * fe_0 + tk_zzzzz_yyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyzz[i] * fz_0;

        tk_yzzzzz_yzzz[i] = tk_zzzzz_zzz[i] * fe_0 + tk_zzzzz_yzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yzzz[i] * fz_0;

        tk_yzzzzz_zzzz[i] = tk_zzzzz_zzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_zzzz[i] * fz_0;
    }

    // Set up 405-420 components of targeted buffer : IG

    auto tk_zzzzzz_xxxx = pbuffer.data(idx_kin_ig + 405);

    auto tk_zzzzzz_xxxy = pbuffer.data(idx_kin_ig + 406);

    auto tk_zzzzzz_xxxz = pbuffer.data(idx_kin_ig + 407);

    auto tk_zzzzzz_xxyy = pbuffer.data(idx_kin_ig + 408);

    auto tk_zzzzzz_xxyz = pbuffer.data(idx_kin_ig + 409);

    auto tk_zzzzzz_xxzz = pbuffer.data(idx_kin_ig + 410);

    auto tk_zzzzzz_xyyy = pbuffer.data(idx_kin_ig + 411);

    auto tk_zzzzzz_xyyz = pbuffer.data(idx_kin_ig + 412);

    auto tk_zzzzzz_xyzz = pbuffer.data(idx_kin_ig + 413);

    auto tk_zzzzzz_xzzz = pbuffer.data(idx_kin_ig + 414);

    auto tk_zzzzzz_yyyy = pbuffer.data(idx_kin_ig + 415);

    auto tk_zzzzzz_yyyz = pbuffer.data(idx_kin_ig + 416);

    auto tk_zzzzzz_yyzz = pbuffer.data(idx_kin_ig + 417);

    auto tk_zzzzzz_yzzz = pbuffer.data(idx_kin_ig + 418);

    auto tk_zzzzzz_zzzz = pbuffer.data(idx_kin_ig + 419);

#pragma omp simd aligned(pa_z,               \
                             tk_zzzz_xxxx,   \
                             tk_zzzz_xxxy,   \
                             tk_zzzz_xxxz,   \
                             tk_zzzz_xxyy,   \
                             tk_zzzz_xxyz,   \
                             tk_zzzz_xxzz,   \
                             tk_zzzz_xyyy,   \
                             tk_zzzz_xyyz,   \
                             tk_zzzz_xyzz,   \
                             tk_zzzz_xzzz,   \
                             tk_zzzz_yyyy,   \
                             tk_zzzz_yyyz,   \
                             tk_zzzz_yyzz,   \
                             tk_zzzz_yzzz,   \
                             tk_zzzz_zzzz,   \
                             tk_zzzzz_xxx,   \
                             tk_zzzzz_xxxx,  \
                             tk_zzzzz_xxxy,  \
                             tk_zzzzz_xxxz,  \
                             tk_zzzzz_xxy,   \
                             tk_zzzzz_xxyy,  \
                             tk_zzzzz_xxyz,  \
                             tk_zzzzz_xxz,   \
                             tk_zzzzz_xxzz,  \
                             tk_zzzzz_xyy,   \
                             tk_zzzzz_xyyy,  \
                             tk_zzzzz_xyyz,  \
                             tk_zzzzz_xyz,   \
                             tk_zzzzz_xyzz,  \
                             tk_zzzzz_xzz,   \
                             tk_zzzzz_xzzz,  \
                             tk_zzzzz_yyy,   \
                             tk_zzzzz_yyyy,  \
                             tk_zzzzz_yyyz,  \
                             tk_zzzzz_yyz,   \
                             tk_zzzzz_yyzz,  \
                             tk_zzzzz_yzz,   \
                             tk_zzzzz_yzzz,  \
                             tk_zzzzz_zzz,   \
                             tk_zzzzz_zzzz,  \
                             tk_zzzzzz_xxxx, \
                             tk_zzzzzz_xxxy, \
                             tk_zzzzzz_xxxz, \
                             tk_zzzzzz_xxyy, \
                             tk_zzzzzz_xxyz, \
                             tk_zzzzzz_xxzz, \
                             tk_zzzzzz_xyyy, \
                             tk_zzzzzz_xyyz, \
                             tk_zzzzzz_xyzz, \
                             tk_zzzzzz_xzzz, \
                             tk_zzzzzz_yyyy, \
                             tk_zzzzzz_yyyz, \
                             tk_zzzzzz_yyzz, \
                             tk_zzzzzz_yzzz, \
                             tk_zzzzzz_zzzz, \
                             ts_zzzz_xxxx,   \
                             ts_zzzz_xxxy,   \
                             ts_zzzz_xxxz,   \
                             ts_zzzz_xxyy,   \
                             ts_zzzz_xxyz,   \
                             ts_zzzz_xxzz,   \
                             ts_zzzz_xyyy,   \
                             ts_zzzz_xyyz,   \
                             ts_zzzz_xyzz,   \
                             ts_zzzz_xzzz,   \
                             ts_zzzz_yyyy,   \
                             ts_zzzz_yyyz,   \
                             ts_zzzz_yyzz,   \
                             ts_zzzz_yzzz,   \
                             ts_zzzz_zzzz,   \
                             ts_zzzzzz_xxxx, \
                             ts_zzzzzz_xxxy, \
                             ts_zzzzzz_xxxz, \
                             ts_zzzzzz_xxyy, \
                             ts_zzzzzz_xxyz, \
                             ts_zzzzzz_xxzz, \
                             ts_zzzzzz_xyyy, \
                             ts_zzzzzz_xyyz, \
                             ts_zzzzzz_xyzz, \
                             ts_zzzzzz_xzzz, \
                             ts_zzzzzz_yyyy, \
                             ts_zzzzzz_yyyz, \
                             ts_zzzzzz_yyzz, \
                             ts_zzzzzz_yzzz, \
                             ts_zzzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzzz_xxxx[i] =
            -10.0 * ts_zzzz_xxxx[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxx[i] * fe_0 + tk_zzzzz_xxxx[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxx[i] * fz_0;

        tk_zzzzzz_xxxy[i] =
            -10.0 * ts_zzzz_xxxy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxy[i] * fe_0 + tk_zzzzz_xxxy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxy[i] * fz_0;

        tk_zzzzzz_xxxz[i] = -10.0 * ts_zzzz_xxxz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxz[i] * fe_0 + tk_zzzzz_xxx[i] * fe_0 +
                            tk_zzzzz_xxxz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxz[i] * fz_0;

        tk_zzzzzz_xxyy[i] =
            -10.0 * ts_zzzz_xxyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyy[i] * fe_0 + tk_zzzzz_xxyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyy[i] * fz_0;

        tk_zzzzzz_xxyz[i] = -10.0 * ts_zzzz_xxyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyz[i] * fe_0 + tk_zzzzz_xxy[i] * fe_0 +
                            tk_zzzzz_xxyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyz[i] * fz_0;

        tk_zzzzzz_xxzz[i] = -10.0 * ts_zzzz_xxzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxzz[i] * fe_0 + 2.0 * tk_zzzzz_xxz[i] * fe_0 +
                            tk_zzzzz_xxzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxzz[i] * fz_0;

        tk_zzzzzz_xyyy[i] =
            -10.0 * ts_zzzz_xyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyy[i] * fe_0 + tk_zzzzz_xyyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyy[i] * fz_0;

        tk_zzzzzz_xyyz[i] = -10.0 * ts_zzzz_xyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyz[i] * fe_0 + tk_zzzzz_xyy[i] * fe_0 +
                            tk_zzzzz_xyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyz[i] * fz_0;

        tk_zzzzzz_xyzz[i] = -10.0 * ts_zzzz_xyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyzz[i] * fe_0 + 2.0 * tk_zzzzz_xyz[i] * fe_0 +
                            tk_zzzzz_xyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyzz[i] * fz_0;

        tk_zzzzzz_xzzz[i] = -10.0 * ts_zzzz_xzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xzzz[i] * fe_0 + 3.0 * tk_zzzzz_xzz[i] * fe_0 +
                            tk_zzzzz_xzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xzzz[i] * fz_0;

        tk_zzzzzz_yyyy[i] =
            -10.0 * ts_zzzz_yyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyy[i] * fe_0 + tk_zzzzz_yyyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyy[i] * fz_0;

        tk_zzzzzz_yyyz[i] = -10.0 * ts_zzzz_yyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyz[i] * fe_0 + tk_zzzzz_yyy[i] * fe_0 +
                            tk_zzzzz_yyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyz[i] * fz_0;

        tk_zzzzzz_yyzz[i] = -10.0 * ts_zzzz_yyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyzz[i] * fe_0 + 2.0 * tk_zzzzz_yyz[i] * fe_0 +
                            tk_zzzzz_yyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyzz[i] * fz_0;

        tk_zzzzzz_yzzz[i] = -10.0 * ts_zzzz_yzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yzzz[i] * fe_0 + 3.0 * tk_zzzzz_yzz[i] * fe_0 +
                            tk_zzzzz_yzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yzzz[i] * fz_0;

        tk_zzzzzz_zzzz[i] = -10.0 * ts_zzzz_zzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_zzzz[i] * fe_0 + 4.0 * tk_zzzzz_zzz[i] * fe_0 +
                            tk_zzzzz_zzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
