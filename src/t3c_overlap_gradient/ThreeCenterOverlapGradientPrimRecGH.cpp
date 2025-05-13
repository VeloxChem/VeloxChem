#include "ThreeCenterOverlapGradientPrimRecGH.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_gh(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_gh,
                              const size_t idx_fh,
                              const size_t idx_gg,
                              const size_t idx_gh,
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

    // Set up components of auxiliary buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_fh + 20);

    auto ts_xxy_xxxxx = pbuffer.data(idx_fh + 21);

    auto ts_xxy_xxxxy = pbuffer.data(idx_fh + 22);

    auto ts_xxy_xxxxz = pbuffer.data(idx_fh + 23);

    auto ts_xxy_xxxyy = pbuffer.data(idx_fh + 24);

    auto ts_xxy_xxxyz = pbuffer.data(idx_fh + 25);

    auto ts_xxy_xxxzz = pbuffer.data(idx_fh + 26);

    auto ts_xxy_xxyyy = pbuffer.data(idx_fh + 27);

    auto ts_xxy_xxyyz = pbuffer.data(idx_fh + 28);

    auto ts_xxy_xxyzz = pbuffer.data(idx_fh + 29);

    auto ts_xxy_xxzzz = pbuffer.data(idx_fh + 30);

    auto ts_xxy_xyyyy = pbuffer.data(idx_fh + 31);

    auto ts_xxy_xyyyz = pbuffer.data(idx_fh + 32);

    auto ts_xxy_xyyzz = pbuffer.data(idx_fh + 33);

    auto ts_xxy_xyzzz = pbuffer.data(idx_fh + 34);

    auto ts_xxy_xzzzz = pbuffer.data(idx_fh + 35);

    auto ts_xxy_yyyyy = pbuffer.data(idx_fh + 36);

    auto ts_xxy_yyyyz = pbuffer.data(idx_fh + 37);

    auto ts_xxy_yyyzz = pbuffer.data(idx_fh + 38);

    auto ts_xxy_yyzzz = pbuffer.data(idx_fh + 39);

    auto ts_xxy_yzzzz = pbuffer.data(idx_fh + 40);

    auto ts_xxy_zzzzz = pbuffer.data(idx_fh + 41);

    auto ts_xxz_xxxxx = pbuffer.data(idx_fh + 42);

    auto ts_xxz_xxxxy = pbuffer.data(idx_fh + 43);

    auto ts_xxz_xxxxz = pbuffer.data(idx_fh + 44);

    auto ts_xxz_xxxyy = pbuffer.data(idx_fh + 45);

    auto ts_xxz_xxxyz = pbuffer.data(idx_fh + 46);

    auto ts_xxz_xxxzz = pbuffer.data(idx_fh + 47);

    auto ts_xxz_xxyyy = pbuffer.data(idx_fh + 48);

    auto ts_xxz_xxyyz = pbuffer.data(idx_fh + 49);

    auto ts_xxz_xxyzz = pbuffer.data(idx_fh + 50);

    auto ts_xxz_xxzzz = pbuffer.data(idx_fh + 51);

    auto ts_xxz_xyyyy = pbuffer.data(idx_fh + 52);

    auto ts_xxz_xyyyz = pbuffer.data(idx_fh + 53);

    auto ts_xxz_xyyzz = pbuffer.data(idx_fh + 54);

    auto ts_xxz_xyzzz = pbuffer.data(idx_fh + 55);

    auto ts_xxz_xzzzz = pbuffer.data(idx_fh + 56);

    auto ts_xxz_yyyyy = pbuffer.data(idx_fh + 57);

    auto ts_xxz_yyyyz = pbuffer.data(idx_fh + 58);

    auto ts_xxz_yyyzz = pbuffer.data(idx_fh + 59);

    auto ts_xxz_yyzzz = pbuffer.data(idx_fh + 60);

    auto ts_xxz_yzzzz = pbuffer.data(idx_fh + 61);

    auto ts_xxz_zzzzz = pbuffer.data(idx_fh + 62);

    auto ts_xyy_xxxxx = pbuffer.data(idx_fh + 63);

    auto ts_xyy_xxxxy = pbuffer.data(idx_fh + 64);

    auto ts_xyy_xxxxz = pbuffer.data(idx_fh + 65);

    auto ts_xyy_xxxyy = pbuffer.data(idx_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_fh + 67);

    auto ts_xyy_xxxzz = pbuffer.data(idx_fh + 68);

    auto ts_xyy_xxyyy = pbuffer.data(idx_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_fh + 71);

    auto ts_xyy_xxzzz = pbuffer.data(idx_fh + 72);

    auto ts_xyy_xyyyy = pbuffer.data(idx_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_fh + 76);

    auto ts_xyy_xzzzz = pbuffer.data(idx_fh + 77);

    auto ts_xyy_yyyyy = pbuffer.data(idx_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_fh + 82);

    auto ts_xyy_zzzzz = pbuffer.data(idx_fh + 83);

    auto ts_xyz_xxxxx = pbuffer.data(idx_fh + 84);

    auto ts_xyz_xxxxy = pbuffer.data(idx_fh + 85);

    auto ts_xyz_xxxxz = pbuffer.data(idx_fh + 86);

    auto ts_xyz_xxxyy = pbuffer.data(idx_fh + 87);

    auto ts_xyz_xxxyz = pbuffer.data(idx_fh + 88);

    auto ts_xyz_xxxzz = pbuffer.data(idx_fh + 89);

    auto ts_xyz_xxyyy = pbuffer.data(idx_fh + 90);

    auto ts_xyz_xxyyz = pbuffer.data(idx_fh + 91);

    auto ts_xyz_xxyzz = pbuffer.data(idx_fh + 92);

    auto ts_xyz_xxzzz = pbuffer.data(idx_fh + 93);

    auto ts_xyz_xyyyy = pbuffer.data(idx_fh + 94);

    auto ts_xyz_xyyyz = pbuffer.data(idx_fh + 95);

    auto ts_xyz_xyyzz = pbuffer.data(idx_fh + 96);

    auto ts_xyz_xyzzz = pbuffer.data(idx_fh + 97);

    auto ts_xyz_xzzzz = pbuffer.data(idx_fh + 98);

    auto ts_xyz_yyyyy = pbuffer.data(idx_fh + 99);

    auto ts_xyz_yyyyz = pbuffer.data(idx_fh + 100);

    auto ts_xyz_yyyzz = pbuffer.data(idx_fh + 101);

    auto ts_xyz_yyzzz = pbuffer.data(idx_fh + 102);

    auto ts_xyz_yzzzz = pbuffer.data(idx_fh + 103);

    auto ts_xyz_zzzzz = pbuffer.data(idx_fh + 104);

    auto ts_xzz_xxxxx = pbuffer.data(idx_fh + 105);

    auto ts_xzz_xxxxy = pbuffer.data(idx_fh + 106);

    auto ts_xzz_xxxxz = pbuffer.data(idx_fh + 107);

    auto ts_xzz_xxxyy = pbuffer.data(idx_fh + 108);

    auto ts_xzz_xxxyz = pbuffer.data(idx_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_fh + 110);

    auto ts_xzz_xxyyy = pbuffer.data(idx_fh + 111);

    auto ts_xzz_xxyyz = pbuffer.data(idx_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_fh + 114);

    auto ts_xzz_xyyyy = pbuffer.data(idx_fh + 115);

    auto ts_xzz_xyyyz = pbuffer.data(idx_fh + 116);

    auto ts_xzz_xyyzz = pbuffer.data(idx_fh + 117);

    auto ts_xzz_xyzzz = pbuffer.data(idx_fh + 118);

    auto ts_xzz_xzzzz = pbuffer.data(idx_fh + 119);

    auto ts_xzz_yyyyy = pbuffer.data(idx_fh + 120);

    auto ts_xzz_yyyyz = pbuffer.data(idx_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_fh + 125);

    auto ts_yyy_xxxxx = pbuffer.data(idx_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_fh + 146);

    auto ts_yyz_xxxxx = pbuffer.data(idx_fh + 147);

    auto ts_yyz_xxxxy = pbuffer.data(idx_fh + 148);

    auto ts_yyz_xxxxz = pbuffer.data(idx_fh + 149);

    auto ts_yyz_xxxyy = pbuffer.data(idx_fh + 150);

    auto ts_yyz_xxxyz = pbuffer.data(idx_fh + 151);

    auto ts_yyz_xxxzz = pbuffer.data(idx_fh + 152);

    auto ts_yyz_xxyyy = pbuffer.data(idx_fh + 153);

    auto ts_yyz_xxyyz = pbuffer.data(idx_fh + 154);

    auto ts_yyz_xxyzz = pbuffer.data(idx_fh + 155);

    auto ts_yyz_xxzzz = pbuffer.data(idx_fh + 156);

    auto ts_yyz_xyyyy = pbuffer.data(idx_fh + 157);

    auto ts_yyz_xyyyz = pbuffer.data(idx_fh + 158);

    auto ts_yyz_xyyzz = pbuffer.data(idx_fh + 159);

    auto ts_yyz_xyzzz = pbuffer.data(idx_fh + 160);

    auto ts_yyz_xzzzz = pbuffer.data(idx_fh + 161);

    auto ts_yyz_yyyyy = pbuffer.data(idx_fh + 162);

    auto ts_yyz_yyyyz = pbuffer.data(idx_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_fh + 167);

    auto ts_yzz_xxxxx = pbuffer.data(idx_fh + 168);

    auto ts_yzz_xxxxy = pbuffer.data(idx_fh + 169);

    auto ts_yzz_xxxxz = pbuffer.data(idx_fh + 170);

    auto ts_yzz_xxxyy = pbuffer.data(idx_fh + 171);

    auto ts_yzz_xxxyz = pbuffer.data(idx_fh + 172);

    auto ts_yzz_xxxzz = pbuffer.data(idx_fh + 173);

    auto ts_yzz_xxyyy = pbuffer.data(idx_fh + 174);

    auto ts_yzz_xxyyz = pbuffer.data(idx_fh + 175);

    auto ts_yzz_xxyzz = pbuffer.data(idx_fh + 176);

    auto ts_yzz_xxzzz = pbuffer.data(idx_fh + 177);

    auto ts_yzz_xyyyy = pbuffer.data(idx_fh + 178);

    auto ts_yzz_xyyyz = pbuffer.data(idx_fh + 179);

    auto ts_yzz_xyyzz = pbuffer.data(idx_fh + 180);

    auto ts_yzz_xyzzz = pbuffer.data(idx_fh + 181);

    auto ts_yzz_xzzzz = pbuffer.data(idx_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_fh + 188);

    auto ts_zzz_xxxxx = pbuffer.data(idx_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_fh + 209);

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

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_gh + 20);

    auto ts_xxxy_xxxxx = pbuffer.data(idx_gh + 21);

    auto ts_xxxy_xxxxy = pbuffer.data(idx_gh + 22);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_gh + 23);

    auto ts_xxxy_xxxyy = pbuffer.data(idx_gh + 24);

    auto ts_xxxy_xxxyz = pbuffer.data(idx_gh + 25);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_gh + 26);

    auto ts_xxxy_xxyyy = pbuffer.data(idx_gh + 27);

    auto ts_xxxy_xxyyz = pbuffer.data(idx_gh + 28);

    auto ts_xxxy_xxyzz = pbuffer.data(idx_gh + 29);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_gh + 30);

    auto ts_xxxy_xyyyy = pbuffer.data(idx_gh + 31);

    auto ts_xxxy_xyyyz = pbuffer.data(idx_gh + 32);

    auto ts_xxxy_xyyzz = pbuffer.data(idx_gh + 33);

    auto ts_xxxy_xyzzz = pbuffer.data(idx_gh + 34);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_gh + 35);

    auto ts_xxxy_yyyyy = pbuffer.data(idx_gh + 36);

    auto ts_xxxy_yyyyz = pbuffer.data(idx_gh + 37);

    auto ts_xxxy_yyyzz = pbuffer.data(idx_gh + 38);

    auto ts_xxxy_yyzzz = pbuffer.data(idx_gh + 39);

    auto ts_xxxy_yzzzz = pbuffer.data(idx_gh + 40);

    auto ts_xxxy_zzzzz = pbuffer.data(idx_gh + 41);

    auto ts_xxxz_xxxxx = pbuffer.data(idx_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_gh + 43);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_gh + 44);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_gh + 45);

    auto ts_xxxz_xxxyz = pbuffer.data(idx_gh + 46);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_gh + 47);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_gh + 48);

    auto ts_xxxz_xxyyz = pbuffer.data(idx_gh + 49);

    auto ts_xxxz_xxyzz = pbuffer.data(idx_gh + 50);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_gh + 51);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_gh + 52);

    auto ts_xxxz_xyyyz = pbuffer.data(idx_gh + 53);

    auto ts_xxxz_xyyzz = pbuffer.data(idx_gh + 54);

    auto ts_xxxz_xyzzz = pbuffer.data(idx_gh + 55);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_gh + 56);

    auto ts_xxxz_yyyyy = pbuffer.data(idx_gh + 57);

    auto ts_xxxz_yyyyz = pbuffer.data(idx_gh + 58);

    auto ts_xxxz_yyyzz = pbuffer.data(idx_gh + 59);

    auto ts_xxxz_yyzzz = pbuffer.data(idx_gh + 60);

    auto ts_xxxz_yzzzz = pbuffer.data(idx_gh + 61);

    auto ts_xxxz_zzzzz = pbuffer.data(idx_gh + 62);

    auto ts_xxyy_xxxxx = pbuffer.data(idx_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_gh + 83);

    auto ts_xxyz_xxxxx = pbuffer.data(idx_gh + 84);

    auto ts_xxyz_xxxxy = pbuffer.data(idx_gh + 85);

    auto ts_xxyz_xxxxz = pbuffer.data(idx_gh + 86);

    auto ts_xxyz_xxxyy = pbuffer.data(idx_gh + 87);

    auto ts_xxyz_xxxyz = pbuffer.data(idx_gh + 88);

    auto ts_xxyz_xxxzz = pbuffer.data(idx_gh + 89);

    auto ts_xxyz_xxyyy = pbuffer.data(idx_gh + 90);

    auto ts_xxyz_xxyyz = pbuffer.data(idx_gh + 91);

    auto ts_xxyz_xxyzz = pbuffer.data(idx_gh + 92);

    auto ts_xxyz_xxzzz = pbuffer.data(idx_gh + 93);

    auto ts_xxyz_xyyyy = pbuffer.data(idx_gh + 94);

    auto ts_xxyz_xyyyz = pbuffer.data(idx_gh + 95);

    auto ts_xxyz_xyyzz = pbuffer.data(idx_gh + 96);

    auto ts_xxyz_xyzzz = pbuffer.data(idx_gh + 97);

    auto ts_xxyz_xzzzz = pbuffer.data(idx_gh + 98);

    auto ts_xxyz_yyyyy = pbuffer.data(idx_gh + 99);

    auto ts_xxyz_yyyyz = pbuffer.data(idx_gh + 100);

    auto ts_xxyz_yyyzz = pbuffer.data(idx_gh + 101);

    auto ts_xxyz_yyzzz = pbuffer.data(idx_gh + 102);

    auto ts_xxyz_yzzzz = pbuffer.data(idx_gh + 103);

    auto ts_xxyz_zzzzz = pbuffer.data(idx_gh + 104);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_gh + 119);

    auto ts_xxzz_yyyyy = pbuffer.data(idx_gh + 120);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_gh + 125);

    auto ts_xyyy_xxxxx = pbuffer.data(idx_gh + 126);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_gh + 127);

    auto ts_xyyy_xxxxz = pbuffer.data(idx_gh + 128);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_gh + 130);

    auto ts_xyyy_xxxzz = pbuffer.data(idx_gh + 131);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_gh + 134);

    auto ts_xyyy_xxzzz = pbuffer.data(idx_gh + 135);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_gh + 139);

    auto ts_xyyy_xzzzz = pbuffer.data(idx_gh + 140);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_gh + 145);

    auto ts_xyyy_zzzzz = pbuffer.data(idx_gh + 146);

    auto ts_xyyz_xxxxx = pbuffer.data(idx_gh + 147);

    auto ts_xyyz_xxxxy = pbuffer.data(idx_gh + 148);

    auto ts_xyyz_xxxxz = pbuffer.data(idx_gh + 149);

    auto ts_xyyz_xxxyy = pbuffer.data(idx_gh + 150);

    auto ts_xyyz_xxxyz = pbuffer.data(idx_gh + 151);

    auto ts_xyyz_xxxzz = pbuffer.data(idx_gh + 152);

    auto ts_xyyz_xxyyy = pbuffer.data(idx_gh + 153);

    auto ts_xyyz_xxyyz = pbuffer.data(idx_gh + 154);

    auto ts_xyyz_xxyzz = pbuffer.data(idx_gh + 155);

    auto ts_xyyz_xxzzz = pbuffer.data(idx_gh + 156);

    auto ts_xyyz_xyyyy = pbuffer.data(idx_gh + 157);

    auto ts_xyyz_xyyyz = pbuffer.data(idx_gh + 158);

    auto ts_xyyz_xyyzz = pbuffer.data(idx_gh + 159);

    auto ts_xyyz_xyzzz = pbuffer.data(idx_gh + 160);

    auto ts_xyyz_xzzzz = pbuffer.data(idx_gh + 161);

    auto ts_xyyz_yyyyy = pbuffer.data(idx_gh + 162);

    auto ts_xyyz_yyyyz = pbuffer.data(idx_gh + 163);

    auto ts_xyyz_yyyzz = pbuffer.data(idx_gh + 164);

    auto ts_xyyz_yyzzz = pbuffer.data(idx_gh + 165);

    auto ts_xyyz_yzzzz = pbuffer.data(idx_gh + 166);

    auto ts_xyyz_zzzzz = pbuffer.data(idx_gh + 167);

    auto ts_xyzz_xxxxx = pbuffer.data(idx_gh + 168);

    auto ts_xyzz_xxxxy = pbuffer.data(idx_gh + 169);

    auto ts_xyzz_xxxxz = pbuffer.data(idx_gh + 170);

    auto ts_xyzz_xxxyy = pbuffer.data(idx_gh + 171);

    auto ts_xyzz_xxxyz = pbuffer.data(idx_gh + 172);

    auto ts_xyzz_xxxzz = pbuffer.data(idx_gh + 173);

    auto ts_xyzz_xxyyy = pbuffer.data(idx_gh + 174);

    auto ts_xyzz_xxyyz = pbuffer.data(idx_gh + 175);

    auto ts_xyzz_xxyzz = pbuffer.data(idx_gh + 176);

    auto ts_xyzz_xxzzz = pbuffer.data(idx_gh + 177);

    auto ts_xyzz_xyyyy = pbuffer.data(idx_gh + 178);

    auto ts_xyzz_xyyyz = pbuffer.data(idx_gh + 179);

    auto ts_xyzz_xyyzz = pbuffer.data(idx_gh + 180);

    auto ts_xyzz_xyzzz = pbuffer.data(idx_gh + 181);

    auto ts_xyzz_xzzzz = pbuffer.data(idx_gh + 182);

    auto ts_xyzz_yyyyy = pbuffer.data(idx_gh + 183);

    auto ts_xyzz_yyyyz = pbuffer.data(idx_gh + 184);

    auto ts_xyzz_yyyzz = pbuffer.data(idx_gh + 185);

    auto ts_xyzz_yyzzz = pbuffer.data(idx_gh + 186);

    auto ts_xyzz_yzzzz = pbuffer.data(idx_gh + 187);

    auto ts_xyzz_zzzzz = pbuffer.data(idx_gh + 188);

    auto ts_xzzz_xxxxx = pbuffer.data(idx_gh + 189);

    auto ts_xzzz_xxxxy = pbuffer.data(idx_gh + 190);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_gh + 191);

    auto ts_xzzz_xxxyy = pbuffer.data(idx_gh + 192);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_gh + 194);

    auto ts_xzzz_xxyyy = pbuffer.data(idx_gh + 195);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_gh + 198);

    auto ts_xzzz_xyyyy = pbuffer.data(idx_gh + 199);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_gh + 203);

    auto ts_xzzz_yyyyy = pbuffer.data(idx_gh + 204);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_gh + 230);

    auto ts_yyyz_xxxxx = pbuffer.data(idx_gh + 231);

    auto ts_yyyz_xxxxy = pbuffer.data(idx_gh + 232);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_gh + 233);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_gh + 234);

    auto ts_yyyz_xxxyz = pbuffer.data(idx_gh + 235);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_gh + 236);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_gh + 237);

    auto ts_yyyz_xxyyz = pbuffer.data(idx_gh + 238);

    auto ts_yyyz_xxyzz = pbuffer.data(idx_gh + 239);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_gh + 240);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_gh + 241);

    auto ts_yyyz_xyyyz = pbuffer.data(idx_gh + 242);

    auto ts_yyyz_xyyzz = pbuffer.data(idx_gh + 243);

    auto ts_yyyz_xyzzz = pbuffer.data(idx_gh + 244);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_gh + 245);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_gh + 246);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_gh + 251);

    auto ts_yyzz_xxxxx = pbuffer.data(idx_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_gh + 272);

    auto ts_yzzz_xxxxx = pbuffer.data(idx_gh + 273);

    auto ts_yzzz_xxxxy = pbuffer.data(idx_gh + 274);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_gh + 275);

    auto ts_yzzz_xxxyy = pbuffer.data(idx_gh + 276);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_gh + 278);

    auto ts_yzzz_xxyyy = pbuffer.data(idx_gh + 279);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_gh + 282);

    auto ts_yzzz_xyyyy = pbuffer.data(idx_gh + 283);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_gh + 314);

    // Set up 0-21 components of targeted buffer : GH

    auto gs_x_xxxx_xxxxx = pbuffer.data(idx_g_gh);

    auto gs_x_xxxx_xxxxy = pbuffer.data(idx_g_gh + 1);

    auto gs_x_xxxx_xxxxz = pbuffer.data(idx_g_gh + 2);

    auto gs_x_xxxx_xxxyy = pbuffer.data(idx_g_gh + 3);

    auto gs_x_xxxx_xxxyz = pbuffer.data(idx_g_gh + 4);

    auto gs_x_xxxx_xxxzz = pbuffer.data(idx_g_gh + 5);

    auto gs_x_xxxx_xxyyy = pbuffer.data(idx_g_gh + 6);

    auto gs_x_xxxx_xxyyz = pbuffer.data(idx_g_gh + 7);

    auto gs_x_xxxx_xxyzz = pbuffer.data(idx_g_gh + 8);

    auto gs_x_xxxx_xxzzz = pbuffer.data(idx_g_gh + 9);

    auto gs_x_xxxx_xyyyy = pbuffer.data(idx_g_gh + 10);

    auto gs_x_xxxx_xyyyz = pbuffer.data(idx_g_gh + 11);

    auto gs_x_xxxx_xyyzz = pbuffer.data(idx_g_gh + 12);

    auto gs_x_xxxx_xyzzz = pbuffer.data(idx_g_gh + 13);

    auto gs_x_xxxx_xzzzz = pbuffer.data(idx_g_gh + 14);

    auto gs_x_xxxx_yyyyy = pbuffer.data(idx_g_gh + 15);

    auto gs_x_xxxx_yyyyz = pbuffer.data(idx_g_gh + 16);

    auto gs_x_xxxx_yyyzz = pbuffer.data(idx_g_gh + 17);

    auto gs_x_xxxx_yyzzz = pbuffer.data(idx_g_gh + 18);

    auto gs_x_xxxx_yzzzz = pbuffer.data(idx_g_gh + 19);

    auto gs_x_xxxx_zzzzz = pbuffer.data(idx_g_gh + 20);

    #pragma omp simd aligned(gc_x, gs_x_xxxx_xxxxx, gs_x_xxxx_xxxxy, gs_x_xxxx_xxxxz, gs_x_xxxx_xxxyy, gs_x_xxxx_xxxyz, gs_x_xxxx_xxxzz, gs_x_xxxx_xxyyy, gs_x_xxxx_xxyyz, gs_x_xxxx_xxyzz, gs_x_xxxx_xxzzz, gs_x_xxxx_xyyyy, gs_x_xxxx_xyyyz, gs_x_xxxx_xyyzz, gs_x_xxxx_xyzzz, gs_x_xxxx_xzzzz, gs_x_xxxx_yyyyy, gs_x_xxxx_yyyyz, gs_x_xxxx_yyyzz, gs_x_xxxx_yyzzz, gs_x_xxxx_yzzzz, gs_x_xxxx_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxzz, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyzz, ts_xxx_xxzzz, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyzz, ts_xxx_xyzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyzz, ts_xxx_yyzzz, ts_xxx_yzzzz, ts_xxx_zzzzz, ts_xxxx_xxxx, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxy, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxz, ts_xxxx_xxxzz, ts_xxxx_xxyy, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyz, ts_xxxx_xxyzz, ts_xxxx_xxzz, ts_xxxx_xxzzz, ts_xxxx_xyyy, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyz, ts_xxxx_xyyzz, ts_xxxx_xyzz, ts_xxxx_xyzzz, ts_xxxx_xzzz, ts_xxxx_xzzzz, ts_xxxx_yyyy, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyz, ts_xxxx_yyyzz, ts_xxxx_yyzz, ts_xxxx_yyzzz, ts_xxxx_yzzz, ts_xxxx_yzzzz, ts_xxxx_zzzz, ts_xxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxx_xxxxx[i] = 8.0 * ts_xxx_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxxy[i] = 8.0 * ts_xxx_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxxz[i] = 8.0 * ts_xxx_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxyy[i] = 8.0 * ts_xxx_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxyz[i] = 8.0 * ts_xxx_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxzz[i] = 8.0 * ts_xxx_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxyyy[i] = 8.0 * ts_xxx_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxyyz[i] = 8.0 * ts_xxx_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxyzz[i] = 8.0 * ts_xxx_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxzzz[i] = 8.0 * ts_xxx_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyyyy[i] = 8.0 * ts_xxx_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyyyz[i] = 8.0 * ts_xxx_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyyzz[i] = 8.0 * ts_xxx_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyzzz[i] = 8.0 * ts_xxx_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xzzzz[i] = 8.0 * ts_xxx_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyyyy[i] = 8.0 * ts_xxx_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyyyz[i] = 8.0 * ts_xxx_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyyzz[i] = 8.0 * ts_xxx_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyzzz[i] = 8.0 * ts_xxx_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yzzzz[i] = 8.0 * ts_xxx_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_zzzzz[i] = 8.0 * ts_xxx_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 21-42 components of targeted buffer : GH

    auto gs_x_xxxy_xxxxx = pbuffer.data(idx_g_gh + 21);

    auto gs_x_xxxy_xxxxy = pbuffer.data(idx_g_gh + 22);

    auto gs_x_xxxy_xxxxz = pbuffer.data(idx_g_gh + 23);

    auto gs_x_xxxy_xxxyy = pbuffer.data(idx_g_gh + 24);

    auto gs_x_xxxy_xxxyz = pbuffer.data(idx_g_gh + 25);

    auto gs_x_xxxy_xxxzz = pbuffer.data(idx_g_gh + 26);

    auto gs_x_xxxy_xxyyy = pbuffer.data(idx_g_gh + 27);

    auto gs_x_xxxy_xxyyz = pbuffer.data(idx_g_gh + 28);

    auto gs_x_xxxy_xxyzz = pbuffer.data(idx_g_gh + 29);

    auto gs_x_xxxy_xxzzz = pbuffer.data(idx_g_gh + 30);

    auto gs_x_xxxy_xyyyy = pbuffer.data(idx_g_gh + 31);

    auto gs_x_xxxy_xyyyz = pbuffer.data(idx_g_gh + 32);

    auto gs_x_xxxy_xyyzz = pbuffer.data(idx_g_gh + 33);

    auto gs_x_xxxy_xyzzz = pbuffer.data(idx_g_gh + 34);

    auto gs_x_xxxy_xzzzz = pbuffer.data(idx_g_gh + 35);

    auto gs_x_xxxy_yyyyy = pbuffer.data(idx_g_gh + 36);

    auto gs_x_xxxy_yyyyz = pbuffer.data(idx_g_gh + 37);

    auto gs_x_xxxy_yyyzz = pbuffer.data(idx_g_gh + 38);

    auto gs_x_xxxy_yyzzz = pbuffer.data(idx_g_gh + 39);

    auto gs_x_xxxy_yzzzz = pbuffer.data(idx_g_gh + 40);

    auto gs_x_xxxy_zzzzz = pbuffer.data(idx_g_gh + 41);

    #pragma omp simd aligned(gc_x, gs_x_xxxy_xxxxx, gs_x_xxxy_xxxxy, gs_x_xxxy_xxxxz, gs_x_xxxy_xxxyy, gs_x_xxxy_xxxyz, gs_x_xxxy_xxxzz, gs_x_xxxy_xxyyy, gs_x_xxxy_xxyyz, gs_x_xxxy_xxyzz, gs_x_xxxy_xxzzz, gs_x_xxxy_xyyyy, gs_x_xxxy_xyyyz, gs_x_xxxy_xyyzz, gs_x_xxxy_xyzzz, gs_x_xxxy_xzzzz, gs_x_xxxy_yyyyy, gs_x_xxxy_yyyyz, gs_x_xxxy_yyyzz, gs_x_xxxy_yyzzz, gs_x_xxxy_yzzzz, gs_x_xxxy_zzzzz, ts_xxxy_xxxx, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxy, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxz, ts_xxxy_xxxzz, ts_xxxy_xxyy, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyz, ts_xxxy_xxyzz, ts_xxxy_xxzz, ts_xxxy_xxzzz, ts_xxxy_xyyy, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyz, ts_xxxy_xyyzz, ts_xxxy_xyzz, ts_xxxy_xyzzz, ts_xxxy_xzzz, ts_xxxy_xzzzz, ts_xxxy_yyyy, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyz, ts_xxxy_yyyzz, ts_xxxy_yyzz, ts_xxxy_yyzzz, ts_xxxy_yzzz, ts_xxxy_yzzzz, ts_xxxy_zzzz, ts_xxxy_zzzzz, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxzz, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyzz, ts_xxy_xxzzz, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyzz, ts_xxy_xyzzz, ts_xxy_xzzzz, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyzz, ts_xxy_yyzzz, ts_xxy_yzzzz, ts_xxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxy_xxxxx[i] = 6.0 * ts_xxy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxxy[i] = 6.0 * ts_xxy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxxz[i] = 6.0 * ts_xxy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxyy[i] = 6.0 * ts_xxy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxyz[i] = 6.0 * ts_xxy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxzz[i] = 6.0 * ts_xxy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxyyy[i] = 6.0 * ts_xxy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxyyz[i] = 6.0 * ts_xxy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxyzz[i] = 6.0 * ts_xxy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxzzz[i] = 6.0 * ts_xxy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyyyy[i] = 6.0 * ts_xxy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyyyz[i] = 6.0 * ts_xxy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyyzz[i] = 6.0 * ts_xxy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyzzz[i] = 6.0 * ts_xxy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xzzzz[i] = 6.0 * ts_xxy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyyyy[i] = 6.0 * ts_xxy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyyyz[i] = 6.0 * ts_xxy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyyzz[i] = 6.0 * ts_xxy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyzzz[i] = 6.0 * ts_xxy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yzzzz[i] = 6.0 * ts_xxy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_zzzzz[i] = 6.0 * ts_xxy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-63 components of targeted buffer : GH

    auto gs_x_xxxz_xxxxx = pbuffer.data(idx_g_gh + 42);

    auto gs_x_xxxz_xxxxy = pbuffer.data(idx_g_gh + 43);

    auto gs_x_xxxz_xxxxz = pbuffer.data(idx_g_gh + 44);

    auto gs_x_xxxz_xxxyy = pbuffer.data(idx_g_gh + 45);

    auto gs_x_xxxz_xxxyz = pbuffer.data(idx_g_gh + 46);

    auto gs_x_xxxz_xxxzz = pbuffer.data(idx_g_gh + 47);

    auto gs_x_xxxz_xxyyy = pbuffer.data(idx_g_gh + 48);

    auto gs_x_xxxz_xxyyz = pbuffer.data(idx_g_gh + 49);

    auto gs_x_xxxz_xxyzz = pbuffer.data(idx_g_gh + 50);

    auto gs_x_xxxz_xxzzz = pbuffer.data(idx_g_gh + 51);

    auto gs_x_xxxz_xyyyy = pbuffer.data(idx_g_gh + 52);

    auto gs_x_xxxz_xyyyz = pbuffer.data(idx_g_gh + 53);

    auto gs_x_xxxz_xyyzz = pbuffer.data(idx_g_gh + 54);

    auto gs_x_xxxz_xyzzz = pbuffer.data(idx_g_gh + 55);

    auto gs_x_xxxz_xzzzz = pbuffer.data(idx_g_gh + 56);

    auto gs_x_xxxz_yyyyy = pbuffer.data(idx_g_gh + 57);

    auto gs_x_xxxz_yyyyz = pbuffer.data(idx_g_gh + 58);

    auto gs_x_xxxz_yyyzz = pbuffer.data(idx_g_gh + 59);

    auto gs_x_xxxz_yyzzz = pbuffer.data(idx_g_gh + 60);

    auto gs_x_xxxz_yzzzz = pbuffer.data(idx_g_gh + 61);

    auto gs_x_xxxz_zzzzz = pbuffer.data(idx_g_gh + 62);

    #pragma omp simd aligned(gc_x, gs_x_xxxz_xxxxx, gs_x_xxxz_xxxxy, gs_x_xxxz_xxxxz, gs_x_xxxz_xxxyy, gs_x_xxxz_xxxyz, gs_x_xxxz_xxxzz, gs_x_xxxz_xxyyy, gs_x_xxxz_xxyyz, gs_x_xxxz_xxyzz, gs_x_xxxz_xxzzz, gs_x_xxxz_xyyyy, gs_x_xxxz_xyyyz, gs_x_xxxz_xyyzz, gs_x_xxxz_xyzzz, gs_x_xxxz_xzzzz, gs_x_xxxz_yyyyy, gs_x_xxxz_yyyyz, gs_x_xxxz_yyyzz, gs_x_xxxz_yyzzz, gs_x_xxxz_yzzzz, gs_x_xxxz_zzzzz, ts_xxxz_xxxx, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxy, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxz, ts_xxxz_xxxzz, ts_xxxz_xxyy, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyz, ts_xxxz_xxyzz, ts_xxxz_xxzz, ts_xxxz_xxzzz, ts_xxxz_xyyy, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyz, ts_xxxz_xyyzz, ts_xxxz_xyzz, ts_xxxz_xyzzz, ts_xxxz_xzzz, ts_xxxz_xzzzz, ts_xxxz_yyyy, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyz, ts_xxxz_yyyzz, ts_xxxz_yyzz, ts_xxxz_yyzzz, ts_xxxz_yzzz, ts_xxxz_yzzzz, ts_xxxz_zzzz, ts_xxxz_zzzzz, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxzz, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyzz, ts_xxz_xxzzz, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyzz, ts_xxz_xyzzz, ts_xxz_xzzzz, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyzz, ts_xxz_yyzzz, ts_xxz_yzzzz, ts_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxz_xxxxx[i] = 6.0 * ts_xxz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxxy[i] = 6.0 * ts_xxz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxxz[i] = 6.0 * ts_xxz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxyy[i] = 6.0 * ts_xxz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxyz[i] = 6.0 * ts_xxz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxzz[i] = 6.0 * ts_xxz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxyyy[i] = 6.0 * ts_xxz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxyyz[i] = 6.0 * ts_xxz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxyzz[i] = 6.0 * ts_xxz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxzzz[i] = 6.0 * ts_xxz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyyyy[i] = 6.0 * ts_xxz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyyyz[i] = 6.0 * ts_xxz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyyzz[i] = 6.0 * ts_xxz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyzzz[i] = 6.0 * ts_xxz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xzzzz[i] = 6.0 * ts_xxz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyyyy[i] = 6.0 * ts_xxz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyyyz[i] = 6.0 * ts_xxz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyyzz[i] = 6.0 * ts_xxz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyzzz[i] = 6.0 * ts_xxz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yzzzz[i] = 6.0 * ts_xxz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_zzzzz[i] = 6.0 * ts_xxz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 63-84 components of targeted buffer : GH

    auto gs_x_xxyy_xxxxx = pbuffer.data(idx_g_gh + 63);

    auto gs_x_xxyy_xxxxy = pbuffer.data(idx_g_gh + 64);

    auto gs_x_xxyy_xxxxz = pbuffer.data(idx_g_gh + 65);

    auto gs_x_xxyy_xxxyy = pbuffer.data(idx_g_gh + 66);

    auto gs_x_xxyy_xxxyz = pbuffer.data(idx_g_gh + 67);

    auto gs_x_xxyy_xxxzz = pbuffer.data(idx_g_gh + 68);

    auto gs_x_xxyy_xxyyy = pbuffer.data(idx_g_gh + 69);

    auto gs_x_xxyy_xxyyz = pbuffer.data(idx_g_gh + 70);

    auto gs_x_xxyy_xxyzz = pbuffer.data(idx_g_gh + 71);

    auto gs_x_xxyy_xxzzz = pbuffer.data(idx_g_gh + 72);

    auto gs_x_xxyy_xyyyy = pbuffer.data(idx_g_gh + 73);

    auto gs_x_xxyy_xyyyz = pbuffer.data(idx_g_gh + 74);

    auto gs_x_xxyy_xyyzz = pbuffer.data(idx_g_gh + 75);

    auto gs_x_xxyy_xyzzz = pbuffer.data(idx_g_gh + 76);

    auto gs_x_xxyy_xzzzz = pbuffer.data(idx_g_gh + 77);

    auto gs_x_xxyy_yyyyy = pbuffer.data(idx_g_gh + 78);

    auto gs_x_xxyy_yyyyz = pbuffer.data(idx_g_gh + 79);

    auto gs_x_xxyy_yyyzz = pbuffer.data(idx_g_gh + 80);

    auto gs_x_xxyy_yyzzz = pbuffer.data(idx_g_gh + 81);

    auto gs_x_xxyy_yzzzz = pbuffer.data(idx_g_gh + 82);

    auto gs_x_xxyy_zzzzz = pbuffer.data(idx_g_gh + 83);

    #pragma omp simd aligned(gc_x, gs_x_xxyy_xxxxx, gs_x_xxyy_xxxxy, gs_x_xxyy_xxxxz, gs_x_xxyy_xxxyy, gs_x_xxyy_xxxyz, gs_x_xxyy_xxxzz, gs_x_xxyy_xxyyy, gs_x_xxyy_xxyyz, gs_x_xxyy_xxyzz, gs_x_xxyy_xxzzz, gs_x_xxyy_xyyyy, gs_x_xxyy_xyyyz, gs_x_xxyy_xyyzz, gs_x_xxyy_xyzzz, gs_x_xxyy_xzzzz, gs_x_xxyy_yyyyy, gs_x_xxyy_yyyyz, gs_x_xxyy_yyyzz, gs_x_xxyy_yyzzz, gs_x_xxyy_yzzzz, gs_x_xxyy_zzzzz, ts_xxyy_xxxx, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxz, ts_xxyy_xxxzz, ts_xxyy_xxyy, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyz, ts_xxyy_xxyzz, ts_xxyy_xxzz, ts_xxyy_xxzzz, ts_xxyy_xyyy, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyz, ts_xxyy_xyyzz, ts_xxyy_xyzz, ts_xxyy_xyzzz, ts_xxyy_xzzz, ts_xxyy_xzzzz, ts_xxyy_yyyy, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyz, ts_xxyy_yyyzz, ts_xxyy_yyzz, ts_xxyy_yyzzz, ts_xxyy_yzzz, ts_xxyy_yzzzz, ts_xxyy_zzzz, ts_xxyy_zzzzz, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxzz, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyzz, ts_xyy_xxzzz, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyzz, ts_xyy_xyzzz, ts_xyy_xzzzz, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyzz, ts_xyy_yyzzz, ts_xyy_yzzzz, ts_xyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyy_xxxxx[i] = 4.0 * ts_xyy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxxy[i] = 4.0 * ts_xyy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxxz[i] = 4.0 * ts_xyy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxyy[i] = 4.0 * ts_xyy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxyz[i] = 4.0 * ts_xyy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxzz[i] = 4.0 * ts_xyy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxyyy[i] = 4.0 * ts_xyy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxyyz[i] = 4.0 * ts_xyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxyzz[i] = 4.0 * ts_xyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxzzz[i] = 4.0 * ts_xyy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyyyy[i] = 4.0 * ts_xyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyyyz[i] = 4.0 * ts_xyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyyzz[i] = 4.0 * ts_xyy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyzzz[i] = 4.0 * ts_xyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xzzzz[i] = 4.0 * ts_xyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyyyy[i] = 4.0 * ts_xyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyyyz[i] = 4.0 * ts_xyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyyzz[i] = 4.0 * ts_xyy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyzzz[i] = 4.0 * ts_xyy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yzzzz[i] = 4.0 * ts_xyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_zzzzz[i] = 4.0 * ts_xyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 84-105 components of targeted buffer : GH

    auto gs_x_xxyz_xxxxx = pbuffer.data(idx_g_gh + 84);

    auto gs_x_xxyz_xxxxy = pbuffer.data(idx_g_gh + 85);

    auto gs_x_xxyz_xxxxz = pbuffer.data(idx_g_gh + 86);

    auto gs_x_xxyz_xxxyy = pbuffer.data(idx_g_gh + 87);

    auto gs_x_xxyz_xxxyz = pbuffer.data(idx_g_gh + 88);

    auto gs_x_xxyz_xxxzz = pbuffer.data(idx_g_gh + 89);

    auto gs_x_xxyz_xxyyy = pbuffer.data(idx_g_gh + 90);

    auto gs_x_xxyz_xxyyz = pbuffer.data(idx_g_gh + 91);

    auto gs_x_xxyz_xxyzz = pbuffer.data(idx_g_gh + 92);

    auto gs_x_xxyz_xxzzz = pbuffer.data(idx_g_gh + 93);

    auto gs_x_xxyz_xyyyy = pbuffer.data(idx_g_gh + 94);

    auto gs_x_xxyz_xyyyz = pbuffer.data(idx_g_gh + 95);

    auto gs_x_xxyz_xyyzz = pbuffer.data(idx_g_gh + 96);

    auto gs_x_xxyz_xyzzz = pbuffer.data(idx_g_gh + 97);

    auto gs_x_xxyz_xzzzz = pbuffer.data(idx_g_gh + 98);

    auto gs_x_xxyz_yyyyy = pbuffer.data(idx_g_gh + 99);

    auto gs_x_xxyz_yyyyz = pbuffer.data(idx_g_gh + 100);

    auto gs_x_xxyz_yyyzz = pbuffer.data(idx_g_gh + 101);

    auto gs_x_xxyz_yyzzz = pbuffer.data(idx_g_gh + 102);

    auto gs_x_xxyz_yzzzz = pbuffer.data(idx_g_gh + 103);

    auto gs_x_xxyz_zzzzz = pbuffer.data(idx_g_gh + 104);

    #pragma omp simd aligned(gc_x, gs_x_xxyz_xxxxx, gs_x_xxyz_xxxxy, gs_x_xxyz_xxxxz, gs_x_xxyz_xxxyy, gs_x_xxyz_xxxyz, gs_x_xxyz_xxxzz, gs_x_xxyz_xxyyy, gs_x_xxyz_xxyyz, gs_x_xxyz_xxyzz, gs_x_xxyz_xxzzz, gs_x_xxyz_xyyyy, gs_x_xxyz_xyyyz, gs_x_xxyz_xyyzz, gs_x_xxyz_xyzzz, gs_x_xxyz_xzzzz, gs_x_xxyz_yyyyy, gs_x_xxyz_yyyyz, gs_x_xxyz_yyyzz, gs_x_xxyz_yyzzz, gs_x_xxyz_yzzzz, gs_x_xxyz_zzzzz, ts_xxyz_xxxx, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxy, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxz, ts_xxyz_xxxzz, ts_xxyz_xxyy, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyz, ts_xxyz_xxyzz, ts_xxyz_xxzz, ts_xxyz_xxzzz, ts_xxyz_xyyy, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyz, ts_xxyz_xyyzz, ts_xxyz_xyzz, ts_xxyz_xyzzz, ts_xxyz_xzzz, ts_xxyz_xzzzz, ts_xxyz_yyyy, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyz, ts_xxyz_yyyzz, ts_xxyz_yyzz, ts_xxyz_yyzzz, ts_xxyz_yzzz, ts_xxyz_yzzzz, ts_xxyz_zzzz, ts_xxyz_zzzzz, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxzz, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyzz, ts_xyz_xxzzz, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyzz, ts_xyz_xyzzz, ts_xyz_xzzzz, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyzz, ts_xyz_yyzzz, ts_xyz_yzzzz, ts_xyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyz_xxxxx[i] = 4.0 * ts_xyz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxxy[i] = 4.0 * ts_xyz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxxz[i] = 4.0 * ts_xyz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxyy[i] = 4.0 * ts_xyz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxyz[i] = 4.0 * ts_xyz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxzz[i] = 4.0 * ts_xyz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxyyy[i] = 4.0 * ts_xyz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxyyz[i] = 4.0 * ts_xyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxyzz[i] = 4.0 * ts_xyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxzzz[i] = 4.0 * ts_xyz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyyyy[i] = 4.0 * ts_xyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyyyz[i] = 4.0 * ts_xyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyyzz[i] = 4.0 * ts_xyz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyzzz[i] = 4.0 * ts_xyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xzzzz[i] = 4.0 * ts_xyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyyyy[i] = 4.0 * ts_xyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyyyz[i] = 4.0 * ts_xyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyyzz[i] = 4.0 * ts_xyz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyzzz[i] = 4.0 * ts_xyz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yzzzz[i] = 4.0 * ts_xyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_zzzzz[i] = 4.0 * ts_xyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-126 components of targeted buffer : GH

    auto gs_x_xxzz_xxxxx = pbuffer.data(idx_g_gh + 105);

    auto gs_x_xxzz_xxxxy = pbuffer.data(idx_g_gh + 106);

    auto gs_x_xxzz_xxxxz = pbuffer.data(idx_g_gh + 107);

    auto gs_x_xxzz_xxxyy = pbuffer.data(idx_g_gh + 108);

    auto gs_x_xxzz_xxxyz = pbuffer.data(idx_g_gh + 109);

    auto gs_x_xxzz_xxxzz = pbuffer.data(idx_g_gh + 110);

    auto gs_x_xxzz_xxyyy = pbuffer.data(idx_g_gh + 111);

    auto gs_x_xxzz_xxyyz = pbuffer.data(idx_g_gh + 112);

    auto gs_x_xxzz_xxyzz = pbuffer.data(idx_g_gh + 113);

    auto gs_x_xxzz_xxzzz = pbuffer.data(idx_g_gh + 114);

    auto gs_x_xxzz_xyyyy = pbuffer.data(idx_g_gh + 115);

    auto gs_x_xxzz_xyyyz = pbuffer.data(idx_g_gh + 116);

    auto gs_x_xxzz_xyyzz = pbuffer.data(idx_g_gh + 117);

    auto gs_x_xxzz_xyzzz = pbuffer.data(idx_g_gh + 118);

    auto gs_x_xxzz_xzzzz = pbuffer.data(idx_g_gh + 119);

    auto gs_x_xxzz_yyyyy = pbuffer.data(idx_g_gh + 120);

    auto gs_x_xxzz_yyyyz = pbuffer.data(idx_g_gh + 121);

    auto gs_x_xxzz_yyyzz = pbuffer.data(idx_g_gh + 122);

    auto gs_x_xxzz_yyzzz = pbuffer.data(idx_g_gh + 123);

    auto gs_x_xxzz_yzzzz = pbuffer.data(idx_g_gh + 124);

    auto gs_x_xxzz_zzzzz = pbuffer.data(idx_g_gh + 125);

    #pragma omp simd aligned(gc_x, gs_x_xxzz_xxxxx, gs_x_xxzz_xxxxy, gs_x_xxzz_xxxxz, gs_x_xxzz_xxxyy, gs_x_xxzz_xxxyz, gs_x_xxzz_xxxzz, gs_x_xxzz_xxyyy, gs_x_xxzz_xxyyz, gs_x_xxzz_xxyzz, gs_x_xxzz_xxzzz, gs_x_xxzz_xyyyy, gs_x_xxzz_xyyyz, gs_x_xxzz_xyyzz, gs_x_xxzz_xyzzz, gs_x_xxzz_xzzzz, gs_x_xxzz_yyyyy, gs_x_xxzz_yyyyz, gs_x_xxzz_yyyzz, gs_x_xxzz_yyzzz, gs_x_xxzz_yzzzz, gs_x_xxzz_zzzzz, ts_xxzz_xxxx, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxy, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxz, ts_xxzz_xxxzz, ts_xxzz_xxyy, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyz, ts_xxzz_xxyzz, ts_xxzz_xxzz, ts_xxzz_xxzzz, ts_xxzz_xyyy, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyz, ts_xxzz_xyyzz, ts_xxzz_xyzz, ts_xxzz_xyzzz, ts_xxzz_xzzz, ts_xxzz_xzzzz, ts_xxzz_yyyy, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyz, ts_xxzz_yyyzz, ts_xxzz_yyzz, ts_xxzz_yyzzz, ts_xxzz_yzzz, ts_xxzz_yzzzz, ts_xxzz_zzzz, ts_xxzz_zzzzz, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxzz, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyzz, ts_xzz_xxzzz, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyzz, ts_xzz_xyzzz, ts_xzz_xzzzz, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyzz, ts_xzz_yyzzz, ts_xzz_yzzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzz_xxxxx[i] = 4.0 * ts_xzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxxy[i] = 4.0 * ts_xzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxxz[i] = 4.0 * ts_xzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxyy[i] = 4.0 * ts_xzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxyz[i] = 4.0 * ts_xzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxzz[i] = 4.0 * ts_xzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxyyy[i] = 4.0 * ts_xzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxyyz[i] = 4.0 * ts_xzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxyzz[i] = 4.0 * ts_xzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxzzz[i] = 4.0 * ts_xzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyyyy[i] = 4.0 * ts_xzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyyyz[i] = 4.0 * ts_xzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyyzz[i] = 4.0 * ts_xzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyzzz[i] = 4.0 * ts_xzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xzzzz[i] = 4.0 * ts_xzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyyyy[i] = 4.0 * ts_xzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyyyz[i] = 4.0 * ts_xzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyyzz[i] = 4.0 * ts_xzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyzzz[i] = 4.0 * ts_xzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yzzzz[i] = 4.0 * ts_xzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_zzzzz[i] = 4.0 * ts_xzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 126-147 components of targeted buffer : GH

    auto gs_x_xyyy_xxxxx = pbuffer.data(idx_g_gh + 126);

    auto gs_x_xyyy_xxxxy = pbuffer.data(idx_g_gh + 127);

    auto gs_x_xyyy_xxxxz = pbuffer.data(idx_g_gh + 128);

    auto gs_x_xyyy_xxxyy = pbuffer.data(idx_g_gh + 129);

    auto gs_x_xyyy_xxxyz = pbuffer.data(idx_g_gh + 130);

    auto gs_x_xyyy_xxxzz = pbuffer.data(idx_g_gh + 131);

    auto gs_x_xyyy_xxyyy = pbuffer.data(idx_g_gh + 132);

    auto gs_x_xyyy_xxyyz = pbuffer.data(idx_g_gh + 133);

    auto gs_x_xyyy_xxyzz = pbuffer.data(idx_g_gh + 134);

    auto gs_x_xyyy_xxzzz = pbuffer.data(idx_g_gh + 135);

    auto gs_x_xyyy_xyyyy = pbuffer.data(idx_g_gh + 136);

    auto gs_x_xyyy_xyyyz = pbuffer.data(idx_g_gh + 137);

    auto gs_x_xyyy_xyyzz = pbuffer.data(idx_g_gh + 138);

    auto gs_x_xyyy_xyzzz = pbuffer.data(idx_g_gh + 139);

    auto gs_x_xyyy_xzzzz = pbuffer.data(idx_g_gh + 140);

    auto gs_x_xyyy_yyyyy = pbuffer.data(idx_g_gh + 141);

    auto gs_x_xyyy_yyyyz = pbuffer.data(idx_g_gh + 142);

    auto gs_x_xyyy_yyyzz = pbuffer.data(idx_g_gh + 143);

    auto gs_x_xyyy_yyzzz = pbuffer.data(idx_g_gh + 144);

    auto gs_x_xyyy_yzzzz = pbuffer.data(idx_g_gh + 145);

    auto gs_x_xyyy_zzzzz = pbuffer.data(idx_g_gh + 146);

    #pragma omp simd aligned(gc_x, gs_x_xyyy_xxxxx, gs_x_xyyy_xxxxy, gs_x_xyyy_xxxxz, gs_x_xyyy_xxxyy, gs_x_xyyy_xxxyz, gs_x_xyyy_xxxzz, gs_x_xyyy_xxyyy, gs_x_xyyy_xxyyz, gs_x_xyyy_xxyzz, gs_x_xyyy_xxzzz, gs_x_xyyy_xyyyy, gs_x_xyyy_xyyyz, gs_x_xyyy_xyyzz, gs_x_xyyy_xyzzz, gs_x_xyyy_xzzzz, gs_x_xyyy_yyyyy, gs_x_xyyy_yyyyz, gs_x_xyyy_yyyzz, gs_x_xyyy_yyzzz, gs_x_xyyy_yzzzz, gs_x_xyyy_zzzzz, ts_xyyy_xxxx, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxy, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxz, ts_xyyy_xxxzz, ts_xyyy_xxyy, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyz, ts_xyyy_xxyzz, ts_xyyy_xxzz, ts_xyyy_xxzzz, ts_xyyy_xyyy, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyz, ts_xyyy_xyyzz, ts_xyyy_xyzz, ts_xyyy_xyzzz, ts_xyyy_xzzz, ts_xyyy_xzzzz, ts_xyyy_yyyy, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyz, ts_xyyy_yyyzz, ts_xyyy_yyzz, ts_xyyy_yyzzz, ts_xyyy_yzzz, ts_xyyy_yzzzz, ts_xyyy_zzzz, ts_xyyy_zzzzz, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxzz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xxzzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_xzzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyy_xxxxx[i] = 2.0 * ts_yyy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxxy[i] = 2.0 * ts_yyy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxxz[i] = 2.0 * ts_yyy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxyy[i] = 2.0 * ts_yyy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxyz[i] = 2.0 * ts_yyy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxzz[i] = 2.0 * ts_yyy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxyyy[i] = 2.0 * ts_yyy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxyyz[i] = 2.0 * ts_yyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxyzz[i] = 2.0 * ts_yyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxzzz[i] = 2.0 * ts_yyy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyyyy[i] = 2.0 * ts_yyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyyyz[i] = 2.0 * ts_yyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyyzz[i] = 2.0 * ts_yyy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyzzz[i] = 2.0 * ts_yyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xzzzz[i] = 2.0 * ts_yyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyyyy[i] = 2.0 * ts_yyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyyyz[i] = 2.0 * ts_yyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyyzz[i] = 2.0 * ts_yyy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyzzz[i] = 2.0 * ts_yyy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yzzzz[i] = 2.0 * ts_yyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_zzzzz[i] = 2.0 * ts_yyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 147-168 components of targeted buffer : GH

    auto gs_x_xyyz_xxxxx = pbuffer.data(idx_g_gh + 147);

    auto gs_x_xyyz_xxxxy = pbuffer.data(idx_g_gh + 148);

    auto gs_x_xyyz_xxxxz = pbuffer.data(idx_g_gh + 149);

    auto gs_x_xyyz_xxxyy = pbuffer.data(idx_g_gh + 150);

    auto gs_x_xyyz_xxxyz = pbuffer.data(idx_g_gh + 151);

    auto gs_x_xyyz_xxxzz = pbuffer.data(idx_g_gh + 152);

    auto gs_x_xyyz_xxyyy = pbuffer.data(idx_g_gh + 153);

    auto gs_x_xyyz_xxyyz = pbuffer.data(idx_g_gh + 154);

    auto gs_x_xyyz_xxyzz = pbuffer.data(idx_g_gh + 155);

    auto gs_x_xyyz_xxzzz = pbuffer.data(idx_g_gh + 156);

    auto gs_x_xyyz_xyyyy = pbuffer.data(idx_g_gh + 157);

    auto gs_x_xyyz_xyyyz = pbuffer.data(idx_g_gh + 158);

    auto gs_x_xyyz_xyyzz = pbuffer.data(idx_g_gh + 159);

    auto gs_x_xyyz_xyzzz = pbuffer.data(idx_g_gh + 160);

    auto gs_x_xyyz_xzzzz = pbuffer.data(idx_g_gh + 161);

    auto gs_x_xyyz_yyyyy = pbuffer.data(idx_g_gh + 162);

    auto gs_x_xyyz_yyyyz = pbuffer.data(idx_g_gh + 163);

    auto gs_x_xyyz_yyyzz = pbuffer.data(idx_g_gh + 164);

    auto gs_x_xyyz_yyzzz = pbuffer.data(idx_g_gh + 165);

    auto gs_x_xyyz_yzzzz = pbuffer.data(idx_g_gh + 166);

    auto gs_x_xyyz_zzzzz = pbuffer.data(idx_g_gh + 167);

    #pragma omp simd aligned(gc_x, gs_x_xyyz_xxxxx, gs_x_xyyz_xxxxy, gs_x_xyyz_xxxxz, gs_x_xyyz_xxxyy, gs_x_xyyz_xxxyz, gs_x_xyyz_xxxzz, gs_x_xyyz_xxyyy, gs_x_xyyz_xxyyz, gs_x_xyyz_xxyzz, gs_x_xyyz_xxzzz, gs_x_xyyz_xyyyy, gs_x_xyyz_xyyyz, gs_x_xyyz_xyyzz, gs_x_xyyz_xyzzz, gs_x_xyyz_xzzzz, gs_x_xyyz_yyyyy, gs_x_xyyz_yyyyz, gs_x_xyyz_yyyzz, gs_x_xyyz_yyzzz, gs_x_xyyz_yzzzz, gs_x_xyyz_zzzzz, ts_xyyz_xxxx, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxy, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxz, ts_xyyz_xxxzz, ts_xyyz_xxyy, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyz, ts_xyyz_xxyzz, ts_xyyz_xxzz, ts_xyyz_xxzzz, ts_xyyz_xyyy, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyz, ts_xyyz_xyyzz, ts_xyyz_xyzz, ts_xyyz_xyzzz, ts_xyyz_xzzz, ts_xyyz_xzzzz, ts_xyyz_yyyy, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyz, ts_xyyz_yyyzz, ts_xyyz_yyzz, ts_xyyz_yyzzz, ts_xyyz_yzzz, ts_xyyz_yzzzz, ts_xyyz_zzzz, ts_xyyz_zzzzz, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxzz, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyzz, ts_yyz_xxzzz, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyzz, ts_yyz_xyzzz, ts_yyz_xzzzz, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyzz, ts_yyz_yyzzz, ts_yyz_yzzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyz_xxxxx[i] = 2.0 * ts_yyz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxxy[i] = 2.0 * ts_yyz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxxz[i] = 2.0 * ts_yyz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxyy[i] = 2.0 * ts_yyz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxyz[i] = 2.0 * ts_yyz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxzz[i] = 2.0 * ts_yyz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxyyy[i] = 2.0 * ts_yyz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxyyz[i] = 2.0 * ts_yyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxyzz[i] = 2.0 * ts_yyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxzzz[i] = 2.0 * ts_yyz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyyyy[i] = 2.0 * ts_yyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyyyz[i] = 2.0 * ts_yyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyyzz[i] = 2.0 * ts_yyz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyzzz[i] = 2.0 * ts_yyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xzzzz[i] = 2.0 * ts_yyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyyyy[i] = 2.0 * ts_yyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyyyz[i] = 2.0 * ts_yyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyyzz[i] = 2.0 * ts_yyz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyzzz[i] = 2.0 * ts_yyz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yzzzz[i] = 2.0 * ts_yyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_zzzzz[i] = 2.0 * ts_yyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 168-189 components of targeted buffer : GH

    auto gs_x_xyzz_xxxxx = pbuffer.data(idx_g_gh + 168);

    auto gs_x_xyzz_xxxxy = pbuffer.data(idx_g_gh + 169);

    auto gs_x_xyzz_xxxxz = pbuffer.data(idx_g_gh + 170);

    auto gs_x_xyzz_xxxyy = pbuffer.data(idx_g_gh + 171);

    auto gs_x_xyzz_xxxyz = pbuffer.data(idx_g_gh + 172);

    auto gs_x_xyzz_xxxzz = pbuffer.data(idx_g_gh + 173);

    auto gs_x_xyzz_xxyyy = pbuffer.data(idx_g_gh + 174);

    auto gs_x_xyzz_xxyyz = pbuffer.data(idx_g_gh + 175);

    auto gs_x_xyzz_xxyzz = pbuffer.data(idx_g_gh + 176);

    auto gs_x_xyzz_xxzzz = pbuffer.data(idx_g_gh + 177);

    auto gs_x_xyzz_xyyyy = pbuffer.data(idx_g_gh + 178);

    auto gs_x_xyzz_xyyyz = pbuffer.data(idx_g_gh + 179);

    auto gs_x_xyzz_xyyzz = pbuffer.data(idx_g_gh + 180);

    auto gs_x_xyzz_xyzzz = pbuffer.data(idx_g_gh + 181);

    auto gs_x_xyzz_xzzzz = pbuffer.data(idx_g_gh + 182);

    auto gs_x_xyzz_yyyyy = pbuffer.data(idx_g_gh + 183);

    auto gs_x_xyzz_yyyyz = pbuffer.data(idx_g_gh + 184);

    auto gs_x_xyzz_yyyzz = pbuffer.data(idx_g_gh + 185);

    auto gs_x_xyzz_yyzzz = pbuffer.data(idx_g_gh + 186);

    auto gs_x_xyzz_yzzzz = pbuffer.data(idx_g_gh + 187);

    auto gs_x_xyzz_zzzzz = pbuffer.data(idx_g_gh + 188);

    #pragma omp simd aligned(gc_x, gs_x_xyzz_xxxxx, gs_x_xyzz_xxxxy, gs_x_xyzz_xxxxz, gs_x_xyzz_xxxyy, gs_x_xyzz_xxxyz, gs_x_xyzz_xxxzz, gs_x_xyzz_xxyyy, gs_x_xyzz_xxyyz, gs_x_xyzz_xxyzz, gs_x_xyzz_xxzzz, gs_x_xyzz_xyyyy, gs_x_xyzz_xyyyz, gs_x_xyzz_xyyzz, gs_x_xyzz_xyzzz, gs_x_xyzz_xzzzz, gs_x_xyzz_yyyyy, gs_x_xyzz_yyyyz, gs_x_xyzz_yyyzz, gs_x_xyzz_yyzzz, gs_x_xyzz_yzzzz, gs_x_xyzz_zzzzz, ts_xyzz_xxxx, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxy, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxz, ts_xyzz_xxxzz, ts_xyzz_xxyy, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyz, ts_xyzz_xxyzz, ts_xyzz_xxzz, ts_xyzz_xxzzz, ts_xyzz_xyyy, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyz, ts_xyzz_xyyzz, ts_xyzz_xyzz, ts_xyzz_xyzzz, ts_xyzz_xzzz, ts_xyzz_xzzzz, ts_xyzz_yyyy, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyz, ts_xyzz_yyyzz, ts_xyzz_yyzz, ts_xyzz_yyzzz, ts_xyzz_yzzz, ts_xyzz_yzzzz, ts_xyzz_zzzz, ts_xyzz_zzzzz, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxzz, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyzz, ts_yzz_xxzzz, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyzz, ts_yzz_xyzzz, ts_yzz_xzzzz, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzz_xxxxx[i] = 2.0 * ts_yzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxxy[i] = 2.0 * ts_yzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxxz[i] = 2.0 * ts_yzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxyy[i] = 2.0 * ts_yzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxyz[i] = 2.0 * ts_yzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxzz[i] = 2.0 * ts_yzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxyyy[i] = 2.0 * ts_yzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxyyz[i] = 2.0 * ts_yzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxyzz[i] = 2.0 * ts_yzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxzzz[i] = 2.0 * ts_yzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyyyy[i] = 2.0 * ts_yzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyyyz[i] = 2.0 * ts_yzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyyzz[i] = 2.0 * ts_yzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyzzz[i] = 2.0 * ts_yzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xzzzz[i] = 2.0 * ts_yzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyyyy[i] = 2.0 * ts_yzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyyyz[i] = 2.0 * ts_yzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyyzz[i] = 2.0 * ts_yzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyzzz[i] = 2.0 * ts_yzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yzzzz[i] = 2.0 * ts_yzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_zzzzz[i] = 2.0 * ts_yzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 189-210 components of targeted buffer : GH

    auto gs_x_xzzz_xxxxx = pbuffer.data(idx_g_gh + 189);

    auto gs_x_xzzz_xxxxy = pbuffer.data(idx_g_gh + 190);

    auto gs_x_xzzz_xxxxz = pbuffer.data(idx_g_gh + 191);

    auto gs_x_xzzz_xxxyy = pbuffer.data(idx_g_gh + 192);

    auto gs_x_xzzz_xxxyz = pbuffer.data(idx_g_gh + 193);

    auto gs_x_xzzz_xxxzz = pbuffer.data(idx_g_gh + 194);

    auto gs_x_xzzz_xxyyy = pbuffer.data(idx_g_gh + 195);

    auto gs_x_xzzz_xxyyz = pbuffer.data(idx_g_gh + 196);

    auto gs_x_xzzz_xxyzz = pbuffer.data(idx_g_gh + 197);

    auto gs_x_xzzz_xxzzz = pbuffer.data(idx_g_gh + 198);

    auto gs_x_xzzz_xyyyy = pbuffer.data(idx_g_gh + 199);

    auto gs_x_xzzz_xyyyz = pbuffer.data(idx_g_gh + 200);

    auto gs_x_xzzz_xyyzz = pbuffer.data(idx_g_gh + 201);

    auto gs_x_xzzz_xyzzz = pbuffer.data(idx_g_gh + 202);

    auto gs_x_xzzz_xzzzz = pbuffer.data(idx_g_gh + 203);

    auto gs_x_xzzz_yyyyy = pbuffer.data(idx_g_gh + 204);

    auto gs_x_xzzz_yyyyz = pbuffer.data(idx_g_gh + 205);

    auto gs_x_xzzz_yyyzz = pbuffer.data(idx_g_gh + 206);

    auto gs_x_xzzz_yyzzz = pbuffer.data(idx_g_gh + 207);

    auto gs_x_xzzz_yzzzz = pbuffer.data(idx_g_gh + 208);

    auto gs_x_xzzz_zzzzz = pbuffer.data(idx_g_gh + 209);

    #pragma omp simd aligned(gc_x, gs_x_xzzz_xxxxx, gs_x_xzzz_xxxxy, gs_x_xzzz_xxxxz, gs_x_xzzz_xxxyy, gs_x_xzzz_xxxyz, gs_x_xzzz_xxxzz, gs_x_xzzz_xxyyy, gs_x_xzzz_xxyyz, gs_x_xzzz_xxyzz, gs_x_xzzz_xxzzz, gs_x_xzzz_xyyyy, gs_x_xzzz_xyyyz, gs_x_xzzz_xyyzz, gs_x_xzzz_xyzzz, gs_x_xzzz_xzzzz, gs_x_xzzz_yyyyy, gs_x_xzzz_yyyyz, gs_x_xzzz_yyyzz, gs_x_xzzz_yyzzz, gs_x_xzzz_yzzzz, gs_x_xzzz_zzzzz, ts_xzzz_xxxx, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxy, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxz, ts_xzzz_xxxzz, ts_xzzz_xxyy, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyz, ts_xzzz_xxyzz, ts_xzzz_xxzz, ts_xzzz_xxzzz, ts_xzzz_xyyy, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyz, ts_xzzz_xyyzz, ts_xzzz_xyzz, ts_xzzz_xyzzz, ts_xzzz_xzzz, ts_xzzz_xzzzz, ts_xzzz_yyyy, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyz, ts_xzzz_yyyzz, ts_xzzz_yyzz, ts_xzzz_yyzzz, ts_xzzz_yzzz, ts_xzzz_yzzzz, ts_xzzz_zzzz, ts_xzzz_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzz_xxxxx[i] = 2.0 * ts_zzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxxy[i] = 2.0 * ts_zzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxxz[i] = 2.0 * ts_zzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxyy[i] = 2.0 * ts_zzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxyz[i] = 2.0 * ts_zzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxzz[i] = 2.0 * ts_zzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxyyy[i] = 2.0 * ts_zzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxyyz[i] = 2.0 * ts_zzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxyzz[i] = 2.0 * ts_zzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxzzz[i] = 2.0 * ts_zzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyyyy[i] = 2.0 * ts_zzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyyyz[i] = 2.0 * ts_zzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyyzz[i] = 2.0 * ts_zzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyzzz[i] = 2.0 * ts_zzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xzzzz[i] = 2.0 * ts_zzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyyyy[i] = 2.0 * ts_zzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyyyz[i] = 2.0 * ts_zzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyyzz[i] = 2.0 * ts_zzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyzzz[i] = 2.0 * ts_zzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yzzzz[i] = 2.0 * ts_zzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_zzzzz[i] = 2.0 * ts_zzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 210-231 components of targeted buffer : GH

    auto gs_x_yyyy_xxxxx = pbuffer.data(idx_g_gh + 210);

    auto gs_x_yyyy_xxxxy = pbuffer.data(idx_g_gh + 211);

    auto gs_x_yyyy_xxxxz = pbuffer.data(idx_g_gh + 212);

    auto gs_x_yyyy_xxxyy = pbuffer.data(idx_g_gh + 213);

    auto gs_x_yyyy_xxxyz = pbuffer.data(idx_g_gh + 214);

    auto gs_x_yyyy_xxxzz = pbuffer.data(idx_g_gh + 215);

    auto gs_x_yyyy_xxyyy = pbuffer.data(idx_g_gh + 216);

    auto gs_x_yyyy_xxyyz = pbuffer.data(idx_g_gh + 217);

    auto gs_x_yyyy_xxyzz = pbuffer.data(idx_g_gh + 218);

    auto gs_x_yyyy_xxzzz = pbuffer.data(idx_g_gh + 219);

    auto gs_x_yyyy_xyyyy = pbuffer.data(idx_g_gh + 220);

    auto gs_x_yyyy_xyyyz = pbuffer.data(idx_g_gh + 221);

    auto gs_x_yyyy_xyyzz = pbuffer.data(idx_g_gh + 222);

    auto gs_x_yyyy_xyzzz = pbuffer.data(idx_g_gh + 223);

    auto gs_x_yyyy_xzzzz = pbuffer.data(idx_g_gh + 224);

    auto gs_x_yyyy_yyyyy = pbuffer.data(idx_g_gh + 225);

    auto gs_x_yyyy_yyyyz = pbuffer.data(idx_g_gh + 226);

    auto gs_x_yyyy_yyyzz = pbuffer.data(idx_g_gh + 227);

    auto gs_x_yyyy_yyzzz = pbuffer.data(idx_g_gh + 228);

    auto gs_x_yyyy_yzzzz = pbuffer.data(idx_g_gh + 229);

    auto gs_x_yyyy_zzzzz = pbuffer.data(idx_g_gh + 230);

    #pragma omp simd aligned(gc_x, gs_x_yyyy_xxxxx, gs_x_yyyy_xxxxy, gs_x_yyyy_xxxxz, gs_x_yyyy_xxxyy, gs_x_yyyy_xxxyz, gs_x_yyyy_xxxzz, gs_x_yyyy_xxyyy, gs_x_yyyy_xxyyz, gs_x_yyyy_xxyzz, gs_x_yyyy_xxzzz, gs_x_yyyy_xyyyy, gs_x_yyyy_xyyyz, gs_x_yyyy_xyyzz, gs_x_yyyy_xyzzz, gs_x_yyyy_xzzzz, gs_x_yyyy_yyyyy, gs_x_yyyy_yyyyz, gs_x_yyyy_yyyzz, gs_x_yyyy_yyzzz, gs_x_yyyy_yzzzz, gs_x_yyyy_zzzzz, ts_yyyy_xxxx, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxz, ts_yyyy_xxxzz, ts_yyyy_xxyy, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyz, ts_yyyy_xxyzz, ts_yyyy_xxzz, ts_yyyy_xxzzz, ts_yyyy_xyyy, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyz, ts_yyyy_xyyzz, ts_yyyy_xyzz, ts_yyyy_xyzzz, ts_yyyy_xzzz, ts_yyyy_xzzzz, ts_yyyy_yyyy, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyz, ts_yyyy_yyyzz, ts_yyyy_yyzz, ts_yyyy_yyzzz, ts_yyyy_yzzz, ts_yyyy_yzzzz, ts_yyyy_zzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyy_xxxxx[i] = 10.0 * ts_yyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxxy[i] = 8.0 * ts_yyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxxz[i] = 8.0 * ts_yyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxyy[i] = 6.0 * ts_yyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxyz[i] = 6.0 * ts_yyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxzz[i] = 6.0 * ts_yyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxyyy[i] = 4.0 * ts_yyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxyyz[i] = 4.0 * ts_yyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxyzz[i] = 4.0 * ts_yyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxzzz[i] = 4.0 * ts_yyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyyyy[i] = 2.0 * ts_yyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyyyz[i] = 2.0 * ts_yyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyyzz[i] = 2.0 * ts_yyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyzzz[i] = 2.0 * ts_yyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xzzzz[i] = 2.0 * ts_yyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyyyy[i] = 2.0 * ts_yyyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyyyz[i] = 2.0 * ts_yyyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyyzz[i] = 2.0 * ts_yyyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyzzz[i] = 2.0 * ts_yyyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yzzzz[i] = 2.0 * ts_yyyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_zzzzz[i] = 2.0 * ts_yyyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 231-252 components of targeted buffer : GH

    auto gs_x_yyyz_xxxxx = pbuffer.data(idx_g_gh + 231);

    auto gs_x_yyyz_xxxxy = pbuffer.data(idx_g_gh + 232);

    auto gs_x_yyyz_xxxxz = pbuffer.data(idx_g_gh + 233);

    auto gs_x_yyyz_xxxyy = pbuffer.data(idx_g_gh + 234);

    auto gs_x_yyyz_xxxyz = pbuffer.data(idx_g_gh + 235);

    auto gs_x_yyyz_xxxzz = pbuffer.data(idx_g_gh + 236);

    auto gs_x_yyyz_xxyyy = pbuffer.data(idx_g_gh + 237);

    auto gs_x_yyyz_xxyyz = pbuffer.data(idx_g_gh + 238);

    auto gs_x_yyyz_xxyzz = pbuffer.data(idx_g_gh + 239);

    auto gs_x_yyyz_xxzzz = pbuffer.data(idx_g_gh + 240);

    auto gs_x_yyyz_xyyyy = pbuffer.data(idx_g_gh + 241);

    auto gs_x_yyyz_xyyyz = pbuffer.data(idx_g_gh + 242);

    auto gs_x_yyyz_xyyzz = pbuffer.data(idx_g_gh + 243);

    auto gs_x_yyyz_xyzzz = pbuffer.data(idx_g_gh + 244);

    auto gs_x_yyyz_xzzzz = pbuffer.data(idx_g_gh + 245);

    auto gs_x_yyyz_yyyyy = pbuffer.data(idx_g_gh + 246);

    auto gs_x_yyyz_yyyyz = pbuffer.data(idx_g_gh + 247);

    auto gs_x_yyyz_yyyzz = pbuffer.data(idx_g_gh + 248);

    auto gs_x_yyyz_yyzzz = pbuffer.data(idx_g_gh + 249);

    auto gs_x_yyyz_yzzzz = pbuffer.data(idx_g_gh + 250);

    auto gs_x_yyyz_zzzzz = pbuffer.data(idx_g_gh + 251);

    #pragma omp simd aligned(gc_x, gs_x_yyyz_xxxxx, gs_x_yyyz_xxxxy, gs_x_yyyz_xxxxz, gs_x_yyyz_xxxyy, gs_x_yyyz_xxxyz, gs_x_yyyz_xxxzz, gs_x_yyyz_xxyyy, gs_x_yyyz_xxyyz, gs_x_yyyz_xxyzz, gs_x_yyyz_xxzzz, gs_x_yyyz_xyyyy, gs_x_yyyz_xyyyz, gs_x_yyyz_xyyzz, gs_x_yyyz_xyzzz, gs_x_yyyz_xzzzz, gs_x_yyyz_yyyyy, gs_x_yyyz_yyyyz, gs_x_yyyz_yyyzz, gs_x_yyyz_yyzzz, gs_x_yyyz_yzzzz, gs_x_yyyz_zzzzz, ts_yyyz_xxxx, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxy, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxz, ts_yyyz_xxxzz, ts_yyyz_xxyy, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyz, ts_yyyz_xxyzz, ts_yyyz_xxzz, ts_yyyz_xxzzz, ts_yyyz_xyyy, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyz, ts_yyyz_xyyzz, ts_yyyz_xyzz, ts_yyyz_xyzzz, ts_yyyz_xzzz, ts_yyyz_xzzzz, ts_yyyz_yyyy, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyz, ts_yyyz_yyyzz, ts_yyyz_yyzz, ts_yyyz_yyzzz, ts_yyyz_yzzz, ts_yyyz_yzzzz, ts_yyyz_zzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyz_xxxxx[i] = 10.0 * ts_yyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxxy[i] = 8.0 * ts_yyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxxz[i] = 8.0 * ts_yyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxyy[i] = 6.0 * ts_yyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxyz[i] = 6.0 * ts_yyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxzz[i] = 6.0 * ts_yyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxyyy[i] = 4.0 * ts_yyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxyyz[i] = 4.0 * ts_yyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxyzz[i] = 4.0 * ts_yyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxzzz[i] = 4.0 * ts_yyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyyyy[i] = 2.0 * ts_yyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyyyz[i] = 2.0 * ts_yyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyyzz[i] = 2.0 * ts_yyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyzzz[i] = 2.0 * ts_yyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xzzzz[i] = 2.0 * ts_yyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyyyy[i] = 2.0 * ts_yyyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyyyz[i] = 2.0 * ts_yyyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyyzz[i] = 2.0 * ts_yyyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyzzz[i] = 2.0 * ts_yyyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yzzzz[i] = 2.0 * ts_yyyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_zzzzz[i] = 2.0 * ts_yyyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 252-273 components of targeted buffer : GH

    auto gs_x_yyzz_xxxxx = pbuffer.data(idx_g_gh + 252);

    auto gs_x_yyzz_xxxxy = pbuffer.data(idx_g_gh + 253);

    auto gs_x_yyzz_xxxxz = pbuffer.data(idx_g_gh + 254);

    auto gs_x_yyzz_xxxyy = pbuffer.data(idx_g_gh + 255);

    auto gs_x_yyzz_xxxyz = pbuffer.data(idx_g_gh + 256);

    auto gs_x_yyzz_xxxzz = pbuffer.data(idx_g_gh + 257);

    auto gs_x_yyzz_xxyyy = pbuffer.data(idx_g_gh + 258);

    auto gs_x_yyzz_xxyyz = pbuffer.data(idx_g_gh + 259);

    auto gs_x_yyzz_xxyzz = pbuffer.data(idx_g_gh + 260);

    auto gs_x_yyzz_xxzzz = pbuffer.data(idx_g_gh + 261);

    auto gs_x_yyzz_xyyyy = pbuffer.data(idx_g_gh + 262);

    auto gs_x_yyzz_xyyyz = pbuffer.data(idx_g_gh + 263);

    auto gs_x_yyzz_xyyzz = pbuffer.data(idx_g_gh + 264);

    auto gs_x_yyzz_xyzzz = pbuffer.data(idx_g_gh + 265);

    auto gs_x_yyzz_xzzzz = pbuffer.data(idx_g_gh + 266);

    auto gs_x_yyzz_yyyyy = pbuffer.data(idx_g_gh + 267);

    auto gs_x_yyzz_yyyyz = pbuffer.data(idx_g_gh + 268);

    auto gs_x_yyzz_yyyzz = pbuffer.data(idx_g_gh + 269);

    auto gs_x_yyzz_yyzzz = pbuffer.data(idx_g_gh + 270);

    auto gs_x_yyzz_yzzzz = pbuffer.data(idx_g_gh + 271);

    auto gs_x_yyzz_zzzzz = pbuffer.data(idx_g_gh + 272);

    #pragma omp simd aligned(gc_x, gs_x_yyzz_xxxxx, gs_x_yyzz_xxxxy, gs_x_yyzz_xxxxz, gs_x_yyzz_xxxyy, gs_x_yyzz_xxxyz, gs_x_yyzz_xxxzz, gs_x_yyzz_xxyyy, gs_x_yyzz_xxyyz, gs_x_yyzz_xxyzz, gs_x_yyzz_xxzzz, gs_x_yyzz_xyyyy, gs_x_yyzz_xyyyz, gs_x_yyzz_xyyzz, gs_x_yyzz_xyzzz, gs_x_yyzz_xzzzz, gs_x_yyzz_yyyyy, gs_x_yyzz_yyyyz, gs_x_yyzz_yyyzz, gs_x_yyzz_yyzzz, gs_x_yyzz_yzzzz, gs_x_yyzz_zzzzz, ts_yyzz_xxxx, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxy, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxz, ts_yyzz_xxxzz, ts_yyzz_xxyy, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyz, ts_yyzz_xxyzz, ts_yyzz_xxzz, ts_yyzz_xxzzz, ts_yyzz_xyyy, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyz, ts_yyzz_xyyzz, ts_yyzz_xyzz, ts_yyzz_xyzzz, ts_yyzz_xzzz, ts_yyzz_xzzzz, ts_yyzz_yyyy, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyz, ts_yyzz_yyyzz, ts_yyzz_yyzz, ts_yyzz_yyzzz, ts_yyzz_yzzz, ts_yyzz_yzzzz, ts_yyzz_zzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzz_xxxxx[i] = 10.0 * ts_yyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxxy[i] = 8.0 * ts_yyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxxz[i] = 8.0 * ts_yyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxyy[i] = 6.0 * ts_yyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxyz[i] = 6.0 * ts_yyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxzz[i] = 6.0 * ts_yyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxyyy[i] = 4.0 * ts_yyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxyyz[i] = 4.0 * ts_yyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxyzz[i] = 4.0 * ts_yyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxzzz[i] = 4.0 * ts_yyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyyyy[i] = 2.0 * ts_yyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyyyz[i] = 2.0 * ts_yyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyyzz[i] = 2.0 * ts_yyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyzzz[i] = 2.0 * ts_yyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xzzzz[i] = 2.0 * ts_yyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyyyy[i] = 2.0 * ts_yyzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyyyz[i] = 2.0 * ts_yyzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyyzz[i] = 2.0 * ts_yyzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyzzz[i] = 2.0 * ts_yyzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yzzzz[i] = 2.0 * ts_yyzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_zzzzz[i] = 2.0 * ts_yyzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 273-294 components of targeted buffer : GH

    auto gs_x_yzzz_xxxxx = pbuffer.data(idx_g_gh + 273);

    auto gs_x_yzzz_xxxxy = pbuffer.data(idx_g_gh + 274);

    auto gs_x_yzzz_xxxxz = pbuffer.data(idx_g_gh + 275);

    auto gs_x_yzzz_xxxyy = pbuffer.data(idx_g_gh + 276);

    auto gs_x_yzzz_xxxyz = pbuffer.data(idx_g_gh + 277);

    auto gs_x_yzzz_xxxzz = pbuffer.data(idx_g_gh + 278);

    auto gs_x_yzzz_xxyyy = pbuffer.data(idx_g_gh + 279);

    auto gs_x_yzzz_xxyyz = pbuffer.data(idx_g_gh + 280);

    auto gs_x_yzzz_xxyzz = pbuffer.data(idx_g_gh + 281);

    auto gs_x_yzzz_xxzzz = pbuffer.data(idx_g_gh + 282);

    auto gs_x_yzzz_xyyyy = pbuffer.data(idx_g_gh + 283);

    auto gs_x_yzzz_xyyyz = pbuffer.data(idx_g_gh + 284);

    auto gs_x_yzzz_xyyzz = pbuffer.data(idx_g_gh + 285);

    auto gs_x_yzzz_xyzzz = pbuffer.data(idx_g_gh + 286);

    auto gs_x_yzzz_xzzzz = pbuffer.data(idx_g_gh + 287);

    auto gs_x_yzzz_yyyyy = pbuffer.data(idx_g_gh + 288);

    auto gs_x_yzzz_yyyyz = pbuffer.data(idx_g_gh + 289);

    auto gs_x_yzzz_yyyzz = pbuffer.data(idx_g_gh + 290);

    auto gs_x_yzzz_yyzzz = pbuffer.data(idx_g_gh + 291);

    auto gs_x_yzzz_yzzzz = pbuffer.data(idx_g_gh + 292);

    auto gs_x_yzzz_zzzzz = pbuffer.data(idx_g_gh + 293);

    #pragma omp simd aligned(gc_x, gs_x_yzzz_xxxxx, gs_x_yzzz_xxxxy, gs_x_yzzz_xxxxz, gs_x_yzzz_xxxyy, gs_x_yzzz_xxxyz, gs_x_yzzz_xxxzz, gs_x_yzzz_xxyyy, gs_x_yzzz_xxyyz, gs_x_yzzz_xxyzz, gs_x_yzzz_xxzzz, gs_x_yzzz_xyyyy, gs_x_yzzz_xyyyz, gs_x_yzzz_xyyzz, gs_x_yzzz_xyzzz, gs_x_yzzz_xzzzz, gs_x_yzzz_yyyyy, gs_x_yzzz_yyyyz, gs_x_yzzz_yyyzz, gs_x_yzzz_yyzzz, gs_x_yzzz_yzzzz, gs_x_yzzz_zzzzz, ts_yzzz_xxxx, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxy, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxz, ts_yzzz_xxxzz, ts_yzzz_xxyy, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyz, ts_yzzz_xxyzz, ts_yzzz_xxzz, ts_yzzz_xxzzz, ts_yzzz_xyyy, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyz, ts_yzzz_xyyzz, ts_yzzz_xyzz, ts_yzzz_xyzzz, ts_yzzz_xzzz, ts_yzzz_xzzzz, ts_yzzz_yyyy, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyz, ts_yzzz_yyyzz, ts_yzzz_yyzz, ts_yzzz_yyzzz, ts_yzzz_yzzz, ts_yzzz_yzzzz, ts_yzzz_zzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzz_xxxxx[i] = 10.0 * ts_yzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxxy[i] = 8.0 * ts_yzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxxz[i] = 8.0 * ts_yzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxyy[i] = 6.0 * ts_yzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxyz[i] = 6.0 * ts_yzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxzz[i] = 6.0 * ts_yzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxyyy[i] = 4.0 * ts_yzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxyyz[i] = 4.0 * ts_yzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxyzz[i] = 4.0 * ts_yzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxzzz[i] = 4.0 * ts_yzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyyyy[i] = 2.0 * ts_yzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyyyz[i] = 2.0 * ts_yzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyyzz[i] = 2.0 * ts_yzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyzzz[i] = 2.0 * ts_yzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xzzzz[i] = 2.0 * ts_yzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyyyy[i] = 2.0 * ts_yzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyyyz[i] = 2.0 * ts_yzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyyzz[i] = 2.0 * ts_yzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyzzz[i] = 2.0 * ts_yzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yzzzz[i] = 2.0 * ts_yzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_zzzzz[i] = 2.0 * ts_yzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 294-315 components of targeted buffer : GH

    auto gs_x_zzzz_xxxxx = pbuffer.data(idx_g_gh + 294);

    auto gs_x_zzzz_xxxxy = pbuffer.data(idx_g_gh + 295);

    auto gs_x_zzzz_xxxxz = pbuffer.data(idx_g_gh + 296);

    auto gs_x_zzzz_xxxyy = pbuffer.data(idx_g_gh + 297);

    auto gs_x_zzzz_xxxyz = pbuffer.data(idx_g_gh + 298);

    auto gs_x_zzzz_xxxzz = pbuffer.data(idx_g_gh + 299);

    auto gs_x_zzzz_xxyyy = pbuffer.data(idx_g_gh + 300);

    auto gs_x_zzzz_xxyyz = pbuffer.data(idx_g_gh + 301);

    auto gs_x_zzzz_xxyzz = pbuffer.data(idx_g_gh + 302);

    auto gs_x_zzzz_xxzzz = pbuffer.data(idx_g_gh + 303);

    auto gs_x_zzzz_xyyyy = pbuffer.data(idx_g_gh + 304);

    auto gs_x_zzzz_xyyyz = pbuffer.data(idx_g_gh + 305);

    auto gs_x_zzzz_xyyzz = pbuffer.data(idx_g_gh + 306);

    auto gs_x_zzzz_xyzzz = pbuffer.data(idx_g_gh + 307);

    auto gs_x_zzzz_xzzzz = pbuffer.data(idx_g_gh + 308);

    auto gs_x_zzzz_yyyyy = pbuffer.data(idx_g_gh + 309);

    auto gs_x_zzzz_yyyyz = pbuffer.data(idx_g_gh + 310);

    auto gs_x_zzzz_yyyzz = pbuffer.data(idx_g_gh + 311);

    auto gs_x_zzzz_yyzzz = pbuffer.data(idx_g_gh + 312);

    auto gs_x_zzzz_yzzzz = pbuffer.data(idx_g_gh + 313);

    auto gs_x_zzzz_zzzzz = pbuffer.data(idx_g_gh + 314);

    #pragma omp simd aligned(gc_x, gs_x_zzzz_xxxxx, gs_x_zzzz_xxxxy, gs_x_zzzz_xxxxz, gs_x_zzzz_xxxyy, gs_x_zzzz_xxxyz, gs_x_zzzz_xxxzz, gs_x_zzzz_xxyyy, gs_x_zzzz_xxyyz, gs_x_zzzz_xxyzz, gs_x_zzzz_xxzzz, gs_x_zzzz_xyyyy, gs_x_zzzz_xyyyz, gs_x_zzzz_xyyzz, gs_x_zzzz_xyzzz, gs_x_zzzz_xzzzz, gs_x_zzzz_yyyyy, gs_x_zzzz_yyyyz, gs_x_zzzz_yyyzz, gs_x_zzzz_yyzzz, gs_x_zzzz_yzzzz, gs_x_zzzz_zzzzz, ts_zzzz_xxxx, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxy, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxz, ts_zzzz_xxxzz, ts_zzzz_xxyy, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyz, ts_zzzz_xxyzz, ts_zzzz_xxzz, ts_zzzz_xxzzz, ts_zzzz_xyyy, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyz, ts_zzzz_xyyzz, ts_zzzz_xyzz, ts_zzzz_xyzzz, ts_zzzz_xzzz, ts_zzzz_xzzzz, ts_zzzz_yyyy, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyz, ts_zzzz_yyyzz, ts_zzzz_yyzz, ts_zzzz_yyzzz, ts_zzzz_yzzz, ts_zzzz_yzzzz, ts_zzzz_zzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzz_xxxxx[i] = 10.0 * ts_zzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxxy[i] = 8.0 * ts_zzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxxz[i] = 8.0 * ts_zzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxyy[i] = 6.0 * ts_zzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxyz[i] = 6.0 * ts_zzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxzz[i] = 6.0 * ts_zzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxyyy[i] = 4.0 * ts_zzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxyyz[i] = 4.0 * ts_zzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxyzz[i] = 4.0 * ts_zzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxzzz[i] = 4.0 * ts_zzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyyyy[i] = 2.0 * ts_zzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyyyz[i] = 2.0 * ts_zzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyyzz[i] = 2.0 * ts_zzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyzzz[i] = 2.0 * ts_zzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xzzzz[i] = 2.0 * ts_zzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyyyy[i] = 2.0 * ts_zzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyyyz[i] = 2.0 * ts_zzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyyzz[i] = 2.0 * ts_zzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyzzz[i] = 2.0 * ts_zzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yzzzz[i] = 2.0 * ts_zzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_zzzzz[i] = 2.0 * ts_zzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 315-336 components of targeted buffer : GH

    auto gs_y_xxxx_xxxxx = pbuffer.data(idx_g_gh + 315);

    auto gs_y_xxxx_xxxxy = pbuffer.data(idx_g_gh + 316);

    auto gs_y_xxxx_xxxxz = pbuffer.data(idx_g_gh + 317);

    auto gs_y_xxxx_xxxyy = pbuffer.data(idx_g_gh + 318);

    auto gs_y_xxxx_xxxyz = pbuffer.data(idx_g_gh + 319);

    auto gs_y_xxxx_xxxzz = pbuffer.data(idx_g_gh + 320);

    auto gs_y_xxxx_xxyyy = pbuffer.data(idx_g_gh + 321);

    auto gs_y_xxxx_xxyyz = pbuffer.data(idx_g_gh + 322);

    auto gs_y_xxxx_xxyzz = pbuffer.data(idx_g_gh + 323);

    auto gs_y_xxxx_xxzzz = pbuffer.data(idx_g_gh + 324);

    auto gs_y_xxxx_xyyyy = pbuffer.data(idx_g_gh + 325);

    auto gs_y_xxxx_xyyyz = pbuffer.data(idx_g_gh + 326);

    auto gs_y_xxxx_xyyzz = pbuffer.data(idx_g_gh + 327);

    auto gs_y_xxxx_xyzzz = pbuffer.data(idx_g_gh + 328);

    auto gs_y_xxxx_xzzzz = pbuffer.data(idx_g_gh + 329);

    auto gs_y_xxxx_yyyyy = pbuffer.data(idx_g_gh + 330);

    auto gs_y_xxxx_yyyyz = pbuffer.data(idx_g_gh + 331);

    auto gs_y_xxxx_yyyzz = pbuffer.data(idx_g_gh + 332);

    auto gs_y_xxxx_yyzzz = pbuffer.data(idx_g_gh + 333);

    auto gs_y_xxxx_yzzzz = pbuffer.data(idx_g_gh + 334);

    auto gs_y_xxxx_zzzzz = pbuffer.data(idx_g_gh + 335);

    #pragma omp simd aligned(gc_y, gs_y_xxxx_xxxxx, gs_y_xxxx_xxxxy, gs_y_xxxx_xxxxz, gs_y_xxxx_xxxyy, gs_y_xxxx_xxxyz, gs_y_xxxx_xxxzz, gs_y_xxxx_xxyyy, gs_y_xxxx_xxyyz, gs_y_xxxx_xxyzz, gs_y_xxxx_xxzzz, gs_y_xxxx_xyyyy, gs_y_xxxx_xyyyz, gs_y_xxxx_xyyzz, gs_y_xxxx_xyzzz, gs_y_xxxx_xzzzz, gs_y_xxxx_yyyyy, gs_y_xxxx_yyyyz, gs_y_xxxx_yyyzz, gs_y_xxxx_yyzzz, gs_y_xxxx_yzzzz, gs_y_xxxx_zzzzz, ts_xxxx_xxxx, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxy, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxz, ts_xxxx_xxxzz, ts_xxxx_xxyy, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyz, ts_xxxx_xxyzz, ts_xxxx_xxzz, ts_xxxx_xxzzz, ts_xxxx_xyyy, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyz, ts_xxxx_xyyzz, ts_xxxx_xyzz, ts_xxxx_xyzzz, ts_xxxx_xzzz, ts_xxxx_xzzzz, ts_xxxx_yyyy, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyz, ts_xxxx_yyyzz, ts_xxxx_yyzz, ts_xxxx_yyzzz, ts_xxxx_yzzz, ts_xxxx_yzzzz, ts_xxxx_zzzz, ts_xxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxx_xxxxx[i] = 2.0 * ts_xxxx_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxxy[i] = 2.0 * ts_xxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxxz[i] = 2.0 * ts_xxxx_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxyy[i] = 4.0 * ts_xxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxyz[i] = 2.0 * ts_xxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxzz[i] = 2.0 * ts_xxxx_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxyyy[i] = 6.0 * ts_xxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxyyz[i] = 4.0 * ts_xxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxyzz[i] = 2.0 * ts_xxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxzzz[i] = 2.0 * ts_xxxx_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyyyy[i] = 8.0 * ts_xxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyyyz[i] = 6.0 * ts_xxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyyzz[i] = 4.0 * ts_xxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyzzz[i] = 2.0 * ts_xxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xzzzz[i] = 2.0 * ts_xxxx_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyyyy[i] = 10.0 * ts_xxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyyyz[i] = 8.0 * ts_xxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyyzz[i] = 6.0 * ts_xxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyzzz[i] = 4.0 * ts_xxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yzzzz[i] = 2.0 * ts_xxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_zzzzz[i] = 2.0 * ts_xxxx_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 336-357 components of targeted buffer : GH

    auto gs_y_xxxy_xxxxx = pbuffer.data(idx_g_gh + 336);

    auto gs_y_xxxy_xxxxy = pbuffer.data(idx_g_gh + 337);

    auto gs_y_xxxy_xxxxz = pbuffer.data(idx_g_gh + 338);

    auto gs_y_xxxy_xxxyy = pbuffer.data(idx_g_gh + 339);

    auto gs_y_xxxy_xxxyz = pbuffer.data(idx_g_gh + 340);

    auto gs_y_xxxy_xxxzz = pbuffer.data(idx_g_gh + 341);

    auto gs_y_xxxy_xxyyy = pbuffer.data(idx_g_gh + 342);

    auto gs_y_xxxy_xxyyz = pbuffer.data(idx_g_gh + 343);

    auto gs_y_xxxy_xxyzz = pbuffer.data(idx_g_gh + 344);

    auto gs_y_xxxy_xxzzz = pbuffer.data(idx_g_gh + 345);

    auto gs_y_xxxy_xyyyy = pbuffer.data(idx_g_gh + 346);

    auto gs_y_xxxy_xyyyz = pbuffer.data(idx_g_gh + 347);

    auto gs_y_xxxy_xyyzz = pbuffer.data(idx_g_gh + 348);

    auto gs_y_xxxy_xyzzz = pbuffer.data(idx_g_gh + 349);

    auto gs_y_xxxy_xzzzz = pbuffer.data(idx_g_gh + 350);

    auto gs_y_xxxy_yyyyy = pbuffer.data(idx_g_gh + 351);

    auto gs_y_xxxy_yyyyz = pbuffer.data(idx_g_gh + 352);

    auto gs_y_xxxy_yyyzz = pbuffer.data(idx_g_gh + 353);

    auto gs_y_xxxy_yyzzz = pbuffer.data(idx_g_gh + 354);

    auto gs_y_xxxy_yzzzz = pbuffer.data(idx_g_gh + 355);

    auto gs_y_xxxy_zzzzz = pbuffer.data(idx_g_gh + 356);

    #pragma omp simd aligned(gc_y, gs_y_xxxy_xxxxx, gs_y_xxxy_xxxxy, gs_y_xxxy_xxxxz, gs_y_xxxy_xxxyy, gs_y_xxxy_xxxyz, gs_y_xxxy_xxxzz, gs_y_xxxy_xxyyy, gs_y_xxxy_xxyyz, gs_y_xxxy_xxyzz, gs_y_xxxy_xxzzz, gs_y_xxxy_xyyyy, gs_y_xxxy_xyyyz, gs_y_xxxy_xyyzz, gs_y_xxxy_xyzzz, gs_y_xxxy_xzzzz, gs_y_xxxy_yyyyy, gs_y_xxxy_yyyyz, gs_y_xxxy_yyyzz, gs_y_xxxy_yyzzz, gs_y_xxxy_yzzzz, gs_y_xxxy_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxzz, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyzz, ts_xxx_xxzzz, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyzz, ts_xxx_xyzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyzz, ts_xxx_yyzzz, ts_xxx_yzzzz, ts_xxx_zzzzz, ts_xxxy_xxxx, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxy, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxz, ts_xxxy_xxxzz, ts_xxxy_xxyy, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyz, ts_xxxy_xxyzz, ts_xxxy_xxzz, ts_xxxy_xxzzz, ts_xxxy_xyyy, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyz, ts_xxxy_xyyzz, ts_xxxy_xyzz, ts_xxxy_xyzzz, ts_xxxy_xzzz, ts_xxxy_xzzzz, ts_xxxy_yyyy, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyz, ts_xxxy_yyyzz, ts_xxxy_yyzz, ts_xxxy_yyzzz, ts_xxxy_yzzz, ts_xxxy_yzzzz, ts_xxxy_zzzz, ts_xxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxy_xxxxx[i] = 2.0 * ts_xxx_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxxy[i] = 2.0 * ts_xxx_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxxz[i] = 2.0 * ts_xxx_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxyy[i] = 2.0 * ts_xxx_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxyz[i] = 2.0 * ts_xxx_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxzz[i] = 2.0 * ts_xxx_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxyyy[i] = 2.0 * ts_xxx_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxyyz[i] = 2.0 * ts_xxx_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxyzz[i] = 2.0 * ts_xxx_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxzzz[i] = 2.0 * ts_xxx_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyyyy[i] = 2.0 * ts_xxx_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyyyz[i] = 2.0 * ts_xxx_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyyzz[i] = 2.0 * ts_xxx_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyzzz[i] = 2.0 * ts_xxx_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xzzzz[i] = 2.0 * ts_xxx_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyyyy[i] = 2.0 * ts_xxx_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyyyz[i] = 2.0 * ts_xxx_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyyzz[i] = 2.0 * ts_xxx_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyzzz[i] = 2.0 * ts_xxx_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yzzzz[i] = 2.0 * ts_xxx_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_zzzzz[i] = 2.0 * ts_xxx_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 357-378 components of targeted buffer : GH

    auto gs_y_xxxz_xxxxx = pbuffer.data(idx_g_gh + 357);

    auto gs_y_xxxz_xxxxy = pbuffer.data(idx_g_gh + 358);

    auto gs_y_xxxz_xxxxz = pbuffer.data(idx_g_gh + 359);

    auto gs_y_xxxz_xxxyy = pbuffer.data(idx_g_gh + 360);

    auto gs_y_xxxz_xxxyz = pbuffer.data(idx_g_gh + 361);

    auto gs_y_xxxz_xxxzz = pbuffer.data(idx_g_gh + 362);

    auto gs_y_xxxz_xxyyy = pbuffer.data(idx_g_gh + 363);

    auto gs_y_xxxz_xxyyz = pbuffer.data(idx_g_gh + 364);

    auto gs_y_xxxz_xxyzz = pbuffer.data(idx_g_gh + 365);

    auto gs_y_xxxz_xxzzz = pbuffer.data(idx_g_gh + 366);

    auto gs_y_xxxz_xyyyy = pbuffer.data(idx_g_gh + 367);

    auto gs_y_xxxz_xyyyz = pbuffer.data(idx_g_gh + 368);

    auto gs_y_xxxz_xyyzz = pbuffer.data(idx_g_gh + 369);

    auto gs_y_xxxz_xyzzz = pbuffer.data(idx_g_gh + 370);

    auto gs_y_xxxz_xzzzz = pbuffer.data(idx_g_gh + 371);

    auto gs_y_xxxz_yyyyy = pbuffer.data(idx_g_gh + 372);

    auto gs_y_xxxz_yyyyz = pbuffer.data(idx_g_gh + 373);

    auto gs_y_xxxz_yyyzz = pbuffer.data(idx_g_gh + 374);

    auto gs_y_xxxz_yyzzz = pbuffer.data(idx_g_gh + 375);

    auto gs_y_xxxz_yzzzz = pbuffer.data(idx_g_gh + 376);

    auto gs_y_xxxz_zzzzz = pbuffer.data(idx_g_gh + 377);

    #pragma omp simd aligned(gc_y, gs_y_xxxz_xxxxx, gs_y_xxxz_xxxxy, gs_y_xxxz_xxxxz, gs_y_xxxz_xxxyy, gs_y_xxxz_xxxyz, gs_y_xxxz_xxxzz, gs_y_xxxz_xxyyy, gs_y_xxxz_xxyyz, gs_y_xxxz_xxyzz, gs_y_xxxz_xxzzz, gs_y_xxxz_xyyyy, gs_y_xxxz_xyyyz, gs_y_xxxz_xyyzz, gs_y_xxxz_xyzzz, gs_y_xxxz_xzzzz, gs_y_xxxz_yyyyy, gs_y_xxxz_yyyyz, gs_y_xxxz_yyyzz, gs_y_xxxz_yyzzz, gs_y_xxxz_yzzzz, gs_y_xxxz_zzzzz, ts_xxxz_xxxx, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxy, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxz, ts_xxxz_xxxzz, ts_xxxz_xxyy, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyz, ts_xxxz_xxyzz, ts_xxxz_xxzz, ts_xxxz_xxzzz, ts_xxxz_xyyy, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyz, ts_xxxz_xyyzz, ts_xxxz_xyzz, ts_xxxz_xyzzz, ts_xxxz_xzzz, ts_xxxz_xzzzz, ts_xxxz_yyyy, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyz, ts_xxxz_yyyzz, ts_xxxz_yyzz, ts_xxxz_yyzzz, ts_xxxz_yzzz, ts_xxxz_yzzzz, ts_xxxz_zzzz, ts_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxz_xxxxx[i] = 2.0 * ts_xxxz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxxy[i] = 2.0 * ts_xxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxxz[i] = 2.0 * ts_xxxz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxyy[i] = 4.0 * ts_xxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxyz[i] = 2.0 * ts_xxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxzz[i] = 2.0 * ts_xxxz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxyyy[i] = 6.0 * ts_xxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxyyz[i] = 4.0 * ts_xxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxyzz[i] = 2.0 * ts_xxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxzzz[i] = 2.0 * ts_xxxz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyyyy[i] = 8.0 * ts_xxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyyyz[i] = 6.0 * ts_xxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyyzz[i] = 4.0 * ts_xxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyzzz[i] = 2.0 * ts_xxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xzzzz[i] = 2.0 * ts_xxxz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyyyy[i] = 10.0 * ts_xxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyyyz[i] = 8.0 * ts_xxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyyzz[i] = 6.0 * ts_xxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyzzz[i] = 4.0 * ts_xxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yzzzz[i] = 2.0 * ts_xxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_zzzzz[i] = 2.0 * ts_xxxz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 378-399 components of targeted buffer : GH

    auto gs_y_xxyy_xxxxx = pbuffer.data(idx_g_gh + 378);

    auto gs_y_xxyy_xxxxy = pbuffer.data(idx_g_gh + 379);

    auto gs_y_xxyy_xxxxz = pbuffer.data(idx_g_gh + 380);

    auto gs_y_xxyy_xxxyy = pbuffer.data(idx_g_gh + 381);

    auto gs_y_xxyy_xxxyz = pbuffer.data(idx_g_gh + 382);

    auto gs_y_xxyy_xxxzz = pbuffer.data(idx_g_gh + 383);

    auto gs_y_xxyy_xxyyy = pbuffer.data(idx_g_gh + 384);

    auto gs_y_xxyy_xxyyz = pbuffer.data(idx_g_gh + 385);

    auto gs_y_xxyy_xxyzz = pbuffer.data(idx_g_gh + 386);

    auto gs_y_xxyy_xxzzz = pbuffer.data(idx_g_gh + 387);

    auto gs_y_xxyy_xyyyy = pbuffer.data(idx_g_gh + 388);

    auto gs_y_xxyy_xyyyz = pbuffer.data(idx_g_gh + 389);

    auto gs_y_xxyy_xyyzz = pbuffer.data(idx_g_gh + 390);

    auto gs_y_xxyy_xyzzz = pbuffer.data(idx_g_gh + 391);

    auto gs_y_xxyy_xzzzz = pbuffer.data(idx_g_gh + 392);

    auto gs_y_xxyy_yyyyy = pbuffer.data(idx_g_gh + 393);

    auto gs_y_xxyy_yyyyz = pbuffer.data(idx_g_gh + 394);

    auto gs_y_xxyy_yyyzz = pbuffer.data(idx_g_gh + 395);

    auto gs_y_xxyy_yyzzz = pbuffer.data(idx_g_gh + 396);

    auto gs_y_xxyy_yzzzz = pbuffer.data(idx_g_gh + 397);

    auto gs_y_xxyy_zzzzz = pbuffer.data(idx_g_gh + 398);

    #pragma omp simd aligned(gc_y, gs_y_xxyy_xxxxx, gs_y_xxyy_xxxxy, gs_y_xxyy_xxxxz, gs_y_xxyy_xxxyy, gs_y_xxyy_xxxyz, gs_y_xxyy_xxxzz, gs_y_xxyy_xxyyy, gs_y_xxyy_xxyyz, gs_y_xxyy_xxyzz, gs_y_xxyy_xxzzz, gs_y_xxyy_xyyyy, gs_y_xxyy_xyyyz, gs_y_xxyy_xyyzz, gs_y_xxyy_xyzzz, gs_y_xxyy_xzzzz, gs_y_xxyy_yyyyy, gs_y_xxyy_yyyyz, gs_y_xxyy_yyyzz, gs_y_xxyy_yyzzz, gs_y_xxyy_yzzzz, gs_y_xxyy_zzzzz, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxzz, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyzz, ts_xxy_xxzzz, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyzz, ts_xxy_xyzzz, ts_xxy_xzzzz, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyzz, ts_xxy_yyzzz, ts_xxy_yzzzz, ts_xxy_zzzzz, ts_xxyy_xxxx, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxz, ts_xxyy_xxxzz, ts_xxyy_xxyy, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyz, ts_xxyy_xxyzz, ts_xxyy_xxzz, ts_xxyy_xxzzz, ts_xxyy_xyyy, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyz, ts_xxyy_xyyzz, ts_xxyy_xyzz, ts_xxyy_xyzzz, ts_xxyy_xzzz, ts_xxyy_xzzzz, ts_xxyy_yyyy, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyz, ts_xxyy_yyyzz, ts_xxyy_yyzz, ts_xxyy_yyzzz, ts_xxyy_yzzz, ts_xxyy_yzzzz, ts_xxyy_zzzz, ts_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyy_xxxxx[i] = 4.0 * ts_xxy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxxy[i] = 4.0 * ts_xxy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxxz[i] = 4.0 * ts_xxy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxyy[i] = 4.0 * ts_xxy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxyz[i] = 4.0 * ts_xxy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxzz[i] = 4.0 * ts_xxy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxyyy[i] = 4.0 * ts_xxy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxyyz[i] = 4.0 * ts_xxy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxyzz[i] = 4.0 * ts_xxy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxzzz[i] = 4.0 * ts_xxy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyyyy[i] = 4.0 * ts_xxy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyyyz[i] = 4.0 * ts_xxy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyyzz[i] = 4.0 * ts_xxy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyzzz[i] = 4.0 * ts_xxy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xzzzz[i] = 4.0 * ts_xxy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyyyy[i] = 4.0 * ts_xxy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyyyz[i] = 4.0 * ts_xxy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyyzz[i] = 4.0 * ts_xxy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyzzz[i] = 4.0 * ts_xxy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yzzzz[i] = 4.0 * ts_xxy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_zzzzz[i] = 4.0 * ts_xxy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 399-420 components of targeted buffer : GH

    auto gs_y_xxyz_xxxxx = pbuffer.data(idx_g_gh + 399);

    auto gs_y_xxyz_xxxxy = pbuffer.data(idx_g_gh + 400);

    auto gs_y_xxyz_xxxxz = pbuffer.data(idx_g_gh + 401);

    auto gs_y_xxyz_xxxyy = pbuffer.data(idx_g_gh + 402);

    auto gs_y_xxyz_xxxyz = pbuffer.data(idx_g_gh + 403);

    auto gs_y_xxyz_xxxzz = pbuffer.data(idx_g_gh + 404);

    auto gs_y_xxyz_xxyyy = pbuffer.data(idx_g_gh + 405);

    auto gs_y_xxyz_xxyyz = pbuffer.data(idx_g_gh + 406);

    auto gs_y_xxyz_xxyzz = pbuffer.data(idx_g_gh + 407);

    auto gs_y_xxyz_xxzzz = pbuffer.data(idx_g_gh + 408);

    auto gs_y_xxyz_xyyyy = pbuffer.data(idx_g_gh + 409);

    auto gs_y_xxyz_xyyyz = pbuffer.data(idx_g_gh + 410);

    auto gs_y_xxyz_xyyzz = pbuffer.data(idx_g_gh + 411);

    auto gs_y_xxyz_xyzzz = pbuffer.data(idx_g_gh + 412);

    auto gs_y_xxyz_xzzzz = pbuffer.data(idx_g_gh + 413);

    auto gs_y_xxyz_yyyyy = pbuffer.data(idx_g_gh + 414);

    auto gs_y_xxyz_yyyyz = pbuffer.data(idx_g_gh + 415);

    auto gs_y_xxyz_yyyzz = pbuffer.data(idx_g_gh + 416);

    auto gs_y_xxyz_yyzzz = pbuffer.data(idx_g_gh + 417);

    auto gs_y_xxyz_yzzzz = pbuffer.data(idx_g_gh + 418);

    auto gs_y_xxyz_zzzzz = pbuffer.data(idx_g_gh + 419);

    #pragma omp simd aligned(gc_y, gs_y_xxyz_xxxxx, gs_y_xxyz_xxxxy, gs_y_xxyz_xxxxz, gs_y_xxyz_xxxyy, gs_y_xxyz_xxxyz, gs_y_xxyz_xxxzz, gs_y_xxyz_xxyyy, gs_y_xxyz_xxyyz, gs_y_xxyz_xxyzz, gs_y_xxyz_xxzzz, gs_y_xxyz_xyyyy, gs_y_xxyz_xyyyz, gs_y_xxyz_xyyzz, gs_y_xxyz_xyzzz, gs_y_xxyz_xzzzz, gs_y_xxyz_yyyyy, gs_y_xxyz_yyyyz, gs_y_xxyz_yyyzz, gs_y_xxyz_yyzzz, gs_y_xxyz_yzzzz, gs_y_xxyz_zzzzz, ts_xxyz_xxxx, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxy, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxz, ts_xxyz_xxxzz, ts_xxyz_xxyy, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyz, ts_xxyz_xxyzz, ts_xxyz_xxzz, ts_xxyz_xxzzz, ts_xxyz_xyyy, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyz, ts_xxyz_xyyzz, ts_xxyz_xyzz, ts_xxyz_xyzzz, ts_xxyz_xzzz, ts_xxyz_xzzzz, ts_xxyz_yyyy, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyz, ts_xxyz_yyyzz, ts_xxyz_yyzz, ts_xxyz_yyzzz, ts_xxyz_yzzz, ts_xxyz_yzzzz, ts_xxyz_zzzz, ts_xxyz_zzzzz, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxzz, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyzz, ts_xxz_xxzzz, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyzz, ts_xxz_xyzzz, ts_xxz_xzzzz, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyzz, ts_xxz_yyzzz, ts_xxz_yzzzz, ts_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyz_xxxxx[i] = 2.0 * ts_xxz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxxy[i] = 2.0 * ts_xxz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxxz[i] = 2.0 * ts_xxz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxyy[i] = 2.0 * ts_xxz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxyz[i] = 2.0 * ts_xxz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxzz[i] = 2.0 * ts_xxz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxyyy[i] = 2.0 * ts_xxz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxyyz[i] = 2.0 * ts_xxz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxyzz[i] = 2.0 * ts_xxz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxzzz[i] = 2.0 * ts_xxz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyyyy[i] = 2.0 * ts_xxz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyyyz[i] = 2.0 * ts_xxz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyyzz[i] = 2.0 * ts_xxz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyzzz[i] = 2.0 * ts_xxz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xzzzz[i] = 2.0 * ts_xxz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyyyy[i] = 2.0 * ts_xxz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyyyz[i] = 2.0 * ts_xxz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyyzz[i] = 2.0 * ts_xxz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyzzz[i] = 2.0 * ts_xxz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yzzzz[i] = 2.0 * ts_xxz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_zzzzz[i] = 2.0 * ts_xxz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 420-441 components of targeted buffer : GH

    auto gs_y_xxzz_xxxxx = pbuffer.data(idx_g_gh + 420);

    auto gs_y_xxzz_xxxxy = pbuffer.data(idx_g_gh + 421);

    auto gs_y_xxzz_xxxxz = pbuffer.data(idx_g_gh + 422);

    auto gs_y_xxzz_xxxyy = pbuffer.data(idx_g_gh + 423);

    auto gs_y_xxzz_xxxyz = pbuffer.data(idx_g_gh + 424);

    auto gs_y_xxzz_xxxzz = pbuffer.data(idx_g_gh + 425);

    auto gs_y_xxzz_xxyyy = pbuffer.data(idx_g_gh + 426);

    auto gs_y_xxzz_xxyyz = pbuffer.data(idx_g_gh + 427);

    auto gs_y_xxzz_xxyzz = pbuffer.data(idx_g_gh + 428);

    auto gs_y_xxzz_xxzzz = pbuffer.data(idx_g_gh + 429);

    auto gs_y_xxzz_xyyyy = pbuffer.data(idx_g_gh + 430);

    auto gs_y_xxzz_xyyyz = pbuffer.data(idx_g_gh + 431);

    auto gs_y_xxzz_xyyzz = pbuffer.data(idx_g_gh + 432);

    auto gs_y_xxzz_xyzzz = pbuffer.data(idx_g_gh + 433);

    auto gs_y_xxzz_xzzzz = pbuffer.data(idx_g_gh + 434);

    auto gs_y_xxzz_yyyyy = pbuffer.data(idx_g_gh + 435);

    auto gs_y_xxzz_yyyyz = pbuffer.data(idx_g_gh + 436);

    auto gs_y_xxzz_yyyzz = pbuffer.data(idx_g_gh + 437);

    auto gs_y_xxzz_yyzzz = pbuffer.data(idx_g_gh + 438);

    auto gs_y_xxzz_yzzzz = pbuffer.data(idx_g_gh + 439);

    auto gs_y_xxzz_zzzzz = pbuffer.data(idx_g_gh + 440);

    #pragma omp simd aligned(gc_y, gs_y_xxzz_xxxxx, gs_y_xxzz_xxxxy, gs_y_xxzz_xxxxz, gs_y_xxzz_xxxyy, gs_y_xxzz_xxxyz, gs_y_xxzz_xxxzz, gs_y_xxzz_xxyyy, gs_y_xxzz_xxyyz, gs_y_xxzz_xxyzz, gs_y_xxzz_xxzzz, gs_y_xxzz_xyyyy, gs_y_xxzz_xyyyz, gs_y_xxzz_xyyzz, gs_y_xxzz_xyzzz, gs_y_xxzz_xzzzz, gs_y_xxzz_yyyyy, gs_y_xxzz_yyyyz, gs_y_xxzz_yyyzz, gs_y_xxzz_yyzzz, gs_y_xxzz_yzzzz, gs_y_xxzz_zzzzz, ts_xxzz_xxxx, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxy, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxz, ts_xxzz_xxxzz, ts_xxzz_xxyy, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyz, ts_xxzz_xxyzz, ts_xxzz_xxzz, ts_xxzz_xxzzz, ts_xxzz_xyyy, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyz, ts_xxzz_xyyzz, ts_xxzz_xyzz, ts_xxzz_xyzzz, ts_xxzz_xzzz, ts_xxzz_xzzzz, ts_xxzz_yyyy, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyz, ts_xxzz_yyyzz, ts_xxzz_yyzz, ts_xxzz_yyzzz, ts_xxzz_yzzz, ts_xxzz_yzzzz, ts_xxzz_zzzz, ts_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzz_xxxxx[i] = 2.0 * ts_xxzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxxy[i] = 2.0 * ts_xxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxxz[i] = 2.0 * ts_xxzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxyy[i] = 4.0 * ts_xxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxyz[i] = 2.0 * ts_xxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxzz[i] = 2.0 * ts_xxzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxyyy[i] = 6.0 * ts_xxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxyyz[i] = 4.0 * ts_xxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxyzz[i] = 2.0 * ts_xxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxzzz[i] = 2.0 * ts_xxzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyyyy[i] = 8.0 * ts_xxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyyyz[i] = 6.0 * ts_xxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyyzz[i] = 4.0 * ts_xxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyzzz[i] = 2.0 * ts_xxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xzzzz[i] = 2.0 * ts_xxzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyyyy[i] = 10.0 * ts_xxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyyyz[i] = 8.0 * ts_xxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyyzz[i] = 6.0 * ts_xxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyzzz[i] = 4.0 * ts_xxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yzzzz[i] = 2.0 * ts_xxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_zzzzz[i] = 2.0 * ts_xxzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 441-462 components of targeted buffer : GH

    auto gs_y_xyyy_xxxxx = pbuffer.data(idx_g_gh + 441);

    auto gs_y_xyyy_xxxxy = pbuffer.data(idx_g_gh + 442);

    auto gs_y_xyyy_xxxxz = pbuffer.data(idx_g_gh + 443);

    auto gs_y_xyyy_xxxyy = pbuffer.data(idx_g_gh + 444);

    auto gs_y_xyyy_xxxyz = pbuffer.data(idx_g_gh + 445);

    auto gs_y_xyyy_xxxzz = pbuffer.data(idx_g_gh + 446);

    auto gs_y_xyyy_xxyyy = pbuffer.data(idx_g_gh + 447);

    auto gs_y_xyyy_xxyyz = pbuffer.data(idx_g_gh + 448);

    auto gs_y_xyyy_xxyzz = pbuffer.data(idx_g_gh + 449);

    auto gs_y_xyyy_xxzzz = pbuffer.data(idx_g_gh + 450);

    auto gs_y_xyyy_xyyyy = pbuffer.data(idx_g_gh + 451);

    auto gs_y_xyyy_xyyyz = pbuffer.data(idx_g_gh + 452);

    auto gs_y_xyyy_xyyzz = pbuffer.data(idx_g_gh + 453);

    auto gs_y_xyyy_xyzzz = pbuffer.data(idx_g_gh + 454);

    auto gs_y_xyyy_xzzzz = pbuffer.data(idx_g_gh + 455);

    auto gs_y_xyyy_yyyyy = pbuffer.data(idx_g_gh + 456);

    auto gs_y_xyyy_yyyyz = pbuffer.data(idx_g_gh + 457);

    auto gs_y_xyyy_yyyzz = pbuffer.data(idx_g_gh + 458);

    auto gs_y_xyyy_yyzzz = pbuffer.data(idx_g_gh + 459);

    auto gs_y_xyyy_yzzzz = pbuffer.data(idx_g_gh + 460);

    auto gs_y_xyyy_zzzzz = pbuffer.data(idx_g_gh + 461);

    #pragma omp simd aligned(gc_y, gs_y_xyyy_xxxxx, gs_y_xyyy_xxxxy, gs_y_xyyy_xxxxz, gs_y_xyyy_xxxyy, gs_y_xyyy_xxxyz, gs_y_xyyy_xxxzz, gs_y_xyyy_xxyyy, gs_y_xyyy_xxyyz, gs_y_xyyy_xxyzz, gs_y_xyyy_xxzzz, gs_y_xyyy_xyyyy, gs_y_xyyy_xyyyz, gs_y_xyyy_xyyzz, gs_y_xyyy_xyzzz, gs_y_xyyy_xzzzz, gs_y_xyyy_yyyyy, gs_y_xyyy_yyyyz, gs_y_xyyy_yyyzz, gs_y_xyyy_yyzzz, gs_y_xyyy_yzzzz, gs_y_xyyy_zzzzz, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxzz, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyzz, ts_xyy_xxzzz, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyzz, ts_xyy_xyzzz, ts_xyy_xzzzz, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyzz, ts_xyy_yyzzz, ts_xyy_yzzzz, ts_xyy_zzzzz, ts_xyyy_xxxx, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxy, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxz, ts_xyyy_xxxzz, ts_xyyy_xxyy, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyz, ts_xyyy_xxyzz, ts_xyyy_xxzz, ts_xyyy_xxzzz, ts_xyyy_xyyy, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyz, ts_xyyy_xyyzz, ts_xyyy_xyzz, ts_xyyy_xyzzz, ts_xyyy_xzzz, ts_xyyy_xzzzz, ts_xyyy_yyyy, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyz, ts_xyyy_yyyzz, ts_xyyy_yyzz, ts_xyyy_yyzzz, ts_xyyy_yzzz, ts_xyyy_yzzzz, ts_xyyy_zzzz, ts_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyy_xxxxx[i] = 6.0 * ts_xyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxxy[i] = 6.0 * ts_xyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxxz[i] = 6.0 * ts_xyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxyy[i] = 6.0 * ts_xyy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxyz[i] = 6.0 * ts_xyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxzz[i] = 6.0 * ts_xyy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxyyy[i] = 6.0 * ts_xyy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxyyz[i] = 6.0 * ts_xyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxyzz[i] = 6.0 * ts_xyy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxzzz[i] = 6.0 * ts_xyy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyyyy[i] = 6.0 * ts_xyy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyyyz[i] = 6.0 * ts_xyy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyyzz[i] = 6.0 * ts_xyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyzzz[i] = 6.0 * ts_xyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xzzzz[i] = 6.0 * ts_xyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyyyy[i] = 6.0 * ts_xyy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyyyz[i] = 6.0 * ts_xyy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyyzz[i] = 6.0 * ts_xyy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyzzz[i] = 6.0 * ts_xyy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yzzzz[i] = 6.0 * ts_xyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_zzzzz[i] = 6.0 * ts_xyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 462-483 components of targeted buffer : GH

    auto gs_y_xyyz_xxxxx = pbuffer.data(idx_g_gh + 462);

    auto gs_y_xyyz_xxxxy = pbuffer.data(idx_g_gh + 463);

    auto gs_y_xyyz_xxxxz = pbuffer.data(idx_g_gh + 464);

    auto gs_y_xyyz_xxxyy = pbuffer.data(idx_g_gh + 465);

    auto gs_y_xyyz_xxxyz = pbuffer.data(idx_g_gh + 466);

    auto gs_y_xyyz_xxxzz = pbuffer.data(idx_g_gh + 467);

    auto gs_y_xyyz_xxyyy = pbuffer.data(idx_g_gh + 468);

    auto gs_y_xyyz_xxyyz = pbuffer.data(idx_g_gh + 469);

    auto gs_y_xyyz_xxyzz = pbuffer.data(idx_g_gh + 470);

    auto gs_y_xyyz_xxzzz = pbuffer.data(idx_g_gh + 471);

    auto gs_y_xyyz_xyyyy = pbuffer.data(idx_g_gh + 472);

    auto gs_y_xyyz_xyyyz = pbuffer.data(idx_g_gh + 473);

    auto gs_y_xyyz_xyyzz = pbuffer.data(idx_g_gh + 474);

    auto gs_y_xyyz_xyzzz = pbuffer.data(idx_g_gh + 475);

    auto gs_y_xyyz_xzzzz = pbuffer.data(idx_g_gh + 476);

    auto gs_y_xyyz_yyyyy = pbuffer.data(idx_g_gh + 477);

    auto gs_y_xyyz_yyyyz = pbuffer.data(idx_g_gh + 478);

    auto gs_y_xyyz_yyyzz = pbuffer.data(idx_g_gh + 479);

    auto gs_y_xyyz_yyzzz = pbuffer.data(idx_g_gh + 480);

    auto gs_y_xyyz_yzzzz = pbuffer.data(idx_g_gh + 481);

    auto gs_y_xyyz_zzzzz = pbuffer.data(idx_g_gh + 482);

    #pragma omp simd aligned(gc_y, gs_y_xyyz_xxxxx, gs_y_xyyz_xxxxy, gs_y_xyyz_xxxxz, gs_y_xyyz_xxxyy, gs_y_xyyz_xxxyz, gs_y_xyyz_xxxzz, gs_y_xyyz_xxyyy, gs_y_xyyz_xxyyz, gs_y_xyyz_xxyzz, gs_y_xyyz_xxzzz, gs_y_xyyz_xyyyy, gs_y_xyyz_xyyyz, gs_y_xyyz_xyyzz, gs_y_xyyz_xyzzz, gs_y_xyyz_xzzzz, gs_y_xyyz_yyyyy, gs_y_xyyz_yyyyz, gs_y_xyyz_yyyzz, gs_y_xyyz_yyzzz, gs_y_xyyz_yzzzz, gs_y_xyyz_zzzzz, ts_xyyz_xxxx, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxy, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxz, ts_xyyz_xxxzz, ts_xyyz_xxyy, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyz, ts_xyyz_xxyzz, ts_xyyz_xxzz, ts_xyyz_xxzzz, ts_xyyz_xyyy, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyz, ts_xyyz_xyyzz, ts_xyyz_xyzz, ts_xyyz_xyzzz, ts_xyyz_xzzz, ts_xyyz_xzzzz, ts_xyyz_yyyy, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyz, ts_xyyz_yyyzz, ts_xyyz_yyzz, ts_xyyz_yyzzz, ts_xyyz_yzzz, ts_xyyz_yzzzz, ts_xyyz_zzzz, ts_xyyz_zzzzz, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxzz, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyzz, ts_xyz_xxzzz, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyzz, ts_xyz_xyzzz, ts_xyz_xzzzz, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyzz, ts_xyz_yyzzz, ts_xyz_yzzzz, ts_xyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyz_xxxxx[i] = 4.0 * ts_xyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxxy[i] = 4.0 * ts_xyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxxz[i] = 4.0 * ts_xyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxyy[i] = 4.0 * ts_xyz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxyz[i] = 4.0 * ts_xyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxzz[i] = 4.0 * ts_xyz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxyyy[i] = 4.0 * ts_xyz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxyyz[i] = 4.0 * ts_xyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxyzz[i] = 4.0 * ts_xyz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxzzz[i] = 4.0 * ts_xyz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyyyy[i] = 4.0 * ts_xyz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyyyz[i] = 4.0 * ts_xyz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyyzz[i] = 4.0 * ts_xyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyzzz[i] = 4.0 * ts_xyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xzzzz[i] = 4.0 * ts_xyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyyyy[i] = 4.0 * ts_xyz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyyyz[i] = 4.0 * ts_xyz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyyzz[i] = 4.0 * ts_xyz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyzzz[i] = 4.0 * ts_xyz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yzzzz[i] = 4.0 * ts_xyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_zzzzz[i] = 4.0 * ts_xyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 483-504 components of targeted buffer : GH

    auto gs_y_xyzz_xxxxx = pbuffer.data(idx_g_gh + 483);

    auto gs_y_xyzz_xxxxy = pbuffer.data(idx_g_gh + 484);

    auto gs_y_xyzz_xxxxz = pbuffer.data(idx_g_gh + 485);

    auto gs_y_xyzz_xxxyy = pbuffer.data(idx_g_gh + 486);

    auto gs_y_xyzz_xxxyz = pbuffer.data(idx_g_gh + 487);

    auto gs_y_xyzz_xxxzz = pbuffer.data(idx_g_gh + 488);

    auto gs_y_xyzz_xxyyy = pbuffer.data(idx_g_gh + 489);

    auto gs_y_xyzz_xxyyz = pbuffer.data(idx_g_gh + 490);

    auto gs_y_xyzz_xxyzz = pbuffer.data(idx_g_gh + 491);

    auto gs_y_xyzz_xxzzz = pbuffer.data(idx_g_gh + 492);

    auto gs_y_xyzz_xyyyy = pbuffer.data(idx_g_gh + 493);

    auto gs_y_xyzz_xyyyz = pbuffer.data(idx_g_gh + 494);

    auto gs_y_xyzz_xyyzz = pbuffer.data(idx_g_gh + 495);

    auto gs_y_xyzz_xyzzz = pbuffer.data(idx_g_gh + 496);

    auto gs_y_xyzz_xzzzz = pbuffer.data(idx_g_gh + 497);

    auto gs_y_xyzz_yyyyy = pbuffer.data(idx_g_gh + 498);

    auto gs_y_xyzz_yyyyz = pbuffer.data(idx_g_gh + 499);

    auto gs_y_xyzz_yyyzz = pbuffer.data(idx_g_gh + 500);

    auto gs_y_xyzz_yyzzz = pbuffer.data(idx_g_gh + 501);

    auto gs_y_xyzz_yzzzz = pbuffer.data(idx_g_gh + 502);

    auto gs_y_xyzz_zzzzz = pbuffer.data(idx_g_gh + 503);

    #pragma omp simd aligned(gc_y, gs_y_xyzz_xxxxx, gs_y_xyzz_xxxxy, gs_y_xyzz_xxxxz, gs_y_xyzz_xxxyy, gs_y_xyzz_xxxyz, gs_y_xyzz_xxxzz, gs_y_xyzz_xxyyy, gs_y_xyzz_xxyyz, gs_y_xyzz_xxyzz, gs_y_xyzz_xxzzz, gs_y_xyzz_xyyyy, gs_y_xyzz_xyyyz, gs_y_xyzz_xyyzz, gs_y_xyzz_xyzzz, gs_y_xyzz_xzzzz, gs_y_xyzz_yyyyy, gs_y_xyzz_yyyyz, gs_y_xyzz_yyyzz, gs_y_xyzz_yyzzz, gs_y_xyzz_yzzzz, gs_y_xyzz_zzzzz, ts_xyzz_xxxx, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxy, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxz, ts_xyzz_xxxzz, ts_xyzz_xxyy, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyz, ts_xyzz_xxyzz, ts_xyzz_xxzz, ts_xyzz_xxzzz, ts_xyzz_xyyy, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyz, ts_xyzz_xyyzz, ts_xyzz_xyzz, ts_xyzz_xyzzz, ts_xyzz_xzzz, ts_xyzz_xzzzz, ts_xyzz_yyyy, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyz, ts_xyzz_yyyzz, ts_xyzz_yyzz, ts_xyzz_yyzzz, ts_xyzz_yzzz, ts_xyzz_yzzzz, ts_xyzz_zzzz, ts_xyzz_zzzzz, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxzz, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyzz, ts_xzz_xxzzz, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyzz, ts_xzz_xyzzz, ts_xzz_xzzzz, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyzz, ts_xzz_yyzzz, ts_xzz_yzzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzz_xxxxx[i] = 2.0 * ts_xzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxxy[i] = 2.0 * ts_xzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxxz[i] = 2.0 * ts_xzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxyy[i] = 2.0 * ts_xzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxyz[i] = 2.0 * ts_xzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxzz[i] = 2.0 * ts_xzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxyyy[i] = 2.0 * ts_xzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxyyz[i] = 2.0 * ts_xzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxyzz[i] = 2.0 * ts_xzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxzzz[i] = 2.0 * ts_xzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyyyy[i] = 2.0 * ts_xzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyyyz[i] = 2.0 * ts_xzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyyzz[i] = 2.0 * ts_xzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyzzz[i] = 2.0 * ts_xzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xzzzz[i] = 2.0 * ts_xzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyyyy[i] = 2.0 * ts_xzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyyyz[i] = 2.0 * ts_xzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyyzz[i] = 2.0 * ts_xzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyzzz[i] = 2.0 * ts_xzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yzzzz[i] = 2.0 * ts_xzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_zzzzz[i] = 2.0 * ts_xzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 504-525 components of targeted buffer : GH

    auto gs_y_xzzz_xxxxx = pbuffer.data(idx_g_gh + 504);

    auto gs_y_xzzz_xxxxy = pbuffer.data(idx_g_gh + 505);

    auto gs_y_xzzz_xxxxz = pbuffer.data(idx_g_gh + 506);

    auto gs_y_xzzz_xxxyy = pbuffer.data(idx_g_gh + 507);

    auto gs_y_xzzz_xxxyz = pbuffer.data(idx_g_gh + 508);

    auto gs_y_xzzz_xxxzz = pbuffer.data(idx_g_gh + 509);

    auto gs_y_xzzz_xxyyy = pbuffer.data(idx_g_gh + 510);

    auto gs_y_xzzz_xxyyz = pbuffer.data(idx_g_gh + 511);

    auto gs_y_xzzz_xxyzz = pbuffer.data(idx_g_gh + 512);

    auto gs_y_xzzz_xxzzz = pbuffer.data(idx_g_gh + 513);

    auto gs_y_xzzz_xyyyy = pbuffer.data(idx_g_gh + 514);

    auto gs_y_xzzz_xyyyz = pbuffer.data(idx_g_gh + 515);

    auto gs_y_xzzz_xyyzz = pbuffer.data(idx_g_gh + 516);

    auto gs_y_xzzz_xyzzz = pbuffer.data(idx_g_gh + 517);

    auto gs_y_xzzz_xzzzz = pbuffer.data(idx_g_gh + 518);

    auto gs_y_xzzz_yyyyy = pbuffer.data(idx_g_gh + 519);

    auto gs_y_xzzz_yyyyz = pbuffer.data(idx_g_gh + 520);

    auto gs_y_xzzz_yyyzz = pbuffer.data(idx_g_gh + 521);

    auto gs_y_xzzz_yyzzz = pbuffer.data(idx_g_gh + 522);

    auto gs_y_xzzz_yzzzz = pbuffer.data(idx_g_gh + 523);

    auto gs_y_xzzz_zzzzz = pbuffer.data(idx_g_gh + 524);

    #pragma omp simd aligned(gc_y, gs_y_xzzz_xxxxx, gs_y_xzzz_xxxxy, gs_y_xzzz_xxxxz, gs_y_xzzz_xxxyy, gs_y_xzzz_xxxyz, gs_y_xzzz_xxxzz, gs_y_xzzz_xxyyy, gs_y_xzzz_xxyyz, gs_y_xzzz_xxyzz, gs_y_xzzz_xxzzz, gs_y_xzzz_xyyyy, gs_y_xzzz_xyyyz, gs_y_xzzz_xyyzz, gs_y_xzzz_xyzzz, gs_y_xzzz_xzzzz, gs_y_xzzz_yyyyy, gs_y_xzzz_yyyyz, gs_y_xzzz_yyyzz, gs_y_xzzz_yyzzz, gs_y_xzzz_yzzzz, gs_y_xzzz_zzzzz, ts_xzzz_xxxx, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxy, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxz, ts_xzzz_xxxzz, ts_xzzz_xxyy, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyz, ts_xzzz_xxyzz, ts_xzzz_xxzz, ts_xzzz_xxzzz, ts_xzzz_xyyy, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyz, ts_xzzz_xyyzz, ts_xzzz_xyzz, ts_xzzz_xyzzz, ts_xzzz_xzzz, ts_xzzz_xzzzz, ts_xzzz_yyyy, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyz, ts_xzzz_yyyzz, ts_xzzz_yyzz, ts_xzzz_yyzzz, ts_xzzz_yzzz, ts_xzzz_yzzzz, ts_xzzz_zzzz, ts_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzz_xxxxx[i] = 2.0 * ts_xzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxxy[i] = 2.0 * ts_xzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxxz[i] = 2.0 * ts_xzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxyy[i] = 4.0 * ts_xzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxyz[i] = 2.0 * ts_xzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxzz[i] = 2.0 * ts_xzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxyyy[i] = 6.0 * ts_xzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxyyz[i] = 4.0 * ts_xzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxyzz[i] = 2.0 * ts_xzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxzzz[i] = 2.0 * ts_xzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyyyy[i] = 8.0 * ts_xzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyyyz[i] = 6.0 * ts_xzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyyzz[i] = 4.0 * ts_xzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyzzz[i] = 2.0 * ts_xzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xzzzz[i] = 2.0 * ts_xzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyyyy[i] = 10.0 * ts_xzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyyyz[i] = 8.0 * ts_xzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyyzz[i] = 6.0 * ts_xzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyzzz[i] = 4.0 * ts_xzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yzzzz[i] = 2.0 * ts_xzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_zzzzz[i] = 2.0 * ts_xzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 525-546 components of targeted buffer : GH

    auto gs_y_yyyy_xxxxx = pbuffer.data(idx_g_gh + 525);

    auto gs_y_yyyy_xxxxy = pbuffer.data(idx_g_gh + 526);

    auto gs_y_yyyy_xxxxz = pbuffer.data(idx_g_gh + 527);

    auto gs_y_yyyy_xxxyy = pbuffer.data(idx_g_gh + 528);

    auto gs_y_yyyy_xxxyz = pbuffer.data(idx_g_gh + 529);

    auto gs_y_yyyy_xxxzz = pbuffer.data(idx_g_gh + 530);

    auto gs_y_yyyy_xxyyy = pbuffer.data(idx_g_gh + 531);

    auto gs_y_yyyy_xxyyz = pbuffer.data(idx_g_gh + 532);

    auto gs_y_yyyy_xxyzz = pbuffer.data(idx_g_gh + 533);

    auto gs_y_yyyy_xxzzz = pbuffer.data(idx_g_gh + 534);

    auto gs_y_yyyy_xyyyy = pbuffer.data(idx_g_gh + 535);

    auto gs_y_yyyy_xyyyz = pbuffer.data(idx_g_gh + 536);

    auto gs_y_yyyy_xyyzz = pbuffer.data(idx_g_gh + 537);

    auto gs_y_yyyy_xyzzz = pbuffer.data(idx_g_gh + 538);

    auto gs_y_yyyy_xzzzz = pbuffer.data(idx_g_gh + 539);

    auto gs_y_yyyy_yyyyy = pbuffer.data(idx_g_gh + 540);

    auto gs_y_yyyy_yyyyz = pbuffer.data(idx_g_gh + 541);

    auto gs_y_yyyy_yyyzz = pbuffer.data(idx_g_gh + 542);

    auto gs_y_yyyy_yyzzz = pbuffer.data(idx_g_gh + 543);

    auto gs_y_yyyy_yzzzz = pbuffer.data(idx_g_gh + 544);

    auto gs_y_yyyy_zzzzz = pbuffer.data(idx_g_gh + 545);

    #pragma omp simd aligned(gc_y, gs_y_yyyy_xxxxx, gs_y_yyyy_xxxxy, gs_y_yyyy_xxxxz, gs_y_yyyy_xxxyy, gs_y_yyyy_xxxyz, gs_y_yyyy_xxxzz, gs_y_yyyy_xxyyy, gs_y_yyyy_xxyyz, gs_y_yyyy_xxyzz, gs_y_yyyy_xxzzz, gs_y_yyyy_xyyyy, gs_y_yyyy_xyyyz, gs_y_yyyy_xyyzz, gs_y_yyyy_xyzzz, gs_y_yyyy_xzzzz, gs_y_yyyy_yyyyy, gs_y_yyyy_yyyyz, gs_y_yyyy_yyyzz, gs_y_yyyy_yyzzz, gs_y_yyyy_yzzzz, gs_y_yyyy_zzzzz, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxzz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xxzzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_xzzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, ts_yyyy_xxxx, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxz, ts_yyyy_xxxzz, ts_yyyy_xxyy, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyz, ts_yyyy_xxyzz, ts_yyyy_xxzz, ts_yyyy_xxzzz, ts_yyyy_xyyy, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyz, ts_yyyy_xyyzz, ts_yyyy_xyzz, ts_yyyy_xyzzz, ts_yyyy_xzzz, ts_yyyy_xzzzz, ts_yyyy_yyyy, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyz, ts_yyyy_yyyzz, ts_yyyy_yyzz, ts_yyyy_yyzzz, ts_yyyy_yzzz, ts_yyyy_yzzzz, ts_yyyy_zzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyy_xxxxx[i] = 8.0 * ts_yyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxxy[i] = 8.0 * ts_yyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxxz[i] = 8.0 * ts_yyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxyy[i] = 8.0 * ts_yyy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxyz[i] = 8.0 * ts_yyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxzz[i] = 8.0 * ts_yyy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxyyy[i] = 8.0 * ts_yyy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxyyz[i] = 8.0 * ts_yyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxyzz[i] = 8.0 * ts_yyy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxzzz[i] = 8.0 * ts_yyy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyyyy[i] = 8.0 * ts_yyy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyyyz[i] = 8.0 * ts_yyy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyyzz[i] = 8.0 * ts_yyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyzzz[i] = 8.0 * ts_yyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xzzzz[i] = 8.0 * ts_yyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyyyy[i] = 8.0 * ts_yyy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyyyz[i] = 8.0 * ts_yyy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyyzz[i] = 8.0 * ts_yyy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyzzz[i] = 8.0 * ts_yyy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yzzzz[i] = 8.0 * ts_yyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_zzzzz[i] = 8.0 * ts_yyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 546-567 components of targeted buffer : GH

    auto gs_y_yyyz_xxxxx = pbuffer.data(idx_g_gh + 546);

    auto gs_y_yyyz_xxxxy = pbuffer.data(idx_g_gh + 547);

    auto gs_y_yyyz_xxxxz = pbuffer.data(idx_g_gh + 548);

    auto gs_y_yyyz_xxxyy = pbuffer.data(idx_g_gh + 549);

    auto gs_y_yyyz_xxxyz = pbuffer.data(idx_g_gh + 550);

    auto gs_y_yyyz_xxxzz = pbuffer.data(idx_g_gh + 551);

    auto gs_y_yyyz_xxyyy = pbuffer.data(idx_g_gh + 552);

    auto gs_y_yyyz_xxyyz = pbuffer.data(idx_g_gh + 553);

    auto gs_y_yyyz_xxyzz = pbuffer.data(idx_g_gh + 554);

    auto gs_y_yyyz_xxzzz = pbuffer.data(idx_g_gh + 555);

    auto gs_y_yyyz_xyyyy = pbuffer.data(idx_g_gh + 556);

    auto gs_y_yyyz_xyyyz = pbuffer.data(idx_g_gh + 557);

    auto gs_y_yyyz_xyyzz = pbuffer.data(idx_g_gh + 558);

    auto gs_y_yyyz_xyzzz = pbuffer.data(idx_g_gh + 559);

    auto gs_y_yyyz_xzzzz = pbuffer.data(idx_g_gh + 560);

    auto gs_y_yyyz_yyyyy = pbuffer.data(idx_g_gh + 561);

    auto gs_y_yyyz_yyyyz = pbuffer.data(idx_g_gh + 562);

    auto gs_y_yyyz_yyyzz = pbuffer.data(idx_g_gh + 563);

    auto gs_y_yyyz_yyzzz = pbuffer.data(idx_g_gh + 564);

    auto gs_y_yyyz_yzzzz = pbuffer.data(idx_g_gh + 565);

    auto gs_y_yyyz_zzzzz = pbuffer.data(idx_g_gh + 566);

    #pragma omp simd aligned(gc_y, gs_y_yyyz_xxxxx, gs_y_yyyz_xxxxy, gs_y_yyyz_xxxxz, gs_y_yyyz_xxxyy, gs_y_yyyz_xxxyz, gs_y_yyyz_xxxzz, gs_y_yyyz_xxyyy, gs_y_yyyz_xxyyz, gs_y_yyyz_xxyzz, gs_y_yyyz_xxzzz, gs_y_yyyz_xyyyy, gs_y_yyyz_xyyyz, gs_y_yyyz_xyyzz, gs_y_yyyz_xyzzz, gs_y_yyyz_xzzzz, gs_y_yyyz_yyyyy, gs_y_yyyz_yyyyz, gs_y_yyyz_yyyzz, gs_y_yyyz_yyzzz, gs_y_yyyz_yzzzz, gs_y_yyyz_zzzzz, ts_yyyz_xxxx, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxy, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxz, ts_yyyz_xxxzz, ts_yyyz_xxyy, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyz, ts_yyyz_xxyzz, ts_yyyz_xxzz, ts_yyyz_xxzzz, ts_yyyz_xyyy, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyz, ts_yyyz_xyyzz, ts_yyyz_xyzz, ts_yyyz_xyzzz, ts_yyyz_xzzz, ts_yyyz_xzzzz, ts_yyyz_yyyy, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyz, ts_yyyz_yyyzz, ts_yyyz_yyzz, ts_yyyz_yyzzz, ts_yyyz_yzzz, ts_yyyz_yzzzz, ts_yyyz_zzzz, ts_yyyz_zzzzz, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxzz, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyzz, ts_yyz_xxzzz, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyzz, ts_yyz_xyzzz, ts_yyz_xzzzz, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyzz, ts_yyz_yyzzz, ts_yyz_yzzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyz_xxxxx[i] = 6.0 * ts_yyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxxy[i] = 6.0 * ts_yyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxxz[i] = 6.0 * ts_yyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxyy[i] = 6.0 * ts_yyz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxyz[i] = 6.0 * ts_yyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxzz[i] = 6.0 * ts_yyz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxyyy[i] = 6.0 * ts_yyz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxyyz[i] = 6.0 * ts_yyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxyzz[i] = 6.0 * ts_yyz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxzzz[i] = 6.0 * ts_yyz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyyyy[i] = 6.0 * ts_yyz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyyyz[i] = 6.0 * ts_yyz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyyzz[i] = 6.0 * ts_yyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyzzz[i] = 6.0 * ts_yyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xzzzz[i] = 6.0 * ts_yyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyyyy[i] = 6.0 * ts_yyz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyyyz[i] = 6.0 * ts_yyz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyyzz[i] = 6.0 * ts_yyz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyzzz[i] = 6.0 * ts_yyz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yzzzz[i] = 6.0 * ts_yyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_zzzzz[i] = 6.0 * ts_yyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 567-588 components of targeted buffer : GH

    auto gs_y_yyzz_xxxxx = pbuffer.data(idx_g_gh + 567);

    auto gs_y_yyzz_xxxxy = pbuffer.data(idx_g_gh + 568);

    auto gs_y_yyzz_xxxxz = pbuffer.data(idx_g_gh + 569);

    auto gs_y_yyzz_xxxyy = pbuffer.data(idx_g_gh + 570);

    auto gs_y_yyzz_xxxyz = pbuffer.data(idx_g_gh + 571);

    auto gs_y_yyzz_xxxzz = pbuffer.data(idx_g_gh + 572);

    auto gs_y_yyzz_xxyyy = pbuffer.data(idx_g_gh + 573);

    auto gs_y_yyzz_xxyyz = pbuffer.data(idx_g_gh + 574);

    auto gs_y_yyzz_xxyzz = pbuffer.data(idx_g_gh + 575);

    auto gs_y_yyzz_xxzzz = pbuffer.data(idx_g_gh + 576);

    auto gs_y_yyzz_xyyyy = pbuffer.data(idx_g_gh + 577);

    auto gs_y_yyzz_xyyyz = pbuffer.data(idx_g_gh + 578);

    auto gs_y_yyzz_xyyzz = pbuffer.data(idx_g_gh + 579);

    auto gs_y_yyzz_xyzzz = pbuffer.data(idx_g_gh + 580);

    auto gs_y_yyzz_xzzzz = pbuffer.data(idx_g_gh + 581);

    auto gs_y_yyzz_yyyyy = pbuffer.data(idx_g_gh + 582);

    auto gs_y_yyzz_yyyyz = pbuffer.data(idx_g_gh + 583);

    auto gs_y_yyzz_yyyzz = pbuffer.data(idx_g_gh + 584);

    auto gs_y_yyzz_yyzzz = pbuffer.data(idx_g_gh + 585);

    auto gs_y_yyzz_yzzzz = pbuffer.data(idx_g_gh + 586);

    auto gs_y_yyzz_zzzzz = pbuffer.data(idx_g_gh + 587);

    #pragma omp simd aligned(gc_y, gs_y_yyzz_xxxxx, gs_y_yyzz_xxxxy, gs_y_yyzz_xxxxz, gs_y_yyzz_xxxyy, gs_y_yyzz_xxxyz, gs_y_yyzz_xxxzz, gs_y_yyzz_xxyyy, gs_y_yyzz_xxyyz, gs_y_yyzz_xxyzz, gs_y_yyzz_xxzzz, gs_y_yyzz_xyyyy, gs_y_yyzz_xyyyz, gs_y_yyzz_xyyzz, gs_y_yyzz_xyzzz, gs_y_yyzz_xzzzz, gs_y_yyzz_yyyyy, gs_y_yyzz_yyyyz, gs_y_yyzz_yyyzz, gs_y_yyzz_yyzzz, gs_y_yyzz_yzzzz, gs_y_yyzz_zzzzz, ts_yyzz_xxxx, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxy, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxz, ts_yyzz_xxxzz, ts_yyzz_xxyy, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyz, ts_yyzz_xxyzz, ts_yyzz_xxzz, ts_yyzz_xxzzz, ts_yyzz_xyyy, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyz, ts_yyzz_xyyzz, ts_yyzz_xyzz, ts_yyzz_xyzzz, ts_yyzz_xzzz, ts_yyzz_xzzzz, ts_yyzz_yyyy, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyz, ts_yyzz_yyyzz, ts_yyzz_yyzz, ts_yyzz_yyzzz, ts_yyzz_yzzz, ts_yyzz_yzzzz, ts_yyzz_zzzz, ts_yyzz_zzzzz, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxzz, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyzz, ts_yzz_xxzzz, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyzz, ts_yzz_xyzzz, ts_yzz_xzzzz, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzz_xxxxx[i] = 4.0 * ts_yzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxxy[i] = 4.0 * ts_yzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxxz[i] = 4.0 * ts_yzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxyy[i] = 4.0 * ts_yzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxyz[i] = 4.0 * ts_yzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxzz[i] = 4.0 * ts_yzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxyyy[i] = 4.0 * ts_yzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxyyz[i] = 4.0 * ts_yzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxyzz[i] = 4.0 * ts_yzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxzzz[i] = 4.0 * ts_yzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyyyy[i] = 4.0 * ts_yzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyyyz[i] = 4.0 * ts_yzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyyzz[i] = 4.0 * ts_yzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyzzz[i] = 4.0 * ts_yzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xzzzz[i] = 4.0 * ts_yzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyyyy[i] = 4.0 * ts_yzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyyyz[i] = 4.0 * ts_yzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyyzz[i] = 4.0 * ts_yzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyzzz[i] = 4.0 * ts_yzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yzzzz[i] = 4.0 * ts_yzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_zzzzz[i] = 4.0 * ts_yzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 588-609 components of targeted buffer : GH

    auto gs_y_yzzz_xxxxx = pbuffer.data(idx_g_gh + 588);

    auto gs_y_yzzz_xxxxy = pbuffer.data(idx_g_gh + 589);

    auto gs_y_yzzz_xxxxz = pbuffer.data(idx_g_gh + 590);

    auto gs_y_yzzz_xxxyy = pbuffer.data(idx_g_gh + 591);

    auto gs_y_yzzz_xxxyz = pbuffer.data(idx_g_gh + 592);

    auto gs_y_yzzz_xxxzz = pbuffer.data(idx_g_gh + 593);

    auto gs_y_yzzz_xxyyy = pbuffer.data(idx_g_gh + 594);

    auto gs_y_yzzz_xxyyz = pbuffer.data(idx_g_gh + 595);

    auto gs_y_yzzz_xxyzz = pbuffer.data(idx_g_gh + 596);

    auto gs_y_yzzz_xxzzz = pbuffer.data(idx_g_gh + 597);

    auto gs_y_yzzz_xyyyy = pbuffer.data(idx_g_gh + 598);

    auto gs_y_yzzz_xyyyz = pbuffer.data(idx_g_gh + 599);

    auto gs_y_yzzz_xyyzz = pbuffer.data(idx_g_gh + 600);

    auto gs_y_yzzz_xyzzz = pbuffer.data(idx_g_gh + 601);

    auto gs_y_yzzz_xzzzz = pbuffer.data(idx_g_gh + 602);

    auto gs_y_yzzz_yyyyy = pbuffer.data(idx_g_gh + 603);

    auto gs_y_yzzz_yyyyz = pbuffer.data(idx_g_gh + 604);

    auto gs_y_yzzz_yyyzz = pbuffer.data(idx_g_gh + 605);

    auto gs_y_yzzz_yyzzz = pbuffer.data(idx_g_gh + 606);

    auto gs_y_yzzz_yzzzz = pbuffer.data(idx_g_gh + 607);

    auto gs_y_yzzz_zzzzz = pbuffer.data(idx_g_gh + 608);

    #pragma omp simd aligned(gc_y, gs_y_yzzz_xxxxx, gs_y_yzzz_xxxxy, gs_y_yzzz_xxxxz, gs_y_yzzz_xxxyy, gs_y_yzzz_xxxyz, gs_y_yzzz_xxxzz, gs_y_yzzz_xxyyy, gs_y_yzzz_xxyyz, gs_y_yzzz_xxyzz, gs_y_yzzz_xxzzz, gs_y_yzzz_xyyyy, gs_y_yzzz_xyyyz, gs_y_yzzz_xyyzz, gs_y_yzzz_xyzzz, gs_y_yzzz_xzzzz, gs_y_yzzz_yyyyy, gs_y_yzzz_yyyyz, gs_y_yzzz_yyyzz, gs_y_yzzz_yyzzz, gs_y_yzzz_yzzzz, gs_y_yzzz_zzzzz, ts_yzzz_xxxx, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxy, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxz, ts_yzzz_xxxzz, ts_yzzz_xxyy, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyz, ts_yzzz_xxyzz, ts_yzzz_xxzz, ts_yzzz_xxzzz, ts_yzzz_xyyy, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyz, ts_yzzz_xyyzz, ts_yzzz_xyzz, ts_yzzz_xyzzz, ts_yzzz_xzzz, ts_yzzz_xzzzz, ts_yzzz_yyyy, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyz, ts_yzzz_yyyzz, ts_yzzz_yyzz, ts_yzzz_yyzzz, ts_yzzz_yzzz, ts_yzzz_yzzzz, ts_yzzz_zzzz, ts_yzzz_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzz_xxxxx[i] = 2.0 * ts_zzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxxy[i] = 2.0 * ts_zzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxxz[i] = 2.0 * ts_zzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxyy[i] = 2.0 * ts_zzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxyz[i] = 2.0 * ts_zzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxzz[i] = 2.0 * ts_zzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxyyy[i] = 2.0 * ts_zzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxyyz[i] = 2.0 * ts_zzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxyzz[i] = 2.0 * ts_zzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxzzz[i] = 2.0 * ts_zzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyyyy[i] = 2.0 * ts_zzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyyyz[i] = 2.0 * ts_zzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyyzz[i] = 2.0 * ts_zzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyzzz[i] = 2.0 * ts_zzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xzzzz[i] = 2.0 * ts_zzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyyyy[i] = 2.0 * ts_zzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyyyz[i] = 2.0 * ts_zzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyyzz[i] = 2.0 * ts_zzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyzzz[i] = 2.0 * ts_zzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yzzzz[i] = 2.0 * ts_zzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_zzzzz[i] = 2.0 * ts_zzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 609-630 components of targeted buffer : GH

    auto gs_y_zzzz_xxxxx = pbuffer.data(idx_g_gh + 609);

    auto gs_y_zzzz_xxxxy = pbuffer.data(idx_g_gh + 610);

    auto gs_y_zzzz_xxxxz = pbuffer.data(idx_g_gh + 611);

    auto gs_y_zzzz_xxxyy = pbuffer.data(idx_g_gh + 612);

    auto gs_y_zzzz_xxxyz = pbuffer.data(idx_g_gh + 613);

    auto gs_y_zzzz_xxxzz = pbuffer.data(idx_g_gh + 614);

    auto gs_y_zzzz_xxyyy = pbuffer.data(idx_g_gh + 615);

    auto gs_y_zzzz_xxyyz = pbuffer.data(idx_g_gh + 616);

    auto gs_y_zzzz_xxyzz = pbuffer.data(idx_g_gh + 617);

    auto gs_y_zzzz_xxzzz = pbuffer.data(idx_g_gh + 618);

    auto gs_y_zzzz_xyyyy = pbuffer.data(idx_g_gh + 619);

    auto gs_y_zzzz_xyyyz = pbuffer.data(idx_g_gh + 620);

    auto gs_y_zzzz_xyyzz = pbuffer.data(idx_g_gh + 621);

    auto gs_y_zzzz_xyzzz = pbuffer.data(idx_g_gh + 622);

    auto gs_y_zzzz_xzzzz = pbuffer.data(idx_g_gh + 623);

    auto gs_y_zzzz_yyyyy = pbuffer.data(idx_g_gh + 624);

    auto gs_y_zzzz_yyyyz = pbuffer.data(idx_g_gh + 625);

    auto gs_y_zzzz_yyyzz = pbuffer.data(idx_g_gh + 626);

    auto gs_y_zzzz_yyzzz = pbuffer.data(idx_g_gh + 627);

    auto gs_y_zzzz_yzzzz = pbuffer.data(idx_g_gh + 628);

    auto gs_y_zzzz_zzzzz = pbuffer.data(idx_g_gh + 629);

    #pragma omp simd aligned(gc_y, gs_y_zzzz_xxxxx, gs_y_zzzz_xxxxy, gs_y_zzzz_xxxxz, gs_y_zzzz_xxxyy, gs_y_zzzz_xxxyz, gs_y_zzzz_xxxzz, gs_y_zzzz_xxyyy, gs_y_zzzz_xxyyz, gs_y_zzzz_xxyzz, gs_y_zzzz_xxzzz, gs_y_zzzz_xyyyy, gs_y_zzzz_xyyyz, gs_y_zzzz_xyyzz, gs_y_zzzz_xyzzz, gs_y_zzzz_xzzzz, gs_y_zzzz_yyyyy, gs_y_zzzz_yyyyz, gs_y_zzzz_yyyzz, gs_y_zzzz_yyzzz, gs_y_zzzz_yzzzz, gs_y_zzzz_zzzzz, ts_zzzz_xxxx, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxy, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxz, ts_zzzz_xxxzz, ts_zzzz_xxyy, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyz, ts_zzzz_xxyzz, ts_zzzz_xxzz, ts_zzzz_xxzzz, ts_zzzz_xyyy, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyz, ts_zzzz_xyyzz, ts_zzzz_xyzz, ts_zzzz_xyzzz, ts_zzzz_xzzz, ts_zzzz_xzzzz, ts_zzzz_yyyy, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyz, ts_zzzz_yyyzz, ts_zzzz_yyzz, ts_zzzz_yyzzz, ts_zzzz_yzzz, ts_zzzz_yzzzz, ts_zzzz_zzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzz_xxxxx[i] = 2.0 * ts_zzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxxy[i] = 2.0 * ts_zzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxxz[i] = 2.0 * ts_zzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxyy[i] = 4.0 * ts_zzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxyz[i] = 2.0 * ts_zzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxzz[i] = 2.0 * ts_zzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxyyy[i] = 6.0 * ts_zzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxyyz[i] = 4.0 * ts_zzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxyzz[i] = 2.0 * ts_zzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxzzz[i] = 2.0 * ts_zzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyyyy[i] = 8.0 * ts_zzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyyyz[i] = 6.0 * ts_zzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyyzz[i] = 4.0 * ts_zzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyzzz[i] = 2.0 * ts_zzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xzzzz[i] = 2.0 * ts_zzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyyyy[i] = 10.0 * ts_zzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyyyz[i] = 8.0 * ts_zzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyyzz[i] = 6.0 * ts_zzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyzzz[i] = 4.0 * ts_zzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yzzzz[i] = 2.0 * ts_zzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_zzzzz[i] = 2.0 * ts_zzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 630-651 components of targeted buffer : GH

    auto gs_z_xxxx_xxxxx = pbuffer.data(idx_g_gh + 630);

    auto gs_z_xxxx_xxxxy = pbuffer.data(idx_g_gh + 631);

    auto gs_z_xxxx_xxxxz = pbuffer.data(idx_g_gh + 632);

    auto gs_z_xxxx_xxxyy = pbuffer.data(idx_g_gh + 633);

    auto gs_z_xxxx_xxxyz = pbuffer.data(idx_g_gh + 634);

    auto gs_z_xxxx_xxxzz = pbuffer.data(idx_g_gh + 635);

    auto gs_z_xxxx_xxyyy = pbuffer.data(idx_g_gh + 636);

    auto gs_z_xxxx_xxyyz = pbuffer.data(idx_g_gh + 637);

    auto gs_z_xxxx_xxyzz = pbuffer.data(idx_g_gh + 638);

    auto gs_z_xxxx_xxzzz = pbuffer.data(idx_g_gh + 639);

    auto gs_z_xxxx_xyyyy = pbuffer.data(idx_g_gh + 640);

    auto gs_z_xxxx_xyyyz = pbuffer.data(idx_g_gh + 641);

    auto gs_z_xxxx_xyyzz = pbuffer.data(idx_g_gh + 642);

    auto gs_z_xxxx_xyzzz = pbuffer.data(idx_g_gh + 643);

    auto gs_z_xxxx_xzzzz = pbuffer.data(idx_g_gh + 644);

    auto gs_z_xxxx_yyyyy = pbuffer.data(idx_g_gh + 645);

    auto gs_z_xxxx_yyyyz = pbuffer.data(idx_g_gh + 646);

    auto gs_z_xxxx_yyyzz = pbuffer.data(idx_g_gh + 647);

    auto gs_z_xxxx_yyzzz = pbuffer.data(idx_g_gh + 648);

    auto gs_z_xxxx_yzzzz = pbuffer.data(idx_g_gh + 649);

    auto gs_z_xxxx_zzzzz = pbuffer.data(idx_g_gh + 650);

    #pragma omp simd aligned(gc_z, gs_z_xxxx_xxxxx, gs_z_xxxx_xxxxy, gs_z_xxxx_xxxxz, gs_z_xxxx_xxxyy, gs_z_xxxx_xxxyz, gs_z_xxxx_xxxzz, gs_z_xxxx_xxyyy, gs_z_xxxx_xxyyz, gs_z_xxxx_xxyzz, gs_z_xxxx_xxzzz, gs_z_xxxx_xyyyy, gs_z_xxxx_xyyyz, gs_z_xxxx_xyyzz, gs_z_xxxx_xyzzz, gs_z_xxxx_xzzzz, gs_z_xxxx_yyyyy, gs_z_xxxx_yyyyz, gs_z_xxxx_yyyzz, gs_z_xxxx_yyzzz, gs_z_xxxx_yzzzz, gs_z_xxxx_zzzzz, ts_xxxx_xxxx, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxy, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxz, ts_xxxx_xxxzz, ts_xxxx_xxyy, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyz, ts_xxxx_xxyzz, ts_xxxx_xxzz, ts_xxxx_xxzzz, ts_xxxx_xyyy, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyz, ts_xxxx_xyyzz, ts_xxxx_xyzz, ts_xxxx_xyzzz, ts_xxxx_xzzz, ts_xxxx_xzzzz, ts_xxxx_yyyy, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyz, ts_xxxx_yyyzz, ts_xxxx_yyzz, ts_xxxx_yyzzz, ts_xxxx_yzzz, ts_xxxx_yzzzz, ts_xxxx_zzzz, ts_xxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxx_xxxxx[i] = 2.0 * ts_xxxx_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxxy[i] = 2.0 * ts_xxxx_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxxz[i] = 2.0 * ts_xxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxyy[i] = 2.0 * ts_xxxx_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxyz[i] = 2.0 * ts_xxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxzz[i] = 4.0 * ts_xxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxyyy[i] = 2.0 * ts_xxxx_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxyyz[i] = 2.0 * ts_xxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxyzz[i] = 4.0 * ts_xxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxzzz[i] = 6.0 * ts_xxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyyyy[i] = 2.0 * ts_xxxx_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyyyz[i] = 2.0 * ts_xxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyyzz[i] = 4.0 * ts_xxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyzzz[i] = 6.0 * ts_xxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xzzzz[i] = 8.0 * ts_xxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyyyy[i] = 2.0 * ts_xxxx_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyyyz[i] = 2.0 * ts_xxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyyzz[i] = 4.0 * ts_xxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyzzz[i] = 6.0 * ts_xxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yzzzz[i] = 8.0 * ts_xxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_zzzzz[i] = 10.0 * ts_xxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 651-672 components of targeted buffer : GH

    auto gs_z_xxxy_xxxxx = pbuffer.data(idx_g_gh + 651);

    auto gs_z_xxxy_xxxxy = pbuffer.data(idx_g_gh + 652);

    auto gs_z_xxxy_xxxxz = pbuffer.data(idx_g_gh + 653);

    auto gs_z_xxxy_xxxyy = pbuffer.data(idx_g_gh + 654);

    auto gs_z_xxxy_xxxyz = pbuffer.data(idx_g_gh + 655);

    auto gs_z_xxxy_xxxzz = pbuffer.data(idx_g_gh + 656);

    auto gs_z_xxxy_xxyyy = pbuffer.data(idx_g_gh + 657);

    auto gs_z_xxxy_xxyyz = pbuffer.data(idx_g_gh + 658);

    auto gs_z_xxxy_xxyzz = pbuffer.data(idx_g_gh + 659);

    auto gs_z_xxxy_xxzzz = pbuffer.data(idx_g_gh + 660);

    auto gs_z_xxxy_xyyyy = pbuffer.data(idx_g_gh + 661);

    auto gs_z_xxxy_xyyyz = pbuffer.data(idx_g_gh + 662);

    auto gs_z_xxxy_xyyzz = pbuffer.data(idx_g_gh + 663);

    auto gs_z_xxxy_xyzzz = pbuffer.data(idx_g_gh + 664);

    auto gs_z_xxxy_xzzzz = pbuffer.data(idx_g_gh + 665);

    auto gs_z_xxxy_yyyyy = pbuffer.data(idx_g_gh + 666);

    auto gs_z_xxxy_yyyyz = pbuffer.data(idx_g_gh + 667);

    auto gs_z_xxxy_yyyzz = pbuffer.data(idx_g_gh + 668);

    auto gs_z_xxxy_yyzzz = pbuffer.data(idx_g_gh + 669);

    auto gs_z_xxxy_yzzzz = pbuffer.data(idx_g_gh + 670);

    auto gs_z_xxxy_zzzzz = pbuffer.data(idx_g_gh + 671);

    #pragma omp simd aligned(gc_z, gs_z_xxxy_xxxxx, gs_z_xxxy_xxxxy, gs_z_xxxy_xxxxz, gs_z_xxxy_xxxyy, gs_z_xxxy_xxxyz, gs_z_xxxy_xxxzz, gs_z_xxxy_xxyyy, gs_z_xxxy_xxyyz, gs_z_xxxy_xxyzz, gs_z_xxxy_xxzzz, gs_z_xxxy_xyyyy, gs_z_xxxy_xyyyz, gs_z_xxxy_xyyzz, gs_z_xxxy_xyzzz, gs_z_xxxy_xzzzz, gs_z_xxxy_yyyyy, gs_z_xxxy_yyyyz, gs_z_xxxy_yyyzz, gs_z_xxxy_yyzzz, gs_z_xxxy_yzzzz, gs_z_xxxy_zzzzz, ts_xxxy_xxxx, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxy, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxz, ts_xxxy_xxxzz, ts_xxxy_xxyy, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyz, ts_xxxy_xxyzz, ts_xxxy_xxzz, ts_xxxy_xxzzz, ts_xxxy_xyyy, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyz, ts_xxxy_xyyzz, ts_xxxy_xyzz, ts_xxxy_xyzzz, ts_xxxy_xzzz, ts_xxxy_xzzzz, ts_xxxy_yyyy, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyz, ts_xxxy_yyyzz, ts_xxxy_yyzz, ts_xxxy_yyzzz, ts_xxxy_yzzz, ts_xxxy_yzzzz, ts_xxxy_zzzz, ts_xxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxy_xxxxx[i] = 2.0 * ts_xxxy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxxy[i] = 2.0 * ts_xxxy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxxz[i] = 2.0 * ts_xxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxyy[i] = 2.0 * ts_xxxy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxyz[i] = 2.0 * ts_xxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxzz[i] = 4.0 * ts_xxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxyyy[i] = 2.0 * ts_xxxy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxyyz[i] = 2.0 * ts_xxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxyzz[i] = 4.0 * ts_xxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxzzz[i] = 6.0 * ts_xxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyyyy[i] = 2.0 * ts_xxxy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyyyz[i] = 2.0 * ts_xxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyyzz[i] = 4.0 * ts_xxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyzzz[i] = 6.0 * ts_xxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xzzzz[i] = 8.0 * ts_xxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyyyy[i] = 2.0 * ts_xxxy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyyyz[i] = 2.0 * ts_xxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyyzz[i] = 4.0 * ts_xxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyzzz[i] = 6.0 * ts_xxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yzzzz[i] = 8.0 * ts_xxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_zzzzz[i] = 10.0 * ts_xxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 672-693 components of targeted buffer : GH

    auto gs_z_xxxz_xxxxx = pbuffer.data(idx_g_gh + 672);

    auto gs_z_xxxz_xxxxy = pbuffer.data(idx_g_gh + 673);

    auto gs_z_xxxz_xxxxz = pbuffer.data(idx_g_gh + 674);

    auto gs_z_xxxz_xxxyy = pbuffer.data(idx_g_gh + 675);

    auto gs_z_xxxz_xxxyz = pbuffer.data(idx_g_gh + 676);

    auto gs_z_xxxz_xxxzz = pbuffer.data(idx_g_gh + 677);

    auto gs_z_xxxz_xxyyy = pbuffer.data(idx_g_gh + 678);

    auto gs_z_xxxz_xxyyz = pbuffer.data(idx_g_gh + 679);

    auto gs_z_xxxz_xxyzz = pbuffer.data(idx_g_gh + 680);

    auto gs_z_xxxz_xxzzz = pbuffer.data(idx_g_gh + 681);

    auto gs_z_xxxz_xyyyy = pbuffer.data(idx_g_gh + 682);

    auto gs_z_xxxz_xyyyz = pbuffer.data(idx_g_gh + 683);

    auto gs_z_xxxz_xyyzz = pbuffer.data(idx_g_gh + 684);

    auto gs_z_xxxz_xyzzz = pbuffer.data(idx_g_gh + 685);

    auto gs_z_xxxz_xzzzz = pbuffer.data(idx_g_gh + 686);

    auto gs_z_xxxz_yyyyy = pbuffer.data(idx_g_gh + 687);

    auto gs_z_xxxz_yyyyz = pbuffer.data(idx_g_gh + 688);

    auto gs_z_xxxz_yyyzz = pbuffer.data(idx_g_gh + 689);

    auto gs_z_xxxz_yyzzz = pbuffer.data(idx_g_gh + 690);

    auto gs_z_xxxz_yzzzz = pbuffer.data(idx_g_gh + 691);

    auto gs_z_xxxz_zzzzz = pbuffer.data(idx_g_gh + 692);

    #pragma omp simd aligned(gc_z, gs_z_xxxz_xxxxx, gs_z_xxxz_xxxxy, gs_z_xxxz_xxxxz, gs_z_xxxz_xxxyy, gs_z_xxxz_xxxyz, gs_z_xxxz_xxxzz, gs_z_xxxz_xxyyy, gs_z_xxxz_xxyyz, gs_z_xxxz_xxyzz, gs_z_xxxz_xxzzz, gs_z_xxxz_xyyyy, gs_z_xxxz_xyyyz, gs_z_xxxz_xyyzz, gs_z_xxxz_xyzzz, gs_z_xxxz_xzzzz, gs_z_xxxz_yyyyy, gs_z_xxxz_yyyyz, gs_z_xxxz_yyyzz, gs_z_xxxz_yyzzz, gs_z_xxxz_yzzzz, gs_z_xxxz_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxzz, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyzz, ts_xxx_xxzzz, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyzz, ts_xxx_xyzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyzz, ts_xxx_yyzzz, ts_xxx_yzzzz, ts_xxx_zzzzz, ts_xxxz_xxxx, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxy, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxz, ts_xxxz_xxxzz, ts_xxxz_xxyy, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyz, ts_xxxz_xxyzz, ts_xxxz_xxzz, ts_xxxz_xxzzz, ts_xxxz_xyyy, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyz, ts_xxxz_xyyzz, ts_xxxz_xyzz, ts_xxxz_xyzzz, ts_xxxz_xzzz, ts_xxxz_xzzzz, ts_xxxz_yyyy, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyz, ts_xxxz_yyyzz, ts_xxxz_yyzz, ts_xxxz_yyzzz, ts_xxxz_yzzz, ts_xxxz_yzzzz, ts_xxxz_zzzz, ts_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxz_xxxxx[i] = 2.0 * ts_xxx_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxxy[i] = 2.0 * ts_xxx_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxxz[i] = 2.0 * ts_xxx_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxyy[i] = 2.0 * ts_xxx_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxyz[i] = 2.0 * ts_xxx_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxzz[i] = 2.0 * ts_xxx_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxyyy[i] = 2.0 * ts_xxx_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxyyz[i] = 2.0 * ts_xxx_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxyzz[i] = 2.0 * ts_xxx_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxzzz[i] = 2.0 * ts_xxx_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyyyy[i] = 2.0 * ts_xxx_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyyyz[i] = 2.0 * ts_xxx_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyyzz[i] = 2.0 * ts_xxx_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyzzz[i] = 2.0 * ts_xxx_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xzzzz[i] = 2.0 * ts_xxx_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyyyy[i] = 2.0 * ts_xxx_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyyyz[i] = 2.0 * ts_xxx_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyyzz[i] = 2.0 * ts_xxx_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyzzz[i] = 2.0 * ts_xxx_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yzzzz[i] = 2.0 * ts_xxx_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_zzzzz[i] = 2.0 * ts_xxx_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 693-714 components of targeted buffer : GH

    auto gs_z_xxyy_xxxxx = pbuffer.data(idx_g_gh + 693);

    auto gs_z_xxyy_xxxxy = pbuffer.data(idx_g_gh + 694);

    auto gs_z_xxyy_xxxxz = pbuffer.data(idx_g_gh + 695);

    auto gs_z_xxyy_xxxyy = pbuffer.data(idx_g_gh + 696);

    auto gs_z_xxyy_xxxyz = pbuffer.data(idx_g_gh + 697);

    auto gs_z_xxyy_xxxzz = pbuffer.data(idx_g_gh + 698);

    auto gs_z_xxyy_xxyyy = pbuffer.data(idx_g_gh + 699);

    auto gs_z_xxyy_xxyyz = pbuffer.data(idx_g_gh + 700);

    auto gs_z_xxyy_xxyzz = pbuffer.data(idx_g_gh + 701);

    auto gs_z_xxyy_xxzzz = pbuffer.data(idx_g_gh + 702);

    auto gs_z_xxyy_xyyyy = pbuffer.data(idx_g_gh + 703);

    auto gs_z_xxyy_xyyyz = pbuffer.data(idx_g_gh + 704);

    auto gs_z_xxyy_xyyzz = pbuffer.data(idx_g_gh + 705);

    auto gs_z_xxyy_xyzzz = pbuffer.data(idx_g_gh + 706);

    auto gs_z_xxyy_xzzzz = pbuffer.data(idx_g_gh + 707);

    auto gs_z_xxyy_yyyyy = pbuffer.data(idx_g_gh + 708);

    auto gs_z_xxyy_yyyyz = pbuffer.data(idx_g_gh + 709);

    auto gs_z_xxyy_yyyzz = pbuffer.data(idx_g_gh + 710);

    auto gs_z_xxyy_yyzzz = pbuffer.data(idx_g_gh + 711);

    auto gs_z_xxyy_yzzzz = pbuffer.data(idx_g_gh + 712);

    auto gs_z_xxyy_zzzzz = pbuffer.data(idx_g_gh + 713);

    #pragma omp simd aligned(gc_z, gs_z_xxyy_xxxxx, gs_z_xxyy_xxxxy, gs_z_xxyy_xxxxz, gs_z_xxyy_xxxyy, gs_z_xxyy_xxxyz, gs_z_xxyy_xxxzz, gs_z_xxyy_xxyyy, gs_z_xxyy_xxyyz, gs_z_xxyy_xxyzz, gs_z_xxyy_xxzzz, gs_z_xxyy_xyyyy, gs_z_xxyy_xyyyz, gs_z_xxyy_xyyzz, gs_z_xxyy_xyzzz, gs_z_xxyy_xzzzz, gs_z_xxyy_yyyyy, gs_z_xxyy_yyyyz, gs_z_xxyy_yyyzz, gs_z_xxyy_yyzzz, gs_z_xxyy_yzzzz, gs_z_xxyy_zzzzz, ts_xxyy_xxxx, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxz, ts_xxyy_xxxzz, ts_xxyy_xxyy, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyz, ts_xxyy_xxyzz, ts_xxyy_xxzz, ts_xxyy_xxzzz, ts_xxyy_xyyy, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyz, ts_xxyy_xyyzz, ts_xxyy_xyzz, ts_xxyy_xyzzz, ts_xxyy_xzzz, ts_xxyy_xzzzz, ts_xxyy_yyyy, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyz, ts_xxyy_yyyzz, ts_xxyy_yyzz, ts_xxyy_yyzzz, ts_xxyy_yzzz, ts_xxyy_yzzzz, ts_xxyy_zzzz, ts_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyy_xxxxx[i] = 2.0 * ts_xxyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxxy[i] = 2.0 * ts_xxyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxxz[i] = 2.0 * ts_xxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxyy[i] = 2.0 * ts_xxyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxyz[i] = 2.0 * ts_xxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxzz[i] = 4.0 * ts_xxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxyyy[i] = 2.0 * ts_xxyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxyyz[i] = 2.0 * ts_xxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxyzz[i] = 4.0 * ts_xxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxzzz[i] = 6.0 * ts_xxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyyyy[i] = 2.0 * ts_xxyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyyyz[i] = 2.0 * ts_xxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyyzz[i] = 4.0 * ts_xxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyzzz[i] = 6.0 * ts_xxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xzzzz[i] = 8.0 * ts_xxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyyyy[i] = 2.0 * ts_xxyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyyyz[i] = 2.0 * ts_xxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyyzz[i] = 4.0 * ts_xxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyzzz[i] = 6.0 * ts_xxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yzzzz[i] = 8.0 * ts_xxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_zzzzz[i] = 10.0 * ts_xxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 714-735 components of targeted buffer : GH

    auto gs_z_xxyz_xxxxx = pbuffer.data(idx_g_gh + 714);

    auto gs_z_xxyz_xxxxy = pbuffer.data(idx_g_gh + 715);

    auto gs_z_xxyz_xxxxz = pbuffer.data(idx_g_gh + 716);

    auto gs_z_xxyz_xxxyy = pbuffer.data(idx_g_gh + 717);

    auto gs_z_xxyz_xxxyz = pbuffer.data(idx_g_gh + 718);

    auto gs_z_xxyz_xxxzz = pbuffer.data(idx_g_gh + 719);

    auto gs_z_xxyz_xxyyy = pbuffer.data(idx_g_gh + 720);

    auto gs_z_xxyz_xxyyz = pbuffer.data(idx_g_gh + 721);

    auto gs_z_xxyz_xxyzz = pbuffer.data(idx_g_gh + 722);

    auto gs_z_xxyz_xxzzz = pbuffer.data(idx_g_gh + 723);

    auto gs_z_xxyz_xyyyy = pbuffer.data(idx_g_gh + 724);

    auto gs_z_xxyz_xyyyz = pbuffer.data(idx_g_gh + 725);

    auto gs_z_xxyz_xyyzz = pbuffer.data(idx_g_gh + 726);

    auto gs_z_xxyz_xyzzz = pbuffer.data(idx_g_gh + 727);

    auto gs_z_xxyz_xzzzz = pbuffer.data(idx_g_gh + 728);

    auto gs_z_xxyz_yyyyy = pbuffer.data(idx_g_gh + 729);

    auto gs_z_xxyz_yyyyz = pbuffer.data(idx_g_gh + 730);

    auto gs_z_xxyz_yyyzz = pbuffer.data(idx_g_gh + 731);

    auto gs_z_xxyz_yyzzz = pbuffer.data(idx_g_gh + 732);

    auto gs_z_xxyz_yzzzz = pbuffer.data(idx_g_gh + 733);

    auto gs_z_xxyz_zzzzz = pbuffer.data(idx_g_gh + 734);

    #pragma omp simd aligned(gc_z, gs_z_xxyz_xxxxx, gs_z_xxyz_xxxxy, gs_z_xxyz_xxxxz, gs_z_xxyz_xxxyy, gs_z_xxyz_xxxyz, gs_z_xxyz_xxxzz, gs_z_xxyz_xxyyy, gs_z_xxyz_xxyyz, gs_z_xxyz_xxyzz, gs_z_xxyz_xxzzz, gs_z_xxyz_xyyyy, gs_z_xxyz_xyyyz, gs_z_xxyz_xyyzz, gs_z_xxyz_xyzzz, gs_z_xxyz_xzzzz, gs_z_xxyz_yyyyy, gs_z_xxyz_yyyyz, gs_z_xxyz_yyyzz, gs_z_xxyz_yyzzz, gs_z_xxyz_yzzzz, gs_z_xxyz_zzzzz, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxzz, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyzz, ts_xxy_xxzzz, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyzz, ts_xxy_xyzzz, ts_xxy_xzzzz, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyzz, ts_xxy_yyzzz, ts_xxy_yzzzz, ts_xxy_zzzzz, ts_xxyz_xxxx, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxy, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxz, ts_xxyz_xxxzz, ts_xxyz_xxyy, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyz, ts_xxyz_xxyzz, ts_xxyz_xxzz, ts_xxyz_xxzzz, ts_xxyz_xyyy, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyz, ts_xxyz_xyyzz, ts_xxyz_xyzz, ts_xxyz_xyzzz, ts_xxyz_xzzz, ts_xxyz_xzzzz, ts_xxyz_yyyy, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyz, ts_xxyz_yyyzz, ts_xxyz_yyzz, ts_xxyz_yyzzz, ts_xxyz_yzzz, ts_xxyz_yzzzz, ts_xxyz_zzzz, ts_xxyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyz_xxxxx[i] = 2.0 * ts_xxy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxxy[i] = 2.0 * ts_xxy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxxz[i] = 2.0 * ts_xxy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxyy[i] = 2.0 * ts_xxy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxyz[i] = 2.0 * ts_xxy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxzz[i] = 2.0 * ts_xxy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxyyy[i] = 2.0 * ts_xxy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxyyz[i] = 2.0 * ts_xxy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxyzz[i] = 2.0 * ts_xxy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxzzz[i] = 2.0 * ts_xxy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyyyy[i] = 2.0 * ts_xxy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyyyz[i] = 2.0 * ts_xxy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyyzz[i] = 2.0 * ts_xxy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyzzz[i] = 2.0 * ts_xxy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xzzzz[i] = 2.0 * ts_xxy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyyyy[i] = 2.0 * ts_xxy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyyyz[i] = 2.0 * ts_xxy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyyzz[i] = 2.0 * ts_xxy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyzzz[i] = 2.0 * ts_xxy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yzzzz[i] = 2.0 * ts_xxy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_zzzzz[i] = 2.0 * ts_xxy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 735-756 components of targeted buffer : GH

    auto gs_z_xxzz_xxxxx = pbuffer.data(idx_g_gh + 735);

    auto gs_z_xxzz_xxxxy = pbuffer.data(idx_g_gh + 736);

    auto gs_z_xxzz_xxxxz = pbuffer.data(idx_g_gh + 737);

    auto gs_z_xxzz_xxxyy = pbuffer.data(idx_g_gh + 738);

    auto gs_z_xxzz_xxxyz = pbuffer.data(idx_g_gh + 739);

    auto gs_z_xxzz_xxxzz = pbuffer.data(idx_g_gh + 740);

    auto gs_z_xxzz_xxyyy = pbuffer.data(idx_g_gh + 741);

    auto gs_z_xxzz_xxyyz = pbuffer.data(idx_g_gh + 742);

    auto gs_z_xxzz_xxyzz = pbuffer.data(idx_g_gh + 743);

    auto gs_z_xxzz_xxzzz = pbuffer.data(idx_g_gh + 744);

    auto gs_z_xxzz_xyyyy = pbuffer.data(idx_g_gh + 745);

    auto gs_z_xxzz_xyyyz = pbuffer.data(idx_g_gh + 746);

    auto gs_z_xxzz_xyyzz = pbuffer.data(idx_g_gh + 747);

    auto gs_z_xxzz_xyzzz = pbuffer.data(idx_g_gh + 748);

    auto gs_z_xxzz_xzzzz = pbuffer.data(idx_g_gh + 749);

    auto gs_z_xxzz_yyyyy = pbuffer.data(idx_g_gh + 750);

    auto gs_z_xxzz_yyyyz = pbuffer.data(idx_g_gh + 751);

    auto gs_z_xxzz_yyyzz = pbuffer.data(idx_g_gh + 752);

    auto gs_z_xxzz_yyzzz = pbuffer.data(idx_g_gh + 753);

    auto gs_z_xxzz_yzzzz = pbuffer.data(idx_g_gh + 754);

    auto gs_z_xxzz_zzzzz = pbuffer.data(idx_g_gh + 755);

    #pragma omp simd aligned(gc_z, gs_z_xxzz_xxxxx, gs_z_xxzz_xxxxy, gs_z_xxzz_xxxxz, gs_z_xxzz_xxxyy, gs_z_xxzz_xxxyz, gs_z_xxzz_xxxzz, gs_z_xxzz_xxyyy, gs_z_xxzz_xxyyz, gs_z_xxzz_xxyzz, gs_z_xxzz_xxzzz, gs_z_xxzz_xyyyy, gs_z_xxzz_xyyyz, gs_z_xxzz_xyyzz, gs_z_xxzz_xyzzz, gs_z_xxzz_xzzzz, gs_z_xxzz_yyyyy, gs_z_xxzz_yyyyz, gs_z_xxzz_yyyzz, gs_z_xxzz_yyzzz, gs_z_xxzz_yzzzz, gs_z_xxzz_zzzzz, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxzz, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyzz, ts_xxz_xxzzz, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyzz, ts_xxz_xyzzz, ts_xxz_xzzzz, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyzz, ts_xxz_yyzzz, ts_xxz_yzzzz, ts_xxz_zzzzz, ts_xxzz_xxxx, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxy, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxz, ts_xxzz_xxxzz, ts_xxzz_xxyy, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyz, ts_xxzz_xxyzz, ts_xxzz_xxzz, ts_xxzz_xxzzz, ts_xxzz_xyyy, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyz, ts_xxzz_xyyzz, ts_xxzz_xyzz, ts_xxzz_xyzzz, ts_xxzz_xzzz, ts_xxzz_xzzzz, ts_xxzz_yyyy, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyz, ts_xxzz_yyyzz, ts_xxzz_yyzz, ts_xxzz_yyzzz, ts_xxzz_yzzz, ts_xxzz_yzzzz, ts_xxzz_zzzz, ts_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzz_xxxxx[i] = 4.0 * ts_xxz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxxy[i] = 4.0 * ts_xxz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxxz[i] = 4.0 * ts_xxz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxyy[i] = 4.0 * ts_xxz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxyz[i] = 4.0 * ts_xxz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxzz[i] = 4.0 * ts_xxz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxyyy[i] = 4.0 * ts_xxz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxyyz[i] = 4.0 * ts_xxz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxyzz[i] = 4.0 * ts_xxz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxzzz[i] = 4.0 * ts_xxz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyyyy[i] = 4.0 * ts_xxz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyyyz[i] = 4.0 * ts_xxz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyyzz[i] = 4.0 * ts_xxz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyzzz[i] = 4.0 * ts_xxz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xzzzz[i] = 4.0 * ts_xxz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyyyy[i] = 4.0 * ts_xxz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyyyz[i] = 4.0 * ts_xxz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyyzz[i] = 4.0 * ts_xxz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyzzz[i] = 4.0 * ts_xxz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yzzzz[i] = 4.0 * ts_xxz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_zzzzz[i] = 4.0 * ts_xxz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 756-777 components of targeted buffer : GH

    auto gs_z_xyyy_xxxxx = pbuffer.data(idx_g_gh + 756);

    auto gs_z_xyyy_xxxxy = pbuffer.data(idx_g_gh + 757);

    auto gs_z_xyyy_xxxxz = pbuffer.data(idx_g_gh + 758);

    auto gs_z_xyyy_xxxyy = pbuffer.data(idx_g_gh + 759);

    auto gs_z_xyyy_xxxyz = pbuffer.data(idx_g_gh + 760);

    auto gs_z_xyyy_xxxzz = pbuffer.data(idx_g_gh + 761);

    auto gs_z_xyyy_xxyyy = pbuffer.data(idx_g_gh + 762);

    auto gs_z_xyyy_xxyyz = pbuffer.data(idx_g_gh + 763);

    auto gs_z_xyyy_xxyzz = pbuffer.data(idx_g_gh + 764);

    auto gs_z_xyyy_xxzzz = pbuffer.data(idx_g_gh + 765);

    auto gs_z_xyyy_xyyyy = pbuffer.data(idx_g_gh + 766);

    auto gs_z_xyyy_xyyyz = pbuffer.data(idx_g_gh + 767);

    auto gs_z_xyyy_xyyzz = pbuffer.data(idx_g_gh + 768);

    auto gs_z_xyyy_xyzzz = pbuffer.data(idx_g_gh + 769);

    auto gs_z_xyyy_xzzzz = pbuffer.data(idx_g_gh + 770);

    auto gs_z_xyyy_yyyyy = pbuffer.data(idx_g_gh + 771);

    auto gs_z_xyyy_yyyyz = pbuffer.data(idx_g_gh + 772);

    auto gs_z_xyyy_yyyzz = pbuffer.data(idx_g_gh + 773);

    auto gs_z_xyyy_yyzzz = pbuffer.data(idx_g_gh + 774);

    auto gs_z_xyyy_yzzzz = pbuffer.data(idx_g_gh + 775);

    auto gs_z_xyyy_zzzzz = pbuffer.data(idx_g_gh + 776);

    #pragma omp simd aligned(gc_z, gs_z_xyyy_xxxxx, gs_z_xyyy_xxxxy, gs_z_xyyy_xxxxz, gs_z_xyyy_xxxyy, gs_z_xyyy_xxxyz, gs_z_xyyy_xxxzz, gs_z_xyyy_xxyyy, gs_z_xyyy_xxyyz, gs_z_xyyy_xxyzz, gs_z_xyyy_xxzzz, gs_z_xyyy_xyyyy, gs_z_xyyy_xyyyz, gs_z_xyyy_xyyzz, gs_z_xyyy_xyzzz, gs_z_xyyy_xzzzz, gs_z_xyyy_yyyyy, gs_z_xyyy_yyyyz, gs_z_xyyy_yyyzz, gs_z_xyyy_yyzzz, gs_z_xyyy_yzzzz, gs_z_xyyy_zzzzz, ts_xyyy_xxxx, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxy, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxz, ts_xyyy_xxxzz, ts_xyyy_xxyy, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyz, ts_xyyy_xxyzz, ts_xyyy_xxzz, ts_xyyy_xxzzz, ts_xyyy_xyyy, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyz, ts_xyyy_xyyzz, ts_xyyy_xyzz, ts_xyyy_xyzzz, ts_xyyy_xzzz, ts_xyyy_xzzzz, ts_xyyy_yyyy, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyz, ts_xyyy_yyyzz, ts_xyyy_yyzz, ts_xyyy_yyzzz, ts_xyyy_yzzz, ts_xyyy_yzzzz, ts_xyyy_zzzz, ts_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyy_xxxxx[i] = 2.0 * ts_xyyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxxy[i] = 2.0 * ts_xyyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxxz[i] = 2.0 * ts_xyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxyy[i] = 2.0 * ts_xyyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxyz[i] = 2.0 * ts_xyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxzz[i] = 4.0 * ts_xyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxyyy[i] = 2.0 * ts_xyyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxyyz[i] = 2.0 * ts_xyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxyzz[i] = 4.0 * ts_xyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxzzz[i] = 6.0 * ts_xyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyyyy[i] = 2.0 * ts_xyyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyyyz[i] = 2.0 * ts_xyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyyzz[i] = 4.0 * ts_xyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyzzz[i] = 6.0 * ts_xyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xzzzz[i] = 8.0 * ts_xyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyyyy[i] = 2.0 * ts_xyyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyyyz[i] = 2.0 * ts_xyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyyzz[i] = 4.0 * ts_xyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyzzz[i] = 6.0 * ts_xyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yzzzz[i] = 8.0 * ts_xyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_zzzzz[i] = 10.0 * ts_xyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 777-798 components of targeted buffer : GH

    auto gs_z_xyyz_xxxxx = pbuffer.data(idx_g_gh + 777);

    auto gs_z_xyyz_xxxxy = pbuffer.data(idx_g_gh + 778);

    auto gs_z_xyyz_xxxxz = pbuffer.data(idx_g_gh + 779);

    auto gs_z_xyyz_xxxyy = pbuffer.data(idx_g_gh + 780);

    auto gs_z_xyyz_xxxyz = pbuffer.data(idx_g_gh + 781);

    auto gs_z_xyyz_xxxzz = pbuffer.data(idx_g_gh + 782);

    auto gs_z_xyyz_xxyyy = pbuffer.data(idx_g_gh + 783);

    auto gs_z_xyyz_xxyyz = pbuffer.data(idx_g_gh + 784);

    auto gs_z_xyyz_xxyzz = pbuffer.data(idx_g_gh + 785);

    auto gs_z_xyyz_xxzzz = pbuffer.data(idx_g_gh + 786);

    auto gs_z_xyyz_xyyyy = pbuffer.data(idx_g_gh + 787);

    auto gs_z_xyyz_xyyyz = pbuffer.data(idx_g_gh + 788);

    auto gs_z_xyyz_xyyzz = pbuffer.data(idx_g_gh + 789);

    auto gs_z_xyyz_xyzzz = pbuffer.data(idx_g_gh + 790);

    auto gs_z_xyyz_xzzzz = pbuffer.data(idx_g_gh + 791);

    auto gs_z_xyyz_yyyyy = pbuffer.data(idx_g_gh + 792);

    auto gs_z_xyyz_yyyyz = pbuffer.data(idx_g_gh + 793);

    auto gs_z_xyyz_yyyzz = pbuffer.data(idx_g_gh + 794);

    auto gs_z_xyyz_yyzzz = pbuffer.data(idx_g_gh + 795);

    auto gs_z_xyyz_yzzzz = pbuffer.data(idx_g_gh + 796);

    auto gs_z_xyyz_zzzzz = pbuffer.data(idx_g_gh + 797);

    #pragma omp simd aligned(gc_z, gs_z_xyyz_xxxxx, gs_z_xyyz_xxxxy, gs_z_xyyz_xxxxz, gs_z_xyyz_xxxyy, gs_z_xyyz_xxxyz, gs_z_xyyz_xxxzz, gs_z_xyyz_xxyyy, gs_z_xyyz_xxyyz, gs_z_xyyz_xxyzz, gs_z_xyyz_xxzzz, gs_z_xyyz_xyyyy, gs_z_xyyz_xyyyz, gs_z_xyyz_xyyzz, gs_z_xyyz_xyzzz, gs_z_xyyz_xzzzz, gs_z_xyyz_yyyyy, gs_z_xyyz_yyyyz, gs_z_xyyz_yyyzz, gs_z_xyyz_yyzzz, gs_z_xyyz_yzzzz, gs_z_xyyz_zzzzz, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxzz, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyzz, ts_xyy_xxzzz, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyzz, ts_xyy_xyzzz, ts_xyy_xzzzz, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyzz, ts_xyy_yyzzz, ts_xyy_yzzzz, ts_xyy_zzzzz, ts_xyyz_xxxx, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxy, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxz, ts_xyyz_xxxzz, ts_xyyz_xxyy, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyz, ts_xyyz_xxyzz, ts_xyyz_xxzz, ts_xyyz_xxzzz, ts_xyyz_xyyy, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyz, ts_xyyz_xyyzz, ts_xyyz_xyzz, ts_xyyz_xyzzz, ts_xyyz_xzzz, ts_xyyz_xzzzz, ts_xyyz_yyyy, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyz, ts_xyyz_yyyzz, ts_xyyz_yyzz, ts_xyyz_yyzzz, ts_xyyz_yzzz, ts_xyyz_yzzzz, ts_xyyz_zzzz, ts_xyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyz_xxxxx[i] = 2.0 * ts_xyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxxy[i] = 2.0 * ts_xyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxxz[i] = 2.0 * ts_xyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxyy[i] = 2.0 * ts_xyy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxyz[i] = 2.0 * ts_xyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxzz[i] = 2.0 * ts_xyy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxyyy[i] = 2.0 * ts_xyy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxyyz[i] = 2.0 * ts_xyy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxyzz[i] = 2.0 * ts_xyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxzzz[i] = 2.0 * ts_xyy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyyyy[i] = 2.0 * ts_xyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyyyz[i] = 2.0 * ts_xyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyyzz[i] = 2.0 * ts_xyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyzzz[i] = 2.0 * ts_xyy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xzzzz[i] = 2.0 * ts_xyy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyyyy[i] = 2.0 * ts_xyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyyyz[i] = 2.0 * ts_xyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyyzz[i] = 2.0 * ts_xyy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyzzz[i] = 2.0 * ts_xyy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yzzzz[i] = 2.0 * ts_xyy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_zzzzz[i] = 2.0 * ts_xyy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 798-819 components of targeted buffer : GH

    auto gs_z_xyzz_xxxxx = pbuffer.data(idx_g_gh + 798);

    auto gs_z_xyzz_xxxxy = pbuffer.data(idx_g_gh + 799);

    auto gs_z_xyzz_xxxxz = pbuffer.data(idx_g_gh + 800);

    auto gs_z_xyzz_xxxyy = pbuffer.data(idx_g_gh + 801);

    auto gs_z_xyzz_xxxyz = pbuffer.data(idx_g_gh + 802);

    auto gs_z_xyzz_xxxzz = pbuffer.data(idx_g_gh + 803);

    auto gs_z_xyzz_xxyyy = pbuffer.data(idx_g_gh + 804);

    auto gs_z_xyzz_xxyyz = pbuffer.data(idx_g_gh + 805);

    auto gs_z_xyzz_xxyzz = pbuffer.data(idx_g_gh + 806);

    auto gs_z_xyzz_xxzzz = pbuffer.data(idx_g_gh + 807);

    auto gs_z_xyzz_xyyyy = pbuffer.data(idx_g_gh + 808);

    auto gs_z_xyzz_xyyyz = pbuffer.data(idx_g_gh + 809);

    auto gs_z_xyzz_xyyzz = pbuffer.data(idx_g_gh + 810);

    auto gs_z_xyzz_xyzzz = pbuffer.data(idx_g_gh + 811);

    auto gs_z_xyzz_xzzzz = pbuffer.data(idx_g_gh + 812);

    auto gs_z_xyzz_yyyyy = pbuffer.data(idx_g_gh + 813);

    auto gs_z_xyzz_yyyyz = pbuffer.data(idx_g_gh + 814);

    auto gs_z_xyzz_yyyzz = pbuffer.data(idx_g_gh + 815);

    auto gs_z_xyzz_yyzzz = pbuffer.data(idx_g_gh + 816);

    auto gs_z_xyzz_yzzzz = pbuffer.data(idx_g_gh + 817);

    auto gs_z_xyzz_zzzzz = pbuffer.data(idx_g_gh + 818);

    #pragma omp simd aligned(gc_z, gs_z_xyzz_xxxxx, gs_z_xyzz_xxxxy, gs_z_xyzz_xxxxz, gs_z_xyzz_xxxyy, gs_z_xyzz_xxxyz, gs_z_xyzz_xxxzz, gs_z_xyzz_xxyyy, gs_z_xyzz_xxyyz, gs_z_xyzz_xxyzz, gs_z_xyzz_xxzzz, gs_z_xyzz_xyyyy, gs_z_xyzz_xyyyz, gs_z_xyzz_xyyzz, gs_z_xyzz_xyzzz, gs_z_xyzz_xzzzz, gs_z_xyzz_yyyyy, gs_z_xyzz_yyyyz, gs_z_xyzz_yyyzz, gs_z_xyzz_yyzzz, gs_z_xyzz_yzzzz, gs_z_xyzz_zzzzz, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxzz, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyzz, ts_xyz_xxzzz, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyzz, ts_xyz_xyzzz, ts_xyz_xzzzz, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyzz, ts_xyz_yyzzz, ts_xyz_yzzzz, ts_xyz_zzzzz, ts_xyzz_xxxx, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxy, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxz, ts_xyzz_xxxzz, ts_xyzz_xxyy, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyz, ts_xyzz_xxyzz, ts_xyzz_xxzz, ts_xyzz_xxzzz, ts_xyzz_xyyy, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyz, ts_xyzz_xyyzz, ts_xyzz_xyzz, ts_xyzz_xyzzz, ts_xyzz_xzzz, ts_xyzz_xzzzz, ts_xyzz_yyyy, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyz, ts_xyzz_yyyzz, ts_xyzz_yyzz, ts_xyzz_yyzzz, ts_xyzz_yzzz, ts_xyzz_yzzzz, ts_xyzz_zzzz, ts_xyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzz_xxxxx[i] = 4.0 * ts_xyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxxy[i] = 4.0 * ts_xyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxxz[i] = 4.0 * ts_xyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxyy[i] = 4.0 * ts_xyz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxyz[i] = 4.0 * ts_xyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxzz[i] = 4.0 * ts_xyz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxyyy[i] = 4.0 * ts_xyz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxyyz[i] = 4.0 * ts_xyz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxyzz[i] = 4.0 * ts_xyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxzzz[i] = 4.0 * ts_xyz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyyyy[i] = 4.0 * ts_xyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyyyz[i] = 4.0 * ts_xyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyyzz[i] = 4.0 * ts_xyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyzzz[i] = 4.0 * ts_xyz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xzzzz[i] = 4.0 * ts_xyz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyyyy[i] = 4.0 * ts_xyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyyyz[i] = 4.0 * ts_xyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyyzz[i] = 4.0 * ts_xyz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyzzz[i] = 4.0 * ts_xyz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yzzzz[i] = 4.0 * ts_xyz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_zzzzz[i] = 4.0 * ts_xyz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 819-840 components of targeted buffer : GH

    auto gs_z_xzzz_xxxxx = pbuffer.data(idx_g_gh + 819);

    auto gs_z_xzzz_xxxxy = pbuffer.data(idx_g_gh + 820);

    auto gs_z_xzzz_xxxxz = pbuffer.data(idx_g_gh + 821);

    auto gs_z_xzzz_xxxyy = pbuffer.data(idx_g_gh + 822);

    auto gs_z_xzzz_xxxyz = pbuffer.data(idx_g_gh + 823);

    auto gs_z_xzzz_xxxzz = pbuffer.data(idx_g_gh + 824);

    auto gs_z_xzzz_xxyyy = pbuffer.data(idx_g_gh + 825);

    auto gs_z_xzzz_xxyyz = pbuffer.data(idx_g_gh + 826);

    auto gs_z_xzzz_xxyzz = pbuffer.data(idx_g_gh + 827);

    auto gs_z_xzzz_xxzzz = pbuffer.data(idx_g_gh + 828);

    auto gs_z_xzzz_xyyyy = pbuffer.data(idx_g_gh + 829);

    auto gs_z_xzzz_xyyyz = pbuffer.data(idx_g_gh + 830);

    auto gs_z_xzzz_xyyzz = pbuffer.data(idx_g_gh + 831);

    auto gs_z_xzzz_xyzzz = pbuffer.data(idx_g_gh + 832);

    auto gs_z_xzzz_xzzzz = pbuffer.data(idx_g_gh + 833);

    auto gs_z_xzzz_yyyyy = pbuffer.data(idx_g_gh + 834);

    auto gs_z_xzzz_yyyyz = pbuffer.data(idx_g_gh + 835);

    auto gs_z_xzzz_yyyzz = pbuffer.data(idx_g_gh + 836);

    auto gs_z_xzzz_yyzzz = pbuffer.data(idx_g_gh + 837);

    auto gs_z_xzzz_yzzzz = pbuffer.data(idx_g_gh + 838);

    auto gs_z_xzzz_zzzzz = pbuffer.data(idx_g_gh + 839);

    #pragma omp simd aligned(gc_z, gs_z_xzzz_xxxxx, gs_z_xzzz_xxxxy, gs_z_xzzz_xxxxz, gs_z_xzzz_xxxyy, gs_z_xzzz_xxxyz, gs_z_xzzz_xxxzz, gs_z_xzzz_xxyyy, gs_z_xzzz_xxyyz, gs_z_xzzz_xxyzz, gs_z_xzzz_xxzzz, gs_z_xzzz_xyyyy, gs_z_xzzz_xyyyz, gs_z_xzzz_xyyzz, gs_z_xzzz_xyzzz, gs_z_xzzz_xzzzz, gs_z_xzzz_yyyyy, gs_z_xzzz_yyyyz, gs_z_xzzz_yyyzz, gs_z_xzzz_yyzzz, gs_z_xzzz_yzzzz, gs_z_xzzz_zzzzz, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxzz, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyzz, ts_xzz_xxzzz, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyzz, ts_xzz_xyzzz, ts_xzz_xzzzz, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyzz, ts_xzz_yyzzz, ts_xzz_yzzzz, ts_xzz_zzzzz, ts_xzzz_xxxx, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxy, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxz, ts_xzzz_xxxzz, ts_xzzz_xxyy, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyz, ts_xzzz_xxyzz, ts_xzzz_xxzz, ts_xzzz_xxzzz, ts_xzzz_xyyy, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyz, ts_xzzz_xyyzz, ts_xzzz_xyzz, ts_xzzz_xyzzz, ts_xzzz_xzzz, ts_xzzz_xzzzz, ts_xzzz_yyyy, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyz, ts_xzzz_yyyzz, ts_xzzz_yyzz, ts_xzzz_yyzzz, ts_xzzz_yzzz, ts_xzzz_yzzzz, ts_xzzz_zzzz, ts_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzz_xxxxx[i] = 6.0 * ts_xzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxxy[i] = 6.0 * ts_xzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxxz[i] = 6.0 * ts_xzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxyy[i] = 6.0 * ts_xzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxyz[i] = 6.0 * ts_xzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxzz[i] = 6.0 * ts_xzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxyyy[i] = 6.0 * ts_xzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxyyz[i] = 6.0 * ts_xzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxyzz[i] = 6.0 * ts_xzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxzzz[i] = 6.0 * ts_xzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyyyy[i] = 6.0 * ts_xzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyyyz[i] = 6.0 * ts_xzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyyzz[i] = 6.0 * ts_xzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyzzz[i] = 6.0 * ts_xzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xzzzz[i] = 6.0 * ts_xzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyyyy[i] = 6.0 * ts_xzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyyyz[i] = 6.0 * ts_xzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyyzz[i] = 6.0 * ts_xzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyzzz[i] = 6.0 * ts_xzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yzzzz[i] = 6.0 * ts_xzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_zzzzz[i] = 6.0 * ts_xzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 840-861 components of targeted buffer : GH

    auto gs_z_yyyy_xxxxx = pbuffer.data(idx_g_gh + 840);

    auto gs_z_yyyy_xxxxy = pbuffer.data(idx_g_gh + 841);

    auto gs_z_yyyy_xxxxz = pbuffer.data(idx_g_gh + 842);

    auto gs_z_yyyy_xxxyy = pbuffer.data(idx_g_gh + 843);

    auto gs_z_yyyy_xxxyz = pbuffer.data(idx_g_gh + 844);

    auto gs_z_yyyy_xxxzz = pbuffer.data(idx_g_gh + 845);

    auto gs_z_yyyy_xxyyy = pbuffer.data(idx_g_gh + 846);

    auto gs_z_yyyy_xxyyz = pbuffer.data(idx_g_gh + 847);

    auto gs_z_yyyy_xxyzz = pbuffer.data(idx_g_gh + 848);

    auto gs_z_yyyy_xxzzz = pbuffer.data(idx_g_gh + 849);

    auto gs_z_yyyy_xyyyy = pbuffer.data(idx_g_gh + 850);

    auto gs_z_yyyy_xyyyz = pbuffer.data(idx_g_gh + 851);

    auto gs_z_yyyy_xyyzz = pbuffer.data(idx_g_gh + 852);

    auto gs_z_yyyy_xyzzz = pbuffer.data(idx_g_gh + 853);

    auto gs_z_yyyy_xzzzz = pbuffer.data(idx_g_gh + 854);

    auto gs_z_yyyy_yyyyy = pbuffer.data(idx_g_gh + 855);

    auto gs_z_yyyy_yyyyz = pbuffer.data(idx_g_gh + 856);

    auto gs_z_yyyy_yyyzz = pbuffer.data(idx_g_gh + 857);

    auto gs_z_yyyy_yyzzz = pbuffer.data(idx_g_gh + 858);

    auto gs_z_yyyy_yzzzz = pbuffer.data(idx_g_gh + 859);

    auto gs_z_yyyy_zzzzz = pbuffer.data(idx_g_gh + 860);

    #pragma omp simd aligned(gc_z, gs_z_yyyy_xxxxx, gs_z_yyyy_xxxxy, gs_z_yyyy_xxxxz, gs_z_yyyy_xxxyy, gs_z_yyyy_xxxyz, gs_z_yyyy_xxxzz, gs_z_yyyy_xxyyy, gs_z_yyyy_xxyyz, gs_z_yyyy_xxyzz, gs_z_yyyy_xxzzz, gs_z_yyyy_xyyyy, gs_z_yyyy_xyyyz, gs_z_yyyy_xyyzz, gs_z_yyyy_xyzzz, gs_z_yyyy_xzzzz, gs_z_yyyy_yyyyy, gs_z_yyyy_yyyyz, gs_z_yyyy_yyyzz, gs_z_yyyy_yyzzz, gs_z_yyyy_yzzzz, gs_z_yyyy_zzzzz, ts_yyyy_xxxx, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxz, ts_yyyy_xxxzz, ts_yyyy_xxyy, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyz, ts_yyyy_xxyzz, ts_yyyy_xxzz, ts_yyyy_xxzzz, ts_yyyy_xyyy, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyz, ts_yyyy_xyyzz, ts_yyyy_xyzz, ts_yyyy_xyzzz, ts_yyyy_xzzz, ts_yyyy_xzzzz, ts_yyyy_yyyy, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyz, ts_yyyy_yyyzz, ts_yyyy_yyzz, ts_yyyy_yyzzz, ts_yyyy_yzzz, ts_yyyy_yzzzz, ts_yyyy_zzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyy_xxxxx[i] = 2.0 * ts_yyyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxxy[i] = 2.0 * ts_yyyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxxz[i] = 2.0 * ts_yyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxyy[i] = 2.0 * ts_yyyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxyz[i] = 2.0 * ts_yyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxzz[i] = 4.0 * ts_yyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxyyy[i] = 2.0 * ts_yyyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxyyz[i] = 2.0 * ts_yyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxyzz[i] = 4.0 * ts_yyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxzzz[i] = 6.0 * ts_yyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyyyy[i] = 2.0 * ts_yyyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyyyz[i] = 2.0 * ts_yyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyyzz[i] = 4.0 * ts_yyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyzzz[i] = 6.0 * ts_yyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xzzzz[i] = 8.0 * ts_yyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyyyy[i] = 2.0 * ts_yyyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyyyz[i] = 2.0 * ts_yyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyyzz[i] = 4.0 * ts_yyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyzzz[i] = 6.0 * ts_yyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yzzzz[i] = 8.0 * ts_yyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_zzzzz[i] = 10.0 * ts_yyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 861-882 components of targeted buffer : GH

    auto gs_z_yyyz_xxxxx = pbuffer.data(idx_g_gh + 861);

    auto gs_z_yyyz_xxxxy = pbuffer.data(idx_g_gh + 862);

    auto gs_z_yyyz_xxxxz = pbuffer.data(idx_g_gh + 863);

    auto gs_z_yyyz_xxxyy = pbuffer.data(idx_g_gh + 864);

    auto gs_z_yyyz_xxxyz = pbuffer.data(idx_g_gh + 865);

    auto gs_z_yyyz_xxxzz = pbuffer.data(idx_g_gh + 866);

    auto gs_z_yyyz_xxyyy = pbuffer.data(idx_g_gh + 867);

    auto gs_z_yyyz_xxyyz = pbuffer.data(idx_g_gh + 868);

    auto gs_z_yyyz_xxyzz = pbuffer.data(idx_g_gh + 869);

    auto gs_z_yyyz_xxzzz = pbuffer.data(idx_g_gh + 870);

    auto gs_z_yyyz_xyyyy = pbuffer.data(idx_g_gh + 871);

    auto gs_z_yyyz_xyyyz = pbuffer.data(idx_g_gh + 872);

    auto gs_z_yyyz_xyyzz = pbuffer.data(idx_g_gh + 873);

    auto gs_z_yyyz_xyzzz = pbuffer.data(idx_g_gh + 874);

    auto gs_z_yyyz_xzzzz = pbuffer.data(idx_g_gh + 875);

    auto gs_z_yyyz_yyyyy = pbuffer.data(idx_g_gh + 876);

    auto gs_z_yyyz_yyyyz = pbuffer.data(idx_g_gh + 877);

    auto gs_z_yyyz_yyyzz = pbuffer.data(idx_g_gh + 878);

    auto gs_z_yyyz_yyzzz = pbuffer.data(idx_g_gh + 879);

    auto gs_z_yyyz_yzzzz = pbuffer.data(idx_g_gh + 880);

    auto gs_z_yyyz_zzzzz = pbuffer.data(idx_g_gh + 881);

    #pragma omp simd aligned(gc_z, gs_z_yyyz_xxxxx, gs_z_yyyz_xxxxy, gs_z_yyyz_xxxxz, gs_z_yyyz_xxxyy, gs_z_yyyz_xxxyz, gs_z_yyyz_xxxzz, gs_z_yyyz_xxyyy, gs_z_yyyz_xxyyz, gs_z_yyyz_xxyzz, gs_z_yyyz_xxzzz, gs_z_yyyz_xyyyy, gs_z_yyyz_xyyyz, gs_z_yyyz_xyyzz, gs_z_yyyz_xyzzz, gs_z_yyyz_xzzzz, gs_z_yyyz_yyyyy, gs_z_yyyz_yyyyz, gs_z_yyyz_yyyzz, gs_z_yyyz_yyzzz, gs_z_yyyz_yzzzz, gs_z_yyyz_zzzzz, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxzz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xxzzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_xzzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, ts_yyyz_xxxx, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxy, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxz, ts_yyyz_xxxzz, ts_yyyz_xxyy, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyz, ts_yyyz_xxyzz, ts_yyyz_xxzz, ts_yyyz_xxzzz, ts_yyyz_xyyy, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyz, ts_yyyz_xyyzz, ts_yyyz_xyzz, ts_yyyz_xyzzz, ts_yyyz_xzzz, ts_yyyz_xzzzz, ts_yyyz_yyyy, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyz, ts_yyyz_yyyzz, ts_yyyz_yyzz, ts_yyyz_yyzzz, ts_yyyz_yzzz, ts_yyyz_yzzzz, ts_yyyz_zzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyz_xxxxx[i] = 2.0 * ts_yyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxxy[i] = 2.0 * ts_yyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxxz[i] = 2.0 * ts_yyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxyy[i] = 2.0 * ts_yyy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxyz[i] = 2.0 * ts_yyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxzz[i] = 2.0 * ts_yyy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxyyy[i] = 2.0 * ts_yyy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxyyz[i] = 2.0 * ts_yyy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxyzz[i] = 2.0 * ts_yyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxzzz[i] = 2.0 * ts_yyy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyyyy[i] = 2.0 * ts_yyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyyyz[i] = 2.0 * ts_yyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyyzz[i] = 2.0 * ts_yyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyzzz[i] = 2.0 * ts_yyy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xzzzz[i] = 2.0 * ts_yyy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyyyy[i] = 2.0 * ts_yyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyyyz[i] = 2.0 * ts_yyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyyzz[i] = 2.0 * ts_yyy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyzzz[i] = 2.0 * ts_yyy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yzzzz[i] = 2.0 * ts_yyy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_zzzzz[i] = 2.0 * ts_yyy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 882-903 components of targeted buffer : GH

    auto gs_z_yyzz_xxxxx = pbuffer.data(idx_g_gh + 882);

    auto gs_z_yyzz_xxxxy = pbuffer.data(idx_g_gh + 883);

    auto gs_z_yyzz_xxxxz = pbuffer.data(idx_g_gh + 884);

    auto gs_z_yyzz_xxxyy = pbuffer.data(idx_g_gh + 885);

    auto gs_z_yyzz_xxxyz = pbuffer.data(idx_g_gh + 886);

    auto gs_z_yyzz_xxxzz = pbuffer.data(idx_g_gh + 887);

    auto gs_z_yyzz_xxyyy = pbuffer.data(idx_g_gh + 888);

    auto gs_z_yyzz_xxyyz = pbuffer.data(idx_g_gh + 889);

    auto gs_z_yyzz_xxyzz = pbuffer.data(idx_g_gh + 890);

    auto gs_z_yyzz_xxzzz = pbuffer.data(idx_g_gh + 891);

    auto gs_z_yyzz_xyyyy = pbuffer.data(idx_g_gh + 892);

    auto gs_z_yyzz_xyyyz = pbuffer.data(idx_g_gh + 893);

    auto gs_z_yyzz_xyyzz = pbuffer.data(idx_g_gh + 894);

    auto gs_z_yyzz_xyzzz = pbuffer.data(idx_g_gh + 895);

    auto gs_z_yyzz_xzzzz = pbuffer.data(idx_g_gh + 896);

    auto gs_z_yyzz_yyyyy = pbuffer.data(idx_g_gh + 897);

    auto gs_z_yyzz_yyyyz = pbuffer.data(idx_g_gh + 898);

    auto gs_z_yyzz_yyyzz = pbuffer.data(idx_g_gh + 899);

    auto gs_z_yyzz_yyzzz = pbuffer.data(idx_g_gh + 900);

    auto gs_z_yyzz_yzzzz = pbuffer.data(idx_g_gh + 901);

    auto gs_z_yyzz_zzzzz = pbuffer.data(idx_g_gh + 902);

    #pragma omp simd aligned(gc_z, gs_z_yyzz_xxxxx, gs_z_yyzz_xxxxy, gs_z_yyzz_xxxxz, gs_z_yyzz_xxxyy, gs_z_yyzz_xxxyz, gs_z_yyzz_xxxzz, gs_z_yyzz_xxyyy, gs_z_yyzz_xxyyz, gs_z_yyzz_xxyzz, gs_z_yyzz_xxzzz, gs_z_yyzz_xyyyy, gs_z_yyzz_xyyyz, gs_z_yyzz_xyyzz, gs_z_yyzz_xyzzz, gs_z_yyzz_xzzzz, gs_z_yyzz_yyyyy, gs_z_yyzz_yyyyz, gs_z_yyzz_yyyzz, gs_z_yyzz_yyzzz, gs_z_yyzz_yzzzz, gs_z_yyzz_zzzzz, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxzz, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyzz, ts_yyz_xxzzz, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyzz, ts_yyz_xyzzz, ts_yyz_xzzzz, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyzz, ts_yyz_yyzzz, ts_yyz_yzzzz, ts_yyz_zzzzz, ts_yyzz_xxxx, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxy, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxz, ts_yyzz_xxxzz, ts_yyzz_xxyy, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyz, ts_yyzz_xxyzz, ts_yyzz_xxzz, ts_yyzz_xxzzz, ts_yyzz_xyyy, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyz, ts_yyzz_xyyzz, ts_yyzz_xyzz, ts_yyzz_xyzzz, ts_yyzz_xzzz, ts_yyzz_xzzzz, ts_yyzz_yyyy, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyz, ts_yyzz_yyyzz, ts_yyzz_yyzz, ts_yyzz_yyzzz, ts_yyzz_yzzz, ts_yyzz_yzzzz, ts_yyzz_zzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzz_xxxxx[i] = 4.0 * ts_yyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxxy[i] = 4.0 * ts_yyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxxz[i] = 4.0 * ts_yyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxyy[i] = 4.0 * ts_yyz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxyz[i] = 4.0 * ts_yyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxzz[i] = 4.0 * ts_yyz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxyyy[i] = 4.0 * ts_yyz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxyyz[i] = 4.0 * ts_yyz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxyzz[i] = 4.0 * ts_yyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxzzz[i] = 4.0 * ts_yyz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyyyy[i] = 4.0 * ts_yyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyyyz[i] = 4.0 * ts_yyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyyzz[i] = 4.0 * ts_yyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyzzz[i] = 4.0 * ts_yyz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xzzzz[i] = 4.0 * ts_yyz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyyyy[i] = 4.0 * ts_yyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyyyz[i] = 4.0 * ts_yyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyyzz[i] = 4.0 * ts_yyz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyzzz[i] = 4.0 * ts_yyz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yzzzz[i] = 4.0 * ts_yyz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_zzzzz[i] = 4.0 * ts_yyz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 903-924 components of targeted buffer : GH

    auto gs_z_yzzz_xxxxx = pbuffer.data(idx_g_gh + 903);

    auto gs_z_yzzz_xxxxy = pbuffer.data(idx_g_gh + 904);

    auto gs_z_yzzz_xxxxz = pbuffer.data(idx_g_gh + 905);

    auto gs_z_yzzz_xxxyy = pbuffer.data(idx_g_gh + 906);

    auto gs_z_yzzz_xxxyz = pbuffer.data(idx_g_gh + 907);

    auto gs_z_yzzz_xxxzz = pbuffer.data(idx_g_gh + 908);

    auto gs_z_yzzz_xxyyy = pbuffer.data(idx_g_gh + 909);

    auto gs_z_yzzz_xxyyz = pbuffer.data(idx_g_gh + 910);

    auto gs_z_yzzz_xxyzz = pbuffer.data(idx_g_gh + 911);

    auto gs_z_yzzz_xxzzz = pbuffer.data(idx_g_gh + 912);

    auto gs_z_yzzz_xyyyy = pbuffer.data(idx_g_gh + 913);

    auto gs_z_yzzz_xyyyz = pbuffer.data(idx_g_gh + 914);

    auto gs_z_yzzz_xyyzz = pbuffer.data(idx_g_gh + 915);

    auto gs_z_yzzz_xyzzz = pbuffer.data(idx_g_gh + 916);

    auto gs_z_yzzz_xzzzz = pbuffer.data(idx_g_gh + 917);

    auto gs_z_yzzz_yyyyy = pbuffer.data(idx_g_gh + 918);

    auto gs_z_yzzz_yyyyz = pbuffer.data(idx_g_gh + 919);

    auto gs_z_yzzz_yyyzz = pbuffer.data(idx_g_gh + 920);

    auto gs_z_yzzz_yyzzz = pbuffer.data(idx_g_gh + 921);

    auto gs_z_yzzz_yzzzz = pbuffer.data(idx_g_gh + 922);

    auto gs_z_yzzz_zzzzz = pbuffer.data(idx_g_gh + 923);

    #pragma omp simd aligned(gc_z, gs_z_yzzz_xxxxx, gs_z_yzzz_xxxxy, gs_z_yzzz_xxxxz, gs_z_yzzz_xxxyy, gs_z_yzzz_xxxyz, gs_z_yzzz_xxxzz, gs_z_yzzz_xxyyy, gs_z_yzzz_xxyyz, gs_z_yzzz_xxyzz, gs_z_yzzz_xxzzz, gs_z_yzzz_xyyyy, gs_z_yzzz_xyyyz, gs_z_yzzz_xyyzz, gs_z_yzzz_xyzzz, gs_z_yzzz_xzzzz, gs_z_yzzz_yyyyy, gs_z_yzzz_yyyyz, gs_z_yzzz_yyyzz, gs_z_yzzz_yyzzz, gs_z_yzzz_yzzzz, gs_z_yzzz_zzzzz, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxzz, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyzz, ts_yzz_xxzzz, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyzz, ts_yzz_xyzzz, ts_yzz_xzzzz, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, ts_yzz_zzzzz, ts_yzzz_xxxx, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxy, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxz, ts_yzzz_xxxzz, ts_yzzz_xxyy, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyz, ts_yzzz_xxyzz, ts_yzzz_xxzz, ts_yzzz_xxzzz, ts_yzzz_xyyy, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyz, ts_yzzz_xyyzz, ts_yzzz_xyzz, ts_yzzz_xyzzz, ts_yzzz_xzzz, ts_yzzz_xzzzz, ts_yzzz_yyyy, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyz, ts_yzzz_yyyzz, ts_yzzz_yyzz, ts_yzzz_yyzzz, ts_yzzz_yzzz, ts_yzzz_yzzzz, ts_yzzz_zzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzz_xxxxx[i] = 6.0 * ts_yzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxxy[i] = 6.0 * ts_yzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxxz[i] = 6.0 * ts_yzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxyy[i] = 6.0 * ts_yzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxyz[i] = 6.0 * ts_yzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxzz[i] = 6.0 * ts_yzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxyyy[i] = 6.0 * ts_yzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxyyz[i] = 6.0 * ts_yzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxyzz[i] = 6.0 * ts_yzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxzzz[i] = 6.0 * ts_yzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyyyy[i] = 6.0 * ts_yzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyyyz[i] = 6.0 * ts_yzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyyzz[i] = 6.0 * ts_yzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyzzz[i] = 6.0 * ts_yzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xzzzz[i] = 6.0 * ts_yzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyyyy[i] = 6.0 * ts_yzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyyyz[i] = 6.0 * ts_yzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyyzz[i] = 6.0 * ts_yzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyzzz[i] = 6.0 * ts_yzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yzzzz[i] = 6.0 * ts_yzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_zzzzz[i] = 6.0 * ts_yzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 924-945 components of targeted buffer : GH

    auto gs_z_zzzz_xxxxx = pbuffer.data(idx_g_gh + 924);

    auto gs_z_zzzz_xxxxy = pbuffer.data(idx_g_gh + 925);

    auto gs_z_zzzz_xxxxz = pbuffer.data(idx_g_gh + 926);

    auto gs_z_zzzz_xxxyy = pbuffer.data(idx_g_gh + 927);

    auto gs_z_zzzz_xxxyz = pbuffer.data(idx_g_gh + 928);

    auto gs_z_zzzz_xxxzz = pbuffer.data(idx_g_gh + 929);

    auto gs_z_zzzz_xxyyy = pbuffer.data(idx_g_gh + 930);

    auto gs_z_zzzz_xxyyz = pbuffer.data(idx_g_gh + 931);

    auto gs_z_zzzz_xxyzz = pbuffer.data(idx_g_gh + 932);

    auto gs_z_zzzz_xxzzz = pbuffer.data(idx_g_gh + 933);

    auto gs_z_zzzz_xyyyy = pbuffer.data(idx_g_gh + 934);

    auto gs_z_zzzz_xyyyz = pbuffer.data(idx_g_gh + 935);

    auto gs_z_zzzz_xyyzz = pbuffer.data(idx_g_gh + 936);

    auto gs_z_zzzz_xyzzz = pbuffer.data(idx_g_gh + 937);

    auto gs_z_zzzz_xzzzz = pbuffer.data(idx_g_gh + 938);

    auto gs_z_zzzz_yyyyy = pbuffer.data(idx_g_gh + 939);

    auto gs_z_zzzz_yyyyz = pbuffer.data(idx_g_gh + 940);

    auto gs_z_zzzz_yyyzz = pbuffer.data(idx_g_gh + 941);

    auto gs_z_zzzz_yyzzz = pbuffer.data(idx_g_gh + 942);

    auto gs_z_zzzz_yzzzz = pbuffer.data(idx_g_gh + 943);

    auto gs_z_zzzz_zzzzz = pbuffer.data(idx_g_gh + 944);

    #pragma omp simd aligned(gc_z, gs_z_zzzz_xxxxx, gs_z_zzzz_xxxxy, gs_z_zzzz_xxxxz, gs_z_zzzz_xxxyy, gs_z_zzzz_xxxyz, gs_z_zzzz_xxxzz, gs_z_zzzz_xxyyy, gs_z_zzzz_xxyyz, gs_z_zzzz_xxyzz, gs_z_zzzz_xxzzz, gs_z_zzzz_xyyyy, gs_z_zzzz_xyyyz, gs_z_zzzz_xyyzz, gs_z_zzzz_xyzzz, gs_z_zzzz_xzzzz, gs_z_zzzz_yyyyy, gs_z_zzzz_yyyyz, gs_z_zzzz_yyyzz, gs_z_zzzz_yyzzz, gs_z_zzzz_yzzzz, gs_z_zzzz_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, ts_zzzz_xxxx, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxy, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxz, ts_zzzz_xxxzz, ts_zzzz_xxyy, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyz, ts_zzzz_xxyzz, ts_zzzz_xxzz, ts_zzzz_xxzzz, ts_zzzz_xyyy, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyz, ts_zzzz_xyyzz, ts_zzzz_xyzz, ts_zzzz_xyzzz, ts_zzzz_xzzz, ts_zzzz_xzzzz, ts_zzzz_yyyy, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyz, ts_zzzz_yyyzz, ts_zzzz_yyzz, ts_zzzz_yyzzz, ts_zzzz_yzzz, ts_zzzz_yzzzz, ts_zzzz_zzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzz_xxxxx[i] = 8.0 * ts_zzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxxy[i] = 8.0 * ts_zzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxxz[i] = 8.0 * ts_zzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxyy[i] = 8.0 * ts_zzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxyz[i] = 8.0 * ts_zzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxzz[i] = 8.0 * ts_zzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxyyy[i] = 8.0 * ts_zzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxyyz[i] = 8.0 * ts_zzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxyzz[i] = 8.0 * ts_zzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxzzz[i] = 8.0 * ts_zzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyyyy[i] = 8.0 * ts_zzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyyyz[i] = 8.0 * ts_zzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyyzz[i] = 8.0 * ts_zzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyzzz[i] = 8.0 * ts_zzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xzzzz[i] = 8.0 * ts_zzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyyyy[i] = 8.0 * ts_zzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyyyz[i] = 8.0 * ts_zzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyyzz[i] = 8.0 * ts_zzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyzzz[i] = 8.0 * ts_zzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yzzzz[i] = 8.0 * ts_zzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_zzzzz[i] = 8.0 * ts_zzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_zzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

