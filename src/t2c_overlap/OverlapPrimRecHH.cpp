#include "OverlapPrimRecHH.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_hh(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_hh,
                     const size_t idx_ovl_fh,
                     const size_t idx_ovl_gg,
                     const size_t idx_ovl_gh,
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

    // Set up components of auxiliary buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_ovl_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_ovl_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_ovl_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_ovl_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_ovl_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_ovl_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_ovl_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_ovl_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_ovl_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_ovl_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_ovl_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_ovl_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_ovl_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_ovl_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_ovl_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_ovl_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_ovl_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_ovl_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_ovl_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_ovl_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_ovl_fh + 20);

    auto ts_xxy_xxxxx = pbuffer.data(idx_ovl_fh + 21);

    auto ts_xxy_xxxxz = pbuffer.data(idx_ovl_fh + 23);

    auto ts_xxy_xxxzz = pbuffer.data(idx_ovl_fh + 26);

    auto ts_xxy_xxzzz = pbuffer.data(idx_ovl_fh + 30);

    auto ts_xxy_xzzzz = pbuffer.data(idx_ovl_fh + 35);

    auto ts_xxy_yyyyy = pbuffer.data(idx_ovl_fh + 36);

    auto ts_xxy_yyyyz = pbuffer.data(idx_ovl_fh + 37);

    auto ts_xxy_yyyzz = pbuffer.data(idx_ovl_fh + 38);

    auto ts_xxy_yyzzz = pbuffer.data(idx_ovl_fh + 39);

    auto ts_xxy_yzzzz = pbuffer.data(idx_ovl_fh + 40);

    auto ts_xxz_xxxxx = pbuffer.data(idx_ovl_fh + 42);

    auto ts_xxz_xxxxy = pbuffer.data(idx_ovl_fh + 43);

    auto ts_xxz_xxxxz = pbuffer.data(idx_ovl_fh + 44);

    auto ts_xxz_xxxyy = pbuffer.data(idx_ovl_fh + 45);

    auto ts_xxz_xxxzz = pbuffer.data(idx_ovl_fh + 47);

    auto ts_xxz_xxyyy = pbuffer.data(idx_ovl_fh + 48);

    auto ts_xxz_xxzzz = pbuffer.data(idx_ovl_fh + 51);

    auto ts_xxz_xyyyy = pbuffer.data(idx_ovl_fh + 52);

    auto ts_xxz_xzzzz = pbuffer.data(idx_ovl_fh + 56);

    auto ts_xxz_yyyyz = pbuffer.data(idx_ovl_fh + 58);

    auto ts_xxz_yyyzz = pbuffer.data(idx_ovl_fh + 59);

    auto ts_xxz_yyzzz = pbuffer.data(idx_ovl_fh + 60);

    auto ts_xxz_yzzzz = pbuffer.data(idx_ovl_fh + 61);

    auto ts_xxz_zzzzz = pbuffer.data(idx_ovl_fh + 62);

    auto ts_xyy_xxxxy = pbuffer.data(idx_ovl_fh + 64);

    auto ts_xyy_xxxyy = pbuffer.data(idx_ovl_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_ovl_fh + 67);

    auto ts_xyy_xxyyy = pbuffer.data(idx_ovl_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_ovl_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_ovl_fh + 71);

    auto ts_xyy_xyyyy = pbuffer.data(idx_ovl_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_ovl_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_ovl_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_ovl_fh + 76);

    auto ts_xyy_yyyyy = pbuffer.data(idx_ovl_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_ovl_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_ovl_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_ovl_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_ovl_fh + 82);

    auto ts_xyy_zzzzz = pbuffer.data(idx_ovl_fh + 83);

    auto ts_xyz_yyyyz = pbuffer.data(idx_ovl_fh + 100);

    auto ts_xyz_yyyzz = pbuffer.data(idx_ovl_fh + 101);

    auto ts_xyz_yyzzz = pbuffer.data(idx_ovl_fh + 102);

    auto ts_xyz_yzzzz = pbuffer.data(idx_ovl_fh + 103);

    auto ts_xzz_xxxxz = pbuffer.data(idx_ovl_fh + 107);

    auto ts_xzz_xxxyz = pbuffer.data(idx_ovl_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_ovl_fh + 110);

    auto ts_xzz_xxyyz = pbuffer.data(idx_ovl_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_ovl_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_ovl_fh + 114);

    auto ts_xzz_xyyyz = pbuffer.data(idx_ovl_fh + 116);

    auto ts_xzz_xyyzz = pbuffer.data(idx_ovl_fh + 117);

    auto ts_xzz_xyzzz = pbuffer.data(idx_ovl_fh + 118);

    auto ts_xzz_xzzzz = pbuffer.data(idx_ovl_fh + 119);

    auto ts_xzz_yyyyy = pbuffer.data(idx_ovl_fh + 120);

    auto ts_xzz_yyyyz = pbuffer.data(idx_ovl_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_ovl_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_ovl_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_ovl_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_ovl_fh + 125);

    auto ts_yyy_xxxxx = pbuffer.data(idx_ovl_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_ovl_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_ovl_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_ovl_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_ovl_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_ovl_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_ovl_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_ovl_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_ovl_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_ovl_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_ovl_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_ovl_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_ovl_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_ovl_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_ovl_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_ovl_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_ovl_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_ovl_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_ovl_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_ovl_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_ovl_fh + 146);

    auto ts_yyz_xxxxy = pbuffer.data(idx_ovl_fh + 148);

    auto ts_yyz_xxxxz = pbuffer.data(idx_ovl_fh + 149);

    auto ts_yyz_xxxyy = pbuffer.data(idx_ovl_fh + 150);

    auto ts_yyz_xxxzz = pbuffer.data(idx_ovl_fh + 152);

    auto ts_yyz_xxyyy = pbuffer.data(idx_ovl_fh + 153);

    auto ts_yyz_xxzzz = pbuffer.data(idx_ovl_fh + 156);

    auto ts_yyz_xyyyy = pbuffer.data(idx_ovl_fh + 157);

    auto ts_yyz_xzzzz = pbuffer.data(idx_ovl_fh + 161);

    auto ts_yyz_yyyyy = pbuffer.data(idx_ovl_fh + 162);

    auto ts_yyz_yyyyz = pbuffer.data(idx_ovl_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_ovl_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_ovl_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_ovl_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_ovl_fh + 167);

    auto ts_yzz_xxxxx = pbuffer.data(idx_ovl_fh + 168);

    auto ts_yzz_xxxxz = pbuffer.data(idx_ovl_fh + 170);

    auto ts_yzz_xxxyz = pbuffer.data(idx_ovl_fh + 172);

    auto ts_yzz_xxxzz = pbuffer.data(idx_ovl_fh + 173);

    auto ts_yzz_xxyyz = pbuffer.data(idx_ovl_fh + 175);

    auto ts_yzz_xxyzz = pbuffer.data(idx_ovl_fh + 176);

    auto ts_yzz_xxzzz = pbuffer.data(idx_ovl_fh + 177);

    auto ts_yzz_xyyyz = pbuffer.data(idx_ovl_fh + 179);

    auto ts_yzz_xyyzz = pbuffer.data(idx_ovl_fh + 180);

    auto ts_yzz_xyzzz = pbuffer.data(idx_ovl_fh + 181);

    auto ts_yzz_xzzzz = pbuffer.data(idx_ovl_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_ovl_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_ovl_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_ovl_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_ovl_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_ovl_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_ovl_fh + 188);

    auto ts_zzz_xxxxx = pbuffer.data(idx_ovl_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_ovl_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_ovl_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_ovl_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_ovl_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_ovl_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_ovl_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_ovl_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_ovl_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_ovl_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_ovl_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_ovl_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_ovl_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_ovl_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_ovl_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_ovl_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_ovl_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_ovl_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_ovl_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_ovl_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_ovl_fh + 209);

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

    auto ts_xxxz_xxxz = pbuffer.data(idx_ovl_gg + 32);

    auto ts_xxxz_xxyz = pbuffer.data(idx_ovl_gg + 34);

    auto ts_xxxz_xxzz = pbuffer.data(idx_ovl_gg + 35);

    auto ts_xxxz_xyyz = pbuffer.data(idx_ovl_gg + 37);

    auto ts_xxxz_xyzz = pbuffer.data(idx_ovl_gg + 38);

    auto ts_xxxz_xzzz = pbuffer.data(idx_ovl_gg + 39);

    auto ts_xxyy_xxxy = pbuffer.data(idx_ovl_gg + 46);

    auto ts_xxyy_xxyy = pbuffer.data(idx_ovl_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_ovl_gg + 49);

    auto ts_xxyy_xyyy = pbuffer.data(idx_ovl_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_ovl_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_ovl_gg + 53);

    auto ts_xxyy_yyyy = pbuffer.data(idx_ovl_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_ovl_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_ovl_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_ovl_gg + 58);

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

    auto ts_xzzz_xxxz = pbuffer.data(idx_ovl_gg + 137);

    auto ts_xzzz_xxyz = pbuffer.data(idx_ovl_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_ovl_gg + 140);

    auto ts_xzzz_xyyz = pbuffer.data(idx_ovl_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_ovl_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_ovl_gg + 144);

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

    auto ts_yyyz_xxxz = pbuffer.data(idx_ovl_gg + 167);

    auto ts_yyyz_xxyz = pbuffer.data(idx_ovl_gg + 169);

    auto ts_yyyz_xxzz = pbuffer.data(idx_ovl_gg + 170);

    auto ts_yyyz_xyyz = pbuffer.data(idx_ovl_gg + 172);

    auto ts_yyyz_xyzz = pbuffer.data(idx_ovl_gg + 173);

    auto ts_yyyz_xzzz = pbuffer.data(idx_ovl_gg + 174);

    auto ts_yyyz_yyyz = pbuffer.data(idx_ovl_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_ovl_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_ovl_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_ovl_gg + 179);

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

    auto ts_yzzz_xxxy = pbuffer.data(idx_ovl_gg + 196);

    auto ts_yzzz_xxxz = pbuffer.data(idx_ovl_gg + 197);

    auto ts_yzzz_xxyy = pbuffer.data(idx_ovl_gg + 198);

    auto ts_yzzz_xxyz = pbuffer.data(idx_ovl_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_ovl_gg + 200);

    auto ts_yzzz_xyyy = pbuffer.data(idx_ovl_gg + 201);

    auto ts_yzzz_xyyz = pbuffer.data(idx_ovl_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_ovl_gg + 203);

    auto ts_yzzz_xzzz = pbuffer.data(idx_ovl_gg + 204);

    auto ts_yzzz_yyyy = pbuffer.data(idx_ovl_gg + 205);

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

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_ovl_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_ovl_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_ovl_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_ovl_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_ovl_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_ovl_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_ovl_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_ovl_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_ovl_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_ovl_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_ovl_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_ovl_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_ovl_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_ovl_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_ovl_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_ovl_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_ovl_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_ovl_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_ovl_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_ovl_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_ovl_gh + 20);

    auto ts_xxxy_xxxxx = pbuffer.data(idx_ovl_gh + 21);

    auto ts_xxxy_xxxxy = pbuffer.data(idx_ovl_gh + 22);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_ovl_gh + 23);

    auto ts_xxxy_xxxyy = pbuffer.data(idx_ovl_gh + 24);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_ovl_gh + 26);

    auto ts_xxxy_xxyyy = pbuffer.data(idx_ovl_gh + 27);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_ovl_gh + 30);

    auto ts_xxxy_xyyyy = pbuffer.data(idx_ovl_gh + 31);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_ovl_gh + 35);

    auto ts_xxxy_yyyyy = pbuffer.data(idx_ovl_gh + 36);

    auto ts_xxxy_yyyyz = pbuffer.data(idx_ovl_gh + 37);

    auto ts_xxxy_yyyzz = pbuffer.data(idx_ovl_gh + 38);

    auto ts_xxxy_yyzzz = pbuffer.data(idx_ovl_gh + 39);

    auto ts_xxxy_yzzzz = pbuffer.data(idx_ovl_gh + 40);

    auto ts_xxxz_xxxxx = pbuffer.data(idx_ovl_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_ovl_gh + 43);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_ovl_gh + 44);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_ovl_gh + 45);

    auto ts_xxxz_xxxyz = pbuffer.data(idx_ovl_gh + 46);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_ovl_gh + 47);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_ovl_gh + 48);

    auto ts_xxxz_xxyyz = pbuffer.data(idx_ovl_gh + 49);

    auto ts_xxxz_xxyzz = pbuffer.data(idx_ovl_gh + 50);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_ovl_gh + 51);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_ovl_gh + 52);

    auto ts_xxxz_xyyyz = pbuffer.data(idx_ovl_gh + 53);

    auto ts_xxxz_xyyzz = pbuffer.data(idx_ovl_gh + 54);

    auto ts_xxxz_xyzzz = pbuffer.data(idx_ovl_gh + 55);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_ovl_gh + 56);

    auto ts_xxxz_yyyyz = pbuffer.data(idx_ovl_gh + 58);

    auto ts_xxxz_yyyzz = pbuffer.data(idx_ovl_gh + 59);

    auto ts_xxxz_yyzzz = pbuffer.data(idx_ovl_gh + 60);

    auto ts_xxxz_yzzzz = pbuffer.data(idx_ovl_gh + 61);

    auto ts_xxxz_zzzzz = pbuffer.data(idx_ovl_gh + 62);

    auto ts_xxyy_xxxxx = pbuffer.data(idx_ovl_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_ovl_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_ovl_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_ovl_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_ovl_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_ovl_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_ovl_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_ovl_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_ovl_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_ovl_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_ovl_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_ovl_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_ovl_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_ovl_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_ovl_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_ovl_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_ovl_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_ovl_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_ovl_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_ovl_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_ovl_gh + 83);

    auto ts_xxyz_xxxxz = pbuffer.data(idx_ovl_gh + 86);

    auto ts_xxyz_xxxzz = pbuffer.data(idx_ovl_gh + 89);

    auto ts_xxyz_xxzzz = pbuffer.data(idx_ovl_gh + 93);

    auto ts_xxyz_xzzzz = pbuffer.data(idx_ovl_gh + 98);

    auto ts_xxyz_yyyyz = pbuffer.data(idx_ovl_gh + 100);

    auto ts_xxyz_yyyzz = pbuffer.data(idx_ovl_gh + 101);

    auto ts_xxyz_yyzzz = pbuffer.data(idx_ovl_gh + 102);

    auto ts_xxyz_yzzzz = pbuffer.data(idx_ovl_gh + 103);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_ovl_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_ovl_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_ovl_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_ovl_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_ovl_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_ovl_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_ovl_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_ovl_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_ovl_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_ovl_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_ovl_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_ovl_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_ovl_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_ovl_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_ovl_gh + 119);

    auto ts_xxzz_yyyyy = pbuffer.data(idx_ovl_gh + 120);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_ovl_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_ovl_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_ovl_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_ovl_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_ovl_gh + 125);

    auto ts_xyyy_xxxxx = pbuffer.data(idx_ovl_gh + 126);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_ovl_gh + 127);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_ovl_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_ovl_gh + 130);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_ovl_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_ovl_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_ovl_gh + 134);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_ovl_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_ovl_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_ovl_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_ovl_gh + 139);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_ovl_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_ovl_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_ovl_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_ovl_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_ovl_gh + 145);

    auto ts_xyyy_zzzzz = pbuffer.data(idx_ovl_gh + 146);

    auto ts_xyyz_yyyyz = pbuffer.data(idx_ovl_gh + 163);

    auto ts_xyyz_yyyzz = pbuffer.data(idx_ovl_gh + 164);

    auto ts_xyyz_yyzzz = pbuffer.data(idx_ovl_gh + 165);

    auto ts_xyyz_yzzzz = pbuffer.data(idx_ovl_gh + 166);

    auto ts_xyyz_zzzzz = pbuffer.data(idx_ovl_gh + 167);

    auto ts_xyzz_yyyyy = pbuffer.data(idx_ovl_gh + 183);

    auto ts_xyzz_yyyyz = pbuffer.data(idx_ovl_gh + 184);

    auto ts_xyzz_yyyzz = pbuffer.data(idx_ovl_gh + 185);

    auto ts_xyzz_yyzzz = pbuffer.data(idx_ovl_gh + 186);

    auto ts_xyzz_yzzzz = pbuffer.data(idx_ovl_gh + 187);

    auto ts_xzzz_xxxxx = pbuffer.data(idx_ovl_gh + 189);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_ovl_gh + 191);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_ovl_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_ovl_gh + 194);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_ovl_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_ovl_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_ovl_gh + 198);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_ovl_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_ovl_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_ovl_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_ovl_gh + 203);

    auto ts_xzzz_yyyyy = pbuffer.data(idx_ovl_gh + 204);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_ovl_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_ovl_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_ovl_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_ovl_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_ovl_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_ovl_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_ovl_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_ovl_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_ovl_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_ovl_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_ovl_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_ovl_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_ovl_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_ovl_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_ovl_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_ovl_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_ovl_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_ovl_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_ovl_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_ovl_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_ovl_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_ovl_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_ovl_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_ovl_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_ovl_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_ovl_gh + 230);

    auto ts_yyyz_xxxxy = pbuffer.data(idx_ovl_gh + 232);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_ovl_gh + 233);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_ovl_gh + 234);

    auto ts_yyyz_xxxyz = pbuffer.data(idx_ovl_gh + 235);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_ovl_gh + 236);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_ovl_gh + 237);

    auto ts_yyyz_xxyyz = pbuffer.data(idx_ovl_gh + 238);

    auto ts_yyyz_xxyzz = pbuffer.data(idx_ovl_gh + 239);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_ovl_gh + 240);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_ovl_gh + 241);

    auto ts_yyyz_xyyyz = pbuffer.data(idx_ovl_gh + 242);

    auto ts_yyyz_xyyzz = pbuffer.data(idx_ovl_gh + 243);

    auto ts_yyyz_xyzzz = pbuffer.data(idx_ovl_gh + 244);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_ovl_gh + 245);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_ovl_gh + 246);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_ovl_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_ovl_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_ovl_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_ovl_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_ovl_gh + 251);

    auto ts_yyzz_xxxxx = pbuffer.data(idx_ovl_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_ovl_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_ovl_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_ovl_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_ovl_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_ovl_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_ovl_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_ovl_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_ovl_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_ovl_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_ovl_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_ovl_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_ovl_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_ovl_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_ovl_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_ovl_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_ovl_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_ovl_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_ovl_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_ovl_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_ovl_gh + 272);

    auto ts_yzzz_xxxxx = pbuffer.data(idx_ovl_gh + 273);

    auto ts_yzzz_xxxxy = pbuffer.data(idx_ovl_gh + 274);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_ovl_gh + 275);

    auto ts_yzzz_xxxyy = pbuffer.data(idx_ovl_gh + 276);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_ovl_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_ovl_gh + 278);

    auto ts_yzzz_xxyyy = pbuffer.data(idx_ovl_gh + 279);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_ovl_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_ovl_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_ovl_gh + 282);

    auto ts_yzzz_xyyyy = pbuffer.data(idx_ovl_gh + 283);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_ovl_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_ovl_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_ovl_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_ovl_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_ovl_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_ovl_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_ovl_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_ovl_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_ovl_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_ovl_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_ovl_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_ovl_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_ovl_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_ovl_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_ovl_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_ovl_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_ovl_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_ovl_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_ovl_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_ovl_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_ovl_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_ovl_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_ovl_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_ovl_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_ovl_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_ovl_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_ovl_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_ovl_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_ovl_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_ovl_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_ovl_gh + 314);

    // Set up 0-21 components of targeted buffer : HH

    auto ts_xxxxx_xxxxx = pbuffer.data(idx_ovl_hh);

    auto ts_xxxxx_xxxxy = pbuffer.data(idx_ovl_hh + 1);

    auto ts_xxxxx_xxxxz = pbuffer.data(idx_ovl_hh + 2);

    auto ts_xxxxx_xxxyy = pbuffer.data(idx_ovl_hh + 3);

    auto ts_xxxxx_xxxyz = pbuffer.data(idx_ovl_hh + 4);

    auto ts_xxxxx_xxxzz = pbuffer.data(idx_ovl_hh + 5);

    auto ts_xxxxx_xxyyy = pbuffer.data(idx_ovl_hh + 6);

    auto ts_xxxxx_xxyyz = pbuffer.data(idx_ovl_hh + 7);

    auto ts_xxxxx_xxyzz = pbuffer.data(idx_ovl_hh + 8);

    auto ts_xxxxx_xxzzz = pbuffer.data(idx_ovl_hh + 9);

    auto ts_xxxxx_xyyyy = pbuffer.data(idx_ovl_hh + 10);

    auto ts_xxxxx_xyyyz = pbuffer.data(idx_ovl_hh + 11);

    auto ts_xxxxx_xyyzz = pbuffer.data(idx_ovl_hh + 12);

    auto ts_xxxxx_xyzzz = pbuffer.data(idx_ovl_hh + 13);

    auto ts_xxxxx_xzzzz = pbuffer.data(idx_ovl_hh + 14);

    auto ts_xxxxx_yyyyy = pbuffer.data(idx_ovl_hh + 15);

    auto ts_xxxxx_yyyyz = pbuffer.data(idx_ovl_hh + 16);

    auto ts_xxxxx_yyyzz = pbuffer.data(idx_ovl_hh + 17);

    auto ts_xxxxx_yyzzz = pbuffer.data(idx_ovl_hh + 18);

    auto ts_xxxxx_yzzzz = pbuffer.data(idx_ovl_hh + 19);

    auto ts_xxxxx_zzzzz = pbuffer.data(idx_ovl_hh + 20);

    #pragma omp simd aligned(pa_x, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxzz, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyzz, ts_xxx_xxzzz, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyzz, ts_xxx_xyzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyzz, ts_xxx_yyzzz, ts_xxx_yzzzz, ts_xxx_zzzzz, ts_xxxx_xxxx, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxy, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxz, ts_xxxx_xxxzz, ts_xxxx_xxyy, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyz, ts_xxxx_xxyzz, ts_xxxx_xxzz, ts_xxxx_xxzzz, ts_xxxx_xyyy, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyz, ts_xxxx_xyyzz, ts_xxxx_xyzz, ts_xxxx_xyzzz, ts_xxxx_xzzz, ts_xxxx_xzzzz, ts_xxxx_yyyy, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyz, ts_xxxx_yyyzz, ts_xxxx_yyzz, ts_xxxx_yyzzz, ts_xxxx_yzzz, ts_xxxx_yzzzz, ts_xxxx_zzzz, ts_xxxx_zzzzz, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxzz, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyyy, ts_xxxxx_yyyyz, ts_xxxxx_yyyzz, ts_xxxxx_yyzzz, ts_xxxxx_yzzzz, ts_xxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxx_xxxxx[i] = 4.0 * ts_xxx_xxxxx[i] * fe_0 + 5.0 * ts_xxxx_xxxx[i] * fe_0 + ts_xxxx_xxxxx[i] * pa_x[i];

        ts_xxxxx_xxxxy[i] = 4.0 * ts_xxx_xxxxy[i] * fe_0 + 4.0 * ts_xxxx_xxxy[i] * fe_0 + ts_xxxx_xxxxy[i] * pa_x[i];

        ts_xxxxx_xxxxz[i] = 4.0 * ts_xxx_xxxxz[i] * fe_0 + 4.0 * ts_xxxx_xxxz[i] * fe_0 + ts_xxxx_xxxxz[i] * pa_x[i];

        ts_xxxxx_xxxyy[i] = 4.0 * ts_xxx_xxxyy[i] * fe_0 + 3.0 * ts_xxxx_xxyy[i] * fe_0 + ts_xxxx_xxxyy[i] * pa_x[i];

        ts_xxxxx_xxxyz[i] = 4.0 * ts_xxx_xxxyz[i] * fe_0 + 3.0 * ts_xxxx_xxyz[i] * fe_0 + ts_xxxx_xxxyz[i] * pa_x[i];

        ts_xxxxx_xxxzz[i] = 4.0 * ts_xxx_xxxzz[i] * fe_0 + 3.0 * ts_xxxx_xxzz[i] * fe_0 + ts_xxxx_xxxzz[i] * pa_x[i];

        ts_xxxxx_xxyyy[i] = 4.0 * ts_xxx_xxyyy[i] * fe_0 + 2.0 * ts_xxxx_xyyy[i] * fe_0 + ts_xxxx_xxyyy[i] * pa_x[i];

        ts_xxxxx_xxyyz[i] = 4.0 * ts_xxx_xxyyz[i] * fe_0 + 2.0 * ts_xxxx_xyyz[i] * fe_0 + ts_xxxx_xxyyz[i] * pa_x[i];

        ts_xxxxx_xxyzz[i] = 4.0 * ts_xxx_xxyzz[i] * fe_0 + 2.0 * ts_xxxx_xyzz[i] * fe_0 + ts_xxxx_xxyzz[i] * pa_x[i];

        ts_xxxxx_xxzzz[i] = 4.0 * ts_xxx_xxzzz[i] * fe_0 + 2.0 * ts_xxxx_xzzz[i] * fe_0 + ts_xxxx_xxzzz[i] * pa_x[i];

        ts_xxxxx_xyyyy[i] = 4.0 * ts_xxx_xyyyy[i] * fe_0 + ts_xxxx_yyyy[i] * fe_0 + ts_xxxx_xyyyy[i] * pa_x[i];

        ts_xxxxx_xyyyz[i] = 4.0 * ts_xxx_xyyyz[i] * fe_0 + ts_xxxx_yyyz[i] * fe_0 + ts_xxxx_xyyyz[i] * pa_x[i];

        ts_xxxxx_xyyzz[i] = 4.0 * ts_xxx_xyyzz[i] * fe_0 + ts_xxxx_yyzz[i] * fe_0 + ts_xxxx_xyyzz[i] * pa_x[i];

        ts_xxxxx_xyzzz[i] = 4.0 * ts_xxx_xyzzz[i] * fe_0 + ts_xxxx_yzzz[i] * fe_0 + ts_xxxx_xyzzz[i] * pa_x[i];

        ts_xxxxx_xzzzz[i] = 4.0 * ts_xxx_xzzzz[i] * fe_0 + ts_xxxx_zzzz[i] * fe_0 + ts_xxxx_xzzzz[i] * pa_x[i];

        ts_xxxxx_yyyyy[i] = 4.0 * ts_xxx_yyyyy[i] * fe_0 + ts_xxxx_yyyyy[i] * pa_x[i];

        ts_xxxxx_yyyyz[i] = 4.0 * ts_xxx_yyyyz[i] * fe_0 + ts_xxxx_yyyyz[i] * pa_x[i];

        ts_xxxxx_yyyzz[i] = 4.0 * ts_xxx_yyyzz[i] * fe_0 + ts_xxxx_yyyzz[i] * pa_x[i];

        ts_xxxxx_yyzzz[i] = 4.0 * ts_xxx_yyzzz[i] * fe_0 + ts_xxxx_yyzzz[i] * pa_x[i];

        ts_xxxxx_yzzzz[i] = 4.0 * ts_xxx_yzzzz[i] * fe_0 + ts_xxxx_yzzzz[i] * pa_x[i];

        ts_xxxxx_zzzzz[i] = 4.0 * ts_xxx_zzzzz[i] * fe_0 + ts_xxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : HH

    auto ts_xxxxy_xxxxx = pbuffer.data(idx_ovl_hh + 21);

    auto ts_xxxxy_xxxxy = pbuffer.data(idx_ovl_hh + 22);

    auto ts_xxxxy_xxxxz = pbuffer.data(idx_ovl_hh + 23);

    auto ts_xxxxy_xxxyy = pbuffer.data(idx_ovl_hh + 24);

    auto ts_xxxxy_xxxyz = pbuffer.data(idx_ovl_hh + 25);

    auto ts_xxxxy_xxxzz = pbuffer.data(idx_ovl_hh + 26);

    auto ts_xxxxy_xxyyy = pbuffer.data(idx_ovl_hh + 27);

    auto ts_xxxxy_xxyyz = pbuffer.data(idx_ovl_hh + 28);

    auto ts_xxxxy_xxyzz = pbuffer.data(idx_ovl_hh + 29);

    auto ts_xxxxy_xxzzz = pbuffer.data(idx_ovl_hh + 30);

    auto ts_xxxxy_xyyyy = pbuffer.data(idx_ovl_hh + 31);

    auto ts_xxxxy_xyyyz = pbuffer.data(idx_ovl_hh + 32);

    auto ts_xxxxy_xyyzz = pbuffer.data(idx_ovl_hh + 33);

    auto ts_xxxxy_xyzzz = pbuffer.data(idx_ovl_hh + 34);

    auto ts_xxxxy_xzzzz = pbuffer.data(idx_ovl_hh + 35);

    auto ts_xxxxy_yyyyy = pbuffer.data(idx_ovl_hh + 36);

    auto ts_xxxxy_yyyyz = pbuffer.data(idx_ovl_hh + 37);

    auto ts_xxxxy_yyyzz = pbuffer.data(idx_ovl_hh + 38);

    auto ts_xxxxy_yyzzz = pbuffer.data(idx_ovl_hh + 39);

    auto ts_xxxxy_yzzzz = pbuffer.data(idx_ovl_hh + 40);

    auto ts_xxxxy_zzzzz = pbuffer.data(idx_ovl_hh + 41);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xxxx, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxy, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxz, ts_xxxx_xxxzz, ts_xxxx_xxyy, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyz, ts_xxxx_xxyzz, ts_xxxx_xxzz, ts_xxxx_xxzzz, ts_xxxx_xyyy, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyz, ts_xxxx_xyyzz, ts_xxxx_xyzz, ts_xxxx_xyzzz, ts_xxxx_xzzz, ts_xxxx_xzzzz, ts_xxxx_zzzzz, ts_xxxxy_xxxxx, ts_xxxxy_xxxxy, ts_xxxxy_xxxxz, ts_xxxxy_xxxyy, ts_xxxxy_xxxyz, ts_xxxxy_xxxzz, ts_xxxxy_xxyyy, ts_xxxxy_xxyyz, ts_xxxxy_xxyzz, ts_xxxxy_xxzzz, ts_xxxxy_xyyyy, ts_xxxxy_xyyyz, ts_xxxxy_xyyzz, ts_xxxxy_xyzzz, ts_xxxxy_xzzzz, ts_xxxxy_yyyyy, ts_xxxxy_yyyyz, ts_xxxxy_yyyzz, ts_xxxxy_yyzzz, ts_xxxxy_yzzzz, ts_xxxxy_zzzzz, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyzz, ts_xxxy_yyzzz, ts_xxxy_yzzzz, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyzz, ts_xxy_yyzzz, ts_xxy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxy_xxxxx[i] = ts_xxxx_xxxxx[i] * pa_y[i];

        ts_xxxxy_xxxxy[i] = ts_xxxx_xxxx[i] * fe_0 + ts_xxxx_xxxxy[i] * pa_y[i];

        ts_xxxxy_xxxxz[i] = ts_xxxx_xxxxz[i] * pa_y[i];

        ts_xxxxy_xxxyy[i] = 2.0 * ts_xxxx_xxxy[i] * fe_0 + ts_xxxx_xxxyy[i] * pa_y[i];

        ts_xxxxy_xxxyz[i] = ts_xxxx_xxxz[i] * fe_0 + ts_xxxx_xxxyz[i] * pa_y[i];

        ts_xxxxy_xxxzz[i] = ts_xxxx_xxxzz[i] * pa_y[i];

        ts_xxxxy_xxyyy[i] = 3.0 * ts_xxxx_xxyy[i] * fe_0 + ts_xxxx_xxyyy[i] * pa_y[i];

        ts_xxxxy_xxyyz[i] = 2.0 * ts_xxxx_xxyz[i] * fe_0 + ts_xxxx_xxyyz[i] * pa_y[i];

        ts_xxxxy_xxyzz[i] = ts_xxxx_xxzz[i] * fe_0 + ts_xxxx_xxyzz[i] * pa_y[i];

        ts_xxxxy_xxzzz[i] = ts_xxxx_xxzzz[i] * pa_y[i];

        ts_xxxxy_xyyyy[i] = 4.0 * ts_xxxx_xyyy[i] * fe_0 + ts_xxxx_xyyyy[i] * pa_y[i];

        ts_xxxxy_xyyyz[i] = 3.0 * ts_xxxx_xyyz[i] * fe_0 + ts_xxxx_xyyyz[i] * pa_y[i];

        ts_xxxxy_xyyzz[i] = 2.0 * ts_xxxx_xyzz[i] * fe_0 + ts_xxxx_xyyzz[i] * pa_y[i];

        ts_xxxxy_xyzzz[i] = ts_xxxx_xzzz[i] * fe_0 + ts_xxxx_xyzzz[i] * pa_y[i];

        ts_xxxxy_xzzzz[i] = ts_xxxx_xzzzz[i] * pa_y[i];

        ts_xxxxy_yyyyy[i] = 3.0 * ts_xxy_yyyyy[i] * fe_0 + ts_xxxy_yyyyy[i] * pa_x[i];

        ts_xxxxy_yyyyz[i] = 3.0 * ts_xxy_yyyyz[i] * fe_0 + ts_xxxy_yyyyz[i] * pa_x[i];

        ts_xxxxy_yyyzz[i] = 3.0 * ts_xxy_yyyzz[i] * fe_0 + ts_xxxy_yyyzz[i] * pa_x[i];

        ts_xxxxy_yyzzz[i] = 3.0 * ts_xxy_yyzzz[i] * fe_0 + ts_xxxy_yyzzz[i] * pa_x[i];

        ts_xxxxy_yzzzz[i] = 3.0 * ts_xxy_yzzzz[i] * fe_0 + ts_xxxy_yzzzz[i] * pa_x[i];

        ts_xxxxy_zzzzz[i] = ts_xxxx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : HH

    auto ts_xxxxz_xxxxx = pbuffer.data(idx_ovl_hh + 42);

    auto ts_xxxxz_xxxxy = pbuffer.data(idx_ovl_hh + 43);

    auto ts_xxxxz_xxxxz = pbuffer.data(idx_ovl_hh + 44);

    auto ts_xxxxz_xxxyy = pbuffer.data(idx_ovl_hh + 45);

    auto ts_xxxxz_xxxyz = pbuffer.data(idx_ovl_hh + 46);

    auto ts_xxxxz_xxxzz = pbuffer.data(idx_ovl_hh + 47);

    auto ts_xxxxz_xxyyy = pbuffer.data(idx_ovl_hh + 48);

    auto ts_xxxxz_xxyyz = pbuffer.data(idx_ovl_hh + 49);

    auto ts_xxxxz_xxyzz = pbuffer.data(idx_ovl_hh + 50);

    auto ts_xxxxz_xxzzz = pbuffer.data(idx_ovl_hh + 51);

    auto ts_xxxxz_xyyyy = pbuffer.data(idx_ovl_hh + 52);

    auto ts_xxxxz_xyyyz = pbuffer.data(idx_ovl_hh + 53);

    auto ts_xxxxz_xyyzz = pbuffer.data(idx_ovl_hh + 54);

    auto ts_xxxxz_xyzzz = pbuffer.data(idx_ovl_hh + 55);

    auto ts_xxxxz_xzzzz = pbuffer.data(idx_ovl_hh + 56);

    auto ts_xxxxz_yyyyy = pbuffer.data(idx_ovl_hh + 57);

    auto ts_xxxxz_yyyyz = pbuffer.data(idx_ovl_hh + 58);

    auto ts_xxxxz_yyyzz = pbuffer.data(idx_ovl_hh + 59);

    auto ts_xxxxz_yyzzz = pbuffer.data(idx_ovl_hh + 60);

    auto ts_xxxxz_yzzzz = pbuffer.data(idx_ovl_hh + 61);

    auto ts_xxxxz_zzzzz = pbuffer.data(idx_ovl_hh + 62);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xxxx, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxy, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxz, ts_xxxx_xxxzz, ts_xxxx_xxyy, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyz, ts_xxxx_xxyzz, ts_xxxx_xxzz, ts_xxxx_xxzzz, ts_xxxx_xyyy, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyz, ts_xxxx_xyyzz, ts_xxxx_xyzz, ts_xxxx_xyzzz, ts_xxxx_xzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxxz_xxxxx, ts_xxxxz_xxxxy, ts_xxxxz_xxxxz, ts_xxxxz_xxxyy, ts_xxxxz_xxxyz, ts_xxxxz_xxxzz, ts_xxxxz_xxyyy, ts_xxxxz_xxyyz, ts_xxxxz_xxyzz, ts_xxxxz_xxzzz, ts_xxxxz_xyyyy, ts_xxxxz_xyyyz, ts_xxxxz_xyyzz, ts_xxxxz_xyzzz, ts_xxxxz_xzzzz, ts_xxxxz_yyyyy, ts_xxxxz_yyyyz, ts_xxxxz_yyyzz, ts_xxxxz_yyzzz, ts_xxxxz_yzzzz, ts_xxxxz_zzzzz, ts_xxxz_yyyyz, ts_xxxz_yyyzz, ts_xxxz_yyzzz, ts_xxxz_yzzzz, ts_xxxz_zzzzz, ts_xxz_yyyyz, ts_xxz_yyyzz, ts_xxz_yyzzz, ts_xxz_yzzzz, ts_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxz_xxxxx[i] = ts_xxxx_xxxxx[i] * pa_z[i];

        ts_xxxxz_xxxxy[i] = ts_xxxx_xxxxy[i] * pa_z[i];

        ts_xxxxz_xxxxz[i] = ts_xxxx_xxxx[i] * fe_0 + ts_xxxx_xxxxz[i] * pa_z[i];

        ts_xxxxz_xxxyy[i] = ts_xxxx_xxxyy[i] * pa_z[i];

        ts_xxxxz_xxxyz[i] = ts_xxxx_xxxy[i] * fe_0 + ts_xxxx_xxxyz[i] * pa_z[i];

        ts_xxxxz_xxxzz[i] = 2.0 * ts_xxxx_xxxz[i] * fe_0 + ts_xxxx_xxxzz[i] * pa_z[i];

        ts_xxxxz_xxyyy[i] = ts_xxxx_xxyyy[i] * pa_z[i];

        ts_xxxxz_xxyyz[i] = ts_xxxx_xxyy[i] * fe_0 + ts_xxxx_xxyyz[i] * pa_z[i];

        ts_xxxxz_xxyzz[i] = 2.0 * ts_xxxx_xxyz[i] * fe_0 + ts_xxxx_xxyzz[i] * pa_z[i];

        ts_xxxxz_xxzzz[i] = 3.0 * ts_xxxx_xxzz[i] * fe_0 + ts_xxxx_xxzzz[i] * pa_z[i];

        ts_xxxxz_xyyyy[i] = ts_xxxx_xyyyy[i] * pa_z[i];

        ts_xxxxz_xyyyz[i] = ts_xxxx_xyyy[i] * fe_0 + ts_xxxx_xyyyz[i] * pa_z[i];

        ts_xxxxz_xyyzz[i] = 2.0 * ts_xxxx_xyyz[i] * fe_0 + ts_xxxx_xyyzz[i] * pa_z[i];

        ts_xxxxz_xyzzz[i] = 3.0 * ts_xxxx_xyzz[i] * fe_0 + ts_xxxx_xyzzz[i] * pa_z[i];

        ts_xxxxz_xzzzz[i] = 4.0 * ts_xxxx_xzzz[i] * fe_0 + ts_xxxx_xzzzz[i] * pa_z[i];

        ts_xxxxz_yyyyy[i] = ts_xxxx_yyyyy[i] * pa_z[i];

        ts_xxxxz_yyyyz[i] = 3.0 * ts_xxz_yyyyz[i] * fe_0 + ts_xxxz_yyyyz[i] * pa_x[i];

        ts_xxxxz_yyyzz[i] = 3.0 * ts_xxz_yyyzz[i] * fe_0 + ts_xxxz_yyyzz[i] * pa_x[i];

        ts_xxxxz_yyzzz[i] = 3.0 * ts_xxz_yyzzz[i] * fe_0 + ts_xxxz_yyzzz[i] * pa_x[i];

        ts_xxxxz_yzzzz[i] = 3.0 * ts_xxz_yzzzz[i] * fe_0 + ts_xxxz_yzzzz[i] * pa_x[i];

        ts_xxxxz_zzzzz[i] = 3.0 * ts_xxz_zzzzz[i] * fe_0 + ts_xxxz_zzzzz[i] * pa_x[i];
    }

    // Set up 63-84 components of targeted buffer : HH

    auto ts_xxxyy_xxxxx = pbuffer.data(idx_ovl_hh + 63);

    auto ts_xxxyy_xxxxy = pbuffer.data(idx_ovl_hh + 64);

    auto ts_xxxyy_xxxxz = pbuffer.data(idx_ovl_hh + 65);

    auto ts_xxxyy_xxxyy = pbuffer.data(idx_ovl_hh + 66);

    auto ts_xxxyy_xxxyz = pbuffer.data(idx_ovl_hh + 67);

    auto ts_xxxyy_xxxzz = pbuffer.data(idx_ovl_hh + 68);

    auto ts_xxxyy_xxyyy = pbuffer.data(idx_ovl_hh + 69);

    auto ts_xxxyy_xxyyz = pbuffer.data(idx_ovl_hh + 70);

    auto ts_xxxyy_xxyzz = pbuffer.data(idx_ovl_hh + 71);

    auto ts_xxxyy_xxzzz = pbuffer.data(idx_ovl_hh + 72);

    auto ts_xxxyy_xyyyy = pbuffer.data(idx_ovl_hh + 73);

    auto ts_xxxyy_xyyyz = pbuffer.data(idx_ovl_hh + 74);

    auto ts_xxxyy_xyyzz = pbuffer.data(idx_ovl_hh + 75);

    auto ts_xxxyy_xyzzz = pbuffer.data(idx_ovl_hh + 76);

    auto ts_xxxyy_xzzzz = pbuffer.data(idx_ovl_hh + 77);

    auto ts_xxxyy_yyyyy = pbuffer.data(idx_ovl_hh + 78);

    auto ts_xxxyy_yyyyz = pbuffer.data(idx_ovl_hh + 79);

    auto ts_xxxyy_yyyzz = pbuffer.data(idx_ovl_hh + 80);

    auto ts_xxxyy_yyzzz = pbuffer.data(idx_ovl_hh + 81);

    auto ts_xxxyy_yzzzz = pbuffer.data(idx_ovl_hh + 82);

    auto ts_xxxyy_zzzzz = pbuffer.data(idx_ovl_hh + 83);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxx_xxxxx, ts_xxx_xxxxz, ts_xxx_xxxzz, ts_xxx_xxzzz, ts_xxx_xzzzz, ts_xxxy_xxxxx, ts_xxxy_xxxxz, ts_xxxy_xxxzz, ts_xxxy_xxzzz, ts_xxxy_xzzzz, ts_xxxyy_xxxxx, ts_xxxyy_xxxxy, ts_xxxyy_xxxxz, ts_xxxyy_xxxyy, ts_xxxyy_xxxyz, ts_xxxyy_xxxzz, ts_xxxyy_xxyyy, ts_xxxyy_xxyyz, ts_xxxyy_xxyzz, ts_xxxyy_xxzzz, ts_xxxyy_xyyyy, ts_xxxyy_xyyyz, ts_xxxyy_xyyzz, ts_xxxyy_xyzzz, ts_xxxyy_xzzzz, ts_xxxyy_yyyyy, ts_xxxyy_yyyyz, ts_xxxyy_yyyzz, ts_xxxyy_yyzzz, ts_xxxyy_yzzzz, ts_xxxyy_zzzzz, ts_xxyy_xxxxy, ts_xxyy_xxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxyy, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyz, ts_xxyy_xxyzz, ts_xxyy_xyyy, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyz, ts_xxyy_xyyzz, ts_xxyy_xyzz, ts_xxyy_xyzzz, ts_xxyy_yyyy, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyz, ts_xxyy_yyyzz, ts_xxyy_yyzz, ts_xxyy_yyzzz, ts_xxyy_yzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, ts_xyy_xxxxy, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyzz, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyzz, ts_xyy_xyzzz, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyzz, ts_xyy_yyzzz, ts_xyy_yzzzz, ts_xyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyy_xxxxx[i] = ts_xxx_xxxxx[i] * fe_0 + ts_xxxy_xxxxx[i] * pa_y[i];

        ts_xxxyy_xxxxy[i] = 2.0 * ts_xyy_xxxxy[i] * fe_0 + 4.0 * ts_xxyy_xxxy[i] * fe_0 + ts_xxyy_xxxxy[i] * pa_x[i];

        ts_xxxyy_xxxxz[i] = ts_xxx_xxxxz[i] * fe_0 + ts_xxxy_xxxxz[i] * pa_y[i];

        ts_xxxyy_xxxyy[i] = 2.0 * ts_xyy_xxxyy[i] * fe_0 + 3.0 * ts_xxyy_xxyy[i] * fe_0 + ts_xxyy_xxxyy[i] * pa_x[i];

        ts_xxxyy_xxxyz[i] = 2.0 * ts_xyy_xxxyz[i] * fe_0 + 3.0 * ts_xxyy_xxyz[i] * fe_0 + ts_xxyy_xxxyz[i] * pa_x[i];

        ts_xxxyy_xxxzz[i] = ts_xxx_xxxzz[i] * fe_0 + ts_xxxy_xxxzz[i] * pa_y[i];

        ts_xxxyy_xxyyy[i] = 2.0 * ts_xyy_xxyyy[i] * fe_0 + 2.0 * ts_xxyy_xyyy[i] * fe_0 + ts_xxyy_xxyyy[i] * pa_x[i];

        ts_xxxyy_xxyyz[i] = 2.0 * ts_xyy_xxyyz[i] * fe_0 + 2.0 * ts_xxyy_xyyz[i] * fe_0 + ts_xxyy_xxyyz[i] * pa_x[i];

        ts_xxxyy_xxyzz[i] = 2.0 * ts_xyy_xxyzz[i] * fe_0 + 2.0 * ts_xxyy_xyzz[i] * fe_0 + ts_xxyy_xxyzz[i] * pa_x[i];

        ts_xxxyy_xxzzz[i] = ts_xxx_xxzzz[i] * fe_0 + ts_xxxy_xxzzz[i] * pa_y[i];

        ts_xxxyy_xyyyy[i] = 2.0 * ts_xyy_xyyyy[i] * fe_0 + ts_xxyy_yyyy[i] * fe_0 + ts_xxyy_xyyyy[i] * pa_x[i];

        ts_xxxyy_xyyyz[i] = 2.0 * ts_xyy_xyyyz[i] * fe_0 + ts_xxyy_yyyz[i] * fe_0 + ts_xxyy_xyyyz[i] * pa_x[i];

        ts_xxxyy_xyyzz[i] = 2.0 * ts_xyy_xyyzz[i] * fe_0 + ts_xxyy_yyzz[i] * fe_0 + ts_xxyy_xyyzz[i] * pa_x[i];

        ts_xxxyy_xyzzz[i] = 2.0 * ts_xyy_xyzzz[i] * fe_0 + ts_xxyy_yzzz[i] * fe_0 + ts_xxyy_xyzzz[i] * pa_x[i];

        ts_xxxyy_xzzzz[i] = ts_xxx_xzzzz[i] * fe_0 + ts_xxxy_xzzzz[i] * pa_y[i];

        ts_xxxyy_yyyyy[i] = 2.0 * ts_xyy_yyyyy[i] * fe_0 + ts_xxyy_yyyyy[i] * pa_x[i];

        ts_xxxyy_yyyyz[i] = 2.0 * ts_xyy_yyyyz[i] * fe_0 + ts_xxyy_yyyyz[i] * pa_x[i];

        ts_xxxyy_yyyzz[i] = 2.0 * ts_xyy_yyyzz[i] * fe_0 + ts_xxyy_yyyzz[i] * pa_x[i];

        ts_xxxyy_yyzzz[i] = 2.0 * ts_xyy_yyzzz[i] * fe_0 + ts_xxyy_yyzzz[i] * pa_x[i];

        ts_xxxyy_yzzzz[i] = 2.0 * ts_xyy_yzzzz[i] * fe_0 + ts_xxyy_yzzzz[i] * pa_x[i];

        ts_xxxyy_zzzzz[i] = 2.0 * ts_xyy_zzzzz[i] * fe_0 + ts_xxyy_zzzzz[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : HH

    auto ts_xxxyz_xxxxx = pbuffer.data(idx_ovl_hh + 84);

    auto ts_xxxyz_xxxxy = pbuffer.data(idx_ovl_hh + 85);

    auto ts_xxxyz_xxxxz = pbuffer.data(idx_ovl_hh + 86);

    auto ts_xxxyz_xxxyy = pbuffer.data(idx_ovl_hh + 87);

    auto ts_xxxyz_xxxyz = pbuffer.data(idx_ovl_hh + 88);

    auto ts_xxxyz_xxxzz = pbuffer.data(idx_ovl_hh + 89);

    auto ts_xxxyz_xxyyy = pbuffer.data(idx_ovl_hh + 90);

    auto ts_xxxyz_xxyyz = pbuffer.data(idx_ovl_hh + 91);

    auto ts_xxxyz_xxyzz = pbuffer.data(idx_ovl_hh + 92);

    auto ts_xxxyz_xxzzz = pbuffer.data(idx_ovl_hh + 93);

    auto ts_xxxyz_xyyyy = pbuffer.data(idx_ovl_hh + 94);

    auto ts_xxxyz_xyyyz = pbuffer.data(idx_ovl_hh + 95);

    auto ts_xxxyz_xyyzz = pbuffer.data(idx_ovl_hh + 96);

    auto ts_xxxyz_xyzzz = pbuffer.data(idx_ovl_hh + 97);

    auto ts_xxxyz_xzzzz = pbuffer.data(idx_ovl_hh + 98);

    auto ts_xxxyz_yyyyy = pbuffer.data(idx_ovl_hh + 99);

    auto ts_xxxyz_yyyyz = pbuffer.data(idx_ovl_hh + 100);

    auto ts_xxxyz_yyyzz = pbuffer.data(idx_ovl_hh + 101);

    auto ts_xxxyz_yyzzz = pbuffer.data(idx_ovl_hh + 102);

    auto ts_xxxyz_yzzzz = pbuffer.data(idx_ovl_hh + 103);

    auto ts_xxxyz_zzzzz = pbuffer.data(idx_ovl_hh + 104);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxy_xxxxy, ts_xxxy_xxxyy, ts_xxxy_xxyyy, ts_xxxy_xyyyy, ts_xxxy_yyyyy, ts_xxxyz_xxxxx, ts_xxxyz_xxxxy, ts_xxxyz_xxxxz, ts_xxxyz_xxxyy, ts_xxxyz_xxxyz, ts_xxxyz_xxxzz, ts_xxxyz_xxyyy, ts_xxxyz_xxyyz, ts_xxxyz_xxyzz, ts_xxxyz_xxzzz, ts_xxxyz_xyyyy, ts_xxxyz_xyyyz, ts_xxxyz_xyyzz, ts_xxxyz_xyzzz, ts_xxxyz_xzzzz, ts_xxxyz_yyyyy, ts_xxxyz_yyyyz, ts_xxxyz_yyyzz, ts_xxxyz_yyzzz, ts_xxxyz_yzzzz, ts_xxxyz_zzzzz, ts_xxxz_xxxxx, ts_xxxz_xxxxz, ts_xxxz_xxxyz, ts_xxxz_xxxz, ts_xxxz_xxxzz, ts_xxxz_xxyyz, ts_xxxz_xxyz, ts_xxxz_xxyzz, ts_xxxz_xxzz, ts_xxxz_xxzzz, ts_xxxz_xyyyz, ts_xxxz_xyyz, ts_xxxz_xyyzz, ts_xxxz_xyzz, ts_xxxz_xyzzz, ts_xxxz_xzzz, ts_xxxz_xzzzz, ts_xxxz_zzzzz, ts_xxyz_yyyyz, ts_xxyz_yyyzz, ts_xxyz_yyzzz, ts_xxyz_yzzzz, ts_xyz_yyyyz, ts_xyz_yyyzz, ts_xyz_yyzzz, ts_xyz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyz_xxxxx[i] = ts_xxxz_xxxxx[i] * pa_y[i];

        ts_xxxyz_xxxxy[i] = ts_xxxy_xxxxy[i] * pa_z[i];

        ts_xxxyz_xxxxz[i] = ts_xxxz_xxxxz[i] * pa_y[i];

        ts_xxxyz_xxxyy[i] = ts_xxxy_xxxyy[i] * pa_z[i];

        ts_xxxyz_xxxyz[i] = ts_xxxz_xxxz[i] * fe_0 + ts_xxxz_xxxyz[i] * pa_y[i];

        ts_xxxyz_xxxzz[i] = ts_xxxz_xxxzz[i] * pa_y[i];

        ts_xxxyz_xxyyy[i] = ts_xxxy_xxyyy[i] * pa_z[i];

        ts_xxxyz_xxyyz[i] = 2.0 * ts_xxxz_xxyz[i] * fe_0 + ts_xxxz_xxyyz[i] * pa_y[i];

        ts_xxxyz_xxyzz[i] = ts_xxxz_xxzz[i] * fe_0 + ts_xxxz_xxyzz[i] * pa_y[i];

        ts_xxxyz_xxzzz[i] = ts_xxxz_xxzzz[i] * pa_y[i];

        ts_xxxyz_xyyyy[i] = ts_xxxy_xyyyy[i] * pa_z[i];

        ts_xxxyz_xyyyz[i] = 3.0 * ts_xxxz_xyyz[i] * fe_0 + ts_xxxz_xyyyz[i] * pa_y[i];

        ts_xxxyz_xyyzz[i] = 2.0 * ts_xxxz_xyzz[i] * fe_0 + ts_xxxz_xyyzz[i] * pa_y[i];

        ts_xxxyz_xyzzz[i] = ts_xxxz_xzzz[i] * fe_0 + ts_xxxz_xyzzz[i] * pa_y[i];

        ts_xxxyz_xzzzz[i] = ts_xxxz_xzzzz[i] * pa_y[i];

        ts_xxxyz_yyyyy[i] = ts_xxxy_yyyyy[i] * pa_z[i];

        ts_xxxyz_yyyyz[i] = 2.0 * ts_xyz_yyyyz[i] * fe_0 + ts_xxyz_yyyyz[i] * pa_x[i];

        ts_xxxyz_yyyzz[i] = 2.0 * ts_xyz_yyyzz[i] * fe_0 + ts_xxyz_yyyzz[i] * pa_x[i];

        ts_xxxyz_yyzzz[i] = 2.0 * ts_xyz_yyzzz[i] * fe_0 + ts_xxyz_yyzzz[i] * pa_x[i];

        ts_xxxyz_yzzzz[i] = 2.0 * ts_xyz_yzzzz[i] * fe_0 + ts_xxyz_yzzzz[i] * pa_x[i];

        ts_xxxyz_zzzzz[i] = ts_xxxz_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : HH

    auto ts_xxxzz_xxxxx = pbuffer.data(idx_ovl_hh + 105);

    auto ts_xxxzz_xxxxy = pbuffer.data(idx_ovl_hh + 106);

    auto ts_xxxzz_xxxxz = pbuffer.data(idx_ovl_hh + 107);

    auto ts_xxxzz_xxxyy = pbuffer.data(idx_ovl_hh + 108);

    auto ts_xxxzz_xxxyz = pbuffer.data(idx_ovl_hh + 109);

    auto ts_xxxzz_xxxzz = pbuffer.data(idx_ovl_hh + 110);

    auto ts_xxxzz_xxyyy = pbuffer.data(idx_ovl_hh + 111);

    auto ts_xxxzz_xxyyz = pbuffer.data(idx_ovl_hh + 112);

    auto ts_xxxzz_xxyzz = pbuffer.data(idx_ovl_hh + 113);

    auto ts_xxxzz_xxzzz = pbuffer.data(idx_ovl_hh + 114);

    auto ts_xxxzz_xyyyy = pbuffer.data(idx_ovl_hh + 115);

    auto ts_xxxzz_xyyyz = pbuffer.data(idx_ovl_hh + 116);

    auto ts_xxxzz_xyyzz = pbuffer.data(idx_ovl_hh + 117);

    auto ts_xxxzz_xyzzz = pbuffer.data(idx_ovl_hh + 118);

    auto ts_xxxzz_xzzzz = pbuffer.data(idx_ovl_hh + 119);

    auto ts_xxxzz_yyyyy = pbuffer.data(idx_ovl_hh + 120);

    auto ts_xxxzz_yyyyz = pbuffer.data(idx_ovl_hh + 121);

    auto ts_xxxzz_yyyzz = pbuffer.data(idx_ovl_hh + 122);

    auto ts_xxxzz_yyzzz = pbuffer.data(idx_ovl_hh + 123);

    auto ts_xxxzz_yzzzz = pbuffer.data(idx_ovl_hh + 124);

    auto ts_xxxzz_zzzzz = pbuffer.data(idx_ovl_hh + 125);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxyy, ts_xxx_xxyyy, ts_xxx_xyyyy, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxyy, ts_xxxz_xxyyy, ts_xxxz_xyyyy, ts_xxxzz_xxxxx, ts_xxxzz_xxxxy, ts_xxxzz_xxxxz, ts_xxxzz_xxxyy, ts_xxxzz_xxxyz, ts_xxxzz_xxxzz, ts_xxxzz_xxyyy, ts_xxxzz_xxyyz, ts_xxxzz_xxyzz, ts_xxxzz_xxzzz, ts_xxxzz_xyyyy, ts_xxxzz_xyyyz, ts_xxxzz_xyyzz, ts_xxxzz_xyzzz, ts_xxxzz_xzzzz, ts_xxxzz_yyyyy, ts_xxxzz_yyyyz, ts_xxxzz_yyyzz, ts_xxxzz_yyzzz, ts_xxxzz_yzzzz, ts_xxxzz_zzzzz, ts_xxzz_xxxxz, ts_xxzz_xxxyz, ts_xxzz_xxxz, ts_xxzz_xxxzz, ts_xxzz_xxyyz, ts_xxzz_xxyz, ts_xxzz_xxyzz, ts_xxzz_xxzz, ts_xxzz_xxzzz, ts_xxzz_xyyyz, ts_xxzz_xyyz, ts_xxzz_xyyzz, ts_xxzz_xyzz, ts_xxzz_xyzzz, ts_xxzz_xzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyz, ts_xxzz_yyyzz, ts_xxzz_yyzz, ts_xxzz_yyzzz, ts_xxzz_yzzz, ts_xxzz_yzzzz, ts_xxzz_zzzz, ts_xxzz_zzzzz, ts_xzz_xxxxz, ts_xzz_xxxyz, ts_xzz_xxxzz, ts_xzz_xxyyz, ts_xzz_xxyzz, ts_xzz_xxzzz, ts_xzz_xyyyz, ts_xzz_xyyzz, ts_xzz_xyzzz, ts_xzz_xzzzz, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyzz, ts_xzz_yyzzz, ts_xzz_yzzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzz_xxxxx[i] = ts_xxx_xxxxx[i] * fe_0 + ts_xxxz_xxxxx[i] * pa_z[i];

        ts_xxxzz_xxxxy[i] = ts_xxx_xxxxy[i] * fe_0 + ts_xxxz_xxxxy[i] * pa_z[i];

        ts_xxxzz_xxxxz[i] = 2.0 * ts_xzz_xxxxz[i] * fe_0 + 4.0 * ts_xxzz_xxxz[i] * fe_0 + ts_xxzz_xxxxz[i] * pa_x[i];

        ts_xxxzz_xxxyy[i] = ts_xxx_xxxyy[i] * fe_0 + ts_xxxz_xxxyy[i] * pa_z[i];

        ts_xxxzz_xxxyz[i] = 2.0 * ts_xzz_xxxyz[i] * fe_0 + 3.0 * ts_xxzz_xxyz[i] * fe_0 + ts_xxzz_xxxyz[i] * pa_x[i];

        ts_xxxzz_xxxzz[i] = 2.0 * ts_xzz_xxxzz[i] * fe_0 + 3.0 * ts_xxzz_xxzz[i] * fe_0 + ts_xxzz_xxxzz[i] * pa_x[i];

        ts_xxxzz_xxyyy[i] = ts_xxx_xxyyy[i] * fe_0 + ts_xxxz_xxyyy[i] * pa_z[i];

        ts_xxxzz_xxyyz[i] = 2.0 * ts_xzz_xxyyz[i] * fe_0 + 2.0 * ts_xxzz_xyyz[i] * fe_0 + ts_xxzz_xxyyz[i] * pa_x[i];

        ts_xxxzz_xxyzz[i] = 2.0 * ts_xzz_xxyzz[i] * fe_0 + 2.0 * ts_xxzz_xyzz[i] * fe_0 + ts_xxzz_xxyzz[i] * pa_x[i];

        ts_xxxzz_xxzzz[i] = 2.0 * ts_xzz_xxzzz[i] * fe_0 + 2.0 * ts_xxzz_xzzz[i] * fe_0 + ts_xxzz_xxzzz[i] * pa_x[i];

        ts_xxxzz_xyyyy[i] = ts_xxx_xyyyy[i] * fe_0 + ts_xxxz_xyyyy[i] * pa_z[i];

        ts_xxxzz_xyyyz[i] = 2.0 * ts_xzz_xyyyz[i] * fe_0 + ts_xxzz_yyyz[i] * fe_0 + ts_xxzz_xyyyz[i] * pa_x[i];

        ts_xxxzz_xyyzz[i] = 2.0 * ts_xzz_xyyzz[i] * fe_0 + ts_xxzz_yyzz[i] * fe_0 + ts_xxzz_xyyzz[i] * pa_x[i];

        ts_xxxzz_xyzzz[i] = 2.0 * ts_xzz_xyzzz[i] * fe_0 + ts_xxzz_yzzz[i] * fe_0 + ts_xxzz_xyzzz[i] * pa_x[i];

        ts_xxxzz_xzzzz[i] = 2.0 * ts_xzz_xzzzz[i] * fe_0 + ts_xxzz_zzzz[i] * fe_0 + ts_xxzz_xzzzz[i] * pa_x[i];

        ts_xxxzz_yyyyy[i] = 2.0 * ts_xzz_yyyyy[i] * fe_0 + ts_xxzz_yyyyy[i] * pa_x[i];

        ts_xxxzz_yyyyz[i] = 2.0 * ts_xzz_yyyyz[i] * fe_0 + ts_xxzz_yyyyz[i] * pa_x[i];

        ts_xxxzz_yyyzz[i] = 2.0 * ts_xzz_yyyzz[i] * fe_0 + ts_xxzz_yyyzz[i] * pa_x[i];

        ts_xxxzz_yyzzz[i] = 2.0 * ts_xzz_yyzzz[i] * fe_0 + ts_xxzz_yyzzz[i] * pa_x[i];

        ts_xxxzz_yzzzz[i] = 2.0 * ts_xzz_yzzzz[i] * fe_0 + ts_xxzz_yzzzz[i] * pa_x[i];

        ts_xxxzz_zzzzz[i] = 2.0 * ts_xzz_zzzzz[i] * fe_0 + ts_xxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : HH

    auto ts_xxyyy_xxxxx = pbuffer.data(idx_ovl_hh + 126);

    auto ts_xxyyy_xxxxy = pbuffer.data(idx_ovl_hh + 127);

    auto ts_xxyyy_xxxxz = pbuffer.data(idx_ovl_hh + 128);

    auto ts_xxyyy_xxxyy = pbuffer.data(idx_ovl_hh + 129);

    auto ts_xxyyy_xxxyz = pbuffer.data(idx_ovl_hh + 130);

    auto ts_xxyyy_xxxzz = pbuffer.data(idx_ovl_hh + 131);

    auto ts_xxyyy_xxyyy = pbuffer.data(idx_ovl_hh + 132);

    auto ts_xxyyy_xxyyz = pbuffer.data(idx_ovl_hh + 133);

    auto ts_xxyyy_xxyzz = pbuffer.data(idx_ovl_hh + 134);

    auto ts_xxyyy_xxzzz = pbuffer.data(idx_ovl_hh + 135);

    auto ts_xxyyy_xyyyy = pbuffer.data(idx_ovl_hh + 136);

    auto ts_xxyyy_xyyyz = pbuffer.data(idx_ovl_hh + 137);

    auto ts_xxyyy_xyyzz = pbuffer.data(idx_ovl_hh + 138);

    auto ts_xxyyy_xyzzz = pbuffer.data(idx_ovl_hh + 139);

    auto ts_xxyyy_xzzzz = pbuffer.data(idx_ovl_hh + 140);

    auto ts_xxyyy_yyyyy = pbuffer.data(idx_ovl_hh + 141);

    auto ts_xxyyy_yyyyz = pbuffer.data(idx_ovl_hh + 142);

    auto ts_xxyyy_yyyzz = pbuffer.data(idx_ovl_hh + 143);

    auto ts_xxyyy_yyzzz = pbuffer.data(idx_ovl_hh + 144);

    auto ts_xxyyy_yzzzz = pbuffer.data(idx_ovl_hh + 145);

    auto ts_xxyyy_zzzzz = pbuffer.data(idx_ovl_hh + 146);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxy_xxxxx, ts_xxy_xxxxz, ts_xxy_xxxzz, ts_xxy_xxzzz, ts_xxy_xzzzz, ts_xxyy_xxxxx, ts_xxyy_xxxxz, ts_xxyy_xxxzz, ts_xxyy_xxzzz, ts_xxyy_xzzzz, ts_xxyyy_xxxxx, ts_xxyyy_xxxxy, ts_xxyyy_xxxxz, ts_xxyyy_xxxyy, ts_xxyyy_xxxyz, ts_xxyyy_xxxzz, ts_xxyyy_xxyyy, ts_xxyyy_xxyyz, ts_xxyyy_xxyzz, ts_xxyyy_xxzzz, ts_xxyyy_xyyyy, ts_xxyyy_xyyyz, ts_xxyyy_xyyzz, ts_xxyyy_xyzzz, ts_xxyyy_xzzzz, ts_xxyyy_yyyyy, ts_xxyyy_yyyyz, ts_xxyyy_yyyzz, ts_xxyyy_yyzzz, ts_xxyyy_yzzzz, ts_xxyyy_zzzzz, ts_xyyy_xxxxy, ts_xyyy_xxxy, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxyy, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyz, ts_xyyy_xxyzz, ts_xyyy_xyyy, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyz, ts_xyyy_xyyzz, ts_xyyy_xyzz, ts_xyyy_xyzzz, ts_xyyy_yyyy, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyz, ts_xyyy_yyyzz, ts_xyyy_yyzz, ts_xyyy_yyzzz, ts_xyyy_yzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, ts_yyy_xxxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyy_xxxxx[i] = 2.0 * ts_xxy_xxxxx[i] * fe_0 + ts_xxyy_xxxxx[i] * pa_y[i];

        ts_xxyyy_xxxxy[i] = ts_yyy_xxxxy[i] * fe_0 + 4.0 * ts_xyyy_xxxy[i] * fe_0 + ts_xyyy_xxxxy[i] * pa_x[i];

        ts_xxyyy_xxxxz[i] = 2.0 * ts_xxy_xxxxz[i] * fe_0 + ts_xxyy_xxxxz[i] * pa_y[i];

        ts_xxyyy_xxxyy[i] = ts_yyy_xxxyy[i] * fe_0 + 3.0 * ts_xyyy_xxyy[i] * fe_0 + ts_xyyy_xxxyy[i] * pa_x[i];

        ts_xxyyy_xxxyz[i] = ts_yyy_xxxyz[i] * fe_0 + 3.0 * ts_xyyy_xxyz[i] * fe_0 + ts_xyyy_xxxyz[i] * pa_x[i];

        ts_xxyyy_xxxzz[i] = 2.0 * ts_xxy_xxxzz[i] * fe_0 + ts_xxyy_xxxzz[i] * pa_y[i];

        ts_xxyyy_xxyyy[i] = ts_yyy_xxyyy[i] * fe_0 + 2.0 * ts_xyyy_xyyy[i] * fe_0 + ts_xyyy_xxyyy[i] * pa_x[i];

        ts_xxyyy_xxyyz[i] = ts_yyy_xxyyz[i] * fe_0 + 2.0 * ts_xyyy_xyyz[i] * fe_0 + ts_xyyy_xxyyz[i] * pa_x[i];

        ts_xxyyy_xxyzz[i] = ts_yyy_xxyzz[i] * fe_0 + 2.0 * ts_xyyy_xyzz[i] * fe_0 + ts_xyyy_xxyzz[i] * pa_x[i];

        ts_xxyyy_xxzzz[i] = 2.0 * ts_xxy_xxzzz[i] * fe_0 + ts_xxyy_xxzzz[i] * pa_y[i];

        ts_xxyyy_xyyyy[i] = ts_yyy_xyyyy[i] * fe_0 + ts_xyyy_yyyy[i] * fe_0 + ts_xyyy_xyyyy[i] * pa_x[i];

        ts_xxyyy_xyyyz[i] = ts_yyy_xyyyz[i] * fe_0 + ts_xyyy_yyyz[i] * fe_0 + ts_xyyy_xyyyz[i] * pa_x[i];

        ts_xxyyy_xyyzz[i] = ts_yyy_xyyzz[i] * fe_0 + ts_xyyy_yyzz[i] * fe_0 + ts_xyyy_xyyzz[i] * pa_x[i];

        ts_xxyyy_xyzzz[i] = ts_yyy_xyzzz[i] * fe_0 + ts_xyyy_yzzz[i] * fe_0 + ts_xyyy_xyzzz[i] * pa_x[i];

        ts_xxyyy_xzzzz[i] = 2.0 * ts_xxy_xzzzz[i] * fe_0 + ts_xxyy_xzzzz[i] * pa_y[i];

        ts_xxyyy_yyyyy[i] = ts_yyy_yyyyy[i] * fe_0 + ts_xyyy_yyyyy[i] * pa_x[i];

        ts_xxyyy_yyyyz[i] = ts_yyy_yyyyz[i] * fe_0 + ts_xyyy_yyyyz[i] * pa_x[i];

        ts_xxyyy_yyyzz[i] = ts_yyy_yyyzz[i] * fe_0 + ts_xyyy_yyyzz[i] * pa_x[i];

        ts_xxyyy_yyzzz[i] = ts_yyy_yyzzz[i] * fe_0 + ts_xyyy_yyzzz[i] * pa_x[i];

        ts_xxyyy_yzzzz[i] = ts_yyy_yzzzz[i] * fe_0 + ts_xyyy_yzzzz[i] * pa_x[i];

        ts_xxyyy_zzzzz[i] = ts_yyy_zzzzz[i] * fe_0 + ts_xyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : HH

    auto ts_xxyyz_xxxxx = pbuffer.data(idx_ovl_hh + 147);

    auto ts_xxyyz_xxxxy = pbuffer.data(idx_ovl_hh + 148);

    auto ts_xxyyz_xxxxz = pbuffer.data(idx_ovl_hh + 149);

    auto ts_xxyyz_xxxyy = pbuffer.data(idx_ovl_hh + 150);

    auto ts_xxyyz_xxxyz = pbuffer.data(idx_ovl_hh + 151);

    auto ts_xxyyz_xxxzz = pbuffer.data(idx_ovl_hh + 152);

    auto ts_xxyyz_xxyyy = pbuffer.data(idx_ovl_hh + 153);

    auto ts_xxyyz_xxyyz = pbuffer.data(idx_ovl_hh + 154);

    auto ts_xxyyz_xxyzz = pbuffer.data(idx_ovl_hh + 155);

    auto ts_xxyyz_xxzzz = pbuffer.data(idx_ovl_hh + 156);

    auto ts_xxyyz_xyyyy = pbuffer.data(idx_ovl_hh + 157);

    auto ts_xxyyz_xyyyz = pbuffer.data(idx_ovl_hh + 158);

    auto ts_xxyyz_xyyzz = pbuffer.data(idx_ovl_hh + 159);

    auto ts_xxyyz_xyzzz = pbuffer.data(idx_ovl_hh + 160);

    auto ts_xxyyz_xzzzz = pbuffer.data(idx_ovl_hh + 161);

    auto ts_xxyyz_yyyyy = pbuffer.data(idx_ovl_hh + 162);

    auto ts_xxyyz_yyyyz = pbuffer.data(idx_ovl_hh + 163);

    auto ts_xxyyz_yyyzz = pbuffer.data(idx_ovl_hh + 164);

    auto ts_xxyyz_yyzzz = pbuffer.data(idx_ovl_hh + 165);

    auto ts_xxyyz_yzzzz = pbuffer.data(idx_ovl_hh + 166);

    auto ts_xxyyz_zzzzz = pbuffer.data(idx_ovl_hh + 167);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxyy, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyz, ts_xxyy_xxyzz, ts_xxyy_xyyy, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyz, ts_xxyy_xyyzz, ts_xxyy_xyzz, ts_xxyy_xyzzz, ts_xxyy_yyyyy, ts_xxyyz_xxxxx, ts_xxyyz_xxxxy, ts_xxyyz_xxxxz, ts_xxyyz_xxxyy, ts_xxyyz_xxxyz, ts_xxyyz_xxxzz, ts_xxyyz_xxyyy, ts_xxyyz_xxyyz, ts_xxyyz_xxyzz, ts_xxyyz_xxzzz, ts_xxyyz_xyyyy, ts_xxyyz_xyyyz, ts_xxyyz_xyyzz, ts_xxyyz_xyzzz, ts_xxyyz_xzzzz, ts_xxyyz_yyyyy, ts_xxyyz_yyyyz, ts_xxyyz_yyyzz, ts_xxyyz_yyzzz, ts_xxyyz_yzzzz, ts_xxyyz_zzzzz, ts_xxyz_xxxxz, ts_xxyz_xxxzz, ts_xxyz_xxzzz, ts_xxyz_xzzzz, ts_xxz_xxxxz, ts_xxz_xxxzz, ts_xxz_xxzzz, ts_xxz_xzzzz, ts_xyyz_yyyyz, ts_xyyz_yyyzz, ts_xyyz_yyzzz, ts_xyyz_yzzzz, ts_xyyz_zzzzz, ts_yyz_yyyyz, ts_yyz_yyyzz, ts_yyz_yyzzz, ts_yyz_yzzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyz_xxxxx[i] = ts_xxyy_xxxxx[i] * pa_z[i];

        ts_xxyyz_xxxxy[i] = ts_xxyy_xxxxy[i] * pa_z[i];

        ts_xxyyz_xxxxz[i] = ts_xxz_xxxxz[i] * fe_0 + ts_xxyz_xxxxz[i] * pa_y[i];

        ts_xxyyz_xxxyy[i] = ts_xxyy_xxxyy[i] * pa_z[i];

        ts_xxyyz_xxxyz[i] = ts_xxyy_xxxy[i] * fe_0 + ts_xxyy_xxxyz[i] * pa_z[i];

        ts_xxyyz_xxxzz[i] = ts_xxz_xxxzz[i] * fe_0 + ts_xxyz_xxxzz[i] * pa_y[i];

        ts_xxyyz_xxyyy[i] = ts_xxyy_xxyyy[i] * pa_z[i];

        ts_xxyyz_xxyyz[i] = ts_xxyy_xxyy[i] * fe_0 + ts_xxyy_xxyyz[i] * pa_z[i];

        ts_xxyyz_xxyzz[i] = 2.0 * ts_xxyy_xxyz[i] * fe_0 + ts_xxyy_xxyzz[i] * pa_z[i];

        ts_xxyyz_xxzzz[i] = ts_xxz_xxzzz[i] * fe_0 + ts_xxyz_xxzzz[i] * pa_y[i];

        ts_xxyyz_xyyyy[i] = ts_xxyy_xyyyy[i] * pa_z[i];

        ts_xxyyz_xyyyz[i] = ts_xxyy_xyyy[i] * fe_0 + ts_xxyy_xyyyz[i] * pa_z[i];

        ts_xxyyz_xyyzz[i] = 2.0 * ts_xxyy_xyyz[i] * fe_0 + ts_xxyy_xyyzz[i] * pa_z[i];

        ts_xxyyz_xyzzz[i] = 3.0 * ts_xxyy_xyzz[i] * fe_0 + ts_xxyy_xyzzz[i] * pa_z[i];

        ts_xxyyz_xzzzz[i] = ts_xxz_xzzzz[i] * fe_0 + ts_xxyz_xzzzz[i] * pa_y[i];

        ts_xxyyz_yyyyy[i] = ts_xxyy_yyyyy[i] * pa_z[i];

        ts_xxyyz_yyyyz[i] = ts_yyz_yyyyz[i] * fe_0 + ts_xyyz_yyyyz[i] * pa_x[i];

        ts_xxyyz_yyyzz[i] = ts_yyz_yyyzz[i] * fe_0 + ts_xyyz_yyyzz[i] * pa_x[i];

        ts_xxyyz_yyzzz[i] = ts_yyz_yyzzz[i] * fe_0 + ts_xyyz_yyzzz[i] * pa_x[i];

        ts_xxyyz_yzzzz[i] = ts_yyz_yzzzz[i] * fe_0 + ts_xyyz_yzzzz[i] * pa_x[i];

        ts_xxyyz_zzzzz[i] = ts_yyz_zzzzz[i] * fe_0 + ts_xyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 168-189 components of targeted buffer : HH

    auto ts_xxyzz_xxxxx = pbuffer.data(idx_ovl_hh + 168);

    auto ts_xxyzz_xxxxy = pbuffer.data(idx_ovl_hh + 169);

    auto ts_xxyzz_xxxxz = pbuffer.data(idx_ovl_hh + 170);

    auto ts_xxyzz_xxxyy = pbuffer.data(idx_ovl_hh + 171);

    auto ts_xxyzz_xxxyz = pbuffer.data(idx_ovl_hh + 172);

    auto ts_xxyzz_xxxzz = pbuffer.data(idx_ovl_hh + 173);

    auto ts_xxyzz_xxyyy = pbuffer.data(idx_ovl_hh + 174);

    auto ts_xxyzz_xxyyz = pbuffer.data(idx_ovl_hh + 175);

    auto ts_xxyzz_xxyzz = pbuffer.data(idx_ovl_hh + 176);

    auto ts_xxyzz_xxzzz = pbuffer.data(idx_ovl_hh + 177);

    auto ts_xxyzz_xyyyy = pbuffer.data(idx_ovl_hh + 178);

    auto ts_xxyzz_xyyyz = pbuffer.data(idx_ovl_hh + 179);

    auto ts_xxyzz_xyyzz = pbuffer.data(idx_ovl_hh + 180);

    auto ts_xxyzz_xyzzz = pbuffer.data(idx_ovl_hh + 181);

    auto ts_xxyzz_xzzzz = pbuffer.data(idx_ovl_hh + 182);

    auto ts_xxyzz_yyyyy = pbuffer.data(idx_ovl_hh + 183);

    auto ts_xxyzz_yyyyz = pbuffer.data(idx_ovl_hh + 184);

    auto ts_xxyzz_yyyzz = pbuffer.data(idx_ovl_hh + 185);

    auto ts_xxyzz_yyzzz = pbuffer.data(idx_ovl_hh + 186);

    auto ts_xxyzz_yzzzz = pbuffer.data(idx_ovl_hh + 187);

    auto ts_xxyzz_zzzzz = pbuffer.data(idx_ovl_hh + 188);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzz_xxxxx, ts_xxyzz_xxxxy, ts_xxyzz_xxxxz, ts_xxyzz_xxxyy, ts_xxyzz_xxxyz, ts_xxyzz_xxxzz, ts_xxyzz_xxyyy, ts_xxyzz_xxyyz, ts_xxyzz_xxyzz, ts_xxyzz_xxzzz, ts_xxyzz_xyyyy, ts_xxyzz_xyyyz, ts_xxyzz_xyyzz, ts_xxyzz_xyzzz, ts_xxyzz_xzzzz, ts_xxyzz_yyyyy, ts_xxyzz_yyyyz, ts_xxyzz_yyyzz, ts_xxyzz_yyzzz, ts_xxyzz_yzzzz, ts_xxyzz_zzzzz, ts_xxzz_xxxx, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxy, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxz, ts_xxzz_xxxzz, ts_xxzz_xxyy, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyz, ts_xxzz_xxyzz, ts_xxzz_xxzz, ts_xxzz_xxzzz, ts_xxzz_xyyy, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyz, ts_xxzz_xyyzz, ts_xxzz_xyzz, ts_xxzz_xyzzz, ts_xxzz_xzzz, ts_xxzz_xzzzz, ts_xxzz_zzzzz, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyzz, ts_xyzz_yyzzz, ts_xyzz_yzzzz, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzz_xxxxx[i] = ts_xxzz_xxxxx[i] * pa_y[i];

        ts_xxyzz_xxxxy[i] = ts_xxzz_xxxx[i] * fe_0 + ts_xxzz_xxxxy[i] * pa_y[i];

        ts_xxyzz_xxxxz[i] = ts_xxzz_xxxxz[i] * pa_y[i];

        ts_xxyzz_xxxyy[i] = 2.0 * ts_xxzz_xxxy[i] * fe_0 + ts_xxzz_xxxyy[i] * pa_y[i];

        ts_xxyzz_xxxyz[i] = ts_xxzz_xxxz[i] * fe_0 + ts_xxzz_xxxyz[i] * pa_y[i];

        ts_xxyzz_xxxzz[i] = ts_xxzz_xxxzz[i] * pa_y[i];

        ts_xxyzz_xxyyy[i] = 3.0 * ts_xxzz_xxyy[i] * fe_0 + ts_xxzz_xxyyy[i] * pa_y[i];

        ts_xxyzz_xxyyz[i] = 2.0 * ts_xxzz_xxyz[i] * fe_0 + ts_xxzz_xxyyz[i] * pa_y[i];

        ts_xxyzz_xxyzz[i] = ts_xxzz_xxzz[i] * fe_0 + ts_xxzz_xxyzz[i] * pa_y[i];

        ts_xxyzz_xxzzz[i] = ts_xxzz_xxzzz[i] * pa_y[i];

        ts_xxyzz_xyyyy[i] = 4.0 * ts_xxzz_xyyy[i] * fe_0 + ts_xxzz_xyyyy[i] * pa_y[i];

        ts_xxyzz_xyyyz[i] = 3.0 * ts_xxzz_xyyz[i] * fe_0 + ts_xxzz_xyyyz[i] * pa_y[i];

        ts_xxyzz_xyyzz[i] = 2.0 * ts_xxzz_xyzz[i] * fe_0 + ts_xxzz_xyyzz[i] * pa_y[i];

        ts_xxyzz_xyzzz[i] = ts_xxzz_xzzz[i] * fe_0 + ts_xxzz_xyzzz[i] * pa_y[i];

        ts_xxyzz_xzzzz[i] = ts_xxzz_xzzzz[i] * pa_y[i];

        ts_xxyzz_yyyyy[i] = ts_yzz_yyyyy[i] * fe_0 + ts_xyzz_yyyyy[i] * pa_x[i];

        ts_xxyzz_yyyyz[i] = ts_yzz_yyyyz[i] * fe_0 + ts_xyzz_yyyyz[i] * pa_x[i];

        ts_xxyzz_yyyzz[i] = ts_yzz_yyyzz[i] * fe_0 + ts_xyzz_yyyzz[i] * pa_x[i];

        ts_xxyzz_yyzzz[i] = ts_yzz_yyzzz[i] * fe_0 + ts_xyzz_yyzzz[i] * pa_x[i];

        ts_xxyzz_yzzzz[i] = ts_yzz_yzzzz[i] * fe_0 + ts_xyzz_yzzzz[i] * pa_x[i];

        ts_xxyzz_zzzzz[i] = ts_xxzz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : HH

    auto ts_xxzzz_xxxxx = pbuffer.data(idx_ovl_hh + 189);

    auto ts_xxzzz_xxxxy = pbuffer.data(idx_ovl_hh + 190);

    auto ts_xxzzz_xxxxz = pbuffer.data(idx_ovl_hh + 191);

    auto ts_xxzzz_xxxyy = pbuffer.data(idx_ovl_hh + 192);

    auto ts_xxzzz_xxxyz = pbuffer.data(idx_ovl_hh + 193);

    auto ts_xxzzz_xxxzz = pbuffer.data(idx_ovl_hh + 194);

    auto ts_xxzzz_xxyyy = pbuffer.data(idx_ovl_hh + 195);

    auto ts_xxzzz_xxyyz = pbuffer.data(idx_ovl_hh + 196);

    auto ts_xxzzz_xxyzz = pbuffer.data(idx_ovl_hh + 197);

    auto ts_xxzzz_xxzzz = pbuffer.data(idx_ovl_hh + 198);

    auto ts_xxzzz_xyyyy = pbuffer.data(idx_ovl_hh + 199);

    auto ts_xxzzz_xyyyz = pbuffer.data(idx_ovl_hh + 200);

    auto ts_xxzzz_xyyzz = pbuffer.data(idx_ovl_hh + 201);

    auto ts_xxzzz_xyzzz = pbuffer.data(idx_ovl_hh + 202);

    auto ts_xxzzz_xzzzz = pbuffer.data(idx_ovl_hh + 203);

    auto ts_xxzzz_yyyyy = pbuffer.data(idx_ovl_hh + 204);

    auto ts_xxzzz_yyyyz = pbuffer.data(idx_ovl_hh + 205);

    auto ts_xxzzz_yyyzz = pbuffer.data(idx_ovl_hh + 206);

    auto ts_xxzzz_yyzzz = pbuffer.data(idx_ovl_hh + 207);

    auto ts_xxzzz_yzzzz = pbuffer.data(idx_ovl_hh + 208);

    auto ts_xxzzz_zzzzz = pbuffer.data(idx_ovl_hh + 209);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxyy, ts_xxz_xxyyy, ts_xxz_xyyyy, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxyy, ts_xxzz_xxyyy, ts_xxzz_xyyyy, ts_xxzzz_xxxxx, ts_xxzzz_xxxxy, ts_xxzzz_xxxxz, ts_xxzzz_xxxyy, ts_xxzzz_xxxyz, ts_xxzzz_xxxzz, ts_xxzzz_xxyyy, ts_xxzzz_xxyyz, ts_xxzzz_xxyzz, ts_xxzzz_xxzzz, ts_xxzzz_xyyyy, ts_xxzzz_xyyyz, ts_xxzzz_xyyzz, ts_xxzzz_xyzzz, ts_xxzzz_xzzzz, ts_xxzzz_yyyyy, ts_xxzzz_yyyyz, ts_xxzzz_yyyzz, ts_xxzzz_yyzzz, ts_xxzzz_yzzzz, ts_xxzzz_zzzzz, ts_xzzz_xxxxz, ts_xzzz_xxxyz, ts_xzzz_xxxz, ts_xzzz_xxxzz, ts_xzzz_xxyyz, ts_xzzz_xxyz, ts_xzzz_xxyzz, ts_xzzz_xxzz, ts_xzzz_xxzzz, ts_xzzz_xyyyz, ts_xzzz_xyyz, ts_xzzz_xyyzz, ts_xzzz_xyzz, ts_xzzz_xyzzz, ts_xzzz_xzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyz, ts_xzzz_yyyzz, ts_xzzz_yyzz, ts_xzzz_yyzzz, ts_xzzz_yzzz, ts_xzzz_yzzzz, ts_xzzz_zzzz, ts_xzzz_zzzzz, ts_zzz_xxxxz, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzz_xxxxx[i] = 2.0 * ts_xxz_xxxxx[i] * fe_0 + ts_xxzz_xxxxx[i] * pa_z[i];

        ts_xxzzz_xxxxy[i] = 2.0 * ts_xxz_xxxxy[i] * fe_0 + ts_xxzz_xxxxy[i] * pa_z[i];

        ts_xxzzz_xxxxz[i] = ts_zzz_xxxxz[i] * fe_0 + 4.0 * ts_xzzz_xxxz[i] * fe_0 + ts_xzzz_xxxxz[i] * pa_x[i];

        ts_xxzzz_xxxyy[i] = 2.0 * ts_xxz_xxxyy[i] * fe_0 + ts_xxzz_xxxyy[i] * pa_z[i];

        ts_xxzzz_xxxyz[i] = ts_zzz_xxxyz[i] * fe_0 + 3.0 * ts_xzzz_xxyz[i] * fe_0 + ts_xzzz_xxxyz[i] * pa_x[i];

        ts_xxzzz_xxxzz[i] = ts_zzz_xxxzz[i] * fe_0 + 3.0 * ts_xzzz_xxzz[i] * fe_0 + ts_xzzz_xxxzz[i] * pa_x[i];

        ts_xxzzz_xxyyy[i] = 2.0 * ts_xxz_xxyyy[i] * fe_0 + ts_xxzz_xxyyy[i] * pa_z[i];

        ts_xxzzz_xxyyz[i] = ts_zzz_xxyyz[i] * fe_0 + 2.0 * ts_xzzz_xyyz[i] * fe_0 + ts_xzzz_xxyyz[i] * pa_x[i];

        ts_xxzzz_xxyzz[i] = ts_zzz_xxyzz[i] * fe_0 + 2.0 * ts_xzzz_xyzz[i] * fe_0 + ts_xzzz_xxyzz[i] * pa_x[i];

        ts_xxzzz_xxzzz[i] = ts_zzz_xxzzz[i] * fe_0 + 2.0 * ts_xzzz_xzzz[i] * fe_0 + ts_xzzz_xxzzz[i] * pa_x[i];

        ts_xxzzz_xyyyy[i] = 2.0 * ts_xxz_xyyyy[i] * fe_0 + ts_xxzz_xyyyy[i] * pa_z[i];

        ts_xxzzz_xyyyz[i] = ts_zzz_xyyyz[i] * fe_0 + ts_xzzz_yyyz[i] * fe_0 + ts_xzzz_xyyyz[i] * pa_x[i];

        ts_xxzzz_xyyzz[i] = ts_zzz_xyyzz[i] * fe_0 + ts_xzzz_yyzz[i] * fe_0 + ts_xzzz_xyyzz[i] * pa_x[i];

        ts_xxzzz_xyzzz[i] = ts_zzz_xyzzz[i] * fe_0 + ts_xzzz_yzzz[i] * fe_0 + ts_xzzz_xyzzz[i] * pa_x[i];

        ts_xxzzz_xzzzz[i] = ts_zzz_xzzzz[i] * fe_0 + ts_xzzz_zzzz[i] * fe_0 + ts_xzzz_xzzzz[i] * pa_x[i];

        ts_xxzzz_yyyyy[i] = ts_zzz_yyyyy[i] * fe_0 + ts_xzzz_yyyyy[i] * pa_x[i];

        ts_xxzzz_yyyyz[i] = ts_zzz_yyyyz[i] * fe_0 + ts_xzzz_yyyyz[i] * pa_x[i];

        ts_xxzzz_yyyzz[i] = ts_zzz_yyyzz[i] * fe_0 + ts_xzzz_yyyzz[i] * pa_x[i];

        ts_xxzzz_yyzzz[i] = ts_zzz_yyzzz[i] * fe_0 + ts_xzzz_yyzzz[i] * pa_x[i];

        ts_xxzzz_yzzzz[i] = ts_zzz_yzzzz[i] * fe_0 + ts_xzzz_yzzzz[i] * pa_x[i];

        ts_xxzzz_zzzzz[i] = ts_zzz_zzzzz[i] * fe_0 + ts_xzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : HH

    auto ts_xyyyy_xxxxx = pbuffer.data(idx_ovl_hh + 210);

    auto ts_xyyyy_xxxxy = pbuffer.data(idx_ovl_hh + 211);

    auto ts_xyyyy_xxxxz = pbuffer.data(idx_ovl_hh + 212);

    auto ts_xyyyy_xxxyy = pbuffer.data(idx_ovl_hh + 213);

    auto ts_xyyyy_xxxyz = pbuffer.data(idx_ovl_hh + 214);

    auto ts_xyyyy_xxxzz = pbuffer.data(idx_ovl_hh + 215);

    auto ts_xyyyy_xxyyy = pbuffer.data(idx_ovl_hh + 216);

    auto ts_xyyyy_xxyyz = pbuffer.data(idx_ovl_hh + 217);

    auto ts_xyyyy_xxyzz = pbuffer.data(idx_ovl_hh + 218);

    auto ts_xyyyy_xxzzz = pbuffer.data(idx_ovl_hh + 219);

    auto ts_xyyyy_xyyyy = pbuffer.data(idx_ovl_hh + 220);

    auto ts_xyyyy_xyyyz = pbuffer.data(idx_ovl_hh + 221);

    auto ts_xyyyy_xyyzz = pbuffer.data(idx_ovl_hh + 222);

    auto ts_xyyyy_xyzzz = pbuffer.data(idx_ovl_hh + 223);

    auto ts_xyyyy_xzzzz = pbuffer.data(idx_ovl_hh + 224);

    auto ts_xyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 225);

    auto ts_xyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 226);

    auto ts_xyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 227);

    auto ts_xyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 228);

    auto ts_xyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 229);

    auto ts_xyyyy_zzzzz = pbuffer.data(idx_ovl_hh + 230);

    #pragma omp simd aligned(pa_x, ts_xyyyy_xxxxx, ts_xyyyy_xxxxy, ts_xyyyy_xxxxz, ts_xyyyy_xxxyy, ts_xyyyy_xxxyz, ts_xyyyy_xxxzz, ts_xyyyy_xxyyy, ts_xyyyy_xxyyz, ts_xyyyy_xxyzz, ts_xyyyy_xxzzz, ts_xyyyy_xyyyy, ts_xyyyy_xyyyz, ts_xyyyy_xyyzz, ts_xyyyy_xyzzz, ts_xyyyy_xzzzz, ts_xyyyy_yyyyy, ts_xyyyy_yyyyz, ts_xyyyy_yyyzz, ts_xyyyy_yyzzz, ts_xyyyy_yzzzz, ts_xyyyy_zzzzz, ts_yyyy_xxxx, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxz, ts_yyyy_xxxzz, ts_yyyy_xxyy, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyz, ts_yyyy_xxyzz, ts_yyyy_xxzz, ts_yyyy_xxzzz, ts_yyyy_xyyy, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyz, ts_yyyy_xyyzz, ts_yyyy_xyzz, ts_yyyy_xyzzz, ts_yyyy_xzzz, ts_yyyy_xzzzz, ts_yyyy_yyyy, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyz, ts_yyyy_yyyzz, ts_yyyy_yyzz, ts_yyyy_yyzzz, ts_yyyy_yzzz, ts_yyyy_yzzzz, ts_yyyy_zzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyy_xxxxx[i] = 5.0 * ts_yyyy_xxxx[i] * fe_0 + ts_yyyy_xxxxx[i] * pa_x[i];

        ts_xyyyy_xxxxy[i] = 4.0 * ts_yyyy_xxxy[i] * fe_0 + ts_yyyy_xxxxy[i] * pa_x[i];

        ts_xyyyy_xxxxz[i] = 4.0 * ts_yyyy_xxxz[i] * fe_0 + ts_yyyy_xxxxz[i] * pa_x[i];

        ts_xyyyy_xxxyy[i] = 3.0 * ts_yyyy_xxyy[i] * fe_0 + ts_yyyy_xxxyy[i] * pa_x[i];

        ts_xyyyy_xxxyz[i] = 3.0 * ts_yyyy_xxyz[i] * fe_0 + ts_yyyy_xxxyz[i] * pa_x[i];

        ts_xyyyy_xxxzz[i] = 3.0 * ts_yyyy_xxzz[i] * fe_0 + ts_yyyy_xxxzz[i] * pa_x[i];

        ts_xyyyy_xxyyy[i] = 2.0 * ts_yyyy_xyyy[i] * fe_0 + ts_yyyy_xxyyy[i] * pa_x[i];

        ts_xyyyy_xxyyz[i] = 2.0 * ts_yyyy_xyyz[i] * fe_0 + ts_yyyy_xxyyz[i] * pa_x[i];

        ts_xyyyy_xxyzz[i] = 2.0 * ts_yyyy_xyzz[i] * fe_0 + ts_yyyy_xxyzz[i] * pa_x[i];

        ts_xyyyy_xxzzz[i] = 2.0 * ts_yyyy_xzzz[i] * fe_0 + ts_yyyy_xxzzz[i] * pa_x[i];

        ts_xyyyy_xyyyy[i] = ts_yyyy_yyyy[i] * fe_0 + ts_yyyy_xyyyy[i] * pa_x[i];

        ts_xyyyy_xyyyz[i] = ts_yyyy_yyyz[i] * fe_0 + ts_yyyy_xyyyz[i] * pa_x[i];

        ts_xyyyy_xyyzz[i] = ts_yyyy_yyzz[i] * fe_0 + ts_yyyy_xyyzz[i] * pa_x[i];

        ts_xyyyy_xyzzz[i] = ts_yyyy_yzzz[i] * fe_0 + ts_yyyy_xyzzz[i] * pa_x[i];

        ts_xyyyy_xzzzz[i] = ts_yyyy_zzzz[i] * fe_0 + ts_yyyy_xzzzz[i] * pa_x[i];

        ts_xyyyy_yyyyy[i] = ts_yyyy_yyyyy[i] * pa_x[i];

        ts_xyyyy_yyyyz[i] = ts_yyyy_yyyyz[i] * pa_x[i];

        ts_xyyyy_yyyzz[i] = ts_yyyy_yyyzz[i] * pa_x[i];

        ts_xyyyy_yyzzz[i] = ts_yyyy_yyzzz[i] * pa_x[i];

        ts_xyyyy_yzzzz[i] = ts_yyyy_yzzzz[i] * pa_x[i];

        ts_xyyyy_zzzzz[i] = ts_yyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : HH

    auto ts_xyyyz_xxxxx = pbuffer.data(idx_ovl_hh + 231);

    auto ts_xyyyz_xxxxy = pbuffer.data(idx_ovl_hh + 232);

    auto ts_xyyyz_xxxxz = pbuffer.data(idx_ovl_hh + 233);

    auto ts_xyyyz_xxxyy = pbuffer.data(idx_ovl_hh + 234);

    auto ts_xyyyz_xxxyz = pbuffer.data(idx_ovl_hh + 235);

    auto ts_xyyyz_xxxzz = pbuffer.data(idx_ovl_hh + 236);

    auto ts_xyyyz_xxyyy = pbuffer.data(idx_ovl_hh + 237);

    auto ts_xyyyz_xxyyz = pbuffer.data(idx_ovl_hh + 238);

    auto ts_xyyyz_xxyzz = pbuffer.data(idx_ovl_hh + 239);

    auto ts_xyyyz_xxzzz = pbuffer.data(idx_ovl_hh + 240);

    auto ts_xyyyz_xyyyy = pbuffer.data(idx_ovl_hh + 241);

    auto ts_xyyyz_xyyyz = pbuffer.data(idx_ovl_hh + 242);

    auto ts_xyyyz_xyyzz = pbuffer.data(idx_ovl_hh + 243);

    auto ts_xyyyz_xyzzz = pbuffer.data(idx_ovl_hh + 244);

    auto ts_xyyyz_xzzzz = pbuffer.data(idx_ovl_hh + 245);

    auto ts_xyyyz_yyyyy = pbuffer.data(idx_ovl_hh + 246);

    auto ts_xyyyz_yyyyz = pbuffer.data(idx_ovl_hh + 247);

    auto ts_xyyyz_yyyzz = pbuffer.data(idx_ovl_hh + 248);

    auto ts_xyyyz_yyzzz = pbuffer.data(idx_ovl_hh + 249);

    auto ts_xyyyz_yzzzz = pbuffer.data(idx_ovl_hh + 250);

    auto ts_xyyyz_zzzzz = pbuffer.data(idx_ovl_hh + 251);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxyy, ts_xyyy_xxyyy, ts_xyyy_xyyyy, ts_xyyyz_xxxxx, ts_xyyyz_xxxxy, ts_xyyyz_xxxxz, ts_xyyyz_xxxyy, ts_xyyyz_xxxyz, ts_xyyyz_xxxzz, ts_xyyyz_xxyyy, ts_xyyyz_xxyyz, ts_xyyyz_xxyzz, ts_xyyyz_xxzzz, ts_xyyyz_xyyyy, ts_xyyyz_xyyyz, ts_xyyyz_xyyzz, ts_xyyyz_xyzzz, ts_xyyyz_xzzzz, ts_xyyyz_yyyyy, ts_xyyyz_yyyyz, ts_xyyyz_yyyzz, ts_xyyyz_yyzzz, ts_xyyyz_yzzzz, ts_xyyyz_zzzzz, ts_yyyz_xxxxz, ts_yyyz_xxxyz, ts_yyyz_xxxz, ts_yyyz_xxxzz, ts_yyyz_xxyyz, ts_yyyz_xxyz, ts_yyyz_xxyzz, ts_yyyz_xxzz, ts_yyyz_xxzzz, ts_yyyz_xyyyz, ts_yyyz_xyyz, ts_yyyz_xyyzz, ts_yyyz_xyzz, ts_yyyz_xyzzz, ts_yyyz_xzzz, ts_yyyz_xzzzz, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyz, ts_yyyz_yyyzz, ts_yyyz_yyzz, ts_yyyz_yyzzz, ts_yyyz_yzzz, ts_yyyz_yzzzz, ts_yyyz_zzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyz_xxxxx[i] = ts_xyyy_xxxxx[i] * pa_z[i];

        ts_xyyyz_xxxxy[i] = ts_xyyy_xxxxy[i] * pa_z[i];

        ts_xyyyz_xxxxz[i] = 4.0 * ts_yyyz_xxxz[i] * fe_0 + ts_yyyz_xxxxz[i] * pa_x[i];

        ts_xyyyz_xxxyy[i] = ts_xyyy_xxxyy[i] * pa_z[i];

        ts_xyyyz_xxxyz[i] = 3.0 * ts_yyyz_xxyz[i] * fe_0 + ts_yyyz_xxxyz[i] * pa_x[i];

        ts_xyyyz_xxxzz[i] = 3.0 * ts_yyyz_xxzz[i] * fe_0 + ts_yyyz_xxxzz[i] * pa_x[i];

        ts_xyyyz_xxyyy[i] = ts_xyyy_xxyyy[i] * pa_z[i];

        ts_xyyyz_xxyyz[i] = 2.0 * ts_yyyz_xyyz[i] * fe_0 + ts_yyyz_xxyyz[i] * pa_x[i];

        ts_xyyyz_xxyzz[i] = 2.0 * ts_yyyz_xyzz[i] * fe_0 + ts_yyyz_xxyzz[i] * pa_x[i];

        ts_xyyyz_xxzzz[i] = 2.0 * ts_yyyz_xzzz[i] * fe_0 + ts_yyyz_xxzzz[i] * pa_x[i];

        ts_xyyyz_xyyyy[i] = ts_xyyy_xyyyy[i] * pa_z[i];

        ts_xyyyz_xyyyz[i] = ts_yyyz_yyyz[i] * fe_0 + ts_yyyz_xyyyz[i] * pa_x[i];

        ts_xyyyz_xyyzz[i] = ts_yyyz_yyzz[i] * fe_0 + ts_yyyz_xyyzz[i] * pa_x[i];

        ts_xyyyz_xyzzz[i] = ts_yyyz_yzzz[i] * fe_0 + ts_yyyz_xyzzz[i] * pa_x[i];

        ts_xyyyz_xzzzz[i] = ts_yyyz_zzzz[i] * fe_0 + ts_yyyz_xzzzz[i] * pa_x[i];

        ts_xyyyz_yyyyy[i] = ts_yyyz_yyyyy[i] * pa_x[i];

        ts_xyyyz_yyyyz[i] = ts_yyyz_yyyyz[i] * pa_x[i];

        ts_xyyyz_yyyzz[i] = ts_yyyz_yyyzz[i] * pa_x[i];

        ts_xyyyz_yyzzz[i] = ts_yyyz_yyzzz[i] * pa_x[i];

        ts_xyyyz_yzzzz[i] = ts_yyyz_yzzzz[i] * pa_x[i];

        ts_xyyyz_zzzzz[i] = ts_yyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 252-273 components of targeted buffer : HH

    auto ts_xyyzz_xxxxx = pbuffer.data(idx_ovl_hh + 252);

    auto ts_xyyzz_xxxxy = pbuffer.data(idx_ovl_hh + 253);

    auto ts_xyyzz_xxxxz = pbuffer.data(idx_ovl_hh + 254);

    auto ts_xyyzz_xxxyy = pbuffer.data(idx_ovl_hh + 255);

    auto ts_xyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 256);

    auto ts_xyyzz_xxxzz = pbuffer.data(idx_ovl_hh + 257);

    auto ts_xyyzz_xxyyy = pbuffer.data(idx_ovl_hh + 258);

    auto ts_xyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 259);

    auto ts_xyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 260);

    auto ts_xyyzz_xxzzz = pbuffer.data(idx_ovl_hh + 261);

    auto ts_xyyzz_xyyyy = pbuffer.data(idx_ovl_hh + 262);

    auto ts_xyyzz_xyyyz = pbuffer.data(idx_ovl_hh + 263);

    auto ts_xyyzz_xyyzz = pbuffer.data(idx_ovl_hh + 264);

    auto ts_xyyzz_xyzzz = pbuffer.data(idx_ovl_hh + 265);

    auto ts_xyyzz_xzzzz = pbuffer.data(idx_ovl_hh + 266);

    auto ts_xyyzz_yyyyy = pbuffer.data(idx_ovl_hh + 267);

    auto ts_xyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 268);

    auto ts_xyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 269);

    auto ts_xyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 270);

    auto ts_xyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 271);

    auto ts_xyyzz_zzzzz = pbuffer.data(idx_ovl_hh + 272);

    #pragma omp simd aligned(pa_x, ts_xyyzz_xxxxx, ts_xyyzz_xxxxy, ts_xyyzz_xxxxz, ts_xyyzz_xxxyy, ts_xyyzz_xxxyz, ts_xyyzz_xxxzz, ts_xyyzz_xxyyy, ts_xyyzz_xxyyz, ts_xyyzz_xxyzz, ts_xyyzz_xxzzz, ts_xyyzz_xyyyy, ts_xyyzz_xyyyz, ts_xyyzz_xyyzz, ts_xyyzz_xyzzz, ts_xyyzz_xzzzz, ts_xyyzz_yyyyy, ts_xyyzz_yyyyz, ts_xyyzz_yyyzz, ts_xyyzz_yyzzz, ts_xyyzz_yzzzz, ts_xyyzz_zzzzz, ts_yyzz_xxxx, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxy, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxz, ts_yyzz_xxxzz, ts_yyzz_xxyy, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyz, ts_yyzz_xxyzz, ts_yyzz_xxzz, ts_yyzz_xxzzz, ts_yyzz_xyyy, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyz, ts_yyzz_xyyzz, ts_yyzz_xyzz, ts_yyzz_xyzzz, ts_yyzz_xzzz, ts_yyzz_xzzzz, ts_yyzz_yyyy, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyz, ts_yyzz_yyyzz, ts_yyzz_yyzz, ts_yyzz_yyzzz, ts_yyzz_yzzz, ts_yyzz_yzzzz, ts_yyzz_zzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzz_xxxxx[i] = 5.0 * ts_yyzz_xxxx[i] * fe_0 + ts_yyzz_xxxxx[i] * pa_x[i];

        ts_xyyzz_xxxxy[i] = 4.0 * ts_yyzz_xxxy[i] * fe_0 + ts_yyzz_xxxxy[i] * pa_x[i];

        ts_xyyzz_xxxxz[i] = 4.0 * ts_yyzz_xxxz[i] * fe_0 + ts_yyzz_xxxxz[i] * pa_x[i];

        ts_xyyzz_xxxyy[i] = 3.0 * ts_yyzz_xxyy[i] * fe_0 + ts_yyzz_xxxyy[i] * pa_x[i];

        ts_xyyzz_xxxyz[i] = 3.0 * ts_yyzz_xxyz[i] * fe_0 + ts_yyzz_xxxyz[i] * pa_x[i];

        ts_xyyzz_xxxzz[i] = 3.0 * ts_yyzz_xxzz[i] * fe_0 + ts_yyzz_xxxzz[i] * pa_x[i];

        ts_xyyzz_xxyyy[i] = 2.0 * ts_yyzz_xyyy[i] * fe_0 + ts_yyzz_xxyyy[i] * pa_x[i];

        ts_xyyzz_xxyyz[i] = 2.0 * ts_yyzz_xyyz[i] * fe_0 + ts_yyzz_xxyyz[i] * pa_x[i];

        ts_xyyzz_xxyzz[i] = 2.0 * ts_yyzz_xyzz[i] * fe_0 + ts_yyzz_xxyzz[i] * pa_x[i];

        ts_xyyzz_xxzzz[i] = 2.0 * ts_yyzz_xzzz[i] * fe_0 + ts_yyzz_xxzzz[i] * pa_x[i];

        ts_xyyzz_xyyyy[i] = ts_yyzz_yyyy[i] * fe_0 + ts_yyzz_xyyyy[i] * pa_x[i];

        ts_xyyzz_xyyyz[i] = ts_yyzz_yyyz[i] * fe_0 + ts_yyzz_xyyyz[i] * pa_x[i];

        ts_xyyzz_xyyzz[i] = ts_yyzz_yyzz[i] * fe_0 + ts_yyzz_xyyzz[i] * pa_x[i];

        ts_xyyzz_xyzzz[i] = ts_yyzz_yzzz[i] * fe_0 + ts_yyzz_xyzzz[i] * pa_x[i];

        ts_xyyzz_xzzzz[i] = ts_yyzz_zzzz[i] * fe_0 + ts_yyzz_xzzzz[i] * pa_x[i];

        ts_xyyzz_yyyyy[i] = ts_yyzz_yyyyy[i] * pa_x[i];

        ts_xyyzz_yyyyz[i] = ts_yyzz_yyyyz[i] * pa_x[i];

        ts_xyyzz_yyyzz[i] = ts_yyzz_yyyzz[i] * pa_x[i];

        ts_xyyzz_yyzzz[i] = ts_yyzz_yyzzz[i] * pa_x[i];

        ts_xyyzz_yzzzz[i] = ts_yyzz_yzzzz[i] * pa_x[i];

        ts_xyyzz_zzzzz[i] = ts_yyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : HH

    auto ts_xyzzz_xxxxx = pbuffer.data(idx_ovl_hh + 273);

    auto ts_xyzzz_xxxxy = pbuffer.data(idx_ovl_hh + 274);

    auto ts_xyzzz_xxxxz = pbuffer.data(idx_ovl_hh + 275);

    auto ts_xyzzz_xxxyy = pbuffer.data(idx_ovl_hh + 276);

    auto ts_xyzzz_xxxyz = pbuffer.data(idx_ovl_hh + 277);

    auto ts_xyzzz_xxxzz = pbuffer.data(idx_ovl_hh + 278);

    auto ts_xyzzz_xxyyy = pbuffer.data(idx_ovl_hh + 279);

    auto ts_xyzzz_xxyyz = pbuffer.data(idx_ovl_hh + 280);

    auto ts_xyzzz_xxyzz = pbuffer.data(idx_ovl_hh + 281);

    auto ts_xyzzz_xxzzz = pbuffer.data(idx_ovl_hh + 282);

    auto ts_xyzzz_xyyyy = pbuffer.data(idx_ovl_hh + 283);

    auto ts_xyzzz_xyyyz = pbuffer.data(idx_ovl_hh + 284);

    auto ts_xyzzz_xyyzz = pbuffer.data(idx_ovl_hh + 285);

    auto ts_xyzzz_xyzzz = pbuffer.data(idx_ovl_hh + 286);

    auto ts_xyzzz_xzzzz = pbuffer.data(idx_ovl_hh + 287);

    auto ts_xyzzz_yyyyy = pbuffer.data(idx_ovl_hh + 288);

    auto ts_xyzzz_yyyyz = pbuffer.data(idx_ovl_hh + 289);

    auto ts_xyzzz_yyyzz = pbuffer.data(idx_ovl_hh + 290);

    auto ts_xyzzz_yyzzz = pbuffer.data(idx_ovl_hh + 291);

    auto ts_xyzzz_yzzzz = pbuffer.data(idx_ovl_hh + 292);

    auto ts_xyzzz_zzzzz = pbuffer.data(idx_ovl_hh + 293);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzz_xxxxx, ts_xyzzz_xxxxy, ts_xyzzz_xxxxz, ts_xyzzz_xxxyy, ts_xyzzz_xxxyz, ts_xyzzz_xxxzz, ts_xyzzz_xxyyy, ts_xyzzz_xxyyz, ts_xyzzz_xxyzz, ts_xyzzz_xxzzz, ts_xyzzz_xyyyy, ts_xyzzz_xyyyz, ts_xyzzz_xyyzz, ts_xyzzz_xyzzz, ts_xyzzz_xzzzz, ts_xyzzz_yyyyy, ts_xyzzz_yyyyz, ts_xyzzz_yyyzz, ts_xyzzz_yyzzz, ts_xyzzz_yzzzz, ts_xyzzz_zzzzz, ts_xzzz_xxxxx, ts_xzzz_xxxxz, ts_xzzz_xxxzz, ts_xzzz_xxzzz, ts_xzzz_xzzzz, ts_yzzz_xxxxy, ts_yzzz_xxxy, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxyy, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyz, ts_yzzz_xxyzz, ts_yzzz_xyyy, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyz, ts_yzzz_xyyzz, ts_yzzz_xyzz, ts_yzzz_xyzzz, ts_yzzz_yyyy, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyz, ts_yzzz_yyyzz, ts_yzzz_yyzz, ts_yzzz_yyzzz, ts_yzzz_yzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzz_xxxxx[i] = ts_xzzz_xxxxx[i] * pa_y[i];

        ts_xyzzz_xxxxy[i] = 4.0 * ts_yzzz_xxxy[i] * fe_0 + ts_yzzz_xxxxy[i] * pa_x[i];

        ts_xyzzz_xxxxz[i] = ts_xzzz_xxxxz[i] * pa_y[i];

        ts_xyzzz_xxxyy[i] = 3.0 * ts_yzzz_xxyy[i] * fe_0 + ts_yzzz_xxxyy[i] * pa_x[i];

        ts_xyzzz_xxxyz[i] = 3.0 * ts_yzzz_xxyz[i] * fe_0 + ts_yzzz_xxxyz[i] * pa_x[i];

        ts_xyzzz_xxxzz[i] = ts_xzzz_xxxzz[i] * pa_y[i];

        ts_xyzzz_xxyyy[i] = 2.0 * ts_yzzz_xyyy[i] * fe_0 + ts_yzzz_xxyyy[i] * pa_x[i];

        ts_xyzzz_xxyyz[i] = 2.0 * ts_yzzz_xyyz[i] * fe_0 + ts_yzzz_xxyyz[i] * pa_x[i];

        ts_xyzzz_xxyzz[i] = 2.0 * ts_yzzz_xyzz[i] * fe_0 + ts_yzzz_xxyzz[i] * pa_x[i];

        ts_xyzzz_xxzzz[i] = ts_xzzz_xxzzz[i] * pa_y[i];

        ts_xyzzz_xyyyy[i] = ts_yzzz_yyyy[i] * fe_0 + ts_yzzz_xyyyy[i] * pa_x[i];

        ts_xyzzz_xyyyz[i] = ts_yzzz_yyyz[i] * fe_0 + ts_yzzz_xyyyz[i] * pa_x[i];

        ts_xyzzz_xyyzz[i] = ts_yzzz_yyzz[i] * fe_0 + ts_yzzz_xyyzz[i] * pa_x[i];

        ts_xyzzz_xyzzz[i] = ts_yzzz_yzzz[i] * fe_0 + ts_yzzz_xyzzz[i] * pa_x[i];

        ts_xyzzz_xzzzz[i] = ts_xzzz_xzzzz[i] * pa_y[i];

        ts_xyzzz_yyyyy[i] = ts_yzzz_yyyyy[i] * pa_x[i];

        ts_xyzzz_yyyyz[i] = ts_yzzz_yyyyz[i] * pa_x[i];

        ts_xyzzz_yyyzz[i] = ts_yzzz_yyyzz[i] * pa_x[i];

        ts_xyzzz_yyzzz[i] = ts_yzzz_yyzzz[i] * pa_x[i];

        ts_xyzzz_yzzzz[i] = ts_yzzz_yzzzz[i] * pa_x[i];

        ts_xyzzz_zzzzz[i] = ts_yzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 294-315 components of targeted buffer : HH

    auto ts_xzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 294);

    auto ts_xzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 295);

    auto ts_xzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 296);

    auto ts_xzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 297);

    auto ts_xzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 298);

    auto ts_xzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 299);

    auto ts_xzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 300);

    auto ts_xzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 301);

    auto ts_xzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 302);

    auto ts_xzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 303);

    auto ts_xzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 304);

    auto ts_xzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 305);

    auto ts_xzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 306);

    auto ts_xzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 307);

    auto ts_xzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 308);

    auto ts_xzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 309);

    auto ts_xzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 310);

    auto ts_xzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 311);

    auto ts_xzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 312);

    auto ts_xzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 313);

    auto ts_xzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 314);

    #pragma omp simd aligned(pa_x, ts_xzzzz_xxxxx, ts_xzzzz_xxxxy, ts_xzzzz_xxxxz, ts_xzzzz_xxxyy, ts_xzzzz_xxxyz, ts_xzzzz_xxxzz, ts_xzzzz_xxyyy, ts_xzzzz_xxyyz, ts_xzzzz_xxyzz, ts_xzzzz_xxzzz, ts_xzzzz_xyyyy, ts_xzzzz_xyyyz, ts_xzzzz_xyyzz, ts_xzzzz_xyzzz, ts_xzzzz_xzzzz, ts_xzzzz_yyyyy, ts_xzzzz_yyyyz, ts_xzzzz_yyyzz, ts_xzzzz_yyzzz, ts_xzzzz_yzzzz, ts_xzzzz_zzzzz, ts_zzzz_xxxx, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxy, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxz, ts_zzzz_xxxzz, ts_zzzz_xxyy, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyz, ts_zzzz_xxyzz, ts_zzzz_xxzz, ts_zzzz_xxzzz, ts_zzzz_xyyy, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyz, ts_zzzz_xyyzz, ts_zzzz_xyzz, ts_zzzz_xyzzz, ts_zzzz_xzzz, ts_zzzz_xzzzz, ts_zzzz_yyyy, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyz, ts_zzzz_yyyzz, ts_zzzz_yyzz, ts_zzzz_yyzzz, ts_zzzz_yzzz, ts_zzzz_yzzzz, ts_zzzz_zzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzz_xxxxx[i] = 5.0 * ts_zzzz_xxxx[i] * fe_0 + ts_zzzz_xxxxx[i] * pa_x[i];

        ts_xzzzz_xxxxy[i] = 4.0 * ts_zzzz_xxxy[i] * fe_0 + ts_zzzz_xxxxy[i] * pa_x[i];

        ts_xzzzz_xxxxz[i] = 4.0 * ts_zzzz_xxxz[i] * fe_0 + ts_zzzz_xxxxz[i] * pa_x[i];

        ts_xzzzz_xxxyy[i] = 3.0 * ts_zzzz_xxyy[i] * fe_0 + ts_zzzz_xxxyy[i] * pa_x[i];

        ts_xzzzz_xxxyz[i] = 3.0 * ts_zzzz_xxyz[i] * fe_0 + ts_zzzz_xxxyz[i] * pa_x[i];

        ts_xzzzz_xxxzz[i] = 3.0 * ts_zzzz_xxzz[i] * fe_0 + ts_zzzz_xxxzz[i] * pa_x[i];

        ts_xzzzz_xxyyy[i] = 2.0 * ts_zzzz_xyyy[i] * fe_0 + ts_zzzz_xxyyy[i] * pa_x[i];

        ts_xzzzz_xxyyz[i] = 2.0 * ts_zzzz_xyyz[i] * fe_0 + ts_zzzz_xxyyz[i] * pa_x[i];

        ts_xzzzz_xxyzz[i] = 2.0 * ts_zzzz_xyzz[i] * fe_0 + ts_zzzz_xxyzz[i] * pa_x[i];

        ts_xzzzz_xxzzz[i] = 2.0 * ts_zzzz_xzzz[i] * fe_0 + ts_zzzz_xxzzz[i] * pa_x[i];

        ts_xzzzz_xyyyy[i] = ts_zzzz_yyyy[i] * fe_0 + ts_zzzz_xyyyy[i] * pa_x[i];

        ts_xzzzz_xyyyz[i] = ts_zzzz_yyyz[i] * fe_0 + ts_zzzz_xyyyz[i] * pa_x[i];

        ts_xzzzz_xyyzz[i] = ts_zzzz_yyzz[i] * fe_0 + ts_zzzz_xyyzz[i] * pa_x[i];

        ts_xzzzz_xyzzz[i] = ts_zzzz_yzzz[i] * fe_0 + ts_zzzz_xyzzz[i] * pa_x[i];

        ts_xzzzz_xzzzz[i] = ts_zzzz_zzzz[i] * fe_0 + ts_zzzz_xzzzz[i] * pa_x[i];

        ts_xzzzz_yyyyy[i] = ts_zzzz_yyyyy[i] * pa_x[i];

        ts_xzzzz_yyyyz[i] = ts_zzzz_yyyyz[i] * pa_x[i];

        ts_xzzzz_yyyzz[i] = ts_zzzz_yyyzz[i] * pa_x[i];

        ts_xzzzz_yyzzz[i] = ts_zzzz_yyzzz[i] * pa_x[i];

        ts_xzzzz_yzzzz[i] = ts_zzzz_yzzzz[i] * pa_x[i];

        ts_xzzzz_zzzzz[i] = ts_zzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : HH

    auto ts_yyyyy_xxxxx = pbuffer.data(idx_ovl_hh + 315);

    auto ts_yyyyy_xxxxy = pbuffer.data(idx_ovl_hh + 316);

    auto ts_yyyyy_xxxxz = pbuffer.data(idx_ovl_hh + 317);

    auto ts_yyyyy_xxxyy = pbuffer.data(idx_ovl_hh + 318);

    auto ts_yyyyy_xxxyz = pbuffer.data(idx_ovl_hh + 319);

    auto ts_yyyyy_xxxzz = pbuffer.data(idx_ovl_hh + 320);

    auto ts_yyyyy_xxyyy = pbuffer.data(idx_ovl_hh + 321);

    auto ts_yyyyy_xxyyz = pbuffer.data(idx_ovl_hh + 322);

    auto ts_yyyyy_xxyzz = pbuffer.data(idx_ovl_hh + 323);

    auto ts_yyyyy_xxzzz = pbuffer.data(idx_ovl_hh + 324);

    auto ts_yyyyy_xyyyy = pbuffer.data(idx_ovl_hh + 325);

    auto ts_yyyyy_xyyyz = pbuffer.data(idx_ovl_hh + 326);

    auto ts_yyyyy_xyyzz = pbuffer.data(idx_ovl_hh + 327);

    auto ts_yyyyy_xyzzz = pbuffer.data(idx_ovl_hh + 328);

    auto ts_yyyyy_xzzzz = pbuffer.data(idx_ovl_hh + 329);

    auto ts_yyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 330);

    auto ts_yyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 331);

    auto ts_yyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 332);

    auto ts_yyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 333);

    auto ts_yyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 334);

    auto ts_yyyyy_zzzzz = pbuffer.data(idx_ovl_hh + 335);

    #pragma omp simd aligned(pa_y, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxzz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xxzzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_xzzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, ts_yyyy_xxxx, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxz, ts_yyyy_xxxzz, ts_yyyy_xxyy, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyz, ts_yyyy_xxyzz, ts_yyyy_xxzz, ts_yyyy_xxzzz, ts_yyyy_xyyy, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyz, ts_yyyy_xyyzz, ts_yyyy_xyzz, ts_yyyy_xyzzz, ts_yyyy_xzzz, ts_yyyy_xzzzz, ts_yyyy_yyyy, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyz, ts_yyyy_yyyzz, ts_yyyy_yyzz, ts_yyyy_yyzzz, ts_yyyy_yzzz, ts_yyyy_yzzzz, ts_yyyy_zzzz, ts_yyyy_zzzzz, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxzz, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyy_xxxxx[i] = 4.0 * ts_yyy_xxxxx[i] * fe_0 + ts_yyyy_xxxxx[i] * pa_y[i];

        ts_yyyyy_xxxxy[i] = 4.0 * ts_yyy_xxxxy[i] * fe_0 + ts_yyyy_xxxx[i] * fe_0 + ts_yyyy_xxxxy[i] * pa_y[i];

        ts_yyyyy_xxxxz[i] = 4.0 * ts_yyy_xxxxz[i] * fe_0 + ts_yyyy_xxxxz[i] * pa_y[i];

        ts_yyyyy_xxxyy[i] = 4.0 * ts_yyy_xxxyy[i] * fe_0 + 2.0 * ts_yyyy_xxxy[i] * fe_0 + ts_yyyy_xxxyy[i] * pa_y[i];

        ts_yyyyy_xxxyz[i] = 4.0 * ts_yyy_xxxyz[i] * fe_0 + ts_yyyy_xxxz[i] * fe_0 + ts_yyyy_xxxyz[i] * pa_y[i];

        ts_yyyyy_xxxzz[i] = 4.0 * ts_yyy_xxxzz[i] * fe_0 + ts_yyyy_xxxzz[i] * pa_y[i];

        ts_yyyyy_xxyyy[i] = 4.0 * ts_yyy_xxyyy[i] * fe_0 + 3.0 * ts_yyyy_xxyy[i] * fe_0 + ts_yyyy_xxyyy[i] * pa_y[i];

        ts_yyyyy_xxyyz[i] = 4.0 * ts_yyy_xxyyz[i] * fe_0 + 2.0 * ts_yyyy_xxyz[i] * fe_0 + ts_yyyy_xxyyz[i] * pa_y[i];

        ts_yyyyy_xxyzz[i] = 4.0 * ts_yyy_xxyzz[i] * fe_0 + ts_yyyy_xxzz[i] * fe_0 + ts_yyyy_xxyzz[i] * pa_y[i];

        ts_yyyyy_xxzzz[i] = 4.0 * ts_yyy_xxzzz[i] * fe_0 + ts_yyyy_xxzzz[i] * pa_y[i];

        ts_yyyyy_xyyyy[i] = 4.0 * ts_yyy_xyyyy[i] * fe_0 + 4.0 * ts_yyyy_xyyy[i] * fe_0 + ts_yyyy_xyyyy[i] * pa_y[i];

        ts_yyyyy_xyyyz[i] = 4.0 * ts_yyy_xyyyz[i] * fe_0 + 3.0 * ts_yyyy_xyyz[i] * fe_0 + ts_yyyy_xyyyz[i] * pa_y[i];

        ts_yyyyy_xyyzz[i] = 4.0 * ts_yyy_xyyzz[i] * fe_0 + 2.0 * ts_yyyy_xyzz[i] * fe_0 + ts_yyyy_xyyzz[i] * pa_y[i];

        ts_yyyyy_xyzzz[i] = 4.0 * ts_yyy_xyzzz[i] * fe_0 + ts_yyyy_xzzz[i] * fe_0 + ts_yyyy_xyzzz[i] * pa_y[i];

        ts_yyyyy_xzzzz[i] = 4.0 * ts_yyy_xzzzz[i] * fe_0 + ts_yyyy_xzzzz[i] * pa_y[i];

        ts_yyyyy_yyyyy[i] = 4.0 * ts_yyy_yyyyy[i] * fe_0 + 5.0 * ts_yyyy_yyyy[i] * fe_0 + ts_yyyy_yyyyy[i] * pa_y[i];

        ts_yyyyy_yyyyz[i] = 4.0 * ts_yyy_yyyyz[i] * fe_0 + 4.0 * ts_yyyy_yyyz[i] * fe_0 + ts_yyyy_yyyyz[i] * pa_y[i];

        ts_yyyyy_yyyzz[i] = 4.0 * ts_yyy_yyyzz[i] * fe_0 + 3.0 * ts_yyyy_yyzz[i] * fe_0 + ts_yyyy_yyyzz[i] * pa_y[i];

        ts_yyyyy_yyzzz[i] = 4.0 * ts_yyy_yyzzz[i] * fe_0 + 2.0 * ts_yyyy_yzzz[i] * fe_0 + ts_yyyy_yyzzz[i] * pa_y[i];

        ts_yyyyy_yzzzz[i] = 4.0 * ts_yyy_yzzzz[i] * fe_0 + ts_yyyy_zzzz[i] * fe_0 + ts_yyyy_yzzzz[i] * pa_y[i];

        ts_yyyyy_zzzzz[i] = 4.0 * ts_yyy_zzzzz[i] * fe_0 + ts_yyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 336-357 components of targeted buffer : HH

    auto ts_yyyyz_xxxxx = pbuffer.data(idx_ovl_hh + 336);

    auto ts_yyyyz_xxxxy = pbuffer.data(idx_ovl_hh + 337);

    auto ts_yyyyz_xxxxz = pbuffer.data(idx_ovl_hh + 338);

    auto ts_yyyyz_xxxyy = pbuffer.data(idx_ovl_hh + 339);

    auto ts_yyyyz_xxxyz = pbuffer.data(idx_ovl_hh + 340);

    auto ts_yyyyz_xxxzz = pbuffer.data(idx_ovl_hh + 341);

    auto ts_yyyyz_xxyyy = pbuffer.data(idx_ovl_hh + 342);

    auto ts_yyyyz_xxyyz = pbuffer.data(idx_ovl_hh + 343);

    auto ts_yyyyz_xxyzz = pbuffer.data(idx_ovl_hh + 344);

    auto ts_yyyyz_xxzzz = pbuffer.data(idx_ovl_hh + 345);

    auto ts_yyyyz_xyyyy = pbuffer.data(idx_ovl_hh + 346);

    auto ts_yyyyz_xyyyz = pbuffer.data(idx_ovl_hh + 347);

    auto ts_yyyyz_xyyzz = pbuffer.data(idx_ovl_hh + 348);

    auto ts_yyyyz_xyzzz = pbuffer.data(idx_ovl_hh + 349);

    auto ts_yyyyz_xzzzz = pbuffer.data(idx_ovl_hh + 350);

    auto ts_yyyyz_yyyyy = pbuffer.data(idx_ovl_hh + 351);

    auto ts_yyyyz_yyyyz = pbuffer.data(idx_ovl_hh + 352);

    auto ts_yyyyz_yyyzz = pbuffer.data(idx_ovl_hh + 353);

    auto ts_yyyyz_yyzzz = pbuffer.data(idx_ovl_hh + 354);

    auto ts_yyyyz_yzzzz = pbuffer.data(idx_ovl_hh + 355);

    auto ts_yyyyz_zzzzz = pbuffer.data(idx_ovl_hh + 356);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxyy, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyz, ts_yyyy_xxyzz, ts_yyyy_xyyy, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyz, ts_yyyy_xyyzz, ts_yyyy_xyzz, ts_yyyy_xyzzz, ts_yyyy_yyyy, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyz, ts_yyyy_yyyzz, ts_yyyy_yyzz, ts_yyyy_yyzzz, ts_yyyy_yzzz, ts_yyyy_yzzzz, ts_yyyyz_xxxxx, ts_yyyyz_xxxxy, ts_yyyyz_xxxxz, ts_yyyyz_xxxyy, ts_yyyyz_xxxyz, ts_yyyyz_xxxzz, ts_yyyyz_xxyyy, ts_yyyyz_xxyyz, ts_yyyyz_xxyzz, ts_yyyyz_xxzzz, ts_yyyyz_xyyyy, ts_yyyyz_xyyyz, ts_yyyyz_xyyzz, ts_yyyyz_xyzzz, ts_yyyyz_xzzzz, ts_yyyyz_yyyyy, ts_yyyyz_yyyyz, ts_yyyyz_yyyzz, ts_yyyyz_yyzzz, ts_yyyyz_yzzzz, ts_yyyyz_zzzzz, ts_yyyz_xxxxz, ts_yyyz_xxxzz, ts_yyyz_xxzzz, ts_yyyz_xzzzz, ts_yyyz_zzzzz, ts_yyz_xxxxz, ts_yyz_xxxzz, ts_yyz_xxzzz, ts_yyz_xzzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyz_xxxxx[i] = ts_yyyy_xxxxx[i] * pa_z[i];

        ts_yyyyz_xxxxy[i] = ts_yyyy_xxxxy[i] * pa_z[i];

        ts_yyyyz_xxxxz[i] = 3.0 * ts_yyz_xxxxz[i] * fe_0 + ts_yyyz_xxxxz[i] * pa_y[i];

        ts_yyyyz_xxxyy[i] = ts_yyyy_xxxyy[i] * pa_z[i];

        ts_yyyyz_xxxyz[i] = ts_yyyy_xxxy[i] * fe_0 + ts_yyyy_xxxyz[i] * pa_z[i];

        ts_yyyyz_xxxzz[i] = 3.0 * ts_yyz_xxxzz[i] * fe_0 + ts_yyyz_xxxzz[i] * pa_y[i];

        ts_yyyyz_xxyyy[i] = ts_yyyy_xxyyy[i] * pa_z[i];

        ts_yyyyz_xxyyz[i] = ts_yyyy_xxyy[i] * fe_0 + ts_yyyy_xxyyz[i] * pa_z[i];

        ts_yyyyz_xxyzz[i] = 2.0 * ts_yyyy_xxyz[i] * fe_0 + ts_yyyy_xxyzz[i] * pa_z[i];

        ts_yyyyz_xxzzz[i] = 3.0 * ts_yyz_xxzzz[i] * fe_0 + ts_yyyz_xxzzz[i] * pa_y[i];

        ts_yyyyz_xyyyy[i] = ts_yyyy_xyyyy[i] * pa_z[i];

        ts_yyyyz_xyyyz[i] = ts_yyyy_xyyy[i] * fe_0 + ts_yyyy_xyyyz[i] * pa_z[i];

        ts_yyyyz_xyyzz[i] = 2.0 * ts_yyyy_xyyz[i] * fe_0 + ts_yyyy_xyyzz[i] * pa_z[i];

        ts_yyyyz_xyzzz[i] = 3.0 * ts_yyyy_xyzz[i] * fe_0 + ts_yyyy_xyzzz[i] * pa_z[i];

        ts_yyyyz_xzzzz[i] = 3.0 * ts_yyz_xzzzz[i] * fe_0 + ts_yyyz_xzzzz[i] * pa_y[i];

        ts_yyyyz_yyyyy[i] = ts_yyyy_yyyyy[i] * pa_z[i];

        ts_yyyyz_yyyyz[i] = ts_yyyy_yyyy[i] * fe_0 + ts_yyyy_yyyyz[i] * pa_z[i];

        ts_yyyyz_yyyzz[i] = 2.0 * ts_yyyy_yyyz[i] * fe_0 + ts_yyyy_yyyzz[i] * pa_z[i];

        ts_yyyyz_yyzzz[i] = 3.0 * ts_yyyy_yyzz[i] * fe_0 + ts_yyyy_yyzzz[i] * pa_z[i];

        ts_yyyyz_yzzzz[i] = 4.0 * ts_yyyy_yzzz[i] * fe_0 + ts_yyyy_yzzzz[i] * pa_z[i];

        ts_yyyyz_zzzzz[i] = 3.0 * ts_yyz_zzzzz[i] * fe_0 + ts_yyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 357-378 components of targeted buffer : HH

    auto ts_yyyzz_xxxxx = pbuffer.data(idx_ovl_hh + 357);

    auto ts_yyyzz_xxxxy = pbuffer.data(idx_ovl_hh + 358);

    auto ts_yyyzz_xxxxz = pbuffer.data(idx_ovl_hh + 359);

    auto ts_yyyzz_xxxyy = pbuffer.data(idx_ovl_hh + 360);

    auto ts_yyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 361);

    auto ts_yyyzz_xxxzz = pbuffer.data(idx_ovl_hh + 362);

    auto ts_yyyzz_xxyyy = pbuffer.data(idx_ovl_hh + 363);

    auto ts_yyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 364);

    auto ts_yyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 365);

    auto ts_yyyzz_xxzzz = pbuffer.data(idx_ovl_hh + 366);

    auto ts_yyyzz_xyyyy = pbuffer.data(idx_ovl_hh + 367);

    auto ts_yyyzz_xyyyz = pbuffer.data(idx_ovl_hh + 368);

    auto ts_yyyzz_xyyzz = pbuffer.data(idx_ovl_hh + 369);

    auto ts_yyyzz_xyzzz = pbuffer.data(idx_ovl_hh + 370);

    auto ts_yyyzz_xzzzz = pbuffer.data(idx_ovl_hh + 371);

    auto ts_yyyzz_yyyyy = pbuffer.data(idx_ovl_hh + 372);

    auto ts_yyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 373);

    auto ts_yyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 374);

    auto ts_yyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 375);

    auto ts_yyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 376);

    auto ts_yyyzz_zzzzz = pbuffer.data(idx_ovl_hh + 377);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyy_xxxxy, ts_yyy_xxxyy, ts_yyy_xxyyy, ts_yyy_xyyyy, ts_yyy_yyyyy, ts_yyyz_xxxxy, ts_yyyz_xxxyy, ts_yyyz_xxyyy, ts_yyyz_xyyyy, ts_yyyz_yyyyy, ts_yyyzz_xxxxx, ts_yyyzz_xxxxy, ts_yyyzz_xxxxz, ts_yyyzz_xxxyy, ts_yyyzz_xxxyz, ts_yyyzz_xxxzz, ts_yyyzz_xxyyy, ts_yyyzz_xxyyz, ts_yyyzz_xxyzz, ts_yyyzz_xxzzz, ts_yyyzz_xyyyy, ts_yyyzz_xyyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzzz, ts_yyyzz_xzzzz, ts_yyyzz_yyyyy, ts_yyyzz_yyyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxz, ts_yyzz_xxxyz, ts_yyzz_xxxz, ts_yyzz_xxxzz, ts_yyzz_xxyyz, ts_yyzz_xxyz, ts_yyzz_xxyzz, ts_yyzz_xxzz, ts_yyzz_xxzzz, ts_yyzz_xyyyz, ts_yyzz_xyyz, ts_yyzz_xyyzz, ts_yyzz_xyzz, ts_yyzz_xyzzz, ts_yyzz_xzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyz, ts_yyzz_yyyz, ts_yyzz_yyyzz, ts_yyzz_yyzz, ts_yyzz_yyzzz, ts_yyzz_yzzz, ts_yyzz_yzzzz, ts_yyzz_zzzz, ts_yyzz_zzzzz, ts_yzz_xxxxx, ts_yzz_xxxxz, ts_yzz_xxxyz, ts_yzz_xxxzz, ts_yzz_xxyyz, ts_yzz_xxyzz, ts_yzz_xxzzz, ts_yzz_xyyyz, ts_yzz_xyyzz, ts_yzz_xyzzz, ts_yzz_xzzzz, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzz_xxxxx[i] = 2.0 * ts_yzz_xxxxx[i] * fe_0 + ts_yyzz_xxxxx[i] * pa_y[i];

        ts_yyyzz_xxxxy[i] = ts_yyy_xxxxy[i] * fe_0 + ts_yyyz_xxxxy[i] * pa_z[i];

        ts_yyyzz_xxxxz[i] = 2.0 * ts_yzz_xxxxz[i] * fe_0 + ts_yyzz_xxxxz[i] * pa_y[i];

        ts_yyyzz_xxxyy[i] = ts_yyy_xxxyy[i] * fe_0 + ts_yyyz_xxxyy[i] * pa_z[i];

        ts_yyyzz_xxxyz[i] = 2.0 * ts_yzz_xxxyz[i] * fe_0 + ts_yyzz_xxxz[i] * fe_0 + ts_yyzz_xxxyz[i] * pa_y[i];

        ts_yyyzz_xxxzz[i] = 2.0 * ts_yzz_xxxzz[i] * fe_0 + ts_yyzz_xxxzz[i] * pa_y[i];

        ts_yyyzz_xxyyy[i] = ts_yyy_xxyyy[i] * fe_0 + ts_yyyz_xxyyy[i] * pa_z[i];

        ts_yyyzz_xxyyz[i] = 2.0 * ts_yzz_xxyyz[i] * fe_0 + 2.0 * ts_yyzz_xxyz[i] * fe_0 + ts_yyzz_xxyyz[i] * pa_y[i];

        ts_yyyzz_xxyzz[i] = 2.0 * ts_yzz_xxyzz[i] * fe_0 + ts_yyzz_xxzz[i] * fe_0 + ts_yyzz_xxyzz[i] * pa_y[i];

        ts_yyyzz_xxzzz[i] = 2.0 * ts_yzz_xxzzz[i] * fe_0 + ts_yyzz_xxzzz[i] * pa_y[i];

        ts_yyyzz_xyyyy[i] = ts_yyy_xyyyy[i] * fe_0 + ts_yyyz_xyyyy[i] * pa_z[i];

        ts_yyyzz_xyyyz[i] = 2.0 * ts_yzz_xyyyz[i] * fe_0 + 3.0 * ts_yyzz_xyyz[i] * fe_0 + ts_yyzz_xyyyz[i] * pa_y[i];

        ts_yyyzz_xyyzz[i] = 2.0 * ts_yzz_xyyzz[i] * fe_0 + 2.0 * ts_yyzz_xyzz[i] * fe_0 + ts_yyzz_xyyzz[i] * pa_y[i];

        ts_yyyzz_xyzzz[i] = 2.0 * ts_yzz_xyzzz[i] * fe_0 + ts_yyzz_xzzz[i] * fe_0 + ts_yyzz_xyzzz[i] * pa_y[i];

        ts_yyyzz_xzzzz[i] = 2.0 * ts_yzz_xzzzz[i] * fe_0 + ts_yyzz_xzzzz[i] * pa_y[i];

        ts_yyyzz_yyyyy[i] = ts_yyy_yyyyy[i] * fe_0 + ts_yyyz_yyyyy[i] * pa_z[i];

        ts_yyyzz_yyyyz[i] = 2.0 * ts_yzz_yyyyz[i] * fe_0 + 4.0 * ts_yyzz_yyyz[i] * fe_0 + ts_yyzz_yyyyz[i] * pa_y[i];

        ts_yyyzz_yyyzz[i] = 2.0 * ts_yzz_yyyzz[i] * fe_0 + 3.0 * ts_yyzz_yyzz[i] * fe_0 + ts_yyzz_yyyzz[i] * pa_y[i];

        ts_yyyzz_yyzzz[i] = 2.0 * ts_yzz_yyzzz[i] * fe_0 + 2.0 * ts_yyzz_yzzz[i] * fe_0 + ts_yyzz_yyzzz[i] * pa_y[i];

        ts_yyyzz_yzzzz[i] = 2.0 * ts_yzz_yzzzz[i] * fe_0 + ts_yyzz_zzzz[i] * fe_0 + ts_yyzz_yzzzz[i] * pa_y[i];

        ts_yyyzz_zzzzz[i] = 2.0 * ts_yzz_zzzzz[i] * fe_0 + ts_yyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 378-399 components of targeted buffer : HH

    auto ts_yyzzz_xxxxx = pbuffer.data(idx_ovl_hh + 378);

    auto ts_yyzzz_xxxxy = pbuffer.data(idx_ovl_hh + 379);

    auto ts_yyzzz_xxxxz = pbuffer.data(idx_ovl_hh + 380);

    auto ts_yyzzz_xxxyy = pbuffer.data(idx_ovl_hh + 381);

    auto ts_yyzzz_xxxyz = pbuffer.data(idx_ovl_hh + 382);

    auto ts_yyzzz_xxxzz = pbuffer.data(idx_ovl_hh + 383);

    auto ts_yyzzz_xxyyy = pbuffer.data(idx_ovl_hh + 384);

    auto ts_yyzzz_xxyyz = pbuffer.data(idx_ovl_hh + 385);

    auto ts_yyzzz_xxyzz = pbuffer.data(idx_ovl_hh + 386);

    auto ts_yyzzz_xxzzz = pbuffer.data(idx_ovl_hh + 387);

    auto ts_yyzzz_xyyyy = pbuffer.data(idx_ovl_hh + 388);

    auto ts_yyzzz_xyyyz = pbuffer.data(idx_ovl_hh + 389);

    auto ts_yyzzz_xyyzz = pbuffer.data(idx_ovl_hh + 390);

    auto ts_yyzzz_xyzzz = pbuffer.data(idx_ovl_hh + 391);

    auto ts_yyzzz_xzzzz = pbuffer.data(idx_ovl_hh + 392);

    auto ts_yyzzz_yyyyy = pbuffer.data(idx_ovl_hh + 393);

    auto ts_yyzzz_yyyyz = pbuffer.data(idx_ovl_hh + 394);

    auto ts_yyzzz_yyyzz = pbuffer.data(idx_ovl_hh + 395);

    auto ts_yyzzz_yyzzz = pbuffer.data(idx_ovl_hh + 396);

    auto ts_yyzzz_yzzzz = pbuffer.data(idx_ovl_hh + 397);

    auto ts_yyzzz_zzzzz = pbuffer.data(idx_ovl_hh + 398);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyz_xxxxy, ts_yyz_xxxyy, ts_yyz_xxyyy, ts_yyz_xyyyy, ts_yyz_yyyyy, ts_yyzz_xxxxy, ts_yyzz_xxxyy, ts_yyzz_xxyyy, ts_yyzz_xyyyy, ts_yyzz_yyyyy, ts_yyzzz_xxxxx, ts_yyzzz_xxxxy, ts_yyzzz_xxxxz, ts_yyzzz_xxxyy, ts_yyzzz_xxxyz, ts_yyzzz_xxxzz, ts_yyzzz_xxyyy, ts_yyzzz_xxyyz, ts_yyzzz_xxyzz, ts_yyzzz_xxzzz, ts_yyzzz_xyyyy, ts_yyzzz_xyyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzzz, ts_yyzzz_xzzzz, ts_yyzzz_yyyyy, ts_yyzzz_yyyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzzz, ts_yzzz_xxxxx, ts_yzzz_xxxxz, ts_yzzz_xxxyz, ts_yzzz_xxxz, ts_yzzz_xxxzz, ts_yzzz_xxyyz, ts_yzzz_xxyz, ts_yzzz_xxyzz, ts_yzzz_xxzz, ts_yzzz_xxzzz, ts_yzzz_xyyyz, ts_yzzz_xyyz, ts_yzzz_xyyzz, ts_yzzz_xyzz, ts_yzzz_xyzzz, ts_yzzz_xzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyz, ts_yzzz_yyyz, ts_yzzz_yyyzz, ts_yzzz_yyzz, ts_yzzz_yyzzz, ts_yzzz_yzzz, ts_yzzz_yzzzz, ts_yzzz_zzzz, ts_yzzz_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxz, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzz_xxxxx[i] = ts_zzz_xxxxx[i] * fe_0 + ts_yzzz_xxxxx[i] * pa_y[i];

        ts_yyzzz_xxxxy[i] = 2.0 * ts_yyz_xxxxy[i] * fe_0 + ts_yyzz_xxxxy[i] * pa_z[i];

        ts_yyzzz_xxxxz[i] = ts_zzz_xxxxz[i] * fe_0 + ts_yzzz_xxxxz[i] * pa_y[i];

        ts_yyzzz_xxxyy[i] = 2.0 * ts_yyz_xxxyy[i] * fe_0 + ts_yyzz_xxxyy[i] * pa_z[i];

        ts_yyzzz_xxxyz[i] = ts_zzz_xxxyz[i] * fe_0 + ts_yzzz_xxxz[i] * fe_0 + ts_yzzz_xxxyz[i] * pa_y[i];

        ts_yyzzz_xxxzz[i] = ts_zzz_xxxzz[i] * fe_0 + ts_yzzz_xxxzz[i] * pa_y[i];

        ts_yyzzz_xxyyy[i] = 2.0 * ts_yyz_xxyyy[i] * fe_0 + ts_yyzz_xxyyy[i] * pa_z[i];

        ts_yyzzz_xxyyz[i] = ts_zzz_xxyyz[i] * fe_0 + 2.0 * ts_yzzz_xxyz[i] * fe_0 + ts_yzzz_xxyyz[i] * pa_y[i];

        ts_yyzzz_xxyzz[i] = ts_zzz_xxyzz[i] * fe_0 + ts_yzzz_xxzz[i] * fe_0 + ts_yzzz_xxyzz[i] * pa_y[i];

        ts_yyzzz_xxzzz[i] = ts_zzz_xxzzz[i] * fe_0 + ts_yzzz_xxzzz[i] * pa_y[i];

        ts_yyzzz_xyyyy[i] = 2.0 * ts_yyz_xyyyy[i] * fe_0 + ts_yyzz_xyyyy[i] * pa_z[i];

        ts_yyzzz_xyyyz[i] = ts_zzz_xyyyz[i] * fe_0 + 3.0 * ts_yzzz_xyyz[i] * fe_0 + ts_yzzz_xyyyz[i] * pa_y[i];

        ts_yyzzz_xyyzz[i] = ts_zzz_xyyzz[i] * fe_0 + 2.0 * ts_yzzz_xyzz[i] * fe_0 + ts_yzzz_xyyzz[i] * pa_y[i];

        ts_yyzzz_xyzzz[i] = ts_zzz_xyzzz[i] * fe_0 + ts_yzzz_xzzz[i] * fe_0 + ts_yzzz_xyzzz[i] * pa_y[i];

        ts_yyzzz_xzzzz[i] = ts_zzz_xzzzz[i] * fe_0 + ts_yzzz_xzzzz[i] * pa_y[i];

        ts_yyzzz_yyyyy[i] = 2.0 * ts_yyz_yyyyy[i] * fe_0 + ts_yyzz_yyyyy[i] * pa_z[i];

        ts_yyzzz_yyyyz[i] = ts_zzz_yyyyz[i] * fe_0 + 4.0 * ts_yzzz_yyyz[i] * fe_0 + ts_yzzz_yyyyz[i] * pa_y[i];

        ts_yyzzz_yyyzz[i] = ts_zzz_yyyzz[i] * fe_0 + 3.0 * ts_yzzz_yyzz[i] * fe_0 + ts_yzzz_yyyzz[i] * pa_y[i];

        ts_yyzzz_yyzzz[i] = ts_zzz_yyzzz[i] * fe_0 + 2.0 * ts_yzzz_yzzz[i] * fe_0 + ts_yzzz_yyzzz[i] * pa_y[i];

        ts_yyzzz_yzzzz[i] = ts_zzz_yzzzz[i] * fe_0 + ts_yzzz_zzzz[i] * fe_0 + ts_yzzz_yzzzz[i] * pa_y[i];

        ts_yyzzz_zzzzz[i] = ts_zzz_zzzzz[i] * fe_0 + ts_yzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 399-420 components of targeted buffer : HH

    auto ts_yzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 399);

    auto ts_yzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 400);

    auto ts_yzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 401);

    auto ts_yzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 402);

    auto ts_yzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 403);

    auto ts_yzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 404);

    auto ts_yzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 405);

    auto ts_yzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 406);

    auto ts_yzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 407);

    auto ts_yzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 408);

    auto ts_yzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 409);

    auto ts_yzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 410);

    auto ts_yzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 411);

    auto ts_yzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 412);

    auto ts_yzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 413);

    auto ts_yzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 414);

    auto ts_yzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 415);

    auto ts_yzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 416);

    auto ts_yzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 417);

    auto ts_yzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 418);

    auto ts_yzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 419);

    #pragma omp simd aligned(pa_y, ts_yzzzz_xxxxx, ts_yzzzz_xxxxy, ts_yzzzz_xxxxz, ts_yzzzz_xxxyy, ts_yzzzz_xxxyz, ts_yzzzz_xxxzz, ts_yzzzz_xxyyy, ts_yzzzz_xxyyz, ts_yzzzz_xxyzz, ts_yzzzz_xxzzz, ts_yzzzz_xyyyy, ts_yzzzz_xyyyz, ts_yzzzz_xyyzz, ts_yzzzz_xyzzz, ts_yzzzz_xzzzz, ts_yzzzz_yyyyy, ts_yzzzz_yyyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzzz, ts_yzzzz_zzzzz, ts_zzzz_xxxx, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxy, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxz, ts_zzzz_xxxzz, ts_zzzz_xxyy, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyz, ts_zzzz_xxyzz, ts_zzzz_xxzz, ts_zzzz_xxzzz, ts_zzzz_xyyy, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyz, ts_zzzz_xyyzz, ts_zzzz_xyzz, ts_zzzz_xyzzz, ts_zzzz_xzzz, ts_zzzz_xzzzz, ts_zzzz_yyyy, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyz, ts_zzzz_yyyzz, ts_zzzz_yyzz, ts_zzzz_yyzzz, ts_zzzz_yzzz, ts_zzzz_yzzzz, ts_zzzz_zzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzz_xxxxx[i] = ts_zzzz_xxxxx[i] * pa_y[i];

        ts_yzzzz_xxxxy[i] = ts_zzzz_xxxx[i] * fe_0 + ts_zzzz_xxxxy[i] * pa_y[i];

        ts_yzzzz_xxxxz[i] = ts_zzzz_xxxxz[i] * pa_y[i];

        ts_yzzzz_xxxyy[i] = 2.0 * ts_zzzz_xxxy[i] * fe_0 + ts_zzzz_xxxyy[i] * pa_y[i];

        ts_yzzzz_xxxyz[i] = ts_zzzz_xxxz[i] * fe_0 + ts_zzzz_xxxyz[i] * pa_y[i];

        ts_yzzzz_xxxzz[i] = ts_zzzz_xxxzz[i] * pa_y[i];

        ts_yzzzz_xxyyy[i] = 3.0 * ts_zzzz_xxyy[i] * fe_0 + ts_zzzz_xxyyy[i] * pa_y[i];

        ts_yzzzz_xxyyz[i] = 2.0 * ts_zzzz_xxyz[i] * fe_0 + ts_zzzz_xxyyz[i] * pa_y[i];

        ts_yzzzz_xxyzz[i] = ts_zzzz_xxzz[i] * fe_0 + ts_zzzz_xxyzz[i] * pa_y[i];

        ts_yzzzz_xxzzz[i] = ts_zzzz_xxzzz[i] * pa_y[i];

        ts_yzzzz_xyyyy[i] = 4.0 * ts_zzzz_xyyy[i] * fe_0 + ts_zzzz_xyyyy[i] * pa_y[i];

        ts_yzzzz_xyyyz[i] = 3.0 * ts_zzzz_xyyz[i] * fe_0 + ts_zzzz_xyyyz[i] * pa_y[i];

        ts_yzzzz_xyyzz[i] = 2.0 * ts_zzzz_xyzz[i] * fe_0 + ts_zzzz_xyyzz[i] * pa_y[i];

        ts_yzzzz_xyzzz[i] = ts_zzzz_xzzz[i] * fe_0 + ts_zzzz_xyzzz[i] * pa_y[i];

        ts_yzzzz_xzzzz[i] = ts_zzzz_xzzzz[i] * pa_y[i];

        ts_yzzzz_yyyyy[i] = 5.0 * ts_zzzz_yyyy[i] * fe_0 + ts_zzzz_yyyyy[i] * pa_y[i];

        ts_yzzzz_yyyyz[i] = 4.0 * ts_zzzz_yyyz[i] * fe_0 + ts_zzzz_yyyyz[i] * pa_y[i];

        ts_yzzzz_yyyzz[i] = 3.0 * ts_zzzz_yyzz[i] * fe_0 + ts_zzzz_yyyzz[i] * pa_y[i];

        ts_yzzzz_yyzzz[i] = 2.0 * ts_zzzz_yzzz[i] * fe_0 + ts_zzzz_yyzzz[i] * pa_y[i];

        ts_yzzzz_yzzzz[i] = ts_zzzz_zzzz[i] * fe_0 + ts_zzzz_yzzzz[i] * pa_y[i];

        ts_yzzzz_zzzzz[i] = ts_zzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 420-441 components of targeted buffer : HH

    auto ts_zzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 420);

    auto ts_zzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 421);

    auto ts_zzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 422);

    auto ts_zzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 423);

    auto ts_zzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 424);

    auto ts_zzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 425);

    auto ts_zzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 426);

    auto ts_zzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 427);

    auto ts_zzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 428);

    auto ts_zzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 429);

    auto ts_zzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 430);

    auto ts_zzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 431);

    auto ts_zzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 432);

    auto ts_zzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 433);

    auto ts_zzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 434);

    auto ts_zzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 435);

    auto ts_zzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 436);

    auto ts_zzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 437);

    auto ts_zzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 438);

    auto ts_zzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 439);

    auto ts_zzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 440);

    #pragma omp simd aligned(pa_z, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, ts_zzzz_xxxx, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxy, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxz, ts_zzzz_xxxzz, ts_zzzz_xxyy, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyz, ts_zzzz_xxyzz, ts_zzzz_xxzz, ts_zzzz_xxzzz, ts_zzzz_xyyy, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyz, ts_zzzz_xyyzz, ts_zzzz_xyzz, ts_zzzz_xyzzz, ts_zzzz_xzzz, ts_zzzz_xzzzz, ts_zzzz_yyyy, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyz, ts_zzzz_yyyzz, ts_zzzz_yyzz, ts_zzzz_yyzzz, ts_zzzz_yzzz, ts_zzzz_yzzzz, ts_zzzz_zzzz, ts_zzzz_zzzzz, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxzz, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzz_xxxxx[i] = 4.0 * ts_zzz_xxxxx[i] * fe_0 + ts_zzzz_xxxxx[i] * pa_z[i];

        ts_zzzzz_xxxxy[i] = 4.0 * ts_zzz_xxxxy[i] * fe_0 + ts_zzzz_xxxxy[i] * pa_z[i];

        ts_zzzzz_xxxxz[i] = 4.0 * ts_zzz_xxxxz[i] * fe_0 + ts_zzzz_xxxx[i] * fe_0 + ts_zzzz_xxxxz[i] * pa_z[i];

        ts_zzzzz_xxxyy[i] = 4.0 * ts_zzz_xxxyy[i] * fe_0 + ts_zzzz_xxxyy[i] * pa_z[i];

        ts_zzzzz_xxxyz[i] = 4.0 * ts_zzz_xxxyz[i] * fe_0 + ts_zzzz_xxxy[i] * fe_0 + ts_zzzz_xxxyz[i] * pa_z[i];

        ts_zzzzz_xxxzz[i] = 4.0 * ts_zzz_xxxzz[i] * fe_0 + 2.0 * ts_zzzz_xxxz[i] * fe_0 + ts_zzzz_xxxzz[i] * pa_z[i];

        ts_zzzzz_xxyyy[i] = 4.0 * ts_zzz_xxyyy[i] * fe_0 + ts_zzzz_xxyyy[i] * pa_z[i];

        ts_zzzzz_xxyyz[i] = 4.0 * ts_zzz_xxyyz[i] * fe_0 + ts_zzzz_xxyy[i] * fe_0 + ts_zzzz_xxyyz[i] * pa_z[i];

        ts_zzzzz_xxyzz[i] = 4.0 * ts_zzz_xxyzz[i] * fe_0 + 2.0 * ts_zzzz_xxyz[i] * fe_0 + ts_zzzz_xxyzz[i] * pa_z[i];

        ts_zzzzz_xxzzz[i] = 4.0 * ts_zzz_xxzzz[i] * fe_0 + 3.0 * ts_zzzz_xxzz[i] * fe_0 + ts_zzzz_xxzzz[i] * pa_z[i];

        ts_zzzzz_xyyyy[i] = 4.0 * ts_zzz_xyyyy[i] * fe_0 + ts_zzzz_xyyyy[i] * pa_z[i];

        ts_zzzzz_xyyyz[i] = 4.0 * ts_zzz_xyyyz[i] * fe_0 + ts_zzzz_xyyy[i] * fe_0 + ts_zzzz_xyyyz[i] * pa_z[i];

        ts_zzzzz_xyyzz[i] = 4.0 * ts_zzz_xyyzz[i] * fe_0 + 2.0 * ts_zzzz_xyyz[i] * fe_0 + ts_zzzz_xyyzz[i] * pa_z[i];

        ts_zzzzz_xyzzz[i] = 4.0 * ts_zzz_xyzzz[i] * fe_0 + 3.0 * ts_zzzz_xyzz[i] * fe_0 + ts_zzzz_xyzzz[i] * pa_z[i];

        ts_zzzzz_xzzzz[i] = 4.0 * ts_zzz_xzzzz[i] * fe_0 + 4.0 * ts_zzzz_xzzz[i] * fe_0 + ts_zzzz_xzzzz[i] * pa_z[i];

        ts_zzzzz_yyyyy[i] = 4.0 * ts_zzz_yyyyy[i] * fe_0 + ts_zzzz_yyyyy[i] * pa_z[i];

        ts_zzzzz_yyyyz[i] = 4.0 * ts_zzz_yyyyz[i] * fe_0 + ts_zzzz_yyyy[i] * fe_0 + ts_zzzz_yyyyz[i] * pa_z[i];

        ts_zzzzz_yyyzz[i] = 4.0 * ts_zzz_yyyzz[i] * fe_0 + 2.0 * ts_zzzz_yyyz[i] * fe_0 + ts_zzzz_yyyzz[i] * pa_z[i];

        ts_zzzzz_yyzzz[i] = 4.0 * ts_zzz_yyzzz[i] * fe_0 + 3.0 * ts_zzzz_yyzz[i] * fe_0 + ts_zzzz_yyzzz[i] * pa_z[i];

        ts_zzzzz_yzzzz[i] = 4.0 * ts_zzz_yzzzz[i] * fe_0 + 4.0 * ts_zzzz_yzzz[i] * fe_0 + ts_zzzz_yzzzz[i] * pa_z[i];

        ts_zzzzz_zzzzz[i] = 4.0 * ts_zzz_zzzzz[i] * fe_0 + 5.0 * ts_zzzz_zzzz[i] * fe_0 + ts_zzzz_zzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

