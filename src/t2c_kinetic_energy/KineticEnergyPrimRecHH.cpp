#include "KineticEnergyPrimRecHH.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_hh(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_hh,
                            const size_t              idx_ovl_fh,
                            const size_t              idx_kin_fh,
                            const size_t              idx_kin_gg,
                            const size_t              idx_kin_gh,
                            const size_t              idx_ovl_hh,
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

    auto ts_xxz_xxxxx = pbuffer.data(idx_ovl_fh + 42);

    auto ts_xxz_xxxxy = pbuffer.data(idx_ovl_fh + 43);

    auto ts_xxz_xxxyy = pbuffer.data(idx_ovl_fh + 45);

    auto ts_xxz_xxyyy = pbuffer.data(idx_ovl_fh + 48);

    auto ts_xxz_xyyyy = pbuffer.data(idx_ovl_fh + 52);

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

    auto ts_yyz_xxxyy = pbuffer.data(idx_ovl_fh + 150);

    auto ts_yyz_xxyyy = pbuffer.data(idx_ovl_fh + 153);

    auto ts_yyz_xyyyy = pbuffer.data(idx_ovl_fh + 157);

    auto ts_yyz_yyyyy = pbuffer.data(idx_ovl_fh + 162);

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

    // Set up components of auxiliary buffer : FH

    auto tk_xxx_xxxxx = pbuffer.data(idx_kin_fh);

    auto tk_xxx_xxxxy = pbuffer.data(idx_kin_fh + 1);

    auto tk_xxx_xxxxz = pbuffer.data(idx_kin_fh + 2);

    auto tk_xxx_xxxyy = pbuffer.data(idx_kin_fh + 3);

    auto tk_xxx_xxxyz = pbuffer.data(idx_kin_fh + 4);

    auto tk_xxx_xxxzz = pbuffer.data(idx_kin_fh + 5);

    auto tk_xxx_xxyyy = pbuffer.data(idx_kin_fh + 6);

    auto tk_xxx_xxyyz = pbuffer.data(idx_kin_fh + 7);

    auto tk_xxx_xxyzz = pbuffer.data(idx_kin_fh + 8);

    auto tk_xxx_xxzzz = pbuffer.data(idx_kin_fh + 9);

    auto tk_xxx_xyyyy = pbuffer.data(idx_kin_fh + 10);

    auto tk_xxx_xyyyz = pbuffer.data(idx_kin_fh + 11);

    auto tk_xxx_xyyzz = pbuffer.data(idx_kin_fh + 12);

    auto tk_xxx_xyzzz = pbuffer.data(idx_kin_fh + 13);

    auto tk_xxx_xzzzz = pbuffer.data(idx_kin_fh + 14);

    auto tk_xxx_yyyyy = pbuffer.data(idx_kin_fh + 15);

    auto tk_xxx_yyyyz = pbuffer.data(idx_kin_fh + 16);

    auto tk_xxx_yyyzz = pbuffer.data(idx_kin_fh + 17);

    auto tk_xxx_yyzzz = pbuffer.data(idx_kin_fh + 18);

    auto tk_xxx_yzzzz = pbuffer.data(idx_kin_fh + 19);

    auto tk_xxx_zzzzz = pbuffer.data(idx_kin_fh + 20);

    auto tk_xxy_xxxxx = pbuffer.data(idx_kin_fh + 21);

    auto tk_xxy_xxxxz = pbuffer.data(idx_kin_fh + 23);

    auto tk_xxy_xxxzz = pbuffer.data(idx_kin_fh + 26);

    auto tk_xxy_xxzzz = pbuffer.data(idx_kin_fh + 30);

    auto tk_xxy_xzzzz = pbuffer.data(idx_kin_fh + 35);

    auto tk_xxz_xxxxx = pbuffer.data(idx_kin_fh + 42);

    auto tk_xxz_xxxxy = pbuffer.data(idx_kin_fh + 43);

    auto tk_xxz_xxxyy = pbuffer.data(idx_kin_fh + 45);

    auto tk_xxz_xxyyy = pbuffer.data(idx_kin_fh + 48);

    auto tk_xxz_xyyyy = pbuffer.data(idx_kin_fh + 52);

    auto tk_xyy_xxxxy = pbuffer.data(idx_kin_fh + 64);

    auto tk_xyy_xxxyy = pbuffer.data(idx_kin_fh + 66);

    auto tk_xyy_xxxyz = pbuffer.data(idx_kin_fh + 67);

    auto tk_xyy_xxyyy = pbuffer.data(idx_kin_fh + 69);

    auto tk_xyy_xxyyz = pbuffer.data(idx_kin_fh + 70);

    auto tk_xyy_xxyzz = pbuffer.data(idx_kin_fh + 71);

    auto tk_xyy_xyyyy = pbuffer.data(idx_kin_fh + 73);

    auto tk_xyy_xyyyz = pbuffer.data(idx_kin_fh + 74);

    auto tk_xyy_xyyzz = pbuffer.data(idx_kin_fh + 75);

    auto tk_xyy_xyzzz = pbuffer.data(idx_kin_fh + 76);

    auto tk_xyy_yyyyy = pbuffer.data(idx_kin_fh + 78);

    auto tk_xyy_yyyyz = pbuffer.data(idx_kin_fh + 79);

    auto tk_xyy_yyyzz = pbuffer.data(idx_kin_fh + 80);

    auto tk_xyy_yyzzz = pbuffer.data(idx_kin_fh + 81);

    auto tk_xyy_yzzzz = pbuffer.data(idx_kin_fh + 82);

    auto tk_xyy_zzzzz = pbuffer.data(idx_kin_fh + 83);

    auto tk_xzz_xxxxz = pbuffer.data(idx_kin_fh + 107);

    auto tk_xzz_xxxyz = pbuffer.data(idx_kin_fh + 109);

    auto tk_xzz_xxxzz = pbuffer.data(idx_kin_fh + 110);

    auto tk_xzz_xxyyz = pbuffer.data(idx_kin_fh + 112);

    auto tk_xzz_xxyzz = pbuffer.data(idx_kin_fh + 113);

    auto tk_xzz_xxzzz = pbuffer.data(idx_kin_fh + 114);

    auto tk_xzz_xyyyz = pbuffer.data(idx_kin_fh + 116);

    auto tk_xzz_xyyzz = pbuffer.data(idx_kin_fh + 117);

    auto tk_xzz_xyzzz = pbuffer.data(idx_kin_fh + 118);

    auto tk_xzz_xzzzz = pbuffer.data(idx_kin_fh + 119);

    auto tk_xzz_yyyyy = pbuffer.data(idx_kin_fh + 120);

    auto tk_xzz_yyyyz = pbuffer.data(idx_kin_fh + 121);

    auto tk_xzz_yyyzz = pbuffer.data(idx_kin_fh + 122);

    auto tk_xzz_yyzzz = pbuffer.data(idx_kin_fh + 123);

    auto tk_xzz_yzzzz = pbuffer.data(idx_kin_fh + 124);

    auto tk_xzz_zzzzz = pbuffer.data(idx_kin_fh + 125);

    auto tk_yyy_xxxxx = pbuffer.data(idx_kin_fh + 126);

    auto tk_yyy_xxxxy = pbuffer.data(idx_kin_fh + 127);

    auto tk_yyy_xxxxz = pbuffer.data(idx_kin_fh + 128);

    auto tk_yyy_xxxyy = pbuffer.data(idx_kin_fh + 129);

    auto tk_yyy_xxxyz = pbuffer.data(idx_kin_fh + 130);

    auto tk_yyy_xxxzz = pbuffer.data(idx_kin_fh + 131);

    auto tk_yyy_xxyyy = pbuffer.data(idx_kin_fh + 132);

    auto tk_yyy_xxyyz = pbuffer.data(idx_kin_fh + 133);

    auto tk_yyy_xxyzz = pbuffer.data(idx_kin_fh + 134);

    auto tk_yyy_xxzzz = pbuffer.data(idx_kin_fh + 135);

    auto tk_yyy_xyyyy = pbuffer.data(idx_kin_fh + 136);

    auto tk_yyy_xyyyz = pbuffer.data(idx_kin_fh + 137);

    auto tk_yyy_xyyzz = pbuffer.data(idx_kin_fh + 138);

    auto tk_yyy_xyzzz = pbuffer.data(idx_kin_fh + 139);

    auto tk_yyy_xzzzz = pbuffer.data(idx_kin_fh + 140);

    auto tk_yyy_yyyyy = pbuffer.data(idx_kin_fh + 141);

    auto tk_yyy_yyyyz = pbuffer.data(idx_kin_fh + 142);

    auto tk_yyy_yyyzz = pbuffer.data(idx_kin_fh + 143);

    auto tk_yyy_yyzzz = pbuffer.data(idx_kin_fh + 144);

    auto tk_yyy_yzzzz = pbuffer.data(idx_kin_fh + 145);

    auto tk_yyy_zzzzz = pbuffer.data(idx_kin_fh + 146);

    auto tk_yyz_xxxxy = pbuffer.data(idx_kin_fh + 148);

    auto tk_yyz_xxxyy = pbuffer.data(idx_kin_fh + 150);

    auto tk_yyz_xxyyy = pbuffer.data(idx_kin_fh + 153);

    auto tk_yyz_xyyyy = pbuffer.data(idx_kin_fh + 157);

    auto tk_yyz_yyyyy = pbuffer.data(idx_kin_fh + 162);

    auto tk_yzz_xxxxx = pbuffer.data(idx_kin_fh + 168);

    auto tk_yzz_xxxxz = pbuffer.data(idx_kin_fh + 170);

    auto tk_yzz_xxxyz = pbuffer.data(idx_kin_fh + 172);

    auto tk_yzz_xxxzz = pbuffer.data(idx_kin_fh + 173);

    auto tk_yzz_xxyyz = pbuffer.data(idx_kin_fh + 175);

    auto tk_yzz_xxyzz = pbuffer.data(idx_kin_fh + 176);

    auto tk_yzz_xxzzz = pbuffer.data(idx_kin_fh + 177);

    auto tk_yzz_xyyyz = pbuffer.data(idx_kin_fh + 179);

    auto tk_yzz_xyyzz = pbuffer.data(idx_kin_fh + 180);

    auto tk_yzz_xyzzz = pbuffer.data(idx_kin_fh + 181);

    auto tk_yzz_xzzzz = pbuffer.data(idx_kin_fh + 182);

    auto tk_yzz_yyyyz = pbuffer.data(idx_kin_fh + 184);

    auto tk_yzz_yyyzz = pbuffer.data(idx_kin_fh + 185);

    auto tk_yzz_yyzzz = pbuffer.data(idx_kin_fh + 186);

    auto tk_yzz_yzzzz = pbuffer.data(idx_kin_fh + 187);

    auto tk_yzz_zzzzz = pbuffer.data(idx_kin_fh + 188);

    auto tk_zzz_xxxxx = pbuffer.data(idx_kin_fh + 189);

    auto tk_zzz_xxxxy = pbuffer.data(idx_kin_fh + 190);

    auto tk_zzz_xxxxz = pbuffer.data(idx_kin_fh + 191);

    auto tk_zzz_xxxyy = pbuffer.data(idx_kin_fh + 192);

    auto tk_zzz_xxxyz = pbuffer.data(idx_kin_fh + 193);

    auto tk_zzz_xxxzz = pbuffer.data(idx_kin_fh + 194);

    auto tk_zzz_xxyyy = pbuffer.data(idx_kin_fh + 195);

    auto tk_zzz_xxyyz = pbuffer.data(idx_kin_fh + 196);

    auto tk_zzz_xxyzz = pbuffer.data(idx_kin_fh + 197);

    auto tk_zzz_xxzzz = pbuffer.data(idx_kin_fh + 198);

    auto tk_zzz_xyyyy = pbuffer.data(idx_kin_fh + 199);

    auto tk_zzz_xyyyz = pbuffer.data(idx_kin_fh + 200);

    auto tk_zzz_xyyzz = pbuffer.data(idx_kin_fh + 201);

    auto tk_zzz_xyzzz = pbuffer.data(idx_kin_fh + 202);

    auto tk_zzz_xzzzz = pbuffer.data(idx_kin_fh + 203);

    auto tk_zzz_yyyyy = pbuffer.data(idx_kin_fh + 204);

    auto tk_zzz_yyyyz = pbuffer.data(idx_kin_fh + 205);

    auto tk_zzz_yyyzz = pbuffer.data(idx_kin_fh + 206);

    auto tk_zzz_yyzzz = pbuffer.data(idx_kin_fh + 207);

    auto tk_zzz_yzzzz = pbuffer.data(idx_kin_fh + 208);

    auto tk_zzz_zzzzz = pbuffer.data(idx_kin_fh + 209);

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

    auto tk_xxxz_xxxz = pbuffer.data(idx_kin_gg + 32);

    auto tk_xxxz_xxyz = pbuffer.data(idx_kin_gg + 34);

    auto tk_xxxz_xxzz = pbuffer.data(idx_kin_gg + 35);

    auto tk_xxxz_xyyz = pbuffer.data(idx_kin_gg + 37);

    auto tk_xxxz_xyzz = pbuffer.data(idx_kin_gg + 38);

    auto tk_xxxz_xzzz = pbuffer.data(idx_kin_gg + 39);

    auto tk_xxxz_yyyz = pbuffer.data(idx_kin_gg + 41);

    auto tk_xxxz_yyzz = pbuffer.data(idx_kin_gg + 42);

    auto tk_xxxz_yzzz = pbuffer.data(idx_kin_gg + 43);

    auto tk_xxxz_zzzz = pbuffer.data(idx_kin_gg + 44);

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

    auto tk_xzzz_xxxz = pbuffer.data(idx_kin_gg + 137);

    auto tk_xzzz_xxyz = pbuffer.data(idx_kin_gg + 139);

    auto tk_xzzz_xxzz = pbuffer.data(idx_kin_gg + 140);

    auto tk_xzzz_xyyz = pbuffer.data(idx_kin_gg + 142);

    auto tk_xzzz_xyzz = pbuffer.data(idx_kin_gg + 143);

    auto tk_xzzz_xzzz = pbuffer.data(idx_kin_gg + 144);

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

    auto tk_yyyz_xxxz = pbuffer.data(idx_kin_gg + 167);

    auto tk_yyyz_xxyz = pbuffer.data(idx_kin_gg + 169);

    auto tk_yyyz_xxzz = pbuffer.data(idx_kin_gg + 170);

    auto tk_yyyz_xyyz = pbuffer.data(idx_kin_gg + 172);

    auto tk_yyyz_xyzz = pbuffer.data(idx_kin_gg + 173);

    auto tk_yyyz_xzzz = pbuffer.data(idx_kin_gg + 174);

    auto tk_yyyz_yyyz = pbuffer.data(idx_kin_gg + 176);

    auto tk_yyyz_yyzz = pbuffer.data(idx_kin_gg + 177);

    auto tk_yyyz_yzzz = pbuffer.data(idx_kin_gg + 178);

    auto tk_yyyz_zzzz = pbuffer.data(idx_kin_gg + 179);

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

    auto tk_yzzz_xxxy = pbuffer.data(idx_kin_gg + 196);

    auto tk_yzzz_xxxz = pbuffer.data(idx_kin_gg + 197);

    auto tk_yzzz_xxyy = pbuffer.data(idx_kin_gg + 198);

    auto tk_yzzz_xxyz = pbuffer.data(idx_kin_gg + 199);

    auto tk_yzzz_xxzz = pbuffer.data(idx_kin_gg + 200);

    auto tk_yzzz_xyyy = pbuffer.data(idx_kin_gg + 201);

    auto tk_yzzz_xyyz = pbuffer.data(idx_kin_gg + 202);

    auto tk_yzzz_xyzz = pbuffer.data(idx_kin_gg + 203);

    auto tk_yzzz_xzzz = pbuffer.data(idx_kin_gg + 204);

    auto tk_yzzz_yyyy = pbuffer.data(idx_kin_gg + 205);

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

    // Set up components of auxiliary buffer : GH

    auto tk_xxxx_xxxxx = pbuffer.data(idx_kin_gh);

    auto tk_xxxx_xxxxy = pbuffer.data(idx_kin_gh + 1);

    auto tk_xxxx_xxxxz = pbuffer.data(idx_kin_gh + 2);

    auto tk_xxxx_xxxyy = pbuffer.data(idx_kin_gh + 3);

    auto tk_xxxx_xxxyz = pbuffer.data(idx_kin_gh + 4);

    auto tk_xxxx_xxxzz = pbuffer.data(idx_kin_gh + 5);

    auto tk_xxxx_xxyyy = pbuffer.data(idx_kin_gh + 6);

    auto tk_xxxx_xxyyz = pbuffer.data(idx_kin_gh + 7);

    auto tk_xxxx_xxyzz = pbuffer.data(idx_kin_gh + 8);

    auto tk_xxxx_xxzzz = pbuffer.data(idx_kin_gh + 9);

    auto tk_xxxx_xyyyy = pbuffer.data(idx_kin_gh + 10);

    auto tk_xxxx_xyyyz = pbuffer.data(idx_kin_gh + 11);

    auto tk_xxxx_xyyzz = pbuffer.data(idx_kin_gh + 12);

    auto tk_xxxx_xyzzz = pbuffer.data(idx_kin_gh + 13);

    auto tk_xxxx_xzzzz = pbuffer.data(idx_kin_gh + 14);

    auto tk_xxxx_yyyyy = pbuffer.data(idx_kin_gh + 15);

    auto tk_xxxx_yyyyz = pbuffer.data(idx_kin_gh + 16);

    auto tk_xxxx_yyyzz = pbuffer.data(idx_kin_gh + 17);

    auto tk_xxxx_yyzzz = pbuffer.data(idx_kin_gh + 18);

    auto tk_xxxx_yzzzz = pbuffer.data(idx_kin_gh + 19);

    auto tk_xxxx_zzzzz = pbuffer.data(idx_kin_gh + 20);

    auto tk_xxxy_xxxxx = pbuffer.data(idx_kin_gh + 21);

    auto tk_xxxy_xxxxy = pbuffer.data(idx_kin_gh + 22);

    auto tk_xxxy_xxxxz = pbuffer.data(idx_kin_gh + 23);

    auto tk_xxxy_xxxyy = pbuffer.data(idx_kin_gh + 24);

    auto tk_xxxy_xxxzz = pbuffer.data(idx_kin_gh + 26);

    auto tk_xxxy_xxyyy = pbuffer.data(idx_kin_gh + 27);

    auto tk_xxxy_xxzzz = pbuffer.data(idx_kin_gh + 30);

    auto tk_xxxy_xyyyy = pbuffer.data(idx_kin_gh + 31);

    auto tk_xxxy_xzzzz = pbuffer.data(idx_kin_gh + 35);

    auto tk_xxxy_yyyyy = pbuffer.data(idx_kin_gh + 36);

    auto tk_xxxz_xxxxx = pbuffer.data(idx_kin_gh + 42);

    auto tk_xxxz_xxxxy = pbuffer.data(idx_kin_gh + 43);

    auto tk_xxxz_xxxxz = pbuffer.data(idx_kin_gh + 44);

    auto tk_xxxz_xxxyy = pbuffer.data(idx_kin_gh + 45);

    auto tk_xxxz_xxxyz = pbuffer.data(idx_kin_gh + 46);

    auto tk_xxxz_xxxzz = pbuffer.data(idx_kin_gh + 47);

    auto tk_xxxz_xxyyy = pbuffer.data(idx_kin_gh + 48);

    auto tk_xxxz_xxyyz = pbuffer.data(idx_kin_gh + 49);

    auto tk_xxxz_xxyzz = pbuffer.data(idx_kin_gh + 50);

    auto tk_xxxz_xxzzz = pbuffer.data(idx_kin_gh + 51);

    auto tk_xxxz_xyyyy = pbuffer.data(idx_kin_gh + 52);

    auto tk_xxxz_xyyyz = pbuffer.data(idx_kin_gh + 53);

    auto tk_xxxz_xyyzz = pbuffer.data(idx_kin_gh + 54);

    auto tk_xxxz_xyzzz = pbuffer.data(idx_kin_gh + 55);

    auto tk_xxxz_xzzzz = pbuffer.data(idx_kin_gh + 56);

    auto tk_xxxz_yyyyz = pbuffer.data(idx_kin_gh + 58);

    auto tk_xxxz_yyyzz = pbuffer.data(idx_kin_gh + 59);

    auto tk_xxxz_yyzzz = pbuffer.data(idx_kin_gh + 60);

    auto tk_xxxz_yzzzz = pbuffer.data(idx_kin_gh + 61);

    auto tk_xxxz_zzzzz = pbuffer.data(idx_kin_gh + 62);

    auto tk_xxyy_xxxxx = pbuffer.data(idx_kin_gh + 63);

    auto tk_xxyy_xxxxy = pbuffer.data(idx_kin_gh + 64);

    auto tk_xxyy_xxxxz = pbuffer.data(idx_kin_gh + 65);

    auto tk_xxyy_xxxyy = pbuffer.data(idx_kin_gh + 66);

    auto tk_xxyy_xxxyz = pbuffer.data(idx_kin_gh + 67);

    auto tk_xxyy_xxxzz = pbuffer.data(idx_kin_gh + 68);

    auto tk_xxyy_xxyyy = pbuffer.data(idx_kin_gh + 69);

    auto tk_xxyy_xxyyz = pbuffer.data(idx_kin_gh + 70);

    auto tk_xxyy_xxyzz = pbuffer.data(idx_kin_gh + 71);

    auto tk_xxyy_xxzzz = pbuffer.data(idx_kin_gh + 72);

    auto tk_xxyy_xyyyy = pbuffer.data(idx_kin_gh + 73);

    auto tk_xxyy_xyyyz = pbuffer.data(idx_kin_gh + 74);

    auto tk_xxyy_xyyzz = pbuffer.data(idx_kin_gh + 75);

    auto tk_xxyy_xyzzz = pbuffer.data(idx_kin_gh + 76);

    auto tk_xxyy_xzzzz = pbuffer.data(idx_kin_gh + 77);

    auto tk_xxyy_yyyyy = pbuffer.data(idx_kin_gh + 78);

    auto tk_xxyy_yyyyz = pbuffer.data(idx_kin_gh + 79);

    auto tk_xxyy_yyyzz = pbuffer.data(idx_kin_gh + 80);

    auto tk_xxyy_yyzzz = pbuffer.data(idx_kin_gh + 81);

    auto tk_xxyy_yzzzz = pbuffer.data(idx_kin_gh + 82);

    auto tk_xxyy_zzzzz = pbuffer.data(idx_kin_gh + 83);

    auto tk_xxzz_xxxxx = pbuffer.data(idx_kin_gh + 105);

    auto tk_xxzz_xxxxy = pbuffer.data(idx_kin_gh + 106);

    auto tk_xxzz_xxxxz = pbuffer.data(idx_kin_gh + 107);

    auto tk_xxzz_xxxyy = pbuffer.data(idx_kin_gh + 108);

    auto tk_xxzz_xxxyz = pbuffer.data(idx_kin_gh + 109);

    auto tk_xxzz_xxxzz = pbuffer.data(idx_kin_gh + 110);

    auto tk_xxzz_xxyyy = pbuffer.data(idx_kin_gh + 111);

    auto tk_xxzz_xxyyz = pbuffer.data(idx_kin_gh + 112);

    auto tk_xxzz_xxyzz = pbuffer.data(idx_kin_gh + 113);

    auto tk_xxzz_xxzzz = pbuffer.data(idx_kin_gh + 114);

    auto tk_xxzz_xyyyy = pbuffer.data(idx_kin_gh + 115);

    auto tk_xxzz_xyyyz = pbuffer.data(idx_kin_gh + 116);

    auto tk_xxzz_xyyzz = pbuffer.data(idx_kin_gh + 117);

    auto tk_xxzz_xyzzz = pbuffer.data(idx_kin_gh + 118);

    auto tk_xxzz_xzzzz = pbuffer.data(idx_kin_gh + 119);

    auto tk_xxzz_yyyyy = pbuffer.data(idx_kin_gh + 120);

    auto tk_xxzz_yyyyz = pbuffer.data(idx_kin_gh + 121);

    auto tk_xxzz_yyyzz = pbuffer.data(idx_kin_gh + 122);

    auto tk_xxzz_yyzzz = pbuffer.data(idx_kin_gh + 123);

    auto tk_xxzz_yzzzz = pbuffer.data(idx_kin_gh + 124);

    auto tk_xxzz_zzzzz = pbuffer.data(idx_kin_gh + 125);

    auto tk_xyyy_xxxxx = pbuffer.data(idx_kin_gh + 126);

    auto tk_xyyy_xxxxy = pbuffer.data(idx_kin_gh + 127);

    auto tk_xyyy_xxxyy = pbuffer.data(idx_kin_gh + 129);

    auto tk_xyyy_xxxyz = pbuffer.data(idx_kin_gh + 130);

    auto tk_xyyy_xxyyy = pbuffer.data(idx_kin_gh + 132);

    auto tk_xyyy_xxyyz = pbuffer.data(idx_kin_gh + 133);

    auto tk_xyyy_xxyzz = pbuffer.data(idx_kin_gh + 134);

    auto tk_xyyy_xyyyy = pbuffer.data(idx_kin_gh + 136);

    auto tk_xyyy_xyyyz = pbuffer.data(idx_kin_gh + 137);

    auto tk_xyyy_xyyzz = pbuffer.data(idx_kin_gh + 138);

    auto tk_xyyy_xyzzz = pbuffer.data(idx_kin_gh + 139);

    auto tk_xyyy_yyyyy = pbuffer.data(idx_kin_gh + 141);

    auto tk_xyyy_yyyyz = pbuffer.data(idx_kin_gh + 142);

    auto tk_xyyy_yyyzz = pbuffer.data(idx_kin_gh + 143);

    auto tk_xyyy_yyzzz = pbuffer.data(idx_kin_gh + 144);

    auto tk_xyyy_yzzzz = pbuffer.data(idx_kin_gh + 145);

    auto tk_xyyy_zzzzz = pbuffer.data(idx_kin_gh + 146);

    auto tk_xzzz_xxxxx = pbuffer.data(idx_kin_gh + 189);

    auto tk_xzzz_xxxxz = pbuffer.data(idx_kin_gh + 191);

    auto tk_xzzz_xxxyz = pbuffer.data(idx_kin_gh + 193);

    auto tk_xzzz_xxxzz = pbuffer.data(idx_kin_gh + 194);

    auto tk_xzzz_xxyyz = pbuffer.data(idx_kin_gh + 196);

    auto tk_xzzz_xxyzz = pbuffer.data(idx_kin_gh + 197);

    auto tk_xzzz_xxzzz = pbuffer.data(idx_kin_gh + 198);

    auto tk_xzzz_xyyyz = pbuffer.data(idx_kin_gh + 200);

    auto tk_xzzz_xyyzz = pbuffer.data(idx_kin_gh + 201);

    auto tk_xzzz_xyzzz = pbuffer.data(idx_kin_gh + 202);

    auto tk_xzzz_xzzzz = pbuffer.data(idx_kin_gh + 203);

    auto tk_xzzz_yyyyy = pbuffer.data(idx_kin_gh + 204);

    auto tk_xzzz_yyyyz = pbuffer.data(idx_kin_gh + 205);

    auto tk_xzzz_yyyzz = pbuffer.data(idx_kin_gh + 206);

    auto tk_xzzz_yyzzz = pbuffer.data(idx_kin_gh + 207);

    auto tk_xzzz_yzzzz = pbuffer.data(idx_kin_gh + 208);

    auto tk_xzzz_zzzzz = pbuffer.data(idx_kin_gh + 209);

    auto tk_yyyy_xxxxx = pbuffer.data(idx_kin_gh + 210);

    auto tk_yyyy_xxxxy = pbuffer.data(idx_kin_gh + 211);

    auto tk_yyyy_xxxxz = pbuffer.data(idx_kin_gh + 212);

    auto tk_yyyy_xxxyy = pbuffer.data(idx_kin_gh + 213);

    auto tk_yyyy_xxxyz = pbuffer.data(idx_kin_gh + 214);

    auto tk_yyyy_xxxzz = pbuffer.data(idx_kin_gh + 215);

    auto tk_yyyy_xxyyy = pbuffer.data(idx_kin_gh + 216);

    auto tk_yyyy_xxyyz = pbuffer.data(idx_kin_gh + 217);

    auto tk_yyyy_xxyzz = pbuffer.data(idx_kin_gh + 218);

    auto tk_yyyy_xxzzz = pbuffer.data(idx_kin_gh + 219);

    auto tk_yyyy_xyyyy = pbuffer.data(idx_kin_gh + 220);

    auto tk_yyyy_xyyyz = pbuffer.data(idx_kin_gh + 221);

    auto tk_yyyy_xyyzz = pbuffer.data(idx_kin_gh + 222);

    auto tk_yyyy_xyzzz = pbuffer.data(idx_kin_gh + 223);

    auto tk_yyyy_xzzzz = pbuffer.data(idx_kin_gh + 224);

    auto tk_yyyy_yyyyy = pbuffer.data(idx_kin_gh + 225);

    auto tk_yyyy_yyyyz = pbuffer.data(idx_kin_gh + 226);

    auto tk_yyyy_yyyzz = pbuffer.data(idx_kin_gh + 227);

    auto tk_yyyy_yyzzz = pbuffer.data(idx_kin_gh + 228);

    auto tk_yyyy_yzzzz = pbuffer.data(idx_kin_gh + 229);

    auto tk_yyyy_zzzzz = pbuffer.data(idx_kin_gh + 230);

    auto tk_yyyz_xxxxy = pbuffer.data(idx_kin_gh + 232);

    auto tk_yyyz_xxxxz = pbuffer.data(idx_kin_gh + 233);

    auto tk_yyyz_xxxyy = pbuffer.data(idx_kin_gh + 234);

    auto tk_yyyz_xxxyz = pbuffer.data(idx_kin_gh + 235);

    auto tk_yyyz_xxxzz = pbuffer.data(idx_kin_gh + 236);

    auto tk_yyyz_xxyyy = pbuffer.data(idx_kin_gh + 237);

    auto tk_yyyz_xxyyz = pbuffer.data(idx_kin_gh + 238);

    auto tk_yyyz_xxyzz = pbuffer.data(idx_kin_gh + 239);

    auto tk_yyyz_xxzzz = pbuffer.data(idx_kin_gh + 240);

    auto tk_yyyz_xyyyy = pbuffer.data(idx_kin_gh + 241);

    auto tk_yyyz_xyyyz = pbuffer.data(idx_kin_gh + 242);

    auto tk_yyyz_xyyzz = pbuffer.data(idx_kin_gh + 243);

    auto tk_yyyz_xyzzz = pbuffer.data(idx_kin_gh + 244);

    auto tk_yyyz_xzzzz = pbuffer.data(idx_kin_gh + 245);

    auto tk_yyyz_yyyyy = pbuffer.data(idx_kin_gh + 246);

    auto tk_yyyz_yyyyz = pbuffer.data(idx_kin_gh + 247);

    auto tk_yyyz_yyyzz = pbuffer.data(idx_kin_gh + 248);

    auto tk_yyyz_yyzzz = pbuffer.data(idx_kin_gh + 249);

    auto tk_yyyz_yzzzz = pbuffer.data(idx_kin_gh + 250);

    auto tk_yyyz_zzzzz = pbuffer.data(idx_kin_gh + 251);

    auto tk_yyzz_xxxxx = pbuffer.data(idx_kin_gh + 252);

    auto tk_yyzz_xxxxy = pbuffer.data(idx_kin_gh + 253);

    auto tk_yyzz_xxxxz = pbuffer.data(idx_kin_gh + 254);

    auto tk_yyzz_xxxyy = pbuffer.data(idx_kin_gh + 255);

    auto tk_yyzz_xxxyz = pbuffer.data(idx_kin_gh + 256);

    auto tk_yyzz_xxxzz = pbuffer.data(idx_kin_gh + 257);

    auto tk_yyzz_xxyyy = pbuffer.data(idx_kin_gh + 258);

    auto tk_yyzz_xxyyz = pbuffer.data(idx_kin_gh + 259);

    auto tk_yyzz_xxyzz = pbuffer.data(idx_kin_gh + 260);

    auto tk_yyzz_xxzzz = pbuffer.data(idx_kin_gh + 261);

    auto tk_yyzz_xyyyy = pbuffer.data(idx_kin_gh + 262);

    auto tk_yyzz_xyyyz = pbuffer.data(idx_kin_gh + 263);

    auto tk_yyzz_xyyzz = pbuffer.data(idx_kin_gh + 264);

    auto tk_yyzz_xyzzz = pbuffer.data(idx_kin_gh + 265);

    auto tk_yyzz_xzzzz = pbuffer.data(idx_kin_gh + 266);

    auto tk_yyzz_yyyyy = pbuffer.data(idx_kin_gh + 267);

    auto tk_yyzz_yyyyz = pbuffer.data(idx_kin_gh + 268);

    auto tk_yyzz_yyyzz = pbuffer.data(idx_kin_gh + 269);

    auto tk_yyzz_yyzzz = pbuffer.data(idx_kin_gh + 270);

    auto tk_yyzz_yzzzz = pbuffer.data(idx_kin_gh + 271);

    auto tk_yyzz_zzzzz = pbuffer.data(idx_kin_gh + 272);

    auto tk_yzzz_xxxxx = pbuffer.data(idx_kin_gh + 273);

    auto tk_yzzz_xxxxy = pbuffer.data(idx_kin_gh + 274);

    auto tk_yzzz_xxxxz = pbuffer.data(idx_kin_gh + 275);

    auto tk_yzzz_xxxyy = pbuffer.data(idx_kin_gh + 276);

    auto tk_yzzz_xxxyz = pbuffer.data(idx_kin_gh + 277);

    auto tk_yzzz_xxxzz = pbuffer.data(idx_kin_gh + 278);

    auto tk_yzzz_xxyyy = pbuffer.data(idx_kin_gh + 279);

    auto tk_yzzz_xxyyz = pbuffer.data(idx_kin_gh + 280);

    auto tk_yzzz_xxyzz = pbuffer.data(idx_kin_gh + 281);

    auto tk_yzzz_xxzzz = pbuffer.data(idx_kin_gh + 282);

    auto tk_yzzz_xyyyy = pbuffer.data(idx_kin_gh + 283);

    auto tk_yzzz_xyyyz = pbuffer.data(idx_kin_gh + 284);

    auto tk_yzzz_xyyzz = pbuffer.data(idx_kin_gh + 285);

    auto tk_yzzz_xyzzz = pbuffer.data(idx_kin_gh + 286);

    auto tk_yzzz_xzzzz = pbuffer.data(idx_kin_gh + 287);

    auto tk_yzzz_yyyyy = pbuffer.data(idx_kin_gh + 288);

    auto tk_yzzz_yyyyz = pbuffer.data(idx_kin_gh + 289);

    auto tk_yzzz_yyyzz = pbuffer.data(idx_kin_gh + 290);

    auto tk_yzzz_yyzzz = pbuffer.data(idx_kin_gh + 291);

    auto tk_yzzz_yzzzz = pbuffer.data(idx_kin_gh + 292);

    auto tk_yzzz_zzzzz = pbuffer.data(idx_kin_gh + 293);

    auto tk_zzzz_xxxxx = pbuffer.data(idx_kin_gh + 294);

    auto tk_zzzz_xxxxy = pbuffer.data(idx_kin_gh + 295);

    auto tk_zzzz_xxxxz = pbuffer.data(idx_kin_gh + 296);

    auto tk_zzzz_xxxyy = pbuffer.data(idx_kin_gh + 297);

    auto tk_zzzz_xxxyz = pbuffer.data(idx_kin_gh + 298);

    auto tk_zzzz_xxxzz = pbuffer.data(idx_kin_gh + 299);

    auto tk_zzzz_xxyyy = pbuffer.data(idx_kin_gh + 300);

    auto tk_zzzz_xxyyz = pbuffer.data(idx_kin_gh + 301);

    auto tk_zzzz_xxyzz = pbuffer.data(idx_kin_gh + 302);

    auto tk_zzzz_xxzzz = pbuffer.data(idx_kin_gh + 303);

    auto tk_zzzz_xyyyy = pbuffer.data(idx_kin_gh + 304);

    auto tk_zzzz_xyyyz = pbuffer.data(idx_kin_gh + 305);

    auto tk_zzzz_xyyzz = pbuffer.data(idx_kin_gh + 306);

    auto tk_zzzz_xyzzz = pbuffer.data(idx_kin_gh + 307);

    auto tk_zzzz_xzzzz = pbuffer.data(idx_kin_gh + 308);

    auto tk_zzzz_yyyyy = pbuffer.data(idx_kin_gh + 309);

    auto tk_zzzz_yyyyz = pbuffer.data(idx_kin_gh + 310);

    auto tk_zzzz_yyyzz = pbuffer.data(idx_kin_gh + 311);

    auto tk_zzzz_yyzzz = pbuffer.data(idx_kin_gh + 312);

    auto tk_zzzz_yzzzz = pbuffer.data(idx_kin_gh + 313);

    auto tk_zzzz_zzzzz = pbuffer.data(idx_kin_gh + 314);

    // Set up components of auxiliary buffer : HH

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

    // Set up 0-21 components of targeted buffer : HH

    auto tk_xxxxx_xxxxx = pbuffer.data(idx_kin_hh);

    auto tk_xxxxx_xxxxy = pbuffer.data(idx_kin_hh + 1);

    auto tk_xxxxx_xxxxz = pbuffer.data(idx_kin_hh + 2);

    auto tk_xxxxx_xxxyy = pbuffer.data(idx_kin_hh + 3);

    auto tk_xxxxx_xxxyz = pbuffer.data(idx_kin_hh + 4);

    auto tk_xxxxx_xxxzz = pbuffer.data(idx_kin_hh + 5);

    auto tk_xxxxx_xxyyy = pbuffer.data(idx_kin_hh + 6);

    auto tk_xxxxx_xxyyz = pbuffer.data(idx_kin_hh + 7);

    auto tk_xxxxx_xxyzz = pbuffer.data(idx_kin_hh + 8);

    auto tk_xxxxx_xxzzz = pbuffer.data(idx_kin_hh + 9);

    auto tk_xxxxx_xyyyy = pbuffer.data(idx_kin_hh + 10);

    auto tk_xxxxx_xyyyz = pbuffer.data(idx_kin_hh + 11);

    auto tk_xxxxx_xyyzz = pbuffer.data(idx_kin_hh + 12);

    auto tk_xxxxx_xyzzz = pbuffer.data(idx_kin_hh + 13);

    auto tk_xxxxx_xzzzz = pbuffer.data(idx_kin_hh + 14);

    auto tk_xxxxx_yyyyy = pbuffer.data(idx_kin_hh + 15);

    auto tk_xxxxx_yyyyz = pbuffer.data(idx_kin_hh + 16);

    auto tk_xxxxx_yyyzz = pbuffer.data(idx_kin_hh + 17);

    auto tk_xxxxx_yyzzz = pbuffer.data(idx_kin_hh + 18);

    auto tk_xxxxx_yzzzz = pbuffer.data(idx_kin_hh + 19);

    auto tk_xxxxx_zzzzz = pbuffer.data(idx_kin_hh + 20);

#pragma omp simd aligned(pa_x,               \
                             tk_xxx_xxxxx,   \
                             tk_xxx_xxxxy,   \
                             tk_xxx_xxxxz,   \
                             tk_xxx_xxxyy,   \
                             tk_xxx_xxxyz,   \
                             tk_xxx_xxxzz,   \
                             tk_xxx_xxyyy,   \
                             tk_xxx_xxyyz,   \
                             tk_xxx_xxyzz,   \
                             tk_xxx_xxzzz,   \
                             tk_xxx_xyyyy,   \
                             tk_xxx_xyyyz,   \
                             tk_xxx_xyyzz,   \
                             tk_xxx_xyzzz,   \
                             tk_xxx_xzzzz,   \
                             tk_xxx_yyyyy,   \
                             tk_xxx_yyyyz,   \
                             tk_xxx_yyyzz,   \
                             tk_xxx_yyzzz,   \
                             tk_xxx_yzzzz,   \
                             tk_xxx_zzzzz,   \
                             tk_xxxx_xxxx,   \
                             tk_xxxx_xxxxx,  \
                             tk_xxxx_xxxxy,  \
                             tk_xxxx_xxxxz,  \
                             tk_xxxx_xxxy,   \
                             tk_xxxx_xxxyy,  \
                             tk_xxxx_xxxyz,  \
                             tk_xxxx_xxxz,   \
                             tk_xxxx_xxxzz,  \
                             tk_xxxx_xxyy,   \
                             tk_xxxx_xxyyy,  \
                             tk_xxxx_xxyyz,  \
                             tk_xxxx_xxyz,   \
                             tk_xxxx_xxyzz,  \
                             tk_xxxx_xxzz,   \
                             tk_xxxx_xxzzz,  \
                             tk_xxxx_xyyy,   \
                             tk_xxxx_xyyyy,  \
                             tk_xxxx_xyyyz,  \
                             tk_xxxx_xyyz,   \
                             tk_xxxx_xyyzz,  \
                             tk_xxxx_xyzz,   \
                             tk_xxxx_xyzzz,  \
                             tk_xxxx_xzzz,   \
                             tk_xxxx_xzzzz,  \
                             tk_xxxx_yyyy,   \
                             tk_xxxx_yyyyy,  \
                             tk_xxxx_yyyyz,  \
                             tk_xxxx_yyyz,   \
                             tk_xxxx_yyyzz,  \
                             tk_xxxx_yyzz,   \
                             tk_xxxx_yyzzz,  \
                             tk_xxxx_yzzz,   \
                             tk_xxxx_yzzzz,  \
                             tk_xxxx_zzzz,   \
                             tk_xxxx_zzzzz,  \
                             tk_xxxxx_xxxxx, \
                             tk_xxxxx_xxxxy, \
                             tk_xxxxx_xxxxz, \
                             tk_xxxxx_xxxyy, \
                             tk_xxxxx_xxxyz, \
                             tk_xxxxx_xxxzz, \
                             tk_xxxxx_xxyyy, \
                             tk_xxxxx_xxyyz, \
                             tk_xxxxx_xxyzz, \
                             tk_xxxxx_xxzzz, \
                             tk_xxxxx_xyyyy, \
                             tk_xxxxx_xyyyz, \
                             tk_xxxxx_xyyzz, \
                             tk_xxxxx_xyzzz, \
                             tk_xxxxx_xzzzz, \
                             tk_xxxxx_yyyyy, \
                             tk_xxxxx_yyyyz, \
                             tk_xxxxx_yyyzz, \
                             tk_xxxxx_yyzzz, \
                             tk_xxxxx_yzzzz, \
                             tk_xxxxx_zzzzz, \
                             ts_xxx_xxxxx,   \
                             ts_xxx_xxxxy,   \
                             ts_xxx_xxxxz,   \
                             ts_xxx_xxxyy,   \
                             ts_xxx_xxxyz,   \
                             ts_xxx_xxxzz,   \
                             ts_xxx_xxyyy,   \
                             ts_xxx_xxyyz,   \
                             ts_xxx_xxyzz,   \
                             ts_xxx_xxzzz,   \
                             ts_xxx_xyyyy,   \
                             ts_xxx_xyyyz,   \
                             ts_xxx_xyyzz,   \
                             ts_xxx_xyzzz,   \
                             ts_xxx_xzzzz,   \
                             ts_xxx_yyyyy,   \
                             ts_xxx_yyyyz,   \
                             ts_xxx_yyyzz,   \
                             ts_xxx_yyzzz,   \
                             ts_xxx_yzzzz,   \
                             ts_xxx_zzzzz,   \
                             ts_xxxxx_xxxxx, \
                             ts_xxxxx_xxxxy, \
                             ts_xxxxx_xxxxz, \
                             ts_xxxxx_xxxyy, \
                             ts_xxxxx_xxxyz, \
                             ts_xxxxx_xxxzz, \
                             ts_xxxxx_xxyyy, \
                             ts_xxxxx_xxyyz, \
                             ts_xxxxx_xxyzz, \
                             ts_xxxxx_xxzzz, \
                             ts_xxxxx_xyyyy, \
                             ts_xxxxx_xyyyz, \
                             ts_xxxxx_xyyzz, \
                             ts_xxxxx_xyzzz, \
                             ts_xxxxx_xzzzz, \
                             ts_xxxxx_yyyyy, \
                             ts_xxxxx_yyyyz, \
                             ts_xxxxx_yyyzz, \
                             ts_xxxxx_yyzzz, \
                             ts_xxxxx_yzzzz, \
                             ts_xxxxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_xxxxx[i] = -8.0 * ts_xxx_xxxxx[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxx[i] * fe_0 + 5.0 * tk_xxxx_xxxx[i] * fe_0 +
                            tk_xxxx_xxxxx[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxx[i] * fz_0;

        tk_xxxxx_xxxxy[i] = -8.0 * ts_xxx_xxxxy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxy[i] * fe_0 + 4.0 * tk_xxxx_xxxy[i] * fe_0 +
                            tk_xxxx_xxxxy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxy[i] * fz_0;

        tk_xxxxx_xxxxz[i] = -8.0 * ts_xxx_xxxxz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxz[i] * fe_0 + 4.0 * tk_xxxx_xxxz[i] * fe_0 +
                            tk_xxxx_xxxxz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxz[i] * fz_0;

        tk_xxxxx_xxxyy[i] = -8.0 * ts_xxx_xxxyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxyy[i] * fe_0 + 3.0 * tk_xxxx_xxyy[i] * fe_0 +
                            tk_xxxx_xxxyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxyy[i] * fz_0;

        tk_xxxxx_xxxyz[i] = -8.0 * ts_xxx_xxxyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxyz[i] * fe_0 + 3.0 * tk_xxxx_xxyz[i] * fe_0 +
                            tk_xxxx_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxyz[i] * fz_0;

        tk_xxxxx_xxxzz[i] = -8.0 * ts_xxx_xxxzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxzz[i] * fe_0 + 3.0 * tk_xxxx_xxzz[i] * fe_0 +
                            tk_xxxx_xxxzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxzz[i] * fz_0;

        tk_xxxxx_xxyyy[i] = -8.0 * ts_xxx_xxyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyyy[i] * fe_0 + 2.0 * tk_xxxx_xyyy[i] * fe_0 +
                            tk_xxxx_xxyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyyy[i] * fz_0;

        tk_xxxxx_xxyyz[i] = -8.0 * ts_xxx_xxyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyyz[i] * fe_0 + 2.0 * tk_xxxx_xyyz[i] * fe_0 +
                            tk_xxxx_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyyz[i] * fz_0;

        tk_xxxxx_xxyzz[i] = -8.0 * ts_xxx_xxyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyzz[i] * fe_0 + 2.0 * tk_xxxx_xyzz[i] * fe_0 +
                            tk_xxxx_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyzz[i] * fz_0;

        tk_xxxxx_xxzzz[i] = -8.0 * ts_xxx_xxzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxzzz[i] * fe_0 + 2.0 * tk_xxxx_xzzz[i] * fe_0 +
                            tk_xxxx_xxzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxzzz[i] * fz_0;

        tk_xxxxx_xyyyy[i] = -8.0 * ts_xxx_xyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyyy[i] * fe_0 + tk_xxxx_yyyy[i] * fe_0 +
                            tk_xxxx_xyyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyyy[i] * fz_0;

        tk_xxxxx_xyyyz[i] = -8.0 * ts_xxx_xyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyyz[i] * fe_0 + tk_xxxx_yyyz[i] * fe_0 +
                            tk_xxxx_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyyz[i] * fz_0;

        tk_xxxxx_xyyzz[i] = -8.0 * ts_xxx_xyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyzz[i] * fe_0 + tk_xxxx_yyzz[i] * fe_0 +
                            tk_xxxx_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyzz[i] * fz_0;

        tk_xxxxx_xyzzz[i] = -8.0 * ts_xxx_xyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyzzz[i] * fe_0 + tk_xxxx_yzzz[i] * fe_0 +
                            tk_xxxx_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyzzz[i] * fz_0;

        tk_xxxxx_xzzzz[i] = -8.0 * ts_xxx_xzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xzzzz[i] * fe_0 + tk_xxxx_zzzz[i] * fe_0 +
                            tk_xxxx_xzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xzzzz[i] * fz_0;

        tk_xxxxx_yyyyy[i] =
            -8.0 * ts_xxx_yyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyyy[i] * fe_0 + tk_xxxx_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyyy[i] * fz_0;

        tk_xxxxx_yyyyz[i] =
            -8.0 * ts_xxx_yyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyyz[i] * fe_0 + tk_xxxx_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyyz[i] * fz_0;

        tk_xxxxx_yyyzz[i] =
            -8.0 * ts_xxx_yyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyzz[i] * fe_0 + tk_xxxx_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyzz[i] * fz_0;

        tk_xxxxx_yyzzz[i] =
            -8.0 * ts_xxx_yyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyzzz[i] * fe_0 + tk_xxxx_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyzzz[i] * fz_0;

        tk_xxxxx_yzzzz[i] =
            -8.0 * ts_xxx_yzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yzzzz[i] * fe_0 + tk_xxxx_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yzzzz[i] * fz_0;

        tk_xxxxx_zzzzz[i] =
            -8.0 * ts_xxx_zzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_zzzzz[i] * fe_0 + tk_xxxx_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_zzzzz[i] * fz_0;
    }

    // Set up 21-42 components of targeted buffer : HH

    auto tk_xxxxy_xxxxx = pbuffer.data(idx_kin_hh + 21);

    auto tk_xxxxy_xxxxy = pbuffer.data(idx_kin_hh + 22);

    auto tk_xxxxy_xxxxz = pbuffer.data(idx_kin_hh + 23);

    auto tk_xxxxy_xxxyy = pbuffer.data(idx_kin_hh + 24);

    auto tk_xxxxy_xxxyz = pbuffer.data(idx_kin_hh + 25);

    auto tk_xxxxy_xxxzz = pbuffer.data(idx_kin_hh + 26);

    auto tk_xxxxy_xxyyy = pbuffer.data(idx_kin_hh + 27);

    auto tk_xxxxy_xxyyz = pbuffer.data(idx_kin_hh + 28);

    auto tk_xxxxy_xxyzz = pbuffer.data(idx_kin_hh + 29);

    auto tk_xxxxy_xxzzz = pbuffer.data(idx_kin_hh + 30);

    auto tk_xxxxy_xyyyy = pbuffer.data(idx_kin_hh + 31);

    auto tk_xxxxy_xyyyz = pbuffer.data(idx_kin_hh + 32);

    auto tk_xxxxy_xyyzz = pbuffer.data(idx_kin_hh + 33);

    auto tk_xxxxy_xyzzz = pbuffer.data(idx_kin_hh + 34);

    auto tk_xxxxy_xzzzz = pbuffer.data(idx_kin_hh + 35);

    auto tk_xxxxy_yyyyy = pbuffer.data(idx_kin_hh + 36);

    auto tk_xxxxy_yyyyz = pbuffer.data(idx_kin_hh + 37);

    auto tk_xxxxy_yyyzz = pbuffer.data(idx_kin_hh + 38);

    auto tk_xxxxy_yyzzz = pbuffer.data(idx_kin_hh + 39);

    auto tk_xxxxy_yzzzz = pbuffer.data(idx_kin_hh + 40);

    auto tk_xxxxy_zzzzz = pbuffer.data(idx_kin_hh + 41);

#pragma omp simd aligned(pa_y,               \
                             tk_xxxx_xxxx,   \
                             tk_xxxx_xxxxx,  \
                             tk_xxxx_xxxxy,  \
                             tk_xxxx_xxxxz,  \
                             tk_xxxx_xxxy,   \
                             tk_xxxx_xxxyy,  \
                             tk_xxxx_xxxyz,  \
                             tk_xxxx_xxxz,   \
                             tk_xxxx_xxxzz,  \
                             tk_xxxx_xxyy,   \
                             tk_xxxx_xxyyy,  \
                             tk_xxxx_xxyyz,  \
                             tk_xxxx_xxyz,   \
                             tk_xxxx_xxyzz,  \
                             tk_xxxx_xxzz,   \
                             tk_xxxx_xxzzz,  \
                             tk_xxxx_xyyy,   \
                             tk_xxxx_xyyyy,  \
                             tk_xxxx_xyyyz,  \
                             tk_xxxx_xyyz,   \
                             tk_xxxx_xyyzz,  \
                             tk_xxxx_xyzz,   \
                             tk_xxxx_xyzzz,  \
                             tk_xxxx_xzzz,   \
                             tk_xxxx_xzzzz,  \
                             tk_xxxx_yyyy,   \
                             tk_xxxx_yyyyy,  \
                             tk_xxxx_yyyyz,  \
                             tk_xxxx_yyyz,   \
                             tk_xxxx_yyyzz,  \
                             tk_xxxx_yyzz,   \
                             tk_xxxx_yyzzz,  \
                             tk_xxxx_yzzz,   \
                             tk_xxxx_yzzzz,  \
                             tk_xxxx_zzzz,   \
                             tk_xxxx_zzzzz,  \
                             tk_xxxxy_xxxxx, \
                             tk_xxxxy_xxxxy, \
                             tk_xxxxy_xxxxz, \
                             tk_xxxxy_xxxyy, \
                             tk_xxxxy_xxxyz, \
                             tk_xxxxy_xxxzz, \
                             tk_xxxxy_xxyyy, \
                             tk_xxxxy_xxyyz, \
                             tk_xxxxy_xxyzz, \
                             tk_xxxxy_xxzzz, \
                             tk_xxxxy_xyyyy, \
                             tk_xxxxy_xyyyz, \
                             tk_xxxxy_xyyzz, \
                             tk_xxxxy_xyzzz, \
                             tk_xxxxy_xzzzz, \
                             tk_xxxxy_yyyyy, \
                             tk_xxxxy_yyyyz, \
                             tk_xxxxy_yyyzz, \
                             tk_xxxxy_yyzzz, \
                             tk_xxxxy_yzzzz, \
                             tk_xxxxy_zzzzz, \
                             ts_xxxxy_xxxxx, \
                             ts_xxxxy_xxxxy, \
                             ts_xxxxy_xxxxz, \
                             ts_xxxxy_xxxyy, \
                             ts_xxxxy_xxxyz, \
                             ts_xxxxy_xxxzz, \
                             ts_xxxxy_xxyyy, \
                             ts_xxxxy_xxyyz, \
                             ts_xxxxy_xxyzz, \
                             ts_xxxxy_xxzzz, \
                             ts_xxxxy_xyyyy, \
                             ts_xxxxy_xyyyz, \
                             ts_xxxxy_xyyzz, \
                             ts_xxxxy_xyzzz, \
                             ts_xxxxy_xzzzz, \
                             ts_xxxxy_yyyyy, \
                             ts_xxxxy_yyyyz, \
                             ts_xxxxy_yyyzz, \
                             ts_xxxxy_yyzzz, \
                             ts_xxxxy_yzzzz, \
                             ts_xxxxy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxy_xxxxx[i] = tk_xxxx_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxx[i] * fz_0;

        tk_xxxxy_xxxxy[i] = tk_xxxx_xxxx[i] * fe_0 + tk_xxxx_xxxxy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxy[i] * fz_0;

        tk_xxxxy_xxxxz[i] = tk_xxxx_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxz[i] * fz_0;

        tk_xxxxy_xxxyy[i] = 2.0 * tk_xxxx_xxxy[i] * fe_0 + tk_xxxx_xxxyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxyy[i] * fz_0;

        tk_xxxxy_xxxyz[i] = tk_xxxx_xxxz[i] * fe_0 + tk_xxxx_xxxyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxyz[i] * fz_0;

        tk_xxxxy_xxxzz[i] = tk_xxxx_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxzz[i] * fz_0;

        tk_xxxxy_xxyyy[i] = 3.0 * tk_xxxx_xxyy[i] * fe_0 + tk_xxxx_xxyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyyy[i] * fz_0;

        tk_xxxxy_xxyyz[i] = 2.0 * tk_xxxx_xxyz[i] * fe_0 + tk_xxxx_xxyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyyz[i] * fz_0;

        tk_xxxxy_xxyzz[i] = tk_xxxx_xxzz[i] * fe_0 + tk_xxxx_xxyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyzz[i] * fz_0;

        tk_xxxxy_xxzzz[i] = tk_xxxx_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxzzz[i] * fz_0;

        tk_xxxxy_xyyyy[i] = 4.0 * tk_xxxx_xyyy[i] * fe_0 + tk_xxxx_xyyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyyy[i] * fz_0;

        tk_xxxxy_xyyyz[i] = 3.0 * tk_xxxx_xyyz[i] * fe_0 + tk_xxxx_xyyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyyz[i] * fz_0;

        tk_xxxxy_xyyzz[i] = 2.0 * tk_xxxx_xyzz[i] * fe_0 + tk_xxxx_xyyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyzz[i] * fz_0;

        tk_xxxxy_xyzzz[i] = tk_xxxx_xzzz[i] * fe_0 + tk_xxxx_xyzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyzzz[i] * fz_0;

        tk_xxxxy_xzzzz[i] = tk_xxxx_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xzzzz[i] * fz_0;

        tk_xxxxy_yyyyy[i] = 5.0 * tk_xxxx_yyyy[i] * fe_0 + tk_xxxx_yyyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyyy[i] * fz_0;

        tk_xxxxy_yyyyz[i] = 4.0 * tk_xxxx_yyyz[i] * fe_0 + tk_xxxx_yyyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyyz[i] * fz_0;

        tk_xxxxy_yyyzz[i] = 3.0 * tk_xxxx_yyzz[i] * fe_0 + tk_xxxx_yyyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyzz[i] * fz_0;

        tk_xxxxy_yyzzz[i] = 2.0 * tk_xxxx_yzzz[i] * fe_0 + tk_xxxx_yyzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyzzz[i] * fz_0;

        tk_xxxxy_yzzzz[i] = tk_xxxx_zzzz[i] * fe_0 + tk_xxxx_yzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yzzzz[i] * fz_0;

        tk_xxxxy_zzzzz[i] = tk_xxxx_zzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_zzzzz[i] * fz_0;
    }

    // Set up 42-63 components of targeted buffer : HH

    auto tk_xxxxz_xxxxx = pbuffer.data(idx_kin_hh + 42);

    auto tk_xxxxz_xxxxy = pbuffer.data(idx_kin_hh + 43);

    auto tk_xxxxz_xxxxz = pbuffer.data(idx_kin_hh + 44);

    auto tk_xxxxz_xxxyy = pbuffer.data(idx_kin_hh + 45);

    auto tk_xxxxz_xxxyz = pbuffer.data(idx_kin_hh + 46);

    auto tk_xxxxz_xxxzz = pbuffer.data(idx_kin_hh + 47);

    auto tk_xxxxz_xxyyy = pbuffer.data(idx_kin_hh + 48);

    auto tk_xxxxz_xxyyz = pbuffer.data(idx_kin_hh + 49);

    auto tk_xxxxz_xxyzz = pbuffer.data(idx_kin_hh + 50);

    auto tk_xxxxz_xxzzz = pbuffer.data(idx_kin_hh + 51);

    auto tk_xxxxz_xyyyy = pbuffer.data(idx_kin_hh + 52);

    auto tk_xxxxz_xyyyz = pbuffer.data(idx_kin_hh + 53);

    auto tk_xxxxz_xyyzz = pbuffer.data(idx_kin_hh + 54);

    auto tk_xxxxz_xyzzz = pbuffer.data(idx_kin_hh + 55);

    auto tk_xxxxz_xzzzz = pbuffer.data(idx_kin_hh + 56);

    auto tk_xxxxz_yyyyy = pbuffer.data(idx_kin_hh + 57);

    auto tk_xxxxz_yyyyz = pbuffer.data(idx_kin_hh + 58);

    auto tk_xxxxz_yyyzz = pbuffer.data(idx_kin_hh + 59);

    auto tk_xxxxz_yyzzz = pbuffer.data(idx_kin_hh + 60);

    auto tk_xxxxz_yzzzz = pbuffer.data(idx_kin_hh + 61);

    auto tk_xxxxz_zzzzz = pbuffer.data(idx_kin_hh + 62);

#pragma omp simd aligned(pa_z,               \
                             tk_xxxx_xxxx,   \
                             tk_xxxx_xxxxx,  \
                             tk_xxxx_xxxxy,  \
                             tk_xxxx_xxxxz,  \
                             tk_xxxx_xxxy,   \
                             tk_xxxx_xxxyy,  \
                             tk_xxxx_xxxyz,  \
                             tk_xxxx_xxxz,   \
                             tk_xxxx_xxxzz,  \
                             tk_xxxx_xxyy,   \
                             tk_xxxx_xxyyy,  \
                             tk_xxxx_xxyyz,  \
                             tk_xxxx_xxyz,   \
                             tk_xxxx_xxyzz,  \
                             tk_xxxx_xxzz,   \
                             tk_xxxx_xxzzz,  \
                             tk_xxxx_xyyy,   \
                             tk_xxxx_xyyyy,  \
                             tk_xxxx_xyyyz,  \
                             tk_xxxx_xyyz,   \
                             tk_xxxx_xyyzz,  \
                             tk_xxxx_xyzz,   \
                             tk_xxxx_xyzzz,  \
                             tk_xxxx_xzzz,   \
                             tk_xxxx_xzzzz,  \
                             tk_xxxx_yyyy,   \
                             tk_xxxx_yyyyy,  \
                             tk_xxxx_yyyyz,  \
                             tk_xxxx_yyyz,   \
                             tk_xxxx_yyyzz,  \
                             tk_xxxx_yyzz,   \
                             tk_xxxx_yyzzz,  \
                             tk_xxxx_yzzz,   \
                             tk_xxxx_yzzzz,  \
                             tk_xxxx_zzzz,   \
                             tk_xxxx_zzzzz,  \
                             tk_xxxxz_xxxxx, \
                             tk_xxxxz_xxxxy, \
                             tk_xxxxz_xxxxz, \
                             tk_xxxxz_xxxyy, \
                             tk_xxxxz_xxxyz, \
                             tk_xxxxz_xxxzz, \
                             tk_xxxxz_xxyyy, \
                             tk_xxxxz_xxyyz, \
                             tk_xxxxz_xxyzz, \
                             tk_xxxxz_xxzzz, \
                             tk_xxxxz_xyyyy, \
                             tk_xxxxz_xyyyz, \
                             tk_xxxxz_xyyzz, \
                             tk_xxxxz_xyzzz, \
                             tk_xxxxz_xzzzz, \
                             tk_xxxxz_yyyyy, \
                             tk_xxxxz_yyyyz, \
                             tk_xxxxz_yyyzz, \
                             tk_xxxxz_yyzzz, \
                             tk_xxxxz_yzzzz, \
                             tk_xxxxz_zzzzz, \
                             ts_xxxxz_xxxxx, \
                             ts_xxxxz_xxxxy, \
                             ts_xxxxz_xxxxz, \
                             ts_xxxxz_xxxyy, \
                             ts_xxxxz_xxxyz, \
                             ts_xxxxz_xxxzz, \
                             ts_xxxxz_xxyyy, \
                             ts_xxxxz_xxyyz, \
                             ts_xxxxz_xxyzz, \
                             ts_xxxxz_xxzzz, \
                             ts_xxxxz_xyyyy, \
                             ts_xxxxz_xyyyz, \
                             ts_xxxxz_xyyzz, \
                             ts_xxxxz_xyzzz, \
                             ts_xxxxz_xzzzz, \
                             ts_xxxxz_yyyyy, \
                             ts_xxxxz_yyyyz, \
                             ts_xxxxz_yyyzz, \
                             ts_xxxxz_yyzzz, \
                             ts_xxxxz_yzzzz, \
                             ts_xxxxz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxz_xxxxx[i] = tk_xxxx_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxx[i] * fz_0;

        tk_xxxxz_xxxxy[i] = tk_xxxx_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxy[i] * fz_0;

        tk_xxxxz_xxxxz[i] = tk_xxxx_xxxx[i] * fe_0 + tk_xxxx_xxxxz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxz[i] * fz_0;

        tk_xxxxz_xxxyy[i] = tk_xxxx_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxyy[i] * fz_0;

        tk_xxxxz_xxxyz[i] = tk_xxxx_xxxy[i] * fe_0 + tk_xxxx_xxxyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxyz[i] * fz_0;

        tk_xxxxz_xxxzz[i] = 2.0 * tk_xxxx_xxxz[i] * fe_0 + tk_xxxx_xxxzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxzz[i] * fz_0;

        tk_xxxxz_xxyyy[i] = tk_xxxx_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyyy[i] * fz_0;

        tk_xxxxz_xxyyz[i] = tk_xxxx_xxyy[i] * fe_0 + tk_xxxx_xxyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyyz[i] * fz_0;

        tk_xxxxz_xxyzz[i] = 2.0 * tk_xxxx_xxyz[i] * fe_0 + tk_xxxx_xxyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyzz[i] * fz_0;

        tk_xxxxz_xxzzz[i] = 3.0 * tk_xxxx_xxzz[i] * fe_0 + tk_xxxx_xxzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxzzz[i] * fz_0;

        tk_xxxxz_xyyyy[i] = tk_xxxx_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyyy[i] * fz_0;

        tk_xxxxz_xyyyz[i] = tk_xxxx_xyyy[i] * fe_0 + tk_xxxx_xyyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyyz[i] * fz_0;

        tk_xxxxz_xyyzz[i] = 2.0 * tk_xxxx_xyyz[i] * fe_0 + tk_xxxx_xyyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyzz[i] * fz_0;

        tk_xxxxz_xyzzz[i] = 3.0 * tk_xxxx_xyzz[i] * fe_0 + tk_xxxx_xyzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyzzz[i] * fz_0;

        tk_xxxxz_xzzzz[i] = 4.0 * tk_xxxx_xzzz[i] * fe_0 + tk_xxxx_xzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xzzzz[i] * fz_0;

        tk_xxxxz_yyyyy[i] = tk_xxxx_yyyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyyy[i] * fz_0;

        tk_xxxxz_yyyyz[i] = tk_xxxx_yyyy[i] * fe_0 + tk_xxxx_yyyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyyz[i] * fz_0;

        tk_xxxxz_yyyzz[i] = 2.0 * tk_xxxx_yyyz[i] * fe_0 + tk_xxxx_yyyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyzz[i] * fz_0;

        tk_xxxxz_yyzzz[i] = 3.0 * tk_xxxx_yyzz[i] * fe_0 + tk_xxxx_yyzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyzzz[i] * fz_0;

        tk_xxxxz_yzzzz[i] = 4.0 * tk_xxxx_yzzz[i] * fe_0 + tk_xxxx_yzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yzzzz[i] * fz_0;

        tk_xxxxz_zzzzz[i] = 5.0 * tk_xxxx_zzzz[i] * fe_0 + tk_xxxx_zzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_zzzzz[i] * fz_0;
    }

    // Set up 63-84 components of targeted buffer : HH

    auto tk_xxxyy_xxxxx = pbuffer.data(idx_kin_hh + 63);

    auto tk_xxxyy_xxxxy = pbuffer.data(idx_kin_hh + 64);

    auto tk_xxxyy_xxxxz = pbuffer.data(idx_kin_hh + 65);

    auto tk_xxxyy_xxxyy = pbuffer.data(idx_kin_hh + 66);

    auto tk_xxxyy_xxxyz = pbuffer.data(idx_kin_hh + 67);

    auto tk_xxxyy_xxxzz = pbuffer.data(idx_kin_hh + 68);

    auto tk_xxxyy_xxyyy = pbuffer.data(idx_kin_hh + 69);

    auto tk_xxxyy_xxyyz = pbuffer.data(idx_kin_hh + 70);

    auto tk_xxxyy_xxyzz = pbuffer.data(idx_kin_hh + 71);

    auto tk_xxxyy_xxzzz = pbuffer.data(idx_kin_hh + 72);

    auto tk_xxxyy_xyyyy = pbuffer.data(idx_kin_hh + 73);

    auto tk_xxxyy_xyyyz = pbuffer.data(idx_kin_hh + 74);

    auto tk_xxxyy_xyyzz = pbuffer.data(idx_kin_hh + 75);

    auto tk_xxxyy_xyzzz = pbuffer.data(idx_kin_hh + 76);

    auto tk_xxxyy_xzzzz = pbuffer.data(idx_kin_hh + 77);

    auto tk_xxxyy_yyyyy = pbuffer.data(idx_kin_hh + 78);

    auto tk_xxxyy_yyyyz = pbuffer.data(idx_kin_hh + 79);

    auto tk_xxxyy_yyyzz = pbuffer.data(idx_kin_hh + 80);

    auto tk_xxxyy_yyzzz = pbuffer.data(idx_kin_hh + 81);

    auto tk_xxxyy_yzzzz = pbuffer.data(idx_kin_hh + 82);

    auto tk_xxxyy_zzzzz = pbuffer.data(idx_kin_hh + 83);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xxx_xxxxx,   \
                             tk_xxx_xxxxz,   \
                             tk_xxx_xxxzz,   \
                             tk_xxx_xxzzz,   \
                             tk_xxx_xzzzz,   \
                             tk_xxxy_xxxxx,  \
                             tk_xxxy_xxxxz,  \
                             tk_xxxy_xxxzz,  \
                             tk_xxxy_xxzzz,  \
                             tk_xxxy_xzzzz,  \
                             tk_xxxyy_xxxxx, \
                             tk_xxxyy_xxxxy, \
                             tk_xxxyy_xxxxz, \
                             tk_xxxyy_xxxyy, \
                             tk_xxxyy_xxxyz, \
                             tk_xxxyy_xxxzz, \
                             tk_xxxyy_xxyyy, \
                             tk_xxxyy_xxyyz, \
                             tk_xxxyy_xxyzz, \
                             tk_xxxyy_xxzzz, \
                             tk_xxxyy_xyyyy, \
                             tk_xxxyy_xyyyz, \
                             tk_xxxyy_xyyzz, \
                             tk_xxxyy_xyzzz, \
                             tk_xxxyy_xzzzz, \
                             tk_xxxyy_yyyyy, \
                             tk_xxxyy_yyyyz, \
                             tk_xxxyy_yyyzz, \
                             tk_xxxyy_yyzzz, \
                             tk_xxxyy_yzzzz, \
                             tk_xxxyy_zzzzz, \
                             tk_xxyy_xxxxy,  \
                             tk_xxyy_xxxy,   \
                             tk_xxyy_xxxyy,  \
                             tk_xxyy_xxxyz,  \
                             tk_xxyy_xxyy,   \
                             tk_xxyy_xxyyy,  \
                             tk_xxyy_xxyyz,  \
                             tk_xxyy_xxyz,   \
                             tk_xxyy_xxyzz,  \
                             tk_xxyy_xyyy,   \
                             tk_xxyy_xyyyy,  \
                             tk_xxyy_xyyyz,  \
                             tk_xxyy_xyyz,   \
                             tk_xxyy_xyyzz,  \
                             tk_xxyy_xyzz,   \
                             tk_xxyy_xyzzz,  \
                             tk_xxyy_yyyy,   \
                             tk_xxyy_yyyyy,  \
                             tk_xxyy_yyyyz,  \
                             tk_xxyy_yyyz,   \
                             tk_xxyy_yyyzz,  \
                             tk_xxyy_yyzz,   \
                             tk_xxyy_yyzzz,  \
                             tk_xxyy_yzzz,   \
                             tk_xxyy_yzzzz,  \
                             tk_xxyy_zzzzz,  \
                             tk_xyy_xxxxy,   \
                             tk_xyy_xxxyy,   \
                             tk_xyy_xxxyz,   \
                             tk_xyy_xxyyy,   \
                             tk_xyy_xxyyz,   \
                             tk_xyy_xxyzz,   \
                             tk_xyy_xyyyy,   \
                             tk_xyy_xyyyz,   \
                             tk_xyy_xyyzz,   \
                             tk_xyy_xyzzz,   \
                             tk_xyy_yyyyy,   \
                             tk_xyy_yyyyz,   \
                             tk_xyy_yyyzz,   \
                             tk_xyy_yyzzz,   \
                             tk_xyy_yzzzz,   \
                             tk_xyy_zzzzz,   \
                             ts_xxx_xxxxx,   \
                             ts_xxx_xxxxz,   \
                             ts_xxx_xxxzz,   \
                             ts_xxx_xxzzz,   \
                             ts_xxx_xzzzz,   \
                             ts_xxxyy_xxxxx, \
                             ts_xxxyy_xxxxy, \
                             ts_xxxyy_xxxxz, \
                             ts_xxxyy_xxxyy, \
                             ts_xxxyy_xxxyz, \
                             ts_xxxyy_xxxzz, \
                             ts_xxxyy_xxyyy, \
                             ts_xxxyy_xxyyz, \
                             ts_xxxyy_xxyzz, \
                             ts_xxxyy_xxzzz, \
                             ts_xxxyy_xyyyy, \
                             ts_xxxyy_xyyyz, \
                             ts_xxxyy_xyyzz, \
                             ts_xxxyy_xyzzz, \
                             ts_xxxyy_xzzzz, \
                             ts_xxxyy_yyyyy, \
                             ts_xxxyy_yyyyz, \
                             ts_xxxyy_yyyzz, \
                             ts_xxxyy_yyzzz, \
                             ts_xxxyy_yzzzz, \
                             ts_xxxyy_zzzzz, \
                             ts_xyy_xxxxy,   \
                             ts_xyy_xxxyy,   \
                             ts_xyy_xxxyz,   \
                             ts_xyy_xxyyy,   \
                             ts_xyy_xxyyz,   \
                             ts_xyy_xxyzz,   \
                             ts_xyy_xyyyy,   \
                             ts_xyy_xyyyz,   \
                             ts_xyy_xyyzz,   \
                             ts_xyy_xyzzz,   \
                             ts_xyy_yyyyy,   \
                             ts_xyy_yyyyz,   \
                             ts_xyy_yyyzz,   \
                             ts_xyy_yyzzz,   \
                             ts_xyy_yzzzz,   \
                             ts_xyy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyy_xxxxx[i] =
            -2.0 * ts_xxx_xxxxx[i] * fbe_0 * fz_0 + tk_xxx_xxxxx[i] * fe_0 + tk_xxxy_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxxx[i] * fz_0;

        tk_xxxyy_xxxxy[i] = -4.0 * ts_xyy_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxxy[i] * fe_0 + 4.0 * tk_xxyy_xxxy[i] * fe_0 +
                            tk_xxyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxxy[i] * fz_0;

        tk_xxxyy_xxxxz[i] =
            -2.0 * ts_xxx_xxxxz[i] * fbe_0 * fz_0 + tk_xxx_xxxxz[i] * fe_0 + tk_xxxy_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxxz[i] * fz_0;

        tk_xxxyy_xxxyy[i] = -4.0 * ts_xyy_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxyy[i] * fe_0 + 3.0 * tk_xxyy_xxyy[i] * fe_0 +
                            tk_xxyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxyy[i] * fz_0;

        tk_xxxyy_xxxyz[i] = -4.0 * ts_xyy_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxyz[i] * fe_0 + 3.0 * tk_xxyy_xxyz[i] * fe_0 +
                            tk_xxyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxyz[i] * fz_0;

        tk_xxxyy_xxxzz[i] =
            -2.0 * ts_xxx_xxxzz[i] * fbe_0 * fz_0 + tk_xxx_xxxzz[i] * fe_0 + tk_xxxy_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxzz[i] * fz_0;

        tk_xxxyy_xxyyy[i] = -4.0 * ts_xyy_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyyy[i] * fe_0 + 2.0 * tk_xxyy_xyyy[i] * fe_0 +
                            tk_xxyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyyy[i] * fz_0;

        tk_xxxyy_xxyyz[i] = -4.0 * ts_xyy_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyyz[i] * fe_0 + 2.0 * tk_xxyy_xyyz[i] * fe_0 +
                            tk_xxyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyyz[i] * fz_0;

        tk_xxxyy_xxyzz[i] = -4.0 * ts_xyy_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyzz[i] * fe_0 + 2.0 * tk_xxyy_xyzz[i] * fe_0 +
                            tk_xxyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyzz[i] * fz_0;

        tk_xxxyy_xxzzz[i] =
            -2.0 * ts_xxx_xxzzz[i] * fbe_0 * fz_0 + tk_xxx_xxzzz[i] * fe_0 + tk_xxxy_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxzzz[i] * fz_0;

        tk_xxxyy_xyyyy[i] = -4.0 * ts_xyy_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyyy[i] * fe_0 + tk_xxyy_yyyy[i] * fe_0 +
                            tk_xxyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyyy[i] * fz_0;

        tk_xxxyy_xyyyz[i] = -4.0 * ts_xyy_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyyz[i] * fe_0 + tk_xxyy_yyyz[i] * fe_0 +
                            tk_xxyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyyz[i] * fz_0;

        tk_xxxyy_xyyzz[i] = -4.0 * ts_xyy_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyzz[i] * fe_0 + tk_xxyy_yyzz[i] * fe_0 +
                            tk_xxyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyzz[i] * fz_0;

        tk_xxxyy_xyzzz[i] = -4.0 * ts_xyy_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyzzz[i] * fe_0 + tk_xxyy_yzzz[i] * fe_0 +
                            tk_xxyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyzzz[i] * fz_0;

        tk_xxxyy_xzzzz[i] =
            -2.0 * ts_xxx_xzzzz[i] * fbe_0 * fz_0 + tk_xxx_xzzzz[i] * fe_0 + tk_xxxy_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xzzzz[i] * fz_0;

        tk_xxxyy_yyyyy[i] =
            -4.0 * ts_xyy_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyyy[i] * fe_0 + tk_xxyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyyy[i] * fz_0;

        tk_xxxyy_yyyyz[i] =
            -4.0 * ts_xyy_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyyz[i] * fe_0 + tk_xxyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyyz[i] * fz_0;

        tk_xxxyy_yyyzz[i] =
            -4.0 * ts_xyy_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyzz[i] * fe_0 + tk_xxyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyzz[i] * fz_0;

        tk_xxxyy_yyzzz[i] =
            -4.0 * ts_xyy_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyzzz[i] * fe_0 + tk_xxyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyzzz[i] * fz_0;

        tk_xxxyy_yzzzz[i] =
            -4.0 * ts_xyy_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yzzzz[i] * fe_0 + tk_xxyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yzzzz[i] * fz_0;

        tk_xxxyy_zzzzz[i] =
            -4.0 * ts_xyy_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_zzzzz[i] * fe_0 + tk_xxyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_zzzzz[i] * fz_0;
    }

    // Set up 84-105 components of targeted buffer : HH

    auto tk_xxxyz_xxxxx = pbuffer.data(idx_kin_hh + 84);

    auto tk_xxxyz_xxxxy = pbuffer.data(idx_kin_hh + 85);

    auto tk_xxxyz_xxxxz = pbuffer.data(idx_kin_hh + 86);

    auto tk_xxxyz_xxxyy = pbuffer.data(idx_kin_hh + 87);

    auto tk_xxxyz_xxxyz = pbuffer.data(idx_kin_hh + 88);

    auto tk_xxxyz_xxxzz = pbuffer.data(idx_kin_hh + 89);

    auto tk_xxxyz_xxyyy = pbuffer.data(idx_kin_hh + 90);

    auto tk_xxxyz_xxyyz = pbuffer.data(idx_kin_hh + 91);

    auto tk_xxxyz_xxyzz = pbuffer.data(idx_kin_hh + 92);

    auto tk_xxxyz_xxzzz = pbuffer.data(idx_kin_hh + 93);

    auto tk_xxxyz_xyyyy = pbuffer.data(idx_kin_hh + 94);

    auto tk_xxxyz_xyyyz = pbuffer.data(idx_kin_hh + 95);

    auto tk_xxxyz_xyyzz = pbuffer.data(idx_kin_hh + 96);

    auto tk_xxxyz_xyzzz = pbuffer.data(idx_kin_hh + 97);

    auto tk_xxxyz_xzzzz = pbuffer.data(idx_kin_hh + 98);

    auto tk_xxxyz_yyyyy = pbuffer.data(idx_kin_hh + 99);

    auto tk_xxxyz_yyyyz = pbuffer.data(idx_kin_hh + 100);

    auto tk_xxxyz_yyyzz = pbuffer.data(idx_kin_hh + 101);

    auto tk_xxxyz_yyzzz = pbuffer.data(idx_kin_hh + 102);

    auto tk_xxxyz_yzzzz = pbuffer.data(idx_kin_hh + 103);

    auto tk_xxxyz_zzzzz = pbuffer.data(idx_kin_hh + 104);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_xxxy_xxxxy,  \
                             tk_xxxy_xxxyy,  \
                             tk_xxxy_xxyyy,  \
                             tk_xxxy_xyyyy,  \
                             tk_xxxy_yyyyy,  \
                             tk_xxxyz_xxxxx, \
                             tk_xxxyz_xxxxy, \
                             tk_xxxyz_xxxxz, \
                             tk_xxxyz_xxxyy, \
                             tk_xxxyz_xxxyz, \
                             tk_xxxyz_xxxzz, \
                             tk_xxxyz_xxyyy, \
                             tk_xxxyz_xxyyz, \
                             tk_xxxyz_xxyzz, \
                             tk_xxxyz_xxzzz, \
                             tk_xxxyz_xyyyy, \
                             tk_xxxyz_xyyyz, \
                             tk_xxxyz_xyyzz, \
                             tk_xxxyz_xyzzz, \
                             tk_xxxyz_xzzzz, \
                             tk_xxxyz_yyyyy, \
                             tk_xxxyz_yyyyz, \
                             tk_xxxyz_yyyzz, \
                             tk_xxxyz_yyzzz, \
                             tk_xxxyz_yzzzz, \
                             tk_xxxyz_zzzzz, \
                             tk_xxxz_xxxxx,  \
                             tk_xxxz_xxxxz,  \
                             tk_xxxz_xxxyz,  \
                             tk_xxxz_xxxz,   \
                             tk_xxxz_xxxzz,  \
                             tk_xxxz_xxyyz,  \
                             tk_xxxz_xxyz,   \
                             tk_xxxz_xxyzz,  \
                             tk_xxxz_xxzz,   \
                             tk_xxxz_xxzzz,  \
                             tk_xxxz_xyyyz,  \
                             tk_xxxz_xyyz,   \
                             tk_xxxz_xyyzz,  \
                             tk_xxxz_xyzz,   \
                             tk_xxxz_xyzzz,  \
                             tk_xxxz_xzzz,   \
                             tk_xxxz_xzzzz,  \
                             tk_xxxz_yyyyz,  \
                             tk_xxxz_yyyz,   \
                             tk_xxxz_yyyzz,  \
                             tk_xxxz_yyzz,   \
                             tk_xxxz_yyzzz,  \
                             tk_xxxz_yzzz,   \
                             tk_xxxz_yzzzz,  \
                             tk_xxxz_zzzz,   \
                             tk_xxxz_zzzzz,  \
                             ts_xxxyz_xxxxx, \
                             ts_xxxyz_xxxxy, \
                             ts_xxxyz_xxxxz, \
                             ts_xxxyz_xxxyy, \
                             ts_xxxyz_xxxyz, \
                             ts_xxxyz_xxxzz, \
                             ts_xxxyz_xxyyy, \
                             ts_xxxyz_xxyyz, \
                             ts_xxxyz_xxyzz, \
                             ts_xxxyz_xxzzz, \
                             ts_xxxyz_xyyyy, \
                             ts_xxxyz_xyyyz, \
                             ts_xxxyz_xyyzz, \
                             ts_xxxyz_xyzzz, \
                             ts_xxxyz_xzzzz, \
                             ts_xxxyz_yyyyy, \
                             ts_xxxyz_yyyyz, \
                             ts_xxxyz_yyyzz, \
                             ts_xxxyz_yyzzz, \
                             ts_xxxyz_yzzzz, \
                             ts_xxxyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyz_xxxxx[i] = tk_xxxz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxxx[i] * fz_0;

        tk_xxxyz_xxxxy[i] = tk_xxxy_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxxxy[i] * fz_0;

        tk_xxxyz_xxxxz[i] = tk_xxxz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxxz[i] * fz_0;

        tk_xxxyz_xxxyy[i] = tk_xxxy_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxxyy[i] * fz_0;

        tk_xxxyz_xxxyz[i] = tk_xxxz_xxxz[i] * fe_0 + tk_xxxz_xxxyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxyz[i] * fz_0;

        tk_xxxyz_xxxzz[i] = tk_xxxz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxzz[i] * fz_0;

        tk_xxxyz_xxyyy[i] = tk_xxxy_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxyyy[i] * fz_0;

        tk_xxxyz_xxyyz[i] = 2.0 * tk_xxxz_xxyz[i] * fe_0 + tk_xxxz_xxyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxyyz[i] * fz_0;

        tk_xxxyz_xxyzz[i] = tk_xxxz_xxzz[i] * fe_0 + tk_xxxz_xxyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxyzz[i] * fz_0;

        tk_xxxyz_xxzzz[i] = tk_xxxz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxzzz[i] * fz_0;

        tk_xxxyz_xyyyy[i] = tk_xxxy_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xyyyy[i] * fz_0;

        tk_xxxyz_xyyyz[i] = 3.0 * tk_xxxz_xyyz[i] * fe_0 + tk_xxxz_xyyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyyyz[i] * fz_0;

        tk_xxxyz_xyyzz[i] = 2.0 * tk_xxxz_xyzz[i] * fe_0 + tk_xxxz_xyyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyyzz[i] * fz_0;

        tk_xxxyz_xyzzz[i] = tk_xxxz_xzzz[i] * fe_0 + tk_xxxz_xyzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyzzz[i] * fz_0;

        tk_xxxyz_xzzzz[i] = tk_xxxz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xzzzz[i] * fz_0;

        tk_xxxyz_yyyyy[i] = tk_xxxy_yyyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_yyyyy[i] * fz_0;

        tk_xxxyz_yyyyz[i] = 4.0 * tk_xxxz_yyyz[i] * fe_0 + tk_xxxz_yyyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyyyz[i] * fz_0;

        tk_xxxyz_yyyzz[i] = 3.0 * tk_xxxz_yyzz[i] * fe_0 + tk_xxxz_yyyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyyzz[i] * fz_0;

        tk_xxxyz_yyzzz[i] = 2.0 * tk_xxxz_yzzz[i] * fe_0 + tk_xxxz_yyzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyzzz[i] * fz_0;

        tk_xxxyz_yzzzz[i] = tk_xxxz_zzzz[i] * fe_0 + tk_xxxz_yzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yzzzz[i] * fz_0;

        tk_xxxyz_zzzzz[i] = tk_xxxz_zzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_zzzzz[i] * fz_0;
    }

    // Set up 105-126 components of targeted buffer : HH

    auto tk_xxxzz_xxxxx = pbuffer.data(idx_kin_hh + 105);

    auto tk_xxxzz_xxxxy = pbuffer.data(idx_kin_hh + 106);

    auto tk_xxxzz_xxxxz = pbuffer.data(idx_kin_hh + 107);

    auto tk_xxxzz_xxxyy = pbuffer.data(idx_kin_hh + 108);

    auto tk_xxxzz_xxxyz = pbuffer.data(idx_kin_hh + 109);

    auto tk_xxxzz_xxxzz = pbuffer.data(idx_kin_hh + 110);

    auto tk_xxxzz_xxyyy = pbuffer.data(idx_kin_hh + 111);

    auto tk_xxxzz_xxyyz = pbuffer.data(idx_kin_hh + 112);

    auto tk_xxxzz_xxyzz = pbuffer.data(idx_kin_hh + 113);

    auto tk_xxxzz_xxzzz = pbuffer.data(idx_kin_hh + 114);

    auto tk_xxxzz_xyyyy = pbuffer.data(idx_kin_hh + 115);

    auto tk_xxxzz_xyyyz = pbuffer.data(idx_kin_hh + 116);

    auto tk_xxxzz_xyyzz = pbuffer.data(idx_kin_hh + 117);

    auto tk_xxxzz_xyzzz = pbuffer.data(idx_kin_hh + 118);

    auto tk_xxxzz_xzzzz = pbuffer.data(idx_kin_hh + 119);

    auto tk_xxxzz_yyyyy = pbuffer.data(idx_kin_hh + 120);

    auto tk_xxxzz_yyyyz = pbuffer.data(idx_kin_hh + 121);

    auto tk_xxxzz_yyyzz = pbuffer.data(idx_kin_hh + 122);

    auto tk_xxxzz_yyzzz = pbuffer.data(idx_kin_hh + 123);

    auto tk_xxxzz_yzzzz = pbuffer.data(idx_kin_hh + 124);

    auto tk_xxxzz_zzzzz = pbuffer.data(idx_kin_hh + 125);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xxx_xxxxx,   \
                             tk_xxx_xxxxy,   \
                             tk_xxx_xxxyy,   \
                             tk_xxx_xxyyy,   \
                             tk_xxx_xyyyy,   \
                             tk_xxxz_xxxxx,  \
                             tk_xxxz_xxxxy,  \
                             tk_xxxz_xxxyy,  \
                             tk_xxxz_xxyyy,  \
                             tk_xxxz_xyyyy,  \
                             tk_xxxzz_xxxxx, \
                             tk_xxxzz_xxxxy, \
                             tk_xxxzz_xxxxz, \
                             tk_xxxzz_xxxyy, \
                             tk_xxxzz_xxxyz, \
                             tk_xxxzz_xxxzz, \
                             tk_xxxzz_xxyyy, \
                             tk_xxxzz_xxyyz, \
                             tk_xxxzz_xxyzz, \
                             tk_xxxzz_xxzzz, \
                             tk_xxxzz_xyyyy, \
                             tk_xxxzz_xyyyz, \
                             tk_xxxzz_xyyzz, \
                             tk_xxxzz_xyzzz, \
                             tk_xxxzz_xzzzz, \
                             tk_xxxzz_yyyyy, \
                             tk_xxxzz_yyyyz, \
                             tk_xxxzz_yyyzz, \
                             tk_xxxzz_yyzzz, \
                             tk_xxxzz_yzzzz, \
                             tk_xxxzz_zzzzz, \
                             tk_xxzz_xxxxz,  \
                             tk_xxzz_xxxyz,  \
                             tk_xxzz_xxxz,   \
                             tk_xxzz_xxxzz,  \
                             tk_xxzz_xxyyz,  \
                             tk_xxzz_xxyz,   \
                             tk_xxzz_xxyzz,  \
                             tk_xxzz_xxzz,   \
                             tk_xxzz_xxzzz,  \
                             tk_xxzz_xyyyz,  \
                             tk_xxzz_xyyz,   \
                             tk_xxzz_xyyzz,  \
                             tk_xxzz_xyzz,   \
                             tk_xxzz_xyzzz,  \
                             tk_xxzz_xzzz,   \
                             tk_xxzz_xzzzz,  \
                             tk_xxzz_yyyyy,  \
                             tk_xxzz_yyyyz,  \
                             tk_xxzz_yyyz,   \
                             tk_xxzz_yyyzz,  \
                             tk_xxzz_yyzz,   \
                             tk_xxzz_yyzzz,  \
                             tk_xxzz_yzzz,   \
                             tk_xxzz_yzzzz,  \
                             tk_xxzz_zzzz,   \
                             tk_xxzz_zzzzz,  \
                             tk_xzz_xxxxz,   \
                             tk_xzz_xxxyz,   \
                             tk_xzz_xxxzz,   \
                             tk_xzz_xxyyz,   \
                             tk_xzz_xxyzz,   \
                             tk_xzz_xxzzz,   \
                             tk_xzz_xyyyz,   \
                             tk_xzz_xyyzz,   \
                             tk_xzz_xyzzz,   \
                             tk_xzz_xzzzz,   \
                             tk_xzz_yyyyy,   \
                             tk_xzz_yyyyz,   \
                             tk_xzz_yyyzz,   \
                             tk_xzz_yyzzz,   \
                             tk_xzz_yzzzz,   \
                             tk_xzz_zzzzz,   \
                             ts_xxx_xxxxx,   \
                             ts_xxx_xxxxy,   \
                             ts_xxx_xxxyy,   \
                             ts_xxx_xxyyy,   \
                             ts_xxx_xyyyy,   \
                             ts_xxxzz_xxxxx, \
                             ts_xxxzz_xxxxy, \
                             ts_xxxzz_xxxxz, \
                             ts_xxxzz_xxxyy, \
                             ts_xxxzz_xxxyz, \
                             ts_xxxzz_xxxzz, \
                             ts_xxxzz_xxyyy, \
                             ts_xxxzz_xxyyz, \
                             ts_xxxzz_xxyzz, \
                             ts_xxxzz_xxzzz, \
                             ts_xxxzz_xyyyy, \
                             ts_xxxzz_xyyyz, \
                             ts_xxxzz_xyyzz, \
                             ts_xxxzz_xyzzz, \
                             ts_xxxzz_xzzzz, \
                             ts_xxxzz_yyyyy, \
                             ts_xxxzz_yyyyz, \
                             ts_xxxzz_yyyzz, \
                             ts_xxxzz_yyzzz, \
                             ts_xxxzz_yzzzz, \
                             ts_xxxzz_zzzzz, \
                             ts_xzz_xxxxz,   \
                             ts_xzz_xxxyz,   \
                             ts_xzz_xxxzz,   \
                             ts_xzz_xxyyz,   \
                             ts_xzz_xxyzz,   \
                             ts_xzz_xxzzz,   \
                             ts_xzz_xyyyz,   \
                             ts_xzz_xyyzz,   \
                             ts_xzz_xyzzz,   \
                             ts_xzz_xzzzz,   \
                             ts_xzz_yyyyy,   \
                             ts_xzz_yyyyz,   \
                             ts_xzz_yyyzz,   \
                             ts_xzz_yyzzz,   \
                             ts_xzz_yzzzz,   \
                             ts_xzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzz_xxxxx[i] =
            -2.0 * ts_xxx_xxxxx[i] * fbe_0 * fz_0 + tk_xxx_xxxxx[i] * fe_0 + tk_xxxz_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxxx[i] * fz_0;

        tk_xxxzz_xxxxy[i] =
            -2.0 * ts_xxx_xxxxy[i] * fbe_0 * fz_0 + tk_xxx_xxxxy[i] * fe_0 + tk_xxxz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxxy[i] * fz_0;

        tk_xxxzz_xxxxz[i] = -4.0 * ts_xzz_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxxz[i] * fe_0 + 4.0 * tk_xxzz_xxxz[i] * fe_0 +
                            tk_xxzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxxz[i] * fz_0;

        tk_xxxzz_xxxyy[i] =
            -2.0 * ts_xxx_xxxyy[i] * fbe_0 * fz_0 + tk_xxx_xxxyy[i] * fe_0 + tk_xxxz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxyy[i] * fz_0;

        tk_xxxzz_xxxyz[i] = -4.0 * ts_xzz_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxyz[i] * fe_0 + 3.0 * tk_xxzz_xxyz[i] * fe_0 +
                            tk_xxzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxyz[i] * fz_0;

        tk_xxxzz_xxxzz[i] = -4.0 * ts_xzz_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxzz[i] * fe_0 + 3.0 * tk_xxzz_xxzz[i] * fe_0 +
                            tk_xxzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxzz[i] * fz_0;

        tk_xxxzz_xxyyy[i] =
            -2.0 * ts_xxx_xxyyy[i] * fbe_0 * fz_0 + tk_xxx_xxyyy[i] * fe_0 + tk_xxxz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxyyy[i] * fz_0;

        tk_xxxzz_xxyyz[i] = -4.0 * ts_xzz_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxyyz[i] * fe_0 + 2.0 * tk_xxzz_xyyz[i] * fe_0 +
                            tk_xxzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxyyz[i] * fz_0;

        tk_xxxzz_xxyzz[i] = -4.0 * ts_xzz_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxyzz[i] * fe_0 + 2.0 * tk_xxzz_xyzz[i] * fe_0 +
                            tk_xxzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxyzz[i] * fz_0;

        tk_xxxzz_xxzzz[i] = -4.0 * ts_xzz_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxzzz[i] * fe_0 + 2.0 * tk_xxzz_xzzz[i] * fe_0 +
                            tk_xxzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxzzz[i] * fz_0;

        tk_xxxzz_xyyyy[i] =
            -2.0 * ts_xxx_xyyyy[i] * fbe_0 * fz_0 + tk_xxx_xyyyy[i] * fe_0 + tk_xxxz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xyyyy[i] * fz_0;

        tk_xxxzz_xyyyz[i] = -4.0 * ts_xzz_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyyyz[i] * fe_0 + tk_xxzz_yyyz[i] * fe_0 +
                            tk_xxzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyyyz[i] * fz_0;

        tk_xxxzz_xyyzz[i] = -4.0 * ts_xzz_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyyzz[i] * fe_0 + tk_xxzz_yyzz[i] * fe_0 +
                            tk_xxzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyyzz[i] * fz_0;

        tk_xxxzz_xyzzz[i] = -4.0 * ts_xzz_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyzzz[i] * fe_0 + tk_xxzz_yzzz[i] * fe_0 +
                            tk_xxzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyzzz[i] * fz_0;

        tk_xxxzz_xzzzz[i] = -4.0 * ts_xzz_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xzzzz[i] * fe_0 + tk_xxzz_zzzz[i] * fe_0 +
                            tk_xxzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xzzzz[i] * fz_0;

        tk_xxxzz_yyyyy[i] =
            -4.0 * ts_xzz_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyyy[i] * fe_0 + tk_xxzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyyy[i] * fz_0;

        tk_xxxzz_yyyyz[i] =
            -4.0 * ts_xzz_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyyz[i] * fe_0 + tk_xxzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyyz[i] * fz_0;

        tk_xxxzz_yyyzz[i] =
            -4.0 * ts_xzz_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyzz[i] * fe_0 + tk_xxzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyzz[i] * fz_0;

        tk_xxxzz_yyzzz[i] =
            -4.0 * ts_xzz_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyzzz[i] * fe_0 + tk_xxzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyzzz[i] * fz_0;

        tk_xxxzz_yzzzz[i] =
            -4.0 * ts_xzz_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yzzzz[i] * fe_0 + tk_xxzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yzzzz[i] * fz_0;

        tk_xxxzz_zzzzz[i] =
            -4.0 * ts_xzz_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_zzzzz[i] * fe_0 + tk_xxzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_zzzzz[i] * fz_0;
    }

    // Set up 126-147 components of targeted buffer : HH

    auto tk_xxyyy_xxxxx = pbuffer.data(idx_kin_hh + 126);

    auto tk_xxyyy_xxxxy = pbuffer.data(idx_kin_hh + 127);

    auto tk_xxyyy_xxxxz = pbuffer.data(idx_kin_hh + 128);

    auto tk_xxyyy_xxxyy = pbuffer.data(idx_kin_hh + 129);

    auto tk_xxyyy_xxxyz = pbuffer.data(idx_kin_hh + 130);

    auto tk_xxyyy_xxxzz = pbuffer.data(idx_kin_hh + 131);

    auto tk_xxyyy_xxyyy = pbuffer.data(idx_kin_hh + 132);

    auto tk_xxyyy_xxyyz = pbuffer.data(idx_kin_hh + 133);

    auto tk_xxyyy_xxyzz = pbuffer.data(idx_kin_hh + 134);

    auto tk_xxyyy_xxzzz = pbuffer.data(idx_kin_hh + 135);

    auto tk_xxyyy_xyyyy = pbuffer.data(idx_kin_hh + 136);

    auto tk_xxyyy_xyyyz = pbuffer.data(idx_kin_hh + 137);

    auto tk_xxyyy_xyyzz = pbuffer.data(idx_kin_hh + 138);

    auto tk_xxyyy_xyzzz = pbuffer.data(idx_kin_hh + 139);

    auto tk_xxyyy_xzzzz = pbuffer.data(idx_kin_hh + 140);

    auto tk_xxyyy_yyyyy = pbuffer.data(idx_kin_hh + 141);

    auto tk_xxyyy_yyyyz = pbuffer.data(idx_kin_hh + 142);

    auto tk_xxyyy_yyyzz = pbuffer.data(idx_kin_hh + 143);

    auto tk_xxyyy_yyzzz = pbuffer.data(idx_kin_hh + 144);

    auto tk_xxyyy_yzzzz = pbuffer.data(idx_kin_hh + 145);

    auto tk_xxyyy_zzzzz = pbuffer.data(idx_kin_hh + 146);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xxy_xxxxx,   \
                             tk_xxy_xxxxz,   \
                             tk_xxy_xxxzz,   \
                             tk_xxy_xxzzz,   \
                             tk_xxy_xzzzz,   \
                             tk_xxyy_xxxxx,  \
                             tk_xxyy_xxxxz,  \
                             tk_xxyy_xxxzz,  \
                             tk_xxyy_xxzzz,  \
                             tk_xxyy_xzzzz,  \
                             tk_xxyyy_xxxxx, \
                             tk_xxyyy_xxxxy, \
                             tk_xxyyy_xxxxz, \
                             tk_xxyyy_xxxyy, \
                             tk_xxyyy_xxxyz, \
                             tk_xxyyy_xxxzz, \
                             tk_xxyyy_xxyyy, \
                             tk_xxyyy_xxyyz, \
                             tk_xxyyy_xxyzz, \
                             tk_xxyyy_xxzzz, \
                             tk_xxyyy_xyyyy, \
                             tk_xxyyy_xyyyz, \
                             tk_xxyyy_xyyzz, \
                             tk_xxyyy_xyzzz, \
                             tk_xxyyy_xzzzz, \
                             tk_xxyyy_yyyyy, \
                             tk_xxyyy_yyyyz, \
                             tk_xxyyy_yyyzz, \
                             tk_xxyyy_yyzzz, \
                             tk_xxyyy_yzzzz, \
                             tk_xxyyy_zzzzz, \
                             tk_xyyy_xxxxy,  \
                             tk_xyyy_xxxy,   \
                             tk_xyyy_xxxyy,  \
                             tk_xyyy_xxxyz,  \
                             tk_xyyy_xxyy,   \
                             tk_xyyy_xxyyy,  \
                             tk_xyyy_xxyyz,  \
                             tk_xyyy_xxyz,   \
                             tk_xyyy_xxyzz,  \
                             tk_xyyy_xyyy,   \
                             tk_xyyy_xyyyy,  \
                             tk_xyyy_xyyyz,  \
                             tk_xyyy_xyyz,   \
                             tk_xyyy_xyyzz,  \
                             tk_xyyy_xyzz,   \
                             tk_xyyy_xyzzz,  \
                             tk_xyyy_yyyy,   \
                             tk_xyyy_yyyyy,  \
                             tk_xyyy_yyyyz,  \
                             tk_xyyy_yyyz,   \
                             tk_xyyy_yyyzz,  \
                             tk_xyyy_yyzz,   \
                             tk_xyyy_yyzzz,  \
                             tk_xyyy_yzzz,   \
                             tk_xyyy_yzzzz,  \
                             tk_xyyy_zzzzz,  \
                             tk_yyy_xxxxy,   \
                             tk_yyy_xxxyy,   \
                             tk_yyy_xxxyz,   \
                             tk_yyy_xxyyy,   \
                             tk_yyy_xxyyz,   \
                             tk_yyy_xxyzz,   \
                             tk_yyy_xyyyy,   \
                             tk_yyy_xyyyz,   \
                             tk_yyy_xyyzz,   \
                             tk_yyy_xyzzz,   \
                             tk_yyy_yyyyy,   \
                             tk_yyy_yyyyz,   \
                             tk_yyy_yyyzz,   \
                             tk_yyy_yyzzz,   \
                             tk_yyy_yzzzz,   \
                             tk_yyy_zzzzz,   \
                             ts_xxy_xxxxx,   \
                             ts_xxy_xxxxz,   \
                             ts_xxy_xxxzz,   \
                             ts_xxy_xxzzz,   \
                             ts_xxy_xzzzz,   \
                             ts_xxyyy_xxxxx, \
                             ts_xxyyy_xxxxy, \
                             ts_xxyyy_xxxxz, \
                             ts_xxyyy_xxxyy, \
                             ts_xxyyy_xxxyz, \
                             ts_xxyyy_xxxzz, \
                             ts_xxyyy_xxyyy, \
                             ts_xxyyy_xxyyz, \
                             ts_xxyyy_xxyzz, \
                             ts_xxyyy_xxzzz, \
                             ts_xxyyy_xyyyy, \
                             ts_xxyyy_xyyyz, \
                             ts_xxyyy_xyyzz, \
                             ts_xxyyy_xyzzz, \
                             ts_xxyyy_xzzzz, \
                             ts_xxyyy_yyyyy, \
                             ts_xxyyy_yyyyz, \
                             ts_xxyyy_yyyzz, \
                             ts_xxyyy_yyzzz, \
                             ts_xxyyy_yzzzz, \
                             ts_xxyyy_zzzzz, \
                             ts_yyy_xxxxy,   \
                             ts_yyy_xxxyy,   \
                             ts_yyy_xxxyz,   \
                             ts_yyy_xxyyy,   \
                             ts_yyy_xxyyz,   \
                             ts_yyy_xxyzz,   \
                             ts_yyy_xyyyy,   \
                             ts_yyy_xyyyz,   \
                             ts_yyy_xyyzz,   \
                             ts_yyy_xyzzz,   \
                             ts_yyy_yyyyy,   \
                             ts_yyy_yyyyz,   \
                             ts_yyy_yyyzz,   \
                             ts_yyy_yyzzz,   \
                             ts_yyy_yzzzz,   \
                             ts_yyy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyy_xxxxx[i] =
            -4.0 * ts_xxy_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxxx[i] * fe_0 + tk_xxyy_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxxx[i] * fz_0;

        tk_xxyyy_xxxxy[i] = -2.0 * ts_yyy_xxxxy[i] * fbe_0 * fz_0 + tk_yyy_xxxxy[i] * fe_0 + 4.0 * tk_xyyy_xxxy[i] * fe_0 +
                            tk_xyyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxxy[i] * fz_0;

        tk_xxyyy_xxxxz[i] =
            -4.0 * ts_xxy_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxxz[i] * fe_0 + tk_xxyy_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxxz[i] * fz_0;

        tk_xxyyy_xxxyy[i] = -2.0 * ts_yyy_xxxyy[i] * fbe_0 * fz_0 + tk_yyy_xxxyy[i] * fe_0 + 3.0 * tk_xyyy_xxyy[i] * fe_0 +
                            tk_xyyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxyy[i] * fz_0;

        tk_xxyyy_xxxyz[i] = -2.0 * ts_yyy_xxxyz[i] * fbe_0 * fz_0 + tk_yyy_xxxyz[i] * fe_0 + 3.0 * tk_xyyy_xxyz[i] * fe_0 +
                            tk_xyyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxyz[i] * fz_0;

        tk_xxyyy_xxxzz[i] =
            -4.0 * ts_xxy_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxzz[i] * fe_0 + tk_xxyy_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxzz[i] * fz_0;

        tk_xxyyy_xxyyy[i] = -2.0 * ts_yyy_xxyyy[i] * fbe_0 * fz_0 + tk_yyy_xxyyy[i] * fe_0 + 2.0 * tk_xyyy_xyyy[i] * fe_0 +
                            tk_xyyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyyy[i] * fz_0;

        tk_xxyyy_xxyyz[i] = -2.0 * ts_yyy_xxyyz[i] * fbe_0 * fz_0 + tk_yyy_xxyyz[i] * fe_0 + 2.0 * tk_xyyy_xyyz[i] * fe_0 +
                            tk_xyyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyyz[i] * fz_0;

        tk_xxyyy_xxyzz[i] = -2.0 * ts_yyy_xxyzz[i] * fbe_0 * fz_0 + tk_yyy_xxyzz[i] * fe_0 + 2.0 * tk_xyyy_xyzz[i] * fe_0 +
                            tk_xyyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyzz[i] * fz_0;

        tk_xxyyy_xxzzz[i] =
            -4.0 * ts_xxy_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxzzz[i] * fe_0 + tk_xxyy_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxzzz[i] * fz_0;

        tk_xxyyy_xyyyy[i] = -2.0 * ts_yyy_xyyyy[i] * fbe_0 * fz_0 + tk_yyy_xyyyy[i] * fe_0 + tk_xyyy_yyyy[i] * fe_0 + tk_xyyy_xyyyy[i] * pa_x[i] +
                            2.0 * ts_xxyyy_xyyyy[i] * fz_0;

        tk_xxyyy_xyyyz[i] = -2.0 * ts_yyy_xyyyz[i] * fbe_0 * fz_0 + tk_yyy_xyyyz[i] * fe_0 + tk_xyyy_yyyz[i] * fe_0 + tk_xyyy_xyyyz[i] * pa_x[i] +
                            2.0 * ts_xxyyy_xyyyz[i] * fz_0;

        tk_xxyyy_xyyzz[i] = -2.0 * ts_yyy_xyyzz[i] * fbe_0 * fz_0 + tk_yyy_xyyzz[i] * fe_0 + tk_xyyy_yyzz[i] * fe_0 + tk_xyyy_xyyzz[i] * pa_x[i] +
                            2.0 * ts_xxyyy_xyyzz[i] * fz_0;

        tk_xxyyy_xyzzz[i] = -2.0 * ts_yyy_xyzzz[i] * fbe_0 * fz_0 + tk_yyy_xyzzz[i] * fe_0 + tk_xyyy_yzzz[i] * fe_0 + tk_xyyy_xyzzz[i] * pa_x[i] +
                            2.0 * ts_xxyyy_xyzzz[i] * fz_0;

        tk_xxyyy_xzzzz[i] =
            -4.0 * ts_xxy_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xzzzz[i] * fe_0 + tk_xxyy_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xzzzz[i] * fz_0;

        tk_xxyyy_yyyyy[i] =
            -2.0 * ts_yyy_yyyyy[i] * fbe_0 * fz_0 + tk_yyy_yyyyy[i] * fe_0 + tk_xyyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyyy[i] * fz_0;

        tk_xxyyy_yyyyz[i] =
            -2.0 * ts_yyy_yyyyz[i] * fbe_0 * fz_0 + tk_yyy_yyyyz[i] * fe_0 + tk_xyyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyyz[i] * fz_0;

        tk_xxyyy_yyyzz[i] =
            -2.0 * ts_yyy_yyyzz[i] * fbe_0 * fz_0 + tk_yyy_yyyzz[i] * fe_0 + tk_xyyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyzz[i] * fz_0;

        tk_xxyyy_yyzzz[i] =
            -2.0 * ts_yyy_yyzzz[i] * fbe_0 * fz_0 + tk_yyy_yyzzz[i] * fe_0 + tk_xyyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyzzz[i] * fz_0;

        tk_xxyyy_yzzzz[i] =
            -2.0 * ts_yyy_yzzzz[i] * fbe_0 * fz_0 + tk_yyy_yzzzz[i] * fe_0 + tk_xyyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yzzzz[i] * fz_0;

        tk_xxyyy_zzzzz[i] =
            -2.0 * ts_yyy_zzzzz[i] * fbe_0 * fz_0 + tk_yyy_zzzzz[i] * fe_0 + tk_xyyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_zzzzz[i] * fz_0;
    }

    // Set up 147-168 components of targeted buffer : HH

    auto tk_xxyyz_xxxxx = pbuffer.data(idx_kin_hh + 147);

    auto tk_xxyyz_xxxxy = pbuffer.data(idx_kin_hh + 148);

    auto tk_xxyyz_xxxxz = pbuffer.data(idx_kin_hh + 149);

    auto tk_xxyyz_xxxyy = pbuffer.data(idx_kin_hh + 150);

    auto tk_xxyyz_xxxyz = pbuffer.data(idx_kin_hh + 151);

    auto tk_xxyyz_xxxzz = pbuffer.data(idx_kin_hh + 152);

    auto tk_xxyyz_xxyyy = pbuffer.data(idx_kin_hh + 153);

    auto tk_xxyyz_xxyyz = pbuffer.data(idx_kin_hh + 154);

    auto tk_xxyyz_xxyzz = pbuffer.data(idx_kin_hh + 155);

    auto tk_xxyyz_xxzzz = pbuffer.data(idx_kin_hh + 156);

    auto tk_xxyyz_xyyyy = pbuffer.data(idx_kin_hh + 157);

    auto tk_xxyyz_xyyyz = pbuffer.data(idx_kin_hh + 158);

    auto tk_xxyyz_xyyzz = pbuffer.data(idx_kin_hh + 159);

    auto tk_xxyyz_xyzzz = pbuffer.data(idx_kin_hh + 160);

    auto tk_xxyyz_xzzzz = pbuffer.data(idx_kin_hh + 161);

    auto tk_xxyyz_yyyyy = pbuffer.data(idx_kin_hh + 162);

    auto tk_xxyyz_yyyyz = pbuffer.data(idx_kin_hh + 163);

    auto tk_xxyyz_yyyzz = pbuffer.data(idx_kin_hh + 164);

    auto tk_xxyyz_yyzzz = pbuffer.data(idx_kin_hh + 165);

    auto tk_xxyyz_yzzzz = pbuffer.data(idx_kin_hh + 166);

    auto tk_xxyyz_zzzzz = pbuffer.data(idx_kin_hh + 167);

#pragma omp simd aligned(pa_z,               \
                             tk_xxyy_xxxx,   \
                             tk_xxyy_xxxxx,  \
                             tk_xxyy_xxxxy,  \
                             tk_xxyy_xxxxz,  \
                             tk_xxyy_xxxy,   \
                             tk_xxyy_xxxyy,  \
                             tk_xxyy_xxxyz,  \
                             tk_xxyy_xxxz,   \
                             tk_xxyy_xxxzz,  \
                             tk_xxyy_xxyy,   \
                             tk_xxyy_xxyyy,  \
                             tk_xxyy_xxyyz,  \
                             tk_xxyy_xxyz,   \
                             tk_xxyy_xxyzz,  \
                             tk_xxyy_xxzz,   \
                             tk_xxyy_xxzzz,  \
                             tk_xxyy_xyyy,   \
                             tk_xxyy_xyyyy,  \
                             tk_xxyy_xyyyz,  \
                             tk_xxyy_xyyz,   \
                             tk_xxyy_xyyzz,  \
                             tk_xxyy_xyzz,   \
                             tk_xxyy_xyzzz,  \
                             tk_xxyy_xzzz,   \
                             tk_xxyy_xzzzz,  \
                             tk_xxyy_yyyy,   \
                             tk_xxyy_yyyyy,  \
                             tk_xxyy_yyyyz,  \
                             tk_xxyy_yyyz,   \
                             tk_xxyy_yyyzz,  \
                             tk_xxyy_yyzz,   \
                             tk_xxyy_yyzzz,  \
                             tk_xxyy_yzzz,   \
                             tk_xxyy_yzzzz,  \
                             tk_xxyy_zzzz,   \
                             tk_xxyy_zzzzz,  \
                             tk_xxyyz_xxxxx, \
                             tk_xxyyz_xxxxy, \
                             tk_xxyyz_xxxxz, \
                             tk_xxyyz_xxxyy, \
                             tk_xxyyz_xxxyz, \
                             tk_xxyyz_xxxzz, \
                             tk_xxyyz_xxyyy, \
                             tk_xxyyz_xxyyz, \
                             tk_xxyyz_xxyzz, \
                             tk_xxyyz_xxzzz, \
                             tk_xxyyz_xyyyy, \
                             tk_xxyyz_xyyyz, \
                             tk_xxyyz_xyyzz, \
                             tk_xxyyz_xyzzz, \
                             tk_xxyyz_xzzzz, \
                             tk_xxyyz_yyyyy, \
                             tk_xxyyz_yyyyz, \
                             tk_xxyyz_yyyzz, \
                             tk_xxyyz_yyzzz, \
                             tk_xxyyz_yzzzz, \
                             tk_xxyyz_zzzzz, \
                             ts_xxyyz_xxxxx, \
                             ts_xxyyz_xxxxy, \
                             ts_xxyyz_xxxxz, \
                             ts_xxyyz_xxxyy, \
                             ts_xxyyz_xxxyz, \
                             ts_xxyyz_xxxzz, \
                             ts_xxyyz_xxyyy, \
                             ts_xxyyz_xxyyz, \
                             ts_xxyyz_xxyzz, \
                             ts_xxyyz_xxzzz, \
                             ts_xxyyz_xyyyy, \
                             ts_xxyyz_xyyyz, \
                             ts_xxyyz_xyyzz, \
                             ts_xxyyz_xyzzz, \
                             ts_xxyyz_xzzzz, \
                             ts_xxyyz_yyyyy, \
                             ts_xxyyz_yyyyz, \
                             ts_xxyyz_yyyzz, \
                             ts_xxyyz_yyzzz, \
                             ts_xxyyz_yzzzz, \
                             ts_xxyyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyz_xxxxx[i] = tk_xxyy_xxxxx[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxx[i] * fz_0;

        tk_xxyyz_xxxxy[i] = tk_xxyy_xxxxy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxy[i] * fz_0;

        tk_xxyyz_xxxxz[i] = tk_xxyy_xxxx[i] * fe_0 + tk_xxyy_xxxxz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxz[i] * fz_0;

        tk_xxyyz_xxxyy[i] = tk_xxyy_xxxyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxyy[i] * fz_0;

        tk_xxyyz_xxxyz[i] = tk_xxyy_xxxy[i] * fe_0 + tk_xxyy_xxxyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxyz[i] * fz_0;

        tk_xxyyz_xxxzz[i] = 2.0 * tk_xxyy_xxxz[i] * fe_0 + tk_xxyy_xxxzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxzz[i] * fz_0;

        tk_xxyyz_xxyyy[i] = tk_xxyy_xxyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyyy[i] * fz_0;

        tk_xxyyz_xxyyz[i] = tk_xxyy_xxyy[i] * fe_0 + tk_xxyy_xxyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyyz[i] * fz_0;

        tk_xxyyz_xxyzz[i] = 2.0 * tk_xxyy_xxyz[i] * fe_0 + tk_xxyy_xxyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyzz[i] * fz_0;

        tk_xxyyz_xxzzz[i] = 3.0 * tk_xxyy_xxzz[i] * fe_0 + tk_xxyy_xxzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxzzz[i] * fz_0;

        tk_xxyyz_xyyyy[i] = tk_xxyy_xyyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyyy[i] * fz_0;

        tk_xxyyz_xyyyz[i] = tk_xxyy_xyyy[i] * fe_0 + tk_xxyy_xyyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyyz[i] * fz_0;

        tk_xxyyz_xyyzz[i] = 2.0 * tk_xxyy_xyyz[i] * fe_0 + tk_xxyy_xyyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyzz[i] * fz_0;

        tk_xxyyz_xyzzz[i] = 3.0 * tk_xxyy_xyzz[i] * fe_0 + tk_xxyy_xyzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyzzz[i] * fz_0;

        tk_xxyyz_xzzzz[i] = 4.0 * tk_xxyy_xzzz[i] * fe_0 + tk_xxyy_xzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xzzzz[i] * fz_0;

        tk_xxyyz_yyyyy[i] = tk_xxyy_yyyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyyy[i] * fz_0;

        tk_xxyyz_yyyyz[i] = tk_xxyy_yyyy[i] * fe_0 + tk_xxyy_yyyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyyz[i] * fz_0;

        tk_xxyyz_yyyzz[i] = 2.0 * tk_xxyy_yyyz[i] * fe_0 + tk_xxyy_yyyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyzz[i] * fz_0;

        tk_xxyyz_yyzzz[i] = 3.0 * tk_xxyy_yyzz[i] * fe_0 + tk_xxyy_yyzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyzzz[i] * fz_0;

        tk_xxyyz_yzzzz[i] = 4.0 * tk_xxyy_yzzz[i] * fe_0 + tk_xxyy_yzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yzzzz[i] * fz_0;

        tk_xxyyz_zzzzz[i] = 5.0 * tk_xxyy_zzzz[i] * fe_0 + tk_xxyy_zzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_zzzzz[i] * fz_0;
    }

    // Set up 168-189 components of targeted buffer : HH

    auto tk_xxyzz_xxxxx = pbuffer.data(idx_kin_hh + 168);

    auto tk_xxyzz_xxxxy = pbuffer.data(idx_kin_hh + 169);

    auto tk_xxyzz_xxxxz = pbuffer.data(idx_kin_hh + 170);

    auto tk_xxyzz_xxxyy = pbuffer.data(idx_kin_hh + 171);

    auto tk_xxyzz_xxxyz = pbuffer.data(idx_kin_hh + 172);

    auto tk_xxyzz_xxxzz = pbuffer.data(idx_kin_hh + 173);

    auto tk_xxyzz_xxyyy = pbuffer.data(idx_kin_hh + 174);

    auto tk_xxyzz_xxyyz = pbuffer.data(idx_kin_hh + 175);

    auto tk_xxyzz_xxyzz = pbuffer.data(idx_kin_hh + 176);

    auto tk_xxyzz_xxzzz = pbuffer.data(idx_kin_hh + 177);

    auto tk_xxyzz_xyyyy = pbuffer.data(idx_kin_hh + 178);

    auto tk_xxyzz_xyyyz = pbuffer.data(idx_kin_hh + 179);

    auto tk_xxyzz_xyyzz = pbuffer.data(idx_kin_hh + 180);

    auto tk_xxyzz_xyzzz = pbuffer.data(idx_kin_hh + 181);

    auto tk_xxyzz_xzzzz = pbuffer.data(idx_kin_hh + 182);

    auto tk_xxyzz_yyyyy = pbuffer.data(idx_kin_hh + 183);

    auto tk_xxyzz_yyyyz = pbuffer.data(idx_kin_hh + 184);

    auto tk_xxyzz_yyyzz = pbuffer.data(idx_kin_hh + 185);

    auto tk_xxyzz_yyzzz = pbuffer.data(idx_kin_hh + 186);

    auto tk_xxyzz_yzzzz = pbuffer.data(idx_kin_hh + 187);

    auto tk_xxyzz_zzzzz = pbuffer.data(idx_kin_hh + 188);

#pragma omp simd aligned(pa_y,               \
                             tk_xxyzz_xxxxx, \
                             tk_xxyzz_xxxxy, \
                             tk_xxyzz_xxxxz, \
                             tk_xxyzz_xxxyy, \
                             tk_xxyzz_xxxyz, \
                             tk_xxyzz_xxxzz, \
                             tk_xxyzz_xxyyy, \
                             tk_xxyzz_xxyyz, \
                             tk_xxyzz_xxyzz, \
                             tk_xxyzz_xxzzz, \
                             tk_xxyzz_xyyyy, \
                             tk_xxyzz_xyyyz, \
                             tk_xxyzz_xyyzz, \
                             tk_xxyzz_xyzzz, \
                             tk_xxyzz_xzzzz, \
                             tk_xxyzz_yyyyy, \
                             tk_xxyzz_yyyyz, \
                             tk_xxyzz_yyyzz, \
                             tk_xxyzz_yyzzz, \
                             tk_xxyzz_yzzzz, \
                             tk_xxyzz_zzzzz, \
                             tk_xxzz_xxxx,   \
                             tk_xxzz_xxxxx,  \
                             tk_xxzz_xxxxy,  \
                             tk_xxzz_xxxxz,  \
                             tk_xxzz_xxxy,   \
                             tk_xxzz_xxxyy,  \
                             tk_xxzz_xxxyz,  \
                             tk_xxzz_xxxz,   \
                             tk_xxzz_xxxzz,  \
                             tk_xxzz_xxyy,   \
                             tk_xxzz_xxyyy,  \
                             tk_xxzz_xxyyz,  \
                             tk_xxzz_xxyz,   \
                             tk_xxzz_xxyzz,  \
                             tk_xxzz_xxzz,   \
                             tk_xxzz_xxzzz,  \
                             tk_xxzz_xyyy,   \
                             tk_xxzz_xyyyy,  \
                             tk_xxzz_xyyyz,  \
                             tk_xxzz_xyyz,   \
                             tk_xxzz_xyyzz,  \
                             tk_xxzz_xyzz,   \
                             tk_xxzz_xyzzz,  \
                             tk_xxzz_xzzz,   \
                             tk_xxzz_xzzzz,  \
                             tk_xxzz_yyyy,   \
                             tk_xxzz_yyyyy,  \
                             tk_xxzz_yyyyz,  \
                             tk_xxzz_yyyz,   \
                             tk_xxzz_yyyzz,  \
                             tk_xxzz_yyzz,   \
                             tk_xxzz_yyzzz,  \
                             tk_xxzz_yzzz,   \
                             tk_xxzz_yzzzz,  \
                             tk_xxzz_zzzz,   \
                             tk_xxzz_zzzzz,  \
                             ts_xxyzz_xxxxx, \
                             ts_xxyzz_xxxxy, \
                             ts_xxyzz_xxxxz, \
                             ts_xxyzz_xxxyy, \
                             ts_xxyzz_xxxyz, \
                             ts_xxyzz_xxxzz, \
                             ts_xxyzz_xxyyy, \
                             ts_xxyzz_xxyyz, \
                             ts_xxyzz_xxyzz, \
                             ts_xxyzz_xxzzz, \
                             ts_xxyzz_xyyyy, \
                             ts_xxyzz_xyyyz, \
                             ts_xxyzz_xyyzz, \
                             ts_xxyzz_xyzzz, \
                             ts_xxyzz_xzzzz, \
                             ts_xxyzz_yyyyy, \
                             ts_xxyzz_yyyyz, \
                             ts_xxyzz_yyyzz, \
                             ts_xxyzz_yyzzz, \
                             ts_xxyzz_yzzzz, \
                             ts_xxyzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzz_xxxxx[i] = tk_xxzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxx[i] * fz_0;

        tk_xxyzz_xxxxy[i] = tk_xxzz_xxxx[i] * fe_0 + tk_xxzz_xxxxy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxy[i] * fz_0;

        tk_xxyzz_xxxxz[i] = tk_xxzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxz[i] * fz_0;

        tk_xxyzz_xxxyy[i] = 2.0 * tk_xxzz_xxxy[i] * fe_0 + tk_xxzz_xxxyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxyy[i] * fz_0;

        tk_xxyzz_xxxyz[i] = tk_xxzz_xxxz[i] * fe_0 + tk_xxzz_xxxyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxyz[i] * fz_0;

        tk_xxyzz_xxxzz[i] = tk_xxzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxzz[i] * fz_0;

        tk_xxyzz_xxyyy[i] = 3.0 * tk_xxzz_xxyy[i] * fe_0 + tk_xxzz_xxyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyyy[i] * fz_0;

        tk_xxyzz_xxyyz[i] = 2.0 * tk_xxzz_xxyz[i] * fe_0 + tk_xxzz_xxyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyyz[i] * fz_0;

        tk_xxyzz_xxyzz[i] = tk_xxzz_xxzz[i] * fe_0 + tk_xxzz_xxyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyzz[i] * fz_0;

        tk_xxyzz_xxzzz[i] = tk_xxzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxzzz[i] * fz_0;

        tk_xxyzz_xyyyy[i] = 4.0 * tk_xxzz_xyyy[i] * fe_0 + tk_xxzz_xyyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyyy[i] * fz_0;

        tk_xxyzz_xyyyz[i] = 3.0 * tk_xxzz_xyyz[i] * fe_0 + tk_xxzz_xyyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyyz[i] * fz_0;

        tk_xxyzz_xyyzz[i] = 2.0 * tk_xxzz_xyzz[i] * fe_0 + tk_xxzz_xyyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyzz[i] * fz_0;

        tk_xxyzz_xyzzz[i] = tk_xxzz_xzzz[i] * fe_0 + tk_xxzz_xyzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyzzz[i] * fz_0;

        tk_xxyzz_xzzzz[i] = tk_xxzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xzzzz[i] * fz_0;

        tk_xxyzz_yyyyy[i] = 5.0 * tk_xxzz_yyyy[i] * fe_0 + tk_xxzz_yyyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyyy[i] * fz_0;

        tk_xxyzz_yyyyz[i] = 4.0 * tk_xxzz_yyyz[i] * fe_0 + tk_xxzz_yyyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyyz[i] * fz_0;

        tk_xxyzz_yyyzz[i] = 3.0 * tk_xxzz_yyzz[i] * fe_0 + tk_xxzz_yyyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyzz[i] * fz_0;

        tk_xxyzz_yyzzz[i] = 2.0 * tk_xxzz_yzzz[i] * fe_0 + tk_xxzz_yyzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyzzz[i] * fz_0;

        tk_xxyzz_yzzzz[i] = tk_xxzz_zzzz[i] * fe_0 + tk_xxzz_yzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yzzzz[i] * fz_0;

        tk_xxyzz_zzzzz[i] = tk_xxzz_zzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_zzzzz[i] * fz_0;
    }

    // Set up 189-210 components of targeted buffer : HH

    auto tk_xxzzz_xxxxx = pbuffer.data(idx_kin_hh + 189);

    auto tk_xxzzz_xxxxy = pbuffer.data(idx_kin_hh + 190);

    auto tk_xxzzz_xxxxz = pbuffer.data(idx_kin_hh + 191);

    auto tk_xxzzz_xxxyy = pbuffer.data(idx_kin_hh + 192);

    auto tk_xxzzz_xxxyz = pbuffer.data(idx_kin_hh + 193);

    auto tk_xxzzz_xxxzz = pbuffer.data(idx_kin_hh + 194);

    auto tk_xxzzz_xxyyy = pbuffer.data(idx_kin_hh + 195);

    auto tk_xxzzz_xxyyz = pbuffer.data(idx_kin_hh + 196);

    auto tk_xxzzz_xxyzz = pbuffer.data(idx_kin_hh + 197);

    auto tk_xxzzz_xxzzz = pbuffer.data(idx_kin_hh + 198);

    auto tk_xxzzz_xyyyy = pbuffer.data(idx_kin_hh + 199);

    auto tk_xxzzz_xyyyz = pbuffer.data(idx_kin_hh + 200);

    auto tk_xxzzz_xyyzz = pbuffer.data(idx_kin_hh + 201);

    auto tk_xxzzz_xyzzz = pbuffer.data(idx_kin_hh + 202);

    auto tk_xxzzz_xzzzz = pbuffer.data(idx_kin_hh + 203);

    auto tk_xxzzz_yyyyy = pbuffer.data(idx_kin_hh + 204);

    auto tk_xxzzz_yyyyz = pbuffer.data(idx_kin_hh + 205);

    auto tk_xxzzz_yyyzz = pbuffer.data(idx_kin_hh + 206);

    auto tk_xxzzz_yyzzz = pbuffer.data(idx_kin_hh + 207);

    auto tk_xxzzz_yzzzz = pbuffer.data(idx_kin_hh + 208);

    auto tk_xxzzz_zzzzz = pbuffer.data(idx_kin_hh + 209);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xxz_xxxxx,   \
                             tk_xxz_xxxxy,   \
                             tk_xxz_xxxyy,   \
                             tk_xxz_xxyyy,   \
                             tk_xxz_xyyyy,   \
                             tk_xxzz_xxxxx,  \
                             tk_xxzz_xxxxy,  \
                             tk_xxzz_xxxyy,  \
                             tk_xxzz_xxyyy,  \
                             tk_xxzz_xyyyy,  \
                             tk_xxzzz_xxxxx, \
                             tk_xxzzz_xxxxy, \
                             tk_xxzzz_xxxxz, \
                             tk_xxzzz_xxxyy, \
                             tk_xxzzz_xxxyz, \
                             tk_xxzzz_xxxzz, \
                             tk_xxzzz_xxyyy, \
                             tk_xxzzz_xxyyz, \
                             tk_xxzzz_xxyzz, \
                             tk_xxzzz_xxzzz, \
                             tk_xxzzz_xyyyy, \
                             tk_xxzzz_xyyyz, \
                             tk_xxzzz_xyyzz, \
                             tk_xxzzz_xyzzz, \
                             tk_xxzzz_xzzzz, \
                             tk_xxzzz_yyyyy, \
                             tk_xxzzz_yyyyz, \
                             tk_xxzzz_yyyzz, \
                             tk_xxzzz_yyzzz, \
                             tk_xxzzz_yzzzz, \
                             tk_xxzzz_zzzzz, \
                             tk_xzzz_xxxxz,  \
                             tk_xzzz_xxxyz,  \
                             tk_xzzz_xxxz,   \
                             tk_xzzz_xxxzz,  \
                             tk_xzzz_xxyyz,  \
                             tk_xzzz_xxyz,   \
                             tk_xzzz_xxyzz,  \
                             tk_xzzz_xxzz,   \
                             tk_xzzz_xxzzz,  \
                             tk_xzzz_xyyyz,  \
                             tk_xzzz_xyyz,   \
                             tk_xzzz_xyyzz,  \
                             tk_xzzz_xyzz,   \
                             tk_xzzz_xyzzz,  \
                             tk_xzzz_xzzz,   \
                             tk_xzzz_xzzzz,  \
                             tk_xzzz_yyyyy,  \
                             tk_xzzz_yyyyz,  \
                             tk_xzzz_yyyz,   \
                             tk_xzzz_yyyzz,  \
                             tk_xzzz_yyzz,   \
                             tk_xzzz_yyzzz,  \
                             tk_xzzz_yzzz,   \
                             tk_xzzz_yzzzz,  \
                             tk_xzzz_zzzz,   \
                             tk_xzzz_zzzzz,  \
                             tk_zzz_xxxxz,   \
                             tk_zzz_xxxyz,   \
                             tk_zzz_xxxzz,   \
                             tk_zzz_xxyyz,   \
                             tk_zzz_xxyzz,   \
                             tk_zzz_xxzzz,   \
                             tk_zzz_xyyyz,   \
                             tk_zzz_xyyzz,   \
                             tk_zzz_xyzzz,   \
                             tk_zzz_xzzzz,   \
                             tk_zzz_yyyyy,   \
                             tk_zzz_yyyyz,   \
                             tk_zzz_yyyzz,   \
                             tk_zzz_yyzzz,   \
                             tk_zzz_yzzzz,   \
                             tk_zzz_zzzzz,   \
                             ts_xxz_xxxxx,   \
                             ts_xxz_xxxxy,   \
                             ts_xxz_xxxyy,   \
                             ts_xxz_xxyyy,   \
                             ts_xxz_xyyyy,   \
                             ts_xxzzz_xxxxx, \
                             ts_xxzzz_xxxxy, \
                             ts_xxzzz_xxxxz, \
                             ts_xxzzz_xxxyy, \
                             ts_xxzzz_xxxyz, \
                             ts_xxzzz_xxxzz, \
                             ts_xxzzz_xxyyy, \
                             ts_xxzzz_xxyyz, \
                             ts_xxzzz_xxyzz, \
                             ts_xxzzz_xxzzz, \
                             ts_xxzzz_xyyyy, \
                             ts_xxzzz_xyyyz, \
                             ts_xxzzz_xyyzz, \
                             ts_xxzzz_xyzzz, \
                             ts_xxzzz_xzzzz, \
                             ts_xxzzz_yyyyy, \
                             ts_xxzzz_yyyyz, \
                             ts_xxzzz_yyyzz, \
                             ts_xxzzz_yyzzz, \
                             ts_xxzzz_yzzzz, \
                             ts_xxzzz_zzzzz, \
                             ts_zzz_xxxxz,   \
                             ts_zzz_xxxyz,   \
                             ts_zzz_xxxzz,   \
                             ts_zzz_xxyyz,   \
                             ts_zzz_xxyzz,   \
                             ts_zzz_xxzzz,   \
                             ts_zzz_xyyyz,   \
                             ts_zzz_xyyzz,   \
                             ts_zzz_xyzzz,   \
                             ts_zzz_xzzzz,   \
                             ts_zzz_yyyyy,   \
                             ts_zzz_yyyyz,   \
                             ts_zzz_yyyzz,   \
                             ts_zzz_yyzzz,   \
                             ts_zzz_yzzzz,   \
                             ts_zzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzz_xxxxx[i] =
            -4.0 * ts_xxz_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxxx[i] * fe_0 + tk_xxzz_xxxxx[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxxx[i] * fz_0;

        tk_xxzzz_xxxxy[i] =
            -4.0 * ts_xxz_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxxy[i] * fe_0 + tk_xxzz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxxy[i] * fz_0;

        tk_xxzzz_xxxxz[i] = -2.0 * ts_zzz_xxxxz[i] * fbe_0 * fz_0 + tk_zzz_xxxxz[i] * fe_0 + 4.0 * tk_xzzz_xxxz[i] * fe_0 +
                            tk_xzzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxxz[i] * fz_0;

        tk_xxzzz_xxxyy[i] =
            -4.0 * ts_xxz_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxyy[i] * fe_0 + tk_xxzz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxyy[i] * fz_0;

        tk_xxzzz_xxxyz[i] = -2.0 * ts_zzz_xxxyz[i] * fbe_0 * fz_0 + tk_zzz_xxxyz[i] * fe_0 + 3.0 * tk_xzzz_xxyz[i] * fe_0 +
                            tk_xzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxyz[i] * fz_0;

        tk_xxzzz_xxxzz[i] = -2.0 * ts_zzz_xxxzz[i] * fbe_0 * fz_0 + tk_zzz_xxxzz[i] * fe_0 + 3.0 * tk_xzzz_xxzz[i] * fe_0 +
                            tk_xzzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxzz[i] * fz_0;

        tk_xxzzz_xxyyy[i] =
            -4.0 * ts_xxz_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxyyy[i] * fe_0 + tk_xxzz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxyyy[i] * fz_0;

        tk_xxzzz_xxyyz[i] = -2.0 * ts_zzz_xxyyz[i] * fbe_0 * fz_0 + tk_zzz_xxyyz[i] * fe_0 + 2.0 * tk_xzzz_xyyz[i] * fe_0 +
                            tk_xzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxyyz[i] * fz_0;

        tk_xxzzz_xxyzz[i] = -2.0 * ts_zzz_xxyzz[i] * fbe_0 * fz_0 + tk_zzz_xxyzz[i] * fe_0 + 2.0 * tk_xzzz_xyzz[i] * fe_0 +
                            tk_xzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxyzz[i] * fz_0;

        tk_xxzzz_xxzzz[i] = -2.0 * ts_zzz_xxzzz[i] * fbe_0 * fz_0 + tk_zzz_xxzzz[i] * fe_0 + 2.0 * tk_xzzz_xzzz[i] * fe_0 +
                            tk_xzzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxzzz[i] * fz_0;

        tk_xxzzz_xyyyy[i] =
            -4.0 * ts_xxz_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xyyyy[i] * fe_0 + tk_xxzz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xyyyy[i] * fz_0;

        tk_xxzzz_xyyyz[i] = -2.0 * ts_zzz_xyyyz[i] * fbe_0 * fz_0 + tk_zzz_xyyyz[i] * fe_0 + tk_xzzz_yyyz[i] * fe_0 + tk_xzzz_xyyyz[i] * pa_x[i] +
                            2.0 * ts_xxzzz_xyyyz[i] * fz_0;

        tk_xxzzz_xyyzz[i] = -2.0 * ts_zzz_xyyzz[i] * fbe_0 * fz_0 + tk_zzz_xyyzz[i] * fe_0 + tk_xzzz_yyzz[i] * fe_0 + tk_xzzz_xyyzz[i] * pa_x[i] +
                            2.0 * ts_xxzzz_xyyzz[i] * fz_0;

        tk_xxzzz_xyzzz[i] = -2.0 * ts_zzz_xyzzz[i] * fbe_0 * fz_0 + tk_zzz_xyzzz[i] * fe_0 + tk_xzzz_yzzz[i] * fe_0 + tk_xzzz_xyzzz[i] * pa_x[i] +
                            2.0 * ts_xxzzz_xyzzz[i] * fz_0;

        tk_xxzzz_xzzzz[i] = -2.0 * ts_zzz_xzzzz[i] * fbe_0 * fz_0 + tk_zzz_xzzzz[i] * fe_0 + tk_xzzz_zzzz[i] * fe_0 + tk_xzzz_xzzzz[i] * pa_x[i] +
                            2.0 * ts_xxzzz_xzzzz[i] * fz_0;

        tk_xxzzz_yyyyy[i] =
            -2.0 * ts_zzz_yyyyy[i] * fbe_0 * fz_0 + tk_zzz_yyyyy[i] * fe_0 + tk_xzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyyy[i] * fz_0;

        tk_xxzzz_yyyyz[i] =
            -2.0 * ts_zzz_yyyyz[i] * fbe_0 * fz_0 + tk_zzz_yyyyz[i] * fe_0 + tk_xzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyyz[i] * fz_0;

        tk_xxzzz_yyyzz[i] =
            -2.0 * ts_zzz_yyyzz[i] * fbe_0 * fz_0 + tk_zzz_yyyzz[i] * fe_0 + tk_xzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyzz[i] * fz_0;

        tk_xxzzz_yyzzz[i] =
            -2.0 * ts_zzz_yyzzz[i] * fbe_0 * fz_0 + tk_zzz_yyzzz[i] * fe_0 + tk_xzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyzzz[i] * fz_0;

        tk_xxzzz_yzzzz[i] =
            -2.0 * ts_zzz_yzzzz[i] * fbe_0 * fz_0 + tk_zzz_yzzzz[i] * fe_0 + tk_xzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yzzzz[i] * fz_0;

        tk_xxzzz_zzzzz[i] =
            -2.0 * ts_zzz_zzzzz[i] * fbe_0 * fz_0 + tk_zzz_zzzzz[i] * fe_0 + tk_xzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_zzzzz[i] * fz_0;
    }

    // Set up 210-231 components of targeted buffer : HH

    auto tk_xyyyy_xxxxx = pbuffer.data(idx_kin_hh + 210);

    auto tk_xyyyy_xxxxy = pbuffer.data(idx_kin_hh + 211);

    auto tk_xyyyy_xxxxz = pbuffer.data(idx_kin_hh + 212);

    auto tk_xyyyy_xxxyy = pbuffer.data(idx_kin_hh + 213);

    auto tk_xyyyy_xxxyz = pbuffer.data(idx_kin_hh + 214);

    auto tk_xyyyy_xxxzz = pbuffer.data(idx_kin_hh + 215);

    auto tk_xyyyy_xxyyy = pbuffer.data(idx_kin_hh + 216);

    auto tk_xyyyy_xxyyz = pbuffer.data(idx_kin_hh + 217);

    auto tk_xyyyy_xxyzz = pbuffer.data(idx_kin_hh + 218);

    auto tk_xyyyy_xxzzz = pbuffer.data(idx_kin_hh + 219);

    auto tk_xyyyy_xyyyy = pbuffer.data(idx_kin_hh + 220);

    auto tk_xyyyy_xyyyz = pbuffer.data(idx_kin_hh + 221);

    auto tk_xyyyy_xyyzz = pbuffer.data(idx_kin_hh + 222);

    auto tk_xyyyy_xyzzz = pbuffer.data(idx_kin_hh + 223);

    auto tk_xyyyy_xzzzz = pbuffer.data(idx_kin_hh + 224);

    auto tk_xyyyy_yyyyy = pbuffer.data(idx_kin_hh + 225);

    auto tk_xyyyy_yyyyz = pbuffer.data(idx_kin_hh + 226);

    auto tk_xyyyy_yyyzz = pbuffer.data(idx_kin_hh + 227);

    auto tk_xyyyy_yyzzz = pbuffer.data(idx_kin_hh + 228);

    auto tk_xyyyy_yzzzz = pbuffer.data(idx_kin_hh + 229);

    auto tk_xyyyy_zzzzz = pbuffer.data(idx_kin_hh + 230);

#pragma omp simd aligned(pa_x,               \
                             tk_xyyyy_xxxxx, \
                             tk_xyyyy_xxxxy, \
                             tk_xyyyy_xxxxz, \
                             tk_xyyyy_xxxyy, \
                             tk_xyyyy_xxxyz, \
                             tk_xyyyy_xxxzz, \
                             tk_xyyyy_xxyyy, \
                             tk_xyyyy_xxyyz, \
                             tk_xyyyy_xxyzz, \
                             tk_xyyyy_xxzzz, \
                             tk_xyyyy_xyyyy, \
                             tk_xyyyy_xyyyz, \
                             tk_xyyyy_xyyzz, \
                             tk_xyyyy_xyzzz, \
                             tk_xyyyy_xzzzz, \
                             tk_xyyyy_yyyyy, \
                             tk_xyyyy_yyyyz, \
                             tk_xyyyy_yyyzz, \
                             tk_xyyyy_yyzzz, \
                             tk_xyyyy_yzzzz, \
                             tk_xyyyy_zzzzz, \
                             tk_yyyy_xxxx,   \
                             tk_yyyy_xxxxx,  \
                             tk_yyyy_xxxxy,  \
                             tk_yyyy_xxxxz,  \
                             tk_yyyy_xxxy,   \
                             tk_yyyy_xxxyy,  \
                             tk_yyyy_xxxyz,  \
                             tk_yyyy_xxxz,   \
                             tk_yyyy_xxxzz,  \
                             tk_yyyy_xxyy,   \
                             tk_yyyy_xxyyy,  \
                             tk_yyyy_xxyyz,  \
                             tk_yyyy_xxyz,   \
                             tk_yyyy_xxyzz,  \
                             tk_yyyy_xxzz,   \
                             tk_yyyy_xxzzz,  \
                             tk_yyyy_xyyy,   \
                             tk_yyyy_xyyyy,  \
                             tk_yyyy_xyyyz,  \
                             tk_yyyy_xyyz,   \
                             tk_yyyy_xyyzz,  \
                             tk_yyyy_xyzz,   \
                             tk_yyyy_xyzzz,  \
                             tk_yyyy_xzzz,   \
                             tk_yyyy_xzzzz,  \
                             tk_yyyy_yyyy,   \
                             tk_yyyy_yyyyy,  \
                             tk_yyyy_yyyyz,  \
                             tk_yyyy_yyyz,   \
                             tk_yyyy_yyyzz,  \
                             tk_yyyy_yyzz,   \
                             tk_yyyy_yyzzz,  \
                             tk_yyyy_yzzz,   \
                             tk_yyyy_yzzzz,  \
                             tk_yyyy_zzzz,   \
                             tk_yyyy_zzzzz,  \
                             ts_xyyyy_xxxxx, \
                             ts_xyyyy_xxxxy, \
                             ts_xyyyy_xxxxz, \
                             ts_xyyyy_xxxyy, \
                             ts_xyyyy_xxxyz, \
                             ts_xyyyy_xxxzz, \
                             ts_xyyyy_xxyyy, \
                             ts_xyyyy_xxyyz, \
                             ts_xyyyy_xxyzz, \
                             ts_xyyyy_xxzzz, \
                             ts_xyyyy_xyyyy, \
                             ts_xyyyy_xyyyz, \
                             ts_xyyyy_xyyzz, \
                             ts_xyyyy_xyzzz, \
                             ts_xyyyy_xzzzz, \
                             ts_xyyyy_yyyyy, \
                             ts_xyyyy_yyyyz, \
                             ts_xyyyy_yyyzz, \
                             ts_xyyyy_yyzzz, \
                             ts_xyyyy_yzzzz, \
                             ts_xyyyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyy_xxxxx[i] = 5.0 * tk_yyyy_xxxx[i] * fe_0 + tk_yyyy_xxxxx[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxx[i] * fz_0;

        tk_xyyyy_xxxxy[i] = 4.0 * tk_yyyy_xxxy[i] * fe_0 + tk_yyyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxy[i] * fz_0;

        tk_xyyyy_xxxxz[i] = 4.0 * tk_yyyy_xxxz[i] * fe_0 + tk_yyyy_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxz[i] * fz_0;

        tk_xyyyy_xxxyy[i] = 3.0 * tk_yyyy_xxyy[i] * fe_0 + tk_yyyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxyy[i] * fz_0;

        tk_xyyyy_xxxyz[i] = 3.0 * tk_yyyy_xxyz[i] * fe_0 + tk_yyyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxyz[i] * fz_0;

        tk_xyyyy_xxxzz[i] = 3.0 * tk_yyyy_xxzz[i] * fe_0 + tk_yyyy_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxzz[i] * fz_0;

        tk_xyyyy_xxyyy[i] = 2.0 * tk_yyyy_xyyy[i] * fe_0 + tk_yyyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyyy[i] * fz_0;

        tk_xyyyy_xxyyz[i] = 2.0 * tk_yyyy_xyyz[i] * fe_0 + tk_yyyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyyz[i] * fz_0;

        tk_xyyyy_xxyzz[i] = 2.0 * tk_yyyy_xyzz[i] * fe_0 + tk_yyyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyzz[i] * fz_0;

        tk_xyyyy_xxzzz[i] = 2.0 * tk_yyyy_xzzz[i] * fe_0 + tk_yyyy_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxzzz[i] * fz_0;

        tk_xyyyy_xyyyy[i] = tk_yyyy_yyyy[i] * fe_0 + tk_yyyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyyy[i] * fz_0;

        tk_xyyyy_xyyyz[i] = tk_yyyy_yyyz[i] * fe_0 + tk_yyyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyyz[i] * fz_0;

        tk_xyyyy_xyyzz[i] = tk_yyyy_yyzz[i] * fe_0 + tk_yyyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyzz[i] * fz_0;

        tk_xyyyy_xyzzz[i] = tk_yyyy_yzzz[i] * fe_0 + tk_yyyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyzzz[i] * fz_0;

        tk_xyyyy_xzzzz[i] = tk_yyyy_zzzz[i] * fe_0 + tk_yyyy_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xzzzz[i] * fz_0;

        tk_xyyyy_yyyyy[i] = tk_yyyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyyy[i] * fz_0;

        tk_xyyyy_yyyyz[i] = tk_yyyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyyz[i] * fz_0;

        tk_xyyyy_yyyzz[i] = tk_yyyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyzz[i] * fz_0;

        tk_xyyyy_yyzzz[i] = tk_yyyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyzzz[i] * fz_0;

        tk_xyyyy_yzzzz[i] = tk_yyyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yzzzz[i] * fz_0;

        tk_xyyyy_zzzzz[i] = tk_yyyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_zzzzz[i] * fz_0;
    }

    // Set up 231-252 components of targeted buffer : HH

    auto tk_xyyyz_xxxxx = pbuffer.data(idx_kin_hh + 231);

    auto tk_xyyyz_xxxxy = pbuffer.data(idx_kin_hh + 232);

    auto tk_xyyyz_xxxxz = pbuffer.data(idx_kin_hh + 233);

    auto tk_xyyyz_xxxyy = pbuffer.data(idx_kin_hh + 234);

    auto tk_xyyyz_xxxyz = pbuffer.data(idx_kin_hh + 235);

    auto tk_xyyyz_xxxzz = pbuffer.data(idx_kin_hh + 236);

    auto tk_xyyyz_xxyyy = pbuffer.data(idx_kin_hh + 237);

    auto tk_xyyyz_xxyyz = pbuffer.data(idx_kin_hh + 238);

    auto tk_xyyyz_xxyzz = pbuffer.data(idx_kin_hh + 239);

    auto tk_xyyyz_xxzzz = pbuffer.data(idx_kin_hh + 240);

    auto tk_xyyyz_xyyyy = pbuffer.data(idx_kin_hh + 241);

    auto tk_xyyyz_xyyyz = pbuffer.data(idx_kin_hh + 242);

    auto tk_xyyyz_xyyzz = pbuffer.data(idx_kin_hh + 243);

    auto tk_xyyyz_xyzzz = pbuffer.data(idx_kin_hh + 244);

    auto tk_xyyyz_xzzzz = pbuffer.data(idx_kin_hh + 245);

    auto tk_xyyyz_yyyyy = pbuffer.data(idx_kin_hh + 246);

    auto tk_xyyyz_yyyyz = pbuffer.data(idx_kin_hh + 247);

    auto tk_xyyyz_yyyzz = pbuffer.data(idx_kin_hh + 248);

    auto tk_xyyyz_yyzzz = pbuffer.data(idx_kin_hh + 249);

    auto tk_xyyyz_yzzzz = pbuffer.data(idx_kin_hh + 250);

    auto tk_xyyyz_zzzzz = pbuffer.data(idx_kin_hh + 251);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xyyy_xxxxx,  \
                             tk_xyyy_xxxxy,  \
                             tk_xyyy_xxxyy,  \
                             tk_xyyy_xxyyy,  \
                             tk_xyyy_xyyyy,  \
                             tk_xyyyz_xxxxx, \
                             tk_xyyyz_xxxxy, \
                             tk_xyyyz_xxxxz, \
                             tk_xyyyz_xxxyy, \
                             tk_xyyyz_xxxyz, \
                             tk_xyyyz_xxxzz, \
                             tk_xyyyz_xxyyy, \
                             tk_xyyyz_xxyyz, \
                             tk_xyyyz_xxyzz, \
                             tk_xyyyz_xxzzz, \
                             tk_xyyyz_xyyyy, \
                             tk_xyyyz_xyyyz, \
                             tk_xyyyz_xyyzz, \
                             tk_xyyyz_xyzzz, \
                             tk_xyyyz_xzzzz, \
                             tk_xyyyz_yyyyy, \
                             tk_xyyyz_yyyyz, \
                             tk_xyyyz_yyyzz, \
                             tk_xyyyz_yyzzz, \
                             tk_xyyyz_yzzzz, \
                             tk_xyyyz_zzzzz, \
                             tk_yyyz_xxxxz,  \
                             tk_yyyz_xxxyz,  \
                             tk_yyyz_xxxz,   \
                             tk_yyyz_xxxzz,  \
                             tk_yyyz_xxyyz,  \
                             tk_yyyz_xxyz,   \
                             tk_yyyz_xxyzz,  \
                             tk_yyyz_xxzz,   \
                             tk_yyyz_xxzzz,  \
                             tk_yyyz_xyyyz,  \
                             tk_yyyz_xyyz,   \
                             tk_yyyz_xyyzz,  \
                             tk_yyyz_xyzz,   \
                             tk_yyyz_xyzzz,  \
                             tk_yyyz_xzzz,   \
                             tk_yyyz_xzzzz,  \
                             tk_yyyz_yyyyy,  \
                             tk_yyyz_yyyyz,  \
                             tk_yyyz_yyyz,   \
                             tk_yyyz_yyyzz,  \
                             tk_yyyz_yyzz,   \
                             tk_yyyz_yyzzz,  \
                             tk_yyyz_yzzz,   \
                             tk_yyyz_yzzzz,  \
                             tk_yyyz_zzzz,   \
                             tk_yyyz_zzzzz,  \
                             ts_xyyyz_xxxxx, \
                             ts_xyyyz_xxxxy, \
                             ts_xyyyz_xxxxz, \
                             ts_xyyyz_xxxyy, \
                             ts_xyyyz_xxxyz, \
                             ts_xyyyz_xxxzz, \
                             ts_xyyyz_xxyyy, \
                             ts_xyyyz_xxyyz, \
                             ts_xyyyz_xxyzz, \
                             ts_xyyyz_xxzzz, \
                             ts_xyyyz_xyyyy, \
                             ts_xyyyz_xyyyz, \
                             ts_xyyyz_xyyzz, \
                             ts_xyyyz_xyzzz, \
                             ts_xyyyz_xzzzz, \
                             ts_xyyyz_yyyyy, \
                             ts_xyyyz_yyyyz, \
                             ts_xyyyz_yyyzz, \
                             ts_xyyyz_yyzzz, \
                             ts_xyyyz_yzzzz, \
                             ts_xyyyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyz_xxxxx[i] = tk_xyyy_xxxxx[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxxx[i] * fz_0;

        tk_xyyyz_xxxxy[i] = tk_xyyy_xxxxy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxxy[i] * fz_0;

        tk_xyyyz_xxxxz[i] = 4.0 * tk_yyyz_xxxz[i] * fe_0 + tk_yyyz_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxxz[i] * fz_0;

        tk_xyyyz_xxxyy[i] = tk_xyyy_xxxyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxyy[i] * fz_0;

        tk_xyyyz_xxxyz[i] = 3.0 * tk_yyyz_xxyz[i] * fe_0 + tk_yyyz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxyz[i] * fz_0;

        tk_xyyyz_xxxzz[i] = 3.0 * tk_yyyz_xxzz[i] * fe_0 + tk_yyyz_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxzz[i] * fz_0;

        tk_xyyyz_xxyyy[i] = tk_xyyy_xxyyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxyyy[i] * fz_0;

        tk_xyyyz_xxyyz[i] = 2.0 * tk_yyyz_xyyz[i] * fe_0 + tk_yyyz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxyyz[i] * fz_0;

        tk_xyyyz_xxyzz[i] = 2.0 * tk_yyyz_xyzz[i] * fe_0 + tk_yyyz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxyzz[i] * fz_0;

        tk_xyyyz_xxzzz[i] = 2.0 * tk_yyyz_xzzz[i] * fe_0 + tk_yyyz_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxzzz[i] * fz_0;

        tk_xyyyz_xyyyy[i] = tk_xyyy_xyyyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xyyyy[i] * fz_0;

        tk_xyyyz_xyyyz[i] = tk_yyyz_yyyz[i] * fe_0 + tk_yyyz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyyyz[i] * fz_0;

        tk_xyyyz_xyyzz[i] = tk_yyyz_yyzz[i] * fe_0 + tk_yyyz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyyzz[i] * fz_0;

        tk_xyyyz_xyzzz[i] = tk_yyyz_yzzz[i] * fe_0 + tk_yyyz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyzzz[i] * fz_0;

        tk_xyyyz_xzzzz[i] = tk_yyyz_zzzz[i] * fe_0 + tk_yyyz_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xzzzz[i] * fz_0;

        tk_xyyyz_yyyyy[i] = tk_yyyz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyyy[i] * fz_0;

        tk_xyyyz_yyyyz[i] = tk_yyyz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyyz[i] * fz_0;

        tk_xyyyz_yyyzz[i] = tk_yyyz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyzz[i] * fz_0;

        tk_xyyyz_yyzzz[i] = tk_yyyz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyzzz[i] * fz_0;

        tk_xyyyz_yzzzz[i] = tk_yyyz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yzzzz[i] * fz_0;

        tk_xyyyz_zzzzz[i] = tk_yyyz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_zzzzz[i] * fz_0;
    }

    // Set up 252-273 components of targeted buffer : HH

    auto tk_xyyzz_xxxxx = pbuffer.data(idx_kin_hh + 252);

    auto tk_xyyzz_xxxxy = pbuffer.data(idx_kin_hh + 253);

    auto tk_xyyzz_xxxxz = pbuffer.data(idx_kin_hh + 254);

    auto tk_xyyzz_xxxyy = pbuffer.data(idx_kin_hh + 255);

    auto tk_xyyzz_xxxyz = pbuffer.data(idx_kin_hh + 256);

    auto tk_xyyzz_xxxzz = pbuffer.data(idx_kin_hh + 257);

    auto tk_xyyzz_xxyyy = pbuffer.data(idx_kin_hh + 258);

    auto tk_xyyzz_xxyyz = pbuffer.data(idx_kin_hh + 259);

    auto tk_xyyzz_xxyzz = pbuffer.data(idx_kin_hh + 260);

    auto tk_xyyzz_xxzzz = pbuffer.data(idx_kin_hh + 261);

    auto tk_xyyzz_xyyyy = pbuffer.data(idx_kin_hh + 262);

    auto tk_xyyzz_xyyyz = pbuffer.data(idx_kin_hh + 263);

    auto tk_xyyzz_xyyzz = pbuffer.data(idx_kin_hh + 264);

    auto tk_xyyzz_xyzzz = pbuffer.data(idx_kin_hh + 265);

    auto tk_xyyzz_xzzzz = pbuffer.data(idx_kin_hh + 266);

    auto tk_xyyzz_yyyyy = pbuffer.data(idx_kin_hh + 267);

    auto tk_xyyzz_yyyyz = pbuffer.data(idx_kin_hh + 268);

    auto tk_xyyzz_yyyzz = pbuffer.data(idx_kin_hh + 269);

    auto tk_xyyzz_yyzzz = pbuffer.data(idx_kin_hh + 270);

    auto tk_xyyzz_yzzzz = pbuffer.data(idx_kin_hh + 271);

    auto tk_xyyzz_zzzzz = pbuffer.data(idx_kin_hh + 272);

#pragma omp simd aligned(pa_x,               \
                             tk_xyyzz_xxxxx, \
                             tk_xyyzz_xxxxy, \
                             tk_xyyzz_xxxxz, \
                             tk_xyyzz_xxxyy, \
                             tk_xyyzz_xxxyz, \
                             tk_xyyzz_xxxzz, \
                             tk_xyyzz_xxyyy, \
                             tk_xyyzz_xxyyz, \
                             tk_xyyzz_xxyzz, \
                             tk_xyyzz_xxzzz, \
                             tk_xyyzz_xyyyy, \
                             tk_xyyzz_xyyyz, \
                             tk_xyyzz_xyyzz, \
                             tk_xyyzz_xyzzz, \
                             tk_xyyzz_xzzzz, \
                             tk_xyyzz_yyyyy, \
                             tk_xyyzz_yyyyz, \
                             tk_xyyzz_yyyzz, \
                             tk_xyyzz_yyzzz, \
                             tk_xyyzz_yzzzz, \
                             tk_xyyzz_zzzzz, \
                             tk_yyzz_xxxx,   \
                             tk_yyzz_xxxxx,  \
                             tk_yyzz_xxxxy,  \
                             tk_yyzz_xxxxz,  \
                             tk_yyzz_xxxy,   \
                             tk_yyzz_xxxyy,  \
                             tk_yyzz_xxxyz,  \
                             tk_yyzz_xxxz,   \
                             tk_yyzz_xxxzz,  \
                             tk_yyzz_xxyy,   \
                             tk_yyzz_xxyyy,  \
                             tk_yyzz_xxyyz,  \
                             tk_yyzz_xxyz,   \
                             tk_yyzz_xxyzz,  \
                             tk_yyzz_xxzz,   \
                             tk_yyzz_xxzzz,  \
                             tk_yyzz_xyyy,   \
                             tk_yyzz_xyyyy,  \
                             tk_yyzz_xyyyz,  \
                             tk_yyzz_xyyz,   \
                             tk_yyzz_xyyzz,  \
                             tk_yyzz_xyzz,   \
                             tk_yyzz_xyzzz,  \
                             tk_yyzz_xzzz,   \
                             tk_yyzz_xzzzz,  \
                             tk_yyzz_yyyy,   \
                             tk_yyzz_yyyyy,  \
                             tk_yyzz_yyyyz,  \
                             tk_yyzz_yyyz,   \
                             tk_yyzz_yyyzz,  \
                             tk_yyzz_yyzz,   \
                             tk_yyzz_yyzzz,  \
                             tk_yyzz_yzzz,   \
                             tk_yyzz_yzzzz,  \
                             tk_yyzz_zzzz,   \
                             tk_yyzz_zzzzz,  \
                             ts_xyyzz_xxxxx, \
                             ts_xyyzz_xxxxy, \
                             ts_xyyzz_xxxxz, \
                             ts_xyyzz_xxxyy, \
                             ts_xyyzz_xxxyz, \
                             ts_xyyzz_xxxzz, \
                             ts_xyyzz_xxyyy, \
                             ts_xyyzz_xxyyz, \
                             ts_xyyzz_xxyzz, \
                             ts_xyyzz_xxzzz, \
                             ts_xyyzz_xyyyy, \
                             ts_xyyzz_xyyyz, \
                             ts_xyyzz_xyyzz, \
                             ts_xyyzz_xyzzz, \
                             ts_xyyzz_xzzzz, \
                             ts_xyyzz_yyyyy, \
                             ts_xyyzz_yyyyz, \
                             ts_xyyzz_yyyzz, \
                             ts_xyyzz_yyzzz, \
                             ts_xyyzz_yzzzz, \
                             ts_xyyzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzz_xxxxx[i] = 5.0 * tk_yyzz_xxxx[i] * fe_0 + tk_yyzz_xxxxx[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxx[i] * fz_0;

        tk_xyyzz_xxxxy[i] = 4.0 * tk_yyzz_xxxy[i] * fe_0 + tk_yyzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxy[i] * fz_0;

        tk_xyyzz_xxxxz[i] = 4.0 * tk_yyzz_xxxz[i] * fe_0 + tk_yyzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxz[i] * fz_0;

        tk_xyyzz_xxxyy[i] = 3.0 * tk_yyzz_xxyy[i] * fe_0 + tk_yyzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxyy[i] * fz_0;

        tk_xyyzz_xxxyz[i] = 3.0 * tk_yyzz_xxyz[i] * fe_0 + tk_yyzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxyz[i] * fz_0;

        tk_xyyzz_xxxzz[i] = 3.0 * tk_yyzz_xxzz[i] * fe_0 + tk_yyzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxzz[i] * fz_0;

        tk_xyyzz_xxyyy[i] = 2.0 * tk_yyzz_xyyy[i] * fe_0 + tk_yyzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyyy[i] * fz_0;

        tk_xyyzz_xxyyz[i] = 2.0 * tk_yyzz_xyyz[i] * fe_0 + tk_yyzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyyz[i] * fz_0;

        tk_xyyzz_xxyzz[i] = 2.0 * tk_yyzz_xyzz[i] * fe_0 + tk_yyzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyzz[i] * fz_0;

        tk_xyyzz_xxzzz[i] = 2.0 * tk_yyzz_xzzz[i] * fe_0 + tk_yyzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxzzz[i] * fz_0;

        tk_xyyzz_xyyyy[i] = tk_yyzz_yyyy[i] * fe_0 + tk_yyzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyyy[i] * fz_0;

        tk_xyyzz_xyyyz[i] = tk_yyzz_yyyz[i] * fe_0 + tk_yyzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyyz[i] * fz_0;

        tk_xyyzz_xyyzz[i] = tk_yyzz_yyzz[i] * fe_0 + tk_yyzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyzz[i] * fz_0;

        tk_xyyzz_xyzzz[i] = tk_yyzz_yzzz[i] * fe_0 + tk_yyzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyzzz[i] * fz_0;

        tk_xyyzz_xzzzz[i] = tk_yyzz_zzzz[i] * fe_0 + tk_yyzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xzzzz[i] * fz_0;

        tk_xyyzz_yyyyy[i] = tk_yyzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyyy[i] * fz_0;

        tk_xyyzz_yyyyz[i] = tk_yyzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyyz[i] * fz_0;

        tk_xyyzz_yyyzz[i] = tk_yyzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyzz[i] * fz_0;

        tk_xyyzz_yyzzz[i] = tk_yyzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyzzz[i] * fz_0;

        tk_xyyzz_yzzzz[i] = tk_yyzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yzzzz[i] * fz_0;

        tk_xyyzz_zzzzz[i] = tk_yyzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_zzzzz[i] * fz_0;
    }

    // Set up 273-294 components of targeted buffer : HH

    auto tk_xyzzz_xxxxx = pbuffer.data(idx_kin_hh + 273);

    auto tk_xyzzz_xxxxy = pbuffer.data(idx_kin_hh + 274);

    auto tk_xyzzz_xxxxz = pbuffer.data(idx_kin_hh + 275);

    auto tk_xyzzz_xxxyy = pbuffer.data(idx_kin_hh + 276);

    auto tk_xyzzz_xxxyz = pbuffer.data(idx_kin_hh + 277);

    auto tk_xyzzz_xxxzz = pbuffer.data(idx_kin_hh + 278);

    auto tk_xyzzz_xxyyy = pbuffer.data(idx_kin_hh + 279);

    auto tk_xyzzz_xxyyz = pbuffer.data(idx_kin_hh + 280);

    auto tk_xyzzz_xxyzz = pbuffer.data(idx_kin_hh + 281);

    auto tk_xyzzz_xxzzz = pbuffer.data(idx_kin_hh + 282);

    auto tk_xyzzz_xyyyy = pbuffer.data(idx_kin_hh + 283);

    auto tk_xyzzz_xyyyz = pbuffer.data(idx_kin_hh + 284);

    auto tk_xyzzz_xyyzz = pbuffer.data(idx_kin_hh + 285);

    auto tk_xyzzz_xyzzz = pbuffer.data(idx_kin_hh + 286);

    auto tk_xyzzz_xzzzz = pbuffer.data(idx_kin_hh + 287);

    auto tk_xyzzz_yyyyy = pbuffer.data(idx_kin_hh + 288);

    auto tk_xyzzz_yyyyz = pbuffer.data(idx_kin_hh + 289);

    auto tk_xyzzz_yyyzz = pbuffer.data(idx_kin_hh + 290);

    auto tk_xyzzz_yyzzz = pbuffer.data(idx_kin_hh + 291);

    auto tk_xyzzz_yzzzz = pbuffer.data(idx_kin_hh + 292);

    auto tk_xyzzz_zzzzz = pbuffer.data(idx_kin_hh + 293);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xyzzz_xxxxx, \
                             tk_xyzzz_xxxxy, \
                             tk_xyzzz_xxxxz, \
                             tk_xyzzz_xxxyy, \
                             tk_xyzzz_xxxyz, \
                             tk_xyzzz_xxxzz, \
                             tk_xyzzz_xxyyy, \
                             tk_xyzzz_xxyyz, \
                             tk_xyzzz_xxyzz, \
                             tk_xyzzz_xxzzz, \
                             tk_xyzzz_xyyyy, \
                             tk_xyzzz_xyyyz, \
                             tk_xyzzz_xyyzz, \
                             tk_xyzzz_xyzzz, \
                             tk_xyzzz_xzzzz, \
                             tk_xyzzz_yyyyy, \
                             tk_xyzzz_yyyyz, \
                             tk_xyzzz_yyyzz, \
                             tk_xyzzz_yyzzz, \
                             tk_xyzzz_yzzzz, \
                             tk_xyzzz_zzzzz, \
                             tk_xzzz_xxxxx,  \
                             tk_xzzz_xxxxz,  \
                             tk_xzzz_xxxzz,  \
                             tk_xzzz_xxzzz,  \
                             tk_xzzz_xzzzz,  \
                             tk_yzzz_xxxxy,  \
                             tk_yzzz_xxxy,   \
                             tk_yzzz_xxxyy,  \
                             tk_yzzz_xxxyz,  \
                             tk_yzzz_xxyy,   \
                             tk_yzzz_xxyyy,  \
                             tk_yzzz_xxyyz,  \
                             tk_yzzz_xxyz,   \
                             tk_yzzz_xxyzz,  \
                             tk_yzzz_xyyy,   \
                             tk_yzzz_xyyyy,  \
                             tk_yzzz_xyyyz,  \
                             tk_yzzz_xyyz,   \
                             tk_yzzz_xyyzz,  \
                             tk_yzzz_xyzz,   \
                             tk_yzzz_xyzzz,  \
                             tk_yzzz_yyyy,   \
                             tk_yzzz_yyyyy,  \
                             tk_yzzz_yyyyz,  \
                             tk_yzzz_yyyz,   \
                             tk_yzzz_yyyzz,  \
                             tk_yzzz_yyzz,   \
                             tk_yzzz_yyzzz,  \
                             tk_yzzz_yzzz,   \
                             tk_yzzz_yzzzz,  \
                             tk_yzzz_zzzzz,  \
                             ts_xyzzz_xxxxx, \
                             ts_xyzzz_xxxxy, \
                             ts_xyzzz_xxxxz, \
                             ts_xyzzz_xxxyy, \
                             ts_xyzzz_xxxyz, \
                             ts_xyzzz_xxxzz, \
                             ts_xyzzz_xxyyy, \
                             ts_xyzzz_xxyyz, \
                             ts_xyzzz_xxyzz, \
                             ts_xyzzz_xxzzz, \
                             ts_xyzzz_xyyyy, \
                             ts_xyzzz_xyyyz, \
                             ts_xyzzz_xyyzz, \
                             ts_xyzzz_xyzzz, \
                             ts_xyzzz_xzzzz, \
                             ts_xyzzz_yyyyy, \
                             ts_xyzzz_yyyyz, \
                             ts_xyzzz_yyyzz, \
                             ts_xyzzz_yyzzz, \
                             ts_xyzzz_yzzzz, \
                             ts_xyzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzz_xxxxx[i] = tk_xzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxxx[i] * fz_0;

        tk_xyzzz_xxxxy[i] = 4.0 * tk_yzzz_xxxy[i] * fe_0 + tk_yzzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxxy[i] * fz_0;

        tk_xyzzz_xxxxz[i] = tk_xzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxxz[i] * fz_0;

        tk_xyzzz_xxxyy[i] = 3.0 * tk_yzzz_xxyy[i] * fe_0 + tk_yzzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxyy[i] * fz_0;

        tk_xyzzz_xxxyz[i] = 3.0 * tk_yzzz_xxyz[i] * fe_0 + tk_yzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxyz[i] * fz_0;

        tk_xyzzz_xxxzz[i] = tk_xzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxzz[i] * fz_0;

        tk_xyzzz_xxyyy[i] = 2.0 * tk_yzzz_xyyy[i] * fe_0 + tk_yzzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyyy[i] * fz_0;

        tk_xyzzz_xxyyz[i] = 2.0 * tk_yzzz_xyyz[i] * fe_0 + tk_yzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyyz[i] * fz_0;

        tk_xyzzz_xxyzz[i] = 2.0 * tk_yzzz_xyzz[i] * fe_0 + tk_yzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyzz[i] * fz_0;

        tk_xyzzz_xxzzz[i] = tk_xzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxzzz[i] * fz_0;

        tk_xyzzz_xyyyy[i] = tk_yzzz_yyyy[i] * fe_0 + tk_yzzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyyy[i] * fz_0;

        tk_xyzzz_xyyyz[i] = tk_yzzz_yyyz[i] * fe_0 + tk_yzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyyz[i] * fz_0;

        tk_xyzzz_xyyzz[i] = tk_yzzz_yyzz[i] * fe_0 + tk_yzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyzz[i] * fz_0;

        tk_xyzzz_xyzzz[i] = tk_yzzz_yzzz[i] * fe_0 + tk_yzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyzzz[i] * fz_0;

        tk_xyzzz_xzzzz[i] = tk_xzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xzzzz[i] * fz_0;

        tk_xyzzz_yyyyy[i] = tk_yzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyyy[i] * fz_0;

        tk_xyzzz_yyyyz[i] = tk_yzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyyz[i] * fz_0;

        tk_xyzzz_yyyzz[i] = tk_yzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyzz[i] * fz_0;

        tk_xyzzz_yyzzz[i] = tk_yzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyzzz[i] * fz_0;

        tk_xyzzz_yzzzz[i] = tk_yzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yzzzz[i] * fz_0;

        tk_xyzzz_zzzzz[i] = tk_yzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_zzzzz[i] * fz_0;
    }

    // Set up 294-315 components of targeted buffer : HH

    auto tk_xzzzz_xxxxx = pbuffer.data(idx_kin_hh + 294);

    auto tk_xzzzz_xxxxy = pbuffer.data(idx_kin_hh + 295);

    auto tk_xzzzz_xxxxz = pbuffer.data(idx_kin_hh + 296);

    auto tk_xzzzz_xxxyy = pbuffer.data(idx_kin_hh + 297);

    auto tk_xzzzz_xxxyz = pbuffer.data(idx_kin_hh + 298);

    auto tk_xzzzz_xxxzz = pbuffer.data(idx_kin_hh + 299);

    auto tk_xzzzz_xxyyy = pbuffer.data(idx_kin_hh + 300);

    auto tk_xzzzz_xxyyz = pbuffer.data(idx_kin_hh + 301);

    auto tk_xzzzz_xxyzz = pbuffer.data(idx_kin_hh + 302);

    auto tk_xzzzz_xxzzz = pbuffer.data(idx_kin_hh + 303);

    auto tk_xzzzz_xyyyy = pbuffer.data(idx_kin_hh + 304);

    auto tk_xzzzz_xyyyz = pbuffer.data(idx_kin_hh + 305);

    auto tk_xzzzz_xyyzz = pbuffer.data(idx_kin_hh + 306);

    auto tk_xzzzz_xyzzz = pbuffer.data(idx_kin_hh + 307);

    auto tk_xzzzz_xzzzz = pbuffer.data(idx_kin_hh + 308);

    auto tk_xzzzz_yyyyy = pbuffer.data(idx_kin_hh + 309);

    auto tk_xzzzz_yyyyz = pbuffer.data(idx_kin_hh + 310);

    auto tk_xzzzz_yyyzz = pbuffer.data(idx_kin_hh + 311);

    auto tk_xzzzz_yyzzz = pbuffer.data(idx_kin_hh + 312);

    auto tk_xzzzz_yzzzz = pbuffer.data(idx_kin_hh + 313);

    auto tk_xzzzz_zzzzz = pbuffer.data(idx_kin_hh + 314);

#pragma omp simd aligned(pa_x,               \
                             tk_xzzzz_xxxxx, \
                             tk_xzzzz_xxxxy, \
                             tk_xzzzz_xxxxz, \
                             tk_xzzzz_xxxyy, \
                             tk_xzzzz_xxxyz, \
                             tk_xzzzz_xxxzz, \
                             tk_xzzzz_xxyyy, \
                             tk_xzzzz_xxyyz, \
                             tk_xzzzz_xxyzz, \
                             tk_xzzzz_xxzzz, \
                             tk_xzzzz_xyyyy, \
                             tk_xzzzz_xyyyz, \
                             tk_xzzzz_xyyzz, \
                             tk_xzzzz_xyzzz, \
                             tk_xzzzz_xzzzz, \
                             tk_xzzzz_yyyyy, \
                             tk_xzzzz_yyyyz, \
                             tk_xzzzz_yyyzz, \
                             tk_xzzzz_yyzzz, \
                             tk_xzzzz_yzzzz, \
                             tk_xzzzz_zzzzz, \
                             tk_zzzz_xxxx,   \
                             tk_zzzz_xxxxx,  \
                             tk_zzzz_xxxxy,  \
                             tk_zzzz_xxxxz,  \
                             tk_zzzz_xxxy,   \
                             tk_zzzz_xxxyy,  \
                             tk_zzzz_xxxyz,  \
                             tk_zzzz_xxxz,   \
                             tk_zzzz_xxxzz,  \
                             tk_zzzz_xxyy,   \
                             tk_zzzz_xxyyy,  \
                             tk_zzzz_xxyyz,  \
                             tk_zzzz_xxyz,   \
                             tk_zzzz_xxyzz,  \
                             tk_zzzz_xxzz,   \
                             tk_zzzz_xxzzz,  \
                             tk_zzzz_xyyy,   \
                             tk_zzzz_xyyyy,  \
                             tk_zzzz_xyyyz,  \
                             tk_zzzz_xyyz,   \
                             tk_zzzz_xyyzz,  \
                             tk_zzzz_xyzz,   \
                             tk_zzzz_xyzzz,  \
                             tk_zzzz_xzzz,   \
                             tk_zzzz_xzzzz,  \
                             tk_zzzz_yyyy,   \
                             tk_zzzz_yyyyy,  \
                             tk_zzzz_yyyyz,  \
                             tk_zzzz_yyyz,   \
                             tk_zzzz_yyyzz,  \
                             tk_zzzz_yyzz,   \
                             tk_zzzz_yyzzz,  \
                             tk_zzzz_yzzz,   \
                             tk_zzzz_yzzzz,  \
                             tk_zzzz_zzzz,   \
                             tk_zzzz_zzzzz,  \
                             ts_xzzzz_xxxxx, \
                             ts_xzzzz_xxxxy, \
                             ts_xzzzz_xxxxz, \
                             ts_xzzzz_xxxyy, \
                             ts_xzzzz_xxxyz, \
                             ts_xzzzz_xxxzz, \
                             ts_xzzzz_xxyyy, \
                             ts_xzzzz_xxyyz, \
                             ts_xzzzz_xxyzz, \
                             ts_xzzzz_xxzzz, \
                             ts_xzzzz_xyyyy, \
                             ts_xzzzz_xyyyz, \
                             ts_xzzzz_xyyzz, \
                             ts_xzzzz_xyzzz, \
                             ts_xzzzz_xzzzz, \
                             ts_xzzzz_yyyyy, \
                             ts_xzzzz_yyyyz, \
                             ts_xzzzz_yyyzz, \
                             ts_xzzzz_yyzzz, \
                             ts_xzzzz_yzzzz, \
                             ts_xzzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzz_xxxxx[i] = 5.0 * tk_zzzz_xxxx[i] * fe_0 + tk_zzzz_xxxxx[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxx[i] * fz_0;

        tk_xzzzz_xxxxy[i] = 4.0 * tk_zzzz_xxxy[i] * fe_0 + tk_zzzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxy[i] * fz_0;

        tk_xzzzz_xxxxz[i] = 4.0 * tk_zzzz_xxxz[i] * fe_0 + tk_zzzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxz[i] * fz_0;

        tk_xzzzz_xxxyy[i] = 3.0 * tk_zzzz_xxyy[i] * fe_0 + tk_zzzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxyy[i] * fz_0;

        tk_xzzzz_xxxyz[i] = 3.0 * tk_zzzz_xxyz[i] * fe_0 + tk_zzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxyz[i] * fz_0;

        tk_xzzzz_xxxzz[i] = 3.0 * tk_zzzz_xxzz[i] * fe_0 + tk_zzzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxzz[i] * fz_0;

        tk_xzzzz_xxyyy[i] = 2.0 * tk_zzzz_xyyy[i] * fe_0 + tk_zzzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyyy[i] * fz_0;

        tk_xzzzz_xxyyz[i] = 2.0 * tk_zzzz_xyyz[i] * fe_0 + tk_zzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyyz[i] * fz_0;

        tk_xzzzz_xxyzz[i] = 2.0 * tk_zzzz_xyzz[i] * fe_0 + tk_zzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyzz[i] * fz_0;

        tk_xzzzz_xxzzz[i] = 2.0 * tk_zzzz_xzzz[i] * fe_0 + tk_zzzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxzzz[i] * fz_0;

        tk_xzzzz_xyyyy[i] = tk_zzzz_yyyy[i] * fe_0 + tk_zzzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyyy[i] * fz_0;

        tk_xzzzz_xyyyz[i] = tk_zzzz_yyyz[i] * fe_0 + tk_zzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyyz[i] * fz_0;

        tk_xzzzz_xyyzz[i] = tk_zzzz_yyzz[i] * fe_0 + tk_zzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyzz[i] * fz_0;

        tk_xzzzz_xyzzz[i] = tk_zzzz_yzzz[i] * fe_0 + tk_zzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyzzz[i] * fz_0;

        tk_xzzzz_xzzzz[i] = tk_zzzz_zzzz[i] * fe_0 + tk_zzzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xzzzz[i] * fz_0;

        tk_xzzzz_yyyyy[i] = tk_zzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyyy[i] * fz_0;

        tk_xzzzz_yyyyz[i] = tk_zzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyyz[i] * fz_0;

        tk_xzzzz_yyyzz[i] = tk_zzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyzz[i] * fz_0;

        tk_xzzzz_yyzzz[i] = tk_zzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyzzz[i] * fz_0;

        tk_xzzzz_yzzzz[i] = tk_zzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yzzzz[i] * fz_0;

        tk_xzzzz_zzzzz[i] = tk_zzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_zzzzz[i] * fz_0;
    }

    // Set up 315-336 components of targeted buffer : HH

    auto tk_yyyyy_xxxxx = pbuffer.data(idx_kin_hh + 315);

    auto tk_yyyyy_xxxxy = pbuffer.data(idx_kin_hh + 316);

    auto tk_yyyyy_xxxxz = pbuffer.data(idx_kin_hh + 317);

    auto tk_yyyyy_xxxyy = pbuffer.data(idx_kin_hh + 318);

    auto tk_yyyyy_xxxyz = pbuffer.data(idx_kin_hh + 319);

    auto tk_yyyyy_xxxzz = pbuffer.data(idx_kin_hh + 320);

    auto tk_yyyyy_xxyyy = pbuffer.data(idx_kin_hh + 321);

    auto tk_yyyyy_xxyyz = pbuffer.data(idx_kin_hh + 322);

    auto tk_yyyyy_xxyzz = pbuffer.data(idx_kin_hh + 323);

    auto tk_yyyyy_xxzzz = pbuffer.data(idx_kin_hh + 324);

    auto tk_yyyyy_xyyyy = pbuffer.data(idx_kin_hh + 325);

    auto tk_yyyyy_xyyyz = pbuffer.data(idx_kin_hh + 326);

    auto tk_yyyyy_xyyzz = pbuffer.data(idx_kin_hh + 327);

    auto tk_yyyyy_xyzzz = pbuffer.data(idx_kin_hh + 328);

    auto tk_yyyyy_xzzzz = pbuffer.data(idx_kin_hh + 329);

    auto tk_yyyyy_yyyyy = pbuffer.data(idx_kin_hh + 330);

    auto tk_yyyyy_yyyyz = pbuffer.data(idx_kin_hh + 331);

    auto tk_yyyyy_yyyzz = pbuffer.data(idx_kin_hh + 332);

    auto tk_yyyyy_yyzzz = pbuffer.data(idx_kin_hh + 333);

    auto tk_yyyyy_yzzzz = pbuffer.data(idx_kin_hh + 334);

    auto tk_yyyyy_zzzzz = pbuffer.data(idx_kin_hh + 335);

#pragma omp simd aligned(pa_y,               \
                             tk_yyy_xxxxx,   \
                             tk_yyy_xxxxy,   \
                             tk_yyy_xxxxz,   \
                             tk_yyy_xxxyy,   \
                             tk_yyy_xxxyz,   \
                             tk_yyy_xxxzz,   \
                             tk_yyy_xxyyy,   \
                             tk_yyy_xxyyz,   \
                             tk_yyy_xxyzz,   \
                             tk_yyy_xxzzz,   \
                             tk_yyy_xyyyy,   \
                             tk_yyy_xyyyz,   \
                             tk_yyy_xyyzz,   \
                             tk_yyy_xyzzz,   \
                             tk_yyy_xzzzz,   \
                             tk_yyy_yyyyy,   \
                             tk_yyy_yyyyz,   \
                             tk_yyy_yyyzz,   \
                             tk_yyy_yyzzz,   \
                             tk_yyy_yzzzz,   \
                             tk_yyy_zzzzz,   \
                             tk_yyyy_xxxx,   \
                             tk_yyyy_xxxxx,  \
                             tk_yyyy_xxxxy,  \
                             tk_yyyy_xxxxz,  \
                             tk_yyyy_xxxy,   \
                             tk_yyyy_xxxyy,  \
                             tk_yyyy_xxxyz,  \
                             tk_yyyy_xxxz,   \
                             tk_yyyy_xxxzz,  \
                             tk_yyyy_xxyy,   \
                             tk_yyyy_xxyyy,  \
                             tk_yyyy_xxyyz,  \
                             tk_yyyy_xxyz,   \
                             tk_yyyy_xxyzz,  \
                             tk_yyyy_xxzz,   \
                             tk_yyyy_xxzzz,  \
                             tk_yyyy_xyyy,   \
                             tk_yyyy_xyyyy,  \
                             tk_yyyy_xyyyz,  \
                             tk_yyyy_xyyz,   \
                             tk_yyyy_xyyzz,  \
                             tk_yyyy_xyzz,   \
                             tk_yyyy_xyzzz,  \
                             tk_yyyy_xzzz,   \
                             tk_yyyy_xzzzz,  \
                             tk_yyyy_yyyy,   \
                             tk_yyyy_yyyyy,  \
                             tk_yyyy_yyyyz,  \
                             tk_yyyy_yyyz,   \
                             tk_yyyy_yyyzz,  \
                             tk_yyyy_yyzz,   \
                             tk_yyyy_yyzzz,  \
                             tk_yyyy_yzzz,   \
                             tk_yyyy_yzzzz,  \
                             tk_yyyy_zzzz,   \
                             tk_yyyy_zzzzz,  \
                             tk_yyyyy_xxxxx, \
                             tk_yyyyy_xxxxy, \
                             tk_yyyyy_xxxxz, \
                             tk_yyyyy_xxxyy, \
                             tk_yyyyy_xxxyz, \
                             tk_yyyyy_xxxzz, \
                             tk_yyyyy_xxyyy, \
                             tk_yyyyy_xxyyz, \
                             tk_yyyyy_xxyzz, \
                             tk_yyyyy_xxzzz, \
                             tk_yyyyy_xyyyy, \
                             tk_yyyyy_xyyyz, \
                             tk_yyyyy_xyyzz, \
                             tk_yyyyy_xyzzz, \
                             tk_yyyyy_xzzzz, \
                             tk_yyyyy_yyyyy, \
                             tk_yyyyy_yyyyz, \
                             tk_yyyyy_yyyzz, \
                             tk_yyyyy_yyzzz, \
                             tk_yyyyy_yzzzz, \
                             tk_yyyyy_zzzzz, \
                             ts_yyy_xxxxx,   \
                             ts_yyy_xxxxy,   \
                             ts_yyy_xxxxz,   \
                             ts_yyy_xxxyy,   \
                             ts_yyy_xxxyz,   \
                             ts_yyy_xxxzz,   \
                             ts_yyy_xxyyy,   \
                             ts_yyy_xxyyz,   \
                             ts_yyy_xxyzz,   \
                             ts_yyy_xxzzz,   \
                             ts_yyy_xyyyy,   \
                             ts_yyy_xyyyz,   \
                             ts_yyy_xyyzz,   \
                             ts_yyy_xyzzz,   \
                             ts_yyy_xzzzz,   \
                             ts_yyy_yyyyy,   \
                             ts_yyy_yyyyz,   \
                             ts_yyy_yyyzz,   \
                             ts_yyy_yyzzz,   \
                             ts_yyy_yzzzz,   \
                             ts_yyy_zzzzz,   \
                             ts_yyyyy_xxxxx, \
                             ts_yyyyy_xxxxy, \
                             ts_yyyyy_xxxxz, \
                             ts_yyyyy_xxxyy, \
                             ts_yyyyy_xxxyz, \
                             ts_yyyyy_xxxzz, \
                             ts_yyyyy_xxyyy, \
                             ts_yyyyy_xxyyz, \
                             ts_yyyyy_xxyzz, \
                             ts_yyyyy_xxzzz, \
                             ts_yyyyy_xyyyy, \
                             ts_yyyyy_xyyyz, \
                             ts_yyyyy_xyyzz, \
                             ts_yyyyy_xyzzz, \
                             ts_yyyyy_xzzzz, \
                             ts_yyyyy_yyyyy, \
                             ts_yyyyy_yyyyz, \
                             ts_yyyyy_yyyzz, \
                             ts_yyyyy_yyzzz, \
                             ts_yyyyy_yzzzz, \
                             ts_yyyyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyy_xxxxx[i] =
            -8.0 * ts_yyy_xxxxx[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxx[i] * fe_0 + tk_yyyy_xxxxx[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxx[i] * fz_0;

        tk_yyyyy_xxxxy[i] = -8.0 * ts_yyy_xxxxy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxy[i] * fe_0 + tk_yyyy_xxxx[i] * fe_0 +
                            tk_yyyy_xxxxy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxy[i] * fz_0;

        tk_yyyyy_xxxxz[i] =
            -8.0 * ts_yyy_xxxxz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxz[i] * fe_0 + tk_yyyy_xxxxz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxz[i] * fz_0;

        tk_yyyyy_xxxyy[i] = -8.0 * ts_yyy_xxxyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxyy[i] * fe_0 + 2.0 * tk_yyyy_xxxy[i] * fe_0 +
                            tk_yyyy_xxxyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxyy[i] * fz_0;

        tk_yyyyy_xxxyz[i] = -8.0 * ts_yyy_xxxyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxyz[i] * fe_0 + tk_yyyy_xxxz[i] * fe_0 +
                            tk_yyyy_xxxyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxyz[i] * fz_0;

        tk_yyyyy_xxxzz[i] =
            -8.0 * ts_yyy_xxxzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxzz[i] * fe_0 + tk_yyyy_xxxzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxzz[i] * fz_0;

        tk_yyyyy_xxyyy[i] = -8.0 * ts_yyy_xxyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyyy[i] * fe_0 + 3.0 * tk_yyyy_xxyy[i] * fe_0 +
                            tk_yyyy_xxyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyyy[i] * fz_0;

        tk_yyyyy_xxyyz[i] = -8.0 * ts_yyy_xxyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyyz[i] * fe_0 + 2.0 * tk_yyyy_xxyz[i] * fe_0 +
                            tk_yyyy_xxyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyyz[i] * fz_0;

        tk_yyyyy_xxyzz[i] = -8.0 * ts_yyy_xxyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyzz[i] * fe_0 + tk_yyyy_xxzz[i] * fe_0 +
                            tk_yyyy_xxyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyzz[i] * fz_0;

        tk_yyyyy_xxzzz[i] =
            -8.0 * ts_yyy_xxzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxzzz[i] * fe_0 + tk_yyyy_xxzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxzzz[i] * fz_0;

        tk_yyyyy_xyyyy[i] = -8.0 * ts_yyy_xyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyyy[i] * fe_0 + 4.0 * tk_yyyy_xyyy[i] * fe_0 +
                            tk_yyyy_xyyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyyy[i] * fz_0;

        tk_yyyyy_xyyyz[i] = -8.0 * ts_yyy_xyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyyz[i] * fe_0 + 3.0 * tk_yyyy_xyyz[i] * fe_0 +
                            tk_yyyy_xyyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyyz[i] * fz_0;

        tk_yyyyy_xyyzz[i] = -8.0 * ts_yyy_xyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyzz[i] * fe_0 + 2.0 * tk_yyyy_xyzz[i] * fe_0 +
                            tk_yyyy_xyyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyzz[i] * fz_0;

        tk_yyyyy_xyzzz[i] = -8.0 * ts_yyy_xyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyzzz[i] * fe_0 + tk_yyyy_xzzz[i] * fe_0 +
                            tk_yyyy_xyzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyzzz[i] * fz_0;

        tk_yyyyy_xzzzz[i] =
            -8.0 * ts_yyy_xzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xzzzz[i] * fe_0 + tk_yyyy_xzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xzzzz[i] * fz_0;

        tk_yyyyy_yyyyy[i] = -8.0 * ts_yyy_yyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyyy[i] * fe_0 + 5.0 * tk_yyyy_yyyy[i] * fe_0 +
                            tk_yyyy_yyyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyyy[i] * fz_0;

        tk_yyyyy_yyyyz[i] = -8.0 * ts_yyy_yyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyyz[i] * fe_0 + 4.0 * tk_yyyy_yyyz[i] * fe_0 +
                            tk_yyyy_yyyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyyz[i] * fz_0;

        tk_yyyyy_yyyzz[i] = -8.0 * ts_yyy_yyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyzz[i] * fe_0 + 3.0 * tk_yyyy_yyzz[i] * fe_0 +
                            tk_yyyy_yyyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyzz[i] * fz_0;

        tk_yyyyy_yyzzz[i] = -8.0 * ts_yyy_yyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyzzz[i] * fe_0 + 2.0 * tk_yyyy_yzzz[i] * fe_0 +
                            tk_yyyy_yyzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyzzz[i] * fz_0;

        tk_yyyyy_yzzzz[i] = -8.0 * ts_yyy_yzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yzzzz[i] * fe_0 + tk_yyyy_zzzz[i] * fe_0 +
                            tk_yyyy_yzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yzzzz[i] * fz_0;

        tk_yyyyy_zzzzz[i] =
            -8.0 * ts_yyy_zzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_zzzzz[i] * fe_0 + tk_yyyy_zzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_zzzzz[i] * fz_0;
    }

    // Set up 336-357 components of targeted buffer : HH

    auto tk_yyyyz_xxxxx = pbuffer.data(idx_kin_hh + 336);

    auto tk_yyyyz_xxxxy = pbuffer.data(idx_kin_hh + 337);

    auto tk_yyyyz_xxxxz = pbuffer.data(idx_kin_hh + 338);

    auto tk_yyyyz_xxxyy = pbuffer.data(idx_kin_hh + 339);

    auto tk_yyyyz_xxxyz = pbuffer.data(idx_kin_hh + 340);

    auto tk_yyyyz_xxxzz = pbuffer.data(idx_kin_hh + 341);

    auto tk_yyyyz_xxyyy = pbuffer.data(idx_kin_hh + 342);

    auto tk_yyyyz_xxyyz = pbuffer.data(idx_kin_hh + 343);

    auto tk_yyyyz_xxyzz = pbuffer.data(idx_kin_hh + 344);

    auto tk_yyyyz_xxzzz = pbuffer.data(idx_kin_hh + 345);

    auto tk_yyyyz_xyyyy = pbuffer.data(idx_kin_hh + 346);

    auto tk_yyyyz_xyyyz = pbuffer.data(idx_kin_hh + 347);

    auto tk_yyyyz_xyyzz = pbuffer.data(idx_kin_hh + 348);

    auto tk_yyyyz_xyzzz = pbuffer.data(idx_kin_hh + 349);

    auto tk_yyyyz_xzzzz = pbuffer.data(idx_kin_hh + 350);

    auto tk_yyyyz_yyyyy = pbuffer.data(idx_kin_hh + 351);

    auto tk_yyyyz_yyyyz = pbuffer.data(idx_kin_hh + 352);

    auto tk_yyyyz_yyyzz = pbuffer.data(idx_kin_hh + 353);

    auto tk_yyyyz_yyzzz = pbuffer.data(idx_kin_hh + 354);

    auto tk_yyyyz_yzzzz = pbuffer.data(idx_kin_hh + 355);

    auto tk_yyyyz_zzzzz = pbuffer.data(idx_kin_hh + 356);

#pragma omp simd aligned(pa_z,               \
                             tk_yyyy_xxxx,   \
                             tk_yyyy_xxxxx,  \
                             tk_yyyy_xxxxy,  \
                             tk_yyyy_xxxxz,  \
                             tk_yyyy_xxxy,   \
                             tk_yyyy_xxxyy,  \
                             tk_yyyy_xxxyz,  \
                             tk_yyyy_xxxz,   \
                             tk_yyyy_xxxzz,  \
                             tk_yyyy_xxyy,   \
                             tk_yyyy_xxyyy,  \
                             tk_yyyy_xxyyz,  \
                             tk_yyyy_xxyz,   \
                             tk_yyyy_xxyzz,  \
                             tk_yyyy_xxzz,   \
                             tk_yyyy_xxzzz,  \
                             tk_yyyy_xyyy,   \
                             tk_yyyy_xyyyy,  \
                             tk_yyyy_xyyyz,  \
                             tk_yyyy_xyyz,   \
                             tk_yyyy_xyyzz,  \
                             tk_yyyy_xyzz,   \
                             tk_yyyy_xyzzz,  \
                             tk_yyyy_xzzz,   \
                             tk_yyyy_xzzzz,  \
                             tk_yyyy_yyyy,   \
                             tk_yyyy_yyyyy,  \
                             tk_yyyy_yyyyz,  \
                             tk_yyyy_yyyz,   \
                             tk_yyyy_yyyzz,  \
                             tk_yyyy_yyzz,   \
                             tk_yyyy_yyzzz,  \
                             tk_yyyy_yzzz,   \
                             tk_yyyy_yzzzz,  \
                             tk_yyyy_zzzz,   \
                             tk_yyyy_zzzzz,  \
                             tk_yyyyz_xxxxx, \
                             tk_yyyyz_xxxxy, \
                             tk_yyyyz_xxxxz, \
                             tk_yyyyz_xxxyy, \
                             tk_yyyyz_xxxyz, \
                             tk_yyyyz_xxxzz, \
                             tk_yyyyz_xxyyy, \
                             tk_yyyyz_xxyyz, \
                             tk_yyyyz_xxyzz, \
                             tk_yyyyz_xxzzz, \
                             tk_yyyyz_xyyyy, \
                             tk_yyyyz_xyyyz, \
                             tk_yyyyz_xyyzz, \
                             tk_yyyyz_xyzzz, \
                             tk_yyyyz_xzzzz, \
                             tk_yyyyz_yyyyy, \
                             tk_yyyyz_yyyyz, \
                             tk_yyyyz_yyyzz, \
                             tk_yyyyz_yyzzz, \
                             tk_yyyyz_yzzzz, \
                             tk_yyyyz_zzzzz, \
                             ts_yyyyz_xxxxx, \
                             ts_yyyyz_xxxxy, \
                             ts_yyyyz_xxxxz, \
                             ts_yyyyz_xxxyy, \
                             ts_yyyyz_xxxyz, \
                             ts_yyyyz_xxxzz, \
                             ts_yyyyz_xxyyy, \
                             ts_yyyyz_xxyyz, \
                             ts_yyyyz_xxyzz, \
                             ts_yyyyz_xxzzz, \
                             ts_yyyyz_xyyyy, \
                             ts_yyyyz_xyyyz, \
                             ts_yyyyz_xyyzz, \
                             ts_yyyyz_xyzzz, \
                             ts_yyyyz_xzzzz, \
                             ts_yyyyz_yyyyy, \
                             ts_yyyyz_yyyyz, \
                             ts_yyyyz_yyyzz, \
                             ts_yyyyz_yyzzz, \
                             ts_yyyyz_yzzzz, \
                             ts_yyyyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyz_xxxxx[i] = tk_yyyy_xxxxx[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxx[i] * fz_0;

        tk_yyyyz_xxxxy[i] = tk_yyyy_xxxxy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxy[i] * fz_0;

        tk_yyyyz_xxxxz[i] = tk_yyyy_xxxx[i] * fe_0 + tk_yyyy_xxxxz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxz[i] * fz_0;

        tk_yyyyz_xxxyy[i] = tk_yyyy_xxxyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxyy[i] * fz_0;

        tk_yyyyz_xxxyz[i] = tk_yyyy_xxxy[i] * fe_0 + tk_yyyy_xxxyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxyz[i] * fz_0;

        tk_yyyyz_xxxzz[i] = 2.0 * tk_yyyy_xxxz[i] * fe_0 + tk_yyyy_xxxzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxzz[i] * fz_0;

        tk_yyyyz_xxyyy[i] = tk_yyyy_xxyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyyy[i] * fz_0;

        tk_yyyyz_xxyyz[i] = tk_yyyy_xxyy[i] * fe_0 + tk_yyyy_xxyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyyz[i] * fz_0;

        tk_yyyyz_xxyzz[i] = 2.0 * tk_yyyy_xxyz[i] * fe_0 + tk_yyyy_xxyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyzz[i] * fz_0;

        tk_yyyyz_xxzzz[i] = 3.0 * tk_yyyy_xxzz[i] * fe_0 + tk_yyyy_xxzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxzzz[i] * fz_0;

        tk_yyyyz_xyyyy[i] = tk_yyyy_xyyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyyy[i] * fz_0;

        tk_yyyyz_xyyyz[i] = tk_yyyy_xyyy[i] * fe_0 + tk_yyyy_xyyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyyz[i] * fz_0;

        tk_yyyyz_xyyzz[i] = 2.0 * tk_yyyy_xyyz[i] * fe_0 + tk_yyyy_xyyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyzz[i] * fz_0;

        tk_yyyyz_xyzzz[i] = 3.0 * tk_yyyy_xyzz[i] * fe_0 + tk_yyyy_xyzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyzzz[i] * fz_0;

        tk_yyyyz_xzzzz[i] = 4.0 * tk_yyyy_xzzz[i] * fe_0 + tk_yyyy_xzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xzzzz[i] * fz_0;

        tk_yyyyz_yyyyy[i] = tk_yyyy_yyyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyyy[i] * fz_0;

        tk_yyyyz_yyyyz[i] = tk_yyyy_yyyy[i] * fe_0 + tk_yyyy_yyyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyyz[i] * fz_0;

        tk_yyyyz_yyyzz[i] = 2.0 * tk_yyyy_yyyz[i] * fe_0 + tk_yyyy_yyyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyzz[i] * fz_0;

        tk_yyyyz_yyzzz[i] = 3.0 * tk_yyyy_yyzz[i] * fe_0 + tk_yyyy_yyzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyzzz[i] * fz_0;

        tk_yyyyz_yzzzz[i] = 4.0 * tk_yyyy_yzzz[i] * fe_0 + tk_yyyy_yzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yzzzz[i] * fz_0;

        tk_yyyyz_zzzzz[i] = 5.0 * tk_yyyy_zzzz[i] * fe_0 + tk_yyyy_zzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_zzzzz[i] * fz_0;
    }

    // Set up 357-378 components of targeted buffer : HH

    auto tk_yyyzz_xxxxx = pbuffer.data(idx_kin_hh + 357);

    auto tk_yyyzz_xxxxy = pbuffer.data(idx_kin_hh + 358);

    auto tk_yyyzz_xxxxz = pbuffer.data(idx_kin_hh + 359);

    auto tk_yyyzz_xxxyy = pbuffer.data(idx_kin_hh + 360);

    auto tk_yyyzz_xxxyz = pbuffer.data(idx_kin_hh + 361);

    auto tk_yyyzz_xxxzz = pbuffer.data(idx_kin_hh + 362);

    auto tk_yyyzz_xxyyy = pbuffer.data(idx_kin_hh + 363);

    auto tk_yyyzz_xxyyz = pbuffer.data(idx_kin_hh + 364);

    auto tk_yyyzz_xxyzz = pbuffer.data(idx_kin_hh + 365);

    auto tk_yyyzz_xxzzz = pbuffer.data(idx_kin_hh + 366);

    auto tk_yyyzz_xyyyy = pbuffer.data(idx_kin_hh + 367);

    auto tk_yyyzz_xyyyz = pbuffer.data(idx_kin_hh + 368);

    auto tk_yyyzz_xyyzz = pbuffer.data(idx_kin_hh + 369);

    auto tk_yyyzz_xyzzz = pbuffer.data(idx_kin_hh + 370);

    auto tk_yyyzz_xzzzz = pbuffer.data(idx_kin_hh + 371);

    auto tk_yyyzz_yyyyy = pbuffer.data(idx_kin_hh + 372);

    auto tk_yyyzz_yyyyz = pbuffer.data(idx_kin_hh + 373);

    auto tk_yyyzz_yyyzz = pbuffer.data(idx_kin_hh + 374);

    auto tk_yyyzz_yyzzz = pbuffer.data(idx_kin_hh + 375);

    auto tk_yyyzz_yzzzz = pbuffer.data(idx_kin_hh + 376);

    auto tk_yyyzz_zzzzz = pbuffer.data(idx_kin_hh + 377);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_yyy_xxxxy,   \
                             tk_yyy_xxxyy,   \
                             tk_yyy_xxyyy,   \
                             tk_yyy_xyyyy,   \
                             tk_yyy_yyyyy,   \
                             tk_yyyz_xxxxy,  \
                             tk_yyyz_xxxyy,  \
                             tk_yyyz_xxyyy,  \
                             tk_yyyz_xyyyy,  \
                             tk_yyyz_yyyyy,  \
                             tk_yyyzz_xxxxx, \
                             tk_yyyzz_xxxxy, \
                             tk_yyyzz_xxxxz, \
                             tk_yyyzz_xxxyy, \
                             tk_yyyzz_xxxyz, \
                             tk_yyyzz_xxxzz, \
                             tk_yyyzz_xxyyy, \
                             tk_yyyzz_xxyyz, \
                             tk_yyyzz_xxyzz, \
                             tk_yyyzz_xxzzz, \
                             tk_yyyzz_xyyyy, \
                             tk_yyyzz_xyyyz, \
                             tk_yyyzz_xyyzz, \
                             tk_yyyzz_xyzzz, \
                             tk_yyyzz_xzzzz, \
                             tk_yyyzz_yyyyy, \
                             tk_yyyzz_yyyyz, \
                             tk_yyyzz_yyyzz, \
                             tk_yyyzz_yyzzz, \
                             tk_yyyzz_yzzzz, \
                             tk_yyyzz_zzzzz, \
                             tk_yyzz_xxxxx,  \
                             tk_yyzz_xxxxz,  \
                             tk_yyzz_xxxyz,  \
                             tk_yyzz_xxxz,   \
                             tk_yyzz_xxxzz,  \
                             tk_yyzz_xxyyz,  \
                             tk_yyzz_xxyz,   \
                             tk_yyzz_xxyzz,  \
                             tk_yyzz_xxzz,   \
                             tk_yyzz_xxzzz,  \
                             tk_yyzz_xyyyz,  \
                             tk_yyzz_xyyz,   \
                             tk_yyzz_xyyzz,  \
                             tk_yyzz_xyzz,   \
                             tk_yyzz_xyzzz,  \
                             tk_yyzz_xzzz,   \
                             tk_yyzz_xzzzz,  \
                             tk_yyzz_yyyyz,  \
                             tk_yyzz_yyyz,   \
                             tk_yyzz_yyyzz,  \
                             tk_yyzz_yyzz,   \
                             tk_yyzz_yyzzz,  \
                             tk_yyzz_yzzz,   \
                             tk_yyzz_yzzzz,  \
                             tk_yyzz_zzzz,   \
                             tk_yyzz_zzzzz,  \
                             tk_yzz_xxxxx,   \
                             tk_yzz_xxxxz,   \
                             tk_yzz_xxxyz,   \
                             tk_yzz_xxxzz,   \
                             tk_yzz_xxyyz,   \
                             tk_yzz_xxyzz,   \
                             tk_yzz_xxzzz,   \
                             tk_yzz_xyyyz,   \
                             tk_yzz_xyyzz,   \
                             tk_yzz_xyzzz,   \
                             tk_yzz_xzzzz,   \
                             tk_yzz_yyyyz,   \
                             tk_yzz_yyyzz,   \
                             tk_yzz_yyzzz,   \
                             tk_yzz_yzzzz,   \
                             tk_yzz_zzzzz,   \
                             ts_yyy_xxxxy,   \
                             ts_yyy_xxxyy,   \
                             ts_yyy_xxyyy,   \
                             ts_yyy_xyyyy,   \
                             ts_yyy_yyyyy,   \
                             ts_yyyzz_xxxxx, \
                             ts_yyyzz_xxxxy, \
                             ts_yyyzz_xxxxz, \
                             ts_yyyzz_xxxyy, \
                             ts_yyyzz_xxxyz, \
                             ts_yyyzz_xxxzz, \
                             ts_yyyzz_xxyyy, \
                             ts_yyyzz_xxyyz, \
                             ts_yyyzz_xxyzz, \
                             ts_yyyzz_xxzzz, \
                             ts_yyyzz_xyyyy, \
                             ts_yyyzz_xyyyz, \
                             ts_yyyzz_xyyzz, \
                             ts_yyyzz_xyzzz, \
                             ts_yyyzz_xzzzz, \
                             ts_yyyzz_yyyyy, \
                             ts_yyyzz_yyyyz, \
                             ts_yyyzz_yyyzz, \
                             ts_yyyzz_yyzzz, \
                             ts_yyyzz_yzzzz, \
                             ts_yyyzz_zzzzz, \
                             ts_yzz_xxxxx,   \
                             ts_yzz_xxxxz,   \
                             ts_yzz_xxxyz,   \
                             ts_yzz_xxxzz,   \
                             ts_yzz_xxyyz,   \
                             ts_yzz_xxyzz,   \
                             ts_yzz_xxzzz,   \
                             ts_yzz_xyyyz,   \
                             ts_yzz_xyyzz,   \
                             ts_yzz_xyzzz,   \
                             ts_yzz_xzzzz,   \
                             ts_yzz_yyyyz,   \
                             ts_yzz_yyyzz,   \
                             ts_yzz_yyzzz,   \
                             ts_yzz_yzzzz,   \
                             ts_yzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzz_xxxxx[i] =
            -4.0 * ts_yzz_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxxx[i] * fe_0 + tk_yyzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxxx[i] * fz_0;

        tk_yyyzz_xxxxy[i] =
            -2.0 * ts_yyy_xxxxy[i] * fbe_0 * fz_0 + tk_yyy_xxxxy[i] * fe_0 + tk_yyyz_xxxxy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxxxy[i] * fz_0;

        tk_yyyzz_xxxxz[i] =
            -4.0 * ts_yzz_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxxz[i] * fe_0 + tk_yyzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxxz[i] * fz_0;

        tk_yyyzz_xxxyy[i] =
            -2.0 * ts_yyy_xxxyy[i] * fbe_0 * fz_0 + tk_yyy_xxxyy[i] * fe_0 + tk_yyyz_xxxyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxxyy[i] * fz_0;

        tk_yyyzz_xxxyz[i] = -4.0 * ts_yzz_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxyz[i] * fe_0 + tk_yyzz_xxxz[i] * fe_0 +
                            tk_yyzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxyz[i] * fz_0;

        tk_yyyzz_xxxzz[i] =
            -4.0 * ts_yzz_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxzz[i] * fe_0 + tk_yyzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxzz[i] * fz_0;

        tk_yyyzz_xxyyy[i] =
            -2.0 * ts_yyy_xxyyy[i] * fbe_0 * fz_0 + tk_yyy_xxyyy[i] * fe_0 + tk_yyyz_xxyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxyyy[i] * fz_0;

        tk_yyyzz_xxyyz[i] = -4.0 * ts_yzz_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxyyz[i] * fe_0 + 2.0 * tk_yyzz_xxyz[i] * fe_0 +
                            tk_yyzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxyyz[i] * fz_0;

        tk_yyyzz_xxyzz[i] = -4.0 * ts_yzz_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxyzz[i] * fe_0 + tk_yyzz_xxzz[i] * fe_0 +
                            tk_yyzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxyzz[i] * fz_0;

        tk_yyyzz_xxzzz[i] =
            -4.0 * ts_yzz_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxzzz[i] * fe_0 + tk_yyzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxzzz[i] * fz_0;

        tk_yyyzz_xyyyy[i] =
            -2.0 * ts_yyy_xyyyy[i] * fbe_0 * fz_0 + tk_yyy_xyyyy[i] * fe_0 + tk_yyyz_xyyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xyyyy[i] * fz_0;

        tk_yyyzz_xyyyz[i] = -4.0 * ts_yzz_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyyyz[i] * fe_0 + 3.0 * tk_yyzz_xyyz[i] * fe_0 +
                            tk_yyzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyyyz[i] * fz_0;

        tk_yyyzz_xyyzz[i] = -4.0 * ts_yzz_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyyzz[i] * fe_0 + 2.0 * tk_yyzz_xyzz[i] * fe_0 +
                            tk_yyzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyyzz[i] * fz_0;

        tk_yyyzz_xyzzz[i] = -4.0 * ts_yzz_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyzzz[i] * fe_0 + tk_yyzz_xzzz[i] * fe_0 +
                            tk_yyzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyzzz[i] * fz_0;

        tk_yyyzz_xzzzz[i] =
            -4.0 * ts_yzz_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xzzzz[i] * fe_0 + tk_yyzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xzzzz[i] * fz_0;

        tk_yyyzz_yyyyy[i] =
            -2.0 * ts_yyy_yyyyy[i] * fbe_0 * fz_0 + tk_yyy_yyyyy[i] * fe_0 + tk_yyyz_yyyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_yyyyy[i] * fz_0;

        tk_yyyzz_yyyyz[i] = -4.0 * ts_yzz_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyyyz[i] * fe_0 + 4.0 * tk_yyzz_yyyz[i] * fe_0 +
                            tk_yyzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyyyz[i] * fz_0;

        tk_yyyzz_yyyzz[i] = -4.0 * ts_yzz_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyyzz[i] * fe_0 + 3.0 * tk_yyzz_yyzz[i] * fe_0 +
                            tk_yyzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyyzz[i] * fz_0;

        tk_yyyzz_yyzzz[i] = -4.0 * ts_yzz_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyzzz[i] * fe_0 + 2.0 * tk_yyzz_yzzz[i] * fe_0 +
                            tk_yyzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyzzz[i] * fz_0;

        tk_yyyzz_yzzzz[i] = -4.0 * ts_yzz_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yzzzz[i] * fe_0 + tk_yyzz_zzzz[i] * fe_0 +
                            tk_yyzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yzzzz[i] * fz_0;

        tk_yyyzz_zzzzz[i] =
            -4.0 * ts_yzz_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_zzzzz[i] * fe_0 + tk_yyzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_zzzzz[i] * fz_0;
    }

    // Set up 378-399 components of targeted buffer : HH

    auto tk_yyzzz_xxxxx = pbuffer.data(idx_kin_hh + 378);

    auto tk_yyzzz_xxxxy = pbuffer.data(idx_kin_hh + 379);

    auto tk_yyzzz_xxxxz = pbuffer.data(idx_kin_hh + 380);

    auto tk_yyzzz_xxxyy = pbuffer.data(idx_kin_hh + 381);

    auto tk_yyzzz_xxxyz = pbuffer.data(idx_kin_hh + 382);

    auto tk_yyzzz_xxxzz = pbuffer.data(idx_kin_hh + 383);

    auto tk_yyzzz_xxyyy = pbuffer.data(idx_kin_hh + 384);

    auto tk_yyzzz_xxyyz = pbuffer.data(idx_kin_hh + 385);

    auto tk_yyzzz_xxyzz = pbuffer.data(idx_kin_hh + 386);

    auto tk_yyzzz_xxzzz = pbuffer.data(idx_kin_hh + 387);

    auto tk_yyzzz_xyyyy = pbuffer.data(idx_kin_hh + 388);

    auto tk_yyzzz_xyyyz = pbuffer.data(idx_kin_hh + 389);

    auto tk_yyzzz_xyyzz = pbuffer.data(idx_kin_hh + 390);

    auto tk_yyzzz_xyzzz = pbuffer.data(idx_kin_hh + 391);

    auto tk_yyzzz_xzzzz = pbuffer.data(idx_kin_hh + 392);

    auto tk_yyzzz_yyyyy = pbuffer.data(idx_kin_hh + 393);

    auto tk_yyzzz_yyyyz = pbuffer.data(idx_kin_hh + 394);

    auto tk_yyzzz_yyyzz = pbuffer.data(idx_kin_hh + 395);

    auto tk_yyzzz_yyzzz = pbuffer.data(idx_kin_hh + 396);

    auto tk_yyzzz_yzzzz = pbuffer.data(idx_kin_hh + 397);

    auto tk_yyzzz_zzzzz = pbuffer.data(idx_kin_hh + 398);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_yyz_xxxxy,   \
                             tk_yyz_xxxyy,   \
                             tk_yyz_xxyyy,   \
                             tk_yyz_xyyyy,   \
                             tk_yyz_yyyyy,   \
                             tk_yyzz_xxxxy,  \
                             tk_yyzz_xxxyy,  \
                             tk_yyzz_xxyyy,  \
                             tk_yyzz_xyyyy,  \
                             tk_yyzz_yyyyy,  \
                             tk_yyzzz_xxxxx, \
                             tk_yyzzz_xxxxy, \
                             tk_yyzzz_xxxxz, \
                             tk_yyzzz_xxxyy, \
                             tk_yyzzz_xxxyz, \
                             tk_yyzzz_xxxzz, \
                             tk_yyzzz_xxyyy, \
                             tk_yyzzz_xxyyz, \
                             tk_yyzzz_xxyzz, \
                             tk_yyzzz_xxzzz, \
                             tk_yyzzz_xyyyy, \
                             tk_yyzzz_xyyyz, \
                             tk_yyzzz_xyyzz, \
                             tk_yyzzz_xyzzz, \
                             tk_yyzzz_xzzzz, \
                             tk_yyzzz_yyyyy, \
                             tk_yyzzz_yyyyz, \
                             tk_yyzzz_yyyzz, \
                             tk_yyzzz_yyzzz, \
                             tk_yyzzz_yzzzz, \
                             tk_yyzzz_zzzzz, \
                             tk_yzzz_xxxxx,  \
                             tk_yzzz_xxxxz,  \
                             tk_yzzz_xxxyz,  \
                             tk_yzzz_xxxz,   \
                             tk_yzzz_xxxzz,  \
                             tk_yzzz_xxyyz,  \
                             tk_yzzz_xxyz,   \
                             tk_yzzz_xxyzz,  \
                             tk_yzzz_xxzz,   \
                             tk_yzzz_xxzzz,  \
                             tk_yzzz_xyyyz,  \
                             tk_yzzz_xyyz,   \
                             tk_yzzz_xyyzz,  \
                             tk_yzzz_xyzz,   \
                             tk_yzzz_xyzzz,  \
                             tk_yzzz_xzzz,   \
                             tk_yzzz_xzzzz,  \
                             tk_yzzz_yyyyz,  \
                             tk_yzzz_yyyz,   \
                             tk_yzzz_yyyzz,  \
                             tk_yzzz_yyzz,   \
                             tk_yzzz_yyzzz,  \
                             tk_yzzz_yzzz,   \
                             tk_yzzz_yzzzz,  \
                             tk_yzzz_zzzz,   \
                             tk_yzzz_zzzzz,  \
                             tk_zzz_xxxxx,   \
                             tk_zzz_xxxxz,   \
                             tk_zzz_xxxyz,   \
                             tk_zzz_xxxzz,   \
                             tk_zzz_xxyyz,   \
                             tk_zzz_xxyzz,   \
                             tk_zzz_xxzzz,   \
                             tk_zzz_xyyyz,   \
                             tk_zzz_xyyzz,   \
                             tk_zzz_xyzzz,   \
                             tk_zzz_xzzzz,   \
                             tk_zzz_yyyyz,   \
                             tk_zzz_yyyzz,   \
                             tk_zzz_yyzzz,   \
                             tk_zzz_yzzzz,   \
                             tk_zzz_zzzzz,   \
                             ts_yyz_xxxxy,   \
                             ts_yyz_xxxyy,   \
                             ts_yyz_xxyyy,   \
                             ts_yyz_xyyyy,   \
                             ts_yyz_yyyyy,   \
                             ts_yyzzz_xxxxx, \
                             ts_yyzzz_xxxxy, \
                             ts_yyzzz_xxxxz, \
                             ts_yyzzz_xxxyy, \
                             ts_yyzzz_xxxyz, \
                             ts_yyzzz_xxxzz, \
                             ts_yyzzz_xxyyy, \
                             ts_yyzzz_xxyyz, \
                             ts_yyzzz_xxyzz, \
                             ts_yyzzz_xxzzz, \
                             ts_yyzzz_xyyyy, \
                             ts_yyzzz_xyyyz, \
                             ts_yyzzz_xyyzz, \
                             ts_yyzzz_xyzzz, \
                             ts_yyzzz_xzzzz, \
                             ts_yyzzz_yyyyy, \
                             ts_yyzzz_yyyyz, \
                             ts_yyzzz_yyyzz, \
                             ts_yyzzz_yyzzz, \
                             ts_yyzzz_yzzzz, \
                             ts_yyzzz_zzzzz, \
                             ts_zzz_xxxxx,   \
                             ts_zzz_xxxxz,   \
                             ts_zzz_xxxyz,   \
                             ts_zzz_xxxzz,   \
                             ts_zzz_xxyyz,   \
                             ts_zzz_xxyzz,   \
                             ts_zzz_xxzzz,   \
                             ts_zzz_xyyyz,   \
                             ts_zzz_xyyzz,   \
                             ts_zzz_xyzzz,   \
                             ts_zzz_xzzzz,   \
                             ts_zzz_yyyyz,   \
                             ts_zzz_yyyzz,   \
                             ts_zzz_yyzzz,   \
                             ts_zzz_yzzzz,   \
                             ts_zzz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzz_xxxxx[i] =
            -2.0 * ts_zzz_xxxxx[i] * fbe_0 * fz_0 + tk_zzz_xxxxx[i] * fe_0 + tk_yzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxxx[i] * fz_0;

        tk_yyzzz_xxxxy[i] =
            -4.0 * ts_yyz_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxxxy[i] * fe_0 + tk_yyzz_xxxxy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxxxy[i] * fz_0;

        tk_yyzzz_xxxxz[i] =
            -2.0 * ts_zzz_xxxxz[i] * fbe_0 * fz_0 + tk_zzz_xxxxz[i] * fe_0 + tk_yzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxxz[i] * fz_0;

        tk_yyzzz_xxxyy[i] =
            -4.0 * ts_yyz_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxxyy[i] * fe_0 + tk_yyzz_xxxyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxxyy[i] * fz_0;

        tk_yyzzz_xxxyz[i] = -2.0 * ts_zzz_xxxyz[i] * fbe_0 * fz_0 + tk_zzz_xxxyz[i] * fe_0 + tk_yzzz_xxxz[i] * fe_0 + tk_yzzz_xxxyz[i] * pa_y[i] +
                            2.0 * ts_yyzzz_xxxyz[i] * fz_0;

        tk_yyzzz_xxxzz[i] =
            -2.0 * ts_zzz_xxxzz[i] * fbe_0 * fz_0 + tk_zzz_xxxzz[i] * fe_0 + tk_yzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxzz[i] * fz_0;

        tk_yyzzz_xxyyy[i] =
            -4.0 * ts_yyz_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxyyy[i] * fe_0 + tk_yyzz_xxyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxyyy[i] * fz_0;

        tk_yyzzz_xxyyz[i] = -2.0 * ts_zzz_xxyyz[i] * fbe_0 * fz_0 + tk_zzz_xxyyz[i] * fe_0 + 2.0 * tk_yzzz_xxyz[i] * fe_0 +
                            tk_yzzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxyyz[i] * fz_0;

        tk_yyzzz_xxyzz[i] = -2.0 * ts_zzz_xxyzz[i] * fbe_0 * fz_0 + tk_zzz_xxyzz[i] * fe_0 + tk_yzzz_xxzz[i] * fe_0 + tk_yzzz_xxyzz[i] * pa_y[i] +
                            2.0 * ts_yyzzz_xxyzz[i] * fz_0;

        tk_yyzzz_xxzzz[i] =
            -2.0 * ts_zzz_xxzzz[i] * fbe_0 * fz_0 + tk_zzz_xxzzz[i] * fe_0 + tk_yzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxzzz[i] * fz_0;

        tk_yyzzz_xyyyy[i] =
            -4.0 * ts_yyz_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xyyyy[i] * fe_0 + tk_yyzz_xyyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xyyyy[i] * fz_0;

        tk_yyzzz_xyyyz[i] = -2.0 * ts_zzz_xyyyz[i] * fbe_0 * fz_0 + tk_zzz_xyyyz[i] * fe_0 + 3.0 * tk_yzzz_xyyz[i] * fe_0 +
                            tk_yzzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyyyz[i] * fz_0;

        tk_yyzzz_xyyzz[i] = -2.0 * ts_zzz_xyyzz[i] * fbe_0 * fz_0 + tk_zzz_xyyzz[i] * fe_0 + 2.0 * tk_yzzz_xyzz[i] * fe_0 +
                            tk_yzzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyyzz[i] * fz_0;

        tk_yyzzz_xyzzz[i] = -2.0 * ts_zzz_xyzzz[i] * fbe_0 * fz_0 + tk_zzz_xyzzz[i] * fe_0 + tk_yzzz_xzzz[i] * fe_0 + tk_yzzz_xyzzz[i] * pa_y[i] +
                            2.0 * ts_yyzzz_xyzzz[i] * fz_0;

        tk_yyzzz_xzzzz[i] =
            -2.0 * ts_zzz_xzzzz[i] * fbe_0 * fz_0 + tk_zzz_xzzzz[i] * fe_0 + tk_yzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xzzzz[i] * fz_0;

        tk_yyzzz_yyyyy[i] =
            -4.0 * ts_yyz_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_yyyyy[i] * fe_0 + tk_yyzz_yyyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_yyyyy[i] * fz_0;

        tk_yyzzz_yyyyz[i] = -2.0 * ts_zzz_yyyyz[i] * fbe_0 * fz_0 + tk_zzz_yyyyz[i] * fe_0 + 4.0 * tk_yzzz_yyyz[i] * fe_0 +
                            tk_yzzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyyyz[i] * fz_0;

        tk_yyzzz_yyyzz[i] = -2.0 * ts_zzz_yyyzz[i] * fbe_0 * fz_0 + tk_zzz_yyyzz[i] * fe_0 + 3.0 * tk_yzzz_yyzz[i] * fe_0 +
                            tk_yzzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyyzz[i] * fz_0;

        tk_yyzzz_yyzzz[i] = -2.0 * ts_zzz_yyzzz[i] * fbe_0 * fz_0 + tk_zzz_yyzzz[i] * fe_0 + 2.0 * tk_yzzz_yzzz[i] * fe_0 +
                            tk_yzzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyzzz[i] * fz_0;

        tk_yyzzz_yzzzz[i] = -2.0 * ts_zzz_yzzzz[i] * fbe_0 * fz_0 + tk_zzz_yzzzz[i] * fe_0 + tk_yzzz_zzzz[i] * fe_0 + tk_yzzz_yzzzz[i] * pa_y[i] +
                            2.0 * ts_yyzzz_yzzzz[i] * fz_0;

        tk_yyzzz_zzzzz[i] =
            -2.0 * ts_zzz_zzzzz[i] * fbe_0 * fz_0 + tk_zzz_zzzzz[i] * fe_0 + tk_yzzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_zzzzz[i] * fz_0;
    }

    // Set up 399-420 components of targeted buffer : HH

    auto tk_yzzzz_xxxxx = pbuffer.data(idx_kin_hh + 399);

    auto tk_yzzzz_xxxxy = pbuffer.data(idx_kin_hh + 400);

    auto tk_yzzzz_xxxxz = pbuffer.data(idx_kin_hh + 401);

    auto tk_yzzzz_xxxyy = pbuffer.data(idx_kin_hh + 402);

    auto tk_yzzzz_xxxyz = pbuffer.data(idx_kin_hh + 403);

    auto tk_yzzzz_xxxzz = pbuffer.data(idx_kin_hh + 404);

    auto tk_yzzzz_xxyyy = pbuffer.data(idx_kin_hh + 405);

    auto tk_yzzzz_xxyyz = pbuffer.data(idx_kin_hh + 406);

    auto tk_yzzzz_xxyzz = pbuffer.data(idx_kin_hh + 407);

    auto tk_yzzzz_xxzzz = pbuffer.data(idx_kin_hh + 408);

    auto tk_yzzzz_xyyyy = pbuffer.data(idx_kin_hh + 409);

    auto tk_yzzzz_xyyyz = pbuffer.data(idx_kin_hh + 410);

    auto tk_yzzzz_xyyzz = pbuffer.data(idx_kin_hh + 411);

    auto tk_yzzzz_xyzzz = pbuffer.data(idx_kin_hh + 412);

    auto tk_yzzzz_xzzzz = pbuffer.data(idx_kin_hh + 413);

    auto tk_yzzzz_yyyyy = pbuffer.data(idx_kin_hh + 414);

    auto tk_yzzzz_yyyyz = pbuffer.data(idx_kin_hh + 415);

    auto tk_yzzzz_yyyzz = pbuffer.data(idx_kin_hh + 416);

    auto tk_yzzzz_yyzzz = pbuffer.data(idx_kin_hh + 417);

    auto tk_yzzzz_yzzzz = pbuffer.data(idx_kin_hh + 418);

    auto tk_yzzzz_zzzzz = pbuffer.data(idx_kin_hh + 419);

#pragma omp simd aligned(pa_y,               \
                             tk_yzzzz_xxxxx, \
                             tk_yzzzz_xxxxy, \
                             tk_yzzzz_xxxxz, \
                             tk_yzzzz_xxxyy, \
                             tk_yzzzz_xxxyz, \
                             tk_yzzzz_xxxzz, \
                             tk_yzzzz_xxyyy, \
                             tk_yzzzz_xxyyz, \
                             tk_yzzzz_xxyzz, \
                             tk_yzzzz_xxzzz, \
                             tk_yzzzz_xyyyy, \
                             tk_yzzzz_xyyyz, \
                             tk_yzzzz_xyyzz, \
                             tk_yzzzz_xyzzz, \
                             tk_yzzzz_xzzzz, \
                             tk_yzzzz_yyyyy, \
                             tk_yzzzz_yyyyz, \
                             tk_yzzzz_yyyzz, \
                             tk_yzzzz_yyzzz, \
                             tk_yzzzz_yzzzz, \
                             tk_yzzzz_zzzzz, \
                             tk_zzzz_xxxx,   \
                             tk_zzzz_xxxxx,  \
                             tk_zzzz_xxxxy,  \
                             tk_zzzz_xxxxz,  \
                             tk_zzzz_xxxy,   \
                             tk_zzzz_xxxyy,  \
                             tk_zzzz_xxxyz,  \
                             tk_zzzz_xxxz,   \
                             tk_zzzz_xxxzz,  \
                             tk_zzzz_xxyy,   \
                             tk_zzzz_xxyyy,  \
                             tk_zzzz_xxyyz,  \
                             tk_zzzz_xxyz,   \
                             tk_zzzz_xxyzz,  \
                             tk_zzzz_xxzz,   \
                             tk_zzzz_xxzzz,  \
                             tk_zzzz_xyyy,   \
                             tk_zzzz_xyyyy,  \
                             tk_zzzz_xyyyz,  \
                             tk_zzzz_xyyz,   \
                             tk_zzzz_xyyzz,  \
                             tk_zzzz_xyzz,   \
                             tk_zzzz_xyzzz,  \
                             tk_zzzz_xzzz,   \
                             tk_zzzz_xzzzz,  \
                             tk_zzzz_yyyy,   \
                             tk_zzzz_yyyyy,  \
                             tk_zzzz_yyyyz,  \
                             tk_zzzz_yyyz,   \
                             tk_zzzz_yyyzz,  \
                             tk_zzzz_yyzz,   \
                             tk_zzzz_yyzzz,  \
                             tk_zzzz_yzzz,   \
                             tk_zzzz_yzzzz,  \
                             tk_zzzz_zzzz,   \
                             tk_zzzz_zzzzz,  \
                             ts_yzzzz_xxxxx, \
                             ts_yzzzz_xxxxy, \
                             ts_yzzzz_xxxxz, \
                             ts_yzzzz_xxxyy, \
                             ts_yzzzz_xxxyz, \
                             ts_yzzzz_xxxzz, \
                             ts_yzzzz_xxyyy, \
                             ts_yzzzz_xxyyz, \
                             ts_yzzzz_xxyzz, \
                             ts_yzzzz_xxzzz, \
                             ts_yzzzz_xyyyy, \
                             ts_yzzzz_xyyyz, \
                             ts_yzzzz_xyyzz, \
                             ts_yzzzz_xyzzz, \
                             ts_yzzzz_xzzzz, \
                             ts_yzzzz_yyyyy, \
                             ts_yzzzz_yyyyz, \
                             ts_yzzzz_yyyzz, \
                             ts_yzzzz_yyzzz, \
                             ts_yzzzz_yzzzz, \
                             ts_yzzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzz_xxxxx[i] = tk_zzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxx[i] * fz_0;

        tk_yzzzz_xxxxy[i] = tk_zzzz_xxxx[i] * fe_0 + tk_zzzz_xxxxy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxy[i] * fz_0;

        tk_yzzzz_xxxxz[i] = tk_zzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxz[i] * fz_0;

        tk_yzzzz_xxxyy[i] = 2.0 * tk_zzzz_xxxy[i] * fe_0 + tk_zzzz_xxxyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxyy[i] * fz_0;

        tk_yzzzz_xxxyz[i] = tk_zzzz_xxxz[i] * fe_0 + tk_zzzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxyz[i] * fz_0;

        tk_yzzzz_xxxzz[i] = tk_zzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxzz[i] * fz_0;

        tk_yzzzz_xxyyy[i] = 3.0 * tk_zzzz_xxyy[i] * fe_0 + tk_zzzz_xxyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyyy[i] * fz_0;

        tk_yzzzz_xxyyz[i] = 2.0 * tk_zzzz_xxyz[i] * fe_0 + tk_zzzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyyz[i] * fz_0;

        tk_yzzzz_xxyzz[i] = tk_zzzz_xxzz[i] * fe_0 + tk_zzzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyzz[i] * fz_0;

        tk_yzzzz_xxzzz[i] = tk_zzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxzzz[i] * fz_0;

        tk_yzzzz_xyyyy[i] = 4.0 * tk_zzzz_xyyy[i] * fe_0 + tk_zzzz_xyyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyyy[i] * fz_0;

        tk_yzzzz_xyyyz[i] = 3.0 * tk_zzzz_xyyz[i] * fe_0 + tk_zzzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyyz[i] * fz_0;

        tk_yzzzz_xyyzz[i] = 2.0 * tk_zzzz_xyzz[i] * fe_0 + tk_zzzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyzz[i] * fz_0;

        tk_yzzzz_xyzzz[i] = tk_zzzz_xzzz[i] * fe_0 + tk_zzzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyzzz[i] * fz_0;

        tk_yzzzz_xzzzz[i] = tk_zzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xzzzz[i] * fz_0;

        tk_yzzzz_yyyyy[i] = 5.0 * tk_zzzz_yyyy[i] * fe_0 + tk_zzzz_yyyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyyy[i] * fz_0;

        tk_yzzzz_yyyyz[i] = 4.0 * tk_zzzz_yyyz[i] * fe_0 + tk_zzzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyyz[i] * fz_0;

        tk_yzzzz_yyyzz[i] = 3.0 * tk_zzzz_yyzz[i] * fe_0 + tk_zzzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyzz[i] * fz_0;

        tk_yzzzz_yyzzz[i] = 2.0 * tk_zzzz_yzzz[i] * fe_0 + tk_zzzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyzzz[i] * fz_0;

        tk_yzzzz_yzzzz[i] = tk_zzzz_zzzz[i] * fe_0 + tk_zzzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yzzzz[i] * fz_0;

        tk_yzzzz_zzzzz[i] = tk_zzzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_zzzzz[i] * fz_0;
    }

    // Set up 420-441 components of targeted buffer : HH

    auto tk_zzzzz_xxxxx = pbuffer.data(idx_kin_hh + 420);

    auto tk_zzzzz_xxxxy = pbuffer.data(idx_kin_hh + 421);

    auto tk_zzzzz_xxxxz = pbuffer.data(idx_kin_hh + 422);

    auto tk_zzzzz_xxxyy = pbuffer.data(idx_kin_hh + 423);

    auto tk_zzzzz_xxxyz = pbuffer.data(idx_kin_hh + 424);

    auto tk_zzzzz_xxxzz = pbuffer.data(idx_kin_hh + 425);

    auto tk_zzzzz_xxyyy = pbuffer.data(idx_kin_hh + 426);

    auto tk_zzzzz_xxyyz = pbuffer.data(idx_kin_hh + 427);

    auto tk_zzzzz_xxyzz = pbuffer.data(idx_kin_hh + 428);

    auto tk_zzzzz_xxzzz = pbuffer.data(idx_kin_hh + 429);

    auto tk_zzzzz_xyyyy = pbuffer.data(idx_kin_hh + 430);

    auto tk_zzzzz_xyyyz = pbuffer.data(idx_kin_hh + 431);

    auto tk_zzzzz_xyyzz = pbuffer.data(idx_kin_hh + 432);

    auto tk_zzzzz_xyzzz = pbuffer.data(idx_kin_hh + 433);

    auto tk_zzzzz_xzzzz = pbuffer.data(idx_kin_hh + 434);

    auto tk_zzzzz_yyyyy = pbuffer.data(idx_kin_hh + 435);

    auto tk_zzzzz_yyyyz = pbuffer.data(idx_kin_hh + 436);

    auto tk_zzzzz_yyyzz = pbuffer.data(idx_kin_hh + 437);

    auto tk_zzzzz_yyzzz = pbuffer.data(idx_kin_hh + 438);

    auto tk_zzzzz_yzzzz = pbuffer.data(idx_kin_hh + 439);

    auto tk_zzzzz_zzzzz = pbuffer.data(idx_kin_hh + 440);

#pragma omp simd aligned(pa_z,               \
                             tk_zzz_xxxxx,   \
                             tk_zzz_xxxxy,   \
                             tk_zzz_xxxxz,   \
                             tk_zzz_xxxyy,   \
                             tk_zzz_xxxyz,   \
                             tk_zzz_xxxzz,   \
                             tk_zzz_xxyyy,   \
                             tk_zzz_xxyyz,   \
                             tk_zzz_xxyzz,   \
                             tk_zzz_xxzzz,   \
                             tk_zzz_xyyyy,   \
                             tk_zzz_xyyyz,   \
                             tk_zzz_xyyzz,   \
                             tk_zzz_xyzzz,   \
                             tk_zzz_xzzzz,   \
                             tk_zzz_yyyyy,   \
                             tk_zzz_yyyyz,   \
                             tk_zzz_yyyzz,   \
                             tk_zzz_yyzzz,   \
                             tk_zzz_yzzzz,   \
                             tk_zzz_zzzzz,   \
                             tk_zzzz_xxxx,   \
                             tk_zzzz_xxxxx,  \
                             tk_zzzz_xxxxy,  \
                             tk_zzzz_xxxxz,  \
                             tk_zzzz_xxxy,   \
                             tk_zzzz_xxxyy,  \
                             tk_zzzz_xxxyz,  \
                             tk_zzzz_xxxz,   \
                             tk_zzzz_xxxzz,  \
                             tk_zzzz_xxyy,   \
                             tk_zzzz_xxyyy,  \
                             tk_zzzz_xxyyz,  \
                             tk_zzzz_xxyz,   \
                             tk_zzzz_xxyzz,  \
                             tk_zzzz_xxzz,   \
                             tk_zzzz_xxzzz,  \
                             tk_zzzz_xyyy,   \
                             tk_zzzz_xyyyy,  \
                             tk_zzzz_xyyyz,  \
                             tk_zzzz_xyyz,   \
                             tk_zzzz_xyyzz,  \
                             tk_zzzz_xyzz,   \
                             tk_zzzz_xyzzz,  \
                             tk_zzzz_xzzz,   \
                             tk_zzzz_xzzzz,  \
                             tk_zzzz_yyyy,   \
                             tk_zzzz_yyyyy,  \
                             tk_zzzz_yyyyz,  \
                             tk_zzzz_yyyz,   \
                             tk_zzzz_yyyzz,  \
                             tk_zzzz_yyzz,   \
                             tk_zzzz_yyzzz,  \
                             tk_zzzz_yzzz,   \
                             tk_zzzz_yzzzz,  \
                             tk_zzzz_zzzz,   \
                             tk_zzzz_zzzzz,  \
                             tk_zzzzz_xxxxx, \
                             tk_zzzzz_xxxxy, \
                             tk_zzzzz_xxxxz, \
                             tk_zzzzz_xxxyy, \
                             tk_zzzzz_xxxyz, \
                             tk_zzzzz_xxxzz, \
                             tk_zzzzz_xxyyy, \
                             tk_zzzzz_xxyyz, \
                             tk_zzzzz_xxyzz, \
                             tk_zzzzz_xxzzz, \
                             tk_zzzzz_xyyyy, \
                             tk_zzzzz_xyyyz, \
                             tk_zzzzz_xyyzz, \
                             tk_zzzzz_xyzzz, \
                             tk_zzzzz_xzzzz, \
                             tk_zzzzz_yyyyy, \
                             tk_zzzzz_yyyyz, \
                             tk_zzzzz_yyyzz, \
                             tk_zzzzz_yyzzz, \
                             tk_zzzzz_yzzzz, \
                             tk_zzzzz_zzzzz, \
                             ts_zzz_xxxxx,   \
                             ts_zzz_xxxxy,   \
                             ts_zzz_xxxxz,   \
                             ts_zzz_xxxyy,   \
                             ts_zzz_xxxyz,   \
                             ts_zzz_xxxzz,   \
                             ts_zzz_xxyyy,   \
                             ts_zzz_xxyyz,   \
                             ts_zzz_xxyzz,   \
                             ts_zzz_xxzzz,   \
                             ts_zzz_xyyyy,   \
                             ts_zzz_xyyyz,   \
                             ts_zzz_xyyzz,   \
                             ts_zzz_xyzzz,   \
                             ts_zzz_xzzzz,   \
                             ts_zzz_yyyyy,   \
                             ts_zzz_yyyyz,   \
                             ts_zzz_yyyzz,   \
                             ts_zzz_yyzzz,   \
                             ts_zzz_yzzzz,   \
                             ts_zzz_zzzzz,   \
                             ts_zzzzz_xxxxx, \
                             ts_zzzzz_xxxxy, \
                             ts_zzzzz_xxxxz, \
                             ts_zzzzz_xxxyy, \
                             ts_zzzzz_xxxyz, \
                             ts_zzzzz_xxxzz, \
                             ts_zzzzz_xxyyy, \
                             ts_zzzzz_xxyyz, \
                             ts_zzzzz_xxyzz, \
                             ts_zzzzz_xxzzz, \
                             ts_zzzzz_xyyyy, \
                             ts_zzzzz_xyyyz, \
                             ts_zzzzz_xyyzz, \
                             ts_zzzzz_xyzzz, \
                             ts_zzzzz_xzzzz, \
                             ts_zzzzz_yyyyy, \
                             ts_zzzzz_yyyyz, \
                             ts_zzzzz_yyyzz, \
                             ts_zzzzz_yyzzz, \
                             ts_zzzzz_yzzzz, \
                             ts_zzzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzz_xxxxx[i] =
            -8.0 * ts_zzz_xxxxx[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxx[i] * fe_0 + tk_zzzz_xxxxx[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxx[i] * fz_0;

        tk_zzzzz_xxxxy[i] =
            -8.0 * ts_zzz_xxxxy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxy[i] * fe_0 + tk_zzzz_xxxxy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxy[i] * fz_0;

        tk_zzzzz_xxxxz[i] = -8.0 * ts_zzz_xxxxz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxz[i] * fe_0 + tk_zzzz_xxxx[i] * fe_0 +
                            tk_zzzz_xxxxz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxz[i] * fz_0;

        tk_zzzzz_xxxyy[i] =
            -8.0 * ts_zzz_xxxyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxyy[i] * fe_0 + tk_zzzz_xxxyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxyy[i] * fz_0;

        tk_zzzzz_xxxyz[i] = -8.0 * ts_zzz_xxxyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxyz[i] * fe_0 + tk_zzzz_xxxy[i] * fe_0 +
                            tk_zzzz_xxxyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxyz[i] * fz_0;

        tk_zzzzz_xxxzz[i] = -8.0 * ts_zzz_xxxzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxzz[i] * fe_0 + 2.0 * tk_zzzz_xxxz[i] * fe_0 +
                            tk_zzzz_xxxzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxzz[i] * fz_0;

        tk_zzzzz_xxyyy[i] =
            -8.0 * ts_zzz_xxyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyyy[i] * fe_0 + tk_zzzz_xxyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyyy[i] * fz_0;

        tk_zzzzz_xxyyz[i] = -8.0 * ts_zzz_xxyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyyz[i] * fe_0 + tk_zzzz_xxyy[i] * fe_0 +
                            tk_zzzz_xxyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyyz[i] * fz_0;

        tk_zzzzz_xxyzz[i] = -8.0 * ts_zzz_xxyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyzz[i] * fe_0 + 2.0 * tk_zzzz_xxyz[i] * fe_0 +
                            tk_zzzz_xxyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyzz[i] * fz_0;

        tk_zzzzz_xxzzz[i] = -8.0 * ts_zzz_xxzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxzzz[i] * fe_0 + 3.0 * tk_zzzz_xxzz[i] * fe_0 +
                            tk_zzzz_xxzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxzzz[i] * fz_0;

        tk_zzzzz_xyyyy[i] =
            -8.0 * ts_zzz_xyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyyy[i] * fe_0 + tk_zzzz_xyyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyyy[i] * fz_0;

        tk_zzzzz_xyyyz[i] = -8.0 * ts_zzz_xyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyyz[i] * fe_0 + tk_zzzz_xyyy[i] * fe_0 +
                            tk_zzzz_xyyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyyz[i] * fz_0;

        tk_zzzzz_xyyzz[i] = -8.0 * ts_zzz_xyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyzz[i] * fe_0 + 2.0 * tk_zzzz_xyyz[i] * fe_0 +
                            tk_zzzz_xyyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyzz[i] * fz_0;

        tk_zzzzz_xyzzz[i] = -8.0 * ts_zzz_xyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyzzz[i] * fe_0 + 3.0 * tk_zzzz_xyzz[i] * fe_0 +
                            tk_zzzz_xyzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyzzz[i] * fz_0;

        tk_zzzzz_xzzzz[i] = -8.0 * ts_zzz_xzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xzzzz[i] * fe_0 + 4.0 * tk_zzzz_xzzz[i] * fe_0 +
                            tk_zzzz_xzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xzzzz[i] * fz_0;

        tk_zzzzz_yyyyy[i] =
            -8.0 * ts_zzz_yyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyyy[i] * fe_0 + tk_zzzz_yyyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyyy[i] * fz_0;

        tk_zzzzz_yyyyz[i] = -8.0 * ts_zzz_yyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyyz[i] * fe_0 + tk_zzzz_yyyy[i] * fe_0 +
                            tk_zzzz_yyyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyyz[i] * fz_0;

        tk_zzzzz_yyyzz[i] = -8.0 * ts_zzz_yyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyzz[i] * fe_0 + 2.0 * tk_zzzz_yyyz[i] * fe_0 +
                            tk_zzzz_yyyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyzz[i] * fz_0;

        tk_zzzzz_yyzzz[i] = -8.0 * ts_zzz_yyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyzzz[i] * fe_0 + 3.0 * tk_zzzz_yyzz[i] * fe_0 +
                            tk_zzzz_yyzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyzzz[i] * fz_0;

        tk_zzzzz_yzzzz[i] = -8.0 * ts_zzz_yzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yzzzz[i] * fe_0 + 4.0 * tk_zzzz_yzzz[i] * fe_0 +
                            tk_zzzz_yzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yzzzz[i] * fz_0;

        tk_zzzzz_zzzzz[i] = -8.0 * ts_zzz_zzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_zzzzz[i] * fe_0 + 5.0 * tk_zzzz_zzzz[i] * fe_0 +
                            tk_zzzz_zzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_zzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
