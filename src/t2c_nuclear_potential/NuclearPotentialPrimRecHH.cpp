#include "NuclearPotentialPrimRecHH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hh(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hh,
                               const size_t              idx_npot_0_fh,
                               const size_t              idx_npot_1_fh,
                               const size_t              idx_npot_0_gg,
                               const size_t              idx_npot_1_gg,
                               const size_t              idx_npot_0_gh,
                               const size_t              idx_npot_1_gh,
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

    // Set up components of auxiliary buffer : FH

    auto ta_xxx_xxxxx_0 = pbuffer.data(idx_npot_0_fh);

    auto ta_xxx_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 1);

    auto ta_xxx_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 2);

    auto ta_xxx_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 3);

    auto ta_xxx_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 4);

    auto ta_xxx_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 5);

    auto ta_xxx_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 6);

    auto ta_xxx_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 7);

    auto ta_xxx_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 8);

    auto ta_xxx_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 9);

    auto ta_xxx_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 10);

    auto ta_xxx_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 11);

    auto ta_xxx_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 12);

    auto ta_xxx_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 13);

    auto ta_xxx_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 14);

    auto ta_xxx_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 15);

    auto ta_xxx_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 16);

    auto ta_xxx_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 17);

    auto ta_xxx_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 18);

    auto ta_xxx_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 19);

    auto ta_xxx_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 20);

    auto ta_xxy_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 21);

    auto ta_xxy_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 23);

    auto ta_xxy_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 26);

    auto ta_xxy_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 30);

    auto ta_xxy_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 35);

    auto ta_xxy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 36);

    auto ta_xxy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 37);

    auto ta_xxy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 38);

    auto ta_xxy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 39);

    auto ta_xxy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 40);

    auto ta_xxz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 42);

    auto ta_xxz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 43);

    auto ta_xxz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 44);

    auto ta_xxz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 45);

    auto ta_xxz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 47);

    auto ta_xxz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 48);

    auto ta_xxz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 51);

    auto ta_xxz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 52);

    auto ta_xxz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 56);

    auto ta_xxz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 58);

    auto ta_xxz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 59);

    auto ta_xxz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 60);

    auto ta_xxz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 61);

    auto ta_xxz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 62);

    auto ta_xyy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 64);

    auto ta_xyy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 66);

    auto ta_xyy_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 67);

    auto ta_xyy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 69);

    auto ta_xyy_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 70);

    auto ta_xyy_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 71);

    auto ta_xyy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 73);

    auto ta_xyy_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 74);

    auto ta_xyy_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 75);

    auto ta_xyy_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 76);

    auto ta_xyy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 78);

    auto ta_xyy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 79);

    auto ta_xyy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 80);

    auto ta_xyy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 81);

    auto ta_xyy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 82);

    auto ta_xyy_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 83);

    auto ta_xyz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 100);

    auto ta_xyz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 101);

    auto ta_xyz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 102);

    auto ta_xyz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 103);

    auto ta_xzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 107);

    auto ta_xzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 109);

    auto ta_xzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 110);

    auto ta_xzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 112);

    auto ta_xzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 113);

    auto ta_xzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 114);

    auto ta_xzz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 116);

    auto ta_xzz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 117);

    auto ta_xzz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 118);

    auto ta_xzz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 119);

    auto ta_xzz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 120);

    auto ta_xzz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 121);

    auto ta_xzz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 122);

    auto ta_xzz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 123);

    auto ta_xzz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 124);

    auto ta_xzz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 125);

    auto ta_yyy_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 126);

    auto ta_yyy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 127);

    auto ta_yyy_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 128);

    auto ta_yyy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 129);

    auto ta_yyy_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 130);

    auto ta_yyy_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 131);

    auto ta_yyy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 132);

    auto ta_yyy_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 133);

    auto ta_yyy_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 134);

    auto ta_yyy_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 135);

    auto ta_yyy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 136);

    auto ta_yyy_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 137);

    auto ta_yyy_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 138);

    auto ta_yyy_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 139);

    auto ta_yyy_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 140);

    auto ta_yyy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 141);

    auto ta_yyy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 142);

    auto ta_yyy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 143);

    auto ta_yyy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 144);

    auto ta_yyy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 145);

    auto ta_yyy_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 146);

    auto ta_yyz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 148);

    auto ta_yyz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 149);

    auto ta_yyz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 150);

    auto ta_yyz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 152);

    auto ta_yyz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 153);

    auto ta_yyz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 156);

    auto ta_yyz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 157);

    auto ta_yyz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 161);

    auto ta_yyz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 162);

    auto ta_yyz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 163);

    auto ta_yyz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 164);

    auto ta_yyz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 165);

    auto ta_yyz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 166);

    auto ta_yyz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 167);

    auto ta_yzz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 168);

    auto ta_yzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 170);

    auto ta_yzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 172);

    auto ta_yzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 173);

    auto ta_yzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 175);

    auto ta_yzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 176);

    auto ta_yzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 177);

    auto ta_yzz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 179);

    auto ta_yzz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 180);

    auto ta_yzz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 181);

    auto ta_yzz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 182);

    auto ta_yzz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 183);

    auto ta_yzz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 184);

    auto ta_yzz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 185);

    auto ta_yzz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 186);

    auto ta_yzz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 187);

    auto ta_yzz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 188);

    auto ta_zzz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 189);

    auto ta_zzz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 190);

    auto ta_zzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 191);

    auto ta_zzz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 192);

    auto ta_zzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 193);

    auto ta_zzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 194);

    auto ta_zzz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 195);

    auto ta_zzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 196);

    auto ta_zzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 197);

    auto ta_zzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 198);

    auto ta_zzz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 199);

    auto ta_zzz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 200);

    auto ta_zzz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 201);

    auto ta_zzz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 202);

    auto ta_zzz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 203);

    auto ta_zzz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 204);

    auto ta_zzz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 205);

    auto ta_zzz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 206);

    auto ta_zzz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 207);

    auto ta_zzz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 208);

    auto ta_zzz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 209);

    // Set up components of auxiliary buffer : FH

    auto ta_xxx_xxxxx_1 = pbuffer.data(idx_npot_1_fh);

    auto ta_xxx_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 1);

    auto ta_xxx_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 2);

    auto ta_xxx_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 3);

    auto ta_xxx_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 4);

    auto ta_xxx_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 5);

    auto ta_xxx_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 6);

    auto ta_xxx_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 7);

    auto ta_xxx_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 8);

    auto ta_xxx_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 9);

    auto ta_xxx_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 10);

    auto ta_xxx_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 11);

    auto ta_xxx_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 12);

    auto ta_xxx_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 13);

    auto ta_xxx_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 14);

    auto ta_xxx_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 15);

    auto ta_xxx_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 16);

    auto ta_xxx_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 17);

    auto ta_xxx_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 18);

    auto ta_xxx_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 19);

    auto ta_xxx_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 20);

    auto ta_xxy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 21);

    auto ta_xxy_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 23);

    auto ta_xxy_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 26);

    auto ta_xxy_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 30);

    auto ta_xxy_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 35);

    auto ta_xxy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 36);

    auto ta_xxy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 37);

    auto ta_xxy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 38);

    auto ta_xxy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 39);

    auto ta_xxy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 40);

    auto ta_xxz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 42);

    auto ta_xxz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 43);

    auto ta_xxz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 44);

    auto ta_xxz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 45);

    auto ta_xxz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 47);

    auto ta_xxz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 48);

    auto ta_xxz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 51);

    auto ta_xxz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 52);

    auto ta_xxz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 56);

    auto ta_xxz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 58);

    auto ta_xxz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 59);

    auto ta_xxz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 60);

    auto ta_xxz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 61);

    auto ta_xxz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 62);

    auto ta_xyy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 64);

    auto ta_xyy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 66);

    auto ta_xyy_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 67);

    auto ta_xyy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 69);

    auto ta_xyy_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 70);

    auto ta_xyy_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 71);

    auto ta_xyy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 73);

    auto ta_xyy_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 74);

    auto ta_xyy_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 75);

    auto ta_xyy_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 76);

    auto ta_xyy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 78);

    auto ta_xyy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 79);

    auto ta_xyy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 80);

    auto ta_xyy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 81);

    auto ta_xyy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 82);

    auto ta_xyy_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 83);

    auto ta_xyz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 100);

    auto ta_xyz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 101);

    auto ta_xyz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 102);

    auto ta_xyz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 103);

    auto ta_xzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 107);

    auto ta_xzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 109);

    auto ta_xzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 110);

    auto ta_xzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 112);

    auto ta_xzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 113);

    auto ta_xzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 114);

    auto ta_xzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 116);

    auto ta_xzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 117);

    auto ta_xzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 118);

    auto ta_xzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 119);

    auto ta_xzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 120);

    auto ta_xzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 121);

    auto ta_xzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 122);

    auto ta_xzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 123);

    auto ta_xzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 124);

    auto ta_xzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 125);

    auto ta_yyy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 126);

    auto ta_yyy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 127);

    auto ta_yyy_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 128);

    auto ta_yyy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 129);

    auto ta_yyy_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 130);

    auto ta_yyy_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 131);

    auto ta_yyy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 132);

    auto ta_yyy_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 133);

    auto ta_yyy_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 134);

    auto ta_yyy_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 135);

    auto ta_yyy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 136);

    auto ta_yyy_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 137);

    auto ta_yyy_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 138);

    auto ta_yyy_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 139);

    auto ta_yyy_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 140);

    auto ta_yyy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 141);

    auto ta_yyy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 142);

    auto ta_yyy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 143);

    auto ta_yyy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 144);

    auto ta_yyy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 145);

    auto ta_yyy_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 146);

    auto ta_yyz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 148);

    auto ta_yyz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 149);

    auto ta_yyz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 150);

    auto ta_yyz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 152);

    auto ta_yyz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 153);

    auto ta_yyz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 156);

    auto ta_yyz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 157);

    auto ta_yyz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 161);

    auto ta_yyz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 162);

    auto ta_yyz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 163);

    auto ta_yyz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 164);

    auto ta_yyz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 165);

    auto ta_yyz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 166);

    auto ta_yyz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 167);

    auto ta_yzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 168);

    auto ta_yzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 170);

    auto ta_yzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 172);

    auto ta_yzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 173);

    auto ta_yzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 175);

    auto ta_yzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 176);

    auto ta_yzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 177);

    auto ta_yzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 179);

    auto ta_yzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 180);

    auto ta_yzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 181);

    auto ta_yzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 182);

    auto ta_yzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 183);

    auto ta_yzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 184);

    auto ta_yzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 185);

    auto ta_yzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 186);

    auto ta_yzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 187);

    auto ta_yzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 188);

    auto ta_zzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 189);

    auto ta_zzz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 190);

    auto ta_zzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 191);

    auto ta_zzz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 192);

    auto ta_zzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 193);

    auto ta_zzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 194);

    auto ta_zzz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 195);

    auto ta_zzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 196);

    auto ta_zzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 197);

    auto ta_zzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 198);

    auto ta_zzz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 199);

    auto ta_zzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 200);

    auto ta_zzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 201);

    auto ta_zzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 202);

    auto ta_zzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 203);

    auto ta_zzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 204);

    auto ta_zzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 205);

    auto ta_zzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 206);

    auto ta_zzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 207);

    auto ta_zzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 208);

    auto ta_zzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 209);

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

    auto ta_xxxz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 32);

    auto ta_xxxz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 34);

    auto ta_xxxz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 35);

    auto ta_xxxz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 37);

    auto ta_xxxz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 38);

    auto ta_xxxz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 39);

    auto ta_xxyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 46);

    auto ta_xxyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 48);

    auto ta_xxyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 49);

    auto ta_xxyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 51);

    auto ta_xxyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 52);

    auto ta_xxyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 53);

    auto ta_xxyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 55);

    auto ta_xxyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 56);

    auto ta_xxyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 57);

    auto ta_xxyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 58);

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

    auto ta_xzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 137);

    auto ta_xzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 139);

    auto ta_xzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 140);

    auto ta_xzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 142);

    auto ta_xzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 143);

    auto ta_xzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 144);

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

    auto ta_yyyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 167);

    auto ta_yyyz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 169);

    auto ta_yyyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 170);

    auto ta_yyyz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 172);

    auto ta_yyyz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 173);

    auto ta_yyyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 174);

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

    auto ta_yzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 196);

    auto ta_yzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 197);

    auto ta_yzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 198);

    auto ta_yzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 199);

    auto ta_yzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 200);

    auto ta_yzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 201);

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

    auto ta_xxxz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 32);

    auto ta_xxxz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 34);

    auto ta_xxxz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 35);

    auto ta_xxxz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 37);

    auto ta_xxxz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 38);

    auto ta_xxxz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 39);

    auto ta_xxyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 46);

    auto ta_xxyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 48);

    auto ta_xxyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 49);

    auto ta_xxyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 51);

    auto ta_xxyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 52);

    auto ta_xxyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 53);

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

    auto ta_xzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 137);

    auto ta_xzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 139);

    auto ta_xzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 140);

    auto ta_xzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 142);

    auto ta_xzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 143);

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

    auto ta_yyyz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 167);

    auto ta_yyyz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 169);

    auto ta_yyyz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 170);

    auto ta_yyyz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 172);

    auto ta_yyyz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 173);

    auto ta_yyyz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 174);

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

    auto ta_yzzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 196);

    auto ta_yzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 197);

    auto ta_yzzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 198);

    auto ta_yzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 199);

    auto ta_yzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 200);

    auto ta_yzzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 201);

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

    // Set up components of auxiliary buffer : GH

    auto ta_xxxx_xxxxx_0 = pbuffer.data(idx_npot_0_gh);

    auto ta_xxxx_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 1);

    auto ta_xxxx_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 2);

    auto ta_xxxx_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 3);

    auto ta_xxxx_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 4);

    auto ta_xxxx_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 5);

    auto ta_xxxx_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 6);

    auto ta_xxxx_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 7);

    auto ta_xxxx_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 8);

    auto ta_xxxx_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 9);

    auto ta_xxxx_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 10);

    auto ta_xxxx_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 11);

    auto ta_xxxx_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 12);

    auto ta_xxxx_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 13);

    auto ta_xxxx_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 14);

    auto ta_xxxx_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 15);

    auto ta_xxxx_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 16);

    auto ta_xxxx_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 17);

    auto ta_xxxx_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 18);

    auto ta_xxxx_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 19);

    auto ta_xxxx_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 20);

    auto ta_xxxy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 21);

    auto ta_xxxy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 22);

    auto ta_xxxy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 23);

    auto ta_xxxy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 24);

    auto ta_xxxy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 26);

    auto ta_xxxy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 27);

    auto ta_xxxy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 30);

    auto ta_xxxy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 31);

    auto ta_xxxy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 35);

    auto ta_xxxy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 36);

    auto ta_xxxy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 37);

    auto ta_xxxy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 38);

    auto ta_xxxy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 39);

    auto ta_xxxy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 40);

    auto ta_xxxz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 42);

    auto ta_xxxz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 43);

    auto ta_xxxz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 44);

    auto ta_xxxz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 45);

    auto ta_xxxz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 46);

    auto ta_xxxz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 47);

    auto ta_xxxz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 48);

    auto ta_xxxz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 49);

    auto ta_xxxz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 50);

    auto ta_xxxz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 51);

    auto ta_xxxz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 52);

    auto ta_xxxz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 53);

    auto ta_xxxz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 54);

    auto ta_xxxz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 55);

    auto ta_xxxz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 56);

    auto ta_xxxz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 58);

    auto ta_xxxz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 59);

    auto ta_xxxz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 60);

    auto ta_xxxz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 61);

    auto ta_xxxz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 62);

    auto ta_xxyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 63);

    auto ta_xxyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 64);

    auto ta_xxyy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 65);

    auto ta_xxyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 66);

    auto ta_xxyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 67);

    auto ta_xxyy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 68);

    auto ta_xxyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 69);

    auto ta_xxyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 70);

    auto ta_xxyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 71);

    auto ta_xxyy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 72);

    auto ta_xxyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 73);

    auto ta_xxyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 74);

    auto ta_xxyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 75);

    auto ta_xxyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 76);

    auto ta_xxyy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 77);

    auto ta_xxyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 78);

    auto ta_xxyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 79);

    auto ta_xxyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 80);

    auto ta_xxyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 81);

    auto ta_xxyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 82);

    auto ta_xxyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 83);

    auto ta_xxyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 86);

    auto ta_xxyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 89);

    auto ta_xxyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 93);

    auto ta_xxyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 98);

    auto ta_xxyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 100);

    auto ta_xxyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 101);

    auto ta_xxyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 102);

    auto ta_xxyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 103);

    auto ta_xxzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 105);

    auto ta_xxzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 106);

    auto ta_xxzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 107);

    auto ta_xxzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 108);

    auto ta_xxzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 109);

    auto ta_xxzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 110);

    auto ta_xxzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 111);

    auto ta_xxzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 112);

    auto ta_xxzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 113);

    auto ta_xxzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 114);

    auto ta_xxzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 115);

    auto ta_xxzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 116);

    auto ta_xxzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 117);

    auto ta_xxzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 118);

    auto ta_xxzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 119);

    auto ta_xxzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 120);

    auto ta_xxzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 121);

    auto ta_xxzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 122);

    auto ta_xxzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 123);

    auto ta_xxzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 124);

    auto ta_xxzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 125);

    auto ta_xyyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 126);

    auto ta_xyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 127);

    auto ta_xyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 129);

    auto ta_xyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 130);

    auto ta_xyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 132);

    auto ta_xyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 133);

    auto ta_xyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 134);

    auto ta_xyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 136);

    auto ta_xyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 137);

    auto ta_xyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 138);

    auto ta_xyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 139);

    auto ta_xyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 141);

    auto ta_xyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 142);

    auto ta_xyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 143);

    auto ta_xyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 144);

    auto ta_xyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 145);

    auto ta_xyyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 146);

    auto ta_xyyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 163);

    auto ta_xyyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 164);

    auto ta_xyyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 165);

    auto ta_xyyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 166);

    auto ta_xyyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 167);

    auto ta_xyzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 183);

    auto ta_xyzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 184);

    auto ta_xyzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 185);

    auto ta_xyzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 186);

    auto ta_xyzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 187);

    auto ta_xzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 189);

    auto ta_xzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 191);

    auto ta_xzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 193);

    auto ta_xzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 194);

    auto ta_xzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 196);

    auto ta_xzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 197);

    auto ta_xzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 198);

    auto ta_xzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 200);

    auto ta_xzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 201);

    auto ta_xzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 202);

    auto ta_xzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 203);

    auto ta_xzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 204);

    auto ta_xzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 205);

    auto ta_xzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 206);

    auto ta_xzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 207);

    auto ta_xzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 208);

    auto ta_xzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 209);

    auto ta_yyyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 210);

    auto ta_yyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 211);

    auto ta_yyyy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 212);

    auto ta_yyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 213);

    auto ta_yyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 214);

    auto ta_yyyy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 215);

    auto ta_yyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 216);

    auto ta_yyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 217);

    auto ta_yyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 218);

    auto ta_yyyy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 219);

    auto ta_yyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 220);

    auto ta_yyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 221);

    auto ta_yyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 222);

    auto ta_yyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 223);

    auto ta_yyyy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 224);

    auto ta_yyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 225);

    auto ta_yyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 226);

    auto ta_yyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 227);

    auto ta_yyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 228);

    auto ta_yyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 229);

    auto ta_yyyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 230);

    auto ta_yyyz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 232);

    auto ta_yyyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 233);

    auto ta_yyyz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 234);

    auto ta_yyyz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 235);

    auto ta_yyyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 236);

    auto ta_yyyz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 237);

    auto ta_yyyz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 238);

    auto ta_yyyz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 239);

    auto ta_yyyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 240);

    auto ta_yyyz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 241);

    auto ta_yyyz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 242);

    auto ta_yyyz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 243);

    auto ta_yyyz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 244);

    auto ta_yyyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 245);

    auto ta_yyyz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 246);

    auto ta_yyyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 247);

    auto ta_yyyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 248);

    auto ta_yyyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 249);

    auto ta_yyyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 250);

    auto ta_yyyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 251);

    auto ta_yyzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 252);

    auto ta_yyzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 253);

    auto ta_yyzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 254);

    auto ta_yyzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 255);

    auto ta_yyzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 256);

    auto ta_yyzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 257);

    auto ta_yyzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 258);

    auto ta_yyzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 259);

    auto ta_yyzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 260);

    auto ta_yyzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 261);

    auto ta_yyzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 262);

    auto ta_yyzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 263);

    auto ta_yyzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 264);

    auto ta_yyzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 265);

    auto ta_yyzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 266);

    auto ta_yyzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 267);

    auto ta_yyzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 268);

    auto ta_yyzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 269);

    auto ta_yyzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 270);

    auto ta_yyzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 271);

    auto ta_yyzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 272);

    auto ta_yzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 273);

    auto ta_yzzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 274);

    auto ta_yzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 275);

    auto ta_yzzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 276);

    auto ta_yzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 277);

    auto ta_yzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 278);

    auto ta_yzzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 279);

    auto ta_yzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 280);

    auto ta_yzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 281);

    auto ta_yzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 282);

    auto ta_yzzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 283);

    auto ta_yzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 284);

    auto ta_yzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 285);

    auto ta_yzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 286);

    auto ta_yzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 287);

    auto ta_yzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 288);

    auto ta_yzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 289);

    auto ta_yzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 290);

    auto ta_yzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 291);

    auto ta_yzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 292);

    auto ta_yzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 293);

    auto ta_zzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 294);

    auto ta_zzzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 295);

    auto ta_zzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 296);

    auto ta_zzzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 297);

    auto ta_zzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 298);

    auto ta_zzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 299);

    auto ta_zzzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 300);

    auto ta_zzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 301);

    auto ta_zzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 302);

    auto ta_zzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 303);

    auto ta_zzzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 304);

    auto ta_zzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 305);

    auto ta_zzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 306);

    auto ta_zzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 307);

    auto ta_zzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 308);

    auto ta_zzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 309);

    auto ta_zzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 310);

    auto ta_zzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 311);

    auto ta_zzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 312);

    auto ta_zzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 313);

    auto ta_zzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 314);

    // Set up components of auxiliary buffer : GH

    auto ta_xxxx_xxxxx_1 = pbuffer.data(idx_npot_1_gh);

    auto ta_xxxx_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 1);

    auto ta_xxxx_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 2);

    auto ta_xxxx_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 3);

    auto ta_xxxx_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 4);

    auto ta_xxxx_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 5);

    auto ta_xxxx_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 6);

    auto ta_xxxx_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 7);

    auto ta_xxxx_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 8);

    auto ta_xxxx_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 9);

    auto ta_xxxx_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 10);

    auto ta_xxxx_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 11);

    auto ta_xxxx_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 12);

    auto ta_xxxx_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 13);

    auto ta_xxxx_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 14);

    auto ta_xxxx_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 15);

    auto ta_xxxx_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 16);

    auto ta_xxxx_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 17);

    auto ta_xxxx_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 18);

    auto ta_xxxx_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 19);

    auto ta_xxxx_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 20);

    auto ta_xxxy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 21);

    auto ta_xxxy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 22);

    auto ta_xxxy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 23);

    auto ta_xxxy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 24);

    auto ta_xxxy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 26);

    auto ta_xxxy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 27);

    auto ta_xxxy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 30);

    auto ta_xxxy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 31);

    auto ta_xxxy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 35);

    auto ta_xxxy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 36);

    auto ta_xxxy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 37);

    auto ta_xxxy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 38);

    auto ta_xxxy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 39);

    auto ta_xxxy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 40);

    auto ta_xxxz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 42);

    auto ta_xxxz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 43);

    auto ta_xxxz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 44);

    auto ta_xxxz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 45);

    auto ta_xxxz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 46);

    auto ta_xxxz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 47);

    auto ta_xxxz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 48);

    auto ta_xxxz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 49);

    auto ta_xxxz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 50);

    auto ta_xxxz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 51);

    auto ta_xxxz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 52);

    auto ta_xxxz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 53);

    auto ta_xxxz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 54);

    auto ta_xxxz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 55);

    auto ta_xxxz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 56);

    auto ta_xxxz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 58);

    auto ta_xxxz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 59);

    auto ta_xxxz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 60);

    auto ta_xxxz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 61);

    auto ta_xxxz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 62);

    auto ta_xxyy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 63);

    auto ta_xxyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 64);

    auto ta_xxyy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 65);

    auto ta_xxyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 66);

    auto ta_xxyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 67);

    auto ta_xxyy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 68);

    auto ta_xxyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 69);

    auto ta_xxyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 70);

    auto ta_xxyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 71);

    auto ta_xxyy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 72);

    auto ta_xxyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 73);

    auto ta_xxyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 74);

    auto ta_xxyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 75);

    auto ta_xxyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 76);

    auto ta_xxyy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 77);

    auto ta_xxyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 78);

    auto ta_xxyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 79);

    auto ta_xxyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 80);

    auto ta_xxyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 81);

    auto ta_xxyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 82);

    auto ta_xxyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 83);

    auto ta_xxyz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 86);

    auto ta_xxyz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 89);

    auto ta_xxyz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 93);

    auto ta_xxyz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 98);

    auto ta_xxyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 100);

    auto ta_xxyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 101);

    auto ta_xxyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 102);

    auto ta_xxyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 103);

    auto ta_xxzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 105);

    auto ta_xxzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 106);

    auto ta_xxzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 107);

    auto ta_xxzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 108);

    auto ta_xxzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 109);

    auto ta_xxzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 110);

    auto ta_xxzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 111);

    auto ta_xxzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 112);

    auto ta_xxzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 113);

    auto ta_xxzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 114);

    auto ta_xxzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 115);

    auto ta_xxzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 116);

    auto ta_xxzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 117);

    auto ta_xxzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 118);

    auto ta_xxzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 119);

    auto ta_xxzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 120);

    auto ta_xxzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 121);

    auto ta_xxzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 122);

    auto ta_xxzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 123);

    auto ta_xxzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 124);

    auto ta_xxzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 125);

    auto ta_xyyy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 126);

    auto ta_xyyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 127);

    auto ta_xyyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 129);

    auto ta_xyyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 130);

    auto ta_xyyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 132);

    auto ta_xyyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 133);

    auto ta_xyyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 134);

    auto ta_xyyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 136);

    auto ta_xyyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 137);

    auto ta_xyyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 138);

    auto ta_xyyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 139);

    auto ta_xyyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 141);

    auto ta_xyyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 142);

    auto ta_xyyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 143);

    auto ta_xyyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 144);

    auto ta_xyyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 145);

    auto ta_xyyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 146);

    auto ta_xyyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 163);

    auto ta_xyyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 164);

    auto ta_xyyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 165);

    auto ta_xyyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 166);

    auto ta_xyyz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 167);

    auto ta_xyzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 183);

    auto ta_xyzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 184);

    auto ta_xyzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 185);

    auto ta_xyzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 186);

    auto ta_xyzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 187);

    auto ta_xzzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 189);

    auto ta_xzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 191);

    auto ta_xzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 193);

    auto ta_xzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 194);

    auto ta_xzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 196);

    auto ta_xzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 197);

    auto ta_xzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 198);

    auto ta_xzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 200);

    auto ta_xzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 201);

    auto ta_xzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 202);

    auto ta_xzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 203);

    auto ta_xzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 204);

    auto ta_xzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 205);

    auto ta_xzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 206);

    auto ta_xzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 207);

    auto ta_xzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 208);

    auto ta_xzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 209);

    auto ta_yyyy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 210);

    auto ta_yyyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 211);

    auto ta_yyyy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 212);

    auto ta_yyyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 213);

    auto ta_yyyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 214);

    auto ta_yyyy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 215);

    auto ta_yyyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 216);

    auto ta_yyyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 217);

    auto ta_yyyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 218);

    auto ta_yyyy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 219);

    auto ta_yyyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 220);

    auto ta_yyyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 221);

    auto ta_yyyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 222);

    auto ta_yyyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 223);

    auto ta_yyyy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 224);

    auto ta_yyyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 225);

    auto ta_yyyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 226);

    auto ta_yyyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 227);

    auto ta_yyyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 228);

    auto ta_yyyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 229);

    auto ta_yyyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 230);

    auto ta_yyyz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 232);

    auto ta_yyyz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 233);

    auto ta_yyyz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 234);

    auto ta_yyyz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 235);

    auto ta_yyyz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 236);

    auto ta_yyyz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 237);

    auto ta_yyyz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 238);

    auto ta_yyyz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 239);

    auto ta_yyyz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 240);

    auto ta_yyyz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 241);

    auto ta_yyyz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 242);

    auto ta_yyyz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 243);

    auto ta_yyyz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 244);

    auto ta_yyyz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 245);

    auto ta_yyyz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 246);

    auto ta_yyyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 247);

    auto ta_yyyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 248);

    auto ta_yyyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 249);

    auto ta_yyyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 250);

    auto ta_yyyz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 251);

    auto ta_yyzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 252);

    auto ta_yyzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 253);

    auto ta_yyzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 254);

    auto ta_yyzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 255);

    auto ta_yyzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 256);

    auto ta_yyzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 257);

    auto ta_yyzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 258);

    auto ta_yyzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 259);

    auto ta_yyzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 260);

    auto ta_yyzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 261);

    auto ta_yyzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 262);

    auto ta_yyzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 263);

    auto ta_yyzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 264);

    auto ta_yyzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 265);

    auto ta_yyzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 266);

    auto ta_yyzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 267);

    auto ta_yyzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 268);

    auto ta_yyzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 269);

    auto ta_yyzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 270);

    auto ta_yyzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 271);

    auto ta_yyzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 272);

    auto ta_yzzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 273);

    auto ta_yzzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 274);

    auto ta_yzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 275);

    auto ta_yzzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 276);

    auto ta_yzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 277);

    auto ta_yzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 278);

    auto ta_yzzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 279);

    auto ta_yzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 280);

    auto ta_yzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 281);

    auto ta_yzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 282);

    auto ta_yzzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 283);

    auto ta_yzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 284);

    auto ta_yzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 285);

    auto ta_yzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 286);

    auto ta_yzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 287);

    auto ta_yzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 288);

    auto ta_yzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 289);

    auto ta_yzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 290);

    auto ta_yzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 291);

    auto ta_yzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 292);

    auto ta_yzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 293);

    auto ta_zzzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 294);

    auto ta_zzzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 295);

    auto ta_zzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 296);

    auto ta_zzzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 297);

    auto ta_zzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 298);

    auto ta_zzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 299);

    auto ta_zzzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 300);

    auto ta_zzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 301);

    auto ta_zzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 302);

    auto ta_zzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 303);

    auto ta_zzzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 304);

    auto ta_zzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 305);

    auto ta_zzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 306);

    auto ta_zzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 307);

    auto ta_zzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 308);

    auto ta_zzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 309);

    auto ta_zzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 310);

    auto ta_zzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 311);

    auto ta_zzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 312);

    auto ta_zzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 313);

    auto ta_zzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 314);

    // Set up 0-21 components of targeted buffer : HH

    auto ta_xxxxx_xxxxx_0 = pbuffer.data(idx_npot_0_hh);

    auto ta_xxxxx_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 1);

    auto ta_xxxxx_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 2);

    auto ta_xxxxx_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 3);

    auto ta_xxxxx_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 4);

    auto ta_xxxxx_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 5);

    auto ta_xxxxx_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 6);

    auto ta_xxxxx_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 7);

    auto ta_xxxxx_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 8);

    auto ta_xxxxx_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 9);

    auto ta_xxxxx_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 10);

    auto ta_xxxxx_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 11);

    auto ta_xxxxx_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 12);

    auto ta_xxxxx_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 13);

    auto ta_xxxxx_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 14);

    auto ta_xxxxx_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 15);

    auto ta_xxxxx_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 16);

    auto ta_xxxxx_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 17);

    auto ta_xxxxx_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 18);

    auto ta_xxxxx_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 19);

    auto ta_xxxxx_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 20);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xxx_xxxxx_0,   \
                             ta_xxx_xxxxx_1,   \
                             ta_xxx_xxxxy_0,   \
                             ta_xxx_xxxxy_1,   \
                             ta_xxx_xxxxz_0,   \
                             ta_xxx_xxxxz_1,   \
                             ta_xxx_xxxyy_0,   \
                             ta_xxx_xxxyy_1,   \
                             ta_xxx_xxxyz_0,   \
                             ta_xxx_xxxyz_1,   \
                             ta_xxx_xxxzz_0,   \
                             ta_xxx_xxxzz_1,   \
                             ta_xxx_xxyyy_0,   \
                             ta_xxx_xxyyy_1,   \
                             ta_xxx_xxyyz_0,   \
                             ta_xxx_xxyyz_1,   \
                             ta_xxx_xxyzz_0,   \
                             ta_xxx_xxyzz_1,   \
                             ta_xxx_xxzzz_0,   \
                             ta_xxx_xxzzz_1,   \
                             ta_xxx_xyyyy_0,   \
                             ta_xxx_xyyyy_1,   \
                             ta_xxx_xyyyz_0,   \
                             ta_xxx_xyyyz_1,   \
                             ta_xxx_xyyzz_0,   \
                             ta_xxx_xyyzz_1,   \
                             ta_xxx_xyzzz_0,   \
                             ta_xxx_xyzzz_1,   \
                             ta_xxx_xzzzz_0,   \
                             ta_xxx_xzzzz_1,   \
                             ta_xxx_yyyyy_0,   \
                             ta_xxx_yyyyy_1,   \
                             ta_xxx_yyyyz_0,   \
                             ta_xxx_yyyyz_1,   \
                             ta_xxx_yyyzz_0,   \
                             ta_xxx_yyyzz_1,   \
                             ta_xxx_yyzzz_0,   \
                             ta_xxx_yyzzz_1,   \
                             ta_xxx_yzzzz_0,   \
                             ta_xxx_yzzzz_1,   \
                             ta_xxx_zzzzz_0,   \
                             ta_xxx_zzzzz_1,   \
                             ta_xxxx_xxxx_0,   \
                             ta_xxxx_xxxx_1,   \
                             ta_xxxx_xxxxx_0,  \
                             ta_xxxx_xxxxx_1,  \
                             ta_xxxx_xxxxy_0,  \
                             ta_xxxx_xxxxy_1,  \
                             ta_xxxx_xxxxz_0,  \
                             ta_xxxx_xxxxz_1,  \
                             ta_xxxx_xxxy_0,   \
                             ta_xxxx_xxxy_1,   \
                             ta_xxxx_xxxyy_0,  \
                             ta_xxxx_xxxyy_1,  \
                             ta_xxxx_xxxyz_0,  \
                             ta_xxxx_xxxyz_1,  \
                             ta_xxxx_xxxz_0,   \
                             ta_xxxx_xxxz_1,   \
                             ta_xxxx_xxxzz_0,  \
                             ta_xxxx_xxxzz_1,  \
                             ta_xxxx_xxyy_0,   \
                             ta_xxxx_xxyy_1,   \
                             ta_xxxx_xxyyy_0,  \
                             ta_xxxx_xxyyy_1,  \
                             ta_xxxx_xxyyz_0,  \
                             ta_xxxx_xxyyz_1,  \
                             ta_xxxx_xxyz_0,   \
                             ta_xxxx_xxyz_1,   \
                             ta_xxxx_xxyzz_0,  \
                             ta_xxxx_xxyzz_1,  \
                             ta_xxxx_xxzz_0,   \
                             ta_xxxx_xxzz_1,   \
                             ta_xxxx_xxzzz_0,  \
                             ta_xxxx_xxzzz_1,  \
                             ta_xxxx_xyyy_0,   \
                             ta_xxxx_xyyy_1,   \
                             ta_xxxx_xyyyy_0,  \
                             ta_xxxx_xyyyy_1,  \
                             ta_xxxx_xyyyz_0,  \
                             ta_xxxx_xyyyz_1,  \
                             ta_xxxx_xyyz_0,   \
                             ta_xxxx_xyyz_1,   \
                             ta_xxxx_xyyzz_0,  \
                             ta_xxxx_xyyzz_1,  \
                             ta_xxxx_xyzz_0,   \
                             ta_xxxx_xyzz_1,   \
                             ta_xxxx_xyzzz_0,  \
                             ta_xxxx_xyzzz_1,  \
                             ta_xxxx_xzzz_0,   \
                             ta_xxxx_xzzz_1,   \
                             ta_xxxx_xzzzz_0,  \
                             ta_xxxx_xzzzz_1,  \
                             ta_xxxx_yyyy_0,   \
                             ta_xxxx_yyyy_1,   \
                             ta_xxxx_yyyyy_0,  \
                             ta_xxxx_yyyyy_1,  \
                             ta_xxxx_yyyyz_0,  \
                             ta_xxxx_yyyyz_1,  \
                             ta_xxxx_yyyz_0,   \
                             ta_xxxx_yyyz_1,   \
                             ta_xxxx_yyyzz_0,  \
                             ta_xxxx_yyyzz_1,  \
                             ta_xxxx_yyzz_0,   \
                             ta_xxxx_yyzz_1,   \
                             ta_xxxx_yyzzz_0,  \
                             ta_xxxx_yyzzz_1,  \
                             ta_xxxx_yzzz_0,   \
                             ta_xxxx_yzzz_1,   \
                             ta_xxxx_yzzzz_0,  \
                             ta_xxxx_yzzzz_1,  \
                             ta_xxxx_zzzz_0,   \
                             ta_xxxx_zzzz_1,   \
                             ta_xxxx_zzzzz_0,  \
                             ta_xxxx_zzzzz_1,  \
                             ta_xxxxx_xxxxx_0, \
                             ta_xxxxx_xxxxy_0, \
                             ta_xxxxx_xxxxz_0, \
                             ta_xxxxx_xxxyy_0, \
                             ta_xxxxx_xxxyz_0, \
                             ta_xxxxx_xxxzz_0, \
                             ta_xxxxx_xxyyy_0, \
                             ta_xxxxx_xxyyz_0, \
                             ta_xxxxx_xxyzz_0, \
                             ta_xxxxx_xxzzz_0, \
                             ta_xxxxx_xyyyy_0, \
                             ta_xxxxx_xyyyz_0, \
                             ta_xxxxx_xyyzz_0, \
                             ta_xxxxx_xyzzz_0, \
                             ta_xxxxx_xzzzz_0, \
                             ta_xxxxx_yyyyy_0, \
                             ta_xxxxx_yyyyz_0, \
                             ta_xxxxx_yyyzz_0, \
                             ta_xxxxx_yyzzz_0, \
                             ta_xxxxx_yzzzz_0, \
                             ta_xxxxx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_xxxxx_0[i] = 4.0 * ta_xxx_xxxxx_0[i] * fe_0 - 4.0 * ta_xxx_xxxxx_1[i] * fe_0 + 5.0 * ta_xxxx_xxxx_0[i] * fe_0 -
                              5.0 * ta_xxxx_xxxx_1[i] * fe_0 + ta_xxxx_xxxxx_0[i] * pa_x[i] - ta_xxxx_xxxxx_1[i] * pc_x[i];

        ta_xxxxx_xxxxy_0[i] = 4.0 * ta_xxx_xxxxy_0[i] * fe_0 - 4.0 * ta_xxx_xxxxy_1[i] * fe_0 + 4.0 * ta_xxxx_xxxy_0[i] * fe_0 -
                              4.0 * ta_xxxx_xxxy_1[i] * fe_0 + ta_xxxx_xxxxy_0[i] * pa_x[i] - ta_xxxx_xxxxy_1[i] * pc_x[i];

        ta_xxxxx_xxxxz_0[i] = 4.0 * ta_xxx_xxxxz_0[i] * fe_0 - 4.0 * ta_xxx_xxxxz_1[i] * fe_0 + 4.0 * ta_xxxx_xxxz_0[i] * fe_0 -
                              4.0 * ta_xxxx_xxxz_1[i] * fe_0 + ta_xxxx_xxxxz_0[i] * pa_x[i] - ta_xxxx_xxxxz_1[i] * pc_x[i];

        ta_xxxxx_xxxyy_0[i] = 4.0 * ta_xxx_xxxyy_0[i] * fe_0 - 4.0 * ta_xxx_xxxyy_1[i] * fe_0 + 3.0 * ta_xxxx_xxyy_0[i] * fe_0 -
                              3.0 * ta_xxxx_xxyy_1[i] * fe_0 + ta_xxxx_xxxyy_0[i] * pa_x[i] - ta_xxxx_xxxyy_1[i] * pc_x[i];

        ta_xxxxx_xxxyz_0[i] = 4.0 * ta_xxx_xxxyz_0[i] * fe_0 - 4.0 * ta_xxx_xxxyz_1[i] * fe_0 + 3.0 * ta_xxxx_xxyz_0[i] * fe_0 -
                              3.0 * ta_xxxx_xxyz_1[i] * fe_0 + ta_xxxx_xxxyz_0[i] * pa_x[i] - ta_xxxx_xxxyz_1[i] * pc_x[i];

        ta_xxxxx_xxxzz_0[i] = 4.0 * ta_xxx_xxxzz_0[i] * fe_0 - 4.0 * ta_xxx_xxxzz_1[i] * fe_0 + 3.0 * ta_xxxx_xxzz_0[i] * fe_0 -
                              3.0 * ta_xxxx_xxzz_1[i] * fe_0 + ta_xxxx_xxxzz_0[i] * pa_x[i] - ta_xxxx_xxxzz_1[i] * pc_x[i];

        ta_xxxxx_xxyyy_0[i] = 4.0 * ta_xxx_xxyyy_0[i] * fe_0 - 4.0 * ta_xxx_xxyyy_1[i] * fe_0 + 2.0 * ta_xxxx_xyyy_0[i] * fe_0 -
                              2.0 * ta_xxxx_xyyy_1[i] * fe_0 + ta_xxxx_xxyyy_0[i] * pa_x[i] - ta_xxxx_xxyyy_1[i] * pc_x[i];

        ta_xxxxx_xxyyz_0[i] = 4.0 * ta_xxx_xxyyz_0[i] * fe_0 - 4.0 * ta_xxx_xxyyz_1[i] * fe_0 + 2.0 * ta_xxxx_xyyz_0[i] * fe_0 -
                              2.0 * ta_xxxx_xyyz_1[i] * fe_0 + ta_xxxx_xxyyz_0[i] * pa_x[i] - ta_xxxx_xxyyz_1[i] * pc_x[i];

        ta_xxxxx_xxyzz_0[i] = 4.0 * ta_xxx_xxyzz_0[i] * fe_0 - 4.0 * ta_xxx_xxyzz_1[i] * fe_0 + 2.0 * ta_xxxx_xyzz_0[i] * fe_0 -
                              2.0 * ta_xxxx_xyzz_1[i] * fe_0 + ta_xxxx_xxyzz_0[i] * pa_x[i] - ta_xxxx_xxyzz_1[i] * pc_x[i];

        ta_xxxxx_xxzzz_0[i] = 4.0 * ta_xxx_xxzzz_0[i] * fe_0 - 4.0 * ta_xxx_xxzzz_1[i] * fe_0 + 2.0 * ta_xxxx_xzzz_0[i] * fe_0 -
                              2.0 * ta_xxxx_xzzz_1[i] * fe_0 + ta_xxxx_xxzzz_0[i] * pa_x[i] - ta_xxxx_xxzzz_1[i] * pc_x[i];

        ta_xxxxx_xyyyy_0[i] = 4.0 * ta_xxx_xyyyy_0[i] * fe_0 - 4.0 * ta_xxx_xyyyy_1[i] * fe_0 + ta_xxxx_yyyy_0[i] * fe_0 - ta_xxxx_yyyy_1[i] * fe_0 +
                              ta_xxxx_xyyyy_0[i] * pa_x[i] - ta_xxxx_xyyyy_1[i] * pc_x[i];

        ta_xxxxx_xyyyz_0[i] = 4.0 * ta_xxx_xyyyz_0[i] * fe_0 - 4.0 * ta_xxx_xyyyz_1[i] * fe_0 + ta_xxxx_yyyz_0[i] * fe_0 - ta_xxxx_yyyz_1[i] * fe_0 +
                              ta_xxxx_xyyyz_0[i] * pa_x[i] - ta_xxxx_xyyyz_1[i] * pc_x[i];

        ta_xxxxx_xyyzz_0[i] = 4.0 * ta_xxx_xyyzz_0[i] * fe_0 - 4.0 * ta_xxx_xyyzz_1[i] * fe_0 + ta_xxxx_yyzz_0[i] * fe_0 - ta_xxxx_yyzz_1[i] * fe_0 +
                              ta_xxxx_xyyzz_0[i] * pa_x[i] - ta_xxxx_xyyzz_1[i] * pc_x[i];

        ta_xxxxx_xyzzz_0[i] = 4.0 * ta_xxx_xyzzz_0[i] * fe_0 - 4.0 * ta_xxx_xyzzz_1[i] * fe_0 + ta_xxxx_yzzz_0[i] * fe_0 - ta_xxxx_yzzz_1[i] * fe_0 +
                              ta_xxxx_xyzzz_0[i] * pa_x[i] - ta_xxxx_xyzzz_1[i] * pc_x[i];

        ta_xxxxx_xzzzz_0[i] = 4.0 * ta_xxx_xzzzz_0[i] * fe_0 - 4.0 * ta_xxx_xzzzz_1[i] * fe_0 + ta_xxxx_zzzz_0[i] * fe_0 - ta_xxxx_zzzz_1[i] * fe_0 +
                              ta_xxxx_xzzzz_0[i] * pa_x[i] - ta_xxxx_xzzzz_1[i] * pc_x[i];

        ta_xxxxx_yyyyy_0[i] =
            4.0 * ta_xxx_yyyyy_0[i] * fe_0 - 4.0 * ta_xxx_yyyyy_1[i] * fe_0 + ta_xxxx_yyyyy_0[i] * pa_x[i] - ta_xxxx_yyyyy_1[i] * pc_x[i];

        ta_xxxxx_yyyyz_0[i] =
            4.0 * ta_xxx_yyyyz_0[i] * fe_0 - 4.0 * ta_xxx_yyyyz_1[i] * fe_0 + ta_xxxx_yyyyz_0[i] * pa_x[i] - ta_xxxx_yyyyz_1[i] * pc_x[i];

        ta_xxxxx_yyyzz_0[i] =
            4.0 * ta_xxx_yyyzz_0[i] * fe_0 - 4.0 * ta_xxx_yyyzz_1[i] * fe_0 + ta_xxxx_yyyzz_0[i] * pa_x[i] - ta_xxxx_yyyzz_1[i] * pc_x[i];

        ta_xxxxx_yyzzz_0[i] =
            4.0 * ta_xxx_yyzzz_0[i] * fe_0 - 4.0 * ta_xxx_yyzzz_1[i] * fe_0 + ta_xxxx_yyzzz_0[i] * pa_x[i] - ta_xxxx_yyzzz_1[i] * pc_x[i];

        ta_xxxxx_yzzzz_0[i] =
            4.0 * ta_xxx_yzzzz_0[i] * fe_0 - 4.0 * ta_xxx_yzzzz_1[i] * fe_0 + ta_xxxx_yzzzz_0[i] * pa_x[i] - ta_xxxx_yzzzz_1[i] * pc_x[i];

        ta_xxxxx_zzzzz_0[i] =
            4.0 * ta_xxx_zzzzz_0[i] * fe_0 - 4.0 * ta_xxx_zzzzz_1[i] * fe_0 + ta_xxxx_zzzzz_0[i] * pa_x[i] - ta_xxxx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : HH

    auto ta_xxxxy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 21);

    auto ta_xxxxy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 22);

    auto ta_xxxxy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 23);

    auto ta_xxxxy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 24);

    auto ta_xxxxy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 25);

    auto ta_xxxxy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 26);

    auto ta_xxxxy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 27);

    auto ta_xxxxy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 28);

    auto ta_xxxxy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 29);

    auto ta_xxxxy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 30);

    auto ta_xxxxy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 31);

    auto ta_xxxxy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 32);

    auto ta_xxxxy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 33);

    auto ta_xxxxy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 34);

    auto ta_xxxxy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 35);

    auto ta_xxxxy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 36);

    auto ta_xxxxy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 37);

    auto ta_xxxxy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 38);

    auto ta_xxxxy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 39);

    auto ta_xxxxy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 40);

    auto ta_xxxxy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 41);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxxx_xxxx_0,   \
                             ta_xxxx_xxxx_1,   \
                             ta_xxxx_xxxxx_0,  \
                             ta_xxxx_xxxxx_1,  \
                             ta_xxxx_xxxxy_0,  \
                             ta_xxxx_xxxxy_1,  \
                             ta_xxxx_xxxxz_0,  \
                             ta_xxxx_xxxxz_1,  \
                             ta_xxxx_xxxy_0,   \
                             ta_xxxx_xxxy_1,   \
                             ta_xxxx_xxxyy_0,  \
                             ta_xxxx_xxxyy_1,  \
                             ta_xxxx_xxxyz_0,  \
                             ta_xxxx_xxxyz_1,  \
                             ta_xxxx_xxxz_0,   \
                             ta_xxxx_xxxz_1,   \
                             ta_xxxx_xxxzz_0,  \
                             ta_xxxx_xxxzz_1,  \
                             ta_xxxx_xxyy_0,   \
                             ta_xxxx_xxyy_1,   \
                             ta_xxxx_xxyyy_0,  \
                             ta_xxxx_xxyyy_1,  \
                             ta_xxxx_xxyyz_0,  \
                             ta_xxxx_xxyyz_1,  \
                             ta_xxxx_xxyz_0,   \
                             ta_xxxx_xxyz_1,   \
                             ta_xxxx_xxyzz_0,  \
                             ta_xxxx_xxyzz_1,  \
                             ta_xxxx_xxzz_0,   \
                             ta_xxxx_xxzz_1,   \
                             ta_xxxx_xxzzz_0,  \
                             ta_xxxx_xxzzz_1,  \
                             ta_xxxx_xyyy_0,   \
                             ta_xxxx_xyyy_1,   \
                             ta_xxxx_xyyyy_0,  \
                             ta_xxxx_xyyyy_1,  \
                             ta_xxxx_xyyyz_0,  \
                             ta_xxxx_xyyyz_1,  \
                             ta_xxxx_xyyz_0,   \
                             ta_xxxx_xyyz_1,   \
                             ta_xxxx_xyyzz_0,  \
                             ta_xxxx_xyyzz_1,  \
                             ta_xxxx_xyzz_0,   \
                             ta_xxxx_xyzz_1,   \
                             ta_xxxx_xyzzz_0,  \
                             ta_xxxx_xyzzz_1,  \
                             ta_xxxx_xzzz_0,   \
                             ta_xxxx_xzzz_1,   \
                             ta_xxxx_xzzzz_0,  \
                             ta_xxxx_xzzzz_1,  \
                             ta_xxxx_zzzzz_0,  \
                             ta_xxxx_zzzzz_1,  \
                             ta_xxxxy_xxxxx_0, \
                             ta_xxxxy_xxxxy_0, \
                             ta_xxxxy_xxxxz_0, \
                             ta_xxxxy_xxxyy_0, \
                             ta_xxxxy_xxxyz_0, \
                             ta_xxxxy_xxxzz_0, \
                             ta_xxxxy_xxyyy_0, \
                             ta_xxxxy_xxyyz_0, \
                             ta_xxxxy_xxyzz_0, \
                             ta_xxxxy_xxzzz_0, \
                             ta_xxxxy_xyyyy_0, \
                             ta_xxxxy_xyyyz_0, \
                             ta_xxxxy_xyyzz_0, \
                             ta_xxxxy_xyzzz_0, \
                             ta_xxxxy_xzzzz_0, \
                             ta_xxxxy_yyyyy_0, \
                             ta_xxxxy_yyyyz_0, \
                             ta_xxxxy_yyyzz_0, \
                             ta_xxxxy_yyzzz_0, \
                             ta_xxxxy_yzzzz_0, \
                             ta_xxxxy_zzzzz_0, \
                             ta_xxxy_yyyyy_0,  \
                             ta_xxxy_yyyyy_1,  \
                             ta_xxxy_yyyyz_0,  \
                             ta_xxxy_yyyyz_1,  \
                             ta_xxxy_yyyzz_0,  \
                             ta_xxxy_yyyzz_1,  \
                             ta_xxxy_yyzzz_0,  \
                             ta_xxxy_yyzzz_1,  \
                             ta_xxxy_yzzzz_0,  \
                             ta_xxxy_yzzzz_1,  \
                             ta_xxy_yyyyy_0,   \
                             ta_xxy_yyyyy_1,   \
                             ta_xxy_yyyyz_0,   \
                             ta_xxy_yyyyz_1,   \
                             ta_xxy_yyyzz_0,   \
                             ta_xxy_yyyzz_1,   \
                             ta_xxy_yyzzz_0,   \
                             ta_xxy_yyzzz_1,   \
                             ta_xxy_yzzzz_0,   \
                             ta_xxy_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxy_xxxxx_0[i] = ta_xxxx_xxxxx_0[i] * pa_y[i] - ta_xxxx_xxxxx_1[i] * pc_y[i];

        ta_xxxxy_xxxxy_0[i] = ta_xxxx_xxxx_0[i] * fe_0 - ta_xxxx_xxxx_1[i] * fe_0 + ta_xxxx_xxxxy_0[i] * pa_y[i] - ta_xxxx_xxxxy_1[i] * pc_y[i];

        ta_xxxxy_xxxxz_0[i] = ta_xxxx_xxxxz_0[i] * pa_y[i] - ta_xxxx_xxxxz_1[i] * pc_y[i];

        ta_xxxxy_xxxyy_0[i] =
            2.0 * ta_xxxx_xxxy_0[i] * fe_0 - 2.0 * ta_xxxx_xxxy_1[i] * fe_0 + ta_xxxx_xxxyy_0[i] * pa_y[i] - ta_xxxx_xxxyy_1[i] * pc_y[i];

        ta_xxxxy_xxxyz_0[i] = ta_xxxx_xxxz_0[i] * fe_0 - ta_xxxx_xxxz_1[i] * fe_0 + ta_xxxx_xxxyz_0[i] * pa_y[i] - ta_xxxx_xxxyz_1[i] * pc_y[i];

        ta_xxxxy_xxxzz_0[i] = ta_xxxx_xxxzz_0[i] * pa_y[i] - ta_xxxx_xxxzz_1[i] * pc_y[i];

        ta_xxxxy_xxyyy_0[i] =
            3.0 * ta_xxxx_xxyy_0[i] * fe_0 - 3.0 * ta_xxxx_xxyy_1[i] * fe_0 + ta_xxxx_xxyyy_0[i] * pa_y[i] - ta_xxxx_xxyyy_1[i] * pc_y[i];

        ta_xxxxy_xxyyz_0[i] =
            2.0 * ta_xxxx_xxyz_0[i] * fe_0 - 2.0 * ta_xxxx_xxyz_1[i] * fe_0 + ta_xxxx_xxyyz_0[i] * pa_y[i] - ta_xxxx_xxyyz_1[i] * pc_y[i];

        ta_xxxxy_xxyzz_0[i] = ta_xxxx_xxzz_0[i] * fe_0 - ta_xxxx_xxzz_1[i] * fe_0 + ta_xxxx_xxyzz_0[i] * pa_y[i] - ta_xxxx_xxyzz_1[i] * pc_y[i];

        ta_xxxxy_xxzzz_0[i] = ta_xxxx_xxzzz_0[i] * pa_y[i] - ta_xxxx_xxzzz_1[i] * pc_y[i];

        ta_xxxxy_xyyyy_0[i] =
            4.0 * ta_xxxx_xyyy_0[i] * fe_0 - 4.0 * ta_xxxx_xyyy_1[i] * fe_0 + ta_xxxx_xyyyy_0[i] * pa_y[i] - ta_xxxx_xyyyy_1[i] * pc_y[i];

        ta_xxxxy_xyyyz_0[i] =
            3.0 * ta_xxxx_xyyz_0[i] * fe_0 - 3.0 * ta_xxxx_xyyz_1[i] * fe_0 + ta_xxxx_xyyyz_0[i] * pa_y[i] - ta_xxxx_xyyyz_1[i] * pc_y[i];

        ta_xxxxy_xyyzz_0[i] =
            2.0 * ta_xxxx_xyzz_0[i] * fe_0 - 2.0 * ta_xxxx_xyzz_1[i] * fe_0 + ta_xxxx_xyyzz_0[i] * pa_y[i] - ta_xxxx_xyyzz_1[i] * pc_y[i];

        ta_xxxxy_xyzzz_0[i] = ta_xxxx_xzzz_0[i] * fe_0 - ta_xxxx_xzzz_1[i] * fe_0 + ta_xxxx_xyzzz_0[i] * pa_y[i] - ta_xxxx_xyzzz_1[i] * pc_y[i];

        ta_xxxxy_xzzzz_0[i] = ta_xxxx_xzzzz_0[i] * pa_y[i] - ta_xxxx_xzzzz_1[i] * pc_y[i];

        ta_xxxxy_yyyyy_0[i] =
            3.0 * ta_xxy_yyyyy_0[i] * fe_0 - 3.0 * ta_xxy_yyyyy_1[i] * fe_0 + ta_xxxy_yyyyy_0[i] * pa_x[i] - ta_xxxy_yyyyy_1[i] * pc_x[i];

        ta_xxxxy_yyyyz_0[i] =
            3.0 * ta_xxy_yyyyz_0[i] * fe_0 - 3.0 * ta_xxy_yyyyz_1[i] * fe_0 + ta_xxxy_yyyyz_0[i] * pa_x[i] - ta_xxxy_yyyyz_1[i] * pc_x[i];

        ta_xxxxy_yyyzz_0[i] =
            3.0 * ta_xxy_yyyzz_0[i] * fe_0 - 3.0 * ta_xxy_yyyzz_1[i] * fe_0 + ta_xxxy_yyyzz_0[i] * pa_x[i] - ta_xxxy_yyyzz_1[i] * pc_x[i];

        ta_xxxxy_yyzzz_0[i] =
            3.0 * ta_xxy_yyzzz_0[i] * fe_0 - 3.0 * ta_xxy_yyzzz_1[i] * fe_0 + ta_xxxy_yyzzz_0[i] * pa_x[i] - ta_xxxy_yyzzz_1[i] * pc_x[i];

        ta_xxxxy_yzzzz_0[i] =
            3.0 * ta_xxy_yzzzz_0[i] * fe_0 - 3.0 * ta_xxy_yzzzz_1[i] * fe_0 + ta_xxxy_yzzzz_0[i] * pa_x[i] - ta_xxxy_yzzzz_1[i] * pc_x[i];

        ta_xxxxy_zzzzz_0[i] = ta_xxxx_zzzzz_0[i] * pa_y[i] - ta_xxxx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : HH

    auto ta_xxxxz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 42);

    auto ta_xxxxz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 43);

    auto ta_xxxxz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 44);

    auto ta_xxxxz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 45);

    auto ta_xxxxz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 46);

    auto ta_xxxxz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 47);

    auto ta_xxxxz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 48);

    auto ta_xxxxz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 49);

    auto ta_xxxxz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 50);

    auto ta_xxxxz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 51);

    auto ta_xxxxz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 52);

    auto ta_xxxxz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 53);

    auto ta_xxxxz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 54);

    auto ta_xxxxz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 55);

    auto ta_xxxxz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 56);

    auto ta_xxxxz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 57);

    auto ta_xxxxz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 58);

    auto ta_xxxxz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 59);

    auto ta_xxxxz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 60);

    auto ta_xxxxz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 61);

    auto ta_xxxxz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 62);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxxx_xxxx_0,   \
                             ta_xxxx_xxxx_1,   \
                             ta_xxxx_xxxxx_0,  \
                             ta_xxxx_xxxxx_1,  \
                             ta_xxxx_xxxxy_0,  \
                             ta_xxxx_xxxxy_1,  \
                             ta_xxxx_xxxxz_0,  \
                             ta_xxxx_xxxxz_1,  \
                             ta_xxxx_xxxy_0,   \
                             ta_xxxx_xxxy_1,   \
                             ta_xxxx_xxxyy_0,  \
                             ta_xxxx_xxxyy_1,  \
                             ta_xxxx_xxxyz_0,  \
                             ta_xxxx_xxxyz_1,  \
                             ta_xxxx_xxxz_0,   \
                             ta_xxxx_xxxz_1,   \
                             ta_xxxx_xxxzz_0,  \
                             ta_xxxx_xxxzz_1,  \
                             ta_xxxx_xxyy_0,   \
                             ta_xxxx_xxyy_1,   \
                             ta_xxxx_xxyyy_0,  \
                             ta_xxxx_xxyyy_1,  \
                             ta_xxxx_xxyyz_0,  \
                             ta_xxxx_xxyyz_1,  \
                             ta_xxxx_xxyz_0,   \
                             ta_xxxx_xxyz_1,   \
                             ta_xxxx_xxyzz_0,  \
                             ta_xxxx_xxyzz_1,  \
                             ta_xxxx_xxzz_0,   \
                             ta_xxxx_xxzz_1,   \
                             ta_xxxx_xxzzz_0,  \
                             ta_xxxx_xxzzz_1,  \
                             ta_xxxx_xyyy_0,   \
                             ta_xxxx_xyyy_1,   \
                             ta_xxxx_xyyyy_0,  \
                             ta_xxxx_xyyyy_1,  \
                             ta_xxxx_xyyyz_0,  \
                             ta_xxxx_xyyyz_1,  \
                             ta_xxxx_xyyz_0,   \
                             ta_xxxx_xyyz_1,   \
                             ta_xxxx_xyyzz_0,  \
                             ta_xxxx_xyyzz_1,  \
                             ta_xxxx_xyzz_0,   \
                             ta_xxxx_xyzz_1,   \
                             ta_xxxx_xyzzz_0,  \
                             ta_xxxx_xyzzz_1,  \
                             ta_xxxx_xzzz_0,   \
                             ta_xxxx_xzzz_1,   \
                             ta_xxxx_xzzzz_0,  \
                             ta_xxxx_xzzzz_1,  \
                             ta_xxxx_yyyyy_0,  \
                             ta_xxxx_yyyyy_1,  \
                             ta_xxxxz_xxxxx_0, \
                             ta_xxxxz_xxxxy_0, \
                             ta_xxxxz_xxxxz_0, \
                             ta_xxxxz_xxxyy_0, \
                             ta_xxxxz_xxxyz_0, \
                             ta_xxxxz_xxxzz_0, \
                             ta_xxxxz_xxyyy_0, \
                             ta_xxxxz_xxyyz_0, \
                             ta_xxxxz_xxyzz_0, \
                             ta_xxxxz_xxzzz_0, \
                             ta_xxxxz_xyyyy_0, \
                             ta_xxxxz_xyyyz_0, \
                             ta_xxxxz_xyyzz_0, \
                             ta_xxxxz_xyzzz_0, \
                             ta_xxxxz_xzzzz_0, \
                             ta_xxxxz_yyyyy_0, \
                             ta_xxxxz_yyyyz_0, \
                             ta_xxxxz_yyyzz_0, \
                             ta_xxxxz_yyzzz_0, \
                             ta_xxxxz_yzzzz_0, \
                             ta_xxxxz_zzzzz_0, \
                             ta_xxxz_yyyyz_0,  \
                             ta_xxxz_yyyyz_1,  \
                             ta_xxxz_yyyzz_0,  \
                             ta_xxxz_yyyzz_1,  \
                             ta_xxxz_yyzzz_0,  \
                             ta_xxxz_yyzzz_1,  \
                             ta_xxxz_yzzzz_0,  \
                             ta_xxxz_yzzzz_1,  \
                             ta_xxxz_zzzzz_0,  \
                             ta_xxxz_zzzzz_1,  \
                             ta_xxz_yyyyz_0,   \
                             ta_xxz_yyyyz_1,   \
                             ta_xxz_yyyzz_0,   \
                             ta_xxz_yyyzz_1,   \
                             ta_xxz_yyzzz_0,   \
                             ta_xxz_yyzzz_1,   \
                             ta_xxz_yzzzz_0,   \
                             ta_xxz_yzzzz_1,   \
                             ta_xxz_zzzzz_0,   \
                             ta_xxz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxz_xxxxx_0[i] = ta_xxxx_xxxxx_0[i] * pa_z[i] - ta_xxxx_xxxxx_1[i] * pc_z[i];

        ta_xxxxz_xxxxy_0[i] = ta_xxxx_xxxxy_0[i] * pa_z[i] - ta_xxxx_xxxxy_1[i] * pc_z[i];

        ta_xxxxz_xxxxz_0[i] = ta_xxxx_xxxx_0[i] * fe_0 - ta_xxxx_xxxx_1[i] * fe_0 + ta_xxxx_xxxxz_0[i] * pa_z[i] - ta_xxxx_xxxxz_1[i] * pc_z[i];

        ta_xxxxz_xxxyy_0[i] = ta_xxxx_xxxyy_0[i] * pa_z[i] - ta_xxxx_xxxyy_1[i] * pc_z[i];

        ta_xxxxz_xxxyz_0[i] = ta_xxxx_xxxy_0[i] * fe_0 - ta_xxxx_xxxy_1[i] * fe_0 + ta_xxxx_xxxyz_0[i] * pa_z[i] - ta_xxxx_xxxyz_1[i] * pc_z[i];

        ta_xxxxz_xxxzz_0[i] =
            2.0 * ta_xxxx_xxxz_0[i] * fe_0 - 2.0 * ta_xxxx_xxxz_1[i] * fe_0 + ta_xxxx_xxxzz_0[i] * pa_z[i] - ta_xxxx_xxxzz_1[i] * pc_z[i];

        ta_xxxxz_xxyyy_0[i] = ta_xxxx_xxyyy_0[i] * pa_z[i] - ta_xxxx_xxyyy_1[i] * pc_z[i];

        ta_xxxxz_xxyyz_0[i] = ta_xxxx_xxyy_0[i] * fe_0 - ta_xxxx_xxyy_1[i] * fe_0 + ta_xxxx_xxyyz_0[i] * pa_z[i] - ta_xxxx_xxyyz_1[i] * pc_z[i];

        ta_xxxxz_xxyzz_0[i] =
            2.0 * ta_xxxx_xxyz_0[i] * fe_0 - 2.0 * ta_xxxx_xxyz_1[i] * fe_0 + ta_xxxx_xxyzz_0[i] * pa_z[i] - ta_xxxx_xxyzz_1[i] * pc_z[i];

        ta_xxxxz_xxzzz_0[i] =
            3.0 * ta_xxxx_xxzz_0[i] * fe_0 - 3.0 * ta_xxxx_xxzz_1[i] * fe_0 + ta_xxxx_xxzzz_0[i] * pa_z[i] - ta_xxxx_xxzzz_1[i] * pc_z[i];

        ta_xxxxz_xyyyy_0[i] = ta_xxxx_xyyyy_0[i] * pa_z[i] - ta_xxxx_xyyyy_1[i] * pc_z[i];

        ta_xxxxz_xyyyz_0[i] = ta_xxxx_xyyy_0[i] * fe_0 - ta_xxxx_xyyy_1[i] * fe_0 + ta_xxxx_xyyyz_0[i] * pa_z[i] - ta_xxxx_xyyyz_1[i] * pc_z[i];

        ta_xxxxz_xyyzz_0[i] =
            2.0 * ta_xxxx_xyyz_0[i] * fe_0 - 2.0 * ta_xxxx_xyyz_1[i] * fe_0 + ta_xxxx_xyyzz_0[i] * pa_z[i] - ta_xxxx_xyyzz_1[i] * pc_z[i];

        ta_xxxxz_xyzzz_0[i] =
            3.0 * ta_xxxx_xyzz_0[i] * fe_0 - 3.0 * ta_xxxx_xyzz_1[i] * fe_0 + ta_xxxx_xyzzz_0[i] * pa_z[i] - ta_xxxx_xyzzz_1[i] * pc_z[i];

        ta_xxxxz_xzzzz_0[i] =
            4.0 * ta_xxxx_xzzz_0[i] * fe_0 - 4.0 * ta_xxxx_xzzz_1[i] * fe_0 + ta_xxxx_xzzzz_0[i] * pa_z[i] - ta_xxxx_xzzzz_1[i] * pc_z[i];

        ta_xxxxz_yyyyy_0[i] = ta_xxxx_yyyyy_0[i] * pa_z[i] - ta_xxxx_yyyyy_1[i] * pc_z[i];

        ta_xxxxz_yyyyz_0[i] =
            3.0 * ta_xxz_yyyyz_0[i] * fe_0 - 3.0 * ta_xxz_yyyyz_1[i] * fe_0 + ta_xxxz_yyyyz_0[i] * pa_x[i] - ta_xxxz_yyyyz_1[i] * pc_x[i];

        ta_xxxxz_yyyzz_0[i] =
            3.0 * ta_xxz_yyyzz_0[i] * fe_0 - 3.0 * ta_xxz_yyyzz_1[i] * fe_0 + ta_xxxz_yyyzz_0[i] * pa_x[i] - ta_xxxz_yyyzz_1[i] * pc_x[i];

        ta_xxxxz_yyzzz_0[i] =
            3.0 * ta_xxz_yyzzz_0[i] * fe_0 - 3.0 * ta_xxz_yyzzz_1[i] * fe_0 + ta_xxxz_yyzzz_0[i] * pa_x[i] - ta_xxxz_yyzzz_1[i] * pc_x[i];

        ta_xxxxz_yzzzz_0[i] =
            3.0 * ta_xxz_yzzzz_0[i] * fe_0 - 3.0 * ta_xxz_yzzzz_1[i] * fe_0 + ta_xxxz_yzzzz_0[i] * pa_x[i] - ta_xxxz_yzzzz_1[i] * pc_x[i];

        ta_xxxxz_zzzzz_0[i] =
            3.0 * ta_xxz_zzzzz_0[i] * fe_0 - 3.0 * ta_xxz_zzzzz_1[i] * fe_0 + ta_xxxz_zzzzz_0[i] * pa_x[i] - ta_xxxz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 63-84 components of targeted buffer : HH

    auto ta_xxxyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 63);

    auto ta_xxxyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 64);

    auto ta_xxxyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 65);

    auto ta_xxxyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 66);

    auto ta_xxxyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 67);

    auto ta_xxxyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 68);

    auto ta_xxxyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 69);

    auto ta_xxxyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 70);

    auto ta_xxxyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 71);

    auto ta_xxxyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 72);

    auto ta_xxxyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 73);

    auto ta_xxxyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 74);

    auto ta_xxxyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 75);

    auto ta_xxxyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 76);

    auto ta_xxxyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 77);

    auto ta_xxxyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 78);

    auto ta_xxxyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 79);

    auto ta_xxxyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 80);

    auto ta_xxxyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 81);

    auto ta_xxxyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 82);

    auto ta_xxxyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 83);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxx_xxxxx_0,   \
                             ta_xxx_xxxxx_1,   \
                             ta_xxx_xxxxz_0,   \
                             ta_xxx_xxxxz_1,   \
                             ta_xxx_xxxzz_0,   \
                             ta_xxx_xxxzz_1,   \
                             ta_xxx_xxzzz_0,   \
                             ta_xxx_xxzzz_1,   \
                             ta_xxx_xzzzz_0,   \
                             ta_xxx_xzzzz_1,   \
                             ta_xxxy_xxxxx_0,  \
                             ta_xxxy_xxxxx_1,  \
                             ta_xxxy_xxxxz_0,  \
                             ta_xxxy_xxxxz_1,  \
                             ta_xxxy_xxxzz_0,  \
                             ta_xxxy_xxxzz_1,  \
                             ta_xxxy_xxzzz_0,  \
                             ta_xxxy_xxzzz_1,  \
                             ta_xxxy_xzzzz_0,  \
                             ta_xxxy_xzzzz_1,  \
                             ta_xxxyy_xxxxx_0, \
                             ta_xxxyy_xxxxy_0, \
                             ta_xxxyy_xxxxz_0, \
                             ta_xxxyy_xxxyy_0, \
                             ta_xxxyy_xxxyz_0, \
                             ta_xxxyy_xxxzz_0, \
                             ta_xxxyy_xxyyy_0, \
                             ta_xxxyy_xxyyz_0, \
                             ta_xxxyy_xxyzz_0, \
                             ta_xxxyy_xxzzz_0, \
                             ta_xxxyy_xyyyy_0, \
                             ta_xxxyy_xyyyz_0, \
                             ta_xxxyy_xyyzz_0, \
                             ta_xxxyy_xyzzz_0, \
                             ta_xxxyy_xzzzz_0, \
                             ta_xxxyy_yyyyy_0, \
                             ta_xxxyy_yyyyz_0, \
                             ta_xxxyy_yyyzz_0, \
                             ta_xxxyy_yyzzz_0, \
                             ta_xxxyy_yzzzz_0, \
                             ta_xxxyy_zzzzz_0, \
                             ta_xxyy_xxxxy_0,  \
                             ta_xxyy_xxxxy_1,  \
                             ta_xxyy_xxxy_0,   \
                             ta_xxyy_xxxy_1,   \
                             ta_xxyy_xxxyy_0,  \
                             ta_xxyy_xxxyy_1,  \
                             ta_xxyy_xxxyz_0,  \
                             ta_xxyy_xxxyz_1,  \
                             ta_xxyy_xxyy_0,   \
                             ta_xxyy_xxyy_1,   \
                             ta_xxyy_xxyyy_0,  \
                             ta_xxyy_xxyyy_1,  \
                             ta_xxyy_xxyyz_0,  \
                             ta_xxyy_xxyyz_1,  \
                             ta_xxyy_xxyz_0,   \
                             ta_xxyy_xxyz_1,   \
                             ta_xxyy_xxyzz_0,  \
                             ta_xxyy_xxyzz_1,  \
                             ta_xxyy_xyyy_0,   \
                             ta_xxyy_xyyy_1,   \
                             ta_xxyy_xyyyy_0,  \
                             ta_xxyy_xyyyy_1,  \
                             ta_xxyy_xyyyz_0,  \
                             ta_xxyy_xyyyz_1,  \
                             ta_xxyy_xyyz_0,   \
                             ta_xxyy_xyyz_1,   \
                             ta_xxyy_xyyzz_0,  \
                             ta_xxyy_xyyzz_1,  \
                             ta_xxyy_xyzz_0,   \
                             ta_xxyy_xyzz_1,   \
                             ta_xxyy_xyzzz_0,  \
                             ta_xxyy_xyzzz_1,  \
                             ta_xxyy_yyyy_0,   \
                             ta_xxyy_yyyy_1,   \
                             ta_xxyy_yyyyy_0,  \
                             ta_xxyy_yyyyy_1,  \
                             ta_xxyy_yyyyz_0,  \
                             ta_xxyy_yyyyz_1,  \
                             ta_xxyy_yyyz_0,   \
                             ta_xxyy_yyyz_1,   \
                             ta_xxyy_yyyzz_0,  \
                             ta_xxyy_yyyzz_1,  \
                             ta_xxyy_yyzz_0,   \
                             ta_xxyy_yyzz_1,   \
                             ta_xxyy_yyzzz_0,  \
                             ta_xxyy_yyzzz_1,  \
                             ta_xxyy_yzzz_0,   \
                             ta_xxyy_yzzz_1,   \
                             ta_xxyy_yzzzz_0,  \
                             ta_xxyy_yzzzz_1,  \
                             ta_xxyy_zzzzz_0,  \
                             ta_xxyy_zzzzz_1,  \
                             ta_xyy_xxxxy_0,   \
                             ta_xyy_xxxxy_1,   \
                             ta_xyy_xxxyy_0,   \
                             ta_xyy_xxxyy_1,   \
                             ta_xyy_xxxyz_0,   \
                             ta_xyy_xxxyz_1,   \
                             ta_xyy_xxyyy_0,   \
                             ta_xyy_xxyyy_1,   \
                             ta_xyy_xxyyz_0,   \
                             ta_xyy_xxyyz_1,   \
                             ta_xyy_xxyzz_0,   \
                             ta_xyy_xxyzz_1,   \
                             ta_xyy_xyyyy_0,   \
                             ta_xyy_xyyyy_1,   \
                             ta_xyy_xyyyz_0,   \
                             ta_xyy_xyyyz_1,   \
                             ta_xyy_xyyzz_0,   \
                             ta_xyy_xyyzz_1,   \
                             ta_xyy_xyzzz_0,   \
                             ta_xyy_xyzzz_1,   \
                             ta_xyy_yyyyy_0,   \
                             ta_xyy_yyyyy_1,   \
                             ta_xyy_yyyyz_0,   \
                             ta_xyy_yyyyz_1,   \
                             ta_xyy_yyyzz_0,   \
                             ta_xyy_yyyzz_1,   \
                             ta_xyy_yyzzz_0,   \
                             ta_xyy_yyzzz_1,   \
                             ta_xyy_yzzzz_0,   \
                             ta_xyy_yzzzz_1,   \
                             ta_xyy_zzzzz_0,   \
                             ta_xyy_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyy_xxxxx_0[i] = ta_xxx_xxxxx_0[i] * fe_0 - ta_xxx_xxxxx_1[i] * fe_0 + ta_xxxy_xxxxx_0[i] * pa_y[i] - ta_xxxy_xxxxx_1[i] * pc_y[i];

        ta_xxxyy_xxxxy_0[i] = 2.0 * ta_xyy_xxxxy_0[i] * fe_0 - 2.0 * ta_xyy_xxxxy_1[i] * fe_0 + 4.0 * ta_xxyy_xxxy_0[i] * fe_0 -
                              4.0 * ta_xxyy_xxxy_1[i] * fe_0 + ta_xxyy_xxxxy_0[i] * pa_x[i] - ta_xxyy_xxxxy_1[i] * pc_x[i];

        ta_xxxyy_xxxxz_0[i] = ta_xxx_xxxxz_0[i] * fe_0 - ta_xxx_xxxxz_1[i] * fe_0 + ta_xxxy_xxxxz_0[i] * pa_y[i] - ta_xxxy_xxxxz_1[i] * pc_y[i];

        ta_xxxyy_xxxyy_0[i] = 2.0 * ta_xyy_xxxyy_0[i] * fe_0 - 2.0 * ta_xyy_xxxyy_1[i] * fe_0 + 3.0 * ta_xxyy_xxyy_0[i] * fe_0 -
                              3.0 * ta_xxyy_xxyy_1[i] * fe_0 + ta_xxyy_xxxyy_0[i] * pa_x[i] - ta_xxyy_xxxyy_1[i] * pc_x[i];

        ta_xxxyy_xxxyz_0[i] = 2.0 * ta_xyy_xxxyz_0[i] * fe_0 - 2.0 * ta_xyy_xxxyz_1[i] * fe_0 + 3.0 * ta_xxyy_xxyz_0[i] * fe_0 -
                              3.0 * ta_xxyy_xxyz_1[i] * fe_0 + ta_xxyy_xxxyz_0[i] * pa_x[i] - ta_xxyy_xxxyz_1[i] * pc_x[i];

        ta_xxxyy_xxxzz_0[i] = ta_xxx_xxxzz_0[i] * fe_0 - ta_xxx_xxxzz_1[i] * fe_0 + ta_xxxy_xxxzz_0[i] * pa_y[i] - ta_xxxy_xxxzz_1[i] * pc_y[i];

        ta_xxxyy_xxyyy_0[i] = 2.0 * ta_xyy_xxyyy_0[i] * fe_0 - 2.0 * ta_xyy_xxyyy_1[i] * fe_0 + 2.0 * ta_xxyy_xyyy_0[i] * fe_0 -
                              2.0 * ta_xxyy_xyyy_1[i] * fe_0 + ta_xxyy_xxyyy_0[i] * pa_x[i] - ta_xxyy_xxyyy_1[i] * pc_x[i];

        ta_xxxyy_xxyyz_0[i] = 2.0 * ta_xyy_xxyyz_0[i] * fe_0 - 2.0 * ta_xyy_xxyyz_1[i] * fe_0 + 2.0 * ta_xxyy_xyyz_0[i] * fe_0 -
                              2.0 * ta_xxyy_xyyz_1[i] * fe_0 + ta_xxyy_xxyyz_0[i] * pa_x[i] - ta_xxyy_xxyyz_1[i] * pc_x[i];

        ta_xxxyy_xxyzz_0[i] = 2.0 * ta_xyy_xxyzz_0[i] * fe_0 - 2.0 * ta_xyy_xxyzz_1[i] * fe_0 + 2.0 * ta_xxyy_xyzz_0[i] * fe_0 -
                              2.0 * ta_xxyy_xyzz_1[i] * fe_0 + ta_xxyy_xxyzz_0[i] * pa_x[i] - ta_xxyy_xxyzz_1[i] * pc_x[i];

        ta_xxxyy_xxzzz_0[i] = ta_xxx_xxzzz_0[i] * fe_0 - ta_xxx_xxzzz_1[i] * fe_0 + ta_xxxy_xxzzz_0[i] * pa_y[i] - ta_xxxy_xxzzz_1[i] * pc_y[i];

        ta_xxxyy_xyyyy_0[i] = 2.0 * ta_xyy_xyyyy_0[i] * fe_0 - 2.0 * ta_xyy_xyyyy_1[i] * fe_0 + ta_xxyy_yyyy_0[i] * fe_0 - ta_xxyy_yyyy_1[i] * fe_0 +
                              ta_xxyy_xyyyy_0[i] * pa_x[i] - ta_xxyy_xyyyy_1[i] * pc_x[i];

        ta_xxxyy_xyyyz_0[i] = 2.0 * ta_xyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xyy_xyyyz_1[i] * fe_0 + ta_xxyy_yyyz_0[i] * fe_0 - ta_xxyy_yyyz_1[i] * fe_0 +
                              ta_xxyy_xyyyz_0[i] * pa_x[i] - ta_xxyy_xyyyz_1[i] * pc_x[i];

        ta_xxxyy_xyyzz_0[i] = 2.0 * ta_xyy_xyyzz_0[i] * fe_0 - 2.0 * ta_xyy_xyyzz_1[i] * fe_0 + ta_xxyy_yyzz_0[i] * fe_0 - ta_xxyy_yyzz_1[i] * fe_0 +
                              ta_xxyy_xyyzz_0[i] * pa_x[i] - ta_xxyy_xyyzz_1[i] * pc_x[i];

        ta_xxxyy_xyzzz_0[i] = 2.0 * ta_xyy_xyzzz_0[i] * fe_0 - 2.0 * ta_xyy_xyzzz_1[i] * fe_0 + ta_xxyy_yzzz_0[i] * fe_0 - ta_xxyy_yzzz_1[i] * fe_0 +
                              ta_xxyy_xyzzz_0[i] * pa_x[i] - ta_xxyy_xyzzz_1[i] * pc_x[i];

        ta_xxxyy_xzzzz_0[i] = ta_xxx_xzzzz_0[i] * fe_0 - ta_xxx_xzzzz_1[i] * fe_0 + ta_xxxy_xzzzz_0[i] * pa_y[i] - ta_xxxy_xzzzz_1[i] * pc_y[i];

        ta_xxxyy_yyyyy_0[i] =
            2.0 * ta_xyy_yyyyy_0[i] * fe_0 - 2.0 * ta_xyy_yyyyy_1[i] * fe_0 + ta_xxyy_yyyyy_0[i] * pa_x[i] - ta_xxyy_yyyyy_1[i] * pc_x[i];

        ta_xxxyy_yyyyz_0[i] =
            2.0 * ta_xyy_yyyyz_0[i] * fe_0 - 2.0 * ta_xyy_yyyyz_1[i] * fe_0 + ta_xxyy_yyyyz_0[i] * pa_x[i] - ta_xxyy_yyyyz_1[i] * pc_x[i];

        ta_xxxyy_yyyzz_0[i] =
            2.0 * ta_xyy_yyyzz_0[i] * fe_0 - 2.0 * ta_xyy_yyyzz_1[i] * fe_0 + ta_xxyy_yyyzz_0[i] * pa_x[i] - ta_xxyy_yyyzz_1[i] * pc_x[i];

        ta_xxxyy_yyzzz_0[i] =
            2.0 * ta_xyy_yyzzz_0[i] * fe_0 - 2.0 * ta_xyy_yyzzz_1[i] * fe_0 + ta_xxyy_yyzzz_0[i] * pa_x[i] - ta_xxyy_yyzzz_1[i] * pc_x[i];

        ta_xxxyy_yzzzz_0[i] =
            2.0 * ta_xyy_yzzzz_0[i] * fe_0 - 2.0 * ta_xyy_yzzzz_1[i] * fe_0 + ta_xxyy_yzzzz_0[i] * pa_x[i] - ta_xxyy_yzzzz_1[i] * pc_x[i];

        ta_xxxyy_zzzzz_0[i] =
            2.0 * ta_xyy_zzzzz_0[i] * fe_0 - 2.0 * ta_xyy_zzzzz_1[i] * fe_0 + ta_xxyy_zzzzz_0[i] * pa_x[i] - ta_xxyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 84-105 components of targeted buffer : HH

    auto ta_xxxyz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 84);

    auto ta_xxxyz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 85);

    auto ta_xxxyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 86);

    auto ta_xxxyz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 87);

    auto ta_xxxyz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 88);

    auto ta_xxxyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 89);

    auto ta_xxxyz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 90);

    auto ta_xxxyz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 91);

    auto ta_xxxyz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 92);

    auto ta_xxxyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 93);

    auto ta_xxxyz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 94);

    auto ta_xxxyz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 95);

    auto ta_xxxyz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 96);

    auto ta_xxxyz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 97);

    auto ta_xxxyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 98);

    auto ta_xxxyz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 99);

    auto ta_xxxyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 100);

    auto ta_xxxyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 101);

    auto ta_xxxyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 102);

    auto ta_xxxyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 103);

    auto ta_xxxyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 104);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta_xxxy_xxxxy_0,  \
                             ta_xxxy_xxxxy_1,  \
                             ta_xxxy_xxxyy_0,  \
                             ta_xxxy_xxxyy_1,  \
                             ta_xxxy_xxyyy_0,  \
                             ta_xxxy_xxyyy_1,  \
                             ta_xxxy_xyyyy_0,  \
                             ta_xxxy_xyyyy_1,  \
                             ta_xxxy_yyyyy_0,  \
                             ta_xxxy_yyyyy_1,  \
                             ta_xxxyz_xxxxx_0, \
                             ta_xxxyz_xxxxy_0, \
                             ta_xxxyz_xxxxz_0, \
                             ta_xxxyz_xxxyy_0, \
                             ta_xxxyz_xxxyz_0, \
                             ta_xxxyz_xxxzz_0, \
                             ta_xxxyz_xxyyy_0, \
                             ta_xxxyz_xxyyz_0, \
                             ta_xxxyz_xxyzz_0, \
                             ta_xxxyz_xxzzz_0, \
                             ta_xxxyz_xyyyy_0, \
                             ta_xxxyz_xyyyz_0, \
                             ta_xxxyz_xyyzz_0, \
                             ta_xxxyz_xyzzz_0, \
                             ta_xxxyz_xzzzz_0, \
                             ta_xxxyz_yyyyy_0, \
                             ta_xxxyz_yyyyz_0, \
                             ta_xxxyz_yyyzz_0, \
                             ta_xxxyz_yyzzz_0, \
                             ta_xxxyz_yzzzz_0, \
                             ta_xxxyz_zzzzz_0, \
                             ta_xxxz_xxxxx_0,  \
                             ta_xxxz_xxxxx_1,  \
                             ta_xxxz_xxxxz_0,  \
                             ta_xxxz_xxxxz_1,  \
                             ta_xxxz_xxxyz_0,  \
                             ta_xxxz_xxxyz_1,  \
                             ta_xxxz_xxxz_0,   \
                             ta_xxxz_xxxz_1,   \
                             ta_xxxz_xxxzz_0,  \
                             ta_xxxz_xxxzz_1,  \
                             ta_xxxz_xxyyz_0,  \
                             ta_xxxz_xxyyz_1,  \
                             ta_xxxz_xxyz_0,   \
                             ta_xxxz_xxyz_1,   \
                             ta_xxxz_xxyzz_0,  \
                             ta_xxxz_xxyzz_1,  \
                             ta_xxxz_xxzz_0,   \
                             ta_xxxz_xxzz_1,   \
                             ta_xxxz_xxzzz_0,  \
                             ta_xxxz_xxzzz_1,  \
                             ta_xxxz_xyyyz_0,  \
                             ta_xxxz_xyyyz_1,  \
                             ta_xxxz_xyyz_0,   \
                             ta_xxxz_xyyz_1,   \
                             ta_xxxz_xyyzz_0,  \
                             ta_xxxz_xyyzz_1,  \
                             ta_xxxz_xyzz_0,   \
                             ta_xxxz_xyzz_1,   \
                             ta_xxxz_xyzzz_0,  \
                             ta_xxxz_xyzzz_1,  \
                             ta_xxxz_xzzz_0,   \
                             ta_xxxz_xzzz_1,   \
                             ta_xxxz_xzzzz_0,  \
                             ta_xxxz_xzzzz_1,  \
                             ta_xxxz_zzzzz_0,  \
                             ta_xxxz_zzzzz_1,  \
                             ta_xxyz_yyyyz_0,  \
                             ta_xxyz_yyyyz_1,  \
                             ta_xxyz_yyyzz_0,  \
                             ta_xxyz_yyyzz_1,  \
                             ta_xxyz_yyzzz_0,  \
                             ta_xxyz_yyzzz_1,  \
                             ta_xxyz_yzzzz_0,  \
                             ta_xxyz_yzzzz_1,  \
                             ta_xyz_yyyyz_0,   \
                             ta_xyz_yyyyz_1,   \
                             ta_xyz_yyyzz_0,   \
                             ta_xyz_yyyzz_1,   \
                             ta_xyz_yyzzz_0,   \
                             ta_xyz_yyzzz_1,   \
                             ta_xyz_yzzzz_0,   \
                             ta_xyz_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyz_xxxxx_0[i] = ta_xxxz_xxxxx_0[i] * pa_y[i] - ta_xxxz_xxxxx_1[i] * pc_y[i];

        ta_xxxyz_xxxxy_0[i] = ta_xxxy_xxxxy_0[i] * pa_z[i] - ta_xxxy_xxxxy_1[i] * pc_z[i];

        ta_xxxyz_xxxxz_0[i] = ta_xxxz_xxxxz_0[i] * pa_y[i] - ta_xxxz_xxxxz_1[i] * pc_y[i];

        ta_xxxyz_xxxyy_0[i] = ta_xxxy_xxxyy_0[i] * pa_z[i] - ta_xxxy_xxxyy_1[i] * pc_z[i];

        ta_xxxyz_xxxyz_0[i] = ta_xxxz_xxxz_0[i] * fe_0 - ta_xxxz_xxxz_1[i] * fe_0 + ta_xxxz_xxxyz_0[i] * pa_y[i] - ta_xxxz_xxxyz_1[i] * pc_y[i];

        ta_xxxyz_xxxzz_0[i] = ta_xxxz_xxxzz_0[i] * pa_y[i] - ta_xxxz_xxxzz_1[i] * pc_y[i];

        ta_xxxyz_xxyyy_0[i] = ta_xxxy_xxyyy_0[i] * pa_z[i] - ta_xxxy_xxyyy_1[i] * pc_z[i];

        ta_xxxyz_xxyyz_0[i] =
            2.0 * ta_xxxz_xxyz_0[i] * fe_0 - 2.0 * ta_xxxz_xxyz_1[i] * fe_0 + ta_xxxz_xxyyz_0[i] * pa_y[i] - ta_xxxz_xxyyz_1[i] * pc_y[i];

        ta_xxxyz_xxyzz_0[i] = ta_xxxz_xxzz_0[i] * fe_0 - ta_xxxz_xxzz_1[i] * fe_0 + ta_xxxz_xxyzz_0[i] * pa_y[i] - ta_xxxz_xxyzz_1[i] * pc_y[i];

        ta_xxxyz_xxzzz_0[i] = ta_xxxz_xxzzz_0[i] * pa_y[i] - ta_xxxz_xxzzz_1[i] * pc_y[i];

        ta_xxxyz_xyyyy_0[i] = ta_xxxy_xyyyy_0[i] * pa_z[i] - ta_xxxy_xyyyy_1[i] * pc_z[i];

        ta_xxxyz_xyyyz_0[i] =
            3.0 * ta_xxxz_xyyz_0[i] * fe_0 - 3.0 * ta_xxxz_xyyz_1[i] * fe_0 + ta_xxxz_xyyyz_0[i] * pa_y[i] - ta_xxxz_xyyyz_1[i] * pc_y[i];

        ta_xxxyz_xyyzz_0[i] =
            2.0 * ta_xxxz_xyzz_0[i] * fe_0 - 2.0 * ta_xxxz_xyzz_1[i] * fe_0 + ta_xxxz_xyyzz_0[i] * pa_y[i] - ta_xxxz_xyyzz_1[i] * pc_y[i];

        ta_xxxyz_xyzzz_0[i] = ta_xxxz_xzzz_0[i] * fe_0 - ta_xxxz_xzzz_1[i] * fe_0 + ta_xxxz_xyzzz_0[i] * pa_y[i] - ta_xxxz_xyzzz_1[i] * pc_y[i];

        ta_xxxyz_xzzzz_0[i] = ta_xxxz_xzzzz_0[i] * pa_y[i] - ta_xxxz_xzzzz_1[i] * pc_y[i];

        ta_xxxyz_yyyyy_0[i] = ta_xxxy_yyyyy_0[i] * pa_z[i] - ta_xxxy_yyyyy_1[i] * pc_z[i];

        ta_xxxyz_yyyyz_0[i] =
            2.0 * ta_xyz_yyyyz_0[i] * fe_0 - 2.0 * ta_xyz_yyyyz_1[i] * fe_0 + ta_xxyz_yyyyz_0[i] * pa_x[i] - ta_xxyz_yyyyz_1[i] * pc_x[i];

        ta_xxxyz_yyyzz_0[i] =
            2.0 * ta_xyz_yyyzz_0[i] * fe_0 - 2.0 * ta_xyz_yyyzz_1[i] * fe_0 + ta_xxyz_yyyzz_0[i] * pa_x[i] - ta_xxyz_yyyzz_1[i] * pc_x[i];

        ta_xxxyz_yyzzz_0[i] =
            2.0 * ta_xyz_yyzzz_0[i] * fe_0 - 2.0 * ta_xyz_yyzzz_1[i] * fe_0 + ta_xxyz_yyzzz_0[i] * pa_x[i] - ta_xxyz_yyzzz_1[i] * pc_x[i];

        ta_xxxyz_yzzzz_0[i] =
            2.0 * ta_xyz_yzzzz_0[i] * fe_0 - 2.0 * ta_xyz_yzzzz_1[i] * fe_0 + ta_xxyz_yzzzz_0[i] * pa_x[i] - ta_xxyz_yzzzz_1[i] * pc_x[i];

        ta_xxxyz_zzzzz_0[i] = ta_xxxz_zzzzz_0[i] * pa_y[i] - ta_xxxz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : HH

    auto ta_xxxzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 105);

    auto ta_xxxzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 106);

    auto ta_xxxzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 107);

    auto ta_xxxzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 108);

    auto ta_xxxzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 109);

    auto ta_xxxzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 110);

    auto ta_xxxzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 111);

    auto ta_xxxzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 112);

    auto ta_xxxzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 113);

    auto ta_xxxzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 114);

    auto ta_xxxzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 115);

    auto ta_xxxzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 116);

    auto ta_xxxzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 117);

    auto ta_xxxzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 118);

    auto ta_xxxzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 119);

    auto ta_xxxzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 120);

    auto ta_xxxzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 121);

    auto ta_xxxzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 122);

    auto ta_xxxzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 123);

    auto ta_xxxzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 124);

    auto ta_xxxzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 125);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxx_xxxxx_0,   \
                             ta_xxx_xxxxx_1,   \
                             ta_xxx_xxxxy_0,   \
                             ta_xxx_xxxxy_1,   \
                             ta_xxx_xxxyy_0,   \
                             ta_xxx_xxxyy_1,   \
                             ta_xxx_xxyyy_0,   \
                             ta_xxx_xxyyy_1,   \
                             ta_xxx_xyyyy_0,   \
                             ta_xxx_xyyyy_1,   \
                             ta_xxxz_xxxxx_0,  \
                             ta_xxxz_xxxxx_1,  \
                             ta_xxxz_xxxxy_0,  \
                             ta_xxxz_xxxxy_1,  \
                             ta_xxxz_xxxyy_0,  \
                             ta_xxxz_xxxyy_1,  \
                             ta_xxxz_xxyyy_0,  \
                             ta_xxxz_xxyyy_1,  \
                             ta_xxxz_xyyyy_0,  \
                             ta_xxxz_xyyyy_1,  \
                             ta_xxxzz_xxxxx_0, \
                             ta_xxxzz_xxxxy_0, \
                             ta_xxxzz_xxxxz_0, \
                             ta_xxxzz_xxxyy_0, \
                             ta_xxxzz_xxxyz_0, \
                             ta_xxxzz_xxxzz_0, \
                             ta_xxxzz_xxyyy_0, \
                             ta_xxxzz_xxyyz_0, \
                             ta_xxxzz_xxyzz_0, \
                             ta_xxxzz_xxzzz_0, \
                             ta_xxxzz_xyyyy_0, \
                             ta_xxxzz_xyyyz_0, \
                             ta_xxxzz_xyyzz_0, \
                             ta_xxxzz_xyzzz_0, \
                             ta_xxxzz_xzzzz_0, \
                             ta_xxxzz_yyyyy_0, \
                             ta_xxxzz_yyyyz_0, \
                             ta_xxxzz_yyyzz_0, \
                             ta_xxxzz_yyzzz_0, \
                             ta_xxxzz_yzzzz_0, \
                             ta_xxxzz_zzzzz_0, \
                             ta_xxzz_xxxxz_0,  \
                             ta_xxzz_xxxxz_1,  \
                             ta_xxzz_xxxyz_0,  \
                             ta_xxzz_xxxyz_1,  \
                             ta_xxzz_xxxz_0,   \
                             ta_xxzz_xxxz_1,   \
                             ta_xxzz_xxxzz_0,  \
                             ta_xxzz_xxxzz_1,  \
                             ta_xxzz_xxyyz_0,  \
                             ta_xxzz_xxyyz_1,  \
                             ta_xxzz_xxyz_0,   \
                             ta_xxzz_xxyz_1,   \
                             ta_xxzz_xxyzz_0,  \
                             ta_xxzz_xxyzz_1,  \
                             ta_xxzz_xxzz_0,   \
                             ta_xxzz_xxzz_1,   \
                             ta_xxzz_xxzzz_0,  \
                             ta_xxzz_xxzzz_1,  \
                             ta_xxzz_xyyyz_0,  \
                             ta_xxzz_xyyyz_1,  \
                             ta_xxzz_xyyz_0,   \
                             ta_xxzz_xyyz_1,   \
                             ta_xxzz_xyyzz_0,  \
                             ta_xxzz_xyyzz_1,  \
                             ta_xxzz_xyzz_0,   \
                             ta_xxzz_xyzz_1,   \
                             ta_xxzz_xyzzz_0,  \
                             ta_xxzz_xyzzz_1,  \
                             ta_xxzz_xzzz_0,   \
                             ta_xxzz_xzzz_1,   \
                             ta_xxzz_xzzzz_0,  \
                             ta_xxzz_xzzzz_1,  \
                             ta_xxzz_yyyyy_0,  \
                             ta_xxzz_yyyyy_1,  \
                             ta_xxzz_yyyyz_0,  \
                             ta_xxzz_yyyyz_1,  \
                             ta_xxzz_yyyz_0,   \
                             ta_xxzz_yyyz_1,   \
                             ta_xxzz_yyyzz_0,  \
                             ta_xxzz_yyyzz_1,  \
                             ta_xxzz_yyzz_0,   \
                             ta_xxzz_yyzz_1,   \
                             ta_xxzz_yyzzz_0,  \
                             ta_xxzz_yyzzz_1,  \
                             ta_xxzz_yzzz_0,   \
                             ta_xxzz_yzzz_1,   \
                             ta_xxzz_yzzzz_0,  \
                             ta_xxzz_yzzzz_1,  \
                             ta_xxzz_zzzz_0,   \
                             ta_xxzz_zzzz_1,   \
                             ta_xxzz_zzzzz_0,  \
                             ta_xxzz_zzzzz_1,  \
                             ta_xzz_xxxxz_0,   \
                             ta_xzz_xxxxz_1,   \
                             ta_xzz_xxxyz_0,   \
                             ta_xzz_xxxyz_1,   \
                             ta_xzz_xxxzz_0,   \
                             ta_xzz_xxxzz_1,   \
                             ta_xzz_xxyyz_0,   \
                             ta_xzz_xxyyz_1,   \
                             ta_xzz_xxyzz_0,   \
                             ta_xzz_xxyzz_1,   \
                             ta_xzz_xxzzz_0,   \
                             ta_xzz_xxzzz_1,   \
                             ta_xzz_xyyyz_0,   \
                             ta_xzz_xyyyz_1,   \
                             ta_xzz_xyyzz_0,   \
                             ta_xzz_xyyzz_1,   \
                             ta_xzz_xyzzz_0,   \
                             ta_xzz_xyzzz_1,   \
                             ta_xzz_xzzzz_0,   \
                             ta_xzz_xzzzz_1,   \
                             ta_xzz_yyyyy_0,   \
                             ta_xzz_yyyyy_1,   \
                             ta_xzz_yyyyz_0,   \
                             ta_xzz_yyyyz_1,   \
                             ta_xzz_yyyzz_0,   \
                             ta_xzz_yyyzz_1,   \
                             ta_xzz_yyzzz_0,   \
                             ta_xzz_yyzzz_1,   \
                             ta_xzz_yzzzz_0,   \
                             ta_xzz_yzzzz_1,   \
                             ta_xzz_zzzzz_0,   \
                             ta_xzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzz_xxxxx_0[i] = ta_xxx_xxxxx_0[i] * fe_0 - ta_xxx_xxxxx_1[i] * fe_0 + ta_xxxz_xxxxx_0[i] * pa_z[i] - ta_xxxz_xxxxx_1[i] * pc_z[i];

        ta_xxxzz_xxxxy_0[i] = ta_xxx_xxxxy_0[i] * fe_0 - ta_xxx_xxxxy_1[i] * fe_0 + ta_xxxz_xxxxy_0[i] * pa_z[i] - ta_xxxz_xxxxy_1[i] * pc_z[i];

        ta_xxxzz_xxxxz_0[i] = 2.0 * ta_xzz_xxxxz_0[i] * fe_0 - 2.0 * ta_xzz_xxxxz_1[i] * fe_0 + 4.0 * ta_xxzz_xxxz_0[i] * fe_0 -
                              4.0 * ta_xxzz_xxxz_1[i] * fe_0 + ta_xxzz_xxxxz_0[i] * pa_x[i] - ta_xxzz_xxxxz_1[i] * pc_x[i];

        ta_xxxzz_xxxyy_0[i] = ta_xxx_xxxyy_0[i] * fe_0 - ta_xxx_xxxyy_1[i] * fe_0 + ta_xxxz_xxxyy_0[i] * pa_z[i] - ta_xxxz_xxxyy_1[i] * pc_z[i];

        ta_xxxzz_xxxyz_0[i] = 2.0 * ta_xzz_xxxyz_0[i] * fe_0 - 2.0 * ta_xzz_xxxyz_1[i] * fe_0 + 3.0 * ta_xxzz_xxyz_0[i] * fe_0 -
                              3.0 * ta_xxzz_xxyz_1[i] * fe_0 + ta_xxzz_xxxyz_0[i] * pa_x[i] - ta_xxzz_xxxyz_1[i] * pc_x[i];

        ta_xxxzz_xxxzz_0[i] = 2.0 * ta_xzz_xxxzz_0[i] * fe_0 - 2.0 * ta_xzz_xxxzz_1[i] * fe_0 + 3.0 * ta_xxzz_xxzz_0[i] * fe_0 -
                              3.0 * ta_xxzz_xxzz_1[i] * fe_0 + ta_xxzz_xxxzz_0[i] * pa_x[i] - ta_xxzz_xxxzz_1[i] * pc_x[i];

        ta_xxxzz_xxyyy_0[i] = ta_xxx_xxyyy_0[i] * fe_0 - ta_xxx_xxyyy_1[i] * fe_0 + ta_xxxz_xxyyy_0[i] * pa_z[i] - ta_xxxz_xxyyy_1[i] * pc_z[i];

        ta_xxxzz_xxyyz_0[i] = 2.0 * ta_xzz_xxyyz_0[i] * fe_0 - 2.0 * ta_xzz_xxyyz_1[i] * fe_0 + 2.0 * ta_xxzz_xyyz_0[i] * fe_0 -
                              2.0 * ta_xxzz_xyyz_1[i] * fe_0 + ta_xxzz_xxyyz_0[i] * pa_x[i] - ta_xxzz_xxyyz_1[i] * pc_x[i];

        ta_xxxzz_xxyzz_0[i] = 2.0 * ta_xzz_xxyzz_0[i] * fe_0 - 2.0 * ta_xzz_xxyzz_1[i] * fe_0 + 2.0 * ta_xxzz_xyzz_0[i] * fe_0 -
                              2.0 * ta_xxzz_xyzz_1[i] * fe_0 + ta_xxzz_xxyzz_0[i] * pa_x[i] - ta_xxzz_xxyzz_1[i] * pc_x[i];

        ta_xxxzz_xxzzz_0[i] = 2.0 * ta_xzz_xxzzz_0[i] * fe_0 - 2.0 * ta_xzz_xxzzz_1[i] * fe_0 + 2.0 * ta_xxzz_xzzz_0[i] * fe_0 -
                              2.0 * ta_xxzz_xzzz_1[i] * fe_0 + ta_xxzz_xxzzz_0[i] * pa_x[i] - ta_xxzz_xxzzz_1[i] * pc_x[i];

        ta_xxxzz_xyyyy_0[i] = ta_xxx_xyyyy_0[i] * fe_0 - ta_xxx_xyyyy_1[i] * fe_0 + ta_xxxz_xyyyy_0[i] * pa_z[i] - ta_xxxz_xyyyy_1[i] * pc_z[i];

        ta_xxxzz_xyyyz_0[i] = 2.0 * ta_xzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xzz_xyyyz_1[i] * fe_0 + ta_xxzz_yyyz_0[i] * fe_0 - ta_xxzz_yyyz_1[i] * fe_0 +
                              ta_xxzz_xyyyz_0[i] * pa_x[i] - ta_xxzz_xyyyz_1[i] * pc_x[i];

        ta_xxxzz_xyyzz_0[i] = 2.0 * ta_xzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xzz_xyyzz_1[i] * fe_0 + ta_xxzz_yyzz_0[i] * fe_0 - ta_xxzz_yyzz_1[i] * fe_0 +
                              ta_xxzz_xyyzz_0[i] * pa_x[i] - ta_xxzz_xyyzz_1[i] * pc_x[i];

        ta_xxxzz_xyzzz_0[i] = 2.0 * ta_xzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xzz_xyzzz_1[i] * fe_0 + ta_xxzz_yzzz_0[i] * fe_0 - ta_xxzz_yzzz_1[i] * fe_0 +
                              ta_xxzz_xyzzz_0[i] * pa_x[i] - ta_xxzz_xyzzz_1[i] * pc_x[i];

        ta_xxxzz_xzzzz_0[i] = 2.0 * ta_xzz_xzzzz_0[i] * fe_0 - 2.0 * ta_xzz_xzzzz_1[i] * fe_0 + ta_xxzz_zzzz_0[i] * fe_0 - ta_xxzz_zzzz_1[i] * fe_0 +
                              ta_xxzz_xzzzz_0[i] * pa_x[i] - ta_xxzz_xzzzz_1[i] * pc_x[i];

        ta_xxxzz_yyyyy_0[i] =
            2.0 * ta_xzz_yyyyy_0[i] * fe_0 - 2.0 * ta_xzz_yyyyy_1[i] * fe_0 + ta_xxzz_yyyyy_0[i] * pa_x[i] - ta_xxzz_yyyyy_1[i] * pc_x[i];

        ta_xxxzz_yyyyz_0[i] =
            2.0 * ta_xzz_yyyyz_0[i] * fe_0 - 2.0 * ta_xzz_yyyyz_1[i] * fe_0 + ta_xxzz_yyyyz_0[i] * pa_x[i] - ta_xxzz_yyyyz_1[i] * pc_x[i];

        ta_xxxzz_yyyzz_0[i] =
            2.0 * ta_xzz_yyyzz_0[i] * fe_0 - 2.0 * ta_xzz_yyyzz_1[i] * fe_0 + ta_xxzz_yyyzz_0[i] * pa_x[i] - ta_xxzz_yyyzz_1[i] * pc_x[i];

        ta_xxxzz_yyzzz_0[i] =
            2.0 * ta_xzz_yyzzz_0[i] * fe_0 - 2.0 * ta_xzz_yyzzz_1[i] * fe_0 + ta_xxzz_yyzzz_0[i] * pa_x[i] - ta_xxzz_yyzzz_1[i] * pc_x[i];

        ta_xxxzz_yzzzz_0[i] =
            2.0 * ta_xzz_yzzzz_0[i] * fe_0 - 2.0 * ta_xzz_yzzzz_1[i] * fe_0 + ta_xxzz_yzzzz_0[i] * pa_x[i] - ta_xxzz_yzzzz_1[i] * pc_x[i];

        ta_xxxzz_zzzzz_0[i] =
            2.0 * ta_xzz_zzzzz_0[i] * fe_0 - 2.0 * ta_xzz_zzzzz_1[i] * fe_0 + ta_xxzz_zzzzz_0[i] * pa_x[i] - ta_xxzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 126-147 components of targeted buffer : HH

    auto ta_xxyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 126);

    auto ta_xxyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 127);

    auto ta_xxyyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 128);

    auto ta_xxyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 129);

    auto ta_xxyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 130);

    auto ta_xxyyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 131);

    auto ta_xxyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 132);

    auto ta_xxyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 133);

    auto ta_xxyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 134);

    auto ta_xxyyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 135);

    auto ta_xxyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 136);

    auto ta_xxyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 137);

    auto ta_xxyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 138);

    auto ta_xxyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 139);

    auto ta_xxyyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 140);

    auto ta_xxyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 141);

    auto ta_xxyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 142);

    auto ta_xxyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 143);

    auto ta_xxyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 144);

    auto ta_xxyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 145);

    auto ta_xxyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 146);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxy_xxxxx_0,   \
                             ta_xxy_xxxxx_1,   \
                             ta_xxy_xxxxz_0,   \
                             ta_xxy_xxxxz_1,   \
                             ta_xxy_xxxzz_0,   \
                             ta_xxy_xxxzz_1,   \
                             ta_xxy_xxzzz_0,   \
                             ta_xxy_xxzzz_1,   \
                             ta_xxy_xzzzz_0,   \
                             ta_xxy_xzzzz_1,   \
                             ta_xxyy_xxxxx_0,  \
                             ta_xxyy_xxxxx_1,  \
                             ta_xxyy_xxxxz_0,  \
                             ta_xxyy_xxxxz_1,  \
                             ta_xxyy_xxxzz_0,  \
                             ta_xxyy_xxxzz_1,  \
                             ta_xxyy_xxzzz_0,  \
                             ta_xxyy_xxzzz_1,  \
                             ta_xxyy_xzzzz_0,  \
                             ta_xxyy_xzzzz_1,  \
                             ta_xxyyy_xxxxx_0, \
                             ta_xxyyy_xxxxy_0, \
                             ta_xxyyy_xxxxz_0, \
                             ta_xxyyy_xxxyy_0, \
                             ta_xxyyy_xxxyz_0, \
                             ta_xxyyy_xxxzz_0, \
                             ta_xxyyy_xxyyy_0, \
                             ta_xxyyy_xxyyz_0, \
                             ta_xxyyy_xxyzz_0, \
                             ta_xxyyy_xxzzz_0, \
                             ta_xxyyy_xyyyy_0, \
                             ta_xxyyy_xyyyz_0, \
                             ta_xxyyy_xyyzz_0, \
                             ta_xxyyy_xyzzz_0, \
                             ta_xxyyy_xzzzz_0, \
                             ta_xxyyy_yyyyy_0, \
                             ta_xxyyy_yyyyz_0, \
                             ta_xxyyy_yyyzz_0, \
                             ta_xxyyy_yyzzz_0, \
                             ta_xxyyy_yzzzz_0, \
                             ta_xxyyy_zzzzz_0, \
                             ta_xyyy_xxxxy_0,  \
                             ta_xyyy_xxxxy_1,  \
                             ta_xyyy_xxxy_0,   \
                             ta_xyyy_xxxy_1,   \
                             ta_xyyy_xxxyy_0,  \
                             ta_xyyy_xxxyy_1,  \
                             ta_xyyy_xxxyz_0,  \
                             ta_xyyy_xxxyz_1,  \
                             ta_xyyy_xxyy_0,   \
                             ta_xyyy_xxyy_1,   \
                             ta_xyyy_xxyyy_0,  \
                             ta_xyyy_xxyyy_1,  \
                             ta_xyyy_xxyyz_0,  \
                             ta_xyyy_xxyyz_1,  \
                             ta_xyyy_xxyz_0,   \
                             ta_xyyy_xxyz_1,   \
                             ta_xyyy_xxyzz_0,  \
                             ta_xyyy_xxyzz_1,  \
                             ta_xyyy_xyyy_0,   \
                             ta_xyyy_xyyy_1,   \
                             ta_xyyy_xyyyy_0,  \
                             ta_xyyy_xyyyy_1,  \
                             ta_xyyy_xyyyz_0,  \
                             ta_xyyy_xyyyz_1,  \
                             ta_xyyy_xyyz_0,   \
                             ta_xyyy_xyyz_1,   \
                             ta_xyyy_xyyzz_0,  \
                             ta_xyyy_xyyzz_1,  \
                             ta_xyyy_xyzz_0,   \
                             ta_xyyy_xyzz_1,   \
                             ta_xyyy_xyzzz_0,  \
                             ta_xyyy_xyzzz_1,  \
                             ta_xyyy_yyyy_0,   \
                             ta_xyyy_yyyy_1,   \
                             ta_xyyy_yyyyy_0,  \
                             ta_xyyy_yyyyy_1,  \
                             ta_xyyy_yyyyz_0,  \
                             ta_xyyy_yyyyz_1,  \
                             ta_xyyy_yyyz_0,   \
                             ta_xyyy_yyyz_1,   \
                             ta_xyyy_yyyzz_0,  \
                             ta_xyyy_yyyzz_1,  \
                             ta_xyyy_yyzz_0,   \
                             ta_xyyy_yyzz_1,   \
                             ta_xyyy_yyzzz_0,  \
                             ta_xyyy_yyzzz_1,  \
                             ta_xyyy_yzzz_0,   \
                             ta_xyyy_yzzz_1,   \
                             ta_xyyy_yzzzz_0,  \
                             ta_xyyy_yzzzz_1,  \
                             ta_xyyy_zzzzz_0,  \
                             ta_xyyy_zzzzz_1,  \
                             ta_yyy_xxxxy_0,   \
                             ta_yyy_xxxxy_1,   \
                             ta_yyy_xxxyy_0,   \
                             ta_yyy_xxxyy_1,   \
                             ta_yyy_xxxyz_0,   \
                             ta_yyy_xxxyz_1,   \
                             ta_yyy_xxyyy_0,   \
                             ta_yyy_xxyyy_1,   \
                             ta_yyy_xxyyz_0,   \
                             ta_yyy_xxyyz_1,   \
                             ta_yyy_xxyzz_0,   \
                             ta_yyy_xxyzz_1,   \
                             ta_yyy_xyyyy_0,   \
                             ta_yyy_xyyyy_1,   \
                             ta_yyy_xyyyz_0,   \
                             ta_yyy_xyyyz_1,   \
                             ta_yyy_xyyzz_0,   \
                             ta_yyy_xyyzz_1,   \
                             ta_yyy_xyzzz_0,   \
                             ta_yyy_xyzzz_1,   \
                             ta_yyy_yyyyy_0,   \
                             ta_yyy_yyyyy_1,   \
                             ta_yyy_yyyyz_0,   \
                             ta_yyy_yyyyz_1,   \
                             ta_yyy_yyyzz_0,   \
                             ta_yyy_yyyzz_1,   \
                             ta_yyy_yyzzz_0,   \
                             ta_yyy_yyzzz_1,   \
                             ta_yyy_yzzzz_0,   \
                             ta_yyy_yzzzz_1,   \
                             ta_yyy_zzzzz_0,   \
                             ta_yyy_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyy_xxxxx_0[i] =
            2.0 * ta_xxy_xxxxx_0[i] * fe_0 - 2.0 * ta_xxy_xxxxx_1[i] * fe_0 + ta_xxyy_xxxxx_0[i] * pa_y[i] - ta_xxyy_xxxxx_1[i] * pc_y[i];

        ta_xxyyy_xxxxy_0[i] = ta_yyy_xxxxy_0[i] * fe_0 - ta_yyy_xxxxy_1[i] * fe_0 + 4.0 * ta_xyyy_xxxy_0[i] * fe_0 - 4.0 * ta_xyyy_xxxy_1[i] * fe_0 +
                              ta_xyyy_xxxxy_0[i] * pa_x[i] - ta_xyyy_xxxxy_1[i] * pc_x[i];

        ta_xxyyy_xxxxz_0[i] =
            2.0 * ta_xxy_xxxxz_0[i] * fe_0 - 2.0 * ta_xxy_xxxxz_1[i] * fe_0 + ta_xxyy_xxxxz_0[i] * pa_y[i] - ta_xxyy_xxxxz_1[i] * pc_y[i];

        ta_xxyyy_xxxyy_0[i] = ta_yyy_xxxyy_0[i] * fe_0 - ta_yyy_xxxyy_1[i] * fe_0 + 3.0 * ta_xyyy_xxyy_0[i] * fe_0 - 3.0 * ta_xyyy_xxyy_1[i] * fe_0 +
                              ta_xyyy_xxxyy_0[i] * pa_x[i] - ta_xyyy_xxxyy_1[i] * pc_x[i];

        ta_xxyyy_xxxyz_0[i] = ta_yyy_xxxyz_0[i] * fe_0 - ta_yyy_xxxyz_1[i] * fe_0 + 3.0 * ta_xyyy_xxyz_0[i] * fe_0 - 3.0 * ta_xyyy_xxyz_1[i] * fe_0 +
                              ta_xyyy_xxxyz_0[i] * pa_x[i] - ta_xyyy_xxxyz_1[i] * pc_x[i];

        ta_xxyyy_xxxzz_0[i] =
            2.0 * ta_xxy_xxxzz_0[i] * fe_0 - 2.0 * ta_xxy_xxxzz_1[i] * fe_0 + ta_xxyy_xxxzz_0[i] * pa_y[i] - ta_xxyy_xxxzz_1[i] * pc_y[i];

        ta_xxyyy_xxyyy_0[i] = ta_yyy_xxyyy_0[i] * fe_0 - ta_yyy_xxyyy_1[i] * fe_0 + 2.0 * ta_xyyy_xyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xyyy_1[i] * fe_0 +
                              ta_xyyy_xxyyy_0[i] * pa_x[i] - ta_xyyy_xxyyy_1[i] * pc_x[i];

        ta_xxyyy_xxyyz_0[i] = ta_yyy_xxyyz_0[i] * fe_0 - ta_yyy_xxyyz_1[i] * fe_0 + 2.0 * ta_xyyy_xyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyz_1[i] * fe_0 +
                              ta_xyyy_xxyyz_0[i] * pa_x[i] - ta_xyyy_xxyyz_1[i] * pc_x[i];

        ta_xxyyy_xxyzz_0[i] = ta_yyy_xxyzz_0[i] * fe_0 - ta_yyy_xxyzz_1[i] * fe_0 + 2.0 * ta_xyyy_xyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyzz_1[i] * fe_0 +
                              ta_xyyy_xxyzz_0[i] * pa_x[i] - ta_xyyy_xxyzz_1[i] * pc_x[i];

        ta_xxyyy_xxzzz_0[i] =
            2.0 * ta_xxy_xxzzz_0[i] * fe_0 - 2.0 * ta_xxy_xxzzz_1[i] * fe_0 + ta_xxyy_xxzzz_0[i] * pa_y[i] - ta_xxyy_xxzzz_1[i] * pc_y[i];

        ta_xxyyy_xyyyy_0[i] = ta_yyy_xyyyy_0[i] * fe_0 - ta_yyy_xyyyy_1[i] * fe_0 + ta_xyyy_yyyy_0[i] * fe_0 - ta_xyyy_yyyy_1[i] * fe_0 +
                              ta_xyyy_xyyyy_0[i] * pa_x[i] - ta_xyyy_xyyyy_1[i] * pc_x[i];

        ta_xxyyy_xyyyz_0[i] = ta_yyy_xyyyz_0[i] * fe_0 - ta_yyy_xyyyz_1[i] * fe_0 + ta_xyyy_yyyz_0[i] * fe_0 - ta_xyyy_yyyz_1[i] * fe_0 +
                              ta_xyyy_xyyyz_0[i] * pa_x[i] - ta_xyyy_xyyyz_1[i] * pc_x[i];

        ta_xxyyy_xyyzz_0[i] = ta_yyy_xyyzz_0[i] * fe_0 - ta_yyy_xyyzz_1[i] * fe_0 + ta_xyyy_yyzz_0[i] * fe_0 - ta_xyyy_yyzz_1[i] * fe_0 +
                              ta_xyyy_xyyzz_0[i] * pa_x[i] - ta_xyyy_xyyzz_1[i] * pc_x[i];

        ta_xxyyy_xyzzz_0[i] = ta_yyy_xyzzz_0[i] * fe_0 - ta_yyy_xyzzz_1[i] * fe_0 + ta_xyyy_yzzz_0[i] * fe_0 - ta_xyyy_yzzz_1[i] * fe_0 +
                              ta_xyyy_xyzzz_0[i] * pa_x[i] - ta_xyyy_xyzzz_1[i] * pc_x[i];

        ta_xxyyy_xzzzz_0[i] =
            2.0 * ta_xxy_xzzzz_0[i] * fe_0 - 2.0 * ta_xxy_xzzzz_1[i] * fe_0 + ta_xxyy_xzzzz_0[i] * pa_y[i] - ta_xxyy_xzzzz_1[i] * pc_y[i];

        ta_xxyyy_yyyyy_0[i] = ta_yyy_yyyyy_0[i] * fe_0 - ta_yyy_yyyyy_1[i] * fe_0 + ta_xyyy_yyyyy_0[i] * pa_x[i] - ta_xyyy_yyyyy_1[i] * pc_x[i];

        ta_xxyyy_yyyyz_0[i] = ta_yyy_yyyyz_0[i] * fe_0 - ta_yyy_yyyyz_1[i] * fe_0 + ta_xyyy_yyyyz_0[i] * pa_x[i] - ta_xyyy_yyyyz_1[i] * pc_x[i];

        ta_xxyyy_yyyzz_0[i] = ta_yyy_yyyzz_0[i] * fe_0 - ta_yyy_yyyzz_1[i] * fe_0 + ta_xyyy_yyyzz_0[i] * pa_x[i] - ta_xyyy_yyyzz_1[i] * pc_x[i];

        ta_xxyyy_yyzzz_0[i] = ta_yyy_yyzzz_0[i] * fe_0 - ta_yyy_yyzzz_1[i] * fe_0 + ta_xyyy_yyzzz_0[i] * pa_x[i] - ta_xyyy_yyzzz_1[i] * pc_x[i];

        ta_xxyyy_yzzzz_0[i] = ta_yyy_yzzzz_0[i] * fe_0 - ta_yyy_yzzzz_1[i] * fe_0 + ta_xyyy_yzzzz_0[i] * pa_x[i] - ta_xyyy_yzzzz_1[i] * pc_x[i];

        ta_xxyyy_zzzzz_0[i] = ta_yyy_zzzzz_0[i] * fe_0 - ta_yyy_zzzzz_1[i] * fe_0 + ta_xyyy_zzzzz_0[i] * pa_x[i] - ta_xyyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 147-168 components of targeted buffer : HH

    auto ta_xxyyz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 147);

    auto ta_xxyyz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 148);

    auto ta_xxyyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 149);

    auto ta_xxyyz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 150);

    auto ta_xxyyz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 151);

    auto ta_xxyyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 152);

    auto ta_xxyyz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 153);

    auto ta_xxyyz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 154);

    auto ta_xxyyz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 155);

    auto ta_xxyyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 156);

    auto ta_xxyyz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 157);

    auto ta_xxyyz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 158);

    auto ta_xxyyz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 159);

    auto ta_xxyyz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 160);

    auto ta_xxyyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 161);

    auto ta_xxyyz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 162);

    auto ta_xxyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 163);

    auto ta_xxyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 164);

    auto ta_xxyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 165);

    auto ta_xxyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 166);

    auto ta_xxyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 167);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta_xxyy_xxxxx_0,  \
                             ta_xxyy_xxxxx_1,  \
                             ta_xxyy_xxxxy_0,  \
                             ta_xxyy_xxxxy_1,  \
                             ta_xxyy_xxxy_0,   \
                             ta_xxyy_xxxy_1,   \
                             ta_xxyy_xxxyy_0,  \
                             ta_xxyy_xxxyy_1,  \
                             ta_xxyy_xxxyz_0,  \
                             ta_xxyy_xxxyz_1,  \
                             ta_xxyy_xxyy_0,   \
                             ta_xxyy_xxyy_1,   \
                             ta_xxyy_xxyyy_0,  \
                             ta_xxyy_xxyyy_1,  \
                             ta_xxyy_xxyyz_0,  \
                             ta_xxyy_xxyyz_1,  \
                             ta_xxyy_xxyz_0,   \
                             ta_xxyy_xxyz_1,   \
                             ta_xxyy_xxyzz_0,  \
                             ta_xxyy_xxyzz_1,  \
                             ta_xxyy_xyyy_0,   \
                             ta_xxyy_xyyy_1,   \
                             ta_xxyy_xyyyy_0,  \
                             ta_xxyy_xyyyy_1,  \
                             ta_xxyy_xyyyz_0,  \
                             ta_xxyy_xyyyz_1,  \
                             ta_xxyy_xyyz_0,   \
                             ta_xxyy_xyyz_1,   \
                             ta_xxyy_xyyzz_0,  \
                             ta_xxyy_xyyzz_1,  \
                             ta_xxyy_xyzz_0,   \
                             ta_xxyy_xyzz_1,   \
                             ta_xxyy_xyzzz_0,  \
                             ta_xxyy_xyzzz_1,  \
                             ta_xxyy_yyyyy_0,  \
                             ta_xxyy_yyyyy_1,  \
                             ta_xxyyz_xxxxx_0, \
                             ta_xxyyz_xxxxy_0, \
                             ta_xxyyz_xxxxz_0, \
                             ta_xxyyz_xxxyy_0, \
                             ta_xxyyz_xxxyz_0, \
                             ta_xxyyz_xxxzz_0, \
                             ta_xxyyz_xxyyy_0, \
                             ta_xxyyz_xxyyz_0, \
                             ta_xxyyz_xxyzz_0, \
                             ta_xxyyz_xxzzz_0, \
                             ta_xxyyz_xyyyy_0, \
                             ta_xxyyz_xyyyz_0, \
                             ta_xxyyz_xyyzz_0, \
                             ta_xxyyz_xyzzz_0, \
                             ta_xxyyz_xzzzz_0, \
                             ta_xxyyz_yyyyy_0, \
                             ta_xxyyz_yyyyz_0, \
                             ta_xxyyz_yyyzz_0, \
                             ta_xxyyz_yyzzz_0, \
                             ta_xxyyz_yzzzz_0, \
                             ta_xxyyz_zzzzz_0, \
                             ta_xxyz_xxxxz_0,  \
                             ta_xxyz_xxxxz_1,  \
                             ta_xxyz_xxxzz_0,  \
                             ta_xxyz_xxxzz_1,  \
                             ta_xxyz_xxzzz_0,  \
                             ta_xxyz_xxzzz_1,  \
                             ta_xxyz_xzzzz_0,  \
                             ta_xxyz_xzzzz_1,  \
                             ta_xxz_xxxxz_0,   \
                             ta_xxz_xxxxz_1,   \
                             ta_xxz_xxxzz_0,   \
                             ta_xxz_xxxzz_1,   \
                             ta_xxz_xxzzz_0,   \
                             ta_xxz_xxzzz_1,   \
                             ta_xxz_xzzzz_0,   \
                             ta_xxz_xzzzz_1,   \
                             ta_xyyz_yyyyz_0,  \
                             ta_xyyz_yyyyz_1,  \
                             ta_xyyz_yyyzz_0,  \
                             ta_xyyz_yyyzz_1,  \
                             ta_xyyz_yyzzz_0,  \
                             ta_xyyz_yyzzz_1,  \
                             ta_xyyz_yzzzz_0,  \
                             ta_xyyz_yzzzz_1,  \
                             ta_xyyz_zzzzz_0,  \
                             ta_xyyz_zzzzz_1,  \
                             ta_yyz_yyyyz_0,   \
                             ta_yyz_yyyyz_1,   \
                             ta_yyz_yyyzz_0,   \
                             ta_yyz_yyyzz_1,   \
                             ta_yyz_yyzzz_0,   \
                             ta_yyz_yyzzz_1,   \
                             ta_yyz_yzzzz_0,   \
                             ta_yyz_yzzzz_1,   \
                             ta_yyz_zzzzz_0,   \
                             ta_yyz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyz_xxxxx_0[i] = ta_xxyy_xxxxx_0[i] * pa_z[i] - ta_xxyy_xxxxx_1[i] * pc_z[i];

        ta_xxyyz_xxxxy_0[i] = ta_xxyy_xxxxy_0[i] * pa_z[i] - ta_xxyy_xxxxy_1[i] * pc_z[i];

        ta_xxyyz_xxxxz_0[i] = ta_xxz_xxxxz_0[i] * fe_0 - ta_xxz_xxxxz_1[i] * fe_0 + ta_xxyz_xxxxz_0[i] * pa_y[i] - ta_xxyz_xxxxz_1[i] * pc_y[i];

        ta_xxyyz_xxxyy_0[i] = ta_xxyy_xxxyy_0[i] * pa_z[i] - ta_xxyy_xxxyy_1[i] * pc_z[i];

        ta_xxyyz_xxxyz_0[i] = ta_xxyy_xxxy_0[i] * fe_0 - ta_xxyy_xxxy_1[i] * fe_0 + ta_xxyy_xxxyz_0[i] * pa_z[i] - ta_xxyy_xxxyz_1[i] * pc_z[i];

        ta_xxyyz_xxxzz_0[i] = ta_xxz_xxxzz_0[i] * fe_0 - ta_xxz_xxxzz_1[i] * fe_0 + ta_xxyz_xxxzz_0[i] * pa_y[i] - ta_xxyz_xxxzz_1[i] * pc_y[i];

        ta_xxyyz_xxyyy_0[i] = ta_xxyy_xxyyy_0[i] * pa_z[i] - ta_xxyy_xxyyy_1[i] * pc_z[i];

        ta_xxyyz_xxyyz_0[i] = ta_xxyy_xxyy_0[i] * fe_0 - ta_xxyy_xxyy_1[i] * fe_0 + ta_xxyy_xxyyz_0[i] * pa_z[i] - ta_xxyy_xxyyz_1[i] * pc_z[i];

        ta_xxyyz_xxyzz_0[i] =
            2.0 * ta_xxyy_xxyz_0[i] * fe_0 - 2.0 * ta_xxyy_xxyz_1[i] * fe_0 + ta_xxyy_xxyzz_0[i] * pa_z[i] - ta_xxyy_xxyzz_1[i] * pc_z[i];

        ta_xxyyz_xxzzz_0[i] = ta_xxz_xxzzz_0[i] * fe_0 - ta_xxz_xxzzz_1[i] * fe_0 + ta_xxyz_xxzzz_0[i] * pa_y[i] - ta_xxyz_xxzzz_1[i] * pc_y[i];

        ta_xxyyz_xyyyy_0[i] = ta_xxyy_xyyyy_0[i] * pa_z[i] - ta_xxyy_xyyyy_1[i] * pc_z[i];

        ta_xxyyz_xyyyz_0[i] = ta_xxyy_xyyy_0[i] * fe_0 - ta_xxyy_xyyy_1[i] * fe_0 + ta_xxyy_xyyyz_0[i] * pa_z[i] - ta_xxyy_xyyyz_1[i] * pc_z[i];

        ta_xxyyz_xyyzz_0[i] =
            2.0 * ta_xxyy_xyyz_0[i] * fe_0 - 2.0 * ta_xxyy_xyyz_1[i] * fe_0 + ta_xxyy_xyyzz_0[i] * pa_z[i] - ta_xxyy_xyyzz_1[i] * pc_z[i];

        ta_xxyyz_xyzzz_0[i] =
            3.0 * ta_xxyy_xyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyzz_1[i] * fe_0 + ta_xxyy_xyzzz_0[i] * pa_z[i] - ta_xxyy_xyzzz_1[i] * pc_z[i];

        ta_xxyyz_xzzzz_0[i] = ta_xxz_xzzzz_0[i] * fe_0 - ta_xxz_xzzzz_1[i] * fe_0 + ta_xxyz_xzzzz_0[i] * pa_y[i] - ta_xxyz_xzzzz_1[i] * pc_y[i];

        ta_xxyyz_yyyyy_0[i] = ta_xxyy_yyyyy_0[i] * pa_z[i] - ta_xxyy_yyyyy_1[i] * pc_z[i];

        ta_xxyyz_yyyyz_0[i] = ta_yyz_yyyyz_0[i] * fe_0 - ta_yyz_yyyyz_1[i] * fe_0 + ta_xyyz_yyyyz_0[i] * pa_x[i] - ta_xyyz_yyyyz_1[i] * pc_x[i];

        ta_xxyyz_yyyzz_0[i] = ta_yyz_yyyzz_0[i] * fe_0 - ta_yyz_yyyzz_1[i] * fe_0 + ta_xyyz_yyyzz_0[i] * pa_x[i] - ta_xyyz_yyyzz_1[i] * pc_x[i];

        ta_xxyyz_yyzzz_0[i] = ta_yyz_yyzzz_0[i] * fe_0 - ta_yyz_yyzzz_1[i] * fe_0 + ta_xyyz_yyzzz_0[i] * pa_x[i] - ta_xyyz_yyzzz_1[i] * pc_x[i];

        ta_xxyyz_yzzzz_0[i] = ta_yyz_yzzzz_0[i] * fe_0 - ta_yyz_yzzzz_1[i] * fe_0 + ta_xyyz_yzzzz_0[i] * pa_x[i] - ta_xyyz_yzzzz_1[i] * pc_x[i];

        ta_xxyyz_zzzzz_0[i] = ta_yyz_zzzzz_0[i] * fe_0 - ta_yyz_zzzzz_1[i] * fe_0 + ta_xyyz_zzzzz_0[i] * pa_x[i] - ta_xyyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 168-189 components of targeted buffer : HH

    auto ta_xxyzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 168);

    auto ta_xxyzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 169);

    auto ta_xxyzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 170);

    auto ta_xxyzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 171);

    auto ta_xxyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 172);

    auto ta_xxyzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 173);

    auto ta_xxyzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 174);

    auto ta_xxyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 175);

    auto ta_xxyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 176);

    auto ta_xxyzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 177);

    auto ta_xxyzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 178);

    auto ta_xxyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 179);

    auto ta_xxyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 180);

    auto ta_xxyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 181);

    auto ta_xxyzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 182);

    auto ta_xxyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 183);

    auto ta_xxyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 184);

    auto ta_xxyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 185);

    auto ta_xxyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 186);

    auto ta_xxyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 187);

    auto ta_xxyzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 188);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xxyzz_xxxxx_0, \
                             ta_xxyzz_xxxxy_0, \
                             ta_xxyzz_xxxxz_0, \
                             ta_xxyzz_xxxyy_0, \
                             ta_xxyzz_xxxyz_0, \
                             ta_xxyzz_xxxzz_0, \
                             ta_xxyzz_xxyyy_0, \
                             ta_xxyzz_xxyyz_0, \
                             ta_xxyzz_xxyzz_0, \
                             ta_xxyzz_xxzzz_0, \
                             ta_xxyzz_xyyyy_0, \
                             ta_xxyzz_xyyyz_0, \
                             ta_xxyzz_xyyzz_0, \
                             ta_xxyzz_xyzzz_0, \
                             ta_xxyzz_xzzzz_0, \
                             ta_xxyzz_yyyyy_0, \
                             ta_xxyzz_yyyyz_0, \
                             ta_xxyzz_yyyzz_0, \
                             ta_xxyzz_yyzzz_0, \
                             ta_xxyzz_yzzzz_0, \
                             ta_xxyzz_zzzzz_0, \
                             ta_xxzz_xxxx_0,   \
                             ta_xxzz_xxxx_1,   \
                             ta_xxzz_xxxxx_0,  \
                             ta_xxzz_xxxxx_1,  \
                             ta_xxzz_xxxxy_0,  \
                             ta_xxzz_xxxxy_1,  \
                             ta_xxzz_xxxxz_0,  \
                             ta_xxzz_xxxxz_1,  \
                             ta_xxzz_xxxy_0,   \
                             ta_xxzz_xxxy_1,   \
                             ta_xxzz_xxxyy_0,  \
                             ta_xxzz_xxxyy_1,  \
                             ta_xxzz_xxxyz_0,  \
                             ta_xxzz_xxxyz_1,  \
                             ta_xxzz_xxxz_0,   \
                             ta_xxzz_xxxz_1,   \
                             ta_xxzz_xxxzz_0,  \
                             ta_xxzz_xxxzz_1,  \
                             ta_xxzz_xxyy_0,   \
                             ta_xxzz_xxyy_1,   \
                             ta_xxzz_xxyyy_0,  \
                             ta_xxzz_xxyyy_1,  \
                             ta_xxzz_xxyyz_0,  \
                             ta_xxzz_xxyyz_1,  \
                             ta_xxzz_xxyz_0,   \
                             ta_xxzz_xxyz_1,   \
                             ta_xxzz_xxyzz_0,  \
                             ta_xxzz_xxyzz_1,  \
                             ta_xxzz_xxzz_0,   \
                             ta_xxzz_xxzz_1,   \
                             ta_xxzz_xxzzz_0,  \
                             ta_xxzz_xxzzz_1,  \
                             ta_xxzz_xyyy_0,   \
                             ta_xxzz_xyyy_1,   \
                             ta_xxzz_xyyyy_0,  \
                             ta_xxzz_xyyyy_1,  \
                             ta_xxzz_xyyyz_0,  \
                             ta_xxzz_xyyyz_1,  \
                             ta_xxzz_xyyz_0,   \
                             ta_xxzz_xyyz_1,   \
                             ta_xxzz_xyyzz_0,  \
                             ta_xxzz_xyyzz_1,  \
                             ta_xxzz_xyzz_0,   \
                             ta_xxzz_xyzz_1,   \
                             ta_xxzz_xyzzz_0,  \
                             ta_xxzz_xyzzz_1,  \
                             ta_xxzz_xzzz_0,   \
                             ta_xxzz_xzzz_1,   \
                             ta_xxzz_xzzzz_0,  \
                             ta_xxzz_xzzzz_1,  \
                             ta_xxzz_zzzzz_0,  \
                             ta_xxzz_zzzzz_1,  \
                             ta_xyzz_yyyyy_0,  \
                             ta_xyzz_yyyyy_1,  \
                             ta_xyzz_yyyyz_0,  \
                             ta_xyzz_yyyyz_1,  \
                             ta_xyzz_yyyzz_0,  \
                             ta_xyzz_yyyzz_1,  \
                             ta_xyzz_yyzzz_0,  \
                             ta_xyzz_yyzzz_1,  \
                             ta_xyzz_yzzzz_0,  \
                             ta_xyzz_yzzzz_1,  \
                             ta_yzz_yyyyy_0,   \
                             ta_yzz_yyyyy_1,   \
                             ta_yzz_yyyyz_0,   \
                             ta_yzz_yyyyz_1,   \
                             ta_yzz_yyyzz_0,   \
                             ta_yzz_yyyzz_1,   \
                             ta_yzz_yyzzz_0,   \
                             ta_yzz_yyzzz_1,   \
                             ta_yzz_yzzzz_0,   \
                             ta_yzz_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzz_xxxxx_0[i] = ta_xxzz_xxxxx_0[i] * pa_y[i] - ta_xxzz_xxxxx_1[i] * pc_y[i];

        ta_xxyzz_xxxxy_0[i] = ta_xxzz_xxxx_0[i] * fe_0 - ta_xxzz_xxxx_1[i] * fe_0 + ta_xxzz_xxxxy_0[i] * pa_y[i] - ta_xxzz_xxxxy_1[i] * pc_y[i];

        ta_xxyzz_xxxxz_0[i] = ta_xxzz_xxxxz_0[i] * pa_y[i] - ta_xxzz_xxxxz_1[i] * pc_y[i];

        ta_xxyzz_xxxyy_0[i] =
            2.0 * ta_xxzz_xxxy_0[i] * fe_0 - 2.0 * ta_xxzz_xxxy_1[i] * fe_0 + ta_xxzz_xxxyy_0[i] * pa_y[i] - ta_xxzz_xxxyy_1[i] * pc_y[i];

        ta_xxyzz_xxxyz_0[i] = ta_xxzz_xxxz_0[i] * fe_0 - ta_xxzz_xxxz_1[i] * fe_0 + ta_xxzz_xxxyz_0[i] * pa_y[i] - ta_xxzz_xxxyz_1[i] * pc_y[i];

        ta_xxyzz_xxxzz_0[i] = ta_xxzz_xxxzz_0[i] * pa_y[i] - ta_xxzz_xxxzz_1[i] * pc_y[i];

        ta_xxyzz_xxyyy_0[i] =
            3.0 * ta_xxzz_xxyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxyy_1[i] * fe_0 + ta_xxzz_xxyyy_0[i] * pa_y[i] - ta_xxzz_xxyyy_1[i] * pc_y[i];

        ta_xxyzz_xxyyz_0[i] =
            2.0 * ta_xxzz_xxyz_0[i] * fe_0 - 2.0 * ta_xxzz_xxyz_1[i] * fe_0 + ta_xxzz_xxyyz_0[i] * pa_y[i] - ta_xxzz_xxyyz_1[i] * pc_y[i];

        ta_xxyzz_xxyzz_0[i] = ta_xxzz_xxzz_0[i] * fe_0 - ta_xxzz_xxzz_1[i] * fe_0 + ta_xxzz_xxyzz_0[i] * pa_y[i] - ta_xxzz_xxyzz_1[i] * pc_y[i];

        ta_xxyzz_xxzzz_0[i] = ta_xxzz_xxzzz_0[i] * pa_y[i] - ta_xxzz_xxzzz_1[i] * pc_y[i];

        ta_xxyzz_xyyyy_0[i] =
            4.0 * ta_xxzz_xyyy_0[i] * fe_0 - 4.0 * ta_xxzz_xyyy_1[i] * fe_0 + ta_xxzz_xyyyy_0[i] * pa_y[i] - ta_xxzz_xyyyy_1[i] * pc_y[i];

        ta_xxyzz_xyyyz_0[i] =
            3.0 * ta_xxzz_xyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyz_1[i] * fe_0 + ta_xxzz_xyyyz_0[i] * pa_y[i] - ta_xxzz_xyyyz_1[i] * pc_y[i];

        ta_xxyzz_xyyzz_0[i] =
            2.0 * ta_xxzz_xyzz_0[i] * fe_0 - 2.0 * ta_xxzz_xyzz_1[i] * fe_0 + ta_xxzz_xyyzz_0[i] * pa_y[i] - ta_xxzz_xyyzz_1[i] * pc_y[i];

        ta_xxyzz_xyzzz_0[i] = ta_xxzz_xzzz_0[i] * fe_0 - ta_xxzz_xzzz_1[i] * fe_0 + ta_xxzz_xyzzz_0[i] * pa_y[i] - ta_xxzz_xyzzz_1[i] * pc_y[i];

        ta_xxyzz_xzzzz_0[i] = ta_xxzz_xzzzz_0[i] * pa_y[i] - ta_xxzz_xzzzz_1[i] * pc_y[i];

        ta_xxyzz_yyyyy_0[i] = ta_yzz_yyyyy_0[i] * fe_0 - ta_yzz_yyyyy_1[i] * fe_0 + ta_xyzz_yyyyy_0[i] * pa_x[i] - ta_xyzz_yyyyy_1[i] * pc_x[i];

        ta_xxyzz_yyyyz_0[i] = ta_yzz_yyyyz_0[i] * fe_0 - ta_yzz_yyyyz_1[i] * fe_0 + ta_xyzz_yyyyz_0[i] * pa_x[i] - ta_xyzz_yyyyz_1[i] * pc_x[i];

        ta_xxyzz_yyyzz_0[i] = ta_yzz_yyyzz_0[i] * fe_0 - ta_yzz_yyyzz_1[i] * fe_0 + ta_xyzz_yyyzz_0[i] * pa_x[i] - ta_xyzz_yyyzz_1[i] * pc_x[i];

        ta_xxyzz_yyzzz_0[i] = ta_yzz_yyzzz_0[i] * fe_0 - ta_yzz_yyzzz_1[i] * fe_0 + ta_xyzz_yyzzz_0[i] * pa_x[i] - ta_xyzz_yyzzz_1[i] * pc_x[i];

        ta_xxyzz_yzzzz_0[i] = ta_yzz_yzzzz_0[i] * fe_0 - ta_yzz_yzzzz_1[i] * fe_0 + ta_xyzz_yzzzz_0[i] * pa_x[i] - ta_xyzz_yzzzz_1[i] * pc_x[i];

        ta_xxyzz_zzzzz_0[i] = ta_xxzz_zzzzz_0[i] * pa_y[i] - ta_xxzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 189-210 components of targeted buffer : HH

    auto ta_xxzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 189);

    auto ta_xxzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 190);

    auto ta_xxzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 191);

    auto ta_xxzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 192);

    auto ta_xxzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 193);

    auto ta_xxzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 194);

    auto ta_xxzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 195);

    auto ta_xxzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 196);

    auto ta_xxzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 197);

    auto ta_xxzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 198);

    auto ta_xxzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 199);

    auto ta_xxzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 200);

    auto ta_xxzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 201);

    auto ta_xxzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 202);

    auto ta_xxzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 203);

    auto ta_xxzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 204);

    auto ta_xxzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 205);

    auto ta_xxzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 206);

    auto ta_xxzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 207);

    auto ta_xxzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 208);

    auto ta_xxzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 209);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xxz_xxxxx_0,   \
                             ta_xxz_xxxxx_1,   \
                             ta_xxz_xxxxy_0,   \
                             ta_xxz_xxxxy_1,   \
                             ta_xxz_xxxyy_0,   \
                             ta_xxz_xxxyy_1,   \
                             ta_xxz_xxyyy_0,   \
                             ta_xxz_xxyyy_1,   \
                             ta_xxz_xyyyy_0,   \
                             ta_xxz_xyyyy_1,   \
                             ta_xxzz_xxxxx_0,  \
                             ta_xxzz_xxxxx_1,  \
                             ta_xxzz_xxxxy_0,  \
                             ta_xxzz_xxxxy_1,  \
                             ta_xxzz_xxxyy_0,  \
                             ta_xxzz_xxxyy_1,  \
                             ta_xxzz_xxyyy_0,  \
                             ta_xxzz_xxyyy_1,  \
                             ta_xxzz_xyyyy_0,  \
                             ta_xxzz_xyyyy_1,  \
                             ta_xxzzz_xxxxx_0, \
                             ta_xxzzz_xxxxy_0, \
                             ta_xxzzz_xxxxz_0, \
                             ta_xxzzz_xxxyy_0, \
                             ta_xxzzz_xxxyz_0, \
                             ta_xxzzz_xxxzz_0, \
                             ta_xxzzz_xxyyy_0, \
                             ta_xxzzz_xxyyz_0, \
                             ta_xxzzz_xxyzz_0, \
                             ta_xxzzz_xxzzz_0, \
                             ta_xxzzz_xyyyy_0, \
                             ta_xxzzz_xyyyz_0, \
                             ta_xxzzz_xyyzz_0, \
                             ta_xxzzz_xyzzz_0, \
                             ta_xxzzz_xzzzz_0, \
                             ta_xxzzz_yyyyy_0, \
                             ta_xxzzz_yyyyz_0, \
                             ta_xxzzz_yyyzz_0, \
                             ta_xxzzz_yyzzz_0, \
                             ta_xxzzz_yzzzz_0, \
                             ta_xxzzz_zzzzz_0, \
                             ta_xzzz_xxxxz_0,  \
                             ta_xzzz_xxxxz_1,  \
                             ta_xzzz_xxxyz_0,  \
                             ta_xzzz_xxxyz_1,  \
                             ta_xzzz_xxxz_0,   \
                             ta_xzzz_xxxz_1,   \
                             ta_xzzz_xxxzz_0,  \
                             ta_xzzz_xxxzz_1,  \
                             ta_xzzz_xxyyz_0,  \
                             ta_xzzz_xxyyz_1,  \
                             ta_xzzz_xxyz_0,   \
                             ta_xzzz_xxyz_1,   \
                             ta_xzzz_xxyzz_0,  \
                             ta_xzzz_xxyzz_1,  \
                             ta_xzzz_xxzz_0,   \
                             ta_xzzz_xxzz_1,   \
                             ta_xzzz_xxzzz_0,  \
                             ta_xzzz_xxzzz_1,  \
                             ta_xzzz_xyyyz_0,  \
                             ta_xzzz_xyyyz_1,  \
                             ta_xzzz_xyyz_0,   \
                             ta_xzzz_xyyz_1,   \
                             ta_xzzz_xyyzz_0,  \
                             ta_xzzz_xyyzz_1,  \
                             ta_xzzz_xyzz_0,   \
                             ta_xzzz_xyzz_1,   \
                             ta_xzzz_xyzzz_0,  \
                             ta_xzzz_xyzzz_1,  \
                             ta_xzzz_xzzz_0,   \
                             ta_xzzz_xzzz_1,   \
                             ta_xzzz_xzzzz_0,  \
                             ta_xzzz_xzzzz_1,  \
                             ta_xzzz_yyyyy_0,  \
                             ta_xzzz_yyyyy_1,  \
                             ta_xzzz_yyyyz_0,  \
                             ta_xzzz_yyyyz_1,  \
                             ta_xzzz_yyyz_0,   \
                             ta_xzzz_yyyz_1,   \
                             ta_xzzz_yyyzz_0,  \
                             ta_xzzz_yyyzz_1,  \
                             ta_xzzz_yyzz_0,   \
                             ta_xzzz_yyzz_1,   \
                             ta_xzzz_yyzzz_0,  \
                             ta_xzzz_yyzzz_1,  \
                             ta_xzzz_yzzz_0,   \
                             ta_xzzz_yzzz_1,   \
                             ta_xzzz_yzzzz_0,  \
                             ta_xzzz_yzzzz_1,  \
                             ta_xzzz_zzzz_0,   \
                             ta_xzzz_zzzz_1,   \
                             ta_xzzz_zzzzz_0,  \
                             ta_xzzz_zzzzz_1,  \
                             ta_zzz_xxxxz_0,   \
                             ta_zzz_xxxxz_1,   \
                             ta_zzz_xxxyz_0,   \
                             ta_zzz_xxxyz_1,   \
                             ta_zzz_xxxzz_0,   \
                             ta_zzz_xxxzz_1,   \
                             ta_zzz_xxyyz_0,   \
                             ta_zzz_xxyyz_1,   \
                             ta_zzz_xxyzz_0,   \
                             ta_zzz_xxyzz_1,   \
                             ta_zzz_xxzzz_0,   \
                             ta_zzz_xxzzz_1,   \
                             ta_zzz_xyyyz_0,   \
                             ta_zzz_xyyyz_1,   \
                             ta_zzz_xyyzz_0,   \
                             ta_zzz_xyyzz_1,   \
                             ta_zzz_xyzzz_0,   \
                             ta_zzz_xyzzz_1,   \
                             ta_zzz_xzzzz_0,   \
                             ta_zzz_xzzzz_1,   \
                             ta_zzz_yyyyy_0,   \
                             ta_zzz_yyyyy_1,   \
                             ta_zzz_yyyyz_0,   \
                             ta_zzz_yyyyz_1,   \
                             ta_zzz_yyyzz_0,   \
                             ta_zzz_yyyzz_1,   \
                             ta_zzz_yyzzz_0,   \
                             ta_zzz_yyzzz_1,   \
                             ta_zzz_yzzzz_0,   \
                             ta_zzz_yzzzz_1,   \
                             ta_zzz_zzzzz_0,   \
                             ta_zzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzz_xxxxx_0[i] =
            2.0 * ta_xxz_xxxxx_0[i] * fe_0 - 2.0 * ta_xxz_xxxxx_1[i] * fe_0 + ta_xxzz_xxxxx_0[i] * pa_z[i] - ta_xxzz_xxxxx_1[i] * pc_z[i];

        ta_xxzzz_xxxxy_0[i] =
            2.0 * ta_xxz_xxxxy_0[i] * fe_0 - 2.0 * ta_xxz_xxxxy_1[i] * fe_0 + ta_xxzz_xxxxy_0[i] * pa_z[i] - ta_xxzz_xxxxy_1[i] * pc_z[i];

        ta_xxzzz_xxxxz_0[i] = ta_zzz_xxxxz_0[i] * fe_0 - ta_zzz_xxxxz_1[i] * fe_0 + 4.0 * ta_xzzz_xxxz_0[i] * fe_0 - 4.0 * ta_xzzz_xxxz_1[i] * fe_0 +
                              ta_xzzz_xxxxz_0[i] * pa_x[i] - ta_xzzz_xxxxz_1[i] * pc_x[i];

        ta_xxzzz_xxxyy_0[i] =
            2.0 * ta_xxz_xxxyy_0[i] * fe_0 - 2.0 * ta_xxz_xxxyy_1[i] * fe_0 + ta_xxzz_xxxyy_0[i] * pa_z[i] - ta_xxzz_xxxyy_1[i] * pc_z[i];

        ta_xxzzz_xxxyz_0[i] = ta_zzz_xxxyz_0[i] * fe_0 - ta_zzz_xxxyz_1[i] * fe_0 + 3.0 * ta_xzzz_xxyz_0[i] * fe_0 - 3.0 * ta_xzzz_xxyz_1[i] * fe_0 +
                              ta_xzzz_xxxyz_0[i] * pa_x[i] - ta_xzzz_xxxyz_1[i] * pc_x[i];

        ta_xxzzz_xxxzz_0[i] = ta_zzz_xxxzz_0[i] * fe_0 - ta_zzz_xxxzz_1[i] * fe_0 + 3.0 * ta_xzzz_xxzz_0[i] * fe_0 - 3.0 * ta_xzzz_xxzz_1[i] * fe_0 +
                              ta_xzzz_xxxzz_0[i] * pa_x[i] - ta_xzzz_xxxzz_1[i] * pc_x[i];

        ta_xxzzz_xxyyy_0[i] =
            2.0 * ta_xxz_xxyyy_0[i] * fe_0 - 2.0 * ta_xxz_xxyyy_1[i] * fe_0 + ta_xxzz_xxyyy_0[i] * pa_z[i] - ta_xxzz_xxyyy_1[i] * pc_z[i];

        ta_xxzzz_xxyyz_0[i] = ta_zzz_xxyyz_0[i] * fe_0 - ta_zzz_xxyyz_1[i] * fe_0 + 2.0 * ta_xzzz_xyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyz_1[i] * fe_0 +
                              ta_xzzz_xxyyz_0[i] * pa_x[i] - ta_xzzz_xxyyz_1[i] * pc_x[i];

        ta_xxzzz_xxyzz_0[i] = ta_zzz_xxyzz_0[i] * fe_0 - ta_zzz_xxyzz_1[i] * fe_0 + 2.0 * ta_xzzz_xyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyzz_1[i] * fe_0 +
                              ta_xzzz_xxyzz_0[i] * pa_x[i] - ta_xzzz_xxyzz_1[i] * pc_x[i];

        ta_xxzzz_xxzzz_0[i] = ta_zzz_xxzzz_0[i] * fe_0 - ta_zzz_xxzzz_1[i] * fe_0 + 2.0 * ta_xzzz_xzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xzzz_1[i] * fe_0 +
                              ta_xzzz_xxzzz_0[i] * pa_x[i] - ta_xzzz_xxzzz_1[i] * pc_x[i];

        ta_xxzzz_xyyyy_0[i] =
            2.0 * ta_xxz_xyyyy_0[i] * fe_0 - 2.0 * ta_xxz_xyyyy_1[i] * fe_0 + ta_xxzz_xyyyy_0[i] * pa_z[i] - ta_xxzz_xyyyy_1[i] * pc_z[i];

        ta_xxzzz_xyyyz_0[i] = ta_zzz_xyyyz_0[i] * fe_0 - ta_zzz_xyyyz_1[i] * fe_0 + ta_xzzz_yyyz_0[i] * fe_0 - ta_xzzz_yyyz_1[i] * fe_0 +
                              ta_xzzz_xyyyz_0[i] * pa_x[i] - ta_xzzz_xyyyz_1[i] * pc_x[i];

        ta_xxzzz_xyyzz_0[i] = ta_zzz_xyyzz_0[i] * fe_0 - ta_zzz_xyyzz_1[i] * fe_0 + ta_xzzz_yyzz_0[i] * fe_0 - ta_xzzz_yyzz_1[i] * fe_0 +
                              ta_xzzz_xyyzz_0[i] * pa_x[i] - ta_xzzz_xyyzz_1[i] * pc_x[i];

        ta_xxzzz_xyzzz_0[i] = ta_zzz_xyzzz_0[i] * fe_0 - ta_zzz_xyzzz_1[i] * fe_0 + ta_xzzz_yzzz_0[i] * fe_0 - ta_xzzz_yzzz_1[i] * fe_0 +
                              ta_xzzz_xyzzz_0[i] * pa_x[i] - ta_xzzz_xyzzz_1[i] * pc_x[i];

        ta_xxzzz_xzzzz_0[i] = ta_zzz_xzzzz_0[i] * fe_0 - ta_zzz_xzzzz_1[i] * fe_0 + ta_xzzz_zzzz_0[i] * fe_0 - ta_xzzz_zzzz_1[i] * fe_0 +
                              ta_xzzz_xzzzz_0[i] * pa_x[i] - ta_xzzz_xzzzz_1[i] * pc_x[i];

        ta_xxzzz_yyyyy_0[i] = ta_zzz_yyyyy_0[i] * fe_0 - ta_zzz_yyyyy_1[i] * fe_0 + ta_xzzz_yyyyy_0[i] * pa_x[i] - ta_xzzz_yyyyy_1[i] * pc_x[i];

        ta_xxzzz_yyyyz_0[i] = ta_zzz_yyyyz_0[i] * fe_0 - ta_zzz_yyyyz_1[i] * fe_0 + ta_xzzz_yyyyz_0[i] * pa_x[i] - ta_xzzz_yyyyz_1[i] * pc_x[i];

        ta_xxzzz_yyyzz_0[i] = ta_zzz_yyyzz_0[i] * fe_0 - ta_zzz_yyyzz_1[i] * fe_0 + ta_xzzz_yyyzz_0[i] * pa_x[i] - ta_xzzz_yyyzz_1[i] * pc_x[i];

        ta_xxzzz_yyzzz_0[i] = ta_zzz_yyzzz_0[i] * fe_0 - ta_zzz_yyzzz_1[i] * fe_0 + ta_xzzz_yyzzz_0[i] * pa_x[i] - ta_xzzz_yyzzz_1[i] * pc_x[i];

        ta_xxzzz_yzzzz_0[i] = ta_zzz_yzzzz_0[i] * fe_0 - ta_zzz_yzzzz_1[i] * fe_0 + ta_xzzz_yzzzz_0[i] * pa_x[i] - ta_xzzz_yzzzz_1[i] * pc_x[i];

        ta_xxzzz_zzzzz_0[i] = ta_zzz_zzzzz_0[i] * fe_0 - ta_zzz_zzzzz_1[i] * fe_0 + ta_xzzz_zzzzz_0[i] * pa_x[i] - ta_xzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 210-231 components of targeted buffer : HH

    auto ta_xyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 210);

    auto ta_xyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 211);

    auto ta_xyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 212);

    auto ta_xyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 213);

    auto ta_xyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 214);

    auto ta_xyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 215);

    auto ta_xyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 216);

    auto ta_xyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 217);

    auto ta_xyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 218);

    auto ta_xyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 219);

    auto ta_xyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 220);

    auto ta_xyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 221);

    auto ta_xyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 222);

    auto ta_xyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 223);

    auto ta_xyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 224);

    auto ta_xyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 225);

    auto ta_xyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 226);

    auto ta_xyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 227);

    auto ta_xyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 228);

    auto ta_xyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 229);

    auto ta_xyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 230);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xyyyy_xxxxx_0, \
                             ta_xyyyy_xxxxy_0, \
                             ta_xyyyy_xxxxz_0, \
                             ta_xyyyy_xxxyy_0, \
                             ta_xyyyy_xxxyz_0, \
                             ta_xyyyy_xxxzz_0, \
                             ta_xyyyy_xxyyy_0, \
                             ta_xyyyy_xxyyz_0, \
                             ta_xyyyy_xxyzz_0, \
                             ta_xyyyy_xxzzz_0, \
                             ta_xyyyy_xyyyy_0, \
                             ta_xyyyy_xyyyz_0, \
                             ta_xyyyy_xyyzz_0, \
                             ta_xyyyy_xyzzz_0, \
                             ta_xyyyy_xzzzz_0, \
                             ta_xyyyy_yyyyy_0, \
                             ta_xyyyy_yyyyz_0, \
                             ta_xyyyy_yyyzz_0, \
                             ta_xyyyy_yyzzz_0, \
                             ta_xyyyy_yzzzz_0, \
                             ta_xyyyy_zzzzz_0, \
                             ta_yyyy_xxxx_0,   \
                             ta_yyyy_xxxx_1,   \
                             ta_yyyy_xxxxx_0,  \
                             ta_yyyy_xxxxx_1,  \
                             ta_yyyy_xxxxy_0,  \
                             ta_yyyy_xxxxy_1,  \
                             ta_yyyy_xxxxz_0,  \
                             ta_yyyy_xxxxz_1,  \
                             ta_yyyy_xxxy_0,   \
                             ta_yyyy_xxxy_1,   \
                             ta_yyyy_xxxyy_0,  \
                             ta_yyyy_xxxyy_1,  \
                             ta_yyyy_xxxyz_0,  \
                             ta_yyyy_xxxyz_1,  \
                             ta_yyyy_xxxz_0,   \
                             ta_yyyy_xxxz_1,   \
                             ta_yyyy_xxxzz_0,  \
                             ta_yyyy_xxxzz_1,  \
                             ta_yyyy_xxyy_0,   \
                             ta_yyyy_xxyy_1,   \
                             ta_yyyy_xxyyy_0,  \
                             ta_yyyy_xxyyy_1,  \
                             ta_yyyy_xxyyz_0,  \
                             ta_yyyy_xxyyz_1,  \
                             ta_yyyy_xxyz_0,   \
                             ta_yyyy_xxyz_1,   \
                             ta_yyyy_xxyzz_0,  \
                             ta_yyyy_xxyzz_1,  \
                             ta_yyyy_xxzz_0,   \
                             ta_yyyy_xxzz_1,   \
                             ta_yyyy_xxzzz_0,  \
                             ta_yyyy_xxzzz_1,  \
                             ta_yyyy_xyyy_0,   \
                             ta_yyyy_xyyy_1,   \
                             ta_yyyy_xyyyy_0,  \
                             ta_yyyy_xyyyy_1,  \
                             ta_yyyy_xyyyz_0,  \
                             ta_yyyy_xyyyz_1,  \
                             ta_yyyy_xyyz_0,   \
                             ta_yyyy_xyyz_1,   \
                             ta_yyyy_xyyzz_0,  \
                             ta_yyyy_xyyzz_1,  \
                             ta_yyyy_xyzz_0,   \
                             ta_yyyy_xyzz_1,   \
                             ta_yyyy_xyzzz_0,  \
                             ta_yyyy_xyzzz_1,  \
                             ta_yyyy_xzzz_0,   \
                             ta_yyyy_xzzz_1,   \
                             ta_yyyy_xzzzz_0,  \
                             ta_yyyy_xzzzz_1,  \
                             ta_yyyy_yyyy_0,   \
                             ta_yyyy_yyyy_1,   \
                             ta_yyyy_yyyyy_0,  \
                             ta_yyyy_yyyyy_1,  \
                             ta_yyyy_yyyyz_0,  \
                             ta_yyyy_yyyyz_1,  \
                             ta_yyyy_yyyz_0,   \
                             ta_yyyy_yyyz_1,   \
                             ta_yyyy_yyyzz_0,  \
                             ta_yyyy_yyyzz_1,  \
                             ta_yyyy_yyzz_0,   \
                             ta_yyyy_yyzz_1,   \
                             ta_yyyy_yyzzz_0,  \
                             ta_yyyy_yyzzz_1,  \
                             ta_yyyy_yzzz_0,   \
                             ta_yyyy_yzzz_1,   \
                             ta_yyyy_yzzzz_0,  \
                             ta_yyyy_yzzzz_1,  \
                             ta_yyyy_zzzz_0,   \
                             ta_yyyy_zzzz_1,   \
                             ta_yyyy_zzzzz_0,  \
                             ta_yyyy_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyy_xxxxx_0[i] =
            5.0 * ta_yyyy_xxxx_0[i] * fe_0 - 5.0 * ta_yyyy_xxxx_1[i] * fe_0 + ta_yyyy_xxxxx_0[i] * pa_x[i] - ta_yyyy_xxxxx_1[i] * pc_x[i];

        ta_xyyyy_xxxxy_0[i] =
            4.0 * ta_yyyy_xxxy_0[i] * fe_0 - 4.0 * ta_yyyy_xxxy_1[i] * fe_0 + ta_yyyy_xxxxy_0[i] * pa_x[i] - ta_yyyy_xxxxy_1[i] * pc_x[i];

        ta_xyyyy_xxxxz_0[i] =
            4.0 * ta_yyyy_xxxz_0[i] * fe_0 - 4.0 * ta_yyyy_xxxz_1[i] * fe_0 + ta_yyyy_xxxxz_0[i] * pa_x[i] - ta_yyyy_xxxxz_1[i] * pc_x[i];

        ta_xyyyy_xxxyy_0[i] =
            3.0 * ta_yyyy_xxyy_0[i] * fe_0 - 3.0 * ta_yyyy_xxyy_1[i] * fe_0 + ta_yyyy_xxxyy_0[i] * pa_x[i] - ta_yyyy_xxxyy_1[i] * pc_x[i];

        ta_xyyyy_xxxyz_0[i] =
            3.0 * ta_yyyy_xxyz_0[i] * fe_0 - 3.0 * ta_yyyy_xxyz_1[i] * fe_0 + ta_yyyy_xxxyz_0[i] * pa_x[i] - ta_yyyy_xxxyz_1[i] * pc_x[i];

        ta_xyyyy_xxxzz_0[i] =
            3.0 * ta_yyyy_xxzz_0[i] * fe_0 - 3.0 * ta_yyyy_xxzz_1[i] * fe_0 + ta_yyyy_xxxzz_0[i] * pa_x[i] - ta_yyyy_xxxzz_1[i] * pc_x[i];

        ta_xyyyy_xxyyy_0[i] =
            2.0 * ta_yyyy_xyyy_0[i] * fe_0 - 2.0 * ta_yyyy_xyyy_1[i] * fe_0 + ta_yyyy_xxyyy_0[i] * pa_x[i] - ta_yyyy_xxyyy_1[i] * pc_x[i];

        ta_xyyyy_xxyyz_0[i] =
            2.0 * ta_yyyy_xyyz_0[i] * fe_0 - 2.0 * ta_yyyy_xyyz_1[i] * fe_0 + ta_yyyy_xxyyz_0[i] * pa_x[i] - ta_yyyy_xxyyz_1[i] * pc_x[i];

        ta_xyyyy_xxyzz_0[i] =
            2.0 * ta_yyyy_xyzz_0[i] * fe_0 - 2.0 * ta_yyyy_xyzz_1[i] * fe_0 + ta_yyyy_xxyzz_0[i] * pa_x[i] - ta_yyyy_xxyzz_1[i] * pc_x[i];

        ta_xyyyy_xxzzz_0[i] =
            2.0 * ta_yyyy_xzzz_0[i] * fe_0 - 2.0 * ta_yyyy_xzzz_1[i] * fe_0 + ta_yyyy_xxzzz_0[i] * pa_x[i] - ta_yyyy_xxzzz_1[i] * pc_x[i];

        ta_xyyyy_xyyyy_0[i] = ta_yyyy_yyyy_0[i] * fe_0 - ta_yyyy_yyyy_1[i] * fe_0 + ta_yyyy_xyyyy_0[i] * pa_x[i] - ta_yyyy_xyyyy_1[i] * pc_x[i];

        ta_xyyyy_xyyyz_0[i] = ta_yyyy_yyyz_0[i] * fe_0 - ta_yyyy_yyyz_1[i] * fe_0 + ta_yyyy_xyyyz_0[i] * pa_x[i] - ta_yyyy_xyyyz_1[i] * pc_x[i];

        ta_xyyyy_xyyzz_0[i] = ta_yyyy_yyzz_0[i] * fe_0 - ta_yyyy_yyzz_1[i] * fe_0 + ta_yyyy_xyyzz_0[i] * pa_x[i] - ta_yyyy_xyyzz_1[i] * pc_x[i];

        ta_xyyyy_xyzzz_0[i] = ta_yyyy_yzzz_0[i] * fe_0 - ta_yyyy_yzzz_1[i] * fe_0 + ta_yyyy_xyzzz_0[i] * pa_x[i] - ta_yyyy_xyzzz_1[i] * pc_x[i];

        ta_xyyyy_xzzzz_0[i] = ta_yyyy_zzzz_0[i] * fe_0 - ta_yyyy_zzzz_1[i] * fe_0 + ta_yyyy_xzzzz_0[i] * pa_x[i] - ta_yyyy_xzzzz_1[i] * pc_x[i];

        ta_xyyyy_yyyyy_0[i] = ta_yyyy_yyyyy_0[i] * pa_x[i] - ta_yyyy_yyyyy_1[i] * pc_x[i];

        ta_xyyyy_yyyyz_0[i] = ta_yyyy_yyyyz_0[i] * pa_x[i] - ta_yyyy_yyyyz_1[i] * pc_x[i];

        ta_xyyyy_yyyzz_0[i] = ta_yyyy_yyyzz_0[i] * pa_x[i] - ta_yyyy_yyyzz_1[i] * pc_x[i];

        ta_xyyyy_yyzzz_0[i] = ta_yyyy_yyzzz_0[i] * pa_x[i] - ta_yyyy_yyzzz_1[i] * pc_x[i];

        ta_xyyyy_yzzzz_0[i] = ta_yyyy_yzzzz_0[i] * pa_x[i] - ta_yyyy_yzzzz_1[i] * pc_x[i];

        ta_xyyyy_zzzzz_0[i] = ta_yyyy_zzzzz_0[i] * pa_x[i] - ta_yyyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 231-252 components of targeted buffer : HH

    auto ta_xyyyz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 231);

    auto ta_xyyyz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 232);

    auto ta_xyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 233);

    auto ta_xyyyz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 234);

    auto ta_xyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 235);

    auto ta_xyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 236);

    auto ta_xyyyz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 237);

    auto ta_xyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 238);

    auto ta_xyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 239);

    auto ta_xyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 240);

    auto ta_xyyyz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 241);

    auto ta_xyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 242);

    auto ta_xyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 243);

    auto ta_xyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 244);

    auto ta_xyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 245);

    auto ta_xyyyz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 246);

    auto ta_xyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 247);

    auto ta_xyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 248);

    auto ta_xyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 249);

    auto ta_xyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 250);

    auto ta_xyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 251);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta_xyyy_xxxxx_0,  \
                             ta_xyyy_xxxxx_1,  \
                             ta_xyyy_xxxxy_0,  \
                             ta_xyyy_xxxxy_1,  \
                             ta_xyyy_xxxyy_0,  \
                             ta_xyyy_xxxyy_1,  \
                             ta_xyyy_xxyyy_0,  \
                             ta_xyyy_xxyyy_1,  \
                             ta_xyyy_xyyyy_0,  \
                             ta_xyyy_xyyyy_1,  \
                             ta_xyyyz_xxxxx_0, \
                             ta_xyyyz_xxxxy_0, \
                             ta_xyyyz_xxxxz_0, \
                             ta_xyyyz_xxxyy_0, \
                             ta_xyyyz_xxxyz_0, \
                             ta_xyyyz_xxxzz_0, \
                             ta_xyyyz_xxyyy_0, \
                             ta_xyyyz_xxyyz_0, \
                             ta_xyyyz_xxyzz_0, \
                             ta_xyyyz_xxzzz_0, \
                             ta_xyyyz_xyyyy_0, \
                             ta_xyyyz_xyyyz_0, \
                             ta_xyyyz_xyyzz_0, \
                             ta_xyyyz_xyzzz_0, \
                             ta_xyyyz_xzzzz_0, \
                             ta_xyyyz_yyyyy_0, \
                             ta_xyyyz_yyyyz_0, \
                             ta_xyyyz_yyyzz_0, \
                             ta_xyyyz_yyzzz_0, \
                             ta_xyyyz_yzzzz_0, \
                             ta_xyyyz_zzzzz_0, \
                             ta_yyyz_xxxxz_0,  \
                             ta_yyyz_xxxxz_1,  \
                             ta_yyyz_xxxyz_0,  \
                             ta_yyyz_xxxyz_1,  \
                             ta_yyyz_xxxz_0,   \
                             ta_yyyz_xxxz_1,   \
                             ta_yyyz_xxxzz_0,  \
                             ta_yyyz_xxxzz_1,  \
                             ta_yyyz_xxyyz_0,  \
                             ta_yyyz_xxyyz_1,  \
                             ta_yyyz_xxyz_0,   \
                             ta_yyyz_xxyz_1,   \
                             ta_yyyz_xxyzz_0,  \
                             ta_yyyz_xxyzz_1,  \
                             ta_yyyz_xxzz_0,   \
                             ta_yyyz_xxzz_1,   \
                             ta_yyyz_xxzzz_0,  \
                             ta_yyyz_xxzzz_1,  \
                             ta_yyyz_xyyyz_0,  \
                             ta_yyyz_xyyyz_1,  \
                             ta_yyyz_xyyz_0,   \
                             ta_yyyz_xyyz_1,   \
                             ta_yyyz_xyyzz_0,  \
                             ta_yyyz_xyyzz_1,  \
                             ta_yyyz_xyzz_0,   \
                             ta_yyyz_xyzz_1,   \
                             ta_yyyz_xyzzz_0,  \
                             ta_yyyz_xyzzz_1,  \
                             ta_yyyz_xzzz_0,   \
                             ta_yyyz_xzzz_1,   \
                             ta_yyyz_xzzzz_0,  \
                             ta_yyyz_xzzzz_1,  \
                             ta_yyyz_yyyyy_0,  \
                             ta_yyyz_yyyyy_1,  \
                             ta_yyyz_yyyyz_0,  \
                             ta_yyyz_yyyyz_1,  \
                             ta_yyyz_yyyz_0,   \
                             ta_yyyz_yyyz_1,   \
                             ta_yyyz_yyyzz_0,  \
                             ta_yyyz_yyyzz_1,  \
                             ta_yyyz_yyzz_0,   \
                             ta_yyyz_yyzz_1,   \
                             ta_yyyz_yyzzz_0,  \
                             ta_yyyz_yyzzz_1,  \
                             ta_yyyz_yzzz_0,   \
                             ta_yyyz_yzzz_1,   \
                             ta_yyyz_yzzzz_0,  \
                             ta_yyyz_yzzzz_1,  \
                             ta_yyyz_zzzz_0,   \
                             ta_yyyz_zzzz_1,   \
                             ta_yyyz_zzzzz_0,  \
                             ta_yyyz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyz_xxxxx_0[i] = ta_xyyy_xxxxx_0[i] * pa_z[i] - ta_xyyy_xxxxx_1[i] * pc_z[i];

        ta_xyyyz_xxxxy_0[i] = ta_xyyy_xxxxy_0[i] * pa_z[i] - ta_xyyy_xxxxy_1[i] * pc_z[i];

        ta_xyyyz_xxxxz_0[i] =
            4.0 * ta_yyyz_xxxz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxz_1[i] * fe_0 + ta_yyyz_xxxxz_0[i] * pa_x[i] - ta_yyyz_xxxxz_1[i] * pc_x[i];

        ta_xyyyz_xxxyy_0[i] = ta_xyyy_xxxyy_0[i] * pa_z[i] - ta_xyyy_xxxyy_1[i] * pc_z[i];

        ta_xyyyz_xxxyz_0[i] =
            3.0 * ta_yyyz_xxyz_0[i] * fe_0 - 3.0 * ta_yyyz_xxyz_1[i] * fe_0 + ta_yyyz_xxxyz_0[i] * pa_x[i] - ta_yyyz_xxxyz_1[i] * pc_x[i];

        ta_xyyyz_xxxzz_0[i] =
            3.0 * ta_yyyz_xxzz_0[i] * fe_0 - 3.0 * ta_yyyz_xxzz_1[i] * fe_0 + ta_yyyz_xxxzz_0[i] * pa_x[i] - ta_yyyz_xxxzz_1[i] * pc_x[i];

        ta_xyyyz_xxyyy_0[i] = ta_xyyy_xxyyy_0[i] * pa_z[i] - ta_xyyy_xxyyy_1[i] * pc_z[i];

        ta_xyyyz_xxyyz_0[i] =
            2.0 * ta_yyyz_xyyz_0[i] * fe_0 - 2.0 * ta_yyyz_xyyz_1[i] * fe_0 + ta_yyyz_xxyyz_0[i] * pa_x[i] - ta_yyyz_xxyyz_1[i] * pc_x[i];

        ta_xyyyz_xxyzz_0[i] =
            2.0 * ta_yyyz_xyzz_0[i] * fe_0 - 2.0 * ta_yyyz_xyzz_1[i] * fe_0 + ta_yyyz_xxyzz_0[i] * pa_x[i] - ta_yyyz_xxyzz_1[i] * pc_x[i];

        ta_xyyyz_xxzzz_0[i] =
            2.0 * ta_yyyz_xzzz_0[i] * fe_0 - 2.0 * ta_yyyz_xzzz_1[i] * fe_0 + ta_yyyz_xxzzz_0[i] * pa_x[i] - ta_yyyz_xxzzz_1[i] * pc_x[i];

        ta_xyyyz_xyyyy_0[i] = ta_xyyy_xyyyy_0[i] * pa_z[i] - ta_xyyy_xyyyy_1[i] * pc_z[i];

        ta_xyyyz_xyyyz_0[i] = ta_yyyz_yyyz_0[i] * fe_0 - ta_yyyz_yyyz_1[i] * fe_0 + ta_yyyz_xyyyz_0[i] * pa_x[i] - ta_yyyz_xyyyz_1[i] * pc_x[i];

        ta_xyyyz_xyyzz_0[i] = ta_yyyz_yyzz_0[i] * fe_0 - ta_yyyz_yyzz_1[i] * fe_0 + ta_yyyz_xyyzz_0[i] * pa_x[i] - ta_yyyz_xyyzz_1[i] * pc_x[i];

        ta_xyyyz_xyzzz_0[i] = ta_yyyz_yzzz_0[i] * fe_0 - ta_yyyz_yzzz_1[i] * fe_0 + ta_yyyz_xyzzz_0[i] * pa_x[i] - ta_yyyz_xyzzz_1[i] * pc_x[i];

        ta_xyyyz_xzzzz_0[i] = ta_yyyz_zzzz_0[i] * fe_0 - ta_yyyz_zzzz_1[i] * fe_0 + ta_yyyz_xzzzz_0[i] * pa_x[i] - ta_yyyz_xzzzz_1[i] * pc_x[i];

        ta_xyyyz_yyyyy_0[i] = ta_yyyz_yyyyy_0[i] * pa_x[i] - ta_yyyz_yyyyy_1[i] * pc_x[i];

        ta_xyyyz_yyyyz_0[i] = ta_yyyz_yyyyz_0[i] * pa_x[i] - ta_yyyz_yyyyz_1[i] * pc_x[i];

        ta_xyyyz_yyyzz_0[i] = ta_yyyz_yyyzz_0[i] * pa_x[i] - ta_yyyz_yyyzz_1[i] * pc_x[i];

        ta_xyyyz_yyzzz_0[i] = ta_yyyz_yyzzz_0[i] * pa_x[i] - ta_yyyz_yyzzz_1[i] * pc_x[i];

        ta_xyyyz_yzzzz_0[i] = ta_yyyz_yzzzz_0[i] * pa_x[i] - ta_yyyz_yzzzz_1[i] * pc_x[i];

        ta_xyyyz_zzzzz_0[i] = ta_yyyz_zzzzz_0[i] * pa_x[i] - ta_yyyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 252-273 components of targeted buffer : HH

    auto ta_xyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 252);

    auto ta_xyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 253);

    auto ta_xyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 254);

    auto ta_xyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 255);

    auto ta_xyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 256);

    auto ta_xyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 257);

    auto ta_xyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 258);

    auto ta_xyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 259);

    auto ta_xyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 260);

    auto ta_xyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 261);

    auto ta_xyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 262);

    auto ta_xyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 263);

    auto ta_xyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 264);

    auto ta_xyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 265);

    auto ta_xyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 266);

    auto ta_xyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 267);

    auto ta_xyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 268);

    auto ta_xyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 269);

    auto ta_xyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 270);

    auto ta_xyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 271);

    auto ta_xyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 272);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xyyzz_xxxxx_0, \
                             ta_xyyzz_xxxxy_0, \
                             ta_xyyzz_xxxxz_0, \
                             ta_xyyzz_xxxyy_0, \
                             ta_xyyzz_xxxyz_0, \
                             ta_xyyzz_xxxzz_0, \
                             ta_xyyzz_xxyyy_0, \
                             ta_xyyzz_xxyyz_0, \
                             ta_xyyzz_xxyzz_0, \
                             ta_xyyzz_xxzzz_0, \
                             ta_xyyzz_xyyyy_0, \
                             ta_xyyzz_xyyyz_0, \
                             ta_xyyzz_xyyzz_0, \
                             ta_xyyzz_xyzzz_0, \
                             ta_xyyzz_xzzzz_0, \
                             ta_xyyzz_yyyyy_0, \
                             ta_xyyzz_yyyyz_0, \
                             ta_xyyzz_yyyzz_0, \
                             ta_xyyzz_yyzzz_0, \
                             ta_xyyzz_yzzzz_0, \
                             ta_xyyzz_zzzzz_0, \
                             ta_yyzz_xxxx_0,   \
                             ta_yyzz_xxxx_1,   \
                             ta_yyzz_xxxxx_0,  \
                             ta_yyzz_xxxxx_1,  \
                             ta_yyzz_xxxxy_0,  \
                             ta_yyzz_xxxxy_1,  \
                             ta_yyzz_xxxxz_0,  \
                             ta_yyzz_xxxxz_1,  \
                             ta_yyzz_xxxy_0,   \
                             ta_yyzz_xxxy_1,   \
                             ta_yyzz_xxxyy_0,  \
                             ta_yyzz_xxxyy_1,  \
                             ta_yyzz_xxxyz_0,  \
                             ta_yyzz_xxxyz_1,  \
                             ta_yyzz_xxxz_0,   \
                             ta_yyzz_xxxz_1,   \
                             ta_yyzz_xxxzz_0,  \
                             ta_yyzz_xxxzz_1,  \
                             ta_yyzz_xxyy_0,   \
                             ta_yyzz_xxyy_1,   \
                             ta_yyzz_xxyyy_0,  \
                             ta_yyzz_xxyyy_1,  \
                             ta_yyzz_xxyyz_0,  \
                             ta_yyzz_xxyyz_1,  \
                             ta_yyzz_xxyz_0,   \
                             ta_yyzz_xxyz_1,   \
                             ta_yyzz_xxyzz_0,  \
                             ta_yyzz_xxyzz_1,  \
                             ta_yyzz_xxzz_0,   \
                             ta_yyzz_xxzz_1,   \
                             ta_yyzz_xxzzz_0,  \
                             ta_yyzz_xxzzz_1,  \
                             ta_yyzz_xyyy_0,   \
                             ta_yyzz_xyyy_1,   \
                             ta_yyzz_xyyyy_0,  \
                             ta_yyzz_xyyyy_1,  \
                             ta_yyzz_xyyyz_0,  \
                             ta_yyzz_xyyyz_1,  \
                             ta_yyzz_xyyz_0,   \
                             ta_yyzz_xyyz_1,   \
                             ta_yyzz_xyyzz_0,  \
                             ta_yyzz_xyyzz_1,  \
                             ta_yyzz_xyzz_0,   \
                             ta_yyzz_xyzz_1,   \
                             ta_yyzz_xyzzz_0,  \
                             ta_yyzz_xyzzz_1,  \
                             ta_yyzz_xzzz_0,   \
                             ta_yyzz_xzzz_1,   \
                             ta_yyzz_xzzzz_0,  \
                             ta_yyzz_xzzzz_1,  \
                             ta_yyzz_yyyy_0,   \
                             ta_yyzz_yyyy_1,   \
                             ta_yyzz_yyyyy_0,  \
                             ta_yyzz_yyyyy_1,  \
                             ta_yyzz_yyyyz_0,  \
                             ta_yyzz_yyyyz_1,  \
                             ta_yyzz_yyyz_0,   \
                             ta_yyzz_yyyz_1,   \
                             ta_yyzz_yyyzz_0,  \
                             ta_yyzz_yyyzz_1,  \
                             ta_yyzz_yyzz_0,   \
                             ta_yyzz_yyzz_1,   \
                             ta_yyzz_yyzzz_0,  \
                             ta_yyzz_yyzzz_1,  \
                             ta_yyzz_yzzz_0,   \
                             ta_yyzz_yzzz_1,   \
                             ta_yyzz_yzzzz_0,  \
                             ta_yyzz_yzzzz_1,  \
                             ta_yyzz_zzzz_0,   \
                             ta_yyzz_zzzz_1,   \
                             ta_yyzz_zzzzz_0,  \
                             ta_yyzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzz_xxxxx_0[i] =
            5.0 * ta_yyzz_xxxx_0[i] * fe_0 - 5.0 * ta_yyzz_xxxx_1[i] * fe_0 + ta_yyzz_xxxxx_0[i] * pa_x[i] - ta_yyzz_xxxxx_1[i] * pc_x[i];

        ta_xyyzz_xxxxy_0[i] =
            4.0 * ta_yyzz_xxxy_0[i] * fe_0 - 4.0 * ta_yyzz_xxxy_1[i] * fe_0 + ta_yyzz_xxxxy_0[i] * pa_x[i] - ta_yyzz_xxxxy_1[i] * pc_x[i];

        ta_xyyzz_xxxxz_0[i] =
            4.0 * ta_yyzz_xxxz_0[i] * fe_0 - 4.0 * ta_yyzz_xxxz_1[i] * fe_0 + ta_yyzz_xxxxz_0[i] * pa_x[i] - ta_yyzz_xxxxz_1[i] * pc_x[i];

        ta_xyyzz_xxxyy_0[i] =
            3.0 * ta_yyzz_xxyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxyy_1[i] * fe_0 + ta_yyzz_xxxyy_0[i] * pa_x[i] - ta_yyzz_xxxyy_1[i] * pc_x[i];

        ta_xyyzz_xxxyz_0[i] =
            3.0 * ta_yyzz_xxyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyz_1[i] * fe_0 + ta_yyzz_xxxyz_0[i] * pa_x[i] - ta_yyzz_xxxyz_1[i] * pc_x[i];

        ta_xyyzz_xxxzz_0[i] =
            3.0 * ta_yyzz_xxzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxzz_1[i] * fe_0 + ta_yyzz_xxxzz_0[i] * pa_x[i] - ta_yyzz_xxxzz_1[i] * pc_x[i];

        ta_xyyzz_xxyyy_0[i] =
            2.0 * ta_yyzz_xyyy_0[i] * fe_0 - 2.0 * ta_yyzz_xyyy_1[i] * fe_0 + ta_yyzz_xxyyy_0[i] * pa_x[i] - ta_yyzz_xxyyy_1[i] * pc_x[i];

        ta_xyyzz_xxyyz_0[i] =
            2.0 * ta_yyzz_xyyz_0[i] * fe_0 - 2.0 * ta_yyzz_xyyz_1[i] * fe_0 + ta_yyzz_xxyyz_0[i] * pa_x[i] - ta_yyzz_xxyyz_1[i] * pc_x[i];

        ta_xyyzz_xxyzz_0[i] =
            2.0 * ta_yyzz_xyzz_0[i] * fe_0 - 2.0 * ta_yyzz_xyzz_1[i] * fe_0 + ta_yyzz_xxyzz_0[i] * pa_x[i] - ta_yyzz_xxyzz_1[i] * pc_x[i];

        ta_xyyzz_xxzzz_0[i] =
            2.0 * ta_yyzz_xzzz_0[i] * fe_0 - 2.0 * ta_yyzz_xzzz_1[i] * fe_0 + ta_yyzz_xxzzz_0[i] * pa_x[i] - ta_yyzz_xxzzz_1[i] * pc_x[i];

        ta_xyyzz_xyyyy_0[i] = ta_yyzz_yyyy_0[i] * fe_0 - ta_yyzz_yyyy_1[i] * fe_0 + ta_yyzz_xyyyy_0[i] * pa_x[i] - ta_yyzz_xyyyy_1[i] * pc_x[i];

        ta_xyyzz_xyyyz_0[i] = ta_yyzz_yyyz_0[i] * fe_0 - ta_yyzz_yyyz_1[i] * fe_0 + ta_yyzz_xyyyz_0[i] * pa_x[i] - ta_yyzz_xyyyz_1[i] * pc_x[i];

        ta_xyyzz_xyyzz_0[i] = ta_yyzz_yyzz_0[i] * fe_0 - ta_yyzz_yyzz_1[i] * fe_0 + ta_yyzz_xyyzz_0[i] * pa_x[i] - ta_yyzz_xyyzz_1[i] * pc_x[i];

        ta_xyyzz_xyzzz_0[i] = ta_yyzz_yzzz_0[i] * fe_0 - ta_yyzz_yzzz_1[i] * fe_0 + ta_yyzz_xyzzz_0[i] * pa_x[i] - ta_yyzz_xyzzz_1[i] * pc_x[i];

        ta_xyyzz_xzzzz_0[i] = ta_yyzz_zzzz_0[i] * fe_0 - ta_yyzz_zzzz_1[i] * fe_0 + ta_yyzz_xzzzz_0[i] * pa_x[i] - ta_yyzz_xzzzz_1[i] * pc_x[i];

        ta_xyyzz_yyyyy_0[i] = ta_yyzz_yyyyy_0[i] * pa_x[i] - ta_yyzz_yyyyy_1[i] * pc_x[i];

        ta_xyyzz_yyyyz_0[i] = ta_yyzz_yyyyz_0[i] * pa_x[i] - ta_yyzz_yyyyz_1[i] * pc_x[i];

        ta_xyyzz_yyyzz_0[i] = ta_yyzz_yyyzz_0[i] * pa_x[i] - ta_yyzz_yyyzz_1[i] * pc_x[i];

        ta_xyyzz_yyzzz_0[i] = ta_yyzz_yyzzz_0[i] * pa_x[i] - ta_yyzz_yyzzz_1[i] * pc_x[i];

        ta_xyyzz_yzzzz_0[i] = ta_yyzz_yzzzz_0[i] * pa_x[i] - ta_yyzz_yzzzz_1[i] * pc_x[i];

        ta_xyyzz_zzzzz_0[i] = ta_yyzz_zzzzz_0[i] * pa_x[i] - ta_yyzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 273-294 components of targeted buffer : HH

    auto ta_xyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 273);

    auto ta_xyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 274);

    auto ta_xyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 275);

    auto ta_xyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 276);

    auto ta_xyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 277);

    auto ta_xyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 278);

    auto ta_xyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 279);

    auto ta_xyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 280);

    auto ta_xyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 281);

    auto ta_xyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 282);

    auto ta_xyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 283);

    auto ta_xyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 284);

    auto ta_xyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 285);

    auto ta_xyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 286);

    auto ta_xyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 287);

    auto ta_xyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 288);

    auto ta_xyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 289);

    auto ta_xyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 290);

    auto ta_xyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 291);

    auto ta_xyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 292);

    auto ta_xyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 293);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta_xyzzz_xxxxx_0, \
                             ta_xyzzz_xxxxy_0, \
                             ta_xyzzz_xxxxz_0, \
                             ta_xyzzz_xxxyy_0, \
                             ta_xyzzz_xxxyz_0, \
                             ta_xyzzz_xxxzz_0, \
                             ta_xyzzz_xxyyy_0, \
                             ta_xyzzz_xxyyz_0, \
                             ta_xyzzz_xxyzz_0, \
                             ta_xyzzz_xxzzz_0, \
                             ta_xyzzz_xyyyy_0, \
                             ta_xyzzz_xyyyz_0, \
                             ta_xyzzz_xyyzz_0, \
                             ta_xyzzz_xyzzz_0, \
                             ta_xyzzz_xzzzz_0, \
                             ta_xyzzz_yyyyy_0, \
                             ta_xyzzz_yyyyz_0, \
                             ta_xyzzz_yyyzz_0, \
                             ta_xyzzz_yyzzz_0, \
                             ta_xyzzz_yzzzz_0, \
                             ta_xyzzz_zzzzz_0, \
                             ta_xzzz_xxxxx_0,  \
                             ta_xzzz_xxxxx_1,  \
                             ta_xzzz_xxxxz_0,  \
                             ta_xzzz_xxxxz_1,  \
                             ta_xzzz_xxxzz_0,  \
                             ta_xzzz_xxxzz_1,  \
                             ta_xzzz_xxzzz_0,  \
                             ta_xzzz_xxzzz_1,  \
                             ta_xzzz_xzzzz_0,  \
                             ta_xzzz_xzzzz_1,  \
                             ta_yzzz_xxxxy_0,  \
                             ta_yzzz_xxxxy_1,  \
                             ta_yzzz_xxxy_0,   \
                             ta_yzzz_xxxy_1,   \
                             ta_yzzz_xxxyy_0,  \
                             ta_yzzz_xxxyy_1,  \
                             ta_yzzz_xxxyz_0,  \
                             ta_yzzz_xxxyz_1,  \
                             ta_yzzz_xxyy_0,   \
                             ta_yzzz_xxyy_1,   \
                             ta_yzzz_xxyyy_0,  \
                             ta_yzzz_xxyyy_1,  \
                             ta_yzzz_xxyyz_0,  \
                             ta_yzzz_xxyyz_1,  \
                             ta_yzzz_xxyz_0,   \
                             ta_yzzz_xxyz_1,   \
                             ta_yzzz_xxyzz_0,  \
                             ta_yzzz_xxyzz_1,  \
                             ta_yzzz_xyyy_0,   \
                             ta_yzzz_xyyy_1,   \
                             ta_yzzz_xyyyy_0,  \
                             ta_yzzz_xyyyy_1,  \
                             ta_yzzz_xyyyz_0,  \
                             ta_yzzz_xyyyz_1,  \
                             ta_yzzz_xyyz_0,   \
                             ta_yzzz_xyyz_1,   \
                             ta_yzzz_xyyzz_0,  \
                             ta_yzzz_xyyzz_1,  \
                             ta_yzzz_xyzz_0,   \
                             ta_yzzz_xyzz_1,   \
                             ta_yzzz_xyzzz_0,  \
                             ta_yzzz_xyzzz_1,  \
                             ta_yzzz_yyyy_0,   \
                             ta_yzzz_yyyy_1,   \
                             ta_yzzz_yyyyy_0,  \
                             ta_yzzz_yyyyy_1,  \
                             ta_yzzz_yyyyz_0,  \
                             ta_yzzz_yyyyz_1,  \
                             ta_yzzz_yyyz_0,   \
                             ta_yzzz_yyyz_1,   \
                             ta_yzzz_yyyzz_0,  \
                             ta_yzzz_yyyzz_1,  \
                             ta_yzzz_yyzz_0,   \
                             ta_yzzz_yyzz_1,   \
                             ta_yzzz_yyzzz_0,  \
                             ta_yzzz_yyzzz_1,  \
                             ta_yzzz_yzzz_0,   \
                             ta_yzzz_yzzz_1,   \
                             ta_yzzz_yzzzz_0,  \
                             ta_yzzz_yzzzz_1,  \
                             ta_yzzz_zzzzz_0,  \
                             ta_yzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzz_xxxxx_0[i] = ta_xzzz_xxxxx_0[i] * pa_y[i] - ta_xzzz_xxxxx_1[i] * pc_y[i];

        ta_xyzzz_xxxxy_0[i] =
            4.0 * ta_yzzz_xxxy_0[i] * fe_0 - 4.0 * ta_yzzz_xxxy_1[i] * fe_0 + ta_yzzz_xxxxy_0[i] * pa_x[i] - ta_yzzz_xxxxy_1[i] * pc_x[i];

        ta_xyzzz_xxxxz_0[i] = ta_xzzz_xxxxz_0[i] * pa_y[i] - ta_xzzz_xxxxz_1[i] * pc_y[i];

        ta_xyzzz_xxxyy_0[i] =
            3.0 * ta_yzzz_xxyy_0[i] * fe_0 - 3.0 * ta_yzzz_xxyy_1[i] * fe_0 + ta_yzzz_xxxyy_0[i] * pa_x[i] - ta_yzzz_xxxyy_1[i] * pc_x[i];

        ta_xyzzz_xxxyz_0[i] =
            3.0 * ta_yzzz_xxyz_0[i] * fe_0 - 3.0 * ta_yzzz_xxyz_1[i] * fe_0 + ta_yzzz_xxxyz_0[i] * pa_x[i] - ta_yzzz_xxxyz_1[i] * pc_x[i];

        ta_xyzzz_xxxzz_0[i] = ta_xzzz_xxxzz_0[i] * pa_y[i] - ta_xzzz_xxxzz_1[i] * pc_y[i];

        ta_xyzzz_xxyyy_0[i] =
            2.0 * ta_yzzz_xyyy_0[i] * fe_0 - 2.0 * ta_yzzz_xyyy_1[i] * fe_0 + ta_yzzz_xxyyy_0[i] * pa_x[i] - ta_yzzz_xxyyy_1[i] * pc_x[i];

        ta_xyzzz_xxyyz_0[i] =
            2.0 * ta_yzzz_xyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyz_1[i] * fe_0 + ta_yzzz_xxyyz_0[i] * pa_x[i] - ta_yzzz_xxyyz_1[i] * pc_x[i];

        ta_xyzzz_xxyzz_0[i] =
            2.0 * ta_yzzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyzz_1[i] * fe_0 + ta_yzzz_xxyzz_0[i] * pa_x[i] - ta_yzzz_xxyzz_1[i] * pc_x[i];

        ta_xyzzz_xxzzz_0[i] = ta_xzzz_xxzzz_0[i] * pa_y[i] - ta_xzzz_xxzzz_1[i] * pc_y[i];

        ta_xyzzz_xyyyy_0[i] = ta_yzzz_yyyy_0[i] * fe_0 - ta_yzzz_yyyy_1[i] * fe_0 + ta_yzzz_xyyyy_0[i] * pa_x[i] - ta_yzzz_xyyyy_1[i] * pc_x[i];

        ta_xyzzz_xyyyz_0[i] = ta_yzzz_yyyz_0[i] * fe_0 - ta_yzzz_yyyz_1[i] * fe_0 + ta_yzzz_xyyyz_0[i] * pa_x[i] - ta_yzzz_xyyyz_1[i] * pc_x[i];

        ta_xyzzz_xyyzz_0[i] = ta_yzzz_yyzz_0[i] * fe_0 - ta_yzzz_yyzz_1[i] * fe_0 + ta_yzzz_xyyzz_0[i] * pa_x[i] - ta_yzzz_xyyzz_1[i] * pc_x[i];

        ta_xyzzz_xyzzz_0[i] = ta_yzzz_yzzz_0[i] * fe_0 - ta_yzzz_yzzz_1[i] * fe_0 + ta_yzzz_xyzzz_0[i] * pa_x[i] - ta_yzzz_xyzzz_1[i] * pc_x[i];

        ta_xyzzz_xzzzz_0[i] = ta_xzzz_xzzzz_0[i] * pa_y[i] - ta_xzzz_xzzzz_1[i] * pc_y[i];

        ta_xyzzz_yyyyy_0[i] = ta_yzzz_yyyyy_0[i] * pa_x[i] - ta_yzzz_yyyyy_1[i] * pc_x[i];

        ta_xyzzz_yyyyz_0[i] = ta_yzzz_yyyyz_0[i] * pa_x[i] - ta_yzzz_yyyyz_1[i] * pc_x[i];

        ta_xyzzz_yyyzz_0[i] = ta_yzzz_yyyzz_0[i] * pa_x[i] - ta_yzzz_yyyzz_1[i] * pc_x[i];

        ta_xyzzz_yyzzz_0[i] = ta_yzzz_yyzzz_0[i] * pa_x[i] - ta_yzzz_yyzzz_1[i] * pc_x[i];

        ta_xyzzz_yzzzz_0[i] = ta_yzzz_yzzzz_0[i] * pa_x[i] - ta_yzzz_yzzzz_1[i] * pc_x[i];

        ta_xyzzz_zzzzz_0[i] = ta_yzzz_zzzzz_0[i] * pa_x[i] - ta_yzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 294-315 components of targeted buffer : HH

    auto ta_xzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 294);

    auto ta_xzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 295);

    auto ta_xzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 296);

    auto ta_xzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 297);

    auto ta_xzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 298);

    auto ta_xzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 299);

    auto ta_xzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 300);

    auto ta_xzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 301);

    auto ta_xzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 302);

    auto ta_xzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 303);

    auto ta_xzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 304);

    auto ta_xzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 305);

    auto ta_xzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 306);

    auto ta_xzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 307);

    auto ta_xzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 308);

    auto ta_xzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 309);

    auto ta_xzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 310);

    auto ta_xzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 311);

    auto ta_xzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 312);

    auto ta_xzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 313);

    auto ta_xzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 314);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta_xzzzz_xxxxx_0, \
                             ta_xzzzz_xxxxy_0, \
                             ta_xzzzz_xxxxz_0, \
                             ta_xzzzz_xxxyy_0, \
                             ta_xzzzz_xxxyz_0, \
                             ta_xzzzz_xxxzz_0, \
                             ta_xzzzz_xxyyy_0, \
                             ta_xzzzz_xxyyz_0, \
                             ta_xzzzz_xxyzz_0, \
                             ta_xzzzz_xxzzz_0, \
                             ta_xzzzz_xyyyy_0, \
                             ta_xzzzz_xyyyz_0, \
                             ta_xzzzz_xyyzz_0, \
                             ta_xzzzz_xyzzz_0, \
                             ta_xzzzz_xzzzz_0, \
                             ta_xzzzz_yyyyy_0, \
                             ta_xzzzz_yyyyz_0, \
                             ta_xzzzz_yyyzz_0, \
                             ta_xzzzz_yyzzz_0, \
                             ta_xzzzz_yzzzz_0, \
                             ta_xzzzz_zzzzz_0, \
                             ta_zzzz_xxxx_0,   \
                             ta_zzzz_xxxx_1,   \
                             ta_zzzz_xxxxx_0,  \
                             ta_zzzz_xxxxx_1,  \
                             ta_zzzz_xxxxy_0,  \
                             ta_zzzz_xxxxy_1,  \
                             ta_zzzz_xxxxz_0,  \
                             ta_zzzz_xxxxz_1,  \
                             ta_zzzz_xxxy_0,   \
                             ta_zzzz_xxxy_1,   \
                             ta_zzzz_xxxyy_0,  \
                             ta_zzzz_xxxyy_1,  \
                             ta_zzzz_xxxyz_0,  \
                             ta_zzzz_xxxyz_1,  \
                             ta_zzzz_xxxz_0,   \
                             ta_zzzz_xxxz_1,   \
                             ta_zzzz_xxxzz_0,  \
                             ta_zzzz_xxxzz_1,  \
                             ta_zzzz_xxyy_0,   \
                             ta_zzzz_xxyy_1,   \
                             ta_zzzz_xxyyy_0,  \
                             ta_zzzz_xxyyy_1,  \
                             ta_zzzz_xxyyz_0,  \
                             ta_zzzz_xxyyz_1,  \
                             ta_zzzz_xxyz_0,   \
                             ta_zzzz_xxyz_1,   \
                             ta_zzzz_xxyzz_0,  \
                             ta_zzzz_xxyzz_1,  \
                             ta_zzzz_xxzz_0,   \
                             ta_zzzz_xxzz_1,   \
                             ta_zzzz_xxzzz_0,  \
                             ta_zzzz_xxzzz_1,  \
                             ta_zzzz_xyyy_0,   \
                             ta_zzzz_xyyy_1,   \
                             ta_zzzz_xyyyy_0,  \
                             ta_zzzz_xyyyy_1,  \
                             ta_zzzz_xyyyz_0,  \
                             ta_zzzz_xyyyz_1,  \
                             ta_zzzz_xyyz_0,   \
                             ta_zzzz_xyyz_1,   \
                             ta_zzzz_xyyzz_0,  \
                             ta_zzzz_xyyzz_1,  \
                             ta_zzzz_xyzz_0,   \
                             ta_zzzz_xyzz_1,   \
                             ta_zzzz_xyzzz_0,  \
                             ta_zzzz_xyzzz_1,  \
                             ta_zzzz_xzzz_0,   \
                             ta_zzzz_xzzz_1,   \
                             ta_zzzz_xzzzz_0,  \
                             ta_zzzz_xzzzz_1,  \
                             ta_zzzz_yyyy_0,   \
                             ta_zzzz_yyyy_1,   \
                             ta_zzzz_yyyyy_0,  \
                             ta_zzzz_yyyyy_1,  \
                             ta_zzzz_yyyyz_0,  \
                             ta_zzzz_yyyyz_1,  \
                             ta_zzzz_yyyz_0,   \
                             ta_zzzz_yyyz_1,   \
                             ta_zzzz_yyyzz_0,  \
                             ta_zzzz_yyyzz_1,  \
                             ta_zzzz_yyzz_0,   \
                             ta_zzzz_yyzz_1,   \
                             ta_zzzz_yyzzz_0,  \
                             ta_zzzz_yyzzz_1,  \
                             ta_zzzz_yzzz_0,   \
                             ta_zzzz_yzzz_1,   \
                             ta_zzzz_yzzzz_0,  \
                             ta_zzzz_yzzzz_1,  \
                             ta_zzzz_zzzz_0,   \
                             ta_zzzz_zzzz_1,   \
                             ta_zzzz_zzzzz_0,  \
                             ta_zzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzz_xxxxx_0[i] =
            5.0 * ta_zzzz_xxxx_0[i] * fe_0 - 5.0 * ta_zzzz_xxxx_1[i] * fe_0 + ta_zzzz_xxxxx_0[i] * pa_x[i] - ta_zzzz_xxxxx_1[i] * pc_x[i];

        ta_xzzzz_xxxxy_0[i] =
            4.0 * ta_zzzz_xxxy_0[i] * fe_0 - 4.0 * ta_zzzz_xxxy_1[i] * fe_0 + ta_zzzz_xxxxy_0[i] * pa_x[i] - ta_zzzz_xxxxy_1[i] * pc_x[i];

        ta_xzzzz_xxxxz_0[i] =
            4.0 * ta_zzzz_xxxz_0[i] * fe_0 - 4.0 * ta_zzzz_xxxz_1[i] * fe_0 + ta_zzzz_xxxxz_0[i] * pa_x[i] - ta_zzzz_xxxxz_1[i] * pc_x[i];

        ta_xzzzz_xxxyy_0[i] =
            3.0 * ta_zzzz_xxyy_0[i] * fe_0 - 3.0 * ta_zzzz_xxyy_1[i] * fe_0 + ta_zzzz_xxxyy_0[i] * pa_x[i] - ta_zzzz_xxxyy_1[i] * pc_x[i];

        ta_xzzzz_xxxyz_0[i] =
            3.0 * ta_zzzz_xxyz_0[i] * fe_0 - 3.0 * ta_zzzz_xxyz_1[i] * fe_0 + ta_zzzz_xxxyz_0[i] * pa_x[i] - ta_zzzz_xxxyz_1[i] * pc_x[i];

        ta_xzzzz_xxxzz_0[i] =
            3.0 * ta_zzzz_xxzz_0[i] * fe_0 - 3.0 * ta_zzzz_xxzz_1[i] * fe_0 + ta_zzzz_xxxzz_0[i] * pa_x[i] - ta_zzzz_xxxzz_1[i] * pc_x[i];

        ta_xzzzz_xxyyy_0[i] =
            2.0 * ta_zzzz_xyyy_0[i] * fe_0 - 2.0 * ta_zzzz_xyyy_1[i] * fe_0 + ta_zzzz_xxyyy_0[i] * pa_x[i] - ta_zzzz_xxyyy_1[i] * pc_x[i];

        ta_xzzzz_xxyyz_0[i] =
            2.0 * ta_zzzz_xyyz_0[i] * fe_0 - 2.0 * ta_zzzz_xyyz_1[i] * fe_0 + ta_zzzz_xxyyz_0[i] * pa_x[i] - ta_zzzz_xxyyz_1[i] * pc_x[i];

        ta_xzzzz_xxyzz_0[i] =
            2.0 * ta_zzzz_xyzz_0[i] * fe_0 - 2.0 * ta_zzzz_xyzz_1[i] * fe_0 + ta_zzzz_xxyzz_0[i] * pa_x[i] - ta_zzzz_xxyzz_1[i] * pc_x[i];

        ta_xzzzz_xxzzz_0[i] =
            2.0 * ta_zzzz_xzzz_0[i] * fe_0 - 2.0 * ta_zzzz_xzzz_1[i] * fe_0 + ta_zzzz_xxzzz_0[i] * pa_x[i] - ta_zzzz_xxzzz_1[i] * pc_x[i];

        ta_xzzzz_xyyyy_0[i] = ta_zzzz_yyyy_0[i] * fe_0 - ta_zzzz_yyyy_1[i] * fe_0 + ta_zzzz_xyyyy_0[i] * pa_x[i] - ta_zzzz_xyyyy_1[i] * pc_x[i];

        ta_xzzzz_xyyyz_0[i] = ta_zzzz_yyyz_0[i] * fe_0 - ta_zzzz_yyyz_1[i] * fe_0 + ta_zzzz_xyyyz_0[i] * pa_x[i] - ta_zzzz_xyyyz_1[i] * pc_x[i];

        ta_xzzzz_xyyzz_0[i] = ta_zzzz_yyzz_0[i] * fe_0 - ta_zzzz_yyzz_1[i] * fe_0 + ta_zzzz_xyyzz_0[i] * pa_x[i] - ta_zzzz_xyyzz_1[i] * pc_x[i];

        ta_xzzzz_xyzzz_0[i] = ta_zzzz_yzzz_0[i] * fe_0 - ta_zzzz_yzzz_1[i] * fe_0 + ta_zzzz_xyzzz_0[i] * pa_x[i] - ta_zzzz_xyzzz_1[i] * pc_x[i];

        ta_xzzzz_xzzzz_0[i] = ta_zzzz_zzzz_0[i] * fe_0 - ta_zzzz_zzzz_1[i] * fe_0 + ta_zzzz_xzzzz_0[i] * pa_x[i] - ta_zzzz_xzzzz_1[i] * pc_x[i];

        ta_xzzzz_yyyyy_0[i] = ta_zzzz_yyyyy_0[i] * pa_x[i] - ta_zzzz_yyyyy_1[i] * pc_x[i];

        ta_xzzzz_yyyyz_0[i] = ta_zzzz_yyyyz_0[i] * pa_x[i] - ta_zzzz_yyyyz_1[i] * pc_x[i];

        ta_xzzzz_yyyzz_0[i] = ta_zzzz_yyyzz_0[i] * pa_x[i] - ta_zzzz_yyyzz_1[i] * pc_x[i];

        ta_xzzzz_yyzzz_0[i] = ta_zzzz_yyzzz_0[i] * pa_x[i] - ta_zzzz_yyzzz_1[i] * pc_x[i];

        ta_xzzzz_yzzzz_0[i] = ta_zzzz_yzzzz_0[i] * pa_x[i] - ta_zzzz_yzzzz_1[i] * pc_x[i];

        ta_xzzzz_zzzzz_0[i] = ta_zzzz_zzzzz_0[i] * pa_x[i] - ta_zzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 315-336 components of targeted buffer : HH

    auto ta_yyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 315);

    auto ta_yyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 316);

    auto ta_yyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 317);

    auto ta_yyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 318);

    auto ta_yyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 319);

    auto ta_yyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 320);

    auto ta_yyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 321);

    auto ta_yyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 322);

    auto ta_yyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 323);

    auto ta_yyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 324);

    auto ta_yyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 325);

    auto ta_yyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 326);

    auto ta_yyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 327);

    auto ta_yyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 328);

    auto ta_yyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 329);

    auto ta_yyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 330);

    auto ta_yyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 331);

    auto ta_yyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 332);

    auto ta_yyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 333);

    auto ta_yyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 334);

    auto ta_yyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 335);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta_yyy_xxxxx_0,   \
                             ta_yyy_xxxxx_1,   \
                             ta_yyy_xxxxy_0,   \
                             ta_yyy_xxxxy_1,   \
                             ta_yyy_xxxxz_0,   \
                             ta_yyy_xxxxz_1,   \
                             ta_yyy_xxxyy_0,   \
                             ta_yyy_xxxyy_1,   \
                             ta_yyy_xxxyz_0,   \
                             ta_yyy_xxxyz_1,   \
                             ta_yyy_xxxzz_0,   \
                             ta_yyy_xxxzz_1,   \
                             ta_yyy_xxyyy_0,   \
                             ta_yyy_xxyyy_1,   \
                             ta_yyy_xxyyz_0,   \
                             ta_yyy_xxyyz_1,   \
                             ta_yyy_xxyzz_0,   \
                             ta_yyy_xxyzz_1,   \
                             ta_yyy_xxzzz_0,   \
                             ta_yyy_xxzzz_1,   \
                             ta_yyy_xyyyy_0,   \
                             ta_yyy_xyyyy_1,   \
                             ta_yyy_xyyyz_0,   \
                             ta_yyy_xyyyz_1,   \
                             ta_yyy_xyyzz_0,   \
                             ta_yyy_xyyzz_1,   \
                             ta_yyy_xyzzz_0,   \
                             ta_yyy_xyzzz_1,   \
                             ta_yyy_xzzzz_0,   \
                             ta_yyy_xzzzz_1,   \
                             ta_yyy_yyyyy_0,   \
                             ta_yyy_yyyyy_1,   \
                             ta_yyy_yyyyz_0,   \
                             ta_yyy_yyyyz_1,   \
                             ta_yyy_yyyzz_0,   \
                             ta_yyy_yyyzz_1,   \
                             ta_yyy_yyzzz_0,   \
                             ta_yyy_yyzzz_1,   \
                             ta_yyy_yzzzz_0,   \
                             ta_yyy_yzzzz_1,   \
                             ta_yyy_zzzzz_0,   \
                             ta_yyy_zzzzz_1,   \
                             ta_yyyy_xxxx_0,   \
                             ta_yyyy_xxxx_1,   \
                             ta_yyyy_xxxxx_0,  \
                             ta_yyyy_xxxxx_1,  \
                             ta_yyyy_xxxxy_0,  \
                             ta_yyyy_xxxxy_1,  \
                             ta_yyyy_xxxxz_0,  \
                             ta_yyyy_xxxxz_1,  \
                             ta_yyyy_xxxy_0,   \
                             ta_yyyy_xxxy_1,   \
                             ta_yyyy_xxxyy_0,  \
                             ta_yyyy_xxxyy_1,  \
                             ta_yyyy_xxxyz_0,  \
                             ta_yyyy_xxxyz_1,  \
                             ta_yyyy_xxxz_0,   \
                             ta_yyyy_xxxz_1,   \
                             ta_yyyy_xxxzz_0,  \
                             ta_yyyy_xxxzz_1,  \
                             ta_yyyy_xxyy_0,   \
                             ta_yyyy_xxyy_1,   \
                             ta_yyyy_xxyyy_0,  \
                             ta_yyyy_xxyyy_1,  \
                             ta_yyyy_xxyyz_0,  \
                             ta_yyyy_xxyyz_1,  \
                             ta_yyyy_xxyz_0,   \
                             ta_yyyy_xxyz_1,   \
                             ta_yyyy_xxyzz_0,  \
                             ta_yyyy_xxyzz_1,  \
                             ta_yyyy_xxzz_0,   \
                             ta_yyyy_xxzz_1,   \
                             ta_yyyy_xxzzz_0,  \
                             ta_yyyy_xxzzz_1,  \
                             ta_yyyy_xyyy_0,   \
                             ta_yyyy_xyyy_1,   \
                             ta_yyyy_xyyyy_0,  \
                             ta_yyyy_xyyyy_1,  \
                             ta_yyyy_xyyyz_0,  \
                             ta_yyyy_xyyyz_1,  \
                             ta_yyyy_xyyz_0,   \
                             ta_yyyy_xyyz_1,   \
                             ta_yyyy_xyyzz_0,  \
                             ta_yyyy_xyyzz_1,  \
                             ta_yyyy_xyzz_0,   \
                             ta_yyyy_xyzz_1,   \
                             ta_yyyy_xyzzz_0,  \
                             ta_yyyy_xyzzz_1,  \
                             ta_yyyy_xzzz_0,   \
                             ta_yyyy_xzzz_1,   \
                             ta_yyyy_xzzzz_0,  \
                             ta_yyyy_xzzzz_1,  \
                             ta_yyyy_yyyy_0,   \
                             ta_yyyy_yyyy_1,   \
                             ta_yyyy_yyyyy_0,  \
                             ta_yyyy_yyyyy_1,  \
                             ta_yyyy_yyyyz_0,  \
                             ta_yyyy_yyyyz_1,  \
                             ta_yyyy_yyyz_0,   \
                             ta_yyyy_yyyz_1,   \
                             ta_yyyy_yyyzz_0,  \
                             ta_yyyy_yyyzz_1,  \
                             ta_yyyy_yyzz_0,   \
                             ta_yyyy_yyzz_1,   \
                             ta_yyyy_yyzzz_0,  \
                             ta_yyyy_yyzzz_1,  \
                             ta_yyyy_yzzz_0,   \
                             ta_yyyy_yzzz_1,   \
                             ta_yyyy_yzzzz_0,  \
                             ta_yyyy_yzzzz_1,  \
                             ta_yyyy_zzzz_0,   \
                             ta_yyyy_zzzz_1,   \
                             ta_yyyy_zzzzz_0,  \
                             ta_yyyy_zzzzz_1,  \
                             ta_yyyyy_xxxxx_0, \
                             ta_yyyyy_xxxxy_0, \
                             ta_yyyyy_xxxxz_0, \
                             ta_yyyyy_xxxyy_0, \
                             ta_yyyyy_xxxyz_0, \
                             ta_yyyyy_xxxzz_0, \
                             ta_yyyyy_xxyyy_0, \
                             ta_yyyyy_xxyyz_0, \
                             ta_yyyyy_xxyzz_0, \
                             ta_yyyyy_xxzzz_0, \
                             ta_yyyyy_xyyyy_0, \
                             ta_yyyyy_xyyyz_0, \
                             ta_yyyyy_xyyzz_0, \
                             ta_yyyyy_xyzzz_0, \
                             ta_yyyyy_xzzzz_0, \
                             ta_yyyyy_yyyyy_0, \
                             ta_yyyyy_yyyyz_0, \
                             ta_yyyyy_yyyzz_0, \
                             ta_yyyyy_yyzzz_0, \
                             ta_yyyyy_yzzzz_0, \
                             ta_yyyyy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyy_xxxxx_0[i] =
            4.0 * ta_yyy_xxxxx_0[i] * fe_0 - 4.0 * ta_yyy_xxxxx_1[i] * fe_0 + ta_yyyy_xxxxx_0[i] * pa_y[i] - ta_yyyy_xxxxx_1[i] * pc_y[i];

        ta_yyyyy_xxxxy_0[i] = 4.0 * ta_yyy_xxxxy_0[i] * fe_0 - 4.0 * ta_yyy_xxxxy_1[i] * fe_0 + ta_yyyy_xxxx_0[i] * fe_0 - ta_yyyy_xxxx_1[i] * fe_0 +
                              ta_yyyy_xxxxy_0[i] * pa_y[i] - ta_yyyy_xxxxy_1[i] * pc_y[i];

        ta_yyyyy_xxxxz_0[i] =
            4.0 * ta_yyy_xxxxz_0[i] * fe_0 - 4.0 * ta_yyy_xxxxz_1[i] * fe_0 + ta_yyyy_xxxxz_0[i] * pa_y[i] - ta_yyyy_xxxxz_1[i] * pc_y[i];

        ta_yyyyy_xxxyy_0[i] = 4.0 * ta_yyy_xxxyy_0[i] * fe_0 - 4.0 * ta_yyy_xxxyy_1[i] * fe_0 + 2.0 * ta_yyyy_xxxy_0[i] * fe_0 -
                              2.0 * ta_yyyy_xxxy_1[i] * fe_0 + ta_yyyy_xxxyy_0[i] * pa_y[i] - ta_yyyy_xxxyy_1[i] * pc_y[i];

        ta_yyyyy_xxxyz_0[i] = 4.0 * ta_yyy_xxxyz_0[i] * fe_0 - 4.0 * ta_yyy_xxxyz_1[i] * fe_0 + ta_yyyy_xxxz_0[i] * fe_0 - ta_yyyy_xxxz_1[i] * fe_0 +
                              ta_yyyy_xxxyz_0[i] * pa_y[i] - ta_yyyy_xxxyz_1[i] * pc_y[i];

        ta_yyyyy_xxxzz_0[i] =
            4.0 * ta_yyy_xxxzz_0[i] * fe_0 - 4.0 * ta_yyy_xxxzz_1[i] * fe_0 + ta_yyyy_xxxzz_0[i] * pa_y[i] - ta_yyyy_xxxzz_1[i] * pc_y[i];

        ta_yyyyy_xxyyy_0[i] = 4.0 * ta_yyy_xxyyy_0[i] * fe_0 - 4.0 * ta_yyy_xxyyy_1[i] * fe_0 + 3.0 * ta_yyyy_xxyy_0[i] * fe_0 -
                              3.0 * ta_yyyy_xxyy_1[i] * fe_0 + ta_yyyy_xxyyy_0[i] * pa_y[i] - ta_yyyy_xxyyy_1[i] * pc_y[i];

        ta_yyyyy_xxyyz_0[i] = 4.0 * ta_yyy_xxyyz_0[i] * fe_0 - 4.0 * ta_yyy_xxyyz_1[i] * fe_0 + 2.0 * ta_yyyy_xxyz_0[i] * fe_0 -
                              2.0 * ta_yyyy_xxyz_1[i] * fe_0 + ta_yyyy_xxyyz_0[i] * pa_y[i] - ta_yyyy_xxyyz_1[i] * pc_y[i];

        ta_yyyyy_xxyzz_0[i] = 4.0 * ta_yyy_xxyzz_0[i] * fe_0 - 4.0 * ta_yyy_xxyzz_1[i] * fe_0 + ta_yyyy_xxzz_0[i] * fe_0 - ta_yyyy_xxzz_1[i] * fe_0 +
                              ta_yyyy_xxyzz_0[i] * pa_y[i] - ta_yyyy_xxyzz_1[i] * pc_y[i];

        ta_yyyyy_xxzzz_0[i] =
            4.0 * ta_yyy_xxzzz_0[i] * fe_0 - 4.0 * ta_yyy_xxzzz_1[i] * fe_0 + ta_yyyy_xxzzz_0[i] * pa_y[i] - ta_yyyy_xxzzz_1[i] * pc_y[i];

        ta_yyyyy_xyyyy_0[i] = 4.0 * ta_yyy_xyyyy_0[i] * fe_0 - 4.0 * ta_yyy_xyyyy_1[i] * fe_0 + 4.0 * ta_yyyy_xyyy_0[i] * fe_0 -
                              4.0 * ta_yyyy_xyyy_1[i] * fe_0 + ta_yyyy_xyyyy_0[i] * pa_y[i] - ta_yyyy_xyyyy_1[i] * pc_y[i];

        ta_yyyyy_xyyyz_0[i] = 4.0 * ta_yyy_xyyyz_0[i] * fe_0 - 4.0 * ta_yyy_xyyyz_1[i] * fe_0 + 3.0 * ta_yyyy_xyyz_0[i] * fe_0 -
                              3.0 * ta_yyyy_xyyz_1[i] * fe_0 + ta_yyyy_xyyyz_0[i] * pa_y[i] - ta_yyyy_xyyyz_1[i] * pc_y[i];

        ta_yyyyy_xyyzz_0[i] = 4.0 * ta_yyy_xyyzz_0[i] * fe_0 - 4.0 * ta_yyy_xyyzz_1[i] * fe_0 + 2.0 * ta_yyyy_xyzz_0[i] * fe_0 -
                              2.0 * ta_yyyy_xyzz_1[i] * fe_0 + ta_yyyy_xyyzz_0[i] * pa_y[i] - ta_yyyy_xyyzz_1[i] * pc_y[i];

        ta_yyyyy_xyzzz_0[i] = 4.0 * ta_yyy_xyzzz_0[i] * fe_0 - 4.0 * ta_yyy_xyzzz_1[i] * fe_0 + ta_yyyy_xzzz_0[i] * fe_0 - ta_yyyy_xzzz_1[i] * fe_0 +
                              ta_yyyy_xyzzz_0[i] * pa_y[i] - ta_yyyy_xyzzz_1[i] * pc_y[i];

        ta_yyyyy_xzzzz_0[i] =
            4.0 * ta_yyy_xzzzz_0[i] * fe_0 - 4.0 * ta_yyy_xzzzz_1[i] * fe_0 + ta_yyyy_xzzzz_0[i] * pa_y[i] - ta_yyyy_xzzzz_1[i] * pc_y[i];

        ta_yyyyy_yyyyy_0[i] = 4.0 * ta_yyy_yyyyy_0[i] * fe_0 - 4.0 * ta_yyy_yyyyy_1[i] * fe_0 + 5.0 * ta_yyyy_yyyy_0[i] * fe_0 -
                              5.0 * ta_yyyy_yyyy_1[i] * fe_0 + ta_yyyy_yyyyy_0[i] * pa_y[i] - ta_yyyy_yyyyy_1[i] * pc_y[i];

        ta_yyyyy_yyyyz_0[i] = 4.0 * ta_yyy_yyyyz_0[i] * fe_0 - 4.0 * ta_yyy_yyyyz_1[i] * fe_0 + 4.0 * ta_yyyy_yyyz_0[i] * fe_0 -
                              4.0 * ta_yyyy_yyyz_1[i] * fe_0 + ta_yyyy_yyyyz_0[i] * pa_y[i] - ta_yyyy_yyyyz_1[i] * pc_y[i];

        ta_yyyyy_yyyzz_0[i] = 4.0 * ta_yyy_yyyzz_0[i] * fe_0 - 4.0 * ta_yyy_yyyzz_1[i] * fe_0 + 3.0 * ta_yyyy_yyzz_0[i] * fe_0 -
                              3.0 * ta_yyyy_yyzz_1[i] * fe_0 + ta_yyyy_yyyzz_0[i] * pa_y[i] - ta_yyyy_yyyzz_1[i] * pc_y[i];

        ta_yyyyy_yyzzz_0[i] = 4.0 * ta_yyy_yyzzz_0[i] * fe_0 - 4.0 * ta_yyy_yyzzz_1[i] * fe_0 + 2.0 * ta_yyyy_yzzz_0[i] * fe_0 -
                              2.0 * ta_yyyy_yzzz_1[i] * fe_0 + ta_yyyy_yyzzz_0[i] * pa_y[i] - ta_yyyy_yyzzz_1[i] * pc_y[i];

        ta_yyyyy_yzzzz_0[i] = 4.0 * ta_yyy_yzzzz_0[i] * fe_0 - 4.0 * ta_yyy_yzzzz_1[i] * fe_0 + ta_yyyy_zzzz_0[i] * fe_0 - ta_yyyy_zzzz_1[i] * fe_0 +
                              ta_yyyy_yzzzz_0[i] * pa_y[i] - ta_yyyy_yzzzz_1[i] * pc_y[i];

        ta_yyyyy_zzzzz_0[i] =
            4.0 * ta_yyy_zzzzz_0[i] * fe_0 - 4.0 * ta_yyy_zzzzz_1[i] * fe_0 + ta_yyyy_zzzzz_0[i] * pa_y[i] - ta_yyyy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 336-357 components of targeted buffer : HH

    auto ta_yyyyz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 336);

    auto ta_yyyyz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 337);

    auto ta_yyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 338);

    auto ta_yyyyz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 339);

    auto ta_yyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 340);

    auto ta_yyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 341);

    auto ta_yyyyz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 342);

    auto ta_yyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 343);

    auto ta_yyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 344);

    auto ta_yyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 345);

    auto ta_yyyyz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 346);

    auto ta_yyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 347);

    auto ta_yyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 348);

    auto ta_yyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 349);

    auto ta_yyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 350);

    auto ta_yyyyz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 351);

    auto ta_yyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 352);

    auto ta_yyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 353);

    auto ta_yyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 354);

    auto ta_yyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 355);

    auto ta_yyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 356);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyyy_xxxxx_0,  \
                             ta_yyyy_xxxxx_1,  \
                             ta_yyyy_xxxxy_0,  \
                             ta_yyyy_xxxxy_1,  \
                             ta_yyyy_xxxy_0,   \
                             ta_yyyy_xxxy_1,   \
                             ta_yyyy_xxxyy_0,  \
                             ta_yyyy_xxxyy_1,  \
                             ta_yyyy_xxxyz_0,  \
                             ta_yyyy_xxxyz_1,  \
                             ta_yyyy_xxyy_0,   \
                             ta_yyyy_xxyy_1,   \
                             ta_yyyy_xxyyy_0,  \
                             ta_yyyy_xxyyy_1,  \
                             ta_yyyy_xxyyz_0,  \
                             ta_yyyy_xxyyz_1,  \
                             ta_yyyy_xxyz_0,   \
                             ta_yyyy_xxyz_1,   \
                             ta_yyyy_xxyzz_0,  \
                             ta_yyyy_xxyzz_1,  \
                             ta_yyyy_xyyy_0,   \
                             ta_yyyy_xyyy_1,   \
                             ta_yyyy_xyyyy_0,  \
                             ta_yyyy_xyyyy_1,  \
                             ta_yyyy_xyyyz_0,  \
                             ta_yyyy_xyyyz_1,  \
                             ta_yyyy_xyyz_0,   \
                             ta_yyyy_xyyz_1,   \
                             ta_yyyy_xyyzz_0,  \
                             ta_yyyy_xyyzz_1,  \
                             ta_yyyy_xyzz_0,   \
                             ta_yyyy_xyzz_1,   \
                             ta_yyyy_xyzzz_0,  \
                             ta_yyyy_xyzzz_1,  \
                             ta_yyyy_yyyy_0,   \
                             ta_yyyy_yyyy_1,   \
                             ta_yyyy_yyyyy_0,  \
                             ta_yyyy_yyyyy_1,  \
                             ta_yyyy_yyyyz_0,  \
                             ta_yyyy_yyyyz_1,  \
                             ta_yyyy_yyyz_0,   \
                             ta_yyyy_yyyz_1,   \
                             ta_yyyy_yyyzz_0,  \
                             ta_yyyy_yyyzz_1,  \
                             ta_yyyy_yyzz_0,   \
                             ta_yyyy_yyzz_1,   \
                             ta_yyyy_yyzzz_0,  \
                             ta_yyyy_yyzzz_1,  \
                             ta_yyyy_yzzz_0,   \
                             ta_yyyy_yzzz_1,   \
                             ta_yyyy_yzzzz_0,  \
                             ta_yyyy_yzzzz_1,  \
                             ta_yyyyz_xxxxx_0, \
                             ta_yyyyz_xxxxy_0, \
                             ta_yyyyz_xxxxz_0, \
                             ta_yyyyz_xxxyy_0, \
                             ta_yyyyz_xxxyz_0, \
                             ta_yyyyz_xxxzz_0, \
                             ta_yyyyz_xxyyy_0, \
                             ta_yyyyz_xxyyz_0, \
                             ta_yyyyz_xxyzz_0, \
                             ta_yyyyz_xxzzz_0, \
                             ta_yyyyz_xyyyy_0, \
                             ta_yyyyz_xyyyz_0, \
                             ta_yyyyz_xyyzz_0, \
                             ta_yyyyz_xyzzz_0, \
                             ta_yyyyz_xzzzz_0, \
                             ta_yyyyz_yyyyy_0, \
                             ta_yyyyz_yyyyz_0, \
                             ta_yyyyz_yyyzz_0, \
                             ta_yyyyz_yyzzz_0, \
                             ta_yyyyz_yzzzz_0, \
                             ta_yyyyz_zzzzz_0, \
                             ta_yyyz_xxxxz_0,  \
                             ta_yyyz_xxxxz_1,  \
                             ta_yyyz_xxxzz_0,  \
                             ta_yyyz_xxxzz_1,  \
                             ta_yyyz_xxzzz_0,  \
                             ta_yyyz_xxzzz_1,  \
                             ta_yyyz_xzzzz_0,  \
                             ta_yyyz_xzzzz_1,  \
                             ta_yyyz_zzzzz_0,  \
                             ta_yyyz_zzzzz_1,  \
                             ta_yyz_xxxxz_0,   \
                             ta_yyz_xxxxz_1,   \
                             ta_yyz_xxxzz_0,   \
                             ta_yyz_xxxzz_1,   \
                             ta_yyz_xxzzz_0,   \
                             ta_yyz_xxzzz_1,   \
                             ta_yyz_xzzzz_0,   \
                             ta_yyz_xzzzz_1,   \
                             ta_yyz_zzzzz_0,   \
                             ta_yyz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyz_xxxxx_0[i] = ta_yyyy_xxxxx_0[i] * pa_z[i] - ta_yyyy_xxxxx_1[i] * pc_z[i];

        ta_yyyyz_xxxxy_0[i] = ta_yyyy_xxxxy_0[i] * pa_z[i] - ta_yyyy_xxxxy_1[i] * pc_z[i];

        ta_yyyyz_xxxxz_0[i] =
            3.0 * ta_yyz_xxxxz_0[i] * fe_0 - 3.0 * ta_yyz_xxxxz_1[i] * fe_0 + ta_yyyz_xxxxz_0[i] * pa_y[i] - ta_yyyz_xxxxz_1[i] * pc_y[i];

        ta_yyyyz_xxxyy_0[i] = ta_yyyy_xxxyy_0[i] * pa_z[i] - ta_yyyy_xxxyy_1[i] * pc_z[i];

        ta_yyyyz_xxxyz_0[i] = ta_yyyy_xxxy_0[i] * fe_0 - ta_yyyy_xxxy_1[i] * fe_0 + ta_yyyy_xxxyz_0[i] * pa_z[i] - ta_yyyy_xxxyz_1[i] * pc_z[i];

        ta_yyyyz_xxxzz_0[i] =
            3.0 * ta_yyz_xxxzz_0[i] * fe_0 - 3.0 * ta_yyz_xxxzz_1[i] * fe_0 + ta_yyyz_xxxzz_0[i] * pa_y[i] - ta_yyyz_xxxzz_1[i] * pc_y[i];

        ta_yyyyz_xxyyy_0[i] = ta_yyyy_xxyyy_0[i] * pa_z[i] - ta_yyyy_xxyyy_1[i] * pc_z[i];

        ta_yyyyz_xxyyz_0[i] = ta_yyyy_xxyy_0[i] * fe_0 - ta_yyyy_xxyy_1[i] * fe_0 + ta_yyyy_xxyyz_0[i] * pa_z[i] - ta_yyyy_xxyyz_1[i] * pc_z[i];

        ta_yyyyz_xxyzz_0[i] =
            2.0 * ta_yyyy_xxyz_0[i] * fe_0 - 2.0 * ta_yyyy_xxyz_1[i] * fe_0 + ta_yyyy_xxyzz_0[i] * pa_z[i] - ta_yyyy_xxyzz_1[i] * pc_z[i];

        ta_yyyyz_xxzzz_0[i] =
            3.0 * ta_yyz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyz_xxzzz_1[i] * fe_0 + ta_yyyz_xxzzz_0[i] * pa_y[i] - ta_yyyz_xxzzz_1[i] * pc_y[i];

        ta_yyyyz_xyyyy_0[i] = ta_yyyy_xyyyy_0[i] * pa_z[i] - ta_yyyy_xyyyy_1[i] * pc_z[i];

        ta_yyyyz_xyyyz_0[i] = ta_yyyy_xyyy_0[i] * fe_0 - ta_yyyy_xyyy_1[i] * fe_0 + ta_yyyy_xyyyz_0[i] * pa_z[i] - ta_yyyy_xyyyz_1[i] * pc_z[i];

        ta_yyyyz_xyyzz_0[i] =
            2.0 * ta_yyyy_xyyz_0[i] * fe_0 - 2.0 * ta_yyyy_xyyz_1[i] * fe_0 + ta_yyyy_xyyzz_0[i] * pa_z[i] - ta_yyyy_xyyzz_1[i] * pc_z[i];

        ta_yyyyz_xyzzz_0[i] =
            3.0 * ta_yyyy_xyzz_0[i] * fe_0 - 3.0 * ta_yyyy_xyzz_1[i] * fe_0 + ta_yyyy_xyzzz_0[i] * pa_z[i] - ta_yyyy_xyzzz_1[i] * pc_z[i];

        ta_yyyyz_xzzzz_0[i] =
            3.0 * ta_yyz_xzzzz_0[i] * fe_0 - 3.0 * ta_yyz_xzzzz_1[i] * fe_0 + ta_yyyz_xzzzz_0[i] * pa_y[i] - ta_yyyz_xzzzz_1[i] * pc_y[i];

        ta_yyyyz_yyyyy_0[i] = ta_yyyy_yyyyy_0[i] * pa_z[i] - ta_yyyy_yyyyy_1[i] * pc_z[i];

        ta_yyyyz_yyyyz_0[i] = ta_yyyy_yyyy_0[i] * fe_0 - ta_yyyy_yyyy_1[i] * fe_0 + ta_yyyy_yyyyz_0[i] * pa_z[i] - ta_yyyy_yyyyz_1[i] * pc_z[i];

        ta_yyyyz_yyyzz_0[i] =
            2.0 * ta_yyyy_yyyz_0[i] * fe_0 - 2.0 * ta_yyyy_yyyz_1[i] * fe_0 + ta_yyyy_yyyzz_0[i] * pa_z[i] - ta_yyyy_yyyzz_1[i] * pc_z[i];

        ta_yyyyz_yyzzz_0[i] =
            3.0 * ta_yyyy_yyzz_0[i] * fe_0 - 3.0 * ta_yyyy_yyzz_1[i] * fe_0 + ta_yyyy_yyzzz_0[i] * pa_z[i] - ta_yyyy_yyzzz_1[i] * pc_z[i];

        ta_yyyyz_yzzzz_0[i] =
            4.0 * ta_yyyy_yzzz_0[i] * fe_0 - 4.0 * ta_yyyy_yzzz_1[i] * fe_0 + ta_yyyy_yzzzz_0[i] * pa_z[i] - ta_yyyy_yzzzz_1[i] * pc_z[i];

        ta_yyyyz_zzzzz_0[i] =
            3.0 * ta_yyz_zzzzz_0[i] * fe_0 - 3.0 * ta_yyz_zzzzz_1[i] * fe_0 + ta_yyyz_zzzzz_0[i] * pa_y[i] - ta_yyyz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 357-378 components of targeted buffer : HH

    auto ta_yyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 357);

    auto ta_yyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 358);

    auto ta_yyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 359);

    auto ta_yyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 360);

    auto ta_yyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 361);

    auto ta_yyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 362);

    auto ta_yyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 363);

    auto ta_yyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 364);

    auto ta_yyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 365);

    auto ta_yyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 366);

    auto ta_yyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 367);

    auto ta_yyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 368);

    auto ta_yyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 369);

    auto ta_yyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 370);

    auto ta_yyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 371);

    auto ta_yyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 372);

    auto ta_yyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 373);

    auto ta_yyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 374);

    auto ta_yyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 375);

    auto ta_yyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 376);

    auto ta_yyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 377);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyy_xxxxy_0,   \
                             ta_yyy_xxxxy_1,   \
                             ta_yyy_xxxyy_0,   \
                             ta_yyy_xxxyy_1,   \
                             ta_yyy_xxyyy_0,   \
                             ta_yyy_xxyyy_1,   \
                             ta_yyy_xyyyy_0,   \
                             ta_yyy_xyyyy_1,   \
                             ta_yyy_yyyyy_0,   \
                             ta_yyy_yyyyy_1,   \
                             ta_yyyz_xxxxy_0,  \
                             ta_yyyz_xxxxy_1,  \
                             ta_yyyz_xxxyy_0,  \
                             ta_yyyz_xxxyy_1,  \
                             ta_yyyz_xxyyy_0,  \
                             ta_yyyz_xxyyy_1,  \
                             ta_yyyz_xyyyy_0,  \
                             ta_yyyz_xyyyy_1,  \
                             ta_yyyz_yyyyy_0,  \
                             ta_yyyz_yyyyy_1,  \
                             ta_yyyzz_xxxxx_0, \
                             ta_yyyzz_xxxxy_0, \
                             ta_yyyzz_xxxxz_0, \
                             ta_yyyzz_xxxyy_0, \
                             ta_yyyzz_xxxyz_0, \
                             ta_yyyzz_xxxzz_0, \
                             ta_yyyzz_xxyyy_0, \
                             ta_yyyzz_xxyyz_0, \
                             ta_yyyzz_xxyzz_0, \
                             ta_yyyzz_xxzzz_0, \
                             ta_yyyzz_xyyyy_0, \
                             ta_yyyzz_xyyyz_0, \
                             ta_yyyzz_xyyzz_0, \
                             ta_yyyzz_xyzzz_0, \
                             ta_yyyzz_xzzzz_0, \
                             ta_yyyzz_yyyyy_0, \
                             ta_yyyzz_yyyyz_0, \
                             ta_yyyzz_yyyzz_0, \
                             ta_yyyzz_yyzzz_0, \
                             ta_yyyzz_yzzzz_0, \
                             ta_yyyzz_zzzzz_0, \
                             ta_yyzz_xxxxx_0,  \
                             ta_yyzz_xxxxx_1,  \
                             ta_yyzz_xxxxz_0,  \
                             ta_yyzz_xxxxz_1,  \
                             ta_yyzz_xxxyz_0,  \
                             ta_yyzz_xxxyz_1,  \
                             ta_yyzz_xxxz_0,   \
                             ta_yyzz_xxxz_1,   \
                             ta_yyzz_xxxzz_0,  \
                             ta_yyzz_xxxzz_1,  \
                             ta_yyzz_xxyyz_0,  \
                             ta_yyzz_xxyyz_1,  \
                             ta_yyzz_xxyz_0,   \
                             ta_yyzz_xxyz_1,   \
                             ta_yyzz_xxyzz_0,  \
                             ta_yyzz_xxyzz_1,  \
                             ta_yyzz_xxzz_0,   \
                             ta_yyzz_xxzz_1,   \
                             ta_yyzz_xxzzz_0,  \
                             ta_yyzz_xxzzz_1,  \
                             ta_yyzz_xyyyz_0,  \
                             ta_yyzz_xyyyz_1,  \
                             ta_yyzz_xyyz_0,   \
                             ta_yyzz_xyyz_1,   \
                             ta_yyzz_xyyzz_0,  \
                             ta_yyzz_xyyzz_1,  \
                             ta_yyzz_xyzz_0,   \
                             ta_yyzz_xyzz_1,   \
                             ta_yyzz_xyzzz_0,  \
                             ta_yyzz_xyzzz_1,  \
                             ta_yyzz_xzzz_0,   \
                             ta_yyzz_xzzz_1,   \
                             ta_yyzz_xzzzz_0,  \
                             ta_yyzz_xzzzz_1,  \
                             ta_yyzz_yyyyz_0,  \
                             ta_yyzz_yyyyz_1,  \
                             ta_yyzz_yyyz_0,   \
                             ta_yyzz_yyyz_1,   \
                             ta_yyzz_yyyzz_0,  \
                             ta_yyzz_yyyzz_1,  \
                             ta_yyzz_yyzz_0,   \
                             ta_yyzz_yyzz_1,   \
                             ta_yyzz_yyzzz_0,  \
                             ta_yyzz_yyzzz_1,  \
                             ta_yyzz_yzzz_0,   \
                             ta_yyzz_yzzz_1,   \
                             ta_yyzz_yzzzz_0,  \
                             ta_yyzz_yzzzz_1,  \
                             ta_yyzz_zzzz_0,   \
                             ta_yyzz_zzzz_1,   \
                             ta_yyzz_zzzzz_0,  \
                             ta_yyzz_zzzzz_1,  \
                             ta_yzz_xxxxx_0,   \
                             ta_yzz_xxxxx_1,   \
                             ta_yzz_xxxxz_0,   \
                             ta_yzz_xxxxz_1,   \
                             ta_yzz_xxxyz_0,   \
                             ta_yzz_xxxyz_1,   \
                             ta_yzz_xxxzz_0,   \
                             ta_yzz_xxxzz_1,   \
                             ta_yzz_xxyyz_0,   \
                             ta_yzz_xxyyz_1,   \
                             ta_yzz_xxyzz_0,   \
                             ta_yzz_xxyzz_1,   \
                             ta_yzz_xxzzz_0,   \
                             ta_yzz_xxzzz_1,   \
                             ta_yzz_xyyyz_0,   \
                             ta_yzz_xyyyz_1,   \
                             ta_yzz_xyyzz_0,   \
                             ta_yzz_xyyzz_1,   \
                             ta_yzz_xyzzz_0,   \
                             ta_yzz_xyzzz_1,   \
                             ta_yzz_xzzzz_0,   \
                             ta_yzz_xzzzz_1,   \
                             ta_yzz_yyyyz_0,   \
                             ta_yzz_yyyyz_1,   \
                             ta_yzz_yyyzz_0,   \
                             ta_yzz_yyyzz_1,   \
                             ta_yzz_yyzzz_0,   \
                             ta_yzz_yyzzz_1,   \
                             ta_yzz_yzzzz_0,   \
                             ta_yzz_yzzzz_1,   \
                             ta_yzz_zzzzz_0,   \
                             ta_yzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzz_xxxxx_0[i] =
            2.0 * ta_yzz_xxxxx_0[i] * fe_0 - 2.0 * ta_yzz_xxxxx_1[i] * fe_0 + ta_yyzz_xxxxx_0[i] * pa_y[i] - ta_yyzz_xxxxx_1[i] * pc_y[i];

        ta_yyyzz_xxxxy_0[i] = ta_yyy_xxxxy_0[i] * fe_0 - ta_yyy_xxxxy_1[i] * fe_0 + ta_yyyz_xxxxy_0[i] * pa_z[i] - ta_yyyz_xxxxy_1[i] * pc_z[i];

        ta_yyyzz_xxxxz_0[i] =
            2.0 * ta_yzz_xxxxz_0[i] * fe_0 - 2.0 * ta_yzz_xxxxz_1[i] * fe_0 + ta_yyzz_xxxxz_0[i] * pa_y[i] - ta_yyzz_xxxxz_1[i] * pc_y[i];

        ta_yyyzz_xxxyy_0[i] = ta_yyy_xxxyy_0[i] * fe_0 - ta_yyy_xxxyy_1[i] * fe_0 + ta_yyyz_xxxyy_0[i] * pa_z[i] - ta_yyyz_xxxyy_1[i] * pc_z[i];

        ta_yyyzz_xxxyz_0[i] = 2.0 * ta_yzz_xxxyz_0[i] * fe_0 - 2.0 * ta_yzz_xxxyz_1[i] * fe_0 + ta_yyzz_xxxz_0[i] * fe_0 - ta_yyzz_xxxz_1[i] * fe_0 +
                              ta_yyzz_xxxyz_0[i] * pa_y[i] - ta_yyzz_xxxyz_1[i] * pc_y[i];

        ta_yyyzz_xxxzz_0[i] =
            2.0 * ta_yzz_xxxzz_0[i] * fe_0 - 2.0 * ta_yzz_xxxzz_1[i] * fe_0 + ta_yyzz_xxxzz_0[i] * pa_y[i] - ta_yyzz_xxxzz_1[i] * pc_y[i];

        ta_yyyzz_xxyyy_0[i] = ta_yyy_xxyyy_0[i] * fe_0 - ta_yyy_xxyyy_1[i] * fe_0 + ta_yyyz_xxyyy_0[i] * pa_z[i] - ta_yyyz_xxyyy_1[i] * pc_z[i];

        ta_yyyzz_xxyyz_0[i] = 2.0 * ta_yzz_xxyyz_0[i] * fe_0 - 2.0 * ta_yzz_xxyyz_1[i] * fe_0 + 2.0 * ta_yyzz_xxyz_0[i] * fe_0 -
                              2.0 * ta_yyzz_xxyz_1[i] * fe_0 + ta_yyzz_xxyyz_0[i] * pa_y[i] - ta_yyzz_xxyyz_1[i] * pc_y[i];

        ta_yyyzz_xxyzz_0[i] = 2.0 * ta_yzz_xxyzz_0[i] * fe_0 - 2.0 * ta_yzz_xxyzz_1[i] * fe_0 + ta_yyzz_xxzz_0[i] * fe_0 - ta_yyzz_xxzz_1[i] * fe_0 +
                              ta_yyzz_xxyzz_0[i] * pa_y[i] - ta_yyzz_xxyzz_1[i] * pc_y[i];

        ta_yyyzz_xxzzz_0[i] =
            2.0 * ta_yzz_xxzzz_0[i] * fe_0 - 2.0 * ta_yzz_xxzzz_1[i] * fe_0 + ta_yyzz_xxzzz_0[i] * pa_y[i] - ta_yyzz_xxzzz_1[i] * pc_y[i];

        ta_yyyzz_xyyyy_0[i] = ta_yyy_xyyyy_0[i] * fe_0 - ta_yyy_xyyyy_1[i] * fe_0 + ta_yyyz_xyyyy_0[i] * pa_z[i] - ta_yyyz_xyyyy_1[i] * pc_z[i];

        ta_yyyzz_xyyyz_0[i] = 2.0 * ta_yzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yzz_xyyyz_1[i] * fe_0 + 3.0 * ta_yyzz_xyyz_0[i] * fe_0 -
                              3.0 * ta_yyzz_xyyz_1[i] * fe_0 + ta_yyzz_xyyyz_0[i] * pa_y[i] - ta_yyzz_xyyyz_1[i] * pc_y[i];

        ta_yyyzz_xyyzz_0[i] = 2.0 * ta_yzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yzz_xyyzz_1[i] * fe_0 + 2.0 * ta_yyzz_xyzz_0[i] * fe_0 -
                              2.0 * ta_yyzz_xyzz_1[i] * fe_0 + ta_yyzz_xyyzz_0[i] * pa_y[i] - ta_yyzz_xyyzz_1[i] * pc_y[i];

        ta_yyyzz_xyzzz_0[i] = 2.0 * ta_yzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzzz_1[i] * fe_0 + ta_yyzz_xzzz_0[i] * fe_0 - ta_yyzz_xzzz_1[i] * fe_0 +
                              ta_yyzz_xyzzz_0[i] * pa_y[i] - ta_yyzz_xyzzz_1[i] * pc_y[i];

        ta_yyyzz_xzzzz_0[i] =
            2.0 * ta_yzz_xzzzz_0[i] * fe_0 - 2.0 * ta_yzz_xzzzz_1[i] * fe_0 + ta_yyzz_xzzzz_0[i] * pa_y[i] - ta_yyzz_xzzzz_1[i] * pc_y[i];

        ta_yyyzz_yyyyy_0[i] = ta_yyy_yyyyy_0[i] * fe_0 - ta_yyy_yyyyy_1[i] * fe_0 + ta_yyyz_yyyyy_0[i] * pa_z[i] - ta_yyyz_yyyyy_1[i] * pc_z[i];

        ta_yyyzz_yyyyz_0[i] = 2.0 * ta_yzz_yyyyz_0[i] * fe_0 - 2.0 * ta_yzz_yyyyz_1[i] * fe_0 + 4.0 * ta_yyzz_yyyz_0[i] * fe_0 -
                              4.0 * ta_yyzz_yyyz_1[i] * fe_0 + ta_yyzz_yyyyz_0[i] * pa_y[i] - ta_yyzz_yyyyz_1[i] * pc_y[i];

        ta_yyyzz_yyyzz_0[i] = 2.0 * ta_yzz_yyyzz_0[i] * fe_0 - 2.0 * ta_yzz_yyyzz_1[i] * fe_0 + 3.0 * ta_yyzz_yyzz_0[i] * fe_0 -
                              3.0 * ta_yyzz_yyzz_1[i] * fe_0 + ta_yyzz_yyyzz_0[i] * pa_y[i] - ta_yyzz_yyyzz_1[i] * pc_y[i];

        ta_yyyzz_yyzzz_0[i] = 2.0 * ta_yzz_yyzzz_0[i] * fe_0 - 2.0 * ta_yzz_yyzzz_1[i] * fe_0 + 2.0 * ta_yyzz_yzzz_0[i] * fe_0 -
                              2.0 * ta_yyzz_yzzz_1[i] * fe_0 + ta_yyzz_yyzzz_0[i] * pa_y[i] - ta_yyzz_yyzzz_1[i] * pc_y[i];

        ta_yyyzz_yzzzz_0[i] = 2.0 * ta_yzz_yzzzz_0[i] * fe_0 - 2.0 * ta_yzz_yzzzz_1[i] * fe_0 + ta_yyzz_zzzz_0[i] * fe_0 - ta_yyzz_zzzz_1[i] * fe_0 +
                              ta_yyzz_yzzzz_0[i] * pa_y[i] - ta_yyzz_yzzzz_1[i] * pc_y[i];

        ta_yyyzz_zzzzz_0[i] =
            2.0 * ta_yzz_zzzzz_0[i] * fe_0 - 2.0 * ta_yzz_zzzzz_1[i] * fe_0 + ta_yyzz_zzzzz_0[i] * pa_y[i] - ta_yyzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 378-399 components of targeted buffer : HH

    auto ta_yyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 378);

    auto ta_yyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 379);

    auto ta_yyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 380);

    auto ta_yyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 381);

    auto ta_yyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 382);

    auto ta_yyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 383);

    auto ta_yyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 384);

    auto ta_yyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 385);

    auto ta_yyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 386);

    auto ta_yyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 387);

    auto ta_yyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 388);

    auto ta_yyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 389);

    auto ta_yyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 390);

    auto ta_yyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 391);

    auto ta_yyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 392);

    auto ta_yyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 393);

    auto ta_yyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 394);

    auto ta_yyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 395);

    auto ta_yyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 396);

    auto ta_yyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 397);

    auto ta_yyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 398);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta_yyz_xxxxy_0,   \
                             ta_yyz_xxxxy_1,   \
                             ta_yyz_xxxyy_0,   \
                             ta_yyz_xxxyy_1,   \
                             ta_yyz_xxyyy_0,   \
                             ta_yyz_xxyyy_1,   \
                             ta_yyz_xyyyy_0,   \
                             ta_yyz_xyyyy_1,   \
                             ta_yyz_yyyyy_0,   \
                             ta_yyz_yyyyy_1,   \
                             ta_yyzz_xxxxy_0,  \
                             ta_yyzz_xxxxy_1,  \
                             ta_yyzz_xxxyy_0,  \
                             ta_yyzz_xxxyy_1,  \
                             ta_yyzz_xxyyy_0,  \
                             ta_yyzz_xxyyy_1,  \
                             ta_yyzz_xyyyy_0,  \
                             ta_yyzz_xyyyy_1,  \
                             ta_yyzz_yyyyy_0,  \
                             ta_yyzz_yyyyy_1,  \
                             ta_yyzzz_xxxxx_0, \
                             ta_yyzzz_xxxxy_0, \
                             ta_yyzzz_xxxxz_0, \
                             ta_yyzzz_xxxyy_0, \
                             ta_yyzzz_xxxyz_0, \
                             ta_yyzzz_xxxzz_0, \
                             ta_yyzzz_xxyyy_0, \
                             ta_yyzzz_xxyyz_0, \
                             ta_yyzzz_xxyzz_0, \
                             ta_yyzzz_xxzzz_0, \
                             ta_yyzzz_xyyyy_0, \
                             ta_yyzzz_xyyyz_0, \
                             ta_yyzzz_xyyzz_0, \
                             ta_yyzzz_xyzzz_0, \
                             ta_yyzzz_xzzzz_0, \
                             ta_yyzzz_yyyyy_0, \
                             ta_yyzzz_yyyyz_0, \
                             ta_yyzzz_yyyzz_0, \
                             ta_yyzzz_yyzzz_0, \
                             ta_yyzzz_yzzzz_0, \
                             ta_yyzzz_zzzzz_0, \
                             ta_yzzz_xxxxx_0,  \
                             ta_yzzz_xxxxx_1,  \
                             ta_yzzz_xxxxz_0,  \
                             ta_yzzz_xxxxz_1,  \
                             ta_yzzz_xxxyz_0,  \
                             ta_yzzz_xxxyz_1,  \
                             ta_yzzz_xxxz_0,   \
                             ta_yzzz_xxxz_1,   \
                             ta_yzzz_xxxzz_0,  \
                             ta_yzzz_xxxzz_1,  \
                             ta_yzzz_xxyyz_0,  \
                             ta_yzzz_xxyyz_1,  \
                             ta_yzzz_xxyz_0,   \
                             ta_yzzz_xxyz_1,   \
                             ta_yzzz_xxyzz_0,  \
                             ta_yzzz_xxyzz_1,  \
                             ta_yzzz_xxzz_0,   \
                             ta_yzzz_xxzz_1,   \
                             ta_yzzz_xxzzz_0,  \
                             ta_yzzz_xxzzz_1,  \
                             ta_yzzz_xyyyz_0,  \
                             ta_yzzz_xyyyz_1,  \
                             ta_yzzz_xyyz_0,   \
                             ta_yzzz_xyyz_1,   \
                             ta_yzzz_xyyzz_0,  \
                             ta_yzzz_xyyzz_1,  \
                             ta_yzzz_xyzz_0,   \
                             ta_yzzz_xyzz_1,   \
                             ta_yzzz_xyzzz_0,  \
                             ta_yzzz_xyzzz_1,  \
                             ta_yzzz_xzzz_0,   \
                             ta_yzzz_xzzz_1,   \
                             ta_yzzz_xzzzz_0,  \
                             ta_yzzz_xzzzz_1,  \
                             ta_yzzz_yyyyz_0,  \
                             ta_yzzz_yyyyz_1,  \
                             ta_yzzz_yyyz_0,   \
                             ta_yzzz_yyyz_1,   \
                             ta_yzzz_yyyzz_0,  \
                             ta_yzzz_yyyzz_1,  \
                             ta_yzzz_yyzz_0,   \
                             ta_yzzz_yyzz_1,   \
                             ta_yzzz_yyzzz_0,  \
                             ta_yzzz_yyzzz_1,  \
                             ta_yzzz_yzzz_0,   \
                             ta_yzzz_yzzz_1,   \
                             ta_yzzz_yzzzz_0,  \
                             ta_yzzz_yzzzz_1,  \
                             ta_yzzz_zzzz_0,   \
                             ta_yzzz_zzzz_1,   \
                             ta_yzzz_zzzzz_0,  \
                             ta_yzzz_zzzzz_1,  \
                             ta_zzz_xxxxx_0,   \
                             ta_zzz_xxxxx_1,   \
                             ta_zzz_xxxxz_0,   \
                             ta_zzz_xxxxz_1,   \
                             ta_zzz_xxxyz_0,   \
                             ta_zzz_xxxyz_1,   \
                             ta_zzz_xxxzz_0,   \
                             ta_zzz_xxxzz_1,   \
                             ta_zzz_xxyyz_0,   \
                             ta_zzz_xxyyz_1,   \
                             ta_zzz_xxyzz_0,   \
                             ta_zzz_xxyzz_1,   \
                             ta_zzz_xxzzz_0,   \
                             ta_zzz_xxzzz_1,   \
                             ta_zzz_xyyyz_0,   \
                             ta_zzz_xyyyz_1,   \
                             ta_zzz_xyyzz_0,   \
                             ta_zzz_xyyzz_1,   \
                             ta_zzz_xyzzz_0,   \
                             ta_zzz_xyzzz_1,   \
                             ta_zzz_xzzzz_0,   \
                             ta_zzz_xzzzz_1,   \
                             ta_zzz_yyyyz_0,   \
                             ta_zzz_yyyyz_1,   \
                             ta_zzz_yyyzz_0,   \
                             ta_zzz_yyyzz_1,   \
                             ta_zzz_yyzzz_0,   \
                             ta_zzz_yyzzz_1,   \
                             ta_zzz_yzzzz_0,   \
                             ta_zzz_yzzzz_1,   \
                             ta_zzz_zzzzz_0,   \
                             ta_zzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzz_xxxxx_0[i] = ta_zzz_xxxxx_0[i] * fe_0 - ta_zzz_xxxxx_1[i] * fe_0 + ta_yzzz_xxxxx_0[i] * pa_y[i] - ta_yzzz_xxxxx_1[i] * pc_y[i];

        ta_yyzzz_xxxxy_0[i] =
            2.0 * ta_yyz_xxxxy_0[i] * fe_0 - 2.0 * ta_yyz_xxxxy_1[i] * fe_0 + ta_yyzz_xxxxy_0[i] * pa_z[i] - ta_yyzz_xxxxy_1[i] * pc_z[i];

        ta_yyzzz_xxxxz_0[i] = ta_zzz_xxxxz_0[i] * fe_0 - ta_zzz_xxxxz_1[i] * fe_0 + ta_yzzz_xxxxz_0[i] * pa_y[i] - ta_yzzz_xxxxz_1[i] * pc_y[i];

        ta_yyzzz_xxxyy_0[i] =
            2.0 * ta_yyz_xxxyy_0[i] * fe_0 - 2.0 * ta_yyz_xxxyy_1[i] * fe_0 + ta_yyzz_xxxyy_0[i] * pa_z[i] - ta_yyzz_xxxyy_1[i] * pc_z[i];

        ta_yyzzz_xxxyz_0[i] = ta_zzz_xxxyz_0[i] * fe_0 - ta_zzz_xxxyz_1[i] * fe_0 + ta_yzzz_xxxz_0[i] * fe_0 - ta_yzzz_xxxz_1[i] * fe_0 +
                              ta_yzzz_xxxyz_0[i] * pa_y[i] - ta_yzzz_xxxyz_1[i] * pc_y[i];

        ta_yyzzz_xxxzz_0[i] = ta_zzz_xxxzz_0[i] * fe_0 - ta_zzz_xxxzz_1[i] * fe_0 + ta_yzzz_xxxzz_0[i] * pa_y[i] - ta_yzzz_xxxzz_1[i] * pc_y[i];

        ta_yyzzz_xxyyy_0[i] =
            2.0 * ta_yyz_xxyyy_0[i] * fe_0 - 2.0 * ta_yyz_xxyyy_1[i] * fe_0 + ta_yyzz_xxyyy_0[i] * pa_z[i] - ta_yyzz_xxyyy_1[i] * pc_z[i];

        ta_yyzzz_xxyyz_0[i] = ta_zzz_xxyyz_0[i] * fe_0 - ta_zzz_xxyyz_1[i] * fe_0 + 2.0 * ta_yzzz_xxyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyz_1[i] * fe_0 +
                              ta_yzzz_xxyyz_0[i] * pa_y[i] - ta_yzzz_xxyyz_1[i] * pc_y[i];

        ta_yyzzz_xxyzz_0[i] = ta_zzz_xxyzz_0[i] * fe_0 - ta_zzz_xxyzz_1[i] * fe_0 + ta_yzzz_xxzz_0[i] * fe_0 - ta_yzzz_xxzz_1[i] * fe_0 +
                              ta_yzzz_xxyzz_0[i] * pa_y[i] - ta_yzzz_xxyzz_1[i] * pc_y[i];

        ta_yyzzz_xxzzz_0[i] = ta_zzz_xxzzz_0[i] * fe_0 - ta_zzz_xxzzz_1[i] * fe_0 + ta_yzzz_xxzzz_0[i] * pa_y[i] - ta_yzzz_xxzzz_1[i] * pc_y[i];

        ta_yyzzz_xyyyy_0[i] =
            2.0 * ta_yyz_xyyyy_0[i] * fe_0 - 2.0 * ta_yyz_xyyyy_1[i] * fe_0 + ta_yyzz_xyyyy_0[i] * pa_z[i] - ta_yyzz_xyyyy_1[i] * pc_z[i];

        ta_yyzzz_xyyyz_0[i] = ta_zzz_xyyyz_0[i] * fe_0 - ta_zzz_xyyyz_1[i] * fe_0 + 3.0 * ta_yzzz_xyyz_0[i] * fe_0 - 3.0 * ta_yzzz_xyyz_1[i] * fe_0 +
                              ta_yzzz_xyyyz_0[i] * pa_y[i] - ta_yzzz_xyyyz_1[i] * pc_y[i];

        ta_yyzzz_xyyzz_0[i] = ta_zzz_xyyzz_0[i] * fe_0 - ta_zzz_xyyzz_1[i] * fe_0 + 2.0 * ta_yzzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyzz_1[i] * fe_0 +
                              ta_yzzz_xyyzz_0[i] * pa_y[i] - ta_yzzz_xyyzz_1[i] * pc_y[i];

        ta_yyzzz_xyzzz_0[i] = ta_zzz_xyzzz_0[i] * fe_0 - ta_zzz_xyzzz_1[i] * fe_0 + ta_yzzz_xzzz_0[i] * fe_0 - ta_yzzz_xzzz_1[i] * fe_0 +
                              ta_yzzz_xyzzz_0[i] * pa_y[i] - ta_yzzz_xyzzz_1[i] * pc_y[i];

        ta_yyzzz_xzzzz_0[i] = ta_zzz_xzzzz_0[i] * fe_0 - ta_zzz_xzzzz_1[i] * fe_0 + ta_yzzz_xzzzz_0[i] * pa_y[i] - ta_yzzz_xzzzz_1[i] * pc_y[i];

        ta_yyzzz_yyyyy_0[i] =
            2.0 * ta_yyz_yyyyy_0[i] * fe_0 - 2.0 * ta_yyz_yyyyy_1[i] * fe_0 + ta_yyzz_yyyyy_0[i] * pa_z[i] - ta_yyzz_yyyyy_1[i] * pc_z[i];

        ta_yyzzz_yyyyz_0[i] = ta_zzz_yyyyz_0[i] * fe_0 - ta_zzz_yyyyz_1[i] * fe_0 + 4.0 * ta_yzzz_yyyz_0[i] * fe_0 - 4.0 * ta_yzzz_yyyz_1[i] * fe_0 +
                              ta_yzzz_yyyyz_0[i] * pa_y[i] - ta_yzzz_yyyyz_1[i] * pc_y[i];

        ta_yyzzz_yyyzz_0[i] = ta_zzz_yyyzz_0[i] * fe_0 - ta_zzz_yyyzz_1[i] * fe_0 + 3.0 * ta_yzzz_yyzz_0[i] * fe_0 - 3.0 * ta_yzzz_yyzz_1[i] * fe_0 +
                              ta_yzzz_yyyzz_0[i] * pa_y[i] - ta_yzzz_yyyzz_1[i] * pc_y[i];

        ta_yyzzz_yyzzz_0[i] = ta_zzz_yyzzz_0[i] * fe_0 - ta_zzz_yyzzz_1[i] * fe_0 + 2.0 * ta_yzzz_yzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yzzz_1[i] * fe_0 +
                              ta_yzzz_yyzzz_0[i] * pa_y[i] - ta_yzzz_yyzzz_1[i] * pc_y[i];

        ta_yyzzz_yzzzz_0[i] = ta_zzz_yzzzz_0[i] * fe_0 - ta_zzz_yzzzz_1[i] * fe_0 + ta_yzzz_zzzz_0[i] * fe_0 - ta_yzzz_zzzz_1[i] * fe_0 +
                              ta_yzzz_yzzzz_0[i] * pa_y[i] - ta_yzzz_yzzzz_1[i] * pc_y[i];

        ta_yyzzz_zzzzz_0[i] = ta_zzz_zzzzz_0[i] * fe_0 - ta_zzz_zzzzz_1[i] * fe_0 + ta_yzzz_zzzzz_0[i] * pa_y[i] - ta_yzzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 399-420 components of targeted buffer : HH

    auto ta_yzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 399);

    auto ta_yzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 400);

    auto ta_yzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 401);

    auto ta_yzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 402);

    auto ta_yzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 403);

    auto ta_yzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 404);

    auto ta_yzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 405);

    auto ta_yzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 406);

    auto ta_yzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 407);

    auto ta_yzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 408);

    auto ta_yzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 409);

    auto ta_yzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 410);

    auto ta_yzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 411);

    auto ta_yzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 412);

    auto ta_yzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 413);

    auto ta_yzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 414);

    auto ta_yzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 415);

    auto ta_yzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 416);

    auto ta_yzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 417);

    auto ta_yzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 418);

    auto ta_yzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 419);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta_yzzzz_xxxxx_0, \
                             ta_yzzzz_xxxxy_0, \
                             ta_yzzzz_xxxxz_0, \
                             ta_yzzzz_xxxyy_0, \
                             ta_yzzzz_xxxyz_0, \
                             ta_yzzzz_xxxzz_0, \
                             ta_yzzzz_xxyyy_0, \
                             ta_yzzzz_xxyyz_0, \
                             ta_yzzzz_xxyzz_0, \
                             ta_yzzzz_xxzzz_0, \
                             ta_yzzzz_xyyyy_0, \
                             ta_yzzzz_xyyyz_0, \
                             ta_yzzzz_xyyzz_0, \
                             ta_yzzzz_xyzzz_0, \
                             ta_yzzzz_xzzzz_0, \
                             ta_yzzzz_yyyyy_0, \
                             ta_yzzzz_yyyyz_0, \
                             ta_yzzzz_yyyzz_0, \
                             ta_yzzzz_yyzzz_0, \
                             ta_yzzzz_yzzzz_0, \
                             ta_yzzzz_zzzzz_0, \
                             ta_zzzz_xxxx_0,   \
                             ta_zzzz_xxxx_1,   \
                             ta_zzzz_xxxxx_0,  \
                             ta_zzzz_xxxxx_1,  \
                             ta_zzzz_xxxxy_0,  \
                             ta_zzzz_xxxxy_1,  \
                             ta_zzzz_xxxxz_0,  \
                             ta_zzzz_xxxxz_1,  \
                             ta_zzzz_xxxy_0,   \
                             ta_zzzz_xxxy_1,   \
                             ta_zzzz_xxxyy_0,  \
                             ta_zzzz_xxxyy_1,  \
                             ta_zzzz_xxxyz_0,  \
                             ta_zzzz_xxxyz_1,  \
                             ta_zzzz_xxxz_0,   \
                             ta_zzzz_xxxz_1,   \
                             ta_zzzz_xxxzz_0,  \
                             ta_zzzz_xxxzz_1,  \
                             ta_zzzz_xxyy_0,   \
                             ta_zzzz_xxyy_1,   \
                             ta_zzzz_xxyyy_0,  \
                             ta_zzzz_xxyyy_1,  \
                             ta_zzzz_xxyyz_0,  \
                             ta_zzzz_xxyyz_1,  \
                             ta_zzzz_xxyz_0,   \
                             ta_zzzz_xxyz_1,   \
                             ta_zzzz_xxyzz_0,  \
                             ta_zzzz_xxyzz_1,  \
                             ta_zzzz_xxzz_0,   \
                             ta_zzzz_xxzz_1,   \
                             ta_zzzz_xxzzz_0,  \
                             ta_zzzz_xxzzz_1,  \
                             ta_zzzz_xyyy_0,   \
                             ta_zzzz_xyyy_1,   \
                             ta_zzzz_xyyyy_0,  \
                             ta_zzzz_xyyyy_1,  \
                             ta_zzzz_xyyyz_0,  \
                             ta_zzzz_xyyyz_1,  \
                             ta_zzzz_xyyz_0,   \
                             ta_zzzz_xyyz_1,   \
                             ta_zzzz_xyyzz_0,  \
                             ta_zzzz_xyyzz_1,  \
                             ta_zzzz_xyzz_0,   \
                             ta_zzzz_xyzz_1,   \
                             ta_zzzz_xyzzz_0,  \
                             ta_zzzz_xyzzz_1,  \
                             ta_zzzz_xzzz_0,   \
                             ta_zzzz_xzzz_1,   \
                             ta_zzzz_xzzzz_0,  \
                             ta_zzzz_xzzzz_1,  \
                             ta_zzzz_yyyy_0,   \
                             ta_zzzz_yyyy_1,   \
                             ta_zzzz_yyyyy_0,  \
                             ta_zzzz_yyyyy_1,  \
                             ta_zzzz_yyyyz_0,  \
                             ta_zzzz_yyyyz_1,  \
                             ta_zzzz_yyyz_0,   \
                             ta_zzzz_yyyz_1,   \
                             ta_zzzz_yyyzz_0,  \
                             ta_zzzz_yyyzz_1,  \
                             ta_zzzz_yyzz_0,   \
                             ta_zzzz_yyzz_1,   \
                             ta_zzzz_yyzzz_0,  \
                             ta_zzzz_yyzzz_1,  \
                             ta_zzzz_yzzz_0,   \
                             ta_zzzz_yzzz_1,   \
                             ta_zzzz_yzzzz_0,  \
                             ta_zzzz_yzzzz_1,  \
                             ta_zzzz_zzzz_0,   \
                             ta_zzzz_zzzz_1,   \
                             ta_zzzz_zzzzz_0,  \
                             ta_zzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzz_xxxxx_0[i] = ta_zzzz_xxxxx_0[i] * pa_y[i] - ta_zzzz_xxxxx_1[i] * pc_y[i];

        ta_yzzzz_xxxxy_0[i] = ta_zzzz_xxxx_0[i] * fe_0 - ta_zzzz_xxxx_1[i] * fe_0 + ta_zzzz_xxxxy_0[i] * pa_y[i] - ta_zzzz_xxxxy_1[i] * pc_y[i];

        ta_yzzzz_xxxxz_0[i] = ta_zzzz_xxxxz_0[i] * pa_y[i] - ta_zzzz_xxxxz_1[i] * pc_y[i];

        ta_yzzzz_xxxyy_0[i] =
            2.0 * ta_zzzz_xxxy_0[i] * fe_0 - 2.0 * ta_zzzz_xxxy_1[i] * fe_0 + ta_zzzz_xxxyy_0[i] * pa_y[i] - ta_zzzz_xxxyy_1[i] * pc_y[i];

        ta_yzzzz_xxxyz_0[i] = ta_zzzz_xxxz_0[i] * fe_0 - ta_zzzz_xxxz_1[i] * fe_0 + ta_zzzz_xxxyz_0[i] * pa_y[i] - ta_zzzz_xxxyz_1[i] * pc_y[i];

        ta_yzzzz_xxxzz_0[i] = ta_zzzz_xxxzz_0[i] * pa_y[i] - ta_zzzz_xxxzz_1[i] * pc_y[i];

        ta_yzzzz_xxyyy_0[i] =
            3.0 * ta_zzzz_xxyy_0[i] * fe_0 - 3.0 * ta_zzzz_xxyy_1[i] * fe_0 + ta_zzzz_xxyyy_0[i] * pa_y[i] - ta_zzzz_xxyyy_1[i] * pc_y[i];

        ta_yzzzz_xxyyz_0[i] =
            2.0 * ta_zzzz_xxyz_0[i] * fe_0 - 2.0 * ta_zzzz_xxyz_1[i] * fe_0 + ta_zzzz_xxyyz_0[i] * pa_y[i] - ta_zzzz_xxyyz_1[i] * pc_y[i];

        ta_yzzzz_xxyzz_0[i] = ta_zzzz_xxzz_0[i] * fe_0 - ta_zzzz_xxzz_1[i] * fe_0 + ta_zzzz_xxyzz_0[i] * pa_y[i] - ta_zzzz_xxyzz_1[i] * pc_y[i];

        ta_yzzzz_xxzzz_0[i] = ta_zzzz_xxzzz_0[i] * pa_y[i] - ta_zzzz_xxzzz_1[i] * pc_y[i];

        ta_yzzzz_xyyyy_0[i] =
            4.0 * ta_zzzz_xyyy_0[i] * fe_0 - 4.0 * ta_zzzz_xyyy_1[i] * fe_0 + ta_zzzz_xyyyy_0[i] * pa_y[i] - ta_zzzz_xyyyy_1[i] * pc_y[i];

        ta_yzzzz_xyyyz_0[i] =
            3.0 * ta_zzzz_xyyz_0[i] * fe_0 - 3.0 * ta_zzzz_xyyz_1[i] * fe_0 + ta_zzzz_xyyyz_0[i] * pa_y[i] - ta_zzzz_xyyyz_1[i] * pc_y[i];

        ta_yzzzz_xyyzz_0[i] =
            2.0 * ta_zzzz_xyzz_0[i] * fe_0 - 2.0 * ta_zzzz_xyzz_1[i] * fe_0 + ta_zzzz_xyyzz_0[i] * pa_y[i] - ta_zzzz_xyyzz_1[i] * pc_y[i];

        ta_yzzzz_xyzzz_0[i] = ta_zzzz_xzzz_0[i] * fe_0 - ta_zzzz_xzzz_1[i] * fe_0 + ta_zzzz_xyzzz_0[i] * pa_y[i] - ta_zzzz_xyzzz_1[i] * pc_y[i];

        ta_yzzzz_xzzzz_0[i] = ta_zzzz_xzzzz_0[i] * pa_y[i] - ta_zzzz_xzzzz_1[i] * pc_y[i];

        ta_yzzzz_yyyyy_0[i] =
            5.0 * ta_zzzz_yyyy_0[i] * fe_0 - 5.0 * ta_zzzz_yyyy_1[i] * fe_0 + ta_zzzz_yyyyy_0[i] * pa_y[i] - ta_zzzz_yyyyy_1[i] * pc_y[i];

        ta_yzzzz_yyyyz_0[i] =
            4.0 * ta_zzzz_yyyz_0[i] * fe_0 - 4.0 * ta_zzzz_yyyz_1[i] * fe_0 + ta_zzzz_yyyyz_0[i] * pa_y[i] - ta_zzzz_yyyyz_1[i] * pc_y[i];

        ta_yzzzz_yyyzz_0[i] =
            3.0 * ta_zzzz_yyzz_0[i] * fe_0 - 3.0 * ta_zzzz_yyzz_1[i] * fe_0 + ta_zzzz_yyyzz_0[i] * pa_y[i] - ta_zzzz_yyyzz_1[i] * pc_y[i];

        ta_yzzzz_yyzzz_0[i] =
            2.0 * ta_zzzz_yzzz_0[i] * fe_0 - 2.0 * ta_zzzz_yzzz_1[i] * fe_0 + ta_zzzz_yyzzz_0[i] * pa_y[i] - ta_zzzz_yyzzz_1[i] * pc_y[i];

        ta_yzzzz_yzzzz_0[i] = ta_zzzz_zzzz_0[i] * fe_0 - ta_zzzz_zzzz_1[i] * fe_0 + ta_zzzz_yzzzz_0[i] * pa_y[i] - ta_zzzz_yzzzz_1[i] * pc_y[i];

        ta_yzzzz_zzzzz_0[i] = ta_zzzz_zzzzz_0[i] * pa_y[i] - ta_zzzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 420-441 components of targeted buffer : HH

    auto ta_zzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 420);

    auto ta_zzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 421);

    auto ta_zzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 422);

    auto ta_zzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 423);

    auto ta_zzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 424);

    auto ta_zzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 425);

    auto ta_zzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 426);

    auto ta_zzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 427);

    auto ta_zzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 428);

    auto ta_zzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 429);

    auto ta_zzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 430);

    auto ta_zzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 431);

    auto ta_zzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 432);

    auto ta_zzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 433);

    auto ta_zzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 434);

    auto ta_zzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 435);

    auto ta_zzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 436);

    auto ta_zzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 437);

    auto ta_zzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 438);

    auto ta_zzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 439);

    auto ta_zzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 440);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta_zzz_xxxxx_0,   \
                             ta_zzz_xxxxx_1,   \
                             ta_zzz_xxxxy_0,   \
                             ta_zzz_xxxxy_1,   \
                             ta_zzz_xxxxz_0,   \
                             ta_zzz_xxxxz_1,   \
                             ta_zzz_xxxyy_0,   \
                             ta_zzz_xxxyy_1,   \
                             ta_zzz_xxxyz_0,   \
                             ta_zzz_xxxyz_1,   \
                             ta_zzz_xxxzz_0,   \
                             ta_zzz_xxxzz_1,   \
                             ta_zzz_xxyyy_0,   \
                             ta_zzz_xxyyy_1,   \
                             ta_zzz_xxyyz_0,   \
                             ta_zzz_xxyyz_1,   \
                             ta_zzz_xxyzz_0,   \
                             ta_zzz_xxyzz_1,   \
                             ta_zzz_xxzzz_0,   \
                             ta_zzz_xxzzz_1,   \
                             ta_zzz_xyyyy_0,   \
                             ta_zzz_xyyyy_1,   \
                             ta_zzz_xyyyz_0,   \
                             ta_zzz_xyyyz_1,   \
                             ta_zzz_xyyzz_0,   \
                             ta_zzz_xyyzz_1,   \
                             ta_zzz_xyzzz_0,   \
                             ta_zzz_xyzzz_1,   \
                             ta_zzz_xzzzz_0,   \
                             ta_zzz_xzzzz_1,   \
                             ta_zzz_yyyyy_0,   \
                             ta_zzz_yyyyy_1,   \
                             ta_zzz_yyyyz_0,   \
                             ta_zzz_yyyyz_1,   \
                             ta_zzz_yyyzz_0,   \
                             ta_zzz_yyyzz_1,   \
                             ta_zzz_yyzzz_0,   \
                             ta_zzz_yyzzz_1,   \
                             ta_zzz_yzzzz_0,   \
                             ta_zzz_yzzzz_1,   \
                             ta_zzz_zzzzz_0,   \
                             ta_zzz_zzzzz_1,   \
                             ta_zzzz_xxxx_0,   \
                             ta_zzzz_xxxx_1,   \
                             ta_zzzz_xxxxx_0,  \
                             ta_zzzz_xxxxx_1,  \
                             ta_zzzz_xxxxy_0,  \
                             ta_zzzz_xxxxy_1,  \
                             ta_zzzz_xxxxz_0,  \
                             ta_zzzz_xxxxz_1,  \
                             ta_zzzz_xxxy_0,   \
                             ta_zzzz_xxxy_1,   \
                             ta_zzzz_xxxyy_0,  \
                             ta_zzzz_xxxyy_1,  \
                             ta_zzzz_xxxyz_0,  \
                             ta_zzzz_xxxyz_1,  \
                             ta_zzzz_xxxz_0,   \
                             ta_zzzz_xxxz_1,   \
                             ta_zzzz_xxxzz_0,  \
                             ta_zzzz_xxxzz_1,  \
                             ta_zzzz_xxyy_0,   \
                             ta_zzzz_xxyy_1,   \
                             ta_zzzz_xxyyy_0,  \
                             ta_zzzz_xxyyy_1,  \
                             ta_zzzz_xxyyz_0,  \
                             ta_zzzz_xxyyz_1,  \
                             ta_zzzz_xxyz_0,   \
                             ta_zzzz_xxyz_1,   \
                             ta_zzzz_xxyzz_0,  \
                             ta_zzzz_xxyzz_1,  \
                             ta_zzzz_xxzz_0,   \
                             ta_zzzz_xxzz_1,   \
                             ta_zzzz_xxzzz_0,  \
                             ta_zzzz_xxzzz_1,  \
                             ta_zzzz_xyyy_0,   \
                             ta_zzzz_xyyy_1,   \
                             ta_zzzz_xyyyy_0,  \
                             ta_zzzz_xyyyy_1,  \
                             ta_zzzz_xyyyz_0,  \
                             ta_zzzz_xyyyz_1,  \
                             ta_zzzz_xyyz_0,   \
                             ta_zzzz_xyyz_1,   \
                             ta_zzzz_xyyzz_0,  \
                             ta_zzzz_xyyzz_1,  \
                             ta_zzzz_xyzz_0,   \
                             ta_zzzz_xyzz_1,   \
                             ta_zzzz_xyzzz_0,  \
                             ta_zzzz_xyzzz_1,  \
                             ta_zzzz_xzzz_0,   \
                             ta_zzzz_xzzz_1,   \
                             ta_zzzz_xzzzz_0,  \
                             ta_zzzz_xzzzz_1,  \
                             ta_zzzz_yyyy_0,   \
                             ta_zzzz_yyyy_1,   \
                             ta_zzzz_yyyyy_0,  \
                             ta_zzzz_yyyyy_1,  \
                             ta_zzzz_yyyyz_0,  \
                             ta_zzzz_yyyyz_1,  \
                             ta_zzzz_yyyz_0,   \
                             ta_zzzz_yyyz_1,   \
                             ta_zzzz_yyyzz_0,  \
                             ta_zzzz_yyyzz_1,  \
                             ta_zzzz_yyzz_0,   \
                             ta_zzzz_yyzz_1,   \
                             ta_zzzz_yyzzz_0,  \
                             ta_zzzz_yyzzz_1,  \
                             ta_zzzz_yzzz_0,   \
                             ta_zzzz_yzzz_1,   \
                             ta_zzzz_yzzzz_0,  \
                             ta_zzzz_yzzzz_1,  \
                             ta_zzzz_zzzz_0,   \
                             ta_zzzz_zzzz_1,   \
                             ta_zzzz_zzzzz_0,  \
                             ta_zzzz_zzzzz_1,  \
                             ta_zzzzz_xxxxx_0, \
                             ta_zzzzz_xxxxy_0, \
                             ta_zzzzz_xxxxz_0, \
                             ta_zzzzz_xxxyy_0, \
                             ta_zzzzz_xxxyz_0, \
                             ta_zzzzz_xxxzz_0, \
                             ta_zzzzz_xxyyy_0, \
                             ta_zzzzz_xxyyz_0, \
                             ta_zzzzz_xxyzz_0, \
                             ta_zzzzz_xxzzz_0, \
                             ta_zzzzz_xyyyy_0, \
                             ta_zzzzz_xyyyz_0, \
                             ta_zzzzz_xyyzz_0, \
                             ta_zzzzz_xyzzz_0, \
                             ta_zzzzz_xzzzz_0, \
                             ta_zzzzz_yyyyy_0, \
                             ta_zzzzz_yyyyz_0, \
                             ta_zzzzz_yyyzz_0, \
                             ta_zzzzz_yyzzz_0, \
                             ta_zzzzz_yzzzz_0, \
                             ta_zzzzz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzz_xxxxx_0[i] =
            4.0 * ta_zzz_xxxxx_0[i] * fe_0 - 4.0 * ta_zzz_xxxxx_1[i] * fe_0 + ta_zzzz_xxxxx_0[i] * pa_z[i] - ta_zzzz_xxxxx_1[i] * pc_z[i];

        ta_zzzzz_xxxxy_0[i] =
            4.0 * ta_zzz_xxxxy_0[i] * fe_0 - 4.0 * ta_zzz_xxxxy_1[i] * fe_0 + ta_zzzz_xxxxy_0[i] * pa_z[i] - ta_zzzz_xxxxy_1[i] * pc_z[i];

        ta_zzzzz_xxxxz_0[i] = 4.0 * ta_zzz_xxxxz_0[i] * fe_0 - 4.0 * ta_zzz_xxxxz_1[i] * fe_0 + ta_zzzz_xxxx_0[i] * fe_0 - ta_zzzz_xxxx_1[i] * fe_0 +
                              ta_zzzz_xxxxz_0[i] * pa_z[i] - ta_zzzz_xxxxz_1[i] * pc_z[i];

        ta_zzzzz_xxxyy_0[i] =
            4.0 * ta_zzz_xxxyy_0[i] * fe_0 - 4.0 * ta_zzz_xxxyy_1[i] * fe_0 + ta_zzzz_xxxyy_0[i] * pa_z[i] - ta_zzzz_xxxyy_1[i] * pc_z[i];

        ta_zzzzz_xxxyz_0[i] = 4.0 * ta_zzz_xxxyz_0[i] * fe_0 - 4.0 * ta_zzz_xxxyz_1[i] * fe_0 + ta_zzzz_xxxy_0[i] * fe_0 - ta_zzzz_xxxy_1[i] * fe_0 +
                              ta_zzzz_xxxyz_0[i] * pa_z[i] - ta_zzzz_xxxyz_1[i] * pc_z[i];

        ta_zzzzz_xxxzz_0[i] = 4.0 * ta_zzz_xxxzz_0[i] * fe_0 - 4.0 * ta_zzz_xxxzz_1[i] * fe_0 + 2.0 * ta_zzzz_xxxz_0[i] * fe_0 -
                              2.0 * ta_zzzz_xxxz_1[i] * fe_0 + ta_zzzz_xxxzz_0[i] * pa_z[i] - ta_zzzz_xxxzz_1[i] * pc_z[i];

        ta_zzzzz_xxyyy_0[i] =
            4.0 * ta_zzz_xxyyy_0[i] * fe_0 - 4.0 * ta_zzz_xxyyy_1[i] * fe_0 + ta_zzzz_xxyyy_0[i] * pa_z[i] - ta_zzzz_xxyyy_1[i] * pc_z[i];

        ta_zzzzz_xxyyz_0[i] = 4.0 * ta_zzz_xxyyz_0[i] * fe_0 - 4.0 * ta_zzz_xxyyz_1[i] * fe_0 + ta_zzzz_xxyy_0[i] * fe_0 - ta_zzzz_xxyy_1[i] * fe_0 +
                              ta_zzzz_xxyyz_0[i] * pa_z[i] - ta_zzzz_xxyyz_1[i] * pc_z[i];

        ta_zzzzz_xxyzz_0[i] = 4.0 * ta_zzz_xxyzz_0[i] * fe_0 - 4.0 * ta_zzz_xxyzz_1[i] * fe_0 + 2.0 * ta_zzzz_xxyz_0[i] * fe_0 -
                              2.0 * ta_zzzz_xxyz_1[i] * fe_0 + ta_zzzz_xxyzz_0[i] * pa_z[i] - ta_zzzz_xxyzz_1[i] * pc_z[i];

        ta_zzzzz_xxzzz_0[i] = 4.0 * ta_zzz_xxzzz_0[i] * fe_0 - 4.0 * ta_zzz_xxzzz_1[i] * fe_0 + 3.0 * ta_zzzz_xxzz_0[i] * fe_0 -
                              3.0 * ta_zzzz_xxzz_1[i] * fe_0 + ta_zzzz_xxzzz_0[i] * pa_z[i] - ta_zzzz_xxzzz_1[i] * pc_z[i];

        ta_zzzzz_xyyyy_0[i] =
            4.0 * ta_zzz_xyyyy_0[i] * fe_0 - 4.0 * ta_zzz_xyyyy_1[i] * fe_0 + ta_zzzz_xyyyy_0[i] * pa_z[i] - ta_zzzz_xyyyy_1[i] * pc_z[i];

        ta_zzzzz_xyyyz_0[i] = 4.0 * ta_zzz_xyyyz_0[i] * fe_0 - 4.0 * ta_zzz_xyyyz_1[i] * fe_0 + ta_zzzz_xyyy_0[i] * fe_0 - ta_zzzz_xyyy_1[i] * fe_0 +
                              ta_zzzz_xyyyz_0[i] * pa_z[i] - ta_zzzz_xyyyz_1[i] * pc_z[i];

        ta_zzzzz_xyyzz_0[i] = 4.0 * ta_zzz_xyyzz_0[i] * fe_0 - 4.0 * ta_zzz_xyyzz_1[i] * fe_0 + 2.0 * ta_zzzz_xyyz_0[i] * fe_0 -
                              2.0 * ta_zzzz_xyyz_1[i] * fe_0 + ta_zzzz_xyyzz_0[i] * pa_z[i] - ta_zzzz_xyyzz_1[i] * pc_z[i];

        ta_zzzzz_xyzzz_0[i] = 4.0 * ta_zzz_xyzzz_0[i] * fe_0 - 4.0 * ta_zzz_xyzzz_1[i] * fe_0 + 3.0 * ta_zzzz_xyzz_0[i] * fe_0 -
                              3.0 * ta_zzzz_xyzz_1[i] * fe_0 + ta_zzzz_xyzzz_0[i] * pa_z[i] - ta_zzzz_xyzzz_1[i] * pc_z[i];

        ta_zzzzz_xzzzz_0[i] = 4.0 * ta_zzz_xzzzz_0[i] * fe_0 - 4.0 * ta_zzz_xzzzz_1[i] * fe_0 + 4.0 * ta_zzzz_xzzz_0[i] * fe_0 -
                              4.0 * ta_zzzz_xzzz_1[i] * fe_0 + ta_zzzz_xzzzz_0[i] * pa_z[i] - ta_zzzz_xzzzz_1[i] * pc_z[i];

        ta_zzzzz_yyyyy_0[i] =
            4.0 * ta_zzz_yyyyy_0[i] * fe_0 - 4.0 * ta_zzz_yyyyy_1[i] * fe_0 + ta_zzzz_yyyyy_0[i] * pa_z[i] - ta_zzzz_yyyyy_1[i] * pc_z[i];

        ta_zzzzz_yyyyz_0[i] = 4.0 * ta_zzz_yyyyz_0[i] * fe_0 - 4.0 * ta_zzz_yyyyz_1[i] * fe_0 + ta_zzzz_yyyy_0[i] * fe_0 - ta_zzzz_yyyy_1[i] * fe_0 +
                              ta_zzzz_yyyyz_0[i] * pa_z[i] - ta_zzzz_yyyyz_1[i] * pc_z[i];

        ta_zzzzz_yyyzz_0[i] = 4.0 * ta_zzz_yyyzz_0[i] * fe_0 - 4.0 * ta_zzz_yyyzz_1[i] * fe_0 + 2.0 * ta_zzzz_yyyz_0[i] * fe_0 -
                              2.0 * ta_zzzz_yyyz_1[i] * fe_0 + ta_zzzz_yyyzz_0[i] * pa_z[i] - ta_zzzz_yyyzz_1[i] * pc_z[i];

        ta_zzzzz_yyzzz_0[i] = 4.0 * ta_zzz_yyzzz_0[i] * fe_0 - 4.0 * ta_zzz_yyzzz_1[i] * fe_0 + 3.0 * ta_zzzz_yyzz_0[i] * fe_0 -
                              3.0 * ta_zzzz_yyzz_1[i] * fe_0 + ta_zzzz_yyzzz_0[i] * pa_z[i] - ta_zzzz_yyzzz_1[i] * pc_z[i];

        ta_zzzzz_yzzzz_0[i] = 4.0 * ta_zzz_yzzzz_0[i] * fe_0 - 4.0 * ta_zzz_yzzzz_1[i] * fe_0 + 4.0 * ta_zzzz_yzzz_0[i] * fe_0 -
                              4.0 * ta_zzzz_yzzz_1[i] * fe_0 + ta_zzzz_yzzzz_0[i] * pa_z[i] - ta_zzzz_yzzzz_1[i] * pc_z[i];

        ta_zzzzz_zzzzz_0[i] = 4.0 * ta_zzz_zzzzz_0[i] * fe_0 - 4.0 * ta_zzz_zzzzz_1[i] * fe_0 + 5.0 * ta_zzzz_zzzz_0[i] * fe_0 -
                              5.0 * ta_zzzz_zzzz_1[i] * fe_0 + ta_zzzz_zzzzz_0[i] * pa_z[i] - ta_zzzz_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
