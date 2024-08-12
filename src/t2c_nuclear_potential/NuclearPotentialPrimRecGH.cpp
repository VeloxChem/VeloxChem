#include "NuclearPotentialPrimRecGH.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_gh(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_gh,
                               const size_t idx_npot_0_dh,
                               const size_t idx_npot_1_dh,
                               const size_t idx_npot_0_fg,
                               const size_t idx_npot_1_fg,
                               const size_t idx_npot_0_fh,
                               const size_t idx_npot_1_fh,
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

    // Set up components of auxiliary buffer : DH

    auto ta_xx_xxxxx_0 = pbuffer.data(idx_npot_0_dh);

    auto ta_xx_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 1);

    auto ta_xx_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 2);

    auto ta_xx_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 3);

    auto ta_xx_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 4);

    auto ta_xx_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 5);

    auto ta_xx_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 6);

    auto ta_xx_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 7);

    auto ta_xx_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 8);

    auto ta_xx_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 9);

    auto ta_xx_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 10);

    auto ta_xx_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 11);

    auto ta_xx_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 12);

    auto ta_xx_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 13);

    auto ta_xx_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 14);

    auto ta_xx_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 15);

    auto ta_xx_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 16);

    auto ta_xx_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 17);

    auto ta_xx_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 18);

    auto ta_xx_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 19);

    auto ta_xx_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 20);

    auto ta_xy_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 36);

    auto ta_xy_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 37);

    auto ta_xy_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 38);

    auto ta_xy_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 39);

    auto ta_xy_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 40);

    auto ta_xz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 58);

    auto ta_xz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 59);

    auto ta_xz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 60);

    auto ta_xz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 61);

    auto ta_xz_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 62);

    auto ta_yy_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 63);

    auto ta_yy_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 64);

    auto ta_yy_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 65);

    auto ta_yy_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 66);

    auto ta_yy_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 67);

    auto ta_yy_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 68);

    auto ta_yy_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 69);

    auto ta_yy_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 70);

    auto ta_yy_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 71);

    auto ta_yy_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 72);

    auto ta_yy_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 73);

    auto ta_yy_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 74);

    auto ta_yy_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 75);

    auto ta_yy_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 76);

    auto ta_yy_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 77);

    auto ta_yy_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 78);

    auto ta_yy_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 79);

    auto ta_yy_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 80);

    auto ta_yy_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 81);

    auto ta_yy_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 82);

    auto ta_yy_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 83);

    auto ta_yz_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 86);

    auto ta_yz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 89);

    auto ta_yz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 93);

    auto ta_yz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 98);

    auto ta_yz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 100);

    auto ta_yz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 101);

    auto ta_yz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 102);

    auto ta_yz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 103);

    auto ta_yz_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 104);

    auto ta_zz_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 105);

    auto ta_zz_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 106);

    auto ta_zz_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 107);

    auto ta_zz_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 108);

    auto ta_zz_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 109);

    auto ta_zz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 110);

    auto ta_zz_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 111);

    auto ta_zz_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 112);

    auto ta_zz_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 113);

    auto ta_zz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 114);

    auto ta_zz_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 115);

    auto ta_zz_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 116);

    auto ta_zz_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 117);

    auto ta_zz_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 118);

    auto ta_zz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 119);

    auto ta_zz_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 120);

    auto ta_zz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 121);

    auto ta_zz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 122);

    auto ta_zz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 123);

    auto ta_zz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 124);

    auto ta_zz_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 125);

    // Set up components of auxiliary buffer : DH

    auto ta_xx_xxxxx_1 = pbuffer.data(idx_npot_1_dh);

    auto ta_xx_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 1);

    auto ta_xx_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 2);

    auto ta_xx_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 3);

    auto ta_xx_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 4);

    auto ta_xx_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 5);

    auto ta_xx_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 6);

    auto ta_xx_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 7);

    auto ta_xx_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 8);

    auto ta_xx_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 9);

    auto ta_xx_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 10);

    auto ta_xx_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 11);

    auto ta_xx_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 12);

    auto ta_xx_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 13);

    auto ta_xx_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 14);

    auto ta_xx_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 15);

    auto ta_xx_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 16);

    auto ta_xx_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 17);

    auto ta_xx_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 18);

    auto ta_xx_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 19);

    auto ta_xx_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 20);

    auto ta_xy_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 36);

    auto ta_xy_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 37);

    auto ta_xy_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 38);

    auto ta_xy_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 39);

    auto ta_xy_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 40);

    auto ta_xz_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 58);

    auto ta_xz_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 59);

    auto ta_xz_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 60);

    auto ta_xz_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 61);

    auto ta_xz_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 62);

    auto ta_yy_xxxxx_1 = pbuffer.data(idx_npot_1_dh + 63);

    auto ta_yy_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 64);

    auto ta_yy_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 65);

    auto ta_yy_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 66);

    auto ta_yy_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 67);

    auto ta_yy_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 68);

    auto ta_yy_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 69);

    auto ta_yy_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 70);

    auto ta_yy_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 71);

    auto ta_yy_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 72);

    auto ta_yy_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 73);

    auto ta_yy_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 74);

    auto ta_yy_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 75);

    auto ta_yy_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 76);

    auto ta_yy_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 77);

    auto ta_yy_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 78);

    auto ta_yy_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 79);

    auto ta_yy_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 80);

    auto ta_yy_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 81);

    auto ta_yy_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 82);

    auto ta_yy_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 83);

    auto ta_yz_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 86);

    auto ta_yz_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 89);

    auto ta_yz_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 93);

    auto ta_yz_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 98);

    auto ta_yz_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 100);

    auto ta_yz_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 101);

    auto ta_yz_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 102);

    auto ta_yz_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 103);

    auto ta_yz_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 104);

    auto ta_zz_xxxxx_1 = pbuffer.data(idx_npot_1_dh + 105);

    auto ta_zz_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 106);

    auto ta_zz_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 107);

    auto ta_zz_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 108);

    auto ta_zz_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 109);

    auto ta_zz_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 110);

    auto ta_zz_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 111);

    auto ta_zz_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 112);

    auto ta_zz_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 113);

    auto ta_zz_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 114);

    auto ta_zz_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 115);

    auto ta_zz_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 116);

    auto ta_zz_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 117);

    auto ta_zz_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 118);

    auto ta_zz_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 119);

    auto ta_zz_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 120);

    auto ta_zz_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 121);

    auto ta_zz_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 122);

    auto ta_zz_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 123);

    auto ta_zz_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 124);

    auto ta_zz_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 125);

    // Set up components of auxiliary buffer : FG

    auto ta_xxx_xxxx_0 = pbuffer.data(idx_npot_0_fg);

    auto ta_xxx_xxxy_0 = pbuffer.data(idx_npot_0_fg + 1);

    auto ta_xxx_xxxz_0 = pbuffer.data(idx_npot_0_fg + 2);

    auto ta_xxx_xxyy_0 = pbuffer.data(idx_npot_0_fg + 3);

    auto ta_xxx_xxyz_0 = pbuffer.data(idx_npot_0_fg + 4);

    auto ta_xxx_xxzz_0 = pbuffer.data(idx_npot_0_fg + 5);

    auto ta_xxx_xyyy_0 = pbuffer.data(idx_npot_0_fg + 6);

    auto ta_xxx_xyyz_0 = pbuffer.data(idx_npot_0_fg + 7);

    auto ta_xxx_xyzz_0 = pbuffer.data(idx_npot_0_fg + 8);

    auto ta_xxx_xzzz_0 = pbuffer.data(idx_npot_0_fg + 9);

    auto ta_xxx_yyyy_0 = pbuffer.data(idx_npot_0_fg + 10);

    auto ta_xxx_yyyz_0 = pbuffer.data(idx_npot_0_fg + 11);

    auto ta_xxx_yyzz_0 = pbuffer.data(idx_npot_0_fg + 12);

    auto ta_xxx_yzzz_0 = pbuffer.data(idx_npot_0_fg + 13);

    auto ta_xxx_zzzz_0 = pbuffer.data(idx_npot_0_fg + 14);

    auto ta_xxz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 32);

    auto ta_xxz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 34);

    auto ta_xxz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 35);

    auto ta_xxz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 37);

    auto ta_xxz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 38);

    auto ta_xxz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 39);

    auto ta_xyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 46);

    auto ta_xyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 48);

    auto ta_xyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 49);

    auto ta_xyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 51);

    auto ta_xyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 52);

    auto ta_xyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 53);

    auto ta_xyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 55);

    auto ta_xyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 56);

    auto ta_xyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 57);

    auto ta_xyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 58);

    auto ta_xzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 77);

    auto ta_xzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 79);

    auto ta_xzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 80);

    auto ta_xzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 82);

    auto ta_xzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 83);

    auto ta_xzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 84);

    auto ta_xzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 86);

    auto ta_xzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 87);

    auto ta_xzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 88);

    auto ta_xzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 89);

    auto ta_yyy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 90);

    auto ta_yyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 91);

    auto ta_yyy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 92);

    auto ta_yyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 93);

    auto ta_yyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 94);

    auto ta_yyy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 95);

    auto ta_yyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 96);

    auto ta_yyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 97);

    auto ta_yyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 98);

    auto ta_yyy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 99);

    auto ta_yyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 100);

    auto ta_yyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 101);

    auto ta_yyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 102);

    auto ta_yyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 103);

    auto ta_yyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 104);

    auto ta_yyz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 107);

    auto ta_yyz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 109);

    auto ta_yyz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 110);

    auto ta_yyz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 112);

    auto ta_yyz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 113);

    auto ta_yyz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 114);

    auto ta_yyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 116);

    auto ta_yyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 117);

    auto ta_yyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 118);

    auto ta_yyz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 119);

    auto ta_yzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 121);

    auto ta_yzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 122);

    auto ta_yzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 123);

    auto ta_yzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 124);

    auto ta_yzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 125);

    auto ta_yzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 126);

    auto ta_yzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 127);

    auto ta_yzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 128);

    auto ta_yzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 129);

    auto ta_yzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 130);

    auto ta_yzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 131);

    auto ta_yzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 132);

    auto ta_yzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 133);

    auto ta_yzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 134);

    auto ta_zzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 135);

    auto ta_zzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 136);

    auto ta_zzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 137);

    auto ta_zzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 138);

    auto ta_zzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 139);

    auto ta_zzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 140);

    auto ta_zzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 141);

    auto ta_zzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 142);

    auto ta_zzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 143);

    auto ta_zzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 144);

    auto ta_zzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 145);

    auto ta_zzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 146);

    auto ta_zzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 147);

    auto ta_zzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 148);

    auto ta_zzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 149);

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

    auto ta_xxz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 32);

    auto ta_xxz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 34);

    auto ta_xxz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 35);

    auto ta_xxz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 37);

    auto ta_xxz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 38);

    auto ta_xxz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 39);

    auto ta_xyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 46);

    auto ta_xyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 48);

    auto ta_xyy_xxyz_1 = pbuffer.data(idx_npot_1_fg + 49);

    auto ta_xyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 51);

    auto ta_xyy_xyyz_1 = pbuffer.data(idx_npot_1_fg + 52);

    auto ta_xyy_xyzz_1 = pbuffer.data(idx_npot_1_fg + 53);

    auto ta_xyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 55);

    auto ta_xyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 56);

    auto ta_xyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 57);

    auto ta_xyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 58);

    auto ta_xzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 77);

    auto ta_xzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 79);

    auto ta_xzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 80);

    auto ta_xzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 82);

    auto ta_xzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 83);

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

    auto ta_yyz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 107);

    auto ta_yyz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 109);

    auto ta_yyz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 110);

    auto ta_yyz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 112);

    auto ta_yyz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 113);

    auto ta_yyz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 114);

    auto ta_yyz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 116);

    auto ta_yyz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 117);

    auto ta_yyz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 118);

    auto ta_yyz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 119);

    auto ta_yzz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 121);

    auto ta_yzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 122);

    auto ta_yzz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 123);

    auto ta_yzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 124);

    auto ta_yzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 125);

    auto ta_yzz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 126);

    auto ta_yzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 127);

    auto ta_yzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 128);

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

    auto ta_xxy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 22);

    auto ta_xxy_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 23);

    auto ta_xxy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 24);

    auto ta_xxy_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 26);

    auto ta_xxy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 27);

    auto ta_xxy_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 30);

    auto ta_xxy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 31);

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

    auto ta_xxz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 46);

    auto ta_xxz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 47);

    auto ta_xxz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 48);

    auto ta_xxz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 49);

    auto ta_xxz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 50);

    auto ta_xxz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 51);

    auto ta_xxz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 52);

    auto ta_xxz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 53);

    auto ta_xxz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 54);

    auto ta_xxz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 55);

    auto ta_xxz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 56);

    auto ta_xxz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 58);

    auto ta_xxz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 59);

    auto ta_xxz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 60);

    auto ta_xxz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 61);

    auto ta_xxz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 62);

    auto ta_xyy_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 63);

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

    auto ta_xzz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 105);

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

    auto ta_yyz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 151);

    auto ta_yyz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 152);

    auto ta_yyz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 153);

    auto ta_yyz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 154);

    auto ta_yyz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 155);

    auto ta_yyz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 156);

    auto ta_yyz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 157);

    auto ta_yyz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 158);

    auto ta_yyz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 159);

    auto ta_yyz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 160);

    auto ta_yyz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 161);

    auto ta_yyz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 162);

    auto ta_yyz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 163);

    auto ta_yyz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 164);

    auto ta_yyz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 165);

    auto ta_yyz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 166);

    auto ta_yyz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 167);

    auto ta_yzz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 168);

    auto ta_yzz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 169);

    auto ta_yzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 170);

    auto ta_yzz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 171);

    auto ta_yzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 172);

    auto ta_yzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 173);

    auto ta_yzz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 174);

    auto ta_yzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 175);

    auto ta_yzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 176);

    auto ta_yzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 177);

    auto ta_yzz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 178);

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

    auto ta_xxy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 22);

    auto ta_xxy_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 23);

    auto ta_xxy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 24);

    auto ta_xxy_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 26);

    auto ta_xxy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 27);

    auto ta_xxy_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 30);

    auto ta_xxy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 31);

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

    auto ta_xxz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 46);

    auto ta_xxz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 47);

    auto ta_xxz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 48);

    auto ta_xxz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 49);

    auto ta_xxz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 50);

    auto ta_xxz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 51);

    auto ta_xxz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 52);

    auto ta_xxz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 53);

    auto ta_xxz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 54);

    auto ta_xxz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 55);

    auto ta_xxz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 56);

    auto ta_xxz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 58);

    auto ta_xxz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 59);

    auto ta_xxz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 60);

    auto ta_xxz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 61);

    auto ta_xxz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 62);

    auto ta_xyy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 63);

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

    auto ta_xzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 105);

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

    auto ta_yyz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 151);

    auto ta_yyz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 152);

    auto ta_yyz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 153);

    auto ta_yyz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 154);

    auto ta_yyz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 155);

    auto ta_yyz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 156);

    auto ta_yyz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 157);

    auto ta_yyz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 158);

    auto ta_yyz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 159);

    auto ta_yyz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 160);

    auto ta_yyz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 161);

    auto ta_yyz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 162);

    auto ta_yyz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 163);

    auto ta_yyz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 164);

    auto ta_yyz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 165);

    auto ta_yyz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 166);

    auto ta_yyz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 167);

    auto ta_yzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 168);

    auto ta_yzz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 169);

    auto ta_yzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 170);

    auto ta_yzz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 171);

    auto ta_yzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 172);

    auto ta_yzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 173);

    auto ta_yzz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 174);

    auto ta_yzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 175);

    auto ta_yzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 176);

    auto ta_yzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 177);

    auto ta_yzz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 178);

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

    // Set up 0-21 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pc_x, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_yyyyy_0, ta_xx_yyyyy_1, ta_xx_yyyyz_0, ta_xx_yyyyz_1, ta_xx_yyyzz_0, ta_xx_yyyzz_1, ta_xx_yyzzz_0, ta_xx_yyzzz_1, ta_xx_yzzzz_0, ta_xx_yzzzz_1, ta_xx_zzzzz_0, ta_xx_zzzzz_1, ta_xxx_xxxx_0, ta_xxx_xxxx_1, ta_xxx_xxxxx_0, ta_xxx_xxxxx_1, ta_xxx_xxxxy_0, ta_xxx_xxxxy_1, ta_xxx_xxxxz_0, ta_xxx_xxxxz_1, ta_xxx_xxxy_0, ta_xxx_xxxy_1, ta_xxx_xxxyy_0, ta_xxx_xxxyy_1, ta_xxx_xxxyz_0, ta_xxx_xxxyz_1, ta_xxx_xxxz_0, ta_xxx_xxxz_1, ta_xxx_xxxzz_0, ta_xxx_xxxzz_1, ta_xxx_xxyy_0, ta_xxx_xxyy_1, ta_xxx_xxyyy_0, ta_xxx_xxyyy_1, ta_xxx_xxyyz_0, ta_xxx_xxyyz_1, ta_xxx_xxyz_0, ta_xxx_xxyz_1, ta_xxx_xxyzz_0, ta_xxx_xxyzz_1, ta_xxx_xxzz_0, ta_xxx_xxzz_1, ta_xxx_xxzzz_0, ta_xxx_xxzzz_1, ta_xxx_xyyy_0, ta_xxx_xyyy_1, ta_xxx_xyyyy_0, ta_xxx_xyyyy_1, ta_xxx_xyyyz_0, ta_xxx_xyyyz_1, ta_xxx_xyyz_0, ta_xxx_xyyz_1, ta_xxx_xyyzz_0, ta_xxx_xyyzz_1, ta_xxx_xyzz_0, ta_xxx_xyzz_1, ta_xxx_xyzzz_0, ta_xxx_xyzzz_1, ta_xxx_xzzz_0, ta_xxx_xzzz_1, ta_xxx_xzzzz_0, ta_xxx_xzzzz_1, ta_xxx_yyyy_0, ta_xxx_yyyy_1, ta_xxx_yyyyy_0, ta_xxx_yyyyy_1, ta_xxx_yyyyz_0, ta_xxx_yyyyz_1, ta_xxx_yyyz_0, ta_xxx_yyyz_1, ta_xxx_yyyzz_0, ta_xxx_yyyzz_1, ta_xxx_yyzz_0, ta_xxx_yyzz_1, ta_xxx_yyzzz_0, ta_xxx_yyzzz_1, ta_xxx_yzzz_0, ta_xxx_yzzz_1, ta_xxx_yzzzz_0, ta_xxx_yzzzz_1, ta_xxx_zzzz_0, ta_xxx_zzzz_1, ta_xxx_zzzzz_0, ta_xxx_zzzzz_1, ta_xxxx_xxxxx_0, ta_xxxx_xxxxy_0, ta_xxxx_xxxxz_0, ta_xxxx_xxxyy_0, ta_xxxx_xxxyz_0, ta_xxxx_xxxzz_0, ta_xxxx_xxyyy_0, ta_xxxx_xxyyz_0, ta_xxxx_xxyzz_0, ta_xxxx_xxzzz_0, ta_xxxx_xyyyy_0, ta_xxxx_xyyyz_0, ta_xxxx_xyyzz_0, ta_xxxx_xyzzz_0, ta_xxxx_xzzzz_0, ta_xxxx_yyyyy_0, ta_xxxx_yyyyz_0, ta_xxxx_yyyzz_0, ta_xxxx_yyzzz_0, ta_xxxx_yzzzz_0, ta_xxxx_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_xxxxx_0[i] = 3.0 * ta_xx_xxxxx_0[i] * fe_0 - 3.0 * ta_xx_xxxxx_1[i] * fe_0 + 5.0 * ta_xxx_xxxx_0[i] * fe_0 - 5.0 * ta_xxx_xxxx_1[i] * fe_0 + ta_xxx_xxxxx_0[i] * pa_x[i] - ta_xxx_xxxxx_1[i] * pc_x[i];

        ta_xxxx_xxxxy_0[i] = 3.0 * ta_xx_xxxxy_0[i] * fe_0 - 3.0 * ta_xx_xxxxy_1[i] * fe_0 + 4.0 * ta_xxx_xxxy_0[i] * fe_0 - 4.0 * ta_xxx_xxxy_1[i] * fe_0 + ta_xxx_xxxxy_0[i] * pa_x[i] - ta_xxx_xxxxy_1[i] * pc_x[i];

        ta_xxxx_xxxxz_0[i] = 3.0 * ta_xx_xxxxz_0[i] * fe_0 - 3.0 * ta_xx_xxxxz_1[i] * fe_0 + 4.0 * ta_xxx_xxxz_0[i] * fe_0 - 4.0 * ta_xxx_xxxz_1[i] * fe_0 + ta_xxx_xxxxz_0[i] * pa_x[i] - ta_xxx_xxxxz_1[i] * pc_x[i];

        ta_xxxx_xxxyy_0[i] = 3.0 * ta_xx_xxxyy_0[i] * fe_0 - 3.0 * ta_xx_xxxyy_1[i] * fe_0 + 3.0 * ta_xxx_xxyy_0[i] * fe_0 - 3.0 * ta_xxx_xxyy_1[i] * fe_0 + ta_xxx_xxxyy_0[i] * pa_x[i] - ta_xxx_xxxyy_1[i] * pc_x[i];

        ta_xxxx_xxxyz_0[i] = 3.0 * ta_xx_xxxyz_0[i] * fe_0 - 3.0 * ta_xx_xxxyz_1[i] * fe_0 + 3.0 * ta_xxx_xxyz_0[i] * fe_0 - 3.0 * ta_xxx_xxyz_1[i] * fe_0 + ta_xxx_xxxyz_0[i] * pa_x[i] - ta_xxx_xxxyz_1[i] * pc_x[i];

        ta_xxxx_xxxzz_0[i] = 3.0 * ta_xx_xxxzz_0[i] * fe_0 - 3.0 * ta_xx_xxxzz_1[i] * fe_0 + 3.0 * ta_xxx_xxzz_0[i] * fe_0 - 3.0 * ta_xxx_xxzz_1[i] * fe_0 + ta_xxx_xxxzz_0[i] * pa_x[i] - ta_xxx_xxxzz_1[i] * pc_x[i];

        ta_xxxx_xxyyy_0[i] = 3.0 * ta_xx_xxyyy_0[i] * fe_0 - 3.0 * ta_xx_xxyyy_1[i] * fe_0 + 2.0 * ta_xxx_xyyy_0[i] * fe_0 - 2.0 * ta_xxx_xyyy_1[i] * fe_0 + ta_xxx_xxyyy_0[i] * pa_x[i] - ta_xxx_xxyyy_1[i] * pc_x[i];

        ta_xxxx_xxyyz_0[i] = 3.0 * ta_xx_xxyyz_0[i] * fe_0 - 3.0 * ta_xx_xxyyz_1[i] * fe_0 + 2.0 * ta_xxx_xyyz_0[i] * fe_0 - 2.0 * ta_xxx_xyyz_1[i] * fe_0 + ta_xxx_xxyyz_0[i] * pa_x[i] - ta_xxx_xxyyz_1[i] * pc_x[i];

        ta_xxxx_xxyzz_0[i] = 3.0 * ta_xx_xxyzz_0[i] * fe_0 - 3.0 * ta_xx_xxyzz_1[i] * fe_0 + 2.0 * ta_xxx_xyzz_0[i] * fe_0 - 2.0 * ta_xxx_xyzz_1[i] * fe_0 + ta_xxx_xxyzz_0[i] * pa_x[i] - ta_xxx_xxyzz_1[i] * pc_x[i];

        ta_xxxx_xxzzz_0[i] = 3.0 * ta_xx_xxzzz_0[i] * fe_0 - 3.0 * ta_xx_xxzzz_1[i] * fe_0 + 2.0 * ta_xxx_xzzz_0[i] * fe_0 - 2.0 * ta_xxx_xzzz_1[i] * fe_0 + ta_xxx_xxzzz_0[i] * pa_x[i] - ta_xxx_xxzzz_1[i] * pc_x[i];

        ta_xxxx_xyyyy_0[i] = 3.0 * ta_xx_xyyyy_0[i] * fe_0 - 3.0 * ta_xx_xyyyy_1[i] * fe_0 + ta_xxx_yyyy_0[i] * fe_0 - ta_xxx_yyyy_1[i] * fe_0 + ta_xxx_xyyyy_0[i] * pa_x[i] - ta_xxx_xyyyy_1[i] * pc_x[i];

        ta_xxxx_xyyyz_0[i] = 3.0 * ta_xx_xyyyz_0[i] * fe_0 - 3.0 * ta_xx_xyyyz_1[i] * fe_0 + ta_xxx_yyyz_0[i] * fe_0 - ta_xxx_yyyz_1[i] * fe_0 + ta_xxx_xyyyz_0[i] * pa_x[i] - ta_xxx_xyyyz_1[i] * pc_x[i];

        ta_xxxx_xyyzz_0[i] = 3.0 * ta_xx_xyyzz_0[i] * fe_0 - 3.0 * ta_xx_xyyzz_1[i] * fe_0 + ta_xxx_yyzz_0[i] * fe_0 - ta_xxx_yyzz_1[i] * fe_0 + ta_xxx_xyyzz_0[i] * pa_x[i] - ta_xxx_xyyzz_1[i] * pc_x[i];

        ta_xxxx_xyzzz_0[i] = 3.0 * ta_xx_xyzzz_0[i] * fe_0 - 3.0 * ta_xx_xyzzz_1[i] * fe_0 + ta_xxx_yzzz_0[i] * fe_0 - ta_xxx_yzzz_1[i] * fe_0 + ta_xxx_xyzzz_0[i] * pa_x[i] - ta_xxx_xyzzz_1[i] * pc_x[i];

        ta_xxxx_xzzzz_0[i] = 3.0 * ta_xx_xzzzz_0[i] * fe_0 - 3.0 * ta_xx_xzzzz_1[i] * fe_0 + ta_xxx_zzzz_0[i] * fe_0 - ta_xxx_zzzz_1[i] * fe_0 + ta_xxx_xzzzz_0[i] * pa_x[i] - ta_xxx_xzzzz_1[i] * pc_x[i];

        ta_xxxx_yyyyy_0[i] = 3.0 * ta_xx_yyyyy_0[i] * fe_0 - 3.0 * ta_xx_yyyyy_1[i] * fe_0 + ta_xxx_yyyyy_0[i] * pa_x[i] - ta_xxx_yyyyy_1[i] * pc_x[i];

        ta_xxxx_yyyyz_0[i] = 3.0 * ta_xx_yyyyz_0[i] * fe_0 - 3.0 * ta_xx_yyyyz_1[i] * fe_0 + ta_xxx_yyyyz_0[i] * pa_x[i] - ta_xxx_yyyyz_1[i] * pc_x[i];

        ta_xxxx_yyyzz_0[i] = 3.0 * ta_xx_yyyzz_0[i] * fe_0 - 3.0 * ta_xx_yyyzz_1[i] * fe_0 + ta_xxx_yyyzz_0[i] * pa_x[i] - ta_xxx_yyyzz_1[i] * pc_x[i];

        ta_xxxx_yyzzz_0[i] = 3.0 * ta_xx_yyzzz_0[i] * fe_0 - 3.0 * ta_xx_yyzzz_1[i] * fe_0 + ta_xxx_yyzzz_0[i] * pa_x[i] - ta_xxx_yyzzz_1[i] * pc_x[i];

        ta_xxxx_yzzzz_0[i] = 3.0 * ta_xx_yzzzz_0[i] * fe_0 - 3.0 * ta_xx_yzzzz_1[i] * fe_0 + ta_xxx_yzzzz_0[i] * pa_x[i] - ta_xxx_yzzzz_1[i] * pc_x[i];

        ta_xxxx_zzzzz_0[i] = 3.0 * ta_xx_zzzzz_0[i] * fe_0 - 3.0 * ta_xx_zzzzz_1[i] * fe_0 + ta_xxx_zzzzz_0[i] * pa_x[i] - ta_xxx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : GH

    auto ta_xxxy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 21);

    auto ta_xxxy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 22);

    auto ta_xxxy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 23);

    auto ta_xxxy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 24);

    auto ta_xxxy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 25);

    auto ta_xxxy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 26);

    auto ta_xxxy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 27);

    auto ta_xxxy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 28);

    auto ta_xxxy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 29);

    auto ta_xxxy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 30);

    auto ta_xxxy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 31);

    auto ta_xxxy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 32);

    auto ta_xxxy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 33);

    auto ta_xxxy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 34);

    auto ta_xxxy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 35);

    auto ta_xxxy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 36);

    auto ta_xxxy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 37);

    auto ta_xxxy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 38);

    auto ta_xxxy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 39);

    auto ta_xxxy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 40);

    auto ta_xxxy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 41);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxx_xxxx_0, ta_xxx_xxxx_1, ta_xxx_xxxxx_0, ta_xxx_xxxxx_1, ta_xxx_xxxxy_0, ta_xxx_xxxxy_1, ta_xxx_xxxxz_0, ta_xxx_xxxxz_1, ta_xxx_xxxy_0, ta_xxx_xxxy_1, ta_xxx_xxxyy_0, ta_xxx_xxxyy_1, ta_xxx_xxxyz_0, ta_xxx_xxxyz_1, ta_xxx_xxxz_0, ta_xxx_xxxz_1, ta_xxx_xxxzz_0, ta_xxx_xxxzz_1, ta_xxx_xxyy_0, ta_xxx_xxyy_1, ta_xxx_xxyyy_0, ta_xxx_xxyyy_1, ta_xxx_xxyyz_0, ta_xxx_xxyyz_1, ta_xxx_xxyz_0, ta_xxx_xxyz_1, ta_xxx_xxyzz_0, ta_xxx_xxyzz_1, ta_xxx_xxzz_0, ta_xxx_xxzz_1, ta_xxx_xxzzz_0, ta_xxx_xxzzz_1, ta_xxx_xyyy_0, ta_xxx_xyyy_1, ta_xxx_xyyyy_0, ta_xxx_xyyyy_1, ta_xxx_xyyyz_0, ta_xxx_xyyyz_1, ta_xxx_xyyz_0, ta_xxx_xyyz_1, ta_xxx_xyyzz_0, ta_xxx_xyyzz_1, ta_xxx_xyzz_0, ta_xxx_xyzz_1, ta_xxx_xyzzz_0, ta_xxx_xyzzz_1, ta_xxx_xzzz_0, ta_xxx_xzzz_1, ta_xxx_xzzzz_0, ta_xxx_xzzzz_1, ta_xxx_zzzzz_0, ta_xxx_zzzzz_1, ta_xxxy_xxxxx_0, ta_xxxy_xxxxy_0, ta_xxxy_xxxxz_0, ta_xxxy_xxxyy_0, ta_xxxy_xxxyz_0, ta_xxxy_xxxzz_0, ta_xxxy_xxyyy_0, ta_xxxy_xxyyz_0, ta_xxxy_xxyzz_0, ta_xxxy_xxzzz_0, ta_xxxy_xyyyy_0, ta_xxxy_xyyyz_0, ta_xxxy_xyyzz_0, ta_xxxy_xyzzz_0, ta_xxxy_xzzzz_0, ta_xxxy_yyyyy_0, ta_xxxy_yyyyz_0, ta_xxxy_yyyzz_0, ta_xxxy_yyzzz_0, ta_xxxy_yzzzz_0, ta_xxxy_zzzzz_0, ta_xxy_yyyyy_0, ta_xxy_yyyyy_1, ta_xxy_yyyyz_0, ta_xxy_yyyyz_1, ta_xxy_yyyzz_0, ta_xxy_yyyzz_1, ta_xxy_yyzzz_0, ta_xxy_yyzzz_1, ta_xxy_yzzzz_0, ta_xxy_yzzzz_1, ta_xy_yyyyy_0, ta_xy_yyyyy_1, ta_xy_yyyyz_0, ta_xy_yyyyz_1, ta_xy_yyyzz_0, ta_xy_yyyzz_1, ta_xy_yyzzz_0, ta_xy_yyzzz_1, ta_xy_yzzzz_0, ta_xy_yzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxy_xxxxx_0[i] = ta_xxx_xxxxx_0[i] * pa_y[i] - ta_xxx_xxxxx_1[i] * pc_y[i];

        ta_xxxy_xxxxy_0[i] = ta_xxx_xxxx_0[i] * fe_0 - ta_xxx_xxxx_1[i] * fe_0 + ta_xxx_xxxxy_0[i] * pa_y[i] - ta_xxx_xxxxy_1[i] * pc_y[i];

        ta_xxxy_xxxxz_0[i] = ta_xxx_xxxxz_0[i] * pa_y[i] - ta_xxx_xxxxz_1[i] * pc_y[i];

        ta_xxxy_xxxyy_0[i] = 2.0 * ta_xxx_xxxy_0[i] * fe_0 - 2.0 * ta_xxx_xxxy_1[i] * fe_0 + ta_xxx_xxxyy_0[i] * pa_y[i] - ta_xxx_xxxyy_1[i] * pc_y[i];

        ta_xxxy_xxxyz_0[i] = ta_xxx_xxxz_0[i] * fe_0 - ta_xxx_xxxz_1[i] * fe_0 + ta_xxx_xxxyz_0[i] * pa_y[i] - ta_xxx_xxxyz_1[i] * pc_y[i];

        ta_xxxy_xxxzz_0[i] = ta_xxx_xxxzz_0[i] * pa_y[i] - ta_xxx_xxxzz_1[i] * pc_y[i];

        ta_xxxy_xxyyy_0[i] = 3.0 * ta_xxx_xxyy_0[i] * fe_0 - 3.0 * ta_xxx_xxyy_1[i] * fe_0 + ta_xxx_xxyyy_0[i] * pa_y[i] - ta_xxx_xxyyy_1[i] * pc_y[i];

        ta_xxxy_xxyyz_0[i] = 2.0 * ta_xxx_xxyz_0[i] * fe_0 - 2.0 * ta_xxx_xxyz_1[i] * fe_0 + ta_xxx_xxyyz_0[i] * pa_y[i] - ta_xxx_xxyyz_1[i] * pc_y[i];

        ta_xxxy_xxyzz_0[i] = ta_xxx_xxzz_0[i] * fe_0 - ta_xxx_xxzz_1[i] * fe_0 + ta_xxx_xxyzz_0[i] * pa_y[i] - ta_xxx_xxyzz_1[i] * pc_y[i];

        ta_xxxy_xxzzz_0[i] = ta_xxx_xxzzz_0[i] * pa_y[i] - ta_xxx_xxzzz_1[i] * pc_y[i];

        ta_xxxy_xyyyy_0[i] = 4.0 * ta_xxx_xyyy_0[i] * fe_0 - 4.0 * ta_xxx_xyyy_1[i] * fe_0 + ta_xxx_xyyyy_0[i] * pa_y[i] - ta_xxx_xyyyy_1[i] * pc_y[i];

        ta_xxxy_xyyyz_0[i] = 3.0 * ta_xxx_xyyz_0[i] * fe_0 - 3.0 * ta_xxx_xyyz_1[i] * fe_0 + ta_xxx_xyyyz_0[i] * pa_y[i] - ta_xxx_xyyyz_1[i] * pc_y[i];

        ta_xxxy_xyyzz_0[i] = 2.0 * ta_xxx_xyzz_0[i] * fe_0 - 2.0 * ta_xxx_xyzz_1[i] * fe_0 + ta_xxx_xyyzz_0[i] * pa_y[i] - ta_xxx_xyyzz_1[i] * pc_y[i];

        ta_xxxy_xyzzz_0[i] = ta_xxx_xzzz_0[i] * fe_0 - ta_xxx_xzzz_1[i] * fe_0 + ta_xxx_xyzzz_0[i] * pa_y[i] - ta_xxx_xyzzz_1[i] * pc_y[i];

        ta_xxxy_xzzzz_0[i] = ta_xxx_xzzzz_0[i] * pa_y[i] - ta_xxx_xzzzz_1[i] * pc_y[i];

        ta_xxxy_yyyyy_0[i] = 2.0 * ta_xy_yyyyy_0[i] * fe_0 - 2.0 * ta_xy_yyyyy_1[i] * fe_0 + ta_xxy_yyyyy_0[i] * pa_x[i] - ta_xxy_yyyyy_1[i] * pc_x[i];

        ta_xxxy_yyyyz_0[i] = 2.0 * ta_xy_yyyyz_0[i] * fe_0 - 2.0 * ta_xy_yyyyz_1[i] * fe_0 + ta_xxy_yyyyz_0[i] * pa_x[i] - ta_xxy_yyyyz_1[i] * pc_x[i];

        ta_xxxy_yyyzz_0[i] = 2.0 * ta_xy_yyyzz_0[i] * fe_0 - 2.0 * ta_xy_yyyzz_1[i] * fe_0 + ta_xxy_yyyzz_0[i] * pa_x[i] - ta_xxy_yyyzz_1[i] * pc_x[i];

        ta_xxxy_yyzzz_0[i] = 2.0 * ta_xy_yyzzz_0[i] * fe_0 - 2.0 * ta_xy_yyzzz_1[i] * fe_0 + ta_xxy_yyzzz_0[i] * pa_x[i] - ta_xxy_yyzzz_1[i] * pc_x[i];

        ta_xxxy_yzzzz_0[i] = 2.0 * ta_xy_yzzzz_0[i] * fe_0 - 2.0 * ta_xy_yzzzz_1[i] * fe_0 + ta_xxy_yzzzz_0[i] * pa_x[i] - ta_xxy_yzzzz_1[i] * pc_x[i];

        ta_xxxy_zzzzz_0[i] = ta_xxx_zzzzz_0[i] * pa_y[i] - ta_xxx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : GH

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

    auto ta_xxxz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 57);

    auto ta_xxxz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 58);

    auto ta_xxxz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 59);

    auto ta_xxxz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 60);

    auto ta_xxxz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 61);

    auto ta_xxxz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 62);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xxx_xxxx_0, ta_xxx_xxxx_1, ta_xxx_xxxxx_0, ta_xxx_xxxxx_1, ta_xxx_xxxxy_0, ta_xxx_xxxxy_1, ta_xxx_xxxxz_0, ta_xxx_xxxxz_1, ta_xxx_xxxy_0, ta_xxx_xxxy_1, ta_xxx_xxxyy_0, ta_xxx_xxxyy_1, ta_xxx_xxxyz_0, ta_xxx_xxxyz_1, ta_xxx_xxxz_0, ta_xxx_xxxz_1, ta_xxx_xxxzz_0, ta_xxx_xxxzz_1, ta_xxx_xxyy_0, ta_xxx_xxyy_1, ta_xxx_xxyyy_0, ta_xxx_xxyyy_1, ta_xxx_xxyyz_0, ta_xxx_xxyyz_1, ta_xxx_xxyz_0, ta_xxx_xxyz_1, ta_xxx_xxyzz_0, ta_xxx_xxyzz_1, ta_xxx_xxzz_0, ta_xxx_xxzz_1, ta_xxx_xxzzz_0, ta_xxx_xxzzz_1, ta_xxx_xyyy_0, ta_xxx_xyyy_1, ta_xxx_xyyyy_0, ta_xxx_xyyyy_1, ta_xxx_xyyyz_0, ta_xxx_xyyyz_1, ta_xxx_xyyz_0, ta_xxx_xyyz_1, ta_xxx_xyyzz_0, ta_xxx_xyyzz_1, ta_xxx_xyzz_0, ta_xxx_xyzz_1, ta_xxx_xyzzz_0, ta_xxx_xyzzz_1, ta_xxx_xzzz_0, ta_xxx_xzzz_1, ta_xxx_xzzzz_0, ta_xxx_xzzzz_1, ta_xxx_yyyyy_0, ta_xxx_yyyyy_1, ta_xxxz_xxxxx_0, ta_xxxz_xxxxy_0, ta_xxxz_xxxxz_0, ta_xxxz_xxxyy_0, ta_xxxz_xxxyz_0, ta_xxxz_xxxzz_0, ta_xxxz_xxyyy_0, ta_xxxz_xxyyz_0, ta_xxxz_xxyzz_0, ta_xxxz_xxzzz_0, ta_xxxz_xyyyy_0, ta_xxxz_xyyyz_0, ta_xxxz_xyyzz_0, ta_xxxz_xyzzz_0, ta_xxxz_xzzzz_0, ta_xxxz_yyyyy_0, ta_xxxz_yyyyz_0, ta_xxxz_yyyzz_0, ta_xxxz_yyzzz_0, ta_xxxz_yzzzz_0, ta_xxxz_zzzzz_0, ta_xxz_yyyyz_0, ta_xxz_yyyyz_1, ta_xxz_yyyzz_0, ta_xxz_yyyzz_1, ta_xxz_yyzzz_0, ta_xxz_yyzzz_1, ta_xxz_yzzzz_0, ta_xxz_yzzzz_1, ta_xxz_zzzzz_0, ta_xxz_zzzzz_1, ta_xz_yyyyz_0, ta_xz_yyyyz_1, ta_xz_yyyzz_0, ta_xz_yyyzz_1, ta_xz_yyzzz_0, ta_xz_yyzzz_1, ta_xz_yzzzz_0, ta_xz_yzzzz_1, ta_xz_zzzzz_0, ta_xz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxz_xxxxx_0[i] = ta_xxx_xxxxx_0[i] * pa_z[i] - ta_xxx_xxxxx_1[i] * pc_z[i];

        ta_xxxz_xxxxy_0[i] = ta_xxx_xxxxy_0[i] * pa_z[i] - ta_xxx_xxxxy_1[i] * pc_z[i];

        ta_xxxz_xxxxz_0[i] = ta_xxx_xxxx_0[i] * fe_0 - ta_xxx_xxxx_1[i] * fe_0 + ta_xxx_xxxxz_0[i] * pa_z[i] - ta_xxx_xxxxz_1[i] * pc_z[i];

        ta_xxxz_xxxyy_0[i] = ta_xxx_xxxyy_0[i] * pa_z[i] - ta_xxx_xxxyy_1[i] * pc_z[i];

        ta_xxxz_xxxyz_0[i] = ta_xxx_xxxy_0[i] * fe_0 - ta_xxx_xxxy_1[i] * fe_0 + ta_xxx_xxxyz_0[i] * pa_z[i] - ta_xxx_xxxyz_1[i] * pc_z[i];

        ta_xxxz_xxxzz_0[i] = 2.0 * ta_xxx_xxxz_0[i] * fe_0 - 2.0 * ta_xxx_xxxz_1[i] * fe_0 + ta_xxx_xxxzz_0[i] * pa_z[i] - ta_xxx_xxxzz_1[i] * pc_z[i];

        ta_xxxz_xxyyy_0[i] = ta_xxx_xxyyy_0[i] * pa_z[i] - ta_xxx_xxyyy_1[i] * pc_z[i];

        ta_xxxz_xxyyz_0[i] = ta_xxx_xxyy_0[i] * fe_0 - ta_xxx_xxyy_1[i] * fe_0 + ta_xxx_xxyyz_0[i] * pa_z[i] - ta_xxx_xxyyz_1[i] * pc_z[i];

        ta_xxxz_xxyzz_0[i] = 2.0 * ta_xxx_xxyz_0[i] * fe_0 - 2.0 * ta_xxx_xxyz_1[i] * fe_0 + ta_xxx_xxyzz_0[i] * pa_z[i] - ta_xxx_xxyzz_1[i] * pc_z[i];

        ta_xxxz_xxzzz_0[i] = 3.0 * ta_xxx_xxzz_0[i] * fe_0 - 3.0 * ta_xxx_xxzz_1[i] * fe_0 + ta_xxx_xxzzz_0[i] * pa_z[i] - ta_xxx_xxzzz_1[i] * pc_z[i];

        ta_xxxz_xyyyy_0[i] = ta_xxx_xyyyy_0[i] * pa_z[i] - ta_xxx_xyyyy_1[i] * pc_z[i];

        ta_xxxz_xyyyz_0[i] = ta_xxx_xyyy_0[i] * fe_0 - ta_xxx_xyyy_1[i] * fe_0 + ta_xxx_xyyyz_0[i] * pa_z[i] - ta_xxx_xyyyz_1[i] * pc_z[i];

        ta_xxxz_xyyzz_0[i] = 2.0 * ta_xxx_xyyz_0[i] * fe_0 - 2.0 * ta_xxx_xyyz_1[i] * fe_0 + ta_xxx_xyyzz_0[i] * pa_z[i] - ta_xxx_xyyzz_1[i] * pc_z[i];

        ta_xxxz_xyzzz_0[i] = 3.0 * ta_xxx_xyzz_0[i] * fe_0 - 3.0 * ta_xxx_xyzz_1[i] * fe_0 + ta_xxx_xyzzz_0[i] * pa_z[i] - ta_xxx_xyzzz_1[i] * pc_z[i];

        ta_xxxz_xzzzz_0[i] = 4.0 * ta_xxx_xzzz_0[i] * fe_0 - 4.0 * ta_xxx_xzzz_1[i] * fe_0 + ta_xxx_xzzzz_0[i] * pa_z[i] - ta_xxx_xzzzz_1[i] * pc_z[i];

        ta_xxxz_yyyyy_0[i] = ta_xxx_yyyyy_0[i] * pa_z[i] - ta_xxx_yyyyy_1[i] * pc_z[i];

        ta_xxxz_yyyyz_0[i] = 2.0 * ta_xz_yyyyz_0[i] * fe_0 - 2.0 * ta_xz_yyyyz_1[i] * fe_0 + ta_xxz_yyyyz_0[i] * pa_x[i] - ta_xxz_yyyyz_1[i] * pc_x[i];

        ta_xxxz_yyyzz_0[i] = 2.0 * ta_xz_yyyzz_0[i] * fe_0 - 2.0 * ta_xz_yyyzz_1[i] * fe_0 + ta_xxz_yyyzz_0[i] * pa_x[i] - ta_xxz_yyyzz_1[i] * pc_x[i];

        ta_xxxz_yyzzz_0[i] = 2.0 * ta_xz_yyzzz_0[i] * fe_0 - 2.0 * ta_xz_yyzzz_1[i] * fe_0 + ta_xxz_yyzzz_0[i] * pa_x[i] - ta_xxz_yyzzz_1[i] * pc_x[i];

        ta_xxxz_yzzzz_0[i] = 2.0 * ta_xz_yzzzz_0[i] * fe_0 - 2.0 * ta_xz_yzzzz_1[i] * fe_0 + ta_xxz_yzzzz_0[i] * pa_x[i] - ta_xxz_yzzzz_1[i] * pc_x[i];

        ta_xxxz_zzzzz_0[i] = 2.0 * ta_xz_zzzzz_0[i] * fe_0 - 2.0 * ta_xz_zzzzz_1[i] * fe_0 + ta_xxz_zzzzz_0[i] * pa_x[i] - ta_xxz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 63-84 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xxy_xxxxx_0, ta_xxy_xxxxx_1, ta_xxy_xxxxz_0, ta_xxy_xxxxz_1, ta_xxy_xxxzz_0, ta_xxy_xxxzz_1, ta_xxy_xxzzz_0, ta_xxy_xxzzz_1, ta_xxy_xzzzz_0, ta_xxy_xzzzz_1, ta_xxyy_xxxxx_0, ta_xxyy_xxxxy_0, ta_xxyy_xxxxz_0, ta_xxyy_xxxyy_0, ta_xxyy_xxxyz_0, ta_xxyy_xxxzz_0, ta_xxyy_xxyyy_0, ta_xxyy_xxyyz_0, ta_xxyy_xxyzz_0, ta_xxyy_xxzzz_0, ta_xxyy_xyyyy_0, ta_xxyy_xyyyz_0, ta_xxyy_xyyzz_0, ta_xxyy_xyzzz_0, ta_xxyy_xzzzz_0, ta_xxyy_yyyyy_0, ta_xxyy_yyyyz_0, ta_xxyy_yyyzz_0, ta_xxyy_yyzzz_0, ta_xxyy_yzzzz_0, ta_xxyy_zzzzz_0, ta_xyy_xxxxy_0, ta_xyy_xxxxy_1, ta_xyy_xxxy_0, ta_xyy_xxxy_1, ta_xyy_xxxyy_0, ta_xyy_xxxyy_1, ta_xyy_xxxyz_0, ta_xyy_xxxyz_1, ta_xyy_xxyy_0, ta_xyy_xxyy_1, ta_xyy_xxyyy_0, ta_xyy_xxyyy_1, ta_xyy_xxyyz_0, ta_xyy_xxyyz_1, ta_xyy_xxyz_0, ta_xyy_xxyz_1, ta_xyy_xxyzz_0, ta_xyy_xxyzz_1, ta_xyy_xyyy_0, ta_xyy_xyyy_1, ta_xyy_xyyyy_0, ta_xyy_xyyyy_1, ta_xyy_xyyyz_0, ta_xyy_xyyyz_1, ta_xyy_xyyz_0, ta_xyy_xyyz_1, ta_xyy_xyyzz_0, ta_xyy_xyyzz_1, ta_xyy_xyzz_0, ta_xyy_xyzz_1, ta_xyy_xyzzz_0, ta_xyy_xyzzz_1, ta_xyy_yyyy_0, ta_xyy_yyyy_1, ta_xyy_yyyyy_0, ta_xyy_yyyyy_1, ta_xyy_yyyyz_0, ta_xyy_yyyyz_1, ta_xyy_yyyz_0, ta_xyy_yyyz_1, ta_xyy_yyyzz_0, ta_xyy_yyyzz_1, ta_xyy_yyzz_0, ta_xyy_yyzz_1, ta_xyy_yyzzz_0, ta_xyy_yyzzz_1, ta_xyy_yzzz_0, ta_xyy_yzzz_1, ta_xyy_yzzzz_0, ta_xyy_yzzzz_1, ta_xyy_zzzzz_0, ta_xyy_zzzzz_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_zzzzz_0, ta_yy_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyy_xxxxx_0[i] = ta_xx_xxxxx_0[i] * fe_0 - ta_xx_xxxxx_1[i] * fe_0 + ta_xxy_xxxxx_0[i] * pa_y[i] - ta_xxy_xxxxx_1[i] * pc_y[i];

        ta_xxyy_xxxxy_0[i] = ta_yy_xxxxy_0[i] * fe_0 - ta_yy_xxxxy_1[i] * fe_0 + 4.0 * ta_xyy_xxxy_0[i] * fe_0 - 4.0 * ta_xyy_xxxy_1[i] * fe_0 + ta_xyy_xxxxy_0[i] * pa_x[i] - ta_xyy_xxxxy_1[i] * pc_x[i];

        ta_xxyy_xxxxz_0[i] = ta_xx_xxxxz_0[i] * fe_0 - ta_xx_xxxxz_1[i] * fe_0 + ta_xxy_xxxxz_0[i] * pa_y[i] - ta_xxy_xxxxz_1[i] * pc_y[i];

        ta_xxyy_xxxyy_0[i] = ta_yy_xxxyy_0[i] * fe_0 - ta_yy_xxxyy_1[i] * fe_0 + 3.0 * ta_xyy_xxyy_0[i] * fe_0 - 3.0 * ta_xyy_xxyy_1[i] * fe_0 + ta_xyy_xxxyy_0[i] * pa_x[i] - ta_xyy_xxxyy_1[i] * pc_x[i];

        ta_xxyy_xxxyz_0[i] = ta_yy_xxxyz_0[i] * fe_0 - ta_yy_xxxyz_1[i] * fe_0 + 3.0 * ta_xyy_xxyz_0[i] * fe_0 - 3.0 * ta_xyy_xxyz_1[i] * fe_0 + ta_xyy_xxxyz_0[i] * pa_x[i] - ta_xyy_xxxyz_1[i] * pc_x[i];

        ta_xxyy_xxxzz_0[i] = ta_xx_xxxzz_0[i] * fe_0 - ta_xx_xxxzz_1[i] * fe_0 + ta_xxy_xxxzz_0[i] * pa_y[i] - ta_xxy_xxxzz_1[i] * pc_y[i];

        ta_xxyy_xxyyy_0[i] = ta_yy_xxyyy_0[i] * fe_0 - ta_yy_xxyyy_1[i] * fe_0 + 2.0 * ta_xyy_xyyy_0[i] * fe_0 - 2.0 * ta_xyy_xyyy_1[i] * fe_0 + ta_xyy_xxyyy_0[i] * pa_x[i] - ta_xyy_xxyyy_1[i] * pc_x[i];

        ta_xxyy_xxyyz_0[i] = ta_yy_xxyyz_0[i] * fe_0 - ta_yy_xxyyz_1[i] * fe_0 + 2.0 * ta_xyy_xyyz_0[i] * fe_0 - 2.0 * ta_xyy_xyyz_1[i] * fe_0 + ta_xyy_xxyyz_0[i] * pa_x[i] - ta_xyy_xxyyz_1[i] * pc_x[i];

        ta_xxyy_xxyzz_0[i] = ta_yy_xxyzz_0[i] * fe_0 - ta_yy_xxyzz_1[i] * fe_0 + 2.0 * ta_xyy_xyzz_0[i] * fe_0 - 2.0 * ta_xyy_xyzz_1[i] * fe_0 + ta_xyy_xxyzz_0[i] * pa_x[i] - ta_xyy_xxyzz_1[i] * pc_x[i];

        ta_xxyy_xxzzz_0[i] = ta_xx_xxzzz_0[i] * fe_0 - ta_xx_xxzzz_1[i] * fe_0 + ta_xxy_xxzzz_0[i] * pa_y[i] - ta_xxy_xxzzz_1[i] * pc_y[i];

        ta_xxyy_xyyyy_0[i] = ta_yy_xyyyy_0[i] * fe_0 - ta_yy_xyyyy_1[i] * fe_0 + ta_xyy_yyyy_0[i] * fe_0 - ta_xyy_yyyy_1[i] * fe_0 + ta_xyy_xyyyy_0[i] * pa_x[i] - ta_xyy_xyyyy_1[i] * pc_x[i];

        ta_xxyy_xyyyz_0[i] = ta_yy_xyyyz_0[i] * fe_0 - ta_yy_xyyyz_1[i] * fe_0 + ta_xyy_yyyz_0[i] * fe_0 - ta_xyy_yyyz_1[i] * fe_0 + ta_xyy_xyyyz_0[i] * pa_x[i] - ta_xyy_xyyyz_1[i] * pc_x[i];

        ta_xxyy_xyyzz_0[i] = ta_yy_xyyzz_0[i] * fe_0 - ta_yy_xyyzz_1[i] * fe_0 + ta_xyy_yyzz_0[i] * fe_0 - ta_xyy_yyzz_1[i] * fe_0 + ta_xyy_xyyzz_0[i] * pa_x[i] - ta_xyy_xyyzz_1[i] * pc_x[i];

        ta_xxyy_xyzzz_0[i] = ta_yy_xyzzz_0[i] * fe_0 - ta_yy_xyzzz_1[i] * fe_0 + ta_xyy_yzzz_0[i] * fe_0 - ta_xyy_yzzz_1[i] * fe_0 + ta_xyy_xyzzz_0[i] * pa_x[i] - ta_xyy_xyzzz_1[i] * pc_x[i];

        ta_xxyy_xzzzz_0[i] = ta_xx_xzzzz_0[i] * fe_0 - ta_xx_xzzzz_1[i] * fe_0 + ta_xxy_xzzzz_0[i] * pa_y[i] - ta_xxy_xzzzz_1[i] * pc_y[i];

        ta_xxyy_yyyyy_0[i] = ta_yy_yyyyy_0[i] * fe_0 - ta_yy_yyyyy_1[i] * fe_0 + ta_xyy_yyyyy_0[i] * pa_x[i] - ta_xyy_yyyyy_1[i] * pc_x[i];

        ta_xxyy_yyyyz_0[i] = ta_yy_yyyyz_0[i] * fe_0 - ta_yy_yyyyz_1[i] * fe_0 + ta_xyy_yyyyz_0[i] * pa_x[i] - ta_xyy_yyyyz_1[i] * pc_x[i];

        ta_xxyy_yyyzz_0[i] = ta_yy_yyyzz_0[i] * fe_0 - ta_yy_yyyzz_1[i] * fe_0 + ta_xyy_yyyzz_0[i] * pa_x[i] - ta_xyy_yyyzz_1[i] * pc_x[i];

        ta_xxyy_yyzzz_0[i] = ta_yy_yyzzz_0[i] * fe_0 - ta_yy_yyzzz_1[i] * fe_0 + ta_xyy_yyzzz_0[i] * pa_x[i] - ta_xyy_yyzzz_1[i] * pc_x[i];

        ta_xxyy_yzzzz_0[i] = ta_yy_yzzzz_0[i] * fe_0 - ta_yy_yzzzz_1[i] * fe_0 + ta_xyy_yzzzz_0[i] * pa_x[i] - ta_xyy_yzzzz_1[i] * pc_x[i];

        ta_xxyy_zzzzz_0[i] = ta_yy_zzzzz_0[i] * fe_0 - ta_yy_zzzzz_1[i] * fe_0 + ta_xyy_zzzzz_0[i] * pa_x[i] - ta_xyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 84-105 components of targeted buffer : GH

    auto ta_xxyz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 84);

    auto ta_xxyz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 85);

    auto ta_xxyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 86);

    auto ta_xxyz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 87);

    auto ta_xxyz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 88);

    auto ta_xxyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 89);

    auto ta_xxyz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 90);

    auto ta_xxyz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 91);

    auto ta_xxyz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 92);

    auto ta_xxyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 93);

    auto ta_xxyz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 94);

    auto ta_xxyz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 95);

    auto ta_xxyz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 96);

    auto ta_xxyz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 97);

    auto ta_xxyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 98);

    auto ta_xxyz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 99);

    auto ta_xxyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 100);

    auto ta_xxyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 101);

    auto ta_xxyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 102);

    auto ta_xxyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 103);

    auto ta_xxyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 104);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxy_xxxxy_0, ta_xxy_xxxxy_1, ta_xxy_xxxyy_0, ta_xxy_xxxyy_1, ta_xxy_xxyyy_0, ta_xxy_xxyyy_1, ta_xxy_xyyyy_0, ta_xxy_xyyyy_1, ta_xxy_yyyyy_0, ta_xxy_yyyyy_1, ta_xxyz_xxxxx_0, ta_xxyz_xxxxy_0, ta_xxyz_xxxxz_0, ta_xxyz_xxxyy_0, ta_xxyz_xxxyz_0, ta_xxyz_xxxzz_0, ta_xxyz_xxyyy_0, ta_xxyz_xxyyz_0, ta_xxyz_xxyzz_0, ta_xxyz_xxzzz_0, ta_xxyz_xyyyy_0, ta_xxyz_xyyyz_0, ta_xxyz_xyyzz_0, ta_xxyz_xyzzz_0, ta_xxyz_xzzzz_0, ta_xxyz_yyyyy_0, ta_xxyz_yyyyz_0, ta_xxyz_yyyzz_0, ta_xxyz_yyzzz_0, ta_xxyz_yzzzz_0, ta_xxyz_zzzzz_0, ta_xxz_xxxxx_0, ta_xxz_xxxxx_1, ta_xxz_xxxxz_0, ta_xxz_xxxxz_1, ta_xxz_xxxyz_0, ta_xxz_xxxyz_1, ta_xxz_xxxz_0, ta_xxz_xxxz_1, ta_xxz_xxxzz_0, ta_xxz_xxxzz_1, ta_xxz_xxyyz_0, ta_xxz_xxyyz_1, ta_xxz_xxyz_0, ta_xxz_xxyz_1, ta_xxz_xxyzz_0, ta_xxz_xxyzz_1, ta_xxz_xxzz_0, ta_xxz_xxzz_1, ta_xxz_xxzzz_0, ta_xxz_xxzzz_1, ta_xxz_xyyyz_0, ta_xxz_xyyyz_1, ta_xxz_xyyz_0, ta_xxz_xyyz_1, ta_xxz_xyyzz_0, ta_xxz_xyyzz_1, ta_xxz_xyzz_0, ta_xxz_xyzz_1, ta_xxz_xyzzz_0, ta_xxz_xyzzz_1, ta_xxz_xzzz_0, ta_xxz_xzzz_1, ta_xxz_xzzzz_0, ta_xxz_xzzzz_1, ta_xxz_zzzzz_0, ta_xxz_zzzzz_1, ta_xyz_yyyyz_0, ta_xyz_yyyyz_1, ta_xyz_yyyzz_0, ta_xyz_yyyzz_1, ta_xyz_yyzzz_0, ta_xyz_yyzzz_1, ta_xyz_yzzzz_0, ta_xyz_yzzzz_1, ta_yz_yyyyz_0, ta_yz_yyyyz_1, ta_yz_yyyzz_0, ta_yz_yyyzz_1, ta_yz_yyzzz_0, ta_yz_yyzzz_1, ta_yz_yzzzz_0, ta_yz_yzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyz_xxxxx_0[i] = ta_xxz_xxxxx_0[i] * pa_y[i] - ta_xxz_xxxxx_1[i] * pc_y[i];

        ta_xxyz_xxxxy_0[i] = ta_xxy_xxxxy_0[i] * pa_z[i] - ta_xxy_xxxxy_1[i] * pc_z[i];

        ta_xxyz_xxxxz_0[i] = ta_xxz_xxxxz_0[i] * pa_y[i] - ta_xxz_xxxxz_1[i] * pc_y[i];

        ta_xxyz_xxxyy_0[i] = ta_xxy_xxxyy_0[i] * pa_z[i] - ta_xxy_xxxyy_1[i] * pc_z[i];

        ta_xxyz_xxxyz_0[i] = ta_xxz_xxxz_0[i] * fe_0 - ta_xxz_xxxz_1[i] * fe_0 + ta_xxz_xxxyz_0[i] * pa_y[i] - ta_xxz_xxxyz_1[i] * pc_y[i];

        ta_xxyz_xxxzz_0[i] = ta_xxz_xxxzz_0[i] * pa_y[i] - ta_xxz_xxxzz_1[i] * pc_y[i];

        ta_xxyz_xxyyy_0[i] = ta_xxy_xxyyy_0[i] * pa_z[i] - ta_xxy_xxyyy_1[i] * pc_z[i];

        ta_xxyz_xxyyz_0[i] = 2.0 * ta_xxz_xxyz_0[i] * fe_0 - 2.0 * ta_xxz_xxyz_1[i] * fe_0 + ta_xxz_xxyyz_0[i] * pa_y[i] - ta_xxz_xxyyz_1[i] * pc_y[i];

        ta_xxyz_xxyzz_0[i] = ta_xxz_xxzz_0[i] * fe_0 - ta_xxz_xxzz_1[i] * fe_0 + ta_xxz_xxyzz_0[i] * pa_y[i] - ta_xxz_xxyzz_1[i] * pc_y[i];

        ta_xxyz_xxzzz_0[i] = ta_xxz_xxzzz_0[i] * pa_y[i] - ta_xxz_xxzzz_1[i] * pc_y[i];

        ta_xxyz_xyyyy_0[i] = ta_xxy_xyyyy_0[i] * pa_z[i] - ta_xxy_xyyyy_1[i] * pc_z[i];

        ta_xxyz_xyyyz_0[i] = 3.0 * ta_xxz_xyyz_0[i] * fe_0 - 3.0 * ta_xxz_xyyz_1[i] * fe_0 + ta_xxz_xyyyz_0[i] * pa_y[i] - ta_xxz_xyyyz_1[i] * pc_y[i];

        ta_xxyz_xyyzz_0[i] = 2.0 * ta_xxz_xyzz_0[i] * fe_0 - 2.0 * ta_xxz_xyzz_1[i] * fe_0 + ta_xxz_xyyzz_0[i] * pa_y[i] - ta_xxz_xyyzz_1[i] * pc_y[i];

        ta_xxyz_xyzzz_0[i] = ta_xxz_xzzz_0[i] * fe_0 - ta_xxz_xzzz_1[i] * fe_0 + ta_xxz_xyzzz_0[i] * pa_y[i] - ta_xxz_xyzzz_1[i] * pc_y[i];

        ta_xxyz_xzzzz_0[i] = ta_xxz_xzzzz_0[i] * pa_y[i] - ta_xxz_xzzzz_1[i] * pc_y[i];

        ta_xxyz_yyyyy_0[i] = ta_xxy_yyyyy_0[i] * pa_z[i] - ta_xxy_yyyyy_1[i] * pc_z[i];

        ta_xxyz_yyyyz_0[i] = ta_yz_yyyyz_0[i] * fe_0 - ta_yz_yyyyz_1[i] * fe_0 + ta_xyz_yyyyz_0[i] * pa_x[i] - ta_xyz_yyyyz_1[i] * pc_x[i];

        ta_xxyz_yyyzz_0[i] = ta_yz_yyyzz_0[i] * fe_0 - ta_yz_yyyzz_1[i] * fe_0 + ta_xyz_yyyzz_0[i] * pa_x[i] - ta_xyz_yyyzz_1[i] * pc_x[i];

        ta_xxyz_yyzzz_0[i] = ta_yz_yyzzz_0[i] * fe_0 - ta_yz_yyzzz_1[i] * fe_0 + ta_xyz_yyzzz_0[i] * pa_x[i] - ta_xyz_yyzzz_1[i] * pc_x[i];

        ta_xxyz_yzzzz_0[i] = ta_yz_yzzzz_0[i] * fe_0 - ta_yz_yzzzz_1[i] * fe_0 + ta_xyz_yzzzz_0[i] * pa_x[i] - ta_xyz_yzzzz_1[i] * pc_x[i];

        ta_xxyz_zzzzz_0[i] = ta_xxz_zzzzz_0[i] * pa_y[i] - ta_xxz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xxz_xxxxx_0, ta_xxz_xxxxx_1, ta_xxz_xxxxy_0, ta_xxz_xxxxy_1, ta_xxz_xxxyy_0, ta_xxz_xxxyy_1, ta_xxz_xxyyy_0, ta_xxz_xxyyy_1, ta_xxz_xyyyy_0, ta_xxz_xyyyy_1, ta_xxzz_xxxxx_0, ta_xxzz_xxxxy_0, ta_xxzz_xxxxz_0, ta_xxzz_xxxyy_0, ta_xxzz_xxxyz_0, ta_xxzz_xxxzz_0, ta_xxzz_xxyyy_0, ta_xxzz_xxyyz_0, ta_xxzz_xxyzz_0, ta_xxzz_xxzzz_0, ta_xxzz_xyyyy_0, ta_xxzz_xyyyz_0, ta_xxzz_xyyzz_0, ta_xxzz_xyzzz_0, ta_xxzz_xzzzz_0, ta_xxzz_yyyyy_0, ta_xxzz_yyyyz_0, ta_xxzz_yyyzz_0, ta_xxzz_yyzzz_0, ta_xxzz_yzzzz_0, ta_xxzz_zzzzz_0, ta_xzz_xxxxz_0, ta_xzz_xxxxz_1, ta_xzz_xxxyz_0, ta_xzz_xxxyz_1, ta_xzz_xxxz_0, ta_xzz_xxxz_1, ta_xzz_xxxzz_0, ta_xzz_xxxzz_1, ta_xzz_xxyyz_0, ta_xzz_xxyyz_1, ta_xzz_xxyz_0, ta_xzz_xxyz_1, ta_xzz_xxyzz_0, ta_xzz_xxyzz_1, ta_xzz_xxzz_0, ta_xzz_xxzz_1, ta_xzz_xxzzz_0, ta_xzz_xxzzz_1, ta_xzz_xyyyz_0, ta_xzz_xyyyz_1, ta_xzz_xyyz_0, ta_xzz_xyyz_1, ta_xzz_xyyzz_0, ta_xzz_xyyzz_1, ta_xzz_xyzz_0, ta_xzz_xyzz_1, ta_xzz_xyzzz_0, ta_xzz_xyzzz_1, ta_xzz_xzzz_0, ta_xzz_xzzz_1, ta_xzz_xzzzz_0, ta_xzz_xzzzz_1, ta_xzz_yyyyy_0, ta_xzz_yyyyy_1, ta_xzz_yyyyz_0, ta_xzz_yyyyz_1, ta_xzz_yyyz_0, ta_xzz_yyyz_1, ta_xzz_yyyzz_0, ta_xzz_yyyzz_1, ta_xzz_yyzz_0, ta_xzz_yyzz_1, ta_xzz_yyzzz_0, ta_xzz_yyzzz_1, ta_xzz_yzzz_0, ta_xzz_yzzz_1, ta_xzz_yzzzz_0, ta_xzz_yzzzz_1, ta_xzz_zzzz_0, ta_xzz_zzzz_1, ta_xzz_zzzzz_0, ta_xzz_zzzzz_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzz_xxxxx_0[i] = ta_xx_xxxxx_0[i] * fe_0 - ta_xx_xxxxx_1[i] * fe_0 + ta_xxz_xxxxx_0[i] * pa_z[i] - ta_xxz_xxxxx_1[i] * pc_z[i];

        ta_xxzz_xxxxy_0[i] = ta_xx_xxxxy_0[i] * fe_0 - ta_xx_xxxxy_1[i] * fe_0 + ta_xxz_xxxxy_0[i] * pa_z[i] - ta_xxz_xxxxy_1[i] * pc_z[i];

        ta_xxzz_xxxxz_0[i] = ta_zz_xxxxz_0[i] * fe_0 - ta_zz_xxxxz_1[i] * fe_0 + 4.0 * ta_xzz_xxxz_0[i] * fe_0 - 4.0 * ta_xzz_xxxz_1[i] * fe_0 + ta_xzz_xxxxz_0[i] * pa_x[i] - ta_xzz_xxxxz_1[i] * pc_x[i];

        ta_xxzz_xxxyy_0[i] = ta_xx_xxxyy_0[i] * fe_0 - ta_xx_xxxyy_1[i] * fe_0 + ta_xxz_xxxyy_0[i] * pa_z[i] - ta_xxz_xxxyy_1[i] * pc_z[i];

        ta_xxzz_xxxyz_0[i] = ta_zz_xxxyz_0[i] * fe_0 - ta_zz_xxxyz_1[i] * fe_0 + 3.0 * ta_xzz_xxyz_0[i] * fe_0 - 3.0 * ta_xzz_xxyz_1[i] * fe_0 + ta_xzz_xxxyz_0[i] * pa_x[i] - ta_xzz_xxxyz_1[i] * pc_x[i];

        ta_xxzz_xxxzz_0[i] = ta_zz_xxxzz_0[i] * fe_0 - ta_zz_xxxzz_1[i] * fe_0 + 3.0 * ta_xzz_xxzz_0[i] * fe_0 - 3.0 * ta_xzz_xxzz_1[i] * fe_0 + ta_xzz_xxxzz_0[i] * pa_x[i] - ta_xzz_xxxzz_1[i] * pc_x[i];

        ta_xxzz_xxyyy_0[i] = ta_xx_xxyyy_0[i] * fe_0 - ta_xx_xxyyy_1[i] * fe_0 + ta_xxz_xxyyy_0[i] * pa_z[i] - ta_xxz_xxyyy_1[i] * pc_z[i];

        ta_xxzz_xxyyz_0[i] = ta_zz_xxyyz_0[i] * fe_0 - ta_zz_xxyyz_1[i] * fe_0 + 2.0 * ta_xzz_xyyz_0[i] * fe_0 - 2.0 * ta_xzz_xyyz_1[i] * fe_0 + ta_xzz_xxyyz_0[i] * pa_x[i] - ta_xzz_xxyyz_1[i] * pc_x[i];

        ta_xxzz_xxyzz_0[i] = ta_zz_xxyzz_0[i] * fe_0 - ta_zz_xxyzz_1[i] * fe_0 + 2.0 * ta_xzz_xyzz_0[i] * fe_0 - 2.0 * ta_xzz_xyzz_1[i] * fe_0 + ta_xzz_xxyzz_0[i] * pa_x[i] - ta_xzz_xxyzz_1[i] * pc_x[i];

        ta_xxzz_xxzzz_0[i] = ta_zz_xxzzz_0[i] * fe_0 - ta_zz_xxzzz_1[i] * fe_0 + 2.0 * ta_xzz_xzzz_0[i] * fe_0 - 2.0 * ta_xzz_xzzz_1[i] * fe_0 + ta_xzz_xxzzz_0[i] * pa_x[i] - ta_xzz_xxzzz_1[i] * pc_x[i];

        ta_xxzz_xyyyy_0[i] = ta_xx_xyyyy_0[i] * fe_0 - ta_xx_xyyyy_1[i] * fe_0 + ta_xxz_xyyyy_0[i] * pa_z[i] - ta_xxz_xyyyy_1[i] * pc_z[i];

        ta_xxzz_xyyyz_0[i] = ta_zz_xyyyz_0[i] * fe_0 - ta_zz_xyyyz_1[i] * fe_0 + ta_xzz_yyyz_0[i] * fe_0 - ta_xzz_yyyz_1[i] * fe_0 + ta_xzz_xyyyz_0[i] * pa_x[i] - ta_xzz_xyyyz_1[i] * pc_x[i];

        ta_xxzz_xyyzz_0[i] = ta_zz_xyyzz_0[i] * fe_0 - ta_zz_xyyzz_1[i] * fe_0 + ta_xzz_yyzz_0[i] * fe_0 - ta_xzz_yyzz_1[i] * fe_0 + ta_xzz_xyyzz_0[i] * pa_x[i] - ta_xzz_xyyzz_1[i] * pc_x[i];

        ta_xxzz_xyzzz_0[i] = ta_zz_xyzzz_0[i] * fe_0 - ta_zz_xyzzz_1[i] * fe_0 + ta_xzz_yzzz_0[i] * fe_0 - ta_xzz_yzzz_1[i] * fe_0 + ta_xzz_xyzzz_0[i] * pa_x[i] - ta_xzz_xyzzz_1[i] * pc_x[i];

        ta_xxzz_xzzzz_0[i] = ta_zz_xzzzz_0[i] * fe_0 - ta_zz_xzzzz_1[i] * fe_0 + ta_xzz_zzzz_0[i] * fe_0 - ta_xzz_zzzz_1[i] * fe_0 + ta_xzz_xzzzz_0[i] * pa_x[i] - ta_xzz_xzzzz_1[i] * pc_x[i];

        ta_xxzz_yyyyy_0[i] = ta_zz_yyyyy_0[i] * fe_0 - ta_zz_yyyyy_1[i] * fe_0 + ta_xzz_yyyyy_0[i] * pa_x[i] - ta_xzz_yyyyy_1[i] * pc_x[i];

        ta_xxzz_yyyyz_0[i] = ta_zz_yyyyz_0[i] * fe_0 - ta_zz_yyyyz_1[i] * fe_0 + ta_xzz_yyyyz_0[i] * pa_x[i] - ta_xzz_yyyyz_1[i] * pc_x[i];

        ta_xxzz_yyyzz_0[i] = ta_zz_yyyzz_0[i] * fe_0 - ta_zz_yyyzz_1[i] * fe_0 + ta_xzz_yyyzz_0[i] * pa_x[i] - ta_xzz_yyyzz_1[i] * pc_x[i];

        ta_xxzz_yyzzz_0[i] = ta_zz_yyzzz_0[i] * fe_0 - ta_zz_yyzzz_1[i] * fe_0 + ta_xzz_yyzzz_0[i] * pa_x[i] - ta_xzz_yyzzz_1[i] * pc_x[i];

        ta_xxzz_yzzzz_0[i] = ta_zz_yzzzz_0[i] * fe_0 - ta_zz_yzzzz_1[i] * fe_0 + ta_xzz_yzzzz_0[i] * pa_x[i] - ta_xzz_yzzzz_1[i] * pc_x[i];

        ta_xxzz_zzzzz_0[i] = ta_zz_zzzzz_0[i] * fe_0 - ta_zz_zzzzz_1[i] * fe_0 + ta_xzz_zzzzz_0[i] * pa_x[i] - ta_xzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 126-147 components of targeted buffer : GH

    auto ta_xyyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 126);

    auto ta_xyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 127);

    auto ta_xyyy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 128);

    auto ta_xyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 129);

    auto ta_xyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 130);

    auto ta_xyyy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 131);

    auto ta_xyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 132);

    auto ta_xyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 133);

    auto ta_xyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 134);

    auto ta_xyyy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 135);

    auto ta_xyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 136);

    auto ta_xyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 137);

    auto ta_xyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 138);

    auto ta_xyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 139);

    auto ta_xyyy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 140);

    auto ta_xyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 141);

    auto ta_xyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 142);

    auto ta_xyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 143);

    auto ta_xyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 144);

    auto ta_xyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 145);

    auto ta_xyyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 146);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyyy_xxxxx_0, ta_xyyy_xxxxy_0, ta_xyyy_xxxxz_0, ta_xyyy_xxxyy_0, ta_xyyy_xxxyz_0, ta_xyyy_xxxzz_0, ta_xyyy_xxyyy_0, ta_xyyy_xxyyz_0, ta_xyyy_xxyzz_0, ta_xyyy_xxzzz_0, ta_xyyy_xyyyy_0, ta_xyyy_xyyyz_0, ta_xyyy_xyyzz_0, ta_xyyy_xyzzz_0, ta_xyyy_xzzzz_0, ta_xyyy_yyyyy_0, ta_xyyy_yyyyz_0, ta_xyyy_yyyzz_0, ta_xyyy_yyzzz_0, ta_xyyy_yzzzz_0, ta_xyyy_zzzzz_0, ta_yyy_xxxx_0, ta_yyy_xxxx_1, ta_yyy_xxxxx_0, ta_yyy_xxxxx_1, ta_yyy_xxxxy_0, ta_yyy_xxxxy_1, ta_yyy_xxxxz_0, ta_yyy_xxxxz_1, ta_yyy_xxxy_0, ta_yyy_xxxy_1, ta_yyy_xxxyy_0, ta_yyy_xxxyy_1, ta_yyy_xxxyz_0, ta_yyy_xxxyz_1, ta_yyy_xxxz_0, ta_yyy_xxxz_1, ta_yyy_xxxzz_0, ta_yyy_xxxzz_1, ta_yyy_xxyy_0, ta_yyy_xxyy_1, ta_yyy_xxyyy_0, ta_yyy_xxyyy_1, ta_yyy_xxyyz_0, ta_yyy_xxyyz_1, ta_yyy_xxyz_0, ta_yyy_xxyz_1, ta_yyy_xxyzz_0, ta_yyy_xxyzz_1, ta_yyy_xxzz_0, ta_yyy_xxzz_1, ta_yyy_xxzzz_0, ta_yyy_xxzzz_1, ta_yyy_xyyy_0, ta_yyy_xyyy_1, ta_yyy_xyyyy_0, ta_yyy_xyyyy_1, ta_yyy_xyyyz_0, ta_yyy_xyyyz_1, ta_yyy_xyyz_0, ta_yyy_xyyz_1, ta_yyy_xyyzz_0, ta_yyy_xyyzz_1, ta_yyy_xyzz_0, ta_yyy_xyzz_1, ta_yyy_xyzzz_0, ta_yyy_xyzzz_1, ta_yyy_xzzz_0, ta_yyy_xzzz_1, ta_yyy_xzzzz_0, ta_yyy_xzzzz_1, ta_yyy_yyyy_0, ta_yyy_yyyy_1, ta_yyy_yyyyy_0, ta_yyy_yyyyy_1, ta_yyy_yyyyz_0, ta_yyy_yyyyz_1, ta_yyy_yyyz_0, ta_yyy_yyyz_1, ta_yyy_yyyzz_0, ta_yyy_yyyzz_1, ta_yyy_yyzz_0, ta_yyy_yyzz_1, ta_yyy_yyzzz_0, ta_yyy_yyzzz_1, ta_yyy_yzzz_0, ta_yyy_yzzz_1, ta_yyy_yzzzz_0, ta_yyy_yzzzz_1, ta_yyy_zzzz_0, ta_yyy_zzzz_1, ta_yyy_zzzzz_0, ta_yyy_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyy_xxxxx_0[i] = 5.0 * ta_yyy_xxxx_0[i] * fe_0 - 5.0 * ta_yyy_xxxx_1[i] * fe_0 + ta_yyy_xxxxx_0[i] * pa_x[i] - ta_yyy_xxxxx_1[i] * pc_x[i];

        ta_xyyy_xxxxy_0[i] = 4.0 * ta_yyy_xxxy_0[i] * fe_0 - 4.0 * ta_yyy_xxxy_1[i] * fe_0 + ta_yyy_xxxxy_0[i] * pa_x[i] - ta_yyy_xxxxy_1[i] * pc_x[i];

        ta_xyyy_xxxxz_0[i] = 4.0 * ta_yyy_xxxz_0[i] * fe_0 - 4.0 * ta_yyy_xxxz_1[i] * fe_0 + ta_yyy_xxxxz_0[i] * pa_x[i] - ta_yyy_xxxxz_1[i] * pc_x[i];

        ta_xyyy_xxxyy_0[i] = 3.0 * ta_yyy_xxyy_0[i] * fe_0 - 3.0 * ta_yyy_xxyy_1[i] * fe_0 + ta_yyy_xxxyy_0[i] * pa_x[i] - ta_yyy_xxxyy_1[i] * pc_x[i];

        ta_xyyy_xxxyz_0[i] = 3.0 * ta_yyy_xxyz_0[i] * fe_0 - 3.0 * ta_yyy_xxyz_1[i] * fe_0 + ta_yyy_xxxyz_0[i] * pa_x[i] - ta_yyy_xxxyz_1[i] * pc_x[i];

        ta_xyyy_xxxzz_0[i] = 3.0 * ta_yyy_xxzz_0[i] * fe_0 - 3.0 * ta_yyy_xxzz_1[i] * fe_0 + ta_yyy_xxxzz_0[i] * pa_x[i] - ta_yyy_xxxzz_1[i] * pc_x[i];

        ta_xyyy_xxyyy_0[i] = 2.0 * ta_yyy_xyyy_0[i] * fe_0 - 2.0 * ta_yyy_xyyy_1[i] * fe_0 + ta_yyy_xxyyy_0[i] * pa_x[i] - ta_yyy_xxyyy_1[i] * pc_x[i];

        ta_xyyy_xxyyz_0[i] = 2.0 * ta_yyy_xyyz_0[i] * fe_0 - 2.0 * ta_yyy_xyyz_1[i] * fe_0 + ta_yyy_xxyyz_0[i] * pa_x[i] - ta_yyy_xxyyz_1[i] * pc_x[i];

        ta_xyyy_xxyzz_0[i] = 2.0 * ta_yyy_xyzz_0[i] * fe_0 - 2.0 * ta_yyy_xyzz_1[i] * fe_0 + ta_yyy_xxyzz_0[i] * pa_x[i] - ta_yyy_xxyzz_1[i] * pc_x[i];

        ta_xyyy_xxzzz_0[i] = 2.0 * ta_yyy_xzzz_0[i] * fe_0 - 2.0 * ta_yyy_xzzz_1[i] * fe_0 + ta_yyy_xxzzz_0[i] * pa_x[i] - ta_yyy_xxzzz_1[i] * pc_x[i];

        ta_xyyy_xyyyy_0[i] = ta_yyy_yyyy_0[i] * fe_0 - ta_yyy_yyyy_1[i] * fe_0 + ta_yyy_xyyyy_0[i] * pa_x[i] - ta_yyy_xyyyy_1[i] * pc_x[i];

        ta_xyyy_xyyyz_0[i] = ta_yyy_yyyz_0[i] * fe_0 - ta_yyy_yyyz_1[i] * fe_0 + ta_yyy_xyyyz_0[i] * pa_x[i] - ta_yyy_xyyyz_1[i] * pc_x[i];

        ta_xyyy_xyyzz_0[i] = ta_yyy_yyzz_0[i] * fe_0 - ta_yyy_yyzz_1[i] * fe_0 + ta_yyy_xyyzz_0[i] * pa_x[i] - ta_yyy_xyyzz_1[i] * pc_x[i];

        ta_xyyy_xyzzz_0[i] = ta_yyy_yzzz_0[i] * fe_0 - ta_yyy_yzzz_1[i] * fe_0 + ta_yyy_xyzzz_0[i] * pa_x[i] - ta_yyy_xyzzz_1[i] * pc_x[i];

        ta_xyyy_xzzzz_0[i] = ta_yyy_zzzz_0[i] * fe_0 - ta_yyy_zzzz_1[i] * fe_0 + ta_yyy_xzzzz_0[i] * pa_x[i] - ta_yyy_xzzzz_1[i] * pc_x[i];

        ta_xyyy_yyyyy_0[i] = ta_yyy_yyyyy_0[i] * pa_x[i] - ta_yyy_yyyyy_1[i] * pc_x[i];

        ta_xyyy_yyyyz_0[i] = ta_yyy_yyyyz_0[i] * pa_x[i] - ta_yyy_yyyyz_1[i] * pc_x[i];

        ta_xyyy_yyyzz_0[i] = ta_yyy_yyyzz_0[i] * pa_x[i] - ta_yyy_yyyzz_1[i] * pc_x[i];

        ta_xyyy_yyzzz_0[i] = ta_yyy_yyzzz_0[i] * pa_x[i] - ta_yyy_yyzzz_1[i] * pc_x[i];

        ta_xyyy_yzzzz_0[i] = ta_yyy_yzzzz_0[i] * pa_x[i] - ta_yyy_yzzzz_1[i] * pc_x[i];

        ta_xyyy_zzzzz_0[i] = ta_yyy_zzzzz_0[i] * pa_x[i] - ta_yyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 147-168 components of targeted buffer : GH

    auto ta_xyyz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 147);

    auto ta_xyyz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 148);

    auto ta_xyyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 149);

    auto ta_xyyz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 150);

    auto ta_xyyz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 151);

    auto ta_xyyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 152);

    auto ta_xyyz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 153);

    auto ta_xyyz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 154);

    auto ta_xyyz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 155);

    auto ta_xyyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 156);

    auto ta_xyyz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 157);

    auto ta_xyyz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 158);

    auto ta_xyyz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 159);

    auto ta_xyyz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 160);

    auto ta_xyyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 161);

    auto ta_xyyz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 162);

    auto ta_xyyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 163);

    auto ta_xyyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 164);

    auto ta_xyyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 165);

    auto ta_xyyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 166);

    auto ta_xyyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 167);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xyy_xxxxx_0, ta_xyy_xxxxx_1, ta_xyy_xxxxy_0, ta_xyy_xxxxy_1, ta_xyy_xxxyy_0, ta_xyy_xxxyy_1, ta_xyy_xxyyy_0, ta_xyy_xxyyy_1, ta_xyy_xyyyy_0, ta_xyy_xyyyy_1, ta_xyyz_xxxxx_0, ta_xyyz_xxxxy_0, ta_xyyz_xxxxz_0, ta_xyyz_xxxyy_0, ta_xyyz_xxxyz_0, ta_xyyz_xxxzz_0, ta_xyyz_xxyyy_0, ta_xyyz_xxyyz_0, ta_xyyz_xxyzz_0, ta_xyyz_xxzzz_0, ta_xyyz_xyyyy_0, ta_xyyz_xyyyz_0, ta_xyyz_xyyzz_0, ta_xyyz_xyzzz_0, ta_xyyz_xzzzz_0, ta_xyyz_yyyyy_0, ta_xyyz_yyyyz_0, ta_xyyz_yyyzz_0, ta_xyyz_yyzzz_0, ta_xyyz_yzzzz_0, ta_xyyz_zzzzz_0, ta_yyz_xxxxz_0, ta_yyz_xxxxz_1, ta_yyz_xxxyz_0, ta_yyz_xxxyz_1, ta_yyz_xxxz_0, ta_yyz_xxxz_1, ta_yyz_xxxzz_0, ta_yyz_xxxzz_1, ta_yyz_xxyyz_0, ta_yyz_xxyyz_1, ta_yyz_xxyz_0, ta_yyz_xxyz_1, ta_yyz_xxyzz_0, ta_yyz_xxyzz_1, ta_yyz_xxzz_0, ta_yyz_xxzz_1, ta_yyz_xxzzz_0, ta_yyz_xxzzz_1, ta_yyz_xyyyz_0, ta_yyz_xyyyz_1, ta_yyz_xyyz_0, ta_yyz_xyyz_1, ta_yyz_xyyzz_0, ta_yyz_xyyzz_1, ta_yyz_xyzz_0, ta_yyz_xyzz_1, ta_yyz_xyzzz_0, ta_yyz_xyzzz_1, ta_yyz_xzzz_0, ta_yyz_xzzz_1, ta_yyz_xzzzz_0, ta_yyz_xzzzz_1, ta_yyz_yyyyy_0, ta_yyz_yyyyy_1, ta_yyz_yyyyz_0, ta_yyz_yyyyz_1, ta_yyz_yyyz_0, ta_yyz_yyyz_1, ta_yyz_yyyzz_0, ta_yyz_yyyzz_1, ta_yyz_yyzz_0, ta_yyz_yyzz_1, ta_yyz_yyzzz_0, ta_yyz_yyzzz_1, ta_yyz_yzzz_0, ta_yyz_yzzz_1, ta_yyz_yzzzz_0, ta_yyz_yzzzz_1, ta_yyz_zzzz_0, ta_yyz_zzzz_1, ta_yyz_zzzzz_0, ta_yyz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyz_xxxxx_0[i] = ta_xyy_xxxxx_0[i] * pa_z[i] - ta_xyy_xxxxx_1[i] * pc_z[i];

        ta_xyyz_xxxxy_0[i] = ta_xyy_xxxxy_0[i] * pa_z[i] - ta_xyy_xxxxy_1[i] * pc_z[i];

        ta_xyyz_xxxxz_0[i] = 4.0 * ta_yyz_xxxz_0[i] * fe_0 - 4.0 * ta_yyz_xxxz_1[i] * fe_0 + ta_yyz_xxxxz_0[i] * pa_x[i] - ta_yyz_xxxxz_1[i] * pc_x[i];

        ta_xyyz_xxxyy_0[i] = ta_xyy_xxxyy_0[i] * pa_z[i] - ta_xyy_xxxyy_1[i] * pc_z[i];

        ta_xyyz_xxxyz_0[i] = 3.0 * ta_yyz_xxyz_0[i] * fe_0 - 3.0 * ta_yyz_xxyz_1[i] * fe_0 + ta_yyz_xxxyz_0[i] * pa_x[i] - ta_yyz_xxxyz_1[i] * pc_x[i];

        ta_xyyz_xxxzz_0[i] = 3.0 * ta_yyz_xxzz_0[i] * fe_0 - 3.0 * ta_yyz_xxzz_1[i] * fe_0 + ta_yyz_xxxzz_0[i] * pa_x[i] - ta_yyz_xxxzz_1[i] * pc_x[i];

        ta_xyyz_xxyyy_0[i] = ta_xyy_xxyyy_0[i] * pa_z[i] - ta_xyy_xxyyy_1[i] * pc_z[i];

        ta_xyyz_xxyyz_0[i] = 2.0 * ta_yyz_xyyz_0[i] * fe_0 - 2.0 * ta_yyz_xyyz_1[i] * fe_0 + ta_yyz_xxyyz_0[i] * pa_x[i] - ta_yyz_xxyyz_1[i] * pc_x[i];

        ta_xyyz_xxyzz_0[i] = 2.0 * ta_yyz_xyzz_0[i] * fe_0 - 2.0 * ta_yyz_xyzz_1[i] * fe_0 + ta_yyz_xxyzz_0[i] * pa_x[i] - ta_yyz_xxyzz_1[i] * pc_x[i];

        ta_xyyz_xxzzz_0[i] = 2.0 * ta_yyz_xzzz_0[i] * fe_0 - 2.0 * ta_yyz_xzzz_1[i] * fe_0 + ta_yyz_xxzzz_0[i] * pa_x[i] - ta_yyz_xxzzz_1[i] * pc_x[i];

        ta_xyyz_xyyyy_0[i] = ta_xyy_xyyyy_0[i] * pa_z[i] - ta_xyy_xyyyy_1[i] * pc_z[i];

        ta_xyyz_xyyyz_0[i] = ta_yyz_yyyz_0[i] * fe_0 - ta_yyz_yyyz_1[i] * fe_0 + ta_yyz_xyyyz_0[i] * pa_x[i] - ta_yyz_xyyyz_1[i] * pc_x[i];

        ta_xyyz_xyyzz_0[i] = ta_yyz_yyzz_0[i] * fe_0 - ta_yyz_yyzz_1[i] * fe_0 + ta_yyz_xyyzz_0[i] * pa_x[i] - ta_yyz_xyyzz_1[i] * pc_x[i];

        ta_xyyz_xyzzz_0[i] = ta_yyz_yzzz_0[i] * fe_0 - ta_yyz_yzzz_1[i] * fe_0 + ta_yyz_xyzzz_0[i] * pa_x[i] - ta_yyz_xyzzz_1[i] * pc_x[i];

        ta_xyyz_xzzzz_0[i] = ta_yyz_zzzz_0[i] * fe_0 - ta_yyz_zzzz_1[i] * fe_0 + ta_yyz_xzzzz_0[i] * pa_x[i] - ta_yyz_xzzzz_1[i] * pc_x[i];

        ta_xyyz_yyyyy_0[i] = ta_yyz_yyyyy_0[i] * pa_x[i] - ta_yyz_yyyyy_1[i] * pc_x[i];

        ta_xyyz_yyyyz_0[i] = ta_yyz_yyyyz_0[i] * pa_x[i] - ta_yyz_yyyyz_1[i] * pc_x[i];

        ta_xyyz_yyyzz_0[i] = ta_yyz_yyyzz_0[i] * pa_x[i] - ta_yyz_yyyzz_1[i] * pc_x[i];

        ta_xyyz_yyzzz_0[i] = ta_yyz_yyzzz_0[i] * pa_x[i] - ta_yyz_yyzzz_1[i] * pc_x[i];

        ta_xyyz_yzzzz_0[i] = ta_yyz_yzzzz_0[i] * pa_x[i] - ta_yyz_yzzzz_1[i] * pc_x[i];

        ta_xyyz_zzzzz_0[i] = ta_yyz_zzzzz_0[i] * pa_x[i] - ta_yyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 168-189 components of targeted buffer : GH

    auto ta_xyzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 168);

    auto ta_xyzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 169);

    auto ta_xyzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 170);

    auto ta_xyzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 171);

    auto ta_xyzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 172);

    auto ta_xyzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 173);

    auto ta_xyzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 174);

    auto ta_xyzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 175);

    auto ta_xyzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 176);

    auto ta_xyzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 177);

    auto ta_xyzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 178);

    auto ta_xyzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 179);

    auto ta_xyzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 180);

    auto ta_xyzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 181);

    auto ta_xyzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 182);

    auto ta_xyzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 183);

    auto ta_xyzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 184);

    auto ta_xyzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 185);

    auto ta_xyzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 186);

    auto ta_xyzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 187);

    auto ta_xyzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 188);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xyzz_xxxxx_0, ta_xyzz_xxxxy_0, ta_xyzz_xxxxz_0, ta_xyzz_xxxyy_0, ta_xyzz_xxxyz_0, ta_xyzz_xxxzz_0, ta_xyzz_xxyyy_0, ta_xyzz_xxyyz_0, ta_xyzz_xxyzz_0, ta_xyzz_xxzzz_0, ta_xyzz_xyyyy_0, ta_xyzz_xyyyz_0, ta_xyzz_xyyzz_0, ta_xyzz_xyzzz_0, ta_xyzz_xzzzz_0, ta_xyzz_yyyyy_0, ta_xyzz_yyyyz_0, ta_xyzz_yyyzz_0, ta_xyzz_yyzzz_0, ta_xyzz_yzzzz_0, ta_xyzz_zzzzz_0, ta_xzz_xxxxx_0, ta_xzz_xxxxx_1, ta_xzz_xxxxz_0, ta_xzz_xxxxz_1, ta_xzz_xxxzz_0, ta_xzz_xxxzz_1, ta_xzz_xxzzz_0, ta_xzz_xxzzz_1, ta_xzz_xzzzz_0, ta_xzz_xzzzz_1, ta_yzz_xxxxy_0, ta_yzz_xxxxy_1, ta_yzz_xxxy_0, ta_yzz_xxxy_1, ta_yzz_xxxyy_0, ta_yzz_xxxyy_1, ta_yzz_xxxyz_0, ta_yzz_xxxyz_1, ta_yzz_xxyy_0, ta_yzz_xxyy_1, ta_yzz_xxyyy_0, ta_yzz_xxyyy_1, ta_yzz_xxyyz_0, ta_yzz_xxyyz_1, ta_yzz_xxyz_0, ta_yzz_xxyz_1, ta_yzz_xxyzz_0, ta_yzz_xxyzz_1, ta_yzz_xyyy_0, ta_yzz_xyyy_1, ta_yzz_xyyyy_0, ta_yzz_xyyyy_1, ta_yzz_xyyyz_0, ta_yzz_xyyyz_1, ta_yzz_xyyz_0, ta_yzz_xyyz_1, ta_yzz_xyyzz_0, ta_yzz_xyyzz_1, ta_yzz_xyzz_0, ta_yzz_xyzz_1, ta_yzz_xyzzz_0, ta_yzz_xyzzz_1, ta_yzz_yyyy_0, ta_yzz_yyyy_1, ta_yzz_yyyyy_0, ta_yzz_yyyyy_1, ta_yzz_yyyyz_0, ta_yzz_yyyyz_1, ta_yzz_yyyz_0, ta_yzz_yyyz_1, ta_yzz_yyyzz_0, ta_yzz_yyyzz_1, ta_yzz_yyzz_0, ta_yzz_yyzz_1, ta_yzz_yyzzz_0, ta_yzz_yyzzz_1, ta_yzz_yzzz_0, ta_yzz_yzzz_1, ta_yzz_yzzzz_0, ta_yzz_yzzzz_1, ta_yzz_zzzzz_0, ta_yzz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzz_xxxxx_0[i] = ta_xzz_xxxxx_0[i] * pa_y[i] - ta_xzz_xxxxx_1[i] * pc_y[i];

        ta_xyzz_xxxxy_0[i] = 4.0 * ta_yzz_xxxy_0[i] * fe_0 - 4.0 * ta_yzz_xxxy_1[i] * fe_0 + ta_yzz_xxxxy_0[i] * pa_x[i] - ta_yzz_xxxxy_1[i] * pc_x[i];

        ta_xyzz_xxxxz_0[i] = ta_xzz_xxxxz_0[i] * pa_y[i] - ta_xzz_xxxxz_1[i] * pc_y[i];

        ta_xyzz_xxxyy_0[i] = 3.0 * ta_yzz_xxyy_0[i] * fe_0 - 3.0 * ta_yzz_xxyy_1[i] * fe_0 + ta_yzz_xxxyy_0[i] * pa_x[i] - ta_yzz_xxxyy_1[i] * pc_x[i];

        ta_xyzz_xxxyz_0[i] = 3.0 * ta_yzz_xxyz_0[i] * fe_0 - 3.0 * ta_yzz_xxyz_1[i] * fe_0 + ta_yzz_xxxyz_0[i] * pa_x[i] - ta_yzz_xxxyz_1[i] * pc_x[i];

        ta_xyzz_xxxzz_0[i] = ta_xzz_xxxzz_0[i] * pa_y[i] - ta_xzz_xxxzz_1[i] * pc_y[i];

        ta_xyzz_xxyyy_0[i] = 2.0 * ta_yzz_xyyy_0[i] * fe_0 - 2.0 * ta_yzz_xyyy_1[i] * fe_0 + ta_yzz_xxyyy_0[i] * pa_x[i] - ta_yzz_xxyyy_1[i] * pc_x[i];

        ta_xyzz_xxyyz_0[i] = 2.0 * ta_yzz_xyyz_0[i] * fe_0 - 2.0 * ta_yzz_xyyz_1[i] * fe_0 + ta_yzz_xxyyz_0[i] * pa_x[i] - ta_yzz_xxyyz_1[i] * pc_x[i];

        ta_xyzz_xxyzz_0[i] = 2.0 * ta_yzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzz_1[i] * fe_0 + ta_yzz_xxyzz_0[i] * pa_x[i] - ta_yzz_xxyzz_1[i] * pc_x[i];

        ta_xyzz_xxzzz_0[i] = ta_xzz_xxzzz_0[i] * pa_y[i] - ta_xzz_xxzzz_1[i] * pc_y[i];

        ta_xyzz_xyyyy_0[i] = ta_yzz_yyyy_0[i] * fe_0 - ta_yzz_yyyy_1[i] * fe_0 + ta_yzz_xyyyy_0[i] * pa_x[i] - ta_yzz_xyyyy_1[i] * pc_x[i];

        ta_xyzz_xyyyz_0[i] = ta_yzz_yyyz_0[i] * fe_0 - ta_yzz_yyyz_1[i] * fe_0 + ta_yzz_xyyyz_0[i] * pa_x[i] - ta_yzz_xyyyz_1[i] * pc_x[i];

        ta_xyzz_xyyzz_0[i] = ta_yzz_yyzz_0[i] * fe_0 - ta_yzz_yyzz_1[i] * fe_0 + ta_yzz_xyyzz_0[i] * pa_x[i] - ta_yzz_xyyzz_1[i] * pc_x[i];

        ta_xyzz_xyzzz_0[i] = ta_yzz_yzzz_0[i] * fe_0 - ta_yzz_yzzz_1[i] * fe_0 + ta_yzz_xyzzz_0[i] * pa_x[i] - ta_yzz_xyzzz_1[i] * pc_x[i];

        ta_xyzz_xzzzz_0[i] = ta_xzz_xzzzz_0[i] * pa_y[i] - ta_xzz_xzzzz_1[i] * pc_y[i];

        ta_xyzz_yyyyy_0[i] = ta_yzz_yyyyy_0[i] * pa_x[i] - ta_yzz_yyyyy_1[i] * pc_x[i];

        ta_xyzz_yyyyz_0[i] = ta_yzz_yyyyz_0[i] * pa_x[i] - ta_yzz_yyyyz_1[i] * pc_x[i];

        ta_xyzz_yyyzz_0[i] = ta_yzz_yyyzz_0[i] * pa_x[i] - ta_yzz_yyyzz_1[i] * pc_x[i];

        ta_xyzz_yyzzz_0[i] = ta_yzz_yyzzz_0[i] * pa_x[i] - ta_yzz_yyzzz_1[i] * pc_x[i];

        ta_xyzz_yzzzz_0[i] = ta_yzz_yzzzz_0[i] * pa_x[i] - ta_yzz_yzzzz_1[i] * pc_x[i];

        ta_xyzz_zzzzz_0[i] = ta_yzz_zzzzz_0[i] * pa_x[i] - ta_yzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 189-210 components of targeted buffer : GH

    auto ta_xzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 189);

    auto ta_xzzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 190);

    auto ta_xzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 191);

    auto ta_xzzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 192);

    auto ta_xzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 193);

    auto ta_xzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 194);

    auto ta_xzzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 195);

    auto ta_xzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 196);

    auto ta_xzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 197);

    auto ta_xzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 198);

    auto ta_xzzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 199);

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

    #pragma omp simd aligned(pa_x, pc_x, ta_xzzz_xxxxx_0, ta_xzzz_xxxxy_0, ta_xzzz_xxxxz_0, ta_xzzz_xxxyy_0, ta_xzzz_xxxyz_0, ta_xzzz_xxxzz_0, ta_xzzz_xxyyy_0, ta_xzzz_xxyyz_0, ta_xzzz_xxyzz_0, ta_xzzz_xxzzz_0, ta_xzzz_xyyyy_0, ta_xzzz_xyyyz_0, ta_xzzz_xyyzz_0, ta_xzzz_xyzzz_0, ta_xzzz_xzzzz_0, ta_xzzz_yyyyy_0, ta_xzzz_yyyyz_0, ta_xzzz_yyyzz_0, ta_xzzz_yyzzz_0, ta_xzzz_yzzzz_0, ta_xzzz_zzzzz_0, ta_zzz_xxxx_0, ta_zzz_xxxx_1, ta_zzz_xxxxx_0, ta_zzz_xxxxx_1, ta_zzz_xxxxy_0, ta_zzz_xxxxy_1, ta_zzz_xxxxz_0, ta_zzz_xxxxz_1, ta_zzz_xxxy_0, ta_zzz_xxxy_1, ta_zzz_xxxyy_0, ta_zzz_xxxyy_1, ta_zzz_xxxyz_0, ta_zzz_xxxyz_1, ta_zzz_xxxz_0, ta_zzz_xxxz_1, ta_zzz_xxxzz_0, ta_zzz_xxxzz_1, ta_zzz_xxyy_0, ta_zzz_xxyy_1, ta_zzz_xxyyy_0, ta_zzz_xxyyy_1, ta_zzz_xxyyz_0, ta_zzz_xxyyz_1, ta_zzz_xxyz_0, ta_zzz_xxyz_1, ta_zzz_xxyzz_0, ta_zzz_xxyzz_1, ta_zzz_xxzz_0, ta_zzz_xxzz_1, ta_zzz_xxzzz_0, ta_zzz_xxzzz_1, ta_zzz_xyyy_0, ta_zzz_xyyy_1, ta_zzz_xyyyy_0, ta_zzz_xyyyy_1, ta_zzz_xyyyz_0, ta_zzz_xyyyz_1, ta_zzz_xyyz_0, ta_zzz_xyyz_1, ta_zzz_xyyzz_0, ta_zzz_xyyzz_1, ta_zzz_xyzz_0, ta_zzz_xyzz_1, ta_zzz_xyzzz_0, ta_zzz_xyzzz_1, ta_zzz_xzzz_0, ta_zzz_xzzz_1, ta_zzz_xzzzz_0, ta_zzz_xzzzz_1, ta_zzz_yyyy_0, ta_zzz_yyyy_1, ta_zzz_yyyyy_0, ta_zzz_yyyyy_1, ta_zzz_yyyyz_0, ta_zzz_yyyyz_1, ta_zzz_yyyz_0, ta_zzz_yyyz_1, ta_zzz_yyyzz_0, ta_zzz_yyyzz_1, ta_zzz_yyzz_0, ta_zzz_yyzz_1, ta_zzz_yyzzz_0, ta_zzz_yyzzz_1, ta_zzz_yzzz_0, ta_zzz_yzzz_1, ta_zzz_yzzzz_0, ta_zzz_yzzzz_1, ta_zzz_zzzz_0, ta_zzz_zzzz_1, ta_zzz_zzzzz_0, ta_zzz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzz_xxxxx_0[i] = 5.0 * ta_zzz_xxxx_0[i] * fe_0 - 5.0 * ta_zzz_xxxx_1[i] * fe_0 + ta_zzz_xxxxx_0[i] * pa_x[i] - ta_zzz_xxxxx_1[i] * pc_x[i];

        ta_xzzz_xxxxy_0[i] = 4.0 * ta_zzz_xxxy_0[i] * fe_0 - 4.0 * ta_zzz_xxxy_1[i] * fe_0 + ta_zzz_xxxxy_0[i] * pa_x[i] - ta_zzz_xxxxy_1[i] * pc_x[i];

        ta_xzzz_xxxxz_0[i] = 4.0 * ta_zzz_xxxz_0[i] * fe_0 - 4.0 * ta_zzz_xxxz_1[i] * fe_0 + ta_zzz_xxxxz_0[i] * pa_x[i] - ta_zzz_xxxxz_1[i] * pc_x[i];

        ta_xzzz_xxxyy_0[i] = 3.0 * ta_zzz_xxyy_0[i] * fe_0 - 3.0 * ta_zzz_xxyy_1[i] * fe_0 + ta_zzz_xxxyy_0[i] * pa_x[i] - ta_zzz_xxxyy_1[i] * pc_x[i];

        ta_xzzz_xxxyz_0[i] = 3.0 * ta_zzz_xxyz_0[i] * fe_0 - 3.0 * ta_zzz_xxyz_1[i] * fe_0 + ta_zzz_xxxyz_0[i] * pa_x[i] - ta_zzz_xxxyz_1[i] * pc_x[i];

        ta_xzzz_xxxzz_0[i] = 3.0 * ta_zzz_xxzz_0[i] * fe_0 - 3.0 * ta_zzz_xxzz_1[i] * fe_0 + ta_zzz_xxxzz_0[i] * pa_x[i] - ta_zzz_xxxzz_1[i] * pc_x[i];

        ta_xzzz_xxyyy_0[i] = 2.0 * ta_zzz_xyyy_0[i] * fe_0 - 2.0 * ta_zzz_xyyy_1[i] * fe_0 + ta_zzz_xxyyy_0[i] * pa_x[i] - ta_zzz_xxyyy_1[i] * pc_x[i];

        ta_xzzz_xxyyz_0[i] = 2.0 * ta_zzz_xyyz_0[i] * fe_0 - 2.0 * ta_zzz_xyyz_1[i] * fe_0 + ta_zzz_xxyyz_0[i] * pa_x[i] - ta_zzz_xxyyz_1[i] * pc_x[i];

        ta_xzzz_xxyzz_0[i] = 2.0 * ta_zzz_xyzz_0[i] * fe_0 - 2.0 * ta_zzz_xyzz_1[i] * fe_0 + ta_zzz_xxyzz_0[i] * pa_x[i] - ta_zzz_xxyzz_1[i] * pc_x[i];

        ta_xzzz_xxzzz_0[i] = 2.0 * ta_zzz_xzzz_0[i] * fe_0 - 2.0 * ta_zzz_xzzz_1[i] * fe_0 + ta_zzz_xxzzz_0[i] * pa_x[i] - ta_zzz_xxzzz_1[i] * pc_x[i];

        ta_xzzz_xyyyy_0[i] = ta_zzz_yyyy_0[i] * fe_0 - ta_zzz_yyyy_1[i] * fe_0 + ta_zzz_xyyyy_0[i] * pa_x[i] - ta_zzz_xyyyy_1[i] * pc_x[i];

        ta_xzzz_xyyyz_0[i] = ta_zzz_yyyz_0[i] * fe_0 - ta_zzz_yyyz_1[i] * fe_0 + ta_zzz_xyyyz_0[i] * pa_x[i] - ta_zzz_xyyyz_1[i] * pc_x[i];

        ta_xzzz_xyyzz_0[i] = ta_zzz_yyzz_0[i] * fe_0 - ta_zzz_yyzz_1[i] * fe_0 + ta_zzz_xyyzz_0[i] * pa_x[i] - ta_zzz_xyyzz_1[i] * pc_x[i];

        ta_xzzz_xyzzz_0[i] = ta_zzz_yzzz_0[i] * fe_0 - ta_zzz_yzzz_1[i] * fe_0 + ta_zzz_xyzzz_0[i] * pa_x[i] - ta_zzz_xyzzz_1[i] * pc_x[i];

        ta_xzzz_xzzzz_0[i] = ta_zzz_zzzz_0[i] * fe_0 - ta_zzz_zzzz_1[i] * fe_0 + ta_zzz_xzzzz_0[i] * pa_x[i] - ta_zzz_xzzzz_1[i] * pc_x[i];

        ta_xzzz_yyyyy_0[i] = ta_zzz_yyyyy_0[i] * pa_x[i] - ta_zzz_yyyyy_1[i] * pc_x[i];

        ta_xzzz_yyyyz_0[i] = ta_zzz_yyyyz_0[i] * pa_x[i] - ta_zzz_yyyyz_1[i] * pc_x[i];

        ta_xzzz_yyyzz_0[i] = ta_zzz_yyyzz_0[i] * pa_x[i] - ta_zzz_yyyzz_1[i] * pc_x[i];

        ta_xzzz_yyzzz_0[i] = ta_zzz_yyzzz_0[i] * pa_x[i] - ta_zzz_yyzzz_1[i] * pc_x[i];

        ta_xzzz_yzzzz_0[i] = ta_zzz_yzzzz_0[i] * pa_x[i] - ta_zzz_yzzzz_1[i] * pc_x[i];

        ta_xzzz_zzzzz_0[i] = ta_zzz_zzzzz_0[i] * pa_x[i] - ta_zzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 210-231 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, pc_y, ta_yy_xxxxx_0, ta_yy_xxxxx_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxxz_0, ta_yy_xxxxz_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxxzz_0, ta_yy_xxxzz_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xxzzz_0, ta_yy_xxzzz_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_xzzzz_0, ta_yy_xzzzz_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_zzzzz_0, ta_yy_zzzzz_1, ta_yyy_xxxx_0, ta_yyy_xxxx_1, ta_yyy_xxxxx_0, ta_yyy_xxxxx_1, ta_yyy_xxxxy_0, ta_yyy_xxxxy_1, ta_yyy_xxxxz_0, ta_yyy_xxxxz_1, ta_yyy_xxxy_0, ta_yyy_xxxy_1, ta_yyy_xxxyy_0, ta_yyy_xxxyy_1, ta_yyy_xxxyz_0, ta_yyy_xxxyz_1, ta_yyy_xxxz_0, ta_yyy_xxxz_1, ta_yyy_xxxzz_0, ta_yyy_xxxzz_1, ta_yyy_xxyy_0, ta_yyy_xxyy_1, ta_yyy_xxyyy_0, ta_yyy_xxyyy_1, ta_yyy_xxyyz_0, ta_yyy_xxyyz_1, ta_yyy_xxyz_0, ta_yyy_xxyz_1, ta_yyy_xxyzz_0, ta_yyy_xxyzz_1, ta_yyy_xxzz_0, ta_yyy_xxzz_1, ta_yyy_xxzzz_0, ta_yyy_xxzzz_1, ta_yyy_xyyy_0, ta_yyy_xyyy_1, ta_yyy_xyyyy_0, ta_yyy_xyyyy_1, ta_yyy_xyyyz_0, ta_yyy_xyyyz_1, ta_yyy_xyyz_0, ta_yyy_xyyz_1, ta_yyy_xyyzz_0, ta_yyy_xyyzz_1, ta_yyy_xyzz_0, ta_yyy_xyzz_1, ta_yyy_xyzzz_0, ta_yyy_xyzzz_1, ta_yyy_xzzz_0, ta_yyy_xzzz_1, ta_yyy_xzzzz_0, ta_yyy_xzzzz_1, ta_yyy_yyyy_0, ta_yyy_yyyy_1, ta_yyy_yyyyy_0, ta_yyy_yyyyy_1, ta_yyy_yyyyz_0, ta_yyy_yyyyz_1, ta_yyy_yyyz_0, ta_yyy_yyyz_1, ta_yyy_yyyzz_0, ta_yyy_yyyzz_1, ta_yyy_yyzz_0, ta_yyy_yyzz_1, ta_yyy_yyzzz_0, ta_yyy_yyzzz_1, ta_yyy_yzzz_0, ta_yyy_yzzz_1, ta_yyy_yzzzz_0, ta_yyy_yzzzz_1, ta_yyy_zzzz_0, ta_yyy_zzzz_1, ta_yyy_zzzzz_0, ta_yyy_zzzzz_1, ta_yyyy_xxxxx_0, ta_yyyy_xxxxy_0, ta_yyyy_xxxxz_0, ta_yyyy_xxxyy_0, ta_yyyy_xxxyz_0, ta_yyyy_xxxzz_0, ta_yyyy_xxyyy_0, ta_yyyy_xxyyz_0, ta_yyyy_xxyzz_0, ta_yyyy_xxzzz_0, ta_yyyy_xyyyy_0, ta_yyyy_xyyyz_0, ta_yyyy_xyyzz_0, ta_yyyy_xyzzz_0, ta_yyyy_xzzzz_0, ta_yyyy_yyyyy_0, ta_yyyy_yyyyz_0, ta_yyyy_yyyzz_0, ta_yyyy_yyzzz_0, ta_yyyy_yzzzz_0, ta_yyyy_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyy_xxxxx_0[i] = 3.0 * ta_yy_xxxxx_0[i] * fe_0 - 3.0 * ta_yy_xxxxx_1[i] * fe_0 + ta_yyy_xxxxx_0[i] * pa_y[i] - ta_yyy_xxxxx_1[i] * pc_y[i];

        ta_yyyy_xxxxy_0[i] = 3.0 * ta_yy_xxxxy_0[i] * fe_0 - 3.0 * ta_yy_xxxxy_1[i] * fe_0 + ta_yyy_xxxx_0[i] * fe_0 - ta_yyy_xxxx_1[i] * fe_0 + ta_yyy_xxxxy_0[i] * pa_y[i] - ta_yyy_xxxxy_1[i] * pc_y[i];

        ta_yyyy_xxxxz_0[i] = 3.0 * ta_yy_xxxxz_0[i] * fe_0 - 3.0 * ta_yy_xxxxz_1[i] * fe_0 + ta_yyy_xxxxz_0[i] * pa_y[i] - ta_yyy_xxxxz_1[i] * pc_y[i];

        ta_yyyy_xxxyy_0[i] = 3.0 * ta_yy_xxxyy_0[i] * fe_0 - 3.0 * ta_yy_xxxyy_1[i] * fe_0 + 2.0 * ta_yyy_xxxy_0[i] * fe_0 - 2.0 * ta_yyy_xxxy_1[i] * fe_0 + ta_yyy_xxxyy_0[i] * pa_y[i] - ta_yyy_xxxyy_1[i] * pc_y[i];

        ta_yyyy_xxxyz_0[i] = 3.0 * ta_yy_xxxyz_0[i] * fe_0 - 3.0 * ta_yy_xxxyz_1[i] * fe_0 + ta_yyy_xxxz_0[i] * fe_0 - ta_yyy_xxxz_1[i] * fe_0 + ta_yyy_xxxyz_0[i] * pa_y[i] - ta_yyy_xxxyz_1[i] * pc_y[i];

        ta_yyyy_xxxzz_0[i] = 3.0 * ta_yy_xxxzz_0[i] * fe_0 - 3.0 * ta_yy_xxxzz_1[i] * fe_0 + ta_yyy_xxxzz_0[i] * pa_y[i] - ta_yyy_xxxzz_1[i] * pc_y[i];

        ta_yyyy_xxyyy_0[i] = 3.0 * ta_yy_xxyyy_0[i] * fe_0 - 3.0 * ta_yy_xxyyy_1[i] * fe_0 + 3.0 * ta_yyy_xxyy_0[i] * fe_0 - 3.0 * ta_yyy_xxyy_1[i] * fe_0 + ta_yyy_xxyyy_0[i] * pa_y[i] - ta_yyy_xxyyy_1[i] * pc_y[i];

        ta_yyyy_xxyyz_0[i] = 3.0 * ta_yy_xxyyz_0[i] * fe_0 - 3.0 * ta_yy_xxyyz_1[i] * fe_0 + 2.0 * ta_yyy_xxyz_0[i] * fe_0 - 2.0 * ta_yyy_xxyz_1[i] * fe_0 + ta_yyy_xxyyz_0[i] * pa_y[i] - ta_yyy_xxyyz_1[i] * pc_y[i];

        ta_yyyy_xxyzz_0[i] = 3.0 * ta_yy_xxyzz_0[i] * fe_0 - 3.0 * ta_yy_xxyzz_1[i] * fe_0 + ta_yyy_xxzz_0[i] * fe_0 - ta_yyy_xxzz_1[i] * fe_0 + ta_yyy_xxyzz_0[i] * pa_y[i] - ta_yyy_xxyzz_1[i] * pc_y[i];

        ta_yyyy_xxzzz_0[i] = 3.0 * ta_yy_xxzzz_0[i] * fe_0 - 3.0 * ta_yy_xxzzz_1[i] * fe_0 + ta_yyy_xxzzz_0[i] * pa_y[i] - ta_yyy_xxzzz_1[i] * pc_y[i];

        ta_yyyy_xyyyy_0[i] = 3.0 * ta_yy_xyyyy_0[i] * fe_0 - 3.0 * ta_yy_xyyyy_1[i] * fe_0 + 4.0 * ta_yyy_xyyy_0[i] * fe_0 - 4.0 * ta_yyy_xyyy_1[i] * fe_0 + ta_yyy_xyyyy_0[i] * pa_y[i] - ta_yyy_xyyyy_1[i] * pc_y[i];

        ta_yyyy_xyyyz_0[i] = 3.0 * ta_yy_xyyyz_0[i] * fe_0 - 3.0 * ta_yy_xyyyz_1[i] * fe_0 + 3.0 * ta_yyy_xyyz_0[i] * fe_0 - 3.0 * ta_yyy_xyyz_1[i] * fe_0 + ta_yyy_xyyyz_0[i] * pa_y[i] - ta_yyy_xyyyz_1[i] * pc_y[i];

        ta_yyyy_xyyzz_0[i] = 3.0 * ta_yy_xyyzz_0[i] * fe_0 - 3.0 * ta_yy_xyyzz_1[i] * fe_0 + 2.0 * ta_yyy_xyzz_0[i] * fe_0 - 2.0 * ta_yyy_xyzz_1[i] * fe_0 + ta_yyy_xyyzz_0[i] * pa_y[i] - ta_yyy_xyyzz_1[i] * pc_y[i];

        ta_yyyy_xyzzz_0[i] = 3.0 * ta_yy_xyzzz_0[i] * fe_0 - 3.0 * ta_yy_xyzzz_1[i] * fe_0 + ta_yyy_xzzz_0[i] * fe_0 - ta_yyy_xzzz_1[i] * fe_0 + ta_yyy_xyzzz_0[i] * pa_y[i] - ta_yyy_xyzzz_1[i] * pc_y[i];

        ta_yyyy_xzzzz_0[i] = 3.0 * ta_yy_xzzzz_0[i] * fe_0 - 3.0 * ta_yy_xzzzz_1[i] * fe_0 + ta_yyy_xzzzz_0[i] * pa_y[i] - ta_yyy_xzzzz_1[i] * pc_y[i];

        ta_yyyy_yyyyy_0[i] = 3.0 * ta_yy_yyyyy_0[i] * fe_0 - 3.0 * ta_yy_yyyyy_1[i] * fe_0 + 5.0 * ta_yyy_yyyy_0[i] * fe_0 - 5.0 * ta_yyy_yyyy_1[i] * fe_0 + ta_yyy_yyyyy_0[i] * pa_y[i] - ta_yyy_yyyyy_1[i] * pc_y[i];

        ta_yyyy_yyyyz_0[i] = 3.0 * ta_yy_yyyyz_0[i] * fe_0 - 3.0 * ta_yy_yyyyz_1[i] * fe_0 + 4.0 * ta_yyy_yyyz_0[i] * fe_0 - 4.0 * ta_yyy_yyyz_1[i] * fe_0 + ta_yyy_yyyyz_0[i] * pa_y[i] - ta_yyy_yyyyz_1[i] * pc_y[i];

        ta_yyyy_yyyzz_0[i] = 3.0 * ta_yy_yyyzz_0[i] * fe_0 - 3.0 * ta_yy_yyyzz_1[i] * fe_0 + 3.0 * ta_yyy_yyzz_0[i] * fe_0 - 3.0 * ta_yyy_yyzz_1[i] * fe_0 + ta_yyy_yyyzz_0[i] * pa_y[i] - ta_yyy_yyyzz_1[i] * pc_y[i];

        ta_yyyy_yyzzz_0[i] = 3.0 * ta_yy_yyzzz_0[i] * fe_0 - 3.0 * ta_yy_yyzzz_1[i] * fe_0 + 2.0 * ta_yyy_yzzz_0[i] * fe_0 - 2.0 * ta_yyy_yzzz_1[i] * fe_0 + ta_yyy_yyzzz_0[i] * pa_y[i] - ta_yyy_yyzzz_1[i] * pc_y[i];

        ta_yyyy_yzzzz_0[i] = 3.0 * ta_yy_yzzzz_0[i] * fe_0 - 3.0 * ta_yy_yzzzz_1[i] * fe_0 + ta_yyy_zzzz_0[i] * fe_0 - ta_yyy_zzzz_1[i] * fe_0 + ta_yyy_yzzzz_0[i] * pa_y[i] - ta_yyy_yzzzz_1[i] * pc_y[i];

        ta_yyyy_zzzzz_0[i] = 3.0 * ta_yy_zzzzz_0[i] * fe_0 - 3.0 * ta_yy_zzzzz_1[i] * fe_0 + ta_yyy_zzzzz_0[i] * pa_y[i] - ta_yyy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 231-252 components of targeted buffer : GH

    auto ta_yyyz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 231);

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

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yyy_xxxxx_0, ta_yyy_xxxxx_1, ta_yyy_xxxxy_0, ta_yyy_xxxxy_1, ta_yyy_xxxy_0, ta_yyy_xxxy_1, ta_yyy_xxxyy_0, ta_yyy_xxxyy_1, ta_yyy_xxxyz_0, ta_yyy_xxxyz_1, ta_yyy_xxyy_0, ta_yyy_xxyy_1, ta_yyy_xxyyy_0, ta_yyy_xxyyy_1, ta_yyy_xxyyz_0, ta_yyy_xxyyz_1, ta_yyy_xxyz_0, ta_yyy_xxyz_1, ta_yyy_xxyzz_0, ta_yyy_xxyzz_1, ta_yyy_xyyy_0, ta_yyy_xyyy_1, ta_yyy_xyyyy_0, ta_yyy_xyyyy_1, ta_yyy_xyyyz_0, ta_yyy_xyyyz_1, ta_yyy_xyyz_0, ta_yyy_xyyz_1, ta_yyy_xyyzz_0, ta_yyy_xyyzz_1, ta_yyy_xyzz_0, ta_yyy_xyzz_1, ta_yyy_xyzzz_0, ta_yyy_xyzzz_1, ta_yyy_yyyy_0, ta_yyy_yyyy_1, ta_yyy_yyyyy_0, ta_yyy_yyyyy_1, ta_yyy_yyyyz_0, ta_yyy_yyyyz_1, ta_yyy_yyyz_0, ta_yyy_yyyz_1, ta_yyy_yyyzz_0, ta_yyy_yyyzz_1, ta_yyy_yyzz_0, ta_yyy_yyzz_1, ta_yyy_yyzzz_0, ta_yyy_yyzzz_1, ta_yyy_yzzz_0, ta_yyy_yzzz_1, ta_yyy_yzzzz_0, ta_yyy_yzzzz_1, ta_yyyz_xxxxx_0, ta_yyyz_xxxxy_0, ta_yyyz_xxxxz_0, ta_yyyz_xxxyy_0, ta_yyyz_xxxyz_0, ta_yyyz_xxxzz_0, ta_yyyz_xxyyy_0, ta_yyyz_xxyyz_0, ta_yyyz_xxyzz_0, ta_yyyz_xxzzz_0, ta_yyyz_xyyyy_0, ta_yyyz_xyyyz_0, ta_yyyz_xyyzz_0, ta_yyyz_xyzzz_0, ta_yyyz_xzzzz_0, ta_yyyz_yyyyy_0, ta_yyyz_yyyyz_0, ta_yyyz_yyyzz_0, ta_yyyz_yyzzz_0, ta_yyyz_yzzzz_0, ta_yyyz_zzzzz_0, ta_yyz_xxxxz_0, ta_yyz_xxxxz_1, ta_yyz_xxxzz_0, ta_yyz_xxxzz_1, ta_yyz_xxzzz_0, ta_yyz_xxzzz_1, ta_yyz_xzzzz_0, ta_yyz_xzzzz_1, ta_yyz_zzzzz_0, ta_yyz_zzzzz_1, ta_yz_xxxxz_0, ta_yz_xxxxz_1, ta_yz_xxxzz_0, ta_yz_xxxzz_1, ta_yz_xxzzz_0, ta_yz_xxzzz_1, ta_yz_xzzzz_0, ta_yz_xzzzz_1, ta_yz_zzzzz_0, ta_yz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyz_xxxxx_0[i] = ta_yyy_xxxxx_0[i] * pa_z[i] - ta_yyy_xxxxx_1[i] * pc_z[i];

        ta_yyyz_xxxxy_0[i] = ta_yyy_xxxxy_0[i] * pa_z[i] - ta_yyy_xxxxy_1[i] * pc_z[i];

        ta_yyyz_xxxxz_0[i] = 2.0 * ta_yz_xxxxz_0[i] * fe_0 - 2.0 * ta_yz_xxxxz_1[i] * fe_0 + ta_yyz_xxxxz_0[i] * pa_y[i] - ta_yyz_xxxxz_1[i] * pc_y[i];

        ta_yyyz_xxxyy_0[i] = ta_yyy_xxxyy_0[i] * pa_z[i] - ta_yyy_xxxyy_1[i] * pc_z[i];

        ta_yyyz_xxxyz_0[i] = ta_yyy_xxxy_0[i] * fe_0 - ta_yyy_xxxy_1[i] * fe_0 + ta_yyy_xxxyz_0[i] * pa_z[i] - ta_yyy_xxxyz_1[i] * pc_z[i];

        ta_yyyz_xxxzz_0[i] = 2.0 * ta_yz_xxxzz_0[i] * fe_0 - 2.0 * ta_yz_xxxzz_1[i] * fe_0 + ta_yyz_xxxzz_0[i] * pa_y[i] - ta_yyz_xxxzz_1[i] * pc_y[i];

        ta_yyyz_xxyyy_0[i] = ta_yyy_xxyyy_0[i] * pa_z[i] - ta_yyy_xxyyy_1[i] * pc_z[i];

        ta_yyyz_xxyyz_0[i] = ta_yyy_xxyy_0[i] * fe_0 - ta_yyy_xxyy_1[i] * fe_0 + ta_yyy_xxyyz_0[i] * pa_z[i] - ta_yyy_xxyyz_1[i] * pc_z[i];

        ta_yyyz_xxyzz_0[i] = 2.0 * ta_yyy_xxyz_0[i] * fe_0 - 2.0 * ta_yyy_xxyz_1[i] * fe_0 + ta_yyy_xxyzz_0[i] * pa_z[i] - ta_yyy_xxyzz_1[i] * pc_z[i];

        ta_yyyz_xxzzz_0[i] = 2.0 * ta_yz_xxzzz_0[i] * fe_0 - 2.0 * ta_yz_xxzzz_1[i] * fe_0 + ta_yyz_xxzzz_0[i] * pa_y[i] - ta_yyz_xxzzz_1[i] * pc_y[i];

        ta_yyyz_xyyyy_0[i] = ta_yyy_xyyyy_0[i] * pa_z[i] - ta_yyy_xyyyy_1[i] * pc_z[i];

        ta_yyyz_xyyyz_0[i] = ta_yyy_xyyy_0[i] * fe_0 - ta_yyy_xyyy_1[i] * fe_0 + ta_yyy_xyyyz_0[i] * pa_z[i] - ta_yyy_xyyyz_1[i] * pc_z[i];

        ta_yyyz_xyyzz_0[i] = 2.0 * ta_yyy_xyyz_0[i] * fe_0 - 2.0 * ta_yyy_xyyz_1[i] * fe_0 + ta_yyy_xyyzz_0[i] * pa_z[i] - ta_yyy_xyyzz_1[i] * pc_z[i];

        ta_yyyz_xyzzz_0[i] = 3.0 * ta_yyy_xyzz_0[i] * fe_0 - 3.0 * ta_yyy_xyzz_1[i] * fe_0 + ta_yyy_xyzzz_0[i] * pa_z[i] - ta_yyy_xyzzz_1[i] * pc_z[i];

        ta_yyyz_xzzzz_0[i] = 2.0 * ta_yz_xzzzz_0[i] * fe_0 - 2.0 * ta_yz_xzzzz_1[i] * fe_0 + ta_yyz_xzzzz_0[i] * pa_y[i] - ta_yyz_xzzzz_1[i] * pc_y[i];

        ta_yyyz_yyyyy_0[i] = ta_yyy_yyyyy_0[i] * pa_z[i] - ta_yyy_yyyyy_1[i] * pc_z[i];

        ta_yyyz_yyyyz_0[i] = ta_yyy_yyyy_0[i] * fe_0 - ta_yyy_yyyy_1[i] * fe_0 + ta_yyy_yyyyz_0[i] * pa_z[i] - ta_yyy_yyyyz_1[i] * pc_z[i];

        ta_yyyz_yyyzz_0[i] = 2.0 * ta_yyy_yyyz_0[i] * fe_0 - 2.0 * ta_yyy_yyyz_1[i] * fe_0 + ta_yyy_yyyzz_0[i] * pa_z[i] - ta_yyy_yyyzz_1[i] * pc_z[i];

        ta_yyyz_yyzzz_0[i] = 3.0 * ta_yyy_yyzz_0[i] * fe_0 - 3.0 * ta_yyy_yyzz_1[i] * fe_0 + ta_yyy_yyzzz_0[i] * pa_z[i] - ta_yyy_yyzzz_1[i] * pc_z[i];

        ta_yyyz_yzzzz_0[i] = 4.0 * ta_yyy_yzzz_0[i] * fe_0 - 4.0 * ta_yyy_yzzz_1[i] * fe_0 + ta_yyy_yzzzz_0[i] * pa_z[i] - ta_yyy_yzzzz_1[i] * pc_z[i];

        ta_yyyz_zzzzz_0[i] = 2.0 * ta_yz_zzzzz_0[i] * fe_0 - 2.0 * ta_yz_zzzzz_1[i] * fe_0 + ta_yyz_zzzzz_0[i] * pa_y[i] - ta_yyz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 252-273 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yyz_xxxxy_0, ta_yyz_xxxxy_1, ta_yyz_xxxyy_0, ta_yyz_xxxyy_1, ta_yyz_xxyyy_0, ta_yyz_xxyyy_1, ta_yyz_xyyyy_0, ta_yyz_xyyyy_1, ta_yyz_yyyyy_0, ta_yyz_yyyyy_1, ta_yyzz_xxxxx_0, ta_yyzz_xxxxy_0, ta_yyzz_xxxxz_0, ta_yyzz_xxxyy_0, ta_yyzz_xxxyz_0, ta_yyzz_xxxzz_0, ta_yyzz_xxyyy_0, ta_yyzz_xxyyz_0, ta_yyzz_xxyzz_0, ta_yyzz_xxzzz_0, ta_yyzz_xyyyy_0, ta_yyzz_xyyyz_0, ta_yyzz_xyyzz_0, ta_yyzz_xyzzz_0, ta_yyzz_xzzzz_0, ta_yyzz_yyyyy_0, ta_yyzz_yyyyz_0, ta_yyzz_yyyzz_0, ta_yyzz_yyzzz_0, ta_yyzz_yzzzz_0, ta_yyzz_zzzzz_0, ta_yzz_xxxxx_0, ta_yzz_xxxxx_1, ta_yzz_xxxxz_0, ta_yzz_xxxxz_1, ta_yzz_xxxyz_0, ta_yzz_xxxyz_1, ta_yzz_xxxz_0, ta_yzz_xxxz_1, ta_yzz_xxxzz_0, ta_yzz_xxxzz_1, ta_yzz_xxyyz_0, ta_yzz_xxyyz_1, ta_yzz_xxyz_0, ta_yzz_xxyz_1, ta_yzz_xxyzz_0, ta_yzz_xxyzz_1, ta_yzz_xxzz_0, ta_yzz_xxzz_1, ta_yzz_xxzzz_0, ta_yzz_xxzzz_1, ta_yzz_xyyyz_0, ta_yzz_xyyyz_1, ta_yzz_xyyz_0, ta_yzz_xyyz_1, ta_yzz_xyyzz_0, ta_yzz_xyyzz_1, ta_yzz_xyzz_0, ta_yzz_xyzz_1, ta_yzz_xyzzz_0, ta_yzz_xyzzz_1, ta_yzz_xzzz_0, ta_yzz_xzzz_1, ta_yzz_xzzzz_0, ta_yzz_xzzzz_1, ta_yzz_yyyyz_0, ta_yzz_yyyyz_1, ta_yzz_yyyz_0, ta_yzz_yyyz_1, ta_yzz_yyyzz_0, ta_yzz_yyyzz_1, ta_yzz_yyzz_0, ta_yzz_yyzz_1, ta_yzz_yyzzz_0, ta_yzz_yyzzz_1, ta_yzz_yzzz_0, ta_yzz_yzzz_1, ta_yzz_yzzzz_0, ta_yzz_yzzzz_1, ta_yzz_zzzz_0, ta_yzz_zzzz_1, ta_yzz_zzzzz_0, ta_yzz_zzzzz_1, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzz_xxxxx_0[i] = ta_zz_xxxxx_0[i] * fe_0 - ta_zz_xxxxx_1[i] * fe_0 + ta_yzz_xxxxx_0[i] * pa_y[i] - ta_yzz_xxxxx_1[i] * pc_y[i];

        ta_yyzz_xxxxy_0[i] = ta_yy_xxxxy_0[i] * fe_0 - ta_yy_xxxxy_1[i] * fe_0 + ta_yyz_xxxxy_0[i] * pa_z[i] - ta_yyz_xxxxy_1[i] * pc_z[i];

        ta_yyzz_xxxxz_0[i] = ta_zz_xxxxz_0[i] * fe_0 - ta_zz_xxxxz_1[i] * fe_0 + ta_yzz_xxxxz_0[i] * pa_y[i] - ta_yzz_xxxxz_1[i] * pc_y[i];

        ta_yyzz_xxxyy_0[i] = ta_yy_xxxyy_0[i] * fe_0 - ta_yy_xxxyy_1[i] * fe_0 + ta_yyz_xxxyy_0[i] * pa_z[i] - ta_yyz_xxxyy_1[i] * pc_z[i];

        ta_yyzz_xxxyz_0[i] = ta_zz_xxxyz_0[i] * fe_0 - ta_zz_xxxyz_1[i] * fe_0 + ta_yzz_xxxz_0[i] * fe_0 - ta_yzz_xxxz_1[i] * fe_0 + ta_yzz_xxxyz_0[i] * pa_y[i] - ta_yzz_xxxyz_1[i] * pc_y[i];

        ta_yyzz_xxxzz_0[i] = ta_zz_xxxzz_0[i] * fe_0 - ta_zz_xxxzz_1[i] * fe_0 + ta_yzz_xxxzz_0[i] * pa_y[i] - ta_yzz_xxxzz_1[i] * pc_y[i];

        ta_yyzz_xxyyy_0[i] = ta_yy_xxyyy_0[i] * fe_0 - ta_yy_xxyyy_1[i] * fe_0 + ta_yyz_xxyyy_0[i] * pa_z[i] - ta_yyz_xxyyy_1[i] * pc_z[i];

        ta_yyzz_xxyyz_0[i] = ta_zz_xxyyz_0[i] * fe_0 - ta_zz_xxyyz_1[i] * fe_0 + 2.0 * ta_yzz_xxyz_0[i] * fe_0 - 2.0 * ta_yzz_xxyz_1[i] * fe_0 + ta_yzz_xxyyz_0[i] * pa_y[i] - ta_yzz_xxyyz_1[i] * pc_y[i];

        ta_yyzz_xxyzz_0[i] = ta_zz_xxyzz_0[i] * fe_0 - ta_zz_xxyzz_1[i] * fe_0 + ta_yzz_xxzz_0[i] * fe_0 - ta_yzz_xxzz_1[i] * fe_0 + ta_yzz_xxyzz_0[i] * pa_y[i] - ta_yzz_xxyzz_1[i] * pc_y[i];

        ta_yyzz_xxzzz_0[i] = ta_zz_xxzzz_0[i] * fe_0 - ta_zz_xxzzz_1[i] * fe_0 + ta_yzz_xxzzz_0[i] * pa_y[i] - ta_yzz_xxzzz_1[i] * pc_y[i];

        ta_yyzz_xyyyy_0[i] = ta_yy_xyyyy_0[i] * fe_0 - ta_yy_xyyyy_1[i] * fe_0 + ta_yyz_xyyyy_0[i] * pa_z[i] - ta_yyz_xyyyy_1[i] * pc_z[i];

        ta_yyzz_xyyyz_0[i] = ta_zz_xyyyz_0[i] * fe_0 - ta_zz_xyyyz_1[i] * fe_0 + 3.0 * ta_yzz_xyyz_0[i] * fe_0 - 3.0 * ta_yzz_xyyz_1[i] * fe_0 + ta_yzz_xyyyz_0[i] * pa_y[i] - ta_yzz_xyyyz_1[i] * pc_y[i];

        ta_yyzz_xyyzz_0[i] = ta_zz_xyyzz_0[i] * fe_0 - ta_zz_xyyzz_1[i] * fe_0 + 2.0 * ta_yzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzz_1[i] * fe_0 + ta_yzz_xyyzz_0[i] * pa_y[i] - ta_yzz_xyyzz_1[i] * pc_y[i];

        ta_yyzz_xyzzz_0[i] = ta_zz_xyzzz_0[i] * fe_0 - ta_zz_xyzzz_1[i] * fe_0 + ta_yzz_xzzz_0[i] * fe_0 - ta_yzz_xzzz_1[i] * fe_0 + ta_yzz_xyzzz_0[i] * pa_y[i] - ta_yzz_xyzzz_1[i] * pc_y[i];

        ta_yyzz_xzzzz_0[i] = ta_zz_xzzzz_0[i] * fe_0 - ta_zz_xzzzz_1[i] * fe_0 + ta_yzz_xzzzz_0[i] * pa_y[i] - ta_yzz_xzzzz_1[i] * pc_y[i];

        ta_yyzz_yyyyy_0[i] = ta_yy_yyyyy_0[i] * fe_0 - ta_yy_yyyyy_1[i] * fe_0 + ta_yyz_yyyyy_0[i] * pa_z[i] - ta_yyz_yyyyy_1[i] * pc_z[i];

        ta_yyzz_yyyyz_0[i] = ta_zz_yyyyz_0[i] * fe_0 - ta_zz_yyyyz_1[i] * fe_0 + 4.0 * ta_yzz_yyyz_0[i] * fe_0 - 4.0 * ta_yzz_yyyz_1[i] * fe_0 + ta_yzz_yyyyz_0[i] * pa_y[i] - ta_yzz_yyyyz_1[i] * pc_y[i];

        ta_yyzz_yyyzz_0[i] = ta_zz_yyyzz_0[i] * fe_0 - ta_zz_yyyzz_1[i] * fe_0 + 3.0 * ta_yzz_yyzz_0[i] * fe_0 - 3.0 * ta_yzz_yyzz_1[i] * fe_0 + ta_yzz_yyyzz_0[i] * pa_y[i] - ta_yzz_yyyzz_1[i] * pc_y[i];

        ta_yyzz_yyzzz_0[i] = ta_zz_yyzzz_0[i] * fe_0 - ta_zz_yyzzz_1[i] * fe_0 + 2.0 * ta_yzz_yzzz_0[i] * fe_0 - 2.0 * ta_yzz_yzzz_1[i] * fe_0 + ta_yzz_yyzzz_0[i] * pa_y[i] - ta_yzz_yyzzz_1[i] * pc_y[i];

        ta_yyzz_yzzzz_0[i] = ta_zz_yzzzz_0[i] * fe_0 - ta_zz_yzzzz_1[i] * fe_0 + ta_yzz_zzzz_0[i] * fe_0 - ta_yzz_zzzz_1[i] * fe_0 + ta_yzz_yzzzz_0[i] * pa_y[i] - ta_yzz_yzzzz_1[i] * pc_y[i];

        ta_yyzz_zzzzz_0[i] = ta_zz_zzzzz_0[i] * fe_0 - ta_zz_zzzzz_1[i] * fe_0 + ta_yzz_zzzzz_0[i] * pa_y[i] - ta_yzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 273-294 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_y, pc_y, ta_yzzz_xxxxx_0, ta_yzzz_xxxxy_0, ta_yzzz_xxxxz_0, ta_yzzz_xxxyy_0, ta_yzzz_xxxyz_0, ta_yzzz_xxxzz_0, ta_yzzz_xxyyy_0, ta_yzzz_xxyyz_0, ta_yzzz_xxyzz_0, ta_yzzz_xxzzz_0, ta_yzzz_xyyyy_0, ta_yzzz_xyyyz_0, ta_yzzz_xyyzz_0, ta_yzzz_xyzzz_0, ta_yzzz_xzzzz_0, ta_yzzz_yyyyy_0, ta_yzzz_yyyyz_0, ta_yzzz_yyyzz_0, ta_yzzz_yyzzz_0, ta_yzzz_yzzzz_0, ta_yzzz_zzzzz_0, ta_zzz_xxxx_0, ta_zzz_xxxx_1, ta_zzz_xxxxx_0, ta_zzz_xxxxx_1, ta_zzz_xxxxy_0, ta_zzz_xxxxy_1, ta_zzz_xxxxz_0, ta_zzz_xxxxz_1, ta_zzz_xxxy_0, ta_zzz_xxxy_1, ta_zzz_xxxyy_0, ta_zzz_xxxyy_1, ta_zzz_xxxyz_0, ta_zzz_xxxyz_1, ta_zzz_xxxz_0, ta_zzz_xxxz_1, ta_zzz_xxxzz_0, ta_zzz_xxxzz_1, ta_zzz_xxyy_0, ta_zzz_xxyy_1, ta_zzz_xxyyy_0, ta_zzz_xxyyy_1, ta_zzz_xxyyz_0, ta_zzz_xxyyz_1, ta_zzz_xxyz_0, ta_zzz_xxyz_1, ta_zzz_xxyzz_0, ta_zzz_xxyzz_1, ta_zzz_xxzz_0, ta_zzz_xxzz_1, ta_zzz_xxzzz_0, ta_zzz_xxzzz_1, ta_zzz_xyyy_0, ta_zzz_xyyy_1, ta_zzz_xyyyy_0, ta_zzz_xyyyy_1, ta_zzz_xyyyz_0, ta_zzz_xyyyz_1, ta_zzz_xyyz_0, ta_zzz_xyyz_1, ta_zzz_xyyzz_0, ta_zzz_xyyzz_1, ta_zzz_xyzz_0, ta_zzz_xyzz_1, ta_zzz_xyzzz_0, ta_zzz_xyzzz_1, ta_zzz_xzzz_0, ta_zzz_xzzz_1, ta_zzz_xzzzz_0, ta_zzz_xzzzz_1, ta_zzz_yyyy_0, ta_zzz_yyyy_1, ta_zzz_yyyyy_0, ta_zzz_yyyyy_1, ta_zzz_yyyyz_0, ta_zzz_yyyyz_1, ta_zzz_yyyz_0, ta_zzz_yyyz_1, ta_zzz_yyyzz_0, ta_zzz_yyyzz_1, ta_zzz_yyzz_0, ta_zzz_yyzz_1, ta_zzz_yyzzz_0, ta_zzz_yyzzz_1, ta_zzz_yzzz_0, ta_zzz_yzzz_1, ta_zzz_yzzzz_0, ta_zzz_yzzzz_1, ta_zzz_zzzz_0, ta_zzz_zzzz_1, ta_zzz_zzzzz_0, ta_zzz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzz_xxxxx_0[i] = ta_zzz_xxxxx_0[i] * pa_y[i] - ta_zzz_xxxxx_1[i] * pc_y[i];

        ta_yzzz_xxxxy_0[i] = ta_zzz_xxxx_0[i] * fe_0 - ta_zzz_xxxx_1[i] * fe_0 + ta_zzz_xxxxy_0[i] * pa_y[i] - ta_zzz_xxxxy_1[i] * pc_y[i];

        ta_yzzz_xxxxz_0[i] = ta_zzz_xxxxz_0[i] * pa_y[i] - ta_zzz_xxxxz_1[i] * pc_y[i];

        ta_yzzz_xxxyy_0[i] = 2.0 * ta_zzz_xxxy_0[i] * fe_0 - 2.0 * ta_zzz_xxxy_1[i] * fe_0 + ta_zzz_xxxyy_0[i] * pa_y[i] - ta_zzz_xxxyy_1[i] * pc_y[i];

        ta_yzzz_xxxyz_0[i] = ta_zzz_xxxz_0[i] * fe_0 - ta_zzz_xxxz_1[i] * fe_0 + ta_zzz_xxxyz_0[i] * pa_y[i] - ta_zzz_xxxyz_1[i] * pc_y[i];

        ta_yzzz_xxxzz_0[i] = ta_zzz_xxxzz_0[i] * pa_y[i] - ta_zzz_xxxzz_1[i] * pc_y[i];

        ta_yzzz_xxyyy_0[i] = 3.0 * ta_zzz_xxyy_0[i] * fe_0 - 3.0 * ta_zzz_xxyy_1[i] * fe_0 + ta_zzz_xxyyy_0[i] * pa_y[i] - ta_zzz_xxyyy_1[i] * pc_y[i];

        ta_yzzz_xxyyz_0[i] = 2.0 * ta_zzz_xxyz_0[i] * fe_0 - 2.0 * ta_zzz_xxyz_1[i] * fe_0 + ta_zzz_xxyyz_0[i] * pa_y[i] - ta_zzz_xxyyz_1[i] * pc_y[i];

        ta_yzzz_xxyzz_0[i] = ta_zzz_xxzz_0[i] * fe_0 - ta_zzz_xxzz_1[i] * fe_0 + ta_zzz_xxyzz_0[i] * pa_y[i] - ta_zzz_xxyzz_1[i] * pc_y[i];

        ta_yzzz_xxzzz_0[i] = ta_zzz_xxzzz_0[i] * pa_y[i] - ta_zzz_xxzzz_1[i] * pc_y[i];

        ta_yzzz_xyyyy_0[i] = 4.0 * ta_zzz_xyyy_0[i] * fe_0 - 4.0 * ta_zzz_xyyy_1[i] * fe_0 + ta_zzz_xyyyy_0[i] * pa_y[i] - ta_zzz_xyyyy_1[i] * pc_y[i];

        ta_yzzz_xyyyz_0[i] = 3.0 * ta_zzz_xyyz_0[i] * fe_0 - 3.0 * ta_zzz_xyyz_1[i] * fe_0 + ta_zzz_xyyyz_0[i] * pa_y[i] - ta_zzz_xyyyz_1[i] * pc_y[i];

        ta_yzzz_xyyzz_0[i] = 2.0 * ta_zzz_xyzz_0[i] * fe_0 - 2.0 * ta_zzz_xyzz_1[i] * fe_0 + ta_zzz_xyyzz_0[i] * pa_y[i] - ta_zzz_xyyzz_1[i] * pc_y[i];

        ta_yzzz_xyzzz_0[i] = ta_zzz_xzzz_0[i] * fe_0 - ta_zzz_xzzz_1[i] * fe_0 + ta_zzz_xyzzz_0[i] * pa_y[i] - ta_zzz_xyzzz_1[i] * pc_y[i];

        ta_yzzz_xzzzz_0[i] = ta_zzz_xzzzz_0[i] * pa_y[i] - ta_zzz_xzzzz_1[i] * pc_y[i];

        ta_yzzz_yyyyy_0[i] = 5.0 * ta_zzz_yyyy_0[i] * fe_0 - 5.0 * ta_zzz_yyyy_1[i] * fe_0 + ta_zzz_yyyyy_0[i] * pa_y[i] - ta_zzz_yyyyy_1[i] * pc_y[i];

        ta_yzzz_yyyyz_0[i] = 4.0 * ta_zzz_yyyz_0[i] * fe_0 - 4.0 * ta_zzz_yyyz_1[i] * fe_0 + ta_zzz_yyyyz_0[i] * pa_y[i] - ta_zzz_yyyyz_1[i] * pc_y[i];

        ta_yzzz_yyyzz_0[i] = 3.0 * ta_zzz_yyzz_0[i] * fe_0 - 3.0 * ta_zzz_yyzz_1[i] * fe_0 + ta_zzz_yyyzz_0[i] * pa_y[i] - ta_zzz_yyyzz_1[i] * pc_y[i];

        ta_yzzz_yyzzz_0[i] = 2.0 * ta_zzz_yzzz_0[i] * fe_0 - 2.0 * ta_zzz_yzzz_1[i] * fe_0 + ta_zzz_yyzzz_0[i] * pa_y[i] - ta_zzz_yyzzz_1[i] * pc_y[i];

        ta_yzzz_yzzzz_0[i] = ta_zzz_zzzz_0[i] * fe_0 - ta_zzz_zzzz_1[i] * fe_0 + ta_zzz_yzzzz_0[i] * pa_y[i] - ta_zzz_yzzzz_1[i] * pc_y[i];

        ta_yzzz_zzzzz_0[i] = ta_zzz_zzzzz_0[i] * pa_y[i] - ta_zzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 294-315 components of targeted buffer : GH

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

    #pragma omp simd aligned(pa_z, pc_z, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, ta_zzz_xxxx_0, ta_zzz_xxxx_1, ta_zzz_xxxxx_0, ta_zzz_xxxxx_1, ta_zzz_xxxxy_0, ta_zzz_xxxxy_1, ta_zzz_xxxxz_0, ta_zzz_xxxxz_1, ta_zzz_xxxy_0, ta_zzz_xxxy_1, ta_zzz_xxxyy_0, ta_zzz_xxxyy_1, ta_zzz_xxxyz_0, ta_zzz_xxxyz_1, ta_zzz_xxxz_0, ta_zzz_xxxz_1, ta_zzz_xxxzz_0, ta_zzz_xxxzz_1, ta_zzz_xxyy_0, ta_zzz_xxyy_1, ta_zzz_xxyyy_0, ta_zzz_xxyyy_1, ta_zzz_xxyyz_0, ta_zzz_xxyyz_1, ta_zzz_xxyz_0, ta_zzz_xxyz_1, ta_zzz_xxyzz_0, ta_zzz_xxyzz_1, ta_zzz_xxzz_0, ta_zzz_xxzz_1, ta_zzz_xxzzz_0, ta_zzz_xxzzz_1, ta_zzz_xyyy_0, ta_zzz_xyyy_1, ta_zzz_xyyyy_0, ta_zzz_xyyyy_1, ta_zzz_xyyyz_0, ta_zzz_xyyyz_1, ta_zzz_xyyz_0, ta_zzz_xyyz_1, ta_zzz_xyyzz_0, ta_zzz_xyyzz_1, ta_zzz_xyzz_0, ta_zzz_xyzz_1, ta_zzz_xyzzz_0, ta_zzz_xyzzz_1, ta_zzz_xzzz_0, ta_zzz_xzzz_1, ta_zzz_xzzzz_0, ta_zzz_xzzzz_1, ta_zzz_yyyy_0, ta_zzz_yyyy_1, ta_zzz_yyyyy_0, ta_zzz_yyyyy_1, ta_zzz_yyyyz_0, ta_zzz_yyyyz_1, ta_zzz_yyyz_0, ta_zzz_yyyz_1, ta_zzz_yyyzz_0, ta_zzz_yyyzz_1, ta_zzz_yyzz_0, ta_zzz_yyzz_1, ta_zzz_yyzzz_0, ta_zzz_yyzzz_1, ta_zzz_yzzz_0, ta_zzz_yzzz_1, ta_zzz_yzzzz_0, ta_zzz_yzzzz_1, ta_zzz_zzzz_0, ta_zzz_zzzz_1, ta_zzz_zzzzz_0, ta_zzz_zzzzz_1, ta_zzzz_xxxxx_0, ta_zzzz_xxxxy_0, ta_zzzz_xxxxz_0, ta_zzzz_xxxyy_0, ta_zzzz_xxxyz_0, ta_zzzz_xxxzz_0, ta_zzzz_xxyyy_0, ta_zzzz_xxyyz_0, ta_zzzz_xxyzz_0, ta_zzzz_xxzzz_0, ta_zzzz_xyyyy_0, ta_zzzz_xyyyz_0, ta_zzzz_xyyzz_0, ta_zzzz_xyzzz_0, ta_zzzz_xzzzz_0, ta_zzzz_yyyyy_0, ta_zzzz_yyyyz_0, ta_zzzz_yyyzz_0, ta_zzzz_yyzzz_0, ta_zzzz_yzzzz_0, ta_zzzz_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzz_xxxxx_0[i] = 3.0 * ta_zz_xxxxx_0[i] * fe_0 - 3.0 * ta_zz_xxxxx_1[i] * fe_0 + ta_zzz_xxxxx_0[i] * pa_z[i] - ta_zzz_xxxxx_1[i] * pc_z[i];

        ta_zzzz_xxxxy_0[i] = 3.0 * ta_zz_xxxxy_0[i] * fe_0 - 3.0 * ta_zz_xxxxy_1[i] * fe_0 + ta_zzz_xxxxy_0[i] * pa_z[i] - ta_zzz_xxxxy_1[i] * pc_z[i];

        ta_zzzz_xxxxz_0[i] = 3.0 * ta_zz_xxxxz_0[i] * fe_0 - 3.0 * ta_zz_xxxxz_1[i] * fe_0 + ta_zzz_xxxx_0[i] * fe_0 - ta_zzz_xxxx_1[i] * fe_0 + ta_zzz_xxxxz_0[i] * pa_z[i] - ta_zzz_xxxxz_1[i] * pc_z[i];

        ta_zzzz_xxxyy_0[i] = 3.0 * ta_zz_xxxyy_0[i] * fe_0 - 3.0 * ta_zz_xxxyy_1[i] * fe_0 + ta_zzz_xxxyy_0[i] * pa_z[i] - ta_zzz_xxxyy_1[i] * pc_z[i];

        ta_zzzz_xxxyz_0[i] = 3.0 * ta_zz_xxxyz_0[i] * fe_0 - 3.0 * ta_zz_xxxyz_1[i] * fe_0 + ta_zzz_xxxy_0[i] * fe_0 - ta_zzz_xxxy_1[i] * fe_0 + ta_zzz_xxxyz_0[i] * pa_z[i] - ta_zzz_xxxyz_1[i] * pc_z[i];

        ta_zzzz_xxxzz_0[i] = 3.0 * ta_zz_xxxzz_0[i] * fe_0 - 3.0 * ta_zz_xxxzz_1[i] * fe_0 + 2.0 * ta_zzz_xxxz_0[i] * fe_0 - 2.0 * ta_zzz_xxxz_1[i] * fe_0 + ta_zzz_xxxzz_0[i] * pa_z[i] - ta_zzz_xxxzz_1[i] * pc_z[i];

        ta_zzzz_xxyyy_0[i] = 3.0 * ta_zz_xxyyy_0[i] * fe_0 - 3.0 * ta_zz_xxyyy_1[i] * fe_0 + ta_zzz_xxyyy_0[i] * pa_z[i] - ta_zzz_xxyyy_1[i] * pc_z[i];

        ta_zzzz_xxyyz_0[i] = 3.0 * ta_zz_xxyyz_0[i] * fe_0 - 3.0 * ta_zz_xxyyz_1[i] * fe_0 + ta_zzz_xxyy_0[i] * fe_0 - ta_zzz_xxyy_1[i] * fe_0 + ta_zzz_xxyyz_0[i] * pa_z[i] - ta_zzz_xxyyz_1[i] * pc_z[i];

        ta_zzzz_xxyzz_0[i] = 3.0 * ta_zz_xxyzz_0[i] * fe_0 - 3.0 * ta_zz_xxyzz_1[i] * fe_0 + 2.0 * ta_zzz_xxyz_0[i] * fe_0 - 2.0 * ta_zzz_xxyz_1[i] * fe_0 + ta_zzz_xxyzz_0[i] * pa_z[i] - ta_zzz_xxyzz_1[i] * pc_z[i];

        ta_zzzz_xxzzz_0[i] = 3.0 * ta_zz_xxzzz_0[i] * fe_0 - 3.0 * ta_zz_xxzzz_1[i] * fe_0 + 3.0 * ta_zzz_xxzz_0[i] * fe_0 - 3.0 * ta_zzz_xxzz_1[i] * fe_0 + ta_zzz_xxzzz_0[i] * pa_z[i] - ta_zzz_xxzzz_1[i] * pc_z[i];

        ta_zzzz_xyyyy_0[i] = 3.0 * ta_zz_xyyyy_0[i] * fe_0 - 3.0 * ta_zz_xyyyy_1[i] * fe_0 + ta_zzz_xyyyy_0[i] * pa_z[i] - ta_zzz_xyyyy_1[i] * pc_z[i];

        ta_zzzz_xyyyz_0[i] = 3.0 * ta_zz_xyyyz_0[i] * fe_0 - 3.0 * ta_zz_xyyyz_1[i] * fe_0 + ta_zzz_xyyy_0[i] * fe_0 - ta_zzz_xyyy_1[i] * fe_0 + ta_zzz_xyyyz_0[i] * pa_z[i] - ta_zzz_xyyyz_1[i] * pc_z[i];

        ta_zzzz_xyyzz_0[i] = 3.0 * ta_zz_xyyzz_0[i] * fe_0 - 3.0 * ta_zz_xyyzz_1[i] * fe_0 + 2.0 * ta_zzz_xyyz_0[i] * fe_0 - 2.0 * ta_zzz_xyyz_1[i] * fe_0 + ta_zzz_xyyzz_0[i] * pa_z[i] - ta_zzz_xyyzz_1[i] * pc_z[i];

        ta_zzzz_xyzzz_0[i] = 3.0 * ta_zz_xyzzz_0[i] * fe_0 - 3.0 * ta_zz_xyzzz_1[i] * fe_0 + 3.0 * ta_zzz_xyzz_0[i] * fe_0 - 3.0 * ta_zzz_xyzz_1[i] * fe_0 + ta_zzz_xyzzz_0[i] * pa_z[i] - ta_zzz_xyzzz_1[i] * pc_z[i];

        ta_zzzz_xzzzz_0[i] = 3.0 * ta_zz_xzzzz_0[i] * fe_0 - 3.0 * ta_zz_xzzzz_1[i] * fe_0 + 4.0 * ta_zzz_xzzz_0[i] * fe_0 - 4.0 * ta_zzz_xzzz_1[i] * fe_0 + ta_zzz_xzzzz_0[i] * pa_z[i] - ta_zzz_xzzzz_1[i] * pc_z[i];

        ta_zzzz_yyyyy_0[i] = 3.0 * ta_zz_yyyyy_0[i] * fe_0 - 3.0 * ta_zz_yyyyy_1[i] * fe_0 + ta_zzz_yyyyy_0[i] * pa_z[i] - ta_zzz_yyyyy_1[i] * pc_z[i];

        ta_zzzz_yyyyz_0[i] = 3.0 * ta_zz_yyyyz_0[i] * fe_0 - 3.0 * ta_zz_yyyyz_1[i] * fe_0 + ta_zzz_yyyy_0[i] * fe_0 - ta_zzz_yyyy_1[i] * fe_0 + ta_zzz_yyyyz_0[i] * pa_z[i] - ta_zzz_yyyyz_1[i] * pc_z[i];

        ta_zzzz_yyyzz_0[i] = 3.0 * ta_zz_yyyzz_0[i] * fe_0 - 3.0 * ta_zz_yyyzz_1[i] * fe_0 + 2.0 * ta_zzz_yyyz_0[i] * fe_0 - 2.0 * ta_zzz_yyyz_1[i] * fe_0 + ta_zzz_yyyzz_0[i] * pa_z[i] - ta_zzz_yyyzz_1[i] * pc_z[i];

        ta_zzzz_yyzzz_0[i] = 3.0 * ta_zz_yyzzz_0[i] * fe_0 - 3.0 * ta_zz_yyzzz_1[i] * fe_0 + 3.0 * ta_zzz_yyzz_0[i] * fe_0 - 3.0 * ta_zzz_yyzz_1[i] * fe_0 + ta_zzz_yyzzz_0[i] * pa_z[i] - ta_zzz_yyzzz_1[i] * pc_z[i];

        ta_zzzz_yzzzz_0[i] = 3.0 * ta_zz_yzzzz_0[i] * fe_0 - 3.0 * ta_zz_yzzzz_1[i] * fe_0 + 4.0 * ta_zzz_yzzz_0[i] * fe_0 - 4.0 * ta_zzz_yzzz_1[i] * fe_0 + ta_zzz_yzzzz_0[i] * pa_z[i] - ta_zzz_yzzzz_1[i] * pc_z[i];

        ta_zzzz_zzzzz_0[i] = 3.0 * ta_zz_zzzzz_0[i] * fe_0 - 3.0 * ta_zz_zzzzz_1[i] * fe_0 + 5.0 * ta_zzz_zzzz_0[i] * fe_0 - 5.0 * ta_zzz_zzzz_1[i] * fe_0 + ta_zzz_zzzzz_0[i] * pa_z[i] - ta_zzz_zzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

